"""
BioCompiler ESMFold Offline Fallback — Sequence-Based Heuristic Structure Prediction
=====================================================================================

Provides a lightweight, offline-capable fallback when the ESMFold API and
local ``esm`` package are both unavailable.  Instead of returning
``success=False`` (which causes structure predicates to return UNCERTAIN),
this module uses simple sequence-based heuristics to produce a
low-confidence estimate of structural properties.

**IMPORTANT — Accuracy Caveat**
    This fallback is *significantly* less accurate than ESMFold.  The
    heuristics here are based on well-known biophysical principles
    (Kyte-Doolittle hydrophobicity, Chou-Fasman secondary structure
    propensity, charge distribution), but they cannot capture the
    long-range interactions and evolutionary constraints that a deep
    learning model like ESMFold learns.  Results from this module should
    be treated as rough estimates only and are flagged with
    ``method="heuristic_fallback"`` and ``confidence < 0.5``.

Heuristic Components
--------------------
1. **Hydrophobicity profile** (Kyte-Doolittle):
   - Sliding-window average hydrophobicity per residue.
   - Hydrophobic regions are assumed to be buried (higher pLDDT).
   - Hydrophilic / charged regions are assumed to be exposed (lower pLDDT).

2. **Charge distribution**:
   - Fraction of positively / negatively charged residues.
   - Charge clustering (surface patches) lowers estimated confidence.
   - Balanced charge distribution is consistent with globular proteins.

3. **Secondary structure propensity** (simplified Chou-Fasman):
   - Per-residue helix, sheet, coil propensity scores.
   - Windowed classification into helix / sheet / coil regions.
   - Proteins with well-defined secondary structure are assigned
     slightly higher estimated pLDDT.

Estimated pLDDT Score
---------------------
The heuristic pLDDT is bounded by ``HEURISTIC_MAX_CONFIDENCE`` (default
40.0, well below the ESMFold "Low confidence" threshold of 50).  This
ensures that downstream predicates never treat heuristic results as
reliable.  The score is computed as::

    base_plddt = 25.0  (starting point for all sequences)

    + hydrophobicity_bonus     (0–8 points, for balanced hydrophobic content)
    + secondary_structure_bonus (0–5 points, for well-defined SS propensity)
    + charge_balance_bonus      (0–4 points, for balanced charge distribution)
    - disorder_penalty          (0–10 points, for long low-complexity runs)
    - charge_clustering_penalty (0–5 points, for concentrated same-charge patches)

    final_plddt = min(base_plddt + bonuses - penalties, HEURISTIC_MAX_CONFIDENCE)

When to Use
-----------
This fallback is invoked automatically by :func:`esmfold.predict_structure`
when both the ESM Atlas API and local ``esm`` package are unavailable.
It should **never** be used as a replacement for ESMFold when the API is
reachable.

References
----------
- Kyte & Doolittle, J. Mol. Biol. 1982; 157:105–132 (hydrophobicity scale)
- Chou & Fasman, Biochemistry 1974; 13:222–245 (secondary structure propensity)
- Lin et al., Science 2023; 379:1043 (ESMFold / ESM-2)
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "predict_structure_heuristic",
    "estimate_plddt_from_sequence",
    "estimate_secondary_structure_from_sequence",
    "compute_hydrophobicity_profile",
    "compute_charge_profile",
    "HEURISTIC_MAX_CONFIDENCE",
]


# ==============================================================================
# Constants
# ==============================================================================

#: Maximum pLDDT score the heuristic fallback will ever return.
#: Intentionally below the ESMFold "Low confidence" threshold of 50.
HEURISTIC_MAX_CONFIDENCE: float = 40.0

#: Kyte-Doolittle hydrophobicity scale (Kyte & Doolittle, 1982).
#: More positive = more hydrophobic.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

#: Simplified Chou-Fasman secondary structure propensity parameters.
#: Each amino acid has (helix_propensity, sheet_propensity).
#: Higher values = stronger propensity for that SS type.
#: Derived from Chou & Fasman (1974) and subsequent refinements.
CHOU_FASMAN_PROPENSITY: dict[str, tuple[float, float]] = {
    "A": (1.42, 0.83), "R": (0.98, 0.93), "N": (0.67, 0.89), "D": (1.01, 0.54),
    "C": (0.70, 1.19), "Q": (1.11, 1.10), "E": (1.51, 0.37), "G": (0.57, 0.75),
    "H": (1.00, 0.87), "I": (1.08, 1.60), "L": (1.21, 1.30), "K": (1.16, 0.74),
    "M": (1.45, 1.05), "F": (1.13, 1.38), "P": (0.57, 0.55), "S": (0.77, 0.75),
    "T": (0.83, 1.19), "W": (1.08, 1.37), "Y": (0.69, 1.47), "V": (1.06, 1.70),
}

POSITIVELY_CHARGED_AAS = set("KRH")
NEGATIVELY_CHARGED_AAS = set("DE")
HYDROPHOBIC_AAS = set("AILMFVW")

#: Threshold for Chou-Fasman helix nucleation (a window of 6 residues
#: with 4 or more above P_helix > 1.00 initiates a helix).
HELIX_NUCLEATION_WINDOW = 6
HELIX_NUCLEATION_THRESHOLD = 4

#: Threshold for Chou-Fasman sheet nucleation (a window of 5 residues
#: with 3 or more above P_sheet > 1.00 initiates a sheet).
SHEET_NUCLEATION_WINDOW = 5
SHEET_NUCLEATION_THRESHOLD = 3

#: Sliding window size for hydrophobicity smoothing.
HYDROPHOBICITY_WINDOW = 9


# ==============================================================================
# Hydrophobicity profile
# ==============================================================================

def compute_hydrophobicity_profile(
    protein: str,
    window: int = HYDROPHOBICITY_WINDOW,
) -> list[float]:
    """Compute a smoothed Kyte-Doolittle hydrophobicity profile.

    For each residue position, the average Kyte-Doolittle score over a
    sliding window centered on that position is computed.  Positions near
    the termini use a smaller window (asymmetric padding).

    Args:
        protein: Amino acid sequence (single-letter codes, uppercase).
        window:  Sliding window size (should be odd; default 9).

    Returns:
        List of smoothed hydrophobicity values, one per residue.
    """
    protein = protein.upper()
    n = len(protein)
    if n == 0:
        return []

    # Raw hydrophobicity per residue
    raw = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in protein]

    # Smooth with sliding window
    half_w = window // 2
    smoothed: list[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        avg = sum(raw[start:end]) / (end - start)
        smoothed.append(round(avg, 3))

    return smoothed


# ==============================================================================
# Charge profile
# ==============================================================================

@dataclass
class ChargeProfile:
    """Summary of charge distribution in a protein sequence.

    Attributes:
        net_charge:           Sum of positive (+1) and negative (-1) charges.
        positive_fraction:    Fraction of positively charged residues (KRH).
        negative_fraction:    Fraction of negatively charged residues (DE).
        charge_balance:       How balanced the charges are (1.0 = perfectly
                              balanced, 0.0 = all one sign).  Computed as
                              ``1.0 - abs(pos_frac - neg_frac) / max(pos_frac + neg_frac, 1e-9)``.
        charge_patch_count:   Number of concentrated same-charge patches
                              (7-residue windows with >=70% same-sign charge).
    """
    net_charge: int
    positive_fraction: float
    negative_fraction: float
    charge_balance: float
    charge_patch_count: int


def compute_charge_profile(protein: str, window: int = 7, min_frac: float = 0.7) -> ChargeProfile:
    """Compute charge distribution statistics for a protein sequence.

    Args:
        protein: Amino acid sequence (single-letter codes).
        window:  Window size for charge patch detection (default 7).
        min_frac: Minimum same-sign charge fraction to count as a patch (default 0.7).

    Returns:
        ChargeProfile with net charge, fractions, balance, and patch count.
    """
    protein = protein.upper()
    n = len(protein)

    pos_count = sum(1 for aa in protein if aa in POSITIVELY_CHARGED_AAS)
    neg_count = sum(1 for aa in protein if aa in NEGATIVELY_CHARGED_AAS)

    pos_frac = pos_count / n if n > 0 else 0.0
    neg_frac = neg_count / n if n > 0 else 0.0
    net_charge = pos_count - neg_count

    total_charged = pos_frac + neg_frac
    if total_charged > 1e-9:
        charge_balance = 1.0 - abs(pos_frac - neg_frac) / total_charged
    else:
        # No charged residues at all — neutral, but not informative
        charge_balance = 0.5

    # Detect charge patches
    patch_count = 0
    if n >= window:
        i = 0
        while i <= n - window:
            window_seq = protein[i:i + window]
            w_pos = sum(1 for aa in window_seq if aa in POSITIVELY_CHARGED_AAS)
            w_neg = sum(1 for aa in window_seq if aa in NEGATIVELY_CHARGED_AAS)
            if w_pos / window >= min_frac or w_neg / window >= min_frac:
                patch_count += 1
                i += window  # skip past this patch
            else:
                i += 1

    return ChargeProfile(
        net_charge=net_charge,
        positive_fraction=round(pos_frac, 4),
        negative_fraction=round(neg_frac, 4),
        charge_balance=round(charge_balance, 4),
        charge_patch_count=patch_count,
    )


# ==============================================================================
# Secondary structure estimation
# ==============================================================================

@dataclass
class SecondaryStructureEstimate:
    """Estimated secondary structure fractions from Chou-Fasman heuristics.

    Attributes:
        helix_fraction:   Fraction of residues predicted to be helical.
        sheet_fraction:   Fraction of residues predicted to be beta-sheet.
        coil_fraction:    Fraction of residues predicted to be coil.
        assignments:      Per-residue SS assignment ('H', 'E', 'C').
    """
    helix_fraction: float
    sheet_fraction: float
    coil_fraction: float
    assignments: list[str]


def estimate_secondary_structure_from_sequence(protein: str) -> SecondaryStructureEstimate:
    """Estimate secondary structure content from amino acid sequence.

    Uses a simplified Chou-Fasman algorithm:
      1. Compute per-residue helix and sheet propensities.
      2. Identify helix nucleation sites (window of 6 with >=4 residues
         above helix threshold 1.00).
      3. Identify sheet nucleation sites (window of 5 with >=3 residues
         above sheet threshold 1.00).
      4. Extend nucleated regions while average propensity stays above 1.00.
      5. Resolve overlaps (helix takes priority).

    This is a rough approximation.  Real Chou-Fasman includes additional
    rules for turn prediction, Pro/Gly breaking, and strand pairing.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        SecondaryStructureEstimate with fractions and per-residue assignments.
    """
    protein = protein.upper()
    n = len(protein)

    if n == 0:
        return SecondaryStructureEstimate(
            helix_fraction=0.0, sheet_fraction=0.0, coil_fraction=1.0,
            assignments=[],
        )

    # Per-residue propensities
    p_helix = [CHOU_FASMAN_PROPENSITY.get(aa, (1.0, 1.0))[0] for aa in protein]
    p_sheet = [CHOU_FASMAN_PROPENSITY.get(aa, (1.0, 1.0))[1] for aa in protein]

    # Assignments: default to coil
    assignments: list[str] = ["C"] * n

    # --- Helix nucleation and extension ---
    helix_regions: list[tuple[int, int]] = []
    i = 0
    while i <= n - HELIX_NUCLEATION_WINDOW:
        window_helix = p_helix[i:i + HELIX_NUCLEATION_WINDOW]
        above_threshold = sum(1 for p in window_helix if p > 1.00)
        if above_threshold >= HELIX_NUCLEATION_THRESHOLD:
            # Nucleate helix, extend in both directions
            start = i
            end = i + HELIX_NUCLEATION_WINDOW

            # Extend N-terminal
            while start > 0:
                avg_helix = sum(p_helix[max(0, start - 1):end]) / (end - start + 1)
                if avg_helix > 1.00:
                    start -= 1
                else:
                    break

            # Extend C-terminal
            while end < n:
                avg_helix = sum(p_helix[start:min(n, end + 1)]) / (end - start + 1)
                if avg_helix > 1.00:
                    end += 1
                else:
                    break

            helix_regions.append((start, end))
            i = end
        else:
            i += 1

    # --- Sheet nucleation and extension ---
    sheet_regions: list[tuple[int, int]] = []
    i = 0
    while i <= n - SHEET_NUCLEATION_WINDOW:
        window_sheet = p_sheet[i:i + SHEET_NUCLEATION_WINDOW]
        above_threshold = sum(1 for p in window_sheet if p > 1.00)
        if above_threshold >= SHEET_NUCLEATION_THRESHOLD:
            start = i
            end = i + SHEET_NUCLEATION_WINDOW

            # Extend N-terminal
            while start > 0:
                avg_sheet = sum(p_sheet[max(0, start - 1):end]) / (end - start + 1)
                if avg_sheet > 1.00:
                    start -= 1
                else:
                    break

            # Extend C-terminal
            while end < n:
                avg_sheet = sum(p_sheet[start:min(n, end + 1)]) / (end - start + 1)
                if avg_sheet > 1.00:
                    end += 1
                else:
                    break

            sheet_regions.append((start, end))
            i = end
        else:
            i += 1

    # Apply helix assignments (priority over sheet)
    for start, end in helix_regions:
        for j in range(start, min(end, n)):
            assignments[j] = "H"

    # Apply sheet assignments (only where not already helix)
    for start, end in sheet_regions:
        for j in range(start, min(end, n)):
            if assignments[j] == "C":
                assignments[j] = "E"

    # Compute fractions
    helix_count = sum(1 for a in assignments if a == "H")
    sheet_count = sum(1 for a in assignments if a == "E")
    coil_count = n - helix_count - sheet_count

    return SecondaryStructureEstimate(
        helix_fraction=round(helix_count / n, 4) if n > 0 else 0.0,
        sheet_fraction=round(sheet_count / n, 4) if n > 0 else 0.0,
        coil_fraction=round(coil_count / n, 4) if n > 0 else 1.0,
        assignments=assignments,
    )


# ==============================================================================
# Heuristic pLDDT estimation
# ==============================================================================

def estimate_plddt_from_sequence(protein: str) -> dict[str, Any]:
    """Estimate a heuristic pLDDT score from protein sequence alone.

    This function combines multiple sequence-based heuristics to produce
    a rough estimate of what ESMFold's mean pLDDT might be.  The estimate
    is capped at :data:`HEURISTIC_MAX_CONFIDENCE` (default 40.0) to
    ensure downstream code never treats it as a reliable prediction.

    Heuristic scoring breakdown:
        - **Base score**: 25.0 (very low starting point)
        - **Hydrophobicity bonus** (0–8): Balanced hydrophobic content
          (around 30–50% hydrophobic residues) suggests a well-folded
          globular protein.
        - **Secondary structure bonus** (0–5): Presence of nucleated
          helix or sheet regions suggests ordered structure.
        - **Charge balance bonus** (0–4): Roughly equal positive and
          negative charges is typical of soluble, folded proteins.
        - **Disorder penalty** (0–10): Long runs of low-complexity or
          disorder-promoting residues reduce confidence.
        - **Charge clustering penalty** (0–5): Concentrated same-charge
          patches suggest surface binding or disorder.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        Dictionary with keys:
            - ``estimated_mean_plddt`` (float): Heuristic pLDDT estimate
              (0–HEURISTIC_MAX_CONFIDENCE).
            - ``confidence`` (float): Confidence in this estimate (0–0.5).
            - ``method`` (str): Always ``"heuristic_fallback"``.
            - ``heuristic_details`` (dict): Breakdown of bonuses/penalties.
            - ``hydrophobicity_profile`` (list[float]): Per-residue smoothed
              Kyte-Doolittle values.
            - ``charge_profile`` (ChargeProfile): Charge distribution summary.
            - ``secondary_structure`` (SecondaryStructureEstimate): SS estimate.
    """
    protein = protein.upper()
    n = len(protein)

    if n == 0:
        return {
            "estimated_mean_plddt": 0.0,
            "confidence": 0.0,
            "method": "heuristic_fallback",
            "heuristic_details": {"base": 0.0},
            "hydrophobicity_profile": [],
            "charge_profile": ChargeProfile(
                net_charge=0, positive_fraction=0.0, negative_fraction=0.0,
                charge_balance=0.0, charge_patch_count=0,
            ),
            "secondary_structure": SecondaryStructureEstimate(
                helix_fraction=0.0, sheet_fraction=0.0, coil_fraction=1.0,
                assignments=[],
            ),
        }

    # --- Compute sub-profiles ---
    hydro_profile = compute_hydrophobicity_profile(protein)
    charge_prof = compute_charge_profile(protein)
    ss_estimate = estimate_secondary_structure_from_sequence(protein)

    # --- Base score ---
    base = 25.0

    # --- Hydrophobicity bonus (0–8) ---
    # Well-folded globular proteins typically have ~30–50% hydrophobic residues.
    hydrophobic_count = sum(1 for aa in protein if aa in HYDROPHOBIC_AAS)
    hydro_frac = hydrophobic_count / n
    # Optimal range is 0.30–0.50; penalize outside that.
    if 0.30 <= hydro_frac <= 0.50:
        hydro_bonus = 8.0
    elif 0.20 <= hydro_frac < 0.30 or 0.50 < hydro_frac <= 0.60:
        hydro_bonus = 5.0
    elif 0.10 <= hydro_frac < 0.20 or 0.60 < hydro_frac <= 0.70:
        hydro_bonus = 2.0
    else:
        hydro_bonus = 0.0

    # Also check mean hydrophobicity — positive mean suggests foldable.
    mean_hydro = sum(hydro_profile) / len(hydro_profile) if hydro_profile else 0.0
    if mean_hydro < -1.0:
        hydro_bonus = min(hydro_bonus, 2.0)  # Very hydrophilic — likely disordered

    # --- Secondary structure bonus (0–5) ---
    # More defined SS → higher bonus.
    defined_ss_frac = ss_estimate.helix_fraction + ss_estimate.sheet_fraction
    if defined_ss_frac > 0.60:
        ss_bonus = 5.0
    elif defined_ss_frac > 0.40:
        ss_bonus = 4.0
    elif defined_ss_frac > 0.25:
        ss_bonus = 3.0
    elif defined_ss_frac > 0.10:
        ss_bonus = 2.0
    elif defined_ss_frac > 0.0:
        ss_bonus = 1.0
    else:
        ss_bonus = 0.0

    # --- Charge balance bonus (0–4) ---
    # Balanced charges suggest soluble, well-folded protein.
    cb = charge_prof.charge_balance
    if cb > 0.8:
        charge_bonus = 4.0
    elif cb > 0.6:
        charge_bonus = 3.0
    elif cb > 0.4:
        charge_bonus = 2.0
    elif cb > 0.2:
        charge_bonus = 1.0
    else:
        charge_bonus = 0.0

    # --- Disorder penalty (0–10) ---
    # Long runs of disorder-promoting residues (P, G, S, Q, E, K).
    # Also penalize low-complexity regions (repeated same amino acid).
    disorder_promoting = set("PGSQEK")
    consecutive_disorder = 0
    max_consecutive_disorder = 0
    for aa in protein:
        if aa in disorder_promoting:
            consecutive_disorder += 1
            max_consecutive_disorder = max(max_consecutive_disorder, consecutive_disorder)
        else:
            consecutive_disorder = 0

    # Long disorder runs (>15 residues) are a strong signal of IDRs.
    if max_consecutive_disorder > 30:
        disorder_penalty = 10.0
    elif max_consecutive_disorder > 20:
        disorder_penalty = 7.0
    elif max_consecutive_disorder > 15:
        disorder_penalty = 5.0
    elif max_consecutive_disorder > 10:
        disorder_penalty = 3.0
    elif max_consecutive_disorder > 5:
        disorder_penalty = 1.0
    else:
        disorder_penalty = 0.0

    # Low-complexity penalty: long runs of same amino acid.
    max_repeat = 1
    current_repeat = 1
    for i in range(1, n):
        if protein[i] == protein[i - 1]:
            current_repeat += 1
            max_repeat = max(max_repeat, current_repeat)
        else:
            current_repeat = 1

    if max_repeat > 10:
        disorder_penalty = min(disorder_penalty + 3.0, 10.0)
    elif max_repeat > 6:
        disorder_penalty = min(disorder_penalty + 1.0, 10.0)

    # --- Charge clustering penalty (0–5) ---
    patch_count = charge_prof.charge_patch_count
    if patch_count >= 3:
        charge_cluster_penalty = 5.0
    elif patch_count >= 2:
        charge_cluster_penalty = 3.0
    elif patch_count >= 1:
        charge_cluster_penalty = 1.0
    else:
        charge_cluster_penalty = 0.0

    # --- Compute final estimated pLDDT ---
    raw_plddt = base + hydro_bonus + ss_bonus + charge_bonus - disorder_penalty - charge_cluster_penalty
    estimated_plddt = max(0.0, min(raw_plddt, HEURISTIC_MAX_CONFIDENCE))

    # --- Confidence in this estimate ---
    # Confidence is inherently low for heuristics.  We scale it based on
    # how many positive signals we found vs. total possible.
    max_bonus = 8.0 + 5.0 + 4.0  # hydro + ss + charge = 17
    max_penalty = 10.0 + 5.0  # disorder + charge_cluster = 15
    total_signal = hydro_bonus + ss_bonus + charge_bonus
    total_noise = disorder_penalty + charge_cluster_penalty

    # More positive signals → slightly higher confidence, but always < 0.5
    confidence = min(0.5, max(0.1, 0.15 + 0.02 * total_signal - 0.01 * total_noise))

    heuristic_details = {
        "base": base,
        "hydrophobicity_bonus": hydro_bonus,
        "hydrophobic_fraction": round(hydro_frac, 4),
        "mean_hydrophobicity": round(mean_hydro, 3),
        "secondary_structure_bonus": ss_bonus,
        "defined_ss_fraction": round(defined_ss_frac, 4),
        "helix_fraction": ss_estimate.helix_fraction,
        "sheet_fraction": ss_estimate.sheet_fraction,
        "charge_balance_bonus": charge_bonus,
        "charge_balance": charge_prof.charge_balance,
        "disorder_penalty": disorder_penalty,
        "max_consecutive_disorder": max_consecutive_disorder,
        "max_repeat": max_repeat,
        "charge_clustering_penalty": charge_cluster_penalty,
        "charge_patch_count": patch_count,
    }

    return {
        "estimated_mean_plddt": round(estimated_plddt, 2),
        "confidence": round(confidence, 4),
        "method": "heuristic_fallback",
        "heuristic_details": heuristic_details,
        "hydrophobicity_profile": hydro_profile,
        "charge_profile": charge_prof,
        "secondary_structure": ss_estimate,
    }


# ==============================================================================
# Main fallback prediction function
# ==============================================================================

def predict_structure_heuristic(protein: str) -> dict[str, Any]:
    """Produce a heuristic structure prediction for offline fallback.

    This function returns a dictionary that can be used to construct an
    :class:`~biocompiler.esmfold.ESMFoldResult` with low-confidence
    estimates.  It is called by :func:`~biocompiler.esmfold.predict_structure`
    when both the ESM Atlas API and the local ``esm`` package are
    unavailable.

    The returned prediction always has:
      - ``success = True`` (a prediction was produced, albeit low-confidence)
      - ``method = "heuristic_fallback"``
      - ``mean_plddt`` capped at :data:`HEURISTIC_MAX_CONFIDENCE`
      - ``confidence < 0.5``

    **Accuracy warning**: This heuristic is NOT a substitute for ESMFold.
    It is based on simple biophysical rules and can be wildly inaccurate
    for any specific protein.  It exists solely to provide a non-UNCERTAIN
    signal when no real structure prediction is available, so that
    downstream predicates can produce a tentative verdict rather than
    defaulting to UNCERTAIN.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        Dictionary with keys needed to build an ESMFoldResult:
            - ``protein`` (str): Input sequence.
            - ``mean_plddt`` (float): Estimated mean pLDDT (capped).
            - ``plddt_scores`` (list[float]): Per-residue estimated pLDDT.
            - ``method`` (str): ``"heuristic_fallback"``.
            - ``model_name`` (str): ``"heuristic_v1"``.
            - ``confidence`` (float): Confidence in the prediction (0–0.5).
            - ``heuristic_details`` (dict): Scoring breakdown.
            - ``secondary_structure`` (dict): Estimated SS fractions.
    """
    protein = protein.upper()
    n = len(protein)

    estimate = estimate_plddt_from_sequence(protein)
    mean_plddt = estimate["estimated_mean_plddt"]

    # Per-residue pLDDT: use hydrophobicity profile to modulate around the mean.
    # Hydrophobic residues get a slight boost, hydrophilic get a slight reduction.
    hydro_profile = estimate["hydrophobicity_profile"]
    ss_assignments = estimate["secondary_structure"].assignments

    plddt_scores: list[float] = []
    for i in range(n):
        # Base per-residue score starts from mean
        residue_plddt = mean_plddt

        # Hydrophobicity modulation: buried residues tend to have higher pLDDT
        if i < len(hydro_profile):
            # Positive hydrophobicity → slight boost; negative → slight reduction
            hydro_mod = hydro_profile[i] * 1.5  # scale factor
            residue_plddt += hydro_mod

        # SS modulation: defined SS gets a slight boost
        if i < len(ss_assignments):
            if ss_assignments[i] == "H":
                residue_plddt += 2.0  # helix residues are often well-predicted
            elif ss_assignments[i] == "E":
                residue_plddt += 1.5  # sheet residues are often well-predicted
            else:
                residue_plddt -= 1.0  # coil residues are harder to predict

        # Clamp each per-residue score to [0, HEURISTIC_MAX_CONFIDENCE]
        residue_plddt = max(0.0, min(residue_plddt, HEURISTIC_MAX_CONFIDENCE))
        plddt_scores.append(round(residue_plddt, 2))

    # Recompute mean from per-residue scores for consistency
    actual_mean = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0
    actual_mean = min(actual_mean, HEURISTIC_MAX_CONFIDENCE)

    ss = estimate["secondary_structure"]
    result = {
        "protein": protein,
        "mean_plddt": round(actual_mean, 2),
        "plddt_scores": plddt_scores,
        "method": "heuristic_fallback",
        "model_name": "heuristic_v1",
        "confidence": estimate["confidence"],
        "heuristic_details": estimate["heuristic_details"],
        "secondary_structure": {
            "helix_fraction": ss.helix_fraction,
            "sheet_fraction": ss.sheet_fraction,
            "coil_fraction": ss.coil_fraction,
            "assignments": ss.assignments,
        },
    }

    logger.info(
        "Heuristic fallback prediction for %d-aa protein: "
        "estimated pLDDT=%.1f, confidence=%.2f, method=%s",
        n, actual_mean, estimate["confidence"], "heuristic_fallback",
    )

    return result
