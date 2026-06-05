"""
BioCompiler ESMFold Offline Fallback — Sequence-Based Heuristic Structure Prediction
=====================================================================================

Provides a lightweight, offline-capable fallback when the ESMFold API and
local ``esm`` package are both unavailable.  Instead of returning
``success=False`` (which causes structure predicates to return UNCERTAIN),
this module uses sequence-based heuristics to produce a
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

3. **Chou-Fasman secondary structure prediction** (full algorithm):
   - Per-residue helix, sheet, turn, coil propensity scores from
     the published Chou & Fasman (1974) parameters.
   - Helix nucleation (window of 6 with ≥4 residues above P_helix > 1.00)
     and bidirectional extension.
   - Sheet nucleation (window of 5 with ≥3 residues above P_sheet > 1.00)
     and bidirectional extension.
   - Turn prediction using per-position turn probability (4-residue windows)
     with published turn propensity parameters.
   - Proline / Glycine helix-breaking rules.
   - Overlap resolution: helix > sheet > turn > coil.
   - Outputs a DSSP-style secondary structure string (H/E/C).
   - Expected Q3 accuracy: ~50-60% on standard benchmarks.

4. **Contact density estimation**:
   - Uses predicted secondary structure to estimate residue contact density.
   - Helices have higher contact density (~1.8 contacts/residue).
   - Sheets have moderate contact density (~1.5 contacts/residue).
   - Coils have lower contact density (~0.8 contacts/residue).
   - Contact density modulates per-residue pLDDT estimates.

Estimated pLDDT Score
---------------------
The heuristic pLDDT is bounded by ``HEURISTIC_MAX_CONFIDENCE`` (default
55.0, still below ESMFold's "Low confidence" band of 50-70).  Per-residue
pLDDT scores are calibrated based on:

    - Predicted secondary structure type:
      - Helix residues: pLDDT 45–55
      - Sheet residues: pLDDT 40–50
      - Coil residues: pLDDT 25–35
    - Sequence length (shorter = slightly higher confidence)
    - Hydrophobic/hydrophilic balance (modulates within SS band)
    - Contact density (higher density → higher within-band pLDDT)
    - Charge distribution regularity

Overall mean pLDDT is typically 35–50 depending on SS content,
never exceeding HEURISTIC_MAX_CONFIDENCE (55).

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
- Chou & Fasman, Biochemistry 1974; 13:211–222 (turn prediction)
- Lin et al., Science 2023; 379:1043 (ESMFold / ESM-2)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "predict_structure_heuristic",
    "estimate_plddt_from_sequence",
    "estimate_secondary_structure_from_sequence",
    "estimate_fold_quality",
    "compute_hydrophobicity_profile",
    "compute_charge_profile",
    "compute_contact_density",
    "HEURISTIC_MAX_CONFIDENCE",
    "HEURISTIC_MIN_CONFIDENCE",
    "ChargeProfile",
    "SecondaryStructureEstimate",
    "ContactDensityProfile",
    "FoldQualityEstimate",
]


# ==============================================================================
# Constants
# ==============================================================================

#: Maximum pLDDT score the heuristic fallback will ever return.
#: Still below the ESMFold "Low confidence" band of 50-70, but allows
#: per-residue variation based on SS prediction quality.
HEURISTIC_MAX_CONFIDENCE: float = 55.0

#: Minimum pLDDT score the heuristic fallback will ever return.
HEURISTIC_MIN_CONFIDENCE: float = 25.0

#: Kyte-Doolittle hydrophobicity scale (Kyte & Doolittle, 1982).
#: More positive = more hydrophobic.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

#: Chou-Fasman secondary structure propensity parameters.
#: Each amino acid has (P_helix, P_sheet, P_turn).
#: P_helix and P_sheet from Chou & Fasman (1974) Table VI.
#: P_turn values from Chou & Fasman (1974) Table IX (position-independent
#: average of the 4 position-specific turn propensities).
CHOU_FASMAN_PROPENSITY: dict[str, tuple[float, float, float]] = {
    "A": (1.42, 0.83, 0.660),   # Ala: strong helix former
    "R": (0.98, 0.93, 0.950),   # Arg: indifferent helix, moderate turn
    "N": (0.67, 0.89, 1.560),   # Asn: helix breaker, strong turn former
    "D": (1.01, 0.54, 1.460),   # Asp: weak helix, strong turn former
    "C": (0.70, 1.19, 1.190),   # Cys: indifferent, moderate turn
    "Q": (1.11, 1.10, 0.980),   # Gln: moderate helix former
    "E": (1.51, 0.37, 0.740),   # Glu: strong helix former, weak turn
    "G": (0.57, 0.75, 1.560),   # Gly: helix breaker, strong turn former
    "H": (1.00, 0.87, 0.950),   # His: indifferent helix
    "I": (1.08, 1.60, 0.470),   # Ile: moderate helix, strong sheet
    "L": (1.21, 1.30, 0.590),   # Leu: strong helix former
    "K": (1.16, 0.74, 1.010),   # Lys: moderate helix former
    "M": (1.45, 1.05, 0.600),   # Met: strong helix former
    "F": (1.13, 1.38, 0.600),   # Phe: moderate helix/sheet former
    "P": (0.57, 0.55, 1.520),   # Pro: helix breaker, strong turn former
    "S": (0.77, 0.75, 1.430),   # Ser: indifferent, strong turn former
    "T": (0.83, 1.19, 0.960),   # Thr: indifferent, moderate turn
    "W": (1.08, 1.37, 0.960),   # Trp: moderate helix/sheet
    "Y": (0.69, 1.47, 1.140),   # Tyr: weak helix, strong sheet
    "V": (1.06, 1.70, 0.500),   # Val: indifferent helix, strong sheet
}

#: Per-position turn propensities from Chou & Fasman (1974) Table IX.
#: Each amino acid has (P_turn_pos1, P_turn_pos2, P_turn_pos3, P_turn_pos4).
#: Position i, i+1, i+2, i+3 in a 4-residue turn window.
TURN_PROPENSITY_BY_POSITION: dict[str, tuple[float, float, float, float]] = {
    "A": (0.660, 1.190, 0.860, 0.570),
    "R": (1.410, 0.890, 1.010, 0.770),
    "N": (1.560, 2.060, 1.780, 0.950),
    "D": (1.460, 1.010, 2.220, 1.010),
    "C": (0.870, 1.350, 1.310, 0.940),
    "Q": (1.110, 0.730, 1.260, 0.900),
    "E": (0.740, 1.350, 0.840, 0.640),
    "G": (1.560, 0.610, 2.640, 1.630),
    "H": (1.130, 0.940, 1.080, 0.700),
    "I": (0.810, 0.300, 0.660, 0.280),
    "L": (0.980, 0.310, 0.710, 0.410),
    "K": (1.010, 0.810, 1.160, 0.780),
    "M": (0.890, 0.350, 0.640, 0.420),
    "F": (0.840, 0.430, 0.830, 0.530),
    "P": (1.520, 1.590, 1.280, 0.630),
    "S": (1.310, 1.270, 1.690, 1.210),
    "T": (0.840, 0.910, 1.070, 0.990),
    "W": (0.860, 0.530, 1.050, 1.040),
    "Y": (0.780, 1.170, 1.360, 0.910),
    "V": (0.680, 0.240, 0.710, 0.280),
}

#: Positively charged amino acids (Lys, Arg, His).
POSITIVELY_CHARGED_AAS: set[str] = set("KRH")

#: Negatively charged amino acids (Asp, Glu).
NEGATIVELY_CHARGED_AAS: set[str] = set("DE")

#: Hydrophobic amino acids (Ala, Ile, Leu, Met, Phe, Val, Trp).
HYDROPHOBIC_AAS: set[str] = set("AILMFVW")

#: Helix-breaking residues (Proline always breaks; Glycine is a weak breaker)
HELIX_BREAKERS: set[str] = set("PG")

#: Amino acids that promote intrinsic disorder (Pro, Gly, Ser, Gln, Glu, Lys).
DISORDER_PROMOTING_AAS: set[str] = set("PGSQEK")

#: Threshold for Chou-Fasman helix nucleation (a window of 6 residues
#: with 4 or more above P_helix > 1.00 initiates a helix).
HELIX_NUCLEATION_WINDOW: int = 6
HELIX_NUCLEATION_THRESHOLD: int = 4

#: Threshold for Chou-Fasman sheet nucleation (a window of 5 residues
#: with 3 or more above P_sheet > 1.00 initiates a sheet).
SHEET_NUCLEATION_WINDOW: int = 5
SHEET_NUCLEATION_THRESHOLD: int = 3

#: Turn nucleation: a 4-residue window is a turn candidate if the
#: sum of position-specific turn propensities (f(i)*f(i+1)*f(i+2)*f(i+3))
#: exceeds TURN_PRODUCT_THRESHOLD and the average P_turn > TURN_AVG_THRESHOLD.
TURN_WINDOW: int = 4
TURN_PRODUCT_THRESHOLD: float = 0.000075  # Chou-Fasman turn threshold
TURN_AVG_THRESHOLD: float = 1.00

#: Sliding window size for hydrophobicity smoothing.
HYDROPHOBICITY_WINDOW: int = 9

#: Contact density estimates by secondary structure type.
#: Based on average contacts per residue from known protein structures.
CONTACT_DENSITY: dict[str, float] = {
    "H": 1.80,  # Alpha-helix: ~3.6 residues/turn, i→i+3, i→i+4 contacts
    "E": 1.50,  # Beta-sheet: H-bonds between strands, moderate contacts
    "T": 1.20,  # Turn: fewer regular contacts
    "C": 0.80,  # Coil: few regular contacts, more solvent-exposed
}

#: Per-residue pLDDT ranges by secondary structure type.
#: These ranges are below ESMFold's "Low confidence" band (50-70).
PLDDT_RANGES: dict[str, tuple[float, float]] = {
    "H": (45.0, 55.0),  # Helix residues: well-ordered, higher pLDDT
    "E": (40.0, 50.0),  # Sheet residues: moderately ordered
    "T": (30.0, 40.0),  # Turn residues: some order but less stable
    "C": (25.0, 35.0),  # Coil residues: disordered, lowest pLDDT
}

# ---- pLDDT calibration parameters ----

#: Sequence length bounds for length-factor bonus.
#: Proteins shorter than LENGTH_SHORT_THRESHOLD get full bonus;
#: proteins longer than LENGTH_LONG_THRESHOLD get zero bonus.
LENGTH_SHORT_THRESHOLD: int = 100
LENGTH_LONG_THRESHOLD: int = 400

#: Epsilon for floating-point comparisons in charge balance.
_CHARGE_BALANCE_EPSILON: float = 1e-9

#: Default charge balance when no charged residues are present.
_NO_CHARGE_BALANCE_DEFAULT: float = 0.5

#: Contact density smoothing window size (residues).
_CONTACT_DENSITY_SMOOTHING_WINDOW: int = 3

#: Default contact density for unknown SS types.
_DEFAULT_CONTACT_DENSITY: float = 0.8

#: Default charge patch detection window size.
DEFAULT_CHARGE_PATCH_WINDOW: int = 7

#: Default minimum same-sign charge fraction to count as a patch.
DEFAULT_CHARGE_PATCH_MIN_FRAC: float = 0.7

#: Chou-Fasman propensity threshold above which an AA is considered a
#: former for that SS type.
PROPENSITY_THRESHOLD: float = 1.00

#: Minimum helix length (residues) after extension; shorter are pruned to coil.
MIN_HELIX_LENGTH: int = 4

#: Minimum sheet length (residues) after extension; shorter are pruned to coil.
MIN_SHEET_LENGTH: int = 3

#: Per-residue pLDDT length bonus scale factor.
_LENGTH_BONUS_SCALE: float = 2.0

#: Per-residue pLDDT charge balance bonus scale factor.
_CHARGE_BONUS_SCALE: float = 2.0

#: Contact density modulation weight (fraction of half-range).
_CONTACT_DENSITY_MODULATION_WEIGHT: float = 0.6

#: Hydrophobicity modulation weight for structured residues (H, E).
_HYDRO_MODULATION_STRUCTURED_WEIGHT: float = 0.4

#: Hydrophobicity modulation weight for coil residues.
_HYDRO_MODULATION_COIL_WEIGHT: float = 0.2

#: Global bonus scaling applied to per-residue length and charge bonuses.
_GLOBAL_BONUS_SCALE: float = 0.5

#: Hydrophobic balance bonus for structured residues (pLDDT points).
_HYDRO_BALANCE_BONUS: float = 1.5

#: Disorder penalty thresholds and factors (consecutive disorder run length → factor).
_DISORDER_PENALTY_THRESHOLDS: list[tuple[int, float]] = [
    (30, 0.75),
    (20, 0.85),
    (15, 0.90),
    (10, 0.95),
]

#: Charge clustering penalty thresholds and factors (patch count → factor).
_CHARGE_PATCH_PENALTY_THRESHOLDS: list[tuple[int, float]] = [
    (3, 0.85),
    (2, 0.90),
    (1, 0.95),
]

#: Low-complexity penalty thresholds and factors (max repeat → factor).
_LOW_COMPLEXITY_PENALTY_THRESHOLDS: list[tuple[int, float]] = [
    (10, 0.80),
    (6, 0.90),
]

#: Base confidence for heuristic estimates.
_HEURISTIC_BASE_CONFIDENCE: float = 0.15

#: Maximum confidence for heuristic estimates.
_HEURISTIC_MAX_CONFIDENCE: float = 0.5

#: Minimum confidence for heuristic estimates.
_HEURISTIC_MIN_CONFIDENCE_VALUE: float = 0.1


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
    raw: list[float] = []
    for aa in protein:
        if aa not in KYTE_DOOLITTLE:
            logger.warning("Unknown amino acid '%s' in hydrophobicity profile; defaulting to 0.0", aa)
        raw.append(KYTE_DOOLITTLE.get(aa, 0.0))

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


def compute_charge_profile(
    protein: str,
    window: int = DEFAULT_CHARGE_PATCH_WINDOW,
    min_frac: float = DEFAULT_CHARGE_PATCH_MIN_FRAC,
) -> ChargeProfile:
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
    if total_charged > _CHARGE_BALANCE_EPSILON:
        charge_balance = 1.0 - abs(pos_frac - neg_frac) / total_charged
    else:
        # No charged residues at all — neutral, but not informative
        charge_balance = _NO_CHARGE_BALANCE_DEFAULT

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
# Full Chou-Fasman secondary structure prediction
# ==============================================================================

@dataclass
class SecondaryStructureEstimate:
    """Estimated secondary structure from the full Chou-Fasman algorithm.

    Attributes:
        helix_fraction:   Fraction of residues predicted to be helical.
        sheet_fraction:   Fraction of residues predicted to be beta-sheet.
        turn_fraction:    Fraction of residues predicted to be turn.
        coil_fraction:    Fraction of residues predicted to be coil.
        assignments:      Per-residue SS assignment ('H', 'E', 'T', 'C').
        ss_string:        DSSP-style secondary structure string (H/E/C).
                          Turn regions are mapped to 'C' in the DSSP string
                          since DSSP uses H/E/C convention.
    """
    helix_fraction: float
    sheet_fraction: float
    turn_fraction: float
    coil_fraction: float
    assignments: list[str]
    ss_string: str


def _compute_turn_probabilities(protein: str) -> list[float]:
    """Compute per-residue turn probability from position-specific propensities.

    For each position i, computes the geometric mean of position-specific
    turn propensities over all 4-residue windows that include position i.
    This gives a smoothed turn propensity per residue.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        List of turn probabilities, one per residue.
    """
    protein = protein.upper()
    n = len(protein)

    if n < TURN_WINDOW:
        return [0.0] * n

    # Compute turn probability for each 4-residue window
    # f(i,i+1,i+2,i+3) = f_t(i) * f_t(i+1) * f_t(i+2) * f_t(i+3)
    # where f_t is the position-specific turn propensity
    window_turn_probs: list[float] = []
    for i in range(n - TURN_WINDOW + 1):
        product = 1.0
        avg_turn = 0.0
        for j in range(TURN_WINDOW):
            aa = protein[i + j]
            pos_props = TURN_PROPENSITY_BY_POSITION.get(aa, (1.0, 1.0, 1.0, 1.0))
            if aa not in TURN_PROPENSITY_BY_POSITION:
                logger.debug("Unknown AA '%s' in turn propensity lookup; using defaults", aa)
            product *= pos_props[j]
            avg_turn += CHOU_FASMAN_PROPENSITY.get(aa, (1.0, 1.0, 1.0))[2]

        avg_turn /= TURN_WINDOW

        # A turn is predicted if product > threshold AND avg P_turn > 1.0
        if product > TURN_PRODUCT_THRESHOLD and avg_turn > TURN_AVG_THRESHOLD:
            window_turn_probs.append(product)
        else:
            window_turn_probs.append(0.0)

    # Convert window-level turn probs to per-residue by averaging
    # over all windows containing each residue
    per_residue: list[float] = [0.0] * n
    for i in range(n):
        count: int = 0
        total: float = 0.0
        # Which windows include position i?
        # Window starting at j covers positions j, j+1, j+2, j+3
        # So j ranges from max(0, i-3) to min(n-4, i)
        for j in range(max(0, i - TURN_WINDOW + 1), min(n - TURN_WINDOW + 1, i + 1)):
            total += window_turn_probs[j]
            count += 1
        per_residue[i] = total / count if count > 0 else 0.0

    return per_residue


def estimate_secondary_structure_from_sequence(protein: str) -> SecondaryStructureEstimate:
    """Estimate secondary structure content using the full Chou-Fasman algorithm.

    The Chou-Fasman algorithm (Chou & Fasman, 1974) predicts helix, sheet,
    turn, and coil regions from amino acid sequence alone.  This implementation
    includes:

    1. **Helix nucleation and extension**: A window of 6 residues with ≥4
       above P_helix > 1.00 nucleates a helix.  The helix is extended in
       both directions while the average P_helix remains > 1.00.

    2. **Sheet nucleation and extension**: A window of 5 residues with ≥3
       above P_sheet > 1.00 nucleates a sheet.  The sheet is extended in
       both directions while the average P_sheet remains > 1.00.

    3. **Turn prediction**: 4-residue windows with high position-specific
       turn propensity products (f(i)*f(i+1)*f(i+2)*f(i+3) > threshold)
       and average P_turn > 1.00 are predicted as turns.

    4. **Proline/Glycine helix-breaking**: Proline breaks helices (except
       at position 1).  Glycine destabilizes helices.

    5. **Overlap resolution**: When regions overlap, priority is:
       helix > sheet > turn > coil.

    6. **Short helix/sheet pruning**: Helices shorter than 4 residues and
       sheets shorter than 3 residues are demoted to coil.

    Expected Q3 accuracy: ~50-60% on standard benchmarks.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        SecondaryStructureEstimate with fractions, per-residue assignments,
        and DSSP-style secondary structure string.
    """
    protein = protein.upper()
    n = len(protein)

    if n == 0:
        return SecondaryStructureEstimate(
            helix_fraction=0.0, sheet_fraction=0.0,
            turn_fraction=0.0, coil_fraction=1.0,
            assignments=[], ss_string="",
        )

    # Per-residue propensities
    p_helix: list[float] = []
    p_sheet: list[float] = []
    for aa in protein:
        if aa not in CHOU_FASMAN_PROPENSITY:
            logger.warning("Unknown amino acid '%s' in SS prediction; using default propensity", aa)
        props = CHOU_FASMAN_PROPENSITY.get(aa, (1.0, 1.0, 1.0))
        p_helix.append(props[0])
        p_sheet.append(props[1])

    # Assignments: default to coil
    assignments: list[str] = ["C"] * n

    # --- Step 1: Turn prediction ---
    # Identify turn regions using position-specific propensities
    turn_regions: list[tuple[int, int]] = []
    for i in range(n - TURN_WINDOW + 1):
        window = protein[i:i + TURN_WINDOW]
        product = 1.0
        avg_turn = 0.0
        for j in range(TURN_WINDOW):
            aa = window[j]
            pos_props = TURN_PROPENSITY_BY_POSITION.get(aa, (1.0, 1.0, 1.0, 1.0))
            if aa not in TURN_PROPENSITY_BY_POSITION:
                logger.debug("Unknown AA '%s' in turn propensity lookup; using defaults", aa)
            product *= pos_props[j]
            avg_turn += CHOU_FASMAN_PROPENSITY.get(aa, (1.0, 1.0, 1.0))[2]

        avg_turn /= TURN_WINDOW

        if product > TURN_PRODUCT_THRESHOLD and avg_turn > TURN_AVG_THRESHOLD:
            turn_regions.append((i, i + TURN_WINDOW))

    # --- Step 2: Helix nucleation and extension ---
    helix_regions: list[tuple[int, int]] = []
    i = 0
    while i <= n - HELIX_NUCLEATION_WINDOW:
        window_helix = p_helix[i:i + HELIX_NUCLEATION_WINDOW]
        above_threshold = sum(1 for p in window_helix if p > PROPENSITY_THRESHOLD)
        if above_threshold >= HELIX_NUCLEATION_THRESHOLD:
            # Nucleate helix, extend in both directions
            start = i
            end = i + HELIX_NUCLEATION_WINDOW

            # Extend N-terminal
            while start > 0:
                # Check for helix breakers at the extension point
                if protein[start - 1] == "P":
                    break
                avg_helix = sum(p_helix[max(0, start - 1):end]) / (end - start + 1)
                if avg_helix > PROPENSITY_THRESHOLD:
                    start -= 1
                else:
                    break

            # Extend C-terminal
            while end < n:
                # Proline can only be at position 1 of a helix
                if protein[end] == "P":
                    break
                avg_helix = sum(p_helix[start:min(n, end + 1)]) / (end - start + 1)
                if avg_helix > PROPENSITY_THRESHOLD:
                    end += 1
                else:
                    break

            # Prune short helices (< MIN_HELIX_LENGTH residues after extension)
            if end - start >= MIN_HELIX_LENGTH:
                helix_regions.append((start, end))
            i = end
        else:
            i += 1

    # --- Step 3: Sheet nucleation and extension ---
    sheet_regions: list[tuple[int, int]] = []
    i = 0
    while i <= n - SHEET_NUCLEATION_WINDOW:
        window_sheet = p_sheet[i:i + SHEET_NUCLEATION_WINDOW]
        above_threshold = sum(1 for p in window_sheet if p > PROPENSITY_THRESHOLD)
        if above_threshold >= SHEET_NUCLEATION_THRESHOLD:
            start = i
            end = i + SHEET_NUCLEATION_WINDOW

            # Extend N-terminal
            while start > 0:
                avg_sheet = sum(p_sheet[max(0, start - 1):end]) / (end - start + 1)
                if avg_sheet > PROPENSITY_THRESHOLD:
                    start -= 1
                else:
                    break

            # Extend C-terminal
            while end < n:
                avg_sheet = sum(p_sheet[start:min(n, end + 1)]) / (end - start + 1)
                if avg_sheet > PROPENSITY_THRESHOLD:
                    end += 1
                else:
                    break

            # Prune short sheets (< MIN_SHEET_LENGTH residues after extension)
            if end - start >= MIN_SHEET_LENGTH:
                sheet_regions.append((start, end))
            i = end
        else:
            i += 1

    # --- Step 4: Merge overlapping helix regions ---
    if helix_regions:
        merged_helix: list[tuple[int, int]] = [helix_regions[0]]
        for start, end in helix_regions[1:]:
            prev_start, prev_end = merged_helix[-1]
            if start <= prev_end:
                # Overlapping or adjacent — merge
                merged_helix[-1] = (prev_start, max(prev_end, end))
            else:
                merged_helix.append((start, end))
        helix_regions = merged_helix

    # --- Step 5: Merge overlapping sheet regions ---
    if sheet_regions:
        merged_sheet: list[tuple[int, int]] = [sheet_regions[0]]
        for start, end in sheet_regions[1:]:
            prev_start, prev_end = merged_sheet[-1]
            if start <= prev_end:
                merged_sheet[-1] = (prev_start, max(prev_end, end))
            else:
                merged_sheet.append((start, end))
        sheet_regions = merged_sheet

    # --- Step 6: Resolve overlaps — helix > sheet > turn > coil ---
    # Apply helix assignments first (highest priority)
    for start, end in helix_regions:
        for j in range(start, min(end, n)):
            assignments[j] = "H"

    # Apply sheet assignments (only where not already helix)
    for start, end in sheet_regions:
        for j in range(start, min(end, n)):
            if assignments[j] in ("C", "T"):
                assignments[j] = "E"

    # Apply turn assignments (only where still coil)
    for start, end in turn_regions:
        for j in range(start, min(end, n)):
            if assignments[j] == "C":
                assignments[j] = "T"

    # --- Step 7: Re-check helix breaking after overlap resolution ---
    # Proline inside a helix (not position 0) should break it
    for i_pos in range(n):
        if assignments[i_pos] == "H" and protein[i_pos] == "P":
            # Proline can be at position 0 of a helix but not internal
            # Find the start of this helix region
            helix_start = i_pos
            while helix_start > 0 and assignments[helix_start - 1] == "H":
                helix_start -= 1
            # If Pro is at helix_start, it's OK (N-cap position)
            if i_pos > helix_start:
                # Break the helix at proline: everything from proline onward
                # becomes coil/turn
                assignments[i_pos] = "C"

    # Compute fractions
    helix_count = sum(1 for a in assignments if a == "H")
    sheet_count = sum(1 for a in assignments if a == "E")
    turn_count = sum(1 for a in assignments if a == "T")
    coil_count = n - helix_count - sheet_count - turn_count

    # Generate DSSP-style string: H=helix, E=sheet, C=coil/turn
    # DSSP convention: H and E are the main SS types; turns map to C
    ss_string = "".join(
        "H" if a == "H" else "E" if a == "E" else "C"
        for a in assignments
    )

    return SecondaryStructureEstimate(
        helix_fraction=round(helix_count / n, 4) if n > 0 else 0.0,
        sheet_fraction=round(sheet_count / n, 4) if n > 0 else 0.0,
        turn_fraction=round(turn_count / n, 4) if n > 0 else 0.0,
        coil_fraction=round(coil_count / n, 4) if n > 0 else 1.0,
        assignments=assignments,
        ss_string=ss_string,
    )


# ==============================================================================
# Contact density estimation
# ==============================================================================

@dataclass
class ContactDensityProfile:
    """Per-residue contact density estimate based on predicted SS.

    Attributes:
        per_residue:  Contact density estimate for each residue.
        mean:         Mean contact density across all residues.
        ss_weighted:  SS-content-weighted average contact density.
    """
    per_residue: list[float]
    mean: float
    ss_weighted: float


def compute_contact_density(ss_assignments: list[str]) -> ContactDensityProfile:
    """Estimate per-residue contact density from predicted secondary structure.

    Contact density is estimated based on the predicted secondary structure:
      - Helix residues: ~1.8 contacts/residue (i→i+3, i→i+4 H-bonds)
      - Sheet residues: ~1.5 contacts/residue (inter-strand H-bonds)
      - Turn residues:  ~1.2 contacts/residue (some local contacts)
      - Coil residues:  ~0.8 contacts/residue (mostly solvent-exposed)

    A sliding window smoothing (window=3) is applied to reduce noise
    at SS boundaries.

    Args:
        ss_assignments: Per-residue SS assignments ('H', 'E', 'T', 'C').

    Returns:
        ContactDensityProfile with per-residue, mean, and SS-weighted densities.
    """
    n = len(ss_assignments)
    if n == 0:
        return ContactDensityProfile(per_residue=[], mean=0.0, ss_weighted=0.0)

    # Raw contact density per residue
    raw = [CONTACT_DENSITY.get(ss, _DEFAULT_CONTACT_DENSITY) for ss in ss_assignments]

    # Smooth with a sliding window to reduce boundary effects
    half_w = _CONTACT_DENSITY_SMOOTHING_WINDOW // 2
    smoothed: list[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        avg = sum(raw[start:end]) / (end - start)
        smoothed.append(round(avg, 3))

    mean_density = sum(smoothed) / n if n > 0 else 0.0

    # SS-weighted density: weight by SS type
    helix_count = sum(1 for a in ss_assignments if a == "H")
    sheet_count = sum(1 for a in ss_assignments if a == "E")
    turn_count = sum(1 for a in ss_assignments if a == "T")
    coil_count = n - helix_count - sheet_count - turn_count

    ss_weighted = (
        (helix_count * CONTACT_DENSITY["H"] +
         sheet_count * CONTACT_DENSITY["E"] +
         turn_count * CONTACT_DENSITY["T"] +
         coil_count * CONTACT_DENSITY["C"]) / n if n > 0 else 0.0
    )

    return ContactDensityProfile(
        per_residue=smoothed,
        mean=round(mean_density, 4),
        ss_weighted=round(ss_weighted, 4),
    )


# ==============================================================================
# Fold quality estimation (hydrophobic burial / polar surface)
# ==============================================================================

@dataclass
class FoldQualityEstimate:
    """Estimate of how well a sequence's hydrophobicity pattern supports folding.

    A well-folded globular protein tends to have:
      - Hydrophobic residues buried in the core (interior of the sequence
        or clustered together)
      - Polar / charged residues on the surface (termini or interspersed)

    This estimate uses the Kyte-Doolittle hydrophobicity profile to assess
    whether the sequence's hydrophobicity pattern is consistent with a
    compact, well-folded structure.

    Attributes:
        hydro_burial_score:  0–1 score indicating how well hydrophobic residues
                             are concentrated in the sequence interior (0 = poor,
                             1 = ideal burial pattern).
        polar_surface_score: 0–1 score indicating how well polar residues are
                             positioned at sequence boundaries / surface (0 = poor,
                             1 = ideal surface distribution).
        hydro_core_detected: Whether a hydrophobic core region is detected
                             (consecutive stretch of above-average hydrophobicity).
        overall_quality:     0–1 composite score (weighted average of burial
                             and surface scores).
        interpretation:      Human-readable interpretation of the quality score.
    """

    hydro_burial_score: float
    polar_surface_score: float
    hydro_core_detected: bool
    overall_quality: float
    interpretation: str


#: Window size for detecting hydrophobic core regions (residues).
_HYDRO_CORE_WINDOW: int = 11

#: Minimum average hydrophobicity in a window to count as "core" region.
#: Positive Kyte-Doolittle values indicate hydrophobicity; this threshold
#: means the window average must be above the scale midpoint.
_HYDRO_CORE_THRESHOLD: float = 0.5

#: Weight for hydrophobic burial score in overall quality.
_BURIAL_WEIGHT: float = 0.6

#: Weight for polar surface score in overall quality.
_SURFACE_WEIGHT: float = 0.4


def estimate_fold_quality(protein: str) -> FoldQualityEstimate:
    """Estimate structure fold quality from hydrophobicity burial patterns.

    Uses the Kyte-Doolittle hydrophobicity profile to assess whether the
    sequence's hydrophobicity pattern is consistent with a compact,
    well-folded globular protein:

    1. **Hydrophobic burial**: In well-folded proteins, hydrophobic residues
       cluster in the core.  We check if the smoothed hydrophobicity profile
       has elevated values in the sequence interior (middle portion) compared
       to the termini.

    2. **Polar surface**: In well-folded proteins, polar/charged residues
       are enriched at the surface.  We check if the N- and C-terminal
       regions have lower hydrophobicity (more polar) than the interior.

    3. **Hydrophobic core detection**: We scan for contiguous stretches
       where the smoothed hydrophobicity exceeds a threshold, indicating
       a potential buried core.

    **Length-aware adjustments**: Small proteins (10–50 residues) may not
    have a traditional hydrophobic core with clear interior/termini
    distinctions, yet they can still be well-folded (e.g., zinc fingers,
    mini-proteins, peptide hormones).  For these proteins, the burial
    score is de-emphasised in favour of overall hydrophobicity and
    secondary-structure content, which are stronger indicators of
    compact folding at short lengths.

    Returns a :class:`FoldQualityEstimate` with scores in 0–1 range.
    Higher scores indicate a hydrophobicity pattern more consistent with
    a well-folded globular protein.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        FoldQualityEstimate with burial, surface, and overall quality scores.
    """
    protein = protein.upper()
    n = len(protein)

    if n < 10:
        return FoldQualityEstimate(
            hydro_burial_score=0.0,
            polar_surface_score=0.0,
            hydro_core_detected=False,
            overall_quality=0.0,
            interpretation="Sequence too short for fold quality assessment",
        )

    hydro_profile = compute_hydrophobicity_profile(protein)
    mean_hydro = sum(hydro_profile) / len(hydro_profile)

    # --- Determine protein size category ---
    # Small proteins (10-50 aa) may lack a traditional hydrophobic core;
    # the interior/termini distinction is less meaningful.  We apply
    # length-aware weighting to the burial and surface scores.
    _SMALL_PROTEIN_THRESHOLD: int = 50
    is_small = n < _SMALL_PROTEIN_THRESHOLD

    # --- 1. Hydrophobic burial score ---
    # In a well-folded globular protein, the interior of the sequence
    # (middle 60%) should be more hydrophobic than the terminal regions
    # (first and last 20%).
    inner_start = max(1, n // 5)
    inner_end = min(n - 1, 4 * n // 5)
    outer_start = 0
    outer_end = inner_start
    outer2_start = inner_end
    outer2_end = n

    inner_hydro = sum(hydro_profile[inner_start:inner_end]) / max(1, inner_end - inner_start)
    outer1_hydro = sum(hydro_profile[outer_start:outer_end]) / max(1, outer_end - outer_start)
    outer2_hydro = sum(hydro_profile[outer2_start:outer2_end]) / max(1, outer2_end - outer2_start)
    outer_hydro = (outer1_hydro + outer2_hydro) / 2.0

    # Burial score: how much more hydrophobic is the interior vs the termini?
    # Normalize: typical difference is 0.5–2.0 on the KD scale for well-folded
    hydro_diff = inner_hydro - outer_hydro
    burial_score = max(0.0, min(1.0, (hydro_diff + 1.0) / 3.0))

    # --- 1b. Small-protein hydrophobicity score ---
    # For short proteins, the burial score is unreliable because there is
    # little distinction between "interior" and "surface" in the sequence.
    # Instead, we assess overall hydrophobicity balance: a mix of
    # hydrophobic and polar residues is consistent with a compact fold,
    # even in the absence of a classic core.
    if is_small:
        hydrophobic_count = sum(1 for aa in protein if aa in HYDROPHOBIC_AAS)
        hydro_frac = hydrophobic_count / n
        # Optimal hydrophobic fraction for small folded proteins: 0.25–0.50
        if 0.25 <= hydro_frac <= 0.50:
            small_hydro_score = 0.7
        elif 0.15 <= hydro_frac < 0.25 or 0.50 < hydro_frac <= 0.60:
            small_hydro_score = 0.5
        else:
            small_hydro_score = 0.2

        # Also consider overall mean hydrophobicity: positive values
        # suggest a hydrophobic tendency consistent with folding.
        if mean_hydro > 0:
            small_hydro_score = min(1.0, small_hydro_score + 0.1)

        # Blend the burial score with the small-protein score, giving
        # more weight to the small-protein heuristic for short sequences.
        # The shorter the protein, the more we trust the heuristic.
        small_weight = max(0.0, min(0.7, (_SMALL_PROTEIN_THRESHOLD - n) / _SMALL_PROTEIN_THRESHOLD))
        burial_score = (1.0 - small_weight) * burial_score + small_weight * small_hydro_score

    # --- 2. Polar surface score ---
    # Terminal regions should be more polar (lower KD score) than interior
    # If outer_hydro < inner_hydro, that's good (polar termini)
    if mean_hydro != 0:
        polarity_ratio = 1.0 - (outer_hydro - mean_hydro) / (abs(mean_hydro) + 1.0)
    else:
        polarity_ratio = 0.5

    # Also check: what fraction of charged residues are in the terminal 20%?
    terminal_count = sum(1 for i, aa in enumerate(protein)
                         if (i < inner_start or i >= inner_end)
                         and aa in POSITIVELY_CHARGED_AAS | NEGATIVELY_CHARGED_AAS)
    terminal_fraction = (inner_start + (n - inner_end)) / n
    expected_terminal_charged = terminal_fraction  # if random
    actual_terminal_charged_frac = terminal_count / max(1, n) / max(0.01, terminal_fraction)

    # Surface score: polar termini + charged residues enriched at termini
    surface_score = max(0.0, min(1.0,
        0.5 * max(0.0, polarity_ratio) +
        0.5 * min(1.0, actual_terminal_charged_frac)
    ))

    # For small proteins, the terminal enrichment signal is weak because
    # the "terminal" and "interior" regions overlap.  We supplement the
    # surface score with overall charge balance — balanced charges are
    # consistent with a well-folded small protein.
    if is_small:
        charge_prof = compute_charge_profile(protein)
        # Blend: replace part of surface_score with charge balance
        small_surface_weight = max(0.0, min(0.5, (_SMALL_PROTEIN_THRESHOLD - n) / _SMALL_PROTEIN_THRESHOLD))
        surface_score = (1.0 - small_surface_weight) * surface_score + small_surface_weight * charge_prof.charge_balance

    # --- 3. Hydrophobic core detection ---
    # Scan for windows with above-average hydrophobicity.
    # For small proteins, use a proportionally smaller window.
    hydro_core_detected = False
    effective_core_window = min(_HYDRO_CORE_WINDOW, max(5, n // 2))
    if n >= effective_core_window:
        half_w = effective_core_window // 2
        for i in range(half_w, n - half_w):
            window_avg = sum(hydro_profile[i - half_w:i + half_w + 1]) / effective_core_window
            if window_avg > _HYDRO_CORE_THRESHOLD:
                hydro_core_detected = True
                break
    else:
        # For very short sequences, just check if overall hydrophobicity is positive
        hydro_core_detected = mean_hydro > _HYDRO_CORE_THRESHOLD

    # --- 4. Overall quality ---
    # For small proteins, reduce the burial weight and increase the surface
    # weight, since the burial score is less reliable.
    if is_small:
        burial_weight = 0.4
        surface_weight = 0.6
    else:
        burial_weight = _BURIAL_WEIGHT
        surface_weight = _SURFACE_WEIGHT

    overall = burial_weight * burial_score + surface_weight * surface_score

    # Bonus for detected hydrophobic core
    if hydro_core_detected:
        overall = min(1.0, overall + 0.05)

    # --- 5. Interpretation (length-aware) ---
    if is_small:
        if overall >= 0.6:
            interpretation = "Good (small protein): hydrophobicity balance consistent with compact fold"
        elif overall >= 0.4:
            interpretation = "Moderate (small protein): partial hydrophobic character, may fold or be partially structured"
        elif overall >= 0.2:
            interpretation = "Weak (small protein): limited hydrophobic content, may be disordered or peptide-like"
        else:
            interpretation = "Poor (small protein): hydrophobicity pattern inconsistent with compact folding"
    else:
        if overall >= 0.7:
            interpretation = "Good: hydrophobic burial pattern consistent with folded globular protein"
        elif overall >= 0.5:
            interpretation = "Moderate: partial hydrophobic core, some polar surface enrichment"
        elif overall >= 0.3:
            interpretation = "Weak: limited hydrophobic burial, may be partially disordered"
        else:
            interpretation = "Poor: hydrophobicity pattern inconsistent with compact folding"

    return FoldQualityEstimate(
        hydro_burial_score=round(burial_score, 4),
        polar_surface_score=round(surface_score, 4),
        hydro_core_detected=hydro_core_detected,
        overall_quality=round(overall, 4),
        interpretation=interpretation,
    )


# ==============================================================================
# Heuristic pLDDT estimation (improved with SS-based calibration)
# ==============================================================================

def _compute_per_residue_plddt(
    ss_assignments: list[str],
    hydro_profile: list[float],
    contact_density: list[float],
    sequence_length: int,
    charge_balance: float,
    hydrophobic_fraction: float,
) -> list[float]:
    """Compute per-residue pLDDT scores based on SS prediction and calibration.

    Each residue's pLDDT is set within the range for its SS type:
      - Helix (H): 45–55
      - Sheet (E): 40–50
      - Turn (T):  30–40
      - Coil (C):  25–35

    Within each range, the score is modulated by:
      - Contact density (higher density → higher within-range score)
      - Hydrophobicity (more hydrophobic → slight boost for structured regions)
      - Sequence length factor (shorter proteins slightly higher)
      - Charge balance (better balance → slight boost)

    The final score is clamped to [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE].

    Args:
        ss_assignments: Per-residue SS assignments ('H', 'E', 'T', 'C').
        hydro_profile: Smoothed Kyte-Doolittle hydrophobicity profile.
        contact_density: Per-residue contact density estimates.
        sequence_length: Length of the protein sequence.
        charge_balance: Charge balance score (0-1).
        hydrophobic_fraction: Fraction of hydrophobic residues.

    Returns:
        List of per-residue pLDDT scores.
    """
    n = len(ss_assignments)
    if n == 0:
        return []

    # Sequence length factor: shorter proteins get a slight boost
    # Proteins < 100 aa: factor = 1.0; proteins > 400 aa: factor = 0.0
    length_factor = max(0.0, min(1.0, (LENGTH_LONG_THRESHOLD - sequence_length) / (LENGTH_LONG_THRESHOLD - LENGTH_SHORT_THRESHOLD)))
    length_bonus = length_factor * _LENGTH_BONUS_SCALE  # 0–2 point bonus for shorter proteins

    # Charge balance bonus: 0–2 points
    charge_bonus = charge_balance * _CHARGE_BONUS_SCALE

    # Hydrophobic balance factor: optimal is 0.30–0.50
    if 0.30 <= hydrophobic_fraction <= 0.50:
        hydro_balance_factor = 1.0
    elif 0.20 <= hydrophobic_fraction < 0.30 or 0.50 < hydrophobic_fraction <= 0.60:
        hydro_balance_factor = 0.5
    else:
        hydro_balance_factor = 0.0

    # Normalize hydrophobicity profile for modulation
    if hydro_profile:
        min_hydro = min(hydro_profile)
        max_hydro = max(hydro_profile)
        hydro_range = max_hydro - min_hydro if max_hydro > min_hydro else 1.0
    else:
        hydro_range = 1.0
        min_hydro = 0.0

    # Normalize contact density for modulation
    if contact_density:
        min_cd = min(contact_density)
        max_cd = max(contact_density)
        cd_range = max_cd - min_cd if max_cd > min_cd else 1.0
    else:
        cd_range = 1.0
        min_cd = 0.0

    plddt_scores: list[float] = []
    for i in range(n):
        ss_type = ss_assignments[i]
        low, high = PLDDT_RANGES.get(ss_type, (25.0, 35.0))

        # Start at midpoint of the SS-type range
        base_score = (low + high) / 2.0
        half_range = (high - low) / 2.0

        # Modulation factor: how far from midpoint (0 = midpoint, -1 = low, +1 = high)
        modulation = 0.0

        # Contact density modulation (0–0.3 of half_range)
        if i < len(contact_density) and cd_range > 0:
            cd_norm = (contact_density[i] - min_cd) / cd_range  # 0–1
            modulation += (cd_norm - 0.5) * _CONTACT_DENSITY_MODULATION_WEIGHT  # -0.3 to +0.3

        # Hydrophobicity modulation (0–0.2 of half_range)
        # More hydrophobic → slight boost for structured regions
        if i < len(hydro_profile) and hydro_range > 0:
            hydro_norm = (hydro_profile[i] - min_hydro) / hydro_range  # 0–1
            if ss_type in ("H", "E"):
                modulation += (hydro_norm - 0.5) * _HYDRO_MODULATION_STRUCTURED_WEIGHT  # -0.2 to +0.2
            elif ss_type == "C":
                # For coil, hydrophobic residues are less well-predicted
                # (they might be buried but we can't know)
                modulation -= abs(hydro_norm - 0.5) * _HYDRO_MODULATION_COIL_WEIGHT

        # Apply modulation within the SS-type range
        residue_score = base_score + modulation * half_range

        # Add global bonuses
        residue_score += length_bonus * _GLOBAL_BONUS_SCALE  # Scale down since it's per-residue
        residue_score += charge_bonus * _GLOBAL_BONUS_SCALE

        # Apply hydrophobic balance bonus for structured residues
        if ss_type in ("H", "E"):
            residue_score += hydro_balance_factor * _HYDRO_BALANCE_BONUS

        # Clamp to [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE]
        residue_score = max(HEURISTIC_MIN_CONFIDENCE, min(residue_score, HEURISTIC_MAX_CONFIDENCE))
        plddt_scores.append(round(residue_score, 2))

    return plddt_scores


def estimate_plddt_from_sequence(protein: str) -> dict[str, Any]:
    """Estimate per-residue pLDDT scores from protein sequence alone.

    This function uses the full Chou-Fasman secondary structure prediction,
    contact density estimation, and calibration factors to produce per-residue
    pLDDT estimates.  The estimates are bounded by
    [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE] and should never
    be treated as reliable predictions.

    Scoring approach:
        - Per-residue pLDDT is primarily determined by predicted SS type:
          - Helix residues: 45–55
          - Sheet residues: 40–50
          - Turn residues: 30–40
          - Coil residues: 25–35
        - Within each range, scores are modulated by:
          - Contact density (higher → higher within range)
          - Hydrophobicity (more hydrophobic → slight boost for structured)
          - Sequence length (shorter → slight boost)
          - Charge balance (better balance → slight boost)
        - Overall mean is typically 35–50 depending on SS content.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        Dictionary with keys:
            - ``estimated_mean_plddt`` (float): Mean of per-residue pLDDT.
            - ``confidence`` (float): Confidence in this estimate (0–0.5).
            - ``method`` (str): Always ``"heuristic_fallback"``.
            - ``heuristic_details`` (dict): Breakdown of calibration factors.
            - ``hydrophobicity_profile`` (list[float]): Per-residue smoothed
              Kyte-Doolittle values.
            - ``charge_profile`` (ChargeProfile): Charge distribution summary.
            - ``secondary_structure`` (SecondaryStructureEstimate): SS estimate.
            - ``contact_density`` (ContactDensityProfile): Contact density estimate.
            - ``ss_prediction`` (str): DSSP-style secondary structure string.
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
                helix_fraction=0.0, sheet_fraction=0.0,
                turn_fraction=0.0, coil_fraction=1.0,
                assignments=[], ss_string="",
            ),
            "contact_density": ContactDensityProfile(
                per_residue=[], mean=0.0, ss_weighted=0.0,
            ),
            "ss_prediction": "",
        }

    # --- Compute sub-profiles ---
    hydro_profile = compute_hydrophobicity_profile(protein)
    charge_prof = compute_charge_profile(protein)
    ss_estimate = estimate_secondary_structure_from_sequence(protein)
    contact_density = compute_contact_density(ss_estimate.assignments)

    # --- Compute hydrophobic fraction ---
    hydrophobic_count = sum(1 for aa in protein if aa in HYDROPHOBIC_AAS)
    hydro_frac = hydrophobic_count / n

    # --- Compute disorder signals ---
    disorder_promoting = DISORDER_PROMOTING_AAS
    consecutive_disorder = 0
    max_consecutive_disorder = 0
    for aa in protein:
        if aa in disorder_promoting:
            consecutive_disorder += 1
            max_consecutive_disorder = max(max_consecutive_disorder, consecutive_disorder)
        else:
            consecutive_disorder = 0

    # Low-complexity: long runs of same amino acid
    max_repeat = 1
    current_repeat = 1
    for i in range(1, n):
        if protein[i] == protein[i - 1]:
            current_repeat += 1
            max_repeat = max(max_repeat, current_repeat)
        else:
            current_repeat = 1

    # --- Compute per-residue pLDDT ---
    plddt_scores = _compute_per_residue_plddt(
        ss_assignments=ss_estimate.assignments,
        hydro_profile=hydro_profile,
        contact_density=contact_density.per_residue,
        sequence_length=n,
        charge_balance=charge_prof.charge_balance,
        hydrophobic_fraction=hydro_frac,
    )

    # Apply disorder penalty: reduce pLDDT for proteins with long disorder runs
    disorder_penalty_factor = 1.0
    for threshold, factor in _DISORDER_PENALTY_THRESHOLDS:
        if max_consecutive_disorder > threshold:
            disorder_penalty_factor = factor
            break

    # Apply charge clustering penalty
    patch_count = charge_prof.charge_patch_count
    charge_penalty_factor = 1.0
    for threshold, factor in _CHARGE_PATCH_PENALTY_THRESHOLDS:
        if patch_count >= threshold:
            charge_penalty_factor = factor
            break

    # Apply low-complexity penalty
    lc_penalty_factor = 1.0
    for threshold, factor in _LOW_COMPLEXITY_PENALTY_THRESHOLDS:
        if max_repeat > threshold:
            lc_penalty_factor = factor
            break

    # Combined penalty (multiplicative)
    combined_penalty = disorder_penalty_factor * charge_penalty_factor * lc_penalty_factor

    # Apply penalty to per-residue scores
    plddt_scores = [
        max(HEURISTIC_MIN_CONFIDENCE,
            min(round(score * combined_penalty, 2), HEURISTIC_MAX_CONFIDENCE))
        for score in plddt_scores
    ]

    # Compute mean
    mean_plddt = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0
    mean_plddt = max(HEURISTIC_MIN_CONFIDENCE, min(mean_plddt, HEURISTIC_MAX_CONFIDENCE))

    # --- Confidence in this estimate ---
    # Confidence is inherently low for heuristics. We scale it based on
    # SS content (more defined SS → slightly higher confidence) and
    # signal quality (more regular → slightly higher confidence).
    defined_ss_frac = ss_estimate.helix_fraction + ss_estimate.sheet_fraction

    # Base confidence: 0.15
    # SS content bonus: 0–0.15 (more defined SS → higher)
    # Charge balance bonus: 0–0.05
    # Disorder penalty: 0–0.10
    # Contact density bonus: 0–0.05
    ss_bonus = min(0.15, defined_ss_frac * 0.25)
    cb_bonus = charge_prof.charge_balance * 0.05
    disorder_pen = min(0.10, max_consecutive_disorder * 0.003) if max_consecutive_disorder > 5 else 0.0
    cd_bonus = min(0.05, (contact_density.ss_weighted - 0.8) * 0.1) if contact_density.ss_weighted > 0.8 else 0.0

    confidence = _HEURISTIC_BASE_CONFIDENCE + ss_bonus + cb_bonus + cd_bonus - disorder_pen
    confidence = min(_HEURISTIC_MAX_CONFIDENCE, max(_HEURISTIC_MIN_CONFIDENCE_VALUE, confidence))

    # Mean hydrophobicity
    mean_hydro = sum(hydro_profile) / len(hydro_profile) if hydro_profile else 0.0

    heuristic_details = {
        "per_residue_method": "ss_based_calibration",
        "helix_plddt_range": PLDDT_RANGES["H"],
        "sheet_plddt_range": PLDDT_RANGES["E"],
        "turn_plddt_range": PLDDT_RANGES["T"],
        "coil_plddt_range": PLDDT_RANGES["C"],
        "hydrophobic_fraction": round(hydro_frac, 4),
        "mean_hydrophobicity": round(mean_hydro, 3),
        "defined_ss_fraction": round(defined_ss_frac, 4),
        "helix_fraction": ss_estimate.helix_fraction,
        "sheet_fraction": ss_estimate.sheet_fraction,
        "turn_fraction": ss_estimate.turn_fraction,
        "charge_balance": charge_prof.charge_balance,
        "contact_density_mean": contact_density.mean,
        "contact_density_ss_weighted": contact_density.ss_weighted,
        "disorder_penalty_factor": round(disorder_penalty_factor, 4),
        "charge_clustering_penalty_factor": round(charge_penalty_factor, 4),
        "low_complexity_penalty_factor": round(lc_penalty_factor, 4),
        "combined_penalty_factor": round(combined_penalty, 4),
        "max_consecutive_disorder": max_consecutive_disorder,
        "max_repeat": max_repeat,
        "charge_patch_count": patch_count,
        "sequence_length_factor": round(max(0.0, min(1.0, (LENGTH_LONG_THRESHOLD - n) / (LENGTH_LONG_THRESHOLD - LENGTH_SHORT_THRESHOLD))), 4),
    }

    return {
        "estimated_mean_plddt": round(mean_plddt, 2),
        "confidence": round(confidence, 4),
        "method": "heuristic_fallback",
        "heuristic_details": heuristic_details,
        "hydrophobicity_profile": hydro_profile,
        "charge_profile": charge_prof,
        "secondary_structure": ss_estimate,
        "contact_density": contact_density,
        "ss_prediction": ss_estimate.ss_string,
        "plddt_scores": plddt_scores,
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
      - ``mean_plddt`` in range [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE]
      - ``confidence < 0.5``
      - Per-residue pLDDT scores calibrated by SS prediction

    **Accuracy warning**: This heuristic is NOT a substitute for ESMFold.
    It is based on biophysical rules (Chou-Fasman SS prediction, contact
    density estimation, hydrophobicity analysis) and can be inaccurate
    for any specific protein.  It exists solely to provide a non-UNCERTAIN
    signal when no real structure prediction is available, so that
    downstream predicates can produce a tentative verdict rather than
    defaulting to UNCERTAIN.

    **Improvements over v1**:
      - Full Chou-Fasman algorithm with turn prediction and Pro/Gly breaking
      - Contact density estimation from predicted SS
      - Per-residue pLDDT calibrated by SS type (helix 45-55, sheet 40-50,
        coil 25-35) instead of uniform modulation around a single mean
      - pLDDT range expanded to 30-55 (from 0-40) for better differentiation
      - DSSP-style secondary structure string output
      - Model name updated to "heuristic_v2"

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        Dictionary with keys needed to build an ESMFoldResult:
            - ``protein`` (str): Input sequence.
            - ``mean_plddt`` (float): Estimated mean pLDDT.
            - ``plddt_scores`` (list[float]): Per-residue estimated pLDDT.
            - ``method`` (str): ``"heuristic_fallback"``.
            - ``model_name`` (str): ``"heuristic_v2"``.
            - ``confidence`` (float): Confidence in the prediction (0–0.5).
            - ``heuristic_details`` (dict): Scoring breakdown.
            - ``secondary_structure`` (dict): Estimated SS fractions.
            - ``ss_prediction`` (str): DSSP-style SS string.
    """
    protein = protein.upper()
    n = len(protein)

    estimate = estimate_plddt_from_sequence(protein)
    plddt_scores = estimate["plddt_scores"]

    # Recompute mean from per-residue scores for consistency
    actual_mean = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0
    actual_mean = max(HEURISTIC_MIN_CONFIDENCE, min(actual_mean, HEURISTIC_MAX_CONFIDENCE))

    ss = estimate["secondary_structure"]
    # --- Fold quality from hydrophobic burial / polar surface analysis ---
    fold_quality = estimate_fold_quality(protein)

    # Modulate pLDDT by fold quality: if the hydrophobicity pattern is
    # inconsistent with folding, reduce the heuristic pLDDT estimates.
    if fold_quality.overall_quality < 0.3:
        fold_modulation = 0.85
    elif fold_quality.overall_quality < 0.5:
        fold_modulation = 0.92
    elif fold_quality.overall_quality >= 0.7:
        fold_modulation = 1.0
    else:
        fold_modulation = 0.97

    if fold_modulation < 1.0:
        plddt_scores = [
            max(HEURISTIC_MIN_CONFIDENCE,
                min(round(s * fold_modulation, 2), HEURISTIC_MAX_CONFIDENCE))
            for s in plddt_scores
        ]
        actual_mean = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0
        actual_mean = max(HEURISTIC_MIN_CONFIDENCE, min(actual_mean, HEURISTIC_MAX_CONFIDENCE))

    result = {
        "protein": protein,
        "mean_plddt": round(actual_mean, 2),
        "plddt_scores": plddt_scores,
        "method": "heuristic_fallback",
        "model_name": "heuristic_v2",
        "confidence": estimate["confidence"],
        "heuristic_details": estimate["heuristic_details"],
        "secondary_structure": {
            "helix_fraction": ss.helix_fraction,
            "sheet_fraction": ss.sheet_fraction,
            "turn_fraction": ss.turn_fraction,
            "coil_fraction": ss.coil_fraction,
            "assignments": ss.assignments,
        },
        "ss_prediction": estimate["ss_prediction"],
        "fold_quality": fold_quality,
    }

    logger.info(
        "Heuristic fallback prediction for %d-aa protein: "
        "estimated pLDDT=%.1f, confidence=%.2f, method=%s, "
        "SS: H=%.2f E=%.2f T=%.2f C=%.2f, "
        "fold_quality=%.2f (%s)",
        n, actual_mean, estimate["confidence"], "heuristic_fallback",
        ss.helix_fraction, ss.sheet_fraction,
        ss.turn_fraction, ss.coil_fraction,
        fold_quality.overall_quality, fold_quality.interpretation,
    )

    return result
