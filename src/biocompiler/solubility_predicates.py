"""
BioCompiler Solubility Predicates v7.2.0
=========================================
Type system predicates for protein solubility assessment.

Uses CamSol-style intrinsic solubility scoring and biophysical heuristics
to predict whether a protein will express in soluble form or aggregate.

Predicates:
  - SolubleExpression: protein predicted soluble via CamSol scoring
  - NoAggregationProneRegion: no long aggregation-prone regions
  - ChargeComposition: adequate charged residue fraction and safe pI
  - NoLongHydrophobicStretch: no excessively long hydrophobic stretches

Helper functions:
  - compute_approximate_pI: estimate isoelectric point via bisection
  - compute_net_charge: net charge at a given pH (Henderson-Hasselbalch)
  - find_hydrophobic_stretches: detect maximal consecutive hydrophobic runs
"""

from __future__ import annotations

import logging
from typing import Optional

from biocompiler.type_system import Verdict, TypeCheckResult

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# pKa values for ionizable groups (charge calculation)
# ────────────────────────────────────────────────────────────
PKA_VALUES: dict[str, float] = {
    "N_term": 9.69,   # N-terminal amino group
    "C_term": 2.34,   # C-terminal carboxyl group
    "K": 10.54,       # Lysine side chain
    "R": 12.48,       # Arginine side chain
    "H": 6.04,        # Histidine side chain
    "D": 3.90,        # Aspartic acid side chain
    "E": 4.07,        # Glutamic acid side chain
}

# ────────────────────────────────────────────────────────────
# CamSol intrinsic solubility profile (per-residue)
# Based on amino acid physicochemical properties.
# Positive values → soluble propensity; negative → aggregation-prone.
# Scaled so that overall score × 3.0 aligns with CamSol thresholds.
# ────────────────────────────────────────────────────────────
_CAMSOL_INTRINSIC: dict[str, float] = {
    "A":  0.10, "R":  1.30, "N":  0.60, "D":  1.20, "C": -0.50,
    "Q":  0.50, "E":  1.40, "G":  0.00, "H":  0.30, "I": -1.50,
    "L": -1.30, "K":  1.50, "M": -0.80, "F": -1.60, "P":  0.20,
    "S":  0.50, "T":  0.10, "W": -1.70, "Y": -0.90, "V": -1.20,
}

# Scaling factor: converts mean smoothed profile to CamSol-score range
# where > 1.5 = highly soluble, < -1.0 = insoluble.
_CAMSOL_SCALE = 3.0

# Default smoothing window for CamSol profile
_CAMSOL_WINDOW = 7

# Default hydrophobic residue set (AILMFWV)
_DEFAULT_HYDROPHOBIC: set[str] = set("AILMFWV")


# ────────────────────────────────────────────────────────────
# Internal CamSol helpers
# ────────────────────────────────────────────────────────────

def _camsol_smoothed_profile(protein: str, window: int = _CAMSOL_WINDOW) -> list[float]:
    """Compute per-residue CamSol smoothed solubility profile.

    Applies a sliding-window average over the raw intrinsic scores,
    then scales by ``_CAMSOL_SCALE`` so the values align with CamSol
    threshold conventions (overall > 1.5 = highly soluble).

    Args:
        protein: Upper-cased amino-acid sequence.
        window: Smoothing window width (odd number preferred).

    Returns:
        List of smoothed, scaled per-residue scores.  Length == len(protein).
    """
    n = len(protein)
    if n == 0:
        return []

    raw = [_CAMSOL_INTRINSIC.get(aa, 0.0) for aa in protein]
    half_w = window // 2

    smoothed: list[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        avg = sum(raw[start:end]) / (end - start)
        smoothed.append(avg * _CAMSOL_SCALE)

    return smoothed


def _camsol_overall_score(protein: str, window: int = _CAMSOL_WINDOW) -> float:
    """Compute overall CamSol intrinsic solubility score.

    Args:
        protein: Upper-cased amino-acid sequence.
        window: Smoothing window width.

    Returns:
        Overall CamSol score (typical range -3 to +3).
        > 1.5 highly soluble, 0–1.5 soluble, -1–0 marginal, < -1 insoluble.
    """
    profile = _camsol_smoothed_profile(protein, window)
    if not profile:
        return 0.0
    return sum(profile) / len(profile)


def _find_aggregation_regions(
    protein: str,
    window: int = _CAMSOL_WINDOW,
    score_threshold: float = -1.0,
) -> list[tuple[int, int, float]]:
    """Find aggregation-prone regions from CamSol per-residue scores.

    A region is a maximal consecutive run of residues whose smoothed
    CamSol score is below *score_threshold*.

    Args:
        protein: Upper-cased amino-acid sequence.
        window: Smoothing window width.
        score_threshold: Per-residue score below which a residue is
            considered aggregation-prone.

    Returns:
        List of (start, end, avg_score) tuples (end is exclusive).
    """
    profile = _camsol_smoothed_profile(protein, window)
    if not profile:
        return []

    regions: list[tuple[int, int, float]] = []
    i = 0
    n = len(profile)
    while i < n:
        if profile[i] < score_threshold:
            start = i
            while i < n and profile[i] < score_threshold:
                i += 1
            end = i
            avg_score = sum(profile[start:end]) / (end - start)
            regions.append((start, end, round(avg_score, 3)))
        else:
            i += 1

    return regions


# ────────────────────────────────────────────────────────────
# Public helper functions
# ────────────────────────────────────────────────────────────

def compute_net_charge(protein: str, pH: float) -> float:
    """Compute the net charge of a protein at a given pH.

    Uses the Henderson-Hasselbalch equation for each ionizable group.

    Positive contributions (protonated at low pH):
        N-terminus, K, R, H  →  count × 1 / (1 + 10^(pH − pKa))

    Negative contributions (deprotonated at high pH):
        C-terminus, D, E     →  −count × 1 / (1 + 10^(pKa − pH))

    Args:
        protein: Amino-acid sequence (single-letter codes).
        pH: pH value at which to compute net charge.

    Returns:
        Net charge (float).  Positive = basic; negative = acidic.
    """
    if not protein:
        return 0.0

    protein = protein.upper()
    charge = 0.0

    # Positive groups
    # N-terminal amino group
    charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["N_term"]))

    for aa in protein:
        if aa == "K":
            charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["K"]))
        elif aa == "R":
            charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["R"]))
        elif aa == "H":
            charge += 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["H"]))

    # Negative groups
    # C-terminal carboxyl group
    charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C_term"] - pH))

    for aa in protein:
        if aa == "D":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["D"] - pH))
        elif aa == "E":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["E"] - pH))

    return charge


def compute_approximate_pI(protein: str) -> float:
    """Compute approximate isoelectric point (pI) of a protein.

    The pI is the pH at which the net charge is zero.  Uses bisection
    over [0, 14] with 100 iterations (~1e-30 precision) to find the
    crossover.

    Args:
        protein: Amino-acid sequence (single-letter codes).

    Returns:
        Approximate pI value (float, in range [0, 14]).
        Returns 7.0 for empty input.
    """
    if not protein:
        return 7.0

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        charge = compute_net_charge(protein, mid)
        if charge > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def find_hydrophobic_stretches(
    protein: str,
    hydrophobic: set[str] | None = None,
) -> list[tuple[int, int]]:
    """Find all maximal consecutive hydrophobic stretches.

    A stretch is a maximal run of residues where every residue belongs
    to the *hydrophobic* set.

    Args:
        protein: Amino-acid sequence (single-letter codes).
        hydrophobic: Set of single-letter hydrophobic residue codes.
            Defaults to ``{'A', 'I', 'L', 'M', 'F', 'W', 'V'}``.

    Returns:
        List of (start, end) tuples (end exclusive) for every maximal
        hydrophobic stretch of length ≥ 1.
    """
    if not protein:
        return []

    protein = protein.upper()
    hydro = hydrophobic if hydrophobic is not None else _DEFAULT_HYDROPHOBIC

    stretches: list[tuple[int, int]] = []
    start: int | None = None

    for i, aa in enumerate(protein):
        if aa in hydro:
            if start is None:
                start = i
        else:
            if start is not None:
                stretches.append((start, i))
                start = None

    # Handle stretch extending to the C-terminus
    if start is not None:
        stretches.append((start, len(protein)))

    return stretches


# ────────────────────────────────────────────────────────────
# Predicate 1: Soluble Expression
# ────────────────────────────────────────────────────────────

def evaluate_soluble_expression(
    sequence: str,
    protein: str,
    organism: str,
    min_solubility_score: float = 0.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check if protein is predicted to be soluble using CamSol scoring.

    Verdict logic (CamSol intrinsic score):
    - Score > 1.5        → PASS  (highly soluble)
    - Score 0.0 to 1.5   → LIKELY_PASS (soluble)
    - Score -1.0 to 0.0  → UNCERTAIN (marginally soluble)
    - Score < -1.0       → LIKELY_FAIL (insoluble)

    If *pdb_string* is provided, a note is added that structural
    correction would improve accuracy.

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        min_solubility_score: Minimum acceptable CamSol score (default 0.0).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and CamSol score derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="SolubleExpression",
            verdict=Verdict.FAIL,
            violation="Empty protein sequence",
        )

    camsol_score = _camsol_overall_score(protein)

    # Identify aggregation-prone regions for derivation
    agg_regions = _find_aggregation_regions(protein)

    # Determine verdict
    if camsol_score > 1.5:
        verdict = Verdict.PASS
        violation = None
    elif camsol_score >= 0.0:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif camsol_score >= -1.0:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Marginal solubility: CamSol score {camsol_score:.3f} "
            f"is in the uncertain range [-1.0, 0.0)"
        )
    else:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Insoluble protein: CamSol score {camsol_score:.3f} < -1.0"
        )

    # If the score is below the user-specified minimum, that's a stronger
    # failure signal
    if camsol_score < min_solubility_score and verdict in (
        Verdict.PASS, Verdict.LIKELY_PASS,
    ):
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"CamSol score {camsol_score:.3f} is below the minimum "
            f"acceptable score {min_solubility_score}"
        )

    # Build derivation
    derivation: list[dict] = [
        {"step": "camsol_intrinsic_score", "value": round(camsol_score, 3)},
        {"step": "min_solubility_score", "value": min_solubility_score},
    ]
    if agg_regions:
        derivation.append({
            "step": "aggregation_prone_regions",
            "value": [
                {"start": s, "end": e, "avg_score": sc}
                for s, e, sc in agg_regions
            ],
        })
    else:
        derivation.append({
            "step": "aggregation_prone_regions",
            "value": [],
        })

    # Knowledge gap when no PDB structure
    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "No PDB structure provided; solubility estimated from "
            "intrinsic sequence properties only.  Structural correction "
            "would improve accuracy."
        )

    return TypeCheckResult(
        predicate=f"SolubleExpression(min={min_solubility_score})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 2: No Aggregation-Prone Region
# ────────────────────────────────────────────────────────────

def evaluate_no_aggregation_prone_region(
    sequence: str,
    protein: str,
    organism: str,
    max_region_length: int = 5,
    score_threshold: float = -1.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check for aggregation-prone regions using CamSol per-residue scoring.

    A region is a consecutive run of residues whose smoothed CamSol score
    is below *score_threshold*.  Verdict depends on the longest region:

    - No region longer than *max_region_length* → PASS
    - Region length max_region_length+1 to 7    → LIKELY_PASS (borderline)
    - Region length 8 to 10                     → UNCERTAIN
    - Region length 11 to 15                    → LIKELY_FAIL
    - Region length > 15                        → FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        max_region_length: Maximum acceptable region length (default 5).
        score_threshold: Per-residue CamSol score threshold (default -1.0).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and aggregation-prone region details.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="NoAggregationProneRegion",
            verdict=Verdict.PASS,
        )

    agg_regions = _find_aggregation_regions(protein, score_threshold=score_threshold)

    if not agg_regions:
        return TypeCheckResult(
            predicate=f"NoAggregationProneRegion(max={max_region_length}, "
                      f"threshold={score_threshold})",
            verdict=Verdict.PASS,
            derivation=[
                {"step": "aggregation_prone_regions", "value": []},
                {"step": "longest_region", "value": 0},
            ],
        )

    # Find the longest region
    longest = max(end - start for start, end, _ in agg_regions)

    # Determine verdict based on longest region length
    if longest <= max_region_length:
        verdict = Verdict.PASS
        violation = None
    elif longest <= 7:
        verdict = Verdict.LIKELY_PASS
        violation = (
            f"Borderline aggregation-prone region of {longest} residues "
            f"(max allowed: {max_region_length})"
        )
    elif longest <= 10:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Aggregation-prone region of {longest} residues detected "
            f"(max allowed: {max_region_length})"
        )
    elif longest <= 15:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Long aggregation-prone region of {longest} residues "
            f"(max allowed: {max_region_length})"
        )
    else:
        verdict = Verdict.FAIL
        violation = (
            f"Very long aggregation-prone region of {longest} residues "
            f"(max allowed: {max_region_length})"
        )

    # Build derivation: list all aggregation-prone regions
    derivation: list[dict] = [
        {
            "step": "aggregation_prone_regions",
            "value": [
                {"start": s, "end": e, "length": e - s, "avg_score": sc}
                for s, e, sc in agg_regions
            ],
        },
        {"step": "longest_region", "value": longest},
        {"step": "max_region_length", "value": max_region_length},
        {"step": "score_threshold", "value": score_threshold},
    ]

    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "Aggregation-prone regions identified from intrinsic "
            "sequence only.  Structural accessibility correction "
            "would refine the prediction."
        )

    return TypeCheckResult(
        predicate=f"NoAggregationProneRegion(max={max_region_length}, "
                  f"threshold={score_threshold})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 3: Charge Composition
# ────────────────────────────────────────────────────────────

def evaluate_charge_composition(
    sequence: str,
    protein: str,
    organism: str,
    min_charged_fraction: float = 0.10,
    max_pI: float = 9.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check charge composition for solubility.

    Computes two metrics:
    1. **Charged fraction** — fraction of K, R, H, D, E residues.
       Below *min_charged_fraction* → LIKELY_FAIL (too few charges
       for solubility).
    2. **Isoelectric point (pI)** — pH at which net charge is zero.
       Above *max_pI* → UNCERTAIN (protein may precipitate near its pI
       in typical buffers).

    Both OK → PASS.

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        min_charged_fraction: Minimum fraction of charged residues (default 0.10).
        max_pI: Maximum acceptable isoelectric point (default 9.0).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict, charged fraction, pI, and residue
        counts in the derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="ChargeComposition",
            verdict=Verdict.FAIL,
            violation="Empty protein sequence",
        )

    n = len(protein)

    # Count charged residues
    charged_set = {"K", "R", "H", "D", "E"}
    pos_count = sum(1 for aa in protein if aa in {"K", "R", "H"})
    neg_count = sum(1 for aa in protein if aa in {"D", "E"})
    charged_count = pos_count + neg_count
    charged_fraction = charged_count / n

    # Compute isoelectric point
    pI = compute_approximate_pI(protein)

    # Evaluate both conditions
    low_charge = charged_fraction < min_charged_fraction
    high_pI = pI > max_pI

    if low_charge and high_pI:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Low charged fraction ({charged_fraction:.1%} < {min_charged_fraction:.0%}) "
            f"and high pI ({pI:.2f} > {max_pI})"
        )
    elif low_charge:
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"Low charged fraction ({charged_fraction:.1%} < {min_charged_fraction:.0%}): "
            f"insufficient surface charges for solubility"
        )
    elif high_pI:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"High isoelectric point (pI={pI:.2f} > {max_pI}): "
            f"protein may precipitate near its pI in typical buffers"
        )
    else:
        verdict = Verdict.PASS
        violation = None

    derivation: list[dict] = [
        {"step": "charged_fraction", "value": round(charged_fraction, 4)},
        {"step": "min_charged_fraction", "value": min_charged_fraction},
        {"step": "isoelectric_point", "value": round(pI, 2)},
        {"step": "max_pI", "value": max_pI},
        {"step": "positive_residues", "value": pos_count},
        {"step": "negative_residues", "value": neg_count},
        {"step": "total_charged", "value": charged_count},
        {"step": "protein_length", "value": n},
    ]

    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "Charge composition assessed from sequence alone.  Surface "
            "accessibility of charged residues (from structure) would "
            "improve solubility prediction."
        )

    return TypeCheckResult(
        predicate=f"ChargeComposition(min_charged={min_charged_fraction}, "
                  f"max_pI={max_pI})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )


# ────────────────────────────────────────────────────────────
# Predicate 4: No Long Hydrophobic Stretch
# ────────────────────────────────────────────────────────────

def evaluate_no_long_hydrophobic_stretch(
    sequence: str,
    protein: str,
    organism: str,
    max_stretch: int = 7,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Check for long consecutive hydrophobic stretches (AILMFWV).

    Stretches longer than *max_stretch* are aggregation-prone.

    Verdict logic (based on the longest stretch found):
    - No stretch > *max_stretch*                     → PASS
    - Longest stretch *max_stretch*+1 to *max_stretch*+3 → LIKELY_PASS
    - Longest stretch *max_stretch*+4 to *max_stretch*+6 → UNCERTAIN
    - Longest stretch > *max_stretch*+6              → FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino-acid sequence (single-letter codes).
        organism: Target organism name.
        max_stretch: Maximum acceptable hydrophobic stretch length (default 7).
        pdb_string: Optional PDB-format structure string.

    Returns:
        TypeCheckResult with verdict and positions of long hydrophobic
        stretches in the derivation.
    """
    protein = protein.upper()

    if not protein:
        return TypeCheckResult(
            predicate="NoLongHydrophobicStretch",
            verdict=Verdict.PASS,
        )

    stretches = find_hydrophobic_stretches(protein)

    # Filter to stretches longer than max_stretch
    long_stretches = [(s, e) for s, e in stretches if (e - s) > max_stretch]

    if not long_stretches:
        # All stretches are within the allowed length
        max_found = max((e - s) for s, e in stretches) if stretches else 0
        return TypeCheckResult(
            predicate=f"NoLongHydrophobicStretch(max={max_stretch})",
            verdict=Verdict.PASS,
            derivation=[
                {"step": "max_stretch_found", "value": max_found},
                {"step": "max_stretch_allowed", "value": max_stretch},
                {"step": "long_stretches", "value": []},
            ],
        )

    # Find the longest stretch
    longest = max(e - s for s, e in long_stretches)
    excess = longest - max_stretch

    # Determine verdict based on how much the longest stretch exceeds the limit
    if excess <= 3:
        verdict = Verdict.LIKELY_PASS
        violation = (
            f"Hydrophobic stretch of {longest} residues slightly exceeds "
            f"limit of {max_stretch} (borderline)"
        )
    elif excess <= 6:
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Hydrophobic stretch of {longest} residues exceeds limit "
            f"of {max_stretch} by {excess}"
        )
    else:
        verdict = Verdict.FAIL
        violation = (
            f"Very long hydrophobic stretch of {longest} residues exceeds "
            f"limit of {max_stretch} by {excess} (aggregation-prone)"
        )

    # Build derivation: list all long hydrophobic stretches
    derivation: list[dict] = [
        {
            "step": "long_stretches",
            "value": [
                {"start": s, "end": e, "length": e - s}
                for s, e in long_stretches
            ],
        },
        {"step": "longest_stretch", "value": longest},
        {"step": "max_stretch_allowed", "value": max_stretch},
    ]

    knowledge_gap: str | None = None
    if pdb_string is None:
        knowledge_gap = (
            "Hydrophobic stretches assessed from sequence alone.  "
            "Structural context (buried vs. exposed) would clarify "
            "whether stretches contribute to aggregation."
        )

    return TypeCheckResult(
        predicate=f"NoLongHydrophobicStretch(max={max_stretch})",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        knowledge_gap=knowledge_gap,
    )
