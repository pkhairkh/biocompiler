"""
BioCompiler CamSol Solubility Module v7.2.0
=============================================
CamSol-inspired solubility prediction algorithm in pure Python.

Predicts protein solubility from sequence using intrinsic physicochemical
properties (hydropathy, charge, secondary structure propensity) with
optional structure-based corrections from PDB data.

Based on the CamSol method (Sormanni et al., J Mol Biol 2015) for
predicting protein solubility from sequence.

Scoring range: -3 to +3 (positive = soluble, negative = aggregation-prone).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

from .constants import BLOSUM62

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# CamSol-specific physicochemical scales
# ────────────────────────────────────────────────────────────

CAMSOL_HYDROPATHY: dict[str, float] = {
    # CamSol-specific hydropathy scale (optimized for solubility prediction)
    # Based on Wimley-White octanol scale with CamSol corrections
    # Positive = hydrophilic (soluble), Negative = hydrophobic (aggregation-prone)
    "I": -0.31,
    "V": -0.07,
    "L": -0.56,
    "F": -0.52,
    "C":  0.24,
    "M": -0.23,
    "A":  0.17,
    "G":  0.27,
    "T":  0.26,
    "S":  0.40,
    "W": -0.60,
    "Y": -0.18,
    "P":  0.99,
    "H":  0.61,
    "E":  1.32,
    "Q":  0.73,
    "D":  1.26,
    "N":  0.64,
    "K":  1.23,
    "R":  0.60,
}

CAMSOL_CHARGE: dict[str, float] = {
    # Net charge contribution at pH 7.4
    "K":  1.0,
    "R":  1.0,
    "H":  0.5,  # partially protonated at physiological pH
    "D": -1.0,
    "E": -1.0,
    "A":  0.0,
    "V":  0.0,
    "I":  0.0,
    "L":  0.0,
    "F":  0.0,
    "C":  0.0,
    "M":  0.0,
    "G":  0.0,
    "T":  0.0,
    "S":  0.0,
    "W":  0.0,
    "Y":  0.0,
    "P":  0.0,
    "Q":  0.0,
    "N":  0.0,
}

CAMSOL_ALPHA_HELIX: dict[str, float] = {
    # Alpha-helix propensity (solubility effect)
    # Positive values indicate helix-forming residues that reduce aggregation
    "A":  0.23,
    "L":  0.11,
    "E":  0.37,
    "M":  0.06,
    "Q":  0.33,
    "K":  0.26,
    "R":  0.21,
    "H":  0.09,
    "V": -0.08,
    "I": -0.07,
    "W": -0.09,
    "F": -0.11,
    "Y": -0.07,
    "C": -0.03,
    "G": -0.12,
    "P": -0.25,
    "S":  0.03,
    "T":  0.04,
    "D":  0.15,
    "N":  0.05,
}

CAMSOL_BETA_STRAND: dict[str, float] = {
    # Beta-strand propensity (aggregation effect)
    # Positive values indicate strand-forming residues that increase aggregation
    "V":  0.34,
    "I":  0.30,
    "T":  0.14,
    "F":  0.28,
    "Y":  0.13,
    "W":  0.19,
    "L":  0.22,
    "C":  0.10,
    "M":  0.08,
    "Q": -0.04,
    "N": -0.05,
    "S": -0.03,
    "R": -0.06,
    "H": -0.05,
    "E": -0.10,
    "D": -0.09,
    "A": -0.01,
    "G": -0.14,
    "P": -0.28,
    "K": -0.11,
}

# ────────────────────────────────────────────────────────────
# Standard amino acids and internal constants
# ────────────────────────────────────────────────────────────

_STANDARD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

# Weights for each component in the intrinsic solubility calculation
_WEIGHT_HYDROPATHY = 0.35
_WEIGHT_CHARGE = 0.25
_WEIGHT_ALPHA_HELIX = 0.15
_WEIGHT_BETA_STRAND = 0.15
_WEIGHT_PROGLY = 0.10

# Aggregation-prone region threshold
_AGGREGATION_THRESHOLD = -0.5

# Minimum length for an aggregation-prone region
_MIN_AGGREGATION_REGION_LENGTH = 3


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class SolubilityResult:
    """Result of CamSol solubility prediction.

    Attributes:
        protein: Input protein sequence (1-letter codes).
        intrinsic_score: CamSol intrinsic solubility score (-3 to +3, >0 = soluble).
        structural_score: Structure-corrected score (if PDB available), else None.
        overall_score: Final solubility score (structural if available, else intrinsic).
        per_residue_scores: Contribution of each residue to solubility.
        aggregation_prone_regions: List of (start, end, avg_score) tuples for
            regions below the aggregation threshold.
        solubility_class: One of "highly_soluble", "soluble",
            "marginally_soluble", "insoluble".
        recommendations: Actionable suggestions for improving solubility.
        method: "camsol_intrinsic" or "camsol_structural".
    """
    protein: str
    intrinsic_score: float
    structural_score: float | None
    overall_score: float
    per_residue_scores: list[float]
    aggregation_prone_regions: list[tuple[int, int, float]]
    solubility_class: str
    recommendations: list[str]
    method: str


# ────────────────────────────────────────────────────────────
# Core computation functions
# ────────────────────────────────────────────────────────────

def compute_intrinsic_solubility(
    protein: str,
    window: int = 7,
    smoothing: int = 3,
) -> SolubilityResult:
    """Compute CamSol intrinsic solubility score from protein sequence.

    The intrinsic solubility is computed as a weighted combination of
    physicochemical properties: hydropathy, charge, secondary structure
    propensity, and proline/glycine content.

    Steps:
        1. Compute per-residue intrinsic score using weighted combination.
        2. Apply sliding window smoothing (default window=7).
        3. Apply additional smoothing pass (3-residue window).
        4. Compute global score as average of per-residue scores.
        5. Identify aggregation-prone regions.
        6. Classify solubility.
        7. Generate recommendations.

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        window: Sliding window size for initial smoothing (default 7).
        smoothing: Secondary smoothing window size (default 3).

    Returns:
        SolubilityResult with intrinsic solubility prediction.

    Raises:
        ValueError: If protein is empty or contains non-standard residues.
    """
    protein = protein.upper().strip()

    if not protein:
        raise ValueError("Protein sequence must not be empty.")

    invalid = set(protein) - set(_STANDARD_AAS)
    if invalid:
        raise ValueError(
            f"Protein contains non-standard residues: {invalid}. "
            f"Only standard 20 amino acids are supported."
        )

    n = len(protein)
    if n == 0:
        raise ValueError("Protein sequence must not be empty.")

    # ── Step 1: Per-residue intrinsic score ──
    raw_scores = []
    for aa in protein:
        hydro = CAMSOL_HYDROPATHY.get(aa, 0.0)
        charge = CAMSOL_CHARGE.get(aa, 0.0)
        alpha = CAMSOL_ALPHA_HELIX.get(aa, 0.0)
        beta = CAMSOL_BETA_STRAND.get(aa, 0.0)

        # Proline/glycine effect: P increases solubility (breaks structure),
        # G decreases solubility (flexible, can adopt aggregation-prone
        # conformations)
        if aa == "P":
            progly = 1.0
        elif aa == "G":
            progly = -0.5
        else:
            progly = 0.0

        score = (
            _WEIGHT_HYDROPATHY * hydro
            + _WEIGHT_CHARGE * charge
            + _WEIGHT_ALPHA_HELIX * alpha
            + _WEIGHT_BETA_STRAND * (-beta)  # strand increases aggregation
            + _WEIGHT_PROGLY * progly
        )
        raw_scores.append(score)

    # ── Step 2: Sliding window smoothing ──
    smoothed = _sliding_window_smooth(raw_scores, window)

    # ── Step 3: Additional smoothing pass ──
    double_smoothed = _sliding_window_smooth(smoothed, smoothing)

    # ── Step 4: Compute global score ──
    intrinsic_score = (
        sum(double_smoothed) / len(double_smoothed) if double_smoothed else 0.0
    )

    # Clamp to [-3, +3]
    intrinsic_score = max(-3.0, min(3.0, intrinsic_score))

    # ── Step 5: Identify aggregation-prone regions ──
    agg_regions = _find_aggregation_prone_regions(
        double_smoothed, _AGGREGATION_THRESHOLD
    )

    # ── Step 6: Classify solubility ──
    solubility_class = classify_solubility(intrinsic_score)

    # ── Step 7: Generate recommendations ──
    result = SolubilityResult(
        protein=protein,
        intrinsic_score=round(intrinsic_score, 4),
        structural_score=None,
        overall_score=round(intrinsic_score, 4),
        per_residue_scores=[round(s, 4) for s in double_smoothed],
        aggregation_prone_regions=[
            (start, end, round(avg, 4)) for start, end, avg in agg_regions
        ],
        solubility_class=solubility_class,
        recommendations=[],  # filled below
        method="camsol_intrinsic",
    )
    result.recommendations = generate_solubility_recommendations(result)

    logger.debug(
        "CamSol intrinsic: protein=%s, len=%d, score=%.4f, class=%s, "
        "agg_regions=%d",
        protein[:10] + "..." if len(protein) > 10 else protein,
        n,
        intrinsic_score,
        solubility_class,
        len(agg_regions),
    )

    return result


def compute_structural_solubility(
    protein: str,
    pdb_string: str,
) -> SolubilityResult:
    """Compute structure-corrected CamSol solubility score.

    Starts with the intrinsic solubility score and applies structure-based
    corrections derived from approximate SASA calculations on the PDB data:
      - Buried residues (low SASA): reduce their aggregation penalty
        (already buried in the structure).
      - Exposed hydrophobic residues: increase aggregation penalty.
      - Disulfide bonds: increase stability, moderate effect on solubility.

    SASA is approximated from PDB CA coordinates using neighbor counting
    (no external dependency required).

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        pdb_string: PDB file content as a string.

    Returns:
        SolubilityResult with structure-corrected solubility prediction.

    Raises:
        ValueError: If protein is empty or PDB string is empty.
    """
    protein = protein.upper().strip()

    if not protein:
        raise ValueError("Protein sequence must not be empty.")
    if not pdb_string.strip():
        raise ValueError("PDB string must not be empty.")

    # Compute intrinsic first
    intrinsic = compute_intrinsic_solubility(protein)

    # Parse PDB coordinates (CA atoms only)
    ca_coords = _parse_pdb_ca_coords(pdb_string)

    if len(ca_coords) < 3:
        logger.warning(
            "Too few CA atoms in PDB (%d). Falling back to intrinsic score.",
            len(ca_coords),
        )
        intrinsic.method = "camsol_structural"
        intrinsic.structural_score = intrinsic.intrinsic_score
        return intrinsic

    # Approximate SASA using CA neighbor counting
    sasa = _approximate_sasa(ca_coords)

    # Detect disulfide bonds from SSBOND records
    disulfide_residues = _parse_disulfide_bonds(pdb_string)

    # Apply structure-based corrections
    corrected_scores = list(intrinsic.per_residue_scores)

    for i in range(min(len(corrected_scores), len(sasa))):
        residue_sasa = sasa[i]

        if residue_sasa < 0.15:
            # Buried residue: reduce aggregation penalty
            # If score is negative (aggregation-prone), reduce the magnitude
            if corrected_scores[i] < 0:
                correction = abs(corrected_scores[i]) * 0.5 * (
                    1.0 - residue_sasa / 0.15
                )
                corrected_scores[i] += correction
        elif residue_sasa > 0.40:
            # Exposed residue
            hydro = CAMSOL_HYDROPATHY.get(protein[i], 0.0)
            if hydro < 0:
                # Exposed hydrophobic: increase aggregation penalty
                penalty = abs(hydro) * 0.3 * (residue_sasa - 0.40) / 0.60
                corrected_scores[i] -= penalty

        # Disulfide bond correction
        if i in disulfide_residues:
            corrected_scores[i] += 0.15  # moderate stabilizing effect

    # Re-smooth after corrections
    corrected_smooth = _sliding_window_smooth(corrected_scores, 3)

    # Recompute global score
    structural_score = (
        sum(corrected_smooth) / len(corrected_smooth)
        if corrected_smooth
        else 0.0
    )
    structural_score = max(-3.0, min(3.0, structural_score))

    # Recompute aggregation-prone regions
    agg_regions = _find_aggregation_prone_regions(
        corrected_smooth, _AGGREGATION_THRESHOLD
    )

    solubility_class = classify_solubility(structural_score)

    result = SolubilityResult(
        protein=protein,
        intrinsic_score=intrinsic.intrinsic_score,
        structural_score=round(structural_score, 4),
        overall_score=round(structural_score, 4),
        per_residue_scores=[round(s, 4) for s in corrected_smooth],
        aggregation_prone_regions=[
            (start, end, round(avg, 4)) for start, end, avg in agg_regions
        ],
        solubility_class=solubility_class,
        recommendations=[],
        method="camsol_structural",
    )
    result.recommendations = generate_solubility_recommendations(result)

    logger.debug(
        "CamSol structural: protein=%s, len=%d, intrinsic=%.4f, "
        "structural=%.4f, class=%s",
        protein[:10] + "..." if len(protein) > 10 else protein,
        len(protein),
        intrinsic.intrinsic_score,
        structural_score,
        solubility_class,
    )

    return result


def compute_solubility(
    protein: str,
    pdb_string: str | None = None,
    **kwargs: object,
) -> SolubilityResult:
    """Unified solubility computation interface.

    Calls compute_structural_solubility if PDB data is provided,
    otherwise compute_intrinsic_solubility.

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        pdb_string: Optional PDB file content as a string.
        **kwargs: Ignored (backward compatibility for callers passing
            extra keyword arguments like organism, structure_correction).

    Returns:
        SolubilityResult with solubility prediction.
    """
    if pdb_string is not None and pdb_string.strip():
        return compute_structural_solubility(protein, pdb_string)
    return compute_intrinsic_solubility(protein)


def find_solubility_mutations(
    protein: str,
    min_score: float = 0.0,
    **kwargs: object,
) -> list[dict]:
    """Find amino acid substitutions to improve solubility.

    For each position in aggregation-prone regions, tries all 19 possible
    substitutions and computes the change in intrinsic solubility score
    (delta-delta-solubility). Only conservative mutations (BLOSUM62 >= 0)
    are returned.

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        min_score: Minimum intrinsic score to target (default 0.0).
        **kwargs: Ignored (backward compatibility for callers passing
            extra keyword arguments like pdb_string).

    Returns:
        List of dicts with keys:
            position (int): 0-based residue position.
            wildtype (str): Original amino acid.
            mutant (str): Substituted amino acid.
            delta_solubility (float): Change in intrinsic score.
            blosum62 (int): BLOSUM62 score for the substitution.
        Sorted by delta_solubility descending (most improvement first).

    Raises:
        ValueError: If protein is empty.
    """
    protein = protein.upper().strip()

    if not protein:
        raise ValueError("Protein sequence must not be empty.")

    # Compute intrinsic solubility to find aggregation-prone regions
    result = compute_intrinsic_solubility(protein)

    # If already above threshold, no mutations needed
    if result.intrinsic_score >= min_score and not result.aggregation_prone_regions:
        return []

    mutations: list[dict] = []

    # Determine positions to target: all residues in aggregation-prone regions
    target_positions: set[int] = set()
    for start, end, _avg in result.aggregation_prone_regions:
        for pos in range(start, end):
            target_positions.add(pos)

    # Also consider positions with the lowest per-residue scores
    if not target_positions and result.intrinsic_score < min_score:
        # Target the worst-scoring residues
        n_target = max(1, len(protein) // 10)
        indexed = sorted(
            enumerate(result.per_residue_scores), key=lambda x: x[1]
        )
        for idx, _score in indexed[:n_target]:
            target_positions.add(idx)

    for pos in sorted(target_positions):
        if pos >= len(protein):
            continue
        wildtype = protein[pos]

        for mutant in _STANDARD_AAS:
            if mutant == wildtype:
                continue

            # Check BLOSUM62 conservation
            blosum = BLOSUM62.get(wildtype, {}).get(mutant, -10)
            if blosum < 0:
                continue

            # Compute delta solubility
            delta = _compute_mutation_delta(protein, pos, mutant)
            if delta > 0:
                mutations.append({
                    "position": pos,
                    "wildtype": wildtype,
                    "mutant": mutant,
                    "delta_solubility": round(delta, 4),
                    "blosum62": blosum,
                })

    # Sort by delta_solubility descending (most improvement first)
    mutations.sort(key=lambda m: m["delta_solubility"], reverse=True)

    return mutations


def classify_solubility(score: float) -> str:
    """Classify a solubility score into a solubility category.

    Args:
        score: CamSol solubility score (-3 to +3).

    Returns:
        One of: "highly_soluble", "soluble", "marginally_soluble", "insoluble".
    """
    if score > 1.5:
        return "highly_soluble"
    elif score > 0.0:
        return "soluble"
    elif score > -1.0:
        return "marginally_soluble"
    else:
        return "insoluble"


def generate_solubility_recommendations(result: SolubilityResult) -> list[str]:
    """Generate actionable recommendations for improving protein solubility.

    Based on aggregation-prone regions and overall score, suggests:
      - Surface mutations to charged residues
      - Proline substitutions at edges of aggregation regions
      - N-terminal tag additions
      - Avoiding specific patterns (e.g., long hydrophobic stretches)

    Args:
        result: SolubilityResult from a solubility computation.

    Returns:
        List of recommendation strings.
    """
    recommendations: list[str] = []
    protein = result.protein
    regions = result.aggregation_prone_regions

    # Overall assessment
    if result.overall_score > 1.5:
        recommendations.append(
            "Protein is highly soluble. No modifications necessary."
        )
        return recommendations

    if not regions:
        if result.overall_score < 0.0:
            recommendations.append(
                "Overall score is below zero but no discrete aggregation-prone "
                "regions detected. Consider increasing overall hydrophilicity "
                "with surface-exposed charged residues (K, E, D, R)."
            )
        return recommendations

    # Report aggregation-prone regions
    for start, end, avg in regions:
        region_seq = protein[start:end]
        hydrophobic_count = sum(
            1 for aa in region_seq
            if CAMSOL_HYDROPATHY.get(aa, 0.0) < 0
        )
        recommendations.append(
            f"Aggregation-prone region at positions {start+1}-{end} "
            f"(avg score: {avg:.2f}): '{region_seq}' — "
            f"{hydrophobic_count}/{len(region_seq)} hydrophobic residues."
        )

        # Suggest charged residue substitutions for hydrophobic residues
        for i in range(start, min(end, len(protein))):
            aa = protein[i]
            if CAMSOL_HYDROPATHY.get(aa, 0.0) < -0.3:
                best_charge = _best_charged_substitution(aa)
                if best_charge:
                    recommendations.append(
                        f"  Position {i+1}: Consider {aa}->{best_charge} "
                        f"substitution (charged residue improves solubility)."
                    )

        # Suggest proline at edges
        if end - start >= _MIN_AGGREGATION_REGION_LENGTH:
            edge_left = max(0, start)
            edge_right = min(len(protein) - 1, end - 1)
            if protein[edge_left] != "P":
                recommendations.append(
                    f"  Consider P substitution at position {edge_left+1} "
                    f"(N-terminal edge of aggregation region) to break "
                    f"beta-sheet propensity."
                )
            if protein[edge_right] != "P":
                recommendations.append(
                    f"  Consider P substitution at position {edge_right+1} "
                    f"(C-terminal edge of aggregation region) to break "
                    f"beta-sheet propensity."
                )

    # Check for long hydrophobic stretches
    _check_hydrophobic_stretches(protein, recommendations)

    # Check net charge
    net_charge = sum(CAMSOL_CHARGE.get(aa, 0.0) for aa in protein)
    if abs(net_charge) < 2 and len(protein) > 30:
        recommendations.append(
            f"Low net charge ({net_charge:+.1f}). Increasing charged residue "
            f"content (K, R, D, E) can improve solubility through charge "
            f"repulsion."
        )

    # N-terminal tag suggestion for marginally soluble or insoluble
    if result.overall_score < 0.0:
        recommendations.append(
            "Consider adding an N-terminal solubility tag (e.g., MBP, GST, "
            "SUMO, or His6) to improve expression solubility."
        )

    # Avoid specific patterns: long hydrophobic tetrapeptide
    for i in range(len(protein) - 3):
        tetrapeptide = protein[i:i + 4]
        if all(CAMSOL_HYDROPATHY.get(aa, 0.0) < -0.2 for aa in tetrapeptide):
            recommendations.append(
                f"Long hydrophobic stretch at positions {i+1}-{i+4} "
                f"('{tetrapeptide}'). Consider breaking this pattern with "
                f"a charged or polar residue."
            )
            break  # only report the first one

    return recommendations


# ────────────────────────────────────────────────────────────
# Internal helper functions
# ────────────────────────────────────────────────────────────

def _sliding_window_smooth(scores: list[float], window: int) -> list[float]:
    """Apply sliding window averaging to a list of scores.

    At the boundaries, uses a truncated window (no zero-padding).

    Args:
        scores: Per-residue scores.
        window: Window size (odd values give symmetric smoothing).

    Returns:
        Smoothed scores of the same length.
    """
    if not scores:
        return []
    if window <= 1:
        return list(scores)
    if len(scores) <= window:
        return [sum(scores) / len(scores)] * len(scores)

    n = len(scores)
    half = window // 2
    result = []

    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        avg = sum(scores[start:end]) / (end - start)
        result.append(avg)

    return result


def _find_aggregation_prone_regions(
    scores: list[float],
    threshold: float,
) -> list[tuple[int, int, float]]:
    """Identify consecutive residues with scores below threshold.

    Args:
        scores: Per-residue solubility scores.
        threshold: Score threshold for aggregation-prone residues.

    Returns:
        List of (start, end, avg_score) tuples. start is inclusive,
        end is exclusive.
    """
    if not scores:
        return []

    regions: list[tuple[int, int, float]] = []
    in_region = False
    region_start = 0

    for i, score in enumerate(scores):
        if score < threshold:
            if not in_region:
                region_start = i
                in_region = True
        else:
            if in_region:
                if i - region_start >= _MIN_AGGREGATION_REGION_LENGTH:
                    avg = sum(scores[region_start:i]) / (i - region_start)
                    regions.append((region_start, i, avg))
                in_region = False

    # Handle region extending to the end
    if in_region:
        end = len(scores)
        if end - region_start >= _MIN_AGGREGATION_REGION_LENGTH:
            avg = sum(scores[region_start:end]) / (end - region_start)
            regions.append((region_start, end, avg))

    # Merge overlapping/adjacent regions
    merged: list[tuple[int, int, float]] = []
    for start, end, avg in regions:
        if merged and start <= merged[-1][1]:
            # Merge with previous
            prev_start, prev_end, _prev_avg = merged.pop()
            new_start = prev_start
            new_end = max(prev_end, end)
            combined_scores = scores[new_start:new_end]
            new_avg = sum(combined_scores) / len(combined_scores)
            merged.append((new_start, new_end, new_avg))
        else:
            merged.append((start, end, avg))

    return merged


def _compute_mutation_delta(
    protein: str,
    position: int,
    mutant: str,
) -> float:
    """Compute the change in intrinsic solubility from a single substitution.

    Computes the delta by recalculating the per-residue score at the
    mutation position and applying local smoothing effects.

    Args:
        protein: Original protein sequence.
        position: 0-based residue position to mutate.
        mutant: Mutant amino acid (1-letter code).

    Returns:
        Delta solubility score (positive = improvement).
    """
    wildtype = protein[position]

    # Compute score for wildtype residue
    wt_score = _compute_single_residue_score(wildtype)

    # Compute score for mutant residue
    mt_score = _compute_single_residue_score(mutant)

    # Direct delta at the position
    direct_delta = mt_score - wt_score

    # Also consider the effect on neighbors through smoothing
    # Check a 7-residue window centered on the position
    n = len(protein)
    window = 7
    half = window // 2

    # Build a local patch with the mutation
    mut_protein = protein[:position] + mutant + protein[position + 1:]

    # Compute local scores for wildtype and mutant
    wt_local = 0.0
    mt_local = 0.0
    count = 0

    for i in range(max(0, position - half), min(n, position + half + 1)):
        wt_local += _compute_single_residue_score(protein[i])
        mt_local += _compute_single_residue_score(mut_protein[i])
        count += 1

    if count > 0:
        wt_local /= count
        mt_local /= count

    # Weighted combination: 60% direct delta, 40% local smoothing delta
    local_delta = mt_local - wt_local
    delta = 0.6 * direct_delta + 0.4 * local_delta

    return delta


def _compute_single_residue_score(aa: str) -> float:
    """Compute the raw intrinsic score contribution for a single amino acid.

    Args:
        aa: Single-letter amino acid code.

    Returns:
        Weighted intrinsic score contribution.
    """
    hydro = CAMSOL_HYDROPATHY.get(aa, 0.0)
    charge = CAMSOL_CHARGE.get(aa, 0.0)
    alpha = CAMSOL_ALPHA_HELIX.get(aa, 0.0)
    beta = CAMSOL_BETA_STRAND.get(aa, 0.0)

    if aa == "P":
        progly = 1.0
    elif aa == "G":
        progly = -0.5
    else:
        progly = 0.0

    return (
        _WEIGHT_HYDROPATHY * hydro
        + _WEIGHT_CHARGE * charge
        + _WEIGHT_ALPHA_HELIX * alpha
        + _WEIGHT_BETA_STRAND * (-beta)
        + _WEIGHT_PROGLY * progly
    )


def _best_charged_substitution(aa: str) -> str | None:
    """Find the best charged residue substitution for a hydrophobic residue.

    Prefers the substitution with the highest BLOSUM62 score among
    charged residues (K, R, D, E, H).

    Args:
        aa: Hydrophobic amino acid (1-letter code).

    Returns:
        Best charged substitution, or None if no acceptable substitution.
    """
    charged = ["K", "R", "D", "E", "H"]
    best = None
    best_blosum = -10

    for candidate in charged:
        blosum = BLOSUM62.get(aa, {}).get(candidate, -10)
        if blosum >= -1 and blosum > best_blosum:
            best_blosum = blosum
            best = candidate

    return best


def _check_hydrophobic_stretches(
    protein: str,
    recommendations: list[str],
) -> None:
    """Check for long consecutive hydrophobic stretches in the protein.

    Appends warnings to the recommendations list.

    Args:
        protein: Protein sequence.
        recommendations: List to append recommendations to (mutated in place).
    """
    hydrophobic = set()
    for aa, val in CAMSOL_HYDROPATHY.items():
        if val < -0.2:
            hydrophobic.add(aa)

    stretch_start = None
    for i, aa in enumerate(protein):
        if aa in hydrophobic:
            if stretch_start is None:
                stretch_start = i
        else:
            if stretch_start is not None:
                length = i - stretch_start
                if length >= 7:
                    recommendations.append(
                        f"Very long hydrophobic stretch at positions "
                        f"{stretch_start+1}-{i} "
                        f"('{protein[stretch_start:i]}', {length} residues). "
                        f"Long hydrophobic stretches strongly promote "
                        f"aggregation. Consider inserting charged/polar "
                        f"residues every 5-6 positions."
                    )
                stretch_start = None

    # Handle stretch at end
    if stretch_start is not None:
        length = len(protein) - stretch_start
        if length >= 7:
            recommendations.append(
                f"Very long hydrophobic stretch at positions "
                f"{stretch_start+1}-{len(protein)} "
                f"('{protein[stretch_start:]}', {length} residues). "
                f"Long hydrophobic stretches strongly promote aggregation. "
                f"Consider inserting charged/polar residues every 5-6 positions."
            )


# ────────────────────────────────────────────────────────────
# PDB parsing helpers
# ────────────────────────────────────────────────────────────

def _parse_pdb_ca_coords(
    pdb_string: str,
) -> list[tuple[float, float, float]]:
    """Extract CA (alpha carbon) coordinates from PDB content.

    Args:
        pdb_string: PDB file content.

    Returns:
        List of (x, y, z) tuples for each CA atom, in order of appearance.
    """
    coords: list[tuple[float, float, float]] = []

    for line in pdb_string.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) >= 54:
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except (ValueError, IndexError):
                        continue

    return coords


def _approximate_sasa(
    ca_coords: list[tuple[float, float, float]],
    cutoff: float = 10.0,
) -> list[float]:
    """Approximate relative SASA using CA neighbor counting.

    For each CA atom, counts the number of neighbor CA atoms within
    the cutoff distance. More neighbors = more buried = lower SASA.

    The raw neighbor count is converted to a relative SASA value in [0, 1]
    using a sigmoid-like transformation.

    Args:
        ca_coords: List of (x, y, z) CA coordinates.
        cutoff: Distance cutoff in Angstroms for neighbor detection.

    Returns:
        List of relative SASA values in [0, 1] for each residue.
    """
    n = len(ca_coords)
    neighbor_counts = [0] * n

    for i in range(n):
        xi, yi, zi = ca_coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = ca_coords[j]
            dist_sq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if dist_sq < cutoff * cutoff:
                neighbor_counts[i] += 1
                neighbor_counts[j] += 1

    # Convert neighbor counts to relative SASA using sigmoid transformation
    # More neighbors -> more buried -> lower SASA
    # Typical CA neighbor count: 0-25 for buried, 0-8 for exposed
    max_neighbors = max(neighbor_counts) if neighbor_counts else 1
    if max_neighbors == 0:
        max_neighbors = 1

    sasa: list[float] = []
    for count in neighbor_counts:
        # Sigmoid: 1 / (1 + exp(k * (count - midpoint)))
        # At count=0, SASA ~ 1.0 (fully exposed)
        # At count=20, SASA ~ 0.0 (fully buried)
        k = 0.4
        midpoint = 12.0
        relative_sasa = 1.0 / (1.0 + math.exp(k * (count - midpoint)))
        sasa.append(relative_sasa)

    return sasa


def _parse_disulfide_bonds(pdb_string: str) -> set[int]:
    """Parse SSBOND records from PDB to find residues in disulfide bonds.

    Args:
        pdb_string: PDB file content.

    Returns:
        Set of 0-based residue indices involved in disulfide bonds.
    """
    residues: set[int] = set()

    for line in pdb_string.splitlines():
        if line.startswith("SSBOND"):
            try:
                # SSBOND record format:
                # COLUMNS 18-21: Sequence number of first residue
                # COLUMNS 32-35: Sequence number of second residue
                res1 = int(line[17:21].strip())
                res2 = int(line[31:35].strip())
                # Convert 1-based PDB numbering to 0-based
                residues.add(res1 - 1)
                residues.add(res2 - 1)
            except (ValueError, IndexError):
                continue

    return residues
