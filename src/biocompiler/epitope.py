"""
BioCompiler B-cell Epitope Prediction
======================================
Linear and conformational B-cell epitope prediction using multiple
classical scales and methods.

Methods implemented:
  - Kolaskar-Tongaonkar antigenicity
  - Parker hydrophilicity
  - Chou-Fasman β-turn propensity
  - Emini surface accessibility (EEA)
  - BepiPred-like composite method
  - Conformational epitope prediction from PDB structure

References:
  - Kolaskar & Tongaonkar, FEBS Lett 1990; 276:172-174
  - Parker et al., Biochemistry 1986; 25:5424-5432
  - Chou & Fasman, Biochemistry 1974; 13:222-245
  - Emini et al., J Virol 1985; 55:836-839
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Amino acid scales
# ────────────────────────────────────────────────────────────

ANTIGENICITY_SCALE: dict[str, float] = {
    # Kolaskar-Tongaonkar antigenicity scale
    "A": 1.064, "R": 1.008, "N": 0.873, "D": 1.026,
    "C": 1.412, "E": 0.895, "Q": 1.091, "G": 0.842,
    "H": 1.105, "I": 1.142, "L": 1.170, "K": 0.933,
    "M": 1.207, "F": 1.279, "P": 0.658, "S": 0.772,
    "T": 0.789, "W": 1.190, "Y": 1.161, "V": 1.132,
}

PARKER_SCALE: dict[str, float] = {
    # Parker hydrophilicity scale (Parker et al. 1986)
    # Positive = hydrophilic (surface/epitope), negative = hydrophobic (buried)
    "A": -1.01, "R": 1.40, "N": 0.81, "D": 1.21,
    "C": -1.20, "E": 1.64, "Q": 0.96, "G": -0.16,
    "H": 0.56, "I": -1.42, "L": -1.42, "K": 1.73,
    "M": -1.27, "F": -1.42, "P": 0.26, "S": 0.52,
    "T": -0.19, "W": -1.07, "Y": -0.31, "V": -1.07,
}

CHOU_FASMAN_TURN: dict[str, float] = {
    # Chou-Fasman β-turn propensity (Chou & Fasman 1974)
    "A": 0.060, "R": 0.095, "N": 0.147, "D": 0.161,
    "C": 0.108, "E": 0.056, "Q": 0.098, "G": 0.102,
    "H": 0.140, "I": 0.043, "L": 0.053, "K": 0.101,
    "M": 0.068, "F": 0.059, "P": 0.301, "S": 0.120,
    "T": 0.086, "W": 0.077, "Y": 0.114, "V": 0.050,
}

EMINI_SCALE: dict[str, float] = {
    # Emini surface probability scale (Emini et al. 1985)
    "A": 0.510, "R": 1.008, "N": 0.849, "D": 0.628,
    "C": 0.358, "E": 0.977, "Q": 0.993, "G": 0.471,
    "H": 0.873, "I": 0.296, "L": 0.332, "K": 1.027,
    "M": 0.411, "F": 0.328, "P": 0.709, "S": 0.643,
    "T": 0.549, "W": 0.307, "Y": 0.361, "V": 0.265,
}

_FLEXIBILITY_SCALE: dict[str, float] = {
    # Flexibility scale used by the BepiPred-like composite method
    # Glycine and proline are most flexible; bulky aromatics least
    "A": 0.360, "R": 0.530, "N": 0.460, "D": 0.510,
    "C": 0.350, "E": 0.500, "Q": 0.490, "G": 0.540,
    "H": 0.320, "I": 0.460, "L": 0.370, "K": 0.470,
    "M": 0.300, "F": 0.310, "P": 0.510, "S": 0.510,
    "T": 0.440, "W": 0.310, "Y": 0.420, "V": 0.390,
}

ALL_SCALES: dict[str, dict[str, float]] = {
    "kolaskar_tongaonkar": ANTIGENICITY_SCALE,
    "parker_hydrophilicity": PARKER_SCALE,
    "chou_fasman": CHOU_FASMAN_TURN,
    "eea": EMINI_SCALE,
    "bepipred_flexibility": _FLEXIBILITY_SCALE,
}

# 3-letter to 1-letter amino acid mapping (for PDB parsing)
_THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class EpitopeRegion:
    """A predicted B-cell epitope region."""

    start: int          # 0-indexed start position
    end: int            # exclusive end position
    peptide: str
    score: float        # 0 to 1
    method: str         # "kolaskar_tongaonkar", "bepipred", "parker_hydrophilicity",
                        # "chou_fasman", "eea", "consensus", "conformational"
    is_linear: bool
    properties: dict = field(default_factory=dict)


@dataclass
class EpitopePredictionResult:
    """Combined B-cell epitope prediction result."""

    protein: str
    linear_epitopes: list[EpitopeRegion]
    conformational_epitopes: list[EpitopeRegion]   # empty if no structure provided
    per_residue_score: list[float]                  # combined epitope propensity per residue
    epitope_coverage: float                         # fraction of residues in predicted epitopes
    methods_used: list[str]


# ────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────

def _sliding_window_average(values: list[float], window: int) -> list[float]:
    """Compute sliding window average of a list of values."""
    n = len(values)
    if n == 0 or window <= 0:
        return []
    half = window // 2
    result: list[float] = []
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        result.append(sum(values[start:end]) / (end - start))
    return result


def _normalize_01(scores: list[float]) -> list[float]:
    """Min-max normalize a list of scores to [0, 1]."""
    if not scores:
        return []
    min_s = min(scores)
    max_s = max(scores)
    if max_s - min_s > 1e-10:
        return [(s - min_s) / (max_s - min_s) for s in scores]
    return [0.5] * len(scores)


def _find_regions(
    norm_scores: list[float],
    protein: str,
    method: str,
    threshold: float = 0.5,
    min_length: int = 6,
    is_linear: bool = True,
    extra_props: dict | None = None,
    raw_scores: list[float] | None = None,
) -> list[EpitopeRegion]:
    """Find contiguous regions where normalized scores exceed threshold.

    Args:
        norm_scores: Per-residue scores normalized to [0, 1].
        protein: Amino acid sequence.
        method: Method name for the EpitopeRegion.
        threshold: Normalized score threshold.
        min_length: Minimum region length in residues.
        is_linear: Whether the epitope is linear.
        extra_props: Additional properties to include.
        raw_scores: Raw (un-normalized) scores for property recording.

    Returns:
        List of EpitopeRegion objects.
    """
    if not norm_scores:
        return []

    regions: list[EpitopeRegion] = []
    in_region = False
    region_start = 0

    for i, s in enumerate(norm_scores):
        if s >= threshold and not in_region:
            in_region = True
            region_start = i
        elif s < threshold and in_region:
            in_region = False
            if i - region_start >= min_length:
                _add_region(
                    regions, norm_scores, raw_scores, protein,
                    region_start, i, method, is_linear, extra_props,
                )

    if in_region:
        i = len(norm_scores)
        if i - region_start >= min_length:
            _add_region(
                regions, norm_scores, raw_scores, protein,
                region_start, i, method, is_linear, extra_props,
            )

    return regions


def _add_region(
    regions: list[EpitopeRegion],
    norm_scores: list[float],
    raw_scores: list[float] | None,
    protein: str,
    start: int,
    end: int,
    method: str,
    is_linear: bool,
    extra_props: dict | None,
) -> None:
    """Append a single EpitopeRegion to the list."""
    length = end - start
    avg_score = sum(norm_scores[start:end]) / length
    props: dict = {}
    if raw_scores is not None:
        props["raw_score_avg"] = round(sum(raw_scores[start:end]) / length, 4)
    if extra_props:
        props.update(extra_props)
    regions.append(EpitopeRegion(
        start=start,
        end=end,
        peptide=protein[start:end],
        score=round(avg_score, 4),
        method=method,
        is_linear=is_linear,
        properties=props,
    ))


# ────────────────────────────────────────────────────────────
# Prediction methods
# ────────────────────────────────────────────────────────────

def predict_kolaskar_tongaonkar(
    protein: str,
    window: int = 7,
    threshold: float = 1.0,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using the Kolaskar-Tongaonkar antigenicity scale.

    Kolaskar & Tongaonkar (1990) developed an antigenicity scale based on
    the frequency of amino acids in known antigenic determinants. Regions
    with higher antigenicity scores are more likely to be epitopes.

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for averaging.
        threshold: Antigenicity threshold (default 1.0, the approximate
                   mean of the scale). Regions above this are predicted
                   as epitopes.

    Returns:
        List of EpitopeRegion predictions.
    """
    if not protein:
        return []

    protein = protein.upper()

    # Compute raw antigenicity scores per residue
    raw_scores = [ANTIGENICITY_SCALE.get(aa, 1.0) for aa in protein]

    # Sliding window average
    smoothed = _sliding_window_average(raw_scores, window)

    # Normalize to [0, 1]
    norm_scores = _normalize_01(smoothed)

    # Convert raw threshold to normalized threshold
    min_s = min(smoothed) if smoothed else 0.0
    max_s = max(smoothed) if smoothed else 1.0
    if max_s - min_s > 1e-10:
        norm_threshold = (threshold - min_s) / (max_s - min_s)
        norm_threshold = max(0.0, min(1.0, norm_threshold))
    else:
        norm_threshold = 0.5

    return _find_regions(
        norm_scores, protein, "kolaskar_tongaonkar",
        threshold=norm_threshold,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": threshold},
        raw_scores=smoothed,
    )


def predict_parker_hydrophilicity(
    protein: str,
    window: int = 7,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using the Parker hydrophilicity scale.

    Parker et al. (1986) showed that hydrophilic regions of proteins
    tend to be on the surface and are more likely to be B-cell epitopes.
    Positive values indicate hydrophilic (surface/epitope) regions;
    negative values indicate hydrophobic (buried) regions.

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for averaging.

    Returns:
        List of EpitopeRegion predictions for hydrophilic regions.
    """
    if not protein:
        return []

    protein = protein.upper()

    raw_scores = [PARKER_SCALE.get(aa, 0.0) for aa in protein]
    smoothed = _sliding_window_average(raw_scores, window)
    norm_scores = _normalize_01(smoothed)

    # Threshold: above the midpoint of the normalized scale
    return _find_regions(
        norm_scores, protein, "parker_hydrophilicity",
        threshold=0.5,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": "above_mean"},
        raw_scores=smoothed,
    )


def predict_chou_fasman_beta_turn(
    protein: str,
    window: int = 7,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using Chou-Fasman β-turn propensity.

    Chou & Fasman (1974) showed that β-turns are often surface-exposed
    and correspond to B-cell epitope regions. Higher turn propensity
    indicates a greater likelihood of being an epitope.

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for averaging.

    Returns:
        List of EpitopeRegion predictions for high turn-propensity regions.
    """
    if not protein:
        return []

    protein = protein.upper()

    raw_scores = [CHOU_FASMAN_TURN.get(aa, 0.0) for aa in protein]
    smoothed = _sliding_window_average(raw_scores, window)
    norm_scores = _normalize_01(smoothed)

    return _find_regions(
        norm_scores, protein, "chou_fasman",
        threshold=0.5,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": "above_mean"},
        raw_scores=smoothed,
    )


def predict_eea(
    protein: str,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using Emini surface accessibility (EEA).

    Emini et al. (1985) developed a surface probability method based on
    the statistical analysis of surface accessibility in protein structures.
    Surface-accessible residues are more likely to be B-cell epitopes.

    The Emini formula computes surface probability as a geometric mean
    over a sliding window:
        S(n) = exp( (1/w) * Σ ln(f(i)) )
    where f(i) is the surface probability for the amino acid at position i.

    Args:
        protein: Amino acid sequence (1-letter codes).

    Returns:
        List of EpitopeRegion predictions for surface-accessible regions.
    """
    if not protein:
        return []

    protein = protein.upper()
    n = len(protein)
    window = 6   # standard EEA window size
    half = window // 2

    # Compute Emini surface probability scores using geometric mean
    raw_scores: list[float] = []
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        log_sum = 0.0
        count = 0
        for j in range(start, end):
            prob = EMINI_SCALE.get(protein[j], 0.5)
            log_sum += math.log(max(prob, 1e-10))
            count += 1
        raw_scores.append(math.exp(log_sum / count) if count > 0 else 0.5)

    # Normalize to [0, 1]
    norm_scores = _normalize_01(raw_scores)

    # Use the original EEA threshold of 1.0 where possible.
    # If the raw threshold of 1.0 is above the max score (common with
    # geometric mean), fall back to the mean-based threshold (0.5 on
    # the normalized scale).
    min_s = min(raw_scores) if raw_scores else 0.0
    max_s = max(raw_scores) if raw_scores else 1.0
    if max_s >= 1.0 and max_s - min_s > 1e-10:
        norm_threshold = (1.0 - min_s) / (max_s - min_s)
        norm_threshold = max(0.0, min(1.0, norm_threshold))
        raw_threshold = 1.0
    else:
        norm_threshold = 0.5
        raw_threshold = round(
            (max_s + min_s) / 2.0, 4
        ) if raw_scores else 1.0

    return _find_regions(
        norm_scores, protein, "eea",
        threshold=norm_threshold,
        min_length=6,
        is_linear=True,
        extra_props={"window": window, "threshold": raw_threshold},
        raw_scores=raw_scores,
    )


def predict_bepipred_like(
    protein: str,
    window: int = 9,
) -> list[EpitopeRegion]:
    """Predict B-cell epitopes using a simplified BepiPred-like composite method.

    Combines three properties into a composite score:
        Score = 0.4 * hydrophilicity + 0.3 * flexibility + 0.3 * surface_accessibility

    Hydrophilicity is from the Parker scale, surface accessibility from the
    Emini scale, and flexibility is estimated from local sequence composition
    (glycine and proline content increases flexibility).

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for smoothing.

    Returns:
        List of EpitopeRegion predictions.
    """
    if not protein:
        return []

    protein = protein.upper()

    # Compute individual property scores
    hydro_raw = [PARKER_SCALE.get(aa, 0.0) for aa in protein]
    flex_raw = [_FLEXIBILITY_SCALE.get(aa, 0.4) for aa in protein]
    surf_raw = [EMINI_SCALE.get(aa, 0.5) for aa in protein]

    # Smooth each property
    hydro_smooth = _sliding_window_average(hydro_raw, window)
    flex_smooth = _sliding_window_average(flex_raw, window)
    surf_smooth = _sliding_window_average(surf_raw, window)

    # Normalize each to [0, 1]
    hydro_norm = _normalize_01(hydro_smooth)
    flex_norm = _normalize_01(flex_smooth)
    surf_norm = _normalize_01(surf_smooth)

    # Composite score
    composite: list[float] = [
        0.4 * h + 0.3 * f + 0.3 * s
        for h, f, s in zip(hydro_norm, flex_norm, surf_norm)
    ]

    # Find regions above 0.5 threshold
    regions: list[EpitopeRegion] = []
    in_region = False
    region_start = 0
    min_length = 6

    for i, s in enumerate(composite):
        if s >= 0.5 and not in_region:
            in_region = True
            region_start = i
        elif s < 0.5 and in_region:
            in_region = False
            if i - region_start >= min_length:
                length = i - region_start
                avg = sum(composite[region_start:i]) / length
                h_avg = sum(hydro_norm[region_start:i]) / length
                f_avg = sum(flex_norm[region_start:i]) / length
                s_avg = sum(surf_norm[region_start:i]) / length
                regions.append(EpitopeRegion(
                    start=region_start,
                    end=i,
                    peptide=protein[region_start:i],
                    score=round(avg, 4),
                    method="bepipred",
                    is_linear=True,
                    properties={
                        "hydrophilicity_avg": round(h_avg, 4),
                        "flexibility_avg": round(f_avg, 4),
                        "surface_avg": round(s_avg, 4),
                        "window": window,
                        "weights": {
                            "hydrophilicity": 0.4,
                            "flexibility": 0.3,
                            "surface": 0.3,
                        },
                    },
                ))

    if in_region:
        i = len(composite)
        if i - region_start >= min_length:
            length = i - region_start
            avg = sum(composite[region_start:i]) / length
            h_avg = sum(hydro_norm[region_start:i]) / length
            f_avg = sum(flex_norm[region_start:i]) / length
            s_avg = sum(surf_norm[region_start:i]) / length
            regions.append(EpitopeRegion(
                start=region_start,
                end=i,
                peptide=protein[region_start:i],
                score=round(avg, 4),
                method="bepipred",
                is_linear=True,
                properties={
                    "hydrophilicity_avg": round(h_avg, 4),
                    "flexibility_avg": round(f_avg, 4),
                    "surface_avg": round(s_avg, 4),
                    "window": window,
                    "weights": {
                        "hydrophilicity": 0.4,
                        "flexibility": 0.3,
                        "surface": 0.3,
                    },
                },
            ))

    return regions


def predict_conformational_epitopes(
    pdb_string: str,
    distance_cutoff: float = 6.0,
) -> list[EpitopeRegion]:
    """Predict conformational B-cell epitopes from a PDB structure.

    Identifies surface patches on the protein structure and scores them
    by hydrophilicity, charge, and flexibility. Surface residues are
    identified by having fewer than 15 Cα neighbors within 12 Å.
    Adjacent surface residues (within distance_cutoff in sequence) are
    clustered into patches.

    Args:
        pdb_string: PDB file content as a string.
        distance_cutoff: Maximum sequence gap to cluster adjacent surface
                        residues into patches (default 6.0).

    Returns:
        List of EpitopeRegion predictions with is_linear=False.
    """
    if not pdb_string:
        return []

    # ── Parse PDB: extract Cα coordinates ──
    ca_atoms: list[tuple[int, float, float, float, str]] = []

    for line in pdb_string.splitlines():
        line = line.rstrip()
        if not line.startswith("ATOM"):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        try:
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except (ValueError, IndexError):
            continue

        resname = line[17:20].strip() if len(line) >= 20 else ""
        aa = _THREE_TO_ONE.get(resname, "X")
        ca_atoms.append((resnum, x, y, z, aa))

    if len(ca_atoms) < 3:
        logger.warning(
            "Too few Cα atoms (%d) in PDB for conformational epitope prediction",
            len(ca_atoms),
        )
        return []

    # ── Identify surface residues: Cα with < 15 neighbors within 12 Å ──
    neighbor_cutoff = 12.0
    max_neighbors = 15

    surface_indices: set[int] = set()
    for i, (_, xi, yi, zi, _) in enumerate(ca_atoms):
        neighbors = 0
        for j, (_, xj, yj, zj, _) in enumerate(ca_atoms):
            if i == j:
                continue
            dist_sq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if dist_sq <= neighbor_cutoff ** 2:
                neighbors += 1
                if neighbors >= max_neighbors:
                    break
        if neighbors < max_neighbors:
            surface_indices.add(i)

    if not surface_indices:
        return []

    # ── Cluster adjacent surface residues into patches ──
    sorted_surface = sorted(surface_indices)
    patches: list[list[int]] = []
    current_patch = [sorted_surface[0]]

    for k in range(1, len(sorted_surface)):
        if sorted_surface[k] - sorted_surface[k - 1] <= distance_cutoff:
            current_patch.append(sorted_surface[k])
        else:
            if len(current_patch) >= 3:
                patches.append(current_patch)
            current_patch = [sorted_surface[k]]

    if len(current_patch) >= 3:
        patches.append(current_patch)

    # ── Score patches by hydrophilicity + charge + flexibility ──
    charged_aas = {"R", "K", "D", "E"}
    flexible_aas = {"G", "P"}

    epitopes: list[EpitopeRegion] = []

    for patch in patches:
        aas: list[str] = []
        hydro_sum = 0.0
        charge_count = 0
        flex_count = 0

        for idx in patch:
            _, _, _, _, aa = ca_atoms[idx]
            aas.append(aa)
            hydro_sum += PARKER_SCALE.get(aa, 0.0)
            if aa in charged_aas:
                charge_count += 1
            if aa in flexible_aas:
                flex_count += 1

        n_res = len(patch)
        hydro_avg = hydro_sum / n_res
        charge_frac = charge_count / n_res
        flex_frac = flex_count / n_res

        # Normalize hydrophilicity: Parker scale ≈ [-1.5, 1.8]
        hydro_norm = max(0.0, min(1.0, (hydro_avg + 1.5) / 3.3))

        # Combined score
        score = 0.4 * hydro_norm + 0.3 * charge_frac + 0.3 * flex_frac

        # Map patch to sequence positions (0-based)
        first_resnum = ca_atoms[patch[0]][0]
        last_resnum = ca_atoms[patch[-1]][0]
        start = first_resnum - 1   # convert to 0-based
        end = last_resnum          # exclusive

        epitopes.append(EpitopeRegion(
            start=start,
            end=end,
            peptide="".join(aas),
            score=round(score, 4),
            method="conformational",
            is_linear=False,
            properties={
                "hydrophilicity_avg": round(hydro_avg, 4),
                "charge_fraction": round(charge_frac, 4),
                "flexibility_fraction": round(flex_frac, 4),
                "surface_residue_count": n_res,
                "pdb_residue_range": (first_resnum, last_resnum),
            },
        ))

    return epitopes


# ────────────────────────────────────────────────────────────
# Method dispatch
# ────────────────────────────────────────────────────────────

_METHOD_MAP: dict[str, object] = {
    "kolaskar_tongaonkar": predict_kolaskar_tongaonkar,
    "parker_hydrophilicity": predict_parker_hydrophilicity,
    "chou_fasman_beta_turn": predict_chou_fasman_beta_turn,
    "eea": predict_eea,
    "bepipred": predict_bepipred_like,
}


# ────────────────────────────────────────────────────────────
# Combined prediction
# ────────────────────────────────────────────────────────────

def predict_epitopes(
    protein: str,
    pdb_string: str | None = None,
    methods: list[str] | None = None,
) -> EpitopePredictionResult:
    """Run multiple B-cell epitope prediction methods and combine results.

    Default methods: kolaskar_tongaonkar, parker_hydrophilicity,
    chou_fasman_beta_turn. If a PDB structure is provided,
    conformational epitope prediction is also run.

    Per-residue scores are computed as the average across all methods
    that predict each residue as an epitope. Consensus epitope regions
    are those predicted by ≥ 2 methods.

    Args:
        protein: Amino acid sequence (1-letter codes).
        pdb_string: Optional PDB file content for conformational prediction.
        methods: List of method names to use. Defaults to
                 ["kolaskar_tongaonkar", "parker_hydrophilicity",
                  "chou_fasman_beta_turn"].

    Returns:
        EpitopePredictionResult with combined predictions.
    """
    if not protein:
        return EpitopePredictionResult(
            protein="",
            linear_epitopes=[],
            conformational_epitopes=[],
            per_residue_score=[],
            epitope_coverage=0.0,
            methods_used=[],
        )

    protein = protein.upper()
    n = len(protein)

    if methods is None:
        methods = [
            "kolaskar_tongaonkar",
            "parker_hydrophilicity",
            "chou_fasman_beta_turn",
        ]

    # Validate method names
    valid_methods = [m for m in methods if m in _METHOD_MAP]
    if not valid_methods:
        logger.warning("No valid methods specified for epitope prediction")
        return EpitopePredictionResult(
            protein=protein,
            linear_epitopes=[],
            conformational_epitopes=[],
            per_residue_score=[0.0] * n,
            epitope_coverage=0.0,
            methods_used=[],
        )

    # ── Run each linear method ──
    all_linear_epitopes: list[EpitopeRegion] = []
    residue_method_count: list[int] = [0] * n
    residue_score_sum: list[float] = [0.0] * n

    for method_name in valid_methods:
        func = _METHOD_MAP[method_name]
        try:
            epitopes = func(protein)  # type: ignore[operator]
            all_linear_epitopes.extend(epitopes)

            # Track per-residue predictions
            for ep in epitopes:
                for pos in range(ep.start, ep.end):
                    if 0 <= pos < n:
                        residue_method_count[pos] += 1
                        residue_score_sum[pos] += ep.score
        except Exception as e:
            logger.warning("Method %s failed: %s", method_name, e)

    # ── Per-residue combined score ──
    per_residue_score: list[float] = []
    for i in range(n):
        if residue_method_count[i] > 0:
            per_residue_score.append(
                residue_score_sum[i] / residue_method_count[i]
            )
        else:
            per_residue_score.append(0.0)

    # ── Consensus epitope regions (≥ 2 methods) ──
    consensus_epitopes: list[EpitopeRegion] = []
    in_region = False
    region_start = 0
    min_consensus = 2
    min_length = 6

    for i in range(n):
        if residue_method_count[i] >= min_consensus and not in_region:
            in_region = True
            region_start = i
        elif residue_method_count[i] < min_consensus and in_region:
            in_region = False
            if i - region_start >= min_length:
                length = i - region_start
                avg_score = sum(per_residue_score[region_start:i]) / length
                method_count_avg = sum(
                    residue_method_count[region_start:i]
                ) / length
                consensus_epitopes.append(EpitopeRegion(
                    start=region_start,
                    end=i,
                    peptide=protein[region_start:i],
                    score=round(avg_score, 4),
                    method="consensus",
                    is_linear=True,
                    properties={
                        "method_count_avg": round(method_count_avg, 4),
                        "contributing_methods": min_consensus,
                    },
                ))

    if in_region:
        i = n
        if i - region_start >= min_length:
            length = i - region_start
            avg_score = sum(per_residue_score[region_start:i]) / length
            method_count_avg = sum(
                residue_method_count[region_start:i]
            ) / length
            consensus_epitopes.append(EpitopeRegion(
                start=region_start,
                end=i,
                peptide=protein[region_start:i],
                score=round(avg_score, 4),
                method="consensus",
                is_linear=True,
                properties={
                    "method_count_avg": round(method_count_avg, 4),
                    "contributing_methods": min_consensus,
                },
            ))

    # Combine individual + consensus epitopes
    final_linear = list(all_linear_epitopes) + consensus_epitopes

    # ── Conformational prediction ──
    conformational_epitopes: list[EpitopeRegion] = []
    if pdb_string:
        try:
            conformational_epitopes = predict_conformational_epitopes(pdb_string)
        except Exception as e:
            logger.warning("Conformational epitope prediction failed: %s", e)

    # ── Epitope coverage ──
    epitope_residues: set[int] = set()
    for ep in all_linear_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)
    for ep in consensus_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)
    for ep in conformational_epitopes:
        for pos in range(ep.start, ep.end):
            epitope_residues.add(pos)

    coverage = len(epitope_residues) / n if n > 0 else 0.0

    methods_used = list(valid_methods)
    if pdb_string:
        methods_used.append("conformational")
    if consensus_epitopes:
        methods_used.append("consensus")

    return EpitopePredictionResult(
        protein=protein,
        linear_epitopes=final_linear,
        conformational_epitopes=conformational_epitopes,
        per_residue_score=per_residue_score,
        epitope_coverage=round(coverage, 4),
        methods_used=methods_used,
    )
