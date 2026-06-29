"""T-cell and linear B-cell epitope prediction.

Classical scale-based methods: Kolaskar-Tongaonkar, Parker hydrophilicity,
Chou-Fasman beta-turn, Emini surface accessibility (EEA), and a
BepiPred-like composite.

Split out of ``core.py`` (W8-a refactor).
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import logging
import math
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, TypedDict

from biocompiler.shared.constants import (
    BLOSUM62,
    DEFAULT_MHC_PEPTIDE_LENGTH,
    HYDROPATHY,
    STANDARD_AAS,
)
from ..engines.base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from biocompiler.shared.exceptions import ImmunogenicityError
from biocompiler.shared.types import Verdict

logger = logging.getLogger(__name__)
from ._constants import *  # noqa: F401,F403
from ._pssm import _get_mhc_i_pssms, _get_mhc_ii_pssms  # noqa: F401
from ._models import *  # noqa: F401,F403
from .mhc_binding import *  # noqa: F401,F403


def _validate_protein(protein: str) -> str:
    """Validate and normalise a protein sequence.

    Uses the shared ``validate_protein_sequence`` from engine_base, then
    raises ``ImmunogenicityError`` on failure instead of ``ValueError``.
    """
    try:
        return validate_protein_sequence(protein, "Immunogenicity")
    except ValueError as exc:
        raise ImmunogenicityError(str(exc)) from exc


def _peptide_hydrophobicity_score(peptide: str) -> float:
    """Score the hydrophobicity of a peptide core (0-1 range).

    Hydrophobic cores favour MHC binding.
    """
    if len(peptide) < 3:
        return 0.0
    core = peptide[1:-1]
    if not core:
        return 0.0
    avg_hydro = sum(HYDROPATHY.get(aa, 0.0) for aa in core) / len(core)
    normalised = (avg_hydro + HYDROPHOBICITY_OFFSET) / HYDROPHOBICITY_RANGE
    return max(0.0, min(1.0, normalised))


def _peptide_charge_score(peptide: str) -> float:
    """Score the charge balance of a peptide (0-1 range).

    A mix of charged and neutral residues favours MHC binding.
    """
    charged = sum(1 for aa in peptide if aa in "DEKRH")
    neutral = len(peptide) - charged
    if len(peptide) == 0:
        return 0.0
    ratio = min(charged, neutral) / max(charged, neutral, 1)
    return min(1.0, ratio * 1.5)


def _score_peptide_for_allele(peptide: str, allele: str) -> float:
    """Score a peptide against an MHC allele using PSSM.

    For MHC-I, the peptide must match the PSSM length.
    For MHC-II, scans all 9-mer cores within the peptide.
    """
    mhc_i = _get_mhc_i_pssms()
    mhc_ii = _get_mhc_ii_pssms()
    if allele in mhc_i:
        return score_peptide_pssm(peptide, allele)
    elif allele in mhc_ii:
        best_score = 0.0
        for offset in range(len(peptide) - 8):
            core = peptide[offset : offset + 9]
            s = score_peptide_pssm(core, allele)
            best_score = max(best_score, s)
        return best_score
    return 0.0


# ═══════════════════════════════════════════════════════════════════════════
# T-cell epitope prediction
# ═══════════════════════════════════════════════════════════════════════════


def predict_t_cell_epitopes(
    protein: str,
    mhc_alleles: list[str] | None = None,
    peptide_length: int = DEFAULT_MHC_PEPTIDE_LENGTH,
) -> list[TCellEpitopeDict]:
    """Predict T-cell epitopes in a protein sequence.

    Uses PSSM-based scoring for MHC-I (9-mers) and MHC-II
    (15-mers with 9-mer core scanning) alleles.

    Parameters
    ----------
    protein : str
        Amino-acid sequence (one-letter codes).
    mhc_alleles : list[str] | None
        MHC alleles to evaluate.  Defaults to all alleles in
        :data:`DEFAULT_MHC_I_ALLELES` + :data:`DEFAULT_MHC_II_ALLELES`.
    peptide_length : int
        Length of the sliding window for MHC-I peptides (default 9).
        MHC-II peptides are always scanned as 15-mers internally.

    Returns
    -------
    list[TCellEpitopeDict]
        Each dict contains: start, end, peptide, score, allele,
        binding_class.
    """
    protein = _validate_protein(protein)

    mhc_i_pssms = _get_mhc_i_pssms()
    mhc_ii_pssms = _get_mhc_ii_pssms()

    if mhc_alleles is not None:
        mhc_i_alleles = [a for a in mhc_alleles if a in mhc_i_pssms]
        mhc_ii_alleles = [a for a in mhc_alleles if a in mhc_ii_pssms]
        unrecognised = set(mhc_alleles) - set(mhc_i_alleles) - set(mhc_ii_alleles)
        for allele in unrecognised:
            logger.warning("Unrecognised MHC allele: %s — skipping", allele)
    else:
        mhc_i_alleles = DEFAULT_MHC_I_ALLELES
        mhc_ii_alleles = DEFAULT_MHC_II_ALLELES

    epitopes: list[TCellEpitopeDict] = []

    # MHC-I predictions
    if mhc_i_alleles:
        mhc_i_results = predict_mhc_i_binding(protein, mhc_i_alleles, peptide_length)
        for r in mhc_i_results:
            epitopes.append(TCellEpitopeDict(
                start=r.start_position,
                end=r.end_position + 1,  # exclusive
                peptide=r.peptide,
                score=round(r.binding_score, 4),
                allele=r.allele,
                binding_class=r.binding_class,
            ))

    # MHC-II predictions
    if mhc_ii_alleles:
        mhc_ii_results = predict_mhc_ii_binding(protein, mhc_ii_alleles)
        for r in mhc_ii_results:
            epitopes.append(TCellEpitopeDict(
                start=r.start_position,
                end=r.end_position + 1,  # exclusive
                peptide=r.peptide,
                score=round(r.binding_score, 4),
                allele=r.allele,
                binding_class=r.binding_class,
            ))

    epitopes.sort(key=lambda e: e["score"], reverse=True)
    return epitopes


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: amino acid scales
# ═══════════════════════════════════════════════════════════════════════════

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
    "A": -1.01, "R": 1.40, "N": 0.81, "D": 1.21,
    "C": -1.20, "E": 1.64, "Q": 0.96, "G": -0.16,
    "H": 0.56, "I": -1.42, "L": -1.42, "K": 1.73,
    "M": -1.27, "F": -1.42, "P": 0.26, "S": 0.52,
    "T": -0.19, "W": -1.07, "Y": -0.31, "V": -1.07,
}

CHOU_FASMAN_TURN: dict[str, float] = {
    # Chou-Fasman beta-turn propensity (Chou & Fasman 1974)
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

# Removed: ANTIGENICITY_PROPENSITY (alias for ANTIGENICITY_SCALE).

# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: data classes
# ═══════════════════════════════════════════════════════════════════════════


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


# ═══════════════════════════════════════════════════════════════════════════
# B-cell epitope: prediction methods
# ═══════════════════════════════════════════════════════════════════════════


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
        threshold: Antigenicity threshold (default 1.0).

    Returns:
        List of EpitopeRegion predictions.
    """
    if not protein:
        return []

    protein = protein.upper()

    raw_scores = [ANTIGENICITY_SCALE.get(aa, 1.0) for aa in protein]
    smoothed = _sliding_window_average(raw_scores, window)
    norm_scores = _normalize_01(smoothed)

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
    """Predict B-cell epitopes using Chou-Fasman beta-turn propensity.

    Chou & Fasman (1974) showed that beta-turns are often surface-exposed
    and correspond to B-cell epitope regions.

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

    Args:
        protein: Amino acid sequence (1-letter codes).

    Returns:
        List of EpitopeRegion predictions for surface-accessible regions.
    """
    if not protein:
        return []

    protein = protein.upper()
    n = len(protein)
    window = 6
    half = window // 2

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

    norm_scores = _normalize_01(raw_scores)

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

    Args:
        protein: Amino acid sequence (1-letter codes).
        window: Sliding window size for smoothing.

    Returns:
        List of EpitopeRegion predictions.
    """
    if not protein:
        return []

    protein = protein.upper()

    hydro_raw = [PARKER_SCALE.get(aa, 0.0) for aa in protein]
    flex_raw = [_FLEXIBILITY_SCALE.get(aa, 0.4) for aa in protein]
    surf_raw = [EMINI_SCALE.get(aa, 0.5) for aa in protein]

    hydro_smooth = _sliding_window_average(hydro_raw, window)
    flex_smooth = _sliding_window_average(flex_raw, window)
    surf_smooth = _sliding_window_average(surf_raw, window)

    hydro_norm = _normalize_01(hydro_smooth)
    flex_norm = _normalize_01(flex_smooth)
    surf_norm = _normalize_01(surf_smooth)

    composite: list[float] = [
        0.4 * h + 0.3 * f + 0.3 * s
        for h, f, s in zip(hydro_norm, flex_norm, surf_norm)
    ]

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



__all__ = [
    "ANTIGENICITY_SCALE",
    "PARKER_SCALE",
    "CHOU_FASMAN_TURN",
    "EMINI_SCALE",
    "ALL_SCALES",
    "predict_t_cell_epitopes",
    "predict_kolaskar_tongaonkar",
    "predict_parker_hydrophilicity",
    "predict_chou_fasman_beta_turn",
    "predict_eea",
    "predict_bepipred_like",
]
