"""MHC binding prediction.

PSSM scoring, IC50 mapping, binding classification, and allele-specific
peptide-MHC binding prediction for MHC-I and MHC-II.

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
from ._constants import _STANDARD_AA_SET  # noqa: F401
from ._pssm import *  # noqa: F401,F403
from ._pssm import _get_mhc_i_pssms, _get_mhc_ii_pssms  # noqa: F401
from ._supertypes import _get_supertype_pssm  # noqa: F401
from ._models import *  # noqa: F401,F403

# ═══════════════════════════════════════════════════════════════════════════
# Prediction cache
# ═══════════════════════════════════════════════════════════════════════════

_prediction_cache: dict[tuple[str, str, int], list[MHCBindingResult]] = {}


def clear_cache() -> None:
    """Clear the MHC binding prediction cache.

    PSSMs are built once and kept across cache clears.
    """
    global _prediction_cache
    _prediction_cache.clear()
    logger.info("Immunogenicity prediction cache cleared")
def score_peptide_pssm(
    peptide: str,
    pssm: list[dict[str, float]] | str,
) -> float:
    """Compute binding score from a PSSM.

    Uses the geometric mean of position-specific scores and normalises
    to the [0, 1] range.

    Parameters
    ----------
    peptide : str
        Amino-acid sequence of the peptide (must equal PSSM length).
    pssm : list[dict[str, float]] or str
        Position-specific scoring matrix, one dict per position.
        Alternatively, an allele name string (e.g. ``"HLA-A*02:01"``)
        which will be looked up in :data:`MHC_I_PSSM` and
        :data:`MHC_II_PSSM`.

    Returns
    -------
    float
        Normalised binding score in [0, 1].
    """
    if isinstance(pssm, str):
        allele = pssm
        lookup = _get_mhc_i_pssms().get(allele) or _get_mhc_ii_pssms().get(allele)
        if lookup is None:
            logger.debug("No PSSM for allele %s — returning 0.0", allele)
            return 0.0
        pssm = lookup

    if len(peptide) != len(pssm):
        logger.warning(
            "Peptide length %d does not match PSSM length %d — returning 0.0",
            len(peptide),
            len(pssm),
        )
        return 0.0

    scores: list[float] = []
    for i, aa in enumerate(peptide):
        aa_upper = aa.upper()
        if aa_upper not in pssm[i]:
            scores.append(PSSM_UNKNOWN_AA_SCORE)
        else:
            scores.append(pssm[i][aa_upper])

    log_sum = sum(math.log(max(s, 1e-10)) for s in scores)
    geo_mean = math.exp(log_sum / len(scores))

    max_scores: list[float] = []
    min_scores: list[float] = []
    for pos_dict in pssm:
        vals = list(pos_dict.values())
        max_scores.append(max(vals))
        min_scores.append(min(vals))
    max_log_sum = sum(math.log(max(s, 1e-10)) for s in max_scores)
    max_geo_mean = math.exp(max_log_sum / len(max_scores))
    min_log_sum = sum(math.log(max(s, 1e-10)) for s in min_scores)
    min_geo_mean = math.exp(min_log_sum / len(min_scores))

    if max_geo_mean <= min_geo_mean:
        return 0.0

    raw = (geo_mean - min_geo_mean) / (max_geo_mean - min_geo_mean)
    raw = max(0.0, min(1.0, raw))

    normalised = raw ** PSSM_CONTRAST_POWER

    return max(0.0, min(1.0, normalised))


def binding_score_to_ic50(score: float, *, from_pssm: bool = True) -> float:
    """Map a binding score to an estimated IC50 (nM) using a log-linear mapping.

    Parameters
    ----------
    score : float
        Normalised binding score in [0, 1].
    from_pssm : bool
        Whether the score was derived from a PSSM heuristic.  If True,
        a DeprecationWarning is emitted because PSSM-derived IC50
        estimates are inaccurate and will be removed in a future version.
        Use MHCflurry or NetMHCpan for reliable IC50 predictions.

    Returns
    -------
    float
        Estimated IC50 in nM.

    Notes
    -----
    Effective formula: IC50 = 10 ** (IC50_LOG_INTERCEPT - IC50_LOG_SLOPE * score), calibrated so:
      - score ~0.9 -> ~50 nM (strong)
      - score ~0.5 -> ~500 nM (moderate)
      - score ~0.1 -> ~5000 nM (weak)
    """
    if from_pssm:
        warnings.warn(
            "binding_score_to_ic50() with PSSM-derived scores is deprecated. "
            "PSSM-based IC50 estimates are inaccurate. Use MHCflurry or "
            "NetMHCpan for reliable IC50 predictions. This function will "
            "be removed in a future version.",
            DeprecationWarning,
            stacklevel=2,
        )
    clamped = max(0.0, min(1.0, score))
    return 10.0 ** (IC50_LOG_INTERCEPT - IC50_LOG_SLOPE * clamped)


def classify_binding(ic50: float) -> str:
    """Classify a peptide by its IC50 value.

    Parameters
    ----------
    ic50 : float
        IC50 in nM.

    Returns
    -------
    str
        One of ``"strong_binder"``, ``"moderate_binder"``,
        ``"weak_binder"``, ``"non_binder"``.
    """
    if ic50 < IC50_STRONG_BINDER_THRESHOLD:
        return "strong_binder"
    elif ic50 <= IC50_MODERATE_BINDER_THRESHOLD:
        return "moderate_binder"
    elif ic50 <= IC50_WEAK_BINDER_THRESHOLD:
        return "weak_binder"
    else:
        return "non_binder"




# ═══════════════════════════════════════════════════════════════════════════
# Offline prediction API
def _is_mhc_ii_allele(allele: str) -> bool:
    """Return True if the allele is an MHC class II allele."""
    return allele.startswith("HLA-DR") or allele.startswith("HLA-DQ") or allele.startswith("HLA-DP")


def predict_immunogenicity(allele: str, peptide: str) -> ImmunogenicityPrediction:
    """Predict immunogenicity for a single peptide-allele pair, offline.

    Prediction hierarchy (never requires external tools):
    1. Precomputed database lookup
    2. PSSM scoring fallback

    Parameters
    ----------
    allele : str
        MHC allele name.
    peptide : str
        Amino-acid sequence of the peptide.

    Returns
    -------
    ImmunogenicityPrediction
    """
    # Step 1: precomputed lookup
    lookup_key = (allele, peptide.upper())
    if lookup_key in PRECOMPUTED_BINDERS:
        ic50 = PRECOMPUTED_BINDERS[lookup_key]
        return ImmunogenicityPrediction(
            allele=allele, peptide=peptide, ic50_nm=ic50,
            binding_class=classify_binding(ic50),
            method="precomputed_lookup", confidence="high",
        )

    # Step 2: PSSM scoring
    mhc_i_pssms = _get_mhc_i_pssms()
    mhc_ii_pssms = _get_mhc_ii_pssms()

    if allele in mhc_i_pssms:
        pssm = mhc_i_pssms[allele]
        if len(peptide) == len(pssm):
            score = score_peptide_pssm(peptide, pssm)
            ic50 = binding_score_to_ic50(score)
            anchor_res, anchor_scores_val = _identify_anchor_positions(peptide, pssm)
            if anchor_res and all(v >= 1.0 for v in anchor_scores_val.values()):
                confidence = "medium"
            else:
                confidence = "low"
            return ImmunogenicityPrediction(
                allele=allele, peptide=peptide,
                ic50_nm=round(ic50, 2), binding_class=classify_binding(ic50),
                method="pssm_fallback", confidence=confidence,
            )
        else:
            return ImmunogenicityPrediction(
                allele=allele, peptide=peptide, ic50_nm=50000.0,
                binding_class="non_binder", method="pssm_fallback", confidence="low",
            )

    if allele in mhc_ii_pssms:
        pssm = mhc_ii_pssms[allele]
        if len(peptide) >= MHC_II_CORE_LENGTH:
            best_score = 0.0
            for offset in range(len(peptide) - MHC_II_CORE_LENGTH + 1):
                core = peptide[offset : offset + MHC_II_CORE_LENGTH]
                s = score_peptide_pssm(core, pssm)
                best_score = max(best_score, s)
            ic50 = binding_score_to_ic50(best_score)
            return ImmunogenicityPrediction(
                allele=allele, peptide=peptide,
                ic50_nm=round(ic50, 2), binding_class=classify_binding(ic50),
                method="pssm_fallback", confidence="low",
            )
        else:
            return ImmunogenicityPrediction(
                allele=allele, peptide=peptide, ic50_nm=50000.0,
                binding_class="non_binder", method="pssm_fallback", confidence="low",
            )

    # Unknown allele
    return ImmunogenicityPrediction(
        allele=allele, peptide=peptide, ic50_nm=50000.0,
        binding_class="non_binder", method="pssm_fallback", confidence="low",
    )


def scan_peptides(protein: str, allele: str, peptide_length: int = 9) -> list[PeptideResult]:
    """Generate all overlapping peptides and score each against an allele.

    Works entirely offline using precomputed database lookups and PSSM
    fallback - never requires MHCflurry or NetMHCpan.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    allele : str
        MHC allele name.
    peptide_length : int
        Length of the sliding window (default 9 for MHC-I).

    Returns
    -------
    list[PeptideResult]
        Peptide results sorted by binding strength (lowest IC50 first).
    """
    if not protein or peptide_length < 1:
        return []
    protein = protein.upper()
    if len(protein) < peptide_length:
        return []
    results: list[PeptideResult] = []
    for start in range(len(protein) - peptide_length + 1):
        peptide = protein[start : start + peptide_length]
        if any(c not in _STANDARD_AA_SET for c in peptide):
            continue
        pred = predict_immunogenicity(allele, peptide)
        results.append(PeptideResult(
            position=start, peptide=peptide,
            ic50_nm=pred.ic50_nm, binding_class=pred.binding_class,
        ))
    results.sort(key=lambda r: r.ic50_nm)
    return results


def _identify_anchor_positions(
    peptide: str,
    pssm: list[dict[str, float]],
    threshold: float = 2.5,
) -> tuple[dict[int, str], dict[int, float]]:
    """Identify anchor residues in a peptide relative to a PSSM.

    An anchor position is one where the position has high selectivity
    (the ratio of max/min score in the PSSM row exceeds *threshold*).

    Returns
    -------
    anchor_residues : dict[int, str]
        Position index -> amino acid at that anchor position.
    anchor_scores : dict[int, float]
        Position index -> the PSSM score at that position.
    """
    anchor_residues: dict[int, str] = {}
    anchor_scores: dict[int, float] = {}

    for i, aa in enumerate(peptide):
        aa_upper = aa.upper()
        row = pssm[i]
        row_values = list(row.values())
        selectivity = max(row_values) / max(min(row_values), 1e-10)
        if selectivity >= threshold:
            score = row.get(aa_upper, 0.5)
            anchor_residues[i] = aa_upper
            anchor_scores[i] = score

    return anchor_residues, anchor_scores


# ═══════════════════════════════════════════════════════════════════════════
# MHC binding: prediction functions
# ═══════════════════════════════════════════════════════════════════════════


def predict_mhc_i_binding(
    protein: str,
    alleles: list[str] | None = None,
    peptide_length: int = DEFAULT_MHC_PEPTIDE_LENGTH,
    use_netmhcpan: bool = False,
    use_mhcflurry: bool = True,
) -> list[MHCBindingResult]:
    """Predict MHC class I binding for overlapping peptides.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    alleles : list[str] or None
        MHC-I alleles to evaluate. Defaults to :data:`DEFAULT_MHC_I_ALLELES`.
    peptide_length : int
        Length of peptides to extract (default 9).
    use_netmhcpan : bool
        If True, try the NetMHCpan web API first for more accurate
        predictions.  Falls back to MHCflurry or PSSM if unavailable.
        Default False.
    use_mhcflurry : bool
        If True, try MHCflurry as an offline neural-network predictor
        (AUC 0.80-0.85, confidence=0.85).  Used as intermediate tier
        between NetMHCpan (0.85-0.95, confidence=0.95) and PSSM
        (0.60-0.75).  Falls back to PSSM if MHCflurry is not installed.
        Default True.

    Returns
    -------
    list[MHCBindingResult]
        Binding predictions for every peptide x allele combination.

    Prediction hierarchy
    --------------------
    NetMHCpan (online, AUC 0.85-0.95)
      → MHCflurry (offline NN, AUC 0.80-0.85)
        → PSSM (offline heuristic, AUC 0.60-0.75)
    """
    if alleles is None:
        alleles = DEFAULT_MHC_I_ALLELES

    if not protein or peptide_length < 1:
        return []

    # Try NetMHCpan if requested
    if use_netmhcpan:
        try:
            from biocompiler.immunogenicity.netmhcpan import NetMHCpanClient
            client = NetMHCpanClient()
            results = client.batch_predict(
                protein, alleles, epitope_lengths=[peptide_length],
            )
            # Convert netmhcpan MHCBindingResult objects to the
            # immunogenicity module's MHCBindingResult format
            converted = []
            for r in results:
                converted.append(MHCBindingResult(
                    allele=r.allele,
                    peptide=r.peptide,
                    start_position=r.start_position,
                    end_position=r.end_position,
                    binding_score=r.binding_score,
                    ic50_nm=r.ic50_nm,
                    binding_class=r.binding_class,
                    anchor_residues=r.anchor_residues,
                    anchor_scores=r.anchor_scores,
                ))
            logger.info(
                "MHC-I prediction via NetMHCpan: %d results for %d alleles, "
                "protein length %d",
                len(converted), len(alleles), len(protein),
            )
            return converted
        except Exception as exc:
            logger.warning(
                "NetMHCpan API failed, falling back to MHCflurry/PSSM: %s", exc,
            )

    # Try MHCflurry adapter if requested (offline NN predictor, AUC 0.80-0.85)
    # The adapter's fallback chain handles MHCflurry → NetMHCpan →
    # precomputed → PSSM gracefully, so it never crashes.
    if use_mhcflurry:
        try:
            from biocompiler.immunogenicity.mhcflurry_adapter import predict_binding as adapter_predict
            results = []
            for start in range(len(protein) - peptide_length + 1):
                peptide = protein[start : start + peptide_length]
                if any(c.upper() not in _STANDARD_AA_SET for c in peptide):
                    continue
                for allele in alleles:
                    r = adapter_predict(allele, peptide)
                    results.append(MHCBindingResult(
                        allele=r.allele,
                        peptide=r.peptide,
                        start_position=start,
                        end_position=start + peptide_length - 1,
                        binding_score=r.binding_score,
                        ic50_nm=r.ic50_nm,
                        binding_class=r.binding_class,
                        anchor_residues=r.anchor_residues,
                        anchor_scores=r.anchor_scores,
                        method=r.method,
                        rank=r.rank,
                        confidence=r.confidence,
                    ))
            logger.info(
                "MHC-I prediction via mhcflurry_adapter: %d results for %d alleles, "
                "protein length %d",
                len(results), len(alleles), len(protein),
            )
            return results
        except Exception as exc:
            logger.warning(
                "MHCflurry adapter prediction failed, falling back to PSSM: %s", exc,
            )

    # PSSM-based prediction (original implementation, also serves as fallback)
    # Check cache
    cache_key = (hashlib.sha256(protein.encode()).hexdigest(), ",".join(alleles), peptide_length)
    if cache_key in _prediction_cache:
        return _prediction_cache[cache_key]

    results: list[MHCBindingResult] = []
    mhc_i_pssms = _get_mhc_i_pssms()

    for allele in alleles:
        pssm = mhc_i_pssms.get(allele)
        method = "pssm_fallback"
        confidence = 0.5

        if pssm is None:
            # Try supertype fallback PSSM
            fallback = _get_supertype_pssm(allele)
            if fallback is not None and len(fallback) == peptide_length:
                pssm = fallback
                method = "pssm_fallback"
                confidence = 0.3
                logger.debug(
                    "Using supertype fallback PSSM for allele %s", allele,
                )
            else:
                logger.debug("No PSSM for allele %s — skipping", allele)
                continue

        if len(pssm) != peptide_length:
            logger.debug(
                "PSSM length %d does not match peptide_length %d for %s — skipping",
                len(pssm),
                peptide_length,
                allele,
            )
            continue

        for start in range(len(protein) - peptide_length + 1):
            peptide = protein[start : start + peptide_length]

            if any(c.upper() not in _STANDARD_AA_SET for c in peptide):
                continue

            score = score_peptide_pssm(peptide, pssm)
            ic50 = binding_score_to_ic50(score)
            binding_class = classify_binding(ic50)
            anchor_residues, anchor_scores = _identify_anchor_positions(peptide, pssm)

            results.append(
                MHCBindingResult(
                    allele=allele,
                    peptide=peptide,
                    start_position=start,
                    end_position=start + peptide_length - 1,
                    binding_score=round(score, 6),
                    ic50_nm=round(ic50, 2),
                    binding_class=binding_class,
                    anchor_residues=anchor_residues,
                    anchor_scores={k: round(v, 4) for k, v in anchor_scores.items()},
                    method=method,
                    confidence=confidence,
                )
            )

    _prediction_cache[cache_key] = results

    logger.info(
        "MHC-I prediction: %d results for %d alleles, protein length %d",
        len(results),
        len(alleles),
        len(protein),
    )
    return results


def predict_mhc_ii_binding(
    protein: str,
    alleles: list[str] | None = None,
    peptide_length: int = 15,
    use_netmhcpan: bool = False,
) -> list[MHCBindingResult]:
    """Predict MHC class II binding for overlapping 15-mer peptides.

    MHC-II binding is evaluated by scanning all possible 9-mer core
    registers within each 15-mer peptide and keeping the best-scoring
    core.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    alleles : list[str] or None
        MHC-II alleles to evaluate. Defaults to :data:`DEFAULT_MHC_II_ALLELES`.
    peptide_length : int
        Length of peptides to extract (default 15).
    use_netmhcpan : bool
        If True, try the NetMHCpan web API first for more accurate
        predictions.  Falls back to PSSM if the API is unavailable or
        fails.  Default False (use PSSM heuristic only).

    Returns
    -------
    list[MHCBindingResult]
        Binding predictions for every peptide x allele combination.
    """
    if alleles is None:
        alleles = DEFAULT_MHC_II_ALLELES

    if not protein or peptide_length < 9:
        return []

    # Try NetMHCpan if requested
    if use_netmhcpan:
        try:
            from biocompiler.immunogenicity.netmhcpan import NetMHCpanClient
            client = NetMHCpanClient()
            results = client.batch_predict(
                protein, alleles, epitope_lengths=[peptide_length],
            )
            converted = []
            for r in results:
                converted.append(MHCBindingResult(
                    allele=r.allele,
                    peptide=r.peptide,
                    start_position=r.start_position,
                    end_position=r.end_position,
                    binding_score=r.binding_score,
                    ic50_nm=r.ic50_nm,
                    binding_class=r.binding_class,
                    anchor_residues=r.anchor_residues,
                    anchor_scores=r.anchor_scores,
                ))
            logger.info(
                "MHC-II prediction via NetMHCpan: %d results for %d alleles, "
                "protein length %d",
                len(converted), len(alleles), len(protein),
            )
            return converted
        except Exception as exc:
            logger.warning(
                "NetMHCpan API failed, falling back to PSSM: %s", exc,
            )

    # PSSM-based prediction (original implementation, also serves as fallback)
    # Check cache
    cache_key = (hashlib.sha256(protein.encode()).hexdigest(), ",".join(alleles), -peptide_length)
    if cache_key in _prediction_cache:
        return _prediction_cache[cache_key]

    results: list[MHCBindingResult] = []
    mhc_ii_pssms = _get_mhc_ii_pssms()

    for allele in alleles:
        pssm = mhc_ii_pssms.get(allele)
        method = "pssm_fallback"
        confidence = 0.5

        if pssm is None:
            # Try supertype fallback PSSM for MHC-II
            fallback = _get_supertype_pssm(allele)
            if fallback is not None and len(fallback) == MHC_II_CORE_LENGTH:
                pssm = fallback
                method = "pssm_fallback"
                confidence = 0.3
                logger.debug(
                    "Using supertype fallback PSSM for MHC-II allele %s", allele,
                )
            else:
                logger.debug("No PSSM for allele %s — skipping", allele)
                continue

        if len(pssm) != MHC_II_CORE_LENGTH:
            logger.debug(
                "PSSM length %d != core length %d for %s — skipping",
                len(pssm),
                MHC_II_CORE_LENGTH,
                allele,
            )
            continue

        for start in range(len(protein) - peptide_length + 1):
            peptide = protein[start : start + peptide_length]

            if any(c.upper() not in _STANDARD_AA_SET for c in peptide):
                continue

            best_score = 0.0
            best_core = peptide[:MHC_II_CORE_LENGTH]
            best_core_offset = 0

            for core_start in range(peptide_length - MHC_II_CORE_LENGTH + 1):
                core = peptide[core_start : core_start + MHC_II_CORE_LENGTH]
                score = score_peptide_pssm(core, pssm)
                if score > best_score:
                    best_score = score
                    best_core = core
                    best_core_offset = core_start

            ic50 = binding_score_to_ic50(best_score)
            binding_class = classify_binding(ic50)

            anchor_residues, anchor_scores = _identify_anchor_positions(
                best_core, pssm
            )

            adjusted_anchors: dict[int, str] = {
                k + best_core_offset: v for k, v in anchor_residues.items()
            }
            adjusted_scores: dict[int, float] = {
                k + best_core_offset: round(v, 4)
                for k, v in anchor_scores.items()
            }

            results.append(
                MHCBindingResult(
                    allele=allele,
                    peptide=peptide,
                    start_position=start,
                    end_position=start + peptide_length - 1,
                    binding_score=round(best_score, 6),
                    ic50_nm=round(ic50, 2),
                    binding_class=binding_class,
                    anchor_residues=adjusted_anchors,
                    anchor_scores=adjusted_scores,
                    method=method,
                    confidence=confidence,
                )
            )

    _prediction_cache[cache_key] = results

    logger.info(
        "MHC-II prediction: %d results for %d alleles, protein length %d",
        len(results),
        len(alleles),
        len(protein),
    )
    return results


def predict_all(
    protein: str,
    mhc_i_alleles: list[str] | None = None,
    mhc_ii_alleles: list[str] | None = None,
    use_netmhcpan: bool = False,
) -> MHCPredictionResult:
    """Run both MHC-I and MHC-II predictions and aggregate results.

    Parameters
    ----------
    protein : str
        Full protein amino-acid sequence.
    mhc_i_alleles : list[str] or None
        MHC-I alleles (defaults to :data:`DEFAULT_MHC_I_ALLELES`).
    mhc_ii_alleles : list[str] or None
        MHC-II alleles (defaults to :data:`DEFAULT_MHC_II_ALLELES`).
    use_netmhcpan : bool
        If True, try the NetMHCpan web API first for more accurate
        predictions.  Falls back to PSSM if the API is unavailable or
        fails.  Default False (use PSSM heuristic only).

    Returns
    -------
    MHCPredictionResult
        Aggregated binding prediction.
    """
    with EngineTimer() as timer:
        mhc_i_results = predict_mhc_i_binding(
            protein, alleles=mhc_i_alleles, use_netmhcpan=use_netmhcpan,
        )
        mhc_ii_results = predict_mhc_ii_binding(
            protein, alleles=mhc_ii_alleles, use_netmhcpan=use_netmhcpan,
        )

        all_results = mhc_i_results + mhc_ii_results

        strong_binders = sum(
            1 for r in all_results if r.binding_class == "strong_binder"
        )
        moderate_binders = sum(
            1 for r in all_results if r.binding_class == "moderate_binder"
        )
        weak_binders = sum(
            1 for r in all_results if r.binding_class == "weak_binder"
        )
        non_binders = sum(
            1 for r in all_results if r.binding_class == "non_binder"
        )

        binding_profile: dict[str, float] = {}
        for r in all_results:
            if r.allele not in binding_profile or r.binding_score > binding_profile[r.allele]:
                binding_profile[r.allele] = round(r.binding_score, 6)

        result = MHCPredictionResult(
            protein=protein,
            mhc_i_results=mhc_i_results,
            mhc_ii_results=mhc_ii_results,
            strong_binders=strong_binders + moderate_binders,
            weak_binders=weak_binders,
            non_binders=non_binders,
            binding_profile=binding_profile,
        )

    result.execution_time_s = round(timer.elapsed, 4)

    logger.info(
        "predict_all: %d MHC-I, %d MHC-II results; "
        "strong+moderate=%d, weak=%d, non=%d (%.2fs)",
        len(mhc_i_results),
        len(mhc_ii_results),
        result.strong_binders,
        result.weak_binders,
        result.non_binders,
        result.execution_time_s,
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════


__all__ = [
    "clear_cache",
    "score_peptide_pssm",
    "binding_score_to_ic50",
    "classify_binding",
    "predict_mhc_i_binding",
    "predict_mhc_ii_binding",
    "predict_all",
    "predict_immunogenicity",
    "scan_peptides",
]
