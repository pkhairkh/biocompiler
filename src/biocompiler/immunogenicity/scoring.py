"""Combined immunogenicity scoring, deimmunization, and batch API.

Top-level ``compute_immunogenicity`` aggregator, epitope-density scoring,
deimmunization mutation suggestions, and the batch computation API.

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
from ._pssm import _get_mhc_i_pssms, _get_mhc_ii_pssms  # noqa: F401
from ._supertypes import _get_supertype_pssm  # noqa: F401
from ._models import *  # noqa: F401,F403
from .mhc_binding import *  # noqa: F401,F403
from .epitopes import *  # noqa: F401,F403
from .epitopes import _validate_protein  # noqa: F401
from .conformational import *  # noqa: F401,F403
from .conformational import _check_real_binding_data_available  # noqa: F401

def compute_immunogenicity(
    protein: str,
    mhc_alleles: list[str] | str | None = None,
    organism: str = "Homo_sapiens",
    is_self_protein: bool = False,
    use_real_data: bool = False,
) -> ImmunogenicityResult:
    """Compute combined immunogenicity score for a protein.

    .. versionchanged:: TIGHTEN-4
       Added an **honesty mode**.  The default PSSMs in this module are
       hand-crafted approximations (NOT real binding data); consequently
       the returned :class:`ImmunogenicityResult` now carries a ``verdict``
       field that is :attr:`~biocompiler.shared.types.Verdict.UNCERTAIN`
       by default, with ``reason="fabricated_scores"`` and
       ``data_source="guessed_pssm"``.  PASS / FAIL verdicts are only
       produced when ``use_real_data=True`` AND a real binding-data
       predictor (NetMHCpan binary on PATH, or MHCflurry with downloaded
       models) is available.

    Runs both T-cell and B-cell epitope prediction and combines
    the results into a single score.

    Scoring formula::

        overall = 0.6 * t_cell_score + 0.4 * b_cell_score

    where:
    - t_cell_score = max epitope score (capped at 1.0)
    - b_cell_score = epitope coverage (fraction of residues in predicted epitopes)

    Classification:
    - low:      overall < 0.3
    - moderate: 0.3 <= overall < 0.6
    - high:     overall >= 0.6

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    mhc_alleles : list[str] | str | None
        MHC alleles for T-cell epitope prediction.  For backward
        compatibility, a bare string is accepted and interpreted as
        the *organism* (so ``compute_immunogenicity(seq, "human")``
        works as a shorthand for
        ``compute_immunogenicity(seq, organism="human")``).
    organism : str
        Host organism (default ``"Homo_sapiens"``).
    is_self_protein : bool
        If True, the protein is from the host organism and
        immunogenicity assessment is not applicable.  Returns a
        result with classification ``"not_applicable"`` and a
        score of 0.0 instead of flagging epitopes.
    use_real_data : bool
        If True, require a real binding-data predictor (NetMHCpan or
        MHCflurry).  If none is available, the result is returned
        with ``verdict=UNCERTAIN`` and ``reason="no_real_predictor"``.
        If False (default), the guessed PSSMs are used to compute
        reference scores, but the verdict is UNCERTAIN with
        ``reason="fabricated_scores"`` — never PASS/FAIL.

    Returns
    -------
    ImmunogenicityResult
        The result carries honesty fields (``verdict``, ``reason``,
        ``message``, ``data_source``, ``scores``) in addition to the
        traditional score / classification fields.
    """
    # Backward-compatible dispatch: allow compute_immunogenicity(seq, "human")
    # by detecting a bare string in the mhc_alleles slot.
    if isinstance(mhc_alleles, str):
        organism = mhc_alleles
        mhc_alleles = None
    # Self-protein short-circuit: skip assessment entirely
    if is_self_protein:
        logger.info(
            "Protein marked as self-protein — skipping immunogenicity assessment"
        )
        return ImmunogenicityResult(
            sequence=protein if protein else "",
            primary_score=0.0,
            classification="not_applicable",
            t_cell_score=0.0,
            b_cell_score=0.0,
            t_cell_epitopes=[],
            b_cell_epitopes=[],
            mutations=[],
            success=True,
            error=None,
            execution_time_s=0.0,
            # Self-protein: no assessment needed.  Verdict is PASS but
            # the data_source remains "guessed_pssm" by default because
            # no real binding prediction was actually performed.
            verdict=Verdict.PASS,
            reason="self_protein_no_assessment",
            message=(
                "Self-protein — immunogenicity assessment not applicable. "
                "Verdict PASS does not rely on any binding-data predictor."
            ),
            data_source="self_protein",
            scores={"overall": 0.0, "t_cell": 0.0, "b_cell": 0.0},
        )

    try:
        protein = _validate_protein(protein)
    except ImmunogenicityError as exc:
        return ImmunogenicityResult(
            sequence=protein if protein else "",
            primary_score=0.0,
            classification="low",
            t_cell_score=0.0,
            b_cell_score=0.0,
            t_cell_epitopes=[],
            b_cell_epitopes=[],
            mutations=[],
            success=False,
            error=str(exc),
            execution_time_s=0.0,
            # Validation failure — cannot make any claim.
            verdict=Verdict.UNCERTAIN,
            reason="invalid_protein",
            message=(
                f"Protein validation failed: {exc}. "
                "No immunogenicity claim is made."
            ),
            data_source="none",
            scores={"overall": 0.0, "t_cell": 0.0, "b_cell": 0.0},
        )

    with EngineTimer() as timer:
        # T-cell prediction
        # NOTE (H18 honesty fix): ``predict_t_cell_epitopes`` uses
        # PSSM-based scoring ONLY — it does NOT accept a
        # ``use_real_data`` flag and does NOT transparently delegate
        # to NetMHCpan / MHCflurry.  Consequently the epitopes /
        # scores computed here are ALWAYS PSSM-derived
        # ("fabricated"), regardless of the ``use_real_data``
        # argument below.  The honesty verdict block further down
        # must NOT issue a PASS/FAIL claim on the basis of these
        # scores, and must NOT label the verdict ``real_data_*``
        # (which would imply real binding data was used).
        t_epitopes = predict_t_cell_epitopes(protein, mhc_alleles)
        if t_epitopes:
            # Use average of top-5 moderate+ binders, not max of all predictions.
            # The previous max-of-all approach was over-sensitive: for any protein
            # >20aa, some 9-mer with good anchor residues would score 0.6-0.8,
            # making t_cell_score ≈ 0.7-1.0 for virtually ALL proteins.  This
            # classified everything as "high immunogenicity" and rendered L6
            # useless.  The average-of-top-5-binders approach is more realistic:
            # only proteins with multiple strong epitopes score high.
            binder_epitopes = [
                e for e in t_epitopes
                if e.get("binding_class") in ("strong_binder", "moderate_binder")
            ]
            if binder_epitopes:
                top_scores = sorted(
                    [e["score"] for e in binder_epitopes], reverse=True
                )[:5]
                t_cell_score = min(1.0, sum(top_scores) / len(top_scores))
            else:
                t_cell_score = 0.0
        else:
            t_cell_score = 0.0

        # Epitope density bonus
        density_score = compute_epitope_density(t_epitopes)
        t_cell_score = min(1.0, t_cell_score + density_score)

        # B-cell prediction using epitope.py's predict_epitopes
        b_result = predict_epitopes(protein)
        b_epitopes_converted: list[BCellEpitopeDict] = [
            BCellEpitopeDict(
                start=ep.start,
                end=ep.end,
                peptide=ep.peptide,
                score=ep.score,
                method=ep.method,
            )
            for ep in b_result.linear_epitopes
        ]
        b_cell_score = b_result.epitope_coverage

        # Combined score
        overall_score = T_CELL_WEIGHT * t_cell_score + B_CELL_WEIGHT * b_cell_score
        overall_score = max(0.0, min(1.0, overall_score))

        # Classification
        if overall_score < IMMUNOGENICITY_LOW_THRESHOLD:
            immuno_class = "low"
        elif overall_score < IMMUNOGENICITY_HIGH_THRESHOLD:
            immuno_class = "moderate"
        else:
            immuno_class = "high"

        # Deimmunization candidates
        deimm_candidates = find_deimmunization_mutations(protein)

        # ── Honesty verdict (TIGHTEN-4) ───────────────────────────
        # The default PSSMs in this module are hand-crafted approximations
        # ("guessed/approximate scores, NOT scores derived from real binding
        # data" — see _build_mhc_i_pssms / _build_mhc_ii_pssms).  We must
        # NOT issue a PASS/FAIL claim on the basis of fabricated scores.
        real_data_available = _check_real_binding_data_available()
        scores_dict = {
            "overall": round(overall_score, 4),
            "t_cell": round(t_cell_score, 4),
            "b_cell": round(b_cell_score, 4),
            "t_cell_epitope_count": len(t_epitopes),
            "b_cell_epitope_count": len(b_epitopes_converted),
        }

        if use_real_data and real_data_available:
            # H18 honesty fix: a real binding-data predictor
            # (NetMHCpan / MHCflurry) IS installed, but the T-cell
            # epitope prediction path above (predict_t_cell_epitopes)
            # uses PSSM-based scoring ONLY — it does NOT accept a
            # ``use_real_data`` flag and does NOT delegate to the
            # real predictor.  The overall_score is therefore STILL
            # PSSM-derived ("fabricated"), and we MUST NOT issue a
            # PASS/FAIL verdict on its basis, nor label the verdict
            # ``real_data_*`` (which would imply real binding data was
            # used).  The verdict is UNCERTAIN; the reason clearly
            # states that the real predictor is available but unused
            # by the PSSM-based epitope path.
            verdict = Verdict.UNCERTAIN
            reason = "pssm_scores_real_predictor_unused"
            message = (
                "use_real_data=True and a real binding-data predictor "
                "(NetMHCpan / MHCflurry) is installed, but the T-cell "
                "epitope prediction path uses PSSM-based scoring only "
                "(predict_t_cell_epitopes does not delegate to the real "
                "predictor). Scores are therefore still PSSM-derived. "
                "Verdict is UNCERTAIN — no PASS/FAIL claim is made. "
                "To obtain real binding-data verdicts, call the "
                "NetMHCpan / MHCflurry adapter directly."
            )
            data_source = "guessed_pssm"
        elif use_real_data and not real_data_available:
            # Caller asked for real data but none is available.
            verdict = Verdict.UNCERTAIN
            reason = "no_real_predictor"
            message = (
                "use_real_data=True but no real MHC binding-data predictor "
                "is installed (NetMHCpan / MHCflurry missing). Verdict is "
                "UNCERTAIN. The guessed-PSSM scores below are provided "
                "for reference only — they are NOT real binding data."
            )
            data_source = "guessed_pssm_no_real_available"
        else:
            # Default: guessed PSSMs only.  Verdict MUST be UNCERTAIN.
            verdict = Verdict.UNCERTAIN
            reason = "fabricated_scores"
            message = (
                "Immunogenicity scores use approximate PSSMs, NOT real "
                "binding data. Verdict is UNCERTAIN. Install NetMHCpan "
                "for verified predictions."
            )
            data_source = "guessed_pssm"

        return ImmunogenicityResult(
            sequence=protein,
            primary_score=round(overall_score, 4),
            classification=immuno_class,
            t_cell_score=round(t_cell_score, 4),
            b_cell_score=round(b_cell_score, 4),
            t_cell_epitopes=t_epitopes,
            b_cell_epitopes=b_epitopes_converted,
            mutations=deimm_candidates,
            execution_time_s=round(timer.elapsed, 4),
            verdict=verdict,
            reason=reason,
            message=message,
            data_source=data_source,
            scores=scores_dict,
        )


# ═══════════════════════════════════════════════════════════════════════════
# Epitope density scoring
# ═══════════════════════════════════════════════════════════════════════════

#: Maximum distance (in amino acids) between epitope centres for them to
#: be considered part of the same cluster.
EPITOPE_DENSITY_CLUSTER_DISTANCE: int = 15


def compute_epitope_density(
    epitopes: list[TCellEpitopeDict],
    cluster_distance: int = EPITOPE_DENSITY_CLUSTER_DISTANCE,
) -> float:
    """Compute epitope density score based on clustering nearby epitopes.

    Epitopes within *cluster_distance* amino acids of each other are
    clustered together.  The density score is scaled by the **total
    number of clustered epitopes** across all clusters (not just the
    largest cluster), providing a more granular assessment of
    immunogenicity risk.  A bonus is applied for clusters containing
    multiple strong binders.

    Parameters
    ----------
    epitopes : list[TCellEpitopeDict]
        Predicted T-cell epitopes (each with ``start``, ``end``, ``score``,
        ``allele``, ``binding_class``).
    cluster_distance : int
        Maximum distance (aa) between epitope centres for clustering
        (default 15).

    Returns
    -------
    float
        Density bonus score in [0, 0.3].  A higher score indicates
        a denser epitope region, suggesting stronger immunogenicity.
    """
    if not epitopes:
        return 0.0

    # Compute centre position for each epitope
    centres: list[tuple[int, TCellEpitopeDict]] = []
    for epi in epitopes:
        centre = (epi["start"] + epi["end"]) // 2
        centres.append((centre, epi))

    # Sort by centre position
    centres.sort(key=lambda x: x[0])

    # Cluster epitopes within cluster_distance of each other
    clusters: list[list[TCellEpitopeDict]] = []
    if not centres:
        return 0.0

    current_cluster: list[TCellEpitopeDict] = [centres[0][1]]
    last_centre = centres[0][0]

    for centre, epi in centres[1:]:
        if centre - last_centre <= cluster_distance:
            current_cluster.append(epi)
        else:
            if len(current_cluster) >= 2:
                clusters.append(current_cluster)
            current_cluster = [epi]
        last_centre = centre

    if len(current_cluster) >= 2:
        clusters.append(current_cluster)

    if not clusters:
        return 0.0

    # Total number of epitopes across all clusters (granular scaling)
    total_clustered = sum(len(c) for c in clusters)
    max_cluster_size = max(len(c) for c in clusters)
    max_strong_in_cluster = max(
        sum(1 for e in c if e.get("binding_class") == "strong_binder")
        for c in clusters
    )

    # Density score scales with the total number of clustered epitopes,
    # weighted by the largest cluster proportion.  This provides a more
    # granular score than a flat bonus — e.g. 3 clusters of 2 epitopes
    # each scores higher than 1 cluster of 2.
    density_score = min(0.2, max_cluster_size * 0.04) * (1 + 0.05 * (total_clustered - max_cluster_size))
    density_score = min(0.2, density_score)

    # Strong-binder bonus also scales with the number of strong binders
    # across all clusters
    total_strong = sum(
        sum(1 for e in c if e.get("binding_class") == "strong_binder")
        for c in clusters
    )
    strong_bonus = min(0.1, max_strong_in_cluster * 0.05) * (1 + 0.05 * max(0, total_strong - max_strong_in_cluster))
    strong_bonus = min(0.1, strong_bonus)

    return round(density_score + strong_bonus, 4)


# ═══════════════════════════════════════════════════════════════════════════
# Deimmunization mutation finding
# ═══════════════════════════════════════════════════════════════════════════


def suggest_mutations(
    protein: str,
    epitope_threshold: float = 0.5,
    blosum62_min: int = 1,
    max_suggestions: int = 20,
    protected_positions: set[int] | None = None,
) -> list[MutationResult]:
    """Suggest conservative mutations to reduce immunogenicity.

    Unlike :func:`find_deimmunization_mutations`, which evaluates all
    19 possible substitutions at each epitope position, this function
    only proposes **conservative** substitutions — those with a
    BLOSUM62 score >= *blosum62_min* (default 1), meaning the mutant
    amino acid is chemically similar to the wildtype.  This is useful
    when protein function must be preserved while reducing
    immunogenicity.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    epitope_threshold : float
        Only consider epitopes with score > this threshold (default 0.5).
    blosum62_min : int
        Minimum BLOSUM62 substitution score for a mutation to be
        considered conservative (default 1, i.e. only similar AAs).
    max_suggestions : int
        Maximum number of mutation suggestions to return (default 20).
    protected_positions : set[int] | None
        Set of 0-based residue positions that should **not** be mutated
        (e.g. catalytic residues, binding-site residues, disulfide
        cysteines).  If ``None``, no positions are protected.

    Returns
    -------
    list[MutationResult]
        Conservative mutation suggestions sorted by largest binding
        score reduction, then by highest BLOSUM62 score.
    """
    protein = _validate_protein(protein)

    # Normalise protected_positions to an empty set when not provided
    _protected: set[int] = protected_positions if protected_positions is not None else set()

    # Get T-cell epitopes
    t_epitopes = predict_t_cell_epitopes(protein)

    # Filter to epitopes above threshold
    strong_epitopes = [
        e for e in t_epitopes if e["score"] > epitope_threshold
    ]

    if not strong_epitopes:
        return []

    # Track which (position, allele) combos we have already scored
    seen: set[tuple[int, str]] = set()
    candidates: list[MutationResult] = []

    mhc_i_pssms = _get_mhc_i_pssms()
    mhc_ii_pssms = _get_mhc_ii_pssms()

    for epi in strong_epitopes:
        allele = epi["allele"]
        original_score = epi["score"]
        start = epi["start"]
        end = epi["end"]  # exclusive
        peptide = epi["peptide"]

        for pos in range(start, end):
            if pos >= len(protein):
                continue
            # Skip protected positions (active sites, catalytic residues, etc.)
            if pos in _protected:
                continue
            key = (pos, allele)
            if key in seen:
                continue
            seen.add(key)

            wildtype = protein[pos]

            for mutant in sorted(_STANDARD_AA_SET):
                if mutant == wildtype:
                    continue

                # Only consider conservative substitutions
                blosum_score = BLOSUM62.get(wildtype, {}).get(mutant, -10)
                if blosum_score < blosum62_min:
                    continue

                # Build mutated protein and re-score
                mutated_protein = protein[:pos] + mutant + protein[pos + 1 :]

                if allele in mhc_i_pssms:
                    if end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:end]
                        new_score = score_peptide_pssm(new_peptide, allele)
                    else:
                        new_score = original_score
                elif allele in mhc_ii_pssms:
                    mhc_ii_window = 15
                    pep_end = start + mhc_ii_window
                    if pep_end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:pep_end]
                        new_score = _score_peptide_for_allele(new_peptide, allele)
                    else:
                        new_score = original_score
                else:
                    # Try supertype fallback PSSM
                    fallback_pssm = _get_supertype_pssm(allele)
                    if fallback_pssm is not None and end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:end]
                        if len(new_peptide) == len(fallback_pssm):
                            new_score = score_peptide_pssm(new_peptide, fallback_pssm)
                        else:
                            new_score = original_score
                    else:
                        new_score = original_score

                score_change = new_score - original_score

                if score_change < 0:
                    candidates.append(
                        MutationResult(
                            position=pos,
                            original=wildtype,
                            mutant=mutant,
                            delta_score=round(-score_change, 4),
                            score_type="immunogenicity",
                            engine="immunogenicity",
                            recommendation="deimmunizing_conservative",
                            description=(
                                f"{wildtype}{pos+1}{mutant}: conservatively reduces "
                                f"{allele} binding by {abs(score_change):.4f} "
                                f"(BLOSUM62={blosum_score})"
                            ),
                            details={
                                "epitope": peptide,
                                "binding_score_change": round(score_change, 4),
                                "blosum62": blosum_score,
                                "conservative": True,
                                "allele": allele,
                            },
                        )
                    )

    # Sort by largest improvement then highest BLOSUM62
    candidates.sort(key=lambda c: (-c.delta_score, -c.details.get("blosum62", 0)))

    return candidates[:max_suggestions]


def find_deimmunization_mutations(
    protein: str,
    epitope_threshold: float = 0.7,
    blosum62_min: int = 0,
    organism: str = "Homo_sapiens",
) -> list[MutationResult]:
    """Find mutations that may reduce immunogenicity.

    For each T-cell epitope scoring above *epitope_threshold*,
    considers every position within the epitope and evaluates
    all 19 possible substitutions.  A substitution that reduces
    the epitope binding score and satisfies the BLOSUM62
    conservation threshold is returned as a candidate.

    Parameters
    ----------
    protein : str
        Amino-acid sequence.
    epitope_threshold : float
        Only consider epitopes with score > this threshold.
    blosum62_min : int
        Minimum BLOSUM62 substitution score (conservative mutations).

    Returns
    -------
    list[MutationResult]
        Mutation suggestions that reduce immunogenicity, sorted by
        largest binding score reduction.
    """
    protein = _validate_protein(protein)

    # Get T-cell epitopes using PSSM-based scoring
    t_epitopes = predict_t_cell_epitopes(protein)

    # Filter to strong epitopes
    strong_epitopes = [
        e for e in t_epitopes if e["score"] > epitope_threshold
    ]

    if not strong_epitopes:
        return []

    # Deduplicate: track which (position, allele) combos we have already scored
    seen: set[tuple[int, str]] = set()
    candidates: list[MutationResult] = []

    mhc_i_pssms = _get_mhc_i_pssms()
    mhc_ii_pssms = _get_mhc_ii_pssms()

    for epi in strong_epitopes:
        allele = epi["allele"]
        original_score = epi["score"]
        start = epi["start"]
        end = epi["end"]  # exclusive
        peptide = epi["peptide"]

        for pos in range(start, end):
            if pos >= len(protein):
                continue
            key = (pos, allele)
            if key in seen:
                continue
            seen.add(key)

            wildtype = protein[pos]

            for mutant in sorted(_STANDARD_AA_SET):
                if mutant == wildtype:
                    continue

                # Check BLOSUM62 conservation
                blosum_score = BLOSUM62.get(wildtype, {}).get(mutant, -10)
                if blosum_score < blosum62_min:
                    continue

                # Build mutated protein and re-score the epitope region
                mutated_protein = protein[:pos] + mutant + protein[pos + 1 :]

                # Re-score using PSSM-based scoring
                if allele in mhc_i_pssms:
                    if end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:end]
                        new_score = score_peptide_pssm(new_peptide, allele)
                    else:
                        new_score = original_score
                elif allele in mhc_ii_pssms:
                    # MHC-II: 15-mer with 9-mer core scanning
                    mhc_ii_window = 15
                    pep_end = start + mhc_ii_window
                    if pep_end <= len(mutated_protein):
                        new_peptide = mutated_protein[start:pep_end]
                        new_score = _score_peptide_for_allele(new_peptide, allele)
                    else:
                        new_score = original_score
                else:
                    new_score = original_score

                score_change = new_score - original_score

                # Only keep substitutions that reduce binding
                if score_change < 0:
                    candidates.append(
                        MutationResult(
                            position=pos,
                            original=wildtype,
                            mutant=mutant,
                            delta_score=round(-score_change, 4),  # positive = improvement
                            score_type="immunogenicity",
                            engine="immunogenicity",
                            recommendation="deimmunizing",
                            description=(
                                f"{wildtype}{pos+1}{mutant}: reduces {allele} "
                                f"binding by {abs(score_change):.4f}"
                            ),
                            details={
                                "epitope": peptide,
                                "binding_score_change": round(score_change, 4),
                                "blosum62": blosum_score,
                                "protein_preserved": blosum_score >= 0,
                                "allele": allele,
                            },
                        )
                    )

    # Sort by largest improvement (highest delta_score), then by BLOSUM62
    candidates.sort(key=lambda c: (-c.delta_score, -c.details.get("blosum62", 0)))

    # Limit to top candidates
    return candidates[:MAX_DEIMMUNIZATION_CANDIDATES]


# ═══════════════════════════════════════════════════════════════════════════
# Batch API
# ═══════════════════════════════════════════════════════════════════════════


def compute_immunogenicity_batch(
    sequences: list[str],
    max_workers: int | None = None,
    **kwargs,
) -> BatchResult[ImmunogenicityResult]:
    """Compute immunogenicity scores for multiple sequences in parallel.

    Uses ``concurrent.futures.ThreadPoolExecutor`` for parallelism.

    Parameters
    ----------
    sequences : list[str]
        List of protein amino-acid sequences.
    **kwargs
        Additional keyword arguments passed to :func:`compute_immunogenicity`
        (e.g. ``mhc_alleles``).

    Returns
    -------
    BatchResult[ImmunogenicityResult]
        Batch result containing one result per input sequence, in the same order.
    """
    logger.info("compute_immunogenicity_batch: processing %d sequences", len(sequences))

    results: list[ImmunogenicityResult] = []
    errors: list[str] = []

    with EngineTimer() as batch_timer:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_idx = {
                executor.submit(compute_immunogenicity, seq, **kwargs): i
                for i, seq in enumerate(sequences)
            }
            result_map: dict[int, ImmunogenicityResult] = {}
            for future in concurrent.futures.as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    result_map[idx] = future.result()
                except Exception as exc:
                    logger.error(
                        "compute_immunogenicity_batch: sequence %d failed: %s",
                        idx, exc,
                    )
                    result_map[idx] = ImmunogenicityResult(
                        sequence=sequences[idx] if sequences[idx] else "",
                        primary_score=0.0,
                        classification="low",
                        success=False,
                        error=str(exc),
                    )
                    errors.append(f"sequence {idx}: {exc}")

            for i in range(len(sequences)):
                results.append(result_map[i])

    logger.info(
        "compute_immunogenicity_batch: completed %d/%d successfully",
        sum(1 for r in results if r.success),
        len(results),
    )

    return BatchResult[ImmunogenicityResult](
        results=results,
        errors=errors,
        total_time_s=round(batch_timer.elapsed, 4),
    )

__all__ = [
    "compute_immunogenicity",
    "compute_epitope_density",
    "suggest_mutations",
    "find_deimmunization_mutations",
    "compute_immunogenicity_batch",
    "EPITOPE_DENSITY_CLUSTER_DISTANCE",
]
