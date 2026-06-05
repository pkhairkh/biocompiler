"""
BioCompiler Immunogenicity Predicates
=======================================

Type-check predicates for immunogenicity assessment.
Each predicate returns a TypeCheckResult with a Verdict.

All predicates work fully offline - they never require MHCflurry or
NetMHCpan.  They use the precomputed database + PSSM fallback via
:func:`immunogenicity.predict_immunogenicity` and related offline APIs.

Context-aware classification
----------------------------
Predicates are **context-aware**: they consider whether the protein is
a self-protein (from the host organism) or a foreign protein, and whether
the protein is intended for therapeutic use.

- **Self-proteins** (organism matches source organism) auto-PASS because
  they are tolerated by the host immune system.
- **Foreign, non-therapeutic proteins** receive an EXPECTED_IMMUNOGENIC
  informational note rather than a hard FAIL, since immunogenicity is
  expected and not necessarily a problem.
- **Therapeutic proteins** have stricter thresholds because immune
  responses can compromise efficacy or cause adverse events.
"""

from __future__ import annotations

from typing import Any, Optional

from .types import Verdict, TypeCheckResult
from .immunogenicity import compute_immunogenicity
from .immunogenicity import predict_all as predict_mhc_all
from .immunogenicity import predict_epitopes
from .immunogenicity import predict_immunogenicity
from .immunogenicity import scan_peptides
from .immunogenicity import (
    DEFAULT_MHC_I_ALLELES,
    DEFAULT_MHC_II_ALLELES,
    IC50_STRONG_BINDER_THRESHOLD,
    binding_score_to_ic50,
)

__all__ = [
    "evaluate_low_immunogenicity",
    "evaluate_no_strong_t_cell_epitope",
    "evaluate_no_dominant_b_cell_epitope",
    "evaluate_population_coverage_safe",
]

# ────────────────────────────────────────────────────────────
# Verdict threshold constants
# ────────────────────────────────────────────────────────────
_IMMO_PASS_THRESHOLD: float = 0.2
_IMMO_LIKELY_PASS_THRESHOLD: float = 0.3
_IMMO_UNCERTAIN_THRESHOLD: float = 0.5
_IMMO_LIKELY_FAIL_THRESHOLD: float = 0.7
_T_CELL_STRONG_SCORE: float = 0.7
_B_CELL_DOMINANT_SCORE: float = 0.7
_POP_COVERAGE_PASS_THRESHOLD: float = 0.2
_POP_COVERAGE_LIKELY_PASS_THRESHOLD: float = 0.3
_POP_COVERAGE_LIKELY_FAIL_THRESHOLD: float = 0.7

# ────────────────────────────────────────────────────────────
# Context-aware classification constants
# ────────────────────────────────────────────────────────────
#: Organisms that are considered "self" (host organisms used as targets).
_SELF_ORGANISMS: frozenset[str] = frozenset({
    "Homo_sapiens", "human", "Human",
    "Mus_musculus", "mouse", "Mouse",
    "Cricetulus_griseus", "CHO", "CHO-K1", "cho",
})

#: IC50 threshold (nM) for "strong" T-cell epitopes in the predicate.
#: Only epitopes with IC50 < 50 nM are classified as "strong" for
#: predicate evaluation (stricter than the general 500 nM moderate threshold).
_T_CELL_STRONG_IC50_THRESHOLD: float = 50.0

#: Immunogenicity score threshold above which a therapeutic protein FAILs.
_THERAPEUTIC_IMMO_FAIL_THRESHOLD: float = 2.0

#: B-cell epitope score threshold above which a therapeutic protein FAILs.
_THERAPEUTIC_B_CELL_FAIL_SCORE: float = 1.5

# ────────────────────────────────────────────────────────────
# Organism-specific epitope databases (placeholder structure)
# ────────────────────────────────────────────────────────────
#: Mapping of organism -> set of known self-epitope peptide sequences
#: that should not be flagged.  Populated from organism databases.
_ORGANISM_SELF_EPITOPES: dict[str, set[str]] = {}


def _is_self_protein(organism: str, source_organism: str | None = None) -> bool:
    """Return True if the protein is from the host (self) organism.

    A protein is considered "self" when the source organism matches the
    target/host organism, meaning the immune system should be tolerant
    to it.

    Parameters
    ----------
    organism : str
        The target/host organism (e.g. ``"Homo_sapiens"``).
    source_organism : str or None
        The organism the protein originates from.  If *None*, defaults
        to *organism* (assumes self unless told otherwise).

    Returns
    -------
    bool
    """
    if source_organism is None:
        # If no source organism is specified, assume it's from the host
        return True
    # Normalise for comparison
    src = source_organism.strip()
    tgt = organism.strip()
    if src == tgt:
        return True
    # Also check against the known self-organism aliases
    return src in _SELF_ORGANISMS and tgt in _SELF_ORGANISMS


def evaluate_low_immunogenicity(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> TypeCheckResult:
    """Evaluate whether a protein has low overall immunogenicity.

    Context-aware classification:
      - **Self-protein** (source_organism matches host): auto-PASS.
      - **Foreign, non-therapeutic**: return UNCERTAIN with an
        ``EXPECTED_IMMUNOGENIC`` informational note instead of FAIL,
        because foreign proteins are *expected* to be immunogenic.
      - **Therapeutic protein**: only FAIL if immunogenicity score > 2.0
        AND ``therapeutic=True``, since immune responses can compromise
        drug efficacy.
      - **Non-therapeutic foreign proteins**: downgrade FAIL/LIKELY_FAIL
        to UNCERTAIN/LIKELY_PASS with an informational note.

    Score bands (when context does not override):
      PASS: score < 0.2
      LIKELY_PASS: < 0.3
      UNCERTAIN: < 0.5
      LIKELY_FAIL: < 0.7
      FAIL: >= 0.7

    Args:
        protein: Amino acid sequence (primary argument).
        sequence: Nucleotide sequence (unused; backward compat).
        organism: Target host organism (e.g. ``"Homo_sapiens"``).
        **kwargs: Optional keyword arguments:
            threshold: Immunogenicity score threshold (default 0.3).
            self_protein: If True, the protein is from the host organism
                (auto-PASS).  Defaults to None (auto-detected from
                source_organism).
            source_organism: Organism the protein originates from.
                If None, defaults to *organism* (assumed self).
            therapeutic: If True, the protein is intended for therapeutic
                use (stricter thresholds).  Default False.

    Returns:
        TypeCheckResult with immunogenicity verdict.
    """
    threshold = kwargs.get('threshold', 0.3)
    self_protein = kwargs.get('self_protein', None)
    source_organism = kwargs.get('source_organism', None)
    therapeutic = kwargs.get('therapeutic', False)

    if not protein:
        return TypeCheckResult(
            predicate="LowImmunogenicity",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    # Determine self-protein status
    is_self = self_protein if self_protein is not None else _is_self_protein(organism, source_organism)

    result = compute_immunogenicity(protein)
    score = result.overall_score

    derivation = [
        {"step": "compute_immunogenicity", "score": score,
         "t_cell_score": result.t_cell_score,
         "b_cell_score": result.b_cell_score,
         "immunogenicity_class": result.immunogenicity_class,
         "self_protein": is_self,
         "therapeutic": therapeutic,
         "source_organism": source_organism},
    ]

    # ── Self-protein auto-PASS ───────────────────────────────
    if is_self:
        return TypeCheckResult(
            predicate="LowImmunogenicity",
            verdict=Verdict.PASS,
            derivation=derivation,
            violation=None,
            knowledge_gap="Self-protein assumed tolerated; actual tolerance depends on expression context",
        )

    # ── Score-based classification ───────────────────────────
    if score < _IMMO_PASS_THRESHOLD:
        verdict = Verdict.PASS
        violation = None
    elif score < _IMMO_LIKELY_PASS_THRESHOLD:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif score < _IMMO_UNCERTAIN_THRESHOLD:
        verdict = Verdict.UNCERTAIN
        violation = f"Immunogenicity score {score:.3f} is moderate"
    elif score < _IMMO_LIKELY_FAIL_THRESHOLD:
        verdict = Verdict.LIKELY_FAIL
        violation = f"Immunogenicity score {score:.3f} is elevated"
    else:
        verdict = Verdict.FAIL
        violation = f"Immunogenicity score {score:.3f} is high (>{_IMMO_LIKELY_FAIL_THRESHOLD})"

    # ── Context-aware overrides for non-self proteins ─────────
    if therapeutic:
        # Therapeutic: only FAIL if score > 2.0 (very high immunogenicity)
        # The score-based bands above are retained for moderate levels,
        # but we upgrade the hard FAIL to only trigger at 2.0+
        if score > _THERAPEUTIC_IMMO_FAIL_THRESHOLD:
            verdict = Verdict.FAIL
            violation = (
                f"Immunogenicity score {score:.3f} is very high "
                f"(>{_THERAPEUTIC_IMMO_FAIL_THRESHOLD}) for therapeutic protein"
            )
        # If score-based verdict was already FAIL but < 2.0, downgrade to LIKELY_FAIL
        elif verdict == Verdict.FAIL:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"Immunogenicity score {score:.3f} is elevated for therapeutic protein "
                f"(threshold for FAIL: {_THERAPEUTIC_IMMO_FAIL_THRESHOLD})"
            )
    else:
        # Non-therapeutic foreign protein: immunogenicity is expected,
        # so downgrade hard FAIL/LIKELY_FAIL to informational
        if verdict == Verdict.FAIL:
            verdict = Verdict.UNCERTAIN
            violation = (
                f"EXPECTED_IMMUNOGENIC: Immunogenicity score {score:.3f} is high, "
                f"which is expected for a foreign (non-self) protein. "
                f"Not a concern unless used therapeutically."
            )
            derivation[-1]["classification"] = "EXPECTED_IMMUNOGENIC"
        elif verdict == Verdict.LIKELY_FAIL:
            verdict = Verdict.LIKELY_PASS
            violation = (
                f"EXPECTED_IMMUNOGENIC: Immunogenicity score {score:.3f} is elevated, "
                f"which is expected for a foreign (non-self) protein."
            )
            derivation[-1]["classification"] = "EXPECTED_IMMUNOGENIC"

    return TypeCheckResult(
        predicate="LowImmunogenicity",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )


def evaluate_no_strong_t_cell_epitope(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> TypeCheckResult:
    """Evaluate whether the protein lacks strong T-cell epitopes.

    Uses the offline scan_peptides API to scan all default MHC-I and
    MHC-II alleles without requiring external tools.

    Context-aware classification:
      - **Self-protein** (source_organism matches host): auto-PASS.
      - **Tiered classification for foreign proteins**:
          PASS: 0 strong epitopes (IC50 < 50 nM)
          WARN: 1-2 strong epitopes (LIKELY_PASS / UNCERTAIN)
          FAIL: >2 strong epitopes **only for therapeutic proteins**
      - **IC50 threshold**: Only epitopes with IC50 < 50 nM are
        classified as "strong" (stricter than the general 500 nM
        moderate threshold used in the immunogenicity module).

    Args:
        protein: Amino acid sequence (primary argument).
        sequence: Nucleotide sequence (unused; backward compat).
        organism: Target host organism (e.g. ``"Homo_sapiens"``).
        **kwargs: Optional keyword arguments:
            max_strong: Maximum allowed strong T-cell epitopes (default 0).
            source_organism: Organism the protein originates from.
                If None, defaults to *organism* (assumed self).
            therapeutic: If True, the protein is intended for therapeutic
                use (stricter FAIL threshold).  Default False.

    Returns:
        TypeCheckResult with T-cell epitope verdict.
    """
    max_strong = kwargs.get('max_strong', 0)
    source_organism = kwargs.get('source_organism', None)
    therapeutic = kwargs.get('therapeutic', False)

    if not protein:
        return TypeCheckResult(
            predicate="NoStrongTCellEpitope",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    # Determine self-protein status
    is_self = _is_self_protein(organism, source_organism)

    # Count strong epitopes using the stricter IC50 < 50 nM threshold
    # (not the broader binding_class which includes moderate_binders at < 500 nM)
    strong_count = 0
    total_epitopes = 0
    strong_epitope_details: list[dict[str, Any]] = []

    for allele in DEFAULT_MHC_I_ALLELES:
        results = scan_peptides(protein, allele, peptide_length=9)
        for r in results:
            total_epitopes += 1
            # Only flag as "strong" if IC50 < 50 nM (the predicate-level threshold)
            if r.ic50_nm < _T_CELL_STRONG_IC50_THRESHOLD:
                strong_count += 1
                strong_epitope_details.append({
                    "allele": allele, "peptide": r.peptide,
                    "ic50_nm": r.ic50_nm, "position": r.position,
                })

    for allele in DEFAULT_MHC_II_ALLELES:
        results = scan_peptides(protein, allele, peptide_length=15)
        for r in results:
            total_epitopes += 1
            if r.ic50_nm < _T_CELL_STRONG_IC50_THRESHOLD:
                strong_count += 1
                strong_epitope_details.append({
                    "allele": allele, "peptide": r.peptide,
                    "ic50_nm": r.ic50_nm, "position": r.position,
                })

    derivation = [
        {"step": "scan_peptides_offline", "total": total_epitopes,
         "strong": strong_count,
         "strong_ic50_threshold_nM": _T_CELL_STRONG_IC50_THRESHOLD,
         "self_protein": is_self,
         "therapeutic": therapeutic,
         "source_organism": source_organism},
    ]

    # ── Self-protein auto-PASS ───────────────────────────────
    if is_self:
        return TypeCheckResult(
            predicate="NoStrongTCellEpitope",
            verdict=Verdict.PASS,
            derivation=derivation,
            violation=None,
            knowledge_gap="Self-protein assumed tolerated; T-cell epitopes may be subject to central tolerance",
        )

    # ── Tiered classification ────────────────────────────────
    if strong_count <= max_strong:
        verdict = Verdict.PASS
        violation = None
    elif strong_count <= 1:
        # WARN tier: 1 strong epitope
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif strong_count <= 2:
        # WARN tier: 2 strong epitopes
        verdict = Verdict.UNCERTAIN
        violation = f"{strong_count} strong T-cell epitope(s) found (IC50 < {_T_CELL_STRONG_IC50_THRESHOLD} nM)"
    else:
        # >2 strong epitopes: FAIL for therapeutic, WARN for non-therapeutic
        if therapeutic:
            verdict = Verdict.FAIL
            violation = (
                f"{strong_count} strong T-cell epitopes found (IC50 < {_T_CELL_STRONG_IC50_THRESHOLD} nM) "
                f"— high immunogenicity risk for therapeutic protein"
            )
        else:
            verdict = Verdict.UNCERTAIN
            violation = (
                f"EXPECTED_IMMUNOGENIC: {strong_count} strong T-cell epitopes found "
                f"(IC50 < {_T_CELL_STRONG_IC50_THRESHOLD} nM). Expected for foreign protein; "
                f"would FAIL if used therapeutically."
            )
            derivation[-1]["classification"] = "EXPECTED_IMMUNOGENIC"

    return TypeCheckResult(
        predicate="NoStrongTCellEpitope",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )


def evaluate_no_dominant_b_cell_epitope(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> TypeCheckResult:
    """Evaluate whether the protein lacks dominant B-cell epitopes.

    Context-aware classification:
      - **Self-protein** (source_organism matches host): auto-PASS.
      - **Therapeutic proteins**: FAIL if B-cell epitope score > 1.5
        AND ``therapeutic=True``.
      - **Non-therapeutic proteins**: return PASS with an informational
        note for any B-cell epitopes detected, since B-cell epitopes
        in foreign proteins are expected.
      - Organism-specific epitope databases can be consulted via
        ``_ORGANISM_SELF_EPITOPES`` to skip known self-epitopes.

    Args:
        protein: Amino acid sequence (primary argument).
        sequence: Nucleotide sequence (unused; backward compat).
        organism: Target host organism (e.g. ``"Homo_sapiens"``).
        **kwargs: Optional keyword arguments:
            score_threshold: Score threshold for "dominant" epitopes
                (default 0.7).
            source_organism: Organism the protein originates from.
                If None, defaults to *organism* (assumed self).
            therapeutic: If True, the protein is intended for therapeutic
                use (stricter FAIL threshold).  Default False.

    Returns:
        TypeCheckResult with B-cell epitope verdict.
    """
    score_threshold = kwargs.get('score_threshold', _B_CELL_DOMINANT_SCORE)
    source_organism = kwargs.get('source_organism', None)
    therapeutic = kwargs.get('therapeutic', False)

    if not protein:
        return TypeCheckResult(
            predicate="NoDominantBCellEpitope",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    # Determine self-protein status
    is_self = _is_self_protein(organism, source_organism)

    result = predict_epitopes(protein)
    b_epitopes = result.linear_epitopes

    # Filter out organism-specific self-epitopes if database available
    self_epitopes = _ORGANISM_SELF_EPITOPES.get(organism, set())
    if self_epitopes:
        b_epitopes = [
            e for e in b_epitopes
            if e.peptide.upper() not in self_epitopes
        ]

    dominant_count = sum(1 for e in b_epitopes if e.score >= score_threshold)
    # Compute aggregate B-cell epitope score (sum of dominant epitope scores)
    dominant_score = sum(e.score for e in b_epitopes if e.score >= score_threshold)

    derivation = [
        {"step": "predict_b_cell_epitopes", "total": len(b_epitopes),
         "dominant": dominant_count, "dominant_score": dominant_score,
         "score_threshold": score_threshold,
         "self_protein": is_self,
         "therapeutic": therapeutic,
         "source_organism": source_organism},
    ]

    # ── Self-protein auto-PASS ───────────────────────────────
    if is_self:
        return TypeCheckResult(
            predicate="NoDominantBCellEpitope",
            verdict=Verdict.PASS,
            derivation=derivation,
            violation=None,
            knowledge_gap="Self-protein assumed tolerated; B-cell epitopes subject to peripheral tolerance",
        )

    # ── Non-self protein classification ──────────────────────
    if therapeutic:
        # Therapeutic: FAIL only if dominant score > 1.5
        if dominant_score > _THERAPEUTIC_B_CELL_FAIL_SCORE:
            verdict = Verdict.FAIL
            violation = (
                f"B-cell epitope score {dominant_score:.3f} exceeds therapeutic "
                f"threshold ({_THERAPEUTIC_B_CELL_FAIL_SCORE}) — "
                f"{dominant_count} dominant epitope(s) found"
            )
        elif dominant_count == 0:
            verdict = Verdict.PASS
            violation = None
        elif dominant_count == 1:
            verdict = Verdict.LIKELY_PASS
            violation = None
        elif dominant_count == 2:
            verdict = Verdict.UNCERTAIN
            violation = f"{dominant_count} dominant B-cell epitope(s) found (therapeutic context)"
        else:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"{dominant_count} dominant B-cell epitopes found — "
                f"elevated risk for therapeutic protein"
            )
    else:
        # Non-therapeutic: PASS with informational note
        if dominant_count == 0:
            verdict = Verdict.PASS
            violation = None
        else:
            verdict = Verdict.PASS
            violation = (
                f"INFORMATIONAL: {dominant_count} dominant B-cell epitope(s) found "
                f"(score threshold: {score_threshold}). Expected for foreign protein; "
                f"not a concern for non-therapeutic use."
            )
            derivation[-1]["classification"] = "EXPECTED_IMMUNOGENIC"

    return TypeCheckResult(
        predicate="NoDominantBCellEpitope",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )


def evaluate_population_coverage_safe(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> TypeCheckResult:
    """Evaluate whether the protein MHC binding profile is safe for
    population coverage.

    Uses the offline predict_all API (PSSM-based only).

    PASS: Binding rate < 0.2
    LIKELY_PASS: < 0.3
    UNCERTAIN: < coverage_threshold
    LIKELY_FAIL: < 0.7
    FAIL: >= 0.7

    Args:
        protein: Amino acid sequence (primary argument).
        sequence: Nucleotide sequence (unused; backward compat).
        organism: Target organism (unused by this predicate).
        **kwargs: Optional keyword arguments:
            coverage_threshold: Maximum allowed binding rate (default 0.5).

    Returns:
        TypeCheckResult with population coverage safety verdict.
    """
    coverage_threshold = kwargs.get('coverage_threshold', 0.5)
    if not protein:
        return TypeCheckResult(
            predicate="PopulationCoverageSafe",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    result = predict_mhc_all(protein)
    binding_rate = result.binding_rate
    num_binders = len(result.binders)
    num_strong = result.strong_binders

    if binding_rate < _POP_COVERAGE_PASS_THRESHOLD:
        verdict = Verdict.PASS
        violation = None
    elif binding_rate < _POP_COVERAGE_LIKELY_PASS_THRESHOLD:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif binding_rate < coverage_threshold:
        verdict = Verdict.UNCERTAIN
        violation = f"MHC binding rate {binding_rate:.3f} is moderate"
    elif binding_rate < _POP_COVERAGE_LIKELY_FAIL_THRESHOLD:
        verdict = Verdict.LIKELY_FAIL
        violation = f"MHC binding rate {binding_rate:.3f} is elevated"
    else:
        verdict = Verdict.FAIL
        violation = f"MHC binding rate {binding_rate:.3f} is high (broad population coverage of binders)"

    derivation = [
        {"step": "compute_population_coverage", "binding_rate": binding_rate,
         "num_binders": num_binders, "num_strong_binders": num_strong,
         "total_predictions": len(result.predictions)},
    ]

    return TypeCheckResult(
        predicate="PopulationCoverageSafe",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )
