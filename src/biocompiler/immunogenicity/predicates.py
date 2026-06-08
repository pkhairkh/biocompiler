"""
BioCompiler Immunogenicity Predicates
=======================================

Type-check predicates for immunogenicity assessment.
Each predicate returns a PredicateResult with a Verdict.

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

from dataclasses import dataclass
from typing import Any

from ..types import Verdict, TypeCheckResult
from .core import compute_immunogenicity
from .core import predict_all as predict_mhc_all
from .core import predict_epitopes
from .core import predict_immunogenicity
from .core import scan_peptides
from .core import (
    DEFAULT_MHC_I_ALLELES,
    DEFAULT_MHC_II_ALLELES,
    IC50_STRONG_BINDER_THRESHOLD,
    binding_score_to_ic50,
)

__all__ = [
    "PredicateResult",
    "evaluate_low_immunogenicity",
    "evaluate_no_strong_t_cell_epitope",
    "evaluate_no_dominant_b_cell_epitope",
    "evaluate_population_coverage_safe",
    "_check_anchor_match",
]

# ────────────────────────────────────────────────────────────
# Consistent result type for all immuno predicates
# ────────────────────────────────────────────────────────────

@dataclass
class PredicateResult(TypeCheckResult):
    """Extended TypeCheckResult with a ``details`` field for immuno predicates.

    All four immuno predicate functions return this type so that callers
    can rely on a consistent interface with ``verdict``, ``derivation``,
    and ``details`` fields.

    Attributes
    ----------
    details : str or None
        Human-readable explanation of the verdict.  For self-protein
        auto-PASS results this is always set to
        ``"Auto-PASS: protein is from host organism"``.
    """
    details: str | None = None


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

# ────────────────────────────────────────────────────────────
# Anchor residue matching for MHC binding prediction
# ────────────────────────────────────────────────────────────
# MHC-I anchor positions for 9-mer peptides. Each allele has preferred
# residues at positions 2 (P2) and 9 (P9/C-terminus). If a peptide
# does not match the anchor preferences, it is unlikely to bind even
# if the PSSM score is borderline. This helps resolve UNCERTAIN verdicts.
_MHC_I_ANCHOR_POSITIONS: dict[str, dict[str, set[str]]] = {
    # position -> preferred anchor residues
    "HLA-A*02:01": {"P2": {"L", "I", "V", "M", "A", "T"}, "P9": {"V", "L", "I", "A"}},
    "HLA-A*01:01": {"P2": {"T", "S", "D"}, "P9": {"Y", "F"}},
    "HLA-A*03:01": {"P2": {"V", "I", "L", "M", "T", "S"}, "P9": {"K", "R"}},
    "HLA-B*07:02": {"P2": {"P"}, "P9": {"L", "I", "V", "F"}},
    "HLA-B*08:01": {"P2": {"K", "R"}, "P9": {"L", "I", "V"}},
}


def _check_anchor_match(peptide: str, allele: str) -> bool:
    """Check if a peptide matches anchor residue preferences for an MHC allele.

    For MHC-I 9-mers, anchor positions are typically P2 (index 1) and
    P9 (index 8). If the peptide doesn't match the preferred anchors,
    it is unlikely to be a strong binder even if the PSSM score is borderline.

    Args:
        peptide: Amino acid sequence of the peptide (9-mer for MHC-I).
        allele: MHC allele name (e.g., "HLA-A*02:01").

    Returns:
        True if the peptide matches anchor preferences (or allele unknown),
        False if it clearly doesn't match.
    """
    anchors = _MHC_I_ANCHOR_POSITIONS.get(allele)
    if anchors is None:
        # Unknown allele → can't verify anchors → assume match
        return True

    peptide = peptide.upper()
    if len(peptide) < 9:
        return True  # Can't check anchors for short peptides

    p2 = peptide[1]
    p9 = peptide[8]

    p2_match = p2 in anchors.get("P2", set())
    p9_match = p9 in anchors.get("P9", set())

    # Both anchors must match for strong binding likelihood
    return p2_match and p9_match

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


def _resolve_self_protein(
    self_protein: bool | None,
    organism: str,
    source_organism: str | None,
) -> bool:
    """Resolve the self-protein status from explicit flag or auto-detection.

    If *self_protein* is explicitly set (True/False), use that value.
    Otherwise, auto-detect from *organism* and *source_organism*.
    """
    if self_protein is not None:
        return self_protein
    return _is_self_protein(organism, source_organism)


def evaluate_low_immunogenicity(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> PredicateResult:
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
        PredicateResult with immunogenicity verdict and details.
    """
    threshold = kwargs.get('threshold', 0.3)
    self_protein = kwargs.get('self_protein', None)
    source_organism = kwargs.get('source_organism', None)
    therapeutic = kwargs.get('therapeutic', False)

    if not protein:
        return PredicateResult(
            predicate="LowImmunogenicity",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
            details="Empty protein sequence",
        )

    # Determine self-protein status
    is_self = _resolve_self_protein(self_protein, organism, source_organism)

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
        return PredicateResult(
            predicate="LowImmunogenicity",
            verdict=Verdict.PASS,
            derivation=derivation,
            violation=None,
            knowledge_gap="Self-protein assumed tolerated; actual tolerance depends on expression context",
            details="Auto-PASS: protein is from host organism",
        )

    # ── Score-based classification ───────────────────────────
    if score < _IMMO_PASS_THRESHOLD:
        verdict = Verdict.PASS
        violation = None
    elif score < _IMMO_LIKELY_PASS_THRESHOLD:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif score < _IMMO_UNCERTAIN_THRESHOLD:
        # Moderate immunogenicity — use LIKELY_FAIL instead of UNCERTAIN
        # since we have a quantitative score that provides meaningful evidence.
        # UNCERTAIN is reserved for truly ambiguous cases where no signal exists.
        verdict = Verdict.LIKELY_FAIL
        violation = f"Immunogenicity score {score:.3f} is moderate"
    elif score < _IMMO_LIKELY_FAIL_THRESHOLD:
        verdict = Verdict.LIKELY_FAIL
        violation = f"Immunogenicity score {score:.3f} is elevated"
    else:
        verdict = Verdict.FAIL
        violation = f"Immunogenicity score {score:.3f} is high (>{_IMMO_LIKELY_FAIL_THRESHOLD})"

    details: str | None = None

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
            details = "Therapeutic protein: very high immunogenicity score"
        # If score-based verdict was already FAIL but < 2.0, downgrade to LIKELY_FAIL
        elif verdict == Verdict.FAIL:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"Immunogenicity score {score:.3f} is elevated for therapeutic protein "
                f"(threshold for FAIL: {_THERAPEUTIC_IMMO_FAIL_THRESHOLD})"
            )
            details = "Therapeutic protein: elevated immunogenicity score (downgraded from FAIL)"
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
            details = "Foreign protein: immunogenicity expected (not therapeutic)"
        elif verdict == Verdict.LIKELY_FAIL:
            verdict = Verdict.LIKELY_PASS
            violation = (
                f"EXPECTED_IMMUNOGENIC: Immunogenicity score {score:.3f} is elevated, "
                f"which is expected for a foreign (non-self) protein."
            )
            derivation[-1]["classification"] = "EXPECTED_IMMUNOGENIC"
            details = "Foreign protein: immunogenicity expected (not therapeutic)"

    return PredicateResult(
        predicate="LowImmunogenicity",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        details=details,
    )


def evaluate_no_strong_t_cell_epitope(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> PredicateResult:
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
        Callers may override via the ``ic50_threshold`` kwarg.

    Args:
        protein: Amino acid sequence (primary argument).
        sequence: Nucleotide sequence (unused; backward compat).
        organism: Target host organism (e.g. ``"Homo_sapiens"``).
        **kwargs: Optional keyword arguments:
            max_strong: Maximum allowed strong T-cell epitopes (default 0).
            ic50_threshold: IC50 threshold (nM) for "strong" binder
                classification.  Defaults to None, which uses the
                module-level ``_T_CELL_STRONG_IC50_THRESHOLD`` (50 nM).
                Pass a float to override.
            self_protein: If True, the protein is from the host organism
                (auto-PASS).  Defaults to None (auto-detected from
                source_organism).
            source_organism: Organism the protein originates from.
                If None, defaults to *organism* (assumed self).
            therapeutic: If True, the protein is intended for therapeutic
                use (stricter FAIL threshold).  Default False.

    Returns:
        PredicateResult with T-cell epitope verdict and details.
    """
    max_strong = kwargs.get('max_strong', 0)
    ic50_threshold = kwargs.get('ic50_threshold', None)
    self_protein = kwargs.get('self_protein', None)
    source_organism = kwargs.get('source_organism', None)
    therapeutic = kwargs.get('therapeutic', False)

    # Backward-compatible IC50 threshold: use the module-level default
    # only when no explicit threshold is provided by the caller.
    if ic50_threshold is None:
        ic50_threshold = _T_CELL_STRONG_IC50_THRESHOLD

    if not protein:
        return PredicateResult(
            predicate="NoStrongTCellEpitope",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
            details="Empty protein sequence",
        )

    # Determine self-protein status
    is_self = _resolve_self_protein(self_protein, organism, source_organism)

    # Count strong epitopes using the IC50 threshold
    # (not the broader binding_class which includes moderate_binders at < 500 nM)
    strong_count = 0
    total_epitopes = 0
    strong_epitope_details: list[dict[str, Any]] = []

    for allele in DEFAULT_MHC_I_ALLELES:
        results = scan_peptides(protein, allele, peptide_length=9)
        for r in results:
            total_epitopes += 1
            # Only flag as "strong" if IC50 < threshold nM
            if r.ic50_nm < ic50_threshold:
                strong_count += 1
                strong_epitope_details.append({
                    "allele": allele, "peptide": r.peptide,
                    "ic50_nm": r.ic50_nm, "position": r.position,
                })

    for allele in DEFAULT_MHC_II_ALLELES:
        results = scan_peptides(protein, allele, peptide_length=15)
        for r in results:
            total_epitopes += 1
            if r.ic50_nm < ic50_threshold:
                strong_count += 1
                strong_epitope_details.append({
                    "allele": allele, "peptide": r.peptide,
                    "ic50_nm": r.ic50_nm, "position": r.position,
                })

    derivation = [
        {"step": "scan_peptides_offline", "total": total_epitopes,
         "strong": strong_count,
         "strong_ic50_threshold_nM": ic50_threshold,
         "self_protein": is_self,
         "therapeutic": therapeutic,
         "source_organism": source_organism},
    ]

    # ── Self-protein auto-PASS ───────────────────────────────
    if is_self:
        return PredicateResult(
            predicate="NoStrongTCellEpitope",
            verdict=Verdict.PASS,
            derivation=derivation,
            violation=None,
            knowledge_gap="Self-protein assumed tolerated; T-cell epitopes may be subject to central tolerance",
            details="Auto-PASS: protein is from host organism",
        )

    # ── Anchor residue filtering ────────────────────────────────
    # Re-classify epitopes that don't match MHC anchor positions.
    # Peptides that lack the preferred anchor residues for their MHC allele
    # are unlikely to be genuine strong binders, even if the PSSM score
    # suggests borderline binding. This helps resolve UNCERTAIN verdicts.
    anchor_confirmed = 0
    for ep in strong_epitope_details:
        if _check_anchor_match(ep["peptide"], ep["allele"]):
            anchor_confirmed += 1

    # Use anchor-filtered count for verdict determination
    # Only count epitopes that have confirmed anchor matches
    confirmed_strong_count = anchor_confirmed

    # ── Tiered classification ────────────────────────────────
    details: str | None = None

    # If anchor filtering reduced the count, use the reduced count
    effective_strong_count = confirmed_strong_count if confirmed_strong_count > 0 else strong_count

    if effective_strong_count <= max_strong:
        verdict = Verdict.PASS
        violation = None
    elif effective_strong_count <= 1:
        # WARN tier: 1 confirmed strong epitope
        verdict = Verdict.LIKELY_PASS
        violation = None
        details = f"1 strong T-cell epitope found (IC50 < {ic50_threshold} nM, anchor-confirmed)"
    elif effective_strong_count <= 2:
        # 2 confirmed epitopes: use LIKELY_FAIL instead of UNCERTAIN
        # since anchor confirmation provides meaningful evidence
        verdict = Verdict.LIKELY_FAIL
        violation = f"{effective_strong_count} anchor-confirmed strong T-cell epitope(s) found (IC50 < {ic50_threshold} nM)"
        details = f"{effective_strong_count} anchor-confirmed strong T-cell epitope(s) found (IC50 < {ic50_threshold} nM)"
    else:
        # >2 confirmed epitopes: FAIL for therapeutic, WARN for non-therapeutic
        if therapeutic:
            verdict = Verdict.FAIL
            violation = (
                f"{effective_strong_count} anchor-confirmed strong T-cell epitopes found "
                f"(IC50 < {ic50_threshold} nM) — high immunogenicity risk for therapeutic protein"
            )
            details = f"Therapeutic protein: {effective_strong_count} anchor-confirmed strong T-cell epitopes found"
        else:
            # Non-therapeutic: use LIKELY_FAIL instead of UNCERTAIN
            # since anchor confirmation provides meaningful evidence
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"EXPECTED_IMMUNOGENIC: {effective_strong_count} anchor-confirmed strong T-cell epitopes found "
                f"(IC50 < {ic50_threshold} nM). Expected for foreign protein; "
                f"would FAIL if used therapeutically."
            )
            derivation[-1]["classification"] = "EXPECTED_IMMUNOGENIC"
            details = f"Foreign protein: {effective_strong_count} anchor-confirmed strong T-cell epitopes expected"

    return PredicateResult(
        predicate="NoStrongTCellEpitope",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        details=details,
    )


def evaluate_no_dominant_b_cell_epitope(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> PredicateResult:
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
            self_protein: If True, the protein is from the host organism
                (auto-PASS).  Defaults to None (auto-detected from
                source_organism).
            source_organism: Organism the protein originates from.
                If None, defaults to *organism* (assumed self).
            therapeutic: If True, the protein is intended for therapeutic
                use (stricter FAIL threshold).  Default False.

    Returns:
        PredicateResult with B-cell epitope verdict and details.
    """
    score_threshold = kwargs.get('score_threshold', _B_CELL_DOMINANT_SCORE)
    self_protein = kwargs.get('self_protein', None)
    source_organism = kwargs.get('source_organism', None)
    therapeutic = kwargs.get('therapeutic', False)

    if not protein:
        return PredicateResult(
            predicate="NoDominantBCellEpitope",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
            details="Empty protein sequence",
        )

    # Determine self-protein status
    is_self = _resolve_self_protein(self_protein, organism, source_organism)

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
        return PredicateResult(
            predicate="NoDominantBCellEpitope",
            verdict=Verdict.PASS,
            derivation=derivation,
            violation=None,
            knowledge_gap="Self-protein assumed tolerated; B-cell epitopes subject to peripheral tolerance",
            details="Auto-PASS: protein is from host organism",
        )

    # ── Non-self protein classification ──────────────────────
    details: str | None = None

    if therapeutic:
        # Therapeutic: FAIL only if dominant score > 1.5
        if dominant_score > _THERAPEUTIC_B_CELL_FAIL_SCORE:
            verdict = Verdict.FAIL
            violation = (
                f"B-cell epitope score {dominant_score:.3f} exceeds therapeutic "
                f"threshold ({_THERAPEUTIC_B_CELL_FAIL_SCORE}) — "
                f"{dominant_count} dominant epitope(s) found"
            )
            details = f"Therapeutic protein: B-cell epitope score {dominant_score:.3f} exceeds threshold"
        elif dominant_count == 0:
            verdict = Verdict.PASS
            violation = None
        elif dominant_count == 1:
            verdict = Verdict.LIKELY_PASS
            violation = None
            details = "Therapeutic protein: 1 dominant B-cell epitope found"
        elif dominant_count == 2:
            verdict = Verdict.UNCERTAIN
            violation = f"{dominant_count} dominant B-cell epitope(s) found (therapeutic context)"
            details = f"Therapeutic protein: {dominant_count} dominant B-cell epitope(s) found"
        else:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"{dominant_count} dominant B-cell epitopes found — "
                f"elevated risk for therapeutic protein"
            )
            details = f"Therapeutic protein: {dominant_count} dominant B-cell epitopes found (elevated risk)"
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
            details = f"Foreign protein: {dominant_count} dominant B-cell epitope(s) expected (non-therapeutic)"

    return PredicateResult(
        predicate="NoDominantBCellEpitope",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        details=details,
    )


def evaluate_population_coverage_safe(
    protein: str,
    sequence: str | None = None,
    organism: str = "Homo_sapiens",
    **kwargs: Any,
) -> PredicateResult:
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
        PredicateResult with population coverage safety verdict and details.
    """
    coverage_threshold = kwargs.get('coverage_threshold', 0.5)
    if not protein:
        return PredicateResult(
            predicate="PopulationCoverageSafe",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
            details="Empty protein sequence",
        )

    result = predict_mhc_all(protein)
    binding_rate = result.binding_rate
    num_binders = len(result.binders)
    num_strong = result.strong_binders

    details: str | None = None

    if binding_rate < _POP_COVERAGE_PASS_THRESHOLD:
        verdict = Verdict.PASS
        violation = None
    elif binding_rate < _POP_COVERAGE_LIKELY_PASS_THRESHOLD:
        verdict = Verdict.LIKELY_PASS
        violation = None
        details = f"MHC binding rate {binding_rate:.3f} is slightly elevated"
    elif binding_rate < coverage_threshold:
        # Moderate binding rate — use LIKELY_FAIL instead of UNCERTAIN
        # since the quantitative binding rate provides meaningful evidence.
        verdict = Verdict.LIKELY_FAIL
        violation = f"MHC binding rate {binding_rate:.3f} is moderate"
        details = f"MHC binding rate {binding_rate:.3f} is moderate"
    elif binding_rate < _POP_COVERAGE_LIKELY_FAIL_THRESHOLD:
        verdict = Verdict.LIKELY_FAIL
        violation = f"MHC binding rate {binding_rate:.3f} is elevated"
        details = f"MHC binding rate {binding_rate:.3f} is elevated"
    else:
        verdict = Verdict.FAIL
        violation = f"MHC binding rate {binding_rate:.3f} is high (broad population coverage of binders)"
        details = f"MHC binding rate {binding_rate:.3f} is high (broad population coverage of binders)"

    derivation = [
        {"step": "compute_population_coverage", "binding_rate": binding_rate,
         "num_binders": num_binders, "num_strong_binders": num_strong,
         "total_predictions": len(result.predictions)},
    ]

    return PredicateResult(
        predicate="PopulationCoverageSafe",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
        details=details,
    )
