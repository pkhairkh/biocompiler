"""
BioCompiler Immunogenicity Predicates
=======================================

Type-check predicates for immunogenicity assessment.
Each predicate returns a TypeCheckResult with a Verdict.
"""

from __future__ import annotations

from typing import Any, Optional

from .types import Verdict, TypeCheckResult
from .immunogenicity import compute_immunogenicity
from .immunogenicity import predict_all as predict_mhc_all
from .immunogenicity import predict_epitopes


def evaluate_low_immunogenicity(sequence: str, protein: str, organism: str, **kwargs: Any) -> TypeCheckResult:
    """Evaluate whether a protein has low overall immunogenicity.

    PASS: Score < 0.2 (very low risk)
    LIKELY_PASS: < 0.3
    UNCERTAIN: < 0.5
    LIKELY_FAIL: < 0.7
    FAIL: >= 0.7 (high risk)

    Args:
        sequence: Nucleotide sequence (unused by this predicate).
        protein: Amino acid sequence.
        organism: Target organism (unused by this predicate).
        **kwargs: Optional keyword arguments:
            threshold: Immunogenicity score threshold (default 0.3).

    Returns:
        TypeCheckResult with immunogenicity verdict.
    """
    threshold = kwargs.get('threshold', 0.3)
    if not protein:
        return TypeCheckResult(
            predicate="LowImmunogenicity",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    result = compute_immunogenicity(protein)
    score = result.overall_score

    if score < 0.2:
        verdict = Verdict.PASS
        violation = None
    elif score < 0.3:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif score < 0.5:
        verdict = Verdict.UNCERTAIN
        violation = f"Immunogenicity score {score:.3f} is moderate"
    elif score < 0.7:
        verdict = Verdict.LIKELY_FAIL
        violation = f"Immunogenicity score {score:.3f} is elevated"
    else:
        verdict = Verdict.FAIL
        violation = f"Immunogenicity score {score:.3f} is high (>{0.7})"

    derivation = [
        {"step": "compute_immunogenicity", "score": score,
         "t_cell_score": result.t_cell_score,
         "b_cell_score": result.b_cell_score,
         "immunogenicity_class": result.immunogenicity_class},
    ]

    return TypeCheckResult(
        predicate="LowImmunogenicity",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )


def evaluate_no_strong_t_cell_epitope(sequence: str, protein: str, organism: str, **kwargs: Any) -> TypeCheckResult:
    """Evaluate whether the protein lacks strong T-cell epitopes.

    PASS: No strong T-cell epitopes
    LIKELY_PASS: 1 strong epitope
    UNCERTAIN: 2 strong epitopes
    LIKELY_FAIL: 3 strong epitopes
    FAIL: >3 strong T-cell epitopes

    Args:
        sequence: Nucleotide sequence (unused by this predicate).
        protein: Amino acid sequence.
        organism: Target organism (unused by this predicate).
        **kwargs: Optional keyword arguments:
            max_strong: Maximum allowed strong T-cell epitopes (default 0).

    Returns:
        TypeCheckResult with T-cell epitope verdict.
    """
    max_strong = kwargs.get('max_strong', 0)
    if not protein:
        return TypeCheckResult(
            predicate="NoStrongTCellEpitope",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    result = predict_epitopes(protein)
    t_epitopes = result.linear_epitopes
    strong_count = sum(1 for e in t_epitopes if e.score >= 0.7)

    if strong_count <= max_strong:
        verdict = Verdict.PASS
        violation = None
    elif strong_count == 1:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif strong_count == 2:
        verdict = Verdict.UNCERTAIN
        violation = f"{strong_count} strong T-cell epitope(s) found"
    elif strong_count == 3:
        verdict = Verdict.LIKELY_FAIL
        violation = f"{strong_count} strong T-cell epitopes found"
    else:
        verdict = Verdict.FAIL
        violation = f"{strong_count} strong T-cell epitopes found (high immunogenicity risk)"

    derivation = [
        {"step": "predict_t_cell_epitopes", "total": len(t_epitopes),
         "strong": strong_count},
    ]

    return TypeCheckResult(
        predicate="NoStrongTCellEpitope",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )


def evaluate_no_dominant_b_cell_epitope(sequence: str, protein: str, organism: str, **kwargs: Any) -> TypeCheckResult:
    """Evaluate whether the protein lacks dominant B-cell epitopes.

    PASS: No B-cell epitopes with score >= threshold
    LIKELY_PASS: 1 dominant epitope
    UNCERTAIN: 2 dominant epitopes
    FAIL: >2 dominant B-cell epitopes

    Args:
        sequence: Nucleotide sequence (unused by this predicate).
        protein: Amino acid sequence.
        organism: Target organism (unused by this predicate).
        **kwargs: Optional keyword arguments:
            score_threshold: Score threshold for "dominant" epitopes (default 0.7).

    Returns:
        TypeCheckResult with B-cell epitope verdict.
    """
    score_threshold = kwargs.get('score_threshold', 0.7)
    if not protein:
        return TypeCheckResult(
            predicate="NoDominantBCellEpitope",
            verdict=Verdict.UNCERTAIN,
            violation="Empty protein sequence",
        )

    result = predict_epitopes(protein)
    b_epitopes = result.linear_epitopes
    dominant_count = sum(1 for e in b_epitopes if e.score >= score_threshold)

    if dominant_count == 0:
        verdict = Verdict.PASS
        violation = None
    elif dominant_count == 1:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif dominant_count == 2:
        verdict = Verdict.UNCERTAIN
        violation = f"{dominant_count} dominant B-cell epitope(s) found"
    else:
        verdict = Verdict.FAIL
        violation = f"{dominant_count} dominant B-cell epitopes found (high immunogenicity risk)"

    derivation = [
        {"step": "predict_b_cell_epitopes", "total": len(b_epitopes),
         "dominant": dominant_count},
    ]

    return TypeCheckResult(
        predicate="NoDominantBCellEpitope",
        verdict=verdict,
        derivation=derivation,
        violation=violation,
    )


def evaluate_population_coverage_safe(
    sequence: str,
    protein: str,
    organism: str,
    **kwargs: Any,
) -> TypeCheckResult:
    """Evaluate whether the protein's MHC binding profile is safe for
    population coverage (i.e., doesn't bind too many common alleles).

    PASS: Binding rate < 0.2 (low population coverage of binders)
    LIKELY_PASS: < 0.3
    UNCERTAIN: < coverage_threshold
    LIKELY_FAIL: < 0.7
    FAIL: >= 0.7 (high population coverage of binders = high immunogenicity risk)

    Args:
        sequence: Nucleotide sequence (unused by this predicate).
        protein: Amino acid sequence.
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

    if binding_rate < 0.2:
        verdict = Verdict.PASS
        violation = None
    elif binding_rate < 0.3:
        verdict = Verdict.LIKELY_PASS
        violation = None
    elif binding_rate < coverage_threshold:
        verdict = Verdict.UNCERTAIN
        violation = f"MHC binding rate {binding_rate:.3f} is moderate"
    elif binding_rate < 0.7:
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
