"""
Strategy path functions for optimize_sequence().

After second-pass cleanup (Task SP2), this module retains only
``evaluate_extended_predicates`` — the only function used by the
default fast path (``pipeline_certification.py`` imports it).

The slow-path helpers (``run_csp_solver_path``, ``run_harmonize_path``,
``run_prokaryote_hybrid_path``) referenced deleted modules (``greedy``,
``codon_harmonization``, the hybrid stack) and have been removed in
lockstep with the file deletions. The corresponding call sites in
``pipeline_core.py`` are being scrubbed by SP1.

Decomposition: Extracted from pipeline_core.py (Task pipeline-decompose).
"""

import logging

from ..type_system import (
    PredicateResult,
    check_mrna_stability as _check_mrna_stability,
    check_no_cryptic_promoter as _check_no_cryptic_promoter,
    check_no_unexpected_tm_domain as _check_no_unexpected_tm_domain,
    check_no_cryptic_orf as _check_no_cryptic_orf,
    check_no_rqc_trigger as _check_no_rqc_trigger,
    check_no_alu_repeat as _check_no_alu_repeat,
    check_no_mirna_binding_site as _check_no_mirna_binding_site,
    check_no_m6a_site as _check_no_m6a_site,
    check_no_polya_signal as _check_no_polya_signal,
    check_mrna_secondary_structure as _check_mrna_secondary_structure,
)

logger = logging.getLogger(__name__)


def evaluate_extended_predicates(
    optimized_seq: str,
    organism: str,
    pred_results: list,
    tissue: str = "",
) -> list:
    """Append extended diagnostic predicates (11-20) to pred_results.

    These were previously only evaluated in the BioOptimizer
    _evaluate_all_predicates path, but the hybrid fast path
    skipped them. Now wired in so that the competitive
    landscape claims are genuine for all paths.

    Args:
        optimized_seq: The optimized DNA sequence.
        organism: Target organism name.
        pred_results: List to append predicate results to.
        tissue: Tissue context for tissue-specific checks (e.g. miRNA).
            Empty string means no tissue filter (genome-wide check).
    """
    # 11. NoInstabilityMotif
    try:
        pred_results.append(_check_mrna_stability(optimized_seq, organism=organism))
    except Exception:
        _attta_pos = [i for i in range(len(optimized_seq) - 4) if optimized_seq[i:i+5] == "ATTTA"]
        pred_results.append(PredicateResult(
            "NoInstabilityMotif", len(_attta_pos) == 0,
            details=f"Found {len(_attta_pos)} ATTTA motifs (fallback)"
        ))

    # 12. NoCrypticPromoter
    try:
        pred_results.append(_check_no_cryptic_promoter(optimized_seq, organism=organism, threshold=0.7))
    except Exception:
        pred_results.append(PredicateResult("NoCrypticPromoter", True, details="Unavailable"))

    # 13. NoUnexpectedTMDomain
    try:
        pred_results.append(_check_no_unexpected_tm_domain(optimized_seq, is_cytosolic=True, organism=organism))
    except Exception:
        pred_results.append(PredicateResult("NoUnexpectedTMDomain", True, details="Unavailable"))

    # 14. NoCrypticORF
    try:
        pred_results.append(_check_no_cryptic_orf(optimized_seq, min_orf_length=30, organism=organism))
    except Exception:
        pred_results.append(PredicateResult("NoCrypticORF", True, details="Unavailable"))

    # 15. NoRQCTrigger
    try:
        pred_results.append(_check_no_rqc_trigger(optimized_seq, organism=organism))
    except Exception:
        pred_results.append(PredicateResult("NoRQCTrigger", True, details="Unavailable"))

    # 16. NoAluRepeat
    try:
        pred_results.append(_check_no_alu_repeat(optimized_seq, organism=organism))
    except Exception:
        pred_results.append(PredicateResult("NoAluRepeat", True, details="Unavailable"))

    # 17. NoMiRNABindingSite (eukaryotes only)
    try:
        mirna_kwargs = {"organism": organism, "min_seed_match": 7}
        if tissue:
            mirna_kwargs["tissue"] = tissue
        pred_results.append(_check_no_mirna_binding_site(optimized_seq, **mirna_kwargs))
    except Exception:
        pred_results.append(PredicateResult("NoMiRNABindingSite", True, details="Unavailable"))

    # 18. NoM6ASite
    try:
        pred_results.append(_check_no_m6a_site(optimized_seq, organism=organism))
    except Exception:
        pred_results.append(PredicateResult("NoM6ASite", True, details="Unavailable"))

    # 19. NoPolyASignal
    try:
        pred_results.append(_check_no_polya_signal(optimized_seq, organism=organism))
    except Exception:
        pred_results.append(PredicateResult("NoPolyASignal", True, details="Unavailable"))

    # 20. mRNASecondaryStructure (ViennaRNA / Nussinov fallback)
    try:
        pred_results.append(_check_mrna_secondary_structure(optimized_seq, organism=organism, use_viennarna=True))
    except Exception:
        pred_results.append(PredicateResult("mRNASecondaryStructure", True, details="Unavailable"))

    return pred_results
