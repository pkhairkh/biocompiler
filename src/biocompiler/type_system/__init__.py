"""
BioCompiler Type System v1.0.0
===============================
Defines the core types, codon tables, BLOSUM62 matrix, and 43 predicate classes
for certified gene optimization: 12 DNA-level + 4 structure + 4 stability +
4 solubility + 4 immunogenicity + 15 extended (miRNA, m6A, PolyA, RQC, Alu,
CrypticORF, BlastMatches, PrimerCompat, NucleosideModGuidance, SlidingGC, etc.).

This package decomposes the original monolith into submodules:
  - codon_tables: Core data tables and data structures
  - checks: Low-level check_* functions returning PredicateResult
  - predicates: High-level evaluate_* functions returning TypeCheckResult
  - logic: Removed (use biocompiler.shared.five_valued_logic and biocompiler.shared.types directly)
  - registry: PredicateRegistry for named dispatch
"""

# Re-export types that other modules import from type_system
# (backward compatibility — these were available in the old monolith)
from biocompiler.shared.types import Verdict, SLOTMode, TypeCheckResult

# Core data tables and types
from .codon_tables import (
    AA_TO_CODONS,
    BLOSUM62,
    CODON_TABLE,
    CertLevel,
    PREDICATE_NAMES,
    PROMOTER_CONSENSUS,
    PredicateResult,
    SpliceVerdict,
    START_CODONS,
    STOP_CODONS,
)

# Low-level predicate checks
from .checks import (
    check_co_translational_folding,
    check_codon_optimality,
    check_conservation_score,
    check_mrna_secondary_structure,
    check_mrna_stability,
    check_no_alu_repeat,
    check_no_avoidable_gt,
    check_no_blast_matches,
    check_no_cpg_island,
    check_no_cryptic_orf,
    check_no_cryptic_promoter,
    check_no_cryptic_splice,
    check_no_m6a_site,
    check_no_mirna_binding_site,
    check_no_polya_signal,
    check_no_rqc_trigger,
    check_nucleoside_modification_guidance,
    check_no_gt_dinucleotide,
    check_no_gt_dinucleotide_soft,
    check_no_restriction_site,
    check_no_stop_codons,
    check_no_unexpected_tm_domain,
    check_primer_compatibility,
    check_valid_coding_seq,
    find_cross_codon_cg,
    find_cross_codon_gt,
    find_cross_codon_restriction,
    # Translation helper (used by mutagenesis tests)
    _translate_dna_to_aa,
    # Organism-aware helpers (public for backward compat)
    _is_prokaryotic_organism,
    _compute_max_gt_count,
    # Additional underscore helpers re-exported for backward compat.
    # Consumed directly by tests/test_uncertain_reduction.py and other
    # call sites that import from `biocompiler.type_system` (the
    # checks.py shim already re-exports these via star-import; this
    # explicit re-export honours the documented public surface).
    _resolve_species_cai,
    _count_dinucs_fast,
    _compute_codon_ramp_score,
    _rna_revcomp_to_dna,
    _mirna_context_score,
)

# High-level evaluate API
from .predicates import (
    analyze_codon_at_position,
    evaluate_all_predicates,
    evaluate_co_translational_folding,
    evaluate_codon_adapted,
    evaluate_codon_optimality,
    evaluate_conservation_score,
    evaluate_gc_in_range,
    evaluate_in_frame,
    evaluate_mrna_secondary_structure,
    evaluate_mrna_stability,
    evaluate_no_alu_repeat,
    evaluate_no_blast_matches,
    evaluate_no_cpg_island,
    evaluate_no_cryptic_orf,
    evaluate_no_cryptic_promoter,
    evaluate_no_cryptic_splice,
    evaluate_no_m6a_site,
    evaluate_no_mirna_binding_site,
    evaluate_no_gt_dinucleotide,
    evaluate_no_instability_motif,
    evaluate_no_polya_signal,
    evaluate_no_restriction_site,
    evaluate_no_rqc_trigger,
    evaluate_no_stop_codons,
    evaluate_no_unexpected_tm_domain,
    evaluate_nucleoside_modification_guidance,
    evaluate_primer_compatibility,
    evaluate_splice_correct,
    evaluate_valid_coding_seq,
)

# Sliding-window GC (re-exported from sliding_gc module)
from biocompiler.sequence.sliding_gc import (
    evaluate_sliding_gc,
    check_sliding_gc,
    SlidingGCResult,
    WindowViolation,
)

# miRNA seed database (multi-organism, tissue-filtered)
from .mirna_seeds import get_mirna_seeds

# Registry
from .registry import PredicateRegistry, registry

# Eukaryote GT constant (public for backward compat)
from .codon_tables import _EUKARYOTE_GT_PER_BP

# Internal constants re-exported for backward compatibility
from .codon_tables import (
    _TM_EUKARYOTIC_MIN_STRETCH,
    _TM_PROKARYOTIC_MIN_STRETCH,
    _MRNA_DG_PROKARYOTE_FAIL,
    _MRNA_DG_EUKARYOTE_FAIL,
    _RESTRICTION_SITE_MIN_LENGTH,
    _TM_BORDERLINE_RATIO,
)


__all__ = [
    # Core data tables
    "CODON_TABLE", "AA_TO_CODONS", "BLOSUM62",
    # Enums
    "CertLevel", "SpliceVerdict",
    # Predicate names
    "PREDICATE_NAMES",
    # Promoter data
    "PROMOTER_CONSENSUS",
    # Data classes
    "PredicateResult",
    # Low-level predicate checks
    "check_no_stop_codons", "check_no_cryptic_splice", "check_no_cpg_island",
    "check_no_restriction_site", "check_no_gt_dinucleotide", "check_no_avoidable_gt",
    "check_no_gt_dinucleotide_soft",
    "check_valid_coding_seq", "check_conservation_score", "check_codon_optimality",
    "check_no_cryptic_promoter", "check_no_unexpected_tm_domain",
    "check_mrna_secondary_structure", "check_co_translational_folding",
    "check_mrna_stability", "evaluate_mrna_stability",
    "check_no_blast_matches", "check_primer_compatibility",
    "check_no_cryptic_orf", "check_no_rqc_trigger", "check_no_alu_repeat",
    "check_no_mirna_binding_site", "check_no_m6a_site", "check_no_polya_signal",
    "check_nucleoside_modification_guidance",
    # High-level evaluate API
    "evaluate_gc_in_range", "evaluate_no_cryptic_splice", "evaluate_splice_correct",
    "evaluate_codon_adapted", "evaluate_no_restriction_site", "evaluate_in_frame",
    "evaluate_no_instability_motif", "evaluate_no_unexpected_tm_domain",
    "evaluate_mrna_secondary_structure", "evaluate_no_cryptic_promoter",
    "evaluate_no_cpg_island", "analyze_codon_at_position",
    "evaluate_co_translational_folding", "evaluate_all_predicates",
    "evaluate_no_stop_codons", "evaluate_no_gt_dinucleotide",
    "evaluate_valid_coding_seq", "evaluate_conservation_score",
    "evaluate_codon_optimality",
    "evaluate_no_blast_matches", "evaluate_primer_compatibility",
    "evaluate_no_cryptic_orf", "evaluate_no_rqc_trigger", "evaluate_no_alu_repeat",
    "evaluate_no_mirna_binding_site", "evaluate_no_m6a_site", "evaluate_no_polya_signal",
    "evaluate_nucleoside_modification_guidance",
    # Sliding-window GC
    "evaluate_sliding_gc", "check_sliding_gc", "SlidingGCResult", "WindowViolation",
    # Cross-codon helpers
    "find_cross_codon_gt", "find_cross_codon_cg", "find_cross_codon_restriction",
    # Translation helper
    "_translate_dna_to_aa",
    # Organism-aware helpers
    "_is_prokaryotic_organism", "_compute_max_gt_count", "_EUKARYOTE_GT_PER_BP",
    # Additional underscore helpers (backward-compat re-exports)
    "_resolve_species_cai", "_count_dinucs_fast", "_compute_codon_ramp_score",
    "_rna_revcomp_to_dna", "_mirna_context_score",
    # Registry
    "PredicateRegistry", "registry",
    # miRNA seed database
    "get_mirna_seeds",
]
# TEST MARKER 1781908261
