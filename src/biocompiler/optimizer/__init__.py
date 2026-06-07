"""
biocompiler.optimizer — Modular optimization subpackage.

This package decomposes the former monolithic optimization.py into
focused submodules while maintaining full backward compatibility.

Submodules:
    cai         — CAI computation and batch scoring
    constraints — Restriction site removal, GT/AG helpers, back-translation
    greedy      — Greedy optimizer loop and splice donor scoring
    pipeline    — BioOptimizer class, optimize_sequence(), batch_optimize()
    utils       — ConvergenceTracker, OptimizationResult, protein_to_aa_list()
"""

# Re-export all public symbols for backward compatibility.
# Users can do:  from biocompiler.optimizer import BioOptimizer
# Or:            from biocompiler.optimization import BioOptimizer  (deprecated shim)

from .utils import (
    ConvergenceTracker,
    OptimizationResult,
    FullConstructResult,
    protein_to_aa_list,
    MAX_RESTRICTION_SITE_ITERATIONS,
    MAX_IUPAC_SITE_ITERATIONS,
    MAX_ATTTA_MOTIF_ITERATIONS,
    MAX_T_RUN_ITERATIONS,
    MAX_GC_ADJUSTMENT_ITERATIONS,
    MAX_SPLICE_ELIMINATION_ITERATIONS,
    MAX_CPG_DISRUPTION_ITERATIONS,
    DEFAULT_MAX_ITERATIONS,
    CONVERGENCE_IMPROVEMENT_THRESHOLD,
    CONVERGENCE_PATIENCE,
    OSCILLATION_WINDOW,
)

from .cai import (
    HAS_NUMBA,
    USE_NUMBA,
    _HAS_NUMBA,
    _adaptiveness_to_array,
    _codon_to_index,
    _dna_to_codon_indices,
    _compute_cai_fast,
    _count_dinucs_fast,
    _BatchSwapScorer,
)

from .constraints import (
    _find_site_in_sequence,
    _get_overlapping_codons,
    _remove_site_multicodon,
    _find_gt_free_codons,
    _find_ag_free_codons,
    _organism_to_species_key,
    _species_key_to_organism,
    _back_translate_protein,
    _back_translate_protein_dp,
    _build_restriction_site_set,
    _contains_restriction_site,
    _count_dinucleotides,
    _count_gts,
    _is_unavoidable_gt,
    _has_gt,
    _codon_creates_boundary_gt,
    _DP_TOP_K,
    _GT_AG_PENALTY,
    _RESTRICTION_PENALTY,
)

from .greedy import (
    score_splice_donor_potential,
    _gt_aware_select_codon,
    _is_in_codon_gt,
    _eukaryote_cai_recovery,
    _eliminate_cpg_dinucleotides,
    _greedy_optimize,
    _expand_iupac_site,
    _check_predicates_via_type_system,
)

from .pipeline import (
    optimize_sequence,
    batch_optimize,
    BioOptimizer,
)

# Also re-export the constants that were part of the original __all__
from .greedy import (
    SPLICE_DONOR_POTENTIAL_THRESHOLD,
    EUKARYOTE_CAI_GT_COST_THRESHOLD,
    GT_BOUNDARY_CAI_TOLERANCE,
    GT_CAI_LOG_ADAPTIVENESS_COST,
)

__all__ = [
    "OptimizationResult",
    "FullConstructResult",
    "optimize_sequence",
    "batch_optimize",
    "protein_to_aa_list",
    "BioOptimizer",
    "ConvergenceTracker",
    "score_splice_donor_potential",
    "SPLICE_DONOR_POTENTIAL_THRESHOLD",
    "EUKARYOTE_CAI_GT_COST_THRESHOLD",
    "GT_BOUNDARY_CAI_TOLERANCE",
    "GT_CAI_LOG_ADAPTIVENESS_COST",
    # Named constants (from utils.py)
    "DEFAULT_MAX_ITERATIONS",
    "CONVERGENCE_IMPROVEMENT_THRESHOLD",
    "CONVERGENCE_PATIENCE",
    "OSCILLATION_WINDOW",
    "MAX_RESTRICTION_SITE_ITERATIONS",
    "MAX_IUPAC_SITE_ITERATIONS",
    "MAX_ATTTA_MOTIF_ITERATIONS",
    "MAX_T_RUN_ITERATIONS",
    "MAX_GC_ADJUSTMENT_ITERATIONS",
    "MAX_SPLICE_ELIMINATION_ITERATIONS",
    "MAX_CPG_DISRUPTION_ITERATIONS",
    # Also expose internal symbols that were part of the original module
    "HAS_NUMBA",
    "USE_NUMBA",
    "_adaptiveness_to_array",
    "_codon_to_index",
    "_dna_to_codon_indices",
    "_compute_cai_fast",
    "_count_dinucs_fast",
    "_BatchSwapScorer",
    "_find_site_in_sequence",
    "_get_overlapping_codons",
    "_remove_site_multicodon",
    "_find_gt_free_codons",
    "_find_ag_free_codons",
    "_organism_to_species_key",
    "_species_key_to_organism",
    "_back_translate_protein",
    "_back_translate_protein_dp",
    "_build_restriction_site_set",
    "_contains_restriction_site",
    "_count_dinucleotides",
    "_count_gts",
    "_is_unavoidable_gt",
    "_has_gt",
    "_codon_creates_boundary_gt",
    "_gt_aware_select_codon",
    "_is_in_codon_gt",
    "_eukaryote_cai_recovery",
    "_eliminate_cpg_dinucleotides",
    "_greedy_optimize",
    "_expand_iupac_site",
    "_check_predicates_via_type_system",
]

# ── Re-export moved submodules for convenient access ────────────────────
# These modules were moved from the top-level biocompiler package into
# the optimizer subpackage for better organization.

from .constraint_helpers import (  # noqa: E402
    _expand_iupac_site,
    _is_unavoidable_gt_aa,
    _gt_free_cai_ratio,
    IUPAC_EXPANSION_CAP,
)

from .gc_adjustment import (  # noqa: E402
    adjust_gc_content,
    MAX_GC_ADJUSTMENT_ITERATIONS as _MAX_GC_ADJUSTMENT_ITERATIONS_GA,
)

from .cpg_disruption import (  # noqa: E402
    disrupt_cpg_dinucleotides,
    reconcile_cpg_sites,
    MAX_CPG_DISRUPTION_ITERATIONS as _MAX_CPG_DISRUPTION_ITERATIONS_CD,
)

from .splice_elimination import (  # noqa: E402
    eliminate_cryptic_splice_sites,
    MAX_SPLICE_ELIMINATION_ITERATIONS as _MAX_SPLICE_ELIMINATION_ITERATIONS_SE,
    ELIMINATED_SITE_SCORE,
    TOP_CAI_ALTERNATIVES,
)

from .objectives import (  # noqa: E402
    ObjectiveFunction,
    cai_objective,
    cai_gc_balanced_objective,
    codon_pair_objective,
    min_max_gc_objective,
    harmonization_objective,
    make_harmonization_objective,
    resolve_objective,
    OBJECTIVE_REGISTRY,
)

from .codon_harmonization import (  # noqa: E402
    compute_rca,
    harmonize_codons,
    harmonize_with_cai_fallback,
    compute_harmonization_score,
)

from .incremental import (  # noqa: E402
    IncrementalSequenceState,
    CodonCache,
    EnzymeSiteCache,
)

from .large_sequence import (  # noqa: E402
    optimize_large_sequence,
    ProteinTooLongError,
    MAX_PROTEIN_LENGTH_DEFAULT,
)

from .multigene import (  # noqa: E402
    GeneSpec,
    MultiGeneResult,
    OperonConfig,
    optimize_multigene,
    optimize_operon,
)

from .numba_kernels import (  # noqa: E402
    count_gc,
    count_dinucleotides,
    compute_cai_kernel,
    scan_restriction_sites,
    find_all_dinucleotide_positions,
    seq_to_bytes,
    compute_cai_incremental,
    batch_codon_swap_score,
    fast_gc_window,
    fast_dinucleotide_count,
    count_gc_parallel,
    scan_restriction_sites_multi,
)

from .grammar_loader import (  # noqa: E402
    load_grammar,
    grammar_to_predicate_params,
    load_builtin_grammar,
    list_builtin_grammars,
)

from .whatif_analysis import (  # noqa: E402
    WhatIfScenario,
    WhatIfAnalyzer,
    WhatIfReport,
)

from .performance import (  # noqa: E402
    should_skip_constraint,
    batch_detect_violations,
    estimate_optimization_complexity,
    get_fast_path_config,
    warm_numba_cache,
    get_organism_data,
    clear_caches,
    should_skip_mrna_stability,
    should_skip_cpg_elimination,
    should_skip_utr_suggestions,
    FastPathConfig,
)
