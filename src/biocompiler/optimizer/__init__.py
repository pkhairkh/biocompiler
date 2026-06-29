"""
biocompiler.optimizer — Modular optimization subpackage.

This package decomposes the former monolithic optimization.py into
focused submodules while maintaining full backward compatibility.

Submodules:
    cai         — CAI computation and batch scoring
    constraints — Restriction site removal, GT/AG helpers, back-translation
    pipeline    — BioOptimizer class, optimize_sequence(), batch_optimize()
    utils       — ConvergenceTracker, OptimizationResult, protein_to_aa_list()

Second-pass cleanup (Task SP2): top-level imports of the slow-path stack
(greedy, hybrid_*, strategy_*, incremental, large_sequence, multigene,
whatif_analysis, performance, codon_harmonization) have been removed in
lockstep with the file deletions. The fast-path modules (integrated_optimizer,
pipeline_core, pipeline_certification, pipeline_paths) and the kept
constraint/CAI/utility modules remain.
"""

# Re-export all public symbols for backward compatibility.
# Users can do:  from biocompiler.optimizer import BioOptimizer

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

from .pipeline_core import (
    optimize_sequence,
    batch_optimize,
    BioOptimizer,
)

from .two_pass import optimize_two_pass

__all__ = [
    "OptimizationResult",
    "FullConstructResult",
    "optimize_sequence",
    "optimize_two_pass",
    "batch_optimize",
    "protein_to_aa_list",
    "BioOptimizer",
    "ConvergenceTracker",
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
    "MAX_MIRNA_AVOIDANCE_ITERATIONS",
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
    "_expand_iupac_site",
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

# gc_adjustment.py, cpg_disruption.py, splice_elimination.py were deleted
# in Task W4 / H15 — the integrated optimizer has its own inline versions
# of GC adjustment, CpG disruption, and splice elimination, so the
# stand-alone modules (906 LOC) were dead production code. The constants
# MAX_GC_ADJUSTMENT_ITERATIONS, MAX_CPG_DISRUPTION_ITERATIONS, and
# MAX_SPLICE_ELIMINATION_ITERATIONS remain defined in utils.py (imported
# at the top of this file). (Task W4 / H15.)

from .mirna_avoidance import (  # noqa: E402
    eliminate_mirna_binding_sites as eliminate_mirna_binding_sites_active,
    MAX_MIRNA_AVOIDANCE_ITERATIONS,
)

from biocompiler.optimizer.objectives import (  # noqa: E402
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

from biocompiler.optimizer.numba_kernels import (  # noqa: E402
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

from biocompiler.optimizer.grammar_loader import (  # noqa: E402
    load_grammar,
    grammar_to_predicate_params,
    load_builtin_grammar,
    list_builtin_grammars,
)

__all__.extend([
    # miRNA active avoidance (mirna_elimination.py was retired in W7-b;
    # mirna_avoidance.py is now the sole implementation)
    "eliminate_mirna_binding_sites_active",
])
