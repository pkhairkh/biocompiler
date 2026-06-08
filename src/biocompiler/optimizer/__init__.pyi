"""
Type stubs for biocompiler.optimizer package.
"""

from typing import Any

from .cai import (
    HAS_NUMBA as HAS_NUMBA,
    _HAS_NUMBA as _HAS_NUMBA,
    _adaptiveness_to_array as _adaptiveness_to_array,
    _codon_to_index as _codon_to_index,
    _dna_to_codon_indices as _dna_to_codon_indices,
    _compute_cai_fast as _compute_cai_fast,
    _count_dinucs_fast as _count_dinucs_fast,
    _BatchSwapScorer as _BatchSwapScorer,
)

from .constraints import (
    _find_site_in_sequence as _find_site_in_sequence,
    _get_overlapping_codons as _get_overlapping_codons,
    _remove_site_multicodon as _remove_site_multicodon,
    _find_gt_free_codons as _find_gt_free_codons,
    _find_ag_free_codons as _find_ag_free_codons,
    _organism_to_species_key as _organism_to_species_key,
    _species_key_to_organism as _species_key_to_organism,
    _back_translate_protein as _back_translate_protein,
    _back_translate_protein_dp as _back_translate_protein_dp,
    _build_restriction_site_set as _build_restriction_site_set,
    _contains_restriction_site as _contains_restriction_site,
    _count_dinucleotides as _count_dinucleotides,
    _count_gts as _count_gts,
    _is_unavoidable_gt as _is_unavoidable_gt,
    _has_gt as _has_gt,
    _codon_creates_boundary_gt as _codon_creates_boundary_gt,
)

from .greedy import (
    score_splice_donor_potential as score_splice_donor_potential,
    _gt_aware_select_codon as _gt_aware_select_codon,
    _is_in_codon_gt as _is_in_codon_gt,
    _eukaryote_cai_recovery as _eukaryote_cai_recovery,
    _eliminate_cpg_dinucleotides as _eliminate_cpg_dinucleotides,
    _greedy_optimize as _greedy_optimize,
    _expand_iupac_site as _expand_iupac_site,
    _check_predicates_via_type_system as _check_predicates_via_type_system,
)

from .pipeline import (
    optimize_sequence as optimize_sequence,
    batch_optimize as batch_optimize,
    BioOptimizer as BioOptimizer,
)

from .utils import (
    ConvergenceTracker as ConvergenceTracker,
    OptimizationResult as OptimizationResult,
    FullConstructResult as FullConstructResult,
    protein_to_aa_list as protein_to_aa_list,
    MAX_RESTRICTION_SITE_ITERATIONS as MAX_RESTRICTION_SITE_ITERATIONS,
    MAX_IUPAC_SITE_ITERATIONS as MAX_IUPAC_SITE_ITERATIONS,
    MAX_ATTTA_MOTIF_ITERATIONS as MAX_ATTTA_MOTIF_ITERATIONS,
    MAX_T_RUN_ITERATIONS as MAX_T_RUN_ITERATIONS,
    MAX_GC_ADJUSTMENT_ITERATIONS as MAX_GC_ADJUSTMENT_ITERATIONS,
    MAX_SPLICE_ELIMINATION_ITERATIONS as MAX_SPLICE_ELIMINATION_ITERATIONS,
    MAX_CPG_DISRUPTION_ITERATIONS as MAX_CPG_DISRUPTION_ITERATIONS,
    DEFAULT_MAX_ITERATIONS as DEFAULT_MAX_ITERATIONS,
    CONVERGENCE_IMPROVEMENT_THRESHOLD as CONVERGENCE_IMPROVEMENT_THRESHOLD,
    CONVERGENCE_PATIENCE as CONVERGENCE_PATIENCE,
    OSCILLATION_WINDOW as OSCILLATION_WINDOW,
)

SPLICE_DONOR_POTENTIAL_THRESHOLD: float
EUKARYOTE_CAI_GT_COST_THRESHOLD: float
GT_CAI_LOG_ADAPTIVENESS_COST: float

__all__: list[str]
