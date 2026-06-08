"""
biocompiler.optimizer — Modular optimization subpackage.

This package decomposes the former monolithic optimization.py into
focused submodules while maintaining full backward compatibility.

Submodules:
    cai                  — CAI computation and batch scoring
    constraints          — Restriction site removal, GT/AG helpers, back-translation
    greedy               — Greedy optimizer loop and splice donor scoring
    mfe_optimization     — MFE minimization (greedy 5', full-gene SA, LinearDesign, IPknot)
    pipeline             — BioOptimizer class, optimize_sequence(), batch_optimize()
    ribosome_simulation  — TASEP simulation, RQC/NGD detection, mRNA structure modulation
    nucleosome           — Nucleosome positioning prediction (Segal PSSM, NuPoP, Teif-Percus, Enformer)
    dna_shape            — DNA shape prediction (MGW, HelT, ProT, Roll) and damage susceptibility (Olson 1998, dnacurve, Deep DNAshape)
    accessibility        — mRNA accessibility (unpaired probability) via ViennaRNA partition function
    chromatin_context    — Chromatin state prediction (euchromatin/heterochromatin, DNase, repair accessibility, Sei, Enformer)
    eukaryotic_decay     — Eukaryotic mRNA decay pathway modeling (XRN1, CCR4-NOT, DCP1/DCP2, NMD, NGD, NSD)
    epitranscriptomics   — Epitranscriptomic modification site detection (m6A, m5C, Ψ, m1A, 2'-O-Me, m6Am)
    utils                — ConvergenceTracker, OptimizationResult, protein_to_aa_list()
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
    # Ribosome simulation (TASEP, RQC/NGD)
    "simulate_tasep_gillespie",
    "simulate_tasep_ensemble",
    "detect_rqc_signals",
    "compute_mrna_structure_elongation_mod",
    "predict_secondary_structure_modern",
    # MFE optimization (LinearDesign, SA, pseudoknots)
    "optimize_with_lineardesign",
    "optimize_mfe_simulated_annealing",
    "detect_pseudoknots",
    # RNA degradation (NMD, decay, epitranscriptomics)
    "detect_nmd_triggers",
    "detect_eukaryotic_decay_signals",
    "detect_epitranscriptomic_marks",
    "compute_mrna_accessibility",
    # DNA damage (mutation risk, COSMIC)
    "compute_net_mutation_risk",
    "compute_cosmic_context_weights",
    "infer_gene_strand",
    # Nucleosome positioning (Kaplan, Segal, NuPoP, Enformer)
    "score_kaplan_pssm",
    "score_segal_pssm",
    "predict_nucleosome_nupop",
    "predict_occupancy_with_exclusion",
    "predict_histone_marks_enformer",
    # Chromatin context (euchromatin/heterochromatin, DNase, repair accessibility, Sei, Enformer)
    "predict_chromatin_state",
    "compute_gc_ishchz_score",
    "compute_cpg_density_score",
    "predict_dnase_accessibility",
    "compute_repair_accessibility",
    "predict_chromatin_sei",
    "predict_chromatin_enformer",
    # Ligand binding v2 (SMILES, 3D conformer, Vina docking)
    "parse_smiles_features_rdkit",
    "generate_3d_conformer",
    "dock_ligand_vina",
    "compute_interaction_fingerprint",
    # Eukaryotic decay pathways (XRN1, CCR4-NOT, DCP1/DCP2, NMD, NGD, NSD)
    "DecayPathway",
    "DecaySignal",
    "predict_decay_rate",
    "detect_decay_signals",
    "optimize_decay_resistance",
    "estimate_halflife",
    # mRNA accessibility (ViennaRNA partition function, GC heuristic fallback)
    "compute_accessibility",
    "compute_accessibility_windows",
    "compute_codon_accessibility",
    "compute_5prime_accessibility",
    "compute_mirna_site_accessibility",
    "adjust_severity_for_accessibility",
    "AccessibilityWindowResult",
    "HAS_VIENNARNA",
    # DNA shape prediction (MGW, HelT, ProT, Roll, damage susceptibility)
    "DNAShapeProfile",
    "compute_minor_groove_width",
    "compute_helix_twist",
    "compute_propeller_twist",
    "compute_roll",
    "compute_dna_shape_profile",
    "compute_damage_susceptibility_from_shape",
    "compute_shape_dnacurve",
    "compute_shape_deep_dnashape",
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

# ── Phase 8: New submodule exports ────────────────────────────────────
# These imports use try/except to allow graceful degradation when the
# new Phase 8 submodules are not yet installed.

try:
    from .rna_degradation import (  # noqa: E402
        detect_nmd_triggers,
        detect_eukaryotic_decay_signals,
        detect_epitranscriptomic_marks,
        compute_mrna_accessibility,
    )
except ImportError:
    pass

try:
    from .dna_damage import (  # noqa: E402
        compute_net_mutation_risk,
        compute_cosmic_context_weights,
        infer_gene_strand,
    )
except ImportError:
    pass

try:
    from .ligand_binding_v2 import (  # noqa: E402
        parse_smiles_features_rdkit,
        generate_3d_conformer,
        dock_ligand_vina,
        compute_interaction_fingerprint,
    )
except ImportError:
    pass

from .ribosome_simulation import (  # noqa: E402
    simulate_tasep_gillespie,
    simulate_tasep_ensemble,
    detect_rqc_signals,
    compute_mrna_structure_elongation_mod,
    predict_secondary_structure_modern,
)

from .mfe_optimization import (  # noqa: E402
    optimize_mfe,
    optimize_mfe_simulated_annealing,
    optimize_with_lineardesign,
    detect_pseudoknots,
    DEFAULT_5PRIME_WINDOW,
)

from .nucleosome import (  # noqa: E402
    NUCLEOSOME_SIZE,
    DEFAULT_STEP,
    FINE_STEP,
    DEFAULT_CHEMICAL_POTENTIAL,
    HELICAL_PERIOD,
    DINUCLEOTIDES,
    NUPOP_SPECIES,
    score_kaplan_pssm,
    score_segal_pssm,  # backward compatibility alias for score_kaplan_pssm
    predict_nucleosome_nupop,
    predict_occupancy_with_exclusion,
    predict_histone_marks_enformer,
    predict_nucleosome_occupancy,
    find_nucleosome_positions,
)

# ── COSMIC mutational signatures ──────────────────────────────────────
try:
    from .cosmic_signatures import (  # noqa: E402
        SBS_SIGNATURES,
        TRINUCLEOTIDE_CONTEXTS,
        get_trinucleotide_context,
        compute_signature_weight,
        scan_signature_hotspots,
        assign_signatures,
        compute_damage_susceptibility_profile,
        assign_signatures_sigprofiler,
    )
except ImportError:
    pass

try:
    from .chromatin_context import (  # noqa: E402
        predict_chromatin_state,
        compute_gc_ishchz_score,
        compute_cpg_density_score,
        predict_dnase_accessibility,
        compute_repair_accessibility,
        predict_chromatin_sei,
        predict_chromatin_enformer,
    )
except ImportError:
    pass

try:
    from .epitranscriptomics import (  # noqa: E402
        EpitranscriptomicSite,
        detect_m6a_sites,
        detect_m5c_sites,
        detect_pseudouridine_sites,
        detect_m1a_sites,
        detect_2om_sites,
        detect_m6am_sites,
        detect_all_epitranscriptomic_marks,
        assess_stability_impact,
        get_modification_function,
        MODIFICATION_FUNCTIONS,
    )
except ImportError:
    pass

try:
    from .eukaryotic_decay import (  # noqa: E402
        DecayPathway,
        DecaySignal,
        predict_decay_rate,
        detect_decay_signals,
        optimize_decay_resistance,
        estimate_halflife,
    )
except ImportError:
    pass

try:
    from .accessibility import (  # noqa: E402
        compute_accessibility,
        compute_accessibility_windows,
        compute_codon_accessibility,
        compute_5prime_accessibility,
        compute_mirna_site_accessibility,
        adjust_severity_for_accessibility,
        AccessibilityWindowResult,
        HAS_VIENNARNA,
    )
except ImportError:
    pass

try:
    from .dna_shape import (  # noqa: E402
        DNAShapeProfile,
        compute_minor_groove_width,
        compute_helix_twist,
        compute_propeller_twist,
        compute_roll,
        compute_dna_shape_profile,
        compute_damage_susceptibility_from_shape,
        compute_shape_dnacurve,
        compute_shape_deep_dnashape,
    )
except ImportError:
    pass
