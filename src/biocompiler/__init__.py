"""
BioCompiler — Machine-Verified Gene Design  (v9.2.0)

A compiler framework for human protein synthesis using intermediate
representations. Pipeline:

  Scanner → NDFST Splicing → Translation → Type Check → Certificate → Verify

All computation is DETERMINISTIC: same input always produces identical output.

v9.2.0 highlights:
  - Unified engine API: BaseEngineResult, MutationResult, BatchResult
    shared across all 6 analysis engines (ESMFold, FoldX, CamSol,
    Immunogenicity, Deimmunization, Protein Design)
  - 28-predicate type system: 12 DNA + 4 structure + 4 stability +
    4 solubility + 4 immunogenicity predicates
  - SLOT architecture: 13 core predicates (PASS/FAIL) + 19 SLOT-dependent
    predicates (always UNCERTAIN); Lean4 proof covers all 28 predicates
  - HBB full pass: all 8 optimizer predicates pass simultaneously
  - CpG reconciliation, CAI reconciliation, cross-codon coordination
"""

__version__ = "9.2.0"

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

# ═══════════════════════════════════════════════════════════════════════
# Core types & exceptions
# ═══════════════════════════════════════════════════════════════════════

from .types import (
    Verdict,
    Token,
    PositionRange,
    SpliceIsoform,
    TypeCheckResult,
    Certificate,
    three_valued_and,
    three_valued_or,
    combined_verdict,
)

# Note: three_valued_and/three_valued_or are kept for backward compatibility.
# For the full 5-valued logic, use five_valued_and/five_valued_or from .types.
from .exceptions import (
    BioCompilerError,
    EngineError,
    ESMFoldError,
    FoldXError,
    CamSolError,
    ImmunogenicityError,
    InvalidSequenceError,
    CertificateGenerationError,
    CertificateVerificationError,
    UnknownPredicateError,
    OptimizationError,
    UnsupportedOrganismError,
    InvalidProteinError,
    FileFormatError,
    SplicingError,
    MutagenesisError,
)

# ═══════════════════════════════════════════════════════════════════════
# Canonical biological constants
# ═══════════════════════════════════════════════════════════════════════

from .constants import (
    BLOSUM62, HYDROPATHY, HYDROPHOBIC_AAS, STANDARD_AAS,
    DEFAULT_ENGINE_TIMEOUT, DEFAULT_BATCH_SIZE,
    DEFAULT_SOLUBILITY_WINDOW, DEFAULT_SOLUBILITY_SMOOTHING,
    DEFAULT_MHC_PEPTIDE_LENGTH,
)

# ═══════════════════════════════════════════════════════════════════════
# Pipeline: Scanner → Splicing → Translation → Type Check
# ═══════════════════════════════════════════════════════════════════════

from .scanner import (
    validate_dna_sequence,
    gc_content,
    scan_sequence,
)

from .splicing import compute_splice_isoforms, maxent_score, score_splice_sites

from .translation import (
    translate,
    compute_cai,
    find_orfs,
)

# Type system imports
from .type_system import (
    evaluate_no_cryptic_splice,
    evaluate_splice_correct,
    evaluate_gc_in_range,
    evaluate_codon_adapted,
    evaluate_no_restriction_site,
    evaluate_in_frame,
    evaluate_no_instability_motif,
    evaluate_no_cpg_island,
    evaluate_all_predicates,
    analyze_codon_at_position,
    registry as predicate_registry,
)

# ═══════════════════════════════════════════════════════════════════════
# Certificate engine (merged from certificates.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .certificate import generate_certificate, verify_certificate, compute_certificate, format_certificate
except ImportError:
    generate_certificate = None
    verify_certificate = None
    compute_certificate = None
    format_certificate = None

# ═══════════════════════════════════════════════════════════════════════
# MaxEntScan splice scoring
# ═══════════════════════════════════════════════════════════════════════

try:
    from .maxentscan import (
        score_donor, score_acceptor, scan_splice_sites,
        max_donor_score, max_acceptor_score,
    )
except ImportError:
    score_donor = None
    score_acceptor = None
    scan_splice_sites = None
    max_donor_score = None
    max_acceptor_score = None

# ═══════════════════════════════════════════════════════════════════════
# Optimization
# ═══════════════════════════════════════════════════════════════════════

from .optimization import BioOptimizer, optimize_sequence, OptimizationResult

# ═══════════════════════════════════════════════════════════════════════
# Grammar, Export, Report
# ═══════════════════════════════════════════════════════════════════════

try:
    from .grammar_loader import load_grammar, grammar_to_predicate_params, load_builtin_grammar, list_builtin_grammars
except ImportError:
    load_grammar = None
    grammar_to_predicate_params = None
    load_builtin_grammar = None
    list_builtin_grammars = None

try:
    from .export import export_fasta, export_genbank, export_genbank_with_certificate, export_multi_fasta
except ImportError:
    export_fasta = None
    export_genbank = None
    export_genbank_with_certificate = None
    export_multi_fasta = None

try:
    from .report import generate_report
except ImportError:
    generate_report = None

# ═══════════════════════════════════════════════════════════════════════
# Benchmark (merged from comprehensive_benchmark.py & tool_comparison.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .benchmark import (
        run_benchmarks, BenchmarkReport,
        # Extended API
        REFERENCE_GENES, run_structured_benchmarks,
        format_benchmark_report_json, format_benchmark_report_text,
        # Head-to-head (merged from tool_comparison)
        run_head_to_head_benchmark, compare_tools, HeadToHeadReport,
        format_head_to_head_text, format_head_to_head_json,
        is_dna_chisel_available,
        # Comprehensive (merged from comprehensive_benchmark)
        run_comprehensive_benchmark,
    )
except ImportError:
    try:
        from .benchmark import run_benchmark as run_benchmarks
        BenchmarkReport = None
    except ImportError:
        run_benchmarks = None
        BenchmarkReport = None
    # Names only available in the extended benchmark API
    REFERENCE_GENES = None
    run_structured_benchmarks = None
    format_benchmark_report_json = None
    format_benchmark_report_text = None
    run_head_to_head_benchmark = None
    compare_tools = None
    HeadToHeadReport = None
    format_head_to_head_text = None
    format_head_to_head_json = None
    is_dna_chisel_available = None
    run_comprehensive_benchmark = None

# ═══════════════════════════════════════════════════════════════════════
# Organism data (merged from species.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .organisms import (
        OrganismDatabase, get_database,
        SPECIES, ECOLI_CAI, HUMAN_CAI, compute_cai_weights,
    )
except ImportError:
    OrganismDatabase = None
    get_database = None
    SPECIES = {}
    ECOLI_CAI = None
    HUMAN_CAI = None
    compute_cai_weights = None

# ═══════════════════════════════════════════════════════════════════════
# Tissue data
# ═══════════════════════════════════════════════════════════════════════

try:
    from .tissue_data import get_tissue_weights, list_available_tissues, add_custom_tissue
except ImportError:
    get_tissue_weights = None
    list_available_tissues = None
    add_custom_tissue = None

# ═══════════════════════════════════════════════════════════════════════
# DNA Chisel compatibility
# ═══════════════════════════════════════════════════════════════════════

try:
    from .dna_chisel_compat import compare_optimizers, run_comparative_benchmark
except ImportError:
    compare_optimizers = None
    run_comparative_benchmark = None

# ═══════════════════════════════════════════════════════════════════════
# Dataset validation
# ═══════════════════════════════════════════════════════════════════════

try:
    from .dataset_validation import run_dataset_validation, DatasetValidationReport
except ImportError:
    run_dataset_validation = None
    DatasetValidationReport = None

# ═══════════════════════════════════════════════════════════════════════
# Sequence import / BioPython compat / Jupyter
# ═══════════════════════════════════════════════════════════════════════

try:
    from .import_seq import import_fasta, import_genbank, import_sequence
except ImportError:
    import_fasta = None
    import_genbank = None
    import_sequence = None

try:
    from .biopython_compat import to_seqrecord, from_seqrecord, optimize_to_seqrecord
except ImportError:
    to_seqrecord = None
    from_seqrecord = None
    optimize_to_seqrecord = None

try:
    from .jupyter import display_sequence, display_optimization_result, display_type_check, plot_gc_content, plot_codon_usage
except ImportError:
    display_sequence = None
    display_optimization_result = None
    display_type_check = None
    plot_gc_content = None
    plot_codon_usage = None

# ═══════════════════════════════════════════════════════════════════════
# Mutagenesis
# ═══════════════════════════════════════════════════════════════════════

try:
    from .mutagenesis import (
        MutagenesisResult, AASubstitution,
        GT_MANDATORY_AAS, AG_MANDATORY_AAS,
        is_gt_mandatory, is_ag_mandatory, diagnose_optimizer_weakness,
        force_gt_free_reoptimization,
        type_directed_mutagenesis, find_unrepairable_cryptic_donors,
        find_unrepairable_cryptic_acceptors, propose_substitutions,
        apply_substitution,
    )
except ImportError:
    MutagenesisResult = None
    AASubstitution = None
    GT_MANDATORY_AAS = None
    AG_MANDATORY_AAS = None
    is_gt_mandatory = None
    is_ag_mandatory = None
    diagnose_optimizer_weakness = None
    force_gt_free_reoptimization = None
    type_directed_mutagenesis = None
    find_unrepairable_cryptic_donors = None
    find_unrepairable_cryptic_acceptors = None
    propose_substitutions = None
    apply_substitution = None

# ═══════════════════════════════════════════════════════════════════════
# Solubility (Camsol + predicates)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .camsol import (
        compute_intrinsic_solubility, compute_solubility,
        compute_structural_solubility,
        classify_solubility, find_solubility_mutations,
        generate_solubility_recommendations,
        compute_solubility_batch,
        clear_cache as camsol_clear_cache,
        SolubilityResult,
        CAMSOL_HYDROPATHY, CAMSOL_CHARGE, CAMSOL_ALPHA_HELIX,
        CAMSOL_BETA_STRAND,
    )
except ImportError:
    compute_intrinsic_solubility = None
    compute_solubility = None
    compute_structural_solubility = None
    classify_solubility = None
    find_solubility_mutations = None
    generate_solubility_recommendations = None
    compute_solubility_batch = None
    camsol_clear_cache = None
    SolubilityResult = None
    CAMSOL_HYDROPATHY = None
    CAMSOL_CHARGE = None
    CAMSOL_ALPHA_HELIX = None
    CAMSOL_BETA_STRAND = None

try:
    from .solubility_predicates import (
        evaluate_soluble_expression, evaluate_no_aggregation_prone_region,
        evaluate_charge_composition, evaluate_no_long_hydrophobic_stretch,
        compute_approximate_pI, compute_net_charge,
        find_hydrophobic_stretches, PKA_VALUES,
    )
except ImportError:
    evaluate_soluble_expression = None
    evaluate_no_aggregation_prone_region = None
    evaluate_charge_composition = None
    evaluate_no_long_hydrophobic_stretch = None
    compute_approximate_pI = None
    compute_net_charge = None
    find_hydrophobic_stretches = None
    PKA_VALUES = None

# ═══════════════════════════════════════════════════════════════════════
# Structure prediction (ESMFold — merged from esmfold_batch.py & esmfold_cache.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .esmfold import (
        ESMFoldResult, ESMFoldError, ESMFoldCache,
        is_esmfold_available, predict_structure, predict_structure_batch,
        predict_batch, predict_proteins, format_batch_report,
        parse_pdb, compute_backbone_dihedrals, classify_plddt,
        estimate_contact_map, analyze_structure,
        validate_batch_input, estimate_batch_time,
        BatchStructureRequest, BatchStructureResult,
        clear_cache as esmfold_clear_cache,
    )
except ImportError:
    ESMFoldResult = None
    # ESMFoldError is already imported unconditionally from .exceptions above
    ESMFoldCache = None
    is_esmfold_available = None
    predict_structure = None
    predict_structure_batch = None
    predict_batch = None
    predict_proteins = None
    format_batch_report = None
    parse_pdb = None
    compute_backbone_dihedrals = None
    classify_plddt = None
    estimate_contact_map = None
    analyze_structure = None
    validate_batch_input = None
    estimate_batch_time = None
    BatchStructureRequest = None
    BatchStructureResult = None
    esmfold_clear_cache = None

# ═══════════════════════════════════════════════════════════════════════
# Immunogenicity (merged from mhc_binding.py & epitope.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .immunogenicity import (
        ImmunogenicityResult,
        predict_t_cell_epitopes, predict_b_cell_epitopes,
        compute_surface_accessibility_approx,
        compute_immunogenicity, find_deimmunization_mutations,
        compute_immunogenicity_batch,
        clear_cache as immunogenicity_clear_cache,
        # MHC binding
        MHCBindingResult, MHCPredictionResult,
        predict_mhc_i_binding, predict_mhc_ii_binding,
        score_peptide_pssm, binding_score_to_ic50, classify_binding,
        predict_all as predict_mhc_binding,
        MHC_I_PSSM, MHC_II_PSSM, POPULATION_COVERAGE,
        DEFAULT_MHC_I_ALLELES, DEFAULT_MHC_II_ALLELES,
        # B-cell epitope
        EpitopeRegion, EpitopePredictionResult,
        predict_kolaskar_tongaonkar, predict_parker_hydrophilicity,
        predict_chou_fasman_beta_turn, predict_eea,
        predict_bepipred_like, predict_conformational_epitopes,
        predict_epitopes, ALL_SCALES,
        ANTIGENICITY_SCALE, PARKER_SCALE,
        CHOU_FASMAN_TURN, EMINI_SCALE,
    )
except ImportError:
    ImmunogenicityResult = None
    predict_t_cell_epitopes = None
    predict_b_cell_epitopes = None
    compute_surface_accessibility_approx = None
    compute_immunogenicity = None
    find_deimmunization_mutations = None
    compute_immunogenicity_batch = None
    immunogenicity_clear_cache = None
    MHCBindingResult = None
    MHCPredictionResult = None
    predict_mhc_i_binding = None
    predict_mhc_ii_binding = None
    score_peptide_pssm = None
    binding_score_to_ic50 = None
    classify_binding = None
    predict_mhc_binding = None
    MHC_I_PSSM = None
    MHC_II_PSSM = None
    POPULATION_COVERAGE = None
    DEFAULT_MHC_I_ALLELES = None
    DEFAULT_MHC_II_ALLELES = None
    EpitopeRegion = None
    EpitopePredictionResult = None
    predict_kolaskar_tongaonkar = None
    predict_parker_hydrophilicity = None
    predict_chou_fasman_beta_turn = None
    predict_eea = None
    predict_bepipred_like = None
    predict_conformational_epitopes = None
    predict_epitopes = None
    ALL_SCALES = None
    ANTIGENICITY_SCALE = None
    PARKER_SCALE = None
    CHOU_FASMAN_TURN = None
    EMINI_SCALE = None

# Deprecated aliases removed in v7.5.0:
# MHC_I_PREFERENCES → use MHC_I_PSSM
# MHC_II_PREFERENCES → use MHC_II_PSSM
# ANTIGENICITY_PROPENSITY → use ANTIGENICITY_SCALE

# ═══════════════════════════════════════════════════════════════════════
# Deimmunization
# ═══════════════════════════════════════════════════════════════════════

try:
    from .deimmunization import (
        DeimmunizationResult,
        EpitopeMutation,
        deimmunize,
        find_epitope_disrupting_mutations,
        rank_deimmunization_mutations,
        validate_deimmunized_protein,
        compute_mutation_impact,
    )
except ImportError:
    DeimmunizationResult = None
    EpitopeMutation = None
    deimmunize = None
    find_epitope_disrupting_mutations = None
    rank_deimmunization_mutations = None
    validate_deimmunized_protein = None
    compute_mutation_impact = None

# ═══════════════════════════════════════════════════════════════════════
# Structure (consolidated from structure/, structure_quality, structure_predicates, structure_report)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .structure import (
        # Parser
        Atom, Residue, Chain, ProteinStructure,
        parse_pdb, parse_pdb_file, compute_dihedral,
        compute_ramachandran, secondary_structure_estimate,
        THREE_TO_ONE, ONE_TO_THREE,
        # Quality
        StructureQualityReport,
        assess_plddt, assess_ramachandran,
        compute_clash_score, compute_packing_density,
        compute_exposed_hydrophobic, compute_structure_quality,
        find_low_confidence_regions, compute_sasa_approximation,
        KYTE_DOOLITTLE, VDW_RADII,
        # Predicates
        evaluate_structure_confidence, evaluate_no_misfolding_risk,
        evaluate_correct_fold_topology, evaluate_no_unexpected_interaction,
        # Report
        ProteinAssessmentReport, assess_protein,
        format_assessment_text, format_assessment_json, format_assessment_html,
        generate_recommendations, compute_overall_verdict,
        plot_plddt_bar_svg, plot_solubility_profile_svg,
    )
except ImportError:
    Atom = None
    Residue = None
    Chain = None
    ProteinStructure = None
    # parse_pdb may already be defined from .esmfold; don't clobber it
    parse_pdb_file = None
    compute_dihedral = None
    compute_ramachandran = None
    secondary_structure_estimate = None
    THREE_TO_ONE = None
    ONE_TO_THREE = None
    StructureQualityReport = None
    assess_plddt = None
    assess_ramachandran = None
    compute_clash_score = None
    compute_packing_density = None
    compute_exposed_hydrophobic = None
    compute_structure_quality = None
    find_low_confidence_regions = None
    compute_sasa_approximation = None
    KYTE_DOOLITTLE = None
    VDW_RADII = None
    evaluate_structure_confidence = None
    evaluate_no_misfolding_risk = None
    evaluate_correct_fold_topology = None
    evaluate_no_unexpected_interaction = None
    ProteinAssessmentReport = None
    assess_protein = None
    format_assessment_text = None
    format_assessment_json = None
    format_assessment_html = None
    generate_recommendations = None
    compute_overall_verdict = None
    plot_plddt_bar_svg = None
    plot_solubility_profile_svg = None

# ═══════════════════════════════════════════════════════════════════════
# FoldX stability (merged from foldx_mutations.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .foldx import (
        FoldXResult, FoldXError, FoldXCache,
        StabilityLandscape, ConservationScore,
        is_foldx_available, run_foldx_stability, run_foldx_repair,
        run_foldx_mutation, empirical_stability,
        run_stability_batch,
        clear_cache as foldx_clear_cache,
        scan_mutations, find_stabilizing_mutations,
        scan_all_mutations, scan_position, compute_conservation,
        find_compensatory_mutations, rank_positions_by_mutability,
        identify_hotspot_regions,
    )
except ImportError:
    FoldXResult = None
    # FoldXError is already imported unconditionally from .exceptions above
    FoldXCache = None
    StabilityLandscape = None
    ConservationScore = None
    is_foldx_available = None
    run_foldx_stability = None
    run_foldx_repair = None
    run_foldx_mutation = None
    empirical_stability = None
    run_stability_batch = None
    foldx_clear_cache = None
    scan_mutations = None
    find_stabilizing_mutations = None
    scan_all_mutations = None
    scan_position = None
    compute_conservation = None
    find_compensatory_mutations = None
    rank_positions_by_mutability = None
    identify_hotspot_regions = None

# ═══════════════════════════════════════════════════════════════════════
# Engine base — unified types for all analysis engines
# ═══════════════════════════════════════════════════════════════════════

try:
    from .engine_base import (
        EngineResult, BaseEngineResult, MutationResult, BatchResult,
        EngineTimer, EngineConfig, validate_protein_sequence,
        classify_score,
    )
except ImportError:
    EngineResult = None
    BaseEngineResult = None
    MutationResult = None
    BatchResult = None
    EngineTimer = None
    EngineConfig = None
    validate_protein_sequence = None
    classify_score = None

# ═══════════════════════════════════════════════════════════════════════
# Protein design
# ═══════════════════════════════════════════════════════════════════════

try:
    from .protein_design import (
        DesignResult, DesignConstraints,
        design_thermostable, design_soluble,
        design_low_immunogenicity, design_multi_objective,
        score_mutation, find_disulfide_opportunities,
        find_proline_substitution_sites,
    )
except ImportError:
    DesignResult = None
    DesignConstraints = None
    design_thermostable = None
    design_soluble = None
    design_low_immunogenicity = None
    design_multi_objective = None
    score_mutation = None
    find_disulfide_opportunities = None
    find_proline_substitution_sites = None

# ═══════════════════════════════════════════════════════════════════════
# CSP/SMT Solver (constraint-based gene optimization)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .solver import CSPSolver, SolverConfig, SolverResult, SolverBackend
    from .solver.dispatch import solve_with_csp, is_solver_available, csp_optimize
except ImportError:
    CSPSolver = None
    SolverConfig = None
    SolverResult = None
    SolverBackend = None
    solve_with_csp = None
    is_solver_available = None
    csp_optimize = None

# ═══════════════════════════════════════════════════════════════════════
# ViennaRNA (mRNA secondary structure prediction)
# ═══════════════════════════════════════════════════════════════════════

try:
    from . import viennarna
    from . import viennarna_fallback
except ImportError:
    viennarna = None  # type: ignore[assignment]
    viennarna_fallback = None  # type: ignore[assignment]

# ═══════════════════════════════════════════════════════════════════════
# MHCflurry (offline MHC-I binding prediction)
# ═══════════════════════════════════════════════════════════════════════

try:
    from . import mhcflurry_adapter
    from . import mhcflurry_population
except ImportError:
    mhcflurry_adapter = None  # type: ignore[assignment]
    mhcflurry_population = None  # type: ignore[assignment]

# ═══════════════════════════════════════════════════════════════════════
# Public API — organized by domain
# ═══════════════════════════════════════════════════════════════════════

__all__ = [
    # ── Version ──────────────────────────────────────────────
    "__version__",

    # ── Core types ───────────────────────────────────────────
    "Verdict", "Token", "PositionRange", "SpliceIsoform",
    "TypeCheckResult", "Certificate",
    "three_valued_and", "three_valued_or", "combined_verdict",

    # ── Exceptions ───────────────────────────────────────────
    "BioCompilerError", "EngineError",
    "ESMFoldError", "FoldXError", "CamSolError", "ImmunogenicityError",
    "InvalidSequenceError",
    "CertificateGenerationError", "CertificateVerificationError",
    "UnknownPredicateError", "OptimizationError",
    "UnsupportedOrganismError", "InvalidProteinError",
    "FileFormatError", "SplicingError", "MutagenesisError",

    # ── Canonical constants ──────────────────────────────────
    "BLOSUM62", "HYDROPATHY", "HYDROPHOBIC_AAS", "STANDARD_AAS",
    "DEFAULT_ENGINE_TIMEOUT", "DEFAULT_BATCH_SIZE",
    "DEFAULT_SOLUBILITY_WINDOW", "DEFAULT_SOLUBILITY_SMOOTHING",
    "DEFAULT_MHC_PEPTIDE_LENGTH",

    # ── Scanner ──────────────────────────────────────────────
    "validate_dna_sequence", "gc_content", "scan_sequence",

    # ── Splicing ─────────────────────────────────────────────
    "compute_splice_isoforms", "maxent_score", "score_splice_sites",

    # ── Translation ──────────────────────────────────────────
    "translate", "compute_cai", "find_orfs",

    # ── Type system & predicates ─────────────────────────────
    "evaluate_no_cryptic_splice", "evaluate_splice_correct",
    "evaluate_gc_in_range", "evaluate_codon_adapted",
    "evaluate_no_restriction_site", "evaluate_in_frame",
    "evaluate_no_instability_motif", "evaluate_no_cpg_island",
    "evaluate_all_predicates", "analyze_codon_at_position",
    "predicate_registry",

    # ── Certificate ──────────────────────────────────────────
    "generate_certificate", "verify_certificate",
    "compute_certificate", "format_certificate",

    # ── MaxEntScan ───────────────────────────────────────────
    "score_donor", "score_acceptor", "scan_splice_sites",
    "max_donor_score", "max_acceptor_score",

    # ── Optimization ─────────────────────────────────────────
    "BioOptimizer", "optimize_sequence", "OptimizationResult",

    # ── Grammar ──────────────────────────────────────────────
    "load_grammar", "grammar_to_predicate_params",
    "load_builtin_grammar", "list_builtin_grammars",

    # ── Export ───────────────────────────────────────────────
    "export_fasta", "export_genbank", "export_genbank_with_certificate", "export_multi_fasta",

    # ── Report ───────────────────────────────────────────────
    "generate_report",

    # ── Benchmark ────────────────────────────────────────────
    "run_benchmarks", "BenchmarkReport",
    "REFERENCE_GENES", "run_structured_benchmarks",
    "format_benchmark_report_json", "format_benchmark_report_text",
    "run_head_to_head_benchmark", "compare_tools", "HeadToHeadReport",
    "format_head_to_head_text", "format_head_to_head_json",
    "is_dna_chisel_available",
    "run_comprehensive_benchmark",

    # ── Organisms ────────────────────────────────────────────
    "OrganismDatabase", "get_database",
    "SPECIES", "ECOLI_CAI", "HUMAN_CAI", "compute_cai_weights",

    # ── Tissue data ──────────────────────────────────────────
    "get_tissue_weights", "list_available_tissues", "add_custom_tissue",

    # ── DNA Chisel compat ────────────────────────────────────
    "compare_optimizers", "run_comparative_benchmark",

    # ── Dataset validation ───────────────────────────────────
    "run_dataset_validation", "DatasetValidationReport",

    # ── Sequence import / BioPython / Jupyter ────────────────
    "import_fasta", "import_genbank", "import_sequence",
    "to_seqrecord", "from_seqrecord", "optimize_to_seqrecord",
    "display_sequence", "display_optimization_result", "display_type_check",
    "plot_gc_content", "plot_codon_usage",

    # ── Mutagenesis ──────────────────────────────────────────
    "MutagenesisResult", "AASubstitution",
    "GT_MANDATORY_AAS", "AG_MANDATORY_AAS",
    "is_gt_mandatory", "is_ag_mandatory", "diagnose_optimizer_weakness",
    "force_gt_free_reoptimization",
    "type_directed_mutagenesis", "find_unrepairable_cryptic_donors",
    "find_unrepairable_cryptic_acceptors", "propose_substitutions",
    "apply_substitution",

    # ── Engine base types ────────────────────────────────────
    "EngineResult", "BaseEngineResult", "MutationResult", "BatchResult",
    "EngineTimer", "EngineConfig", "validate_protein_sequence",
    "classify_score",

    # ── Solubility ───────────────────────────────────────────
    "compute_intrinsic_solubility", "compute_solubility",
    "compute_structural_solubility",
    "classify_solubility", "find_solubility_mutations",
    "generate_solubility_recommendations",
    "compute_solubility_batch", "camsol_clear_cache",
    "SolubilityResult",
    "CAMSOL_HYDROPATHY", "CAMSOL_CHARGE", "CAMSOL_ALPHA_HELIX",
    "CAMSOL_BETA_STRAND",
    "evaluate_soluble_expression", "evaluate_no_aggregation_prone_region",
    "evaluate_charge_composition", "evaluate_no_long_hydrophobic_stretch",
    "compute_approximate_pI", "compute_net_charge",
    "find_hydrophobic_stretches", "PKA_VALUES",

    # ── ESMFold structure prediction ─────────────────────────
    "ESMFoldResult", "ESMFoldCache",
    "is_esmfold_available", "predict_structure", "predict_structure_batch",
    "predict_batch", "predict_proteins", "format_batch_report",
    "compute_backbone_dihedrals", "classify_plddt",
    "estimate_contact_map", "analyze_structure",
    "validate_batch_input", "estimate_batch_time",
    "BatchStructureRequest", "BatchStructureResult",
    "esmfold_clear_cache",

    # ── Immunogenicity & MHC binding ─────────────────────────
    "ImmunogenicityResult",
    "predict_t_cell_epitopes", "predict_b_cell_epitopes",
    "compute_surface_accessibility_approx",
    "compute_immunogenicity", "find_deimmunization_mutations",
    "compute_immunogenicity_batch", "immunogenicity_clear_cache",
    "MHCBindingResult", "MHCPredictionResult",
    "predict_mhc_i_binding", "predict_mhc_ii_binding",
    "score_peptide_pssm", "binding_score_to_ic50", "classify_binding",
    "predict_mhc_binding",
    "MHC_I_PSSM", "MHC_II_PSSM", "POPULATION_COVERAGE",
    "DEFAULT_MHC_I_ALLELES", "DEFAULT_MHC_II_ALLELES",
    "EpitopeRegion", "EpitopePredictionResult",
    "predict_kolaskar_tongaonkar", "predict_parker_hydrophilicity",
    "predict_chou_fasman_beta_turn", "predict_eea",
    "predict_bepipred_like", "predict_conformational_epitopes",
    "predict_epitopes", "ALL_SCALES",
    "ANTIGENICITY_SCALE", "PARKER_SCALE",
    "CHOU_FASMAN_TURN", "EMINI_SCALE",

    # ── Deimmunization ───────────────────────────────────────
    "DeimmunizationResult", "EpitopeMutation",
    "deimmunize", "find_epitope_disrupting_mutations",
    "rank_deimmunization_mutations", "validate_deimmunized_protein",
    "compute_mutation_impact",

    # ── Structure (consolidated) ──────────────────────────────
    # Parser
    "Atom", "Residue", "Chain", "ProteinStructure",
    "parse_pdb", "parse_pdb_file", "compute_dihedral",
    "compute_ramachandran", "secondary_structure_estimate",
    "THREE_TO_ONE", "ONE_TO_THREE",
    # Quality
    "StructureQualityReport",
    "assess_plddt", "assess_ramachandran",
    "compute_clash_score", "compute_packing_density",
    "compute_exposed_hydrophobic", "compute_structure_quality",
    "find_low_confidence_regions", "compute_sasa_approximation",
    "KYTE_DOOLITTLE", "VDW_RADII",
    # Predicates
    "evaluate_structure_confidence", "evaluate_no_misfolding_risk",
    "evaluate_correct_fold_topology", "evaluate_no_unexpected_interaction",
    # Report
    "ProteinAssessmentReport", "assess_protein",
    "format_assessment_text", "format_assessment_json", "format_assessment_html",
    "generate_recommendations", "compute_overall_verdict",
    "plot_plddt_bar_svg", "plot_solubility_profile_svg",

    # ── FoldX stability ──────────────────────────────────────
    "FoldXResult", "FoldXCache",
    "StabilityLandscape", "ConservationScore",
    "is_foldx_available", "run_foldx_stability", "run_foldx_repair",
    "run_foldx_mutation", "empirical_stability",
    "run_stability_batch", "foldx_clear_cache",
    "scan_mutations", "find_stabilizing_mutations",
    "scan_all_mutations", "scan_position", "compute_conservation",
    "find_compensatory_mutations", "rank_positions_by_mutability",
    "identify_hotspot_regions",

    # ── Protein design ───────────────────────────────────────
    "DesignResult", "DesignConstraints",
    "design_thermostable", "design_soluble",
    "design_low_immunogenicity", "design_multi_objective",
    "score_mutation", "find_disulfide_opportunities",
    "find_proline_substitution_sites",

    # ── CSP/SMT Solver ───────────────────────────────────────────
    "CSPSolver", "SolverConfig", "SolverResult", "SolverBackend",
    "solve_with_csp", "is_solver_available", "csp_optimize",

    # ── ViennaRNA ────────────────────────────────────────────────
    "viennarna", "viennarna_fallback",

    # ── MHCflurry ────────────────────────────────────────────────
    "mhcflurry_adapter", "mhcflurry_population",
]
