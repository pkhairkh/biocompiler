"""
BioCompiler — Machine-Verified Gene Design

A compiler framework for human protein synthesis using intermediate
representations. Pipeline:

  Scanner → NDFST Splicing → Translation → Type Check → Certificate → Verify

All computation is DETERMINISTIC: same input always produces identical output.
"""

__version__ = "7.4.0"

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
from .exceptions import (
    BioCompilerError,
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

from .constants import BLOSUM62, HYDROPATHY, HYDROPHOBIC_AAS, STANDARD_AAS

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
    pass

# ═══════════════════════════════════════════════════════════════════════
# MaxEntScan splice scoring
# ═══════════════════════════════════════════════════════════════════════

try:
    from .maxentscan import (
        score_donor, score_acceptor, scan_splice_sites,
        max_donor_score, max_acceptor_score,
    )
except ImportError:
    pass

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
    pass

try:
    from .export import export_fasta, export_genbank, export_genbank_with_certificate, export_multi_fasta
except ImportError:
    pass

try:
    from .report import generate_report
except ImportError:
    pass

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

# ═══════════════════════════════════════════════════════════════════════
# Organism data (merged from species.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .organisms import (
        OrganismDatabase, get_database,
        SPECIES, ECOLI_CAI, HUMAN_CAI, compute_cai_weights,
    )
except ImportError:
    pass

# ═══════════════════════════════════════════════════════════════════════
# Tissue data
# ═══════════════════════════════════════════════════════════════════════

try:
    from .tissue_data import get_tissue_weights, list_available_tissues, add_custom_tissue
except ImportError:
    pass

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
    pass

# ═══════════════════════════════════════════════════════════════════════
# Sequence import / BioPython compat / Jupyter
# ═══════════════════════════════════════════════════════════════════════

try:
    from .import_seq import import_fasta, import_genbank, import_sequence
except ImportError:
    pass

try:
    from .biopython_compat import to_seqrecord, from_seqrecord, optimize_to_seqrecord
except ImportError:
    pass

try:
    from .jupyter import display_sequence, display_optimization_result, display_type_check, plot_gc_content, plot_codon_usage
except ImportError:
    pass

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
    pass

# ═══════════════════════════════════════════════════════════════════════
# Solubility (Camsol + predicates)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .camsol import (
        compute_intrinsic_solubility, compute_solubility,
        compute_structural_solubility,
        classify_solubility, find_solubility_mutations,
        generate_solubility_recommendations,
        SolubilityResult,
        CAMSOL_HYDROPATHY, CAMSOL_CHARGE, CAMSOL_ALPHA_HELIX,
        CAMSOL_BETA_STRAND,
    )
except ImportError:
    pass

try:
    from .solubility_predicates import (
        evaluate_soluble_expression, evaluate_no_aggregation_prone_region,
        evaluate_charge_composition, evaluate_no_long_hydrophobic_stretch,
        compute_approximate_pI, compute_net_charge,
        find_hydrophobic_stretches, PKA_VALUES,
    )
except ImportError:
    pass

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
    )
except ImportError:
    pass

# ═══════════════════════════════════════════════════════════════════════
# Immunogenicity (merged from mhc_binding.py & epitope.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .immunogenicity import (
        ImmunogenicityResult,
        predict_t_cell_epitopes, predict_b_cell_epitopes,
        compute_surface_accessibility_approx,
        compute_immunogenicity, find_deimmunization_mutations,
        MHC_I_PREFERENCES, MHC_II_PREFERENCES,
        ANTIGENICITY_PROPENSITY,
        # MHC binding (merged from mhc_binding.py)
        MHCBindingResult, MHCPredictionResult,
        predict_mhc_i_binding, predict_mhc_ii_binding,
        score_peptide_pssm, binding_score_to_ic50, classify_binding,
        predict_all as predict_mhc_binding,
        MHC_I_PSSM, MHC_II_PSSM, POPULATION_COVERAGE,
        DEFAULT_MHC_I_ALLELES, DEFAULT_MHC_II_ALLELES,
        # B-cell epitope (merged from epitope.py)
        EpitopeRegion, EpitopePredictionResult,
        predict_kolaskar_tongaonkar, predict_parker_hydrophilicity,
        predict_chou_fasman_beta_turn, predict_eea,
        predict_bepipred_like, predict_conformational_epitopes,
        predict_epitopes, ALL_SCALES,
        ANTIGENICITY_SCALE, PARKER_SCALE,
        CHOU_FASMAN_TURN, EMINI_SCALE,
    )
except ImportError:
    pass

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
    pass

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
    pass

# ═══════════════════════════════════════════════════════════════════════
# FoldX stability (merged from foldx_mutations.py)
# ═══════════════════════════════════════════════════════════════════════

try:
    from .foldx import (
        FoldXResult, MutationResult, FoldXError,
        StabilityLandscape, ConservationScore,
        is_foldx_available, run_foldx_stability, run_foldx_repair,
        run_foldx_mutation, empirical_stability,
        scan_mutations, find_stabilizing_mutations,
        scan_all_mutations, scan_position, compute_conservation,
        find_compensatory_mutations, rank_positions_by_mutability,
        identify_hotspot_regions,
    )
except ImportError:
    pass

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
    pass

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
    "BioCompilerError", "InvalidSequenceError",
    "CertificateGenerationError", "CertificateVerificationError",
    "UnknownPredicateError", "OptimizationError",
    "UnsupportedOrganismError", "InvalidProteinError",
    "FileFormatError", "SplicingError", "MutagenesisError",

    # ── Canonical constants ──────────────────────────────────
    "BLOSUM62", "HYDROPATHY", "HYDROPHOBIC_AAS", "STANDARD_AAS",

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

    # ── Solubility ───────────────────────────────────────────
    "compute_intrinsic_solubility", "compute_solubility",
    "compute_structural_solubility",
    "classify_solubility", "find_solubility_mutations",
    "generate_solubility_recommendations",
    "SolubilityResult",
    "CAMSOL_HYDROPATHY", "CAMSOL_CHARGE", "CAMSOL_ALPHA_HELIX",
    "CAMSOL_BETA_STRAND",
    "evaluate_soluble_expression", "evaluate_no_aggregation_prone_region",
    "evaluate_charge_composition", "evaluate_no_long_hydrophobic_stretch",
    "compute_approximate_pI", "compute_net_charge",
    "find_hydrophobic_stretches", "PKA_VALUES",

    # ── ESMFold structure prediction ─────────────────────────
    "ESMFoldResult", "ESMFoldError", "ESMFoldCache",
    "is_esmfold_available", "predict_structure", "predict_structure_batch",
    "predict_batch", "predict_proteins", "format_batch_report",
    "parse_pdb", "compute_backbone_dihedrals", "classify_plddt",
    "estimate_contact_map", "analyze_structure",
    "validate_batch_input", "estimate_batch_time",
    "BatchStructureRequest", "BatchStructureResult",

    # ── Immunogenicity & MHC binding ─────────────────────────
    "ImmunogenicityResult",
    "predict_t_cell_epitopes", "predict_b_cell_epitopes",
    "compute_surface_accessibility_approx",
    "compute_immunogenicity", "find_deimmunization_mutations",
    "MHC_I_PREFERENCES", "MHC_II_PREFERENCES",
    "ANTIGENICITY_PROPENSITY",
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
    "FoldXResult", "MutationResult", "FoldXError",
    "StabilityLandscape", "ConservationScore",
    "is_foldx_available", "run_foldx_stability", "run_foldx_repair",
    "run_foldx_mutation", "empirical_stability",
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
]
