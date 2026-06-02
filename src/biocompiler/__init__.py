"""
BioCompiler — Machine-Verified Gene Design

A compiler framework for human protein synthesis using intermediate
representations. Pipeline:

  Scanner → NDFST Splicing → Translation → Type Check → Certificate → Verify

All computation is DETERMINISTIC: same input always produces identical output.
"""

__version__ = "7.3.0"

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

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
from .scanner import (
    validate_dna_sequence,
    gc_content,
    scan_sequence,
)
from .translation import (
    translate,
    compute_cai,
    find_orfs,
)
from .splicing import compute_splice_isoforms

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

try:
    from .certificate import generate_certificate, verify_certificate
except ImportError:
    pass

try:
    from .maxentscan import (
        score_donor, score_acceptor, scan_splice_sites,
        max_donor_score, max_acceptor_score,
    )
except ImportError:
    pass

from .optimization import BioOptimizer, optimize_sequence, OptimizationResult

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

try:
    from .benchmark import run_benchmarks, BenchmarkReport
except ImportError:
    try:
        from .benchmark import run_benchmark as run_benchmarks
        BenchmarkReport = None
    except ImportError:
        run_benchmarks = None
        BenchmarkReport = None

try:
    from .organism_db import OrganismDatabase, get_database
except ImportError:
    pass

try:
    from .tissue_data import get_tissue_weights, list_available_tissues, add_custom_tissue
except ImportError:
    pass

try:
    from .dna_chisel_compat import compare_optimizers, run_comparative_benchmark
except ImportError:
    compare_optimizers = None
    run_comparative_benchmark = None

try:
    from .tool_comparison import (
        run_head_to_head_benchmark, compare_tools, HeadToHeadReport,
        format_head_to_head_text, format_head_to_head_json,
        is_dna_chisel_available,
    )
except ImportError:
    pass

try:
    from .dataset_validation import run_dataset_validation, DatasetValidationReport
except ImportError:
    pass

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

try:
    from .mutagenesis import (
        MutagenesisResult, AASubstitution, BLOSUM62,
        GT_MANDATORY_AAS, AG_MANDATORY_AAS,
        is_gt_mandatory, is_ag_mandatory, diagnose_optimizer_weakness,
        force_gt_free_reoptimization,
        type_directed_mutagenesis, find_unrepairable_cryptic_donors,
        find_unrepairable_cryptic_acceptors, propose_substitutions,
        apply_substitution,
    )
except ImportError:
    pass

try:
    from .camsol import (
        compute_intrinsic_solubility, compute_solubility,
        compute_structural_solubility,
        classify_solubility, find_solubility_mutations,
        generate_solubility_recommendations,
        SolubilityResult,
        CAMSOL_HYDROPATHY, CAMSOL_CHARGE, CAMSOL_ALPHA_HELIX,
        CAMSOL_BETA_STRAND, BLOSUM62,
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

try:
    from .esmfold import (
        ESMFoldResult, ESMFoldError,
        is_esmfold_available, predict_structure, predict_structure_batch,
        parse_pdb, compute_backbone_dihedrals, classify_plddt,
        estimate_contact_map, analyze_structure,
    )
except ImportError:
    pass

try:
    from .immunogenicity import (
        ImmunogenicityResult,
        predict_t_cell_epitopes,
        predict_b_cell_epitopes,
        compute_surface_accessibility_approx,
        compute_immunogenicity,
        find_deimmunization_mutations,
        MHC_I_PREFERENCES,
        MHC_II_PREFERENCES,
        ANTIGENICITY_PROPENSITY,
        BLOSUM62 as ImmunogenicityBLOSUM62,
    )
except ImportError:
    pass

try:
    from .deimmunization import (
        DeimmunizationResult,
        EpitopeMutation,
        deimmunize,
        find_epitope_disrupting_mutations,
        rank_deimmunization_mutations,
        validate_deimmunized_protein,
        compute_mutation_impact,
        BLOSUM62 as DeimmunizationBLOSUM62,
        HYDROPATHY,
    )
except ImportError:
    pass

try:
    from .structure_quality import (
        StructureQualityReport,
        assess_plddt,
        assess_ramachandran,
        compute_clash_score,
        compute_packing_density,
        compute_exposed_hydrophobic,
        compute_structure_quality,
        find_low_confidence_regions,
        compute_sasa_approximation,
        KYTE_DOOLITTLE,
        VDW_RADII,
    )
except ImportError:
    pass

try:
    from .structure_report import (
        ProteinAssessmentReport,
        assess_protein,
        format_assessment_text,
        format_assessment_json,
        format_assessment_html,
        generate_recommendations,
        compute_overall_verdict,
        plot_plddt_bar_svg,
        plot_solubility_profile_svg,
    )
except ImportError:
    pass

try:
    from .foldx import (
        FoldXResult, MutationResult, FoldXError,
        is_foldx_available, run_foldx_stability, run_foldx_repair,
        run_foldx_mutation, empirical_stability,
        scan_mutations, find_stabilizing_mutations,
        BLOSUM62 as FoldXBLOSUM62,
        HYDROPATHY as FoldXHYDROPATHY,
    )
except ImportError:
    pass

try:
    from .protein_design import (
        DesignResult, DesignConstraints,
        BLOSUM62 as ProteinDesignBLOSUM62,
        HYDROPATHY as ProteinDesignHYDROPATHY,
        design_thermostable, design_soluble,
        design_low_immunogenicity, design_multi_objective,
        score_mutation, find_disulfide_opportunities,
        find_proline_substitution_sites,
    )
except ImportError:
    pass

try:
    from .mhc_binding import (
        MHCBindingResult, MHCPredictionResult,
        predict_mhc_i_binding, predict_mhc_ii_binding,
        score_peptide_pssm, binding_score_to_ic50, classify_binding,
        predict_all as predict_mhc_binding,
        MHC_I_PSSM, MHC_II_PSSM, POPULATION_COVERAGE,
        DEFAULT_MHC_I_ALLELES, DEFAULT_MHC_II_ALLELES,
    )
except ImportError:
    pass

try:
    from .epitope import (
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

__all__ = [
    "__version__",
    "Verdict", "Token", "PositionRange", "SpliceIsoform",
    "TypeCheckResult", "Certificate",
    "three_valued_and", "three_valued_or", "combined_verdict",
    "BioCompilerError", "InvalidSequenceError",
    "CertificateGenerationError", "CertificateVerificationError",
    "UnknownPredicateError", "OptimizationError",
    "UnsupportedOrganismError", "InvalidProteinError",
    "FileFormatError", "SplicingError", "MutagenesisError",
    "validate_dna_sequence", "gc_content", "scan_sequence",
    "translate", "compute_cai", "find_orfs",
    "compute_splice_isoforms",
    "generate_certificate", "verify_certificate",
    "score_donor", "score_acceptor", "scan_splice_sites",
    "max_donor_score", "max_acceptor_score",
    "BioOptimizer",
    "load_grammar", "grammar_to_predicate_params",
    "load_builtin_grammar", "list_builtin_grammars",
    "export_fasta", "export_genbank", "export_genbank_with_certificate", "export_multi_fasta",
    "generate_report",
    "run_benchmarks", "BenchmarkReport",
    "OrganismDatabase", "get_database",
    "get_tissue_weights", "list_available_tissues", "add_custom_tissue",
    "compare_optimizers", "run_comparative_benchmark",
    "run_head_to_head_benchmark", "compare_tools", "HeadToHeadReport",
    "format_head_to_head_text", "format_head_to_head_json",
    "is_dna_chisel_available",
    "run_dataset_validation", "DatasetValidationReport",
    "import_fasta", "import_genbank", "import_sequence",
    "to_seqrecord", "from_seqrecord", "optimize_to_seqrecord",
    "display_sequence", "display_optimization_result", "display_type_check",
    "plot_gc_content", "plot_codon_usage",
    "MutagenesisResult", "AASubstitution", "BLOSUM62",
    "GT_MANDATORY_AAS", "AG_MANDATORY_AAS",
    "is_gt_mandatory", "is_ag_mandatory", "diagnose_optimizer_weakness",
    "force_gt_free_reoptimization",
    "type_directed_mutagenesis", "find_unrepairable_cryptic_donors",
    "find_unrepairable_cryptic_acceptors", "propose_substitutions",
    "apply_substitution",
    "ESMFoldResult", "ESMFoldError",
    "is_esmfold_available", "predict_structure", "predict_structure_batch",
    "parse_pdb", "compute_backbone_dihedrals", "classify_plddt",
    "estimate_contact_map", "analyze_structure",
    "ImmunogenicityResult",
    "predict_t_cell_epitopes",
    "predict_b_cell_epitopes",
    "compute_surface_accessibility_approx",
    "compute_immunogenicity",
    "find_deimmunization_mutations",
    "MHC_I_PREFERENCES",
    "MHC_II_PREFERENCES",
    "ANTIGENICITY_PROPENSITY",
    "ImmunogenicityBLOSUM62",
    "compute_intrinsic_solubility", "compute_solubility",
    "compute_structural_solubility",
    "classify_solubility", "find_solubility_mutations",
    "generate_solubility_recommendations",
    "SolubilityResult",
    "CAMSOL_HYDROPATHY", "CAMSOL_CHARGE", "CAMSOL_ALPHA_HELIX",
    "CAMSOL_BETA_STRAND", "BLOSUM62",
    "evaluate_soluble_expression", "evaluate_no_aggregation_prone_region",
    "evaluate_charge_composition", "evaluate_no_long_hydrophobic_stretch",
    "compute_approximate_pI", "compute_net_charge",
    "find_hydrophobic_stretches", "PKA_VALUES",
    "DeimmunizationResult",
    "EpitopeMutation",
    "deimmunize",
    "find_epitope_disrupting_mutations",
    "rank_deimmunization_mutations",
    "validate_deimmunized_protein",
    "compute_mutation_impact",
    "DeimmunizationBLOSUM62",
    "HYDROPATHY",
    "StructureQualityReport",
    "assess_plddt",
    "assess_ramachandran",
    "compute_clash_score",
    "compute_packing_density",
    "compute_exposed_hydrophobic",
    "compute_structure_quality",
    "find_low_confidence_regions",
    "compute_sasa_approximation",
    "KYTE_DOOLITTLE",
    "VDW_RADII",
    "ProteinAssessmentReport",
    "assess_protein",
    "format_assessment_text",
    "format_assessment_json",
    "format_assessment_html",
    "generate_recommendations",
    "compute_overall_verdict",
    "plot_plddt_bar_svg",
    "plot_solubility_profile_svg",
    "FoldXResult", "MutationResult", "FoldXError",
    "is_foldx_available", "run_foldx_stability", "run_foldx_repair",
    "run_foldx_mutation", "empirical_stability",
    "scan_mutations", "find_stabilizing_mutations",
    "FoldXBLOSUM62", "FoldXHYDROPATHY",
    "DesignResult", "DesignConstraints",
    "ProteinDesignBLOSUM62", "ProteinDesignHYDROPATHY",
    "design_thermostable", "design_soluble",
    "design_low_immunogenicity", "design_multi_objective",
    "score_mutation", "find_disulfide_opportunities",
    "find_proline_substitution_sites",
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
]
