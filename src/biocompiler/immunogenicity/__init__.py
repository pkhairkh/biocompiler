"""BioCompiler Immunogenicity Subpackage.

Provides immunogenicity scoring, MHC binding prediction, deimmunization,
and immunogenicity type-check predicates.
"""
from .core import *  # noqa: F401,F403
from biocompiler.immunogenicity.deimmunization import *  # noqa: F401,F403
from biocompiler.immunogenicity.mhcflurry_adapter import *  # noqa: F401,F403
from biocompiler.immunogenicity.mhcflurry_population import *  # noqa: F401,F403
from biocompiler.immunogenicity.netmhcpan import *  # noqa: F401,F403
from .predicates import *  # noqa: F401,F403

__all__ = [  # noqa: F405
    "ALLELE_CLASSIFICATION", "ALL_SCALES", "ANTIGENICITY_SCALE",
    "AppliedMutation", "BCellEpitopeDict", "BindingImpactEntry",
    "DeimmunizationResult", "EpitopeMutation", "EpitopePredictionResult",
    "EpitopeRegion", "ImmunogenicityError", "ImmunogenicityPrediction",
    "ImmunogenicityResult", "MHCBindingResult", "MHCPredictionResult",
    "MHCflurryClient", "MutationImpactResult", "MutationValidation",
    "NetMHCpanCache", "NetMHCpanClient", "NetMHCpanError",
    "PeptideResult", "PopulationCoverageResult", "PredicateResult",
    "ValidationReport",
    "batch_predict", "batch_predict_binding_netmhcpan",
    "binding_score_to_ic50", "classify_binding", "classify_binding_rank",
    "clear_cache", "compute_epitope_density", "compute_immunogenicity",
    "compute_immunogenicity_batch", "compute_mutation_impact",
    "compute_population_coverage", "deimmunize", "download_models",
    "evaluate_low_immunogenicity",
    "evaluate_no_dominant_b_cell_epitope",
    "evaluate_no_strong_t_cell_epitope",
    "evaluate_population_coverage_safe",
    "find_coverage_optimizing_alleles", "find_deimmunization_mutations",
    "find_epitope_disrupting_mutations", "get_allele_frequency",
    "get_population_coverage", "ic50_to_binding_score",
    "is_mhcflurry_available", "is_netmhcpan_available",
    "is_netmhcpan_installed", "parse_netmhcpan_output",
    "predict_all", "predict_b_cell_epitopes", "predict_batch",
    "predict_binding", "predict_binding_netmhcpan",
    "predict_immunogenicity", "predict_mhc_i_binding",
    "predict_mhc_ii_binding", "predict_t_cell_epitopes",
    "rank_deimmunization_mutations", "scan_peptides",
    "score_peptide_pssm", "suggest_mutations",
    "validate_deimmunized_protein",
]
