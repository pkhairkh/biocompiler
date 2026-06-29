"""
BioCompiler Application Layer — Service modules for business logic.

This package contains service functions that implement the core business
logic of the API and CLI, decoupled from HTTP/CLI concerns. Route handlers
and CLI command handlers are thin wrappers that parse input, call these
services, and format responses.
"""

# API-layer services (Wave 4a)
from .optimization_service import (
    resolve_organism_domain,
    type_check_sequence,
    optimize_protein,
    type_check_batch_item,
    optimize_batch_item,
)
from .assessment_service import (
    predict_structure,
    assess_structure_quality,
    analyze_stability,
    scan_stability_mutations,
    analyze_solubility,
    find_solubility_mutations,
    analyze_immunogenicity,
    deimmunize_protein,
    full_assessment,
    structure_batch_item,
    stability_batch_item,
    solubility_batch_item,
    immunogenicity_batch_item,
    assessment_verdict_to_verdict,
)
from .export_service import (
    export_fasta,
    export_genbank,
    export_sbol3,
    export_batch_item,
)
from .provenance_service import (
    store_provenance,
    retrieve_provenance,
    query_provenance,
    get_provenance_store,
)

# CLI-layer services (Wave 4b)
from .cli_services import (
    run_optimization,
    run_batch_optimization,
    run_check_predicates,
    resolve_organism_arg,
    resolve_source_organism_arg,
    clear_engine_caches,
    format_optimization_json,
    format_batch_json,
    OptimizationResult,
    BatchOptimizationResult,
    CheckResult,
)

__all__ = [
    # API: Optimization
    "resolve_organism_domain",
    "type_check_sequence",
    "optimize_protein",
    "type_check_batch_item",
    "optimize_batch_item",
    # API: Assessment
    "predict_structure",
    "assess_structure_quality",
    "analyze_stability",
    "scan_stability_mutations",
    "analyze_solubility",
    "find_solubility_mutations",
    "analyze_immunogenicity",
    "deimmunize_protein",
    "full_assessment",
    "structure_batch_item",
    "stability_batch_item",
    "solubility_batch_item",
    "immunogenicity_batch_item",
    "assessment_verdict_to_verdict",
    # API: Export
    "export_fasta",
    "export_genbank",
    "export_sbol3",
    "export_batch_item",
    # API: Provenance
    "store_provenance",
    "retrieve_provenance",
    "query_provenance",
    "get_provenance_store",
    # CLI: Services
    "run_optimization",
    "run_batch_optimization",
    "run_check_predicates",
    "resolve_organism_arg",
    "resolve_source_organism_arg",
    "clear_engine_caches",
    "format_optimization_json",
    "format_batch_json",
    "OptimizationResult",
    "BatchOptimizationResult",
    "CheckResult",
]
