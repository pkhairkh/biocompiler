"""
BioCompiler API — FastAPI Web Service (refactored package).

This package provides the BioCompiler REST API, split into:
- :mod:`biocompiler.api.models` — Pydantic request/response models and constants
- :mod:`biocompiler.api.auth` — API key authentication and rate limiting
- :mod:`biocompiler.api.routes` — thin HTTP route handlers
- :mod:`biocompiler.api.app` — FastAPI app factory (``create_app``)

Business logic lives in :mod:`biocompiler.application`:
- :mod:`biocompiler.application.optimization_service` — optimization orchestration
- :mod:`biocompiler.application.assessment_service` — protein structure/stability/solubility/immunogenicity
- :mod:`biocompiler.application.export_service` — FASTA, GenBank, SBOL3 export
- :mod:`biocompiler.application.provenance_service` — provenance storage and retrieval

Backward compatibility: ``from biocompiler.api import create_app`` still works.
"""

from .app import create_app
from .models import (
    # Constants
    MAX_PROTEIN_LENGTH,
    MAX_PROTEIN_SEQUENCE_LENGTH,
    MAX_BATCH_SIZE,
    MAX_REQUEST_SIZE,
    OPTIMIZE_TIMEOUT_S,
    MAX_DNA_LENGTH,
    MAX_DNA_SEQUENCE_LENGTH,
    BATCH_CHECK_MAX,
    BATCH_OPTIMIZE_MAX,
    BATCH_EXPORT_MAX,
    BATCH_ITEM_TIMEOUT_S,
    ESMFOLD_TIMEOUT_S,
    # Input models
    SequenceInput,
    ProteinInput,
    CertificateInput,
    ExportFastaInput,
    ExportGenbankInput,
    ScanInput,
    ExportSbol3Input,
    StructurePredictInput,
    QualityAssessInput,
    StabilityInput,
    MutationScanInput,
    SolubilityInput,
    ImmunogenicityInput,
    DeimmunizeInput,
    FullAssessmentInput,
    # Response models
    TypeCheckResponse,
    OptimizeResponse,
    VerifyResponse,
    ScanResponse,
    OrganismResponse,
    PredicateResponse,
    HealthResponse,
    InfoResponse,
    EnzymeListResponse,
    ProvenanceDetailResponse,
    ProvenanceListResponse,
    ProvenanceRecordSummary,
    ExportFastaResponse,
    ExportGenbankResponse,
    ExportSbol3Response,
    DatasetValidationResponse,
    DatasetValidationResult,
    ProvenanceResponse,
    ProvenanceExplainResponse,
    ProvenanceReportResponse,
    StructurePredictResponse,
    QualityAssessResponse,
    StabilityResponse,
    MutationScanResponse,
    SolubilityResponse,
    ImmunogenicityResponse,
    DeimmunizeResponse,
    FullAssessmentResponse,
    # Batch models
    BatchCheckItem,
    BatchCheckInput,
    BatchCheckSummary,
    BatchCheckResponse,
    BatchOptimizeItem,
    BatchOptimizeInput,
    BatchOptimizeSummary,
    BatchOptimizeResponse,
    BatchExportItem,
    BatchExportInput,
    BatchExportResultItem,
    BatchExportResponse,
    BatchStructureItem,
    BatchStructureInput,
    BatchStructureResultItem,
    BatchStructureResponse,
    BatchStabilityItem,
    BatchStabilityInput,
    BatchStabilityResponse,
    BatchSolubilityItem,
    BatchSolubilityInput,
    BatchSolubilityResponse,
    BatchImmunogenicityItem,
    BatchImmunogenicityInput,
    BatchImmunogenicityResponse,
    # Domain resolution
    resolve_organism_domain,
    # Validation helpers
    validate_protein_input,
    validate_organism_input,
)
from .auth import (
    verify_api_key,
    API_KEY_NAME,
    RATE_LIMIT_RPM,
    set_no_auth_flag,
    get_auth_mode,
    get_configured_api_keys,
    is_auth_enabled,
    _rate_limiter,
    _check_rate_limit,
    _check_batch_rate_limit,
    _generate_and_persist_api_key,
    _get_api_key_file_path,
    _AUTH_MODE,
    _CONFIGURED_API_KEYS,
    _api_key_header,
    _NO_AUTH_CLI_FLAG,
)
from biocompiler.optimizer import batch_optimize

# Create the default app instance
app = create_app()


# ── Lazy attributes ────────────────────────────────────────────────────
# Deferred to first access to keep `import biocompiler` side-effect-free.
# Importing them eagerly would transitively load routes.py and
# provenance_service.py, both of which create files/directories in
# ~/.biocompiler/ at module-load time.

_LAZY_NAMES = {"create_app", "app", "_provenance_store"}


def __getattr__(name: str):
    """Lazily resolve deferred attributes to avoid import-time side effects."""
    if name == "create_app":
        from .app import create_app as _create_app
        globals()["create_app"] = _create_app
        return _create_app
    if name == "app":
        from .app import create_app as _create_app
        _app = _create_app()
        globals()["app"] = _app
        return _app
    if name == "_provenance_store":
        from ..application.provenance_service import _provenance_store as _store
        globals()["_provenance_store"] = _store
        return _store
    raise AttributeError(f"module 'biocompiler.api' has no attribute {name!r}")


def __dir__() -> list[str]:
    """Include lazy attributes in dir(biocompiler.api)."""
    return sorted(list(globals().keys()) + list(_LAZY_NAMES))


__all__ = [
    "create_app",
    "app",
    "verify_api_key",
    "validate_protein_input",
    "validate_organism_input",
    "API_KEY_NAME",
    "RATE_LIMIT_RPM",
    "BATCH_CHECK_MAX",
    "BATCH_OPTIMIZE_MAX",
    "BATCH_EXPORT_MAX",
    "BATCH_ITEM_TIMEOUT_S",
    "ESMFOLD_TIMEOUT_S",
    "MAX_PROTEIN_LENGTH",
    "MAX_PROTEIN_SEQUENCE_LENGTH",
    "MAX_BATCH_SIZE",
    "MAX_REQUEST_SIZE",
    "MAX_DNA_LENGTH",
    "MAX_DNA_SEQUENCE_LENGTH",
    "OPTIMIZE_TIMEOUT_S",
    "batch_optimize",
    "SequenceInput",
    "ProteinInput",
    "CertificateInput",
    "ExportFastaInput",
    "ExportGenbankInput",
    "ScanInput",
    "TypeCheckResponse",
    "OptimizeResponse",
    "VerifyResponse",
    "ScanResponse",
    "OrganismResponse",
    "PredicateResponse",
    "HealthResponse",
    "InfoResponse",
    "EnzymeListResponse",
    "ProvenanceDetailResponse",
    "ProvenanceListResponse",
    "ProvenanceRecordSummary",
    "ExportFastaResponse",
    "ExportGenbankResponse",
    "ExportSbol3Input",
    "ExportSbol3Response",
    "DatasetValidationResponse",
    "DatasetValidationResult",
    "StructurePredictInput",
    "QualityAssessInput",
    "StabilityInput",
    "MutationScanInput",
    "SolubilityInput",
    "ImmunogenicityInput",
    "DeimmunizeInput",
    "FullAssessmentInput",
    "StructurePredictResponse",
    "QualityAssessResponse",
    "StabilityResponse",
    "MutationScanResponse",
    "SolubilityResponse",
    "ImmunogenicityResponse",
    "DeimmunizeResponse",
    "FullAssessmentResponse",
    # Batch models
    "BatchCheckItem",
    "BatchCheckInput",
    "BatchCheckSummary",
    "BatchCheckResponse",
    "BatchOptimizeItem",
    "BatchOptimizeInput",
    "BatchOptimizeSummary",
    "BatchOptimizeResponse",
    "BatchExportItem",
    "BatchExportInput",
    "BatchExportResultItem",
    "BatchExportResponse",
    "BatchStructureItem",
    "BatchStructureInput",
    "BatchStructureResultItem",
    "BatchStructureResponse",
    "BatchStabilityItem",
    "BatchStabilityInput",
    "BatchStabilityResponse",
    "BatchSolubilityItem",
    "BatchSolubilityInput",
    "BatchSolubilityResponse",
    "BatchImmunogenicityItem",
    "BatchImmunogenicityInput",
    "BatchImmunogenicityResponse",
    "resolve_organism_domain",
    "set_no_auth_flag",
    "get_auth_mode",
    "get_configured_api_keys",
    "is_auth_enabled",
    "ProvenanceResponse",
    "ProvenanceExplainResponse",
    "ProvenanceReportResponse",
    "ProvenanceStore",
    "_provenance_store",
    "_rate_limiter",
    "_check_rate_limit",
    "_check_batch_rate_limit",
]

# Re-export ProvenanceStore from decision_provenance for backward compat
from biocompiler.provenance.decision_provenance import ProvenanceStore
