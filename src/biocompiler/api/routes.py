"""
BioCompiler API — Thin HTTP route handlers.

Route handlers parse input, call application service functions,
and return results. All business logic lives in
:mod:`biocompiler.application`.
"""

import asyncio
import json
import logging
from datetime import datetime, timezone
from typing import Any, Optional

from fastapi import APIRouter, Depends, HTTPException

from .auth import verify_api_key, _check_rate_limit, _check_batch_rate_limit
from .models import (
    # Constants
    BATCH_CHECK_MAX,
    BATCH_OPTIMIZE_MAX,
    BATCH_EXPORT_MAX,
    BATCH_ITEM_TIMEOUT_S,
    ESMFOLD_TIMEOUT_S,
    OPTIMIZE_TIMEOUT_S,
    MAX_REQUEST_SIZE,
    BATCH_PROTEIN_MAX,
    # Input models
    SequenceInput,
    ProteinInput,
    CertificateInput,
    ExportFastaInput,
    ExportGenbankInput,
    ScanInput,
    ExportSbol3Input,
    BatchCheckInput,
    BatchOptimizeInput,
    FastBatchOptimizeInput,
    BatchExportInput,
    BatchStructureInput,
    BatchStabilityInput,
    BatchSolubilityInput,
    BatchImmunogenicityInput,
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
    InfoResponse,
    PredicateResponse,
    HealthResponse,
    BatchCheckResponse,
    BatchCheckSummary,
    BatchOptimizeResponse,
    BatchOptimizeSummary,
    BatchExportResponse,
    BatchExportResultItem,
    BatchStructureResponse,
    BatchStructureResultItem,
    BatchStabilityResponse,
    BatchSolubilityResponse,
    BatchImmunogenicityResponse,
    StructurePredictResponse,
    QualityAssessResponse,
    StabilityResponse,
    MutationScanResponse,
    SolubilityResponse,
    ImmunogenicityResponse,
    DeimmunizeResponse,
    FullAssessmentResponse,
)
from biocompiler.sequence.scanner import scan_sequence, validate_dna_sequence, gc_content
from biocompiler.expression.translation import find_orfs
from biocompiler.sequence.splicing import compute_splice_isoforms
from biocompiler.provenance.certificate import verify_certificate
from biocompiler.optimizer import batch_optimize
from biocompiler.export.core import export_fasta as _export_fasta, export_genbank as _export_genbank
from biocompiler.shared.types import Certificate
from biocompiler.shared.constants import RESTRICTION_ENZYMES
from ..organisms import SUPPORTED_ORGANISMS, CODON_USAGE_TABLES
from biocompiler.shared.exceptions import (
    BioCompilerError, BiosecurityError, InvalidSequenceError,
    CertificateGenerationError, CertificateVerificationError,
    UnsupportedOrganismError, InvalidProteinError,
    OptimizationConstraintError,
)
from ..application.optimization_service import (
    type_check_sequence,
    optimize_protein,
    type_check_batch_item,
    optimize_batch_item,
)
from ..application.assessment_service import (
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
from ..application.export_service import (
    export_fasta as svc_export_fasta,
    export_genbank as svc_export_genbank,
    export_sbol3,
    export_batch_item,
)
from ..application.provenance_service import (
    store_provenance,
    retrieve_provenance,
    query_provenance,
)

logger = logging.getLogger(__name__)

# ─── Routers ──────────────────────────────────────────────────────

_main_router = APIRouter()
_protein_router = APIRouter(tags=["Protein Analysis"])


# ═══════════════════════════════════════════════════════════════════
# Health / Info / Organisms / Predicates / Enzymes
# ═══════════════════════════════════════════════════════════════════


@_main_router.get("/health", response_model=HealthResponse)
async def health_check() -> HealthResponse:
    """Health check endpoint (no auth required)."""
    from .auth import _is_auth_enabled, RATE_LIMIT_RPM
    from .app import get_cors_config
    import biocompiler

    cors_origins, _allow_creds = get_cors_config()

    return HealthResponse(
        status="healthy",
        version=biocompiler.__version__,
        timestamp=datetime.now(timezone.utc).isoformat(),
        auth_enabled=_is_auth_enabled(),
        rate_limit_rpm=RATE_LIMIT_RPM,
        cors_origins=cors_origins,
        cors_allow_credentials=_allow_creds,
    )


@_main_router.get("/info", response_model=InfoResponse)
async def api_info() -> InfoResponse:
    """Return API configuration and limits."""
    from .models import (
        MAX_PROTEIN_LENGTH, MAX_DNA_LENGTH, MAX_BATCH_SIZE,
        MAX_REQUEST_SIZE, OPTIMIZE_TIMEOUT_S,
    )
    import biocompiler

    # Try to get safety version from biosecurity module
    try:
        from ..biosecurity import __version__ as safety_ver
    except (ImportError, AttributeError):
        safety_ver = biocompiler.__version__

    # Get the unique canonical organism names (filter out aliases)
    canonical_organisms = sorted(set(
        name for name in SUPPORTED_ORGANISMS
        if not name.islower()
    ))

    return InfoResponse(
        max_protein_length=MAX_PROTEIN_LENGTH,
        max_dna_length=MAX_DNA_LENGTH,
        max_batch_size=MAX_BATCH_SIZE,
        max_request_size=MAX_REQUEST_SIZE,
        optimize_timeout_s=OPTIMIZE_TIMEOUT_S,
        supported_organisms=canonical_organisms,
        api_version=biocompiler.__version__,
        safety_version=safety_ver,
    )


@_main_router.get("/organisms", response_model=OrganismResponse)
async def list_organisms() -> OrganismResponse:
    """List all supported organisms with codon usage data."""
    organisms = []
    for name in SUPPORTED_ORGANISMS:
        usage = CODON_USAGE_TABLES.get(name, {})
        organisms.append({
            "name": name,
            "codon_count": len(usage),
            "available": bool(usage),
        })
    return OrganismResponse(organisms=organisms)


@_main_router.get("/predicates", response_model=PredicateResponse)
async def list_predicates() -> PredicateResponse:
    """List all registered type predicates."""
    from ..type_system import registry
    return PredicateResponse(predicates=registry.names())


@_main_router.get("/enzymes")
async def list_enzymes() -> dict[str, dict[str, str]]:
    """List all known restriction enzymes."""
    return {
        "enzymes": {
            name: site for name, site in RESTRICTION_ENZYMES.items()
        }
    }


# ═══════════════════════════════════════════════════════════════════
# Provenance
# ═══════════════════════════════════════════════════════════════════


@_main_router.get("/provenance/{record_id}")
async def get_provenance(record_id: str, client_id: str = Depends(verify_api_key)) -> dict:
    """Retrieve a stored provenance record by ID."""
    try:
        return retrieve_provenance(record_id)
    except (FileNotFoundError, ValueError):
        raise HTTPException(
            status_code=404,
            detail=f"Provenance record not found: {record_id}",
        )


@_main_router.get("/provenance")
async def query_provenance_endpoint(
    protein_name: Optional[str] = None,
    organism: Optional[str] = None,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
    client_id: str = Depends(verify_api_key),
) -> dict:
    """Query/list stored provenance records."""
    date_range = None
    if date_from and date_to:
        date_range = (date_from, date_to)
    return query_provenance(
        protein_name=protein_name,
        organism=organism,
        date_range=date_range,
    )


# ═══════════════════════════════════════════════════════════════════
# Core Endpoints: Check, Optimize, Verify, Scan
# ═══════════════════════════════════════════════════════════════════


@_main_router.post("/check", response_model=TypeCheckResponse)
async def check_sequence(input_data: SequenceInput, client_id: str = Depends(verify_api_key)) -> TypeCheckResponse:
    """Type-check a DNA sequence against all registered predicates."""
    try:
        result = type_check_sequence(
            seq=input_data.sequence,
            organism=input_data.organism,
            exon_boundaries=input_data.exon_boundaries,
            gc_lo=input_data.gc_lo,
            gc_hi=input_data.gc_hi,
            cai_threshold=input_data.cai_threshold,
            enzymes=input_data.enzymes,
            cellular_context=input_data.cellular_context,
        )
        return TypeCheckResponse(**result)
    except (InvalidSequenceError, UnsupportedOrganismError) as e:
        raise HTTPException(status_code=400, detail=str(e))
    except BioCompilerError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        logger.exception("Unexpected error during type-check")
        raise HTTPException(status_code=500, detail=f"Type-check failed: {e}")


@_main_router.post("/optimize", response_model=OptimizeResponse)
async def optimize(input_data: ProteinInput, client_id: str = Depends(verify_api_key)) -> OptimizeResponse:
    """Optimize a DNA sequence for a target protein."""
    try:
        result = await asyncio.wait_for(
            asyncio.to_thread(
                optimize_protein,
                protein=input_data.protein,
                organism=input_data.organism,
                gc_lo=input_data.gc_lo,
                gc_hi=input_data.gc_hi,
                cai_threshold=input_data.cai_threshold,
                restriction_sites=input_data.enzymes,
                cryptic_splice_threshold=input_data.cryptic_splice_threshold,
                organism_domain=input_data.organism_domain,
                source_organism=input_data.source_organism,
                therapeutic=input_data.therapeutic,
                self_protein=input_data.self_protein,
                track_provenance=input_data.track_provenance,
                strict_mode=input_data.strict_mode,
            ),
            timeout=OPTIMIZE_TIMEOUT_S,
        )

        # Persist provenance if tracking enabled
        provenance_id = None
        if input_data.track_provenance and result.get("decision_trail") is not None:
            try:
                provenance_id = store_provenance(result["decision_trail"])
            except Exception as exc:
                logger.warning("Failed to persist provenance trail: %s", exc)

        return OptimizeResponse(
            sequence=result["sequence"],
            protein=result["protein"],
            cai=result["cai"],
            gc_content=result["gc_content"],
            satisfied_predicates=result["satisfied_predicates"],
            failed_predicates=result["failed_predicates"],
            fallback_used=result["fallback_used"],
            provenance_id=provenance_id,
            organism_domain=result["organism_domain"],
            source_organism=result["source_organism"],
            therapeutic=result["therapeutic"],
            self_protein=result["self_protein"],
        )
    except (InvalidProteinError, UnsupportedOrganismError) as e:
        raise HTTPException(status_code=400, detail=str(e))
    except OptimizationConstraintError as e:
        raise HTTPException(
            status_code=422,
            detail={
                "error": "OptimizationConstraintError",
                "message": str(e),
                "failed_predicates": e.failed_predicates,
                "suggestion": "Set strict_mode=False to allow partial results, "
                              "or relax constraints (e.g., widen GC range, "
                              "remove restriction enzymes).",
                "partial_result": {
                    "sequence": e.partial_result.sequence if e.partial_result else None,
                    "protein": e.partial_result.protein if e.partial_result else None,
                    "cai": e.partial_result.cai if e.partial_result else None,
                    "gc_content": e.partial_result.gc_content if e.partial_result else None,
                    "failed_predicates": e.partial_result.failed_predicates if e.partial_result else None,
                    "satisfied_predicates": e.partial_result.satisfied_predicates if e.partial_result else None,
                } if e.partial_result else None,
            },
        )
    except BiosecurityError as e:
        detail: dict[str, Any] = {
            "error": "BiosecurityError",
            "message": str(e),
            "risk_level": e.risk_level,
            "flagged_categories": e.flagged_categories,
            "matches": [
                {
                    "category": m.category,
                    "name": m.name,
                    "position": m.position,
                    "matched_sequence": m.matched_sequence,
                    "confidence": m.confidence,
                }
                for m in (e.matches or [])
            ],
            "recommendations": (
                e.report.recommendations
                if e.report is not None and hasattr(e.report, "recommendations")
                else []
            ),
        }
        raise HTTPException(status_code=403, detail=detail)
    except BioCompilerError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except asyncio.TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=(
                f"Optimization timed out after {OPTIMIZE_TIMEOUT_S} seconds. "
                f"Try a shorter protein or relax constraints."
            ),
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("Unexpected error during optimization")
        raise HTTPException(status_code=500, detail=f"Optimization failed: {e}")


@_main_router.post("/verify", response_model=VerifyResponse)
async def verify(input_data: CertificateInput, client_id: str = Depends(verify_api_key)) -> VerifyResponse:
    """Independently verify a guarantee certificate."""
    try:
        status, failures = verify_certificate(input_data.certificate)
        return VerifyResponse(status=status, failure_reasons=failures)
    except (CertificateVerificationError, ValueError, KeyError) as e:
        raise HTTPException(status_code=400, detail=f"Verification error: {e}")
    except Exception as e:
        logger.exception("Unexpected error during certificate verification")
        raise HTTPException(status_code=500, detail="Internal verification error")


@_main_router.post("/scan", response_model=ScanResponse)
async def scan(input_data: ScanInput, client_id: str = Depends(verify_api_key)) -> ScanResponse:
    """Scan a DNA sequence for biological motifs."""
    try:
        seq = input_data.sequence.upper()
        seq = validate_dna_sequence(seq)
        tokens = scan_sequence(seq, input_data.enzymes)

        token_dicts = [
            {
                "position": t.position,
                "element_type": t.element_type,
                "match_sequence": t.match_sequence,
                "score": t.score,
                "frame": t.frame,
                "strand": t.strand,
            }
            for t in tokens
        ]

        orf_dicts = None
        if input_data.find_orfs:
            orf_dicts = find_orfs(seq)

        return ScanResponse(
            sequence_length=len(seq),
            tokens=token_dicts,
            orfs=orf_dicts,
        )
    except InvalidSequenceError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Unexpected error during sequence scan")
        raise HTTPException(status_code=500, detail=f"Scan failed: {e}")


# ═══════════════════════════════════════════════════════════════════
# Export Endpoints
# ═══════════════════════════════════════════════════════════════════


@_main_router.post("/export/fasta")
async def export_fasta_endpoint(input_data: ExportFastaInput, client_id: str = Depends(verify_api_key)) -> dict[str, str]:
    """Export a sequence in FASTA format."""
    try:
        fasta = svc_export_fasta(
            sequence=input_data.sequence,
            identifier=input_data.identifier,
            description=input_data.description,
            organism=input_data.organism,
        )
        return {"format": "fasta", "content": fasta}
    except (InvalidSequenceError, ValueError) as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        logger.exception("FASTA export failed unexpectedly")
        raise HTTPException(status_code=500, detail=f"FASTA export failed: {e}")


@_main_router.post("/export/genbank")
async def export_genbank_endpoint(input_data: ExportGenbankInput, client_id: str = Depends(verify_api_key)) -> dict[str, str]:
    """Export a sequence in GenBank format with optional certificate embedding."""
    try:
        cert = None
        if input_data.certificate:
            cert = Certificate.from_dict(input_data.certificate)

        genbank = svc_export_genbank(
            sequence=input_data.sequence,
            locus_name=input_data.locus_name,
            definition=input_data.definition,
            organism=input_data.organism,
            gene_name=input_data.gene_name,
            exon_boundaries=input_data.exon_boundaries,
            certificate=cert,
        )
        return {"format": "genbank", "content": genbank}
    except (InvalidSequenceError, ValueError, KeyError) as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        logger.exception("GenBank export failed unexpectedly")
        raise HTTPException(status_code=500, detail=f"GenBank export failed: {e}")


@_main_router.post("/export/sbol3")
async def export_sbol3_endpoint(input_data: ExportSbol3Input, client_id: str = Depends(verify_api_key)) -> dict:
    """Export a sequence in SBOL3 format."""
    try:
        content = export_sbol3(
            sequence=input_data.sequence,
            organism=input_data.organism,
            gene_name=input_data.gene_name,
            fmt=input_data.format,
        )
        return {"format": input_data.format, "content": content}
    except Exception as e:
        logger.exception("SBOL3 export failed unexpectedly")
        raise HTTPException(status_code=500, detail=f"SBOL3 export failed: {e}")


# ═══════════════════════════════════════════════════════════════════
# Batch Endpoints
# ═══════════════════════════════════════════════════════════════════


@_main_router.post("/batch/check", response_model=BatchCheckResponse)
async def batch_check(
    input_data: BatchCheckInput,
    client_id: str = Depends(verify_api_key),
) -> BatchCheckResponse:
    """Type-check multiple DNA sequences in a single request."""
    n = len(input_data.sequences)
    if n > BATCH_CHECK_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_CHECK_MAX} sequences",
        )

    _check_batch_rate_limit(client_id, n)

    results: list[dict[str, Any]] = []
    pass_count = 0
    fail_count = 0
    uncertain_count = 0
    error_count = 0

    for idx, item in enumerate(input_data.sequences):
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(type_check_batch_item,
                    sequence=item.sequence,
                    organism=item.organism,
                    exon_boundaries=item.exon_boundaries,
                    gc_lo=item.gc_lo,
                    gc_hi=item.gc_hi,
                    cai_threshold=item.cai_threshold,
                    enzymes=item.enzymes,
                    cellular_context=item.cellular_context,
                ),
                timeout=BATCH_ITEM_TIMEOUT_S,
            )
            results.append({"index": idx, "status": "success", **result})
            verdict = result.get("overall_verdict", "ERROR")
            if verdict == "PASS":
                pass_count += 1
            elif verdict == "FAIL":
                fail_count += 1
            elif verdict == "UNCERTAIN":
                uncertain_count += 1
            else:
                error_count += 1
        except asyncio.TimeoutError:
            logger.warning("Batch check item %d timed out", idx)
            results.append({
                "index": idx, "status": "error",
                "error": f"Item timed out after {BATCH_ITEM_TIMEOUT_S}s",
            })
            error_count += 1
        except Exception as e:
            logger.warning("Batch check item %d failed: %s", idx, e)
            results.append({
                "index": idx, "status": "error", "error": str(e),
            })
            error_count += 1

    return BatchCheckResponse(
        results=results,
        summary=BatchCheckSummary(
            total=n,
            **{"pass": pass_count},
            fail=fail_count,
            uncertain=uncertain_count,
            errors=error_count,
        ),
    )


@_main_router.post("/batch/optimize", response_model=BatchOptimizeResponse)
async def batch_optimize_endpoint(
    input_data: BatchOptimizeInput,
    client_id: str = Depends(verify_api_key),
) -> BatchOptimizeResponse:
    """Optimize multiple proteins in a single request."""
    n = len(input_data.proteins)
    if n > BATCH_OPTIMIZE_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_OPTIMIZE_MAX} proteins",
        )

    _check_batch_rate_limit(client_id, n)

    optimize_timeout = max(BATCH_ITEM_TIMEOUT_S, 60)

    results: list[dict[str, Any]] = []
    all_satisfied_count = 0
    partial_count = 0
    error_count = 0

    for idx, item in enumerate(input_data.proteins):
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(optimize_batch_item,
                    protein=item.protein,
                    organism=item.organism,
                    gc_lo=item.gc_lo,
                    gc_hi=item.gc_hi,
                    cai_threshold=item.cai_threshold,
                    enzymes=item.enzymes,
                    cryptic_splice_threshold=item.cryptic_splice_threshold,
                    organism_domain=item.organism_domain,
                    source_organism=item.source_organism,
                    therapeutic=item.therapeutic,
                    self_protein=item.self_protein,
                ),
                timeout=optimize_timeout,
            )
            results.append({"index": idx, "status": "success", **result})
            if not result.get("failed_predicates", []):
                all_satisfied_count += 1
            else:
                partial_count += 1
        except asyncio.TimeoutError:
            logger.warning("Batch optimize item %d timed out", idx)
            results.append({
                "index": idx, "status": "error",
                "error": f"Item timed out after {optimize_timeout}s",
            })
            error_count += 1
        except Exception as e:
            logger.warning("Batch optimize item %d failed: %s", idx, e)
            results.append({
                "index": idx, "status": "error", "error": str(e),
            })
            error_count += 1

    return BatchOptimizeResponse(
        results=results,
        summary=BatchOptimizeSummary(
            total=n,
            all_satisfied=all_satisfied_count,
            partial=partial_count,
            errors=error_count,
        ),
    )


@_main_router.post("/batch/optimize/fast", response_model=BatchOptimizeResponse)
async def fast_batch_optimize(
    input_data: FastBatchOptimizeInput,
    client_id: str = Depends(verify_api_key),
) -> BatchOptimizeResponse:
    """Fast batch optimization with shared organism/parameters."""
    n = len(input_data.proteins)
    if n > BATCH_OPTIMIZE_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_OPTIMIZE_MAX} proteins",
        )

    _check_batch_rate_limit(client_id, n)

    try:
        opt_results = batch_optimize(
            proteins=input_data.proteins,
            organism=input_data.organism,
            gc_lo=input_data.gc_lo,
            gc_hi=input_data.gc_hi,
            cai_threshold=input_data.cai_threshold,
            enzymes=input_data.enzymes,
            organism_domain=input_data.organism_domain,
            source_organism=input_data.source_organism,
            therapeutic=input_data.therapeutic,
            self_protein=input_data.self_protein,
        )
    except Exception as e:
        logger.warning("Fast batch optimize failed: %s", e)
        raise HTTPException(status_code=500, detail=str(e))

    results: list[dict[str, Any]] = []
    all_satisfied_count = 0
    partial_count = 0
    error_count = 0

    for idx, result in enumerate(opt_results):
        result_dict = {
            "index": idx,
            "status": "success",
            "sequence": result.sequence,
            "protein": result.protein,
            "cai": result.cai,
            "gc_content": result.gc_content,
            "satisfied_predicates": result.satisfied_predicates,
            "failed_predicates": result.failed_predicates,
            "fallback_used": result.fallback_used,
            "organism_domain": input_data.organism_domain,
        }
        results.append(result_dict)
        if not result.failed_predicates:
            all_satisfied_count += 1
        else:
            partial_count += 1

    return BatchOptimizeResponse(
        results=results,
        summary=BatchOptimizeSummary(
            total=n,
            all_satisfied=all_satisfied_count,
            partial=partial_count,
            errors=error_count,
        ),
    )


@_main_router.post("/batch/export", response_model=BatchExportResponse)
async def batch_export(
    input_data: BatchExportInput,
    client_id: str = Depends(verify_api_key),
) -> BatchExportResponse:
    """Export multiple sequences in a single request."""
    n = len(input_data.sequences)
    if n > BATCH_EXPORT_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_EXPORT_MAX} sequences",
        )

    _check_batch_rate_limit(client_id, n)

    results: list[BatchExportResultItem] = []

    for idx, item in enumerate(input_data.sequences):
        try:
            content = await asyncio.wait_for(
                asyncio.to_thread(
                    export_batch_item,
                    sequence=item.sequence,
                    fmt=item.format,
                    identifier=item.identifier,
                    description=item.description,
                    organism=item.organism,
                    locus_name=item.locus_name,
                    definition=item.definition,
                    gene_name=item.gene_name,
                    exon_boundaries=item.exon_boundaries,
                ),
                timeout=BATCH_ITEM_TIMEOUT_S,
            )
            results.append(BatchExportResultItem(
                format=item.format,
                content=content,
                error=None,
            ))
        except asyncio.TimeoutError:
            logger.warning("Batch export item %d timed out", idx)
            results.append(BatchExportResultItem(
                format=item.format,
                content=None,
                error=f"Item timed out after {BATCH_ITEM_TIMEOUT_S}s",
            ))
        except Exception as e:
            logger.warning("Batch export item %d failed: %s", idx, e)
            results.append(BatchExportResultItem(
                format=item.format,
                content=None,
                error=str(e),
            ))

    return BatchExportResponse(results=results)


@_main_router.post("/validate-datasets")
async def validate_datasets(
    client_id: str = Depends(verify_api_key),
    datasets: str | None = None,
    include_cross_organism: bool = True,
    include_optimization_improvement: bool = True,
) -> dict[str, Any]:
    """Validate the optimizer against common biological datasets."""
    from biocompiler.validation.dataset_validation import run_dataset_validation, format_dataset_report_json

    ds_list = datasets.split(",") if datasets else None
    report = await asyncio.to_thread(
        run_dataset_validation,
        datasets=ds_list,
        include_cross_organism=include_cross_organism,
        include_optimization_improvement=include_optimization_improvement,
    )

    return json.loads(format_dataset_report_json(report))


# ═══════════════════════════════════════════════════════════════════
# Protein Analysis Endpoints
# ═══════════════════════════════════════════════════════════════════


@_protein_router.post(
    "/structure/predict",
    response_model=StructurePredictResponse,
    summary="Predict protein structure via ESMFold",
)
async def structure_predict(input_data: StructurePredictInput) -> StructurePredictResponse:
    """Predict protein 3D structure using ESMFold."""
    try:
        from biocompiler.engines.esmfold import is_esmfold_available
    except ImportError:
        raise HTTPException(status_code=503, detail="ESMFold module not available.")

    if not is_esmfold_available():
        raise HTTPException(status_code=503, detail="ESMFold model is not available.")

    try:
        result = predict_structure(protein=input_data.protein)
    except RuntimeError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except TimeoutError:
        raise HTTPException(status_code=504, detail=f"ESMFold prediction timed out after {ESMFOLD_TIMEOUT_S}s.")
    except Exception as e:
        logger.exception("ESMFold prediction failed")
        raise HTTPException(status_code=500, detail=f"Structure prediction failed: {e}")

    return StructurePredictResponse(**result)


@_protein_router.post(
    "/structure/quality",
    response_model=QualityAssessResponse,
    summary="Assess structure quality from PDB",
)
async def structure_quality(input_data: QualityAssessInput) -> QualityAssessResponse:
    """Assess the quality of a protein structure from a PDB string."""
    try:
        from ..structure.quality import compute_structure_quality
    except ImportError:
        raise HTTPException(status_code=503, detail="Structure quality module not available.")

    try:
        result = assess_structure_quality(pdb_string=input_data.pdb_string)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid PDB input: {e}")
    except Exception as e:
        logger.exception("Structure quality assessment failed")
        raise HTTPException(status_code=500, detail=f"Quality assessment failed: {e}")

    return QualityAssessResponse(**result)


@_protein_router.post(
    "/structure/batch",
    response_model=BatchStructureResponse,
    summary="Batch protein structure prediction",
)
async def structure_batch(input_data: BatchStructureInput) -> BatchStructureResponse:
    """Predict protein structures for multiple proteins in one request."""
    try:
        from biocompiler.engines.esmfold import is_esmfold_available
    except ImportError:
        raise HTTPException(status_code=503, detail="ESMFold module not available.")

    if not is_esmfold_available():
        raise HTTPException(status_code=503, detail="ESMFold model is not available.")

    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(status_code=400, detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX}")

    _check_batch_rate_limit("batch_structure", n)

    results: list[BatchStructureResultItem] = []
    total = n
    high_count = 0
    medium_count = 0
    low_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(structure_batch_item, item.protein),
                timeout=ESMFOLD_TIMEOUT_S,
            )
            if "error" in result:
                results.append(BatchStructureResultItem(error=result["error"]))
                error_count += 1
            else:
                results.append(BatchStructureResultItem(
                    pdb_string=result.get("pdb_string"),
                    mean_plddt=result.get("mean_plddt", 0.0),
                    quality_class=result.get("quality_class"),
                ))
                qc = result.get("quality_class")
                if qc in ("very_high", "high"):
                    high_count += 1
                elif qc == "medium":
                    medium_count += 1
                else:
                    low_count += 1
        except asyncio.TimeoutError:
            results.append(BatchStructureResultItem(error=f"Timeout after {ESMFOLD_TIMEOUT_S}s"))
            error_count += 1
        except Exception as e:
            results.append(BatchStructureResultItem(error=str(e)))
            error_count += 1

    return BatchStructureResponse(
        results=results,
        summary={
            "total": total,
            "high_confidence": high_count,
            "medium_confidence": medium_count,
            "low_confidence": low_count,
            "errors": error_count,
        },
    )


# ─── Stability Endpoints ────────────────────────────────────────────

@_protein_router.post(
    "/stability/analyze",
    response_model=StabilityResponse,
    summary="Analyze protein stability",
)
async def stability_analyze(input_data: StabilityInput) -> StabilityResponse:
    """Analyze protein stability using FoldX."""
    try:
        from biocompiler.engines.foldx import empirical_stability
    except ImportError:
        raise HTTPException(status_code=503, detail="FoldX stability module not available.")

    try:
        result = analyze_stability(protein=input_data.protein)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Stability analysis failed")
        raise HTTPException(status_code=500, detail=f"Stability analysis failed: {e}")

    # Serialize Verdict enum value for JSON compatibility
    if "verdict_enum" in result and result["verdict_enum"] is not None:
        result["verdict_enum"] = result["verdict_enum"].value
    return StabilityResponse(**result)


@_protein_router.post(
    "/stability/mutations",
    response_model=MutationScanResponse,
    summary="Scan mutations for stability",
)
async def stability_mutations(input_data: MutationScanInput) -> MutationScanResponse:
    """Scan single-point mutations for stability effects."""
    try:
        from biocompiler.engines.foldx import scan_mutations
    except ImportError:
        raise HTTPException(status_code=503, detail="FoldX mutation scanning module not available.")

    if input_data.positions is not None:
        for p in input_data.positions:
            if p > len(input_data.protein):
                raise HTTPException(status_code=400, detail=f"Position {p} exceeds protein length ({len(input_data.protein)}).")

    try:
        result = scan_stability_mutations(protein=input_data.protein, positions=input_data.positions)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Mutation scanning failed")
        raise HTTPException(status_code=500, detail=f"Mutation scanning failed: {e}")

    return MutationScanResponse(**result)


@_protein_router.post(
    "/stability/batch",
    response_model=BatchStabilityResponse,
    summary="Batch protein stability analysis",
)
async def stability_batch(input_data: BatchStabilityInput) -> BatchStabilityResponse:
    """Analyze protein stability for multiple proteins in one request."""
    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(status_code=400, detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX}")

    _check_batch_rate_limit("batch_stability", n)

    results: list[dict] = []
    total = n
    stable_count = 0
    marginal_count = 0
    unstable_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(stability_batch_item, item.protein),
                timeout=BATCH_ITEM_TIMEOUT_S,
            )
            results.append(result)
            verdict = result.get("verdict", "")
            verdict_enum = assessment_verdict_to_verdict(verdict) if verdict else None
            if verdict_enum is not None:
                result["verdict_enum"] = verdict_enum.value
            if verdict == "STABLE":
                stable_count += 1
            elif verdict == "MARGINAL":
                marginal_count += 1
            elif verdict == "UNSTABLE":
                unstable_count += 1
        except asyncio.TimeoutError:
            results.append({"status": "error", "error": f"Timeout after {BATCH_ITEM_TIMEOUT_S}s"})
            error_count += 1
        except Exception as e:
            results.append({"status": "error", "error": str(e)})
            error_count += 1

    return BatchStabilityResponse(
        results=results,
        summary={
            "total": total,
            "stable": stable_count,
            "marginal": marginal_count,
            "unstable": unstable_count,
            "errors": error_count,
        },
    )


# ─── Solubility Endpoints ───────────────────────────────────────────

@_protein_router.post(
    "/solubility/analyze",
    response_model=SolubilityResponse,
    summary="Analyze protein solubility",
)
async def solubility_analyze(input_data: SolubilityInput) -> SolubilityResponse:
    """Analyze protein solubility using CamSol algorithm."""
    try:
        from biocompiler.engines.camsol import compute_solubility
    except ImportError:
        raise HTTPException(status_code=503, detail="CamSol solubility module not available.")

    try:
        result = analyze_solubility(protein=input_data.protein, pdb_string=input_data.pdb_string)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Solubility analysis failed")
        raise HTTPException(status_code=500, detail=f"Solubility analysis failed: {e}")

    return SolubilityResponse(**result)


@_protein_router.post(
    "/solubility/mutations",
    response_model=MutationScanResponse,
    summary="Find solubility-improving mutations",
)
async def solubility_mutations(input_data: SolubilityInput) -> MutationScanResponse:
    """Find mutations that improve protein solubility."""
    try:
        from biocompiler.engines.camsol import find_solubility_mutations as _fsm
    except ImportError:
        raise HTTPException(status_code=503, detail="CamSol solubility module not available.")

    try:
        result = find_solubility_mutations(protein=input_data.protein)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Solubility mutation scan failed")
        raise HTTPException(status_code=500, detail=f"Solubility mutation scan failed: {e}")

    return MutationScanResponse(**result)


@_protein_router.post(
    "/solubility/batch",
    response_model=BatchSolubilityResponse,
    summary="Batch protein solubility analysis",
)
async def solubility_batch(input_data: BatchSolubilityInput) -> BatchSolubilityResponse:
    """Analyze protein solubility for multiple proteins in one request."""
    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(status_code=400, detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX}")

    _check_batch_rate_limit("batch_solubility", n)

    results: list[dict] = []
    total = n
    high_count = 0
    medium_count = 0
    low_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(solubility_batch_item, item.protein, item.pdb_string),
                timeout=BATCH_ITEM_TIMEOUT_S,
            )
            results.append(result)
            sol_class = result.get("solubility_class", "")
            if sol_class == "high":
                high_count += 1
            elif sol_class == "medium":
                medium_count += 1
            elif sol_class in ("low", "very_low"):
                low_count += 1
        except asyncio.TimeoutError:
            results.append({"status": "error", "error": f"Timeout after {BATCH_ITEM_TIMEOUT_S}s"})
            error_count += 1
        except Exception as e:
            results.append({"status": "error", "error": str(e)})
            error_count += 1

    return BatchSolubilityResponse(
        results=results,
        summary={
            "total": total,
            "high": high_count,
            "medium": medium_count,
            "low": low_count,
            "errors": error_count,
        },
    )


# ─── Immunogenicity Endpoints ───────────────────────────────────────

@_protein_router.post(
    "/immunogenicity/analyze",
    response_model=ImmunogenicityResponse,
    summary="Analyze protein immunogenicity",
)
async def immunogenicity_analyze(input_data: ImmunogenicityInput) -> ImmunogenicityResponse:
    """Analyze protein immunogenicity by predicting T-cell and B-cell epitopes."""
    try:
        from biocompiler.immunogenicity.core import compute_immunogenicity
    except ImportError:
        raise HTTPException(status_code=503, detail="Immunogenicity module not available.")

    try:
        result = analyze_immunogenicity(protein=input_data.protein, mhc_alleles=input_data.mhc_alleles)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Immunogenicity analysis failed")
        raise HTTPException(status_code=500, detail=f"Immunogenicity analysis failed: {e}")

    return ImmunogenicityResponse(**result)


@_protein_router.post(
    "/immunogenicity/deimmunize",
    response_model=DeimmunizeResponse,
    summary="Deimmunize a protein",
)
async def immunogenicity_deimmunize(input_data: DeimmunizeInput) -> DeimmunizeResponse:
    """Deimmunize a protein by introducing conservative amino acid substitutions."""
    try:
        from biocompiler.immunogenicity.deimmunization import deimmunize
    except ImportError:
        raise HTTPException(status_code=503, detail="Deimmunization module not available.")

    try:
        result = deimmunize_protein(
            protein=input_data.protein,
            organism=input_data.organism,
            target_score=input_data.target_score,
            max_mutations=input_data.max_mutations,
            blosum62_min=input_data.blosum62_min,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Deimmunization failed")
        raise HTTPException(status_code=500, detail=f"Deimmunization failed: {e}")

    return DeimmunizeResponse(**result)


@_protein_router.post(
    "/immunogenicity/batch",
    response_model=BatchImmunogenicityResponse,
    summary="Batch protein immunogenicity analysis",
)
async def immunogenicity_batch(input_data: BatchImmunogenicityInput) -> BatchImmunogenicityResponse:
    """Analyze protein immunogenicity for multiple proteins in one request."""
    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(status_code=400, detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX}")

    _check_batch_rate_limit("batch_immunogenicity", n)

    results: list[dict] = []
    total = n
    low_count = 0
    moderate_count = 0
    high_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(immunogenicity_batch_item, item.protein, item.mhc_alleles),
                timeout=BATCH_ITEM_TIMEOUT_S,
            )
            results.append(result)
            imm_class = result.get("immunogenicity_class", "")
            if imm_class == "low":
                low_count += 1
            elif imm_class == "moderate":
                moderate_count += 1
            elif imm_class in ("high", "very_high"):
                high_count += 1
        except asyncio.TimeoutError:
            results.append({"status": "error", "error": f"Timeout after {BATCH_ITEM_TIMEOUT_S}s"})
            error_count += 1
        except Exception as e:
            results.append({"status": "error", "error": str(e)})
            error_count += 1

    return BatchImmunogenicityResponse(
        results=results,
        summary={
            "total": total,
            "low": low_count,
            "moderate": moderate_count,
            "high": high_count,
            "errors": error_count,
        },
    )


# ─── Full Assessment Endpoint ────────────────────────────────────────

@_protein_router.post(
    "/assessment/full",
    response_model=FullAssessmentResponse,
    summary="Full protein assessment",
)
async def assessment_full(input_data: FullAssessmentInput) -> FullAssessmentResponse:
    """Run a full protein assessment combining all analyses."""
    try:
        result = full_assessment(
            protein=input_data.protein,
            organism=input_data.organism,
            pdb_string=input_data.pdb_string,
            run_structure=input_data.run_structure,
            run_stability=input_data.run_stability,
            run_solubility=input_data.run_solubility,
            run_immunogenicity=input_data.run_immunogenicity,
        )
        # Serialize Verdict enum values for JSON compatibility
        response_data = dict(result)
        if "overall_verdict_enum" in response_data and response_data["overall_verdict_enum"] is not None:
            response_data["overall_verdict_enum"] = response_data["overall_verdict_enum"].value
        if response_data.get("stability") and "verdict_enum" in response_data["stability"]:
            response_data["stability"]["verdict_enum"] = response_data["stability"]["verdict_enum"].value
        return FullAssessmentResponse(**response_data)
    except Exception as e:
        logger.exception("Full assessment failed")
        raise HTTPException(status_code=500, detail=f"Full assessment failed: {e}")
