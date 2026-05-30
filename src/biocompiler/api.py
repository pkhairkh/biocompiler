"""
BioCompiler REST API — FastAPI Web Service

Production-grade REST API for integration with other bioinformatics tools.

Endpoints:
- POST /check            — Type-check a DNA sequence
- POST /optimize         — Optimize DNA for a target protein
- POST /verify           — Verify a certificate
- POST /scan             — Scan a sequence for motifs
- POST /export/fasta     — Export sequence as FASTA
- POST /export/genbank   — Export sequence as GenBank
- GET  /organisms        — List supported organisms
- GET  /predicates       — List registered predicates
- GET  /health           — Health check

Batch Endpoints:
- POST /batch/check      — Type-check multiple sequences in one request (max 50)
- POST /batch/optimize   — Optimize multiple proteins in one request (max 20)
- POST /batch/export     — Export multiple sequences in one request (max 50)

Security:
- API key authentication (optional, set BIOCOMPILER_API_KEY env var)
- Rate limiting (60 requests/minute by default, configurable)
  Batch requests consume rate-limit units per item (e.g., 50 sequences = 50 units)
- CORS with configurable origins

All endpoints accept and return JSON. Certificate data is embedded
directly in responses for seamless pipeline integration.

Batch Processing:
- Each batch item is processed independently; one failure does not affect others
- Per-item timeout prevents a single slow item from blocking the entire batch
- Batch size limits are enforced with 400 responses for oversized requests
"""

import asyncio
import logging
import os
import time
from collections import defaultdict
from datetime import datetime, timezone
from typing import Any, Optional

from fastapi import FastAPI, HTTPException, Query, Request, Security, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import APIKeyHeader
from pydantic import BaseModel, Field, field_validator

from .scanner import scan_sequence, gc_content, validate_dna_sequence
from .translation import translate, compute_cai, find_orfs
from .splicing import compute_splice_isoforms
from .type_system import evaluate_all_predicates
from .certificate import generate_certificate, verify_certificate
from .optimization import optimize_sequence, OptimizationResult
from .export import export_fasta, export_genbank, export_genbank_with_certificate
from .types import Verdict, Certificate
from .constants import RESTRICTION_ENZYMES
from .organisms import SUPPORTED_ORGANISMS, CODON_USAGE_TABLES
from .exceptions import (
    BioCompilerError, InvalidSequenceError, CertificateGenerationError,
    CertificateVerificationError, UnsupportedOrganismError, InvalidProteinError,
)

logger = logging.getLogger(__name__)

# ─── API Key Authentication ─────────────────────────────────────────

API_KEY_NAME = "X-API-Key"
_api_key_header = APIKeyHeader(name=API_KEY_NAME, auto_error=False)

# The API key is optional: if BIOCOMPILER_API_KEY is not set, auth is disabled.
# Supports multiple keys via comma-separated BIOCOMPILER_API_KEYS env var,
# or single key via BIOCOMPILER_API_KEY for backward compatibility.
_CONFIGURED_API_KEYS: set[str] = set()

_single_key = os.environ.get("BIOCOMPILER_API_KEY", "")
_multi_keys = os.environ.get("BIOCOMPILER_API_KEYS", "")

if _multi_keys:
    _CONFIGURED_API_KEYS = {k.strip() for k in _multi_keys.split(",") if k.strip()}
elif _single_key:
    _CONFIGURED_API_KEYS = {_single_key}

# ─── Rate Limiting ──────────────────────────────────────────────────

RATE_LIMIT_RPM = int(os.environ.get("BIOCOMPILER_RATE_LIMIT", "60"))  # requests per minute
_rate_limit_store: dict[str, list[float]] = defaultdict(list)


def _check_rate_limit(client_id: str) -> None:
    """Enforce per-client rate limiting (sliding window)."""
    now = time.monotonic()
    window = _rate_limit_store[client_id]
    # Remove timestamps older than 60 seconds
    _rate_limit_store[client_id] = [t for t in window if now - t < 60.0]
    window = _rate_limit_store[client_id]
    if len(window) >= RATE_LIMIT_RPM:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded: {RATE_LIMIT_RPM} requests/minute. "
                   f"Retry after {60.0 - (now - window[0]):.0f} seconds.",
        )
    window.append(now)


def _check_batch_rate_limit(client_id: str, item_count: int) -> None:
    """
    Enforce per-client rate limiting for batch requests.

    Each item in the batch consumes one rate-limit unit. We pre-check
    whether the entire batch can be accommodated, then record all units.
    This prevents partial batch execution when rate-limited mid-way.
    """
    now = time.monotonic()
    window = _rate_limit_store[client_id]
    # Remove timestamps older than 60 seconds
    _rate_limit_store[client_id] = [t for t in window if now - t < 60.0]
    window = _rate_limit_store[client_id]
    remaining = RATE_LIMIT_RPM - len(window)
    if item_count > remaining:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded: batch has {item_count} items but only "
                   f"{remaining} rate-limit units remaining this minute. "
                   f"Limit: {RATE_LIMIT_RPM} requests/minute.",
        )
    # Record all items
    window.extend([now] * item_count)


async def verify_api_key(api_key: str = Security(_api_key_header)) -> str:
    """
    Verify the API key if authentication is configured.

    If BIOCOMPILER_API_KEY / BIOCOMPILER_API_KEYS is not set, auth is disabled.
    If set, the request must include a matching X-API-Key header.
    Supports multiple API keys for key rotation scenarios.
    """
    if not _CONFIGURED_API_KEYS:
        return "anonymous"  # Auth disabled

    if api_key is None or api_key not in _CONFIGURED_API_KEYS:
        raise HTTPException(
            status_code=401,
            detail="Invalid or missing API key. Set X-API-Key header.",
        )
    return api_key


# ─── Pydantic Models ──────────────────────────────────────────────

class SequenceInput(BaseModel):
    """DNA sequence input."""
    sequence: str = Field(..., description="DNA sequence (ACGTN)")
    organism: str = Field("Homo_sapiens", description="Target organism")
    exon_boundaries: Optional[list[tuple[int, int]]] = Field(
        None, description="Exon boundaries as [(start, end), ...]"
    )
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.5, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to check")
    cellular_context: str = Field("HEK293T", description="Cellular context for splicing")

    @field_validator("sequence")
    @classmethod
    def validate_seq(cls, v: str) -> str:
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        if v not in SUPPORTED_ORGANISMS:
            raise ValueError(f"Unsupported organism: {v}. Supported: {SUPPORTED_ORGANISMS}")
        return v


class ProteinInput(BaseModel):
    """Protein sequence input for optimization."""
    protein: str = Field(..., description="Target protein sequence (single-letter codes)")
    organism: str = Field("Homo_sapiens", description="Target organism")
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.2, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to avoid")
    cryptic_splice_threshold: float = Field(3.0, description="Cryptic splice site threshold")

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        v = v.upper()
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(v) - valid_aas
        if invalid:
            raise ValueError(f"Invalid amino acids: {invalid}")
        return v


class CertificateInput(BaseModel):
    """Certificate verification input."""
    certificate: dict = Field(..., description="Certificate as JSON dict")


class ExportFastaInput(BaseModel):
    """FASTA export input."""
    sequence: str = Field(..., description="DNA sequence")
    identifier: str = Field("BioCompiler_design", description="Sequence identifier")
    description: str = Field("", description="Description line")
    organism: str = Field("Homo_sapiens", description="Source organism")


class ExportGenbankInput(BaseModel):
    """GenBank export input."""
    sequence: str = Field(..., description="DNA sequence")
    locus_name: str = Field("BIOCOMPILER", description="LOCUS name (max 16 chars)")
    definition: str = Field("BioCompiler designed sequence", description="DEFINITION line")
    organism: str = Field("Homo_sapiens", description="Source organism")
    gene_name: Optional[str] = Field(None, description="Gene name")
    exon_boundaries: Optional[list[tuple[int, int]]] = Field(None)
    certificate: Optional[dict] = Field(None, description="Certificate dict to embed")


class ScanInput(BaseModel):
    """Sequence scan input."""
    sequence: str = Field(..., description="DNA sequence")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to scan for")
    find_orfs: bool = Field(False, description="Find open reading frames")


# ─── Response Models ──────────────────────────────────────────────

class TypeCheckResponse(BaseModel):
    sequence_length: int
    gc_content: float
    protein: str
    results: list[dict]
    overall_verdict: str
    certificate: Optional[dict] = None


class OptimizeResponse(BaseModel):
    sequence: str
    protein: str
    cai: float
    gc_content: float
    satisfied_predicates: list[str]
    failed_predicates: list[str]
    fallback_used: bool


class VerifyResponse(BaseModel):
    status: str
    failure_reasons: list[str]


class ScanResponse(BaseModel):
    sequence_length: int
    tokens: list[dict]
    orfs: Optional[list[dict]] = None


class OrganismResponse(BaseModel):
    organisms: list[dict]


class PredicateResponse(BaseModel):
    predicates: list[str]


class HealthResponse(BaseModel):
    status: str
    version: str
    timestamp: str
    auth_enabled: bool = False
    rate_limit_rpm: int = 60


# ─── Batch Pydantic Models ────────────────────────────────────────

BATCH_CHECK_MAX = 50
BATCH_OPTIMIZE_MAX = 20
BATCH_EXPORT_MAX = 50
BATCH_ITEM_TIMEOUT_S = int(os.environ.get("BIOCOMPILER_BATCH_ITEM_TIMEOUT", "30"))


class BatchCheckItem(BaseModel):
    """Single item in a batch type-check request."""
    sequence: str = Field(..., description="DNA sequence (ACGTN)")
    organism: str = Field("Homo_sapiens", description="Target organism")
    exon_boundaries: Optional[list[tuple[int, int]]] = Field(
        None, description="Exon boundaries as [(start, end), ...]"
    )
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.5, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to check")
    cellular_context: str = Field("HEK293T", description="Cellular context for splicing")

    @field_validator("sequence")
    @classmethod
    def validate_seq(cls, v: str) -> str:
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        if v not in SUPPORTED_ORGANISMS:
            raise ValueError(f"Unsupported organism: {v}. Supported: {SUPPORTED_ORGANISMS}")
        return v


class BatchCheckInput(BaseModel):
    """Input for batch type-check endpoint."""
    sequences: list[BatchCheckItem] = Field(
        ..., description=f"List of sequences to type-check (max {BATCH_CHECK_MAX})"
    )


class BatchCheckSummary(BaseModel):
    """Summary statistics for batch type-check results."""
    total: int = Field(..., description="Total number of sequences processed")
    pass_: int = Field(..., alias="pass", description="Number of sequences with PASS verdict")
    fail: int = Field(..., description="Number of sequences with FAIL verdict")
    uncertain: int = Field(..., description="Number of sequences with UNCERTAIN verdict")
    errors: int = Field(..., description="Number of sequences that encountered errors")

    model_config = {"populate_by_name": True}


class BatchCheckResponse(BaseModel):
    """Response for batch type-check endpoint."""
    results: list[dict] = Field(..., description="Per-item results")
    summary: BatchCheckSummary = Field(..., description="Aggregate summary")


class BatchOptimizeItem(BaseModel):
    """Single item in a batch optimize request."""
    protein: str = Field(..., description="Target protein sequence (single-letter codes)")
    organism: str = Field("Homo_sapiens", description="Target organism")
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.2, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to avoid")
    cryptic_splice_threshold: float = Field(3.0, description="Cryptic splice site threshold")

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        v = v.upper()
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(v) - valid_aas
        if invalid:
            raise ValueError(f"Invalid amino acids: {invalid}")
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        if v not in SUPPORTED_ORGANISMS:
            raise ValueError(f"Unsupported organism: {v}. Supported: {SUPPORTED_ORGANISMS}")
        return v


class BatchOptimizeInput(BaseModel):
    """Input for batch optimize endpoint."""
    proteins: list[BatchOptimizeItem] = Field(
        ..., description=f"List of proteins to optimize (max {BATCH_OPTIMIZE_MAX})"
    )


class BatchOptimizeSummary(BaseModel):
    """Summary statistics for batch optimize results."""
    total: int = Field(..., description="Total number of proteins processed")
    all_satisfied: int = Field(..., description="Proteins where all predicates are satisfied")
    partial: int = Field(..., description="Proteins where some predicates failed")
    errors: int = Field(..., description="Proteins that encountered errors")


class BatchOptimizeResponse(BaseModel):
    """Response for batch optimize endpoint."""
    results: list[dict] = Field(..., description="Per-item results")
    summary: BatchOptimizeSummary = Field(..., description="Aggregate summary")


class BatchExportItem(BaseModel):
    """Single item in a batch export request."""
    sequence: str = Field(..., description="DNA sequence")
    format: str = Field("fasta", description="Export format: 'fasta' or 'genbank'")
    identifier: str = Field("BioCompiler_design", description="Sequence identifier (FASTA)")
    description: str = Field("", description="Description line (FASTA)")
    organism: str = Field("Homo_sapiens", description="Source organism")
    locus_name: str = Field("BIOCOMPILER", description="LOCUS name (GenBank, max 16 chars)")
    definition: str = Field("BioCompiler designed sequence", description="DEFINITION line (GenBank)")
    gene_name: Optional[str] = Field(None, description="Gene name (GenBank)")
    exon_boundaries: Optional[list[tuple[int, int]]] = Field(None, description="Exon boundaries (GenBank)")

    @field_validator("format")
    @classmethod
    def validate_format(cls, v: str) -> str:
        v = v.lower()
        if v not in ("fasta", "genbank"):
            raise ValueError(f"Unsupported format: {v}. Supported: fasta, genbank")
        return v


class BatchExportInput(BaseModel):
    """Input for batch export endpoint."""
    sequences: list[BatchExportItem] = Field(
        ..., description=f"List of sequences to export (max {BATCH_EXPORT_MAX})"
    )


class BatchExportResultItem(BaseModel):
    """Result for a single export item."""
    format: str = Field(..., description="Export format used")
    content: Optional[str] = Field(None, description="Exported content (None on error)")
    error: Optional[str] = Field(None, description="Error message if export failed")


class BatchExportResponse(BaseModel):
    """Response for batch export endpoint."""
    results: list[BatchExportResultItem] = Field(..., description="Per-item export results")


# ─── Batch Helper Functions ────────────────────────────────────────

def _type_check_single(item: BatchCheckItem) -> dict[str, Any]:
    """
    Process a single type-check item for batch processing.

    Returns a dict matching the TypeCheckResponse structure.
    Raises exceptions on error (caller handles isolation).
    """
    seq = item.sequence.upper()
    exon_boundaries = item.exon_boundaries or [(0, len(seq))]

    # Type check
    results = evaluate_all_predicates(
        seq=seq,
        known_exon_boundaries=exon_boundaries,
        organism=item.organism,
        gc_lo=item.gc_lo,
        gc_hi=item.gc_hi,
        cai_threshold=item.cai_threshold,
        enzymes=item.enzymes,
    )

    # Overall verdict
    verdicts = [r.verdict for r in results]
    overall = Verdict.PASS
    for v in verdicts:
        if overall == Verdict.FAIL or v == Verdict.FAIL:
            overall = Verdict.FAIL
        elif overall == Verdict.UNCERTAIN or v == Verdict.UNCERTAIN:
            overall = Verdict.UNCERTAIN

    # Build response
    result_dicts = [
        {
            "predicate": r.predicate,
            "verdict": r.verdict.value,
            "violation": r.violation,
            "knowledge_gap": r.knowledge_gap,
        }
        for r in results
    ]

    # Protein translation
    coding_seq = "".join(seq[start:end] for start, end in exon_boundaries)
    protein = translate(coding_seq)

    # Generate certificate if all PASS
    cert_dict = None
    if overall == Verdict.PASS:
        try:
            cert = generate_certificate(
                coding_seq, results,
                {
                    "organism": item.organism,
                    "exon_boundaries": exon_boundaries,
                    "gc_lo": item.gc_lo,
                    "gc_hi": item.gc_hi,
                    "cai_threshold": item.cai_threshold,
                    "enzymes": item.enzymes or list(RESTRICTION_ENZYMES.keys()),
                    "cell_type": item.cellular_context,
                },
            )
            cert_dict = cert.to_dict()
        except CertificateGenerationError:
            pass

    return {
        "sequence_length": len(seq),
        "gc_content": gc_content(seq),
        "protein": protein[:100] + ("..." if len(protein) > 100 else ""),
        "results": result_dicts,
        "overall_verdict": overall.value,
        "certificate": cert_dict,
    }


def _optimize_single(item: BatchOptimizeItem) -> dict[str, Any]:
    """
    Process a single optimization item for batch processing.

    Returns a dict matching the OptimizeResponse structure.
    Raises exceptions on error (caller handles isolation).
    """
    result = optimize_sequence(
        target_protein=item.protein,
        organism=item.organism,
        gc_lo=item.gc_lo,
        gc_hi=item.gc_hi,
        cai_threshold=item.cai_threshold,
        restriction_sites=item.enzymes,
        cryptic_splice_threshold=item.cryptic_splice_threshold,
    )

    return {
        "sequence": result.sequence,
        "protein": result.protein,
        "cai": result.cai,
        "gc_content": result.gc_content,
        "satisfied_predicates": result.satisfied_predicates,
        "failed_predicates": result.failed_predicates,
        "fallback_used": result.fallback_used,
    }


def _export_single(item: BatchExportItem) -> str:
    """
    Process a single export item for batch processing.

    Returns the exported content string.
    Raises exceptions on error (caller handles isolation).
    """
    if item.format == "fasta":
        return export_fasta(
            sequence=item.sequence,
            identifier=item.identifier,
            description=item.description,
            organism=item.organism,
        )
    elif item.format == "genbank":
        return export_genbank(
            sequence=item.sequence,
            locus_name=item.locus_name,
            definition=item.definition,
            organism=item.organism,
            gene_name=item.gene_name,
            exon_boundaries=item.exon_boundaries,
        )
    else:
        raise ValueError(f"Unsupported export format: {item.format}")


# ─── FastAPI Application ──────────────────────────────────────────

def create_app() -> FastAPI:
    """
    Create and configure the FastAPI application.

    Single-item endpoints:
    - POST /check            — Type-check a DNA sequence
    - POST /optimize         — Optimize DNA for a target protein
    - POST /verify           — Verify a certificate
    - POST /scan             — Scan a sequence for motifs
    - POST /export/fasta     — Export sequence as FASTA
    - POST /export/genbank   — Export sequence as GenBank
    - GET  /organisms        — List supported organisms
    - GET  /predicates       — List registered predicates
    - GET  /health           — Health check

    Batch endpoints:
    - POST /batch/check      — Type-check multiple sequences (max 50)
    - POST /batch/optimize   — Optimize multiple proteins (max 20)
    - POST /batch/export     — Export multiple sequences (max 50)
    """
    from . import __version__

    # CORS origins from environment (default: allow all for development)
    cors_origins = os.environ.get("BIOCOMPILER_CORS_ORIGINS", "*").split(",")

    app = FastAPI(
        title="BioCompiler API",
        description="Machine-verified gene design REST API with authentication and rate limiting",
        version=__version__,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # CORS middleware for cross-origin access from web tools
    app.add_middleware(
        CORSMiddleware,
        allow_origins=cors_origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Rate limiting middleware
    @app.middleware("http")
    async def rate_limit_middleware(request: Request, call_next):
        client_id = request.client.host if request.client else "unknown"
        try:
            _check_rate_limit(client_id)
        except HTTPException as e:
            from fastapi.responses import JSONResponse
            return JSONResponse(status_code=e.status_code, content={"detail": e.detail})
        response = await call_next(request)
        return response

    # ─── Endpoints ──────────────────────────────────────────────

    @app.get("/health", response_model=HealthResponse)
    async def health_check():
        """Health check endpoint (no auth required)."""
        return HealthResponse(
            status="healthy",
            version=__version__,
            timestamp=datetime.now(timezone.utc).isoformat(),
            auth_enabled=bool(_CONFIGURED_API_KEYS),
            rate_limit_rpm=RATE_LIMIT_RPM,
        )

    @app.get("/organisms", response_model=OrganismResponse)
    async def list_organisms():
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

    @app.get("/predicates", response_model=PredicateResponse)
    async def list_predicates():
        """List all registered type predicates."""
        from .type_system import registry
        return PredicateResponse(predicates=registry.names())

    @app.get("/enzymes")
    async def list_enzymes():
        """List all known restriction enzymes."""
        return {
            "enzymes": {
                name: site for name, site in RESTRICTION_ENZYMES.items()
            }
        }

    @app.post("/check", response_model=TypeCheckResponse)
    async def check_sequence(input_data: SequenceInput, client_id: str = Depends(verify_api_key)):
        """
        Type-check a DNA sequence against all registered predicates.

        Optionally generates a guarantee certificate if all predicates pass.
        """
        try:
            seq = input_data.sequence.upper()
            exon_boundaries = input_data.exon_boundaries or [(0, len(seq))]

            # Type check
            results = evaluate_all_predicates(
                seq=seq,
                known_exon_boundaries=exon_boundaries,
                organism=input_data.organism,
                gc_lo=input_data.gc_lo,
                gc_hi=input_data.gc_hi,
                cai_threshold=input_data.cai_threshold,
                enzymes=input_data.enzymes,
            )

            # Overall verdict
            verdicts = [r.verdict for r in results]
            overall = Verdict.PASS
            for v in verdicts:
                if overall == Verdict.FAIL or v == Verdict.FAIL:
                    overall = Verdict.FAIL
                elif overall == Verdict.UNCERTAIN or v == Verdict.UNCERTAIN:
                    overall = Verdict.UNCERTAIN

            # Build response
            result_dicts = [
                {
                    "predicate": r.predicate,
                    "verdict": r.verdict.value,
                    "violation": r.violation,
                    "knowledge_gap": r.knowledge_gap,
                }
                for r in results
            ]

            # Protein translation
            coding_seq = "".join(seq[start:end] for start, end in exon_boundaries)
            protein = translate(coding_seq)

            # Generate certificate if all PASS
            cert_dict = None
            if overall == Verdict.PASS:
                try:
                    cert = generate_certificate(
                        coding_seq, results,
                        {
                            "organism": input_data.organism,
                            "exon_boundaries": exon_boundaries,
                            "gc_lo": input_data.gc_lo,
                            "gc_hi": input_data.gc_hi,
                            "cai_threshold": input_data.cai_threshold,
                            "enzymes": input_data.enzymes or list(RESTRICTION_ENZYMES.keys()),
                            "cell_type": input_data.cellular_context,
                        },
                    )
                    cert_dict = cert.to_dict()
                except CertificateGenerationError:
                    pass

            return TypeCheckResponse(
                sequence_length=len(seq),
                gc_content=gc_content(seq),
                protein=protein[:100] + ("..." if len(protein) > 100 else ""),
                results=result_dicts,
                overall_verdict=overall.value,
                certificate=cert_dict,
            )
        except (InvalidSequenceError, UnsupportedOrganismError) as e:
            raise HTTPException(status_code=400, detail=str(e))
        except BioCompilerError as e:
            raise HTTPException(status_code=422, detail=str(e))

    @app.post("/optimize", response_model=OptimizeResponse)
    async def optimize(input_data: ProteinInput, client_id: str = Depends(verify_api_key)):
        """
        Optimize a DNA sequence for a target protein.

        Uses z3 constraint solver with greedy fallback to find
        a sequence that satisfies all type predicates.
        """
        try:
            result = optimize_sequence(
                target_protein=input_data.protein,
                organism=input_data.organism,
                gc_lo=input_data.gc_lo,
                gc_hi=input_data.gc_hi,
                cai_threshold=input_data.cai_threshold,
                restriction_sites=input_data.enzymes,
                cryptic_splice_threshold=input_data.cryptic_splice_threshold,
            )

            return OptimizeResponse(
                sequence=result.sequence,
                protein=result.protein,
                cai=result.cai,
                gc_content=result.gc_content,
                satisfied_predicates=result.satisfied_predicates,
                failed_predicates=result.failed_predicates,
                fallback_used=result.fallback_used,
            )
        except (InvalidProteinError, UnsupportedOrganismError) as e:
            raise HTTPException(status_code=400, detail=str(e))
        except BioCompilerError as e:
            raise HTTPException(status_code=422, detail=str(e))

    @app.post("/verify", response_model=VerifyResponse)
    async def verify(input_data: CertificateInput, client_id: str = Depends(verify_api_key)):
        """
        Independently verify a guarantee certificate.

        Re-evaluates every predicate from scratch using only the
        sequence and parameters embedded in the certificate.
        """
        try:
            status, failures = verify_certificate(input_data.certificate)
            return VerifyResponse(status=status, failure_reasons=failures)
        except (CertificateVerificationError, ValueError, KeyError) as e:
            raise HTTPException(status_code=400, detail=f"Verification error: {e}")
        except Exception as e:
            logger.exception("Unexpected error during certificate verification")
            raise HTTPException(status_code=500, detail="Internal verification error")

    @app.post("/scan", response_model=ScanResponse)
    async def scan(input_data: ScanInput, client_id: str = Depends(verify_api_key)):
        """
        Scan a DNA sequence for biological motifs.

        Detects splice sites, start/stop codons, Kozak contexts,
        instability motifs, and restriction enzyme sites.
        """
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
                orfs = find_orfs(seq)
                orf_dicts = orfs

            return ScanResponse(
                sequence_length=len(seq),
                tokens=token_dicts,
                orfs=orf_dicts,
            )
        except InvalidSequenceError as e:
            raise HTTPException(status_code=400, detail=str(e))

    @app.post("/export/fasta")
    async def export_fasta_endpoint(input_data: ExportFastaInput, client_id: str = Depends(verify_api_key)):
        """Export a sequence in FASTA format."""
        try:
            fasta = export_fasta(
                sequence=input_data.sequence,
                identifier=input_data.identifier,
                description=input_data.description,
                organism=input_data.organism,
            )
            return {"format": "fasta", "content": fasta}
        except Exception as e:
            raise HTTPException(status_code=400, detail=str(e))

    @app.post("/export/genbank")
    async def export_genbank_endpoint(input_data: ExportGenbankInput, client_id: str = Depends(verify_api_key)):
        """Export a sequence in GenBank format with optional certificate embedding."""
        try:
            cert = None
            if input_data.certificate:
                cert = Certificate.from_dict(input_data.certificate)

            genbank = export_genbank(
                sequence=input_data.sequence,
                locus_name=input_data.locus_name,
                definition=input_data.definition,
                organism=input_data.organism,
                gene_name=input_data.gene_name,
                exon_boundaries=input_data.exon_boundaries,
                certificate=cert,
            )
            return {"format": "genbank", "content": genbank}
        except Exception as e:
            raise HTTPException(status_code=400, detail=str(e))

    # ─── Batch Endpoints ─────────────────────────────────────────

    @app.post("/batch/check", response_model=BatchCheckResponse)
    async def batch_check(
        input_data: BatchCheckInput,
        client_id: str = Depends(verify_api_key),
    ):
        """
        Type-check multiple DNA sequences in a single request.

        Each sequence is processed independently — one failure does not
        affect others. Per-item timeout prevents a single slow item from
        blocking the entire batch.

        Rate limiting: each sequence consumes one rate-limit unit.
        Maximum batch size: 50 sequences.
        """
        n = len(input_data.sequences)
        if n > BATCH_CHECK_MAX:
            raise HTTPException(
                status_code=400,
                detail=f"Batch size {n} exceeds maximum of {BATCH_CHECK_MAX} sequences",
            )

        # Per-item rate limiting
        request_client = client_id
        _check_batch_rate_limit(request_client, n)

        results: list[dict[str, Any]] = []
        pass_count = 0
        fail_count = 0
        uncertain_count = 0
        error_count = 0

        for idx, item in enumerate(input_data.sequences):
            try:
                result = await asyncio.wait_for(
                    asyncio.to_thread(_type_check_single, item),
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
                    "index": idx,
                    "status": "error",
                    "error": f"Item timed out after {BATCH_ITEM_TIMEOUT_S}s",
                })
                error_count += 1
            except Exception as e:
                logger.warning("Batch check item %d failed: %s", idx, e)
                results.append({
                    "index": idx,
                    "status": "error",
                    "error": str(e),
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

    @app.post("/batch/optimize", response_model=BatchOptimizeResponse)
    async def batch_optimize(
        input_data: BatchOptimizeInput,
        client_id: str = Depends(verify_api_key),
    ):
        """
        Optimize multiple proteins in a single request.

        Each protein is optimized independently — one failure does not
        affect others. Per-item timeout prevents a single slow optimization
        from blocking the entire batch.

        Rate limiting: each protein consumes one rate-limit unit.
        Maximum batch size: 20 proteins (optimization is expensive).
        """
        n = len(input_data.proteins)
        if n > BATCH_OPTIMIZE_MAX:
            raise HTTPException(
                status_code=400,
                detail=f"Batch size {n} exceeds maximum of {BATCH_OPTIMIZE_MAX} proteins",
            )

        # Per-item rate limiting
        request_client = client_id
        _check_batch_rate_limit(request_client, n)

        # Use a longer timeout for optimization
        optimize_timeout = max(BATCH_ITEM_TIMEOUT_S, 60)

        results: list[dict[str, Any]] = []
        all_satisfied_count = 0
        partial_count = 0
        error_count = 0

        for idx, item in enumerate(input_data.proteins):
            try:
                result = await asyncio.wait_for(
                    asyncio.to_thread(_optimize_single, item),
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
                    "index": idx,
                    "status": "error",
                    "error": f"Item timed out after {optimize_timeout}s",
                })
                error_count += 1
            except Exception as e:
                logger.warning("Batch optimize item %d failed: %s", idx, e)
                results.append({
                    "index": idx,
                    "status": "error",
                    "error": str(e),
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

    @app.post("/batch/export", response_model=BatchExportResponse)
    async def batch_export(
        input_data: BatchExportInput,
        client_id: str = Depends(verify_api_key),
    ):
        """
        Export multiple sequences in a single request.

        Each sequence is exported independently — one failure does not
        affect others. Supports both FASTA and GenBank formats per item.

        Rate limiting: each sequence consumes one rate-limit unit.
        Maximum batch size: 50 sequences.
        """
        n = len(input_data.sequences)
        if n > BATCH_EXPORT_MAX:
            raise HTTPException(
                status_code=400,
                detail=f"Batch size {n} exceeds maximum of {BATCH_EXPORT_MAX} sequences",
            )

        # Per-item rate limiting
        request_client = client_id
        _check_batch_rate_limit(request_client, n)

        results: list[BatchExportResultItem] = []

        for idx, item in enumerate(input_data.sequences):
            try:
                content = await asyncio.wait_for(
                    asyncio.to_thread(_export_single, item),
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

    return app


# Create the default app instance
app = create_app()
