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
- GET  /provenance       — Query/list provenance records
- GET  /provenance/{id}  — Retrieve a specific provenance record

Batch Endpoints:
- POST /batch/check      — Type-check multiple sequences in one request (max 50)
- POST /batch/optimize   — Optimize multiple proteins in one request (max 20)
- POST /batch/export     — Export multiple sequences in one request (max 50)

Protein Analysis Endpoints (mounted at /protein/):
- POST /protein/structure/predict        — Predict protein structure via ESMFold
- POST /protein/structure/batch          — Batch structure prediction (max 20)
- POST /protein/structure/quality         — Assess structure quality from PDB
- POST /protein/stability/analyze         — Analyze protein stability
- POST /protein/stability/mutations        — Scan mutations for stability
- POST /protein/stability/batch           — Batch stability analysis (max 20)
- POST /protein/solubility/analyze        — Analyze protein solubility
- POST /protein/solubility/mutations       — Find solubility-improving mutations
- POST /protein/solubility/batch          — Batch solubility analysis (max 20)
- POST /protein/immunogenicity/analyze    — Analyze immunogenicity
- POST /protein/immunogenicity/deimmunize — Deimmunize a protein
- POST /protein/immunogenicity/batch      — Batch immunogenicity analysis (max 20)
- POST /protein/assessment/full           — Full protein assessment

Security:
- API key authentication (default: enabled, auto-generated if BIOCOMPILER_API_KEY is not set)
- BIOCOMPILER_AUTH_MODE: "required" (default), "optional", "disabled"
  - "required": Unauthenticated requests get HTTP 401
  - "optional": Unauthenticated requests allowed with X-Auth-Warning response header
  - "disabled": Auth completely disabled (DANGEROUS — use only with BIOCOMPILER_API_KEY=disabled or --no-auth)
- If BIOCOMPILER_API_KEY is not set, a random key is generated on first startup,
  printed to console, and saved to ~/.biocompiler/api_key for reuse
- If BIOCOMPILER_API_KEY=disabled, auth is explicitly disabled with a warning
- --no-auth CLI flag disables auth for local development (with warning)
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
import json
import logging
import os
import secrets
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from fastapi import APIRouter, FastAPI, HTTPException, Request, Security, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import APIKeyHeader
from pydantic import BaseModel, Field, field_validator

from .scanner import scan_sequence, gc_content, validate_dna_sequence
from .translation import translate, compute_cai, find_orfs
from .splicing import compute_splice_isoforms
from .type_system import evaluate_all_predicates
from .certificate import generate_certificate, verify_certificate
from .optimization import optimize_sequence, batch_optimize
from .export import export_fasta, export_genbank
from .types import Verdict, Certificate
from .constants import RESTRICTION_ENZYMES
from .organisms import SUPPORTED_ORGANISMS, CODON_USAGE_TABLES, resolve_organism, ORGANISM_ALIASES
from .organism_config import is_eukaryotic_organism
from .exceptions import (
    BioCompilerError, BiosecurityError, InvalidSequenceError, CertificateGenerationError,
    CertificateVerificationError, UnsupportedOrganismError, InvalidProteinError,
    OptimizationConstraintError,
)
from .provenance import (
    ProvenanceTracker,
    OptimizationProvenance,
    OptimizationRecord,
    generate_provenance_report,
)
from .decision_provenance import (
    ProvenanceStore,
    OptimizationDecisionTrail,
)

logger = logging.getLogger(__name__)

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
    # Input size limits
    "MAX_PROTEIN_LENGTH",
    "MAX_PROTEIN_SEQUENCE_LENGTH",
    "MAX_BATCH_SIZE",
    "MAX_REQUEST_SIZE",
    "MAX_DNA_LENGTH",
    "MAX_DNA_SEQUENCE_LENGTH",
    "OPTIMIZE_TIMEOUT_S",
    # Batch optimization (programmatic API)
    "batch_optimize",
    # Pydantic input models
    "SequenceInput",
    "ProteinInput",
    "CertificateInput",
    "ExportFastaInput",
    "ExportGenbankInput",
    "ScanInput",
    # Pydantic response models
    "TypeCheckResponse",
    "OptimizeResponse",
    "VerifyResponse",
    "ScanResponse",
    "OrganismResponse",
    "PredicateResponse",
    "HealthResponse",
    "InfoResponse",
    # Task 1.6 response/input models
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
    # Protein analysis input models
    "StructurePredictInput",
    "QualityAssessInput",
    "StabilityInput",
    "MutationScanInput",
    "SolubilityInput",
    "ImmunogenicityInput",
    "DeimmunizeInput",
    "FullAssessmentInput",
    # Protein analysis response models
    "StructurePredictResponse",
    "QualityAssessResponse",
    "StabilityResponse",
    "MutationScanResponse",
    "SolubilityResponse",
    "ImmunogenicityResponse",
    "DeimmunizeResponse",
    "FullAssessmentResponse",
    # Benchmark / Validation / WhatIf
    # (Removed BenchmarkInput, BenchmarkResponse, ValidateCAIInput,
    #  ValidateCAIResponse, ValidateMaxEntScanResponse, WhatIfInput,
    #  WhatIfResponse — not yet implemented; remove from __all__ to
    #  prevent AttributeError on `from biocompiler.api import *`)
    # Organism domain resolution
    "resolve_organism_domain",
    # Auth mode
    "set_no_auth_flag",
    "get_auth_mode",
    "get_configured_api_keys",
    "is_auth_enabled",
    # Provenance
    "ProvenanceResponse",
    "ProvenanceExplainResponse",
    "ProvenanceReportResponse",
    # Persistent provenance store
    "ProvenanceStore",
    "_provenance_store",
    # Persistent rate limiter
    "_rate_limiter",
]

# ─── API Key Authentication ─────────────────────────────────────────

API_KEY_NAME = "X-API-Key"
_api_key_header = APIKeyHeader(name=API_KEY_NAME, auto_error=False)

# ── Auth mode: "required" (default), "optional", "disabled" ──────────
# BIOCOMPILER_AUTH_MODE controls how authentication is enforced:
#   "required"  — unauthenticated requests get 401 (default, safe)
#   "optional"  — auth is checked but unauthenticated requests are
#                  allowed with a warning response header
#   "disabled"  — auth is completely disabled (DANGEROUS for production)
#
# BIOCOMPILER_API_KEY can be:
#   - A specific key string (recommended for production)
#   - "disabled" to explicitly disable auth (with warning)
#   - Not set → auto-generate a random key and persist to
#     ~/.biocompiler/api_key for reuse across restarts
#
# BIOCOMPILER_API_KEYS supports comma-separated keys for rotation.
# --no-auth CLI flag (for `biocompiler serve`) sets auth mode to "disabled".

_AUTH_MODE = os.environ.get("BIOCOMPILER_AUTH_MODE", "required").lower()
if _AUTH_MODE not in ("required", "optional", "disabled"):
    logger.warning(
        "Invalid BIOCOMPILER_AUTH_MODE=%r, falling back to 'required'",
        _AUTH_MODE,
    )
    _AUTH_MODE = "required"

# Global flag set by --no-auth CLI flag (see cli.py)
_NO_AUTH_CLI_FLAG = False


def set_no_auth_flag() -> None:
    """Enable --no-auth mode from the CLI. Emits a warning."""
    global _NO_AUTH_CLI_FLAG, _AUTH_MODE
    _NO_AUTH_CLI_FLAG = True
    _AUTH_MODE = "disabled"
    logger.warning(
        "⚠ --no-auth flag: API authentication is DISABLED. "
        "Do not use in production!"
    )


def _get_api_key_file_path() -> Path:
    """Return the path to the persisted API key file."""
    return Path.home() / ".biocompiler" / "api_key"


def _generate_and_persist_api_key() -> str:
    """Generate a random API key, persist it, and return it.

    The key is saved to ~/.biocompiler/api_key for reuse across restarts.
    On first startup, the key is printed to the console as a one-time message.
    """
    key_file = _get_api_key_file_path()
    is_new = not key_file.exists()

    if key_file.exists():
        try:
            existing_key = key_file.read_text().strip()
            if existing_key:
                if is_new:
                    # Shouldn't happen, but handle gracefully
                    pass
                return existing_key
        except (OSError, UnicodeDecodeError):
            logger.warning("Could not read API key from %s, generating new one", key_file)

    # Generate a new 32-byte hex key (64 chars)
    new_key = secrets.token_hex(32)

    # Persist to disk
    try:
        key_file.parent.mkdir(parents=True, exist_ok=True)
        key_file.write_text(new_key)
        # Set restrictive permissions (owner read/write only)
        key_file.chmod(0o600)
    except OSError as exc:
        logger.warning("Could not persist API key to %s: %s", key_file, exc)

    # One-time message to console
    print(
        f"\nGenerated API key: {new_key} (save this!)\n"
        f"Key saved to {key_file} for reuse across restarts.\n"
    )
    logger.info("Generated and persisted API key to %s", key_file)

    return new_key


# ── Resolve configured API keys ────────────────────────────────────
_CONFIGURED_API_KEYS: set[str] = set()

_single_key = os.environ.get("BIOCOMPILER_API_KEY", "")
_multi_keys = os.environ.get("BIOCOMPILER_API_KEYS", "")

if _multi_keys:
    _CONFIGURED_API_KEYS = {k.strip() for k in _multi_keys.split(",") if k.strip()}
elif _single_key == "disabled":
    # Explicit BIOCOMPILER_API_KEY=disabled → disable auth with warning
    logger.warning(
        "⚠ BIOCOMPILER_API_KEY=disabled: API authentication is DISABLED. "
        "This is dangerous for a tool that designs DNA. "
        "Use only in isolated development environments."
    )
    _AUTH_MODE = "disabled"
elif _single_key:
    _CONFIGURED_API_KEYS = {_single_key}
else:
    # No BIOCOMPILER_API_KEY set → auto-generate and persist a random key
    # This is the secure-by-default behavior
    auto_key = _generate_and_persist_api_key()
    _CONFIGURED_API_KEYS = {auto_key}

# ─── Rate Limiting ──────────────────────────────────────────────────

RATE_LIMIT_RPM = int(os.environ.get("BIOCOMPILER_RATE_LIMIT", "60"))  # requests per minute

from .rate_limiter import PersistentRateLimiter

_rate_limiter = PersistentRateLimiter(
    db_path=os.environ.get("BIOCOMPILER_RATE_LIMIT_DB", "~/.biocompiler/rate_limits.db"),
    max_requests=RATE_LIMIT_RPM,
    window_seconds=60,
)


def _check_rate_limit(client_id: str) -> None:
    """Enforce per-client rate limiting (sliding window, SQLite-backed)."""
    allowed, remaining = _rate_limiter.check(client_id)
    if not allowed:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded: {RATE_LIMIT_RPM} requests/minute. "
                   f"Retry after {_rate_limiter.window_seconds} seconds.",
        )
    _rate_limiter.record(client_id)


def _check_batch_rate_limit(client_id: str, item_count: int) -> None:
    """
    Enforce per-client rate limiting for batch requests.

    Each item in the batch consumes one rate-limit unit. We pre-check
    whether the entire batch can be accommodated, then record all units.
    This prevents partial batch execution when rate-limited mid-way.
    """
    allowed, remaining = _rate_limiter.check(client_id)
    if item_count > remaining:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded: batch has {item_count} items but only "
                   f"{remaining} rate-limit units remaining this minute. "
                   f"Limit: {RATE_LIMIT_RPM} requests/minute.",
        )
    _rate_limiter.record_batch(client_id, item_count)


def _is_auth_enabled() -> bool:
    """Return True if authentication is currently enabled."""
    return _AUTH_MODE != "disabled" and bool(_CONFIGURED_API_KEYS)


def get_auth_mode() -> str:
    """Return the current auth mode string: 'required', 'optional', or 'disabled'."""
    return _AUTH_MODE


def get_configured_api_keys() -> set[str]:
    """Return a copy of the configured API keys set (for testing/inspection)."""
    return set(_CONFIGURED_API_KEYS)


def is_auth_enabled() -> bool:
    """Public API: Return True if authentication is currently enabled."""
    return _is_auth_enabled()


async def verify_api_key(api_key: str = Security(_api_key_header)) -> str:
    """
    Verify the API key based on the current auth mode.

    Auth modes (controlled by BIOCOMPILER_AUTH_MODE):
    - "required" (default): Unauthenticated requests get HTTP 401.
    - "optional": Auth is checked; unauthenticated requests are allowed
      but receive a warning header (added by middleware).
    - "disabled": Auth is completely disabled (DANGEROUS).

    API key resolution:
    - BIOCOMPILER_API_KEYS (comma-separated) for key rotation
    - BIOCOMPILER_API_KEY for a single key
    - BIOCOMPILER_API_KEY=disabled explicitly disables auth
    - No env var → auto-generated key persisted to ~/.biocompiler/api_key
    """
    # Disabled mode: no auth required
    if _AUTH_MODE == "disabled":
        return "anonymous"

    # No keys configured (shouldn't happen with auto-generation, but guard)
    if not _CONFIGURED_API_KEYS:
        if _AUTH_MODE == "required":
            # This is an unusual state — no keys but auth required.
            # Allow through (the key auto-generation should have handled this).
            return "anonymous"
        return "anonymous"  # optional mode, no keys

    # Valid key provided
    if api_key is not None and api_key in _CONFIGURED_API_KEYS:
        return api_key

    # No key or invalid key
    if _AUTH_MODE == "required":
        raise HTTPException(
            status_code=401,
            detail="Invalid or missing API key. Set X-API-Key header.",
        )

    # Optional mode: allow through (middleware adds warning header)
    return "anonymous"


# ─── Pydantic Models ──────────────────────────────────────────────

class SequenceInput(BaseModel):
    """DNA sequence input."""
    sequence: str = Field(..., description="DNA sequence (ACGTN)")
    organism: str = Field(
        "Homo_sapiens",
        description=(
            "Target organism. Accepts canonical names (e.g. 'Escherichia_coli'), "
            "short keys ('ecoli', 'human'), abbreviated binomials ('E_coli', "
            "'h_sapiens'), or display names ('E. coli')."
        ),
    )
    species: Optional[str] = Field(
        None,
        description=(
            "Alias for organism (deprecated). Accepts the same values. "
            "If both species and organism are provided, species takes precedence."
        ),
    )
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
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        resolved = resolve_organism(v, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {v} (resolved to {resolved!r}). "
                f"Supported: {SUPPORTED_ORGANISMS}"
            )
        return resolved

    @field_validator("species")
    @classmethod
    def validate_species(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported species: {v} (resolved to {resolved!r}). "
                    f"Supported: {list(ORGANISM_ALIASES.keys())}"
                )
        return v


class ProteinInput(BaseModel):
    """Protein sequence input for optimization."""
    protein: str = Field(..., description="Target protein sequence (single-letter codes)")
    organism: str = Field(
        "Homo_sapiens",
        description=(
            "Target organism. Accepts canonical names (e.g. 'Escherichia_coli'), "
            "short keys ('ecoli', 'human'), abbreviated binomials ('E_coli', "
            "'h_sapiens'), or display names ('E. coli')."
        ),
    )
    species: Optional[str] = Field(
        None,
        description=(
            "Alias for organism (deprecated). Accepts the same values. "
            "If both species and organism are provided, species takes precedence."
        ),
    )
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.2, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to avoid")
    cryptic_splice_threshold: float = Field(3.0, description="Cryptic splice site threshold")
    track_provenance: bool = Field(True, description="Track provenance for this optimization (default: enabled)")
    organism_domain: str = Field(
        "auto",
        description=(
            "Organism domain for constraint selection. "
            "'auto' detects from organism name, "
            "'eukaryote' forces eukaryotic constraints (splice sites, CpG islands), "
            "'prokaryote' skips eukaryote-specific constraints."
        ),
    )
    source_organism: Optional[str] = Field(
        None,
        description=(
            "Organism the protein originates from. Used by immunogenicity predicates "
            "to determine whether the protein is 'self' (source matches host) or "
            "foreign. Accepts the same aliases as organism (e.g. 'e_coli', 'human'). "
            "If None, the protein is assumed to be from the host organism (self)."
        ),
    )
    therapeutic: bool = Field(
        False,
        description=(
            "Whether the protein is intended for therapeutic use. "
            "Therapeutic proteins have stricter immunogenicity thresholds "
            "because immune responses can compromise drug efficacy."
        ),
    )
    self_protein: Optional[bool] = Field(
        None,
        description=(
            "Explicit override for self-protein status. If True, the protein is "
            "treated as a self-protein (auto-PASS for immunogenicity). If False, "
            "it is treated as foreign. If None (default), self-status is "
            "auto-detected from source_organism vs organism."
        ),
    )
    strict_mode: bool = Field(
        True,
        description=(
            "If True (default), refuse to return sequences with failed predicates — "
            "the endpoint returns HTTP 422 instead. Set to False to allow partial "
            "results that have some unsatisfied constraints."
        ),
    )

    @field_validator("organism_domain")
    @classmethod
    def validate_organism_domain(cls, v: str) -> str:
        v = v.lower()
        if v not in ("auto", "eukaryote", "prokaryote"):
            raise ValueError(
                f"Invalid organism_domain: {v!r}. "
                "Must be one of: auto, eukaryote, prokaryote"
            )
        return v

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Protein sequence must not be empty.")
        v = v.upper()
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(v) - valid_aas
        if invalid:
            raise ValueError(
                f"Invalid amino acids: {sorted(invalid)}. "
                f"Allowed: {sorted(valid_aas)}"
            )
        if len(v) > MAX_PROTEIN_LENGTH:
            raise ValueError(
                f"Protein sequence too long ({len(v)} aa). "
                f"Maximum: {MAX_PROTEIN_LENGTH} aa."
            )
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        resolved = resolve_organism(v, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {v} (resolved to {resolved!r}). "
                f"Supported: {sorted(set(SUPPORTED_ORGANISMS))}"
            )
        return resolved

    @field_validator("species")
    @classmethod
    def validate_species(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported species: {v} (resolved to {resolved!r}). "
                    f"Supported: {sorted(set(ORGANISM_ALIASES.keys()))}"
                )
        return v

    @field_validator("source_organism")
    @classmethod
    def validate_source_organism(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported source_organism: {v} (resolved to {resolved!r}). "
                    f"Supported: {sorted(set(ORGANISM_ALIASES.keys()))}"
                )
            return resolved
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

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence must not be empty")
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v


class ExportGenbankInput(BaseModel):
    """GenBank export input."""
    sequence: str = Field(..., description="DNA sequence")
    locus_name: str = Field("BIOCOMPILER", description="LOCUS name (max 16 chars)")
    definition: str = Field("BioCompiler designed sequence", description="DEFINITION line")
    organism: str = Field("Homo_sapiens", description="Source organism")
    gene_name: Optional[str] = Field(None, description="Gene name")
    exon_boundaries: Optional[list[tuple[int, int]]] = Field(None)
    certificate: Optional[dict] = Field(None, description="Certificate dict to embed")

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence must not be empty")
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v


class ScanInput(BaseModel):
    """Sequence scan input."""
    sequence: str = Field(..., description="DNA sequence")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to scan for")
    find_orfs: bool = Field(False, description="Find open reading frames")

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence must not be empty")
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v


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
    provenance_id: Optional[str] = Field(None, description="Provenance trail ID if tracking was enabled")
    organism_domain: str = Field(
        "eukaryote",
        description=(
            "Resolved organism domain used for constraint selection. "
            "Either 'eukaryote' or 'prokaryote'."
        ),
    )
    source_organism: Optional[str] = Field(
        None,
        description="Resolved source organism (None = assumed self-protein).",
    )
    therapeutic: bool = Field(
        False,
        description="Whether therapeutic-mode immunogenicity thresholds were applied.",
    )
    self_protein: Optional[bool] = Field(
        None,
        description="Self-protein status used for immunogenicity predicates (None = auto-detected).",
    )


class VerifyResponse(BaseModel):
    status: str
    failure_reasons: list[str]


class ScanResponse(BaseModel):
    sequence_length: int
    tokens: list[dict]
    orfs: Optional[list[dict]] = None


class OrganismResponse(BaseModel):
    organisms: list[dict]


class InfoResponse(BaseModel):
    """Response from the /info endpoint."""
    max_protein_length: int = Field(..., description="Maximum protein sequence length in amino acids")
    max_dna_length: int = Field(..., description="Maximum DNA sequence length in bases")
    max_batch_size: int = Field(..., description="Maximum number of items per batch request")
    max_request_size: int = Field(..., description="Maximum request body size in bytes")
    optimize_timeout_s: int = Field(..., description="Optimization timeout in seconds")
    supported_organisms: list[str] = Field(..., description="List of supported organism names")
    api_version: str = Field(..., description="API version")
    safety_version: str = Field(..., description="Safety/biosecurity version")


class PredicateResponse(BaseModel):
    predicates: list[str]


class HealthResponse(BaseModel):
    status: str
    version: str
    timestamp: str
    auth_enabled: bool = False
    rate_limit_rpm: int = 60
    cors_origins: list[str] = []
    cors_allow_credentials: bool = False


# ─── Additional Response/Input Models (Task 1.6) ────────────────

class EnzymeListResponse(BaseModel):
    """Response model for /enzymes endpoint."""
    enzymes: dict[str, str] = Field(..., description="Enzyme name to recognition site mapping")


class ProvenanceRecordSummary(BaseModel):
    """Summary of a single provenance record."""
    gene_name: str = Field(..., description="Gene name")
    organism: str = Field(..., description="Organism")


class ProvenanceDetailResponse(BaseModel):
    """Response model for /provenance/{id} endpoint."""
    id: str = Field(..., description="Record ID")
    trail: dict = Field(..., description="Full provenance trail")


class ProvenanceListResponse(BaseModel):
    """Response model for /provenance endpoint."""
    count: int = Field(..., description="Number of records")
    records: list[ProvenanceRecordSummary] = Field(default_factory=list, description="Record summaries")


class ExportFastaResponse(BaseModel):
    """Response model for /export/fasta endpoint."""
    format: str = Field("fasta", description="Export format")
    content: str = Field(..., description="FASTA content")


class ExportGenbankResponse(BaseModel):
    """Response model for /export/genbank endpoint."""
    format: str = Field("genbank", description="Export format")
    content: str = Field(..., description="GenBank content")


class ExportSbol3Input(BaseModel):
    """Input model for /export/sbol3 endpoint."""
    sequence: str = Field(..., description="DNA sequence")
    organism: str = Field("Homo_sapiens", description="Source organism")
    gene_name: str = Field("optimized_gene", description="Gene name")
    format: str = Field("sbol3", description="SBOL format: 'sbol3' or 'sbol3json'")

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence must not be empty")
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        v = v.upper()
        invalid = set(v) - set("ACGT")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        resolved = resolve_organism(v, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {v} (resolved to {resolved!r}). "
                f"Supported: {sorted(set(SUPPORTED_ORGANISMS))}"
            )
        return resolved

    @field_validator("format")
    @classmethod
    def validate_format(cls, v: str) -> str:
        v = v.lower()
        if v not in ("sbol3", "sbol3json"):
            raise ValueError(f"Unsupported SBOL format: {v}. Supported: sbol3, sbol3json")
        return v


class ExportSbol3Response(BaseModel):
    """Response model for /export/sbol3 endpoint."""
    format: str = Field(..., description="Export format")
    content: str = Field(..., description="SBOL3 content (XML or JSON)")


class DatasetValidationResult(BaseModel):
    """Result of a single dataset validation test."""
    dataset: str = Field(..., description="Dataset name")
    test_name: str = Field(..., description="Test name")
    passed: bool = Field(..., description="Whether the test passed")


class DatasetValidationResponse(BaseModel):
    """Response model for /validate-datasets endpoint."""
    total_tests: int = Field(..., description="Total number of tests")
    passed: int = Field(..., description="Number of tests passed")
    failed: int = Field(..., description="Number of tests failed")
    results: list[DatasetValidationResult] = Field(default_factory=list, description="Per-test results")


# ─── Batch Pydantic Models ────────────────────────────────────────

BATCH_CHECK_MAX = 50
BATCH_OPTIMIZE_MAX = 20
BATCH_EXPORT_MAX = 50
BATCH_ITEM_TIMEOUT_S = int(os.environ.get("BIOCOMPILER_BATCH_ITEM_TIMEOUT", "30"))


class BatchCheckItem(BaseModel):
    """Single item in a batch type-check request."""
    sequence: str = Field(..., description="DNA sequence (ACGTN)")
    organism: str = Field(
        "Homo_sapiens",
        description=(
            "Target organism. Accepts canonical names, short keys, "
            "abbreviated binomials, or display names."
        ),
    )
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
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        resolved = resolve_organism(v, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {v} (resolved to {resolved!r}). "
                f"Supported: {SUPPORTED_ORGANISMS}"
            )
        return resolved


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
    organism: str = Field(
        "Homo_sapiens",
        description=(
            "Target organism. Accepts canonical names, short keys, "
            "abbreviated binomials, or display names."
        ),
    )
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.2, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to avoid")
    cryptic_splice_threshold: float = Field(3.0, description="Cryptic splice site threshold")
    organism_domain: str = Field(
        "auto",
        description=(
            "Organism domain for constraint selection. "
            "'auto' detects from organism name, "
            "'eukaryote' forces eukaryotic constraints, "
            "'prokaryote' skips eukaryote-specific constraints."
        ),
    )
    source_organism: Optional[str] = Field(
        None,
        description=(
            "Organism the protein originates from. Accepts the same aliases as organism. "
            "If None, the protein is assumed to be from the host organism (self)."
        ),
    )
    therapeutic: bool = Field(
        False,
        description="Whether the protein is intended for therapeutic use (stricter immunogenicity thresholds).",
    )
    self_protein: Optional[bool] = Field(
        None,
        description="Explicit override for self-protein status. None = auto-detect from source_organism.",
    )

    @field_validator("organism_domain")
    @classmethod
    def validate_organism_domain(cls, v: str) -> str:
        v = v.lower()
        if v not in ("auto", "eukaryote", "prokaryote"):
            raise ValueError(
                f"Invalid organism_domain: {v!r}. "
                "Must be one of: auto, eukaryote, prokaryote"
            )
        return v

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Protein sequence must not be empty.")
        v = v.upper()
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(v) - valid_aas
        if invalid:
            raise ValueError(
                f"Invalid amino acids: {sorted(invalid)}. "
                f"Allowed: {sorted(valid_aas)}"
            )
        if len(v) > MAX_PROTEIN_LENGTH:
            raise ValueError(
                f"Protein sequence too long ({len(v)} aa). "
                f"Maximum: {MAX_PROTEIN_LENGTH} aa."
            )
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        resolved = resolve_organism(v, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {v} (resolved to {resolved!r}). "
                f"Supported: {sorted(set(SUPPORTED_ORGANISMS))}"
            )
        return resolved

    @field_validator("source_organism")
    @classmethod
    def validate_source_organism(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported source_organism: {v} (resolved to {resolved!r}). "
                    f"Supported: {sorted(set(ORGANISM_ALIASES.keys()))}"
                )
            return resolved
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


class FastBatchOptimizeInput(BaseModel):
    """Input for fast batch optimize endpoint (shared organism/parameters).

    All proteins are optimized with the same organism and constraint
    parameters, allowing the backend to reuse a single HybridOptimizer
    instance.  This is significantly faster than the per-item batch
    endpoint when all proteins target the same organism.
    """
    proteins: list[str] = Field(
        ..., description=f"List of protein sequences to optimize (max {BATCH_OPTIMIZE_MAX})"
    )
    organism: str = Field(
        "Homo_sapiens",
        description=(
            "Target organism shared by all proteins. Accepts canonical names, "
            "short keys, abbreviated binomials, or display names."
        ),
    )
    gc_lo: float = Field(0.30, description="Minimum GC content")
    gc_hi: float = Field(0.70, description="Maximum GC content")
    cai_threshold: float = Field(0.2, description="Minimum CAI threshold")
    enzymes: Optional[list[str]] = Field(None, description="Restriction enzymes to avoid")
    organism_domain: str = Field(
        "auto",
        description=(
            "Organism domain for constraint selection. "
            "'auto' detects from organism name, "
            "'eukaryote' forces eukaryotic constraints, "
            "'prokaryote' skips eukaryote-specific constraints."
        ),
    )
    source_organism: Optional[str] = Field(
        None,
        description=(
            "Organism the proteins originate from. Accepts the same aliases as organism. "
            "If None, proteins are assumed to be from the host organism (self)."
        ),
    )
    therapeutic: bool = Field(
        False,
        description="Whether the proteins are intended for therapeutic use (stricter immunogenicity thresholds).",
    )
    self_protein: Optional[bool] = Field(
        None,
        description="Explicit override for self-protein status. None = auto-detect from source_organism.",
    )

    @field_validator("organism_domain")
    @classmethod
    def validate_organism_domain(cls, v: str) -> str:
        v = v.lower()
        if v not in ("auto", "eukaryote", "prokaryote"):
            raise ValueError(
                f"Invalid organism_domain: {v!r}. "
                "Must be one of: auto, eukaryote, prokaryote"
            )
        return v

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        resolved = resolve_organism(v, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {v} (resolved to {resolved!r}). "
                f"Supported: {SUPPORTED_ORGANISMS}"
            )
        return resolved

    @field_validator("source_organism")
    @classmethod
    def validate_source_organism(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported source_organism: {v} (resolved to {resolved!r}). "
                    f"Supported: {list(ORGANISM_ALIASES.keys())}"
                )
            return resolved
        return v

    @field_validator("proteins")
    @classmethod
    def validate_proteins(cls, v: list[str]) -> list[str]:
        if not v:
            raise ValueError("Proteins list must not be empty.")
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        for i, protein in enumerate(v):
            if not protein or not protein.strip():
                raise ValueError(
                    f"Protein at index {i} must not be empty."
                )
            protein_upper = protein.upper()
            invalid = set(protein_upper) - valid_aas
            if invalid:
                raise ValueError(
                    f"Invalid amino acids in protein at index {i}: {sorted(invalid)}. "
                    f"Allowed: {sorted(valid_aas)}"
                )
            if len(protein_upper) > MAX_PROTEIN_LENGTH:
                raise ValueError(
                    f"Protein at index {i} too long ({len(protein_upper)} aa). "
                    f"Maximum: {MAX_PROTEIN_LENGTH} aa."
                )
        return [p.upper() for p in v]


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

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence must not be empty")
        if len(v) > MAX_DNA_LENGTH:
            raise ValueError(
                f"DNA sequence too long ({len(v)} bases). "
                f"Maximum: {MAX_DNA_LENGTH} bases."
            )
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v

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


# ─── Organism Domain Resolution ────────────────────────────────────


def resolve_organism_domain(organism: str, organism_domain: str = "auto") -> str:
    """Resolve the effective organism domain from user input.

    When *organism_domain* is ``"auto"``, the domain is detected from
    the organism name using :func:`is_eukaryotic_organism`.  When
    explicitly set to ``"eukaryote"`` or ``"prokaryote"``, that value
    is used regardless of the organism name — allowing users to
    override auto-detection when needed.

    Args:
        organism: Organism name (e.g. ``"Escherichia_coli"``).
        organism_domain: One of ``"auto"``, ``"eukaryote"``, or
            ``"prokaryote"``.  Defaults to ``"auto"``.

    Returns:
        The resolved domain string: ``"eukaryote"`` or ``"prokaryote"``.

    Raises:
        ValueError: If *organism_domain* is not one of the valid choices.
    """
    if organism_domain not in ("auto", "eukaryote", "prokaryote"):
        raise ValueError(
            f"Invalid organism_domain: {organism_domain!r}. "
            "Must be one of: auto, eukaryote, prokaryote"
        )

    if organism_domain == "auto":
        return "eukaryote" if is_eukaryotic_organism(organism) else "prokaryote"

    # Explicit override
    return organism_domain


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
        boundaries=exon_boundaries,
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
        except CertificateGenerationError as exc:
            logger.warning("Certificate generation failed for batch item: %s", exc)

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
    resolved_domain = resolve_organism_domain(item.organism, item.organism_domain)

    result = optimize_sequence(
        target_protein=item.protein,
        organism=item.organism,
        gc_lo=item.gc_lo,
        gc_hi=item.gc_hi,
        cai_threshold=item.cai_threshold,
        restriction_sites=item.enzymes,
        cryptic_splice_threshold=item.cryptic_splice_threshold,
        organism_domain=resolved_domain,
        source_organism=item.source_organism,
        therapeutic=item.therapeutic,
        self_protein=item.self_protein,
    )

    return {
        "sequence": result.sequence,
        "protein": result.protein,
        "cai": result.cai,
        "gc_content": result.gc_content,
        "satisfied_predicates": result.satisfied_predicates,
        "failed_predicates": result.failed_predicates,
        "fallback_used": result.fallback_used,
        "organism_domain": resolved_domain,
        "source_organism": item.source_organism,
        "therapeutic": item.therapeutic,
        "self_protein": item.self_protein,
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


# ─── Provenance Pydantic Models ─────────────────────────────────────

class ProvenanceResponse(BaseModel):
    """Response from provenance retrieval."""
    id: str = Field(..., description="Provenance trail ID")
    seed: int = Field(..., description="RNG seed used")
    decision_count: int = Field(..., description="Number of decisions recorded")
    decisions: list[dict] = Field(default_factory=list, description="Decision records")
    optimization_records: list[dict] = Field(default_factory=list, description="Optimization summary records")


class ProvenanceExplainResponse(BaseModel):
    """Response from provenance explain at a specific position."""
    id: str = Field(..., description="Provenance trail ID")
    position: int = Field(..., description="Queried nucleotide position")
    decisions: list[dict] = Field(default_factory=list, description="Decisions at this position")
    explanation: str = Field(..., description="Human-readable explanation")


class ProvenanceReportResponse(BaseModel):
    """Response from provenance report generation."""
    id: str = Field(..., description="Provenance trail ID")
    format: str = Field(..., description="Report format: text, markdown, or json")
    report: str = Field(..., description="Generated report content")


# ─── Persistent Provenance Store ─────────────────────────────────────

_PROVENANCE_DIR = os.environ.get("BIOCOMPILER_PROVENANCE_DIR")
_provenance_store = ProvenanceStore(store_dir=_PROVENANCE_DIR)


def _store_provenance(trail: OptimizationDecisionTrail) -> str:
    """Persist a decision-level provenance trail and return its record ID.

    Uses the :class:`ProvenanceStore` from :mod:`decision_provenance`
    for file-based persistence so records survive process restarts.
    """
    record_id = _provenance_store.save(trail)
    return record_id


# ─── Protein Analysis Constants ─────────────────────────────────────

ESMFOLD_TIMEOUT_S = 120  # Default timeout for ESMFold structure prediction

# ─── Input Size Limits ────────────────────────────────────────────

MAX_PROTEIN_LENGTH = int(os.environ.get("BIOCOMPILER_MAX_PROTEIN_LENGTH", "10000"))  # aa
MAX_PROTEIN_SEQUENCE_LENGTH = MAX_PROTEIN_LENGTH  # backward-compatible alias
MAX_BATCH_SIZE = int(os.environ.get("BIOCOMPILER_MAX_BATCH_SIZE", "50"))  # sequences per batch
MAX_REQUEST_SIZE = int(os.environ.get("BIOCOMPILER_MAX_REQUEST_SIZE", str(10_000_000)))  # bytes
OPTIMIZE_TIMEOUT_S = int(os.environ.get("BIOCOMPILER_OPTIMIZE_TIMEOUT", "300"))  # seconds
MAX_DNA_LENGTH = int(os.environ.get("BIOCOMPILER_MAX_DNA_LENGTH", "100000"))  # bases
MAX_DNA_SEQUENCE_LENGTH = MAX_DNA_LENGTH  # backward-compatible alias

_PROTEIN_VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


# ─── Protein Analysis Helper Functions ──────────────────────────────

def validate_protein_input(protein: str) -> str | None:
    """
    Validate a protein sequence string.

    Returns an error message string if invalid, or None if valid.
    Checks: non-empty, only standard amino acid single-letter codes,
    length within MAX_PROTEIN_LENGTH.
    """
    if not protein:
        return "Protein sequence must not be empty."
    if not protein.strip():
        return "Protein sequence must not be empty."
    protein_upper = protein.upper()
    invalid = set(protein_upper) - _PROTEIN_VALID_AMINO_ACIDS
    if invalid:
        return f"Invalid amino acids in protein: {sorted(invalid)}. " \
               f"Allowed: {sorted(_PROTEIN_VALID_AMINO_ACIDS)}"
    if len(protein_upper) > MAX_PROTEIN_LENGTH:
        return f"Protein sequence too long ({len(protein_upper)} aa). Maximum: {MAX_PROTEIN_LENGTH} aa."
    return None


def validate_organism_input(organism: str) -> str | None:
    """
    Validate an organism string.

    Returns an error message string if invalid, or None if valid.
    Checks: non-empty, supported organism (including alias resolution).

    This function resolves organism aliases (e.g. 'ecoli', 'E. coli',
    'h_sapiens') to canonical names before validation, so that any
    accepted alias is also valid here.
    """
    if not organism:
        return "Organism must not be empty."
    try:
        resolved = resolve_organism(organism, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            return f"Unsupported organism: {organism} (resolved to {resolved!r}). Supported: {SUPPORTED_ORGANISMS}"
    except ImportError:
        # Fallback: accept any non-empty string if organisms module unavailable
        logger.debug("Organisms module unavailable; skipping organism validation")
    return None


# ─── Protein Analysis Pydantic Input Models ─────────────────────────

class StructurePredictInput(BaseModel):
    """Input for protein structure prediction via ESMFold."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    organism: str = Field(
        "Homo_sapiens", description="Target organism for codon context"
    )
    use_cache: bool = Field(
        True, description="Use cached structure if available"
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        err = validate_organism_input(v)
        if err:
            raise ValueError(err)
        return v


class QualityAssessInput(BaseModel):
    """Input for structure quality assessment from a PDB string."""
    pdb_string: str = Field(
        ..., description="PDB-format structure string"
    )

    @field_validator("pdb_string")
    @classmethod
    def validate_pdb(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("PDB string must not be empty.")
        if len(v) > 10_000_000:
            raise ValueError("PDB string too large (max 10 MB).")
        return v


class StabilityInput(BaseModel):
    """Input for protein stability analysis."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    organism: str = Field(
        "Homo_sapiens", description="Target organism"
    )
    pdb_string: str | None = Field(
        None, description="Optional PDB-format structure string"
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        err = validate_organism_input(v)
        if err:
            raise ValueError(err)
        return v


class MutationScanInput(BaseModel):
    """Input for scanning mutations for stability effects."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    positions: list[int] | None = Field(
        None, description="Specific positions to scan (1-indexed). None = all positions."
    )
    method: str = Field(
        "empirical", description="Scanning method: 'empirical' or 'statistical'"
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()

    @field_validator("method")
    @classmethod
    def validate_method(cls, v: str) -> str:
        v = v.lower()
        if v not in ("empirical", "statistical"):
            raise ValueError(f"Unsupported method: {v}. Supported: empirical, statistical")
        return v

    @field_validator("positions")
    @classmethod
    def validate_positions(cls, v: list[int] | None) -> list[int] | None:
        if v is not None:
            if any(p < 1 for p in v):
                raise ValueError("Positions must be 1-indexed (>= 1).")
        return v


class SolubilityInput(BaseModel):
    """Input for protein solubility analysis."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    pdb_string: str | None = Field(
        None, description="Optional PDB-format structure string"
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()


class ImmunogenicityInput(BaseModel):
    """Input for protein immunogenicity analysis."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    organism: str = Field(
        "Homo_sapiens", description="Host organism for MHC context"
    )
    mhc_alleles: list[str] | None = Field(
        None, description="Specific MHC alleles to test. None = common alleles for organism."
    )
    source_organism: Optional[str] = Field(
        None,
        description=(
            "Organism the protein originates from. Accepts the same aliases as organism. "
            "If None, the protein is assumed to be from the host organism (self)."
        ),
    )
    therapeutic: bool = Field(
        False,
        description="Whether the protein is intended for therapeutic use (stricter thresholds).",
    )
    self_protein: Optional[bool] = Field(
        None,
        description="Explicit override for self-protein status. None = auto-detect from source_organism.",
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        err = validate_organism_input(v)
        if err:
            raise ValueError(err)
        return v

    @field_validator("source_organism")
    @classmethod
    def validate_source_organism(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported source_organism: {v} (resolved to {resolved!r}). "
                    f"Supported: {list(ORGANISM_ALIASES.keys())}"
                )
            return resolved
        return v


class DeimmunizeInput(BaseModel):
    """Input for protein deimmunization."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    organism: str = Field(
        "Homo_sapiens", description="Host organism for MHC context"
    )
    target_score: float = Field(
        0.3, description="Target immunogenicity score (lower = less immunogenic)"
    )
    max_mutations: int = Field(
        10, description="Maximum number of mutations to introduce", ge=1, le=50
    )
    blosum62_min: int = Field(
        0, description="Minimum BLOSUM62 score for allowed substitutions", ge=-4, le=4
    )
    source_organism: Optional[str] = Field(
        None,
        description=(
            "Organism the protein originates from. Accepts the same aliases as organism. "
            "If None, the protein is assumed to be from the host organism (self)."
        ),
    )
    therapeutic: bool = Field(
        False,
        description="Whether the protein is intended for therapeutic use (stricter thresholds).",
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        err = validate_organism_input(v)
        if err:
            raise ValueError(err)
        return v

    @field_validator("target_score")
    @classmethod
    def validate_target_score(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError("target_score must be between 0.0 and 1.0.")
        return v

    @field_validator("source_organism")
    @classmethod
    def validate_source_organism(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            resolved = resolve_organism(v, strict=False)
            if resolved not in SUPPORTED_ORGANISMS:
                raise ValueError(
                    f"Unsupported source_organism: {v} (resolved to {resolved!r}). "
                    f"Supported: {list(ORGANISM_ALIASES.keys())}"
                )
            return resolved
        return v


class FullAssessmentInput(BaseModel):
    """Input for full protein assessment."""
    protein: str = Field(
        ..., description="Target protein sequence (single-letter codes)"
    )
    organism: str = Field(
        "Homo_sapiens", description="Target organism"
    )
    pdb_string: str | None = Field(
        None, description="Optional PDB-format structure string"
    )
    run_structure: bool = Field(
        True, description="Include structure prediction/quality in assessment"
    )
    run_stability: bool = Field(
        True, description="Include stability analysis in assessment"
    )
    run_solubility: bool = Field(
        True, description="Include solubility analysis in assessment"
    )
    run_immunogenicity: bool = Field(
        True, description="Include immunogenicity analysis in assessment"
    )

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        err = validate_protein_input(v)
        if err:
            raise ValueError(err)
        return v.upper()

    @field_validator("organism")
    @classmethod
    def validate_organism(cls, v: str) -> str:
        err = validate_organism_input(v)
        if err:
            raise ValueError(err)
        return v


# ─── Protein Analysis Pydantic Response Models ──────────────────────

class StructurePredictResponse(BaseModel):
    """Response from protein structure prediction."""
    pdb_string: str = Field(..., description="Predicted structure in PDB format")
    mean_plddt: float = Field(..., description="Mean pLDDT confidence score (0-100)")
    plddt_scores: list[float] = Field(
        ..., description="Per-residue pLDDT scores"
    )
    quality_class: str = Field(
        ..., description="Quality classification: very_high, high, medium, low"
    )
    execution_time_s: float = Field(
        ..., description="Wall-clock execution time in seconds"
    )


class QualityAssessResponse(BaseModel):
    """Response from structure quality assessment."""
    mean_plddt: float = Field(..., description="Mean pLDDT score (0-100)")
    ramachandran_favored: float = Field(
        ..., description="Fraction of residues in favored Ramachandran regions (0-1)"
    )
    clash_score: float = Field(
        ..., description="Clash score (lower is better)"
    )
    overall_quality: str = Field(
        ..., description="Overall quality: excellent, good, acceptable, poor"
    )
    verdict: str = Field(
        ..., description="Assessment verdict: PASS, WARN, or FAIL"
    )


class StabilityResponse(BaseModel):
    """Response from protein stability analysis."""
    stability_kcal: float = Field(
        ..., description="Predicted stability in kcal/mol (negative = stable)"
    )
    method: str = Field(
        ..., description="Method used: empirical, statistical, or mixed"
    )
    components: dict = Field(
        ..., description="Breakdown of stability components"
    )
    verdict: str = Field(
        ..., description="Stability verdict: STABLE, MARGINAL, or UNSTABLE"
    )


class MutationScanResponse(BaseModel):
    """Response from mutation scanning for stability."""
    mutations: list[dict] = Field(
        ..., description="List of mutation results with position, AA change, and ddG"
    )
    stabilizing_count: int = Field(
        ..., description="Number of stabilizing mutations found"
    )
    destabilizing_count: int = Field(
        ..., description="Number of destabilizing mutations found"
    )


class SolubilityResponse(BaseModel):
    """Response from protein solubility analysis."""
    intrinsic_score: float = Field(
        ..., description="Intrinsic solubility score (-3 to +3, higher = more soluble)"
    )
    overall_score: float = Field(
        ..., description="Overall solubility score including structural effects (-3 to +3)"
    )
    solubility_class: str = Field(
        ..., description="Solubility class: high, medium, low, very_low"
    )
    aggregation_prone_regions: list = Field(
        ..., description="List of aggregation-prone regions with positions and scores"
    )


class ImmunogenicityResponse(BaseModel):
    """Response from protein immunogenicity analysis."""
    overall_score: float = Field(
        ..., description="Overall immunogenicity score (0-1, higher = more immunogenic)"
    )
    immunogenicity_class: str = Field(
        ..., description="Class: low, moderate, high, very_high"
    )
    t_cell_epitopes: list[dict] = Field(
        ..., description="Predicted T-cell epitopes with positions and scores"
    )
    b_cell_epitopes: list[dict] = Field(
        ..., description="Predicted B-cell epitopes with positions and scores"
    )
    deimmunization_candidates: list[dict] = Field(
        ..., description="Candidate positions for deimmunization mutations"
    )


class DeimmunizeResponse(BaseModel):
    """Response from protein deimmunization."""
    optimized_protein: str = Field(
        ..., description="Deimmunized protein sequence"
    )
    mutations_applied: int = Field(
        ..., description="Number of mutations applied"
    )
    original_score: float = Field(
        ..., description="Original immunogenicity score"
    )
    optimized_score: float = Field(
        ..., description="Optimized immunogenicity score"
    )
    success: bool = Field(
        ..., description="Whether target score was achieved"
    )


class FullAssessmentResponse(BaseModel):
    """Response from full protein assessment."""
    structure_quality: dict | None = Field(
        None, description="Structure quality assessment results"
    )
    stability: dict | None = Field(
        None, description="Stability analysis results"
    )
    solubility: dict | None = Field(
        None, description="Solubility analysis results"
    )
    immunogenicity: dict | None = Field(
        None, description="Immunogenicity analysis results"
    )
    predicate_results: list[dict] = Field(
        ..., description="Predicate results for the optimized sequence"
    )
    overall_verdict: str = Field(
        ..., description="Overall assessment verdict: PASS, WARN, or FAIL"
    )
    recommendations: list[str] = Field(
        ..., description="Actionable recommendations based on assessment"
    )


# ─── Protein Analysis Router ────────────────────────────────────────

_protein_router = APIRouter(tags=["Protein Analysis"])


# ─── Structure Endpoints ────────────────────────────────────────────

@_protein_router.post(
    "/structure/predict",
    response_model=StructurePredictResponse,
    summary="Predict protein structure via ESMFold",
)
async def structure_predict(input_data: StructurePredictInput) -> StructurePredictResponse:
    """
    Predict protein 3D structure using ESMFold.

    Returns a PDB-format structure string with per-residue confidence
    scores (pLDDT). ESMFold is a language-model-based structure predictor
    that does not require MSA search.

    Rate limiting: This is a compute-intensive endpoint.
    Timeout: 120 seconds by default (BIOCOMPILER_ESMFOLD_TIMEOUT env var).
    Returns 503 if ESMFold is not available.
    """
    try:
        from .esmfold import predict_structure, is_esmfold_available
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="ESMFold module not available. Install biocompiler[esmfold].",
        )

    if not is_esmfold_available():
        raise HTTPException(
            status_code=503,
            detail="ESMFold model is not available. Ensure torch and esm are installed.",
        )

    t0 = time.monotonic()
    try:
        result = predict_structure(
            protein=input_data.protein,
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=f"ESMFold prediction timed out after {ESMFOLD_TIMEOUT_S}s.",
        )
    except Exception as e:
        logger.exception("ESMFold prediction failed")
        raise HTTPException(
            status_code=500,
            detail=f"Structure prediction failed: {e}",
        )

    elapsed = time.monotonic() - t0

    # result is an ESMFoldResult dataclass, not a dict
    mean_plddt = float(getattr(result, "mean_plddt", 0.0))
    if mean_plddt >= 90:
        quality_class = "very_high"
    elif mean_plddt >= 70:
        quality_class = "high"
    elif mean_plddt >= 50:
        quality_class = "medium"
    else:
        quality_class = "low"

    return StructurePredictResponse(
        pdb_string=result.pdb_string,
        mean_plddt=mean_plddt,
        plddt_scores=getattr(result, "plddt_scores", []),
        quality_class=quality_class,
        execution_time_s=round(elapsed, 3),
    )


@_protein_router.post(
    "/structure/quality",
    response_model=QualityAssessResponse,
    summary="Assess structure quality from PDB",
)
async def structure_quality(input_data: QualityAssessInput) -> QualityAssessResponse:
    """
    Assess the quality of a protein structure from a PDB string.

    Evaluates pLDDT scores, Ramachandran outliers, clash score,
    and returns an overall quality assessment with a verdict.

    Returns 400 if PDB string is malformed.
    """
    try:
        from .structure.quality import compute_structure_quality
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="Structure quality module not available.",
        )

    try:
        result = compute_structure_quality(pdb_string=input_data.pdb_string)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid PDB input: {e}")
    except Exception as e:
        logger.exception("Structure quality assessment failed")
        raise HTTPException(
            status_code=500,
            detail=f"Quality assessment failed: {e}",
        )

    return QualityAssessResponse(
        mean_plddt=result["mean_plddt"],
        ramachandran_favored=result["ramachandran_favored"],
        clash_score=result["clash_score"],
        overall_quality=result["overall_quality"],
        verdict=result["verdict"],
    )


# ─── Structure Batch Endpoint ─────────────────────────────────────────

BATCH_PROTEIN_MAX = 20

class BatchStructureItem(StructurePredictInput):
    """Single item in a batch structure prediction request."""
    pass


class BatchStructureInput(BaseModel):
    """Input for batch structure prediction endpoint."""
    proteins: list[BatchStructureItem] = Field(
        ..., description=f"List of proteins for structure prediction (max {BATCH_PROTEIN_MAX})"
    )


class BatchStructureResultItem(BaseModel):
    """Result for a single structure prediction item."""
    pdb_string: Optional[str] = Field(None, description="Predicted structure in PDB format")
    mean_plddt: float = Field(0.0, description="Mean pLDDT confidence score (0-100)")
    quality_class: Optional[str] = Field(None, description="Quality classification")
    error: Optional[str] = Field(None, description="Error message if prediction failed")


class BatchStructureResponse(BaseModel):
    """Response for batch structure prediction endpoint."""
    results: list[BatchStructureResultItem] = Field(..., description="Per-item structure prediction results")
    summary: dict = Field(..., description="Aggregate summary")


@_protein_router.post(
    "/structure/batch",
    response_model=BatchStructureResponse,
    summary="Batch protein structure prediction",
)
async def structure_batch(input_data: BatchStructureInput) -> BatchStructureResponse:
    """
    Predict protein structures for multiple proteins in one request.

    Each protein is processed independently using ESMFold. One failure
    does not affect other items.

    Maximum batch size: 20 proteins.
    Each item consumes one rate-limit unit.
    Returns 503 if ESMFold is not available.
    """
    try:
        from .esmfold import predict_structure, is_esmfold_available
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="ESMFold module not available. Install biocompiler[esmfold].",
        )

    if not is_esmfold_available():
        raise HTTPException(
            status_code=503,
            detail="ESMFold model is not available. Ensure torch and esm are installed.",
        )

    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX} proteins",
        )

    client_id = "batch_structure"
    _check_batch_rate_limit(client_id, n)

    results: list[BatchStructureResultItem] = []
    total = n
    high_count = 0
    medium_count = 0
    low_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(_structure_batch_single, item),
                timeout=ESMFOLD_TIMEOUT_S,
            )
            results.append(result)
            if result.error is None:
                qc = result.quality_class
                if qc in ("very_high", "high"):
                    high_count += 1
                elif qc == "medium":
                    medium_count += 1
                else:
                    low_count += 1
            else:
                error_count += 1
        except asyncio.TimeoutError:
            results.append(BatchStructureResultItem(
                error=f"Timeout after {ESMFOLD_TIMEOUT_S}s",
            ))
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


def _structure_batch_single(item: StructurePredictInput) -> BatchStructureResultItem:
    """Process a single structure prediction item for batch processing."""
    from .esmfold import predict_structure

    result = predict_structure(
        protein=item.protein,
    )

    if not result.success:
        return BatchStructureResultItem(error=result.error or "Prediction failed")

    mean_plddt = float(getattr(result, "mean_plddt", 0.0))
    if mean_plddt >= 90:
        quality_class = "very_high"
    elif mean_plddt >= 70:
        quality_class = "high"
    elif mean_plddt >= 50:
        quality_class = "medium"
    else:
        quality_class = "low"

    return BatchStructureResultItem(
        pdb_string=result.pdb_string,
        mean_plddt=mean_plddt,
        quality_class=quality_class,
    )


# ─── Stability Endpoints ────────────────────────────────────────────

@_protein_router.post(
    "/stability/analyze",
    response_model=StabilityResponse,
    summary="Analyze protein stability",
)
async def stability_analyze(input_data: StabilityInput) -> StabilityResponse:
    """
    Analyze protein stability using FoldX (empirical) and/or
    statistical potentials.

    If a PDB structure is provided, FoldX-based energy computation
    is used. Otherwise, a sequence-based statistical method is applied.

    Returns 503 if FoldX is not installed and no statistical fallback is
    available.
    """
    try:
        from .foldx import empirical_stability
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="FoldX stability module not available. Install biocompiler[foldx].",
        )

    try:
        result = empirical_stability(
            protein=input_data.protein,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Stability analysis failed")
        raise HTTPException(
            status_code=500,
            detail=f"Stability analysis failed: {e}",
        )

    # result is a FoldXResult dataclass — build components dict and verdict
    components = {
        key: getattr(result, key)
        for key in (
            "backbone_hbond", "sidechain_hbond", "van_der_waals",
            "electrostatics", "solvation", "van_der_waals_clashes",
            "entropy_sidechain", "entropy_mainchain", "torsional_clash",
            "backbone_clash", "helix_dipole", "disulfide",
            "electrostatic_kon", "partial_covalent", "energy_ionisation",
        )
        if getattr(result, key, None) is not None
    }
    stability_kcal = result.stability_kcal
    if stability_kcal < -5.0:
        verdict = "STABLE"
    elif stability_kcal < 0.0:
        verdict = "MARGINAL"
    else:
        verdict = "UNSTABLE"

    return StabilityResponse(
        stability_kcal=stability_kcal,
        method=result.method,
        components=components,
        verdict=verdict,
    )


@_protein_router.post(
    "/stability/mutations",
    response_model=MutationScanResponse,
    summary="Scan mutations for stability",
)
async def stability_mutations(input_data: MutationScanInput) -> MutationScanResponse:
    """
    Scan single-point mutations for stability effects.

    For each position (or specified positions), all possible amino acid
    substitutions are evaluated and ranked by predicted ddG.

    Method 'empirical' uses FoldX BuildModel; 'statistical' uses
    sequence-based energy functions.

    Rate limiting: This is compute-intensive; full scans on large proteins
    may take several minutes.
    """
    try:
        from .foldx import scan_mutations
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="FoldX mutation scanning module not available.",
        )

    # Validate positions against protein length
    if input_data.positions is not None:
        for p in input_data.positions:
            if p > len(input_data.protein):
                raise HTTPException(
                    status_code=400,
                    detail=f"Position {p} exceeds protein length ({len(input_data.protein)}).",
                )

    try:
        result = scan_mutations(
            protein=input_data.protein,
            positions=input_data.positions,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Mutation scanning failed")
        raise HTTPException(
            status_code=500,
            detail=f"Mutation scanning failed: {e}",
        )

    # result is a list of MutationResult dataclasses
    mutations = [
        {
            "position": m.position,
            "original": m.original,
            "mutant": m.mutant,
            "score": m.score,
            "description": m.description,
            "details": m.details,
        }
        for m in result
    ]
    stabilizing_count = sum(
        1 for m in result if m.details and m.details.get("stabilizing")
    )
    destabilizing_count = sum(
        1 for m in result if m.details and m.details.get("destabilizing")
    )

    return MutationScanResponse(
        mutations=mutations,
        stabilizing_count=stabilizing_count,
        destabilizing_count=destabilizing_count,
    )


# ─── Solubility Endpoints ───────────────────────────────────────────

@_protein_router.post(
    "/solubility/analyze",
    response_model=SolubilityResponse,
    summary="Analyze protein solubility",
)
async def solubility_analyze(input_data: SolubilityInput) -> SolubilityResponse:
    """
    Analyze protein solubility using CamSol algorithm.

    Computes intrinsic solubility from sequence features and, if a PDB
    structure is provided, an overall solubility score incorporating
    structural corrections. Identifies aggregation-prone regions.

    Returns 400 for invalid protein or PDB input.
    """
    try:
        from .camsol import compute_solubility
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="CamSol solubility module not available.",
        )

    try:
        result = compute_solubility(
            protein=input_data.protein,
            pdb_string=input_data.pdb_string,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Solubility analysis failed")
        raise HTTPException(
            status_code=500,
            detail=f"Solubility analysis failed: {e}",
        )

    # result is a SolubilityResult dataclass, not a dict
    return SolubilityResponse(
        intrinsic_score=result.intrinsic_score,
        overall_score=result.overall_score,
        solubility_class=result.solubility_class,
        aggregation_prone_regions=result.aggregation_prone_regions,
    )


@_protein_router.post(
    "/solubility/mutations",
    response_model=MutationScanResponse,
    summary="Find solubility-improving mutations",
)
async def solubility_mutations(input_data: SolubilityInput) -> MutationScanResponse:
    """
    Find mutations that improve protein solubility.

    Scans positions in aggregation-prone regions and proposes
    amino acid substitutions that increase solubility score while
    preserving structural stability.

    Returns a list of proposed mutations with predicted solubility impact.
    """
    try:
        from .camsol import find_solubility_mutations
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="CamSol solubility module not available.",
        )

    try:
        result = find_solubility_mutations(
            protein=input_data.protein,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Solubility mutation scan failed")
        raise HTTPException(
            status_code=500,
            detail=f"Solubility mutation scan failed: {e}",
        )

    # result is a list of MutationResult dataclasses
    mutations = [
        {
            "position": m.position,
            "original": m.original,
            "mutant": m.mutant,
            "score": m.score,
            "description": m.description,
            "details": m.details,
        }
        for m in result
    ]
    stabilizing = sum(1 for m in result if m.score > 0)
    destabilizing = sum(1 for m in result if m.score <= 0)

    return MutationScanResponse(
        mutations=mutations,
        stabilizing_count=stabilizing,
        destabilizing_count=destabilizing,
    )


# ─── Immunogenicity Endpoints ───────────────────────────────────────

@_protein_router.post(
    "/immunogenicity/analyze",
    response_model=ImmunogenicityResponse,
    summary="Analyze protein immunogenicity",
)
async def immunogenicity_analyze(input_data: ImmunogenicityInput) -> ImmunogenicityResponse:
    """
    Analyze protein immunogenicity by predicting T-cell and B-cell epitopes.

    Uses MHC binding prediction for T-cell epitope identification and
    surface accessibility / physicochemical properties for B-cell epitope
    prediction. Returns an overall immunogenicity score and candidate
    positions for deimmunization.

    If specific MHC alleles are not provided, common alleles for the
    specified organism are used.
    """
    try:
        from .immunogenicity import compute_immunogenicity
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="Immunogenicity module not available.",
        )

    try:
        result = compute_immunogenicity(
            protein=input_data.protein,
            mhc_alleles=input_data.mhc_alleles,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Immunogenicity analysis failed")
        raise HTTPException(
            status_code=500,
            detail=f"Immunogenicity analysis failed: {e}",
        )

    # result is an ImmunogenicityResult dataclass
    deimmunization_candidates = [
        {
            "position": m.position,
            "original": m.original,
            "mutant": m.mutant,
            "score": m.score,
            "description": m.description,
        }
        for m in result.deimmunization_candidates
    ]

    return ImmunogenicityResponse(
        overall_score=result.overall_score,
        immunogenicity_class=result.immunogenicity_class,
        t_cell_epitopes=result.t_cell_epitopes,
        b_cell_epitopes=result.b_cell_epitopes,
        deimmunization_candidates=deimmunization_candidates,
    )


@_protein_router.post(
    "/immunogenicity/deimmunize",
    response_model=DeimmunizeResponse,
    summary="Deimmunize a protein",
)
async def immunogenicity_deimmunize(input_data: DeimmunizeInput) -> DeimmunizeResponse:
    """
    Deimmunize a protein by introducing conservative amino acid
    substitutions that reduce immunogenicity while preserving
    protein function.

    Constraints:
    - Maximum number of mutations controlled by max_mutations
    - Only substitutions with BLOSUM62 score >= blosum62_min are considered
    - Process continues until target_score is reached or max_mutations exceeded

    Returns the optimized protein, applied mutations, and before/after scores.
    """
    try:
        from .deimmunization import deimmunize
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="Deimmunization module not available.",
        )

    try:
        result = deimmunize(
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
        raise HTTPException(
            status_code=500,
            detail=f"Deimmunization failed: {e}",
        )

    # result is a DeimmunizationResult dataclass
    return DeimmunizeResponse(
        optimized_protein=result.optimized_protein,
        mutations_applied=len(result.mutations_applied),
        original_score=result.original_immunogenicity,
        optimized_score=result.optimized_immunogenicity,
        success=result.success,
    )


# ─── Full Assessment Endpoint ────────────────────────────────────────

@_protein_router.post(
    "/assessment/full",
    response_model=FullAssessmentResponse,
    summary="Full protein assessment",
)
async def assessment_full(input_data: FullAssessmentInput) -> FullAssessmentResponse:
    """
    Run a full protein assessment combining all analyses.

    Optionally includes structure prediction, stability analysis,
    solubility analysis, and immunogenicity assessment. Also runs
    predicates on the optimized sequence for completeness.

    This is the most comprehensive endpoint — suitable for final
    protein design validation before synthesis.

    Compute time: May take 2-5 minutes depending on protein length
    and which analyses are enabled.

    Rate limiting: Consumes multiple rate-limit units (one per analysis).
    """
    structure_quality_result: dict | None = None
    stability_result: dict | None = None
    solubility_result: dict | None = None
    immunogenicity_result: dict | None = None
    predicate_results: list[dict] = []
    recommendations: list[str] = []

    # ── Structure Prediction & Quality ─────────────────────────────
    if input_data.run_structure:
        pdb_str = input_data.pdb_string
        if pdb_str is None:
            # Predict structure first
            try:
                from .esmfold import predict_structure, is_esmfold_available

                if is_esmfold_available():
                    pred = predict_structure(
                        protein=input_data.protein,
                    )
                    pdb_str = pred.pdb_string
                    mean_plddt = float(getattr(pred, "mean_plddt", 0.0))

                    if mean_plddt < 50:
                        recommendations.append(
                            f"Predicted structure has low confidence "
                            f"(mean pLDDT={mean_plddt:.1f}). Consider experimental "
                            f"structure determination."
                        )
                else:
                    recommendations.append(
                        "ESMFold is not available. Structure prediction skipped. "
                        "Provide a PDB string or install biocompiler[esmfold]."
                    )
            except ImportError:
                recommendations.append(
                    "ESMFold module not available. Structure prediction skipped."
                )
            except Exception as e:
                logger.warning("Structure prediction failed: %s", e)
                recommendations.append(
                    f"Structure prediction failed: {e}"
                )

        # Assess structure quality if we have a PDB string
        if pdb_str is not None:
            try:
                from .structure.quality import compute_structure_quality

                quality = compute_structure_quality(pdb_string=pdb_str)
                structure_quality_result = quality

                if quality["verdict"] == "FAIL":
                    recommendations.append(
                        "Structure quality is poor. Results from stability and "
                        "solubility analyses may be unreliable."
                    )
                elif quality["verdict"] == "WARN":
                    recommendations.append(
                        "Structure quality is acceptable but not ideal. "
                        "Consider verifying with experimental data."
                    )
            except ImportError:
                recommendations.append(
                    "Structure quality module not available. Quality assessment skipped."
                )
            except Exception as e:
                logger.warning("Structure quality assessment failed: %s", e)
                recommendations.append(
                    f"Structure quality assessment failed: {e}"
                )

    # ── Stability Analysis ─────────────────────────────────────────
    if input_data.run_stability:
        try:
            from .foldx import empirical_stability

            stab = empirical_stability(
                protein=input_data.protein,
            )
            # stab is a FoldXResult dataclass — build dict for response
            stability_result = {
                "stability_kcal": stab.stability_kcal,
                "method": stab.method,
                "success": stab.success,
            }
            # Derive verdict from stability_kcal
            if stab.stability_kcal < -5.0:
                stability_result["verdict"] = "STABLE"
            elif stab.stability_kcal < 0.0:
                stability_result["verdict"] = "MARGINAL"
            else:
                stability_result["verdict"] = "UNSTABLE"

            if stability_result["verdict"] == "UNSTABLE":
                recommendations.append(
                    f"Protein is predicted to be unstable "
                    f"(ΔG={stab.stability_kcal:.1f} kcal/mol). "
                    f"Consider stability-enhancing mutations."
                )
            elif stability_result["verdict"] == "MARGINAL":
                recommendations.append(
                    f"Protein stability is marginal "
                    f"(ΔG={stab.stability_kcal:.1f} kcal/mol). "
                    f"Monitor during expression."
                )
        except ImportError:
            recommendations.append(
                "FoldX stability module not available. Stability analysis skipped."
            )
        except Exception as e:
            logger.warning("Stability analysis failed: %s", e)
            recommendations.append(f"Stability analysis failed: {e}")

    # ── Solubility Analysis ────────────────────────────────────────
    if input_data.run_solubility:
        try:
            from .camsol import compute_solubility

            sol = compute_solubility(
                protein=input_data.protein,
                pdb_string=input_data.pdb_string,
            )
            # sol is a SolubilityResult dataclass — build dict for response
            solubility_result = {
                "intrinsic_score": sol.intrinsic_score,
                "overall_score": sol.overall_score,
                "solubility_class": sol.solubility_class,
                "aggregation_prone_regions": sol.aggregation_prone_regions,
            }

            if sol.solubility_class in ("low", "very_low", "marginally_soluble", "insoluble"):
                recommendations.append(
                    f"Protein solubility is {sol.solubility_class} "
                    f"(score={sol.overall_score:.2f}). "
                    f"Consider solubility-enhancing mutations or fusion tags."
                )
            if sol.aggregation_prone_regions:
                recommendations.append(
                    f"Found {len(sol.aggregation_prone_regions)} "
                    f"aggregation-prone region(s). Review for potential redesign."
                )
        except ImportError:
            recommendations.append(
                "CamSol solubility module not available. Solubility analysis skipped."
            )
        except Exception as e:
            logger.warning("Solubility analysis failed: %s", e)
            recommendations.append(f"Solubility analysis failed: {e}")

    # ── Immunogenicity Analysis ────────────────────────────────────
    if input_data.run_immunogenicity:
        try:
            from .immunogenicity import compute_immunogenicity

            imm = compute_immunogenicity(
                protein=input_data.protein,
                mhc_alleles=None,
            )
            # imm is an ImmunogenicityResult dataclass — build dict for response
            immunogenicity_result = {
                "overall_score": imm.overall_score,
                "immunogenicity_class": imm.immunogenicity_class,
                "t_cell_epitopes": imm.t_cell_epitopes,
                "b_cell_epitopes": imm.b_cell_epitopes,
                "deimmunization_candidates": [
                    {
                        "position": m.position,
                        "original": m.original,
                        "mutant": m.mutant,
                        "score": m.score,
                    }
                    for m in imm.deimmunization_candidates
                ],
            }

            if imm.immunogenicity_class in ("high", "very_high"):
                recommendations.append(
                    f"Protein immunogenicity is {imm.immunogenicity_class} "
                    f"(score={imm.overall_score:.2f}). "
                    f"Consider deimmunization for therapeutic applications."
                )
            if imm.deimmunization_candidates:
                recommendations.append(
                    f"Found {len(imm.deimmunization_candidates)} "
                    f"deimmunization candidate position(s)."
                )
        except ImportError:
            recommendations.append(
                "Immunogenicity module not available. Immunogenicity analysis skipped."
            )
        except Exception as e:
            logger.warning("Immunogenicity analysis failed: %s", e)
            recommendations.append(f"Immunogenicity analysis failed: {e}")

    # ── Predicate Results ──────────────────────────────────────────
    try:
        from .structure.report import assess_protein

        report = assess_protein(
            protein=input_data.protein,
            organism=input_data.organism,
            pdb_string=input_data.pdb_string,
        )
        predicate_results = report.get("predicate_results", [])
    except ImportError:
        # Fallback: try type system directly
        try:
            from .translation import translate, compute_cai
            from .optimization import optimize_sequence

            opt_result = optimize_sequence(
                target_protein=input_data.protein,
                organism=input_data.organism,
            )
            predicate_results = [
                {"predicate": p, "verdict": "PASS"}
                for p in opt_result.satisfied_predicates
            ] + [
                {"predicate": p, "verdict": "FAIL"}
                for p in opt_result.failed_predicates
            ]
        except Exception as exc:
            logger.warning("Predicate fallback evaluation failed: %s", exc)
            predicate_results = []
    except Exception as e:
        logger.warning("Predicate evaluation failed: %s", e)
        predicate_results = []

    # ── Overall Verdict ────────────────────────────────────────────
    verdicts: list[str] = []

    if structure_quality_result is not None:
        verdicts.append(structure_quality_result.get("verdict", "WARN"))
    if stability_result is not None:
        verdicts.append(stability_result.get("verdict", "WARN"))
    if solubility_result is not None:
        sol_class = solubility_result.get("solubility_class", "medium")
        if sol_class in ("high",):
            verdicts.append("PASS")
        elif sol_class == "medium":
            verdicts.append("WARN")
        else:
            verdicts.append("FAIL")
    if immunogenicity_result is not None:
        imm_class = immunogenicity_result.get("immunogenicity_class", "moderate")
        if imm_class in ("low",):
            verdicts.append("PASS")
        elif imm_class == "moderate":
            verdicts.append("WARN")
        else:
            verdicts.append("FAIL")

    # Combine verdicts: any FAIL -> FAIL, any WARN -> WARN, else PASS
    if any(v == "FAIL" for v in verdicts):
        overall_verdict = "FAIL"
    elif any(v == "WARN" for v in verdicts):
        overall_verdict = "WARN"
    elif verdicts:
        overall_verdict = "PASS"
    else:
        overall_verdict = "WARN"

    if not recommendations:
        recommendations.append("No issues detected. Protein design looks good.")

    return FullAssessmentResponse(
        structure_quality=structure_quality_result,
        stability=stability_result,
        solubility=solubility_result,
        immunogenicity=immunogenicity_result,
        predicate_results=predicate_results,
        overall_verdict=overall_verdict,
        recommendations=recommendations,
    )


# ─── Protein Analysis Batch Endpoints ────────────────────────────────


class BatchStabilityItem(StabilityInput):
    """Single item in a batch stability analysis request."""
    pass


class BatchStabilityInput(BaseModel):
    """Input for batch stability analysis endpoint."""
    proteins: list[BatchStabilityItem] = Field(
        ..., description=f"List of proteins for stability analysis (max {BATCH_PROTEIN_MAX})"
    )


class BatchStabilityResponse(BaseModel):
    """Response for batch stability analysis endpoint."""
    results: list[dict] = Field(..., description="Per-item stability analysis results")
    summary: dict = Field(..., description="Aggregate summary")


class BatchSolubilityItem(SolubilityInput):
    """Single item in a batch solubility analysis request."""
    pass


class BatchSolubilityInput(BaseModel):
    """Input for batch solubility analysis endpoint."""
    proteins: list[BatchSolubilityItem] = Field(
        ..., description=f"List of proteins for solubility analysis (max {BATCH_PROTEIN_MAX})"
    )


class BatchSolubilityResponse(BaseModel):
    """Response for batch solubility analysis endpoint."""
    results: list[dict] = Field(..., description="Per-item solubility analysis results")
    summary: dict = Field(..., description="Aggregate summary")


class BatchImmunogenicityItem(ImmunogenicityInput):
    """Single item in a batch immunogenicity analysis request."""
    pass


class BatchImmunogenicityInput(BaseModel):
    """Input for batch immunogenicity analysis endpoint."""
    proteins: list[BatchImmunogenicityItem] = Field(
        ..., description=f"List of proteins for immunogenicity analysis (max {BATCH_PROTEIN_MAX})"
    )


class BatchImmunogenicityResponse(BaseModel):
    """Response for batch immunogenicity analysis endpoint."""
    results: list[dict] = Field(..., description="Per-item immunogenicity analysis results")
    summary: dict = Field(..., description="Aggregate summary")


@_protein_router.post(
    "/stability/batch",
    response_model=BatchStabilityResponse,
    summary="Batch protein stability analysis",
)
async def stability_batch(input_data: BatchStabilityInput) -> BatchStabilityResponse:
    """
    Analyze protein stability for multiple proteins in one request.

    Each protein is analyzed independently using FoldX (empirical) and/or
    statistical potentials. One failure does not affect other items.

    Maximum batch size: 20 proteins.
    Each item consumes one rate-limit unit.
    """
    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX} proteins",
        )

    client_id = "batch_stability"  # Rate limit key for batch operations
    _check_batch_rate_limit(client_id, n)

    results: list[dict] = []
    total = n
    stable_count = 0
    marginal_count = 0
    unstable_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(_stability_batch_single, item),
                timeout=BATCH_ITEM_TIMEOUT_S,
            )
            results.append(result)
            verdict = result.get("verdict", "")
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


def _stability_batch_single(item: StabilityInput) -> dict[str, Any]:
    """Process a single stability analysis item for batch processing."""
    try:
        from .foldx import empirical_stability
        result = empirical_stability(
            protein=item.protein,
        )
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="FoldX stability module not available.",
        )
    # result is a FoldXResult dataclass — build dict for API response
    stability_kcal = result.stability_kcal
    if stability_kcal < -5.0:
        verdict = "STABLE"
    elif stability_kcal < 0.0:
        verdict = "MARGINAL"
    else:
        verdict = "UNSTABLE"
    components = {
        key: getattr(result, key)
        for key in (
            "backbone_hbond", "sidechain_hbond", "van_der_waals",
            "electrostatics", "solvation", "van_der_waals_clashes",
            "entropy_sidechain", "entropy_mainchain", "torsional_clash",
            "backbone_clash", "helix_dipole", "disulfide",
            "electrostatic_kon", "partial_covalent", "energy_ionisation",
        )
        if getattr(result, key, None) is not None
    }
    return {
        "stability_kcal": stability_kcal,
        "method": result.method,
        "components": components,
        "verdict": verdict,
    }


@_protein_router.post(
    "/solubility/batch",
    response_model=BatchSolubilityResponse,
    summary="Batch protein solubility analysis",
)
async def solubility_batch(input_data: BatchSolubilityInput) -> BatchSolubilityResponse:
    """
    Analyze protein solubility for multiple proteins in one request.

    Each protein is analyzed independently using the CamSol algorithm.
    One failure does not affect other items.

    Maximum batch size: 20 proteins.
    Each item consumes one rate-limit unit.
    """
    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX} proteins",
        )

    client_id = "batch_solubility"
    _check_batch_rate_limit(client_id, n)

    results: list[dict] = []
    total = n
    high_count = 0
    medium_count = 0
    low_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(_solubility_batch_single, item),
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


def _solubility_batch_single(item: SolubilityInput) -> dict[str, Any]:
    """Process a single solubility analysis item for batch processing."""
    try:
        from .camsol import compute_solubility
        result = compute_solubility(
            protein=item.protein,
            pdb_string=item.pdb_string,
        )
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="CamSol solubility module not available.",
        )
    # result is a SolubilityResult dataclass — build dict for API response
    return {
        "intrinsic_score": result.intrinsic_score,
        "overall_score": result.overall_score,
        "solubility_class": result.solubility_class,
        "aggregation_prone_regions": result.aggregation_prone_regions,
    }


@_protein_router.post(
    "/immunogenicity/batch",
    response_model=BatchImmunogenicityResponse,
    summary="Batch protein immunogenicity analysis",
)
async def immunogenicity_batch(input_data: BatchImmunogenicityInput) -> BatchImmunogenicityResponse:
    """
    Analyze protein immunogenicity for multiple proteins in one request.

    Each protein is analyzed independently for T-cell and B-cell epitopes.
    One failure does not affect other items.

    Maximum batch size: 20 proteins.
    Each item consumes one rate-limit unit.
    """
    n = len(input_data.proteins)
    if n > BATCH_PROTEIN_MAX:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size {n} exceeds maximum of {BATCH_PROTEIN_MAX} proteins",
        )

    client_id = "batch_immunogenicity"
    _check_batch_rate_limit(client_id, n)

    results: list[dict] = []
    total = n
    low_count = 0
    moderate_count = 0
    high_count = 0
    error_count = 0

    for item in input_data.proteins:
        try:
            result = await asyncio.wait_for(
                asyncio.to_thread(_immunogenicity_batch_single, item),
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


def _immunogenicity_batch_single(item: ImmunogenicityInput) -> dict[str, Any]:
    """Process a single immunogenicity analysis item for batch processing."""
    try:
        from .immunogenicity import compute_immunogenicity
        result = compute_immunogenicity(
            protein=item.protein,
            mhc_alleles=item.mhc_alleles,
        )
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="Immunogenicity module not available.",
        )
    # result is an ImmunogenicityResult dataclass — build dict for API response
    return {
        "overall_score": result.overall_score,
        "immunogenicity_class": result.immunogenicity_class,
        "t_cell_epitopes": result.t_cell_epitopes,
        "b_cell_epitopes": result.b_cell_epitopes,
    }


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

    Protein Analysis endpoints:
    - POST /protein/structure/predict       — Predict protein structure
    - POST /protein/structure/batch         — Batch structure prediction
    - POST /protein/structure/quality       — Assess structure quality
    - POST /protein/stability/analyze       — Analyze stability
    - POST /protein/stability/mutations     — Scan stability mutations
    - POST /protein/stability/batch         — Batch stability analysis
    - POST /protein/solubility/analyze      — Analyze solubility
    - POST /protein/solubility/mutations    — Find solubility mutations
    - POST /protein/solubility/batch        — Batch solubility analysis
    - POST /protein/immunogenicity/analyze  — Analyze immunogenicity
    - POST /protein/immunogenicity/deimmunize — Deimmunize a protein
    - POST /protein/immunogenicity/batch    — Batch immunogenicity analysis
    - POST /protein/assessment/full         — Full protein assessment
    """
    from . import __version__

    # CORS origins from environment (default: no CORS — restrictive)
    # SECURITY: Never default to "*" with credentials. If no origins
    # are explicitly configured, CORS is disabled entirely.
    _raw_cors = os.environ.get("BIOCOMPILER_CORS_ORIGINS", "").strip()
    cors_origins = [o.strip() for o in _raw_cors.split(",") if o.strip()] if _raw_cors else []

    # Credentials only when explicitly enabled AND origins are not "*"
    # (CORS spec forbids allow_credentials=True with allow_origins=["*"])
    _raw_creds = os.environ.get("BIOCOMPILER_CORS_ALLOW_CREDENTIALS", "").strip().lower()
    _allow_creds = _raw_creds in ("true", "1", "yes") and cors_origins != ["*"]

    app = FastAPI(
        title="BioCompiler API",
        description="Machine-verified gene design REST API with authentication and rate limiting",
        version=__version__,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # CORS middleware: only added when explicit origins are configured.
    # SECURITY: Never allow_credentials=True with allow_origins=["*"].
    if cors_origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=cors_origins,
            allow_credentials=_allow_creds,
            allow_methods=["*"],
            allow_headers=["*"],
        )
    # else: no CORS middleware added (most restrictive)

    # Request body size middleware
    @app.middleware("http")
    async def request_size_limit_middleware(request: Request, call_next) -> Any:
        """Reject requests whose body exceeds MAX_REQUEST_SIZE bytes."""
        if request.method in ("POST", "PUT", "PATCH"):
            content_length = request.headers.get("content-length")
            if content_length is not None:
                try:
                    if int(content_length) > MAX_REQUEST_SIZE:
                        from fastapi.responses import JSONResponse
                        return JSONResponse(
                            status_code=413,
                            content={
                                "detail": (
                                    f"Request body too large ({int(content_length)} bytes). "
                                    f"Maximum: {MAX_REQUEST_SIZE} bytes."
                                ),
                            },
                        )
                except (ValueError, TypeError):
                    pass
        response = await call_next(request)
        return response

    # Auth-mode middleware: add warning header when auth is optional
    # and the request is unauthenticated
    @app.middleware("http")
    async def auth_mode_middleware(request: Request, call_next) -> Any:
        response = await call_next(request)
        # In optional mode, add a warning header if no API key was provided
        if _AUTH_MODE == "optional":
            api_key = request.headers.get(API_KEY_NAME)
            if not api_key or api_key not in _CONFIGURED_API_KEYS:
                response.headers["X-Auth-Warning"] = (
                    "Authentication is optional but recommended. "
                    "Unauthenticated access may be restricted in future versions."
                )
        return response

    # Rate limiting middleware
    @app.middleware("http")
    async def rate_limit_middleware(request: Request, call_next) -> Any:
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
    async def health_check() -> HealthResponse:
        """Health check endpoint (no auth required)."""
        # Compute CORS configuration for health response
        _cors_origins_list = cors_origins  # Report actual configured origins
        _cors_allow_creds = _allow_creds   # Report actual credential setting

        return HealthResponse(
            status="healthy",
            version=__version__,
            timestamp=datetime.now(timezone.utc).isoformat(),
            auth_enabled=_is_auth_enabled(),
            rate_limit_rpm=RATE_LIMIT_RPM,
            cors_origins=_cors_origins_list,
            cors_allow_credentials=_cors_allow_creds,
        )

    @app.get("/info", response_model=InfoResponse)
    async def api_info() -> InfoResponse:
        """
        Return API configuration and limits.

        Provides runtime limits, supported organisms, and version info
        so clients can validate inputs before making requests.
        """
        # Try to get safety version from biosecurity module
        try:
            from .biosecurity import __version__ as safety_ver
        except (ImportError, AttributeError):
            safety_ver = __version__

        # Get the unique canonical organism names (filter out aliases)
        canonical_organisms = sorted(set(
            name for name in SUPPORTED_ORGANISMS
            if not name.islower()  # aliases are lowercase shortcuts
        ))

        return InfoResponse(
            max_protein_length=MAX_PROTEIN_LENGTH,
            max_dna_length=MAX_DNA_LENGTH,
            max_batch_size=MAX_BATCH_SIZE,
            max_request_size=MAX_REQUEST_SIZE,
            optimize_timeout_s=OPTIMIZE_TIMEOUT_S,
            supported_organisms=canonical_organisms,
            api_version=__version__,
            safety_version=safety_ver,
        )

    @app.get("/organisms", response_model=OrganismResponse)
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

    @app.get("/predicates", response_model=PredicateResponse)
    async def list_predicates() -> PredicateResponse:
        """List all registered type predicates."""
        from .type_system import registry
        return PredicateResponse(predicates=registry.names())

    @app.get("/enzymes")
    async def list_enzymes() -> dict[str, dict[str, str]]:
        """List all known restriction enzymes."""
        return {
            "enzymes": {
                name: site for name, site in RESTRICTION_ENZYMES.items()
            }
        }

    # ─── Provenance Endpoints ────────────────────────────────────

    @app.get("/provenance/{record_id}")
    async def get_provenance(record_id: str, client_id: str = Depends(verify_api_key)) -> dict:
        """
        Retrieve a stored provenance record by ID.

        Returns the full decision-level provenance trail for a previously
        run optimization.  Records are persisted on disk and survive
        process restarts.

        The ``record_id`` is the UUID returned in the ``provenance_id``
        field of an ``/optimize`` response.
        """
        try:
            trail = _provenance_store.load(record_id)
        except (FileNotFoundError, ValueError):
            raise HTTPException(
                status_code=404,
                detail=f"Provenance record not found: {record_id}",
            )
        return {
            "id": record_id,
            "trail": trail.to_dict(),
        }

    @app.get("/provenance")
    async def query_provenance(
        protein_name: Optional[str] = None,
        organism: Optional[str] = None,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None,
        client_id: str = Depends(verify_api_key),
    ) -> dict:
        """
        Query/list stored provenance records.

        Supports filtering by gene/protein name, organism, and date range.
        All filters are optional and combined with AND logic.

        Query parameters:
        - ``protein_name`` — match records whose gene_name equals this value
        - ``organism`` — match records whose organism equals this value
        - ``date_from`` — ISO 8601 start timestamp (inclusive)
        - ``date_to`` — ISO 8601 end timestamp (inclusive)
        """
        date_range = None
        if date_from and date_to:
            date_range = (date_from, date_to)

        trails = _provenance_store.query(
            protein_name=protein_name,
            organism=organism,
            date_range=date_range,
        )
        return {
            "count": len(trails),
            "records": [
                {
                    "gene_name": t.gene_name,
                    "organism": t.organism,
                    "timestamp": t.timestamp,
                    "total_cai": t.total_cai,
                    "total_gc": t.total_gc,
                    "codon_decision_count": len(t.codon_decisions),
                    "constraint_decision_count": len(t.constraint_decisions),
                }
                for t in trails
            ],
        }

    @app.post("/check", response_model=TypeCheckResponse)
    async def check_sequence(input_data: SequenceInput, client_id: str = Depends(verify_api_key)) -> TypeCheckResponse:
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
                boundaries=exon_boundaries,
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
                except CertificateGenerationError as exc:
                    logger.warning("Certificate generation failed during type-check: %s", exc)

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
        except Exception as e:
            logger.exception("Unexpected error during type-check")
            raise HTTPException(status_code=500, detail=f"Type-check failed: {e}")

    @app.post("/optimize", response_model=OptimizeResponse)
    async def optimize(input_data: ProteinInput, client_id: str = Depends(verify_api_key)) -> OptimizeResponse:
        """
        Optimize a DNA sequence for a target protein.

        Uses z3 constraint solver with greedy fallback to find
        a sequence that satisfies all type predicates.

        The ``organism_domain`` parameter controls whether eukaryotic
        constraints (cryptic splice sites, CpG islands) are applied.
        When set to ``"auto"`` (default), the domain is detected from
        the organism name.  Set to ``"eukaryote"`` or ``"prokaryote"``
        to override auto-detection.

        Immunogenicity predicates (LowImmunogenicity, NoStrongTCellEpitope,
        NoDominantBCellEpitope, PopulationCoverageSafe) are context-aware:
        they consider ``source_organism``, ``therapeutic``, and
        ``self_protein`` to determine self-protein status and apply
        appropriate thresholds.
        """
        try:
            resolved_domain = resolve_organism_domain(
                input_data.organism, input_data.organism_domain
            )

            # Run optimization with timeout
            try:
                result = await asyncio.wait_for(
                    asyncio.to_thread(
                        optimize_sequence,
                        target_protein=input_data.protein,
                        organism=input_data.organism,
                        gc_lo=input_data.gc_lo,
                        gc_hi=input_data.gc_hi,
                        cai_threshold=input_data.cai_threshold,
                        restriction_sites=input_data.enzymes,
                        cryptic_splice_threshold=input_data.cryptic_splice_threshold,
                        organism_domain=resolved_domain,
                        source_organism=input_data.source_organism,
                        therapeutic=input_data.therapeutic,
                        self_protein=input_data.self_protein,
                        track_provenance=input_data.track_provenance,
                        strict_mode=input_data.strict_mode,
                    ),
                    timeout=OPTIMIZE_TIMEOUT_S,
                )
            except asyncio.TimeoutError:
                raise HTTPException(
                    status_code=504,
                    detail=(
                        f"Optimization timed out after {OPTIMIZE_TIMEOUT_S} seconds. "
                        f"Try a shorter protein or relax constraints."
                    ),
                )

            # Persist decision-level provenance trail if tracking is enabled
            provenance_id = None
            if input_data.track_provenance and result.decision_trail is not None:
                try:
                    provenance_id = _store_provenance(result.decision_trail)
                except Exception as exc:
                    logger.warning(
                        "Failed to persist provenance trail: %s", exc
                    )

            return OptimizeResponse(
                sequence=result.sequence,
                protein=result.protein,
                cai=result.cai,
                gc_content=result.gc_content,
                satisfied_predicates=result.satisfied_predicates,
                failed_predicates=result.failed_predicates,
                fallback_used=result.fallback_used,
                provenance_id=provenance_id,
                organism_domain=resolved_domain,
                source_organism=input_data.source_organism,
                therapeutic=input_data.therapeutic,
                self_protein=input_data.self_protein,
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
        except HTTPException:
            raise  # Re-raise HTTPExceptions (including our 504 timeout)
        except Exception as e:
            logger.exception("Unexpected error during optimization")
            raise HTTPException(status_code=500, detail=f"Optimization failed: {e}")

    @app.post("/verify", response_model=VerifyResponse)
    async def verify(input_data: CertificateInput, client_id: str = Depends(verify_api_key)) -> VerifyResponse:
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
    async def scan(input_data: ScanInput, client_id: str = Depends(verify_api_key)) -> ScanResponse:
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
        except Exception as e:
            logger.exception("Unexpected error during sequence scan")
            raise HTTPException(status_code=500, detail=f"Scan failed: {e}")

    @app.post("/export/fasta")
    async def export_fasta_endpoint(input_data: ExportFastaInput, client_id: str = Depends(verify_api_key)) -> dict[str, str]:
        """Export a sequence in FASTA format."""
        try:
            fasta = export_fasta(
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

    @app.post("/export/genbank")
    async def export_genbank_endpoint(input_data: ExportGenbankInput, client_id: str = Depends(verify_api_key)) -> dict[str, str]:
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
        except (InvalidSequenceError, ValueError, KeyError) as e:
            raise HTTPException(status_code=422, detail=str(e))
        except Exception as e:
            logger.exception("GenBank export failed unexpectedly")
            raise HTTPException(status_code=500, detail=f"GenBank export failed: {e}")

    # ─── SBOL3 Export Endpoint ────────────────────────────────────

    @app.post("/export/sbol3")
    async def export_sbol3_endpoint(input_data: ExportSbol3Input, client_id: str = Depends(verify_api_key)) -> dict:
        """Export a sequence in SBOL3 format."""
        try:
            # Generate SBOL3 XML content
            seq = input_data.sequence.upper()
            organism = input_data.organism
            gene_name = input_data.gene_name
            fmt = input_data.format

            if fmt == "sbol3json":
                import json as _json
                sbol3_content = _json.dumps({
                    "SBOL3": True,
                    "type": "http://sbols.org/v3#Component",
                    "identity": f"https://biocompiler.org/{gene_name}",
                    "displayId": gene_name,
                    "sequence": {
                        "type": "http://sbols.org/v3#Sequence",
                        "identity": f"https://biocompiler.org/{gene_name}_seq",
                        "displayId": f"{gene_name}_seq",
                        "elements": seq,
                        "encoding": "http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html",
                    },
                    "organism": organism,
                }, indent=2)
            else:
                sbol3_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:sbol="http://sbols.org/v3#"
         xmlns:dcterms="http://purl.org/dc/terms/">
  <sbol:Component rdf:about="https://biocompiler.org/{gene_name}">
    <dcterms:title>{gene_name}</dcterms:title>
    <sbol:displayId>{gene_name}</sbol:displayId>
    <sbol:hasSequence rdf:resource="https://biocoder.org/{gene_name}_seq"/>
    <sbol:organism>{organism}</sbol:organism>
  </sbol:Component>
  <sbol:Sequence rdf:about="https://biocompiler.org/{gene_name}_seq">
    <sbol:displayId>{gene_name}_seq</sbol:displayId>
    <sbol:elements>{seq}</sbol:elements>
    <sbol:encoding rdf:resource="http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"/>
  </sbol:Sequence>
</rdf:RDF>"""
            return {"format": fmt, "content": sbol3_content}
        except Exception as e:
            logger.exception("SBOL3 export failed unexpectedly")
            raise HTTPException(status_code=500, detail=f"SBOL3 export failed: {e}")

    # ─── Batch Endpoints ─────────────────────────────────────────

    @app.post("/batch/check", response_model=BatchCheckResponse)
    async def batch_check(
        input_data: BatchCheckInput,
        client_id: str = Depends(verify_api_key),
    ) -> BatchCheckResponse:
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
    ) -> BatchOptimizeResponse:
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

    @app.post("/batch/optimize/fast", response_model=BatchOptimizeResponse)
    async def fast_batch_optimize(
        input_data: FastBatchOptimizeInput,
        client_id: str = Depends(verify_api_key),
    ) -> BatchOptimizeResponse:
        """
        Fast batch optimization with shared organism/parameters.

        All proteins are optimized for the same organism with the same
        GC range and restriction enzyme settings.  This allows the
        backend to reuse a single HybridOptimizer instance, avoiding
        redundant precomputation of codon tables, restriction site
        lookups, and GC/GT/AG data structures for each protein.

        This endpoint is significantly faster than POST /batch/optimize
        for high-throughput workflows where all proteins target the
        same organism.

        Rate limiting: each protein consumes one rate-limit unit.
        Maximum batch size: 20 proteins.
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

        # Build response
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

    @app.post("/batch/export", response_model=BatchExportResponse)
    async def batch_export(
        input_data: BatchExportInput,
        client_id: str = Depends(verify_api_key),
    ) -> BatchExportResponse:
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

    @app.post("/validate-datasets")
    async def validate_datasets(
        client_id: str = Depends(verify_api_key),
        datasets: str | None = None,
        include_cross_organism: bool = True,
        include_optimization_improvement: bool = True,
    ) -> dict[str, Any]:
        """
        Validate the optimizer against common biological datasets.

        Tests the optimizer on well-known gene sequences from human, E. coli,
        yeast, and synthetic benchmark proteins. Validates translation fidelity,
        GC content, CAI bounds, protein length, and cross-organism consistency.

        Returns a structured report with pass/fail status for each test.
        """
        from .dataset_validation import run_dataset_validation, format_dataset_report_json

        ds_list = datasets.split(",") if datasets else None
        report = await asyncio.to_thread(
            run_dataset_validation,
            datasets=ds_list,
            include_cross_organism=include_cross_organism,
            include_optimization_improvement=include_optimization_improvement,
        )

        return json.loads(format_dataset_report_json(report))

    # ─── Protein Analysis Router ───────────────────────────────────
    # Mount at /protein/ for clean URLs (primary)
    app.include_router(_protein_router, prefix="/protein")
    return app


# Create the default app instance
app = create_app()
