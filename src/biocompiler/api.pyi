"""Type stubs for biocompiler.api — REST API endpoints, Pydantic models, and auth."""

from __future__ import annotations

from typing import Any, Optional

from fastapi import FastAPI
from fastapi.security import APIKeyHeader
from pydantic import BaseModel, Field


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

API_KEY_NAME: str
RATE_LIMIT_RPM: int
BATCH_CHECK_MAX: int
BATCH_OPTIMIZE_MAX: int
BATCH_EXPORT_MAX: int
BATCH_ITEM_TIMEOUT_S: int
ESMFOLD_TIMEOUT_S: int
MAX_PROTEIN_LENGTH: int
MAX_BATCH_SIZE: int
MAX_REQUEST_SIZE: int
OPTIMIZE_TIMEOUT_S: int


# ────────────────────────────────────────────────────────────
# Auth
# ────────────────────────────────────────────────────────────

async def verify_api_key(api_key: str = ...) -> str: ...
def set_no_auth_flag() -> None: ...
def get_auth_mode() -> str: ...
def get_configured_api_keys() -> set[str]: ...
def is_auth_enabled() -> bool: ...


# ────────────────────────────────────────────────────────────
# Pydantic Input Models
# ────────────────────────────────────────────────────────────

class SequenceInput(BaseModel):
    sequence: str
    organism: str
    species: Optional[str]
    exon_boundaries: Optional[list[tuple[int, int]]]
    gc_lo: float
    gc_hi: float
    cai_threshold: float
    enzymes: Optional[list[str]]
    cellular_context: str


class ProteinInput(BaseModel):
    protein: str
    organism: str
    species: Optional[str]
    gc_lo: float
    gc_hi: float
    cai_threshold: float
    enzymes: Optional[list[str]]
    cryptic_splice_threshold: float
    track_provenance: bool
    organism_domain: str
    source_organism: Optional[str]
    therapeutic: bool
    self_protein: Optional[bool]
    strict_mode: bool


class CertificateInput(BaseModel):
    certificate: dict


class ExportFastaInput(BaseModel):
    sequence: str
    identifier: str
    description: str
    organism: str


class ExportGenbankInput(BaseModel):
    sequence: str
    locus_name: str
    definition: str
    organism: str
    gene_name: Optional[str]
    exon_boundaries: Optional[list[tuple[int, int]]]
    certificate: Optional[dict]


class ScanInput(BaseModel):
    sequence: str
    enzymes: Optional[list[str]]
    find_orfs: bool


class BatchCheckItem(BaseModel):
    sequence: str
    organism: str
    exon_boundaries: Optional[list[tuple[int, int]]]
    gc_lo: float
    gc_hi: float
    cai_threshold: float
    enzymes: Optional[list[str]]
    cellular_context: str


class BatchCheckInput(BaseModel):
    sequences: list[BatchCheckItem]


class BatchOptimizeItem(BaseModel):
    protein: str
    organism: str
    gc_lo: float
    gc_hi: float
    cai_threshold: float
    enzymes: Optional[list[str]]
    cryptic_splice_threshold: float
    organism_domain: str
    source_organism: Optional[str]
    therapeutic: bool
    self_protein: Optional[bool]


class StructurePredictInput(BaseModel):
    protein: str


class QualityAssessInput(BaseModel):
    pdb_content: str


class StabilityInput(BaseModel):
    protein: str


class MutationScanInput(BaseModel):
    protein: str


class SolubilityInput(BaseModel):
    protein: str


class ImmunogenicityInput(BaseModel):
    protein: str
    alleles: Optional[list[str]]


class DeimmunizeInput(BaseModel):
    protein: str
    alleles: Optional[list[str]]


class FullAssessmentInput(BaseModel):
    protein: str


# ────────────────────────────────────────────────────────────
# Pydantic Response Models
# ────────────────────────────────────────────────────────────

class TypeCheckResponse(BaseModel):
    sequence_length: int
    gc_content: float
    protein: str
    results: list[dict]
    overall_verdict: str
    certificate: Optional[dict]


class OptimizeResponse(BaseModel):
    sequence: str
    protein: str
    cai: float
    gc_content: float
    satisfied_predicates: list[str]
    failed_predicates: list[str]
    fallback_used: bool
    provenance_id: Optional[str]
    organism_domain: str
    source_organism: Optional[str]
    therapeutic: bool
    self_protein: Optional[bool]


class VerifyResponse(BaseModel):
    status: str
    failure_reasons: list[str]


class ScanResponse(BaseModel):
    sequence_length: int
    tokens: list[dict]
    orfs: Optional[list[dict]]


class OrganismResponse(BaseModel):
    organisms: list[dict]


class PredicateResponse(BaseModel):
    predicates: list[str]


class HealthResponse(BaseModel):
    status: str
    version: str
    timestamp: str
    auth_enabled: bool
    rate_limit_rpm: int


class InfoResponse(BaseModel):
    max_protein_length: int
    max_batch_size: int
    max_request_size: int
    optimize_timeout_s: int
    supported_organisms: list[str]
    api_version: str
    safety_version: str


class BatchCheckSummary(BaseModel):
    total: int
    pass_: int
    fail: int
    uncertain: int
    errors: int


class BatchCheckResponse(BaseModel):
    results: list[dict]
    summary: BatchCheckSummary


class ProvenanceResponse(BaseModel):
    count: int
    records: list[dict]


class ProvenanceExplainResponse(BaseModel):
    explanation: str


class ProvenanceReportResponse(BaseModel):
    report: str


class ExportFastaResponse(BaseModel):
    format: str
    content: str


class ExportGenbankResponse(BaseModel):
    format: str
    content: str


class EnzymeListResponse(BaseModel):
    enzymes: dict[str, str]


class ProvenanceDetailResponse(BaseModel):
    id: str
    trail: dict


class ProvenanceRecordSummary(BaseModel):
    gene_name: Optional[str]
    organism: Optional[str]
    timestamp: Optional[str]
    total_cai: Optional[float]
    total_gc: Optional[float]
    codon_decision_count: int
    constraint_decision_count: int


class ProvenanceListResponse(BaseModel):
    count: int
    records: list[ProvenanceRecordSummary]


class ExportSbol3Input(BaseModel):
    sequence: str
    organism: str
    gene_name: str
    base_uri: str
    format: str


class ExportSbol3Response(BaseModel):
    format: str
    content: str


class DatasetValidationResult(BaseModel):
    dataset: str
    test_name: str
    passed: bool
    expected: Optional[str]
    actual: Optional[str]


class DatasetValidationResponse(BaseModel):
    total_tests: int
    passed: int
    failed: int
    results: list[DatasetValidationResult]


# ────────────────────────────────────────────────────────────
# Organism domain resolution
# ────────────────────────────────────────────────────────────

def resolve_organism_domain(organism: str, domain: str = ...) -> str: ...


# ────────────────────────────────────────────────────────────
# FastAPI app factory
# ────────────────────────────────────────────────────────────

def create_app() -> FastAPI: ...

app: FastAPI
