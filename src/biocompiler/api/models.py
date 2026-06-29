"""
BioCompiler API — Pydantic request/response models and related constants.

This module contains ALL Pydantic models used by the REST API, along with
the size-limit constants they reference in field validators.
"""

import os
from typing import Optional

from pydantic import BaseModel, Field, field_validator

from ..organisms import SUPPORTED_ORGANISMS, CODON_USAGE_TABLES, resolve_organism, ORGANISM_ALIASES
from biocompiler.organisms.config import is_eukaryotic_organism

# ─── Input Size Limits ────────────────────────────────────────────

MAX_PROTEIN_LENGTH = int(os.environ.get("BIOCOMPILER_MAX_PROTEIN_LENGTH", "10000"))  # aa
MAX_PROTEIN_SEQUENCE_LENGTH = MAX_PROTEIN_LENGTH  # backward-compatible alias
MAX_BATCH_SIZE = int(os.environ.get("BIOCOMPILER_MAX_BATCH_SIZE", "50"))  # sequences per batch
MAX_REQUEST_SIZE = int(os.environ.get("BIOCOMPILER_MAX_REQUEST_SIZE", str(10_000_000)))  # bytes
OPTIMIZE_TIMEOUT_S = int(os.environ.get("BIOCOMPILER_OPTIMIZE_TIMEOUT", "300"))  # seconds
MAX_DNA_LENGTH = int(os.environ.get("BIOCOMPILER_MAX_DNA_LENGTH", "100000"))  # bases
MAX_DNA_SEQUENCE_LENGTH = MAX_DNA_LENGTH  # backward-compatible alias

ESMFOLD_TIMEOUT_S = 120  # Default timeout for ESMFold structure prediction

_PROTEIN_VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# ─── Batch Limits ─────────────────────────────────────────────────

BATCH_CHECK_MAX = 50
BATCH_OPTIMIZE_MAX = 20
BATCH_EXPORT_MAX = 50
BATCH_ITEM_TIMEOUT_S = int(os.environ.get("BIOCOMPILER_BATCH_ITEM_TIMEOUT", "30"))
BATCH_PROTEIN_MAX = 20


# ─── Validation Helpers ───────────────────────────────────────────

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
    """
    if not organism:
        return "Organism must not be empty."
    try:
        resolved = resolve_organism(organism, strict=False)
        if resolved not in SUPPORTED_ORGANISMS:
            return f"Unsupported organism: {organism} (resolved to {resolved!r}). Supported: {SUPPORTED_ORGANISMS}"
    except ImportError:
        pass
    return None


# ─── Organism Domain Resolution ────────────────────────────────────

# Imported from application service to avoid duplication.
# The canonical definition lives in biocompiler.application.optimization_service.
from ..application.optimization_service import resolve_organism_domain  # noqa: E402


# ═══════════════════════════════════════════════════════════════════
# Pydantic Input Models
# ═══════════════════════════════════════════════════════════════════


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
                    f"Supported: {sorted(set(SUPPORTED_ORGANISMS))}"
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


# ═══════════════════════════════════════════════════════════════════
# Pydantic Response Models
# ═══════════════════════════════════════════════════════════════════


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


# ═══════════════════════════════════════════════════════════════════
# Batch Pydantic Models
# ═══════════════════════════════════════════════════════════════════


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
    """Input for fast batch optimize endpoint (shared organism/parameters)."""
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


# ═══════════════════════════════════════════════════════════════════
# Protein Analysis Input Models
# ═══════════════════════════════════════════════════════════════════


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


# ═══════════════════════════════════════════════════════════════════
# Protein Analysis Response Models
# ═══════════════════════════════════════════════════════════════════


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
    verdict_enum: str | None = Field(
        None, description="Stability verdict as Verdict enum value (PASS, UNCERTAIN, or FAIL)"
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
    overall_verdict_enum: str | None = Field(
        None, description="Overall verdict as Verdict enum value (PASS, LIKELY_PASS, UNCERTAIN, LIKELY_FAIL, FAIL)"
    )
    recommendations: list[str] = Field(
        ..., description="Actionable recommendations based on assessment"
    )


# ═══════════════════════════════════════════════════════════════════
# Protein Analysis Batch Models
# ═══════════════════════════════════════════════════════════════════


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
