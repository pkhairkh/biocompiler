"""
BioCompiler REST API — FastAPI Web Service

Production-grade REST API for integration with other bioinformatics tools.

Endpoints:
- POST /check       — Type-check a DNA sequence
- POST /optimize    — Optimize DNA for a target protein
- POST /verify      — Verify a certificate
- POST /scan        — Scan a sequence for motifs
- POST /export/fasta   — Export sequence as FASTA
- POST /export/genbank — Export sequence as GenBank
- GET  /organisms   — List supported organisms
- GET  /predicates  — List registered predicates
- GET  /health      — Health check

All endpoints accept and return JSON. Certificate data is embedded
directly in responses for seamless pipeline integration.
"""

import logging
from datetime import datetime, timezone
from typing import Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, validator

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
    UnsupportedOrganismError, InvalidProteinError,
)

logger = logging.getLogger(__name__)

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

    @validator("sequence")
    def validate_seq(cls, v):
        v = v.upper()
        invalid = set(v) - set("ACGTN")
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v

    @validator("organism")
    def validate_organism(cls, v):
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

    @validator("protein")
    def validate_protein(cls, v):
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


# ─── FastAPI Application ──────────────────────────────────────────

def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    from . import __version__

    app = FastAPI(
        title="BioCompiler API",
        description="Machine-verified gene design REST API",
        version=__version__,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # CORS middleware for cross-origin access from web tools
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # ─── Endpoints ──────────────────────────────────────────────

    @app.get("/health", response_model=HealthResponse)
    async def health_check():
        """Health check endpoint."""
        return HealthResponse(
            status="healthy",
            version=__version__,
            timestamp=datetime.now(timezone.utc).isoformat(),
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
    async def check_sequence(input_data: SequenceInput):
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
    async def optimize(input_data: ProteinInput):
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
    async def verify(input_data: CertificateInput):
        """
        Independently verify a guarantee certificate.

        Re-evaluates every predicate from scratch using only the
        sequence and parameters embedded in the certificate.
        """
        try:
            status, failures = verify_certificate(input_data.certificate)
            return VerifyResponse(status=status, failure_reasons=failures)
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Verification error: {e}")

    @app.post("/scan", response_model=ScanResponse)
    async def scan(input_data: ScanInput):
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
    async def export_fasta_endpoint(input_data: ExportFastaInput):
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
    async def export_genbank_endpoint(input_data: ExportGenbankInput):
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

    return app


# Create the default app instance
app = create_app()
