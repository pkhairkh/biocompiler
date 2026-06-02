"""
BioCompiler REST API — Phase 2+3 Endpoints

REST API endpoints for protein structure prediction (ESMFold),
structure quality assessment, stability analysis, solubility analysis,
immunogenicity assessment, deimmunization, and full protein assessment.

Endpoints:
- POST /phase2_3/structure/predict        — Predict protein structure via ESMFold
- POST /phase2_3/structure/quality         — Assess structure quality from PDB
- POST /phase2_3/stability/analyze         — Analyze protein stability
- POST /phase2_3/stability/mutations        — Scan mutations for stability
- POST /phase2_3/solubility/analyze        — Analyze protein solubility
- POST /phase2_3/solubility/mutations       — Find solubility-improving mutations
- POST /phase2_3/immunogenicity/analyze    — Analyze immunogenicity
- POST /phase2_3/immunogenicity/deimmunize — Deimmunize a protein
- POST /phase2_3/assessment/full           — Full protein assessment (all Phase 2+3)

Rate limiting: These endpoints are compute-intensive. They share the same
rate-limit pool as Phase 1 endpoints (60 req/min by default). ESMFold
calls have a 120-second default timeout.

HTTP Status Codes:
- 200: Success
- 400: Bad input (invalid protein, organism, PDB string)
- 503: Service unavailable (ESMFold offline, FoldX not installed, etc.)
- 422: Validation error
- 500: Internal error
"""

from __future__ import annotations

import logging
import time
from typing import Any

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, field_validator

logger = logging.getLogger(__name__)

# ─── Constants ────────────────────────────────────────────────────────

ESMFOLD_TIMEOUT_S = 120  # Default timeout for ESMFold structure prediction

_VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


# ─── Helper Functions ─────────────────────────────────────────────────

def validate_protein_input(protein: str) -> str | None:
    """
    Validate a protein sequence string.

    Returns an error message string if invalid, or None if valid.
    Checks: non-empty, only standard amino acid single-letter codes.
    """
    if not protein:
        return "Protein sequence must not be empty."
    protein_upper = protein.upper()
    invalid = set(protein_upper) - _VALID_AMINO_ACIDS
    if invalid:
        return f"Invalid amino acids in protein: {invalid}. " \
               f"Allowed: {sorted(_VALID_AMINO_ACIDS)}"
    if len(protein_upper) > 5000:
        return f"Protein sequence too long ({len(protein_upper)} aa). Maximum: 5000 aa."
    return None


def validate_organism_input(organism: str) -> str | None:
    """
    Validate an organism string.

    Returns an error message string if invalid, or None if valid.
    Checks: non-empty, supported organism.
    """
    if not organism:
        return "Organism must not be empty."
    try:
        from .organisms import SUPPORTED_ORGANISMS
        if organism not in SUPPORTED_ORGANISMS:
            return f"Unsupported organism: {organism}. Supported: {SUPPORTED_ORGANISMS}"
    except ImportError:
        # Fallback: accept any non-empty string if organisms module unavailable
        pass
    return None


# ─── Pydantic Input Models ────────────────────────────────────────────

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


class FullAssessmentInput(BaseModel):
    """Input for full protein assessment (all Phase 2+3 analyses)."""
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


# ─── Pydantic Response Models ─────────────────────────────────────────

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
        ..., description="Intrinsic solubility score (0-1, higher = more soluble)"
    )
    overall_score: float = Field(
        ..., description="Overall solubility score including structural effects (0-1)"
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
        ..., description="Phase 1 predicate results for the optimized sequence"
    )
    overall_verdict: str = Field(
        ..., description="Overall assessment verdict: PASS, WARN, or FAIL"
    )
    recommendations: list[str] = Field(
        ..., description="Actionable recommendations based on assessment"
    )


# ─── Router ────────────────────────────────────────────────────────────

router = APIRouter(prefix="/phase2_3", tags=["Phase 2+3"])


# ─── Structure Endpoints ──────────────────────────────────────────────

@router.post(
    "/structure/predict",
    response_model=StructurePredictResponse,
    summary="Predict protein structure via ESMFold",
)
async def structure_predict(input_data: StructurePredictInput):
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
        from biocompiler.esmfold import predict_structure, is_esmfold_available
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
            organism=input_data.organism,
            use_cache=input_data.use_cache,
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

    # Classify quality from mean pLDDT
    mean_plddt = float(result.get("mean_plddt", 0.0))
    if mean_plddt >= 90:
        quality_class = "very_high"
    elif mean_plddt >= 70:
        quality_class = "high"
    elif mean_plddt >= 50:
        quality_class = "medium"
    else:
        quality_class = "low"

    return StructurePredictResponse(
        pdb_string=result["pdb_string"],
        mean_plddt=mean_plddt,
        plddt_scores=result.get("plddt_scores", []),
        quality_class=quality_class,
        execution_time_s=round(elapsed, 3),
    )


@router.post(
    "/structure/quality",
    response_model=QualityAssessResponse,
    summary="Assess structure quality from PDB",
)
async def structure_quality(input_data: QualityAssessInput):
    """
    Assess the quality of a protein structure from a PDB string.

    Evaluates pLDDT scores, Ramachandran outliers, clash score,
    and returns an overall quality assessment with a verdict.

    Returns 400 if PDB string is malformed.
    """
    try:
        from biocompiler.structure_quality import compute_structure_quality
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


# ─── Stability Endpoints ──────────────────────────────────────────────

@router.post(
    "/stability/analyze",
    response_model=StabilityResponse,
    summary="Analyze protein stability",
)
async def stability_analyze(input_data: StabilityInput):
    """
    Analyze protein stability using FoldX (empirical) and/or
    statistical potentials.

    If a PDB structure is provided, FoldX-based energy computation
    is used. Otherwise, a sequence-based statistical method is applied.

    Returns 503 if FoldX is not installed and no statistical fallback is
    available.
    """
    try:
        from biocompiler.foldx import empirical_stability
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="FoldX stability module not available. Install biocompiler[foldx].",
        )

    try:
        result = empirical_stability(
            protein=input_data.protein,
            organism=input_data.organism,
            pdb_string=input_data.pdb_string,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Stability analysis failed")
        raise HTTPException(
            status_code=500,
            detail=f"Stability analysis failed: {e}",
        )

    return StabilityResponse(
        stability_kcal=result["stability_kcal"],
        method=result["method"],
        components=result["components"],
        verdict=result["verdict"],
    )


@router.post(
    "/stability/mutations",
    response_model=MutationScanResponse,
    summary="Scan mutations for stability",
)
async def stability_mutations(input_data: MutationScanInput):
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
        from biocompiler.foldx import scan_mutations
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
            method=input_data.method,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Mutation scanning failed")
        raise HTTPException(
            status_code=500,
            detail=f"Mutation scanning failed: {e}",
        )

    return MutationScanResponse(
        mutations=result["mutations"],
        stabilizing_count=result["stabilizing_count"],
        destabilizing_count=result["destabilizing_count"],
    )


# ─── Solubility Endpoints ─────────────────────────────────────────────

@router.post(
    "/solubility/analyze",
    response_model=SolubilityResponse,
    summary="Analyze protein solubility",
)
async def solubility_analyze(input_data: SolubilityInput):
    """
    Analyze protein solubility using CamSol algorithm.

    Computes intrinsic solubility from sequence features and, if a PDB
    structure is provided, an overall solubility score incorporating
    structural corrections. Identifies aggregation-prone regions.

    Returns 400 for invalid protein or PDB input.
    """
    try:
        from biocompiler.camsol import compute_solubility
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

    return SolubilityResponse(
        intrinsic_score=result["intrinsic_score"],
        overall_score=result["overall_score"],
        solubility_class=result["solubility_class"],
        aggregation_prone_regions=result["aggregation_prone_regions"],
    )


@router.post(
    "/solubility/mutations",
    response_model=MutationScanResponse,
    summary="Find solubility-improving mutations",
)
async def solubility_mutations(input_data: SolubilityInput):
    """
    Find mutations that improve protein solubility.

    Scans positions in aggregation-prone regions and proposes
    amino acid substitutions that increase solubility score while
    preserving structural stability.

    Returns a list of proposed mutations with predicted solubility impact.
    """
    try:
        from biocompiler.camsol import find_solubility_mutations
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="CamSol solubility module not available.",
        )

    try:
        result = find_solubility_mutations(
            protein=input_data.protein,
            pdb_string=input_data.pdb_string,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Solubility mutation scan failed")
        raise HTTPException(
            status_code=500,
            detail=f"Solubility mutation scan failed: {e}",
        )

    return MutationScanResponse(
        mutations=result["mutations"],
        stabilizing_count=result.get("stabilizing_count", 0),
        destabilizing_count=result.get("destabilizing_count", 0),
    )


# ─── Immunogenicity Endpoints ─────────────────────────────────────────

@router.post(
    "/immunogenicity/analyze",
    response_model=ImmunogenicityResponse,
    summary="Analyze protein immunogenicity",
)
async def immunogenicity_analyze(input_data: ImmunogenicityInput):
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
        from biocompiler.immunogenicity import compute_immunogenicity
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="Immunogenicity module not available.",
        )

    try:
        result = compute_immunogenicity(
            protein=input_data.protein,
            organism=input_data.organism,
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

    return ImmunogenicityResponse(
        overall_score=result["overall_score"],
        immunogenicity_class=result["immunogenicity_class"],
        t_cell_epitopes=result["t_cell_epitopes"],
        b_cell_epitopes=result["b_cell_epitopes"],
        deimmunization_candidates=result["deimmunization_candidates"],
    )


@router.post(
    "/immunogenicity/deimmunize",
    response_model=DeimmunizeResponse,
    summary="Deimmunize a protein",
)
async def immunogenicity_deimmunize(input_data: DeimmunizeInput):
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
        from biocompiler.deimmunization import deimmunize
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

    return DeimmunizeResponse(
        optimized_protein=result["optimized_protein"],
        mutations_applied=result["mutations_applied"],
        original_score=result["original_score"],
        optimized_score=result["optimized_score"],
        success=result["success"],
    )


# ─── Full Assessment Endpoint ──────────────────────────────────────────

@router.post(
    "/assessment/full",
    response_model=FullAssessmentResponse,
    summary="Full protein assessment (all Phase 2+3)",
)
async def assessment_full(input_data: FullAssessmentInput):
    """
    Run a full protein assessment combining all Phase 2+3 analyses.

    Optionally includes structure prediction, stability analysis,
    solubility analysis, and immunogenicity assessment. Also runs
    Phase 1 predicates on the optimized sequence for completeness.

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
                from biocompiler.esmfold import predict_structure, is_esmfold_available

                if is_esmfold_available():
                    pred = predict_structure(
                        protein=input_data.protein,
                        organism=input_data.organism,
                        use_cache=True,
                    )
                    pdb_str = pred["pdb_string"]
                    mean_plddt = float(pred.get("mean_plddt", 0.0))

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
                from biocompiler.structure_quality import compute_structure_quality

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
            from biocompiler.foldx import empirical_stability

            stab = empirical_stability(
                protein=input_data.protein,
                organism=input_data.organism,
                pdb_string=input_data.pdb_string,
            )
            stability_result = stab

            if stab["verdict"] == "UNSTABLE":
                recommendations.append(
                    f"Protein is predicted to be unstable "
                    f"(ΔG={stab['stability_kcal']:.1f} kcal/mol). "
                    f"Consider stability-enhancing mutations."
                )
            elif stab["verdict"] == "MARGINAL":
                recommendations.append(
                    f"Protein stability is marginal "
                    f"(ΔG={stab['stability_kcal']:.1f} kcal/mol). "
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
            from biocompiler.camsol import compute_solubility

            sol = compute_solubility(
                protein=input_data.protein,
                pdb_string=input_data.pdb_string,
            )
            solubility_result = sol

            if sol["solubility_class"] in ("low", "very_low"):
                recommendations.append(
                    f"Protein solubility is {sol['solubility_class']} "
                    f"(score={sol['overall_score']:.2f}). "
                    f"Consider solubility-enhancing mutations or fusion tags."
                )
            if sol["aggregation_prone_regions"]:
                recommendations.append(
                    f"Found {len(sol['aggregation_prone_regions'])} "
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
            from biocompiler.immunogenicity import compute_immunogenicity

            imm = compute_immunogenicity(
                protein=input_data.protein,
                organism=input_data.organism,
                mhc_alleles=None,
            )
            immunogenicity_result = imm

            if imm["immunogenicity_class"] in ("high", "very_high"):
                recommendations.append(
                    f"Protein immunogenicity is {imm['immunogenicity_class']} "
                    f"(score={imm['overall_score']:.2f}). "
                    f"Consider deimmunization for therapeutic applications."
                )
            if imm["deimmunization_candidates"]:
                recommendations.append(
                    f"Found {len(imm['deimmunization_candidates'])} "
                    f"deimmunization candidate position(s)."
                )
        except ImportError:
            recommendations.append(
                "Immunogenicity module not available. Immunogenicity analysis skipped."
            )
        except Exception as e:
            logger.warning("Immunogenicity analysis failed: %s", e)
            recommendations.append(f"Immunogenicity analysis failed: {e}")

    # ── Phase 1 Predicate Results ──────────────────────────────────
    try:
        from biocompiler.structure_report import assess_protein

        report = assess_protein(
            protein=input_data.protein,
            organism=input_data.organism,
            pdb_string=input_data.pdb_string,
        )
        predicate_results = report.get("predicate_results", [])
    except ImportError:
        # Fallback: try Phase 1 type system directly
        try:
            from biocompiler.translation import translate, compute_cai
            from biocompiler.optimization import optimize_sequence

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
        except Exception:
            predicate_results = []
    except Exception as e:
        logger.warning("Phase 1 predicate evaluation failed: %s", e)
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


# ─── Registration Function ─────────────────────────────────────────────

def register_phase2_3_routes(app: Any) -> None:
    """
    Register all Phase 2+3 routes on an existing FastAPI application.

    Usage::

        from fastapi import FastAPI
        from biocompiler.api_phase2_3 import register_phase2_3_routes

        app = FastAPI()
        register_phase2_3_routes(app)

    This adds all /phase2_3/* endpoints to the application.
    """
    app.include_router(router)
    logger.info("Phase 2+3 routes registered on FastAPI app at /phase2_3/")
