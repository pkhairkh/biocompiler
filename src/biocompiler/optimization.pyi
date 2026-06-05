"""Type stubs for biocompiler.optimization — public API surface."""

from __future__ import annotations

from typing import Any

from .type_system import PredicateResult
from .decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
)

# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

class OptimizationResult:
    """Result of optimizing a protein sequence."""
    sequence: str
    gc_content: float
    cai: float
    failed_predicates: list[str]
    predicate_results: list[PredicateResult]
    certificate_text: str
    protein: str
    fallback_used: bool
    satisfied_predicates: list[str]
    aa_substitutions: list[dict[str, Any]]
    mutagenesis_applied: bool
    mrna_stability_score: float | None
    destabilizing_motifs_removed: int
    stability_improvement: float | None
    provenance: Any
    codon_pair_bias: float | None
    suggested_5utr: str | None
    suggested_3utr: str | None
    utr_score_5: float | None
    utr_score_3: float | None
    decision_trail: OptimizationDecisionTrail | None

    def __post_init__(self) -> None: ...


class FullConstructResult:
    """Complete expression construct: 5' UTR + CDS + 3' UTR."""
    utr5: str
    cds: str
    utr3: str
    full_construct: str
    organism: str
    gc_content: float
    cai: float
    utr_score_5: float | None
    utr_score_3: float | None
    protein: str

    def __post_init__(self) -> None: ...


# ────────────────────────────────────────────────────────────
# Public functions
# ────────────────────────────────────────────────────────────

def optimize_sequence(
    target_protein: str | None = ...,
    organism: str = ...,
    species: str | None = ...,
    gc_lo: float = ...,
    gc_hi: float = ...,
    cai_threshold: float = ...,
    enzymes: list[str] | None = ...,
    strategy: str = ...,
    use_csp_solver: bool = ...,
    optimize_mrna_stability: bool = ...,
    seed: int | None = ...,
    include_utr: bool = ...,
    consider_codon_pair_bias: bool = ...,
    track_provenance: bool = ...,
    **kwargs: Any,
) -> OptimizationResult:
    """Optimize a protein sequence for expression in the target organism.

    Args:
        target_protein: Amino acid sequence (1-letter codes, no stop).
        organism: Target organism (canonical, short key, or alias).
        species: Alias for ``organism`` (deprecated).
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        enzymes: List of restriction enzyme names to avoid.
        strategy: Optimization strategy ('hybrid', 'constraint_first', 'cai_first').
        use_csp_solver: If True, try CSP/SMT solver before greedy.
        optimize_mrna_stability: If True, run mRNA stability improvement pass.
        seed: Deterministic seed for reproducibility.
        include_utr: If True, generate UTR suggestions.
        consider_codon_pair_bias: If True, run codon pair bias pass.
        track_provenance: If True, collect decision-level provenance.
        **kwargs: Additional arguments.

    Returns:
        OptimizationResult with optimized sequence and metrics.

    Raises:
        InvalidProteinError: If the protein contains invalid amino acid codes.
        UnsupportedOrganismError: If the organism is not supported.
    """
    ...


def protein_to_aa_list(protein: str) -> list[str]:
    """Convert protein string to list of amino acid codes."""
    ...


# ────────────────────────────────────────────────────────────
# BioOptimizer class (legacy optimizer)
# ────────────────────────────────────────────────────────────

class BioOptimizer:
    """Certified gene sequence optimizer with multi-step CAI-maximizing pipeline."""
    species: str
    organism_name: str
    species_cai: dict[str, float]
    enzymes: list[str]
    splice_low: float
    splice_high: float
    cpg_window: int
    cpg_threshold: float
    min_blosum: int
    min_cai: float
    avoid_gt: bool
    strategy: str
    optimize_mrna_stability: bool
    organism_domain: str
    is_prokaryote: bool

    def __init__(
        self,
        species: str = ...,
        enzymes: list[str] | None = ...,
        splice_low: float = ...,
        splice_high: float = ...,
        cpg_window: int = ...,
        cpg_threshold: float = ...,
        min_blosum: int = ...,
        min_cai: float = ...,
        avoid_gt: bool = ...,
        strategy: str = ...,
        optimize_mrna_stability: bool = ...,
        **kwargs: Any,
    ) -> None: ...

    def optimize(
        self,
        seq: str,
        strategy: str | None = ...,
    ) -> tuple[str, list[PredicateResult], str]:
        """Run the full optimization pipeline."""
        ...
