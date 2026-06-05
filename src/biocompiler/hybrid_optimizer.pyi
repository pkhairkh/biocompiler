"""Type stubs for biocompiler.hybrid_optimizer — public API surface."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class Violation:
    """A single constraint violation with severity score."""
    violation_type: str
    position: int
    severity: float
    codon_indices: list[int]
    details: str = ""

    def __lt__(self, other: Violation) -> bool: ...


@dataclass
class HybridResult:
    """Result from the hybrid optimizer."""
    sequence: str
    cai: float
    gc_content: float
    violations_fixed: int = 0
    hill_climb_improvements: int = 0
    iterations_used: int = 0
    phase1_cai: float = 0.0
    phase2_cai: float = 0.0
    phase3_cai: float = 0.0
    warnings: list[str] = field(default_factory=list)


# ────────────────────────────────────────────────────────────
# HybridOptimizer class
# ────────────────────────────────────────────────────────────

class HybridOptimizer:
    """Hybrid gene optimizer combining greedy initialization with
    priority-based local search and CAI hill climbing.

    Architecture:
    1. Phase 1: Greedy CAI maximization (best codon per position)
    2. Phase 2: Priority-based constraint satisfaction
    3. Phase 3: CAI hill climbing (upgrade codons while maintaining constraints)
    """

    species: str
    organism: str
    species_cai: dict[str, float]
    enzymes: list[str]
    gc_lo: float
    gc_hi: float
    avoid_gt: bool
    splice_threshold: float
    cpg_window: int
    cpg_threshold: float
    max_local_search_iterations: int
    max_hill_climb_iterations: int
    cai_weight: float
    provenance_collector: Any
    is_prokaryote: bool
    sorted_codons: dict[str, list[str]]
    optimal_codon: dict[str, str]
    gt_free: dict[str, list[str]]
    ag_free: dict[str, list[str]]
    codon_gc: dict[str, int]

    def __init__(
        self,
        species: str = ...,
        organism: str | None = ...,
        enzymes: list[str] | None = ...,
        gc_lo: float = ...,
        gc_hi: float = ...,
        avoid_gt: bool = ...,
        splice_threshold: float = ...,
        cpg_window: int = ...,
        cpg_threshold: float = ...,
        max_local_search_iterations: int = ...,
        max_hill_climb_iterations: int = ...,
        cai_weight: float = ...,
        provenance_collector: Any = ...,
    ) -> None: ...

    def optimize(
        self,
        protein: str,
        is_prokaryote: bool = ...,
    ) -> HybridResult:
        """Run the hybrid optimization pipeline.

        Args:
            protein: Amino acid sequence (single-letter codes, no stop).
            is_prokaryote: When True, skip eukaryote-specific constraint steps.

        Returns:
            HybridResult with the optimized DNA sequence and metrics.
        """
        ...
