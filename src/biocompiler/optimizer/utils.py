"""
Utility functions and data classes for the optimizer.

Contains convergence tracking, result dataclasses, and protein validation helpers.
"""

from typing import Any

import logging
from dataclasses import dataclass, field

from ..type_system import AA_TO_CODONS, PredicateResult
from ..decision_provenance import OptimizationDecisionTrail
from ..exceptions import InvalidProteinError

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Named Constants
# ────────────────────────────────────────────────────────────

# Iteration limits for each optimization step
MAX_RESTRICTION_SITE_ITERATIONS: int = 100
MAX_IUPAC_SITE_ITERATIONS: int = 100
MAX_ATTTA_MOTIF_ITERATIONS: int = 100
MAX_T_RUN_ITERATIONS: int = 100
MAX_GC_ADJUSTMENT_ITERATIONS: int = 200
MAX_SPLICE_ELIMINATION_ITERATIONS: int = 300
MAX_CPG_DISRUPTION_ITERATIONS: int = 200

# Main optimization loop convergence settings
DEFAULT_MAX_ITERATIONS: int = 1000
CONVERGENCE_IMPROVEMENT_THRESHOLD: float = 1e-6
CONVERGENCE_PATIENCE: int = 3  # consecutive non-improving iterations before convergence
OSCILLATION_WINDOW: int = 10   # window size for oscillation detection

# Thresholds and sentinel values
# NOTE (Task 1.8): These are now imported from .constants via aliased imports
# at the top of this file.  The original names (T_RUN_LENGTH_THRESHOLD,
# ELIMINATED_SITE_SCORE, TOP_CAI_ALTERNATIVES, IUPAC_EXPANSION_CAP) are
# preserved as import aliases, so all downstream code continues to work.


# ────────────────────────────────────────────────────────────
# Convergence Tracker
# ────────────────────────────────────────────────────────────

class ConvergenceTracker:
    """Track optimization objective convergence and detect stagnation or oscillation.

    The objective is CAI * constraint_satisfaction_score.  CAI is in [0, 1]
    and constraint_satisfaction_score is the fraction of predicates that pass,
    also in [0, 1], so the combined objective is in [0, 1].

    Convergence is declared when the objective has not improved by more than
    ``improvement_threshold`` for ``patience`` consecutive iterations.

    Oscillation is declared when the objective has both increased and decreased
    within the last ``oscillation_window`` iterations — this indicates the
    optimizer is cycling and should stop at the best point seen so far.
    """

    def __init__(
        self,
        improvement_threshold: float = CONVERGENCE_IMPROVEMENT_THRESHOLD,
        patience: int = CONVERGENCE_PATIENCE,
        oscillation_window: int = OSCILLATION_WINDOW,
    ) -> None:
        self.improvement_threshold = improvement_threshold
        self.patience = patience
        self.oscillation_window = oscillation_window
        self.history: list[float] = []
        self.best_objective: float = float('-inf')
        self.best_iteration: int = -1
        self._no_improvement_count: int = 0

    def record(self, objective: float) -> None:
        """Record the objective value for the current iteration."""
        self.history.append(objective)
        if objective > self.best_objective:
            self.best_objective = objective
            self.best_iteration = len(self.history) - 1

    def check_convergence(self) -> str | None:
        """Check if the optimizer has converged, hit max iterations, or is oscillating.

        Returns:
            None if optimization should continue,
            "converged" if the objective has plateaued,
            "oscillating" if the objective is cycling,
        """
        if len(self.history) < 2:
            return None

        latest = self.history[-1]
        # Check plateau: no significant improvement for `patience` iterations
        recent = self.history[-(self.patience + 1):]
        if len(recent) >= self.patience + 1:
            improvement = max(recent) - min(recent)
            if improvement < self.improvement_threshold:
                return "converged"

        # Check oscillation: both increases and decreases within the window
        window = self.history[-self.oscillation_window:]
        if len(window) >= 3:
            has_increase = any(window[i + 1] > window[i] for i in range(len(window) - 1))
            has_decrease = any(window[i + 1] < window[i] for i in range(len(window) - 1))
            if has_increase and has_decrease:
                # Additional check: the objective hasn't improved overall in the window
                if window[-1] <= window[0] + self.improvement_threshold:
                    return "oscillating"

        return None

    @property
    def iterations(self) -> int:
        """Number of iterations recorded."""
        return len(self.history)

    @property
    def best(self) -> float:
        """Best objective value seen so far."""
        return self.best_objective

    @property
    def best_iteration_index(self) -> int:
        """Index of the iteration where the best objective was achieved."""
        return self.best_iteration


# ────────────────────────────────────────────────────────────
# High-level OptimizationResult and optimize_sequence API
# ────────────────────────────────────────────────────────────

@dataclass
class OptimizationResult:
    """Result of optimizing a protein sequence.

    Provides the optimized DNA sequence along with quality metrics
    and a list of predicates that the result fails to satisfy.
    """
    sequence: str
    gc_content: float
    cai: float
    failed_predicates: list[str] = field(default_factory=list)
    predicate_results: list[PredicateResult] = field(default_factory=list)
    certificate_text: str = ""
    # Extended attributes for API/visualization compatibility
    protein: str = ""
    fallback_used: bool = False
    satisfied_predicates: list[str] = field(default_factory=list)
    aa_substitutions: list[dict[str, Any]] = field(default_factory=list)
    mutagenesis_applied: bool = False
    # mRNA stability metrics (populated when optimize_mrna_stability=True)
    mrna_stability_score: float | None = None
    destabilizing_motifs_removed: int = 0
    stability_improvement: float | None = None
    # Provenance: OptimizationRecord for this run (populated by optimize_sequence)
    provenance: Any = field(default=None, repr=False)
    # Codon pair bias metric (populated when consider_codon_pair_bias=True)
    codon_pair_bias: float | None = None
    # UTR suggestions (populated when include_utr=True)
    suggested_5utr: str | None = None
    suggested_3utr: str | None = None
    utr_score_5: float | None = None
    utr_score_3: float | None = None
    # Decision-level provenance trail (populated when track_provenance=True)
    decision_trail: OptimizationDecisionTrail | None = None
    # Convergence tracking (populated by optimize_sequence)
    convergence_status: str | None = None  # "converged" | "max_iterations" | "oscillating" | None
    iterations_used: int = 0
    # Optimization warnings (cap-exceeded notices, convergence issues)
    warnings: list[str] = field(default_factory=list)
    # Custom objective score (populated when a non-default objective is used)
    objective_score: float | None = None

    def __post_init__(self):
        """Validate OptimizationResult invariants."""
        if self.protein and self.sequence:
            if len(self.sequence) != len(self.protein) * 3:
                raise ValueError(
                    f"Sequence length ({len(self.sequence)}) must equal "
                    f"protein length * 3 ({len(self.protein) * 3})"
                )
        if not (0.0 <= self.cai <= 1.0):
            raise ValueError(f"CAI must be in [0, 1], got {self.cai}")
        if not (0.0 <= self.gc_content <= 1.0):
            raise ValueError(f"GC content must be in [0, 1], got {self.gc_content}")
        if self.mutagenesis_applied:
            if self.aa_substitutions is None or len(self.aa_substitutions) == 0:
                raise ValueError(
                    "Mutagenesis applied but no substitutions recorded"
                )
        if self.utr_score_5 is not None:
            if not (0.0 <= self.utr_score_5 <= 1.0):
                raise ValueError(
                    f"UTR 5' score must be in [0, 1], got {self.utr_score_5}"
                )
        if self.utr_score_3 is not None:
            if not (0.0 <= self.utr_score_3 <= 1.0):
                raise ValueError(
                    f"UTR 3' score must be in [0, 1], got {self.utr_score_3}"
                )


@dataclass
class FullConstructResult:
    """Complete expression construct: 5' UTR + CDS + 3' UTR.

    This represents what a biologist would actually order from a gene
    synthesis company — the full DNA construct ready for cloning or
    direct expression in the target organism.

    The CDS is the optimized coding sequence. The UTRs are suggested
    (not enforced) and should be evaluated by the user before ordering.

    Attributes:
        utr5: 5' UTR sequence (empty string if not provided).
        cds: Optimized coding sequence (starts with ATG, ends with stop codon).
        utr3: 3' UTR sequence (empty string if not provided).
        full_construct: Concatenated 5'UTR + CDS + 3'UTR.
        organism: Target organism for this construct.
        gc_content: GC fraction of the full construct.
        cai: Codon Adaptation Index of the CDS.
        utr_score_5: Expression suitability score for the 5' UTR (0.0–1.0).
        utr_score_3: Expression suitability score for the 3' UTR (0.0–1.0).
        protein: Amino acid sequence encoded by the CDS.
    """
    utr5: str
    cds: str
    utr3: str
    full_construct: str
    organism: str
    gc_content: float
    cai: float
    utr_score_5: float | None = None
    utr_score_3: float | None = None
    protein: str = ""

    def __post_init__(self):
        """Validate FullConstructResult invariants."""
        if self.full_construct != self.utr5 + self.cds + self.utr3:
            raise ValueError(
                "full_construct must equal utr5 + cds + utr3"
            )
        if not (0.0 <= self.gc_content <= 1.0):
            raise ValueError(
                f"GC content must be in [0, 1], got {self.gc_content}"
            )
        if not (0.0 <= self.cai <= 1.0):
            raise ValueError(f"CAI must be in [0, 1], got {self.cai}")
        if self.utr_score_5 is not None:
            if not (0.0 <= self.utr_score_5 <= 1.0):
                raise ValueError(
                    f"UTR 5' score must be in [0, 1], got {self.utr_score_5}"
                )
        if self.utr_score_3 is not None:
            if not (0.0 <= self.utr_score_3 <= 1.0):
                raise ValueError(
                    f"UTR 3' score must be in [0, 1], got {self.utr_score_3}"
                )


# ==============================================================================
# Input Validation
# ==============================================================================

def protein_to_aa_list(protein: str) -> list[str]:
    """Convert protein string to list of amino acid codes. Raises InvalidProteinError for bad input.

    Pre-conditions:
    - protein must be a non-empty string of standard amino acid codes

    Post-conditions:
    - result is a list of valid single-letter amino acid codes
    - len(result) == len(protein.strip())
    """
    if not protein or not protein.strip():
        raise InvalidProteinError(protein, set())
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    invalid = set(ch for ch in protein if ch not in valid_aas)
    if invalid:
        raise InvalidProteinError(protein, invalid)
    return list(protein)


__all__ = [
    "ConvergenceTracker",
    "OptimizationResult",
    "FullConstructResult",
    "protein_to_aa_list",
    # Named constants
    "MAX_RESTRICTION_SITE_ITERATIONS",
    "MAX_IUPAC_SITE_ITERATIONS",
    "MAX_ATTTA_MOTIF_ITERATIONS",
    "MAX_T_RUN_ITERATIONS",
    "MAX_GC_ADJUSTMENT_ITERATIONS",
    "MAX_SPLICE_ELIMINATION_ITERATIONS",
    "MAX_CPG_DISRUPTION_ITERATIONS",
    "DEFAULT_MAX_ITERATIONS",
    "CONVERGENCE_IMPROVEMENT_THRESHOLD",
    "CONVERGENCE_PATIENCE",
    "OSCILLATION_WINDOW",
]
