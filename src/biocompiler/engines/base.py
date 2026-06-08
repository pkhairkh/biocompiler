"""BioCompiler Engine Base — Unified base types for all analysis engines.

v9.2.0 — Complete API unification foundation

All analysis engines (ESMFold, FoldX, CamSol, Immunogenicity) share:
  - BaseEngineResult: concrete base class with unified field names
  - MutationResult: unified mutation suggestion type
  - BatchResult[T]: generic batch result
  - EngineTimer: timing context manager
  - EngineConfig: unified configuration
  - validate_protein_sequence(): input validation
  - classify_score(): unified score→classification helper

Design principles:
  1. Every engine result inherits from BaseEngineResult
  2. Unified field names: primary_score, classification, mutations
  3. Domain-specific names preserved as properties (plddt, ddg, etc.)
  4. All batch functions return BatchResult[SpecificResult]
  5. All mutation-finding functions return list[MutationResult]
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import (
    Any,
    Generic,
    Protocol,
    TypeVar,
    runtime_checkable,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Protein input validation
# ---------------------------------------------------------------------------

__all__ = [
    "STANDARD_AMINO_ACIDS",
    "DEFAULT_TIMEOUT_S",
    "DEFAULT_MAX_WORKERS",
    "DEFAULT_CONFIDENCE",
    "validate_protein_sequence",
    "EngineResult",
    "BaseEngineResult",
    "MutationResult",
    "BatchResult",
    "EngineTimer",
    "EngineConfig",
    "classify_score",
]

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: Default wall-clock timeout for engine operations (seconds).
DEFAULT_TIMEOUT_S: float = 300.0

#: Default parallelism for batch engine operations.
DEFAULT_MAX_WORKERS: int = 4

#: Default confidence score for mutation suggestions.
DEFAULT_CONFIDENCE: float = 1.0


def validate_protein_sequence(protein: str, engine_name: str = "engine") -> str:
    """Validate and normalize a protein sequence for engine use.

    Args:
        protein: amino acid sequence string
        engine_name: name of the engine (for error messages)

    Returns:
        Normalized (uppercase, whitespace-stripped) protein sequence

    Raises:
        ValueError: if the sequence is empty or contains invalid characters
    """
    if not protein or not protein.strip():
        raise ValueError(f"{engine_name}: protein sequence must not be empty")

    protein = protein.strip().upper()

    invalid = set(protein) - STANDARD_AMINO_ACIDS
    if invalid:
        raise ValueError(
            f"{engine_name}: protein sequence contains invalid amino acids: {sorted(invalid)}"
        )

    return protein


# ---------------------------------------------------------------------------
# Engine result protocol (structural typing)
# ---------------------------------------------------------------------------

@runtime_checkable
class EngineResult(Protocol):
    """Protocol that all engine result types must satisfy.

    Every analysis engine returns a result that includes:
    - success: whether the analysis completed without errors
    - error: error message if the analysis failed
    - execution_time_s: wall-clock time for the computation
    """
    success: bool
    error: str | None
    execution_time_s: float


# ---------------------------------------------------------------------------
# Base engine result (concrete base class)
# ---------------------------------------------------------------------------

@dataclass
class BaseEngineResult:
    """Concrete base class for all engine results.

    Provides unified field names that work across all engines:
      - primary_score: the main metric (aliases: plddt, ddg, score, immunogenicity_score)
      - classification: the categorical label (aliases: confidence_class, stability_class, etc.)
      - mutations: list of suggested mutations (aliases: stabilizing_mutations, etc.)

    Subclasses add engine-specific fields and property aliases.
    """
    sequence: str
    primary_score: float
    classification: str
    success: bool
    error: str | None = None
    execution_time_s: float = 0.0
    engine_name: str = ""
    primary_score_label: str = "score"

    @property
    def passed(self) -> bool:
        """Whether the analysis completed successfully."""
        return self.success


# ---------------------------------------------------------------------------
# Unified mutation result
# ---------------------------------------------------------------------------

@dataclass
class MutationResult:
    """Unified mutation suggestion from any engine.

    All engines that suggest mutations (FoldX, CamSol, Immunogenicity,
    Deimmunization) return this type.  Replaces the old ad-hoc Mutation
    class from mutations.py and the various per-engine mutation types.

    Fields:
      position: 0-based residue index in the protein sequence
      original: original amino acid (single letter)
      mutant: suggested replacement (single letter)
      delta_score: predicted change in the relevant metric
        (negative ΔΔG = stabilizing, positive solubility = improving,
         negative immunogenicity = deimmunizing)
      score_type: which metric the delta_score refers to
        ('ddg', 'solubility', 'immunogenicity')
      engine: which engine suggested this
        ('foldx', 'camsol', 'immunogenicity', 'deimmunization')
      recommendation: human-readable category
        ('stabilizing', 'solubility_improving', 'deimmunizing')
      confidence: confidence score (0.0-1.0) for this suggestion
      details: engine-specific extra data
    """
    position: int
    original: str
    mutant: str
    delta_score: float = 0.0
    score_type: str = ""  # 'ddg', 'solubility', 'immunogenicity'
    engine: str = ""  # 'foldx', 'camsol', 'immunogenicity', 'deimmunization'
    recommendation: str = ""  # 'stabilizing', 'solubility_improving', 'deimmunizing'
    description: str = ""
    confidence: float = DEFAULT_CONFIDENCE
    details: dict[str, Any] = field(default_factory=dict)

    def __init__(
        self,
        *,
        score: float | None = None,
        delta_score: float | None = None,
        position: int = 0,
        original: str = "",
        mutant: str = "",
        score_type: str = "",
        engine: str = "",
        recommendation: str = "",
        description: str = "",
        confidence: float = DEFAULT_CONFIDENCE,
        details: dict[str, Any] | None = None,
        # Old field name aliases for backward compat
        original_aa: str | None = None,
        mutant_aa: str | None = None,
        **kwargs: Any,
    ) -> None:
        """Initialize MutationResult with backward-compatible score= alias."""
        # Handle score alias: score= maps to delta_score
        effective_delta = 0.0
        if delta_score is not None:
            effective_delta = delta_score
        elif score is not None:
            effective_delta = score
        # Handle old field name aliases: original_aa/mutant_aa override original/mutant if original/mutant are empty
        effective_original = original_aa if (not original and original_aa) else original
        effective_mutant = mutant_aa if (not mutant and mutant_aa) else mutant
        # Set fields using object.__setattr__ to bypass dataclass restrictions
        self.position = position
        self.original = effective_original or ""
        self.mutant = effective_mutant or ""
        self.delta_score = effective_delta
        self.score_type = score_type
        self.engine = engine
        self.recommendation = recommendation
        self.description = description
        self.confidence = confidence
        self.details = details if details is not None else {}

        if kwargs:
            logger.warning(
                "MutationResult received unexpected keyword arguments: %s",
                ", ".join(sorted(kwargs.keys())),
            )

    # Backward compatibility: 'score' property alias for delta_score
    @property
    def score(self) -> float:
        """Alias for delta_score (backward compatibility)."""
        return self.delta_score

    @score.setter
    def score(self, value: float) -> None:
        """Set delta_score via the score alias."""
        self.delta_score = value

    # Old field name aliases
    @property
    def original_aa(self) -> str:
        """Alias for original (backward compatibility)."""
        return self.original

    @property
    def mutant_aa(self) -> str:
        """Alias for mutant (backward compatibility)."""
        return self.mutant

    def __str__(self) -> str:
        """Return human-readable mutation string (e.g. 'A15G')."""
        return (
            f"{self.original}{self.position + 1}{self.mutant} "
            f"({self.engine}: {self.score_type}={self.delta_score:.2f})"
        )


# ---------------------------------------------------------------------------
# Batch result
# ---------------------------------------------------------------------------

T = TypeVar("T")


@dataclass
class BatchResult(Generic[T]):
    """Unified result type for batch processing across all engines.

    Generic in the element type T so callers get proper typing:
        BatchResult[ESMFoldResult], BatchResult[FoldXResult], etc.
    """
    results: list[T] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    total_time_s: float = 0.0
    successful: int = 0
    failed: int = 0

    def __post_init__(self) -> None:
        """Auto-compute successful/failed counts when not explicitly provided."""
        if self.successful == 0 and self.failed == 0 and self.results:
            for r in self.results:
                if not hasattr(r, "success"):
                    logger.debug(
                        "BatchResult item %s lacks 'success' attribute; "
                        "assuming successful",
                        type(r).__name__,
                    )
            self.successful = sum(
                1 for r in self.results if getattr(r, "success", True)
            )
            self.failed = len(self.results) - self.successful

    @property
    def success_count(self) -> int:
        """Alias for successful — matches older code."""
        return self.successful

    @property
    def failure_count(self) -> int:
        """Alias for failed — matches older code."""
        return self.failed

    @property
    def total(self) -> int:
        """Total number of results (successful + failed)."""
        return len(self.results)


# ---------------------------------------------------------------------------
# Engine timer
# ---------------------------------------------------------------------------

class EngineTimer:
    """Context manager for tracking engine execution time.

    Usage:
        with EngineTimer() as timer:
            result = do_analysis(...)
        result.execution_time_s = timer.elapsed
    """

    def __init__(self) -> None:
        """Initialize timer with zero start and elapsed times."""
        self.start: float = 0.0
        self.elapsed: float = 0.0

    def __enter__(self) -> EngineTimer:
        """Start timing when entering context."""
        self.start = time.perf_counter()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: Any | None,
    ) -> None:
        """Stop timing when exiting context."""
        self.elapsed = time.perf_counter() - self.start


# ---------------------------------------------------------------------------
# Engine configuration
# ---------------------------------------------------------------------------

@dataclass
class EngineConfig:
    """Unified configuration for all analysis engines.

    Every engine function accepts these as keyword arguments:
        use_cache: whether to use result caching
        timeout_s: wall-clock timeout in seconds
        verbose: whether to log detailed progress
        max_workers: parallelism for batch operations
    """
    use_cache: bool = True
    timeout_s: float = DEFAULT_TIMEOUT_S
    verbose: bool = False
    max_workers: int = DEFAULT_MAX_WORKERS


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def classify_score(
    score: float,
    thresholds: list[tuple[float, str]],
    fallback: str = "unknown",
) -> str:
    """Classify a numeric score into a category.

    Args:
        score: the numeric value to classify
        thresholds: list of (threshold, label) pairs in descending order.
            The first threshold where score >= threshold wins.
        fallback: label if no threshold matches

    Returns:
        The classification label.

    Example:
        >>> classify_score(85.0, [(90, "very_high"), (70, "high"), (50, "medium")])
        'high'
    """
    for threshold, label in thresholds:
        if score >= threshold:
            return label
    return fallback
