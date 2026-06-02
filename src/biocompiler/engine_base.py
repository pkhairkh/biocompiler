"""BioCompiler Engine Base — Unified base types for analysis engines.

Provides common protocols and dataclasses that all analysis engines
(ESMFold, FoldX, CamSol, Immunogenicity) implement consistently.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Protocol, Sequence, runtime_checkable


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


@dataclass
class MutationResult:
    """Unified mutation suggestion from any engine.

    All engines that suggest mutations (FoldX, CamSol, Immunogenicity)
    return this type instead of ad-hoc dicts.
    """
    position: int  # 0-based residue index
    original: str  # original amino acid (single letter)
    mutant: str    # suggested replacement (single letter)
    score: float   # engine-specific score (higher = better improvement)
    engine: str    # which engine suggested this ("foldx", "camsol", "immunogenicity")
    description: str = ""  # human-readable description
    details: dict = field(default_factory=dict)  # engine-specific extra data

    def __str__(self) -> str:
        return f"{self.original}{self.position+1}{self.mutant} ({self.engine}: {self.score:.2f})"


@dataclass
class BatchResult:
    """Unified result type for batch processing across all engines.

    All engines that support batch processing return this type.
    """
    results: list  # list of engine-specific result objects
    errors: list[str] = field(default_factory=list)
    total_time_s: float = 0.0
    successful: int = 0
    failed: int = 0

    def __post_init__(self):
        if self.successful == 0 and self.failed == 0 and self.results:
            self.successful = sum(1 for r in self.results if getattr(r, 'success', True))
            self.failed = len(self.results) - self.successful


class EngineTimer:
    """Context manager for tracking engine execution time.

    Usage:
        with EngineTimer() as timer:
            result = do_analysis(...)
        result.execution_time_s = timer.elapsed
    """
    def __init__(self):
        self.start: float = 0.0
        self.elapsed: float = 0.0

    def __enter__(self) -> EngineTimer:
        self.start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.elapsed = time.perf_counter() - self.start


def validate_protein_sequence(protein: str, engine_name: str) -> str:
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

    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(protein) - valid_aas
    if invalid:
        raise ValueError(
            f"{engine_name}: protein sequence contains invalid amino acids: {sorted(invalid)}"
        )

    return protein
