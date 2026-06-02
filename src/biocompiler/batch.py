"""BioCompiler Batch Processing Module v7.6.0

Unified batch processing for all BioCompiler analysis engines.

Provides:
  - BatchProcessor: Configurable batch processor for running analyses
  - batch_analyze(): Convenience function for batch analysis
  - batch_compare(): Compare results across engines

Uses unified types from engine_base:
  - BatchResult[T]: Generic batch result container
  - EngineTimer: Timing context manager
  - EngineConfig: Unified configuration
  - BaseEngineResult: Base class for all engine results
"""

from __future__ import annotations

import logging
from typing import Callable, List, Optional, Sequence, TypeVar

from .engine_base import BaseEngineResult, BatchResult, EngineTimer, EngineConfig

logger = logging.getLogger(__name__)

T = TypeVar("T", bound=BaseEngineResult)


class BatchProcessor:
    """Configurable batch processor for running analyses across multiple sequences.

    Uses unified types from engine_base:
      - Returns BatchResult[T] from all batch operations
      - Uses EngineTimer for timing
      - Accepts EngineConfig for configuration

    Example:
        processor = BatchProcessor(config=EngineConfig(max_workers=2))
        result = processor.run(sequences, analyze_fn)
        print(f"Successful: {result.successful}, Failed: {result.failed}")
    """

    def __init__(self, config: Optional[EngineConfig] = None):
        self.config = config or EngineConfig()

    def run(
        self,
        sequences: Sequence[str],
        analyze_fn: Callable[[str], T],
        engine_name: str = "batch",
    ) -> BatchResult[T]:
        """Run batch analysis on multiple sequences.

        Args:
            sequences: Protein sequences to analyze
            analyze_fn: Analysis function that takes a sequence and returns T
            engine_name: Name of the engine for error messages

        Returns:
            BatchResult[T] containing all results, errors, and timing
        """
        results: List[T] = []
        errors: List[str] = []
        timer = EngineTimer()
        timer.__enter__()

        for seq in sequences:
            try:
                result = analyze_fn(seq)
                results.append(result)
            except Exception as exc:
                error_msg = f"{engine_name}: failed on sequence ({len(seq)} aa): {exc}"
                errors.append(error_msg)
                logger.error(error_msg)

        timer.__exit__(None, None, None)
        return BatchResult[T](
            results=results,
            errors=errors,
            total_time_s=round(timer.elapsed, 4),
        )


def batch_analyze(
    sequences: Sequence[str],
    analyze_fn: Callable[[str], T],
    engine_name: str = "batch",
    config: Optional[EngineConfig] = None,
) -> BatchResult[T]:
    """Convenience function for batch analysis.

    Analyzes multiple protein sequences using the provided analysis function
    and returns a unified BatchResult.

    Args:
        sequences: Protein sequences to analyze
        analyze_fn: Analysis function that takes a sequence and returns T
        engine_name: Name of the engine for error messages
        config: Optional configuration (max_workers, timeout, etc.)

    Returns:
        BatchResult[T] with results, errors, and timing information

    Example:
        from biocompiler.camsol import compute_solubility
        result = batch_analyze(proteins, compute_solubility, "camsol")
        print(f"Completed: {result.successful}/{result.total}")
    """
    processor = BatchProcessor(config=config)
    return processor.run(sequences, analyze_fn, engine_name=engine_name)


def batch_compare(
    sequences: Sequence[str],
    *analyze_fns: Callable[[str], BaseEngineResult],
    engine_names: Optional[List[str]] = None,
    config: Optional[EngineConfig] = None,
) -> dict[str, BatchResult[BaseEngineResult]]:
    """Compare results from multiple analysis engines on the same sequences.

    Args:
        sequences: Protein sequences to analyze
        *analyze_fns: Analysis functions to compare
        engine_names: Optional names for each analysis function
        config: Optional configuration

    Returns:
        Dict mapping engine name to BatchResult

    Example:
        results = batch_compare(
            proteins,
            compute_solubility, compute_immunogenicity,
            engine_names=["camsol", "immunogenicity"],
        )
    """
    if engine_names is None:
        engine_names = [f"engine_{i}" for i in range(len(analyze_fns))]

    output: dict[str, BatchResult[BaseEngineResult]] = {}
    for name, fn in zip(engine_names, analyze_fns):
        output[name] = batch_analyze(
            sequences, fn, engine_name=name, config=config,
        )
    return output
