"""Deprecated: use biocompiler.benchmarking.cai_benchmarking instead."""
import warnings

warnings.warn(
    "biocompiler.cai_benchmarking is deprecated — use biocompiler.benchmarking.cai_benchmarking instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.benchmarking.cai_benchmarking import *  # noqa: F401,F403

__all__ = [
    "BenchmarkSummary",
]
