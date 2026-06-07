"""Deprecated: use biocompiler.benchmarking.core instead."""
import warnings

warnings.warn(
    "biocompiler.benchmark is deprecated — use biocompiler.benchmarking.core instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.benchmarking.core import *  # noqa: F401,F403

__all__ = [
    "_build_best_codon_sequence",
    "_compute_metrics",
    "_compute_cai",
    "_count_cpg_ratio",
    "_count_gt",
    "_count_restriction_sites",
]
