"""Deprecated: use biocompiler.engines.viennarna instead."""
import warnings

warnings.warn(
    "biocompiler.viennarna is deprecated — use biocompiler.engines.viennarna instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.engines.viennarna import *  # noqa: F401,F403

__all__ = [
    "DEFAULT_5PRIME_WINDOW",
    "DEFAULT_FOLD_TIMEOUT_SECONDS",
    "DEFAULT_FULL_LENGTH_CUTOFF",
    "DEFAULT_OVERLAP_THRESHOLD",
    "DEFAULT_STEP",
    "DEFAULT_WINDOW_SIZE",
    "EXPECTED_VIENNARNA_VERSION",
    "MAX_FOLD_TIMEOUT_SECONDS",
    "NEAREST_NEIGHBOR_AU",
    "NEAREST_NEIGHBOR_GC",
    "NEAREST_NEIGHBOR_GU",
]
