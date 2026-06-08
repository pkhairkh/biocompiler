"""Deprecated: use biocompiler.sequence.sliding_gc instead."""
import warnings

warnings.warn(
    "biocompiler.sliding_gc is deprecated — use biocompiler.sequence.sliding_gc instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.sliding_gc import *  # noqa: F401,F403

__all__ = [
    "WindowViolation",
    "SlidingGCResult",
    "check_sliding_gc",
    "fix_sliding_gc_violations",
    "evaluate_sliding_gc",
    "_check_sliding_gc_python",
    "_FORCE_PYTHON_GC_WINDOW",
    "_HAS_NUMBA",
    "_check_sliding_gc_numba",
]
