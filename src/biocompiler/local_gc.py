"""Deprecated: use biocompiler.sequence.local_gc instead."""
import warnings

warnings.warn(
    "biocompiler.local_gc is deprecated — use biocompiler.sequence.local_gc instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.local_gc import *  # noqa: F401,F403

__all__ = [
    "LocalGCConstraint",
    "LocalGCResult",
    "check_local_gc",
    "optimize_local_gc",
]
