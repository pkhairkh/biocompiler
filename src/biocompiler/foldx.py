"""Deprecated: use biocompiler.engines.foldx instead."""
import warnings

warnings.warn(
    "biocompiler.foldx is deprecated — use biocompiler.engines.foldx instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.engines.foldx import *  # noqa: F401,F403

__all__ = [
    "STANDARD_AAS",
]
