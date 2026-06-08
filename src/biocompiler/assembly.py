"""Deprecated: use biocompiler.optimizer.assembly instead."""
import warnings

warnings.warn(
    "biocompiler.assembly is deprecated — use biocompiler.optimizer.assembly instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.optimizer.assembly import *  # noqa: F401,F403
