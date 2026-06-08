"""Deprecated: use biocompiler.optimizer.parts instead."""
import warnings

warnings.warn(
    "biocompiler.parts is deprecated — use biocompiler.optimizer.parts instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.optimizer.parts import *  # noqa: F401,F403
