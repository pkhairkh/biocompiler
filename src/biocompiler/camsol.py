"""Deprecated: use biocompiler.engines.camsol instead."""
import warnings

warnings.warn(
    "biocompiler.camsol is deprecated — use biocompiler.engines.camsol instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.engines.camsol import *  # noqa: F401,F403
