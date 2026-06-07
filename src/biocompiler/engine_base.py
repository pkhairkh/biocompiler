"""Deprecated: use biocompiler.engines.base instead."""
import warnings

warnings.warn(
    "biocompiler.engine_base is deprecated — use biocompiler.engines.base instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.engines.base import *  # noqa: F401,F403
