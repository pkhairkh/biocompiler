"""Deprecated: use biocompiler.export.core instead."""
import warnings

warnings.warn(
    "biocompiler.export is deprecated — use biocompiler.export.core instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.core import *  # noqa: F401,F403
