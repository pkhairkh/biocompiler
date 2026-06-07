"""Deprecated: use biocompiler.organisms.config instead."""
import warnings

warnings.warn(
    "biocompiler.organism_config is deprecated — use biocompiler.organisms.config instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.organisms.config import *  # noqa: F401,F403
