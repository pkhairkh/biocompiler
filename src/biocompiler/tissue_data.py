"""Deprecated: use biocompiler.organisms.tissue_data instead."""
import warnings

warnings.warn(
    "biocompiler.tissue_data is deprecated — use biocompiler.organisms.tissue_data instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.organisms.tissue_data import *  # noqa: F401,F403
