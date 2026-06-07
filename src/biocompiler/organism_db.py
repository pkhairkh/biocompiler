"""Deprecated: use biocompiler.organisms.db instead."""
import warnings

warnings.warn(
    "biocompiler.organism_db is deprecated — use biocompiler.organisms.db instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.organisms.db import *  # noqa: F401,F403
