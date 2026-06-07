"""Deprecated: use biocompiler.export.sbol_import instead."""
import warnings

warnings.warn(
    "biocompiler.sbol_import is deprecated — use biocompiler.export.sbol_import instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.sbol_import import *  # noqa: F401,F403
