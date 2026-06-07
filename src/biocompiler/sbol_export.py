"""Deprecated: use biocompiler.export.sbol_export instead."""
import warnings

warnings.warn(
    "biocompiler.sbol_export is deprecated — use biocompiler.export.sbol_export instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.sbol_export import *  # noqa: F401,F403

__all__ = [
    "SBOL3_NS",
    "RDF_NS",
    "SO_CDS",
    "SO_NS",
    "SO_PROMOTER",
    "SO_TERMINATOR",
    "_resolve_role_uri",
]
