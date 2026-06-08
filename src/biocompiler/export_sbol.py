"""Deprecated: use biocompiler.export.sbol_legacy instead."""
import warnings

warnings.warn(
    "biocompiler.export_sbol is deprecated — use biocompiler.export.sbol_legacy instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.sbol_legacy import *  # noqa: F401,F403

__all__ = [
    "_build_rdf_xml",
    "_generate_sbol_identity",
]
