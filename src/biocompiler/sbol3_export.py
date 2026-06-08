"""Deprecated: use biocompiler.export.sbol3_export instead."""
import warnings

warnings.warn(
    "biocompiler.sbol3_export is deprecated — use biocompiler.export.sbol3_export instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.sbol3_export import *  # noqa: F401,F403

__all__ = [
    "_generate_identity",
    "BC_NS",
    "BIOPAX_NS",
    "DCT_NS",
    "DNAREGION_TYPE",
    "IUPAC_DNA_ENCODING",
    "RDF_NS",
    "SBOL3_NS",
    "SO_CDS",
    "SO_GENE",
]
