"""Deprecated: use biocompiler.export.genbank_annotations instead."""
import warnings

warnings.warn(
    "biocompiler.genbank_annotations is deprecated — use biocompiler.export.genbank_annotations instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.genbank_annotations import *  # noqa: F401,F403
