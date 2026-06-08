"""Deprecated: use biocompiler.export.genbank_roundtrip instead."""
import warnings

warnings.warn(
    "biocompiler.genbank_roundtrip is deprecated — use biocompiler.export.genbank_roundtrip instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.genbank_roundtrip import *  # noqa: F401,F403
