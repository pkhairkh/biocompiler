"""Deprecated: use biocompiler.provenance.reporting instead."""
import warnings

warnings.warn(
    "biocompiler.provenance_reporting is deprecated — use biocompiler.provenance.reporting instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance.reporting import *  # noqa: F401,F403

__all__ = [
    "ProvenanceQuery",
    "ProvenanceReport",
    "explain_position",
]
