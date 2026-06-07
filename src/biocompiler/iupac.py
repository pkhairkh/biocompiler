"""Deprecated: use biocompiler.sequence.iupac instead."""
import warnings

warnings.warn(
    "biocompiler.iupac is deprecated — use biocompiler.sequence.iupac instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.iupac import *  # noqa: F401,F403

__all__ = [
    "IUPAC_DNA",
    "resolve_ambiguous",
    "is_ambiguous",
    "expand_ambiguous",
    "has_ambiguous",
    "validate_iupac_sequence",
    "VALID_IUPAC_BASES",
]
