"""Deprecated: use biocompiler.expression.translation instead."""
import warnings

warnings.warn(
    "biocompiler.translation is deprecated — use biocompiler.expression.translation instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.translation import *  # noqa: F401,F403

__all__ = [
    "translate",
    "translate_with_confidence",
    "compute_cai",
    "find_orfs",
    "ORFResult",
    "DEFAULT_MIN_ORF_LENGTH_AA",
    "BACTERIAL_START_CODONS",
    "STANDARD_START_CODON",
    "PartialCodonError",
]
