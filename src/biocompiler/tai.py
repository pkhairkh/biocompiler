"""Deprecated: use biocompiler.expression.tai instead."""
import warnings

warnings.warn(
    "biocompiler.tai is deprecated — use biocompiler.expression.tai instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.tai import *  # noqa: F401,F403

__all__ = [
    "compute_tai",
    "calculate_tai",
    "compute_tai_and_cai",
    "optimize_for_tai",
    "TRNA_GENE_COPIES",
    "WOBBLE_RULES",
    "WOBBLE_EFFICIENCY",
    "SUPPORTED_ORGANISMS_TAI",
    "compute_codon_weights",
]
