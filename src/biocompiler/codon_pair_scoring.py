"""Deprecated: use biocompiler.expression.codon_pair_scoring instead."""
import warnings

warnings.warn(
    "biocompiler.codon_pair_scoring is deprecated — use biocompiler.expression.codon_pair_scoring instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.codon_pair_scoring import *  # noqa: F401,F403

__all__ = [
    "compute_cpb",
    "compute_cpb_score",
    "estimate_cpb_from_codon_freq",
    "get_codon_pair_data",
    "score_codon_pair",
    "suggest_better_pair",
]
