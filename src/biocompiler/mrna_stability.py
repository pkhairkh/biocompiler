"""Deprecated: use biocompiler.expression.mrna_stability instead."""
import warnings

warnings.warn(
    "biocompiler.mrna_stability is deprecated — use biocompiler.expression.mrna_stability instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.mrna_stability import *  # noqa: F401,F403

__all__ = [
    "STABILITY_MOTIFS",
    "MRNAStabilityScore",
    "score_mrna_stability",
    "compute_mrna_half_life_score",
    "predict_mrna_stability",
    "suggest_mutations_for_stability",
]
