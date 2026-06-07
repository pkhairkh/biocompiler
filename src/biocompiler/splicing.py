"""Deprecated: use biocompiler.sequence.splicing instead."""
import warnings

warnings.warn(
    "biocompiler.splicing is deprecated — use biocompiler.sequence.splicing instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.splicing import *  # noqa: F401,F403

__all__ = [
    "maxent_score",
    "maxent_score_v2",
    "score_splice_sites",
    "compute_splice_isoforms",
    "_MAXENT_PWM",
    "_BASE_INDEX",
    "_MIN_CONTEXT_LEN",
    "_PWM_CONTEXT_LEN",
    "_PWM_DOWNSTREAM",
    "_PWM_UPSTREAM",
]
