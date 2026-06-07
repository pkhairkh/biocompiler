"""Deprecated: use biocompiler.sequence.maxentscan instead."""
import warnings

warnings.warn(
    "biocompiler.maxentscan is deprecated — use biocompiler.sequence.maxentscan instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.maxentscan import *  # noqa: F401,F403

__all__ = [
    "BASE_TO_INDEX",
    "BG_PROB",
    "DONOR_PWM",
    "DONOR_PWM_SCORE",
    "ACCEPTOR_PWM",
    "ACCEPTOR_PWM_SCORE",
    "score_donor",
    "score_acceptor",
    "scan_splice_sites",
    "max_donor_score",
    "max_acceptor_score",
    "validate_against_published",
    "check_consistency",
    "CRYPTIC_SPLICE_THRESHOLD",
    "_EDGE_CASE_SCORE",
    "_IMPOSSIBLE_SCORE",
    "_ACCEPTOR_DOWNSTREAM",
    "_ACCEPTOR_UPSTREAM",
    "_DONOR_DOWNSTREAM",
    "_DONOR_UPSTREAM",
]
