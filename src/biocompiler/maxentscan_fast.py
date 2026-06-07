"""Deprecated: use biocompiler.sequence.maxentscan_fast instead."""
import warnings

warnings.warn(
    "biocompiler.maxentscan_fast is deprecated — use biocompiler.sequence.maxentscan_fast instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.maxentscan_fast import *  # noqa: F401,F403

__all__ = [
    "HAS_NUMBA_MAXENT",
    "scan_splice_sites_fast",
    "scan_splice_sites_fast_str",
    "score_donor_fast",
    "score_acceptor_fast",
    "check_consistency",
    "score_donor",
    "score_acceptor",
    "scan_splice_sites",
    "max_donor_score",
    "max_acceptor_score",
    "validate_against_published",
    "CRYPTIC_SPLICE_THRESHOLD",
    "_DONOR_PWM_NP",
    "_ACCEPTOR_PWM_NP",
    "_BASE_IDX_NP",
    "_DONOR_LOG_NP",
    "_ACCEPTOR_LOG_NP",
]
