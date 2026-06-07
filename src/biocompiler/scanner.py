"""Deprecated: use biocompiler.sequence.scanner instead."""
import warnings

warnings.warn(
    "biocompiler.scanner is deprecated — use biocompiler.sequence.scanner instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.scanner import *  # noqa: F401,F403

__all__ = [
    "validate_dna_sequence",
    "gc_content",
    "scan_sequence",
    "KOZAK_POSITION_WEIGHTS",
    "SPLICE_DONOR_MIN_SCORE",
    "SPLICE_ACCEPTOR_MIN_SCORE",
    "DONOR_FALLBACK_SCORE",
    "POLYPYRIMIDINE_MIN_FRACTION",
    "ACCEPTOR_SCORE_MULTIPLIER",
    "KOZAK_REPORT_THRESHOLD",
    "NUM_READING_FRAMES",
    "CODON_LENGTH",
    "DEFAULT_MOTIF_SCORE",
    "SCORE_ROUND_DIGITS",
    "KOZAK_UPSTREAM_CONTEXT",
    "KOZAK_DOWNSTREAM_CONTEXT",
    "_iupac_match",
    "_score_kozak",
]
