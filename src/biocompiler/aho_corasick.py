"""Deprecated: use biocompiler.sequence.aho_corasick instead."""
import warnings

warnings.warn(
    "biocompiler.aho_corasick is deprecated — use biocompiler.sequence.aho_corasick instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.aho_corasick import *  # noqa: F401,F403

__all__ = [
    "AhoCorasickNode",
    "AhoCorasickScanner",
    "build_scanner_from_enzymes",
    "build_scanner_from_sites",
]
