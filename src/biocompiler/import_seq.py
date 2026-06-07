"""Deprecated: use biocompiler.sequence.import_seq instead."""
import warnings

warnings.warn(
    "biocompiler.import_seq is deprecated — use biocompiler.sequence.import_seq instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.import_seq import *  # noqa: F401,F403

__all__ = [
    "import_fasta",
    "import_genbank",
    "import_sequence",
    "_parse_exon_boundaries",
    "_clean_qualifier_value",
    "_looks_like_path",
    "_resolve_input",
]
