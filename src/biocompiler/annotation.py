"""Deprecated: use biocompiler.export.annotation instead."""
import warnings

warnings.warn(
    "biocompiler.annotation is deprecated — use biocompiler.export.annotation instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.export.annotation import *  # noqa: F401,F403

__all__ = [
    "SequenceAnnotation",
    "annotate_sequence",
    "annotate_to_genbank",
    "_find_orfs",
    "_find_restriction_sites",
    "_find_cpg_islands",
    "_find_gc_at_rich_regions",
    "_find_rbs",
    "_find_simple_repeats",
    "_find_splice_sites",
]
