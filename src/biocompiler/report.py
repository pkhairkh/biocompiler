"""Deprecated: use biocompiler.provenance.report instead."""
import warnings

warnings.warn(
    "biocompiler.report is deprecated — use biocompiler.provenance.report instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance.report import *  # noqa: F401,F403

__all__ = [
    "generate_report",
    "_build_organism_aware_constraints",
    "_build_cai_comparison_section",
    "_EUKARYOTE_ONLY_CONSTRAINTS",
    "_compute_gc_windows",
    "_DEFAULT_GC_WINDOW_SIZE",
    "_DESIGN_ID_DISPLAY_LENGTH",
    "_MAX_ISOFORM_DISPLAY",
    "_MAX_TOKEN_DISPLAY",
    "_VERDICT_SYMBOLS",
    "_compute_codon_usage",
    "_element_type_description",
]
