"""Deprecated: use biocompiler.expression.sharp_li_tables instead."""
import warnings

warnings.warn(
    "biocompiler.sharp_li_tables is deprecated — use biocompiler.expression.sharp_li_tables instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.sharp_li_tables import *  # noqa: F401,F403

__all__ = [
    # Thread-safe accessors
    "get_sharp_li_table",
    "set_sharp_li_table",
    # State class (for advanced/testing use)
    "_SharpLiState",
    # Re-exported constants (backward compat)
    "_CODON_TABLE",
    "_AA_TO_CODONS",
    "_STOP_CODONS",
    "ECOLI_SHARP_LI_REFERENCE_GENES",
    "ECOLI_SHARP_LI_CODON_USAGE",
    "ECOLI_SHARP_LI_CAI_WEIGHTS",
    "YEAST_SHARP_LI_REFERENCE_GENES",
    "YEAST_SHARP_LI_CODON_USAGE",
    "YEAST_SHARP_LI_CAI_WEIGHTS",
    "SHARP_LI_PUBLISHED_CAI",
    "_REFERENCE_WEIGHTS",
    "SHARP_LI_REFERENCE_GENES",
    "SHARP_LI_CODON_USAGE",
    "SHARP_LI_CAI_WEIGHTS",
    # Re-exported functions
    "compute_cai_with_reference",
    "get_sharp_li_cai_weights",
]
