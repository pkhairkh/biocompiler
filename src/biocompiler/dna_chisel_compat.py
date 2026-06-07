"""Deprecated: use biocompiler.shared.dna_chisel_compat instead."""
import warnings

warnings.warn(
    "biocompiler.dna_chisel_compat is deprecated — use biocompiler.shared.dna_chisel_compat instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.dna_chisel_compat import *  # noqa: F401,F403

__all__ = [
    "GC_ENFORCEMENT_WINDOW",
    "_DNA_CHISEL_AVAILABLE",
    "_build_initial_sequence",
    "CAI_COMPARISON_EPSILON",
    "DEFAULT_BACTERIAL_PROMOTER_LENGTH",
    "DEFAULT_HAIRPIN_BOOST",
    "DEFAULT_HAIRPIN_STEM_SIZE",
    "DEFAULT_KMER_UNIQUIFY_SIZE",
    "MAX_RESTRICTION_ENZYMES",
    "_AA_ONE_TO_THREE",
    "_CONSTRAINT_BUILDERS",
    "_DNA_CHISEL_CONSTRAINTS",
    "_ORGANISM_MAP",
    "_build_dna_chisel_spec",
    "_compute_comparative_summary",
    "_compute_winners",
    "_count_restriction_sites",
    "avoid_bacterial_promoter",
    "avoid_changes",
    "avoid_hairpins",
    "avoid_pattern",
    "enforce_gc_content",
    "enforce_gc_content_local",
    "enforce_sequence",
    "enforce_start_codon",
    "enforce_stop_codon",
    "enforce_translation",
    "uniquify_all_kmers",
]
