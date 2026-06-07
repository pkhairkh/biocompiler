"""Deprecated: use biocompiler.immunogenicity.mhcflurry_adapter instead."""
import warnings

warnings.warn(
    "biocompiler.mhcflurry_adapter is deprecated — use biocompiler.immunogenicity.mhcflurry_adapter instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.immunogenicity.mhcflurry_adapter import *  # noqa: F401,F403

__all__ = [
    "MHCBindingResult",
    "_LRUCache",
    "_extract_overlapping_peptides",
    "_mhcflurry_result_to_binding_result",
    "_validate_peptide",
    "_validate_protein",
    "_validate_sequence",
    "ic50_to_binding_score",
    "classify_binding",
]
