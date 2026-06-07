"""Deprecated: use biocompiler.immunogenicity.deimmunization instead."""
import warnings

warnings.warn(
    "biocompiler.deimmunization is deprecated — use biocompiler.immunogenicity.deimmunization instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.immunogenicity.deimmunization import *  # noqa: F401,F403

__all__ = [
    "_estimate_ddg",
    "_detect_mhc_class",
    "_estimate_solubility_impact",
    "_filter_binder_epitopes",
    "_filter_strong_binder_epitopes",
    "_get_mhc_alleles",
    "_is_anchor_position",
    "_is_structurally_dangerous",
    "_organism_to_species",
]
