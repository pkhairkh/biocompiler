"""Deprecated: use biocompiler.engines.protein_design instead."""
import warnings

warnings.warn(
    "biocompiler.protein_design is deprecated — use biocompiler.engines.protein_design instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.engines.protein_design import *  # noqa: F401,F403

__all__ = [
    "_base_immunogenicity",
    "_estimate_ddg",
    "_base_solubility",
    "_base_stability",
    "_check_constraints",
    "_estimate_immunogenicity_delta",
    "_estimate_solubility_delta",
    "_is_preserved",
    "_predict_secondary_structure_simple",
]
