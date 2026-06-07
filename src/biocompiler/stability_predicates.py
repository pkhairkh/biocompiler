"""Deprecated: use biocompiler.type_system.stability_predicates instead."""
import warnings

warnings.warn(
    "biocompiler.stability_predicates is deprecated — use biocompiler.type_system.stability_predicates instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.type_system.stability_predicates import *  # noqa: F401,F403

__all__ = [
    "_parse_pdb_coords",
    "_BLOSUM62_DDG_FACTOR",
    "_BLOSUM62_UNKNOWN_SCORE",
    "_CLEARLY_UNSTABLE_DG",
    "_DISULFIDE_BOND_KCAL",
    "_DISULFIDE_CB_DIST_THRESHOLD",
    "_ENTROPY_PENALTY_COEFF",
    "_HYDRO_CONTRIBUTION_WEIGHT",
    "_HYDRO_FRAC_HI",
    "_HYDRO_FRAC_LO",
    "_HYDRO_PEAK_FRAC",
    "_PRO_GLY_CONFIDENCE_THRESHOLD",
    "_PRO_GLY_PENALTY_THRESHOLD",
    "_PRO_GLY_PENALTY_WEIGHT",
    "_SALT_BRIDGE_KCAL_PER_PAIR",
    "_euclidean",
    "_get_cb_coords",
]
