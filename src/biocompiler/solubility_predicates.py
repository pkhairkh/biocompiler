"""Deprecated: use biocompiler.type_system.solubility_predicates instead."""
import warnings

warnings.warn(
    "biocompiler.solubility_predicates is deprecated — use biocompiler.type_system.solubility_predicates instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.type_system.solubility_predicates import *  # noqa: F401,F403

__all__ = [
    "_CAMSOL_INTRINSIC",
    "_AGG_BORDERLINE_MAX",
    "_AGG_LIKELY_FAIL_MAX",
    "_AGG_UNCERTAIN_MAX",
    "_CAMSOL_HIGHLY_SOLUBLE",
    "_CAMSOL_MARGINAL",
    "_CAMSOL_SCALE",
    "_CAMSOL_WINDOW",
    "_DEFAULT_HYDROPHOBIC",
    "_HYDRO_EXCESS_BORDERLINE",
    "_HYDRO_EXCESS_UNCERTAIN",
    "_camsol_overall_score",
    "_camsol_smoothed_profile",
    "_find_aggregation_regions",
]
