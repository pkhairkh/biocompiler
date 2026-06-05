"""BioCompiler Engines Subpackage.

Contains engine-specific modules and offline fallbacks.
"""

from .esmfold_fallback import (
    predict_structure_heuristic,
    estimate_plddt_from_sequence,
    estimate_secondary_structure_from_sequence,
    compute_hydrophobicity_profile,
    compute_charge_profile,
    compute_contact_density,
    HEURISTIC_MAX_CONFIDENCE,
    HEURISTIC_MIN_CONFIDENCE,
    ChargeProfile,
    SecondaryStructureEstimate,
    ContactDensityProfile,
)
from .viennarna_fallback import (
    nussinov_fold,
    compute_approx_dg,
    compute_gc_dg_estimate,
    compute_nntm_dg,
    predict_mfe_fallback,
    predict_accessibility_fallback,
    find_stable_structures_fallback,
)

# Re-export types from viennarna (or local fallbacks)
try:
    from ..viennarna import MFEResult, StemLoop, AccessibilityResult
except ImportError:
    from .viennarna_fallback import MFEResult, StemLoop, AccessibilityResult  # type: ignore[no-redef]

__all__ = [
    # esmfold_fallback
    "predict_structure_heuristic",
    "estimate_plddt_from_sequence",
    "estimate_secondary_structure_from_sequence",
    "compute_hydrophobicity_profile",
    "compute_charge_profile",
    "compute_contact_density",
    "HEURISTIC_MAX_CONFIDENCE",
    "HEURISTIC_MIN_CONFIDENCE",
    "ChargeProfile",
    "SecondaryStructureEstimate",
    "ContactDensityProfile",
    # viennarna_fallback
    "nussinov_fold",
    "compute_approx_dg",
    "compute_gc_dg_estimate",
    "compute_nntm_dg",
    "predict_mfe_fallback",
    "predict_accessibility_fallback",
    "find_stable_structures_fallback",
    # types
    "MFEResult",
    "StemLoop",
    "AccessibilityResult",
]
