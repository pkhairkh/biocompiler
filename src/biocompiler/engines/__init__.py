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
    predict_mfe_fallback,
    predict_accessibility_fallback,
    find_stable_structures_fallback,
    MFEResult,
    AccessibilityResult,
    StableStructure,
    MIN_LOOP_LENGTH,
    PAIR_ENERGIES,
    DEFAULT_WINDOW_SIZE,
    DEFAULT_STABLE_DG_THRESHOLD,
)

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
    "predict_mfe_fallback",
    "predict_accessibility_fallback",
    "find_stable_structures_fallback",
    "MFEResult",
    "AccessibilityResult",
    "StableStructure",
    "MIN_LOOP_LENGTH",
    "PAIR_ENERGIES",
    "DEFAULT_WINDOW_SIZE",
    "DEFAULT_STABLE_DG_THRESHOLD",
]
