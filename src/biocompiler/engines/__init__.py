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
)

__all__ = [
    "predict_structure_heuristic",
    "estimate_plddt_from_sequence",
    "estimate_secondary_structure_from_sequence",
    "compute_hydrophobicity_profile",
    "compute_charge_profile",
    "compute_contact_density",
    "HEURISTIC_MAX_CONFIDENCE",
]
