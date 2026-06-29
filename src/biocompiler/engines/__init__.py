"""BioCompiler Engines Subpackage.

Contains engine-specific modules and offline fallbacks.
"""

from .esmfold_fallback import (
    predict_structure_heuristic,
    estimate_plddt_from_sequence,
    estimate_secondary_structure_from_sequence,
    estimate_fold_quality,
    compute_hydrophobicity_profile,
    compute_charge_profile,
    compute_contact_density,
    HEURISTIC_MAX_CONFIDENCE,
    HEURISTIC_MIN_CONFIDENCE,
    ChargeProfile,
    SecondaryStructureEstimate,
    ContactDensityProfile,
    FoldQualityEstimate,
)
from biocompiler.engines.viennarna_fallback import (
    nussinov_fold,
    compute_approx_dg,
    compute_gc_dg_estimate,
    compute_nntm_dg,
    predict_mfe_fallback,
    predict_accessibility_fallback,
    find_stable_structures_fallback,
)
from .base import *  # noqa: F401,F403
from biocompiler.engines.viennarna import *  # noqa: F401,F403
from biocompiler.engines.esmfold import *  # noqa: F401,F403
from biocompiler.engines.foldx import *  # noqa: F401,F403
from biocompiler.engines.camsol import *  # noqa: F401,F403
from biocompiler.engines.protein_design import *  # noqa: F401,F403

# Re-export types from viennarna (or local fallbacks).
# Use __getattr__ to avoid circular imports at module load time: if
# viennarna.py ever imports from engines (directly or transitively),
# the lazy lookup below will still resolve correctly because by the
# time any attribute is actually accessed the module graph is fully
# initialised.
_LAZY_TYPE_NAMES = ("MFEResult", "StemLoop", "AccessibilityResult")


def __getattr__(name: str):  # type: ignore[no-untyped-def]
    if name in _LAZY_TYPE_NAMES:
        try:
            import biocompiler.engines.viennarna as _vr
            return getattr(_vr, name)
        except (ImportError, AttributeError):
            from biocompiler.engines.viennarna_fallback import MFEResult, StemLoop, AccessibilityResult
            return {"MFEResult": MFEResult, "StemLoop": StemLoop,
                    "AccessibilityResult": AccessibilityResult}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    # esmfold_fallback
    "predict_structure_heuristic",
    "estimate_plddt_from_sequence",
    "estimate_secondary_structure_from_sequence",
    "estimate_fold_quality",
    "compute_hydrophobicity_profile",
    "compute_charge_profile",
    "compute_contact_density",
    "HEURISTIC_MAX_CONFIDENCE",
    "HEURISTIC_MIN_CONFIDENCE",
    "ChargeProfile",
    "SecondaryStructureEstimate",
    "ContactDensityProfile",
    "FoldQualityEstimate",
    # viennarna_fallback
    "nussinov_fold",
    "compute_approx_dg",
    "compute_gc_dg_estimate",
    "compute_nntm_dg",
    "predict_mfe_fallback",
    "predict_accessibility_fallback",
    "find_stable_structures_fallback",
    # types (resolved lazily via __getattr__)
    "MFEResult",
    "StemLoop",
    "AccessibilityResult",
]
