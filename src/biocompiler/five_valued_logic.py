"""Deprecated: use biocompiler.shared.five_valued_logic instead."""
import warnings

warnings.warn(
    "biocompiler.five_valued_logic is deprecated — use biocompiler.shared.five_valued_logic instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.five_valued_logic import *  # noqa: F401,F403

__all__ = [
    "FiveValuedVerdict",
    "to_three_valued",
    "confidence_score",
    "combine_verdicts",
    "verify_five_valued_soundness",
    "five_valued_and",
    "five_valued_or",
    "refinement_is_sound",
]
