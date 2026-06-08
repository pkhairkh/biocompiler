"""Deprecated: use biocompiler.validation.wetlab_validation instead."""
import warnings

warnings.warn(
    "biocompiler.wet_lab_validation is deprecated — use biocompiler.validation.wetlab_validation instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.validation.wetlab_validation import *  # noqa: F401,F403

__all__ = [
    "ExperimentalResult",
    "ValidationComparison",
    "WetLabValidator",
    "_global_validator",
]
