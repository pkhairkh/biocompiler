"""Deprecated: use biocompiler.validation.dataset_validation instead."""
import warnings

warnings.warn(
    "biocompiler.dataset_validation is deprecated — use biocompiler.validation.dataset_validation instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.validation.dataset_validation import *  # noqa: F401,F403
