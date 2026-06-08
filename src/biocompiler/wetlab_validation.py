"""Deprecated: use biocompiler.validation.wetlab_validation instead."""
import warnings

warnings.warn(
    "biocompiler.wetlab_validation is deprecated — use biocompiler.validation.wetlab_validation instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.validation.wetlab_validation import *  # noqa: F401,F403
