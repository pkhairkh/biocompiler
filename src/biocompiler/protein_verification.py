"""Deprecated: use biocompiler.validation.protein_verification instead."""
import warnings

warnings.warn(
    "biocompiler.protein_verification is deprecated — use biocompiler.validation.protein_verification instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.validation.protein_verification import *  # noqa: F401,F403
