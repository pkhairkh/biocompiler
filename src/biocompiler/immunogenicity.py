"""Deprecated: use biocompiler.immunogenicity.core instead."""
import warnings

warnings.warn(
    "biocompiler.immunogenicity is deprecated — use biocompiler.immunogenicity.core instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.immunogenicity.core import *  # noqa: F401,F403
