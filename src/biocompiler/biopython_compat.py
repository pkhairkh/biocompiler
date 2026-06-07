"""Deprecated: use biocompiler.shared.biopython_compat instead."""
import warnings

warnings.warn(
    "biocompiler.biopython_compat is deprecated — use biocompiler.shared.biopython_compat instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.biopython_compat import *  # noqa: F401,F403
