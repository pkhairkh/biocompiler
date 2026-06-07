"""Deprecated: use biocompiler.shared.exceptions instead."""
import warnings

warnings.warn(
    "biocompiler.exceptions is deprecated — use biocompiler.shared.exceptions instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.exceptions import *  # noqa: F401,F403
