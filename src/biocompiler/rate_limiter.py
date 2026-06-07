"""Deprecated: use biocompiler.shared.rate_limiter instead."""
import warnings

warnings.warn(
    "biocompiler.rate_limiter is deprecated — use biocompiler.shared.rate_limiter instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.rate_limiter import *  # noqa: F401,F403
