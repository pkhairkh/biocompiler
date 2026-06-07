"""Deprecated: use biocompiler.optimizer.hybrid_optimizer instead."""
import warnings

warnings.warn(
    "biocompiler.hybrid_optimizer is deprecated — use biocompiler.optimizer.hybrid_optimizer instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.optimizer.hybrid_optimizer import *  # noqa: F401,F403
