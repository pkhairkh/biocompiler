"""Deprecated: use biocompiler.sequence.pattern_enforcement instead."""
import warnings

warnings.warn(
    "biocompiler.pattern_enforcement is deprecated — use biocompiler.sequence.pattern_enforcement instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.sequence.pattern_enforcement import *  # noqa: F401,F403

__all__ = [
    "PatternConstraint",
    "PatternResult",
    "check_pattern",
    "check_patterns",
    "enforce_pattern",
    "enforce_patterns",
    "build_avoidance_scanner",
]
