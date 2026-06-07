"""Deprecated: use biocompiler.shared.types instead."""
import warnings

warnings.warn(
    "biocompiler.types is deprecated — use biocompiler.shared.types instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.types import *  # noqa: F401,F403

__all__ = [
    "SLOTMode",
    "Verdict",
    "five_valued_and",
    "five_valued_or",
    "three_valued_and",
    "three_valued_or",
    "combined_verdict",
    "PositionRange",
    "Token",
    "SpliceIsoform",
    "TypeCheckResult",
    "Certificate",
    "_VERDICT_ORDER",
]
