"""Deprecated: use biocompiler.expression.utr_models instead."""
import warnings

warnings.warn(
    "biocompiler.utr_models is deprecated — use biocompiler.expression.utr_models instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.utr_models import *  # noqa: F401,F403

__all__ = [
    "UTRConfig",
    "ORGANISM_UTR_CONFIGS",
    "AVAILABLE_ORGANISMS",
    "score_5utr",
    "score_3utr",
    "suggest_5utr",
    "suggest_3utr",
]
