"""Deprecated: use biocompiler.provenance instead."""
import warnings

warnings.warn(
    "biocompiler.provenance is deprecated — use biocompiler.provenance instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance import *  # noqa: F401,F403

__all__ = [
    "DecisionRecord",
    "ProvenanceTracker",
    "OptimizationProvenance",
    "OptimizationRecord",
    "ProvenanceStore",
    "generate_provenance_report",
    "DECISION_CATEGORY_CAI",
    "DECISION_CATEGORY_GT_AVOIDANCE",
    "DECISION_CATEGORY_GC_CONTENT",
    "DECISION_CATEGORY_RESTRICTION_SITE",
    "DECISION_CATEGORY_SPLICE_PREVENTION",
    "DECISION_CATEGORY_MUTATION",
    "DECISION_CATEGORY_CONSTRAINT_RELAXATION",
    "DECISION_CATEGORY_OTHER",
    "ALL_DECISION_CATEGORIES",
]
