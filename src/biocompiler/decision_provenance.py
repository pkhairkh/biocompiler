"""Deprecated: use biocompiler.provenance.decision_provenance instead."""
import warnings

warnings.warn(
    "biocompiler.decision_provenance is deprecated — use biocompiler.provenance.decision_provenance instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.provenance.decision_provenance import *  # noqa: F401,F403

__all__ = [
    "CodonDecision",
    "ConstraintDecision",
    "OptimizationDecisionTrail",
    "DecisionProvenanceCollector",
    "ProvenanceStore",
]
