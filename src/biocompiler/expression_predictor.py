"""Deprecated: use biocompiler.expression.expression_predictor instead."""
import warnings

warnings.warn(
    "biocompiler.expression_predictor is deprecated — use biocompiler.expression.expression_predictor instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.expression.expression_predictor import *  # noqa: F401,F403

__all__ = [
    "ExpressionPrediction",
    "predict_expression",
    "ExpressionPredictor",
    "_compute_gc_optimality",
    "_FACTOR_WEIGHTS",
    "_GC_SWEET_SPOT",
    "_estimate_mrna_stability",
]
