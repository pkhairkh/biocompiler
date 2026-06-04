"""
BioCompiler Validation Sub-package
====================================
Provides benchmarking modules that validate BioCompiler's heuristic
engines against experimental data.
"""

from .ground_truth import (
    GroundTruthEntry as GroundTruthEntry,
    ValidationResult as ValidationResult,
    GROUND_TRUTH_DATA as GROUND_TRUTH_DATA,
    validate_against_ground_truth as validate_against_ground_truth,
)

from .iedb_comparison import (
    IEDBBenchmarkEntry as IEDBBenchmarkEntry,
    benchmark_mhc_predictions as benchmark_mhc_predictions,
)

__all__ = [
    # ground_truth
    "GroundTruthEntry",
    "ValidationResult",
    "GROUND_TRUTH_DATA",
    "validate_against_ground_truth",
    # iedb_comparison
    "IEDBBenchmarkEntry",
    "benchmark_mhc_predictions",
]
