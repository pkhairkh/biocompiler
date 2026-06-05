"""
BioCompiler Validation Sub-package
====================================
Provides benchmarking modules that validate BioCompiler's heuristic
engines against experimental data.
"""

from .ground_truth import (
    GroundTruthEntry as GroundTruthEntry,
    ValidationResult as ValidationResult,
    GroundTruthResult as GroundTruthResult,
    GROUND_TRUTH_DATA as GROUND_TRUTH_DATA,
    validate_against_ground_truth as validate_against_ground_truth,
    validate_optimization_result as validate_optimization_result,
)

from .iedb_comparison import (
    IEDBBenchmarkEntry as IEDBBenchmarkEntry,
    IEDBComparisonResult as IEDBComparisonResult,
    benchmark_mhc_predictions as benchmark_mhc_predictions,
    compare_with_iedb as compare_with_iedb,
    get_available_alleles as get_available_alleles,
    get_known_binders as get_known_binders,
    get_known_non_binders as get_known_non_binders,
)

__all__ = [
    # ground_truth
    "GroundTruthEntry",
    "ValidationResult",
    "GroundTruthResult",
    "GROUND_TRUTH_DATA",
    "validate_against_ground_truth",
    "validate_optimization_result",
    # iedb_comparison
    "IEDBBenchmarkEntry",
    "IEDBComparisonResult",
    "benchmark_mhc_predictions",
    "compare_with_iedb",
    "get_available_alleles",
    "get_known_binders",
    "get_known_non_binders",
]
