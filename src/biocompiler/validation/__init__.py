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

from .published_expression_data import (
    PublishedExpressionResult as PublishedExpressionResult,
    KUDLA_2009_GFP_DATASET as KUDLA_2009_GFP_DATASET,
    WELCH_2009_DATASET as WELCH_2009_DATASET,
    PUIGBO_2008_DATASET as PUIGBO_2008_DATASET,
    SHARP_LI_1987_DATASET as SHARP_LI_1987_DATASET,
    ALL_PUBLISHED_DATASETS as ALL_PUBLISHED_DATASETS,
)

from .benchmark_runner import (
    BenchmarkResult as BenchmarkResult,
    BenchmarkReport as BenchmarkReport,
    run_single_benchmark as run_single_benchmark,
    run_benchmark_suite as run_benchmark_suite,
    format_report_text as format_report_text,
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
    # published_expression_data
    "PublishedExpressionResult",
    "KUDLA_2009_GFP_DATASET",
    "WELCH_2009_DATASET",
    "PUIGBO_2008_DATASET",
    "SHARP_LI_1987_DATASET",
    "ALL_PUBLISHED_DATASETS",
    # benchmark_runner
    "BenchmarkResult",
    "BenchmarkReport",
    "run_single_benchmark",
    "run_benchmark_suite",
    "format_report_text",
]
