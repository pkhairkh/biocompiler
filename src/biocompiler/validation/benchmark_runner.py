"""Automated benchmark runner for wet-lab validation.

Runs the biocompiler optimizer on published protein sequences and compares
the output CAI/expression predictions against measured experimental data.

Usage:
    python -m biocompiler.validation.benchmark_runner [--output report.json]

Or programmatically:
    from biocompiler.validation.benchmark_runner import run_benchmark_suite
    report = run_benchmark_suite()
"""
from __future__ import annotations

import json
import logging
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "BenchmarkResult",
    "BenchmarkReport",
    "run_single_benchmark",
    "run_benchmark_suite",
    "format_report_text",
]


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run.

    Attributes:
        gene_name: Name of the gene benchmarked.
        organism: Target organism.
        measured_expression: Measured expression level from the published data.
        predicted_cai: CAI of the sequence produced by biocompiler.
        expected_cai: CAI reported in the published dataset.
        cai_error: Difference between predicted and expected CAI.
        cai_relative_error: Relative CAI error (|cai_error| / max(expected, 0.01)).
        optimization_time_seconds: Wall-clock time for the optimization.
        passed: Whether CAI is within tolerance.
        notes: Additional notes about this benchmark.
    """

    gene_name: str
    organism: str
    measured_expression: float
    predicted_cai: float
    expected_cai: float
    cai_error: float  # predicted - expected
    cai_relative_error: float  # |cai_error| / max(expected, 0.01)
    optimization_time_seconds: float
    passed: bool  # whether CAI is within tolerance
    notes: str = ""


@dataclass
class BenchmarkReport:
    """Full benchmark report across all datasets.

    Attributes:
        timestamp: ISO-8601 timestamp of when the benchmark was run.
        biocompiler_version: Version string of biocompiler.
        results: Per-benchmark results.
        summary: Aggregate summary statistics.
    """

    timestamp: str
    biocompiler_version: str
    results: list[BenchmarkResult] = field(default_factory=list)
    summary: dict = field(default_factory=dict)


def run_single_benchmark(
    protein: str,
    organism: str,
    expected_cai: float,
    measured_expression: float = 0.0,
    cai_tolerance: float = 0.05,
    gene_name: str = "",
) -> BenchmarkResult:
    """Run a single benchmark: optimize the protein and compare CAI.

    Parameters
    ----------
    protein : str
        Amino acid sequence to optimize.
    organism : str
        Target organism (e.g. "Escherichia_coli", "Homo_sapiens").
    expected_cai : float
        CAI value reported in the published dataset.
    measured_expression : float
        Measured expression level from the publication (for reference).
    cai_tolerance : float
        Maximum acceptable absolute CAI deviation for the benchmark
        to be considered "passed" (default 0.05).
    gene_name : str
        Name of the gene (for reporting).

    Returns
    -------
    BenchmarkResult
        The benchmark result with predicted CAI and pass/fail status.
    """
    from ..optimizer.pipeline_core import optimize_sequence

    start_time = time.perf_counter()
    try:
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            strategy="hybrid",
        )
        predicted_cai = result.cai
    except Exception as e:
        logger.error("Benchmark failed for %s: %s", gene_name, e)
        predicted_cai = 0.0

    opt_time = time.perf_counter() - start_time
    cai_error = predicted_cai - expected_cai
    cai_rel_error = abs(cai_error) / max(expected_cai, 0.01)

    return BenchmarkResult(
        gene_name=gene_name,
        organism=organism,
        measured_expression=measured_expression,
        predicted_cai=round(predicted_cai, 4),
        expected_cai=round(expected_cai, 4),
        cai_error=round(cai_error, 4),
        cai_relative_error=round(cai_rel_error, 4),
        optimization_time_seconds=round(opt_time, 3),
        passed=abs(cai_error) <= cai_tolerance,
    )


def run_benchmark_suite(
    cai_tolerance: float = 0.05,
    output_path: Optional[str] = None,
) -> BenchmarkReport:
    """Run the full benchmark suite against published data.

    Parameters
    ----------
    cai_tolerance : float
        Maximum acceptable absolute CAI deviation (default 0.05).
    output_path : str, optional
        If provided, write the JSON report to this file path.

    Returns
    -------
    BenchmarkReport
        The complete benchmark report with per-gene results and summary.
    """
    from .published_expression_data import ALL_PUBLISHED_DATASETS

    report = BenchmarkReport(
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        biocompiler_version=_get_version(),
    )

    for dataset_name, entries in ALL_PUBLISHED_DATASETS.items():
        logger.info("Running dataset: %s (%d entries)", dataset_name, len(entries))
        for entry in entries:
            if not entry.protein_sequence:
                continue
            result = run_single_benchmark(
                protein=entry.protein_sequence,
                organism=entry.organism,
                expected_cai=entry.cai_predicted,
                measured_expression=entry.measured_expression_level,
                cai_tolerance=cai_tolerance,
                gene_name=entry.gene_name,
            )
            report.results.append(result)

    # Compute summary
    total = len(report.results)
    passed = sum(1 for r in report.results if r.passed)
    report.summary = {
        "total_benchmarks": total,
        "passed": passed,
        "failed": total - passed,
        "pass_rate": round(passed / max(total, 1), 4),
        "mean_cai_error": round(
            sum(r.cai_error for r in report.results) / max(total, 1), 4
        ),
        "max_cai_error": round(
            max((abs(r.cai_error) for r in report.results), default=0.0), 4
        ),
        "cai_tolerance": cai_tolerance,
    }

    if output_path:
        with open(output_path, "w") as f:
            json.dump(asdict(report), f, indent=2, default=str)
        logger.info("Benchmark report written to: %s", output_path)

    return report


def format_report_text(report: BenchmarkReport) -> str:
    """Format a BenchmarkReport as a human-readable text string.

    Parameters
    ----------
    report : BenchmarkReport
        The report to format.

    Returns
    -------
    str
        Human-readable report text.
    """
    lines: list[str] = []
    lines.append("=" * 70)
    lines.append("BIOCOMPILER WET-LAB VALIDATION BENCHMARK REPORT")
    lines.append("=" * 70)
    lines.append(f"Timestamp:    {report.timestamp}")
    lines.append(f"Version:      {report.biocompiler_version}")
    lines.append("")

    # Per-result table
    lines.append(f"{'Gene':<12} {'Organism':<22} {'Pred CAI':>9} {'Exp CAI':>9} "
                 f"{'Error':>8} {'Pass':>5} {'Time(s)':>8}")
    lines.append("-" * 70)
    for r in report.results:
        lines.append(
            f"{r.gene_name:<12} {r.organism:<22} {r.predicted_cai:>9.4f} "
            f"{r.expected_cai:>9.4f} {r.cai_error:>+8.4f} "
            f"{'YES' if r.passed else 'NO':>5} {r.optimization_time_seconds:>8.3f}"
        )
    lines.append("")

    # Summary
    s = report.summary
    lines.append("=" * 70)
    lines.append("SUMMARY")
    lines.append("=" * 70)
    lines.append(f"Total benchmarks: {s.get('total_benchmarks', 0)}")
    lines.append(f"Passed:           {s.get('passed', 0)}")
    lines.append(f"Failed:           {s.get('failed', 0)}")
    lines.append(f"Pass rate:        {s.get('pass_rate', 0.0):.1%}")
    lines.append(f"Mean CAI error:   {s.get('mean_cai_error', 0.0):+.4f}")
    lines.append(f"Max |CAI error|:  {s.get('max_cai_error', 0.0):.4f}")
    lines.append(f"CAI tolerance:    {s.get('cai_tolerance', 0.05):.4f}")
    lines.append("=" * 70)

    return "\n".join(lines)


def _get_version() -> str:
    """Get biocompiler version string."""
    try:
        from .. import __version__
        return __version__
    except Exception:
        return "unknown"


# ────────────────────────────────────────────────────────────
# CLI entry point
# ────────────────────────────────────────────────────────────
def main() -> None:
    """CLI entry point for running the benchmark suite."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run biocompiler wet-lab validation benchmarks"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output path for JSON report (default: stdout only)",
    )
    parser.add_argument(
        "--tolerance", "-t",
        type=float,
        default=0.05,
        help="CAI tolerance for pass/fail (default: 0.05)",
    )
    parser.add_argument(
        "--json-only",
        action="store_true",
        help="Only output JSON (no text report)",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    report = run_benchmark_suite(
        cai_tolerance=args.tolerance,
        output_path=args.output,
    )

    if not args.json_only:
        print(format_report_text(report))

    # Exit with non-zero code if any benchmarks failed
    if report.summary.get("failed", 0) > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
