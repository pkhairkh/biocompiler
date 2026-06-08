"""
BioCompiler Performance Benchmark Runner
==========================================
Runs comprehensive optimization benchmarks across the standard protein x
organism matrix, with warmup iterations, JSON result output, and baseline
comparison.

This runner focuses on **BioCompiler's own performance** across releases,
unlike the head-to-head ``runner.py`` which compares against DNAchisel.

Usage::

    from biocompiler.benchmarking.perf_runner import PerformanceBenchmarkRunner

    runner = PerformanceBenchmarkRunner(suite_name="core")
    results = runner.run()
    runner.save_json(results, "benchmark_results.json")

    # Compare against a baseline from a previous release
    regression_report = runner.compare_against_baseline(
        results, "baseline_v11.json"
    )
"""

from __future__ import annotations

import json
import logging
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .registry import (
    BENCHMARK_SUITES,
    BenchmarkSuite,
    StandardProtein,
    TargetOrganism,
    get_suite,
)

logger = logging.getLogger(__name__)

__all__ = [
    "PerformanceBenchmarkRunner",
    "BenchmarkRunResult",
    "ProteinOrganismResult",
    "load_baseline",
]


# ---------------------------------------------------------------------------
# ProteinOrganismResult — result for a single (protein, organism) pair
# ---------------------------------------------------------------------------

class ProteinOrganismResult:
    """Result of benchmarking a single protein on a single organism.

    Attributes
    ----------
    protein_name : str
        Gene symbol.
    organism_name : str
        Target organism (normalised).
    protein_length : int
        Protein length in amino acids.
    cai : float
        Codon Adaptation Index of the optimized sequence.
    gc_mean : float
        Mean GC content fraction.
    gc_std : float
        Standard deviation of GC content.
    gc_min : float
        Minimum GC in any sliding window.
    gc_max : float
        Maximum GC in any sliding window.
    restriction_site_total : int
        Total restriction enzyme sites found.
    cryptic_splice_sites : int
        Number of cryptic splice-site dinucleotides.
    cpg_islands : int
        Number of CpG islands.
    mrna_stability : float
        Composite mRNA stability score.
    constraint_satisfaction_rate : float
        Fraction of constraints satisfied (0.0 - 1.0).
    optimization_time_s : float
        Median wall-clock optimization time in seconds.
    optimization_time_std : float
        Standard deviation of optimization times.
    num_iterations : int
        Number of timed iterations (excluding warmup).
    error : str | None
        Error message if optimization failed, else ``None``.
    """

    def __init__(
        self,
        protein_name: str,
        organism_name: str,
        protein_length: int,
        *,
        cai: float = 0.0,
        gc_mean: float = 0.0,
        gc_std: float = 0.0,
        gc_min: float = 0.0,
        gc_max: float = 0.0,
        restriction_site_total: int = 0,
        cryptic_splice_sites: int = 0,
        cpg_islands: int = 0,
        mrna_stability: float = 0.0,
        constraint_satisfaction_rate: float = 0.0,
        optimization_time_s: float = 0.0,
        optimization_time_std: float = 0.0,
        num_iterations: int = 0,
        error: str | None = None,
    ) -> None:
        self.protein_name = protein_name
        self.organism_name = organism_name
        self.protein_length = protein_length
        self.cai = cai
        self.gc_mean = gc_mean
        self.gc_std = gc_std
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.restriction_site_total = restriction_site_total
        self.cryptic_splice_sites = cryptic_splice_sites
        self.cpg_islands = cpg_islands
        self.mrna_stability = mrna_stability
        self.constraint_satisfaction_rate = constraint_satisfaction_rate
        self.optimization_time_s = optimization_time_s
        self.optimization_time_std = optimization_time_std
        self.num_iterations = num_iterations
        self.error = error

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "protein_name": self.protein_name,
            "organism_name": self.organism_name,
            "protein_length": self.protein_length,
            "cai": round(self.cai, 6),
            "gc_mean": round(self.gc_mean, 6),
            "gc_std": round(self.gc_std, 6),
            "gc_min": round(self.gc_min, 6),
            "gc_max": round(self.gc_max, 6),
            "restriction_site_total": self.restriction_site_total,
            "cryptic_splice_sites": self.cryptic_splice_sites,
            "cpg_islands": self.cpg_islands,
            "mrna_stability": round(self.mrna_stability, 6),
            "constraint_satisfaction_rate": round(self.constraint_satisfaction_rate, 6),
            "optimization_time_s": round(self.optimization_time_s, 6),
            "optimization_time_std": round(self.optimization_time_std, 6),
            "num_iterations": self.num_iterations,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ProteinOrganismResult:
        """Deserialize from a dict (e.g. loaded from JSON)."""
        return cls(
            protein_name=data["protein_name"],
            organism_name=data["organism_name"],
            protein_length=data["protein_length"],
            cai=data.get("cai", 0.0),
            gc_mean=data.get("gc_mean", 0.0),
            gc_std=data.get("gc_std", 0.0),
            gc_min=data.get("gc_min", 0.0),
            gc_max=data.get("gc_max", 0.0),
            restriction_site_total=data.get("restriction_site_total", 0),
            cryptic_splice_sites=data.get("cryptic_splice_sites", 0),
            cpg_islands=data.get("cpg_islands", 0),
            mrna_stability=data.get("mrna_stability", 0.0),
            constraint_satisfaction_rate=data.get("constraint_satisfaction_rate", 0.0),
            optimization_time_s=data.get("optimization_time_s", 0.0),
            optimization_time_std=data.get("optimization_time_std", 0.0),
            num_iterations=data.get("num_iterations", 0),
            error=data.get("error"),
        )


# ---------------------------------------------------------------------------
# BenchmarkRunResult — aggregate result for a full benchmark run
# ---------------------------------------------------------------------------

class BenchmarkRunResult:
    """Aggregate result for a complete benchmark run.

    Attributes
    ----------
    suite_name : str
        Name of the benchmark suite that was run.
    version : str
        BioCompiler version at time of run.
    timestamp : str
        ISO 8601 timestamp of when the run started.
    results : list[ProteinOrganismResult]
        Individual (protein, organism) results.
    total_time_s : float
        Total wall-clock time for the entire run.
    summary : dict
        Aggregate summary statistics.
    """

    def __init__(
        self,
        suite_name: str,
        version: str = "",
        timestamp: str = "",
        results: list[ProteinOrganismResult] | None = None,
        total_time_s: float = 0.0,
        summary: dict[str, Any] | None = None,
    ) -> None:
        self.suite_name = suite_name
        self.version = version
        self.timestamp = timestamp
        self.results = results or []
        self.total_time_s = total_time_s
        self.summary = summary or {}

    def compute_summary(self) -> dict[str, Any]:
        """Compute and store aggregate summary statistics.

        Returns
        -------
        dict
            Summary with means, medians, and success rates.
        """
        successful = [r for r in self.results if r.error is None]
        failed = [r for r in self.results if r.error is not None]

        if not successful:
            self.summary = {
                "total_benchmarks": len(self.results),
                "successful": 0,
                "failed": len(failed),
                "errors": [r.error for r in failed],
            }
            return self.summary

        n = len(successful)
        summary: dict[str, Any] = {
            "total_benchmarks": len(self.results),
            "successful": n,
            "failed": len(failed),
            "errors": [r.error for r in failed],
            # CAI stats
            "mean_cai": round(statistics.mean(r.cai for r in successful), 4),
            "median_cai": round(statistics.median(r.cai for r in successful), 4),
            "min_cai": round(min(r.cai for r in successful), 4),
            "max_cai": round(max(r.cai for r in successful), 4),
            # GC stats
            "mean_gc": round(statistics.mean(r.gc_mean for r in successful), 4),
            "median_gc": round(statistics.median(r.gc_mean for r in successful), 4),
            # Time stats
            "mean_optimization_time_s": round(
                statistics.mean(r.optimization_time_s for r in successful), 4
            ),
            "median_optimization_time_s": round(
                statistics.median(r.optimization_time_s for r in successful), 4
            ),
            "max_optimization_time_s": round(
                max(r.optimization_time_s for r in successful), 4
            ),
            # Constraint satisfaction
            "mean_constraint_satisfaction_rate": round(
                statistics.mean(r.constraint_satisfaction_rate for r in successful), 4
            ),
            # Per-organism breakdown
            "per_organism": {},
        }

        # Per-organism breakdown
        by_organism: dict[str, list[ProteinOrganismResult]] = {}
        for r in successful:
            by_organism.setdefault(r.organism_name, []).append(r)

        for org_name, org_results in sorted(by_organism.items()):
            summary["per_organism"][org_name] = {
                "count": len(org_results),
                "mean_cai": round(statistics.mean(r.cai for r in org_results), 4),
                "mean_gc": round(statistics.mean(r.gc_mean for r in org_results), 4),
                "mean_optimization_time_s": round(
                    statistics.mean(r.optimization_time_s for r in org_results), 4
                ),
                "mean_constraint_satisfaction_rate": round(
                    statistics.mean(r.constraint_satisfaction_rate for r in org_results), 4
                ),
            }

        self.summary = summary
        return summary

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "suite_name": self.suite_name,
            "version": self.version,
            "timestamp": self.timestamp,
            "total_time_s": round(self.total_time_s, 4),
            "num_results": len(self.results),
            "summary": self.summary,
            "results": [r.to_dict() for r in self.results],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> BenchmarkRunResult:
        """Deserialize from a dict."""
        results = [
            ProteinOrganismResult.from_dict(r)
            for r in data.get("results", [])
        ]
        return cls(
            suite_name=data.get("suite_name", ""),
            version=data.get("version", ""),
            timestamp=data.get("timestamp", ""),
            results=results,
            total_time_s=data.get("total_time_s", 0.0),
            summary=data.get("summary", {}),
        )


# ---------------------------------------------------------------------------
# PerformanceBenchmarkRunner
# ---------------------------------------------------------------------------

class PerformanceBenchmarkRunner:
    """Runs performance benchmarks across the protein x organism matrix.

    Features:
      - Configurable warmup iterations to prime JIT / caches
      - Multiple timed iterations with median / std reporting
      - JSON result serialization
      - Baseline comparison for regression detection

    Parameters
    ----------
    suite_name : str
        Name of the benchmark suite from the registry.
    warmup_iterations : int
        Number of warmup iterations (not timed). Default 1.
    timed_iterations : int
        Number of timed iterations. Median time is reported. Default 3.
    gc_lo : float
        Minimum GC content fraction for optimization. Default 0.30.
    gc_hi : float
        Maximum GC content fraction for optimization. Default 0.70.

    Examples
    --------
    >>> runner = PerformanceBenchmarkRunner(suite_name="core")
    >>> results = runner.run()
    >>> runner.save_json(results, "bench_v12.json")
    """

    def __init__(
        self,
        suite_name: str = "core",
        warmup_iterations: int = 1,
        timed_iterations: int = 3,
        gc_lo: float = 0.30,
        gc_hi: float = 0.70,
    ) -> None:
        self.suite_name = suite_name
        self.suite = get_suite(suite_name)
        self.warmup_iterations = warmup_iterations
        self.timed_iterations = timed_iterations
        self.gc_lo = gc_lo
        self.gc_hi = gc_hi

    def run(self) -> BenchmarkRunResult:
        """Run all benchmarks in the suite.

        Returns
        -------
        BenchmarkRunResult
            Complete results including per-(protein, organism) metrics
            and aggregate summary.
        """
        version = self._get_version()
        timestamp = datetime.now(timezone.utc).isoformat()

        logger.info(
            "Starting benchmark suite '%s' (%d proteins x %d organisms = %d combos)",
            self.suite_name,
            len(self.suite.proteins),
            len(self.suite.organisms),
            self.suite.num_benchmarks,
        )

        t_start = time.perf_counter()
        results: list[ProteinOrganismResult] = []

        for protein_name, protein in self.suite.proteins.items():
            for organism in self.suite.organisms:
                result = self._run_single(protein, organism)
                results.append(result)

        total_time = time.perf_counter() - t_start

        run_result = BenchmarkRunResult(
            suite_name=self.suite_name,
            version=version,
            timestamp=timestamp,
            results=results,
            total_time_s=total_time,
        )
        run_result.compute_summary()

        logger.info(
            "Benchmark suite '%s' completed: %d/%d successful in %.1fs",
            self.suite_name,
            run_result.summary.get("successful", 0),
            len(results),
            total_time,
        )

        return run_result

    def _run_single(
        self,
        protein: StandardProtein,
        organism: TargetOrganism,
    ) -> ProteinOrganismResult:
        """Run benchmark for a single (protein, organism) pair.

        Performs warmup iterations (untimed) followed by timed iterations.
        Returns the median metrics across timed iterations.
        """
        logger.info(
            "Benchmarking %s (%d aa) on %s",
            protein.name, protein.length, organism.display_name,
        )

        # Warmup iterations (not timed)
        for _ in range(self.warmup_iterations):
            try:
                self._optimize_once(protein, organism)
            except Exception as exc:
                logger.warning(
                    "Warmup failed for %s/%s: %s",
                    protein.name, organism.name, exc,
                )

        # Timed iterations
        iteration_data: list[dict[str, Any]] = []
        times: list[float] = []

        for _ in range(self.timed_iterations):
            try:
                t0 = time.perf_counter()
                metrics = self._optimize_once(protein, organism)
                elapsed = time.perf_counter() - t0
                times.append(elapsed)
                metrics["optimization_time_s"] = elapsed
                iteration_data.append(metrics)
            except Exception as exc:
                logger.error(
                    "Benchmark failed for %s/%s: %s",
                    protein.name, organism.name, exc,
                )
                return ProteinOrganismResult(
                    protein_name=protein.name,
                    organism_name=organism.name,
                    protein_length=protein.length,
                    num_iterations=0,
                    error=str(exc),
                )

        if not iteration_data:
            return ProteinOrganismResult(
                protein_name=protein.name,
                organism_name=organism.name,
                protein_length=protein.length,
                num_iterations=0,
                error="All iterations failed",
            )

        # Use the last iteration's quality metrics (optimization is deterministic
        # for a given input, so quality doesn't vary across iterations).
        # Timing uses median across iterations.
        last_metrics = iteration_data[-1]
        median_time = statistics.median(times) if times else 0.0
        std_time = statistics.stdev(times) if len(times) > 1 else 0.0

        return ProteinOrganismResult(
            protein_name=protein.name,
            organism_name=organism.name,
            protein_length=protein.length,
            cai=last_metrics.get("cai", 0.0),
            gc_mean=last_metrics.get("gc_mean", 0.0),
            gc_std=last_metrics.get("gc_std", 0.0),
            gc_min=last_metrics.get("gc_min", 0.0),
            gc_max=last_metrics.get("gc_max", 0.0),
            restriction_site_total=last_metrics.get("restriction_site_total", 0),
            cryptic_splice_sites=last_metrics.get("cryptic_splice_sites", 0),
            cpg_islands=last_metrics.get("cpg_islands", 0),
            mrna_stability=last_metrics.get("mrna_stability", 0.0),
            constraint_satisfaction_rate=last_metrics.get(
                "constraint_satisfaction_rate", 0.0
            ),
            optimization_time_s=median_time,
            optimization_time_std=std_time,
            num_iterations=len(times),
        )

    def _optimize_once(
        self,
        protein: StandardProtein,
        organism: TargetOrganism,
    ) -> dict[str, Any]:
        """Run a single optimization and compute all metrics.

        Returns
        -------
        dict
            All computed metrics for this optimization run.
        """
        from ..optimizer import optimize_sequence

        result = optimize_sequence(
            target_protein=protein.protein_sequence,
            organism=organism.name,
            gc_lo=self.gc_lo,
            gc_hi=self.gc_hi,
        )

        dna = result.sequence

        # Compute comprehensive metrics
        from .metrics import compute_all_metrics, compute_cai_validated
        from .metrics import (
            compute_gc_distribution,
            count_restriction_sites,
            count_cryptic_splice_sites,
            count_cpg_islands,
            compute_mrna_stability_score,
        )
        from .metrics import STANDARD_ENZYME_PANEL

        # Use compute_all_metrics for comprehensive data
        bench_metrics = compute_all_metrics(
            dna=dna,
            protein=protein.protein_sequence,
            organism=organism.name,
            enzymes=list(STANDARD_ENZYME_PANEL),
        )

        # Compute constraint satisfaction rate
        constraints = {
            "gc_in_range": self.gc_lo <= bench_metrics.gc_profile.mean <= self.gc_hi,
            "no_restriction_sites": bench_metrics.restriction_site_total == 0,
            "cai_above_threshold": bench_metrics.cai >= 0.5,
            "low_cryptic_splice_sites": bench_metrics.cryptic_splice_sites == 0,
            "no_cpg_islands": bench_metrics.cpg_islands == 0,
            "good_mrna_stability": bench_metrics.mrna_stability >= 0.5,
        }
        constraint_rate = sum(constraints.values()) / len(constraints)

        return {
            "cai": bench_metrics.cai,
            "gc_mean": bench_metrics.gc_profile.mean,
            "gc_std": bench_metrics.gc_profile.std,
            "gc_min": bench_metrics.gc_profile.min_,
            "gc_max": bench_metrics.gc_profile.max_,
            "restriction_site_total": bench_metrics.restriction_site_total,
            "cryptic_splice_sites": bench_metrics.cryptic_splice_sites,
            "cpg_islands": bench_metrics.cpg_islands,
            "mrna_stability": bench_metrics.mrna_stability,
            "constraint_satisfaction_rate": constraint_rate,
        }

    @staticmethod
    def _get_version() -> str:
        """Get the current BioCompiler version."""
        try:
            from .. import __version__
            return __version__
        except (ImportError, AttributeError):
            pass
        try:
            from importlib.metadata import version
            return version("biocompiler")
        except Exception:
            return "unknown"

    # ── Serialization ──────────────────────────────────────────────────

    @staticmethod
    def save_json(
        result: BenchmarkRunResult,
        filepath: str | Path,
    ) -> Path:
        """Save benchmark results to a JSON file.

        Parameters
        ----------
        result : BenchmarkRunResult
            Benchmark results to save.
        filepath : str | Path
            Output file path.

        Returns
        -------
        Path
            The path the file was written to.
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(result.to_dict(), f, indent=2, ensure_ascii=False)
        logger.info("Benchmark results saved to %s", filepath)
        return filepath

    # ── Baseline comparison ────────────────────────────────────────────

    def compare_against_baseline(
        self,
        current: BenchmarkRunResult,
        baseline_path: str | Path,
    ) -> dict[str, Any]:
        """Compare current results against a baseline from a previous release.

        Parameters
        ----------
        current : BenchmarkRunResult
            Current benchmark results.
        baseline_path : str | Path
            Path to the baseline JSON file.

        Returns
        -------
        dict
            Comparison report with regressions highlighted.
        """
        from .regression import compare_results, generate_regression_report

        baseline = load_baseline(baseline_path)
        comparison = compare_results(current, baseline)
        report = generate_regression_report(comparison)
        return report


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def load_baseline(filepath: str | Path) -> BenchmarkRunResult:
    """Load a baseline benchmark result from a JSON file.

    Parameters
    ----------
    filepath : str | Path
        Path to the JSON baseline file.

    Returns
    -------
    BenchmarkRunResult
        Deserialized baseline results.
    """
    filepath = Path(filepath)
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
    return BenchmarkRunResult.from_dict(data)
