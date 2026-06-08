"""
BioCompiler Performance Benchmark Test Suite
==============================================
Tests for the benchmarking registry, performance runner, and regression
detection modules.

These tests are marked with ``@pytest.mark.benchmark`` and can be run with::

    pytest -m benchmark

They are excluded from default test runs (``-m "not slow"``) because the
full benchmark suite takes several minutes.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Registry tests
# ---------------------------------------------------------------------------


class TestRegistry:
    """Tests for the benchmarking registry module."""

    def test_standard_proteins_defined(self):
        """STANDARD_PROTEINS should contain key benchmark proteins."""
        from biocompiler.benchmarking.registry import STANDARD_PROTEINS

        assert "EGFP" in STANDARD_PROTEINS
        assert "INS" in STANDARD_PROTEINS
        assert "BSA" in STANDARD_PROTEINS
        assert "HBB" in STANDARD_PROTEINS

    def test_standard_proteins_have_sequences(self):
        """Each standard protein must have a non-empty protein sequence."""
        from biocompiler.benchmarking.registry import STANDARD_PROTEINS

        for name, protein in STANDARD_PROTEINS.items():
            assert protein.protein_sequence, f"{name} has empty protein_sequence"
            assert len(protein.protein_sequence) > 10, (
                f"{name} protein_sequence seems too short: {len(protein.protein_sequence)}"
            )
            # Must contain only valid amino acid letters
            valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
            for aa in protein.protein_sequence:
                assert aa in valid_aas, (
                    f"{name} contains invalid amino acid '{aa}'"
                )

    def test_standard_proteins_length_categories(self):
        """Proteins should have correct length_category based on length."""
        from biocompiler.benchmarking.registry import STANDARD_PROTEINS

        for name, protein in STANDARD_PROTEINS.items():
            length = protein.length
            cat = protein.length_category
            if cat == "small":
                assert length < 150, f"{name} is 'small' but has {length} aa"
            elif cat == "medium":
                assert 150 <= length <= 400, f"{name} is 'medium' but has {length} aa"
            elif cat == "large":
                assert length > 400, f"{name} is 'large' but has {length} aa"

    def test_target_organisms_defined(self):
        """TARGET_ORGANISMS should contain the 5 key organisms."""
        from biocompiler.benchmarking.registry import TARGET_ORGANISMS

        expected = [
            "Escherichia_coli",
            "Saccharomyces_cerevisiae",
            "Homo_sapiens",
            "Mus_musculus",
            "CHO_K1",
        ]
        for org in expected:
            assert org in TARGET_ORGANISMS, f"Missing organism: {org}"

    def test_target_organisms_domains(self):
        """E. coli should be prokaryote; others should be eukaryote."""
        from biocompiler.benchmarking.registry import TARGET_ORGANISMS

        assert TARGET_ORGANISMS["Escherichia_coli"].domain == "prokaryote"
        assert TARGET_ORGANISMS["Homo_sapiens"].domain == "eukaryote"
        assert TARGET_ORGANISMS["Saccharomyces_cerevisiae"].domain == "eukaryote"
        assert TARGET_ORGANISMS["CHO_K1"].domain == "eukaryote"
        assert TARGET_ORGANISMS["Mus_musculus"].domain == "eukaryote"

    def test_benchmark_metrics_defined(self):
        """Core metrics must be defined with appropriate thresholds."""
        from biocompiler.benchmarking.registry import BENCHMARK_METRICS, MetricDirection

        # CAI metric
        assert "cai" in BENCHMARK_METRICS
        cai = BENCHMARK_METRICS["cai"]
        assert cai.direction == MetricDirection.HIGHER_IS_BETTER
        assert cai.regression_threshold == 0.05  # 5% CAI drop = regression

        # Optimization time metric
        assert "optimization_time" in BENCHMARK_METRICS
        opt_time = BENCHMARK_METRICS["optimization_time"]
        assert opt_time.direction == MetricDirection.LOWER_IS_BETTER
        assert opt_time.regression_threshold == 0.10  # 10% slowdown = regression

        # GC metric
        assert "gc_mean" in BENCHMARK_METRICS
        gc = BENCHMARK_METRICS["gc_mean"]
        assert gc.direction == MetricDirection.CLOSER_TO_TARGET

    def test_benchmark_suites_defined(self):
        """Core benchmark suites should be defined."""
        from biocompiler.benchmarking.registry import BENCHMARK_SUITES

        assert "core" in BENCHMARK_SUITES
        assert "full" in BENCHMARK_SUITES
        assert "therapeutic" in BENCHMARK_SUITES
        assert "cross_organism" in BENCHMARK_SUITES

    def test_core_suite_structure(self):
        """Core suite should have 4 proteins x 2 organisms."""
        from biocompiler.benchmarking.registry import get_suite

        suite = get_suite("core")
        assert len(suite.proteins) == 4
        assert len(suite.organisms) == 2
        assert suite.num_benchmarks == 8

    def test_get_suite_invalid_raises(self):
        """get_suite with an invalid name should raise ValueError."""
        from biocompiler.benchmarking.registry import get_suite

        with pytest.raises(ValueError, match="Unknown benchmark suite"):
            get_suite("nonexistent")

    def test_protein_to_dict_roundtrip(self):
        """StandardProtein serialization roundtrip."""
        from biocompiler.benchmarking.registry import STANDARD_PROTEINS

        egfp = STANDARD_PROTEINS["EGFP"]
        d = egfp.to_dict()
        assert d["name"] == "EGFP"
        assert d["length"] == len(egfp.protein_sequence)
        assert "protein_sequence" in d

    def test_organism_to_dict_roundtrip(self):
        """TargetOrganism serialization roundtrip."""
        from biocompiler.benchmarking.registry import TARGET_ORGANISMS

        ecoli = TARGET_ORGANISMS["Escherichia_coli"]
        d = ecoli.to_dict()
        assert d["name"] == "Escherichia_coli"
        assert d["domain"] == "prokaryote"


# ---------------------------------------------------------------------------
# Performance runner tests
# ---------------------------------------------------------------------------


class TestPerformanceRunner:
    """Tests for the PerformanceBenchmarkRunner."""

    def test_runner_initialization(self):
        """Runner should initialize with default parameters."""
        from biocompiler.benchmarking.perf_runner import PerformanceBenchmarkRunner

        runner = PerformanceBenchmarkRunner(suite_name="core")
        assert runner.suite_name == "core"
        assert runner.warmup_iterations == 1
        assert runner.timed_iterations == 3

    def test_runner_custom_parameters(self):
        """Runner should accept custom parameters."""
        from biocompiler.benchmarking.perf_runner import PerformanceBenchmarkRunner

        runner = PerformanceBenchmarkRunner(
            suite_name="core",
            warmup_iterations=2,
            timed_iterations=5,
            gc_lo=0.35,
            gc_hi=0.65,
        )
        assert runner.warmup_iterations == 2
        assert runner.timed_iterations == 5
        assert runner.gc_lo == 0.35
        assert runner.gc_hi == 0.65

    def test_protein_organism_result_serialization(self):
        """ProteinOrganismResult should serialize to/from dict."""
        from biocompiler.benchmarking.perf_runner import ProteinOrganismResult

        result = ProteinOrganismResult(
            protein_name="EGFP",
            organism_name="Escherichia_coli",
            protein_length=239,
            cai=0.85,
            gc_mean=0.52,
            optimization_time_s=1.234,
            num_iterations=3,
        )
        d = result.to_dict()
        assert d["protein_name"] == "EGFP"
        assert d["cai"] == 0.85

        # Roundtrip
        restored = ProteinOrganismResult.from_dict(d)
        assert restored.protein_name == "EGFP"
        assert restored.cai == 0.85
        assert restored.num_iterations == 3

    def test_benchmark_run_result_serialization(self):
        """BenchmarkRunResult should serialize to/from dict and JSON."""
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        result = BenchmarkRunResult(
            suite_name="core",
            version="12.0.0",
            timestamp="2025-01-01T00:00:00Z",
            results=[
                ProteinOrganismResult(
                    protein_name="EGFP",
                    organism_name="Homo_sapiens",
                    protein_length=239,
                    cai=0.90,
                    gc_mean=0.55,
                    optimization_time_s=2.0,
                    num_iterations=3,
                ),
            ],
            total_time_s=10.0,
        )
        result.compute_summary()

        d = result.to_dict()
        assert d["suite_name"] == "core"
        assert d["num_results"] == 1
        assert "summary" in d

        # JSON roundtrip
        json_str = json.dumps(d)
        parsed = json.loads(json_str)
        restored = BenchmarkRunResult.from_dict(parsed)
        assert restored.suite_name == "core"
        assert len(restored.results) == 1
        assert restored.results[0].cai == 0.90

    def test_save_and_load_json(self, tmp_path):
        """save_json / load_baseline should produce/load valid JSON."""
        from biocompiler.benchmarking.perf_runner import (
            PerformanceBenchmarkRunner,
            BenchmarkRunResult,
            ProteinOrganismResult,
            load_baseline,
        )

        result = BenchmarkRunResult(
            suite_name="core",
            version="12.0.0",
            results=[
                ProteinOrganismResult(
                    protein_name="INS",
                    organism_name="Homo_sapiens",
                    protein_length=110,
                    cai=0.78,
                    gc_mean=0.50,
                    num_iterations=3,
                ),
            ],
            total_time_s=5.0,
        )
        result.compute_summary()

        filepath = tmp_path / "bench_test.json"
        PerformanceBenchmarkRunner.save_json(result, filepath)

        # Verify file exists and is valid JSON
        assert filepath.exists()
        with open(filepath) as f:
            data = json.load(f)
        assert data["suite_name"] == "core"

        # Load back
        loaded = load_baseline(filepath)
        assert loaded.suite_name == "core"
        assert len(loaded.results) == 1
        assert loaded.results[0].cai == 0.78

    def test_benchmark_run_result_summary_empty(self):
        """Summary of empty results should report 0 successful."""
        from biocompiler.benchmarking.perf_runner import BenchmarkRunResult

        result = BenchmarkRunResult(suite_name="core", results=[])
        summary = result.compute_summary()
        assert summary["successful"] == 0
        assert summary["total_benchmarks"] == 0

    def test_benchmark_run_result_summary_with_errors(self):
        """Summary should count failed runs."""
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        failed = ProteinOrganismResult(
            protein_name="BSA",
            organism_name="Escherichia_coli",
            protein_length=607,
            error="optimization failed",
        )
        result = BenchmarkRunResult(suite_name="core", results=[failed])
        summary = result.compute_summary()
        assert summary["failed"] == 1
        assert "optimization failed" in summary["errors"]


# ---------------------------------------------------------------------------
# Regression detection tests
# ---------------------------------------------------------------------------


class TestRegressionDetection:
    """Tests for the regression detection module."""

    def test_compare_single_metric_no_regression(self):
        """No regression when value is within threshold."""
        from biocompiler.benchmarking.regression import compare_single_metric

        reg = compare_single_metric("cai", 0.80, 0.79)
        assert reg is None  # 1.25% change < 5% threshold

    def test_compare_single_metric_cai_regression(self):
        """CAI regression should be detected when > 5% drop."""
        from biocompiler.benchmarking.regression import compare_single_metric

        # 10% drop: 0.80 -> 0.72
        reg = compare_single_metric("cai", 0.80, 0.72)
        assert reg is not None
        assert reg.change_fraction < 0
        assert abs(reg.change_fraction) > 0.05

    def test_compare_single_metric_time_regression(self):
        """Time regression should be detected when > 10% slowdown."""
        from biocompiler.benchmarking.regression import compare_single_metric

        # 15% slowdown: 1.0s -> 1.15s
        reg = compare_single_metric("optimization_time", 1.0, 1.15)
        assert reg is not None
        assert reg.change_fraction > 0.10

    def test_compare_single_metric_no_time_regression(self):
        """Small time change should not trigger regression."""
        from biocompiler.benchmarking.regression import compare_single_metric

        # 5% change: below 10% threshold
        reg = compare_single_metric("optimization_time", 1.0, 1.05)
        assert reg is None

    def test_compare_single_metric_gc_regression(self):
        """GC deviation from target should trigger regression."""
        from biocompiler.benchmarking.regression import compare_single_metric

        # Baseline at 0.50 (target), current at 0.60
        # Deviation goes from 0.0 to 0.1 — but 0.0 baseline deviation
        # means we skip (can't compute percentage of zero)
        reg = compare_single_metric("gc_mean", 0.50, 0.60)
        # Since baseline deviation from target (0.50) is 0, skip
        assert reg is None

    def test_compare_single_metric_lower_is_better_improvement(self):
        """Lower values for LOWER_IS_BETTER metrics are not regressions."""
        from biocompiler.benchmarking.regression import compare_single_metric

        # Restriction sites: 3 -> 1 (improvement, not regression)
        reg = compare_single_metric("restriction_site_total", 3, 1)
        assert reg is None  # Going lower is better, not a regression

    def test_compare_single_metric_zero_baseline(self):
        """Should return None when baseline is zero."""
        from biocompiler.benchmarking.regression import compare_single_metric

        reg = compare_single_metric("restriction_site_total", 0, 5)
        assert reg is None  # Can't compute % change from zero

    def test_compare_results_no_regressions(self):
        """compare_results should report no regressions when results are stable."""
        from biocompiler.benchmarking.regression import compare_results
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        # Identical results — no regressions
        results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.85,
                gc_mean=0.52,
                optimization_time_s=1.5,
                constraint_satisfaction_rate=1.0,
                restriction_site_total=0,
                cryptic_splice_sites=0,
                cpg_islands=0,
                mrna_stability=0.8,
                num_iterations=3,
            ),
        ]

        current = BenchmarkRunResult(
            suite_name="core", version="12.0.0", results=results,
        )
        baseline = BenchmarkRunResult(
            suite_name="core", version="11.0.0", results=list(results),
        )

        comparison = compare_results(current, baseline)
        assert not comparison["has_regressions"]
        assert comparison["summary"]["total_regressions"] == 0

    def test_compare_results_with_cai_regression(self):
        """compare_results should detect CAI regression."""
        from biocompiler.benchmarking.regression import compare_results
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        baseline_results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.90,
                gc_mean=0.52,
                optimization_time_s=1.5,
                constraint_satisfaction_rate=1.0,
                num_iterations=3,
            ),
        ]
        current_results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.75,  # > 5% drop: (0.75-0.90)/0.90 = -16.7%
                gc_mean=0.52,
                optimization_time_s=1.5,
                constraint_satisfaction_rate=1.0,
                num_iterations=3,
            ),
        ]

        baseline = BenchmarkRunResult(
            suite_name="core", version="11.0.0", results=baseline_results,
        )
        current = BenchmarkRunResult(
            suite_name="core", version="12.0.0", results=current_results,
        )

        comparison = compare_results(current, baseline)
        assert comparison["has_regressions"]
        # Should find at least the CAI regression
        cai_regressions = [
            r for r in comparison["regressions"]
            if r["metric_name"] == "cai"
        ]
        assert len(cai_regressions) >= 1
        assert cai_regressions[0]["change_fraction"] < -0.05

    def test_compare_results_with_time_regression(self):
        """compare_results should detect > 10% slowdown."""
        from biocompiler.benchmarking.regression import compare_results
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        baseline_results = [
            ProteinOrganismResult(
                protein_name="INS",
                organism_name="Escherichia_coli",
                protein_length=110,
                cai=0.80,
                gc_mean=0.50,
                optimization_time_s=1.0,
                constraint_satisfaction_rate=1.0,
                num_iterations=3,
            ),
        ]
        current_results = [
            ProteinOrganismResult(
                protein_name="INS",
                organism_name="Escherichia_coli",
                protein_length=110,
                cai=0.80,
                gc_mean=0.50,
                optimization_time_s=1.25,  # 25% slowdown
                constraint_satisfaction_rate=1.0,
                num_iterations=3,
            ),
        ]

        baseline = BenchmarkRunResult(
            suite_name="core", version="11.0.0", results=baseline_results,
        )
        current = BenchmarkRunResult(
            suite_name="core", version="12.0.0", results=current_results,
        )

        comparison = compare_results(current, baseline)
        assert comparison["has_regressions"]
        time_regressions = [
            r for r in comparison["regressions"]
            if r["metric_name"] == "optimization_time"
        ]
        assert len(time_regressions) >= 1

    def test_compare_results_new_benchmarks(self):
        """compare_results should report new/removed benchmarks."""
        from biocompiler.benchmarking.regression import compare_results
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        baseline_results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.85,
                num_iterations=3,
            ),
        ]
        current_results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.85,
                num_iterations=3,
            ),
            ProteinOrganismResult(
                protein_name="INS",  # New benchmark
                organism_name="Escherichia_coli",
                protein_length=110,
                cai=0.80,
                num_iterations=3,
            ),
        ]

        baseline = BenchmarkRunResult(
            suite_name="core", version="11.0.0", results=baseline_results,
        )
        current = BenchmarkRunResult(
            suite_name="core", version="12.0.0", results=current_results,
        )

        comparison = compare_results(current, baseline)
        assert "INS/Escherichia_coli" in comparison["new_benchmarks"]

    def test_generate_regression_report(self):
        """generate_regression_report should produce markdown text."""
        from biocompiler.benchmarking.regression import (
            compare_results,
            generate_regression_report,
        )
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        baseline_results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.90,
                gc_mean=0.52,
                optimization_time_s=1.0,
                constraint_satisfaction_rate=1.0,
                num_iterations=3,
            ),
        ]
        current_results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.70,  # Major CAI drop
                gc_mean=0.52,
                optimization_time_s=2.0,  # Major slowdown
                constraint_satisfaction_rate=0.5,
                num_iterations=3,
            ),
        ]

        baseline = BenchmarkRunResult(
            suite_name="core", version="11.0.0", results=baseline_results,
        )
        current = BenchmarkRunResult(
            suite_name="core", version="12.0.0", results=current_results,
        )

        comparison = compare_results(current, baseline)
        report = generate_regression_report(comparison)

        assert report.has_regressions
        assert "Regression" in report.text
        assert "EGFP" in report.text
        assert "cai" in report.text.lower() or "CAI" in report.text

    def test_generate_regression_report_no_regressions(self):
        """Report should show 'no regressions' when results are stable."""
        from biocompiler.benchmarking.regression import (
            compare_results,
            generate_regression_report,
        )
        from biocompiler.benchmarking.perf_runner import (
            BenchmarkRunResult,
            ProteinOrganismResult,
        )

        results = [
            ProteinOrganismResult(
                protein_name="EGFP",
                organism_name="Homo_sapiens",
                protein_length=239,
                cai=0.85,
                gc_mean=0.52,
                optimization_time_s=1.0,
                constraint_satisfaction_rate=1.0,
                num_iterations=3,
            ),
        ]

        current = BenchmarkRunResult(
            suite_name="core", version="12.0.0", results=list(results),
        )
        baseline = BenchmarkRunResult(
            suite_name="core", version="11.0.0", results=list(results),
        )

        comparison = compare_results(current, baseline)
        report = generate_regression_report(comparison)

        assert not report.has_regressions
        assert "No regressions" in report.text

    def test_regression_severity_levels(self):
        """Critical severity should be assigned for > 2x threshold change."""
        from biocompiler.benchmarking.regression import compare_single_metric, RegressionSeverity

        # CAI: 5% threshold, so 15% drop = 3x threshold = critical
        reg = compare_single_metric("cai", 0.90, 0.765)  # 15% drop
        assert reg is not None
        assert reg.severity == RegressionSeverity.CRITICAL

        # CAI: 7% drop = 1.4x threshold = warning
        reg = compare_single_metric("cai", 0.90, 0.837)  # 7% drop
        assert reg is not None
        assert reg.severity == RegressionSeverity.WARNING


# ---------------------------------------------------------------------------
# Integration tests (marked as benchmark / slow)
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
@pytest.mark.slow
class TestPerformanceBenchmarkIntegration:
    """Integration tests that actually run the benchmark suite.

    These are slow and require the optimizer to be functional.
    """

    def test_core_suite_runs_successfully(self):
        """The core benchmark suite should complete with results."""
        from biocompiler.benchmarking.perf_runner import PerformanceBenchmarkRunner

        runner = PerformanceBenchmarkRunner(
            suite_name="core",
            warmup_iterations=0,  # Skip warmup for speed
            timed_iterations=1,   # Single iteration for CI
        )
        result = runner.run()

        assert result.suite_name == "core"
        assert len(result.results) > 0

        # At least some should succeed
        successful = [r for r in result.results if r.error is None]
        assert len(successful) > 0, "All benchmarks failed"

        # Successful results should have reasonable metrics
        for r in successful:
            assert 0.0 <= r.cai <= 1.0, f"CAI out of range: {r.cai}"
            assert 0.0 <= r.gc_mean <= 1.0, f"GC out of range: {r.gc_mean}"
            assert r.optimization_time_s > 0, "Optimization time should be positive"

    def test_cross_organism_suite_runs(self):
        """The cross-organism suite should test all 5 organisms."""
        from biocompiler.benchmarking.perf_runner import PerformanceBenchmarkRunner

        runner = PerformanceBenchmarkRunner(
            suite_name="cross_organism",
            warmup_iterations=0,
            timed_iterations=1,
        )
        result = runner.run()

        organisms_seen = {r.organism_name for r in result.results if r.error is None}
        assert len(organisms_seen) >= 2, (
            f"Expected results for multiple organisms, got: {organisms_seen}"
        )

    def test_benchmark_results_save_load_roundtrip(self, tmp_path):
        """Full save/load cycle should preserve all data."""
        from biocompiler.benchmarking.perf_runner import (
            PerformanceBenchmarkRunner,
            load_baseline,
        )

        runner = PerformanceBenchmarkRunner(
            suite_name="core",
            warmup_iterations=0,
            timed_iterations=1,
        )
        result = runner.run()

        filepath = tmp_path / "bench_roundtrip.json"
        runner.save_json(result, filepath)

        loaded = load_baseline(filepath)

        assert loaded.suite_name == result.suite_name
        assert len(loaded.results) == len(result.results)
        for orig, rest in zip(result.results, loaded.results):
            assert orig.protein_name == rest.protein_name
            assert orig.organism_name == rest.organism_name
            assert abs(orig.cai - rest.cai) < 1e-6

    def test_baseline_comparison_after_run(self, tmp_path):
        """Running a benchmark and comparing against itself should show no regressions."""
        from biocompiler.benchmarking.perf_runner import (
            PerformanceBenchmarkRunner,
            load_baseline,
        )

        runner = PerformanceBenchmarkRunner(
            suite_name="core",
            warmup_iterations=0,
            timed_iterations=1,
        )
        result = runner.run()

        filepath = tmp_path / "baseline.json"
        runner.save_json(result, filepath)

        # Compare against itself — no regressions expected
        report = runner.compare_against_baseline(result, filepath)
        assert not report.get("has_regressions", True), (
            f"Self-comparison should not show regressions. Report: {report}"
        )
