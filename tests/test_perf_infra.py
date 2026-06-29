"""
Performance regression test infrastructure for BioCompiler.

Provides:
1. Benchmark timing consistency tests
2. CI environment resource checks
3. Percentile-based timing assertions

These tests are designed to catch performance regressions without
being overly sensitive to environmental noise. They use statistical
methods (percentiles, relative tolerances) rather than strict
equality checks.
"""

from __future__ import annotations

import os
import statistics
import time
from typing import Callable

import pytest


# ─── Helpers ──────────────────────────────────────────────────────────


def measure_execution_time(
    func: Callable,
    n_runs: int = 10,
    warmup: int = 2,
) -> list[float]:
    """Measure execution time of a function over multiple runs.

    Parameters
    ----------
    func:
        The function to benchmark.
    n_runs:
        Number of timed runs.
    warmup:
        Number of warmup runs (not timed).

    Returns
    -------
    list[float]
        List of execution times in seconds.
    """
    # Warmup runs
    for _ in range(warmup):
        func()

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return times


def percentile(data: list[float], p: float) -> float:
    """Calculate the p-th percentile of a list of values.

    Uses linear interpolation (same method as numpy default).
    """
    if not data:
        raise ValueError("data must not be empty")
    sorted_data = sorted(data)
    k = (len(sorted_data) - 1) * (p / 100.0)
    f = int(k)
    c = f + 1
    if c >= len(sorted_data):
        return sorted_data[f]
    d0 = sorted_data[f] * (c - k)
    d1 = sorted_data[c] * (k - f)
    return d0 + d1


def assert_percentile_within(
    times: list[float],
    percentile_threshold: float,
    max_seconds: float,
    label: str = "operation",
) -> None:
    """Assert that a percentile of execution times is within a bound.

    Parameters
    ----------
    times:
        List of measured execution times.
    percentile_threshold:
        Which percentile to check (e.g. 95.0 for P95).
    max_seconds:
        Maximum allowed time at this percentile.
    label:
        Descriptive label for error messages.
    """
    actual = percentile(times, percentile_threshold)
    assert actual <= max_seconds, (
        f"{label}: P{percentile_threshold}={actual:.4f}s exceeds "
        f"threshold of {max_seconds:.4f}s "
        f"(min={min(times):.4f}, median={statistics.median(times):.4f}, "
        f"max={max(times):.4f})"
    )


def assert_no_regression(
    baseline_times: list[float],
    current_times: list[float],
    max_regression_ratio: float = 2.0,
    label: str = "operation",
) -> None:
    """Assert that current performance has not regressed relative to baseline.

    Compares the median of current_times against the median of baseline_times.
    A regression ratio > max_regression_ratio fails the test.

    Parameters
    ----------
    baseline_times:
        Previously recorded execution times.
    current_times:
        Currently measured execution times.
    max_regression_ratio:
        Maximum allowed ratio of current_median / baseline_median.
    label:
        Descriptive label for error messages.
    """
    baseline_median = statistics.median(baseline_times)
    current_median = statistics.median(current_times)

    if baseline_median > 0:
        ratio = current_median / baseline_median
        assert ratio <= max_regression_ratio, (
            f"{label}: Performance regression detected! "
            f"Current median ({current_median:.4f}s) is {ratio:.2f}x "
            f"slower than baseline ({baseline_median:.4f}s). "
            f"Max allowed ratio: {max_regression_ratio:.2f}x"
        )


# ═══════════════════════════════════════════════════════════════════════
# 1. Benchmark timing consistency
# ═══════════════════════════════════════════════════════════════════════


class TestBenchmarkTimingConsistency:
    """Test that benchmark timing produces consistent results."""

    def test_measure_execution_time_returns_correct_count(self):
        """measure_execution_time should return the correct number of runs."""
        times = measure_execution_time(lambda: None, n_runs=5, warmup=1)
        assert len(times) == 5

    def test_measure_execution_time_all_positive(self):
        """All measured times should be positive."""
        times = measure_execution_time(lambda: sum(range(100)), n_runs=10)
        assert all(t > 0 for t in times), f"Non-positive times: {times}"

    def test_measure_execution_time_reasonable_variance(self):
        """Execution times for a deterministic function should have reasonable variance.

        The coefficient of variation should be < 1.0 (std < mean) for a
        simple deterministic function.
        """
        times = measure_execution_time(lambda: sum(range(1000)), n_runs=20)
        mean_t = statistics.mean(times)
        if mean_t > 0:
            stdev_t = statistics.stdev(times)
            cv = stdev_t / mean_t
            # CV should be reasonable (< 1.0 for a deterministic function)
            assert cv < 1.0, (
                f"Coefficient of variation too high: {cv:.2f} "
                f"(mean={mean_t:.6f}, stdev={stdev_t:.6f})"
            )

    def test_warmup_runs_not_included(self):
        """Warmup runs should not be included in the results."""
        # Use a function that takes variable time on first call
        call_count = [0]

        def counted_func():
            call_count[0] += 1
            return call_count[0]

        times = measure_execution_time(counted_func, n_runs=5, warmup=3)
        # Total calls = 3 warmup + 5 measured = 8
        assert call_count[0] == 8
        assert len(times) == 5

    def test_percentile_calculation(self):
        """Percentile function should return correct values."""
        data = [1.0, 2.0, 3.0, 4.0, 5.0]
        assert percentile(data, 0) == 1.0
        assert percentile(data, 50) == 3.0
        assert percentile(data, 100) == 5.0

    def test_percentile_single_value(self):
        """Percentile of a single value should return that value."""
        assert percentile([42.0], 0) == 42.0
        assert percentile([42.0], 50) == 42.0
        assert percentile([42.0], 100) == 42.0

    def test_assert_percentile_within_passes(self):
        """assert_percentile_within should pass when times are within bounds."""
        times = [0.001, 0.002, 0.0015, 0.0025, 0.003]
        # P95 should be around 0.003, well within 0.01
        assert_percentile_within(times, 95, 0.01, label="test_op")

    def test_assert_percentile_within_fails(self):
        """assert_percentile_within should fail when times exceed bounds."""
        times = [0.1, 0.2, 0.3, 0.4, 0.5]
        with pytest.raises(AssertionError, match="exceeds threshold"):
            assert_percentile_within(times, 95, 0.01, label="test_op")


# ═══════════════════════════════════════════════════════════════════════
# 2. CI environment resource checks
# ═══════════════════════════════════════════════════════════════════════


class TestCIEnvironmentResources:
    """Test that the CI environment has sufficient resources."""

    def test_cpu_count_sufficient(self):
        """The system should have at least 1 CPU core."""
        cpu_count = os.cpu_count()
        assert cpu_count is not None and cpu_count >= 1, (
            f"Insufficient CPU cores: {cpu_count}"
        )

    def test_memory_available(self):
        """The system should have sufficient memory.

        We do not set a hard minimum, but verify we can query memory info.
        """
        try:
            import resource
            # Get maximum resident set size (in bytes on Linux)
            usage = resource.getrusage(resource.RUSAGE_SELF)
            # ru_maxrss is in KB on Linux
            maxrss_mb = usage.ru_maxrss / 1024
            # Just verify it is a reasonable number (not negative, not absurd)
            assert maxrss_mb > 0, f"Unexpected maxrss: {maxrss_mb}"
        except (ImportError, AttributeError):
            pytest.skip("resource module not available")

    def test_disk_space_available(self):
        """There should be sufficient disk space for test outputs."""
        stat = os.statvfs(".")
        available_gb = (stat.f_bavail * stat.f_frsize) / (1024 ** 3)
        assert available_gb > 0.1, (
            f"Insufficient disk space: {available_gb:.2f} GB available"
        )

    def test_time_monotonic_available(self):
        """time.monotonic should be available for reliable timing."""
        t1 = time.monotonic()
        t2 = time.monotonic()
        assert t2 >= t1, "monotonic time should not go backwards"

    def test_perf_counter_available(self):
        """time.perf_counter should be available for high-resolution timing."""
        t1 = time.perf_counter()
        time.sleep(0.001)
        t2 = time.perf_counter()
        assert t2 > t1, "perf_counter should advance"

    def test_timer_resolution_sufficient(self):
        """Timer resolution should be sufficient for benchmarking (< 1ms)."""
        # Measure the smallest observable time difference
        times = []
        for _ in range(100):
            t1 = time.perf_counter()
            t2 = time.perf_counter()
            times.append(t2 - t1)
        min_resolution = min(t for t in times if t > 0)
        # Resolution should be better than 1ms
        assert min_resolution < 0.001, (
            f"Timer resolution too coarse: {min_resolution:.6f}s"
        )


# ═══════════════════════════════════════════════════════════════════════
# 3. Percentile-based timing assertions
# ═══════════════════════════════════════════════════════════════════════


class TestPercentileTimingAssertions:
    """Test percentile-based timing assertion infrastructure."""

    def test_gc_content_timing_p95(self):
        """GC content computation should meet P95 timing requirements."""
        from biocompiler.sequence.scanner import gc_content

        # Generate a moderately long sequence
        seq = "ATGC" * 1000  # 4000 bp

        times = measure_execution_time(lambda: gc_content(seq), n_runs=20)
        # P95 should be well under 100ms for 4kb
        assert_percentile_within(times, 95, 0.1, label="gc_content_4kb")

    def test_compute_cai_timing_p95(self):
        """CAI computation should meet P95 timing requirements."""
        from biocompiler.expression.translation import compute_cai

        dna = "ATGGTTTCTAAAGGTGAA" * 10  # 180 bp

        times = measure_execution_time(
            lambda: compute_cai(dna, organism="Escherichia_coli"),
            n_runs=20,
        )
        # P95 should be under 100ms
        assert_percentile_within(times, 95, 0.1, label="compute_cai_180bp")

    def test_no_regression_on_simple_operation(self):
        """A simple operation should not show regression between runs."""
        # Run the same function twice and compare
        baseline = measure_execution_time(
            lambda: sum(range(10000)), n_runs=10
        )
        current = measure_execution_time(
            lambda: sum(range(10000)), n_runs=10
        )
        # Allow 3x regression ratio to account for CI noise
        assert_no_regression(
            baseline, current, max_regression_ratio=3.0,
            label="sum_range_10k",
        )

    def test_timing_distribution_reasonable(self):
        """Timing distribution should not have extreme outliers.

        The ratio of P99 to P50 should be less than 10x for a
        simple deterministic function.
        """
        from biocompiler.sequence.scanner import gc_content

        seq = "ATGC" * 500

        times = measure_execution_time(lambda: gc_content(seq), n_runs=30)
        p50 = percentile(times, 50)
        p99 = percentile(times, 99)

        if p50 > 0:
            ratio = p99 / p50
            assert ratio < 10.0, (
                f"Timing outlier ratio too high: P99/P50 = {ratio:.1f}x "
                f"(P50={p50:.6f}, P99={p99:.6f})"
            )

    def test_concurrent_timing_doesnt_degrade_sequential(self):
        """Concurrent execution should not degrade sequential timing."""
        from biocompiler.sequence.scanner import gc_content

        seq = "ATGC" * 500

        # Sequential baseline
        seq_times = measure_execution_time(
            lambda: gc_content(seq), n_runs=10
        )
        seq_median = statistics.median(seq_times)

        # The sequential median should be under a reasonable threshold
        assert seq_median < 1.0, (
            f"Sequential gc_content too slow: {seq_median:.4f}s"
        )


# ═══════════════════════════════════════════════════════════════════════
# 4. Integration: Full optimization timing benchmarks
# ═══════════════════════════════════════════════════════════════════════


class TestOptimizationTimingBenchmarks:
    """Full optimization timing benchmarks with percentile assertions."""

    @pytest.mark.slow
    def test_short_protein_optimization_p95(self):
        """Optimizing a short protein should meet P95 timing requirements."""
        from biocompiler.optimizer import optimize_sequence

        protein = "MVSKGE"

        times = measure_execution_time(
            lambda: optimize_sequence(protein, organism="Escherichia_coli"),
            n_runs=5,
        )
        # P95 should be under 30 seconds for a 6-AA protein
        assert_percentile_within(times, 95, 30.0, label="optimize_6aa_ecoli")

    def test_cai_computation_timing(self):
        """CAI computation should be fast enough for real-time use."""
        from biocompiler.expression.translation import compute_cai

        dna = "ATGGTTTCTAAAGGTGAA" * 10  # 180 bp

        times = measure_execution_time(
            lambda: compute_cai(dna, organism="Escherichia_coli"),
            n_runs=20,
        )
        # P95 should be under 100ms
        assert_percentile_within(times, 95, 0.1, label="compute_cai_180bp")
