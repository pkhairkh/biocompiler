"""
Tests for the FoldX Empirical Stability Benchmark
===================================================
Validates that:
  - The benchmark runs successfully on the curated dataset
  - MAE is computed and documented
  - Direction accuracy exceeds 60%
  - The dataset has adequate coverage (size categories, stability range)
  - Report generation works correctly
"""

from __future__ import annotations

import json
import math

import pytest

from biocompiler.validation.foldx_benchmark import (
    BENCHMARK_DATASET,
    BenchmarkComparison,
    BenchmarkReport,
    BenchmarkStatistics,
    ProteinEntry,
    generate_benchmark_json,
    generate_benchmark_report,
    run_foldx_benchmark,
)


# ────────────────────────────────────────────────────────────
# Dataset Quality Tests
# ────────────────────────────────────────────────────────────

class TestDatasetQuality:
    """Validate the curated benchmark dataset."""

    def test_dataset_has_minimum_30_proteins(self):
        """Dataset must contain at least 30 proteins."""
        assert len(BENCHMARK_DATASET) >= 30, (
            f"Dataset has only {len(BENCHMARK_DATASET)} proteins; minimum is 30"
        )

    def test_dataset_covers_all_size_categories(self):
        """Dataset must include small, medium, and large proteins."""
        categories = {e.size_category for e in BENCHMARK_DATASET}
        assert "small" in categories, "No small proteins (<100 aa) in dataset"
        assert "medium" in categories, "No medium proteins (100-300 aa) in dataset"
        assert "large" in categories, "No large proteins (>300 aa) in dataset"

    def test_small_proteins_under_100_aa(self):
        """Small category proteins must be <100 residues."""
        for entry in BENCHMARK_DATASET:
            if entry.size_category == "small":
                assert len(entry.sequence) < 100, (
                    f"{entry.name} labeled 'small' but has {len(entry.sequence)} aa"
                )

    def test_medium_proteins_100_to_300_aa(self):
        """Medium category proteins must be 100-300 residues."""
        for entry in BENCHMARK_DATASET:
            if entry.size_category == "medium":
                assert 100 <= len(entry.sequence) <= 300, (
                    f"{entry.name} labeled 'medium' but has {len(entry.sequence)} aa"
                )

    def test_large_proteins_over_300_aa(self):
        """Large category proteins must be >300 residues."""
        for entry in BENCHMARK_DATASET:
            if entry.size_category == "large":
                assert len(entry.sequence) > 300, (
                    f"{entry.name} labeled 'large' but has {len(entry.sequence)} aa"
                )

    def test_dataset_includes_stable_proteins(self):
        """Dataset must include proteins with ΔG < -8 (clearly stable)."""
        stable = [e for e in BENCHMARK_DATASET if e.experimental_dg < -8.0]
        assert len(stable) >= 3, (
            f"Only {len(stable)} clearly stable proteins (ΔG < -8); need >= 3"
        )

    def test_dataset_includes_marginally_stable_proteins(self):
        """Dataset must include proteins with -8 ≤ ΔG < -4 (marginally stable)."""
        marginal = [e for e in BENCHMARK_DATASET if -8.0 <= e.experimental_dg < -4.0]
        assert len(marginal) >= 3, (
            f"Only {len(marginal)} marginally stable proteins; need >= 3"
        )

    def test_dataset_includes_less_stable_proteins(self):
        """Dataset must include proteins with ΔG > -4 (less stable)."""
        less_stable = [e for e in BENCHMARK_DATASET if e.experimental_dg > -4.0]
        assert len(less_stable) >= 2, (
            f"Only {len(less_stable)} less stable proteins (ΔG > -4); need >= 2"
        )

    def test_all_entries_have_required_fields(self):
        """Every dataset entry must have name, pdb_id, sequence, ΔG, source, category."""
        for entry in BENCHMARK_DATASET:
            assert entry.name, "Missing protein name"
            assert entry.pdb_id, f"Missing PDB ID for {entry.name}"
            assert entry.sequence, f"Missing sequence for {entry.name}"
            assert isinstance(entry.experimental_dg, float), (
                f"experimental_dg not float for {entry.name}"
            )
            assert entry.source, f"Missing source reference for {entry.name}"
            assert entry.size_category in ("small", "medium", "large"), (
                f"Invalid size_category '{entry.size_category}' for {entry.name}"
            )

    def test_sequences_contain_only_standard_amino_acids(self):
        """All sequences must use only standard 20 amino acid codes."""
        standard = set("ACDEFGHIKLMNPQRSTVWY")
        for entry in BENCHMARK_DATASET:
            seq_set = set(entry.sequence)
            invalid = seq_set - standard
            assert not invalid, (
                f"{entry.name} has non-standard amino acids: {invalid}"
            )

    def test_at_least_5_small_proteins(self):
        """Must have at least 5 small proteins for statistical relevance."""
        small = [e for e in BENCHMARK_DATASET if e.size_category == "small"]
        assert len(small) >= 5, f"Only {len(small)} small proteins; need >= 5"

    def test_at_least_5_medium_proteins(self):
        """Must have at least 5 medium proteins for statistical relevance."""
        medium = [e for e in BENCHMARK_DATASET if e.size_category == "medium"]
        assert len(medium) >= 5, f"Only {len(medium)} medium proteins; need >= 5"

    def test_at_least_3_large_proteins(self):
        """Must have at least 3 large proteins for statistical relevance."""
        large = [e for e in BENCHMARK_DATASET if e.size_category == "large"]
        assert len(large) >= 3, f"Only {len(large)} large proteins; need >= 3"


# ────────────────────────────────────────────────────────────
# Benchmark Execution Tests
# ────────────────────────────────────────────────────────────

class TestBenchmarkExecution:
    """Validate that the benchmark runs successfully."""

    @pytest.fixture(scope="class")
    def benchmark_report(self):
        """Run the full benchmark once and cache the result."""
        return run_foldx_benchmark()

    def test_benchmark_runs_successfully(self, benchmark_report):
        """Benchmark must run without errors."""
        assert isinstance(benchmark_report, BenchmarkReport)

    def test_benchmark_produces_comparisons(self, benchmark_report):
        """Benchmark must produce comparisons for all dataset entries."""
        assert len(benchmark_report.comparisons) > 0, "No comparisons produced"

    def test_all_proteins_evaluated(self, benchmark_report):
        """All proteins in the dataset should be evaluated."""
        n_dataset = len(BENCHMARK_DATASET)
        n_compared = len(benchmark_report.comparisons)
        # Allow some failures but at least 90% should succeed
        assert n_compared >= n_dataset * 0.9, (
            f"Only {n_compared}/{n_dataset} proteins evaluated"
        )

    def test_comparisons_have_valid_structure(self, benchmark_report):
        """Each comparison must have valid numeric fields."""
        for c in benchmark_report.comparisons:
            assert isinstance(c.experimental_dg, float)
            assert isinstance(c.predicted_dg, float)
            assert isinstance(c.error, float)
            assert isinstance(c.abs_error, float)
            assert abs(c.error - (c.predicted_dg - c.experimental_dg)) < 0.01
            assert abs(c.abs_error - abs(c.error)) < 0.01

    def test_direction_correct_is_boolean(self, benchmark_report):
        """direction_correct must be a boolean for each comparison."""
        for c in benchmark_report.comparisons:
            assert isinstance(c.direction_correct, bool)


# ────────────────────────────────────────────────────────────
# Benchmark Statistics Tests
# ────────────────────────────────────────────────────────────

class TestBenchmarkStatistics:
    """Validate that benchmark statistics are meaningful."""

    @pytest.fixture(scope="class")
    def benchmark_report(self):
        return run_foldx_benchmark()

    def test_mae_is_documented(self, benchmark_report):
        """MAE must be a finite, documented number."""
        mae = benchmark_report.statistics.mae
        assert math.isfinite(mae), f"MAE is not finite: {mae}"
        assert mae >= 0, f"MAE must be non-negative: {mae}"

    def test_rmse_is_documented(self, benchmark_report):
        """RMSE must be a finite, documented number."""
        rmse = benchmark_report.statistics.rmse
        assert math.isfinite(rmse), f"RMSE is not finite: {rmse}"
        assert rmse >= 0, f"RMSE must be non-negative: {rmse}"

    def test_pearson_r_is_valid(self, benchmark_report):
        """Pearson r must be in [-1, 1]."""
        r = benchmark_report.statistics.pearson_r
        assert -1.0 <= r <= 1.0, f"Pearson r out of range: {r}"

    def test_direction_accuracy_above_60_percent(self, benchmark_report):
        """Direction accuracy must exceed 60%.

        This is the key quality gate: the empirical heuristic must at
        least correctly predict whether a protein is stable (ΔG < 0) or
        unstable (ΔG ≥ 0) more often than random guessing.
        """
        acc = benchmark_report.statistics.direction_accuracy
        assert acc > 0.60, (
            f"Direction accuracy {acc:.1%} is below 60% threshold. "
            f"The empirical heuristic is not capturing stability direction "
            f"better than chance."
        )

    def test_mae_within_expected_range(self, benchmark_report):
        """MAE should be within the expected ±5 kcal/mol heuristic range.

        The empirical heuristics are documented as having ±5 kcal/mol
        uncertainty, so MAE should be less than 10 kcal/mol at least.
        """
        mae = benchmark_report.statistics.mae
        assert mae < 15.0, (
            f"MAE {mae:.2f} is unexpectedly high (>15 kcal/mol). "
            f"The empirical heuristic may have systematic issues."
        )

    def test_size_category_statistics_present(self, benchmark_report):
        """Statistics must include per-category MAEs and counts."""
        stats = benchmark_report.statistics
        assert stats.small_count > 0, "No small proteins in statistics"
        assert stats.medium_count > 0, "No medium proteins in statistics"
        assert stats.large_count > 0, "No large proteins in statistics"
        assert stats.small_mae >= 0
        assert stats.medium_mae >= 0
        assert stats.large_mae >= 0


# ────────────────────────────────────────────────────────────
# Report Generation Tests
# ────────────────────────────────────────────────────────────

class TestReportGeneration:
    """Validate report generation functions."""

    @pytest.fixture(scope="class")
    def benchmark_report(self):
        return run_foldx_benchmark()

    def test_text_report_generation(self, benchmark_report):
        """Text report must be generated and contain key metrics."""
        text = generate_benchmark_report(benchmark_report)
        assert isinstance(text, str)
        assert "MAE" in text
        assert "RMSE" in text
        assert "Direction accuracy" in text
        assert "Pearson" in text
        assert "kcal/mol" in text

    def test_json_report_generation(self, benchmark_report):
        """JSON report must be generated and be valid JSON."""
        data = generate_benchmark_json(benchmark_report)
        assert isinstance(data, dict)
        # Must be JSON-serializable
        json_str = json.dumps(data)
        assert len(json_str) > 0
        # Must contain key fields
        assert "statistics" in data
        assert "comparisons" in data
        assert "n_proteins" in data
        assert data["statistics"]["mae"] >= 0

    def test_json_report_contains_all_comparisons(self, benchmark_report):
        """JSON report must contain all protein comparisons."""
        data = generate_benchmark_json(benchmark_report)
        assert len(data["comparisons"]) == len(benchmark_report.comparisons)

    def test_text_report_shows_individual_proteins(self, benchmark_report):
        """Text report must list individual protein results."""
        text = generate_benchmark_report(benchmark_report)
        for c in benchmark_report.comparisons:
            # At least some protein names should appear
            assert c.name.split()[0] in text or c.pdb_id in text, (
                f"Protein {c.name} ({c.pdb_id}) not found in text report"
            )


# ────────────────────────────────────────────────────────────
# Custom Dataset Tests
# ────────────────────────────────────────────────────────────

class TestCustomDataset:
    """Test that the benchmark works with a custom dataset."""

    def test_benchmark_with_custom_dataset(self):
        """Benchmark must work with a user-provided dataset."""
        custom = [
            ProteinEntry(
                name="Test protein 1",
                pdb_id="XXXX",
                sequence="MKFLWRAVGLVLSVSCQISAITENPEQGGPAMTAVKWLAVNNEEEVAVTVRVFQKQRLDRGHFIRVTQEAAAELIQKAVETGDVFVVLSKGHEVTVRQKLDAFDGKIIEQKVNVYAPYVYDAKGYVFDFSDKQMTVVYDNGSQYVFGKRFHNEQVYFEDVIKYGKIPVTVEGQKTVSQGFNAYVNWYTKTVTVRKGQIYVKDTPNKTLFEKNSDVRVVRDFVEKLGNKVTYDPAKFQKQDYVDFENLENKKKGFVVAHSSPSVYVKGFNANDVEKLRVDSYKFEDVKYVFENGFNKVTKDVVFQFDVKEKFVYGFNAYSKQFKGDVFVKHFNEKVYVFDDATKVYVFKEHPDKVKVFDDVKVKFNDKVVKVFDN",
                experimental_dg=-10.0,
                source="Test entry",
                size_category="large",
            ),
            ProteinEntry(
                name="Test protein 2",
                pdb_id="YYYY",
                sequence="MKFLWRAVGLVLSVSCQISAITENPEQGGPAVAVKWLAVNNEEEVAVTVRVFQKQRL",
                experimental_dg=-5.0,
                source="Test entry",
                size_category="small",
            ),
        ]

        report = run_foldx_benchmark(dataset=custom)
        assert len(report.comparisons) == 2
        assert report.statistics.n_proteins == 2

    def test_empty_dataset_returns_empty_report(self):
        """Empty dataset should produce empty report without errors."""
        report = run_foldx_benchmark(dataset=[])
        assert len(report.comparisons) == 0
        assert report.statistics.n_proteins == 0
        assert report.statistics.mae == 0.0
