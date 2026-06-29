"""
Tests for the publication benchmark suite, reproducibility utilities,
and DNAchisel comparison module.

Covers:
  - PUBLICATION_GENES contain valid protein sequences
  - benchmark_cai_quality runs without errors
  - generate_latex_table produces valid LaTeX
  - capture_environment returns required fields
  - save/load benchmark results round-trips correctly
  - ComparisonResult dataclass construction
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# PUBLICATION_GENES validation
# ---------------------------------------------------------------------------


class TestPublicationGenes:
    """Test that PUBLICATION_GENES are valid proteins."""

    def test_publication_genes_exist(self):
        """PUBLICATION_GENES dict should be non-empty."""
        from biocompiler.benchmarking.publication_benchmark import PUBLICATION_GENES
        assert len(PUBLICATION_GENES) > 0

    def test_publication_genes_are_valid_proteins(self):
        """All PUBLICATION_GENES sequences should contain only standard amino acids."""
        from biocompiler.benchmarking.publication_benchmark import PUBLICATION_GENES
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for name, seq in PUBLICATION_GENES.items():
            invalid_chars = set(seq) - valid_aa
            assert not invalid_chars, (
                f"Gene '{name}' contains non-standard amino acid codes: {invalid_chars}"
            )

    def test_publication_genes_minimum_length(self):
        """All PUBLICATION_GENES should be at least 20 amino acids long."""
        from biocompiler.benchmarking.publication_benchmark import PUBLICATION_GENES
        for name, seq in PUBLICATION_GENES.items():
            assert len(seq) >= 20, (
                f"Gene '{name}' is too short ({len(seq)} aa); minimum is 20."
            )

    def test_publication_genes_has_expected_keys(self):
        """PUBLICATION_GENES should contain the expected gene names."""
        from biocompiler.benchmarking.publication_benchmark import PUBLICATION_GENES
        expected = {"GFP", "Insulin", "HBB", "mCherry"}
        assert expected.issubset(set(PUBLICATION_GENES.keys())), (
            f"Missing expected genes: {expected - set(PUBLICATION_GENES.keys())}"
        )


# ---------------------------------------------------------------------------
# benchmark_cai_quality
# ---------------------------------------------------------------------------


class TestBenchmarkCAIQuality:
    """Test that benchmark_cai_quality runs without errors."""

    def test_benchmark_cai_quality_small_gene_set(self):
        """benchmark_cai_quality should run on a small gene set without errors."""
        from biocompiler.benchmarking.publication_benchmark import benchmark_cai_quality

        # Use a small subset for fast testing
        small_genes = {"GFP": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYQ"}
        organisms = ["Escherichia_coli"]

        results = benchmark_cai_quality(gene_set=small_genes, organisms=organisms)

        # Verify structure
        assert "per_gene" in results
        assert "mean_cai_biocompiler" in results
        assert "organisms" in results
        assert "genes" in results
        assert "GFP" in results["per_gene"]
        assert "Escherichia_coli" in results["per_gene"]["GFP"]
        assert "cai_biocompiler" in results["per_gene"]["GFP"]["Escherichia_coli"]

    def test_benchmark_cai_quality_returns_float_cai(self):
        """CAI values should be floats in [0, 1]."""
        from biocompiler.benchmarking.publication_benchmark import benchmark_cai_quality

        small_genes = {"GFP": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYQ"}
        results = benchmark_cai_quality(gene_set=small_genes, organisms=["Escherichia_coli"])

        cai = results["per_gene"]["GFP"]["Escherichia_coli"]["cai_biocompiler"]
        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0


# ---------------------------------------------------------------------------
# generate_latex_table
# ---------------------------------------------------------------------------


class TestGenerateLatexTable:
    """Test that generate_latex_table produces valid LaTeX."""

    def test_generate_latex_table_cai_type(self):
        """CAI table type should produce valid LaTeX with proper structure."""
        from biocompiler.benchmarking.publication_benchmark import generate_latex_table

        results = {
            "cai_quality": {
                "per_gene": {
                    "GFP": {
                        "Escherichia_coli": {
                            "cai_biocompiler": 0.85,
                            "cai_dnachisel": 0.78,
                            "delta": 0.07,
                        }
                    }
                }
            }
        }
        latex = generate_latex_table(results, table_type="cai")

        assert r"\begin{table}" in latex
        assert r"\end{table}" in latex
        assert r"\toprule" in latex
        assert r"\bottomrule" in latex
        assert "GFP" in latex
        assert "0.85" in latex

    def test_generate_latex_table_speed_type(self):
        """Speed table type should produce valid LaTeX with speed columns."""
        from biocompiler.benchmarking.publication_benchmark import generate_latex_table

        results = {
            "speed": {
                "per_gene": {
                    "GFP": {
                        "Escherichia_coli": {
                            "biocompiler_mean_ms": 45.2,
                            "biocompiler_std_ms": 3.1,
                            "dnachisel_mean_ms": 120.5,
                        }
                    }
                }
            }
        }
        latex = generate_latex_table(results, table_type="speed")

        assert r"\begin{table}" in latex
        assert "GFP" in latex
        assert "45.2" in latex

    def test_generate_latex_table_summary_type(self):
        """Summary table type should produce valid LaTeX with summary metrics."""
        from biocompiler.benchmarking.publication_benchmark import generate_latex_table

        results = {
            "summary": {
                "mean_cai_biocompiler": 0.85,
                "mean_cai_dnachisel": 0.78,
                "cai_advantage": 0.07,
                "mean_speed_biocompiler_ms": 45.0,
                "mean_speed_dnachisel_ms": 120.0,
                "speed_ratio": 0.38,
                "constraint_satisfaction_rate_biocompiler": 0.95,
                "constraint_satisfaction_rate_dnachisel": 0.80,
            }
        }
        latex = generate_latex_table(results, table_type="summary")

        assert r"\begin{table}" in latex
        assert "CAI" in latex
        assert "Speed" in latex

    def test_generate_latex_table_default_type_is_summary(self):
        """Default table_type should produce a summary-style table."""
        from biocompiler.benchmarking.publication_benchmark import generate_latex_table

        results = {
            "summary": {
                "mean_cai_biocompiler": 0.85,
                "mean_cai_dnachisel": 0.78,
                "cai_advantage": 0.07,
                "mean_speed_biocompiler_ms": 45.0,
                "mean_speed_dnachisel_ms": 120.0,
                "speed_ratio": 0.38,
                "constraint_satisfaction_rate_biocompiler": 0.95,
                "constraint_satisfaction_rate_dnachisel": 0.80,
            }
        }
        latex = generate_latex_table(results)  # default table_type
        assert r"\begin{table}" in latex


# ---------------------------------------------------------------------------
# capture_environment
# ---------------------------------------------------------------------------


class TestCaptureEnvironment:
    """Test that capture_environment returns required fields."""

    def test_capture_environment_returns_dict(self):
        """capture_environment should return a dict."""
        from biocompiler.benchmarking.reproducibility import capture_environment
        env = capture_environment()
        assert isinstance(env, dict)

    def test_capture_environment_has_required_fields(self):
        """capture_environment should contain all required fields."""
        from biocompiler.benchmarking.reproducibility import capture_environment
        env = capture_environment()

        required_fields = [
            "python_version",
            "platform_system",
            "platform_machine",
            "cpu_count",
            "package_versions",
            "timestamp",
        ]
        for field_name in required_fields:
            assert field_name in env, f"Missing required field: {field_name}"

    def test_capture_environment_python_version_is_string(self):
        """python_version should be a non-empty string."""
        from biocompiler.benchmarking.reproducibility import capture_environment
        env = capture_environment()
        assert isinstance(env["python_version"], str)
        assert len(env["python_version"]) > 0

    def test_capture_environment_cpu_count_is_positive(self):
        """cpu_count should be a positive integer."""
        from biocompiler.benchmarking.reproducibility import capture_environment
        env = capture_environment()
        assert isinstance(env["cpu_count"], int)
        assert env["cpu_count"] > 0

    def test_capture_environment_package_versions_is_dict(self):
        """package_versions should be a dict."""
        from biocompiler.benchmarking.reproducibility import capture_environment
        env = capture_environment()
        assert isinstance(env["package_versions"], dict)


# ---------------------------------------------------------------------------
# save/load benchmark results
# ---------------------------------------------------------------------------


class TestSaveLoadBenchmarkResults:
    """Test that save/load benchmark results round-trips correctly."""

    def test_save_and_load_roundtrip(self):
        """Saving and loading benchmark results should preserve data."""
        from biocompiler.benchmarking.reproducibility import (
            save_benchmark_results,
            load_benchmark_results,
        )

        original_results = {
            "mean_cai_biocompiler": 0.85,
            "mean_cai_dnachisel": 0.78,
            "mean_speed_biocompiler_ms": 45.0,
            "constraint_satisfaction_rate_biocompiler": 0.95,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_results.json")
            save_benchmark_results(original_results, path)

            # File should exist
            assert Path(path).exists()

            # Load and verify
            loaded = load_benchmark_results(path)
            assert "results" in loaded
            assert "environment" in loaded
            assert loaded["results"]["mean_cai_biocompiler"] == 0.85
            assert loaded["results"]["mean_cai_dnachisel"] == 0.78
            assert loaded["results"]["mean_speed_biocompiler_ms"] == 45.0

    def test_saved_file_contains_environment(self):
        """Saved file should include environment metadata."""
        from biocompiler.benchmarking.reproducibility import (
            save_benchmark_results,
            load_benchmark_results,
        )

        results = {"test_metric": 42.0}

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "env_test.json")
            save_benchmark_results(results, path)

            loaded = load_benchmark_results(path)
            assert "environment" in loaded
            assert "python_version" in loaded["environment"]
            assert "timestamp" in loaded["environment"]

    def test_save_creates_parent_directories(self):
        """save_benchmark_results should create parent directories if needed."""
        from biocompiler.benchmarking.reproducibility import (
            save_benchmark_results,
            load_benchmark_results,
        )

        results = {"test": True}

        with tempfile.TemporaryDirectory() as tmpdir:
            nested_path = os.path.join(tmpdir, "subdir1", "subdir2", "results.json")
            save_benchmark_results(results, nested_path)

            assert Path(nested_path).exists()
            loaded = load_benchmark_results(nested_path)
            assert loaded["results"]["test"] is True


# ---------------------------------------------------------------------------
# validate_benchmark_reproducibility
# ---------------------------------------------------------------------------


class TestValidateBenchmarkReproducibility:
    """Test benchmark reproducibility validation."""

    def test_identical_results_pass(self):
        """Identical results should pass validation."""
        from biocompiler.benchmarking.reproducibility import (
            save_benchmark_results,
            validate_benchmark_reproducibility,
        )

        results = {
            "mean_cai_biocompiler": 0.85,
            "mean_speed_biocompiler_ms": 45.0,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "baseline.json")
            save_benchmark_results(results, path)

            is_reproducible = validate_benchmark_reproducibility(
                path, new_results=results, tolerance=0.01
            )
            assert is_reproducible is True

    def test_within_tolerance_passes(self):
        """Results within tolerance should pass validation."""
        from biocompiler.benchmarking.reproducibility import (
            save_benchmark_results,
            validate_benchmark_reproducibility,
        )

        baseline = {"mean_cai_biocompiler": 0.85}
        new = {"mean_cai_biocompiler": 0.855}  # ~0.6% difference

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "baseline.json")
            save_benchmark_results(baseline, path)

            is_reproducible = validate_benchmark_reproducibility(
                path, new_results=new, tolerance=0.01
            )
            assert is_reproducible is True

    def test_outside_tolerance_fails(self):
        """Results outside tolerance should fail validation."""
        from biocompiler.benchmarking.reproducibility import (
            save_benchmark_results,
            validate_benchmark_reproducibility,
        )

        baseline = {"mean_cai_biocompiler": 0.85}
        new = {"mean_cai_biocompiler": 0.70}  # ~17.6% difference

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "baseline.json")
            save_benchmark_results(baseline, path)

            is_reproducible = validate_benchmark_reproducibility(
                path, new_results=new, tolerance=0.01
            )
            assert is_reproducible is False


# ---------------------------------------------------------------------------
# ComparisonResult dataclass
# ---------------------------------------------------------------------------


class TestComparisonResult:
    """Test ComparisonResult dataclass construction."""

    def test_basic_construction(self):
        """ComparisonResult should be constructable with required fields."""
        from biocompiler.benchmarking.dnachisel_comparison import ComparisonResult

        result = ComparisonResult(
            category="CAI",
            biocompiler_metrics={"mean": 0.85, "wins": 5},
            dnachisel_metrics={"mean": 0.78, "wins": 2},
            winner="biocompiler",
            details="BioCompiler leads on CAI by +0.07.",
        )

        assert result.category == "CAI"
        assert result.biocompiler_metrics["mean"] == 0.85
        assert result.dnachisel_metrics["mean"] == 0.78
        assert result.winner == "biocompiler"
        assert "CAI" in result.details

    def test_default_values(self):
        """ComparisonResult defaults should be sensible."""
        from biocompiler.benchmarking.dnachisel_comparison import ComparisonResult

        result = ComparisonResult(category="Test")
        assert result.biocompiler_metrics == {}
        assert result.dnachisel_metrics == {}
        assert result.winner == "tie"
        assert result.details == ""

    def test_winner_values(self):
        """ComparisonResult winner should accept valid values."""
        from biocompiler.benchmarking.dnachisel_comparison import ComparisonResult

        for winner in ["biocompiler", "dnachisel", "tie"]:
            result = ComparisonResult(
                category="Test",
                winner=winner,
            )
            assert result.winner == winner

    def test_comparison_result_with_eukaryotic_metrics(self):
        """ComparisonResult should work with eukaryotic feature metrics."""
        from biocompiler.benchmarking.dnachisel_comparison import ComparisonResult

        result = ComparisonResult(
            category="Eukaryotic Features",
            biocompiler_metrics={"mean_rate": 0.85, "wins": 4, "n": 8},
            dnachisel_metrics={"mean_rate": 0.55, "wins": 1, "n": 8},
            winner="biocompiler",
            details="BioCompiler handles eukaryotic features better.",
        )

        assert result.category == "Eukaryotic Features"
        assert result.biocompiler_metrics["n"] == 8
        assert result.winner == "biocompiler"


# ---------------------------------------------------------------------------
# generate_comparison_table
# ---------------------------------------------------------------------------


class TestGenerateComparisonTable:
    """Test LaTeX table generation from ComparisonResult objects."""

    def test_generates_valid_latex(self):
        """generate_comparison_table should produce valid LaTeX."""
        from biocompiler.benchmarking.dnachisel_comparison import (
            ComparisonResult,
            generate_comparison_table,
        )

        results = [
            ComparisonResult(
                category="CAI",
                biocompiler_metrics={"mean": 0.85, "wins": 5},
                dnachisel_metrics={"mean": 0.78, "wins": 2},
                winner="biocompiler",
                details="BC leads on CAI.",
            ),
        ]

        latex = generate_comparison_table(results)
        assert r"\begin{table}" in latex
        assert r"\end{table}" in latex
        assert "CAI" in latex
        assert "biocompiler" in latex

    def test_handles_speed_metrics(self):
        """generate_comparison_table should handle speed metrics with ms units."""
        from biocompiler.benchmarking.dnachisel_comparison import (
            ComparisonResult,
            generate_comparison_table,
        )

        results = [
            ComparisonResult(
                category="Speed",
                biocompiler_metrics={"mean_ms": 45.2, "wins": 3},
                dnachisel_metrics={"mean_ms": 120.5, "wins": 1},
                winner="biocompiler",
                details="BC is faster.",
            ),
        ]

        latex = generate_comparison_table(results)
        assert "ms" in latex


# ---------------------------------------------------------------------------
# generate_comparison_report
# ---------------------------------------------------------------------------


class TestGenerateComparisonReport:
    """Test Markdown report generation from ComparisonResult objects."""

    def test_generates_markdown_report(self):
        """generate_comparison_report should produce valid Markdown."""
        from biocompiler.benchmarking.dnachisel_comparison import (
            ComparisonResult,
            generate_comparison_report,
        )

        results = [
            ComparisonResult(
                category="CAI",
                biocompiler_metrics={"mean": 0.85, "wins": 5},
                dnachisel_metrics={"mean": 0.78, "wins": 2},
                winner="biocompiler",
                details="BioCompiler leads on CAI by +0.07.",
            ),
            ComparisonResult(
                category="Speed",
                biocompiler_metrics={"mean_ms": 45.2, "wins": 3},
                dnachisel_metrics={"mean_ms": 120.5, "wins": 1},
                winner="biocompiler",
                details="BioCompiler is 2.7x faster.",
            ),
        ]

        report = generate_comparison_report(results)
        assert "# BioCompiler vs DNAchisel" in report
        assert "## Summary" in report
        assert "CAI" in report
        assert "Speed" in report
        assert "## Details" in report
        assert "## Overall Verdict" in report
        assert "biocompiler" in report.lower()

    def test_report_with_mixed_winners(self):
        """Report should correctly reflect mixed winner categories."""
        from biocompiler.benchmarking.dnachisel_comparison import (
            ComparisonResult,
            generate_comparison_report,
        )

        results = [
            ComparisonResult(
                category="CAI",
                biocompiler_metrics={"mean": 0.85, "wins": 5},
                dnachisel_metrics={"mean": 0.78, "wins": 2},
                winner="biocompiler",
                details="BC leads on CAI.",
            ),
            ComparisonResult(
                category="Speed",
                biocompiler_metrics={"mean_ms": 200.0, "wins": 1},
                dnachisel_metrics={"mean_ms": 45.0, "wins": 5},
                winner="dnachisel",
                details="DC is faster.",
            ),
        ]

        report = generate_comparison_report(results)
        assert "Tie" in report  # 1-1 is a tie


# ---------------------------------------------------------------------------
# benchmark_speed
# ---------------------------------------------------------------------------


class TestBenchmarkSpeed:
    """Test benchmark_speed function."""

    def test_benchmark_speed_returns_stats(self):
        """benchmark_speed should return statistical summary."""
        from biocompiler.benchmarking.publication_benchmark import benchmark_speed

        small_genes = {"GFP": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYQ"}
        results = benchmark_speed(
            gene_set=small_genes,
            organisms=["Escherichia_coli"],
            warmup=1,
            repeats=3,
        )

        assert "per_gene" in results
        assert "mean_speed_biocompiler_ms" in results
        assert "warmup" in results
        assert results["warmup"] == 1
        assert results["repeats"] == 3

        # Check per-gene stats include std, min, max
        gene_stats = results["per_gene"]["GFP"]["Escherichia_coli"]
        assert "biocompiler_mean_ms" in gene_stats
        assert "biocompiler_std_ms" in gene_stats
        assert "biocompiler_min_ms" in gene_stats
        assert "biocompiler_max_ms" in gene_stats


# ---------------------------------------------------------------------------
# benchmark_constraint_satisfaction
# ---------------------------------------------------------------------------


class TestBenchmarkConstraintSatisfaction:
    """Test benchmark_constraint_satisfaction function."""

    def test_benchmark_constraint_satisfaction_returns_rates(self):
        """benchmark_constraint_satisfaction should return satisfaction rates."""
        from biocompiler.benchmarking.publication_benchmark import (
            benchmark_constraint_satisfaction,
        )

        small_genes = {"GFP": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYQ"}
        constraint_sets = [
            {
                "name": "gc_range_only",
                "constraints": [
                    {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                ],
            },
        ]

        results = benchmark_constraint_satisfaction(
            gene_set=small_genes,
            constraint_sets=constraint_sets,
        )

        assert "constraint_satisfaction_rate_biocompiler" in results
        assert "constraint_satisfaction_rate_dnachisel" in results
        assert "constraint_sets_used" in results
        assert "gc_range_only" in results["constraint_sets_used"]
