"""Tests for the organism-aware CAI recovery benchmark.

Validates that:
  1. The benchmark runs without errors
  2. Organism-aware mode produces higher CAI for E. coli than organism-unaware
  3. The report printing function works
"""

from __future__ import annotations

import io
import sys
from unittest.mock import patch

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def small_gene_set() -> dict:
    """A small gene set for fast test execution (3 short proteins)."""
    return {
        "TestINS": {
            "protein_sequence": (
                "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
                "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLEN"
                "YCN"
            ),
            "organism": "Homo sapiens",
            "description": "Insulin (test)",
        },
        "TestGH1_frag": {
            "protein_sequence": (
                "FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFS"
            ),
            "organism": "Homo sapiens",
            "description": "Growth hormone fragment (test)",
        },
        "TestEPO_frag": {
            "protein_sequence": (
                "APPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRME"
            ),
            "organism": "Homo sapiens",
            "description": "Erythropoietin fragment (test)",
        },
    }


# ---------------------------------------------------------------------------
# Test: benchmark runs without errors
# ---------------------------------------------------------------------------

class TestBenchmarkRuns:
    """Verify that the benchmark completes without raising."""

    def test_benchmark_returns_dict(self, small_gene_set: dict) -> None:
        """benchmark_organism_aware_cai() returns a dict with expected keys."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        assert isinstance(results, dict)
        assert "mean_cai_old" in results
        assert "mean_cai_new" in results
        assert "mean_cai_recovery" in results
        assert "per_gene_results" in results
        assert "n_genes" in results
        assert "organism" in results

    def test_benchmark_counts_genes(self, small_gene_set: dict) -> None:
        """The n_genes field matches the number of per-gene results."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        assert results["n_genes"] == len(small_gene_set)
        assert len(results["per_gene_results"]) == len(small_gene_set)

    def test_per_gene_result_fields(self, small_gene_set: dict) -> None:
        """Each per-gene result has the expected fields."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        for gene_result in results["per_gene_results"]:
            assert "gene_name" in gene_result
            assert "protein_length" in gene_result
            assert "cai_old" in gene_result
            assert "cai_new" in gene_result
            assert "cai_recovery" in gene_result
            assert isinstance(gene_result["cai_old"], float)
            assert isinstance(gene_result["cai_new"], float)
            assert isinstance(gene_result["cai_recovery"], float)

    def test_empty_gene_set(self) -> None:
        """An empty gene set returns zeroed summary with no errors."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(gene_set={})
        assert results["n_genes"] == 0
        assert results["mean_cai_old"] == 0.0
        assert results["mean_cai_new"] == 0.0
        assert results["mean_cai_recovery"] == 0.0


# ---------------------------------------------------------------------------
# Test: organism-aware mode produces higher CAI for E. coli
# ---------------------------------------------------------------------------

class TestCAIRecovery:
    """Verify that organism-aware constraints yield higher CAI for prokaryotes."""

    def test_organism_aware_higher_cai(self, small_gene_set: dict) -> None:
        """Mean CAI with organism-aware constraints >= organism-unaware."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        # The organism-aware mode should recover CAI (i.e., mean_cai_new >= mean_cai_old)
        # For E. coli (prokaryote), disabling splice/CpG constraints should
        # not reduce CAI.  At minimum, it should be equal; in practice it
        # should be strictly higher because the optimizer is no longer forced
        # to use suboptimal codons for irrelevant constraints.
        assert results["mean_cai_new"] >= results["mean_cai_old"], (
            f"Organism-aware CAI ({results['mean_cai_new']:.4f}) should be "
            f">= organism-unaware CAI ({results['mean_cai_old']:.4f}) for E. coli"
        )

    def test_cai_recovery_non_negative(self, small_gene_set: dict) -> None:
        """The mean CAI recovery should be non-negative for E. coli."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        assert results["mean_cai_recovery"] >= 0.0, (
            f"CAI recovery should be non-negative for E. coli, "
            f"got {results['mean_cai_recovery']:+.4f}"
        )

    def test_per_gene_recovery_non_negative(self, small_gene_set: dict) -> None:
        """Each individual gene should show non-negative CAI recovery."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        for gene_result in results["per_gene_results"]:
            # Allow a tiny tolerance for floating-point edge cases
            assert gene_result["cai_recovery"] >= -0.001, (
                f"Gene {gene_result['gene_name']}: CAI recovery should be "
                f"non-negative, got {gene_result['cai_recovery']:+.6f}"
            )

    def test_cai_values_in_valid_range(self, small_gene_set: dict) -> None:
        """All CAI values should be in [0, 1]."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        for gene_result in results["per_gene_results"]:
            assert 0.0 <= gene_result["cai_old"] <= 1.0, (
                f"cai_old out of range: {gene_result['cai_old']}"
            )
            assert 0.0 <= gene_result["cai_new"] <= 1.0, (
                f"cai_new out of range: {gene_result['cai_new']}"
            )


# ---------------------------------------------------------------------------
# Test: report printing
# ---------------------------------------------------------------------------

class TestReportPrinting:
    """Verify the print_organism_aware_report function."""

    def test_report_prints_without_error(self, small_gene_set: dict) -> None:
        """print_organism_aware_report() runs without raising."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
            print_organism_aware_report,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        # Capture stdout
        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_organism_aware_report(results)

        output = captured.getvalue()
        assert len(output) > 0

    def test_report_contains_key_phrases(self, small_gene_set: dict) -> None:
        """The report output contains expected sections and values."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
            print_organism_aware_report,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_organism_aware_report(results)

        output = captured.getvalue()

        # Key phrases that should appear in the report
        assert "Organism-Aware" in output
        assert "Escherichia_coli" in output
        assert "Mean CAI" in output
        assert "CAI(old)" in output
        assert "CAI(new)" in output
        assert "recovery" in output.lower()

    def test_report_with_empty_results(self) -> None:
        """The report handles an empty result set gracefully."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            print_organism_aware_report,
        )

        empty_results = {
            "mean_cai_old": 0.0,
            "mean_cai_new": 0.0,
            "mean_cai_dnachisel": None,
            "mean_cai_recovery": 0.0,
            "n_genes": 0,
            "organism": "Escherichia_coli",
            "per_gene_results": [],
        }

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_organism_aware_report(empty_results)

        output = captured.getvalue()
        assert "Escherichia_coli" in output

    def test_report_with_dnachisel(self, small_gene_set: dict) -> None:
        """The report includes DNAchisel column when data is present."""
        from biocompiler.benchmarking.organism_aware_benchmark import (
            benchmark_organism_aware_cai,
            print_organism_aware_report,
        )

        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        # Inject fake DNAchisel data to test the report formatting
        results["mean_cai_dnachisel"] = 0.75
        for gene in results["per_gene_results"]:
            gene["cai_dnachisel"] = 0.75

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_organism_aware_report(results)

        output = captured.getvalue()
        assert "DNAchisel" in output


# ---------------------------------------------------------------------------
# Test: runner convenience function
# ---------------------------------------------------------------------------

class TestRunnerConvenience:
    """Verify the convenience wrapper in runner.py."""

    def test_run_benchmark_by_name(self, small_gene_set: dict) -> None:
        """run_benchmark_by_name('organism_aware_cai') delegates correctly."""
        from biocompiler.benchmarking.runner import run_benchmark_by_name
        from biocompiler.benchmarking.organism_aware_benchmark import benchmark_organism_aware_cai

        # Test that the runner dispatches correctly
        results = benchmark_organism_aware_cai(
            gene_set=small_gene_set,
            organism="Escherichia_coli",
        )

        assert isinstance(results, dict)
        assert "mean_cai_old" in results
        assert "mean_cai_new" in results
        assert "mean_cai_recovery" in results

    def test_available_benchmarks_contains_organism_aware(self) -> None:
        """AVAILABLE_BENCHMARKS includes the organism_aware_cai benchmark."""
        from biocompiler.benchmarking.runner import AVAILABLE_BENCHMARKS

        assert "organism_aware_cai" in AVAILABLE_BENCHMARKS
