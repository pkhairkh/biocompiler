"""Tests for the head-to-head benchmark module.

Tests that:
  1. The benchmark runs without errors (even with DNAchisel absent)
  2. HeadToHeadResult is populated correctly
  3. Report printing produces output without errors
"""

from __future__ import annotations

import io
import sys
from unittest.mock import patch

import pytest

from biocompiler.benchmarking.head_to_head_benchmark import (
    HeadToHeadResult,
    run_head_to_head,
    print_head_to_head_report,
)


# ─── Fixtures ─────────────────────────────────────────────────────────


class _FakeOptimizationResult:
    """Minimal stand-in for biocompiler.optimization.OptimizationResult."""

    def __init__(
        self,
        sequence: str = "ATGAAAGCGTTTTGA",
        cai: float = 0.85,
        gc_content: float = 0.40,
        provenance: object = True,
    ) -> None:
        self.sequence = sequence
        self.cai = cai
        self.gc_content = gc_content
        self.provenance = provenance
        self.failed_predicates: list = []
        self.predicate_results: list = []


# A small test gene panel (2 genes)
_TEST_GENES = ["lacZ", "recA"]

# Patch target for DNAchisel availability (imported inside run_head_to_head)
_DNACHISEL_PATCH = "biocompiler.benchmarking.dnachisel_adapter.is_dnachisel_available"
_OPTIMIZE_PATCH = "biocompiler.optimization.optimize_sequence"


# ─── Test: benchmark runs without errors ─────────────────────────────


class TestBenchmarkRuns:
    """Test that the benchmark runs to completion under various conditions."""

    def test_run_with_no_genes(self) -> None:
        """Benchmark with an empty gene list returns zeroed result."""
        result = run_head_to_head(genes=[], organism="Escherichia_coli")
        assert result.num_genes == 0
        assert result.biocompiler_mean_cai == 0.0
        assert result.dnachisel_mean_cai == 0.0
        assert result.per_gene_results == []

    def test_run_with_nonexistent_genes(self) -> None:
        """Benchmark with gene names not in the panel produces empty results."""
        result = run_head_to_head(
            genes=["NONEXISTENT_GENE_XYZ"], organism="Escherichia_coli",
        )
        # The gene won't be found, so panel is empty
        assert result.num_genes == 0

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_run_with_biocompiler_only(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """Benchmark runs with only BioCompiler (DNAchisel unavailable)."""
        # Return a valid fake result for each gene
        mock_optimize.return_value = _FakeOptimizationResult(
            sequence="ATGAAAGCGTTTTGA",
            cai=0.78,
            gc_content=0.40,
        )

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        assert result.num_genes == 2
        assert result.biocompiler_mean_cai > 0.0
        # DNAchisel should be 0.0 since unavailable
        assert result.dnachisel_mean_cai == 0.0
        assert result.dnachisel_mean_time_ms == 0.0
        assert result.speed_ratio == 0.0
        # Per-gene results should have 2 entries
        assert len(result.per_gene_results) == 2
        # Each entry should have required keys
        for g in result.per_gene_results:
            assert "gene_name" in g
            assert "biocompiler_cai" in g
            assert "biocompiler_gc" in g
            assert "biocompiler_time_ms" in g
            assert "biocompiler_violations" in g
            assert "biocompiler_success" in g
            assert g["biocompiler_success"] is True

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_run_handles_biocompiler_failure(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """Benchmark gracefully handles BioCompiler optimization failure."""
        mock_optimize.side_effect = RuntimeError("optimization failed")

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        assert result.num_genes == 2
        # All genes should have failed
        for g in result.per_gene_results:
            assert g["biocompiler_success"] is False
            assert g["biocompiler_cai"] == 0.0
        # Mean CAI should be 0 since all failed
        assert result.biocompiler_mean_cai == 0.0


# ─── Test: HeadToHeadResult population ────────────────────────────────


class TestHeadToHeadResultPopulation:
    """Test that HeadToHeadResult fields are correctly populated."""

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_result_fields_populated(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """All HeadToHeadResult fields are populated with reasonable values."""
        mock_optimize.return_value = _FakeOptimizationResult(
            sequence="ATGAAAGCGTTTTGA",
            cai=0.82,
            gc_content=0.38,
        )

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        # Required fields exist and have correct types
        assert isinstance(result.organism, str)
        assert result.organism == "Escherichia_coli"
        assert isinstance(result.num_genes, int)
        assert result.num_genes == 2
        assert isinstance(result.biocompiler_mean_cai, float)
        assert isinstance(result.dnachisel_mean_cai, float)
        assert isinstance(result.biocompiler_mean_time_ms, float)
        assert isinstance(result.dnachisel_mean_time_ms, float)
        assert isinstance(result.speed_ratio, float)
        assert isinstance(result.cai_gap, float)
        assert isinstance(result.per_gene_results, list)
        assert isinstance(result.biocompiler_provenance_available, bool)
        assert isinstance(result.dnachisel_provenance_available, bool)

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_cai_gap_calculation(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """CAI gap is correctly computed as BC - DC."""
        mock_optimize.return_value = _FakeOptimizationResult(cai=0.90)

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        # DNAchisel unavailable → mean_cai=0 → gap = BC_mean
        assert result.cai_gap == result.biocompiler_mean_cai

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_per_gene_structure(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """Each per-gene result dict has the expected keys and value types."""
        mock_optimize.return_value = _FakeOptimizationResult(cai=0.75)

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        expected_keys = {
            "gene_name", "protein_length", "biocompiler_cai",
            "dnachisel_cai", "biocompiler_gc", "dnachisel_gc",
            "biocompiler_time_ms", "dnachisel_time_ms",
            "biocompiler_violations", "dnachisel_violations",
            "biocompiler_success", "dnachisel_success", "sequence_length",
        }

        for g in result.per_gene_results:
            assert expected_keys <= set(g.keys()), (
                f"Missing keys: {expected_keys - set(g.keys())}"
            )
            assert isinstance(g["gene_name"], str)
            assert isinstance(g["protein_length"], int)
            assert isinstance(g["biocompiler_cai"], float)
            assert isinstance(g["biocompiler_time_ms"], float)
            assert isinstance(g["biocompiler_success"], bool)
            # DNAchisel fields should be None since unavailable
            assert g["dnachisel_cai"] is None
            assert g["dnachisel_gc"] is None
            assert g["dnachisel_time_ms"] is None
            assert g["dnachisel_violations"] is None
            assert g["dnachisel_success"] is None

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_provenance_tracking(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """Provenance availability is correctly reported."""
        mock_optimize.return_value = _FakeOptimizationResult(
            provenance={"decisions": ["codon_choice_1"]},
        )

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        assert result.biocompiler_provenance_available is True
        assert result.dnachisel_provenance_available is False

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_mean_cai_calculation(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """Mean CAI is correctly averaged across successful genes."""
        call_count = [0]

        def side_effect(**kwargs):
            call_count[0] += 1
            return _FakeOptimizationResult(cai=0.80 + call_count[0] * 0.05)

        mock_optimize.side_effect = side_effect

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        # Gene 1: cai=0.85, Gene 2: cai=0.90, mean=0.875
        expected_mean = 0.875
        assert abs(result.biocompiler_mean_cai - expected_mean) < 0.01

    @patch(_DNACHISEL_PATCH, return_value=False)
    @patch(_OPTIMIZE_PATCH)
    def test_constraint_violations_counted(
        self,
        mock_optimize: patch,
        mock_dnachisel_avail: patch,
    ) -> None:
        """Constraint violations are counted per gene."""
        # Use a short valid sequence
        mock_optimize.return_value = _FakeOptimizationResult(
            sequence="ATGAAAGCGTTTTGA",
            cai=0.75,
            gc_content=0.38,
        )

        result = run_head_to_head(
            genes=_TEST_GENES, organism="Escherichia_coli",
        )

        for g in result.per_gene_results:
            assert isinstance(g["biocompiler_violations"], int)
            assert g["biocompiler_violations"] >= 0


# ─── Test: Report printing ────────────────────────────────────────────


class TestReportPrinting:
    """Test that print_head_to_head_report produces output without errors."""

    def _make_result(self, **overrides) -> HeadToHeadResult:
        """Create a HeadToHeadResult with sensible defaults."""
        defaults = dict(
            organism="Escherichia_coli",
            num_genes=2,
            biocompiler_mean_cai=0.85,
            dnachisel_mean_cai=0.0,
            biocompiler_mean_time_ms=12.5,
            dnachisel_mean_time_ms=0.0,
            speed_ratio=0.0,
            cai_gap=0.85,
            per_gene_results=[
                {
                    "gene_name": "lacZ",
                    "protein_length": 1024,
                    "biocompiler_cai": 0.82,
                    "dnachisel_cai": None,
                    "biocompiler_gc": 0.42,
                    "dnachisel_gc": None,
                    "biocompiler_time_ms": 15.0,
                    "dnachisel_time_ms": None,
                    "biocompiler_violations": 0,
                    "dnachisel_violations": None,
                    "biocompiler_success": True,
                    "dnachisel_success": None,
                    "sequence_length": 3072,
                },
                {
                    "gene_name": "recA",
                    "protein_length": 353,
                    "biocompiler_cai": 0.88,
                    "dnachisel_cai": None,
                    "biocompiler_gc": 0.45,
                    "dnachisel_gc": None,
                    "biocompiler_time_ms": 10.0,
                    "dnachisel_time_ms": None,
                    "biocompiler_violations": 0,
                    "dnachisel_violations": None,
                    "biocompiler_success": True,
                    "dnachisel_success": None,
                    "sequence_length": 1059,
                },
            ],
            biocompiler_provenance_available=True,
            dnachisel_provenance_available=False,
        )
        defaults.update(overrides)
        return HeadToHeadResult(**defaults)

    def test_print_biocompiler_only(self) -> None:
        """Report prints successfully for BioCompiler-only results."""
        result = self._make_result()

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_head_to_head_report(result)

        output = captured.getvalue()
        assert "HEAD-TO-HEAD" in output
        assert "Escherichia_coli" in output
        assert "BioCompiler" in output
        assert "PROVENANCE" in output

    def test_print_with_dnachisel_results(self) -> None:
        """Report prints successfully when DNAchisel results are present."""
        per_gene = [
            {
                "gene_name": "lacZ",
                "protein_length": 1024,
                "biocompiler_cai": 0.82,
                "dnachisel_cai": 0.78,
                "biocompiler_gc": 0.42,
                "dnachisel_gc": 0.44,
                "biocompiler_time_ms": 15.0,
                "dnachisel_time_ms": 25.0,
                "biocompiler_violations": 0,
                "dnachisel_violations": 0,
                "biocompiler_success": True,
                "dnachisel_success": True,
                "sequence_length": 3072,
            },
        ]
        result = self._make_result(
            num_genes=1,
            biocompiler_mean_cai=0.82,
            dnachisel_mean_cai=0.78,
            biocompiler_mean_time_ms=15.0,
            dnachisel_mean_time_ms=25.0,
            speed_ratio=0.6,
            cai_gap=0.04,
            per_gene_results=per_gene,
            dnachisel_provenance_available=True,
        )

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_head_to_head_report(result)

        output = captured.getvalue()
        assert "DNAchisel" in output
        assert "0.78" in output
        assert "VERDICT" in output

    def test_print_empty_results(self) -> None:
        """Report handles empty results gracefully."""
        result = self._make_result(
            num_genes=0,
            per_gene_results=[],
        )

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_head_to_head_report(result)

        output = captured.getvalue()
        assert "Genes: 0" in output

    def test_print_bc_leads_on_cai(self) -> None:
        """Report highlights BioCompiler leading on CAI."""
        per_gene = [
            {
                "gene_name": "groEL",
                "protein_length": 547,
                "biocompiler_cai": 0.90,
                "dnachisel_cai": 0.80,
                "biocompiler_gc": 0.48,
                "dnachisel_gc": 0.46,
                "biocompiler_time_ms": 8.0,
                "dnachisel_time_ms": 20.0,
                "biocompiler_violations": 0,
                "dnachisel_violations": 0,
                "biocompiler_success": True,
                "dnachisel_success": True,
                "sequence_length": 1641,
            },
        ]
        result = self._make_result(
            num_genes=1,
            biocompiler_mean_cai=0.90,
            dnachisel_mean_cai=0.80,
            cai_gap=0.10,
            per_gene_results=per_gene,
        )

        captured = io.StringIO()
        with patch("sys.stdout", captured):
            print_head_to_head_report(result)

        output = captured.getvalue()
        assert "BC leads" in output
