"""
Comprehensive tests for the wet-lab validation framework.

Tests cover:
- ExperimentalResult construction, serialization, deserialization
- ValidationComparison construction and serialization
- WetLabValidator: add/remove experimental results
- WetLabValidator: compare_with_prediction
- WetLabValidator: compute_correlation (Pearson)
- WetLabValidator: compute_rank_correlation (Spearman)
- WetLabValidator: validation_report
- WetLabValidator: load_from_csv
- WetLabValidator: save_report
- Edge cases (empty data, single point, duplicate prevention)
"""

import json
import os
import tempfile

import pytest

from biocompiler.validation.wetlab_validation import (
    ExperimentalResult,
    ValidationComparison,
    WetLabValidator,
)
from biocompiler.optimizer import OptimizationResult


# ═══════════════════════════════════════════════════════════════════════
# ExperimentalResult Tests
# ═══════════════════════════════════════════════════════════════════════

class TestExperimentalResult:
    """Test ExperimentalResult dataclass."""

    def test_construction_basic(self):
        result = ExperimentalResult(
            gene_name="INS",
            organism="Escherichia_coli",
            measured_expression_level=1.2e6,
        )
        assert result.gene_name == "INS"
        assert result.organism == "Escherichia_coli"
        assert result.measured_expression_level == 1.2e6
        assert result.measured_cai is None
        assert result.sequence_used == ""
        assert result.notes == ""

    def test_construction_full(self):
        result = ExperimentalResult(
            gene_name="GFP",
            organism="Escherichia_coli",
            measured_expression_level=5.0e5,
            measured_cai=0.95,
            sequence_used="ATGGCTAGCAAAGGAGAAGAACTTTTCACTGG",
            notes="Shake flask, 37°C, 16h",
        )
        assert result.measured_cai == 0.95
        assert result.sequence_used.startswith("ATG")
        assert result.notes == "Shake flask, 37°C, 16h"

    def test_to_dict(self):
        result = ExperimentalResult(
            gene_name="INS",
            organism="Escherichia_coli",
            measured_expression_level=1.2e6,
            measured_cai=0.88,
        )
        d = result.to_dict()
        assert d["gene_name"] == "INS"
        assert d["organism"] == "Escherichia_coli"
        assert d["measured_expression_level"] == 1.2e6
        assert d["measured_cai"] == 0.88

    def test_from_dict(self):
        data = {
            "gene_name": "GFP",
            "organism": "Saccharomyces_cerevisiae",
            "measured_expression_level": 3.0e5,
            "measured_cai": None,
            "sequence_used": "ATGAAACGT",
            "notes": "Test condition",
        }
        result = ExperimentalResult.from_dict(data)
        assert result.gene_name == "GFP"
        assert result.measured_cai is None
        assert result.notes == "Test condition"

    def test_from_dict_missing_optional_fields(self):
        data = {
            "gene_name": "GFP",
            "organism": "Escherichia_coli",
            "measured_expression_level": 1.0,
        }
        result = ExperimentalResult.from_dict(data)
        assert result.sequence_used == ""
        assert result.notes == ""

    def test_round_trip_serialization(self):
        result = ExperimentalResult(
            gene_name="EPO",
            organism="Homo_sapiens",
            measured_expression_level=2.5e4,
            measured_cai=0.72,
            sequence_used="ATGGAG...",
            notes="HEK293T cells",
        )
        d = result.to_dict()
        restored = ExperimentalResult.from_dict(d)
        assert restored.gene_name == result.gene_name
        assert restored.organism == result.organism
        assert restored.measured_expression_level == result.measured_expression_level
        assert restored.measured_cai == result.measured_cai
        assert restored.sequence_used == result.sequence_used
        assert restored.notes == result.notes


# ═══════════════════════════════════════════════════════════════════════
# ValidationComparison Tests
# ═══════════════════════════════════════════════════════════════════════

class TestValidationComparison:
    """Test ValidationComparison dataclass."""

    def test_construction(self):
        comp = ValidationComparison(
            gene_name="INS",
            predicted_cai=0.95,
            measured_expression=1.2e6,
            correlation=0.85,
            rank_order_match=True,
        )
        assert comp.gene_name == "INS"
        assert comp.predicted_cai == 0.95
        assert comp.measured_expression == 1.2e6
        assert comp.correlation == 0.85
        assert comp.rank_order_match is True

    def test_to_dict(self):
        comp = ValidationComparison(
            gene_name="GFP",
            predicted_cai=0.88,
            measured_expression=5.0e5,
            correlation=-0.3,
            rank_order_match=False,
        )
        d = comp.to_dict()
        assert d["gene_name"] == "GFP"
        assert d["rank_order_match"] is False

    def test_from_dict(self):
        data = {
            "gene_name": "EPO",
            "predicted_cai": 0.72,
            "measured_expression": 2.5e4,
            "correlation": 1.0,
            "rank_order_match": True,
        }
        comp = ValidationComparison.from_dict(data)
        assert comp.gene_name == "EPO"
        assert comp.predicted_cai == 0.72

    def test_round_trip(self):
        comp = ValidationComparison(
            gene_name="HGH",
            predicted_cai=0.90,
            measured_expression=8.0e5,
            correlation=0.95,
            rank_order_match=True,
        )
        d = comp.to_dict()
        restored = ValidationComparison.from_dict(d)
        assert restored.gene_name == comp.gene_name
        assert restored.correlation == comp.correlation


# ═══════════════════════════════════════════════════════════════════════
# WetLabValidator Core Tests
# ═══════════════════════════════════════════════════════════════════════

class TestWetLabValidatorAddRemove:
    """Test adding and removing experimental results."""

    def test_add_result(self):
        validator = WetLabValidator()
        result = ExperimentalResult("INS", "Escherichia_coli", 1.2e6)
        validator.add_experimental_result(result)
        assert len(validator.results) == 1
        assert validator.results[0].gene_name == "INS"

    def test_add_multiple_results(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        validator.add_experimental_result(ExperimentalResult("GFP", "Escherichia_coli", 5.0e5))
        validator.add_experimental_result(ExperimentalResult("EPO", "Homo_sapiens", 2.5e4))
        assert len(validator.results) == 3

    def test_prevent_duplicate_gene_organism(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        with pytest.raises(ValueError, match="already exists"):
            validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 2.0e6))

    def test_allow_same_gene_different_organism(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        validator.add_experimental_result(ExperimentalResult("INS", "Homo_sapiens", 5.0e4))
        assert len(validator.results) == 2

    def test_remove_result(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        validator.add_experimental_result(ExperimentalResult("GFP", "Escherichia_coli", 5.0e5))
        assert validator.remove_experimental_result("INS", "Escherichia_coli") is True
        assert len(validator.results) == 1
        assert validator.results[0].gene_name == "GFP"

    def test_remove_nonexistent_result(self):
        validator = WetLabValidator()
        assert validator.remove_experimental_result("NONEXISTENT", "Escherichia_coli") is False

    def test_results_property_returns_copy(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        results = validator.results
        results.clear()  # Should not affect internal state
        assert len(validator.results) == 1


class TestWetLabValidatorCompare:
    """Test compare_with_prediction method."""

    def _make_opt_result(self, cai: float = 0.95) -> OptimizationResult:
        """Create a mock OptimizationResult with a given CAI."""
        return OptimizationResult(
            sequence="ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTGGTCC",
            gc_content=0.50,
            cai=cai,
        )

    def test_compare_with_prediction(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        opt = self._make_opt_result(cai=0.95)
        comp = validator.compare_with_prediction(opt, "INS")
        assert comp.gene_name == "INS"
        assert comp.predicted_cai == 0.95
        assert comp.measured_expression == 1.2e6

    def test_compare_raises_for_unknown_gene(self):
        validator = WetLabValidator()
        opt = self._make_opt_result()
        with pytest.raises(ValueError, match="No experimental result found"):
            validator.compare_with_prediction(opt, "UNKNOWN_GENE")

    def test_compare_updates_existing_comparison(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        opt1 = self._make_opt_result(cai=0.80)
        comp1 = validator.compare_with_prediction(opt1, "INS")
        assert comp1.predicted_cai == 0.80

        # Second comparison should update the existing one
        opt2 = self._make_opt_result(cai=0.95)
        comp2 = validator.compare_with_prediction(opt2, "INS")
        assert comp2.predicted_cai == 0.95
        assert len(validator.comparisons) == 1  # Should not duplicate


class TestWetLabValidatorCorrelation:
    """Test Pearson and Spearman correlation methods."""

    def _make_opt_result(self, cai: float = 0.95) -> OptimizationResult:
        return OptimizationResult(
            sequence="ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTGGTCC",
            gc_content=0.50,
            cai=cai,
        )

    def test_correlation_with_no_comparisons_raises(self):
        validator = WetLabValidator()
        with pytest.raises(ValueError):
            validator.compute_correlation()

    def test_correlation_single_point(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        opt = self._make_opt_result(cai=0.95)
        validator.compare_with_prediction(opt, "INS")
        r, p = validator.compute_correlation()
        assert r == 1.0  # Single point trivially correlated

    def test_correlation_two_points_positive(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        validator.add_experimental_result(ExperimentalResult("GFP", "Escherichia_coli", 5.0e5))
        opt1 = self._make_opt_result(cai=0.95)
        opt2 = self._make_opt_result(cai=0.80)
        validator.compare_with_prediction(opt1, "INS")
        validator.compare_with_prediction(opt2, "GFP")
        r, p = validator.compute_correlation()
        assert r == 1.0  # Higher CAI → higher expression

    def test_correlation_two_points_negative(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 5.0e5))
        validator.add_experimental_result(ExperimentalResult("GFP", "Escherichia_coli", 1.2e6))
        opt1 = self._make_opt_result(cai=0.95)
        opt2 = self._make_opt_result(cai=0.80)
        validator.compare_with_prediction(opt1, "INS")
        validator.compare_with_prediction(opt2, "GFP")
        r, p = validator.compute_correlation()
        assert r == -1.0  # Higher CAI → lower expression (anti-correlated)

    def test_correlation_three_points(self):
        validator = WetLabValidator()
        # Create data with a positive correlation
        validator.add_experimental_result(ExperimentalResult("A", "Escherichia_coli", 100))
        validator.add_experimental_result(ExperimentalResult("B", "Escherichia_coli", 200))
        validator.add_experimental_result(ExperimentalResult("C", "Escherichia_coli", 300))

        for gene, cai in [("A", 0.70), ("B", 0.85), ("C", 0.95)]:
            opt = self._make_opt_result(cai=cai)
            validator.compare_with_prediction(opt, gene)

        r, p = validator.compute_correlation()
        assert r > 0.9  # Strong positive correlation

    def test_rank_correlation_with_no_comparisons_raises(self):
        validator = WetLabValidator()
        with pytest.raises(ValueError):
            validator.compute_rank_correlation()

    def test_rank_correlation_perfect(self):
        validator = WetLabValidator()
        # Perfect rank ordering
        validator.add_experimental_result(ExperimentalResult("A", "Escherichia_coli", 100))
        validator.add_experimental_result(ExperimentalResult("B", "Escherichia_coli", 200))
        validator.add_experimental_result(ExperimentalResult("C", "Escherichia_coli", 300))

        for gene, cai in [("A", 0.70), ("B", 0.85), ("C", 0.95)]:
            opt = self._make_opt_result(cai=cai)
            validator.compare_with_prediction(opt, gene)

        rho, p = validator.compute_rank_correlation()
        assert rho == 1.0  # Perfect rank correlation

    def test_rank_correlation_reversed(self):
        validator = WetLabValidator()
        # Reversed rank ordering
        validator.add_experimental_result(ExperimentalResult("A", "Escherichia_coli", 300))
        validator.add_experimental_result(ExperimentalResult("B", "Escherichia_coli", 200))
        validator.add_experimental_result(ExperimentalResult("C", "Escherichia_coli", 100))

        for gene, cai in [("A", 0.70), ("B", 0.85), ("C", 0.95)]:
            opt = self._make_opt_result(cai=cai)
            validator.compare_with_prediction(opt, gene)

        rho, p = validator.compute_rank_correlation()
        assert rho == -1.0  # Perfect inverse rank correlation


# ═══════════════════════════════════════════════════════════════════════
# Validation Report Tests
# ═══════════════════════════════════════════════════════════════════════

class TestValidationReport:
    """Test validation report generation."""

    def _make_opt_result(self, cai: float = 0.95) -> OptimizationResult:
        return OptimizationResult(
            sequence="ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTGGTCC",
            gc_content=0.50,
            cai=cai,
        )

    def test_report_no_data(self):
        validator = WetLabValidator()
        report = validator.validation_report()
        assert report["status"] == "no_data"
        assert report["num_comparisons"] == 0

    def test_report_with_data(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        validator.add_experimental_result(ExperimentalResult("GFP", "Escherichia_coli", 5.0e5))
        validator.add_experimental_result(ExperimentalResult("EPO", "Homo_sapiens", 2.5e4))

        for gene, cai in [("INS", 0.95), ("GFP", 0.88), ("EPO", 0.72)]:
            opt = self._make_opt_result(cai=cai)
            validator.compare_with_prediction(opt, gene)

        report = validator.validation_report()
        assert report["status"] == "ok"
        assert report["num_results"] == 3
        assert report["num_comparisons"] == 3
        assert "pearson_correlation" in report
        assert "spearman_rank_correlation" in report
        assert isinstance(report["pearson_correlation"]["r"], float)
        assert isinstance(report["spearman_rank_correlation"]["rho"], float)

    def test_report_contains_comparison_details(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        opt = self._make_opt_result(cai=0.95)
        validator.compare_with_prediction(opt, "INS")

        report = validator.validation_report()
        assert len(report["comparisons"]) == 1
        assert report["comparisons"][0]["gene_name"] == "INS"
        assert report["comparisons"][0]["predicted_cai"] == 0.95

    def test_report_contains_stats(self):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))
        validator.add_experimental_result(ExperimentalResult("GFP", "Escherichia_coli", 5.0e5))
        opt1 = self._make_opt_result(cai=0.95)
        opt2 = self._make_opt_result(cai=0.88)
        validator.compare_with_prediction(opt1, "INS")
        validator.compare_with_prediction(opt2, "GFP")

        report = validator.validation_report()
        assert "predicted_cai_stats" in report
        assert "measured_expression_stats" in report
        assert report["predicted_cai_stats"]["min"] == 0.88
        assert report["predicted_cai_stats"]["max"] == 0.95


# ═══════════════════════════════════════════════════════════════════════
# CSV Loading Tests
# ═══════════════════════════════════════════════════════════════════════

class TestCSVLoading:
    """Test loading experimental results from CSV."""

    def test_load_from_csv(self, tmp_path):
        csv_content = (
            "gene_name,organism,measured_expression_level,measured_cai,notes\n"
            "INS,Escherichia_coli,1200000,0.88,Shake flask\n"
            "GFP,Escherichia_coli,500000,,Fluorescence\n"
            "EPO,Homo_sapiens,25000,0.72,HEK293T\n"
        )
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(csv_content)

        validator = WetLabValidator()
        count = validator.load_from_csv(str(csv_path))
        assert count == 3
        assert len(validator.results) == 3

    def test_load_from_csv_missing_columns(self, tmp_path):
        csv_content = "gene_name,organism\nINS,Escherichia_coli\n"
        csv_path = tmp_path / "bad.csv"
        csv_path.write_text(csv_content)

        validator = WetLabValidator()
        with pytest.raises(ValueError, match="Missing required columns"):
            validator.load_from_csv(str(csv_path))

    def test_load_from_csv_file_not_found(self):
        validator = WetLabValidator()
        with pytest.raises(FileNotFoundError):
            validator.load_from_csv("/nonexistent/path.csv")

    def test_load_from_csv_optional_fields(self, tmp_path):
        csv_content = (
            "gene_name,organism,measured_expression_level,measured_cai,sequence_used,notes\n"
            "INS,Escherichia_coli,1200000,0.88,ATGGCT,Test\n"
        )
        csv_path = tmp_path / "full.csv"
        csv_path.write_text(csv_content)

        validator = WetLabValidator()
        count = validator.load_from_csv(str(csv_path))
        assert count == 1
        result = validator.results[0]
        assert result.measured_cai == 0.88
        assert result.sequence_used == "ATGGCT"
        assert result.notes == "Test"

    def test_load_from_csv_skip_invalid_rows(self, tmp_path):
        csv_content = (
            "gene_name,organism,measured_expression_level\n"
            "INS,Escherichia_coli,1200000\n"
            "BAD,Escherichia_coli,not_a_number\n"
            "GFP,Escherichia_coli,500000\n"
        )
        csv_path = tmp_path / "mixed.csv"
        csv_path.write_text(csv_content)

        validator = WetLabValidator()
        count = validator.load_from_csv(str(csv_path))
        # BAD row should be skipped, but INS and GFP should load
        assert count == 2


# ═══════════════════════════════════════════════════════════════════════
# Report Saving Tests
# ═══════════════════════════════════════════════════════════════════════

class TestSaveReport:
    """Test saving validation reports as JSON."""

    def test_save_report(self, tmp_path):
        validator = WetLabValidator()
        validator.add_experimental_result(ExperimentalResult("INS", "Escherichia_coli", 1.2e6))

        opt = OptimizationResult(
            sequence="ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTGGTCC",
            gc_content=0.50,
            cai=0.95,
        )
        validator.compare_with_prediction(opt, "INS")

        output_path = tmp_path / "report.json"
        validator.save_report(str(output_path))

        assert output_path.exists()
        with open(output_path) as f:
            report = json.load(f)
        assert report["status"] == "ok"
        assert report["num_comparisons"] == 1

    def test_save_report_creates_directories(self, tmp_path):
        validator = WetLabValidator()
        output_path = tmp_path / "nested" / "dir" / "report.json"
        validator.save_report(str(output_path))
        assert output_path.exists()

    def test_save_report_no_data(self, tmp_path):
        validator = WetLabValidator()
        output_path = tmp_path / "empty_report.json"
        validator.save_report(str(output_path))

        with open(output_path) as f:
            report = json.load(f)
        assert report["status"] == "no_data"


# ═══════════════════════════════════════════════════════════════════════
# Integration / End-to-End Tests
# ═══════════════════════════════════════════════════════════════════════

class TestWetLabIntegration:
    """End-to-end integration tests."""

    def _make_opt_result(self, cai: float = 0.95) -> OptimizationResult:
        return OptimizationResult(
            sequence="ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTGGTCC",
            gc_content=0.50,
            cai=cai,
        )

    def test_full_pipeline(self, tmp_path):
        """Test the complete pipeline: add data → compare → report → save."""
        validator = WetLabValidator()

        # Step 1: Add experimental results
        validator.add_experimental_result(ExperimentalResult(
            "INS", "Escherichia_coli", 1.2e6, 0.90, notes="High expression",
        ))
        validator.add_experimental_result(ExperimentalResult(
            "GFP", "Escherichia_coli", 5.0e5, 0.85, notes="Moderate expression",
        ))
        validator.add_experimental_result(ExperimentalResult(
            "EPO", "Homo_sapiens", 2.5e4, 0.65, notes="Low expression",
        ))

        # Step 2: Compare with predictions
        for gene, cai in [("INS", 0.95), ("GFP", 0.88), ("EPO", 0.72)]:
            opt = self._make_opt_result(cai=cai)
            comp = validator.compare_with_prediction(opt, gene)
            assert comp.gene_name == gene
            assert comp.predicted_cai == cai

        # Step 3: Compute correlations
        r, p = validator.compute_correlation()
        assert r > 0  # Positive correlation expected

        rho, p_rho = validator.compute_rank_correlation()
        assert rho > 0  # Positive rank correlation expected

        # Step 4: Generate and save report
        report = validator.validation_report()
        assert report["status"] == "ok"
        assert report["num_comparisons"] == 3

        output_path = tmp_path / "validation_report.json"
        validator.save_report(str(output_path))
        assert output_path.exists()

    def test_csv_pipeline(self, tmp_path):
        """Test loading from CSV and comparing."""
        # Create CSV
        csv_content = (
            "gene_name,organism,measured_expression_level,measured_cai\n"
            "INS,Escherichia_coli,1200000,0.90\n"
            "GFP,Escherichia_coli,500000,0.85\n"
        )
        csv_path = tmp_path / "exp_results.csv"
        csv_path.write_text(csv_content)

        # Load
        validator = WetLabValidator()
        count = validator.load_from_csv(str(csv_path))
        assert count == 2

        # Compare
        for gene, cai in [("INS", 0.95), ("GFP", 0.88)]:
            opt = self._make_opt_result(cai=cai)
            comp = validator.compare_with_prediction(opt, gene)
            assert comp.predicted_cai == cai

        # Report
        report = validator.validation_report()
        assert report["num_comparisons"] == 2


# ═══════════════════════════════════════════════════════════════════════
# Global Validator Tests
# ═══════════════════════════════════════════════════════════════════════

class TestGlobalValidator:
    """Test the module-level global validator instance."""

    def test_global_validator_exists(self):
        from biocompiler.validation.wetlab_validation import _global_validator
        assert isinstance(_global_validator, WetLabValidator)

    def test_global_validator_is_independent(self):
        """Each import gets the same singleton."""
        from biocompiler.validation.wetlab_validation import _global_validator as v1
        from biocompiler.validation.wetlab_validation import _global_validator as v2
        assert v1 is v2
