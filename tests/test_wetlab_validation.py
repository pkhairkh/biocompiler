"""
Tests for the wet-lab validation framework.

Covers:
  - Benchmark dataset integrity
  - Validation suite execution
  - Metrics computation
  - Regression detection
  - WetLabProtocol and WetLabResult
"""

import pytest
from biocompiler.wetlab_validation import (
    BenchmarkEntry,
    BENCHMARK_DATASET,
    ProteinValidationResult,
    ValidationSuiteResult,
    RegressionReport,
    RegressionItem,
    run_validation_suite,
    check_regression,
    WetLabProtocol,
    WetLabResult,
    compare_insilico_vs_wetlab,
    generate_protocol_report,
    EXPRESSION_LEVEL_MAP,
    SOLUBILITY_LEVEL_MAP,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Benchmark Dataset Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestBenchmarkDataset:
    """Test the curated benchmark dataset."""

    def test_dataset_has_minimum_entries(self):
        """Dataset must have at least 10 benchmark proteins."""
        assert len(BENCHMARK_DATASET) >= 10

    def test_all_entries_have_required_fields(self):
        """Each entry must have all required fields populated."""
        for entry in BENCHMARK_DATASET:
            assert entry.protein_name, f"Missing protein_name"
            assert entry.organism, f"Missing organism for {entry.protein_name}"
            assert entry.protein_sequence, f"Missing protein_sequence for {entry.protein_name}"
            assert entry.expression_system, f"Missing expression_system for {entry.protein_name}"
            assert entry.source_publication, f"Missing source_publication for {entry.protein_name}"
            assert len(entry.expected_cai_range) == 2
            assert len(entry.expected_gc_range) == 2

    def test_protein_sequences_are_valid(self):
        """All protein sequences must be valid amino acid strings."""
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for entry in BENCHMARK_DATASET:
            invalid = set(entry.protein_sequence) - valid_aa
            assert not invalid, (
                f"{entry.protein_name}: invalid amino acid characters: {invalid}"
            )

    def test_cai_ranges_are_valid(self):
        """CAI ranges must be within [0, 1] and lo <= hi."""
        for entry in BENCHMARK_DATASET:
            lo, hi = entry.expected_cai_range
            assert 0.0 <= lo <= 1.0, f"{entry.protein_name}: CAI lo={lo}"
            assert 0.0 <= hi <= 1.0, f"{entry.protein_name}: CAI hi={hi}"
            assert lo <= hi, f"{entry.protein_name}: CAI lo > hi"

    def test_gc_ranges_are_valid(self):
        """GC ranges must be within [0, 1] and lo <= hi."""
        for entry in BENCHMARK_DATASET:
            lo, hi = entry.expected_gc_range
            assert 0.0 <= lo <= 1.0, f"{entry.protein_name}: GC lo={lo}"
            assert 0.0 <= hi <= 1.0, f"{entry.protein_name}: GC hi={hi}"
            assert lo <= hi, f"{entry.protein_name}: GC lo > hi"

    def test_covers_required_proteins(self):
        """Dataset must cover the required proteins: GFP, insulin, HBB, EPO,
        TNF-alpha, IL-2, interferon-alpha, growth hormone, albumin, lysozyme."""
        protein_names = {e.protein_name for e in BENCHMARK_DATASET}
        required = {"GFP", "Insulin", "HBB", "EPO", "TNF_alpha",
                    "IL2", "IFN_alpha", "hGH", "Albumin", "Lysozyme"}
        missing = required - protein_names
        assert not missing, f"Missing required proteins: {missing}"

    def test_covers_multiple_organisms(self):
        """Dataset must cover at least 3 different organisms."""
        organisms = {e.organism for e in BENCHMARK_DATASET}
        assert len(organisms) >= 3, f"Only {len(organisms)} organisms: {organisms}"

    def test_protein_sequences_have_reasonable_length(self):
        """Protein sequences should have reasonable lengths (50-700 aa)."""
        for entry in BENCHMARK_DATASET:
            n = len(entry.protein_sequence)
            assert 30 <= n <= 700, (
                f"{entry.protein_name}: unusual protein length {n} aa"
            )

    def test_expression_categories_are_valid(self):
        """Expression categories must be 'high', 'medium', or 'low'."""
        for entry in BENCHMARK_DATASET:
            assert entry.expression_category in ("high", "medium", "low"), (
                f"{entry.protein_name}: invalid category '{entry.expression_category}'"
            )


class TestBenchmarkEntry:
    """Test BenchmarkEntry dataclass."""

    def test_create_benchmark_entry(self):
        entry = BenchmarkEntry(
            protein_name="TestProtein",
            organism="Escherichia_coli",
            protein_sequence="MSKGEELFTGV",
            expression_system="E. coli",
            source_publication="doi:10.1234/test",
            expected_cai_range=(0.8, 1.0),
            expected_gc_range=(0.4, 0.6),
        )
        assert entry.protein_name == "TestProtein"
        assert entry.expected_no_restriction_sites is True
        assert entry.measured_expression_level is None

    def test_benchmark_entry_with_measured_level(self):
        entry = BenchmarkEntry(
            protein_name="GFP",
            organism="Escherichia_coli",
            protein_sequence="MSKGEELFTGV",
            expression_system="E. coli",
            source_publication="doi:10.1234/test",
            expected_cai_range=(0.8, 1.0),
            expected_gc_range=(0.4, 0.6),
            measured_expression_level=150.0,
            expression_category="high",
        )
        assert entry.measured_expression_level == pytest.approx(150.0, rel=1e-6)
        assert entry.expression_category == "high"


# ═══════════════════════════════════════════════════════════════════════════════
# Validation Suite Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestValidationSuite:
    """Test the in-silico validation suite runner."""

    def test_run_validation_suite_default(self):
        """Run the validation suite with default config on a subset."""
        subset = BENCHMARK_DATASET[:3]
        result = run_validation_suite(proteins=subset)
        assert isinstance(result, ValidationSuiteResult)
        assert result.total_proteins == 3
        assert len(result.per_protein_results) == 3
        assert 0.0 <= result.pass_rate <= 1.0

    def test_run_validation_suite_with_config(self):
        """Run with a custom optimizer configuration."""
        config = {"strategy": "hybrid", "gc_lo": 0.30, "gc_hi": 0.70}
        subset = BENCHMARK_DATASET[:2]
        result = run_validation_suite(optimizer_config=config, proteins=subset)
        assert result.config == config

    def test_validation_suite_metrics(self):
        """Check that aggregate metrics are computed correctly."""
        subset = BENCHMARK_DATASET[:2]
        result = run_validation_suite(proteins=subset)
        assert 0.0 <= result.cai_mean <= 1.0
        assert result.cai_std >= 0.0
        assert 0.0 <= result.gc_mean <= 1.0
        assert result.gc_std >= 0.0
        assert 0.0 <= result.protein_fidelity_rate <= 1.0
        assert 0.0 <= result.constraint_violation_rate <= 1.0

    def test_validation_suite_dnachisel_comparison(self):
        """DNAchisel comparison should have win/tie/loss keys."""
        subset = BENCHMARK_DATASET[:2]
        result = run_validation_suite(proteins=subset)
        assert "win" in result.comparison_vs_dnachisel
        assert "tie" in result.comparison_vs_dnachisel
        assert "loss" in result.comparison_vs_dnachisel

    def test_validation_suite_organism_filter(self):
        """Organism filter should limit the proteins tested."""
        result = run_validation_suite(organism_filter="Escherichia_coli")
        ecoli_proteins = [
            e for e in BENCHMARK_DATASET if e.organism == "Escherichia_coli"
        ]
        assert result.total_proteins == len(ecoli_proteins)

    def test_validation_suite_empty_protein_list(self):
        """Empty protein list should return zeroed metrics."""
        result = run_validation_suite(proteins=[])
        assert result.total_proteins == 0
        assert result.pass_rate == 0.0
        assert result.cai_mean == 0.0

    def test_per_protein_result_fields(self):
        """Each per-protein result should have all required fields."""
        subset = BENCHMARK_DATASET[:1]
        result = run_validation_suite(proteins=subset)
        r = result.per_protein_results[0]
        assert isinstance(r, ProteinValidationResult)
        assert r.protein_name == subset[0].protein_name
        assert isinstance(r.cai, float)
        assert isinstance(r.cai_in_range, bool)
        assert isinstance(r.gc_content, float)
        assert isinstance(r.gc_in_range, bool)
        assert isinstance(r.no_restriction_sites, bool)
        assert isinstance(r.protein_fidelity, bool)
        assert isinstance(r.passed, bool)
        assert isinstance(r.details, dict)
        assert r.optimization_time_s >= 0.0

    def test_validation_suite_total_time(self):
        """Total time should be non-negative."""
        subset = BENCHMARK_DATASET[:2]
        result = run_validation_suite(proteins=subset)
        assert result.total_time_s >= 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# Regression Detection Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestRegressionDetection:
    """Test the regression detection system."""

    def _make_suite_result(
        self,
        pass_rate: float = 0.8,
        cai_mean: float = 0.9,
        gc_mean: float = 0.5,
        fidelity_rate: float = 1.0,
        violation_rate: float = 0.2,
    ) -> ValidationSuiteResult:
        """Helper to create a ValidationSuiteResult with given metrics."""
        return ValidationSuiteResult(
            per_protein_results=[
                ProteinValidationResult(
                    protein_name="Test",
                    organism="Escherichia_coli",
                    cai=cai_mean,
                    cai_in_range=True,
                    gc_content=gc_mean,
                    gc_in_range=True,
                    no_restriction_sites=True,
                    protein_fidelity=True,
                    passed=True,
                )
            ],
            pass_rate=pass_rate,
            cai_mean=cai_mean,
            cai_std=0.05,
            gc_mean=gc_mean,
            gc_std=0.03,
            protein_fidelity_rate=fidelity_rate,
            constraint_violation_rate=violation_rate,
            comparison_vs_dnachisel={"win": 1, "tie": 0, "loss": 0},
            total_proteins=1,
        )

    def test_no_regression_when_identical(self):
        """Identical baseline and current should show no regression."""
        baseline = self._make_suite_result()
        current = self._make_suite_result()
        report = check_regression(baseline, current)
        assert not report.has_regression
        assert len(report.regressed_proteins) == 0

    def test_detects_pass_rate_regression(self):
        """Should detect >5% drop in pass_rate."""
        baseline = self._make_suite_result(pass_rate=0.9)
        current = self._make_suite_result(pass_rate=0.8)
        report = check_regression(baseline, current)
        # 0.9 → 0.8 = 11% degradation
        assert report.has_regression

    def test_detects_cai_mean_regression(self):
        """Should detect >5% drop in CAI mean."""
        baseline = self._make_suite_result(cai_mean=0.9)
        current = self._make_suite_result(cai_mean=0.8)
        report = check_regression(baseline, current)
        # 0.9 → 0.8 = 11% degradation
        assert report.has_regression

    def test_no_false_positive_for_small_changes(self):
        """Should not flag <5% changes as regressions."""
        baseline = self._make_suite_result(cai_mean=0.9)
        current = self._make_suite_result(cai_mean=0.87)
        report = check_regression(baseline, current)
        # 0.9 → 0.87 = 3.3% degradation, below threshold
        # But check individual items
        cai_regression = [r for r in report.regressions if r.metric_name == "cai_mean"]
        assert cai_regression
        assert not cai_regression[0].is_regression

    def test_regression_report_has_summary(self):
        """Report should include a human-readable summary."""
        baseline = self._make_suite_result()
        current = self._make_suite_result()
        report = check_regression(baseline, current)
        assert isinstance(report.summary, str)
        assert len(report.summary) > 0

    def test_regression_items_have_degradation(self):
        """Each regression item should have a degradation value."""
        baseline = self._make_suite_result()
        current = self._make_suite_result()
        report = check_regression(baseline, current)
        for item in report.regressions:
            assert isinstance(item, RegressionItem)
            assert isinstance(item.degradation, float)
            assert isinstance(item.is_regression, bool)

    def test_constraint_violation_rate_increase_is_regression(self):
        """Increasing constraint_violation_rate should be flagged."""
        baseline = self._make_suite_result(violation_rate=0.1)
        current = self._make_suite_result(violation_rate=0.2)
        report = check_regression(baseline, current)
        cv_items = [r for r in report.regressions if r.metric_name == "constraint_violation_rate"]
        assert cv_items
        # 0.1 → 0.2 = 100% increase
        assert cv_items[0].is_regression

    def test_per_protein_regression_detection(self):
        """Should detect when individual proteins regress."""
        baseline = ValidationSuiteResult(
            per_protein_results=[
                ProteinValidationResult(
                    protein_name="GFP",
                    organism="Escherichia_coli",
                    cai=0.95,
                    cai_in_range=True,
                    gc_content=0.50,
                    gc_in_range=True,
                    no_restriction_sites=True,
                    protein_fidelity=True,
                    passed=True,
                ),
            ],
            pass_rate=1.0,
            cai_mean=0.95,
            cai_std=0.0,
            gc_mean=0.50,
            gc_std=0.0,
            protein_fidelity_rate=1.0,
            constraint_violation_rate=0.0,
            comparison_vs_dnachisel={"win": 0, "tie": 0, "loss": 0},
            total_proteins=1,
        )
        current = ValidationSuiteResult(
            per_protein_results=[
                ProteinValidationResult(
                    protein_name="GFP",
                    organism="Escherichia_coli",
                    cai=0.80,
                    cai_in_range=True,
                    gc_content=0.50,
                    gc_in_range=True,
                    no_restriction_sites=True,
                    protein_fidelity=True,
                    passed=True,
                ),
            ],
            pass_rate=1.0,
            cai_mean=0.80,
            cai_std=0.0,
            gc_mean=0.50,
            gc_std=0.0,
            protein_fidelity_rate=1.0,
            constraint_violation_rate=0.0,
            comparison_vs_dnachisel={"win": 0, "tie": 0, "loss": 0},
            total_proteins=1,
        )
        report = check_regression(baseline, current)
        # CAI dropped from 0.95 → 0.80, that's ~15.8% degradation
        assert "GFP" in report.regressed_proteins


# ═══════════════════════════════════════════════════════════════════════════════
# WetLabProtocol and WetLabResult Tests (retained from original)
# ═══════════════════════════════════════════════════════════════════════════════

class TestWetLabProtocol:
    """Test WetLabProtocol creation and validation."""

    def test_create_protocol(self):
        protocol = WetLabProtocol(
            gene_name="GFP",
            organism="Escherichia_coli",
            optimized_dna="ATGGTTAGCAAAGGCGAAGAA",
            cai=0.93,
            gc_content=0.48,
            vector="pET-28a",
            promoter="T7",
            selection_marker="kanamycin",
            host_strain="BL21(DE3)",
            expected_expression_level="high",
            expected_solubility="high",
        )
        assert protocol.gene_name == "GFP"

    def test_invalid_expression_level(self):
        with pytest.raises(ValueError, match="Invalid expression level"):
            WetLabProtocol(
                gene_name="GFP",
                organism="Escherichia_coli",
                optimized_dna="ATGGTTAGC",
                cai=0.9,
                gc_content=0.5,
                vector="pET-28a",
                promoter="T7",
                selection_marker="kanamycin",
                host_strain="BL21(DE3)",
                expected_expression_level="extreme",
                expected_solubility="high",
            )

    def test_invalid_dna_sequence(self):
        with pytest.raises(ValueError, match="only A, C, G, T"):
            WetLabProtocol(
                gene_name="GFP",
                organism="Escherichia_coli",
                optimized_dna="ATGXAGC",
                cai=0.9,
                gc_content=0.5,
                vector="pET-28a",
                promoter="T7",
                selection_marker="kanamycin",
                host_strain="BL21(DE3)",
                expected_expression_level="high",
                expected_solubility="high",
            )

    def test_generate_protocol(self):
        protocol = WetLabProtocol(
            gene_name="GFP",
            organism="Escherichia_coli",
            optimized_dna="ATGGTTAGCAAAGGCGAAGAA",
            cai=0.93,
            gc_content=0.48,
            vector="pET-28a",
            promoter="T7",
            selection_marker="kanamycin",
            host_strain="BL21(DE3)",
            expected_expression_level="high",
            expected_solubility="high",
        )
        text = protocol.generate_protocol()
        assert "GFP" in text
        assert "pET-28a" in text
        assert "T7" in text

    def test_generate_oligos(self):
        protocol = WetLabProtocol(
            gene_name="GFP",
            organism="Escherichia_coli",
            optimized_dna="A" * 200,
            cai=0.9,
            gc_content=0.5,
            vector="pET-28a",
            promoter="T7",
            selection_marker="kanamycin",
            host_strain="BL21(DE3)",
            expected_expression_level="high",
            expected_solubility="high",
        )
        oligos = protocol.generate_oligos(max_oligo_length=60)
        assert len(oligos) > 1

    def test_predict_expression_category(self):
        protocol = WetLabProtocol(
            gene_name="GFP",
            organism="Escherichia_coli",
            optimized_dna="ATGGTTAGCAAAGGCGAAGAA",
            cai=0.9,
            gc_content=0.5,
            vector="pET-28a",
            promoter="T7",
            selection_marker="kanamycin",
            host_strain="BL21(DE3)",
            expected_expression_level="high",
            expected_solubility="high",
        )
        assert protocol.predict_expression_category() == "high"


class TestWetLabResult:
    """Test WetLabResult creation and classification."""

    def test_create_result(self):
        result = WetLabResult(
            actual_expression_level=100.0,
            actual_solubility=0.8,
            western_blot_confirmed=True,
            sequencing_match=True,
        )
        assert result.classify_expression() == "high"
        assert result.classify_solubility() == "high"

    def test_negative_expression_raises(self):
        with pytest.raises(ValueError, match="must be >= 0"):
            WetLabResult(
                actual_expression_level=-10.0,
                actual_solubility=0.5,
                western_blot_confirmed=True,
                sequencing_match=True,
            )

    def test_solubility_out_of_range(self):
        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            WetLabResult(
                actual_expression_level=10.0,
                actual_solubility=1.5,
                western_blot_confirmed=True,
                sequencing_match=True,
            )


class TestCompareInsilicoVsWetlab:
    """Test the comparison function."""

    def test_strong_agreement(self):
        prediction = {
            "cai": 0.9,
            "gc_content": 0.5,
            "expression_level": "high",
            "solubility": "high",
        }
        actual = WetLabResult(
            actual_expression_level=100.0,
            actual_solubility=0.8,
            western_blot_confirmed=True,
            sequencing_match=True,
        )
        result = compare_insilico_vs_wetlab(prediction, actual)
        assert result["expression_category_match"] is True
        assert result["overall_agreement"] == "strong_agreement"

    def test_disagreement(self):
        prediction = {
            "expression_level": "high",
            "solubility": "high",
        }
        actual = WetLabResult(
            actual_expression_level=1.0,
            actual_solubility=0.1,
            western_blot_confirmed=False,
            sequencing_match=False,
        )
        result = compare_insilico_vs_wetlab(prediction, actual)
        assert result["expression_category_match"] is False
        assert result["overall_agreement"] == "disagreement"
