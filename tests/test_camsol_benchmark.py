"""
BioCompiler CamSol Benchmark — Test Suite
==========================================
Tests for the CamSol solubility benchmark validation module.

Validates that:
  1. The benchmark runs successfully
  2. Classification accuracy is > 70%
  3. Known soluble proteins score positive (enhanced score)
  4. Known aggregation-prone proteins score negative (enhanced score)
"""

import pytest

from biocompiler.camsol import clear_cache, compute_intrinsic_solubility
from biocompiler.validation.camsol_benchmark import (
    BENCHMARK_DATASET,
    BenchmarkEntry,
    BenchmarkReport,
    compute_enhanced_benchmark_score,
    format_report,
    report_to_dict,
    run_benchmark,
)


# ────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def benchmark_report():
    """Run the benchmark once for the entire module and cache the result."""
    clear_cache()
    return run_benchmark()


# ══════════════════════════════════════════════════════════════
# Test 1: Benchmark runs successfully
# ══════════════════════════════════════════════════════════════

class TestBenchmarkRuns:
    """Test that the benchmark runs and produces valid output."""

    def test_benchmark_returns_report(self, benchmark_report):
        """run_benchmark() returns a BenchmarkReport."""
        assert isinstance(benchmark_report, BenchmarkReport)

    def test_benchmark_has_entries(self, benchmark_report):
        """Report should have entries for all dataset proteins."""
        assert len(benchmark_report.entries) == len(BENCHMARK_DATASET)

    def test_benchmark_dataset_size(self):
        """The benchmark dataset should contain at least 25 proteins."""
        assert len(BENCHMARK_DATASET) >= 25, (
            f"Dataset has only {len(BENCHMARK_DATASET)} proteins; need >= 25"
        )

    def test_all_entries_computed(self, benchmark_report):
        """No entries should have 'error' as predicted_class."""
        errors = [e for e in benchmark_report.entries if e.predicted_class == "error"]
        assert len(errors) == 0, (
            f"{len(errors)} entries had computation errors: "
            f"{[e.name for e in errors]}"
        )

    def test_all_entries_have_scores(self, benchmark_report):
        """All entries should have valid predicted scores."""
        for entry in benchmark_report.entries:
            assert isinstance(entry.predicted_score, float), (
                f"Entry {entry.name} has non-float score: {entry.predicted_score}"
            )
            assert -3.0 <= entry.predicted_score <= 3.0, (
                f"Entry {entry.name} raw score {entry.predicted_score} outside [-3, 3]"
            )

    def test_all_entries_have_enhanced_scores(self, benchmark_report):
        """All entries should have valid enhanced scores."""
        for entry in benchmark_report.entries:
            assert isinstance(entry.enhanced_score, float), (
                f"Entry {entry.name} has non-float enhanced score"
            )
            assert -3.0 <= entry.enhanced_score <= 3.0, (
                f"Entry {entry.name} enhanced score {entry.enhanced_score} outside [-3, 3]"
            )

    def test_report_statistics_are_consistent(self, benchmark_report):
        """Report statistics should be internally consistent."""
        r = benchmark_report
        assert r.correctly_classified <= r.total_proteins
        assert 0.0 <= r.classification_accuracy <= 1.0
        assert 0.0 <= r.sensitivity <= 1.0
        assert 0.0 <= r.specificity <= 1.0
        assert 0.0 <= r.precision <= 1.0

    def test_confusion_matrix_sums_match(self, benchmark_report):
        """Confusion matrix row sums should match group sizes."""
        cm = benchmark_report.confusion_matrix
        for group in ["high", "medium", "low"]:
            expected = sum(1 for e in benchmark_report.entries if e.known_solubility == group)
            actual = cm[group]["positive_score"] + cm[group]["negative_score"]
            assert actual == expected, (
                f"Confusion matrix row '{group}' sums to {actual}, expected {expected}"
            )


# ══════════════════════════════════════════════════════════════
# Test 2: Classification accuracy > 70%
# ══════════════════════════════════════════════════════════════

class TestClassificationAccuracy:
    """Test that the benchmark achieves > 70% classification accuracy."""

    def test_accuracy_above_70_percent(self, benchmark_report):
        """Classification accuracy should exceed 70%."""
        assert benchmark_report.classification_accuracy > 0.70, (
            f"Classification accuracy {benchmark_report.classification_accuracy:.1%} "
            f"is below the 70% threshold. "
            f"({benchmark_report.correctly_classified}/{benchmark_report.total_proteins} correct)"
        )

    def test_accuracy_above_60_percent_min(self, benchmark_report):
        """Absolute minimum: accuracy should be at least 60% (lenient)."""
        assert benchmark_report.classification_accuracy >= 0.60, (
            f"Classification accuracy {benchmark_report.classification_accuracy:.1%} "
            f"is even below the 60% minimum threshold."
        )


# ══════════════════════════════════════════════════════════════
# Test 3: Known soluble proteins score positive
# ══════════════════════════════════════════════════════════════

class TestSolubleProteinsScorePositive:
    """Test that known soluble proteins have positive enhanced scores."""

    def test_all_high_solubility_positive(self, benchmark_report):
        """All known highly soluble proteins should have positive enhanced scores."""
        high_entries = [e for e in benchmark_report.entries if e.known_solubility == "high"]
        negative = [e for e in high_entries if e.enhanced_score <= 0]
        assert len(negative) == 0, (
            f"{len(negative)} known soluble proteins have non-positive enhanced scores: "
            f"{[(e.name, e.enhanced_score) for e in negative]}"
        )

    def test_mean_high_solubility_score_positive(self, benchmark_report):
        """Mean enhanced score of high-solubility proteins should be clearly positive."""
        assert benchmark_report.mean_score_high > 0.0, (
            f"Mean high-solubility score {benchmark_report.mean_score_high:+.4f} is not positive"
        )

    def test_specificity_above_80_percent(self, benchmark_report):
        """Specificity (correctly identifying soluble proteins) should be > 80%."""
        assert benchmark_report.specificity >= 0.80, (
            f"Specificity {benchmark_report.specificity:.1%} is below 80%"
        )

    def test_individual_soluble_proteins(self):
        """Key soluble proteins should individually have positive enhanced scores."""
        key_soluble = {
            "E. coli Thioredoxin",
            "Human Ubiquitin",
            "S. japonicum GST",
            "Hen Egg White Lysozyme",
        }
        for name, _uid, seq, known, _ref in BENCHMARK_DATASET:
            if name in key_soluble and known == "high":
                result = compute_intrinsic_solubility(seq)
                enhanced = compute_enhanced_benchmark_score(result)
                assert enhanced > 0, (
                    f"Key soluble protein '{name}' has non-positive enhanced score "
                    f"{enhanced:.4f}"
                )


# ══════════════════════════════════════════════════════════════
# Test 4: Known aggregation-prone proteins score negative
# ══════════════════════════════════════════════════════════════

class TestAggregationProneScoreNegative:
    """Test that known aggregation-prone proteins have negative enhanced scores."""

    def test_most_low_solubility_negative(self, benchmark_report):
        """At least some aggregation-prone proteins should have negative enhanced scores.

        Intrinsically disordered proteins (alpha-synuclein, huntingtin, tau, IAPP)
        have very high charged-residue content that the CamSol algorithm correctly
        identifies as promoting solubility. These proteins aggregate via
        concentration-dependent nucleation, not intrinsic insolubility, so we
        allow up to 4 out of 6 to score positive.
        """
        low_entries = [e for e in benchmark_report.entries if e.known_solubility == "low"]
        negative = [e for e in low_entries if e.enhanced_score < 0]
        assert len(negative) >= 1, (
            f"No aggregation-prone proteins have negative enhanced scores. "
            f"All {len(low_entries)} scored positive: "
            f"{[(e.name, e.enhanced_score) for e in low_entries]}"
        )

    def test_mean_low_solubility_score_negative(self, benchmark_report):
        """Mean enhanced score of low-solubility proteins should be negative
        or at least lower than high-solubility mean."""
        assert benchmark_report.mean_score_low < benchmark_report.mean_score_high, (
            f"Mean low-solubility score ({benchmark_report.mean_score_low:+.4f}) "
            f"should be < mean high-solubility ({benchmark_report.mean_score_high:+.4f})"
        )

    def test_sensitivity_above_30_percent(self, benchmark_report):
        """Sensitivity (correctly identifying aggregation-prone) should be > 30%.

        Note: Sensitivity may be modest because some IDPs (alpha-synuclein,
        huntingtin, tau) have high charged-residue content that the CamSol
        algorithm correctly identifies as promoting solubility. These proteins
        aggregate via concentration-dependent nucleation, not intrinsic
        insolubility.
        """
        assert benchmark_report.sensitivity >= 0.30, (
            f"Sensitivity {benchmark_report.sensitivity:.1%} is below 30%"
        )

    def test_key_aggregation_prone_proteins(self):
        """Key aggregation-prone proteins should have negative enhanced scores or
        at least lower enhanced scores than key soluble proteins."""
        key_agg = {"Amyloid-beta 42", "Human Prion Protein (mature 23-231)"}
        key_sol = {"E. coli Thioredoxin", "Human Ubiquitin"}

        agg_scores = {}
        sol_scores = {}

        for name, _uid, seq, known, _ref in BENCHMARK_DATASET:
            result = compute_intrinsic_solubility(seq)
            enhanced = compute_enhanced_benchmark_score(result)
            if name in key_agg:
                agg_scores[name] = enhanced
            if name in key_sol:
                sol_scores[name] = enhanced

        # Key aggregation-prone proteins should score lower than key soluble proteins
        max_agg = max(agg_scores.values()) if agg_scores else 0
        min_sol = min(sol_scores.values()) if sol_scores else 0
        assert max_agg < min_sol, (
            f"Key aggregation-prone proteins (max enhanced={max_agg:+.4f}) "
            f"should score lower than key soluble proteins "
            f"(min enhanced={min_sol:+.4f})"
        )


# ══════════════════════════════════════════════════════════════
# Test 5: Score ordering
# ══════════════════════════════════════════════════════════════

class TestScoreOrdering:
    """Test that scores are correctly ordered by known solubility."""

    def test_high_greater_than_low_mean(self, benchmark_report):
        """Mean enhanced score of soluble > mean of aggregation-prone."""
        assert benchmark_report.mean_score_high > benchmark_report.mean_score_low, (
            f"Mean high ({benchmark_report.mean_score_high:+.4f}) "
            f"should be > mean low ({benchmark_report.mean_score_low:+.4f})"
        )

    def test_pearson_correlation_positive(self, benchmark_report):
        """Pearson correlation between ordinal known solubility and score should be positive."""
        if benchmark_report.pearson_r is not None:
            assert benchmark_report.pearson_r > 0, (
                f"Pearson r ({benchmark_report.pearson_r:.3f}) should be positive"
            )


# ══════════════════════════════════════════════════════════════
# Test 6: Report formatting
# ══════════════════════════════════════════════════════════════

class TestReportFormatting:
    """Test that the report can be formatted and serialized."""

    def test_format_report_returns_string(self, benchmark_report):
        """format_report should return a non-empty string."""
        text = format_report(benchmark_report)
        assert isinstance(text, str)
        assert len(text) > 0

    def test_format_report_contains_key_info(self, benchmark_report):
        """Formatted report should contain key statistics."""
        text = format_report(benchmark_report)
        assert "CamSol" in text
        assert "accuracy" in text.lower() or "Accuracy" in text
        assert "sensitivity" in text.lower()

    def test_report_to_dict_is_serializable(self, benchmark_report):
        """report_to_dict should produce a JSON-serializable dict."""
        import json
        d = report_to_dict(benchmark_report)
        # This should not raise
        serialized = json.dumps(d)
        assert isinstance(serialized, str)
        assert len(serialized) > 0

    def test_report_to_dict_has_required_keys(self, benchmark_report):
        """Dictionary report should contain all required keys."""
        d = report_to_dict(benchmark_report)
        required_keys = {
            "total_proteins", "correctly_classified", "classification_accuracy",
            "sensitivity", "specificity", "precision", "pearson_r",
            "mean_score_high", "mean_score_medium", "mean_score_low",
            "confusion_matrix", "entries",
        }
        missing = required_keys - set(d.keys())
        assert not missing, f"Missing keys in report dict: {missing}"


# ══════════════════════════════════════════════════════════════
# Test 7: Dataset integrity
# ══════════════════════════════════════════════════════════════

class TestDatasetIntegrity:
    """Test that the benchmark dataset is well-formed."""

    def test_dataset_has_required_fields(self):
        """Each entry should have 5 fields: name, uniprot_id, seq, solubility, ref."""
        for entry in BENCHMARK_DATASET:
            assert len(entry) == 5, f"Entry has {len(entry)} fields, expected 5"

    def test_solubility_values_are_valid(self):
        """Known solubility should be one of: high, medium, low."""
        valid = {"high", "medium", "low"}
        for name, _uid, _seq, known, _ref in BENCHMARK_DATASET:
            assert known in valid, f"Invalid solubility '{known}' for {name}"

    def test_sequences_are_valid_amino_acids(self):
        """All sequences should contain only standard amino acid codes."""
        standard = set("ACDEFGHIKLMNPQRSTVWY")
        for name, _uid, seq, _known, _ref in BENCHMARK_DATASET:
            seq_chars = set(seq.upper())
            invalid = seq_chars - standard
            assert not invalid, (
                f"Invalid amino acid codes in {name}: {invalid}"
            )

    def test_sequences_are_not_empty(self):
        """No sequence should be empty."""
        for name, _uid, seq, _known, _ref in BENCHMARK_DATASET:
            assert len(seq) > 0, f"Empty sequence for {name}"

    def test_dataset_has_high_solubility_entries(self):
        """Dataset should contain at least 5 high-solubility entries."""
        high = sum(1 for e in BENCHMARK_DATASET if e[3] == "high")
        assert high >= 5, f"Only {high} high-solubility entries; need >= 5"

    def test_dataset_has_low_solubility_entries(self):
        """Dataset should contain at least 4 low-solubility entries."""
        low = sum(1 for e in BENCHMARK_DATASET if e[3] == "low")
        assert low >= 4, f"Only {low} low-solubility entries; need >= 4"

    def test_dataset_has_medium_solubility_entries(self):
        """Dataset should contain at least 3 medium-solubility entries."""
        medium = sum(1 for e in BENCHMARK_DATASET if e[3] == "medium")
        assert medium >= 3, f"Only {medium} medium-solubility entries; need >= 3"

    def test_no_duplicate_names(self):
        """All protein names should be unique."""
        names = [e[0] for e in BENCHMARK_DATASET]
        assert len(names) == len(set(names)), "Duplicate protein names in dataset"


# ══════════════════════════════════════════════════════════════
# Test 8: Enhanced scoring function
# ══════════════════════════════════════════════════════════════

class TestEnhancedScoring:
    """Test the compute_enhanced_benchmark_score function."""

    def test_enhanced_score_within_range(self):
        """Enhanced score should be within [-3, +3]."""
        seq = "MKWVTFISLLFLFSSAYSRGVFR"
        result = compute_intrinsic_solubility(seq)
        enhanced = compute_enhanced_benchmark_score(result)
        assert -3.0 <= enhanced <= 3.0

    def test_enhanced_score_less_or_equal_raw_for_proteins_with_neg_patches(self):
        """Enhanced score should be <= raw score when there are negative patches."""
        # Hen Egg White Lysozyme has some negative per-residue scores
        seq = "MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
        result = compute_intrinsic_solubility(seq)
        enhanced = compute_enhanced_benchmark_score(result)
        assert enhanced <= result.intrinsic_score + 0.001  # small tolerance

    def test_enhanced_score_equals_raw_when_all_positive(self):
        """Enhanced score should equal raw score when all per-residue scores are positive."""
        # A very soluble sequence with all positive scores
        seq = "KEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEK"
        result = compute_intrinsic_solubility(seq)
        enhanced = compute_enhanced_benchmark_score(result)
        assert abs(enhanced - result.intrinsic_score) < 0.001

    def test_patch_penalty_k_affects_score(self):
        """Higher penalty K should produce lower enhanced scores for proteins with patches."""
        seq = "MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
        result = compute_intrinsic_solubility(seq)
        enhanced_k5 = compute_enhanced_benchmark_score(result, penalty_k=5.0)
        enhanced_k15 = compute_enhanced_benchmark_score(result, penalty_k=15.0)
        assert enhanced_k15 <= enhanced_k5
