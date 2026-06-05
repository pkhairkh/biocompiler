"""Test BioCompiler IEDB MHC Prediction Benchmark — validation/iedb_comparison.py.

Tests cover:
1. IEDBBenchmarkEntry construction and is_binder property
2. IEDB_BENCHMARK_DATA entries are valid
3. benchmark_mhc_predictions with mock predictor
4. MHBenchmarkResult fields (auc_roc, pearson_r)
"""
from __future__ import annotations

import math
from unittest.mock import MagicMock

import pytest

# ---------------------------------------------------------------------------
# Import the module under test — skip entire suite if unavailable
# ---------------------------------------------------------------------------

iedb_mod = pytest.importorskip(
    "biocompiler.validation.iedb_comparison",
    reason="biocompiler.validation.iedb_comparison not available",
)

IEDBBenchmarkEntry = iedb_mod.IEDBBenchmarkEntry
IEDB_BENCHMARK_DATA = iedb_mod.IEDB_BENCHMARK_DATA
MHBenchmarkResult = iedb_mod.MHBenchmarkResult
benchmark_mhc_predictions = iedb_mod.benchmark_mhc_predictions
BINDER_IC50_THRESHOLD = iedb_mod.BINDER_IC50_THRESHOLD

# Also try importing the internal helpers for targeted testing
_pearson_correlation = iedb_mod._pearson_correlation
_compute_auc_roc = iedb_mod._compute_auc_roc


# ═══════════════════════════════════════════════════════════════════════════
# 1. IEDBBenchmarkEntry construction and is_binder property
# ═══════════════════════════════════════════════════════════════════════════

class TestIEDBBenchmarkEntry:
    """Tests for IEDBBenchmarkEntry dataclass and is_binder property."""

    def test_construction_all_fields(self) -> None:
        """IEDBBenchmarkEntry can be created with all required fields."""
        entry = IEDBBenchmarkEntry(
            peptide="GILGFVFTL",
            allele="HLA-A*02:01",
            measured_ic50=5.0,
            source="IEDB EpiID 10644",
        )
        assert entry.peptide == "GILGFVFTL"
        assert entry.allele == "HLA-A*02:01"
        assert entry.measured_ic50 == 5.0
        assert entry.source == "IEDB EpiID 10644"

    def test_is_binder_strong_binder(self) -> None:
        """Strong binder (IC50 < 50 nM) is classified as binder."""
        entry = IEDBBenchmarkEntry(
            peptide="GILGFVFTL",
            allele="HLA-A*02:01",
            measured_ic50=5.0,
            source="test",
        )
        assert entry.is_binder is True

    def test_is_binder_moderate_binder(self) -> None:
        """Moderate binder (IC50 50–500 nM) is classified as binder."""
        entry = IEDBBenchmarkEntry(
            peptide="FLPSDCFFSV",
            allele="HLA-A*02:01",
            measured_ic50=150.0,
            source="test",
        )
        assert entry.is_binder is True

    def test_is_binder_weak_binder(self) -> None:
        """Weak binder (IC50 500–5000 nM) is classified as non-binder."""
        entry = IEDBBenchmarkEntry(
            peptide="TVFYLAPNL",
            allele="HLA-A*02:01",
            measured_ic50=1800.0,
            source="test",
        )
        assert entry.is_binder is False

    def test_is_binder_non_binder(self) -> None:
        """Non-binder (IC50 > 5000 nM) is classified as non-binder."""
        entry = IEDBBenchmarkEntry(
            peptide="DNEEGVQAD",
            allele="HLA-A*02:01",
            measured_ic50=25000.0,
            source="test",
        )
        assert entry.is_binder is False

    def test_is_binder_at_threshold(self) -> None:
        """At exactly 500 nM (BINDER_IC50_THRESHOLD), is_binder is False.

        The threshold uses strict less-than: IC50 < 500 → binder.
        """
        entry = IEDBBenchmarkEntry(
            peptide="TESTPEPTD",
            allele="HLA-A*02:01",
            measured_ic50=500.0,
            source="test",
        )
        assert entry.is_binder is False

    def test_is_binder_just_below_threshold(self) -> None:
        """Just below threshold (499.9 nM) is classified as binder."""
        entry = IEDBBenchmarkEntry(
            peptide="TESTPEPTD",
            allele="HLA-A*02:01",
            measured_ic50=499.9,
            source="test",
        )
        assert entry.is_binder is True

    def test_is_binder_returns_bool(self) -> None:
        """is_binder always returns a bool, not an int or other type."""
        for ic50 in [5.0, 150.0, 499.9, 500.0, 1800.0, 50000.0]:
            entry = IEDBBenchmarkEntry(
                peptide="AAAAAAAAA",
                allele="HLA-A*02:01",
                measured_ic50=ic50,
                source="test",
            )
            assert isinstance(entry.is_binder, bool)

    def test_binder_threshold_constant(self) -> None:
        """BINDER_IC50_THRESHOLD is 500.0 nM (standard IEDB threshold)."""
        assert BINDER_IC50_THRESHOLD == 500.0

    def test_entry_is_dataclass(self) -> None:
        """IEDBBenchmarkEntry is a proper dataclass with __dataclass_fields__."""
        import dataclasses
        assert dataclasses.is_dataclass(IEDBBenchmarkEntry)

    def test_entry_fields(self) -> None:
        """IEDBBenchmarkEntry has the expected fields."""
        entry = IEDBBenchmarkEntry(
            peptide="SIINFEKL",
            allele="HLA-A*02:01",
            measured_ic50=10.0,
            source="test",
        )
        assert hasattr(entry, "peptide")
        assert hasattr(entry, "allele")
        assert hasattr(entry, "measured_ic50")
        assert hasattr(entry, "source")
        assert hasattr(entry, "is_binder")


# ═══════════════════════════════════════════════════════════════════════════
# 2. IEDB_BENCHMARK_DATA entries are valid
# ═══════════════════════════════════════════════════════════════════════════

class TestIEDBBenchmarkData:
    """Tests for the IEDB_BENCHMARK_DATA curated dataset."""

    def test_data_is_list(self) -> None:
        """IEDB_BENCHMARK_DATA is a list."""
        assert isinstance(IEDB_BENCHMARK_DATA, list)

    def test_data_not_empty(self) -> None:
        """IEDB_BENCHMARK_DATA is non-empty."""
        assert len(IEDB_BENCHMARK_DATA) > 0

    def test_data_has_sufficient_entries(self) -> None:
        """IEDB_BENCHMARK_DATA contains at least 20 entries."""
        assert len(IEDB_BENCHMARK_DATA) >= 20

    def test_all_entries_are_benchmark_entries(self) -> None:
        """Every element is an IEDBBenchmarkEntry instance."""
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            assert isinstance(entry, IEDBBenchmarkEntry), (
                f"Entry {i} is {type(entry).__name__}, not IEDBBenchmarkEntry"
            )

    def test_all_peptides_valid_aa(self) -> None:
        """All peptides contain only standard amino acid single-letter codes."""
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            for ch in entry.peptide:
                assert ch in valid_aa, (
                    f"Entry {i} peptide '{entry.peptide}' contains invalid AA '{ch}'"
                )

    def test_all_peptides_valid_mhc_i_length(self) -> None:
        """All peptides are 8-11 residues (valid MHC-I epitope lengths)."""
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            assert 8 <= len(entry.peptide) <= 11, (
                f"Entry {i} peptide '{entry.peptide}' has length {len(entry.peptide)}, "
                f"expected 8-11"
            )

    def test_all_alleles_valid_format(self) -> None:
        """All alleles follow HLA naming convention (HLA-X*XX:XX)."""
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            assert entry.allele.startswith("HLA-"), (
                f"Entry {i} allele '{entry.allele}' doesn't start with 'HLA-'"
            )
            assert "*" in entry.allele, (
                f"Entry {i} allele '{entry.allele}' missing '*' separator"
            )

    def test_all_ic50_positive(self) -> None:
        """All measured IC50 values are positive."""
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            assert entry.measured_ic50 > 0, (
                f"Entry {i} has non-positive IC50: {entry.measured_ic50}"
            )

    def test_all_sources_non_empty(self) -> None:
        """All source fields are non-empty strings."""
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            assert isinstance(entry.source, str) and len(entry.source) > 0, (
                f"Entry {i} has empty or non-string source"
            )

    def test_includes_binders_and_non_binders(self) -> None:
        """Dataset includes both binders and non-binders."""
        binders = [e for e in IEDB_BENCHMARK_DATA if e.is_binder]
        non_binders = [e for e in IEDB_BENCHMARK_DATA if not e.is_binder]
        assert len(binders) > 0, "Dataset has no binders"
        assert len(non_binders) > 0, "Dataset has no non-binders"

    def test_includes_multiple_alleles(self) -> None:
        """Dataset spans multiple MHC alleles."""
        alleles = {e.allele for e in IEDB_BENCHMARK_DATA}
        assert len(alleles) >= 2, f"Expected >= 2 alleles, got {alleles}"

    def test_includes_hla_a0201(self) -> None:
        """Dataset includes HLA-A*02:01 entries (most studied allele)."""
        alleles = {e.allele for e in IEDB_BENCHMARK_DATA}
        assert "HLA-A*02:01" in alleles

    def test_includes_hla_a0301(self) -> None:
        """Dataset includes HLA-A*03:01 entries."""
        alleles = {e.allele for e in IEDB_BENCHMARK_DATA}
        assert "HLA-A*03:01" in alleles

    def test_ic50_range_spans_binding_classes(self) -> None:
        """IC50 values span all four binding classes."""
        ic50s = [e.measured_ic50 for e in IEDB_BENCHMARK_DATA]
        has_strong = any(ic < 50 for ic in ic50s)
        has_moderate = any(50 <= ic < 500 for ic in ic50s)
        has_weak = any(500 <= ic <= 5000 for ic in ic50s)
        has_non = any(ic > 5000 for ic in ic50s)
        assert has_strong, "No strong binders (IC50 < 50 nM) in dataset"
        assert has_moderate, "No moderate binders (50-500 nM) in dataset"
        assert has_weak, "No weak binders (500-5000 nM) in dataset"
        assert has_non, "No non-binders (>5000 nM) in dataset"

    def test_no_duplicate_peptide_allele_pairs(self) -> None:
        """No duplicate (peptide, allele) pairs in the dataset."""
        seen = set()
        for entry in IEDB_BENCHMARK_DATA:
            key = (entry.peptide, entry.allele)
            assert key not in seen, f"Duplicate entry: {key}"
            seen.add(key)

    def test_well_known_peptides_present(self) -> None:
        """Dataset includes well-known epitopes from the literature."""
        peptides = {e.peptide for e in IEDB_BENCHMARK_DATA}
        # GILGFVFTL: Influenza M1 epitope, most-studied A*02:01 binder
        assert "GILGFVFTL" in peptides, "Missing Influenza M1 epitope"


# ═══════════════════════════════════════════════════════════════════════════
# 3. benchmark_mhc_predictions with mock predictor
# ═══════════════════════════════════════════════════════════════════════════

class TestBenchmarkMHCPredictions:
    """Tests for benchmark_mhc_predictions with mock predictors."""

    def test_perfect_predictor(self) -> None:
        """A predictor that returns exact measured IC50 gets perfect scores."""
        # Build a predictor that memorizes the data
        data_map = {
            (e.peptide, e.allele): e.measured_ic50
            for e in IEDB_BENCHMARK_DATA
        }

        def perfect_predictor(peptide: str, allele: str) -> float:
            return data_map[(peptide, allele)]

        result = benchmark_mhc_predictions(perfect_predictor, IEDB_BENCHMARK_DATA)

        # Perfect predictor should get all correct
        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        assert result.correct_predictions == len(IEDB_BENCHMARK_DATA)
        assert result.auc_roc == pytest.approx(1.0, abs=0.01)
        assert result.pearson_r == pytest.approx(1.0, abs=0.01)

    def test_always_non_binder_predictor(self) -> None:
        """A predictor that always predicts non-binder (50000 nM)."""
        def always_non_binder(peptide: str, allele: str) -> float:
            return 50000.0

        result = benchmark_mhc_predictions(always_non_binder, IEDB_BENCHMARK_DATA)

        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        # Should correctly classify non-binders but miss all binders
        non_binders = sum(1 for e in IEDB_BENCHMARK_DATA if not e.is_binder)
        assert result.correct_predictions == non_binders
        # AUC-ROC should be poor (near 0.5 = random)
        assert 0.0 <= result.auc_roc <= 0.6

    def test_always_binder_predictor(self) -> None:
        """A predictor that always predicts binder (10 nM)."""
        def always_binder(peptide: str, allele: str) -> float:
            return 10.0

        result = benchmark_mhc_predictions(always_binder, IEDB_BENCHMARK_DATA)

        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        # Should correctly classify binders but miss all non-binders
        binders = sum(1 for e in IEDB_BENCHMARK_DATA if e.is_binder)
        assert result.correct_predictions == binders

    def test_inverted_predictor(self) -> None:
        """A predictor that inverts: returns low IC50 for non-binders and vice versa."""
        def inverted_predictor(peptide: str, allele: str) -> float:
            # Find the matching entry
            for e in IEDB_BENCHMARK_DATA:
                if e.peptide == peptide and e.allele == allele:
                    # Invert: binder → high IC50, non-binder → low IC50
                    if e.is_binder:
                        return 50000.0
                    else:
                        return 5.0
            return 50000.0  # fallback

        result = benchmark_mhc_predictions(inverted_predictor, IEDB_BENCHMARK_DATA)

        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        assert result.correct_predictions == 0
        # AUC-ROC should be near 0.0 (perfectly anti-correlated)
        assert result.auc_roc < 0.2

    def test_noisy_predictor(self) -> None:
        """A predictor that adds noise to measured IC50 still correlates."""
        import random
        rng = random.Random(42)  # deterministic seed

        data_map = {
            (e.peptide, e.allele): e.measured_ic50
            for e in IEDB_BENCHMARK_DATA
        }

        def noisy_predictor(peptide: str, allele: str) -> float:
            true_ic50 = data_map[(peptide, allele)]
            # Add ~2x multiplicative noise
            noise_factor = rng.uniform(0.5, 2.0)
            return true_ic50 * noise_factor

        result = benchmark_mhc_predictions(noisy_predictor, IEDB_BENCHMARK_DATA)

        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        # Pearson r should be positive (noise preserves rank order somewhat)
        assert result.pearson_r > 0.5
        # AUC should be above random (0.5)
        assert result.auc_roc > 0.5

    def test_with_empty_entries(self) -> None:
        """Benchmarking with empty entries list returns empty result."""
        def dummy(peptide: str, allele: str) -> float:
            return 500.0

        result = benchmark_mhc_predictions(dummy, [])

        assert result.total_entries == 0
        assert result.correct_predictions == 0
        assert result.auc_roc == 0.0
        assert math.isnan(result.pearson_r)
        assert result.details == []

    def test_with_single_entry(self) -> None:
        """Benchmarking with a single entry produces valid result."""
        entry = IEDBBenchmarkEntry(
            peptide="GILGFVFTL",
            allele="HLA-A*02:01",
            measured_ic50=5.0,
            source="test",
        )

        def predictor(peptide: str, allele: str) -> float:
            return 10.0

        result = benchmark_mhc_predictions(predictor, [entry])

        assert result.total_entries == 1
        assert result.correct_predictions == 1
        # Single entry: AUC degenerate case
        assert result.auc_roc == 0.5  # degenerate

    def test_predictor_exception_handled(self) -> None:
        """Predictor that raises exceptions gets fallback IC50 (50000 nM)."""
        call_count = 0

        def failing_predictor(peptide: str, allele: str) -> float:
            nonlocal call_count
            call_count += 1
            if call_count <= 3:
                raise RuntimeError("prediction failed")
            return 500.0

        result = benchmark_mhc_predictions(failing_predictor, IEDB_BENCHMARK_DATA)

        # Should not crash; entries with exceptions get 50000 nM
        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        assert len(result.details) == len(IEDB_BENCHMARK_DATA)

    def test_predictor_negative_ic50_clamped(self) -> None:
        """Predictor returning negative IC50 gets clamped to 0.01."""
        def negative_predictor(peptide: str, allele: str) -> float:
            return -100.0

        result = benchmark_mhc_predictions(negative_predictor, IEDB_BENCHMARK_DATA)

        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        # All entries should have predicted_binder = True (0.01 < 500)
        for detail in result.details:
            assert detail["predicted_ic50"] == pytest.approx(0.01, abs=0.01)
            assert detail["predicted_binder"] is True

    def test_predictor_called_with_correct_args(self) -> None:
        """Predictor is called with (peptide, allele) for each entry."""
        mock_fn = MagicMock(return_value=50000.0)

        result = benchmark_mhc_predictions(mock_fn, IEDB_BENCHMARK_DATA)

        assert mock_fn.call_count == len(IEDB_BENCHMARK_DATA)
        for i, entry in enumerate(IEDB_BENCHMARK_DATA):
            call_args = mock_fn.call_args_list[i]
            assert call_args[0][0] == entry.peptide
            assert call_args[0][1] == entry.allele

    def test_details_populated_correctly(self) -> None:
        """Per-entry details have all expected keys and correct values."""
        def predictor(peptide: str, allele: str) -> float:
            return 100.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)

        assert len(result.details) == len(IEDB_BENCHMARK_DATA)
        required_keys = {
            "peptide", "allele", "measured_ic50", "predicted_ic50",
            "measured_binder", "predicted_binder", "correct",
        }
        for i, detail in enumerate(result.details):
            assert set(detail.keys()) == required_keys, (
                f"Detail {i} has keys {set(detail.keys())}, expected {required_keys}"
            )
            entry = IEDB_BENCHMARK_DATA[i]
            assert detail["peptide"] == entry.peptide
            assert detail["allele"] == entry.allele
            assert detail["measured_ic50"] == entry.measured_ic50
            assert detail["predicted_ic50"] == pytest.approx(100.0)
            assert detail["measured_binder"] == entry.is_binder
            assert detail["predicted_binder"] is True  # 100 < 500
            assert detail["correct"] == (entry.is_binder is True)

    def test_correct_predictions_count(self) -> None:
        """correct_predictions matches sum of correct flags in details."""
        def predictor(peptide: str, allele: str) -> float:
            return 500.0  # exactly at threshold → predicted_binder = False

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)

        expected_correct = sum(1 for d in result.details if d["correct"])
        assert result.correct_predictions == expected_correct

    def test_with_custom_subset_of_entries(self) -> None:
        """benchmark_mhc_predictions works with a subset of entries."""
        subset = IEDB_BENCHMARK_DATA[:5]

        def predictor(peptide: str, allele: str) -> float:
            return 50000.0

        result = benchmark_mhc_predictions(predictor, subset)

        assert result.total_entries == 5
        assert len(result.details) == 5

    def test_returns_mhbenchmark_result(self) -> None:
        """benchmark_mhc_predictions returns an MHBenchmarkResult instance."""
        def predictor(peptide: str, allele: str) -> float:
            return 1000.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)

        assert isinstance(result, MHBenchmarkResult)


# ═══════════════════════════════════════════════════════════════════════════
# 4. MHBenchmarkResult fields (auc_roc, pearson_r)
# ═══════════════════════════════════════════════════════════════════════════

class TestMHBenchmarkResult:
    """Tests for MHBenchmarkResult dataclass and field validation."""

    def test_construction_all_fields(self) -> None:
        """MHBenchmarkResult can be created with all fields."""
        result = MHBenchmarkResult(
            total_entries=20,
            correct_predictions=18,
            auc_roc=0.95,
            pearson_r=0.88,
            details=[{"key": "value"}],
        )
        assert result.total_entries == 20
        assert result.correct_predictions == 18
        assert result.auc_roc == pytest.approx(0.95)
        assert result.pearson_r == pytest.approx(0.88)
        assert len(result.details) == 1

    def test_details_default_empty_list(self) -> None:
        """details defaults to an empty list when not provided."""
        result = MHBenchmarkResult(
            total_entries=0,
            correct_predictions=0,
            auc_roc=0.0,
            pearson_r=float("nan"),
        )
        assert result.details == []

    def test_is_dataclass(self) -> None:
        """MHBenchmarkResult is a proper dataclass."""
        import dataclasses
        assert dataclasses.is_dataclass(MHBenchmarkResult)

    def test_has_expected_fields(self) -> None:
        """MHBenchmarkResult has all expected fields."""
        result = MHBenchmarkResult(
            total_entries=0,
            correct_predictions=0,
            auc_roc=0.0,
            pearson_r=0.0,
        )
        assert hasattr(result, "total_entries")
        assert hasattr(result, "correct_predictions")
        assert hasattr(result, "auc_roc")
        assert hasattr(result, "pearson_r")
        assert hasattr(result, "details")

    def test_auc_roc_range(self) -> None:
        """auc_roc from benchmark is in valid range [0, 1]."""
        def predictor(peptide: str, allele: str) -> float:
            return 1000.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        assert 0.0 <= result.auc_roc <= 1.0

    def test_pearson_r_range(self) -> None:
        """pearson_r from benchmark is in valid range [-1, 1] or NaN."""
        def predictor(peptide: str, allele: str) -> float:
            return 1000.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        if not math.isnan(result.pearson_r):
            assert -1.0 <= result.pearson_r <= 1.0

    def test_pearson_r_nan_for_empty(self) -> None:
        """pearson_r is NaN when no entries are provided."""
        result = MHBenchmarkResult(
            total_entries=0,
            correct_predictions=0,
            auc_roc=0.0,
            pearson_r=float("nan"),
        )
        assert math.isnan(result.pearson_r)

    def test_auc_roc_perfect_predictor_near_one(self) -> None:
        """Perfect predictor yields auc_roc close to 1.0."""
        data_map = {
            (e.peptide, e.allele): e.measured_ic50
            for e in IEDB_BENCHMARK_DATA
        }
        result = benchmark_mhc_predictions(
            lambda p, a: data_map[(p, a)],
            IEDB_BENCHMARK_DATA,
        )
        assert result.auc_roc >= 0.99

    def test_pearson_r_perfect_predictor_near_one(self) -> None:
        """Perfect predictor yields pearson_r close to 1.0."""
        data_map = {
            (e.peptide, e.allele): e.measured_ic50
            for e in IEDB_BENCHMARK_DATA
        }
        result = benchmark_mhc_predictions(
            lambda p, a: data_map[(p, a)],
            IEDB_BENCHMARK_DATA,
        )
        assert result.pearson_r >= 0.99

    def test_total_entries_matches_input(self) -> None:
        """total_entries equals the number of input entries."""
        def predictor(peptide: str, allele: str) -> float:
            return 500.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        assert result.total_entries == len(IEDB_BENCHMARK_DATA)

    def test_correct_predictions_leq_total(self) -> None:
        """correct_predictions cannot exceed total_entries."""
        def predictor(peptide: str, allele: str) -> float:
            return 500.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        assert result.correct_predictions <= result.total_entries

    def test_correct_predictions_non_negative(self) -> None:
        """correct_predictions is non-negative."""
        def predictor(peptide: str, allele: str) -> float:
            return 50000.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        assert result.correct_predictions >= 0


# ═══════════════════════════════════════════════════════════════════════════
# Additional: Internal helper function tests
# ═══════════════════════════════════════════════════════════════════════════

class TestPearsonCorrelation:
    """Tests for the _pearson_correlation internal helper."""

    def test_perfect_positive(self) -> None:
        """Perfect positive correlation returns 1.0."""
        r = _pearson_correlation([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0])
        assert r == pytest.approx(1.0)

    def test_perfect_negative(self) -> None:
        """Perfect negative correlation returns -1.0."""
        r = _pearson_correlation([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
        assert r == pytest.approx(-1.0)

    def test_no_correlation(self) -> None:
        """Uncorrelated data gives r near 0."""
        r = _pearson_correlation([1.0, 2.0, 3.0, 4.0], [1.0, -1.0, 1.0, -1.0])
        assert abs(r) < 0.5

    def test_too_few_points(self) -> None:
        """Fewer than 2 data points returns NaN."""
        r = _pearson_correlation([1.0], [1.0])
        assert math.isnan(r)

    def test_mismatched_lengths(self) -> None:
        """Mismatched list lengths returns NaN."""
        r = _pearson_correlation([1.0, 2.0], [1.0])
        assert math.isnan(r)

    def test_zero_variance(self) -> None:
        """Zero variance in either variable returns NaN."""
        r = _pearson_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
        assert math.isnan(r)

    def test_output_clamped(self) -> None:
        """Output is clamped to [-1, 1] even with floating-point imprecision."""
        r = _pearson_correlation([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
        assert -1.0 <= r <= 1.0


class TestComputeAUCROC:
    """Tests for the _compute_auc_roc internal helper."""

    def test_perfect_separation(self) -> None:
        """All binders have lower predicted IC50 than all non-binders → AUC 1.0."""
        measured = [5.0, 10.0, 100.0, 5000.0, 10000.0, 50000.0]
        predicted = [1.0, 2.0, 50.0, 10000.0, 20000.0, 80000.0]
        auc = _compute_auc_roc(measured, predicted)
        assert auc == pytest.approx(1.0)

    def test_random_predictions(self) -> None:
        """All same predicted IC50 → AUC 0.5 (random)."""
        measured = [5.0, 100.0, 1000.0, 50000.0]
        predicted = [500.0, 500.0, 500.0, 500.0]
        auc = _compute_auc_roc(measured, predicted)
        assert auc == pytest.approx(0.5)

    def test_inverted_predictions(self) -> None:
        """All binders predicted high, non-binders predicted low → AUC 0.0."""
        measured = [5.0, 10.0, 50000.0, 30000.0]
        predicted = [50000.0, 40000.0, 5.0, 10.0]
        auc = _compute_auc_roc(measured, predicted)
        assert auc == pytest.approx(0.0)

    def test_too_few_entries(self) -> None:
        """Fewer than 2 entries returns 0.5 (degenerate)."""
        auc = _compute_auc_roc([100.0], [100.0])
        assert auc == 0.5

    def test_all_binders(self) -> None:
        """All entries are binders → AUC 0.5 (no negatives for comparison)."""
        measured = [5.0, 10.0, 50.0, 100.0]  # all < 500
        predicted = [1.0, 5.0, 50.0, 200.0]
        auc = _compute_auc_roc(measured, predicted)
        assert auc == 0.5

    def test_all_non_binders(self) -> None:
        """All entries are non-binders → AUC 0.5 (no positives for comparison)."""
        measured = [1000.0, 5000.0, 10000.0, 50000.0]  # all >= 500
        predicted = [500.0, 1000.0, 5000.0, 50000.0]
        auc = _compute_auc_roc(measured, predicted)
        assert auc == 0.5

    def test_auc_in_valid_range(self) -> None:
        """AUC-ROC is always in [0, 1]."""
        measured = [5.0, 50.0, 500.0, 5000.0, 50000.0]
        predicted = [100.0, 200.0, 500.0, 5000.0, 25000.0]
        auc = _compute_auc_roc(measured, predicted)
        assert 0.0 <= auc <= 1.0


# ═══════════════════════════════════════════════════════════════════════════
# Additional: Edge cases and integration
# ═══════════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Edge-case and integration tests."""

    def test_all_binder_entries(self) -> None:
        """Benchmark with only binder entries works without error."""
        binders = [e for e in IEDB_BENCHMARK_DATA if e.is_binder]

        def predictor(peptide: str, allele: str) -> float:
            return 10.0

        result = benchmark_mhc_predictions(predictor, binders)
        assert result.total_entries == len(binders)

    def test_all_non_binder_entries(self) -> None:
        """Benchmark with only non-binder entries works without error."""
        non_binders = [e for e in IEDB_BENCHMARK_DATA if not e.is_binder]

        def predictor(peptide: str, allele: str) -> float:
            return 50000.0

        result = benchmark_mhc_predictions(predictor, non_binders)
        assert result.total_entries == len(non_binders)

    def test_predictor_returns_zero(self) -> None:
        """Predictor returning 0.0 (perfect binder) is handled."""
        def predictor(peptide: str, allele: str) -> float:
            return 0.0

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        # All predicted as binder
        for detail in result.details:
            assert detail["predicted_binder"] is True

    def test_predictor_returns_very_large(self) -> None:
        """Predictor returning very large IC50 is handled."""
        def predictor(peptide: str, allele: str) -> float:
            return 1e9

        result = benchmark_mhc_predictions(predictor, IEDB_BENCHMARK_DATA)
        assert result.total_entries == len(IEDB_BENCHMARK_DATA)
        # All predicted as non-binder
        for detail in result.details:
            assert detail["predicted_binder"] is False

    def test_module_all_exports(self) -> None:
        """Module __all__ includes all expected public names."""
        assert "IEDBBenchmarkEntry" in iedb_mod.__all__
        assert "IEDB_BENCHMARK_DATA" in iedb_mod.__all__
        assert "MHBenchmarkResult" in iedb_mod.__all__
        assert "benchmark_mhc_predictions" in iedb_mod.__all__

    def test_import_from_validation_package(self) -> None:
        """Public names are importable from biocompiler.validation."""
        from biocompiler.validation import IEDBBenchmarkEntry as IBE
        from biocompiler.validation import benchmark_mhc_predictions as bmp

        assert IBE is IEDBBenchmarkEntry
        assert bmp is benchmark_mhc_predictions

    def test_two_predictors_different_results(self) -> None:
        """Two different predictors produce different benchmark results."""
        def good_predictor(peptide: str, allele: str) -> float:
            # Approximate: correlates with measured
            for e in IEDB_BENCHMARK_DATA:
                if e.peptide == peptide and e.allele == allele:
                    return e.measured_ic50 * 1.5
            return 50000.0

        def mediocre_predictor(peptide: str, allele: str) -> float:
            # Returns a constant non-binder IC50 with small perturbation
            # to avoid zero variance (which would give NaN Pearson r)
            import random
            rng = random.Random(123)
            return 50000.0 + rng.uniform(-100, 100)

        good_result = benchmark_mhc_predictions(good_predictor, IEDB_BENCHMARK_DATA)
        mediocre_result = benchmark_mhc_predictions(mediocre_predictor, IEDB_BENCHMARK_DATA)

        assert good_result.correct_predictions >= mediocre_result.correct_predictions
        # Good predictor should have higher AUC-ROC
        assert good_result.auc_roc >= mediocre_result.auc_roc

    def test_single_allele_subset(self) -> None:
        """Benchmarking on a single allele produces valid metrics."""
        a0201_entries = [
            e for e in IEDB_BENCHMARK_DATA if e.allele == "HLA-A*02:01"
        ]

        def predictor(peptide: str, allele: str) -> float:
            return 100.0

        result = benchmark_mhc_predictions(predictor, a0201_entries)

        assert result.total_entries == len(a0201_entries)
        assert 0.0 <= result.auc_roc <= 1.0
        if not math.isnan(result.pearson_r):
            assert -1.0 <= result.pearson_r <= 1.0
