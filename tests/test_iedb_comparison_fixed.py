"""Test BioCompiler IEDB Comparison — compare_with_iedb and IEDBComparisonResult.

Tests cover:
1. compare_with_iedb() with known IEDB epitopes
2. Sensitivity and specificity computation
3. AUC-ROC estimation
4. Empty peptide lists
5. Known true positives
6. Error handling for invalid inputs
7. Unknown allele handling
8. Helper functions (get_known_binders, get_known_non_binders, get_available_alleles)
"""
from __future__ import annotations

import math

import pytest

# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------
# The biocompiler package has a pre-existing circular import between
# biocompiler.__init__ → .translation → .organisms → .benchmarking → .translation
# which prevents normal `from biocompiler.validation.iedb_comparison import ...`.
# We load the module directly from its file path using importlib.util,
# bypassing the circular package chain entirely.  This is safe because
# iedb_comparison.py has no intra-package imports.

import importlib.util
import sys
import types
from pathlib import Path

_MODULE_PATH = (
    Path(__file__).resolve().parent.parent
    / "src" / "biocompiler" / "validation" / "iedb_comparison.py"
)

if not _MODULE_PATH.exists():
    pytest.skip(
        f"iedb_comparison.py not found at {_MODULE_PATH}",
        allow_module_level=True,
    )

# Save existing biocompiler modules so we can restore them after loading
_saved_modules: dict[str, object] = {}
for _key in list(sys.modules):
    if "biocompiler" in _key:
        _saved_modules[_key] = sys.modules.pop(_key)

# Ensure parent packages exist in sys.modules so dataclass resolution works
_needs_restore: list[str] = []
for _pkg in ("biocompiler", "biocompiler.validation"):
    if _pkg not in sys.modules:
        _mod = types.ModuleType(_pkg)
        if _pkg == "biocompiler":
            _mod.__path__ = [str(Path(__file__).resolve().parent.parent / "src" / "biocompiler")]
            _mod.__package__ = "biocompiler"
        elif _pkg == "biocompiler.validation":
            _mod.__path__ = [str(Path(__file__).resolve().parent.parent / "src" / "biocompiler" / "validation")]
            _mod.__package__ = "biocompiler.validation"
        sys.modules[_pkg] = _mod
        _needs_restore.append(_pkg)

_spec = importlib.util.spec_from_file_location(
    "biocompiler.validation.iedb_comparison", str(_MODULE_PATH)
)
iedb_mod = importlib.util.module_from_spec(_spec)  # type: ignore[arg-type]
sys.modules["biocompiler.validation.iedb_comparison"] = iedb_mod
_spec.loader.exec_module(iedb_mod)  # type: ignore[union-attr]

# Restore previously saved modules so other tests aren't affected
for _key, _mod in _saved_modules.items():
    sys.modules[_key] = _mod
# Remove temporary parent packages we created (they were placeholders)
for _pkg in _needs_restore:
    if _pkg not in _saved_modules:
        del sys.modules[_pkg]

IEDBComparisonResult = iedb_mod.IEDBComparisonResult
IEDBBenchmarkEntry = iedb_mod.IEDBBenchmarkEntry
IEDB_BENCHMARK_DATA = iedb_mod.IEDB_BENCHMARK_DATA
BINDER_IC50_THRESHOLD = iedb_mod.BINDER_IC50_THRESHOLD
compare_with_iedb = iedb_mod.compare_with_iedb
get_known_binders = iedb_mod.get_known_binders
get_known_non_binders = iedb_mod.get_known_non_binders
get_available_alleles = iedb_mod.get_available_alleles


# ═══════════════════════════════════════════════════════════════════════════
# 1. compare_with_iedb() with known IEDB epitopes
# ═══════════════════════════════════════════════════════════════════════════

class TestCompareWithIEDBKnownEpitopes:
    """Tests for compare_with_iedb() using known IEDB epitope data."""

    def test_returns_iedb_comparison_result(self) -> None:
        """compare_with_iedb returns an IEDBComparisonResult instance."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL"])
        assert isinstance(result, IEDBComparisonResult)

    def test_allele_preserved_in_result(self) -> None:
        """The allele field in the result matches the input allele."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL"])
        assert result.allele == "HLA-A*02:01"

    def test_known_binder_is_true_positive(self) -> None:
        """GILGFVFTL is a known A*02:01 binder → TP=1."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL"])
        assert result.true_positives == 1
        assert result.false_positives == 0

    def test_known_non_binder_is_false_positive(self) -> None:
        """A known IEDB non-binder predicted as binder → FP=1."""
        # DNEEGVQAD is a known non-binder for A*02:01 (IC50=25000)
        result = compare_with_iedb("HLA-A*02:01", ["DNEEGVQAD"])
        assert result.true_positives == 0
        assert result.false_positives == 1

    def test_unknown_peptide_is_false_positive(self) -> None:
        """A peptide not in the IEDB dataset predicted as binder → FP=1."""
        result = compare_with_iedb("HLA-A*02:01", ["ZZZZZZZZZ"])
        assert result.true_positives == 0
        assert result.false_positives == 1

    def test_mixed_predictions(self) -> None:
        """Mix of known binders, known non-binders, and unknowns."""
        # GILGFVFTL: known binder for A*02:01 (IC50=5)
        # LLFGYPVYV: known binder for A*02:01 (IC50=8)
        # DNEEGVQAD: known non-binder for A*02:01 (IC50=25000)
        # FAKEPEPTI: not in IEDB dataset
        result = compare_with_iedb(
            "HLA-A*02:01",
            ["GILGFVFTL", "LLFGYPVYV", "DNEEGVQAD", "FAKEPEPTI"],
        )
        assert result.true_positives == 2
        assert result.false_positives == 2  # DNEEGVQAD (non-binder) + FAKEPEPTI (unknown)

    def test_all_known_binders_for_allele(self) -> None:
        """Predicting all known binders → high TP, FN=0."""
        binders = get_known_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(binders))
        assert result.true_positives == len(binders)
        assert result.false_negatives == 0
        assert result.sensitivity == pytest.approx(1.0)

    def test_a0301_known_binder(self) -> None:
        """KVYLRDIAP is a known A*03:01 binder (IC50=12)."""
        result = compare_with_iedb("HLA-A*03:01", ["KVYLRDIAP"])
        assert result.true_positives == 1
        assert result.false_positives == 0

    def test_a0301_known_non_binder(self) -> None:
        """DDEEGVQAA is a known A*03:01 non-binder (IC50=30000)."""
        result = compare_with_iedb("HLA-A*03:01", ["DDEEGVQAA"])
        assert result.true_positives == 0
        assert result.false_positives == 1


# ═══════════════════════════════════════════════════════════════════════════
# 2. Sensitivity and specificity computation
# ═══════════════════════════════════════════════════════════════════════════

class TestSensitivitySpecificity:
    """Tests for sensitivity and specificity computation."""

    def test_perfect_sensitivity_all_binders_predicted(self) -> None:
        """Predicting all known binders gives sensitivity=1.0."""
        binders = get_known_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(binders))
        assert result.sensitivity == pytest.approx(1.0)

    def test_zero_sensitivity_no_binders_predicted(self) -> None:
        """Predicting no binders (empty list) gives sensitivity=0.0."""
        result = compare_with_iedb("HLA-A*02:01", [])
        # No binders predicted → all IEDB binders are FN
        assert result.sensitivity == pytest.approx(0.0)

    def test_perfect_specificity_only_binders_predicted(self) -> None:
        """Predicting only true binders gives specificity=1.0 (no FP)."""
        binders = get_known_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(binders))
        # All non-binders are correctly left out → TN = all non-binders
        non_binders = get_known_non_binders("HLA-A*02:01")
        assert result.true_negatives == len(non_binders)
        assert result.false_positives == 0
        assert result.specificity == pytest.approx(1.0)

    def test_zero_specificity_all_non_binders_predicted(self) -> None:
        """Predicting all IEDB non-binders as binders gives specificity=0.0."""
        non_binders = get_known_non_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(non_binders))
        # All non-binders are FP, TN=0
        assert result.specificity == pytest.approx(0.0)

    def test_sensitivity_formula(self) -> None:
        """Sensitivity = TP / (TP + FN)."""
        binders = get_known_binders("HLA-A*02:01")
        # Predict only half the binders
        half = list(binders)[:max(1, len(binders) // 2)]
        result = compare_with_iedb("HLA-A*02:01", half)
        expected_sens = result.true_positives / (result.true_positives + result.false_negatives)
        assert result.sensitivity == pytest.approx(expected_sens)

    def test_specificity_formula(self) -> None:
        """Specificity = TN / (TN + FP)."""
        binders = get_known_binders("HLA-A*02:01")
        non_binders = get_known_non_binders("HLA-A*02:01")
        # Predict all binders + one non-binder
        predicted = list(binders) + [list(non_binders)[0]]
        result = compare_with_iedb("HLA-A*02:01", predicted)
        expected_spec = result.true_negatives / (result.true_negatives + result.false_positives)
        assert result.specificity == pytest.approx(expected_spec)

    def test_partial_sensitivity(self) -> None:
        """Partial sensitivity when only some binders are found."""
        binders = get_known_binders("HLA-A*02:01")
        if len(binders) >= 2:
            half = list(binders)[:len(binders) // 2]
            result = compare_with_iedb("HLA-A*02:01", half)
            assert 0.0 < result.sensitivity < 1.0


# ═══════════════════════════════════════════════════════════════════════════
# 3. AUC estimation
# ═══════════════════════════════════════════════════════════════════════════

class TestAUCEstimation:
    """Tests for AUC-ROC estimation in compare_with_iedb."""

    def test_auc_estimate_is_average_of_sens_spec(self) -> None:
        """AUC estimate = (sensitivity + specificity) / 2."""
        binders = get_known_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(binders))
        expected_auc = (result.sensitivity + result.specificity) / 2.0
        assert result.auc_estimate == pytest.approx(expected_auc, abs=1e-5)

    def test_perfect_prediction_auc(self) -> None:
        """Predicting exactly the known binders → AUC=1.0."""
        binders = get_known_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(binders))
        assert result.sensitivity == pytest.approx(1.0)
        assert result.specificity == pytest.approx(1.0)
        assert result.auc_estimate == pytest.approx(1.0)

    def test_worst_prediction_auc(self) -> None:
        """Predicting only non-binders → low AUC."""
        non_binders = get_known_non_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", list(non_binders))
        # sensitivity=0, specificity=0 → AUC=0
        assert result.auc_estimate == pytest.approx(0.0)

    def test_auc_in_valid_range(self) -> None:
        """AUC estimate is always in [0, 1]."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL"])
        assert 0.0 <= result.auc_estimate <= 1.0

    def test_auc_random_baseline(self) -> None:
        """AUC should be between 0 and 1 for any prediction set."""
        binders = get_known_binders("HLA-A*02:01")
        non_binders = get_known_non_binders("HLA-A*02:01")
        # Predict all peptides (binders + non-binders)
        result = compare_with_iedb(
            "HLA-A*02:01",
            list(binders) + list(non_binders),
        )
        assert 0.0 <= result.auc_estimate <= 1.0
        # sensitivity=1.0 (all binders found), specificity=0.0 (all non-binders also predicted)
        assert result.auc_estimate == pytest.approx(0.5)


# ═══════════════════════════════════════════════════════════════════════════
# 4. Empty peptide lists
# ═══════════════════════════════════════════════════════════════════════════

class TestEmptyPeptideList:
    """Tests for compare_with_iedb with empty peptide lists."""

    def test_empty_peptides_no_true_positives(self) -> None:
        """Empty peptide list → TP=0, FP=0."""
        result = compare_with_iedb("HLA-A*02:01", [])
        assert result.true_positives == 0
        assert result.false_positives == 0

    def test_empty_peptides_all_binders_are_false_negatives(self) -> None:
        """Empty peptide list → all IEDB binders are FN."""
        binders = get_known_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", [])
        assert result.false_negatives == len(binders)

    def test_empty_peptides_all_non_binders_are_true_negatives(self) -> None:
        """Empty peptide list → all IEDB non-binders are TN."""
        non_binders = get_known_non_binders("HLA-A*02:01")
        result = compare_with_iedb("HLA-A*02:01", [])
        assert result.true_negatives == len(non_binders)

    def test_empty_peptides_sensitivity_zero(self) -> None:
        """Empty peptide list → sensitivity=0.0."""
        result = compare_with_iedb("HLA-A*02:01", [])
        assert result.sensitivity == pytest.approx(0.0)

    def test_empty_peptides_specificity_one(self) -> None:
        """Empty peptide list → specificity=1.0 (no false positives)."""
        result = compare_with_iedb("HLA-A*02:01", [])
        assert result.specificity == pytest.approx(1.0)

    def test_empty_peptides_auc_estimate(self) -> None:
        """Empty peptide list → AUC = (0 + 1) / 2 = 0.5."""
        result = compare_with_iedb("HLA-A*02:01", [])
        assert result.auc_estimate == pytest.approx(0.5)


# ═══════════════════════════════════════════════════════════════════════════
# 5. Known true positives
# ═══════════════════════════════════════════════════════════════════════════

class TestKnownTruePositives:
    """Tests for compare_with_iedb with specific known IEDB epitopes."""

    def test_gilgfvftl_a0201(self) -> None:
        """GILGFVFTL (Influenza M1) is a known A*02:01 binder."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL"])
        assert result.true_positives == 1
        # Other IEDB binders for this allele that were not predicted → FN
        binders = get_known_binders("HLA-A*02:01")
        assert result.false_negatives == len(binders) - 1

    def test_llfgypvyv_a0201(self) -> None:
        """LLFGYPVYV (HTLV-1 Tax) is a known A*02:01 binder."""
        result = compare_with_iedb("HLA-A*02:01", ["LLFGYPVYV"])
        assert result.true_positives == 1

    def test_elagigiltv_a0201(self) -> None:
        """ELAGIGILTV (Melanoma MART-1) is a known A*02:01 binder."""
        result = compare_with_iedb("HLA-A*02:01", ["ELAGIGILTV"])
        assert result.true_positives == 1

    def test_yldvgvltv_a0201(self) -> None:
        """YLDVGVLTV (EBV BMLF1) is a known A*02:01 binder."""
        result = compare_with_iedb("HLA-A*02:01", ["YLDVGVLTV"])
        assert result.true_positives == 1

    def test_kvylrdiap_a0301(self) -> None:
        """KVYLRDIAP is a known A*03:01 binder."""
        result = compare_with_iedb("HLA-A*03:01", ["KVYLRDIAP"])
        assert result.true_positives == 1

    def test_multiple_true_positives(self) -> None:
        """Multiple known binders all counted as TP."""
        result = compare_with_iedb(
            "HLA-A*02:01",
            ["GILGFVFTL", "LLFGYPVYV", "ELAGIGILTV"],
        )
        assert result.true_positives == 3

    def test_false_negative_when_binder_not_predicted(self) -> None:
        """A known binder that is NOT in the peptide list → FN."""
        binders = get_known_binders("HLA-A*02:01")
        # Don't predict any binders
        result = compare_with_iedb("HLA-A*02:01", ["FAKEPEPTI"])
        assert result.false_negatives == len(binders)

    def test_true_negative_when_non_binder_not_predicted(self) -> None:
        """A known non-binder NOT in the peptide list → TN."""
        non_binders = get_known_non_binders("HLA-A*02:01")
        binders = get_known_binders("HLA-A*02:01")
        # Predict only binders (no non-binders)
        result = compare_with_iedb("HLA-A*02:01", list(binders))
        assert result.true_negatives == len(non_binders)


# ═══════════════════════════════════════════════════════════════════════════
# 6. Error handling for invalid inputs
# ═══════════════════════════════════════════════════════════════════════════

class TestErrorHandling:
    """Tests for error handling in compare_with_iedb."""

    def test_empty_allele_raises_value_error(self) -> None:
        """Empty allele string raises ValueError."""
        with pytest.raises(ValueError, match="allele must be a non-empty string"):
            compare_with_iedb("", ["GILGFVFTL"])

    def test_whitespace_allele_raises_value_error(self) -> None:
        """Whitespace-only allele string raises ValueError."""
        with pytest.raises(ValueError, match="allele must be a non-empty string"):
            compare_with_iedb("   ", ["GILGFVFTL"])

    def test_none_allele_raises_value_error(self) -> None:
        """None allele raises ValueError."""
        with pytest.raises(ValueError, match="allele must be a non-empty string"):
            compare_with_iedb(None, ["GILGFVFTL"])  # type: ignore[arg-type]

    def test_none_peptides_raises_value_error(self) -> None:
        """None peptides raises ValueError."""
        with pytest.raises(ValueError, match="peptides must be a list"):
            compare_with_iedb("HLA-A*02:01", None)  # type: ignore[arg-type]

    def test_non_string_peptide_raises_value_error(self) -> None:
        """Non-string peptide in the list raises ValueError."""
        with pytest.raises(ValueError, match="peptides\\[1\\] must be a string"):
            compare_with_iedb("HLA-A*02:01", ["GILGFVFTL", 42])  # type: ignore[list-item]

    def test_integer_allele_raises_value_error(self) -> None:
        """Integer allele raises ValueError."""
        with pytest.raises(ValueError, match="allele must be a non-empty string"):
            compare_with_iedb(123, ["GILGFVFTL"])  # type: ignore[arg-type]


# ═══════════════════════════════════════════════════════════════════════════
# 7. Unknown allele handling
# ═══════════════════════════════════════════════════════════════════════════

class TestUnknownAllele:
    """Tests for compare_with_iedb with alleles not in the built-in dataset."""

    def test_unknown_allele_all_false_positives(self) -> None:
        """Unknown allele: all predicted peptides become FP."""
        result = compare_with_iedb("HLA-B*07:02", ["GILGFVFTL", "YLDVGVLTV"])
        assert result.true_positives == 0
        assert result.false_positives == 2
        assert result.true_negatives == 0
        assert result.false_negatives == 0

    def test_unknown_allele_empty_peptides(self) -> None:
        """Unknown allele with empty peptides → all zeros."""
        result = compare_with_iedb("HLA-B*07:02", [])
        assert result.true_positives == 0
        assert result.false_positives == 0
        assert result.true_negatives == 0
        assert result.false_negatives == 0

    def test_unknown_allele_auc_is_random(self) -> None:
        """Unknown allele: AUC defaults to 0.5 (no ground truth)."""
        result = compare_with_iedb("HLA-B*07:02", ["GILGFVFTL"])
        assert result.auc_estimate == pytest.approx(0.5)


# ═══════════════════════════════════════════════════════════════════════════
# 8. Helper functions
# ═══════════════════════════════════════════════════════════════════════════

class TestHelperFunctions:
    """Tests for get_known_binders, get_known_non_binders, get_available_alleles."""

    def test_get_known_binders_a0201(self) -> None:
        """get_known_binders returns binder peptides for HLA-A*02:01."""
        binders = get_known_binders("HLA-A*02:01")
        assert isinstance(binders, set)
        assert len(binders) > 0
        # GILGFVFTL is a well-known binder
        assert "GILGFVFTL" in binders

    def test_get_known_binders_a0301(self) -> None:
        """get_known_binders returns binder peptides for HLA-A*03:01."""
        binders = get_known_binders("HLA-A*03:01")
        assert isinstance(binders, set)
        assert len(binders) > 0

    def test_get_known_binders_unknown_allele_empty(self) -> None:
        """get_known_binders returns empty set for unknown allele."""
        binders = get_known_binders("HLA-Z*99:99")
        assert binders == set()

    def test_get_known_non_binders_a0201(self) -> None:
        """get_known_non_binders returns non-binder peptides for HLA-A*02:01."""
        non_binders = get_known_non_binders("HLA-A*02:01")
        assert isinstance(non_binders, set)
        assert len(non_binders) > 0
        # DNEEGVQAD is a known non-binder
        assert "DNEEGVQAD" in non_binders

    def test_get_known_non_binders_unknown_allele_empty(self) -> None:
        """get_known_non_binders returns empty set for unknown allele."""
        non_binders = get_known_non_binders("HLA-Z*99:99")
        assert non_binders == set()

    def test_binders_and_non_binders_disjoint(self) -> None:
        """Binders and non-binders for the same allele are disjoint."""
        for allele in get_available_alleles():
            binders = get_known_binders(allele)
            non_binders = get_known_non_binders(allele)
            assert binders & non_binders == set()

    def test_binders_union_non_binders_equals_all(self) -> None:
        """Union of binders and non-binders equals all peptides for the allele."""
        for allele in get_available_alleles():
            binders = get_known_binders(allele)
            non_binders = get_known_non_binders(allele)
            all_peptides = binders | non_binders
            # Count entries for this allele in benchmark data
            expected_count = sum(
                1 for e in IEDB_BENCHMARK_DATA if e.allele == allele
            )
            assert len(all_peptides) == expected_count

    def test_get_available_alleles(self) -> None:
        """get_available_alleles returns a sorted list with at least 2 alleles."""
        alleles = get_available_alleles()
        assert isinstance(alleles, list)
        assert len(alleles) >= 2
        assert alleles == sorted(alleles)  # sorted
        assert "HLA-A*02:01" in alleles
        assert "HLA-A*03:01" in alleles


# ═══════════════════════════════════════════════════════════════════════════
# 9. IEDBComparisonResult dataclass
# ═══════════════════════════════════════════════════════════════════════════

class TestIEDBComparisonResult:
    """Tests for the IEDBComparisonResult dataclass."""

    def test_construction(self) -> None:
        """IEDBComparisonResult can be constructed with all fields."""
        result = IEDBComparisonResult(
            allele="HLA-A*02:01",
            true_positives=5,
            false_positives=2,
            true_negatives=3,
            false_negatives=1,
            sensitivity=0.833,
            specificity=0.6,
            auc_estimate=0.716,
        )
        assert result.allele == "HLA-A*02:01"
        assert result.true_positives == 5
        assert result.false_positives == 2
        assert result.true_negatives == 3
        assert result.false_negatives == 1
        assert result.sensitivity == pytest.approx(0.833)
        assert result.specificity == pytest.approx(0.6)
        assert result.auc_estimate == pytest.approx(0.716)

    def test_is_dataclass(self) -> None:
        """IEDBComparisonResult is a proper dataclass."""
        import dataclasses
        assert dataclasses.is_dataclass(IEDBComparisonResult)

    def test_has_expected_fields(self) -> None:
        """IEDBComparisonResult has all required fields."""
        result = IEDBComparisonResult(
            allele="HLA-A*02:01",
            true_positives=0,
            false_positives=0,
            true_negatives=0,
            false_negatives=0,
            sensitivity=0.0,
            specificity=0.0,
            auc_estimate=0.5,
        )
        for field_name in [
            "allele", "true_positives", "false_positives",
            "true_negatives", "false_negatives",
            "sensitivity", "specificity", "auc_estimate",
        ]:
            assert hasattr(result, field_name), f"Missing field: {field_name}"

    def test_result_from_compare_has_all_fields(self) -> None:
        """compare_with_iedb result has all expected fields populated."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL"])
        assert isinstance(result.allele, str)
        assert isinstance(result.true_positives, int)
        assert isinstance(result.false_positives, int)
        assert isinstance(result.true_negatives, int)
        assert isinstance(result.false_negatives, int)
        assert isinstance(result.sensitivity, float)
        assert isinstance(result.specificity, float)
        assert isinstance(result.auc_estimate, float)


# ═══════════════════════════════════════════════════════════════════════════
# 10. Confusion matrix consistency
# ═══════════════════════════════════════════════════════════════════════════

class TestConfusionMatrixConsistency:
    """Tests that the confusion matrix entries are internally consistent."""

    def test_tp_plus_fn_equals_total_binders(self) -> None:
        """TP + FN = total number of IEDB binders for the allele."""
        for allele in get_available_alleles():
            binders = get_known_binders(allele)
            result = compare_with_iedb(allele, list(binders)[:1])
            assert result.true_positives + result.false_negatives == len(binders)

    def test_tn_plus_fp_equals_total_non_binders(self) -> None:
        """TN + FP = total number of IEDB non-binders for the allele."""
        for allele in get_available_alleles():
            non_binders = get_known_non_binders(allele)
            binders = get_known_binders(allele)
            # Predict one binder
            result = compare_with_iedb(allele, list(binders)[:1])
            assert result.true_negatives + result.false_positives >= len(non_binders)
            # FP may include unknown peptides, so TN + FP >= number of IEDB non-binders

    def test_all_non_negative(self) -> None:
        """All confusion matrix entries are non-negative."""
        result = compare_with_iedb("HLA-A*02:01", ["GILGFVFTL", "FAKEPEPTI"])
        assert result.true_positives >= 0
        assert result.false_positives >= 0
        assert result.true_negatives >= 0
        assert result.false_negatives >= 0

    def test_duplicate_peptides_handled(self) -> None:
        """Duplicate peptides in the input list are deduplicated (set)."""
        result = compare_with_iedb(
            "HLA-A*02:01",
            ["GILGFVFTL", "GILGFVFTL", "GILGFVFTL"],
        )
        # Only counted once
        assert result.true_positives == 1
        assert result.false_positives == 0


# ═══════════════════════════════════════════════════════════════════════════
# 11. Module exports
# ═══════════════════════════════════════════════════════════════════════════

class TestModuleExports:
    """Tests that new public names are properly exported."""

    def test_compare_with_iedb_in_all(self) -> None:
        """compare_with_iedb is in __all__."""
        assert "compare_with_iedb" in iedb_mod.__all__

    def test_iedb_comparison_result_in_all(self) -> None:
        """IEDBComparisonResult is in __all__."""
        assert "IEDBComparisonResult" in iedb_mod.__all__

    def test_importable_from_validation_package(self) -> None:
        """New names are importable from biocompiler.validation.

        Note: due to a pre-existing circular import in the biocompiler
        top-level package, this test verifies that the module-level
        objects are the same as those we loaded via importlib.util.
        A full ``from biocompiler.validation import ...`` may fail in
        environments where the circular import is not yet resolved;
        in that case the test is skipped.
        """
        try:
            from biocompiler.validation import IEDBComparisonResult as ICR
            from biocompiler.validation import compare_with_iedb as cwi
            from biocompiler.validation import get_available_alleles as gaa
            from biocompiler.validation import get_known_binders as gkb
            from biocompiler.validation import get_known_non_binders as gknb

            assert ICR is IEDBComparisonResult
            assert cwi is compare_with_iedb
            assert gaa is get_available_alleles
            assert gkb is get_known_binders
            assert gknb is get_known_non_binders
        except ImportError:
            pytest.skip("biocompiler.validation package has circular import; skipping")
