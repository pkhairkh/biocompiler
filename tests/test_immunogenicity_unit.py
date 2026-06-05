"""Unit tests for biocompiler.immunogenicity — focused on pure-function coverage.

Covers four areas:
  1. classify_binding() — correct IC50 → binding-class mapping
  2. MHCBindingResult — dataclass construction & field validation
  3. binding_score_to_ic50() — forward/inverse round-trip & boundary behaviour
  4. Basic scanning functions — predict_mhc_i_binding / predict_mhc_ii_binding /
     predict_t_cell_epitopes / predict_b_cell_epitopes return reasonable results
"""

from __future__ import annotations

import math

import pytest

from biocompiler.immunogenicity import (
    MHCBindingResult,
    MHCPredictionResult,
    binding_score_to_ic50,
    classify_binding,
    clear_cache,
    predict_b_cell_epitopes,
    predict_mhc_i_binding,
    predict_mhc_ii_binding,
    predict_t_cell_epitopes,
    score_peptide_pssm,
)

# ---------------------------------------------------------------------------
# Test proteins
# ---------------------------------------------------------------------------

# Known HLA-A*02:01 binder from influenza M1 (GILGFVFTL)
KNOWN_BINDER_9MER = "GILGFVFTL"

# Extended version with flanking residues for 15-mer scanning
SHORT_PROTEIN = "AGILGFVFTLAGILGF"

# Human hemoglobin beta chain
HBB_HUMAN = (
    "MVHLTPEEKSAVTALWGKVNVADIVGHALSDLHAKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. classify_binding
# ═══════════════════════════════════════════════════════════════════════════

class TestClassifyBinding:
    """Test classify_binding() returns the correct category for every IC50 tier."""

    # -- Strong binders: IC50 < 50 nM --
    @pytest.mark.parametrize("ic50", [0.1, 1, 10, 25, 49, 49.99])
    def test_strong_binder(self, ic50: float):
        assert classify_binding(ic50) == "strong_binder", (
            f"IC50={ic50} should be strong_binder"
        )

    # -- Moderate binders: 50 <= IC50 <= 500 nM --
    @pytest.mark.parametrize("ic50", [50, 100, 250, 499, 500])
    def test_moderate_binder(self, ic50: float):
        assert classify_binding(ic50) == "moderate_binder", (
            f"IC50={ic50} should be moderate_binder"
        )

    # -- Weak binders: 500 < IC50 <= 5000 nM --
    @pytest.mark.parametrize("ic50", [501, 1000, 2500, 5000])
    def test_weak_binder(self, ic50: float):
        assert classify_binding(ic50) == "weak_binder", (
            f"IC50={ic50} should be weak_binder"
        )

    # -- Non-binders: IC50 > 5000 nM --
    @pytest.mark.parametrize("ic50", [5001, 10000, 100000])
    def test_non_binder(self, ic50: float):
        assert classify_binding(ic50) == "non_binder", (
            f"IC50={ic50} should be non_binder"
        )

    # -- Boundary: exact thresholds --
    def test_exact_threshold_50(self):
        """IC50 = 50 is the lower bound of moderate_binder."""
        assert classify_binding(50) == "moderate_binder"

    def test_exact_threshold_500(self):
        """IC50 = 500 is the upper bound of moderate_binder."""
        assert classify_binding(500) == "moderate_binder"

    def test_exact_threshold_5000(self):
        """IC50 = 5000 is the upper bound of weak_binder."""
        assert classify_binding(5000) == "weak_binder"

    # -- Off-by-one at boundaries --
    def test_boundary_just_below_50(self):
        assert classify_binding(49.999) == "strong_binder"

    def test_boundary_just_above_500(self):
        assert classify_binding(500.001) == "weak_binder"

    def test_boundary_just_above_5000(self):
        assert classify_binding(5000.001) == "non_binder"

    # -- Return type is always a string from the valid set --
    def test_return_value_is_valid_string(self):
        valid = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        for ic50 in [0.01, 1, 49, 50, 200, 500, 1000, 5000, 5001, 1e6]:
            assert classify_binding(ic50) in valid


# ═══════════════════════════════════════════════════════════════════════════
# 2. MHCBindingResult construction
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCBindingResult:
    """Test MHCBindingResult dataclass construction and field access."""

    def test_basic_construction(self):
        """Construct with all required fields; verify attribute access."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="GILGFVFTL",
            start_position=0,
            end_position=8,
            binding_score=0.85,
            ic50_nm=45.3,
            binding_class="strong_binder",
            anchor_residues={1: "I", 8: "L"},
            anchor_scores={1: 2.0, 8: 1.5},
        )
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "GILGFVFTL"
        assert result.start_position == 0
        assert result.end_position == 8
        assert result.binding_score == 0.85
        assert result.ic50_nm == 45.3
        assert result.binding_class == "strong_binder"
        assert result.anchor_residues == {1: "I", 8: "L"}
        assert result.anchor_scores == {1: 2.0, 8: 1.5}

    def test_construction_with_none_ic50(self):
        """ic50_nm can be None when score is not computable."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="AAAAAAAAA",
            start_position=0,
            end_position=8,
            binding_score=0.0,
            ic50_nm=None,
            binding_class="non_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.ic50_nm is None

    def test_construction_empty_anchor_dicts(self):
        """Peptides with no high-selectivity positions should have empty anchor dicts."""
        result = MHCBindingResult(
            allele="HLA-B*08:01",
            peptide="AAAAAAAAA",
            start_position=5,
            end_position=13,
            binding_score=0.1,
            ic50_nm=7000.0,
            binding_class="non_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.anchor_residues == {}
        assert result.anchor_scores == {}

    def test_is_dataclass(self):
        """MHCBindingResult should be a dataclass."""
        from dataclasses import is_dataclass
        assert is_dataclass(MHCBindingResult)

    def test_binding_class_is_string(self):
        """binding_class should always be a string."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="GILGFVFTL",
            start_position=0,
            end_position=8,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert isinstance(result.binding_class, str)

    def test_binding_score_range_semantic(self):
        """binding_score should be in [0, 1] by convention."""
        # We construct both extremes to verify the dataclass accepts them
        for score in (0.0, 0.5, 1.0):
            result = MHCBindingResult(
                allele="HLA-A*02:01",
                peptide="GILGFVFTL",
                start_position=0,
                end_position=8,
                binding_score=score,
                ic50_nm=500.0,
                binding_class="moderate_binder",
                anchor_residues={},
                anchor_scores={},
            )
            assert 0.0 <= result.binding_score <= 1.0

    def test_position_indices_non_negative(self):
        """start_position and end_position should be non-negative integers."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="GILGFVFTL",
            start_position=0,
            end_position=8,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.start_position >= 0
        assert result.end_position >= 0
        assert result.end_position >= result.start_position


# ═══════════════════════════════════════════════════════════════════════════
# 3. binding_score_to_ic50 and inverse
# ═══════════════════════════════════════════════════════════════════════════

class TestBindingScoreIc50:
    """Test binding_score_to_ic50() forward mapping and inverse round-trip."""

    # -- Forward mapping calibration points (from docstring) --
    def test_high_score_gives_low_ic50(self):
        """Score ~0.9 should give IC50 near 50 nM (strong binder range)."""
        ic50 = binding_score_to_ic50(0.9)
        # 10^(3.949 - 2.5*0.9) = 10^(3.949 - 2.25) = 10^1.699 ≈ 50.0
        assert 20 < ic50 < 100, f"Score 0.9 → IC50={ic50:.1f} nM, expected ~50"

    def test_mid_score_gives_mid_ic50(self):
        """Score ~0.5 should give IC50 near 500 nM (moderate binder range)."""
        ic50 = binding_score_to_ic50(0.5)
        # 10^(3.949 - 2.5*0.5) = 10^(3.949 - 1.25) = 10^2.699 ≈ 500
        assert 300 < ic50 < 800, f"Score 0.5 → IC50={ic50:.1f} nM, expected ~500"

    def test_low_score_gives_high_ic50(self):
        """Score ~0.1 should give IC50 near 5000 nM (weak binder range)."""
        ic50 = binding_score_to_ic50(0.1)
        # 10^(3.949 - 2.5*0.1) = 10^(3.949 - 0.25) = 10^3.699 ≈ 5000
        assert 3000 < ic50 < 8000, f"Score 0.1 → IC50={ic50:.1f} nM, expected ~5000"

    # -- Monotonicity --
    def test_monotonically_decreasing(self):
        """Higher binding scores should produce lower IC50 values."""
        scores = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]
        ic50s = [binding_score_to_ic50(s) for s in scores]
        for i in range(len(ic50s) - 1):
            assert ic50s[i] > ic50s[i + 1], (
                f"IC50 not decreasing: score={scores[i]} → {ic50s[i]:.1f}, "
                f"score={scores[i+1]} → {ic50s[i+1]:.1f}"
            )

    # -- Boundary clamping --
    def test_score_below_zero_clamped(self):
        """Score < 0 should be clamped to 0 before mapping."""
        assert binding_score_to_ic50(-1.0) == binding_score_to_ic50(0.0)

    def test_score_above_one_clamped(self):
        """Score > 1 should be clamped to 1 before mapping."""
        assert binding_score_to_ic50(2.0) == binding_score_to_ic50(1.0)

    # -- Exact formula verification --
    def test_exact_formula(self):
        """Verify the formula IC50 = 10^(3.949 - 2.5 * score)."""
        for score in [0.0, 0.25, 0.5, 0.75, 1.0]:
            clamped = max(0.0, min(1.0, score))
            expected = 10.0 ** (3.949 - 2.5 * clamped)
            actual = binding_score_to_ic50(score)
            assert abs(actual - expected) < 0.01, (
                f"Formula mismatch at score={score}: expected={expected:.2f}, actual={actual:.2f}"
            )

    # -- Inverse round-trip --
    @staticmethod
    def _ic50_to_binding_score(ic50: float) -> float:
        """Inverse of binding_score_to_ic50: score = (3.949 - log10(IC50)) / 2.5."""
        if ic50 <= 0:
            return 1.0
        raw = (3.949 - math.log10(ic50)) / 2.5
        return max(0.0, min(1.0, raw))

    @pytest.mark.parametrize("score", [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    def test_inverse_round_trip(self, score: float):
        """score → ic50 → score should return the original value (within float tolerance)."""
        ic50 = binding_score_to_ic50(score)
        recovered = self._ic50_to_binding_score(ic50)
        assert abs(recovered - score) < 1e-10, (
            f"Round-trip failed: {score} → IC50={ic50:.2f} → {recovered}"
        )

    @pytest.mark.parametrize("ic50", [50, 500, 5000])
    def test_inverse_round_trip_from_ic50(self, ic50: float):
        """ic50 → score → ic50 should return the original value."""
        score = self._ic50_to_binding_score(ic50)
        recovered_ic50 = binding_score_to_ic50(score)
        # Allow 0.1% tolerance for floating-point rounding
        assert abs(recovered_ic50 - ic50) / ic50 < 1e-3, (
            f"Round-trip from IC50 failed: {ic50} → score={score:.4f} → IC50={recovered_ic50:.2f}"
        )

    # -- Consistency with classify_binding --
    def test_classify_binding_consistent_with_score(self):
        """binding_score_to_ic50 → classify_binding should match the score tier."""
        # score=0.9 → IC50 near 50 → should be moderate_binder (since IC50 >= 50)
        ic50 = binding_score_to_ic50(0.9)
        cls = classify_binding(ic50)
        assert cls in ("strong_binder", "moderate_binder"), (
            f"Score 0.9 → IC50={ic50:.1f} classified as {cls}"
        )

        # score=0.0 → IC50 very high → should be non_binder
        ic50_low = binding_score_to_ic50(0.0)
        assert classify_binding(ic50_low) == "non_binder"


# ═══════════════════════════════════════════════════════════════════════════
# 4. Basic scanning functions return reasonable results
# ═══════════════════════════════════════════════════════════════════════════

class TestScanningFunctions:
    """Test that predict_mhc_i_binding, predict_mhc_ii_binding,
    predict_t_cell_epitopes, and predict_b_cell_epitopes return
    structurally reasonable results."""

    @pytest.fixture(autouse=True)
    def _clear_prediction_cache(self):
        """Ensure a fresh cache for each test."""
        clear_cache()
        yield
        clear_cache()

    # -- predict_mhc_i_binding --

    def test_mhc_i_returns_list_of_binding_results(self):
        results = predict_mhc_i_binding(SHORT_PROTEIN)
        assert isinstance(results, list)
        if results:
            assert all(isinstance(r, MHCBindingResult) for r in results)

    def test_mhc_i_binding_scores_in_range(self):
        results = predict_mhc_i_binding(SHORT_PROTEIN)
        for r in results:
            assert 0.0 <= r.binding_score <= 1.0, (
                f"binding_score {r.binding_score} out of [0,1]"
            )

    def test_mhc_i_ic50_positive(self):
        results = predict_mhc_i_binding(SHORT_PROTEIN)
        for r in results:
            if r.ic50_nm is not None:
                assert r.ic50_nm > 0, f"IC50 should be positive, got {r.ic50_nm}"

    def test_mhc_i_binding_class_valid(self):
        valid = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        results = predict_mhc_i_binding(SHORT_PROTEIN)
        for r in results:
            assert r.binding_class in valid, (
                f"Invalid binding_class: {r.binding_class}"
            )

    def test_mhc_i_peptide_length(self):
        """All MHC-I peptides should be 9-mers (default peptide_length)."""
        results = predict_mhc_i_binding(SHORT_PROTEIN)
        for r in results:
            assert len(r.peptide) == 9, (
                f"Peptide length {len(r.peptide)} != 9: {r.peptide}"
            )

    def test_mhc_i_position_indices_consistent(self):
        """end_position should equal start_position + len(peptide) - 1."""
        results = predict_mhc_i_binding(SHORT_PROTEIN)
        for r in results:
            expected_end = r.start_position + len(r.peptide) - 1
            assert r.end_position == expected_end, (
                f"Position mismatch: start={r.start_position}, end={r.end_position}, "
                f"peptide_len={len(r.peptide)}"
            )

    def test_mhc_i_empty_protein(self):
        assert predict_mhc_i_binding("") == []

    def test_mhc_i_known_binder_scores_high(self):
        """Known HLA-A*02:01 binder GILGFVFTL should score positively."""
        results = predict_mhc_i_binding(
            KNOWN_BINDER_9MER, alleles=["HLA-A*02:01"]
        )
        assert len(results) > 0, "No results for known binder"
        best = max(results, key=lambda r: r.binding_score)
        assert best.binding_score > 0.0, (
            f"Known binder scored {best.binding_score}, expected > 0"
        )

    # -- predict_mhc_ii_binding --

    def test_mhc_ii_returns_list_of_binding_results(self):
        results = predict_mhc_ii_binding(SHORT_PROTEIN)
        assert isinstance(results, list)
        if results:
            assert all(isinstance(r, MHCBindingResult) for r in results)

    def test_mhc_ii_binding_scores_in_range(self):
        results = predict_mhc_ii_binding(SHORT_PROTEIN)
        for r in results:
            assert 0.0 <= r.binding_score <= 1.0, (
                f"binding_score {r.binding_score} out of [0,1]"
            )

    def test_mhc_ii_binding_class_valid(self):
        valid = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        results = predict_mhc_ii_binding(SHORT_PROTEIN)
        for r in results:
            assert r.binding_class in valid, (
                f"Invalid binding_class: {r.binding_class}"
            )

    def test_mhc_ii_peptide_length(self):
        """Default MHC-II peptides should be 15-mers."""
        results = predict_mhc_ii_binding(SHORT_PROTEIN)
        for r in results:
            assert len(r.peptide) == 15, (
                f"Peptide length {len(r.peptide)} != 15: {r.peptide}"
            )

    def test_mhc_ii_empty_protein(self):
        assert predict_mhc_ii_binding("") == []

    def test_mhc_ii_short_protein_returns_empty(self):
        """Protein shorter than peptide_length should yield no results."""
        assert predict_mhc_ii_binding("ACDE") == []

    # -- predict_t_cell_epitopes --

    def test_t_cell_epitopes_returns_list(self):
        epitopes = predict_t_cell_epitopes(SHORT_PROTEIN)
        assert isinstance(epitopes, list)

    def test_t_cell_epitope_dict_keys(self):
        """Each T-cell epitope dict should contain the required keys."""
        epitopes = predict_t_cell_epitopes(SHORT_PROTEIN)
        required_keys = {"start", "end", "peptide", "score", "allele", "binding_class"}
        for ep in epitopes:
            assert required_keys.issubset(ep.keys()), (
                f"Missing keys: {required_keys - set(ep.keys())}"
            )

    def test_t_cell_epitope_scores_in_range(self):
        epitopes = predict_t_cell_epitopes(SHORT_PROTEIN)
        for ep in epitopes:
            assert 0.0 <= ep["score"] <= 1.0, (
                f"T-cell epitope score {ep['score']} out of [0,1]"
            )

    def test_t_cell_epitope_binding_class_valid(self):
        valid = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        epitopes = predict_t_cell_epitopes(SHORT_PROTEIN)
        for ep in epitopes:
            assert ep["binding_class"] in valid, (
                f"Invalid binding_class: {ep['binding_class']}"
            )

    def test_t_cell_epitopes_sorted_by_score(self):
        """Epitopes should be sorted by score descending."""
        epitopes = predict_t_cell_epitopes(HBB_HUMAN)
        if len(epitopes) >= 2:
            for i in range(len(epitopes) - 1):
                assert epitopes[i]["score"] >= epitopes[i + 1]["score"], (
                    f"Epitopes not sorted: {epitopes[i]['score']} < {epitopes[i+1]['score']}"
                )

    # -- predict_b_cell_epitopes --

    def test_b_cell_epitopes_returns_list(self):
        epitopes = predict_b_cell_epitopes(SHORT_PROTEIN)
        assert isinstance(epitopes, list)

    def test_b_cell_epitope_dict_keys(self):
        """Each B-cell epitope dict should contain at least start, end, peptide, score."""
        epitopes = predict_b_cell_epitopes(HBB_HUMAN)
        required_keys = {"start", "end", "peptide", "score"}
        for ep in epitopes:
            assert required_keys.issubset(ep.keys()), (
                f"Missing keys: {required_keys - set(ep.keys())}"
            )

    def test_b_cell_epitope_scores_in_range(self):
        epitopes = predict_b_cell_epitopes(HBB_HUMAN)
        for ep in epitopes:
            assert 0.0 <= ep["score"] <= 1.0, (
                f"B-cell epitope score {ep['score']} out of [0,1]"
            )

    # -- score_peptide_pssm edge cases --

    def test_score_peptide_wrong_length_returns_zero(self):
        """Peptide length != PSSM length should return 0.0."""
        score = score_peptide_pssm("GI", "HLA-A*02:01")
        assert score == 0.0

    def test_score_peptide_unknown_allele_returns_zero(self):
        """Unknown allele string should return 0.0."""
        score = score_peptide_pssm("GILGFVFTL", "HLA-FAKE*99:99")
        assert score == 0.0

    def test_score_peptide_known_binder_positive(self):
        """Known binder should score > 0."""
        score = score_peptide_pssm(KNOWN_BINDER_9MER, "HLA-A*02:01")
        assert score > 0.0

    def test_score_peptide_result_in_range(self):
        """score_peptide_pssm should always return a value in [0, 1]."""
        score = score_peptide_pssm(KNOWN_BINDER_9MER, "HLA-A*02:01")
        assert 0.0 <= score <= 1.0
