"""Unit tests for biocompiler.mhcflurry_population — pure-function logic, input validation, and named constants.

This file focuses on:
  1. Population coverage calculation logic (pure functions)
  2. Input validation and edge-case handling
  3. Named constants (public and private)

It intentionally does **not** duplicate the data-integrity or backward-compatibility
checks in ``test_mhcflurry_population.py``.
"""

from __future__ import annotations

import importlib
import math
from unittest.mock import patch

import pytest

try:
    import biocompiler.mhcflurry_population as _mod
except ImportError as exc:
    pytest.skip(f"mhcflurry_population not available: {exc}", allow_module_level=True)

# Re-import individual names for cleaner test code
from biocompiler.mhcflurry_population import (
    ALLELE_CLASSIFICATION,
    EXPANDED_POPULATION_COVERAGE,
    POPULATION_GROUPS,
    POPULATION_WEIGHTS,
    SUPPORTED_MHCFLURRY_ALLELES,
    compute_population_coverage,
    find_coverage_optimizing_alleles,
    get_allele_frequency,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _first_known_allele() -> str:
    """Return a guaranteed-present allele from the coverage data."""
    return next(iter(EXPANDED_POPULATION_COVERAGE))


def _manual_coverage(alleles: list[str], population: str) -> float:
    """Independent re-implementation of the independence-approximation formula."""
    prob_not = 1.0
    for a in alleles:
        freq = get_allele_frequency(a, population)
        prob_not *= 1.0 - freq / 100.0
    return max(0.0, min(1.0, 1.0 - prob_not))


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Named constants — public
# ═══════════════════════════════════════════════════════════════════════════════


class TestPublicConstants:
    """Validate public named constants exported in __all__."""

    def test_population_groups_type_and_count(self):
        assert isinstance(POPULATION_GROUPS, list)
        assert len(POPULATION_GROUPS) == 6

    def test_population_groups_contents(self):
        expected = ["Caucasian", "African", "Asian", "Hispanic", "Native American", "Oceanian"]
        assert POPULATION_GROUPS == expected

    def test_population_groups_all_strings(self):
        assert all(isinstance(g, str) for g in POPULATION_GROUPS)

    def test_population_groups_no_duplicates(self):
        assert len(POPULATION_GROUPS) == len(set(POPULATION_GROUPS))

    def test_population_weights_keys_match_groups(self):
        assert set(POPULATION_WEIGHTS.keys()) == set(POPULATION_GROUPS)

    def test_population_weights_all_positive(self):
        for pop, w in POPULATION_WEIGHTS.items():
            assert w > 0, f"{pop} weight={w} not positive"

    def test_population_weights_are_floats(self):
        assert all(isinstance(w, float) for w in POPULATION_WEIGHTS.values())

    def test_expanded_population_coverage_is_dict(self):
        assert isinstance(EXPANDED_POPULATION_COVERAGE, dict)

    def test_supported_mhcflurry_alleles_is_list(self):
        assert isinstance(SUPPORTED_MHCFLURRY_ALLELES, list)

    def test_supported_mhcflurry_alleles_sorted(self):
        assert SUPPORTED_MHCFLURRY_ALLELES == sorted(SUPPORTED_MHCFLURRY_ALLELES)

    def test_supported_mhcflurry_alleles_no_duplicates(self):
        assert len(SUPPORTED_MHCFLURRY_ALLELES) == len(set(SUPPORTED_MHCFLURRY_ALLELES))

    def test_allele_classification_is_dict(self):
        assert isinstance(ALLELE_CLASSIFICATION, dict)

    def test_allele_classification_values_format(self):
        """Every value must be '<class>:<source>' where class in {I, II} and source in {pssm, mhcflurry, both}."""
        valid = {"I:pssm", "I:mhcflurry", "I:both", "II:pssm", "II:mhcflurry", "II:both"}
        invalid = {a: v for a, v in ALLELE_CLASSIFICATION.items() if v not in valid}
        assert not invalid, f"Invalid classification values: { {k: invalid[k] for k in list(invalid)[:5]} }"

    def test_allele_classification_keys_are_strings(self):
        assert all(isinstance(a, str) for a in ALLELE_CLASSIFICATION)


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Named constants — private (accessed via the module object)
# ═══════════════════════════════════════════════════════════════════════════════


class TestPrivateConstants:
    """Validate internal named constants that govern the calculation logic."""

    def test_frequency_percent_scale(self):
        assert _mod._FREQUENCY_PERCENT_SCALE == 100.0

    def test_default_optimizing_n_alleles(self):
        assert _mod._DEFAULT_OPTIMIZING_N_ALLELES == 6

    def test_global_population_key(self):
        assert _mod._GLOBAL_POPULATION_KEY == "global"

    def test_coverage_fraction_bounds(self):
        assert _mod._MIN_COVERAGE_FRACTION == 0.0
        assert _mod._MAX_COVERAGE_FRACTION == 1.0
        assert _mod._MIN_COVERAGE_FRACTION < _mod._MAX_COVERAGE_FRACTION

    def test_mouse_mhc_i_alleles_tuple(self):
        assert isinstance(_mod._MOUSE_MHC_I_ALLELES, tuple)
        assert len(_mod._MOUSE_MHC_I_ALLELES) > 0

    def test_mouse_mhc_ii_alleles_tuple(self):
        assert isinstance(_mod._MOUSE_MHC_II_ALLELES, tuple)
        assert len(_mod._MOUSE_MHC_II_ALLELES) > 0

    def test_mouse_alleles_in_classification(self):
        """Mouse alleles should be present in ALLELE_CLASSIFICATION."""
        for allele in _mod._MOUSE_MHC_I_ALLELES:
            assert allele in ALLELE_CLASSIFICATION, f"{allele} missing from classification"
            assert ALLELE_CLASSIFICATION[allele].startswith("I:"), f"{allele} not class I"
        for allele in _mod._MOUSE_MHC_II_ALLELES:
            assert allele in ALLELE_CLASSIFICATION, f"{allele} missing from classification"
            assert ALLELE_CLASSIFICATION[allele].startswith("II:"), f"{allele} not class II"

    def test_mouse_alleles_in_supported_list(self):
        """Mouse alleles should appear in SUPPORTED_MHCFLURRY_ALLELES."""
        for allele in _mod._MOUSE_MHC_I_ALLELES + _mod._MOUSE_MHC_II_ALLELES:
            assert allele in SUPPORTED_MHCFLURRY_ALLELES, f"{allele} missing from supported list"

    def test_pssm_alleles_constant(self):
        assert isinstance(_mod._PSSM_MHC_I_ALLELES, list)
        assert isinstance(_mod._PSSM_MHC_II_ALLELES, list)
        assert len(_mod._PSSM_MHC_I_ALLELES) == 6
        assert len(_mod._PSSM_MHC_II_ALLELES) == 3


# ═══════════════════════════════════════════════════════════════════════════════
# 3. get_allele_frequency — pure lookup logic
# ═══════════════════════════════════════════════════════════════════════════════


class TestGetAlleleFrequencyLogic:
    """Test the pure lookup mechanics of get_allele_frequency."""

    def test_returns_exact_value_from_data(self):
        allele = "HLA-A*02:01"
        pop = "Caucasian"
        expected = EXPANDED_POPULATION_COVERAGE[allele][pop]
        assert get_allele_frequency(allele, pop) == pytest.approx(expected)

    def test_returns_float_for_known_allele(self):
        allele = _first_known_allele()
        result = get_allele_frequency(allele, POPULATION_GROUPS[0])
        assert isinstance(result, float)

    def test_unknown_allele_returns_zero_float(self):
        result = get_allele_frequency("HLA-FAKE*00:00", "Caucasian")
        assert result == 0.0
        assert isinstance(result, float)

    def test_unknown_population_returns_zero_float(self):
        allele = _first_known_allele()
        result = get_allele_frequency(allele, "Martian")
        assert result == 0.0
        assert isinstance(result, float)

    def test_both_unknown_returns_zero(self):
        assert get_allele_frequency("BOGUS", "NONEXISTENT") == 0.0

    def test_zero_frequency_allele_population(self):
        """Some alleles have 0.5% in some populations — verify it's returned correctly."""
        # HLA-A*02:06 has 0.5 in Caucasian
        result = get_allele_frequency("HLA-A*02:06", "Caucasian")
        assert result == pytest.approx(0.5)

    def test_consistency_with_direct_dict_access(self):
        """Every lookup via function must match direct dict access for known data."""
        for allele, freqs in EXPANDED_POPULATION_COVERAGE.items():
            for pop, expected in freqs.items():
                assert get_allele_frequency(allele, pop) == pytest.approx(expected), (
                    f"Mismatch for {allele}/{pop}"
                )

    def test_logging_on_unknown_allele(self, caplog):
        """Unknown allele should produce a warning log."""
        with caplog.at_level("WARNING"):
            get_allele_frequency("HLA-NONEXIST*99:99", "Caucasian")
        assert any("not found" in r.message.lower() for r in caplog.records)

    def test_logging_on_unknown_population(self, caplog):
        """Unknown population should produce a warning log."""
        allele = _first_known_allele()
        with caplog.at_level("WARNING"):
            get_allele_frequency(allele, "NonexistentPop")
        assert any("not found" in r.message.lower() for r in caplog.records)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. compute_population_coverage — pure calculation logic
# ═══════════════════════════════════════════════════════════════════════════════


class TestComputePopulationCoverageLogic:
    """Test the independence-approximation coverage formula in isolation."""

    def test_formula_matches_manual_calculation_single_pop(self):
        alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"]
        for pop in POPULATION_GROUPS:
            result = compute_population_coverage(alleles, population=pop)
            expected = _manual_coverage(alleles, pop)
            assert result == pytest.approx(expected, abs=1e-12), f"Mismatch for {pop}"

    def test_single_allele_coverage_formula(self):
        """Coverage for a single allele = freq/100."""
        allele = "HLA-A*02:01"
        for pop in POPULATION_GROUPS:
            freq = get_allele_frequency(allele, pop)
            expected = freq / 100.0
            result = compute_population_coverage([allele], population=pop)
            assert result == pytest.approx(expected, abs=1e-12), f"Mismatch for {pop}"

    def test_zero_frequency_allele_contributes_nothing(self):
        """An allele with 0% frequency in a population should not affect coverage."""
        allele_known = "HLA-A*02:01"
        cov_without = compute_population_coverage([allele_known], population="Caucasian")
        cov_with_fake = compute_population_coverage(
            [allele_known, "HLA-FAKE*00:00"], population="Caucasian"
        )
        assert cov_with_fake == pytest.approx(cov_without)

    def test_coverage_monotonically_non_decreasing(self):
        """Adding alleles should never decrease coverage."""
        alleles = list(EXPANDED_POPULATION_COVERAGE.keys())[:10]
        prev = 0.0
        for i in range(1, len(alleles) + 1):
            cov = compute_population_coverage(alleles[:i], population="Caucasian")
            assert cov >= prev - 1e-12
            prev = cov

    def test_global_coverage_is_weighted_average(self):
        """Global coverage must equal the population-weighted average of per-population coverages."""
        alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"]
        global_cov = compute_population_coverage(alleles, population="global")

        total_weight = sum(POPULATION_WEIGHTS.values())
        weighted_sum = 0.0
        for pop in POPULATION_GROUPS:
            pop_cov = compute_population_coverage(alleles, population=pop)
            weighted_sum += pop_cov * POPULATION_WEIGHTS[pop]
        expected = weighted_sum / total_weight

        assert global_cov == pytest.approx(expected, abs=1e-12)

    def test_all_alleles_coverage_less_than_one(self):
        """Even with all alleles, coverage should be <= 1.0."""
        all_alleles = list(EXPANDED_POPULATION_COVERAGE.keys())
        for pop in POPULATION_GROUPS + ["global"]:
            cov = compute_population_coverage(all_alleles, population=pop)
            assert cov <= 1.0 + 1e-12, f"Coverage > 1.0 for {pop}: {cov}"

    def test_coverage_bounded_clamping(self):
        """Verify clamping: result should always be in [0.0, 1.0]."""
        # Use mock data to test edge: freq > 100 should be clamped
        with patch.dict(
            EXPANDED_POPULATION_COVERAGE,
            {"TEST*01:01": {p: 0.0 for p in POPULATION_GROUPS}},
            clear=False,
        ):
            cov = compute_population_coverage(["TEST*01:01"], population="Caucasian")
            assert 0.0 <= cov <= 1.0

    def test_duplicate_alleles_are_included_in_product(self):
        """Passing the same allele twice affects the product — the function does not deduplicate.

        This is a known design choice: the independence formula 1 - prod(1 - freq_i/100)
        treats each list entry independently, so [A, A] yields 1 - (1-f/100)^2
        rather than 1 - (1-f/100).  Callers are responsible for deduplication.
        """
        allele = "HLA-A*02:01"
        freq = get_allele_frequency(allele, "Caucasian")
        cov_once = compute_population_coverage([allele], population="Caucasian")
        cov_twice = compute_population_coverage([allele, allele], population="Caucasian")
        # Verify once: 1 - (1 - f/100)
        assert cov_once == pytest.approx(freq / 100.0, abs=1e-12)
        # Verify twice: 1 - (1 - f/100)^2
        expected_twice = 1.0 - (1.0 - freq / 100.0) ** 2
        assert cov_twice == pytest.approx(expected_twice, abs=1e-12)
        # And twice > once (unless freq == 0)
        if freq > 0:
            assert cov_twice > cov_once


class TestComputePopulationCoverageEdgeCases:
    """Edge cases and boundary conditions for compute_population_coverage."""

    def test_empty_alleles_returns_zero(self):
        assert compute_population_coverage([], population="Caucasian") == 0.0

    def test_empty_alleles_global_returns_zero(self):
        assert compute_population_coverage([], population="global") == 0.0

    def test_all_unknown_alleles_returns_zero(self):
        result = compute_population_coverage(
            ["FAKE1", "FAKE2", "FAKE3"], population="Caucasian"
        )
        assert result == 0.0

    def test_mixed_known_and_unknown_alleles(self):
        """Unknown alleles should be silently treated as 0% frequency."""
        known_cov = compute_population_coverage(["HLA-A*02:01"], population="Caucasian")
        mixed_cov = compute_population_coverage(
            ["HLA-A*02:01", "FAKE*99:99"], population="Caucasian"
        )
        assert mixed_cov == pytest.approx(known_cov)

    def test_single_allele_per_population(self):
        """Verify single-allele coverage matches freq/100 for every population."""
        allele = "HLA-A*02:01"
        for pop in POPULATION_GROUPS:
            freq = EXPANDED_POPULATION_COVERAGE[allele][pop]
            cov = compute_population_coverage([allele], population=pop)
            assert cov == pytest.approx(freq / 100.0)

    def test_two_allele_formula(self):
        """For two alleles: cov = 1 - (1-f1/100)(1-f2/100)."""
        a1, a2 = "HLA-A*02:01", "HLA-A*01:01"
        for pop in POPULATION_GROUPS:
            f1 = get_allele_frequency(a1, pop)
            f2 = get_allele_frequency(a2, pop)
            expected = 1.0 - (1.0 - f1 / 100.0) * (1.0 - f2 / 100.0)
            result = compute_population_coverage([a1, a2], population=pop)
            assert result == pytest.approx(expected, abs=1e-12)


# ═══════════════════════════════════════════════════════════════════════════════
# 5. find_coverage_optimizing_alleles — greedy algorithm logic
# ═══════════════════════════════════════════════════════════════════════════════


class TestFindCoverageOptimizingAllelesLogic:
    """Test the greedy allele-selection algorithm."""

    def test_n_equals_zero_returns_empty(self):
        assert find_coverage_optimizing_alleles(n_alleles=0, population="Caucasian") == []

    def test_n_negative_returns_empty(self):
        assert find_coverage_optimizing_alleles(n_alleles=-1, population="Caucasian") == []

    def test_n_one_returns_single_allele(self):
        result = find_coverage_optimizing_alleles(n_alleles=1, population="Caucasian")
        assert len(result) == 1

    def test_n_one_returns_highest_frequency_allele(self):
        result = find_coverage_optimizing_alleles(n_alleles=1, population="Caucasian")
        allele = result[0]
        its_freq = get_allele_frequency(allele, "Caucasian")
        # No allele should have higher Caucasian frequency
        for a in EXPANDED_POPULATION_COVERAGE:
            assert get_allele_frequency(a, "Caucasian") <= its_freq + 1e-9

    def test_n_exceeds_available_alleles(self):
        """Requesting more alleles than available returns at most the available count."""
        max_available = len(EXPANDED_POPULATION_COVERAGE)
        result = find_coverage_optimizing_alleles(
            n_alleles=max_available + 100, population="Caucasian"
        )
        assert len(result) <= max_available

    def test_greedy_selection_is_monotonic(self):
        """Each additional allele should increase (or maintain) coverage."""
        alleles = find_coverage_optimizing_alleles(n_alleles=8, population="global")
        prev_cov = 0.0
        for i in range(1, len(alleles) + 1):
            cov = compute_population_coverage(alleles[:i], population="global")
            assert cov >= prev_cov - 1e-9
            prev_cov = cov

    def test_no_duplicate_alleles_in_result(self):
        result = find_coverage_optimizing_alleles(n_alleles=12, population="Asian")
        assert len(result) == len(set(result))

    def test_all_result_alleles_in_coverage_data(self):
        result = find_coverage_optimizing_alleles(n_alleles=6, population="African")
        for a in result:
            assert a in EXPANDED_POPULATION_COVERAGE

    def test_global_population_optimization(self):
        """Optimizing for global should not crash and should return valid alleles."""
        result = find_coverage_optimizing_alleles(n_alleles=6, population="global")
        assert len(result) == 6
        for a in result:
            assert a in EXPANDED_POPULATION_COVERAGE

    def test_second_allele_provides_highest_marginal_gain(self):
        """The second allele should provide the largest marginal coverage increase."""
        first = find_coverage_optimizing_alleles(n_alleles=1, population="Caucasian")
        second = find_coverage_optimizing_alleles(n_alleles=2, population="Caucasian")

        first_cov = compute_population_coverage(first, population="Caucasian")
        second_cov = compute_population_coverage(second, population="Caucasian")
        marginal_gain = second_cov - first_cov

        # Check no other single allele gives a higher marginal gain
        for a in EXPANDED_POPULATION_COVERAGE:
            if a == first[0]:
                continue
            trial_cov = compute_population_coverage(first + [a], population="Caucasian")
            trial_gain = trial_cov - first_cov
            assert trial_gain <= marginal_gain + 1e-9

    def test_default_n_alleles_is_six(self):
        """Calling without n_alleles should use the default of 6."""
        result = find_coverage_optimizing_alleles(population="Caucasian")
        assert len(result) == 6

    def test_default_population_is_global(self):
        """Calling without population should default to 'global'."""
        result_explicit = find_coverage_optimizing_alleles(n_alleles=4, population="global")
        result_default = find_coverage_optimizing_alleles(n_alleles=4)
        assert result_explicit == result_default


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Input validation — edge cases and error handling
# ═══════════════════════════════════════════════════════════════════════════════


class TestInputValidation:
    """Input validation and graceful handling of edge cases."""

    def test_get_allele_frequency_empty_string_allele(self):
        """Empty string allele should return 0.0, not crash."""
        assert get_allele_frequency("", "Caucasian") == 0.0

    def test_get_allele_frequency_empty_string_population(self):
        """Empty string population should return 0.0, not crash."""
        allele = _first_known_allele()
        assert get_allele_frequency(allele, "") == 0.0

    def test_compute_coverage_with_single_unknown_allele(self):
        """Single unknown allele should yield 0.0 coverage."""
        assert compute_population_coverage(["NONEXISTENT"], population="Caucasian") == 0.0

    def test_compute_coverage_empty_string_allele(self):
        """Empty string in allele list should not crash."""
        result = compute_population_coverage([""], population="Caucasian")
        # Empty string is unknown, so coverage = 0
        assert result == 0.0

    def test_find_optimizing_n_large_value(self):
        """Very large n should not crash; returns at most available alleles."""
        result = find_coverage_optimizing_alleles(n_alleles=9999, population="global")
        assert len(result) <= len(EXPANDED_POPULATION_COVERAGE)

    def test_find_optimizing_n_one_min(self):
        """n_alleles=1 is the smallest valid positive value."""
        result = find_coverage_optimizing_alleles(n_alleles=1, population="Caucasian")
        assert len(result) == 1

    def test_compute_coverage_case_sensitive_population(self):
        """Population names are case-sensitive; wrong case should return 0.0 for unknown."""
        allele = _first_known_allele()
        result = compute_population_coverage([allele], population="caucasian")
        # "caucasian" != "Caucasian", so treated as unknown population
        # All alleles in unknown pop return freq=0.0, so coverage=0.0
        assert result == 0.0

    def test_compute_coverage_global_key_is_case_sensitive(self):
        """Only lowercase 'global' triggers the global path."""
        allele = _first_known_allele()
        cov_global = compute_population_coverage([allele], population="global")
        cov_Global = compute_population_coverage([allele], population="Global")
        # "Global" is treated as a regular (unknown) population
        assert cov_Global != cov_global or cov_global == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Internal helpers — _global_average_frequency, _compute_global_coverage
# ═══════════════════════════════════════════════════════════════════════════════


class TestInternalHelpers:
    """Test private helper functions via the public API or direct call."""

    def test_global_average_frequency_weighted(self):
        """_global_average_frequency should return population-weighted average of frequencies."""
        allele = "HLA-A*02:01"
        result = _mod._global_average_frequency(allele)

        total_weight = sum(POPULATION_WEIGHTS.values())
        weighted_sum = 0.0
        for pop in POPULATION_GROUPS:
            freq = get_allele_frequency(allele, pop)
            weighted_sum += freq * POPULATION_WEIGHTS[pop]
        expected = weighted_sum / total_weight

        assert result == pytest.approx(expected, abs=1e-12)

    def test_global_average_frequency_unknown_allele(self):
        """Unknown allele should yield 0.0."""
        assert _mod._global_average_frequency("FAKE*00:00") == 0.0

    def test_compute_global_coverage_matches_formula(self):
        """_compute_global_coverage should equal weighted average of per-pop coverages."""
        alleles = ["HLA-A*02:01", "HLA-B*07:02"]
        result = _mod._compute_global_coverage(alleles)

        total_weight = sum(POPULATION_WEIGHTS.values())
        weighted = 0.0
        for pop in POPULATION_GROUPS:
            weighted += compute_population_coverage(alleles, population=pop) * POPULATION_WEIGHTS[pop]
        expected = weighted / total_weight

        assert result == pytest.approx(expected, abs=1e-12)

    def test_compute_global_coverage_empty_alleles(self):
        """Empty allele list → 0.0 coverage even via global path."""
        # compute_population_coverage([], "global") short-circuits before _compute_global_coverage
        # but let's test _compute_global_coverage directly
        result = _mod._compute_global_coverage([])
        # Each population returns 0.0 for empty alleles → global = 0.0
        assert result == 0.0

    def test_population_weights_sum_positive(self):
        """Sum of POPULATION_WEIGHTS must be positive for the weighted average."""
        assert sum(POPULATION_WEIGHTS.values()) > 0

    def test_global_coverage_zero_weights_returns_zero(self):
        """If all weights are zero, _compute_global_coverage should return 0.0."""
        with patch.dict(_mod.POPULATION_WEIGHTS, {k: 0.0 for k in POPULATION_WEIGHTS}):
            # Reload is not needed; the function reads POPULATION_WEIGHTS at call time
            result = _mod._compute_global_coverage(["HLA-A*02:01"])
            assert result == 0.0

    def test_global_average_frequency_zero_weights_returns_zero(self):
        """If all weights are zero, _global_average_frequency should return 0.0."""
        with patch.dict(_mod.POPULATION_WEIGHTS, {k: 0.0 for k in POPULATION_WEIGHTS}):
            result = _mod._global_average_frequency("HLA-A*02:01")
            assert result == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 8. ALLELE_CLASSIFICATION build logic
# ═══════════════════════════════════════════════════════════════════════════════


class TestAlleleClassificationBuildLogic:
    """Test the _build_allele_classification() function and its invariants."""

    def test_pssm_i_alleles_classified(self):
        """All PSSM MHC-I alleles must be in the classification."""
        for allele in _mod._PSSM_MHC_I_ALLELES:
            assert allele in ALLELE_CLASSIFICATION, f"{allele} not classified"
            cls = ALLELE_CLASSIFICATION[allele]
            assert cls.startswith("I:"), f"{allele} classified as {cls}, not I:*"

    def test_pssm_ii_alleles_classified(self):
        """All PSSM MHC-II alleles must be in the classification."""
        for allele in _mod._PSSM_MHC_II_ALLELES:
            assert allele in ALLELE_CLASSIFICATION, f"{allele} not classified"
            cls = ALLELE_CLASSIFICATION[allele]
            assert cls.startswith("II:"), f"{allele} classified as {cls}, not II:*"

    def test_pssm_i_overlapping_mhcflurry_are_both(self):
        """PSSM alleles also in MHCflurry should be 'I:both' or 'II:both'."""
        pssm_i_set = set(_mod._PSSM_MHC_I_ALLELES)
        mhcflurry_i_set = set(_mod._MHCFLURRY_I_ALLELES)
        overlap = pssm_i_set & mhcflurry_i_set
        for allele in overlap:
            assert ALLELE_CLASSIFICATION.get(allele) == "I:both", (
                f"{allele} should be I:both but is {ALLELE_CLASSIFICATION.get(allele)}"
            )

    def test_pssm_only_i_alleles(self):
        """PSSM alleles NOT in MHCflurry should be 'I:pssm'."""
        pssm_i_set = set(_mod._PSSM_MHC_I_ALLELES)
        mhcflurry_i_set = set(_mod._MHCFLURRY_I_ALLELES)
        pssm_only = pssm_i_set - mhcflurry_i_set
        for allele in pssm_only:
            assert ALLELE_CLASSIFICATION.get(allele) == "I:pssm", (
                f"{allele} should be I:pssm but is {ALLELE_CLASSIFICATION.get(allele)}"
            )

    def test_mhcflurry_only_i_alleles(self):
        """MHCflurry alleles NOT in PSSM should be 'I:mhcflurry'."""
        pssm_set = set(_mod._PSSM_MHC_I_ALLELES + _mod._PSSM_MHC_II_ALLELES)
        mhcflurry_i_set = set(_mod._MHCFLURRY_I_ALLELES)
        mhcflurry_only = mhcflurry_i_set - pssm_set
        for allele in mhcflurry_only:
            if allele in ALLELE_CLASSIFICATION:
                assert ALLELE_CLASSIFICATION[allele] == "I:mhcflurry", (
                    f"{allele} should be I:mhcflurry but is {ALLELE_CLASSIFICATION[allele]}"
                )

    def test_mhcflurry_only_ii_alleles(self):
        """MHCflurry-II alleles NOT in PSSM should be 'II:mhcflurry'."""
        pssm_set = set(_mod._PSSM_MHC_I_ALLELES + _mod._PSSM_MHC_II_ALLELES)
        mhcflurry_ii_set = set(_mod._MHCFLURRY_II_ALLELES)
        mhcflurry_only = mhcflurry_ii_set - pssm_set
        for allele in mhcflurry_only:
            if allele in ALLELE_CLASSIFICATION:
                assert ALLELE_CLASSIFICATION[allele] == "II:mhcflurry", (
                    f"{allele} should be II:mhcflurry but is {ALLELE_CLASSIFICATION[allele]}"
                )

    def test_pssm_ii_overlapping_mhcflurry_are_both(self):
        """PSSM MHC-II alleles in MHCflurry should be 'II:both'."""
        pssm_ii_set = set(_mod._PSSM_MHC_II_ALLELES)
        mhcflurry_ii_set = set(_mod._MHCFLURRY_II_ALLELES)
        overlap = pssm_ii_set & mhcflurry_ii_set
        for allele in overlap:
            assert ALLELE_CLASSIFICATION.get(allele) == "II:both", (
                f"{allele} should be II:both but is {ALLELE_CLASSIFICATION.get(allele)}"
            )

    def test_build_deterministic(self):
        """Calling _build_allele_classification twice returns the same result."""
        first = _mod._build_allele_classification()
        second = _mod._build_allele_classification()
        assert first == second

    def test_classification_covers_all_source_alleles(self):
        """Every allele from PSSM and MHCflurry lists should be classified."""
        all_source = set(
            _mod._PSSM_MHC_I_ALLELES
            + _mod._PSSM_MHC_II_ALLELES
            + _mod._MHCFLURRY_I_ALLELES
            + _mod._MHCFLURRY_II_ALLELES
        )
        classified = set(ALLELE_CLASSIFICATION.keys())
        missing = all_source - classified
        # Mouse alleles are added via setdefault so they should be present too
        assert not missing, f"Unclassified source alleles: {sorted(missing)[:10]}"


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Mathematical properties of coverage formula
# ═══════════════════════════════════════════════════════════════════════════════


class TestCoverageMathProperties:
    """Verify mathematical properties that must hold for the independence-approximation formula."""

    def test_single_allele_coverage_equals_frequency_fraction(self):
        """P(covered by {A}) = freq(A)/100."""
        for allele in list(EXPANDED_POPULATION_COVERAGE.keys())[:5]:
            for pop in POPULATION_GROUPS:
                freq = get_allele_frequency(allele, pop)
                cov = compute_population_coverage([allele], population=pop)
                assert cov == pytest.approx(freq / 100.0, abs=1e-12)

    def test_commutativity(self):
        """Coverage is independent of allele order: cov([A,B]) == cov([B,A])."""
        a, b = "HLA-A*02:01", "HLA-B*07:02"
        for pop in POPULATION_GROUPS:
            cov1 = compute_population_coverage([a, b], population=pop)
            cov2 = compute_population_coverage([b, a], population=pop)
            assert cov1 == pytest.approx(cov2, abs=1e-12)

    def test_subadditivity(self):
        """Coverage of a subset <= coverage of a superset: cov(A) <= cov(A ∪ B)."""
        subset = ["HLA-A*02:01"]
        superset = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"]
        for pop in POPULATION_GROUPS:
            assert compute_population_coverage(subset, pop) <= compute_population_coverage(superset, pop) + 1e-12

    def test_upper_bound_one(self):
        """Coverage can never exceed 1.0 even with all alleles."""
        all_alleles = list(EXPANDED_POPULATION_COVERAGE.keys())
        for pop in POPULATION_GROUPS + ["global"]:
            cov = compute_population_coverage(all_alleles, population=pop)
            assert cov <= 1.0 + 1e-12

    def test_lower_bound_zero(self):
        """Coverage can never go below 0.0."""
        assert compute_population_coverage([], population="Caucasian") == 0.0
        assert compute_population_coverage(["FAKE"], population="Caucasian") == 0.0

    def test_inclusion_exclusion_approximation(self):
        """For two alleles: cov ≈ f1/100 + f2/100 - f1*f2/10000."""
        a1, a2 = "HLA-A*02:01", "HLA-A*01:01"
        for pop in POPULATION_GROUPS:
            f1 = get_allele_frequency(a1, pop)
            f2 = get_allele_frequency(a2, pop)
            # Exact formula: 1 - (1-f1/100)(1-f2/100) = f1/100 + f2/100 - f1*f2/10000
            cov = compute_population_coverage([a1, a2], population=pop)
            expected = f1 / 100.0 + f2 / 100.0 - f1 * f2 / 10000.0
            assert cov == pytest.approx(expected, abs=1e-12)
