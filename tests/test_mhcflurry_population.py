"""Tests for biocompiler.mhcflurry_population — expanded population coverage for MHCflurry.

Covers:
  1. Data integrity of EXPANDED_POPULATION_COVERAGE
  2. get_allele_frequency lookups
  3. compute_population_coverage calculations
  4. find_coverage_optimizing_alleles greedy optimiser
  5. ALLELE_CLASSIFICATION consistency
  6. SUPPORTED_MHCFLURRY_ALLELES validity
  7. Backward-compatibility with original POPULATION_COVERAGE
"""

from __future__ import annotations

import re

import pytest

from biocompiler.immunogenicity import DEFAULT_MHC_I_ALLELES, POPULATION_COVERAGE

try:
    from biocompiler.mhcflurry_population import (
        ALLELE_CLASSIFICATION,
        EXPANDED_POPULATION_COVERAGE,
        POPULATION_GROUPS,
        SUPPORTED_MHCFLURRY_ALLELES,
        compute_population_coverage,
        find_coverage_optimizing_alleles,
        get_allele_frequency,
    )
except ImportError as exc:
    pytest.skip(f"mhcflurry_population not available: {exc}", allow_module_level=True)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Data integrity
# ═══════════════════════════════════════════════════════════════════════════


class TestDataIntegrity:
    """Validate structural and range invariants of EXPANDED_POPULATION_COVERAGE."""

    def test_at_least_50_alleles(self):
        assert len(EXPANDED_POPULATION_COVERAGE) >= 50

    def test_every_allele_has_all_population_groups(self):
        groups = set(POPULATION_GROUPS)
        for allele, freqs in EXPANDED_POPULATION_COVERAGE.items():
            absent = groups - set(freqs.keys())
            assert not absent, f"{allele} missing groups: {absent}"

    def test_all_frequencies_in_valid_range(self):
        for allele, freqs in EXPANDED_POPULATION_COVERAGE.items():
            for pop, freq in freqs.items():
                assert 0.0 <= freq <= 100.0, f"{allele}/{pop}={freq} outside [0,100]"

    def test_no_allele_has_zero_for_all_populations(self):
        all_zero = [a for a, f in EXPANDED_POPULATION_COVERAGE.items() if all(v == 0.0 for v in f.values())]
        assert not all_zero, f"Zero-frequency alleles: {all_zero}"

    def test_frequency_values_are_numeric(self):
        for allele, freqs in EXPANDED_POPULATION_COVERAGE.items():
            for pop, freq in freqs.items():
                assert isinstance(freq, (int, float)), f"{allele}/{pop}: {freq!r}"

    def test_population_groups_non_empty(self):
        assert len(POPULATION_GROUPS) > 0


# ═══════════════════════════════════════════════════════════════════════════
# 2. get_allele_frequency
# ═══════════════════════════════════════════════════════════════════════════


class TestGetAlleleFrequency:

    def test_known_allele_returns_correct_frequency(self):
        allele = next(iter(EXPANDED_POPULATION_COVERAGE))
        pop = POPULATION_GROUPS[0]
        expected = EXPANDED_POPULATION_COVERAGE[allele][pop]
        assert get_allele_frequency(allele, pop) == pytest.approx(expected)

    def test_unknown_allele_returns_zero(self):
        assert get_allele_frequency("HLA-FAKE*99:99", POPULATION_GROUPS[0]) == 0.0

    def test_unknown_population_returns_zero(self):
        allele = next(iter(EXPANDED_POPULATION_COVERAGE))
        assert get_allele_frequency(allele, "NonExistentPopulation") == 0.0

    def test_returns_float(self):
        allele = next(iter(EXPANDED_POPULATION_COVERAGE))
        assert isinstance(get_allele_frequency(allele, POPULATION_GROUPS[0]), float)

    def test_hla_a0201_caucasian(self):
        """HLA-A*02:01 Caucasian = 28.0 (original data)."""
        assert get_allele_frequency("HLA-A*02:01", "Caucasian") == pytest.approx(28.0)


# ═══════════════════════════════════════════════════════════════════════════
# 3. compute_population_coverage
# ═══════════════════════════════════════════════════════════════════════════


class TestComputePopulationCoverage:

    def test_single_common_allele_high_coverage(self):
        cov = compute_population_coverage(["HLA-A*02:01"], population="Caucasian")
        assert cov > 0.20, f"Single common allele coverage={cov:.3f} too low"

    def test_multiple_alleles_higher_coverage(self):
        single = compute_population_coverage(["HLA-A*02:01"], population="Caucasian")
        multi = compute_population_coverage(
            ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"], population="Caucasian"
        )
        assert multi >= single

    def test_global_population_reasonable(self):
        cov = compute_population_coverage(DEFAULT_MHC_I_ALLELES, population="global")
        assert 0.0 < cov <= 1.0

    def test_empty_allele_list_returns_zero(self):
        assert compute_population_coverage([], population="Caucasian") == 0.0

    def test_coverage_bounded_zero_to_one(self):
        all_alleles = list(EXPANDED_POPULATION_COVERAGE.keys())
        for pop in POPULATION_GROUPS:
            cov = compute_population_coverage(all_alleles, population=pop)
            assert 0.0 <= cov <= 1.0, f"{pop} coverage={cov}"

    def test_coverage_increases_with_more_alleles(self):
        cumulative = []
        running: list[str] = []
        for allele in DEFAULT_MHC_I_ALLELES:
            running.append(allele)
            cumulative.append(compute_population_coverage(running, population="Caucasian"))
        for i in range(1, len(cumulative)):
            assert cumulative[i] >= cumulative[i - 1] - 1e-9


# ═══════════════════════════════════════════════════════════════════════════
# 4. find_coverage_optimizing_alleles
# ═══════════════════════════════════════════════════════════════════════════


class TestFindCoverageOptimizingAlleles:

    def test_n6_returns_six_alleles(self):
        assert len(find_coverage_optimizing_alleles(n=6, population="Caucasian")) == 6

    def test_first_allele_is_most_common(self):
        alleles = find_coverage_optimizing_alleles(n=6, population="Caucasian")
        first_cov = compute_population_coverage([alleles[0]], population="Caucasian")
        for a in EXPANDED_POPULATION_COVERAGE:
            c = compute_population_coverage([a], population="Caucasian")
            assert c <= first_cov + 1e-9, f"{a} ({c:.4f}) > {alleles[0]} ({first_cov:.4f})"

    def test_coverage_monotonically_increases(self):
        alleles = find_coverage_optimizing_alleles(n=8, population="global")
        prev = 0.0
        for i in range(1, len(alleles) + 1):
            cov = compute_population_coverage(alleles[:i], population="global")
            assert cov >= prev - 1e-9
            prev = cov

    def test_n1_returns_single_most_common(self):
        alleles = find_coverage_optimizing_alleles(n=1, population="Caucasian")
        assert len(alleles) == 1
        freq = get_allele_frequency(alleles[0], "Caucasian")
        max_freq = max(EXPANDED_POPULATION_COVERAGE[a]["Caucasian"] for a in EXPANDED_POPULATION_COVERAGE)
        assert freq == pytest.approx(max_freq)

    def test_no_duplicate_alleles(self):
        alleles = find_coverage_optimizing_alleles(n=10, population="Caucasian")
        assert len(alleles) == len(set(alleles))

    def test_all_returned_in_coverage_data(self):
        alleles = find_coverage_optimizing_alleles(n=10, population="global")
        for a in alleles:
            assert a in EXPANDED_POPULATION_COVERAGE, f"{a!r} not in coverage data"


# ═══════════════════════════════════════════════════════════════════════════
# 5. ALLELE_CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════════


class TestAlleleClassification:

    def test_each_allele_classified_as_i_or_ii(self):
        invalid = {a: c for a, c in ALLELE_CLASSIFICATION.items() if c not in ("I", "II")}
        assert not invalid, f"Invalid classifications: {invalid}"

    def test_mhc_i_default_alleles_are_class_i(self):
        for allele in DEFAULT_MHC_I_ALLELES:
            if allele in ALLELE_CLASSIFICATION:
                assert ALLELE_CLASSIFICATION[allele] == "I", f"{allele} != 'I'"

    def test_mhc_ii_alleles_are_class_ii(self):
        for allele in ("HLA-DRB1*01:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01"):
            if allele in ALLELE_CLASSIFICATION:
                assert ALLELE_CLASSIFICATION[allele] == "II", f"{allele} != 'II'"

    def test_hla_a_b_c_are_class_i(self):
        for allele, cls in ALLELE_CLASSIFICATION.items():
            if allele.startswith(("HLA-A", "HLA-B", "HLA-C")):
                assert cls == "I", f"{allele} classified as {cls}, not 'I'"

    def test_hla_drb1_are_class_ii(self):
        for allele, cls in ALLELE_CLASSIFICATION.items():
            if allele.startswith("HLA-DRB1"):
                assert cls == "II", f"{allele} classified as {cls}, not 'II'"

    def test_classified_alleles_in_coverage_data(self):
        missing = set(ALLELE_CLASSIFICATION) - set(EXPANDED_POPULATION_COVERAGE)
        assert not missing, f"Classified alleles not in coverage: {sorted(missing)}"


# ═══════════════════════════════════════════════════════════════════════════
# 6. SUPPORTED_MHCFLURRY_ALLELES
# ═══════════════════════════════════════════════════════════════════════════


class TestSupportedMhcflurryAlleles:

    def test_list_is_non_empty(self):
        assert len(SUPPORTED_MHCFLURRY_ALLELES) > 0

    def test_contains_default_mhc_i_alleles(self):
        supported = set(SUPPORTED_MHCFLURRY_ALLELES)
        for allele in DEFAULT_MHC_I_ALLELES:
            assert allele in supported, f"Missing default allele {allele}"

    def test_all_match_hla_naming_pattern(self):
        pat = re.compile(r"^HLA-[A-Z]+\*\d{2}:\d{2}$")
        invalid = [a for a in SUPPORTED_MHCFLURRY_ALLELES if not pat.match(a)]
        assert not invalid, f"Non-HLA-pattern alleles: {invalid[:10]}"

    def test_no_duplicates(self):
        assert len(SUPPORTED_MHCFLURRY_ALLELES) == len(set(SUPPORTED_MHCFLURRY_ALLELES))

    def test_supported_in_coverage_data(self):
        missing = set(SUPPORTED_MHCFLURRY_ALLELES) - set(EXPANDED_POPULATION_COVERAGE)
        assert not missing, f"Supported alleles not in coverage: {sorted(missing)[:10]}"


# ═══════════════════════════════════════════════════════════════════════════
# 7. Comparison with original POPULATION_COVERAGE
# ═══════════════════════════════════════════════════════════════════════════


class TestBackwardCompatibility:

    def test_original_alleles_still_present(self):
        missing = set(POPULATION_COVERAGE) - set(EXPANDED_POPULATION_COVERAGE)
        assert not missing, f"Original alleles missing: {sorted(missing)}"

    def test_original_frequencies_match(self):
        mismatches = []
        for allele, orig_freqs in POPULATION_COVERAGE.items():
            exp_freqs = EXPANDED_POPULATION_COVERAGE[allele]
            for pop, orig in orig_freqs.items():
                exp = exp_freqs.get(pop)
                if exp is None:
                    mismatches.append(f"{allele}/{pop}: MISSING in expanded")
                elif exp != pytest.approx(orig):
                    mismatches.append(f"{allele}/{pop}: orig={orig} exp={exp}")
        assert not mismatches, "Frequency mismatches:\n" + "\n".join(mismatches[:10])

    def test_expanded_is_superset(self):
        assert len(EXPANDED_POPULATION_COVERAGE) > len(POPULATION_COVERAGE)

    def test_expanded_has_original_population_keys(self):
        orig_pops = {"Caucasian", "Asian", "African", "Hispanic"}
        for allele in POPULATION_COVERAGE:
            missing = orig_pops - set(EXPANDED_POPULATION_COVERAGE[allele])
            assert not missing, f"{allele}: missing original pops {missing}"

    def test_original_coverage_recomputed_correctly(self):
        """Cross-check: expanded module coverage matches direct computation from original."""
        alleles = list(POPULATION_COVERAGE.keys())[:3]
        for pop in ("Caucasian", "Asian"):
            coverage = compute_population_coverage(alleles, population=pop)
            # Direct independent calculation from original data
            product = 1.0
            for a in alleles:
                product *= (1.0 - POPULATION_COVERAGE[a].get(pop, 0.0) / 100.0)
            expected = 1.0 - product
            assert coverage == pytest.approx(expected, abs=0.01), (
                f"{pop}: module={coverage:.4f}, direct={expected:.4f}"
            )
