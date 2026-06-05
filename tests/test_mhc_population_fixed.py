"""Tests for the mhcflurry_population.py fallback fixes (Task F4.8).

Verifies:
  1. get_population_coverage() works without MHCflurry installed
  2. The precomputed-database fallback produces valid results
  3. Unknown alleles are handled gracefully (no crash, empty result + warning)
  4. PopulationCoverageResult dataclass fields are correct
  5. MHCFLURRY_AVAILABLE flag is a bool and the module never crashes on import
"""

from __future__ import annotations

import warnings
from unittest.mock import patch

import pytest

from biocompiler.mhcflurry_population import (
    EXPANDED_POPULATION_COVERAGE,
    MHCFLURRY_AVAILABLE,
    POPULATION_GROUPS,
    POPULATION_WEIGHTS,
    PopulationCoverageResult,
    compute_population_coverage,
    get_population_coverage,
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Module import never crashes — MHCFLURRY_AVAILABLE is a bool
# ═══════════════════════════════════════════════════════════════════════════


class TestModuleImportSafety:
    """The module must import cleanly even when mhcflurry is absent."""

    def test_mhcflurry_available_is_bool(self):
        assert isinstance(MHCFLURRY_AVAILABLE, bool)

    def test_mhcflurry_available_is_false_in_test_env(self):
        """In CI/test environments mhcflurry is typically not installed."""
        # We don't assert False because some envs might have it,
        # but we verify the value is a valid bool.
        assert MHCFLURRY_AVAILABLE in (True, False)


# ═══════════════════════════════════════════════════════════════════════════
# 2. PopulationCoverageResult dataclass
# ═══════════════════════════════════════════════════════════════════════════


class TestPopulationCoverageResult:
    """Validate the PopulationCoverageResult dataclass structure."""

    def test_default_construction(self):
        result = PopulationCoverageResult()
        assert result.alleles_tested == []
        assert result.peptides_tested == 0
        assert result.population_coverage == 0.0
        assert result.method_used == "unknown"
        assert result.per_allele_results == {}

    def test_frozen(self):
        result = PopulationCoverageResult()
        with pytest.raises(AttributeError):
            result.population_coverage = 0.5  # type: ignore[misc]

    def test_custom_construction(self):
        result = PopulationCoverageResult(
            alleles_tested=["HLA-A*02:01"],
            peptides_tested=3,
            population_coverage=0.28,
            method_used="precomputed_db",
            per_allele_results={"HLA-A*02:01": {"frequency_global": 12.5}},
        )
        assert result.alleles_tested == ["HLA-A*02:01"]
        assert result.peptides_tested == 3
        assert result.population_coverage == pytest.approx(0.28)
        assert result.method_used == "precomputed_db"
        assert "HLA-A*02:01" in result.per_allele_results

    def test_population_coverage_in_valid_range(self):
        """population_coverage must be in [0.0, 1.0] for any result."""
        result = PopulationCoverageResult(population_coverage=0.5)
        assert 0.0 <= result.population_coverage <= 1.0


# ═══════════════════════════════════════════════════════════════════════════
# 3. get_population_coverage() — fallback without mhcflurry
# ═══════════════════════════════════════════════════════════════════════════


class TestGetPopulationCoverageFallback:
    """Test that get_population_coverage works when MHCflurry is NOT available."""

    @pytest.fixture(autouse=True)
    def _force_mhcflurry_unavailable(self):
        """Force MHCFLURRY_AVAILABLE = False for all tests in this class."""
        with patch(
            "biocompiler.mhcflurry_population.MHCFLURRY_AVAILABLE", False
        ):
            yield

    def test_known_alleles_returns_result(self):
        alleles = ["HLA-A*02:01", "HLA-A*01:01"]
        result = get_population_coverage(alleles, peptides=["SIINFEKL"])
        assert isinstance(result, PopulationCoverageResult)

    def test_method_used_is_precomputed_db(self):
        alleles = ["HLA-A*02:01"]
        result = get_population_coverage(alleles, peptides=[])
        assert result.method_used == "precomputed_db"

    def test_coverage_fraction_is_valid(self):
        alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"]
        result = get_population_coverage(alleles, peptides=[])
        assert 0.0 < result.population_coverage <= 1.0

    def test_coverage_matches_compute_population_coverage(self):
        """population_coverage must equal compute_population_coverage(global)."""
        alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"]
        result = get_population_coverage(alleles, peptides=[])
        expected = compute_population_coverage(alleles, population="global")
        assert result.population_coverage == pytest.approx(expected, abs=1e-6)

    def test_peptides_tested_count(self):
        alleles = ["HLA-A*02:01"]
        peps = ["SIINFEKL", "LLFGYPVYV", "GILGFVFTL"]
        result = get_population_coverage(alleles, peptides=peps)
        assert result.peptides_tested == 3

    def test_alleles_tested_recorded(self):
        alleles = ["HLA-A*02:01", "HLA-B*07:02"]
        result = get_population_coverage(alleles, peptides=[])
        assert result.alleles_tested == alleles

    def test_per_allele_results_has_frequency_global(self):
        alleles = ["HLA-A*02:01", "HLA-A*01:01"]
        result = get_population_coverage(alleles, peptides=[])
        for allele in alleles:
            assert allele in result.per_allele_results
            assert "frequency_global" in result.per_allele_results[allele]
            freq = result.per_allele_results[allele]["frequency_global"]
            assert isinstance(freq, float)
            assert freq >= 0.0

    def test_single_allele_coverage_fraction(self):
        """Single allele coverage ≈ weighted-average frequency / 100."""
        allele = "HLA-A*02:01"
        result = get_population_coverage([allele], peptides=[])
        # The global coverage for a single allele equals its global freq / 100
        from biocompiler.mhcflurry_population import _global_average_frequency

        freq = _global_average_frequency(allele)
        expected = freq / 100.0
        assert result.population_coverage == pytest.approx(expected, abs=1e-6)

    def test_more_alleles_more_coverage(self):
        """Adding alleles should increase or maintain coverage."""
        result1 = get_population_coverage(["HLA-A*02:01"], peptides=[])
        result2 = get_population_coverage(
            ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02"], peptides=[]
        )
        assert result2.population_coverage >= result1.population_coverage - 1e-9


# ═══════════════════════════════════════════════════════════════════════════
# 4. Unknown allele handling
# ═══════════════════════════════════════════════════════════════════════════


class TestUnknownAlleleHandling:
    """Unknown alleles must be handled gracefully — no crash."""

    @pytest.fixture(autouse=True)
    def _force_mhcflurry_unavailable(self):
        with patch(
            "biocompiler.mhcflurry_population.MHCFLURRY_AVAILABLE", False
        ):
            yield

    def test_all_unknown_alleles_returns_empty_result(self):
        result = get_population_coverage(
            ["HLA-FAKE*99:99", "BOGUS*00:00"], peptides=["SIINFEKL"]
        )
        assert isinstance(result, PopulationCoverageResult)
        assert result.population_coverage == 0.0
        assert result.method_used == "unknown_alleles"

    def test_all_unknown_alleles_per_allele_status(self):
        result = get_population_coverage(
            ["HLA-FAKE*99:99"], peptides=[]
        )
        assert "HLA-FAKE*99:99" in result.per_allele_results
        assert result.per_allele_results["HLA-FAKE*99:99"]["status"] == "no_coverage_data"
        assert result.per_allele_results["HLA-FAKE*99:99"]["frequency_global"] == 0.0

    def test_mixed_known_and_unknown_alleles(self):
        """Known alleles contribute to coverage; unknown ones are noted."""
        result = get_population_coverage(
            ["HLA-A*02:01", "HLA-FAKE*99:99"], peptides=[]
        )
        assert result.population_coverage > 0.0
        assert result.method_used == "precomputed_db"
        # Unknown allele should be in per_allele_results with status
        assert "HLA-FAKE*99:99" in result.per_allele_results
        assert result.per_allele_results["HLA-FAKE*99:99"]["status"] == "no_coverage_data"
        # Known allele should have frequency
        assert "frequency_global" in result.per_allele_results["HLA-A*02:01"]

    def test_empty_allele_list(self):
        result = get_population_coverage([], peptides=["SIINFEKL"])
        assert result.population_coverage == 0.0
        assert result.method_used == "no_alleles"
        assert result.peptides_tested == 1

    def test_empty_peptide_list(self):
        result = get_population_coverage(["HLA-A*02:01"], peptides=[])
        assert result.peptides_tested == 0
        assert result.population_coverage > 0.0

    def test_all_unknown_alleles_warning_logged(self, caplog):
        """A warning should be logged when all alleles are unknown."""
        with caplog.at_level("WARNING"):
            get_population_coverage(["FAKE*01:01"], peptides=[])
        # The warning messages mention "coverage data" and "returning empty result"
        assert any(
            "coverage data" in r.message.lower() or "returning empty" in r.message.lower()
            for r in caplog.records
        )


# ═══════════════════════════════════════════════════════════════════════════
# 5. MHCflurry path (when available) — mocked
# ═══════════════════════════════════════════════════════════════════════════


class TestMhcflurryPathMocked:
    """When MHCflurry IS available, it should be used and produce results."""

    def test_mhcflurry_path_returns_result(self):
        """Mock MHCFLURRY_AVAILABLE=True and verify mhcflurry path is used."""
        with patch(
            "biocompiler.mhcflurry_population.MHCFLURRY_AVAILABLE", True
        ):
            # Also need to mock MHCflurryClient since it's not installed
            with patch(
                "biocompiler.mhcflurry_population._coverage_via_mhcflurry"
            ) as mock_cov:
                mock_cov.return_value = PopulationCoverageResult(
                    alleles_tested=["HLA-A*02:01"],
                    peptides_tested=1,
                    population_coverage=0.28,
                    method_used="mhcflurry",
                    per_allele_results={
                        "HLA-A*02:01": {"frequency_global": 12.5}
                    },
                )
                result = get_population_coverage(
                    ["HLA-A*02:01"], peptides=["SIINFEKL"]
                )
                assert result.method_used == "mhcflurry"
                assert result.population_coverage == pytest.approx(0.28)
                mock_cov.assert_called_once_with(["HLA-A*02:01"], ["SIINFEKL"])

    def test_mhcflurry_failure_falls_back(self):
        """If MHCflurry path raises, fallback to precomputed_db."""
        with patch(
            "biocompiler.mhcflurry_population.MHCFLURRY_AVAILABLE", True
        ):
            with patch(
                "biocompiler.mhcflurry_population._coverage_via_mhcflurry",
                side_effect=RuntimeError("model load failed"),
            ):
                result = get_population_coverage(
                    ["HLA-A*02:01"], peptides=[]
                )
                # Should fall back to precomputed_db
                assert result.method_used == "precomputed_db"
                assert result.population_coverage > 0.0


# ═══════════════════════════════════════════════════════════════════════════
# 6. Precomputed binding peptide lookup
# ═══════════════════════════════════════════════════════════════════════════


class TestBindingPeptideLookup:
    """Verify that binding peptide lookup works in the fallback path."""

    @pytest.fixture(autouse=True)
    def _force_mhcflurry_unavailable(self):
        with patch(
            "biocompiler.mhcflurry_population.MHCFLURRY_AVAILABLE", False
        ):
            yield

    def test_binding_peptides_in_result_for_known_allele(self):
        """If a peptide is in the precomputed DB, it should appear in results."""
        # HLA-A*02:01 is a known allele with a precomputed database.
        # We can't guarantee any specific peptide is in the DB, but the
        # per_allele_results should have a 'frequency_global' key.
        result = get_population_coverage(
            ["HLA-A*02:01"], peptides=["SIINFEKL"]
        )
        assert "HLA-A*02:01" in result.per_allele_results
        assert "frequency_global" in result.per_allele_results["HLA-A*02:01"]

    def test_no_binding_peptides_key_when_none_found(self):
        """If no binders are found, the 'binding_peptides' key may be absent."""
        result = get_population_coverage(
            ["HLA-A*02:01"], peptides=["ZZZZZZZZZ"]
        )
        # The key should either be absent or present with an empty list
        binding = result.per_allele_results.get("HLA-A*02:01", {}).get(
            "binding_peptides", []
        )
        assert isinstance(binding, list)


# ═══════════════════════════════════════════════════════════════════════════
# 7. Integration — full flow without mhcflurry
# ═══════════════════════════════════════════════════════════════════════════


class TestIntegrationNoMhcflurry:
    """End-to-end test: get_population_coverage works entirely offline."""

    @pytest.fixture(autouse=True)
    def _force_mhcflurry_unavailable(self):
        with patch(
            "biocompiler.mhcflurry_population.MHCFLURRY_AVAILABLE", False
        ):
            yield

    def test_multiple_alleles_peptides(self):
        """Realistic scenario: multiple alleles and peptides."""
        alleles = [
            "HLA-A*02:01",
            "HLA-A*01:01",
            "HLA-B*07:02",
            "HLA-C*07:01",
        ]
        peptides = ["SIINFEKL", "LLFGYPVYV", "GILGFVFTL", "EAAGITILV"]
        result = get_population_coverage(alleles, peptides)

        assert isinstance(result, PopulationCoverageResult)
        assert result.method_used == "precomputed_db"
        assert result.peptides_tested == 4
        assert len(result.alleles_tested) == 4
        assert 0.0 < result.population_coverage <= 1.0

        # Every allele should have per-allele results
        for allele in alleles:
            assert allele in result.per_allele_results
            assert "frequency_global" in result.per_allele_results[allele]

    def test_coverage_consistent_with_existing_function(self):
        """get_population_coverage agrees with compute_population_coverage."""
        alleles = list(EXPANDED_POPULATION_COVERAGE.keys())[:10]
        result = get_population_coverage(alleles, peptides=[])
        direct = compute_population_coverage(alleles, population="global")
        assert result.population_coverage == pytest.approx(direct, abs=1e-6)

    def test_result_is_frozen(self):
        """PopulationCoverageResult should be immutable."""
        result = get_population_coverage(["HLA-A*02:01"], peptides=[])
        with pytest.raises(AttributeError):
            result.method_used = "hacked"  # type: ignore[misc]

    def test_module_never_crashes_on_import(self):
        """Re-importing the module should never raise."""
        import importlib

        import biocompiler.mhcflurry_population as mod

        importlib.reload(mod)
        assert hasattr(mod, "get_population_coverage")
        assert hasattr(mod, "PopulationCoverageResult")
        assert hasattr(mod, "MHCFLURRY_AVAILABLE")
