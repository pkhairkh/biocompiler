"""
Tests for the Multi-Constraint Stress Test
===========================================

Validates that:
  1. Stress tests run without errors
  2. cai_constrained <= cai_unconstrained (constraints always reduce CAI)
  3. Provenance records capture constraint-codon mappings
  4. The report printing function works
"""

from __future__ import annotations

import pytest

from biocompiler.benchmarking.multi_constraint_stress import (
    StressTestResult,
    StressScenario,
    STRESS_SCENARIOS,
    run_stress_test,
    run_all_stress_tests,
    print_stress_test_report,
    _HBB_PROTEIN,
    _GFP_PROTEIN,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def hbb_restriction_result() -> StressTestResult:
    """Run the restriction-site stress test once for HBB (module-scoped)."""
    return run_stress_test(
        protein=_HBB_PROTEIN,
        organism="Homo_sapiens",
        scenario="a_restriction_sites",
    )


@pytest.fixture(scope="module")
def hbb_no_cpg_result() -> StressTestResult:
    """Run the no-CpG stress test once for HBB (module-scoped)."""
    return run_stress_test(
        protein=_HBB_PROTEIN,
        organism="Homo_sapiens",
        scenario="b_no_cpg",
    )


@pytest.fixture(scope="module")
def hbb_everything_result() -> StressTestResult:
    """Run the everything-conflicts stress test once for HBB."""
    return run_stress_test(
        protein=_HBB_PROTEIN,
        organism="Homo_sapiens",
        scenario="c_everything",
    )


@pytest.fixture(scope="module")
def hbb_splice_result() -> StressTestResult:
    """Run the no-cryptic-splice stress test once for HBB."""
    return run_stress_test(
        protein=_HBB_PROTEIN,
        organism="Homo_sapiens",
        scenario="d_no_cryptic_splice",
    )


# ---------------------------------------------------------------------------
# Test 1: Stress tests run without errors
# ---------------------------------------------------------------------------

class TestStressTestRuns:
    """Verify that all stress-test scenarios execute without raising."""

    @pytest.mark.parametrize("scenario_key", list(STRESS_SCENARIOS.keys()))
    def test_hbb_scenario_runs(self, scenario_key: str) -> None:
        """Each HBB scenario should complete without error."""
        result = run_stress_test(
            protein=_HBB_PROTEIN,
            organism="Homo_sapiens",
            scenario=scenario_key,
        )
        assert isinstance(result, StressTestResult)
        assert result.cai_constrained > 0.0, (
            f"Scenario {scenario_key}: constrained CAI should be > 0"
        )

    def test_gfp_restriction_runs(self) -> None:
        """GFP restriction-site scenario should complete without error."""
        result = run_stress_test(
            protein=_GFP_PROTEIN,
            organism="Homo_sapiens",
            scenario="a_restriction_sites",
        )
        assert isinstance(result, StressTestResult)

    def test_run_all_stress_tests(self) -> None:
        """run_all_stress_tests should return 8 results (4 scenarios × 2 proteins)."""
        results = run_all_stress_tests()
        assert len(results) == 8, f"Expected 8 results, got {len(results)}"
        for r in results:
            assert isinstance(r, StressTestResult)


# ---------------------------------------------------------------------------
# Test 2: cai_constrained <= cai_unconstrained
# ---------------------------------------------------------------------------

class TestCAIReduction:
    """Verify that adding constraints never increases CAI."""

    def test_restriction_sites_reduce_cai(self, hbb_restriction_result: StressTestResult) -> None:
        """Restriction site avoidance should reduce (or maintain) CAI."""
        assert hbb_restriction_result.cai_constrained <= hbb_restriction_result.cai_unconstrained + 0.005, (
            f"Constrained CAI ({hbb_restriction_result.cai_constrained:.4f}) should be "
            f"≤ unconstrained CAI ({hbb_restriction_result.cai_unconstrained:.4f}) + tolerance"
        )

    def test_no_cpg_reduces_cai(self, hbb_no_cpg_result: StressTestResult) -> None:
        """CpG avoidance should reduce (or maintain) CAI."""
        assert hbb_no_cpg_result.cai_constrained <= hbb_no_cpg_result.cai_unconstrained + 0.005, (
            f"Constrained CAI ({hbb_no_cpg_result.cai_constrained:.4f}) should be "
            f"≤ unconstrained CAI ({hbb_no_cpg_result.cai_unconstrained:.4f}) + tolerance"
        )

    def test_everything_reduces_cai(self, hbb_everything_result: StressTestResult) -> None:
        """The 'everything conflicts' scenario should reduce CAI."""
        assert hbb_everything_result.cai_constrained <= hbb_everything_result.cai_unconstrained + 0.005, (
            f"Constrained CAI ({hbb_everything_result.cai_constrained:.4f}) should be "
            f"≤ unconstrained CAI ({hbb_everything_result.cai_unconstrained:.4f}) + tolerance"
        )

    def test_splice_reduces_cai(self, hbb_splice_result: StressTestResult) -> None:
        """Cryptic splice avoidance should reduce (or maintain) CAI."""
        assert hbb_splice_result.cai_constrained <= hbb_splice_result.cai_unconstrained + 0.005, (
            f"Constrained CAI ({hbb_splice_result.cai_constrained:.4f}) should be "
            f"≤ unconstrained CAI ({hbb_splice_result.cai_unconstrained:.4f}) + tolerance"
        )


# ---------------------------------------------------------------------------
# Test 3: Provenance records capture constraint-codon mappings
# ---------------------------------------------------------------------------

class TestProvenanceCapture:
    """Verify that provenance records contain constraint-codon attribution."""

    def test_restriction_result_has_provenance(self, hbb_restriction_result: StressTestResult) -> None:
        """Restriction-site result should have non-empty provenance records."""
        assert len(hbb_restriction_result.provenance_records) > 0, (
            "Provenance records should not be empty"
        )

    def test_cai_loss_per_constraint_populated(self, hbb_everything_result: StressTestResult) -> None:
        """The 'everything' scenario should have per-constraint CAI loss data."""
        # Even if we can't attribute to individual constraints, we should
        # at least have aggregate data
        assert len(hbb_everything_result.cai_loss_per_constraint) > 0 or \
               hbb_everything_result.total_cai_loss != 0.0, (
            "Should have some CAI loss data for the 'everything' scenario"
        )

    def test_codons_changed_per_constraint_populated(self, hbb_restriction_result: StressTestResult) -> None:
        """If CAI changed, codon changes should be tracked."""
        if hbb_restriction_result.total_cai_loss < -0.001:
            # There was a meaningful CAI loss → we should track which codons changed
            total_changed = sum(hbb_restriction_result.codons_changed_per_constraint.values())
            assert total_changed > 0, (
                "If CAI decreased, some codons must have changed"
            )

    def test_most_expensive_constraint_set(self, hbb_everything_result: StressTestResult) -> None:
        """The most_expensive_constraint field should be populated."""
        assert hbb_everything_result.most_expensive_constraint, (
            "most_expensive_constraint should be a non-empty string"
        )

    def test_provenance_contains_decision_trail(self, hbb_restriction_result: StressTestResult) -> None:
        """Provenance records should contain the decision trail object."""
        from biocompiler.decision_provenance import OptimizationDecisionTrail
        has_trail = any(
            isinstance(r, OptimizationDecisionTrail)
            for r in hbb_restriction_result.provenance_records
        )
        assert has_trail, "Provenance should contain an OptimizationDecisionTrail"


# ---------------------------------------------------------------------------
# Test 4: Report printing function
# ---------------------------------------------------------------------------

class TestReportPrinting:
    """Verify that print_stress_test_report produces formatted output."""

    def test_report_returns_string(self) -> None:
        """Report function should return a non-empty string."""
        # Run a minimal test with just the restriction-site scenario on HBB
        results = [
            run_stress_test(
                protein=_HBB_PROTEIN,
                organism="Homo_sapiens",
                scenario="a_restriction_sites",
            ),
        ]
        report = print_stress_test_report(results)
        assert isinstance(report, str)
        assert len(report) > 100, "Report should be substantial"
        assert "MULTI-CONSTRAINT STRESS TEST" in report
        assert "CAI Loss Breakdown" in report or "No per-constraint" in report

    def test_report_with_all_results(self) -> None:
        """Report with all results should contain scenario summaries."""
        results = run_all_stress_tests()
        report = print_stress_test_report(results)
        assert "SUMMARY" in report
        assert "Most Expensive" in report
        # Should mention at least one constraint type
        assert any(
            name in report
            for name in ["NoRestrictionSite", "NoCpG", "GCInRange", "NoCrypticSplice", "all_constraints"]
        ), f"Report should mention constraint types; got:\n{report[:500]}"


# ---------------------------------------------------------------------------
# Test 5: StressScenario and StressTestResult data classes
# ---------------------------------------------------------------------------

class TestDataClasses:
    """Test the data class constructors and properties."""

    def test_stress_scenario_creation(self) -> None:
        """StressScenario should be constructable and frozen."""
        s = StressScenario(
            name="Test Scenario",
            description="Test description",
            enzymes=["EcoRI"],
            avoid_cpg=True,
            gc_lo=0.40,
            gc_hi=0.60,
        )
        assert s.name == "Test Scenario"
        assert s.avoid_cpg is True
        with pytest.raises(AttributeError):
            s.name = "Changed"  # type: ignore[misc]

    def test_stress_test_result_total_loss(self) -> None:
        """total_cai_loss property should compute correctly."""
        r = StressTestResult(
            scenario_name="Test",
            protein="MV",
            organism="Homo_sapiens",
            cai_unconstrained=0.95,
            cai_constrained=0.88,
            cai_loss_per_constraint={"NoCpG": -0.07},
            codons_changed_per_constraint={"NoCpG": 5},
            provenance_records=[],
            most_expensive_constraint="NoCpG",
        )
        assert abs(r.total_cai_loss - (-0.07)) < 0.001

    def test_predefined_scenarios(self) -> None:
        """All predefined scenarios should have required attributes."""
        for key, scenario in STRESS_SCENARIOS.items():
            assert scenario.name, f"Scenario {key} missing name"
            assert scenario.description, f"Scenario {key} missing description"
            assert 0.0 < scenario.gc_lo < scenario.gc_hi <= 1.0, (
                f"Scenario {key} has invalid GC bounds: [{scenario.gc_lo}, {scenario.gc_hi}]"
            )


# ---------------------------------------------------------------------------
# Test 6: Edge case — custom scenario
# ---------------------------------------------------------------------------

class TestCustomScenario:
    """Test with a custom StressScenario passed directly."""

    def test_custom_scenario(self) -> None:
        """A custom StressScenario should work with run_stress_test."""
        custom = StressScenario(
            name="Custom: Tight GC",
            description="GC [0.48, 0.52] only",
            gc_lo=0.48,
            gc_hi=0.52,
        )
        result = run_stress_test(
            protein=_HBB_PROTEIN,
            organism="Homo_sapiens",
            scenario=custom,
        )
        assert isinstance(result, StressTestResult)
        assert result.cai_constrained > 0.0

    def test_invalid_scenario_key(self) -> None:
        """An unknown scenario key should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown scenario"):
            run_stress_test(
                protein=_HBB_PROTEIN,
                organism="Homo_sapiens",
                scenario="nonexistent_scenario",
            )
