"""Comprehensive tests for uncertainty enforcement in BioCompiler.

Validates all uncertainty enforcement features including:
- Certificate uncertainty capping (GOLD/SILVER/BRONZE based on UNCERTAIN verdicts)
- Uncertainty warnings emitted when predicates return UNCERTAIN
- Uncertainty summary statistics
- Solver verdict mapping (hard/soft constraints → Verdict)
- Assessment verdict mapping (STABLE/MARGINAL/UNSTABLE → Verdict)
"""

from __future__ import annotations

import logging

import pytest

from biocompiler.provenance.certificate import (
    compute_certificate,
    compute_uncertainty_summary,
)
from biocompiler.shared.types import Verdict
from biocompiler.type_system.codon_tables import CertLevel, PredicateResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_result(
    predicate: str, verdict: Verdict, passed: bool = True
) -> PredicateResult:
    """Create a PredicateResult with a specific verdict.

    Args:
        predicate: Name of the predicate.
        verdict: The Verdict to assign.
        passed: Whether the predicate is considered passed (default True).

    Returns:
        A PredicateResult with the specified verdict.
    """
    return PredicateResult(
        predicate=predicate,
        passed=passed,
        verdict=verdict,
        details="test",
    )


# ===========================================================================
# TestCertificateUncertaintyCapping
# ===========================================================================

class TestCertificateUncertaintyCapping:
    """Tests for how UNCERTAIN verdicts cap certificate levels.

    Uncertainty capping rules (from compute_certificate docstring):
      1. Any Verdict.UNCERTAIN → certificate cannot exceed SILVER
      2. Multiple UNCERTAIN verdicts → certificate cannot exceed BRONZE
      3. LIKELY_PASS/LIKELY_FAIL do NOT cap the certificate
    """

    def test_no_uncertain_gold(self) -> None:
        """All PASS → GOLD."""
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("GCInRange", Verdict.PASS),
            _make_result("NoCpGIsland", Verdict.PASS),
        ]
        assert compute_certificate(results) == CertLevel.GOLD

    def test_one_uncertain_silver(self) -> None:
        """One UNCERTAIN → SILVER (even if all pass)."""
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result("GCInRange", Verdict.PASS),
        ]
        assert compute_certificate(results) == CertLevel.SILVER

    def test_two_uncertain_bronze(self) -> None:
        """Two UNCERTAIN → BRONZE."""
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result("SolubleExpression", Verdict.UNCERTAIN),
            _make_result("GCInRange", Verdict.PASS),
        ]
        assert compute_certificate(results) == CertLevel.BRONZE

    def test_uncertain_overrides_mutagenesis(self) -> None:
        """One UNCERTAIN + mutagenesis → SILVER (not BRONZE).

        When there is exactly one UNCERTAIN and mutagenesis applied,
        the certificate should still be SILVER because:
        - uncertain_count == 1 → SILVER cap
        - mutagenesis alone also → SILVER
        The combined result is SILVER, not BRONZE.
        """
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result(
                "NoGTDinucleotide",
                Verdict.PASS,
                passed=True,
            ),
        ]
        # Add mutagenesis via the details string
        results[-1] = PredicateResult(
            predicate="NoGTDinucleotide",
            passed=True,
            verdict=Verdict.PASS,
            details="mutagenesis applied: pos 3:V→I",
            mutagenesis_applied=True,
        )
        # 1 UNCERTAIN → SILVER; mutagenesis also → SILVER; combined → SILVER
        assert compute_certificate(results) == CertLevel.SILVER

    def test_likely_pass_does_not_cap(self) -> None:
        """LIKELY_PASS should NOT cap the certificate."""
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("StableFolding", Verdict.LIKELY_PASS),
            _make_result("GCInRange", Verdict.PASS),
        ]
        assert compute_certificate(results) == CertLevel.GOLD

    def test_likely_fail_does_not_cap(self) -> None:
        """LIKELY_FAIL should NOT cap the certificate (but it makes the
        predicate fail, which results in BRONZE for unsatisfied).

        Note: LIKELY_FAIL does not cap via the uncertainty rule, but
        if passed=False (unsatisfied), it triggers BRONZE via the
        has_unsatisfied path. To test that LIKELY_FAIL alone does not
        cap via uncertainty, we need passed=True.
        """
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("CodonOptimality", Verdict.LIKELY_FAIL, passed=True),
            _make_result("GCInRange", Verdict.PASS),
        ]
        # LIKELY_FAIL does not trigger the uncertainty capping rule,
        # and since passed=True, there is no unsatisfied predicate.
        # So the certificate should be GOLD (no mutagenesis, no unavoidable).
        assert compute_certificate(results) == CertLevel.GOLD

    def test_uncertain_and_fail_bronze(self) -> None:
        """FAIL + UNCERTAIN → BRONZE.

        A failed predicate (has_unsatisfied=True) already caps at BRONZE.
        Adding an UNCERTAIN does not change that.
        """
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("NoRestrictionSite", Verdict.FAIL, passed=False),
            _make_result("StableFolding", Verdict.UNCERTAIN),
        ]
        assert compute_certificate(results) == CertLevel.BRONZE

    def test_all_uncertain_bronze(self) -> None:
        """All UNCERTAIN → BRONZE (2+ UNCERTAIN)."""
        results = [
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result("SolubleExpression", Verdict.UNCERTAIN),
            _make_result("StructureConfidence", Verdict.UNCERTAIN),
        ]
        assert compute_certificate(results) == CertLevel.BRONZE


# ===========================================================================
# TestUncertaintyWarnings
# ===========================================================================

class TestUncertaintyWarnings:
    """Tests for uncertainty warning emission.

    The _emit_uncertainty_warnings function in pipeline_core.py logs
    warnings when predicates have UNCERTAIN verdicts. These tests
    verify the warning behavior.
    """

    def test_no_uncertain_no_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        """No UNCERTAIN → no warning logged."""
        from biocompiler.optimizer.pipeline_core import _emit_uncertainty_warnings

        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("GCInRange", Verdict.PASS),
        ]
        log = logging.getLogger("test_uncertainty_warnings")
        with caplog.at_level(logging.WARNING):
            count = _emit_uncertainty_warnings(results, log)
        assert count == 0
        # No UNCERTAIN-related warnings should be emitted
        assert not any(
            "UNCERTAIN" in record.message for record in caplog.records
        )

    def test_uncertain_emits_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        """UNCERTAIN → warning logged with predicate names."""
        from biocompiler.optimizer.pipeline_core import _emit_uncertainty_warnings

        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result("GCInRange", Verdict.PASS),
        ]
        log = logging.getLogger("test_uncertainty_warnings")
        with caplog.at_level(logging.WARNING):
            count = _emit_uncertainty_warnings(results, log)
        assert count == 1
        # Warning should mention the predicate name
        assert any(
            "StableFolding" in record.message for record in caplog.records
        )

    def test_multiple_uncertain_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        """6+ UNCERTAIN → truncated warning with '...'."""
        from biocompiler.optimizer.pipeline_core import _emit_uncertainty_warnings

        results = [
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result("SolubleExpression", Verdict.UNCERTAIN),
            _make_result("StructureConfidence", Verdict.UNCERTAIN),
            _make_result("NoMisfoldingRisk", Verdict.UNCERTAIN),
            _make_result("CorrectFoldTopology", Verdict.UNCERTAIN),
            _make_result("HydrophobicCoreQuality", Verdict.UNCERTAIN),
        ]
        log = logging.getLogger("test_uncertainty_warnings")
        with caplog.at_level(logging.WARNING):
            count = _emit_uncertainty_warnings(results, log)
        assert count == 6
        # With 6 predicates (more than 5), the warning should contain "..."
        assert any(
            "..." in record.message for record in caplog.records
        )


# ===========================================================================
# TestUncertaintySummary
# ===========================================================================

class TestUncertaintySummary:
    """Tests for compute_uncertainty_summary statistics."""

    def test_summary_all_pass(self) -> None:
        """All PASS → confidence=1.0, uncertain_count=0."""
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("GCInRange", Verdict.PASS),
            _make_result("NoCpGIsland", Verdict.PASS),
        ]
        summary = compute_uncertainty_summary(results)
        assert summary["uncertain_count"] == 0
        assert summary["confidence_score"] == 1.0
        assert summary["definite_count"] == 3
        assert summary["likely_pass_count"] == 0
        assert summary["likely_fail_count"] == 0

    def test_summary_mixed(self) -> None:
        """Mixed verdicts → verify all counts."""
        results = [
            _make_result("NoStopCodons", Verdict.PASS),         # definite +1, confidence +1.0
            _make_result("StableFolding", Verdict.UNCERTAIN),   # uncertain +1, confidence +0.5
            _make_result("CodonOptimality", Verdict.LIKELY_PASS),  # likely_pass +1, confidence +0.75
            _make_result("SolubleExpression", Verdict.FAIL, passed=False),  # definite +1, confidence +0.0
            _make_result("StructureConfidence", Verdict.LIKELY_FAIL),  # likely_fail +1, confidence +0.25
        ]
        summary = compute_uncertainty_summary(results)
        assert summary["total_predicates"] == 5
        assert summary["uncertain_count"] == 1
        assert summary["likely_pass_count"] == 1
        assert summary["likely_fail_count"] == 1
        assert summary["definite_count"] == 2  # PASS + FAIL
        assert summary["uncertain_predicates"] == ["StableFolding"]

    def test_summary_confidence_score(self) -> None:
        """Verify weighted average confidence score.

        2 PASS (1.0 each) + 1 UNCERTAIN (0.5) + 1 LIKELY_PASS (0.75) = 3.25 / 4
        """
        results = [
            _make_result("NoStopCodons", Verdict.PASS),
            _make_result("GCInRange", Verdict.PASS),
            _make_result("StableFolding", Verdict.UNCERTAIN),
            _make_result("CodonOptimality", Verdict.LIKELY_PASS),
        ]
        summary = compute_uncertainty_summary(results)
        expected_confidence = (1.0 + 1.0 + 0.5 + 0.75) / 4.0
        assert abs(summary["confidence_score"] - round(expected_confidence, 3)) < 0.001


# ===========================================================================
# TestSolverVerdictMapping
# ===========================================================================

class TestSolverVerdictMapping:
    """Tests for solver constraint enforcement → type-system verdict mapping.

    Verifies the SOLVER_VERDICT_MAP:
    - (HARD, satisfied)  → PASS
    - (HARD, violated)   → FAIL
    - (SOFT, satisfied)  → PASS
    - (SOFT, violated)   → UNCERTAIN
    """

    def test_hard_satisfied_pass(self) -> None:
        """Hard constraint satisfied → PASS."""
        from biocompiler.solver.types import (
            ConstraintStrictness,
            SOLVER_VERDICT_MAP,
        )

        verdict_str = SOLVER_VERDICT_MAP[(ConstraintStrictness.HARD, True)]
        assert verdict_str == "PASS"

    def test_hard_violated_fail(self) -> None:
        """Hard constraint violated → FAIL."""
        from biocompiler.solver.types import (
            ConstraintStrictness,
            SOLVER_VERDICT_MAP,
        )

        verdict_str = SOLVER_VERDICT_MAP[(ConstraintStrictness.HARD, False)]
        assert verdict_str == "FAIL"

    def test_soft_satisfied_pass(self) -> None:
        """Soft constraint satisfied → PASS."""
        from biocompiler.solver.types import (
            ConstraintStrictness,
            SOLVER_VERDICT_MAP,
        )

        verdict_str = SOLVER_VERDICT_MAP[(ConstraintStrictness.SOFT, True)]
        assert verdict_str == "PASS"

    def test_soft_violated_uncertain(self) -> None:
        """Soft constraint violated → UNCERTAIN."""
        from biocompiler.solver.types import (
            ConstraintStrictness,
            SOLVER_VERDICT_MAP,
        )

        verdict_str = SOLVER_VERDICT_MAP[(ConstraintStrictness.SOFT, False)]
        assert verdict_str == "UNCERTAIN"


# ===========================================================================
# TestAssessmentVerdictMapping
# ===========================================================================

class TestAssessmentVerdictMapping:
    """Tests for assessment service verdict → type-system Verdict mapping.

    Verifies assessment_verdict_to_verdict:
    - "STABLE"   → Verdict.PASS
    - "MARGINAL" → Verdict.UNCERTAIN
    - "UNSTABLE" → Verdict.FAIL
    - unknown    → Verdict.UNCERTAIN
    """

    def test_stable_to_pass(self) -> None:
        """'STABLE' → PASS."""
        from biocompiler.application.assessment_service import (
            assessment_verdict_to_verdict,
        )

        assert assessment_verdict_to_verdict("STABLE") == Verdict.PASS

    def test_marginal_to_uncertain(self) -> None:
        """'MARGINAL' → UNCERTAIN."""
        from biocompiler.application.assessment_service import (
            assessment_verdict_to_verdict,
        )

        assert assessment_verdict_to_verdict("MARGINAL") == Verdict.UNCERTAIN

    def test_unstable_to_fail(self) -> None:
        """'UNSTABLE' → FAIL."""
        from biocompiler.application.assessment_service import (
            assessment_verdict_to_verdict,
        )

        assert assessment_verdict_to_verdict("UNSTABLE") == Verdict.FAIL

    def test_unknown_to_uncertain(self) -> None:
        """Unknown verdict string → UNCERTAIN (safe default)."""
        from biocompiler.application.assessment_service import (
            assessment_verdict_to_verdict,
        )

        assert assessment_verdict_to_verdict("UNKNOWN") == Verdict.UNCERTAIN
        assert assessment_verdict_to_verdict("") == Verdict.UNCERTAIN
        assert assessment_verdict_to_verdict("something_else") == Verdict.UNCERTAIN
