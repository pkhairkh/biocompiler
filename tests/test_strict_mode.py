"""
Tests for the strict_mode hard-stop feature.

When strict_mode=True (the default), optimize_sequence() raises
OptimizationConstraintError instead of returning a result with
failed predicates.  This guarantees that callers never silently
receive a sequence that violates constraints.

Tests:
1. strict_mode=True raises OptimizationConstraintError on failures
2. strict_mode=False returns result with failed_predicates populated
3. API returns HTTP 422 with failure details when strict_mode=True
4. Exception contains useful diagnostic information
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from biocompiler.optimizer import optimize_sequence, OptimizationResult
from biocompiler.shared.exceptions import (
    OptimizationConstraintError,
    BioCompilerError,
)


# ─── Test Data ───────────────────────────────────────────────────────

# A short protein that is easy to optimize
_SIMPLE_PROTEIN = "MVSKGE"

# A protein with many Val (V) codons — Val codons all contain GT,
# which triggers GT-dinucleotide / cryptic splice site predicates for eukaryotes.
_VAL_HEAVY_PROTEIN = "MVVVLVVFVV"


# ═══════════════════════════════════════════════════════════════════════
# 1. optimize_sequence — strict_mode=True raises on failures
# ═══════════════════════════════════════════════════════════════════════

class TestStrictModeRaises:
    """Verify that strict_mode=True refuses to return failed results."""

    def test_strict_mode_raises_on_impossible_gc(self):
        """An impossible GC range (0.99–1.0) must raise OptimizationConstraintError."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        # The exception must list at least one failed predicate
        assert len(exc_info.value.failed_predicates) >= 1

    def test_strict_mode_default_is_true(self):
        """When strict_mode is not specified, it defaults to True, so
        impossible constraints should raise."""
        with pytest.raises(OptimizationConstraintError):
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                # strict_mode NOT specified — should default to True
            )

    def test_strict_mode_raises_is_biocompiler_error(self):
        """OptimizationConstraintError must be a subclass of BioCompilerError."""
        with pytest.raises(BioCompilerError):
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )

    def test_strict_mode_passes_when_all_predicates_satisfied(self):
        """With reasonable constraints and strict_mode=True, a valid
        optimization should succeed without raising."""
        result = optimize_sequence(
            target_protein=_SIMPLE_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.20,
            gc_hi=0.80,
            strict_mode=True,
        )
        assert isinstance(result, OptimizationResult)
        assert result.failed_predicates == []


# ═══════════════════════════════════════════════════════════════════════
# 2. optimize_sequence — strict_mode=False returns partial results
# ═══════════════════════════════════════════════════════════════════════

class TestStrictModeFalseReturnsPartial:
    """Verify that strict_mode=False returns results even with failures."""

    def test_strict_mode_false_returns_result_with_failures(self):
        """With impossible GC and strict_mode=False, we get a result with
        failed_predicates populated."""
        result = optimize_sequence(
            target_protein=_SIMPLE_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.99,
            gc_hi=1.0,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert len(result.failed_predicates) >= 1

    def test_strict_mode_false_result_has_sequence(self):
        """The partial result should still contain a valid sequence."""
        result = optimize_sequence(
            target_protein=_SIMPLE_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.99,
            gc_hi=1.0,
            strict_mode=False,
        )
        assert result.sequence is not None
        assert len(result.sequence) > 0

    def test_strict_mode_false_result_encodes_protein(self):
        """The partial result's sequence should still encode the original protein."""
        from biocompiler.expression.translation import translate
        result = optimize_sequence(
            target_protein=_SIMPLE_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.99,
            gc_hi=1.0,
            strict_mode=False,
        )
        translated = translate(result.sequence)
        assert translated == _SIMPLE_PROTEIN

    def test_strict_mode_false_no_failures_returns_success(self):
        """With reasonable constraints and strict_mode=False, a fully
        satisfied result has empty failed_predicates."""
        result = optimize_sequence(
            target_protein=_SIMPLE_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.20,
            gc_hi=0.80,
            strict_mode=False,
        )
        assert isinstance(result, OptimizationResult)
        assert result.failed_predicates == []


# ═══════════════════════════════════════════════════════════════════════
# 3. API endpoint — HTTP 422 on strict mode failure
# ═══════════════════════════════════════════════════════════════════════

class TestAPIStrictMode:
    """Test the /optimize API endpoint strict_mode behavior."""

    @pytest.fixture()
    def client(self):
        """Create a fresh TestClient for each test with auth disabled."""
        import os
        # Must set env var before importing api module
        old_auth = os.environ.get("BIOCOMPILER_AUTH_MODE")
        os.environ["BIOCOMPILER_AUTH_MODE"] = "disabled"
        # Force re-import to pick up env var
        import importlib
        import biocompiler.api as api_mod
        importlib.reload(api_mod)
        from biocompiler.api import create_app
        api_mod._rate_limiter.clear()
        client = TestClient(create_app())
        yield client
        if old_auth is not None:
            os.environ["BIOCOMPILER_AUTH_MODE"] = old_auth
        else:
            os.environ.pop("BIOCOMPILER_AUTH_MODE", None)

    def test_api_strict_mode_true_returns_422_on_failure(self, client):
        """POST /optimize with strict_mode=True and impossible constraints
        returns HTTP 422."""
        resp = client.post("/optimize", json={
            "protein": _SIMPLE_PROTEIN,
            "gc_lo": 0.99,
            "gc_hi": 1.0,
            "strict_mode": True,
        })
        assert resp.status_code == 422

    def test_api_strict_mode_false_returns_200_on_failure(self, client):
        """POST /optimize with strict_mode=False and impossible constraints
        returns HTTP 200 with failed_predicates populated."""
        resp = client.post("/optimize", json={
            "protein": _SIMPLE_PROTEIN,
            "gc_lo": 0.99,
            "gc_hi": 1.0,
            "strict_mode": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["failed_predicates"]) >= 1

    def test_api_422_response_contains_failure_details(self, client):
        """The 422 response body should contain failed_predicates list
        and a helpful suggestion."""
        resp = client.post("/optimize", json={
            "protein": _SIMPLE_PROTEIN,
            "gc_lo": 0.99,
            "gc_hi": 1.0,
            "strict_mode": True,
        })
        assert resp.status_code == 422
        detail = resp.json()["detail"]
        # FastAPI wraps the detail; it may be a dict or string depending
        # on how the HTTPException is serialized
        if isinstance(detail, dict):
            assert "failed_predicates" in detail
            assert len(detail["failed_predicates"]) >= 1
            assert "suggestion" in detail
            assert "strict_mode=False" in detail["suggestion"]
        elif isinstance(detail, str):
            # If detail is a string, it should at least mention the failure
            assert "strict mode" in detail.lower() or "predicate" in detail.lower()

    def test_api_422_response_has_partial_result(self, client):
        """The 422 response should include a partial_result for inspection."""
        resp = client.post("/optimize", json={
            "protein": _SIMPLE_PROTEIN,
            "gc_lo": 0.99,
            "gc_hi": 1.0,
            "strict_mode": True,
        })
        assert resp.status_code == 422
        detail = resp.json()["detail"]
        if isinstance(detail, dict):
            assert "partial_result" in detail
            partial = detail["partial_result"]
            assert partial is not None
            assert partial["sequence"] is not None
            assert partial["protein"] == _SIMPLE_PROTEIN

    def test_api_strict_mode_default_true(self, client):
        """When strict_mode is not specified in the API, it defaults to True,
        so impossible constraints should trigger 422."""
        resp = client.post("/optimize", json={
            "protein": _SIMPLE_PROTEIN,
            "gc_lo": 0.99,
            "gc_hi": 1.0,
            # strict_mode NOT specified — default is True
        })
        assert resp.status_code == 422

    def test_api_succeeds_with_reasonable_constraints(self, client):
        """With reasonable constraints and strict_mode=True, optimization
        succeeds and returns HTTP 200."""
        resp = client.post("/optimize", json={
            "protein": _SIMPLE_PROTEIN,
            "organism": "ecoli",
            "gc_lo": 0.20,
            "gc_hi": 0.80,
            "strict_mode": True,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["sequence"] is not None


# ═══════════════════════════════════════════════════════════════════════
# 4. Exception diagnostics
# ═══════════════════════════════════════════════════════════════════════

class TestExceptionDiagnostics:
    """Verify OptimizationConstraintError contains useful information."""

    def test_exception_has_failed_predicates(self):
        """The exception must store the list of failed predicates."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        err = exc_info.value
        assert isinstance(err.failed_predicates, list)
        assert len(err.failed_predicates) >= 1
        # Each predicate name should be a non-empty string
        for pred in err.failed_predicates:
            assert isinstance(pred, str)
            assert len(pred) > 0

    def test_exception_has_partial_result(self):
        """The exception must provide the partial result for inspection."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        err = exc_info.value
        assert err.partial_result is not None
        assert isinstance(err.partial_result, OptimizationResult)
        assert err.partial_result.sequence is not None
        assert len(err.partial_result.sequence) > 0

    def test_exception_message_contains_count(self):
        """The exception message should mention the number of failed predicates."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        msg = str(exc_info.value)
        assert "predicate(s)" in msg
        # Should mention the count
        n = len(exc_info.value.failed_predicates)
        assert str(n) in msg

    def test_exception_message_contains_predicate_names(self):
        """The exception message should list the failed predicate names."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        msg = str(exc_info.value)
        for pred in exc_info.value.failed_predicates:
            assert pred in msg

    def test_exception_message_suggests_relaxing(self):
        """The exception message should suggest setting strict_mode=False."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        msg = str(exc_info.value).lower()
        assert "strict_mode=false" in msg or "strict_mode = false" in msg

    def test_partial_result_has_cai_and_gc(self):
        """The partial result should have CAI and GC metrics for inspection."""
        with pytest.raises(OptimizationConstraintError) as exc_info:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        partial = exc_info.value.partial_result
        assert partial.cai is not None
        assert 0.0 <= partial.cai <= 1.0
        assert partial.gc_content is not None
        assert 0.0 <= partial.gc_content <= 1.0

    def test_exception_is_catchable_as_biocompiler_error(self):
        """OptimizationConstraintError should be catchable as BioCompilerError
        for generic error handlers."""
        caught = False
        try:
            optimize_sequence(
                target_protein=_SIMPLE_PROTEIN,
                organism="Homo_sapiens",
                gc_lo=0.99,
                gc_hi=1.0,
                strict_mode=True,
            )
        except BioCompilerError:
            caught = True
        assert caught
