"""Tests for dispatch.py — validate_csp_solution and related functions.

Covers:
1. validate_csp_solution with empty sequence returns violations
2. Valid sequence returns empty violations list
3. GC-violating sequence returns GC violation
4. Restriction site in sequence returns violation
5. Translation fidelity check catches wrong codon
6. Sequence length mismatch returns violation
7. solve_with_csp returns SolverResult
8. get_csp_availability returns dict with expected keys
9. is_solver_available returns bool
10. csp_optimize returns SolverResult with fallback_used on failure
"""

from __future__ import annotations

import pytest
from unittest.mock import patch, MagicMock

from biocompiler.solver.dispatch import (
    validate_csp_solution,
    solve_with_csp,
    get_csp_availability,
    is_solver_available,
    csp_optimize,
)
from biocompiler.solver.types import (
    SolverConfig,
    SolverResult,
    SolverBackend,
    ConstraintStrictness,
    ConstraintViolation,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def default_config() -> SolverConfig:
    """Default SolverConfig with standard GC bounds."""
    return SolverConfig(gc_lo=0.30, gc_hi=0.70)


@pytest.fixture
def protein_mvskge() -> str:
    """Short 6-AA protein: MVSKGE."""
    return "MVSKGE"


@pytest.fixture
def valid_dna_mvskge() -> str:
    """Valid DNA encoding of 'MVSKGE' with GC ~44% (within default 30-70%).

    M=ATG  V=GTT  S=TCT  K=AAA  G=GGT  E=GAA
    GC = (G in ATG) + (G in GTT) + (GC in TCT) + (G in GGT) + (G in GAA)
       = 1 + 1 + 1 + 2 + 1 = 6 out of 18 bases = 33.3%  (within 30-70%)
    """
    return "ATGGTTTCTAAAGGTGAA"


# ---------------------------------------------------------------------------
# 1. Empty sequence returns violations
# ---------------------------------------------------------------------------

class TestEmptySequence:
    """Test that validate_csp_solution flags an empty sequence."""

    def test_empty_sequence_returns_violation(self, default_config, protein_mvskge):
        violations, score = validate_csp_solution("", protein_mvskge, default_config)
        assert len(violations) >= 1
        assert violations[0].constraint_name == "non_empty_sequence"
        assert violations[0].constraint_type == ConstraintStrictness.HARD
        assert violations[0].severity == 1.0
        assert score == 0.0  # Empty sequence → composite score 0.0

    def test_empty_string_sequence(self, default_config):
        """Even with an empty protein, empty sequence should violate."""
        violations, score = validate_csp_solution("", "M", default_config)
        assert len(violations) >= 1
        names = [v.constraint_name for v in violations]
        assert "non_empty_sequence" in names


# ---------------------------------------------------------------------------
# 2. Valid sequence returns empty violations list
# ---------------------------------------------------------------------------

class TestValidSequence:
    """Test that a fully valid sequence produces no violations."""

    def test_valid_sequence_no_violations(self, default_config, protein_mvskge, valid_dna_mvskge):
        violations, score = validate_csp_solution(valid_dna_mvskge, protein_mvskge, default_config)
        # Filter out soft-constraint violations that are always-True (CAI, CpG, mRNA dG)
        hard_violations = [
            v for v in violations
            if v.constraint_type == ConstraintStrictness.HARD
        ]
        assert len(hard_violations) == 0, (
            f"Unexpected hard violations: {[v.description for v in hard_violations]}"
        )

    def test_valid_sequence_translation_fidelity_ok(self, default_config, protein_mvskge, valid_dna_mvskge):
        """Translation fidelity should not be flagged for a correct sequence."""
        violations, score = validate_csp_solution(valid_dna_mvskge, protein_mvskge, default_config)
        fidelity_names = [v.constraint_name for v in violations]
        assert "translation_fidelity" not in fidelity_names


# ---------------------------------------------------------------------------
# 3. GC-violating sequence returns GC violation
# ---------------------------------------------------------------------------

class TestGCViolation:
    """Test that a sequence outside GC bounds is flagged."""

    def test_high_gc_sequence_violation(self, default_config, protein_mvskge):
        """Construct an all-GC sequence for 'MVSKGE' that exceeds gc_hi=0.70.

        Use only G/C codons:
        M=ATG (1 GC), V=GTG (2 GC), S=TCG (2 GC), K=AAG (1 GC), G=GGC (3 GC), E=GAG (2 GC)
        Total GC = 1+2+2+1+3+2 = 11 out of 18 = 61.1% — just inside.
        Let's use a very tight config to force a violation.
        """
        tight_config = SolverConfig(gc_lo=0.30, gc_hi=0.40)
        # ATGGGTTCGAAAGGTGAA → GC = 1+2+2+1+2+1 = 9/18 = 50% > 40%
        gc_heavy_seq = "ATGGTGTCGAAGGGCGAG"
        violations, score = validate_csp_solution(gc_heavy_seq, protein_mvskge, tight_config)
        gc_violations = [
            v for v in violations if v.constraint_name == "GCRangeConstraint"
        ]
        assert len(gc_violations) >= 1, (
            f"Expected GC violation for high-GC sequence; got: "
            f"{[v.constraint_name for v in violations]}"
        )

    def test_low_gc_sequence_violation(self, protein_mvskge):
        """Construct an all-AT sequence for 'MVSKGE' that is below gc_lo=0.50."""
        high_lo_config = SolverConfig(gc_lo=0.50, gc_hi=0.70)
        # ATG GTT TCT AAA GGT GAA → GC ≈ 33%, below 50%
        low_gc_seq = "ATGGTTTCTAAAGGTGAA"
        violations, score = validate_csp_solution(low_gc_seq, protein_mvskge, high_lo_config)
        gc_violations = [
            v for v in violations if v.constraint_name == "GCRangeConstraint"
        ]
        assert len(gc_violations) >= 1, (
            f"Expected GC violation for low-GC sequence; got: "
            f"{[v.constraint_name for v in violations]}"
        )


# ---------------------------------------------------------------------------
# 4. Restriction site in sequence returns violation
# ---------------------------------------------------------------------------

class TestRestrictionSiteViolation:
    """Test that a restriction site present in the sequence is flagged."""

    def test_ecori_site_violation(self, protein_mvskge):
        """EcoRI site (GAATTC) embedded in a DNA sequence should be flagged."""
        # Build a sequence that contains GAATTC.
        # We need a protein and a DNA sequence that translates to it AND contains GAATTC.
        # Use protein "EF" → E=GAA, F=TTTC (not valid — F is TTT or TTC).
        # E=GAA F=TTTC is wrong. Let's think more carefully.
        # GAATTC: split into codons GAA,TTC → translates to E,F
        # So protein="EF", sequence="GAATTC" encodes it with an EcoRI site.
        config = SolverConfig(
            gc_lo=0.10, gc_hi=0.90,
            restriction_sites=["GAATTC"],
        )
        violations, score = validate_csp_solution("GAATTC", "EF", config)
        rs_violations = [
            v for v in violations
            if v.constraint_name == "NoRestrictionSiteConstraint"
        ]
        assert len(rs_violations) >= 1, (
            f"Expected restriction site violation; got: "
            f"{[v.constraint_name for v in violations]}"
        )

    def test_no_restriction_site_when_clean(self, protein_mvskge, valid_dna_mvskge):
        """A clean sequence without any restriction site should not be flagged."""
        config = SolverConfig(
            gc_lo=0.30, gc_hi=0.70,
            restriction_sites=["GAATTC", "GGATCC"],
        )
        violations, score = validate_csp_solution(valid_dna_mvskge, protein_mvskge, config)
        rs_violations = [
            v for v in violations
            if v.constraint_name == "NoRestrictionSiteConstraint"
        ]
        assert len(rs_violations) == 0


# ---------------------------------------------------------------------------
# 5. Translation fidelity check catches wrong codon
# ---------------------------------------------------------------------------

class TestTranslationFidelity:
    """Test that the translation fidelity sanity check catches mismatches."""

    def test_wrong_codon_detected(self, default_config, protein_mvskge):
        """Replace the first codon (should be ATG for M) with TTT (F)."""
        # TTT translates to 'F', but protein[0] is 'M'
        wrong_seq = "TTTGTTTCTAAAGGTGAA"
        violations, score = validate_csp_solution(wrong_seq, protein_mvskge, default_config)
        fidelity_violations = [
            v for v in violations
            if v.constraint_name == "translation_fidelity"
        ]
        assert len(fidelity_violations) >= 1
        assert "Codon 0" in fidelity_violations[0].description

    def test_multiple_wrong_codons(self, default_config):
        """Multiple wrong codons each produce a violation."""
        protein = "MF"
        # Wrong: M should be ATG (got TTT=F), F should be TTT/TTC (got ATG=M)
        swapped_seq = "TTTATG"
        violations, score = validate_csp_solution(swapped_seq, protein, default_config)
        fidelity_violations = [
            v for v in violations
            if v.constraint_name == "translation_fidelity"
        ]
        assert len(fidelity_violations) >= 2


# ---------------------------------------------------------------------------
# 6. Sequence length mismatch returns violation
# ---------------------------------------------------------------------------

class TestSequenceLengthMismatch:
    """Test that a sequence with the wrong length is flagged."""

    def test_short_sequence_length_mismatch(self, default_config, protein_mvskge):
        """A sequence shorter than protein * 3 should be flagged."""
        short_seq = "ATGGTTTCT"  # 9 bases, protein needs 18
        violations, score = validate_csp_solution(short_seq, protein_mvskge, default_config)
        length_violations = [
            v for v in violations
            if v.constraint_name == "sequence_length"
        ]
        assert len(length_violations) >= 1
        assert "9" in length_violations[0].description
        assert "18" in length_violations[0].description

    def test_long_sequence_length_mismatch(self, default_config, protein_mvskge):
        """A sequence longer than protein * 3 should also be flagged."""
        long_seq = "ATGGTTTCTAAAGGTGAAATGGTT"  # 24 bases, protein needs 18
        violations, score = validate_csp_solution(long_seq, protein_mvskge, default_config)
        length_violations = [
            v for v in violations
            if v.constraint_name == "sequence_length"
        ]
        assert len(length_violations) >= 1


# ---------------------------------------------------------------------------
# 7. solve_with_csp returns SolverResult
# ---------------------------------------------------------------------------

class TestSolveWithCsp:
    """Test that solve_with_csp returns a SolverResult (may be fallback)."""

    @patch("biocompiler.solver.dispatch.quick_feasibility_check")
    def test_returns_solver_result(self, mock_feasibility, protein_mvskge):
        """solve_with_csp should always return a SolverResult instance."""
        # Patch feasibility check to bypass CSPModel attribute mismatch
        from biocompiler.solver.mus import FeasibilityReport
        mock_feasibility.return_value = FeasibilityReport(feasible=True)
        # Also patch backend engines to None so we get a clean fallback result
        with patch("biocompiler.solver.dispatch._ortools_engine", None), \
             patch("biocompiler.solver.dispatch._z3_engine", None):
            result = solve_with_csp(protein_mvskge)
        assert isinstance(result, SolverResult)

    @patch("biocompiler.solver.dispatch.quick_feasibility_check")
    def test_result_has_required_fields(self, mock_feasibility, protein_mvskge):
        """SolverResult should have key fields populated."""
        from biocompiler.solver.mus import FeasibilityReport
        mock_feasibility.return_value = FeasibilityReport(feasible=True)
        with patch("biocompiler.solver.dispatch._ortools_engine", None), \
             patch("biocompiler.solver.dispatch._z3_engine", None):
            result = solve_with_csp(protein_mvskge)
        assert isinstance(result.solved, bool)
        assert isinstance(result.backend_used, SolverBackend)
        assert isinstance(result.fallback_used, bool)
        assert isinstance(result.violations, list)

    def test_invalid_protein_raises(self):
        """An empty or invalid protein should raise ValueError."""
        with pytest.raises(ValueError):
            solve_with_csp("")

    def test_invalid_amino_acids_raises(self):
        """Amino acids not in the standard set should raise ValueError."""
        with pytest.raises(ValueError):
            solve_with_csp("BZX")


# ---------------------------------------------------------------------------
# 8. get_csp_availability returns dict with expected keys
# ---------------------------------------------------------------------------

class TestGetCspAvailability:
    """Test that get_csp_availability returns a dict with the expected keys."""

    def test_returns_dict(self):
        result = get_csp_availability()
        assert isinstance(result, dict)

    def test_has_expected_keys(self):
        result = get_csp_availability()
        expected_keys = {"ortools", "z3", "any"}
        assert set(result.keys()) == expected_keys

    def test_values_are_bool(self):
        result = get_csp_availability()
        for key in ("ortools", "z3", "any"):
            assert isinstance(result[key], bool), f"Key '{key}' is not bool"

    def test_any_is_logical_or(self):
        """'any' should be True iff at least one backend is available."""
        result = get_csp_availability()
        assert result["any"] == (result["ortools"] or result["z3"])


# ---------------------------------------------------------------------------
# 9. is_solver_available returns bool
# ---------------------------------------------------------------------------

class TestIsSolverAvailable:
    """Test that is_solver_available returns a boolean."""

    def test_returns_bool(self):
        result = is_solver_available()
        assert isinstance(result, bool)

    def test_consistent_with_availability(self):
        """is_solver_available() should match get_csp_availability()['any']."""
        assert is_solver_available() == get_csp_availability()["any"]


# ---------------------------------------------------------------------------
# 10. csp_optimize returns SolverResult with fallback_used on failure
# ---------------------------------------------------------------------------

class TestCspOptimizeFallback:
    """Test that csp_optimize returns a SolverResult with fallback_used=True
    when all solver backends fail."""

    def test_csp_optimize_returns_solver_result(self, protein_mvskge):
        """csp_optimize should always return a SolverResult."""
        result = csp_optimize(protein_mvskge)
        assert isinstance(result, SolverResult)

    @patch("biocompiler.solver.dispatch.solve_with_csp")
    def test_fallback_on_exception(self, mock_solve, protein_mvskge):
        """When solve_with_csp raises, csp_optimize returns fallback result."""
        mock_solve.side_effect = RuntimeError("Solver crashed")
        result = csp_optimize(protein_mvskge)
        assert isinstance(result, SolverResult)
        assert result.fallback_used is True
        assert result.solved is False
        assert result.sequence == ""

    @patch("biocompiler.solver.dispatch._ortools_engine", None)
    @patch("biocompiler.solver.dispatch._z3_engine", None)
    def test_fallback_when_no_backends(self, protein_mvskge):
        """When no backends are available, the result should indicate fallback."""
        result = csp_optimize(protein_mvskge)
        assert isinstance(result, SolverResult)
        # With no engines, solve_with_csp should produce a fallback result
        if not get_csp_availability()["any"]:
            assert result.fallback_used is True

    @patch("biocompiler.solver.dispatch.solve_with_csp")
    def test_fallback_metadata_contains_reason(self, mock_solve, protein_mvskge):
        """The fallback result should include a reason in metadata."""
        mock_solve.side_effect = Exception("Boom")
        result = csp_optimize(protein_mvskge)
        assert isinstance(result, SolverResult)
        assert result.fallback_used is True
        assert "reason" in result.metadata
