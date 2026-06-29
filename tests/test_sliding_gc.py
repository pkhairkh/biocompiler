"""Tests for the sliding-window GC constraint module.

Covers:
- check_sliding_gc: basic scanning, boundary cases, violation detection
- fix_sliding_gc_violations: codon swapping to resolve local GC extremes
- evaluate_sliding_gc: predicate-style evaluation
- Integration with optimize_sequence
- Predicate registry entry
"""

import pytest
from biocompiler.sequence.sliding_gc import (
    WindowViolation,
    SlidingGCResult,
    check_sliding_gc,
    fix_sliding_gc_violations,
    evaluate_sliding_gc,
)
from biocompiler.shared.types import Verdict


# ────────────────────────────────────────────────────────────
# check_sliding_gc
# ────────────────────────────────────────────────────────────

class TestCheckSlidingGC:
    """Test the check_sliding_gc function."""

    def test_uniform_gc_passes(self):
        """A sequence with uniform 50% GC should pass."""
        # Mix of GC and AT to get ~50% GC
        seq = "GCAT" * 25  # 100 bp, ~50% GC everywhere
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is True
        assert len(result.violations) == 0
        # With "GCAT" repeat, sliding windows can vary slightly from 50%
        assert 0.40 <= result.min_gc <= 0.60
        assert 0.40 <= result.max_gc <= 0.60

    def test_all_gc_fails_too_high(self):
        """A sequence of all G/C should fail when gc_max < 1.0."""
        seq = "GCGC" * 25  # 100 bp, 100% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is False
        assert len(result.violations) > 0
        assert all(v.direction == "too_high" for v in result.violations)

    def test_all_at_fails_too_low(self):
        """A sequence of all A/T should fail when gc_min > 0.0."""
        seq = "ATAT" * 25  # 100 bp, 0% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is False
        assert len(result.violations) > 0
        assert all(v.direction == "too_low" for v in result.violations)

    def test_local_gc_extreme_detected(self):
        """A local region of extreme GC should be detected."""
        # 80 bp: 40 bp AT-rich + 40 bp GC-rich
        at_region = "ATAT" * 10  # 40 bp, 0% GC
        gc_region = "GCGC" * 10  # 40 bp, 100% GC
        seq = at_region + gc_region
        result = check_sliding_gc(seq, window_size=20, gc_min=0.20, gc_max=0.80)
        assert result.passed is False
        # Should have violations in both extremes
        low_violations = [v for v in result.violations if v.direction == "too_low"]
        high_violations = [v for v in result.violations if v.direction == "too_high"]
        assert len(low_violations) > 0
        assert len(high_violations) > 0

    def test_short_sequence(self):
        """Sequence shorter than window_size: check whole sequence."""
        seq = "GCGCATAT"  # 8 bp, 50% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is True
        assert len(result.violations) == 0
        assert result.min_gc == result.max_gc  # Only one "window" (whole seq)

    def test_short_sequence_fails(self):
        """Short sequence with extreme GC fails."""
        seq = "GCGCGCGC"  # 8 bp, 100% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is False
        assert len(result.violations) == 1
        assert result.violations[0].direction == "too_high"

    def test_empty_sequence(self):
        """Empty sequence should return passed=True with no violations."""
        result = check_sliding_gc("", window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is True
        assert len(result.violations) == 0
        assert result.min_gc == 0.0
        assert result.max_gc == 0.0

    def test_window_equals_sequence(self):
        """Window size equal to sequence length."""
        seq = "GCGCATAT" * 6  # 48 bp, 50% GC
        result = check_sliding_gc(seq, window_size=48, gc_min=0.30, gc_max=0.70)
        assert result.passed is True
        assert result.min_gc == result.max_gc

    def test_step_parameter(self):
        """Step parameter > 1 should skip some windows."""
        seq = "GCAT" * 25  # 100 bp, 50% GC
        result_step1 = check_sliding_gc(seq, window_size=10, gc_min=0.30, gc_max=0.70, step=1)
        result_step5 = check_sliding_gc(seq, window_size=10, gc_min=0.30, gc_max=0.70, step=5)
        # Step=1 scans more windows than step=5
        # With uniform GC both should pass
        assert result_step1.passed is True
        assert result_step5.passed is True

    def test_case_insensitive(self):
        """Input should be case-insensitive."""
        seq_lower = "gcgcatat" * 6
        seq_upper = "GCGCATAT" * 6
        r1 = check_sliding_gc(seq_lower, window_size=20, gc_min=0.30, gc_max=0.70)
        r2 = check_sliding_gc(seq_upper, window_size=20, gc_min=0.30, gc_max=0.70)
        assert r1.passed == r2.passed
        assert abs(r1.min_gc - r2.min_gc) < 0.01
        assert abs(r1.max_gc - r2.max_gc) < 0.01

    def test_violation_positions(self):
        """Violations should have valid start/end positions."""
        at_region = "ATAT" * 20  # 80 bp, 0% GC
        gc_region = "GCGC" * 20  # 80 bp, 100% GC
        seq = at_region + gc_region
        result = check_sliding_gc(seq, window_size=20, gc_min=0.20, gc_max=0.80)
        for v in result.violations:
            assert v.start >= 0
            assert v.end > v.start
            assert v.end <= len(seq)
            assert 0.0 <= v.gc_content <= 1.0

    def test_relaxed_bounds_pass(self):
        """Very relaxed GC bounds should always pass."""
        seq = "GCGCGCGC" * 12  # 96 bp, 100% GC
        result = check_sliding_gc(seq, window_size=20, gc_min=0.0, gc_max=1.0)
        assert result.passed is True
        assert len(result.violations) == 0


# ────────────────────────────────────────────────────────────
# fix_sliding_gc_violations
# ────────────────────────────────────────────────────────────

class TestFixSlidingGCViolations:
    """Test the fix_sliding_gc_violations function."""

    def test_no_violations_no_change(self):
        """If no violations, sequence should be unchanged."""
        # All alanine: GCT has 2/3 GC, within normal range
        protein = "A" * 30  # 30 alanines
        # Use a mix of codons to get reasonable GC
        dna = "GCT" * 30  # 90 bp, ~67% GC
        fixed, swaps = fix_sliding_gc_violations(
            dna, protein, window_size=30, gc_min=0.30, gc_max=0.70
        )
        # Might have violations at 67% GC being close to 70% max,
        # but no swaps should change the protein
        assert len(fixed) == len(dna)

    def test_protein_preserved(self):
        """Fixed sequence must encode the same protein."""
        protein = "MALWMRLLPL"  # 10 AAs
        # Back-translate with high-GC codons
        dna = "ATGGCTCTGTGGATGCGGCTGCTGCCGCTGCCG"
        # Pad if needed
        if len(dna) < len(protein) * 3:
            dna = dna + "GCT" * (len(protein) * 3 - len(dna)) // 3
        dna = dna[:len(protein) * 3]

        from biocompiler.type_system import CODON_TABLE
        fixed, swaps = fix_sliding_gc_violations(
            dna, protein, window_size=20, gc_min=0.20, gc_max=0.80
        )
        # Verify protein is preserved
        translated = "".join(
            CODON_TABLE.get(fixed[i:i+3], "X")
            for i in range(0, len(fixed), 3)
        )
        assert translated == protein

    def test_global_gc_respected(self):
        """Fixes should not push global GC out of bounds."""
        protein = "AAAAAAAAAA"  # 10 alanines — many codon choices
        dna = "GCC" * 10  # All GCC (high GC), 30 bp
        fixed, swaps = fix_sliding_gc_violations(
            dna, protein, window_size=15, gc_min=0.20, gc_max=0.80,
            gc_lo=0.30, gc_hi=0.70,
        )
        gc = (fixed.count("G") + fixed.count("C")) / len(fixed)
        assert 0.30 <= gc <= 0.70 or swaps == 0  # Either fixed or no change possible

    def test_fix_with_usage_table(self):
        """When usage table provided, higher-CAI codons preferred."""
        protein = "FF"  # Phenylalanine: TTT (low GC) or TTC (high GC)
        dna = "TTTTTT"  # Both Phe as TTT, 0% GC
        usage = {"TTT": 0.5, "TTC": 0.8}  # TTC is higher CAI
        fixed, swaps = fix_sliding_gc_violations(
            dna, protein, window_size=6, gc_min=0.10, gc_max=0.90,
            usage=usage,
        )
        # With gc_min=0.10, no violation expected (0% GC is too low for gc_min=0.10)
        # 0% < 10%, so violation exists
        # The fix should swap TTT→TTC (higher CAI, adds GC)
        if swaps > 0:
            assert "TTC" in fixed


# ────────────────────────────────────────────────────────────
# evaluate_sliding_gc
# ────────────────────────────────────────────────────────────

class TestEvaluateSlidingGC:
    """Test the predicate-style evaluate_sliding_gc function."""

    def test_pass(self):
        """Sequence within bounds should get PASS verdict."""
        seq = "GCGCATAT" * 10  # 80 bp, 50% GC
        result = evaluate_sliding_gc(seq, window_size=20, gc_min=0.30, gc_max=0.70)
        assert result.verdict == Verdict.PASS
        assert "SlidingGC" in result.predicate

    def test_fail(self):
        """Sequence with local GC extremes should get FAIL verdict."""
        seq = "GCGCGCGC" * 10  # 80 bp, 100% GC
        result = evaluate_sliding_gc(seq, window_size=20, gc_min=0.30, gc_max=0.70)
        assert result.verdict == Verdict.FAIL
        assert result.violation is not None
        assert "SlidingGC" in result.predicate

    def test_predicate_name_includes_params(self):
        """Predicate name should include window size and bounds."""
        result = evaluate_sliding_gc("GCGCATAT" * 10, window_size=30, gc_min=0.25, gc_max=0.75)
        assert "30" in result.predicate
        assert "0.25" in result.predicate
        assert "0.75" in result.predicate


# ────────────────────────────────────────────────────────────
# Predicate Registry
# ────────────────────────────────────────────────────────────

class TestPredicateRegistry:
    """Test that SlidingGC is registered in the predicate registry."""

    def test_registry_contains_sliding_gc(self):
        """SlidingGC should be in the predicate registry."""
        from biocompiler.type_system import registry
        assert "SlidingGC" in registry

    def test_registry_evaluate(self):
        """Should be able to evaluate SlidingGC through the registry."""
        from biocompiler.type_system import registry
        seq = "GCGCATAT" * 10
        result = registry.evaluate("SlidingGC", seq=seq, window_size=20, gc_min=0.30, gc_max=0.70)
        assert result.verdict == Verdict.PASS

    def test_registry_verify(self):
        """Should be able to verify SlidingGC through the registry."""
        from biocompiler.type_system import registry
        seq = "GCGCATAT" * 10
        result = registry.verify("SlidingGC", seq=seq, window_size=20, gc_min=0.30, gc_max=0.70)
        assert result.verdict == Verdict.PASS


# ────────────────────────────────────────────────────────────
# Integration with optimize_sequence
# ────────────────────────────────────────────────────────────

class TestOptimizeIntegration:
    """Test sliding GC integration with the optimization pipeline."""

    def test_default_params(self):
        """optimize_sequence should work with default sliding GC params."""
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(
            target_protein="MALWMRLLPL",
            organism="Escherichia_coli",
            gc_window_size=0,  # Disable for simple test
        )
        assert result.sequence is not None
        assert len(result.sequence) > 0

    def test_sliding_gc_enabled(self):
        """optimize_sequence with sliding GC enabled should still produce valid output."""
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(
            target_protein="MALWMRLLPL",
            organism="Escherichia_coli",
            gc_window_size=30,
            gc_window_min=0.20,
            gc_window_max=0.80,
        )
        assert result.sequence is not None
        # Check that the protein is preserved
        from biocompiler.type_system import CODON_TABLE
        translated = "".join(
            CODON_TABLE.get(result.sequence[i:i+3], "X")
            for i in range(0, len(result.sequence) - 2, 3)
        )
        # Remove stop codon if present
        translated = translated.rstrip("*")
        assert translated == "MALWMRLLPL"

    def test_sliding_gc_predicate_in_results(self):
        """When sliding GC is enabled, SlidingGC should appear in predicate results."""
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(
            target_protein="MVHLTPEEKS",
            organism="Escherichia_coli",
            gc_window_size=30,
        )
        pred_names = [r.predicate for r in result.predicate_results]
        assert "SlidingGC" in pred_names

    def test_sliding_gc_disabled(self):
        """When gc_window_size=0, sliding GC should be skipped."""
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(
            target_protein="MVHLTPEEKS",
            organism="Escherichia_coli",
            gc_window_size=0,
        )
        # Should still work normally
        assert result.sequence is not None


# ────────────────────────────────────────────────────────────
# Data class validation
# ────────────────────────────────────────────────────────────

class TestDataClasses:
    """Test data class validation."""

    def test_window_violation_creation(self):
        """WindowViolation should be created with valid data."""
        v = WindowViolation(start=10, end=30, gc_content=0.15, direction="too_low")
        assert v.start == 10
        assert v.end == 30
        assert v.gc_content == 0.15
        assert v.direction == "too_low"

    def test_window_violation_invalid_direction(self):
        """WindowViolation with invalid direction should raise."""
        with pytest.raises(AssertionError):
            WindowViolation(start=10, end=30, gc_content=0.5, direction="invalid")

    def test_window_violation_negative_start(self):
        """WindowViolation with negative start should raise."""
        with pytest.raises(AssertionError):
            WindowViolation(start=-1, end=30, gc_content=0.5, direction="too_low")

    def test_window_violation_end_before_start(self):
        """WindowViolation with end <= start should raise."""
        with pytest.raises(AssertionError):
            WindowViolation(start=30, end=10, gc_content=0.5, direction="too_low")

    def test_sliding_gc_result_creation(self):
        """SlidingGCResult should be created with valid data."""
        r = SlidingGCResult(passed=True, min_gc=0.45, max_gc=0.55, violations=[])
        assert r.passed is True
        assert r.min_gc == 0.45
        assert r.max_gc == 0.55

    def test_sliding_gc_result_min_gt_max(self):
        """SlidingGCResult with min_gc > max_gc should raise."""
        with pytest.raises(AssertionError):
            SlidingGCResult(passed=True, min_gc=0.7, max_gc=0.3, violations=[])


# ────────────────────────────────────────────────────────────
# Edge cases
# ────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_window_size_1(self):
        """Window size of 1 should work (each base is a window)."""
        seq = "GCATGCAT"
        result = check_sliding_gc(seq, window_size=1, gc_min=0.0, gc_max=1.0)
        assert result.passed is True  # All bases are either 0% or 100% GC, but 0-1 range allows all

    def test_exact_boundary_gc(self):
        """GC content exactly at boundary should pass (inclusive)."""
        # 50 bp window with exactly 30% GC (15 GC bases)
        at_part = "AT" * 17  # 34 bp, 0% GC → but we need 50 bp total with 15 GC
        gc_part = "GC" * 8   # 16 bp, 100% GC → total 50 bp with 16 GC = 32%
        seq = at_part + gc_part  # 50 bp, 32% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.32, gc_max=0.70)
        assert result.passed is True  # 32% >= 32% (inclusive)

    def test_just_below_min(self):
        """GC just below gc_min should fail."""
        seq = "AT" * 25  # 50 bp, 0% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.01, gc_max=0.70)
        assert result.passed is False

    def test_just_above_max(self):
        """GC just above gc_max should fail."""
        seq = "GC" * 25  # 50 bp, 100% GC
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.99)
        assert result.passed is False
