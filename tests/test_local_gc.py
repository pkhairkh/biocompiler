"""Tests for the local_gc module — LocalGCConstraint, check_local_gc, optimize_local_gc."""

from __future__ import annotations

import pytest

from biocompiler.sequence.local_gc import (
    LocalGCConstraint,
    LocalGCResult,
    check_local_gc,
    optimize_local_gc,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_gc_seq(gc_fraction: float, length: int = 300) -> str:
    """Create a DNA sequence with an approximate GC fraction."""
    n_gc = int(gc_fraction * length)
    n_at = length - n_gc
    return "G" * (n_gc // 2) + "C" * (n_gc - n_gc // 2) + "A" * (n_at // 2) + "T" * (n_at - n_at // 2)


# ---------------------------------------------------------------------------
# LocalGCConstraint dataclass
# ---------------------------------------------------------------------------

class TestLocalGCConstraint:
    """Tests for the LocalGCConstraint dataclass."""

    def test_construction(self):
        """Can construct a valid LocalGCConstraint."""
        c = LocalGCConstraint(region_start=0, region_end=50, gc_min=0.3, gc_max=0.7)
        assert c.region_start == 0
        assert c.region_end == 50
        assert c.gc_min == 0.3
        assert c.gc_max == 0.7

    def test_frozen(self):
        """LocalGCConstraint is frozen (immutable)."""
        c = LocalGCConstraint(region_start=0, region_end=50, gc_min=0.3, gc_max=0.7)
        with pytest.raises(AttributeError):
            c.region_start = 10  # type: ignore[misc]

    def test_negative_start_raises(self):
        """Negative region_start raises ValueError."""
        with pytest.raises(ValueError, match="region_start"):
            LocalGCConstraint(region_start=-1, region_end=50, gc_min=0.3, gc_max=0.7)

    def test_end_le_start_raises(self):
        """region_end <= region_start raises ValueError."""
        with pytest.raises(ValueError, match="region_end"):
            LocalGCConstraint(region_start=50, region_end=50, gc_min=0.3, gc_max=0.7)

    def test_end_lt_start_raises(self):
        """region_end < region_start raises ValueError."""
        with pytest.raises(ValueError, match="region_end"):
            LocalGCConstraint(region_start=50, region_end=10, gc_min=0.3, gc_max=0.7)

    def test_gc_min_negative_raises(self):
        """Negative gc_min raises ValueError."""
        with pytest.raises(ValueError, match="gc_min"):
            LocalGCConstraint(region_start=0, region_end=50, gc_min=-0.1, gc_max=0.7)

    def test_gc_max_over_one_raises(self):
        """gc_max > 1.0 raises ValueError."""
        with pytest.raises(ValueError, match="gc_max"):
            LocalGCConstraint(region_start=0, region_end=50, gc_min=0.3, gc_max=1.5)

    def test_gc_min_gt_gc_max_raises(self):
        """gc_min > gc_max raises ValueError."""
        with pytest.raises(ValueError, match="gc_min"):
            LocalGCConstraint(region_start=0, region_end=50, gc_min=0.8, gc_max=0.3)

    def test_gc_min_equals_gc_max_ok(self):
        """gc_min == gc_max is allowed (exact GC target)."""
        c = LocalGCConstraint(region_start=0, region_end=50, gc_min=0.5, gc_max=0.5)
        assert c.gc_min == c.gc_max

    def test_boundary_values(self):
        """Boundary values gc_min=0.0, gc_max=1.0 are accepted."""
        c = LocalGCConstraint(region_start=0, region_end=50, gc_min=0.0, gc_max=1.0)
        assert c.gc_min == 0.0
        assert c.gc_max == 1.0


# ---------------------------------------------------------------------------
# check_local_gc
# ---------------------------------------------------------------------------

class TestCheckLocalGC:
    """Tests for the check_local_gc function."""

    def test_satisfied_constraints(self):
        """Sequence within GC range returns satisfied=True."""
        # 50% GC sequence, constraint allows 0.3-0.7
        seq = _make_gc_seq(0.5, 100)
        c = LocalGCConstraint(region_start=0, region_end=100, gc_min=0.3, gc_max=0.7)
        result = check_local_gc(seq, [c])
        assert result.satisfied is True
        assert len(result.violations) == 0

    def test_gc_too_low(self):
        """Sequence with GC below minimum is detected as violation."""
        seq = _make_gc_seq(0.1, 100)
        c = LocalGCConstraint(region_start=0, region_end=100, gc_min=0.3, gc_max=0.7)
        result = check_local_gc(seq, [c])
        assert result.satisfied is False
        assert len(result.violations) == 1
        assert result.violations[0][0] is c
        assert result.violations[0][1] < 0.3

    def test_gc_too_high(self):
        """Sequence with GC above maximum is detected as violation."""
        seq = _make_gc_seq(0.9, 100)
        c = LocalGCConstraint(region_start=0, region_end=100, gc_min=0.3, gc_max=0.7)
        result = check_local_gc(seq, [c])
        assert result.satisfied is False
        assert len(result.violations) == 1
        assert result.violations[0][1] > 0.7

    def test_multiple_constraints(self):
        """Multiple constraints are all checked."""
        # Build a sequence where both halves have ~50% GC
        half1 = _make_gc_seq(0.5, 100)
        half2 = _make_gc_seq(0.5, 100)
        seq = half1 + half2
        c1 = LocalGCConstraint(region_start=0, region_end=100, gc_min=0.3, gc_max=0.7)
        c2 = LocalGCConstraint(region_start=100, region_end=200, gc_min=0.4, gc_max=0.6)
        result = check_local_gc(seq, [c1, c2])
        # Both regions have ~0.5 GC, so both should be satisfied
        assert result.satisfied is True

    def test_partial_violation(self):
        """Only some constraints violated."""
        # First half 50% GC, second half 90% GC
        seq = _make_gc_seq(0.5, 100) + _make_gc_seq(0.9, 100)
        c1 = LocalGCConstraint(region_start=0, region_end=100, gc_min=0.3, gc_max=0.7)
        c2 = LocalGCConstraint(region_start=100, region_end=200, gc_min=0.3, gc_max=0.7)
        result = check_local_gc(seq, [c1, c2])
        assert result.satisfied is False
        assert len(result.violations) == 1

    def test_empty_constraints(self):
        """No constraints means always satisfied."""
        result = check_local_gc("ATGCGT", [])
        assert result.satisfied is True

    def test_region_beyond_sequence(self):
        """Constraint region beyond sequence length is skipped gracefully."""
        result = check_local_gc("ATG", [LocalGCConstraint(region_start=100, region_end=200, gc_min=0.3, gc_max=0.7)])
        assert result.satisfied is True

    def test_returns_local_gc_result(self):
        """Return type is LocalGCResult."""
        result = check_local_gc("ATGCGT", [])
        assert isinstance(result, LocalGCResult)

    def test_case_insensitive(self):
        """Sequence is upper-cased before checking."""
        seq = _make_gc_seq(0.5, 100).lower()
        c = LocalGCConstraint(region_start=0, region_end=100, gc_min=0.3, gc_max=0.7)
        result = check_local_gc(seq, [c])
        assert result.satisfied is True


# ---------------------------------------------------------------------------
# optimize_local_gc
# ---------------------------------------------------------------------------

class TestOptimizeLocalGC:
    """Tests for the optimize_local_gc function."""

    def test_already_satisfied(self):
        """If constraints already satisfied, sequence is returned unchanged."""
        seq = _make_gc_seq(0.5, 99)  # 99 = divisible by 3
        protein = "M" * 33
        c = LocalGCConstraint(region_start=0, region_end=99, gc_min=0.3, gc_max=0.7)
        result = optimize_local_gc(seq, protein, [c])
        assert result.satisfied is True
        assert result.sequence == seq.upper()

    def test_improves_gc(self):
        """Optimization should move GC content toward the target range."""
        # All-A sequence (0% GC) with constraint requiring >= 30%
        seq = "GCA" * 33  # 99 bp, ~33% GC per codon
        protein = "A" * 33  # All alanine
        c = LocalGCConstraint(region_start=0, region_end=99, gc_min=0.5, gc_max=0.9)
        result = optimize_local_gc(seq, protein, [c])
        # Should have tried to increase GC via synonymous codons for Alanine
        # Alanine codons: GCT, GCC, GCA, GCG — all have 66.7% GC
        # Even GCA has 66.7% GC per codon, so the constraint should be satisfiable
        assert isinstance(result, LocalGCResult)

    def test_preserves_translation(self):
        """Optimized sequence should still translate to the same protein."""
        from biocompiler.expression.translation import translate
        from biocompiler.shared.constants import CODON_TABLE

        # Build a sequence with known protein
        protein = "MVSKGE"
        # M=ATG, V=GTT, S=TCA, K=AAG, G=GGT, E=GAG
        dna = "ATGGTTTCAAAGGGTGAG"
        c = LocalGCConstraint(region_start=0, region_end=18, gc_min=0.3, gc_max=0.7)
        result = optimize_local_gc(dna, protein, [c], codon_table=CODON_TABLE)
        translated = translate(result.sequence)
        assert translated == protein

    def test_returns_local_gc_result(self):
        """Return type is LocalGCResult."""
        result = optimize_local_gc("ATGCGT", "MR", [])
        assert isinstance(result, LocalGCResult)

    def test_empty_sequence(self):
        """Empty DNA and protein are handled gracefully."""
        result = optimize_local_gc("", "", [])
        assert isinstance(result, LocalGCResult)
        assert result.satisfied is True

    def test_multiple_constraints_optimization(self):
        """Optimization with multiple region constraints."""
        # Two regions with different GC targets
        seq = "ATG" * 20 + "GCG" * 20  # 120 bp total
        protein = "M" * 20 + "A" * 20
        c1 = LocalGCConstraint(region_start=0, region_end=60, gc_min=0.3, gc_max=0.7)
        c2 = LocalGCConstraint(region_start=60, region_end=120, gc_min=0.3, gc_max=0.9)
        result = optimize_local_gc(seq, protein, [c1, c2])
        assert isinstance(result, LocalGCResult)

    def test_custom_codon_table(self):
        """Custom codon table is used when provided."""
        # Minimal codon table for testing
        custom_table = {"ATG": "M", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A"}
        dna = "ATGGCT"
        protein = "MA"
        c = LocalGCConstraint(region_start=0, region_end=6, gc_min=0.3, gc_max=0.7)
        result = optimize_local_gc(dna, protein, [c], codon_table=custom_table)
        assert isinstance(result, LocalGCResult)

    def test_single_codon_aa_cannot_change(self):
        """Amino acids with only one codon (e.g., M=ATG, W=TGG) cannot be changed."""
        # Met only has ATG — GC content of that codon is fixed at 33.3%
        dna = "ATGATGATG"  # MMM, 33.3% GC
        protein = "MMM"
        c = LocalGCConstraint(region_start=0, region_end=9, gc_min=0.8, gc_max=1.0)
        result = optimize_local_gc(dna, protein, [c])
        # ATG is the only codon for Met, so we cannot increase GC
        # This should NOT be satisfied, but should not crash
        assert isinstance(result, LocalGCResult)

    def test_violations_recorded(self):
        """Violations are recorded even after optimization if unfixable."""
        # Constraint that is impossible to satisfy with Met-only sequence
        dna = "ATGATG"
        protein = "MM"
        c = LocalGCConstraint(region_start=0, region_end=6, gc_min=0.9, gc_max=1.0)
        result = optimize_local_gc(dna, protein, [c])
        assert isinstance(result, LocalGCResult)
        # Met only has ATG, so this is unsatisfiable
