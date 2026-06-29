"""Tests for the assembly module — AssemblyPlan, plan_golden_gate, plan_gibson."""

from __future__ import annotations

import pytest

from biocompiler.optimizer.assembly import (
    AssemblyPlan,
    plan_golden_gate,
    plan_gibson,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Short DNA sequences for testing
FRAG1 = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
FRAG2 = "TTGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
FRAG3 = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"


# ---------------------------------------------------------------------------
# AssemblyPlan dataclass
# ---------------------------------------------------------------------------

class TestAssemblyPlan:
    """Tests for the AssemblyPlan dataclass."""

    def test_construction_golden_gate(self):
        """Can construct a valid Golden Gate AssemblyPlan."""
        plan = AssemblyPlan(
            method="golden_gate",
            fragments=[FRAG1, FRAG2],
            enzymes=["BsaI"],
            overlap_sequences=["GGAG"],
            total_length=84,
        )
        assert plan.method == "golden_gate"
        assert len(plan.fragments) == 2
        assert plan.enzymes == ["BsaI"]
        assert plan.total_length == 84

    def test_construction_gibson(self):
        """Can construct a valid Gibson AssemblyPlan."""
        plan = AssemblyPlan(
            method="gibson",
            fragments=[FRAG1, FRAG2],
            enzymes=[],
            overlap_sequences=["TACT"],
            total_length=64,
        )
        assert plan.method == "gibson"
        assert plan.enzymes == []

    def test_invalid_method_raises(self):
        """Invalid assembly method raises ValueError."""
        with pytest.raises(ValueError, match="Invalid assembly method"):
            AssemblyPlan(
                method="invalid",
                fragments=[],
                enzymes=[],
                overlap_sequences=[],
                total_length=0,
            )

    def test_all_valid_methods(self):
        """All three valid methods are accepted."""
        for method in ("golden_gate", "gibson", "restriction_ligation"):
            plan = AssemblyPlan(
                method=method,
                fragments=[],
                enzymes=[],
                overlap_sequences=[],
                total_length=0,
            )
            assert plan.method == method


# ---------------------------------------------------------------------------
# plan_golden_gate
# ---------------------------------------------------------------------------

class TestPlanGoldenGate:
    """Tests for the plan_golden_gate function."""

    def test_basic_plan(self):
        """Basic Golden Gate plan with two fragments."""
        plan = plan_golden_gate([FRAG1, FRAG2])
        assert plan.method == "golden_gate"
        assert len(plan.fragments) == 2
        assert len(plan.enzymes) >= 1
        assert plan.total_length == len(FRAG1) + len(FRAG2)

    def test_default_enzymes(self):
        """Default enzymes include BsaI and BsmBI."""
        plan = plan_golden_gate([FRAG1, FRAG2])
        assert "BsaI" in plan.enzymes
        assert "BsmBI" in plan.enzymes

    def test_custom_enzymes(self):
        """Custom enzyme list is used."""
        plan = plan_golden_gate([FRAG1, FRAG2], enzymes=["BsaI"])
        assert plan.enzymes == ["BsaI"]

    def test_three_fragments(self):
        """Three fragments produce two overlaps."""
        plan = plan_golden_gate([FRAG1, FRAG2, FRAG3])
        assert len(plan.fragments) == 3
        assert len(plan.overlap_sequences) == 2

    def test_single_fragment(self):
        """Single fragment produces no overlaps."""
        plan = plan_golden_gate([FRAG1])
        assert len(plan.overlap_sequences) == 0
        assert plan.total_length == len(FRAG1)

    def test_empty_sequences_raises(self):
        """Empty sequences list raises ValueError."""
        with pytest.raises(ValueError, match="At least one"):
            plan_golden_gate([])

    def test_total_length_correct(self):
        """Total length is sum of all fragment lengths."""
        plan = plan_golden_gate([FRAG1, FRAG2, FRAG3])
        expected = len(FRAG1) + len(FRAG2) + len(FRAG3)
        assert plan.total_length == expected

    def test_overhang_assignment(self):
        """Overhangs are assigned from the MoClo standard set."""
        plan = plan_golden_gate([FRAG1, FRAG2])
        assert len(plan.overlap_sequences) == 1
        # Should be a 4bp overhang from MoClo standard
        assert len(plan.overlap_sequences[0]) == 4

    def test_many_fragments_cycle_overhangs(self):
        """More junctions than MoClo overhangs cycles through them."""
        fragments = [FRAG1] * 6  # 6 fragments = 5 junctions
        plan = plan_golden_gate(fragments)
        assert len(plan.overlap_sequences) == 5

    def test_invalid_enzyme_skipped(self):
        """Unknown enzyme names are skipped with a warning."""
        plan = plan_golden_gate([FRAG1, FRAG2], enzymes=["BsaI", "FakeEnzyme"])
        assert "BsaI" in plan.enzymes
        assert "FakeEnzyme" not in plan.enzymes

    def test_all_invalid_enzymes_raises(self):
        """All invalid enzymes raises ValueError."""
        with pytest.raises(ValueError, match="No valid enzymes"):
            plan_golden_gate([FRAG1, FRAG2], enzymes=["Fake1", "Fake2"])

    def test_returns_assembly_plan(self):
        """Return type is AssemblyPlan."""
        plan = plan_golden_gate([FRAG1, FRAG2])
        assert isinstance(plan, AssemblyPlan)


# ---------------------------------------------------------------------------
# plan_gibson
# ---------------------------------------------------------------------------

class TestPlanGibson:
    """Tests for the plan_gibson function."""

    def test_basic_plan(self):
        """Basic Gibson plan with two fragments."""
        plan = plan_gibson([FRAG1, FRAG2])
        assert plan.method == "gibson"
        assert len(plan.fragments) == 2
        assert plan.enzymes == []  # Gibson does not use enzymes

    def test_default_overlap_length(self):
        """Default overlap length is 20 bp."""
        plan = plan_gibson([FRAG1, FRAG2])
        for overlap in plan.overlap_sequences:
            assert len(overlap) == 20

    def test_custom_overlap_length(self):
        """Custom overlap length is used."""
        plan = plan_gibson([FRAG1, FRAG2], overlap_length=30)
        for overlap in plan.overlap_sequences:
            assert len(overlap) == 30

    def test_overlap_sequences_derived_from_fragments(self):
        """Overlap sequences come from fragment ends."""
        plan = plan_gibson([FRAG1, FRAG2], overlap_length=10)
        # First overlap should be last 10bp of FRAG1
        assert plan.overlap_sequences[0] == FRAG1[-10:].upper()

    def test_three_fragments_two_overlaps(self):
        """Three fragments produce two overlaps."""
        plan = plan_gibson([FRAG1, FRAG2, FRAG3])
        assert len(plan.overlap_sequences) == 2

    def test_single_fragment_no_overlaps(self):
        """Single fragment produces no overlaps."""
        plan = plan_gibson([FRAG1])
        assert len(plan.overlap_sequences) == 0

    def test_total_length_subtracts_overlaps(self):
        """Total length subtracts overlapping regions."""
        plan = plan_gibson([FRAG1, FRAG2], overlap_length=20)
        expected = len(FRAG1) + len(FRAG2) - 20
        assert plan.total_length == expected

    def test_total_length_three_fragments(self):
        """Total length for 3 fragments subtracts 2 overlaps."""
        overlap = 15
        plan = plan_gibson([FRAG1, FRAG2, FRAG3], overlap_length=overlap)
        expected = len(FRAG1) + len(FRAG2) + len(FRAG3) - 2 * overlap
        assert plan.total_length == expected

    def test_empty_sequences_raises(self):
        """Empty sequences list raises ValueError."""
        with pytest.raises(ValueError, match="At least one"):
            plan_gibson([])

    def test_overlap_too_short_raises(self):
        """Overlap length < 4 raises ValueError."""
        with pytest.raises(ValueError, match="overlap_length"):
            plan_gibson([FRAG1, FRAG2], overlap_length=2)

    def test_overlap_length_4_minimum(self):
        """Overlap length of 4 is the minimum accepted."""
        plan = plan_gibson([FRAG1, FRAG2], overlap_length=4)
        assert len(plan.overlap_sequences[0]) == 4

    def test_short_fragment_uses_full_as_overlap(self):
        """Fragment shorter than overlap_length uses full fragment as overlap."""
        short_frag = "ATGC"
        plan = plan_gibson([short_frag, FRAG2], overlap_length=20)
        # Short fragment is only 4bp, so overlap should be the full fragment
        assert plan.overlap_sequences[0] == short_frag.upper()

    def test_returns_assembly_plan(self):
        """Return type is AssemblyPlan."""
        plan = plan_gibson([FRAG1, FRAG2])
        assert isinstance(plan, AssemblyPlan)

    def test_case_insensitive_sequences(self):
        """Sequences are uppercased in overlap computation."""
        plan = plan_gibson([FRAG1.lower(), FRAG2.lower()], overlap_length=10)
        assert plan.overlap_sequences[0] == FRAG1[-10:].upper()
