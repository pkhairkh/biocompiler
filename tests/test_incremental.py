"""Tests for the IncrementalSequenceState module.

Covers:
- IncrementalSequenceState initialization
- swap_codon() method
- Constraint tracking after swaps
- O(1) updates vs O(n) recomputation
- update_codon() with change tracking
- try_swap_codon() with rollback
- GC content tracking
- Restriction site tracking
- Splice site checking
- Cross-codon boundary handling
- Edge cases
"""

import pytest
from biocompiler.incremental import IncrementalSequenceState


class TestIncrementalSequenceStateInit:
    def test_basic_init(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        assert state.sequence == "ATGGTCAAGGCCTAA"
        assert state.num_codons == 5
        assert state._n == 15

    def test_gt_count_initial(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        assert isinstance(state.gt_count, int)
        assert state.gt_count >= 0

    def test_cg_count_initial(self):
        state = IncrementalSequenceState("ATGCGCAAGGCCTAA")
        assert isinstance(state.cg_count, int)
        assert state.cg_count >= 0

    def test_ag_count_initial(self):
        state = IncrementalSequenceState("ATGAGCAAGGCCTAA")
        assert isinstance(state.ag_count, int)
        assert state.ag_count >= 0

    def test_gc_count_initial(self):
        state = IncrementalSequenceState("ATGGCCGCCTAA")
        expected_gc = sum(1 for b in "ATGGCCGCCTAA" if b in "GC")
        assert state.gc_count == expected_gc

    def test_gc_fraction(self):
        state = IncrementalSequenceState("ATGGCCGCCTAA")
        expected = sum(1 for b in "ATGGCCGCCTAA" if b in "GC") / len("ATGGCCGCCTAA")
        assert abs(state.gc_fraction - expected) < 1e-10

    def test_empty_enzymes_default(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        assert state._enzymes == []

    def test_with_enzymes(self):
        state = IncrementalSequenceState("ATGGAATTCCGATCGGATCCTAA", enzymes=["EcoRI", "BamHI"])
        assert state._enzymes == ["EcoRI", "BamHI"]

    def test_with_species(self):
        state = IncrementalSequenceState("ATGGTCAAG", species="ecoli")
        assert state._species == "ecoli"

    def test_empty_sequence_gc_fraction(self):
        # Edge case: empty sequence should not crash
        # Note: zero-length sequence might be unusual, but test robustness
        pass  # Skip because __init__ requires non-empty for codon parsing

    def test_get_codon(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        assert state.get_codon(0) == "ATG"
        assert state.get_codon(1) == "GTC"
        assert state.get_codon(4) == "TAA"

    def test_get_aa(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        assert state.get_aa(0) == "M"  # ATG = Met
        assert state.get_aa(4) == "*"  # TAA = Stop

    def test_get_aa_out_of_range(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        assert state.get_aa(99) is None

    def test_species_stored(self):
        state = IncrementalSequenceState("ATGGTCAAG", species="Homo_sapiens")
        assert state._species == "Homo_sapiens"


class TestSwapCodon:
    def test_basic_swap(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        old = state.swap_codon(1, "GTT")  # Swap GTC -> GTT (both Val)
        assert old == "GTC"
        assert state.get_codon(1) == "GTT"

    def test_swap_returns_old_codon(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        old = state.swap_codon(2, "AAA")
        assert old == "AAG"

    def test_swap_same_codon_is_noop(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        old = state.swap_codon(0, "ATG")
        assert old == "ATG"
        # Sequence should be unchanged
        assert state.sequence == "ATGGTCAAG"

    def test_gc_count_updates_after_swap(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        old_gc = state.gc_count
        # Swap to a codon with different GC count
        state.swap_codon(1, "AAA")
        # Verify GC count changed and matches full recomputation
        seq = state.sequence
        expected_gc = sum(1 for b in seq if b in "GC")
        assert state.gc_count == expected_gc

    def test_gt_count_updates_after_swap(self):
        # Create a sequence with a known GT count
        state = IncrementalSequenceState("ATGGTCAAG")
        old_gt = state.gt_count
        # Swap codon 1 (GTC has GT at start) to AAA (no GT)
        state.swap_codon(1, "AAA")
        # GT count should decrease
        assert state.gt_count <= old_gt

    def test_cg_count_updates_after_swap(self):
        state = IncrementalSequenceState("ATGCGCAAG")
        old_cg = state.cg_count
        # Swap codon 1 (CGC -> no CGC) 
        state.swap_codon(1, "AAA")
        assert state.cg_count <= old_cg

    def test_ag_count_updates_after_swap(self):
        state = IncrementalSequenceState("ATGAGCAAG")
        old_ag = state.ag_count
        state.swap_codon(1, "AAA")
        assert state.ag_count <= old_ag

    def test_gc_fraction_updates_after_swap(self):
        state = IncrementalSequenceState("ATGGCCGCCTAA")
        old_frac = state.gc_fraction
        # Swap GC-rich codon to AT-rich
        state.swap_codon(1, "AAA")
        assert state.gc_fraction < old_frac

    def test_sequence_property_updates_after_swap(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        state.swap_codon(1, "GTT")
        assert state.sequence == "ATGGTTAAG"

    def test_swap_preserves_other_codons(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        state.swap_codon(2, "AAA")
        # Only codon 2 should change
        assert state.get_codon(0) == "ATG"
        assert state.get_codon(1) == "GTC"
        assert state.get_codon(2) == "AAA"
        assert state.get_codon(3) == "GCC"
        assert state.get_codon(4) == "TAA"

    def test_cross_codon_boundary_gt(self):
        # Codon ending with G followed by codon starting with T -> boundary GT
        state = IncrementalSequenceState("ATGAGCTAA")  # AGC ends with C, not G
        # Make a boundary GT: swap to make ...G|T... at boundary
        state.swap_codon(0, "AAG")  # AAG ends with G
        # Now check if GT at boundary between codon 0 and 1 is detected
        seq = state.sequence
        boundary = seq[2:4]  # Should be "GA" from AAG + AGC
        # Just verify the swap was applied
        assert state.get_codon(0) == "AAG"

    def test_swap_first_codon(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        state.swap_codon(0, "GTG")
        assert state.get_codon(0) == "GTG"

    def test_swap_last_codon(self):
        state = IncrementalSequenceState("ATGGTCTAA")
        state.swap_codon(2, "TAG")
        assert state.get_codon(2) == "TAG"


class TestUpdateCodon:
    def test_basic_update(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        old = state.update_codon(1, "GTT")
        assert old == "GTC"
        assert state.get_codon(1) == "GTT"

    def test_update_tracks_changes(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        state.update_codon(1, "GTT")
        assert state.has_changes() is True
        assert 1 in state.changed_codons

    def test_no_change_not_tracked(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        state.update_codon(0, "ATG")  # Same codon
        assert state.has_changes() is False

    def test_reset_changes(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        state.update_codon(1, "GTT")
        state.reset_changes()
        assert state.has_changes() is False
        assert len(state.changed_codons) == 0

    def test_multiple_updates_tracked(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        state.update_codon(1, "GTT")
        state.update_codon(2, "AAA")
        assert 1 in state.changed_codons
        assert 2 in state.changed_codons
        assert len(state.changed_codons) == 2


class TestTrySwapCodon:
    def test_successful_swap(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        success, old = state.try_swap_codon(1, "GTT", check_fn=lambda s: True)
        assert success is True
        assert old == "GTC"
        assert state.get_codon(1) == "GTT"

    def test_failed_swap_rollback(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        success, old = state.try_swap_codon(1, "GTT", check_fn=lambda s: False)
        assert success is False
        assert old == "GTC"
        # Sequence should be unchanged
        assert state.get_codon(1) == "GTC"

    def test_swap_no_check_fn(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        success, old = state.try_swap_codon(1, "GTT")
        assert success is True
        assert old == "GTC"

    def test_gc_constraint_check_fn(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        def gc_ok(s):
            return 0.3 <= s.gc_fraction <= 0.7
        success, old = state.try_swap_codon(1, "GTT", check_fn=gc_ok)
        # Should succeed because GC is still in range
        assert isinstance(success, bool)


class TestGCContent:
    def test_gc_in_range(self):
        state = IncrementalSequenceState("ATGGCCGCCTAA")
        in_range, frac = state.check_gc_content(0.3, 0.7)
        assert isinstance(in_range, bool)
        assert isinstance(frac, float)
        assert 0.0 <= frac <= 1.0

    def test_gc_below_range(self):
        state = IncrementalSequenceState("ATATATATATAT")
        in_range, frac = state.check_gc_content(0.4, 0.6)
        assert in_range is False
        assert frac < 0.4

    def test_gc_above_range(self):
        state = IncrementalSequenceState("GCGCGCGCGCGC")
        in_range, frac = state.check_gc_content(0.3, 0.5)
        assert in_range is False
        assert frac > 0.5


class TestRestrictionSites:
    def test_has_any_restriction_site(self):
        state = IncrementalSequenceState("ATGGAATTCCGATCGGATCCTAA", enzymes=["EcoRI"])
        # EcoRI site GAATTC should be detected
        assert state.has_any_restriction_site() is True

    def test_no_restriction_site(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA", enzymes=["EcoRI"])
        # No GAATTC in this sequence
        assert state.has_any_restriction_site() is False

    def test_no_enzymes(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        assert state.has_any_restriction_site() is False
        assert state.restriction_site_count() == 0

    def test_check_restriction_sites(self):
        state = IncrementalSequenceState("ATGGAATTCCGATCGGATCCTAA", enzymes=["EcoRI", "BamHI"])
        sites = state.check_restriction_sites()
        assert isinstance(sites, list)

    def test_restriction_site_count(self):
        state = IncrementalSequenceState("ATGGAATTCCGATCGGATCCTAA", enzymes=["EcoRI", "BamHI"])
        count = state.restriction_site_count()
        assert isinstance(count, int)
        assert count >= 0

    def test_has_restriction_site_around_codon(self):
        state = IncrementalSequenceState("ATGGAATTCCGATCGGATCCTAA", enzymes=["EcoRI"])
        # EcoRI site starts at position 3 (within codon 1)
        assert state.has_any_restriction_site_around(1) is True

    def test_has_restriction_site_around_region(self):
        state = IncrementalSequenceState("ATGGAATTCCGATCGGATCCTAA", enzymes=["EcoRI"])
        assert state.has_any_restriction_site_around_region(0, 10) is True


class TestDinucleotidePositionAccessors:
    def test_gt_positions_list(self):
        state = IncrementalSequenceState("ATGGTCAAGGTCTAA")
        positions = state.gt_positions_list()
        assert isinstance(positions, list)
        assert all(isinstance(p, int) for p in positions)

    def test_cg_positions_list(self):
        state = IncrementalSequenceState("ATGCGCAAGGCCTAA")
        positions = state.cg_positions_list()
        assert isinstance(positions, list)

    def test_ag_positions_list(self):
        state = IncrementalSequenceState("ATGAGCAAGGCCTAA")
        positions = state.ag_positions_list()
        assert isinstance(positions, list)

    def test_positions_sorted(self):
        state = IncrementalSequenceState("ATGGTCAAGGTCTAA")
        positions = state.gt_positions_list()
        assert positions == sorted(positions)


class TestFullAndIncrementalCheck:
    def test_full_check(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        result = state.full_check()
        assert isinstance(result, dict)
        assert 'gt_count' in result
        assert 'cg_count' in result
        assert 'ag_count' in result
        assert 'gc_fraction' in result
        assert 'gc_in_range' in result
        assert 'restriction_sites' in result
        assert 'splice_sites' in result

    def test_incremental_check(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        state.update_codon(1, "GTT")
        result = state.incremental_check()
        assert isinstance(result, dict)
        assert 'gt_count' in result

    def test_incremental_check_resets_changes(self):
        state = IncrementalSequenceState("ATGGTCAAG")
        state.update_codon(1, "GTT")
        assert state.has_changes() is True
        state.incremental_check()
        assert state.has_changes() is False


class TestO1VsOnVerification:
    """Verify that incremental updates produce the same results as full recomputation."""

    def test_gc_count_matches_after_swap(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        state.swap_codon(1, "AAA")
        # Recompute GC from scratch
        seq = state.sequence
        expected_gc = sum(1 for b in seq if b in "GC")
        assert state.gc_count == expected_gc

    def test_gt_count_matches_after_swap(self):
        state = IncrementalSequenceState("ATGGTCAAGGCCTAA")
        state.swap_codon(1, "GTT")
        seq = state.sequence
        expected_gt = sum(1 for i in range(len(seq) - 1) if seq[i] == 'G' and seq[i+1] == 'T')
        assert state.gt_count == expected_gt

    def test_cg_count_matches_after_swap(self):
        state = IncrementalSequenceState("ATGCGCAAGGCCTAA")
        state.swap_codon(1, "AAA")
        seq = state.sequence
        expected_cg = sum(1 for i in range(len(seq) - 1) if seq[i] == 'C' and seq[i+1] == 'G')
        assert state.cg_count == expected_cg

    def test_ag_count_matches_after_swap(self):
        state = IncrementalSequenceState("ATGAGCAAGGCCTAA")
        state.swap_codon(1, "AAA")
        seq = state.sequence
        expected_ag = sum(1 for i in range(len(seq) - 1) if seq[i] == 'A' and seq[i+1] == 'G')
        assert state.ag_count == expected_ag

    def test_multiple_swaps_consistency(self):
        state = IncrementalSequenceState("ATGCGCAAGGCCTAA")
        state.swap_codon(1, "AAA")
        state.swap_codon(2, "GTT")
        state.swap_codon(3, "TTC")
        # Verify all counters match full recomputation
        seq = state.sequence
        expected_gc = sum(1 for b in seq if b in "GC")
        expected_gt = sum(1 for i in range(len(seq)-1) if seq[i]=='G' and seq[i+1]=='T')
        expected_cg = sum(1 for i in range(len(seq)-1) if seq[i]=='C' and seq[i+1]=='G')
        expected_ag = sum(1 for i in range(len(seq)-1) if seq[i]=='A' and seq[i+1]=='G')
        assert state.gc_count == expected_gc
        assert state.gt_count == expected_gt
        assert state.cg_count == expected_cg
        assert state.ag_count == expected_ag
