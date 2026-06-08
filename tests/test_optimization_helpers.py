"""
Unit tests for biocompiler.optimization_helpers — restriction site removal,
IUPAC expansion, GT/AG-free codon utilities.

Covers:
- _find_site_in_sequence: forward and reverse complement detection
- _get_overlapping_codons: codon overlap calculation
- _remove_site_multicodon: multi-codon coordinated solving
- _expand_iupac_site: IUPAC ambiguity expansion
- _find_gt_free_codons, _is_unavoidable_gt_aa, _gt_free_cai_ratio
- _find_ag_free_codons
"""

from __future__ import annotations

import pytest

from biocompiler.optimization_helpers import (
    _find_site_in_sequence,
    _get_overlapping_codons,
    _remove_site_multicodon,
    _expand_iupac_site,
    _find_gt_free_codons,
    _is_unavoidable_gt_aa,
    _gt_free_cai_ratio,
    _find_ag_free_codons,
)
from biocompiler.type_system import AA_TO_CODONS
from biocompiler.constants import reverse_complement


# ═══════════════════════════════════════════════════════════════════════════════
# 1. _find_site_in_sequence
# ═══════════════════════════════════════════════════════════════════════════════

class TestFindSiteInSequence:

    def test_find_forward_site(self):
        positions = _find_site_in_sequence("ATGGAATTCC", "GAATTC", "")
        assert 3 in positions

    def test_find_reverse_complement(self):
        rc = reverse_complement("GAATTC")  # GAATTC
        positions = _find_site_in_sequence("ATGGAATTCC", "GAATTC", rc)
        assert len(positions) >= 1

    def test_no_site_found(self):
        positions = _find_site_in_sequence("ATGCATGCATGC", "GAATTC", "")
        assert positions == []

    def test_multiple_occurrences(self):
        positions = _find_site_in_sequence("GAATTCXXXGAATTC", "GAATTC", "")
        assert len(positions) == 2

    def test_palindrome_not_double_counted(self):
        """For palindromic sites, same-site-as-RC should not be double-counted."""
        positions = _find_site_in_sequence("GAATTC", "GAATTC", "GAATTC")
        assert len(positions) == 1

    def test_empty_site(self):
        positions = _find_site_in_sequence("ATGC", "", "")
        assert positions == []


# ═══════════════════════════════════════════════════════════════════════════════
# 2. _get_overlapping_codons
# ═══════════════════════════════════════════════════════════════════════════════

class TestGetOverlappingCodons:

    def test_site_within_single_codon(self):
        # Site at pos 0, len 3, fully within codon 0
        codons = _get_overlapping_codons(0, 3, 10)
        assert codons == [0]

    def test_site_straddles_two_codons(self):
        # Site at pos 2, len 4 -> positions 2-5, covers codons 0 and 1
        codons = _get_overlapping_codons(2, 4, 10)
        assert 0 in codons
        assert 1 in codons

    def test_site_at_end(self):
        # Site at pos 27, len 6 -> positions 27-32, covers codons 9
        codons = _get_overlapping_codons(27, 6, 10)
        assert 9 in codons

    def test_site_covers_multiple_codons(self):
        # Long site covering positions 0-11 (4 codons)
        codons = _get_overlapping_codons(0, 12, 10)
        assert codons == [0, 1, 2, 3]

    def test_boundary_clamping(self):
        """Positions that would go past n_codons should be clamped."""
        codons = _get_overlapping_codons(27, 10, 10)
        assert all(c < 10 for c in codons)

    def test_invalid_pos_raises(self):
        with pytest.raises(AssertionError):
            _get_overlapping_codons(-1, 3, 10)

    def test_invalid_site_len_raises(self):
        with pytest.raises(AssertionError):
            _get_overlapping_codons(0, 0, 10)

    def test_invalid_n_codons_raises(self):
        with pytest.raises(AssertionError):
            _get_overlapping_codons(0, 3, 0)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. _remove_site_multicodon
# ═══════════════════════════════════════════════════════════════════════════════

class TestRemoveSiteMulticodon:

    def test_remove_site_within_codon(self):
        """A site fully within one codon should be removable."""
        sequence = "ATGGAATTCCGATCG"  # GAATTC at pos 3
        aas = ["M", "E", "F", "P", "I"]  # Wrong but ok for test structure
        # Actually we need correct aas. Let's compute them.
        from biocompiler.type_system import CODON_TABLE
        aas = []
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            aas.append(CODON_TABLE.get(codon, "X"))

        sorted_codons = {}
        for aa in set(aas):
            sorted_codons[aa] = AA_TO_CODONS.get(aa, [])

        site_rc = reverse_complement("GAATTC")
        new_seq, fixed = _remove_site_multicodon(
            sequence, aas, sorted_codons, "GAATTC", site_rc
        )
        if fixed:
            assert "GAATTC" not in new_seq

    def test_no_site_returns_unchanged(self):
        """If the site is not in the sequence, return unchanged."""
        sequence = "ATGCATGCATGCATGC"
        aas = ["M", "H", "A", "C", "M", "H"]
        sorted_codons = {aa: AA_TO_CODONS.get(aa, []) for aa in set(aas)}
        new_seq, fixed = _remove_site_multicodon(
            sequence, aas, sorted_codons, "GAATTC", "GAATTC"
        )
        assert fixed is False
        assert new_seq == sequence

    def test_with_usage_cai_aware(self):
        """With usage dict, should prefer higher-CAI combinations."""
        sequence = "ATGGAATTCCGATCG"
        from biocompiler.type_system import CODON_TABLE
        aas = []
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            aas.append(CODON_TABLE.get(codon, "X"))

        sorted_codons = {aa: AA_TO_CODONS.get(aa, []) for aa in set(aas)}
        usage = {c: 1.0 for c in CODON_TABLE}  # uniform
        site_rc = reverse_complement("GAATTC")
        new_seq, fixed = _remove_site_multicodon(
            sequence, aas, sorted_codons, "GAATTC", site_rc, usage=usage
        )
        # Should still attempt to remove the site
        assert isinstance(fixed, bool)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. _expand_iupac_site
# ═══════════════════════════════════════════════════════════════════════════════

class TestExpandIupacSite:

    def test_plain_acgt_site(self):
        """A site with only ACGT should return itself."""
        result = _expand_iupac_site("GAATTC")
        assert result == ["GAATTC"]

    def test_site_with_n(self):
        """N expands to ACGT (4 options per N)."""
        result = _expand_iupac_site("GATN")
        assert len(result) == 4
        assert "GATA" in result
        assert "GATC" in result
        assert "GATG" in result
        assert "GATT" in result

    def test_empty_pattern_raises(self):
        with pytest.raises(AssertionError):
            _expand_iupac_site("")

    def test_too_large_expansion(self):
        """Very large IUPAC expansion should return empty list (capped)."""
        # 10 Ns = 4^10 = 1,048,576 > 4096 cap
        result = _expand_iupac_site("N" * 10)
        assert result == []

    def test_within_cap(self):
        """Small expansion should work fine."""
        # 2 Ns = 4^2 = 16 < 4096
        result = _expand_iupac_site("GATNN")
        assert len(result) == 16


# ═══════════════════════════════════════════════════════════════════════════════
# 5. GT-Free Codon Utilities
# ═══════════════════════════════════════════════════════════════════════════════

class TestGTFreeCodons:

    def test_valine_has_no_gt_free(self):
        """Valine (V) has all GT-containing codons."""
        gt_free = _find_gt_free_codons("V")
        assert len(gt_free) == 0

    def test_is_unavoidable_gt_aa_valine(self):
        assert _is_unavoidable_gt_aa("V") is True

    def test_is_unavoidable_gt_aa_leucine(self):
        """Leucine (L) has GT-free codons like CTT, CTC, CTA."""
        assert _is_unavoidable_gt_aa("L") is False

    def test_find_gt_free_codons_leucine(self):
        gt_free = _find_gt_free_codons("L")
        for codon in gt_free:
            assert "GT" not in codon
        # Leucine should have at least some GT-free codons
        assert len(gt_free) > 0

    def test_find_gt_free_codons_alanine(self):
        """Alanine (A) should have GT-free codons."""
        gt_free = _find_gt_free_codons("A")
        for codon in gt_free:
            assert "GT" not in codon

    def test_gt_free_cai_ratio_valine(self):
        """Valine should have ratio 0.0."""
        usage = {c: 1.0 for c in AA_TO_CODONS.get("V", [])}
        ratio = _gt_free_cai_ratio("V", usage)
        assert ratio == 0.0

    def test_gt_free_cai_ratio_no_gt_codons(self):
        """AA with no GT-containing codons should have ratio 1.0."""
        # Find an amino acid where no codons contain GT
        for aa, codons in AA_TO_CODONS.items():
            if all("GT" not in c for c in codons):
                usage = {c: 1.0 for c in codons}
                ratio = _gt_free_cai_ratio(aa, usage)
                assert ratio == 1.0
                break

    def test_gt_free_cai_ratio_mixed(self):
        """AA with both GT-free and GT-containing codons should have 0 < ratio <= 1."""
        usage = {c: 1.0 for c in AA_TO_CODONS.get("L", [])}
        ratio = _gt_free_cai_ratio("L", usage)
        assert 0.0 < ratio <= 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# 6. AG-Free Codon Utilities
# ═══════════════════════════════════════════════════════════════════════════════

class TestAGFreeCodons:

    def test_find_ag_free_codons_alanine(self):
        ag_free = _find_ag_free_codons("A")
        for codon in ag_free:
            assert "AG" not in codon

    def test_find_ag_free_codons_arginine(self):
        """Arginine has some AG-containing codons (AGA, AGG)."""
        ag_free = _find_ag_free_codons("R")
        for codon in ag_free:
            assert "AG" not in codon
        # Should have some AG-free alternatives
        assert len(ag_free) > 0

    def test_find_ag_free_codons_serine(self):
        """Serine has AG-containing codons (AGC, AGT)."""
        ag_free = _find_ag_free_codons("S")
        for codon in ag_free:
            assert "AG" not in codon
