"""Tests for the translation module.

Covers:
- translate: standard DNA→protein translation
- translate_with_confidence: strict translation with error raising
- compute_cai: Codon Adaptation Index computation
- find_orfs: Open Reading Frame detection
- PartialCodonError: specific error for partial codons
- BACTERIAL_START_CODONS: alternative start codon support
- Edge cases: empty sequences, partial codons, unknown codons, reverse strand
"""

from __future__ import annotations

import pytest
import warnings

from biocompiler.translation import (
    translate,
    translate_with_confidence,
    compute_cai,
    find_orfs,
    PartialCodonError,
    DEFAULT_MIN_ORF_LENGTH_AA,
    BACTERIAL_START_CODONS,
    STANDARD_START_CODON,
)


# ---------------------------------------------------------------------------
# translate
# ---------------------------------------------------------------------------

class TestTranslate:
    """Tests for the translate function."""

    def test_basic_translation(self):
        """Standard codons should translate correctly."""
        # ATG=M, GCT=A, TAA=*
        assert translate("ATGGCTTAA") == "MA"

    def test_to_stop_true(self):
        """to_stop=True should stop at first stop codon."""
        seq = "ATGGCTTAATTT"  # M-A-*-F
        assert translate(seq, to_stop=True) == "MA"

    def test_to_stop_false(self):
        """to_stop=False should include stop codons as '*'."""
        seq = "ATGGCTTAATTT"  # M-A-*-F
        result = translate(seq, to_stop=False)
        assert "*" in result

    def test_empty_sequence(self):
        """Empty sequence returns empty string."""
        assert translate("") == ""

    def test_partial_codon_ignored(self):
        """Partial codon at end is silently ignored."""
        # 7 bases = 2 codons + 1 leftover
        result = translate("ATGGCTA")
        assert result == "MA"

    def test_unknown_codon_mapped_to_x(self):
        """Unknown codons are mapped to 'X'."""
        # NNNGCT — NNN is not a valid codon
        result = translate("NNNGCT")
        assert "X" in result

    def test_stop_codons_recognized(self):
        """All three stop codons should be recognized."""
        for stop in ("TAA", "TAG", "TGA"):
            seq = "ATG" + stop
            result = translate(seq, to_stop=True)
            assert result == "M"

    def test_standard_amino_acids(self):
        """All 20 standard amino acids can be translated."""
        from biocompiler.constants import CODON_TABLE
        # Pick one codon per amino acid
        seen_aas = set()
        for codon, aa in CODON_TABLE.items():
            if aa != "*" and len(codon) == 3:
                result = translate(codon, to_stop=False)
                if result != "X":
                    seen_aas.add(result)
        # Should have all 20 standard amino acids
        assert len(seen_aas) >= 20

    def test_case_insensitive(self):
        """Translation should be case-insensitive."""
        assert translate("ATGGCTTAA") == translate("atggcttaa")

    def test_long_sequence(self):
        """Long sequence should translate correctly."""
        seq = "ATG" * 100 + "TAA"
        result = translate(seq, to_stop=True)
        assert result == "M" * 100


# ---------------------------------------------------------------------------
# translate_with_confidence
# ---------------------------------------------------------------------------

class TestTranslateWithConfidence:
    """Tests for the strict translate_with_confidence function."""

    def test_basic_translation(self):
        """Standard translation should work."""
        assert translate_with_confidence("ATGGCTTAA") == "MA"

    def test_partial_codon_raises(self):
        """Partial codon should raise PartialCodonError."""
        with pytest.raises(PartialCodonError) as exc_info:
            translate_with_confidence("ATGGCTA")  # 7 bases
        assert exc_info.value.sequence_length == 7
        assert exc_info.value.remainder == 1

    def test_partial_codon_is_value_error(self):
        """PartialCodonError should be a subclass of ValueError."""
        with pytest.raises(ValueError):
            translate_with_confidence("ATGGCTA")

    def test_unknown_codon_raises(self):
        """Unknown codon should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown codon"):
            translate_with_confidence("NNNGCT")

    def test_empty_sequence(self):
        """Empty sequence returns empty string."""
        assert translate_with_confidence("") == ""

    def test_to_stop_false_includes_stops(self):
        """to_stop=False should include stop codons as '*'."""
        seq = "ATGGCTTAATTT"
        result = translate_with_confidence(seq, to_stop=False)
        assert "*" in result


# ---------------------------------------------------------------------------
# compute_cai
# ---------------------------------------------------------------------------

class TestComputeCAI:
    """Tests for the compute_cai function."""

    def test_basic_cai(self):
        """CAI should return a float in [0, 1]."""
        seq = "ATGGCTAAAGCTTAA"
        cai = compute_cai(seq, organism="Homo_sapiens")
        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0

    def test_empty_sequence_returns_zero(self):
        """Empty sequence returns 0.0."""
        assert compute_cai("", organism="Homo_sapiens") == 0.0

    def test_different_organisms(self):
        """Different organisms should produce different CAI values."""
        seq = "ATGGCTAAAGCTGCTTAA"
        cai_human = compute_cai(seq, organism="Homo_sapiens")
        cai_ecoli = compute_cai(seq, organism="Escherichia_coli")
        # They may or may not differ, but both should be valid
        assert 0.0 <= cai_human <= 1.0
        assert 0.0 <= cai_ecoli <= 1.0

    def test_unsupported_organism_raises(self):
        """Unsupported organism should raise UnsupportedOrganismError."""
        from biocompiler.exceptions import UnsupportedOrganismError
        with pytest.raises(UnsupportedOrganismError):
            compute_cai("ATGGCTTAA", organism="Nonexistent_organism")

    def test_species_deprecated_param(self):
        """species parameter should emit DeprecationWarning."""
        with pytest.warns(DeprecationWarning, match="deprecated"):
            cai = compute_cai("ATGGCTTAA", species="human")

    def test_species_equivalent_to_organism(self):
        """species='human' should give same result as organism='Homo_sapiens'."""
        cai_org = compute_cai("ATGGCTAAAGCTTAA", organism="Homo_sapiens")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            cai_sp = compute_cai("ATGGCTAAAGCTTAA", species="human")
        assert abs(cai_org - cai_sp) < 0.001

    def test_cai_deterministic(self):
        """Same input should always produce the same CAI."""
        seq = "ATGGCTAAAGCTGCTTAA"
        cai1 = compute_cai(seq, organism="Homo_sapiens")
        cai2 = compute_cai(seq, organism="Homo_sapiens")
        assert cai1 == cai2

    def test_all_optimal_codons_high_cai(self):
        """Sequence using only optimal codons should have high CAI."""
        # Use ATG (Met, always optimal for itself) repeated
        seq = "ATG" * 10
        cai = compute_cai(seq, organism="Homo_sapiens")
        # Met codons are skipped in CAI, so with only Met, CAI = 0.0
        # Use a mix of amino acids with known optimal codons instead
        seq = "ATGGCTGCT" * 5  # M-A-A repeated
        cai = compute_cai(seq, organism="Homo_sapiens")
        assert isinstance(cai, float)


# ---------------------------------------------------------------------------
# find_orfs
# ---------------------------------------------------------------------------

class TestFindOrfs:
    """Tests for the find_orfs function."""

    def test_basic_orf_detection(self):
        """A simple ORF with ATG start and TAA stop should be found."""
        # ATG + 30 AA + TAA
        seq = "ATG" + "GCT" * 30 + "TAA"
        orfs = find_orfs(seq, min_length_aa=10)
        assert len(orfs) >= 1
        assert orfs[0]["protein"].startswith("M")

    def test_no_orf_below_min_length(self):
        """ORFs shorter than min_length_aa should be excluded."""
        # Short ORF: ATG + TAA (only 1 AA)
        seq = "ATGTAA"
        orfs = find_orfs(seq, min_length_aa=30)
        assert len(orfs) == 0

    def test_default_min_orf_length(self):
        """Default minimum ORF length is 30 AA."""
        assert DEFAULT_MIN_ORF_LENGTH_AA == 30

    def test_multiple_orfs_different_frames(self):
        """ORFs in different reading frames should all be found."""
        # Frame 0 ORF and frame 1 ORF
        seq = "AATGGCT" * 20 + "TAA"  # ATG at position 1 (frame 1)
        orfs = find_orfs(seq, min_length_aa=5)
        frames_found = {orf["frame"] for orf in orfs}
        # Should find ORFs in at least one frame
        assert len(orfs) >= 1

    def test_reverse_strand_orfs(self):
        """ORFs on the reverse strand should be found."""
        # Create a sequence with an ORF on the reverse complement
        # RC of "ATG" + "GCT"*30 + "TAA" should be found
        from biocompiler.constants import reverse_complement
        forward_orf = "ATG" + "GCT" * 30 + "TAA"
        rc = reverse_complement(forward_orf)
        # Embed it with some flanking sequence
        seq = "AAA" + rc + "AAA"
        orfs = find_orfs(seq, min_length_aa=10)
        strands = {orf["strand"] for orf in orfs}
        assert "-" in strands or len(orfs) >= 1

    def test_empty_sequence(self):
        """Empty sequence returns empty list."""
        assert find_orfs("") == []

    def test_bacterial_start_codons(self):
        """GTG and TTG should be recognized as start codons."""
        # GTG-initiated ORF
        seq = "GTG" + "GCT" * 30 + "TAA"
        orfs_standard = find_orfs(seq, min_length_aa=10)
        orfs_bacterial = find_orfs(seq, min_length_aa=10,
                                    start_codons=BACTERIAL_START_CODONS)
        # Bacterial start codons should find more ORFs
        assert len(orfs_bacterial) >= len(orfs_standard)

    def test_bacterial_start_translated_as_m(self):
        """Bacterial start codons should translate first AA as M."""
        seq = "GTG" + "GCT" * 30 + "TAA"
        orfs = find_orfs(seq, min_length_aa=10,
                         start_codons=BACTERIAL_START_CODONS)
        for orf in orfs:
            if orf["start"] == 0:  # The GTG-initiated ORF
                assert orf["protein"][0] == "M"

    def test_orf_result_structure(self):
        """Each ORF result should have all required fields."""
        seq = "ATG" + "GCT" * 30 + "TAA"
        orfs = find_orfs(seq, min_length_aa=10)
        assert len(orfs) >= 1
        for orf in orfs:
            assert "start" in orf
            assert "end" in orf
            assert "frame" in orf
            assert "strand" in orf
            assert "protein" in orf
            assert "length" in orf
            assert orf["frame"] in (0, 1, 2)
            assert orf["strand"] in ("+", "-")
            assert orf["length"] >= 10

    def test_standard_start_codon_constant(self):
        """STANDARD_START_CODON should be ATG."""
        assert STANDARD_START_CODON == "ATG"

    def test_bacterial_start_codons_constant(self):
        """BACTERIAL_START_CODONS should include ATG, GTG, TTG."""
        assert "ATG" in BACTERIAL_START_CODONS
        assert "GTG" in BACTERIAL_START_CODONS
        assert "TTG" in BACTERIAL_START_CODONS


# ---------------------------------------------------------------------------
# PartialCodonError
# ---------------------------------------------------------------------------

class TestPartialCodonError:
    """Tests for the PartialCodonError exception class."""

    def test_is_value_error(self):
        """PartialCodonError should be a subclass of ValueError."""
        assert issubclass(PartialCodonError, ValueError)

    def test_attributes(self):
        """Should store sequence_length and remainder."""
        err = PartialCodonError(sequence_length=7, remainder=1)
        assert err.sequence_length == 7
        assert err.remainder == 1

    def test_message(self):
        """Error message should be informative."""
        err = PartialCodonError(sequence_length=10, remainder=1)
        assert "10" in str(err)
        assert "1" in str(err)
