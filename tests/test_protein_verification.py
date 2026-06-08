"""
BioCompiler Protein Verification Tests
========================================

Tests for the protein_verification module that double-checks optimized DNA
actually encodes the input protein.

Test categories:
  1. Correct translation passes verification
  2. Single mismatch detected at correct position
  3. Premature stop codon detected
  4. Missing stop codon detected (when expected)
  5. Wrong length DNA detected
  6. Full end-to-end: optimize insulin, verify translation
  7. verify_and_raise raises TranslationVerificationError on failure
  8. Edge cases: empty sequences, single amino acid, sequences with stop codon
"""

from __future__ import annotations

import pytest

from biocompiler.protein_verification import (
    PositionMismatch,
    VerificationResult,
    verify_translation,
    verify_and_raise,
)
from biocompiler.exceptions import TranslationVerificationError
from biocompiler.constants import CODON_TABLE, AA_TO_CODONS


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _protein_to_dna(protein: str, codon_choice: str = "first") -> str:
    """Convert a protein sequence to DNA using a deterministic codon choice.

    Args:
        protein: Amino acid sequence (single-letter codes).
        codon_choice: Which codon to use — "first" for the first listed
            codon for each amino acid, "last" for the last.

    Returns:
        DNA sequence (no stop codon).
    """
    codons = []
    for aa in protein:
        available = AA_TO_CODONS.get(aa, [])
        if not available:
            raise ValueError(f"No codons for amino acid '{aa}'")
        if codon_choice == "first":
            codons.append(available[0])
        elif codon_choice == "last":
            codons.append(available[-1])
        else:
            codons.append(available[0])
    return "".join(codons)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Correct translation passes verification
# ═══════════════════════════════════════════════════════════════════════════════

class TestCorrectTranslation:
    """Verify that a correctly translated DNA sequence passes all checks."""

    def test_simple_protein_passes(self):
        """A simple protein (MALL) translated correctly should pass."""
        # M=ATG, A=GCT, L=CTT, L=CTT
        dna = "ATGGCTCTTCTT"
        protein = "MALL"
        result = verify_translation(dna, protein)
        assert result.is_valid is True
        assert result.matches_expected is True
        assert result.mismatches == []
        assert result.has_premature_stop is False
        assert result.length_correct is True
        assert result.translated_protein == protein

    def test_all_twenty_amino_acids(self):
        """DNA encoding all 20 standard amino acids should pass."""
        # Use first codon for each amino acid
        protein = "ACDEFGHIKLMNPQRSTVWY"
        dna = _protein_to_dna(protein)
        result = verify_translation(dna, protein)
        assert result.is_valid is True
        assert result.matches_expected is True
        assert result.mismatches == []
        assert result.translated_protein == protein

    def test_single_methionine(self):
        """Single amino acid (M) should pass."""
        dna = "ATG"
        result = verify_translation(dna, "M")
        assert result.is_valid is True
        assert result.matches_expected is True
        assert result.translated_protein == "M"

    def test_with_terminal_stop_codon(self):
        """DNA with a terminal stop codon should pass and report has_stop_codon."""
        # M=ATG, A=GCT, stop=TAA
        dna = "ATGGCTTAA"
        protein = "MA"
        result = verify_translation(dna, protein)
        assert result.is_valid is True
        assert result.matches_expected is True
        assert result.has_stop_codon is True
        assert result.length_correct is True
        assert result.translated_protein == "MA"

    def test_without_terminal_stop_codon(self):
        """DNA without a terminal stop codon should also pass."""
        dna = "ATGGCT"
        protein = "MA"
        result = verify_translation(dna, protein)
        assert result.is_valid is True
        assert result.matches_expected is True
        assert result.has_stop_codon is False
        assert result.length_correct is True

    def test_different_codon_choices_same_protein(self):
        """Different codon choices for the same protein should all pass."""
        protein = "MALWMR"
        dna_first = _protein_to_dna(protein, "first")
        dna_last = _protein_to_dna(protein, "last")
        for dna in [dna_first, dna_last]:
            result = verify_translation(dna, protein)
            assert result.is_valid is True
            assert result.translated_protein == protein

    def test_insulin_bchain(self):
        """Insulin B-chain should pass verification."""
        insulin_b = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
        dna = _protein_to_dna(insulin_b)
        result = verify_translation(dna, insulin_b)
        assert result.is_valid is True
        assert result.translated_protein == insulin_b


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Single mismatch detected
# ═══════════════════════════════════════════════════════════════════════════════

class TestSingleMismatch:
    """Verify that a single amino acid mismatch is detected at the correct position."""

    def test_first_position_mismatch(self):
        """Mismatch at position 0 should be detected."""
        # Expected: M=ATG, A=GCT. Actual: K=AAG, A=GCT
        dna = "AAGGCT"  # AAG=K instead of ATG=M
        result = verify_translation(dna, "MA")
        assert result.is_valid is False
        assert result.matches_expected is False
        assert len(result.mismatches) == 1
        mm = result.mismatches[0]
        assert mm.position == 0
        assert mm.expected == "M"
        assert mm.actual == "K"
        assert mm.codon_used == "AAG"

    def test_middle_position_mismatch(self):
        """Mismatch in the middle of the sequence."""
        # Expected: M A L -> ATG GCT CTT
        # Actual:   M K L -> ATG AAG CTT
        dna = "ATGAAGCTT"
        result = verify_translation(dna, "MAL")
        assert result.is_valid is False
        assert len(result.mismatches) == 1
        mm = result.mismatches[0]
        assert mm.position == 1
        assert mm.expected == "A"
        assert mm.actual == "K"
        assert mm.codon_used == "AAG"

    def test_last_position_mismatch(self):
        """Mismatch at the last position."""
        # Expected: M A L -> ATG GCT CTT
        # Actual:   M A V -> ATG GCT GTT
        dna = "ATGGCTGTT"
        result = verify_translation(dna, "MAL")
        assert result.is_valid is False
        assert len(result.mismatches) == 1
        mm = result.mismatches[0]
        assert mm.position == 2
        assert mm.expected == "L"
        assert mm.actual == "V"
        assert mm.codon_used == "GTT"

    def test_mismatch_with_stop_codon_sequence(self):
        """Mismatch detection works when DNA has a terminal stop codon."""
        # Expected: MA + stop. Actual: MK + stop
        dna = "ATGAAGTAA"  # K at position 1 instead of A
        result = verify_translation(dna, "MA")
        assert result.is_valid is False
        assert len(result.mismatches) == 1
        mm = result.mismatches[0]
        assert mm.position == 1
        assert mm.expected == "A"
        assert mm.actual == "K"

    def test_position_mismatch_str_representation(self):
        """PositionMismatch __str__ should produce readable output."""
        mm = PositionMismatch(position=5, expected="A", actual="K", codon_used="AAG")
        s = str(mm)
        assert "5" in s
        assert "A" in s
        assert "K" in s
        assert "AAG" in s


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Premature stop codon detected
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrematureStopCodon:
    """Verify that premature stop codons are detected."""

    def test_stop_codon_in_middle(self):
        """A stop codon in the middle of the CDS should be detected."""
        # M * L -> ATG TAA CTT (TAA is stop at position 1)
        dna = "ATGTAACTT"
        result = verify_translation(dna, "MAL")
        assert result.is_valid is False
        assert result.has_premature_stop is True
        assert result.matches_expected is False
        # Should have mismatches starting at position 1
        assert len(result.mismatches) >= 1

    def test_stop_codon_at_first_position(self):
        """Stop codon at the very first amino acid position."""
        # TAA = stop at position 0
        dna = "TAAGCTCTT"
        result = verify_translation(dna, "MAL")
        assert result.is_valid is False
        assert result.has_premature_stop is True

    def test_stop_codon_just_before_end(self):
        """Stop codon one position before the end."""
        # M A * -> ATG GCT TGA
        dna = "ATGGCTTGA"
        result = verify_translation(dna, "MAL")
        assert result.is_valid is False
        assert result.has_premature_stop is True

    def test_tag_stop_codon_detected(self):
        """TAG stop codon should also be detected as premature."""
        # M * A -> ATG TAG GCT
        dna = "ATGTAGGCT"
        result = verify_translation(dna, "MAA")
        assert result.has_premature_stop is True

    def test_tga_stop_codon_detected(self):
        """TGA stop codon should also be detected as premature."""
        # M * A -> ATG TGA GCT
        dna = "ATGTGAGCT"
        result = verify_translation(dna, "MAA")
        assert result.has_premature_stop is True

    def test_terminal_stop_not_premature(self):
        """Terminal stop codon should NOT be flagged as premature."""
        # MA + TAA stop
        dna = "ATGGCTTAA"
        result = verify_translation(dna, "MA")
        assert result.has_premature_stop is False
        assert result.has_stop_codon is True
        assert result.is_valid is True


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Missing stop codon detected
# ═══════════════════════════════════════════════════════════════════════════════

class TestMissingStopCodon:
    """Verify that missing stop codons are reported correctly."""

    def test_no_stop_codon_reported(self):
        """Sequence without terminal stop codon should report has_stop_codon=False."""
        dna = "ATGGCT"
        result = verify_translation(dna, "MA")
        assert result.has_stop_codon is False
        # This is still valid (biocompiler doesn't include stop in results)
        assert result.is_valid is True

    def test_with_stop_codon_reported(self):
        """Sequence with terminal stop codon should report has_stop_codon=True."""
        dna = "ATGGCTTAA"
        result = verify_translation(dna, "MA")
        assert result.has_stop_codon is True
        assert result.is_valid is True

    def test_non_stop_at_expected_stop_position(self):
        """If DNA has length protein*3+3 but last codon is NOT a stop, report correctly."""
        # MA + GCT (not a stop codon at the terminal position)
        dna = "ATGGCTGCT"  # 3 codons: M, A, A (not stop)
        protein = "MA"
        result = verify_translation(dna, protein)
        # Length is 9 = 2*3 + 3, so length_correct=True
        # But the last codon is NOT a stop, so has_stop_codon=False
        # And the translated protein has an extra AA
        assert result.has_stop_codon is False
        # The translated protein should be "MAA" (3 codons, none is stop)
        # which doesn't match "MA" — so it's invalid
        assert result.matches_expected is False


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Wrong length DNA detected
# ═══════════════════════════════════════════════════════════════════════════════

class TestWrongLength:
    """Verify that incorrect DNA lengths are detected."""

    def test_dna_too_short(self):
        """DNA shorter than protein*3 should be flagged."""
        # Protein "MA" needs 6 bases, we only provide 3
        dna = "ATG"
        result = verify_translation(dna, "MA")
        assert result.length_correct is False
        assert result.is_valid is False

    def test_dna_too_long(self):
        """DNA longer than protein*3+3 should be flagged."""
        # Protein "MA" needs 6 or 9 bases, we provide 12
        dna = "ATGGCTCTTTAA"
        result = verify_translation(dna, "MA")
        assert result.length_correct is False
        assert result.is_valid is False

    def test_dna_off_by_one(self):
        """DNA that is not a multiple of 3 should be handled."""
        # Protein "MA" needs 6 bases, we provide 7
        dna = "ATGGCTC"
        result = verify_translation(dna, "MA")
        assert result.length_correct is False
        assert result.is_valid is False

    def test_correct_length_no_stop(self):
        """DNA with exactly protein*3 length should have length_correct=True."""
        dna = "ATGGCT"
        result = verify_translation(dna, "MA")
        assert result.length_correct is True

    def test_correct_length_with_stop(self):
        """DNA with exactly protein*3+3 length should have length_correct=True."""
        dna = "ATGGCTTAA"
        result = verify_translation(dna, "MA")
        assert result.length_correct is True

    def test_empty_dna_nonempty_protein(self):
        """Empty DNA with non-empty protein should fail."""
        result = verify_translation("", "MA")
        assert result.length_correct is False
        assert result.is_valid is False

    def test_empty_protein_empty_dna(self):
        """Empty DNA and empty protein should pass."""
        result = verify_translation("", "")
        assert result.is_valid is True
        assert result.matches_expected is True
        assert result.mismatches == []


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Full end-to-end: optimize insulin, verify translation
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.integration
class TestInsulinE2E:
    """End-to-end: optimize insulin for E. coli, verify translation matches."""

    INSULIN = "FVNQHLCGSHLVEALYLVCGERGFFYTPKTGIVEQCCTSICSLYQLENYCN"

    def test_ecoli_insulin_verification(self):
        """Optimize insulin for E. coli, then verify the translation."""
        from biocompiler.optimization import optimize_sequence
        result = optimize_sequence(
            self.INSULIN,
            organism="e_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=["EcoRI", "BamHI", "HindIII"],
        )
        # The optimizer now automatically verifies, but we also do it explicitly
        vr = verify_translation(result.sequence, self.INSULIN)
        assert vr.is_valid is True
        assert vr.matches_expected is True
        assert vr.mismatches == []
        assert vr.has_premature_stop is False
        assert vr.translated_protein == self.INSULIN

    def test_human_insulin_verification(self):
        """Optimize insulin for human, then verify the translation."""
        from biocompiler.optimization import optimize_sequence
        from biocompiler.exceptions import OptimizationConstraintError
        try:
            result = optimize_sequence(
                self.INSULIN,
                organism="human",
                gc_lo=0.30,
                gc_hi=0.70,
                enzymes=["EcoRI", "BamHI", "HindIII"],
            )
        except OptimizationConstraintError:
            # Some organisms/constraint combinations may be infeasible;
            # skip the test rather than fail — the verification module
            # is not responsible for optimizer feasibility.
            pytest.skip("Human insulin optimization hit GT dinucleotide constraint")
        vr = verify_translation(result.sequence, self.INSULIN)
        assert vr.is_valid is True
        assert vr.matches_expected is True
        assert vr.translated_protein == self.INSULIN

    def test_yeast_insulin_verification(self):
        """Optimize insulin for yeast, then verify the translation."""
        from biocompiler.optimization import optimize_sequence
        from biocompiler.exceptions import OptimizationConstraintError
        try:
            result = optimize_sequence(
                self.INSULIN,
                organism="yeast",
                gc_lo=0.30,
                gc_hi=0.70,
                enzymes=["EcoRI", "BamHI", "HindIII"],
            )
        except OptimizationConstraintError:
            pytest.skip("Yeast insulin optimization hit constraint")
        vr = verify_translation(result.sequence, self.INSULIN)
        assert vr.is_valid is True
        assert vr.matches_expected is True
        assert vr.translated_protein == self.INSULIN

    def test_verify_and_raise_on_valid_result(self):
        """verify_and_raise should return the VerificationResult for valid DNA."""
        from biocompiler.optimization import optimize_sequence
        result = optimize_sequence(
            self.INSULIN,
            organism="e_coli",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        vr = verify_and_raise(result.sequence, self.INSULIN)
        assert isinstance(vr, VerificationResult)
        assert vr.is_valid is True


# ═══════════════════════════════════════════════════════════════════════════════
# 7. verify_and_raise raises TranslationVerificationError on failure
# ═══════════════════════════════════════════════════════════════════════════════

class TestVerifyAndRaise:
    """Test that verify_and_raise raises TranslationVerificationError on failure."""

    def test_raises_on_mismatch(self):
        """Should raise TranslationVerificationError when DNA doesn't match protein."""
        dna = "AAGGCT"  # K, A instead of M, A
        with pytest.raises(TranslationVerificationError) as exc_info:
            verify_and_raise(dna, "MA")
        err = exc_info.value
        assert len(err.mismatches) > 0
        assert err.has_premature_stop is False

    def test_raises_on_premature_stop(self):
        """Should raise on premature stop codon."""
        dna = "ATGTAAGCT"  # M, *, A
        with pytest.raises(TranslationVerificationError) as exc_info:
            verify_and_raise(dna, "MAA")
        err = exc_info.value
        assert err.has_premature_stop is True

    def test_raises_on_wrong_length(self):
        """Should raise on wrong length DNA."""
        dna = "ATG"  # Only 1 codon for 2-AA protein
        with pytest.raises(TranslationVerificationError) as exc_info:
            verify_and_raise(dna, "MA")
        err = exc_info.value
        assert err.length_correct is False

    def test_no_raise_on_valid(self):
        """Should NOT raise for valid DNA."""
        dna = "ATGGCT"
        result = verify_and_raise(dna, "MA")
        assert result.is_valid is True

    def test_error_contains_position_details(self):
        """TranslationVerificationError should contain position-level details."""
        dna = "AAGGCT"  # K at pos 0 instead of M
        with pytest.raises(TranslationVerificationError) as exc_info:
            verify_and_raise(dna, "MA")
        err = exc_info.value
        assert len(err.mismatches) >= 1
        mm = err.mismatches[0]
        assert mm.position == 0
        assert mm.expected == "M"
        assert mm.actual == "K"
        assert mm.codon_used == "AAG"

    def test_error_message_is_informative(self):
        """The error message should contain useful diagnostic info."""
        dna = "AAGGCT"
        with pytest.raises(TranslationVerificationError) as exc_info:
            verify_and_raise(dna, "MA")
        msg = str(exc_info.value)
        assert "mismatch" in msg.lower()
        assert "Translated" in msg
        assert "Expected" in msg


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Edge cases for the protein verification module."""

    def test_lowercase_dna(self):
        """Lowercase DNA should be handled (normalized to uppercase)."""
        dna = "atggct"
        result = verify_translation(dna, "MA")
        assert result.is_valid is True

    def test_whitespace_in_dna(self):
        """Whitespace in DNA should be stripped."""
        dna = " ATGGCT "
        result = verify_translation(dna, "MA")
        assert result.is_valid is True

    def test_lowercase_protein(self):
        """Lowercase protein should be handled (normalized to uppercase)."""
        dna = "ATGGCT"
        result = verify_translation(dna, "ma")
        assert result.is_valid is True

    def test_multiple_mismatches(self):
        """Multiple mismatches should all be reported."""
        # M->K (AAG), A->R (CGT), L->I (ATT)
        dna = "AAGCGTATT"
        result = verify_translation(dna, "MAL")
        assert result.is_valid is False
        assert len(result.mismatches) == 3
        positions = [mm.position for mm in result.mismatches]
        assert positions == [0, 1, 2]

    def test_dna_longer_than_expected_protein(self):
        """Extra codons beyond expected protein should be flagged as mismatches."""
        # DNA encodes MALL but expected is MA
        dna = "ATGGCTCTTCTT"
        result = verify_translation(dna, "MA")
        assert result.matches_expected is False
        assert len(result.mismatches) > 0  # Extra positions flagged

    def test_dna_shorter_than_expected_protein(self):
        """Missing codons beyond available DNA should be flagged as mismatches."""
        # DNA encodes M but expected is MA
        dna = "ATG"
        result = verify_translation(dna, "MA")
        assert result.matches_expected is False
        assert len(result.mismatches) > 0
        # Position 1 should be marked as <missing>
        mm = result.mismatches[0]
        assert mm.position == 1
        assert mm.expected == "A"
        assert mm.actual == "<missing>"

    def test_organism_parameter_accepted(self):
        """The organism parameter should be accepted without error."""
        dna = "ATGGCT"
        result = verify_translation(dna, "MA", organism="Escherichia_coli")
        assert result.is_valid is True

    def test_all_stop_codons_as_terminal(self):
        """All three stop codons should be recognized as terminal stop codons."""
        protein = "MA"
        base_dna = "ATGGCT"
        for stop_codon in ["TAA", "TAG", "TGA"]:
            dna = base_dna + stop_codon
            result = verify_translation(dna, protein)
            assert result.has_stop_codon is True
            assert result.is_valid is True
