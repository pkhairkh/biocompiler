"""Test BioCompiler Type System v7.0.0 — CODON_TABLE, AA_TO_CODONS, BLOSUM62, and 8 predicates."""

import pytest
from biocompiler.type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_gt_dinucleotide, check_no_avoidable_gt,
    check_valid_coding_seq, check_conservation_score, check_codon_optimality,
    SpliceVerdict, PREDICATE_NAMES,
)


# ────────────────────────────────────────────────────────────
# CODON_TABLE tests
# ────────────────────────────────────────────────────────────

class TestCodonTable:
    """Tests for the standard genetic code CODON_TABLE."""

    def test_codon_table_has_64_entries(self):
        """Verify exactly 64 codons in the table."""
        assert len(CODON_TABLE) == 64

    def test_codon_table_no_invalid_entries(self):
        """All keys are 3-character DNA strings (A/T/G/C only)."""
        valid_bases = {"A", "T", "G", "C"}
        for codon in CODON_TABLE:
            assert len(codon) == 3, f"Codon {codon!r} is not 3 chars"
            assert all(b in valid_bases for b in codon), f"Codon {codon!r} has invalid bases"

    def test_codon_table_stop_codons(self):
        """TAA, TAG, TGA all map to '*' (stop)."""
        assert CODON_TABLE["TAA"] == "*"
        assert CODON_TABLE["TAG"] == "*"
        assert CODON_TABLE["TGA"] == "*"


# ────────────────────────────────────────────────────────────
# AA_TO_CODONS reverse mapping
# ────────────────────────────────────────────────────────────

class TestAAtoCodons:
    """Tests for the reverse amino acid to codon mapping."""

    def test_aa_to_codons_reverse(self):
        """Every codon in AA_TO_CODONS maps back to the correct amino acid in CODON_TABLE."""
        for aa, codons in AA_TO_CODONS.items():
            for codon in codons:
                assert CODON_TABLE[codon] == aa, (
                    f"Codon {codon} maps to {CODON_TABLE[codon]!r}, not {aa!r}"
                )


# ────────────────────────────────────────────────────────────
# BLOSUM62 tests
# ────────────────────────────────────────────────────────────

class TestBlosum62:
    """Tests for the BLOSUM62 substitution matrix."""

    def test_blosum62_diagonal_positive(self):
        """All (X, X) diagonal entries are positive (conservative self-substitution)."""
        for (a1, a2), score in BLOSUM62.items():
            if a1 == a2:
                assert score > 0, f"BLOSUM62({a1},{a1})={score} is not positive"


# ────────────────────────────────────────────────────────────
# Predicate check tests
# ────────────────────────────────────────────────────────────

class TestCheckNoStopCodons:
    """Predicate 1: No internal stop codons."""

    def test_check_no_stop_codons_pass(self):
        """Clean sequence with no internal stops passes."""
        # ATG (M) GCT (A) GCT (A) TAA (stop, last codon is allowed)
        seq = "ATGGCTGCTTAA"
        result = check_no_stop_codons(seq)
        assert result.passed is True

    def test_check_no_stop_codons_fail(self):
        """Sequence with internal TAA fails."""
        # ATG (M) TAA (stop!) GCT (A) TAA (stop)
        seq = "ATGTAAGCTTAA"
        result = check_no_stop_codons(seq)
        assert result.passed is False
        assert len(result.positions) > 0


class TestCheckNoCrypticSplice:
    """Predicate 2: No cryptic splice sites."""

    def test_check_no_cryptic_splice_pass(self):
        """Sequence without GT dinucleotides passes."""
        # A sequence with no GT at all
        seq = "ATGGCTAAGCCTAA"
        result = check_no_cryptic_splice(seq)
        assert result.passed is True


class TestCheckNoCpGIsland:
    """Predicate 3: No CpG islands."""

    def test_check_no_cpg_island_pass(self):
        """Low CG content sequence passes."""
        # Mostly A/T with very few CG dinucleotides
        seq = "ATGAAAATTTTTAAAATTTTTAAAATTTTTAAA" * 10
        result = check_no_cpg_island(seq)
        assert result.passed is True


class TestCheckValidCodingSeq:
    """Predicate 6: Valid coding sequence."""

    def test_check_valid_coding_seq_pass(self):
        """Valid sequence (length divisible by 3, all valid codons) passes."""
        seq = "ATGGCTGCCTAA"
        result = check_valid_coding_seq(seq)
        assert result.passed is True

    def test_check_valid_coding_seq_fail_bad_length(self):
        """Sequence whose length is not divisible by 3 fails."""
        seq = "ATGGCTGC"  # 8 bases, not divisible by 3
        result = check_valid_coding_seq(seq)
        assert result.passed is False
        assert "not divisible by 3" in result.details


class TestCheckNoAvoidableGT:
    """Predicate 5 (relaxed): No avoidable GT dinucleotides."""

    def test_check_no_avoidable_gt_valine_pass(self):
        """Valine codons (unavoidable GT) should pass the relaxed check.

        All Valine codons start with GT (GTT, GTC, GTA, GTG), so the
        GT is unavoidable and the relaxed predicate should pass.
        """
        # A simple sequence containing only Valine and a stop codon
        # ATG (M) GTT (V) GTT (V) TAA (stop)
        seq = "ATGGTTGTTTAA"
        result = check_no_avoidable_gt(seq)
        assert result.passed is True


class TestPredicateNames:
    """Test that all predicate names are defined."""

    def test_28_predicates(self):
        assert len(PREDICATE_NAMES) == 28
        assert "NoStopCodons" in PREDICATE_NAMES
        assert "NoCrypticSplice" in PREDICATE_NAMES
        assert "NoCpGIsland" in PREDICATE_NAMES
        assert "NoRestrictionSite" in PREDICATE_NAMES
        assert "NoGTDinucleotide" in PREDICATE_NAMES
        assert "ValidCodingSeq" in PREDICATE_NAMES
        assert "ConservationScore" in PREDICATE_NAMES
        assert "CodonOptimality" in PREDICATE_NAMES
