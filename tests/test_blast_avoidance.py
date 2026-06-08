"""Tests for BLAST match avoidance module."""

import pytest
from biocompiler.optimizer.blast_avoidance import (
    check_kmer_overlap,
    eliminate_kmer_overlaps,
    avoid_blast_matches,
    BLASTMatch,
)


class TestKmerOverlap:
    """Test k-mer overlap detection."""

    def test_no_overlap(self):
        """Sequences with no shared k-mers should return empty list."""
        seq = "ATCGATCGATCGATCG"
        refs = ["TTTTTTTTTTTTTTTT"]
        overlaps = check_kmer_overlap(seq, refs, k=15)
        assert overlaps == []

    def test_exact_match_overlap(self):
        """A sequence that is a substring of the reference should be detected."""
        seq = "ATCGATCGATCGATCG"
        refs = ["AAAA" + seq + "TTTT"]
        overlaps = check_kmer_overlap(seq, refs, k=15)
        assert len(overlaps) > 0

    def test_short_kmers(self):
        """Short k values should find more overlaps."""
        seq = "ATCG"
        refs = ["ATCG"]
        overlaps = check_kmer_overlap(seq, refs, k=3)
        assert len(overlaps) > 0

    def test_empty_reference(self):
        """Empty reference list should return no overlaps."""
        seq = "ATCGATCG"
        overlaps = check_kmer_overlap(seq, [], k=5)
        assert overlaps == []

    def test_sequence_shorter_than_k(self):
        """Sequence shorter than k should return no overlaps."""
        seq = "ATC"
        refs = ["ATCGATCG"]
        overlaps = check_kmer_overlap(seq, refs, k=15)
        assert overlaps == []


class TestEliminateKmerOverlaps:
    """Test k-mer overlap elimination."""

    def test_preserves_protein(self):
        """Eliminated sequence should translate to the same protein."""
        seq = "ATGAAAGCGTGA"
        protein = "MKA"
        organism = "Escherichia_coli"
        # Create overlaps that span the whole sequence
        overlaps = [(0, 12, "ATGAAAGCGTGA")]
        result = eliminate_kmer_overlaps(seq, protein, organism, overlaps, max_iterations=10)
        # Result should still translate to MKA
        from biocompiler.translation import translate
        assert translate(result) == protein

    def test_no_overlaps_returns_same(self):
        """With no overlaps, the sequence should be returned unchanged."""
        seq = "ATGAAAGCGTGA"
        protein = "MKA"
        organism = "Escherichia_coli"
        result = eliminate_kmer_overlaps(seq, protein, organism, [], max_iterations=10)
        assert result == seq


class TestAvoidBlastMatches:
    """Test full BLAST avoidance pipeline."""

    def test_no_references(self):
        """With no references, sequence should be returned unchanged."""
        seq = "ATGAAAGCGTGA"
        protein = "MKA"
        organism = "Escherichia_coli"
        result = avoid_blast_matches(seq, protein, organism, reference_sequences=[], word_size=15)
        assert result == seq

    def test_preserves_protein(self):
        """Optimized sequence should translate to the same protein."""
        seq = "ATGAAAGCGTGA"
        protein = "MKA"
        organism = "Escherichia_coli"
        # Use a reference that shares k-mers
        ref = "ATGAAAGCGTGAXXXX"
        result = avoid_blast_matches(seq, protein, organism, reference_sequences=[ref], word_size=8)
        from biocompiler.translation import translate
        assert translate(result) == protein


class TestBLASTMatch:
    """Test BLASTMatch dataclass."""

    def test_creation(self):
        """Test basic dataclass creation."""
        match = BLASTMatch(
            query_start=10, query_end=25,
            reference_start=30, reference_end=45,
            identity=95.0, e_value=1e-10,
            reference_id="ref1",
        )
        assert match.query_start == 10
        assert match.identity == 95.0

    def test_default_reference_id(self):
        """Test default reference_id is empty string."""
        match = BLASTMatch(
            query_start=0, query_end=10,
            reference_start=0, reference_end=10,
            identity=100.0, e_value=0.0,
        )
        assert match.reference_id == ""


class TestTypeSystemIntegration:
    """Test that BLAST avoidance integrates with the type system."""

    def test_check_no_blast_matches_pass(self):
        """No BLAST matches should result in PASS."""
        from biocompiler.type_system.checks import check_no_blast_matches
        seq = "ATCGATCGATCGATCG"
        result = check_no_blast_matches(seq, reference_sequences=["TTTTTTTTTTTTTTTT"], k=15)
        assert result.passed is True

    def test_check_no_blast_matches_fail(self):
        """BLAST matches should result in FAIL."""
        from biocompiler.type_system.checks import check_no_blast_matches
        seq = "ATCGATCGATCGATCG"
        ref = ["XXXX" + seq + "YYYY"]
        result = check_no_blast_matches(seq, reference_sequences=ref, k=10)
        assert result.passed is False

    def test_evaluate_no_blast_matches(self):
        """Test the high-level evaluate function."""
        from biocompiler.type_system.predicates import evaluate_no_blast_matches
        seq = "ATCGATCGATCGATCG"
        result = evaluate_no_blast_matches(seq, reference_sequences=[], k=15)
        assert result.verdict.value == "PASS"
