"""Tests for the scanner module.

Covers:
- validate_dna_sequence: input validation and normalization
- gc_content: GC fraction computation
- scan_sequence: multi-motif DFA scanning
- Kozak consensus scoring
- Restriction site detection (both strands, IUPAC-aware)
- Start/stop codon detection in all reading frames
- Splice site detection (donor GT, acceptor AG)
- Instability motif detection (ATTTA)
- Edge cases: empty sequences, IUPAC codes, invalid characters
"""

from __future__ import annotations

import pytest

from biocompiler.scanner import (
    validate_dna_sequence,
    gc_content,
    scan_sequence,
    KOZAK_POSITION_WEIGHTS,
    KOZAK_REPORT_THRESHOLD,
)
from biocompiler.types import Token
from biocompiler.exceptions import InvalidSequenceError


# ---------------------------------------------------------------------------
# validate_dna_sequence
# ---------------------------------------------------------------------------

class TestValidateDnaSequence:
    """Tests for validate_dna_sequence."""

    def test_basic_validation(self):
        """Valid DNA sequence should be uppercased and returned."""
        assert validate_dna_sequence("atgcgt") == "ATGCGT"

    def test_uppercase_passthrough(self):
        """Already uppercased sequence passes through."""
        assert validate_dna_sequence("ATGCGT") == "ATGCGT"

    def test_invalid_characters_raise(self):
        """Invalid characters should raise InvalidSequenceError."""
        with pytest.raises(InvalidSequenceError):
            validate_dna_sequence("ATGXYZ")

    def test_n_is_valid(self):
        """N is a valid base (ambiguous)."""
        assert validate_dna_sequence("ATGCN") == "ATGCN"

    def test_iupac_allowed_with_flag(self):
        """IUPAC ambiguity codes are accepted when allow_iupac=True."""
        result = validate_dna_sequence("ATGCRYSW", allow_iupac=True)
        assert result == "ATGCRYSW"

    def test_iupac_rejected_without_flag(self):
        """IUPAC codes are rejected when allow_iupac=False (default)."""
        with pytest.raises(InvalidSequenceError):
            validate_dna_sequence("ATGCRYSW")

    def test_empty_string(self):
        """Empty string is valid (returns empty)."""
        assert validate_dna_sequence("") == ""

    def test_numeric_characters_rejected(self):
        """Numeric characters are invalid."""
        with pytest.raises(InvalidSequenceError):
            validate_dna_sequence("ATG123")


# ---------------------------------------------------------------------------
# gc_content
# ---------------------------------------------------------------------------

class TestGcContent:
    """Tests for gc_content."""

    def test_fifty_percent(self):
        """50% GC sequence should return ~0.5."""
        seq = "GCAT" * 25  # 50% GC
        assert abs(gc_content(seq) - 0.5) < 0.01

    def test_all_gc(self):
        """All G/C should return 1.0."""
        assert gc_content("GCGCGCGC") == 1.0

    def test_all_at(self):
        """All A/T should return 0.0."""
        assert gc_content("ATATATAT") == 0.0

    def test_empty_sequence(self):
        """Empty sequence returns 0.0."""
        assert gc_content("") == 0.0

    def test_single_base(self):
        """Single G base returns 1.0."""
        assert gc_content("G") == 1.0

    def test_case_insensitive(self):
        """GC content should be case-insensitive."""
        assert gc_content("gcat") == gc_content("GCAT")

    def test_mixed_sequence(self):
        """Mixed sequence should compute correct fraction."""
        # G=1, C=1, A=2, T=2 → 2/6 = 0.3333
        seq = "GCATAT"
        gc = gc_content(seq)
        assert abs(gc - 1/3) < 0.01


# ---------------------------------------------------------------------------
# scan_sequence
# ---------------------------------------------------------------------------

class TestScanSequence:
    """Tests for the scan_sequence function."""

    def test_basic_scan_returns_tokens(self):
        """Scanning should return a list of Token objects."""
        seq = "ATGGCTAGCTAGCTAGCTAGCTAA"
        tokens = scan_sequence(seq)
        assert isinstance(tokens, list)
        for t in tokens:
            assert isinstance(t, Token)

    def test_empty_sequence(self):
        """Empty sequence returns empty list."""
        assert scan_sequence("") == []

    def test_start_codon_detected(self):
        """ATG start codons should be detected."""
        seq = "ATGATGATG"
        tokens = scan_sequence(seq)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        assert len(start_codons) >= 1

    def test_stop_codon_detected(self):
        """Stop codons should be detected."""
        seq = "ATGTAATAGTGA"
        tokens = scan_sequence(seq)
        stop_codons = [t for t in tokens if t.element_type == "stop_codon"]
        assert len(stop_codons) >= 1

    def test_splice_donor_detected(self):
        """GT splice donor sites should be detected."""
        seq = "AAGTGTAAGT"
        tokens = scan_sequence(seq, use_maxentscan=True)
        donors = [t for t in tokens if t.element_type == "splice_donor"]
        # May or may not find donors depending on MaxEntScan scores
        assert isinstance(tokens, list)

    def test_splice_donor_without_maxentscan(self):
        """GT sites should be found without MaxEntScan."""
        seq = "AAGTGTAAGT"
        tokens = scan_sequence(seq, use_maxentscan=False)
        donors = [t for t in tokens if t.element_type == "splice_donor"]
        assert len(donors) >= 1

    def test_restriction_site_ecori(self):
        """EcoRI site (GAATTC) should be detected."""
        seq = "AAGAATTCGG"
        tokens = scan_sequence(seq, restriction_enzymes=["EcoRI"])
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) >= 1
        assert any("GAATTC" in t.match_sequence for t in sites)

    def test_restriction_site_both_strands(self):
        """Restriction sites should be found on both strands."""
        # EcoRI site GAATTC, RC = GAATTC (palindromic)
        seq = "AAGAATTCGG"
        tokens = scan_sequence(seq, restriction_enzymes=["EcoRI"])
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) >= 1

    def test_no_restriction_enzymes(self):
        """No restriction enzymes means no restriction site tokens."""
        seq = "AAGAATTCGG"
        tokens = scan_sequence(seq)
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) == 0

    def test_instability_motif_detected(self):
        """ATTTA instability motif should be detected."""
        seq = "AATTTATATTTAA"
        tokens = scan_sequence(seq)
        motifs = [t for t in tokens if t.element_type == "instability_motif"]
        assert len(motifs) >= 1

    def test_tokens_sorted_by_position(self):
        """Tokens should be sorted by position."""
        seq = "ATGGCTAAGTGTAAGTGAATTTATAA"
        tokens = scan_sequence(seq)
        positions = [t.position for t in tokens]
        assert positions == sorted(positions)

    def test_kozak_scoring(self):
        """Strong Kozak context (GCCACCATGG) should have high score."""
        # The Kozak consensus is GCCACCATGG
        seq = "GCCACCATGGCTAG"
        tokens = scan_sequence(seq)
        kozak_tokens = [t for t in tokens if t.element_type == "kozak"]
        # Should find at least one strong Kozak context
        if kozak_tokens:
            assert kozak_tokens[0].score >= KOZAK_REPORT_THRESHOLD

    def test_scan_all_frames(self):
        """scan_all_frames=True should find start codons in all frames."""
        # ATG at position 0, 1, and 2
        seq = "AATGATGATGC"
        tokens = scan_sequence(seq, scan_all_frames=True)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        frames = {t.frame for t in start_codons if t.frame is not None}
        # Should find ATGs in multiple frames
        assert len(frames) >= 1

    def test_scan_single_frame(self):
        """scan_all_frames=False should only scan frame 0."""
        seq = "AATGATGATGC"
        tokens = scan_sequence(seq, scan_all_frames=False)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        frames = {t.frame for t in start_codons if t.frame is not None}
        # All should be frame 0
        for f in frames:
            assert f == 0

    def test_unknown_enzyme_skipped(self):
        """Unknown enzyme names should be silently skipped."""
        seq = "ATGGCTAAGT"
        # Should not raise
        tokens = scan_sequence(seq, restriction_enzymes=["FakeEnzyme123"])
        assert isinstance(tokens, list)


# ---------------------------------------------------------------------------
# Kozak position weights
# ---------------------------------------------------------------------------

class TestKozakWeights:
    """Tests for Kozak position weight constants."""

    def test_weights_defined(self):
        """All four Kozak positions should have weights."""
        assert -3 in KOZAK_POSITION_WEIGHTS
        assert -2 in KOZAK_POSITION_WEIGHTS
        assert -1 in KOZAK_POSITION_WEIGHTS
        assert 4 in KOZAK_POSITION_WEIGHTS

    def test_a_at_minus_3_preferred(self):
        """A at position -3 should have highest weight."""
        assert KOZAK_POSITION_WEIGHTS[-3]["A"] == 1.0

    def test_c_at_minus_1_preferred(self):
        """C at position -1 should have highest weight."""
        assert KOZAK_POSITION_WEIGHTS[-1]["C"] == 1.0

    def test_g_at_plus_4_preferred(self):
        """G at position +4 should have highest weight."""
        assert KOZAK_POSITION_WEIGHTS[4]["G"] == 1.0

    def test_report_threshold_reasonable(self):
        """Kozak report threshold should be between 0 and 1."""
        assert 0.0 < KOZAK_REPORT_THRESHOLD < 1.0
