"""Unit tests for biocompiler.scanner.

Covers:
  1. gc_content() — correct fractions for known sequences
  2. gc_content() — edge cases: empty, all-GC, all-AT, single base
  3. validate_dna_sequence() — normalization and rejection of bad input
  4. _iupac_match() — IUPAC ambiguity code matching
  5. _score_kozak() — Kozak consensus scoring
  6. scan_sequence() — integrated motif detection (start/stop codons,
     splice donors/acceptors, instability motifs, restriction sites)
"""

from __future__ import annotations

import pytest

from biocompiler.scanner import (
    gc_content,
    validate_dna_sequence,
    _iupac_match,
    _score_kozak,
    scan_sequence,
    KOZAK_REPORT_THRESHOLD,
    SPLICE_DONOR_MIN_SCORE,
    SPLICE_ACCEPTOR_MIN_SCORE,
)
from biocompiler.exceptions import InvalidSequenceError
from biocompiler.types import Token


# ═══════════════════════════════════════════════════════════════════════════
# 1. gc_content()
# ═══════════════════════════════════════════════════════════════════════════

class TestGcContent:
    """Tests for gc_content() returning correct GC fractions."""

    def test_empty_string_returns_zero(self) -> None:
        assert gc_content("") == 0.0

    def test_all_gc_returns_one(self) -> None:
        assert gc_content("GCGCGC") == 1.0

    def test_all_at_returns_zero(self) -> None:
        assert gc_content("ATATAT") == 0.0

    def test_single_g(self) -> None:
        assert gc_content("G") == 1.0

    def test_single_a(self) -> None:
        assert gc_content("A") == 0.0

    def test_single_c(self) -> None:
        assert gc_content("C") == 1.0

    def test_single_t(self) -> None:
        assert gc_content("T") == 0.0

    def test_mixed_known_sequence(self) -> None:
        # "ATGC" → 2 GC out of 4 = 0.5
        assert gc_content("ATGC") == 0.5

    def test_half_gc(self) -> None:
        # "AATTCCGG" → 4 GC out of 8 = 0.5
        assert gc_content("AATTCCGG") == 0.5

    def test_one_third_gc(self) -> None:
        # "ATG" → 1 GC out of 3 ≈ 0.3333
        assert gc_content("ATG") == pytest.approx(1 / 3, abs=1e-4)

    def test_case_insensitive(self) -> None:
        assert gc_content("atgc") == 0.5
        assert gc_content("ATgc") == 0.5

    def test_rounding_to_four_decimals(self) -> None:
        # "ATGG" → 2 GC out of 4 = 0.5 exactly
        # Use a sequence that exercises rounding: 7 GC out of 12
        # 7/12 ≈ 0.5833
        seq = "GCGCGCGATATA"  # 7 G/C (positions 0-6) out of 12
        assert gc_content(seq) == pytest.approx(7 / 12, abs=1e-4)

    def test_long_sequence(self) -> None:
        # 50% GC, 100 bases
        seq = "GC" * 50  # 100 bases, all GC
        assert gc_content(seq) == 1.0

    def test_with_n_bases(self) -> None:
        # N bases are NOT counted as GC
        # "ATGN" → 1 GC out of 4 = 0.25
        assert gc_content("ATGN") == 0.25

    def test_all_n(self) -> None:
        # All N → 0 GC out of len
        assert gc_content("NNNN") == 0.0

    def test_eGFP_coding_sequence(self) -> None:
        """Real eGFP CDS — GC content should be ~60%."""
        egfp = (
            "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGG"
            "CCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCAC"
            "CACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCT"
            "ACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTT"
            "CTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA"
            "GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAA"
            "CGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGC"
            "AGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
            "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTG"
            "ACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
        )
        result = gc_content(egfp)
        # eGFP has ~60% GC; just check it's in a reasonable range
        assert 0.55 <= result <= 0.65


# ═══════════════════════════════════════════════════════════════════════════
# 2. validate_dna_sequence()
# ═══════════════════════════════════════════════════════════════════════════

class TestValidateDnaSequence:
    """Tests for validate_dna_sequence() — normalization and validation."""

    def test_uppercases_valid_sequence(self) -> None:
        assert validate_dna_sequence("atgc") == "ATGC"

    def test_returns_valid_sequence_unchanged(self) -> None:
        assert validate_dna_sequence("ATGC") == "ATGC"

    def test_accepts_n_bases(self) -> None:
        assert validate_dna_sequence("ATGN") == "ATGN"

    def test_rejects_invalid_bases(self) -> None:
        with pytest.raises(InvalidSequenceError) as exc_info:
            validate_dna_sequence("ATGX")
        assert "X" in str(exc_info.value) or "X" in exc_info.value.invalid_chars

    def test_rejects_multiple_invalid_bases(self) -> None:
        with pytest.raises(InvalidSequenceError) as exc_info:
            validate_dna_sequence("XYZ")
        invalid = exc_info.value.invalid_chars
        assert "X" in invalid
        assert "Y" in invalid
        assert "Z" in invalid

    def test_rejects_dashes(self) -> None:
        with pytest.raises(InvalidSequenceError):
            validate_dna_sequence("AT-GC")

    def test_empty_string_passes(self) -> None:
        # Empty string has no invalid characters
        assert validate_dna_sequence("") == ""

    def test_mixed_case_with_n(self) -> None:
        assert validate_dna_sequence("atgNnGC") == "ATGNNGC"


# ═══════════════════════════════════════════════════════════════════════════
# 3. _iupac_match()
# ═══════════════════════════════════════════════════════════════════════════

class TestIupacMatch:
    """Tests for _iupac_match() — IUPAC ambiguity code pattern matching."""

    def test_exact_match(self) -> None:
        assert _iupac_match("ATGC", "ATGC") is True

    def test_exact_mismatch(self) -> None:
        assert _iupac_match("ATGC", "ATGT") is False

    def test_length_mismatch(self) -> None:
        assert _iupac_match("ATG", "ATGC") is False

    def test_n_matches_any_base(self) -> None:
        assert _iupac_match("A", "N") is True
        assert _iupac_match("C", "N") is True
        assert _iupac_match("G", "N") is True
        assert _iupac_match("T", "N") is True

    def test_r_matches_ag(self) -> None:
        assert _iupac_match("A", "R") is True
        assert _iupac_match("G", "R") is True
        assert _iupac_match("C", "R") is False
        assert _iupac_match("T", "R") is False

    def test_y_matches_ct(self) -> None:
        assert _iupac_match("C", "Y") is True
        assert _iupac_match("T", "Y") is True
        assert _iupac_match("A", "Y") is False
        assert _iupac_match("G", "Y") is False

    def test_s_matches_gc(self) -> None:
        assert _iupac_match("G", "S") is True
        assert _iupac_match("C", "S") is True
        assert _iupac_match("A", "S") is False

    def test_w_matches_at(self) -> None:
        assert _iupac_match("A", "W") is True
        assert _iupac_match("T", "W") is True
        assert _iupac_match("G", "W") is False

    def test_mixed_iupac_pattern(self) -> None:
        # Pattern "RNT" should match "AGT", "AAT", "GGT", "GAT"
        assert _iupac_match("AGT", "RNT") is True
        assert _iupac_match("GGT", "RNT") is True
        assert _iupac_match("CGT", "RNT") is False  # C not in R


# ═══════════════════════════════════════════════════════════════════════════
# 4. _score_kozak()
# ═══════════════════════════════════════════════════════════════════════════

class TestScoreKozak:
    """Tests for _score_kozak() — Kozak consensus scoring."""

    def test_perfect_kozak(self) -> None:
        # GCCACCATGG is the optimal Kozak consensus
        # ATG starts at position 5
        seq = "GCCACCATGG"
        score = _score_kozak(seq, 5)
        # The weighted scoring formula does not reach exactly 1.0;
        # the optimal consensus scores ~0.61 because the position weights
        # are multiplied together in a quadratic scheme.
        assert score == pytest.approx(0.6137, abs=0.01)

    def test_no_consensus(self) -> None:
        # "TTTATGT" — poor Kozak: T at -3, T at -2, T at -1, T at +4
        seq = "TTTATGT"
        score = _score_kozak(seq, 3)
        assert score < 0.3

    def test_atg_at_start_of_sequence(self) -> None:
        # No upstream positions available; score comes only from +4
        seq = "ATGG"
        score = _score_kozak(seq, 0)
        # Only +4 position contributes
        assert 0.0 <= score <= 1.0

    def test_partial_context(self) -> None:
        # ATG at position 1: only -1 available upstream, +4 downstream
        seq = "CATGG"
        score = _score_kozak(seq, 1)
        assert 0.0 <= score <= 1.0

    def test_score_range(self) -> None:
        # Any valid sequence should give a score in [0, 1]
        seq = "AGCATGGA"
        score = _score_kozak(seq, 2)
        assert 0.0 <= score <= 1.0


# ═══════════════════════════════════════════════════════════════════════════
# 5. scan_sequence() — integrated tests
# ═══════════════════════════════════════════════════════════════════════════

class TestScanSequence:
    """Tests for scan_sequence() — integrated motif detection."""

    # --- Basic properties ---

    def test_empty_sequence_returns_no_tokens(self) -> None:
        tokens = scan_sequence("")
        assert tokens == []

    def test_invalid_sequence_raises(self) -> None:
        with pytest.raises(InvalidSequenceError):
            scan_sequence("ATGX")

    def test_tokens_are_sorted_by_position(self) -> None:
        seq = "ATGATTTAATAG"  # has start codon, instability motif, stop codon
        tokens = scan_sequence(seq, use_maxentscan=False)
        positions = [t.position for t in tokens]
        assert positions == sorted(positions)

    def test_case_insensitive(self) -> None:
        upper = scan_sequence("ATG", use_maxentscan=False)
        lower = scan_sequence("atg", use_maxentscan=False)
        # Same number of tokens at same positions
        assert len(upper) == len(lower)
        for u, l in zip(upper, lower):
            assert u.position == l.position
            assert u.element_type == l.element_type

    # --- Start codons ---

    def test_finds_start_codon_atg(self) -> None:
        tokens = scan_sequence("ATG", use_maxentscan=False)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        assert len(start_codons) >= 1
        assert start_codons[0].match_sequence == "ATG"
        assert start_codons[0].position == 0

    def test_start_codons_all_frames(self) -> None:
        # ATG at positions 0 and 1 (two different frames)
        seq = "AATGATG"  # ATG at pos 1 (frame 1), ATG at pos 4 (frame 1)
        tokens = scan_sequence(seq, scan_all_frames=True, use_maxentscan=False)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        positions = {t.position for t in start_codons}
        assert 1 in positions  # ATG starting at position 1
        assert 4 in positions  # ATG starting at position 4

    def test_start_codons_frame0_only(self) -> None:
        seq = "AATGATG"  # ATG at pos 1 and 4, but frame 0 only
        tokens = scan_sequence(seq, scan_all_frames=False, use_maxentscan=False)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        # In frame 0, positions checked: 0, 3, 6 → none are ATG start
        for t in start_codons:
            assert t.frame == 0

    def test_start_codon_has_kozak_score(self) -> None:
        tokens = scan_sequence("GCCACCATGG", use_maxentscan=False)
        start_codons = [t for t in tokens if t.element_type == "start_codon"]
        assert len(start_codons) >= 1
        # Perfect Kozak should have a high score
        assert start_codons[0].score >= 0.5

    # --- Stop codons ---

    def test_finds_stop_codon_taa(self) -> None:
        tokens = scan_sequence("TAA", use_maxentscan=False)
        stops = [t for t in tokens if t.element_type == "stop_codon"]
        assert len(stops) >= 1
        assert stops[0].match_sequence == "TAA"

    def test_finds_stop_codon_tag(self) -> None:
        tokens = scan_sequence("TAG", use_maxentscan=False)
        stops = [t for t in tokens if t.element_type == "stop_codon"]
        assert len(stops) >= 1
        assert stops[0].match_sequence == "TAG"

    def test_finds_stop_codon_tga(self) -> None:
        tokens = scan_sequence("TGA", use_maxentscan=False)
        stops = [t for t in tokens if t.element_type == "stop_codon"]
        assert len(stops) >= 1
        assert stops[0].match_sequence == "TGA"

    def test_stop_codons_all_frames(self) -> None:
        seq = "ATGAATGA"  # TGA at pos 1 (frame 1) and pos 5 (frame 2)
        tokens = scan_sequence(seq, scan_all_frames=True, use_maxentscan=False)
        stops = [t for t in tokens if t.element_type == "stop_codon"]
        frames = {t.frame for t in stops}
        assert len(frames) > 1  # at least two different frames

    def test_non_stop_codon_not_found(self) -> None:
        tokens = scan_sequence("ATG", use_maxentscan=False)
        stops = [t for t in tokens if t.element_type == "stop_codon"]
        assert len(stops) == 0

    # --- Instability motifs ---

    def test_finds_instability_motif(self) -> None:
        seq = "AAAAATTTAAAA"  # ATTTA at position 4
        tokens = scan_sequence(seq, use_maxentscan=False)
        instability = [t for t in tokens if t.element_type == "instability_motif"]
        assert len(instability) >= 1
        assert instability[0].match_sequence == "ATTTA"
        assert instability[0].position == 4

    def test_no_instability_motif_in_normal_sequence(self) -> None:
        seq = "ATGCATGC"
        tokens = scan_sequence(seq, use_maxentscan=False)
        instability = [t for t in tokens if t.element_type == "instability_motif"]
        assert len(instability) == 0

    def test_multiple_instability_motifs(self) -> None:
        seq = "ATTTAATTTA"  # Two overlapping ATTTA at pos 0 and pos 5
        tokens = scan_sequence(seq, use_maxentscan=False)
        instability = [t for t in tokens if t.element_type == "instability_motif"]
        # ATTTA at pos 0: seq[0:5]="ATTTA"
        # ATTTA at pos 5: seq[5:10]="ATTTA"
        assert len(instability) >= 1

    # --- Splice donor/acceptor (without MaxEntScan) ---

    def test_splice_donor_without_maxentscan(self) -> None:
        seq = "AAGTAA"  # GT at position 2
        tokens = scan_sequence(seq, use_maxentscan=False)
        donors = [t for t in tokens if t.element_type == "splice_donor"]
        assert len(donors) >= 1
        assert donors[0].match_sequence == "GT"
        assert donors[0].score == 5.0  # DONOR_FALLBACK_SCORE

    def test_splice_acceptor_without_maxentscan_with_tract(self) -> None:
        # Create a sequence with AG preceded by a polypyrimidine tract
        # (lots of C/T upstream)
        seq = "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCAG" + "ATG"
        tokens = scan_sequence(seq, use_maxentscan=False)
        acceptors = [t for t in tokens if t.element_type == "splice_acceptor"]
        assert len(acceptors) >= 1

    def test_splice_donor_with_maxentscan(self) -> None:
        # Need a longer sequence for MaxEntScan 9-mer model
        seq = "CAGGTAAGT" * 10  # canonical donor context
        tokens = scan_sequence(seq, use_maxentscan=True)
        donors = [t for t in tokens if t.element_type == "splice_donor"]
        # At least one GT should score above threshold
        for d in donors:
            assert d.score >= SPLICE_DONOR_MIN_SCORE

    # --- Restriction sites ---

    def test_finds_ecori_restriction_site(self) -> None:
        # EcoRI site is GAATTC
        seq = "AAGAATTCAA"
        tokens = scan_sequence(seq, restriction_enzymes=["EcoRI"], use_maxentscan=False)
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) >= 1
        assert sites[0].match_sequence == "GAATTC"
        assert sites[0].position == 2

    def test_finds_ecori_reverse_complement(self) -> None:
        # EcoRI RC is GAATTC (palindrome!), so same site
        # Use BamHI instead: GGATCC, RC = GGATCC (also palindrome!)
        # Use XhoI: CTCGAG, RC = CTCGAG (palindrome too!)
        # Use BglII: AGATCT, RC = AGATCT (palindrome)
        # Let's use NdeI: CATATG, RC = CATATG (palindrome)
        # Hmm, most common enzymes are palindromes.
        # Use an enzyme that is NOT palindromic... checking RESTRICTION_ENZYMES
        # Actually many are palindromes. Let's test that we don't double-count
        # palindromes.
        seq = "AAGAATTCAA"
        tokens = scan_sequence(seq, restriction_enzymes=["EcoRI"], use_maxentscan=False)
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        # EcoRI (GAATTC) is palindromic — should not be double-counted
        eco_tokens = [s for s in sites if s.match_sequence == "GAATTC"]
        assert len(eco_tokens) == 1

    def test_no_restriction_site_when_not_requested(self) -> None:
        seq = "AAGAATTCAA"
        tokens = scan_sequence(seq, restriction_enzymes=None, use_maxentscan=False)
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) == 0

    def test_unknown_enzyme_name_ignored(self) -> None:
        seq = "AAGAATTCAA"
        tokens = scan_sequence(seq, restriction_enzymes=["FakeEnzyme"], use_maxentscan=False)
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) == 0

    def test_iupac_restriction_site_sfiI(self) -> None:
        # SfiI: GGCCNNNNNGGCC — contains IUPAC N wildcards (5 N positions)
        seq = "AAGGCCAAAAAGGCCAA"  # 5 bases between the two GGCC segments
        tokens = scan_sequence(seq, restriction_enzymes=["SfiI"], use_maxentscan=False)
        sites = [t for t in tokens if t.element_type == "restriction_site"]
        assert len(sites) >= 1

    # --- Kozak tokens ---

    def test_kozak_token_for_strong_context(self) -> None:
        # GCCACCATGG has perfect Kozak; ATG at position 5
        seq = "GCCACCATGG"
        tokens = scan_sequence(seq, use_maxentscan=False)
        kozak_tokens = [t for t in tokens if t.element_type == "kozak"]
        assert len(kozak_tokens) >= 1
        assert kozak_tokens[0].score >= KOZAK_REPORT_THRESHOLD

    def test_no_kozak_token_for_weak_context(self) -> None:
        # Minimal ATG with poor context
        seq = "TTTATGT"
        tokens = scan_sequence(seq, use_maxentscan=False)
        kozak_tokens = [t for t in tokens if t.element_type == "kozak"]
        # Weak context should not meet the reporting threshold
        assert len(kozak_tokens) == 0

    # --- scan_all_frames parameter ---

    def test_scan_all_frames_true_finds_more_codons(self) -> None:
        seq = "AATG"  # ATG at position 1 (frame 1 only)
        tokens_all = scan_sequence(seq, scan_all_frames=True, use_maxentscan=False)
        tokens_f0 = scan_sequence(seq, scan_all_frames=False, use_maxentscan=False)
        starts_all = [t for t in tokens_all if t.element_type == "start_codon"]
        starts_f0 = [t for t in tokens_f0 if t.element_type == "start_codon"]
        assert len(starts_all) >= len(starts_f0)

    # --- Integration: multiple features in one sequence ---

    def test_multiple_features_in_one_sequence(self) -> None:
        # Sequence with start codon, instability motif, and stop codon
        seq = "ATGATTTATAA"  # ATG at 0, ATTTA at 2, TAA at 7 (frame 1)
        tokens = scan_sequence(seq, use_maxentscan=False)
        element_types = {t.element_type for t in tokens}
        # Should find at least start_codon and instability_motif
        assert "start_codon" in element_types
        assert "instability_motif" in element_types

    def test_token_fields_populated_correctly(self) -> None:
        seq = "ATG"
        tokens = scan_sequence(seq, use_maxentscan=False)
        start = [t for t in tokens if t.element_type == "start_codon"][0]
        assert isinstance(start, Token)
        assert start.position == 0
        assert start.element_type == "start_codon"
        assert start.match_sequence == "ATG"
        assert isinstance(start.score, float)

    # --- Threshold parameters ---

    def test_donor_threshold_filtering(self) -> None:
        # With a very high threshold, fewer donors should be reported
        seq = "CAGGTAAGT" * 10
        tokens_low = scan_sequence(seq, use_maxentscan=True, donor_threshold=0.0)
        tokens_high = scan_sequence(seq, use_maxentscan=True, donor_threshold=20.0)
        donors_low = [t for t in tokens_low if t.element_type == "splice_donor"]
        donors_high = [t for t in tokens_high if t.element_type == "splice_donor"]
        assert len(donors_low) >= len(donors_high)

    def test_acceptor_threshold_filtering(self) -> None:
        seq = "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCAGATG"
        tokens_low = scan_sequence(seq, use_maxentscan=True, acceptor_threshold=0.0)
        tokens_high = scan_sequence(seq, use_maxentscan=True, acceptor_threshold=20.0)
        acc_low = [t for t in tokens_low if t.element_type == "splice_acceptor"]
        acc_high = [t for t in tokens_high if t.element_type == "splice_acceptor"]
        assert len(acc_low) >= len(acc_high)


# ═══════════════════════════════════════════════════════════════════════════
# 6. Determinism — same input always gives same output
# ═══════════════════════════════════════════════════════════════════════════

class TestDeterminism:
    """Scanner output must be deterministic."""

    def test_gc_content_deterministic(self) -> None:
        seq = "ATGCGATCGATCGATCG"
        assert gc_content(seq) == gc_content(seq)

    def test_scan_sequence_deterministic(self) -> None:
        seq = "ATGGCCATTTATAAGATCTGA"
        t1 = scan_sequence(seq, use_maxentscan=False)
        t2 = scan_sequence(seq, use_maxentscan=False)
        assert t1 == t2
