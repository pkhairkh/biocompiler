"""Comprehensive pytest tests for biocompiler.solver.automaton.

Tests the DFA (Deterministic Finite Automaton) construction and acceptance
functions for forbidden-pattern avoidance in DNA sequences.

Test categories:
  1. Single-pattern DFA — build and verify rejection/acceptance
  2. Composite DFA — multiple forbidden patterns
  3. Reverse-complement DFA — palindromic and non-palindromic enzymes
  4. DFA acceptance testing — dfa_accepts function
  5. Edge cases — short/empty sequences, overlapping patterns, single-char
  6. Alphabet encoding — BASE_MAP correctness
  7. DFA structural properties
  8. Integration / real-world scenarios
"""

from __future__ import annotations

import pytest

from biocompiler.solver.automaton import (
    BASE_MAP,
    build_composite_dfa,
    build_forbidden_pattern_dfa,
    build_reverse_complement_dfa,
    build_trun_dfa,
    dfa_accepts,
    negate_dfa,
)

from biocompiler.constants import BASE_REV, reverse_complement


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def ecori_dfa():
    """DFA that forbids the EcoRI recognition site GAATTC."""
    return build_forbidden_pattern_dfa("GAATTC")


@pytest.fixture
def bamhi_dfa():
    """DFA that forbids the BamHI recognition site GGATCC."""
    return build_forbidden_pattern_dfa("GGATCC")


@pytest.fixture
def composite_dfa():
    """Composite DFA forbidding GAATTC, GGATCC, and CTCGAG."""
    return build_composite_dfa(["GAATTC", "GGATCC", "CTCGAG"])


@pytest.fixture
def composite_dfa_with_rc():
    """Composite DFA forbidding GAATTC, GGATCC, CTCGAG and their RCs."""
    patterns = ["GAATTC", "GGATCC", "CTCGAG"]
    all_patterns = list(set(patterns + [reverse_complement(p) for p in patterns]))
    return build_composite_dfa(all_patterns)


@pytest.fixture
def single_char_dfa():
    """DFA that forbids the single character 'G'."""
    return build_forbidden_pattern_dfa("G")


@pytest.fixture
def aaa_dfa():
    """DFA that forbids the triple-A pattern 'AAA'."""
    return build_forbidden_pattern_dfa("AAA")


# ═══════════════════════════════════════════════════════════════════════════
# 1. Single-pattern DFA
# ═══════════════════════════════════════════════════════════════════════════

class TestSinglePatternDFA:
    """Tests for build_forbidden_pattern_dfa with a single forbidden pattern."""

    # --- EcoRI (GAATTC) rejection -----------------------------------------

    def test_ecori_rejects_exact_pattern(self, ecori_dfa):
        """DFA should reject a string that IS the forbidden pattern."""
        assert not dfa_accepts(ecori_dfa, "GAATTC")

    def test_ecori_rejects_pattern_at_start(self, ecori_dfa):
        """DFA should reject when pattern appears at the beginning."""
        assert not dfa_accepts(ecori_dfa, "GAATTCATCGATCG")

    def test_ecori_rejects_pattern_at_middle(self, ecori_dfa):
        """DFA should reject when pattern appears in the middle."""
        assert not dfa_accepts(ecori_dfa, "ATCGAATTCGATCG")

    def test_ecori_rejects_pattern_at_end(self, ecori_dfa):
        """DFA should reject when pattern appears at the end."""
        assert not dfa_accepts(ecori_dfa, "ATCGATCGAATTC")

    def test_ecori_accepts_clean_sequence(self, ecori_dfa):
        """DFA should accept a sequence that does NOT contain GAATTC."""
        assert dfa_accepts(ecori_dfa, "ATCGATCG")

    def test_ecori_accepts_all_a_sequence(self, ecori_dfa):
        """An all-A sequence does not contain GAATTC."""
        assert dfa_accepts(ecori_dfa, "AAAAAAAAAA")

    def test_ecori_accepts_partial_prefix_only(self, ecori_dfa):
        """Partial match of the pattern prefix should NOT be rejected."""
        # 'GAATT' is a prefix of 'GAATTC' but not the full pattern
        assert dfa_accepts(ecori_dfa, "GAATT")

    def test_ecori_accepts_partial_prefix_with_mismatch(self, ecori_dfa):
        """Prefix followed by a mismatch should NOT be rejected."""
        # 'GAATTA' — 5 out of 6 characters match, but last is 'A' not 'C'
        assert dfa_accepts(ecori_dfa, "GAATTA")

    def test_ecori_rejects_multiple_occurrences(self, ecori_dfa):
        """DFA should reject even if the pattern appears more than once."""
        assert not dfa_accepts(ecori_dfa, "GAATTCGAATTC")

    def test_ecori_rejects_overlapping_occurrence(self, ecori_dfa):
        """Overlapping partial matches should still detect the full pattern.

        For GAATTC there are no true overlaps (the pattern is not
        self-overlapping), but this is a sanity check.
        """
        # Embed the pattern in a longer string
        assert not dfa_accepts(ecori_dfa, "AAAGAATTCAAA")

    # --- BamHI (GGATCC) rejection -----------------------------------------

    def test_bamhi_rejects_exact_pattern(self, bamhi_dfa):
        """DFA should reject the BamHI site GGATCC."""
        assert not dfa_accepts(bamhi_dfa, "GGATCC")

    def test_bamhi_rejects_pattern_in_context(self, bamhi_dfa):
        """DFA should reject when BamHI site is embedded."""
        assert not dfa_accepts(bamhi_dfa, "ATCGGATCCGATC")

    def test_bamhi_accepts_clean_sequence(self, bamhi_dfa):
        """DFA should accept sequences without GGATCC."""
        assert dfa_accepts(bamhi_dfa, "ATCGATCG")

    def test_bamhi_accepts_partial_gg(self, bamhi_dfa):
        """'GG' alone is not a BamHI site."""
        assert dfa_accepts(bamhi_dfa, "GGATCA")

    # --- General single-pattern properties --------------------------------

    def test_dfa_has_correct_initial_state(self, ecori_dfa):
        """The initial state is always 0 by convention."""
        transitions, accepting = ecori_dfa
        # Initial state is 0 — it must be an accepting state
        # (no pattern matched yet)
        assert len(transitions) > 0  # state 0 exists
        assert 0 in accepting

    def test_dfa_has_transitions(self, ecori_dfa):
        """The DFA should have at least one transition row."""
        transitions, _ = ecori_dfa
        assert len(transitions) > 0

    def test_dfa_has_forbidden_states(self, ecori_dfa):
        """The DFA should have at least one forbidden (non-accepting/trap) state."""
        transitions, accepting = ecori_dfa
        all_states = set(range(len(transitions)))
        forbidden = all_states - set(accepting)
        assert len(forbidden) > 0

    def test_pattern_length_equals_forbidden_state(self, ecori_dfa):
        """For a single pattern of length n, state n is the forbidden (trap) state."""
        transitions, accepting = ecori_dfa
        assert 6 not in accepting  # GAATTC has length 6; state 6 is the trap

    def test_build_forbidden_pattern_dfa_raises_on_invalid_base(self):
        """build_forbidden_pattern_dfa should raise ValueError for invalid bases."""
        with pytest.raises(ValueError, match="Invalid base"):
            build_forbidden_pattern_dfa("GAXTC")

    def test_build_forbidden_pattern_dfa_accepts_lowercase(self):
        """The DFA builder should accept lowercase input (uppercased internally)."""
        dfa = build_forbidden_pattern_dfa("gaattc")
        assert not dfa_accepts(dfa, "GAATTC")

    def test_ecori_near_miss_gaatta(self, ecori_dfa):
        """GAATTA differs from GAATTC in the last position — should be accepted."""
        assert dfa_accepts(ecori_dfa, "GAATTA")

    def test_ecori_near_miss_gaattc_with_prefix(self, ecori_dfa):
        """GAATTC preceded by its own prefix should still be rejected."""
        assert not dfa_accepts(ecori_dfa, "GAGAATTC")


# ═══════════════════════════════════════════════════════════════════════════
# 2. Composite DFA (multiple forbidden patterns)
# ═══════════════════════════════════════════════════════════════════════════

class TestCompositeDFA:
    """Tests for build_composite_dfa with multiple forbidden patterns."""

    def test_rejects_ecori(self, composite_dfa):
        """Composite DFA should reject sequences containing GAATTC."""
        assert not dfa_accepts(composite_dfa, "GAATTC")

    def test_rejects_bamhi(self, composite_dfa):
        """Composite DFA should reject sequences containing GGATCC."""
        assert not dfa_accepts(composite_dfa, "GGATCC")

    def test_rejects_xhoi(self, composite_dfa):
        """Composite DFA should reject sequences containing CTCGAG."""
        assert not dfa_accepts(composite_dfa, "CTCGAG")

    def test_accepts_clean_sequence(self, composite_dfa):
        """Composite DFA should accept sequences with none of the patterns."""
        assert dfa_accepts(composite_dfa, "ATCGATCG")

    def test_accepts_all_a_sequence(self, composite_dfa):
        """An all-A sequence contains none of the forbidden patterns."""
        assert dfa_accepts(composite_dfa, "AAAAAAAAAA")

    def test_rejects_first_pattern_in_long_sequence(self, composite_dfa):
        """Rejects a long sequence containing GAATTC."""
        assert not dfa_accepts(composite_dfa, "ATCGAATTCGATCGATCG")

    def test_rejects_second_pattern_in_long_sequence(self, composite_dfa):
        """Rejects a long sequence containing GGATCC."""
        assert not dfa_accepts(composite_dfa, "ATCGATCGGATCCGATC")

    def test_rejects_third_pattern_in_long_sequence(self, composite_dfa):
        """Rejects a long sequence containing CTCGAG."""
        assert not dfa_accepts(composite_dfa, "ATCCTCGAGGATCG")

    def test_rejects_sequence_with_multiple_patterns(self, composite_dfa):
        """Rejects a sequence containing multiple forbidden patterns."""
        assert not dfa_accepts(composite_dfa, "GAATTCGGATCC")

    def test_accepts_partial_match(self, composite_dfa):
        """Partial matches of any forbidden pattern should NOT be rejected."""
        # GAATT is a prefix of GAATTC but not the full pattern
        assert dfa_accepts(composite_dfa, "GAATT")
        assert dfa_accepts(composite_dfa, "GGATC")
        assert dfa_accepts(composite_dfa, "CTCGA")

    def test_composite_without_rc_excludes_rc(self):
        """When only the forward pattern is given, RC patterns are not forbidden."""
        # 'AACG' → RC = 'CGTT'
        dfa = build_composite_dfa(["AACG"])
        assert not dfa_accepts(dfa, "AACG")  # original is rejected
        assert dfa_accepts(dfa, "CGTT")       # RC is accepted (not forbidden)

    def test_composite_with_rc_includes_rc(self):
        """When both pattern and RC are included, both are forbidden."""
        # 'AACG' → RC = 'CGTT'
        dfa = build_reverse_complement_dfa("AACG")
        assert not dfa_accepts(dfa, "AACG")   # original is rejected
        assert not dfa_accepts(dfa, "CGTT")    # RC is also rejected

    def test_composite_palindrome_rc_not_duplicated(self):
        """For palindromic patterns, RC = pattern, so no duplication effect."""
        # GAATTC RC = GAATTC (palindrome)
        dfa = build_reverse_complement_dfa("GAATTC")
        assert not dfa_accepts(dfa, "GAATTC")

    def test_empty_patterns_gives_trivial_dfa(self):
        """An empty pattern list should produce a DFA that accepts everything."""
        dfa = build_composite_dfa([])
        assert dfa_accepts(dfa, "GAATTC")
        assert dfa_accepts(dfa, "GGATCC")

    def test_iupac_patterns_raise_error(self):
        """Patterns with IUPAC ambiguity codes (e.g. 'N') should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid base"):
            build_composite_dfa(["GAATTC", "GATNNN"])


# ═══════════════════════════════════════════════════════════════════════════
# 3. Reverse-complement DFA
# ═══════════════════════════════════════════════════════════════════════════

class TestReverseComplementDFA:
    """Tests for reverse complement handling in DFAs."""

    # --- Palindromic enzymes ----------------------------------------------

    def test_ecori_rc_is_same(self):
        """EcoRI site GAATTC is a palindrome: RC(GAATTC) == GAATTC."""
        assert reverse_complement("GAATTC") == "GAATTC"

    def test_xhoi_rc_is_same(self):
        """XhoI site CTCGAG is a palindrome: RC(CTCGAG) == CTCGAG."""
        assert reverse_complement("CTCGAG") == "CTCGAG"

    def test_ecori_rc_dfa_rejects_original(self):
        """RC DFA for EcoRI rejects the original pattern (palindrome)."""
        rc_dfa = build_reverse_complement_dfa("GAATTC")
        assert not dfa_accepts(rc_dfa, "GAATTC")

    def test_ecori_rc_dfa_rejects_rc(self):
        """RC DFA for EcoRI rejects its RC (same as original)."""
        rc_dfa = build_reverse_complement_dfa("GAATTC")
        rc = reverse_complement("GAATTC")
        assert not dfa_accepts(rc_dfa, rc)

    def test_xhoi_rc_dfa_rejects_original(self):
        """RC DFA for XhoI rejects the original pattern (palindrome)."""
        rc_dfa = build_reverse_complement_dfa("CTCGAG")
        assert not dfa_accepts(rc_dfa, "CTCGAG")

    def test_xhoi_rc_dfa_rejects_rc(self):
        """RC DFA for XhoI rejects its RC (same as original)."""
        rc_dfa = build_reverse_complement_dfa("CTCGAG")
        rc = reverse_complement("CTCGAG")
        assert not dfa_accepts(rc_dfa, rc)

    # --- Non-palindromic enzymes ------------------------------------------

    def test_non_palindromic_rc_differs(self):
        """For non-palindromic patterns, RC != original."""
        assert reverse_complement("AACG") == "CGTT"
        assert reverse_complement("AACG") != "AACG"

    def test_rc_dfa_rejects_rc_pattern(self):
        """RC DFA for a non-palindromic pattern rejects the RC."""
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert not dfa_accepts(rc_dfa, "CGTT")

    def test_rc_dfa_rejects_original_for_non_palindrome(self):
        """RC DFA for a non-palindromic pattern also rejects the original.

        build_reverse_complement_dfa builds a composite DFA that rejects
        BOTH the pattern and its reverse complement.
        """
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert not dfa_accepts(rc_dfa, "AACG")

    def test_composite_rejects_both_original_and_rc(self):
        """RC DFA should reject both the pattern and its RC."""
        dfa = build_reverse_complement_dfa("AACG")
        assert not dfa_accepts(dfa, "AACG")
        assert not dfa_accepts(dfa, "CGTT")

    def test_composite_accepts_neither_original_nor_rc(self):
        """RC DFA accepts sequences with neither pattern nor RC."""
        dfa = build_reverse_complement_dfa("AACG")
        assert dfa_accepts(dfa, "TTTTTTTT")

    # --- Other restriction enzyme RCs -------------------------------------

    def test_hindiii_palindrome(self):
        """HindIII AAGCTT → RC AAGCTT (palindrome)."""
        assert reverse_complement("AAGCTT") == "AAGCTT"

    def test_sali_palindrome(self):
        """SalI GTCGAC → RC GTCGAC (palindrome)."""
        assert reverse_complement("GTCGAC") == "GTCGAC"

    def test_multiple_restriction_sites_with_rc(self):
        """Composite DFA for multiple enzymes with RC rejects all variants."""
        patterns = ["GAATTC", "AACG"]
        # Include reverse complements explicitly
        all_patterns = list(set(patterns + [reverse_complement(p) for p in patterns]))
        dfa = build_composite_dfa(all_patterns)
        # Palindrome
        assert not dfa_accepts(dfa, "GAATTC")
        # Non-palindrome + RC
        assert not dfa_accepts(dfa, "AACG")
        assert not dfa_accepts(dfa, "CGTT")
        # Clean sequence
        assert dfa_accepts(dfa, "TTTTAAAA")


# ═══════════════════════════════════════════════════════════════════════════
# 4. DFA acceptance testing
# ═══════════════════════════════════════════════════════════════════════════

class TestDfaAccepts:
    """Tests for the dfa_accepts function."""

    def test_accepts_clean_sequence(self, ecori_dfa):
        """dfa_accepts returns True for a sequence without the forbidden pattern."""
        assert dfa_accepts(ecori_dfa, "ATCGATCG") is True

    def test_rejects_contaminated_sequence(self, ecori_dfa):
        """dfa_accepts returns False for a sequence with the forbidden pattern."""
        assert dfa_accepts(ecori_dfa, "GAATTC") is False

    def test_rejects_embedded_pattern(self, ecori_dfa):
        """dfa_accepts returns False when the pattern is embedded."""
        assert dfa_accepts(ecori_dfa, "AAGAATTCGG") is False

    def test_accepts_after_near_miss(self, ecori_dfa):
        """dfa_accepts returns True after a near-miss that doesn't complete the pattern."""
        # 'GAATTA' almost matches GAATTC but diverges at position 6
        assert dfa_accepts(ecori_dfa, "GAATTA") is True

    def test_accepts_with_composite_dfa(self, composite_dfa):
        """dfa_accepts works correctly with composite DFAs."""
        assert dfa_accepts(composite_dfa, "ATCGATCG") is True
        assert dfa_accepts(composite_dfa, "GAATTC") is False
        assert dfa_accepts(composite_dfa, "GGATCC") is False
        assert dfa_accepts(composite_dfa, "CTCGAG") is False

    def test_accepts_returns_bool(self, ecori_dfa):
        """dfa_accepts should return a boolean, not a truthy/falsy value."""
        result = dfa_accepts(ecori_dfa, "ATCG")
        assert isinstance(result, bool)

    def test_case_insensitive_sequence(self, ecori_dfa):
        """dfa_accepts should handle lowercase sequences (uppercased internally)."""
        assert dfa_accepts(ecori_dfa, "atcgatcg") is True
        assert dfa_accepts(ecori_dfa, "gaattc") is False

    def test_long_clean_sequence(self, ecori_dfa):
        """A long random-like sequence without the pattern should be accepted."""
        seq = "ATCGATCGATCGATCGATCGATCGATCGATCG"
        assert dfa_accepts(ecori_dfa, seq) is True

    def test_pattern_spans_codon_boundary(self, ecori_dfa):
        """Pattern spanning a codon boundary should still be detected."""
        # GAATTC could span: ...XXG|AAT|TCX...
        assert dfa_accepts(ecori_dfa, "AAGAATTCGA") is False

    def test_rejects_early_return(self, ecori_dfa):
        """dfa_accepts should return False as soon as the pattern is found."""
        # Even though the rest of the string is fine, the DFA should have
        # already transitioned to the forbidden state
        assert dfa_accepts(ecori_dfa, "GAATTCATCGATCG") is False


# ═══════════════════════════════════════════════════════════════════════════
# 5. Edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    """Edge-case tests for DFA construction and acceptance."""

    # --- Single-character pattern -----------------------------------------

    def test_single_char_pattern_rejects_matching(self, single_char_dfa):
        """A single-character pattern should reject all strings containing that char."""
        assert not dfa_accepts(single_char_dfa, "G")

    def test_single_char_pattern_rejects_in_context(self, single_char_dfa):
        """Any string with 'G' should be rejected."""
        assert not dfa_accepts(single_char_dfa, "ATCG")
        assert not dfa_accepts(single_char_dfa, "GATC")
        assert not dfa_accepts(single_char_dfa, "ATCGATCG")

    def test_single_char_pattern_accepts_without(self, single_char_dfa):
        """Strings without 'G' should be accepted."""
        assert dfa_accepts(single_char_dfa, "ATC")
        assert dfa_accepts(single_char_dfa, "AAAA")
        assert dfa_accepts(single_char_dfa, "TCTC")

    # --- Very short sequences --------------------------------------------

    def test_sequence_shorter_than_pattern(self, ecori_dfa):
        """A sequence shorter than the pattern length should be accepted."""
        assert dfa_accepts(ecori_dfa, "GAAT")     # length 4 < 6
        assert dfa_accepts(ecori_dfa, "GAA")       # length 3 < 6
        assert dfa_accepts(ecori_dfa, "GA")         # length 2 < 6

    def test_sequence_length_one(self, ecori_dfa):
        """A single-character sequence cannot contain a 6-char pattern."""
        assert dfa_accepts(ecori_dfa, "G")
        assert dfa_accepts(ecori_dfa, "A")
        assert dfa_accepts(ecori_dfa, "T")
        assert dfa_accepts(ecori_dfa, "C")

    # --- Empty sequence ---------------------------------------------------

    def test_empty_sequence_accepted(self, ecori_dfa):
        """An empty sequence should be accepted (no pattern can be found)."""
        assert dfa_accepts(ecori_dfa, "")

    def test_empty_sequence_composite(self, composite_dfa):
        """Empty sequence should be accepted by composite DFA."""
        assert dfa_accepts(composite_dfa, "")

    # --- Pattern equals sequence -----------------------------------------

    def test_pattern_equals_sequence(self, ecori_dfa):
        """When the sequence IS the pattern, it should be rejected."""
        assert not dfa_accepts(ecori_dfa, "GAATTC")

    def test_single_char_equals_sequence(self):
        """When pattern and sequence are both 'G', should be rejected."""
        dfa = build_forbidden_pattern_dfa("G")
        assert not dfa_accepts(dfa, "G")

    # --- Pattern appears multiple times -----------------------------------

    def test_pattern_appears_twice(self, ecori_dfa):
        """Two occurrences of the pattern should still be rejected."""
        assert not dfa_accepts(ecori_dfa, "GAATTCGAATTC")

    def test_pattern_appears_three_times(self, ecori_dfa):
        """Three occurrences should still be rejected."""
        assert not dfa_accepts(ecori_dfa, "GAATTCGAATTCGAATTC")

    def test_pattern_separated_by_other_bases(self, ecori_dfa):
        """Pattern occurrences separated by other bases should be rejected."""
        assert not dfa_accepts(ecori_dfa, "GAATTCATCGAATTC")

    # --- Overlapping patterns ---------------------------------------------

    def test_overlapping_aaa_in_aaaaaa(self, aaa_dfa):
        """'AAAAAA' contains overlapping 'AAA' — should be rejected."""
        assert not dfa_accepts(aaa_dfa, "AAAAAA")

    def test_overlapping_aaa_in_aaaa(self, aaa_dfa):
        """'AAAA' contains 'AAA' starting at position 0 and 1."""
        assert not dfa_accepts(aaa_dfa, "AAAA")

    def test_overlapping_aaa_in_aaa(self, aaa_dfa):
        """'AAA' is exactly the pattern."""
        assert not dfa_accepts(aaa_dfa, "AAA")

    def test_overlapping_aaa_accepts_aa(self, aaa_dfa):
        """'AA' is too short to contain 'AAA'."""
        assert dfa_accepts(aaa_dfa, "AA")

    def test_overlapping_gatc_in_ggatcc(self):
        """'GGATCC' contains 'GATC' starting at position 1."""
        dfa = build_forbidden_pattern_dfa("GATC")
        assert not dfa_accepts(dfa, "GGATCC")

    def test_overlapping_gatc_at_boundaries(self):
        """'AGATCG' contains 'GATC' starting at position 1."""
        dfa = build_forbidden_pattern_dfa("GATC")
        assert not dfa_accepts(dfa, "AGATCG")

    # --- Pattern with repeated subpatterns --------------------------------

    def test_pattern_with_repeated_subpattern(self):
        """Pattern 'GATAGATC' contains 'GATA' twice — verify DFA handles it."""
        dfa = build_forbidden_pattern_dfa("GATAGATC")
        assert not dfa_accepts(dfa, "GATAGATC")
        assert dfa_accepts(dfa, "GATAGATA")  # missing last 'C'

    def test_ttt_run_dfa(self):
        """T-run DFA should reject sequences with 6+ consecutive T's."""
        dfa = build_trun_dfa(max_t=5)
        assert dfa_accepts(dfa, "TTTTT")       # 5 T's = OK
        assert not dfa_accepts(dfa, "TTTTTT")  # 6 T's = forbidden

    def test_ttt_run_with_interruption(self):
        """T-run should reset on non-T character."""
        dfa = build_trun_dfa(max_t=5)
        assert dfa_accepts(dfa, "TTTTTATTTTT")  # 5+5 T's separated by A


# ═══════════════════════════════════════════════════════════════════════════
# 6. Alphabet encoding
# ═══════════════════════════════════════════════════════════════════════════

class TestAlphabetEncoding:
    """Tests for BASE_MAP and BASE_REV encoding consistency."""

    def test_base_map_0_is_a(self):
        """0 should map to 'A'."""
        assert BASE_MAP[0] == "A"

    def test_base_map_1_is_c(self):
        """1 should map to 'C'."""
        assert BASE_MAP[1] == "C"

    def test_base_map_2_is_g(self):
        """2 should map to 'G'."""
        assert BASE_MAP[2] == "G"

    def test_base_map_3_is_t(self):
        """3 should map to 'T'."""
        assert BASE_MAP[3] == "T"

    def test_base_map_has_four_entries(self):
        """BASE_MAP should have exactly 4 entries."""
        assert len(BASE_MAP) == 4

    def test_base_map_keys_are_zero_to_three(self):
        """BASE_MAP keys should be exactly {0, 1, 2, 3}."""
        assert set(BASE_MAP.keys()) == {0, 1, 2, 3}

    def test_base_map_values_are_acgt(self):
        """BASE_MAP values should be exactly {'A', 'C', 'G', 'T'}."""
        assert set(BASE_MAP.values()) == {"A", "C", "G", "T"}

    def test_base_rev_a_is_0(self):
        """'A' should encode as 0."""
        assert BASE_REV["A"] == 0

    def test_base_rev_c_is_1(self):
        """'C' should encode as 1."""
        assert BASE_REV["C"] == 1

    def test_base_rev_g_is_2(self):
        """'G' should encode as 2."""
        assert BASE_REV["G"] == 2

    def test_base_rev_t_is_3(self):
        """'T' should encode as 3."""
        assert BASE_REV["T"] == 3

    def test_automaton_uses_correct_encoding(self):
        """The DFA transitions should use the BASE_REV encoding.

        Verifying that the DFA correctly encodes 'G' → 2 and 'A' → 0
        by checking a specific transition for the EcoRI pattern GAATTC.
        """
        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, accepting = dfa

        # From state 0, reading 'G' (index 2) should advance to state 1
        # because GAATTC starts with 'G'
        assert transitions[0][BASE_REV["G"]] == 1

        # From state 0, reading 'A' (index 0) should stay at state 0
        # because GAATTC doesn't start with 'A'
        assert transitions[0][BASE_REV["A"]] == 0

    def test_transition_from_state_1_reads_a(self):
        """From state 1 (matched 'G'), reading 'A' should advance to state 2."""
        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, _ = dfa
        # State 1 means we've matched 'G'; next char is 'A' (index 0)
        assert transitions[1][BASE_REV["A"]] == 2

    def test_transition_from_state_1_reads_wrong_char(self):
        """From state 1 (matched 'G'), reading 'C' should fall back."""
        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, _ = dfa
        # State 1 means we've matched 'G'; reading 'C' doesn't match 'A'
        next_state = transitions[1][BASE_REV["C"]]
        # Should NOT advance to state 2
        assert next_state != 2

    def test_reverse_complement_encoding(self):
        """Reverse complement should use correct base pairings."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"

    def test_reverse_complement_double(self):
        """Double reverse complement should be the identity."""
        seq = "GAATTC"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_reverse_complement_longer(self):
        """Reverse complement of a longer sequence."""
        assert reverse_complement("AACG") == "CGTT"
        assert reverse_complement("ATCG") == "CGAT"


# ═══════════════════════════════════════════════════════════════════════════
# 7. DFA structural properties
# ═══════════════════════════════════════════════════════════════════════════

class TestDFAStructure:
    """Tests verifying structural properties of the constructed DFAs."""

    def test_single_pattern_num_states(self):
        """A single-pattern DFA should have len(pattern)+1 states."""
        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, accepting = dfa
        # States: 0, 1, 2, 3, 4, 5, 6 (7 states for a 6-char pattern)
        # Number of states = number of rows in transition table
        assert len(transitions) == 7

    def test_every_state_has_four_transitions(self):
        """Every state in the DFA should have transitions for all 4 bases."""
        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, _ = dfa
        for state in range(len(transitions)):
            assert len(transitions[state]) == 4

    def test_forbidden_state_is_trap(self):
        """The pattern-matching DFA's forbidden state is a trap.

        Once the forbidden state is entered, the DFA never leaves.
        dfa_accepts() returns False because the final state is not accepting.
        """
        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, accepting = dfa
        # State 6 is the trap (forbidden) state
        assert 6 not in accepting
        # All transitions from state 6 loop back to itself
        for c in range(4):
            assert transitions[6][c] == 6
        assert not dfa_accepts(dfa, "GAATTC")

    def test_composite_dfa_has_transitions(self, composite_dfa):
        """Composite DFA should have a non-empty transition table."""
        transitions, _ = composite_dfa
        assert len(transitions) > 0

    def test_composite_dfa_has_forbidden_states(self, composite_dfa):
        """Composite DFA should have at least one non-accepting (trap) state."""
        transitions, accepting = composite_dfa
        all_states = set(range(len(transitions)))
        forbidden = all_states - set(accepting)
        assert len(forbidden) > 0

    def test_dfa_is_deterministic(self, ecori_dfa):
        """Each state should have exactly one transition per symbol."""
        transitions, _ = ecori_dfa
        for state in range(len(transitions)):
            assert len(transitions[state]) == 4
            for base in range(4):
                assert isinstance(transitions[state][base], int)

    def test_negate_dfa_swaps_acceptance(self):
        """negate_dfa should swap accepting and non-accepting states."""
        transitions, accepting = build_forbidden_pattern_dfa("AA")
        # Original: accepts strings WITHOUT "AA"
        assert dfa_accepts((transitions, accepting), "AC")
        assert not dfa_accepts((transitions, accepting), "AAC")
        # Negated: accepts strings WITH "AA"
        n_trans, n_acc = negate_dfa(transitions, accepting)
        assert dfa_accepts((n_trans, n_acc), "AAC")


# ═══════════════════════════════════════════════════════════════════════════
# 8. Integration / real-world scenarios
# ═══════════════════════════════════════════════════════════════════════════

class TestRealWorldScenarios:
    """Tests mimicking real codon-optimization use cases."""

    def test_codon_sequence_without_restriction_sites(self, composite_dfa):
        """A typical codon-optimized sequence should be accepted."""
        # ATG GCT GCT GCT = M A A A (using GC-rich codons)
        seq = "ATGGCTGCTGCT"
        assert dfa_accepts(composite_dfa, seq)

    def test_sequence_with_ecori_site_detected(self, composite_dfa):
        """A sequence that happens to contain an EcoRI site should be rejected."""
        # ATG GAA TTC GCT = M E F A — but GAA|TTC contains GAATTC
        seq = "ATGGAATTCGCT"
        assert not dfa_accepts(composite_dfa, seq)

    def test_multiple_enzymes_simultaneously(self):
        """DFA for many restriction enzymes should reject all their sites."""
        sites = ["GAATTC", "GGATCC", "CTCGAG", "AAGCTT", "GTCGAC"]
        # Include reverse complements explicitly
        all_sites = list(set(sites + [reverse_complement(s) for s in sites]))
        dfa = build_composite_dfa(all_sites)
        for site in sites:
            assert not dfa_accepts(dfa, site), f"Should reject {site}"
        # Clean sequence
        assert dfa_accepts(dfa, "ATGGCTGCTGCT")

    def test_very_long_sequence(self, ecori_dfa):
        """DFA should handle very long sequences efficiently."""
        # 10000 bases, no GAATTC
        seq = "ATCG" * 2500
        assert dfa_accepts(ecori_dfa, seq)

    def test_very_long_sequence_with_pattern(self, ecori_dfa):
        """DFA should find pattern in a very long sequence."""
        seq = "ATCG" * 2500 + "GAATTC" + "ATCG" * 100
        assert not dfa_accepts(ecori_dfa, seq)

    def test_or_tools_format_conversion(self):
        """DFA should be convertible to OR-Tools format."""
        from biocompiler.solver.automaton import dfa_to_ortools_format

        dfa = build_forbidden_pattern_dfa("GAATTC")
        transitions, accepting = dfa
        init, trans_list, final = dfa_to_ortools_format(transitions, accepting)

        assert init == 0
        assert len(trans_list) > 0
        assert len(final) > 0
        # Each transition is [from_state, symbol, to_state]
        for t in trans_list:
            assert len(t) == 3
            assert isinstance(t[0], int)
            assert isinstance(t[1], int)
            assert isinstance(t[2], int)
