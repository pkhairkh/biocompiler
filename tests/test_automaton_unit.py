"""Unit tests for biocompiler.solver.automaton.

Covers:
1. DFA construction for forbidden patterns  (build_forbidden_pattern_dfa)
2. DFA acceptance/rejection                  (dfa_accepts)
3. Reverse complement DFA                    (build_reverse_complement_dfa)
4. T-run DFA                                 (build_trun_dfa)
5. DFA to OR-Tools format conversion         (dfa_to_ortools_format)

Also exercises: build_composite_dfa, negate_dfa, dfa_to_dot, edge cases.
"""

from __future__ import annotations

import pytest

from biocompiler.solver.automaton import (
    build_composite_dfa,
    build_forbidden_pattern_dfa,
    build_reverse_complement_dfa,
    build_trun_dfa,
    dfa_accepts,
    dfa_to_dot,
    dfa_to_ortools_format,
    negate_dfa,
)
from biocompiler.constants import BASE_MAP, BASE_REV, reverse_complement


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def ecori_dfa():
    """DFA that forbids the EcoRI site GAATTC."""
    return build_forbidden_pattern_dfa("GAATTC")


@pytest.fixture
def bamhi_dfa():
    """DFA that forbids the BamHI site GGATCC."""
    return build_forbidden_pattern_dfa("GGATCC")


@pytest.fixture
def single_base_dfa():
    """DFA that forbids the single character 'G'."""
    return build_forbidden_pattern_dfa("G")


@pytest.fixture
def aaa_dfa():
    """DFA that forbids the triple-A pattern 'AAA'."""
    return build_forbidden_pattern_dfa("AAA")


@pytest.fixture
def composite_ecori_bamhi_dfa():
    """Composite DFA forbidding both GAATTC and GGATCC."""
    return build_composite_dfa(["GAATTC", "GGATCC"])


@pytest.fixture
def trun5_dfa():
    """T-run DFA that forbids 6+ consecutive T's (max_t=5)."""
    return build_trun_dfa(max_t=5)


# ═══════════════════════════════════════════════════════════════════════════
# 1. DFA construction for forbidden patterns
# ═══════════════════════════════════════════════════════════════════════════

class TestBuildForbiddenPatternDFA:
    """Tests for build_forbidden_pattern_dfa — single-pattern DFA construction."""

    # --- Basic structure ---------------------------------------------------

    def test_returns_transition_table_and_accepting_states(self, ecori_dfa):
        """DFA should return a (transition_table, accepting_states) tuple."""
        delta, accepting = ecori_dfa
        assert isinstance(delta, list)
        assert isinstance(accepting, list)
        # Transition table is a list of lists
        for row in delta:
            assert isinstance(row, list)
            assert len(row) == 4  # one entry per base (A, C, G, T)

    def test_num_states_equals_pattern_length_plus_one(self, ecori_dfa):
        """For pattern GAATTC (length 6), there should be 7 states (0..6)."""
        delta, _ = ecori_dfa
        assert len(delta) == 7

    def test_accepting_states_exclude_trap(self, ecori_dfa):
        """The trap state (state len(pattern)) should NOT be accepting."""
        delta, accepting = ecori_dfa
        trap = len(delta) - 1  # state 6 for GAATTC
        assert trap not in accepting
        # All other states should be accepting
        for s in range(trap):
            assert s in accepting

    def test_trap_state_is_self_loop(self, ecori_dfa):
        """The trap state should transition to itself on every input."""
        delta, _ = ecori_dfa
        trap = len(delta) - 1
        for c in range(4):
            assert delta[trap][c] == trap

    def test_single_char_pattern_has_two_states(self, single_base_dfa):
        """A single-char pattern should produce 2 states: start and trap."""
        delta, _ = single_base_dfa
        assert len(delta) == 2

    def test_single_char_trap_not_accepting(self, single_base_dfa):
        """For pattern 'G', state 1 (trap) should not be accepting."""
        _, accepting = single_base_dfa
        assert 1 not in accepting
        assert 0 in accepting

    def test_empty_pattern_accepts_everything(self):
        """An empty pattern means nothing is forbidden — single state, always accepting."""
        delta, accepting = build_forbidden_pattern_dfa("")
        assert len(delta) == 1
        assert accepting == [0]
        # Every symbol loops back to state 0
        for c in range(4):
            assert delta[0][c] == 0

    # --- Transition correctness for EcoRI (GAATTC) ------------------------

    def test_ecori_state0_on_G_goes_to_state1(self, ecori_dfa):
        """From state 0, reading 'G' (first char of GAATTC) should go to state 1."""
        delta, _ = ecori_dfa
        assert delta[0][BASE_REV["G"]] == 1

    def test_ecori_state0_on_A_stays_at_0(self, ecori_dfa):
        """From state 0, reading 'A' (not first char) should stay at 0."""
        delta, _ = ecori_dfa
        assert delta[0][BASE_REV["A"]] == 0

    def test_ecori_state1_on_A_goes_to_state2(self, ecori_dfa):
        """State 1 (matched 'G'), reading 'A' → state 2."""
        delta, _ = ecori_dfa
        assert delta[1][BASE_REV["A"]] == 2

    def test_ecori_full_match_goes_to_trap(self, ecori_dfa):
        """Matching the full GAATTC should land in the trap state."""
        delta, _ = ecori_dfa
        state = 0
        for ch in "GAATTC":
            state = delta[state][BASE_REV[ch]]
        assert state == 6  # trap

    def test_ecori_mismatch_after_prefix_falls_back(self, ecori_dfa):
        """After partial match GAAT, a 'G' should fall back, not go to trap."""
        delta, _ = ecori_dfa
        state = 0
        for ch in "GAATG":  # GAAT then G instead of T
            state = delta[state][BASE_REV[ch]]
        # Should NOT be in trap state 6
        assert state != 6

    # --- KMP failure link: self-overlapping pattern ------------------------

    def test_overlapping_pattern_aaa_state_transitions(self, aaa_dfa):
        """For pattern AAA, verify KMP failure links work correctly."""
        delta, _ = aaa_dfa
        # State 0 + A → 1
        assert delta[0][BASE_REV["A"]] == 1
        # State 1 + A → 2
        assert delta[1][BASE_REV["A"]] == 2
        # State 2 + A → 3 (trap)
        assert delta[2][BASE_REV["A"]] == 3
        # State 2 + non-A should fall back based on KMP
        # After "AA" + C, the longest prefix of "AAA" that matches suffix of "AAC" is 0
        assert delta[2][BASE_REV["C"]] == 0

    # --- Invalid input -----------------------------------------------------

    def test_invalid_base_raises_valueerror(self):
        """Passing an invalid character in the pattern should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid base"):
            build_forbidden_pattern_dfa("GAXTC")

    def test_lowercase_pattern_is_accepted(self):
        """Lowercase pattern should be uppercased internally and work correctly."""
        dfa = build_forbidden_pattern_dfa("gaattc")
        assert not dfa_accepts(dfa, "GAATTC")
        assert dfa_accepts(dfa, "GAATTG")

    # --- Completeness: every state has 4 transitions ----------------------

    def test_every_state_has_four_transitions(self, ecori_dfa):
        """Every state should have a defined transition for each of the 4 bases."""
        delta, _ = ecori_dfa
        for state_idx, row in enumerate(delta):
            assert len(row) == 4, f"State {state_idx} has {len(row)} transitions, expected 4"
            for c in range(4):
                assert isinstance(row[c], int), f"State {state_idx}, symbol {c} is not int"
                assert 0 <= row[c] < len(delta), (
                    f"State {state_idx}, symbol {c} → {row[c]} out of range [0, {len(delta)})"
                )


# ═══════════════════════════════════════════════════════════════════════════
# 2. DFA acceptance / rejection
# ═══════════════════════════════════════════════════════════════════════════

class TestDfaAccepts:
    """Tests for the dfa_accepts function."""

    # --- Exact pattern -----------------------------------------------------

    def test_rejects_exact_forbidden_pattern(self, ecori_dfa):
        """DFA should reject the exact forbidden pattern."""
        assert dfa_accepts(ecori_dfa, "GAATTC") is False

    def test_accepts_sequence_without_pattern(self, ecori_dfa):
        """DFA should accept a clean sequence."""
        assert dfa_accepts(ecori_dfa, "ATCGATCG") is True

    # --- Pattern embedded in context ---------------------------------------

    def test_rejects_pattern_at_start(self, ecori_dfa):
        """Pattern at the start of a sequence should be rejected."""
        assert dfa_accepts(ecori_dfa, "GAATTCATCG") is False

    def test_rejects_pattern_at_middle(self, ecori_dfa):
        """Pattern in the middle should be rejected."""
        assert dfa_accepts(ecori_dfa, "ATCGAATTCGAT") is False

    def test_rejects_pattern_at_end(self, ecori_dfa):
        """Pattern at the end should be rejected."""
        assert dfa_accepts(ecori_dfa, "ATCGATGAATTC") is False

    # --- Near misses -------------------------------------------------------

    def test_accepts_near_miss_last_char(self, ecori_dfa):
        """GAATTG differs in last char — should be accepted."""
        assert dfa_accepts(ecori_dfa, "GAATTG") is True

    def test_accepts_near_miss_middle_char(self, ecori_dfa):
        """GAATAC differs in middle — should be accepted."""
        assert dfa_accepts(ecori_dfa, "GAATAC") is True

    def test_accepts_partial_prefix(self, ecori_dfa):
        """GAATT is a prefix of GAATTC but not the full pattern — accepted."""
        assert dfa_accepts(ecori_dfa, "GAATT") is True

    # --- Short / empty sequences -------------------------------------------

    def test_accepts_empty_sequence(self, ecori_dfa):
        """An empty sequence cannot contain any forbidden pattern."""
        assert dfa_accepts(ecori_dfa, "") is True

    def test_accepts_sequence_shorter_than_pattern(self, ecori_dfa):
        """Sequences shorter than the pattern length are always accepted."""
        assert dfa_accepts(ecori_dfa, "G") is True
        assert dfa_accepts(ecori_dfa, "GA") is True
        assert dfa_accepts(ecori_dfa, "GAA") is True
        assert dfa_accepts(ecori_dfa, "GAAT") is True
        assert dfa_accepts(ecori_dfa, "GAATT") is True

    def test_single_char_all_bases_accepted(self, ecori_dfa):
        """Any single base is accepted (pattern has length 6)."""
        for base in "ACGT":
            assert dfa_accepts(ecori_dfa, base) is True

    # --- Single-char forbidden pattern -------------------------------------

    def test_single_char_rejects_matching_char(self, single_base_dfa):
        """For pattern 'G', any string containing G is rejected."""
        assert dfa_accepts(single_base_dfa, "G") is False
        assert dfa_accepts(single_base_dfa, "ATCG") is False
        assert dfa_accepts(single_base_dfa, "GATC") is False

    def test_single_char_accepts_without_char(self, single_base_dfa):
        """Strings without G should be accepted."""
        assert dfa_accepts(single_base_dfa, "A") is True
        assert dfa_accepts(single_base_dfa, "ATC") is True
        assert dfa_accepts(single_base_dfa, "TTTT") is True

    # --- Multiple occurrences / overlapping --------------------------------

    def test_rejects_multiple_occurrences(self, ecori_dfa):
        """Two occurrences should still be rejected."""
        assert dfa_accepts(ecori_dfa, "GAATTCGAATTC") is False

    def test_overlapping_aaa_in_aaaa(self, aaa_dfa):
        """AAAA contains AAA starting at positions 0 and 1 — rejected."""
        assert dfa_accepts(aaa_dfa, "AAAA") is False

    def test_overlapping_aaa_accepts_aa(self, aaa_dfa):
        """AA is too short to contain AAA — accepted."""
        assert dfa_accepts(aaa_dfa, "AA") is True

    # --- Case handling -----------------------------------------------------

    def test_lowercase_sequence(self, ecori_dfa):
        """dfa_accepts should handle lowercase sequences (uppercased internally)."""
        assert dfa_accepts(ecori_dfa, "gaattc") is False
        assert dfa_accepts(ecori_dfa, "atcgatcg") is True

    # --- Return type -------------------------------------------------------

    def test_returns_bool(self, ecori_dfa):
        """dfa_accepts should return a bool, not a truthy/falsy value."""
        result = dfa_accepts(ecori_dfa, "ATCG")
        assert isinstance(result, bool)

    # --- Very long sequences -----------------------------------------------

    def test_long_clean_sequence(self, ecori_dfa):
        """DFA should handle long clean sequences efficiently."""
        seq = "ATCG" * 2500
        assert dfa_accepts(ecori_dfa, seq) is True

    def test_long_sequence_with_pattern_at_end(self, ecori_dfa):
        """DFA should detect pattern even in a very long sequence."""
        seq = "ATCG" * 2500 + "GAATTC"
        assert dfa_accepts(ecori_dfa, seq) is False


# ═══════════════════════════════════════════════════════════════════════════
# 3. Reverse complement DFA
# ═══════════════════════════════════════════════════════════════════════════

class TestReverseComplementDFA:
    """Tests for build_reverse_complement_dfa."""

    # --- Palindromic patterns (pattern == revcomp) -------------------------

    def test_palindromic_ecori_uses_single_dfa(self):
        """For palindromic GAATTC, RC DFA should be equivalent to single-pattern DFA."""
        rc_dfa = build_reverse_complement_dfa("GAATTC")
        single_dfa = build_forbidden_pattern_dfa("GAATTC")
        # Same number of states
        assert len(rc_dfa[0]) == len(single_dfa[0])

    def test_palindromic_rejects_original(self):
        """Palindromic RC DFA rejects the original pattern."""
        rc_dfa = build_reverse_complement_dfa("GAATTC")
        assert dfa_accepts(rc_dfa, "GAATTC") is False

    def test_palindromic_rejects_rc(self):
        """Palindromic RC DFA rejects the reverse complement (same as original)."""
        rc_dfa = build_reverse_complement_dfa("GAATTC")
        rc = reverse_complement("GAATTC")
        assert rc == "GAATTC"  # confirm it's palindromic
        assert dfa_accepts(rc_dfa, rc) is False

    def test_palindromic_hindiii(self):
        """HindIII (AAGCTT) is palindromic — RC DFA rejects it."""
        assert reverse_complement("AAGCTT") == "AAGCTT"
        rc_dfa = build_reverse_complement_dfa("AAGCTT")
        assert dfa_accepts(rc_dfa, "AAGCTT") is False

    # --- Non-palindromic patterns (pattern != revcomp) ---------------------

    def test_non_palindromic_rc_differs(self):
        """For non-palindromic pattern, RC != original."""
        assert reverse_complement("AACG") != "AACG"
        assert reverse_complement("AACG") == "CGTT"

    def test_non_palindromic_rejects_original(self):
        """RC DFA should reject the original pattern."""
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert dfa_accepts(rc_dfa, "AACG") is False

    def test_non_palindromic_rejects_rc(self):
        """RC DFA should reject the reverse complement."""
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert dfa_accepts(rc_dfa, "CGTT") is False

    def test_non_palindromic_accepts_neither(self):
        """RC DFA should accept sequences containing neither pattern nor RC."""
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert dfa_accepts(rc_dfa, "TTTTTT") is True

    def test_non_palindromic_rejects_original_embedded(self):
        """RC DFA rejects sequences with the original pattern embedded."""
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert dfa_accepts(rc_dfa, "ATAACGAT") is False

    def test_non_palindromic_rejects_rc_embedded(self):
        """RC DFA rejects sequences with the RC embedded."""
        rc_dfa = build_reverse_complement_dfa("AACG")
        assert dfa_accepts(rc_dfa, "ATCGTTAT") is False

    # --- Composite DFA with RC vs build_reverse_complement_dfa ------------

    def test_composite_with_rc_rejects_both(self):
        """build_composite_dfa with both pattern and RC rejects both."""
        dfa = build_composite_dfa(["AACG", "CGTT"])
        assert dfa_accepts(dfa, "AACG") is False
        assert dfa_accepts(dfa, "CGTT") is False

    # --- Edge case: single-char palindromic --------------------------------

    def test_single_char_not_palindromic(self):
        """A single base like 'A' has RC 'T' — not palindromic."""
        assert reverse_complement("A") == "T"
        rc_dfa = build_reverse_complement_dfa("A")
        assert dfa_accepts(rc_dfa, "A") is False
        assert dfa_accepts(rc_dfa, "T") is False
        assert dfa_accepts(rc_dfa, "C") is True


# ═══════════════════════════════════════════════════════════════════════════
# 4. T-run DFA
# ═══════════════════════════════════════════════════════════════════════════

class TestTrunDFA:
    """Tests for build_trun_dfa — poly-T run length constraint."""

    # --- Basic acceptance / rejection --------------------------------------

    def test_accepts_max_t_consecutive_ts(self, trun5_dfa):
        """5 consecutive T's should be accepted (max_t=5)."""
        assert dfa_accepts(trun5_dfa, "TTTTT") is True

    def test_rejects_max_t_plus_one_ts(self, trun5_dfa):
        """6 consecutive T's should be rejected (max_t=5 → 6 forbidden)."""
        assert dfa_accepts(trun5_dfa, "TTTTTT") is False

    def test_accepts_shorter_t_run(self, trun5_dfa):
        """3 consecutive T's should be accepted."""
        assert dfa_accepts(trun5_dfa, "TTT") is True

    def test_accepts_single_t(self, trun5_dfa):
        """A single T should be accepted."""
        assert dfa_accepts(trun5_dfa, "T") is True

    # --- Run interruption --------------------------------------------------

    def test_accepts_broken_t_run(self, trun5_dfa):
        """Two T-runs of length 5 separated by a non-T base — accepted."""
        assert dfa_accepts(trun5_dfa, "TTTTTATTTTT") is True

    def test_rejects_long_run_after_interruption(self, trun5_dfa):
        """A too-long run after interruption should still be rejected."""
        assert dfa_accepts(trun5_dfa, "TTTTTATTTTTT") is False

    def test_non_t_base_resets_counter(self, trun5_dfa):
        """Any non-T base should reset the T counter to 0."""
        assert dfa_accepts(trun5_dfa, "TTTTTCTTTTT") is True
        assert dfa_accepts(trun5_dfa, "TTTTTGTTTTT") is True

    # --- All-non-T sequences -----------------------------------------------

    def test_accepts_all_a_sequence(self, trun5_dfa):
        """An all-A sequence has no T runs — accepted."""
        assert dfa_accepts(trun5_dfa, "AAAAAAAAAA") is True

    def test_accepts_no_t_at_all(self, trun5_dfa):
        """A sequence with no T's at all — accepted."""
        assert dfa_accepts(trun5_dfa, "ACGACGACG") is True

    # --- Empty / short sequences -------------------------------------------

    def test_accepts_empty_sequence(self, trun5_dfa):
        """Empty sequence is accepted."""
        assert dfa_accepts(trun5_dfa, "") is True

    # --- Structure ---------------------------------------------------------

    def test_num_states_is_max_t_plus_two(self, trun5_dfa):
        """For max_t=5, there should be 7 states (0..5 + trap=6)."""
        delta, _ = trun5_dfa
        assert len(delta) == 7  # states 0,1,2,3,4,5,6

    def test_trap_state_is_last(self, trun5_dfa):
        """The trap state should be state max_t+1 = 6."""
        delta, _ = trun5_dfa
        trap = 6
        # Trap is self-loop
        for c in range(4):
            assert delta[trap][c] == trap

    def test_trap_not_in_accepting(self, trun5_dfa):
        """The trap state should not be in the accepting list."""
        _, accepting = trun5_dfa
        assert 6 not in accepting

    def test_accepting_states_are_0_through_max_t(self, trun5_dfa):
        """Accepting states should be [0, 1, 2, 3, 4, 5]."""
        _, accepting = trun5_dfa
        assert accepting == [0, 1, 2, 3, 4, 5]

    def test_non_t_transition_goes_to_state0(self, trun5_dfa):
        """From any non-trap state, a non-T base should go to state 0."""
        delta, _ = trun5_dfa
        t_idx = BASE_REV["T"]
        for state in range(6):  # states 0..5
            for c in range(4):
                if c != t_idx:
                    assert delta[state][c] == 0, (
                        f"State {state}, base {c} (non-T) should go to 0, got {delta[state][c]}"
                    )

    def test_t_transition_increments(self, trun5_dfa):
        """From state s (< max_t), T should go to state s+1."""
        delta, _ = trun5_dfa
        t_idx = BASE_REV["T"]
        for state in range(5):  # states 0..4
            assert delta[state][t_idx] == state + 1

    def test_t_at_max_t_goes_to_trap(self, trun5_dfa):
        """From state 5 (max_t), T should go to trap state 6."""
        delta, _ = trun5_dfa
        t_idx = BASE_REV["T"]
        assert delta[5][t_idx] == 6

    # --- Different max_t values -------------------------------------------

    def test_max_t_3_rejects_4_ts(self):
        """With max_t=3, 4 consecutive T's are forbidden."""
        dfa = build_trun_dfa(max_t=3)
        assert dfa_accepts(dfa, "TTT") is True
        assert dfa_accepts(dfa, "TTTT") is False

    def test_max_t_1_rejects_2_ts(self):
        """With max_t=1, 2 consecutive T's are forbidden."""
        dfa = build_trun_dfa(max_t=1)
        assert dfa_accepts(dfa, "T") is True
        assert dfa_accepts(dfa, "TT") is False
        assert dfa_accepts(dfa, "AT") is True

    def test_max_t_0_rejects_any_t(self):
        """With max_t=0, any T at all is forbidden."""
        dfa = build_trun_dfa(max_t=0)
        assert dfa_accepts(dfa, "T") is False
        assert dfa_accepts(dfa, "A") is True
        assert dfa_accepts(dfa, "ATCG") is False

    def test_default_max_t_is_5(self):
        """Default max_t should be 5, so 6 T's are forbidden."""
        dfa = build_trun_dfa()
        assert dfa_accepts(dfa, "TTTTT") is True
        assert dfa_accepts(dfa, "TTTTTT") is False


# ═══════════════════════════════════════════════════════════════════════════
# 5. DFA to OR-Tools format conversion
# ═══════════════════════════════════════════════════════════════════════════

class TestDfaToOrToolsFormat:
    """Tests for dfa_to_ortools_format conversion."""

    def test_returns_triple(self, ecori_dfa):
        """Should return (initial_state, transition_list, final_states)."""
        delta, accepting = ecori_dfa
        result = dfa_to_ortools_format(delta, accepting)
        assert len(result) == 3

    def test_initial_state_is_zero(self, ecori_dfa):
        """Initial state should always be 0."""
        delta, accepting = ecori_dfa
        init, _, _ = dfa_to_ortools_format(delta, accepting)
        assert init == 0

    def test_transition_list_length(self, ecori_dfa):
        """Number of transitions = num_states × alphabet_size."""
        delta, accepting = ecori_dfa
        _, trans_list, _ = dfa_to_ortools_format(delta, accepting)
        expected = len(delta) * 4  # 7 states × 4 symbols
        assert len(trans_list) == expected

    def test_transition_list_entry_format(self, ecori_dfa):
        """Each transition entry should be [from_state, symbol, to_state]."""
        delta, accepting = ecori_dfa
        _, trans_list, _ = dfa_to_ortools_format(delta, accepting)
        for entry in trans_list:
            assert len(entry) == 3
            assert isinstance(entry[0], int)  # from_state
            assert isinstance(entry[1], int)  # symbol (0-3)
            assert isinstance(entry[2], int)  # to_state

    def test_symbol_range(self, ecori_dfa):
        """Symbol values should be in {0, 1, 2, 3}."""
        delta, accepting = ecori_dfa
        _, trans_list, _ = dfa_to_ortools_format(delta, accepting)
        for entry in trans_list:
            assert 0 <= entry[1] <= 3

    def test_final_states_are_sorted(self, ecori_dfa):
        """Final states should be returned sorted."""
        delta, accepting = ecori_dfa
        _, _, final = dfa_to_ortools_format(delta, accepting)
        assert final == sorted(final)

    def test_final_states_match_accepting(self, ecori_dfa):
        """Final states should match the accepting states (sorted)."""
        delta, accepting = ecori_dfa
        _, _, final = dfa_to_ortools_format(delta, accepting)
        assert final == sorted(accepting)

    def test_transition_consistency_with_original_table(self, ecori_dfa):
        """Each transition in the OR-Tools list should match the original table."""
        delta, accepting = ecori_dfa
        _, trans_list, _ = dfa_to_ortools_format(delta, accepting)
        for from_state, symbol, to_state in trans_list:
            assert delta[from_state][symbol] == to_state

    def test_single_char_pattern_conversion(self, single_base_dfa):
        """Single-char pattern should produce 2×4=8 transitions."""
        delta, accepting = single_base_dfa
        _, trans_list, _ = dfa_to_ortools_format(delta, accepting)
        assert len(trans_list) == 8  # 2 states × 4 symbols

    def test_trun_dfa_conversion(self, trun5_dfa):
        """T-run DFA with max_t=5 should produce 7×4=28 transitions."""
        delta, accepting = trun5_dfa
        _, trans_list, _ = dfa_to_ortools_format(delta, accepting)
        assert len(trans_list) == 28  # 7 states × 4 symbols

    def test_empty_pattern_conversion(self):
        """Empty pattern DFA (1 state) should produce 4 transitions."""
        delta, accepting = build_forbidden_pattern_dfa("")
        init, trans_list, final = dfa_to_ortools_format(delta, accepting)
        assert init == 0
        assert len(trans_list) == 4  # 1 state × 4 symbols
        assert final == [0]

    def test_composite_dfa_conversion(self, composite_ecori_bamhi_dfa):
        """Composite DFA conversion should produce valid OR-Tools format."""
        delta, accepting = composite_ecori_bamhi_dfa
        init, trans_list, final = dfa_to_ortools_format(delta, accepting)
        assert init == 0
        assert len(trans_list) == len(delta) * 4
        assert len(final) > 0
        # Verify all accepting states are present
        for s in accepting:
            assert s in final


# ═══════════════════════════════════════════════════════════════════════════
# Additional: Composite DFA, negate_dfa, dfa_to_dot
# ═══════════════════════════════════════════════════════════════════════════

class TestCompositeDFA:
    """Tests for build_composite_dfa — Aho-Corasick multi-pattern DFA."""

    def test_rejects_first_pattern(self, composite_ecori_bamhi_dfa):
        """Composite DFA rejects sequences containing GAATTC."""
        assert dfa_accepts(composite_ecori_bamhi_dfa, "GAATTC") is False

    def test_rejects_second_pattern(self, composite_ecori_bamhi_dfa):
        """Composite DFA rejects sequences containing GGATCC."""
        assert dfa_accepts(composite_ecori_bamhi_dfa, "GGATCC") is False

    def test_accepts_clean_sequence(self, composite_ecori_bamhi_dfa):
        """Composite DFA accepts sequences with neither pattern."""
        assert dfa_accepts(composite_ecori_bamhi_dfa, "ATCGATCG") is True

    def test_rejects_sequence_with_both_patterns(self, composite_ecori_bamhi_dfa):
        """Sequence containing both forbidden patterns is rejected."""
        assert dfa_accepts(composite_ecori_bamhi_dfa, "GAATTCGGATCC") is False

    def test_empty_patterns_accepts_everything(self):
        """Empty pattern list should produce a DFA that accepts everything."""
        dfa = build_composite_dfa([])
        assert dfa_accepts(dfa, "GAATTC") is True
        assert dfa_accepts(dfa, "GGATCC") is True

    def test_single_pattern_composite_equivalent_to_single_dfa(self):
        """Composite DFA with one pattern should behave like single-pattern DFA."""
        composite = build_composite_dfa(["GAATTC"])
        single = build_forbidden_pattern_dfa("GAATTC")
        test_seqs = ["GAATTC", "GAATTG", "ATCG", "", "GAATTCGAATTC"]
        for seq in test_seqs:
            assert dfa_accepts(composite, seq) == dfa_accepts(single, seq), (
                f"Disagreement on '{seq}'"
            )

    def test_overlapping_patterns_in_composite(self):
        """Patterns that share prefixes (e.g. GAT and GATC) should both be rejected."""
        dfa = build_composite_dfa(["GAT", "GATC"])
        assert dfa_accepts(dfa, "GAT") is False
        assert dfa_accepts(dfa, "GATC") is False
        assert dfa_accepts(dfa, "GA") is True

    def test_one_pattern_is_prefix_of_another(self):
        """If one pattern is a prefix of another, both should be detected."""
        dfa = build_composite_dfa(["AA", "AAA"])
        assert dfa_accepts(dfa, "AA") is False
        assert dfa_accepts(dfa, "AAA") is False
        assert dfa_accepts(dfa, "A") is True


class TestNegateDFA:
    """Tests for negate_dfa — swapping accepting and non-accepting states."""

    def test_negated_dfa_accepts_forbidden_sequence(self, ecori_dfa):
        """Negating should make the DFA accept the previously forbidden pattern."""
        delta, accepting = ecori_dfa
        n_delta, n_accepting = negate_dfa(delta, accepting)
        assert dfa_accepts((n_delta, n_accepting), "GAATTC") is True

    def test_negated_dfa_rejects_clean_sequence(self, ecori_dfa):
        """Negating should make the DFA reject previously accepted sequences."""
        delta, accepting = ecori_dfa
        n_delta, n_accepting = negate_dfa(delta, accepting)
        # A clean sequence like "ATCGATCG" was accepted; now it's rejected
        assert dfa_accepts((n_delta, n_accepting), "ATCGATCG") is False

    def test_double_negation_restores_original(self, ecori_dfa):
        """Negating twice should restore the original accepting set."""
        delta, accepting = ecori_dfa
        _, n_acc = negate_dfa(delta, accepting)
        _, nn_acc = negate_dfa(delta, n_acc)
        assert set(nn_acc) == set(accepting)

    def test_negation_preserves_transition_table(self, ecori_dfa):
        """The transition table should be the same object after negation."""
        delta, accepting = ecori_dfa
        n_delta, _ = negate_dfa(delta, accepting)
        assert n_delta is delta  # same object

    def test_negated_accepting_complement(self, ecori_dfa):
        """Negated accepting states should be the complement of the original."""
        delta, accepting = ecori_dfa
        _, n_accepting = negate_dfa(delta, accepting)
        all_states = set(range(len(delta)))
        assert set(n_accepting) == all_states - set(accepting)

    def test_negate_trun_dfa(self, trun5_dfa):
        """Negating T-run DFA should accept sequences with 6+ T's."""
        delta, accepting = trun5_dfa
        n_delta, n_accepting = negate_dfa(delta, accepting)
        assert dfa_accepts((n_delta, n_accepting), "TTTTTT") is True
        assert dfa_accepts((n_delta, n_accepting), "TTTTT") is False


class TestDfaToDot:
    """Tests for dfa_to_dot — Graphviz DOT export."""

    def test_produces_dot_string(self, ecori_dfa):
        """Should produce a non-empty DOT string."""
        dot = dfa_to_dot(ecori_dfa)
        assert isinstance(dot, str)
        assert len(dot) > 0

    def test_contains_digraph_header(self, ecori_dfa):
        """DOT string should contain the digraph header."""
        dot = dfa_to_dot(ecori_dfa, name="ecori")
        assert "digraph ecori" in dot

    def test_contains_rankdir(self, ecori_dfa):
        """DOT string should contain rankdir=LR."""
        dot = dfa_to_dot(ecori_dfa)
        assert "rankdir=LR" in dot

    def test_custom_name(self, ecori_dfa):
        """Custom name should appear in the digraph header."""
        dot = dfa_to_dot(ecori_dfa, name="my_dfa")
        assert "digraph my_dfa" in dot

    def test_accepting_states_double_circle(self, ecori_dfa):
        """Accepting states should be drawn as doublecircle."""
        dot = dfa_to_dot(ecori_dfa)
        assert "doublecircle" in dot

    def test_start_arrow(self, ecori_dfa):
        """DOT should have a start arrow pointing to state 0."""
        dot = dfa_to_dot(ecori_dfa)
        assert '"" -> 0' in dot or '"" ->0' in dot or 'start' in dot

    def test_single_state_dfa_dot(self):
        """Empty-pattern DFA (1 state) should produce valid DOT."""
        dfa = build_forbidden_pattern_dfa("")
        dot = dfa_to_dot(dfa, name="trivial")
        assert "digraph trivial" in dot
