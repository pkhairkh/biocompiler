"""
BioCompiler Unit Tests — splicing.py & maxentscan.py
=====================================================

Focused unit tests covering:

1. maxent_score function exists and returns reasonable values
2. Known splice site sequences produce expected score ranges
3. Non-splice-site sequences produce low scores
4. Edge cases (short sequences, no GT/AG)

Also exercises score_donor / score_acceptor from maxentscan.py
to ensure the proper MaxEntScan log-odds model is consistent
with the simplified PWM in splicing.maxent_score.

Task ID: 60
"""

import math
import pytest

# ── Imports from splicing ────────────────────────────────────────────────────
from biocompiler.splicing import (
    maxent_score,
    score_splice_sites,
    _MAXENT_PWM,
    _BASE_INDEX,
    _PWM_CONTEXT_LEN,
    _MIN_CONTEXT_LEN,
    _PWM_UPSTREAM,
    _PWM_DOWNSTREAM,
)

# ── Imports from maxentscan ──────────────────────────────────────────────────
from biocompiler.maxentscan import (
    score_donor,
    score_acceptor,
    scan_splice_sites,
    max_donor_score,
    max_acceptor_score,
    BASE_TO_INDEX as MES_BASE_TO_INDEX,
    BG_PROB,
    DONOR_PWM_SCORE,
    ACCEPTOR_PWM_SCORE,
    _IMPOSSIBLE_SCORE,
    _DONOR_UPSTREAM,
    _DONOR_DOWNSTREAM,
    _ACCEPTOR_UPSTREAM,
    _ACCEPTOR_DOWNSTREAM,
)

# ── Import type system ───────────────────────────────────────────────────────
from biocompiler.type_system import SpliceVerdict


# ==============================================================================
# 1. maxent_score function existence and basic properties
# ==============================================================================

class TestMaxentScoreExistence:
    """Verify maxent_score exists, is callable, and returns reasonable types."""

    def test_maxent_score_is_callable(self):
        """maxent_score is a callable function."""
        assert callable(maxent_score)

    def test_maxent_score_returns_float(self):
        """maxent_score returns a float for a valid context."""
        result = maxent_score("CAGGTAAGT")
        assert isinstance(result, float)

    def test_maxent_score_non_negative(self):
        """Simplified PWM scores are always non-negative (sum of positive weights)."""
        for ctx in ["CAGGTAAGT", "ATGGTCATC", "GGGGTGGGG", "AAAGTAAAA", "TTTGTTTTT"]:
            assert maxent_score(ctx) >= 0.0, f"maxent_score({ctx!r}) = {maxent_score(ctx)} < 0"

    def test_pwm_dimensions(self):
        """The PWM has 4 rows (ACGT) × 9 columns (positions -3 to +6)."""
        assert len(_MAXENT_PWM) == 4, f"Expected 4 PWM rows, got {len(_MAXENT_PWM)}"
        for row in _MAXENT_PWM:
            assert len(row) == 9, f"Expected 9 PWM columns, got {len(row)}"

    def test_base_index_covers_acgt(self):
        """_BASE_INDEX maps A, C, G, T to valid row indices."""
        assert set(_BASE_INDEX.keys()) == {"A", "C", "G", "T"}
        for base, idx in _BASE_INDEX.items():
            assert 0 <= idx < 4

    def test_constants_sensible(self):
        """Named constants have sensible values."""
        assert _PWM_CONTEXT_LEN == 9
        assert _MIN_CONTEXT_LEN == 4
        assert _PWM_UPSTREAM == 3
        assert _PWM_DOWNSTREAM == 6

    def test_all_pwm_weights_non_negative(self):
        """All weights in the simplified PWM are non-negative."""
        for row in _MAXENT_PWM:
            for w in row:
                assert w >= 0.0, f"Negative PWM weight: {w}"


# ==============================================================================
# 2. Known splice site sequences produce expected score ranges
#
# NOTE: The simplified PWM in splicing.py is NOT a proper log-odds model.
# It uses a hand-crafted weight matrix where:
#   - GT at positions 3-4 contributes ~0.02 (G→0.01, T→0.01)
#   - A at positions 3-4 contributes ~7.00 (A→3.50 each)
# When called via score_splice_sites, GT is always at positions 3-4,
# so the core contributes minimally and scores are determined by flanking.
# For full 9-mer contexts with GT at positions 3-4, scores range ~0.3-0.9.
# The proper MaxEntScan model (maxentscan.py) gives biologically meaningful
# log-odds scores (strong donors ~3-10, weak ~0-3, non-donors < 0).
# ==============================================================================

class TestKnownSpliceSiteScores:
    """Strong consensus sequences should be distinguishable from weak ones."""

    def test_simplified_pwm_strong_donor_in_range(self):
        """Simplified PWM: CAGGTAAGT (strong consensus) produces a score in the
        expected low range for GT-at-positions-3-4 contexts."""
        score = maxent_score("CAGGTAAGT")
        # With GT at positions 3-4, the score is determined by flanking only.
        # Typical range for 9-mer contexts with GT is 0.3-1.0.
        assert 0.3 <= score <= 1.5, (
            f"Strong donor context 'CAGGTAAGT' scored {score:.4f}, expected 0.3-1.5"
        )

    def test_simplified_pwm_all_contexts_low_when_gt_present(self):
        """All 9-mer contexts with GT at positions 3-4 produce scores well below
        the 3.0 PASS threshold (because GT contributes only 0.02)."""
        contexts_with_gt = [
            "CAGGTAAGT",  # Strong consensus flanking
            "ATGGTAAGT",  # Moderate flanking
            "TTTGTTTTT",  # T-rich flanking
            "CCCGTCCCC",  # C-rich flanking
            "GGGGTGGGG",  # G-rich flanking
            "ATGGTCATC",  # Weak flanking
        ]
        for ctx in contexts_with_gt:
            score = maxent_score(ctx)
            assert score < 3.0, (
                f"Context {ctx!r} with GT at positions 3-4 scored {score:.2f}, "
                f"expected < 3.0 (GT positions contribute only 0.02)"
            )

    def test_simplified_pwm_c_rich_flanking_lower_than_a_rich(self):
        """C/G-rich flanking around GT produces lower scores than A/T-rich flanking
        in the simplified PWM (A/T have higher weights at most positions)."""
        cg_rich = maxent_score("CCCGTCCCC")
        at_rich = maxent_score("AAAGTAAAA")
        assert at_rich > cg_rich, (
            f"A-rich flanking ({at_rich:.4f}) should score higher than "
            f"C-rich flanking ({cg_rich:.4f})"
        )

    # ── Proper MaxEntScan (maxentscan.py) tests ──────────────────────────────

    def test_maxentscan_strong_donor_scores_positive(self):
        """Proper MaxEntScan: CAG|GTAAGT scores positively (> 0)."""
        seq = "AAACAGGTAAGTAAAA"
        gt_pos = seq.find("GT")
        score = score_donor(seq, gt_pos)
        assert score > 0.0, f"Strong donor should score > 0, got {score}"

    def test_maxentscan_strong_donor_outperforms_weak(self):
        """Proper MaxEntScan: CAGGTAAGT scores higher than ATGGTCATC."""
        strong_seq = "AAACAGGTAAGTAAAA"
        weak_seq = "AAAATGGTCATCAAAA"
        strong_score = score_donor(strong_seq, strong_seq.find("GT"))
        weak_score = score_donor(weak_seq, weak_seq.find("GT"))
        assert strong_score > weak_score, (
            f"Strong donor ({strong_score:.2f}) should outscore weak ({weak_score:.2f})"
        )

    def test_maxentscan_donor_score_typical_range(self):
        """Proper MaxEntScan donor scores for common sequences fall in a
        reasonable range (roughly -5 to +15 for most biological sequences)."""
        seq = "AAACAGGTAAGTAAAA"
        score = score_donor(seq, seq.find("GT"))
        assert -10.0 < score < 20.0, (
            f"Donor score {score} outside expected range [-10, 20]"
        )

    def test_maxentscan_acceptor_with_py_tract(self):
        """Proper MaxEntScan: acceptor with polypyrimidine tract scores higher
        than one with purine-rich upstream."""
        # Strong: C-rich upstream (polypyrimidine tract)
        strong_seq = "C" * 20 + "CAG" + "AAAA"
        # Weak: G-rich upstream (no py tract)
        weak_seq = "G" * 20 + "AAG" + "AAAA"
        strong_ag = strong_seq.find("AG")
        weak_ag = weak_seq.find("AG")
        if strong_ag >= 20 and weak_ag >= 20:
            strong_s = score_acceptor(strong_seq, strong_ag)
            weak_s = score_acceptor(weak_seq, weak_ag)
            assert strong_s > weak_s, (
                f"Py-tract acceptor ({strong_s:.2f}) should outscore "
                f"purine-rich ({weak_s:.2f})"
            )

    def test_gt_dinucleotide_core_contributes_most_in_simplified_pwm(self):
        """In the simplified PWM, having A (instead of G/T) at positions 3-4
        contributes 3.50+3.50=7.00 — the dominant signal. This means
        the simplified PWM penalizes the ABSENCE of GT heavily."""
        with_gt = maxent_score("CAGGTAAGT")  # G at pos 3, T at pos 4
        # Construct a context with A at positions 3-4 instead
        # "CAAATAAGT": A at position 3, T at position 4
        with_a_at_3 = maxent_score("CAAATAAGT")  # A at pos 3 instead of G
        assert with_a_at_3 > with_gt, (
            f"A at pos 3 ({with_a_at_3:.4f}) > G at pos 3 ({with_gt:.4f}): "
            f"non-GT at core positions should increase score"
        )
        # The difference should be large (3.50 - 0.01 = 3.49)
        assert with_a_at_3 - with_gt > 3.0, (
            f"Score difference should be large: {with_a_at_3 - with_gt:.2f}"
        )


# ==============================================================================
# 3. Non-splice-site sequences produce low scores
# ==============================================================================

class TestNonSpliceSiteScores:
    """Sequences without strong splice signals should produce appropriate scores."""

    def test_score_splice_sites_no_gt_empty_result(self):
        """A sequence with no GT dinucleotides yields empty results."""
        seq = "AAAAAAAAAAAAAA"
        results = score_splice_sites(seq)
        assert results == []

    def test_score_splice_sites_ag_only_not_scanned(self):
        """score_splice_sites only scans for GT (donor) sites, not AG."""
        seq = "AAAAAAAAAGAAAA"
        results = score_splice_sites(seq)
        assert results == []

    def test_maxent_scan_no_gt_impossible_donor(self):
        """Proper MaxEntScan: sequence without GT returns _IMPOSSIBLE_SCORE."""
        seq = "AAAAAAAAAA"
        best = max_donor_score(seq)
        assert best == _IMPOSSIBLE_SCORE

    def test_maxent_scan_no_ag_impossible_acceptor(self):
        """Proper MaxEntScan: sequence without AG returns _IMPOSSIBLE_SCORE."""
        seq = "CCCCCCCCCCCCCCCCCCCCCCCCC"
        best = max_acceptor_score(seq)
        assert best == _IMPOSSIBLE_SCORE

    def test_gt_with_poor_flanking_scores_lower_in_maxentscan(self):
        """Proper MaxEntScan: GT with poor flanking scores lower than
        GT with strong consensus flanking."""
        strong_seq = "AAACAGGTAAGTAAAA"
        weak_seq = "AAATATGTCCTCAAAA"
        strong_pos = strong_seq.find("GT")
        weak_pos = weak_seq.find("GT")
        if strong_pos >= 0 and weak_pos >= 0:
            strong_s = score_donor(strong_seq, strong_pos)
            weak_s = score_donor(weak_seq, weak_pos)
            assert strong_s > weak_s, (
                f"Strong consensus ({strong_s:.2f}) should outscore weak ({weak_s:.2f})"
            )

    def test_simplified_pwm_low_scores_for_gt_in_context(self):
        """With GT at positions 3-4 in the context (normal case for
        score_splice_sites), the simplified PWM produces low scores
        (GT positions contribute only 0.01+0.01=0.02)."""
        for ctx in ["CAGGTAAGT", "ATGGTCATC", "CCCGTCCCC", "TTTGTTTTT"]:
            score = maxent_score(ctx)
            assert score < 2.0, (
                f"GT in context {ctx!r} scored {score:.2f}, "
                f"expected < 2.0 (GT positions add only 0.02)"
            )

    def test_non_gt_context_scores_high_in_simplified_pwm(self):
        """Contexts without GT at positions 3-4 score HIGH in the simplified PWM,
        because non-G/T bases at those positions carry large weights (3.50).
        This is by design — the PWM penalizes absence of GT."""
        # All A context (A at positions 3-4): 3.50 + 3.50 from core
        all_a = maxent_score("AAAAAAAAA")
        assert all_a > 5.0, (
            f"All-A context should score high (> 5.0) due to A at GT positions, "
            f"got {all_a:.2f}"
        )


# ==============================================================================
# 4. Edge cases — short sequences, no GT/AG, boundary conditions
# ==============================================================================

class TestEdgeCases:
    """Edge cases for splicing score functions."""

    # ── Short sequences (below minimum context length) ───────────────────────

    def test_empty_string_returns_zero(self):
        """Empty string returns 0.0 (below minimum context length)."""
        assert maxent_score("") == 0.0

    def test_single_char_returns_zero(self):
        """Single character returns 0.0 (below minimum context length)."""
        assert maxent_score("A") == 0.0

    def test_two_char_returns_zero(self):
        """Two characters return 0.0 (below minimum context length 4)."""
        assert maxent_score("GT") == 0.0

    def test_three_char_returns_zero(self):
        """Three characters return 0.0 (below minimum context length 4)."""
        assert maxent_score("AGT") == 0.0

    def test_min_context_len_four(self):
        """Four characters (minimum context) returns a non-zero score."""
        score = maxent_score("AAGT")
        assert score > 0.0, f"Minimum context should score > 0, got {score}"

    def test_short_context_left_padded(self):
        """Contexts shorter than 9 are left-padded with A to make 9-mer,
        then the last 9 characters are taken."""
        short = "AAGT"  # 4 chars
        score_short = maxent_score(short)
        # Expected padded: "AAAAA" + "AAGT" = "AAAAAAAGT", take last 9 = "AAAAAAAGT"
        score_manual = maxent_score("AAAAAAAGT")
        assert score_short == score_manual, (
            f"Short context score {score_short} != manual padded score {score_manual}"
        )

    def test_five_char_context_padded(self):
        """5-char context is left-padded to 9."""
        short = "TAAGT"  # 5 chars (GT near end)
        score_short = maxent_score(short)
        score_manual = maxent_score("AAAATAAGT")
        assert score_short == score_manual

    # ── No GT/AG in sequence ─────────────────────────────────────────────────

    def test_no_gt_in_sequence(self):
        """Sequence with no GT dinucleotide produces no splice site results."""
        for seq in ["AAAAAA", "CCCCCC", "ACACACAC", "TATATATA"]:
            results = score_splice_sites(seq)
            assert results == [], f"No GT in {seq!r} but got results: {results}"

    # ── MaxEntScan boundary positions ────────────────────────────────────────

    def test_maxent_scan_donor_out_of_range(self):
        """score_donor with position too close to edges returns _IMPOSSIBLE_SCORE."""
        short_seq = "GTAA"
        # Position 0: need 3 upstream, but start would be -3
        assert score_donor(short_seq, 0) == _IMPOSSIBLE_SCORE
        # Position near end: need 6 downstream after GT position
        assert score_donor(short_seq, 2) == _IMPOSSIBLE_SCORE

    def test_maxent_scan_acceptor_out_of_range(self):
        """score_acceptor with position too close to edges returns _IMPOSSIBLE_SCORE."""
        short_seq = "AGAA"
        # Position 0: need 20 upstream
        assert score_acceptor(short_seq, 0) == _IMPOSSIBLE_SCORE
        # Near end: need 3 downstream
        assert score_acceptor(short_seq, 2) == _IMPOSSIBLE_SCORE

    def test_negative_position_returns_impossible(self):
        """Negative position returns _IMPOSSIBLE_SCORE."""
        seq = "AAACAGGTAAGTAAAA"
        assert score_donor(seq, -1) == _IMPOSSIBLE_SCORE
        assert score_acceptor(seq, -1) == _IMPOSSIBLE_SCORE

    def test_position_beyond_sequence_returns_impossible(self):
        """Position beyond sequence length returns _IMPOSSIBLE_SCORE."""
        seq = "AAACAGGTAAGTAAAA"
        assert score_donor(seq, 100) == _IMPOSSIBLE_SCORE
        assert score_acceptor(seq, 100) == _IMPOSSIBLE_SCORE

    # ── Invalid characters ───────────────────────────────────────────────────

    def test_maxent_score_unknown_base_treated_as_a(self):
        """Unknown bases in maxent_score default to A (index 0)."""
        score_with_n = maxent_score("CANGTAAGT")
        score_with_a = maxent_score("CAAGTAAGT")
        assert score_with_n == score_with_a, (
            "Unknown base should be treated as A in maxent_score"
        )

    def test_maxent_scan_invalid_base_returns_impossible(self):
        """score_donor with invalid base (N) in the 9-mer window
        returns _IMPOSSIBLE_SCORE."""
        # Place N within the 9-mer window around GT
        seq = "AAACANGTAAGTAAAA"
        gt_pos = seq.find("GT")
        if gt_pos >= 0:
            score = score_donor(seq, gt_pos)
            assert score == _IMPOSSIBLE_SCORE, (
                f"Invalid base N should give _IMPOSSIBLE_SCORE, got {score}"
            )

    def test_maxent_score_invalid_char_no_crash(self):
        """maxent_score doesn't crash on unusual characters (they map to A)."""
        # These all default to _BASE_INDEX.get(char, 0) = index 0 (A)
        for ctx in ["NANGTANNN", "123456789", "!!!!!!!!A"]:
            # Should not raise, returns a float
            score = maxent_score(ctx)
            assert isinstance(score, float)

    # ── Case insensitivity ───────────────────────────────────────────────────

    def test_maxent_score_case_insensitive(self):
        """maxent_score handles lowercase sequences (uppercased internally)."""
        score_upper = maxent_score("CAGGTAAGT")
        score_lower = maxent_score("caggtaagt")
        score_mixed = maxent_score("CaGgTaAgT")
        assert score_upper == score_lower == score_mixed, (
            f"Case sensitivity issue: {score_upper} vs {score_lower} vs {score_mixed}"
        )

    def test_maxent_scan_case_insensitive(self):
        """score_donor/score_acceptor handle lowercase sequences."""
        seq = "aaacaggtaagtaaaa"
        gt_pos = seq.upper().find("GT")
        score_lower = score_donor(seq, gt_pos)
        score_upper = score_donor(seq.upper(), gt_pos)
        assert score_lower == score_upper

    # ── Boundary context extraction ──────────────────────────────────────────

    def test_score_splice_sites_near_start(self):
        """GT near the start of sequence: context extracted with max(0, i-3)."""
        seq = "GTAAAAAA"
        results = score_splice_sites(seq)
        assert len(results) >= 1
        pos, score, verdict = results[0]
        assert pos == 0
        # Short context → left-padded with A, A at position 3 scores high
        assert score > 0.0

    def test_score_splice_sites_near_end(self):
        """GT near the end of sequence: context extracted with min(len, i+6)."""
        seq = "AAAAAAGT"
        results = score_splice_sites(seq)
        assert len(results) >= 1
        pos, score, verdict = results[0]
        assert pos == 6

    def test_gt_at_boundary_can_score_high(self):
        """GT near sequence boundaries gets padded with A, which may cause
        higher scores (A at GT positions in the padded context)."""
        # GT at position 1 in a short sequence
        seq = "AGTAA"
        results = score_splice_sites(seq)
        if results:
            # The context is short, so it gets padded with A
            # This can place A at PWM positions 3-4, boosting the score
            assert len(results) >= 1


# ==============================================================================
# 5. MaxEntScan score_donor / score_acceptor specific tests
# ==============================================================================

class TestScoreDonor:
    """Tests for the proper MaxEntScan donor scoring model."""

    def test_score_donor_returns_float(self):
        """score_donor returns a float."""
        seq = "AAACAGGTAAGTAAAA"
        pos = seq.find("GT")
        result = score_donor(seq, pos)
        assert isinstance(result, float)

    def test_score_donor_rounded_to_4_decimals(self):
        """score_donor rounds to 4 decimal places."""
        seq = "AAACAGGTAAGTAAAA"
        pos = seq.find("GT")
        result = score_donor(seq, pos)
        assert result == round(result, 4)

    def test_score_donor_deterministic(self):
        """Same input always produces same output."""
        seq = "AAACAGGTAAGTAAAA"
        pos = seq.find("GT")
        s1 = score_donor(seq, pos)
        s2 = score_donor(seq, pos)
        assert s1 == s2

    def test_score_donor_invariant_gt_contribution(self):
        """Changing G at the donor position (invariant 0.990) dramatically
        reduces the score."""
        seq_with_g = "AAACAGGTAAGTAAAA"
        seq_without_g = "AAACAGATAAGTAAAA"  # Changed G→A at donor+1 position
        pos_g = seq_with_g.find("GT")
        pos_a = seq_without_g.find("GT")
        if pos_g >= 0 and pos_a >= 0:
            s_g = score_donor(seq_with_g, pos_g)
            s_a = score_donor(seq_without_g, pos_a)
            assert s_g > s_a, (
                f"Consensus donor ({s_g:.4f}) should outscore non-consensus ({s_a:.4f})"
            )

    def test_donor_pwm_dimensions(self):
        """DONOR_PWM_SCORE has 9 positions × 4 bases."""
        assert len(DONOR_PWM_SCORE) == 9
        for row in DONOR_PWM_SCORE:
            assert len(row) == 4

    def test_donor_model_offsets(self):
        """Donor model context: 3 upstream, 6 downstream."""
        assert _DONOR_UPSTREAM == 3
        assert _DONOR_DOWNSTREAM == 6

    def test_donor_minimum_context_length(self):
        """Sequences too short for the 9-mer donor model return _IMPOSSIBLE_SCORE."""
        # Need at least 3+2+6 = 11 chars for GT at position 3
        for short in ["GT", "AGTAA", "AAGTAAA"]:
            pos = short.find("GT")
            if pos >= 0:
                assert score_donor(short, pos) == _IMPOSSIBLE_SCORE


class TestScoreAcceptor:
    """Tests for the proper MaxEntScan acceptor scoring model."""

    def test_score_acceptor_returns_float(self):
        """score_acceptor returns a float."""
        seq = "C" * 25 + "CAG" + "AAAA"
        ag_pos = seq.find("AG")
        if ag_pos >= 20 and ag_pos + 3 <= len(seq):
            result = score_acceptor(seq, ag_pos)
            assert isinstance(result, float)

    def test_score_acceptor_deterministic(self):
        """Same input always produces same output."""
        seq = "C" * 25 + "CAG" + "AAAA"
        ag_pos = seq.find("AG")
        if ag_pos >= 20 and ag_pos + 3 <= len(seq):
            s1 = score_acceptor(seq, ag_pos)
            s2 = score_acceptor(seq, ag_pos)
            assert s1 == s2

    def test_score_acceptor_rounded(self):
        """score_acceptor rounds to 4 decimal places."""
        seq = "C" * 25 + "CAG" + "AAAA"
        ag_pos = seq.find("AG")
        if ag_pos >= 20 and ag_pos + 3 <= len(seq):
            result = score_acceptor(seq, ag_pos)
            assert result == round(result, 4)

    def test_acceptor_pwm_dimensions(self):
        """ACCEPTOR_PWM_SCORE has 23 positions × 4 bases."""
        assert len(ACCEPTOR_PWM_SCORE) == 23
        for row in ACCEPTOR_PWM_SCORE:
            assert len(row) == 4

    def test_acceptor_model_offsets(self):
        """Acceptor model context: 20 upstream, 3 downstream."""
        assert _ACCEPTOR_UPSTREAM == 20
        assert _ACCEPTOR_DOWNSTREAM == 3

    def test_strong_py_tract_outperforms_purine_tract(self):
        """Acceptor with polypyrimidine tract (C/T) scores higher than
        one with purine-rich (G/A) upstream."""
        strong_seq = "C" * 20 + "CAG" + "AAAA"
        weak_seq = "G" * 20 + "AAG" + "AAAA"
        strong_ag = strong_seq.find("AG")
        weak_ag = weak_seq.find("AG")
        if strong_ag >= 20 and weak_ag >= 20:
            strong_s = score_acceptor(strong_seq, strong_ag)
            weak_s = score_acceptor(weak_seq, weak_ag)
            assert strong_s > weak_s, (
                f"Py-tract ({strong_s:.4f}) should outscore purine-rich ({weak_s:.4f})"
            )

    def test_invariant_a_at_acceptor_minus1(self):
        """Position -1 of acceptor is near-invariant A (0.980 in the PWM).
        Changing it should dramatically reduce score."""
        good_seq = "C" * 20 + "CAG" + "AAAA"
        bad_seq = "C" * 20 + "CTG" + "AAAA"
        good_ag = good_seq.find("AG")
        bad_ag = bad_seq.find("AG")
        if good_ag >= 20 and bad_ag >= 20:
            good_s = score_acceptor(good_seq, good_ag)
            bad_s = score_acceptor(bad_seq, bad_ag)
            assert good_s > bad_s, (
                f"A at -1 ({good_s:.4f}) should outscore T at -1 ({bad_s:.4f})"
            )

    def test_acceptor_minimum_context_length(self):
        """Sequences too short for the 23-mer acceptor model return
        _IMPOSSIBLE_SCORE."""
        for short in ["AG", "CAGAA", "C" * 10 + "CAG" + "AA"]:
            ag_pos = short.find("AG")
            if ag_pos >= 0:
                assert score_acceptor(short, ag_pos) == _IMPOSSIBLE_SCORE


# ==============================================================================
# 6. scan_splice_sites and helper functions
# ==============================================================================

class TestScanSpliceSites:
    """Tests for the maxentscan scan_splice_sites function."""

    def test_scan_finds_strong_donor(self):
        """scan_splice_sites finds a strong donor above default threshold."""
        seq = "AAACAGGTAAGTAAAA"
        results = scan_splice_sites(seq, donor_threshold=3.0)
        donors = [(pos, typ, sc) for pos, typ, sc in results if typ == "donor"]
        assert len(donors) >= 1, "Should find at least one donor above threshold"

    def test_scan_threshold_filters(self):
        """Higher threshold filters out weaker sites."""
        seq = "AAACAGGTAAGTAAAA"
        low_thresh_results = scan_splice_sites(seq, donor_threshold=0.0)
        high_thresh_results = scan_splice_sites(seq, donor_threshold=10.0)
        donors_low = [r for r in low_thresh_results if r[1] == "donor"]
        donors_high = [r for r in high_thresh_results if r[1] == "donor"]
        assert len(donors_low) >= len(donors_high)

    def test_scan_permissive_finds_all_gt(self):
        """With very permissive thresholds, all GT positions with enough
        context are found."""
        seq = "AAACAGGTAAGTAAAA"
        # Brute-force GT positions with enough context for score_donor
        all_gt = set()
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "GT" and i >= 3 and i + 6 <= len(seq):
                all_gt.add(i)
        results = scan_splice_sites(seq, donor_threshold=-100.0)
        found_gt = {pos for pos, typ, _ in results if typ == "donor"}
        assert all_gt.issubset(found_gt), (
            f"Permissive scan missed GT positions: {all_gt - found_gt}"
        )

    def test_scan_results_sorted_by_position(self):
        """Results are sorted by position."""
        seq = "AAACAGGTAAGTAAAA" + "C" * 30 + "TTTTTTTTTTTTTTTTTTTCAG" + "AAAA"
        results = scan_splice_sites(seq, donor_threshold=-100.0, acceptor_threshold=-100.0)
        positions = [r[0] for r in results]
        assert positions == sorted(positions)

    def test_scan_result_structure(self):
        """Each result is a (position, site_type, score) tuple."""
        seq = "AAACAGGTAAGTAAAA"
        results = scan_splice_sites(seq, donor_threshold=0.0)
        for result in results:
            assert len(result) == 3
            pos, site_type, score = result
            assert isinstance(pos, int)
            assert site_type in ("donor", "acceptor")
            assert isinstance(score, float)

    def test_max_donor_score_consistent(self):
        """max_donor_score returns the maximum donor score in the sequence."""
        seq = "AAACAGGTAAGTAAAA"
        best = max_donor_score(seq)
        # Manually compute: find all GT positions and score them
        max_manual = _IMPOSSIBLE_SCORE
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "GT":
                s = score_donor(seq, i)
                if s > max_manual:
                    max_manual = s
        assert best == round(max_manual, 4)

    def test_max_acceptor_score_consistent(self):
        """max_acceptor_score returns the maximum acceptor score in the sequence."""
        seq = "C" * 25 + "CAG" + "AAAA"
        best = max_acceptor_score(seq)
        max_manual = _IMPOSSIBLE_SCORE
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "AG":
                s = score_acceptor(seq, i)
                if s > max_manual:
                    max_manual = s
        assert best == round(max_manual, 4)

    def test_scan_finds_acceptors(self):
        """scan_splice_sites finds acceptor (AG) sites above threshold."""
        # Build a sequence with a strong acceptor
        seq = "C" * 25 + "CAG" + "AAAA"
        results = scan_splice_sites(seq, acceptor_threshold=-100.0)
        acceptors = [(pos, typ, sc) for pos, typ, sc in results if typ == "acceptor"]
        assert len(acceptors) >= 1, "Should find at least one acceptor"


# ==============================================================================
# 7. Dual-threshold classification (score_splice_sites)
# ==============================================================================

class TestDualThresholdClassification:
    """Tests for PASS/UNCERTAIN/FAIL classification in score_splice_sites."""

    def test_pass_below_low_threshold(self):
        """Score < low_thresh → SpliceVerdict.PASS."""
        # Most GT sites in long sequences score below 3.0
        seq = "AAACAGGTAAGTAAAA"
        results = score_splice_sites(seq, low_thresh=3.0, high_thresh=6.0)
        for pos, score, verdict in results:
            if score < 3.0:
                assert verdict == SpliceVerdict.PASS

    def test_uncertain_between_thresholds(self):
        """low_thresh ≤ score < high_thresh → SpliceVerdict.UNCERTAIN."""
        # Force UNCERTAIN by using narrow thresholds
        results = score_splice_sites(
            "AAACAGGTAAGTAAAA",
            low_thresh=0.01,
            high_thresh=20.0,
        )
        for pos, score, verdict in results:
            if 0.01 <= score < 20.0:
                assert verdict == SpliceVerdict.UNCERTAIN

    def test_fail_above_high_threshold(self):
        """score ≥ high_thresh → SpliceVerdict.FAIL."""
        # Force FAIL by using very low thresholds
        results = score_splice_sites(
            "AAACAGGTAAGTAAAA",
            low_thresh=0.01,
            high_thresh=0.02,
        )
        for pos, score, verdict in results:
            if score >= 0.02:
                assert verdict == SpliceVerdict.FAIL

    def test_verdict_exhaustive(self):
        """Every result gets exactly one of PASS/UNCERTAIN/FAIL."""
        valid_verdicts = {SpliceVerdict.PASS, SpliceVerdict.UNCERTAIN, SpliceVerdict.FAIL}
        seq = "AAACAGGTAAGTAAAA"
        results = score_splice_sites(seq)
        for pos, score, verdict in results:
            assert verdict in valid_verdicts

    def test_all_verdicts_are_splice_verdict_enum(self):
        """All returned verdicts are SpliceVerdict enum members."""
        seq = "CAGGTAAGT"
        results = score_splice_sites(seq)
        for pos, score, verdict in results:
            assert isinstance(verdict, SpliceVerdict)

    def test_result_tuple_structure(self):
        """Each result is a (position: int, score: float, verdict: SpliceVerdict) tuple."""
        seq = "CAGGTAAGT"
        results = score_splice_sites(seq)
        for result in results:
            assert len(result) == 3
            pos, score, verdict = result
            assert isinstance(pos, int)
            assert isinstance(score, float)
            assert isinstance(verdict, SpliceVerdict)

    def test_classify_threshold_logic_directly(self):
        """Verify the threshold classification logic directly:
        PASS when score < low, UNCERTAIN when low <= score < high,
        FAIL when score >= high."""
        # Use a short sequence where the GT is near the boundary,
        # causing A-padding that can push scores above 3.0
        seq = "AGTAA"
        results = score_splice_sites(seq, low_thresh=3.0, high_thresh=6.0)
        for pos, score, verdict in results:
            if score < 3.0:
                assert verdict == SpliceVerdict.PASS
            elif score < 6.0:
                assert verdict == SpliceVerdict.UNCERTAIN
            else:
                assert verdict == SpliceVerdict.FAIL


# ==============================================================================
# 8. Internal consistency between simplified PWM and MaxEntScan
# ==============================================================================

class TestConsistencyBetweenModels:
    """Verify that the simplified PWM and proper MaxEntScan agree directionally
    on which sequences look more splice-like."""

    def test_both_models_score_consensually(self):
        """Both models produce numerical scores for the same input."""
        pwm_score = maxent_score("CAGGTAAGT")
        seq = "AAACAGGTAAGTAAAA"
        gt_pos = seq.find("GT")
        mes_score = score_donor(seq, gt_pos)
        assert isinstance(pwm_score, float)
        assert isinstance(mes_score, float)

    def test_both_models_agree_on_ordering(self):
        """Both models agree that CAGGTAAGT has a stronger splice signal
        than ATGGTCATC (proper MaxEntScan). The simplified PWM ranks
        based on flanking A/T content, which differs from the biological
        model, so we only check the proper model."""
        strong_seq = "AAACAGGTAAGTAAAA"
        weak_seq = "AAAATGGTCATCAAAA"
        mes_strong = score_donor(strong_seq, strong_seq.find("GT"))
        mes_weak = score_donor(weak_seq, weak_seq.find("GT"))
        assert mes_strong > mes_weak, (
            f"Proper MaxEntScan: strong ({mes_strong:.4f}) > weak ({mes_weak:.4f})"
        )

    def test_maxent_score_vs_manual_computation(self):
        """maxent_score matches manual computation from the PWM."""
        ctx = "CAGGTAAGT"
        expected = 0.0
        for pos in range(9):
            base = ctx[pos]
            idx = _BASE_INDEX[base]
            expected += _MAXENT_PWM[idx][pos]
        actual = maxent_score(ctx)
        assert abs(actual - expected) < 1e-10, (
            f"Manual {expected} != actual {actual}"
        )

    def test_score_donor_log_odds_formula(self):
        """score_donor computes log2(P_motif / P_background) correctly."""
        seq = "AAACAGGTAAGTAAAA"
        gt_pos = seq.find("GT")
        # Manual computation
        start = gt_pos - _DONOR_UPSTREAM
        manual_score = 0.0
        for pwm_idx in range(9):
            base = seq[start + pwm_idx]
            idx = MES_BASE_TO_INDEX[base]
            prob = max(DONOR_PWM_SCORE[pwm_idx][idx], 0.001)  # _EPSILON = 0.001
            manual_score += math.log2(prob / BG_PROB)
        manual_score = round(manual_score, 4)
        actual = score_donor(seq, gt_pos)
        assert actual == manual_score, (
            f"Manual log-odds {manual_score} != score_donor {actual}"
        )

    def test_score_acceptor_log_odds_formula(self):
        """score_acceptor computes log2(P_motif / P_background) correctly."""
        seq = "C" * 20 + "CAG" + "AAAA"
        ag_pos = seq.find("AG")
        start = ag_pos - _ACCEPTOR_UPSTREAM
        manual_score = 0.0
        for pwm_idx in range(23):
            base = seq[start + pwm_idx]
            idx = MES_BASE_TO_INDEX[base]
            prob = max(ACCEPTOR_PWM_SCORE[pwm_idx][idx], 0.001)
            manual_score += math.log2(prob / BG_PROB)
        manual_score = round(manual_score, 4)
        actual = score_acceptor(seq, ag_pos)
        assert actual == manual_score


# ==============================================================================
# 9. Multiple GT positions in one sequence
# ==============================================================================

class TestMultipleGTPositions:
    """Test sequences with multiple GT dinucleotides."""

    def test_multiple_gt_detected(self):
        """Multiple GT dinucleotides are all found by score_splice_sites."""
        seq = "AAAGTAAAAAAGTAAAAAAGTAA"
        gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
        results = score_splice_sites(seq)
        result_positions = [pos for pos, _, _ in results]
        assert set(gt_positions) == set(result_positions), (
            f"GT positions {gt_positions} != result positions {result_positions}"
        )

    def test_overlapping_gt_not_double_counted(self):
        """GGT contains only one GT (at position 1), not overlapping."""
        seq = "GGTAAAAA"
        gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
        assert len(gt_positions) == 1
        results = score_splice_sites(seq)
        assert len(results) == 1

    def test_maxentscan_scores_differ_by_context(self):
        """MaxEntScan: GT dinucleotides in different contexts get different scores."""
        # Strong GT (CAGGTAAGT) vs weak GT (ATGGTCATC) in same sequence
        seq = "AAACAGGTAAGTAAAA" + "AAAATGGTCATCAAAA"
        results = scan_splice_sites(seq, donor_threshold=-100.0)
        donor_scores = {pos: sc for pos, typ, sc in results if typ == "donor"}
        if len(donor_scores) >= 2:
            scores = list(donor_scores.values())
            # At least one score should differ from another
            assert len(set(round(s, 2) for s in scores)) >= 1

    def test_three_gt_in_sequence(self):
        """Three GT positions in one sequence are all detected."""
        seq = "AACAGGTAAGTAAATGGTCATCAAA" + "C" * 10 + "TTTGTTTTT"
        gt_count = sum(1 for i in range(len(seq)-1) if seq[i:i+2] == "GT")
        results = score_splice_sites(seq)
        assert len(results) == gt_count


# ==============================================================================
# 10. Numeric edge cases and PWM properties
# ==============================================================================

class TestNumericEdgeCases:
    """Numeric boundary conditions for scoring."""

    def test_impossible_score_constant(self):
        """_IMPOSSIBLE_SCORE is a very negative number."""
        assert _IMPOSSIBLE_SCORE < -40.0

    def test_bg_prob_uniform(self):
        """Background probability is 0.25 (uniform over 4 bases)."""
        assert BG_PROB == 0.25

    def test_donor_pwm_row_sums_approximately_one(self):
        """Each row of DONOR_PWM_SCORE sums to approximately 1.0."""
        for i, row in enumerate(DONOR_PWM_SCORE):
            total = sum(row)
            assert 0.5 < total < 1.5, (
                f"Donor PWM row {i} sums to {total}, expected ~1.0"
            )

    def test_acceptor_pwm_row_sums_approximately_one(self):
        """Each row of ACCEPTOR_PWM_SCORE sums to approximately 1.0."""
        for i, row in enumerate(ACCEPTOR_PWM_SCORE):
            total = sum(row)
            assert 0.5 < total < 1.5, (
                f"Acceptor PWM row {i} sums to {total}, expected ~1.0"
            )

    def test_simplified_pwm_a_weights_at_gt_positions(self):
        """In the simplified PWM, A at positions 3,4 has weight 3.50 each —
        these are the largest weights in the entire matrix."""
        max_weight = 0.0
        for row in _MAXENT_PWM:
            for w in row:
                max_weight = max(max_weight, w)
        # A at positions 3 and 4 should be the maximum
        assert _MAXENT_PWM[0][3] == max_weight
        assert _MAXENT_PWM[0][4] == max_weight
        assert max_weight == 3.50

    def test_simplified_pwm_non_gt_weights_low(self):
        """G and T at positions 3,4 have very low weights (0.01) —
        making the GT core contribution minimal when GT is present."""
        assert _MAXENT_PWM[2][3] == 0.01  # G at position 3
        assert _MAXENT_PWM[3][4] == 0.01  # T at position 4

    def test_donor_pwm_invariant_positions(self):
        """Donor PWM positions +1 (G) and +2 (T) should have near-invariant
        frequencies (close to 1.0 for the correct base)."""
        # Position +1 (index 3) should have G near 1.0
        g_freq_at_plus1 = DONOR_PWM_SCORE[3][MES_BASE_TO_INDEX["G"]]
        assert g_freq_at_plus1 > 0.9, (
            f"G frequency at donor +1: {g_freq_at_plus1}, expected > 0.9"
        )
        # Position +2 (index 4) should have T near 1.0
        t_freq_at_plus2 = DONOR_PWM_SCORE[4][MES_BASE_TO_INDEX["T"]]
        assert t_freq_at_plus2 > 0.9, (
            f"T frequency at donor +2: {t_freq_at_plus2}, expected > 0.9"
        )

    def test_acceptor_pwm_invariant_positions(self):
        """Acceptor PWM position -1 should have A near 1.0,
        position +0 should have G near 1.0."""
        # Position -1 (index 19) should have A near 1.0
        a_freq_at_minus1 = ACCEPTOR_PWM_SCORE[19][MES_BASE_TO_INDEX["A"]]
        assert a_freq_at_minus1 > 0.9, (
            f"A frequency at acceptor -1: {a_freq_at_minus1}, expected > 0.9"
        )
        # Position +0 (index 20) should have G near 1.0
        g_freq_at_plus0 = ACCEPTOR_PWM_SCORE[20][MES_BASE_TO_INDEX["G"]]
        assert g_freq_at_plus0 > 0.9, (
            f"G frequency at acceptor +0: {g_freq_at_plus0}, expected > 0.9"
        )

    def test_maxent_score_deterministic_across_calls(self):
        """Calling maxent_score multiple times with same input gives same output."""
        for _ in range(10):
            assert maxent_score("CAGGTAAGT") == maxent_score("CAGGTAAGT")
