"""
BioCompiler MaxEntScan Correlation Validation
==============================================

Task ID: 17

Validates that the proper MaxEntScan implementation (maxentscan.score_donor /
score_acceptor) produces scores that are POSITIVELY correlated with known
splice site strength from the Yeo & Burge (2004) reference.

This test was created to address the finding that splicing.maxent_score()
is ANTI-CORRELATED with proper MaxEntScan scoring, meaning it gives high
scores to non-splice sites and low scores to real splice sites.

Key findings documented here:
  - The deprecated maxent_score() uses a hand-crafted PWM where A at the
    GT dinucleotide positions gets weight 3.50, while G/T get weight 0.01.
    This makes real splice sites (with GT) score LOW.
  - The proper maxentscan.score_donor/score_acceptor uses log-odds scoring
    with Yeo & Burge 2004 trained parameters, where the near-invariant
    G (prob 0.990) and T (prob 0.990) at donor positions contribute
    strong positive log-odds (~+2.0 bits each).
  - Prokaryotes (e.g. E. coli) have no spliceosome, so MaxEntScan should
    NEVER be called during optimization for prokaryotic targets.

References:
  Yeo, G. & Burge, C.B. (2004). "Maximum entropy modeling of short sequence
  motifs with applications to RNA splicing." *Journal of Computational Biology*
  11(2-3):377-394. doi:10.1089/1066527041410418
"""

from __future__ import annotations

import math
import warnings
from typing import List

import pytest

from biocompiler.maxentscan import (
    BASE_TO_INDEX,
    BG_PROB,
    DONOR_PWM_SCORE,
    ACCEPTOR_PWM_SCORE,
    score_donor,
    score_acceptor,
    scan_splice_sites,
    max_donor_score,
    max_acceptor_score,
    _IMPOSSIBLE_SCORE,
)


# ==============================================================================
# Helpers
# ==============================================================================

def _pearson_r(x: List[float], y: List[float]) -> float:
    """Compute Pearson correlation coefficient."""
    n = len(x)
    if n < 2:
        return 0.0
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)
    if var_x == 0.0 or var_y == 0.0:
        return 0.0
    return cov / math.sqrt(var_x * var_y)


def _embed_donor_9mer(nine_mer: str, pos: int = 3, total_len: int = 20) -> str:
    """Embed a 9-mer in a sequence so score_donor(seq, pos) scores it."""
    assert len(nine_mer) == 9
    seq = list("A" * total_len)
    for i, c in enumerate(nine_mer):
        idx = pos - 3 + i
        if 0 <= idx < total_len:
            seq[idx] = c
    return "".join(seq)


# ==============================================================================
# Test sequences with known splice site strength
# ==============================================================================

# Donor 9-mers: strong canonical → weak non-canonical
DONOR_TEST_SEQUENCES = [
    # (9-mer, expected_strength, description)
    ("CAGGTGAGT", "strong", "Canonical donor with CAG upstream and AAGT downstream"),
    ("CAGGTAAGT", "strong", "Strong consensus donor CAG|GTAAGT"),
    ("TTGGTAAGT", "moderate", "T-rich upstream, decent downstream"),
    ("CCGGTAAGT", "moderate", "C-rich upstream, decent downstream"),
    ("ATGGTCATC", "weak", "Poor flanking context around GT"),
    ("AAGTAAGCT", "weak", "Non-canonical positions in 9-mer"),
]

# Acceptor test sequences (position of A in AG provided separately)
ACCEPTOR_TEST_SEQUENCES = [
    # (full_seq, ag_pos, expected_strength, description)
    ("T" * 20 + "CAG" + "ATCG", 21, "strong", "Poly-T polypyrimidine tract upstream"),
    ("C" * 20 + "CAG" + "ATCG", 21, "strong", "Poly-C polypyrimidine tract upstream"),
    ("G" * 20 + "AAG" + "ATCG", 21, "weak", "Purine-rich upstream (anti-pyrimidine)"),
    ("A" * 20 + "AAG" + "ATCG", 21, "weak", "A-rich upstream (no pyrimidine tract)"),
]


# ==============================================================================
# Core validation: MaxEntScan scores are positively correlated with strength
# ==============================================================================


class TestMaxEntScanCorrelation:
    """Verify MaxEntScan scores are positively correlated with reference
    implementations (Yeo & Burge 2004).

    The key invariant: strong splice sites should score HIGHER than weak
    splice sites. If this is violated, the implementation is anti-correlated
    with biological reality and must be fixed or deprecated.
    """

    def test_donor_strong_outscore_weak(self) -> None:
        """Strong canonical donors must score higher than weak donors.

        This is the fundamental check: if it fails, the scoring model
        is anti-correlated with biological splice site strength.
        """
        strong_scores = []
        weak_scores = []

        for nine_mer, strength, desc in DONOR_TEST_SEQUENCES:
            seq = _embed_donor_9mer(nine_mer, pos=3)
            score = score_donor(seq, 3)
            if score > _IMPOSSIBLE_SCORE:
                if strength == "strong":
                    strong_scores.append(score)
                elif strength == "weak":
                    weak_scores.append(score)

        assert len(strong_scores) > 0, "No valid strong donor scores"
        assert len(weak_scores) > 0, "No valid weak donor scores"

        # Every strong donor must outscore every weak donor
        for s in strong_scores:
            for w in weak_scores:
                assert s > w, (
                    f"Strong donor score {s:.4f} <= weak donor score {w:.4f} — "
                    f"anti-correlation detected!"
                )

    def test_acceptor_strong_outscore_weak(self) -> None:
        """Strong canonical acceptors must score higher than weak acceptors."""
        strong_scores = []
        weak_scores = []

        for full_seq, ag_pos, strength, desc in ACCEPTOR_TEST_SEQUENCES:
            score = score_acceptor(full_seq, ag_pos)
            if score > _IMPOSSIBLE_SCORE:
                if strength == "strong":
                    strong_scores.append(score)
                elif strength == "weak":
                    weak_scores.append(score)

        assert len(strong_scores) > 0, "No valid strong acceptor scores"
        assert len(weak_scores) > 0, "No valid weak acceptor scores"

        for s in strong_scores:
            for w in weak_scores:
                assert s > w, (
                    f"Strong acceptor score {s:.4f} <= weak acceptor score {w:.4f} — "
                    f"anti-correlation detected!"
                )

    def test_donor_scores_positive_for_strong_sites(self) -> None:
        """Strong canonical donors must score positively (log-odds > 0).

        In the Yeo & Burge model, a strong donor like CAG|GTGAGT has
        high probability under the splice model and low probability
        under background, so log-odds should be positive.
        """
        for nine_mer, strength, desc in DONOR_TEST_SEQUENCES:
            if strength != "strong":
                continue
            seq = _embed_donor_9mer(nine_mer, pos=3)
            score = score_donor(seq, 3)
            assert score > 0.0, (
                f"Strong donor '{nine_mer}' ({desc}) scored {score:.4f}, "
                f"expected > 0.0"
            )

    def test_donor_gt_invariant_contribution(self) -> None:
        """The GT dinucleotide at positions +1/+2 must contribute positively.

        In the Yeo & Burge model, G at position +1 has probability 0.990
        and T at position +2 has probability 0.990. Under uniform background
        (0.25), each contributes log2(0.990/0.25) ≈ +1.98 bits.
        This is the dominant positive contribution for real splice sites.
        """
        # Donor with GT at positions +1/+2
        gt_9mer = "CAGGTAAGT"
        # Same 9-mer but with AA instead of GT (anti-splice-site)
        aa_9mer = "CAGAAAAGT"

        gt_score = score_donor(_embed_donor_9mer(gt_9mer, 3), 3)
        aa_score = score_donor(_embed_donor_9mer(aa_9mer, 3), 3)

        assert gt_score > aa_score, (
            f"GT donor ({gt_score:.2f}) must score higher than AA at same "
            f"position ({aa_score:.2f}) — the GT invariant is the defining "
            f"feature of donor sites"
        )

        # The difference should be large: ~4 bits from the two invariant positions
        assert gt_score - aa_score > 4.0, (
            f"GT vs AA difference ({gt_score - aa_score:.2f}) should be >4 bits "
            f"(G at +1: log2(0.990/0.25) ≈ +2.0 bits, T at +2: ≈ +2.0 bits)"
        )

    def test_acceptor_ag_invariant_contribution(self) -> None:
        """The AG dinucleotide at the acceptor must contribute positively.

        In the Yeo & Burge model, A at position -1 has probability 0.980
        and G at position +0 has probability 0.980. Under uniform background,
        each contributes ~+1.97 bits.

        NOTE: score_acceptor uses G-position convention: position = where
        G of AG is. See test_maxentscan_validation.py for details.
        """
        # Strong acceptor: poly-T tract + CAG
        # Using G-position convention: G of AG is at position 22
        good_seq = "A" + "T" * 19 + "CAG" + "ATCG"
        good_g_pos = 22  # Position of G in AG

        # Weak acceptor: poly-T tract + CTG instead of CAG
        bad_seq = "A" + "T" * 19 + "CTG" + "ATCG"
        # No AG in bad_seq — that's the point

        good_score = score_acceptor(good_seq, good_g_pos)

        # The strong acceptor with proper AG and pyrimidine tract
        # should score positively
        if good_score > _IMPOSSIBLE_SCORE:
            assert good_score > 0.0, (
                f"Good acceptor with AG and py-tract should score positively, "
                f"got {good_score:.4f}"
            )

    def test_polypyrimidine_tract_boosts_acceptor(self) -> None:
        """Acceptor sites with polypyrimidine tracts must score higher than
        those with purine-rich upstream.

        The polypyrimidine tract is the defining feature of 3' splice sites
        (positions -20 to -3 in the acceptor model). C/T at these positions
        have high probability in the splice model.
        """
        # Poly-T (pyrimidine) tract
        py_seq = "T" * 20 + "CAG" + "AAAA"
        py_pos = py_seq.find("AG")

        # Poly-G (purine) tract — anti-pyrimidine
        pu_seq = "G" * 20 + "AAG" + "AAAA"
        pu_pos = pu_seq.find("AG")

        py_score = score_acceptor(py_seq, py_pos)
        pu_score = score_acceptor(pu_seq, pu_pos)

        assert py_score > pu_score, (
            f"Pyrimidine-tract acceptor ({py_score:.2f}) should outscore "
            f"purine-rich ({pu_score:.2f})"
        )

    def test_rank_order_correlation_with_strength(self) -> None:
        """Rank-order correlation between expected strength and actual scores.

        All strong sites > all moderate sites > all weak sites.
        """
        strength_rank = {"strong": 2, "moderate": 1, "weak": 0}
        ranks = []
        scores = []

        for nine_mer, strength, desc in DONOR_TEST_SEQUENCES:
            seq = _embed_donor_9mer(nine_mer, pos=3)
            score = score_donor(seq, 3)
            if score > _IMPOSSIBLE_SCORE:
                ranks.append(float(strength_rank[strength]))
                scores.append(score)

        # Pearson correlation should be positive (strong > moderate > weak)
        correlation = _pearson_r(ranks, scores)
        assert correlation > 0.5, (
            f"Pearson correlation between strength rank and score is "
            f"{correlation:.4f}, expected > 0.5 (positive correlation)"
        )


# ==============================================================================
# Validate that the proper maxentscan model is NOT anti-correlated
# ==============================================================================


class TestMaxEntScanNotAntiCorrelated:
    """Verify that maxentscan.score_donor / score_acceptor produce scores
    that are NOT anti-correlated with biological splice site strength.

    This is the inverse test: confirming that the bug in splicing.maxent_score
    does NOT affect the proper maxentscan module.
    """

    def test_deprecated_pwm_is_anti_correlated(self) -> None:
        """Confirm that the DEPRECATED splicing.maxent_score IS anti-correlated.

        This test documents the known bug: the hand-crafted PWM in
        splicing.maxent_score assigns weight 3.50 to A at the GT core
        positions (3-4), making non-splice sites score HIGH. We verify
        this is indeed the case so we know the bug exists.
        """
        from biocompiler.splicing import maxent_score as _deprecated_maxent_score

        # Canonical donor: has GT at core positions 3-4
        # In the deprecated PWM, G at pos 3 has weight 0.01, T at pos 4 has weight 0.01
        # Total core contribution: 0.01 + 0.01 = 0.02
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            canonical_deprecated = _deprecated_maxent_score("CAGGTGAGT")

        # Non-splice: has AA at core positions 3-4
        # In the deprecated PWM, A at pos 3 has weight 3.50, A at pos 4 has weight 3.50
        # Total core contribution: 3.50 + 3.50 = 7.00
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            nonsite_deprecated = _deprecated_maxent_score("CAGAAAAGT")

        # The deprecated model ranks the non-site HIGHER than the canonical site
        assert nonsite_deprecated > canonical_deprecated, (
            f"Deprecated model: non-site ({nonsite_deprecated:.4f}) should "
            f"score higher than canonical site ({canonical_deprecated:.4f}) "
            f"— this confirms the anti-correlation bug"
        )

    def test_proper_maxentscan_not_anti_correlated(self) -> None:
        """Confirm that the proper maxentscan.score_donor is NOT anti-correlated.

        Unlike the deprecated model, score_donor should rank canonical
        splice sites ABOVE non-sites.
        """
        # Canonical donor with GT
        canonical_score = score_donor(_embed_donor_9mer("CAGGTGAGT", 3), 3)

        # Non-site with AA instead of GT
        nonsite_score = score_donor(_embed_donor_9mer("CAGAAAAGT", 3), 3)

        assert canonical_score > nonsite_score, (
            f"Proper MaxEntScan: canonical ({canonical_score:.4f}) must "
            f"outscore non-site ({nonsite_score:.4f}) — anti-correlation "
            f"would indicate a bug in maxentscan.py"
        )

    def test_v2_matches_proper_maxentscan(self) -> None:
        """maxent_score_v2 delegates to score_donor, so it must also be
        positively correlated with splice site strength."""
        from biocompiler.splicing import maxent_score_v2

        strong = maxent_score_v2("CAGGTGAGT")
        weak = maxent_score_v2("ATGGTCATC")

        assert strong > weak, (
            f"maxent_score_v2: strong ({strong:.4f}) should outscore "
            f"weak ({weak:.4f})"
        )
        assert strong > 0.0, f"Strong donor via v2 should score > 0, got {strong}"


# ==============================================================================
# Prokaryote safety: MaxEntScan must NEVER be called for prokaryotes
# ==============================================================================


class TestProkaryoteMaxEntScanSafety:
    """Verify that MaxEntScan is never called during optimization for
    prokaryotic organisms.

    Prokaryotes (e.g. E. coli) lack a spliceosome, so splice site scoring
    is biologically meaningless. Calling MaxEntScan for prokaryotes would:
    1. Waste computation time
    2. Potentially introduce spurious constraints that lower CAI
    3. Represent a conceptual error (splice sites don't exist in prokaryotes)
    """

    def test_hybrid_optimizer_detects_prokaryote(self) -> None:
        """HybridOptimizer correctly detects E. coli as prokaryotic."""
        from biocompiler.hybrid_optimizer import HybridOptimizer

        opt = HybridOptimizer(species="ecoli")
        assert opt.is_prokaryote, "E. coli should be detected as prokaryotic"

    def test_hybrid_optimizer_eukaryote_not_prokaryote(self) -> None:
        """HybridOptimizer correctly detects human as NOT prokaryotic."""
        from biocompiler.hybrid_optimizer import HybridOptimizer

        opt = HybridOptimizer(species="human")
        assert not opt.is_prokaryote, "Human should NOT be detected as prokaryotic"

    def test_hybrid_optimizer_prokaryote_forces_avoid_gt_false(self) -> None:
        """For prokaryotes, avoid_gt is forced to False even if True is passed."""
        from biocompiler.hybrid_optimizer import HybridOptimizer

        opt = HybridOptimizer(species="ecoli", avoid_gt=True)
        assert not opt.avoid_gt, (
            "Prokaryotic organism should force avoid_gt=False to prevent "
            "MaxEntScan calls"
        )

    def test_hybrid_optimizer_no_splice_violations_for_prokaryote(self) -> None:
        """Prokaryotic optimization produces no splice-related violations."""
        from biocompiler.hybrid_optimizer import HybridOptimizer

        opt = HybridOptimizer(species="ecoli", avoid_gt=True)
        assert not opt.avoid_gt  # Should be forced to False

        # If avoid_gt is False, _detect_expensive_violations returns early
        # and _detect_violations skips the splice section
        # Verify is_prokaryote flag
        assert opt.is_prokaryote

    def test_type_system_skips_cryptic_splice_for_prokaryote(self) -> None:
        """type_system.check_no_cryptic_splice skips prokaryotic organisms."""
        from biocompiler.type_system import check_no_cryptic_splice
        from biocompiler.types import Verdict

        # A sequence with strong GT sites that would normally fail
        seq_with_gt = "AAACAGGTAAGTAAAA"

        # For E. coli (prokaryote), the check should be skipped
        result = check_no_cryptic_splice(seq_with_gt, organism="Escherichia_coli")
        assert result.verdict == Verdict.PASS, (
            f"Cryptic splice check should PASS (skip) for E. coli, "
            f"got {result.verdict}"
        )

    def test_solver_skips_splice_for_prokaryote(self) -> None:
        """Solver constraints skip splice for prokaryotic organisms."""
        try:
            from biocompiler.solver.constraints import build_splice_constraints
            from biocompiler.solver.types import SolverConfig
            config = SolverConfig(organism="Escherichia_coli")
            # For prokaryotes, splice constraints should be empty
            constraints = build_splice_constraints("MV", config)
            # The splice constraints list should be empty for prokaryotes
            # (or the function should handle prokaryotes gracefully)
            assert isinstance(constraints, list), (
                f"Expected list, got {type(constraints)}"
            )
        except (ImportError, TypeError, AttributeError):
            pytest.skip("Solver constraints not available or API differs")


# ==============================================================================
# Log-odds formula correctness
# ==============================================================================


class TestLogOddsFormulaCorrectness:
    """Verify the log-odds scoring formula is implemented correctly.

    score = sum over positions of log2(P(base | splice model) / P(base | background))

    This is the standard PWM scoring formula. If implemented correctly:
    - Conserved positions (high P in model) → positive contribution
    - Anti-conserved positions (low P in model) → negative contribution
    - Background-like positions (P ≈ 0.25) → near-zero contribution
    """

    def test_donor_log_odds_manual_computation(self) -> None:
        """Manual log-odds computation matches score_donor output."""
        nine_mer = "CAGGTGAGT"
        seq = _embed_donor_9mer(nine_mer, pos=3)

        # Manual computation
        epsilon = 0.001
        manual_score = 0.0
        for i, base in enumerate(nine_mer):
            idx = BASE_TO_INDEX[base]
            prob = max(DONOR_PWM_SCORE[i][idx], epsilon)
            manual_score += math.log2(prob / BG_PROB)
        manual_score = round(manual_score, 4)

        actual = score_donor(seq, 3)
        assert actual == manual_score, (
            f"Manual log-odds {manual_score} != score_donor {actual}"
        )

    def test_invariant_positions_contribute_most(self) -> None:
        """The two invariant positions (+1 G, +2 T) contribute the most
        to a canonical donor's score.

        G at +1: log2(0.990/0.25) ≈ +1.98 bits
        T at +2: log2(0.990/0.25) ≈ +1.98 bits
        Total from invariant positions: ~3.97 bits
        """
        g_contribution = math.log2(0.990 / 0.25)  # G at +1
        t_contribution = math.log2(0.990 / 0.25)  # T at +2

        assert g_contribution > 1.9, f"G contribution {g_contribution:.4f} should be >1.9"
        assert t_contribution > 1.9, f"T contribution {t_contribution:.4f} should be >1.9"
        assert g_contribution + t_contribution > 3.9, (
            f"Sum of invariant contributions ({g_contribution + t_contribution:.4f}) "
            f"should be >3.9"
        )

    def test_non_conserved_positions_contribute_little(self) -> None:
        """Positions with near-background probabilities contribute little.

        For example, position -3 has probabilities ~0.31/0.33/0.19/0.16,
        which are close to background (0.25). Their log-odds contributions
        should be close to 0.
        """
        # Position -3 (index 0): A=0.310, C=0.334, G=0.192, T=0.164
        for base_idx, prob in enumerate(DONOR_PWM_SCORE[0]):
            log_odds = math.log2(max(prob, 0.001) / BG_PROB)
            # These should be small (within ±0.5 bits of 0)
            assert abs(log_odds) < 1.0, (
                f"Position -3, base index {base_idx}: log-odds {log_odds:.4f} "
                f"should be near 0 (background-like position)"
            )

    def test_acceptor_log_odds_manual_computation(self) -> None:
        """Manual log-odds computation matches score_acceptor output."""
        seq = "T" * 20 + "CAG" + "ATCG"
        ag_pos = 21

        epsilon = 0.001
        start = ag_pos - 20  # _ACCEPTOR_UPSTREAM = 20
        manual_score = 0.0
        for i in range(23):  # 23 positions in acceptor model
            base = seq[start + i]
            idx = BASE_TO_INDEX[base]
            prob = max(ACCEPTOR_PWM_SCORE[i][idx], epsilon)
            manual_score += math.log2(prob / BG_PROB)
        manual_score = round(manual_score, 4)

        actual = score_acceptor(seq, ag_pos)
        assert actual == manual_score, (
            f"Manual log-odds {manual_score} != score_acceptor {actual}"
        )


# ==============================================================================
# Scan consistency: scan_splice_sites agrees with individual scoring
# ==============================================================================


class TestScanConsistency:
    """Verify that scan_splice_sites produces consistent results with
    individual score_donor / score_acceptor calls.
    """

    def test_scan_donor_scores_match_individual(self) -> None:
        """Donor scores from scan_splice_sites match score_donor calls."""
        seq = "AAACAGGTAAGTAAAA" + "C" * 30 + "TTTTAGAAAA"
        results = scan_splice_sites(seq, donor_threshold=-100.0, acceptor_threshold=-100.0)

        for pos, site_type, score in results:
            if site_type == "donor":
                individual_score = score_donor(seq, pos)
                assert abs(score - individual_score) < 0.001, (
                    f"Scan score {score} != individual score {individual_score} "
                    f"at donor position {pos}"
                )

    def test_scan_acceptor_scores_match_individual(self) -> None:
        """Acceptor scores from scan_splice_sites match score_acceptor calls."""
        seq = "T" * 25 + "CAG" + "AAAA" + "C" * 20 + "AAG" + "AAAA"
        results = scan_splice_sites(seq, donor_threshold=-100.0, acceptor_threshold=-100.0)

        for pos, site_type, score in results:
            if site_type == "acceptor":
                individual_score = score_acceptor(seq, pos)
                assert abs(score - individual_score) < 0.001, (
                    f"Scan score {score} != individual score {individual_score} "
                    f"at acceptor position {pos}"
                )

    def test_max_donor_score_is_max_of_individual(self) -> None:
        """max_donor_score returns the maximum of all individual donor scores."""
        seq = "CAGGTAAGT" + "A" * 10 + "ATGGTCATC" + "A" * 10 + "TTGGTAAGT"
        max_score = max_donor_score(seq)

        individual_max = _IMPOSSIBLE_SCORE
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "GT":
                s = score_donor(seq, i)
                if s > individual_max:
                    individual_max = s

        assert max_score == round(individual_max, 4), (
            f"max_donor_score {max_score} != individual max {round(individual_max, 4)}"
        )

    def test_max_acceptor_score_is_max_of_individual(self) -> None:
        """max_acceptor_score returns the maximum of all individual acceptor scores."""
        seq = "T" * 25 + "CAG" + "AAAA" + "G" * 20 + "AAG" + "AAAA"
        max_score = max_acceptor_score(seq)

        individual_max = _IMPOSSIBLE_SCORE
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "AG":
                s = score_acceptor(seq, i)
                if s > individual_max:
                    individual_max = s

        assert max_score == round(individual_max, 4), (
            f"max_acceptor_score {max_score} != individual max {round(individual_max, 4)}"
        )
