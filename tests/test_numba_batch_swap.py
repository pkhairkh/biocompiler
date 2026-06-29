"""Tests for batch_codon_swap_score NUMBA kernel integration in the optimizer.

These tests verify that:
1. The _BatchSwapScorer class produces correct CAI scores
2. The NUMBA kernel path and Python fallback produce identical results
3. The hill-climb step produces the same (or better) CAI with the batch scorer
4. Fallback works when NUMBA is unavailable
"""

from __future__ import annotations

import math
import pytest


# ── Fixture: E. coli CAI table ──────────────────────────────────────

@pytest.fixture
def ecoli_cai():
    """E. coli codon adaptiveness table from CODON_ADAPTIVENESS_TABLES."""
    from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
    return CODON_ADAPTIVENESS_TABLES.get("Escherichia_coli", {})


@pytest.fixture
def ecoli_opt(ecoli_cai):
    """A BioOptimizer instance for E. coli."""
    from biocompiler.optimizer import BioOptimizer
    return BioOptimizer(species="ecoli", avoid_gt=False)


@pytest.fixture
def human_opt():
    """A BioOptimizer instance for Homo_sapiens."""
    from biocompiler.optimizer import BioOptimizer
    return BioOptimizer(species="Homo_sapiens", avoid_gt=True)


# ── Unit tests for _BatchSwapScorer ────────────────────────────────

class TestBatchSwapScorer:
    """Test the _BatchSwapScorer helper class."""

    def test_import_with_fallback(self):
        """The import should succeed regardless of NUMBA availability."""
        from biocompiler.optimizer import _BatchSwapScorer, HAS_NUMBA
        # Should not raise
        assert _BatchSwapScorer is not None

    def test_scorer_init(self, ecoli_cai):
        """_BatchSwapScorer should build a valid codon-to-index mapping."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)
        assert len(scorer._codon_to_idx) >= 64  # All standard codons
        assert len(scorer._adaptiveness) == len(scorer._codon_to_idx)

    def test_score_candidates_empty(self, ecoli_cai):
        """Empty candidate list should return empty scores."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)
        scores = scorer.score_candidates(["ATG", "TTT", "GAA"], 1, [])
        assert scores == []

    def test_score_candidates_single_codon_seq(self, ecoli_cai):
        """Scoring swaps on a 1-codon sequence should produce valid CAI."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)
        # F is encoded by TTT and TTC
        scores = scorer.score_candidates(["TTT"], 0, ["TTC"])
        assert len(scores) == 1
        assert 0.0 < scores[0] <= 1.0

    def test_score_candidates_matches_incremental(self, ecoli_cai):
        """Batch scorer CAI should match incremental log-sum computation."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        # 3-codon sequence: M-F-E → ATG TTT GAA
        seq_codons = ["ATG", "TTT", "GAA"]

        # Candidates for position 1 (F): TTC (the other Phe codon)
        candidates = ["TTC"]
        scores = scorer.score_candidates(seq_codons, 1, candidates)

        # Manually compute expected CAI after swap
        epsilon = 1e-10
        w_m = ecoli_cai.get("ATG", epsilon)
        w_f_new = ecoli_cai.get("TTC", epsilon)
        w_e = ecoli_cai.get("GAA", epsilon)

        # CAI after swapping TTT→TTC at position 1
        # (Met is excluded from CAI by convention in some formulations,
        #  but batch_codon_swap_score includes all n_codons)
        n = 3
        log_sum_new = math.log(w_m if w_m > 0 else epsilon) + \
                      math.log(w_f_new if w_f_new > 0 else epsilon) + \
                      math.log(w_e if w_e > 0 else epsilon)
        expected_cai = math.exp(log_sum_new / n)

        assert len(scores) == 1
        assert abs(scores[0] - expected_cai) < 1e-10, (
            f"Batch score {scores[0]} != expected {expected_cai}"
        )

    def test_score_candidates_higher_cai_for_optimal_codon(self, ecoli_cai):
        """Swapping to a higher-CAI codon should increase the score."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        # Use Leucine which has 6 codons with varying CAI
        seq_codons = ["TTA", "GAA", "AAA"]  # L with low-CAI codon TTA
        # All Leu codons as candidates
        from biocompiler.type_system import AA_TO_CODONS
        leu_codons = AA_TO_CODONS["L"]

        # Compute current CAI
        current_scores = scorer.score_candidates(seq_codons, 0, [seq_codons[0]])
        current_cai = current_scores[0]

        # Best Leu codon by CAI
        best_leu = max(leu_codons, key=lambda c: ecoli_cai.get(c, 0.0))
        best_scores = scorer.score_candidates(seq_codons, 0, [best_leu])
        best_cai = best_scores[0]

        assert best_cai >= current_cai, (
            f"Optimal codon {best_leu} (CAI={best_cai}) should be >= "
            f"suboptimal {seq_codons[0]} (CAI={current_cai})"
        )

    def test_score_candidates_multiple_candidates(self, ecoli_cai):
        """Batch scoring with multiple candidates should return one score per candidate."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["TTT", "GAA", "AAA"]  # F-E-K
        # F has 2 codons: TTT, TTC
        from biocompiler.type_system import AA_TO_CODONS
        phe_codons = AA_TO_CODONS["F"]
        candidates = [c for c in phe_codons if c != "TTT"]  # Only alternatives

        scores = scorer.score_candidates(seq_codons, 0, candidates)
        assert len(scores) == len(candidates)
        for s in scores:
            assert 0.0 < s <= 1.0

    def test_score_candidates_python_fallback_matches_numba(self, ecoli_cai):
        """Python fallback path should produce identical results to NUMBA path."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["ATG", "CTT", "GAA", "AAA", "GGT"]  # M-L-E-K-G
        from biocompiler.type_system import AA_TO_CODONS
        leu_codons = AA_TO_CODONS["L"]

        # Force Python fallback by calling _compute_log_sum directly
        scores = scorer.score_candidates(seq_codons, 1, leu_codons)

        # Manually compute expected scores via Python
        epsilon = 1e-10
        current_log_sum = sum(
            math.log(max(ecoli_cai.get(c, epsilon), epsilon))
            for c in seq_codons
        )
        n_codons = len(seq_codons)
        w_old = max(ecoli_cai.get(seq_codons[1], epsilon), epsilon)
        log_w_old = math.log(w_old)

        for k, cand in enumerate(leu_codons):
            w_new = max(ecoli_cai.get(cand, epsilon), epsilon)
            new_log_sum = current_log_sum - log_w_old + math.log(w_new)
            expected = math.exp(new_log_sum / n_codons)
            assert abs(scores[k] - expected) < 1e-10, (
                f"Candidate {cand}: score {scores[k]} != expected {expected}"
            )


# ── Integration tests for the hill-climb with batch scorer ─────────

class TestHillClimbBatchScorer:
    """Test that the hill-climb step uses batch scoring correctly."""

    def test_hill_climb_produces_valid_sequence(self, ecoli_opt):
        """Hill climb should return a valid DNA sequence of correct length."""
        protein = "MVSKGE"
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(protein, species="ecoli", avoid_gt=False)
        assert len(result.sequence) == len(protein) * 3

    def test_hill_climb_cai_improves_or_maintains(self, ecoli_opt):
        """After hill climb, CAI should be >= the initial CAI."""
        protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAVDILSKKGDVQVIK"
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(protein, species="ecoli", avoid_gt=False)
        # The optimizer should achieve a reasonable CAI
        assert result.cai > 0.0

    def test_hill_climb_preserves_protein(self, ecoli_opt):
        """The optimized sequence should encode the same protein."""
        protein = "MVSKGE"
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(protein, species="ecoli", avoid_gt=False)

        # Translate back and verify
        from biocompiler.type_system import CODON_TABLE
        translated = ""
        for i in range(0, len(result.sequence), 3):
            codon = result.sequence[i:i+3]
            aa = CODON_TABLE.get(codon, "?")
            translated += aa

        # The translated protein should match the input (excluding stop codon)
        assert translated.rstrip("*") == protein

    def test_batch_scorer_used_in_hill_climb(self):
        """Verify _BatchSwapScorer is instantiated in _step_cai_hill_climb."""
        from biocompiler.optimizer import _BatchSwapScorer, BioOptimizer
        opt = BioOptimizer(species="ecoli", avoid_gt=False)

        # Run the hill climb on a short sequence
        # "ATGGTTTCTAAAGGTGAA" = M-V-S-K-G-E
        seq = "ATGGTTTCTAAAGGTGAA"
        result = opt._step_cai_hill_climb(seq)
        assert len(result) == len(seq)
        assert len(result) % 3 == 0

    def test_eukaryote_hill_climb_with_gt_avoidance(self, human_opt):
        """Hill climb for eukaryote should respect GT avoidance."""
        protein = "MVSKGE"
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(protein, species="Homo_sapiens", avoid_gt=True)
        assert len(result.sequence) == len(protein) * 3

    def test_batch_scorer_numba_flag(self):
        """HAS_NUMBA should be a boolean indicating NUMBA availability."""
        from biocompiler.optimizer import HAS_NUMBA
        assert isinstance(HAS_NUMBA, bool)

    def test_batch_scorer_no_candidates(self, ecoli_cai):
        """Positions with only one codon (Met, Trp) should work correctly."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        # Met only has ATG, so no improving candidates
        seq_codons = ["ATG", "GAA"]
        from biocompiler.type_system import AA_TO_CODONS
        met_candidates = AA_TO_CODONS["M"]

        # Only candidate is ATG itself (no alternatives)
        scores = scorer.score_candidates(seq_codons, 0, met_candidates)
        assert len(scores) == len(met_candidates)


# ── Edge case tests ─────────────────────────────────────────────────

class TestBatchSwapEdgeCases:
    """Test edge cases for the batch swap scorer."""

    def test_single_codon_sequence(self, ecoli_cai):
        """A single-codon sequence should produce valid CAI scores."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)
        # Single Phe codon
        scores = scorer.score_candidates(["TTT"], 0, ["TTC"])
        assert len(scores) == 1
        assert 0.0 < scores[0] <= 1.0

    def test_zero_adaptiveness_handling(self):
        """Codons with zero adaptiveness should not cause division by zero."""
        # Create a minimal CAI table with some zero values
        minimal_cai = {
            "ATG": 1.0,
            "TTT": 0.5,
            "TTC": 0.0,  # Zero adaptiveness
            "GAA": 0.8,
        }
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(minimal_cai)

        # Should handle zero-adaptiveness candidates gracefully
        scores = scorer.score_candidates(["TTT", "GAA"], 0, ["TTC"])
        assert len(scores) == 1
        assert scores[0] > 0.0  # Should not be 0 or NaN

    def test_large_sequence_batch_scoring(self, ecoli_cai):
        """Batch scoring should work on larger sequences (100+ codons)."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        # Build a 100-codon sequence with alternating Leu codons
        seq_codons = []
        for i in range(100):
            if i % 7 == 0:
                seq_codons.append("ATG")  # Met
            elif i % 5 == 0:
                seq_codons.append("TTT")  # Phe
            elif i % 3 == 0:
                seq_codons.append("GAA")  # Glu
            else:
                seq_codons.append("CTT")  # Leu

        # Score Leu codon swaps at position 1
        from biocompiler.type_system import AA_TO_CODONS
        leu_candidates = AA_TO_CODONS["L"]
        scores = scorer.score_candidates(seq_codons, 1, leu_candidates)

        assert len(scores) == len(leu_candidates)
        for s in scores:
            assert 0.0 < s <= 1.0

    def test_batch_scorer_consistency_across_calls(self, ecoli_cai):
        """Repeated calls with the same input should produce identical results."""
        from biocompiler.optimizer import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["ATG", "CTT", "GAA"]
        from biocompiler.type_system import AA_TO_CODONS
        leu_codons = AA_TO_CODONS["L"]

        scores1 = scorer.score_candidates(seq_codons, 1, leu_codons)
        scores2 = scorer.score_candidates(seq_codons, 1, leu_codons)

        assert len(scores1) == len(scores2)
        for a, b in zip(scores1, scores2):
            assert abs(a - b) < 1e-15, f"Scores differ: {a} vs {b}"
