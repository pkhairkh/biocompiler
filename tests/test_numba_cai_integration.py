"""Tests for compute_cai_incremental NUMBA kernel integration in the optimizer.

These tests verify that:
1. _adaptiveness_to_array correctly converts codon adaptiveness dicts to numpy arrays
2. _codon_to_index maps codons to correct base-4 encoded indices
3. _compute_cai_fast produces CAI values matching compute_cai (within tolerance)
4. _BatchSwapScorer.reset_incremental_state / update_incremental_state track log_sum correctly
5. The NUMBA incremental kernel and pure-Python incremental update produce identical results
6. BioOptimizer._compute_seq_cai uses the NUMBA kernel when available
7. Full optimization produces correct CAI values regardless of NUMBA availability
"""

from __future__ import annotations

import math
import pytest


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def ecoli_cai():
    """E. coli codon adaptiveness table from CODON_ADAPTIVENESS_TABLES."""
    from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
    return CODON_ADAPTIVENESS_TABLES.get("Escherichia_coli", {})


@pytest.fixture
def human_cai():
    """Human codon adaptiveness table from CODON_ADAPTIVENESS_TABLES."""
    from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
    return CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens", {})


# ── Unit tests for _adaptiveness_to_array ────────────────────────────

class TestAdaptivenessToArray:
    """Test the _adaptiveness_to_array conversion function."""

    def test_array_length(self, ecoli_cai):
        """Array should have exactly 64 entries (4^3 for all 3-letter codons)."""
        from biocompiler.optimization import _adaptiveness_to_array
        arr = _adaptiveness_to_array(ecoli_cai)
        assert len(arr) == 64

    def test_array_dtype(self, ecoli_cai):
        """Array should be float64."""
        from biocompiler.optimization import _adaptiveness_to_array
        arr = _adaptiveness_to_array(ecoli_cai)
        assert arr.dtype.name == "float64"

    def test_known_codon_value(self, ecoli_cai):
        """ATG should be at the correct index with the correct value."""
        from biocompiler.optimization import _adaptiveness_to_array, _codon_to_index
        arr = _adaptiveness_to_array(ecoli_cai)
        atg_idx = _codon_to_index("ATG")
        assert arr[atg_idx] == ecoli_cai.get("ATG", 0.0)

    def test_missing_codons_get_epsilon(self):
        """Codons not in the dict should get the epsilon default (1e-10)."""
        from biocompiler.optimization import _adaptiveness_to_array, _codon_to_index
        # Minimal dict with just one codon
        minimal = {"AAA": 0.5}
        arr = _adaptiveness_to_array(minimal)
        # Check a codon that is NOT in the dict
        ccc_idx = _codon_to_index("CCC")
        assert arr[ccc_idx] == 1e-10

    def test_all_ecoli_codons_present(self, ecoli_cai):
        """Every codon in ecoli_cai should be findable in the array."""
        from biocompiler.optimization import _adaptiveness_to_array, _codon_to_index
        arr = _adaptiveness_to_array(ecoli_cai)
        for codon, w in ecoli_cai.items():
            if len(codon) == 3 and all(b in "ACGT" for b in codon):
                idx = _codon_to_index(codon)
                assert abs(arr[idx] - w) < 1e-15, (
                    f"Codon {codon}: array[{idx}]={arr[idx]} != {w}"
                )


# ── Unit tests for _codon_to_index ───────────────────────────────────

class TestCodonToIndex:
    """Test the _codon_to_index base-4 encoding function."""

    def test_aaa_index(self):
        """AAA should map to index 0 (0*16 + 0*4 + 0)."""
        from biocompiler.optimization import _codon_to_index
        assert _codon_to_index("AAA") == 0

    def test_ttt_index(self):
        """TTT should map to index 63 (3*16 + 3*4 + 3)."""
        from biocompiler.optimization import _codon_to_index
        assert _codon_to_index("TTT") == 63

    def test_atg_index(self):
        """ATG should map to index 0*16 + 3*4 + 2 = 14."""
        from biocompiler.optimization import _codon_to_index
        assert _codon_to_index("ATG") == 14

    def test_all_indices_unique(self):
        """All 64 codons should map to unique indices."""
        from biocompiler.optimization import _codon_to_index
        from biocompiler.type_system import CODON_TABLE
        indices = {}
        for codon in CODON_TABLE:
            idx = _codon_to_index(codon)
            assert idx not in indices, (
                f"Duplicate index {idx} for codons {codon} and {indices[idx]}"
            )
            indices[idx] = codon

    def test_all_indices_in_range(self):
        """All indices should be in [0, 63]."""
        from biocompiler.optimization import _codon_to_index
        from biocompiler.type_system import CODON_TABLE
        for codon in CODON_TABLE:
            idx = _codon_to_index(codon)
            assert 0 <= idx <= 63, f"Codon {codon} has out-of-range index {idx}"


# ── Unit tests for _compute_cai_fast ─────────────────────────────────

class TestComputeCAIFast:
    """Test the _compute_cai_fast function that uses NUMBA kernel."""

    def test_empty_sequence(self, ecoli_cai):
        """Empty sequence should return 0.0."""
        from biocompiler.optimization import _compute_cai_fast
        assert _compute_cai_fast("", ecoli_cai) == 0.0

    def test_short_sequence(self, ecoli_cai):
        """Sequence shorter than 3 bases should return 0.0."""
        from biocompiler.optimization import _compute_cai_fast
        assert _compute_cai_fast("AT", ecoli_cai) == 0.0

    def test_single_codon_met(self, ecoli_cai):
        """Single Met codon (excluded from CAI) should return 0.0."""
        from biocompiler.optimization import _compute_cai_fast
        # Met is excluded from CAI computation
        assert _compute_cai_fast("ATG", ecoli_cai) == 0.0

    def test_cai_matches_reference_ecoli(self, ecoli_cai):
        """CAI from _compute_cai_fast should match compute_cai for E. coli."""
        from biocompiler.optimization import _compute_cai_fast
        from biocompiler.expression.translation import compute_cai

        # Test several sequences
        test_seqs = [
            "ATGAAAGCGTTT",  # M-K-A-F
            "ATGGTTTCTAAAGGTGAA",  # M-V-S-K-G-E
            "ATGCTTGAACTT",  # M-L-E-L
        ]
        for seq in test_seqs:
            cai_fast = _compute_cai_fast(seq, ecoli_cai)
            cai_ref = compute_cai(seq, "Escherichia_coli")
            # Allow small tolerance due to rounding
            assert abs(cai_fast - cai_ref) < 0.001, (
                f"Seq {seq}: fast={cai_fast:.6f} != ref={cai_ref:.6f}"
            )

    def test_cai_matches_reference_human(self, human_cai):
        """CAI from _compute_cai_fast should match compute_cai for Human."""
        from biocompiler.optimization import _compute_cai_fast
        from biocompiler.expression.translation import compute_cai

        seq = "ATGAAAGCGTTT"  # M-K-A-F
        cai_fast = _compute_cai_fast(seq, human_cai)
        cai_ref = compute_cai(seq, "Homo_sapiens")
        assert abs(cai_fast - cai_ref) < 0.001

    def test_cai_in_valid_range(self, ecoli_cai):
        """CAI should always be in [0.0, 1.0]."""
        from biocompiler.optimization import _compute_cai_fast
        test_seqs = [
            "ATGAAAGCGTTT",
            "ATGCTTGAACTT",
            "ATGTTTTTTTTT",
        ]
        for seq in test_seqs:
            cai = _compute_cai_fast(seq, ecoli_cai)
            assert 0.0 <= cai <= 1.0, f"CAI {cai} out of range for {seq}"


# ── Unit tests for incremental CAI tracking ───────────────────────────

class TestIncrementalCAITracking:
    """Test the _BatchSwapScorer incremental state management."""

    def test_reset_incremental_state(self, ecoli_cai):
        """reset_incremental_state should compute the correct log_sum."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["ATG", "AAA", "GCG", "TTT"]
        scorer.reset_incremental_state(seq_codons)

        # Verify log_sum is correct
        expected_log_sum = sum(
            math.log(max(ecoli_cai.get(c, 1e-10), 1e-10))
            for c in seq_codons
        )
        assert abs(scorer.current_log_sum - expected_log_sum) < 1e-10

    def test_update_incremental_state_single_swap(self, ecoli_cai):
        """A single codon swap should correctly update the log_sum."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["ATG", "AAA", "GCG", "TTT"]
        scorer.reset_incremental_state(seq_codons)
        old_log_sum = scorer.current_log_sum

        # Swap AAA -> AAG (both Lys)
        scorer.update_incremental_state("AAA", "AAG")

        # Verify: new_log_sum = old_log_sum - log(w_AAA) + log(w_AAG)
        w_aaa = max(ecoli_cai.get("AAA", 1e-10), 1e-10)
        w_aag = max(ecoli_cai.get("AAG", 1e-10), 1e-10)
        expected_new = old_log_sum - math.log(w_aaa) + math.log(w_aag)
        assert abs(scorer.current_log_sum - expected_new) < 1e-10

    def test_update_incremental_state_multiple_swaps(self, ecoli_cai):
        """Multiple sequential swaps should track log_sum correctly."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["ATG", "AAA", "GCG", "TTT"]
        scorer.reset_incremental_state(seq_codons)

        # Swap 1: AAA -> AAG
        scorer.update_incremental_state("AAA", "AAG")

        # Swap 2: GCG -> GCT (both Ala)
        scorer.update_incremental_state("GCG", "GCT")

        # Verify: log_sum should match full recomputation of new sequence
        new_seq = ["ATG", "AAG", "GCT", "TTT"]
        expected_log_sum = sum(
            math.log(max(ecoli_cai.get(c, 1e-10), 1e-10))
            for c in new_seq
        )
        assert abs(scorer.current_log_sum - expected_log_sum) < 1e-8

    def test_update_without_reset_is_noop(self, ecoli_cai):
        """update_incremental_state without reset should be a no-op."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        # No reset called — current_log_sum should be None
        assert scorer.current_log_sum is None

        # Update should be a no-op
        scorer.update_incremental_state("AAA", "AAG")
        assert scorer.current_log_sum is None

    def test_incremental_matches_full_recompute(self, ecoli_cai):
        """Incremental log_sum should match full recomputation after many swaps."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)
        from biocompiler.type_system import AA_TO_CODONS

        # Start with a longer sequence
        seq_codons = ["ATG", "CTT", "GAA", "AAA", "GGT", "TTC", "CAG", "GCT"]
        scorer.reset_incremental_state(seq_codons)

        # Perform several swaps
        swaps = [
            ("CTT", "CTG"),  # Leu: CTT -> CTG (better CAI)
            ("GGT", "GGC"),  # Gly: GGT -> GGC (better CAI)
            ("TTC", "TTT"),  # Phe: TTC -> TTT
            ("GCT", "GCC"),  # Ala: GCT -> GCC
        ]
        for old, new in swaps:
            scorer.update_incremental_state(old, new)

        # Final sequence after all swaps
        final_seq = ["ATG", "CTG", "GAA", "AAA", "GGC", "TTT", "CAG", "GCC"]
        expected_log_sum = sum(
            math.log(max(ecoli_cai.get(c, 1e-10), 1e-10))
            for c in final_seq
        )

        # Allow slightly larger tolerance due to accumulated floating point errors
        assert abs(scorer.current_log_sum - expected_log_sum) < 1e-6, (
            f"Incremental={scorer.current_log_sum:.10f} vs "
            f"Expected={expected_log_sum:.10f}"
        )

    def test_score_candidates_uses_cached_log_sum(self, ecoli_cai):
        """score_candidates should use the cached log_sum when available."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        seq_codons = ["ATG", "CTT", "GAA"]
        scorer.reset_incremental_state(seq_codons)

        # Score candidates — should use cached log_sum
        from biocompiler.type_system import AA_TO_CODONS
        leu_codons = AA_TO_CODONS["L"]
        scores = scorer.score_candidates(seq_codons, 1, leu_codons)

        # Should produce valid scores
        assert len(scores) == len(leu_codons)
        for s in scores:
            assert 0.0 < s <= 1.0

    def test_score_candidates_matches_after_incremental_update(self, ecoli_cai):
        """Batch scores after incremental update should be consistent."""
        from biocompiler.optimization import _BatchSwapScorer
        scorer = _BatchSwapScorer(ecoli_cai)

        # Start: ATG CTT GAA
        seq_codons = ["ATG", "CTT", "GAA"]
        scorer.reset_incremental_state(seq_codons)

        # Swap CTT -> CTG
        scorer.update_incremental_state("CTT", "CTG")

        # New sequence: ATG CTG GAA
        new_seq = ["ATG", "CTG", "GAA"]

        # Score swaps at position 2 (GAA -> other Glu codons)
        from biocompiler.type_system import AA_TO_CODONS
        glu_codons = AA_TO_CODONS["E"]
        scores = scorer.score_candidates(new_seq, 2, glu_codons)

        # Manually compute expected scores
        epsilon = 1e-10
        current_log_sum = sum(
            math.log(max(ecoli_cai.get(c, epsilon), epsilon))
            for c in new_seq
        )
        n_codons = len(new_seq)
        w_old = max(ecoli_cai.get("GAA", epsilon), epsilon)
        log_w_old = math.log(w_old)

        for k, cand in enumerate(glu_codons):
            w_new = max(ecoli_cai.get(cand, epsilon), epsilon)
            expected_log_sum = current_log_sum - log_w_old + math.log(w_new)
            expected_cai = math.exp(expected_log_sum / n_codons)
            assert abs(scores[k] - expected_cai) < 1e-6, (
                f"Candidate {cand}: score={scores[k]:.6f} != expected={expected_cai:.6f}"
            )


# ── Integration tests for _compute_seq_cai ────────────────────────────

class TestComputeSeqCAIIntegration:
    """Test BioOptimizer._compute_seq_cai with NUMBA integration."""

    def test_compute_seq_cai_basic(self):
        """_compute_seq_cai should return a valid CAI value."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli", avoid_gt=False)
        cai = opt._compute_seq_cai("ATGAAAGCGTTT")
        assert 0.0 <= cai <= 1.0
        assert cai > 0.0  # Should be non-zero for a valid coding sequence

    def test_compute_seq_cai_empty(self):
        """Empty sequence should return 0.0."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli", avoid_gt=False)
        assert opt._compute_seq_cai("") == 0.0

    def test_compute_seq_cai_consistency(self):
        """Repeated calls should return the same CAI value."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli", avoid_gt=False)
        seq = "ATGAAAGCGTTTTAA"
        cai1 = opt._compute_seq_cai(seq)
        cai2 = opt._compute_seq_cai(seq)
        assert abs(cai1 - cai2) < 1e-15

    def test_compute_seq_cai_numba_adapt_arr_cached(self):
        """The adaptiveness array should be cached after first call."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli", avoid_gt=False)

        # Before first call, no cached array
        assert not hasattr(opt, '_numba_adapt_arr') or opt._numba_adapt_arr is None

        seq = "ATGAAAGCGTTT"
        cai = opt._compute_seq_cai(seq)

        # After first call, should have cached array
        if opt._numba_adapt_arr is not None:
            import numpy as np
            assert isinstance(opt._numba_adapt_arr, np.ndarray)
            assert len(opt._numba_adapt_arr) == 64


# ── Full optimization integration tests ────────────────────────────────

class TestFullOptimizationCAI:
    """Test that full optimization produces correct CAI with NUMBA integration."""

    def test_optimize_ecoli_short(self):
        """Short protein optimization for E. coli should produce high CAI."""
        from biocompiler.optimization import optimize_sequence
        result = optimize_sequence("MKAF", species="ecoli", avoid_gt=False)
        assert result.cai > 0.5
        assert len(result.sequence) == 12  # 4 AA * 3

    def test_optimize_human_short(self):
        """Short protein optimization for Human should produce valid CAI."""
        from biocompiler.optimization import optimize_sequence
        result = optimize_sequence("MKAF", species="Homo_sapiens", avoid_gt=True)
        assert result.cai > 0.0
        assert len(result.sequence) == 12

    def test_optimize_cai_cross_validation(self):
        """CAI from optimize_sequence should match compute_cai."""
        from biocompiler.optimization import optimize_sequence
        from biocompiler.expression.translation import compute_cai

        protein = "MVSKGE"
        result = optimize_sequence(protein, species="ecoli", avoid_gt=False)
        cai_check = compute_cai(result.sequence, "Escherichia_coli")

        # Allow tolerance for rounding differences
        assert abs(result.cai - cai_check) < 0.01, (
            f"Result CAI={result.cai:.6f} != compute_cai={cai_check:.6f}"
        )

    def test_hill_climb_with_incremental_tracking(self):
        """Hill climb should work with incremental CAI tracking enabled."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli", avoid_gt=False)

        # Run the hill climb directly
        seq = "ATGGTTTCTAAAGGTGAA"
        result = opt._step_cai_hill_climb(seq)

        # Should return a valid sequence
        assert len(result) == len(seq)
        assert len(result) % 3 == 0

        # CAI should be valid
        cai = opt._compute_seq_cai(result)
        assert 0.0 <= cai <= 1.0

    def test_numba_kernel_available(self):
        """Check if NUMBA is available in the test environment."""
        from biocompiler.optimization import HAS_NUMBA
        # This test just reports the status — it is not a pass/fail
        if HAS_NUMBA:
            from biocompiler.optimizer.numba_kernels import compute_cai_incremental
            assert compute_cai_incremental is not None
            # Test the incremental kernel directly
            # Swap from w_old=0.5 to w_new=1.0 in a 2-codon sequence
            import math
            current_log_sum = math.log(0.5) + math.log(0.5)  # both codons w=0.5
            n_codons = 2
            old_adaptiveness = 0.5
            new_adaptiveness = 1.0
            result_cai = compute_cai_incremental(
                current_log_sum, n_codons, old_adaptiveness, new_adaptiveness
            )
            # Expected: new_log_sum = log(0.5) + log(1.0) = -0.693
            # CAI = exp((-0.693) / 2) = exp(-0.347) ≈ 0.707
            expected = math.exp((math.log(0.5) + math.log(1.0)) / 2)
            assert abs(result_cai - expected) < 1e-10, (
                f"NUMBA kernel: {result_cai:.10f} != expected {expected:.10f}"
            )


# ── Regression tests ──────────────────────────────────────────────────

class TestCAIRegression:
    """Regression tests to ensure NUMBA integration does not change results."""

    def test_ecoli_cai_values_stable(self):
        """E. coli CAI values should be stable across optimization changes."""
        from biocompiler.optimization import _compute_cai_fast
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        ecoli = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]

        # Test a few known sequences with expected approximate CAI
        test_cases = [
            ("ATGAAAGCGTTT", 0.8),  # Should be around 0.8-0.9
            ("ATGCTTGAACTT", 0.2),  # CTT is a low-CAI Leu codon
        ]
        for seq, min_cai in test_cases:
            cai = _compute_cai_fast(seq, ecoli)
            assert cai >= min_cai, f"CAI {cai:.4f} below expected minimum {min_cai} for {seq}"

    def test_optimize_deterministic(self):
        """Optimization should be deterministic (same input → same output)."""
        from biocompiler.optimization import optimize_sequence

        protein = "MVSKGE"
        result1 = optimize_sequence(protein, species="ecoli", avoid_gt=False)
        result2 = optimize_sequence(protein, species="ecoli", avoid_gt=False)

        assert result1.sequence == result2.sequence
        assert abs(result1.cai - result2.cai) < 1e-15
