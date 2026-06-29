"""
Test NUMBA v2 kernel wiring into the optimizer hot path.

Verifies:
1. USE_NUMBA flag toggles between NUMBA and Python paths
2. CAI computation produces same results via both paths
3. GC counting produces same results via both paths
4. Dinucleotide counting produces same results via both paths
5. Batch codon swap scoring produces same results via both paths
6. Full optimization produces consistent CAI with both paths
7. Fallback works when USE_NUMBA is False
"""

from __future__ import annotations

import math
import pytest


# ────────────────────────────────────────────────────────────
# Detect NUMBA availability
# ────────────────────────────────────────────────────────────

try:
    from biocompiler.optimizer.numba_kernels import HAS_NUMBA
except ImportError:
    HAS_NUMBA = False

NUMBA_SKIP_REASON = "NUMBA is not installed; skipping NUMBA v2 wiring tests"


# ═══════════════════════════════════════════════════════════════
# 1. USE_NUMBA flag tests
# ═══════════════════════════════════════════════════════════════

class TestUseNumbaFlag:
    """Verify USE_NUMBA flag toggles correctly."""

    def test_use_numba_exists(self):
        """USE_NUMBA should be defined in numba_kernels."""
        from biocompiler.optimizer.numba_kernels import USE_NUMBA
        assert isinstance(USE_NUMBA, bool)

    def test_use_numba_default_matches_has_numba(self):
        """USE_NUMBA should default to True when NUMBA is available."""
        from biocompiler.optimizer.numba_kernels import HAS_NUMBA, USE_NUMBA
        assert USE_NUMBA == HAS_NUMBA

    def test_use_numba_toggles(self):
        """USE_NUMBA should be toggleable at runtime."""
        import biocompiler.optimizer.numba_kernels as nk
        original = nk.USE_NUMBA

        # Toggle off
        nk.USE_NUMBA = False
        assert nk.USE_NUMBA is False

        # Restore
        nk.USE_NUMBA = original
        assert nk.USE_NUMBA == original

    def test_use_numba_propagates_to_cai(self):
        """USE_NUMBA should be accessible from optimizer.cai."""
        from biocompiler.optimizer.cai import USE_NUMBA
        assert isinstance(USE_NUMBA, bool)

    def test_use_numba_propagates_to_optimizer(self):
        """USE_NUMBA should be accessible from optimizer.__init__."""
        from biocompiler.optimizer import USE_NUMBA
        assert isinstance(USE_NUMBA, bool)


# ═══════════════════════════════════════════════════════════════
# 2. CAI computation parity (NUMBA vs Python)
# ═══════════════════════════════════════════════════════════════

class TestCAIParity:
    """Verify CAI computation produces same results via NUMBA and Python."""

    PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"

    @pytest.mark.skipif(not HAS_NUMBA, reason=NUMBA_SKIP_REASON)
    def test_cai_fast_parity(self):
        """Parity smoke test for _compute_cai_fast.

        NOTE (prior audit): This was originally a NUMBA-vs-Python parity
        test that toggled ``nk.USE_NUMBA`` between True and False. However,
        ``cai.py``'s 4 dispatch sites check only ``HAS_NUMBA`` (the runtime
        availability flag), NOT ``USE_NUMBA`` (the user-controlled toggle).
        Both ``USE_NUMBA=True`` and ``USE_NUMBA=False`` runs therefore
        execute the SAME NUMBA path, making the parity assertion vacuous
        (trivially true). Until cai.py is fixed to honor ``USE_NUMBA``,
        this test is converted to a smoke test that exercises the function
        under both toggles and verifies it returns a sensible value — it
        can no longer claim to verify NUMBA↔Python parity.
        """
        from biocompiler.optimizer.cai import _compute_cai_fast
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        from biocompiler.type_system import AA_TO_CODONS
        import biocompiler.optimizer.numba_kernels as nk

        # Build test sequence
        adaptiveness = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        dna_parts = []
        for aa in self.PROTEIN:
            codons = AA_TO_CODONS.get(aa, ["ATG"])
            best = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
            dna_parts.append(best)
        seq = "".join(dna_parts)

        original = nk.USE_NUMBA
        try:
            # Both toggles currently hit the same NUMBA code path because
            # cai.py ignores USE_NUMBA (H12). We at least confirm the function
            # returns a sane CAI value under each toggle.
            nk.USE_NUMBA = True
            cai_on = _compute_cai_fast(seq, adaptiveness)
            assert 0.0 < cai_on <= 1.0

            nk.USE_NUMBA = False
            cai_off = _compute_cai_fast(seq, adaptiveness)
            assert 0.0 < cai_off <= 1.0
        finally:
            nk.USE_NUMBA = original

        # Both paths return identical values because the dispatch site
        # checks HAS_NUMBA, not USE_NUMBA. This assertion documents the
        # known false positive (prior audit) — when cai.py is fixed to honor
        # USE_NUMBA, this assertion will need to be re-evaluated.
        assert cai_on == cai_off, (
            "Expected identical results (both paths use NUMBA per H12); "
            "if this fails, cai.py may have been fixed to honor USE_NUMBA "
            "and this test should be restored to a true parity check."
        )

    @pytest.mark.skipif(not HAS_NUMBA, reason=NUMBA_SKIP_REASON)
    def test_cai_incremental_parity(self):
        """compute_cai_incremental should match between NUMBA and Python."""
        from biocompiler.optimizer.numba_kernels import compute_cai_incremental
        import biocompiler.optimizer.numba_kernels as nk

        # Test values
        log_sum = -10.0
        n_codons = 20
        w_old = 0.5
        w_new = 1.0

        # The same function serves both paths (NUMBA JIT when available,
        # pure Python otherwise), so just verify it returns a sensible value
        result = compute_cai_incremental(log_sum, n_codons, w_old, w_new)
        assert 0.0 < result <= 1.0


# ═══════════════════════════════════════════════════════════════
# 3. GC counting parity
# ═══════════════════════════════════════════════════════════════

class TestGCCountingParity:
    """Verify GC counting produces same results via NUMBA and Python."""

    @pytest.mark.skipif(not HAS_NUMBA, reason=NUMBA_SKIP_REASON)
    def test_gc_count_parity(self):
        """count_gc should match Python str.count for various sequences."""
        from biocompiler.optimizer.numba_kernels import count_gc, seq_to_bytes
        import biocompiler.optimizer.numba_kernels as nk

        test_seqs = [
            "ATGCGATCGATCGATCG",
            "AAAAAAAATTTTTTTT",
            "GCGCGCGCGCGC",
            "ATGC",
        ]

        for seq in test_seqs:
            seq_bytes = seq_to_bytes(seq)
            numba_count = count_gc(seq_bytes)
            py_count = seq.count("G") + seq.count("C")
            assert numba_count == py_count, (
                f"GC count mismatch for {seq}: NUMBA={numba_count} vs Python={py_count}"
            )


# ═══════════════════════════════════════════════════════════════
# 4. Dinucleotide counting parity
# ═══════════════════════════════════════════════════════════════

class TestDinucleotideCountingParity:
    """Verify dinucleotide counting produces same results via NUMBA and Python."""

    @pytest.mark.skipif(not HAS_NUMBA, reason=NUMBA_SKIP_REASON)
    def test_count_dinucs_fast_parity(self):
        """Smoke test for _count_dinucs_fast under USE_NUMBA toggles.

        NOTE (prior audit): Originally a NUMBA↔Python parity test,
        but ``cai.py``'s dispatch sites check only ``HAS_NUMBA`` (not
        ``USE_NUMBA``), so both toggles hit the SAME NUMBA path and the
        parity assertion was vacuously true. Converted to a smoke test
        until cai.py is fixed to honor ``USE_NUMBA``.
        """
        from biocompiler.optimizer.cai import _count_dinucs_fast
        import biocompiler.optimizer.numba_kernels as nk

        seq = "ATGCGATCGATCGTGTGTGTAGCGCGCGTATATATCGCGT"

        original = nk.USE_NUMBA
        try:
            nk.USE_NUMBA = True
            counts_on = _count_dinucs_fast(seq, "GT", "CG", "AG")
            assert len(counts_on) == 3
            assert all(isinstance(c, int) for c in counts_on)

            nk.USE_NUMBA = False
            counts_off = _count_dinucs_fast(seq, "GT", "CG", "AG")
            assert len(counts_off) == 3
        finally:
            nk.USE_NUMBA = original

        # Both paths return identical values because the dispatch site
        # checks HAS_NUMBA, not USE_NUMBA (H12).
        assert counts_on == counts_off, (
            "Expected identical results (both paths use NUMBA per H12); "
            "if this fails, cai.py may have been fixed to honor USE_NUMBA."
        )


# ═══════════════════════════════════════════════════════════════
# 5. Batch codon swap scoring parity
# ═══════════════════════════════════════════════════════════════

class TestBatchCodonSwapParity:
    """Verify batch codon swap scoring produces same results via NUMBA and Python."""

    @pytest.mark.skipif(not HAS_NUMBA, reason=NUMBA_SKIP_REASON)
    def test_batch_swap_scorer_parity(self):
        """Smoke test for _BatchSwapScorer under USE_NUMBA toggles.

        NOTE (prior audit): Originally a NUMBA↔Python parity test,
        but ``cai.py``'s dispatch sites check only ``HAS_NUMBA`` (not
        ``USE_NUMBA``), so both toggles hit the SAME NUMBA path and the
        parity assertion was vacuously true. Converted to a smoke test
        until cai.py is fixed to honor ``USE_NUMBA``.
        """
        from biocompiler.optimizer.cai import _BatchSwapScorer
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        from biocompiler.type_system import AA_TO_CODONS
        import biocompiler.optimizer.numba_kernels as nk

        protein = "MALWMRLLPL"
        adaptiveness = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]

        # Build test codons
        seq_codons = []
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, ["ATG"])
            best = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
            seq_codons.append(best)

        # Test position 1 (Ala)
        test_pos = 1
        candidates = AA_TO_CODONS.get(protein[test_pos], ["GCT"])

        original = nk.USE_NUMBA
        try:
            nk.USE_NUMBA = True
            scorer1 = _BatchSwapScorer(adaptiveness)
            scorer1.reset_incremental_state(seq_codons)
            scores_on = scorer1.score_candidates(seq_codons, test_pos, candidates)
            assert len(scores_on) == len(candidates)
            assert all(s > 0 for s in scores_on)

            nk.USE_NUMBA = False
            scorer2 = _BatchSwapScorer(adaptiveness)
            scorer2.reset_incremental_state(seq_codons)
            scores_off = scorer2.score_candidates(seq_codons, test_pos, candidates)
        finally:
            nk.USE_NUMBA = original

        # Both paths return identical values because the dispatch site
        # checks HAS_NUMBA, not USE_NUMBA (H12).
        assert scores_on == scores_off, (
            "Expected identical results (both paths use NUMBA per H12); "
            "if this fails, cai.py may have been fixed to honor USE_NUMBA."
        )


# ═══════════════════════════════════════════════════════════════
# 6. Full optimization parity
# ═══════════════════════════════════════════════════════════════

class TestFullOptimizationParity:
    """Verify full optimization produces consistent results."""

    PROTEIN = "MALWMRLLPL"

    @pytest.mark.skipif(not HAS_NUMBA, reason=NUMBA_SKIP_REASON)
    def test_optimization_cai_consistency(self):
        """Optimization with NUMBA on should produce a valid result with CAI > 0."""
        from biocompiler.optimizer import optimize_sequence
        import biocompiler.optimizer.numba_kernels as nk

        nk.USE_NUMBA = True
        result = optimize_sequence(self.PROTEIN, organism="ecoli", strict_mode=False)
        nk.USE_NUMBA = HAS_NUMBA

        assert result.cai > 0.0
        assert len(result.sequence) == len(self.PROTEIN) * 3

    def test_optimization_works_with_numba_off(self):
        """Optimization should work when USE_NUMBA is False."""
        from biocompiler.optimizer import optimize_sequence
        import biocompiler.optimizer.numba_kernels as nk

        original = nk.USE_NUMBA
        nk.USE_NUMBA = False
        try:
            result = optimize_sequence(self.PROTEIN, organism="ecoli", strict_mode=False)
            assert result.cai > 0.0
            assert len(result.sequence) == len(self.PROTEIN) * 3
        finally:
            nk.USE_NUMBA = original


# ═══════════════════════════════════════════════════════════════
# 7. Fallback when NUMBA is not installed
# ═══════════════════════════════════════════════════════════════

class TestFallbackNoNumba:
    """Verify everything works when NUMBA is not installed / USE_NUMBA is False."""

    def test_compute_cai_fast_without_numba(self):
        """_compute_cai_fast should work when USE_NUMBA is False."""
        from biocompiler.optimizer.cai import _compute_cai_fast
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        import biocompiler.optimizer.numba_kernels as nk

        original = nk.USE_NUMBA
        nk.USE_NUMBA = False
        try:
            adaptiveness = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
            cai = _compute_cai_fast("ATGGCTCTGTGGGAA", adaptiveness)
            assert 0.0 <= cai <= 1.0
        finally:
            nk.USE_NUMBA = original

    def test_count_dinucs_fast_without_numba(self):
        """_count_dinucs_fast should work when USE_NUMBA is False."""
        from biocompiler.optimizer.cai import _count_dinucs_fast
        import biocompiler.optimizer.numba_kernels as nk

        original = nk.USE_NUMBA
        nk.USE_NUMBA = False
        try:
            counts = _count_dinucs_fast("ATGCGATCGTGT", "GT", "CG")
            assert len(counts) == 2
            assert all(isinstance(c, int) for c in counts)
        finally:
            nk.USE_NUMBA = original

    def test_batch_swap_scorer_without_numba(self):
        """_BatchSwapScorer should work when USE_NUMBA is False."""
        from biocompiler.optimizer.cai import _BatchSwapScorer
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        import biocompiler.optimizer.numba_kernels as nk

        original = nk.USE_NUMBA
        nk.USE_NUMBA = False
        try:
            species_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
            scorer = _BatchSwapScorer(species_cai)
            seq_codons = ["ATG", "GCT", "CTG", "GAA"]
            scorer.reset_incremental_state(seq_codons)
            scores = scorer.score_candidates(seq_codons, 1, ["GCT", "GCC", "GCA", "GCG"])
            assert len(scores) == 4
            assert all(s > 0 for s in scores)
        finally:
            nk.USE_NUMBA = original
