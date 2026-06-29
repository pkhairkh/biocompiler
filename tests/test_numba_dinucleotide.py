"""Tests for fast_dinucleotide_count NUMBA kernel wiring.

Verifies that:
1. The _count_dinucs_fast helper correctly counts dinucleotides (with/without NUMBA)
2. The NUMBA kernel matches pure-Python results for GT, CG, AG, AT dinucleotides
3. NoGTDinucleotide checks use the fast kernel (short-circuit)
4. AvoidCpG checks use the fast kernel (short-circuit)
5. _count_dinucleotides() and _count_gts() in optimization.py use the kernel
6. End-to-end optimizer still produces correct results
"""

import pytest
import importlib


# Use direct module import to avoid __init__ chain issues from other agents
def _get_optimization():
    """Get optimization module directly (avoids __init__ import chain issues)."""
    return importlib.import_module('biocompiler.optimizer')


def _get_type_system():
    """Get type_system module directly."""
    return importlib.import_module('biocompiler.type_system')


# ── Helper: _count_dinucs_fast ──────────────────────────────────────────

class TestCountDinucsFast:
    """Tests for the _count_dinucs_fast wrapper in both optimization and type_system."""

    def test_import_from_optimization(self):
        """_count_dinucs_fast is importable from optimization."""
        mod = _get_optimization()
        assert callable(mod._count_dinucs_fast)

    def test_import_from_type_system(self):
        """_count_dinucs_fast is importable from type_system."""
        mod = _get_type_system()
        assert callable(mod._count_dinucs_fast)

    def test_single_gt(self):
        """Count GT dinucleotides."""
        mod = _get_optimization()
        result = mod._count_dinucs_fast("AACGTCAAGT", "GT")
        assert result == (2,)  # positions 3 and 8

    def test_single_cg(self):
        """Count CG dinucleotides."""
        mod = _get_optimization()
        result = mod._count_dinucs_fast("AACGTCGAA", "CG")
        assert result == (2,)  # positions 2 and 5

    def test_multiple_dinucleotides(self):
        """Count multiple dinucleotides in one pass."""
        mod = _get_optimization()
        seq = "AACGTCGAAGTAT"
        gt, cg, ag = mod._count_dinucs_fast(seq, "GT", "CG", "AG")
        assert gt == 2
        assert cg == 2
        assert ag == 1

    def test_empty_dinucleotides(self):
        """No dinucleotides requested returns empty tuple."""
        mod = _get_optimization()
        result = mod._count_dinucs_fast("AACGT")
        assert result == ()

    def test_no_matches(self):
        """Sequence with no matching dinucleotides."""
        mod = _get_optimization()
        result = mod._count_dinucs_fast("AACCCC", "GT")
        assert result == (0,)

    def test_short_sequence(self):
        """Single character sequence has no dinucleotides."""
        mod = _get_optimization()
        result = mod._count_dinucs_fast("A", "GT")
        assert result == (0,)

    def test_empty_sequence(self):
        """Empty sequence has no dinucleotides."""
        mod = _get_optimization()
        result = mod._count_dinucs_fast("", "GT")
        assert result == (0,)

    def test_overlapping_dinucleotides(self):
        """Overlapping dinucleotides are counted correctly."""
        mod = _get_optimization()
        # CGCG has CG at positions 0,2
        result = mod._count_dinucs_fast("CGCG", "CG")
        assert result == (2,)

    def test_type_system_matches_optimization(self):
        """Both modules produce the same results."""
        opt_mod = _get_optimization()
        ts_mod = _get_type_system()
        seq = "AACGTCGAAGTATCGGT"
        for dinucs in [("GT",), ("CG",), ("GT", "CG", "AG")]:
            assert opt_mod._count_dinucs_fast(seq, *dinucs) == ts_mod._count_dinucs_fast(seq, *dinucs)

    def test_matches_pure_python(self):
        """NUMBA kernel matches pure-Python counting for various sequences."""
        mod = _get_optimization()
        test_seqs = [
            "ATGCGTACGTTGCGATCGATCG",
            "GGGGTTTTCCCCAAAA",
            "GTCGCGCGTGTGT",
            "A",
            "",
            "GTGTGTGTGT",
        ]
        for seq in test_seqs:
            for di in ["GT", "CG", "AG", "AT", "TA", "GC"]:
                # Pure-Python reference
                expected = 0
                pos = 0
                while True:
                    pos = seq.find(di, pos)
                    if pos == -1:
                        break
                    expected += 1
                    pos += 1
                result = mod._count_dinucs_fast(seq, di)[0]
                assert result == expected, (
                    f"Mismatch for {di!r} in {seq!r}: "
                    f"expected {expected}, got {result}"
                )


# ── NoGTDinucleotide predicate with fast kernel ────────────────────────

class TestNoGTDinucleotideFastKernel:
    """Tests that NoGTDinucleotide checks use the fast kernel."""

    def test_check_no_gt_passes_no_gts(self):
        """Sequences without GT pass the check quickly."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_gt_dinucleotide("AACCCCAAAA")
        assert result.passed is True

    def test_check_no_gt_fails_with_gts(self):
        """Sequences with GT fail the strict check."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_gt_dinucleotide("AAGTCC")
        assert result.passed is False
        assert result.positions == [2]

    def test_check_no_gt_soft_eukaryote_no_gts(self):
        """Soft check passes for eukaryote with no GTs."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_gt_dinucleotide_soft("AACCCCAAAA", organism="Homo_sapiens")
        assert result.passed is True

    def test_check_no_gt_soft_eukaryote_with_gts(self):
        """Soft check for eukaryote with GTs returns appropriate verdict."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_gt_dinucleotide_soft("AAGTCC", organism="Homo_sapiens")
        # Should not crash; verdict depends on max_gt_count
        assert result.predicate == "NoGTDinucleotide"

    def test_check_no_avoidable_gt_no_gts(self):
        """Avoidable GT check passes with no GTs."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_avoidable_gt("AACCCCAAAA", organism="Homo_sapiens")
        assert result.passed is True

    def test_check_no_avoidable_gt_prokaryote(self):
        """Avoidable GT check is skipped for prokaryotes."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_avoidable_gt("AAGTCC", organism="Escherichia_coli")
        assert result.passed is True


# ── AvoidCpG predicate with fast kernel ─────────────────────────────────

class TestAvoidCpGFastKernel:
    """Tests that AvoidCpG checks use the fast kernel."""

    def test_check_no_cpg_no_cg(self):
        """Sequence with no CG dinucleotides passes immediately."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_cpg_island("AATTTTAAAA", organism="Homo_sapiens")
        assert result.passed is True

    def test_check_no_cpg_prokaryote_skip(self):
        """CpG check is skipped for prokaryotes."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_cpg_island("CGCGCG", organism="Escherichia_coli")
        assert result.passed is True

    def test_check_no_cpg_with_cg_short_seq(self):
        """Short sequence with CG does not crash."""
        ts_mod = _get_type_system()
        result = ts_mod.check_no_cpg_island("AACGTT", organism="Homo_sapiens")
        assert result.predicate == "NoCpGIsland"


# ── Optimization-level function tests ───────────────────────────────────

class TestOptimizationDinucleotideFunctions:
    """Tests that optimization.py functions use the fast kernel."""

    def test_count_dinucleotides(self):
        """_count_dinucleotides uses the kernel."""
        mod = _get_optimization()
        assert mod._count_dinucleotides("AACGTCCGT", "GT") == 2
        assert mod._count_dinucleotides("AACGTCCGT", "CG") == 2  # CG at positions 2, 6
        assert mod._count_dinucleotides("AAAAAAAA", "GT") == 0

    def test_count_gts(self):
        """_count_gts uses the kernel."""
        mod = _get_optimization()
        assert mod._count_gts("AACGTCCGT") == 2
        assert mod._count_gts("AAAAAAAA") == 0
        assert mod._count_gts("GTGTGT") == 3

    def test_count_dinucleotides_matches_str_find(self):
        """_count_dinucleotides matches pure Python str.find-based counting."""
        mod = _get_optimization()
        seq = "ATGCGTACGTTGCGATCGATCG"
        for di in ["GT", "CG", "AG", "AT", "TA", "GC"]:
            expected = 0
            pos = 0
            while True:
                pos = seq.find(di, pos)
                if pos == -1:
                    break
                expected += 1
                pos += 1
            assert mod._count_dinucleotides(seq, di) == expected

    def test_count_gts_matches_manual(self):
        """_count_gts matches manual counting."""
        mod = _get_optimization()
        seq = "ATGGGTACGGTTGATCG"
        expected = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "GT")
        assert mod._count_gts(seq) == expected


# ── End-to-end optimizer test ───────────────────────────────────────────

class TestOptimizerWithFastKernel:
    """End-to-end tests ensuring the optimizer works with the wired kernel."""

    def test_optimize_ecoli_insulin(self):
        """E. coli insulin optimization produces valid result."""
        mod = _get_optimization()
        result = mod.optimize_sequence(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT",
            organism="e_coli",
            strict_mode=False,
        )
        assert result.cai > 0.9
        assert len(result.sequence) == 55 * 3

    def test_optimize_human_insulin(self):
        """Human insulin optimization produces valid result."""
        mod = _get_optimization()
        result = mod.optimize_sequence(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT",
            organism="human",
            strict_mode=False,
        )
        assert result.cai > 0.9
        assert len(result.sequence) == 55 * 3


# ── NUMBA availability test ─────────────────────────────────────────────

class TestNumbaAvailability:
    """Tests for NUMBA availability detection."""

    def test_has_numba_flag_optimization(self):
        """optimization.HAS_NUMBA is a bool."""
        mod = _get_optimization()
        assert isinstance(mod.HAS_NUMBA, bool)

    def test_has_numba_flag_type_system(self):
        """type_system._HAS_NUMBA is a bool."""
        mod = _get_type_system()
        assert isinstance(mod._HAS_NUMBA, bool)

    def test_kernel_import_with_fallback(self):
        """fast_dinucleotide_count import falls back gracefully."""
        from biocompiler.optimizer.numba_kernels import fast_dinucleotide_count
        assert callable(fast_dinucleotide_count)
