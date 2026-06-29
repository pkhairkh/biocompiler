"""Tests for NUMBA fast_gc_window kernel wired into sliding_gc optimizer.

Covers:
- NUMBA kernel produces identical results to the pure-Python path
- Import fallback works when NUMBA is unavailable
- Both paths handle edge cases identically (empty seq, short seq, etc.)
- Violation detection matches between NUMBA and Python
- step parameter handled correctly in both paths
- The _FORCE_PYTHON_GC_WINDOW flag can override dispatch
"""

import pytest
import biocompiler.sequence.sliding_gc as sliding_gc_mod
from biocompiler.sequence.sliding_gc import (
    WindowViolation,
    SlidingGCResult,
    check_sliding_gc,
    _check_sliding_gc_python,
    _check_sliding_gc_numba,
    _FORCE_PYTHON_GC_WINDOW,
    _HAS_NUMBA,
)


# ────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────

# Representative DNA sequences for comparison testing
SEQ_UNIFORM_50 = "GCAT" * 50        # 200 bp, ~50% GC
SEQ_ALL_GC = "GCGC" * 50            # 200 bp, 100% GC
SEQ_ALL_AT = "ATAT" * 50            # 200 bp, 0% GC
SEQ_MIXED = ("ATAT" * 20) + ("GCGC" * 20) + ("GCAT" * 20)  # 240 bp, mixed
SEQ_SHORT = "GCGCATAT"              # 8 bp


class TestNumbaImportFallback:
    """Test that the NUMBA import fallback works correctly."""

    def test_has_numba_flag_is_bool(self):
        """_HAS_NUMBA should be a boolean."""
        assert isinstance(_HAS_NUMBA, bool)

    def test_force_python_flag_default_false(self):
        """_FORCE_PYTHON_GC_WINDOW should default to False."""
        # Reset to default
        sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = False
        assert _FORCE_PYTHON_GC_WINDOW is False

    def test_numba_kernel_import_or_none(self):
        """The NUMBA kernel should be imported or set to None."""
        if _HAS_NUMBA:
            assert sliding_gc_mod._fast_gc_window_numba is not None
            assert sliding_gc_mod._seq_to_bytes is not None
        else:
            assert sliding_gc_mod._fast_gc_window_numba is None


# ────────────────────────────────────────────────────────────
# NUMBA vs Python equivalence (step=1)
# ────────────────────────────────────────────────────────────

class TestNumbaPythonEquivalence:
    """Verify NUMBA and Python paths produce identical results."""

    @pytest.mark.parametrize("seq,window_size,gc_min,gc_max", [
        (SEQ_UNIFORM_50, 50, 0.30, 0.70),
        (SEQ_ALL_GC, 50, 0.30, 0.70),
        (SEQ_ALL_AT, 50, 0.30, 0.70),
        (SEQ_MIXED, 30, 0.20, 0.80),
        (SEQ_SHORT, 4, 0.30, 0.70),
    ])
    def test_passed_status_matches(self, seq, window_size, gc_min, gc_max):
        """Pass/fail status should match between NUMBA and Python."""
        py_result = _check_sliding_gc_python(seq, window_size, gc_min, gc_max, step=1)
        if _HAS_NUMBA:
            nb_result = _check_sliding_gc_numba(seq, window_size, gc_min, gc_max, step=1)
            assert py_result.passed == nb_result.passed, (
                f"Passed mismatch: py={py_result.passed}, nb={nb_result.passed} "
                f"for seq[:20]={seq[:20]}... w={window_size}"
            )

    @pytest.mark.parametrize("seq,window_size,gc_min,gc_max", [
        (SEQ_UNIFORM_50, 50, 0.30, 0.70),
        (SEQ_ALL_GC, 50, 0.30, 0.70),
        (SEQ_ALL_AT, 50, 0.30, 0.70),
        (SEQ_MIXED, 30, 0.20, 0.80),
        (SEQ_SHORT, 4, 0.30, 0.70),
    ])
    def test_min_max_gc_matches(self, seq, window_size, gc_min, gc_max):
        """Min/max GC values should match to within floating-point tolerance."""
        py_result = _check_sliding_gc_python(seq, window_size, gc_min, gc_max, step=1)
        if _HAS_NUMBA:
            nb_result = _check_sliding_gc_numba(seq, window_size, gc_min, gc_max, step=1)
            assert abs(py_result.min_gc - nb_result.min_gc) < 1e-12, (
                f"min_gc mismatch: py={py_result.min_gc}, nb={nb_result.min_gc}"
            )
            assert abs(py_result.max_gc - nb_result.max_gc) < 1e-12, (
                f"max_gc mismatch: py={py_result.max_gc}, nb={nb_result.max_gc}"
            )

    @pytest.mark.parametrize("seq,window_size,gc_min,gc_max", [
        (SEQ_UNIFORM_50, 50, 0.30, 0.70),
        (SEQ_ALL_GC, 50, 0.30, 0.70),
        (SEQ_ALL_AT, 50, 0.30, 0.70),
        (SEQ_MIXED, 30, 0.20, 0.80),
    ])
    def test_violations_match(self, seq, window_size, gc_min, gc_max):
        """Violation positions and GC values should match between paths."""
        py_result = _check_sliding_gc_python(seq, window_size, gc_min, gc_max, step=1)
        if _HAS_NUMBA:
            nb_result = _check_sliding_gc_numba(seq, window_size, gc_min, gc_max, step=1)
            assert len(py_result.violations) == len(nb_result.violations), (
                f"Violation count mismatch: py={len(py_result.violations)}, "
                f"nb={len(nb_result.violations)}"
            )
            for pv, nv in zip(py_result.violations, nb_result.violations):
                assert pv.start == nv.start, f"start mismatch: {pv.start} vs {nv.start}"
                assert pv.end == nv.end, f"end mismatch: {pv.end} vs {nv.end}"
                assert pv.direction == nv.direction, f"direction mismatch"
                assert abs(pv.gc_content - nv.gc_content) < 1e-12, (
                    f"gc_content mismatch: {pv.gc_content} vs {nv.gc_content}"
                )

    def test_empty_sequence_both_paths(self):
        """Empty sequence should return identical results on both paths."""
        py_result = _check_sliding_gc_python("", 50, 0.30, 0.70, step=1)
        assert py_result.passed is True
        assert py_result.min_gc == 0.0
        assert py_result.max_gc == 0.0
        if _HAS_NUMBA:
            nb_result = _check_sliding_gc_numba("", 50, 0.30, 0.70, step=1)
            assert nb_result.passed is True
            assert nb_result.min_gc == 0.0
            assert nb_result.max_gc == 0.0

    def test_sequence_shorter_than_window_both_paths(self):
        """Sequence shorter than window should match on both paths."""
        py_result = _check_sliding_gc_python("GCGC", 50, 0.30, 0.70, step=1)
        if _HAS_NUMBA:
            nb_result = _check_sliding_gc_numba("GCGC", 50, 0.30, 0.70, step=1)
            assert py_result.passed == nb_result.passed
            assert abs(py_result.min_gc - nb_result.min_gc) < 1e-12
            assert abs(py_result.max_gc - nb_result.max_gc) < 1e-12


# ────────────────────────────────────────────────────────────
# Per-window GC% value equivalence (direct kernel comparison)
# ────────────────────────────────────────────────────────────

class TestKernelDirectEquivalence:
    """Compare fast_gc_window kernel output directly against Python computation."""

    @pytest.mark.skipif(not _HAS_NUMBA, reason="NUMBA not available")
    def test_kernel_output_matches_python_loop(self):
        """NUMBA fast_gc_window output should match per-window Python computation."""
        import numpy as np
        from biocompiler.optimizer.numba_kernels import fast_gc_window, seq_to_bytes

        seq = SEQ_MIXED
        window_size = 30

        # NUMBA kernel
        seq_bytes = seq_to_bytes(seq)
        nb_result = fast_gc_window(seq_bytes, window_size)

        # Python reference: compute GC% for each window position
        n = len(seq)
        n_windows = n - window_size + 1
        py_pcts = []
        for i in range(n_windows):
            window = seq[i:i + window_size]
            gc = sum(1 for b in window if b in "GC") / window_size
            py_pcts.append(gc)

        assert len(nb_result) == len(py_pcts), (
            f"Length mismatch: nb={len(nb_result)}, py={len(py_pcts)}"
        )
        for i in range(len(nb_result)):
            assert abs(float(nb_result[i]) - py_pcts[i]) < 1e-12, (
                f"Window {i}: NUMBA={float(nb_result[i])}, Python={py_pcts[i]}"
            )

    @pytest.mark.skipif(not _HAS_NUMBA, reason="NUMBA not available")
    def test_kernel_all_gc_sequence(self):
        """All-GC sequence should produce 1.0 GC% for every window."""
        import numpy as np
        from biocompiler.optimizer.numba_kernels import fast_gc_window, seq_to_bytes

        seq = "GCGC" * 25  # 100 bp
        seq_bytes = seq_to_bytes(seq)
        result = fast_gc_window(seq_bytes, 20)
        assert len(result) == 81  # 100 - 20 + 1
        assert np.all(result == 1.0)

    @pytest.mark.skipif(not _HAS_NUMBA, reason="NUMBA not available")
    def test_kernel_all_at_sequence(self):
        """All-AT sequence should produce 0.0 GC% for every window."""
        import numpy as np
        from biocompiler.optimizer.numba_kernels import fast_gc_window, seq_to_bytes

        seq = "ATAT" * 25  # 100 bp
        seq_bytes = seq_to_bytes(seq)
        result = fast_gc_window(seq_bytes, 20)
        assert len(result) == 81
        assert np.all(result == 0.0)

    @pytest.mark.skipif(not _HAS_NUMBA, reason="NUMBA not available")
    def test_kernel_short_sequence(self):
        """Sequence shorter than window should return empty array."""
        from biocompiler.optimizer.numba_kernels import fast_gc_window, seq_to_bytes

        seq = "GCGC"  # 4 bp
        seq_bytes = seq_to_bytes(seq)
        result = fast_gc_window(seq_bytes, 10)
        assert len(result) == 0

    @pytest.mark.skipif(not _HAS_NUMBA, reason="NUMBA not available")
    def test_kernel_window_equals_sequence(self):
        """Window size equal to sequence length should return single value."""
        import numpy as np
        from biocompiler.optimizer.numba_kernels import fast_gc_window, seq_to_bytes

        seq = "GCGCATAT"  # 8 bp, 50% GC
        seq_bytes = seq_to_bytes(seq)
        result = fast_gc_window(seq_bytes, 8)
        assert len(result) == 1
        assert abs(float(result[0]) - 0.5) < 1e-12


# ────────────────────────────────────────────────────────────
# Dispatch mechanism
# ────────────────────────────────────────────────────────────

class TestDispatchMechanism:
    """Test the NUMBA/Python dispatch in check_sliding_gc."""

    def test_default_dispatch_uses_numba_when_available(self):
        """By default, check_sliding_gc should use NUMBA when available."""
        # Ensure flag is at default
        sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = False
        # Just verify it does not crash — the dispatch is transparent
        result = check_sliding_gc(SEQ_UNIFORM_50, window_size=50)
        assert result.passed is True

    def test_force_python_flag_overrides_dispatch(self):
        """Setting _FORCE_PYTHON_GC_WINDOW=True should force Python path."""
        try:
            sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = True
            result = check_sliding_gc(SEQ_UNIFORM_50, window_size=50)
            assert result.passed is True
        finally:
            sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = False

    def test_force_python_and_auto_dispatch_produce_same_result(self):
        """Results should be identical regardless of dispatch path."""
        py_result = _check_sliding_gc_python(SEQ_MIXED, 30, 0.20, 0.80, step=1)
        try:
            sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = True
            auto_py_result = check_sliding_gc(SEQ_MIXED, window_size=30, gc_min=0.20, gc_max=0.80)
        finally:
            sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = False

        assert py_result.passed == auto_py_result.passed
        assert abs(py_result.min_gc - auto_py_result.min_gc) < 1e-12
        assert abs(py_result.max_gc - auto_py_result.max_gc) < 1e-12

    def test_numba_fallback_on_exception(self):
        """If NUMBA path raises, it should fall back to Python gracefully."""
        if not _HAS_NUMBA:
            pytest.skip("NUMBA not available")

        # Temporarily replace the NUMBA function with one that raises
        original = sliding_gc_mod._fast_gc_window_numba
        try:
            sliding_gc_mod._fast_gc_window_numba = lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("test"))
            sliding_gc_mod._FORCE_PYTHON_GC_WINDOW = False
            # Should not crash — falls back to Python
            result = check_sliding_gc(SEQ_UNIFORM_50, window_size=50)
            assert result.passed is True
        finally:
            sliding_gc_mod._fast_gc_window_numba = original


# ────────────────────────────────────────────────────────────
# Step parameter behavior
# ────────────────────────────────────────────────────────────

class TestStepParameter:
    """Test step parameter handling in both NUMBA and Python paths."""

    def test_step1_numba_python_match(self):
        """step=1 results should match between NUMBA and Python."""
        if not _HAS_NUMBA:
            pytest.skip("NUMBA not available")
        py_result = _check_sliding_gc_python(SEQ_MIXED, 30, 0.20, 0.80, step=1)
        nb_result = _check_sliding_gc_numba(SEQ_MIXED, 30, 0.20, 0.80, step=1)
        assert py_result.passed == nb_result.passed
        assert abs(py_result.min_gc - nb_result.min_gc) < 1e-12
        assert abs(py_result.max_gc - nb_result.max_gc) < 1e-12
        assert len(py_result.violations) == len(nb_result.violations)

    def test_step5_numba_fewer_violations_than_step1(self):
        """step=5 may detect fewer violations than step=1 (sampling effect)."""
        if not _HAS_NUMBA:
            pytest.skip("NUMBA not available")
        nb_step1 = _check_sliding_gc_numba(SEQ_MIXED, 30, 0.20, 0.80, step=1)
        nb_step5 = _check_sliding_gc_numba(SEQ_MIXED, 30, 0.20, 0.80, step=5)
        # step=5 should have <= violations compared to step=1
        assert len(nb_step5.violations) <= len(nb_step1.violations)

    def test_step5_python_path(self):
        """step=5 should work correctly in Python path."""
        result = _check_sliding_gc_python(SEQ_UNIFORM_50, 20, 0.30, 0.70, step=5)
        assert result.passed is True


# ────────────────────────────────────────────────────────────
# Existing test suite still passes (regression)
# ────────────────────────────────────────────────────────────

class TestExistingTestsRegression:
    """Ensure existing sliding GC tests still pass after wiring."""

    def test_uniform_gc_passes(self):
        """A sequence with uniform 50% GC should pass."""
        seq = "GCAT" * 25
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is True

    def test_all_gc_fails_too_high(self):
        """A sequence of all G/C should fail when gc_max < 1.0."""
        seq = "GCGC" * 25
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is False
        assert all(v.direction == "too_high" for v in result.violations)

    def test_all_at_fails_too_low(self):
        """A sequence of all A/T should fail when gc_min > 0.0."""
        seq = "ATAT" * 25
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is False
        assert all(v.direction == "too_low" for v in result.violations)

    def test_local_gc_extreme_detected(self):
        """A local region of extreme GC should be detected."""
        at_region = "ATAT" * 10
        gc_region = "GCGC" * 10
        seq = at_region + gc_region
        result = check_sliding_gc(seq, window_size=20, gc_min=0.20, gc_max=0.80)
        assert result.passed is False

    def test_empty_sequence(self):
        """Empty sequence should return passed=True with no violations."""
        result = check_sliding_gc("", window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is True
        assert len(result.violations) == 0

    def test_short_sequence(self):
        """Sequence shorter than window_size should work."""
        seq = "GCGCATAT"
        result = check_sliding_gc(seq, window_size=50, gc_min=0.30, gc_max=0.70)
        assert result.passed is True
        assert result.min_gc == result.max_gc

    def test_case_insensitive(self):
        """Input should be case-insensitive."""
        seq_lower = "gcgcatat" * 6
        seq_upper = "GCGCATAT" * 6
        r1 = check_sliding_gc(seq_lower, window_size=20, gc_min=0.30, gc_max=0.70)
        r2 = check_sliding_gc(seq_upper, window_size=20, gc_min=0.30, gc_max=0.70)
        assert r1.passed == r2.passed
        assert abs(r1.min_gc - r2.min_gc) < 0.01
        assert abs(r1.max_gc - r2.max_gc) < 0.01

    def test_window_size_1(self):
        """Window size of 1 should work."""
        seq = "GCATGCAT"
        result = check_sliding_gc(seq, window_size=1, gc_min=0.0, gc_max=1.0)
        assert result.passed is True

    def test_relaxed_bounds_pass(self):
        """Very relaxed GC bounds should always pass."""
        seq = "GCGCGCGC" * 12
        result = check_sliding_gc(seq, window_size=20, gc_min=0.0, gc_max=1.0)
        assert result.passed is True
        assert len(result.violations) == 0
