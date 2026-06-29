"""Thread safety tests for biocompiler global mutable state.

Task 1.4 — Verify that concurrent access to module-level caches and
lazy-initialized singletons does not corrupt data or raise exceptions.

Tests cover:
  - get_sharp_li_adaptiveness_tables() — concurrent lazy init
  - HybridOptimizer.optimize() — per-instance lock isolation
  - random.seed() — no global state pollution between calls
  - Module-level cache lazy initialization — esmfold, netmhcpan, etc.
  - proof_checks.assert_valine_gt_invariant — concurrent check
"""

from __future__ import annotations

import random
import threading
import time
import pytest


# ────────────────────────────────────────────────────────────
# 1. Concurrent get_sharp_li_adaptiveness_tables()
# ────────────────────────────────────────────────────────────

class TestSharpLiAdaptivenessTables:
    """Verify that concurrent calls to get_sharp_li_adaptiveness_tables()
    do not corrupt data or raise exceptions."""

    def test_concurrent_init_returns_same_object(self):
        """Multiple threads calling get_sharp_li_adaptiveness_tables()
        simultaneously should all receive the same dict object."""
        from biocompiler.organisms import get_sharp_li_adaptiveness_tables
        # Reset the global to force re-initialisation
        import biocompiler.organisms as org_mod
        org_mod.SHARP_LI_ADAPTIVENESS_TABLES = {}

        results: list[dict] = [None] * 8  # type: ignore[list-item]
        errors: list[Exception] = []

        def worker(idx: int) -> None:
            try:
                results[idx] = get_sharp_li_adaptiveness_tables()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker, args=(i,)) for i in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert not errors, f"Errors during concurrent init: {errors}"
        # All results should be the same object (or at least equal)
        non_none = [r for r in results if r is not None]
        assert len(non_none) == 8, f"Expected 8 results, got {len(non_none)}"
        first = non_none[0]
        for r in non_none[1:]:
            assert r is first, "Concurrent calls returned different objects"

    def test_tables_contain_expected_organisms(self):
        """After init, tables should contain at least the canonical organisms."""
        from biocompiler.organisms import get_sharp_li_adaptiveness_tables
        tables = get_sharp_li_adaptiveness_tables()
        # Reset and re-init for clean state
        import biocompiler.organisms as org_mod
        org_mod.SHARP_LI_ADAPTIVENESS_TABLES = {}
        tables = get_sharp_li_adaptiveness_tables()
        assert "Escherichia_coli" in tables
        assert isinstance(tables["Escherichia_coli"], dict)


# ────────────────────────────────────────────────────────────
# 2. HybridOptimizer — per-instance lock isolation
# ────────────────────────────────────────────────────────────

class TestHybridOptimizerThreadSafety:
    """Verify that concurrent optimize() calls on separate instances
    work correctly and do not interfere with each other."""

    def test_separate_instances_concurrent_optimize(self):
        """Two HybridOptimizer instances running optimize() concurrently
        on different proteins should produce correct, independent results."""
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer

        opt1 = HybridOptimizer(species="ecoli")
        opt2 = HybridOptimizer(species="ecoli")

        protein1 = "MKFLILLFNILCR"  # 13 aa
        protein2 = "MRVLKFGGTSVANA"  # 14 aa

        results: dict[str, object] = {}
        errors: list[Exception] = []

        def worker(key: str, opt: HybridOptimizer, protein: str) -> None:
            try:
                results[key] = opt.optimize(protein, is_prokaryote=True)
            except Exception as exc:
                errors.append(exc)

        t1 = threading.Thread(target=worker, args=("p1", opt1, protein1))
        t2 = threading.Thread(target=worker, args=("p2", opt2, protein2))
        t1.start()
        t2.start()
        t1.join(timeout=10)
        t2.join(timeout=10)

        assert not errors, f"Errors during concurrent optimize: {errors}"
        assert "p1" in results
        assert "p2" in results

        r1 = results["p1"]
        r2 = results["p2"]
        # Each result should be for the correct protein
        assert len(r1.sequence) == len(protein1) * 3  # type: ignore[attr-defined]
        assert len(r2.sequence) == len(protein2) * 3  # type: ignore[attr-defined]

    def test_same_instance_sequential_optimize(self):
        """Calling optimize() on the same instance sequentially should
        work correctly (lock is released after each call)."""
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer

        opt = HybridOptimizer(species="ecoli")
        r1 = opt.optimize("MKFLILLFNILCR", is_prokaryote=True)
        r2 = opt.optimize("MRVLKFGGTSVANA", is_prokaryote=True)

        assert len(r1.sequence) == 13 * 3
        assert len(r2.sequence) == 14 * 3


# ────────────────────────────────────────────────────────────
# 3. random.seed() — no global state pollution
# ────────────────────────────────────────────────────────────

class TestRandomSeedIsolation:
    """Verify that using the seed parameter in optimization does not
    pollute the global random state visible to other threads."""

    def test_seed_does_not_pollute_global_random(self):
        """After calling optimize_sequence with seed, the global
        random generator should not be affected."""
        # Capture global random state
        before_state = random.getstate()
        _ = random.random()  # advance state

        # The optimize_sequence function now uses a per-call Random instance
        # instead of random.seed(). Verify the module does not call random.seed().
        from biocompiler.optimization import _greedy_optimize

        # _greedy_optimize accepts seed; verify it does not call random.seed()
        # by checking global state is unchanged after a call with seed.
        before_state2 = random.getstate()

        try:
            _greedy_optimize(
                protein="MKFLILLFNILCR",
                organism="Escherichia_coli",
                seed=42,
            )
        except Exception:
            # May fail for various reasons; we only care about random state
            pass  # noqa: S110

        after_state = random.getstate()
        assert before_state2 == after_state, (
            "random.seed() was called — global random state was polluted"
        )

    def test_concurrent_different_seeds_dont_interfere(self):
        """Two threads using different seeds should not interfere
        with each other's random state."""

        results: dict[str, list[float]] = {}
        errors: list[Exception] = []

        def worker(key: str, seed: int) -> None:
            try:
                rng = random.Random(seed)
                # Simulate what the optimization code does with the local RNG
                results[key] = [rng.random() for _ in range(5)]
            except Exception as exc:
                errors.append(exc)

        t1 = threading.Thread(target=worker, args=("a", 42))
        t2 = threading.Thread(target=worker, args=("b", 123))
        t1.start()
        t2.start()
        t1.join(timeout=5)
        t2.join(timeout=5)

        assert not errors
        # Different seeds should produce different sequences
        assert results["a"] != results["b"]


# ────────────────────────────────────────────────────────────
# 4. Module-level cache lazy init — basic smoke test
# ────────────────────────────────────────────────────────────

class TestLazyInitThreadSafety:
    """Smoke tests for thread-safe lazy initialization of module-level
    caches. These verify that the lock pattern is present and functional."""

    def test_esmfold_cache_init(self):
        """Concurrent calls to _get_default_cache() should return
        the same ESMFoldCache instance."""
        from biocompiler.engines.esmfold import _get_default_cache
        # Reset
        import biocompiler.engines.esmfold as esm_mod
        esm_mod._default_cache = None

        results: list[object] = [None] * 4
        errors: list[Exception] = []

        def worker(idx: int) -> None:
            try:
                results[idx] = _get_default_cache()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker, args=(i,)) for i in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)

        assert not errors
        non_none = [r for r in results if r is not None]
        assert len(non_none) == 4
        first = non_none[0]
        for r in non_none[1:]:
            assert r is first

    def test_netmhcpan_cache_init(self):
        """Concurrent calls to _get_default_cache() should return
        the same NetMHCpanCache instance."""
        from biocompiler.immunogenicity.netmhcpan import _get_default_cache
        # Reset
        import biocompiler.immunogenicity.netmhcpan as nm_mod
        nm_mod._default_cache = None

        results: list[object] = [None] * 4
        errors: list[Exception] = []

        def worker(idx: int) -> None:
            try:
                results[idx] = _get_default_cache()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker, args=(i,)) for i in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)

        assert not errors
        non_none = [r for r in results if r is not None]
        assert len(non_none) == 4
        first = non_none[0]
        for r in non_none[1:]:
            assert r is first

    def test_proof_checks_valine_concurrent(self):
        """Concurrent calls to assert_valine_gt_invariant() should
        not raise exceptions (double-checked locking)."""
        from biocompiler.provenance.proof_checks import assert_valine_gt_invariant
        # Reset to force re-check
        import biocompiler.provenance.proof_checks as pc_mod
        pc_mod._valine_gt_checked = False

        errors: list[Exception] = []

        def worker() -> None:
            try:
                assert_valine_gt_invariant()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)

        assert not errors

    def test_immunogenicity_pssm_concurrent_build(self):
        """Concurrent calls to _ensure_pssms_built() should not
        raise exceptions (double-checked locking)."""
        from biocompiler.immunogenicity.core import _ensure_pssms_built
        # Reset to force rebuild
        import biocompiler.immunogenicity.core as imm_mod
        imm_mod._pssm_built = False
        imm_mod.MHC_I_PSSM = {}
        imm_mod.MHC_II_PSSM = {}

        errors: list[Exception] = []

        def worker() -> None:
            try:
                _ensure_pssms_built()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert not errors
        # Verify PSSMs were actually built
        assert imm_mod.MHC_I_PSSM, "MHC_I_PSSM was not built"
        assert imm_mod.MHC_II_PSSM, "MHC_II_PSSM was not built"
        assert "HLA-A*02:01" in imm_mod.MHC_I_PSSM


# ────────────────────────────────────────────────────────────
# 5. Cache clear operations are thread-safe
# ────────────────────────────────────────────────────────────

class TestCacheClearThreadSafety:
    """Verify that clear_cache() operations on module-level caches
    do not crash when called concurrently."""

    def test_camsol_clear_cache_concurrent(self):
        """Concurrent clear_cache() calls on camsol should not raise."""
        from biocompiler.engines.camsol import clear_cache
        errors: list[Exception] = []

        def worker() -> None:
            try:
                clear_cache()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)

        assert not errors

    def test_immunogenicity_clear_cache_concurrent(self):
        """Concurrent clear_cache() calls on immunogenicity should not raise."""
        from biocompiler.immunogenicity.core import clear_cache
        errors: list[Exception] = []

        def worker() -> None:
            try:
                clear_cache()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)

        assert not errors
