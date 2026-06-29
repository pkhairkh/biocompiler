"""Test that the optimizer never hangs — it times out gracefully.

Task TIGHTEN-2: the optimizer's greedy constraint-fixing loops can
iterate forever on pathological inputs (e.g. EC_lacI, 399 aa, killed
the process with OOM / infinite loop in testing).  ``optimize_sequence``
now wraps its entire body in a ``signal.alarm``-based timeout so it
either completes or returns a partial result with ``timed_out=True``.

These tests verify:
    1. Normal short inputs still complete (timeout doesn't break things).
    2. Even with a short timeout on a large input, the function returns
       within a reasonable wall-clock bound (it never hangs).
    3. A 1-second timeout on a 2000-aa input does not crash the process.
    4. The SIGALRM mechanism actually fires when the optimizer is
       blocked in pure-Python code (positive test of the timeout).
    5. The SIGALRM handler is restored to its prior state after each
       call (no global side effects).
"""

import signal
import time

import pytest

from biocompiler.optimizer.pipeline_core import (
    optimize_sequence,
    OptimizerTimeout,
    _timeout_handler,
)
from biocompiler.optimizer.utils import OPTIMIZER_TIMEOUT_SECONDS


class TestOptimizerTimeout:
    """Verify the optimizer's wall-clock timeout behaves correctly."""

    def test_normal_completes(self):
        """A short, well-behaved input should complete normally.

        The timeout must not interfere with ordinary optimization.
        """
        result = optimize_sequence(
            "MVHLTPEEK", organism="human", strict_mode=False, timeout_seconds=10
        )
        assert result.timed_out is False
        assert result.optimized_sequence is not None
        # Length invariant: DNA length == protein length * 3
        assert len(result.optimized_sequence) == 9 * 3

    def test_short_timeout_returns(self):
        """A short timeout on a large input must not hang the process.

        Uses a 1000-aa input with a 5-second timeout.  The optimizer
        may complete (in which case ``timed_out`` is False) or it may
        time out (in which case ``timed_out`` is True) — either way,
        the function MUST return within 15 seconds (3x the timeout,
        generous buffer for partial-result construction).
        """
        large = "M" * 500 + "A" * 500
        start = time.time()
        result = optimize_sequence(
            large, organism="human", strict_mode=False, timeout_seconds=5
        )
        elapsed = time.time() - start
        assert elapsed < 15, (
            f"Optimizer did not return within 15s wall-clock "
            f"(elapsed={elapsed:.2f}s) — possible hang"
        )
        assert result is not None

    def test_timeout_does_not_crash(self):
        """A very short timeout on a very large input must not crash.

        2000-aa input with a 1-second timeout.  The function must
        return *something* (either a completed result or a partial
        timed-out result) — never raise an unhandled exception or
        kill the process.
        """
        result = optimize_sequence(
            "M" * 2000, organism="human", strict_mode=False, timeout_seconds=1
        )
        assert result is not None
        # Either completed normally or timed out — both are acceptable.
        assert result.timed_out is True or result.optimized_sequence is not None

    def test_timeout_actually_fires(self):
        """Positive test: the SIGALRM timeout MUST fire when the
        optimizer is blocked in pure-Python code.

        Monkey-patches ``HybridOptimizer.optimize`` to insert a long
        ``time.sleep`` (pure-Python, so the signal is delivered at the
        next bytecode boundary).  With ``timeout_seconds=1`` and a
        10-second sleep, the timeout MUST fire after ~1 second and
        return a partial result with ``timed_out=True``.

        This is the test that actually exercises the timeout mechanism
        — the previous tests verify the optimizer doesn't hang, but
        this one verifies the timeout fires when it should.
        """
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer

        original_optimize = HybridOptimizer.optimize

        def slow_optimize(self, protein, *args, **kwargs):
            # Pure-Python sleep — SIGALRM is delivered here.
            time.sleep(10)
            return original_optimize(self, protein, *args, **kwargs)

        start = time.time()
        try:
            HybridOptimizer.optimize = slow_optimize
            result = optimize_sequence(
                "MVHLTPEEK",
                organism="human",
                strict_mode=False,
                timeout_seconds=1,
            )
        finally:
            HybridOptimizer.optimize = original_optimize

        elapsed = time.time() - start

        # The timeout MUST have fired (not completed normally).
        assert result.timed_out is True, (
            f"Expected timed_out=True but got timed_out={result.timed_out}; "
            f"elapsed={elapsed:.2f}s"
        )
        assert result.convergence_status == "timeout"
        # Elapsed should be slightly more than 1s (the timeout) and
        # much less than 10s (the injected sleep).
        assert elapsed < 5, (
            f"Timeout did not fire promptly (elapsed={elapsed:.2f}s); "
            f"the optimizer hung past the 1s budget"
        )
        # The partial result must still satisfy the length invariant.
        assert len(result.sequence) == len(result.protein) * 3
        # A warning must be recorded explaining the partial result.
        assert any("timed out" in w.lower() for w in result.warnings)

    def test_signal_handler_restored(self):
        """The SIGALRM handler must be restored to its prior state
        after ``optimize_sequence`` returns.

        Without this, a timeout installed by one call would leak into
        subsequent unrelated code (or be missing when a later call
        needs it).
        """
        before = signal.getsignal(signal.SIGALRM)
        try:
            optimize_sequence(
                "MVHLTPEEK",
                organism="human",
                strict_mode=False,
                timeout_seconds=5,
            )
        except Exception:
            pass
        after = signal.getsignal(signal.SIGALRM)
        assert before == after, (
            f"SIGALRM handler leaked: before={before!r}, after={after!r}"
        )

    def test_timeout_zero_disables_alarm(self):
        """``timeout_seconds=0`` must disable the alarm entirely
        (no SIGALRM installed, no handler swap).
        """
        # Should complete normally without installing any alarm.
        result = optimize_sequence(
            "MVHLTPEEK", organism="human", strict_mode=False, timeout_seconds=0
        )
        assert result.timed_out is False
        assert result.optimized_sequence is not None

    def test_default_timeout_constant_exists(self):
        """The ``OPTIMIZER_TIMEOUT_SECONDS`` constant must be exported
        from ``utils.py`` and have a sensible default value.
        """
        assert OPTIMIZER_TIMEOUT_SECONDS == 30
        # The default parameter on optimize_sequence should match.
        import inspect

        sig = inspect.signature(optimize_sequence)
        default = sig.parameters["timeout_seconds"].default
        assert default == OPTIMIZER_TIMEOUT_SECONDS

    def test_optimizer_timeout_exception_class(self):
        """The ``OptimizerTimeout`` exception class and
        ``_timeout_handler`` function must be importable from
        ``pipeline_core``.
        """
        assert issubclass(OptimizerTimeout, Exception)
        assert callable(_timeout_handler)
        # _timeout_handler must raise OptimizerTimeout when invoked.
        with pytest.raises(OptimizerTimeout):
            _timeout_handler(signal.SIGALRM, None)
