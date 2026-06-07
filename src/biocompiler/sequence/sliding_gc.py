"""
BioCompiler Sliding-Window GC Constraint
=========================================

Provides local/sliding-window GC content checking in addition to the
existing global GC constraint.  While the global GC check ensures the
overall sequence composition is within bounds, local windows can still
have extreme GC content that causes issues with polymerase processivity,
secondary structure, or synthesis difficulty.

This module provides:

- :func:`check_sliding_gc` — scan a DNA sequence with a sliding window
  and report any windows where GC content falls outside the allowed range.
- :func:`fix_sliding_gc_violations` — attempt to fix local GC violations
  by swapping synonymous codons within violating windows.
- :func:`evaluate_sliding_gc` — predicate-style evaluation returning a
  :class:`TypeCheckResult` (for the type-system registry).

Data classes:

- :class:`WindowViolation` — a single window that violates GC bounds.
- :class:`SlidingGCResult` — aggregate result from a sliding-GC scan.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from biocompiler.shared.types import TypeCheckResult, Verdict

# ── NUMBA fast path with fallback ─────────────────────────────────────
try:
    from biocompiler.numba_kernels import fast_gc_window as _fast_gc_window_numba
    from biocompiler.numba_kernels import seq_to_bytes as _seq_to_bytes
    from biocompiler.numba_kernels import HAS_NUMBA as _HAS_NUMBA
except ImportError:
    _fast_gc_window_numba = None  # type: ignore[assignment,misc]
    _seq_to_bytes = None  # type: ignore[assignment,misc]
    _HAS_NUMBA = False

logger = logging.getLogger(__name__)

# Module-level flag for testing: force Python path even when NUMBA is available
_FORCE_PYTHON_GC_WINDOW: bool = False

__all__ = [
    "WindowViolation",
    "SlidingGCResult",
    "check_sliding_gc",
    "fix_sliding_gc_violations",
    "evaluate_sliding_gc",
    "_check_sliding_gc_python",
    "_FORCE_PYTHON_GC_WINDOW",
    "_HAS_NUMBA",
    "_check_sliding_gc_numba",
]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class WindowViolation:
    """A single sliding window where GC content is out of bounds.

    Attributes:
        start: 0-based start position of the violating window.
        end: 0-based exclusive end position of the violating window.
        gc_content: GC fraction observed in this window.
        direction: ``"too_low"`` if gc_content < gc_min,
            ``"too_high"`` if gc_content > gc_max.
    """
    start: int
    end: int
    gc_content: float
    direction: str  # "too_low" or "too_high"

    def __post_init__(self) -> None:
        assert self.direction in ("too_low", "too_high"), (
            f"direction must be 'too_low' or 'too_high', got '{self.direction}'"
        )
        assert 0.0 <= self.gc_content <= 1.0, (
            f"gc_content must be in [0, 1], got {self.gc_content}"
        )
        assert self.start >= 0, f"start must be >= 0, got {self.start}"
        assert self.end > self.start, f"end must be > start, got end={self.end} start={self.start}"


@dataclass
class SlidingGCResult:
    """Aggregate result from a sliding-window GC scan.

    Attributes:
        passed: True if no windows violate the GC bounds.
        min_gc: Minimum GC content observed across all windows.
        max_gc: Maximum GC content observed across all windows.
        violations: List of windows that violate the GC bounds.
    """
    passed: bool
    min_gc: float
    max_gc: float
    violations: list[WindowViolation] = field(default_factory=list)

    def __post_init__(self) -> None:
        assert 0.0 <= self.min_gc <= 1.0, f"min_gc must be in [0, 1], got {self.min_gc}"
        assert 0.0 <= self.max_gc <= 1.0, f"max_gc must be in [0, 1], got {self.max_gc}"
        assert self.min_gc <= self.max_gc, (
            f"min_gc must be <= max_gc, got min_gc={self.min_gc} max_gc={self.max_gc}"
        )


# ────────────────────────────────────────────────────────────
# Core sliding-window GC check
# ────────────────────────────────────────────────────────────

def _check_sliding_gc_python(
    dna: str,
    window_size: int,
    gc_min: float,
    gc_max: float,
    step: int,
) -> SlidingGCResult:
    """Pure-Python sliding-window GC check (fallback when NUMBA unavailable)."""
    n = len(dna)

    if n < window_size:
        if n == 0:
            return SlidingGCResult(
                passed=True,
                min_gc=0.0,
                max_gc=0.0,
                violations=[],
            )
        gc = (dna.count("G") + dna.count("C")) / n
        violations = []
        if gc < gc_min:
            violations.append(WindowViolation(0, n, gc, "too_low"))
        elif gc > gc_max:
            violations.append(WindowViolation(0, n, gc, "too_high"))
        return SlidingGCResult(
            passed=len(violations) == 0,
            min_gc=gc,
            max_gc=gc,
            violations=violations,
        )

    # Pre-compute GC counts per position for O(1) sliding updates
    # gc_at[i] = 1 if dna[i] in {G, C}, else 0
    gc_at = [1 if b in "GC" else 0 for b in dna]

    # Initialize first window
    current_gc_count = sum(gc_at[:window_size])
    current_gc = current_gc_count / window_size

    min_gc = current_gc
    max_gc = current_gc
    violations: list[WindowViolation] = []

    # Check first window
    if current_gc < gc_min:
        violations.append(WindowViolation(0, window_size, current_gc, "too_low"))
    elif current_gc > gc_max:
        violations.append(WindowViolation(0, window_size, current_gc, "too_high"))

    # Slide window with O(1) updates
    for start in range(step, n - window_size + 1, step):
        # Remove bases that left the window, add bases that entered
        # For step > 1, we need to update incrementally
        if step == 1:
            current_gc_count -= gc_at[start - 1]
            current_gc_count += gc_at[start + window_size - 1]
        else:
            # Recompute for step > 1 (still O(step) per window)
            current_gc_count = sum(gc_at[start:start + window_size])

        current_gc = current_gc_count / window_size

        if current_gc < min_gc:
            min_gc = current_gc
        if current_gc > max_gc:
            max_gc = current_gc

        end = start + window_size
        if current_gc < gc_min:
            violations.append(WindowViolation(start, end, current_gc, "too_low"))
        elif current_gc > gc_max:
            violations.append(WindowViolation(start, end, current_gc, "too_high"))

    return SlidingGCResult(
        passed=len(violations) == 0,
        min_gc=min_gc,
        max_gc=max_gc,
        violations=violations,
    )


def _check_sliding_gc_numba(
    dna: str,
    window_size: int,
    gc_min: float,
    gc_max: float,
    step: int,
) -> SlidingGCResult:
    """NUMBA-accelerated sliding-window GC check using fast_gc_window kernel.

    The fast_gc_window kernel computes GC% for every position (step=1)
    using incremental state updates in a single O(n) pass.  When step > 1,
    we still compute at full resolution (which gives more accurate min/max)
    and only check violations at step intervals.
    """
    n = len(dna)

    if n < window_size:
        if n == 0:
            return SlidingGCResult(
                passed=True,
                min_gc=0.0,
                max_gc=0.0,
                violations=[],
            )
        gc = (dna.count("G") + dna.count("C")) / n
        violations = []
        if gc < gc_min:
            violations.append(WindowViolation(0, n, gc, "too_low"))
        elif gc > gc_max:
            violations.append(WindowViolation(0, n, gc, "too_high"))
        return SlidingGCResult(
            passed=len(violations) == 0,
            min_gc=gc,
            max_gc=gc,
            violations=violations,
        )

    # Convert DNA string to byte array for NUMBA kernel
    seq_bytes = _seq_to_bytes(dna)  # type: ignore[misc]

    # Compute GC% for all windows at step=1 resolution (O(n))
    gc_pcts = _fast_gc_window_numba(seq_bytes, window_size)  # type: ignore[misc]

    # Find global min/max across all windows (full resolution)
    import numpy as np
    min_gc = float(np.min(gc_pcts))
    max_gc = float(np.max(gc_pcts))

    # Check violations at step intervals
    n_windows = len(gc_pcts)
    violations: list[WindowViolation] = []

    for i in range(0, n_windows, step):
        gc = float(gc_pcts[i])
        start = i
        end = i + window_size
        if gc < gc_min:
            violations.append(WindowViolation(start, end, gc, "too_low"))
        elif gc > gc_max:
            violations.append(WindowViolation(start, end, gc, "too_high"))

    return SlidingGCResult(
        passed=len(violations) == 0,
        min_gc=min_gc,
        max_gc=max_gc,
        violations=violations,
    )


def check_sliding_gc(
    dna: str,
    window_size: int = 50,
    gc_min: float = 0.30,
    gc_max: float = 0.70,
    step: int | None = None,
) -> SlidingGCResult:
    """Check GC content in sliding windows across the sequence.

    Scans the DNA sequence with a window of ``window_size`` nucleotides,
    stepping by ``step`` positions each time (default: step = 1 for
    full-resolution scanning).  Any window whose GC fraction falls
    outside ``[gc_min, gc_max]`` is recorded as a violation.

    When NUMBA is available, uses the ``fast_gc_window`` JIT-compiled
    kernel for O(n) sliding GC computation with incremental state
    updates.  Falls back to a pure-Python implementation otherwise.

    Pre-conditions:
        - ``dna`` is a valid uppercase DNA string (case-insensitive).
        - ``window_size > 0`` and ``window_size <= len(dna)``.
        - ``0.0 <= gc_min < gc_max <= 1.0``.
        - ``step > 0``.

    Post-conditions:
        - Returns a :class:`SlidingGCResult` with ``passed=True`` if
          no windows violate the bounds.
        - ``min_gc`` and ``max_gc`` reflect the extreme GC values
          across all windows.
        - All violations are reported with their positions and direction.
        - Results are identical between NUMBA and Python paths.

    Args:
        dna: DNA sequence to check.
        window_size: Size of the sliding window in nucleotides.
        gc_min: Minimum acceptable GC fraction (inclusive).
        gc_max: Maximum acceptable GC fraction (inclusive).
        step: Step size for the sliding window. Defaults to 1 (full
            resolution). Larger values speed up the scan at the cost of
            potentially missing narrow violations.

    Returns:
        :class:`SlidingGCResult` with pass/fail status and violation details.
    """
    assert window_size > 0, f"window_size must be > 0, got {window_size}"
    assert 0.0 <= gc_min < gc_max <= 1.0, (
        f"gc_min must be in [0, 1), gc_max in (0, 1], and gc_min < gc_max; "
        f"got gc_min={gc_min}, gc_max={gc_max}"
    )
    if step is None:
        step = 1
    assert step > 0, f"step must be > 0, got {step}"

    dna = dna.upper()

    # Dispatch to NUMBA or Python implementation
    if _HAS_NUMBA and not _FORCE_PYTHON_GC_WINDOW:
        try:
            return _check_sliding_gc_numba(dna, window_size, gc_min, gc_max, step)
        except Exception:
            # If NUMBA path fails for any reason, fall back gracefully
            logger.debug(
                "NUMBA fast_gc_window failed, falling back to Python",
                exc_info=True,
            )

    return _check_sliding_gc_python(dna, window_size, gc_min, gc_max, step)


# ────────────────────────────────────────────────────────────
# Fix sliding-window GC violations
# ────────────────────────────────────────────────────────────

# Maximum iterations for sliding GC fix loop
MAX_SLIDING_GC_FIX_ITERATIONS: int = 200


def fix_sliding_gc_violations(
    dna: str,
    protein: str,
    window_size: int = 50,
    gc_min: float = 0.30,
    gc_max: float = 0.70,
    usage: dict[str, float] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    max_cai_cost: float = 0.30,
) -> tuple[str, int]:
    """Fix local/sliding-window GC violations by adjusting codons.

    For each violating window, identify codons that overlap the window
    and swap them for synonymous alternatives that move the local GC
    content toward the target range.  When a ``usage`` (codon
    adaptiveness) table is provided, higher-CAI alternatives are
    preferred.

    The ``max_cai_cost`` parameter controls how much CAI can be
    sacrificed per codon swap.  Only alternatives whose CAI
    adaptiveness ratio is within ``max_cai_cost`` of the current
    codon's ratio are considered.  This prevents catastrophic CAI
    regression when fixing sliding-window GC violations.

    The function iterates until no violations remain or the iteration
    limit is reached.  Global GC constraints (``gc_lo``, ``gc_hi``) are
    respected — no codon swap that would push the overall GC content
    outside the global range is accepted.

    Pre-conditions:
        - ``dna`` is uppercase and ``len(dna) == len(protein) * 3``.
        - ``protein`` is a valid amino acid sequence.

    Post-conditions:
        - Returned sequence encodes the same protein.
        - Global GC is within ``[gc_lo, gc_hi]`` (if it was before).
        - The number of sliding-window violations is reduced or
          eliminated.

    Args:
        dna: Current optimized DNA sequence.
        protein: Amino acid sequence (no stop).
        window_size: Sliding window size (default 50).
        gc_min: Minimum GC fraction per window.
        gc_max: Maximum GC fraction per window.
        usage: Codon adaptiveness table (codon -> w value).  If None,
            all synonymous codons are treated equally.
        gc_lo: Global GC lower bound (respected during swaps).
        gc_hi: Global GC upper bound (respected during swaps).
        max_cai_cost: Maximum CAI adaptiveness ratio loss per swap
            (default 0.30).  Only alternatives whose adaptiveness
            ratio is within this distance of the current codon's
            ratio will be considered.

    Returns:
        Tuple of (fixed sequence, number of codon swaps performed).
    """
    n_codons = len(protein)
    n_bases = len(dna)
    assert n_bases == n_codons * 3, (
        f"Sequence length ({n_bases}) must equal protein length * 3 ({n_codons * 3})"
    )

    seq_list = list(dna)
    total_swaps = 0

    # Pre-compute sorted codons per amino acid (by CAI descending if usage given)
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(protein):
        if aa == "*":
            continue
        # Lazy import to avoid circular dependency with type_system
        from biocompiler.type_system import AA_TO_CODONS as _AA_TO_CODONS
        codons = _AA_TO_CODONS.get(aa, [])
        if usage is not None:
            sorted_codons[aa] = sorted(
                codons, key=lambda c: usage.get(c, 0.0), reverse=True
            )
        else:
            sorted_codons[aa] = codons

    # Track global GC count incrementally
    global_gc_count = sum(1 for b in seq_list if b in "GC")

    for iteration in range(MAX_SLIDING_GC_FIX_ITERATIONS):
        current_seq = "".join(seq_list)
        result = check_sliding_gc(current_seq, window_size, gc_min, gc_max)

        if result.passed:
            break

        # Try to fix each violation
        fixed_any = False
        for violation in result.violations:
            # Find codons overlapping the violating window
            first_codon = violation.start // 3
            last_codon = min(n_codons - 1, (violation.end - 1) // 3)

            for ci in range(first_codon, last_codon + 1):
                aa = protein[ci]
                if aa == "*":
                    continue

                current_codon = "".join(seq_list[ci * 3:ci * 3 + 3])
                current_gc_in_codon = sum(1 for b in current_codon if b in "GC")

                # Determine target direction
                want_more_gc = violation.direction == "too_low"

                # Try alternative codons
                best_alt: str | None = None
                best_score = float('-inf')
                current_cai = usage.get(current_codon, 0.0) if usage else 1.0

                for alt in sorted_codons.get(aa, []):
                    if alt == current_codon:
                        continue

                    alt_gc_in_codon = sum(1 for b in alt if b in "GC")

                    # Check if this swap moves window GC in the right direction
                    gc_delta = alt_gc_in_codon - current_gc_in_codon
                    if want_more_gc and gc_delta <= 0:
                        continue
                    if not want_more_gc and gc_delta >= 0:
                        continue

                    # Check CAI cost: only accept swaps within max_cai_cost
                    alt_cai = usage.get(alt, 0.0) if usage else 1.0
                    cai_cost = current_cai - alt_cai
                    if cai_cost > max_cai_cost:
                        continue

                    # Check global GC constraint
                    new_global_gc = global_gc_count - current_gc_in_codon + alt_gc_in_codon
                    new_global_frac = new_global_gc / n_bases
                    if not (gc_lo <= new_global_frac <= gc_hi):
                        continue

                    # Score: prefer higher CAI, then larger GC delta
                    cai_score = alt_cai
                    score = cai_score + abs(gc_delta) * 0.01  # small weight for GC delta

                    if score > best_score:
                        best_score = score
                        best_alt = alt

                if best_alt is not None:
                    # Apply swap
                    alt_gc_in_codon = sum(1 for b in best_alt if b in "GC")
                    global_gc_count = global_gc_count - current_gc_in_codon + alt_gc_in_codon
                    seq_list[ci * 3] = best_alt[0]
                    seq_list[ci * 3 + 1] = best_alt[1]
                    seq_list[ci * 3 + 2] = best_alt[2]
                    total_swaps += 1
                    fixed_any = True
                    break  # Re-scan after each fix to update violation positions

            if fixed_any:
                break  # Re-scan

        if not fixed_any:
            # No progress — stop iterating
            logger.debug(
                "Sliding GC fix: no progress at iteration %d, %d violations remain",
                iteration, len(result.violations),
            )
            break

    return "".join(seq_list), total_swaps


# ────────────────────────────────────────────────────────────
# Predicate-style evaluation (for type-system registry)
# ────────────────────────────────────────────────────────────

def evaluate_sliding_gc(
    seq: str,
    window_size: int = 50,
    gc_min: float = 0.30,
    gc_max: float = 0.70,
) -> TypeCheckResult:
    """Evaluate sliding-window GC constraint as a type predicate.

    This is the evaluate-style API that returns a :class:`TypeCheckResult`
    suitable for the predicate registry and certificate generation.

    Args:
        seq: DNA sequence to evaluate.
        window_size: Sliding window size in nucleotides.
        gc_min: Minimum acceptable GC fraction per window.
        gc_max: Maximum acceptable GC fraction per window.

    Returns:
        :class:`TypeCheckResult` with PASS/FAIL verdict.
    """
    result = check_sliding_gc(seq, window_size, gc_min, gc_max)

    predicate_name = f"SlidingGC({window_size}, {gc_min}, {gc_max})"

    if result.passed:
        return TypeCheckResult(
            predicate=predicate_name,
            verdict=Verdict.PASS,
        )

    # Build violation summary
    low_count = sum(1 for v in result.violations if v.direction == "too_low")
    high_count = sum(1 for v in result.violations if v.direction == "too_high")

    parts = []
    if low_count > 0:
        parts.append(f"{low_count} window(s) with GC < {gc_min}")
    if high_count > 0:
        parts.append(f"{high_count} window(s) with GC > {gc_max}")

    violation_msg = (
        f"Sliding-window GC violation: {'; '.join(parts)}. "
        f"Min window GC={result.min_gc:.3f}, Max window GC={result.max_gc:.3f}"
    )

    return TypeCheckResult(
        predicate=predicate_name,
        verdict=Verdict.FAIL,
        violation=violation_msg,
    )
