"""
BioCompiler NUMBA JIT-Compiled Kernels
=======================================

NUMBA-accelerated inner loops for codon optimization hot paths:

Existing kernels:
1. count_gc                — Count G and C characters in a DNA sequence
2. count_dinucleotides     — Count occurrences of a dinucleotide (GT, CG, AG)
3. compute_cai_kernel      — Compute CAI from pre-indexed adaptiveness values (log-sum-exp)
4. scan_restriction_sites  — Find all positions of a restriction site pattern
5. find_all_dinucleotide_positions — Find all positions of a dinucleotide

New kernels (v2 — incremental / batch):
6. compute_cai_incremental — Update CAI after a single codon swap without full recompute
7. batch_codon_swap_score  — Score all possible codon swaps at one position (vectorized)
8. fast_gc_window          — Sliding-window GC% with incremental state updates
9. fast_dinucleotide_count — Count multiple dinucleotides in a single pass

New kernels (v3 — mRNA structure / premature stops):
10. _numba_count_dinucleotides — Count a 2-char pattern (CpG counting primitive)
11. _numba_has_premature_stop  — Scan codons 0..n-2 for TAA/TAG/TGA stop codons
12. _numba_compute_approx_dg   — Nussinov-like hairpin ΔG (eliminates 125K+ helper calls)

All kernels use @numba.njit(cache=True, nogil=True) for maximum performance
with cached compilation and thread safety. Pure-Python fallbacks are provided
when NUMBA is unavailable.

Architecture:
    - Input: byte arrays (np.frombuffer or np.array) for zero-copy NUMBA access
    - Output: scalar int/float or numpy arrays of positions
    - All kernels are stateless and side-effect free
    - Compilation cache avoids re-JIT overhead on repeated imports
    - prange used for data-parallel loops where dependencies allow

Usage:
    from biocompiler.optimizer.numba_kernels import HAS_NUMBA, count_gc, count_dinucleotides

    if HAS_NUMBA:
        seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        gc = count_gc(seq_bytes)
        gt_count = count_dinucleotides(seq_bytes, b'GT')
"""

from __future__ import annotations

import math
from typing import List

# ── NUMBA availability detection ──────────────────────────────────────
try:
    import numba
    import numpy as np
    from numba import njit, prange

    HAS_NUMBA: bool = True
    _NUMBA_VERSION: str = numba.__version__
except ImportError:
    HAS_NUMBA: bool = False
    _NUMBA_VERSION = ""
    numba = None  # type: ignore[assignment]
    prange = range  # type: ignore[assignment]

# Runtime toggle: set to False to disable NUMBA even when available
import os as _os
USE_NUMBA: bool = HAS_NUMBA and _os.environ.get("BIOCOMPILER_USE_NUMBA", "1") not in ("0", "false", "False")

__all__ = [
    "HAS_NUMBA",
    "USE_NUMBA",
    "count_gc",
    "count_dinucleotides",
    "compute_cai_kernel",
    "scan_restriction_sites",
    "find_all_dinucleotide_positions",
    "seq_to_bytes",
    # New v2 kernels
    "compute_cai_incremental",
    "batch_codon_swap_score",
    "fast_gc_window",
    "fast_dinucleotide_count",
    # Additional optimized utility kernels
    "count_gc_parallel",
    "scan_restriction_sites_multi",
    # v3 kernels (mRNA ΔG, premature stops, dinucleotide count)
    "_numba_count_dinucleotides",
    "_numba_has_premature_stop",
    "_numba_compute_approx_dg",
]

# ── Byte constants for DNA characters ─────────────────────────────────
_G_BYTE = ord('G')
_C_BYTE = ord('C')
_A_BYTE = ord('A')
_T_BYTE = ord('T')


# ══════════════════════════════════════════════════════════════════════
# Helper: convert DNA string to numpy byte array for zero-copy access
# ══════════════════════════════════════════════════════════════════════

def seq_to_bytes(seq: str) -> "np.ndarray":
    """Convert a DNA string to a numpy uint8 array for NUMBA kernel access.

    This is a zero-copy operation when using np.frombuffer on an ASCII
    encoded string. The resulting array can be passed directly to NUMBA
    njit kernels.

    Args:
        seq: DNA sequence string (uppercase ACGT).

    Returns:
        numpy uint8 array of byte values.
    """
    import numpy as np
    return np.frombuffer(seq.encode('ascii'), dtype=np.uint8)


# ══════════════════════════════════════════════════════════════════════
# NUMBA JIT-compiled kernels
# ══════════════════════════════════════════════════════════════════════

if HAS_NUMBA:
    import numpy as np

    # ── Existing kernels (optimized) ──────────────────────────────────

    @njit(cache=True, nogil=True)
    def count_gc(seq_bytes: np.ndarray) -> int:
        """Count G and C characters in a DNA byte sequence.

        Uses a simple byte comparison loop, which NUMBA compiles to
        tight machine code with branch prediction.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.

        Returns:
            Number of G and C characters.
        """
        count = 0
        for i in range(len(seq_bytes)):
            b = seq_bytes[i]
            if b == _G_BYTE or b == _C_BYTE:
                count += 1
        return count

    @njit(cache=True, nogil=True)
    def count_dinucleotides(seq_bytes: np.ndarray, dinuc_bytes: np.ndarray) -> int:
        """Count occurrences of a dinucleotide in a DNA byte sequence.

        Scans the sequence for the 2-byte pattern. The dinuc_bytes
        must be exactly 2 bytes long.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.
            dinuc_bytes: numpy uint8 array of exactly 2 bytes (e.g. b'GT').

        Returns:
            Number of occurrences of the dinucleotide.
        """
        count = 0
        b0 = dinuc_bytes[0]
        b1 = dinuc_bytes[1]
        n = len(seq_bytes)
        for i in range(n - 1):
            if seq_bytes[i] == b0 and seq_bytes[i + 1] == b1:
                count += 1
        return count

    @njit(cache=True, nogil=True)
    def compute_cai_kernel(
        adaptiveness_array: np.ndarray,
        codon_indices: np.ndarray,
        n_codons: int,
    ) -> float:
        """Compute CAI from pre-indexed adaptiveness values.

        Uses log-sum for numerical stability (geometric mean of
        relative adaptiveness values). This replaces the pure-Python
        loop in compute_cai() for the hottest part.

        The caller must pre-compute:
        - adaptiveness_array: flat array of adaptiveness values indexed by codon
        - codon_indices: array mapping each codon position to its adaptiveness index
        - n_codons: number of codons contributing to CAI (excludes Met and stop)

        The geometric mean is computed as:
            CAI = exp(sum(log(w_i)) / n)

        where w_i is the adaptiveness value for codon i.

        Args:
            adaptiveness_array: numpy float64 array of adaptiveness values.
            codon_indices: numpy int64 array of codon adaptiveness indices.
            n_codons: number of codons to include in the CAI computation.

        Returns:
            CAI value (float64). Returns 0.0 if n_codons == 0.
        """
        if n_codons == 0:
            return 0.0

        log_sum = 0.0
        epsilon = 1e-10
        for i in range(n_codons):
            idx = codon_indices[i]
            w = adaptiveness_array[idx]
            if w <= 0.0:
                w = epsilon
            log_sum += math.log(w)

        return math.exp(log_sum / n_codons)

    @njit(cache=True, nogil=True)
    def scan_restriction_sites(
        seq_bytes: np.ndarray,
        pattern_bytes: np.ndarray,
        pattern_len: int,
    ) -> np.ndarray:
        """Find all positions of a restriction site pattern in a DNA sequence.

        Simple sliding-window pattern matching. For multi-pattern
        matching, use AhoCorasickScanner instead. This kernel is
        optimized for single-pattern scans where the pattern is short
        (4-8 bases typical for restriction enzymes).

        Returns a numpy int64 array of positions where the pattern
        starts. The array length equals the number of matches found.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.
            pattern_bytes: numpy uint8 array of the pattern.
            pattern_len: length of the pattern.

        Returns:
            numpy int64 array of match positions (may be empty).
        """
        n = len(seq_bytes)
        max_matches = n  # worst case: every position matches
        # Pre-allocate result buffer
        results = np.empty(max_matches, dtype=np.int64)
        count = 0

        if pattern_len > n:
            # Return empty array
            return results[:0].copy()

        for i in range(n - pattern_len + 1):
            match = True
            for j in range(pattern_len):
                if seq_bytes[i + j] != pattern_bytes[j]:
                    match = False
                    break
            if match:
                results[count] = i
                count += 1

        return results[:count].copy()

    @njit(cache=True, nogil=True)
    def find_all_dinucleotide_positions(
        seq_bytes: np.ndarray,
        dinuc_bytes: np.ndarray,
    ) -> np.ndarray:
        """Find all positions of a dinucleotide in a DNA sequence.

        Returns a numpy int64 array of positions where the dinucleotide
        starts. This is the NUMBA equivalent of scanning for GT, CG, AG
        positions in the incremental sequence state.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.
            dinuc_bytes: numpy uint8 array of exactly 2 bytes.

        Returns:
            numpy int64 array of match positions (may be empty).
        """
        n = len(seq_bytes)
        results = np.empty(n, dtype=np.int64)
        count = 0
        b0 = dinuc_bytes[0]
        b1 = dinuc_bytes[1]

        for i in range(n - 1):
            if seq_bytes[i] == b0 and seq_bytes[i + 1] == b1:
                results[count] = i
                count += 1

        return results[:count].copy()

    # ── New v2 kernels ────────────────────────────────────────────────

    @njit(cache=True, nogil=True)
    def compute_cai_incremental(
        current_log_sum: float,
        n_codons: int,
        old_adaptiveness: float,
        new_adaptiveness: float,
    ) -> float:
        """Update CAI after a single codon swap without full recompute.

        Instead of recomputing the entire log-sum from scratch, this
        kernel adjusts the running sum by subtracting the old codon's
        log-adaptiveness and adding the new one. This is O(1) instead
        of O(n_codons).

        The geometric mean is:
            CAI = exp(log_sum / n_codons)

        After swapping codon at position i:
            new_log_sum = old_log_sum - log(w_old) + log(w_new)

        Args:
            current_log_sum: The current sum of log(w_i) values.
            n_codons: Number of codons in the CAI computation (unchanged).
            old_adaptiveness: Adaptiveness value of the old codon being replaced.
            new_adaptiveness: Adaptiveness value of the new codon.

        Returns:
            Updated CAI value (float64). Returns 0.0 if n_codons == 0.
        """
        # Edge cases: empty sequence or invalid codon count
        if n_codons <= 0:
            return 0.0

        epsilon = 1e-10

        w_old = old_adaptiveness
        if w_old <= 0.0 or w_old != w_old:  # w_old != w_old catches NaN
            w_old = epsilon

        w_new = new_adaptiveness
        if w_new <= 0.0 or w_new != w_new:  # w_new != w_new catches NaN
            w_new = epsilon

        # Single codon: CAI is simply the adaptiveness of that codon
        # (after the swap).  The general formula still works for n_codons==1
        # but we make it explicit for clarity.
        new_log_sum = current_log_sum - math.log(w_old) + math.log(w_new)
        return math.exp(new_log_sum / n_codons)

    @njit(cache=True, nogil=True)
    def batch_codon_swap_score(
        adaptiveness_array: np.ndarray,
        codon_indices: np.ndarray,
        n_codons: int,
        swap_position: int,
        candidate_indices: np.ndarray,
        n_candidates: int,
        current_log_sum: float,
    ) -> np.ndarray:
        """Score all possible codon swaps at a single position.

        Given the current CAI log-sum and the index of the codon being
        swapped, compute the resulting CAI for each candidate replacement
        codon in a single vectorized pass.

        This avoids calling compute_cai_incremental N times from Python,
        eliminating per-call overhead and enabling SIMD-friendly access
        patterns.

        Args:
            adaptiveness_array: numpy float64 array of adaptiveness values.
            codon_indices: numpy int64 array of codon adaptiveness indices.
            n_codons: Total number of codons in CAI computation.
            swap_position: Index into codon_indices of the codon being swapped.
            candidate_indices: numpy int64 array of adaptiveness indices for
                candidate replacement codons.
            n_candidates: Number of candidate codons.
            current_log_sum: Current sum of log(w_i) for all codons.

        Returns:
            numpy float64 array of CAI scores, one per candidate.
            Length equals n_candidates.
        """
        # No alternative codons (e.g. Met/Trp): return empty scores array
        if n_candidates == 0:
            return np.empty(0, dtype=np.float64)

        scores = np.empty(n_candidates, dtype=np.float64)
        epsilon = 1e-10

        if n_codons <= 0:
            for k in range(n_candidates):
                scores[k] = 0.0
            return scores

        # Get old codon's adaptiveness
        old_idx = codon_indices[swap_position]
        w_old = adaptiveness_array[old_idx]
        if w_old <= 0.0 or w_old != w_old:  # catch NaN
            w_old = epsilon
        log_w_old = math.log(w_old)

        for k in range(n_candidates):
            w_new = adaptiveness_array[candidate_indices[k]]
            if w_new <= 0.0 or w_new != w_new:  # catch NaN
                w_new = epsilon
            new_log_sum = current_log_sum - log_w_old + math.log(w_new)
            scores[k] = math.exp(new_log_sum / n_codons)

        return scores

    @njit(cache=True, nogil=True, parallel=True)
    def fast_gc_window(
        seq_bytes: np.ndarray,
        window_size: int,
    ) -> np.ndarray:
        """Sliding-window GC% computation with incremental state updates.

        Computes GC content for each window of `window_size` bases using
        an O(n) incremental approach: maintain a running count of G+C in
        the window, subtract the outgoing base and add the incoming base
        as the window slides.

        The parallel=True flag enables prange for the initial GC count
        computation, while the sliding window itself is sequential due
        to data dependencies.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.
            window_size: Size of the sliding window in bases.

        Returns:
            numpy float64 array of GC% values for each window position.
            Length is max(len(seq_bytes) - window_size + 1, 0).
            Returns empty array if sequence is shorter than window_size.
        """
        n = len(seq_bytes)
        n_windows = n - window_size + 1

        if n_windows <= 0:
            return np.empty(0, dtype=np.float64)

        results = np.empty(n_windows, dtype=np.float64)

        # Compute GC count for the first window using prange
        gc_count = 0
        for i in prange(window_size):
            b = seq_bytes[i]
            if b == _G_BYTE or b == _C_BYTE:
                gc_count += 1

        results[0] = gc_count / window_size

        # Slide the window: subtract outgoing, add incoming
        for i in range(1, n_windows):
            outgoing = seq_bytes[i - 1]
            incoming = seq_bytes[i + window_size - 1]

            if outgoing == _G_BYTE or outgoing == _C_BYTE:
                gc_count -= 1
            if incoming == _G_BYTE or incoming == _C_BYTE:
                gc_count += 1

            results[i] = gc_count / window_size

        return results

    @njit(cache=True, nogil=True)
    def fast_dinucleotide_count(
        seq_bytes: np.ndarray,
        dinuc_keys: np.ndarray,
        n_dinucs: int,
    ) -> np.ndarray:
        """Count multiple dinucleotides in a single pass over the sequence.

        Instead of calling count_dinucleotides() separately for each
        dinucleotide (which scans the entire sequence each time), this
        kernel scans the sequence once and counts all requested
        dinucleotides simultaneously.

        Each dinucleotide is encoded as a pair of bytes in dinuc_keys:
            dinuc_keys[i, 0] = first byte of dinucleotide i
            dinuc_keys[i, 1] = second byte of dinucleotide i

        This is O(n * d) where d is the number of dinucleotides, but
        with far better cache behavior than d separate O(n) scans.

        Common use case: counting CpG (CG), GT, and AG dinucleotides
        for splice site and CpG island checks.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.
            dinuc_keys: numpy uint8 2D array of shape (n_dinucs, 2),
                where each row is a dinucleotide's two bytes.
            n_dinucs: Number of dinucleotides to count.

        Returns:
            numpy int64 array of counts, one per dinucleotide.
            Length equals n_dinucs.
        """
        counts = np.zeros(n_dinucs, dtype=np.int64)
        n = len(seq_bytes)

        # Sequence shorter than 2 bases cannot contain any dinucleotide
        if n < 2:
            return counts

        for i in range(n - 1):
            b0 = seq_bytes[i]
            b1 = seq_bytes[i + 1]
            for d in range(n_dinucs):
                if b0 == dinuc_keys[d, 0] and b1 == dinuc_keys[d, 1]:
                    counts[d] += 1

        return counts

    # ── Additional optimized utility kernel ────────────────────────────

    @njit(cache=True, nogil=True, parallel=True)
    def count_gc_parallel(seq_bytes: np.ndarray) -> int:
        """Count G and C characters in a DNA byte sequence using parallel reduction.

        This is the parallel version of count_gc, suitable for very long
        sequences (e.g., whole chromosomes). For sequences shorter than
        ~10KB, the single-threaded count_gc is faster due to threading
        overhead.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.

        Returns:
            Number of G and C characters.
        """
        n = len(seq_bytes)
        count = 0
        # prange with += reduction is supported by numba
        for i in prange(n):
            b = seq_bytes[i]
            if b == _G_BYTE or b == _C_BYTE:
                count += 1
        return count

    @njit(cache=True, nogil=True, parallel=True)
    def scan_restriction_sites_multi(
        seq_bytes: np.ndarray,
        pattern_bytes: np.ndarray,
        pattern_offsets: np.ndarray,
        pattern_lens: np.ndarray,
        n_patterns: int,
    ) -> np.ndarray:
        """Find positions of multiple restriction site patterns in parallel.

        Each pattern is scanned independently, making this suitable for
        prange parallelism. Results for all patterns are concatenated
        with pattern_id encoded as (position << 16) | pattern_id to
        allow the caller to separate them.

        Args:
            seq_bytes: numpy uint8 array of DNA characters.
            pattern_bytes: numpy uint8 flat array of all pattern bytes concatenated.
            pattern_offsets: numpy int64 array of start offsets for each pattern
                in pattern_bytes.
            pattern_lens: numpy int64 array of lengths for each pattern.
            n_patterns: Number of patterns to scan.

        Returns:
            numpy int64 array of encoded matches. Each match is encoded as
            (position << 16) | pattern_id. Caller should decode to separate
            positions from pattern IDs. Array may be larger than actual
            matches; trailing entries should be ignored up to returned count.
        """
        n = len(seq_bytes)
        max_total = n * n_patterns  # worst case
        # We'll use a simpler approach: count per pattern, then merge
        # For efficiency, pre-allocate per-pattern counts
        match_counts = np.zeros(n_patterns, dtype=np.int64)

        for p in prange(n_patterns):
            plen = pattern_lens[p]
            poff = pattern_offsets[p]
            cnt = 0
            for i in range(n - plen + 1):
                match = True
                for j in range(plen):
                    if seq_bytes[i + j] != pattern_bytes[poff + j]:
                        match = False
                        break
                if match:
                    cnt += 1
            match_counts[p] = cnt

        # Total matches
        total = 0
        for p in range(n_patterns):
            total += match_counts[p]

        results = np.empty(total, dtype=np.int64)
        # Second pass: fill results
        # (Cannot easily parallelize the fill due to write ordering,
        #  but the expensive pattern matching was parallelized above)
        offsets = np.zeros(n_patterns, dtype=np.int64)
        for p in range(1, n_patterns):
            offsets[p] = offsets[p - 1] + match_counts[p - 1]

        for p in range(n_patterns):
            plen = pattern_lens[p]
            poff = pattern_offsets[p]
            idx = offsets[p]
            for i in range(n - plen + 1):
                match = True
                for j in range(plen):
                    if seq_bytes[i + j] != pattern_bytes[poff + j]:
                        match = False
                        break
                if match:
                    results[idx] = (i << 16) | p
                    idx += 1

        return results

    # ── New v3 kernels (mRNA ΔG, premature stops, dinucleotide count) ─
    # These target the top 3 bottlenecks in the optimization profile:
    #   1. compute_approx_dg (0.042s self + 0.019s helpers = 0.061s)
    #   2. _has_premature_stop (0.013s, 312 calls)
    #   3. CpG dinucleotide counting (called 1642 times)
    #
    # All three are pure-Python hot loops that benefit from @njit.

    # Byte constants for nucleotides (uppercase ASCII)
    _A_BYTE_V3 = 65  # 'A'
    _C_BYTE_V3 = 67  # 'C'
    _G_BYTE_V3 = 71  # 'G'
    _T_BYTE_V3 = 84  # 'T'
    _U_BYTE_V3 = 85  # 'U'

    # Energetic constants — must mirror viennarna_fallback.py exactly
    _STACKING_BONUS_V3 = -0.5
    _HAIRPIN_INIT_V3 = 3.4
    _LOOP_CLOSURE_COEFF_V3 = 1.75
    _RT_37C_V3 = 0.616
    _MIN_LOOP_PENALTY_V3 = 5.4
    _DEFAULT_MIN_LOOP_V3 = 3

    @njit(cache=True, nogil=True)
    def _numba_count_dinucleotides(seq_array, pattern):
        """Count occurrences of a 2-character pattern in a DNA byte sequence.

        NUMBA-JIT kernel that scans ``seq_array`` for the 2-byte ``pattern``.
        Used as the building block for CpG dinucleotide counting
        (``b'CG'``) — replacing ``str.count`` calls in the CAI recovery
        passes.  While ``str.count`` is itself a C-level fast path, this
        kernel enables zero-copy reuse of the same byte array across
        multiple dinucleotide checks (CpG, ATTTA-precheck, etc.) and
        avoids re-encoding the sequence for every check.

        Args:
            seq_array: numpy uint8 array of DNA characters (ASCII bytes).
            pattern:   numpy uint8 array of exactly 2 bytes (e.g. b'CG').

        Returns:
            Number of (possibly overlapping) occurrences of the pattern.
            Returns 0 for sequences shorter than 2 bytes.
        """
        n = len(seq_array)
        if n < 2 or len(pattern) < 2:
            return 0
        b0 = pattern[0]
        b1 = pattern[1]
        count = 0
        for i in range(n - 1):
            if seq_array[i] == b0 and seq_array[i + 1] == b1:
                count += 1
        return count

    @njit(cache=True, nogil=True)
    def _numba_has_premature_stop(seq_array, n_codons):
        """Check if a DNA sequence contains a premature in-frame stop codon.

        Scans codons 0 .. n_codons-2 (i.e., excludes the final codon)
        for the stop codons TAA, TAG, TGA.  Returns ``True`` on the
        first match (early-exit).

        This mirrors the pure-Python ``_has_premature_stop`` in
        ``pipeline_cai.py`` which uses ``range(0, len(seq) - 5, 3)`` —
        that range excludes the last codon when ``len(seq) = 3 * n_codons``.

        Args:
            seq_array: numpy uint8 array of DNA characters.
            n_codons:  Total number of in-frame codons (including the
                       final codon).  Codons 0 .. n_codons-2 are checked.

        Returns:
            True if any of codons 0 .. n_codons-2 is a stop codon
            (TAA, TAG, or TGA); False otherwise (including empty
            sequences and sequences with fewer than 2 codons).
        """
        n = len(seq_array)
        if n_codons < 2:
            return False

        T = _T_BYTE_V3
        A = _A_BYTE_V3
        G = _G_BYTE_V3

        # Check codons 0 .. n_codons-2 (exclude last codon).
        last_ci = n_codons - 1
        for ci in range(last_ci):
            i = ci * 3
            # Bounds check: need 3 bytes available.
            if i + 3 > n:
                break
            if seq_array[i] == T:
                b1 = seq_array[i + 1]
                b2 = seq_array[i + 2]
                # TAA
                if b1 == A and b2 == A:
                    return True
                # TAG
                if b1 == A and b2 == G:
                    return True
                # TGA
                if b1 == G and b2 == A:
                    return True
        return False

    @njit(cache=True, nogil=True)
    def _numba_compute_approx_dg(seq_array):
        """Nussinov-like mRNA secondary-structure ΔG approximation.

        Ports the pure-Python ``compute_approx_dg`` from
        ``biocompiler.engines.viennarna_fallback`` to a single @njit
        kernel.  This eliminates 85K+ calls to ``_can_pair``,
        32K+ calls to ``_pair_energy`` and 7.8K+ calls to
        ``_loop_penalty`` per optimization — currently the #1
        bottleneck (0.061s cumulative, ~20% of total runtime).

        Algorithm (identical to pure Python):

        1. Two O(n²) scans for hairpin stems:
           - Pass 1: iterate over every (stem_start, loop_size) and
             extend the stem outward while base-pairing holds.
           - Pass 2: iterate over every (loop_center, loop_half) and
             extend the stem outward from the symmetric loop.

        2. For each candidate stem, sum:
           - Base-pair energies (GC: -2.4, AU: -1.5, GU: -0.8 kcal/mol;
             T is treated as U).
           - Stacking bonus (-0.5 kcal/mol) for each consecutive pair
             after the first.
           - Hairpin loop penalty (Turner-inspired: max(5.4,
             3.4 + 1.75 * 0.616 * ln(loop_size))).

        3. Return the most negative ΔG found, rounded to 2 decimal places.

        Args:
            seq_array: numpy uint8 array of DNA characters (T or U both
                       accepted; uppercase ASCII).

        Returns:
            Estimated ΔG in kcal/mol (negative = stable structure).
            Returns 0.0 for sequences too short to form a hairpin
            (< 7 nt) or for empty input.
        """
        n = len(seq_array)

        # Min hairpin: 2 * DEFAULT_MIN_LOOP + 1 = 7
        if n < 2 * _DEFAULT_MIN_LOOP_V3 + 1:
            return 0.0

        best_dg = 0.0

        # ── Pass 1: scan by (stem_start, loop_size) ────────────────
        for stem_start in range(n):
            max_loop = 20
            if n - stem_start < max_loop:
                max_loop = n - stem_start
            for loop_size in range(_DEFAULT_MIN_LOOP_V3, max_loop):
                loop_end = stem_start + loop_size
                if loop_end >= n:
                    break

                stem_pairs_count = 0
                stem_dg = 0.0

                # k iterates outward from (stem_start-1, loop_end)
                upper_k = stem_start
                if n - loop_end < upper_k:
                    upper_k = n - loop_end
                for k in range(upper_k):
                    i = stem_start - 1 - k
                    j = loop_end + k
                    if i < 0 or j >= n:
                        break

                    # Read bytes, treating T (84) as U (85)
                    a = seq_array[i]
                    if a == _T_BYTE_V3:
                        a = _U_BYTE_V3
                    b = seq_array[j]
                    if b == _T_BYTE_V3:
                        b = _U_BYTE_V3

                    # Inline _can_pair + _pair_energy:
                    # GC/CG -> -2.4, AU/UA -> -1.5, GU/UG -> -0.8
                    pair_dg = 0.0
                    can_pair = False
                    if a == _G_BYTE_V3 and b == _C_BYTE_V3:
                        pair_dg = -2.4
                        can_pair = True
                    elif a == _C_BYTE_V3 and b == _G_BYTE_V3:
                        pair_dg = -2.4
                        can_pair = True
                    elif a == _A_BYTE_V3 and b == _U_BYTE_V3:
                        pair_dg = -1.5
                        can_pair = True
                    elif a == _U_BYTE_V3 and b == _A_BYTE_V3:
                        pair_dg = -1.5
                        can_pair = True
                    elif a == _G_BYTE_V3 and b == _U_BYTE_V3:
                        pair_dg = -0.8
                        can_pair = True
                    elif a == _U_BYTE_V3 and b == _G_BYTE_V3:
                        pair_dg = -0.8
                        can_pair = True

                    if not can_pair:
                        break

                    stem_dg += pair_dg
                    stem_pairs_count += 1
                    if k > 0:
                        stem_dg += _STACKING_BONUS_V3

                if stem_pairs_count < 2:
                    continue

                # Inline _loop_penalty(loop_size):
                # loop_size >= 3 here (range starts at DEFAULT_MIN_LOOP=3).
                penalty = (_HAIRPIN_INIT_V3
                           + _LOOP_CLOSURE_COEFF_V3 * _RT_37C_V3
                           * math.log(float(loop_size)))
                if penalty < _MIN_LOOP_PENALTY_V3:
                    penalty = _MIN_LOOP_PENALTY_V3

                total_dg = stem_dg + penalty
                if total_dg < best_dg:
                    best_dg = total_dg

        # ── Pass 2: scan by (loop_center, loop_half) ───────────────
        for center in range(_DEFAULT_MIN_LOOP_V3 + 1, n - _DEFAULT_MIN_LOOP_V3 - 1):
            max_half = 15
            if center < max_half:
                max_half = center
            if n - center < max_half:
                max_half = n - center
            for loop_half in range(_DEFAULT_MIN_LOOP_V3, max_half):
                loop_start = center - loop_half
                loop_end = center + loop_half
                if loop_start < 0 or loop_end >= n:
                    break
                actual_loop_size = loop_end - loop_start + 1

                stem_dg = 0.0
                stem_pairs_count = 0

                upper_k = loop_start
                if n - loop_end - 1 < upper_k:
                    upper_k = n - loop_end - 1
                for k in range(upper_k):
                    i = loop_start - 1 - k
                    j = loop_end + 1 + k
                    if i < 0 or j >= n:
                        break

                    a = seq_array[i]
                    if a == _T_BYTE_V3:
                        a = _U_BYTE_V3
                    b = seq_array[j]
                    if b == _T_BYTE_V3:
                        b = _U_BYTE_V3

                    pair_dg = 0.0
                    can_pair = False
                    if a == _G_BYTE_V3 and b == _C_BYTE_V3:
                        pair_dg = -2.4
                        can_pair = True
                    elif a == _C_BYTE_V3 and b == _G_BYTE_V3:
                        pair_dg = -2.4
                        can_pair = True
                    elif a == _A_BYTE_V3 and b == _U_BYTE_V3:
                        pair_dg = -1.5
                        can_pair = True
                    elif a == _U_BYTE_V3 and b == _A_BYTE_V3:
                        pair_dg = -1.5
                        can_pair = True
                    elif a == _G_BYTE_V3 and b == _U_BYTE_V3:
                        pair_dg = -0.8
                        can_pair = True
                    elif a == _U_BYTE_V3 and b == _G_BYTE_V3:
                        pair_dg = -0.8
                        can_pair = True

                    if not can_pair:
                        break

                    stem_dg += pair_dg
                    stem_pairs_count += 1
                    if k > 0:
                        stem_dg += _STACKING_BONUS_V3

                if stem_pairs_count < 2:
                    continue

                # actual_loop_size = 2*loop_half+1, always >= 7 here.
                penalty = (_HAIRPIN_INIT_V3
                           + _LOOP_CLOSURE_COEFF_V3 * _RT_37C_V3
                           * math.log(float(actual_loop_size)))
                if penalty < _MIN_LOOP_PENALTY_V3:
                    penalty = _MIN_LOOP_PENALTY_V3

                total_dg = stem_dg + penalty
                if total_dg < best_dg:
                    best_dg = total_dg

        # Round to 2 decimal places (matches pure-Python `round(x, 2)`).
        return round(best_dg, 2)

    # ── Pre-warm: trigger NUMBA compilation at import time ──────────
    # This adds ~1-2s to first import but avoids JIT latency during
    # optimization. We use tiny dummy arrays to minimize compile time.
    def _warmup() -> None:
        """Pre-compile all NUMBA kernels with tiny dummy inputs."""
        try:
            _dummy_seq = np.array([65, 84, 71, 67], dtype=np.uint8)  # b'ATGC'
            _dummy_dinuc = np.array([71, 84], dtype=np.uint8)  # b'GT'
            _dummy_adapt = np.array([0.5, 1.0], dtype=np.float64)
            _dummy_indices = np.array([0, 1], dtype=np.int64)

            # Existing kernels
            count_gc(_dummy_seq)
            count_dinucleotides(_dummy_seq, _dummy_dinuc)
            compute_cai_kernel(_dummy_adapt, _dummy_indices, 2)
            scan_restriction_sites(_dummy_seq, _dummy_dinuc, 2)
            find_all_dinucleotide_positions(_dummy_seq, _dummy_dinuc)

            # New v2 kernels
            compute_cai_incremental(-1.3862943611198906, 2, 0.5, 1.0)
            batch_codon_swap_score(
                _dummy_adapt, _dummy_indices, 2, 0,
                np.array([0, 1], dtype=np.int64), 2, -1.3862943611198906
            )
            fast_gc_window(_dummy_seq, 2)
            fast_dinucleotide_count(
                _dummy_seq,
                np.array([[71, 84], [67, 71]], dtype=np.uint8), 2
            )

            # Additional utility kernels
            count_gc_parallel(_dummy_seq)
            scan_restriction_sites_multi(
                _dummy_seq, _dummy_dinuc,
                np.array([0], dtype=np.int64),
                np.array([2], dtype=np.int64), 1
            )

            # v3 kernels
            _numba_count_dinucleotides(_dummy_seq, _dummy_dinuc)
            _numba_has_premature_stop(_dummy_seq, 1)
            # Use a slightly longer dummy for the ΔG kernel so both
            # passes execute (need >= 7 nt for any structure).
            _dummy_seq_long = np.array(
                [71, 71, 71, 65, 65, 65, 67, 67, 67], dtype=np.uint8
            )  # b'GGGAAACCC' — canonical hairpin
            _numba_compute_approx_dg(_dummy_seq_long)
        except Exception:
            pass  # Warmup failure is non-fatal

    _warmup()

else:
    # ════════════════════════════════════════════════════════════════
    # Pure-Python fallbacks (when NUMBA is not available)
    # ════════════════════════════════════════════════════════════════
    import array

    # numpy may still be available even if numba is not
    try:
        import numpy as _np
    except ImportError:
        _np = None  # type: ignore[assignment]

    def _as_array(data, dtype=None):
        """Return a numpy array if numpy is available, otherwise the raw list."""
        if _np is not None:
            return _np.array(data, dtype=dtype)
        return data

    def count_gc(seq_bytes) -> int:  # type: ignore[misc]
        """Pure-Python fallback: Count G and C characters."""
        count = 0
        for b in seq_bytes:
            if b == _G_BYTE or b == _C_BYTE:
                count += 1
        return count

    def count_dinucleotides(seq_bytes, dinuc_bytes) -> int:  # type: ignore[misc]
        """Pure-Python fallback: Count dinucleotide occurrences."""
        count = 0
        b0 = dinuc_bytes[0]
        b1 = dinuc_bytes[1]
        n = len(seq_bytes)
        for i in range(n - 1):
            if seq_bytes[i] == b0 and seq_bytes[i + 1] == b1:
                count += 1
        return count

    def compute_cai_kernel(adaptiveness_array, codon_indices, n_codons) -> float:  # type: ignore[misc]
        """Pure-Python fallback: Compute CAI via geometric mean of log ratios."""
        if n_codons == 0:
            return 0.0
        log_sum = 0.0
        epsilon = 1e-10
        for i in range(n_codons):
            idx = codon_indices[i]
            w = adaptiveness_array[idx]
            if w <= 0.0:
                w = epsilon
            log_sum += math.log(w)
        return math.exp(log_sum / n_codons)

    def scan_restriction_sites(seq_bytes, pattern_bytes, pattern_len):  # type: ignore[misc]
        """Pure-Python fallback: Find all positions of a pattern."""
        n = len(seq_bytes)
        results = []
        if pattern_len > n:
            return _as_array(results, dtype=_np.int64 if _np is not None else None)
        for i in range(n - pattern_len + 1):
            match = True
            for j in range(pattern_len):
                if seq_bytes[i + j] != pattern_bytes[j]:
                    match = False
                    break
            if match:
                results.append(i)
        return _as_array(results, dtype=_np.int64 if _np is not None else None)

    def find_all_dinucleotide_positions(seq_bytes, dinuc_bytes):  # type: ignore[misc]
        """Pure-Python fallback: Find all dinucleotide positions."""
        n = len(seq_bytes)
        results = []
        b0 = dinuc_bytes[0]
        b1 = dinuc_bytes[1]
        for i in range(n - 1):
            if seq_bytes[i] == b0 and seq_bytes[i + 1] == b1:
                results.append(i)
        return _as_array(results, dtype=_np.int64 if _np is not None else None)

    def seq_to_bytes(seq: str):  # type: ignore[misc]
        """Pure-Python fallback: Convert DNA string to byte array."""
        if _np is not None:
            return _np.frombuffer(seq.encode('ascii'), dtype=_np.uint8)
        return array.array('B', seq.encode('ascii'))

    # ── New v2 kernel fallbacks ─────────────────────────────────────

    def compute_cai_incremental(
        current_log_sum: float,
        n_codons: int,
        old_adaptiveness: float,
        new_adaptiveness: float,
    ) -> float:
        """Pure-Python fallback: Update CAI after a single codon swap.

        Instead of recomputing the entire log-sum, adjust by subtracting
        the old codon's log-adaptiveness and adding the new one.
        """
        # Edge cases: empty sequence or invalid codon count
        if n_codons <= 0:
            return 0.0
        epsilon = 1e-10
        w_old = old_adaptiveness
        if w_old <= 0.0 or w_old != w_old:  # catch NaN
            w_old = epsilon
        w_new = new_adaptiveness
        if w_new <= 0.0 or w_new != w_new:  # catch NaN
            w_new = epsilon
        new_log_sum = current_log_sum - math.log(w_old) + math.log(w_new)
        return math.exp(new_log_sum / n_codons)

    def batch_codon_swap_score(
        adaptiveness_array,
        codon_indices,
        n_codons: int,
        swap_position: int,
        candidate_indices,
        n_candidates: int,
        current_log_sum: float,
    ):
        """Pure-Python fallback: Score all candidate codon swaps at a position."""
        # No alternative codons (e.g. Met/Trp): return empty array
        if n_candidates == 0:
            return _as_array([], dtype=_np.float64 if _np is not None else None)
        scores = []
        if n_codons <= 0:
            return _as_array([0.0] * n_candidates, dtype=_np.float64 if _np is not None else None)
        epsilon = 1e-10
        old_idx = codon_indices[swap_position]
        w_old = adaptiveness_array[old_idx]
        if w_old <= 0.0 or w_old != w_old:  # catch NaN
            w_old = epsilon
        log_w_old = math.log(w_old)
        for k in range(n_candidates):
            w_new = adaptiveness_array[candidate_indices[k]]
            if w_new <= 0.0 or w_new != w_new:  # catch NaN
                w_new = epsilon
            new_log_sum = current_log_sum - log_w_old + math.log(w_new)
            scores.append(math.exp(new_log_sum / n_codons))
        return _as_array(scores, dtype=_np.float64 if _np is not None else None)

    def fast_gc_window(seq_bytes, window_size: int):
        """Pure-Python fallback: Sliding-window GC% with incremental updates."""
        n = len(seq_bytes)
        n_windows = n - window_size + 1
        if n_windows <= 0:
            return _as_array([], dtype=_np.float64 if _np is not None else None)
        results = []
        gc_count = 0
        for i in range(window_size):
            b = seq_bytes[i]
            if b == _G_BYTE or b == _C_BYTE:
                gc_count += 1
        results.append(gc_count / window_size)
        for i in range(1, n_windows):
            outgoing = seq_bytes[i - 1]
            incoming = seq_bytes[i + window_size - 1]
            if outgoing == _G_BYTE or outgoing == _C_BYTE:
                gc_count -= 1
            if incoming == _G_BYTE or incoming == _C_BYTE:
                gc_count += 1
            results.append(gc_count / window_size)
        return _as_array(results, dtype=_np.float64 if _np is not None else None)

    def fast_dinucleotide_count(seq_bytes, dinuc_keys, n_dinucs: int):
        """Pure-Python fallback: Count multiple dinucleotides in one pass."""
        counts = [0] * n_dinucs
        n = len(seq_bytes)
        # Sequence shorter than 2 bases cannot contain any dinucleotide
        if n < 2:
            return _as_array(counts, dtype=_np.int64 if _np is not None else None)
        for i in range(n - 1):
            b0 = seq_bytes[i]
            b1 = seq_bytes[i + 1]
            for d in range(n_dinucs):
                if b0 == dinuc_keys[d][0] and b1 == dinuc_keys[d][1]:
                    counts[d] += 1
        return _as_array(counts, dtype=_np.int64 if _np is not None else None)

    # ── Additional utility fallbacks ────────────────────────────────

    def count_gc_parallel(seq_bytes) -> int:
        """Pure-Python fallback: Same as count_gc (no parallelism available)."""
        return count_gc(seq_bytes)

    def scan_restriction_sites_multi(
        seq_bytes, pattern_bytes, pattern_offsets, pattern_lens, n_patterns: int
    ):
        """Pure-Python fallback: Find positions of multiple patterns."""
        results = []
        for p in range(n_patterns):
            plen = pattern_lens[p]
            poff = pattern_offsets[p]
            for i in range(len(seq_bytes) - plen + 1):
                match = True
                for j in range(plen):
                    if seq_bytes[i + j] != pattern_bytes[poff + j]:
                        match = False
                        break
                if match:
                    results.append((i << 16) | p)
        return _as_array(results, dtype=_np.int64 if _np is not None else None)

    # ── v3 kernel fallbacks ──────────────────────────────────────────

    def _numba_count_dinucleotides(seq_array, pattern) -> int:
        """Pure-Python fallback: Count occurrences of a 2-character pattern.

        Mirrors the NUMBA ``_numba_count_dinucleotides`` kernel for
        environments where NUMBA is unavailable.  Used for CpG
        dinucleotide counting in the CAI recovery pass.
        """
        n = len(seq_array)
        if n < 2 or len(pattern) < 2:
            return 0
        b0 = pattern[0]
        b1 = pattern[1]
        count = 0
        for i in range(n - 1):
            if seq_array[i] == b0 and seq_array[i + 1] == b1:
                count += 1
        return count

    def _numba_has_premature_stop(seq_array, n_codons: int) -> bool:
        """Pure-Python fallback: Detect in-frame premature stop codons.

        Mirrors the NUMBA ``_numba_has_premature_stop`` kernel.  Scans
        codons 0 .. n_codons-2 for TAA, TAG, TGA.
        """
        n = len(seq_array)
        if n_codons < 2:
            return False
        T = 84  # ord('T')
        A = 65  # ord('A')
        G = 71  # ord('G')
        last_ci = n_codons - 1
        for ci in range(last_ci):
            i = ci * 3
            if i + 3 > n:
                break
            if seq_array[i] == T:
                b1 = seq_array[i + 1]
                b2 = seq_array[i + 2]
                if b1 == A and (b2 == A or b2 == G):
                    return True
                if b1 == G and b2 == A:
                    return True
        return False

    def _numba_compute_approx_dg(seq_array) -> float:
        """Pure-Python fallback: Nussinov-like mRNA ΔG approximation.

        Mirrors the NUMBA ``_numba_compute_approx_dg`` kernel.  This is
        the byte-array equivalent of
        ``biocompiler.engines.viennarna_fallback.compute_approx_dg``,
        kept here so callers can use a single API regardless of NUMBA
        availability.
        """
        n = len(seq_array)
        DEFAULT_MIN_LOOP = 3
        STACKING_BONUS = -0.5
        HAIRPIN_INIT = 3.4
        LOOP_CLOSURE_COEFF = 1.75
        RT_37C = 0.616
        MIN_LOOP_PENALTY = 5.4
        A = 65; C = 67; G = 71; T = 84; U = 85

        if n < 2 * DEFAULT_MIN_LOOP + 1:
            return 0.0

        def _pair_dg(a, b):
            # Treat T as U
            if a == T: a = U
            if b == T: b = U
            if a == G and b == C: return -2.4, True
            if a == C and b == G: return -2.4, True
            if a == A and b == U: return -1.5, True
            if a == U and b == A: return -1.5, True
            if a == G and b == U: return -0.8, True
            if a == U and b == G: return -0.8, True
            return 0.0, False

        def _loop_pen(loop_size):
            if loop_size <= 0:
                return 0.0
            if loop_size >= 3:
                pen = HAIRPIN_INIT + LOOP_CLOSURE_COEFF * RT_37C * math.log(loop_size)
                return max(pen, MIN_LOOP_PENALTY)
            return MIN_LOOP_PENALTY + 2.0 * (3 - loop_size)

        best_dg = 0.0

        # Pass 1
        for stem_start in range(n):
            max_loop = min(20, n - stem_start)
            for loop_size in range(DEFAULT_MIN_LOOP, max_loop):
                loop_end = stem_start + loop_size
                if loop_end >= n:
                    break
                stem_pairs_count = 0
                stem_dg = 0.0
                upper_k = min(stem_start, n - loop_end)
                for k in range(upper_k):
                    i = stem_start - 1 - k
                    j = loop_end + k
                    if i < 0 or j >= n:
                        break
                    pdg, can = _pair_dg(seq_array[i], seq_array[j])
                    if not can:
                        break
                    stem_dg += pdg
                    stem_pairs_count += 1
                    if k > 0:
                        stem_dg += STACKING_BONUS
                if stem_pairs_count < 2:
                    continue
                total = stem_dg + _loop_pen(loop_size)
                if total < best_dg:
                    best_dg = total

        # Pass 2
        for center in range(DEFAULT_MIN_LOOP + 1, n - DEFAULT_MIN_LOOP - 1):
            max_half = min(15, center, n - center)
            for loop_half in range(DEFAULT_MIN_LOOP, max_half):
                loop_start = center - loop_half
                loop_end = center + loop_half
                if loop_start < 0 or loop_end >= n:
                    break
                actual = loop_end - loop_start + 1
                stem_dg = 0.0
                stem_pairs_count = 0
                upper_k = min(loop_start, n - loop_end - 1)
                for k in range(upper_k):
                    i = loop_start - 1 - k
                    j = loop_end + 1 + k
                    if i < 0 or j >= n:
                        break
                    pdg, can = _pair_dg(seq_array[i], seq_array[j])
                    if not can:
                        break
                    stem_dg += pdg
                    stem_pairs_count += 1
                    if k > 0:
                        stem_dg += STACKING_BONUS
                if stem_pairs_count < 2:
                    continue
                total = stem_dg + _loop_pen(actual)
                if total < best_dg:
                    best_dg = total

        return round(best_dg, 2)


# ─── Biosecurity fuzzy matching kernels ─────────────────────────────

if HAS_NUMBA:
    @numba.njit(cache=True)
    def _numba_has_shared_kmer(window_bytes, motif_kmer_hashes, k=4):
        """Check if any k-mer in window appears in motif_kmer_hashes.
        
        Args:
            window_bytes: bytes object (the window as ASCII)
            motif_kmer_hashes: sorted numpy array (uint64) of k-mer hashes
            k: k-mer length (default 4)
        
        Returns True if any k-mer is shared.
        """
        n = len(window_bytes)
        if n < k or len(motif_kmer_hashes) == 0:
            return True
        for i in range(n - k + 1):
            # Simple hash: sum of bytes * position weight
            h = 0
            for j in range(k):
                h = h * 2654435761 + window_bytes[i + j]
            # Binary search in sorted array
            lo, hi = 0, len(motif_kmer_hashes)
            while lo < hi:
                mid = (lo + hi) // 2
                if motif_kmer_hashes[mid] < h:
                    lo = mid + 1
                elif motif_kmer_hashes[mid] > h:
                    hi = mid
                else:
                    return True
        return False
    
    @numba.njit(cache=True)
    def _numba_hamming_distance(s1_bytes, s2_bytes):
        """Compute Hamming distance between two byte strings."""
        n = len(s1_bytes)
        dist = 0
        for i in range(n):
            if s1_bytes[i] != s2_bytes[i]:
                dist += 1
        return dist
    
    @numba.njit(cache=True)
    def _numba_fuzzy_match_hamming(seq_bytes, motif_bytes, max_distance=2, k=4):
        """Find fuzzy matches of motif in sequence using Hamming distance.
        
        Returns array of (position, distance) pairs.
        """
        slen = len(seq_bytes)
        mlen = len(motif_bytes)
        if mlen == 0 or slen < mlen:
            return np.zeros((0, 2), dtype=np.int64)
        
        # Build motif k-mer hashes
        n_motif_kmers = mlen - k + 1
        if n_motif_kmers <= 0:
            return np.zeros((0, 2), dtype=np.int64)
        
        motif_hashes = np.empty(n_motif_kmers, dtype=np.int64)
        for i in range(n_motif_kmers):
            h = 0
            for j in range(k):
                h = h * 2654435761 + motif_bytes[i + j]
            motif_hashes[i] = h
        motif_hashes.sort()
        
        # Scan sequence
        max_results = slen - mlen + 1
        positions = np.empty(max_results, dtype=np.int64)
        distances = np.empty(max_results, dtype=np.int64)
        count = 0
        
        for i in range(slen - mlen + 1):
            # Check shared k-mer
            has_shared = False
            for j in range(mlen - k + 1):
                h = 0
                for jj in range(k):
                    h = h * 2654435761 + seq_bytes[i + j + jj]
                # Binary search
                lo, hi = 0, len(motif_hashes)
                while lo < hi:
                    mid = (lo + hi) // 2
                    if motif_hashes[mid] < h:
                        lo = mid + 1
                    elif motif_hashes[mid] > h:
                        hi = mid
                    else:
                        has_shared = True
                        break
                if has_shared:
                    break
            
            if not has_shared:
                continue
            
            # Compute Hamming distance
            dist = 0
            for j in range(mlen):
                if seq_bytes[i + j] != motif_bytes[j]:
                    dist += 1
                    if dist > max_distance:
                        break
            
            if 1 <= dist <= max_distance:
                positions[count] = i
                distances[count] = dist
                count += 1
        
        return np.column_stack((positions[:count], distances[:count]))
else:
    _numba_has_shared_kmer = None
    _numba_hamming_distance = None
    _numba_fuzzy_match_hamming = None


# ─── Levenshtein distance kernel ────────────────────────────────────

if HAS_NUMBA:
    @numba.njit(cache=True)
    def _numba_levenshtein_within(s1_bytes, s2_bytes, max_dist):
        """Check if Levenshtein distance <= max_dist (early exit)."""
        n = len(s1_bytes)
        m = len(s2_bytes)
        if abs(n - m) > max_dist:
            return False
        # Use two rows for DP
        prev = np.empty(m + 1, dtype=np.int64)
        curr = np.empty(m + 1, dtype=np.int64)
        for j in range(m + 1):
            prev[j] = j
        for i in range(1, n + 1):
            curr[0] = i
            min_in_row = i
            for j in range(1, m + 1):
                cost = 0 if s1_bytes[i-1] == s2_bytes[j-1] else 1
                curr[j] = min(prev[j] + 1, curr[j-1] + 1, prev[j-1] + cost)
                if curr[j] < min_in_row:
                    min_in_row = curr[j]
            if min_in_row > max_dist:
                return False
            prev, curr = curr, prev
        return prev[m] <= max_dist
    
    @numba.njit(cache=True)
    def _numba_levenshtein_distance(s1_bytes, s2_bytes):
        """Full Levenshtein distance."""
        n = len(s1_bytes)
        m = len(s2_bytes)
        if n == 0: return m
        if m == 0: return n
        prev = np.empty(m + 1, dtype=np.int64)
        curr = np.empty(m + 1, dtype=np.int64)
        for j in range(m + 1):
            prev[j] = j
        for i in range(1, n + 1):
            curr[0] = i
            for j in range(1, m + 1):
                cost = 0 if s1_bytes[i-1] == s2_bytes[j-1] else 1
                curr[j] = min(prev[j] + 1, curr[j-1] + 1, prev[j-1] + cost)
            prev, curr = curr, prev
        return prev[m]
    
    @numba.njit(cache=True)
    def _numba_has_shared_kmer_fast(window_bytes, motif_kmer_hashes, k=4):
        """Fast k-mer check using pre-computed sorted hash array."""
        n = len(window_bytes)
        if n < k or len(motif_kmer_hashes) == 0:
            return True
        nkmers = n - k + 1
        for i in range(nkmers):
            h = np.int64(0)
            for j in range(k):
                h = h * np.int64(2654435761) + np.int64(window_bytes[i + j])
            lo, hi = 0, len(motif_kmer_hashes)
            while lo < hi:
                mid = (lo + hi) // 2
                if motif_kmer_hashes[mid] < h:
                    lo = mid + 1
                elif motif_kmer_hashes[mid] > h:
                    hi = mid
                else:
                    return True
        return False
else:
    _numba_levenshtein_within = None
    _numba_levenshtein_distance = None
    _numba_has_shared_kmer_fast = None
