"""
BioCompiler NUMBA JIT-Compiled Kernels
=======================================

NUMBA-accelerated inner loops for codon optimization hot paths:

1. count_dinucleotides  — Count occurrences of a dinucleotide (GT, CG, AG) in a DNA sequence
2. count_gc             — Count G and C characters in a DNA sequence
3. compute_cai_kernel   — Compute CAI from pre-indexed adaptiveness values (log-sum-exp)
4. scan_restriction_sites — Find all positions of a restriction site pattern
5. find_all_dinucleotide_positions — Find all positions of a dinucleotide

All kernels use @numba.njit(cache=True) for maximum performance with cached compilation.
Pure-Python fallbacks are provided when NUMBA is unavailable.

Architecture:
    - Input: byte arrays (np.frombuffer or np.array) for zero-copy NUMBA access
    - Output: scalar int/float or numpy arrays of positions
    - All kernels are stateless and side-effect free
    - Compilation cache avoids re-JIT overhead on repeated imports

Usage:
    from biocompiler.numba_kernels import HAS_NUMBA, count_gc, count_dinucleotides

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
    from numba import njit

    HAS_NUMBA: bool = True
    _NUMBA_VERSION: str = numba.__version__
except ImportError:
    HAS_NUMBA = False
    _NUMBA_VERSION = ""
    numba = None  # type: ignore[assignment]

__all__ = [
    "HAS_NUMBA",
    "count_gc",
    "count_dinucleotides",
    "compute_cai_kernel",
    "scan_restriction_sites",
    "find_all_dinucleotide_positions",
    "seq_to_bytes",
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

    @njit(cache=True)
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

    @njit(cache=True)
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

    @njit(cache=True)
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

    @njit(cache=True)
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

    @njit(cache=True)
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

            count_gc(_dummy_seq)
            count_dinucleotides(_dummy_seq, _dummy_dinuc)
            compute_cai_kernel(_dummy_adapt, _dummy_indices, 2)
            scan_restriction_sites(_dummy_seq, _dummy_dinuc, 2)
            find_all_dinucleotide_positions(_dummy_seq, _dummy_dinuc)
        except Exception:
            pass  # Warmup failure is non-fatal

    _warmup()

else:
    # ════════════════════════════════════════════════════════════════
    # Pure-Python fallbacks (when NUMBA is not available)
    # ════════════════════════════════════════════════════════════════
    import array

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

    def scan_restriction_sites(seq_bytes, pattern_bytes, pattern_len) -> list:  # type: ignore[misc]
        """Pure-Python fallback: Find all positions of a pattern."""
        n = len(seq_bytes)
        results = []
        if pattern_len > n:
            return results
        for i in range(n - pattern_len + 1):
            match = True
            for j in range(pattern_len):
                if seq_bytes[i + j] != pattern_bytes[j]:
                    match = False
                    break
            if match:
                results.append(i)
        return results

    def find_all_dinucleotide_positions(seq_bytes, dinuc_bytes) -> list:  # type: ignore[misc]
        """Pure-Python fallback: Find all dinucleotide positions."""
        n = len(seq_bytes)
        results = []
        b0 = dinuc_bytes[0]
        b1 = dinuc_bytes[1]
        for i in range(n - 1):
            if seq_bytes[i] == b0 and seq_bytes[i + 1] == b1:
                results.append(i)
        return results

    def seq_to_bytes(seq: str):  # type: ignore[misc]
        """Pure-Python fallback: Convert DNA string to byte array."""
        return array.array('B', seq.encode('ascii'))
