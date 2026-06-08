"""Type stubs for biocompiler.numba_kernels — NUMBA-accelerated computation kernels."""

from __future__ import annotations

from typing import Any

import numpy as _np


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

HAS_NUMBA: bool
USE_NUMBA: bool


# ────────────────────────────────────────────────────────────
# Sequence utilities
# ────────────────────────────────────────────────────────────

def seq_to_bytes(seq: str) -> _np.ndarray: ...


# ────────────────────────────────────────────────────────────
# Counting kernels
# ────────────────────────────────────────────────────────────

def count_gc(seq_bytes: _np.ndarray) -> int: ...
def fast_dinucleotide_count(
    seq_bytes: _np.ndarray,
    dinuc_keys: _np.ndarray,
    n_dinucs: int,
) -> _np.ndarray: ...


# ────────────────────────────────────────────────────────────
# CAI kernels
# ────────────────────────────────────────────────────────────

def compute_cai_kernel(
    adapt_array: _np.ndarray,
    indices: _np.ndarray,
    n_codons: int,
) -> float: ...
def compute_cai_incremental(
    log_sum: float,
    n_codons: int,
    old_adapt: float,
    new_adapt: float,
) -> tuple[float, int]: ...


# ────────────────────────────────────────────────────────────
# GC window kernel
# ────────────────────────────────────────────────────────────

def fast_gc_window(
    seq_bytes: _np.ndarray,
    window_size: int,
) -> _np.ndarray: ...
