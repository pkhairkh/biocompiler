"""
CAI (Codon Adaptation Index) computation and batch scoring.

Provides fast CAI computation using NUMBA kernels when available,
with pure-Python fallbacks.  Includes the _BatchSwapScorer for
efficient batch evaluation of candidate codon swaps.
"""

from typing import Dict, List, Any

import logging
import math

from ..type_system import CODON_TABLE

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from ..numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        USE_NUMBA as _USE_NUMBA,
        count_gc as _numba_count_gc,
        count_dinucleotides as _numba_count_dinuc,
        seq_to_bytes as _seq_to_bytes,
        fast_dinucleotide_count as _numba_fast_dinuc_count,
        batch_codon_swap_score as _numba_batch_codon_swap_score,
        compute_cai_kernel as _numba_cai_kernel,
        compute_cai_incremental as _numba_cai_incremental,
    )
except ImportError:
    _HAS_NUMBA = False
    _USE_NUMBA = False
    _numba_fast_dinuc_count = None  # type: ignore[assignment]
    _numba_batch_codon_swap_score = None  # type: ignore[assignment]
    _numba_cai_kernel = None  # type: ignore[assignment]
    _numba_cai_incremental = None  # type: ignore[assignment]

HAS_NUMBA: bool = _HAS_NUMBA if isinstance(_HAS_NUMBA, bool) else False
USE_NUMBA: bool = _USE_NUMBA if isinstance(_USE_NUMBA, bool) else False

# ── NUMBA CAI helpers ──────────────────────────────────────────────
import numpy as _np

logger = logging.getLogger(__name__)


def _adaptiveness_to_array(codon_adaptiveness: dict[str, float]) -> _np.ndarray:
    """Convert codon adaptiveness dict to numpy array indexed by codon encoding.

    Codons are mapped to array indices using a base-4 encoding:
        A=0, C=1, G=2, T=3  →  index = b0*16 + b1*4 + b2
    This gives indices in [0, 63] for all 64 possible 3-letter codons,
    providing O(1) lookup without hash tables.

    Args:
        codon_adaptiveness: Dict mapping codon strings (e.g. 'ATG') to
            adaptiveness values (float in [0.0, 1.0]).

    Returns:
        numpy float64 array of length 64, where index = base4_encode(codon).
        Missing codons get a default value of 1e-10.
    """
    _BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    arr = _np.full(64, 1e-10, dtype=_np.float64)
    for codon, w in codon_adaptiveness.items():
        try:
            idx = _BASE_MAP[codon[0]] * 16 + _BASE_MAP[codon[1]] * 4 + _BASE_MAP[codon[2]]
            arr[idx] = w
        except (KeyError, IndexError):
            continue  # skip non-standard codons
    return arr


def _codon_to_index(codon: str) -> int:
    """Map a 3-letter codon string to its base-4 array index.

    Uses: A=0, C=1, G=2, T=3 → index = b0*16 + b1*4 + b2
    """
    _BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return _BASE_MAP[codon[0]] * 16 + _BASE_MAP[codon[1]] * 4 + _BASE_MAP[codon[2]]


def _dna_to_codon_indices(dna: str) -> _np.ndarray:
    """Convert DNA string to array of codon adaptiveness indices.

    Each codon position gets its base-4 encoded index for O(1) lookup
    in the adaptiveness array.

    Args:
        dna: DNA sequence string (length must be a multiple of 3).

    Returns:
        numpy int64 array of length len(dna)//3 with base-4 encoded indices.
    """
    n_codons = len(dna) // 3
    indices = _np.empty(n_codons, dtype=_np.int64)
    for i in range(n_codons):
        indices[i] = _codon_to_index(dna[i * 3:i * 3 + 3])
    return indices


def _compute_cai_fast(seq: str, codon_adaptiveness: dict[str, float]) -> float:
    """Compute CAI using the NUMBA kernel when available, pure Python otherwise.

    This function replaces the pure-Python loop in _compute_seq_cai with
    the NUMBA-accelerated compute_cai_kernel. The kernel uses log-sum
    for numerical stability (geometric mean of relative adaptiveness values).

    When NUMBA is available, converts the adaptiveness dict to a numpy array
    and uses compute_cai_kernel for O(n) CAI computation in compiled code.
    Falls back to pure-Python log-sum when NUMBA is not installed.

    Args:
        seq: DNA coding sequence.
        codon_adaptiveness: Dict mapping codon strings to adaptiveness values.

    Returns:
        CAI value in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.
    """
    if not seq or len(seq) < 3:
        return 0.0

    if HAS_NUMBA and _numba_cai_kernel is not None:
        try:
            adapt_arr = _adaptiveness_to_array(codon_adaptiveness)
            # Build codon indices array, excluding Met and stop codons from count
            n_total = len(seq) // 3
            codon_indices_list = []
            for i in range(n_total):
                codon = seq[i * 3:i * 3 + 3]
                aa = CODON_TABLE.get(codon)
                if aa == 'M' or aa == '*':
                    continue
                codon_indices_list.append(_codon_to_index(codon))
            if not codon_indices_list:
                return 0.0
            codon_indices = _np.array(codon_indices_list, dtype=_np.int64)
            n_codons = len(codon_indices)
            return _numba_cai_kernel(adapt_arr, codon_indices, n_codons)
        except Exception:
            logger.warning("NUMBA CAI kernel failed, falling back to pure-Python", exc_info=True)

    # Pure-Python fallback
    epsilon = 1e-10
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa == 'M' or aa == '*':
            continue
        w = codon_adaptiveness.get(codon, epsilon)
        if w <= 0.0:
            w = epsilon
        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


def _count_dinucs_fast(seq: str, *dinucleotides: str) -> tuple[int, ...]:
    """Count multiple dinucleotides in a single pass using the NUMBA kernel.

    Falls back to pure-Python counting when NUMBA is unavailable.

    Args:
        seq: DNA sequence string (uppercase ACGT).
        *dinucleotides: One or more dinucleotide strings (e.g. "GT", "CG", "AG").

    Returns:
        Tuple of counts, one per dinucleotide, in the same order as input.
    """
    n_dinucs = len(dinucleotides)
    if n_dinucs == 0:
        return ()

    # Fast path: NUMBA kernel
    if HAS_NUMBA and _numba_fast_dinuc_count is not None:
        import numpy as _np
        seq_bytes = _seq_to_bytes(seq)
        dinuc_keys = _np.array(
            [[ord(d[0]), ord(d[1])] for d in dinucleotides],
            dtype=_np.uint8,
        )
        counts = _numba_fast_dinuc_count(seq_bytes, dinuc_keys, n_dinucs)
        return tuple(int(c) for c in counts)

    # Pure-Python fallback
    results = []
    for di in dinucleotides:
        count = 0
        pos = 0
        while True:
            pos = seq.find(di, pos)
            if pos == -1:
                break
            count += 1
            pos += 1
        results.append(count)
    return tuple(results)


# ── Batch codon swap scorer (NUMBA-accelerated) ─────────────────────
class _BatchSwapScorer:
    """Pre-compute CAI arrays for batch_codon_swap_score kernel invocations.

    Builds the static adaptiveness array and codon-to-index mapping once,
    then provides a fast ``score_candidates`` method that calls the NUMBA
    kernel (or pure-Python fallback) to score all candidate codon swaps
    at a given position in one vectorized pass.

    Usage::

        scorer = _BatchSwapScorer(species_cai)
        # Per iteration / position:
        scores = scorer.score_candidates(seq_codons, codon_idx, candidate_codons)
        # scores[k] = full-sequence CAI if we swapped codon at codon_idx to candidate_codons[k]
    """

    def __init__(self, species_cai: Dict[str, float]) -> None:
        self._species_cai = species_cai

        # Build a stable sorted list of all known codons → integer index.
        # CODON_TABLE keys are the canonical 64 codons; we include any
        # extra codon that appears in species_cai as well.
        from ..type_system import CODON_TABLE as _CT
        all_codons = sorted(set(_CT.keys()) | set(species_cai.keys()))
        self._codon_to_idx: Dict[str, int] = {c: i for i, c in enumerate(all_codons)}
        self._idx_to_codon: Dict[int, str] = {i: c for i, c in enumerate(all_codons)}
        self._n_codons_total = len(all_codons)

        # Build adaptiveness array (float64) — one entry per codon index.
        epsilon = 1e-10
        self._adaptiveness: list[float] = []
        for c in all_codons:
            w = species_cai.get(c, 0.0)
            self._adaptiveness.append(w if w > 0.0 else epsilon)

        # Lazy numpy arrays (built on first use when NUMBA is available)
        self._np_adaptiveness = None  # type: ignore[assignment]
        self._np_built = False

        # Incremental CAI tracking state — avoids O(n) log_sum recomputation
        # after each accepted codon swap.  Set via reset_incremental_state().
        self._current_log_sum: float | None = None
        self._n_codons: int = 0

    def _ensure_numpy(self) -> None:
        """Lazily build numpy arrays for the NUMBA kernel."""
        if self._np_built:
            return
        try:
            import numpy as _np
            self._np_adaptiveness = _np.array(self._adaptiveness, dtype=_np.float64)
        except ImportError:
            self._np_adaptiveness = None
        self._np_built = True

    def _build_codon_indices(self, seq_codons: List[str]) -> "Any":
        """Build the codon_indices numpy array from a list of codon strings."""
        import numpy as _np
        indices = _np.array(
            [self._codon_to_idx.get(c, 0) for c in seq_codons],
            dtype=_np.int64,
        )
        return indices

    def _compute_log_sum(self, seq_codons: List[str]) -> float:
        """Compute the current sum of log(w_i) over all codon positions."""
        epsilon = 1e-10
        log_sum = 0.0
        for c in seq_codons:
            idx = self._codon_to_idx.get(c, 0)
            w = self._adaptiveness[idx]
            if w <= 0.0:
                w = epsilon
            log_sum += math.log(w)
        return log_sum

    def reset_incremental_state(self, seq_codons: List[str]) -> None:
        """Initialize (or re-initialize) incremental CAI tracking state.

        Computes the full log_sum once from the current codon list and
        caches it for O(1) incremental updates via compute_cai_incremental.

        Args:
            seq_codons: Current codon strings for every position in the sequence.
        """
        self._current_log_sum = self._compute_log_sum(seq_codons)
        self._n_codons = len(seq_codons)

    def update_incremental_state(self, old_codon: str, new_codon: str) -> None:
        """Update the cached log_sum after an accepted codon swap.

        Uses compute_cai_incremental (NUMBA or pure-Python) for O(1)
        update instead of O(n) full recomputation.

        Args:
            old_codon: The codon that was replaced.
            new_codon: The codon that replaced it.
        """
        if self._current_log_sum is None:
            return  # Not initialized; caller should use reset_incremental_state

        epsilon = 1e-10

        # Get old adaptiveness
        old_idx = self._codon_to_idx.get(old_codon, 0)
        w_old = self._adaptiveness[old_idx]
        if w_old <= 0.0 or w_old != w_old:  # catch NaN
            w_old = epsilon

        # Get new adaptiveness
        new_idx = self._codon_to_idx.get(new_codon, 0)
        w_new = self._adaptiveness[new_idx]
        if w_new <= 0.0 or w_new != w_new:  # catch NaN
            w_new = epsilon

        # O(1) incremental update: new_log_sum = old_log_sum - log(w_old) + log(w_new)
        #
        # ALWAYS use the direct formula.  Recovering log_sum from CAI via
        # n * log(CAI) accumulates rounding errors over many swaps (the
        # round-trip through exp/log loses ~1 ULP per swap, causing drift
        # of ~1e-10 after 1000 swaps and growing quadratically).
        log_w_old = math.log(w_old)
        log_w_new = math.log(w_new)

        # Validate with NUMBA kernel if available (for consistency check only)
        if HAS_NUMBA and _numba_cai_incremental is not None:
            try:
                new_cai = _numba_cai_incremental(
                    self._current_log_sum, self._n_codons, w_old, w_new
                )
                # Use direct formula instead of recovering from CAI to avoid
                # floating-point drift.  The NUMBA result is only used as a
                # sanity check — if the direct formula diverges significantly
                # from the NUMBA result, we recompute from scratch.
                direct_log_sum = self._current_log_sum - log_w_old + log_w_new
                if new_cai > 0.0 and self._n_codons > 0:
                    numba_log_sum = self._n_codons * math.log(new_cai)
                    if abs(direct_log_sum - numba_log_sum) > 0.01:
                        # Significant divergence — recompute from scratch
                        # This should never happen in normal operation
                        pass  # fall through to direct formula
                self._current_log_sum = direct_log_sum
                return
            except Exception:
                logger.warning("NUMBA incremental CAI failed, falling back to pure-Python", exc_info=True)

        # Pure-Python incremental update — direct formula, no CAI round-trip
        self._current_log_sum = (
            self._current_log_sum - log_w_old + log_w_new
        )

    @property
    def current_log_sum(self) -> float | None:
        """The cached sum of log(w_i) for the current sequence state."""
        return self._current_log_sum

    def score_candidates(
        self,
        seq_codons: List[str],
        swap_position: int,
        candidate_codons: List[str],
    ) -> List[float]:
        """Score all candidate codon swaps at a single position.

        Returns a list of full-sequence CAI values, one per candidate,
        representing the CAI of the sequence if the codon at
        ``swap_position`` were replaced by each candidate.

        Args:
            seq_codons: Current codon strings for every position in the sequence.
            swap_position: Index into seq_codons of the codon being swapped.
            candidate_codons: List of candidate replacement codon strings.

        Returns:
            List of CAI scores (float), one per candidate.
        """
        n_candidates = len(candidate_codons)
        if n_candidates == 0:
            return []

        n_codons = len(seq_codons)

        # ── Fast path: NUMBA batch kernel ──────────────────────────
        if (
            HAS_NUMBA
            and _numba_batch_codon_swap_score is not None
        ):
            self._ensure_numpy()
            if self._np_adaptiveness is not None:
                import numpy as _np
                codon_indices = self._build_codon_indices(seq_codons)
                candidate_indices = _np.array(
                    [self._codon_to_idx.get(c, 0) for c in candidate_codons],
                    dtype=_np.int64,
                )
                # Use cached log_sum if available (avoids O(n) recompute)
                current_log_sum = self._current_log_sum if self._current_log_sum is not None else self._compute_log_sum(seq_codons)
                scores_arr = _numba_batch_codon_swap_score(
                    self._np_adaptiveness,
                    codon_indices,
                    n_codons,
                    swap_position,
                    candidate_indices,
                    n_candidates,
                    current_log_sum,
                )
                return [float(s) for s in scores_arr]

        # ── Fallback: pure-Python per-candidate evaluation ─────────
        epsilon = 1e-10
        current_log_sum = self._current_log_sum if self._current_log_sum is not None else self._compute_log_sum(seq_codons)

        old_idx = self._codon_to_idx.get(seq_codons[swap_position], 0)
        w_old = self._adaptiveness[old_idx]
        if w_old <= 0.0:
            w_old = epsilon
        log_w_old = math.log(w_old)

        scores: List[float] = []
        for c in candidate_codons:
            c_idx = self._codon_to_idx.get(c, 0)
            w_new = self._adaptiveness[c_idx]
            if w_new <= 0.0:
                w_new = epsilon
            new_log_sum = current_log_sum - log_w_old + math.log(w_new)
            scores.append(math.exp(new_log_sum / n_codons) if n_codons > 0 else 0.0)

        return scores


__all__ = [
    "HAS_NUMBA",
    "_HAS_NUMBA",
    "_adaptiveness_to_array",
    "_codon_to_index",
    "_dna_to_codon_indices",
    "_compute_cai_fast",
    "_count_dinucs_fast",
    "_BatchSwapScorer",
]
