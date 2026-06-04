"""
BioCompiler ViennaRNA Offline Fallback — Nussinov-Based RNA Structure Prediction
==================================================================================

Provides a lightweight, offline-capable fallback when the ViennaRNA Python
bindings (``import RNA``) are not installed.  Instead of failing silently,
this module uses the classic **Nussinov–Jacobson** dynamic programming
algorithm to predict RNA secondary structure and approximates free-energy
values using simplified nearest-neighbor parameters.

**IMPORTANT — Accuracy Caveat**
    This fallback is *significantly* less accurate than ViennaRNA.
    The Nussinov algorithm maximizes the *number* of base pairs rather than
    minimizing free energy.  The ΔG estimates are crude approximations based
    on per-pair enthalpy/entropy contributions and lack the stacking and
    loop parameters that make ViennaRNA accurate.  All results are flagged
    with ``method="nussinov_fallback"`` and include a ``warning`` field.

Algorithm Components
--------------------
1. **Nussinov–Jacobson DP**: O(n³) dynamic programming that maximises
   base-pair count with a minimum hairpin loop length of 3 nt.
   Produces a dot-bracket secondary structure string via traceback.

2. **Approximate ΔG**: Per-pair energy contributions:
   - GC pair: −3.4 kcal/mol
   - AU pair: −2.1 kcal/mol
   - GU wobble: −1.4 kcal/mol
   These are rough averages of nearest-neighbor stacking parameters;
   stacking bonuses and loop penalties are *not* included.

3. **MFE prediction**: Wraps Nussinov + ΔG into an MFEResult dataclass
   matching the interface of ``viennarna.predict_mfe``.

4. **Accessibility prediction**: Counts base-pairing partners for each
   nucleotide across a sliding window; unpaired positions get higher
   accessibility scores.

5. **Stable structure detection**: Identifies stem-loop structures whose
   estimated ΔG falls below a user-specified threshold.

When to Use
-----------
This fallback is invoked automatically by the ViennaRNA wrapper when
``import RNA`` fails.  It should **never** be used as a replacement for
ViennaRNA when the library is available.

References
----------
- Nussinov & Jacobson, PNAS 1980; 77:6309–6313 (base-pair maximisation)
- Zuker, Nucleic Acids Res 2003; 31:3406–3415 (MFE / nearest-neighbor)
- Lorenz et al., J Chem Inf Model 2011; 51:2547–2557 (ViennaRNA Package)
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    "nussinov_fold",
    "compute_approx_dg",
    "predict_mfe_fallback",
    "predict_accessibility_fallback",
    "find_stable_structures_fallback",
    "MFEResult",
    "AccessibilityResult",
    "StableStructure",
]


# ==============================================================================
# Constants
# ==============================================================================

#: Minimum hairpin loop length (Nussinov constraint).
MIN_LOOP_LENGTH: int = 3

#: Approximate per-pair free-energy contributions (kcal/mol).
#: These are simplified averages — real nearest-neighbor parameters
#: depend on stacking context, loop type, and sequence.
PAIR_ENERGIES: dict[tuple[str, str], float] = {
    ("G", "C"): -3.4,
    ("C", "G"): -3.4,
    ("A", "U"): -2.1,
    ("U", "A"): -2.1,
    ("G", "U"): -1.4,
    ("U", "G"): -1.4,
}

#: RNA complement map (DNA T is treated as U).
_COMPLEMENT: dict[str, str] = {
    "A": "U", "U": "A", "G": "C", "C": "G",
    "T": "A",  # DNA input support
}

#: Default sliding window size for accessibility computation.
DEFAULT_WINDOW_SIZE: int = 120

#: Default ΔG threshold for "stable" structures (kcal/mol).
DEFAULT_STABLE_DG_THRESHOLD: float = -4.0

#: Maximum length for Nussinov O(n³) DP before switching to heuristic.
MAX_NUSSINOV_LENGTH: int = 2000


# ==============================================================================
# Data classes
# ==============================================================================

@dataclass
class MFEResult:
    """Result of minimum free energy structure prediction.

    Mirrors the interface of ``viennarna.predict_mfe`` so that downstream
    code can treat both paths uniformly.

    Attributes:
        sequence:  Input RNA/DNA sequence (upper-case).
        structure: Dot-bracket secondary structure string.
        mfe:       Estimated minimum free energy (kcal/mol).
        method:    Always ``"nussinov_fallback"``.
        warning:   User-facing warning about approximate results.
    """
    sequence: str
    structure: str
    mfe: float
    method: str = "nussinov_fallback"
    warning: str = "Approximate results from Nussinov fallback; install ViennaRNA for accurate prediction"


@dataclass
class AccessibilityResult:
    """Result of accessibility (unpaired probability) prediction.

    Accessibility is the probability that a nucleotide is unpaired,
    estimated by the fraction of windows in which the position is
    not base-paired in the Nussinov prediction.

    Attributes:
        sequence:       Input RNA/DNA sequence (upper-case).
        accessibility:  Per-nucleotide accessibility values in [0, 1].
        window_size:    Sliding window length used.
        method:         Always ``"nussinov_fallback"``.
        warning:        User-facing warning about approximate results.
    """
    sequence: str
    accessibility: List[float]
    window_size: int
    method: str = "nussinov_fallback"
    warning: str = "Approximate accessibility from Nussinov fallback; install ViennaRNA for accurate prediction"


@dataclass
class StableStructure:
    """A predicted stable stem-loop structure.

    Attributes:
        start:      0-based start position of the structure.
        end:        0-based end position (exclusive).
        structure:  Dot-bracket string for this region.
        dg:         Estimated free energy (kcal/mol).
        stem_length: Number of base pairs in the stem.
    """
    start: int
    end: int
    structure: str
    dg: float
    stem_length: int


# ==============================================================================
# Helper: normalise input
# ==============================================================================

def _normalise_rna(seq: str) -> str:
    """Upper-case and convert T → U for RNA pairing logic."""
    return seq.upper().replace("T", "U")


def _can_pair(a: str, b: str) -> bool:
    """Check if two RNA bases can form a canonical or wobble pair."""
    return (a, b) in PAIR_ENERGIES


# ==============================================================================
# Nussinov–Jacobson DP
# ==============================================================================

def nussinov_fold(seq: str, min_loop: int = MIN_LOOP_LENGTH) -> Tuple[str, float]:
    """Predict RNA secondary structure using the Nussinov–Jacobson algorithm.

    The Nussinov algorithm fills an n×n DP table where ``dp[i][j]`` is the
    maximum number of base pairs achievable in the subsequence ``seq[i:j+1]``.
    A traceback then recovers the optimal pairing, which is encoded as a
    dot-bracket string.

    Free energy is estimated from the identified pairs using the simplified
    per-pair energies in ``PAIR_ENERGIES``.

    Args:
        seq:      RNA or DNA sequence (case-insensitive; T is treated as U).
        min_loop: Minimum number of unpaired nucleotides in a hairpin loop
                  (default 3).

    Returns:
        ``(structure, dg)`` where *structure* is a dot-bracket string of the
        same length as *seq* and *dg* is the approximate free energy in
        kcal/mol (negative = stable).
    """
    rna = _normalise_rna(seq)
    n = len(rna)

    if n == 0:
        return "", 0.0

    # --- DP table ---
    dp = [[0] * n for _ in range(n)]

    for length in range(min_loop + 1, n):  # span length
        for i in range(n - length):
            j = i + length
            # Option 1: j is unpaired
            best = dp[i][j - 1]
            # Option 2: j pairs with some k in [i, j-min_loop)
            for k in range(i, j - min_loop + 1):
                if _can_pair(rna[k], rna[j]):
                    left = dp[i][k - 1] if k > i else 0
                    pair_score = left + 1 + (dp[k + 1][j - 1] if k + 1 <= j - 1 else 0)
                    if pair_score > best:
                        best = pair_score
            dp[i][j] = best

    # --- Traceback ---
    structure = ["."] * n
    _traceback(dp, rna, 0, n - 1, structure, min_loop)

    dot_bracket = "".join(structure)

    # --- Approximate ΔG ---
    dg = compute_approx_dg(seq, dot_bracket)

    return dot_bracket, dg


def _traceback(
    dp: List[List[int]],
    rna: str,
    i: int,
    j: int,
    structure: List[str],
    min_loop: int,
) -> None:
    """Recursive traceback for Nussinov DP table."""
    if i >= j or j - i < min_loop + 1:
        return
    if dp[i][j] == dp[i][j - 1]:
        # j is unpaired
        _traceback(dp, rna, i, j - 1, structure, min_loop)
        return
    # j pairs with some k
    for k in range(i, j - min_loop + 1):
        if _can_pair(rna[k], rna[j]):
            left = dp[i][k - 1] if k > i else 0
            inner = dp[k + 1][j - 1] if k + 1 <= j - 1 else 0
            if dp[i][j] == left + 1 + inner:
                structure[k] = "("
                structure[j] = ")"
                if k > i:
                    _traceback(dp, rna, i, k - 1, structure, min_loop)
                if k + 1 <= j - 1:
                    _traceback(dp, rna, k + 1, j - 1, structure, min_loop)
                return


# ==============================================================================
# Approximate ΔG
# ==============================================================================

def compute_approx_dg(seq: str, structure: str) -> float:
    """Compute an approximate free energy from a sequence and dot-bracket structure.

    Sums per-pair energy contributions for each base pair identified in the
    dot-bracket notation.  Stacking bonuses and loop penalties are *not*
    included — this is a very rough estimate.

    Args:
        seq:        RNA or DNA sequence (case-insensitive; T → U).
        structure:  Dot-bracket secondary structure string of same length.

    Returns:
        Approximate free energy in kcal/mol (negative = stable).
    """
    rna = _normalise_rna(seq)
    if len(rna) != len(structure):
        raise ValueError(
            f"Sequence length ({len(rna)}) != structure length ({len(structure)})"
        )

    # Identify base pairs from dot-bracket using a stack
    stack: List[int] = []
    pairs: List[Tuple[int, int]] = []
    for idx, char in enumerate(structure):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if stack:
                open_idx = stack.pop()
                pairs.append((open_idx, idx))

    # Sum per-pair energies
    dg = 0.0
    for open_idx, close_idx in pairs:
        a = rna[open_idx]
        b = rna[close_idx]
        pair_energy = PAIR_ENERGIES.get((a, b), 0.0)
        dg += pair_energy

    return dg


# ==============================================================================
# predict_mfe_fallback
# ==============================================================================

def predict_mfe_fallback(seq: str) -> MFEResult:
    """Fallback MFE prediction matching the ``viennarna.predict_mfe`` interface.

    Uses the Nussinov algorithm internally and wraps the result in an
    :class:`MFEResult` dataclass with appropriate warning fields.

    Args:
        seq: RNA or DNA sequence (case-insensitive).

    Returns:
        MFEResult with structure, mfe, method, and warning fields.
    """
    seq_upper = seq.upper()
    structure, dg = nussinov_fold(seq_upper)
    return MFEResult(
        sequence=seq_upper,
        structure=structure,
        mfe=dg,
    )


# ==============================================================================
# predict_accessibility_fallback
# ==============================================================================

def predict_accessibility_fallback(
    seq: str,
    window_size: int = DEFAULT_WINDOW_SIZE,
) -> AccessibilityResult:
    """Fallback accessibility prediction using Nussinov folding.

    For each nucleotide position *i*, the accessibility is estimated as:

        acc[i] = 1.0 − (paired_count[i] / total_windows_containing_i)

    where ``paired_count[i]`` is the number of sliding windows in which
    position *i* is paired (``(`` or ``)`` in the dot-bracket) and
    ``total_windows_containing_i`` is the total number of windows that
    include position *i*.

    For sequences shorter than *window_size*, a single fold of the entire
    sequence is used, so accessibility is 1.0 for unpaired and 0.0 for
    paired positions.

    Args:
        seq:         RNA or DNA sequence (case-insensitive).
        window_size: Sliding window length (default 120).

    Returns:
        AccessibilityResult with per-nucleotide values in [0, 1].
    """
    seq_upper = seq.upper()
    n = len(seq_upper)

    if n == 0:
        return AccessibilityResult(
            sequence="",
            accessibility=[],
            window_size=window_size,
        )

    # For short sequences, fold the whole thing once
    if n <= window_size:
        structure, _ = nussinov_fold(seq_upper)
        accessibility = [
            0.0 if ch in ("(", ")") else 1.0
            for ch in structure
        ]
        return AccessibilityResult(
            sequence=seq_upper,
            accessibility=accessibility,
            window_size=window_size,
        )

    # Sliding window approach
    paired_count = [0] * n
    window_count = [0] * n

    for start in range(0, n - window_size + 1, max(1, window_size // 4)):
        end = start + window_size
        window_seq = seq_upper[start:end]
        structure, _ = nussinov_fold(window_seq)

        for local_idx, ch in enumerate(structure):
            global_idx = start + local_idx
            window_count[global_idx] += 1
            if ch in ("(", ")"):
                paired_count[global_idx] += 1

    # Fill positions that were not covered by any window
    for i in range(n):
        if window_count[i] == 0:
            # Edge positions — fold the nearest window
            ws = max(0, min(i - window_size // 2, n - window_size))
            we = ws + window_size
            window_seq = seq_upper[ws:we]
            structure, _ = nussinov_fold(window_seq)
            local_idx = i - ws
            if 0 <= local_idx < len(structure):
                accessibility_val = 0.0 if structure[local_idx] in ("(", ")") else 1.0
            else:
                accessibility_val = 1.0
            # Store directly
            paired_count[i] = 0 if accessibility_val == 1.0 else 1
            window_count[i] = 1

    accessibility = [
        round(1.0 - (paired_count[i] / window_count[i]), 4)
        for i in range(n)
    ]

    return AccessibilityResult(
        sequence=seq_upper,
        accessibility=accessibility,
        window_size=window_size,
    )


# ==============================================================================
# find_stable_structures_fallback
# ==============================================================================

def find_stable_structures_fallback(
    seq: str,
    dg_threshold: float = DEFAULT_STABLE_DG_THRESHOLD,
    window_size: int = 60,
    min_stem: int = 3,
) -> List[StableStructure]:
    """Find stable stem-loop structures in an RNA sequence.

    Slides a window across the sequence, folds each window with the
    Nussinov algorithm, and reports regions whose estimated ΔG is
    more negative than *dg_threshold*.

    Args:
        seq:          RNA or DNA sequence (case-insensitive).
        dg_threshold: Maximum ΔG (most negative) for a structure to be
                      considered "stable" (default −4.0 kcal/mol).
        window_size:  Sliding window length (default 60 nt).
        min_stem:     Minimum number of base pairs in the stem for
                      reporting (default 3).

    Returns:
        List of :class:`StableStructure` instances sorted by position.
    """
    seq_upper = seq.upper()
    n = len(seq_upper)

    if n == 0:
        return []

    results: List[StableStructure] = []

    # If sequence is shorter than window_size, fold the whole thing
    if n <= window_size:
        structure, dg = nussinov_fold(seq_upper)
        stem_length = _count_stem_pairs(structure)
        if dg <= dg_threshold and stem_length >= min_stem:
            results.append(StableStructure(
                start=0,
                end=n,
                structure=structure,
                dg=dg,
                stem_length=stem_length,
            ))
        return results

    # Slide window with 50% overlap
    step = max(1, window_size // 2)
    seen_ranges: set[tuple[int, int]] = set()

    for start in range(0, n - window_size + 1, step):
        end = start + window_size
        window_seq = seq_upper[start:end]
        structure, dg = nussinov_fold(window_seq)
        stem_length = _count_stem_pairs(structure)

        if dg <= dg_threshold and stem_length >= min_stem:
            # Avoid duplicating overlapping hits
            range_key = (start, end)
            if range_key not in seen_ranges:
                seen_ranges.add(range_key)
                results.append(StableStructure(
                    start=start,
                    end=end,
                    structure=structure,
                    dg=dg,
                    stem_length=stem_length,
                ))

    # Sort by start position
    results.sort(key=lambda s: s.start)
    return results


def _count_stem_pairs(structure: str) -> int:
    """Count the number of base pairs (open parentheses) in a dot-bracket."""
    return structure.count("(")
