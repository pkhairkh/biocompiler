"""
BioCompiler Splice Site Scoring
================================
MaxEntScan-based scoring for cryptic splice site detection.
"""

import math
from typing import Optional, List, Tuple

# Position weight matrix for 5' splice site (positions -3 to +6 relative to GT)
_MAXENT_PWM = [
    # pos -3   -2   -1    G    T   +1   +2   +3   +4
    [0.10, 0.07, 0.12, 3.50, 3.50, 0.16, 0.20, 0.12, 0.08],  # A
    [0.06, 0.04, 0.06, 0.01, 0.01, 0.06, 0.04, 0.06, 0.05],  # C
    [0.06, 0.04, 0.06, 0.01, 0.01, 0.06, 0.04, 0.06, 0.05],  # G
    [0.08, 0.15, 0.10, 0.01, 0.01, 0.14, 0.10, 0.12, 0.08],  # T
]

_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def maxent_score(context: str) -> float:
    """Compute MaxEntScan score for a potential splice site context.

    Args:
        context: DNA sequence around a GT dinucleotide (ideally 9-mer)

    Returns:
        MaxEntScan log-odds score. Higher = stronger splice signal.
        Thresholds: < 3.0 PASS, 3.0-6.0 UNCERTAIN, >= 6.0 FAIL.
    """
    if len(context) < 4:
        return 0.0

    if len(context) < 9:
        context = "A" * (9 - len(context)) + context
        context = context[-9:]

    score = 0.0
    for pos in range(min(len(context), 9)):
        base = context[pos].upper() if pos < len(context) else "A"
        idx = _BASE_INDEX.get(base, 0)
        score += _MAXENT_PWM[idx][pos]

    return score


def score_splice_sites(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0):
    """Score all potential splice sites in a sequence."""
    from .type_system import SpliceVerdict
    results: List[Tuple[int, float, SpliceVerdict]] = []
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            start = max(0, i - 3)
            end = min(len(seq), i + 6)
            context = seq[start:end]
            sc = maxent_score(context)
            if sc < low_thresh:
                verdict = SpliceVerdict.PASS
            elif sc < high_thresh:
                verdict = SpliceVerdict.UNCERTAIN
            else:
                verdict = SpliceVerdict.FAIL
            results.append((i, sc, verdict))
    return results
