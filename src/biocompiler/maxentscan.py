"""
BioCompiler MaxEntScan — Splice Site Scoring

Position weight matrix implementation based on Yeo & Burge (2004).
Deterministic: same input always produces same output.
"""

import math
import logging
from typing import List, Tuple, Dict

logger = logging.getLogger(__name__)

BASE_TO_INDEX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
BG_PROB = 0.25

# Donor site PWM (9 positions: -3 to +6 relative to GT)
DONOR_PWM: List[List[float]] = [
    [0.30, 0.30, 0.20, 0.20],  # -3
    [0.30, 0.30, 0.20, 0.20],  # -2
    [0.10, 0.50, 0.10, 0.30],  # -1: C/T preferred
    [0.01, 0.01, 0.98, 0.00],  # +1: G invariant
    [0.00, 0.00, 0.01, 0.99],  # +2: T invariant
    [0.30, 0.20, 0.30, 0.20],  # +3
    [0.25, 0.25, 0.25, 0.25],  # +4
    [0.15, 0.15, 0.40, 0.30],  # +5: G/T preference
    [0.20, 0.20, 0.30, 0.30],  # +6
]

# Acceptor site PWM (23 positions)
ACCEPTOR_PWM: List[List[float]] = [
    [0.20, 0.30, 0.20, 0.30],  # 0
    [0.20, 0.30, 0.15, 0.35],  # 1
    [0.20, 0.25, 0.20, 0.35],  # 2
    [0.20, 0.30, 0.15, 0.35],  # 3
    [0.15, 0.35, 0.10, 0.40],  # 4
    [0.15, 0.30, 0.15, 0.40],  # 5
    [0.15, 0.35, 0.10, 0.40],  # 6
    [0.10, 0.40, 0.10, 0.40],  # 7
    [0.10, 0.40, 0.10, 0.40],  # 8
    [0.10, 0.40, 0.10, 0.40],  # 9
    [0.10, 0.45, 0.05, 0.40],  # 10
    [0.10, 0.45, 0.05, 0.40],  # 11
    [0.10, 0.45, 0.05, 0.40],  # 12
    [0.10, 0.45, 0.05, 0.40],  # 13
    [0.05, 0.50, 0.05, 0.40],  # 14
    [0.05, 0.50, 0.05, 0.40],  # 15
    [0.05, 0.55, 0.05, 0.35],  # 16
    [0.05, 0.55, 0.05, 0.35],  # 17
    [0.08, 0.60, 0.02, 0.30],  # 18
    [0.02, 0.65, 0.01, 0.32],  # 19
    [0.98, 0.01, 0.005, 0.005], # 20: A invariant
    [0.005, 0.01, 0.98, 0.005], # 21: G invariant
    [0.25, 0.30, 0.30, 0.15],  # 22: exonic
]


def _log2(x: float) -> float:
    if x <= 0:
        return -50.0
    return math.log2(x)


def score_donor(seq: str, position: int) -> float:
    """Score a 9-mer at position as a donor splice site."""
    seq = seq.upper()
    start = position - 3
    end = position + 6
    if start < 0 or end > len(seq):
        return -50.0
    score = 0.0
    for pwm_idx in range(9):
        base = seq[start + pwm_idx]
        if base not in BASE_TO_INDEX:
            return -50.0
        prob = DONOR_PWM[pwm_idx][BASE_TO_INDEX[base]]
        score += _log2(prob / BG_PROB)
    return score


def score_acceptor(seq: str, position: int) -> float:
    """Score a 23-mer at position as an acceptor splice site."""
    seq = seq.upper()
    start = position - 20
    end = position + 3
    if start < 0 or end > len(seq):
        return -50.0
    score = 0.0
    for pwm_idx in range(23):
        base = seq[start + pwm_idx]
        if base not in BASE_TO_INDEX:
            return -50.0
        prob = ACCEPTOR_PWM[pwm_idx][BASE_TO_INDEX[base]]
        score += _log2(prob / BG_PROB)
    return score


def scan_splice_sites(
    seq: str,
    donor_threshold: float = 3.0,
    acceptor_threshold: float = 3.0,
) -> List[Tuple[int, str, float]]:
    """Scan entire sequence for donor and acceptor splice sites above threshold."""
    seq = seq.upper()
    results: List[Tuple[int, str, float]] = []
    for i in range(len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "T":
            s = score_donor(seq, i)
            if s >= donor_threshold:
                results.append((i, "donor", round(s, 4)))
    for i in range(len(seq) - 1):
        if seq[i] == "A" and seq[i + 1] == "G":
            s = score_acceptor(seq, i)
            if s >= acceptor_threshold:
                results.append((i, "acceptor", round(s, 4)))
    results.sort(key=lambda x: x[0])
    return results


def max_donor_score(seq: str) -> float:
    """Return the maximum donor splice site score across the sequence."""
    seq = seq.upper()
    best = -50.0
    for i in range(len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "T":
            s = score_donor(seq, i)
            if s > best:
                best = s
    return round(best, 4)


def max_acceptor_score(seq: str) -> float:
    """Return the maximum acceptor splice site score across the sequence."""
    seq = seq.upper()
    best = -50.0
    for i in range(len(seq) - 1):
        if seq[i] == "A" and seq[i + 1] == "G":
            s = score_acceptor(seq, i)
            if s > best:
                best = s
    return round(best, 4)
