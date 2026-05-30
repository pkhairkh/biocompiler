"""
MaxEntScan Splice Site Scoring Module

Implements the Maximum Entropy Modeling approach for scoring donor and acceptor
splice sites, based on Yeo & Burge (2004) "Maximum Entropy Modeling of Short
Sequence Motifs".

Uses position weight matrices (PWMs) with approximate but biologically
realistic values that capture the consensus patterns of human splice sites.

All scoring is DETERMINISTIC: same input always produces same output.
"""

import math
from typing import List, Tuple, Dict

# Base encoding
BASE_TO_INDEX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
INDEX_TO_BASE: Dict[int, str] = {0: "A", 1: "C", 2: "G", 3: "T"}

# Background probability (uniform for DNA)
BG_PROB = 0.25

# ---------------------------------------------------------------------------
# Donor site PWM (9 positions)
# Positions correspond to: -3, -2, -1, +1, +2, +3, +4, +5, +6
# relative to the intron start (the GT dinucleotide is at +1, +2).
# Each row is [P(A), P(C), P(G), P(T)] for that position.
# ---------------------------------------------------------------------------
DONOR_PWM: List[List[float]] = [
    # Position -3: weak consensus
    [0.30, 0.30, 0.20, 0.20],
    # Position -2: weak consensus
    [0.30, 0.30, 0.20, 0.20],
    # Position -1: C or T preferred
    [0.10, 0.50, 0.10, 0.30],
    # Position +1: almost always G (first base of GT)
    [0.01, 0.01, 0.98, 0.00],
    # Position +2: almost always T (second base of GT)
    [0.00, 0.00, 0.01, 0.99],
    # Position +3: mild G preference
    [0.30, 0.20, 0.30, 0.20],
    # Position +4: near uniform
    [0.25, 0.25, 0.25, 0.25],
    # Position +5: G and T preference
    [0.15, 0.15, 0.40, 0.30],
    # Position +6: mild G preference
    [0.20, 0.20, 0.30, 0.30],
]

# ---------------------------------------------------------------------------
# Acceptor site PWM (23 positions)
# Positions 0-19: upstream context (polypyrimidine tract, intronic)
# Positions 20, 21: AG dinucleotide (end of intron)
# Position 22: first exonic base
# ---------------------------------------------------------------------------
ACCEPTOR_PWM: List[List[float]] = [
    # Positions 0-3: distant upstream (weak pyrimidine bias)
    [0.20, 0.30, 0.20, 0.30],
    [0.20, 0.30, 0.15, 0.35],
    [0.20, 0.25, 0.20, 0.35],
    [0.20, 0.30, 0.15, 0.35],
    # Positions 4-7: mid polypyrimidine tract (stronger C/T bias)
    [0.15, 0.35, 0.10, 0.40],
    [0.15, 0.30, 0.15, 0.40],
    [0.15, 0.35, 0.10, 0.40],
    [0.10, 0.40, 0.10, 0.40],
    # Positions 8-11: strong polypyrimidine tract
    [0.10, 0.40, 0.10, 0.40],
    [0.10, 0.40, 0.10, 0.40],
    [0.10, 0.45, 0.05, 0.40],
    [0.10, 0.45, 0.05, 0.40],
    # Positions 12-15: proximal polypyrimidine tract
    [0.10, 0.45, 0.05, 0.40],
    [0.10, 0.45, 0.05, 0.40],
    [0.05, 0.50, 0.05, 0.40],
    [0.05, 0.50, 0.05, 0.40],
    # Positions 16-18: just upstream of AG
    [0.05, 0.55, 0.05, 0.35],
    [0.05, 0.55, 0.05, 0.35],
    [0.08, 0.60, 0.02, 0.30],
    # Position 19: the "Y" before AG (pyrimidine)
    [0.02, 0.65, 0.01, 0.32],
    # Position 20: A of AG (almost always A)
    [0.98, 0.01, 0.005, 0.005],
    # Position 21: G of AG (almost always G)
    [0.005, 0.01, 0.98, 0.005],
    # Position 22: first exonic base (weak preference for G/C)
    [0.25, 0.30, 0.30, 0.15],
]


def _log2(x: float) -> float:
    """Compute log2, returning a large negative number for x <= 0."""
    if x <= 0:
        return -50.0
    return math.log2(x)


def score_donor(seq: str, position: int) -> float:
    """
    Score a 9-mer at the given position as a donor splice site.

    The donor site is a 9-mer spanning positions -3 to +6 relative to the
    GT dinucleotide. If GT starts at `position`, then the 9-mer is
    seq[position-3 : position+6].

    The score is the log-odds ratio: sum of log2(P_i / 0.25) for each
    position, where P_i is the PWM probability for the observed base.

    Args:
        seq: DNA sequence (uppercase, A/C/G/T only)
        position: index of the G in the GT dinucleotide

    Returns:
        Log-odds score. Real MaxEntScan scores range from ~-10 to ~+12.
        Scores above ~3 indicate strong donor sites.
    """
    seq = seq.upper()
    start = position - 3
    end = position + 6  # exclusive; 9 bases from start to end

    if start < 0 or end > len(seq):
        return -50.0  # insufficient context

    score = 0.0
    for pwm_idx in range(9):
        seq_idx = start + pwm_idx
        if seq_idx >= len(seq):
            return -50.0
        base = seq[seq_idx]
        if base not in BASE_TO_INDEX:
            return -50.0  # invalid character
        base_idx = BASE_TO_INDEX[base]
        prob = DONOR_PWM[pwm_idx][base_idx]
        score += _log2(prob / BG_PROB)

    return score


def score_acceptor(seq: str, position: int) -> float:
    """
    Score a 23-mer at the given position as an acceptor splice site.

    The acceptor site is a 23-mer where the AG dinucleotide is at positions
    20-21 of the model. If A of AG is at `position`, the 23-mer spans
    seq[position-20 : position+3].

    The score is the log-odds ratio: sum of log2(P_i / 0.25) for each
    position.

    Args:
        seq: DNA sequence (uppercase, A/C/G/T only)
        position: index of the A in the AG dinucleotide

    Returns:
        Log-odds score. Real MaxEntScan scores range from ~-12 to ~+12.
        Scores above ~3 indicate strong acceptor sites.
    """
    seq = seq.upper()
    start = position - 20
    end = position + 3  # exclusive; 23 bases from start to end

    if start < 0 or end > len(seq):
        return -50.0  # insufficient context

    score = 0.0
    for pwm_idx in range(23):
        seq_idx = start + pwm_idx
        if seq_idx >= len(seq):
            return -50.0
        base = seq[seq_idx]
        if base not in BASE_TO_INDEX:
            return -50.0  # invalid character
        base_idx = BASE_TO_INDEX[base]
        prob = ACCEPTOR_PWM[pwm_idx][base_idx]
        score += _log2(prob / BG_PROB)

    return score


def scan_splice_sites(
    seq: str,
    donor_threshold: float = 3.0,
    acceptor_threshold: float = 3.0,
) -> List[Tuple[int, str, float]]:
    """
    Scan an entire DNA sequence for donor and acceptor splice sites.

    For donor sites: searches for GT dinucleotides and scores each as a
    potential 5' splice site using the 9-mer donor model.

    For acceptor sites: searches for AG dinucleotides and scores each as a
    potential 3' splice site using the 23-mer acceptor model.

    Args:
        seq: DNA sequence (uppercase or lowercase A/C/G/T)
        donor_threshold: minimum score to report a donor site (default 3.0)
        acceptor_threshold: minimum score to report an acceptor site (default 3.0)

    Returns:
        List of (position, site_type, score) tuples, sorted by position.
        site_type is 'donor' or 'acceptor'.
        position is the index of the first base of the dinucleotide (G for GT, A for AG).
    """
    seq = seq.upper()
    results: List[Tuple[int, str, float]] = []

    if len(seq) < 2:
        return results

    # Validate sequence
    valid_bases = set("ACGT")
    if not all(b in valid_bases for b in seq):
        # Skip invalid characters but still scan what we can
        pass

    # Scan for donor sites (GT dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "T":
            score = score_donor(seq, i)
            if score >= donor_threshold:
                results.append((i, "donor", round(score, 4)))

    # Scan for acceptor sites (AG dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i] == "A" and seq[i + 1] == "G":
            score = score_acceptor(seq, i)
            if score >= acceptor_threshold:
                results.append((i, "acceptor", round(score, 4)))

    # Sort by position
    results.sort(key=lambda x: x[0])
    return results


def max_donor_score(seq: str) -> float:
    """
    Return the maximum donor splice site score across the entire sequence.

    Useful for quickly checking if any cryptic donor site exists.
    Returns -50.0 if no scorable position exists.
    """
    seq = seq.upper()
    best = -50.0
    for i in range(len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "T":
            s = score_donor(seq, i)
            if s > best:
                best = s
    return round(best, 4)


def max_acceptor_score(seq: str) -> float:
    """
    Return the maximum acceptor splice site score across the entire sequence.

    Useful for quickly checking if any cryptic acceptor site exists.
    Returns -50.0 if no scorable position exists.
    """
    seq = seq.upper()
    best = -50.0
    for i in range(len(seq) - 1):
        if seq[i] == "A" and seq[i + 1] == "G":
            s = score_acceptor(seq, i)
            if s > best:
                best = s
    return round(best, 4)


if __name__ == "__main__":
    # Quick self-test with known splice site sequences
    # Canonical donor: CAG|GTAAGT -> 9-mer around GT
    test_donor = "CAGGTAAGT"
    print(f"Test donor sequence: {test_donor}")
    print(f"  Donor score (GT at pos 3): {score_donor(test_donor, 3):.4f}")

    # Canonical acceptor: (poly-py tract)TTTTTTTTTTTTTTTTTCAG|G
    test_acceptor = "TTTTTTTTTTTTTTTTTCAGG"
    print(f"\nTest acceptor sequence: {test_acceptor}")
    print(f"  Acceptor score (AG at pos 18): {score_acceptor(test_acceptor, 18):.4f}")

    # Scan a test sequence with both sites
    test_seq = "CAGGTAAGTNNNNNNNTTTTTTTTTTTTTTTTTCAGGATGG"
    print(f"\nScanning test sequence: {test_seq}")
    sites = scan_splice_sites(test_seq)
    for pos, stype, score in sites:
        print(f"  Position {pos}: {stype} site, score={score:.4f}")
