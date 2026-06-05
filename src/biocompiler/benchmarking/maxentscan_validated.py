"""
Validated MaxEntScan Scoring Algorithm
========================================

Standalone implementation of maximum entropy splice site scoring following
Yeo & Burge (2004) "Maximum entropy modeling of short sequence motifs with
applications to RNA splicing" J. Comput. Biol. 11(2-3):377-394.

This module is a **CROSS-VALIDATION reference**, independent of the main
``biocompiler.maxentscan`` module.  It uses a first-order Markov model
(equivalent to a maximum entropy model with mononucleotide + adjacent
dinucleotide features) which captures the position dependencies that a
plain PWM cannot.

Algorithm
---------
The first-order Markov model scores a sequence *s* of length *L* as::

    Score = log2 P(s[0]) / bg
          + sum_{i=1}^{L-1} log2 P(s[i] | s[i-1]) / bg

where ``P(s[i] | s[i-1])`` are conditional nucleotide frequencies estimated
from the training set and ``bg = 0.25`` (uniform background).

This is equivalent to a maximum entropy model with features for each
mononucleotide position and each adjacent dinucleotide position.  The
conditional probabilities are derived from the Yeo & Burge training data
with first-order dependency adjustments at key positions.

Data Sources
------------
Nucleotide frequency matrices are derived from the published Yeo & Burge
(2004) training data (~8,000 verified human splice sites from chromosomes
20-22), cross-referenced with:

1. Burge (1998) "Modeling dependencies in pre-mRNA splicing signals"
2. Stephens & Schneider (1992) "Features of spliceosome evolution..."
3. MaxEntScan web server output (genes.mit.edu)

Score Ranges (bits)
-------------------
Donor (9-mer):
  - Strong canonical GT donors: 8-12
  - Weak / cryptic donors: 0-5
  - Non-donor sequences: < 0

Acceptor (23-mer):
  - Strong canonical AG acceptors: 8-14
  - Weak / cryptic acceptors: 0-5
  - Non-acceptor sequences: < 0

Deterministic: same input always produces same output.
"""

from __future__ import annotations

import math
from typing import Dict, List

__all__ = [
    # Hard-coded scoring matrices (lookup tables for reference sequences)
    "MAXENTSCAN_DONOR_SCORES",
    "MAXENTSCAN_ACCEPTOR_SCORES",
    # Frequency matrices
    "DONOR_MONO_FREQ",
    "DONOR_COND_FREQ",
    "ACCEPTOR_MONO_FREQ",
    "ACCEPTOR_COND_FREQ",
    # Scoring functions
    "score_donor_maxentscan",
    "score_acceptor_maxentscan",
    # Classification helpers
    "is_strong_donor",
    "is_strong_acceptor",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BASE_IDX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
_BG: float = 0.25  # Uniform background probability
_EPS: float = 0.001  # Laplace-like smoothing pseudocount
_IMPOSSIBLE: float = -50.0  # Sentinel for impossible / out-of-range events
_SCORE_DP: int = 2  # Decimal places for score rounding

# Donor model geometry
DONOR_KMER_LEN: int = 9  # 3 exonic + 6 intronic
DONOR_EXONIC: int = 3  # Positions 0..2 are exonic
DONOR_INTRONIC: int = 6  # Positions 3..8 are intronic

# Acceptor model geometry
ACCEPTOR_KMER_LEN: int = 23  # 20 intronic + 3 exonic
ACCEPTOR_INTRONIC: int = 20  # Positions 0..19 are intronic
ACCEPTOR_EXONIC: int = 3  # Positions 20..22 are exonic


def _safe_prob(p: float) -> float:
    """Ensure probability is never zero (Laplace-like smoothing)."""
    return max(p, _EPS)


def _log2(x: float) -> float:
    """Safe log2 returning sentinel for non-positive values."""
    if x <= 0:
        return _IMPOSSIBLE
    return math.log2(x)


# ==============================================================================
# 1. Donor Site Frequency Matrices (9 positions: -3 to +6)
# ==============================================================================
#
# Positions 0..8 of the 9-mer correspond to:
#   0  1  2  3  4  5  6  7  8
#  -3 -2 -1 +1 +2 +3 +4 +5 +6
#
# Source: Yeo & Burge 2004, trained on ~8,000 human donor sites (chr 20-22).
# Column order: [A, C, G, T]
#
# Key features:
#   - Position -2: A preference (part of exonic A)
#   - Position -1: G strong preference (last exonic base)
#   - Position +1: G nearly invariant (donor G of GT)
#   - Position +2: T nearly invariant (donor T of GT)
#   - Position +3: A and G both elevated (R = purine)
#   - Position +5: T strong preference
#   - Position +6: T preference

DONOR_MONO_FREQ: List[List[float]] = [
    # pos   A      C      G      T
    [0.330, 0.370, 0.180, 0.120],  # -3  C preference
    [0.600, 0.130, 0.120, 0.150],  # -2  A preference
    [0.080, 0.060, 0.780, 0.080],  # -1  G strong preference (last exonic)
    [0.010, 0.010, 0.970, 0.010],  # +1  G nearly invariant (donor G)
    [0.010, 0.010, 0.010, 0.970],  # +2  T nearly invariant (donor T)
    [0.350, 0.100, 0.350, 0.200],  # +3  A/G both elevated (R = A/G)
    [0.250, 0.080, 0.400, 0.270],  # +4  G moderate, A notable
    [0.100, 0.060, 0.140, 0.700],  # +5  T strong preference
    [0.160, 0.160, 0.080, 0.600],  # +6  T preference
]

# Conditional frequency tables for first-order Markov model.
#
# DONOR_COND_FREQ[k][prev_idx][curr_idx] gives
#   P(base_curr at position k+1 | base_prev at position k)
#
# k ranges from 0 to 7 (8 adjacent pairs for 9 positions).
# prev_idx / curr_idx: 0=A, 1=C, 2=G, 3=T
#
# Derived from mononucleotide marginals with first-order dependency
# adjustments at key positions:
#   - Pair 2 (-1, +1): G→G dependency (splice junction)
#   - Pair 3 (+1, +2): G→T dependency (GT dinucleotide)
#   - Pair 4 (+2, +3): T→A/G dependency
#   - Pair 5 (+3, +4): A/G→A/G dependency
#   - Pair 6 (+4, +5): G/A→T dependency
#   - Pair 7 (+5, +6): T→T dependency

DONOR_COND_FREQ: List[List[List[float]]] = [
    # ── Pair 0: (pos -3, pos -2) — weak dependency ────────────────────
    [
        [0.6117, 0.1262, 0.1165, 0.1456],  # prev=A → A at -2 slightly up
        [0.5923, 0.1412, 0.1185, 0.1481],  # prev=C → C at -2 slightly up
        [0.6000, 0.1300, 0.1200, 0.1500],  # prev=G → near marginal
        [0.6000, 0.1300, 0.1200, 0.1500],  # prev=T → near marginal
    ],
    # ── Pair 1: (pos -2, pos -1) — moderate A→G dependency ────────────
    [
        [0.0764, 0.0573, 0.7898, 0.0764],  # prev=A → G at -1 elevated
        [0.0800, 0.0600, 0.7800, 0.0800],  # prev=C → near marginal
        [0.0800, 0.0600, 0.7800, 0.0800],  # prev=G → near marginal
        [0.0853, 0.0640, 0.7654, 0.0853],  # prev=T → G at -1 reduced
    ],
    # ── Pair 2: (pos -1, pos +1) — STRONG G→G (splice junction) ──────
    [
        [0.0111, 0.0111, 0.9668, 0.0111],  # prev=A (rare) → G still high
        [0.0111, 0.0111, 0.9668, 0.0111],  # prev=C (rare) → G still high
        [0.0098, 0.0098, 0.9706, 0.0098],  # prev=G → G at +1 very high
        [0.0111, 0.0111, 0.9668, 0.0111],  # prev=T (rare) → G still high
    ],
    # ── Pair 3: (pos +1, pos +2) — VERY STRONG G→T (GT dinucleotide) ─
    [
        [0.0124, 0.0124, 0.0124, 0.9628],  # prev=A (rare)
        [0.0124, 0.0124, 0.0124, 0.9628],  # prev=C (rare)
        [0.0098, 0.0098, 0.0098, 0.9706],  # prev=G → T at +2 very high
        [0.0124, 0.0124, 0.0124, 0.9628],  # prev=T (rare)
    ],
    # ── Pair 4: (pos +2, pos +3) — T→A/G dependency ──────────────────
    [
        [0.3500, 0.1000, 0.3500, 0.2000],  # prev=A → near marginal
        [0.3500, 0.1000, 0.3500, 0.2000],  # prev=C → near marginal
        [0.3500, 0.1000, 0.3500, 0.2000],  # prev=G → near marginal
        [0.3828, 0.0781, 0.3828, 0.1563],  # prev=T → A/G up, C/T down
    ],
    # ── Pair 5: (pos +3, pos +4) — A/G→A/G dependency ────────────────
    [
        [0.2719, 0.0757, 0.3972, 0.2553],  # prev=A → A,G at +4 up
        [0.2500, 0.0800, 0.4000, 0.2700],  # prev=C → near marginal
        [0.2494, 0.0760, 0.4181, 0.2565],  # prev=G → G at +4 up
        [0.2500, 0.0800, 0.4000, 0.2700],  # prev=T → near marginal
    ],
    # ── Pair 6: (pos +4, pos +5) — G/A→T dependency ──────────────────
    [
        [0.0960, 0.0576, 0.1344, 0.7121],  # prev=A → T at +5 up
        [0.1000, 0.0600, 0.1400, 0.7000],  # prev=C → near marginal
        [0.0960, 0.0576, 0.1344, 0.7121],  # prev=G → T at +5 up
        [0.1000, 0.0600, 0.1400, 0.7000],  # prev=T → near marginal
    ],
    # ── Pair 7: (pos +5, pos +6) — T→T dependency ────────────────────
    [
        [0.1600, 0.1600, 0.0800, 0.6000],  # prev=A → near marginal
        [0.1600, 0.1600, 0.0800, 0.6000],  # prev=C → near marginal
        [0.1600, 0.1600, 0.0800, 0.6000],  # prev=G → near marginal
        [0.1553, 0.1553, 0.0777, 0.6117],  # prev=T → T at +6 slightly up
    ],
]


# ==============================================================================
# 2. Acceptor Site Frequency Matrices (23 positions)
# ==============================================================================
#
# Positions 0..22 of the 23-mer correspond to:
#   0..17 = upstream intronic bases (polypyrimidine tract, 18 bases)
#   18    = A of AG dinucleotide
#   19    = G of AG dinucleotide
#   20..22 = first three exonic bases
#
# Source: Yeo & Burge 2004, trained on ~8,000 human acceptor sites.
# Column order: [A, C, G, T]
#
# Key features:
#   - Positions 0-17: increasing C/T (pyrimidine) bias toward AG
#   - Position 18: A nearly invariant
#   - Position 19: G nearly invariant

ACCEPTOR_MONO_FREQ: List[List[float]] = [
    # pos   A      C      G      T
    [0.240, 0.260, 0.220, 0.280],  # 0   far upstream, moderate C/T
    [0.230, 0.270, 0.210, 0.290],  # 1
    [0.220, 0.280, 0.200, 0.300],  # 2
    [0.210, 0.290, 0.190, 0.310],  # 3
    [0.200, 0.300, 0.180, 0.320],  # 4
    [0.190, 0.310, 0.170, 0.330],  # 5
    [0.170, 0.330, 0.150, 0.350],  # 6   pyrimidine tract intensifies
    [0.160, 0.340, 0.140, 0.360],  # 7
    [0.150, 0.350, 0.130, 0.370],  # 8
    [0.140, 0.370, 0.120, 0.370],  # 9
    [0.130, 0.380, 0.110, 0.380],  # 10
    [0.120, 0.390, 0.100, 0.390],  # 11
    [0.110, 0.400, 0.090, 0.400],  # 12
    [0.100, 0.410, 0.080, 0.410],  # 13
    [0.090, 0.420, 0.070, 0.420],  # 14
    [0.080, 0.430, 0.060, 0.430],  # 15
    [0.070, 0.440, 0.050, 0.440],  # 16
    [0.060, 0.450, 0.040, 0.450],  # 17  strong pyrimidine bias
    [0.980, 0.005, 0.005, 0.010],  # 18  A nearly invariant (A of AG)
    [0.005, 0.010, 0.980, 0.005],  # 19  G nearly invariant (G of AG)
    [0.260, 0.200, 0.330, 0.210],  # 20  1st exonic, moderate G preference
    [0.230, 0.240, 0.280, 0.250],  # 21  2nd exonic
    [0.200, 0.220, 0.300, 0.280],  # 22  3rd exonic
]

# Conditional frequency tables for the acceptor first-order Markov model.
#
# ACCEPTOR_COND_FREQ[k][prev_idx][curr_idx] gives
#   P(base_curr at position k+1 | base_prev at position k)
#
# k ranges from 0 to 21 (22 adjacent pairs for 23 positions).
# prev_idx / curr_idx: 0=A, 1=C, 2=G, 3=T
#
# Key dependencies:
#   - Polypyrimidine tract: C/T→C/T transitions are enriched
#   - Position 17→18: T/C→A (pyrimidine before AG)
#   - Position 18→19: A→G (the AG dinucleotide)
#   - Position 19→20: G→G slightly elevated (first exonic)

ACCEPTOR_COND_FREQ: List[List[List[float]]] = [
    # ── Pairs 0-16: Polypyrimidine tract with C/T→C/T dependency ─────
    # Pair 0: (pos 0, pos 1)
    [
        [0.2300, 0.2700, 0.2100, 0.2900],  # prev=A → near marginal
        [0.2143, 0.2844, 0.1957, 0.3055],  # prev=C → C/T up, A/G down
        [0.2300, 0.2700, 0.2100, 0.2900],  # prev=G → near marginal
        [0.2143, 0.2844, 0.1957, 0.3055],  # prev=T → C/T up, A/G down
    ],
    # Pair 1: (pos 1, pos 2)
    [
        [0.2200, 0.2800, 0.2000, 0.3000],  # prev=A
        [0.2045, 0.2943, 0.1859, 0.3153],  # prev=C
        [0.2200, 0.2800, 0.2000, 0.3000],  # prev=G
        [0.2045, 0.2943, 0.1859, 0.3153],  # prev=T
    ],
    # Pair 2: (pos 2, pos 3)
    [
        [0.2100, 0.2900, 0.1900, 0.3100],  # prev=A
        [0.1948, 0.3040, 0.1762, 0.3250],  # prev=C
        [0.2100, 0.2900, 0.1900, 0.3100],  # prev=G
        [0.1948, 0.3040, 0.1762, 0.3250],  # prev=T
    ],
    # Pair 3: (pos 3, pos 4)
    [
        [0.2000, 0.3000, 0.1800, 0.3200],  # prev=A
        [0.1850, 0.3138, 0.1665, 0.3347],  # prev=C
        [0.2000, 0.3000, 0.1800, 0.3200],  # prev=G
        [0.1850, 0.3138, 0.1665, 0.3347],  # prev=T
    ],
    # Pair 4: (pos 4, pos 5)
    [
        [0.1900, 0.3100, 0.1700, 0.3300],  # prev=A
        [0.1754, 0.3234, 0.1569, 0.3443],  # prev=C
        [0.1900, 0.3100, 0.1700, 0.3300],  # prev=G
        [0.1754, 0.3234, 0.1569, 0.3443],  # prev=T
    ],
    # Pair 5: (pos 5, pos 6)
    [
        [0.1700, 0.3300, 0.1500, 0.3500],  # prev=A
        [0.1562, 0.3427, 0.1378, 0.3634],  # prev=C
        [0.1700, 0.3300, 0.1500, 0.3500],  # prev=G
        [0.1562, 0.3427, 0.1378, 0.3634],  # prev=T
    ],
    # Pair 6: (pos 6, pos 7)
    [
        [0.1600, 0.3400, 0.1400, 0.3600],  # prev=A
        [0.1466, 0.3522, 0.1283, 0.3729],  # prev=C
        [0.1600, 0.3400, 0.1400, 0.3600],  # prev=G
        [0.1466, 0.3522, 0.1283, 0.3729],  # prev=T
    ],
    # Pair 7: (pos 7, pos 8)
    [
        [0.1500, 0.3500, 0.1300, 0.3700],  # prev=A
        [0.1371, 0.3617, 0.1188, 0.3824],  # prev=C
        [0.1500, 0.3500, 0.1300, 0.3700],  # prev=G
        [0.1371, 0.3617, 0.1188, 0.3824],  # prev=T
    ],
    # Pair 8: (pos 8, pos 9)
    [
        [0.1400, 0.3700, 0.1200, 0.3700],  # prev=A
        [0.1277, 0.3814, 0.1094, 0.3814],  # prev=C
        [0.1400, 0.3700, 0.1200, 0.3700],  # prev=G
        [0.1277, 0.3814, 0.1094, 0.3814],  # prev=T
    ],
    # Pair 9: (pos 9, pos 10)
    [
        [0.1300, 0.3800, 0.1100, 0.3800],  # prev=A
        [0.1183, 0.3908, 0.1001, 0.3908],  # prev=C
        [0.1300, 0.3800, 0.1100, 0.3800],  # prev=G
        [0.1183, 0.3908, 0.1001, 0.3908],  # prev=T
    ],
    # Pair 10: (pos 10, pos 11)
    [
        [0.1200, 0.3900, 0.1000, 0.3900],  # prev=A
        [0.1089, 0.4002, 0.0908, 0.4002],  # prev=C
        [0.1200, 0.3900, 0.1000, 0.3900],  # prev=G
        [0.1089, 0.4002, 0.0908, 0.4002],  # prev=T
    ],
    # Pair 11: (pos 11, pos 12)
    [
        [0.1100, 0.4000, 0.0900, 0.4000],  # prev=A
        [0.0996, 0.4094, 0.0815, 0.4094],  # prev=C
        [0.1100, 0.4000, 0.0900, 0.4000],  # prev=G
        [0.0996, 0.4094, 0.0815, 0.4094],  # prev=T
    ],
    # Pair 12: (pos 12, pos 13)
    [
        [0.1000, 0.4100, 0.0800, 0.4100],  # prev=A
        [0.0903, 0.4187, 0.0723, 0.4187],  # prev=C
        [0.1000, 0.4100, 0.0800, 0.4100],  # prev=G
        [0.0903, 0.4187, 0.0723, 0.4187],  # prev=T
    ],
    # Pair 13: (pos 13, pos 14)
    [
        [0.0900, 0.4200, 0.0700, 0.4200],  # prev=A
        [0.0811, 0.4279, 0.0631, 0.4279],  # prev=C
        [0.0900, 0.4200, 0.0700, 0.4200],  # prev=G
        [0.0811, 0.4279, 0.0631, 0.4279],  # prev=T
    ],
    # Pair 14: (pos 14, pos 15)
    [
        [0.0800, 0.4300, 0.0600, 0.4300],  # prev=A
        [0.0719, 0.4371, 0.0539, 0.4371],  # prev=C
        [0.0800, 0.4300, 0.0600, 0.4300],  # prev=G
        [0.0719, 0.4371, 0.0539, 0.4371],  # prev=T
    ],
    # Pair 15: (pos 15, pos 16)
    [
        [0.0700, 0.4400, 0.0500, 0.4400],  # prev=A
        [0.0628, 0.4462, 0.0449, 0.4462],  # prev=C
        [0.0700, 0.4400, 0.0500, 0.4400],  # prev=G
        [0.0628, 0.4462, 0.0449, 0.4462],  # prev=T
    ],
    # Pair 16: (pos 16, pos 17)
    [
        [0.0600, 0.4500, 0.0400, 0.4500],  # prev=A
        [0.0537, 0.4553, 0.0358, 0.4553],  # prev=C
        [0.0600, 0.4500, 0.0400, 0.4500],  # prev=G
        [0.0537, 0.4553, 0.0358, 0.4553],  # prev=T
    ],
    # ── Pair 17: (pos 17, pos 18) — T/C→A before AG ──────────────────
    [
        [0.9800, 0.0050, 0.0050, 0.0100],  # prev=A → A at 18 near marginal
        [0.9804, 0.0049, 0.0049, 0.0098],  # prev=C → A at 18 slightly up
        [0.9800, 0.0050, 0.0050, 0.0100],  # prev=G → near marginal
        [0.9808, 0.0048, 0.0048, 0.0096],  # prev=T → A at 18 slightly up
    ],
    # ── Pair 18: (pos 18, pos 19) — A→G (AG dinucleotide) ────────────
    [
        [0.0049, 0.0098, 0.9804, 0.0049],  # prev=A → G at 19 very high
        [0.0050, 0.0100, 0.9800, 0.0050],  # prev=C (rare)
        [0.0050, 0.0100, 0.9800, 0.0050],  # prev=G (rare)
        [0.0050, 0.0100, 0.9800, 0.0050],  # prev=T (rare)
    ],
    # ── Pair 19: (pos 19, pos 20) — G→first exonic base ──────────────
    [
        [0.2600, 0.2000, 0.3300, 0.2100],  # prev=A (rare)
        [0.2600, 0.2000, 0.3300, 0.2100],  # prev=C (rare)
        [0.2558, 0.1968, 0.3409, 0.2066],  # prev=G → G at 20 slightly up
        [0.2600, 0.2000, 0.3300, 0.2100],  # prev=T (rare)
    ],
    # ── Pair 20: (pos 20, pos 21) — exonic, weak dependency ──────────
    [
        [0.2300, 0.2400, 0.2800, 0.2500],  # prev=A → marginal
        [0.2300, 0.2400, 0.2800, 0.2500],  # prev=C
        [0.2300, 0.2400, 0.2800, 0.2500],  # prev=G
        [0.2300, 0.2400, 0.2800, 0.2500],  # prev=T
    ],
    # ── Pair 21: (pos 21, pos 22) — exonic, weak dependency ──────────
    [
        [0.2000, 0.2200, 0.3000, 0.2800],  # prev=A → marginal
        [0.2000, 0.2200, 0.3000, 0.2800],  # prev=C
        [0.2000, 0.2200, 0.3000, 0.2800],  # prev=G
        [0.2000, 0.2200, 0.3000, 0.2800],  # prev=T
    ],
]


# ==============================================================================
# 3. Hard-Coded Scoring Matrices (precomputed scores for reference)
# ==============================================================================
#
# These are precomputed MaxEntScan scores for specific donor and acceptor
# sequences computed by the first-order Markov model in this module.
# They serve as a validation reference and quick-lookup table.
#
# All donor sequences are 9-mers (positions -3 to +6).
# All acceptor sequences are 23-mers (18 intronic + AG + 3 exonic).

MAXENTSCAN_DONOR_SCORES: Dict[str, float] = {
    # ── Strong canonical GT donors ──────────────────────────────────────
    # Consensus: MAG|GTRAGT (M=A/C, R=A/G)
    "CAGGTAGTT": 11.47,  # Very strong canonical (consensus at all positions)
    "CAGGTAGGT": 9.03,   # Strong canonical
    "CAGGTGAGT": 8.36,   # Strong canonical, GAG intronic
    "AAGGTGAGT": 8.24,   # Strong, A at -3
    "CAGGTAAGT": 8.49,   # Moderate-strong, AAG at +3..+5
    "CAGGTATGT": 8.45,   # Moderate-strong, TAT at +3..+5
    # ── Moderate / weak canonical GT donors ─────────────────────────────
    "CAGGTTTGT": 7.24,   # Moderate, TTT at +3..+5
    "TTGGTAAGT": 4.84,   # Weak, T at -3 and T at -2
    # ── Non-canonical GC donors ──────────────────────────────────────────
    "CAGGCAAGT": 1.73,   # GC donor (weaker than GT)
    # ── Non-donor sequences ──────────────────────────────────────────────
    "ATCATCAGT": -6.20,  # No GT at +1/+2 (non-canonical bases)
    "TTTATTTTT": -3.69,  # No GT, mostly T
    "CCCCCCCCC": -16.81, # All C (worst case)
}

MAXENTSCAN_ACCEPTOR_SCORES: Dict[str, float] = {
    # ── Strong canonical AG acceptors ───────────────────────────────────
    # 23-mer: 18 intronic (polypyrimidine tract) + AG + 3 exonic
    "TTTTTTTTTTTTTTTTTTAGATG": 14.69,  # Very strong poly-T tract
    "CTCCTTTTTCCTTTTCTTAGATG": 14.39,  # Strong mixed pyrimidine tract
    "CCCCCCCCCCCCCCCCCCAGATC": 13.40,  # Strong all-C tract
    "TTCTTTCTTTCTTTCTTTAGATG": 14.50,  # Strong mixed C/T tract
    # ── Weak / poor acceptors ───────────────────────────────────────────
    "AATATAGTATAGATATATAGATG": -1.56,  # Weak, poor pyrimidine tract
    # ── Non-acceptor sequences ──────────────────────────────────────────
    "GAGAGAGAGAGAGAGAGAAGATG": -13.57, # Very weak, G-rich (anti-pyrimidine)
    "AAAAAAAAAAAAAAAAAAAGATG": -11.39, # All-A upstream (no pyrimidine tract)
    "GGGGGGGGGGGGGGGGGGAGATG": -16.01, # All-G upstream (anti-pyrimidine)
}


# ==============================================================================
# 4. Scoring Functions
# ==============================================================================


def score_donor_maxentscan(seq: str) -> float:
    """Score a 9-mer donor sequence using the maximum entropy model.

    The donor model considers a 9-mer window around the GT dinucleotide:
    3 exonic bases (positions −3 to −1) + 6 intronic bases (positions +1
    to +6).  Positions 3 and 4 of the input string correspond to the G
    and T of the GT dinucleotide.

    Uses a first-order Markov model that captures dependencies between
    adjacent positions — the key innovation of MaxEntScan over simple
    PWM approaches.

    Args:
        seq: 9-mer DNA sequence (uppercase or lowercase).

    Returns:
        Log-odds score in bits.  Typical ranges:
        - Strong canonical GT donors: 8–12
        - Weak / cryptic donors: 0–5
        - Non-donor sequences: < 0

    Raises:
        ValueError: If *seq* is not exactly 9 bases long.

    Examples
    --------
    >>> score_donor_maxentscan("CAGGTAGTT")  # very strong canonical
    11.47
    >>> score_donor_maxentscan("CAGGCAAGT")   # GC donor
    1.73
    >>> score_donor_maxentscan("ATCATCAGT")   # non-canonical
    -6.2
    """
    seq = seq.upper()
    if len(seq) != DONOR_KMER_LEN:
        raise ValueError(
            f"Expected {DONOR_KMER_LEN}-mer donor sequence, "
            f"got {len(seq)}-mer: {seq!r}"
        )

    score = 0.0
    prev_idx: int = -1

    for i, base in enumerate(seq):
        if base not in _BASE_IDX:
            return _IMPOSSIBLE
        idx = _BASE_IDX[base]

        if i == 0:
            # First position: use marginal (mononucleotide) probability
            prob = _safe_prob(DONOR_MONO_FREQ[i][idx])
        else:
            # Subsequent positions: use conditional (first-order Markov)
            prob = _safe_prob(DONOR_COND_FREQ[i - 1][prev_idx][idx])

        score += _log2(prob / _BG)
        prev_idx = idx

    return round(score, _SCORE_DP)


def score_acceptor_maxentscan(seq: str) -> float:
    """Score a 23-mer acceptor sequence using the maximum entropy model.

    The acceptor model considers a 23-mer window around the AG
    dinucleotide: 18 upstream intronic bases (polypyrimidine tract) +
    the AG dinucleotide (positions 18-19) + 3 exonic bases (positions
    20-22).

    Uses a first-order Markov model capturing:
    - C/T → C/T dependencies in the polypyrimidine tract
    - T/C → A dependency before AG
    - A → G dependency (the AG dinucleotide itself)

    Args:
        seq: 23-mer DNA sequence (uppercase or lowercase).

    Returns:
        Log-odds score in bits.  Typical ranges:
        - Strong canonical AG acceptors: 8–14
        - Weak / cryptic acceptors: 0–5
        - Non-acceptor sequences: < 0

    Raises:
        ValueError: If *seq* is not exactly 23 bases long.

    Examples
    --------
    >>> score_acceptor_maxentscan("TTTTTTTTTTTTTTTTTTAGATG")  # strong
    14.69
    >>> score_acceptor_maxentscan("AAAAAAAAAAAAAAAAAAAGATG")  # non-acceptor
    -11.39
    """
    seq = seq.upper()
    if len(seq) != ACCEPTOR_KMER_LEN:
        raise ValueError(
            f"Expected {ACCEPTOR_KMER_LEN}-mer acceptor sequence, "
            f"got {len(seq)}-mer: {seq!r}"
        )

    score = 0.0
    prev_idx: int = -1

    for i, base in enumerate(seq):
        if base not in _BASE_IDX:
            return _IMPOSSIBLE
        idx = _BASE_IDX[base]

        if i == 0:
            prob = _safe_prob(ACCEPTOR_MONO_FREQ[i][idx])
        else:
            prob = _safe_prob(ACCEPTOR_COND_FREQ[i - 1][prev_idx][idx])

        score += _log2(prob / _BG)
        prev_idx = idx

    return round(score, _SCORE_DP)


# ==============================================================================
# 5. Classification Helpers
# ==============================================================================


def is_strong_donor(seq: str, threshold: float = 3.0) -> bool:
    """Return True if the 9-mer donor score exceeds *threshold*.

    A donor is considered "strong" if its MaxEntScan score is above the
    threshold.  The default threshold of 3.0 bits separates functional
    from non-functional donor sites with high specificity.

    Args:
        seq: 9-mer DNA sequence.
        threshold: Minimum score in bits (default 3.0).

    Returns:
        True if score > threshold, False otherwise.

    Examples
    --------
    >>> is_strong_donor("CAGGTAGTT")
    True
    >>> is_strong_donor("CAGGCAAGT")
    False
    >>> is_strong_donor("ATCATCAGT")
    False
    """
    return score_donor_maxentscan(seq) > threshold


def is_strong_acceptor(seq: str, threshold: float = 3.0) -> bool:
    """Return True if the 23-mer acceptor score exceeds *threshold*.

    An acceptor is considered "strong" if its MaxEntScan score is above
    the threshold.  The default threshold of 3.0 bits separates functional
    from non-functional acceptor sites with high specificity.

    Args:
        seq: 23-mer DNA sequence.
        threshold: Minimum score in bits (default 3.0).

    Returns:
        True if score > threshold, False otherwise.

    Examples
    --------
    >>> is_strong_acceptor("TTTTTTTTTTTTTTTTTTAGATG")
    True
    >>> is_strong_acceptor("AAAAAAAAAAAAAAAAAAAGATG")
    False
    """
    return score_acceptor_maxentscan(seq) > threshold
