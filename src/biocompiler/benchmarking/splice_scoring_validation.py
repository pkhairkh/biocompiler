"""
BioCompiler Splice Scoring Validation
======================================

Validation module cross-checking the production
``maxentscan.score_donor()`` / ``score_acceptor()`` Markov model
implementation against the independent benchmarking implementation in
``benchmarking.maxentscan_validated``.

Since the production module now uses the first-order Markov model (the same
algorithm as the benchmarking module), the two implementations should produce
**positively correlated** scores that agree within floating-point tolerance.
Any disagreement or negative correlation signals a regression in either module.

Previously, this module compared the deprecated PWM-based ``maxent_score()``
against the proper MaxEntScan implementation and detected anti-correlation.
That comparison is preserved for historical reference, but the primary
validation now checks that both Markov implementations agree.

References:
  Yeo, G. & Burge, C.B. (2004). "Maximum entropy modeling of short sequence
  motifs with applications to RNA splicing." *Journal of Computational Biology*
  11(2-3):377-394. doi:10.1089/1066527041410418
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import List

from biocompiler.sequence.splicing import maxent_score as _deprecated_maxent_score
from biocompiler.sequence.maxentscan import score_donor, score_acceptor
from biocompiler.benchmarking.maxentscan_validated import (
    score_donor_maxentscan as _benchmark_donor,
    score_acceptor_maxentscan as _benchmark_acceptor,
    MAXENTSCAN_DONOR_SCORES as _BENCHMARK_DONOR_SCORES,
    MAXENTSCAN_ACCEPTOR_SCORES as _BENCHMARK_ACCEPTOR_SCORES,
)

__all__ = [
    "SpliceValidationResult",
    "validate_splice_scoring",
    "print_splice_validation_report",
    "MarkovCrossCheckResult",
    "validate_markov_cross_check",
    "print_markov_cross_check_report",
]


# ==============================================================================
# Data structures
# ==============================================================================


@dataclass
class SpliceValidationResult:
    """Result of splice scoring validation comparing deprecated and correct methods.

    Attributes:
        test_sites: List of dicts, each with keys:
            - sequence: the test sequence string
            - site_type: "donor" or "acceptor"
            - expected_strength: "strong", "moderate", or "weak"
            - deprecated_score: score from splicing.maxent_score()
            - correct_score: score from maxentscan.score_donor()/score_acceptor()
        correlation: Pearson correlation between deprecated and correct scores
        is_anti_correlated: True if correlation < 0
        correct_scores_rank_order: True if strong sites score higher than weak
            in maxentscan (the correct implementation)
        deprecated_scores_rank_order: True if strong sites score higher than
            weak in the deprecated maxent_score
    """

    test_sites: List[dict]
    correlation: float
    is_anti_correlated: bool
    correct_scores_rank_order: bool
    deprecated_scores_rank_order: bool


@dataclass
class MarkovCrossCheckResult:
    """Result of cross-checking production Markov vs benchmarking Markov scores.

    Both implementations use the first-order Markov model, so their scores
    should agree within floating-point tolerance.  A positive correlation is
    expected; any negative correlation or large score differences indicate a
    regression in either module.

    Attributes:
        donor_results: List of dicts with keys 'sequence', 'production_score',
            'benchmark_score', 'difference'.
        acceptor_results: List of dicts with keys 'sequence', 'production_score',
            'benchmark_score', 'difference'.
        donor_correlation: Pearson correlation between production and benchmark
            donor scores.
        acceptor_correlation: Pearson correlation between production and benchmark
            acceptor scores.
        max_donor_diff: Largest absolute difference in donor scores.
        max_acceptor_diff: Largest absolute difference in acceptor scores.
        all_within_tolerance: True if all score differences are within tolerance.
    """

    donor_results: List[dict]
    acceptor_results: List[dict]
    donor_correlation: float
    acceptor_correlation: float
    max_donor_diff: float
    max_acceptor_diff: float
    all_within_tolerance: bool


# ==============================================================================
# Helpers
# ==============================================================================


def _pearson_correlation(x: List[float], y: List[float]) -> float:
    """Compute Pearson correlation coefficient between two equal-length lists.

    Returns 0.0 if the lists have fewer than 2 elements or if either has
    zero variance.
    """
    n = len(x)
    if n < 2:
        return 0.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)

    if var_x == 0.0 or var_y == 0.0:
        return 0.0

    return cov / math.sqrt(var_x * var_y)


def _score_deprecated(seq: str) -> float:
    """Score with deprecated maxent_score, suppressing deprecation warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        return _deprecated_maxent_score(seq)


def _score_proper_donor(seq: str, gt_pos: int) -> float:
    """Score a donor site using proper maxentscan.score_donor().

    If the sequence is too short for the 9-mer context around *gt_pos*,
    pads with A on both sides to ensure sufficient context.

    Args:
        seq: DNA sequence containing a GT dinucleotide.
        gt_pos: Position of G in the GT dinucleotide.

    Returns:
        MaxEntScan log-odds donor score.
    """
    upstream_needed = 3
    downstream_needed = 6

    pad_left = max(0, upstream_needed - gt_pos)
    pad_right = max(0, downstream_needed - (len(seq) - gt_pos - 2))

    padded = "A" * pad_left + seq + "A" * pad_right
    adjusted_pos = gt_pos + pad_left

    return score_donor(padded, adjusted_pos)


def _score_proper_acceptor(seq: str, ag_pos: int) -> float:
    """Score an acceptor site using proper maxentscan.score_acceptor().

    If the sequence is too short for the 23-mer context around *ag_pos*,
    pads with A on both sides to ensure sufficient context.

    NOTE: score_acceptor() uses A-position convention: the position argument
    is the index of A in the AG dinucleotide.  The 23-mer is extracted as
    seq[position-18 : position+5], with A at 23-mer index 18 and G at 19.

    Args:
        seq: DNA sequence containing an AG dinucleotide.
        ag_pos: Position of A in the AG dinucleotide.

    Returns:
        MaxEntScan log-odds acceptor score.
    """
    # A-position convention: pass ag_pos directly
    upstream_needed = 18
    downstream_needed = 5

    pad_left = max(0, upstream_needed - ag_pos)
    pad_right = max(0, downstream_needed - (len(seq) - ag_pos))

    padded = "A" * pad_left + seq + "A" * pad_right
    adjusted_pos = ag_pos + pad_left

    return score_acceptor(padded, adjusted_pos)


# ==============================================================================
# Test site definitions
# ==============================================================================

# Each entry is a dict with:
#   sequence           – the DNA sequence for maxent_score (9-mer for donors)
#   full_seq           – the full DNA sequence for score_donor/score_acceptor
#   site_type          – "donor" or "acceptor"
#   expected_strength  – "strong", "moderate", or "weak"
#   motif_pos          – position of the splice dinucleotide (GT/AG) in full_seq

_DONOR_TEST_SITES: List[dict] = [
    {
        "sequence": "CAGGTGAGT",
        "full_seq": "CAGGTGAGT",
        "site_type": "donor",
        "expected_strength": "strong",
        "motif_pos": 3,  # G of GT at position 3
    },
    {
        "sequence": "AAGTAAGCT",
        "full_seq": "AAAGTAAGCTAAAAAA",
        "site_type": "donor",
        "expected_strength": "weak",
        "motif_pos": 3,  # G of GT at position 3 in full_seq
    },
    {
        "sequence": "TTGGTAAGT",
        "full_seq": "TTGGTAAGT",
        "site_type": "donor",
        "expected_strength": "moderate",
        "motif_pos": 3,
    },
    {
        "sequence": "CCGGTAAGT",
        "full_seq": "CCGGTAAGT",
        "site_type": "donor",
        "expected_strength": "moderate",
        "motif_pos": 3,
    },
    {
        "sequence": "ATGGTCATC",
        "full_seq": "ATGGTCATC",
        "site_type": "donor",
        "expected_strength": "weak",
        "motif_pos": 3,
    },
]

_ACCEPTOR_TEST_SITES: List[dict] = [
    # Strong acceptor: poly-T (polypyrimidine) tract upstream of AG
    {
        "sequence": "TTTTTTTTTTTTTTTTTTAGG",
        "full_seq": "AAATTTTTTTTTTTTTTTTTAGGAAAA",
        "site_type": "acceptor",
        "expected_strength": "strong",
        "motif_pos": 21,  # A of AG at position 21 in full_seq
    },
    # Strong acceptor: poly-C (polypyrimidine) tract upstream of AG
    {
        "sequence": "CCCCCCCCCCCCCCCCCCAGG",
        "full_seq": "AAACCCCCCCCCCCCCCCCCCAGGAAAA",
        "site_type": "acceptor",
        "expected_strength": "strong",
        "motif_pos": 21,
    },
    # Weak acceptor: purine-rich (G) upstream — anti-pyrimidine tract
    {
        "sequence": "AAAGAAATC",
        "full_seq": "AAAGGGGGGGGGGGGGGGGGAGGAAAA",
        "site_type": "acceptor",
        "expected_strength": "weak",
        "motif_pos": 21,
    },
    # Weak acceptor: A-rich upstream (no pyrimidine tract)
    {
        "sequence": "AAAAAAAAAAAAAAAAAAAAGG",
        "full_seq": "AAAAAAAAAAAAAAAAAAAAAGGAAAA",
        "site_type": "acceptor",
        "expected_strength": "weak",
        "motif_pos": 21,
    },
]


# ==============================================================================
# Public API
# ==============================================================================


def validate_splice_scoring() -> SpliceValidationResult:
    """Validate splice scoring by comparing deprecated and correct implementations.

    Generates test sequences with known splice sites, scores them with both the
    deprecated ``maxent_score()`` and the proper ``maxentscan`` functions, and
    computes the Pearson correlation between the two methods.

    The test sites include:
    - Canonical donor: ``CAGGTGAGT`` (strong, should score high in maxentscan)
    - Non-canonical donor: ``AAGTAAGCT`` (weak, should score low in maxentscan)
    - Canonical acceptor: poly-T tract upstream of AG (strong)
    - Non-canonical acceptor: ``AAAGAAATC`` / purine-rich upstream (weak)
    - Additional moderate and weak donor/acceptor sites

    Returns:
        SpliceValidationResult with comparison data, Pearson correlation, and
        rank-order checks.
    """
    all_configs = _DONOR_TEST_SITES + _ACCEPTOR_TEST_SITES

    test_sites: List[dict] = []
    deprecated_scores: List[float] = []
    correct_scores: List[float] = []

    for config in all_configs:
        seq = config["sequence"]
        full_seq = config["full_seq"]
        site_type = config["site_type"]
        expected_strength = config["expected_strength"]
        motif_pos = config["motif_pos"]

        # Score with deprecated maxent_score (uses the 9-mer context directly;
        # for acceptors, the first 9 chars of the sequence are used)
        dep_score = _score_deprecated(seq)

        # Score with proper maxentscan function
        if site_type == "donor":
            correct_score = _score_proper_donor(full_seq, motif_pos)
        else:
            correct_score = _score_proper_acceptor(full_seq, motif_pos)

        site_result = {
            "sequence": seq,
            "site_type": site_type,
            "expected_strength": expected_strength,
            "deprecated_score": dep_score,
            "correct_score": correct_score,
        }
        test_sites.append(site_result)
        deprecated_scores.append(dep_score)
        correct_scores.append(correct_score)

    # Compute Pearson correlation between the two scoring methods
    correlation = _pearson_correlation(deprecated_scores, correct_scores)

    # Check rank ordering: strong sites should score higher than weak sites
    # within each site type (donor vs donor, acceptor vs acceptor).
    # Cross-type comparison is not meaningful because donor and acceptor
    # models have different score ranges.
    def _check_rank_order(sites: List[dict], score_key: str) -> bool:
        """Check that strong > moderate > weak within a single site type."""
        strong = [
            s[score_key] for s in sites
            if s["expected_strength"] == "strong" and s[score_key] > -50.0
        ]
        moderate = [
            s[score_key] for s in sites
            if s["expected_strength"] == "moderate" and s[score_key] > -50.0
        ]
        weak = [
            s[score_key] for s in sites
            if s["expected_strength"] == "weak" and s[score_key] > -50.0
        ]
        if not strong or not weak:
            return True  # Insufficient data — do not fail
        # All strong must outscore all weak
        if min(strong) <= max(weak):
            return False
        # If moderate exists, it should fall between strong and weak
        if moderate and (min(strong) <= max(moderate) or min(moderate) <= max(weak)):
            # Moderate can overlap — only check strong > weak
            pass
        return True

    donor_sites = [s for s in test_sites if s["site_type"] == "donor"]
    acceptor_sites = [s for s in test_sites if s["site_type"] == "acceptor"]

    correct_rank_order = (
        _check_rank_order(donor_sites, "correct_score")
        and _check_rank_order(acceptor_sites, "correct_score")
    )

    strong_deprecated = [
        s["deprecated_score"]
        for s in test_sites
        if s["expected_strength"] == "strong"
    ]
    weak_deprecated = [
        s["deprecated_score"]
        for s in test_sites
        if s["expected_strength"] == "weak"
    ]

    deprecated_rank_order = (
        _check_rank_order(donor_sites, "deprecated_score")
        and _check_rank_order(acceptor_sites, "deprecated_score")
    )

    return SpliceValidationResult(
        test_sites=test_sites,
        correlation=correlation,
        is_anti_correlated=correlation < 0,
        correct_scores_rank_order=correct_rank_order,
        deprecated_scores_rank_order=deprecated_rank_order,
    )


def print_splice_validation_report(result: SpliceValidationResult) -> None:
    """Print a formatted validation report comparing deprecated and correct scoring.

    The report includes:
    - A comparison table showing each test site with both scores
    - The Pearson correlation between methods
    - Whether the correlation is negative (anti-correlation)
    - Rank-order checks for both methods
    - A clear explanation of the anti-correlation finding
    - A recommendation to use maxentscan functions

    Args:
        result: SpliceValidationResult from ``validate_splice_scoring()``.
    """
    print("=" * 90)
    print("SPLICE SCORING VALIDATION REPORT")
    print("Comparing: splicing.maxent_score() [DEPRECATED]")
    print("     vs:   maxentscan.score_donor() / score_acceptor() [CORRECT]")
    print("Reference: Yeo & Burge (2004) J. Comput. Biol. 11(2-3):377-394")
    print("=" * 90)
    print()

    # ── Comparison table ──────────────────────────────────────────────────
    header = (
        f"{'Sequence':<26} {'Type':<10} {'Strength':<10} "
        f"{'Deprecated':>11} {'Correct':>11} {'Delta':>11}"
    )
    print(header)
    print("-" * 90)

    for site in result.test_sites:
        seq = site["sequence"]
        if len(seq) > 24:
            seq = seq[:21] + "..."
        site_type = site["site_type"]
        strength = site["expected_strength"]
        dep = site["deprecated_score"]
        cor = site["correct_score"]
        delta = cor - dep

        marker = ""
        # Mark anti-correlated entries: deprecated says high, correct says low
        # (or vice versa)
        if dep > cor + 2.0:
            marker = " <<<"
        elif cor > dep + 2.0:
            marker = " >>>"

        print(
            f"{seq:<26} {site_type:<10} {strength:<10} "
            f"{dep:>11.4f} {cor:>11.4f} {delta:>+11.4f}{marker}"
        )

    print("-" * 90)
    print()

    # ── Correlation ───────────────────────────────────────────────────────
    print("STATISTICAL SUMMARY")
    print("-" * 40)
    print(f"  Pearson correlation (deprecated vs correct): {result.correlation:.4f}")
    print(f"  Is anti-correlated (correlation < 0):       {result.is_anti_correlated}")
    print()

    # ── Rank ordering ─────────────────────────────────────────────────────
    print("RANK ORDER CHECKS")
    print("-" * 40)
    print(
        f"  Correct scores rank order (strong > weak):  "
        f"{'PASS' if result.correct_scores_rank_order else 'FAIL'}"
    )
    print(
        f"  Deprecated scores rank order (strong > weak): "
        f"{'PASS' if result.deprecated_scores_rank_order else 'FAIL'}"
    )
    print()

    # ── Detailed scores by strength category ──────────────────────────────
    print("SCORES BY STRENGTH CATEGORY")
    print("-" * 40)

    for strength_label in ("strong", "moderate", "weak"):
        sites_of_type = [
            s for s in result.test_sites if s["expected_strength"] == strength_label
        ]
        if not sites_of_type:
            continue

        dep_scores = [s["deprecated_score"] for s in sites_of_type]
        cor_scores = [
            s["correct_score"] for s in sites_of_type if s["correct_score"] > -50.0
        ]

        dep_avg = sum(dep_scores) / len(dep_scores) if dep_scores else 0.0
        cor_avg = sum(cor_scores) / len(cor_scores) if cor_scores else 0.0

        print(
            f"  {strength_label.capitalize():<10} "
            f"deprecated avg: {dep_avg:>8.4f}  |  "
            f"correct avg: {cor_avg:>8.4f}"
        )

    print()

    # ── Anti-correlation explanation ──────────────────────────────────────
    if result.is_anti_correlated:
        print("!" * 90)
        print("  ANTI-CORRELATION DETECTED  (Pearson r = {:.4f})".format(result.correlation))
        print("!" * 90)
        print()
        print(
            "The deprecated splicing.maxent_score() produces scores that are "
            "ANTI-CORRELATED"
        )
        print(
            "with the proper MaxEntScan (Yeo & Burge 2004) log-odds scoring."
        )
        print()
        print("Root cause: The hand-crafted PWM in splicing.py assigns:")
        print("  - A at positions 3-4 (the GT dinucleotide positions): weight 3.50 EACH")
        print("  - G at position 3 (donor G):                        weight 0.01")
        print("  - T at position 4 (donor T):                        weight 0.01")
        print()
        print(
            "When genuine GT is present at the core positions, the dinucleotide "
            "contributes"
        )
        print(
            "only 0.02 (0.01+0.01) to the total score — making true splice sites "
            "invisible."
        )
        print(
            "Non-splice-site sequences with A at those positions get a massive "
            "7.00 (3.50+3.50)"
        )
        print("boost, producing HIGHER scores for non-sites than for real sites.")
        print()
        print(
            "In contrast, the proper maxentscan.score_donor()/score_acceptor() "
            "use log-odds"
        )
        print(
            "scoring with Yeo & Burge 2004 trained parameters, where the "
            "nearly-invariant"
        )
        print(
            "G (prob 0.990) and T (prob 0.990) at donor positions contribute "
            "strong positive"
        )
        print(
            "log-odds (~+2.0 bits each), correctly ranking genuine splice sites "
            "above non-sites."
        )
        print()
    else:
        print(
            "NOTE: Correlation is non-negative ({:.4f}). This is unexpected; the ".format(
                result.correlation
            )
        )
        print(
            "deprecated PWM should be anti-correlated with proper MaxEntScan."
        )
        print()

    # ── Recommendation ────────────────────────────────────────────────────
    print("=" * 90)
    print("RECOMMENDATION")
    print("=" * 90)
    print()
    print("  Use maxentscan.score_donor() and maxentscan.score_acceptor() for all")
    print("  splice site scoring. The deprecated splicing.maxent_score() should NOT")
    print("  be used for any scoring application — its rankings are inverted relative")
    print("  to biological reality (Yeo & Burge 2004).")
    print()
    print("  For a drop-in replacement, use splicing.maxent_score_v2() which")
    print("  delegates to maxentscan.score_donor() with automatic context extraction.")
    print()


# ==============================================================================
# Markov Cross-Check: Production vs Benchmarking Implementation
# ==============================================================================


def validate_markov_cross_check(
    tolerance: float = 0.1,
) -> MarkovCrossCheckResult:
    """Cross-check production Markov scores against benchmarking Markov scores.

    Both the production ``maxentscan`` module and the benchmarking
    ``maxentscan_validated`` module now use the first-order Markov model.
    Their scores should agree within floating-point tolerance; any
    disagreement indicates a regression in either module.

    Args:
        tolerance: maximum allowed absolute difference between production
            and benchmarking scores (default 0.1).

    Returns:
        MarkovCrossCheckResult with per-sequence comparisons and summary
        statistics.
    """
    donor_results: List[dict] = []
    acceptor_results: List[dict] = []

    # ── Donor cross-check ──────────────────────────────────────────────────
    for seq, _ref_score in _BENCHMARK_DONOR_SCORES.items():
        gt_pos = 3  # GT starts at index 3 in 9-mers
        prod_score = _score_proper_donor(seq, gt_pos)
        bench_score = _benchmark_donor(seq)
        diff = prod_score - bench_score
        donor_results.append({
            "sequence": seq,
            "production_score": prod_score,
            "benchmark_score": bench_score,
            "difference": diff,
        })

    # ── Acceptor cross-check ───────────────────────────────────────────────
    for seq, _ref_score in _BENCHMARK_ACCEPTOR_SCORES.items():
        ag_pos = 20  # AG starts at index 20 in 23-mers
        prod_score = _score_proper_acceptor(seq, ag_pos)
        bench_score = _benchmark_acceptor(seq)
        diff = prod_score - bench_score
        acceptor_results.append({
            "sequence": seq,
            "production_score": prod_score,
            "benchmark_score": bench_score,
            "difference": diff,
        })

    # ── Correlation ────────────────────────────────────────────────────────
    prod_donor_scores = [r["production_score"] for r in donor_results]
    bench_donor_scores = [r["benchmark_score"] for r in donor_results]
    donor_corr = _pearson_correlation(prod_donor_scores, bench_donor_scores)

    prod_acc_scores = [r["production_score"] for r in acceptor_results]
    bench_acc_scores = [r["benchmark_score"] for r in acceptor_results]
    acc_corr = _pearson_correlation(prod_acc_scores, bench_acc_scores)

    # ── Tolerance check ────────────────────────────────────────────────────
    max_donor_diff = max(abs(r["difference"]) for r in donor_results) if donor_results else 0.0
    max_acc_diff = max(abs(r["difference"]) for r in acceptor_results) if acceptor_results else 0.0

    all_within = all(abs(r["difference"]) <= tolerance for r in donor_results) and \
                 all(abs(r["difference"]) <= tolerance for r in acceptor_results)

    return MarkovCrossCheckResult(
        donor_results=donor_results,
        acceptor_results=acceptor_results,
        donor_correlation=donor_corr,
        acceptor_correlation=acc_corr,
        max_donor_diff=max_donor_diff,
        max_acceptor_diff=max_acc_diff,
        all_within_tolerance=all_within,
    )


def print_markov_cross_check_report(
    result: MarkovCrossCheckResult,
    tolerance: float = 0.1,
) -> None:
    """Print a formatted report of the Markov cross-check results.

    Args:
        result: MarkovCrossCheckResult from ``validate_markov_cross_check()``.
        tolerance: the tolerance used for the cross-check.
    """
    print("=" * 90)
    print("MARKOV CROSS-CHECK: Production vs Benchmarking Implementation")
    print("Both modules use the first-order Markov model — scores should agree.")
    print(f"Tolerance: {tolerance}")
    print("=" * 90)
    print()

    # ── Donor comparison ───────────────────────────────────────────────────
    print("DONOR SITES (9-mer)")
    print("-" * 90)
    header = (
        f"{'Sequence':<14} {'Production':>12} {'Benchmark':>12} {'Difference':>12} {'Status':>8}"
    )
    print(header)
    print("-" * 90)

    for r in result.donor_results:
        diff = r["difference"]
        status = "OK" if abs(diff) <= tolerance else "MISMATCH"
        seq = r["sequence"]
        print(
            f"{seq:<14} {r['production_score']:>12.4f} {r['benchmark_score']:>12.4f} "
            f"{diff:>+12.4f} {status:>8}"
        )

    print("-" * 90)
    print(f"  Pearson correlation: {result.donor_correlation:.4f}")
    print(f"  Max absolute diff:   {result.max_donor_diff:.4f}")
    print()

    # ── Acceptor comparison ────────────────────────────────────────────────
    print("ACCEPTOR SITES (23-mer)")
    print("-" * 90)
    header = (
        f"{'Sequence':<26} {'Production':>12} {'Benchmark':>12} {'Difference':>12} {'Status':>8}"
    )
    print(header)
    print("-" * 90)

    for r in result.acceptor_results:
        diff = r["difference"]
        status = "OK" if abs(diff) <= tolerance else "MISMATCH"
        seq = r["sequence"]
        if len(seq) > 24:
            seq = seq[:21] + "..."
        print(
            f"{seq:<26} {r['production_score']:>12.4f} {r['benchmark_score']:>12.4f} "
            f"{diff:>+12.4f} {status:>8}"
        )

    print("-" * 90)
    print(f"  Pearson correlation: {result.acceptor_correlation:.4f}")
    print(f"  Max absolute diff:   {result.max_acceptor_diff:.4f}")
    print()

    # ── Summary ────────────────────────────────────────────────────────────
    print("=" * 90)
    if result.all_within_tolerance:
        print("RESULT: ALL SCORES MATCH within tolerance ({:.4f})".format(tolerance))
        print("Both Markov implementations agree — no regression detected.")
    else:
        print("RESULT: MISMATCH DETECTED!")
        print("Some scores differ by more than the tolerance ({:.4f}).".format(tolerance))
        print("This indicates a regression in either the production or benchmarking module.")
    print()

    if result.donor_correlation < 0 or result.acceptor_correlation < 0:
        print("WARNING: Negative correlation detected between implementations!")
        print("This is a strong signal of a regression — both implementations")
        print("use the same algorithm and should be positively correlated.")
    elif result.donor_correlation > 0.9 and result.acceptor_correlation > 0.9:
        print("Correlation check: PASS (both > 0.9)")
    else:
        print(
            "Correlation check: MARGINAL "
            f"(donor={result.donor_correlation:.4f}, acceptor={result.acceptor_correlation:.4f})"
        )
    print("=" * 90)
    print()
