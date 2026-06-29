"""
BioCompiler Benchmarking — Head-to-Head Comparison Module
==========================================================

Core data structures and comparison functions for benchmarking
BioCompiler against DNAchisel on a per-gene basis, with statistical
analysis across multiple genes.

Provides:
  - ``BenchmarkResult``: Per-gene comparison dataclass with metrics
    for both BioCompiler and DNAchisel (including sequences and
    restriction site counts)
  - ``ComparisonSummary``: Aggregate statistical analysis with mean,
    median, and p-value for each metric
  - ``run_comparison``: Single-gene head-to-head comparison
  - ``compare_results``: Multi-gene statistical analysis
  - ``ALL_CONSTRAINTS``: Standard constraint names for visualization

Design:
    All metrics are computed using BioCompiler's validated evaluators for
    fairness — both tools' outputs are evaluated with the same CAI
    computation (``compute_cai_validated``), GC measurement, and restriction
    site scanner.  DNAchisel's own CAI output is NOT trusted.
"""

from __future__ import annotations

import logging
import math
import time
from dataclasses import dataclass, field
from typing import Optional

from ..organisms import resolve_organism
from biocompiler.sequence.scanner import gc_content
from .metrics import compute_cai_validated

logger = logging.getLogger(__name__)

__all__ = [
    "BenchmarkResult",
    "ComparisonSummary",
    "run_comparison",
    "compare_results",
    "ALL_CONSTRAINTS",
]

# ─── Constraint Names ─────────────────────────────────────────────────
# Used by the visualization module for constraint satisfaction plots.

CONSTRAINT_TRANSLATION = "Translation"
CONSTRAINT_GC_RANGE = "GC Range"
CONSTRAINT_NO_RS = "No Restriction Sites"
CONSTRAINT_NO_GT = "No GT Dinucleotides"
CONSTRAINT_LOW_CPG = "Low CpG Ratio"

ALL_CONSTRAINTS: list[str] = [
    CONSTRAINT_TRANSLATION,
    CONSTRAINT_GC_RANGE,
    CONSTRAINT_NO_RS,
    CONSTRAINT_NO_GT,
    CONSTRAINT_LOW_CPG,
]

# Default GC range for constraint checking
_DEFAULT_GC_LO = 0.30
_DEFAULT_GC_HI = 0.70

# Default restriction enzymes
_DEFAULT_ENZYMES = ["EcoRI", "BamHI", "HindIII", "XhoI"]

# CpG ratio threshold
_CPG_THRESHOLD = 0.6

# GC midpoint for "closer to target" comparison
_GC_MIDPOINT = 0.50

# Epsilon for comparing floating-point differences
_EPSILON = 0.001

# Runtime noise margin (10%)
_RUNTIME_MARGIN = 0.10


# ─── Helper: count GT dinucleotides ──────────────────────────────────

def _count_gt(seq: str) -> int:
    """Count GT dinucleotides in a sequence."""
    return sum(1 for i in range(len(seq) - 1) if seq[i:i + 2] == "GT")


# ─── Helper: compute CpG obs/exp ratio ──────────────────────────────

def _compute_cpg_ratio(seq: str) -> float:
    """Compute CpG Obs/Exp ratio for the full sequence."""
    c = seq.count("C")
    g = seq.count("G")
    cg = sum(1 for i in range(len(seq) - 1) if seq[i:i + 2] == "CG")
    expected = (c * g) / len(seq) if len(seq) > 0 else 0
    return cg / expected if expected > 0 else 0.0


# ─── Helper: count restriction sites (both strands) ─────────────────

def _count_restriction_sites_both_strands(
    seq: str,
    enzymes: list[str] | None = None,
) -> int:
    """Count total restriction sites in sequence (forward + reverse complement)."""
    from biocompiler.shared.constants import RESTRICTION_ENZYMES, reverse_complement

    if enzymes is None:
        enzymes = _DEFAULT_ENZYMES

    count = 0
    seq_upper = seq.upper()
    for enz in enzymes:
        site = RESTRICTION_ENZYMES.get(enz, "")
        if not site:
            continue
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            continue
        # Forward strand
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
        # Reverse complement strand
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:  # Avoid double-counting palindromes
            start = 0
            while True:
                pos = seq_upper.find(site_rc, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1
    return count


# ─── Helper: count restriction sites (forward strand only) ──────────

def _count_restriction_sites(seq: str, enzymes: list[str] | None = None) -> int:
    """Count total restriction sites in sequence (forward strand only)."""
    from biocompiler.sequence.restriction_sites import get_recognition_site

    if enzymes is None:
        enzymes = _DEFAULT_ENZYMES
    total = 0
    for enz in enzymes:
        site = get_recognition_site(enz)
        if site is None:
            continue
        pos = seq.find(site)
        while pos != -1:
            total += 1
            pos = seq.find(site, pos + 1)
    return total


# ─── Helper: evaluate constraint satisfaction ────────────────────────

def _evaluate_constraints(
    sequence: str,
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
    gc_lo: float = _DEFAULT_GC_LO,
    gc_hi: float = _DEFAULT_GC_HI,
) -> dict[str, bool]:
    """Evaluate which constraints are satisfied by a sequence.

    Returns a dict mapping constraint name to bool (True = satisfied).
    """
    from biocompiler.expression.translation import translate

    constraints: dict[str, bool] = {}

    # 1. Translation fidelity
    translated = translate(sequence, to_stop=True)
    constraints[CONSTRAINT_TRANSLATION] = (translated == protein)

    # 2. GC content in range
    gc = gc_content(sequence)
    constraints[CONSTRAINT_GC_RANGE] = (gc_lo <= gc <= gc_hi)

    # 3. No restriction sites
    rs_count = _count_restriction_sites(sequence, enzymes)
    constraints[CONSTRAINT_NO_RS] = (rs_count == 0)

    # 4. No GT dinucleotides
    gt_count = _count_gt(sequence)
    constraints[CONSTRAINT_NO_GT] = (gt_count == 0)

    # 5. Low CpG ratio
    cpg = _compute_cpg_ratio(sequence)
    constraints[CONSTRAINT_LOW_CPG] = (cpg < _CPG_THRESHOLD)

    return constraints


# ─── Core Data Structures ────────────────────────────────────────────


@dataclass
class BenchmarkResult:
    """Per-gene benchmark result comparing BioCompiler vs DNAchisel.

    Captures metrics from both BioCompiler and DNAchisel for a single
    protein/organism pair, enabling direct comparison on:
      - CAI (Codon Adaptation Index — higher is better)
      - GC content (target range compliance)
      - Restriction site count (fewer is better)
      - Runtime (faster is better)
      - Raw sequences for further analysis
      - Constraint satisfaction for each tool

    Attributes:
        gene_name: Name of the gene (e.g., "HBB", "EGFP").
        organism: Target organism name (e.g., "Homo_sapiens").
        biocompiler_cai: CAI score from BioCompiler optimization.
        dnachisel_cai: CAI score from DNAchisel optimization (None if unavailable).
        biocompiler_gc: GC content fraction from BioCompiler.
        dnachisel_gc: GC content fraction from DNAchisel (None if unavailable).
        biocompiler_restriction_sites: Restriction site count from BioCompiler.
        dnachisel_restriction_sites: Restriction site count from DNAchisel
            (None if unavailable).
        biocompiler_runtime: Optimization runtime in seconds (BioCompiler).
        dnachisel_runtime: Optimization runtime in seconds (DNAchisel,
            None if unavailable).
        biocompiler_sequence: Optimized DNA sequence from BioCompiler.
        dnachisel_sequence: Optimized DNA sequence from DNAchisel
            (None if unavailable).
        biocompiler_success: Whether BioCompiler optimization succeeded.
        dnachisel_success: Whether DNAchisel optimization succeeded
            (None if unavailable).
        biocompiler_error: Error message if BioCompiler failed, else None.
        dnachisel_error: Error message if DNAchisel failed, else None.
        enzyme_list: List of restriction enzyme names used in the comparison.
        protein_length: Length of the protein in amino acids.
        constraints_biocompiler: Dict mapping constraint name → satisfied (bool).
        constraints_dnachisel: Dict mapping constraint name → satisfied (bool),
            or None if DNAchisel was unavailable/failed.
    """

    gene_name: str
    organism: str
    biocompiler_cai: float
    dnachisel_cai: float | None
    biocompiler_gc: float
    dnachisel_gc: float | None
    biocompiler_restriction_sites: int
    dnachisel_restriction_sites: int | None
    biocompiler_runtime: float
    dnachisel_runtime: float | None
    biocompiler_sequence: str
    dnachisel_sequence: str | None
    biocompiler_success: bool = True
    dnachisel_success: bool | None = None
    biocompiler_error: str | None = None
    dnachisel_error: str | None = None
    enzyme_list: list[str] = field(default_factory=list)
    protein_length: int = 0
    constraints_biocompiler: dict[str, bool] = field(default_factory=dict)
    constraints_dnachisel: dict[str, bool] | None = None

    # ── Aliases for visualization module compatibility ──
    # The visualization module uses cai_biocompiler, gc_biocompiler,
    # runtime_biocompiler naming convention.

    @property
    def cai_biocompiler(self) -> float:
        """Alias for biocompiler_cai (visualization compatibility)."""
        return self.biocompiler_cai

    @property
    def cai_dnachisel(self) -> float | None:
        """Alias for dnachisel_cai (visualization compatibility)."""
        return self.dnachisel_cai

    @property
    def gc_biocompiler(self) -> float:
        """Alias for biocompiler_gc (visualization compatibility)."""
        return self.biocompiler_gc

    @property
    def gc_dnachisel(self) -> float | None:
        """Alias for dnachisel_gc (visualization compatibility)."""
        return self.dnachisel_gc

    @property
    def runtime_biocompiler(self) -> float:
        """Alias for biocompiler_runtime (visualization compatibility)."""
        return self.biocompiler_runtime

    @property
    def runtime_dnachisel(self) -> float | None:
        """Alias for dnachisel_runtime (visualization compatibility)."""
        return self.dnachisel_runtime


@dataclass
class ComparisonSummary:
    """Aggregate statistical analysis across multiple BenchmarkResults.

    For each metric (CAI, GC, restriction sites, runtime), provides:
      - Mean, median for both BioCompiler and DNAchisel
      - Paired Wilcoxon signed-rank p-value (or t-test p-value if n >= 30)
      - Win count (how many genes each tool wins on that metric)

    Also provides high-level counts and constraint satisfaction rates.

    Attributes:
        total_comparisons: Number of BenchmarkResults analyzed.
        successful_comparisons: Number where both tools succeeded.
        genes_with_dnachisel: Number of genes where DNAchisel results are available.
        cai_mean_biocompiler: Mean CAI across genes (BioCompiler).
        cai_mean_dnachisel: Mean CAI across genes (DNAchisel).
        cai_median_biocompiler: Median CAI across genes (BioCompiler).
        cai_median_dnachisel: Median CAI across genes (DNAchisel).
        cai_p_value: Statistical significance of CAI difference.
        cai_biocompiler_wins: Number of genes where BioCompiler has higher CAI.
        cai_dnachisel_wins: Number of genes where DNAchisel has higher CAI.
        gc_mean_biocompiler: Mean GC content (BioCompiler).
        gc_mean_dnachisel: Mean GC content (DNAchisel).
        gc_median_biocompiler: Median GC content (BioCompiler).
        gc_median_dnachisel: Median GC content (DNAchisel).
        gc_p_value: Statistical significance of GC difference.
        gc_biocompiler_wins: Number of genes where BioCompiler GC is closer to target.
        gc_dnachisel_wins: Number of genes where DNAchisel GC is closer to target.
        rs_mean_biocompiler: Mean restriction site count (BioCompiler).
        rs_mean_dnachisel: Mean restriction site count (DNAchisel).
        rs_median_biocompiler: Median restriction site count (BioCompiler).
        rs_median_dnachisel: Median restriction site count (DNAchisel).
        rs_p_value: Statistical significance of restriction site difference.
        rs_biocompiler_wins: Number of genes where BioCompiler has fewer sites.
        rs_dnachisel_wins: Number of genes where DNAchisel has fewer sites.
        runtime_mean_biocompiler: Mean runtime in seconds (BioCompiler).
        runtime_mean_dnachisel: Mean runtime in seconds (DNAchisel).
        runtime_median_biocompiler: Median runtime (BioCompiler).
        runtime_median_dnachisel: Median runtime (DNAchisel).
        runtime_p_value: Statistical significance of runtime difference.
        runtime_biocompiler_wins: Number of genes where BioCompiler is faster.
        runtime_dnachisel_wins: Number of genes where DNAchisel is faster.
        constraint_satisfaction_rate_biocompiler: Overall constraint satisfaction
            rate for BioCompiler across all genes and constraints.
        constraint_satisfaction_rate_dnachisel: Overall constraint satisfaction
            rate for DNAchisel across all genes and constraints.
    """

    total_comparisons: int = 0
    successful_comparisons: int = 0
    genes_with_dnachisel: int = 0
    # CAI
    cai_mean_biocompiler: float = 0.0
    cai_mean_dnachisel: float = 0.0
    cai_median_biocompiler: float = 0.0
    cai_median_dnachisel: float = 0.0
    cai_p_value: float = 1.0
    cai_biocompiler_wins: int = 0
    cai_dnachisel_wins: int = 0
    # GC
    gc_mean_biocompiler: float = 0.0
    gc_mean_dnachisel: float = 0.0
    gc_median_biocompiler: float = 0.0
    gc_median_dnachisel: float = 0.0
    gc_p_value: float = 1.0
    gc_biocompiler_wins: int = 0
    gc_dnachisel_wins: int = 0
    # Restriction sites
    rs_mean_biocompiler: float = 0.0
    rs_mean_dnachisel: float = 0.0
    rs_median_biocompiler: float = 0.0
    rs_median_dnachisel: float = 0.0
    rs_p_value: float = 1.0
    rs_biocompiler_wins: int = 0
    rs_dnachisel_wins: int = 0
    # Runtime
    runtime_mean_biocompiler: float = 0.0
    runtime_mean_dnachisel: float = 0.0
    runtime_median_biocompiler: float = 0.0
    runtime_median_dnachisel: float = 0.0
    runtime_p_value: float = 1.0
    runtime_biocompiler_wins: int = 0
    runtime_dnachisel_wins: int = 0
    # Constraint satisfaction rates
    constraint_satisfaction_rate_biocompiler: float = 0.0
    constraint_satisfaction_rate_dnachisel: float = 0.0

    # Backward-compatible aliases for the older ComparisonSummary API
    @property
    def total_genes(self) -> int:
        """Alias for total_comparisons."""
        return self.total_comparisons

    @property
    def avg_cai_biocompiler(self) -> float:
        """Alias for cai_mean_biocompiler."""
        return self.cai_mean_biocompiler

    @property
    def avg_cai_dnachisel(self) -> float:
        """Alias for cai_mean_dnachisel."""
        return self.cai_mean_dnachisel

    @property
    def avg_gc_biocompiler(self) -> float:
        """Alias for gc_mean_biocompiler."""
        return self.gc_mean_biocompiler

    @property
    def avg_gc_dnachisel(self) -> float:
        """Alias for gc_mean_dnachisel."""
        return self.gc_mean_dnachisel

    @property
    def avg_runtime_biocompiler(self) -> float:
        """Alias for runtime_mean_biocompiler."""
        return self.runtime_mean_biocompiler

    @property
    def avg_runtime_dnachisel(self) -> float:
        """Alias for runtime_mean_dnachisel."""
        return self.runtime_mean_dnachisel

    @property
    def cai_wins(self) -> int:
        """Alias for cai_biocompiler_wins."""
        return self.cai_biocompiler_wins

    @property
    def cai_losses(self) -> int:
        """Alias for cai_dnachisel_wins."""
        return self.cai_dnachisel_wins

    @property
    def cai_ties(self) -> int:
        """Number of CAI ties."""
        n = self.total_comparisons
        return n - self.cai_biocompiler_wins - self.cai_dnachisel_wins


# ─── Statistical helpers ─────────────────────────────────────────────


def _median(values: list[float]) -> float:
    """Compute the median of a list of floats."""
    if not values:
        return 0.0
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    if n % 2 == 1:
        return sorted_vals[n // 2]
    return (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2.0


def _mean(values: list[float]) -> float:
    """Compute the mean of a list of floats."""
    if not values:
        return 0.0
    return sum(values) / len(values)


def _normal_cdf(x: float) -> float:
    """Approximate the standard normal CDF (Abramowitz & Stegun)."""
    if x < -8.0:
        return 0.0
    if x > 8.0:
        return 1.0
    t = 1.0 / (1.0 + 0.2316419 * abs(x))
    d = 0.3989422804014327  # 1/sqrt(2*pi)
    p = d * math.exp(-x * x / 2.0) * (
        t * (0.319381530 + t * (-0.356563782 + t * (
            1.781477937 + t * (-1.821255978 + t * 1.330274429)
        )))
    )
    if x > 0:
        return 1.0 - p
    return p


def _wilcoxon_signed_rank_p_value(differences: list[float]) -> float:
    """Approximate Wilcoxon signed-rank test p-value (two-sided).

    Uses a normal approximation suitable for n >= 10. For smaller
    samples, returns 1.0 (not enough data for significance).
    """
    diffs = [d for d in differences if abs(d) > 1e-12]
    n = len(diffs)

    if n < 10:
        return 1.0

    abs_diffs = sorted(enumerate(diffs), key=lambda x: abs(x[1]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i + 1
        while j < n and abs(abs_diffs[j][1]) - abs(abs_diffs[i][1]) < 1e-12:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[abs_diffs[k][0]] = avg_rank
        i = j

    w_plus = sum(ranks[i] for i in range(n) if diffs[i] > 0)
    mean_w = n * (n + 1) / 4.0
    var_w = n * (n + 1) * (2 * n + 1) / 24.0

    if var_w <= 0:
        return 1.0

    z = abs(w_plus - mean_w) / math.sqrt(var_w)
    p = 2.0 * (1.0 - _normal_cdf(abs(z)))
    return max(0.0, min(1.0, p))


def _paired_t_test_p_value(differences: list[float]) -> float:
    """Compute paired t-test p-value (two-sided)."""
    n = len(differences)
    if n < 2:
        return 1.0

    mean_d = _mean(differences)
    var_d = sum((d - mean_d) ** 2 for d in differences) / (n - 1)
    if var_d <= 1e-15:
        return 1.0

    se = math.sqrt(var_d / n)
    t_stat = abs(mean_d / se)
    p = 2.0 * (1.0 - _normal_cdf(t_stat))
    return max(0.0, min(1.0, p))


# ─── Single-Gene Comparison ──────────────────────────────────────────


def run_comparison(
    gene_name: str,
    protein: str,
    organism: str = "Homo_sapiens",
    enzymes: list[str] | None = None,
    gc_lo: float = _DEFAULT_GC_LO,
    gc_hi: float = _DEFAULT_GC_HI,
) -> BenchmarkResult:
    """Run both BioCompiler and DNAchisel on a single gene and collect metrics.

    This is the primary entry point for per-gene head-to-head comparison.
    All metrics are computed using BioCompiler's own evaluators for fairness.

    Args:
        gene_name: Name of the gene (e.g., "HBB", "EGFP").
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism for codon usage (e.g., "Homo_sapiens").
        enzymes: Restriction enzymes to avoid. If None, uses defaults.
        gc_lo: Minimum GC content fraction.
        gc_hi: Maximum GC content fraction.

    Returns:
        BenchmarkResult with metrics from both optimizers.
        DNAchisel fields will be None if DNAchisel is unavailable.
    """
    from biocompiler.optimizer import optimize_sequence

    enzymes = enzymes or _DEFAULT_ENZYMES

    # Resolve organism name once for consistent alias handling
    resolved_org = resolve_organism(organism)

    # ── BioCompiler ──
    bc_cai = 0.0
    bc_gc = 0.0
    bc_rs = 0
    bc_runtime = 0.0
    bc_seq = ""
    bc_success = True
    bc_error: str | None = None
    bc_constraints: dict[str, bool] = {}

    t0 = time.perf_counter()
    try:
        bc_result = optimize_sequence(
            target_protein=protein,
            organism=resolved_org,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
        )
        bc_seq = bc_result.sequence
        # Use validated CAI for fair comparison — do NOT trust optimizer's CAI
        bc_cai = compute_cai_validated(bc_seq, resolved_org)
        # Use validated gc_content for fair comparison
        bc_gc = gc_content(bc_seq)
        bc_rs = _count_restriction_sites_both_strands(bc_seq, enzymes)
    except Exception as exc:
        bc_success = False
        bc_error = str(exc)
        logger.error("BioCompiler failed for %s: %s", gene_name, exc)
    bc_runtime = time.perf_counter() - t0

    # Evaluate BioCompiler constraints
    if bc_success and bc_seq:
        bc_constraints = _evaluate_constraints(
            bc_seq, protein, resolved_org, enzymes, gc_lo, gc_hi,
        )
    else:
        bc_constraints = {c: False for c in ALL_CONSTRAINTS}

    # ── DNAchisel ──
    dc_cai: float | None = None
    dc_gc: float | None = None
    dc_rs: int | None = None
    dc_runtime: float | None = None
    dc_seq: str | None = None
    dc_success: bool | None = None
    dc_error: str | None = None
    dc_constraints: dict[str, bool] | None = None

    try:
        from .dnachisel_adapter import DNAchiselAdapter, is_dnachisel_available

        if is_dnachisel_available():
            adapter = DNAchiselAdapter()
            t1 = time.perf_counter()
            dc_result = adapter.optimize(
                protein=protein,
                organism=resolved_org,
                constraints=[
                    {"type": "gc_range", "gc_lo": gc_lo, "gc_hi": gc_hi},
                    {"type": "avoid_restriction", "enzymes": enzymes},
                ],
            )
            dc_runtime = time.perf_counter() - t1

            if dc_result.success:
                dc_seq = dc_result.sequence
                # Use BioCompiler's validated evaluators for fairness —
                # do NOT trust DNAchisel's own CAI / GC / RS output.
                # This matches the module design principle: both tools'
                # outputs are evaluated with the same metrics.
                dc_cai = compute_cai_validated(dc_seq, resolved_org)
                dc_gc = gc_content(dc_seq)
                dc_rs = _count_restriction_sites_both_strands(dc_seq, enzymes)
                dc_success = True
                dc_constraints = _evaluate_constraints(
                    dc_seq, protein, resolved_org, enzymes, gc_lo, gc_hi,
                )
            else:
                dc_success = False
                dc_error = dc_result.error
                dc_constraints = {c: False for c in ALL_CONSTRAINTS}
    except ImportError as exc:
        logger.debug(
            "DNAchisel not available for %s comparison: %s", gene_name, exc
        )
    except Exception as exc:
        dc_success = False
        dc_error = str(exc)
        logger.error("DNAchisel failed for %s: %s", gene_name, exc)

    return BenchmarkResult(
        gene_name=gene_name,
        organism=organism,
        biocompiler_cai=bc_cai,
        dnachisel_cai=dc_cai,
        biocompiler_gc=bc_gc,
        dnachisel_gc=dc_gc,
        biocompiler_restriction_sites=bc_rs,
        dnachisel_restriction_sites=dc_rs,
        biocompiler_runtime=bc_runtime,
        dnachisel_runtime=dc_runtime,
        biocompiler_sequence=bc_seq,
        dnachisel_sequence=dc_seq,
        biocompiler_success=bc_success,
        dnachisel_success=dc_success,
        biocompiler_error=bc_error,
        dnachisel_error=dc_error,
        enzyme_list=enzymes,
        protein_length=len(protein),
        constraints_biocompiler=bc_constraints,
        constraints_dnachisel=dc_constraints,
    )


# ─── Multi-Gene Statistical Comparison ───────────────────────────────


def compare_results(
    results: list[BenchmarkResult],
) -> ComparisonSummary:
    """Perform statistical analysis across multiple BenchmarkResults.

    Computes mean, median, and paired significance tests (Wilcoxon
    signed-rank for n < 30, paired t-test for n >= 30) for each
    metric across all provided BenchmarkResults.

    Only results where DNAchisel data is available are included in
    statistical comparisons (BioCompiler-only results are counted
    as BioCompiler wins).

    Args:
        results: List of BenchmarkResult objects from individual comparisons.

    Returns:
        ComparisonSummary with aggregate statistics and significance tests.
    """
    if not results:
        return ComparisonSummary()

    n = len(results)
    summary = ComparisonSummary(total_comparisons=n)

    # Genes with DNAchisel data
    dc_available = [
        r for r in results
        if r.dnachisel_cai is not None
    ]
    summary.genes_with_dnachisel = len(dc_available)

    # Paired results (both tools succeeded with data)
    paired = [
        r for r in results
        if r.biocompiler_success
        and r.dnachisel_cai is not None
        and r.dnachisel_success is True
    ]
    summary.successful_comparisons = len(paired)

    # ── Constraint satisfaction rates ──
    bc_total_constraints = 0
    bc_satisfied = 0
    dc_total_constraints = 0
    dc_satisfied = 0

    for r in results:
        for satisfied in r.constraints_biocompiler.values():
            bc_total_constraints += 1
            if satisfied:
                bc_satisfied += 1
        if r.constraints_dnachisel is not None:
            for satisfied in r.constraints_dnachisel.values():
                dc_total_constraints += 1
                if satisfied:
                    dc_satisfied += 1

    summary.constraint_satisfaction_rate_biocompiler = (
        round(bc_satisfied / bc_total_constraints, 4)
        if bc_total_constraints > 0
        else 0.0
    )
    summary.constraint_satisfaction_rate_dnachisel = (
        round(dc_satisfied / dc_total_constraints, 4)
        if dc_total_constraints > 0
        else 0.0
    )

    if not paired:
        # Still compute BioCompiler-only stats
        bc_cais = [r.biocompiler_cai for r in results if r.biocompiler_success]
        bc_gcs = [r.biocompiler_gc for r in results if r.biocompiler_success]
        bc_rss = [
            float(r.biocompiler_restriction_sites)
            for r in results
            if r.biocompiler_success
        ]
        bc_times = [
            r.biocompiler_runtime for r in results if r.biocompiler_success
        ]

        summary.cai_mean_biocompiler = round(_mean(bc_cais), 4)
        summary.cai_median_biocompiler = round(_median(bc_cais), 4)
        summary.gc_mean_biocompiler = round(_mean(bc_gcs), 4)
        summary.gc_median_biocompiler = round(_median(bc_gcs), 4)
        summary.rs_mean_biocompiler = round(_mean(bc_rss), 2)
        summary.rs_median_biocompiler = round(_median(bc_rss), 2)
        summary.runtime_mean_biocompiler = round(_mean(bc_times), 4)
        summary.runtime_median_biocompiler = round(_median(bc_times), 4)

        # All wins go to BioCompiler by default
        summary.cai_biocompiler_wins = n
        summary.gc_biocompiler_wins = n
        summary.rs_biocompiler_wins = n
        summary.runtime_biocompiler_wins = n
        return summary

    # ── CAI comparison ──
    bc_cais = [r.biocompiler_cai for r in paired]
    dc_cais = [r.dnachisel_cai for r in paired]  # type: ignore[misc]

    summary.cai_mean_biocompiler = round(_mean(bc_cais), 4)
    summary.cai_mean_dnachisel = round(_mean(dc_cais), 4)
    summary.cai_median_biocompiler = round(_median(bc_cais), 4)
    summary.cai_median_dnachisel = round(_median(dc_cais), 4)

    # CAI wins (higher is better) — include BioCompiler-only as BC wins
    for r in results:
        if r.dnachisel_cai is None:
            summary.cai_biocompiler_wins += 1
        elif r.biocompiler_cai > r.dnachisel_cai + _EPSILON:
            summary.cai_biocompiler_wins += 1
        elif r.dnachisel_cai > r.biocompiler_cai + _EPSILON:
            summary.cai_dnachisel_wins += 1

    # CAI significance
    cai_diffs = [bc - dc for bc, dc in zip(bc_cais, dc_cais)]
    npaired = len(paired)
    if npaired >= 30:
        summary.cai_p_value = round(_paired_t_test_p_value(cai_diffs), 6)
    else:
        summary.cai_p_value = round(_wilcoxon_signed_rank_p_value(cai_diffs), 6)

    # ── GC content comparison ──
    bc_gcs = [r.biocompiler_gc for r in paired]
    dc_gcs = [r.dnachisel_gc for r in paired]  # type: ignore[misc]

    summary.gc_mean_biocompiler = round(_mean(bc_gcs), 4)
    summary.gc_mean_dnachisel = round(_mean(dc_gcs), 4)
    summary.gc_median_biocompiler = round(_median(bc_gcs), 4)
    summary.gc_median_dnachisel = round(_median(dc_gcs), 4)

    # GC wins (closer to midpoint is better)
    for r in results:
        if r.dnachisel_gc is None:
            summary.gc_biocompiler_wins += 1
        else:
            bc_diff = abs(r.biocompiler_gc - _GC_MIDPOINT)
            dc_diff = abs(r.dnachisel_gc - _GC_MIDPOINT)
            if bc_diff < dc_diff - _EPSILON:
                summary.gc_biocompiler_wins += 1
            elif dc_diff < bc_diff - _EPSILON:
                summary.gc_dnachisel_wins += 1

    # GC significance
    gc_diffs = [
        abs(bc - _GC_MIDPOINT) - abs(dc - _GC_MIDPOINT)
        for bc, dc in zip(bc_gcs, dc_gcs)
    ]
    if npaired >= 30:
        summary.gc_p_value = round(_paired_t_test_p_value(gc_diffs), 6)
    else:
        summary.gc_p_value = round(_wilcoxon_signed_rank_p_value(gc_diffs), 6)

    # ── Restriction site comparison ──
    bc_rss = [float(r.biocompiler_restriction_sites) for r in paired]
    dc_rss = [
        float(r.dnachisel_restriction_sites)
        for r in paired
        if r.dnachisel_restriction_sites is not None
    ]

    summary.rs_mean_biocompiler = round(_mean(bc_rss), 2)
    summary.rs_mean_dnachisel = round(_mean(dc_rss), 2)
    summary.rs_median_biocompiler = round(_median(bc_rss), 2)
    summary.rs_median_dnachisel = round(_median(dc_rss), 2)

    # RS wins (fewer is better)
    for r in results:
        if r.dnachisel_restriction_sites is None:
            summary.rs_biocompiler_wins += 1
        elif r.biocompiler_restriction_sites < r.dnachisel_restriction_sites:
            summary.rs_biocompiler_wins += 1
        elif r.dnachisel_restriction_sites < r.biocompiler_restriction_sites:
            summary.rs_dnachisel_wins += 1

    # RS significance
    if dc_rss:
        rs_diffs = [bc - dc for bc, dc in zip(bc_rss, dc_rss)]
        if npaired >= 30:
            summary.rs_p_value = round(_paired_t_test_p_value(rs_diffs), 6)
        else:
            summary.rs_p_value = round(
                _wilcoxon_signed_rank_p_value(rs_diffs), 6
            )

    # ── Runtime comparison ──
    bc_times = [r.biocompiler_runtime for r in paired]
    dc_times = [
        r.dnachisel_runtime
        for r in paired
        if r.dnachisel_runtime is not None
    ]

    summary.runtime_mean_biocompiler = round(_mean(bc_times), 4)
    summary.runtime_mean_dnachisel = round(_mean(dc_times), 4)
    summary.runtime_median_biocompiler = round(_median(bc_times), 4)
    summary.runtime_median_dnachisel = round(_median(dc_times), 4)

    # Runtime wins (faster is better, with 10% noise margin)
    for r in results:
        if r.dnachisel_runtime is None:
            summary.runtime_biocompiler_wins += 1
        elif r.biocompiler_runtime < r.dnachisel_runtime * (
            1.0 - _RUNTIME_MARGIN
        ):
            summary.runtime_biocompiler_wins += 1
        elif r.dnachisel_runtime < r.biocompiler_runtime * (
            1.0 - _RUNTIME_MARGIN
        ):
            summary.runtime_dnachisel_wins += 1

    # Runtime significance
    if dc_times:
        time_diffs = [bc - dc for bc, dc in zip(bc_times, dc_times)]
        if npaired >= 30:
            summary.runtime_p_value = round(
                _paired_t_test_p_value(time_diffs), 6
            )
        else:
            summary.runtime_p_value = round(
                _wilcoxon_signed_rank_p_value(time_diffs), 6
            )

    return summary
