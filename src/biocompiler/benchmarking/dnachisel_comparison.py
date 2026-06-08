"""
Head-to-head comparison between BioCompiler and DNAchisel.

Provides a fair, reproducible comparison across all specification
types supported by both tools. Results can be formatted as LaTeX
tables or Markdown reports suitable for publication.

ComparisonResult encapsulates the outcome of a single comparison
dimension (e.g., CAI quality, speed, constraint satisfaction),
providing structured data for both tools along with the determined
winner and a human-readable details string.

Usage::

    from biocompiler.benchmarking.dnachisel_comparison import (
        ComparisonResult,
        compare_on_cai,
        compare_on_speed,
        compare_on_constraints,
        compare_on_eukaryotic_features,
        generate_comparison_table,
        generate_comparison_report,
    )

    genes = {"GFP": "MSKGEELFTGVVPILVELDGDVNG...", "HBB": "MVHLTPEEK..."}
    organisms = ["Escherichia_coli", "Homo_sapiens"]

    cai_result = compare_on_cai(genes, organisms)
    speed_result = compare_on_speed(genes, organisms, repeats=10)
    constraint_result = compare_on_constraints(genes, enzymes=["EcoRI", "BamHI"])
    euk_result = compare_on_eukaryotic_features(genes)

    latex = generate_comparison_table([cai_result, speed_result])
    report = generate_comparison_report([cai_result, speed_result, constraint_result, euk_result])
"""

from __future__ import annotations

import logging
import math
import statistics
import time
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "ComparisonResult",
    "compare_on_cai",
    "compare_on_speed",
    "compare_on_constraints",
    "compare_on_eukaryotic_features",
    "generate_comparison_table",
    "generate_comparison_report",
]


# ---------------------------------------------------------------------------
# ComparisonResult — structured comparison outcome
# ---------------------------------------------------------------------------


@dataclass
class ComparisonResult:
    """Structured result from a head-to-head comparison.

    Attributes:
        category: Comparison dimension (e.g., "CAI", "Speed",
            "Constraint Satisfaction", "Eukaryotic Features").
        biocompiler_metrics: Dict of metric values for BioCompiler
            (e.g., ``{"mean": 0.85, "std": 0.02, "wins": 5}``).
        dnachisel_metrics: Dict of metric values for DNAchisel
            (e.g., ``{"mean": 0.78, "std": 0.03, "wins": 2}``).
            Empty dict if DNAchisel was unavailable.
        winner: One of ``"biocompiler"``, ``"dnachisel"``, or ``"tie"``.
        details: Human-readable summary string.
    """

    category: str
    biocompiler_metrics: dict = field(default_factory=dict)
    dnachisel_metrics: dict = field(default_factory=dict)
    winner: str = "tie"
    details: str = ""


# ---------------------------------------------------------------------------
# Internal helper: resolve DNAchisel adapter
# ---------------------------------------------------------------------------


def _get_adapter():
    """Get a DNAchiselAdapter instance if available, else None."""
    try:
        from .dnachisel_adapter import DNAchiselAdapter, is_dnachisel_available
        if is_dnachisel_available():
            return DNAchiselAdapter()
    except ImportError:
        pass
    return None


# ---------------------------------------------------------------------------
# Internal helper: optimize with BioCompiler
# ---------------------------------------------------------------------------


def _optimize_biocompiler(protein: str, organism: str) -> str:
    """Optimize a protein with BioCompiler and return the DNA sequence.

    Returns empty string on failure.
    """
    try:
        from ..optimization import optimize_sequence
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=0.30,
            gc_hi=0.70,
        )
        return result.sequence
    except Exception as exc:
        logger.debug("BioCompiler optimization failed: %s", exc)
        return ""


# ---------------------------------------------------------------------------
# CAI comparison
# ---------------------------------------------------------------------------


def compare_on_cai(
    genes: dict[str, str],
    organisms: list[str],
) -> ComparisonResult:
    """Compare BioCompiler vs DNAchisel on CAI quality.

    Optimizes each gene for each organism with both tools and
    computes CAI using BioCompiler's validated evaluator for fairness.

    Args:
        genes: Dict mapping gene name to protein sequence.
        organisms: List of target organism names.

    Returns:
        ComparisonResult with CAI comparison data.
    """
    from .metrics import compute_cai_validated

    adapter = _get_adapter()
    epsilon = 0.001

    bc_cais: list[float] = []
    dc_cais: list[float] = []
    bc_wins = 0
    dc_wins = 0
    ties = 0
    per_gene_details: list[str] = []

    for gene_name, protein in genes.items():
        for organism in organisms:
            bc_seq = _optimize_biocompiler(protein, organism)
            bc_cai = compute_cai_validated(bc_seq, organism) if bc_seq else 0.0

            dc_cai = 0.0
            if adapter is not None:
                try:
                    dc_result = adapter.optimize(
                        protein=protein,
                        organism=organism,
                        constraints=[
                            {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                        ],
                    )
                    if dc_result.success:
                        dc_cai = compute_cai_validated(
                            dc_result.sequence, organism
                        )
                except Exception as exc:
                    logger.debug(
                        "DNAchisel CAI failed for %s/%s: %s",
                        gene_name, organism, exc,
                    )

            bc_cais.append(bc_cai)
            if dc_cai > 0:
                dc_cais.append(dc_cai)

            # Determine per-comparison winner
            if dc_cai > 0:
                if bc_cai > dc_cai + epsilon:
                    bc_wins += 1
                elif dc_cai > bc_cai + epsilon:
                    dc_wins += 1
                else:
                    ties += 1
                per_gene_details.append(
                    f"{gene_name}/{organism}: BC={bc_cai:.4f} DC={dc_cai:.4f}"
                )
            else:
                per_gene_details.append(
                    f"{gene_name}/{organism}: BC={bc_cai:.4f} DC=N/A"
                )

    bc_mean = statistics.mean(bc_cais) if bc_cais else 0.0
    bc_std = statistics.stdev(bc_cais) if len(bc_cais) > 1 else 0.0
    dc_mean = statistics.mean(dc_cais) if dc_cais else 0.0
    dc_std = statistics.stdev(dc_cais) if len(dc_cais) > 1 else 0.0

    # Determine overall winner
    if bc_mean > dc_mean + epsilon:
        winner = "biocompiler"
    elif dc_mean > bc_mean + epsilon:
        winner = "dnachisel"
    else:
        winner = "tie"

    delta = bc_mean - dc_mean

    details = (
        f"CAI comparison across {len(bc_cais)} gene/organism pairs. "
        f"BioCompiler mean={bc_mean:.4f}±{bc_std:.4f}, "
        f"DNAchisel mean={dc_mean:.4f}±{dc_std:.4f}, "
        f"Δ={delta:+.4f}. "
        f"Wins: BC={bc_wins}, DC={dc_wins}, Ties={ties}."
    )

    return ComparisonResult(
        category="CAI",
        biocompiler_metrics={
            "mean": round(bc_mean, 4),
            "std": round(bc_std, 4),
            "wins": bc_wins,
            "n": len(bc_cais),
        },
        dnachisel_metrics={
            "mean": round(dc_mean, 4),
            "std": round(dc_std, 4),
            "wins": dc_wins,
            "n": len(dc_cais),
        },
        winner=winner,
        details=details,
    )


# ---------------------------------------------------------------------------
# Speed comparison
# ---------------------------------------------------------------------------


def compare_on_speed(
    genes: dict[str, str],
    organisms: list[str],
    repeats: int = 10,
) -> ComparisonResult:
    """Compare BioCompiler vs DNAchisel on execution speed.

    Runs each optimization multiple times (after warmup) and reports
    median time per gene.

    Args:
        genes: Dict mapping gene name to protein sequence.
        organisms: List of target organism names.
        repeats: Number of measured iterations per gene/organism (default 10).

    Returns:
        ComparisonResult with speed comparison data.
    """
    adapter = _get_adapter()
    warmup = 1
    epsilon = 0.10  # 10% noise margin

    bc_times: list[float] = []
    dc_times: list[float] = []
    bc_wins = 0
    dc_wins = 0
    ties = 0

    for gene_name, protein in genes.items():
        for organism in organisms:
            # BioCompiler speed
            for _ in range(warmup):
                try:
                    from ..optimization import optimize_sequence
                    optimize_sequence(
                        target_protein=protein,
                        organism=organism,
                        gc_lo=0.30,
                        gc_hi=0.70,
                    )
                except Exception:
                    pass

            bc_run_times: list[float] = []
            for _ in range(repeats):
                t0 = time.perf_counter()
                try:
                    from ..optimization import optimize_sequence
                    optimize_sequence(
                        target_protein=protein,
                        organism=organism,
                        gc_lo=0.30,
                        gc_hi=0.70,
                    )
                except Exception:
                    pass
                bc_run_times.append((time.perf_counter() - t0) * 1000.0)

            bc_median = statistics.median(bc_run_times) if bc_run_times else 0.0

            # DNAchisel speed
            dc_median = 0.0
            if adapter is not None:
                for _ in range(warmup):
                    try:
                        adapter.optimize(
                            protein=protein,
                            organism=organism,
                            constraints=[
                                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                            ],
                        )
                    except Exception:
                        pass

                dc_run_times: list[float] = []
                for _ in range(repeats):
                    t0 = time.perf_counter()
                    try:
                        adapter.optimize(
                            protein=protein,
                            organism=organism,
                            constraints=[
                                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                            ],
                        )
                    except Exception:
                        pass
                    dc_run_times.append((time.perf_counter() - t0) * 1000.0)

                dc_median = statistics.median(dc_run_times) if dc_run_times else 0.0

            bc_times.append(bc_median)
            if dc_median > 0:
                dc_times.append(dc_median)

            # Determine winner (lower is better)
            if dc_median > 0:
                if bc_median < dc_median * (1.0 - epsilon):
                    bc_wins += 1
                elif dc_median < bc_median * (1.0 - epsilon):
                    dc_wins += 1
                else:
                    ties += 1

    bc_mean = statistics.mean(bc_times) if bc_times else 0.0
    bc_std = statistics.stdev(bc_times) if len(bc_times) > 1 else 0.0
    dc_mean = statistics.mean(dc_times) if dc_times else 0.0
    dc_std = statistics.stdev(dc_times) if len(dc_times) > 1 else 0.0

    # Determine overall winner (lower is better for speed)
    if dc_mean > 0 and bc_mean < dc_mean * (1.0 - epsilon):
        winner = "biocompiler"
    elif dc_mean > 0 and dc_mean < bc_mean * (1.0 - epsilon):
        winner = "dnachisel"
    else:
        winner = "tie"

    speed_ratio = bc_mean / dc_mean if dc_mean > 0 else 0.0

    details = (
        f"Speed comparison (median of {repeats} runs). "
        f"BioCompiler mean={bc_mean:.2f}±{bc_std:.2f}ms, "
        f"DNAchisel mean={dc_mean:.2f}±{dc_std:.2f}ms, "
        f"ratio={speed_ratio:.2f}x. "
        f"Wins: BC={bc_wins}, DC={dc_wins}, Ties={ties}."
    )

    return ComparisonResult(
        category="Speed",
        biocompiler_metrics={
            "mean_ms": round(bc_mean, 2),
            "std_ms": round(bc_std, 2),
            "wins": bc_wins,
            "n": len(bc_times),
        },
        dnachisel_metrics={
            "mean_ms": round(dc_mean, 2),
            "std_ms": round(dc_std, 2),
            "wins": dc_wins,
            "n": len(dc_times),
        },
        winner=winner,
        details=details,
    )


# ---------------------------------------------------------------------------
# Constraint satisfaction comparison
# ---------------------------------------------------------------------------


def compare_on_constraints(
    genes: dict[str, str],
    enzymes: list[str],
) -> ComparisonResult:
    """Compare BioCompiler vs DNAchisel on constraint satisfaction.

    For each gene, optimizes the sequence and checks whether
    translation fidelity, GC range compliance, and restriction
    site avoidance constraints are satisfied.

    Args:
        genes: Dict mapping gene name to protein sequence.
        enzymes: List of restriction enzyme names to check for.

    Returns:
        ComparisonResult with constraint satisfaction data.
    """
    adapter = _get_adapter()
    organisms = ["Escherichia_coli", "Homo_sapiens"]

    bc_rates: list[float] = []
    dc_rates: list[float] = []
    bc_wins = 0
    dc_wins = 0
    ties = 0

    for gene_name, protein in genes.items():
        for organism in organisms:
            # BioCompiler constraints
            bc_satisfied = 0
            bc_total = 0
            bc_seq = _optimize_biocompiler(protein, organism)
            if bc_seq:
                # Translation
                bc_total += 1
                try:
                    from ..translation import translate
                    if translate(bc_seq, to_stop=True) == protein:
                        bc_satisfied += 1
                except Exception:
                    pass
                # GC
                bc_total += 1
                try:
                    from ..scanner import gc_content
                    gc = gc_content(bc_seq)
                    if 0.30 <= gc <= 0.70:
                        bc_satisfied += 1
                except Exception:
                    pass
                # Restriction sites
                bc_total += 1
                try:
                    from .metrics import count_restriction_sites
                    rs = count_restriction_sites(bc_seq, enzymes)
                    if sum(rs.values()) == 0:
                        bc_satisfied += 1
                except Exception:
                    pass

            # DNAchisel constraints
            dc_satisfied = 0
            dc_total = 0
            if adapter is not None:
                try:
                    dc_result = adapter.optimize(
                        protein=protein,
                        organism=organism,
                        constraints=[
                            {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                            {"type": "avoid_restriction", "enzymes": enzymes},
                        ],
                    )
                    if dc_result.success:
                        dc_seq = dc_result.sequence
                        # Translation
                        dc_total += 1
                        try:
                            from ..translation import translate
                            if translate(dc_seq, to_stop=True) == protein:
                                dc_satisfied += 1
                        except Exception:
                            pass
                        # GC
                        dc_total += 1
                        try:
                            from ..scanner import gc_content
                            gc = gc_content(dc_seq)
                            if 0.30 <= gc <= 0.70:
                                dc_satisfied += 1
                        except Exception:
                            pass
                        # Restriction sites
                        dc_total += 1
                        try:
                            from .metrics import count_restriction_sites
                            rs = count_restriction_sites(dc_seq, enzymes)
                            if sum(rs.values()) == 0:
                                dc_satisfied += 1
                        except Exception:
                            pass
                except Exception:
                    pass

            bc_rate = bc_satisfied / bc_total if bc_total > 0 else 0.0
            dc_rate = dc_satisfied / dc_total if dc_total > 0 else 0.0

            bc_rates.append(bc_rate)
            if dc_total > 0:
                dc_rates.append(dc_rate)

                if bc_rate > dc_rate + 0.001:
                    bc_wins += 1
                elif dc_rate > bc_rate + 0.001:
                    dc_wins += 1
                else:
                    ties += 1

    bc_mean = statistics.mean(bc_rates) if bc_rates else 0.0
    dc_mean = statistics.mean(dc_rates) if dc_rates else 0.0

    if bc_mean > dc_mean + 0.001:
        winner = "biocompiler"
    elif dc_mean > bc_mean + 0.001:
        winner = "dnachisel"
    else:
        winner = "tie"

    details = (
        f"Constraint satisfaction comparison (enzymes: {', '.join(enzymes)}). "
        f"BioCompiler mean={bc_mean:.1%}, DNAchisel mean={dc_mean:.1%}. "
        f"Wins: BC={bc_wins}, DC={dc_wins}, Ties={ties}."
    )

    return ComparisonResult(
        category="Constraint Satisfaction",
        biocompiler_metrics={
            "mean_rate": round(bc_mean, 4),
            "wins": bc_wins,
            "n": len(bc_rates),
        },
        dnachisel_metrics={
            "mean_rate": round(dc_mean, 4),
            "wins": dc_wins,
            "n": len(dc_rates),
        },
        winner=winner,
        details=details,
    )


# ---------------------------------------------------------------------------
# Eukaryotic feature comparison
# ---------------------------------------------------------------------------


def compare_on_eukaryotic_features(
    genes: dict[str, str],
) -> ComparisonResult:
    """Compare BioCompiler vs DNAchisel on eukaryotic feature handling.

    Tests CpG avoidance, splice site avoidance, and other eukaryotic-
    specific constraints that BioCompiler supports natively.

    Args:
        genes: Dict mapping gene name to protein sequence.

    Returns:
        ComparisonResult with eukaryotic feature comparison data.
    """
    adapter = _get_adapter()
    organisms = ["Homo_sapiens", "CHO_K1"]

    bc_rates: list[float] = []
    dc_rates: list[float] = []
    bc_wins = 0
    dc_wins = 0
    ties = 0

    for gene_name, protein in genes.items():
        for organism in organisms:
            bc_features_satisfied = 0
            bc_total_features = 3  # CpG, splice sites, GC compliance

            try:
                bc_seq = _optimize_biocompiler(protein, organism)
                if bc_seq:
                    # Check GC compliance
                    from ..scanner import gc_content
                    gc = gc_content(bc_seq)
                    if 0.30 <= gc <= 0.70:
                        bc_features_satisfied += 1

                    # Check CpG ratio
                    seq_upper = bc_seq.upper()
                    c_count = seq_upper.count("C")
                    g_count = seq_upper.count("G")
                    cpg_count = sum(
                        1 for i in range(len(seq_upper) - 1)
                        if seq_upper[i:i + 2] == "CG"
                    )
                    if c_count > 0 and g_count > 0:
                        cpg_ratio = (cpg_count * len(seq_upper)) / (c_count * g_count)
                        if cpg_ratio < 0.60:
                            bc_features_satisfied += 1

                    # Check cryptic splice sites (GT/AG in codons)
                    gt_count = sum(
                        1 for i in range(0, len(seq_upper) - 2, 3)
                        for j in range(2)
                        if seq_upper[i + j:i + j + 2] in ("GT", "AG")
                    )
                    expected_gt = len(protein) * 0.3
                    if gt_count <= expected_gt:
                        bc_features_satisfied += 1

            except Exception as exc:
                logger.debug(
                    "BioCompiler eukaryotic check failed for %s: %s",
                    gene_name, exc,
                )

            # DNAchisel: check same features
            dc_features_satisfied = 0
            dc_total_features = 3

            if adapter is not None:
                try:
                    dc_result = adapter.optimize(
                        protein=protein,
                        organism=organism,
                        constraints=[
                            {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                        ],
                    )
                    if dc_result.success:
                        dc_seq = dc_result.sequence
                        # GC compliance
                        from ..scanner import gc_content
                        gc = gc_content(dc_seq)
                        if 0.30 <= gc <= 0.70:
                            dc_features_satisfied += 1

                        # CpG ratio
                        seq_upper = dc_seq.upper()
                        c_count = seq_upper.count("C")
                        g_count = seq_upper.count("G")
                        cpg_count = sum(
                            1 for i in range(len(seq_upper) - 1)
                            if seq_upper[i:i + 2] == "CG"
                        )
                        if c_count > 0 and g_count > 0:
                            cpg_ratio = (cpg_count * len(seq_upper)) / (c_count * g_count)
                        else:
                            cpg_ratio = 0.0
                        if cpg_ratio < 0.60:
                            dc_features_satisfied += 1

                        # Splice sites
                        gt_count = sum(
                            1 for i in range(0, len(seq_upper) - 2, 3)
                            for j in range(2)
                            if seq_upper[i + j:i + j + 2] in ("GT", "AG")
                        )
                        expected_gt = len(protein) * 0.3
                        if gt_count <= expected_gt:
                            dc_features_satisfied += 1

                except Exception as exc:
                    logger.debug(
                        "DNAchisel eukaryotic check failed for %s: %s",
                        gene_name, exc,
                    )

            bc_rate = bc_features_satisfied / bc_total_features if bc_total_features > 0 else 0.0
            dc_rate = dc_features_satisfied / dc_total_features if dc_total_features > 0 else 0.0

            bc_rates.append(bc_rate)
            dc_rates.append(dc_rate)

            if bc_rate > dc_rate + 0.001:
                bc_wins += 1
            elif dc_rate > bc_rate + 0.001:
                dc_wins += 1
            else:
                ties += 1

    bc_mean = statistics.mean(bc_rates) if bc_rates else 0.0
    dc_mean = statistics.mean(dc_rates) if dc_rates else 0.0

    if bc_mean > dc_mean + 0.001:
        winner = "biocompiler"
    elif dc_mean > bc_mean + 0.001:
        winner = "dnachisel"
    else:
        winner = "tie"

    details = (
        f"Eukaryotic feature comparison (CpG, splice sites, GC). "
        f"BioCompiler mean={bc_mean:.1%}, DNAchisel mean={dc_mean:.1%}. "
        f"Wins: BC={bc_wins}, DC={dc_wins}, Ties={ties}."
    )

    return ComparisonResult(
        category="Eukaryotic Features",
        biocompiler_metrics={
            "mean_rate": round(bc_mean, 4),
            "wins": bc_wins,
            "n": len(bc_rates),
        },
        dnachisel_metrics={
            "mean_rate": round(dc_mean, 4),
            "wins": dc_wins,
            "n": len(dc_rates),
        },
        winner=winner,
        details=details,
    )


# ---------------------------------------------------------------------------
# LaTeX table generation
# ---------------------------------------------------------------------------


def generate_comparison_table(results: list[ComparisonResult]) -> str:
    r"""Generate a LaTeX table from a list of ComparisonResult objects.

    Args:
        results: List of ComparisonResult from comparison functions.

    Returns:
        LaTeX source string for the comparison table.
    """
    lines: list[str] = []

    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(
        r"  \caption{Head-to-head comparison: BioCompiler vs DNAchisel}"
    )
    lines.append(r"  \label{tab:head-to-head}")
    lines.append(r"  \begin{tabular}{lrrrrr}")
    lines.append(r"    \toprule")

    lines.append(
        r"    Category & BC Mean & DC Mean & Winner & BC Wins & DC Wins \\"
    )
    lines.append(r"    \midrule")

    for r in results:
        bc_metrics = r.biocompiler_metrics
        dc_metrics = r.dnachisel_metrics

        # Get mean value from whichever key is present
        bc_mean = bc_metrics.get("mean", bc_metrics.get("mean_ms", bc_metrics.get("mean_rate", 0.0)))
        dc_mean = dc_metrics.get("mean", dc_metrics.get("mean_ms", dc_metrics.get("mean_rate", 0.0)))

        # Determine unit suffix for display
        if "mean_ms" in bc_metrics:
            bc_str = f"{bc_mean:.2f} ms"
            dc_str = f"{dc_mean:.2f} ms" if dc_mean > 0 else "---"
        elif "mean_rate" in bc_metrics:
            bc_str = f"{bc_mean:.1%}"
            dc_str = f"{dc_mean:.1%}" if dc_mean > 0 else "---"
        else:
            bc_str = f"{bc_mean:.4f}"
            dc_str = f"{dc_mean:.4f}" if dc_mean > 0 else "---"

        bc_wins = bc_metrics.get("wins", 0)
        dc_wins = dc_metrics.get("wins", 0)

        lines.append(
            f"    {r.category} & {bc_str} & {dc_str} & "
            f"{r.winner} & {bc_wins} & {dc_wins} \\\\"
        )

    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Markdown report generation
# ---------------------------------------------------------------------------


def generate_comparison_report(results: list[ComparisonResult]) -> str:
    """Generate a Markdown report from a list of ComparisonResult objects.

    Args:
        results: List of ComparisonResult from comparison functions.

    Returns:
        Markdown string with the comparison report.
    """
    lines: list[str] = []

    lines.append("# BioCompiler v12 vs DNAchisel — Head-to-Head Comparison")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("| Category | BioCompiler | DNAchisel | Winner |")
    lines.append("|----------|-------------|-----------|--------|")

    for r in results:
        bc_metrics = r.biocompiler_metrics
        dc_metrics = r.dnachisel_metrics

        bc_mean = bc_metrics.get("mean", bc_metrics.get("mean_ms", bc_metrics.get("mean_rate", 0.0)))
        dc_mean = dc_metrics.get("mean", dc_metrics.get("mean_ms", dc_metrics.get("mean_rate", 0.0)))

        # Format display value
        if "mean_ms" in bc_metrics:
            bc_display = f"{bc_mean:.2f} ms"
            dc_display = f"{dc_mean:.2f} ms" if dc_mean > 0 else "N/A"
        elif "mean_rate" in bc_metrics:
            bc_display = f"{bc_mean:.1%}"
            dc_display = f"{dc_mean:.1%}" if dc_mean > 0 else "N/A"
        else:
            bc_display = f"{bc_mean:.4f}"
            dc_display = f"{dc_mean:.4f}" if dc_mean > 0 else "N/A"

        lines.append(
            f"| {r.category} | {bc_display} | {dc_display} | {r.winner} |"
        )

    lines.append("")
    lines.append("## Details")
    lines.append("")

    for r in results:
        lines.append(f"### {r.category}")
        lines.append("")
        lines.append(r.details)
        lines.append("")

        # Per-metric breakdown
        lines.append("**BioCompiler metrics:**")
        for key, val in r.biocompiler_metrics.items():
            lines.append(f"- {key}: {val}")
        lines.append("")

        if r.dnachisel_metrics:
            lines.append("**DNAchisel metrics:**")
            for key, val in r.dnachisel_metrics.items():
                lines.append(f"- {key}: {val}")
        else:
            lines.append("**DNAchisel metrics:** Not available")
        lines.append("")

    # Overall verdict
    bc_overall_wins = sum(
        1 for r in results if r.winner == "biocompiler"
    )
    dc_overall_wins = sum(
        1 for r in results if r.winner == "dnachisel"
    )
    overall_ties = sum(
        1 for r in results if r.winner == "tie"
    )

    lines.append("## Overall Verdict")
    lines.append("")
    if bc_overall_wins > dc_overall_wins:
        lines.append(
            f"**BioCompiler wins** in {bc_overall_wins}/{len(results)} "
            f"categories (DNAchisel: {dc_overall_wins}, Ties: {overall_ties})."
        )
    elif dc_overall_wins > bc_overall_wins:
        lines.append(
            f"**DNAchisel wins** in {dc_overall_wins}/{len(results)} "
            f"categories (BioCompiler: {bc_overall_wins}, Ties: {overall_ties})."
        )
    else:
        lines.append(
            f"**Tie** — each tool wins {bc_overall_wins} categories "
            f"with {overall_ties} ties."
        )
    lines.append("")

    return "\n".join(lines)
