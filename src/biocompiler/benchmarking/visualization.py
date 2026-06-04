"""
BioCompiler Benchmarking — Visualization Module
=================================================

Publication-quality visualization functions for benchmark results
comparing BioCompiler against DNAchisel.

Provides:
  - ``plot_cai_comparison``: Scatter plot of CAI values (biocompiler vs DNAchisel)
  - ``plot_gc_comparison``: Paired bar chart of GC content
  - ``plot_runtime_comparison``: Box plot of runtimes
  - ``plot_constraint_satisfaction``: Stacked bar chart of constraint satisfaction
  - ``plot_summary_dashboard``: 4-panel dashboard combining all plots
  - ``generate_latex_table``: LaTeX table for publication

All plots use Noto Sans SC font for CJK support and degrade gracefully
when DNAchisel results are ``None``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm

try:
    fm.fontManager.addfont("/usr/share/fonts/truetype/chinese/NotoSansSC[wght].ttf")
except Exception:
    logging.getLogger(__name__).debug("NotoSansSC font loading failed", exc_info=True)

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = ["Noto Sans SC", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

if TYPE_CHECKING:
    from .comparison import BenchmarkResult

__all__ = [
    "plot_cai_comparison",
    "plot_gc_comparison",
    "plot_runtime_comparison",
    "plot_constraint_satisfaction",
    "plot_summary_dashboard",
    "generate_latex_table",
]

# ─── Color Palette ────────────────────────────────────────────────────

_COLOR_BIOCOMPILER = "#2563EB"   # Blue
_COLOR_DNACHISEL = "#DC2626"     # Red
_COLOR_DIAGONAL = "#9CA3AF"      # Gray
_COLOR_SATISFIED = "#16A34A"     # Green
_COLOR_VIOLATED = "#EF4444"      # Red
_COLOR_NEUTRAL = "#F59E0B"       # Amber

# ─── Constraint display names (consistent ordering) ──────────────────

_CONSTRAINT_DISPLAY: dict[str, str] = {
    "Translation": "Translation",
    "GC Range": "GC Range",
    "No Restriction Sites": "No RS",
    "No GT Dinucleotides": "No GT",
    "Low CpG Ratio": "Low CpG",
}

# Figure settings
_DPI = 150
_FIGSIZE_SINGLE = (8, 6)
_FIGSIZE_DASHBOARD = (16, 12)


# ─── Helper ───────────────────────────────────────────────────────────

def _ensure_dir(path: str | Path) -> Path:
    """Ensure the parent directory of *path* exists and return a Path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _filter_dnachisel_available(
    results: list[BenchmarkResult],
) -> list[BenchmarkResult]:
    """Return only results that have DNAchisel data."""
    return [r for r in results if r.cai_dnachisel is not None]


# ─── 1. CAI Scatter Plot ─────────────────────────────────────────────

def plot_cai_comparison(
    results: list[BenchmarkResult],
    output_path: str,
) -> None:
    """Scatter plot of CAI values: BioCompiler (x) vs DNAchisel (y).

    Each gene is labeled.  A diagonal reference line shows where
    both tools would have equal CAI.  Points above the line indicate
    DNAchisel is better; points below indicate BioCompiler is better.

    Genes where DNAchisel results are ``None`` are plotted on the
    x-axis at ``y = 0`` with a distinct marker.

    Args:
        results: List of BenchmarkResult objects.
        output_path: File path to save the figure (e.g., ``"cai.png"``).
    """
    out = _ensure_dir(output_path)
    fig, ax = plt.subplots(figsize=_FIGSIZE_SINGLE)

    dc_results = _filter_dnachisel_available(results)
    no_dc = [r for r in results if r.cai_dnachisel is None]

    # Diagonal reference line
    if dc_results:
        all_cai = (
            [r.cai_biocompiler for r in dc_results]
            + [r.cai_dnachisel for r in dc_results]
        )
        lo = min(all_cai) - 0.05
        hi = max(all_cai) + 0.05
    else:
        lo = 0.0
        hi = 1.0
    lo = max(0.0, lo)
    hi = min(1.0, hi)
    ax.plot([lo, hi], [lo, hi], "--", color=_COLOR_DIAGONAL, linewidth=1, label="Equal CAI")

    # Plot paired points
    if dc_results:
        x = [r.cai_biocompiler for r in dc_results]
        y = [r.cai_dnachisel for r in dc_results]
        ax.scatter(x, y, s=80, color=_COLOR_BIOCOMPILER, edgecolors="white",
                   linewidths=0.5, zorder=3, label="BioCompiler vs DNAchisel")
        for r in dc_results:
            ax.annotate(
                r.gene_name,
                (r.cai_biocompiler, r.cai_dnachisel),
                textcoords="offset points",
                xytext=(6, 6),
                fontsize=8,
                alpha=0.85,
            )

    # Plot BioCompiler-only results on x-axis
    if no_dc:
        x_no = [r.cai_biocompiler for r in no_dc]
        y_no = [0.0] * len(no_dc)
        ax.scatter(x_no, y_no, s=60, color=_COLOR_NEUTRAL, marker="D",
                   edgecolors="white", linewidths=0.5, zorder=3,
                   label="BioCompiler only (DNAchisel N/A)")
        for r in no_dc:
            ax.annotate(
                r.gene_name,
                (r.cai_biocompiler, 0.0),
                textcoords="offset points",
                xytext=(6, -12),
                fontsize=8,
                alpha=0.85,
            )

    ax.set_xlabel("BioCompiler CAI", fontsize=12)
    ax.set_ylabel("DNAchisel CAI", fontsize=12)
    ax.set_title("CAI Comparison: BioCompiler vs DNAchisel", fontsize=14)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(out), dpi=_DPI, bbox_inches="tight")
    plt.close(fig)


# ─── 2. GC Content Paired Bar Chart ──────────────────────────────────

def plot_gc_comparison(
    results: list[BenchmarkResult],
    output_path: str,
) -> None:
    """Paired bar chart of GC content for each gene.

    BioCompiler bars are shown in blue; DNAchisel bars (when available)
    in red.  Genes without DNAchisel data show only the BioCompiler bar.

    Args:
        results: List of BenchmarkResult objects.
        output_path: File path to save the figure.
    """
    out = _ensure_dir(output_path)
    fig, ax = plt.subplots(figsize=(max(8, len(results) * 0.8), 6))

    gene_names = [r.gene_name for r in results]
    n = len(gene_names)
    x = np.arange(n)
    bar_width = 0.35

    bc_gc = [r.gc_biocompiler for r in results]
    dc_gc = [r.gc_dnachisel if r.gc_dnachisel is not None else 0.0 for r in results]
    has_dc = [r.gc_dnachisel is not None for r in results]

    bars_bc = ax.bar(x - bar_width / 2, bc_gc, bar_width,
                     color=_COLOR_BIOCOMPILER, label="BioCompiler", edgecolor="white")

    # Only show DNAchisel bars where data exists
    dc_vals = [r.gc_dnachisel if r.gc_dnachisel is not None else np.nan for r in results]
    bars_dc = ax.bar(x + bar_width / 2, dc_vals, bar_width,
                     color=_COLOR_DNACHISEL, label="DNAchisel", edgecolor="white")

    # Add value labels on bars
    for bar in bars_bc:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, height + 0.005,
                    f"{height:.2f}", ha="center", va="bottom", fontsize=7)

    for bar, val in zip(bars_dc, dc_vals):
        if not np.isnan(val) and val > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, val + 0.005,
                    f"{val:.2f}", ha="center", va="bottom", fontsize=7)

    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("GC Content", fontsize=12)
    ax.set_title("GC Content Comparison: BioCompiler vs DNAchisel", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(gene_names, rotation=45, ha="right", fontsize=9)
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(max(bc_gc), max(v for v in dc_vals if not np.isnan(v)) if any(not np.isnan(v) for v in dc_vals) else 0) + 0.1)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(out), dpi=_DPI, bbox_inches="tight")
    plt.close(fig)


# ─── 3. Runtime Box Plot ─────────────────────────────────────────────

def plot_runtime_comparison(
    results: list[BenchmarkResult],
    output_path: str,
) -> None:
    """Box plot of runtimes for BioCompiler and DNAchisel.

    Two boxes are drawn side-by-side.  If no DNAchisel data is
    available, only the BioCompiler box is shown with a note.

    Args:
        results: List of BenchmarkResult objects.
        output_path: File path to save the figure.
    """
    out = _ensure_dir(output_path)
    fig, ax = plt.subplots(figsize=_FIGSIZE_SINGLE)

    bc_runtimes = [r.runtime_biocompiler for r in results]
    dc_results = _filter_dnachisel_available(results)
    dc_runtimes = [r.runtime_dnachisel for r in dc_results if r.runtime_dnachisel is not None]  # type: ignore[union-attr]

    data = [bc_runtimes]
    labels = ["BioCompiler"]
    colors = [_COLOR_BIOCOMPILER]

    if dc_runtimes:
        data.append(dc_runtimes)
        labels.append("DNAchisel")
        colors.append(_COLOR_DNACHISEL)

    bp = ax.boxplot(
        data,
        labels=labels,
        patch_artist=True,
        widths=0.5,
        showmeans=True,
        meanprops=dict(marker="D", markerfacecolor="white", markeredgecolor="black", markersize=6),
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Scatter individual points
    for i, times in enumerate(data, start=1):
        jitter = np.random.default_rng(42).uniform(-0.08, 0.08, size=len(times))
        ax.scatter(np.full(len(times), i) + jitter, times,
                   alpha=0.5, s=20, color=colors[i - 1], edgecolors="white", linewidths=0.3)

    ax.set_ylabel("Runtime (seconds)", fontsize=12)
    ax.set_title("Runtime Comparison: BioCompiler vs DNAchisel", fontsize=14)

    if not dc_runtimes:
        ax.text(0.5, 0.95, "DNAchisel results unavailable",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=10, color="gray", style="italic")

    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(out), dpi=_DPI, bbox_inches="tight")
    plt.close(fig)


# ─── 4. Constraint Satisfaction Stacked Bar Chart ────────────────────

def plot_constraint_satisfaction(
    results: list[BenchmarkResult],
    output_path: str,
) -> None:
    """Stacked bar chart showing percentage of genes satisfying each constraint.

    For each constraint type, two bars are shown (BioCompiler and DNAchisel).
    Each bar is split into satisfied (green) and violated (red) segments.

    If DNAchisel results are ``None``, only BioCompiler bars are shown.

    Args:
        results: List of BenchmarkResult objects.
        output_path: File path to save the figure.
    """
    out = _ensure_dir(output_path)

    # Collect constraint names from the first result that has them
    constraint_names: list[str] = []
    for r in results:
        if r.constraints_biocompiler:
            constraint_names = list(r.constraints_biocompiler.keys())
            break
    if not constraint_names:
        # No constraint data available — write an empty figure
        fig, ax = plt.subplots(figsize=_FIGSIZE_SINGLE)
        ax.text(0.5, 0.5, "No constraint data available",
                transform=ax.transAxes, ha="center", va="center", fontsize=14)
        fig.savefig(str(out), dpi=_DPI, bbox_inches="tight")
        plt.close(fig)
        return

    n_constraints = len(constraint_names)
    fig, ax = plt.subplots(figsize=(max(8, n_constraints * 1.5), 6))

    x = np.arange(n_constraints)
    bar_width = 0.35

    # Compute satisfaction percentages for BioCompiler
    bc_satisfied_pct: list[float] = []
    bc_violated_pct: list[float] = []
    for c in constraint_names:
        total = sum(1 for r in results if c in r.constraints_biocompiler)
        sat = sum(1 for r in results if r.constraints_biocompiler.get(c, False))
        pct = (sat / total * 100) if total > 0 else 0.0
        bc_satisfied_pct.append(pct)
        bc_violated_pct.append(100.0 - pct)

    # BioCompiler stacked bar
    ax.bar(x - bar_width / 2, bc_satisfied_pct, bar_width,
           color=_COLOR_SATISFIED, label="Satisfied (BC)")
    ax.bar(x - bar_width / 2, bc_violated_pct, bar_width,
           bottom=bc_satisfied_pct, color=_COLOR_VIOLATED, label="Violated (BC)")

    # DNAchisel stacked bar (only if data is available)
    dc_results = _filter_dnachisel_available(results)
    if dc_results:
        dc_satisfied_pct: list[float] = []
        dc_violated_pct: list[float] = []
        for c in constraint_names:
            total = sum(1 for r in dc_results if r.constraints_dnachisel and c in r.constraints_dnachisel)
            sat = sum(1 for r in dc_results if r.constraints_dnachisel and r.constraints_dnachisel.get(c, False))
            pct = (sat / total * 100) if total > 0 else 0.0
            dc_satisfied_pct.append(pct)
            dc_violated_pct.append(100.0 - pct)

        ax.bar(x + bar_width / 2, dc_satisfied_pct, bar_width,
               color=_COLOR_SATISFIED, alpha=0.6, label="Satisfied (DC)")
        ax.bar(x + bar_width / 2, dc_violated_pct, bar_width,
               bottom=dc_satisfied_pct, color=_COLOR_VIOLATED, alpha=0.6, label="Violated (DC)")

    # Labels
    display_names = [_CONSTRAINT_DISPLAY.get(c, c) for c in constraint_names]
    ax.set_xlabel("Constraint", fontsize=12)
    ax.set_ylabel("% of Genes", fontsize=12)
    ax.set_title("Constraint Satisfaction Rate", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(display_names, rotation=30, ha="right", fontsize=9)
    ax.set_ylim(0, 105)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(out), dpi=_DPI, bbox_inches="tight")
    plt.close(fig)


# ─── 5. Summary Dashboard ────────────────────────────────────────────

def plot_summary_dashboard(
    results: list[BenchmarkResult],
    output_path: str,
) -> None:
    """4-panel dashboard combining all comparison plots.

    Panels:
      1. Top-left: CAI scatter comparison
      2. Top-right: GC content paired bars
      3. Bottom-left: Runtime box plot
      4. Bottom-right: Constraint satisfaction stacked bars

    Args:
        results: List of BenchmarkResult objects.
        output_path: File path to save the figure.
    """
    out = _ensure_dir(output_path)
    fig, axes = plt.subplots(2, 2, figsize=_FIGSIZE_DASHBOARD)

    # ── Panel 1: CAI Scatter ──
    ax1 = axes[0, 0]
    dc_results = _filter_dnachisel_available(results)
    no_dc = [r for r in results if r.cai_dnachisel is None]

    if dc_results:
        all_cai = (
            [r.cai_biocompiler for r in dc_results]
            + [r.cai_dnachisel for r in dc_results]  # type: ignore[union-attr]
        )
        lo = max(0.0, min(all_cai) - 0.05)
        hi = min(1.0, max(all_cai) + 0.05)
    else:
        lo, hi = 0.0, 1.0
    ax1.plot([lo, hi], [lo, hi], "--", color=_COLOR_DIAGONAL, linewidth=1)

    if dc_results:
        ax1.scatter(
            [r.cai_biocompiler for r in dc_results],
            [r.cai_dnachisel for r in dc_results],  # type: ignore[union-attr]
            s=60, color=_COLOR_BIOCOMPILER, edgecolors="white", linewidths=0.5, zorder=3,
        )
        for r in dc_results:
            ax1.annotate(r.gene_name, (r.cai_biocompiler, r.cai_dnachisel),  # type: ignore[union-attr]
                         textcoords="offset points", xytext=(4, 4), fontsize=7, alpha=0.8)

    if no_dc:
        ax1.scatter(
            [r.cai_biocompiler for r in no_dc], [0.0] * len(no_dc),
            s=40, color=_COLOR_NEUTRAL, marker="D", edgecolors="white",
            linewidths=0.5, zorder=3,
        )
        for r in no_dc:
            ax1.annotate(r.gene_name, (r.cai_biocompiler, 0.0),
                         textcoords="offset points", xytext=(4, -10), fontsize=7, alpha=0.8)

    ax1.set_xlabel("BioCompiler CAI")
    ax1.set_ylabel("DNAchisel CAI")
    ax1.set_title("CAI Comparison")
    ax1.set_xlim(left=0.0)
    ax1.set_ylim(bottom=0.0)
    ax1.grid(True, alpha=0.3)

    # ── Panel 2: GC Paired Bars ──
    ax2 = axes[0, 1]
    gene_names = [r.gene_name for r in results]
    n = len(gene_names)
    x_pos = np.arange(n)
    bar_w = 0.35

    bc_gc = [r.gc_biocompiler for r in results]
    ax2.bar(x_pos - bar_w / 2, bc_gc, bar_w,
            color=_COLOR_BIOCOMPILER, label="BioCompiler", edgecolor="white")

    dc_gc_vals = [r.gc_dnachisel if r.gc_dnachisel is not None else np.nan for r in results]
    if any(not np.isnan(v) for v in dc_gc_vals):
        ax2.bar(x_pos + bar_w / 2, dc_gc_vals, bar_w,
                color=_COLOR_DNACHISEL, label="DNAchisel", edgecolor="white")

    ax2.set_xlabel("Gene")
    ax2.set_ylabel("GC Content")
    ax2.set_title("GC Content Comparison")
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(gene_names, rotation=45, ha="right", fontsize=7)
    ax2.legend(fontsize=7)
    ax2.grid(axis="y", alpha=0.3)

    # ── Panel 3: Runtime Box Plot ──
    ax3 = axes[1, 0]
    bc_runtimes = [r.runtime_biocompiler for r in results]
    dc_runtimes = [r.runtime_dnachisel for r in dc_results if r.runtime_dnachisel is not None]  # type: ignore[union-attr]

    data = [bc_runtimes]
    labels_rt = ["BioCompiler"]
    colors_rt = [_COLOR_BIOCOMPILER]

    if dc_runtimes:
        data.append(dc_runtimes)
        labels_rt.append("DNAchisel")
        colors_rt.append(_COLOR_DNACHISEL)

    bp = ax3.boxplot(data, labels=labels_rt, patch_artist=True, widths=0.5,
                     showmeans=True,
                     meanprops=dict(marker="D", markerfacecolor="white",
                                    markeredgecolor="black", markersize=5))
    for patch, color in zip(bp["boxes"], colors_rt):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax3.set_ylabel("Runtime (s)")
    ax3.set_title("Runtime Comparison")
    ax3.grid(axis="y", alpha=0.3)

    # ── Panel 4: Constraint Satisfaction ──
    ax4 = axes[1, 1]
    constraint_names: list[str] = []
    for r in results:
        if r.constraints_biocompiler:
            constraint_names = list(r.constraints_biocompiler.keys())
            break

    if constraint_names:
        nc = len(constraint_names)
        cx = np.arange(nc)

        bc_sat = []
        for c in constraint_names:
            total = sum(1 for r in results if c in r.constraints_biocompiler)
            sat = sum(1 for r in results if r.constraints_biocompiler.get(c, False))
            bc_sat.append((sat / total * 100) if total > 0 else 0.0)
        bc_viol = [100.0 - s for s in bc_sat]

        ax4.bar(cx - bar_w / 2, bc_sat, bar_w, color=_COLOR_SATISFIED, label="Satisfied (BC)")
        ax4.bar(cx - bar_w / 2, bc_viol, bar_w, bottom=bc_sat, color=_COLOR_VIOLATED, label="Violated (BC)")

        if dc_results:
            dc_sat = []
            for c in constraint_names:
                total = sum(1 for r in dc_results if r.constraints_dnachisel and c in r.constraints_dnachisel)
                sat = sum(1 for r in dc_results if r.constraints_dnachisel and r.constraints_dnachisel.get(c, False))
                dc_sat.append((sat / total * 100) if total > 0 else 0.0)
            dc_viol = [100.0 - s for s in dc_sat]

            ax4.bar(cx + bar_w / 2, dc_sat, bar_w, color=_COLOR_SATISFIED, alpha=0.6, label="Satisfied (DC)")
            ax4.bar(cx + bar_w / 2, dc_viol, bar_w, bottom=dc_sat, color=_COLOR_VIOLATED, alpha=0.6, label="Violated (DC)")

        display_names = [_CONSTRAINT_DISPLAY.get(c, c) for c in constraint_names]
        ax4.set_xticks(cx)
        ax4.set_xticklabels(display_names, rotation=30, ha="right", fontsize=8)
        ax4.set_ylim(0, 105)
    else:
        ax4.text(0.5, 0.5, "No constraint data",
                 transform=ax4.transAxes, ha="center", va="center", fontsize=12)

    ax4.set_xlabel("Constraint")
    ax4.set_ylabel("% of Genes")
    ax4.set_title("Constraint Satisfaction")
    ax4.legend(fontsize=6, loc="lower right")
    ax4.grid(axis="y", alpha=0.3)

    fig.suptitle("BioCompiler vs DNAchisel — Benchmark Dashboard", fontsize=16, y=1.01)
    fig.tight_layout()
    fig.savefig(str(out), dpi=_DPI, bbox_inches="tight")
    plt.close(fig)


# ─── 6. LaTeX Table ──────────────────────────────────────────────────

def generate_latex_table(results: list[BenchmarkResult]) -> str:
    r"""Generate a LaTeX table for publication from benchmark results.

    The table includes columns for gene name, organism, CAI (both tools),
    GC content (both tools), runtime (both tools), and constraint
    satisfaction counts.

    DNAchisel values are shown as "---" when unavailable.

    Args:
        results: List of BenchmarkResult objects.

    Returns:
        LaTeX source string for the table.
    """
    from .comparison import ALL_CONSTRAINTS

    lines: list[str] = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{Benchmark comparison: BioCompiler vs DNAchisel}")
    lines.append(r"  \label{tab:benchmark-comparison}")
    lines.append(r"  \begin{tabular}{llcccccccc}")
    lines.append(r"    \toprule")

    # Header
    lines.append(
        r"    Gene & Organism & "
        r"CAI$_{\text{BC}}$ & CAI$_{\text{DC}}$ & "
        r"GC$_{\text{BC}}$ & GC$_{\text{DC}}$ & "
        r"Time$_{\text{BC}}$ (s) & Time$_{\text{DC}}$ (s) & "
        r"Const$_{\text{BC}}$ & Const$_{\text{DC}}$ \\"
    )
    lines.append(r"    \midrule")

    for r in results:
        # Format CAI
        cai_bc = f"{r.cai_biocompiler:.3f}"
        cai_dc = f"{r.cai_dnachisel:.3f}" if r.cai_dnachisel is not None else "---"

        # Format GC
        gc_bc = f"{r.gc_biocompiler:.3f}"
        gc_dc = f"{r.gc_dnachisel:.3f}" if r.gc_dnachisel is not None else "---"

        # Format runtime
        rt_bc = f"{r.runtime_biocompiler:.3f}"
        rt_dc = f"{r.runtime_dnachisel:.3f}" if r.runtime_dnachisel is not None else "---"

        # Constraint satisfaction counts
        bc_sat = sum(1 for v in r.constraints_biocompiler.values() if v)
        bc_total = len(r.constraints_biocompiler) if r.constraints_biocompiler else 0
        const_bc = f"{bc_sat}/{bc_total}" if bc_total > 0 else "---"

        if r.constraints_dnachisel is not None:
            dc_sat = sum(1 for v in r.constraints_dnachisel.values() if v)
            dc_total = len(r.constraints_dnachisel)
            const_dc = f"{dc_sat}/{dc_total}"
        else:
            const_dc = "---"

        # Escape underscores in organism names for LaTeX
        org_latex = r.organism.replace("_", r"\_")

        lines.append(
            f"    {r.gene_name} & {org_latex} & "
            f"{cai_bc} & {cai_dc} & "
            f"{gc_bc} & {gc_dc} & "
            f"{rt_bc} & {rt_dc} & "
            f"{const_bc} & {const_dc} \\\\"
        )

    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")

    return "\n".join(lines)
