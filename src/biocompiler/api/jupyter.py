"""
BioCompiler Jupyter Notebook Integration

Provides Jupyter notebook-friendly visualization and interaction.
Since scientists use Jupyter for everything, BioCompiler needs
first-class Jupyter support.

IPython, matplotlib, and ipywidgets are OPTIONAL. All functions that
depend on them raise ImportError with helpful messages if not installed.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "display_sequence",
    "display_optimization_result",
    "display_type_check",
    "plot_gc_content",
    "plot_codon_usage",
    "interactive_optimize",
]

# ==============================================================================
# Color Constants for HTML Display
# ==============================================================================

_VERDICT_COLORS: dict[str, str] = {
    "PASS": "#28a745",
    "FAIL": "#dc3545",
    "UNCERTAIN": "#ffc107",
}

_BASE_COLORS: dict[str, str] = {
    "A": "#4CAF50",  # Green
    "T": "#F44336",  # Red
    "G": "#FF9800",  # Orange
    "C": "#2196F3",  # Blue
}

# GC content display thresholds
_GC_LOW_THRESHOLD: float = 0.30
_GC_HIGH_THRESHOLD: float = 0.70
_CAI_DISPLAY_THRESHOLD: float = 0.5


def _check_ipython() -> None:
    """Check that IPython is installed, raise ImportError with helpful message if not."""
    try:
        import IPython  # noqa: F401
    except ImportError:
        raise ImportError(
            "IPython is required for BioCompiler Jupyter display functions but is not installed. "
            "Install it with: pip install ipython>=8.0  "
            "or: pip install biocompiler[jupyter]"
        )


def _check_matplotlib() -> None:
    """Check that matplotlib is installed, raise ImportError with helpful message if not."""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        raise ImportError(
            "matplotlib is required for BioCompiler Jupyter plotting functions but is not installed. "
            "Install it with: pip install matplotlib>=3.5  "
            "or: pip install biocompiler[jupyter]"
        )


def _check_ipywidgets() -> None:
    """Check that ipywidgets is installed, raise ImportError with helpful message if not."""
    try:
        import ipywidgets  # noqa: F401
    except ImportError:
        raise ImportError(
            "ipywidgets is required for BioCompiler interactive widgets but is not installed. "
            "Install it with: pip install ipywidgets>=8.0  "
            "or: pip install biocompiler[jupyter]"
        )


# ==============================================================================
# 1. display_sequence — Rich HTML Sequence Display
# ==============================================================================

def display_sequence(
    sequence: str,
    organism: str = "Homo_sapiens",
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
) -> None:
    """
    Display a DNA sequence as rich HTML in a Jupyter notebook.

    Shows the sequence with:
    - Color-coded nucleotides (A=green, T=red, G=orange, C=blue)
    - Exon regions highlighted with a background shade
    - GC content statistics
    - Codon grouping (every 10 bases)

    Args:
        sequence: DNA sequence string
        organism: Organism name for display
        exon_boundaries: Optional list of (start, end) tuples for exon highlighting

    Raises:
        ImportError: If IPython is not installed
    """
    _check_ipython()

    from IPython.display import HTML, display
    from biocompiler.sequence.scanner import gc_content as _gc_content, validate_dna_sequence

    seq = validate_dna_sequence(sequence)
    gc = _gc_content(seq)
    seq_len = len(seq)

    # Build exon position set for highlighting
    exon_positions = set()
    if exon_boundaries:
        for start, end in exon_boundaries:
            for pos in range(start, min(end, seq_len)):
                exon_positions.add(pos)

    # Color-coded nucleotide HTML
    nucleotide_spans = []
    for i, base in enumerate(seq):
        color = _BASE_COLORS.get(base, "#666666")
        bg_style = ""
        if i in exon_positions:
            bg_style = "background-color: #e8f5e9;"
        nucleotide_spans.append(
            f'<span style="color:{color};font-family:monospace;font-size:14px;{bg_style}">{base}</span>'
        )

    # Group into codons (every 3) with spacing every 10 bases
    formatted_parts = []
    for i, span in enumerate(nucleotide_spans):
        if i > 0 and i % 10 == 0:
            formatted_parts.append(" ")
        formatted_parts.append(span)

    nucleotide_html = "".join(formatted_parts)

    # Stats row
    a_count = seq.count("A")
    t_count = seq.count("T")
    g_count = seq.count("G")
    c_count = seq.count("C")

    stats_html = f"""
    <div style="display:flex;gap:20px;margin-top:8px;font-family:sans-serif;font-size:13px;">
        <span><b>Length:</b> {seq_len} bp</span>
        <span><b>GC:</b> {gc:.1%}</span>
        <span style="color:{_BASE_COLORS['A']}"><b>A:</b> {a_count}</span>
        <span style="color:{_BASE_COLORS['T']}"><b>T:</b> {t_count}</span>
        <span style="color:{_BASE_COLORS['G']}"><b>G:</b> {g_count}</span>
        <span style="color:{_BASE_COLORS['C']}"><b>C:</b> {c_count}</span>
        <span><b>Organism:</b> {organism}</span>
    </div>
    """

    # Exon info
    exon_info = ""
    if exon_boundaries:
        exon_strs = [f"[{s}, {e})" for s, e in exon_boundaries]
        exon_info = f'<div style="margin-top:4px;font-family:sans-serif;font-size:12px;color:#555;"><b>Exons:</b> {", ".join(exon_strs)}</div>'

    html = f"""
    <div style="border:1px solid #ddd;border-radius:8px;padding:12px;background:#fafafa;margin:8px 0;">
        <div style="font-weight:bold;font-size:15px;margin-bottom:6px;font-family:sans-serif;">
            [DNA] Sequence View
        </div>
        <div style="white-space:pre-wrap;word-break:break-all;line-height:1.6;">
            {nucleotide_html}
        </div>
        {stats_html}
        {exon_info}
    </div>
    """

    display(HTML(html))


# ==============================================================================
# 2. display_optimization_result — HTML Table with Color Coding
# ==============================================================================

def display_optimization_result(result) -> None:
    """
    Display an optimization result as a styled HTML table in Jupyter.

    Shows:
    - Optimization metrics (CAI, GC, length) in a card layout
    - Satisfied predicates in green
    - Failed predicates in red
    - Warnings if greedy fallback was used

    Args:
        result: OptimizationResult object from optimize_sequence()

    Raises:
        ImportError: If IPython is not installed
    """
    _check_ipython()

    from IPython.display import HTML, display

    # Metric cards
    cai = result.cai
    gc = result.gc_content
    seq_len = len(result.sequence)

    cai_color = _VERDICT_COLORS["PASS"] if cai >= _CAI_DISPLAY_THRESHOLD else _VERDICT_COLORS["FAIL"]
    gc_color = _VERDICT_COLORS["PASS"] if _GC_LOW_THRESHOLD <= gc <= _GC_HIGH_THRESHOLD else _VERDICT_COLORS["FAIL"]

    metrics_html = f"""
    <div style="display:flex;gap:16px;margin-bottom:12px;">
        <div style="border:1px solid #ddd;border-radius:8px;padding:12px 20px;background:white;text-align:center;min-width:120px;">
            <div style="font-size:12px;color:#888;font-family:sans-serif;">CAI</div>
            <div style="font-size:24px;font-weight:bold;color:{cai_color};">{cai:.4f}</div>
        </div>
        <div style="border:1px solid #ddd;border-radius:8px;padding:12px 20px;background:white;text-align:center;min-width:120px;">
            <div style="font-size:12px;color:#888;font-family:sans-serif;">GC Content</div>
            <div style="font-size:24px;font-weight:bold;color:{gc_color};">{gc:.1%}</div>
        </div>
        <div style="border:1px solid #ddd;border-radius:8px;padding:12px 20px;background:white;text-align:center;min-width:120px;">
            <div style="font-size:12px;color:#888;font-family:sans-serif;">Length</div>
            <div style="font-size:24px;font-weight:bold;color:#333;">{seq_len} bp</div>
        </div>
    </div>
    """

    # Predicate table
    satisfied = result.satisfied_predicates or []
    failed = result.failed_predicates or []

    predicate_rows = ""
    for pred in satisfied:
        predicate_rows += f"""
        <tr>
            <td style="padding:6px 12px;border-bottom:1px solid #eee;font-family:monospace;font-size:13px;">{pred}</td>
            <td style="padding:6px 12px;border-bottom:1px solid #eee;">
                <span style="background-color:{_VERDICT_COLORS['PASS']};color:white;padding:2px 8px;border-radius:4px;font-size:12px;font-weight:bold;">PASS</span>
            </td>
        </tr>
        """
    for pred in failed:
        predicate_rows += f"""
        <tr>
            <td style="padding:6px 12px;border-bottom:1px solid #eee;font-family:monospace;font-size:13px;">{pred}</td>
            <td style="padding:6px 12px;border-bottom:1px solid #eee;">
                <span style="background-color:{_VERDICT_COLORS['FAIL']};color:white;padding:2px 8px;border-radius:4px;font-size:12px;font-weight:bold;">FAIL</span>
            </td>
        </tr>
        """

    predicate_table = f"""
    <table style="border-collapse:collapse;width:100%;font-family:sans-serif;">
        <thead>
            <tr style="background:#f5f5f5;">
                <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Predicate</th>
                <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Verdict</th>
            </tr>
        </thead>
        <tbody>
            {predicate_rows}
        </tbody>
    </table>
    """

    # Fallback warning
    fallback_html = ""
    if result.fallback_used:
        fallback_html = """
        <div style="background:#fff3cd;border:1px solid #ffc107;border-radius:4px;padding:8px 12px;margin-top:8px;font-family:sans-serif;font-size:13px;">
            [WARN] Greedy fallback optimizer was used (z3 not available or sequence too long)
        </div>
        """

    html = f"""
    <div style="border:1px solid #ddd;border-radius:8px;padding:16px;background:#fafafa;margin:8px 0;">
        <div style="font-weight:bold;font-size:16px;margin-bottom:12px;font-family:sans-serif;">
            [SCIENCE] Optimization Result
        </div>
        {metrics_html}
        {predicate_table}
        {fallback_html}
    </div>
    """

    display(HTML(html))


# ==============================================================================
# 3. display_type_check — Styled HTML Type-Check Results
# ==============================================================================

def display_type_check(results: list) -> None:
    """
    Display type-check results as a styled HTML table in Jupyter.

    Shows each predicate with its verdict, violation (if any),
    and knowledge gap (if any), using color-coded badges.

    Args:
        results: List of TypeCheckResult objects

    Raises:
        ImportError: If IPython is not installed
    """
    _check_ipython()

    from IPython.display import HTML, display

    rows_html = ""
    for result in results:
        verdict = result.verdict.value if hasattr(result.verdict, "value") else str(result.verdict)
        verdict_color = _VERDICT_COLORS.get(verdict, "#6c757d")

        violation_text = result.violation or "—"
        gap_text = result.knowledge_gap or "—"

        # Derivation summary
        derivation_text = ""
        if result.derivation:
            steps = [d.get("step", "") for d in result.derivation[:3]]
            derivation_text = " → ".join(steps)
            if len(result.derivation) > 3:
                derivation_text += f" (+{len(result.derivation) - 3} more)"

        rows_html += f"""
        <tr>
            <td style="padding:8px 12px;border-bottom:1px solid #eee;font-family:monospace;font-size:13px;font-weight:500;">
                {result.predicate}
            </td>
            <td style="padding:8px 12px;border-bottom:1px solid #eee;">
                <span style="background-color:{verdict_color};color:white;padding:3px 10px;border-radius:4px;font-size:12px;font-weight:bold;">
                    {verdict}
                </span>
            </td>
            <td style="padding:8px 12px;border-bottom:1px solid #eee;font-size:12px;color:#666;max-width:300px;word-break:break-word;">
                {violation_text}
            </td>
            <td style="padding:8px 12px;border-bottom:1px solid #eee;font-size:12px;color:#888;max-width:250px;word-break:break-word;">
                {gap_text}
            </td>
            <td style="padding:8px 12px;border-bottom:1px solid #eee;font-size:11px;color:#aaa;max-width:200px;word-break:break-word;">
                {derivation_text}
            </td>
        </tr>
        """

    # Overall verdict
    from biocompiler.shared.types import combined_verdict
    verdicts = [r.verdict for r in results]
    overall = combined_verdict(verdicts)
    overall_color = _VERDICT_COLORS.get(overall.value, "#6c757d")

    html = f"""
    <div style="border:1px solid #ddd;border-radius:8px;padding:16px;background:#fafafa;margin:8px 0;">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:12px;">
            <div style="font-weight:bold;font-size:16px;font-family:sans-serif;">
                [PASS] Type-Check Results
            </div>
            <div>
                <span style="background-color:{overall_color};color:white;padding:4px 14px;border-radius:6px;font-size:14px;font-weight:bold;">
                    Overall: {overall.value}
                </span>
            </div>
        </div>
        <table style="border-collapse:collapse;width:100%;font-family:sans-serif;">
            <thead>
                <tr style="background:#f5f5f5;">
                    <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Predicate</th>
                    <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Verdict</th>
                    <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Violation</th>
                    <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Knowledge Gap</th>
                    <th style="padding:8px 12px;text-align:left;border-bottom:2px solid #ddd;">Derivation</th>
                </tr>
            </thead>
            <tbody>
                {rows_html}
            </tbody>
        </table>
    </div>
    """

    display(HTML(html))


# ==============================================================================
# 4. plot_gc_content — GC Content Plot
# ==============================================================================

def plot_gc_content(
    sequence: str,
    window_size: int = 100,
) -> "matplotlib.figure.Figure":
    """
    Plot GC content along a DNA sequence using a sliding window.

    Returns a matplotlib Figure object that can be displayed in Jupyter
    or saved to a file.

    Args:
        sequence: DNA sequence string
        window_size: Size of the sliding window in base pairs (default 100)

    Returns:
        matplotlib.figure.Figure with the GC content plot

    Raises:
        ImportError: If matplotlib is not installed
    """
    _check_matplotlib()

    import matplotlib.pyplot as plt
    from biocompiler.sequence.scanner import validate_dna_sequence

    seq = validate_dna_sequence(sequence)
    seq_len = len(seq)

    if window_size > seq_len:
        window_size = max(seq_len, 1)

    # Compute GC content in sliding windows
    positions = []
    gc_values = []

    for i in range(seq_len - window_size + 1):
        window = seq[i:i + window_size]
        gc = (window.count("G") + window.count("C")) / window_size
        positions.append(i + window_size // 2)
        gc_values.append(gc)

    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(positions, gc_values, color="#2196F3", linewidth=1.5, alpha=0.9)
    ax.fill_between(positions, gc_values, alpha=0.15, color="#2196F3")

    # Reference lines
    ax.axhline(y=0.50, color="#666666", linestyle="--", linewidth=0.8, label="50% GC")
    ax.axhline(y=0.30, color="#FF9800", linestyle=":", linewidth=0.8, label="30% GC (low)")
    ax.axhline(y=0.70, color="#FF9800", linestyle=":", linewidth=0.8, label="70% GC (high)")

    ax.set_xlabel("Position (bp)", fontsize=11)
    ax.set_ylabel("GC Content", fontsize=11)
    ax.set_title(f"GC Content Sliding Window (window={window_size} bp, seq_len={seq_len})", fontsize=13)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, seq_len)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ==============================================================================
# 5. plot_codon_usage — Codon Usage Bar Chart
# ==============================================================================

def plot_codon_usage(
    sequence: str,
    organism: str = "Homo_sapiens",
) -> "matplotlib.figure.Figure":
    """
    Plot codon usage as a grouped bar chart, one group per amino acid.

    Shows the frequency of each codon within its amino acid family,
    colored by relative adaptiveness.

    Args:
        sequence: DNA coding sequence
        organism: Organism name for reference adaptiveness data

    Returns:
        matplotlib.figure.Figure with the codon usage plot

    Raises:
        ImportError: If matplotlib is not installed
    """
    _check_matplotlib()

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from biocompiler.sequence.scanner import validate_dna_sequence
    from biocompiler.shared.constants import AA_TO_CODONS
    from .organisms import CODON_ADAPTIVENESS_TABLES

    seq = validate_dna_sequence(sequence)

    # Count codons in the sequence
    codon_counts: dict[str, int] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1

    # Get adaptiveness data for coloring
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})

    # Group codons by amino acid
    aa_order = ["M", "W", "F", "L", "I", "V", "S", "P", "T", "A",
                "Y", "H", "Q", "N", "K", "D", "E", "C", "R", "G"]

    # Only include amino acids that appear in the sequence
    aa_present = []
    for aa in aa_order:
        codons = AA_TO_CODONS.get(aa, [])
        if any(codon_counts.get(c, 0) > 0 for c in codons):
            aa_present.append(aa)

    if not aa_present:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No codons found in sequence", ha="center", va="center", fontsize=14)
        ax.set_title("Codon Usage")
        plt.tight_layout()
        return fig

    # Build plot data
    n_aa = len(aa_present)
    fig, ax = plt.subplots(figsize=(max(12, n_aa * 1.5), 6))

    group_width = 0.8
    cmap = cm.get_cmap("RdYlGn")

    x_tick_positions = []
    x_tick_labels = []

    for idx, aa in enumerate(aa_present):
        codons = AA_TO_CODONS[aa]
        n_codons = len(codons)
        bar_width = group_width / n_codons

        for j, codon in enumerate(codons):
            count = codon_counts.get(codon, 0)
            # Normalize within AA family
            total_aa = sum(codon_counts.get(c, 0) for c in codons)
            freq = count / total_aa if total_aa > 0 else 0.0

            # Color by adaptiveness
            adapt = adaptiveness.get(codon, 0.5)
            color = cmap(adapt)

            x = idx + j * bar_width
            ax.bar(x, freq, width=bar_width * 0.9, color=color, edgecolor="white", linewidth=0.5)

            # Label the codon below the bar if frequency > 0
            if freq > 0:
                ax.text(x + bar_width * 0.45, freq + 0.02, codon,
                        ha="center", va="bottom", fontsize=7, rotation=90,
                        fontfamily="monospace")

        x_tick_positions.append(idx + group_width / 2)
        x_tick_labels.append(aa)

    ax.set_xticks(x_tick_positions)
    ax.set_xticklabels(x_tick_labels, fontsize=11, fontweight="bold")
    ax.set_ylabel("Frequency within AA family", fontsize=11)
    ax.set_title(f"Codon Usage — {organism} (green=high adaptiveness, red=low)", fontsize=13)
    ax.set_ylim(0, 1.2)
    ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    return fig


# ==============================================================================
# 6. interactive_optimize — Interactive Widget for Parameter Tuning
# ==============================================================================

def interactive_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
) -> "ipywidgets.Widget":
    """
    Create an interactive Jupyter widget for optimization parameter tuning.

    Provides sliders for:
    - GC low bound (0.1–0.5)
    - GC high bound (0.5–0.9)
    - CAI threshold (0.1–0.9)

    When parameters change, the optimization is re-run and results
    displayed in real-time.

    Args:
        protein: Amino acid sequence to optimize
        organism: Target organism for codon optimization

    Returns:
        ipywidgets interactive widget

    Raises:
        ImportError: If ipywidgets or IPython is not installed
    """
    _check_ipywidgets()
    _check_ipython()

    import ipywidgets as widgets
    from IPython.display import display, clear_output

    from biocompiler.optimizer import optimize_sequence
    from biocompiler.sequence.scanner import gc_content as _gc_content

    # Create sliders
    gc_lo_slider = widgets.FloatSlider(
        value=0.30, min=0.10, max=0.50, step=0.01,
        description="GC Low:", style={"description_width": "80px"},
        layout=widgets.Layout(width="400px"),
    )
    gc_hi_slider = widgets.FloatSlider(
        value=0.70, min=0.50, max=0.90, step=0.01,
        description="GC High:", style={"description_width": "80px"},
        layout=widgets.Layout(width="400px"),
    )
    cai_slider = widgets.FloatSlider(
        value=0.20, min=0.10, max=0.90, step=0.01,
        description="CAI Thresh:", style={"description_width": "80px"},
        layout=widgets.Layout(width="400px"),
    )

    optimize_button = widgets.Button(
        description="[DNA] Optimize",
        button_style="success",
        layout=widgets.Layout(width="150px", margin="10px 0"),
    )

    output = widgets.Output()

    def on_optimize_click(b):
        with output:
            clear_output(wait=True)
            gc_lo = gc_lo_slider.value
            gc_hi = gc_hi_slider.value
            cai_threshold = cai_slider.value

            print(f"Optimizing {len(protein)} aa protein for {organism}...")
            print(f"  GC range: [{gc_lo:.2f}, {gc_hi:.2f}]")
            print(f"  CAI threshold: {cai_threshold:.2f}")

            try:
                result = optimize_sequence(
                    target_protein=protein,
                    organism=organism,
                    gc_lo=gc_lo,
                    gc_hi=gc_hi,
                    cai_threshold=cai_threshold,
                )

                # Display results
                print(f"\n[PASS] Optimization complete!")
                print(f"  Sequence length: {len(result.sequence)} bp")
                print(f"  CAI: {result.cai:.4f}")
                print(f"  GC: {result.gc_content:.1%}")
                print(f"  Satisfied: {len(result.satisfied_predicates)} predicates")
                print(f"  Failed: {len(result.failed_predicates)} predicates")

                if result.satisfied_predicates:
                    print(f"\n  [OK] {', '.join(result.satisfied_predicates)}")
                if result.failed_predicates:
                    print(f"\n  [FAIL] {', '.join(result.failed_predicates)}")

                if result.fallback_used:
                    print("\n  [WARN] Greedy fallback was used")

                # Try to show as rich HTML
                try:
                    display_optimization_result(result)
                except ImportError:
                    logger.debug("IPython not available for rich HTML display")

            except Exception as e:
                print(f"\n[FAIL] Optimization failed: {e}")

    optimize_button.on_click(on_optimize_click)

    # Layout
    controls = widgets.VBox([
        widgets.HTML("<b>[TOOL] Optimization Parameters</b>"),
        gc_lo_slider,
        gc_hi_slider,
        cai_slider,
        optimize_button,
    ])

    widget = widgets.VBox([
        controls,
        output,
    ])

    return widget
