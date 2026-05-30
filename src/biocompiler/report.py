"""
BioCompiler Report Engine — Interactive HTML Report Generation

Production-grade report generation with:
- Self-contained HTML with embedded CSS and JavaScript
- Interactive verdict dashboard with color-coded predicate results
- Sequence visualization with exon/intron annotations
- Splice isoform diagrams
- GC content sliding-window plot (SVG-based, no external dependencies)
- Codon usage heatmap (SVG-based)
- Certificate provenance display
- Print-friendly layout

Reports are self-contained single HTML files that can be shared,
archived, or opened offline. No external JavaScript or CSS required.
"""

import html
import json
import logging
from datetime import datetime, timezone
from typing import Optional

from .types import Certificate, TypeCheckResult, Verdict, Token, SpliceIsoform
from .scanner import gc_content, scan_sequence
from .translation import translate, compute_cai
from .splicing import compute_splice_isoforms
from .constants import CODON_TABLE, AA_TO_CODONS

logger = logging.getLogger(__name__)


def generate_report(
    sequence: str,
    type_results: Optional[list[TypeCheckResult]] = None,
    certificate: Optional[Certificate] = None,
    organism: str = "Homo_sapiens",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    tokens: Optional[list[Token]] = None,
    isoforms: Optional[list[SpliceIsoform]] = None,
    optimization_result: Optional[dict] = None,
) -> str:
    """
    Generate a self-contained interactive HTML report.

    The report includes:
    - Summary dashboard with overall verdict
    - Predicate results table with PASS/FAIL/UNCERTAIN indicators
    - Sequence statistics (GC, CAI, length)
    - GC content sliding window plot (SVG)
    - Codon usage heatmap (SVG)
    - Exon/intron structure diagram
    - Splice isoform table
    - Certificate details (if available)
    - Token scan results table

    Args:
        sequence: DNA sequence
        type_results: Type-check predicate results
        certificate: Optional certificate
        organism: Target organism
        gene_name: Optional gene name
        exon_boundaries: Optional exon boundary positions
        tokens: Optional scan tokens
        isoforms: Optional splice isoforms
        optimization_result: Optional optimization result dict

    Returns:
        Self-contained HTML string
    """
    seq = sequence.upper()
    gc = gc_content(seq)
    protein = translate(seq)
    cai = compute_cai(seq, organism)

    if exon_boundaries is None:
        exon_boundaries = [(0, len(seq))]

    if tokens is None:
        tokens = scan_sequence(seq)

    if isoforms is None and len(exon_boundaries) > 1:
        isoforms = compute_splice_isoforms(seq, exon_boundaries)

    # Compute GC sliding window data
    gc_window_data = _compute_gc_windows(seq, window_size=50)

    # Compute codon usage data
    codon_usage_data = _compute_codon_usage(seq)

    # Build HTML
    report_html = _build_html(
        seq=seq,
        gc=gc,
        cai=cai,
        protein=protein,
        organism=organism,
        gene_name=gene_name,
        exon_boundaries=exon_boundaries,
        type_results=type_results or [],
        certificate=certificate,
        tokens=tokens,
        isoforms=isoforms or [],
        gc_window_data=gc_window_data,
        codon_usage_data=codon_usage_data,
        optimization_result=optimization_result,
    )

    return report_html


def _compute_gc_windows(seq: str, window_size: int = 50) -> list[tuple[int, float]]:
    """Compute GC content in sliding windows."""
    results = []
    for i in range(0, len(seq) - window_size + 1, max(1, window_size // 4)):
        window = seq[i:i + window_size]
        gc_val = gc_content(window)
        results.append((i, gc_val))
    return results


def _compute_codon_usage(seq: str) -> dict[str, dict[str, float]]:
    """Compute codon frequency per amino acid."""
    aa_codons: dict[str, dict[str, int]] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa and aa != "*":
            if aa not in aa_codons:
                aa_codons[aa] = {}
            aa_codons[aa][codon] = aa_codons[aa].get(codon, 0) + 1

    # Convert to frequencies
    aa_freqs: dict[str, dict[str, float]] = {}
    for aa, codons in aa_codons.items():
        total = sum(codons.values())
        aa_freqs[aa] = {c: count / total for c, count in codons.items()}

    return aa_freqs


def _build_html(
    seq: str, gc: float, cai: float, protein: str,
    organism: str, gene_name: Optional[str],
    exon_boundaries: list[tuple[int, int]],
    type_results: list[TypeCheckResult],
    certificate: Optional[Certificate],
    tokens: list[Token],
    isoforms: list[SpliceIsoform],
    gc_window_data: list[tuple[int, float]],
    codon_usage_data: dict[str, dict[str, float]],
    optimization_result: Optional[dict],
) -> str:
    """Build the complete HTML report."""

    # Overall verdict
    overall_verdict = "N/A"
    overall_class = "uncertain"
    if type_results:
        verdicts = [r.verdict for r in type_results]
        overall = Verdict.PASS
        for v in verdicts:
            if overall == Verdict.FAIL or v == Verdict.FAIL:
                overall = Verdict.FAIL
            elif overall == Verdict.UNCERTAIN or v == Verdict.UNCERTAIN:
                overall = Verdict.UNCERTAIN
        overall_verdict = overall.value
        overall_class = overall.value.lower()

    # GC plot SVG
    gc_plot_svg = _generate_gc_plot_svg(gc_window_data, gc)

    # Codon heatmap SVG
    codon_heatmap_svg = _generate_codon_heatmap_svg(codon_usage_data)

    # Sequence structure SVG
    structure_svg = _generate_structure_svg(seq, exon_boundaries, tokens)

    # Build predicate rows
    predicate_rows = ""
    for r in type_results:
        css_class = r.verdict.value.lower()
        symbol = {"PASS": "&#10003;", "FAIL": "&#10007;", "UNCERTAIN": "?"}[r.verdict.value]
        violation_html = f'<div class="violation">{html.escape(r.violation or "")}</div>' if r.violation else ""
        gap_html = f'<div class="knowledge-gap">Gap: {html.escape(r.knowledge_gap or "")}</div>' if r.knowledge_gap else ""
        predicate_rows += f"""
            <tr class="{css_class}">
                <td class="verdict-cell"><span class="badge {css_class}">{symbol}</span></td>
                <td>{html.escape(r.predicate)}</td>
                <td>{violation_html}{gap_html}</td>
            </tr>"""

    # Build token rows
    token_rows = ""
    for t in tokens[:50]:  # Limit display
        token_rows += f"""
            <tr>
                <td>{t.position}</td>
                <td>{html.escape(t.element_type)}</td>
                <td class="mono">{html.escape(t.match_sequence)}</td>
                <td>{t.score:.2f}</td>
                <td>{t.frame if t.frame is not None else "-"}</td>
                <td>{html.escape(t.strand)}</td>
            </tr>"""

    # Build isoform rows
    isoform_rows = ""
    for i, iso in enumerate(isoforms[:20]):
        isoform_rows += f"""
            <tr>
                <td>{i + 1}</td>
                <td>{len(iso.sequence)} nt</td>
                <td>{len(iso.exon_boundaries)}</td>
                <td>{html.escape(str(iso.parse_path))}</td>
                <td>{iso.score:.2f}</td>
            </tr>"""

    # Certificate section
    cert_section = ""
    if certificate:
        cert_section = f"""
        <section class="section">
            <h2>Certificate</h2>
            <div class="cert-grid">
                <div class="cert-field">
                    <span class="cert-label">Design ID</span>
                    <span class="cert-value mono">{html.escape(certificate.design_id[:32])}...</span>
                </div>
                <div class="cert-field">
                    <span class="cert-label">Version</span>
                    <span class="cert-value">{html.escape(certificate.version)}</span>
                </div>
                <div class="cert-field">
                    <span class="cert-label">Timestamp</span>
                    <span class="cert-value">{html.escape(certificate.provenance.get('timestamp', 'N/A'))}</span>
                </div>
                <div class="cert-field">
                    <span class="cert-label">Tool</span>
                    <span class="cert-value">{html.escape(certificate.provenance.get('tool', 'N/A'))}</span>
                </div>
            </div>
        </section>"""

    # Optimization section
    opt_section = ""
    if optimization_result:
        opt_section = f"""
        <section class="section">
            <h2>Optimization Results</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{optimization_result.get('cai', 0):.4f}</div>
                    <div class="stat-label">CAI</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{optimization_result.get('gc_content', 0):.1%}</div>
                    <div class="stat-label">GC Content</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(optimization_result.get('satisfied', []))}</div>
                    <div class="stat-label">Satisfied</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(optimization_result.get('failed', []))}</div>
                    <div class="stat-label">Failed</div>
                </div>
            </div>
        </section>"""

    # Build full HTML
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>BioCompiler Report{f' — {html.escape(gene_name)}' if gene_name else ''}</title>
<style>
:root {{
    --pass: #16a34a;
    --fail: #dc2626;
    --uncertain: #d97706;
    --bg: #ffffff;
    --surface: #f8fafc;
    --border: #e2e8f0;
    --text: #1e293b;
    --text-secondary: #64748b;
    --accent: #2563eb;
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background: var(--surface);
    color: var(--text);
    line-height: 1.6;
    padding: 2rem;
    max-width: 1200px;
    margin: 0 auto;
}}
.header {{
    background: linear-gradient(135deg, #1e293b, #334155);
    color: white;
    padding: 2rem;
    border-radius: 12px;
    margin-bottom: 2rem;
    display: flex;
    justify-content: space-between;
    align-items: center;
}}
.header h1 {{ font-size: 1.75rem; font-weight: 700; }}
.header .subtitle {{ color: #94a3b8; margin-top: 0.25rem; }}
.verdict-badge {{
    font-size: 1.5rem;
    font-weight: 800;
    padding: 0.75rem 2rem;
    border-radius: 8px;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}}
.verdict-badge.pass {{ background: var(--pass); color: white; }}
.verdict-badge.fail {{ background: var(--fail); color: white; }}
.verdict-badge.uncertain {{ background: var(--uncertain); color: white; }}
.section {{
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.5rem;
    margin-bottom: 1.5rem;
}}
.section h2 {{
    font-size: 1.25rem;
    font-weight: 600;
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--border);
}}
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 1rem;
}}
.stat-card {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1rem;
    text-align: center;
}}
.stat-value {{
    font-size: 1.75rem;
    font-weight: 700;
    color: var(--accent);
}}
.stat-label {{
    font-size: 0.8rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin-top: 0.5rem;
}}
th, td {{
    padding: 0.5rem 0.75rem;
    text-align: left;
    border-bottom: 1px solid var(--border);
}}
th {{
    background: var(--surface);
    font-weight: 600;
    font-size: 0.85rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    color: var(--text-secondary);
}}
tr.pass {{ background: #f0fdf4; }}
tr.fail {{ background: #fef2f2; }}
tr.uncertain {{ background: #fffbeb; }}
.badge {{
    display: inline-block;
    padding: 0.15rem 0.5rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.9rem;
}}
.badge.pass {{ background: #dcfce7; color: var(--pass); }}
.badge.fail {{ background: #fee2e2; color: var(--fail); }}
.badge.uncertain {{ background: #fef3c7; color: var(--uncertain); }}
.violation {{ color: var(--fail); font-size: 0.85rem; margin-top: 0.25rem; }}
.knowledge-gap {{ color: var(--uncertain); font-size: 0.85rem; }}
.mono {{ font-family: 'SF Mono', 'Fira Code', monospace; font-size: 0.9rem; }}
.cert-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 1rem;
}}
.cert-field {{
    padding: 0.75rem;
    background: var(--surface);
    border-radius: 6px;
}}
.cert-label {{
    display: block;
    font-size: 0.75rem;
    color: var(--text-secondary);
    text-transform: uppercase;
}}
.cert-value {{ font-weight: 500; }}
.svg-container {{
    overflow-x: auto;
    margin: 1rem 0;
}}
.svg-container svg {{ max-width: 100%; height: auto; }}
.sequence-display {{
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.8rem;
    line-height: 1.4;
    background: var(--surface);
    padding: 1rem;
    border-radius: 6px;
    overflow-x: auto;
    word-break: break-all;
    max-height: 200px;
    overflow-y: auto;
}}
.footer {{
    text-align: center;
    color: var(--text-secondary);
    font-size: 0.85rem;
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
}}
@media print {{
    body {{ padding: 0; }}
    .section {{ break-inside: avoid; }}
}}
</style>
</head>
<body>

<div class="header">
    <div>
        <h1>BioCompiler Report</h1>
        <div class="subtitle">
            {html.escape(gene_name or 'Designed Sequence')} &mdash;
            {html.escape(organism)} &mdash;
            {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
        </div>
    </div>
    <div class="verdict-badge {overall_class}">{overall_verdict}</div>
</div>

<!-- Summary Statistics -->
<section class="section">
    <h2>Sequence Summary</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{len(seq):,}</div>
            <div class="stat-label">Length (bp)</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{gc:.1%}</div>
            <div class="stat-label">GC Content</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{cai:.4f}</div>
            <div class="stat-label">CAI</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(protein)}</div>
            <div class="stat-label">Protein (aa)</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(exon_boundaries)}</div>
            <div class="stat-label">Exons</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(tokens)}</div>
            <div class="stat-label">Tokens</div>
        </div>
    </div>
</section>

{opt_section}

<!-- Predicate Results -->
<section class="section">
    <h2>Type-Check Results</h2>
    <table>
        <thead>
            <tr>
                <th>Verdict</th>
                <th>Predicate</th>
                <th>Details</th>
            </tr>
        </thead>
        <tbody>
            {predicate_rows if predicate_rows else '<tr><td colspan="3">No type-check results</td></tr>'}
        </tbody>
    </table>
</section>

<!-- GC Content Plot -->
<section class="section">
    <h2>GC Content Profile</h2>
    <div class="svg-container">
        {gc_plot_svg}
    </div>
</section>

<!-- Sequence Structure -->
<section class="section">
    <h2>Sequence Structure</h2>
    <div class="svg-container">
        {structure_svg}
    </div>
</section>

<!-- Codon Usage -->
<section class="section">
    <h2>Codon Usage</h2>
    <div class="svg-container">
        {codon_heatmap_svg}
    </div>
</section>

<!-- Splice Isoforms -->
{f'''<section class="section">
    <h2>Splice Isoforms ({len(isoforms)})</h2>
    <table>
        <thead>
            <tr><th>#</th><th>Length</th><th>Exons</th><th>Path</th><th>Score</th></tr>
        </thead>
        <tbody>
            {isoform_rows if isoform_rows else '<tr><td colspan="5">No alternative isoforms</td></tr>'}
        </tbody>
    </table>
</section>''' if isoforms else ''}

<!-- Scan Tokens -->
<section class="section">
    <h2>Scan Tokens ({len(tokens)})</h2>
    <table>
        <thead>
            <tr><th>Position</th><th>Type</th><th>Match</th><th>Score</th><th>Frame</th><th>Strand</th></tr>
        </thead>
        <tbody>
            {token_rows if token_rows else '<tr><td colspan="6">No tokens found</td></tr>'}
        </tbody>
    </table>
</section>

{cert_section}

<!-- Sequence -->
<section class="section">
    <h2>Sequence</h2>
    <div class="sequence-display">{_format_sequence_html(seq, exon_boundaries)}</div>
</section>

<!-- Protein -->
<section class="section">
    <h2>Protein Translation</h2>
    <div class="sequence-display">{html.escape(protein)}</div>
</section>

<div class="footer">
    Generated by BioCompiler v2.2.0 &mdash; Machine-Verified Gene Design
</div>

</body>
</html>"""


def _format_sequence_html(seq: str, exon_boundaries: list[tuple[int, int]]) -> str:
    """Format sequence with exon highlighting."""
    exon_positions = set()
    for start, end in exon_boundaries:
        for i in range(start, min(end, len(seq))):
            exon_positions.add(i)

    result = []
    for i, base in enumerate(seq):
        if i in exon_positions:
            result.append(f'<span style="color:var(--accent);font-weight:600">{html.escape(base)}</span>')
        else:
            result.append(f'<span style="color:var(--text-secondary)">{html.escape(base)}</span>')

    # Add line breaks every 60 characters
    formatted = ""
    for i, char in enumerate(result):
        if i > 0 and i % 60 == 0:
            formatted += "\n"
        formatted += char
    return formatted


def _generate_gc_plot_svg(data: list[tuple[int, float]], avg_gc: float) -> str:
    """Generate SVG line chart for GC content."""
    if not data:
        return "<p>No GC data available</p>"

    width = 800
    height = 200
    margin = 40
    plot_w = width - 2 * margin
    plot_h = height - 2 * margin

    max_pos = max(d[0] for d in data)
    min_pos = min(d[0] for d in data)
    pos_range = max(max_pos - min_pos, 1)

    # Build path
    points = []
    for pos, gc_val in data:
        x = margin + (pos - min_pos) / pos_range * plot_w
        y = margin + plot_h - gc_val * plot_h
        points.append(f"{x:.1f},{y:.1f}")

    path_d = "M " + " L ".join(points)

    # Average line
    avg_y = margin + plot_h - avg_gc * plot_h

    # Threshold lines (0.3 and 0.7)
    y30 = margin + plot_h - 0.3 * plot_h
    y70 = margin + plot_h - 0.7 * plot_h

    return f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    <!-- Grid lines -->
    <line x1="{margin}" y1="{y30}" x2="{width - margin}" y2="{y30}" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="4,4"/>
    <line x1="{margin}" y1="{y70}" x2="{width - margin}" y2="{y70}" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="4,4"/>
    <!-- Average line -->
    <line x1="{margin}" y1="{avg_y:.1f}" x2="{width - margin}" y2="{avg_y:.1f}" stroke="#2563eb" stroke-width="1.5" stroke-dasharray="8,4"/>
    <!-- GC curve -->
    <path d="{path_d}" fill="none" stroke="#16a34a" stroke-width="2"/>
    <!-- Axes -->
    <line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="#64748b" stroke-width="1"/>
    <line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="#64748b" stroke-width="1"/>
    <!-- Labels -->
    <text x="5" y="{margin + 4}" font-size="10" fill="#64748b">1.0</text>
    <text x="5" y="{height - margin}" font-size="10" fill="#64748b">0.0</text>
    <text x="{margin}" y="{height - 5}" font-size="10" fill="#64748b">0</text>
    <text x="{width - margin - 30}" y="{height - 5}" font-size="10" fill="#64748b">{max_pos}</text>
    <text x="{width // 2}" y="15" font-size="11" fill="#1e293b" text-anchor="middle" font-weight="600">GC Content (sliding window)</text>
    <text x="{margin - 5}" y="{y30 + 3}" font-size="9" fill="#94a3b8" text-anchor="end">0.30</text>
    <text x="{margin - 5}" y="{y70 + 3}" font-size="9" fill="#94a3b8" text-anchor="end">0.70</text>
</svg>"""


def _generate_codon_heatmap_svg(codon_usage: dict[str, dict[str, float]]) -> str:
    """Generate SVG codon usage heatmap."""
    if not codon_usage:
        return "<p>No codon usage data</p>"

    aas = sorted(codon_usage.keys())
    cell_w = 40
    cell_h = 25
    margin_left = 30
    margin_top = 40
    width = margin_left + 64 * cell_w + 20
    height = margin_top + len(aas) * cell_h + 20

    cells = ""
    for row, aa in enumerate(aas):
        codons = codon_usage[aa]
        y = margin_top + row * cell_h
        # AA label
        cells += f'<text x="{margin_left - 5}" y="{y + cell_h // 2 + 4}" font-size="11" fill="#1e293b" text-anchor="end" font-weight="600">{aa}</text>'
        for col, (codon, freq) in enumerate(sorted(codons.items())):
            x = margin_left + col * cell_w
            # Color intensity based on frequency
            intensity = min(freq * 255, 255)
            color = f"rgb({255 - int(intensity * 0.8)}, {255 - int(intensity * 0.3)}, {255 - int(intensity * 0.8)})"
            cells += f'<rect x="{x}" y="{y}" width="{cell_w - 1}" height="{cell_h - 1}" fill="{color}" rx="2"/>'
            cells += f'<text x="{x + cell_w // 2}" y="{y + cell_h // 2 + 4}" font-size="8" fill="#1e293b" text-anchor="middle">{codon}</text>'

    return f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    <text x="{width // 2}" y="15" font-size="11" fill="#1e293b" text-anchor="middle" font-weight="600">Codon Usage Heatmap</text>
    {cells}
</svg>"""


def _generate_structure_svg(
    seq: str,
    exon_boundaries: list[tuple[int, int]],
    tokens: list[Token],
) -> str:
    """Generate SVG diagram of exon/intron structure."""
    width = 800
    height = 80
    margin = 40
    bar_y = 30
    bar_h = 20

    seq_len = max(len(seq), 1)
    scale = (width - 2 * margin) / seq_len

    elements = ""

    # Draw full sequence bar (introns)
    elements += f'<rect x="{margin}" y="{bar_y}" width="{width - 2 * margin}" height="{bar_h}" fill="#e2e8f0" rx="3"/>'

    # Draw exons
    colors = ["#2563eb", "#7c3aed", "#0891b2", "#059669", "#d97706"]
    for i, (start, end) in enumerate(exon_boundaries):
        x = margin + start * scale
        w = (end - start) * scale
        color = colors[i % len(colors)]
        elements += f'<rect x="{x:.1f}" y="{bar_y}" width="{w:.1f}" height="{bar_h}" fill="{color}" rx="2"/>'
        # Label
        label_x = x + w / 2
        elements += f'<text x="{label_x:.1f}" y="{bar_y - 5}" font-size="9" fill="#1e293b" text-anchor="middle">Exon {i + 1}</text>'

    # Draw splice donor/acceptor markers
    for t in tokens:
        if t.element_type in ("splice_donor", "splice_acceptor"):
            x = margin + t.position * scale
            marker_color = "#dc2626" if t.element_type == "splice_donor" else "#16a34a"
            elements += f'<line x1="{x:.1f}" y1="{bar_y + bar_h}" x2="{x:.1f}" y2="{bar_y + bar_h + 8}" stroke="{marker_color}" stroke-width="1.5"/>'

    # Position markers
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        pos = int(frac * seq_len)
        x = margin + frac * (width - 2 * margin)
        elements += f'<text x="{x:.1f}" y="{height - 5}" font-size="9" fill="#64748b" text-anchor="middle">{pos}</text>'

    return f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    {elements}
</svg>"""
