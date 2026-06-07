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
import logging
from datetime import datetime, timezone
from typing import Optional

from ..types import Certificate, TypeCheckResult, Token, SpliceIsoform, combined_verdict
from ..scanner import gc_content, scan_sequence
from ..translation import translate, compute_cai
from ..splicing import compute_splice_isoforms
from ..constants import CODON_TABLE
from ..engine_base import BaseEngineResult
from ..organism_config import get_organism_config, is_eukaryotic_organism
from .. import __version__

logger = logging.getLogger(__name__)

__all__ = ["generate_report", "_build_organism_aware_constraints", "_build_cai_comparison_section"]

# ── Display & layout constants ──────────────────────────────────────────
_MAX_TOKEN_DISPLAY: int = 50
_MAX_ISOFORM_DISPLAY: int = 20
_SEQUENCE_LINE_WIDTH: int = 60
_STOP_CODON_MARKER: str = "*"
_DESIGN_ID_DISPLAY_LENGTH: int = 32

# ── GC window computation constants ─────────────────────────────────────
_DEFAULT_GC_WINDOW_SIZE: int = 50
_GC_WINDOW_STEP_DIVISOR: int = 4

# ── GC plot SVG constants ───────────────────────────────────────────────
_GC_PLOT_WIDTH: int = 800
_GC_PLOT_HEIGHT: int = 200
_GC_PLOT_MARGIN: int = 40
_GC_LOW_THRESHOLD: float = 0.3
_GC_HIGH_THRESHOLD: float = 0.7

# ── Verdict badge symbols ──────────────────────────────────────────────
_VERDICT_SYMBOLS: dict[str, str] = {
    "PASS": "&#10003;",
    "LIKELY_PASS": "&#10003;~",
    "UNCERTAIN": "?",
    "LIKELY_FAIL": "&#10007;~",
    "FAIL": "&#10007;",
}

# ── Structure SVG constants ────────────────────────────────────────────
_STRUCT_SVG_WIDTH: int = 800
_STRUCT_SVG_HEIGHT: int = 80
_STRUCT_SVG_MARGIN: int = 40
_STRUCT_SVG_BAR_Y: int = 30
_STRUCT_SVG_BAR_H: int = 20
_STRUCT_EXON_COLORS: list[str] = ["#2563eb", "#7c3aed", "#0891b2", "#059669", "#d97706"]
_POSITION_TICK_FRACTIONS: list[float] = [0, 0.25, 0.5, 0.75, 1.0]

# ── Codon heatmap SVG constants ────────────────────────────────────────
_HEATMAP_CELL_W: int = 40
_HEATMAP_CELL_H: int = 25
_HEATMAP_MARGIN_LEFT: int = 30
_HEATMAP_MARGIN_TOP: int = 40
_HEATMAP_MAX_CODON_COLS: int = 64

# ── Organism-aware constraint mapping ──────────────────────────────────
# Constraints that are only relevant for eukaryotic organisms.
_EUKARYOTE_ONLY_CONSTRAINTS: dict[str, str] = {
    "NoCrypticSpliceConstraint": "eukaryotic organism (splice sites are eukaryote-specific)",
    "NoCpGIslandConstraint": "eukaryotic organism (CpG islands are primarily a mammalian concern)",
    "NoTRunConstraint": "eukaryotic organism (poly-T termination signals are eukaryote-specific)",
}

_ALL_CONSTRAINTS: dict[str, str] = {
    "TranslationConstraint": "Ensures every codon translates to the correct amino acid",
    "NoRestrictionSiteConstraint": "Eliminates restriction enzyme recognition sites",
    "GCRangeConstraint": "Keeps GC content within organism-specific range",
    "NoCrypticSpliceConstraint": "Eliminates cryptic splice donor/acceptor sites",
    "NoCpGIslandConstraint": "Prevents CpG island formation (methylation risk)",
    "NoATTTAMotifConstraint": "Removes ATTTA mRNA instability motifs",
    "NoTRunConstraint": "Prevents poly-T transcription termination signals",
    "MaximizeCAI": "Maximize Codon Adaptation Index (soft)",
    "MinimizeCpG": "Minimize CpG dinucleotide count (soft)",
    "MinimizeMRNADG": "Minimize mRNA 5' folding free energy (soft)",
}


def _build_organism_aware_constraints(
    organism: str,
    optimization_result: Optional[dict[str, object]] = None,
) -> dict[str, object]:
    """Build organism-aware constraint status information."""
    config = get_organism_config(organism)
    is_euk = config.is_eukaryote

    applied_names: list[str] = []
    if optimization_result and isinstance(optimization_result.get("applied_constraints"), list):
        applied_names = [c if isinstance(c, str) else str(c) for c in optimization_result["applied_constraints"]]

    skipped_names: list[str] = []
    if optimization_result and isinstance(optimization_result.get("skipped_constraints"), list):
        skipped_names = [c if isinstance(c, str) else str(c) for c in optimization_result["skipped_constraints"]]

    if not applied_names and not skipped_names:
        for cname in _ALL_CONSTRAINTS:
            if cname in _EUKARYOTE_ONLY_CONSTRAINTS and not is_euk:
                skipped_names.append(cname)
            else:
                applied_names.append(cname)

    active_constraints = [(name, _ALL_CONSTRAINTS.get(name, "Optimization constraint")) for name in applied_names]
    skipped_constraints = []
    for name in skipped_names:
        if name in _EUKARYOTE_ONLY_CONSTRAINTS:
            reason = f"{name} \u2014 skipped for {_EUKARYOTE_ONLY_CONSTRAINTS[name]}"
        else:
            reason = f"{name} \u2014 skipped for {config.domain} organism"
        skipped_constraints.append((name, reason))

    return {
        "organism": organism,
        "domain": config.domain,
        "is_eukaryote": is_euk,
        "config_name": config.name,
        "active_constraints": active_constraints,
        "skipped_constraints": skipped_constraints,
    }


def _build_cai_comparison_section(
    organism: str,
    optimization_result: Optional[dict[str, object]] = None,
) -> Optional[dict[str, float]]:
    """Build CAI comparison data between all-constraints and organism-aware runs."""
    if not optimization_result:
        return None
    cai_all = optimization_result.get("cai_all_constraints")
    cai_aware = optimization_result.get("cai_organism_aware_constraints")
    if cai_aware is None and cai_all is not None:
        cai_aware = optimization_result.get("cai", cai_all)
    if cai_all is not None and cai_aware is not None:
        try:
            cai_all_f = float(cai_all)
            cai_aware_f = float(cai_aware)
            recovery = cai_aware_f - cai_all_f
            return {"cai_all_constraints": cai_all_f, "cai_organism_aware": cai_aware_f, "recovery": recovery}
        except (ValueError, TypeError):
            return None
    return None


def generate_report(
    sequence: str,
    type_results: Optional[list[TypeCheckResult]] = None,
    certificate: Optional[Certificate] = None,
    organism: str = "Homo_sapiens",
    gene_name: Optional[str] = None,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    tokens: Optional[list[Token]] = None,
    isoforms: Optional[list[SpliceIsoform]] = None,
    optimization_result: Optional[dict[str, object]] = None,
    engine_results: Optional[list[BaseEngineResult]] = None,
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
    - Engine analysis results (if provided via engine_results)

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
        engine_results: Optional list of BaseEngineResult objects from
            analysis engines (ESMFold, FoldX, CamSol, Immunogenicity, etc.)
            Uses unified field names (primary_score, classification) for
            display.

    Returns:
        Self-contained HTML string
    """
    seq = sequence.upper()

    try:
        gc = gc_content(seq)
    except Exception:
        logger.warning("Failed to compute GC content", exc_info=True)
        gc = 0.0

    try:
        protein = translate(seq)
    except Exception:
        logger.warning("Failed to translate sequence", exc_info=True)
        protein = ""

    try:
        cai = compute_cai(seq, organism)
    except Exception:
        logger.warning("Failed to compute CAI for organism=%s", organism, exc_info=True)
        cai = 0.0

    logger.info("Generating report: %s bp, GC=%.1f%%, CAI=%.4f, organism=%s",
                len(seq), gc * 100, cai, organism)

    if exon_boundaries is None:
        exon_boundaries = [(0, len(seq))]

    if tokens is None:
        try:
            tokens = scan_sequence(seq)
        except Exception:
            logger.warning("Failed to scan sequence for tokens", exc_info=True)
            tokens = []

    if isoforms is None and len(exon_boundaries) > 1:
        try:
            isoforms = compute_splice_isoforms(seq, exon_boundaries)
        except Exception:
            logger.warning("Failed to compute splice isoforms", exc_info=True)
            isoforms = []

    # Compute GC sliding window data
    gc_window_data = _compute_gc_windows(seq, window_size=_DEFAULT_GC_WINDOW_SIZE)

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
        engine_results=engine_results or [],
    )

    return report_html


def _compute_gc_windows(seq: str, window_size: int = 50) -> list[tuple[int, float]]:
    """Compute GC content in sliding windows.

    Args:
        seq: DNA sequence to analyze.
        window_size: Width of the sliding window in base pairs.

    Returns:
        List of (position, gc_fraction) tuples, one per window step.
    """
    results = []
    for i in range(0, len(seq) - window_size + 1, max(1, window_size // _GC_WINDOW_STEP_DIVISOR)):
        window = seq[i:i + window_size]
        gc_val = gc_content(window)
        results.append((i, gc_val))
    return results


def _compute_codon_usage(seq: str) -> dict[str, dict[str, float]]:
    """Compute codon frequency per amino acid.

    Args:
        seq: DNA sequence to analyze.

    Returns:
        Nested dict mapping amino-acid → {codon → frequency}.
    """
    aa_codons: dict[str, dict[str, int]] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa and aa != _STOP_CODON_MARKER:
            if aa not in aa_codons:
                aa_codons[aa] = {}
            aa_codons[aa][codon] = aa_codons[aa].get(codon, 0) + 1

    # Convert to frequencies
    aa_freqs: dict[str, dict[str, float]] = {}
    for aa, codons in aa_codons.items():
        total = sum(codons.values())
        aa_freqs[aa] = {c: count / total for c, count in codons.items()}

    return aa_freqs


def _element_type_description(element_type: str) -> str:
    """Return a human-readable description for a token element type.

    Args:
        element_type: Internal element type identifier string.

    Returns:
        Human-readable description string.
    """
    descriptions = {
        "start_codon": "Start codon — initiates translation (ATG)",
        "stop_codon": "Stop codon — terminates translation (TAA, TAG, TGA)",
        "splice_donor": "Splice donor site — 5' end of intron (GT consensus)",
        "splice_acceptor": "Splice acceptor site — 3' end of intron (AG consensus)",
        "kozak": "Kozak consensus sequence — enhances translation initiation",
        "cpg_island": "CpG island — region of high C-G dinucleotide frequency, often regulatory",
        "u_rich_region": "U-rich region — potential mRNA instability element",
        "tata_box": "TATA box — core promoter element for transcription initiation",
        "poly_a_signal": "Polyadenylation signal — directs 3' end processing of mRNA",
        "restriction_site": "Restriction enzyme recognition site",
        "repeat": "Repetitive sequence element",
        "promoter": "Promoter region — binds RNA polymerase for transcription",
        "enhancer": "Enhancer element — regulatory DNA that increases transcription",
        "silencer": "Silencer element — regulatory DNA that decreases transcription",
        "insulator": "Insulator element — blocks interaction between regulatory regions",
    }
    return descriptions.get(element_type, f"Biological element of type: {element_type}")


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
    optimization_result: Optional[dict[str, object]],
    engine_results: Optional[list[BaseEngineResult]] = None,
) -> str:
    """Build the complete self-contained HTML report.

    Args:
        seq: Uppercase DNA sequence.
        gc: Overall GC fraction.
        cai: Codon Adaptation Index value.
        protein: Translated amino-acid sequence.
        organism: Target organism name.
        gene_name: Optional gene name for the header.
        exon_boundaries: Exon start/end positions.
        type_results: Type-check predicate results.
        certificate: Optional design certificate.
        tokens: Scanned biological elements.
        isoforms: Computed splice isoforms.
        gc_window_data: GC sliding-window data points.
        codon_usage_data: Per-amino-acid codon frequencies.
        optimization_result: Optional optimization metrics dict.
        engine_results: Optional analysis engine results.

    Returns:
        Complete HTML document as a string.
    """

    # Overall verdict
    overall_verdict = "N/A"
    overall_class = "uncertain"
    if type_results:
        verdicts = [r.verdict for r in type_results]
        overall = combined_verdict(verdicts)
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
        symbol = _VERDICT_SYMBOLS.get(r.verdict.value, "?")
        violation_html = f'<div class="violation">{html.escape(r.violation or "")}</div>' if r.violation else ""
        gap_html = f'<div class="knowledge-gap">Gap: {html.escape(r.knowledge_gap or "")}</div>' if r.knowledge_gap else ""
        # Build derivation tooltip text
        derivation_parts = []
        derivation_parts.append(f"Predicate: {r.predicate}")
        derivation_parts.append(f"Verdict: {r.verdict.value}")
        if r.violation:
            derivation_parts.append(f"Violation: {r.violation}")
        if r.knowledge_gap:
            derivation_parts.append(f"Knowledge gap: {r.knowledge_gap}")
        tooltip_text = html.escape(" | ".join(derivation_parts))
        predicate_rows += f"""
            <tr class="{css_class}" data-verdict="{css_class}">
                <td class="verdict-cell"><span class="badge {css_class}" data-tooltip="{tooltip_text}">{symbol}</span></td>
                <td>{html.escape(r.predicate)}</td>
                <td>{violation_html}{gap_html}</td>
            </tr>"""

    # Build token rows
    token_rows = ""
    for t in tokens[:_MAX_TOKEN_DISPLAY]:  # Limit display
        token_rows += f"""
            <tr data-element-type="{html.escape(t.element_type.lower())}" data-match-sequence="{html.escape(t.match_sequence.lower())}">
                <td>{t.position}</td>
                <td class="token-type-cell" data-tooltip="{html.escape(_element_type_description(t.element_type))}">{html.escape(t.element_type)}</td>
                <td class="mono">{html.escape(t.match_sequence)}</td>
                <td>{t.score:.2f}</td>
                <td>{t.frame if t.frame is not None else "-"}</td>
                <td>{html.escape(t.strand)}</td>
            </tr>"""

    # Build isoform rows
    isoform_rows = ""
    for i, iso in enumerate(isoforms[:_MAX_ISOFORM_DISPLAY]):
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
                    <span class="cert-value mono">{html.escape(certificate.design_id[:_DESIGN_ID_DISPLAY_LENGTH])}...</span>
                </div>
                <div class="cert-field">
                    <span class="cert-label">Version</span>
                    <span class="cert-value">{html.escape(certificate.version)}</span>
                </div>
                <div class="cert-field">
                    <span class="cert-label">Timestamp</span>
                    <span class="cert-value">{html.escape(certificate.provenance.get('timestamp', 'N/A') if isinstance(certificate.provenance, dict) else 'N/A')}</span>
                </div>
                <div class="cert-field">
                    <span class="cert-label">Tool</span>
                    <span class="cert-value">{html.escape(certificate.provenance.get('tool', 'N/A') if isinstance(certificate.provenance, dict) else 'N/A')}</span>
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

    # Organism-aware constraints section
    org_constraints = _build_organism_aware_constraints(organism, optimization_result)
    cai_comparison = _build_cai_comparison_section(organism, optimization_result)

    org_constraints_section = f"""
        <section class="section">
            <h2>Organism-Aware Constraints</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{html.escape(org_constraints['config_name'])}</div>
                    <div class="stat-label">Organism</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{html.escape(str(org_constraints['domain']))}</div>
                    <div class="stat-label">Domain</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(org_constraints['active_constraints'])}</div>
                    <div class="stat-label">Active Constraints</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(org_constraints['skipped_constraints'])}</div>
                    <div class="stat-label">Skipped Constraints</div>
                </div>
            </div>
            <h3>Constraints Applied</h3>
            <table>
                <thead>
                    <tr><th>Status</th><th>Constraint</th><th>Description / Reason</th></tr>
                </thead>
                <tbody>"""

    for name, desc in org_constraints["active_constraints"]:
        org_constraints_section += f"""
                    <tr class="pass">
                        <td><span class="badge pass">&#10003; Active</span></td>
                        <td><code>{html.escape(name)}</code></td>
                        <td>{html.escape(desc)}</td>
                    </tr>"""

    for name, reason in org_constraints["skipped_constraints"]:
        org_constraints_section += f"""
                    <tr class="uncertain">
                        <td><span class="badge uncertain">&#9888; Skipped</span></td>
                        <td><code>{html.escape(name)}</code></td>
                        <td>{html.escape(reason)}</td>
                    </tr>"""

    org_constraints_section += """
                </tbody>
            </table>"""

    if cai_comparison:
        recovery_sign = "+" if cai_comparison["recovery"] >= 0 else ""
        org_constraints_section += f"""
            <h3>CAI Comparison</h3>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{cai_comparison['cai_all_constraints']:.4f}</div>
                    <div class="stat-label">CAI with all constraints</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{cai_comparison['cai_organism_aware']:.4f}</div>
                    <div class="stat-label">CAI with organism-aware constraints</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{recovery_sign}{cai_comparison['recovery']:.4f}</div>
                    <div class="stat-label">Recovery</div>
                </div>
            </div>"""

    org_constraints_section += """
        </section>"""

    # Engine analysis section (using unified BaseEngineResult fields)
    engine_section = ""
    if engine_results:
        engine_rows = ""
        for r in engine_results:
            status_class = "pass" if r.success else "fail"
            status_text = "PASS" if r.success else "FAIL"
            score_label = r.primary_score_label or "Score"
            engine_label = r.engine_name or "Unknown"
            error_html = f'<div class="violation">{html.escape(r.error)}</div>' if r.error else ""
            engine_rows += f"""
                <tr class="{status_class}">
                    <td><span class="badge {status_class}">{status_text}</span></td>
                    <td>{html.escape(engine_label)}</td>
                    <td>{html.escape(score_label)}</td>
                    <td>{r.primary_score:.4f}</td>
                    <td>{html.escape(r.classification)}</td>
                    <td>{r.execution_time_s:.3f}s</td>
                    <td>{error_html}</td>
                </tr>"""
        engine_section = f"""
        <section class="section">
            <h2>Engine Analysis Results</h2>
            <table>
                <thead>
                    <tr><th>Status</th><th>Engine</th><th>Metric</th><th>Score</th><th>Classification</th><th>Time</th><th>Details</th></tr>
                </thead>
                <tbody>
                    {engine_rows}
                </tbody>
            </table>
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
.verdict-badge.likely_pass {{ background: #22c55e; color: white; }}
.verdict-badge.fail {{ background: var(--fail); color: white; }}
.verdict-badge.likely_fail {{ background: #f87171; color: white; }}
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
tr.likely_pass {{ background: #f0fdf4; }}
tr.fail {{ background: #fef2f2; }}
tr.likely_fail {{ background: #fef2f2; }}
tr.uncertain {{ background: #fffbeb; }}
.badge {{
    display: inline-block;
    padding: 0.15rem 0.5rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.9rem;
}}
.badge.pass {{ background: #dcfce7; color: var(--pass); }}
.badge.likely_pass {{ background: #dcfce7; color: #16a34a; }}
.badge.fail {{ background: #fee2e2; color: var(--fail); }}
.badge.likely_fail {{ background: #fee2e2; color: #dc2626; }}
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
/* === Interactive Feature Styles === */
/* Tooltip */
.tooltip {{
    position: fixed;
    background: #1e293b;
    color: #f8fafc;
    padding: 0.5rem 0.75rem;
    border-radius: 6px;
    font-size: 0.8rem;
    line-height: 1.4;
    max-width: 350px;
    z-index: 9999;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.15s ease;
    box-shadow: 0 4px 12px rgba(0,0,0,0.3);
}}
.tooltip.visible {{ opacity: 1; }}
.tooltip::after {{
    content: '';
    position: absolute;
    bottom: -5px;
    left: 50%;
    transform: translateX(-50%);
    border-width: 5px 5px 0;
    border-style: solid;
    border-color: #1e293b transparent transparent;
}}
/* Stat card hover */
.stat-card[data-tooltip] {{
    cursor: help;
    transition: box-shadow 0.2s ease, transform 0.2s ease;
}}
.stat-card[data-tooltip]:hover {{
    box-shadow: 0 2px 8px rgba(37,99,235,0.15);
    transform: translateY(-1px);
}}
/* Verdict badge tooltip */
.badge[data-tooltip] {{ cursor: help; }}
/* Token type cell tooltip */
.token-type-cell[data-tooltip] {{
    cursor: help;
    border-bottom: 1px dashed var(--text-secondary);
}}
/* Filter buttons */
.filter-bar {{
    display: flex;
    gap: 0.5rem;
    margin-bottom: 1rem;
    flex-wrap: wrap;
}}
.filter-btn {{
    padding: 0.35rem 0.85rem;
    border: 1px solid var(--border);
    border-radius: 6px;
    background: var(--surface);
    color: var(--text);
    font-size: 0.85rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.15s ease;
}}
.filter-btn:hover {{
    border-color: var(--accent);
    color: var(--accent);
}}
.filter-btn.active {{
    background: var(--accent);
    color: white;
    border-color: var(--accent);
}}
.filter-btn[data-filter="pass"].active {{ background: var(--pass); border-color: var(--pass); }}
.filter-btn[data-filter="likely_pass"].active {{ background: #22c55e; border-color: #22c55e; }}
.filter-btn[data-filter="fail"].active {{ background: var(--fail); border-color: var(--fail); }}
.filter-btn[data-filter="likely_fail"].active {{ background: #f87171; border-color: #f87171; }}
.filter-btn[data-filter="uncertain"].active {{ background: var(--uncertain); border-color: var(--uncertain); }}
/* Search input */
.search-bar {{
    margin-bottom: 1rem;
}}
.search-input {{
    width: 100%;
    max-width: 400px;
    padding: 0.5rem 0.75rem;
    border: 1px solid var(--border);
    border-radius: 6px;
    font-size: 0.9rem;
    background: var(--surface);
    color: var(--text);
    outline: none;
    transition: border-color 0.15s ease;
}}
.search-input:focus {{
    border-color: var(--accent);
    box-shadow: 0 0 0 3px rgba(37,99,235,0.1);
}}
.search-input::placeholder {{
    color: var(--text-secondary);
}}
/* Copy button */
.copy-btn {{
    font-size: 0.75rem;
    padding: 0.25rem 0.6rem;
    border: 1px solid var(--border);
    border-radius: 4px;
    background: var(--surface);
    color: var(--text-secondary);
    cursor: pointer;
    margin-left: 0.75rem;
    vertical-align: middle;
    transition: all 0.15s ease;
}}
.copy-btn:hover {{
    border-color: var(--accent);
    color: var(--accent);
}}
.copy-btn.copied {{
    background: var(--pass);
    color: white;
    border-color: var(--pass);
}}
/* Collapsible sections */
.section h2 {{
    cursor: pointer;
    user-select: none;
    position: relative;
    padding-right: 2rem;
}}
.section h2::after {{
    content: '−';
    position: absolute;
    right: 0;
    top: 50%;
    transform: translateY(-50%);
    font-size: 1.5rem;
    font-weight: 300;
    color: var(--text-secondary);
    transition: transform 0.2s ease;
}}
.section h2.collapsed::after {{
    content: '+';
}}
.section-content {{
    overflow: hidden;
    transition: max-height 0.3s ease;
    max-height: 5000px;
}}
.section-content.collapsed {{
    max-height: 0;
}}
/* Dark mode toggle */
.toggle-btn {{
    background: rgba(255,255,255,0.15);
    border: 1px solid rgba(255,255,255,0.3);
    color: white;
    font-size: 1.25rem;
    padding: 0.35rem 0.65rem;
    border-radius: 6px;
    cursor: pointer;
    transition: all 0.2s ease;
    line-height: 1;
}}
.toggle-btn:hover {{
    background: rgba(255,255,255,0.25);
}}
/* Dark mode variables */
body.dark-mode {{
    --bg: #1e293b;
    --surface: #0f172a;
    --border: #334155;
    --text: #f1f5f9;
    --text-secondary: #94a3b8;
    --accent: #60a5fa;
}}
body.dark-mode .header {{
    background: linear-gradient(135deg, #0f172a, #1e293b);
}}
body.dark-mode tr.pass {{ background: #052e16; }}
body.dark-mode tr.likely_pass {{ background: #052e16; }}
body.dark-mode tr.fail {{ background: #450a0a; }}
body.dark-mode tr.likely_fail {{ background: #450a0a; }}
body.dark-mode tr.uncertain {{ background: #451a03; }}
body.dark-mode .badge.pass {{ background: #052e16; color: #4ade80; }}
body.dark-mode .badge.likely_pass {{ background: #052e16; color: #4ade80; }}
body.dark-mode .badge.fail {{ background: #450a0a; color: #f87171; }}
body.dark-mode .badge.likely_fail {{ background: #450a0a; color: #f87171; }}
body.dark-mode .badge.uncertain {{ background: #451a03; color: #fbbf24; }}
body.dark-mode .violation {{ color: #f87171; }}
body.dark-mode .knowledge-gap {{ color: #fbbf24; }}
body.dark-mode .filter-btn {{
    background: #334155;
    border-color: #475569;
    color: #e2e8f0;
}}
body.dark-mode .filter-btn:hover {{
    border-color: var(--accent);
    color: var(--accent);
}}
body.dark-mode .search-input {{
    background: #334155;
    border-color: #475569;
    color: #f1f5f9;
}}
body.dark-mode .copy-btn {{
    background: #334155;
    border-color: #475569;
    color: #94a3b8;
}}
body.dark-mode .copy-btn:hover {{
    border-color: var(--accent);
    color: var(--accent);
}}
body.dark-mode .sequence-display {{
    background: #0f172a;
}}
body.dark-mode .stat-card {{
    background: #0f172a;
    border-color: #334155;
}}
body.dark-mode .tooltip {{
    background: #f8fafc;
    color: #1e293b;
}}
body.dark-mode .tooltip::after {{
    border-color: #f8fafc transparent transparent;
}}
@media print {{
    body {{ padding: 0; }}
    .section {{ break-inside: avoid; }}
    .filter-bar, .search-bar, .copy-btn, .toggle-btn {{ display: none; }}
    .section-content {{ max-height: none !important; }}
    .section h2 {{ cursor: default; }}
}}
</style>
</head>
<body>

<div class="header">
    <div>
        <h1>BioCompiler Report</h1>
        <div class="subtitle">
            {html.escape(gene_name or 'Designed Sequence')} &mdash;
            {html.escape(organism)} ({html.escape(org_constraints['domain'])}) &mdash;
            {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
        </div>
    </div>
    <div style="display:flex;align-items:center;gap:1rem;">
        <button id="dark-mode-toggle" class="toggle-btn" title="Toggle dark mode">&#9790;</button>
        <div class="verdict-badge {overall_class}">{overall_verdict}</div>
    </div>
</div>

<!-- Summary Statistics -->
<section class="section">
    <h2>Sequence Summary</h2>
    <div class="stats-grid">
        <div class="stat-card" data-tooltip="Total number of nucleotide base pairs in the sequence">
            <div class="stat-value">{len(seq):,}</div>
            <div class="stat-label">Length (bp)</div>
        </div>
        <div class="stat-card" data-tooltip="GC Content = (G+C)/total, measured over the full sequence">
            <div class="stat-value">{gc:.1%}</div>
            <div class="stat-label">GC Content</div>
        </div>
        <div class="stat-card" data-tooltip="Codon Adaptation Index — geometric mean of relative codon usage, compared to {html.escape(organism)} reference">
            <div class="stat-value">{cai:.4f}</div>
            <div class="stat-label">CAI</div>
        </div>
        <div class="stat-card" data-tooltip="Number of amino acids in the translated protein product">
            <div class="stat-value">{len(protein)}</div>
            <div class="stat-label">Protein (aa)</div>
        </div>
        <div class="stat-card" data-tooltip="Number of exon regions in the gene structure">
            <div class="stat-value">{len(exon_boundaries)}</div>
            <div class="stat-label">Exons</div>
        </div>
        <div class="stat-card" data-tooltip="Number of biological elements detected by the sequence scanner">
            <div class="stat-value">{len(tokens)}</div>
            <div class="stat-label">Tokens</div>
        </div>
    </div>
</section>

{opt_section}

{org_constraints_section}

{engine_section}

<!-- Predicate Results -->
<section class="section">
    <h2>Type-Check Results</h2>
    <div class="filter-bar">
        <button class="filter-btn active" data-filter="all">All</button>
        <button class="filter-btn" data-filter="pass">&#10003; PASS</button>
        <button class="filter-btn" data-filter="fail">&#10007; FAIL</button>
        <button class="filter-btn" data-filter="likely_pass">&#10003;~ LIKELY_PASS</button>
        <button class="filter-btn" data-filter="uncertain">? UNCERTAIN</button>
        <button class="filter-btn" data-filter="likely_fail">&#10007;~ LIKELY_FAIL</button>
    </div>
    <table id="predicate-table">
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
    <div class="search-bar">
        <input type="text" id="token-search" placeholder="Search tokens by type or sequence..." class="search-input" />
    </div>
    <table id="token-table">
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
    <h2>Sequence <button id="copy-seq-btn" class="copy-btn" title="Copy sequence to clipboard">Copy to Clipboard</button></h2>
    <div class="sequence-display" id="sequence-display">{_format_sequence_html(seq, exon_boundaries)}</div>
    <textarea id="raw-sequence" style="position:absolute;left:-9999px;">{html.escape(seq)}</textarea>
</section>

<!-- Protein -->
<section class="section">
    <h2>Protein Translation</h2>
    <div class="sequence-display">{html.escape(protein)}</div>
</section>

<div class="footer">
    Generated by BioCompiler v{__version__} &mdash; Machine-Verified Gene Design
</div>

<script>
(function() {{
    'use strict';

    // === Tooltip System ===
    var tooltipEl = document.createElement('div');
    tooltipEl.className = 'tooltip';
    document.body.appendChild(tooltipEl);

    function showTooltip(e) {{
        var target = e.currentTarget;
        var text = target.getAttribute('data-tooltip');
        if (!text) return;
        tooltipEl.textContent = text;
        tooltipEl.classList.add('visible');
        positionTooltip(e);
    }}

    function positionTooltip(e) {{
        var rect = tooltipEl.getBoundingClientRect();
        var x = e.clientX - rect.width / 2;
        var y = e.clientY - rect.height - 12;
        if (x < 8) x = 8;
        if (x + rect.width > window.innerWidth - 8) x = window.innerWidth - rect.width - 8;
        if (y < 8) y = e.clientY + 18;
        tooltipEl.style.left = x + 'px';
        tooltipEl.style.top = y + 'px';
    }}

    function hideTooltip() {{
        tooltipEl.classList.remove('visible');
    }}

    // Attach tooltip events to all elements with data-tooltip
    var tooltipTargets = document.querySelectorAll('[data-tooltip]');
    for (var i = 0; i < tooltipTargets.length; i++) {{
        tooltipTargets[i].addEventListener('mouseenter', showTooltip);
        tooltipTargets[i].addEventListener('mousemove', positionTooltip);
        tooltipTargets[i].addEventListener('mouseleave', hideTooltip);
    }}

    // === Collapsible Sections ===
    var sections = document.querySelectorAll('.section');
    for (var i = 0; i < sections.length; i++) {{
        var section = sections[i];
        var h2 = section.querySelector('h2');
        if (!h2) continue;

        // Wrap all content after h2 in a section-content div
        var wrapper = document.createElement('div');
        wrapper.className = 'section-content';
        var sibling = h2.nextSibling;
        while (sibling) {{
            var next = sibling.nextSibling;
            wrapper.appendChild(sibling);
            sibling = next;
        }}
        section.appendChild(wrapper);

        // Click handler
        (function(header, content) {{
            header.addEventListener('click', function(e) {{
                // Don't collapse if clicking the copy button
                if (e.target.id === 'copy-seq-btn' || e.target.closest('.copy-btn')) return;
                header.classList.toggle('collapsed');
                content.classList.toggle('collapsed');
            }});
        }})(h2, wrapper);
    }}

    // === Predicate Filter Buttons ===
    var filterBtns = document.querySelectorAll('.filter-btn');
    var predicateTable = document.getElementById('predicate-table');
    if (predicateTable) {{
        for (var i = 0; i < filterBtns.length; i++) {{
            (function(btn) {{
                btn.addEventListener('click', function() {{
                    // Update active state
                    for (var j = 0; j < filterBtns.length; j++) {{
                        filterBtns[j].classList.remove('active');
                    }}
                    btn.classList.add('active');

                    // Filter rows
                    var filter = btn.getAttribute('data-filter');
                    var tbody = predicateTable.querySelector('tbody');
                    var rows = tbody.querySelectorAll('tr');
                    for (var k = 0; k < rows.length; k++) {{
                        var row = rows[k];
                        var verdict = row.getAttribute('data-verdict');
                        if (filter === 'all' || verdict === filter) {{
                            row.style.display = '';
                        }} else {{
                            row.style.display = 'none';
                        }}
                    }}
                }});
            }})(filterBtns[i]);
        }}
    }}

    // === Token Search/Filter ===
    var tokenSearch = document.getElementById('token-search');
    var tokenTable = document.getElementById('token-table');
    if (tokenSearch && tokenTable) {{
        tokenSearch.addEventListener('input', function() {{
            var query = this.value.toLowerCase().trim();
            var tbody = tokenTable.querySelector('tbody');
            var rows = tbody.querySelectorAll('tr');
            for (var k = 0; k < rows.length; k++) {{
                var row = rows[k];
                var elemType = (row.getAttribute('data-element-type') || '').toLowerCase();
                var matchSeq = (row.getAttribute('data-match-sequence') || '').toLowerCase();
                if (!query || elemType.indexOf(query) !== -1 || matchSeq.indexOf(query) !== -1) {{
                    row.style.display = '';
                }} else {{
                    row.style.display = 'none';
                }}
            }}
        }});
    }}

    // === Copy Sequence Button ===
    var copyBtn = document.getElementById('copy-seq-btn');
    var rawSeq = document.getElementById('raw-sequence');
    if (copyBtn && rawSeq) {{
        copyBtn.addEventListener('click', function() {{
            var text = rawSeq.value;
            if (navigator.clipboard && navigator.clipboard.writeText) {{
                navigator.clipboard.writeText(text).then(function() {{
                    showCopied(copyBtn);
                }});
            }} else {{
                rawSeq.style.position = 'fixed';
                rawSeq.style.left = '0';
                rawSeq.style.top = '0';
                rawSeq.style.opacity = '0';
                rawSeq.select();
                rawSeq.setSelectionRange(0, text.length);
                try {{
                    document.execCommand('copy');
                    showCopied(copyBtn);
                }} catch(err) {{}}
                rawSeq.style.position = 'absolute';
                rawSeq.style.left = '-9999px';
                rawSeq.style.opacity = '';
            }}
        }});
    }}

    function showCopied(btn) {{
        var orig = btn.textContent;
        btn.textContent = 'Copied!';
        btn.classList.add('copied');
        setTimeout(function() {{
            btn.textContent = orig;
            btn.classList.remove('copied');
        }}, 2000);
    }}

    // === Dark Mode Toggle ===
    var darkToggle = document.getElementById('dark-mode-toggle');
    if (darkToggle) {{
        // Check for saved preference
        var saved = localStorage.getItem('biocompiler-dark-mode');
        if (saved === 'true') {{
            document.body.classList.add('dark-mode');
            darkToggle.innerHTML = '&#9788;';
        }}

        darkToggle.addEventListener('click', function() {{
            document.body.classList.toggle('dark-mode');
            var isDark = document.body.classList.contains('dark-mode');
            darkToggle.innerHTML = isDark ? '&#9788;' : '&#9790;';
            try {{
                localStorage.setItem('biocompiler-dark-mode', isDark);
            }} catch(e) {{}}
        }});
    }}

}})();
</script>

</body>
</html>"""


def _format_sequence_html(seq: str, exon_boundaries: list[tuple[int, int]]) -> str:
    """Format sequence with exon highlighting.

    Args:
        seq: DNA sequence to format.
        exon_boundaries: List of (start, end) exon positions.

    Returns:
        HTML string with exon bases highlighted and line breaks.
    """
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
        if i > 0 and i % _SEQUENCE_LINE_WIDTH == 0:
            formatted += "\n"
        formatted += char
    return formatted


def _generate_gc_plot_svg(data: list[tuple[int, float]], avg_gc: float) -> str:
    """Generate SVG line chart for GC content.

    Args:
        data: List of (position, gc_fraction) tuples from sliding window.
        avg_gc: Overall average GC fraction for the reference line.

    Returns:
        SVG markup string, or a fallback <p> element when data is empty.
    """
    if not data:
        return "<p>No GC data available</p>"

    width = _GC_PLOT_WIDTH
    height = _GC_PLOT_HEIGHT
    margin = _GC_PLOT_MARGIN
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

    # Threshold lines
    y30 = margin + plot_h - _GC_LOW_THRESHOLD * plot_h
    y70 = margin + plot_h - _GC_HIGH_THRESHOLD * plot_h

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
    """Generate SVG codon usage heatmap.

    Args:
        codon_usage: Nested dict mapping amino-acid → {codon → frequency}.

    Returns:
        SVG markup string, or a fallback <p> element when data is empty.
    """
    if not codon_usage:
        return "<p>No codon usage data</p>"

    aas = sorted(codon_usage.keys())
    cell_w = _HEATMAP_CELL_W
    cell_h = _HEATMAP_CELL_H
    margin_left = _HEATMAP_MARGIN_LEFT
    margin_top = _HEATMAP_MARGIN_TOP
    width = margin_left + _HEATMAP_MAX_CODON_COLS * cell_w + 20
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
    """Generate SVG diagram of exon/intron structure.

    Args:
        seq: DNA sequence (used for length scaling).
        exon_boundaries: List of (start, end) exon positions.
        tokens: Scan tokens (splice donor/acceptor markers are drawn).

    Returns:
        SVG markup string.
    """
    width = _STRUCT_SVG_WIDTH
    height = _STRUCT_SVG_HEIGHT
    margin = _STRUCT_SVG_MARGIN
    bar_y = _STRUCT_SVG_BAR_Y
    bar_h = _STRUCT_SVG_BAR_H

    seq_len = max(len(seq), 1)
    scale = (width - 2 * margin) / seq_len

    elements = ""

    # Draw full sequence bar (introns)
    elements += f'<rect x="{margin}" y="{bar_y}" width="{width - 2 * margin}" height="{bar_h}" fill="#e2e8f0" rx="3"/>'

    # Draw exons
    colors = _STRUCT_EXON_COLORS
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
    for frac in _POSITION_TICK_FRACTIONS:
        pos = int(frac * seq_len)
        x = margin + frac * (width - 2 * margin)
        elements += f'<text x="{x:.1f}" y="{height - 5}" font-size="9" fill="#64748b" text-anchor="middle">{pos}</text>'

    return f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    {elements}
</svg>"""
