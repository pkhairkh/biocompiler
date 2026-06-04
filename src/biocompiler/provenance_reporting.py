"""
BioCompiler Provenance Query & Reporting — Decision Trail Analysis and Visualization
=====================================================================================

Provides query, analysis, and visualization tools for
:class:`~biocompiler.decision_provenance.OptimizationDecisionTrail` objects.

Every gene optimization produces a rich decision trail recording *why* each
codon was chosen, what alternatives were considered, and which constraints
drove each tradeoff.  This module makes that trail **queryable** (via
:class:`ProvenanceQuery`) and **visible** (via :class:`ProvenanceReport` and
:func:`explain_position`), enabling users to understand and audit every
optimization decision.

Key components:

- **ProvenanceQuery** — structured queries over a decision trail (filter by
  amino acid, find sub-optimal positions, quantify constraint costs, etc.)
- **ProvenanceReport** — generate Markdown, HTML, or JSON reports from a trail
- **explain_position** — human-readable explanation of a single codon decision

Usage::

    from biocompiler.provenance_reporting import (
        ProvenanceQuery, ProvenanceReport, explain_position,
    )
    from biocompiler.decision_provenance import OptimizationDecisionTrail

    # Load a trail (e.g., from JSON or produced by DecisionProvenanceCollector)
    trail = OptimizationDecisionTrail.from_json(open("trail.json").read())

    # Query the trail
    query = ProvenanceQuery(trail)
    leucine_decisions = query.decisions_for_amino_acid("L")
    constrained = query.constraints_that_reduced_cai()
    cai_lost = query.cai_lost_to_constraints()
    top_positions = query.most_constrained_positions(n=5)
    usage = query.codon_usage_summary()
    alternatives = query.alternatives_not_chosen()

    # Generate reports
    md_report = ProvenanceReport.generate_markdown(trail)
    html_report = ProvenanceReport.generate_html(trail)
    json_report = ProvenanceReport.generate_json(trail)

    # Explain a specific position
    explanation = explain_position(trail, position=42)
    print(explanation)
"""

from __future__ import annotations

import html as html_module
import json
import logging
from datetime import datetime, timezone
from typing import Any

from .decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
)

logger = logging.getLogger(__name__)

__all__ = [
    "ProvenanceQuery",
    "ProvenanceReport",
    "explain_position",
]


# ---------------------------------------------------------------------------
# Amino acid name mapping for human-readable explanations
# ---------------------------------------------------------------------------

_AMINO_ACID_FULL_NAMES: dict[str, str] = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
}


def _aa_full_name(aa: str) -> str:
    """Return the full name of an amino acid given its 1-letter code.

    Falls back to ``"Amino acid <code>"`` if the code is unknown.
    """
    return _AMINO_ACID_FULL_NAMES.get(aa.upper(), f"Amino acid {aa}")


def _humanize_reason(reason: str) -> str:
    """Convert a machine-readable constraint_reason into a human-readable phrase.

    Examples::

        "maximize_cai" -> "maximize CAI"
        "avoid_restriction_site:EcoRI" -> "avoid the EcoRI restriction site"
        "gc_content" -> "satisfy GC content constraints"
    """
    if reason == "maximize_cai":
        return "maximize CAI"
    if reason.startswith("avoid_restriction_site:"):
        enzyme = reason.split(":", 1)[1]
        return f"avoid the {enzyme} restriction site"
    if reason == "gc_content":
        return "satisfy GC content constraints"
    if reason == "avoid_cpg":
        return "avoid CpG dinucleotides"
    if reason == "no_cryptic_splice":
        return "eliminate cryptic splice sites"
    if reason == "codon_pair_bias":
        return "optimize codon pair bias"
    # Generic fallback: replace underscores with spaces
    return reason.replace("_", " ")


# ---------------------------------------------------------------------------
# ProvenanceQuery
# ---------------------------------------------------------------------------

class ProvenanceQuery:
    """Structured queries over an :class:`OptimizationDecisionTrail`.

    Provides methods to filter codon decisions by amino acid or position,
    identify constraints that forced sub-optimal codon choices, quantify
    CAI losses, and summarize codon usage.

    Example::

        query = ProvenanceQuery(trail)
        print(query.decisions_for_amino_acid("L"))
        print(query.cai_lost_to_constraints())
        print(query.most_constrained_positions(n=5))
    """

    def __init__(self, trail: OptimizationDecisionTrail) -> None:
        if not isinstance(trail, OptimizationDecisionTrail):
            raise TypeError(
                f"Expected OptimizationDecisionTrail, got {type(trail).__name__}"
            )
        self._trail = trail
        # Build position index for O(1) lookup
        self._position_index: dict[int, CodonDecision] = {}
        for cd in trail.codon_decisions:
            self._position_index[cd.position] = cd

    # -- Amino acid filtering ------------------------------------------------

    def decisions_for_amino_acid(self, aa: str) -> list[CodonDecision]:
        """Return all codon decisions for a given amino acid.

        Args:
            aa: 1-letter amino acid code (e.g. ``"L"`` for Leucine).

        Returns:
            List of :class:`CodonDecision` instances whose ``amino_acid``
            matches *aa*, in trail order.  Empty list if none found.
        """
        aa_upper = aa.upper()
        return [
            cd for cd in self._trail.codon_decisions
            if cd.amino_acid.upper() == aa_upper
        ]

    # -- Position lookup -----------------------------------------------------

    def decisions_at_position(self, pos: int) -> CodonDecision:
        """Return the codon decision at a specific position.

        Args:
            pos: 0-indexed codon position in the protein sequence.

        Returns:
            The :class:`CodonDecision` at *pos*.

        Raises:
            KeyError: If no decision exists for the given position.
        """
        if pos in self._position_index:
            return self._position_index[pos]
        raise KeyError(
            f"No codon decision found at position {pos}. "
            f"Available positions: {sorted(self._position_index.keys())}"
        )

    # -- Constraint analysis -------------------------------------------------

    def constraints_that_reduced_cai(self) -> list[ConstraintDecision]:
        """Return constraint decisions that forced sub-optimal codon choices.

        A constraint reduced CAI if its ``impact_on_cai`` is negative
        (i.e. it cost CAI to satisfy the constraint).

        Returns:
            List of :class:`ConstraintDecision` instances with
            ``impact_on_cai < 0``, sorted by CAI impact (most negative first).
        """
        costly = [
            cd for cd in self._trail.constraint_decisions
            if cd.impact_on_cai < 0
        ]
        return sorted(costly, key=lambda c: c.impact_on_cai)

    def cai_lost_to_constraints(self) -> float:
        """Return the total CAI points lost due to all constraints.

        Sums the ``impact_on_cai`` of all constraint decisions.  Because
        ``impact_on_cai`` is negative when a constraint costs CAI, the
        returned value is also negative (or zero if no CAI was lost).

        Returns:
            Total CAI loss as a negative float (or 0.0).
        """
        return sum(cd.impact_on_cai for cd in self._trail.constraint_decisions)

    def most_constrained_positions(self, n: int = 5) -> list[tuple[int, str]]:
        """Return positions with the most constraint pressure.

        Constraint pressure at a position is measured by counting how many
        constraint decisions list that position in their
        ``positions_affected``.

        Args:
            n: Number of top positions to return.

        Returns:
            List of ``(position, amino_acid)`` tuples sorted by constraint
            pressure (most constrained first).  Ties are broken by position
            order.
        """
        pressure: dict[int, int] = {}
        for cd in self._trail.constraint_decisions:
            for pos in cd.positions_affected:
                pressure[pos] = pressure.get(pos, 0) + 1

        # Sort by pressure descending, then position ascending for ties
        sorted_positions = sorted(
            pressure.keys(),
            key=lambda p: (-pressure[p], p),
        )

        result: list[tuple[int, str]] = []
        for pos in sorted_positions[:n]:
            aa = "?"
            if pos in self._position_index:
                aa = self._position_index[pos].amino_acid
            result.append((pos, aa))
        return result

    # -- Usage summary -------------------------------------------------------

    def codon_usage_summary(self) -> dict[str, int]:
        """Return a summary of how often each codon was chosen.

        Returns:
            Dict mapping codon strings to their selection count, sorted by
            count descending.
        """
        usage: dict[str, int] = {}
        for cd in self._trail.codon_decisions:
            usage[cd.chosen_codon] = usage.get(cd.chosen_codon, 0) + 1
        return dict(sorted(usage.items(), key=lambda x: x[1], reverse=True))

    def alternatives_not_chosen(self) -> dict[str, list[str]]:
        """Return codons available but not chosen at each position.

        For each codon position, lists the alternative codons that were
        considered but rejected.  The ``alternatives_considered`` field on
        each :class:`CodonDecision` is a list of dicts; each dict should
        contain a ``"codon"`` key (falling back to string representation).

        Returns:
            Dict mapping ``"{position}:{amino_acid}"`` to a list of codon
            strings that were available but not chosen.
        """
        result: dict[str, list[str]] = {}
        for cd in self._trail.codon_decisions:
            alt_codons: list[str] = []
            for alt in cd.alternatives_considered:
                if isinstance(alt, dict) and "codon" in alt:
                    alt_codons.append(str(alt["codon"]))
                elif isinstance(alt, str):
                    alt_codons.append(alt)
                else:
                    alt_codons.append(str(alt))
            key = f"{cd.position}:{cd.amino_acid}"
            result[key] = alt_codons
        return result


# ---------------------------------------------------------------------------
# ProvenanceReport
# ---------------------------------------------------------------------------

class ProvenanceReport:
    """Generate formatted reports from an :class:`OptimizationDecisionTrail`.

    Supports three output formats:

    - **Markdown** — human-readable, suitable for documentation or review
    - **HTML** — styled, with syntax highlighting for code and tables
    - **JSON** — complete machine-readable export

    Example::

        md = ProvenanceReport.generate_markdown(trail)
        html = ProvenanceReport.generate_html(trail)
        json_str = ProvenanceReport.generate_json(trail)
    """

    @staticmethod
    def generate_markdown(trail: OptimizationDecisionTrail) -> str:
        """Generate a Markdown report from a decision trail.

        Sections:
        1. Executive summary (CAI, GC, constraints, key tradeoffs)
        2. Constraint impact table (constraint -> CAI impact)
        3. Top 5 positions where constraints forced sub-optimal choices
        4. Codon usage comparison (chosen vs available)
        5. What-if scenarios summary

        Args:
            trail: The optimization decision trail to report on.

        Returns:
            Multi-line Markdown string.
        """
        query = ProvenanceQuery(trail)
        lines: list[str] = []

        # -- 1. Executive summary --------------------------------------------
        lines.append("# BioCompiler Provenance Report")
        lines.append("")
        lines.append(f"**Gene:** {trail.gene_name or '(unnamed)'}")
        lines.append(f"**Organism:** {trail.organism}")
        lines.append(f"**Solver backend:** {trail.solver_backend}")
        lines.append(f"**Seed:** {trail.seed}")
        lines.append(f"**Timestamp:** {trail.timestamp}")
        lines.append(f"**Version:** {trail.version}")
        lines.append("")

        lines.append("## Executive Summary")
        lines.append("")
        lines.append(f"- **CAI:** {trail.total_cai:.4f}")
        lines.append(f"- **GC content:** {trail.total_gc:.3f}")
        lines.append(f"- **Protein length:** {len(trail.input_protein)} aa")
        lines.append(f"- **Codon decisions:** {len(trail.codon_decisions)}")
        lines.append(f"- **Constraint decisions:** {len(trail.constraint_decisions)}")

        cai_lost = query.cai_lost_to_constraints()
        constraints_reducing = query.constraints_that_reduced_cai()
        lines.append(
            f"- **Total CAI lost to constraints:** "
            f"{abs(cai_lost):.4f} ({len(constraints_reducing)} constraints)"
        )

        # CAI impact from codon decisions (per-position CAI cost)
        codon_cai_cost = sum(
            cd.cai_impact for cd in trail.codon_decisions if cd.cai_impact < 0
        )
        constrained_positions = [
            cd for cd in trail.codon_decisions if cd.cai_impact < 0
        ]
        lines.append(
            f"- **Total CAI cost of constraint satisfaction:** "
            f"{abs(codon_cai_cost):.4f} ({len(constrained_positions)} positions)"
        )

        # Key tradeoffs
        if constraints_reducing:
            lines.append("")
            lines.append("**Key tradeoffs:**")
            for cd in constraints_reducing[:3]:
                lines.append(
                    f"  - {cd.constraint_name}: {cd.tradeoff_description} "
                    f"(CAI impact: {cd.impact_on_cai:+.4f})"
                )
        lines.append("")

        # -- 2. Constraint impact table --------------------------------------
        lines.append("## Constraint Impact")
        lines.append("")
        if trail.constraint_decisions:
            lines.append(
                "| Constraint | Type | Action | CAI Impact | Positions Affected |"
            )
            lines.append(
                "|------------|------|--------|------------|--------------------|"
            )
            for cd in trail.constraint_decisions:
                positions_str = (
                    ", ".join(str(p) for p in cd.positions_affected[:5])
                )
                if len(cd.positions_affected) > 5:
                    positions_str += (
                        f", ... (+{len(cd.positions_affected) - 5})"
                    )
                lines.append(
                    f"| {cd.constraint_name} | {cd.constraint_type} | "
                    f"{cd.action_taken} | {cd.impact_on_cai:+.4f} | "
                    f"{positions_str or '(none)'} |"
                )
        else:
            lines.append("*No constraint decisions recorded.*")
        lines.append("")

        # -- 3. Top 5 most constrained positions ----------------------------
        lines.append("## Most Constrained Positions (Top 5)")
        lines.append("")

        # Top 5 most CAI-expensive constraint fixes
        lines.append("## Top 5 Most CAI-Expensive Constraint Fixes")
        lines.append("")
        if constrained_positions:
            # Sort by CAI impact (most negative = most expensive)
            sorted_by_cost = sorted(
                constrained_positions, key=lambda cd: cd.cai_impact
            )
            lines.append(
                "| Position | Amino Acid | Chosen Codon | Constraint Reason | CAI Impact |"
            )
            lines.append(
                "|----------|------------|-------------|-------------------|------------|"
            )
            for cd in sorted_by_cost[:5]:
                lines.append(
                    f"| {cd.position} | {cd.amino_acid} | {cd.chosen_codon} | "
                    f"{_humanize_reason(cd.constraint_reason)} | "
                    f"{cd.cai_impact:+.4f} |"
                )
        else:
            lines.append("*No positions incurred CAI cost due to constraints.*")
        lines.append("")

        top_positions = query.most_constrained_positions(n=5)
        if top_positions:
            lines.append("| Position | Amino Acid | Full Name | Constraint Pressure |")
            lines.append("|----------|------------|-----------|---------------------|")
            # Compute pressure for display
            pressure_map: dict[int, int] = {}
            for cd in trail.constraint_decisions:
                for pos in cd.positions_affected:
                    pressure_map[pos] = pressure_map.get(pos, 0) + 1

            for pos, aa in top_positions:
                lines.append(
                    f"| {pos} | {aa} | {_aa_full_name(aa)} | "
                    f"{pressure_map.get(pos, 0)} constraints |"
                )
        else:
            lines.append("*No positions affected by constraints.*")
        lines.append("")

        # -- 4. Codon usage comparison ---------------------------------------
        lines.append("## Codon Usage Comparison")
        lines.append("")
        usage = query.codon_usage_summary()

        if usage:
            lines.append("| Codon | Times Chosen | Amino Acid | Available Alternatives |")
            lines.append("|-------|-------------|------------|------------------------|")
            # Group by amino acid for richer display
            aa_for_codon: dict[str, str] = {}
            for cd in trail.codon_decisions:
                aa_for_codon[cd.chosen_codon] = cd.amino_acid

            for codon, count in usage.items():
                aa = aa_for_codon.get(codon, "?")
                # Collect alternatives for this codon across all positions
                alt_set: set[str] = set()
                for cd in trail.codon_decisions:
                    if cd.chosen_codon == codon:
                        for alt in cd.alternatives_considered:
                            if isinstance(alt, dict) and "codon" in alt:
                                alt_set.add(str(alt["codon"]))
                            elif isinstance(alt, str):
                                alt_set.add(alt)
                alt_str = ", ".join(sorted(alt_set)) if alt_set else "(none)"
                lines.append(
                    f"| {codon} | {count} | {aa} | {alt_str} |"
                )
        else:
            lines.append("*No codon decisions recorded.*")
        lines.append("")

        # -- 5. What-if scenarios summary ------------------------------------
        lines.append("## What-If Scenarios Summary")
        lines.append("")
        if trail.iteration_log:
            lines.append(
                f"The optimization completed in "
                f"{len(trail.iteration_log)} iteration(s)."
            )
            lines.append("")
            for i, entry in enumerate(trail.iteration_log):
                action = entry.get("action", "(no action)")
                score = entry.get("score", "N/A")
                lines.append(f"  {i + 1}. action={action}, score={score}")
        else:
            lines.append("*No iteration log recorded.*")

        # Summarize constraint relaxation potential
        if constraints_reducing:
            lines.append("")
            lines.append("### Potential CAI Recovery")
            lines.append("")
            total_recoverable = sum(
                abs(cd.impact_on_cai)
                for cd in constraints_reducing
            )
            lines.append(
                f"If all constraints that reduced CAI were removed, "
                f"up to **{total_recoverable:.4f}** CAI points could "
                f"potentially be recovered."
            )
        lines.append("")

        return "\n".join(lines)

    @staticmethod
    def generate_html(trail: OptimizationDecisionTrail) -> str:
        """Generate an HTML report from a decision trail.

        Same content as :meth:`generate_markdown` but formatted as styled
        HTML with syntax highlighting.

        Args:
            trail: The optimization decision trail to report on.

        Returns:
            Complete HTML document as a string.
        """
        query = ProvenanceQuery(trail)
        cai_lost = query.cai_lost_to_constraints()
        constraints_reducing = query.constraints_that_reduced_cai()
        top_positions = query.most_constrained_positions(n=5)
        usage = query.codon_usage_summary()

        # Build pressure map for top positions display
        pressure_map: dict[int, int] = {}
        for cd in trail.constraint_decisions:
            for pos in cd.positions_affected:
                pressure_map[pos] = pressure_map.get(pos, 0) + 1

        # Escape helper
        def _esc(s: str) -> str:
            return html_module.escape(str(s))

        parts: list[str] = []

        # HTML header
        parts.append("<!DOCTYPE html>")
        parts.append("<html lang='en'>")
        parts.append("<head>")
        parts.append("<meta charset='utf-8'>")
        parts.append(
            f"<title>Provenance Report — "
            f"{_esc(trail.gene_name or 'unnamed')}</title>"
        )
        parts.append("<style>")
        parts.append("body { font-family: -apple-system, BlinkMacSystemFont, "
                      "'Segoe UI', Roboto, sans-serif; margin: 2em auto; "
                      "max-width: 960px; color: #1a1a1a; background: #fafafa; }")
        parts.append("h1 { color: #2c3e50; border-bottom: 2px solid #3498db; "
                      "padding-bottom: 0.3em; }")
        parts.append("h2 { color: #2c3e50; margin-top: 1.5em; }")
        parts.append("h3 { color: #34495e; }")
        parts.append("table { border-collapse: collapse; width: 100%; "
                      "margin: 1em 0; }")
        parts.append("th, td { border: 1px solid #ddd; padding: 8px 12px; "
                      "text-align: left; }")
        parts.append("th { background: #3498db; color: white; }")
        parts.append("tr:nth-child(even) { background: #f2f2f2; }")
        parts.append("tr:hover { background: #e8f4fd; }")
        parts.append("code { background: #f0f0f0; padding: 2px 6px; "
                      "border-radius: 3px; font-size: 0.9em; }")
        parts.append(".metric { font-size: 1.2em; font-weight: bold; "
                      "color: #2c3e50; }")
        parts.append(".negative { color: #e74c3c; }")
        parts.append(".positive { color: #27ae60; }")
        parts.append(".summary-grid { display: grid; grid-template-columns: "
                      "repeat(auto-fit, minmax(200px, 1fr)); gap: 1em; "
                      "margin: 1em 0; }")
        parts.append(".summary-card { background: white; border: 1px solid "
                      "#ddd; border-radius: 8px; padding: 1em; "
                      "text-align: center; }")
        parts.append(".summary-card .label { font-size: 0.85em; "
                      "color: #7f8c8d; }")
        parts.append(".summary-card .value { font-size: 1.5em; "
                      "font-weight: bold; color: #2c3e50; }")
        parts.append("</style>")
        parts.append("</head>")
        parts.append("<body>")

        # Title
        parts.append(
            f"<h1>Provenance Report: "
            f"{_esc(trail.gene_name or 'unnamed')}</h1>"
        )
        parts.append(
            f"<p><strong>Organism:</strong> {_esc(trail.organism)} &nbsp;|&nbsp; "
            f"<strong>Backend:</strong> {_esc(trail.solver_backend)} &nbsp;|&nbsp; "
            f"<strong>Seed:</strong> {_esc(str(trail.seed))} &nbsp;|&nbsp; "
            f"<strong>Timestamp:</strong> {_esc(trail.timestamp)}</p>"
        )

        # Executive summary as cards
        parts.append("<h2>Executive Summary</h2>")
        parts.append("<div class='summary-grid'>")
        parts.append(
            f"<div class='summary-card'><div class='label'>CAI</div>"
            f"<div class='value'>{trail.total_cai:.4f}</div></div>"
        )
        parts.append(
            f"<div class='summary-card'><div class='label'>GC Content</div>"
            f"<div class='value'>{trail.total_gc:.3f}</div></div>"
        )
        parts.append(
            f"<div class='summary-card'><div class='label'>Protein Length</div>"
            f"<div class='value'>{len(trail.input_protein)} aa</div></div>"
        )
        parts.append(
            f"<div class='summary-card'><div class='label'>Codon Decisions</div>"
            f"<div class='value'>{len(trail.codon_decisions)}</div></div>"
        )
        parts.append(
            f"<div class='summary-card'><div class='label'>CAI Lost to Constraints</div>"
            f"<div class='value negative'>{abs(cai_lost):.4f}</div></div>"
        )
        parts.append(
            f"<div class='summary-card'><div class='label'>Constraints Reducing CAI</div>"
            f"<div class='value'>{len(constraints_reducing)}</div></div>"
        )
        parts.append("</div>")

        # CAI cost of constraint satisfaction
        codon_cai_cost = sum(
            cd.cai_impact for cd in trail.codon_decisions if cd.cai_impact < 0
        )
        constrained_positions = [
            cd for cd in trail.codon_decisions if cd.cai_impact < 0
        ]
        parts.append("<div class='summary-grid'>")
        parts.append(
            f"<div class='summary-card'><div class='label'>"
            f"Total CAI Cost of Constraint Satisfaction</div>"
            f"<div class='value negative'>{abs(codon_cai_cost):.4f}</div></div>"
        )
        parts.append(
            f"<div class='summary-card'><div class='label'>"
            f"Positions with CAI Cost</div>"
            f"<div class='value'>{len(constrained_positions)}</div></div>"
        )
        parts.append("</div>")

        # Key tradeoffs
        if constraints_reducing:
            parts.append("<h3>Key Tradeoffs</h3>")
            parts.append("<ul>")
            for cd in constraints_reducing[:3]:
                parts.append(
                    f"<li><strong>{_esc(cd.constraint_name)}:</strong> "
                    f"{_esc(cd.tradeoff_description)} "
                    f"(CAI impact: <span class='negative'>"
                    f"{cd.impact_on_cai:+.4f}</span>)</li>"
                )
            parts.append("</ul>")

        # Constraint impact table
        parts.append("<h2>Constraint Impact</h2>")
        if trail.constraint_decisions:
            parts.append("<table>")
            parts.append(
                "<tr><th>Constraint</th><th>Type</th><th>Action</th>"
                "<th>CAI Impact</th><th>Positions Affected</th></tr>"
            )
            for cd in trail.constraint_decisions:
                positions_str = ", ".join(
                    str(p) for p in cd.positions_affected[:5]
                )
                if len(cd.positions_affected) > 5:
                    positions_str += (
                        f", ... (+{len(cd.positions_affected) - 5})"
                    )
                impact_class = (
                    "negative" if cd.impact_on_cai < 0 else "positive"
                )
                parts.append(
                    f"<tr>"
                    f"<td><code>{_esc(cd.constraint_name)}</code></td>"
                    f"<td>{_esc(cd.constraint_type)}</td>"
                    f"<td>{_esc(cd.action_taken)}</td>"
                    f"<td class='{impact_class}'>"
                    f"{cd.impact_on_cai:+.4f}</td>"
                    f"<td>{_esc(positions_str or '(none)')}</td>"
                    f"</tr>"
                )
            parts.append("</table>")
        else:
            parts.append("<p><em>No constraint decisions recorded.</em></p>")

        # Most constrained positions
        parts.append("<h2>Most Constrained Positions (Top 5)</h2>")

        # Top 5 most CAI-expensive constraint fixes
        parts.append("<h2>Top 5 Most CAI-Expensive Constraint Fixes</h2>")
        if constrained_positions:
            sorted_by_cost = sorted(
                constrained_positions, key=lambda cd: cd.cai_impact
            )
            parts.append("<table>")
            parts.append(
                "<tr><th>Position</th><th>Amino Acid</th>"
                "<th>Chosen Codon</th><th>Constraint Reason</th>"
                "<th>CAI Impact</th></tr>"
            )
            for cd in sorted_by_cost[:5]:
                impact_class = "negative"
                parts.append(
                    f"<tr><td>{cd.position}</td>"
                    f"<td><code>{_esc(cd.amino_acid)}</code></td>"
                    f"<td><code>{_esc(cd.chosen_codon)}</code></td>"
                    f"<td>{_esc(_humanize_reason(cd.constraint_reason))}</td>"
                    f"<td class='{impact_class}'>"
                    f"{cd.cai_impact:+.4f}</td></tr>"
                )
            parts.append("</table>")
        else:
            parts.append(
                "<p><em>No positions incurred CAI cost due to constraints.</em></p>"
            )

        if top_positions:
            parts.append("<table>")
            parts.append(
                "<tr><th>Position</th><th>Amino Acid</th>"
                "<th>Full Name</th><th>Constraint Pressure</th></tr>"
            )
            for pos, aa in top_positions:
                parts.append(
                    f"<tr><td>{pos}</td><td><code>{_esc(aa)}</code></td>"
                    f"<td>{_esc(_aa_full_name(aa))}</td>"
                    f"<td>{pressure_map.get(pos, 0)} constraints</td></tr>"
                )
            parts.append("</table>")
        else:
            parts.append(
                "<p><em>No positions affected by constraints.</em></p>"
            )

        # Codon usage comparison
        parts.append("<h2>Codon Usage Comparison</h2>")
        if usage:
            # Build amino acid lookup for chosen codons
            aa_for_codon: dict[str, str] = {}
            for cd in trail.codon_decisions:
                aa_for_codon[cd.chosen_codon] = cd.amino_acid

            parts.append("<table>")
            parts.append(
                "<tr><th>Codon</th><th>Times Chosen</th>"
                "<th>Amino Acid</th><th>Available Alternatives</th></tr>"
            )
            for codon, count in usage.items():
                aa = aa_for_codon.get(codon, "?")
                alt_set: set[str] = set()
                for cd in trail.codon_decisions:
                    if cd.chosen_codon == codon:
                        for alt in cd.alternatives_considered:
                            if isinstance(alt, dict) and "codon" in alt:
                                alt_set.add(str(alt["codon"]))
                            elif isinstance(alt, str):
                                alt_set.add(alt)
                alt_str = (
                    ", ".join(sorted(alt_set)) if alt_set else "(none)"
                )
                parts.append(
                    f"<tr><td><code>{_esc(codon)}</code></td>"
                    f"<td>{count}</td><td>{_esc(aa)}</td>"
                    f"<td>{_esc(alt_str)}</td></tr>"
                )
            parts.append("</table>")
        else:
            parts.append("<p><em>No codon decisions recorded.</em></p>")

        # What-if scenarios summary
        parts.append("<h2>What-If Scenarios Summary</h2>")
        if trail.iteration_log:
            parts.append(
                f"<p>The optimization completed in "
                f"{len(trail.iteration_log)} iteration(s).</p>"
            )
            parts.append("<ol>")
            for entry in trail.iteration_log:
                action = entry.get("action", "(no action)")
                score = entry.get("score", "N/A")
                parts.append(
                    f"<li>action=<code>{_esc(str(action))}</code>, "
                    f"score=<code>{_esc(str(score))}</code></li>"
                )
            parts.append("</ol>")
        else:
            parts.append("<p><em>No iteration log recorded.</em></p>")

        # Potential CAI recovery
        if constraints_reducing:
            total_recoverable = sum(
                abs(cd.impact_on_cai) for cd in constraints_reducing
            )
            parts.append("<h3>Potential CAI Recovery</h3>")
            parts.append(
                f"<p>If all constraints that reduced CAI were removed, "
                f"up to <strong>{total_recoverable:.4f}</strong> CAI points "
                f"could potentially be recovered.</p>"
            )

        parts.append("</body>")
        parts.append("</html>")

        return "\n".join(parts)

    @staticmethod
    def generate_json(trail: OptimizationDecisionTrail) -> str:
        """Generate a complete machine-readable JSON export.

        Includes the full trail data plus derived analysis metrics.

        Args:
            trail: The optimization decision trail to export.

        Returns:
            Indented JSON string.
        """
        query = ProvenanceQuery(trail)

        report_data: dict[str, Any] = {
            "report_type": "provenance_report",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "trail": trail.to_dict(),
            "analysis": {
                "cai_lost_to_constraints": query.cai_lost_to_constraints(),
                "constraints_reducing_cai": [
                    cd.to_dict()
                    for cd in query.constraints_that_reduced_cai()
                ],
                "most_constrained_positions": [
                    {"position": pos, "amino_acid": aa}
                    for pos, aa in query.most_constrained_positions(n=10)
                ],
                "codon_usage_summary": query.codon_usage_summary(),
                "alternatives_not_chosen": query.alternatives_not_chosen(),
                "cai_impact_summary": {
                    "total_cai_cost_of_constraint_satisfaction": sum(
                        cd.cai_impact for cd in trail.codon_decisions
                        if cd.cai_impact < 0
                    ),
                    "positions_with_cai_cost": len([
                        cd for cd in trail.codon_decisions
                        if cd.cai_impact < 0
                    ]),
                    "top_5_most_expensive_fixes": [
                        {
                            "position": cd.position,
                            "amino_acid": cd.amino_acid,
                            "chosen_codon": cd.chosen_codon,
                            "constraint_reason": cd.constraint_reason,
                            "cai_impact": cd.cai_impact,
                        }
                        for cd in sorted(
                            [cd for cd in trail.codon_decisions if cd.cai_impact < 0],
                            key=lambda cd: cd.cai_impact,
                        )[:5]
                    ],
                },
            },
        }
        return json.dumps(report_data, indent=2, sort_keys=True, default=str)


# ---------------------------------------------------------------------------
# explain_position
# ---------------------------------------------------------------------------

def explain_position(trail: OptimizationDecisionTrail, position: int) -> str:
    """Produce a human-readable explanation of the codon decision at *position*.

    The explanation covers what codon was chosen, why it was chosen (CAI,
    constraints), what alternatives were rejected and why, and whether any
    constraints were violated or satisfied by the chosen codon.

    Example output::

        "At position 42 (Leucine), codon CTG was chosen because it has the
        highest CAI contribution (0.95) among the 6 synonymous codons. TTA
        was rejected because it would lower CAI to 0.14. No constraints were
        violated by CTG."

    Args:
        trail: The optimization decision trail.
        position: 0-indexed codon position to explain.

    Returns:
        Human-readable explanation string.

    Raises:
        KeyError: If no decision exists for the given position.
    """
    query = ProvenanceQuery(trail)
    decision = query.decisions_at_position(position)

    aa_name = _aa_full_name(decision.amino_acid)
    reason = _humanize_reason(decision.constraint_reason)
    sentences: list[str] = []

    # --- Collect CAI and rejection data from alternatives ---
    chosen_cai: float | None = None
    alt_details: list[str] = []
    num_alternatives = len(decision.alternatives_considered)
    total_synonymous = 1 + num_alternatives  # chosen + alternatives

    for alt in decision.alternatives_considered:
        if isinstance(alt, dict):
            alt_codon = alt.get("codon", "?")
            alt_cai = alt.get("cai", None)
            alt_rejection = alt.get("rejection_reason", None)
            alt_violations = alt.get("violations", [])

            if alt_codon == decision.chosen_codon and alt_cai is not None:
                chosen_cai = float(alt_cai)

            rejection_parts: list[str] = []
            if alt_rejection:
                rejection_parts.append(str(alt_rejection))
            if alt_violations:
                if isinstance(alt_violations, list):
                    rejection_parts.extend(str(v) for v in alt_violations)
                else:
                    rejection_parts.append(str(alt_violations))

            if rejection_parts:
                alt_details.append(
                    f"{alt_codon} was rejected because "
                    + "; ".join(rejection_parts)
                )

    # --- Opening sentence ---
    if chosen_cai is not None:
        cai_phrase = (
            f" because it has the highest CAI contribution "
            f"({chosen_cai:.2f})"
        )
        if num_alternatives > 0:
            cai_phrase += f" among the {total_synonymous} synonymous codons"
        sentences.append(
            f"At position {position} ({aa_name}), codon "
            f"{decision.chosen_codon} was chosen{cai_phrase}."
        )
    else:
        sentences.append(
            f"At position {position} ({aa_name}), codon "
            f"{decision.chosen_codon} was chosen to {reason}."
        )

    # --- Confidence ---
    if decision.confidence < 1.0:
        sentences.append(
            f"Confidence in this choice: {decision.confidence:.2f} "
            f"(lower than 1.0 indicates conflicting constraints)."
        )

    # --- Rejected alternatives ---
    if alt_details:
        for detail in alt_details[:5]:  # Limit to top 5
            sentences.append(f"{detail}.")
    elif num_alternatives > 0 and not alt_details:
        # Generic statement about alternatives
        alt_codons: list[str] = []
        for alt in decision.alternatives_considered:
            if isinstance(alt, dict) and "codon" in alt:
                alt_codons.append(str(alt["codon"]))
            elif isinstance(alt, str):
                alt_codons.append(alt)
        if alt_codons:
            sentences.append(
                f"Alternative codons considered but not chosen: "
                f"{', '.join(alt_codons)}."
            )

    # --- Constraints affecting this position ---
    constraints_at_pos = [
        cd for cd in trail.constraint_decisions
        if position in cd.positions_affected
    ]
    if constraints_at_pos:
        constraint_names = [cd.constraint_name for cd in constraints_at_pos]
        sentences.append(
            f"Constraints active at this position: "
            f"{', '.join(constraint_names)}."
        )
        for cd in constraints_at_pos:
            sentences.append(
                f"  - {cd.constraint_name}: {cd.tradeoff_description} "
                f"(CAI impact: {cd.impact_on_cai:+.4f})."
            )
    else:
        sentences.append(
            f"No constraints were violated by {decision.chosen_codon}."
        )

    # --- Original codon comparison (if applicable) ---
    if decision.original_codon and decision.original_codon != decision.chosen_codon:
        sentences.append(
            f"Original codon was {decision.original_codon}, "
            f"changed to {decision.chosen_codon}."
        )

    return " ".join(sentences)
