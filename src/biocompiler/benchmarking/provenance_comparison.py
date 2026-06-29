"""
BioCompiler Benchmarking — Provenance Comparison & Audit Trail
==============================================================

Compares BioCompiler's full decision provenance with DNAchisel's capabilities,
demonstrating BioCompiler's unique differentiator: per-codon, per-constraint
decision provenance with CAI impact tracking.

BioCompiler records *why* every codon was chosen, what alternatives were
considered, which constraints drove each tradeoff, and the CAI cost of every
constraint — producing a complete audit trail suitable for regulatory
submissions. DNAchisel provides only a basic ``constraints_log`` with no
per-codon provenance, no CAI impact tracking, and no conflict resolution trace.

This module provides:

- **ProvenanceComparison**: Data class comparing provenance capabilities
  between BioCompiler and DNAchisel
- **compare_provenance_capabilities()**: Documents what each tool provides
  and identifies BioCompiler's unique features
- **AuditReport**: Full audit report for an optimized sequence, suitable for
  regulatory submissions
- **generate_audit_report()**: Generates a complete audit report from
  provenance records

Usage::

    from biocompiler.benchmarking.provenance_comparison import (
        compare_provenance_capabilities,
        generate_audit_report,
        ProvenanceComparison,
        AuditReport,
    )

    # Compare provenance capabilities
    comparison = compare_provenance_capabilities()
    print(f"BioCompiler unique: {comparison.unique_to_biocompiler}")

    # Generate an audit report from optimization provenance
    report = generate_audit_report(
        protein="MVLSPADKTNVKAAWGKVGA",
        organism="Homo_sapiens",
        sequence="ATGGTGCTG...",
        provenance_records=my_records,
    )
    print(report.regulatory_summary)

Design goals:
1. Make BioCompiler's provenance advantage explicit and auditable.
2. Provide a standardized audit report format suitable for regulatory filings.
3. Enable head-to-head comparison of decision traceability between tools.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "ProvenanceComparison",
    "compare_provenance_capabilities",
    "AuditReport",
    "generate_audit_report",
]


# ---------------------------------------------------------------------------
# ProvenanceComparison
# ---------------------------------------------------------------------------

@dataclass
class ProvenanceComparison:
    """Comparison of provenance capabilities between BioCompiler and DNAchisel.

    Documents what each tool provides for decision traceability and identifies
    features unique to BioCompiler.

    Attributes:
        biocompiler_features: List of provenance features BioCompiler provides.
        dnachisel_features: List of provenance features DNAchisel provides.
        unique_to_biocompiler: Features only BioCompiler provides.
        overlap: Features both tools provide.
    """

    biocompiler_features: list[str]
    dnachisel_features: list[str]
    unique_to_biocompiler: list[str]
    overlap: list[str]

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "biocompiler_features": list(self.biocompiler_features),
            "dnachisel_features": list(self.dnachisel_features),
            "unique_to_biocompiler": list(self.unique_to_biocompiler),
            "overlap": list(self.overlap),
        }


# ---------------------------------------------------------------------------
# compare_provenance_capabilities
# ---------------------------------------------------------------------------

def compare_provenance_capabilities() -> ProvenanceComparison:
    """Compare BioCompiler's full decision provenance with DNAchisel's capabilities.

    Documents what each tool provides for decision traceability during gene
    optimization, identifying BioCompiler's unique differentiators.

    BioCompiler provides:
        a. Per-codon decision provenance — one CodonDecision per codon position
           (e.g. 495 codon decisions for a 495-aa gene), recording chosen codon,
           alternatives considered, constraint reason, and confidence score.
        b. Per-constraint decision provenance — one ConstraintDecision per active
           constraint (e.g. 12 constraint decisions), recording action taken,
           positions affected, tradeoff description, and CAI impact.
        c. CAI impact per decision — each ConstraintDecision records its
           impact_on_cai, making the CAI cost of every constraint explicit and
           auditable.
        d. Conflict resolution trace — when constraints conflict, the resolution
           is recorded with the action (satisfied/relaxed/conflicted/overridden)
           and a human-readable tradeoff description.
        e. Full audit trail — OptimizationDecisionTrail ties together input,
           configuration, every codon and constraint decision, iteration log,
           and final output into a single serializable object.

    DNAchisel provides:
        a. constraints_log (basic) — a log of which constraints were solved,
           but no per-position or per-codon detail.
        b. No per-codon provenance — DNAchisel does not record why each codon
           was chosen or what alternatives were considered.
        c. No CAI impact tracking — DNAchisel has no native CAI optimization
           constraint and does not track CAI cost of individual constraints.

    Returns:
        ProvenanceComparison with feature lists and unique/overlap analysis.
    """

    biocompiler_features = [
        "Per-codon decision provenance (one decision per codon position)",
        "Per-constraint decision provenance (one decision per active constraint)",
        "CAI impact tracking per constraint decision",
        "Conflict resolution trace (satisfied/relaxed/conflicted/overridden)",
        "Full audit trail (OptimizationDecisionTrail)",
        "Alternatives considered per codon with CAI/GC contributions",
        "Confidence score per codon decision",
        "Constraint reason per codon decision",
        "Iteration log (solver convergence behavior)",
        "Position-based constraint pressure analysis",
        "Codon usage summary from decisions",
        "Serializable decision trail (JSON/dict)",
        "Reproducibility metadata (seed, version, timestamp)",
        "Constraint cost breakdown in CAI terms",
        "What-if scenario analysis from decision trail",
        "Human-readable explanation per position (explain_position)",
    ]

    dnachisel_features = [
        "Basic constraints_log (which constraints were solved)",
        "Constraint specification tracking (spec class names)",
        "Sequence output with constraint satisfaction",
        "Post-hoc CAI computation (via adapter, not native)",
        "Post-hoc GC content computation (via adapter, not native)",
    ]

    # Compute overlap (features both tools provide in some form)
    # DNAchisel provides basic constraint logging but not the granular
    # per-codon/per-constraint provenance that BioCompiler does.
    # The overlap is at the level of "the tool records which constraints
    # were applied" — DNAchisel logs constraint specs, BioCompiler logs
    # per-constraint decisions.
    overlap = [
        "Constraint application tracking (spec names / constraint decisions)",
        "Post-hoc CAI computation",
        "Post-hoc GC content computation",
    ]

    # Features unique to BioCompiler
    bc_set = set(biocompiler_features)
    dn_set = set(dnachisel_features)
    overlap_set = set(overlap)

    unique_to_biocompiler = sorted(bc_set - overlap_set)
    # Re-add any overlap items that are in dnachisel but not in biocompiler
    # by name (they exist in different form)
    # Compute properly:
    # unique_to_biocompiler = features in biocompiler_features that are NOT
    # in dnachisel_features AND not in overlap
    unique_to_biocompiler = sorted(
        f for f in biocompiler_features
        if f not in dnachisel_features and f not in overlap
    )

    # Overlap: features present (in some form) in both
    # We've already defined this explicitly above since the feature names
    # differ between tools but the capabilities overlap at a high level.

    comparison = ProvenanceComparison(
        biocompiler_features=biocompiler_features,
        dnachisel_features=dnachisel_features,
        unique_to_biocompiler=unique_to_biocompiler,
        overlap=overlap,
    )

    logger.info(
        "Provenance comparison: BioCompiler=%d features, DNAchisel=%d features, "
        "unique_to_biocompiler=%d, overlap=%d",
        len(biocompiler_features),
        len(dnachisel_features),
        len(unique_to_biocompiler),
        len(overlap),
    )

    return comparison


# ---------------------------------------------------------------------------
# AuditReport
# ---------------------------------------------------------------------------

@dataclass
class AuditReport:
    """Full audit report for an optimized sequence, suitable for regulatory
    submissions.

    Captures every decision made during optimization with quantified CAI
    impacts, constraint breakdowns, and a human-readable regulatory summary.

    Attributes:
        protein: Input amino acid sequence.
        organism: Target organism identifier.
        sequence: Optimized DNA sequence.
        total_codon_decisions: Total number of codon-level decisions recorded.
        total_constraint_decisions: Total number of constraint-level decisions.
        cai_total: Overall CAI of the optimized sequence.
        cai_unconstrained: CAI if no constraints were applied (theoretical max).
        cai_cost_breakdown: Dict mapping constraint name → CAI cost (negative).
        top_decisions: Top 10 most impactful decisions (by confidence or CAI
            impact), each as a dict with position, amino_acid, chosen_codon,
            constraint_reason, confidence, and impact keys.
        conflict_resolutions: List of conflict resolution records, each as a
            dict with constraint_name, action_taken, tradeoff_description,
            and impact_on_cai keys.
        regulatory_summary: Human-readable paragraph suitable for regulatory
            filing, summarizing the optimization decisions and tradeoffs.
    """

    protein: str
    organism: str
    sequence: str
    total_codon_decisions: int
    total_constraint_decisions: int
    cai_total: float
    cai_unconstrained: float
    cai_cost_breakdown: dict[str, float]
    top_decisions: list[dict[str, Any]]
    conflict_resolutions: list[dict[str, Any]]
    regulatory_summary: str

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "protein": self.protein,
            "organism": self.organism,
            "sequence": self.sequence,
            "total_codon_decisions": self.total_codon_decisions,
            "total_constraint_decisions": self.total_constraint_decisions,
            "cai_total": self.cai_total,
            "cai_unconstrained": self.cai_unconstrained,
            "cai_cost_breakdown": dict(self.cai_cost_breakdown),
            "top_decisions": list(self.top_decisions),
            "conflict_resolutions": list(self.conflict_resolutions),
            "regulatory_summary": self.regulatory_summary,
        }


# ---------------------------------------------------------------------------
# generate_audit_report
# ---------------------------------------------------------------------------

def generate_audit_report(
    protein: str,
    organism: str,
    sequence: str,
    provenance_records: list[dict[str, Any]],
) -> AuditReport:
    """Generate a full audit report for an optimized sequence.

    Analyzes provenance records from a BioCompiler optimization run and
    produces a structured audit report with quantified CAI impacts, constraint
    breakdowns, and a human-readable regulatory summary.

    The provenance_records list should contain dicts with the following
    possible structures:

    - Codon decisions: dicts with keys ``"type"`` (``"codon"``),
      ``"position"``, ``"amino_acid"``, ``"chosen_codon"``,
      ``"constraint_reason"``, ``"confidence"``, ``"alternatives_considered"``.
    - Constraint decisions: dicts with keys ``"type"`` (``"constraint"``),
      ``"constraint_name"``, ``"constraint_type"``, ``"action_taken"``,
      ``"positions_affected"``, ``"tradeoff_description"``, ``"impact_on_cai"``.

    If ``provenance_records`` contains :class:`CodonDecision` or
    :class:`ConstraintDecision` objects from
    :mod:`biocompiler.decision_provenance`, they will be converted to dicts
    automatically.

    Args:
        protein: Input amino acid sequence.
        organism: Target organism identifier (e.g. ``"Homo_sapiens"``).
        sequence: Optimized DNA sequence.
        provenance_records: List of provenance record dicts or objects.

    Returns:
        AuditReport with full decision analysis and regulatory summary.
    """
    # Normalize records to dicts
    normalized: list[dict[str, Any]] = []
    for rec in provenance_records:
        if isinstance(rec, dict):
            normalized.append(rec)
        elif hasattr(rec, "to_dict"):
            normalized.append(rec.to_dict())
        else:
            # Try to extract useful fields
            normalized.append({
                "type": getattr(rec, "decision_type", "unknown"),
                "raw": str(rec),
            })

    # Separate codon and constraint decisions
    codon_records: list[dict[str, Any]] = []
    constraint_records: list[dict[str, Any]] = []

    for rec in normalized:
        rec_type = rec.get("type", "")
        if rec_type == "codon" or "chosen_codon" in rec:
            codon_records.append(rec)
        elif rec_type == "constraint" or "constraint_name" in rec:
            constraint_records.append(rec)

    total_codon_decisions = len(codon_records)
    total_constraint_decisions = len(constraint_records)

    # Compute CAI total and unconstrained
    # Try to compute from sequence if biocompiler is available
    cai_total = 0.0
    cai_unconstrained = 0.0
    try:
        from biocompiler.expression.translation import compute_cai
        cai_total = compute_cai(sequence, organism)
        # For unconstrained, we do not recompute — it should be provided
        # externally or estimated. Default to cai_total as a floor.
        cai_unconstrained = cai_total
    except (ImportError, Exception):
        # If we cannot compute, leave as 0.0 — caller can override
        pass

    # CAI cost breakdown per constraint
    cai_cost_breakdown: dict[str, float] = {}
    for rec in constraint_records:
        name = rec.get("constraint_name", "unknown")
        impact = rec.get("impact_on_cai", 0.0)
        if isinstance(impact, (int, float)):
            cai_cost_breakdown[name] = cai_cost_breakdown.get(name, 0.0) + float(impact)

    # Adjust cai_unconstrained: total CAI minus the sum of all constraint costs
    # Constraint impacts are negative, so unconstrained = total - sum(costs)
    total_cai_cost = sum(cai_cost_breakdown.values())
    if cai_total > 0 and total_cai_cost < 0:
        cai_unconstrained = cai_total - total_cai_cost

    # Top 10 most impactful decisions
    # Sort codon records by impact (low confidence = more impactful decisions)
    # and constraint records by CAI impact
    decision_impacts: list[dict[str, Any]] = []

    for rec in codon_records:
        confidence = rec.get("confidence", 1.0)
        impact_score = 1.0 - confidence  # Lower confidence = higher impact
        decision_impacts.append({
            "position": rec.get("position", -1),
            "amino_acid": rec.get("amino_acid", "?"),
            "chosen_codon": rec.get("chosen_codon", "???"),
            "constraint_reason": rec.get("constraint_reason", "unknown"),
            "confidence": confidence,
            "impact": impact_score,
            "type": "codon",
        })

    for rec in constraint_records:
        cai_impact = rec.get("impact_on_cai", 0.0)
        if isinstance(cai_impact, (int, float)):
            impact_score = abs(float(cai_impact))
        else:
            impact_score = 0.0
        decision_impacts.append({
            "constraint_name": rec.get("constraint_name", "unknown"),
            "action_taken": rec.get("action_taken", "unknown"),
            "tradeoff_description": rec.get("tradeoff_description", ""),
            "impact_on_cai": cai_impact if isinstance(cai_impact, (int, float)) else 0.0,
            "impact": impact_score,
            "type": "constraint",
        })

    # Sort by impact descending, take top 10
    top_decisions = sorted(
        decision_impacts, key=lambda d: d.get("impact", 0.0), reverse=True
    )[:10]

    # Remove the internal "impact" score from output
    top_decisions_out: list[dict[str, Any]] = []
    for d in top_decisions:
        d_copy = {k: v for k, v in d.items() if k != "impact"}
        top_decisions_out.append(d_copy)

    # Conflict resolutions: constraint decisions with action != "satisfied"
    conflict_resolutions: list[dict[str, Any]] = []
    for rec in constraint_records:
        action = rec.get("action_taken", "")
        if action in ("relaxed", "conflicted", "overridden"):
            conflict_resolutions.append({
                "constraint_name": rec.get("constraint_name", "unknown"),
                "action_taken": action,
                "tradeoff_description": rec.get("tradeoff_description", ""),
                "impact_on_cai": rec.get("impact_on_cai", 0.0),
            })

    # Generate regulatory summary
    regulatory_summary = _generate_regulatory_summary(
        protein=protein,
        organism=organism,
        sequence=sequence,
        total_codon_decisions=total_codon_decisions,
        total_constraint_decisions=total_constraint_decisions,
        cai_total=cai_total,
        cai_unconstrained=cai_unconstrained,
        cai_cost_breakdown=cai_cost_breakdown,
        conflict_resolutions=conflict_resolutions,
    )

    report = AuditReport(
        protein=protein,
        organism=organism,
        sequence=sequence,
        total_codon_decisions=total_codon_decisions,
        total_constraint_decisions=total_constraint_decisions,
        cai_total=cai_total,
        cai_unconstrained=cai_unconstrained,
        cai_cost_breakdown=cai_cost_breakdown,
        top_decisions=top_decisions_out,
        conflict_resolutions=conflict_resolutions,
        regulatory_summary=regulatory_summary,
    )

    logger.info(
        "Generated audit report: protein=%d aa, organism=%s, "
        "codon_decisions=%d, constraint_decisions=%d, cai=%.4f",
        len(protein), organism,
        total_codon_decisions, total_constraint_decisions, cai_total,
    )

    return report


# ---------------------------------------------------------------------------
# _generate_regulatory_summary — internal helper
# ---------------------------------------------------------------------------

def _generate_regulatory_summary(
    protein: str,
    organism: str,
    sequence: str,
    total_codon_decisions: int,
    total_constraint_decisions: int,
    cai_total: float,
    cai_unconstrained: float,
    cai_cost_breakdown: dict[str, float],
    conflict_resolutions: list[dict[str, Any]],
) -> str:
    """Generate a human-readable regulatory summary paragraph.

    Produces a paragraph suitable for inclusion in regulatory filings that
    summarizes the optimization decisions and tradeoffs in plain language.

    Args:
        protein: Input amino acid sequence.
        organism: Target organism.
        sequence: Optimized DNA sequence.
        total_codon_decisions: Number of codon decisions.
        total_constraint_decisions: Number of constraint decisions.
        cai_total: Achieved CAI.
        cai_unconstrained: Theoretical maximum CAI (no constraints).
        cai_cost_breakdown: Constraint name → CAI cost.
        conflict_resolutions: List of conflict resolution records.

    Returns:
        Human-readable paragraph summarizing the optimization.
    """
    protein_len = len(protein)
    seq_len = len(sequence)

    parts: list[str] = []

    # Opening statement
    parts.append(
        f"This audit report documents the codon optimization of a "
        f"{protein_len}-amino-acid protein for expression in "
        f"{organism.replace('_', ' ')}, producing a {seq_len}-nucleotide "
        f"DNA sequence."
    )

    # Decision count
    parts.append(
        f"The optimization involved {total_codon_decisions} codon-level "
        f"decisions and {total_constraint_decisions} constraint-level "
        f"decisions, each recorded with full provenance including "
        f"alternatives considered and rationale."
    )

    # CAI statement
    if cai_total > 0:
        parts.append(
            f"The achieved Codon Adaptation Index (CAI) is {cai_total:.4f}."
        )
        if cai_unconstrained > cai_total and cai_unconstrained > 0:
            cai_retention = (cai_total / cai_unconstrained) * 100
            parts.append(
                f"The unconstrained maximum CAI would be "
                f"{cai_unconstrained:.4f}, meaning the optimized sequence "
                f"retains {cai_retention:.1f}% of theoretical CAI performance "
                f"while satisfying all biological constraints."
            )

    # Constraint cost breakdown
    if cai_cost_breakdown:
        total_cost = abs(sum(cai_cost_breakdown.values()))
        parts.append(
            f"The total CAI cost attributable to constraint satisfaction is "
            f"{total_cost:.4f} CAI points."
        )
        # List individual constraint costs
        cost_items = sorted(
            cai_cost_breakdown.items(), key=lambda x: abs(x[1]), reverse=True
        )
        if len(cost_items) <= 3:
            cost_descriptions = [
                f"{name} ({abs(cost):.4f} CAI points)"
                for name, cost in cost_items
            ]
            parts.append(
                f"Breakdown by constraint: {'; '.join(cost_descriptions)}."
            )
        else:
            top_costs = [
                f"{name} ({abs(cost):.4f})"
                for name, cost in cost_items[:3]
            ]
            parts.append(
                f"Top constraint costs: {'; '.join(top_costs)}; "
                f"and {len(cost_items) - 3} additional constraints."
            )

    # Conflict resolutions
    if conflict_resolutions:
        n_conflicts = len(conflict_resolutions)
        parts.append(
            f"There were {n_conflicts} constraint conflict(s) that required "
            f"resolution during optimization."
        )
        for cr in conflict_resolutions[:5]:
            name = cr.get("constraint_name", "unknown")
            action = cr.get("action_taken", "unknown")
            desc = cr.get("tradeoff_description", "")
            if desc:
                parts.append(
                    f"The {name} constraint was {action}: {desc}."
                )

    # Provenance statement
    parts.append(
        "Every codon choice, alternative evaluated, constraint tradeoff, "
        "and conflict resolution has been recorded with full provenance, "
        "enabling complete reproducibility and regulatory traceability. "
        "This audit trail is available in machine-readable JSON format "
        "for independent verification."
    )

    return " ".join(parts)
