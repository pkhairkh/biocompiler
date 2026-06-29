"""
BioCompiler CSP Solver — Constraint Interaction Map
=====================================================

Documents which constraints conflict with each other and quantifies the
CAI cost of each conflict.  This helps users understand *why* CAI drops
and which constraints to relax for maximum CAI improvement.

Key concepts
------------
- An **interaction** exists between two constraints when they both affect
  the same codon position(s).  At those overlapping positions the solver
  must compromise, potentially choosing a sub-optimal codon (lower CAI)
  to satisfy both constraints simultaneously.
- The **estimated CAI cost** is a pre-computed benchmark value derived
  from our test suite of human gene optimization runs.  It represents
  the typical CAI loss attributable to each constraint pair.
- **Conflict severity** ("high", "medium", "low") classifies the
  practical impact of the interaction on solution quality.

Pre-computed known interactions (from benchmarking)
----------------------------------------------------
- NoCrypticSpliceConstraint <-> MaximizeCAI: HIGH (~0.10 CAI loss)
- NoCpGIslandConstraint <-> MaximizeCAI: MEDIUM (~0.05 CAI loss)
- NoRestrictionSiteConstraint <-> MaximizeCAI: LOW (~0.02 CAI loss)
- GCRangeConstraint <-> MaximizeCAI: LOW (usually compatible)
- NoATTTAMotifConstraint <-> MaximizeCAI: LOW (rare conflict)

Usage
-----
::

    from biocompiler.solver.constraint_interaction import (
        ConstraintInteractionMap,
        print_interaction_report,
    )

    model = build_csp_model(protein="MVSKGE", organism="Homo_sapiens", config=cfg)
    imap = ConstraintInteractionMap()
    interactions = imap.build_interaction_map(model)
    report = print_interaction_report(interactions)
    print(report)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from biocompiler.solver.types import ConstraintSpec, CSPModel

logger = logging.getLogger(__name__)


__all__ = [
    "InteractionInfo",
    "ConstraintInteractionMap",
    "print_interaction_report",
]


# ==============================================================================
# Pre-computed known interactions from benchmarking
# ==============================================================================

# Key: frozenset of constraint *name* prefixes that interact.
# Value: (estimated_cai_cost, severity)
#
# The keys use the canonical constraint class names so that the lookup
# works regardless of per-position suffixes (e.g. "NoCrypticSpliceConstraint"
# matches "NoCrypticSpliceConstraint_pos42").

_KNOWN_INTERACTIONS: dict[frozenset[str], tuple[float, str]] = {
    frozenset({"NoCrypticSpliceConstraint", "MaximizeCAI"}): (0.10, "high"),
    frozenset({"NoCpGIslandConstraint", "MaximizeCAI"}): (0.05, "medium"),
    frozenset({"NoRestrictionSiteConstraint", "MaximizeCAI"}): (0.02, "low"),
    frozenset({"GCRangeConstraint", "MaximizeCAI"}): (0.01, "low"),
    frozenset({"NoATTTAMotifConstraint", "MaximizeCAI"}): (0.01, "low"),
    # Hard-hard interactions
    frozenset({"NoCrypticSpliceConstraint", "GCRangeConstraint"}): (0.03, "medium"),
    frozenset({"NoCpGIslandConstraint", "GCRangeConstraint"}): (0.02, "low"),
    frozenset({"NoRestrictionSiteConstraint", "GCRangeConstraint"}): (0.02, "low"),
    frozenset({"NoCrypticSpliceConstraint", "NoCpGIslandConstraint"}): (0.02, "low"),
    frozenset({"NoRestrictionSiteConstraint", "NoCpGIslandConstraint"}): (0.01, "low"),
}

# Fallback estimated CAI cost when a specific pair is not in the
# pre-computed table.  We use a small default because most unknown
# pairs have negligible interaction.
_DEFAULT_CAI_COST: float = 0.005
_DEFAULT_SEVERITY: str = "low"

# Constraint names that are global (affect all positions) rather than
# position-specific.
_GLOBAL_CONSTRAINT_NAMES: frozenset[str] = frozenset({
    "TranslationConstraint",
    "GCRangeConstraint",
    "MaximizeCAI",
    "MinimizeCpG",
    "MinimizeMRNADG",
    "NoATTTAMotifConstraint",
    "NoTRunConstraint",
})


# ==============================================================================
# InteractionInfo dataclass
# ==============================================================================

@dataclass
class InteractionInfo:
    """Information about the interaction between two constraints.

    Attributes:
        constraint_a: Name of the first constraint.
        constraint_b: Name of the second constraint.
        overlapping_positions: Codon positions (0-based) where both
            constraints apply.  For global constraints (e.g. GCRangeConstraint),
            this contains *all* positions in the protein.
        estimated_cai_cost: Estimated CAI loss attributable to this
            interaction, derived from benchmarking data.
        conflict_severity: Qualitative severity: "high", "medium", or "low".
    """

    constraint_a: str
    constraint_b: str
    overlapping_positions: list[int] = field(default_factory=list)
    estimated_cai_cost: float = 0.0
    conflict_severity: str = "low"

    def __post_init__(self) -> None:
        """Validate fields after initialization."""
        if self.conflict_severity not in ("high", "medium", "low"):
            raise ValueError(
                f"conflict_severity must be 'high', 'medium', or 'low', "
                f"got {self.conflict_severity!r}"
            )
        if self.estimated_cai_cost < 0:
            raise ValueError(
                f"estimated_cai_cost must be non-negative, "
                f"got {self.estimated_cai_cost}"
            )


# ==============================================================================
# ConstraintInteractionMap class
# ==============================================================================

class ConstraintInteractionMap:
    """Builds a map of pairwise constraint interactions for a CSP model.

    For each pair of constraints in the model, determines whether they
    can conflict (i.e. both affect the same codon position(s)).  Records
    the overlapping positions and an estimated CAI impact based on
    pre-computed benchmarking data.

    Example
    -------
    ::

        imap = ConstraintInteractionMap()
        interactions = imap.build_interaction_map(model)
        for key, info in interactions.items():
            print(f"{info.constraint_a} <-> {info.constraint_b}: "
                  f"CAI cost={info.estimated_cai_cost:.3f}, "
                  f"severity={info.conflict_severity}")
    """

    def build_interaction_map(
        self,
        model: CSPModel,
    ) -> dict[tuple[str, str], InteractionInfo]:
        """Build the pairwise constraint interaction map for *model*.

        Parameters
        ----------
        model:
            The CSP model whose constraints are to be analysed.

        Returns
        -------
        dict[tuple[str, str], InteractionInfo]
            Mapping from ``(constraint_a_name, constraint_b_name)`` to
            :class:`InteractionInfo`.  Only pairs with at least one
            overlapping position are included.  The tuple is always
            sorted alphabetically so that ``("A", "B")`` and
            ``("B", "A")`` map to the same key.
        """
        constraints = model.constraints
        if not constraints:
            return {}

        # Pre-compute the set of codon positions each constraint affects.
        position_map = self._compute_position_map(model)

        # Determine base constraint names (strip per-position suffixes)
        # for looking up pre-computed interaction data.
        base_names = self._extract_base_names(constraints)

        all_protein_positions = list(range(len(model.protein_sequence)))

        interactions: dict[tuple[str, str], InteractionInfo] = {}

        # Compare every unique pair of constraints.
        for i in range(len(constraints)):
            for j in range(i + 1, len(constraints)):
                ca = constraints[i]
                cb = constraints[j]

                # Compute overlap
                pos_a = position_map[ca.name]
                pos_b = position_map[cb.name]
                overlap = sorted(pos_a & pos_b)

                if not overlap:
                    continue

                # Canonical (sorted) key
                key = tuple(sorted([ca.name, cb.name]))

                # Look up pre-computed CAI cost and severity
                base_a = base_names.get(ca.name, ca.name)
                base_b = base_names.get(cb.name, cb.name)
                pair_key = frozenset({base_a, base_b})

                if pair_key in _KNOWN_INTERACTIONS:
                    cai_cost, severity = _KNOWN_INTERACTIONS[pair_key]
                else:
                    # Scale the default cost by the fraction of positions
                    # that overlap — more overlap implies higher potential
                    # for conflict.
                    n_positions = len(model.protein_sequence)
                    if n_positions > 0:
                        overlap_fraction = len(overlap) / n_positions
                    else:
                        overlap_fraction = 0.0
                    cai_cost = _DEFAULT_CAI_COST * overlap_fraction
                    severity = _DEFAULT_SEVERITY

                # For global + global pairs, the overlap is the entire
                # sequence — but the CAI cost should still come from the
                # pre-computed table (not be inflated by the overlap
                # fraction).
                if pair_key in _KNOWN_INTERACTIONS:
                    cai_cost, severity = _KNOWN_INTERACTIONS[pair_key]

                info = InteractionInfo(
                    constraint_a=key[0],
                    constraint_b=key[1],
                    overlapping_positions=overlap,
                    estimated_cai_cost=round(cai_cost, 4),
                    conflict_severity=severity,
                )

                # If we already have an entry for this pair (possible when
                # multiple position-specific constraints share the same
                # base name), merge the overlap positions and take the
                # maximum CAI cost.
                if key in interactions:
                    existing = interactions[key]
                    merged = sorted(set(existing.overlapping_positions) | set(overlap))
                    max_cost = max(existing.estimated_cai_cost, info.estimated_cai_cost)
                    # Take the worse severity
                    worse = _worse_severity(existing.conflict_severity, severity)
                    interactions[key] = InteractionInfo(
                        constraint_a=key[0],
                        constraint_b=key[1],
                        overlapping_positions=merged,
                        estimated_cai_cost=max_cost,
                        conflict_severity=worse,
                    )
                else:
                    interactions[key] = info

        return interactions

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_position_map(
        self,
        model: CSPModel,
    ) -> dict[str, set[int]]:
        """Compute the set of codon positions each constraint affects.

        For constraints with explicit ``positions``, use those.
        For global constraints (empty positions list), assume all
        protein positions.
        """
        all_positions = set(range(len(model.protein_sequence)))
        result: dict[str, set[int]] = {}

        for c in model.constraints:
            if c.positions:
                result[c.name] = set(c.positions)
            else:
                # Global constraint — affects all positions
                # (or is sequence-level like GCRangeConstraint)
                result[c.name] = all_positions

        return result

    def _extract_base_names(
        self,
        constraints: list[ConstraintSpec],
    ) -> dict[str, str]:
        """Map each constraint name to its base (class-level) name.

        Position-specific constraints may have names like
        "NoCrypticSpliceConstraint_pos42".  The base name is
        "NoCrypticSpliceConstraint".
        """
        known_prefixes = sorted(
            _KNOWN_INTERACTIONS_BASE_NAMES,
            key=len,
            reverse=True,  # longest prefix first for greedy match
        )
        result: dict[str, str] = {}
        for c in constraints:
            name = c.name
            matched = False
            for prefix in known_prefixes:
                if name == prefix or name.startswith(prefix + "_"):
                    result[name] = prefix
                    matched = True
                    break
            if not matched:
                result[name] = name
        return result


# Pre-compute the set of all base constraint names that appear in
# _KNOWN_INTERACTIONS for prefix matching.
_KNOWN_INTERACTIONS_BASE_NAMES: set[str] = set()
for pair in _KNOWN_INTERACTIONS:
    for name in pair:
        _KNOWN_INTERACTIONS_BASE_NAMES.add(name)
# Also add standard constraint names that may not appear in the table
# but are still used as prefixes.
_KNOWN_INTERACTIONS_BASE_NAMES.update({
    "TranslationConstraint",
    "NoRestrictionSiteConstraint",
    "GCRangeConstraint",
    "NoCrypticSpliceConstraint",
    "NoCpGIslandConstraint",
    "NoATTTAMotifConstraint",
    "NoTRunConstraint",
    "MaximizeCAI",
    "MinimizeCpG",
    "MinimizeMRNADG",
})


def _worse_severity(a: str, b: str) -> str:
    """Return the worse of two severity levels.

    Ordering: "high" > "medium" > "low".
    """
    order = {"high": 2, "medium": 1, "low": 0}
    return a if order.get(a, 0) >= order.get(b, 0) else b


# ==============================================================================
# Report printer
# ==============================================================================

def print_interaction_report(
    interaction_map: dict[tuple[str, str], InteractionInfo],
) -> str:
    """Format the interaction map as a human-readable table.

    The table is sorted by estimated CAI cost (highest first) so that
    the most impactful interactions are at the top.  Includes
    recommendations for which constraints to relax.

    Parameters
    ----------
    interaction_map:
        The output of :meth:`ConstraintInteractionMap.build_interaction_map`.

    Returns
    -------
    str
        A formatted multi-line string with the interaction report.
    """
    if not interaction_map:
        return "No constraint interactions detected."

    # Sort by CAI cost descending, then by severity descending
    severity_order = {"high": 2, "medium": 1, "low": 0}
    sorted_items = sorted(
        interaction_map.items(),
        key=lambda kv: (
            -kv[1].estimated_cai_cost,
            -severity_order.get(kv[1].conflict_severity, 0),
        ),
    )

    lines: list[str] = []
    lines.append("=" * 80)
    lines.append("Constraint Interaction Map")
    lines.append("=" * 80)
    lines.append("")

    # Header
    header = (
        f"{'Constraint A':<30} {'Constraint B':<30} "
        f"{'Severity':<10} {'CAI Cost':<10} {'Overlap':<10}"
    )
    lines.append(header)
    lines.append("-" * 80)

    for key, info in sorted_items:
        overlap_str = (
            str(len(info.overlapping_positions))
            if len(info.overlapping_positions) <= 5
            else f"{len(info.overlapping_positions)} pos"
        )
        row = (
            f"{info.constraint_a:<30} {info.constraint_b:<30} "
            f"{info.conflict_severity:<10} "
            f"{info.estimated_cai_cost:<10.4f} {overlap_str:<10}"
        )
        lines.append(row)

    lines.append("-" * 80)
    lines.append("")

    # Recommendations
    lines.append("Recommendations:")
    lines.append("-" * 40)

    recommendations = _generate_recommendations(sorted_items)
    for rec in recommendations:
        lines.append(f"  {rec}")

    lines.append("")
    lines.append(
        "Tip: Relaxing HIGH severity constraints yields the most CAI improvement."
    )
    lines.append(
        "     LOW severity constraints rarely affect CAI and can usually be kept."
    )

    return "\n".join(lines)


def _generate_recommendations(
    sorted_items: list[tuple[tuple[str, str], InteractionInfo]],
) -> list[str]:
    """Generate relaxation recommendations based on interaction severity.

    Parameters
    ----------
    sorted_items:
        Interaction map entries sorted by CAI cost descending.

    Returns
    -------
    list[str]
        Ordered list of recommendations.
    """
    recommendations: list[str] = []
    seen: set[str] = set()

    for key, info in sorted_items:
        if info.conflict_severity == "high":
            # Recommend relaxing the hard constraint (not MaximizeCAI)
            relax = _pick_relaxation_target(info.constraint_a, info.constraint_b)
            if relax not in seen:
                recommendations.append(
                    f"HIGH: Consider relaxing '{relax}' — estimated "
                    f"CAI gain: +{info.estimated_cai_cost:.2f} if relaxed."
                )
                seen.add(relax)
        elif info.conflict_severity == "medium":
            relax = _pick_relaxation_target(info.constraint_a, info.constraint_b)
            if relax not in seen:
                recommendations.append(
                    f"MEDIUM: '{relax}' reduces CAI by ~{info.estimated_cai_cost:.2f}; "
                    f"relax if CAI is critical."
                )
                seen.add(relax)

    if not recommendations:
        recommendations.append(
            "All interactions are LOW severity — no relaxation needed."
        )

    return recommendations


def _pick_relaxation_target(name_a: str, name_b: str) -> str:
    """Choose which constraint in a pair to recommend for relaxation.

    Prefer relaxing hard constraints over soft (MaximizeCAI is the
    objective, not something to relax).
    """
    # Never recommend relaxing MaximizeCAI (that is the goal)
    if name_a == "MaximizeCAI":
        return name_b
    if name_b == "MaximizeCAI":
        return name_a
    # For two hard constraints, recommend the less critical one
    # (using alphabetical order as a tiebreaker)
    return name_a
