"""BioCompiler CSP Solver — Constraint Conflict Resolution Provenance.

This module extends the conflict resolution system with **provenance tracking**:
every conflict detected and every resolution applied is recorded as a
:class:`ConflictProvenance` instance so that downstream consumers can
understand *why* a particular constraint was relaxed and what the impact
was.

Provenance is opt-in: when ``track_provenance=False`` (the default) the
resolver behaves identically to the base :class:`ConflictResolver`.  When
enabled, a list of :class:`ConflictProvenance` records is returned
alongside the adjusted sequence.

Public API
----------
- ``ConflictProvenance``           — dataclass recording one conflict resolution
- ``ConflictResolverWithProvenance`` — provenance-aware conflict resolver
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field
from typing import Optional

from biocompiler.solver.types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
    SolverConfig,
)
from .conflict_resolution import (
    ConstraintConflict,
    ConflictResolver,
    prioritize_constraints,
)

logger = logging.getLogger(__name__)


__all__ = [
    "ConflictProvenance",
    "ConflictResolverWithProvenance",
]


# =====================================================================
# ConflictProvenance dataclass
# =====================================================================

_VALID_RESOLUTION_METHODS = frozenset({
    "priority_based",
    "weight_based",
    "manual",
    "csp_backtrack",
    "cai_aware",
})


@dataclass
class ConflictProvenance:
    """A record of one constraint conflict and its resolution.

    Captures the full provenance of a single conflict resolution event:
    which constraints conflicted, how the conflict was resolved, which
    constraint won and which was relaxed, the estimated impact on the
    optimisation objective, and the positions affected.

    Attributes:
        conflicting_constraints: Names of the constraints that conflicted.
        resolution_method: How the conflict was resolved.  One of
            ``"priority_based"``, ``"weight_based"``, ``"manual"``,
            ``"csp_backtrack"``, or ``"cai_aware"``.
        winner: Name of the constraint that was satisfied (won).
        loser: Name of the constraint that was relaxed (lost).
        impact: Human-readable description of the resolution's impact.
        positions_affected: Codon positions where the conflict manifested
            and the resolution was applied.
        cai_impact: Estimated CAI impact of the resolution.  Positive
            values mean the resolution *helped* CAI; negative values
            mean CAI was sacrificed.  Zero means no estimated impact.
            This is a heuristic estimate based on constraint type.
        cai_delta: Actual measured CAI delta of the chosen resolution.
            Unlike ``cai_impact`` (which is heuristic), this is the
            computed difference in CAI contribution at the affected
            position(s).  Positive = CAI improved; negative = CAI lost.
            Defaults to ``None`` when not computed.
        codon_changes: List of ``(position, old_codon, new_codon, cai_delta)``
            tuples recording the actual codon substitutions made to resolve
            this conflict.  Each tuple captures a single position where a
            codon was changed, the original codon, the replacement, and the
            measured CAI delta for that specific change.  Empty list when no
            codon changes were needed or the information was not available.
    """

    conflicting_constraints: list[str]
    resolution_method: str
    winner: str
    loser: str
    impact: str
    positions_affected: list[int]
    cai_impact: float = 0.0
    cai_delta: Optional[float] = None
    codon_changes: list[tuple[int, str, str, float]] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate ``resolution_method`` after initialization."""
        if self.resolution_method not in _VALID_RESOLUTION_METHODS:
            raise ValueError(
                f"Invalid resolution_method {self.resolution_method!r}; "
                f"must be one of {sorted(_VALID_RESOLUTION_METHODS)}"
            )
        # Ensure codon_changes is always a list
        if self.codon_changes is None:
            object.__setattr__(self, 'codon_changes', [])

    def __repr__(self) -> str:
        delta_str = (
            f", cai_delta={self.cai_delta:+.6f}"
            if self.cai_delta is not None
            else ""
        )
        changes_str = f", codon_changes={len(self.codon_changes)}" if self.codon_changes else ""
        return (
            f"ConflictProvenance("
            f"{self.winner!r} > {self.loser!r}, "
            f"method={self.resolution_method!r}, "
            f"positions={self.positions_affected}, "
            f"cai_impact={self.cai_impact:+.4f}{delta_str}{changes_str})"
        )


# =====================================================================
# ConflictResolverWithProvenance
# =====================================================================

class ConflictResolverWithProvenance:
    """Provenance-aware conflict resolver.

    Extends the existing :class:`ConflictResolver` logic with full
    provenance tracking.  Every conflict detected and every resolution
    applied is recorded as a :class:`ConflictProvenance` instance.

    When ``track_provenance=False`` (default) the resolver delegates to
    :class:`ConflictResolver` without recording any provenance, keeping
    behaviour backward-compatible.

    Usage::

        resolver = ConflictResolverWithProvenance(track_provenance=True)
        adjusted_sequence, records = resolver.resolve_conflicts(
            constraints, sequence,
        )
        for rec in records:
            print(f"{rec.winner} beat {rec.loser} via {rec.resolution_method}")

    Parameters
    ----------
    track_provenance : bool
        Whether to record provenance.  Defaults to ``False`` for backward
        compatibility.
    organism : str
        Target organism name (used for CAI impact estimation).
    """

    def __init__(
        self,
        track_provenance: bool = False,
        organism: str = "Homo_sapiens",
        cai_aware: bool = True,
    ) -> None:
        self._track_provenance = track_provenance
        self._organism = organism
        self._cai_aware = cai_aware
        self._base_resolver = ConflictResolver(cai_aware=cai_aware)

        logger.info(
            "ConflictResolverWithProvenance initialized: "
            "track_provenance=%s, cai_aware=%s",
            track_provenance, cai_aware,
        )

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def resolve_conflicts(
        self,
        constraints: list[ConstraintSpec],
        sequence: str,
    ) -> tuple[str, list[ConflictProvenance]]:
        """Resolve conflicts among constraints and record provenance.

        Detects pairwise conflicts among hard constraints, resolves them
        using the priority-based strategy from :class:`ConflictResolver`,
        and returns the adjusted sequence together with a list of
        :class:`ConflictProvenance` records.

        If ``track_provenance`` is ``False``, no provenance is recorded
        and the returned list is always empty (backward compatible).

        Parameters
        ----------
        constraints:
            List of constraint specifications to check for conflicts.
        sequence:
            The current candidate DNA sequence.

        Returns
        -------
        tuple[str, list[ConflictProvenance]]
            A tuple of (adjusted_sequence, provenance_records).
            - adjusted_sequence: The input sequence, possibly modified
              if conflicts required constraint relaxation.  Currently
              the sequence itself is not modified here — the returned
              value is the same as the input, and the provenance records
              indicate which constraints were relaxed so the caller can
              re-solve if needed.
            - provenance_records: List of :class:`ConflictProvenance`
              instances, one per conflict resolution.  Empty if
              ``track_provenance`` is ``False`` or if no conflicts were
              detected.
        """
        provenance_records: list[ConflictProvenance] = []

        if not self._track_provenance:
            return sequence, provenance_records

        # Build a minimal CSPModel to use the base resolver's detect_conflicts
        model = self._build_model(constraints, sequence)

        # Detect conflicts using the base resolver
        conflicts = self._base_resolver.detect_conflicts(model)

        if not conflicts:
            logger.debug("No conflicts detected; returning sequence unchanged.")
            return sequence, provenance_records

        # Build a name -> ConstraintSpec lookup
        constraint_map: dict[str, ConstraintSpec] = {
            c.name: c for c in constraints
        }

        # Resolve each conflict and record provenance
        for conflict in conflicts:
            record = self._resolve_one_conflict(
                conflict, constraint_map, sequence,
            )
            if record is not None:
                provenance_records.append(record)

        logger.info(
            "Recorded %d conflict provenance record(s)",
            len(provenance_records),
        )

        return sequence, provenance_records

    def record_violation_provenance(
        self,
        violations: list,
        sequence: str,
        cai_impact_override: float | None = None,
        codon_changes: list[tuple[int, str, str, float]] | None = None,
    ) -> list[ConflictProvenance]:
        """Record provenance for constraint violations found post-solve.

        When :func:`validate_csp_solution` finds violations, this method
        converts them into :class:`ConflictProvenance` records so that
        every violation's tradeoff is tracked.

        Parameters
        ----------
        violations:
            List of :class:`ConstraintViolation` objects from validation.
        sequence:
            The candidate DNA sequence.
        cai_impact_override:
            If provided, use this CAI impact value for all violation records
            instead of the heuristic estimate.
        codon_changes:
            List of ``(position, old_codon, new_codon, cai_delta)`` tuples
            documenting the actual codon substitutions associated with the
            violations.

        Returns
        -------
        list[ConflictProvenance]
            Provenance records for each violation that represents a
            conflict tradeoff.  Empty if ``track_provenance`` is ``False``.
        """
        if not self._track_provenance or not violations:
            return []

        records: list[ConflictProvenance] = []
        for violation in violations:
            # Each violation represents a constraint that was implicitly
            # relaxed because it conflicted with a higher-priority one.
            # We record this as a csp_backtrack provenance entry.
            cai_value = (
                cai_impact_override
                if cai_impact_override is not None
                else self._estimate_cai_impact(
                    violation.constraint_name, sequence,
                )
            )
            record = ConflictProvenance(
                conflicting_constraints=[violation.constraint_name, "<higher_priority_constraint>"],
                resolution_method="csp_backtrack",
                winner="<higher_priority_constraint>",
                loser=violation.constraint_name,
                impact=(
                    f"Constraint '{violation.constraint_name}' was violated "
                    f"(severity={violation.severity:.2f}, "
                    f"priority={violation.priority.name}) in favor of a "
                    f"higher-priority constraint."
                ),
                positions_affected=list(violation.positions) if violation.positions else [],
                cai_impact=cai_value,
                cai_delta=cai_impact_override,
                codon_changes=codon_changes or [],
            )
            records.append(record)

        logger.info(
            "Recorded %d violation provenance record(s)",
            len(records),
        )
        return records

    def record_cai_aware_provenance(
        self,
        constraint_name: str,
        codon_position: int,
        old_codon: str,
        new_codon: str,
        cai_delta: float,
        positions_affected: list[int] | None = None,
    ) -> ConflictProvenance:
        """Record provenance for a CAI-aware constraint fix.

        When :meth:`ConstraintEnforcer.enforce_with_cai_awareness` applies
        a fix that minimizes CAI loss, this method creates a provenance
        record documenting the tradeoff.

        Parameters
        ----------
        constraint_name:
            Name of the constraint that was fixed.
        codon_position:
            Codon position where the fix was applied.
        old_codon:
            The codon that was replaced.
        new_codon:
            The replacement codon.
        cai_delta:
            The measured CAI delta of the fix (negative = CAI lost).
        positions_affected:
            Nucleotide positions affected by the fix.

        Returns
        -------
        ConflictProvenance
            The provenance record with ``cai_aware`` resolution method
            and the measured ``cai_delta``.
        """
        if positions_affected is None:
            positions_affected = [
                codon_position * 3,
                codon_position * 3 + 1,
                codon_position * 3 + 2,
            ]

        impact_desc = (
            f"Constraint '{constraint_name}' was fixed at codon position "
            f"{codon_position} by replacing {old_codon} with {new_codon}. "
            f"CAI delta: {cai_delta:+.6f}. "
        )
        if cai_delta >= 0:
            impact_desc += "Fix improved or maintained CAI."
        else:
            impact_desc += (
                "Fix sacrificed CAI — this was the least-damaging "
                "alternative among valid codon choices."
            )

        return ConflictProvenance(
            conflicting_constraints=[constraint_name, "MaximizeCAI"],
            resolution_method="cai_aware",
            winner=constraint_name,
            loser="MaximizeCAI",
            impact=impact_desc,
            positions_affected=positions_affected,
            cai_impact=self._estimate_cai_impact(constraint_name, ""),
            cai_delta=cai_delta,
            codon_changes=[(codon_position, old_codon, new_codon, cai_delta)],
        )

    def record_relaxation_provenance(
        self,
        relaxed_constraint_name: str,
        kept_constraint_name: str,
        positions_affected: list[int],
        sequence: str,
        resolution_method: str = "priority_based",
        cai_impact_override: float | None = None,
        codon_changes: list[tuple[int, str, str, float]] | None = None,
    ) -> ConflictProvenance:
        """Record provenance for an explicit constraint relaxation.

        When the solver has to relax a constraint, this creates a
        provenance record documenting the tradeoff.

        Parameters
        ----------
        relaxed_constraint_name:
            Name of the constraint that was relaxed (the loser).
        kept_constraint_name:
            Name of the constraint that was kept (the winner).
        positions_affected:
            Codon positions where the relaxation has effect.
        sequence:
            The candidate DNA sequence.
        resolution_method:
            How the relaxation was decided.  Defaults to ``"priority_based"``.
        cai_impact_override:
            If provided, use this CAI impact value instead of the heuristic
            estimate.  This should be the measured CAI delta from the
            constraint resolution.
        codon_changes:
            List of ``(position, old_codon, new_codon, cai_delta)`` tuples
            documenting the actual codon substitutions associated with the
            relaxation.

        Returns
        -------
        ConflictProvenance
            The provenance record.
        """
        cai_value = (
            cai_impact_override
            if cai_impact_override is not None
            else self._estimate_cai_impact(
                relaxed_constraint_name, sequence,
            )
        )
        return ConflictProvenance(
            conflicting_constraints=[relaxed_constraint_name, kept_constraint_name],
            resolution_method=resolution_method,
            winner=kept_constraint_name,
            loser=relaxed_constraint_name,
            impact=(
                f"Constraint '{relaxed_constraint_name}' was relaxed to "
                f"satisfy '{kept_constraint_name}' at position(s) "
                f"{positions_affected}."
            ),
            positions_affected=positions_affected,
            cai_impact=cai_value,
            cai_delta=cai_impact_override,
            codon_changes=codon_changes or [],
        )

    # -----------------------------------------------------------------
    # Private helpers
    # -----------------------------------------------------------------

    def _build_model(
        self,
        constraints: list[ConstraintSpec],
        sequence: str,
    ) -> CSPModel:
        """Build a minimal CSPModel from constraints and sequence.

        Needed because :class:`ConflictResolver.detect_conflicts`
        requires a :class:`CSPModel`.
        """
        # Infer protein length from sequence
        protein_len = len(sequence) // 3
        # Build a dummy protein sequence (all A — just for length)
        protein_sequence = "A" * protein_len

        # Build codon_domains — not used by conflict detection but required
        from biocompiler.shared.constants import AA_TO_CODONS
        codon_domains: dict[int, list[str]] = {}
        for i in range(protein_len):
            aa = protein_sequence[i] if i < len(protein_sequence) else "A"
            codon_domains[i] = list(AA_TO_CODONS.get(aa, ["GCT"]))

        return CSPModel(
            protein_sequence=protein_sequence,
            codon_domains=codon_domains,
            constraints=list(constraints),
            config=SolverConfig(),
        )

    def _resolve_one_conflict(
        self,
        conflict: ConstraintConflict,
        constraint_map: dict[str, ConstraintSpec],
        sequence: str,
    ) -> Optional[ConflictProvenance]:
        """Convert a :class:`ConstraintConflict` into a :class:`ConflictProvenance`.

        Determines the winner/loser and resolution method based on the
        conflict's ``resolution_strategy``.
        """
        strategy = conflict.resolution_strategy
        name_a = conflict.constraint_a
        name_b = conflict.constraint_b
        positions = conflict.conflict_positions

        spec_a = constraint_map.get(name_a)
        spec_b = constraint_map.get(name_b)
        prio_a = spec_a.priority if spec_a else ConstraintPriority.MEDIUM
        prio_b = spec_b.priority if spec_b else ConstraintPriority.MEDIUM

        # Determine winner/loser and method from strategy
        if strategy == "relax_a":
            # A is relaxed, B wins
            winner = name_b
            loser = name_a
            method = self._determine_resolution_method(spec_a, spec_b)
            if self._cai_aware:
                method = "cai_aware"
            impact = (
                f"Constraint '{name_a}' was relaxed to satisfy '{name_b}' "
                f"at position(s) {positions}. "
                f"'{name_b}' has higher priority ({prio_b.name}) than "
                f"'{name_a}' ({prio_a.name})."
            )
        elif strategy == "relax_b":
            # B is relaxed, A wins
            winner = name_a
            loser = name_b
            method = self._determine_resolution_method(spec_a, spec_b)
            if self._cai_aware:
                method = "cai_aware"
            impact = (
                f"Constraint '{name_b}' was relaxed to satisfy '{name_a}' "
                f"at position(s) {positions}. "
                f"'{name_a}' has higher priority ({prio_a.name}) than "
                f"'{name_b}' ({prio_b.name})."
            )
        elif strategy == "compromise":
            # Both partially relaxed — pick the one with higher priority
            # as "winner" for provenance, but mark method as weight_based
            # (or cai_aware if CAI-aware resolution is active)
            if prio_a < prio_b:
                winner, loser = name_a, name_b
            elif prio_b < prio_a:
                winner, loser = name_b, name_a
            else:
                winner, loser = name_a, name_b  # arbitrary but consistent
            method = "cai_aware" if self._cai_aware else "weight_based"
            impact = (
                f"Both constraints '{name_a}' and '{name_b}' were partially "
                f"relaxed (compromise) at position(s) {positions}. "
                f"Both have equal priority ({prio_a.name})."
            )
        elif strategy == "infeasible":
            # Cannot auto-resolve — record as manual intervention needed
            winner = "<none>"
            loser = "<both>"
            method = "manual"
            impact = (
                f"Conflict between '{name_a}' and '{name_b}' at position(s) "
                f"{positions} is structurally infeasible. Manual intervention "
                f"or redesign is required."
            )
        else:
            logger.warning(
                "Unknown conflict resolution strategy %r — skipping provenance",
                strategy,
            )
            return None

        return ConflictProvenance(
            conflicting_constraints=[name_a, name_b],
            resolution_method=method,
            winner=winner,
            loser=loser,
            impact=impact,
            positions_affected=positions,
            cai_impact=self._compute_cai_impact(loser, constraint_map, sequence),
            cai_delta=None,
        )

    @staticmethod
    def _determine_resolution_method(
        spec_a: Optional[ConstraintSpec],
        spec_b: Optional[ConstraintSpec],
    ) -> str:
        """Determine the resolution method label for a conflict.

        If the weight difference is significant (>2x), label as
        ``weight_based``; otherwise ``priority_based``.
        """
        weight_a = spec_a.weight if spec_a else 1.0
        weight_b = spec_b.weight if spec_b else 1.0

        if weight_a > 0 and weight_b > 0:
            ratio = max(weight_a, weight_b) / min(weight_a, weight_b)
            if ratio > 2.0:
                return "weight_based"

        return "priority_based"

    def _compute_cai_impact(
        self,
        loser: str,
        constraint_map: dict[str, ConstraintSpec],
        sequence: str,
    ) -> float:
        """Compute the CAI impact for a conflict resolution.

        If CAI-aware resolution is active, uses the CAIAwareConstraintResolver
        for a more precise estimate. Otherwise falls back to the heuristic
        name-based estimation.
        """
        cai_impact = self._estimate_cai_impact(loser, sequence)

        if self._cai_aware:
            from .conflict_resolution import CAIImpactEstimator
            cai_resolver = CAIImpactEstimator()
            loser_spec = constraint_map.get(loser)
            if loser_spec is not None:
                cai_impact = cai_resolver.estimate_cai_impact(loser_spec)

        return cai_impact

    @staticmethod
    def _estimate_cai_impact(
        constraint_name: str,
        sequence: str,
    ) -> float:
        """Estimate the CAI impact of relaxing a constraint.

        Uses a heuristic based on constraint type: constraints that directly
        affect codon choice (GC, codon usage, CpG) tend to have higher
        CAI impact, while positional constraints (restriction sites, splice
        sites) have lower impact since they only affect specific positions.

        Returns a float in [-1.0, +1.0]:
        - Negative: relaxing this constraint *hurts* CAI (less optimal codons
          were chosen to satisfy it, so relaxing would actually help CAI).
        - Positive: relaxing this constraint *helps* CAI.
        - Near zero: minimal CAI impact.

        Parameters
        ----------
        constraint_name:
            Name of the constraint that was relaxed.
        sequence:
            The candidate DNA sequence.

        Returns
        -------
        float
            Estimated CAI impact.
        """
        # Heuristic impact estimates based on constraint name patterns
        name_lower = constraint_name.lower()

        # GC constraints: significant CAI impact because they constrain
        # codon choice globally
        if "gc" in name_lower:
            return -0.05  # Relaxing GC range slightly hurts CAI

        # CpG constraints: moderate CAI impact
        if "cpg" in name_lower or "cpg_island" in name_lower:
            return 0.03  # Relaxing CpG can slightly help CAI

        # Restriction site constraints: low CAI impact (local)
        if "restriction" in name_lower or "eco" in name_lower:
            return 0.01

        # Splice site constraints: low CAI impact (local)
        if "splice" in name_lower or "cryptic" in name_lower:
            return 0.01

        # ATTTA motif constraints: low CAI impact
        if "attta" in name_lower or "instability" in name_lower:
            return 0.005

        # T-run / homopolymer constraints: very low CAI impact
        if "t_run" in name_lower or "homopolymer" in name_lower:
            return 0.002

        # Codon usage / CAI: direct impact
        if "codon_usage" in name_lower or "cai" in name_lower:
            return -0.10

        # Translation constraints: should never be relaxed, but if so,
        # high impact
        if "translation" in name_lower or "amino_acid" in name_lower:
            return -0.50

        # Default: small unknown impact
        return 0.0
