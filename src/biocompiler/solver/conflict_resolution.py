"""BioCompiler CSP Solver — Hard Constraint Conflict Resolution.

When the CSP solver determines that hard constraints cannot all be satisfied
simultaneously, this module provides structured conflict detection and
resolution strategies. It builds on top of the MUS (Minimal Unsatisfiable
Subset) analysis from :mod:`biocompiler.solver.mus` to identify *pairs* of
conflicting constraints and propose targeted relaxations.

Conflict resolution strategies
------------------------------
- **relax_a**: Relax constraint A (lower priority) to satisfy constraint B.
- **relax_b**: Relax constraint B (lower priority) to satisfy constraint A.
- **compromise**: Partially relax both constraints to find a middle ground.
- **infeasible**: The conflict cannot be resolved by relaxation alone; the
  problem is structurally infeasible and requires a design change.

Auto-resolution strategies
--------------------------
- **min_relaxation**: Choose the resolution that relaxes the fewest
  constraints with the smallest total priority delta.
- **max_priority**: Always relax the lowest-priority constraint(s).
- **compromise_first**: Prefer compromise strategies when available.

Public API
----------
- ``ConstraintConflict``    — dataclass representing a pairwise constraint conflict
- ``ConflictResolver``      — main conflict detection and resolution class
- ``prioritize_constraints`` — sort constraints by priority field
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
    MUSReport,
    SolverConfig,
)
from .constraints import (
    HardConstraint,
    SoftConstraint,
    GCRangeConstraint,
    NoCrypticSpliceConstraint,
    NoRestrictionSiteConstraint,
    NoCpGIslandConstraint,
    NoATTTAMotifConstraint,
    NoTRunConstraint,
)

logger = logging.getLogger(__name__)


__all__ = [
    "ConstraintConflict",
    "ConflictResolver",
    "CAIImpactEstimator",
    "prioritize_constraints",
]


# =====================================================================
# ConstraintConflict dataclass
# =====================================================================

@dataclass
class ConstraintConflict:
    """A pairwise conflict between two hard constraints.

    Represents a situation where two constraints cannot both be satisfied
    simultaneously within the current model. The ``resolution_strategy``
    field suggests how to resolve the conflict, and ``conflict_positions``
    identifies the codon positions where the conflict manifests.

    Attributes:
        constraint_a: Name of the first conflicting constraint.
        constraint_b: Name of the second conflicting constraint.
        conflict_positions: Codon positions where both constraints
            are active and in conflict.
        resolution_strategy: Suggested resolution approach. One of
            ``"relax_a"``, ``"relax_b"``, ``"compromise"``, or
            ``"infeasible"``.
    """

    constraint_a: str
    constraint_b: str
    conflict_positions: list[int]
    resolution_strategy: str  # "relax_a" | "relax_b" | "compromise" | "infeasible"

    def __post_init__(self) -> None:
        """Validate resolution_strategy after initialization."""
        valid_strategies = {"relax_a", "relax_b", "compromise", "infeasible"}
        if self.resolution_strategy not in valid_strategies:
            raise ValueError(
                f"Invalid resolution_strategy {self.resolution_strategy!r}; "
                f"must be one of {valid_strategies}"
            )

    def __repr__(self) -> str:
        return (
            f"ConstraintConflict("
            f"{self.constraint_a!r} vs {self.constraint_b!r}, "
            f"positions={self.conflict_positions}, "
            f"strategy={self.resolution_strategy!r})"
        )


# =====================================================================
# prioritize_constraints function
# =====================================================================

def prioritize_constraints(constraints: list[ConstraintSpec]) -> list[ConstraintSpec]:
    """Sort constraints by their ``priority`` field (ascending).

    Lower priority values indicate harder-to-relax (more critical)
    constraints. After sorting, the most critical constraints appear
    first, making it easy to iterate from most to least important.

    Parameters
    ----------
    constraints:
        List of :class:`ConstraintSpec` instances to sort.

    Returns
    -------
    list[ConstraintSpec]
        New list sorted by ``priority`` ascending (1 = most critical first).

    Examples
    --------
    >>> from biocompiler.solver.types import ConstraintSpec, ConstraintType
    >>> c1 = ConstraintSpec(ctype=ConstraintType.GC_CONTENT, name="gc", priority=3)
    >>> c2 = ConstraintSpec(ctype=ConstraintType.NO_CPG, name="cpg", priority=7)
    >>> c3 = ConstraintSpec(ctype=ConstraintType.RESTRICTION_SITE, name="eco", priority=1)
    >>> result = prioritize_constraints([c1, c2, c3])
    >>> [c.name for c in result]
    ['eco', 'gc', 'cpg']
    """
    return sorted(constraints, key=lambda c: c.priority)




# =====================================================================
# CAIImpactEstimator class
# =====================================================================

class CAIImpactEstimator:
    """Rank conflict resolution options by estimated CAI impact."""

    _CAI_IMPACT_TABLE: dict[ConstraintType, float] = {
        ConstraintType.GC_CONTENT: -0.05,
        ConstraintType.CODON_USAGE: -0.10,
        ConstraintType.NO_CPG: 0.03,
        ConstraintType.NO_GT_DINUCLEOTIDE: 0.04,
        ConstraintType.NO_CRYPTIC_SPLICE: 0.03,
        ConstraintType.SPLICE_DONOR_AVOIDANCE: 0.03,
        ConstraintType.RESTRICTION_SITE: 0.01,
        ConstraintType.MRNA_STABILITY: 0.002,
        ConstraintType.NO_INSTABILITY_MOTIF: 0.005,
        ConstraintType.AMINO_ACID_IDENTITY: -0.50,
        ConstraintType.MHC_BINDING: 0.0,
        ConstraintType.TCELL_EPITOPE: 0.0,
        ConstraintType.PROTEIN_STABILITY: 0.0,
        ConstraintType.CUSTOM: 0.0,
    }

    def estimate_cai_impact(self, spec: ConstraintSpec) -> float:
        """Estimate the CAI impact of relaxing a constraint."""
        return self._CAI_IMPACT_TABLE.get(spec.ctype, 0.0)

    def rank_resolution(self, ca: ConstraintSpec, cb: ConstraintSpec, priority_strategy: str) -> str:
        """Adjust the resolution strategy based on CAI impact."""
        cai_a = self.estimate_cai_impact(ca)
        cai_b = self.estimate_cai_impact(cb)
        cai_diff = abs(cai_a - cai_b)
        if cai_diff < 0.03:
            return priority_strategy
        if ca.priority != cb.priority:
            return priority_strategy
        if priority_strategy == "compromise":
            if cai_a > cai_b:
                return "relax_a"
            elif cai_b > cai_a:
                return "relax_b"
        elif priority_strategy == "relax_a":
            if cai_b > cai_a and ca.priority == cb.priority:
                return "relax_b"
        elif priority_strategy == "relax_b":
            if cai_a > cai_b and ca.priority == cb.priority:
                return "relax_a"
        return priority_strategy

# =====================================================================
# ConflictResolver class
# =====================================================================

class ConflictResolver:
    """Detect and resolve hard constraint conflicts in a CSP model.

    This class analyzes a :class:`CSPModel` for pairwise conflicts between
    hard constraints, suggests resolution strategies, and can auto-resolve
    conflicts by relaxing constraints.

    Usage::

        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        for conflict in conflicts:
            suggestion = resolver.suggest_resolution(conflict)
            print(suggestion)
        resolved_model = resolver.auto_resolve(model)

    The resolver uses constraint priority values and position overlap to
    determine which constraints conflict and how to resolve them.
    """

    def __init__(
        self,
        cai_aware: bool = True,
        resolution_strategy: str = "cai_aware",
    ) -> None:
        """Initialize the ConflictResolver."""
        # resolution_strategy overrides cai_aware only when explicitly
        # set to "priority_only". Otherwise, the cai_aware parameter
        # takes precedence (allows cai_aware=False to work).
        if resolution_strategy == "priority_only":
            cai_aware = False
        # Note: we do NOT force cai_aware=True when resolution_strategy
        # is "cai_aware" — that would override an explicit cai_aware=False.

        self._cai_aware = cai_aware
        self._resolution_strategy = resolution_strategy
        self._conflict_cache: dict[int, list[ConstraintConflict]] = {}
        self._cai_resolver: Optional[CAIImpactEstimator] = None

        if self._cai_aware:
            self._cai_resolver = CAIImpactEstimator()
            logger.info(
                "ConflictResolver initialized with resolution_strategy=%r (cai_aware=True)",
                self._resolution_strategy,
            )
        else:
            logger.info(
                "ConflictResolver initialized with resolution_strategy=%r (cai_aware=False)",
                self._resolution_strategy,
            )

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def detect_conflicts(self, model: CSPModel) -> list[ConstraintConflict]:
        """Detect pairwise conflicts between hard constraints.

        Analyzes the model's constraint list to find pairs of hard constraints
        that have overlapping positions and cannot be satisfied simultaneously.
        Uses position overlap as a proxy for conflict: two constraints that
        both apply to the same codon positions and are of different types
        may conflict because they constrain the same variables differently.

        The resolution strategy is determined by comparing the priorities
        of the two constraints:

        - If constraint A has strictly higher priority (lower number), the
          strategy is ``"relax_b"``.
        - If constraint B has strictly higher priority, the strategy is
          ``"relax_a"``.
        - If priorities are equal, the strategy is ``"compromise"``.
        - If both constraints are translation constraints (fundamentally
          non-negotiable), the strategy is ``"infeasible"``.

        Parameters
        ----------
        model:
            The CSP model to analyze for conflicts.

        Returns
        -------
        list[ConstraintConflict]
            List of detected conflicts with suggested resolution strategies.
        """
        logger.info(
            "Detecting conflicts in model with %d constraints",
            len(model.constraints),
        )

        hard_constraints = [
            c for c in model.constraints
            if c.strictness == ConstraintStrictness.HARD
        ]

        conflicts: list[ConstraintConflict] = []
        seen_pairs: set[tuple[str, str]] = set()

        for i, ca in enumerate(hard_constraints):
            for cb in hard_constraints[i + 1:]:
                # Create a canonical pair key to avoid duplicates
                pair_key = tuple(sorted([ca.name, cb.name]))
                if pair_key in seen_pairs:
                    continue

                # Check position overlap
                overlap = self._compute_position_overlap(ca, cb, model)
                if not overlap:
                    continue

                seen_pairs.add(pair_key)

                # Determine resolution strategy
                strategy = self._determine_strategy(ca, cb)

                # If CAI-aware, adjust strategy based on CAI impact
                if self._cai_aware and self._cai_resolver is not None:
                    strategy = self._cai_resolver.rank_resolution(ca, cb, strategy)

                conflict = ConstraintConflict(
                    constraint_a=ca.name,
                    constraint_b=cb.name,
                    conflict_positions=overlap,
                    resolution_strategy=strategy,
                )
                conflicts.append(conflict)
                logger.debug(
                    "Detected conflict: %s vs %s at positions %s (strategy=%s, cai_aware=%s)",
                    ca.name, cb.name, overlap, strategy, self._cai_aware,
                )

        logger.info("Detected %d pairwise conflicts", len(conflicts))

        # Cache result keyed by model id
        self._conflict_cache[id(model)] = conflicts
        return conflicts

    def suggest_resolution(self, conflict: ConstraintConflict) -> str:
        """Generate a human-readable resolution suggestion for a conflict.

        Provides a concrete, actionable suggestion based on the conflict's
        resolution strategy.

        Parameters
        ----------
        conflict:
            The conflict to generate a suggestion for.

        Returns
        -------
        str
            A human-readable suggestion string.
        """
        a = conflict.constraint_a
        b = conflict.constraint_b
        positions = conflict.conflict_positions
        pos_str = f" at position(s) {positions}" if positions else ""

        strategy = conflict.resolution_strategy

        if strategy == "relax_a":
            return (
                f"Relax constraint '{a}'{pos_str} to satisfy '{b}'. "
                f"Constraint '{a}' has equal or lower priority than '{b}', "
                f"so relaxing it has less impact on the design."
            )
        elif strategy == "relax_b":
            return (
                f"Relax constraint '{b}'{pos_str} to satisfy '{a}'. "
                f"Constraint '{b}' has equal or lower priority than '{a}', "
                f"so relaxing it has less impact on the design."
            )
        elif strategy == "compromise":
            return (
                f"Partially relax both '{a}' and '{b}'{pos_str}. "
                f"Both constraints have equal priority — a compromise "
                f"relaxation minimizes total violation across both."
            )
        elif strategy == "infeasible":
            return (
                f"Conflict between '{a}' and '{b}'{pos_str} is structurally "
                f"infeasible. No amount of relaxation can resolve this; "
                f"redesign the target sequence or remove one constraint entirely."
            )
        else:
            return (
                f"Conflict between '{a}' and '{b}'{pos_str} with "
                f"unknown strategy '{strategy}'."
            )

    def auto_resolve(
        self,
        model: CSPModel,
        strategy: str = "min_relaxation",
    ) -> CSPModel:
        """Automatically resolve conflicts by relaxing constraints.

        Applies an auto-resolution strategy to the model's constraint list
        and returns a new :class:`CSPModel` with relaxed constraints.

        Supported strategies:

        - ``"min_relaxation"``: Remove the minimum number of constraints
          (preferring low-priority ones) to resolve all conflicts.
        - ``"max_priority"``: Always remove the lowest-priority constraint
          from each conflict pair.
        - ``"compromise_first"``: For conflicts with ``"compromise"``
          strategy, convert both constraints from HARD to SOFT; for
          others, apply ``"max_priority"``.

        Parameters
        ----------
        model:
            The CSP model to resolve.
        strategy:
            Auto-resolution strategy. One of ``"min_relaxation"``,
            ``"max_priority"``, or ``"compromise_first"``.

        Returns
        -------
        CSPModel
            A new CSP model with conflicts resolved.

        Raises
        ------
        ValueError
            If *strategy* is not a recognized strategy name.
        """
        valid_strategies = {"min_relaxation", "max_priority", "compromise_first"}
        if strategy not in valid_strategies:
            raise ValueError(
                f"Unknown auto-resolution strategy {strategy!r}; "
                f"must be one of {valid_strategies}"
            )

        logger.info(
            "Auto-resolving conflicts with strategy=%s, cai_aware=%s (model has %d constraints)",
            strategy, self._cai_aware, len(model.constraints),
        )

        # Detect conflicts if not cached
        conflicts = self._conflict_cache.get(id(model))
        if conflicts is None:
            conflicts = self.detect_conflicts(model)

        if not conflicts:
            logger.info("No conflicts detected; returning model unchanged")
            return model

        # Build a lookup from constraint name to ConstraintSpec
        constraint_map: dict[str, ConstraintSpec] = {
            c.name: c for c in model.constraints
        }

        # Determine which constraints to relax based on strategy
        names_to_remove: set[str] = set()
        names_to_downgrade: set[str] = set()  # HARD -> SOFT

        if strategy == "min_relaxation":
            if self._cai_aware and self._cai_resolver is not None:
                names_to_remove = self._cai_aware_min_relaxation_select(conflicts, constraint_map)
            else:
                names_to_remove = self._min_relaxation_select(conflicts, constraint_map)
        elif strategy == "max_priority":
            if self._cai_aware and self._cai_resolver is not None:
                names_to_remove = self._cai_aware_max_priority_select(conflicts, constraint_map)
            else:
                names_to_remove = self._max_priority_select(conflicts, constraint_map)
        elif strategy == "compromise_first":
            if self._cai_aware and self._cai_resolver is not None:
                names_to_remove, names_to_downgrade = self._cai_aware_compromise_first_select(conflicts, constraint_map)
            else:
                names_to_remove, names_to_downgrade = self._compromise_first_select(conflicts, constraint_map)

        # Build new constraint list
        new_constraints: list[ConstraintSpec] = []
        for c in model.constraints:
            if c.name in names_to_remove:
                logger.info("  Removing constraint %r (auto-resolved)", c.name)
                continue
            if c.name in names_to_downgrade:
                logger.info("  Downgrading constraint %r from HARD to SOFT", c.name)
                new_c = ConstraintSpec(
                    ctype=c.ctype,
                    name=c.name,
                    strictness=ConstraintStrictness.SOFT,
                    params=dict(c.params),
                    positions=list(c.positions),
                    priority=c.priority,
                    weight=c.weight,
                )
                new_constraints.append(new_c)
            else:
                new_constraints.append(c)

        # Build the resolved model
        resolved_model = CSPModel(
            protein_sequence=model.protein_sequence,
            codon_domains=model.codon_domains,
            constraints=new_constraints,
            config=model.config,
        )

        logger.info(
            "Auto-resolution complete: removed %d constraints, "
            "downgraded %d constraints, %d remaining",
            len(names_to_remove), len(names_to_downgrade),
            len(new_constraints),
        )

        # Invalidate cache for old model since we have resolved it
        self._conflict_cache.pop(id(model), None)

        return resolved_model

    # -----------------------------------------------------------------
    # Private helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _compute_position_overlap(
        ca: ConstraintSpec,
        cb: ConstraintSpec,
        model: CSPModel,
    ) -> list[int]:
        """Compute overlapping positions between two constraints.

        If both constraints have empty ``positions`` lists (meaning they
        apply globally), the overlap is the full sequence length.

        Parameters
        ----------
        ca:
            First constraint.
        cb:
            Second constraint.
        model:
            The CSP model (for sequence length context).

        Returns
        -------
        list[int]
            Sorted list of overlapping codon positions.
        """
        pos_a = set(ca.positions) if ca.positions else set(range(model.length))
        pos_b = set(cb.positions) if cb.positions else set(range(model.length))
        overlap = sorted(pos_a & pos_b)
        return overlap

    @staticmethod
    def _determine_strategy(ca: ConstraintSpec, cb: ConstraintSpec) -> str:
        """Determine the resolution strategy for a pair of conflicting constraints.

        Compares the priority values of the two constraints to decide
        which one to relax. Lower priority values are harder to relax
        (more critical).

        Parameters
        ----------
        ca:
            First constraint (constraint_a in the conflict).
        cb:
            Second constraint (constraint_b in the conflict).

        Returns
        -------
        str
            One of ``"relax_a"``, ``"relax_b"``, ``"compromise"``, or
            ``"infeasible"``.
        """
        # Translation constraints are fundamentally non-negotiable
        both_translation = (
            ca.ctype == ConstraintType.AMINO_ACID_IDENTITY
            and cb.ctype == ConstraintType.AMINO_ACID_IDENTITY
        )
        if both_translation:
            return "infeasible"

        # Compare priorities — lower number = higher priority = harder to relax
        if ca.priority < cb.priority:
            # A is more critical → relax B
            return "relax_b"
        elif cb.priority < ca.priority:
            # B is more critical → relax A
            return "relax_a"
        else:
            # Equal priority → compromise
            return "compromise"

    @staticmethod
    def _min_relaxation_select(
        conflicts: list[ConstraintConflict],
        constraint_map: dict[str, ConstraintSpec],
    ) -> set[str]:
        """Select the minimum set of constraints to remove.

        Uses a greedy set-cover approach: iteratively pick the constraint
        whose removal resolves the most conflicts, breaking ties by
        priority (lower priority = easier to relax).

        Parameters
        ----------
        conflicts:
            Detected conflicts.
        constraint_map:
            Mapping from constraint name to spec.

        Returns
        -------
        set[str]
            Names of constraints to remove.
        """
        # Build a mapping: constraint name -> set of conflict indices it appears in
        remaining: list[tuple[int, ConstraintConflict]] = list(enumerate(conflicts))
        to_remove: set[str] = set()

        while remaining:
            # Count how many conflicts each constraint participates in
            involvement: dict[str, int] = {}
            for _, conflict in remaining:
                involvement[conflict.constraint_a] = (
                    involvement.get(conflict.constraint_a, 0) + 1
                )
                involvement[conflict.constraint_b] = (
                    involvement.get(conflict.constraint_b, 0) + 1
                )

            # Pick the constraint with highest involvement (most conflicts resolved)
            # Break ties by priority rank (higher rank = easier to relax = prefer)
            best_name: Optional[str] = None
            best_score: tuple[int, int] = (-1, -1)  # (involvement, priority_rank)
            for name, count in involvement.items():
                spec = constraint_map.get(name)
                priority = spec.priority if spec else ConstraintPriority.MEDIUM
                score = (count, priority.rank)
                if score > best_score:
                    best_score = score
                    best_name = name

            if best_name is None:
                break

            to_remove.add(best_name)

            # Remove all conflicts involving the chosen constraint
            remaining = [
                (idx, c) for idx, c in remaining
                if c.constraint_a != best_name and c.constraint_b != best_name
            ]

        return to_remove

    @staticmethod
    def _max_priority_select(
        conflicts: list[ConstraintConflict],
        constraint_map: dict[str, ConstraintSpec],
    ) -> set[str]:
        """Select constraints to remove using max-priority (relax lowest priority).

        For each conflict, remove the constraint with the highest priority
        number (easiest to relax).

        Parameters
        ----------
        conflicts:
            Detected conflicts.
        constraint_map:
            Mapping from constraint name to spec.

        Returns
        -------
        set[str]
            Names of constraints to remove.
        """
        to_remove: set[str] = set()
        for conflict in conflicts:
            spec_a = constraint_map.get(conflict.constraint_a)
            spec_b = constraint_map.get(conflict.constraint_b)
            prio_a = spec_a.priority if spec_a else ConstraintPriority.MEDIUM
            prio_b = spec_b.priority if spec_b else ConstraintPriority.MEDIUM
            # Remove the one with higher priority rank (easier to relax)
            if prio_a.rank >= prio_b.rank:
                to_remove.add(conflict.constraint_a)
            else:
                to_remove.add(conflict.constraint_b)
        return to_remove

    @staticmethod
    def _compromise_first_select(
        conflicts: list[ConstraintConflict],
        constraint_map: dict[str, ConstraintSpec],
    ) -> tuple[set[str], set[str]]:
        """Select constraints for compromise-first strategy.

        For conflicts with ``"compromise"`` strategy, downgrade both
        constraints from HARD to SOFT. For other conflicts, apply
        max-priority removal.

        Parameters
        ----------
        conflicts:
            Detected conflicts.
        constraint_map:
            Mapping from constraint name to spec.

        Returns
        -------
        tuple[set[str], set[str]]
            (names_to_remove, names_to_downgrade)
        """
        to_remove: set[str] = set()
        to_downgrade: set[str] = set()

        for conflict in conflicts:
            if conflict.resolution_strategy == "compromise":
                # Downgrade both constraints from HARD to SOFT
                to_downgrade.add(conflict.constraint_a)
                to_downgrade.add(conflict.constraint_b)
            elif conflict.resolution_strategy == "relax_a":
                to_remove.add(conflict.constraint_a)
            elif conflict.resolution_strategy == "relax_b":
                to_remove.add(conflict.constraint_b)
            elif conflict.resolution_strategy == "infeasible":
                # Cannot auto-resolve infeasible conflicts — skip with warning
                logger.warning(
                    "Skipping infeasible conflict: %s vs %s — "
                    "cannot auto-resolve; redesign required",
                    conflict.constraint_a, conflict.constraint_b,
                )
            else:
                logger.warning(
                    "Unknown strategy %r for conflict %s vs %s — skipping",
                    conflict.resolution_strategy,
                    conflict.constraint_a,
                    conflict.constraint_b,
                )

        return to_remove, to_downgrade

    # -----------------------------------------------------------------
    # CAI-aware private helpers
    # -----------------------------------------------------------------

    def _cai_aware_min_relaxation_select(self, conflicts: list[ConstraintConflict], constraint_map: dict[str, ConstraintSpec]) -> set[str]:
        """Select the minimum set of constraints to remove, CAI-aware."""
        remaining: list[tuple[int, ConstraintConflict]] = list(enumerate(conflicts))
        to_remove: set[str] = set()
        while remaining:
            involvement: dict[str, int] = {}
            for _, conflict in remaining:
                involvement[conflict.constraint_a] = involvement.get(conflict.constraint_a, 0) + 1
                involvement[conflict.constraint_b] = involvement.get(conflict.constraint_b, 0) + 1
            best_name: Optional[str] = None
            best_score: tuple[int, int, float] = (-1, -1, -float("inf"))
            for name, count in involvement.items():
                spec = constraint_map.get(name)
                priority = spec.priority if spec else ConstraintPriority.MEDIUM
                cai_impact = self._cai_resolver.estimate_cai_impact(spec) if self._cai_resolver and spec else 0.0
                score = (count, priority.rank, cai_impact)
                if score > best_score:
                    best_score = score
                    best_name = name
            if best_name is None:
                break
            to_remove.add(best_name)
            remaining = [(idx, c) for idx, c in remaining if c.constraint_a != best_name and c.constraint_b != best_name]
        return to_remove

    def _cai_aware_max_priority_select(self, conflicts: list[ConstraintConflict], constraint_map: dict[str, ConstraintSpec]) -> set[str]:
        """Select constraints to remove, CAI-aware max-priority."""
        to_remove: set[str] = set()
        for conflict in conflicts:
            spec_a = constraint_map.get(conflict.constraint_a)
            spec_b = constraint_map.get(conflict.constraint_b)
            prio_a = spec_a.priority if spec_a else ConstraintPriority.MEDIUM
            prio_b = spec_b.priority if spec_b else ConstraintPriority.MEDIUM
            if prio_a.rank > prio_b.rank:
                to_remove.add(conflict.constraint_a)
            elif prio_b.rank > prio_a.rank:
                to_remove.add(conflict.constraint_b)
            else:
                cai_a = self._cai_resolver.estimate_cai_impact(spec_a) if self._cai_resolver and spec_a else 0.0
                cai_b = self._cai_resolver.estimate_cai_impact(spec_b) if self._cai_resolver and spec_b else 0.0
                if cai_a >= cai_b:
                    to_remove.add(conflict.constraint_a)
                else:
                    to_remove.add(conflict.constraint_b)
        return to_remove

    def _cai_aware_compromise_first_select(self, conflicts: list[ConstraintConflict], constraint_map: dict[str, ConstraintSpec]) -> tuple[set[str], set[str]]:
        """Select constraints for compromise-first strategy, CAI-aware."""
        to_remove: set[str] = set()
        to_downgrade: set[str] = set()
        for conflict in conflicts:
            if conflict.resolution_strategy == "compromise":
                to_downgrade.add(conflict.constraint_a)
                to_downgrade.add(conflict.constraint_b)
            elif conflict.resolution_strategy == "relax_a":
                to_remove.add(conflict.constraint_a)
            elif conflict.resolution_strategy == "relax_b":
                to_remove.add(conflict.constraint_b)
            elif conflict.resolution_strategy == "infeasible":
                logger.warning("Skipping infeasible conflict: %s vs %s", conflict.constraint_a, conflict.constraint_b)
            else:
                logger.warning("Unknown strategy %r for conflict %s vs %s", conflict.resolution_strategy, conflict.constraint_a, conflict.constraint_b)
        return to_remove, to_downgrade
