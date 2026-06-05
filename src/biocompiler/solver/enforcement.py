"""Constraint enforcement mechanics for the BioCompiler CSP solver.

This module implements the missing constraint enforcement layer that bridges
the gap between the declarative constraint model (``constraints.py``) and
the solver backends.  It provides:

- **ConstraintEnforcer**: Runs all hard constraints against a candidate
  sequence, collects violations with severity, computes weighted soft
  scores, and identifies conflicting hard constraints with suggested
  relaxations.
- **EnforcementResult**: Detailed result of enforcement including
  violation list, soft score, and optional conflict resolution.
- **ConflictResolution**: Diagnosis of conflicting hard constraints with
  human-readable relaxation suggestions.

Typical usage::

    from biocompiler.solver.enforcement import ConstraintEnforcer

    enforcer = ConstraintEnforcer()
    result = enforcer.enforce(model, candidate_sequence)
    if not result.all_hard_satisfied:
        print("Conflicts:", result.conflict_resolution.conflicting_constraints)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from .types import ConstraintStrictness, ConstraintViolation
from .constraints import CSPModel, HardConstraint, SoftConstraint, MaximizeCAI
from .scoring import SoftConstraintScorer
from ..constants import AA_TO_CODONS, CODON_TABLE

logger = logging.getLogger(__name__)


__all__ = [
    "ConflictResolution",
    "EnforcementResult",
    "CAIAwareFix",
    "CAIAwareEnforcementResult",
    "ConstraintEnforcer",
]


# ==============================================================================
# Data structures
# ==============================================================================


@dataclass
class ConflictResolution:
    """Diagnosis of conflicting hard constraints and suggested relaxations.

    When hard constraints cannot all be satisfied simultaneously, this
    dataclass identifies the subset of constraints that are in conflict
    and proposes relaxations the user can apply to restore feasibility.

    Attributes:
        conflicting_constraints: Names of hard constraints that conflict
            with each other (i.e. cannot all be satisfied at once for the
            given sequence).
        suggested_relaxations: Ordered list of human-readable suggestions
            for resolving the conflict, from least impactful to most
            impactful.  Each suggestion references a specific constraint
            and a proposed parameter change.
        is_resolvable: Whether the conflict can be resolved by relaxing
            one or more constraints (as opposed to fundamental
            incompatibility such as an empty codon domain).
    """

    conflicting_constraints: list[str]
    suggested_relaxations: list[str]
    is_resolvable: bool


@dataclass
class CAIAwareFix:
    """A single CAI-aware constraint fix applied to a violated position.

    Records which codon was changed, what it was changed to, the CAI
    delta of the change, and which constraint violation it resolved.

    Attributes:
        codon_position: 0-based codon position that was fixed.
        old_codon: The codon that was replaced.
        new_codon: The replacement codon (chosen for minimal CAI loss).
        cai_delta: The CAI impact of the change (negative = CAI lost).
        constraint_name: Name of the constraint whose violation was fixed.
    """

    codon_position: int
    old_codon: str
    new_codon: str
    cai_delta: float
    constraint_name: str

    def __repr__(self) -> str:
        return (
            f"CAIAwareFix(pos={self.codon_position}, "
            f"{self.old_codon}->{self.new_codon}, "
            f"cai_delta={self.cai_delta:+.6f}, "
            f"constraint={self.constraint_name!r})"
        )


@dataclass
class CAIAwareEnforcementResult:
    """Result of CAI-aware constraint enforcement.

    Contains the fixed sequence, the list of CAI-aware fixes applied,
    remaining violations (if any), and aggregate CAI impact.

    Attributes:
        sequence: The input sequence with fixes applied.
        fixes: List of :class:`CAIAwareFix` records, one per applied fix.
        total_cai_delta: Sum of all CAI deltas from applied fixes.
            Negative means net CAI loss; positive means net gain.
        remaining_violations: Violations that could not be fixed
            (no valid alternative codon found).
        all_hard_satisfied: Whether all hard constraints are now
            satisfied after applying fixes.
    """

    sequence: str
    fixes: list[CAIAwareFix]
    total_cai_delta: float
    remaining_violations: list[ConstraintViolation]
    all_hard_satisfied: bool

    def __repr__(self) -> str:
        return (
            f"CAIAwareEnforcementResult(fixes={len(self.fixes)}, "
            f"total_cai_delta={self.total_cai_delta:+.6f}, "
            f"remaining={len(self.remaining_violations)}, "
            f"hard_satisfied={self.all_hard_satisfied})"
        )


@dataclass
class EnforcementResult:
    """Detailed result of constraint enforcement on a candidate sequence.

    Produced by :meth:`ConstraintEnforcer.enforce`.  Contains the full
    violation breakdown, the aggregate soft-constraint score, and an
    optional conflict resolution if hard constraints are violated.

    Attributes:
        all_hard_satisfied: True if every hard constraint passes.
        violations: List of all constraint violations (hard and soft).
        soft_score: Weighted aggregate score of all soft constraints
            (higher is better).  Zero if the model has no soft constraints.
        conflict_resolution: If hard constraints are violated, a
            :class:`ConflictResolution` describing the conflict and
            suggested relaxations; ``None`` when all hard constraints are
            satisfied.
    """

    all_hard_satisfied: bool
    violations: list[ConstraintViolation]
    soft_score: float
    conflict_resolution: Optional[ConflictResolution] = None


# ==============================================================================
# ConstraintEnforcer
# ==============================================================================


class ConstraintEnforcer:
    """Runs constraint enforcement against candidate DNA sequences.

    This class provides three complementary operations:

    1. :meth:`enforce` — Run all hard constraints, collect violations,
       compute soft scores, and optionally diagnose conflicts.
    2. :meth:`score_soft_constraints` — Compute the weighted soft
       constraint objective value for a sequence.
    3. :meth:`resolve_conflicts` — Identify which hard constraints
       conflict and suggest parameter relaxations.

    The enforcer is stateless; a single instance can be reused across
    multiple models and sequences.

    Example::

        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGGC...")
        if result.all_hard_satisfied:
            print(f"Soft score: {result.soft_score:.4f}")
        else:
            for v in result.violations:
                print(f"  {v.constraint_name}: {v.description}")
    """

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def enforce(self, model: CSPModel, sequence: str) -> EnforcementResult:
        """Run all hard constraints and return detailed enforcement results.

        Evaluates every hard constraint in *model* against *sequence*.
        For each violated constraint, a :class:`ConstraintViolation` is
        created with the violated positions and a severity estimate.

        If any hard constraint is violated, :meth:`resolve_conflicts` is
        called automatically to diagnose the conflict and suggest
        relaxations.

        Args:
            model: The CSP model containing hard/soft constraints and
                configuration (including objective weights).
            sequence: Candidate DNA sequence to validate.

        Returns:
            An :class:`EnforcementResult` with the full violation list,
            soft score, and optional conflict resolution.
        """
        logger.debug(
            "Enforcing %d hard constraints on sequence of length %d",
            len(model.hard_constraints),
            len(sequence),
        )

        violations: list[ConstraintViolation] = []

        # --- Hard constraints ---
        for hc in model.hard_constraints:
            if hc.check(sequence):
                continue

            positions = hc.violated_positions(sequence)
            severity = self._estimate_severity(positions, len(sequence))
            description = (
                f"Hard constraint '{hc.name}' violated at "
                f"{len(positions)} position(s)"
            )

            violation = ConstraintViolation(
                constraint_name=hc.name,
                constraint_type=ConstraintStrictness.HARD,
                description=description,
                positions=positions,
                severity=severity,
            )
            violations.append(violation)
            logger.debug("  VIOLATED: %s (severity %.2f)", hc.name, severity)

        all_hard_satisfied = len(violations) == 0

        # --- Soft constraint score ---
        soft_score = self.score_soft_constraints(model, sequence)

        # --- Conflict resolution (only when hard constraints fail) ---
        conflict_resolution: Optional[ConflictResolution] = None
        if not all_hard_satisfied:
            conflict_resolution = self.resolve_conflicts(
                model.hard_constraints, sequence
            )

        result = EnforcementResult(
            all_hard_satisfied=all_hard_satisfied,
            violations=violations,
            soft_score=soft_score,
            conflict_resolution=conflict_resolution,
        )

        logger.info(
            "Enforcement complete: hard_satisfied=%s, violations=%d, "
            "soft_score=%.4f",
            all_hard_satisfied,
            len(violations),
            soft_score,
        )
        return result

    def enforce_with_cai_awareness(
        self,
        model: CSPModel,
        sequence: str,
    ) -> CAIAwareEnforcementResult:
        """Enforce constraints considering CAI impact of different fixes.

        Like :meth:`enforce`, but when fixing constraint violations the
        method prefers fixes that minimize CAI loss.  For each violated
        hard constraint:

        1. Identify the violated codon positions.
        2. For each position, enumerate alternative codons from the
           model's codon domains.
        3. Score each alternative by its CAI impact using
           :class:`SoftConstraintScorer.score_cai_impact`.
        4. Apply the fix with the lowest CAI cost (smallest negative
           delta, or the most positive delta).

        Args:
            model: The CSP model containing hard/soft constraints
                and codon domains.
            sequence: Candidate DNA sequence to fix.

        Returns:
            A :class:`CAIAwareEnforcementResult` with the fixed sequence,
            applied fixes, and remaining violations.
        """
        logger.debug(
            "CAI-aware enforcement on sequence of length %d with %d "
            "hard constraints",
            len(sequence), len(model.hard_constraints),
        )

        # Build the scorer for CAI impact analysis
        scorer = SoftConstraintScorer(model.config)

        # Work on a mutable copy of the sequence
        seq_list = list(sequence)
        fixes: list[CAIAwareFix] = []
        total_cai_delta = 0.0

        # Identify all violations
        violations: list[ConstraintViolation] = []
        for hc in model.hard_constraints:
            if hc.check(sequence):
                continue
            positions = hc.violated_positions(sequence)
            severity = self._estimate_severity(positions, len(sequence))
            violation = ConstraintViolation(
                constraint_name=hc.name,
                constraint_type=ConstraintStrictness.HARD,
                description=f"Hard constraint '{hc.name}' violated",
                positions=positions,
                severity=severity,
            )
            violations.append(violation)

        if not violations:
            return CAIAwareEnforcementResult(
                sequence=sequence,
                fixes=[],
                total_cai_delta=0.0,
                remaining_violations=[],
                all_hard_satisfied=True,
            )

        # For each violation, try CAI-aware fixes at violated positions
        remaining_violations: list[ConstraintViolation] = []
        fixed_violations: set[str] = set()

        for violation in violations:
            # Convert nucleotide positions to codon positions
            codon_positions = sorted(
                set(p // 3 for p in violation.positions)
            )

            violation_fixed = False

            for codon_pos in codon_positions:
                if codon_pos < 0 or codon_pos >= len(model.protein):
                    continue

                # Get the current codon
                nt_start = codon_pos * 3
                current_codon = "".join(seq_list[nt_start : nt_start + 3])

                # Get alternative codons from model's codon domains
                alternatives = model.codon_domains.get(codon_pos, [])
                # Filter out the current codon
                alternatives = [c for c in alternatives if c != current_codon]

                if not alternatives:
                    continue

                # Score alternatives by CAI impact
                current_seq = "".join(seq_list)
                cai_deltas = scorer.score_cai_impact(
                    model, current_seq, codon_pos, alternatives,
                )

                # Try each alternative in order of CAI impact (best first)
                sorted_alts = sorted(
                    cai_deltas.items(), key=lambda kv: kv[1], reverse=True,
                )

                for alt_codon, cai_delta in sorted_alts:
                    # Apply the fix tentatively
                    for j, base in enumerate(alt_codon):
                        seq_list[nt_start + j] = base

                    test_seq = "".join(seq_list)

                    # Check if this fix resolves the specific violation
                    hc = None
                    for h in model.hard_constraints:
                        if h.name == violation.constraint_name:
                            hc = h
                            break

                    if hc is not None and hc.check(test_seq):
                        # Fix succeeded — record it
                        fix = CAIAwareFix(
                            codon_position=codon_pos,
                            old_codon=current_codon,
                            new_codon=alt_codon,
                            cai_delta=cai_delta,
                            constraint_name=violation.constraint_name,
                        )
                        fixes.append(fix)
                        total_cai_delta += cai_delta
                        violation_fixed = True
                        break
                    else:
                        # Revert the change
                        for j, base in enumerate(current_codon):
                            seq_list[nt_start + j] = base

                if violation_fixed:
                    fixed_violations.add(violation.constraint_name)
                    break  # Move to next violation

            if not violation_fixed:
                remaining_violations.append(violation)

        # Build final sequence and check all hard constraints
        final_sequence = "".join(seq_list)
        all_hard_satisfied = all(
            hc.check(final_sequence) for hc in model.hard_constraints
        )

        logger.info(
            "CAI-aware enforcement: applied %d fixes, "
            "total_cai_delta=%+.6f, remaining=%d, hard_satisfied=%s",
            len(fixes), total_cai_delta,
            len(remaining_violations), all_hard_satisfied,
        )

        return CAIAwareEnforcementResult(
            sequence=final_sequence,
            fixes=fixes,
            total_cai_delta=total_cai_delta,
            remaining_violations=remaining_violations,
            all_hard_satisfied=all_hard_satisfied,
        )

    def score_soft_constraints(self, model: CSPModel, sequence: str) -> float:
        """Compute the weighted soft constraint score for a sequence.

        Uses the objective weights from ``model.config`` to combine
        individual soft constraint scores into a single scalar.  The
        convention is that higher is better (minimization objectives
        return negated costs).

        If the model has no soft constraints, returns 0.0.

        Args:
            model: The CSP model containing soft constraints and config.
            sequence: Candidate DNA sequence to score.

        Returns:
            Weighted aggregate soft score (higher = better).
        """
        if not model.soft_constraints:
            return 0.0

        # Build weight map from config — mirrors CSPModel.objective_value()
        weight_map: dict[str, float] = {
            "MaximizeCAI": model.config.cai_weight,
            "MinimizeCpG": model.config.cpg_weight,
            "MinimizeMRNADG": model.config.mrna_dg_weight,
        }

        total = 0.0
        for sc in model.soft_constraints:
            weight = weight_map.get(sc.name, 0.0)
            # Default weight for custom/unknown soft constraints: 0.1
            if weight == 0.0 and sc.name not in weight_map:
                weight = 0.1
                logger.debug(
                    "Using default weight 0.1 for unknown soft "
                    "constraint '%s'",
                    sc.name,
                )
            total += weight * sc.score(sequence)

        logger.debug("Soft constraint score: %.4f", total)
        return total

    def resolve_conflicts(
        self, hard_constraints: list[HardConstraint], sequence: str
    ) -> ConflictResolution:
        """Identify conflicting hard constraints and suggest relaxations.

        A set of hard constraints *conflict* when they cannot all be
        satisfied simultaneously for the given sequence.  This method
        detects pairwise and group conflicts by analysing violated
        constraints and their positions, then proposes targeted
        parameter relaxations.

        The relaxation heuristic prefers relaxing constraints that:

        1. Have the widest GC range (easiest to widen further).
        2. Involve restriction sites (can remove individual enzymes).
        3. Involve splice thresholds (can be raised).
        4. Have the most violated positions (suggesting they are hardest
           to satisfy with current parameters).

        Args:
            hard_constraints: List of hard constraint instances.
            sequence: Candidate DNA sequence that violates one or more
                constraints.

        Returns:
            A :class:`ConflictResolution` identifying the conflicting
            constraints and suggesting relaxations.
        """
        # Identify which constraints are violated
        violated: list[HardConstraint] = []
        for hc in hard_constraints:
            if not hc.check(sequence):
                violated.append(hc)

        if not violated:
            # No conflicts — all constraints satisfied
            return ConflictResolution(
                conflicting_constraints=[],
                suggested_relaxations=[],
                is_resolvable=True,
            )

        conflicting_names = [hc.name for hc in violated]

        # --- Detect positional overlap (constraints that fight over the
        #     same nucleotide positions are more likely to conflict) ---
        position_map: dict[str, set[int]] = {}
        for hc in violated:
            positions = hc.violated_positions(sequence)
            position_map[hc.name] = set(positions)

        # Check pairwise positional overlap
        overlapping_pairs: list[tuple[str, str]] = []
        for i, name_i in enumerate(conflicting_names):
            for name_j in conflicting_names[i + 1 :]:
                if position_map[name_i] & position_map[name_j]:
                    overlapping_pairs.append((name_i, name_j))

        # --- Generate relaxation suggestions ---
        suggestions = self._generate_relaxations(violated, overlapping_pairs)

        # A conflict is considered resolvable if at least one relaxation
        # suggestion exists and the conflict is not due to an empty model.
        is_resolvable = len(suggestions) > 0

        logger.info(
            "Conflict resolution: %d conflicting constraints, "
            "%d overlapping pairs, resolvable=%s",
            len(conflicting_names),
            len(overlapping_pairs),
            is_resolvable,
        )

        return ConflictResolution(
            conflicting_constraints=conflicting_names,
            suggested_relaxations=suggestions,
            is_resolvable=is_resolvable,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _estimate_severity(positions: list[int], seq_length: int) -> float:
        """Estimate violation severity as the fraction of affected positions.

        Args:
            positions: Nucleotide positions where the violation occurs.
            seq_length: Total length of the sequence.

        Returns:
            Severity in [0.0, 1.0].  Returns 1.0 for an empty sequence
            (degenerate case) and 0.0 if no positions are affected.
        """
        if seq_length == 0:
            return 1.0
        if not positions:
            return 0.0
        # Clamp to [0.0, 1.0]
        severity = min(len(set(positions)) / seq_length, 1.0)
        return severity

    @staticmethod
    def _generate_relaxations(
        violated: list[HardConstraint],
        overlapping_pairs: list[tuple[str, str]],
    ) -> list[str]:
        """Generate human-readable relaxation suggestions for violated constraints.

        Each suggestion references a specific constraint and proposes a
        concrete parameter change.  Suggestions are ordered from least to
        most impactful.

        Args:
            violated: List of violated hard constraints.
            overlapping_pairs: Pairs of constraint names that share
                violated positions (indicating positional conflict).

        Returns:
            Ordered list of relaxation suggestion strings.
        """
        suggestions: list[str] = []

        for hc in violated:
            name = hc.name

            if name == "GCRangeConstraint":
                suggestions.append(
                    "Widen GC range (e.g. increase gc_hi or decrease gc_lo)"
                )
            elif name == "NoRestrictionSiteConstraint":
                suggestions.append(
                    "Remove one or more restriction enzymes from the "
                    "avoidance list"
                )
            elif name == "NoCrypticSpliceConstraint":
                suggestions.append(
                    "Raise cryptic splice threshold (e.g. from 3.0 to 5.0) "
                    "to tolerate weaker splice sites"
                )
            elif name == "NoCpGIslandConstraint":
                suggestions.append(
                    "Increase CpG Obs/Exp threshold or increase the "
                    "sliding window size"
                )
            elif name == "NoATTTAMotifConstraint":
                suggestions.append(
                    "Disable ATTTA motif avoidance (set avoid_attta=False) "
                    "if mRNA stability is not critical"
                )
            elif name == "NoTRunConstraint":
                suggestions.append(
                    "Increase the max T-run length (e.g. from 5 to 8) "
                    "to allow longer poly-T stretches"
                )
            elif name == "TranslationConstraint":
                suggestions.append(
                    "Verify the input protein sequence matches the codon "
                    "domains — translation errors usually indicate a "
                    "model construction bug"
                )
            else:
                suggestions.append(
                    f"Consider relaxing or disabling constraint "
                    f"'{name}'"
                )

        # Add pairwise overlap warnings
        for name_i, name_j in overlapping_pairs:
            suggestions.append(
                f"Constraints '{name_i}' and '{name_j}' overlap at "
                f"shared positions — relaxing one may resolve both"
            )

        return suggestions
