"""BioCompiler CSP Solver — CAI-Aware Constraint Conflict Resolution.

This module extends the conflict resolution system with **CAI-awareness**:
when two hard constraints conflict, the resolver evaluates which resolution
option has the smallest impact on the Codon Adaptation Index (CAI) and
prefers that resolution.  This avoids the situation where a naive
priority-based resolution unnecessarily sacrifices codon optimality.

Key concepts
------------
- **CAI impact estimation**: For each constraint involved in a conflict,
  we estimate how much satisfying (or relaxing) that constraint would
  change the overall CAI score.  The estimate is per-position and
  constraint-type-aware:
    - *Codon-usage constraints*: direct CAI weight impact — the
      difference between the best codon's adaptiveness and the codon
      required by the constraint.
    - *GC constraints*: estimated from the CAI weight of GC-rich vs
      GC-poor codons for the amino acid at each position.
    - *Splice constraints*: estimated from the CAI penalty of excluding
      GT/AG-containing codons from the domain.

- **Min-CAI-loss resolution**: Given two conflicting constraints A and B,
  we compute ``cai_loss(A_relaxed)`` and ``cai_loss(B_relaxed)`` and
  choose the option with the smaller CAI loss.

Public API
----------
- ``ResolutionChoice``          — dataclass representing a conflict resolution decision
- ``CAIAwareConstraintResolver`` — main CAI-aware conflict resolver class
"""

from __future__ import annotations

import copy
import logging
import math
from dataclasses import dataclass, field
from typing import Optional

from biocompiler.solver.types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
)
from .constraints import (
    MaximizeCAI,
    codon_gc_count,
    codon_contains_gt,
    codon_contains_ag,
    CAI_LOG_EPSILON,
)
# NOTE: We intentionally do NOT import from .conflict_resolution at module
# level to avoid a circular import (conflict_resolution imports from this
# module).  Instead, we use lazy imports inside methods that need them.

logger = logging.getLogger(__name__)


__all__ = [
    "ResolutionChoice",
    "CAIAwareConstraintResolver",
]


# =====================================================================
# ResolutionChoice dataclass
# =====================================================================

@dataclass
class ResolutionChoice:
    """The result of a CAI-aware conflict resolution decision.

    Encapsulates which constraint to relax, the estimated CAI impact
    of that decision, and the alternative that was rejected.

    Attributes:
        relax_constraint: Name of the constraint that should be relaxed.
        keep_constraint: Name of the constraint that should be kept.
        cai_loss: Estimated CAI loss from relaxing the chosen constraint.
            Non-negative — the CAI lost by choosing this resolution.
        alternative_relax: Name of the constraint that was *not* relaxed
            (the alternative).
        alternative_cai_loss: Estimated CAI loss if the *other* constraint
            had been relaxed instead.
        conflict_positions: Codon positions where the conflict manifests.
        strategy: The resolution strategy chosen (``"relax_a"``,
            ``"relax_b"``, ``"compromise"``, or ``"cai_aware"``).
    """

    relax_constraint: str
    keep_constraint: str
    cai_loss: float
    alternative_relax: str
    alternative_cai_loss: float
    conflict_positions: list[int]
    strategy: str = "cai_aware"

    def __post_init__(self) -> None:
        """Validate strategy after initialization."""
        valid_strategies = {"relax_a", "relax_b", "compromise", "cai_aware"}
        if self.strategy not in valid_strategies:
            raise ValueError(
                f"Invalid strategy {self.strategy!r}; "
                f"must be one of {valid_strategies}"
            )

    @property
    def cai_savings(self) -> float:
        """How much CAI was saved by choosing this resolution over the alternative.

        Positive means this resolution is better for CAI than the alternative.
        """
        return self.alternative_cai_loss - self.cai_loss

    def __repr__(self) -> str:
        return (
            f"ResolutionChoice("
            f"relax={self.relax_constraint!r}, keep={self.keep_constraint!r}, "
            f"cai_loss={self.cai_loss:.4f}, "
            f"alt_loss={self.alternative_cai_loss:.4f}, "
            f"savings={self.cai_savings:+.4f}, "
            f"strategy={self.strategy!r})"
        )


# =====================================================================
# CAIAwareConstraintResolver
# =====================================================================

class CAIAwareConstraintResolver:
    """Resolve hard constraint conflicts with minimal CAI loss.

    When constraints conflict, this resolver evaluates the CAI impact
    of each possible resolution (relax A or relax B) and chooses the
    one that minimizes CAI loss.  This is in contrast to the base
    :class:`ConflictResolver` which resolves purely by priority.

    The resolver uses the :class:`MaximizeCAI` scoring methods
    (``score()`` and ``cai()``) to compute actual CAI values and
    estimates per-position CAI impact by examining how each constraint
    restricts the codon domain.

    Usage::

        from biocompiler.solver.cai_aware_resolver import CAIAwareConstraintResolver
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        adaptiveness = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        resolver = CAIAwareConstraintResolver(model, adaptiveness)

        # Resolve a specific conflict
        choice = resolver.resolve_with_min_cai_loss(
            conflicting_constraints=[spec_a, spec_b],
            current_sequence="ATGGGCTGA...",
            codon_domains=model.codon_domains,
        )
        print(f"Relax {choice.relax_constraint}, CAI loss = {choice.cai_loss:.4f}")

        # Or auto-resolve all conflicts
        resolved_model = resolver.auto_resolve_cai_aware(model)

    Parameters
    ----------
    model : CSPModel
        The CSP model containing constraints and codon domains.
    adaptiveness : dict[str, float]
        Codon adaptiveness table for the target organism (codon → 0.0–1.0).
        Used for CAI impact estimation.
    """

    def __init__(
        self,
        model: CSPModel,
        adaptiveness: dict[str, float],
    ) -> None:
        self._model = model
        self._adaptiveness = dict(adaptiveness)
        # Lazy import to avoid circular dependency
        from .conflict_resolution import ConflictResolver
        self._base_resolver = ConflictResolver()

        # Build a MaximizeCAI instance for CAI computation
        self._cai_constraint = MaximizeCAI(
            adaptiveness=self._adaptiveness,
            protein=model.protein_sequence,
        )

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def resolve_with_min_cai_loss(
        self,
        conflicting_constraints: list[ConstraintSpec],
        current_sequence: str,
        codon_domains: dict[int, list[str]],
    ) -> ResolutionChoice:
        """Choose a resolution for conflicting constraints that minimises CAI loss.

        For each of the two possible resolutions (relax A, keep B or
        relax B, keep A), compute the estimated CAI loss.  Return the
        resolution with the smaller loss.

        Parameters
        ----------
        conflicting_constraints:
            A list of exactly two :class:`ConstraintSpec` instances that
            conflict.  Must have length 2.
        current_sequence:
            The current candidate DNA sequence (used for baseline CAI
            computation and position-specific analysis).
        codon_domains:
            Mapping from codon position index to the list of allowed
            codons at that position.

        Returns
        -------
        ResolutionChoice
            The resolution that minimises CAI loss, with impact estimates
            for both alternatives.

        Raises
        ------
        ValueError
            If *conflicting_constraints* does not contain exactly 2
            constraints.
        """
        if len(conflicting_constraints) != 2:
            raise ValueError(
                f"Expected exactly 2 conflicting constraints, "
                f"got {len(conflicting_constraints)}"
            )

        spec_a = conflicting_constraints[0]
        spec_b = conflicting_constraints[1]

        # Compute CAI impact of keeping each constraint (i.e. the CAI loss
        # if the *other* constraint is relaxed and this one is enforced).
        # Relaxing A means B is kept → CAI loss from removing A's influence
        # Relaxing B means A is kept → CAI loss from removing B's influence

        # Find overlap positions
        positions_a = set(spec_a.positions) if spec_a.positions else set(range(self._model.length))
        positions_b = set(spec_b.positions) if spec_b.positions else set(range(self._model.length))
        overlap = sorted(positions_a & positions_b)

        # Compute CAI impact of each constraint at overlapping positions
        cai_impact_a = sum(
            self.estimate_cai_impact(spec_a, pos, codon_domains, self._adaptiveness)
            for pos in overlap
        )
        cai_impact_b = sum(
            self.estimate_cai_impact(spec_b, pos, codon_domains, self._adaptiveness)
            for pos in overlap
        )

        # If we relax A (keep B), the CAI loss is the impact of A that
        # we lose minus the impact of B that we gain.
        # Think of it this way:
        #   - A was constraining codon choice → relaxing A frees up codons
        #     which may improve CAI → negative CAI loss (CAI gain)
        #   - B is still constraining codon choice → B's CAI impact persists
        #
        # Net CAI loss when relaxing A = impact of keeping B
        # Net CAI loss when relaxing B = impact of keeping A
        #
        # We want to relax the constraint whose *keeping* causes more CAI loss.
        # So we relax the one with the higher "keep impact" = higher cai_impact.

        # cai_loss_relax_a = CAI loss incurred by keeping B enforced
        # cai_loss_relax_b = CAI loss incurred by keeping A enforced
        cai_loss_relax_a = cai_impact_b  # Keeping B → B's CAI impact persists
        cai_loss_relax_b = cai_impact_a  # Keeping A → A's CAI impact persists

        # Choose the resolution with lower CAI loss
        if cai_loss_relax_a <= cai_loss_relax_b:
            # Relaxing A is better for CAI (or equal)
            return ResolutionChoice(
                relax_constraint=spec_a.name,
                keep_constraint=spec_b.name,
                cai_loss=cai_loss_relax_a,
                alternative_relax=spec_b.name,
                alternative_cai_loss=cai_loss_relax_b,
                conflict_positions=overlap,
                strategy="cai_aware",
            )
        else:
            # Relaxing B is better for CAI
            return ResolutionChoice(
                relax_constraint=spec_b.name,
                keep_constraint=spec_a.name,
                cai_loss=cai_loss_relax_b,
                alternative_relax=spec_a.name,
                alternative_cai_loss=cai_loss_relax_a,
                conflict_positions=overlap,
                strategy="cai_aware",
            )

    def estimate_cai_impact(
        self,
        constraint: ConstraintSpec,
        position: int,
        codon_domains: dict[int, list[str]],
        adaptiveness: dict[str, float],
    ) -> float:
        """Estimate the CAI impact of satisfying a constraint at a given position.

        Returns a non-negative float representing how much CAI is lost by
        this constraint at this position.  Higher values mean the constraint
        forces a bigger deviation from the optimal codon.

        For different constraint types:

        - **Codon-usage constraints** (CODON_USAGE): Direct CAI weight
          impact — the difference between the best codon's adaptiveness
          and the average adaptiveness of remaining domain codons.
        - **GC constraints** (GC_CONTENT): Estimated from the GC content
          of codons — high-GC codons may have different CAI weights than
          low-GC codons for the same amino acid.
        - **Splice constraints** (NO_CRYPTIC_SPLICE,
          SPLICE_DONOR_AVOIDANCE): Estimated from the CAI penalty of
          excluding GT/AG-containing codons from the domain.
        - **Other constraints**: A default small impact is returned based
          on the constraint type's typical effect on codon choice.

        Parameters
        ----------
        constraint:
            The constraint specification to evaluate.
        position:
            The codon position (0-based index into the protein) where
            the constraint applies.
        codon_domains:
            Mapping from position index to list of allowed codons.
        adaptiveness:
            Codon adaptiveness table (codon → 0.0–1.0).

        Returns
        -------
        float
            Estimated CAI impact at this position (non-negative).
            Zero means the constraint has no estimated CAI impact.
        """
        domain = codon_domains.get(position, [])
        if not domain:
            return 0.0

        ctype = constraint.ctype

        # ── Codon-usage constraints: direct CAI weight impact ────────
        if ctype == ConstraintType.CODON_USAGE:
            return self._cai_impact_codon_usage(domain, adaptiveness)

        # ── GC constraints: estimate from codon GC vs CAI ───────────
        if ctype == ConstraintType.GC_CONTENT:
            return self._cai_impact_gc(domain, adaptiveness, constraint)

        # ── Splice constraints: estimate from removing GT/AG codons ──
        if ctype in (ConstraintType.NO_CRYPTIC_SPLICE, ConstraintType.SPLICE_DONOR_AVOIDANCE):
            return self._cai_impact_splice(domain, adaptiveness)

        # ── CpG constraints: estimate from removing CG-containing codons ─
        if ctype == ConstraintType.NO_CPG:
            return self._cai_impact_cpg(domain, adaptiveness)

        # ── Restriction site constraints: typically local, low impact ─
        if ctype == ConstraintType.RESTRICTION_SITE:
            return self._cai_impact_restriction(domain, adaptiveness)

        # ── Instability motif (ATTTA): low CAI impact ───────────────
        if ctype == ConstraintType.NO_INSTABILITY_MOTIF:
            return self._cai_impact_motif(domain, adaptiveness)

        # ── T-run / mRNA stability: very low CAI impact ─────────────
        if ctype == ConstraintType.MRNA_STABILITY:
            return self._cai_impact_trun(domain, adaptiveness)

        # ── GT dinucleotide avoidance ────────────────────────────────
        if ctype == ConstraintType.NO_GT_DINUCLEOTIDE:
            return self._cai_impact_splice(domain, adaptiveness)

        # ── Amino acid identity: should never be relaxed ────────────
        if ctype == ConstraintType.AMINO_ACID_IDENTITY:
            return 1.0  # Maximum impact — translation must be preserved

        # ── Unknown / custom: small default impact ───────────────────
        return 0.0

    def auto_resolve_cai_aware(
        self,
        model: CSPModel,
        current_sequence: str = "",
    ) -> CSPModel:
        """Auto-resolve all conflicts using CAI-aware strategy.

        Detects all pairwise conflicts in the model and resolves each
        one by choosing the resolution that minimises CAI loss.

        Parameters
        ----------
        model:
            The CSP model to resolve.
        current_sequence:
            The current candidate DNA sequence.  If empty, a default
            sequence is built from the model's codon domains using
            highest-CAI codons.

        Returns
        -------
        CSPModel
            A new CSP model with conflicts resolved (some constraints
            removed or downgraded).
        """
        logger.info(
            "CAI-aware auto-resolve: model has %d constraints",
            len(model.constraints),
        )

        # Detect conflicts
        conflicts = self._base_resolver.detect_conflicts(model)

        if not conflicts:
            logger.info("No conflicts detected; returning model unchanged")
            return model

        # Build sequence if not provided
        if not current_sequence:
            current_sequence = self._build_default_sequence(model)

        # Build constraint name → spec lookup
        constraint_map: dict[str, ConstraintSpec] = {
            c.name: c for c in model.constraints
        }

        # Resolve each conflict and collect constraints to relax
        names_to_remove: set[str] = set()
        names_to_downgrade: set[str] = set()

        for conflict in conflicts:
            spec_a = constraint_map.get(conflict.constraint_a)
            spec_b = constraint_map.get(conflict.constraint_b)

            if spec_a is None or spec_b is None:
                logger.warning(
                    "Conflict references unknown constraint: %s vs %s — skipping",
                    conflict.constraint_a, conflict.constraint_b,
                )
                continue

            # Infeasible conflicts cannot be auto-resolved
            if conflict.resolution_strategy == "infeasible":
                logger.warning(
                    "Skipping infeasible conflict: %s vs %s",
                    conflict.constraint_a, conflict.constraint_b,
                )
                continue

            choice = self.resolve_with_min_cai_loss(
                conflicting_constraints=[spec_a, spec_b],
                current_sequence=current_sequence,
                codon_domains=model.codon_domains,
            )

            logger.info(
                "CAI-aware resolution: relax %r (cai_loss=%.4f) "
                "instead of %r (cai_loss=%.4f), savings=%.4f",
                choice.relax_constraint, choice.cai_loss,
                choice.alternative_relax, choice.alternative_cai_loss,
                choice.cai_savings,
            )

            # For compromise situations, downgrade instead of remove
            if conflict.resolution_strategy == "compromise":
                names_to_downgrade.add(choice.relax_constraint)
            else:
                names_to_remove.add(choice.relax_constraint)

        # Build new constraint list
        new_constraints: list[ConstraintSpec] = []
        for c in model.constraints:
            if c.name in names_to_remove:
                logger.info("  Removing constraint %r (CAI-aware auto-resolved)", c.name)
                continue
            if c.name in names_to_downgrade:
                logger.info("  Downgrading constraint %r from HARD to SOFT (CAI-aware)", c.name)
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

        resolved_model = CSPModel(
            protein_sequence=model.protein_sequence,
            codon_domains=model.codon_domains,
            constraints=new_constraints,
            config=model.config,
        )

        logger.info(
            "CAI-aware auto-resolve complete: removed %d, downgraded %d, "
            "%d remaining",
            len(names_to_remove), len(names_to_downgrade),
            len(new_constraints),
        )

        return resolved_model

    # -----------------------------------------------------------------
    # Private helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _best_cai_in_domain(
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """Return the highest adaptiveness value in the domain."""
        if not domain:
            return 0.0
        return max(adaptiveness.get(c, 0.0) for c in domain)

    @staticmethod
    def _avg_cai_in_domain(
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """Return the average adaptiveness value across the domain."""
        if not domain:
            return 0.0
        values = [adaptiveness.get(c, 0.0) for c in domain]
        return sum(values) / len(values)

    def _cai_impact_codon_usage(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """CAI impact for codon-usage constraints.

        Measures the difference between the best possible CAI and the
        average CAI across the domain.  A larger gap means the constraint
        is more likely to force a sub-optimal codon choice.
        """
        best = self._best_cai_in_domain(domain, adaptiveness)
        avg = self._avg_cai_in_domain(domain, adaptiveness)
        # Impact is how much CAI we lose on average vs. optimal
        return max(0.0, best - avg)

    def _cai_impact_gc(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
        constraint: ConstraintSpec,
    ) -> float:
        """CAI impact for GC content constraints.

        Estimates the CAI penalty of GC constraining codon choice.
        If GC needs to be raised, AT-rich codons (which may be
        sub-optimal) are excluded; if GC needs to be lowered, GC-rich
        codons (which often have higher CAI) are excluded.

        The impact is the CAI difference between GC-rich and GC-poor
        codons in the domain, weighted by the GC constraint tightness.
        """
        gc_lo = constraint.params.get("gc_lo", 0.30)
        gc_hi = constraint.params.get("gc_hi", 0.70)

        # Tightness: how narrow the GC range is
        gc_range = gc_hi - gc_lo
        # Normalise: a range of 0.4 (0.30–0.70) is "normal"; tighter ranges
        # have higher impact.  Clamp to [0, 1].
        tightness = max(0.0, min(1.0, 1.0 - gc_range))

        # Separate domain codons by GC content
        gc_rich = [c for c in domain if codon_gc_count(c) >= 2]
        gc_poor = [c for c in domain if codon_gc_count(c) <= 1]

        if not gc_rich or not gc_poor:
            # All codons have similar GC — minimal impact
            return 0.0

        avg_cai_rich = self._avg_cai_in_domain(gc_rich, adaptiveness)
        avg_cai_poor = self._avg_cai_in_domain(gc_poor, adaptiveness)

        # The CAI impact is the difference between GC-rich and GC-poor
        # codon adaptiveness, modulated by constraint tightness
        impact = abs(avg_cai_rich - avg_cai_poor) * tightness
        return impact

    def _cai_impact_splice(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """CAI impact for splice site avoidance constraints.

        Measures the CAI penalty of excluding GT/AG-containing codons
        from the domain.  These codons might have been optimal for CAI.
        """
        # Find codons containing GT or AG dinucleotides
        gt_ag_codons = [c for c in domain if codon_contains_gt(c) or codon_contains_ag(c)]
        safe_codons = [c for c in domain if not codon_contains_gt(c) and not codon_contains_ag(c)]

        if not gt_ag_codons:
            # No GT/AG codons in domain → constraint has no CAI impact
            return 0.0

        if not safe_codons:
            # All codons contain GT/AG → constraint cannot be satisfied
            # without changing the amino acid (infeasible)
            return 1.0

        # Impact = best CAI among excluded codons - best CAI among safe codons
        best_excluded = self._best_cai_in_domain(gt_ag_codons, adaptiveness)
        best_safe = self._best_cai_in_domain(safe_codons, adaptiveness)

        # If excluded codons have higher CAI, the impact is positive
        return max(0.0, best_excluded - best_safe)

    def _cai_impact_cpg(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """CAI impact for CpG island avoidance constraints.

        Measures the CAI penalty of excluding CpG-containing codons.
        CpG dinucleotides (CG) are found in arginine codons (CGN) and
        some others; these often have high CAI in GC-rich organisms.
        """
        from .constraints import codon_contains_cpg

        cpg_codons = [c for c in domain if codon_contains_cpg(c)]
        non_cpg_codons = [c for c in domain if not codon_contains_cpg(c)]

        if not cpg_codons:
            return 0.0

        if not non_cpg_codons:
            return 0.5  # Moderate impact — all codons have CpG

        best_cpg = self._best_cai_in_domain(cpg_codons, adaptiveness)
        best_non_cpg = self._best_cai_in_domain(non_cpg_codons, adaptiveness)

        return max(0.0, best_cpg - best_non_cpg)

    def _cai_impact_restriction(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """CAI impact for restriction site avoidance constraints.

        Restriction site constraints are typically local (affecting only
        specific codon positions within the site).  The CAI impact is
        estimated as a small fraction of the maximum domain CAI
        difference, since restriction sites span multiple codons and
        often only one codon needs to change.
        """
        best = self._best_cai_in_domain(domain, adaptiveness)
        avg = self._avg_cai_in_domain(domain, adaptiveness)
        # Restriction sites are local — scale down the impact
        return max(0.0, (best - avg) * 0.3)

    def _cai_impact_motif(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """CAI impact for instability motif (ATTTA) avoidance.

        ATTTA motifs span codon boundaries and may force specific
        codon choices at adjacent positions.  The impact is typically
        small since many alternative codons avoid the motif.
        """
        best = self._best_cai_in_domain(domain, adaptiveness)
        avg = self._avg_cai_in_domain(domain, adaptiveness)
        return max(0.0, (best - avg) * 0.2)

    def _cai_impact_trun(
        self,
        domain: list[str],
        adaptiveness: dict[str, float],
    ) -> float:
        """CAI impact for T-run / homopolymer constraints.

        T-run constraints only affect codons with multiple T bases.
        The impact is typically very small since the constraint only
        excludes specific codons that create long T runs.
        """
        # Count T-rich codons (2+ T bases)
        t_rich = [c for c in domain if c.count("T") >= 2]
        if not t_rich:
            return 0.0

        best_t = self._best_cai_in_domain(t_rich, adaptiveness)
        best_all = self._best_cai_in_domain(domain, adaptiveness)

        # If T-rich codons are the best, there is some impact
        return max(0.0, best_t - best_all) * 0.1

    def _build_default_sequence(self, model: CSPModel) -> str:
        """Build a default DNA sequence using highest-CAI codons.

        Used when no current sequence is provided for CAI impact
        estimation.
        """
        codons: list[str] = []
        for pos in range(model.length):
            domain = model.codon_domains.get(pos, [])
            if not domain:
                codons.append("ATG")  # Fallback
                continue
            # Pick the codon with highest adaptiveness
            best_codon = max(
                domain,
                key=lambda c: self._adaptiveness.get(c, 0.0),
            )
            codons.append(best_codon)
        return "".join(codons)

    def compute_cai(self, sequence: str) -> float:
        """Compute the CAI value for a sequence using the resolver's adaptiveness table.

        Parameters
        ----------
        sequence:
            DNA sequence to compute CAI for.

        Returns
        -------
        float
            CAI value in [0.0, 1.0].
        """
        return self._cai_constraint.cai(sequence)

    def compute_cai_loss(
        self,
        sequence_before: str,
        sequence_after: str,
    ) -> float:
        """Compute the CAI loss between two sequences.

        Parameters
        ----------
        sequence_before:
            The original sequence (before resolution).
        sequence_after:
            The modified sequence (after resolution).

        Returns
        -------
        float
            CAI loss (non-negative).  Zero means no CAI change.
        """
        cai_before = self.compute_cai(sequence_before)
        cai_after = self.compute_cai(sequence_after)
        return max(0.0, cai_before - cai_after)
