"""
BioCompiler Constraint Scoring & Enforcement
==============================================

Implements scoring mechanics for soft constraint evaluation, constraint
enforcement scoring with priority-weighted penalties, and multi-objective
optimization of codon-optimized gene sequences.

The constraint enforcement scoring pipeline:
    1. Each constraint is checked via its ``.check()`` method
    2. Violations are penalised according to their ``ConstraintPriority``
       (CRITICAL → score 0.0; HIGH/MEDIUM/LOW → weighted penalties)
    3. Constraint weights (``weight`` field) modulate penalty magnitude
    4. A composite score in [0.0, 1.0] summarises overall satisfaction

The soft constraint scoring pipeline:
    1. Each soft constraint (CAI, CpG, mRNA dG) produces a raw score
    2. Scores are normalized to [0.0, 1.0] for comparability
    3. A weighted sum produces the overall optimization objective
    4. Pareto frontier analysis identifies trade-off solutions

Usage::

    from biocompiler.solver.scoring import ConstraintScorer

    scorer = ConstraintScorer()
    score = scorer.score_solution("ATGGGCTGA...", constraints)
    violations = scorer.rank_violations("ATGGGCTGA...", constraints)

    # Soft constraint scoring
    from biocompiler.solver.scoring import SoftConstraintScorer, ScoringResult

    soft_scorer = SoftConstraintScorer(config)
    result = soft_scorer.score_sequence(model, "ATGGGCTGA...")
    print(result.weighted_total, result.cai_score, result.cpg_score)

    # Pareto analysis across multiple candidate sequences
    candidates = [soft_scorer.score_sequence(model, seq) for seq in sequences]
    pareto = compute_pareto_frontier(candidates)
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

from .types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintViolation,
    SolverConfig,
)
from .constraints import (
    CSPModel,
    HardConstraint,
    MaximizeCAI,
    MinimizeCodonPairBias,
    MinimizeCpG,
    MinimizeMRNADG,
    SoftConstraint,
    CAI_LOG_EPSILON,
)
from ..constants import AA_TO_CODONS

logger = logging.getLogger(__name__)


__all__ = [
    # Data classes
    "ScoringResult",
    "CAIImpactResult",
    # Enforcement scorer
    "ConstraintScorer",
    # Soft constraint scorer
    "SoftConstraintScorer",
    # Pareto analysis
    "compute_pareto_frontier",
]


# ==============================================================================
# ConstraintScorer — priority-weighted enforcement scoring
# ==============================================================================

# Default priority mapping for well-known constraint types.
# Used by ConstraintScorer when a ConstraintSpec doesn't carry an explicit
# priority (e.g. when passed as a plain data object without .check()).
_DEFAULT_CONSTRAINT_PRIORITIES: dict[str, ConstraintPriority] = {
    "TranslationConstraint": ConstraintPriority.CRITICAL,
    "NoRestrictionSiteConstraint": ConstraintPriority.HIGH,
    "GCRangeConstraint": ConstraintPriority.MEDIUM,
    "NoCrypticSpliceConstraint": ConstraintPriority.MEDIUM,
    "NoCpGIslandConstraint": ConstraintPriority.MEDIUM,
    "NoATTTAMotifConstraint": ConstraintPriority.MEDIUM,
    "NoTRunConstraint": ConstraintPriority.LOW,
    "MaximizeCAI": ConstraintPriority.LOW,
    "MinimizeCpG": ConstraintPriority.LOW,
    "MinimizeMRNADG": ConstraintPriority.LOW,
}


class ConstraintScorer:
    """Priority-weighted constraint enforcement scorer.

    Scores a candidate DNA sequence against a list of constraints and
    produces:

    1. A **composite score** in [0.0, 1.0] reflecting how well all
       constraints are satisfied, weighted by each constraint's priority
       and weight.
    2. A **ranked violation list** sorted by severity (CRITICAL first,
       then by priority descending, then by severity descending).

    Scoring formula
    ---------------
    The composite score starts at 1.0.  For each violated constraint,
    a penalty is deducted:

        penalty = severity × weight × priority.penalty_weight

    If any **CRITICAL** constraint is violated, the composite score is
    immediately 0.0 regardless of other results — critical constraints
    represent infeasible solutions.

    For non-critical violations, the score is:

        score = max(0.0, 1.0 - total_penalty / normalisation_factor)

    where the normalisation factor ensures the result stays in [0, 1].

    Usage::

        scorer = ConstraintScorer()
        score = scorer.score_solution(sequence, constraint_specs)
        if score < 1.0:
            violations = scorer.rank_violations(sequence, constraint_specs)
            for v in violations:
                print(f"  {v.constraint_name} ({v.priority.name}): {v.description}")
    """

    # Normalisation: the maximum possible penalty per constraint when
    # severity=1.0, weight=1.0, and priority=MEDIUM.  Used to keep
    # the composite score in [0, 1] for typical violation counts.
    _PENALTY_NORMALISATION: float = 3.0  # matches MEDIUM.penalty_weight

    def score_solution(
        self,
        sequence: str,
        constraints: list[ConstraintSpec],
    ) -> float:
        """Score a candidate sequence against constraints, returning 0.0–1.0.

        Each constraint is checked via its ``.check()`` method if available
        (on ``HardConstraint`` / ``SoftConstraint`` instances).  For plain
        ``ConstraintSpec`` objects that lack ``.check()``, the constraint is
        assumed satisfied (a warning is logged).

        Parameters
        ----------
        sequence:
            Candidate DNA sequence to evaluate.
        constraints:
            List of :class:`ConstraintSpec` (or subclass) instances to
            check against.

        Returns
        -------
        float
            Composite score in [0.0, 1.0].
            - 1.0 = all constraints satisfied
            - 0.0 = at least one CRITICAL constraint violated
            - Between = partial satisfaction with weighted penalties
        """
        if not sequence:
            return 0.0

        if not constraints:
            return 1.0

        has_critical_violation = False
        total_penalty = 0.0

        for spec in constraints:
            satisfied = self._check_constraint(spec, sequence)
            if satisfied:
                continue

            priority = spec.priority
            weight = spec.weight
            severity = self._estimate_severity(spec, sequence)

            if priority == ConstraintPriority.CRITICAL:
                has_critical_violation = True
                # Don't break — we want to log all violations

            penalty = severity * weight * priority.penalty_weight
            total_penalty += penalty

            logger.debug(
                "Constraint '%s' violated: priority=%s, weight=%.2f, "
                "severity=%.2f, penalty=%.4f",
                spec.name, priority.name, weight, severity, penalty,
            )

        if has_critical_violation:
            logger.info(
                "Composite score = 0.0 (CRITICAL constraint violated)"
            )
            return 0.0

        # Normalise: divide by number of constraints × normalisation factor
        # so that a single MEDIUM violation with severity=1.0 reduces
        # score by ~1/N (proportional impact).
        n = len(constraints)
        normalisation = n * self._PENALTY_NORMALISATION if n > 0 else 1.0
        score = max(0.0, 1.0 - total_penalty / normalisation)

        logger.info("Composite score = %.4f (total_penalty=%.4f)", score, total_penalty)
        return score

    def rank_violations(
        self,
        sequence: str,
        constraints: list[ConstraintSpec],
    ) -> list[ConstraintViolation]:
        """Return violations ranked by severity (most severe first).

        Checks each constraint against the sequence and returns a list
        of :class:`ConstraintViolation` objects for those that are not
        satisfied.  The list is sorted by:

        1. Priority (CRITICAL > HIGH > MEDIUM > LOW)
        2. Severity (higher severity first)
        3. Constraint name (alphabetical, for deterministic ordering)

        Parameters
        ----------
        sequence:
            Candidate DNA sequence to evaluate.
        constraints:
            List of :class:`ConstraintSpec` instances to check.

        Returns
        -------
        list[ConstraintViolation]
            Ranked violations, most severe first.  Empty if all
            constraints are satisfied.
        """
        if not sequence or not constraints:
            return []

        violations: list[ConstraintViolation] = []

        for spec in constraints:
            satisfied = self._check_constraint(spec, sequence)
            if satisfied:
                continue

            severity = self._estimate_severity(spec, sequence)
            strictness = spec.strictness
            description = (
                f"{strictness.value.capitalize()} constraint "
                f"'{spec.name}' is not satisfied "
                f"(priority={spec.priority.name}, weight={spec.weight:.2f})."
            )

            violation = ConstraintViolation(
                constraint_name=spec.name,
                constraint_type=strictness,
                description=description,
                severity=severity,
                priority=spec.priority,
                weight=spec.weight,
            )
            violations.append(violation)

        # Sort: CRITICAL first, then HIGH, MEDIUM, LOW;
        # within same priority, higher severity first; then by name.
        violations.sort(key=lambda v: (v.priority.rank, -v.severity, v.constraint_name))

        logger.info(
            "Ranked %d violation(s): %s",
            len(violations),
            [f"{v.constraint_name}({v.priority.name})" for v in violations],
        )
        return violations

    # ── Private helpers ───────────────────────────────────────────────

    @staticmethod
    def _check_constraint(spec: ConstraintSpec, sequence: str) -> bool:
        """Check whether a constraint is satisfied for the given sequence.

        Uses the constraint's ``.check()`` method if available (duck-typed).
        For plain ``ConstraintSpec`` objects without ``.check()``, logs a
        warning and assumes the constraint is satisfied.

        Parameters
        ----------
        spec:
            Constraint specification to check.
        sequence:
            Candidate DNA sequence.

        Returns
        -------
        bool
            True if the constraint is satisfied or cannot be checked.
        """
        check_method = getattr(spec, "check", None)
        if check_method is not None and callable(check_method):
            try:
                return bool(check_method(sequence))
            except Exception as exc:
                logger.warning(
                    "Constraint '%s' .check() raised %s: %s",
                    spec.name, type(exc).__name__, exc,
                )
                return False  # Treat errors as violations

        # Plain ConstraintSpec without .check() — assume satisfied
        logger.debug(
            "ConstraintSpec '%s' lacks .check() — assuming satisfied",
            spec.name,
        )
        return True

    @staticmethod
    def _estimate_severity(spec: ConstraintSpec, sequence: str) -> float:
        """Estimate violation severity for a constraint.

        Uses the constraint's ``.violated_positions()`` method if available
        to compute the fraction of the sequence affected.  Falls back to a
        fixed severity based on constraint strictness.

        Parameters
        ----------
        spec:
            The violated constraint specification.
        sequence:
            The candidate DNA sequence.

        Returns
        -------
        float
            Severity in [0.0, 1.0].
        """
        violated_positions_method = getattr(spec, "violated_positions", None)
        if violated_positions_method is not None and callable(violated_positions_method):
            try:
                positions = violated_positions_method(sequence)
                if not sequence:
                    return 1.0
                if not positions:
                    return 0.5  # Violated but no positions reported
                return min(len(set(positions)) / len(sequence), 1.0)
            except Exception:
                pass  # Fall through to default

        # Default severity based on strictness
        if spec.strictness == ConstraintStrictness.HARD:
            return 0.8
        return 0.5


# ==============================================================================
# Normalization constants
# ==============================================================================

# CAI is already in [0.0, 1.0] by definition (geometric mean of [0,1] values).
# For CpG, we normalize the negated count using a sigmoid-like transform
# so that 0 CpG → 1.0 and many CpG → 0.0.
# For mRNA dG, we normalize the negated |dG| using a sigmoid-like transform
# so that weak structure (|dG| near 0) → 1.0 and strong structure → 0.0.

_CPG_NORMALIZATION_SCALE: float = 10.0
"""Sigmoid scale for CpG normalization: controls how quickly score drops with CpG count."""

_MRNA_DG_NORMALIZATION_SCALE: float = 30.0
"""Sigmoid scale for mRNA dG normalization: controls normalization sensitivity."""

_CPB_NORMALIZATION_SCALE: float = 10.0
"""Sigmoid scale for codon pair bias normalization: maps CPB to [0, 1]."""

# ==============================================================================
# Calibrated default weights
# ==============================================================================

# CAI is the dominant soft objective — codon optimality is the primary driver
# of expression level.  Weight 10.0 ensures the solver prioritises codon
# adaptation over secondary concerns.
_DEFAULT_CAI_WEIGHT: float = 10.0

# CpG weight is GC-dependent.  High-GC targets inherently have more CpG
# dinucleotides, so the weight should be reduced to avoid over-penalising
# unavoidable CpG occurrences.  Low-GC targets can afford stronger CpG
# avoidance.
_DEFAULT_CPG_WEIGHT_HIGH_GC: float = 0.5   # For target GC > 60%
_DEFAULT_CPG_WEIGHT_LOW_GC: float = 2.0    # For target GC < 40%
_DEFAULT_CPG_WEIGHT_MID_GC: float = 1.0   # For 40% <= GC <= 60%

# Codon pair bias has a moderate but real effect on translation efficiency.
# Weight 3.0 reflects that CPB is less dominant than CAI but more important
# than the old default of 0.2.
_DEFAULT_CPB_WEIGHT: float = 3.0

# mRNA stability weight differs by domain of life:
# - Prokaryotes: mRNA is short-lived regardless; moderate weight (1.0)
#   acknowledges the role of 5' structure in ribosome binding without
#   over-optimising for stability that doesn't persist.
# - Eukaryotes: mRNA stability is a major determinant of expression level;
#   high weight (5.0) ensures the solver prioritises 5' accessibility.
_DEFAULT_PROKARYOTE_MRNA_DG_WEIGHT: float = 1.0
_DEFAULT_EUKARYOTE_MRNA_DG_WEIGHT: float = 5.0


# ==============================================================================
# ScoringResult dataclass
# ==============================================================================

@dataclass
class ScoringResult:
    """Result of scoring a candidate DNA sequence against soft constraints.

    Contains both individual normalized scores (each in [0.0, 1.0]) and the
    weighted total.  Individual scores are normalized for cross-objective
    comparability — raw scores from different constraint types have
    different scales and cannot be compared directly.

    Attributes:
        cai_score: Normalized CAI score in [0.0, 1.0].  Higher = better
            codon adaptation.  CAI is naturally in [0, 1] so no additional
            normalization is needed.
        cpg_score: Normalized CpG avoidance score in [0.0, 1.0].  Higher =
            fewer CpG dinucleotides.  Derived from sigmoid transform of
            the negated CpG count.
        mrna_dg_score: Normalized mRNA stability score in [0.0, 1.0].
            Higher = weaker 5' secondary structure (better for translation
            initiation).  Derived from sigmoid transform of |dG|.
        cpb_score: Normalized codon pair bias score in [0.0, 1.0].  Higher =
            more over-represented (favoured) codon pairs.  Derived from
            sigmoid transform of the mean CPB score.
        weighted_total: Weighted combination of individual scores using
            weights from SolverConfig.  Higher = better overall.
        individual_scores: Dictionary mapping constraint names to their
            raw (un-normalized) score values as returned by the constraint's
            ``score()`` method.
    """

    cai_score: float
    cpg_score: float
    mrna_dg_score: float
    cpb_score: float = 0.0
    weighted_total: float = 0.0
    individual_scores: dict[str, float] = field(default_factory=dict)

    # Backward compatibility: allow construction without cpb_score
    def __init__(
        self,
        cai_score: float = 0.0,
        cpg_score: float = 0.0,
        mrna_dg_score: float = 0.0,
        cpb_score: float = 0.0,
        weighted_total: float = 0.0,
        individual_scores: Optional[dict[str, float]] = None,
    ) -> None:
        self.cai_score = cai_score
        self.cpg_score = cpg_score
        self.mrna_dg_score = mrna_dg_score
        self.cpb_score = cpb_score
        self.weighted_total = weighted_total
        self.individual_scores = individual_scores if individual_scores is not None else {}

    def __repr__(self) -> str:
        return (
            f"ScoringResult(cai={self.cai_score:.4f}, cpg={self.cpg_score:.4f}, "
            f"mrna_dg={self.mrna_dg_score:.4f}, cpb={self.cpb_score:.4f}, "
            f"total={self.weighted_total:.4f})"
        )


# ==============================================================================
# CAIImpactResult dataclass
# ==============================================================================

@dataclass
class CAIImpactResult:
    """Result of CAI impact analysis for alternative codons at a position.

    Captures the CAI delta for each alternative codon, the current codon,
    and the best (least-CAI-damaging) alternative.

    Attributes:
        position: Codon position (0-based) that was analyzed.
        current_codon: The codon currently at this position.
        cai_deltas: Dict mapping alternative codon -> CAI delta
            (positive = improvement, negative = loss).
        best_codon: The alternative codon with the highest CAI delta
            (least damaging, or most beneficial).  None if empty.
        best_delta: The CAI delta of *best_codon*.
    """

    position: int
    current_codon: str
    cai_deltas: dict[str, float]
    best_codon: Optional[str] = None
    best_delta: float = 0.0

    def __post_init__(self) -> None:
        if self.best_codon is None and self.cai_deltas:
            best = max(self.cai_deltas.items(), key=lambda kv: kv[1])
            self.best_codon = best[0]
            self.best_delta = best[1]

    def __repr__(self) -> str:
        return (
            f"CAIImpactResult(pos={self.position}, current={self.current_codon}, "
            f"best={self.best_codon}, delta={self.best_delta:+.6f})"
        )


# ==============================================================================
# GC-aware CpG weight helper
# ==============================================================================

def _gc_aware_cpg_weight(target_gc: float) -> float:
    """Compute CpG weight scaled inversely with target GC content.

    High-GC targets (>60%) inherently contain more CpG dinucleotides,
    making aggressive CpG avoidance counterproductive.  Low-GC targets
    (<40%) have fewer unavoidable CpG occurrences, so stronger avoidance
    is feasible.

    Edge cases:
    - GC = 0: No G or C nucleotides → CpG dinucleotides are impossible.
      Weight is 0.0 since there is nothing to minimize.
    - GC = 1: Only G and C nucleotides → CpG dinucleotides are maximally
      present and largely unavoidable.  Weight is set to the high-GC
      value (0.5) to avoid over-penalising unavoidable occurrences.

    Args:
        target_gc: Target GC content as a fraction (e.g. 0.55 for 55%).

    Returns:
        CpG weight: 0.5 for GC > 60%, 2.0 for GC < 40%, 1.0 otherwise.
        Returns 0.0 for the degenerate cases GC <= 0 or GC >= 1.
    """
    # Edge cases: degenerate GC content
    if target_gc <= 0.0:
        # No G or C nucleotides → CpG impossible, no need to minimize
        return 0.0
    if target_gc >= 1.0:
        # Only G and C → CpG maximally present and unavoidable;
        # use high-GC weight to avoid penalising the inevitable
        return _DEFAULT_CPG_WEIGHT_HIGH_GC   # 0.5

    if target_gc > 0.60:
        return _DEFAULT_CPG_WEIGHT_HIGH_GC   # 0.5
    elif target_gc < 0.40:
        return _DEFAULT_CPG_WEIGHT_LOW_GC    # 2.0
    else:
        return _DEFAULT_CPG_WEIGHT_MID_GC    # 1.0


# ==============================================================================
# SoftConstraintScorer
# ==============================================================================

class SoftConstraintScorer:
    """Scores candidate DNA sequences against soft constraints.

    Evaluates a candidate sequence using the soft constraints (CAI, CpG,
    mRNA dG, and optionally codon pair bias) from the CSP model, normalizes
    each score to [0, 1], and computes a weighted total using calibrated
    SolverConfig weights.

    Weight calibration is applied during initialization to ensure proper
    relative importance of objectives:

    - CAI weight = 10.0 (dominant soft objective)
    - CpG weight = GC-dependent (0.5–2.0, inversely scaled with target GC)
    - CPB weight = 3.0 (moderate, affects translation efficiency)
    - mRNA dG weight = organism-dependent (1.0 prokaryote / 5.0 eukaryote)

    The scorer is stateless with respect to sequences — it can score any
    number of candidate sequences against the same or different models.

    Attributes:
        config: Solver configuration providing objective weights.

    Usage::

        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, dna_sequence)
        if result.weighted_total > best_total:
            best = result
    """

    def __init__(self, config: SolverConfig) -> None:
        """Initialize the scorer with solver configuration.

        Applies weight calibration on top of the SolverConfig defaults
        to ensure proper relative weighting of soft objectives:

        - CAI weight is set to 10.0 (dominant soft objective)
        - CpG weight is scaled by target GC content
        - CPB weight is set to 3.0 (moderate, affects translation efficiency)
        - mRNA dG weight is set based on organism domain (prokaryote vs eukaryote)

        Args:
            config: Solver configuration providing CAI, CpG, mRNA dG, and
                CPB weights.  The weights are calibrated upon initialization
                to ensure proper relative importance.
        """
        self._config = config
        self._calibrate_weights()
        logger.debug(
            "SoftConstraintScorer initialized with calibrated weights: "
            "cai=%.2f, cpg=%.2f, mrna_dg=%.2f, cpb=%.2f",
            self._config.cai_weight,
            self._config.cpg_weight,
            self._config.mrna_dg_weight,
            self._config.codon_pair_bias_weight,
        )

    def adjust_weights_for_organism(self, organism: str) -> None:
        """Dynamically adjust scoring weights based on the target organism.

        Inspects the organism's domain (prokaryote vs eukaryote) and
        adjusts **all** scorer weights accordingly:

        - **CAI weight** → 10.0 for both domains (dominant objective).
        - **CPB weight** → 3.0 for both domains (moderate, affects
          translation efficiency).
        - **Prokaryotes**: CpG weight is zeroed out (no DNA methylation),
          and mRNA dG weight is set to a moderate value (1.0) since
          prokaryotic mRNA is short-lived regardless.
        - **Eukaryotes**: CpG weight is scaled by target GC content
          (high GC → low weight, low GC → high weight), and mRNA dG
          weight is set to a high value (5.0) since mRNA stability is
          a major determinant of expression level.

        This method mutates the scorer's config in-place so that
        subsequent calls to :meth:`score_sequence` use the adjusted
        weights.

        Parameters
        ----------
        organism:
            Organism identifier string (e.g. ``"E_coli_K12"``,
            ``"human"``).  Resolved via
            :func:`~biocompiler.organism_config.get_organism_config`.

        Examples
        --------
        >>> scorer = SoftConstraintScorer(SolverConfig())
        >>> scorer.adjust_weights_for_organism("E_coli_K12")
        >>> scorer.config.cpg_weight
        0.0
        >>> scorer.config.mrna_dg_weight  # moderate for prokaryotes
        1.0
        >>> scorer.config.cai_weight  # dominant for all organisms
        10.0
        """
        from ..organism_config import get_organism_config

        org_config = get_organism_config(organism)

        old_cai = self._config.cai_weight
        old_cpg = self._config.cpg_weight
        old_mrna = self._config.mrna_dg_weight
        old_cpb = self._config.codon_pair_bias_weight

        # ── Domain-independent weights ──────────────────────────────
        # CAI is the dominant soft objective for ALL organisms.
        self._config.cai_weight = _DEFAULT_CAI_WEIGHT
        # CPB has moderate weight for ALL organisms.
        self._config.codon_pair_bias_weight = _DEFAULT_CPB_WEIGHT

        # ── Prokaryote-specific adjustments ────────────────────────────
        if not org_config.is_eukaryote:
            # CpG weight: zero for prokaryotes (no DNA methylation)
            self._config.cpg_weight = 0.0
            self._config.avoid_cpg = False
            logger.info(
                "Prokaryotic organism %r detected: zeroing cpg_weight "
                "(was %.2f) and disabling avoid_cpg",
                organism, old_cpg,
            )

            # mRNA dG weight: moderate for prokaryotes (short-lived mRNA)
            self._config.mrna_dg_weight = _DEFAULT_PROKARYOTE_MRNA_DG_WEIGHT
            logger.info(
                "Prokaryotic organism %r: setting mrna_dg_weight to "
                "%.2f (moderate, short-lived mRNA) (was %.2f)",
                organism, _DEFAULT_PROKARYOTE_MRNA_DG_WEIGHT, old_mrna,
            )

        # ── Eukaryote-specific adjustments ─────────────────────────────
        else:
            # mRNA dG weight: high for eukaryotes (stability matters)
            self._config.mrna_dg_weight = _DEFAULT_EUKARYOTE_MRNA_DG_WEIGHT
            logger.info(
                "Eukaryotic organism %r: setting mrna_dg_weight to "
                "%.2f (high, stability matters) (was %.2f)",
                organism, _DEFAULT_EUKARYOTE_MRNA_DG_WEIGHT, old_mrna,
            )

            # CpG weight: scale by target GC content
            target_gc = (org_config.gc_target_lo + org_config.gc_target_hi) / 2.0
            cpg_weight = _gc_aware_cpg_weight(target_gc)
            self._config.cpg_weight = cpg_weight
            logger.info(
                "Eukaryotic organism %r: GC-aware cpg_weight=%.2f "
                "(target GC=%.2f) (was %.2f)",
                organism, cpg_weight, target_gc, old_cpg,
            )

        # Adjust mRNA dG weight based on the organism's degradation model
        if org_config.mrna_degradation_model == "none":
            self._config.mrna_dg_weight = 0.0
            logger.info(
                "Organism %r has no mRNA degradation model: zeroing "
                "mrna_dg_weight (was %.2f)",
                organism, old_mrna,
            )

        # Update GC bounds from organism config
        if org_config.gc_target_lo != org_config.gc_target_hi:
            self._config.gc_lo = org_config.gc_target_lo
            self._config.gc_hi = org_config.gc_target_hi
            logger.debug(
                "Organism %r: updated GC range to [%.2f, %.2f]",
                organism, org_config.gc_target_lo, org_config.gc_target_hi,
            )

        logger.info(
            "Weight adjustment for organism %r: "
            "cai_weight %.2f -> %.2f, cpg_weight %.2f -> %.2f, "
            "mrna_dg_weight %.2f -> %.2f, cpb_weight %.2f -> %.2f",
            organism,
            old_cai, self._config.cai_weight,
            old_cpg, self._config.cpg_weight,
            old_mrna, self._config.mrna_dg_weight,
            old_cpb, self._config.codon_pair_bias_weight,
        )

    @property
    def config(self) -> SolverConfig:
        """The solver configuration used for weighting."""
        return self._config

    def score_sequence(self, model: CSPModel, sequence: str) -> ScoringResult:
        """Score a candidate DNA sequence against the model's soft constraints.

        Extracts each soft constraint from the model, computes its raw score
        and normalized score, then combines them into a weighted total.

        Args:
            model: The CSP model containing soft constraints and configuration.
            sequence: The candidate DNA sequence to score.  Must be the same
                length as the model's protein sequence * 3.

        Returns:
            ScoringResult with normalized individual scores and weighted total.

        Raises:
            ValueError: If the sequence length doesn't match the model's
                expected length (protein_length * 3).
        """
        expected_len = len(model.protein) * 3
        if len(sequence) != expected_len:
            raise ValueError(
                f"Sequence length ({len(sequence)}) doesn't match "
                f"expected length ({expected_len} = {len(model.protein)} aa * 3)"
            )

        # ── Compute raw scores from model's soft constraints ──────────
        individual_scores: dict[str, float] = {}
        for sc in model.soft_constraints:
            individual_scores[sc.name] = sc.score(sequence)

        logger.debug(
            "Raw scores for sequence (len=%d): %s",
            len(sequence),
            individual_scores,
        )

        # ── Normalize each score to [0.0, 1.0] ───────────────────────
        cai_score = self._normalize_cai(individual_scores, model, sequence)
        cpg_score = self._normalize_cpg(individual_scores)
        mrna_dg_score = self._normalize_mrna_dg(individual_scores, model, sequence)
        cpb_score = self._normalize_cpb(individual_scores)

        # ── Compute weighted total ────────────────────────────────────
        weights = self._build_weight_map()
        normalized_scores: dict[str, float] = {
            "MaximizeCAI": cai_score,
            "MinimizeCpG": cpg_score,
            "MinimizeMRNADG": mrna_dg_score,
            "MinimizeCodonPairBias": cpb_score,
        }
        weighted_total = self.compute_weighted_score(normalized_scores, weights)

        result = ScoringResult(
            cai_score=cai_score,
            cpg_score=cpg_score,
            mrna_dg_score=mrna_dg_score,
            cpb_score=cpb_score,
            weighted_total=weighted_total,
            individual_scores=individual_scores,
        )

        logger.debug("ScoringResult: %s", result)
        return result

    def compute_weighted_score(
        self,
        scores: dict[str, float],
        weights: dict[str, float],
    ) -> float:
        """Compute a weighted sum of normalized scores.

        Each score is multiplied by its corresponding weight, and the
        products are summed.  Constraint names not present in the weights
        dict are assigned a weight of 0.0 (excluded from the total).

        Args:
            scores: Dictionary mapping constraint names to normalized
                score values (typically in [0.0, 1.0]).
            weights: Dictionary mapping constraint names to weight values.
                Higher weight = more important.

        Returns:
            Weighted sum (float).  Higher = better.  The value depends on
            the magnitude of the weights and scores; for normalized scores
            in [0, 1] and weights summing to ~1.8 (default config), the
            total will typically be in [0.0, 1.8].
        """
        total = 0.0
        for name, score_val in scores.items():
            w = weights.get(name, 0.0)
            total += w * score_val
        return total

    # ── CAI-aware scoring ───────────────────────────────────────────────

    def score_cai_impact(
        self,
        model: CSPModel,
        sequence: str,
        position: int,
        alternative_codons: list[str],
    ) -> dict[str, float]:
        """Compute the CAI impact of each alternative codon at a position.

        Positive delta = CAI improvement; negative = CAI loss.
        This enables the constraint resolver to pick the least-CAI-damaging fix.

        Args:
            model: CSP model with soft constraints (including MaximizeCAI).
            sequence: Current candidate DNA sequence.
            position: 0-based codon position.
            alternative_codons: List of candidate codons to evaluate.

        Returns:
            Dict mapping codon -> CAI delta.
        """
        n = len(model.protein)
        if position < 0 or position >= n:
            raise ValueError(
                f"Codon position {position} out of range "
                f"[0, {n - 1}] for protein of length {n}"
            )

        cai_constraint = None
        for sc in model.soft_constraints:
            if isinstance(sc, MaximizeCAI):
                cai_constraint = sc
                break

        if cai_constraint is None:
            return {codon: 0.0 for codon in alternative_codons}

        adaptiveness = cai_constraint.adaptiveness
        current_codon = sequence[position * 3 : position * 3 + 3]
        current_w = adaptiveness.get(current_codon, CAI_LOG_EPSILON)
        current_log_w = math.log(max(current_w, CAI_LOG_EPSILON))

        result: dict[str, float] = {}
        for codon in alternative_codons:
            alt_w = adaptiveness.get(codon, CAI_LOG_EPSILON)
            alt_log_w = math.log(max(alt_w, CAI_LOG_EPSILON))
            delta = (alt_log_w - current_log_w) / n
            result[codon] = delta

        logger.debug(
            "CAI impact at codon position %d (current=%s): %s",
            position, current_codon,
            {c: f"{d:+.6f}" for c, d in result.items()},
        )
        return result

    # ── Private helpers ───────────────────────────────────────────────

    def _build_weight_map(self) -> dict[str, float]:
        """Build the constraint-name → weight mapping from config.

        Returns:
            Dictionary with keys 'MaximizeCAI', 'MinimizeCpG',
            'MinimizeMRNADG', 'MinimizeCodonPairBias' and their
            respective config weights.
        """
        weights = {
            "MaximizeCAI": self._config.cai_weight,
            "MinimizeCpG": self._config.cpg_weight,
            "MinimizeMRNADG": self._config.mrna_dg_weight,
        }
        # Include CPB weight only if codon pair bias optimisation is enabled
        if self._config.optimize_codon_pair_bias:
            weights["MinimizeCodonPairBias"] = self._config.codon_pair_bias_weight
        return weights

    def _normalize_cai(
        self,
        individual_scores: dict[str, float],
        model: CSPModel,
        sequence: str,
    ) -> float:
        """Normalize the CAI score to [0.0, 1.0].

        CAI is the geometric mean of relative adaptiveness values and is
        naturally in [0.0, 1.0].  We extract it directly from the
        MaximizeCAI constraint's ``cai()`` method rather than converting
        the log-CAI raw score.

        Args:
            individual_scores: Raw scores dict (not used for CAI directly).
            model: CSP model (to find the MaximizeCAI constraint).
            sequence: DNA sequence to compute CAI for.

        Returns:
            CAI value in [0.0, 1.0].
        """
        for sc in model.soft_constraints:
            if isinstance(sc, MaximizeCAI):
                return max(0.0, min(1.0, sc.cai(sequence)))

        # Fallback: if no MaximizeCAI constraint found, compute from raw score
        raw = individual_scores.get("MaximizeCAI", 0.0)
        n = len(model.protein)
        if n == 0:
            return 0.0
        cai = math.exp(raw / n)
        return max(0.0, min(1.0, cai))

    def _normalize_cpg(self, individual_scores: dict[str, float]) -> float:
        """Normalize the CpG score to [0.0, 1.0].

        The raw MinimizeCpG score is the negated CpG count (higher = better).
        We convert this to [0, 1] using a sigmoid-like transform:

            normalized = 1.0 / (1.0 + exp(cpg_count / scale))

        This maps 0 CpG → ~1.0 and many CpG → ~0.0.

        Args:
            individual_scores: Raw scores dict containing 'MinimizeCpG'.

        Returns:
            Normalized CpG avoidance score in [0.0, 1.0].
        """
        raw = individual_scores.get("MinimizeCpG", 0.0)
        # raw = -cpg_count, so cpg_count = -raw
        cpg_count = -raw

        if cpg_count <= 0:
            return 1.0

        # Sigmoid normalization: 1 / (1 + exp(count / scale))
        try:
            normalized = 1.0 / (1.0 + math.exp(cpg_count / _CPG_NORMALIZATION_SCALE))
        except OverflowError:
            # Very large cpg_count → exp overflows → normalized ≈ 0
            normalized = 0.0

        return max(0.0, min(1.0, normalized))

    def _normalize_mrna_dg(
        self,
        individual_scores: dict[str, float],
        model: CSPModel,
        sequence: str,
    ) -> float:
        """Normalize the mRNA dG score to [0.0, 1.0].

        The raw MinimizeMRNADG score is the negated |dG| (higher = better =
        weaker structure).  We normalize using a sigmoid-like transform:

            normalized = 1.0 / (1.0 + exp(|dG| / scale))

        This maps weak structure (|dG| near 0) → ~1.0 and strong structure
        (very negative dG, large |dG|) → ~0.0.

        Args:
            individual_scores: Raw scores dict containing 'MinimizeMRNADG'.
            model: CSP model (to find MinimizeMRNADG for dG extraction).
            sequence: DNA sequence being scored.

        Returns:
            Normalized mRNA stability score in [0.0, 1.0].
        """
        raw = individual_scores.get("MinimizeMRNADG", 0.0)
        # raw = -|dG|, so |dG| = -raw
        abs_dg = -raw

        if abs_dg <= 0:
            return 1.0

        # Sigmoid normalization
        try:
            normalized = 1.0 / (
                1.0 + math.exp(abs_dg / _MRNA_DG_NORMALIZATION_SCALE)
            )
        except OverflowError:
            normalized = 0.0

        return max(0.0, min(1.0, normalized))

    def _normalize_cpb(self, individual_scores: dict[str, float]) -> float:
        """Normalize the codon pair bias score to [0.0, 1.0].

        The raw MinimizeCodonPairBias score is the mean CPB across all
        codon pairs.  CPB values are typically in [-0.5, 0.5], where
        positive = over-represented (favoured) pairs and negative =
        under-represented (disfavoured) pairs.

        We normalize using a sigmoid-like transform:

            normalized = 1.0 / (1.0 + exp(-cpb * scale))

        This maps negative CPB (disfavoured) → ~0.0, zero CPB → 0.5,
        and positive CPB (favoured) → ~1.0.

        To prevent overflow for extreme CPB values (e.g. due to
        corrupted input or very short sequences), the raw CPB is
        clamped to [-2, 2] before the sigmoid is applied.  With
        ``_CPB_NORMALIZATION_SCALE = 10.0``, this bounds the exponent
        to [-20, 20], well within ``float`` range.

        Args:
            individual_scores: Raw scores dict containing
                'MinimizeCodonPairBias'.

        Returns:
            Normalized codon pair bias score in [0.0, 1.0].
        """
        raw = individual_scores.get("MinimizeCodonPairBias", 0.0)

        # Clamp CPB to [-2, 2] to prevent sigmoid overflow
        clamped_cpb = max(-2.0, min(2.0, raw))

        # Sigmoid normalization: maps CPB to (0, 1)
        # positive CPB → high score, negative CPB → low score
        normalized = 1.0 / (1.0 + math.exp(-clamped_cpb * _CPB_NORMALIZATION_SCALE))

        return max(0.0, min(1.0, normalized))

    def _calibrate_weights(self) -> None:
        """Calibrate scoring weights to ensure proper relative importance.

        Applies the following calibrations to the scorer's config:

        - **CAI weight** → 10.0: CAI is the dominant soft objective; codon
          optimality is the primary driver of expression level.
        - **CPG weight** → GC-dependent: high GC (>60%) = 0.5 (less
          avoidance, more unavoidable CpG), low GC (<40%) = 2.0 (strong
          avoidance feasible), mid-range = 1.0.
        - **CPB weight** → 3.0: Codon pair bias has a moderate but real
          effect on translation efficiency.
        - **mRNA dG weight** → organism-dependent: moderate (1.0) for
          prokaryotes (short-lived mRNA), high (5.0) for eukaryotes
          (stability critical for expression).

        This method is called during ``__init__`` and ensures the scorer
        uses calibrated defaults regardless of the SolverConfig's initial
        values.  Subsequent calls to :meth:`adjust_weights_for_organism`
        may further refine these weights.
        """
        # CAI: dominant soft objective
        old_cai = self._config.cai_weight
        self._config.cai_weight = _DEFAULT_CAI_WEIGHT

        # CpG: GC-dependent weight
        old_cpg = self._config.cpg_weight
        target_gc = (self._config.gc_lo + self._config.gc_hi) / 2.0
        self._config.cpg_weight = _gc_aware_cpg_weight(target_gc)

        # CPB: moderate weight for translation efficiency
        old_cpb = self._config.codon_pair_bias_weight
        self._config.codon_pair_bias_weight = _DEFAULT_CPB_WEIGHT

        # mRNA dG: organism-dependent weight
        old_mrna = self._config.mrna_dg_weight
        if self._config.is_eukaryotic:
            self._config.mrna_dg_weight = _DEFAULT_EUKARYOTE_MRNA_DG_WEIGHT
        else:
            self._config.mrna_dg_weight = _DEFAULT_PROKARYOTE_MRNA_DG_WEIGHT

        logger.info(
            "Weight calibration: cai %.2f -> %.2f, cpg %.2f -> %.2f (target GC=%.2f), "
            "cpb %.2f -> %.2f, mrna_dg %.2f -> %.2f (%s)",
            old_cai, self._config.cai_weight,
            old_cpg, self._config.cpg_weight, target_gc,
            old_cpb, self._config.codon_pair_bias_weight,
            old_mrna, self._config.mrna_dg_weight,
            "eukaryote" if self._config.is_eukaryotic else "prokaryote",
        )


# ==============================================================================
# Pareto frontier analysis
# ==============================================================================

def compute_pareto_frontier(results: list[ScoringResult]) -> list[ScoringResult]:
    """Identify Pareto-optimal solutions from a list of scoring results.

    A solution is Pareto-optimal (non-dominated) if no other solution is
    at least as good in ALL objectives and strictly better in at least one.
    In other words, a Pareto-optimal solution represents a trade-off where
    improving one objective would require sacrificing another.

    The four objectives (all to be maximized) are:
    - cai_score: Codon adaptation quality
    - cpg_score: CpG avoidance
    - mrna_dg_score: mRNA structural accessibility
    - cpb_score: Codon pair bias (over-represented pairs)

    Args:
        results: List of ScoringResult instances to analyze.  May be empty.

    Returns:
        List of Pareto-optimal ScoringResult instances, preserving the
        original order from the input.  Returns an empty list if the
        input is empty.

    Examples:
        >>> r1 = ScoringResult(cai_score=0.9, cpg_score=0.3, mrna_dg_score=0.5, cpb_score=0.6)
        >>> r2 = ScoringResult(cai_score=0.7, cpg_score=0.8, mrna_dg_score=0.4, cpb_score=0.7)
        >>> r3 = ScoringResult(cai_score=0.6, cpg_score=0.4, mrna_dg_score=0.3, cpb_score=0.4)
        >>> pareto = compute_pareto_frontier([r1, r2, r3])
        >>> len(pareto)
        2
        >>> # r3 is dominated by both r1 and r2, so only r1 and r2 survive
    """
    if not results:
        return []

    if len(results) == 1:
        return list(results)

    # Objectives to maximize (name → type hint)
    objectives: list[tuple[str, type[ScoringResult]]] = [
        ("cai_score", float),
        ("cpg_score", float),
        ("mrna_dg_score", float),
        ("cpb_score", float),
    ]

    def _dominates(a: ScoringResult, b: ScoringResult) -> bool:
        """Return True if *a* dominates *b*.

        a dominates b iff:
        - a is >= b in ALL objectives, AND
        - a is > b in AT LEAST ONE objective.
        """
        at_least_one_strictly_better = False
        for obj_name, _ in objectives:
            a_val = getattr(a, obj_name)
            b_val = getattr(b, obj_name)
            if a_val < b_val:
                # a is worse in this objective → does not dominate b
                return False
            if a_val > b_val:
                at_least_one_strictly_better = True
        return at_least_one_strictly_better

    # Identify Pareto-optimal solutions
    pareto_indices: list[int] = []
    n = len(results)

    for i in range(n):
        is_dominated = False
        for j in range(n):
            if i == j:
                continue
            if _dominates(results[j], results[i]):
                is_dominated = True
                break
        if not is_dominated:
            pareto_indices.append(i)

    pareto_results = [results[i] for i in pareto_indices]

    logger.debug(
        "Pareto frontier: %d of %d solutions are Pareto-optimal",
        len(pareto_results),
        n,
    )

    return pareto_results
