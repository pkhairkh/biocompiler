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
from typing import Sequence

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
    MinimizeCpG,
    MinimizeMRNADG,
    SoftConstraint,
    CAI_LOG_EPSILON,
)

logger = logging.getLogger(__name__)


__all__ = [
    # Data classes
    "ScoringResult",
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
        weighted_total: Weighted combination of individual scores using
            weights from SolverConfig.  Higher = better overall.
        individual_scores: Dictionary mapping constraint names to their
            raw (un-normalized) score values as returned by the constraint's
            ``score()`` method.
    """

    cai_score: float
    cpg_score: float
    mrna_dg_score: float
    weighted_total: float
    individual_scores: dict[str, float]

    def __repr__(self) -> str:
        return (
            f"ScoringResult(cai={self.cai_score:.4f}, cpg={self.cpg_score:.4f}, "
            f"mrna_dg={self.mrna_dg_score:.4f}, total={self.weighted_total:.4f})"
        )


# ==============================================================================
# SoftConstraintScorer
# ==============================================================================

class SoftConstraintScorer:
    """Scores candidate DNA sequences against soft constraints.

    Evaluates a candidate sequence using the three standard soft constraints
    (CAI, CpG, mRNA dG) from the CSP model, normalizes each score to [0, 1],
    and computes a weighted total using the SolverConfig weights.

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

        Args:
            config: Solver configuration providing CAI, CpG, and mRNA dG
                weights.  The weights control the relative importance of
                each objective in the weighted total.
        """
        self._config = config
        logger.debug(
            "SoftConstraintScorer initialized with weights: "
            "cai=%.2f, cpg=%.2f, mrna_dg=%.2f",
            config.cai_weight,
            config.cpg_weight,
            config.mrna_dg_weight,
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

        # ── Compute weighted total ────────────────────────────────────
        weights = self._build_weight_map()
        normalized_scores: dict[str, float] = {
            "MaximizeCAI": cai_score,
            "MinimizeCpG": cpg_score,
            "MinimizeMRNADG": mrna_dg_score,
        }
        weighted_total = self.compute_weighted_score(normalized_scores, weights)

        result = ScoringResult(
            cai_score=cai_score,
            cpg_score=cpg_score,
            mrna_dg_score=mrna_dg_score,
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

    # ── Private helpers ───────────────────────────────────────────────

    def _build_weight_map(self) -> dict[str, float]:
        """Build the constraint-name → weight mapping from config.

        Returns:
            Dictionary with keys 'MaximizeCAI', 'MinimizeCpG',
            'MinimizeMRNADG' and their respective config weights.
        """
        return {
            "MaximizeCAI": self._config.cai_weight,
            "MinimizeCpG": self._config.cpg_weight,
            "MinimizeMRNADG": self._config.mrna_dg_weight,
        }

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


# ==============================================================================
# Pareto frontier analysis
# ==============================================================================

def compute_pareto_frontier(results: list[ScoringResult]) -> list[ScoringResult]:
    """Identify Pareto-optimal solutions from a list of scoring results.

    A solution is Pareto-optimal (non-dominated) if no other solution is
    at least as good in ALL objectives and strictly better in at least one.
    In other words, a Pareto-optimal solution represents a trade-off where
    improving one objective would require sacrificing another.

    The three objectives (all to be maximized) are:
    - cai_score: Codon adaptation quality
    - cpg_score: CpG avoidance
    - mrna_dg_score: mRNA structural accessibility

    Args:
        results: List of ScoringResult instances to analyze.  May be empty.

    Returns:
        List of Pareto-optimal ScoringResult instances, preserving the
        original order from the input.  Returns an empty list if the
        input is empty.

    Examples:
        >>> r1 = ScoringResult(cai_score=0.9, cpg_score=0.3, mrna_dg_score=0.5, weighted_total=0.0, individual_scores={})
        >>> r2 = ScoringResult(cai_score=0.7, cpg_score=0.8, mrna_dg_score=0.4, weighted_total=0.0, individual_scores={})
        >>> r3 = ScoringResult(cai_score=0.6, cpg_score=0.4, mrna_dg_score=0.3, weighted_total=0.0, individual_scores={})
        >>> pareto = compute_pareto_frontier([r1, r2, r3])
        >>> len(pareto)
        2
        >>> # r3 is dominated by both r1 and r2, so only r1 and r2 survive
    """
    if not results:
        return []

    if len(results) == 1:
        return list(results)

    # Objectives to maximize (name → accessor function)
    objectives: list[tuple[str, type[ScoringResult]]] = [
        ("cai_score", float),
        ("cpg_score", float),
        ("mrna_dg_score", float),
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
