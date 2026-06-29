"""
BioCompiler Decision-Level Provenance — Granular Codon & Constraint Audit Trail
================================================================================

Provides fine-grained, decision-level provenance that captures *why* every codon
was chosen, what alternatives were considered, which constraints drove each
choice, and what tradeoffs were made.

This module complements the coarser-grained :mod:`biocompiler.provenance` module
(which records solver-level decisions) by operating at the individual codon
position level.  Each ``CodonDecision`` records the full set of alternatives
evaluated with their CAI/GC contributions and constraint violations, enabling
full reproducibility and debugging.  ``ConstraintDecision`` entries capture the
cost of each constraint in CAI terms, making tradeoffs explicit and auditable.

The ``DecisionProvenanceCollector`` is a builder that accumulates decisions
during an optimization run and produces an ``OptimizationDecisionTrail`` — a
complete, serializable audit trail tying input, configuration, every codon and
constraint decision, and final output together.  This trail can be persisted as
JSON alongside certificates so that any third party can audit *why* every codon
was chosen.

Usage::

    from biocompiler.provenance.decision_provenance import (
        CodonDecision, ConstraintDecision, DecisionProvenanceCollector,
    )

    collector = DecisionProvenanceCollector()
    collector.start_optimization(
        protein="MVLSPADKTN",
        organism="Homo_sapiens",
        constraints=["GCInRange", "NoCrypticSplice"],
        gene_name="eGFP",
        solver_backend="greedy",
    )

    # During optimization, record each codon choice
    collector.record_codon_decision(CodonDecision(
        position=0, amino_acid="M", original_codon=None,
        chosen_codon="ATG", alternatives_considered=[],
        constraint_reason="maximize_cai", confidence=1.0,
    ))

    # Record constraint outcomes
    collector.record_constraint_decision(ConstraintDecision(
        constraint_name="NoCrypticSplice", constraint_type="hard",
        action_taken="satisfied", positions_affected=[3, 7],
        tradeoff_description="Chose GTC over GTG to avoid cryptic donor",
        impact_on_cai=-0.003,
    ))

    # Finalize and serialize
    trail = collector.finalize(output_dna="ATGGTG...", cai=0.91, gc=0.54)
    print(collector.summary())
    collector.to_json("decision_trail.json")

Key components:
- CodonDecision: why a specific codon was chosen at a position, with full
  alternative analysis
- ConstraintDecision: how each constraint influenced the optimization, including
  tradeoff costs in CAI
- OptimizationDecisionTrail: complete end-to-end audit trail tying input,
  configuration, every codon and constraint decision, and final output together
- DecisionProvenanceCollector: builder that accumulates decisions during an
  optimization run and produces an OptimizationDecisionTrail

Design goals:
1. Every codon choice is recorded with the full set of alternatives evaluated,
   including CAI/GC contributions and constraint violations — enabling full
   reproducibility and debugging.
2. Constraint decisions capture the cost of each constraint in CAI terms,
   making tradeoffs explicit and auditable.
3. The iteration log records what happened at each solver iteration, enabling
   understanding of convergence behavior.
4. Serialization (to_dict / to_json) enables persistence, comparison across
   runs, and integration with certificate provenance metadata.
"""

from __future__ import annotations

import json
import logging
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "CodonDecision",
    "ConstraintDecision",
    "OptimizationDecisionTrail",
    "DecisionProvenanceCollector",
    "ProvenanceStore",
]


# ---------------------------------------------------------------------------
# CodonDecision
# ---------------------------------------------------------------------------

@dataclass
class CodonDecision:
    """Record of why a specific codon was chosen at a given position.

    Captures the full decision context: what was chosen, what alternatives
    were considered and why they were rejected, which constraint or
    optimization goal drove the choice, and how confident we are.

    Attributes:
        position: 0-indexed codon position in the protein sequence.
        amino_acid: 1-letter amino acid code (e.g. ``"V"`` for Valine).
        original_codon: The codon from the input DNA sequence, if this
            optimization started from an existing sequence.  ``None`` for
            de novo optimization.
        chosen_codon: The codon that was selected (e.g. ``"GTC"``).
        alternatives_considered: Ordered list of alternative codons that were
            evaluated, each with its CAI/GC contributions, constraint
            violations, and rejection reason.  Ordering reflects solver
            preference (best first).
        constraint_reason: Which constraint or optimization goal primarily
            drove this choice (e.g. ``"maximize_cai"``,
            ``"avoid_restriction_site:EcoRI"``, ``"gc_content"``).
        confidence: How confident we are that this is the right codon,
            on a 0–1 scale.  1.0 means no doubt; lower values indicate
            uncertainty due to conflicting constraints or marginal
            improvements.
        cai_impact: CAI delta caused by this codon choice versus the
            best-CAI codon for this amino acid.  Negative values indicate
            CAI was sacrificed to satisfy a constraint (e.g. GC content,
            restriction site avoidance).  Zero means no CAI cost (either
            the best-CAI codon was chosen, or CAI was not relevant).
    """

    position: int
    amino_acid: str
    original_codon: str | None
    chosen_codon: str
    alternatives_considered: list[dict[str, Any]]
    constraint_reason: str
    confidence: float
    # CAI-aware provenance: records the CAI cost when a codon is chosen
    # due to a constraint rather than pure CAI optimization.
    # Negative = CAI lost vs best-CAI codon; 0.0 = no CAI cost.
    cai_impact: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        """Serialize this decision to a JSON-compatible dict."""
        return {
            "position": self.position,
            "amino_acid": self.amino_acid,
            "original_codon": self.original_codon,
            "chosen_codon": self.chosen_codon,
            "alternatives_considered": [
                dict(alt) for alt in self.alternatives_considered
            ],
            "constraint_reason": self.constraint_reason,
            "confidence": self.confidence,
            "cai_impact": self.cai_impact,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> CodonDecision:
        """Deserialize a CodonDecision from a plain dict.

        Args:
            data: Dict with the same keys as the dataclass fields.

        Returns:
            A CodonDecision instance.

        Raises:
            ValueError: If required keys are missing.
        """
        required = {
            "position", "amino_acid", "original_codon", "chosen_codon",
            "alternatives_considered", "constraint_reason", "confidence",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize CodonDecision: missing keys {missing}"
            )
        return cls(
            position=int(data["position"]),
            amino_acid=data["amino_acid"],
            original_codon=data["original_codon"],
            chosen_codon=data["chosen_codon"],
            alternatives_considered=[
                dict(alt) for alt in data["alternatives_considered"]
            ],
            constraint_reason=data["constraint_reason"],
            confidence=float(data["confidence"]),
            cai_impact=float(data.get("cai_impact", 0.0)),
        )


# ---------------------------------------------------------------------------
# ConstraintDecision
# ---------------------------------------------------------------------------

_VALID_ACTIONS = {"satisfied", "relaxed", "conflicted", "overridden"}


@dataclass
class ConstraintDecision:
    """Record of how a constraint influenced the optimization.

    Captures the outcome of applying (or failing to fully apply) a constraint,
    including the positions it affected, a human-readable tradeoff description,
    and the cost in CAI terms.

    Attributes:
        constraint_name: Name of the constraint (e.g. ``"NoRestrictionSite"``,
            ``"GCInRange"``, ``"NoCrypticSplice"``).
        constraint_type: Category of the constraint (e.g. ``"hard"``,
            ``"soft"``, ``"preference"``).
        action_taken: What happened with this constraint — one of
            ``"satisfied"`` (fully met), ``"relaxed"`` (partially relaxed to
            find a feasible solution), ``"conflicted"`` (in conflict with
            another constraint), or ``"overridden"`` (explicitly overridden by
            user or higher-priority constraint).
        positions_affected: List of 0-indexed codon positions where this
            constraint influenced the codon choice.
        tradeoff_description: Human-readable explanation of the tradeoff,
            e.g. ``"Chose codon X over Y to avoid restriction site at
            position 42"``.
        impact_on_cai: How much this constraint cost us in CAI.  Negative
            values indicate a CAI penalty; zero means no CAI impact.
    """

    constraint_name: str
    constraint_type: str
    action_taken: str
    positions_affected: list[int]
    tradeoff_description: str
    impact_on_cai: float

    def __post_init__(self) -> None:
        """Validate action_taken is one of the allowed values."""
        if self.action_taken not in _VALID_ACTIONS:
            raise ValueError(
                f"Invalid action_taken '{self.action_taken}'. "
                f"Must be one of {_VALID_ACTIONS}"
            )

    def to_dict(self) -> dict[str, Any]:
        """Serialize this decision to a JSON-compatible dict."""
        return {
            "constraint_name": self.constraint_name,
            "constraint_type": self.constraint_type,
            "action_taken": self.action_taken,
            "positions_affected": list(self.positions_affected),
            "tradeoff_description": self.tradeoff_description,
            "impact_on_cai": self.impact_on_cai,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ConstraintDecision:
        """Deserialize a ConstraintDecision from a plain dict.

        Args:
            data: Dict with the same keys as the dataclass fields.

        Returns:
            A ConstraintDecision instance.

        Raises:
            ValueError: If required keys are missing or action_taken is invalid.
        """
        required = {
            "constraint_name", "constraint_type", "action_taken",
            "positions_affected", "tradeoff_description", "impact_on_cai",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize ConstraintDecision: missing keys {missing}"
            )
        return cls(
            constraint_name=data["constraint_name"],
            constraint_type=data["constraint_type"],
            action_taken=data["action_taken"],
            positions_affected=list(data["positions_affected"]),
            tradeoff_description=data["tradeoff_description"],
            impact_on_cai=float(data["impact_on_cai"]),
        )


# ---------------------------------------------------------------------------
# OptimizationDecisionTrail
# ---------------------------------------------------------------------------

def _get_biocompiler_version() -> str:
    """Return the current biocompiler version string."""
    try:
        from .. import __version__
        return __version__
    except ImportError:
        return "unknown"


@dataclass
class OptimizationDecisionTrail:
    """Complete end-to-end audit trail for a single optimization run.

    Ties together the input protein, solver configuration, every codon-level
    and constraint-level decision, iteration log, and final output.  Designed
    for persistence alongside certificates so that any third party can audit
    *why* every codon was chosen.

    Attributes:
        gene_name: Optional gene name (e.g. ``"eGFP"``, ``"HBB"``).
        input_protein: The original amino acid sequence submitted for
            optimization.
        output_dna: The optimized DNA sequence that was output.
        organism: Target organism identifier (e.g. ``"Homo_sapiens"``).
        solver_backend: Name of the solver backend used
            (e.g. ``"greedy"``, ``"z3"``, ``"ortools"``).
        seed: RNG seed used for this run.  ``None`` if no seed was provided.
        total_cai: Overall Codon Adaptation Index of the output sequence.
        total_gc: Overall GC content of the output sequence (0–1).
        codon_decisions: Ordered list of every CodonDecision produced during
            the run, one per codon position.
        constraint_decisions: Ordered list of every ConstraintDecision produced
            during the run, one per constraint that was active.
        iteration_log: Ordered list of dicts describing what happened at each
            optimization iteration (solver state, score changes, etc.).
        timestamp: ISO 8601 timestamp of when the optimization completed.
        version: Version of biocompiler that produced this trail.
    """

    gene_name: str | None
    input_protein: str
    output_dna: str
    organism: str
    solver_backend: str
    seed: int | None
    total_cai: float
    total_gc: float
    codon_decisions: list[CodonDecision]
    constraint_decisions: list[ConstraintDecision]
    iteration_log: list[dict[str, Any]]
    timestamp: str
    version: str

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict.

        Returns:
            Dict representation of the full optimization decision trail.
        """
        return {
            "gene_name": self.gene_name,
            "input_protein": self.input_protein,
            "output_dna": self.output_dna,
            "organism": self.organism,
            "solver_backend": self.solver_backend,
            "seed": self.seed,
            "total_cai": self.total_cai,
            "total_gc": self.total_gc,
            "codon_decisions": [d.to_dict() for d in self.codon_decisions],
            "constraint_decisions": [
                d.to_dict() for d in self.constraint_decisions
            ],
            "iteration_log": [dict(entry) for entry in self.iteration_log],
            "timestamp": self.timestamp,
            "version": self.version,
        }

    def to_json(self) -> str:
        """Serialize to an indented JSON string.

        Returns:
            JSON string representation.
        """
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> OptimizationDecisionTrail:
        """Deserialize from a plain dict.

        Args:
            data: Dict with the same keys as the dataclass fields.

        Returns:
            An OptimizationDecisionTrail instance.

        Raises:
            ValueError: If required keys are missing.
        """
        required = {
            "gene_name", "input_protein", "output_dna", "organism",
            "solver_backend", "seed", "total_cai", "total_gc",
            "codon_decisions", "constraint_decisions", "iteration_log",
            "timestamp", "version",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize OptimizationDecisionTrail: missing keys {missing}"
            )
        return cls(
            gene_name=data["gene_name"],
            input_protein=data["input_protein"],
            output_dna=data["output_dna"],
            organism=data["organism"],
            solver_backend=data["solver_backend"],
            seed=data["seed"],
            total_cai=float(data["total_cai"]),
            total_gc=float(data["total_gc"]),
            codon_decisions=[
                CodonDecision.from_dict(d) for d in data["codon_decisions"]
            ],
            constraint_decisions=[
                ConstraintDecision.from_dict(d)
                for d in data["constraint_decisions"]
            ],
            iteration_log=[dict(entry) for entry in data["iteration_log"]],
            timestamp=data["timestamp"],
            version=data["version"],
        )

    @classmethod
    def from_json(cls, json_str: str) -> OptimizationDecisionTrail:
        """Deserialize from a JSON string.

        Args:
            json_str: JSON string produced by :meth:`to_json`.

        Returns:
            An OptimizationDecisionTrail instance.
        """
        return cls.from_dict(json.loads(json_str))

    def __repr__(self) -> str:
        return (
            f"OptimizationDecisionTrail("
            f"gene={self.gene_name!r}, "
            f"organism={self.organism!r}, "
            f"backend={self.solver_backend!r}, "
            f"codons={len(self.codon_decisions)}, "
            f"constraints={len(self.constraint_decisions)}, "
            f"cai={self.total_cai:.4f}, "
            f"gc={self.total_gc:.3f})"
        )


# ---------------------------------------------------------------------------
# DecisionProvenanceCollector
# ---------------------------------------------------------------------------

class DecisionProvenanceCollector:
    """Builder that accumulates codon and constraint decisions during
    optimization and produces an OptimizationDecisionTrail.

    Typical usage::

        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            constraints=["GCInRange", "NoCrypticSplice"],
        )
        # ... during optimization:
        collector.record_codon_decision(CodonDecision(...))
        collector.record_constraint_decision(ConstraintDecision(...))
        collector.record_iteration({"step": 1, "score": 0.85})
        # ... after optimization:
        trail = collector.finalize(output_dna="ATGGTG...", cai=0.91, gc=0.54)
        print(collector.summary())
        collector.to_json("trail.json")
    """

    def __init__(self) -> None:
        self._protein: str = ""
        self._organism: str = ""
        self._constraints: list[str] = []
        self._gene_name: str | None = None
        self._solver_backend: str = ""
        self._seed: int | None = None
        self._codon_decisions: list[CodonDecision] = []
        self._constraint_decisions: list[ConstraintDecision] = []
        self._iteration_log: list[dict[str, Any]] = []
        self._started: bool = False
        self._finalized: bool = False
        self._trail: OptimizationDecisionTrail | None = None
        logger.debug("DecisionProvenanceCollector initialized")

    def start_optimization(
        self,
        protein: str,
        organism: str,
        constraints: list[str] | None = None,
        *,
        gene_name: str | None = None,
        solver_backend: str = "greedy",
        seed: int | None = None,
    ) -> None:
        """Begin recording a new optimization run.

        Args:
            protein: The input amino acid sequence.
            organism: Target organism identifier.
            constraints: List of constraint/predicate names that will be active.
            gene_name: Optional gene name.
            solver_backend: Solver backend that will be used.
            seed: RNG seed for reproducibility.

        Raises:
            RuntimeError: If an optimization is already in progress.
        """
        if self._started and not self._finalized:
            raise RuntimeError(
                "Optimization already in progress. "
                "Call finalize() before starting a new one."
            )
        self._protein = protein
        self._organism = organism
        self._constraints = list(constraints) if constraints else []
        self._gene_name = gene_name
        self._solver_backend = solver_backend
        self._seed = seed
        self._codon_decisions = []
        self._constraint_decisions = []
        self._iteration_log = []
        self._started = True
        self._finalized = False
        logger.debug(
            "Optimization started: protein_len=%d, organism=%s, "
            "constraints=%d, backend=%s, seed=%s",
            len(protein), organism, len(self._constraints),
            solver_backend, seed,
        )

    def record_codon_decision(self, decision: CodonDecision) -> None:
        """Record a codon-level decision made by the optimizer.

        Args:
            decision: The CodonDecision to store.

        Raises:
            RuntimeError: If optimization has not been started.
            TypeError: If *decision* is not a CodonDecision.
        """
        if not self._started:
            raise RuntimeError(
                "Call start_optimization() before recording decisions."
            )
        if not isinstance(decision, CodonDecision):
            raise TypeError(
                f"Expected CodonDecision, got {type(decision).__name__}"
            )
        self._codon_decisions.append(decision)
        logger.debug(
            "Recorded codon decision: position=%d, aa=%s, chosen=%s, "
            "confidence=%.3f",
            decision.position, decision.amino_acid, decision.chosen_codon,
            decision.confidence,
        )

    def record_constraint_decision(self, decision: ConstraintDecision) -> None:
        """Record a constraint-level decision.

        Args:
            decision: The ConstraintDecision to store.

        Raises:
            RuntimeError: If optimization has not been started.
            TypeError: If *decision* is not a ConstraintDecision.
        """
        if not self._started:
            raise RuntimeError(
                "Call start_optimization() before recording decisions."
            )
        if not isinstance(decision, ConstraintDecision):
            raise TypeError(
                f"Expected ConstraintDecision, got {type(decision).__name__}"
            )
        self._constraint_decisions.append(decision)
        logger.debug(
            "Recorded constraint decision: %s (%s) → %s, CAI impact=%.4f",
            decision.constraint_name, decision.constraint_type,
            decision.action_taken, decision.impact_on_cai,
        )

    def record_iteration(self, iteration_data: dict[str, Any]) -> None:
        """Record what happened at an optimization iteration.

        Args:
            iteration_data: Dict describing the iteration (solver state,
                score changes, etc.).  Typical keys include ``"step"``,
                ``"score"``, ``"violations"``, ``"action"``.

        Raises:
            RuntimeError: If optimization has not been started.
        """
        if not self._started:
            raise RuntimeError(
                "Call start_optimization() before recording iterations."
            )
        self._iteration_log.append(dict(iteration_data))
        logger.debug(
            "Recorded iteration %d: %s",
            len(self._iteration_log),
            iteration_data.get("action", "(no action key)"),
        )

    def finalize(
        self,
        output_dna: str,
        cai: float,
        gc: float,
    ) -> OptimizationDecisionTrail:
        """Finalize the optimization and produce the complete decision trail.

        Args:
            output_dna: The optimized DNA sequence.
            cai: Overall Codon Adaptation Index.
            gc: Overall GC content (0–1).

        Returns:
            An OptimizationDecisionTrail with all recorded decisions.

        Raises:
            RuntimeError: If optimization has not been started.
        """
        if not self._started:
            raise RuntimeError(
                "Call start_optimization() before finalizing."
            )
        timestamp = datetime.now(timezone.utc).isoformat()
        version = _get_biocompiler_version()

        self._trail = OptimizationDecisionTrail(
            gene_name=self._gene_name,
            input_protein=self._protein,
            output_dna=output_dna,
            organism=self._organism,
            solver_backend=self._solver_backend,
            seed=self._seed,
            total_cai=cai,
            total_gc=gc,
            codon_decisions=list(self._codon_decisions),
            constraint_decisions=list(self._constraint_decisions),
            iteration_log=list(self._iteration_log),
            timestamp=timestamp,
            version=version,
        )
        self._finalized = True
        logger.info(
            "Optimization finalized: gene=%s, cai=%.4f, gc=%.3f, "
            "codon_decisions=%d, constraint_decisions=%d, iterations=%d",
            self._gene_name, cai, gc,
            len(self._codon_decisions), len(self._constraint_decisions),
            len(self._iteration_log),
        )
        return self._trail

    def to_json(self, filepath: str) -> None:
        """Serialize the finalized decision trail to a JSON file.

        Args:
            filepath: Path to write the JSON file.

        Raises:
            RuntimeError: If :meth:`finalize` has not been called.
        """
        if not self._finalized or self._trail is None:
            raise RuntimeError(
                "Call finalize() before to_json(). "
                "The trail needs output_dna, cai, and gc values."
            )
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(self._trail.to_json())
        logger.info("Decision trail written to %s", filepath)

    def summary(self) -> str:
        """Produce a human-readable summary of all recorded decisions.

        Returns:
            A formatted multi-line string summarizing the optimization.

        Raises:
            RuntimeError: If optimization has not been started.
        """
        if not self._started:
            return "No optimization has been started."

        lines: list[str] = []
        lines.append("=" * 72)
        lines.append("BioCompiler Decision Provenance Summary")
        lines.append("=" * 72)
        lines.append(f"  Gene:                 {self._gene_name or '(unnamed)'}")
        lines.append(f"  Organism:             {self._organism}")
        lines.append(f"  Solver backend:       {self._solver_backend}")
        lines.append(f"  Seed:                 {self._seed}")
        lines.append(f"  Input protein length: {len(self._protein)} aa")
        lines.append(f"  Constraints active:   "
                     f"{', '.join(self._constraints) or '(none)'}")
        lines.append("")

        # Codon decisions summary
        lines.append(f"  Codon decisions:      {len(self._codon_decisions)}")
        if self._codon_decisions:
            avg_confidence = (
                sum(d.confidence for d in self._codon_decisions)
                / len(self._codon_decisions)
            )
            low_confidence = [
                d for d in self._codon_decisions if d.confidence < 0.5
            ]
            lines.append(f"  Average confidence:   {avg_confidence:.3f}")
            lines.append(f"  Low-confidence (<0.5): {len(low_confidence)}")

            # CAI impact summary for codon decisions
            total_cai_cost_codons = sum(d.cai_impact for d in self._codon_decisions)
            constrained_codons = [
                d for d in self._codon_decisions if d.cai_impact < 0
            ]
            lines.append(f"  Total CAI cost (codon decisions): {total_cai_cost_codons:.4f}")
            lines.append(f"  Constrained codon positions:      {len(constrained_codons)}")

            # Show constraint reasons breakdown
            reason_counts: dict[str, int] = {}
            for d in self._codon_decisions:
                reason_counts[d.constraint_reason] = (
                    reason_counts.get(d.constraint_reason, 0) + 1
                )
            lines.append("  Codon choice drivers:")
            for reason, count in sorted(
                reason_counts.items(), key=lambda x: x[1], reverse=True
            ):
                lines.append(f"    {reason}: {count}")
        lines.append("")

        # Constraint decisions summary
        lines.append(f"  Constraint decisions: {len(self._constraint_decisions)}")
        if self._constraint_decisions:
            total_cai_cost = sum(
                d.impact_on_cai for d in self._constraint_decisions
            )
            satisfied = sum(
                1 for d in self._constraint_decisions
                if d.action_taken == "satisfied"
            )
            relaxed = sum(
                1 for d in self._constraint_decisions
                if d.action_taken == "relaxed"
            )
            conflicted = sum(
                1 for d in self._constraint_decisions
                if d.action_taken == "conflicted"
            )
            overridden = sum(
                1 for d in self._constraint_decisions
                if d.action_taken == "overridden"
            )
            lines.append(f"  Total CAI cost:       {total_cai_cost:.4f}")
            lines.append(f"  Satisfied:            {satisfied}")
            lines.append(f"  Relaxed:              {relaxed}")
            lines.append(f"  Conflicted:           {conflicted}")
            lines.append(f"  Overridden:           {overridden}")

            # Per-constraint breakdown
            lines.append("  Per-constraint breakdown:")
            for d in self._constraint_decisions:
                positions_str = (
                    ", ".join(str(p) for p in d.positions_affected[:5])
                )
                if len(d.positions_affected) > 5:
                    positions_str += f", ... (+{len(d.positions_affected) - 5})"
                lines.append(
                    f"    {d.constraint_name} ({d.constraint_type}): "
                    f"{d.action_taken}, CAI impact={d.impact_on_cai:.4f}, "
                    f"positions=[{positions_str}]"
                )
                lines.append(f"      {d.tradeoff_description}")
        lines.append("")

        # Iteration log summary
        lines.append(f"  Iterations:           {len(self._iteration_log)}")
        if self._iteration_log:
            lines.append("  Iteration log:")
            for i, entry in enumerate(self._iteration_log):
                action = entry.get("action", "(no action)")
                score = entry.get("score", "N/A")
                lines.append(f"    [{i}] action={action}, score={score}")
        lines.append("")

        # Finalized state
        if self._finalized:
            lines.append("  Status:               FINALIZED")
        else:
            lines.append("  Status:               IN PROGRESS (not yet finalized)")

        lines.append("=" * 72)
        return "\n".join(lines)

    def __repr__(self) -> str:
        state = "finalized" if self._finalized else (
            "active" if self._started else "idle"
        )
        return (
            f"DecisionProvenanceCollector("
            f"state={state}, "
            f"codon_decisions={len(self._codon_decisions)}, "
            f"constraint_decisions={len(self._constraint_decisions)}, "
            f"iterations={len(self._iteration_log)})"
        )


# ---------------------------------------------------------------------------
# ProvenanceStore — Persistent storage for optimization provenance records
# ---------------------------------------------------------------------------

_logger = logging.getLogger(__name__)


class ProvenanceStore:
    """Persistent storage for optimization provenance records.

    Saves and loads :class:`OptimizationDecisionTrail` instances as JSON
    files on disk, one per optimization run.  Supports querying by
    protein name, organism, and date range, as well as exporting a
    full audit trail.

    The default store directory is ``~/.biocompiler/provenance/``.
    Each record is stored as a JSON file named ``<record_id>.json``.

    Usage::

        store = ProvenanceStore()
        record_id = store.save(trail)
        loaded = store.load(record_id)
        results = store.query(organism="Homo_sapiens")
        audit_json = store.export_audit_trail(record_id)
    """

    # UUID format regex: strict lowercase hex with hyphens
    _UUID_PATTERN = __import__("re").compile(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
    )

    def __init__(self, store_dir: str | None = None) -> None:
        self._store_dir = Path(store_dir or Path.home() / ".biocompiler" / "provenance").resolve()
        self._store_dir.mkdir(parents=True, exist_ok=True)
        _logger.debug("ProvenanceStore initialized at %s", self._store_dir)

    @property
    def store_dir(self) -> Path:
        """Return the store directory path."""
        return self._store_dir

    def save(self, record: OptimizationDecisionTrail) -> str:
        """Save a provenance record, return its ID.

        Generates a UUID-based record ID and persists the record as a
        JSON file in the store directory.

        Args:
            record: An :class:`OptimizationDecisionTrail` to persist.

        Returns:
            The generated record ID string.
        """
        record_id = str(uuid.uuid4())
        filepath = self._store_dir / f"{record_id}.json"
        data = record.to_dict()
        # Inject the record_id so it is self-describing
        data["_record_id"] = record_id
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, sort_keys=True)
        _logger.info(
            "Saved provenance record %s (%d codon decisions, organism=%s)",
            record_id, len(record.codon_decisions), record.organism,
        )
        return record_id

    @staticmethod
    def _validate_uuid(record_id: str) -> None:
        """Validate that *record_id* is a well-formed UUID string.

        Security: prevents path-traversal attacks by rejecting anything
        that is not a strict ``xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx``
        lowercase-hex UUID.

        Raises:
            ValueError: If *record_id* is not a valid UUID format.
        """
        # Layer 1: Reject path separators
        if "/" in record_id or "\\" in record_id:
            raise ValueError(
                f"Invalid provenance record ID (contains path separators): {record_id!r}"
            )
        # Layer 2: Reject parent-directory traversal
        if ".." in record_id:
            raise ValueError(
                f"Invalid provenance record ID (contains '..'): {record_id!r}"
            )
        # Layer 3: Strict UUID format check (lowercase hex only)
        if not ProvenanceStore._UUID_PATTERN.match(record_id):
            raise ValueError(
                f"Invalid provenance record ID (not a valid UUID): {record_id!r}"
            )
        # Layer 4: uuid.UUID parsing for definitive validation
        try:
            uuid.UUID(record_id)
        except ValueError:
            raise ValueError(
                f"Invalid provenance record ID (UUID parse failure): {record_id!r}"
            )

    def _safe_record_path(self, record_id: str) -> Path:
        """Return a validated, path-contained file path for *record_id*.

        Calls :meth:`_validate_uuid` first, then resolves the resulting
        path and verifies it is still inside ``self._store_dir``.

        Raises:
            ValueError: If *record_id* is invalid or path escapes store.
        """
        self._validate_uuid(record_id)
        filepath = (self._store_dir / f"{record_id}.json").resolve()
        # Defense-in-depth: verify resolved path stays within store dir
        if not str(filepath).startswith(str(self._store_dir)):
            raise ValueError(
                f"Invalid provenance record ID (path escapes store): {record_id!r}"
            )
        return filepath

    def load(self, record_id: str) -> OptimizationDecisionTrail:
        """Load a provenance record by ID.

        Args:
            record_id: The record ID returned by :meth:`save`.
                Must be a valid UUID string (no path traversal).

        Returns:
            The deserialized :class:`OptimizationDecisionTrail`.

        Raises:
            FileNotFoundError: If no record with the given ID exists.
            ValueError: If the record_id is invalid or the stored file
                is corrupted or cannot be deserialized.
        """
        filepath = self._safe_record_path(record_id)
        if not filepath.exists():
            raise FileNotFoundError(
                f"Provenance record not found: {record_id}"
            )
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
            # Remove the injected _record_id before deserialization
            data.pop("_record_id", None)
            return OptimizationDecisionTrail.from_dict(data)
        except (json.JSONDecodeError, ValueError) as exc:
            _logger.error("Failed to load provenance record %s: %s", record_id, exc)
            raise ValueError(
                f"Corrupted provenance record: {record_id}"
            ) from exc

    def query(
        self,
        protein_name: str | None = None,
        organism: str | None = None,
        date_range: tuple[str, str] | None = None,
    ) -> list[OptimizationDecisionTrail]:
        """Query provenance records.

        Scans all stored records and returns those matching all
        specified filters.  Filters are combined with AND logic —
        only records matching *all* non-None criteria are returned.

        Args:
            protein_name: If provided, match records whose
                ``gene_name`` equals this value (case-sensitive).
            organism: If provided, match records whose ``organism``
                equals this value (case-sensitive).
            date_range: If provided, a ``(start_iso, end_iso)`` tuple.
                Only records whose ``timestamp`` falls within the
                inclusive range are returned.  ISO 8601 format
                (e.g. ``"2026-01-01T00:00:00+00:00"``).

        Returns:
            List of matching :class:`OptimizationDecisionTrail` instances.
        """
        results: list[OptimizationDecisionTrail] = []

        for filepath in sorted(self._store_dir.glob("*.json")):
            try:
                with open(filepath, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except (json.JSONDecodeError, OSError) as exc:
                _logger.warning("Skipping corrupted file %s: %s", filepath, exc)
                continue

            # Apply protein_name filter
            if protein_name is not None:
                if data.get("gene_name") != protein_name:
                    continue

            # Apply organism filter
            if organism is not None:
                if data.get("organism") != organism:
                    continue

            # Apply date_range filter
            if date_range is not None:
                start_iso, end_iso = date_range
                ts = data.get("timestamp", "")
                if not ts:
                    continue
                try:
                    ts_dt = datetime.fromisoformat(ts)
                    start_dt = datetime.fromisoformat(start_iso)
                    end_dt = datetime.fromisoformat(end_iso)
                    if not (start_dt <= ts_dt <= end_dt):
                        continue
                except (ValueError, TypeError) as exc:
                    _logger.warning(
                        "Skipping record with unparseable timestamp %s: %s",
                        ts, exc,
                    )
                    continue

            # Remove injected _record_id before deserialization
            data.pop("_record_id", None)
            try:
                results.append(OptimizationDecisionTrail.from_dict(data))
            except ValueError as exc:
                _logger.warning("Skipping invalid record %s: %s", filepath, exc)

        _logger.debug(
            "Query returned %d records (protein_name=%s, organism=%s, "
            "date_range=%s)",
            len(results), protein_name, organism, date_range,
        )
        return results

    def export_audit_trail(self, record_id: str, format: str = "json") -> str:
        """Export full audit trail for a record.

        Loads the record by ID and returns it in the specified format.

        Args:
            record_id: The record ID returned by :meth:`save`.
            format: Output format — ``"json"`` (default) returns an
                indented JSON string.  ``"text"`` returns a compact
                human-readable summary.

        Returns:
            The audit trail in the requested format as a string.

        Raises:
            FileNotFoundError: If no record with the given ID exists.
            ValueError: If the format is not supported.
        """
        record = self.load(record_id)

        if format == "json":
            return record.to_json()
        elif format == "text":
            lines: list[str] = []
            lines.append(f"Audit Trail: {record.gene_name or '(unnamed)'}")
            lines.append(f"Organism:    {record.organism}")
            lines.append(f"Backend:     {record.solver_backend}")
            lines.append(f"Timestamp:   {record.timestamp}")
            lines.append(f"CAI:         {record.total_cai:.4f}")
            lines.append(f"GC:          {record.total_gc:.3f}")
            lines.append(f"Codon decisions: {len(record.codon_decisions)}")
            for cd in record.codon_decisions:
                lines.append(
                    f"  [{cd.position}] {cd.amino_acid} -> {cd.chosen_codon} "
                    f"(reason={cd.constraint_reason}, confidence={cd.confidence:.3f})"
                )
            lines.append(f"Constraint decisions: {len(record.constraint_decisions)}")
            for cd in record.constraint_decisions:
                lines.append(
                    f"  {cd.constraint_name}: {cd.action_taken} "
                    f"(CAI impact={cd.impact_on_cai:.4f})"
                )
            return "\n".join(lines)
        else:
            raise ValueError(
                f"Unsupported audit trail format: {format!r}. "
                "Use 'json' or 'text'."
            )
