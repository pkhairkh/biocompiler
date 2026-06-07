"""
BioCompiler Provenance — Decision Audit Trail & Output Provenance
==================================================================

Provides a structured record of every decision made during gene optimization,
enabling full reproducibility, debugging, and regulatory traceability.

Every gene optimization involves a sequence of decisions — which codon to
choose at each position, which mutations to apply, when to relax constraints.
This module captures all of those decisions as ``DecisionRecord`` instances
and organizes them via the ``ProvenanceTracker``.  The tracker supports
position-based querying, full audit trail retrieval, and serialization to
dict/JSON for persistence.

At the run level, ``OptimizationProvenance`` ties together the input protein,
solver configuration, all decisions, and the final output into a single
snapshot.  ``OptimizationRecord`` provides a lightweight summary with
reproducibility metadata (seed, biocompiler version, timestamps).  The
``generate_provenance_report`` function produces a human-readable text report
from a list of optimization records.

Usage::

    from biocompiler.provenance import (
        DecisionRecord, ProvenanceTracker, OptimizationRecord,
        generate_provenance_report,
    )
    from datetime import datetime, timezone

    # Create a tracker for an optimization run
    tracker = ProvenanceTracker(seed=42)

    # Record decisions during optimization
    tracker.record_decision(DecisionRecord(
        timestamp=datetime.now(timezone.utc).isoformat(),
        decision_type="codon_selected",
        position=12,
        chosen_value="GTC",
        alternatives_considered=["GTG", "GTA", "GTT"],
        rationale="Highest CAI codon for Valine in Homo_sapiens",
        constraint_context={"cai": 0.92, "gc": 0.54},
    ))

    # Query decisions by position
    decisions_at_12 = tracker.get_decisions_for_position(12)

    # Get the full audit trail
    trail = tracker.get_full_audit_trail()

    # Serialize for persistence
    json_str = tracker.to_json()

    # Generate a report from optimization records
    report = generate_provenance_report(tracker.get_optimization_records())
    print(report)

Key components:
- DecisionRecord: immutable record of a single solver decision
- ProvenanceTracker: accumulates decisions per position and OptimizationRecords
  per run, supports querying and serialization (dict / JSON)
- OptimizationProvenance: end-to-end provenance snapshot tying input, config,
  all decisions, and final output together
- OptimizationRecord: summary record of a single optimization run with
  reproducibility metadata (seed, version, timestamps)
- generate_provenance_report: human-readable report from a list of
  OptimizationRecords

Design goals:
1. Every codon choice, mutation, and constraint relaxation is recorded with
   the alternatives considered and the rationale — no silent decisions.
2. The ``seed`` on ProvenanceTracker ensures reproducibility: the same seed +
   same input + same solver must yield an identical decision trail.
3. Serialization (to_dict / to_json) enables persistence, comparison across
   runs, and integration with certificate provenance metadata.
4. OptimizationRecord provides a lightweight, serializable summary that
   captures all information needed to reproduce or audit an optimization run,
   including the biocompiler version and seeded randomness.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "DecisionRecord",
    "ProvenanceTracker",
    "OptimizationProvenance",
    "OptimizationRecord",
    "ProvenanceStore",
    "generate_provenance_report",
    # Decision category constants
    "DECISION_CATEGORY_CAI",
    "DECISION_CATEGORY_GT_AVOIDANCE",
    "DECISION_CATEGORY_GC_CONTENT",
    "DECISION_CATEGORY_RESTRICTION_SITE",
    "DECISION_CATEGORY_SPLICE_PREVENTION",
    "DECISION_CATEGORY_MUTATION",
    "DECISION_CATEGORY_CONSTRAINT_RELAXATION",
    "DECISION_CATEGORY_OTHER",
    "ALL_DECISION_CATEGORIES",
]

# UUID format regex — enforces standard lowercase hex UUID format to prevent
# path traversal attacks when IDs are used as file path components.
_UUID_PATTERN = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
)


# ---------------------------------------------------------------------------
# DecisionRecord
# ---------------------------------------------------------------------------

# Standard decision categories
DECISION_CATEGORY_CAI = "cai_optimization"
DECISION_CATEGORY_GT_AVOIDANCE = "gt_avoidance"
DECISION_CATEGORY_GC_CONTENT = "gc_content"
DECISION_CATEGORY_RESTRICTION_SITE = "restriction_site"
DECISION_CATEGORY_SPLICE_PREVENTION = "splice_prevention"
DECISION_CATEGORY_MUTATION = "mutation"
DECISION_CATEGORY_CONSTRAINT_RELAXATION = "constraint_relaxation"
DECISION_CATEGORY_OTHER = "other"

ALL_DECISION_CATEGORIES: list[str] = [
    DECISION_CATEGORY_CAI,
    DECISION_CATEGORY_GT_AVOIDANCE,
    DECISION_CATEGORY_GC_CONTENT,
    DECISION_CATEGORY_RESTRICTION_SITE,
    DECISION_CATEGORY_SPLICE_PREVENTION,
    DECISION_CATEGORY_MUTATION,
    DECISION_CATEGORY_CONSTRAINT_RELAXATION,
    DECISION_CATEGORY_OTHER,
]


@dataclass
class DecisionRecord:
    """Immutable record of a single decision made by the optimizer.

    Attributes:
        timestamp: ISO 8601 timestamp of when the decision was made.
        decision_type: Category of the decision, e.g. ``"codon_selected"``,
            ``"mutation_applied"``, ``"constraint_relaxed"``.
        position: 0-based nucleotide position in the output sequence to which
            this decision applies.
        chosen_value: The value that was selected (e.g. a codon like ``"GTC"``
            or a mutation like ``"V42I"``).
        alternatives_considered: Ordered list of values that were evaluated but
            not chosen. Ordering reflects solver preference (best first).
        rationale: Human-readable explanation of *why* ``chosen_value`` was
            preferred over the alternatives.
        constraint_context: Snapshot of active constraints / scores at the time
            of the decision.  Free-form dict — typical keys include
            ``"cai"``, ``"gc"``, ``"max_donor_score"``, ``"active_predicates"``.
        cai_impact: CAI delta caused by this decision.  Negative values
            indicate CAI was sacrificed (e.g. to satisfy a constraint);
            positive values indicate CAI was gained; zero means no CAI
            impact or not applicable.
        codon_before: The codon at this position *before* the decision was
            applied.  Empty string if not applicable (e.g. de novo).
        codon_after: The codon at this position *after* the decision was
            applied (typically same as ``chosen_value`` for codon decisions).
            Empty string if not applicable.
        decision_category: High-level category for this decision, used to
            separate GT avoidance decisions from CAI decisions for eukaryotes.
            One of the ``DECISION_CATEGORY_*`` constants.  Defaults to
            ``DECISION_CATEGORY_OTHER`` for backward compatibility.
        gt_avoidance_impact: For eukaryotic organisms, records whether this
            decision was influenced by GT dinucleotide avoidance (splice donor
            prevention).  Non-zero when a codon was chosen specifically to
            avoid creating a GT at a splice-relevant position; zero otherwise.
            This is separate from ``cai_impact`` — a decision may sacrifice
            CAI for GT avoidance, and both impacts should be recorded.
    """

    timestamp: str
    decision_type: str
    position: int
    chosen_value: str
    alternatives_considered: list[str]
    rationale: str
    constraint_context: dict[str, Any]
    # CAI-aware provenance fields
    cai_impact: float = 0.0
    codon_before: str = ""
    codon_after: str = ""
    # Category and GT-avoidance fields (eukaryote-aware provenance)
    decision_category: str = DECISION_CATEGORY_OTHER
    gt_avoidance_impact: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        """Serialize this record to a JSON-compatible dict."""
        return {
            "timestamp": self.timestamp,
            "decision_type": self.decision_type,
            "position": self.position,
            "chosen_value": self.chosen_value,
            "alternatives_considered": list(self.alternatives_considered),
            "rationale": self.rationale,
            "constraint_context": dict(self.constraint_context),
            "cai_impact": self.cai_impact,
            "codon_before": self.codon_before,
            "codon_after": self.codon_after,
            "decision_category": self.decision_category,
            "gt_avoidance_impact": self.gt_avoidance_impact,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DecisionRecord":
        """Deserialize a DecisionRecord from a plain dict.

        Args:
            data: Dict with the same keys as the dataclass fields.

        Returns:
            A DecisionRecord instance.

        Raises:
            ValueError: If required keys are missing.
        """
        required = {
            "timestamp", "decision_type", "position",
            "chosen_value", "alternatives_considered",
            "rationale", "constraint_context",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(f"Cannot deserialize DecisionRecord: missing keys {missing}")
        return cls(
            timestamp=data["timestamp"],
            decision_type=data["decision_type"],
            position=int(data["position"]),
            chosen_value=data["chosen_value"],
            alternatives_considered=list(data["alternatives_considered"]),
            rationale=data["rationale"],
            constraint_context=dict(data["constraint_context"]),
            cai_impact=float(data.get("cai_impact", 0.0)),
            codon_before=str(data.get("codon_before", "")),
            codon_after=str(data.get("codon_after", "")),
            decision_category=str(data.get("decision_category", DECISION_CATEGORY_OTHER)),
            gt_avoidance_impact=float(data.get("gt_avoidance_impact", 0.0)),
        )


# ---------------------------------------------------------------------------
# ProvenanceTracker
# ---------------------------------------------------------------------------

class ProvenanceTracker:
    """Accumulates DecisionRecords and provides position-based querying.

    The tracker is typically instantiated once per optimization run and passed
    through the solver pipeline so that every decision point can call
    :meth:`record_decision`.

    Attributes:
        seed: Integer seed used for the solver's RNG.  Recording this inside
            the tracker ensures the seed is always co-serialized with the
            decisions it influenced.

    Example::

        tracker = ProvenanceTracker(seed=42)
        tracker.record_decision(DecisionRecord(
            timestamp=datetime.now(timezone.utc).isoformat(),
            decision_type="codon_selected",
            position=12,
            chosen_value="GTC",
            alternatives_considered=["GTG", "GTA", "GTT"],
            rationale="Highest CAI codon for Valine in Homo_sapiens",
            constraint_context={"cai": 0.92, "gc": 0.54},
        ))
        decisions_at_12 = tracker.get_decisions_for_position(12)
    """

    def __init__(self, seed: int = 0) -> None:
        self.seed: int = seed
        self._decisions: list[DecisionRecord] = []
        # Position index for O(1) lookup by position
        self._position_index: dict[int, list[DecisionRecord]] = {}
        # Optimization records for run-level provenance
        self._optimization_records: list[OptimizationRecord] = []
        logger.debug("ProvenanceTracker initialized with seed=%d", seed)

    # -- Mutation -----------------------------------------------------------

    def record_decision(self, decision: DecisionRecord) -> None:
        """Append a decision record to the audit trail.

        Args:
            decision: The DecisionRecord to store.

        Raises:
            TypeError: If *decision* is not a DecisionRecord.
        """
        if not isinstance(decision, DecisionRecord):
            raise TypeError(
                f"Expected DecisionRecord, got {type(decision).__name__}"
            )
        self._decisions.append(decision)
        pos = decision.position
        if pos not in self._position_index:
            self._position_index[pos] = []
        self._position_index[pos].append(decision)
        logger.debug(
            "Recorded decision: type=%s position=%d chosen=%s",
            decision.decision_type,
            decision.position,
            decision.chosen_value,
        )

    # -- Query --------------------------------------------------------------

    def get_decisions_for_position(self, position: int) -> list[DecisionRecord]:
        """Return all decisions recorded for the given position.

        Args:
            position: 0-based nucleotide position.

        Returns:
            List of DecisionRecords at *position*, in chronological order.
            Returns an empty list if no decisions exist for the position.
        """
        return list(self._position_index.get(position, []))

    def get_full_audit_trail(self) -> list[DecisionRecord]:
        """Return the complete ordered list of all recorded decisions.

        Returns:
            Chronologically ordered list of all DecisionRecords.
        """
        return list(self._decisions)

    def get_decisions_by_category(self, category: str) -> list[DecisionRecord]:
        """Return all decisions matching the given decision category.

        This is the primary way to separate GT avoidance decisions from CAI
        decisions for eukaryotic organisms.

        Args:
            category: One of the ``DECISION_CATEGORY_*`` constants, e.g.
                ``DECISION_CATEGORY_GT_AVOIDANCE`` or
                ``DECISION_CATEGORY_CAI``.

        Returns:
            List of DecisionRecords whose ``decision_category`` matches
            *category*, in chronological order.
        """
        return [d for d in self._decisions if d.decision_category == category]

    def get_gt_avoidance_decisions(self) -> list[DecisionRecord]:
        """Return all GT avoidance decisions (convenience method for eukaryotes).

        For eukaryotic organisms, GT dinucleotide avoidance is a distinct
        concern from CAI optimization.  This method returns only the decisions
        that were driven by GT avoidance, allowing separate analysis of CAI
        vs. splice-donor-prevention tradeoffs.

        Returns:
            List of DecisionRecords with
            ``decision_category == DECISION_CATEGORY_GT_AVOIDANCE``.
        """
        return self.get_decisions_by_category(DECISION_CATEGORY_GT_AVOIDANCE)

    def get_cai_decisions(self) -> list[DecisionRecord]:
        """Return all CAI optimization decisions (convenience method).

        Returns:
            List of DecisionRecords with
            ``decision_category == DECISION_CATEGORY_CAI``.
        """
        return self.get_decisions_by_category(DECISION_CATEGORY_CAI)

    # -- Optimization records ------------------------------------------------

    def add_optimization_record(self, record: OptimizationRecord) -> None:
        """Append an OptimizationRecord summarizing a completed optimization run.

        Args:
            record: The OptimizationRecord to store.

        Raises:
            TypeError: If *record* is not an OptimizationRecord.
        """
        if not isinstance(record, OptimizationRecord):
            raise TypeError(
                f"Expected OptimizationRecord, got {type(record).__name__}"
            )
        self._optimization_records.append(record)
        logger.debug(
            "Recorded optimization: organism=%s backend=%s seed=%s",
            record.organism,
            record.solver_backend,
            record.seed_used,
        )

    def get_optimization_records(self) -> list[OptimizationRecord]:
        """Return all stored OptimizationRecords.

        Returns:
            List of OptimizationRecords in chronological order.
        """
        return list(self._optimization_records)

    # -- Serialization ------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Serialize the tracker to a JSON-compatible dict.

        The returned dict includes the seed, all decision records, and all
        optimization records, making it self-contained for persistence.

        Returns:
            Dict with keys ``"seed"``, ``"decision_count"``, ``"decisions"``,
            ``"optimization_records"``.
        """
        return {
            "seed": self.seed,
            "decision_count": len(self._decisions),
            "decisions": [d.to_dict() for d in self._decisions],
            "optimization_records": [r.to_dict() for r in self._optimization_records],
        }

    def to_json(self) -> str:
        """Serialize the tracker to a JSON string.

        Returns:
            Indented JSON string representation of the tracker.
        """
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ProvenanceTracker":
        """Reconstruct a ProvenanceTracker from a serialized dict.

        Args:
            data: Dict produced by :meth:`to_dict`.

        Returns:
            A ProvenanceTracker with all decisions restored.

        Raises:
            ValueError: If required keys are missing or records are invalid.
        """
        if "seed" not in data:
            raise ValueError("Cannot deserialize ProvenanceTracker: missing key 'seed'")
        tracker = cls(seed=data["seed"])
        for rec_data in data.get("decisions", []):
            tracker.record_decision(DecisionRecord.from_dict(rec_data))
        for opt_data in data.get("optimization_records", []):
            tracker.add_optimization_record(OptimizationRecord.from_dict(opt_data))
        return tracker

    @classmethod
    def from_json(cls, json_str: str) -> "ProvenanceTracker":
        """Reconstruct a ProvenanceTracker from a JSON string.

        Args:
            json_str: JSON string produced by :meth:`to_json`.

        Returns:
            A ProvenanceTracker with all decisions restored.
        """
        return cls.from_dict(json.loads(json_str))

    # -- Introspection ------------------------------------------------------

    def __len__(self) -> int:
        """Return the number of recorded decisions."""
        return len(self._decisions)

    def __repr__(self) -> str:
        return (
            f"ProvenanceTracker(seed={self.seed}, "
            f"decisions={len(self._decisions)}, "
            f"optimization_records={len(self._optimization_records)})"
        )


# ---------------------------------------------------------------------------
# OptimizationProvenance
# ---------------------------------------------------------------------------

@dataclass
class OptimizationProvenance:
    """End-to-end provenance snapshot for a single optimization run.

    Ties together the input protein, solver configuration, every decision
    made, and the final output sequence.  Designed for persistence alongside
    certificates so that any third party can audit *why* a particular output
    was produced.

    Attributes:
        input_protein: The original amino acid sequence submitted for
            optimization.
        organism: Target organism identifier (e.g. ``"Homo_sapiens"``).
        solver_backend: Name of the solver backend used
            (e.g. ``"greedy"``, ``"z3"``, ``"ortools"``).
        config_snapshot: Deep copy of the optimizer configuration at the start
            of the run.  Includes GC bounds, restriction enzymes, thresholds,
            etc.
        decisions: Ordered list of every DecisionRecord produced during the
            run.
        final_sequence: The optimized DNA sequence that was output.
        solve_time_seconds: Wall-clock time of the full solve, in seconds.
        constraints_active: List of constraint/predicate names that were
            active during the run (e.g. ``["NoCrypticSplice", "GCInRange"]``).
    """

    input_protein: str
    organism: str
    solver_backend: str
    config_snapshot: dict[str, Any]
    decisions: list[DecisionRecord]
    final_sequence: str
    solve_time_seconds: float
    constraints_active: list[str]

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict.

        Returns:
            Dict representation of the full optimization provenance.
        """
        return {
            "input_protein": self.input_protein,
            "organism": self.organism,
            "solver_backend": self.solver_backend,
            "config_snapshot": dict(self.config_snapshot),
            "decisions": [d.to_dict() for d in self.decisions],
            "final_sequence": self.final_sequence,
            "solve_time_seconds": self.solve_time_seconds,
            "constraints_active": list(self.constraints_active),
        }

    def to_json(self) -> str:
        """Serialize to an indented JSON string.

        Returns:
            JSON string representation.
        """
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "OptimizationProvenance":
        """Deserialize from a plain dict.

        Args:
            data: Dict with the same keys as the dataclass fields.

        Returns:
            An OptimizationProvenance instance.

        Raises:
            ValueError: If required keys are missing.
        """
        required = {
            "input_protein", "organism", "solver_backend",
            "config_snapshot", "decisions", "final_sequence",
            "solve_time_seconds", "constraints_active",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize OptimizationProvenance: missing keys {missing}"
            )
        decisions = [
            DecisionRecord.from_dict(d) for d in data["decisions"]
        ]
        return cls(
            input_protein=data["input_protein"],
            organism=data["organism"],
            solver_backend=data["solver_backend"],
            config_snapshot=dict(data["config_snapshot"]),
            decisions=decisions,
            final_sequence=data["final_sequence"],
            solve_time_seconds=float(data["solve_time_seconds"]),
            constraints_active=list(data["constraints_active"]),
        )

    @classmethod
    def from_json(cls, json_str: str) -> "OptimizationProvenance":
        """Deserialize from a JSON string.

        Args:
            json_str: JSON string produced by :meth:`to_json`.

        Returns:
            An OptimizationProvenance instance.
        """
        return cls.from_dict(json.loads(json_str))

    def __repr__(self) -> str:
        return (
            f"OptimizationProvenance("
            f"organism={self.organism!r}, "
            f"backend={self.solver_backend!r}, "
            f"decisions={len(self.decisions)}, "
            f"solve_time={self.solve_time_seconds:.3f}s)"
        )


# ---------------------------------------------------------------------------
# OptimizationRecord — lightweight run summary for reproducibility
# ---------------------------------------------------------------------------

def _get_biocompiler_version() -> str:
    """Return the current biocompiler version string."""
    try:
        from .. import __version__
        return __version__
    except ImportError:
        return "unknown"


def _coerce_stability_score(value: Any) -> float | None:
    """Coerce an MRNAStabilityScore object (or plain float) to a plain float.

    The optimizer may populate ``mrna_stability_score`` with an
    ``MRNAStabilityScore`` dataclass instead of a plain float.  This helper
    normalises it for JSON serialization and report formatting.
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    # MRNAStabilityScore or similar objects with an overall_score attribute
    if hasattr(value, "overall_score"):
        return float(value.overall_score)
    # Fallback: try to cast
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


@dataclass
class OptimizationRecord:
    """Summary record of a single optimization run for reproducibility auditing.

    Unlike :class:`OptimizationProvenance` (which captures every fine-grained
    decision), this record is a lightweight summary that captures *all*
    information needed to reproduce or audit an optimization run — including
    the biocompiler version, the seed used for randomness, and the list of
    mutations applied.

    Attributes:
        input_sequence: The original input (protein or DNA) submitted.
        output_sequence: The optimized DNA sequence produced.
        organism: Target organism identifier (e.g. ``"Homo_sapiens"``).
        constraints_applied: List of constraint/predicate names active during
            the run (e.g. ``["NoCrypticSplice", "GCInRange"]``).
        mutations_made: List of amino acid substitutions applied during
            optimization (e.g. ``["V42I", "V98I"]``).
        solver_backend: Solver backend used (e.g. ``"greedy"``, ``"z3"``).
        solve_time: Wall-clock solve time in seconds.
        seed_used: The RNG seed used for this run. ``None`` if no seed was
            provided (non-reproducible).
        timestamp: ISO 8601 timestamp of when the optimization completed.
        biocompiler_version: Version of biocompiler that produced this record.
    """

    input_sequence: str
    output_sequence: str
    organism: str
    constraints_applied: list[str]
    mutations_made: list[str]
    solver_backend: str
    solve_time: float
    seed_used: int | None
    timestamp: str
    biocompiler_version: str
    # mRNA stability provenance (optional — populated when optimize_mrna_stability=True)
    mrna_stability_score: float | None = None
    destabilizing_motifs_removed: int = 0
    stability_improvement: float | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize this record to a JSON-compatible dict."""
        return {
            "input_sequence": self.input_sequence,
            "output_sequence": self.output_sequence,
            "organism": self.organism,
            "constraints_applied": list(self.constraints_applied),
            "mutations_made": list(self.mutations_made),
            "solver_backend": self.solver_backend,
            "solve_time": self.solve_time,
            "seed_used": self.seed_used,
            "timestamp": self.timestamp,
            "biocompiler_version": self.biocompiler_version,
            "mrna_stability_score": _coerce_stability_score(self.mrna_stability_score),
            "destabilizing_motifs_removed": self.destabilizing_motifs_removed,
            "stability_improvement": _coerce_stability_score(self.stability_improvement),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "OptimizationRecord":
        """Deserialize an OptimizationRecord from a plain dict.

        Args:
            data: Dict with the same keys as the dataclass fields.

        Returns:
            An OptimizationRecord instance.

        Raises:
            ValueError: If required keys are missing.
        """
        required = {
            "input_sequence", "output_sequence", "organism",
            "constraints_applied", "mutations_made", "solver_backend",
            "solve_time", "seed_used", "timestamp", "biocompiler_version",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize OptimizationRecord: missing keys {missing}"
            )
        return cls(
            input_sequence=data["input_sequence"],
            output_sequence=data["output_sequence"],
            organism=data["organism"],
            constraints_applied=list(data["constraints_applied"]),
            mutations_made=list(data["mutations_made"]),
            solver_backend=data["solver_backend"],
            solve_time=float(data["solve_time"]),
            seed_used=data["seed_used"],
            timestamp=data["timestamp"],
            biocompiler_version=data["biocompiler_version"],
            mrna_stability_score=data.get("mrna_stability_score"),
            destabilizing_motifs_removed=data.get("destabilizing_motifs_removed", 0),
            stability_improvement=data.get("stability_improvement"),
        )

    def to_json(self) -> str:
        """Serialize to an indented JSON string."""
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_json(cls, json_str: str) -> "OptimizationRecord":
        """Deserialize from a JSON string."""
        return cls.from_dict(json.loads(json_str))

    def __repr__(self) -> str:
        return (
            f"OptimizationRecord("
            f"organism={self.organism!r}, "
            f"backend={self.solver_backend!r}, "
            f"mutations={len(self.mutations_made)}, "
            f"seed={self.seed_used}, "
            f"solve_time={self.solve_time:.3f}s)"
        )


# ---------------------------------------------------------------------------
# generate_provenance_report
# ---------------------------------------------------------------------------

def generate_provenance_report(records: list[OptimizationRecord]) -> str:
    """Produce a human-readable provenance report from optimization records.

    The report includes a summary table of all runs, details of mutations,
    and reproducibility information (seed, version).

    Args:
        records: List of OptimizationRecord instances to summarize.

    Returns:
        A formatted multi-line string report.
    """
    if not records:
        return "No optimization records to report.\n"

    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("BioCompiler Provenance Report")
    lines.append("=" * 72)
    lines.append(f"Total optimization runs: {len(records)}")
    lines.append("")

    for idx, rec in enumerate(records, 1):
        lines.append(f"--- Run {idx} ---")
        lines.append(f"  Timestamp:            {rec.timestamp}")
        lines.append(f"  BioCompiler version:  {rec.biocompiler_version}")
        lines.append(f"  Organism:             {rec.organism}")
        lines.append(f"  Solver backend:       {rec.solver_backend}")
        lines.append(f"  Seed:                 {rec.seed_used}")
        lines.append(f"  Solve time:           {rec.solve_time:.3f}s")
        lines.append(f"  Input length:         {len(rec.input_sequence)} chars")
        lines.append(f"  Output length:        {len(rec.output_sequence)} chars")
        lines.append(f"  Constraints applied:  {', '.join(rec.constraints_applied) or '(none)'}")
        lines.append(f"  Mutations made:       {', '.join(rec.mutations_made) or '(none)'}")
        if rec.mrna_stability_score is not None:
            score_val = _coerce_stability_score(rec.mrna_stability_score)
            if score_val is not None:
                lines.append(f"  mRNA stability score: {score_val:.4f}")
            lines.append(f"  Motifs removed:       {rec.destabilizing_motifs_removed}")
            if rec.stability_improvement is not None:
                imp_val = _coerce_stability_score(rec.stability_improvement)
                if imp_val is not None:
                    lines.append(f"  Stability improvement: {imp_val:+.4f}")
        # Report GT avoidance vs CAI decision breakdown for eukaryotes
        gt_predicates = {"NoGTDinucleotide", "NoCrypticSplice", "SpliceCorrect"}
        has_gt_constraint = bool(gt_predicates & set(rec.constraints_applied))
        if has_gt_constraint:
            lines.append("  GT avoidance:         active (eukaryote splice-donor prevention)")
        lines.append("")

    # Summary statistics
    total_mutations = sum(len(r.mutations_made) for r in records)
    total_time = sum(r.solve_time for r in records)
    backends_used = sorted(set(r.solver_backend for r in records))
    organisms_targeted = sorted(set(r.organism for r in records))
    seeded_runs = sum(1 for r in records if r.seed_used is not None)

    lines.append("-" * 72)
    lines.append("Summary")
    lines.append("-" * 72)
    lines.append(f"  Total runs:           {len(records)}")
    lines.append(f"  Total solve time:     {total_time:.3f}s")
    lines.append(f"  Total mutations:      {total_mutations}")
    lines.append(f"  Backends used:        {', '.join(backends_used)}")
    lines.append(f"  Organisms targeted:   {', '.join(organisms_targeted)}")
    lines.append(f"  Seeded (reproducible): {seeded_runs}/{len(records)}")
    if seeded_runs < len(records):
        lines.append("  WARNING: Some runs lack a seed and are NOT reproducible.")
    lines.append("=" * 72)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# ProvenanceStore — file-backed provenance storage with UUID validation
# ---------------------------------------------------------------------------

class ProvenanceStore:
    """File-backed store for :class:`ProvenanceTracker` objects.

    Each provenance record is persisted as a JSON file named ``<uuid>.json``
    under a configurable base directory.  The UUID format is strictly
    validated before any filesystem access to prevent **path traversal**
    attacks (e.g. ``"../../../etc/passwd"``).

    Example::

        store = ProvenanceStore(base_dir=Path("/data/provenance"))

        # Save a tracker
        tracker = ProvenanceTracker(seed=42)
        record_id = store.save(tracker)

        # Load it back
        restored = store.load(record_id)

    Security note:
        The ``load()`` method validates that the *record_id* matches the
        standard UUID format (lowercase hexadecimal,
        ``xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx``) before constructing any
        file path.  This prevents directory traversal attacks where an
        attacker could supply ``"../../../etc/passwd"`` to read arbitrary
        files on the system.
    """

    def __init__(self, base_dir: Path | str | None = None) -> None:
        if base_dir is None:
            base_dir = Path(".") / "provenance_store"
        self._base_dir = Path(base_dir)
        self._base_dir.mkdir(parents=True, exist_ok=True)
        logger.debug("ProvenanceStore initialized at %s", self._base_dir)

    @staticmethod
    def _validate_uuid(record_id: str) -> None:
        """Validate that *record_id* is a properly formatted UUID.

        Args:
            record_id: The ID string to validate.

        Raises:
            ValueError: If *record_id* does not match the standard UUID format
                (lowercase hexadecimal, ``xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx``).
        """
        if not _UUID_PATTERN.match(record_id):
            raise ValueError(
                f"Invalid provenance record ID: {record_id!r}. "
                f"ID must be a valid UUID in lowercase hexadecimal format "
                f"(e.g. 'a1b2c3d4-e5f6-7890-abcd-ef1234567890'). "
                f"This check prevents path traversal attacks."
            )

    def save(self, tracker: ProvenanceTracker) -> str:
        """Persist a :class:`ProvenanceTracker` to a JSON file.

        Args:
            tracker: The provenance tracker to persist.

        Returns:
            The UUID record ID under which the tracker was stored.
        """
        import uuid
        record_id = str(uuid.uuid4())
        file_path = self._base_dir / f"{record_id}.json"
        file_path.write_text(tracker.to_json(), encoding="utf-8")
        logger.debug("Saved provenance tracker: id=%s path=%s", record_id, file_path)
        return record_id

    def load(self, record_id: str) -> ProvenanceTracker:
        """Load a :class:`ProvenanceTracker` from a JSON file by its UUID.

        The *record_id* **must** be a valid UUID string in lowercase
        hexadecimal format.  This validation prevents path traversal attacks
        where a malicious caller could supply ``"../../../etc/passwd"`` to
        read arbitrary files.

        Args:
            record_id: UUID string identifying the provenance record.

        Returns:
            The deserialized :class:`ProvenanceTracker`.

        Raises:
            ValueError: If *record_id* is not a valid UUID format.
            FileNotFoundError: If no provenance file exists for *record_id*.
        """
        self._validate_uuid(record_id)
        file_path = self._base_dir / f"{record_id}.json"
        if not file_path.exists():
            raise FileNotFoundError(
                f"Provenance record not found: {record_id} "
                f"(expected at {file_path})"
            )
        json_str = file_path.read_text(encoding="utf-8")
        tracker = ProvenanceTracker.from_json(json_str)
        logger.debug("Loaded provenance tracker: id=%s", record_id)
        return tracker

    def delete(self, record_id: str) -> None:
        """Delete a persisted provenance record by its UUID.

        Args:
            record_id: UUID string identifying the provenance record.

        Raises:
            ValueError: If *record_id* is not a valid UUID format.
            FileNotFoundError: If no provenance file exists for *record_id*.
        """
        self._validate_uuid(record_id)
        file_path = self._base_dir / f"{record_id}.json"
        if not file_path.exists():
            raise FileNotFoundError(
                f"Provenance record not found: {record_id} "
                f"(expected at {file_path})"
            )
        file_path.unlink()
        logger.debug("Deleted provenance record: id=%s", record_id)
