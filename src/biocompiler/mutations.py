"""BioCompiler Mutations — Unified mutation types and utilities.

v9.0.0 — Unified with MutationResult from engine_base.py

Provides:
  - Mutation: backward-compatible subclass of MutationResult
  - MutationSuggestion: enriched mutation with ranking and confidence
  - rank_mutations(): sort mutations by impact (works with both types)
  - combine_mutations(): merge mutation lists from multiple engines

All engines that produce mutation suggestions should use MutationResult
(from engine_base) as the canonical type.  This module provides the
legacy Mutation class (a subclass) for backward compatibility, plus
utility functions that accept both Mutation and MutationResult.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Union

from .engine_base import MutationResult


# ---------------------------------------------------------------------------
# Mutation — backward-compatible subclass of MutationResult
# ---------------------------------------------------------------------------

@dataclass
class Mutation(MutationResult):
    """Legacy mutation type, now a subclass of MutationResult.

    Preserved for backward compatibility.  All new code should prefer
    MutationResult directly.  Mutation IS-A MutationResult, so any
    function accepting MutationResult also accepts Mutation.

    Since v9.0.0, MutationResult has its own `confidence` field, so
    Mutation no longer adds a separate one — it simply inherits it.

    Conversion:
      - to_mutation_result() → plain MutationResult copy
      - Mutation.from_mutation_result(mr) → Mutation copy from any MutationResult
    """

    # source_engine is a backward-compat alias for the inherited `engine` field

    @property
    def source_engine(self) -> str:
        """Alias for the inherited `engine` field (backward compat)."""
        return self.engine

    @source_engine.setter
    def source_engine(self, value: str) -> None:
        self.engine = value

    # -- old field name aliases -------------------------------------------------

    @property
    def original_aa(self) -> str:
        """Alias for `original` (backward compat)."""
        return self.original

    @original_aa.setter
    def original_aa(self, value: str) -> None:
        self.original = value

    @property
    def mutant_aa(self) -> str:
        """Alias for `mutant` (backward compat)."""
        return self.mutant

    @mutant_aa.setter
    def mutant_aa(self, value: str) -> None:
        self.mutant = value

    # -- conversion methods ----------------------------------------------------

    def to_mutation_result(self) -> MutationResult:
        """Convert to a plain MutationResult (strips Mutation-specific fields)."""
        return MutationResult(
            position=self.position,
            original=self.original,
            mutant=self.mutant,
            delta_score=self.delta_score,
            score_type=self.score_type,
            engine=self.engine,
            recommendation=self.recommendation,
            description=self.description,
            confidence=self.confidence,
            details=self.details,
        )

    @classmethod
    def from_mutation_result(cls, mr: MutationResult, **overrides) -> Mutation:
        """Create a Mutation from any MutationResult instance.

        Extra keyword arguments override fields on the new Mutation.
        """
        return cls(
            position=mr.position,
            original=mr.original,
            mutant=mr.mutant,
            delta_score=mr.delta_score,
            score_type=mr.score_type,
            engine=mr.engine,
            recommendation=mr.recommendation,
            description=mr.description,
            confidence=mr.confidence,
            **overrides,
        )


# ---------------------------------------------------------------------------
# MutationSuggestion — enriched mutation with ranking and metadata
# ---------------------------------------------------------------------------

@dataclass
class MutationSuggestion:
    """A mutation suggestion with ranking metadata.

    Wraps a MutationResult (or Mutation) with additional context:
      - rank: ordinal rank within a ranked list (1 = best)
      - confidence: confidence score (0.0–1.0)
      - rationale: human-readable explanation of why this mutation is suggested
      - supporting_engines: list of engine names that independently suggest this mutation

    This is NOT a subclass of MutationResult; it composes one instead,
    because a suggestion is metadata *about* a mutation, not a mutation itself.
    """

    mutation: MutationResult
    rank: int = 0
    confidence: float = 1.0
    rationale: str = ""
    supporting_engines: List[str] = field(default_factory=list)

    # -- convenience delegates -------------------------------------------------

    @property
    def position(self) -> int:
        return self.mutation.position

    @property
    def original(self) -> str:
        return self.mutation.original

    @property
    def mutant(self) -> str:
        return self.mutation.mutant

    @property
    def delta_score(self) -> float:
        return self.mutation.delta_score

    @property
    def score(self) -> float:
        """Alias for delta_score (backward compat)."""
        return self.mutation.delta_score

    @property
    def score_type(self) -> str:
        return self.mutation.score_type

    @property
    def engine(self) -> str:
        return self.mutation.engine

    @property
    def recommendation(self) -> str:
        return self.mutation.recommendation

    @property
    def description(self) -> str:
        return self.mutation.description

    def to_mutation_result(self) -> MutationResult:
        """Extract the underlying MutationResult (strips suggestion metadata)."""
        if isinstance(self.mutation, MutationResult):
            # If it's a Mutation (subclass), strip to plain MutationResult
            if type(self.mutation) is not MutationResult:
                return self.mutation.to_mutation_result()
            return self.mutation
        # Shouldn't happen, but be safe
        return self.mutation  # type: ignore[return-value]

    def __str__(self) -> str:
        eng = self.mutation.engine
        return (
            f"#{self.rank} {self.mutation.original}{self.mutation.position + 1}"
            f"{self.mutation.mutant} ({eng}: {self.score_type}="
            f"{self.delta_score:.2f}, conf={self.confidence:.2f})"
        )


# ---------------------------------------------------------------------------
# Type alias for functions accepting either Mutation or MutationResult
# ---------------------------------------------------------------------------

AnyMutation = Union[MutationResult, Mutation]


# ---------------------------------------------------------------------------
# rank_mutations — sort by impact
# ---------------------------------------------------------------------------

def rank_mutations(
    mutations: Sequence[AnyMutation],
    *,
    descending: bool = True,
    score_type: Optional[str] = None,
    engine: Optional[str] = None,
) -> List[MutationSuggestion]:
    """Rank mutations by delta_score magnitude.

    Works with both Mutation and MutationResult (and any subclass).

    Args:
        mutations: sequence of Mutation / MutationResult objects
        descending: if True, highest delta_score first (default);
                    if False, lowest first (e.g. most stabilizing ΔΔG)
        score_type: if given, only include mutations matching this score_type
        engine: if given, only include mutations matching this engine

    Returns:
        List of MutationSuggestion objects, ordered by rank.
    """
    # Filter
    filtered: List[MutationResult] = []
    for m in mutations:
        # Ensure we have a MutationResult-compatible object
        mr = m if isinstance(m, MutationResult) else m  # type: ignore[redundant-expr]
        if score_type is not None and mr.score_type != score_type:
            continue
        if engine is not None and mr.engine != engine:
            continue
        filtered.append(mr)

    # Sort by delta_score
    filtered.sort(key=lambda m: m.delta_score, reverse=descending)

    # Wrap in MutationSuggestion with rank
    results: List[MutationSuggestion] = []
    for i, mr in enumerate(filtered):
        confidence = 1.0
        # If it's a Mutation, use its confidence field
        if isinstance(mr, Mutation):
            confidence = mr.confidence

        results.append(MutationSuggestion(
            mutation=mr,
            rank=i + 1,
            confidence=confidence,
            rationale=mr.recommendation or mr.description or "",
            supporting_engines=[mr.engine],
        ))

    return results


# ---------------------------------------------------------------------------
# combine_mutations — merge from multiple engines
# ---------------------------------------------------------------------------

def combine_mutations(
    *mutation_lists: Sequence[AnyMutation],
    deduplicate: bool = True,
    descending: bool = True,
) -> List[MutationSuggestion]:
    """Combine mutation lists from multiple engines, deduplicating by position.

    Works with both Mutation and MutationResult (and any subclass).

    When the same position+original+mutant appears from multiple engines:
      - The delta_scores are averaged
      - All engine names are collected in supporting_engines
      - The confidence is the maximum confidence across engines

    Args:
        *mutation_lists: one or more sequences of Mutation/MutationResult
        deduplicate: if True, merge duplicates (same position+original+mutant)
        descending: if True, rank highest delta_score first

    Returns:
        List of MutationSuggestion objects, ordered by rank.
    """
    # Flatten all lists into MutationResult objects
    all_mutations: List[MutationResult] = []
    for ml in mutation_lists:
        for m in ml:
            if isinstance(m, MutationResult):
                all_mutations.append(m)
            else:
                # Shouldn't happen, but handle gracefully
                continue

    if not deduplicate:
        return rank_mutations(all_mutations, descending=descending)

    # Deduplicate: group by (position, original, mutant)
    from collections import defaultdict

    groups: dict[tuple[int, str, str], List[MutationResult]] = defaultdict(list)
    for mr in all_mutations:
        key = (mr.position, mr.original, mr.mutant)
        groups[key].append(mr)

    # Merge duplicates
    merged: List[MutationResult] = []
    for key, group in groups.items():
        if len(group) == 1:
            merged.append(group[0])
        else:
            # Average delta_score
            avg_delta = sum(m.delta_score for m in group) / len(group)
            # Collect all engines
            engines = list(dict.fromkeys(m.engine for m in group))  # unique, order-preserving
            # Use first as template, override merged fields
            template = group[0]
            merged_details: dict = dict(template.details)
            if len(engines) > 1:
                merged_details["merged_engines"] = engines
                merged_details["original_delta_scores"] = {m.engine: m.delta_score for m in group}

            merged.append(MutationResult(
                position=template.position,
                original=template.original,
                mutant=template.mutant,
                delta_score=round(avg_delta, 4),
                score_type=template.score_type,
                engine=engines[0] if len(engines) == 1 else "combined",
                recommendation=template.recommendation,
                description=template.description,
                details=merged_details,
            ))

    # Rank the merged list
    ranked = rank_mutations(merged, descending=descending)

    # Update supporting_engines for deduplicated entries
    if all_mutations:
        for suggestion in ranked:
            key = (suggestion.position, suggestion.original, suggestion.mutant)
            group = groups.get(key, [])
            if len(group) > 1:
                engines = list(dict.fromkeys(m.engine for m in group))
                suggestion.supporting_engines = engines
                # Max confidence across engines
                max_conf = 1.0
                for m in group:
                    if isinstance(m, Mutation):
                        max_conf = max(max_conf, m.confidence)
                suggestion.confidence = max_conf

    return ranked


# ---------------------------------------------------------------------------
# Helper: ensure a value is a MutationResult
# ---------------------------------------------------------------------------

def as_mutation_result(m: AnyMutation) -> MutationResult:
    """Convert any mutation-like object to a plain MutationResult.

    If already a MutationResult (including Mutation), returns as-is.
    """
    if isinstance(m, Mutation) and type(m) is Mutation:
        return m.to_mutation_result()
    if isinstance(m, MutationResult):
        return m
    raise TypeError(f"Expected MutationResult or Mutation, got {type(m).__name__}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    "Mutation",
    "MutationSuggestion",
    "AnyMutation",
    "rank_mutations",
    "combine_mutations",
    "as_mutation_result",
]
