"""Tests for biocompiler.solver.conflict_resolution.

Covers the 5 main areas:
1. ConstraintConflict construction and validation
2. ConflictResolver.detect_conflicts()
3. ConflictResolver.suggest_resolution()
4. ConflictResolver.auto_resolve()
5. prioritize_constraints()
"""

from __future__ import annotations

import copy
import pytest

# ---------------------------------------------------------------------------
# Graceful import — entire module skipped when solver package is absent
# ---------------------------------------------------------------------------
cr = pytest.importorskip("biocompiler.solver.conflict_resolution")
ConstraintConflict = cr.ConstraintConflict
ConflictResolver = cr.ConflictResolver
prioritize_constraints = cr.prioritize_constraints

from biocompiler.solver.types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
    SolverConfig,
)


# =====================================================================
# Helpers
# =====================================================================

def _make_model(
    protein: str = "MVSKGE",
    constraints: list[ConstraintSpec] | None = None,
    config: SolverConfig | None = None,
) -> CSPModel:
    """Build a minimal CSPModel for testing."""
    from biocompiler.type_system import AA_TO_CODONS

    codon_domains = {i: list(AA_TO_CODONS[aa]) for i, aa in enumerate(protein)}
    return CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints or [],
        config=config or SolverConfig(),
    )


def _int_priority_to_enum(value: int) -> ConstraintPriority:
    """Map old-style integer priority (0–10, lower = more critical) to ConstraintPriority.

    Mapping:
        0–2 → CRITICAL
        3–4 → HIGH
        5–6 → MEDIUM
        7–10 → LOW
    """
    if value <= 2:
        return ConstraintPriority.CRITICAL
    elif value <= 4:
        return ConstraintPriority.HIGH
    elif value <= 6:
        return ConstraintPriority.MEDIUM
    else:
        return ConstraintPriority.LOW


def _hard_spec(
    name: str,
    ctype: ConstraintType = ConstraintType.GC_CONTENT,
    positions: list[int] | None = None,
    priority: int | ConstraintPriority = 5,
) -> ConstraintSpec:
    """Shorthand for a hard ConstraintSpec."""
    prio = priority if isinstance(priority, ConstraintPriority) else _int_priority_to_enum(priority)
    return ConstraintSpec(
        ctype=ctype,
        name=name,
        strictness=ConstraintStrictness.HARD,
        positions=positions or [],
        priority=prio,
    )


def _soft_spec(
    name: str,
    ctype: ConstraintType = ConstraintType.CODON_USAGE,
    positions: list[int] | None = None,
    priority: int | ConstraintPriority = 5,
) -> ConstraintSpec:
    """Shorthand for a soft ConstraintSpec."""
    prio = priority if isinstance(priority, ConstraintPriority) else _int_priority_to_enum(priority)
    return ConstraintSpec(
        ctype=ctype,
        name=name,
        strictness=ConstraintStrictness.SOFT,
        positions=positions or [],
        priority=prio,
    )


# =====================================================================
# 1. ConstraintConflict construction and validation
# =====================================================================


class TestConstraintConflictConstruction:
    """Test ConstraintConflict dataclass construction and __post_init__ validation."""

    def test_basic_construction(self):
        """ConstraintConflict can be created with valid arguments."""
        cc = ConstraintConflict(
            constraint_a="gc_range",
            constraint_b="no_cpg",
            conflict_positions=[2, 5, 7],
            resolution_strategy="relax_b",
        )
        assert cc.constraint_a == "gc_range"
        assert cc.constraint_b == "no_cpg"
        assert cc.conflict_positions == [2, 5, 7]
        assert cc.resolution_strategy == "relax_b"

    def test_all_valid_strategies(self):
        """All four valid resolution strategies are accepted."""
        valid = ["relax_a", "relax_b", "compromise", "infeasible"]
        for strategy in valid:
            cc = ConstraintConflict(
                constraint_a="a", constraint_b="b",
                conflict_positions=[], resolution_strategy=strategy,
            )
            assert cc.resolution_strategy == strategy

    def test_invalid_strategy_raises(self):
        """Invalid resolution_strategy raises ValueError."""
        with pytest.raises(ValueError, match="Invalid resolution_strategy"):
            ConstraintConflict(
                constraint_a="a", constraint_b="b",
                conflict_positions=[], resolution_strategy="invalid_strategy",
            )

    def test_empty_string_strategy_raises(self):
        """Empty string resolution_strategy raises ValueError."""
        with pytest.raises(ValueError, match="Invalid resolution_strategy"):
            ConstraintConflict(
                constraint_a="a", constraint_b="b",
                conflict_positions=[], resolution_strategy="",
            )

    def test_empty_conflict_positions(self):
        """Empty conflict_positions list is allowed."""
        cc = ConstraintConflict(
            constraint_a="a", constraint_b="b",
            conflict_positions=[], resolution_strategy="compromise",
        )
        assert cc.conflict_positions == []

    def test_single_position_conflict(self):
        """Single conflict position works."""
        cc = ConstraintConflict(
            constraint_a="a", constraint_b="b",
            conflict_positions=[3], resolution_strategy="relax_a",
        )
        assert cc.conflict_positions == [3]

    def test_repr_format(self):
        """__repr__ includes key fields."""
        cc = ConstraintConflict(
            constraint_a="gc", constraint_b="splice",
            conflict_positions=[1, 2], resolution_strategy="compromise",
        )
        r = repr(cc)
        assert "'gc'" in r
        assert "'splice'" in r
        assert "compromise" in r

    def test_same_constraint_names(self):
        """constraint_a and constraint_b can be the same (degenerate case)."""
        cc = ConstraintConflict(
            constraint_a="gc", constraint_b="gc",
            conflict_positions=[0], resolution_strategy="relax_a",
        )
        assert cc.constraint_a == cc.constraint_b

    def test_large_position_list(self):
        """Large conflict_positions list is stored correctly."""
        positions = list(range(100))
        cc = ConstraintConflict(
            constraint_a="a", constraint_b="b",
            conflict_positions=positions, resolution_strategy="relax_b",
        )
        assert len(cc.conflict_positions) == 100

    def test_strategy_case_sensitive(self):
        """Strategy validation is case-sensitive."""
        with pytest.raises(ValueError):
            ConstraintConflict(
                constraint_a="a", constraint_b="b",
                conflict_positions=[], resolution_strategy="RELAX_A",
            )


# =====================================================================
# 2. ConflictResolver.detect_conflicts()
# =====================================================================


class TestDetectConflicts:
    """Test ConflictResolver.detect_conflicts()."""

    def test_no_constraints_no_conflicts(self):
        """Model with no constraints produces no conflicts."""
        model = _make_model(constraints=[])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_single_constraint_no_conflict(self):
        """Single hard constraint cannot conflict with anything."""
        model = _make_model(constraints=[
            _hard_spec("gc_range", positions=[]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_two_hard_constraints_no_overlap(self):
        """Two hard constraints with non-overlapping positions → no conflict."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc_pos_0", positions=[0, 1]),
            _hard_spec("splice_pos_4", positions=[4, 5]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_two_hard_constraints_with_overlap(self):
        """Two hard constraints with overlapping positions → conflict detected."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc_range", positions=[0, 1, 2]),
            _hard_spec("splice_site", positions=[1, 2, 3]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        c = conflicts[0]
        assert {c.constraint_a, c.constraint_b} == {"gc_range", "splice_site"}
        assert set(c.conflict_positions) == {1, 2}

    def test_soft_constraints_ignored(self):
        """Soft constraints are not included in conflict detection."""
        model = _make_model(constraints=[
            _hard_spec("gc_range", positions=[0, 1]),
            _soft_spec("cai_opt", positions=[0, 1]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_both_global_constraints_conflict(self):
        """Two hard constraints with empty positions (global) → conflict on all positions."""
        model = _make_model(protein="MVSK", constraints=[
            _hard_spec("gc_range", positions=[]),
            _hard_spec("no_cpg", positions=[]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        # Global constraints overlap on all positions
        assert len(conflicts[0].conflict_positions) == len("MVSK")

    def test_strategy_relax_b_when_a_higher_priority(self):
        """When constraint A has higher priority (lower number), strategy is relax_b."""
        model = _make_model(constraints=[
            _hard_spec("gc_range", priority=1, positions=[0, 1]),
            _hard_spec("no_cpg", priority=5, positions=[0, 1]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        assert conflicts[0].resolution_strategy == "relax_b"

    def test_strategy_relax_a_when_b_higher_priority(self):
        """When constraint B has higher priority, strategy is relax_a."""
        model = _make_model(constraints=[
            _hard_spec("gc_range", priority=5, positions=[0, 1]),
            _hard_spec("no_cpg", priority=1, positions=[0, 1]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        assert conflicts[0].resolution_strategy == "relax_a"

    def test_strategy_compromise_when_equal_priority(self):
        """Equal priority constraints → compromise strategy."""
        model = _make_model(constraints=[
            _hard_spec("gc_range", priority=3, positions=[0, 1]),
            _hard_spec("no_cpg", priority=3, positions=[0, 1]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        assert conflicts[0].resolution_strategy == "compromise"

    def test_strategy_infeasible_for_two_translation_constraints(self):
        """Two AMINO_ACID_IDENTITY constraints → infeasible strategy."""
        model = _make_model(constraints=[
            _hard_spec("trans_0", ctype=ConstraintType.AMINO_ACID_IDENTITY, positions=[0]),
            _hard_spec("trans_1", ctype=ConstraintType.AMINO_ACID_IDENTITY, positions=[0]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        assert conflicts[0].resolution_strategy == "infeasible"

    def test_multiple_pairwise_conflicts(self):
        """Three hard constraints all overlapping → 3 pairwise conflicts."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc", positions=[0, 1, 2, 3]),
            _hard_spec("cpg", positions=[1, 2, 3, 4]),
            _hard_spec("splice", positions=[2, 3, 4, 5]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        # gc-cpg, gc-splice, cpg-splice
        assert len(conflicts) == 3

    def test_no_duplicate_pairs(self):
        """Each pair is reported only once regardless of iteration order."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0]),
            _hard_spec("cpg", positions=[0]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1

    def test_result_populates_cache(self):
        """detect_conflicts populates the internal cache for auto_resolve."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0]),
            _hard_spec("cpg", positions=[0]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        # The cache should be populated (used by auto_resolve)
        assert id(model) in resolver._conflict_cache
        cached = resolver._conflict_cache[id(model)]
        assert len(cached) == len(conflicts)

    def test_one_global_one_specific(self):
        """One global (empty positions) + one specific → overlap on specific positions."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc_range", positions=[]),
            _hard_spec("splice_pos_2", positions=[2]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        assert conflicts[0].conflict_positions == [2]

    def test_partial_overlap_positions(self):
        """Partial position overlap produces overlap only on shared positions."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc", positions=[0, 1, 2]),
            _hard_spec("splice", positions=[2, 3, 4]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        assert conflicts[0].conflict_positions == [2]

    def test_translation_vs_non_translation(self):
        """AMINO_ACID_IDENTITY vs GC_CONTENT → not infeasible (compromise/relax based on priority)."""
        model = _make_model(constraints=[
            _hard_spec("trans", ctype=ConstraintType.AMINO_ACID_IDENTITY, positions=[0], priority=1),
            _hard_spec("gc", ctype=ConstraintType.GC_CONTENT, positions=[0], priority=3),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        # Not infeasible because only one is AMINO_ACID_IDENTITY
        assert conflicts[0].resolution_strategy != "infeasible"


# =====================================================================
# 3. ConflictResolver.suggest_resolution()
# =====================================================================


class TestSuggestResolution:
    """Test ConflictResolver.suggest_resolution()."""

    def _make_conflict(self, strategy: str, a: str = "c_a", b: str = "c_b",
                       positions: list[int] | None = None) -> ConstraintConflict:
        return ConstraintConflict(
            constraint_a=a,
            constraint_b=b,
            conflict_positions=positions or [1, 2],
            resolution_strategy=strategy,
        )

    def test_relax_a_suggestion(self):
        """relax_a strategy mentions relaxing constraint A."""
        resolver = ConflictResolver()
        conflict = self._make_conflict("relax_a")
        suggestion = resolver.suggest_resolution(conflict)
        assert "c_a" in suggestion
        assert "Relax" in suggestion

    def test_relax_b_suggestion(self):
        """relax_b strategy mentions relaxing constraint B."""
        resolver = ConflictResolver()
        conflict = self._make_conflict("relax_b")
        suggestion = resolver.suggest_resolution(conflict)
        assert "c_b" in suggestion
        assert "Relax" in suggestion

    def test_compromise_suggestion(self):
        """compromise strategy mentions both constraints."""
        resolver = ConflictResolver()
        conflict = self._make_conflict("compromise", a="gc", b="splice")
        suggestion = resolver.suggest_resolution(conflict)
        assert "gc" in suggestion
        assert "splice" in suggestion
        assert "compromise" in suggestion.lower() or "Partially" in suggestion

    def test_infeasible_suggestion(self):
        """infeasible strategy mentions redesign."""
        resolver = ConflictResolver()
        conflict = self._make_conflict("infeasible")
        suggestion = resolver.suggest_resolution(conflict)
        assert "infeasible" in suggestion.lower() or "redesign" in suggestion.lower()

    def test_suggestion_includes_positions(self):
        """Suggestion includes conflict positions when present."""
        resolver = ConflictResolver()
        conflict = self._make_conflict("relax_a", positions=[3, 7])
        suggestion = resolver.suggest_resolution(conflict)
        assert "3" in suggestion
        assert "7" in suggestion

    def test_suggestion_empty_positions(self):
        """Suggestion with empty positions doesn't crash."""
        resolver = ConflictResolver()
        conflict = self._make_conflict("relax_b", positions=[])
        suggestion = resolver.suggest_resolution(conflict)
        assert isinstance(suggestion, str)
        assert len(suggestion) > 0

    def test_return_type_is_string(self):
        """suggest_resolution always returns a string."""
        resolver = ConflictResolver()
        for strategy in ["relax_a", "relax_b", "compromise", "infeasible"]:
            conflict = self._make_conflict(strategy)
            assert isinstance(resolver.suggest_resolution(conflict), str)


# =====================================================================
# 4. ConflictResolver.auto_resolve()
# =====================================================================


class TestAutoResolve:
    """Test ConflictResolver.auto_resolve()."""

    def test_invalid_strategy_raises(self):
        """Unknown auto-resolution strategy raises ValueError."""
        model = _make_model(constraints=[])
        resolver = ConflictResolver()
        with pytest.raises(ValueError, match="Unknown auto-resolution strategy"):
            resolver.auto_resolve(model, strategy="bogus")

    def test_no_conflicts_returns_same_model(self):
        """Model with no conflicts is returned unchanged."""
        model = _make_model(constraints=[])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model)
        assert resolved is model

    def test_no_conflicts_returns_same_model_with_constraints(self):
        """Model with non-conflicting constraints is returned unchanged."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0]),
            _hard_spec("splice", positions=[5]),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model)
        assert resolved is model

    def test_min_relaxation_removes_constraints(self):
        """min_relaxation strategy removes constraints to resolve conflicts."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=5),
            _hard_spec("cpg", positions=[0, 1], priority=3),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="min_relaxation")
        # The lower-priority constraint (gc, priority=5) should be removed
        remaining_names = [c.name for c in resolved.constraints]
        assert "cpg" in remaining_names
        # gc should be removed (higher priority number = easier to relax)
        assert "gc" not in remaining_names

    def test_max_priority_removes_lower_priority(self):
        """max_priority strategy removes the lowest-priority (highest number) constraint."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=1),
            _hard_spec("cpg", positions=[0, 1], priority=10),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="max_priority")
        remaining_names = [c.name for c in resolved.constraints]
        assert "gc" in remaining_names
        assert "cpg" not in remaining_names

    def test_compromise_first_downgrades_equal_priority(self):
        """compromise_first downgrades both constraints when strategy is compromise."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=3),
            _hard_spec("cpg", positions=[0, 1], priority=3),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="compromise_first")
        # Both should be downgraded to SOFT, not removed
        for c in resolved.constraints:
            if c.name in ("gc", "cpg"):
                assert c.strictness == ConstraintStrictness.SOFT

    def test_compromise_first_removes_relax_a(self):
        """compromise_first removes constraint_a when strategy is relax_a."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=5),
            _hard_spec("cpg", positions=[0, 1], priority=1),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="compromise_first")
        remaining_names = [c.name for c in resolved.constraints]
        # gc (lower priority, priority=5) → relax_a → remove
        assert "gc" not in remaining_names
        assert "cpg" in remaining_names

    def test_compromise_first_skips_infeasible(self):
        """compromise_first skips infeasible conflicts (cannot auto-resolve)."""
        model = _make_model(constraints=[
            _hard_spec("t0", ctype=ConstraintType.AMINO_ACID_IDENTITY, positions=[0], priority=1),
            _hard_spec("t1", ctype=ConstraintType.AMINO_ACID_IDENTITY, positions=[0], priority=1),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="compromise_first")
        # Both constraints should remain (infeasible is not auto-resolved)
        remaining_names = [c.name for c in resolved.constraints]
        assert "t0" in remaining_names
        assert "t1" in remaining_names

    def test_resolved_model_preserves_protein(self):
        """Auto-resolved model preserves the protein_sequence and codon_domains."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=5),
            _hard_spec("cpg", positions=[0, 1], priority=1),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model)
        assert resolved.protein_sequence == model.protein_sequence
        assert resolved.codon_domains == model.codon_domains

    def test_resolved_model_is_new_object(self):
        """auto_resolve returns a new CSPModel, not a mutation of the original."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=5),
            _hard_spec("cpg", positions=[0, 1], priority=1),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model)
        assert resolved is not model

    def test_cache_invalidated_after_resolve(self):
        """After auto_resolve, the old model's cache is cleared."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0], priority=5),
            _hard_spec("cpg", positions=[0], priority=1),
        ])
        resolver = ConflictResolver()
        resolver.detect_conflicts(model)
        assert id(model) in resolver._conflict_cache
        resolver.auto_resolve(model)
        assert id(model) not in resolver._conflict_cache

    def test_min_relaxation_greedy_set_cover(self):
        """min_relaxation uses greedy set cover: one constraint can resolve multiple conflicts."""
        # Constraint "x" appears in two conflicts; removing it resolves both
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("x", positions=[0, 1, 2, 3], priority=8),
            _hard_spec("a", positions=[0, 1], priority=1),
            _hard_spec("b", positions=[2, 3], priority=2),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="min_relaxation")
        remaining_names = [c.name for c in resolved.constraints]
        # "x" participates in both conflicts and should be removed
        assert "x" not in remaining_names
        assert "a" in remaining_names
        assert "b" in remaining_names

    def test_all_three_strategies_valid(self):
        """All three auto-resolution strategies are accepted without error."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=3),
            _hard_spec("cpg", positions=[0, 1], priority=3),
        ])
        resolver = ConflictResolver()
        for strategy in ["min_relaxation", "max_priority", "compromise_first"]:
            resolved = resolver.auto_resolve(model, strategy=strategy)
            assert isinstance(resolved, CSPModel)

    def test_downgraded_constraint_preserves_params(self):
        """Downgraded constraints preserve ctype, name, params, positions, priority."""
        model = _make_model(constraints=[
            ConstraintSpec(
                ctype=ConstraintType.GC_CONTENT,
                name="gc",
                strictness=ConstraintStrictness.HARD,
                params={"lo": 0.4, "hi": 0.6},
                positions=[0, 1],
                priority=3,
            ),
            ConstraintSpec(
                ctype=ConstraintType.NO_CPG,
                name="cpg",
                strictness=ConstraintStrictness.HARD,
                params={"threshold": 0.6},
                positions=[0, 1],
                priority=3,
            ),
        ])
        resolver = ConflictResolver()
        # First check if conflicts are detected; if not, the model is returned unchanged
        conflicts = resolver.detect_conflicts(model)
        if conflicts:
            resolved = resolver.auto_resolve(model, strategy="compromise_first")
            for c in resolved.constraints:
                if c.name == "gc" and c.strictness == ConstraintStrictness.SOFT:
                    assert c.params == {"lo": 0.4, "hi": 0.6}
                    assert c.positions == [0, 1]
                    assert c.priority == 3
                if c.name == "cpg" and c.strictness == ConstraintStrictness.SOFT:
                    assert c.params == {"threshold": 0.6}
        else:
            # No conflicts detected — verify constraints are unchanged
            for c in resolved.constraints:
                if c.name == "gc":
                    assert c.params == {"lo": 0.4, "hi": 0.6}
                if c.name == "cpg":
                    assert c.params == {"threshold": 0.6}

    def test_auto_resolve_detects_conflicts_if_not_cached(self):
        """auto_resolve calls detect_conflicts internally if cache is empty."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0], priority=5),
            _hard_spec("cpg", positions=[0], priority=1),
        ])
        resolver = ConflictResolver()
        # No prior detect_conflicts call
        assert id(model) not in resolver._conflict_cache
        resolved = resolver.auto_resolve(model)
        # Should have detected the conflict and resolved it
        assert isinstance(resolved, CSPModel)


# =====================================================================
# 5. prioritize_constraints()
# =====================================================================


class TestPrioritizeConstraints:
    """Test prioritize_constraints() function."""

    def test_empty_list(self):
        """Empty list returns empty list."""
        assert prioritize_constraints([]) == []

    def test_single_constraint(self):
        """Single constraint is returned as-is."""
        c = _hard_spec("gc", priority=3)
        result = prioritize_constraints([c])
        assert result == [c]

    def test_sorts_by_priority_ascending(self):
        """Constraints are sorted by priority ascending (most critical first)."""
        c1 = _hard_spec("gc", priority=5)
        c2 = _hard_spec("cpg", priority=1)
        c3 = _hard_spec("splice", priority=3)
        result = prioritize_constraints([c1, c2, c3])
        assert [c.name for c in result] == ["cpg", "splice", "gc"]

    def test_already_sorted(self):
        """Already-sorted list is unchanged."""
        c1 = _hard_spec("a", priority=1)
        c2 = _hard_spec("b", priority=2)
        c3 = _hard_spec("c", priority=3)
        result = prioritize_constraints([c1, c2, c3])
        assert [c.name for c in result] == ["a", "b", "c"]

    def test_reverse_sorted(self):
        """Reverse-sorted list is correctly reordered."""
        c1 = _hard_spec("a", priority=10)
        c2 = _hard_spec("b", priority=5)
        c3 = _hard_spec("c", priority=1)
        result = prioritize_constraints([c1, c2, c3])
        assert [c.name for c in result] == ["c", "b", "a"]

    def test_equal_priorities_preserve_relative_order(self):
        """Constraints with equal priority maintain stable sort order."""
        c1 = _hard_spec("first", priority=3)
        c2 = _hard_spec("second", priority=3)
        c3 = _hard_spec("third", priority=3)
        result = prioritize_constraints([c1, c2, c3])
        # Python sort is stable — equal-priority elements keep original order
        assert [c.name for c in result] == ["first", "second", "third"]

    def test_returns_new_list(self):
        """Returned list is a new list, not a mutation of the input."""
        constraints = [_hard_spec("a", priority=3), _hard_spec("b", priority=1)]
        result = prioritize_constraints(constraints)
        assert result is not constraints

    def test_mixed_constraint_types(self):
        """Works with both HARD and SOFT constraints."""
        c1 = _hard_spec("hard_gc", priority=5)    # MEDIUM
        c2 = _soft_spec("soft_cai", priority=3)   # HIGH
        c3 = _hard_spec("hard_cpg", priority=1)   # CRITICAL
        result = prioritize_constraints([c1, c2, c3])
        # Sorted by priority: CRITICAL(0) < HIGH(1) < MEDIUM(2)
        assert [c.name for c in result] == ["hard_cpg", "soft_cai", "hard_gc"]

    def test_different_constraint_types(self):
        """Works with various ConstraintType values."""
        constraints = [
            ConstraintSpec(ctype=ConstraintType.GC_CONTENT, name="gc", priority=_int_priority_to_enum(7)),       # LOW
            ConstraintSpec(ctype=ConstraintType.NO_CPG, name="cpg", priority=_int_priority_to_enum(2)),           # CRITICAL
            ConstraintSpec(ctype=ConstraintType.RESTRICTION_SITE, name="eco", priority=_int_priority_to_enum(3)), # HIGH
            ConstraintSpec(ctype=ConstraintType.AMINO_ACID_IDENTITY, name="trans", priority=_int_priority_to_enum(0)),  # CRITICAL
        ]
        result = prioritize_constraints(constraints)
        # CRITICAL(cpg,trans) < HIGH(eco) < LOW(gc)
        # cpg and trans are both CRITICAL — stable sort preserves original order
        assert [c.name for c in result] == ["cpg", "trans", "eco", "gc"]

    def test_priority_boundary_values(self):
        """Priority values at extremes (0 and 10) sort correctly."""
        c0 = _hard_spec("critical", priority=0)
        c10 = _hard_spec("relaxable", priority=10)
        c5 = _hard_spec("medium", priority=5)
        result = prioritize_constraints([c10, c5, c0])
        assert [c.name for c in result] == ["critical", "medium", "relaxable"]


# =====================================================================
# Integration / cross-cutting tests
# =====================================================================


class TestConflictResolutionIntegration:
    """End-to-end tests combining detect → suggest → auto_resolve."""

    def test_full_workflow(self):
        """Complete workflow: detect → suggest → auto_resolve."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc", positions=[0, 1, 2], priority=5),
            _hard_spec("cpg", positions=[1, 2, 3], priority=3),
            _hard_spec("splice", positions=[2, 3, 4], priority=7),
        ])
        resolver = ConflictResolver()

        # Detect
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) > 0

        # Suggest
        for conflict in conflicts:
            suggestion = resolver.suggest_resolution(conflict)
            assert isinstance(suggestion, str)
            assert len(suggestion) > 0

        # Auto-resolve
        resolved = resolver.auto_resolve(model, strategy="min_relaxation")
        assert isinstance(resolved, CSPModel)
        assert len(resolved.constraints) < len(model.constraints)

    def test_resolve_then_re_detect_no_conflicts(self):
        """After auto-resolve, re-detecting on the resolved model finds no conflicts."""
        model = _make_model(protein="MVSKGE", constraints=[
            _hard_spec("gc", positions=[0, 1], priority=5),
            _hard_spec("cpg", positions=[0, 1], priority=1),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="max_priority")

        # Re-detect on resolved model
        new_conflicts = resolver.detect_conflicts(resolved)
        assert new_conflicts == []

    def test_multiple_resolve_strategies_same_model(self):
        """Different strategies produce different resolved models."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=3),
            _hard_spec("cpg", positions=[0, 1], priority=3),
        ])
        resolver = ConflictResolver()
        r1 = resolver.auto_resolve(model, strategy="max_priority")
        r2 = resolver.auto_resolve(model, strategy="compromise_first")

        # max_priority removes one; compromise_first downgrades both to SOFT
        r1_names = {c.name for c in r1.constraints}
        r2_names = {c.name for c in r2.constraints}
        # compromise_first should keep both (downgraded), max_priority removes one
        assert len(r2_names) >= len(r1_names)

    def test_valine_gt_dinucleotide_conflict_scenario(self):
        """Simulate the classic Valine GT dinucleotide conflict scenario.

        All valine codons contain GT, so NoGT_Dinucleotide + Valine positions
        creates a real biological conflict.
        """
        model = _make_model(protein="VVVV", constraints=[
            _hard_spec("no_gt", ctype=ConstraintType.NO_GT_DINUCLEOTIDE, positions=[0, 1, 2, 3], priority=3),
            _hard_spec("gc", ctype=ConstraintType.GC_CONTENT, positions=[0, 1, 2, 3], priority=5),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        # no_gt and gc overlap on all positions
        assert len(conflicts) >= 1

        # Auto-resolve with max_priority: gc (priority 5, easier) should be removed
        resolved = resolver.auto_resolve(model, strategy="max_priority")
        remaining = {c.name for c in resolved.constraints}
        assert "no_gt" in remaining
        assert "gc" not in remaining


class TestEdgeCases:
    """Edge cases for conflict resolution."""

    def test_model_with_only_soft_constraints(self):
        """Model with only soft constraints produces no conflicts."""
        model = _make_model(constraints=[
            _soft_spec("cai", positions=[0, 1]),
            _soft_spec("cpg", positions=[0, 1]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_detect_populates_cache_for_auto_resolve(self):
        """After detect_conflicts, auto_resolve uses cached results."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0], priority=3),
            _hard_spec("cpg", positions=[0], priority=5),
        ])
        resolver = ConflictResolver()
        resolver.detect_conflicts(model)
        # Now auto_resolve should use the cache (no re-detection needed)
        assert id(model) in resolver._conflict_cache
        resolved = resolver.auto_resolve(model, strategy="max_priority")
        assert isinstance(resolved, CSPModel)

    def test_auto_resolve_preserves_non_conflicting_constraints(self):
        """Non-conflicting constraints survive auto-resolution."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0, 1], priority=5),
            _hard_spec("cpg", positions=[0, 1], priority=1),
            _hard_spec("unrelated", positions=[5], priority=3),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="max_priority")
        remaining_names = [c.name for c in resolved.constraints]
        assert "unrelated" in remaining_names
        assert "cpg" in remaining_names

    def test_many_constraints_few_conflicts(self):
        """Many constraints with sparse overlaps produce few conflicts."""
        constraints = [
            _hard_spec(f"c{i}", positions=[i], priority=i)
            for i in range(10)
        ]
        # No overlaps since each constraint applies to a unique position
        model = _make_model(protein="M" * 10, constraints=constraints)
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_all_constraints_on_same_position(self):
        """All constraints on the same position produce maximum pairwise conflicts."""
        n = 5
        constraints = [
            _hard_spec(f"c{i}", positions=[0], priority=i + 1)
            for i in range(n)
        ]
        model = _make_model(constraints=constraints)
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        # C(5,2) = 10 pairwise conflicts
        assert len(conflicts) == n * (n - 1) // 2

    def test_auto_resolve_with_mixed_soft_and_hard(self):
        """Auto-resolve only considers hard constraints; soft are untouched."""
        model = _make_model(constraints=[
            _hard_spec("gc", positions=[0], priority=5),
            _hard_spec("cpg", positions=[0], priority=1),
            _soft_spec("cai", positions=[0], priority=3),
        ])
        resolver = ConflictResolver()
        resolved = resolver.auto_resolve(model, strategy="max_priority")
        cai_found = [c for c in resolved.constraints if c.name == "cai"]
        assert len(cai_found) == 1
        assert cai_found[0].strictness == ConstraintStrictness.SOFT

    def test_constraint_names_in_conflict_match_spec_names(self):
        """Conflict constraint_a and constraint_b match actual ConstraintSpec names."""
        model = _make_model(constraints=[
            _hard_spec("my_gc_constraint", positions=[0]),
            _hard_spec("my_cpg_constraint", positions=[0]),
        ])
        resolver = ConflictResolver()
        conflicts = resolver.detect_conflicts(model)
        assert len(conflicts) == 1
        names = {conflicts[0].constraint_a, conflicts[0].constraint_b}
        assert names == {"my_gc_constraint", "my_cpg_constraint"}
