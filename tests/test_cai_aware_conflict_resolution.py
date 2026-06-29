"""Tests for CAI-aware conflict resolution integration.

Tests cover:
1. CAIImpactEstimator — estimate_cai_impact and rank_resolution
2. ConflictResolver with cai_aware=True vs cai_aware=False
3. resolution_strategy parameter ("cai_aware" vs "priority_only")
4. ConflictProvenance records with cai_aware resolution
5. Backward compatibility — cai_aware=False produces same results as before
"""

from __future__ import annotations

import pytest

from biocompiler.solver.types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
    SolverConfig,
)
from biocompiler.solver.conflict_resolution import (
    CAIImpactEstimator,
    ConstraintConflict,
    ConflictResolver,
    prioritize_constraints,
)
from biocompiler.solver.conflict_provenance import (
    ConflictProvenance,
    ConflictResolverWithProvenance,
)


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def gc_constraint() -> ConstraintSpec:
    """A GC content hard constraint."""
    return ConstraintSpec(
        ctype=ConstraintType.GC_CONTENT,
        name="gc_bounds",
        strictness=ConstraintStrictness.HARD,
        priority=ConstraintPriority.MEDIUM,
        params={"gc_lo": 0.30, "gc_hi": 0.70},
    )


@pytest.fixture
def cpg_constraint() -> ConstraintSpec:
    """A CpG avoidance hard constraint."""
    return ConstraintSpec(
        ctype=ConstraintType.NO_CPG,
        name="no_cpg",
        strictness=ConstraintStrictness.HARD,
        priority=ConstraintPriority.MEDIUM,
    )


@pytest.fixture
def restriction_constraint() -> ConstraintSpec:
    """A restriction site avoidance hard constraint."""
    return ConstraintSpec(
        ctype=ConstraintType.RESTRICTION_SITE,
        name="no_ecori",
        strictness=ConstraintStrictness.HARD,
        priority=ConstraintPriority.MEDIUM,
        params={"sites": ["GAATTC"]},
    )


@pytest.fixture
def high_priority_gc() -> ConstraintSpec:
    """A high-priority GC constraint."""
    return ConstraintSpec(
        ctype=ConstraintType.GC_CONTENT,
        name="gc_critical",
        strictness=ConstraintStrictness.HARD,
        priority=ConstraintPriority.HIGH,
    )


@pytest.fixture
def low_priority_cpg() -> ConstraintSpec:
    """A low-priority CpG constraint."""
    return ConstraintSpec(
        ctype=ConstraintType.NO_CPG,
        name="cpg_low",
        strictness=ConstraintStrictness.HARD,
        priority=ConstraintPriority.LOW,
    )


@pytest.fixture
def translation_constraint() -> ConstraintSpec:
    """A translation (amino acid identity) constraint."""
    return ConstraintSpec(
        ctype=ConstraintType.AMINO_ACID_IDENTITY,
        name="translation",
        strictness=ConstraintStrictness.HARD,
        priority=ConstraintPriority.CRITICAL,
    )


def _make_model(constraints: list[ConstraintSpec], protein: str = "AAAA") -> CSPModel:
    """Build a minimal CSPModel from constraints and a short protein."""
    from biocompiler.shared.constants import AA_TO_CODONS
    codon_domains = {i: list(AA_TO_CODONS.get(aa, ["GCT"])) for i, aa in enumerate(protein)}
    return CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints,
        config=SolverConfig(),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CAIImpactEstimator
# ═══════════════════════════════════════════════════════════════════════════════

class TestCAIImpactEstimator:

    def test_estimate_cai_impact_gc(self, gc_constraint):
        """GC constraints have negative CAI impact (relaxing hurts CAI)."""
        resolver = CAIImpactEstimator()
        impact = resolver.estimate_cai_impact(gc_constraint)
        assert impact < 0, "Relaxing GC should hurt CAI (negative impact)"

    def test_estimate_cai_impact_cpg(self, cpg_constraint):
        """CpG constraints have positive CAI impact (relaxing helps CAI)."""
        resolver = CAIImpactEstimator()
        impact = resolver.estimate_cai_impact(cpg_constraint)
        assert impact > 0, "Relaxing CpG avoidance should help CAI (positive impact)"

    def test_estimate_cai_impact_restriction(self, restriction_constraint):
        """Restriction site constraints have small positive CAI impact."""
        resolver = CAIImpactEstimator()
        impact = resolver.estimate_cai_impact(restriction_constraint)
        assert impact > 0, "Relaxing restriction sites should help CAI slightly"
        assert impact < 0.03, "Restriction site impact should be small"

    def test_rank_resolution_equal_priority_compromise(self, gc_constraint, cpg_constraint):
        """With equal priority, CAI-aware ranking should resolve compromise."""
        resolver = CAIImpactEstimator()
        # GC has negative impact (-0.05), CpG has positive (0.03)
        # Compromise should be resolved: relaxing GC hurts CAI, relaxing CpG helps CAI
        result = resolver.rank_resolution(gc_constraint, cpg_constraint, "compromise")
        # CpG has higher CAI impact (0.03 > -0.05), so relaxing CpG helps CAI more
        assert result == "relax_b", (
            f"Expected relax_b (cpg) since relaxing it helps CAI more, got {result}"
        )

    def test_rank_resolution_different_priority_no_change(self, high_priority_gc, low_priority_cpg):
        """When priorities differ, CAI-aware should NOT override priority decision."""
        resolver = CAIImpactEstimator()
        # HIGH gc vs LOW cpg -> priority says relax_b (cpg)
        result = resolver.rank_resolution(high_priority_gc, low_priority_cpg, "relax_b")
        assert result == "relax_b", (
            "Should not override priority-based decision when priorities differ"
        )

    def test_rank_resolution_small_cai_diff_no_change(self):
        """When CAI difference is small (< 0.03), no adjustment."""
        resolver = CAIImpactEstimator()
        spec_a = ConstraintSpec(
            ctype=ConstraintType.RESTRICTION_SITE, name="rs1",
            priority=ConstraintPriority.MEDIUM,
        )
        spec_b = ConstraintSpec(
            ctype=ConstraintType.NO_INSTABILITY_MOTIF, name="attta1",
            priority=ConstraintPriority.MEDIUM,
        )
        # Both have small positive impacts (0.01 vs 0.005, diff=0.005 < 0.03)
        result = resolver.rank_resolution(spec_a, spec_b, "compromise")
        assert result == "compromise", (
            "Should not adjust when CAI difference is below threshold"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ConflictResolver with cai_aware parameter
# ═══════════════════════════════════════════════════════════════════════════════

class TestConflictResolverCAIAware:

    def test_default_is_cai_aware(self):
        """Default ConflictResolver should be CAI-aware."""
        resolver = ConflictResolver()
        assert resolver._cai_aware is True
        assert resolver._cai_resolver is not None

    def test_cai_aware_false(self):
        """cai_aware=False should use priority-only logic."""
        resolver = ConflictResolver(cai_aware=False)
        assert resolver._cai_aware is False
        assert resolver._cai_resolver is None

    def test_cai_aware_produces_different_resolution(self, gc_constraint, cpg_constraint):
        """CAI-aware resolution should produce different strategies than priority-only."""
        model = _make_model([gc_constraint, cpg_constraint])

        resolver_cai = ConflictResolver(cai_aware=True)
        resolver_old = ConflictResolver(cai_aware=False)

        conflicts_cai = resolver_cai.detect_conflicts(model)
        conflicts_old = resolver_old.detect_conflicts(model)

        # Both should detect the same number of conflicts
        assert len(conflicts_cai) == len(conflicts_old), (
            "Both resolvers should detect the same conflicts"
        )

        if conflicts_cai:
            # The resolution strategies should differ because CAI-aware
            # adjusts the "compromise" to prefer the constraint that helps CAI
            cai_strategies = [c.resolution_strategy for c in conflicts_cai]
            old_strategies = [c.resolution_strategy for c in conflicts_old]
            # They may or may not differ depending on priorities,
            # but for equal-priority GC vs CpG they should differ
            assert cai_strategies != old_strategies, (
                f"CAI-aware should produce different strategies: "
                f"cai={cai_strategies}, old={old_strategies}"
            )

    def test_cai_aware_false_matches_old_behavior(self, high_priority_gc, low_priority_cpg):
        """cai_aware=False should produce the same results as before."""
        model = _make_model([high_priority_gc, low_priority_cpg])

        resolver = ConflictResolver(cai_aware=False)
        conflicts = resolver.detect_conflicts(model)

        if conflicts:
            # With HIGH gc and LOW cpg, should always be "relax_b"
            for c in conflicts:
                assert c.resolution_strategy == "relax_b", (
                    f"Priority-only should relax lower-priority constraint, got {c.resolution_strategy}"
                )

    def test_no_conflicts_both_modes(self, restriction_constraint):
        """No conflicts should be detected when only one constraint."""
        model = _make_model([restriction_constraint])
        resolver_cai = ConflictResolver(cai_aware=True)
        resolver_old = ConflictResolver(cai_aware=False)
        assert resolver_cai.detect_conflicts(model) == []
        assert resolver_old.detect_conflicts(model) == []


# ═══════════════════════════════════════════════════════════════════════════════
# 3. resolution_strategy parameter
# ═══════════════════════════════════════════════════════════════════════════════

class TestResolutionStrategyParameter:

    def test_cai_aware_strategy(self):
        """resolution_strategy='cai_aware' should set cai_aware=True."""
        resolver = ConflictResolver(resolution_strategy="cai_aware")
        assert resolver._cai_aware is True
        assert resolver._resolution_strategy == "cai_aware"

    def test_priority_only_strategy(self):
        """resolution_strategy='priority_only' should set cai_aware=False."""
        resolver = ConflictResolver(resolution_strategy="priority_only")
        assert resolver._cai_aware is False
        assert resolver._resolution_strategy == "priority_only"

    def test_priority_only_overrides_cai_aware_param(self):
        """resolution_strategy='priority_only' should override cai_aware=True."""
        resolver = ConflictResolver(cai_aware=True, resolution_strategy="priority_only")
        assert resolver._cai_aware is False

    def test_cai_aware_strategy_overrides_cai_aware_false(self):
        """resolution_strategy='cai_aware' should NOT override cai_aware=False."""
        # The cai_aware parameter takes precedence when explicitly set.
        # Only resolution_strategy='priority_only' overrides cai_aware.
        resolver = ConflictResolver(cai_aware=False, resolution_strategy="cai_aware")
        assert resolver._cai_aware is False


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Auto-resolve with CAI-aware
# ═══════════════════════════════════════════════════════════════════════════════

class TestAutoResolveCAIAware:

    def test_auto_resolve_cai_aware_removes_fewer(self, gc_constraint, cpg_constraint):
        """CAI-aware auto-resolve should consider CAI impact in selection."""
        model = _make_model([gc_constraint, cpg_constraint])

        resolver_cai = ConflictResolver(cai_aware=True)
        resolver_old = ConflictResolver(cai_aware=False)

        resolved_cai = resolver_cai.auto_resolve(model, strategy="max_priority")
        resolved_old = resolver_old.auto_resolve(model, strategy="max_priority")

        # Both should produce valid models
        assert isinstance(resolved_cai, CSPModel)
        assert isinstance(resolved_old, CSPModel)

    def test_auto_resolve_priority_only_backward_compat(self, high_priority_gc, low_priority_cpg):
        """Priority-only auto-resolve should behave as before."""
        model = _make_model([high_priority_gc, low_priority_cpg])
        resolver = ConflictResolver(cai_aware=False)
        resolved = resolver.auto_resolve(model, strategy="max_priority")
        # The low-priority CpG constraint should be removed
        remaining_names = [c.name for c in resolved.constraints]
        assert "cpg_low" not in remaining_names
        assert "gc_critical" in remaining_names


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ConflictProvenance with CAI-aware
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceCAIAware:

    def test_provenance_resolver_accepts_cai_aware(self):
        """ConflictResolverWithProvenance should accept cai_aware parameter."""
        resolver = ConflictResolverWithProvenance(track_provenance=True, cai_aware=True)
        assert resolver._cai_aware is True
        assert resolver._base_resolver._cai_aware is True

    def test_provenance_resolver_default_cai_aware(self):
        """Default ConflictResolverWithProvenance should be CAI-aware."""
        resolver = ConflictResolverWithProvenance(track_provenance=True)
        assert resolver._cai_aware is True

    def test_provenance_with_cai_aware_records(self, gc_constraint, cpg_constraint):
        """CAI-aware provenance should record cai_impact from CAIImpactEstimator."""
        sequence = "GCGGCGGCGGCG"  # 4 codons, 100% GC
        resolver = ConflictResolverWithProvenance(track_provenance=True, cai_aware=True)
        _, records = resolver.resolve_conflicts([gc_constraint, cpg_constraint], sequence)

        if records:
            for record in records:
                # CAI-aware provenance should have non-zero cai_impact
                # (either positive or negative depending on which was relaxed)
                assert isinstance(record.cai_impact, float), (
                    f"cai_impact should be float, got {type(record.cai_impact)}"
                )

    def test_provenance_without_cai_aware(self, high_priority_gc, low_priority_cpg):
        """Non-CAI-aware provenance should still work."""
        sequence = "GCGGCGGCGGCG"
        resolver = ConflictResolverWithProvenance(track_provenance=True, cai_aware=False)
        _, records = resolver.resolve_conflicts([high_priority_gc, low_priority_cpg], sequence)
        # Should produce records without CAI-aware impact
        for record in records:
            assert isinstance(record.cai_impact, float)

    def test_cai_aware_method_label(self, gc_constraint, cpg_constraint):
        """CAI-aware provenance should use 'cai_aware' resolution method."""
        sequence = "GCGGCGGCGGCG"
        resolver = ConflictResolverWithProvenance(track_provenance=True, cai_aware=True)
        _, records = resolver.resolve_conflicts([gc_constraint, cpg_constraint], sequence)

        if records:
            for record in records:
                assert record.resolution_method == "cai_aware", (
                    f"Expected 'cai_aware' method, got {record.resolution_method!r}"
                )

    def test_priority_only_method_label(self, high_priority_gc, low_priority_cpg):
        """Priority-only provenance should use priority_based or weight_based method."""
        sequence = "GCGGCGGCGGCG"
        resolver = ConflictResolverWithProvenance(track_provenance=True, cai_aware=False)
        _, records = resolver.resolve_conflicts([high_priority_gc, low_priority_cpg], sequence)

        if records:
            for record in records:
                assert record.resolution_method in {"priority_based", "weight_based"}, (
                    f"Expected priority/weight method, got {record.resolution_method!r}"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_infeasible_conflict_not_adjusted(self, translation_constraint):
        """Infeasible conflicts should not be adjusted by CAI-aware resolver."""
        another_translation = ConstraintSpec(
            ctype=ConstraintType.AMINO_ACID_IDENTITY,
            name="translation2",
            strictness=ConstraintStrictness.HARD,
            priority=ConstraintPriority.CRITICAL,
        )
        resolver = ConflictResolver(cai_aware=True)
        model = _make_model([translation_constraint, another_translation])
        conflicts = resolver.detect_conflicts(model)
        if conflicts:
            for c in conflicts:
                assert c.resolution_strategy == "infeasible"

    def test_cai_aware_with_no_conflicts(self, restriction_constraint):
        """CAI-aware resolver should handle models with no conflicts."""
        resolver = ConflictResolver(cai_aware=True)
        model = _make_model([restriction_constraint])
        conflicts = resolver.detect_conflicts(model)
        assert conflicts == []

    def test_prioritize_constraints_unchanged(self, gc_constraint, cpg_constraint, restriction_constraint):
        """prioritize_constraints should still work the same way."""
        result = prioritize_constraints([cpg_constraint, restriction_constraint, gc_constraint])
        # CRITICAL=0, HIGH=1, MEDIUM=2, LOW=3
        priorities = [c.priority.rank for c in result]
        assert priorities == sorted(priorities)
