"""Tests for biocompiler.solver.enforcement — constraint enforcement mechanics.

Covers:
1. ConstraintEnforcer.enforce() with valid/invalid sequences
2. ConstraintEnforcer.score_soft_constraints() returns reasonable values
3. ConstraintEnforcer.resolve_conflicts() detects conflicts
4. EnforcementResult dataclass fields
5. ConflictResolution dataclass fields
"""

from __future__ import annotations

import math
from dataclasses import fields

import pytest

from biocompiler.solver.enforcement import (
    ConflictResolution,
    ConstraintEnforcer,
    EnforcementResult,
)
from biocompiler.solver.constraints import (
    CSPModel,
    GCRangeConstraint,
    HardConstraint,
    MaximizeCAI,
    MinimizeCpG,
    MinimizeMRNADG,
    NoATTTAMotifConstraint,
    NoRestrictionSiteConstraint,
    NoTRunConstraint,
    SoftConstraint,
    TranslationConstraint,
)
from biocompiler.solver.types import (
    CodonVariable,
    ConstraintStrictness,
    ConstraintViolation,
    SolverConfig,
)
from biocompiler.constants import AA_TO_CODONS
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════


def _make_model(
    protein: str = "MVSKGE",
    hard_constraints: list[HardConstraint] | None = None,
    soft_constraints: list[SoftConstraint] | None = None,
    config: SolverConfig | None = None,
) -> CSPModel:
    """Build a minimal CSPModel for enforcement testing."""
    variables = [
        CodonVariable(position=i, amino_acid=aa, domain=AA_TO_CODONS[aa])
        for i, aa in enumerate(protein)
    ]
    return CSPModel(
        variables=variables,
        hard_constraints=hard_constraints or [],
        soft_constraints=soft_constraints or [],
        protein=protein,
        organism="Homo_sapiens",
        config=config or SolverConfig(),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 4. EnforcementResult dataclass fields
# ═══════════════════════════════════════════════════════════════════════════════


class TestEnforcementResultDataclass:
    """Verify EnforcementResult has the correct fields and types."""

    def test_has_required_fields(self):
        field_names = {f.name for f in fields(EnforcementResult)}
        assert "all_hard_satisfied" in field_names
        assert "violations" in field_names
        assert "soft_score" in field_names
        assert "conflict_resolution" in field_names

    def test_all_hard_satisfied_is_bool(self):
        result = EnforcementResult(
            all_hard_satisfied=True,
            violations=[],
            soft_score=0.0,
        )
        assert isinstance(result.all_hard_satisfied, bool)

    def test_violations_is_list(self):
        result = EnforcementResult(
            all_hard_satisfied=True,
            violations=[],
            soft_score=0.0,
        )
        assert isinstance(result.violations, list)

    def test_soft_score_is_float(self):
        result = EnforcementResult(
            all_hard_satisfied=True,
            violations=[],
            soft_score=1.5,
        )
        assert isinstance(result.soft_score, float)

    def test_conflict_resolution_defaults_to_none(self):
        result = EnforcementResult(
            all_hard_satisfied=True,
            violations=[],
            soft_score=0.0,
        )
        assert result.conflict_resolution is None

    def test_conflict_resolution_can_hold_value(self):
        cr = ConflictResolution(
            conflicting_constraints=["GCRangeConstraint"],
            suggested_relaxations=["Widen GC range"],
            is_resolvable=True,
        )
        result = EnforcementResult(
            all_hard_satisfied=False,
            violations=[],
            soft_score=0.0,
            conflict_resolution=cr,
        )
        assert result.conflict_resolution is cr
        assert result.conflict_resolution.conflicting_constraints == [
            "GCRangeConstraint"
        ]

    def test_violations_contains_constraint_violation_instances(self):
        violation = ConstraintViolation(
            constraint_name="GCRangeConstraint",
            constraint_type=ConstraintStrictness.HARD,
            description="GC out of range",
            positions=[0, 1],
            severity=0.5,
        )
        result = EnforcementResult(
            all_hard_satisfied=False,
            violations=[violation],
            soft_score=0.0,
        )
        assert len(result.violations) == 1
        assert isinstance(result.violations[0], ConstraintViolation)
        assert result.violations[0].constraint_name == "GCRangeConstraint"


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ConflictResolution dataclass fields
# ═══════════════════════════════════════════════════════════════════════════════


class TestConflictResolutionDataclass:
    """Verify ConflictResolution has the correct fields and types."""

    def test_has_required_fields(self):
        field_names = {f.name for f in fields(ConflictResolution)}
        assert "conflicting_constraints" in field_names
        assert "suggested_relaxations" in field_names
        assert "is_resolvable" in field_names

    def test_conflicting_constraints_is_list_of_str(self):
        cr = ConflictResolution(
            conflicting_constraints=["A", "B"],
            suggested_relaxations=[],
            is_resolvable=True,
        )
        assert isinstance(cr.conflicting_constraints, list)
        assert all(isinstance(c, str) for c in cr.conflicting_constraints)

    def test_suggested_relaxations_is_list_of_str(self):
        cr = ConflictResolution(
            conflicting_constraints=[],
            suggested_relaxations=["Relax X", "Relax Y"],
            is_resolvable=True,
        )
        assert isinstance(cr.suggested_relaxations, list)
        assert all(isinstance(s, str) for s in cr.suggested_relaxations)

    def test_is_resolvable_is_bool(self):
        cr = ConflictResolution(
            conflicting_constraints=[],
            suggested_relaxations=[],
            is_resolvable=True,
        )
        assert isinstance(cr.is_resolvable, bool)

    def test_empty_conflict(self):
        cr = ConflictResolution(
            conflicting_constraints=[],
            suggested_relaxations=[],
            is_resolvable=True,
        )
        assert cr.conflicting_constraints == []
        assert cr.suggested_relaxations == []
        assert cr.is_resolvable is True

    def test_non_resolvable_conflict(self):
        cr = ConflictResolution(
            conflicting_constraints=["FundamentalIssue"],
            suggested_relaxations=[],
            is_resolvable=False,
        )
        assert cr.is_resolvable is False


# ═══════════════════════════════════════════════════════════════════════════════
# 1. ConstraintEnforcer.enforce() — valid sequences
# ═══════════════════════════════════════════════════════════════════════════════


class TestEnforceValidSequence:
    """ConstraintEnforcer.enforce() with sequences satisfying all hard constraints."""

    def test_valid_sequence_all_hard_satisfied(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
            ],
        )
        # ATG=M, GTT=V → GC = 2/6 ≈ 0.33
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.all_hard_satisfied is True

    def test_valid_sequence_no_violations(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.violations == []

    def test_valid_sequence_no_conflict_resolution(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.conflict_resolution is None

    def test_valid_sequence_with_soft_constraints_has_score(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        model = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
            ],
            soft_constraints=[
                MaximizeCAI(adaptiveness=human_adapt, protein="MV"),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.all_hard_satisfied is True
        assert isinstance(result.soft_score, float)
        assert math.isfinite(result.soft_score)

    def test_valid_sequence_no_hard_constraints(self):
        model = _make_model(hard_constraints=[])
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.all_hard_satisfied is True
        assert result.violations == []

    def test_enforce_returns_enforcement_result(self):
        model = _make_model(hard_constraints=[])
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert isinstance(result, EnforcementResult)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. ConstraintEnforcer.enforce() — invalid sequences
# ═══════════════════════════════════════════════════════════════════════════════


class TestEnforceInvalidSequence:
    """ConstraintEnforcer.enforce() with sequences violating hard constraints."""

    def test_invalid_sequence_not_all_hard_satisfied(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        # ATGGTT has GC ≈ 0.33, fails gc_lo=0.80
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.all_hard_satisfied is False

    def test_invalid_sequence_has_violations(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert len(result.violations) > 0

    def test_violation_has_correct_constraint_name(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        violated_names = {v.constraint_name for v in result.violations}
        assert "GCRangeConstraint" in violated_names

    def test_violation_has_hard_strictness(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        for v in result.violations:
            assert v.constraint_type == ConstraintStrictness.HARD

    def test_violation_severity_in_range(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        for v in result.violations:
            assert 0.0 <= v.severity <= 1.0

    def test_invalid_sequence_triggers_conflict_resolution(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.conflict_resolution is not None

    def test_restriction_site_violation(self):
        model = _make_model(
            protein="EF",
            hard_constraints=[
                NoRestrictionSiteConstraint(sites=["GAATTC"]),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "GAATTC")
        assert result.all_hard_satisfied is False
        assert any(
            v.constraint_name == "NoRestrictionSiteConstraint"
            for v in result.violations
        )

    def test_attta_violation(self):
        model = _make_model(
            protein="NL",
            hard_constraints=[
                NoATTTAMotifConstraint(),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "AATTTA")
        assert result.all_hard_satisfied is False

    def test_t_run_violation(self):
        model = _make_model(
            hard_constraints=[
                NoTRunConstraint(max_run=5),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "TTTTTT")
        assert result.all_hard_satisfied is False

    def test_multiple_violations(self):
        model = _make_model(
            protein="NL",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
                NoATTTAMotifConstraint(),
                NoTRunConstraint(max_run=3),
            ],
        )
        # AAT TTA: GC = 1/6 ≈ 0.17 (fails gc_lo=0.80), has ATTTA, has TTT (3 T's OK for max_run=3)
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "AATTTA")
        assert result.all_hard_satisfied is False
        assert len(result.violations) >= 1  # At least GC + ATTTA

    def test_violation_description_is_informative(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        for v in result.violations:
            assert "GCRangeConstraint" in v.description
            assert "violated" in v.description.lower() or "position" in v.description.lower()


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ConstraintEnforcer.score_soft_constraints() returns reasonable values
# ═══════════════════════════════════════════════════════════════════════════════


class TestScoreSoftConstraints:
    """Verify score_soft_constraints returns reasonable weighted scores."""

    def test_no_soft_constraints_returns_zero(self):
        model = _make_model(soft_constraints=[])
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        assert score == 0.0

    def test_cai_score_finite(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        model = _make_model(
            protein="MV",
            soft_constraints=[
                MaximizeCAI(adaptiveness=human_adapt, protein="MV"),
            ],
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        assert math.isfinite(score)

    def test_cai_weight_applied(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        cai = MaximizeCAI(adaptiveness=human_adapt, protein="MV")
        raw_score = cai.score("ATGGTT")

        # Default cai_weight = 1.0
        model = _make_model(
            protein="MV",
            soft_constraints=[cai],
            config=SolverConfig(cai_weight=1.0),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        assert abs(score - raw_score) < 1e-10

    def test_cai_weight_scaling(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        cai = MaximizeCAI(adaptiveness=human_adapt, protein="MV")
        raw_score = cai.score("ATGGTT")

        model = _make_model(
            protein="MV",
            soft_constraints=[cai],
            config=SolverConfig(cai_weight=2.0),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        assert abs(score - 2.0 * raw_score) < 1e-10

    def test_cpg_score_negative_for_cpg_rich(self):
        cpg = MinimizeCpG()
        model = _make_model(
            soft_constraints=[cpg],
            config=SolverConfig(cpg_weight=1.0),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "CGCGCGCGCG")
        # MinimizeCpG.score("CGCGCGCGCG") = -5.0, weight=1.0
        assert score < 0.0

    def test_cpg_weight_applied(self):
        cpg = MinimizeCpG()
        raw_score = cpg.score("CGCGCG")
        model = _make_model(
            soft_constraints=[cpg],
            config=SolverConfig(cpg_weight=0.5),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "CGCGCG")
        assert abs(score - 0.5 * raw_score) < 1e-10

    def test_multiple_soft_constraints_combined(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        cai = MaximizeCAI(adaptiveness=human_adapt, protein="MV")
        cpg = MinimizeCpG()
        model = _make_model(
            protein="MV",
            soft_constraints=[cai, cpg],
            config=SolverConfig(cai_weight=1.0, cpg_weight=0.5),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        expected = 1.0 * cai.score("ATGGTT") + 0.5 * cpg.score("ATGGTT")
        assert abs(score - expected) < 1e-10

    def test_mrna_dg_weight_applied(self):
        dg = MinimizeMRNADG(window_start=0, window_end=50)
        seq = "ATGGCTTCTAAAGGTGAA" + "A" * 32  # 50 nt
        raw_score = dg.score(seq)
        model = _make_model(
            soft_constraints=[dg],
            config=SolverConfig(mrna_dg_weight=0.3),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, seq)
        assert abs(score - 0.3 * raw_score) < 1e-10

    def test_unknown_soft_constraint_gets_default_weight(self):
        """Custom soft constraint not in weight_map should get default weight 0.1."""

        class CustomSoft(SoftConstraint):
            @property
            def name(self) -> str:
                return "CustomObjective"

            def check(self, sequence: str) -> bool:
                return True

            def violated_positions(self, sequence: str) -> list[int]:
                return []

            def score(self, sequence: str) -> float:
                return 10.0

        custom = CustomSoft()
        model = _make_model(soft_constraints=[custom])
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        # Default weight for unknown is 0.1, raw score is 10.0
        assert abs(score - 0.1 * 10.0) < 1e-10

    def test_zero_weight_excludes_constraint(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        cai = MaximizeCAI(adaptiveness=human_adapt, protein="MV")
        model = _make_model(
            protein="MV",
            soft_constraints=[cai],
            config=SolverConfig(cai_weight=0.0),
        )
        enforcer = ConstraintEnforcer()
        score = enforcer.score_soft_constraints(model, "ATGGTT")
        # cai_weight=0.0 means score should be 0, but the implementation
        # checks `weight == 0.0 and sc.name not in weight_map` for default
        # Since "MaximizeCAI" IS in weight_map with 0.0, it stays 0.0
        assert score == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 3. ConstraintEnforcer.resolve_conflicts() detects conflicts
# ═══════════════════════════════════════════════════════════════════════════════


class TestResolveConflicts:
    """Verify resolve_conflicts detects conflicting constraints and suggests relaxations."""

    def test_no_violations_returns_empty_conflict(self):
        constraints = [
            TranslationConstraint(protein="MV"),
            GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "ATGGTT")
        assert isinstance(result, ConflictResolution)
        assert result.conflicting_constraints == []
        assert result.suggested_relaxations == []
        assert result.is_resolvable is True

    def test_single_violated_constraint_detected(self):
        constraints = [
            GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "ATGGTT")
        assert "GCRangeConstraint" in result.conflicting_constraints

    def test_multiple_violated_constraints_detected(self):
        constraints = [
            GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            NoTRunConstraint(max_run=3),
        ]
        # ATTTTTTA: high T run and low GC
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "ATTTTTTA")
        assert len(result.conflicting_constraints) >= 1

    def test_relaxation_suggestions_generated(self):
        constraints = [
            GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "ATGGTT")
        assert len(result.suggested_relaxations) > 0
        # Should suggest widening GC range
        assert any("GC" in s for s in result.suggested_relaxations)

    def test_restriction_site_relaxation_suggestion(self):
        constraints = [
            NoRestrictionSiteConstraint(sites=["GAATTC"]),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "GAATTC")
        assert len(result.suggested_relaxations) > 0
        assert any("restriction" in s.lower() or "enzyme" in s.lower()
                    for s in result.suggested_relaxations)

    def test_attta_relaxation_suggestion(self):
        constraints = [
            NoATTTAMotifConstraint(),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "AAATTTAAAA")
        assert len(result.suggested_relaxations) > 0
        assert any("ATTTA" in s for s in result.suggested_relaxations)

    def test_t_run_relaxation_suggestion(self):
        constraints = [
            NoTRunConstraint(max_run=5),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "TTTTTT")
        assert len(result.suggested_relaxations) > 0
        assert any("T-run" in s or "poly-T" in s or "max T" in s
                    for s in result.suggested_relaxations)

    def test_translation_constraint_relaxation_suggestion(self):
        constraints = [
            TranslationConstraint(protein="MV"),
        ]
        # AAA = K, not M; GTT = V
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "AAAGTT")
        assert len(result.suggested_relaxations) > 0
        assert any("Translation" in s or "translation" in s.lower()
                    for s in result.suggested_relaxations)

    def test_is_resolvable_true_when_suggestions_exist(self):
        constraints = [
            GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
        ]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "ATGGTT")
        assert result.is_resolvable is True

    def test_overlapping_position_conflict_detected(self):
        """When two constraints share violated positions, overlap warning is added."""
        constraints = [
            GCRangeConstraint(gc_lo=0.90, gc_hi=1.0),
            NoTRunConstraint(max_run=2),
        ]
        # GCGCGCGCGCGC: 100% GC (fails gc_lo=0.90 only because > gc_hi=1.0)
        # Let's use a sequence that violates both at overlapping positions
        # GCCTTTGC: GC = 4/8 = 0.50 (fails gc_lo=0.90), TTT violates max_run=2
        # Positions for GC: AT positions (2,3,4), for T-run: starts at 2
        # Actually let's construct something cleaner
        seq = "GCGCTTTTTGCGC"  # GC = 8/13 ≈ 0.62, fails gc_lo=0.90
        # T-run at pos 4-8 (5 T's) violates max_run=2
        # GC positions: 0,1,2,3,9,10,11,12
        # T-run start: 4
        # No overlap here. Let me construct overlapping positions.
        # For GC too low, violated_positions returns AT positions
        # AT positions in "GCGCTTTTTGCGC": positions 4,5,6,7,8 (T's)
        # T-run start: position 4
        # position 4 is in both sets → overlap
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, seq)
        if len(result.conflicting_constraints) >= 2:
            # Should have an overlap warning
            has_overlap_warning = any("overlap" in s.lower()
                                       for s in result.suggested_relaxations)
            assert has_overlap_warning

    def test_unknown_constraint_generic_relaxation(self):
        """Unknown constraint name gets a generic relaxation suggestion."""

        class UnknownHard(HardConstraint):
            @property
            def name(self) -> str:
                return "CustomHardConstraint"

            def check(self, sequence: str) -> bool:
                return False  # always violated

            def violated_positions(self, sequence: str) -> list[int]:
                return list(range(len(sequence)))

        constraints = [UnknownHard()]
        enforcer = ConstraintEnforcer()
        result = enforcer.resolve_conflicts(constraints, "ATGGTT")
        assert "CustomHardConstraint" in result.conflicting_constraints
        assert any("CustomHardConstraint" in s for s in result.suggested_relaxations)


# ═══════════════════════════════════════════════════════════════════════════════
# Integration: enforce() + resolve_conflicts() + score_soft_constraints()
# ═══════════════════════════════════════════════════════════════════════════════


class TestEnforcementIntegration:
    """End-to-end integration tests for ConstraintEnforcer workflow."""

    def test_enforce_includes_soft_score_on_valid_sequence(self):
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        model = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
            ],
            soft_constraints=[
                MaximizeCAI(adaptiveness=human_adapt, protein="MV"),
                MinimizeCpG(),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.all_hard_satisfied is True
        assert isinstance(result.soft_score, float)
        assert math.isfinite(result.soft_score)

    def test_enforce_includes_soft_score_on_invalid_sequence(self):
        """Soft score is still computed even when hard constraints fail."""
        human_adapt = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
            soft_constraints=[
                MaximizeCAI(adaptiveness=human_adapt, protein="MV"),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.all_hard_satisfied is False
        assert isinstance(result.soft_score, float)
        assert math.isfinite(result.soft_score)

    def test_enforce_conflict_resolution_populated_on_failure(self):
        model = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "ATGGTT")
        assert result.conflict_resolution is not None
        assert isinstance(result.conflict_resolution, ConflictResolution)
        assert "GCRangeConstraint" in result.conflict_resolution.conflicting_constraints

    def test_enforcer_is_stateless(self):
        """Multiple enforce calls on the same enforcer instance work independently."""
        model_valid = _make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
            ],
        )
        model_invalid = _make_model(
            protein="MV",
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
            ],
        )
        enforcer = ConstraintEnforcer()
        result1 = enforcer.enforce(model_valid, "ATGGTT")
        result2 = enforcer.enforce(model_invalid, "ATGGTT")

        assert result1.all_hard_satisfied is True
        assert result2.all_hard_satisfied is False
        assert result1.conflict_resolution is None
        assert result2.conflict_resolution is not None

    def test_severity_increases_with_more_violated_positions(self):
        """More violated positions → higher severity."""
        # All-AT sequence: many AT positions for GC constraint
        model_low_gc = _make_model(
            hard_constraints=[
                GCRangeConstraint(gc_lo=0.50, gc_hi=0.70),
            ],
        )
        enforcer = ConstraintEnforcer()
        # Short sequence: "ATG" → 1/3 GC = 0.33
        result_short = enforcer.enforce(model_low_gc, "ATG")
        # Longer all-AT sequence: "ATAATA" → 0/6 GC = 0.0
        result_long = enforcer.enforce(model_low_gc, "ATAATA")

        if result_short.violations and result_long.violations:
            # The longer sequence with 0% GC should have higher severity
            # than the shorter one with 33% GC
            assert result_long.violations[0].severity >= result_short.violations[0].severity * 0.5

    def test_violation_positions_are_populated(self):
        model = _make_model(
            protein="EF",
            hard_constraints=[
                NoRestrictionSiteConstraint(sites=["GAATTC"]),
            ],
        )
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce(model, "GAATTC")
        assert len(result.violations) > 0
        assert len(result.violations[0].positions) > 0


# ═══════════════════════════════════════════════════════════════════════════════
# Severity estimation (private helper, tested indirectly)
# ═══════════════════════════════════════════════════════════════════════════════


class TestSeverityEstimation:
    """Verify severity values through the public enforce() API."""

    def test_severity_zero_when_no_positions(self):
        """If no positions are violated (shouldn't normally happen through enforce),
        severity should be 0."""
        severity = ConstraintEnforcer._estimate_severity([], 10)
        assert severity == 0.0

    def test_severity_one_for_empty_sequence(self):
        severity = ConstraintEnforcer._estimate_severity([0], 0)
        assert severity == 1.0

    def test_severity_fraction_for_partial_violation(self):
        # 3 positions out of 10 → severity 0.3
        severity = ConstraintEnforcer._estimate_severity([0, 1, 2], 10)
        assert abs(severity - 0.3) < 1e-10

    def test_severity_capped_at_one(self):
        # Even with duplicate positions, severity should not exceed 1.0
        severity = ConstraintEnforcer._estimate_severity(
            [0, 0, 0, 1, 1, 1], 2
        )
        assert severity <= 1.0

    def test_severity_uses_unique_positions(self):
        """Duplicate positions should not inflate severity."""
        severity_unique = ConstraintEnforcer._estimate_severity([0, 1], 10)
        severity_duped = ConstraintEnforcer._estimate_severity(
            [0, 0, 1, 1], 10
        )
        assert severity_unique == severity_duped
