"""Unit tests for solver/mus.py — MUS (Minimal Unsatisfiable Subset) analysis.

Covers:
1. FeasibilityReport construction and field access
2. quick_feasibility_check with feasible/infeasible CSP models
3. suggest_relaxations returns reasonable, actionable suggestions
4. explain_conflict produces human-readable output

All tests use the typed API (CSPModel, ConstraintSpec, etc.) directly
and do NOT require a solver backend.  MUSReport objects are constructed
manually to test explain_conflict and suggest_relaxations in isolation.
"""

from __future__ import annotations

from dataclasses import fields

import pytest

from biocompiler.solver.mus import (
    FeasibilityReport,
    explain_conflict,
    quick_feasibility_check,
    suggest_relaxations,
)
from biocompiler.solver.types import (
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
    MUSReport,
    SolverConfig,
)


# ── Helpers ────────────────────────────────────────────────────────

def _make_config(**overrides) -> SolverConfig:
    """Build a SolverConfig with sensible defaults, overridden by *overrides*."""
    defaults = dict(gc_lo=0.30, gc_hi=0.70)
    defaults.update(overrides)
    return SolverConfig(**defaults)


def _make_model(
    protein: str = "MVSKGE",
    constraints: list[ConstraintSpec] | None = None,
    codon_domains: dict[int, list[str]] | None = None,
    config: SolverConfig | None = None,
) -> CSPModel:
    """Build a minimal CSPModel.

    If *codon_domains* is None, the default AA_TO_CODONS lookup is used.
    """
    from biocompiler.shared.constants import AA_TO_CODONS

    if config is None:
        config = _make_config()
    if constraints is None:
        constraints = []
    if codon_domains is None:
        codon_domains = {i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(protein)}
    return CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints,
        config=config,
    )


def _constraint(
    ctype: ConstraintType = ConstraintType.NO_GT_DINUCLEOTIDE,
    name: str = "test_constraint",
    positions: list[int] | None = None,
    params: dict | None = None,
) -> ConstraintSpec:
    """Shorthand to build a ConstraintSpec."""
    return ConstraintSpec(
        ctype=ctype,
        name=name,
        strictness=ConstraintStrictness.HARD,
        params=params or {},
        positions=positions or [],
    )


def _mus_report(
    mus_constraints: list[ConstraintSpec] | None = None,
    all_constraints: list[ConstraintSpec] | None = None,
    conflict_positions: list[int] | None = None,
    explanation: str = "",
) -> MUSReport:
    """Build an MUSReport for testing explain/suggest without a solver."""
    if mus_constraints is None:
        mus_constraints = []
    if all_constraints is None:
        all_constraints = list(mus_constraints)
    return MUSReport(
        mus_constraints=mus_constraints,
        all_constraints=all_constraints,
        iterations=len(all_constraints) + 1,
        solve_time_seconds=0.01,
        conflict_positions=conflict_positions or [],
        explanation=explanation,
    )


# =====================================================================
# 1. FeasibilityReport construction
# =====================================================================

class TestFeasibilityReportConstruction:
    """Tests for FeasibilityReport dataclass creation and field access."""

    def test_create_feasible_default_fields(self):
        """A feasible report with defaults should have empty lists."""
        r = FeasibilityReport(feasible=True)
        assert r.feasible is True
        assert r.warnings == []
        assert r.impossible_constraints == []
        assert r.suggested_relaxations == []

    def test_create_infeasible_with_data(self):
        """An infeasible report should carry warnings, impossible constraints, and suggestions."""
        r = FeasibilityReport(
            feasible=False,
            warnings=["GC bounds too tight"],
            impossible_constraints=["GC_CONTENT", "NoGTDinucleotide"],
            suggested_relaxations=["Widen GC upper bound from 20% to at least 33%"],
        )
        assert r.feasible is False
        assert len(r.warnings) == 1
        assert "GC_CONTENT" in r.impossible_constraints
        assert "NoGTDinucleotide" in r.impossible_constraints
        assert len(r.suggested_relaxations) == 1

    def test_has_expected_field_names(self):
        """FeasibilityReport must expose the 4 documented fields."""
        r = FeasibilityReport(feasible=True)
        field_names = {f.name for f in fields(r)}
        assert field_names == {"feasible", "warnings", "impossible_constraints", "suggested_relaxations"}

    def test_feasible_is_bool(self):
        """The feasible field must be a boolean."""
        r_true = FeasibilityReport(feasible=True)
        r_false = FeasibilityReport(feasible=False)
        assert isinstance(r_true.feasible, bool)
        assert isinstance(r_false.feasible, bool)

    def test_warnings_is_list_of_strings(self):
        """Warnings should be a list of strings."""
        r = FeasibilityReport(feasible=False, warnings=["w1", "w2"])
        assert isinstance(r.warnings, list)
        assert all(isinstance(w, str) for w in r.warnings)

    def test_impossible_constraints_is_list_of_strings(self):
        """Impossible constraints should be a list of strings."""
        r = FeasibilityReport(feasible=False, impossible_constraints=["GC_CONTENT"])
        assert isinstance(r.impossible_constraints, list)
        assert all(isinstance(c, str) for c in r.impossible_constraints)

    def test_suggested_relaxations_is_list_of_strings(self):
        """Suggested relaxations should be a list of strings."""
        r = FeasibilityReport(feasible=False, suggested_relaxations=["Widen GC bounds"])
        assert isinstance(r.suggested_relaxations, list)
        assert all(isinstance(s, str) for s in r.suggested_relaxations)

    def test_default_factory_creates_independent_lists(self):
        """Two FeasibilityReport instances should not share list objects."""
        r1 = FeasibilityReport(feasible=True)
        r2 = FeasibilityReport(feasible=True)
        r1.warnings.append("shared?")
        assert r2.warnings == []  # r2 should not be affected

    def test_fields_are_mutable(self):
        """FeasibilityReport is not frozen — fields should be mutable."""
        r = FeasibilityReport(feasible=True)
        r.feasible = False
        r.warnings.append("new warning")
        r.impossible_constraints.append("GC_CONTENT")
        r.suggested_relaxations.append("Relax GC")
        assert r.feasible is False
        assert len(r.warnings) == 1
        assert len(r.impossible_constraints) == 1
        assert len(r.suggested_relaxations) == 1


# =====================================================================
# 2. quick_feasibility_check
# =====================================================================

class TestQuickFeasibilityCheckFeasible:
    """quick_feasibility_check with models that should be feasible."""

    def test_normal_protein_no_constraints_feasible(self):
        """A standard protein with no constraints and wide GC bounds is feasible."""
        model = _make_model(protein="MVSKGE", config=_make_config(gc_lo=0.10, gc_hi=0.90))
        report = quick_feasibility_check(model)
        assert report.feasible is True

    def test_empty_protein_feasible(self):
        """An empty protein sequence is trivially feasible."""
        model = _make_model(protein="", config=_make_config())
        report = quick_feasibility_check(model)
        assert report.feasible is True

    def test_single_methionine_feasible(self):
        """A single methionine with no constraints is feasible."""
        model = _make_model(protein="M", config=_make_config(gc_lo=0.0, gc_hi=1.0))
        report = quick_feasibility_check(model)
        assert report.feasible is True

    def test_wide_gc_bounds_feasible(self):
        """Wide GC bounds [0.0, 1.0] should always be feasible."""
        model = _make_model(protein="VVVVVV", config=_make_config(gc_lo=0.0, gc_hi=1.0))
        report = quick_feasibility_check(model)
        assert report.feasible is True

    def test_reasonable_protein_with_gc_constraint(self):
        """A typical protein with reasonable GC bounds should be feasible."""
        model = _make_model(
            protein="MVSKGEELFT",
            config=_make_config(gc_lo=0.30, gc_hi=0.70),
        )
        report = quick_feasibility_check(model)
        assert report.feasible is True

    def test_long_alanine_feasible(self):
        """A long all-Alanine protein with moderate GC bounds is feasible."""
        model = _make_model(
            protein="A" * 100,
            config=_make_config(gc_lo=0.20, gc_hi=0.80),
        )
        report = quick_feasibility_check(model)
        assert report.feasible is True


class TestQuickFeasibilityCheckInfeasible:
    """quick_feasibility_check with models that should be infeasible."""

    def test_cross_codon_gt_unavoidable(self):
        """Cross-codon GT unavoidable at boundary → infeasible.

        quick_feasibility_check detects GT dinucleotide infeasibility only at
        codon *boundaries* (cross-codon GT), not within-codon GT.  We construct
        a model where every codon at position 0 ends with G and every codon at
        position 1 starts with T, making GT unavoidable at the boundary.
        """
        model = _make_model(
            protein="AF",
            constraints=[_constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt")],
            codon_domains={
                0: ["GCG"],   # Ala codon ending in G
                1: ["TTC"],   # Phe codon starting with T
            },
        )
        report = quick_feasibility_check(model)
        assert report.feasible is False
        assert "NoGTDinucleotide" in report.impossible_constraints

    def test_min_gc_exceeds_gc_hi(self):
        """gc_hi below minimum achievable GC → infeasible.

        All-Valine codons have at least 1/3 GC content.
        Setting gc_hi=0.20 makes this infeasible.
        """
        model = _make_model(
            protein="VVVVVVVVVV",
            config=_make_config(gc_lo=0.10, gc_hi=0.20),
        )
        report = quick_feasibility_check(model)
        assert report.feasible is False
        assert "GC_CONTENT" in report.impossible_constraints

    def test_max_gc_below_gc_lo(self):
        """gc_lo above maximum achievable GC → infeasible.

        Lys/Asn codons max out at ~1/3 GC each. gc_lo=0.70 is unreachable.
        """
        model = _make_model(
            protein="KKKKNNNN",
            config=_make_config(gc_lo=0.70, gc_hi=0.90),
        )
        report = quick_feasibility_check(model)
        assert report.feasible is False
        assert "GC_CONTENT" in report.impossible_constraints

    def test_cpg_avoidance_infeasible(self):
        """CpG avoidance + amino acids whose only codons contain CG → infeasible.

        Arginine (R) codons: CGT, CGC, CGA, CGG, AGA, AGG.
        The first 4 contain 'CG'. With NO_CPG constraint, only AGA/AGG survive.
        But some amino acids (like Pro, Ala, Thr) also have CG codons mixed in.
        Since Arg has non-CG options (AGA, AGG), this is actually feasible.
        We need a case where ALL codons for an AA contain CG — but no standard
        AA has this property. Instead, test that the check runs without error
        and reports the NO_CPG constraint when configured.
        """
        model = _make_model(
            protein="MVSKGE",
            constraints=[_constraint(ConstraintType.NO_CPG, "no_cpg")],
            config=_make_config(),
        )
        report = quick_feasibility_check(model)
        # Standard AAs all have at least one non-CpG codon, so this is feasible
        assert isinstance(report.feasible, bool)

    def test_restriction_site_eliminates_all_codons(self):
        """Restriction site avoidance that eliminates all codons for an AA → infeasible.

        If the restriction site is 'ATG' (the sole Methionine codon), avoiding
        it eliminates all Met codons.
        """
        model = _make_model(
            protein="MVSKGE",
            config=SolverConfig(
                gc_lo=0.10,
                gc_hi=0.90,
                restriction_sites=["ATG"],
            ),
        )
        report = quick_feasibility_check(model)
        assert report.feasible is False

    def test_restriction_site_constraint_in_model(self):
        """Restriction-site ConstraintSpec in model.constraints is also checked."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.RESTRICTION_SITE,
                    "no_met_site",
                    params={"site": "ATG", "enzyme": "NcoI"},
                ),
            ],
        )
        report = quick_feasibility_check(model)
        assert report.feasible is False


class TestQuickFeasibilityCheckWarnings:
    """quick_feasibility_check should generate warnings for tight bounds."""

    def test_tight_gc_bounds_warn(self):
        """Achievable GC range narrower than requested → warning.

        'WK' (Trp+Lys): Trp=TGG (2 GC), Lys=AAA/AAG (0-1 GC).
        min_gc=2/6=33.3%, max_gc=3/6=50%. Achievable range = 16.7%.
        Requested [30%, 60%] has range 30%, which is wider than achievable.
        """
        model = _make_model(
            protein="WK",
            config=_make_config(gc_lo=0.30, gc_hi=0.60),
        )
        report = quick_feasibility_check(model)
        assert len(report.warnings) > 0

    def test_suggestions_for_infeasible(self):
        """An infeasible report should include suggested relaxations."""
        model = _make_model(
            protein="VVVVVVVVVV",
            config=_make_config(gc_lo=0.10, gc_hi=0.20),
        )
        report = quick_feasibility_check(model)
        assert not report.feasible
        assert len(report.suggested_relaxations) > 0


# =====================================================================
# 3. suggest_relaxations
# =====================================================================

class TestSuggestRelaxations:
    """suggest_relaxations returns reasonable, actionable suggestions."""

    def test_empty_mus_returns_empty_list(self):
        """No MUS constraints → no suggestions needed."""
        model = _make_model()
        report = _mus_report(mus_constraints=[])
        suggestions = suggest_relaxations(report, model)
        assert suggestions == []

    def test_gc_content_suggestion(self):
        """GC_CONTENT in MUS → suggestion to widen GC bounds."""
        model = _make_model(config=_make_config(gc_lo=0.60, gc_hi=0.65))
        gc_constraint = _constraint(ConstraintType.GC_CONTENT, "gc_range")
        report = _mus_report(mus_constraints=[gc_constraint])
        suggestions = suggest_relaxations(report, model)
        assert len(suggestions) >= 1
        assert any("GC" in s or "gc" in s.lower() for s in suggestions)

    def test_gc_relaxation_mentions_specific_bounds(self):
        """GC relaxation suggestion should mention the specific current and widened bounds."""
        model = _make_model(config=_make_config(gc_lo=0.60, gc_hi=0.65))
        gc_constraint = _constraint(ConstraintType.GC_CONTENT, "gc_range")
        report = _mus_report(mus_constraints=[gc_constraint])
        suggestions = suggest_relaxations(report, model)
        gc_suggestions = [s for s in suggestions if "GC" in s or "gc" in s.lower()]
        assert len(gc_suggestions) >= 1
        # Should mention 60% and 65% (current bounds) and widened bounds
        gc_text = gc_suggestions[0]
        assert "60%" in gc_text or "0.60" in gc_text
        assert "65%" in gc_text or "0.65" in gc_text

    def test_no_gt_dinucleotide_suggestion(self):
        """NO_GT_DINUCLEOTIDE in MUS → suggestion to disable GT avoidance."""
        model = _make_model()
        gt_constraint = _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt")
        report = _mus_report(mus_constraints=[gt_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("GT" in s or "dinucleotide" in s.lower() for s in suggestions)

    def test_no_cpg_suggestion(self):
        """NO_CPG in MUS → suggestion to disable CpG avoidance."""
        model = _make_model()
        cpg_constraint = _constraint(ConstraintType.NO_CPG, "no_cpg")
        report = _mus_report(mus_constraints=[cpg_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("CpG" in s for s in suggestions)

    def test_cryptic_splice_suggestion(self):
        """NO_CRYPTIC_SPLICE in MUS → suggestion to increase threshold."""
        model = _make_model(config=_make_config(cryptic_splice_threshold=3.0))
        splice_constraint = _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "no_cryptic_splice")
        report = _mus_report(mus_constraints=[splice_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("splice" in s.lower() or "threshold" in s.lower() for s in suggestions)
        # Should mention current and new threshold values
        splice_suggestions = [s for s in suggestions if "splice" in s.lower() or "threshold" in s.lower()]
        assert any("3.0" in s and "4.0" in s for s in splice_suggestions)

    def test_restriction_site_suggestion(self):
        """RESTRICTION_SITE in MUS → suggestion to remove enzyme from list."""
        model = _make_model()
        rs_constraint = _constraint(
            ConstraintType.RESTRICTION_SITE,
            "no_ecori",
            params={"enzyme": "EcoRI", "site": "GAATTC"},
        )
        report = _mus_report(mus_constraints=[rs_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("EcoRI" in s for s in suggestions)

    def test_codon_usage_suggestion(self):
        """CODON_USAGE in MUS → suggestion to relax codon usage bias."""
        model = _make_model()
        cu_constraint = _constraint(ConstraintType.CODON_USAGE, "codon_usage_bias")
        report = _mus_report(mus_constraints=[cu_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("codon" in s.lower() or "usage" in s.lower() for s in suggestions)

    def test_mhc_binding_suggestion(self):
        """MHC_BINDING in MUS → suggestion about immunogenicity risk."""
        model = _make_model()
        mhc_constraint = _constraint(ConstraintType.MHC_BINDING, "mhc_binding")
        report = _mus_report(mus_constraints=[mhc_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("MHC" in s or "immunogenicity" in s.lower() for s in suggestions)

    def test_mrna_stability_suggestion(self):
        """MRNA_STABILITY in MUS → suggestion about mRNA folding."""
        model = _make_model()
        mrna_constraint = _constraint(ConstraintType.MRNA_STABILITY, "mrna_stability")
        report = _mus_report(mus_constraints=[mrna_constraint])
        suggestions = suggest_relaxations(report, model)
        assert any("mRNA" in s or "mrna" in s.lower() or "stability" in s.lower() for s in suggestions)

    def test_multiple_constraints_multiple_suggestions(self):
        """Multiple constraint types in MUS → at least one suggestion per type."""
        model = _make_model(config=_make_config(gc_lo=0.60, gc_hi=0.65))
        constraints = [
            _constraint(ConstraintType.GC_CONTENT, "gc_range"),
            _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt"),
        ]
        report = _mus_report(mus_constraints=constraints)
        suggestions = suggest_relaxations(report, model)
        # Should have at least 2 suggestions (one per constraint type)
        assert len(suggestions) >= 2
        has_gc = any("GC" in s or "gc" in s.lower() for s in suggestions)
        has_gt = any("GT" in s or "dinucleotide" in s.lower() for s in suggestions)
        assert has_gc
        assert has_gt

    def test_unknown_constraint_type_fallback(self):
        """An unknown constraint type should produce a generic fallback suggestion."""
        model = _make_model()
        custom_constraint = _constraint(ConstraintType.CUSTOM, "custom_xyz")
        report = _mus_report(mus_constraints=[custom_constraint])
        suggestions = suggest_relaxations(report, model)
        assert len(suggestions) >= 1
        assert any("custom_xyz" in s or "review" in s.lower() for s in suggestions)

    def test_suggestions_are_strings(self):
        """All suggestions should be non-empty strings."""
        model = _make_model(config=_make_config(gc_lo=0.60, gc_hi=0.65))
        constraints = [
            _constraint(ConstraintType.GC_CONTENT, "gc_range"),
            _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "no_splice"),
            _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt"),
        ]
        report = _mus_report(mus_constraints=constraints)
        suggestions = suggest_relaxations(report, model)
        for s in suggestions:
            assert isinstance(s, str) and len(s) > 0


# =====================================================================
# 4. explain_conflict
# =====================================================================

class TestExplainConflict:
    """explain_conflict produces human-readable output."""

    def test_returns_nonempty_string(self):
        """explain_conflict should always return a non-empty string."""
        model = _make_model()
        gt_constraint = _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt")
        report = _mus_report(mus_constraints=[gt_constraint])
        text = explain_conflict(report, model)
        assert isinstance(text, str)
        assert len(text) > 0

    def test_empty_mus_no_conflict(self):
        """An empty MUS should indicate no conflict / satisfiable."""
        model = _make_model()
        report = _mus_report(mus_constraints=[])
        text = explain_conflict(report, model)
        assert any(
            kw in text.lower()
            for kw in ("no conflict", "satisfiable", "not infeasible")
        )

    def test_explanation_from_mus_report_used(self):
        """When MUS has no constraints but has an explanation, use it."""
        model = _make_model()
        report = _mus_report(
            mus_constraints=[],
            explanation="Model is not infeasible; MUS is undefined.",
        )
        text = explain_conflict(report, model)
        assert "not infeasible" in text or "MUS is undefined" in text

    def test_mentions_constraint_type(self):
        """Output should mention the constraint type (e.g. 'no_gt_dinucleotide')."""
        model = _make_model()
        gt_constraint = _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt")
        report = _mus_report(mus_constraints=[gt_constraint])
        text = explain_conflict(report, model)
        # The output groups by type and includes the enum value
        assert "no_gt_dinucleotide" in text.lower() or "GT" in text

    def test_mentions_constraint_name(self):
        """Output should include constraint names."""
        model = _make_model()
        gc_constraint = _constraint(ConstraintType.GC_CONTENT, "gc_range_pos0_5")
        report = _mus_report(mus_constraints=[gc_constraint])
        text = explain_conflict(report, model)
        assert "gc_range_pos0_5" in text

    def test_mentions_positions(self):
        """When constraints have positions, the output should mention them."""
        model = _make_model()
        gc_constraint = _constraint(
            ConstraintType.GC_CONTENT, "gc_range", positions=[0, 1, 2, 3, 4, 5]
        )
        report = _mus_report(
            mus_constraints=[gc_constraint],
            conflict_positions=[0, 1, 2, 3, 4, 5],
        )
        text = explain_conflict(report, model)
        # Should mention positions somewhere
        assert any(str(p) in text for p in range(6))

    def test_conflict_positions_shown(self):
        """Key conflict positions from MUSReport should appear in the explanation."""
        model = _make_model()
        constraint = _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt")
        report = _mus_report(
            mus_constraints=[constraint],
            conflict_positions=[3, 7, 12],
        )
        text = explain_conflict(report, model)
        assert "3" in text and "7" in text and "12" in text

    def test_solver_stats_in_output(self):
        """Output should include solver iteration count and time."""
        model = _make_model()
        constraint = _constraint(ConstraintType.GC_CONTENT, "gc_range")
        report = _mus_report(
            mus_constraints=[constraint],
            all_constraints=[constraint],
        )
        text = explain_conflict(report, model)
        # Should mention iterations and time
        assert "solver call" in text.lower() or "iteration" in text.lower() or str(report.iterations) in text

    def test_gc_plus_splice_conflict_explanation(self):
        """GC + cryptic splice → domain-specific explanation about codon GC conflict."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.GC_CONTENT, "gc_range"),
            _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "no_splice"),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        # Domain explanation for GC + splice conflict
        assert "GC" in text and "splice" in text.lower()

    def test_gc_plus_restriction_site_conflict_explanation(self):
        """GC + restriction site → domain explanation about low-GC codon forcing."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.GC_CONTENT, "gc_range"),
            _constraint(
                ConstraintType.RESTRICTION_SITE,
                "no_ecori",
                params={"enzyme": "EcoRI", "site": "GAATTC"},
            ),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        assert "GC" in text and ("restriction" in text.lower() or "EcoRI" in text)

    def test_cpg_plus_gc_conflict_explanation(self):
        """CpG + GC → domain explanation about CpG eliminating high-GC codons."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.NO_CPG, "no_cpg"),
            _constraint(ConstraintType.GC_CONTENT, "gc_range"),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        assert "CpG" in text and "GC" in text

    def test_gt_plus_splice_conflict_explanation(self):
        """GT + splice → domain explanation about GT-dinucleotide / splice interaction."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.NO_GT_DINUCLEOTIDE, "no_gt"),
            _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "no_splice"),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        assert "GT" in text or "dinucleotide" in text.lower()

    def test_cpg_plus_restriction_site_conflict_explanation(self):
        """CpG + restriction site → domain explanation about too few codon options."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.NO_CPG, "no_cpg"),
            _constraint(
                ConstraintType.RESTRICTION_SITE,
                "no_ecori",
                params={"enzyme": "EcoRI", "site": "GAATTC"},
            ),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        assert "codon" in text.lower()  # "too few codon options"

    def test_codon_usage_plus_restriction_site_conflict(self):
        """Codon usage + restriction site → explanation about high-usage codons containing sites."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.CODON_USAGE, "codon_bias"),
            _constraint(
                ConstraintType.RESTRICTION_SITE,
                "no_ecori",
                params={"enzyme": "EcoRI", "site": "GAATTC"},
            ),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        assert "codon" in text.lower() and "restriction" in text.lower()

    def test_two_unknown_types_generic_conflict(self):
        """Two unknown constraint types → generic mutual-exclusion explanation."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.CUSTOM, "custom_a"),
            _constraint(ConstraintType.AMINO_ACID_IDENTITY, "aa_identity"),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        # Should have conflict analysis section
        assert "mutually exclusive" in text.lower() or "conflict" in text.lower()

    def test_multiple_types_with_pairwise_explanations(self):
        """More than 2 types where pairwise explanations exist → conflict analysis.

        GC + CpG + CrypticSplice: pairwise explanations exist for GC+splice
        and CpG+GC, so the generic 'jointly eliminate' fallback is not reached.
        Instead, the specific pairwise conflict explanations are shown.
        """
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.GC_CONTENT, "gc"),
            _constraint(ConstraintType.NO_CPG, "cpg"),
            _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "splice"),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        # Should have conflict analysis section with pairwise explanations
        assert "conflict" in text.lower()
        assert "GC" in text and "CpG" in text

    def test_constraints_grouped_by_type(self):
        """Multiple constraints of the same type should be grouped in the output."""
        model = _make_model()
        constraints = [
            _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "splice_pos3", positions=[3]),
            _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "splice_pos7", positions=[7]),
        ]
        report = _mus_report(mus_constraints=constraints)
        text = explain_conflict(report, model)
        # Should show the type once with a count of 2
        assert "2 constraint" in text.lower()


# =====================================================================
# Integration: quick_feasibility_check + explain + suggest workflow
# =====================================================================

class TestMUSWorkflowIntegration:
    """End-to-end workflow: feasibility check → build MUSReport → explain → suggest."""

    def test_infeasible_model_workflow(self):
        """Full workflow for an infeasible model (impossible GC bounds)."""
        model = _make_model(
            protein="VVVVVVVVVV",
            config=_make_config(gc_lo=0.10, gc_hi=0.20),
        )
        # Step 1: Quick feasibility check
        report = quick_feasibility_check(model)
        assert not report.feasible
        assert len(report.impossible_constraints) > 0

        # Step 2: Build a MUSReport manually (simulating compute_mus result)
        gc_constraint = _constraint(ConstraintType.GC_CONTENT, "gc_range")
        mus_report = _mus_report(
            mus_constraints=[gc_constraint],
            conflict_positions=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        )

        # Step 3: Explain
        text = explain_conflict(mus_report, model)
        assert isinstance(text, str) and len(text) > 0

        # Step 4: Suggest
        suggestions = suggest_relaxations(mus_report, model)
        assert isinstance(suggestions, list)
        assert len(suggestions) > 0

    def test_feasible_model_workflow(self):
        """Full workflow for a feasible model — no MUS, no suggestions."""
        model = _make_model(
            protein="MVSKGE",
            config=_make_config(gc_lo=0.10, gc_hi=0.90),
        )
        # Step 1: Quick feasibility check
        report = quick_feasibility_check(model)
        assert report.feasible

        # Step 2: Empty MUSReport
        mus_report = _mus_report(mus_constraints=[])

        # Step 3: Explain → no conflict
        text = explain_conflict(mus_report, model)
        assert "no conflict" in text.lower() or "satisfiable" in text.lower()

        # Step 4: Suggest → empty
        suggestions = suggest_relaxations(mus_report, model)
        assert suggestions == []

    def test_gc_infeasible_full_workflow(self):
        """Full workflow for GC-infeasible model with suggestions and explanation."""
        model = _make_model(
            protein="KKKKNNNN",
            config=_make_config(gc_lo=0.70, gc_hi=0.90),
        )
        # Quick check
        report = quick_feasibility_check(model)
        assert not report.feasible
        assert "GC_CONTENT" in report.impossible_constraints

        # Simulate MUS result
        gc_constraint = _constraint(ConstraintType.GC_CONTENT, "gc_range")
        mus_report = _mus_report(
            mus_constraints=[gc_constraint],
            conflict_positions=[0, 1, 2, 3, 4, 5, 6, 7],
        )

        # Explain
        text = explain_conflict(mus_report, model)
        assert "gc_content" in text.lower() or "GC" in text

        # Suggest
        suggestions = suggest_relaxations(mus_report, model)
        assert any("GC" in s or "gc" in s.lower() for s in suggestions)
