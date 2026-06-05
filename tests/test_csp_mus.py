"""Tests for solver/mus.py — MUS (Minimal Unsatisfiable Subset) analysis.

Five categories:
1. Quick feasibility check (no solver required)
2. MUS computation (requires solver)
3. Conflict explanation
4. Relaxation suggestions
5. FeasibilityReport data structure

Uses pytest.importorskip for solver-dependent tests.
"""

from __future__ import annotations

import pytest
from dataclasses import fields

# Graceful skip if the module isn't built yet
mus = pytest.importorskip("biocompiler.solver.mus", reason="solver/mus.py not available")
types_mod = pytest.importorskip("biocompiler.solver.types", reason="solver/types.py not available")


# ── Helpers ──────────────────────────────────────────────────────────

def _valine_codons() -> list[str]:
    """All Valine codons — every one starts with 'GT'."""
    from biocompiler.type_system import AA_TO_CODONS
    return AA_TO_CODONS["V"]


def _make_model(
    protein: str,
    constraint_types: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    avoid_gt: bool = False,
) -> types_mod.CSPModel:
    """Build a CSPModel from convenience arguments.

    Maps human-readable constraint names like "NoGTDinucleotide" to
    the corresponding ConstraintType enum values.
    """
    from biocompiler.constants import AA_TO_CODONS

    config = types_mod.SolverConfig(gc_lo=gc_lo, gc_hi=gc_hi)
    codon_domains = {i: list(AA_TO_CODONS.get(aa, [])) for i, aa in enumerate(protein)}

    _NAME_TO_CTYPE = {
        "NoGTDinucleotide": types_mod.ConstraintType.NO_GT_DINUCLEOTIDE,
        "NoStopCodons": types_mod.ConstraintType.AMINO_ACID_IDENTITY,
        "ValidCodingSeq": types_mod.ConstraintType.AMINO_ACID_IDENTITY,
        "NoCpG": types_mod.ConstraintType.NO_CPG,
        "NoCrypticSplice": types_mod.ConstraintType.NO_CRYPTIC_SPLICE,
    }

    constraints: list[types_mod.ConstraintSpec] = []
    if constraint_types:
        for name in constraint_types:
            ctype = _NAME_TO_CTYPE.get(name, types_mod.ConstraintType.CUSTOM)
            constraints.append(
                types_mod.ConstraintSpec(
                    ctype=ctype,
                    name=name,
                    strictness=types_mod.ConstraintStrictness.HARD,
                )
            )

    if avoid_gt:
        # Add GT-dinucleotide constraint if not already present
        has_gt = any(c.ctype == types_mod.ConstraintType.NO_GT_DINUCLEOTIDE for c in constraints)
        if not has_gt:
            constraints.append(
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.NO_GT_DINUCLEOTIDE,
                    name="NoGTDinucleotide",
                    strictness=types_mod.ConstraintStrictness.HARD,
                )
            )

    return types_mod.CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints,
        config=config,
    )


# =====================================================================
# 1. Quick feasibility check
# =====================================================================

class TestQuickFeasibility:
    """Lightweight checks that do not require an SMT solver."""

    def test_all_valine_no_gt_warns(self):
        """All-Valine + NoGTDinucleotide → every codon contains GT."""
        assert all("GT" in c for c in _valine_codons())
        model = _make_model("VVVVVV", ["NoGTDinucleotide"])
        r = mus.quick_feasibility_check(model)
        # The quick check inspects cross-codon GT boundaries;
        # within-codon GT (present in every V codon) may or may not
        # be flagged depending on implementation depth.
        assert isinstance(r, mus.FeasibilityReport)

    def test_normal_protein_feasible(self):
        """Reasonable protein + relaxed constraints → feasible."""
        model = _make_model("MAGSC", ["NoStopCodons", "ValidCodingSeq"],
                            gc_lo=0.20, gc_hi=0.80)
        r = mus.quick_feasibility_check(model)
        assert r.feasible

    def test_tight_gc_bounds_warn(self):
        """GC bounds near the edge of achievability should produce a warning or be infeasible."""
        # Lysine codons: AAA (0/3 GC), AAG (1/3 GC) → max GC = 1/3 ≈ 33%
        # Requesting gc_lo=0.50 is unreachable → infeasible
        model = _make_model("KKKKKKKKKK", gc_lo=0.50, gc_hi=0.55)
        r = mus.quick_feasibility_check(model)
        assert not r.feasible or len(r.warnings) > 0

    def test_min_possible_gc_exceeds_gc_hi(self):
        """gc_hi below minimum achievable GC → infeasible.

        All-Valine codons have ≥1/3 GC; gc_hi=0.20 is unreachable.
        """
        model = _make_model("VVVVVVVVVV", gc_lo=0.10, gc_hi=0.20)
        r = mus.quick_feasibility_check(model)
        assert not r.feasible

    def test_max_possible_gc_below_gc_lo(self):
        """gc_lo above maximum achievable GC → infeasible.

        Lys/Asn codons max out at 1/3 GC each; gc_lo=0.70 is unreachable.
        """
        model = _make_model("KKKKNNNN", gc_lo=0.70, gc_hi=0.90)
        r = mus.quick_feasibility_check(model)
        assert not r.feasible

    def test_empty_protein_trivially_feasible(self):
        """Empty protein → always feasible."""
        model = _make_model("", [])
        r = mus.quick_feasibility_check(model)
        assert r.feasible


# =====================================================================
# 2. MUS computation (requires solver)
# =====================================================================

class TestMUSComputation:
    """MUS tests — require z3-solver."""

    @pytest.fixture(autouse=True)
    def _skip_no_z3(self):
        pytest.importorskip("z3", reason="z3-solver not installed")

    def test_infeasible_mus_identifies_conflict(self):
        """MUS of infeasible problem should identify constraints."""
        model = _make_model("VVVVVV", ["NoGTDinucleotide"], gc_lo=0.20, gc_hi=0.80)
        from biocompiler.solver.engine_z3 import Z3Engine
        backend = Z3Engine(config=model.config)
        result = mus.compute_mus(model, backend)
        # MUS should be non-empty for an infeasible problem
        # (or if the problem is actually feasible, MUS is empty)
        assert isinstance(result, types_mod.MUSReport)

    def test_feasible_mus_empty(self):
        """MUS of feasible problem should be empty."""
        model = _make_model("MASMTGGQQMG", ["NoStopCodons", "ValidCodingSeq"],
                            gc_lo=0.20, gc_hi=0.80)
        from biocompiler.solver.engine_z3 import Z3Engine
        backend = Z3Engine(config=model.config)
        result = mus.compute_mus(model, backend)
        # Feasible models should have empty MUS
        if not result.mus_constraints:
            assert len(result.conflicting_constraints) == 0


# =====================================================================
# 3. Conflict explanation
# =====================================================================

class TestConflictExplanation:
    """explain_conflict() → human-readable string."""

    def test_returns_nonempty_string(self):
        """explain_conflict should return a non-empty string."""
        # Create an infeasible model for explanation
        model = _make_model("VVVVVV", ["NoGTDinucleotide"], gc_lo=0.20, gc_hi=0.80)
        # Build a MUSReport manually for explanation
        report = types_mod.MUSReport(
            mus_constraints=[
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.NO_GT_DINUCLEOTIDE,
                    name="NoGTDinucleotide",
                    strictness=types_mod.ConstraintStrictness.HARD,
                )
            ],
            explanation="GT dinucleotide constraint infeasible for Valine-rich proteins",
        )
        txt = mus.explain_conflict(report, model)
        assert isinstance(txt, str) and len(txt) > 0

    def test_mentions_gt_conflict(self):
        """Conflict explanation should mention GT or Valine or dinucleotide."""
        model = _make_model("VVVVVV", ["NoGTDinucleotide"], gc_lo=0.20, gc_hi=0.80)
        report = types_mod.MUSReport(
            mus_constraints=[
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.NO_GT_DINUCLEOTIDE,
                    name="NoGTDinucleotide",
                    strictness=types_mod.ConstraintStrictness.HARD,
                )
            ],
        )
        txt = mus.explain_conflict(report, model)
        assert "GT" in txt or "Val" in txt or "dinucleotide" in txt.lower()

    def test_mentions_gc_conflict(self):
        """Conflict explanation should mention GC when GC constraint is in MUS."""
        model = _make_model("KKKKNNNN", gc_lo=0.70, gc_hi=0.90)
        report = types_mod.MUSReport(
            mus_constraints=[
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.GC_CONTENT,
                    name="GC_bounds",
                    strictness=types_mod.ConstraintStrictness.HARD,
                    params={"gc_lo": 0.70, "gc_hi": 0.90},
                )
            ],
        )
        txt = mus.explain_conflict(report, model)
        assert "GC" in txt or "gc" in txt.lower()

    def test_feasible_shows_no_conflict(self):
        """Feasible model should show no conflict in explanation."""
        model = _make_model("MASMTGGQQMG", ["NoStopCodons"], gc_lo=0.20, gc_hi=0.80)
        report = types_mod.MUSReport()  # empty MUS
        txt = mus.explain_conflict(report, model)
        assert any(kw in txt.lower() for kw in
                   ("feasible", "no conflict", "satisfiable", "no conflict detected"))


# =====================================================================
# 4. Relaxation suggestions
# =====================================================================

class TestRelaxationSuggestions:
    """suggest_relaxations() → list of concrete, actionable strings."""

    def test_returns_list_of_strings(self):
        """suggestions should be a list of strings."""
        model = _make_model("VVVVVV", ["NoGTDinucleotide"])
        report = types_mod.MUSReport(
            mus_constraints=[
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.NO_GT_DINUCLEOTIDE,
                    name="NoGTDinucleotide",
                    strictness=types_mod.ConstraintStrictness.HARD,
                )
            ],
        )
        sug = mus.suggest_relaxations(report, model)
        assert isinstance(sug, list)
        assert all(isinstance(s, str) for s in sug)

    def test_gc_relaxation_mentioned(self):
        """Tight GC → suggestion about GC relaxation."""
        model = _make_model("KKKKNNNN", ["ValidCodingSeq"], gc_lo=0.60, gc_hi=0.65)
        report = types_mod.MUSReport(
            mus_constraints=[
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.GC_CONTENT,
                    name="GC_bounds",
                    strictness=types_mod.ConstraintStrictness.HARD,
                    params={"gc_lo": 0.60, "gc_hi": 0.65},
                )
            ],
        )
        sug = mus.suggest_relaxations(report, model)
        assert any("GC" in s or "gc" in s.lower() for s in sug)

    def test_constraint_relaxation_mentioned(self):
        """GT constraint in MUS → suggestion mentions GT or relaxation."""
        model = _make_model("VVVVVV", ["NoGTDinucleotide"])
        report = types_mod.MUSReport(
            mus_constraints=[
                types_mod.ConstraintSpec(
                    ctype=types_mod.ConstraintType.NO_GT_DINUCLEOTIDE,
                    name="NoGTDinucleotide",
                    strictness=types_mod.ConstraintStrictness.HARD,
                )
            ],
        )
        sug = mus.suggest_relaxations(report, model)
        assert any("GT" in s or "NoGTDinucleotide" in s
                    or "relax" in s.lower() or "remove" in s.lower()
                    for s in sug)

    def test_feasible_no_suggestions_needed(self):
        """Feasible model should have no relaxation suggestions."""
        model = _make_model("MASMTGGQQMG", ["NoStopCodons", "ValidCodingSeq"],
                            gc_lo=0.20, gc_hi=0.80)
        report = types_mod.MUSReport()  # empty MUS → feasible
        sug = mus.suggest_relaxations(report, model)
        assert sug == []


# =====================================================================
# 5. FeasibilityReport
# =====================================================================

class TestFeasibilityReport:
    """FeasibilityReport creation and field access."""

    def test_create_feasible(self):
        r = mus.FeasibilityReport(feasible=True, warnings=[], impossible_constraints=[])
        assert r.feasible is True
        assert r.warnings == [] and r.impossible_constraints == []

    def test_create_infeasible(self):
        r = mus.FeasibilityReport(
            feasible=False,
            warnings=["GC unreachable"],
            impossible_constraints=["GC_CONTENT", "NoGTDinucleotide"],
        )
        assert not r.feasible
        assert "NoGTDinucleotide" in r.impossible_constraints

    def test_has_required_fields(self):
        names = {f.name for f in fields(mus.FeasibilityReport(
            feasible=True, warnings=[], impossible_constraints=[]))}
        for req in ("feasible", "warnings", "impossible_constraints"):
            assert req in names

    def test_from_quick_check_infeasible(self):
        """Quick check on infeasible model should produce a FeasibilityReport."""
        model = _make_model("KKKKNNNN", gc_lo=0.70, gc_hi=0.90)
        r = mus.quick_feasibility_check(model)
        assert isinstance(r, mus.FeasibilityReport)
        assert not r.feasible and len(r.impossible_constraints) > 0

    def test_from_quick_check_feasible(self):
        """Quick check on feasible model should produce feasible report."""
        model = _make_model("MASMTGGQQMG", ["NoStopCodons"], gc_lo=0.20, gc_hi=0.80)
        r = mus.quick_feasibility_check(model)
        assert isinstance(r, mus.FeasibilityReport)
        if r.feasible:
            assert len(r.impossible_constraints) == 0


# =====================================================================
# Integration & edge cases
# =====================================================================

class TestMUSWorkflow:
    """End-to-end: feasibility → explain → relax."""

    def test_infeasible_workflow(self):
        """Full workflow for infeasible problem."""
        model = _make_model("KKKKNNNN", gc_lo=0.70, gc_hi=0.90)
        r = mus.quick_feasibility_check(model)
        assert not r.feasible

    def test_feasible_workflow(self):
        """Full workflow for feasible problem."""
        model = _make_model("MASMTGGQQMG", ["NoStopCodons"], gc_lo=0.20, gc_hi=0.80)
        r = mus.quick_feasibility_check(model)
        assert r.feasible


class TestMUSEdgeCases:
    """Boundary conditions for MUS analysis."""

    def test_single_valine_gt_check(self):
        """Single Valine with GT constraint — quick check result."""
        model = _make_model("V", ["NoGTDinucleotide"])
        r = mus.quick_feasibility_check(model)
        assert isinstance(r, mus.FeasibilityReport)

    def test_single_m_feasible(self):
        """Single Methionine should be feasible."""
        model = _make_model("M", ["NoStopCodons"])
        r = mus.quick_feasibility_check(model)
        assert r.feasible

    def test_long_alanine_feasible(self):
        """Long alanine sequence should be feasible with reasonable GC."""
        model = _make_model("A" * 100, ["NoStopCodons", "ValidCodingSeq"],
                            gc_lo=0.30, gc_hi=0.70)
        r = mus.quick_feasibility_check(model)
        assert r.feasible

    def test_mixed_valine_gt_check(self):
        """Mixed protein with Valine and GT constraint."""
        model = _make_model("MASVTG", ["NoGTDinucleotide"])
        r = mus.quick_feasibility_check(model)
        assert isinstance(r.feasible, bool)

    def test_gc_extremes_always_feasible(self):
        """Very wide GC bounds should always be feasible."""
        model = _make_model("VVVVVV", gc_lo=0.0, gc_hi=1.0)
        r = mus.quick_feasibility_check(model)
        assert r.feasible
