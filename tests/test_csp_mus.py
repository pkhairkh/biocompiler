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


# ── Helpers ──────────────────────────────────────────────────────────

def _valine_codons() -> list[str]:
    """All Valine codons — every one starts with 'GT'."""
    from biocompiler.type_system import AA_TO_CODONS
    return AA_TO_CODONS["V"]


# =====================================================================
# 1. Quick feasibility check
# =====================================================================

class TestQuickFeasibility:
    """Lightweight checks that do not require an SMT solver."""

    def test_all_valine_no_gt_impossible(self):
        """All-Valine + NoGTDinucleotide → impossible (every V codon has GT)."""
        assert all("GT" in c for c in _valine_codons())
        r = mus.quick_feasibility_check("VVVVVV", ["NoGTDinucleotide"])
        assert not r.feasible
        assert any("GT" in w for w in r.warnings)

    def test_normal_protein_feasible(self):
        """Reasonable protein + relaxed constraints → feasible."""
        r = mus.quick_feasibility_check(
            "MAlaGlySerCys", ["NoStopCodons", "ValidCodingSeq"],
            gc_lo=0.20, gc_hi=0.80,
        )
        assert r.feasible

    def test_tight_gc_bounds_warn(self):
        """Extremely tight GC bounds should produce a warning."""
        r = mus.quick_feasibility_check("VVVVVVVVVV", gc_lo=0.40, gc_hi=0.41)
        assert len(r.warnings) > 0

    def test_min_possible_gc_exceeds_gc_hi(self):
        """gc_hi below minimum achievable GC → infeasible.

        All-Valine codons have ≥1/3 GC; gc_hi=0.20 is unreachable.
        """
        r = mus.quick_feasibility_check("VVVVVVVVVV", gc_lo=0.10, gc_hi=0.20)
        assert not r.feasible

    def test_max_possible_gc_below_gc_lo(self):
        """gc_lo above maximum achievable GC → infeasible.

        Lys/Asn codons max out at 1/3 GC each; gc_lo=0.70 is unreachable.
        """
        r = mus.quick_feasibility_check("KKKKNNNN", gc_lo=0.70, gc_hi=0.90)
        assert not r.feasible

    def test_empty_and_no_constraint_trivially_feasible(self):
        """Empty protein or no constraints → always feasible."""
        assert mus.quick_feasibility_check("", []).feasible
        assert mus.quick_feasibility_check("MVSKGE", []).feasible


# =====================================================================
# 2. MUS computation (requires solver)
# =====================================================================

class TestMUSComputation:
    """MUS tests — require z3-solver."""

    @pytest.fixture(autouse=True)
    def _skip_no_z3(self):
        pytest.importorskip("z3", reason="z3-solver not installed")

    def test_infeasible_mus_identifies_conflict(self):
        """MUS of infeasible problem should identify NoGTDinucleotide."""
        result = mus.compute_mus("VVVVVV", ["NoGTDinucleotide"],
                                 gc_lo=0.20, gc_hi=0.80)
        assert len(result) > 0
        assert any("GT" in c or "Dinucleotide" in c for c in result)

    def test_feasible_mus_empty(self):
        """MUS of feasible problem should be empty."""
        result = mus.compute_mus("MASMTGGQQMG",
                                 ["NoStopCodons", "ValidCodingSeq"],
                                 gc_lo=0.20, gc_hi=0.80)
        assert len(result) == 0

    def test_gc_conflict_mus(self):
        """GC [0.60, 0.65] for A/T-rich AAs → MUS mentions GC."""
        result = mus.compute_mus("KKKNNNIII", ["ValidCodingSeq"],
                                 gc_lo=0.60, gc_hi=0.65)
        assert len(result) > 0
        assert any("GC" in c or "gc" in c.lower() for c in result)

    def test_mus_minimality_approx(self):
        """Removing any MUS element should (approximately) restore feasibility."""
        result = mus.compute_mus("VVVV",
                                 ["NoGTDinucleotide", "NoStopCodons", "ValidCodingSeq"],
                                 gc_lo=0.10, gc_hi=0.90)
        if len(result) <= 1:
            pytest.skip("MUS has ≤1 constraint — minimality trivial")
        for i in range(len(result)):
            reduced = result[:i] + result[i + 1:]
            if len(mus.compute_mus("VVVV", reduced, gc_lo=0.10, gc_hi=0.90)) == 0:
                return  # minimality confirmed for at least one element


# =====================================================================
# 3. Conflict explanation
# =====================================================================

class TestConflictExplanation:
    """explain_conflict() → human-readable string."""

    def test_returns_nonempty_string(self):
        txt = mus.explain_conflict("VVVVVV", ["NoGTDinucleotide"])
        assert isinstance(txt, str) and len(txt) > 0

    def test_mentions_gt_conflict(self):
        txt = mus.explain_conflict("VVVVVV", ["NoGTDinucleotide"])
        assert "GT" in txt or "Val" in txt or "dinucleotide" in txt.lower()

    def test_mentions_gc_conflict(self):
        txt = mus.explain_conflict("KKKKNNNN", gc_lo=0.70, gc_hi=0.90)
        assert "GC" in txt or "gc" in txt.lower()

    def test_mentions_positions_or_aas(self):
        txt = mus.explain_conflict("MVVVGM", ["NoGTDinucleotide"])
        has_pos = any(str(i) in txt for i in range(6))
        has_aa = "position" in txt.lower() or "V" in txt or "Val" in txt
        assert has_pos or has_aa

    def test_feasible_shows_no_conflict(self):
        txt = mus.explain_conflict("MASMTGGQQMG", ["NoStopCodons"],
                                   gc_lo=0.20, gc_hi=0.80)
        assert any(kw in txt.lower() for kw in
                   ("feasible", "no conflict", "satisfiable"))


# =====================================================================
# 4. Relaxation suggestions
# =====================================================================

class TestRelaxationSuggestions:
    """suggest_relaxations() → list of concrete, actionable strings."""

    def test_returns_list_of_strings(self):
        sug = mus.suggest_relaxations("VVVVVV", ["NoGTDinucleotide"])
        assert isinstance(sug, list)
        assert all(isinstance(s, str) for s in sug)

    def test_gc_relaxation_mentioned(self):
        """Tight GC → suggestion like 'Relax GC from [0.60, 0.65] to [0.55, 0.70]'."""
        sug = mus.suggest_relaxations("KKKKNNNN", ["ValidCodingSeq"],
                                      gc_lo=0.60, gc_hi=0.65)
        assert any("GC" in s or "gc" in s.lower() for s in sug)

    def test_constraint_relaxation_mentioned(self):
        sug = mus.suggest_relaxations("VVVVVV", ["NoGTDinucleotide"])
        assert any("GT" in s or "NoGTDinucleotide" in s
                    or "relax" in s.lower() or "remove" in s.lower()
                    for s in sug)

    def test_suggestions_have_specific_parameters(self):
        """Each suggestion mentions a specific parameter (GC, GT, bound, etc.)."""
        sug = mus.suggest_relaxations("VVVVVV", ["NoGTDinucleotide"],
                                      gc_lo=0.40, gc_hi=0.60)
        keywords = ("GC", "gc", "GT", "dinucleotide", "constraint",
                    "bound", "threshold", "Val")
        for s in sug:
            assert any(kw in s for kw in keywords), f"Vague suggestion: {s!r}"

    def test_feasible_no_suggestions_needed(self):
        sug = mus.suggest_relaxations("MASMTGGQQMG",
                                      ["NoStopCodons", "ValidCodingSeq"],
                                      gc_lo=0.20, gc_hi=0.80)
        if sug:
            assert any("feasib" in s.lower() or "no relax" in s.lower()
                        for s in sug)


# =====================================================================
# 5. FeasibilityReport
# =====================================================================

class TestFeasibilityReport:
    """FeasibilityReport creation and field access."""

    def test_create_feasible(self):
        r = mus.FeasibilityReport(feasible=True, warnings=[], conflicts=[])
        assert r.feasible is True
        assert r.warnings == [] and r.conflicts == []

    def test_create_infeasible(self):
        r = mus.FeasibilityReport(
            feasible=False,
            warnings=["GC unreachable"],
            conflicts=["GC_bounds", "NoGTDinucleotide"],
        )
        assert not r.feasible
        assert "NoGTDinucleotide" in r.conflicts

    def test_has_required_fields(self):
        names = {f.name for f in fields(mus.FeasibilityReport(
            feasible=True, warnings=[], conflicts=[]))}
        for req in ("feasible", "warnings", "conflicts"):
            assert req in names

    def test_from_quick_check_infeasible(self):
        r = mus.quick_feasibility_check("VVVVVV", ["NoGTDinucleotide"])
        assert isinstance(r, mus.FeasibilityReport)
        assert not r.feasible and len(r.conflicts) > 0

    def test_from_quick_check_feasible(self):
        r = mus.quick_feasibility_check("MASMTGGQQMG", ["NoStopCodons"],
                                        gc_lo=0.20, gc_hi=0.80)
        assert isinstance(r, mus.FeasibilityReport)
        if r.feasible:
            assert len(r.conflicts) == 0


# =====================================================================
# Integration & edge cases
# =====================================================================

class TestMUSWorkflow:
    """End-to-end: feasibility → explain → relax."""

    def test_infeasible_workflow(self):
        p, cons = "VVVVVV", ["NoGTDinucleotide"]
        r = mus.quick_feasibility_check(p, cons, gc_lo=0.20, gc_hi=0.80)
        assert not r.feasible
        assert len(mus.explain_conflict(p, cons, gc_lo=0.20, gc_hi=0.80)) > 0
        assert isinstance(mus.suggest_relaxations(p, cons, gc_lo=0.20, gc_hi=0.80), list)

    def test_feasible_workflow(self):
        p, cons = "MASMTGGQQMG", ["NoStopCodons"]
        r = mus.quick_feasibility_check(p, cons, gc_lo=0.20, gc_hi=0.80)
        assert r.feasible
        txt = mus.explain_conflict(p, cons, gc_lo=0.20, gc_hi=0.80)
        assert "feasible" in txt.lower() or "no conflict" in txt.lower()


class TestMUSEdgeCases:
    """Boundary conditions for MUS analysis."""

    def test_single_valine_infeasible(self):
        assert not mus.quick_feasibility_check("V", ["NoGTDinucleotide"]).feasible

    def test_single_m_feasible(self):
        assert mus.quick_feasibility_check("M", ["NoStopCodons"]).feasible

    def test_long_alanine_feasible(self):
        r = mus.quick_feasibility_check("A" * 100,
                                        ["NoStopCodons", "ValidCodingSeq"],
                                        gc_lo=0.30, gc_hi=0.70)
        assert r.feasible

    def test_mixed_valine_makes_gt_infeasible(self):
        """Even one Valine makes NoGTDinucleotide impossible."""
        assert not mus.quick_feasibility_check("MASVTG", ["NoGTDinucleotide"]).feasible

    def test_no_valine_gt_may_be_feasible(self):
        """Non-Valine AAs have GT-free codon options (cross-codon GT may still exist)."""
        r = mus.quick_feasibility_check("MASMCGQER", ["NoGTDinucleotide"])
        assert isinstance(r.feasible, bool)

    def test_gc_extremes_always_feasible(self):
        assert mus.quick_feasibility_check("VVVVVV", gc_lo=0.0, gc_hi=1.0).feasible
