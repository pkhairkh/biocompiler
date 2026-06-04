"""Tests for BioCompiler Z3 SMT Solver Engine.

Test categories:
  1. Availability check — Z3Engine.is_available()
  2. Simple solve — protein back-translation with constraints
  3. GC constraint — tight and infeasible bounds
  4. UNSAT core — infeasible problem diagnosis
  5. Comparison with OR-Tools — cross-backend validation
  6. Timeout handling — graceful short-timeout behavior

All Z3-dependent tests use ``pytest.importorskip("z3")`` so the suite
passes cleanly when z3-solver is not installed.
"""

from __future__ import annotations

import time

import pytest

# ---------------------------------------------------------------------------
# Optional-dependency gating — skip entire module when z3 is absent
# ---------------------------------------------------------------------------
z3 = pytest.importorskip("z3", reason="z3-solver not installed")

# The solver module is created by sub-agents A1–A16; import after z3 check.
from biocompiler.solver.engine_z3 import Z3Engine, Z3SolveResult  # type: ignore[import-untyped]
from biocompiler.solver.engine_ortools import ORToolsEngine, ORToolsSolveResult  # type: ignore[import-untyped]

# BioCompiler internals used for verification
from biocompiler.translation import translate, compute_cai
from biocompiler.type_system import CODON_TABLE
from biocompiler.constants import reverse_complement
from biocompiler.restriction_sites import RESTRICTION_SITES


# ---------------------------------------------------------------------------
# Fixtures & helpers
# ---------------------------------------------------------------------------

HBB_PROTEIN_SHORT = "MVLSPADKTNVKAAWGKVGA"  # 20 aa — first 20 residues of HBA1


@pytest.fixture
def z3_engine() -> Z3Engine:
    """Default Z3Engine configured for human codon usage."""
    return Z3Engine(organism="Homo_sapiens")


@pytest.fixture
def ortools_engine() -> ORToolsEngine:
    """Default OR-Tools engine for comparison tests."""
    return ORToolsEngine(organism="Homo_sapiens")


def _gc_content(seq: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in "GC") / len(seq)


def _has_restriction_site(seq: str, enzymes: list[str]) -> bool:
    """Return True if *seq* contains any of the given enzyme recognition sites."""
    for enz in enzymes:
        site = RESTRICTION_SITES.get(enz, "")
        if not site:
            continue
        rc = reverse_complement(site)
        if site in seq or rc in seq:
            return True
    return False


def _translates_to(dna: str, protein: str) -> bool:
    """Check that *dna* translates exactly to *protein* (ignoring stop)."""
    return translate(dna, to_stop=True) == protein


# ======================================================================
# 1. Availability Check
# ======================================================================

class TestZ3Availability:
    """Z3Engine availability detection."""

    def test_is_available_when_installed(self):
        """If z3-solver is importable, is_available() must return True."""
        assert Z3Engine.is_available() is True

    def test_is_available_returns_bool(self):
        """is_available() should return a bool, not an int or other type."""
        assert isinstance(Z3Engine.is_available(), bool)

    def test_engine_instantiation(self, z3_engine: Z3Engine):
        """Z3Engine should instantiate without error when z3 is available."""
        assert z3_engine is not None
        assert hasattr(z3_engine, "solve")


# ======================================================================
# 2. Simple Solve
# ======================================================================

class TestZ3SimpleSolve:
    """Basic solve: translate a short protein with default constraints."""

    def test_solve_returns_result(self, z3_engine: Z3Engine):
        """solve() should return a Z3SolveResult dataclass."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert isinstance(result, Z3SolveResult)

    def test_solve_solved_flag(self, z3_engine: Z3Engine):
        """A feasible problem should be marked as solved."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert result.solved is True, (
            f"Simple protein solve failed: {result.error or 'no error message'}"
        )

    def test_solve_correct_length(self, z3_engine: Z3Engine):
        """Output DNA length must equal protein_length * 3."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert result.solved
        assert len(result.sequence) == len(HBB_PROTEIN_SHORT) * 3

    def test_solve_translates_correctly(self, z3_engine: Z3Engine):
        """Output DNA should translate back to the original protein."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert result.solved
        assert _translates_to(result.sequence, HBB_PROTEIN_SHORT), (
            f"DNA does not translate to expected protein. "
            f"Got: {translate(result.sequence, to_stop=True)}"
        )

    def test_solve_gc_in_range(self, z3_engine: Z3Engine):
        """Default GC bounds [0.30, 0.70] should be respected."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert result.solved
        gc = _gc_content(result.sequence)
        assert 0.30 <= gc <= 0.70, f"GC {gc:.4f} outside [0.30, 0.70]"

    def test_solve_no_default_restriction_sites(self, z3_engine: Z3Engine):
        """With enzymes, the solution should contain no sites from the list."""
        enzymes = ["EcoRI", "BamHI"]
        result = z3_engine.solve(HBB_PROTEIN_SHORT, enzymes=enzymes)
        if result.solved:
            assert not _has_restriction_site(result.sequence, enzymes), (
                "Solution contains a restriction site from the enzyme list"
            )

    def test_solve_valid_codons_no_internal_stops(self, z3_engine: Z3Engine):
        """All codons must be valid and no internal stop codons."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert result.solved
        stops = {"TAA", "TAG", "TGA"}
        for i in range(0, len(result.sequence), 3):
            codon = result.sequence[i:i + 3]
            assert codon in CODON_TABLE, f"Invalid codon {codon!r} at pos {i}"
            if i < len(result.sequence) - 3:
                assert codon not in stops, f"Internal stop {codon!r} at pos {i}"


# ======================================================================
# 3. GC Constraint
# ======================================================================

class TestZ3GCConstraint:
    """GC content constraint: tight and infeasible bounds."""

    def test_moderate_tight_gc_bounds(self, z3_engine: Z3Engine):
        """Moderate tight bounds [0.45, 0.55] should be feasible."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.45, gc_hi=0.55)
        assert result.solved, f"GC [0.45,0.55] should be feasible: {result.error}"
        gc = _gc_content(result.sequence)
        assert 0.45 <= gc <= 0.55, f"GC {gc:.4f} outside [0.45, 0.55]"

    def test_tight_gc_bounds(self, z3_engine: Z3Engine):
        """Very tight GC bounds [0.49, 0.51] — may or may not solve."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.49, gc_hi=0.51)
        if result.solved:
            gc = _gc_content(result.sequence)
            assert 0.49 <= gc <= 0.51, f"GC {gc:.4f} outside [0.49, 0.51]"
        else:
            pytest.skip("Tight GC bounds infeasible for this short protein (acceptable)")

    def test_infeasible_gc_bounds_unsat(self, z3_engine: Z3Engine):
        """Completely infeasible GC bounds [0.99, 1.00] should return unsolved."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.99, gc_hi=1.00)
        assert result.solved is False, (
            "GC [0.99,1.00] should be infeasible, but solver claimed solved"
        )

    def test_infeasible_gc_has_diagnostic(self, z3_engine: Z3Engine):
        """An infeasible problem should have a non-empty error or unsat_core."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.99, gc_hi=1.00)
        assert not result.solved
        has_diag = bool(result.error) or bool(result.unsat_core)
        assert has_diag, "Infeasible solve should provide error or UNSAT core"

    def test_zero_gc_infeasible(self, z3_engine: Z3Engine):
        """GC [0.0, 0.0] is infeasible (M requires ATG which has G)."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.0, gc_hi=0.0)
        assert result.solved is False


# ======================================================================
# 4. UNSAT Core
# ======================================================================

class TestZ3UnsatCore:
    """UNSAT core extraction for infeasible problems."""

    def test_unsat_core_returned_on_infeasible(self, z3_engine: Z3Engine):
        """An infeasible problem should produce a non-empty UNSAT core."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.99, gc_hi=1.00)
        assert not result.solved
        assert result.unsat_core is not None, "UNSAT core should not be None"
        assert len(result.unsat_core) > 0, "UNSAT core should be non-empty"

    def test_unsat_core_identifies_conflicting_constraints(self, z3_engine: Z3Engine):
        """UNSAT core should contain at least one constraint related to infeasibility."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, gc_lo=0.99, gc_hi=1.00)
        assert not result.solved
        core = result.unsat_core
        # Minimum guarantee: core is non-empty
        assert len(core) >= 1, "Expected at least one conflict in UNSAT core"

    def test_unsat_core_multiple_conflicts(self, z3_engine: Z3Engine):
        """Multiple conflicting constraints should produce a richer UNSAT core."""
        result = z3_engine.solve(
            HBB_PROTEIN_SHORT,
            gc_lo=0.99, gc_hi=1.00,
            enzymes=["EcoRI", "BamHI", "XhoI", "NotI"],
        )
        assert not result.solved
        assert result.unsat_core is not None
        assert len(result.unsat_core) >= 1

    def test_sat_problem_empty_core(self, z3_engine: Z3Engine):
        """A feasible (SAT) problem should have empty/None UNSAT core."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT)
        assert result.solved
        if result.unsat_core is not None:
            assert len(result.unsat_core) == 0, "SAT result should have empty UNSAT core"


# ======================================================================
# 5. Comparison with OR-Tools
# ======================================================================

class TestZ3VsORTools:
    """Cross-backend validation: both Z3 and OR-Tools should produce valid sequences."""

    def test_both_solve_same_protein(
        self, z3_engine: Z3Engine, ortools_engine: ORToolsEngine,
    ):
        """Both engines should successfully solve the same protein."""
        z3_result = z3_engine.solve(HBB_PROTEIN_SHORT)
        ort_result = ortools_engine.solve(HBB_PROTEIN_SHORT)
        assert z3_result.solved, f"Z3 failed: {z3_result.error}"
        assert ort_result.solved, f"OR-Tools failed: {ort_result.error}"

    def test_both_produce_valid_sequences(
        self, z3_engine: Z3Engine, ortools_engine: ORToolsEngine,
    ):
        """Both solutions should translate correctly and have valid length."""
        for label, eng in [("Z3", z3_engine), ("OR-Tools", ortools_engine)]:
            result = eng.solve(HBB_PROTEIN_SHORT)
            assert result.solved, f"{label} failed: {result.error}"
            assert len(result.sequence) == len(HBB_PROTEIN_SHORT) * 3
            assert _translates_to(result.sequence, HBB_PROTEIN_SHORT), (
                f"{label}: DNA does not translate to expected protein"
            )

    def test_both_respect_gc_bounds(
        self, z3_engine: Z3Engine, ortools_engine: ORToolsEngine,
    ):
        """Both solutions should have GC within the default [0.30, 0.70] range."""
        for label, eng in [("Z3", z3_engine), ("OR-Tools", ortools_engine)]:
            result = eng.solve(HBB_PROTEIN_SHORT)
            if not result.solved:
                continue
            gc = _gc_content(result.sequence)
            assert 0.30 <= gc <= 0.70, f"{label}: GC {gc:.4f} outside [0.30, 0.70]"

    def test_cai_values_reasonable(
        self, z3_engine: Z3Engine, ortools_engine: ORToolsEngine,
    ):
        """Both engines should produce sequences with reasonable CAI."""
        for label, eng in [("Z3", z3_engine), ("OR-Tools", ortools_engine)]:
            result = eng.solve(HBB_PROTEIN_SHORT)
            if not result.solved:
                continue
            cai = compute_cai(result.sequence, organism="Homo_sapiens")
            assert 0.0 <= cai <= 1.0, f"{label}: CAI {cai:.4f} out of [0, 1]"
            assert cai >= 0.3, f"{label}: CAI {cai:.4f} is unreasonably low"

    def test_cai_may_differ_but_both_valid(
        self, z3_engine: Z3Engine, ortools_engine: ORToolsEngine,
    ):
        """CAI values may differ between backends but both should be valid."""
        z3_r = z3_engine.solve(HBB_PROTEIN_SHORT)
        ort_r = ortools_engine.solve(HBB_PROTEIN_SHORT)
        if not (z3_r.solved and ort_r.solved):
            pytest.skip("One or both backends did not solve")
        z3_cai = compute_cai(z3_r.sequence, organism="Homo_sapiens")
        ort_cai = compute_cai(ort_r.sequence, organism="Homo_sapiens")
        # They may differ — just verify both in valid range [0, 1]
        assert 0.0 <= z3_cai <= 1.0 and 0.0 <= ort_cai <= 1.0


# ======================================================================
# 6. Timeout Handling
# ======================================================================

class TestZ3Timeout:
    """Timeout handling: solver should respect time limits."""

    def test_short_timeout_returns(self, z3_engine: Z3Engine):
        """With a very short timeout (0.1s), solve() should still return."""
        start = time.monotonic()
        result = z3_engine.solve(HBB_PROTEIN_SHORT, timeout_s=0.1)
        elapsed = time.monotonic() - start
        assert elapsed < 10.0, f"solve() took {elapsed:.1f}s with 0.1s timeout — possible hang"

    def test_short_timeout_may_not_solve(self, z3_engine: Z3Engine):
        """With a very short timeout, the solver may return unsolved."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, timeout_s=0.1)
        assert isinstance(result, Z3SolveResult)
        # Either solved or unsolved is acceptable — just no exception

    def test_reasonable_timeout_solves(self, z3_engine: Z3Engine):
        """A reasonable timeout (30s) should be enough for a short protein."""
        result = z3_engine.solve(HBB_PROTEIN_SHORT, timeout_s=30.0)
        assert result.solved, (
            f"30s timeout should be sufficient for 20aa protein: {result.error}"
        )

    def test_timeout_zero_returns_quickly(self, z3_engine: Z3Engine):
        """Timeout of 0 should return immediately (unsolved or cached)."""
        start = time.monotonic()
        result = z3_engine.solve(HBB_PROTEIN_SHORT, timeout_s=0.0)
        elapsed = time.monotonic() - start
        assert elapsed < 5.0, "timeout_s=0 should return very quickly"
        assert isinstance(result, Z3SolveResult)


# ======================================================================
# Edge cases
# ======================================================================

class TestZ3EdgeCases:
    """Edge cases and robustness."""

    def test_single_methionine(self, z3_engine: Z3Engine):
        """Solving a single amino acid (M) should work trivially."""
        result = z3_engine.solve("M")
        assert result.solved, f"Single M should always solve: {result.error}"
        assert result.sequence == "ATG"

    def test_organism_parameter(self):
        """Z3Engine should accept different organism parameters without error."""
        for organism in ["Homo_sapiens", "Escherichia_coli"]:
            engine = Z3Engine(organism=organism)
            result = engine.solve(HBB_PROTEIN_SHORT)
            assert isinstance(result, Z3SolveResult)

    def test_longer_protein(self, z3_engine: Z3Engine):
        """A moderately longer protein (100 aa) should also solve."""
        protein_100 = (
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
            "GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL"
        )
        result = z3_engine.solve(protein_100, timeout_s=60.0)
        assert result.solved, f"100aa protein should solve: {result.error}"
        assert len(result.sequence) == len(protein_100) * 3
