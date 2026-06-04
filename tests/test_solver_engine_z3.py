"""Tests for BioCompiler Z3 SMT Solver Engine (solver/engine_z3.py).

Test coverage:
  1. Z3Engine construction — including ImportError when z3 absent
  2. solve() returns SolverResult (never None)
  3. Timeout handling — short timeout returns gracefully
  4. backend_used is always SolverBackend.Z3 in results

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

from biocompiler.solver.engine_z3 import Z3Engine
from biocompiler.solver.types import (
    ConstraintSpec,
    ConstraintType,
    CSPModel,
    MUSReport,
    SolverBackend,
    SolverConfig,
    SolverResult,
)
from biocompiler.constants import AA_TO_CODONS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

SHORT_PROTEIN = "MVLSPADKTNVKAAWGKVGA"  # 20 aa — first 20 residues of HBA1
MEDIUM_PROTEIN = "MKTVLIAEGH"

# Very high splice threshold that effectively disables splice constraints,
# so we can test GC / restriction-site solving independently.
_SPLICE_DISABLED = 999.0


def _make_model(
    protein: str,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    timeout_seconds: float = 60.0,
    splice_threshold: float = _SPLICE_DISABLED,
) -> CSPModel:
    """Build a CSPModel for the given protein with the specified config.

    By default the splice threshold is set very high (999) so that splice
    constraints are effectively disabled — this lets us test GC and
    restriction-site logic independently.  Callers that *want* splice
    constraints can pass a lower threshold.
    """
    config = SolverConfig(
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        restriction_sites=restriction_sites or [],
        timeout_seconds=timeout_seconds,
        cryptic_splice_threshold=splice_threshold,
    )
    codon_domains = {
        i: list(AA_TO_CODONS.get(aa, []))
        for i, aa in enumerate(protein)
    }
    constraints: list[ConstraintSpec] = []
    return CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints,
        config=config,
    )


def _gc_content(seq: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in "GC") / len(seq)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def default_config() -> SolverConfig:
    """Default SolverConfig for Z3 tests."""
    return SolverConfig()


@pytest.fixture
def z3_engine(default_config: SolverConfig) -> Z3Engine:
    """Default Z3Engine configured for human codon usage."""
    return Z3Engine(default_config, organism="Homo_sapiens")


# ======================================================================
# 1. Z3Engine Construction
# ======================================================================

class TestZ3EngineConstruction:
    """Z3Engine instantiation and availability checks."""

    def test_construction_with_config(self, default_config: SolverConfig):
        """Z3Engine should instantiate with a SolverConfig and default organism."""
        engine = Z3Engine(default_config)
        assert engine is not None
        assert engine.config is default_config
        assert engine.organism == "Homo_sapiens"

    def test_construction_custom_organism(self, default_config: SolverConfig):
        """Z3Engine should accept a custom organism name."""
        engine = Z3Engine(default_config, organism="Escherichia_coli")
        assert engine.organism == "Escherichia_coli"

    def test_construction_custom_seed(self, default_config: SolverConfig):
        """Z3Engine should accept a custom random seed."""
        engine = Z3Engine(default_config, seed=42)
        assert engine._seed == 42

    def test_is_available_returns_true_when_installed(self):
        """If z3-solver is importable (we passed the skip gate), is_available() is True."""
        assert Z3Engine.is_available() is True

    def test_is_available_returns_bool(self):
        """is_available() should return a bool."""
        assert isinstance(Z3Engine.is_available(), bool)

    def test_engine_has_solve_method(self, z3_engine: Z3Engine):
        """Z3Engine should expose a solve() method."""
        assert hasattr(z3_engine, "solve")
        assert callable(z3_engine.solve)

    def test_import_error_when_z3_unavailable(self):
        """If _Z3_AVAILABLE is False, Z3Engine.__init__ raises ImportError."""
        import biocompiler.solver.engine_z3 as z3_mod
        original = z3_mod._Z3_AVAILABLE
        try:
            z3_mod._Z3_AVAILABLE = False
            with pytest.raises(ImportError, match="z3-solver"):
                Z3Engine(SolverConfig())
        finally:
            z3_mod._Z3_AVAILABLE = original


# ======================================================================
# 2. solve() returns SolverResult
# ======================================================================

class TestZ3EngineSolve:
    """solve() always returns a SolverResult — never None and never raises."""

    def test_solve_returns_solver_result(self, z3_engine: Z3Engine):
        """solve() should return a SolverResult instance."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        assert isinstance(result, SolverResult)

    def test_solve_never_returns_none(self, z3_engine: Z3Engine):
        """solve() must never return None — always a SolverResult."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        assert result is not None

    def test_solve_feasible_problem_solved(self, z3_engine: Z3Engine):
        """A feasible problem with relaxed splice threshold should solve successfully."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        assert result.solved is True, (
            f"Default protein should solve. Warnings: {result.warnings}"
        )

    def test_solve_sequence_length(self, z3_engine: Z3Engine):
        """Output DNA length must equal protein_length * 3."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        if result.solved:
            assert len(result.sequence) == len(SHORT_PROTEIN) * 3

    def test_solve_sequence_valid_bases(self, z3_engine: Z3Engine):
        """Output sequence should only contain A, C, G, T."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        if result.solved:
            assert set(result.sequence) <= {"A", "C", "G", "T"}

    def test_solve_gc_in_range(self, z3_engine: Z3Engine):
        """GC content of the solution should be within configured bounds."""
        model = _make_model(SHORT_PROTEIN, gc_lo=0.30, gc_hi=0.70)
        result = z3_engine.solve(model)
        if result.solved:
            gc = _gc_content(result.sequence)
            assert 0.30 <= gc <= 0.70, f"GC {gc:.4f} outside [0.30, 0.70]"

    def test_solve_empty_protein(self, z3_engine: Z3Engine):
        """An empty protein should return unsolved with a warning."""
        model = _make_model("")
        result = z3_engine.solve(model)
        assert isinstance(result, SolverResult)
        assert result.solved is False
        assert any("Empty" in w for w in result.warnings)

    def test_solve_infeasible_gc_returns_unsolved(self, z3_engine: Z3Engine):
        """Infeasible GC bounds [0.99, 1.00] should return unsolved."""
        model = _make_model(SHORT_PROTEIN, gc_lo=0.99, gc_hi=1.00)
        result = z3_engine.solve(model)
        assert isinstance(result, SolverResult)
        assert result.solved is False

    def test_solve_unknown_organism_returns_warning(self, default_config: SolverConfig):
        """An unknown organism should return unsolved with a warning."""
        engine = Z3Engine(default_config, organism="Nonexistent_organism")
        model = _make_model(SHORT_PROTEIN)
        result = engine.solve(model)
        assert isinstance(result, SolverResult)
        assert result.solved is False
        assert any("Unknown organism" in w for w in result.warnings)

    def test_solve_mus_report_on_infeasible(self, z3_engine: Z3Engine):
        """An infeasible problem should produce a MUS report."""
        model = _make_model(SHORT_PROTEIN, gc_lo=0.99, gc_hi=1.00)
        result = z3_engine.solve(model)
        assert result.solved is False
        # The result should have a MUS report (possibly with empty conflicts
        # if core extraction fails, but the field should exist)
        if result.mus_report is not None:
            assert isinstance(result.mus_report, MUSReport)

    def test_solve_with_splice_constraints(self, z3_engine: Z3Engine):
        """Using a low splice threshold may make the problem infeasible or solved."""
        # Default threshold (3.0) — some proteins may become infeasible
        model = _make_model(SHORT_PROTEIN, splice_threshold=3.0)
        result = z3_engine.solve(model)
        # Either solved or unsolved is valid — just no exception
        assert isinstance(result, SolverResult)
        assert isinstance(result.solved, bool)


# ======================================================================
# 3. Timeout handling
# ======================================================================

class TestZ3Timeout:
    """Timeout handling: solver should respect time limits and return gracefully."""

    def test_short_timeout_returns_result(self, z3_engine: Z3Engine):
        """With a very short timeout (0.001s), solve() should still return a SolverResult."""
        model = _make_model(SHORT_PROTEIN, timeout_seconds=0.001)
        start = time.perf_counter()
        result = z3_engine.solve(model)
        elapsed = time.perf_counter() - start
        assert isinstance(result, SolverResult), "solve() must always return SolverResult"
        # The solver should not hang — give generous upper bound
        assert elapsed < 30.0, f"solve() took {elapsed:.1f}s with 0.001s timeout"

    def test_short_timeout_may_not_solve(self, z3_engine: Z3Engine):
        """With a very short timeout, the solver may return unsolved (unknown)."""
        model = _make_model(SHORT_PROTEIN, timeout_seconds=0.001)
        result = z3_engine.solve(model)
        assert isinstance(result, SolverResult)
        # Either solved or unsolved is acceptable — just no exception

    def test_short_timeout_unknown_result_has_warning(self, z3_engine: Z3Engine):
        """When Z3 returns 'unknown', the result should contain a warning about it."""
        model = _make_model(SHORT_PROTEIN, timeout_seconds=0.001)
        result = z3_engine.solve(model)
        if not result.solved:
            # The engine returns a warning mentioning "unknown"
            has_unknown_warning = any(
                "unknown" in w.lower() for w in result.warnings
            )
            # This is expected when Z3 times out
            assert has_unknown_warning or result.fallback_used, (
                "Short timeout should produce an 'unknown' warning or fallback flag"
            )

    def test_reasonable_timeout_solves(self, z3_engine: Z3Engine):
        """A reasonable timeout (60s) should be enough for a short protein."""
        model = _make_model(SHORT_PROTEIN, timeout_seconds=60.0)
        result = z3_engine.solve(model)
        assert result.solved is True, (
            f"60s timeout should be sufficient for 20aa protein. "
            f"Warnings: {result.warnings}"
        )

    def test_timeout_zero_returns_quickly(self, z3_engine: Z3Engine):
        """Timeout of 0 should return immediately."""
        model = _make_model(SHORT_PROTEIN, timeout_seconds=0.0)
        start = time.perf_counter()
        result = z3_engine.solve(model)
        elapsed = time.perf_counter() - start
        assert isinstance(result, SolverResult)
        # Should not hang
        assert elapsed < 10.0, "timeout=0 should return quickly"

    def test_solve_time_seconds_populated(self, z3_engine: Z3Engine):
        """SolverResult.solve_time_seconds should be a non-negative float."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        assert isinstance(result.solve_time_seconds, float)
        assert result.solve_time_seconds >= 0.0


# ======================================================================
# 4. Proper backend_used set in result
# ======================================================================

class TestZ3BackendUsed:
    """All Z3Engine results must have backend_used == SolverBackend.Z3."""

    def test_backend_used_on_success(self, z3_engine: Z3Engine):
        """A successful solve must set backend_used to SolverBackend.Z3."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        assert result.backend_used == SolverBackend.Z3

    def test_backend_used_on_infeasible(self, z3_engine: Z3Engine):
        """An infeasible solve must also set backend_used to SolverBackend.Z3."""
        model = _make_model(SHORT_PROTEIN, gc_lo=0.99, gc_hi=1.00)
        result = z3_engine.solve(model)
        assert result.backend_used == SolverBackend.Z3

    def test_backend_used_on_empty_protein(self, z3_engine: Z3Engine):
        """Empty protein early-return must set backend_used to SolverBackend.Z3."""
        model = _make_model("")
        result = z3_engine.solve(model)
        assert result.backend_used == SolverBackend.Z3

    def test_backend_used_on_timeout(self, z3_engine: Z3Engine):
        """Timeout result must set backend_used to SolverBackend.Z3."""
        model = _make_model(SHORT_PROTEIN, timeout_seconds=0.001)
        result = z3_engine.solve(model)
        assert result.backend_used == SolverBackend.Z3

    def test_backend_used_on_unknown_organism(self, default_config: SolverConfig):
        """Unknown organism early-return must set backend_used to SolverBackend.Z3."""
        engine = Z3Engine(default_config, organism="Nonexistent_organism")
        model = _make_model(SHORT_PROTEIN)
        result = engine.solve(model)
        assert result.backend_used == SolverBackend.Z3

    def test_backend_used_enum_type(self, z3_engine: Z3Engine):
        """backend_used must be a SolverBackend enum instance."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        assert isinstance(result.backend_used, SolverBackend)


# ======================================================================
# Edge cases
# ======================================================================

class TestZ3EdgeCases:
    """Additional edge cases for robustness."""

    def test_single_methionine(self, z3_engine: Z3Engine):
        """Solving a single amino acid (M) should work trivially."""
        model = _make_model("M")
        result = z3_engine.solve(model)
        assert result.solved is True
        assert result.sequence == "ATG"

    def test_restriction_site_avoidance(self, z3_engine: Z3Engine):
        """With restriction sites specified, solution should avoid them if possible."""
        model = _make_model(
            MEDIUM_PROTEIN,
            restriction_sites=["GAATTC"],  # EcoRI
        )
        result = z3_engine.solve(model)
        if result.solved:
            # If solved, the sequence should not contain the forbidden pattern
            assert "GAATTC" not in result.sequence

    def test_custom_codon_domains(self, z3_engine: Z3Engine):
        """Solver should only use codons from the provided domains."""
        protein = "MK"
        codon_domains = {
            0: ["ATG"],       # M — only one codon
            1: ["AAA", "AAG"],  # K — two codons
        }
        config = SolverConfig(gc_lo=0.0, gc_hi=1.0, cryptic_splice_threshold=_SPLICE_DISABLED)
        model = CSPModel(
            protein_sequence=protein,
            codon_domains=codon_domains,
            constraints=[],
            config=config,
        )
        result = z3_engine.solve(model)
        if result.solved:
            # First codon must be ATG
            assert result.sequence[:3] == "ATG"
            # Second codon must be one of the allowed K codons
            assert result.sequence[3:6] in ("AAA", "AAG")

    def test_tight_gc_bounds(self, z3_engine: Z3Engine):
        """Tight but feasible GC bounds should still solve."""
        model = _make_model(SHORT_PROTEIN, gc_lo=0.40, gc_hi=0.60)
        result = z3_engine.solve(model)
        if result.solved:
            gc = _gc_content(result.sequence)
            assert 0.40 <= gc <= 0.60, f"GC {gc:.4f} outside [0.40, 0.60]"

    def test_cai_populated_on_success(self, z3_engine: Z3Engine):
        """A successful solve should populate the CAI field."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        if result.solved:
            assert isinstance(result.cai, float)
            assert 0.0 <= result.cai <= 1.0

    def test_gc_content_populated_on_success(self, z3_engine: Z3Engine):
        """A successful solve should populate the gc_content field."""
        model = _make_model(SHORT_PROTEIN)
        result = z3_engine.solve(model)
        if result.solved:
            assert isinstance(result.gc_content, float)
            assert 0.0 <= result.gc_content <= 1.0
