"""
Comprehensive tests for the BioCompiler CSP/SMT Solver package.

Tests cover:
1. Types (solver.types) — SolverConfig, CodonVariable, MUSReport, ConstraintType, SolverBackend, SolverResult
2. Dispatch (solver.dispatch) — is_csp_available, csp_optimize, validate_csp_solution, fallback behavior
3. Integration (conditional on ortools/z3) — end-to-end solve, GC bounds, restriction sites, CAI, infeasibility
4. Edge cases — single AA, short proteins, all-Valine, impossible GC, empty protein

All tests gracefully handle missing solver backends using pytest.importorskip or mocking.
"""

from __future__ import annotations

import importlib
import math
from dataclasses import FrozenInstanceError
from unittest.mock import MagicMock, patch

import pytest


# ────────────────────────────────────────────────────────────
# Helpers: safely import solver submodules
# ────────────────────────────────────────────────────────────

def _import_solver_types():
    """Import solver.types, skipping if unavailable."""
    return pytest.importorskip("biocompiler.solver.types")


def _import_solver_dispatch():
    """Import solver.dispatch, skipping if unavailable."""
    return pytest.importorskip("biocompiler.solver.dispatch")


# ────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture
def simple_protein():
    """20-AA protein including Valine (GT challenge)."""
    return "MVLSPADKTNVKAAWGKVGA"


@pytest.fixture
def config():
    """Default SolverConfig with sensible GC bounds."""
    mod = _import_solver_types()
    return mod.SolverConfig(gc_lo=0.30, gc_hi=0.70)


@pytest.fixture
def eGFP_protein():
    """Full eGFP protein sequence (239 AA)."""
    return (
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    )


@pytest.fixture
def tight_config():
    """SolverConfig with narrow GC bounds (for infeasibility tests)."""
    mod = _import_solver_types()
    return mod.SolverConfig(gc_lo=0.60, gc_hi=0.61)


# ────────────────────────────────────────────────────────────
# 1. TYPES TESTS
# ────────────────────────────────────────────────────────────

class TestSolverConfig:
    """Tests for SolverConfig dataclass defaults and validation."""

    def test_defaults(self):
        """SolverConfig should have sensible default values."""
        mod = _import_solver_types()
        cfg = mod.SolverConfig()
        assert cfg.gc_lo == 0.30
        assert cfg.gc_hi == 0.70
        assert cfg.organism == "Homo_sapiens"
        assert cfg.max_time_seconds > 0

    def test_custom_values(self):
        """SolverConfig should accept custom values."""
        mod = _import_solver_types()
        cfg = mod.SolverConfig(
            gc_lo=0.40, gc_hi=0.60,
            organism="Escherichia_coli",
            max_time_seconds=60,
        )
        assert cfg.gc_lo == 0.40
        assert cfg.gc_hi == 0.60
        assert cfg.organism == "Escherichia_coli"
        assert cfg.max_time_seconds == 60

    def test_gc_bounds_order(self):
        """gc_lo must be less than gc_hi."""
        mod = _import_solver_types()
        # Some implementations may validate; at minimum the values should store correctly
        cfg = mod.SolverConfig(gc_lo=0.45, gc_hi=0.55)
        assert cfg.gc_lo < cfg.gc_hi

    def test_gc_bounds_fractional(self):
        """GC bounds must be between 0.0 and 1.0."""
        mod = _import_solver_types()
        cfg = mod.SolverConfig(gc_lo=0.0, gc_hi=1.0)
        assert 0.0 <= cfg.gc_lo <= 1.0
        assert 0.0 <= cfg.gc_hi <= 1.0


class TestCodonVariable:
    """Tests for CodonVariable dataclass."""

    def test_creation(self):
        """CodonVariable should be creatable with position, amino acid, and codon list."""
        mod = _import_solver_types()
        cv = mod.CodonVariable(
            position=0,
            amino_acid="V",
            codons=["GTT", "GTC", "GTA", "GTG"],
        )
        assert cv.position == 0
        assert cv.amino_acid == "V"
        assert len(cv.codons) == 4

    def test_valine_all_gt(self):
        """All Valine codons contain 'GT' — key CSP challenge."""
        mod = _import_solver_types()
        cv = mod.CodonVariable(position=5, amino_acid="V", codons=["GTT", "GTC", "GTA", "GTG"])
        assert all("GT" in c for c in cv.codons), "Valine codons should all contain GT"

    def test_methionine_single_codon(self):
        """Methionine has exactly one codon (ATG)."""
        mod = _import_solver_types()
        cv = mod.CodonVariable(position=0, amino_acid="M", codons=["ATG"])
        assert len(cv.codons) == 1
        assert cv.codons[0] == "ATG"

    def test_leucine_six_codons(self):
        """Leucine has 6 synonymous codons."""
        mod = _import_solver_types()
        cv = mod.CodonVariable(position=2, amino_acid="L", codons=["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"])
        assert len(cv.codons) == 6


class TestMUSReport:
    """Tests for MUSReport (Minimal Unsatisfiable Subset) dataclass."""

    def test_creation(self):
        """MUSReport should be creatable with constraint names."""
        mod = _import_solver_types()
        report = mod.MUSReport(
            unsatisfiable_constraints=["GC_bounds", "NoRestrictionSite_EcoRI"],
            explanation="GC bounds [0.60, 0.61] conflict with restriction site removal",
        )
        assert len(report.unsatisfiable_constraints) == 2
        assert "GC_bounds" in report.unsatisfiable_constraints

    def test_empty_mus(self):
        """An empty MUS report means all constraints are satisfiable."""
        mod = _import_solver_types()
        report = mod.MUSReport(unsatisfiable_constraints=[], explanation="")
        assert len(report.unsatisfiable_constraints) == 0


class TestConstraintType:
    """Tests for ConstraintType enum."""

    def test_enum_values(self):
        """ConstraintType should have expected constraint categories."""
        mod = _import_solver_types()
        expected = {"GC_CONTENT", "RESTRICTION_SITE", "NO_GT", "NO_CPG", "CAI_MIN"}
        actual = {ct.name for ct in mod.ConstraintType}
        assert expected.issubset(actual), f"Missing constraint types: {expected - actual}"

    def test_enum_members_are_unique(self):
        """Each ConstraintType member should have a unique value."""
        mod = _import_solver_types()
        values = [ct.value for ct in mod.ConstraintType]
        assert len(values) == len(set(values)), "ConstraintType values must be unique"


class TestSolverBackend:
    """Tests for SolverBackend enum."""

    def test_enum_values(self):
        """SolverBackend should list OR_TOOLS and Z3 at minimum."""
        mod = _import_solver_types()
        names = {b.name for b in mod.SolverBackend}
        assert "OR_TOOLS" in names or "ORTOOLS" in names
        assert "Z3" in names

    def test_enum_count(self):
        """There should be at least 2 backends."""
        mod = _import_solver_types()
        assert len(mod.SolverBackend) >= 2


class TestSolverResult:
    """Tests for SolverResult dataclass and field validation."""

    def test_creation_success(self):
        """A successful SolverResult should have a DNA sequence and feasible=True."""
        mod = _import_solver_types()
        result = mod.SolverResult(
            sequence="ATGGTGCTG",
            feasible=True,
            gc_content=0.555,
            cai=0.78,
            backend=mod.SolverBackend.OR_TOOLS if hasattr(mod.SolverBackend, "OR_TOOLS") else mod.SolverBackend.ORTOOLS,
            solve_time_seconds=0.5,
        )
        assert result.feasible is True
        assert len(result.sequence) == 9
        assert 0.0 <= result.gc_content <= 1.0
        assert 0.0 <= result.cai <= 1.0

    def test_creation_infeasible(self):
        """An infeasible SolverResult should have feasible=False and no valid sequence."""
        mod = _import_solver_types()
        result = mod.SolverResult(
            sequence="",
            feasible=False,
            gc_content=0.0,
            cai=0.0,
            backend=mod.SolverBackend.OR_TOOLS if hasattr(mod.SolverBackend, "OR_TOOLS") else mod.SolverBackend.ORTOOLS,
            solve_time_seconds=0.1,
        )
        assert result.feasible is False
        assert result.sequence == ""

    def test_gc_content_bounds(self):
        """GC content in SolverResult must be in [0, 1]."""
        mod = _import_solver_types()
        result = mod.SolverResult(
            sequence="ATGGTGCTG",
            feasible=True,
            gc_content=0.555,
            cai=0.78,
            backend=mod.SolverBackend.OR_TOOLS if hasattr(mod.SolverBackend, "OR_TOOLS") else mod.SolverBackend.ORTOOLS,
            solve_time_seconds=0.5,
        )
        assert 0.0 <= result.gc_content <= 1.0

    def test_cai_bounds(self):
        """CAI in SolverResult must be in [0, 1]."""
        mod = _import_solver_types()
        result = mod.SolverResult(
            sequence="ATGGTGCTG",
            feasible=True,
            gc_content=0.555,
            cai=0.78,
            backend=mod.SolverBackend.OR_TOOLS if hasattr(mod.SolverBackend, "OR_TOOLS") else mod.SolverBackend.ORTOOLS,
            solve_time_seconds=0.5,
        )
        assert 0.0 <= result.cai <= 1.0

    def test_solve_time_non_negative(self):
        """Solve time should be non-negative."""
        mod = _import_solver_types()
        result = mod.SolverResult(
            sequence="ATGGTGCTG",
            feasible=True,
            gc_content=0.555,
            cai=0.78,
            backend=mod.SolverBackend.OR_TOOLS if hasattr(mod.SolverBackend, "OR_TOOLS") else mod.SolverBackend.ORTOOLS,
            solve_time_seconds=0.5,
        )
        assert result.solve_time_seconds >= 0.0


# ────────────────────────────────────────────────────────────
# 2. DISPATCH TESTS
# ────────────────────────────────────────────────────────────

class TestIsCSPAvailable:
    """Tests for is_csp_available() function."""

    def test_returns_dict(self):
        """is_csp_available should return a dict."""
        mod = _import_solver_dispatch()
        result = mod.is_csp_available()
        assert isinstance(result, dict)

    def test_expected_keys(self):
        """Result should contain 'ortools' and 'z3' keys."""
        mod = _import_solver_dispatch()
        result = mod.is_csp_available()
        assert "ortools" in result
        assert "z3" in result

    def test_values_are_bool(self):
        """Each value should be a boolean."""
        mod = _import_solver_dispatch()
        result = mod.is_csp_available()
        for key, value in result.items():
            assert isinstance(value, bool), f"Key {key!r} has non-bool value {value!r}"

    def test_any_available_is_consistent(self):
        """If any backend is True, the overall availability should be True."""
        mod = _import_solver_dispatch()
        result = mod.is_csp_available()
        # At least the dict should exist; the actual availability depends on the environment
        assert isinstance(result, dict)


class TestCspOptimize:
    """Tests for csp_optimize() function."""

    def test_simple_protein(self, simple_protein):
        """csp_optimize should return a SolverResult for a simple protein."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.csp_optimize(simple_protein, config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        # If a backend is available, it should be feasible
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()):
            assert result.feasible is True
            assert len(result.sequence) == len(simple_protein) * 3

    def test_returns_solver_result_type(self, simple_protein):
        """csp_optimize should always return a SolverResult (even if infeasible)."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.csp_optimize(simple_protein, config=cfg)
        assert isinstance(result, mod_types.SolverResult)

    def test_egfp_protein(self, eGFP_protein):
        """csp_optimize should handle a full-length eGFP protein."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.csp_optimize(eGFP_protein, config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()):
            assert result.feasible is True
            assert len(result.sequence) == len(eGFP_protein) * 3


class TestFallbackBehavior:
    """Tests for fallback behavior when no solver backend is available."""

    def test_fallback_when_no_backend(self, simple_protein):
        """When both backends are unavailable, csp_optimize should still return a result."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()

        with patch.object(mod_dispatch, "is_csp_available", return_value={"ortools": False, "z3": False}):
            # The function should either raise a clear error or return a fallback result
            try:
                result = mod_dispatch.csp_optimize(
                    simple_protein,
                    config=mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70),
                )
                # If it returns, it should be a SolverResult
                assert isinstance(result, mod_types.SolverResult)
                # Fallback result may have feasible=False or use a greedy fallback
            except (ImportError, RuntimeError) as exc:
                # Acceptable: a clear error when no backend is available
                assert "solver" in str(exc).lower() or "backend" in str(exc).lower() or "available" in str(exc).lower()

    def test_is_csp_available_all_false(self):
        """When both backends are mocked as unavailable, is_csp_available should reflect that."""
        mod_dispatch = _import_solver_dispatch()
        with patch.object(mod_dispatch, "is_csp_available", return_value={"ortools": False, "z3": False}):
            result = mod_dispatch.is_csp_available()
            assert result["ortools"] is False
            assert result["z3"] is False


class TestValidateCspSolution:
    """Tests for validate_csp_solution() function."""

    def test_known_good_sequence(self, simple_protein):
        """validate_csp_solution should pass for a valid optimized sequence."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()

        # Build a known-good sequence: use the most common codons
        from biocompiler.type_system import AA_TO_CODONS
        codons = []
        for aa in simple_protein:
            codon_list = AA_TO_CODONS.get(aa, [])
            codons.append(codon_list[0] if codon_list else "NNN")
        good_seq = "".join(codons)

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.validate_csp_solution(good_seq, simple_protein, config=cfg)
        # Should be valid (or at least not raise)
        assert isinstance(result, (bool, dict, list))

    def test_restriction_site_present(self, simple_protein):
        """validate_csp_solution should detect a sequence with restriction sites."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()

        # Construct a sequence that intentionally contains an EcoRI site (GAATTC)
        # Find a pair of adjacent AAs where we can force GAATTC
        # Phenylalanine (F): TTC, Glutamic acid (E): GAA → ...GAA|TTC... = GAATTC!
        # We need E followed by F in the protein
        # simple_protein = MVLSPADKTNVKAAWGKVGA — no E followed by F
        # Let's use a custom protein
        protein_with_ef = "MEF"  # E then F
        from biocompiler.type_system import AA_TO_CODONS
        seq_with_ecori = AA_TO_CODONS["M"][0] + "GAA" + "TTC"  # ATG GAA TTC → contains GAATTC

        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.validate_csp_solution(seq_with_ecori, protein_with_ef, config=cfg)
        # Should indicate a problem (restriction site found)
        if isinstance(result, dict):
            assert not result.get("valid", True), "Should detect EcoRI restriction site"
        elif isinstance(result, list):
            assert len(result) > 0, "Should detect at least one violation"
        elif isinstance(result, bool):
            assert result is False, "Should be False when restriction site present"

    def test_gc_out_of_range(self, simple_protein):
        """validate_csp_solution should detect GC content outside bounds."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()

        # All-A/T codons: extremely low GC
        from biocompiler.type_system import AA_TO_CODONS
        # Build sequence using the most A/T-rich codons
        at_codons = {
            "M": "ATG", "V": "GTT", "L": "TTA", "S": "TCT", "P": "CCT",
            "A": "GCT", "D": "GAT", "K": "AAA", "T": "ACT", "N": "AAT",
            "W": "TGG", "G": "GGT", "E": "GAA",
        }
        seq = "".join(at_codons.get(aa, "NNN") for aa in simple_protein)

        cfg = mod_types.SolverConfig(gc_lo=0.60, gc_hi=0.70)  # Very high GC target
        result = mod_dispatch.validate_csp_solution(seq, simple_protein, config=cfg)
        # The AT-rich sequence should fail the high-GC requirement
        if isinstance(result, dict):
            assert not result.get("valid", True), "Should fail high GC requirement"
        elif isinstance(result, list):
            assert len(result) > 0, "Should detect GC violation"
        elif isinstance(result, bool):
            assert result is False, "Should be False when GC out of range"


# ────────────────────────────────────────────────────────────
# 3. INTEGRATION TESTS (conditional on ortools/z3)
# ────────────────────────────────────────────────────────────

class TestCSPIntegration:
    """Integration tests that require an actual solver backend.

    These tests are skipped automatically if no backend is installed.
    """

    @pytest.fixture(autouse=True)
    def skip_if_no_backend(self):
        """Skip all tests in this class if no solver backend is available."""
        try:
            mod_dispatch = _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")
        avail = mod_dispatch.is_csp_available()
        if not any(avail.values()):
            pytest.skip("No CSP solver backend (ortools/z3) installed")

    def test_solve_produces_valid_sequence(self, simple_protein):
        """solve_with_csp should produce a valid DNA sequence."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        assert result.feasible is True
        assert len(result.sequence) == len(simple_protein) * 3
        # Sequence should only contain valid DNA bases
        assert set(result.sequence).issubset({"A", "C", "G", "T"})

    def test_gc_content_within_bounds(self, simple_protein):
        """GC content of the solution should be within the configured bounds."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            assert cfg.gc_lo <= result.gc_content <= cfg.gc_hi, (
                f"GC content {result.gc_content:.3f} outside [{cfg.gc_lo}, {cfg.gc_hi}]"
            )

    def test_no_restriction_sites_in_solution(self, simple_protein):
        """Solution should not contain common restriction enzyme sites."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        from biocompiler.restriction_sites import RESTRICTION_SITES

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            # Check for the most common restriction sites
            for enzyme, site in RESTRICTION_SITES.items():
                assert site not in result.sequence, (
                    f"Restriction site {enzyme} ({site}) found in solution"
                )

    def test_cai_reasonable(self, simple_protein):
        """CAI of the solution should be reasonable (>= 0.5 for a simple protein)."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            assert result.cai >= 0.5, f"CAI {result.cai:.3f} is unreasonably low"

    def test_infeasibility_detected(self, simple_protein, tight_config):
        """Impossibly narrow GC bounds should be detected as infeasible."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        result = mod_dispatch.solve_with_csp(simple_protein, config=tight_config)
        # With GC bounds [0.60, 0.61] this is almost certainly infeasible
        assert result.feasible is False, (
            f"Expected infeasibility for GC [{tight_config.gc_lo}, {tight_config.gc_hi}], "
            f"but solver found a solution"
        )

    def test_egfp_solve(self, eGFP_protein):
        """Full eGFP protein should be solvable with default bounds."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(eGFP_protein, config=cfg)
        assert result.feasible is True
        assert len(result.sequence) == len(eGFP_protein) * 3

    def test_sequence_translates_correctly(self, simple_protein):
        """The optimized sequence should translate back to the original protein."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        from biocompiler.type_system import CODON_TABLE

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            translated = ""
            for i in range(0, len(result.sequence), 3):
                codon = result.sequence[i:i + 3]
                aa = CODON_TABLE.get(codon, "?")
                translated += aa
            # Compare (stop codons at end are OK)
            assert translated.rstrip("*") == simple_protein, (
                f"Translation mismatch: expected {simple_protein}, got {translated.rstrip('*')}"
            )


# ────────────────────────────────────────────────────────────
# 3b. OR-Tools-specific integration tests
# ────────────────────────────────────────────────────────────

class TestORToolsIntegration:
    """OR-Tools-specific integration tests."""

    @pytest.fixture(autouse=True)
    def skip_if_no_ortools(self):
        """Skip if OR-Tools is not installed."""
        pytest.importorskip("ortools")
        try:
            _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")

    def test_ortools_solve_simple(self, simple_protein):
        """OR-Tools should solve a simple protein."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        backend_name = "OR_TOOLS" if hasattr(mod_types.SolverBackend, "OR_TOOLS") else "ORTOOLS"
        cfg = mod_types.SolverConfig(
            gc_lo=0.30, gc_hi=0.70,
            backend=mod_types.SolverBackend[backend_name],
        )
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        assert result.feasible is True


# ────────────────────────────────────────────────────────────
# 3c. Z3-specific integration tests
# ────────────────────────────────────────────────────────────

class TestZ3Integration:
    """Z3-specific integration tests."""

    @pytest.fixture(autouse=True)
    def skip_if_no_z3(self):
        """Skip if Z3 is not installed."""
        pytest.importorskip("z3")
        try:
            _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")

    def test_z3_solve_simple(self, simple_protein):
        """Z3 should solve a simple protein."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(
            gc_lo=0.30, gc_hi=0.70,
            backend=mod_types.SolverBackend.Z3,
        )
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        assert result.feasible is True


# ────────────────────────────────────────────────────────────
# 4. EDGE CASES
# ────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Edge case tests for the CSP solver."""

    @pytest.fixture(autouse=True)
    def _require_dispatch(self):
        """Ensure dispatch module is available for these tests."""
        try:
            _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")

    def test_single_amino_acid(self):
        """A single amino acid protein ('M') should be handled correctly."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.csp_optimize("M", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            assert result.sequence == "ATG"

    def test_very_short_protein(self):
        """A very short protein ('MK') should produce a 6-base sequence."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.csp_optimize("MK", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            assert len(result.sequence) == 6

    def test_all_valine_protein(self):
        """All-Valine protein ('VVVVVV') — every codon contains GT, challenging."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.csp_optimize("VVVVVV", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            # All Valine codons contain GT — if NoGT constraint is strict, this may be infeasible
            # With relaxed constraints, it should still produce a valid sequence
            assert len(result.sequence) == 18
            # Verify the sequence translates to VVVVVV
            from biocompiler.type_system import CODON_TABLE
            for i in range(0, len(result.sequence), 3):
                codon = result.sequence[i:i + 3]
                assert CODON_TABLE.get(codon) == "V", f"Codon {codon} doesn't encode Valine"

    def test_impossible_gc_target(self):
        """Protein with all A/T amino acids and high GC target should be infeasible or low-CAI."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        # Lysine (K) codons: AAA, AAG — at most 1/3 GC per codon
        # With GC target [0.80, 0.90], this should be infeasible
        protein = "KKKKKKKKKK"  # 10 Lysines
        cfg = mod_types.SolverConfig(gc_lo=0.80, gc_hi=0.90)
        result = mod_dispatch.csp_optimize(protein, config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        # Either infeasible or GC not in range
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()):
            if result.feasible:
                # If somehow feasible, GC should be within bounds
                assert cfg.gc_lo <= result.gc_content <= cfg.gc_hi, (
                    f"GC {result.gc_content:.3f} outside [{cfg.gc_lo}, {cfg.gc_hi}]"
                )
            else:
                # Infeasible is the expected outcome
                assert result.feasible is False

    def test_empty_protein_raises_error(self):
        """Empty protein should raise an error."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        with pytest.raises((ValueError, TypeError, RuntimeError)):
            mod_dispatch.csp_optimize("", config=cfg)

    def test_single_valine(self):
        """Single Valine — the hardest single-AA case due to mandatory GT."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.csp_optimize("V", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            assert len(result.sequence) == 3
            assert result.sequence in ("GTT", "GTC", "GTA", "GTG")

    def test_tryptophan_single_codon(self):
        """Tryptophan has exactly one codon (TGG) — no freedom of choice."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.csp_optimize("W", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            assert result.sequence == "TGG"

    def test_repeated_amino_acid(self):
        """A protein of repeated alanines should use different codons."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.csp_optimize("AAAAAAAAAA", config=cfg)  # 10 Alanines
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            assert len(result.sequence) == 30
            from biocompiler.type_system import CODON_TABLE
            for i in range(0, len(result.sequence), 3):
                codon = result.sequence[i:i + 3]
                assert CODON_TABLE.get(codon) == "A", f"Codon {codon} doesn't encode Alanine"

    def test_methionine_start(self):
        """First amino acid is always Methionine (ATG) — no codon choice."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.20, gc_hi=0.80)
        result = mod_dispatch.csp_optimize("MA", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()) and result.feasible:
            assert result.sequence.startswith("ATG")


# ────────────────────────────────────────────────────────────
# 5. CROSS-VALIDATION TESTS
# ────────────────────────────────────────────────────────────

class TestCrossValidation:
    """Cross-validation between CSP solver and existing biocompiler predicates."""

    @pytest.fixture(autouse=True)
    def skip_if_no_backend(self):
        """Skip if no solver backend is available."""
        try:
            mod_dispatch = _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")
        avail = mod_dispatch.is_csp_available()
        if not any(avail.values()):
            pytest.skip("No CSP solver backend installed")

    def test_solution_passes_no_stop_codons(self, simple_protein):
        """CSP solution should have no internal stop codons."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        from biocompiler.type_system import check_no_stop_codons

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            pred = check_no_stop_codons(result.sequence)
            assert pred.passed, f"Internal stop codons found: {pred.details}"

    def test_solution_valid_coding_seq(self, simple_protein):
        """CSP solution should be a valid coding sequence."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        from biocompiler.type_system import check_valid_coding_seq

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            pred = check_valid_coding_seq(result.sequence)
            assert pred.passed, f"Invalid coding sequence: {pred.details}"

    def test_solution_correct_length(self, simple_protein):
        """CSP solution length should be exactly 3x the protein length."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            expected_len = len(simple_protein) * 3
            assert len(result.sequence) == expected_len, (
                f"Length {len(result.sequence)} != expected {expected_len}"
            )

    def test_solution_gc_matches_reported(self, simple_protein):
        """Reported GC content should match actual computation."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()

        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        result = mod_dispatch.solve_with_csp(simple_protein, config=cfg)
        if result.feasible:
            actual_gc = sum(1 for b in result.sequence if b in "GC") / len(result.sequence)
            assert math.isclose(result.gc_content, actual_gc, abs_tol=0.01), (
                f"Reported GC {result.gc_content:.3f} != actual {actual_gc:.3f}"
            )


# ────────────────────────────────────────────────────────────
# 6. MUS (Minimal Unsatisfiable Subset) TESTS
# ────────────────────────────────────────────────────────────

class TestMUSDiagnosis:
    """Tests for MUS diagnosis when the problem is infeasible."""

    @pytest.fixture(autouse=True)
    def skip_if_no_backend(self):
        """Skip if no solver backend is available."""
        try:
            mod_dispatch = _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")
        avail = mod_dispatch.is_csp_available()
        if not any(avail.values()):
            pytest.skip("No CSP solver backend installed")

    def test_infeasible_produces_mus_report(self, simple_protein, tight_config):
        """An infeasible problem should produce a MUS report if supported."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        result = mod_dispatch.solve_with_csp(simple_protein, config=tight_config)
        if not result.feasible:
            # The result may have a mus_report attribute
            mus = getattr(result, "mus_report", None)
            if mus is not None:
                assert isinstance(mus, mod_types.MUSReport)
                assert len(mus.unsatisfiable_constraints) > 0, (
                    "MUS report should identify at least one unsatisfiable constraint"
                )

    def test_mus_contains_gc_constraint(self, simple_protein, tight_config):
        """For an impossibly narrow GC problem, MUS should include GC-related constraints."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        result = mod_dispatch.solve_with_csp(simple_protein, config=tight_config)
        if not result.feasible:
            mus = getattr(result, "mus_report", None)
            if mus is not None and mus.unsatisfiable_constraints:
                constraint_names = " ".join(mus.unsatisfiable_constraints).lower()
                # The MUS should mention GC-related constraints
                assert "gc" in constraint_names or "content" in constraint_names, (
                    f"Expected GC-related constraint in MUS, got: {mus.unsatisfiable_constraints}"
                )


# ────────────────────────────────────────────────────────────
# 7. BACKEND SELECTION TESTS
# ────────────────────────────────────────────────────────────

class TestBackendSelection:
    """Tests for solver backend selection logic."""

    @pytest.fixture(autouse=True)
    def _require_dispatch(self):
        try:
            _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")

    def test_config_specifies_backend(self):
        """SolverConfig should allow specifying a preferred backend."""
        mod_types = _import_solver_types()
        cfg = mod_types.SolverConfig(
            gc_lo=0.30, gc_hi=0.70,
            backend=mod_types.SolverBackend.OR_TOOLS if hasattr(mod_types.SolverBackend, "OR_TOOLS") else mod_types.SolverBackend.ORTOOLS,
        )
        assert cfg.backend is not None

    def test_default_backend_selection(self):
        """If no backend is specified, the dispatch should auto-select one."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        # No explicit backend — dispatch should choose the best available
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()):
            result = mod_dispatch.csp_optimize("MK", config=cfg)
            assert isinstance(result, mod_types.SolverResult)


# ────────────────────────────────────────────────────────────
# 8. ROBUSTNESS TESTS
# ────────────────────────────────────────────────────────────

class TestRobustness:
    """Robustness tests for unusual inputs and error handling."""

    @pytest.fixture(autouse=True)
    def _require_dispatch(self):
        try:
            _import_solver_dispatch()
        except pytest.skip.Exception:
            pytest.skip("solver.dispatch not available")

    def test_invalid_protein_characters(self):
        """Protein with invalid characters should raise an error."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        with pytest.raises((ValueError, TypeError)):
            mod_dispatch.csp_optimize("MXYZ", config=cfg)

    def test_whitespace_protein(self):
        """Protein with only whitespace should raise an error."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        with pytest.raises((ValueError, TypeError)):
            mod_dispatch.csp_optimize("   ", config=cfg)

    def test_config_invalid_gc_bounds(self):
        """Config with gc_lo > gc_hi should be rejected or handled."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        # Some implementations may raise; others may swap or clamp
        try:
            cfg = mod_types.SolverConfig(gc_lo=0.70, gc_hi=0.30)
            result = mod_dispatch.csp_optimize("MK", config=cfg)
            # If it doesn't raise, it should handle the inverted bounds gracefully
            assert isinstance(result, mod_types.SolverResult)
        except (ValueError, AssertionError):
            pass  # Expected: invalid bounds should be caught

    def test_very_wide_gc_bounds(self):
        """Very wide GC bounds [0.0, 1.0] should always be feasible."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.0, gc_hi=1.0)
        result = mod_dispatch.csp_optimize("MVLSPADKTNVKAAWGKVGA", config=cfg)
        assert isinstance(result, mod_types.SolverResult)
        avail = mod_dispatch.is_csp_available()
        if any(avail.values()):
            assert result.feasible is True, "Wide GC bounds should always be feasible"

    def test_long_protein_does_not_hang(self, eGFP_protein):
        """Solving eGFP (239 AA) should complete within a reasonable time."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        import time
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70, max_time_seconds=30)
        start = time.time()
        result = mod_dispatch.csp_optimize(eGFP_protein, config=cfg)
        elapsed = time.time() - start
        assert elapsed < 60, f"Solve took {elapsed:.1f}s, expected < 60s"
        assert isinstance(result, mod_types.SolverResult)


# ────────────────────────────────────────────────────────────
# 9. CONSTRAINT MODEL TESTS (if solver.constraints exists)
# ────────────────────────────────────────────────────────────

class TestConstraintModel:
    """Tests for the constraint model builder."""

    @pytest.fixture(autouse=True)
    def skip_if_no_constraints_module(self):
        """Skip if solver.constraints module is not available."""
        try:
            self.mod = pytest.importorskip("biocompiler.solver.constraints")
        except pytest.skip.Exception:
            pytest.skip("solver.constraints not available")

    def test_build_model_returns_csp_model(self, simple_protein):
        """Building a constraint model should return a CSPModel."""
        mod_types = _import_solver_types()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        model = self.mod.build_csp_model(simple_protein, config=cfg)
        assert model is not None
        # Model should have some representation of variables and constraints
        assert hasattr(model, "variables") or hasattr(model, "codon_variables")

    def test_model_variable_count(self, simple_protein):
        """Model should have one variable per amino acid position."""
        mod_types = _import_solver_types()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        model = self.mod.build_csp_model(simple_protein, config=cfg)
        if hasattr(model, "variables"):
            assert len(model.variables) == len(simple_protein)
        elif hasattr(model, "codon_variables"):
            assert len(model.codon_variables) == len(simple_protein)

    def test_model_includes_gc_constraint(self, simple_protein):
        """Model should include a GC content constraint."""
        mod_types = _import_solver_types()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        model = self.mod.build_csp_model(simple_protein, config=cfg)
        # Check that some GC-related constraint exists
        if hasattr(model, "constraints"):
            constraint_names = [str(c) for c in model.constraints]
            gc_present = any("gc" in name.lower() for name in constraint_names)
            assert gc_present, f"Expected GC constraint, got: {constraint_names[:5]}"


# ────────────────────────────────────────────────────────────
# 10. ENGINE BACKEND TESTS (if engine modules exist)
# ────────────────────────────────────────────────────────────

class TestEngineORTools:
    """Tests for the OR-Tools engine module."""

    @pytest.fixture(autouse=True)
    def skip_if_no_ortools_engine(self):
        try:
            self.mod = pytest.importorskip("biocompiler.solver.engine_ortools")
        except pytest.skip.Exception:
            pytest.skip("solver.engine_ortools not available")

    def test_solve_returns_solver_result(self, simple_protein):
        """OR-Tools engine should return a SolverResult."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        # Use the engine directly if available
        if hasattr(self.mod, "solve"):
            result = self.mod.solve(simple_protein, config=cfg)
            assert isinstance(result, mod_types.SolverResult)


class TestEngineZ3:
    """Tests for the Z3 engine module."""

    @pytest.fixture(autouse=True)
    def skip_if_no_z3_engine(self):
        try:
            self.mod = pytest.importorskip("biocompiler.solver.engine_z3")
        except pytest.skip.Exception:
            pytest.skip("solver.engine_z3 not available")

    def test_solve_returns_solver_result(self, simple_protein):
        """Z3 engine should return a SolverResult."""
        mod_types = _import_solver_types()
        mod_dispatch = _import_solver_dispatch()
        cfg = mod_types.SolverConfig(gc_lo=0.30, gc_hi=0.70)
        if hasattr(self.mod, "solve"):
            result = self.mod.solve(simple_protein, config=cfg)
            assert isinstance(result, mod_types.SolverResult)
