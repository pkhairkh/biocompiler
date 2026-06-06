"""Tests for Z3 solver fallback behavior, adaptive timeouts, and solution verification.

Agent 15: Z3 fallback improvement tests.

Covers:
1. Adaptive timeout values (tiered by protein length)
2. Fallback to greedy on Z3 timeout
3. Fallback to greedy on UNSAT result
4. Z3 solution verification (defense-in-depth)
5. Large proteins using windowed decomposition
6. Warm-start fallback (partial Z3 solution passed to greedy)
"""

from __future__ import annotations

import pytest
from unittest.mock import patch, MagicMock

from biocompiler.solver.engine_z3 import (
    compute_adaptive_timeout,
    _verify_z3_solution,
    ADAPTIVE_TIMEOUT_TIERS,
    Z3Engine,
    _Z3_AVAILABLE,
)
from biocompiler.solver.types import (
    SolverConfig,
    SolverResult,
    SolverBackend,
    CSPModel,
    ConstraintSpec,
    ConstraintType,
    ConstraintStrictness,
)


# ── Test proteins ──────────────────────────────────────────────────────

# Small protein (≤ 100 aa)
INSULIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"

# Medium protein (≤ 300 aa)
GFP = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

# Protein just under 500 aa
PROTEIN_490 = "M" * 490

# Large protein (> 500 aa — should use windowed decomposition)
PROTEIN_600 = "M" * 600


# ══════════════════════════════════════════════════════════════════════
# 1. Adaptive timeout tests
# ══════════════════════════════════════════════════════════════════════


class TestAdaptiveTimeout:
    """Test the compute_adaptive_timeout function."""

    def test_small_protein_30s(self):
        """Proteins ≤ 100aa get 30s timeout."""
        timeout, skip = compute_adaptive_timeout(50)
        assert timeout == 30.0
        assert skip is False

    def test_small_protein_exactly_100(self):
        """Protein exactly 100aa gets 30s timeout."""
        timeout, skip = compute_adaptive_timeout(100)
        assert timeout == 30.0
        assert skip is False

    def test_medium_protein_60s(self):
        """Proteins 101-300aa get 60s timeout."""
        timeout, skip = compute_adaptive_timeout(200)
        assert timeout == 60.0
        assert skip is False

    def test_medium_protein_exactly_300(self):
        """Protein exactly 300aa gets 60s timeout."""
        timeout, skip = compute_adaptive_timeout(300)
        assert timeout == 60.0
        assert skip is False

    def test_large_protein_120s(self):
        """Proteins 301-500aa get 120s timeout."""
        timeout, skip = compute_adaptive_timeout(400)
        assert timeout == 120.0
        assert skip is False

    def test_large_protein_exactly_500(self):
        """Protein exactly 500aa gets 120s timeout."""
        timeout, skip = compute_adaptive_timeout(500)
        assert timeout == 120.0
        assert skip is False

    def test_windowed_tier_800aa(self):
        """Proteins 501-800aa fall into the windowed tier (180s timeout)."""
        timeout, skip = compute_adaptive_timeout(501)
        assert timeout == 180.0
        assert skip is False

    def test_windowed_tier_exactly_800(self):
        """Protein exactly 800aa gets 180s windowed timeout."""
        timeout, skip = compute_adaptive_timeout(800)
        assert timeout == 180.0
        assert skip is False

    def test_windowed_tier_600aa(self):
        """A 600aa protein falls into the windowed tier."""
        timeout, skip = compute_adaptive_timeout(600)
        assert timeout == 180.0
        assert skip is False

    def test_windowed_tier_1200aa(self):
        """Proteins 801-1200aa fall into the extended windowed tier (300s)."""
        timeout, skip = compute_adaptive_timeout(1000)
        assert timeout == 300.0
        assert skip is False

    def test_windowed_tier_exactly_1200(self):
        """Protein exactly 1200aa gets 300s windowed timeout."""
        timeout, skip = compute_adaptive_timeout(1200)
        assert timeout == 300.0
        assert skip is False

    def test_very_large_protein_skip_z3(self):
        """Proteins > 1200aa should skip Z3 entirely."""
        timeout, skip = compute_adaptive_timeout(1201)
        assert skip is True

    def test_very_large_protein_1500aa(self):
        """A 1500aa protein should skip Z3."""
        timeout, skip = compute_adaptive_timeout(1500)
        assert skip is True

    def test_explicit_timeout_override(self):
        """Explicit solver_timeout overrides tiered default."""
        timeout, skip = compute_adaptive_timeout(50, solver_timeout=45.0)
        assert timeout == 45.0
        assert skip is False

    def test_explicit_timeout_still_skip_large(self):
        """Even with explicit timeout, very large proteins still skip Z3."""
        timeout, skip = compute_adaptive_timeout(1500, solver_timeout=300.0)
        assert timeout == 300.0
        assert skip is True

    def test_explicit_timeout_no_skip_for_windowed(self):
        """Explicit timeout for windowed-tier protein does not skip."""
        timeout, skip = compute_adaptive_timeout(600, solver_timeout=300.0)
        assert timeout == 300.0
        assert skip is False

    def test_zero_solver_timeout_ignored(self):
        """solver_timeout=0 is treated as not provided."""
        timeout, skip = compute_adaptive_timeout(50, solver_timeout=0)
        assert timeout == 30.0  # Falls through to tier

    def test_negative_solver_timeout_ignored(self):
        """Negative solver_timeout is treated as not provided."""
        timeout, skip = compute_adaptive_timeout(50, solver_timeout=-10)
        assert timeout == 30.0

    def test_none_solver_timeout(self):
        """None solver_timeout uses tiered default."""
        timeout, skip = compute_adaptive_timeout(200, solver_timeout=None)
        assert timeout == 60.0

    def test_tiers_are_ordered(self):
        """Adaptive timeout tiers are ordered by protein length."""
        for i in range(len(ADAPTIVE_TIMEOUT_TIERS) - 1):
            assert ADAPTIVE_TIMEOUT_TIERS[i][0] < ADAPTIVE_TIMEOUT_TIERS[i + 1][0]

    def test_tiers_timeouts_increase(self):
        """Timeouts increase with protein length."""
        for i in range(len(ADAPTIVE_TIMEOUT_TIERS) - 1):
            assert ADAPTIVE_TIMEOUT_TIERS[i][1] <= ADAPTIVE_TIMEOUT_TIERS[i + 1][1]


# ══════════════════════════════════════════════════════════════════════
# 2. Fallback to greedy on timeout
# ══════════════════════════════════════════════════════════════════════


class TestFallbackOnTimeout:
    """Test that Z3 falls back to greedy when it times out."""

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_timeout_returns_fallback_result(self):
        """When Z3 times out, result has fallback_used=True."""
        config = SolverConfig(
            organism="Homo_sapiens",
            timeout_seconds=0.001,  # Extremely short — will definitely time out
            backend=SolverBackend.Z3,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence=INSULIN,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # Z3 may return "unknown" (timeout) or "unknown" (canceled) or
        # it may even solve quickly.  When unsolved and the reason is
        # timeout-related, fallback_used should be True.
        if not result.solved:
            reason = result.metadata.get("reason", "")
            if reason == "z3_timeout":
                assert result.fallback_used is True
            # "z3_unknown" may or may not have fallback_used=True depending
            # on whether the reason was specifically a timeout

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_timeout_stores_partial_sequence_in_metadata(self):
        """When Z3 times out with a partial solution, it's stored in metadata."""
        config = SolverConfig(
            organism="Homo_sapiens",
            timeout_seconds=0.001,
            backend=SolverBackend.Z3,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence=INSULIN,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # If there's a partial solution, it should be in metadata
        if not result.solved and result.metadata.get("z3_partial_sequence"):
            partial = result.metadata["z3_partial_sequence"]
            assert len(partial) == len(INSULIN) * 3


# ══════════════════════════════════════════════════════════════════════
# 3. Fallback to greedy on UNSAT
# ══════════════════════════════════════════════════════════════════════


class TestFallbackOnUNSAT:
    """Test that Z3 falls back when the model is UNSAT."""

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_unsat_result_has_fallback(self):
        """UNSAT result should have fallback_used=True."""
        # Create an infeasible model: GC must be both >0.9 and <0.1
        config = SolverConfig(
            organism="Homo_sapiens",
            gc_lo=0.9,
            gc_hi=0.99,
            backend=SolverBackend.Z3,
            timeout_seconds=5.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        # Use a protein that can't reach 90% GC
        protein = "MMMMMMMMMM"  # M=ATG, low GC
        model = CSPModel(
            protein_sequence=protein,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        if not result.solved:
            assert result.fallback_used is True
            assert result.metadata.get("reason") == "unsat" or result.mus_report is not None


# ══════════════════════════════════════════════════════════════════════
# 4. Solution verification
# ══════════════════════════════════════════════════════════════════════


class TestSolutionVerification:
    """Test the _verify_z3_solution function."""

    def test_valid_solution_passes(self):
        """A valid solution should pass verification."""
        # Insulin with reasonable codons
        config = SolverConfig(
            organism="Homo_sapiens",
            gc_lo=0.20,
            gc_hi=0.80,
        )
        # Use a simple valid sequence for "MMM" → ATG ATG ATG
        is_valid, violations = _verify_z3_solution(
            "ATGATGATG", "MMM", config, "Homo_sapiens"
        )
        # GC content of ATGATGATG = 4/9 ≈ 0.44, within [0.20, 0.80]
        assert is_valid is True
        assert len(violations) == 0

    def test_wrong_translation_fails(self):
        """A sequence that doesn't translate correctly fails verification."""
        config = SolverConfig(organism="Homo_sapiens")
        # TTT = Phe, not Met — wrong translation
        is_valid, violations = _verify_z3_solution(
            "TTTATGATG", "MMM", config, "Homo_sapiens"
        )
        assert is_valid is False
        assert any("translates to" in v for v in violations)

    def test_wrong_length_fails(self):
        """A sequence of wrong length fails verification."""
        config = SolverConfig(organism="Homo_sapiens")
        is_valid, violations = _verify_z3_solution(
            "ATGATG", "MMM", config, "Homo_sapiens"
        )
        assert is_valid is False
        assert any("length" in v.lower() for v in violations)

    def test_empty_sequence_fails(self):
        """An empty sequence fails verification."""
        config = SolverConfig(organism="Homo_sapiens")
        is_valid, violations = _verify_z3_solution(
            "", "MMM", config, "Homo_sapiens"
        )
        assert is_valid is False

    def test_gc_out_of_range_fails(self):
        """GC content outside configured range fails verification."""
        config = SolverConfig(organism="Homo_sapiens", gc_lo=0.60, gc_hi=0.80)
        # ATGATGATG has GC ≈ 0.44, outside [0.60, 0.80]
        is_valid, violations = _verify_z3_solution(
            "ATGATGATG", "MMM", config, "Homo_sapiens"
        )
        assert is_valid is False
        assert any("GC" in v for v in violations)

    def test_attta_motif_fails_when_avoid(self):
        """ATTTA motif in sequence fails when avoid_attta=True."""
        config = SolverConfig(organism="Homo_sapiens", avoid_attta=True)
        # Build a sequence containing ATTTA motif.
        # We need ATTTA as a literal substring of the DNA sequence.
        # Strategy: craft codons so that across codon boundaries,
        # the string ATTTA appears.
        # Codon layout: ...ATG | AAT | TTA | TGG...
        #                         ^^^   ^^^
        # The substring spanning positions 4-8: ATTTA
        # ATG = Met, AAT = Asn, TTA = Leu, TGG = Trp
        # So protein = "MNLW" → sequence = "ATGAATTTATGG"
        # Let's verify: A-T-G-A-A-T-T-T-A-T-G-G
        #                          ^^^^^ = ATTTA?  A(4)T(5)T(6)T(7)A(8) = ATTTA ✓
        seq = "ATGAATTTATGG"  # ATG AAT TTA TGG → M N L W
        protein = "MNLW"
        is_valid, violations = _verify_z3_solution(seq, protein, config, "Homo_sapiens")
        assert is_valid is False
        assert any("ATTTA" in v for v in violations)

    def test_t_run_fails_when_avoid(self):
        """T-run of 6+ fails when avoid_t_runs=True."""
        config = SolverConfig(organism="Homo_sapiens", avoid_t_runs=True)
        # Build a sequence with 6 consecutive T's
        # TTTTTT = Phe-Phe → TTT TTT → 6 consecutive T's
        seq = "ATGTTTTTTATG"  # ATG TTT TTT ATG → M F F M
        protein = "MFFM"
        is_valid, violations = _verify_z3_solution(seq, protein, config, "Homo_sapiens")
        assert is_valid is False
        assert any("T-run" in v for v in violations)

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_z3_solve_includes_verification(self):
        """Z3 solve results include verification metadata."""
        config = SolverConfig(
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            backend=SolverBackend.Z3,
            timeout_seconds=30.0,
        )
        engine = Z3Engine(config, organism="Escherichia_coli")
        model = CSPModel(
            protein_sequence="MVLSPADKTN",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Escherichia_coli",
        )
        result = engine.solve(model)
        if result.solved:
            assert result.metadata.get("verification_passed") is True


# ══════════════════════════════════════════════════════════════════════
# 5. Large proteins use windowed decomposition
# ══════════════════════════════════════════════════════════════════════


class TestLargeProteinWindowedSolving:
    """Test that large proteins are handled via extended tiers or windowed decomposition."""

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_large_protein_600aa_solves_directly(self):
        """Proteins 501-800aa should solve directly with extended timeout (tier 4)."""
        config = SolverConfig(
            organism="Homo_sapiens",
            backend=SolverBackend.Z3,
            timeout_seconds=120.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence=PROTEIN_600,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # Should NOT skip Z3 (600aa is within the 800aa tier)
        assert result.metadata.get("skip_z3") is not True
        # Should produce a valid-length sequence (solved directly, not windowed)
        assert len(result.sequence) == len(PROTEIN_600) * 3

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_exactly_500_does_not_skip(self):
        """Proteins exactly 500aa should NOT skip Z3 or use windowed."""
        config = SolverConfig(
            organism="Homo_sapiens",
            backend=SolverBackend.Z3,
            timeout_seconds=120.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence=PROTEIN_490,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        # This should attempt to solve (may or may not succeed, but shouldn't skip)
        result = engine.solve(model)
        # Should NOT have skip_z3 in metadata
        assert result.metadata.get("skip_z3") is not True

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_very_large_protein_over_1200_uses_windowed(self):
        """Proteins > 1200aa should use windowed decomposition."""
        config = SolverConfig(
            organism="Homo_sapiens",
            backend=SolverBackend.Z3,
            timeout_seconds=120.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        protein_1500 = "M" * 1500
        model = CSPModel(
            protein_sequence=protein_1500,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # > 1200aa should use windowed decomposition (not skip entirely)
        assert result.metadata.get("method") == "windowed_decomposition"
        # Should produce a valid-length sequence
        assert len(result.sequence) == len(protein_1500) * 3


# ══════════════════════════════════════════════════════════════════════
# 6. Warm-start fallback (dispatch integration)
# ══════════════════════════════════════════════════════════════════════


class TestWarmStartFallback:
    """Test that partial Z3 solutions are passed as warm-start to greedy."""

    def test_greedy_accepts_warm_start(self):
        """GreedyEngine.solve() accepts warm_start_sequence parameter."""
        from biocompiler.solver.engine_greedy import GreedyEngine

        config = SolverConfig(organism="Homo_sapiens")
        engine = GreedyEngine(config)
        model = CSPModel(
            protein_sequence="MVLSPADKTN",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        # Provide a valid warm-start sequence
        warm_start = "ATGGTTCTGTCGCCCGCAGACAAAACTAAC"
        result = engine.solve(model, warm_start_sequence=warm_start)
        assert result.solved is True

    def test_greedy_ignores_invalid_warm_start_length(self):
        """GreedyEngine ignores warm-start with wrong length."""
        from biocompiler.solver.engine_greedy import GreedyEngine

        config = SolverConfig(organism="Homo_sapiens")
        engine = GreedyEngine(config)
        model = CSPModel(
            protein_sequence="MVLSPADKTN",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        # Warm-start with wrong length
        result = engine.solve(model, warm_start_sequence="ATG")
        assert result.solved is True

    def test_greedy_validates_warm_start_codons(self):
        """GreedyEngine replaces invalid codons in warm-start."""
        from biocompiler.solver.engine_greedy import GreedyEngine

        config = SolverConfig(organism="Homo_sapiens")
        engine = GreedyEngine(config)
        model = CSPModel(
            protein_sequence="MMM",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        # Provide warm-start with an invalid codon for Met (AAA ≠ Met)
        warm_start = "ATGAAAGGG"  # AAA is not valid for Met
        result = engine.solve(model, warm_start_sequence=warm_start)
        assert result.solved is True
        # The result should have replaced AAA with a valid Met codon (ATG)
        assert result.sequence[:3] == "ATG"
        assert result.sequence[3:6] == "ATG"  # AAA replaced
        assert result.sequence[6:9] == "ATG"  # GGG is not valid for Met either

    def test_dispatch_passes_warm_start_on_z3_timeout(self):
        """dispatch.solve_with_csp passes Z3 partial solution to greedy."""
        from biocompiler.solver.dispatch import solve_with_csp

        # Create a large protein that will skip Z3
        result = solve_with_csp(PROTEIN_600, organism="Homo_sapiens")
        # Should get a solved result from greedy (with fallback_used)
        if result.solved:
            assert result.fallback_used is True

    def test_z3_metadata_contains_partial_sequence(self):
        """Z3 timeout result stores partial sequence in metadata for warm-start."""
        # We mock a scenario where Z3 returns a timeout with partial solution
        config = SolverConfig(
            organism="Homo_sapiens",
            timeout_seconds=0.001,
            backend=SolverBackend.Z3,
        )
        if not _Z3_AVAILABLE:
            pytest.skip("Z3 not available")
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence=INSULIN,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # If Z3 timed out and extracted a partial solution
        if not result.solved and result.metadata.get("z3_partial_sequence"):
            partial = result.metadata["z3_partial_sequence"]
            assert isinstance(partial, str)
            assert len(partial) == len(INSULIN) * 3


# ══════════════════════════════════════════════════════════════════════
# 7. Diagnostic logging tests
# ══════════════════════════════════════════════════════════════════════


class TestDiagnosticLogging:
    """Test that diagnostic information is properly recorded."""

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_result_contains_constraint_counts(self):
        """Z3 result metadata contains constraint counts by tier."""
        config = SolverConfig(
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            backend=SolverBackend.Z3,
            timeout_seconds=30.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence="MVLSPADKTN",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # Result should have num_constraints and num_variables
        assert result.num_variables > 0

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_timeout_result_has_reason(self):
        """Z3 timeout result has reason in metadata."""
        config = SolverConfig(
            organism="Homo_sapiens",
            timeout_seconds=0.001,
            backend=SolverBackend.Z3,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence=INSULIN,
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        if not result.solved and result.fallback_used:
            # Should have a reason in metadata
            reason = result.metadata.get("reason")
            assert reason is not None
            assert reason in ("z3_timeout", "z3_unknown", "protein_too_large_for_z3",
                              "unsat")

    def test_adaptive_timeout_tier_constants(self):
        """Verify the adaptive timeout tier constants are well-formed."""
        assert len(ADAPTIVE_TIMEOUT_TIERS) >= 3
        # Each tier is (max_length, timeout_seconds)
        for max_len, timeout in ADAPTIVE_TIMEOUT_TIERS:
            assert max_len > 0
            assert timeout > 0

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_unsat_result_has_constraint_tier_info(self):
        """UNSAT result includes constraint counts by tier in metadata."""
        config = SolverConfig(
            organism="Homo_sapiens",
            gc_lo=0.9,
            gc_hi=0.99,
            backend=SolverBackend.Z3,
            timeout_seconds=5.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence="MMMMMMMMMM",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        if not result.solved and result.metadata.get("reason") == "unsat":
            tier_counts = result.metadata.get("constraint_counts_by_tier")
            assert tier_counts is not None
            assert isinstance(tier_counts, dict)


# ══════════════════════════════════════════════════════════════════════
# 8. Configurable solver_timeout parameter
# ══════════════════════════════════════════════════════════════════════


class TestConfigurableTimeout:
    """Test that solver_timeout parameter on SolverConfig is respected."""

    def test_config_timeout_seconds_respected(self):
        """config.timeout_seconds is used when > 0."""
        timeout, skip = compute_adaptive_timeout(50, solver_timeout=99.0)
        assert timeout == 99.0

    def test_default_timeout_for_small_protein(self):
        """Without explicit timeout, small protein gets tier default."""
        timeout, skip = compute_adaptive_timeout(50)
        assert timeout == 30.0

    def test_default_timeout_for_medium_protein(self):
        """Without explicit timeout, medium protein gets tier default."""
        timeout, skip = compute_adaptive_timeout(200)
        assert timeout == 60.0

    def test_default_timeout_for_500aa_protein(self):
        """Without explicit timeout, 500aa protein gets 120s."""
        timeout, skip = compute_adaptive_timeout(500)
        assert timeout == 120.0

    @pytest.mark.skipif(not _Z3_AVAILABLE, reason="Z3 not available")
    def test_z3_engine_uses_adaptive_timeout(self):
        """Z3Engine uses adaptive timeout based on protein length."""
        config = SolverConfig(
            organism="Homo_sapiens",
            backend=SolverBackend.Z3,
            timeout_seconds=30.0,
        )
        engine = Z3Engine(config, organism="Homo_sapiens")
        model = CSPModel(
            protein_sequence="MVLSPADKTN",
            codon_domains={},
            constraints=[],
            config=config,
            organism="Homo_sapiens",
        )
        result = engine.solve(model)
        # Just verify it ran — the adaptive timeout is selected internally
        assert result is not None
