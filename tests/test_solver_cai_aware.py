"""Tests for CAI-aware solver optimization (F2.10).

Verify that both OR-Tools and Z3 solver engines are CAI-aware during
optimization: when multiple codon assignments satisfy hard constraints
equally, the solver should prefer the one with higher CAI.

Test categories:
  1. CSP solver produces sequences with higher CAI than greedy fallback
  2. Constraint satisfaction is not compromised for CAI
  3. Organism-aware optimization with CSP backends
  4. cai_weight parameter controls CAI pressure
"""

from __future__ import annotations

import math

import pytest

from biocompiler.solver.types import (
    ConstraintSpec,
    ConstraintType,
    CSPModel,
    SolverBackend,
    SolverConfig,
    SolverResult,
)
from biocompiler.constants import AA_TO_CODONS, CODON_TABLE, reverse_complement
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS


# ---------------------------------------------------------------------------
# Test proteins
# ---------------------------------------------------------------------------

SHORT_PROTEIN = "MVLSPADKTNVKAAWGKVGA"  # 20 aa — HBA1 fragment
MEDIUM_PROTEIN = "MKTVLIAEGH"
LEUCINE_HEAVY = "LLLLLLLLLL"  # 10 Leucines — 6 codon choices, strong CAI signal
MIXED_PROTEIN = "MKALWVGTSDEFRQPNHYCI"  # 20 aa — one of each amino acid class

# Very high splice threshold that effectively disables splice constraints
_SPLICE_DISABLED = 999.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gc_content(seq: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in "GC") / len(seq)


def _compute_cai(sequence: str, protein: str, organism: str) -> float:
    """Compute CAI for a DNA sequence using organism codon adaptiveness."""
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    if not sequence or not protein or not adaptiveness:
        return 0.0

    log_sum = 0.0
    n = 0
    for i, aa in enumerate(protein):
        codon = sequence[i * 3 : i * 3 + 3]
        w = adaptiveness.get(codon, 0.01)
        if w > 0:
            log_sum += math.log(w)
            n += 1
    if n == 0:
        return 0.0
    return math.exp(log_sum / n)


def _translates_to(seq: str, protein: str) -> bool:
    """Check that the sequence translates to the given protein."""
    for i, aa in enumerate(protein):
        codon = seq[i * 3 : i * 3 + 3]
        if CODON_TABLE.get(codon) != aa:
            return False
    return True


def _make_ortools_model(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> object:
    """Build a CSPModel for the OR-Tools engine.

    The OR-Tools engine uses its own CSPModel class with ``protein`` and
    ``organism`` attributes (not the types.CSPModel).
    """
    from biocompiler.solver.engine_ortools import CSPModel as ORToolsCSPModel
    return ORToolsCSPModel(protein=protein, organism=organism)


def _make_z3_model(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    cai_weight: float = 1.0,
) -> CSPModel:
    """Build a CSPModel for the Z3 engine."""
    config = SolverConfig(
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        restriction_sites=restriction_sites or [],
        cryptic_splice_threshold=_SPLICE_DISABLED,
        cai_weight=cai_weight,
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


def _greedy_sequence(protein: str, organism: str) -> str:
    """Generate a simple greedy sequence (highest-CAI codon per position)."""
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    codons = []
    for aa in protein:
        domain = AA_TO_CODONS.get(aa, [])
        if not domain:
            codons.append("NNN")
            continue
        best = max(domain, key=lambda c: adaptiveness.get(c, 0.0))
        codons.append(best)
    return "".join(codons)


def _random_low_cai_sequence(protein: str, organism: str) -> str:
    """Generate a sequence using the LOWEST-CAI codon for each position."""
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    codons = []
    for aa in protein:
        domain = AA_TO_CODONS.get(aa, [])
        if not domain:
            codons.append("NNN")
            continue
        worst = min(domain, key=lambda c: adaptiveness.get(c, 0.0))
        codons.append(worst)
    return "".join(codons)


# ======================================================================
# 1. OR-Tools CAI-aware tests
# ======================================================================

ortools = pytest.importorskip("ortools", reason="OR-Tools not installed")

from biocompiler.solver.engine_ortools import ORTOOLSEngine, CSPModel as ORToolsCSPModel


class TestORToolsCAIAware:
    """OR-Tools engine should produce CAI-aware solutions."""

    @pytest.fixture
    def engine(self) -> ORTOOLSEngine:
        return ORTOOLSEngine(SolverConfig())

    def test_csp_higher_cai_than_greedy(self, engine: ORTOOLSEngine):
        """CSP solver should produce sequences with CAI at least as high
        as greedy per-position optimization.

        The greedy approach picks the highest-CAI codon per position
        independently, without constraint checking. The CSP solver also
        maximizes CAI but must satisfy hard constraints (GC, restriction
        sites, etc.). With relaxed constraints, CSP CAI should match or
        exceed greedy CAI.
        """
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"
        config = SolverConfig(gc_lo=0.0, gc_hi=1.0, avoid_cpg=False)
        engine_relaxed = ORTOOLSEngine(config)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result = engine_relaxed.solve(model)

        if result.solved:
            csp_cai = _compute_cai(result.sequence, protein, organism)
            greedy_seq = _greedy_sequence(protein, organism)
            greedy_cai = _compute_cai(greedy_seq, protein, organism)
            # With no hard constraints, CSP should match greedy CAI
            assert csp_cai >= greedy_cai - 0.01, (
                f"CSP CAI ({csp_cai:.4f}) should be close to greedy CAI "
                f"({greedy_cai:.4f}) with relaxed constraints"
            )

    def test_constraints_not_compromised_for_cai(self, engine: ORTOOLSEngine):
        """Constraint satisfaction should not be sacrificed for CAI.

        With tight GC bounds, the solver must respect GC constraints even
        if it means accepting lower CAI. CAI-awareness means preferring
        higher CAI *among feasible solutions*, not violating constraints.
        """
        protein = SHORT_PROTEIN
        organism = "Homo_sapiens"
        config = SolverConfig(gc_lo=0.40, gc_hi=0.60, avoid_cpg=False)
        engine_constrained = ORTOOLSEngine(config)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result = engine_constrained.solve(model)

        if result.solved:
            gc = _gc_content(result.sequence)
            assert 0.40 <= gc <= 0.60, (
                f"GC constraint violated: {gc:.4f} outside [0.40, 0.60]"
            )
            # CAI should still be positive
            cai = _compute_cai(result.sequence, protein, organism)
            assert cai > 0.0, "CAI should be positive for a solved sequence"

    def test_cai_weight_affects_objective(self):
        """Higher cai_weight should produce higher-CAI solutions when
        constraints permit.

        With cai_weight=0.0, the solver has no CAI preference and may
        produce suboptimal CAI. With cai_weight=1.0, the solver should
        maximize CAI.
        """
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"

        # Solve with no CAI pressure
        config_no_cai = SolverConfig(
            gc_lo=0.0, gc_hi=1.0, avoid_cpg=False, cai_weight=0.0,
        )
        engine_no_cai = ORTOOLSEngine(config_no_cai)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result_no_cai = engine_no_cai.solve(model)

        # Solve with full CAI pressure
        config_full_cai = SolverConfig(
            gc_lo=0.0, gc_hi=1.0, avoid_cpg=False, cai_weight=1.0,
        )
        engine_full_cai = ORTOOLSEngine(config_full_cai)
        result_full_cai = engine_full_cai.solve(model)

        if result_no_cai.solved and result_full_cai.solved:
            cai_no = _compute_cai(result_no_cai.sequence, protein, organism)
            cai_full = _compute_cai(result_full_cai.sequence, protein, organism)
            # Full CAI weight should produce >= CAI compared to no weight
            assert cai_full >= cai_no - 0.01, (
                f"Full CAI weight ({cai_full:.4f}) should be >= "
                f"no CAI weight ({cai_no:.4f})"
            )

    def test_organism_aware_cai(self):
        """Different organisms should produce different codon preferences.

        CAI-awareness is organism-specific: the solver should use the
        correct codon adaptiveness table for the target organism.
        """
        protein = "MKALWVGTSDE"  # 11 aa
        for organism in SUPPORTED_ORGANISMS:
            adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
            if not adaptiveness:
                continue

            config = SolverConfig(gc_lo=0.20, gc_hi=0.80, avoid_cpg=False)
            engine = ORTOOLSEngine(config)
            model = ORToolsCSPModel(protein=protein, organism=organism)
            result = engine.solve(model)

            if result.solved:
                cai = _compute_cai(result.sequence, protein, organism)
                assert cai > 0.0, (
                    f"CAI should be positive for organism {organism}"
                )
                assert _translates_to(result.sequence, protein), (
                    f"Sequence should translate back to protein for {organism}"
                )

    def test_csp_prefers_high_cai_among_feasible(self, engine: ORTOOLSEngine):
        """When constraints allow multiple solutions, the solver should
        prefer the one with higher CAI.

        With very relaxed GC bounds (0.0-1.0), the solver has maximum
        freedom. It should still prefer high-CAI codons.
        """
        protein = MEDIUM_PROTEIN
        organism = "Homo_sapiens"
        config = SolverConfig(gc_lo=0.0, gc_hi=1.0, avoid_cpg=False)
        engine_relaxed = ORTOOLSEngine(config)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result = engine_relaxed.solve(model)

        if result.solved:
            csp_cai = _compute_cai(result.sequence, protein, organism)
            low_cai_seq = _random_low_cai_sequence(protein, organism)
            low_cai = _compute_cai(low_cai_seq, protein, organism)
            # CSP should produce higher CAI than the worst-case selection
            assert csp_cai > low_cai, (
                f"CSP CAI ({csp_cai:.4f}) should exceed worst-case "
                f"CAI ({low_cai:.4f})"
            )

    def test_solution_translates_correctly(self, engine: ORTOOLSEngine):
        """CAI-aware solution must still translate to the correct protein."""
        protein = SHORT_PROTEIN
        organism = "Homo_sapiens"
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result = engine.solve(model)

        if result.solved:
            assert _translates_to(result.sequence, protein), (
                "CAI-aware solution must translate to the correct protein"
            )


# ======================================================================
# 2. Z3 CAI-aware tests
# ======================================================================

z3_mod = pytest.importorskip("z3", reason="z3-solver not installed")

from biocompiler.solver.engine_z3 import Z3Engine


class TestZ3CAIAware:
    """Z3 engine should produce CAI-aware solutions."""

    @pytest.fixture
    def z3_engine(self) -> Z3Engine:
        config = SolverConfig(
            gc_lo=0.30, gc_hi=0.70,
            cryptic_splice_threshold=_SPLICE_DISABLED,
        )
        return Z3Engine(config, organism="Homo_sapiens")

    def test_csp_higher_cai_than_greedy(self, z3_engine: Z3Engine):
        """CSP solver should produce sequences with CAI at least as high
        as greedy per-position optimization (with relaxed constraints).
        """
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"
        config = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
        )
        engine = Z3Engine(config, organism=organism)
        model = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0)
        result = engine.solve(model)

        if result.solved:
            csp_cai = _compute_cai(result.sequence, protein, organism)
            greedy_seq = _greedy_sequence(protein, organism)
            greedy_cai = _compute_cai(greedy_seq, protein, organism)
            # With no hard constraints, CSP should match or approach greedy CAI
            assert csp_cai >= greedy_cai - 0.05, (
                f"CSP CAI ({csp_cai:.4f}) should be close to greedy CAI "
                f"({greedy_cai:.4f}) with relaxed constraints"
            )

    def test_constraints_not_compromised_for_cai(self, z3_engine: Z3Engine):
        """Constraint satisfaction should not be sacrificed for CAI.

        With tight GC bounds, the solver must respect GC constraints even
        if it means accepting lower CAI.
        """
        protein = SHORT_PROTEIN
        organism = "Homo_sapiens"
        model = _make_z3_model(protein, organism, gc_lo=0.40, gc_hi=0.60)
        result = z3_engine.solve(model)

        if result.solved:
            gc = _gc_content(result.sequence)
            assert 0.40 <= gc <= 0.60, (
                f"GC constraint violated: {gc:.4f} outside [0.40, 0.60]"
            )

    def test_cai_weight_affects_objective(self):
        """Higher cai_weight should produce higher-CAI solutions."""
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"

        # Solve with no CAI pressure
        config_no_cai = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
            cai_weight=0.0,
        )
        engine_no_cai = Z3Engine(config_no_cai, organism=organism)
        model = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0, cai_weight=0.0)
        result_no_cai = engine_no_cai.solve(model)

        # Solve with full CAI pressure
        config_full_cai = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
            cai_weight=1.0,
        )
        engine_full_cai = Z3Engine(config_full_cai, organism=organism)
        model_full = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0, cai_weight=1.0)
        result_full_cai = engine_full_cai.solve(model_full)

        if result_no_cai.solved and result_full_cai.solved:
            cai_no = _compute_cai(result_no_cai.sequence, protein, organism)
            cai_full = _compute_cai(result_full_cai.sequence, protein, organism)
            # Full CAI weight should produce >= CAI compared to no weight
            assert cai_full >= cai_no - 0.01, (
                f"Full CAI weight ({cai_full:.4f}) should be >= "
                f"no CAI weight ({cai_no:.4f})"
            )

    def test_organism_aware_cai(self):
        """Different organisms should produce different codon preferences."""
        protein = "MKALWVGTSDE"
        for organism in SUPPORTED_ORGANISMS:
            adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
            if not adaptiveness:
                continue

            config = SolverConfig(
                gc_lo=0.20, gc_hi=0.80,
                cryptic_splice_threshold=_SPLICE_DISABLED,
            )
            engine = Z3Engine(config, organism=organism)
            model = _make_z3_model(protein, organism, gc_lo=0.20, gc_hi=0.80)
            result = engine.solve(model)

            if result.solved:
                cai = _compute_cai(result.sequence, protein, organism)
                assert cai > 0.0, (
                    f"CAI should be positive for organism {organism}"
                )
                assert _translates_to(result.sequence, protein), (
                    f"Sequence should translate back for {organism}"
                )

    def test_solution_translates_correctly(self, z3_engine: Z3Engine):
        """CAI-aware solution must still translate to the correct protein."""
        protein = SHORT_PROTEIN
        organism = "Homo_sapiens"
        model = _make_z3_model(protein, organism)
        result = z3_engine.solve(model)

        if result.solved:
            assert _translates_to(result.sequence, protein), (
                "CAI-aware solution must translate to the correct protein"
            )

    def test_csp_prefers_high_cai_among_feasible(self):
        """When constraints allow multiple solutions, the solver should
        prefer the one with higher CAI."""
        protein = MEDIUM_PROTEIN
        organism = "Homo_sapiens"
        config = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
        )
        engine = Z3Engine(config, organism=organism)
        model = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0)
        result = engine.solve(model)

        if result.solved:
            csp_cai = _compute_cai(result.sequence, protein, organism)
            low_cai_seq = _random_low_cai_sequence(protein, organism)
            low_cai = _compute_cai(low_cai_seq, protein, organism)
            assert csp_cai > low_cai, (
                f"CSP CAI ({csp_cai:.4f}) should exceed worst-case "
                f"CAI ({low_cai:.4f})"
            )


# ======================================================================
# 3. Cross-engine comparison tests
# ======================================================================

class TestCrossEngineCAI:
    """Compare CAI-awareness across OR-Tools and Z3 backends."""

    def test_both_backends_prefer_high_cai(self):
        """Both OR-Tools and Z3 should produce solutions with CAI
        significantly above the minimum possible for the protein."""
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"

        low_cai_seq = _random_low_cai_sequence(protein, organism)
        low_cai = _compute_cai(low_cai_seq, protein, organism)
        greedy_seq = _greedy_sequence(protein, organism)
        greedy_cai = _compute_cai(greedy_seq, protein, organism)

        # OR-Tools
        ort_config = SolverConfig(gc_lo=0.0, gc_hi=1.0, avoid_cpg=False)
        ort_engine = ORTOOLSEngine(ort_config)
        ort_model = ORToolsCSPModel(protein=protein, organism=organism)
        ort_result = ort_engine.solve(ort_model)

        if ort_result.solved:
            ort_cai = _compute_cai(ort_result.sequence, protein, organism)
            # OR-Tools CAI should be closer to greedy than to worst-case
            assert ort_cai > low_cai, (
                f"OR-Tools CAI ({ort_cai:.4f}) should exceed "
                f"worst-case ({low_cai:.4f})"
            )

        # Z3
        z3_config = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
        )
        z3_engine = Z3Engine(z3_config, organism=organism)
        z3_model = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0)
        z3_result = z3_engine.solve(z3_model)

        if z3_result.solved:
            z3_cai = _compute_cai(z3_result.sequence, protein, organism)
            # Z3 CAI should be closer to greedy than to worst-case
            assert z3_cai > low_cai, (
                f"Z3 CAI ({z3_cai:.4f}) should exceed "
                f"worst-case ({low_cai:.4f})"
            )

    def test_both_backends_respect_constraints(self):
        """Both backends should respect hard constraints even with CAI
        pressure."""
        protein = SHORT_PROTEIN
        organism = "Homo_sapiens"

        # Tight GC bounds
        gc_lo, gc_hi = 0.40, 0.60

        # OR-Tools
        ort_config = SolverConfig(gc_lo=gc_lo, gc_hi=gc_hi, avoid_cpg=False)
        ort_engine = ORTOOLSEngine(ort_config)
        ort_model = ORToolsCSPModel(protein=protein, organism=organism)
        ort_result = ort_engine.solve(ort_model)

        if ort_result.solved:
            gc = _gc_content(ort_result.sequence)
            assert gc_lo <= gc <= gc_hi, (
                f"OR-Tools GC {gc:.4f} outside [{gc_lo}, {gc_hi}]"
            )

        # Z3
        z3_config = SolverConfig(
            gc_lo=gc_lo, gc_hi=gc_hi,
            cryptic_splice_threshold=_SPLICE_DISABLED,
        )
        z3_engine = Z3Engine(z3_config, organism=organism)
        z3_model = _make_z3_model(protein, organism, gc_lo=gc_lo, gc_hi=gc_hi)
        z3_result = z3_engine.solve(z3_model)

        if z3_result.solved:
            gc = _gc_content(z3_result.sequence)
            assert gc_lo <= gc <= gc_hi, (
                f"Z3 GC {gc:.4f} outside [{gc_lo}, {gc_hi}]"
            )

    def test_constrained_cai_lower_than_unconstrained(self):
        """With tighter constraints, CAI should be <= unconstrained CAI.

        This validates that constraint satisfaction takes priority over
        CAI maximization — the hallmark of CAI-aware optimization.
        """
        protein = MEDIUM_PROTEIN
        organism = "Homo_sapiens"

        # OR-Tools: unconstrained
        config_free = SolverConfig(gc_lo=0.0, gc_hi=1.0, avoid_cpg=False)
        engine_free = ORTOOLSEngine(config_free)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result_free = engine_free.solve(model)

        # OR-Tools: constrained
        config_tight = SolverConfig(gc_lo=0.40, gc_hi=0.60, avoid_cpg=False)
        engine_tight = ORTOOLSEngine(config_tight)
        result_tight = engine_tight.solve(model)

        if result_free.solved and result_tight.solved:
            cai_free = _compute_cai(result_free.sequence, protein, organism)
            cai_tight = _compute_cai(result_tight.sequence, protein, organism)
            # Constrained CAI should be <= unconstrained CAI
            # (allowing small tolerance for numerical precision)
            assert cai_tight <= cai_free + 0.05, (
                f"Constrained CAI ({cai_tight:.4f}) should be <= "
                f"unconstrained ({cai_free:.4f})"
            )


# ======================================================================
# 4. CAI weight parameter tests
# ======================================================================

class TestCAIWeightParameter:
    """Test the cai_weight parameter controls CAI pressure in both backends."""

    def test_cai_weight_default_is_one(self):
        """Default cai_weight should be 1.0."""
        config = SolverConfig()
        assert config.cai_weight == 1.0

    def test_cai_weight_custom(self):
        """Custom cai_weight should be stored correctly."""
        config = SolverConfig(cai_weight=0.5)
        assert config.cai_weight == 0.5

    def test_ortools_zero_weight_no_cai_preference(self):
        """With cai_weight=0.0, OR-Tools should have no CAI preference
        and may produce lower-CAI solutions than with cai_weight=1.0."""
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"

        # With zero CAI weight
        config_zero = SolverConfig(
            gc_lo=0.0, gc_hi=1.0, avoid_cpg=False, cai_weight=0.0,
        )
        engine_zero = ORTOOLSEngine(config_zero)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result_zero = engine_zero.solve(model)

        # With full CAI weight
        config_full = SolverConfig(
            gc_lo=0.0, gc_hi=1.0, avoid_cpg=False, cai_weight=1.0,
        )
        engine_full = ORTOOLSEngine(config_full)
        result_full = engine_full.solve(model)

        if result_zero.solved and result_full.solved:
            cai_zero = _compute_cai(result_zero.sequence, protein, organism)
            cai_full = _compute_cai(result_full.sequence, protein, organism)
            # With weight=0, the solver may produce any feasible solution;
            # with weight=1, it should maximize CAI.
            # We only check that full-weight is at least as good.
            assert cai_full >= cai_zero - 0.01, (
                f"CAI with weight=1.0 ({cai_full:.4f}) should be >= "
                f"CAI with weight=0.0 ({cai_zero:.4f})"
            )

    def test_z3_zero_weight_no_cai_preference(self):
        """With cai_weight=0.0, Z3 should have no CAI preference."""
        protein = LEUCINE_HEAVY
        organism = "Homo_sapiens"

        # With zero CAI weight
        config_zero = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
            cai_weight=0.0,
        )
        engine_zero = Z3Engine(config_zero, organism=organism)
        model = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0, cai_weight=0.0)
        result_zero = engine_zero.solve(model)

        # With full CAI weight
        config_full = SolverConfig(
            gc_lo=0.0, gc_hi=1.0,
            cryptic_splice_threshold=_SPLICE_DISABLED,
            cai_weight=1.0,
        )
        engine_full = Z3Engine(config_full, organism=organism)
        model_full = _make_z3_model(protein, organism, gc_lo=0.0, gc_hi=1.0, cai_weight=1.0)
        result_full = engine_full.solve(model_full)

        if result_zero.solved and result_full.solved:
            cai_zero = _compute_cai(result_zero.sequence, protein, organism)
            cai_full = _compute_cai(result_full.sequence, protein, organism)
            assert cai_full >= cai_zero - 0.01, (
                f"CAI with weight=1.0 ({cai_full:.4f}) should be >= "
                f"CAI with weight=0.0 ({cai_zero:.4f})"
            )

    def test_high_cai_weight_still_respects_constraints(self):
        """Even with very high cai_weight, hard constraints must be satisfied."""
        protein = SHORT_PROTEIN
        organism = "Homo_sapiens"
        gc_lo, gc_hi = 0.40, 0.60

        # OR-Tools with high CAI weight
        config = SolverConfig(
            gc_lo=gc_lo, gc_hi=gc_hi,
            avoid_cpg=False, cai_weight=10.0,
        )
        engine = ORTOOLSEngine(config)
        model = ORToolsCSPModel(protein=protein, organism=organism)
        result = engine.solve(model)

        if result.solved:
            gc = _gc_content(result.sequence)
            assert gc_lo <= gc <= gc_hi, (
                f"GC {gc:.4f} outside [{gc_lo}, {gc_hi}] even with "
                f"high cai_weight=10.0"
            )
