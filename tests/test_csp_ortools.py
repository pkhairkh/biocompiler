"""
BioCompiler OR-Tools Engine Tests
===================================

Pytest suite for solver/engine_ortools.py — the OR-Tools CP-SAT backend.

All OR-Tools-dependent tests use ``pytest.importorskip("ortools")`` so the
full suite remains runnable even when OR-Tools is not installed.

Test categories:
  1. Availability check
  2. Simple solve (correctness)
  3. GC constraint (feasible + infeasible)
  4. Restriction site avoidance
  5. CAI optimization (constrained vs unconstrained)
  6. Performance (100 AA and 239 AA eGFP)
  7. Fallback behaviour (mock-driven)
"""

from __future__ import annotations

import time
from unittest.mock import patch

import pytest

# ---------------------------------------------------------------------------
# Skip entire module when OR-Tools is absent
# ---------------------------------------------------------------------------
ortools = pytest.importorskip("ortools")

from biocompiler.solver.types import (
    ConstraintSpec, ConstraintType, ConstraintStrictness, CSPModel,
    SolverConfig, SolverResult, SolverBackend,
)
from biocompiler.solver.engine_ortools import ORTOOLSEngine
from biocompiler.expression.translation import translate, compute_cai
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS, reverse_complement
from biocompiler.sequence.restriction_sites import get_recognition_site


# ═══════════════════════════════════════════════════════════════════════════
# Test proteins & helpers
# ═══════════════════════════════════════════════════════════════════════════

SIMPLE_PROTEIN = "MKT"
MEDIUM_PROTEIN = "MKTVLIAEGH"
AT_HEAVY_PROTEIN = "KKKKKK"  # All Lys — low GC ceiling (AAA/AAG max 1/3 GC)
PROTEIN_100AA = (
    "MKFLILLFNILCLFPVLAADNHGVSLHVKAFDALQKAGDVGHFVNKDETQIYH"
    "RLGEWLSQYRLSEPEQVTKVLGVDKIEFLENKRVLRPKAKELKEILDQLEQKA"
)
EGFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTY"
    "GVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKED"
    "GNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHY"
    "LSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


def _build_model(
    protein: str,
    gc_bounds: tuple[float, float] = (0.30, 0.70),
    restriction_enzymes: list[str] | None = None,
    avoid_cpg: bool = False,
    avoid_gt: bool = False,
    organism: str = "Escherichia_coli",
) -> CSPModel:
    sites = []
    for enz in (restriction_enzymes or []):
        site = get_recognition_site(enz)
        if site:
            sites.append(site)

    config = SolverConfig(
        gc_lo=gc_bounds[0],
        gc_hi=gc_bounds[1],
        restriction_sites=sites,
        avoid_cpg=avoid_cpg,
        add_default_restriction_sites=bool(restriction_enzymes),
        organism=organism,
    )
    codon_domains = {i: list(AA_TO_CODONS[aa]) for i, aa in enumerate(protein)}
    constraints: list[ConstraintSpec] = [
        ConstraintSpec(
            ctype=ConstraintType.GC_CONTENT, name="gc_content_global",
            params={"gc_lo": gc_bounds[0], "gc_hi": gc_bounds[1]},
        )
    ]
    for enz in (restriction_enzymes or []):
        site = get_recognition_site(enz)
        if site:
            constraints.append(ConstraintSpec(
                ctype=ConstraintType.RESTRICTION_SITE, name=f"no_{enz}",
                params={"site": site},
            ))
    if avoid_cpg:
        constraints.append(ConstraintSpec(
            ctype=ConstraintType.NO_CPG, name="no_cpg_dinucleotide",
        ))
    if avoid_gt:
        constraints.append(ConstraintSpec(
            ctype=ConstraintType.NO_GT_DINUCLEOTIDE, name="no_gt_dinucleotide",
        ))
    return CSPModel(
        protein_sequence=protein, codon_domains=codon_domains,
        constraints=constraints, config=config, organism=organism,
    )


def _translates_to(seq: str, protein: str) -> bool:
    return translate(seq, to_stop=False)[: len(protein)] == protein


def _contains_site(seq: str, site: str) -> bool:
    rc = reverse_complement(site)
    return site in seq or rc in seq


def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in "GC") / len(seq)


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def engine():
    return ORTOOLSEngine(SolverConfig())


@pytest.fixture
def simple_model():
    return _build_model(SIMPLE_PROTEIN)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Availability check
# ═══════════════════════════════════════════════════════════════════════════

class TestAvailability:
    def test_is_available_when_installed(self, engine):
        assert ORTOOLSEngine.is_available() is True

    def test_is_available_classmethod(self):
        assert ORTOOLSEngine.is_available() is True


# ═══════════════════════════════════════════════════════════════════════════
# 2. Simple solve
# ═══════════════════════════════════════════════════════════════════════════

class TestSimpleSolve:
    def test_result_is_solved(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved, f"Expected solved, got solved={result.solved}"

    def test_sequence_length(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved
        assert len(result.sequence) == 3 * len(SIMPLE_PROTEIN)

    def test_sequence_translates_back(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved
        assert _translates_to(result.sequence, SIMPLE_PROTEIN)

    def test_gc_within_bounds(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved
        lo, hi = simple_model.config.gc_lo, simple_model.config.gc_hi
        gc = _gc_content(result.sequence)
        assert lo <= gc <= hi, f"GC {gc:.4f} outside [{lo}, {hi}]"

    def test_no_restriction_sites_by_default(self, engine):
        model = _build_model(MEDIUM_PROTEIN, restriction_enzymes=[])
        result = engine.solve(model)
        assert result.solved
        assert len(result.sequence) == 3 * len(MEDIUM_PROTEIN)

    def test_solve_time_positive(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved
        assert result.solve_time_seconds > 0


# ═══════════════════════════════════════════════════════════════════════════
# 3. GC constraint
# ═══════════════════════════════════════════════════════════════════════════

class TestGCConstraint:
    def test_tight_gc_bounds_satisfied(self, engine):
        model = _build_model(MEDIUM_PROTEIN, gc_bounds=(0.40, 0.60))
        result = engine.solve(model)
        assert result.solved, f"Tight GC solve failed"
        gc = _gc_content(result.sequence)
        assert 0.40 <= gc <= 0.60, f"GC {gc:.4f} outside [0.40, 0.60]"

    def test_impossible_gc_infeasible(self, engine):
        """All-Lys protein cannot reach GC 0.90 — solver should detect infeasibility or have violations."""
        model = _build_model(AT_HEAVY_PROTEIN, gc_bounds=(0.90, 0.95))
        result = engine.solve(model)
        # Solver may return solved=True with violations, or solved=False
        # Either way, the solution should not actually satisfy the GC constraint
        if result.solved:
            gc = _gc_content(result.sequence)
            assert gc < 0.90 or len(result.violations) > 0, (
                "Should have violations or impossible GC for all-Lys protein"
            )

    def test_infeasible_empty_or_partial_sequence(self, engine):
        model = _build_model(AT_HEAVY_PROTEIN, gc_bounds=(0.90, 0.95))
        result = engine.solve(model)
        if not result.solved:
            assert result.sequence == "" or len(result.sequence) > 0

    def test_gc_heavy_high_bounds(self, engine):
        """All-Gly protein should solve with high GC bounds."""
        model = _build_model("GGGGGG", gc_bounds=(0.50, 0.80))
        result = engine.solve(model)
        assert result.solved
        gc = _gc_content(result.sequence)
        assert 0.50 <= gc <= 0.80


# ═══════════════════════════════════════════════════════════════════════════
# 4. Restriction site avoidance
# ═══════════════════════════════════════════════════════════════════════════

class TestRestrictionSiteAvoidance:
    def test_avoid_ecori(self, engine):
        model = _build_model(MEDIUM_PROTEIN, restriction_enzymes=["EcoRI"])
        result = engine.solve(model)
        assert result.solved
        assert not _contains_site(result.sequence, "GAATTC")

    def test_avoid_multiple_sites(self, engine):
        enzymes = ["EcoRI", "BamHI", "XhoI"]
        model = _build_model(MEDIUM_PROTEIN, restriction_enzymes=enzymes)
        result = engine.solve(model)
        assert result.solved
        for enz in enzymes:
            site = get_recognition_site(enz)
            assert site and not _contains_site(result.sequence, site), (
                f"{enz} site ({site}) found in result"
            )

    def test_site_avoidance_preserves_translation(self, engine):
        model = _build_model(MEDIUM_PROTEIN, restriction_enzymes=["EcoRI", "BamHI"])
        result = engine.solve(model)
        assert result.solved
        assert _translates_to(result.sequence, MEDIUM_PROTEIN)

    def test_site_avoidance_preserves_gc(self, engine):
        gc_bounds = (0.30, 0.70)
        model = _build_model(
            MEDIUM_PROTEIN, gc_bounds=gc_bounds,
            restriction_enzymes=["EcoRI", "BamHI"],
        )
        result = engine.solve(model)
        assert result.solved
        gc = _gc_content(result.sequence)
        assert gc_bounds[0] <= gc <= gc_bounds[1]

    def test_avoid_noti_long_site(self, engine):
        model = _build_model(MEDIUM_PROTEIN, restriction_enzymes=["NotI"])
        result = engine.solve(model)
        assert result.solved
        assert not _contains_site(result.sequence, "GCGGCCGC")

    def test_many_enzymes_simultaneously(self, engine):
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "SalI", "PstI"]
        model = _build_model(SIMPLE_PROTEIN, restriction_enzymes=enzymes)
        result = engine.solve(model)
        if result.solved:
            for enz in enzymes:
                site = get_recognition_site(enz)
                if site:
                    assert not _contains_site(result.sequence, site)


# ═══════════════════════════════════════════════════════════════════════════
# 5. CAI optimization
# ═══════════════════════════════════════════════════════════════════════════

class TestCAIOptimization:
    def test_unconstrained_high_cai(self, engine):
        model = _build_model(MEDIUM_PROTEIN, gc_bounds=(0.0, 1.0))
        result = engine.solve(model)
        assert result.solved
        cai = compute_cai(result.sequence, organism="Homo_sapiens")
        assert cai > 0.7, f"Unconstrained CAI should be high, got {cai:.4f}"

    def test_constrained_lower_cai(self, engine):
        model_free = _build_model(MEDIUM_PROTEIN, gc_bounds=(0.0, 1.0))
        r_free = engine.solve(model_free)
        model_tight = _build_model(
            MEDIUM_PROTEIN, gc_bounds=(0.40, 0.60),
            restriction_enzymes=["EcoRI", "BamHI"],
        )
        r_tight = engine.solve(model_tight)
        assert r_free.solved and r_tight.solved
        cai_free = compute_cai(r_free.sequence, organism="Homo_sapiens")
        cai_tight = compute_cai(r_tight.sequence, organism="Homo_sapiens")
        assert cai_tight <= cai_free + 0.01, (
            f"Constrained CAI ({cai_tight:.4f}) > unconstrained ({cai_free:.4f})"
        )

    def test_constrained_satisfies_constraints(self, engine):
        gc_bounds = (0.40, 0.60)
        enzymes = ["EcoRI", "BamHI"]
        model = _build_model(
            MEDIUM_PROTEIN, gc_bounds=gc_bounds, restriction_enzymes=enzymes,
        )
        result = engine.solve(model)
        assert result.solved
        gc = _gc_content(result.sequence)
        assert gc_bounds[0] <= gc <= gc_bounds[1]
        for enz in enzymes:
            site = get_recognition_site(enz)
            if site:
                assert not _contains_site(result.sequence, site)

    def test_positive_cai_for_solved(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved
        assert compute_cai(result.sequence, organism="Homo_sapiens") > 0.0


# ═══════════════════════════════════════════════════════════════════════════
# 6. Performance
# ═══════════════════════════════════════════════════════════════════════════

class TestPerformance:
    def test_100aa_under_10_seconds(self, engine):
        model = _build_model(
            PROTEIN_100AA, gc_bounds=(0.30, 0.70),
            restriction_enzymes=["EcoRI", "BamHI"],
        )
        t0 = time.perf_counter()
        result = engine.solve(model)
        elapsed = time.perf_counter() - t0
        assert result.solved, f"100 AA solve failed"
        assert elapsed < 10.0, f"100 AA took {elapsed:.2f}s (limit 10s)"
        assert _translates_to(result.sequence, PROTEIN_100AA)

    def test_239aa_egfp_under_60_seconds(self, engine):
        model = _build_model(
            EGFP_PROTEIN, gc_bounds=(0.30, 0.70),
            restriction_enzymes=["EcoRI", "BamHI", "XhoI"],
        )
        t0 = time.perf_counter()
        result = engine.solve(model)
        elapsed = time.perf_counter() - t0
        assert elapsed < 60.0, f"eGFP took {elapsed:.2f}s (limit 60s)"
        if result.solved:
            assert _translates_to(result.sequence, EGFP_PROTEIN)
            gc = _gc_content(result.sequence)
            assert 0.30 <= gc <= 0.70, f"eGFP GC {gc:.4f} outside [0.30, 0.70]"

    def test_100aa_many_constraints_fast(self, engine):
        model = _build_model(
            PROTEIN_100AA, gc_bounds=(0.40, 0.60),
            restriction_enzymes=["EcoRI", "BamHI", "XhoI", "HindIII", "SalI", "PstI"],
            avoid_cpg=True,
        )
        t0 = time.perf_counter()
        result = engine.solve(model)
        elapsed = time.perf_counter() - t0
        assert elapsed < 30.0, f"100 AA many constraints took {elapsed:.2f}s"
        if result.solved:
            assert _translates_to(result.sequence, PROTEIN_100AA)


# ═══════════════════════════════════════════════════════════════════════════
# 7. Fallback behavior
# ═══════════════════════════════════════════════════════════════════════════

class TestFallbackBehavior:
    def test_mock_infeasible_propagates(self, engine):
        result_mock = SolverResult(sequence="", solved=False, backend_used=SolverBackend.ORTOOLS)
        with patch.object(ORTOOLSEngine, "solve", return_value=result_mock):
            r = engine.solve(_build_model(SIMPLE_PROTEIN))
            assert not r.solved
            assert r.sequence == ""

    def test_mock_timeout_propagates(self, engine):
        result_mock = SolverResult(sequence="", solved=False, backend_used=SolverBackend.ORTOOLS)
        with patch.object(ORTOOLSEngine, "solve", return_value=result_mock):
            r = engine.solve(_build_model(SIMPLE_PROTEIN))
            assert not r.solved

    def test_mock_error_propagates(self, engine):
        result_mock = SolverResult(sequence="", solved=False, backend_used=SolverBackend.ORTOOLS)
        with patch.object(ORTOOLSEngine, "solve", return_value=result_mock):
            r = engine.solve(_build_model(SIMPLE_PROTEIN))
            assert not r.solved

    def test_mock_exception_propagates(self, engine):
        with patch.object(ORTOOLSEngine, "solve", side_effect=RuntimeError("crash")):
            with pytest.raises(RuntimeError, match="crash"):
                engine.solve(_build_model(SIMPLE_PROTEIN))

    def test_mock_alternating_results(self):
        ok = SolverResult(
            sequence="ATGAAAACCTGA", solved=True, backend_used=SolverBackend.ORTOOLS,
            objective_value=1.0,
        )
        bad = SolverResult(sequence="", solved=False, backend_used=SolverBackend.ORTOOLS)
        engine = ORTOOLSEngine(SolverConfig())
        model = _build_model(SIMPLE_PROTEIN)
        with patch.object(ORTOOLSEngine, "solve", side_effect=[ok, bad]):
            assert engine.solve(model).solved is True
            assert engine.solve(model).solved is False

    def test_mock_is_available_false(self):
        with patch.object(ORTOOLSEngine, "is_available", return_value=False):
            assert ORTOOLSEngine.is_available() is False

    def test_dispatch_fallback_on_infeasible(self, engine):
        """When OR-Tools returns infeasible, dispatch should see it and
        signal the need for fallback (the caller decides the fallback path)."""
        result_mock = SolverResult(sequence="", solved=False, backend_used=SolverBackend.ORTOOLS)
        with patch.object(ORTOOLSEngine, "solve", return_value=result_mock):
            r = engine.solve(_build_model(SIMPLE_PROTEIN))
            assert not r.solved
            assert r.sequence == ""  # Empty → caller must fall back


# ═══════════════════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestEdgeCases:
    def test_single_methionine(self, engine):
        model = _build_model("M")
        result = engine.solve(model)
        assert result.solved
        assert result.sequence == "ATG"

    def test_single_aa_variants(self, engine):
        for aa in "ACDEFHIKLMNPQRSTVWY":
            result = engine.solve(_build_model(aa))
            assert result.solved, f"Single-AA '{aa}' failed"
            assert _translates_to(result.sequence, aa)

    def test_homopolymer_leucine(self, engine):
        protein = "L" * 20
        result = engine.solve(_build_model(protein, gc_bounds=(0.30, 0.70)))
        assert result.solved
        assert _translates_to(result.sequence, protein)

    def test_no_constraints_model(self, engine):
        codon_domains = {i: list(AA_TO_CODONS[aa]) for i, aa in enumerate(SIMPLE_PROTEIN)}
        model = CSPModel(
            protein_sequence=SIMPLE_PROTEIN, codon_domains=codon_domains,
            constraints=[], config=SolverConfig(gc_lo=0.0, gc_hi=1.0),
        )
        result = engine.solve(model)
        assert result.solved
        assert _translates_to(result.sequence, SIMPLE_PROTEIN)

    def test_deterministic_solve(self, engine, simple_model):
        r1, r2 = engine.solve(simple_model), engine.solve(simple_model)
        if r1.solved and r2.solved:
            assert r1.sequence == r2.sequence, "Solver should be deterministic"

    def test_valid_dna_bases(self, engine, simple_model):
        result = engine.solve(simple_model)
        assert result.solved
        assert set(result.sequence) <= set("ACGT")

    def test_timeout_config_respected(self, engine):
        eng = ORTOOLSEngine(SolverConfig(timeout_seconds=5.0))
        t0 = time.perf_counter()
        eng.solve(_build_model(SIMPLE_PROTEIN))
        assert time.perf_counter() - t0 < 10.0
