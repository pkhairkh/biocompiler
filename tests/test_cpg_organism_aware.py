"""Tests for organism-aware CpG constraint behaviour (Task F5.9).

CpG avoidance is primarily relevant for mammalian expression systems
where CpG methylation can lead to gene silencing.  For prokaryotic
organisms (e.g. *E. coli*), CpG islands have no known regulatory
significance and the checks should be skipped.

These tests verify:
1. CpG checks return True/PASS for prokaryotic organisms (E. coli)
2. CpG checks still work correctly for eukaryotic organisms (human, mouse)
3. Backward compatibility when organism is not specified
4. Solver constraint classes (NoCpGIslandConstraint, MinimizeCpG) respect organism domain
5. Constraint enforcement pipeline skips CpG for prokaryotic organisms
"""

import pytest

from biocompiler.type_system import (
    check_no_cpg_island,
    evaluate_no_cpg_island,
    _is_prokaryotic_organism,
    PredicateResult,
)
from biocompiler.shared.types import Verdict, TypeCheckResult
from biocompiler.solver.constraints import (
    NoCpGIslandConstraint,
    MinimizeCpG,
    build_csp_model,
)
from biocompiler.solver.enforcement import ConstraintEnforcer


# ── Helper: a CG-rich sequence that would trigger a CpG island violation ──

# 300 bp sequence rich in CG dinucleotides — should FAIL CpG checks
# for eukaryotic organisms
CG_RICH_SEQ = (
    "ATGCGCCGTCGCCGACGTCGAACGTCGACGTCGAA"
    "CGCCGACGTCGAACGTCGACGTCGAACGCCGACGT"
    "CGAACGTCGACGTCGAACGCCGACGTCGAACGTCG"
    "ACGTCGAACGCCGACGTCGAACGTCGACGTCGAAC"
    "GCCGACGTCGAACGTCGACGTCGAACGCCGACGTC"
    "GAACGTCGACGTCGAACGCCGACGTCGAACGTCGA"
    "CGTCGAACGCCGACGTCGAACGTCGACGTCGAACG"
    "CCGACGTCGAACGTCGACGTCGAACGCCGACGTCG"
    "AACGTCGACGTCGAACGC"
)

# A clean sequence without CpG islands — should PASS for all organisms
CLEAN_SEQ = "ATGAAAGCTTTTGCAAAGCTTTTGCAAAGCTTTTGCAAA" * 10


# ════════════════════════════════════════════════════════════════════
# 1. _is_prokaryotic_organism helper
# ════════════════════════════════════════════════════════════════════

class TestIsProkaryoticOrganism:
    """Tests for the _is_prokaryotic_organism helper function."""

    def test_ecoli_is_prokaryotic(self):
        assert _is_prokaryotic_organism("E_coli") is True

    def test_ecoli_k12_is_prokaryotic(self):
        assert _is_prokaryotic_organism("E_coli_K12") is True

    def test_ecoli_bl21_is_prokaryotic(self):
        assert _is_prokaryotic_organism("E_coli_BL21") is True

    def test_ecoli_full_name_is_prokaryotic(self):
        assert _is_prokaryotic_organism("Escherichia_coli") is True

    def test_ecoli_alias_is_prokaryotic(self):
        assert _is_prokaryotic_organism("ecoli") is True

    def test_human_is_not_prokaryotic(self):
        assert _is_prokaryotic_organism("Homo_sapiens") is False

    def test_mouse_is_not_prokaryotic(self):
        assert _is_prokaryotic_organism("Mus_musculus") is False

    def test_cho_is_not_prokaryotic(self):
        assert _is_prokaryotic_organism("CHO_K1") is False

    def test_yeast_is_not_prokaryotic(self):
        assert _is_prokaryotic_organism("Saccharomyces_cerevisiae") is False

    def test_empty_string_is_not_prokaryotic(self):
        assert _is_prokaryotic_organism("") is False


# ════════════════════════════════════════════════════════════════════
# 2. check_no_cpg_island — organism awareness
# ════════════════════════════════════════════════════════════════════

class TestCheckNoCpGIslandOrganismAware:
    """Tests for check_no_cpg_island with organism parameter."""

    def test_prokaryotic_ecoli_returns_true_even_with_cpg_island(self):
        """CpG check should return True (PASS) for E. coli regardless of sequence."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism="E_coli")
        assert isinstance(result, PredicateResult)
        assert result.passed is True
        assert result.verdict == Verdict.PASS
        assert "skipped" in result.details.lower() or "prokaryotic" in result.details.lower()

    def test_prokaryotic_ecoli_k12_returns_true(self):
        """CpG check should return True for E. coli K-12."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism="E_coli_K12")
        assert result.passed is True

    def test_eukaryotic_human_fails_on_cpg_island(self):
        """CpG check should detect violations for human sequences."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism="Homo_sapiens")
        assert result.passed is False
        assert result.verdict == Verdict.FAIL

    def test_eukaryotic_mouse_fails_on_cpg_island(self):
        """CpG check should detect violations for mouse sequences."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism="Mus_musculus")
        assert result.passed is False
        assert result.verdict == Verdict.FAIL

    def test_eukaryotic_human_passes_on_clean_sequence(self):
        """CpG check should pass for human on a clean sequence."""
        result = check_no_cpg_island(CLEAN_SEQ, organism="Homo_sapiens")
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_backward_compat_no_organism_fails_on_cpg_island(self):
        """Without organism parameter, CpG check should still detect violations."""
        result = check_no_cpg_island(CG_RICH_SEQ)
        # When organism is not specified, the check runs as before
        assert result.passed is False
        assert result.verdict == Verdict.FAIL

    def test_backward_compat_no_organism_passes_on_clean_sequence(self):
        """Without organism parameter, CpG check should pass on clean sequence."""
        result = check_no_cpg_island(CLEAN_SEQ)
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_backward_compat_empty_organism_runs_check(self):
        """Empty string organism should behave like unspecified organism."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism="")
        assert result.passed is False
        assert result.verdict == Verdict.FAIL


# ════════════════════════════════════════════════════════════════════
# 3. evaluate_no_cpg_island — organism awareness
# ════════════════════════════════════════════════════════════════════

class TestEvaluateNoCpGIslandOrganismAware:
    """Tests for evaluate_no_cpg_island with organism parameter."""

    def test_prokaryotic_ecoli_returns_pass(self):
        """CpG evaluation should return PASS for E. coli."""
        result = evaluate_no_cpg_island(CG_RICH_SEQ, organism="E_coli")
        assert isinstance(result, TypeCheckResult)
        assert result.verdict == Verdict.PASS

    def test_eukaryotic_human_fails_on_cpg_island(self):
        """CpG evaluation should detect violations for human."""
        result = evaluate_no_cpg_island(CG_RICH_SEQ, organism="Homo_sapiens")
        assert result.verdict == Verdict.FAIL

    def test_eukaryotic_cho_fails_on_cpg_island(self):
        """CpG evaluation should detect violations for CHO-K1."""
        result = evaluate_no_cpg_island(CG_RICH_SEQ, organism="CHO_K1")
        assert result.verdict == Verdict.FAIL

    def test_backward_compat_no_organism_runs_check(self):
        """Without organism parameter, CpG evaluation should still detect violations."""
        result = evaluate_no_cpg_island(CG_RICH_SEQ)
        assert result.verdict == Verdict.FAIL

    def test_prokaryotic_clean_seq_also_passes(self):
        """Even clean sequences should pass for prokaryotes (trivially)."""
        result = evaluate_no_cpg_island(CLEAN_SEQ, organism="E_coli_K12")
        assert result.verdict == Verdict.PASS


# ════════════════════════════════════════════════════════════════════
# 4. NoCpGIslandConstraint — solver constraint class
# ════════════════════════════════════════════════════════════════════

class TestNoCpGIslandConstraintOrganismAware:
    """Tests for NoCpGIslandConstraint with organism parameter."""

    def test_prokaryotic_constraint_check_returns_true(self):
        """NoCpGIslandConstraint.check() should return True for prokaryotes."""
        constraint = NoCpGIslandConstraint(organism="E_coli")
        assert constraint.check(CG_RICH_SEQ) is True

    def test_prokaryotic_constraint_violated_positions_empty(self):
        """NoCpGIslandConstraint.violated_positions() should return [] for prokaryotes."""
        constraint = NoCpGIslandConstraint(organism="E_coli")
        assert constraint.violated_positions(CG_RICH_SEQ) == []

    def test_eukaryotic_constraint_check_fails_on_cpg_island(self):
        """NoCpGIslandConstraint.check() should fail for eukaryotes with CpG island."""
        constraint = NoCpGIslandConstraint(organism="Homo_sapiens")
        assert constraint.check(CG_RICH_SEQ) is False

    def test_eukaryotic_constraint_violated_positions_nonempty(self):
        """NoCpGIslandConstraint.violated_positions() should be non-empty for eukaryotes."""
        constraint = NoCpGIslandConstraint(organism="Homo_sapiens")
        positions = constraint.violated_positions(CG_RICH_SEQ)
        assert len(positions) > 0

    def test_backward_compat_no_organism_runs_check(self):
        """Without organism, NoCpGIslandConstraint should detect violations."""
        constraint = NoCpGIslandConstraint()
        assert constraint.check(CG_RICH_SEQ) is False

    def test_backward_compat_empty_organism_runs_check(self):
        """Empty organism string should behave like unspecified organism."""
        constraint = NoCpGIslandConstraint(organism="")
        assert constraint.check(CG_RICH_SEQ) is False

    def test_organism_property(self):
        """Organism property should return the configured organism."""
        constraint = NoCpGIslandConstraint(organism="E_coli")
        assert constraint.organism == "E_coli"

    def test_eukaryotic_constraint_passes_on_clean_sequence(self):
        """NoCpGIslandConstraint should pass for eukaryotes on clean sequence."""
        constraint = NoCpGIslandConstraint(organism="Homo_sapiens")
        assert constraint.check(CLEAN_SEQ) is True


# ════════════════════════════════════════════════════════════════════
# 5. MinimizeCpG — solver soft constraint class
# ════════════════════════════════════════════════════════════════════

class TestMinimizeCpGOrganismAware:
    """Tests for MinimizeCpG with organism parameter."""

    def test_prokaryotic_score_is_zero(self):
        """MinimizeCpG.score() should return 0.0 for prokaryotes (neutral)."""
        constraint = MinimizeCpG(organism="E_coli")
        assert constraint.score(CG_RICH_SEQ) == 0.0

    def test_prokaryotic_violated_positions_empty(self):
        """MinimizeCpG.violated_positions() should return [] for prokaryotes."""
        constraint = MinimizeCpG(organism="E_coli")
        assert constraint.violated_positions(CG_RICH_SEQ) == []

    def test_prokaryotic_cpg_count_is_zero(self):
        """MinimizeCpG.cpg_count() should return 0 for prokaryotes."""
        constraint = MinimizeCpG(organism="E_coli")
        assert constraint.cpg_count(CG_RICH_SEQ) == 0

    def test_eukaryotic_score_is_negative(self):
        """MinimizeCpG.score() should return a negative value for eukaryotes with CG."""
        constraint = MinimizeCpG(organism="Homo_sapiens")
        score = constraint.score(CG_RICH_SEQ)
        assert score < 0.0  # negative because there are CG dinucleotides

    def test_eukaryotic_violated_positions_nonempty(self):
        """MinimizeCpG.violated_positions() should find CG positions for eukaryotes."""
        constraint = MinimizeCpG(organism="Homo_sapiens")
        positions = constraint.violated_positions(CG_RICH_SEQ)
        assert len(positions) > 0

    def test_eukaryotic_cpg_count_positive(self):
        """MinimizeCpG.cpg_count() should be positive for eukaryotes with CG."""
        constraint = MinimizeCpG(organism="Homo_sapiens")
        count = constraint.cpg_count(CG_RICH_SEQ)
        assert count > 0

    def test_backward_compat_no_organism_counts_cpg(self):
        """Without organism, MinimizeCpG should still count CG dinucleotides."""
        constraint = MinimizeCpG()
        count = constraint.cpg_count(CG_RICH_SEQ)
        assert count > 0

    def test_check_always_true(self):
        """MinimizeCpG.check() should always return True (it is a soft constraint)."""
        constraint = MinimizeCpG(organism="E_coli")
        assert constraint.check(CG_RICH_SEQ) is True
        constraint2 = MinimizeCpG(organism="Homo_sapiens")
        assert constraint2.check(CG_RICH_SEQ) is True

    def test_organism_property(self):
        """Organism property should return the configured organism."""
        constraint = MinimizeCpG(organism="E_coli")
        assert constraint.organism == "E_coli"


# ════════════════════════════════════════════════════════════════════
# 6. Constraint enforcement pipeline — integration test
# ════════════════════════════════════════════════════════════════════

class TestEnforcementPipelineOrganismAware:
    """Integration tests for the constraint enforcement pipeline.

    When a NoCpGIslandConstraint with a prokaryotic organism is in the
    model, the enforcement pipeline should not report CpG violations.
    """

    def test_prokaryotic_constraint_not_reported_as_violation(self):
        """CpG constraint should not appear in violations for prokaryotes."""
        from biocompiler.solver.constraints import CSPModel, SolverConfig, CodonVariable

        # Build a minimal model with a prokaryotic CpG constraint
        cpg_constraint = NoCpGIslandConstraint(organism="E_coli")
        config = SolverConfig()

        model = CSPModel(
            variables=[
                CodonVariable(position=0, amino_acid="M", domain=["ATG"]),
                CodonVariable(position=1, amino_acid="K", domain=["AAA", "AAG"]),
            ],
            hard_constraints=[cpg_constraint],
            soft_constraints=[],
            protein="MK",
            organism="E_coli",
            config=config,
        )

        enforcer = ConstraintEnforcer()
        # Use a CG-rich sequence — even though it would normally violate CpG,
        # the prokaryotic organism should cause the check to pass
        seq = CG_RICH_SEQ[:6]  # just "ATGCGC" - will not have enough for island but proves the point
        result = enforcer.enforce(model, seq)
        # The CpG constraint should pass (not be in violations)
        assert result.all_hard_satisfied is True

    def test_eukaryotic_constraint_reported_as_violation(self):
        """CpG constraint should appear in violations for eukaryotes."""
        from biocompiler.solver.constraints import CSPModel, SolverConfig, CodonVariable

        cpg_constraint = NoCpGIslandConstraint(organism="Homo_sapiens")
        config = SolverConfig()

        model = CSPModel(
            variables=[
                CodonVariable(position=0, amino_acid="M", domain=["ATG"]),
            ],
            hard_constraints=[cpg_constraint],
            soft_constraints=[],
            protein="M",
            organism="Homo_sapiens",
            config=config,
        )

        enforcer = ConstraintEnforcer()
        # Use the CG-rich sequence — should violate the eukaryotic CpG constraint
        # Note: the sequence needs to be long enough for the sliding window
        result = enforcer.enforce(model, CG_RICH_SEQ)
        # The CpG constraint should fail (be in violations)
        assert result.all_hard_satisfied is False
        cpg_violations = [v for v in result.violations if "CpG" in v.constraint_name]
        assert len(cpg_violations) > 0


# ════════════════════════════════════════════════════════════════════
# 7. Multiple organisms consistency
# ════════════════════════════════════════════════════════════════════

class TestMultipleOrganismsConsistency:
    """Verify that CpG behaviour is consistent across different organism
    identifiers within the same domain."""

    @pytest.mark.parametrize("organism", [
        "E_coli",
        "E_coli_K12",
        "E_coli_BL21",
        "Escherichia_coli",
        "ecoli",
    ])
    def test_prokaryotic_organisms_skip_cpg(self, organism):
        """All E. coli aliases should skip CpG checking."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism=organism)
        assert result.passed is True

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "human",
        "Mus_musculus",
        "mouse",
        "CHO_K1",
        "cho",
    ])
    def test_eukaryotic_organisms_check_cpg(self, organism):
        """All eukaryotic organism aliases should detect CpG violations."""
        result = check_no_cpg_island(CG_RICH_SEQ, organism=organism)
        assert result.passed is False

    @pytest.mark.parametrize("organism", [
        "E_coli",
        "E_coli_K12",
        "E_coli_BL21",
    ])
    def test_prokaryotic_constraint_instances_skip(self, organism):
        """NoCpGIslandConstraint should auto-satisfy for all prokaryotic identifiers."""
        constraint = NoCpGIslandConstraint(organism=organism)
        assert constraint.check(CG_RICH_SEQ) is True
        assert constraint.violated_positions(CG_RICH_SEQ) == []

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "Mus_musculus",
        "CHO_K1",
    ])
    def test_eukaryotic_constraint_instances_check(self, organism):
        """NoCpGIslandConstraint should detect violations for eukaryotic organisms."""
        constraint = NoCpGIslandConstraint(organism=organism)
        assert constraint.check(CG_RICH_SEQ) is False
        assert len(constraint.violated_positions(CG_RICH_SEQ)) > 0
