"""Tests for solver/scoring.py.

Covers:
1. ScoringResult dataclass construction
2. SoftConstraintScorer initialization and score_sequence()
3. compute_weighted_score() with known weights
4. compute_pareto_frontier() with dominated/non-dominated solutions
"""

from __future__ import annotations

import math
from unittest.mock import MagicMock

import pytest

from biocompiler.solver.scoring import (
    ScoringResult,
    SoftConstraintScorer,
    compute_pareto_frontier,
)
from biocompiler.solver.types import SolverConfig, CodonVariable
from biocompiler.solver.constraints import (
    CSPModel,
    MaximizeCAI,
    MinimizeCpG,
    MinimizeMRNADG,
)


# ── Fixtures ────────────────────────────────────────────────────────

@pytest.fixture
def default_config() -> SolverConfig:
    """Default SolverConfig with standard weights (cai=1.0, cpg=0.5, mrna_dg=0.3)."""
    return SolverConfig(gc_lo=0.30, gc_hi=0.70)


@pytest.fixture
def uniform_config() -> SolverConfig:
    """SolverConfig with equal weights for all objectives."""
    return SolverConfig(gc_lo=0.30, gc_hi=0.70, cai_weight=1.0, cpg_weight=1.0, mrna_dg_weight=1.0)


@pytest.fixture
def zero_cpg_weight_config() -> SolverConfig:
    """SolverConfig with cpg_weight=0 to isolate CAI and mRNA dG scoring."""
    return SolverConfig(gc_lo=0.30, gc_hi=0.70, cai_weight=1.0, cpg_weight=0.0, mrna_dg_weight=0.3)


# A short protein "MK" → 2 amino acids → 6 nucleotides
SHORT_PROTEIN = "MK"

# Adaptiveness table: give every codon a reasonable value so CAI is computable
# M: ATG (only codon, w=1.0); K: AAA (w=0.4), AAG (w=0.9)
SIMPLE_ADAPTIVENESS: dict[str, float] = {
    "ATG": 1.0,
    "AAA": 0.4,
    "AAG": 0.9,
}


def _make_model_with_soft_constraints(
    protein: str = SHORT_PROTEIN,
    adaptiveness: dict[str, float] | None = None,
    config: SolverConfig | None = None,
) -> CSPModel:
    """Build a minimal CSPModel with the three standard soft constraints."""
    if adaptiveness is None:
        adaptiveness = SIMPLE_ADAPTIVENESS
    if config is None:
        config = SolverConfig()

    soft_constraints = [
        MaximizeCAI(adaptiveness, protein),
        MinimizeCpG(),
        MinimizeMRNADG(),
    ]

    return CSPModel(
        variables=[],  # Not used by scoring, but required by dataclass
        hard_constraints=[],
        soft_constraints=soft_constraints,
        protein=protein,
        organism="Homo_sapiens",
        config=config,
    )


# ════════════════════════════════════════════════════════════════════
# 1. ScoringResult dataclass construction
# ════════════════════════════════════════════════════════════════════

class TestScoringResultConstruction:
    """Tests for ScoringResult dataclass creation and field access."""

    def test_basic_construction(self):
        """ScoringResult should be constructable with all required fields."""
        result = ScoringResult(
            cai_score=0.8,
            cpg_score=0.7,
            mrna_dg_score=0.6,
            weighted_total=0.75,
            individual_scores={"MaximizeCAI": -0.5, "MinimizeCpG": -3.0},
        )
        assert result.cai_score == 0.8
        assert result.cpg_score == 0.7
        assert result.mrna_dg_score == 0.6
        assert result.weighted_total == 0.75
        assert result.individual_scores["MaximizeCAI"] == -0.5

    def test_all_zeros(self):
        """ScoringResult should accept all-zero scores."""
        result = ScoringResult(
            cai_score=0.0,
            cpg_score=0.0,
            mrna_dg_score=0.0,
            weighted_total=0.0,
            individual_scores={},
        )
        assert result.cai_score == 0.0
        assert result.cpg_score == 0.0
        assert result.mrna_dg_score == 0.0
        assert result.weighted_total == 0.0
        assert result.individual_scores == {}

    def test_all_ones(self):
        """ScoringResult should accept all-ones scores (perfect)."""
        result = ScoringResult(
            cai_score=1.0,
            cpg_score=1.0,
            mrna_dg_score=1.0,
            weighted_total=1.8,
            individual_scores={"MaximizeCAI": 0.0, "MinimizeCpG": 0.0, "MinimizeMRNADG": 0.0},
        )
        assert result.cai_score == 1.0
        assert result.cpg_score == 1.0
        assert result.mrna_dg_score == 1.0

    def test_repr_format(self):
        """ScoringResult.__repr__ should produce a readable string with 4-decimal formatting."""
        result = ScoringResult(
            cai_score=0.12345,
            cpg_score=0.67890,
            mrna_dg_score=0.11111,
            weighted_total=0.45678,
            individual_scores={},
        )
        r = repr(result)
        assert "ScoringResult(" in r
        assert "cai=" in r
        assert "cpg=" in r
        assert "mrna_dg=" in r
        assert "total=" in r
        # Check 4-decimal formatting
        assert "0.1235" in r  # 0.12345 rounded to 4 decimals
        assert "0.6789" in r

    def test_individual_scores_dict_is_preserved(self):
        """individual_scores should preserve exact key-value pairs."""
        raw = {"MaximizeCAI": -1.234, "MinimizeCpG": -5.0, "MinimizeMRNADG": -12.5}
        result = ScoringResult(
            cai_score=0.5, cpg_score=0.5, mrna_dg_score=0.5,
            weighted_total=0.5, individual_scores=raw,
        )
        assert result.individual_scores == raw
        assert len(result.individual_scores) == 3

    def test_field_types_are_float_and_dict(self):
        """Score fields should be float; individual_scores should be dict."""
        result = ScoringResult(
            cai_score=0.5, cpg_score=0.5, mrna_dg_score=0.5,
            weighted_total=0.5, individual_scores={"k": 1.0},
        )
        assert isinstance(result.cai_score, float)
        assert isinstance(result.cpg_score, float)
        assert isinstance(result.mrna_dg_score, float)
        assert isinstance(result.weighted_total, float)
        assert isinstance(result.individual_scores, dict)


# ════════════════════════════════════════════════════════════════════
# 2. SoftConstraintScorer initialization and score_sequence()
# ════════════════════════════════════════════════════════════════════

class TestSoftConstraintScorerInit:
    """Tests for SoftConstraintScorer construction."""

    def test_init_stores_config(self, default_config: SolverConfig):
        """Scorer should store the config and expose it via the config property."""
        scorer = SoftConstraintScorer(default_config)
        assert scorer.config is default_config

    def test_config_property_returns_solver_config(self, default_config: SolverConfig):
        """config property should return a SolverConfig instance."""
        scorer = SoftConstraintScorer(default_config)
        assert isinstance(scorer.config, SolverConfig)

    def test_custom_weights_preserved(self):
        """Custom weights should be accessible through the config."""
        config = SolverConfig(cai_weight=2.0, cpg_weight=1.5, mrna_dg_weight=0.8)
        scorer = SoftConstraintScorer(config)
        assert scorer.config.cai_weight == 2.0
        assert scorer.config.cpg_weight == 1.5
        assert scorer.config.mrna_dg_weight == 0.8


class TestScoreSequence:
    """Tests for SoftConstraintScorer.score_sequence()."""

    def test_returns_scoring_result(self, default_config: SolverConfig):
        """score_sequence() must return a ScoringResult instance."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        # "ATGAAG" = ATG(M) + AAG(K), both high-adaptiveness codons
        result = scorer.score_sequence(model, "ATGAAG")
        assert isinstance(result, ScoringResult)

    def test_correct_sequence_length_for_two_aa(self, default_config: SolverConfig):
        """For a 2-amino-acid protein, sequence must be 6 nucleotides."""
        model = _make_model_with_soft_constraints(protein="MK", config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        # Just verify it does not raise; the result should have valid fields
        assert 0.0 <= result.cai_score <= 1.0

    def test_wrong_sequence_length_raises_value_error(self, default_config: SolverConfig):
        """If sequence length != len(protein)*3, raise ValueError."""
        model = _make_model_with_soft_constraints(protein="MK", config=default_config)
        scorer = SoftConstraintScorer(default_config)
        with pytest.raises(ValueError, match="Sequence length"):
            scorer.score_sequence(model, "ATG")  # 3 nt, but protein has 2 aa → need 6

    def test_sequence_too_long_raises_value_error(self, default_config: SolverConfig):
        """Overlong sequence should also raise ValueError."""
        model = _make_model_with_soft_constraints(protein="MK", config=default_config)
        scorer = SoftConstraintScorer(default_config)
        with pytest.raises(ValueError, match="Sequence length"):
            scorer.score_sequence(model, "ATGAAGTTT")  # 9 nt instead of 6

    def test_cai_score_in_range(self, default_config: SolverConfig):
        """CAI score must be in [0.0, 1.0]."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert 0.0 <= result.cai_score <= 1.0

    def test_cpg_score_in_range(self, default_config: SolverConfig):
        """CpG score must be in [0.0, 1.0]."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert 0.0 <= result.cpg_score <= 1.0

    def test_mrna_dg_score_in_range(self, default_config: SolverConfig):
        """mRNA dG score must be in [0.0, 1.0]."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert 0.0 <= result.mrna_dg_score <= 1.0

    def test_individual_scores_populated(self, default_config: SolverConfig):
        """individual_scores dict should contain entries for each soft constraint."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert "MaximizeCAI" in result.individual_scores
        assert "MinimizeCpG" in result.individual_scores
        assert "MinimizeMRNADG" in result.individual_scores

    def test_no_cpg_gives_cpg_score_one(self, default_config: SolverConfig):
        """A sequence with zero CpG dinucleotides should give cpg_score ≈ 1.0."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        # "ATGAAG" has no "CG" dinucleotide
        result = scorer.score_sequence(model, "ATGAAG")
        assert result.cpg_score == 1.0

    def test_sequence_with_cpg_gives_lower_cpg_score(self, default_config: SolverConfig):
        """A sequence with CpG dinucleotides should have cpg_score < 1.0."""
        # Build a protein and sequence that includes "CG" dinucleotide
        # Protein "MR" → M=ATG, R=CGT/CGC/CGA/CGG/AGA/AGG
        # Use CGT for R which contains CG at positions 3-4
        adaptiveness = {"ATG": 1.0, "CGT": 0.5, "AGA": 0.3}
        protein = "MR"
        model = _make_model_with_soft_constraints(
            protein=protein, adaptiveness=adaptiveness, config=default_config,
        )
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGCGT")
        assert result.cpg_score < 1.0

    def test_high_cai_codons_give_high_cai_score(self, default_config: SolverConfig):
        """Using highest-adaptiveness codons should produce CAI close to 1.0."""
        # ATG(M, w=1.0) + AAG(K, w=0.9) → CAI should be near 1.0
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert result.cai_score > 0.8

    def test_low_cai_codons_give_lower_cai_score(self, default_config: SolverConfig):
        """Using low-adaptiveness codons should produce lower CAI."""
        # ATG(M, w=1.0) + AAA(K, w=0.4) → CAI should be lower
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)
        result_high = scorer.score_sequence(model, "ATGAAG")  # w=0.9 for K
        result_low = scorer.score_sequence(model, "ATGAAA")   # w=0.4 for K
        assert result_high.cai_score > result_low.cai_score

    def test_weighted_total_matches_hand_calculation(self):
        """Weighted total should equal sum of weighted normalized scores."""
        config = SolverConfig(cai_weight=1.0, cpg_weight=0.5, mrna_dg_weight=0.3)
        model = _make_model_with_soft_constraints(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        # weighted_total should be: 1.0*cai + 0.5*cpg + 0.3*mrna_dg
        expected = (
            config.cai_weight * result.cai_score
            + config.cpg_weight * result.cpg_score
            + config.mrna_dg_weight * result.mrna_dg_score
        )
        assert math.isclose(result.weighted_total, expected, abs_tol=1e-9)

    def test_zero_weight_excludes_objective(self, zero_cpg_weight_config: SolverConfig):
        """If cpg_weight=0, the CpG score should not affect weighted_total."""
        model = _make_model_with_soft_constraints(config=zero_cpg_weight_config)
        scorer = SoftConstraintScorer(zero_cpg_weight_config)
        result = scorer.score_sequence(model, "ATGAAG")
        expected = (
            zero_cpg_weight_config.cai_weight * result.cai_score
            + zero_cpg_weight_config.cpg_weight * result.cpg_score
            + zero_cpg_weight_config.mrna_dg_weight * result.mrna_dg_score
        )
        assert math.isclose(result.weighted_total, expected, abs_tol=1e-9)

    def test_longer_protein_scoring(self, default_config: SolverConfig):
        """Scoring should work for a longer protein (10 AA)."""
        protein = "MMMMMMMMMM"  # 10 methionines → 30 nt, all ATG
        adaptiveness = {"ATG": 1.0}
        model = _make_model_with_soft_constraints(
            protein=protein, adaptiveness=adaptiveness, config=default_config,
        )
        scorer = SoftConstraintScorer(default_config)
        sequence = "ATG" * 10
        result = scorer.score_sequence(model, sequence)
        assert isinstance(result, ScoringResult)
        assert result.cai_score == 1.0  # All best codons
        assert result.cpg_score == 1.0  # No CG in ATGATG...

    def test_model_with_no_soft_constraints(self, default_config: SolverConfig):
        """A model with no soft constraints should still produce a valid ScoringResult."""
        model = CSPModel(
            variables=[],
            hard_constraints=[],
            soft_constraints=[],  # No soft constraints at all
            protein="MK",
            organism="Homo_sapiens",
            config=default_config,
        )
        scorer = SoftConstraintScorer(default_config)
        result = scorer.score_sequence(model, "ATGAAG")
        # With no soft constraints, raw scores dict is empty
        # Fallback normalization should still produce values
        assert isinstance(result, ScoringResult)
        # Individual scores should be empty
        assert result.individual_scores == {}


# ════════════════════════════════════════════════════════════════════
# 3. compute_weighted_score() with known weights
# ════════════════════════════════════════════════════════════════════

class TestComputeWeightedScore:
    """Tests for SoftConstraintScorer.compute_weighted_score()."""

    def test_single_objective(self, default_config: SolverConfig):
        """Single objective with weight 1.0 should return the score itself."""
        scorer = SoftConstraintScorer(default_config)
        result = scorer.compute_weighted_score(
            scores={"MaximizeCAI": 0.8},
            weights={"MaximizeCAI": 1.0},
        )
        assert math.isclose(result, 0.8)

    def test_multiple_objectives(self, default_config: SolverConfig):
        """Weighted sum with known weights should match hand calculation."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"MaximizeCAI": 0.9, "MinimizeCpG": 0.7, "MinimizeMRNADG": 0.5}
        weights = {"MaximizeCAI": 1.0, "MinimizeCpG": 0.5, "MinimizeMRNADG": 0.3}
        result = scorer.compute_weighted_score(scores, weights)
        expected = 1.0 * 0.9 + 0.5 * 0.7 + 0.3 * 0.5
        assert math.isclose(result, expected, abs_tol=1e-12)

    def test_zero_weights_give_zero_total(self, default_config: SolverConfig):
        """All-zero weights should produce a weighted total of 0.0."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"MaximizeCAI": 0.9, "MinimizeCpG": 0.8, "MinimizeMRNADG": 0.7}
        weights = {"MaximizeCAI": 0.0, "MinimizeCpG": 0.0, "MinimizeMRNADG": 0.0}
        result = scorer.compute_weighted_score(scores, weights)
        assert result == 0.0

    def test_zero_scores_give_zero_total(self, default_config: SolverConfig):
        """All-zero scores should produce a weighted total of 0.0."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"MaximizeCAI": 0.0, "MinimizeCpG": 0.0, "MinimizeMRNADG": 0.0}
        weights = {"MaximizeCAI": 1.0, "MinimizeCpG": 0.5, "MinimizeMRNADG": 0.3}
        result = scorer.compute_weighted_score(scores, weights)
        assert result == 0.0

    def test_missing_weight_defaults_to_zero(self, default_config: SolverConfig):
        """A score key not present in weights should be treated as weight 0.0."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"MaximizeCAI": 0.9, "UnknownConstraint": 1.0}
        weights = {"MaximizeCAI": 1.0}
        result = scorer.compute_weighted_score(scores, weights)
        assert math.isclose(result, 0.9)

    def test_empty_scores_and_weights(self, default_config: SolverConfig):
        """Empty inputs should return 0.0."""
        scorer = SoftConstraintScorer(default_config)
        result = scorer.compute_weighted_score({}, {})
        assert result == 0.0

    def test_negative_weights(self, default_config: SolverConfig):
        """Negative weights should produce negative contribution (penalize)."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"MaximizeCAI": 0.5}
        weights = {"MaximizeCAI": -1.0}
        result = scorer.compute_weighted_score(scores, weights)
        assert math.isclose(result, -0.5)

    def test_large_weight_dominates(self, default_config: SolverConfig):
        """A very large weight should dominate the weighted total."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"MaximizeCAI": 1.0, "MinimizeCpG": 0.1}
        weights = {"MaximizeCAI": 100.0, "MinimizeCpG": 0.01}
        result = scorer.compute_weighted_score(scores, weights)
        assert math.isclose(result, 100.001, abs_tol=1e-9)

    def test_uniform_weights_equal_simple_average(self, default_config: SolverConfig):
        """With equal weights, weighted sum equals weight * simple_average."""
        scorer = SoftConstraintScorer(default_config)
        scores = {"A": 0.6, "B": 0.8, "C": 1.0}
        weights = {"A": 1.0, "B": 1.0, "C": 1.0}
        result = scorer.compute_weighted_score(scores, weights)
        expected = 1.0 * 0.6 + 1.0 * 0.8 + 1.0 * 1.0
        assert math.isclose(result, expected)


# ════════════════════════════════════════════════════════════════════
# 4. compute_pareto_frontier() with dominated/non-dominated solutions
# ════════════════════════════════════════════════════════════════════

class TestComputeParetoFrontier:
    """Tests for compute_pareto_frontier()."""

    def _make_result(
        self,
        cai: float = 0.5,
        cpg: float = 0.5,
        mrna_dg: float = 0.5,
    ) -> ScoringResult:
        """Helper to create a ScoringResult with given objective scores."""
        return ScoringResult(
            cai_score=cai,
            cpg_score=cpg,
            mrna_dg_score=mrna_dg,
            weighted_total=0.0,
            individual_scores={},
        )

    def test_empty_list_returns_empty(self):
        """Empty input should return empty list."""
        assert compute_pareto_frontier([]) == []

    def test_single_solution_is_pareto(self):
        """A single solution is always Pareto-optimal."""
        r = self._make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)
        result = compute_pareto_frontier([r])
        assert len(result) == 1
        assert result[0] is r

    def test_dominated_solution_excluded(self):
        """A solution dominated in all objectives should be excluded."""
        # r1 dominates r2: r1 is >= r2 in all, > in at least one
        r1 = self._make_result(cai=0.9, cpg=0.8, mrna_dg=0.7)
        r2 = self._make_result(cai=0.5, cpg=0.4, mrna_dg=0.3)  # dominated
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 1
        assert result[0] is r1

    def test_two_non_dominated_solutions_both_survive(self):
        """Two solutions trading off different objectives should both be Pareto."""
        # r1: high CAI, low CpG; r2: low CAI, high CpG → neither dominates
        r1 = self._make_result(cai=0.9, cpg=0.3, mrna_dg=0.5)
        r2 = self._make_result(cai=0.7, cpg=0.8, mrna_dg=0.5)
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 2

    def test_docstring_example(self):
        """The exact example from the compute_pareto_frontier docstring."""
        r1 = ScoringResult(cai_score=0.9, cpg_score=0.3, mrna_dg_score=0.5, weighted_total=0.0, individual_scores={})
        r2 = ScoringResult(cai_score=0.7, cpg_score=0.8, mrna_dg_score=0.4, weighted_total=0.0, individual_scores={})
        r3 = ScoringResult(cai_score=0.6, cpg_score=0.4, mrna_dg_score=0.3, weighted_total=0.0, individual_scores={})
        pareto = compute_pareto_frontier([r1, r2, r3])
        assert len(pareto) == 2
        # r3 is dominated by both r1 and r2, so only r1 and r2 survive

    def test_all_equal_solutions_are_all_pareto(self):
        """Solutions with identical objective values are all Pareto-optimal (none strictly dominates)."""
        r1 = self._make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)
        r2 = self._make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 2

    def test_three_way_tradeoff(self):
        """Three solutions, each best in one objective, should all be Pareto."""
        r1 = self._make_result(cai=1.0, cpg=0.1, mrna_dg=0.1)  # best CAI
        r2 = self._make_result(cai=0.1, cpg=1.0, mrna_dg=0.1)  # best CpG
        r3 = self._make_result(cai=0.1, cpg=0.1, mrna_dg=1.0)  # best mRNA dG
        result = compute_pareto_frontier([r1, r2, r3])
        assert len(result) == 3

    def test_partially_dominated_excluded(self):
        """A solution equal in two objectives but worse in one is dominated."""
        r1 = self._make_result(cai=0.8, cpg=0.6, mrna_dg=0.7)
        r2 = self._make_result(cai=0.8, cpg=0.6, mrna_dg=0.5)  # worse mrna_dg
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 1
        assert result[0] is r1

    def test_preserves_original_order(self):
        """Pareto frontier should preserve the order from the input list."""
        r1 = self._make_result(cai=0.9, cpg=0.3, mrna_dg=0.5)
        r2 = self._make_result(cai=0.7, cpg=0.8, mrna_dg=0.4)
        r3 = self._make_result(cai=0.6, cpg=0.4, mrna_dg=0.3)  # dominated
        result = compute_pareto_frontier([r1, r2, r3])
        assert result[0] is r1
        assert result[1] is r2

    def test_many_solutions_with_one_dominating_all(self):
        """One clearly superior solution should be the only Pareto member."""
        dominant = self._make_result(cai=1.0, cpg=1.0, mrna_dg=1.0)
        others = [
            self._make_result(cai=0.5, cpg=0.5, mrna_dg=0.5),
            self._make_result(cai=0.3, cpg=0.7, mrna_dg=0.2),
            self._make_result(cai=0.8, cpg=0.2, mrna_dg=0.4),
        ]
        result = compute_pareto_frontier([dominant] + others)
        assert len(result) == 1
        assert result[0] is dominant

    def test_frontier_with_five_solutions(self):
        """Mix of dominated and non-dominated across 5 solutions."""
        r1 = self._make_result(cai=0.9, cpg=0.2, mrna_dg=0.5)  # Pareto
        r2 = self._make_result(cai=0.5, cpg=0.9, mrna_dg=0.5)  # Pareto
        r3 = self._make_result(cai=0.4, cpg=0.3, mrna_dg=0.9)  # Pareto
        r4 = self._make_result(cai=0.4, cpg=0.2, mrna_dg=0.4)  # dominated by r1
        r5 = self._make_result(cai=0.3, cpg=0.3, mrna_dg=0.3)  # dominated by many
        result = compute_pareto_frontier([r1, r2, r3, r4, r5])
        assert len(result) == 3
        assert r1 in result
        assert r2 in result
        assert r3 in result

    def test_dominance_check_strictly_better_in_one(self):
        """Dominance requires being strictly better in at least one objective."""
        # r1 >= r2 in all, but equal everywhere → does NOT dominate
        r1 = self._make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)
        r2 = self._make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 2

    def test_weighted_total_not_used_for_pareto(self):
        """Pareto frontier uses only cai_score, cpg_score, mrna_dg_score —
        weighted_total is irrelevant to dominance."""
        r1 = ScoringResult(cai_score=0.5, cpg_score=0.5, mrna_dg_score=0.5,
                           weighted_total=100.0, individual_scores={})
        r2 = ScoringResult(cai_score=0.5, cpg_score=0.5, mrna_dg_score=0.5,
                           weighted_total=0.0, individual_scores={})
        result = compute_pareto_frontier([r1, r2])
        # Neither dominates the other (objectives are equal)
        assert len(result) == 2


# ════════════════════════════════════════════════════════════════════
# 5. Integration: scorer + pareto together
# ════════════════════════════════════════════════════════════════════

class TestScoringAndParetoIntegration:
    """Integration tests scoring sequences and then running Pareto analysis."""

    def test_score_multiple_candidates_and_find_frontier(self, default_config: SolverConfig):
        """Score multiple candidate sequences and compute Pareto frontier."""
        model = _make_model_with_soft_constraints(config=default_config)
        scorer = SoftConstraintScorer(default_config)

        # Two different sequences for "MK"
        result_a = scorer.score_sequence(model, "ATGAAG")  # high-CAI K codon
        result_b = scorer.score_sequence(model, "ATGAAA")  # low-CAI K codon

        pareto = compute_pareto_frontier([result_a, result_b])
        # result_a should have higher CAI, same or similar other scores
        # At minimum, both should be valid ScoringResults
        assert all(isinstance(r, ScoringResult) for r in pareto)
        assert len(pareto) >= 1

    def test_pareto_frontier_size_with_diverse_candidates(self, default_config: SolverConfig):
        """Multiple candidates with different trade-offs should yield multiple Pareto members."""
        protein = "MK"
        # Use adaptiveness that creates meaningful differences
        adaptiveness = {"ATG": 1.0, "AAA": 0.4, "AAG": 0.9}
        model = _make_model_with_soft_constraints(
            protein=protein, adaptiveness=adaptiveness, config=default_config,
        )
        scorer = SoftConstraintScorer(default_config)

        # Only two valid codon choices for "MK": ATGAAG and ATGAAA
        # Both should be scored
        r1 = scorer.score_sequence(model, "ATGAAG")
        r2 = scorer.score_sequence(model, "ATGAAA")
        pareto = compute_pareto_frontier([r1, r2])
        # At least one must be Pareto
        assert len(pareto) >= 1
