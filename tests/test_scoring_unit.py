"""Unit tests for solver/scoring.py — beyond the basic test.

Covers:
1. SoftConstraintScorer with custom weights
2. compute_pareto_frontier edge cases (empty, single, all-dominated)
3. ScoringResult field validation
"""

from __future__ import annotations

import math
from dataclasses import fields

import pytest

from biocompiler.solver.scoring import (
    ScoringResult,
    SoftConstraintScorer,
    compute_pareto_frontier,
)
from biocompiler.solver.types import SolverConfig
from biocompiler.solver.constraints import (
    CSPModel,
    MaximizeCAI,
    MinimizeCpG,
    MinimizeMRNADG,
)


# ── Helpers ────────────────────────────────────────────────────────

SHORT_PROTEIN = "MK"

SIMPLE_ADAPTIVENESS: dict[str, float] = {
    "ATG": 1.0,
    "AAA": 0.4,
    "AAG": 0.9,
}


def _make_model(
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
        variables=[],
        hard_constraints=[],
        soft_constraints=soft_constraints,
        protein=protein,
        organism="Homo_sapiens",
        config=config,
    )


def _make_result(
    cai: float = 0.5,
    cpg: float = 0.5,
    mrna_dg: float = 0.5,
    weighted_total: float = 0.0,
    individual_scores: dict[str, float] | None = None,
) -> ScoringResult:
    """Shorthand to create a ScoringResult for testing."""
    return ScoringResult(
        cai_score=cai,
        cpg_score=cpg,
        mrna_dg_score=mrna_dg,
        weighted_total=weighted_total,
        individual_scores=individual_scores if individual_scores is not None else {},
    )


# ════════════════════════════════════════════════════════════════════
# 1. SoftConstraintScorer with custom weights
# ════════════════════════════════════════════════════════════════════

class TestSoftConstraintScorerCustomWeights:
    """Tests for SoftConstraintScorer under non-default weight configurations."""

    def test_cai_only_weight(self):
        """When only cai_weight is nonzero, weighted_total should equal cai_score."""
        config = SolverConfig(cai_weight=1.0, cpg_weight=0.0, mrna_dg_weight=0.0)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        expected = 1.0 * result.cai_score
        assert math.isclose(result.weighted_total, expected, abs_tol=1e-9)

    def test_cpg_only_weight(self):
        """When only cpg_weight is nonzero, weighted_total should equal cpg_score * weight."""
        config = SolverConfig(cai_weight=0.0, cpg_weight=2.5, mrna_dg_weight=0.0)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        expected = 2.5 * result.cpg_score
        assert math.isclose(result.weighted_total, expected, abs_tol=1e-9)

    def test_mrna_dg_only_weight(self):
        """When only mrna_dg_weight is nonzero, weighted_total tracks that objective."""
        config = SolverConfig(cai_weight=0.0, cpg_weight=0.0, mrna_dg_weight=3.0)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        expected = 3.0 * result.mrna_dg_score
        assert math.isclose(result.weighted_total, expected, abs_tol=1e-9)

    def test_all_weights_zero(self):
        """All weights zero should produce weighted_total == 0.0 regardless of scores."""
        config = SolverConfig(cai_weight=0.0, cpg_weight=0.0, mrna_dg_weight=0.0)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert result.weighted_total == 0.0

    def test_large_weight_dominates_total(self):
        """A very large weight should dominate the weighted total."""
        config = SolverConfig(cai_weight=1000.0, cpg_weight=0.01, mrna_dg_weight=0.01)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        # The total should be dominated by 1000 * cai_score
        cai_contribution = 1000.0 * result.cai_score
        other_contribution = 0.01 * result.cpg_score + 0.01 * result.mrna_dg_score
        assert cai_contribution > 100 * other_contribution  # cai dominates by >100x

    def test_negative_weight_reduces_total(self):
        """A negative weight should subtract from the weighted total."""
        config_pos = SolverConfig(cai_weight=1.0, cpg_weight=0.5, mrna_dg_weight=0.3)
        config_neg = SolverConfig(cai_weight=1.0, cpg_weight=-0.5, mrna_dg_weight=0.3)
        model = _make_model()
        scorer_pos = SoftConstraintScorer(config_pos)
        scorer_neg = SoftConstraintScorer(config_neg)
        result_pos = scorer_pos.score_sequence(model, "ATGAAG")
        result_neg = scorer_neg.score_sequence(model, "ATGAAG")
        # Negative cpg_weight should reduce total (cpg_score is always >= 0)
        assert result_neg.weighted_total < result_pos.weighted_total

    def test_uniform_weights_proportional_to_sum(self):
        """Equal weights should make weighted_total = w * (cai + cpg + mrna_dg)."""
        w = 0.7
        config = SolverConfig(cai_weight=w, cpg_weight=w, mrna_dg_weight=w)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        expected = w * (result.cai_score + result.cpg_score + result.mrna_dg_score)
        assert math.isclose(result.weighted_total, expected, abs_tol=1e-9)

    def test_custom_weights_with_different_sequences(self):
        """Custom weights should correctly differentiate sequences."""
        config = SolverConfig(cai_weight=2.0, cpg_weight=0.0, mrna_dg_weight=0.0)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        # ATGAAG has high CAI; ATGAAA has low CAI
        result_high = scorer.score_sequence(model, "ATGAAG")
        result_low = scorer.score_sequence(model, "ATGAAA")
        assert result_high.weighted_total > result_low.weighted_total

    def test_weight_map_keys_match_constraint_names(self):
        """_build_weight_map() should return keys matching soft constraint names."""
        config = SolverConfig()
        scorer = SoftConstraintScorer(config)
        weight_map = scorer._build_weight_map()
        assert set(weight_map.keys()) == {"MaximizeCAI", "MinimizeCpG", "MinimizeMRNADG"}
        assert weight_map["MaximizeCAI"] == config.cai_weight
        assert weight_map["MinimizeCpG"] == config.cpg_weight
        assert weight_map["MinimizeMRNADG"] == config.mrna_dg_weight

    def test_compute_weighted_score_with_extra_keys_in_scores(self):
        """Extra keys in scores dict (not in weights) should be ignored (weight=0)."""
        config = SolverConfig()
        scorer = SoftConstraintScorer(config)
        scores = {"MaximizeCAI": 0.8, "UnknownObjective": 0.9}
        weights = {"MaximizeCAI": 1.0}
        result = scorer.compute_weighted_score(scores, weights)
        assert math.isclose(result, 0.8)

    def test_scorer_reused_across_sequences(self):
        """A single scorer instance should correctly score multiple sequences."""
        config = SolverConfig(cai_weight=1.0, cpg_weight=0.5, mrna_dg_weight=0.3)
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        r1 = scorer.score_sequence(model, "ATGAAG")
        r2 = scorer.score_sequence(model, "ATGAAA")
        assert isinstance(r1, ScoringResult)
        assert isinstance(r2, ScoringResult)
        assert r1.cai_score > r2.cai_score


# ════════════════════════════════════════════════════════════════════
# 2. compute_pareto_frontier edge cases
# ════════════════════════════════════════════════════════════════════

class TestParetoFrontierEdgeCases:
    """Edge-case tests for compute_pareto_frontier beyond basic dominance."""

    def test_empty_input_returns_empty_list(self):
        """Empty list input should return an empty list."""
        result = compute_pareto_frontier([])
        assert result == []
        assert isinstance(result, list)

    def test_single_solution_is_always_pareto(self):
        """A single solution can never be dominated — it is Pareto-optimal."""
        r = _make_result(cai=0.1, cpg=0.1, mrna_dg=0.1)
        result = compute_pareto_frontier([r])
        assert len(result) == 1
        assert result[0] is r

    def test_all_dominated_except_one(self):
        """When one solution dominates all others, only it survives."""
        dominant = _make_result(cai=1.0, cpg=1.0, mrna_dg=1.0)
        weak = [
            _make_result(cai=0.5, cpg=0.5, mrna_dg=0.5),
            _make_result(cai=0.3, cpg=0.7, mrna_dg=0.2),
            _make_result(cai=0.9, cpg=0.1, mrna_dg=0.3),
        ]
        result = compute_pareto_frontier([dominant] + weak)
        assert len(result) == 1
        assert result[0] is dominant

    def test_solution_dominated_by_one_but_not_another(self):
        """A solution may be dominated by one solution but not another.
        Only the dominance relationship matters, not how many dominate it."""
        # r1 dominates r3 (>= in all, > in CAI and mrna_dG)
        # r2 does NOT dominate r3 (r2.cai < r3.cai)
        # But r3 is still excluded because r1 dominates it
        r1 = _make_result(cai=0.9, cpg=0.5, mrna_dg=0.8)
        r2 = _make_result(cai=0.3, cpg=0.9, mrna_dg=0.8)
        r3 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)  # dominated by r1
        result = compute_pareto_frontier([r1, r2, r3])
        assert r3 not in result
        assert r1 in result
        assert r2 in result

    def test_all_solutions_dominated_by_at_least_one_other(self):
        """When every solution is dominated by at least one other, the
        frontier contains only the non-dominated ones."""
        # Chain: r1 > r2 > r3 > r4 in all objectives
        r1 = _make_result(cai=1.0, cpg=1.0, mrna_dg=1.0)
        r2 = _make_result(cai=0.8, cpg=0.8, mrna_dg=0.8)
        r3 = _make_result(cai=0.6, cpg=0.6, mrna_dg=0.6)
        r4 = _make_result(cai=0.4, cpg=0.4, mrna_dg=0.4)
        result = compute_pareto_frontier([r1, r2, r3, r4])
        assert len(result) == 1
        assert result[0] is r1

    def test_two_identical_solutions_both_pareto(self):
        """Identical solutions don't dominate each other (need strictly better)."""
        r1 = _make_result(cai=0.7, cpg=0.7, mrna_dg=0.7)
        r2 = _make_result(cai=0.7, cpg=0.7, mrna_dg=0.7)
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 2

    def test_pareto_with_boundary_scores_zero_and_one(self):
        """Solutions at score boundaries (0.0 and 1.0) should be handled correctly."""
        r_perfect = _make_result(cai=1.0, cpg=1.0, mrna_dg=1.0)
        r_zero = _make_result(cai=0.0, cpg=0.0, mrna_dg=0.0)
        result = compute_pareto_frontier([r_perfect, r_zero])
        assert len(result) == 1
        assert result[0] is r_perfect

    def test_pareto_with_mixed_boundary_scores(self):
        """Solutions with some 0.0 and some 1.0 scores — trade-off scenario."""
        r1 = _make_result(cai=1.0, cpg=0.0, mrna_dg=0.5)
        r2 = _make_result(cai=0.0, cpg=1.0, mrna_dg=0.5)
        r3 = _make_result(cai=0.0, cpg=0.0, mrna_dg=0.5)  # dominated by both
        result = compute_pareto_frontier([r1, r2, r3])
        assert len(result) == 2
        assert r1 in result
        assert r2 in result
        assert r3 not in result

    def test_frontier_preserves_input_order(self):
        """Pareto-optimal solutions should appear in their original input order."""
        r_a = _make_result(cai=0.9, cpg=0.2, mrna_dg=0.5)
        r_b = _make_result(cai=0.2, cpg=0.9, mrna_dg=0.5)
        r_c = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.9)
        # All three trade off, none dominates another
        result = compute_pareto_frontier([r_a, r_b, r_c])
        assert result == [r_a, r_b, r_c]

    def test_many_duplicates_all_pareto(self):
        """Multiple copies of the same solution should all be Pareto-optimal."""
        r = _make_result(cai=0.6, cpg=0.6, mrna_dg=0.6)
        result = compute_pareto_frontier([r, r, r])
        assert len(result) == 3

    def test_dominance_tie_in_two_objectives_strictly_better_in_third(self):
        """If equal in two objectives but strictly better in the third, dominates."""
        r1 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.9)
        r2 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.3)
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 1
        assert result[0] is r1

    def test_all_equal_on_two_objectives_vary_third(self):
        """Solutions equal in CAI and CpG but varying in mrna_dG: only the best dG survives."""
        r1 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.3)
        r2 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.7)
        r3 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.5)
        result = compute_pareto_frontier([r1, r2, r3])
        assert len(result) == 1
        assert result[0] is r2

    def test_weighted_total_does_not_affect_pareto(self):
        """Pareto frontier should ignore weighted_total — only the three
        individual objective scores matter for dominance."""
        r1 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.5, weighted_total=999.0)
        r2 = _make_result(cai=0.5, cpg=0.5, mrna_dg=0.5, weighted_total=0.0)
        result = compute_pareto_frontier([r1, r2])
        assert len(result) == 2  # neither dominates the other


# ════════════════════════════════════════════════════════════════════
# 3. ScoringResult field validation
# ════════════════════════════════════════════════════════════════════

class TestScoringResultFieldValidation:
    """Tests for ScoringResult dataclass field constraints and validation."""

    def test_has_exactly_five_fields(self):
        """ScoringResult should have exactly the 5 expected fields."""
        field_names = {f.name for f in fields(ScoringResult)}
        expected = {"cai_score", "cpg_score", "mrna_dg_score", "weighted_total", "individual_scores"}
        assert field_names == expected

    def test_all_score_fields_are_float(self):
        """Score fields should accept and store float values."""
        result = _make_result(cai=0.1, cpg=0.2, mrna_dg=0.3, weighted_total=0.6)
        assert isinstance(result.cai_score, float)
        assert isinstance(result.cpg_score, float)
        assert isinstance(result.mrna_dg_score, float)
        assert isinstance(result.weighted_total, float)

    def test_individual_scores_is_dict(self):
        """individual_scores should be a dict."""
        result = _make_result(individual_scores={"MaximizeCAI": -0.5})
        assert isinstance(result.individual_scores, dict)

    def test_individual_scores_can_be_empty(self):
        """individual_scores should accept an empty dict."""
        result = _make_result(individual_scores={})
        assert result.individual_scores == {}

    def test_individual_scores_preserves_all_entries(self):
        """individual_scores should preserve all key-value pairs exactly."""
        scores = {
            "MaximizeCAI": -1.234,
            "MinimizeCpG": -5.0,
            "MinimizeMRNADG": -12.5,
        }
        result = _make_result(individual_scores=scores)
        assert result.individual_scores == scores
        assert len(result.individual_scores) == 3

    def test_score_fields_accept_zero(self):
        """All score fields should accept 0.0."""
        result = _make_result(cai=0.0, cpg=0.0, mrna_dg=0.0, weighted_total=0.0)
        assert result.cai_score == 0.0
        assert result.cpg_score == 0.0
        assert result.mrna_dg_score == 0.0
        assert result.weighted_total == 0.0

    def test_score_fields_accept_one(self):
        """All score fields should accept 1.0."""
        result = _make_result(cai=1.0, cpg=1.0, mrna_dg=1.0, weighted_total=1.8)
        assert result.cai_score == 1.0
        assert result.cpg_score == 1.0
        assert result.mrna_dg_score == 1.0

    def test_score_fields_accept_negative_values(self):
        """Score fields should accept negative values (dataclass has no validation)."""
        result = _make_result(cai=-0.5, cpg=-0.3, mrna_dg=-0.1, weighted_total=-0.9)
        assert result.cai_score == -0.5
        assert result.cpg_score == -0.3
        assert result.mrna_dg_score == -0.1
        assert result.weighted_total == -0.9

    def test_score_fields_accept_values_above_one(self):
        """Score fields should accept values > 1.0 (dataclass has no range check)."""
        result = _make_result(cai=2.0, cpg=3.0, mrna_dg=1.5, weighted_total=6.5)
        assert result.cai_score == 2.0
        assert result.cpg_score == 3.0

    def test_repr_contains_all_short_names(self):
        """ScoringResult.__repr__ should include cai, cpg, mrna_dg, total."""
        result = _make_result(cai=0.12, cpg=0.34, mrna_dg=0.56, weighted_total=0.78)
        r = repr(result)
        assert "cai=" in r
        assert "cpg=" in r
        assert "mrna_dg=" in r
        assert "total=" in r

    def test_repr_four_decimal_precision(self):
        """ScoringResult.__repr__ should format scores to 4 decimal places."""
        result = _make_result(cai=0.123456, cpg=0.654321, mrna_dg=0.111119, weighted_total=0.999999)
        r = repr(result)
        assert "0.1235" in r  # rounded from 0.123456
        assert "0.6543" in r  # rounded from 0.654321

    def test_repr_starts_with_class_name(self):
        """ScoringResult.__repr__ should start with 'ScoringResult('."""
        result = _make_result()
        assert repr(result).startswith("ScoringResult(")

    def test_score_fields_are_mutable(self):
        """ScoringResult fields should be mutable (it's not frozen)."""
        result = _make_result(cai=0.5)
        result.cai_score = 0.9
        assert result.cai_score == 0.9

    def test_individual_scores_dict_is_mutable(self):
        """The individual_scores dict should be mutable after construction."""
        result = _make_result(individual_scores={"A": 1.0})
        result.individual_scores["B"] = 2.0
        assert "B" in result.individual_scores
        assert result.individual_scores["B"] == 2.0

    def test_float_integer_coercion(self):
        """Integer values should be stored as-is (Python dataclass, no coercion)."""
        result = ScoringResult(
            cai_score=1, cpg_score=0, mrna_dg_score=1,
            weighted_total=2, individual_scores={},
        )
        # In Python, 1 == 1.0 but type(1) != float
        # Dataclass doesn't coerce, but the value should still compare equal
        assert result.cai_score == 1.0
        assert result.weighted_total == 2.0

    def test_construction_with_keyword_only(self):
        """ScoringResult should be constructable with keyword arguments."""
        result = ScoringResult(
            cai_score=0.5,
            cpg_score=0.6,
            mrna_dg_score=0.7,
            weighted_total=0.8,
            individual_scores={},
        )
        assert result.cai_score == 0.5

    def test_large_individual_scores_dict(self):
        """individual_scores should handle many entries without error."""
        scores = {f"constraint_{i}": float(i) for i in range(100)}
        result = _make_result(individual_scores=scores)
        assert len(result.individual_scores) == 100
        assert result.individual_scores["constraint_50"] == 50.0

    def test_nan_score_value(self):
        """ScoringResult should accept NaN values (no validation constraint)."""
        result = _make_result(cai=float("nan"))
        assert math.isnan(result.cai_score)

    def test_inf_score_value(self):
        """ScoringResult should accept infinity values (no validation constraint)."""
        result = _make_result(cpg=float("inf"))
        assert math.isinf(result.cpg_score)
        assert result.cpg_score > 0


# ════════════════════════════════════════════════════════════════════
# 4. Integration: scorer with real model produces valid ScoringResults
# ════════════════════════════════════════════════════════════════════

class TestScoringIntegration:
    """Integration tests combining scorer, model, and Pareto analysis."""

    def test_scored_results_are_valid_for_pareto(self):
        """ScoringResults from score_sequence() should work as input to
        compute_pareto_frontier without error."""
        config = SolverConfig()
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        r1 = scorer.score_sequence(model, "ATGAAG")
        r2 = scorer.score_sequence(model, "ATGAAA")
        pareto = compute_pareto_frontier([r1, r2])
        assert len(pareto) >= 1
        for p in pareto:
            assert isinstance(p, ScoringResult)
            assert 0.0 <= p.cai_score <= 1.0
            assert 0.0 <= p.cpg_score <= 1.0
            assert 0.0 <= p.mrna_dg_score <= 1.0

    def test_normalized_scores_bounded_after_scoring(self):
        """All normalized scores from score_sequence() should be in [0.0, 1.0]."""
        config = SolverConfig()
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert 0.0 <= result.cai_score <= 1.0
        assert 0.0 <= result.cpg_score <= 1.0
        assert 0.0 <= result.mrna_dg_score <= 1.0

    def test_individual_scores_keys_match_constraint_names(self):
        """individual_scores dict keys should match the constraint .name attributes."""
        config = SolverConfig()
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        result = scorer.score_sequence(model, "ATGAAG")
        assert "MaximizeCAI" in result.individual_scores
        assert "MinimizeCpG" in result.individual_scores
        assert "MinimizeMRNADG" in result.individual_scores

    def test_different_sequences_produce_different_results(self):
        """Two different codon assignments for the same protein should yield
        different ScoringResults (at least in CAI)."""
        config = SolverConfig()
        model = _make_model(config=config)
        scorer = SoftConstraintScorer(config)
        r1 = scorer.score_sequence(model, "ATGAAG")
        r2 = scorer.score_sequence(model, "ATGAAA")
        assert r1.cai_score != r2.cai_score
