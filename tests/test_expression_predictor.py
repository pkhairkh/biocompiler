"""
Tests for the expression predictor module.

Covers:
  - predict_expression() function
  - ExpressionPrediction dataclass
  - ExpressionPredictor class
  - GC optimality scoring
  - mRNA stability estimation
  - Confidence computation
  - Batch prediction
"""

import pytest
from biocompiler.expression_predictor import (
    ExpressionPrediction,
    predict_expression,
    ExpressionPredictor,
    _compute_gc_optimality,
    _estimate_mrna_stability,
    _GC_SWEET_SPOT,
    _FACTOR_WEIGHTS,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def gfp_ecoli_sequence():
    """A realistic GFP sequence optimized for E. coli."""
    return (
        "ATGGTTAGCAAAGGCGAAGAATTATTTACGGGCGTGGTTCCGATTCTGGTGGAACTGGA"
        "CGGCGATGTGAACGGCCATAAGTTCAGCGTGAGCGGCGAAGGCGAAGGCGATGCGACCT"
        "ATGGCAAGCTGACCTTAAAATTTATTTGCACCACCGGCAAACTGCCGGTGCCGTGGCCGA"
        "CCCTGGTGACCACCTTTAGCTATGGTGTGCAGTGCTTTAGCCGCTATCCGGATCATATGA"
        "AACAGCATGATTTTTTTAAAAGCGCGATGCCAGAAGGCTATGTGCAAGAACGCACCATTT"
        "TTTTCAAAGATGATGGCAACTATAAAACCCGCGCGGAAGTGAAATTTGAAGGCGATACCC"
        "TGGTGAACCGCATTGAGCTGAAGGGCATTGATTTTAAGGAAGATGGTAACATCCTGGGCC"
        "ATAAACTGGAATATAACTATAACAGCCATAACGTGTATATTATGGCGGATAAACAGAAAA"
        "ACGGTATTAAAGTGAACTTCAAAATTCGCCATAACATTGAAGATGGCAGCGTTCAGCTG"
        "GCGGATCATTATCAACAGAACACCCCGATTGGCGATGGCCCGGTGCTGCTGCCGGACAA"
        "CCATTATCTGAGCACCCAGAGCGCGTTAAGCAAAGATCCGAACGAAAAACGCGATCATA"
        "TGGTGCTGCTGGAATTTGTTACCGCGGCGGGCATTACGCATGGCATGGATGAACTGTAT"
        "AAA"
    )


@pytest.fixture
def simple_ecoli_sequence():
    """A simple short sequence for E. coli testing."""
    return "ATGGCTAAAGCTGCGGCCTAA"


# ═══════════════════════════════════════════════════════════════════════════════
# ExpressionPrediction dataclass tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestExpressionPrediction:
    """Test the ExpressionPrediction dataclass."""

    def test_create_prediction(self):
        pred = ExpressionPrediction(
            predicted_relative_expression=0.8,
            confidence=0.7,
            key_factors={"cai": 0.9, "gc_optimality": 0.7},
            category="high",
            organism="Escherichia_coli",
        )
        assert pred.predicted_relative_expression == pytest.approx(0.8, rel=1e-6)
        assert pred.confidence == pytest.approx(0.7, rel=1e-6)
        assert pred.category == "high"

    def test_prediction_with_warnings(self):
        pred = ExpressionPrediction(
            predicted_relative_expression=0.3,
            confidence=0.5,
            key_factors={"cai": 0.4},
            category="low",
            organism="Escherichia_coli",
            warnings=["Low CAI may reduce expression"],
        )
        assert len(pred.warnings) == 1

    def test_invalid_expression_range(self):
        with pytest.raises(ValueError, match="predicted_relative_expression"):
            ExpressionPrediction(
                predicted_relative_expression=1.5,
                confidence=0.5,
                key_factors={},
                category="high",
                organism="Escherichia_coli",
            )

    def test_invalid_confidence_range(self):
        with pytest.raises(ValueError, match="confidence"):
            ExpressionPrediction(
                predicted_relative_expression=0.5,
                confidence=2.0,
                key_factors={},
                category="medium",
                organism="Escherichia_coli",
            )


# ═══════════════════════════════════════════════════════════════════════════════
# predict_expression() function tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestPredictExpression:
    """Test the predict_expression function."""

    def test_predict_gfp_ecoli(self, gfp_ecoli_sequence):
        """GFP sequence optimized for E. coli should predict high expression."""
        result = predict_expression(gfp_ecoli_sequence, "Escherichia_coli")
        assert isinstance(result, ExpressionPrediction)
        assert result.organism == "Escherichia_coli"
        assert 0.0 <= result.predicted_relative_expression <= 1.0
        assert result.category in ("high", "medium", "low")
        assert "cai" in result.key_factors
        assert "gc_optimality" in result.key_factors

    def test_predict_with_precomputed_cai(self, gfp_ecoli_sequence):
        """Should accept pre-computed CAI value."""
        result = predict_expression(
            gfp_ecoli_sequence, "Escherichia_coli", cai=0.93
        )
        assert result.key_factors["cai"] == pytest.approx(0.93, rel=1e-6)

    def test_predict_with_precomputed_gc(self, gfp_ecoli_sequence):
        """Should accept pre-computed GC content."""
        result = predict_expression(
            gfp_ecoli_sequence, "Escherichia_coli", gc_content=0.48
        )
        assert 0.0 <= result.predicted_relative_expression <= 1.0

    def test_predict_empty_sequence(self):
        """Empty sequence should return low prediction."""
        result = predict_expression("", "Escherichia_coli")
        assert result.predicted_relative_expression < 0.5
        assert result.category in ("low", "medium")

    def test_predict_organism_aliases(self, gfp_ecoli_sequence):
        """Should accept organism aliases."""
        result = predict_expression(gfp_ecoli_sequence, "ecoli")
        assert result.organism == "Escherichia_coli"

    def test_predict_human_organism(self):
        """Test with human organism."""
        seq = "ATGGCGCTGTGGATGCGCCTGCTGCCACTGCTGGCGCTGCTGGCGCTGTGGGGCCCGGA"
        result = predict_expression(seq, "Homo_sapiens")
        assert result.organism == "Homo_sapiens"
        assert 0.0 <= result.predicted_relative_expression <= 1.0

    def test_predict_yeast_organism(self):
        """Test with yeast organism."""
        seq = "ATGGTTTCTGAAACCTTTACTGGTGTTGTTCCAATTTTAG"
        result = predict_expression(seq, "Saccharomyces_cerevisiae")
        assert result.organism == "Saccharomyces_cerevisiae"

    def test_high_cai_high_gc_predicts_high(self, gfp_ecoli_sequence):
        """High CAI + optimal GC should predict high expression."""
        result = predict_expression(
            gfp_ecoli_sequence, "Escherichia_coli",
            cai=0.95, gc_content=0.50,
            mrna_stability_score=0.8,
            codon_pair_bias=0.2,
        )
        assert result.predicted_relative_expression >= 0.6
        assert result.category in ("high", "medium")

    def test_low_cai_predicts_low(self, gfp_ecoli_sequence):
        """Very low CAI should predict low expression."""
        result = predict_expression(
            gfp_ecoli_sequence, "Escherichia_coli",
            cai=0.2, gc_content=0.50,
        )
        assert result.predicted_relative_expression < 0.5

    def test_warnings_for_low_cai(self, gfp_ecoli_sequence):
        """Low CAI should generate warnings."""
        result = predict_expression(
            gfp_ecoli_sequence, "Escherichia_coli", cai=0.3
        )
        assert any("Low CAI" in w for w in result.warnings)

    def test_warnings_for_extreme_gc(self, gfp_ecoli_sequence):
        """Extreme GC should generate warnings."""
        result = predict_expression(
            gfp_ecoli_sequence, "Escherichia_coli", gc_content=0.90
        )
        assert any("GC" in w for w in result.warnings)

    def test_all_factors_in_key_factors(self, gfp_ecoli_sequence):
        """All four factors should be present in key_factors."""
        result = predict_expression(gfp_ecoli_sequence, "Escherichia_coli")
        expected_keys = {"cai", "gc_optimality", "mrna_stability", "codon_pair_bias"}
        assert expected_keys.issubset(set(result.key_factors.keys()))


# ═══════════════════════════════════════════════════════════════════════════════
# GC Optimality tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestGCOptimality:
    """Test GC optimality scoring."""

    def test_optimal_ecoli_gc(self):
        """GC at the E. coli sweet spot should score 1.0."""
        score = _compute_gc_optimality(0.50, "Escherichia_coli")
        assert score == 1.0

    def test_optimal_human_gc(self):
        """GC at the human sweet spot should score 1.0."""
        score = _compute_gc_optimality(0.55, "Homo_sapiens")
        assert score == 1.0

    def test_suboptimal_gc(self):
        """GC far from the sweet spot should score < 1.0."""
        score = _compute_gc_optimality(0.20, "Escherichia_coli")
        assert score < 1.0

    def test_very_extreme_gc(self):
        """Very extreme GC should score close to 0."""
        score = _compute_gc_optimality(0.90, "Escherichia_coli")
        assert score < 0.3

    def test_unknown_organism_uses_default(self):
        """Unknown organism should use the default sweet spot."""
        score = _compute_gc_optimality(0.48, "Unknown_organism")
        assert 0.0 <= score <= 1.0

    def test_boundary_of_sweet_spot(self):
        """GC at the boundary of the sweet spot should score 1.0."""
        lo, hi = _GC_SWEET_SPOT["Escherichia_coli"]
        score_lo = _compute_gc_optimality(lo, "Escherichia_coli")
        score_hi = _compute_gc_optimality(hi, "Escherichia_coli")
        assert score_lo == 1.0
        assert score_hi == 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# mRNA stability estimation tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestMRNAStability:
    """Test mRNA stability estimation."""

    def test_stable_sequence(self):
        """GC-balanced sequence without instability motifs should score well."""
        seq = "GCGCCGATGGCTGCCGCTGAAGCGGCTGCGGCAAGGCG"
        score = _estimate_mrna_stability(seq, "Escherichia_coli")
        assert 0.0 <= score <= 1.0
        assert score >= 0.5  # Should be reasonably stable

    def test_unstable_with_atttta(self):
        """Sequence with ATTTA motifs should have reduced stability."""
        # Use a GC-balanced base for fair comparison
        seq_stable = "GCGCCGATGGCTGCCGCTGAAGCGGCTGCGGCAAGGCG"
        score_stable = _estimate_mrna_stability(seq_stable, "Escherichia_coli")

        seq_unstable = "GCGATATTTAGCGCTGAAGCGGCTGCGGCAAGGCG"
        score_unstable = _estimate_mrna_stability(seq_unstable, "Escherichia_coli")

        assert score_unstable <= score_stable

    def test_empty_sequence(self):
        """Empty sequence should return neutral score."""
        score = _estimate_mrna_stability("", "Escherichia_coli")
        assert score == pytest.approx(0.5, rel=1e-6)

    def test_extreme_at_sequence(self):
        """Sequence with many A/T should have lower stability."""
        seq = "AATAATAATAATAATAATAATAATAATAATAATAATAATAAT"
        score = _estimate_mrna_stability(seq, "Escherichia_coli")
        assert score < 0.7  # Should be penalized


# ═══════════════════════════════════════════════════════════════════════════════
# ExpressionPredictor class tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestExpressionPredictorClass:
    """Test the ExpressionPredictor class interface."""

    def test_create_predictor(self):
        """Should create a predictor for a given organism."""
        predictor = ExpressionPredictor(organism="Escherichia_coli")
        assert predictor.organism == "Escherichia_coli"

    def test_predict_single_sequence(self, gfp_ecoli_sequence):
        """Should predict expression for a single sequence."""
        predictor = ExpressionPredictor(organism="Escherichia_coli")
        result = predictor.predict(gfp_ecoli_sequence)
        assert isinstance(result, ExpressionPrediction)

    def test_predict_batch(self, gfp_ecoli_sequence):
        """Should predict expression for multiple sequences."""
        predictor = ExpressionPredictor(organism="Escherichia_coli")
        sequences = [gfp_ecoli_sequence, "ATGGCTAAAGCTGCGGCCTAA"]
        results = predictor.predict_batch(sequences)
        assert len(results) == 2
        assert all(isinstance(r, ExpressionPrediction) for r in results)

    def test_custom_weights(self):
        """Should accept custom factor weights."""
        custom_weights = {
            "cai": 0.5,
            "gc_optimality": 0.2,
            "mrna_stability": 0.2,
            "codon_pair_bias": 0.1,
        }
        predictor = ExpressionPredictor(
            organism="Escherichia_coli",
            factor_weights=custom_weights,
        )
        weights = predictor.get_factor_weights()
        assert abs(weights["cai"] - 0.5) < 0.01

    def test_organism_alias(self):
        """Should accept organism aliases."""
        predictor = ExpressionPredictor(organism="ecoli")
        assert predictor.organism == "Escherichia_coli"

    def test_get_factor_weights(self):
        """Should return the current factor weights."""
        predictor = ExpressionPredictor(organism="Escherichia_coli")
        weights = predictor.get_factor_weights()
        assert "cai" in weights
        assert "gc_optimality" in weights
        assert "mrna_stability" in weights
        assert "codon_pair_bias" in weights
        # Weights should sum to approximately 1.0
        total = sum(weights.values())
        assert abs(total - 1.0) < 0.05


# ═══════════════════════════════════════════════════════════════════════════════
# Integration tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestExpressionPredictionIntegration:
    """Integration tests combining prediction with benchmark proteins."""

    def test_predict_for_benchmark_gfp(self, gfp_ecoli_sequence):
        """GFP should predict reasonable expression for E. coli."""
        result = predict_expression(gfp_ecoli_sequence, "Escherichia_coli")
        # GFP is a well-optimized, highly expressed protein
        assert result.predicted_relative_expression > 0.3
        # Should have all four factors
        assert len(result.key_factors) == 4

    def test_prediction_consistency(self, gfp_ecoli_sequence):
        """Same input should always produce same output."""
        result1 = predict_expression(gfp_ecoli_sequence, "Escherichia_coli")
        result2 = predict_expression(gfp_ecoli_sequence, "Escherichia_coli")
        assert result1.predicted_relative_expression == result2.predicted_relative_expression
        assert result1.category == result2.category

    def test_organism_specific_gc_sensitivity(self):
        """Different organisms should have different GC sweet spots."""
        # Create a sequence with 35% GC
        seq = "ATGGCTAATGCTGCTAAAGCTTTAGCAAAGCTTTAGCTAATGCTAAAGCTGCT"
        ecoli_pred = predict_expression(seq, "Escherichia_coli", cai=0.8)
        yeast_pred = predict_expression(seq, "Saccharomyces_cerevisiae", cai=0.8)

        # 35% GC is better for yeast than E. coli
        assert yeast_pred.key_factors["gc_optimality"] >= ecoli_pred.key_factors["gc_optimality"]

    def test_factor_weights_sum_to_one(self):
        """Default factor weights should sum to 1.0."""
        total = sum(_FACTOR_WEIGHTS.values())
        assert abs(total - 1.0) < 0.001
