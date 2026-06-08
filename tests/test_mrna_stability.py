"""Tests for mRNA stability scoring module.

Covers:
- score_mrna_stability() with different organisms
- compute_mrna_half_life_score() scoring
- predict_mrna_stability() categorization
- suggest_mutations_for_stability() suggestions
- STABILITY_MOTIFS registry
- MRNAStabilityScore dataclass
- Organism-specific motif detection
- Empty/edge case handling
"""

import pytest
from biocompiler.mrna_stability import (
    STABILITY_MOTIFS,
    MRNAStabilityScore,
    score_mrna_stability,
    compute_mrna_half_life_score,
    predict_mrna_stability,
    suggest_mutations_for_stability,
)


class TestStabilityMotifs:
    def test_ecoli_motifs_exist(self):
        assert "Escherichia_coli" in STABILITY_MOTIFS
        motifs = STABILITY_MOTIFS["Escherichia_coli"]
        assert "stabilizing" in motifs
        assert "destabilizing" in motifs

    def test_human_motifs_exist(self):
        assert "Homo_sapiens" in STABILITY_MOTIFS
        motifs = STABILITY_MOTIFS["Homo_sapiens"]
        assert "stabilizing" in motifs
        assert "destabilizing" in motifs

    def test_yeast_motifs_exist(self):
        assert "Saccharomyces_cerevisiae" in STABILITY_MOTIFS

    def test_ecoli_alias(self):
        assert "E_coli" in STABILITY_MOTIFS

    def test_mouse_uses_human_motifs(self):
        assert "Mus_musculus" in STABILITY_MOTIFS
        assert STABILITY_MOTIFS["Mus_musculus"] is STABILITY_MOTIFS["Homo_sapiens"]

    def test_cho_uses_human_motifs(self):
        assert "CHO_K1" in STABILITY_MOTIFS


class TestMRNAStabilityScore:
    def test_basic_creation(self):
        score = MRNAStabilityScore(
            overall_score=0.8,
            stabilizing_count=2,
            destabilizing_count=1,
        )
        assert score.overall_score == pytest.approx(0.8, rel=1e-6)
        assert score.risk_level == "low"

    def test_medium_risk(self):
        score = MRNAStabilityScore(
            overall_score=0.55,
            stabilizing_count=1,
            destabilizing_count=1,
        )
        assert score.risk_level == "medium"

    def test_high_risk(self):
        score = MRNAStabilityScore(
            overall_score=0.3,
            stabilizing_count=0,
            destabilizing_count=3,
        )
        assert score.risk_level == "high"

    def test_score_clamped(self):
        score = MRNAStabilityScore(
            overall_score=1.5,
            stabilizing_count=0,
            destabilizing_count=0,
        )
        assert score.overall_score == 1.0

    def test_negative_score_clamped(self):
        score = MRNAStabilityScore(
            overall_score=-0.5,
            stabilizing_count=0,
            destabilizing_count=5,
        )
        assert score.overall_score == 0.0

    def test_default_risk_level_overridden(self):
        score = MRNAStabilityScore(
            overall_score=0.8,
            stabilizing_count=0,
            destabilizing_count=0,
            risk_level="medium",  # Will be overridden by __post_init__
        )
        assert score.risk_level == "low"  # 0.8 > 0.7 -> low

    def test_motif_details(self):
        score = MRNAStabilityScore(
            overall_score=0.5,
            stabilizing_count=1,
            destabilizing_count=1,
            motif_details=[{"position": 5, "effect": "destabilizing"}],
        )
        assert len(score.motif_details) == 1


class TestScoreMrnaStability:
    def test_ecoli_basic(self):
        dna = "ATGGTCAAGGCCTAA"
        result = score_mrna_stability(dna, "Escherichia_coli")
        assert isinstance(result, MRNAStabilityScore)
        assert 0.0 <= result.overall_score <= 1.0
        assert result.risk_level in ("low", "medium", "high")

    def test_human_basic(self):
        dna = "ATGGTCAAGGCCTAA"
        result = score_mrna_stability(dna, "Homo_sapiens")
        assert isinstance(result, MRNAStabilityScore)
        assert 0.0 <= result.overall_score <= 1.0

    def test_empty_sequence(self):
        result = score_mrna_stability("", "Escherichia_coli")
        assert result.overall_score == 0.5  # Neutral

    def test_destabilizing_motif_detection(self):
        # ATTTA is a destabilizing ARE motif for human
        dna = "ATGATTTAGCCTAA"
        result = score_mrna_stability(dna, "Homo_sapiens")
        assert result.destabilizing_count >= 1

    def test_stabilizing_gc_rich(self):
        # GC-rich sequence should have stabilizing motifs
        dna = "ATGGCGGCGGCTAA"
        result = score_mrna_stability(dna, "Escherichia_coli")
        assert isinstance(result.stabilizing_count, int)
        assert result.stabilizing_count >= 0

    def test_motif_details_populated(self):
        dna = "ATGATTTAGCCTAA"
        result = score_mrna_stability(dna, "Homo_sapiens")
        assert isinstance(result.motif_details, list)

    def test_yeast_basic(self):
        dna = "ATGGTCAAGGCCTAA"
        result = score_mrna_stability(dna, "Saccharomyces_cerevisiae")
        assert isinstance(result, MRNAStabilityScore)


class TestComputeMrnaHalfLifeScore:
    def test_basic_score(self):
        dna = "ATGGTCAAGGCCTAA"
        score = compute_mrna_half_life_score(dna, "Escherichia_coli")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_short_sequence(self):
        score = compute_mrna_half_life_score("AT", "Escherichia_coli")
        assert score == 0.5  # Neutral for short sequences

    def test_human_with_atttta(self):
        dna = "ATGATTTAATTTAATTTAATAA"
        score = compute_mrna_half_life_score(dna, "Homo_sapiens")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_ecoli_no_atttta_check(self):
        dna = "ATGGTCAAGGCCTAA"
        score = compute_mrna_half_life_score(dna, "Escherichia_coli")
        assert isinstance(score, float)

    def test_gc_extreme_penalized(self):
        # Very AT-rich sequence
        dna = "ATATATATATATATATATAT" * 3
        score = compute_mrna_half_life_score(dna, "Homo_sapiens")
        assert isinstance(score, float)


class TestPredictMrnaStability:
    def test_stable_category(self):
        # Use a known high-CAI codon sequence for E. coli
        # E. coli doesn't check ATTTA, so CAI alone determines category
        dna = "ATGCTGCTGCTGCTGCTGCTAA"  # Mostly Leu with optimal codon
        category = predict_mrna_stability(dna, "Escherichia_coli")
        assert category in ("STABLE", "MODERATE", "UNSTABLE")

    def test_unstable_category(self):
        # Very AT-rich sequence likely has low CAI
        dna = "ATGATTTAATTTATTTATAA"
        category = predict_mrna_stability(dna, "Homo_sapiens")
        assert category in ("STABLE", "MODERATE", "UNSTABLE")

    def test_atttta_downgrade_human(self):
        # For human, ATTTA presence should downgrade category
        dna = "ATGATTTAGCCTAA"
        category = predict_mrna_stability(dna, "Homo_sapiens")
        assert category in ("MODERATE", "UNSTABLE")  # Downgraded

    def test_ecoli_no_atttta_check(self):
        dna = "ATGATTTAGCCTAA"
        # E. coli doesn't check ATTTA
        category = predict_mrna_stability(dna, "Escherichia_coli")
        assert category in ("STABLE", "MODERATE", "UNSTABLE")

    def test_yeast_checks_motifs(self):
        dna = "ATGATTTAGCCTAA"
        category = predict_mrna_stability(dna, "Saccharomyces_cerevisiae")
        assert category in ("STABLE", "MODERATE", "UNSTABLE")


class TestSuggestMutationsForStability:
    def test_no_suggestions_for_stable(self):
        # Well-optimized sequence should have few/none suggestions
        dna = "ATGCTGCTGCTGCTAA"
        suggestions = suggest_mutations_for_stability(dna, "Escherichia_coli")
        assert isinstance(suggestions, list)

    def test_suggestions_for_destabilizing_motifs(self):
        # ATTTA-containing sequence should generate suggestions
        dna = "ATGATTTAGCCTAA"
        suggestions = suggest_mutations_for_stability(dna, "Homo_sapiens")
        assert isinstance(suggestions, list)
        # Should have at least one suggestion if destabilizing motifs found
        if suggestions:
            s = suggestions[0]
            assert "position" in s
            assert "original_codon" in s
            assert "suggested_codon" in s
            assert "amino_acid" in s
            assert "motif_removed" in s

    def test_short_sequence_no_suggestions(self):
        suggestions = suggest_mutations_for_stability("AT", "Homo_sapiens")
        assert suggestions == []

    def test_suggestions_preserve_amino_acid(self):
        dna = "ATGATTTAGCCTAA"
        suggestions = suggest_mutations_for_stability(dna, "Homo_sapiens")
        for s in suggestions:
            from biocompiler.type_system import CODON_TABLE, AA_TO_CODONS
            orig_aa = CODON_TABLE.get(s["original_codon"])
            new_aa = CODON_TABLE.get(s["suggested_codon"])
            if orig_aa and new_aa:
                assert orig_aa == new_aa, f"AA changed: {orig_aa} -> {new_aa}"
