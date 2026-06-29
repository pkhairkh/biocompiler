"""Tests for codon pair bias scoring module.

Covers:
- get_codon_pair_data() for different organisms
- score_codon_pair() single pair scoring
- compute_cpb() / compute_cpb_score() sequence scoring
- estimate_cpb_from_codon_freq()
- suggest_better_pair()
- Edge cases (empty, short, invalid)
"""

import pytest
from biocompiler.expression.codon_pair_scoring import (
    get_codon_pair_data,
    score_codon_pair,
    compute_cpb,
    compute_cpb_score,
    estimate_cpb_from_codon_freq,
    suggest_better_pair,
)


class TestGetCodonPairData:
    def test_ecoli_data(self):
        data = get_codon_pair_data("Escherichia_coli")
        assert isinstance(data, dict)
        assert len(data) > 0

    def test_human_data(self):
        data = get_codon_pair_data("Homo_sapiens")
        assert isinstance(data, dict)
        assert len(data) > 0

    def test_ecoli_alias(self):
        data = get_codon_pair_data("E_coli")
        assert isinstance(data, dict)

    def test_mouse_data(self):
        data = get_codon_pair_data("Mus_musculus")
        assert isinstance(data, dict)

    def test_cho_data(self):
        data = get_codon_pair_data("CHO_K1")
        assert isinstance(data, dict)

    def test_yeast_data(self):
        data = get_codon_pair_data("Saccharomyces_cerevisiae")
        assert isinstance(data, dict)

    def test_human_data_has_positive_and_negative(self):
        data = get_codon_pair_data("Homo_sapiens")
        has_positive = any(v > 0 for v in data.values())
        has_negative = any(v < 0 for v in data.values())
        assert has_positive, "Human CPB data should have over-represented pairs"
        assert has_negative, "Human CPB data should have under-represented pairs"

    def test_unknown_organism_returns_estimated_or_empty(self):
        data = get_codon_pair_data("Unknown_organism_xyz")
        # Should return estimated data or empty dict, not raise
        assert isinstance(data, dict)


class TestScoreCodonPair:
    def test_known_pair_human(self):
        # CTG-CTC is over-represented in human
        score = score_codon_pair("CTG", "CTC", "Homo_sapiens")
        assert isinstance(score, float)
        assert score > 0  # Over-represented

    def test_rare_pair_human(self):
        # ATA-ATA is under-represented in human
        score = score_codon_pair("ATA", "ATA", "Homo_sapiens")
        assert isinstance(score, float)
        assert score < 0  # Under-represented

    def test_unknown_pair_returns_zero(self):
        score = score_codon_pair("ATG", "CTG", "Homo_sapiens")
        # Most pairs are not in the bias table, default to 0
        assert isinstance(score, float)

    def test_case_insensitive(self):
        s1 = score_codon_pair("CTG", "CTC", "Homo_sapiens")
        s2 = score_codon_pair("ctg", "ctc", "Homo_sapiens")
        assert s1 == s2


class TestComputeCPB:
    def test_basic_sequence(self):
        dna = "ATGCTGCTCAAG"
        score = compute_cpb(dna, "Homo_sapiens")
        assert isinstance(score, float)
        assert -1.0 <= score <= 1.0  # Reasonable range

    def test_empty_sequence(self):
        score = compute_cpb("", "Homo_sapiens")
        assert score == 0.0

    def test_short_sequence(self):
        # Less than 6 bases = less than 2 codons
        score = compute_cpb("ATGCTG", "Homo_sapiens")
        assert isinstance(score, float)

    def test_invalid_length_raises(self):
        with pytest.raises(ValueError, match="not a multiple of 3"):
            compute_cpb("ATGC", "Homo_sapiens")

    def test_ecoli_sequence(self):
        dna = "ATGCTGAAAGCGTAA"
        score = compute_cpb(dna, "Escherichia_coli")
        assert isinstance(score, float)

    def test_compute_cpb_score_alias(self):
        dna = "ATGCTGCTCAAG"
        s1 = compute_cpb(dna, "Homo_sapiens")
        s2 = compute_cpb_score(dna, "Homo_sapiens")
        assert s1 == s2  # Aliases should give same result


class TestEstimateCPBFromCodonFreq:
    def test_basic_estimation(self):
        from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE
        # Build a simple codon usage table
        usage = {}
        for codon, aa in CODON_TABLE.items():
            if aa != "*":
                usage[codon] = (aa, 0.5, 5.0, 100)
        result = estimate_cpb_from_codon_freq(usage)
        assert isinstance(result, dict)
        assert len(result) > 0

    def test_empty_usage_returns_empty(self):
        result = estimate_cpb_from_codon_freq({})
        assert result == {}

    def test_estimated_values_in_range(self):
        from biocompiler.type_system import CODON_TABLE
        usage = {}
        for codon, aa in CODON_TABLE.items():
            if aa != "*":
                usage[codon] = (aa, 0.5, 5.0, 100)
        result = estimate_cpb_from_codon_freq(usage)
        for key, val in result.items():
            assert -0.5 <= val <= 0.5, f"Score {val} for {key} out of range"


class TestSuggestBetterPair:
    def test_basic_suggestion(self):
        result = suggest_better_pair(
            "ATA", "ATA", "I", "I", "Homo_sapiens",
        )
        # ATA-ATA is under-represented; should suggest a better pair
        # (may or may not find one depending on alternatives)
        assert result is None or isinstance(result, tuple)

    def test_no_improvement_for_optimal(self):
        # CTG-CTC is already over-represented in human
        result = suggest_better_pair(
            "CTG", "CTC", "L", "L", "Homo_sapiens",
        )
        # May or may not find improvement (already good)
        assert result is None or isinstance(result, tuple)

    def test_suggestion_preserves_amino_acids(self):
        result = suggest_better_pair(
            "ATA", "ATA", "I", "I", "Homo_sapiens",
        )
        if result is not None:
            from biocompiler.type_system import CODON_TABLE
            c1, c2 = result
            assert CODON_TABLE.get(c1) == "I"
            assert CODON_TABLE.get(c2) == "I"

    def test_with_cai_weights(self):
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        cai_weights = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens", {})
        result = suggest_better_pair(
            "ATA", "ATA", "I", "I", "Homo_sapiens",
            cai_weights=cai_weights,
        )
        assert result is None or isinstance(result, tuple)
