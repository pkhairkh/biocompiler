"""Tests for UTR models module.

Covers:
- UTRConfig dataclass
- ORGANISM_UTR_CONFIGS registry
- score_5utr() for different organisms
- score_3utr() for different organisms
- suggest_5utr() for different organisms
- suggest_3utr() for different organisms
- Alias resolution
- Edge cases
"""

import pytest
from biocompiler.utr_models import (
    UTRConfig,
    ORGANISM_UTR_CONFIGS,
    AVAILABLE_ORGANISMS,
    score_5utr,
    score_3utr,
    suggest_5utr,
    suggest_3utr,
)


class TestUTRConfig:
    def test_ecoli_config_exists(self):
        assert "Escherichia_coli" in ORGANISM_UTR_CONFIGS
        config = ORGANISM_UTR_CONFIGS["Escherichia_coli"]
        assert config.shine_dalgarno is not None
        assert config.kozak_sequence is None
        assert config.splicing_signals is False

    def test_human_config_exists(self):
        assert "Homo_sapiens" in ORGANISM_UTR_CONFIGS
        config = ORGANISM_UTR_CONFIGS["Homo_sapiens"]
        assert config.kozak_sequence is not None
        assert config.shine_dalgarno is None
        assert config.splicing_signals is True

    def test_yeast_config_exists(self):
        assert "Saccharomyces_cerevisiae" in ORGANISM_UTR_CONFIGS
        config = ORGANISM_UTR_CONFIGS["Saccharomyces_cerevisiae"]
        assert config.polya_signal is not None

    def test_mouse_config_exists(self):
        assert "Mus_musculus" in ORGANISM_UTR_CONFIGS

    def test_cho_config_exists(self):
        assert "CHO_K1" in ORGANISM_UTR_CONFIGS

    def test_available_organisms_list(self):
        assert isinstance(AVAILABLE_ORGANISMS, list)
        assert len(AVAILABLE_ORGANISMS) >= 5
        assert "Escherichia_coli" in AVAILABLE_ORGANISMS


class TestScore5utr:
    def test_ecoli_good_utr(self):
        # Good Shine-Dalgarno + proper spacing
        score = score_5utr("TAAGGAGGTAAAAAAAATG", "Escherichia_coli")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0
        assert score > 0.5  # Should score well

    def test_ecoli_no_sd(self):
        # No Shine-Dalgarno sequence
        score = score_5utr("AAAAAAAAAATG", "Escherichia_coli")
        assert isinstance(score, float)
        assert score < 0.8  # Should score poorly

    def test_human_kozak_perfect(self):
        # Perfect Kozak consensus
        score = score_5utr("GCCACCATGG", "Homo_sapiens")
        assert isinstance(score, float)
        assert score > 0.5

    def test_human_no_kozak(self):
        # Poor Kozak context
        score = score_5utr("TTTTTTATG", "Homo_sapiens")
        assert isinstance(score, float)

    def test_yeast_utr(self):
        score = score_5utr("TATATAAAAAAAAACAATG", "Saccharomyces_cerevisiae")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_empty_sequence(self):
        score = score_5utr("", "Homo_sapiens")
        assert score == 0.0

    def test_alias_ecoli(self):
        score = score_5utr("TAAGGAGGTAAAAAAAATG", "E_coli")
        assert isinstance(score, float)

    def test_alias_human(self):
        score = score_5utr("GCCACCATGG", "human")
        assert isinstance(score, float)

    def test_unknown_organism_raises(self):
        with pytest.raises(ValueError, match="Unknown organism"):
            score_5utr("ATGATG", "Alien_martian")


class TestScore3utr:
    def test_ecoli_terminator(self):
        # Good Rho-independent terminator
        score = score_3utr("GCGCCGCTTTTTT", "Escherichia_coli")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_human_polya_signal(self):
        # Good polyA signal
        score = score_3utr("AAAAAAAAAAAAAAATAAAAATAAA", "Homo_sapiens")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_yeast_3utr(self):
        score = score_3utr("TAATAAATAA", "Saccharomyces_cerevisiae")
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_empty_sequence(self):
        score = score_3utr("", "Homo_sapiens")
        assert score == 0.0

    def test_alias_ecoli(self):
        score = score_3utr("GCGCCGCTTTTTT", "E_coli")
        assert isinstance(score, float)


class TestSuggest5utr:
    def test_ecoli_suggestion(self):
        utr = suggest_5utr("Escherichia_coli")
        assert isinstance(utr, str)
        assert len(utr) > 0
        assert "ATG" in utr  # Should contain start codon

    def test_human_suggestion(self):
        utr = suggest_5utr("Homo_sapiens")
        assert isinstance(utr, str)
        assert len(utr) > 0
        assert "ATG" in utr or "GCCACC" in utr

    def test_yeast_suggestion(self):
        utr = suggest_5utr("Saccharomyces_cerevisiae")
        assert isinstance(utr, str)
        assert len(utr) > 0

    def test_ecoli_suggestion_has_sd(self):
        utr = suggest_5utr("Escherichia_coli")
        assert "AGGAGG" in utr or "GGAG" in utr  # SD-like sequence

    def test_human_suggestion_has_kozak(self):
        utr = suggest_5utr("Homo_sapiens")
        assert "GCCACC" in utr  # Kozak consensus part

    def test_alias(self):
        utr = suggest_5utr("ecoli")
        assert isinstance(utr, str)

    def test_unknown_organism_raises(self):
        with pytest.raises(ValueError, match="Unknown organism"):
            suggest_5utr("Alien_martian")


class TestSuggest3utr:
    def test_ecoli_suggestion(self):
        utr = suggest_3utr("Escherichia_coli")
        assert isinstance(utr, str)
        assert len(utr) > 0

    def test_human_suggestion(self):
        utr = suggest_3utr("Homo_sapiens")
        assert isinstance(utr, str)
        assert len(utr) > 0
        assert "AATAAA" in utr  # PolyA signal

    def test_yeast_suggestion(self):
        utr = suggest_3utr("Saccharomyces_cerevisiae")
        assert isinstance(utr, str)
        assert len(utr) > 0

    def test_ecoli_has_terminator(self):
        utr = suggest_3utr("Escherichia_coli")
        # Should have GC-rich region and poly-T tract
        assert "GCG" in utr or "TTT" in utr

    def test_unknown_organism_raises(self):
        with pytest.raises(ValueError, match="Unknown organism"):
            suggest_3utr("Alien_martian")
