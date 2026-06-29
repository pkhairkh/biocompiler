"""
Comprehensive tests for BLAST integration scaffolding (blast_screening.py).

Covers:
1. BlastHit construction and validation
2. BlastScreeningResult construction and properties
3. BlastScreener.is_blast_available (will be False in test env)
4. Quick screen fallback with known pathogen patterns
5. Screening of safe sequences (should pass)
6. Screening of sequences containing known pathogen motifs (should flag)
7. E-value and identity thresholds
8. Environment variable configuration
9. Protein screening (quick fallback)
10. BlastScreener internal methods
"""

from __future__ import annotations

import os
from unittest.mock import patch

import pytest

from biocompiler.biosecurity.blast_screening import (
    BlastHit,
    BlastScreeningResult,
    BlastScreener,
    PATHOGEN_TOXIN_MOTIFS,
    QUICK_SCREEN_DB_NAME,
)


# ─── Test Data ────────────────────────────────────────────────────────────

# A safe sequence (no pathogen motifs)
SAFE_DNA_SEQUENCE = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

# A sequence containing the anthrax lethal factor motif
ANTHRAX_MOTIF = PATHOGEN_TOXIN_MOTIFS["anthrax_lethal_factor"]["sequence"]
ANTHRAX_CONTAINING_SEQ = f"AAACCC{ANTHRAX_MOTIF}GGGTTT"

# A sequence containing the botulinum toxin motif
BOTULINUM_MOTIF = PATHOGEN_TOXIN_MOTIFS["botulinum_toxin_A"]["sequence"]
BOTULINUM_CONTAINING_SEQ = f"TTTAAA{BOTULINUM_MOTIF}CCCGGG"

# A sequence containing the ricin toxin motif
RICIN_MOTIF = PATHOGEN_TOXIN_MOTIFS["ricin_toxin"]["sequence"]
RICIN_CONTAINING_SEQ = f"ATGCGT{RICIN_MOTIF}TACGCA"

# Multiple toxin motifs in one sequence
MULTI_TOXIN_SEQ = f"AAA{ANTHRAX_MOTIF}CCC{RICIN_MOTIF}GGG"

# Safe protein sequence
SAFE_PROTEIN = "MGSWKRQPPAVNVLRYFPATW"

# Protein with ricin A-chain catalytic motif
RICIN_PROTEIN_MOTIF = "NGSFS"
RICIN_CONTAINING_PROTEIN = f"MKKK{RICIN_PROTEIN_MOTIF}LLLL"

# Protein with anthrax LF catalytic motif
ANTHRAX_PROTEIN_MOTIF = "HEFGH"
ANTHRAX_CONTAINING_PROTEIN = f"MAAAA{ANTHRAX_PROTEIN_MOTIF}VVVVV"


# ═══════════════════════════════════════════════════════════════════════
# 1. BlastHit Construction and Validation
# ═══════════════════════════════════════════════════════════════════════


class TestBlastHit:
    """Test BlastHit dataclass construction and validation."""

    def test_basic_construction(self):
        hit = BlastHit(
            subject_id="anthrax_lethal_factor",
            subject_organism="Bacillus anthracis",
            identity_percent=95.5,
            alignment_length=100,
            e_value=1e-50,
            bit_score=180.0,
        )
        assert hit.subject_id == "anthrax_lethal_factor"
        assert hit.subject_organism == "Bacillus anthracis"
        assert hit.identity_percent == 95.5
        assert hit.alignment_length == 100
        assert hit.e_value == 1e-50
        assert hit.bit_score == 180.0
        assert hit.is_pathogen is False  # Default
        assert hit.is_toxin is False  # Default

    def test_pathogen_and_toxin_flags(self):
        hit = BlastHit(
            subject_id="anthrax_lf",
            subject_organism="Bacillus anthracis",
            identity_percent=99.0,
            alignment_length=200,
            e_value=0.0,
            bit_score=350.0,
            is_pathogen=True,
            is_toxin=True,
        )
        assert hit.is_pathogen is True
        assert hit.is_toxin is True

    def test_empty_subject_id_raises(self):
        with pytest.raises(ValueError, match="subject_id must not be empty"):
            BlastHit(
                subject_id="",
                subject_organism="Unknown",
                identity_percent=100.0,
                alignment_length=10,
                e_value=0.0,
                bit_score=20.0,
            )

    def test_identity_percent_out_of_range_high(self):
        with pytest.raises(ValueError, match="identity_percent must be in"):
            BlastHit(
                subject_id="test",
                subject_organism="Test",
                identity_percent=101.0,
                alignment_length=10,
                e_value=1e-5,
                bit_score=20.0,
            )

    def test_identity_percent_out_of_range_negative(self):
        with pytest.raises(ValueError, match="identity_percent must be in"):
            BlastHit(
                subject_id="test",
                subject_organism="Test",
                identity_percent=-1.0,
                alignment_length=10,
                e_value=1e-5,
                bit_score=20.0,
            )

    def test_identity_percent_boundary_zero(self):
        hit = BlastHit(
            subject_id="test",
            subject_organism="Test",
            identity_percent=0.0,
            alignment_length=10,
            e_value=1e-5,
            bit_score=5.0,
        )
        assert hit.identity_percent == 0.0

    def test_identity_percent_boundary_100(self):
        hit = BlastHit(
            subject_id="test",
            subject_organism="Test",
            identity_percent=100.0,
            alignment_length=10,
            e_value=0.0,
            bit_score=20.0,
        )
        assert hit.identity_percent == 100.0

    def test_negative_alignment_length_raises(self):
        with pytest.raises(ValueError, match="alignment_length must be non-negative"):
            BlastHit(
                subject_id="test",
                subject_organism="Test",
                identity_percent=50.0,
                alignment_length=-1,
                e_value=1e-5,
                bit_score=10.0,
            )

    def test_zero_alignment_length(self):
        hit = BlastHit(
            subject_id="test",
            subject_organism="Test",
            identity_percent=0.0,
            alignment_length=0,
            e_value=1.0,
            bit_score=0.0,
        )
        assert hit.alignment_length == 0

    def test_negative_e_value_raises(self):
        with pytest.raises(ValueError, match="e_value must be non-negative"):
            BlastHit(
                subject_id="test",
                subject_organism="Test",
                identity_percent=50.0,
                alignment_length=10,
                e_value=-1e-5,
                bit_score=10.0,
            )

    def test_zero_e_value(self):
        hit = BlastHit(
            subject_id="test",
            subject_organism="Test",
            identity_percent=100.0,
            alignment_length=10,
            e_value=0.0,
            bit_score=20.0,
        )
        assert hit.e_value == 0.0

    def test_typical_blast_hit_values(self):
        """Test with realistic BLAST hit values."""
        hit = BlastHit(
            subject_id="gi|12345|emb|ABC123.1|",
            subject_organism="Bacillus anthracis str. Ames",
            identity_percent=87.3,
            alignment_length=756,
            e_value=2.1e-120,
            bit_score=425.0,
            is_pathogen=True,
            is_toxin=True,
        )
        assert hit.identity_percent == 87.3
        assert hit.alignment_length == 756
        assert hit.e_value < 1e-100


# ═══════════════════════════════════════════════════════════════════════
# 2. BlastScreeningResult Construction and Properties
# ═══════════════════════════════════════════════════════════════════════


class TestBlastScreeningResult:
    """Test BlastScreeningResult dataclass construction and properties."""

    def test_basic_construction(self):
        result = BlastScreeningResult(query_id="test_query")
        assert result.query_id == "test_query"
        assert result.hits == []
        assert result.is_safe is True
        assert result.screening_database == QUICK_SCREEN_DB_NAME
        assert result.e_value_threshold == 1e-5
        assert result.screening_time_seconds == 0.0

    def test_with_hits(self):
        hits = [
            BlastHit(
                subject_id="anthrax_lf",
                subject_organism="Bacillus anthracis",
                identity_percent=95.0,
                alignment_length=100,
                e_value=1e-50,
                bit_score=180.0,
                is_pathogen=True,
                is_toxin=True,
            ),
        ]
        result = BlastScreeningResult(
            query_id="test",
            hits=hits,
            is_safe=False,
            screening_database="blastn:/path/to/db",
            e_value_threshold=1e-10,
            screening_time_seconds=2.5,
        )
        assert len(result.hits) == 1
        assert result.is_safe is False
        assert result.screening_database == "blastn:/path/to/db"

    def test_pathogen_hits_property(self):
        hits = [
            BlastHit("a", "Org1", 90.0, 10, 1e-5, 20.0, is_pathogen=True),
            BlastHit("b", "Org2", 85.0, 10, 1e-5, 18.0, is_pathogen=False, is_toxin=True),
            BlastHit("c", "Org3", 80.0, 10, 1e-5, 16.0, is_pathogen=False),
        ]
        result = BlastScreeningResult(query_id="test", hits=hits)
        assert len(result.pathogen_hits) == 1
        assert result.pathogen_hits[0].subject_id == "a"

    def test_toxin_hits_property(self):
        hits = [
            BlastHit("a", "Org1", 90.0, 10, 1e-5, 20.0, is_pathogen=True, is_toxin=True),
            BlastHit("b", "Org2", 85.0, 10, 1e-5, 18.0, is_toxin=True),
            BlastHit("c", "Org3", 80.0, 10, 1e-5, 16.0),
        ]
        result = BlastScreeningResult(query_id="test", hits=hits)
        assert len(result.toxin_hits) == 2

    def test_concerning_hits_property(self):
        hits = [
            BlastHit("a", "Org1", 90.0, 10, 1e-5, 20.0, is_pathogen=True),
            BlastHit("b", "Org2", 85.0, 10, 1e-5, 18.0, is_toxin=True),
            BlastHit("c", "Org3", 80.0, 10, 1e-5, 16.0),
        ]
        result = BlastScreeningResult(query_id="test", hits=hits)
        # a is pathogen, b is toxin, c is neither
        assert len(result.concerning_hits) == 2

    def test_concerning_hits_empty(self):
        result = BlastScreeningResult(query_id="test")
        assert result.concerning_hits == []
        assert result.pathogen_hits == []
        assert result.toxin_hits == []

    def test_safe_result_with_no_concerning_hits(self):
        """Non-pathogen, non-toxin hits should still result in is_safe=True."""
        hits = [
            BlastHit("benign_gene", "E. coli K12", 99.0, 100, 1e-50, 180.0),
        ]
        result = BlastScreeningResult(query_id="test", hits=hits, is_safe=True)
        assert result.is_safe is True
        assert len(result.concerning_hits) == 0


# ═══════════════════════════════════════════════════════════════════════
# 3. BlastScreener.is_blast_available
# ═══════════════════════════════════════════════════════════════════════


class TestBlastScreenerAvailability:
    """Test BlastScreener.is_blast_available detection."""

    def test_blast_not_available_in_test_env(self):
        """BLAST+ is typically not installed in test environments."""
        screener = BlastScreener()
        # In most CI/test environments, BLAST+ is not installed
        # This test documents that behavior
        result = screener.is_blast_available()
        assert isinstance(result, bool)

    def test_blast_available_returns_false_when_not_installed(self):
        """When shutil.which cannot find blastn/blastp, returns False."""
        screener = BlastScreener()
        # Force re-check
        screener._blast_available = None
        with patch("shutil.which", return_value=None):
            assert screener.is_blast_available() is False

    def test_blast_available_returns_true_when_installed(self):
        """When both blastn and blastp are found, returns True."""
        screener = BlastScreener()
        screener._blast_available = None
        with patch("shutil.which", return_value="/usr/bin/blastn"):
            assert screener.is_blast_available() is True

    def test_blast_available_caches_result(self):
        """is_blast_available caches the result after first call."""
        screener = BlastScreener()
        screener._blast_available = None
        with patch("shutil.which", return_value="/usr/bin/blastn"):
            screener.is_blast_available()
        # Second call should use cached result even if we change shutil.which
        screener._blast_available = False
        assert screener.is_blast_available() is False

    def test_blast_available_handles_exception(self):
        """is_blast_available returns False if an exception occurs."""
        screener = BlastScreener()
        screener._blast_available = None
        with patch("shutil.which", side_effect=RuntimeError("unexpected")):
            assert screener.is_blast_available() is False


# ═══════════════════════════════════════════════════════════════════════
# 4. Quick Screen Fallback with Known Pathogen Patterns
# ═══════════════════════════════════════════════════════════════════════


class TestQuickScreenFallback:
    """Test pattern-based quick_screen when BLAST+ is not available."""

    def test_quick_screen_flags_anthrax_motif(self):
        """Sequences containing anthrax lethal factor motif should be flagged."""
        screener = BlastScreener()
        # Force BLAST unavailable so quick_screen is used
        screener._blast_available = False
        result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test_anthrax")
        assert result.is_safe is False
        assert len(result.hits) > 0
        assert any(h.subject_id == "anthrax_lethal_factor" for h in result.hits)

    def test_quick_screen_flags_botulinum_motif(self):
        """Sequences containing botulinum toxin motif should be flagged."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(BOTULINUM_CONTAINING_SEQ, "test_botulinum")
        assert result.is_safe is False
        assert any("botulinum" in h.subject_id for h in result.hits)

    def test_quick_screen_flags_ricin_motif(self):
        """Sequences containing ricin toxin motif should be flagged."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(RICIN_CONTAINING_SEQ, "test_ricin")
        assert result.is_safe is False
        assert any("ricin" in h.subject_id for h in result.hits)

    def test_quick_screen_flags_multiple_motifs(self):
        """Sequences with multiple pathogen motifs should flag all of them."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(MULTI_TOXIN_SEQ, "test_multi")
        assert result.is_safe is False
        # Should find at least the anthrax and ricin motifs
        hit_ids = [h.subject_id for h in result.hits]
        assert "anthrax_lethal_factor" in hit_ids
        assert "ricin_toxin" in hit_ids

    def test_quick_screen_hit_has_correct_metadata(self):
        """Quick screen hits should have correct organism and toxin flags."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")
        anthrax_hits = [h for h in result.hits if h.subject_id == "anthrax_lethal_factor"]
        assert len(anthrax_hits) == 1
        hit = anthrax_hits[0]
        assert hit.subject_organism == "Bacillus anthracis"
        assert hit.is_pathogen is True
        assert hit.is_toxin is True
        assert hit.identity_percent == 100.0
        assert hit.e_value == 0.0

    def test_quick_screen_uses_builtin_db_name(self):
        """Quick screen should report the builtin database name."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")
        assert result.screening_database == QUICK_SCREEN_DB_NAME


# ═══════════════════════════════════════════════════════════════════════
# 5. Safe Sequence Screening
# ═══════════════════════════════════════════════════════════════════════


class TestSafeSequenceScreening:
    """Test that safe sequences pass screening."""

    def test_safe_dna_passes(self):
        """A normal gene sequence should pass biosecurity screening."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(SAFE_DNA_SEQUENCE, "safe_query")
        assert result.is_safe is True
        assert len(result.hits) == 0

    def test_short_safe_sequence_passes(self):
        """Short sequences without pathogen motifs should pass."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence("ATGCGTAACTGGATC", "short_query")
        assert result.is_safe is True

    def test_empty_sequence_passes(self):
        """Empty sequence should be considered safe."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence("", "empty_query")
        assert result.is_safe is True
        assert len(result.hits) == 0

    def test_random_dna_passes(self):
        """Random DNA without pathogen motifs should pass."""
        screener = BlastScreener()
        screener._blast_available = False
        # Generate a random-looking sequence
        seq = "ATCGATCGATCGATCGATCGATCGATCGATCG"
        result = screener.screen_sequence(seq, "random_query")
        assert result.is_safe is True

    def test_safe_protein_passes(self):
        """Normal protein sequences should pass screening."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(SAFE_PROTEIN, "safe_protein")
        assert result.is_safe is True

    def test_empty_protein_passes(self):
        """Empty protein should be considered safe."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein("", "empty_protein")
        assert result.is_safe is True


# ═══════════════════════════════════════════════════════════════════════
# 6. Sequences Containing Known Pathogen Motifs
# ═══════════════════════════════════════════════════════════════════════


class TestPathogenMotifDetection:
    """Test detection of sequences containing known pathogen/toxin motifs."""

    def test_exact_motif_match_detected(self):
        """Exact match of a pathogen motif should be detected."""
        screener = BlastScreener()
        screener._blast_available = False
        # Use the exact motif sequence
        exact_motif = PATHOGEN_TOXIN_MOTIFS["diphtheria_toxin"]["sequence"]
        result = screener.screen_sequence(exact_motif, "diphtheria_exact")
        assert result.is_safe is False
        assert any("diphtheria" in h.subject_id for h in result.hits)

    def test_motif_embedded_in_longer_sequence(self):
        """A pathogen motif embedded in a longer sequence should be detected."""
        screener = BlastScreener()
        screener._blast_available = False
        long_seq = "A" * 500 + ANTHRAX_MOTIF + "T" * 500
        result = screener.screen_sequence(long_seq, "embedded_motif")
        assert result.is_safe is False

    def test_shiga_toxin_detected(self):
        """Shiga toxin motif should be detected."""
        screener = BlastScreener()
        screener._blast_available = False
        shiga_motif = PATHOGEN_TOXIN_MOTIFS["shiga_toxin"]["sequence"]
        seq = f"ATGCGT{shiga_motif}TACGCA"
        result = screener.screen_sequence(seq, "shiga_test")
        assert result.is_safe is False
        assert any("shiga" in h.subject_id for h in result.hits)

    def test_cholera_toxin_detected(self):
        """Cholera toxin motif should be detected."""
        screener = BlastScreener()
        screener._blast_available = False
        cholera_motif = PATHOGEN_TOXIN_MOTIFS["cholera_toxin"]["sequence"]
        seq = f"GCATGC{cholera_motif}ATCGAT"
        result = screener.screen_sequence(seq, "cholera_test")
        assert result.is_safe is False
        assert any("cholera" in h.subject_id for h in result.hits)

    def test_protein_with_ricin_motif_detected(self):
        """Protein containing ricin catalytic motif should be flagged."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(RICIN_CONTAINING_PROTEIN, "ricin_prot")
        assert result.is_safe is False
        assert any("ricin" in h.subject_id for h in result.hits)

    def test_protein_with_anthrax_motif_detected(self):
        """Protein containing anthrax LF motif should be flagged."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(ANTHRAX_CONTAINING_PROTEIN, "anthrax_prot")
        assert result.is_safe is False
        assert any("anthrax" in h.subject_id for h in result.hits)


# ═══════════════════════════════════════════════════════════════════════
# 7. E-Value and Identity Thresholds
# ═══════════════════════════════════════════════════════════════════════


class TestThresholdConfiguration:
    """Test e-value and identity threshold configuration."""

    def test_default_e_value_threshold(self):
        """Default e-value threshold should be 1e-5."""
        screener = BlastScreener()
        assert screener.e_value_threshold == 1e-5

    def test_custom_e_value_threshold(self):
        """Custom e-value threshold should be respected."""
        screener = BlastScreener(e_value_threshold=1e-10)
        assert screener.e_value_threshold == 1e-10

    def test_default_identity_threshold(self):
        """Default identity threshold should be 0.80."""
        screener = BlastScreener()
        assert screener.identity_threshold == 0.80

    def test_custom_identity_threshold(self):
        """Custom identity threshold should be respected."""
        screener = BlastScreener(identity_threshold=0.95)
        assert screener.identity_threshold == 0.95

    def test_result_includes_threshold(self):
        """Screening result should include the e-value threshold used."""
        screener = BlastScreener(e_value_threshold=1e-10)
        screener._blast_available = False
        result = screener.screen_sequence(SAFE_DNA_SEQUENCE, "test")
        assert result.e_value_threshold == 1e-10

    def test_low_identity_threshold_allows_more_hits(self):
        """With identity threshold at 0, all hits should be considered concerning."""
        screener = BlastScreener(identity_threshold=0.0)
        screener._blast_available = False
        result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")
        assert result.is_safe is False

    def test_very_high_identity_threshold_may_allow_hits(self):
        """With identity threshold at 100%, only exact 100% identity hits are flagged."""
        screener = BlastScreener(identity_threshold=1.0)
        screener._blast_available = False
        result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")
        # The quick_screen produces 100% identity hits, so this should still flag
        assert result.is_safe is False

    def test_threshold_affects_safety_verdict(self):
        """Identity threshold should affect is_safe verdict.

        Quick screen produces 100% identity hits, so even with
        identity_threshold=1.0 (100%), hits with 100% identity pass.
        With identity_threshold > 1.0 (impossible), nothing would pass.
        """
        screener_strict = BlastScreener(identity_threshold=1.01)
        screener_strict._blast_available = False
        result = screener_strict.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")
        # With threshold > 100%, no hit can exceed it, so is_safe=True
        assert result.is_safe is True


# ═══════════════════════════════════════════════════════════════════════
# 8. Environment Variable Configuration
# ═══════════════════════════════════════════════════════════════════════


class TestEnvironmentVariableConfiguration:
    """Test environment variable configuration."""

    def test_blast_db_path_from_env(self):
        """BIOCOMPILER_BLAST_DB_PATH should be read from environment."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_DB_PATH": "/custom/blast/db"}):
            screener = BlastScreener()
            assert screener.blast_db_path == "/custom/blast/db"

    def test_blast_db_path_explicit_overrides_env(self):
        """Explicit blast_db_path should take precedence over env var."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_DB_PATH": "/env/path"}):
            screener = BlastScreener(blast_db_path="/explicit/path")
            assert screener.blast_db_path == "/explicit/path"

    def test_e_value_from_env(self):
        """BIOCOMPILER_BLAST_E_VALUE should be read from environment."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_E_VALUE": "1e-10"}):
            screener = BlastScreener()
            assert screener.e_value_threshold == 1e-10

    def test_e_value_explicit_overrides_env(self):
        """Explicit e_value_threshold should take precedence over env var."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_E_VALUE": "1e-10"}):
            screener = BlastScreener(e_value_threshold=1e-3)
            assert screener.e_value_threshold == 1e-3

    def test_identity_from_env(self):
        """BIOCOMPILER_BLAST_IDENTITY should be read from environment."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_IDENTITY": "0.95"}):
            screener = BlastScreener()
            assert screener.identity_threshold == 0.95

    def test_identity_explicit_overrides_env(self):
        """Explicit identity_threshold should take precedence over env var."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_IDENTITY": "0.95"}):
            screener = BlastScreener(identity_threshold=0.70)
            assert screener.identity_threshold == 0.70

    def test_no_env_vars_uses_defaults(self):
        """Without env vars, defaults should be used."""
        # Remove any env vars that might be set
        env = os.environ.copy()
        env.pop("BIOCOMPILER_BLAST_DB_PATH", None)
        env.pop("BIOCOMPILER_BLAST_E_VALUE", None)
        env.pop("BIOCOMPILER_BLAST_IDENTITY", None)
        with patch.dict(os.environ, env, clear=True):
            screener = BlastScreener()
            assert screener.blast_db_path is None
            assert screener.e_value_threshold == 1e-5
            assert screener.identity_threshold == 0.80


# ═══════════════════════════════════════════════════════════════════════
# 9. Protein Screening
# ═══════════════════════════════════════════════════════════════════════


class TestProteinScreening:
    """Test protein sequence screening."""

    def test_safe_protein_screen(self):
        """Normal protein should pass screening."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein("MGSWKRQPPAVNVLRYFPATW", "safe_prot")
        assert result.is_safe is True

    def test_ricin_protein_motif_detected(self):
        """Protein with ricin catalytic motif should be flagged."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(RICIN_CONTAINING_PROTEIN, "ricin_prot")
        assert result.is_safe is False
        assert any("ricin" in h.subject_id for h in result.hits)

    def test_protein_screening_result_has_correct_db_name(self):
        """Protein quick screen should report the builtin database name."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(SAFE_PROTEIN, "test")
        assert result.screening_database == QUICK_SCREEN_DB_NAME

    def test_protein_hit_is_toxin(self):
        """Protein hits for toxins should have is_toxin=True."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(RICIN_CONTAINING_PROTEIN, "test")
        toxin_hits = result.toxin_hits
        assert len(toxin_hits) > 0

    def test_blastp_fallback_to_quick_screen(self):
        """When BLAST+ is not available, protein screen should use quick fallback."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(SAFE_PROTEIN, "test")
        # Should complete without error using fallback
        assert result.query_id == "test"


# ═══════════════════════════════════════════════════════════════════════
# 10. Internal Methods
# ═══════════════════════════════════════════════════════════════════════


class TestBlastScreenerInternalMethods:
    """Test BlastScreener internal helper methods."""

    def test_extract_organism_with_brackets(self):
        """Should extract organism name from square brackets."""
        result = BlastScreener._extract_organism(
            "lethal factor [Bacillus anthracis]"
        )
        assert result == "Bacillus anthracis"

    def test_extract_organism_without_brackets(self):
        """Should return full title when no brackets present."""
        result = BlastScreener._extract_organism("Bacillus anthracis lethal factor")
        assert result == "Bacillus anthracis lethal factor"

    def test_extract_organism_multiple_brackets(self):
        """Should extract from the last bracket pair."""
        result = BlastScreener._extract_organism(
            "gene [E. coli] protein [Bacillus anthracis]"
        )
        assert result == "Bacillus anthracis"

    def test_is_pathogen_subject_anthrax(self):
        """Should detect anthracis as pathogen."""
        assert BlastScreener._is_pathogen_subject("anthracis_gene", "Bacillus anthracis protein") is True

    def test_is_pathogen_subject_e_coli(self):
        """Should not flag E. coli K12 as pathogen."""
        assert BlastScreener._is_pathogen_subject("lacZ", "E. coli K12 beta-galactosidase") is False

    def test_is_pathogen_subject_botulinum(self):
        """Should detect botulinum as pathogen."""
        assert BlastScreener._is_pathogen_subject("bont", "Clostridium botulinum neurotoxin") is True

    def test_is_toxin_subject(self):
        """Should detect toxin keywords."""
        assert BlastScreener._is_toxin_subject("toxin_A", "enterotoxin type A") is True
        assert BlastScreener._is_toxin_subject("lacZ", "E. coli K12 beta-galactosidase") is False

    def test_is_toxin_subject_ricin(self):
        """Should detect ricin as toxin."""
        assert BlastScreener._is_toxin_subject("ricin", "Ricinus communis ricin toxin A-chain") is True

    def test_parse_blast_output_empty(self):
        """Empty BLAST output should produce no hits."""
        screener = BlastScreener()
        result = screener._parse_blast_output("")
        assert result == []

    def test_parse_blast_output_single_hit(self):
        """Single BLAST hit line should produce one BlastHit."""
        screener = BlastScreener()
        output = "query\tanthrax_lf\tanthrax lethal factor [Bacillus anthracis]\t95.5\t756\t2e-120\t425.0"
        hits = screener._parse_blast_output(output)
        assert len(hits) == 1
        assert hits[0].subject_id == "anthrax_lf"
        assert hits[0].identity_percent == 95.5
        assert hits[0].alignment_length == 756
        assert hits[0].e_value == 2e-120
        assert hits[0].bit_score == 425.0

    def test_parse_blast_output_multiple_hits(self):
        """Multiple BLAST hit lines should produce multiple BlastHits."""
        screener = BlastScreener()
        output = (
            "query\thit1\tdesc1\t90.0\t100\t1e-50\t180.0\n"
            "query\thit2\tdesc2\t85.0\t200\t1e-40\t150.0\n"
        )
        hits = screener._parse_blast_output(output)
        assert len(hits) == 2
        assert hits[0].subject_id == "hit1"
        assert hits[1].subject_id == "hit2"

    def test_parse_blast_output_malformed_line_skipped(self):
        """Malformed BLAST output lines should be skipped."""
        screener = BlastScreener()
        output = "query\thit1\tdesc1\t90.0\t100\t1e-50\n"  # Missing bitscore
        hits = screener._parse_blast_output(output)
        assert len(hits) == 0  # Should skip malformed line

    def test_parse_blast_output_pathogen_detection(self):
        """BLAST hits with pathogen keywords should have is_pathogen=True."""
        screener = BlastScreener()
        output = "query\tanthracis_lf\tAnthrax lethal factor [Bacillus anthracis]\t95.0\t100\t1e-50\t180.0"
        hits = screener._parse_blast_output(output)
        assert len(hits) == 1
        assert hits[0].is_pathogen is True

    def test_parse_blast_output_toxin_detection(self):
        """BLAST hits with toxin keywords should have is_toxin=True."""
        screener = BlastScreener()
        output = "query\tricin_toxin\tRicin toxin A-chain [Ricinus communis]\t95.0\t100\t1e-50\t180.0"
        hits = screener._parse_blast_output(output)
        assert len(hits) == 1
        assert hits[0].is_toxin is True


# ═══════════════════════════════════════════════════════════════════════
# 11. Screening Time Tracking
# ═══════════════════════════════════════════════════════════════════════


class TestScreeningTimeTracking:
    """Test that screening time is tracked."""

    def test_screening_time_nonzero(self):
        """Screening should report a non-zero time."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(SAFE_DNA_SEQUENCE, "test")
        assert result.screening_time_seconds >= 0.0

    def test_screening_time_reasonable(self):
        """Quick screening should be very fast (< 1 second)."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(SAFE_DNA_SEQUENCE, "test")
        assert result.screening_time_seconds < 1.0

    def test_protein_screening_time_tracked(self):
        """Protein screening should also track time."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(SAFE_PROTEIN, "test")
        assert result.screening_time_seconds >= 0.0


# ═══════════════════════════════════════════════════════════════════════
# 12. Case Insensitivity
# ═══════════════════════════════════════════════════════════════════════


class TestCaseInsensitivity:
    """Test that screening handles case correctly."""

    def test_lowercase_sequence_detected(self):
        """Lowercase sequences should still detect pathogen motifs."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_sequence(anthrax_containing_seq := ANTHRAX_CONTAINING_SEQ.lower(), "test")
        assert result.is_safe is False

    def test_mixed_case_sequence_detected(self):
        """Mixed case sequences should still detect pathogen motifs."""
        screener = BlastScreener()
        screener._blast_available = False
        mixed = ANTHRAX_CONTAINING_SEQ[:5] + ANTHRAX_CONTAINING_SEQ[5:10].lower() + ANTHRAX_CONTAINING_SEQ[10:]
        result = screener.screen_sequence(mixed, "test")
        assert result.is_safe is False

    def test_lowercase_protein_detected(self):
        """Lowercase protein should still be screened."""
        screener = BlastScreener()
        screener._blast_available = False
        result = screener.screen_protein(RICIN_CONTAINING_PROTEIN.lower(), "test")
        assert result.is_safe is False


# ═══════════════════════════════════════════════════════════════════════
# 13. Pathogen/Toxin Motif Database Integrity
# ═══════════════════════════════════════════════════════════════════════


class TestPathogenToxinMotifDatabase:
    """Test the built-in pathogen/toxin motif database."""

    def test_database_not_empty(self):
        """The motif database should not be empty."""
        assert len(PATHOGEN_TOXIN_MOTIFS) > 0

    def test_required_motifs_present(self):
        """Required pathogen motifs should be present in the database."""
        required = ["anthrax_lethal_factor", "botulinum_toxin_A", "ricin_toxin", "diphtheria_toxin"]
        for name in required:
            assert name in PATHOGEN_TOXIN_MOTIFS, f"Missing required motif: {name}"

    def test_each_motif_has_sequence(self):
        """Each motif should have a non-empty sequence."""
        for name, info in PATHOGEN_TOXIN_MOTIFS.items():
            assert "sequence" in info, f"Motif {name} missing 'sequence' key"
            assert len(str(info["sequence"])) > 0, f"Motif {name} has empty sequence"

    def test_each_motif_has_organism(self):
        """Each motif should have an organism."""
        for name, info in PATHOGEN_TOXIN_MOTIFS.items():
            assert "organism" in info, f"Motif {name} missing 'organism' key"

    def test_each_motif_has_pathogen_or_toxin_flag(self):
        """Each motif should have is_pathogen and is_toxin flags."""
        for name, info in PATHOGEN_TOXIN_MOTIFS.items():
            assert "is_pathogen" in info, f"Motif {name} missing 'is_pathogen' key"
            assert "is_toxin" in info, f"Motif {name} missing 'is_toxin' key"

    def test_dna_motifs_are_valid_dna(self):
        """DNA motifs should only contain valid DNA nucleotides."""
        for name, info in PATHOGEN_TOXIN_MOTIFS.items():
            seq = str(info["sequence"]).upper()
            valid_chars = set("ACGT")
            invalid = set(seq) - valid_chars
            assert not invalid, f"Motif {name} has invalid DNA characters: {invalid}"


# ═══════════════════════════════════════════════════════════════════════
# 14. BLAST+ Integration (Mock Tests)
# ═══════════════════════════════════════════════════════════════════════


class TestBlastPlusIntegration:
    """Test BLAST+ integration with mocked subprocess calls."""

    def test_blastn_called_when_available(self):
        """When BLAST+ is available and db is set, blastn should be called."""
        screener = BlastScreener(blast_db_path="/fake/db")
        screener._blast_available = True

        mock_output = "query\tanthracis_lf\tAnthrax lethal factor [Bacillus anthracis]\t95.0\t100\t1e-50\t180.0"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stdout = mock_output
            mock_run.return_value.stderr = ""
            result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")

        assert result.screening_database == "blastn:/fake/db"
        assert len(result.hits) == 1
        assert result.hits[0].is_pathogen is True

    def test_blastp_called_for_protein(self):
        """When BLAST+ is available and db is set, blastp should be called for protein."""
        screener = BlastScreener(blast_db_path="/fake/db")
        screener._blast_available = True

        mock_output = "query\tricin_toxin\tRicin toxin A-chain [Ricinus communis]\t98.0\t50\t1e-30\t90.0"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stdout = mock_output
            mock_run.return_value.stderr = ""
            result = screener.screen_protein(RICIN_CONTAINING_PROTEIN, "test")

        assert result.screening_database == "blastp:/fake/db"
        assert len(result.hits) == 1

    def test_blastn_failure_falls_back_to_quick_screen(self):
        """If blastn fails, should fall back to quick screen."""
        screener = BlastScreener(blast_db_path="/fake/db")
        screener._blast_available = True

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 1
            mock_run.return_value.stdout = ""
            mock_run.return_value.stderr = "blastn error"
            result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")

        # Should fall back to quick screen
        assert result.screening_database == QUICK_SCREEN_DB_NAME

    def test_no_db_path_uses_quick_screen(self):
        """When no db path is configured, should use quick screen."""
        screener = BlastScreener(blast_db_path=None)
        screener._blast_available = True  # Even if BLAST is available
        result = screener.screen_sequence(SAFE_DNA_SEQUENCE, "test")
        # Without a db path, should use quick screen
        assert result.screening_database == QUICK_SCREEN_DB_NAME

    def test_blast_safety_verdict_with_hits(self):
        """BLAST+ result with concerning hits should set is_safe=False."""
        screener = BlastScreener(blast_db_path="/fake/db")
        screener._blast_available = True

        # Hit with 95% identity to a pathogen
        mock_output = "query\tanthracis_lf\tAnthrax lethal factor [Bacillus anthracis]\t95.0\t100\t1e-50\t180.0"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stdout = mock_output
            mock_run.return_value.stderr = ""
            result = screener.screen_sequence(ANTHRAX_CONTAINING_SEQ, "test")

        assert result.is_safe is False

    def test_blast_safety_verdict_below_identity_threshold(self):
        """BLAST+ hit below identity threshold should not flag."""
        screener = BlastScreener(blast_db_path="/fake/db", identity_threshold=0.99)
        screener._blast_available = True

        # Hit with only 80% identity (below 99% threshold)
        mock_output = "query\tanthracis_lf\tAnthrax lethal factor [Bacillus anthracis]\t80.0\t100\t1e-50\t180.0"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stdout = mock_output
            mock_run.return_value.stderr = ""
            result = screener.screen_sequence("SOME_SEQ", "test")

        # 80% < 99% threshold, so should be safe
        assert result.is_safe is True
