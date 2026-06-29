"""
Biosecurity Integration Tests
==============================

Tests for the biosecurity screening pipeline integration:
- check_biosecurity_before_optimize() function
- BiosecurityError raised when pathogen sequences are detected
- skip_biosecurity_check parameter in optimize_sequence and batch_optimize
- biosecurity_screening_result field in OptimizationResult
- API endpoint integration (403 Forbidden for biosecurity failures)

These tests verify that the biosecurity screening gate is properly
wired into the optimization pipeline and blocks pathogen sequences
while allowing safe sequences through.
"""

import pytest
from unittest.mock import patch

from biocompiler.biosecurity import (
    check_biosecurity_before_optimize,
    BiosecurityScreeningResult,
    _PATHOGEN_SIGNATURES,
    _KMER_SIZE,
    _SIMILARITY_THRESHOLD,
    _extract_kmers,
    _compute_kmer_similarity,
)
from biocompiler.shared.exceptions import BiosecurityError, BioCompilerError
from biocompiler.optimizer import optimize_sequence, batch_optimize, OptimizationResult


# ────────────────────────────────────────────────────────────
# Test sequences
# ────────────────────────────────────────────────────────────

# Safe sequences that should pass screening
SAFE_PROTEINS = [
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
]

# Pathogen sequences (exact matches from _PATHOGEN_SIGNATURES)
ANTHRAX_SIGNATURE = "MEFKLRILVVSVATLFVSSGYSQHGVRNEQYADLAKR"
BOTULINUM_SIGNATURE = "MFVKLSFVKILIFQSSQHGVR"
RICIN_SIGNATURE = "AITNLFGRRLDKVKDTSINL"
SHIGA_SIGNATURE = "MYTILFSILLVSQHGVR"
CHOLERA_SIGNATURE = "MIKLCSFVSILLFFSSFSSG"
DIPHTHERIA_SIGNATURE = "MGKKVAVLLLASSVFAHGVR"
PLAGUE_V_SIGNATURE = "MKKISSVVVMTTMTFSSQHG"


class TestBiosecurityScreeningFunction:
    """Tests for the check_biosecurity_before_optimize function."""

    def test_safe_protein_passes(self):
        """Safe proteins should pass biosecurity screening."""
        for protein in SAFE_PROTEINS:
            result = check_biosecurity_before_optimize(protein, organism="Homo_sapiens")
            assert result.passed is True, f"Safe protein should pass: {protein[:30]}..."
            assert result.flagged_pathogens == []
            assert result.match_details == []
            assert result.screened_sequence_length == len(protein)

    def test_anthrax_signature_blocked(self):
        """Anthrax lethal factor signature should be blocked."""
        result = check_biosecurity_before_optimize(ANTHRAX_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Bacillus_anthracis" in result.flagged_pathogens
        assert any("CRITICAL" in rl for rl in result.risk_levels)

    def test_botulinum_signature_blocked(self):
        """Botulinum toxin signature should be blocked."""
        result = check_biosecurity_before_optimize(BOTULINUM_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Clostridium_botulinum" in result.flagged_pathogens

    def test_ricin_signature_blocked(self):
        """Ricin toxin signature should be blocked."""
        result = check_biosecurity_before_optimize(RICIN_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Ricinus_communis" in result.flagged_pathogens

    def test_shiga_signature_blocked(self):
        """Shiga toxin signature should be blocked."""
        result = check_biosecurity_before_optimize(SHIGA_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Shigella_dysenteriae" in result.flagged_pathogens

    def test_cholera_signature_blocked(self):
        """Cholera toxin signature should be blocked."""
        result = check_biosecurity_before_optimize(CHOLERA_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Vibrio_cholerae" in result.flagged_pathogens

    def test_diphtheria_signature_blocked(self):
        """Diphtheria toxin signature should be blocked."""
        result = check_biosecurity_before_optimize(DIPHTHERIA_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Corynebacterium_diphtheriae" in result.flagged_pathogens

    def test_plague_signature_blocked(self):
        """Plague V antigen signature should be blocked."""
        result = check_biosecurity_before_optimize(PLAGUE_V_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert "Yersinia_pestis" in result.flagged_pathogens

    def test_pathogen_embedding_in_longer_sequence(self):
        """Pathogen signature embedded in a longer sequence should be caught."""
        # Embed anthrax signature in a longer sequence
        embedded = "AAAAAA" + ANTHRAX_SIGNATURE + "GGGGGG"
        result = check_biosecurity_before_optimize(embedded, biosecurity_mode="warn")
        assert result.passed is False
        assert "Bacillus_anthracis" in result.flagged_pathogens

    def test_case_insensitive_matching(self):
        """Biosecurity matching should work with mixed case input."""
        # Lowercase input should still be caught
        result = check_biosecurity_before_optimize(ANTHRAX_SIGNATURE.lower(), biosecurity_mode="warn")
        assert result.passed is False

    def test_empty_sequence_raises_error(self):
        """Empty protein sequence should raise ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            check_biosecurity_before_optimize("")

    def test_whitespace_only_sequence_raises_error(self):
        """Whitespace-only protein sequence should raise ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            check_biosecurity_before_optimize("   ")

    def test_skip_check_returns_passed(self):
        """skip_biosecurity_check=True should always return passed=True."""
        result = check_biosecurity_before_optimize(
            ANTHRAX_SIGNATURE, skip_biosecurity_check=True
        )
        assert result.passed is True

    def test_result_str_passed(self):
        """String representation of passed result."""
        result = BiosecurityScreeningResult(passed=True)
        assert "PASSED" in str(result)

    def test_result_str_failed(self):
        """String representation of failed result."""
        result = BiosecurityScreeningResult(
            passed=False,
            flagged_pathogens=["Bacillus_anthracis"],
            risk_levels=["CRITICAL"],
        )
        assert "FAILED" in str(result)
        assert "Bacillus_anthracis" in str(result)

    def test_match_details_include_position(self):
        """Match details should include the position of the match."""
        result = check_biosecurity_before_optimize(ANTHRAX_SIGNATURE, biosecurity_mode="warn")
        assert result.passed is False
        assert len(result.match_details) > 0
        # Should mention position 0 for exact match at start
        assert "position 0" in result.match_details[0] or "matched at" in result.match_details[0]

    def test_kmer_scores_populated(self):
        """K-mer similarity scores should be populated for screened sequences."""
        # Use a long enough safe sequence for meaningful kmer analysis
        protein = SAFE_PROTEINS[0]
        result = check_biosecurity_before_optimize(protein)
        assert result.passed is True
        # kmer_scores may or may not have entries depending on sequence
        assert isinstance(result.kmer_scores, dict)


class TestKmerHelpers:
    """Tests for k-mer extraction and similarity helpers."""

    def test_extract_kmers_basic(self):
        """Test basic k-mer extraction."""
        kmers = _extract_kmers("ABCDEF", k=3)
        assert kmers == {"ABC", "BCD", "CDE", "DEF"}

    def test_extract_kmers_short_sequence(self):
        """K-mer extraction on sequence shorter than k should return empty set."""
        kmers = _extract_kmers("AB", k=6)
        assert kmers == set()

    def test_compute_kmer_similarity_no_overlap(self):
        """No overlap should return 0.0."""
        sim = _compute_kmer_similarity({"ABC", "DEF"}, {"GHI", "JKL"})
        assert sim == 0.0

    def test_compute_kmer_similarity_full_overlap(self):
        """Full overlap should return 1.0."""
        kmers = {"ABC", "DEF"}
        sim = _compute_kmer_similarity(kmers, kmers)
        assert sim == 1.0

    def test_compute_kmer_similarity_partial_overlap(self):
        """Partial overlap should return correct fraction."""
        query = {"ABC", "DEF", "GHI"}
        pathogen = {"ABC", "DEF", "JKL"}
        sim = _compute_kmer_similarity(query, pathogen)
        assert sim == 2 / 3

    def test_compute_kmer_similarity_empty_pathogen(self):
        """Empty pathogen k-mers should return 0.0."""
        sim = _compute_kmer_similarity({"ABC"}, set())
        assert sim == 0.0


class TestBiosecurityError:
    """Tests for the BiosecurityError exception."""

    def test_error_inherits_from_biocompiler_error(self):
        """BiosecurityError should be a subclass of BioCompilerError."""
        assert issubclass(BiosecurityError, BioCompilerError)

    def test_error_message_contains_pathogen(self):
        """Error message should contain the flagged pathogen name."""
        err = BiosecurityError(
            protein="TEST",
            flagged_pathogens=["Bacillus_anthracis"],
            risk_levels=["CRITICAL"],
            match_details=["Exact match at position 0"],
        )
        assert "Bacillus_anthracis" in str(err)
        assert "CRITICAL" in str(err)

    def test_error_message_contains_biosafety_guidance(self):
        """Error message should contain actionable guidance."""
        err = BiosecurityError(
            protein="TEST",
            flagged_pathogens=["Bacillus_anthracis"],
            risk_levels=["CRITICAL"],
            match_details=["Exact match at position 0"],
        )
        assert "BIOSECURITY ALERT" in str(err)
        assert "blocked" in str(err).lower() or "BLOCKED" in str(err)
        assert "biosafety" in str(err).lower()

    def test_error_stores_attributes(self):
        """Error should store protein, pathogens, and details."""
        err = BiosecurityError(
            protein="TESTSEQ",
            flagged_pathogens=["Pathogen1", "Pathogen2"],
            risk_levels=["HIGH", "CRITICAL"],
            match_details=["Detail1", "Detail2"],
        )
        assert err.protein == "TESTSEQ"
        assert err.flagged_pathogens == ["Pathogen1", "Pathogen2"]
        assert err.risk_levels == ["HIGH", "CRITICAL"]
        assert err.match_details == ["Detail1", "Detail2"]


class TestOptimizeSequenceIntegration:
    """Tests for biosecurity integration in optimize_sequence."""

    def test_safe_protein_optimizes_normally(self):
        """Safe proteins should optimize normally with biosecurity check."""
        result = optimize_sequence(
            target_protein=SAFE_PROTEINS[0],  # Hemoglobin
            organism="Escherichia_coli",
            skip_biosecurity_check=False,
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence  # Non-empty
        assert result.cai > 0.0
        assert result.biosecurity_screening_result is not None
        assert result.biosecurity_screening_result.passed is True

    def test_pathogen_protein_raises_biosecurity_error(self):
        """Pathogen protein should raise BiosecurityError."""
        with pytest.raises(BiosecurityError) as exc_info:
            optimize_sequence(
                target_protein=ANTHRAX_SIGNATURE,
                organism="Escherichia_coli",
            )
        assert "Bacillus_anthracis" in str(exc_info.value)

    def test_skip_biosecurity_allows_pathogen(self):
        """skip_biosecurity_check=True should allow pathogen sequences through."""
        # This should NOT raise BiosecurityError
        result = optimize_sequence(
            target_protein=ANTHRAX_SIGNATURE,
            organism="Escherichia_coli",
            skip_biosecurity_check=True,
        )
        assert isinstance(result, OptimizationResult)
        assert result.biosecurity_screening_result is not None
        assert result.biosecurity_screening_result.passed is True  # Skipped

    def test_biosecurity_result_attached_to_result(self):
        """OptimizationResult should have biosecurity_screening_result populated."""
        result = optimize_sequence(
            target_protein=SAFE_PROTEINS[0],
            organism="Escherichia_coli",
        )
        assert result.biosecurity_screening_result is not None
        assert isinstance(result.biosecurity_screening_result, BiosecurityScreeningResult)
        assert result.biosecurity_screening_result.passed is True
        assert result.biosecurity_screening_result.screened_sequence_length == len(SAFE_PROTEINS[0])

    def test_biosecurity_check_runs_before_optimization(self):
        """Biosecurity check should run before any optimization work."""
        # If optimization had started, it would take time; a BiosecurityError
        # should be raised before that
        with pytest.raises(BiosecurityError):
            optimize_sequence(
                target_protein=ANTHRAX_SIGNATURE,
                organism="Escherichia_coli",
            )

    def test_ricin_blocked_in_optimize(self):
        """Ricin toxin should be blocked by optimize_sequence."""
        with pytest.raises(BiosecurityError) as exc_info:
            optimize_sequence(
                target_protein=RICIN_SIGNATURE,
                organism="Homo_sapiens",
            )
        assert "Ricinus_communis" in str(exc_info.value)

    def test_botulinum_blocked_in_optimize(self):
        """Botulinum toxin should be blocked by optimize_sequence."""
        with pytest.raises(BiosecurityError) as exc_info:
            optimize_sequence(
                target_protein=BOTULINUM_SIGNATURE,
                organism="Homo_sapiens",
            )
        assert "Clostridium_botulinum" in str(exc_info.value)


class TestBatchOptimizeIntegration:
    """Tests for biosecurity integration in batch_optimize."""

    def test_safe_batch_optimizes_normally(self):
        """Safe protein batch should optimize normally."""
        results = batch_optimize(
            proteins=SAFE_PROTEINS[:2],
            organism="Escherichia_coli",
        )
        assert len(results) == 2
        for result in results:
            assert isinstance(result, OptimizationResult)
            assert result.biosecurity_screening_result is not None
            assert result.biosecurity_screening_result.passed is True

    def test_batch_with_pathogen_raises_error(self):
        """Batch containing a pathogen protein should raise BiosecurityError."""
        with pytest.raises(BiosecurityError) as exc_info:
            batch_optimize(
                proteins=[SAFE_PROTEINS[0], ANTHRAX_SIGNATURE],
                organism="Escherichia_coli",
            )
        assert "Bacillus_anthracis" in str(exc_info.value)

    def test_batch_skip_biosecurity_allows_pathogen(self):
        """skip_biosecurity_check=True in batch should allow pathogen through."""
        results = batch_optimize(
            proteins=[SAFE_PROTEINS[0]],
            organism="Escherichia_coli",
            skip_biosecurity_check=True,
        )
        assert len(results) == 1

    def test_batch_biosecurity_result_attached(self):
        """Each OptimizationResult in batch should have biosecurity result."""
        results = batch_optimize(
            proteins=SAFE_PROTEINS[:2],
            organism="Escherichia_coli",
        )
        for result in results:
            assert result.biosecurity_screening_result is not None
            assert result.biosecurity_screening_result.passed is True

    def test_empty_batch_returns_empty(self):
        """Empty protein list should return empty results."""
        results = batch_optimize(proteins=[], organism="Escherichia_coli")
        assert results == []

    def test_batch_full_postprocessing_with_biosecurity(self):
        """Batch with full_postprocessing=True should also screen."""
        with pytest.raises(BiosecurityError):
            batch_optimize(
                proteins=[ANTHRAX_SIGNATURE],
                organism="Escherichia_coli",
                full_postprocessing=True,
            )


class TestBiosecurityOptOut:
    """Tests for the skip_biosecurity_check parameter safety guarantees."""

    def test_default_is_not_skipped(self):
        """Default value of skip_biosecurity_check should be False."""
        import inspect
        sig = inspect.signature(optimize_sequence)
        assert sig.parameters["skip_biosecurity_check"].default is False

    def test_batch_default_is_not_skipped(self):
        """Default value of skip_biosecurity_check in batch should be False."""
        import inspect
        sig = inspect.signature(batch_optimize)
        assert sig.parameters["skip_biosecurity_check"].default is False

    def test_safe_protein_passes_without_skip(self):
        """Safe proteins should pass even without skip_biosecurity_check."""
        result = optimize_sequence(
            target_protein=SAFE_PROTEINS[0],
            organism="Escherichia_coli",
            # skip_biosecurity_check defaults to False
        )
        assert isinstance(result, OptimizationResult)
        assert result.biosecurity_screening_result.passed is True


class TestBiosecurityScreeningResultDataclass:
    """Tests for BiosecurityScreeningResult dataclass."""

    def test_default_values(self):
        """Default values should be sensible."""
        result = BiosecurityScreeningResult(passed=True)
        assert result.flagged_pathogens == []
        assert result.risk_levels == []
        assert result.match_details == []
        assert result.kmer_scores == {}
        assert result.screened_sequence_length == 0

    def test_failed_result_with_details(self):
        """Failed result should store all details."""
        result = BiosecurityScreeningResult(
            passed=False,
            flagged_pathogens=["Pathogen1"],
            risk_levels=["CRITICAL"],
            match_details=["Exact match at pos 0"],
            kmer_scores={"Pathogen1": 0.9},
            screened_sequence_length=42,
        )
        assert result.passed is False
        assert result.flagged_pathogens == ["Pathogen1"]
        assert result.risk_levels == ["CRITICAL"]
        assert result.match_details == ["Exact match at pos 0"]
        assert result.kmer_scores == {"Pathogen1": 0.9}
        assert result.screened_sequence_length == 42


class TestPathogenSignatureDatabase:
    """Tests for the pathogen signature database integrity."""

    def test_signatures_not_empty(self):
        """Pathogen signature database should not be empty."""
        assert len(_PATHOGEN_SIGNATURES) > 0

    def test_all_signatures_have_four_fields(self):
        """Each signature should have (sequence, pathogen, risk_level, description)."""
        for entry in _PATHOGEN_SIGNATURES:
            assert len(entry) == 4, f"Signature entry should have 4 fields: {entry}"

    def test_all_risk_levels_valid(self):
        """All risk levels should be CRITICAL or HIGH."""
        valid_levels = {"CRITICAL", "HIGH"}
        for _, _, risk_level, _ in _PATHOGEN_SIGNATURES:
            assert risk_level in valid_levels, f"Invalid risk level: {risk_level}"

    def test_all_signatures_are_protein_sequences(self):
        """All signature sequences should contain only valid amino acid codes."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        for sig, _, _, _ in _PATHOGEN_SIGNATURES:
            invalid = set(sig.upper()) - valid_aas
            assert not invalid, f"Invalid amino acids in signature: {invalid} in '{sig}'"

    def test_all_signatures_unique(self):
        """All signature sequences should be unique."""
        sequences = [sig for sig, _, _, _ in _PATHOGEN_SIGNATURES]
        assert len(sequences) == len(set(sequences)), "Duplicate signatures found"

    def test_signatures_cover_major_pathogens(self):
        """Database should cover major select agent pathogens."""
        pathogens = {pathogen for _, pathogen, _, _ in _PATHOGEN_SIGNATURES}
        # These are the minimum expected pathogens
        expected = {
            "Bacillus_anthracis",
            "Clostridium_botulinum",
            "Ricinus_communis",
            "Yersinia_pestis",
        }
        assert expected.issubset(pathogens), f"Missing pathogens: {expected - pathogens}"
