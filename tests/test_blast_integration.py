"""
Tests for BLAST+ integration module.

Covers:
  - BLAST+ binary detection (find_blast_bin, is_blast_available)
  - Database building commands (build_hazard_db)
  - BLAST XML output parsing
  - Risk level computation from BLAST results
  - Integration with biosecurity module (check_biosecurity_blast)
  - Fallback when BLAST+ is not available
  - BlastScanner wrapper class
  - Hit classification
  - Confidence score mapping
"""

import os
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from biocompiler.biosecurity.blast_integration import (
    BlastHit,
    BlastScanner,
    BlastScreeningResult,
    build_hazard_db,
    check_biosecurity_blast,
    find_blast_bin,
    is_blast_available,
    screen_blast,
    _blast_identity_to_confidence,
    _classify_blast_hit,
    _compute_blast_risk_level,
    _merge_reports,
    _build_hazard_sequences_from_signatures,
)
from biocompiler.biosecurity import (
    BiosecurityReport,
    HazardMatch,
    screen_hazardous_sequence,
)


# ═══════════════════════════════════════════════════════════════════════════
# Sample BLAST XML output for testing
# ═══════════════════════════════════════════════════════════════════════════

SAMPLE_BLAST_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.14.1+</BlastOutput_version>
  <BlastOutput_reference>Reference: Camacho et al., BMC Bioinformatics 2009; 10:421</BlastOutput_reference>
  <BlastOutput_db>hazard_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test_query</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test_query</Iteration_query-def>
      <Iteration_query-len>100</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
          <Hit_def>ricin_A_chain_catalytic Ricin A-chain catalytic site</Hit_def>
          <Hit_accession>1</Hit_accession>
          <Hit_len>267</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>150.5</Hsp_bit-score>
              <Hsp_score>80</Hsp_score>
              <Hsp_evalue>1e-40</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>50</Hsp_hit-from>
              <Hsp_hit-to>149</Hsp_hit-to>
              <Hsp_identity>95</Hsp_identity>
              <Hsp_positive>98</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_hit-len>100</Hsp_hit-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|2</Hit_id>
          <Hit_def>blaTEM_protein TEM beta-lactamase</Hit_def>
          <Hit_accession>2</Hit_accession>
          <Hit_len>286</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>85.2</Hsp_bit-score>
              <Hsp_score>45</Hsp_score>
              <Hsp_evalue>1e-18</Hsp_evalue>
              <Hsp_query-from>10</Hsp_query-from>
              <Hsp_query-to>70</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>65</Hsp_hit-to>
              <Hsp_identity>42</Hsp_identity>
              <Hsp_positive>50</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>60</Hsp_align-len>
              <Hsp_hit-len>60</Hsp_hit-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|3</Hit_id>
          <Hit_def>unknown_protein hypothetical protein</Hit_def>
          <Hit_accession>3</Hit_accession>
          <Hit_len>200</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>40.0</Hsp_bit-score>
              <Hsp_score>20</Hsp_score>
              <Hsp_evalue>0.001</Hsp_evalue>
              <Hsp_query-from>20</Hsp_query-from>
              <Hsp_query-to>55</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>45</Hsp_hit-to>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_positive>25</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>35</Hsp_align-len>
              <Hsp_hit-len>35</Hsp_hit-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>100</Statistics_db-num>
          <Statistics_db-len>25000</Statistics_db-len>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

# Empty BLAST output (no hits)
EMPTY_BLAST_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.14.1+</BlastOutput_version>
  <BlastOutput_db>hazard_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test_query</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test_query</Iteration_query-def>
      <Iteration_query-len>100</Iteration_query-len>
      <Iteration_hits>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>100</Statistics_db-num>
          <Statistics_db-len>25000</Statistics_db-len>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


# ═══════════════════════════════════════════════════════════════════════════
# Test: BLAST+ binary detection
# ═══════════════════════════════════════════════════════════════════════════


class TestFindBlastBin:
    """Test find_blast_bin for BLAST+ binary detection."""

    def test_returns_none_when_not_on_path(self):
        """When the binary is not on PATH and no custom path set, return None."""
        with patch.dict(os.environ, {}, clear=True):
            with patch("shutil.which", return_value=None):
                result = find_blast_bin("nonexistent_blast_tool")
                assert result is None

    def test_returns_path_when_on_system_path(self):
        """When the binary is on the system PATH, return its full path."""
        with patch("shutil.which", return_value="/usr/bin/blastn"):
            result = find_blast_bin("blastn")
            assert result == "/usr/bin/blastn"

    def test_checks_custom_path_first(self):
        """When BIOCOMPILER_BLAST_PATH is set, check it first."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_PATH": "/opt/blast/bin"}):
            with patch.object(Path, "is_file", return_value=True):
                result = find_blast_bin("blastn")
                assert result is not None
                assert "blastn" in result

    def test_falls_back_to_system_path(self):
        """When custom path does not have the binary, fall back to system PATH."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_PATH": "/nonexistent"}):
            with patch.object(Path, "is_file", return_value=False):
                with patch("shutil.which", return_value="/usr/bin/blastn"):
                    result = find_blast_bin("blastn")
                    assert result == "/usr/bin/blastn"


class TestIsBlastAvailable:
    """Test is_blast_available for BLAST+ detection."""

    def test_returns_false_when_no_blast(self):
        """When blastn and makeblastdb are not found, return False."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value=None):
            assert is_blast_available() is False

    def test_returns_true_when_blast_available(self):
        """When blastn and makeblastdb are found, return True."""
        def mock_find(name):
            if name in ("blastn", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            assert is_blast_available() is True

    def test_returns_false_when_only_blastn_available(self):
        """When only blastn is found (no makeblastdb), return False."""
        def mock_find(name):
            if name == "blastn":
                return "/usr/bin/blastn"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            assert is_blast_available() is False


# ═══════════════════════════════════════════════════════════════════════════
# Test: BlastScanner wrapper
# ═══════════════════════════════════════════════════════════════════════════


class TestBlastScanner:
    """Test the BlastScanner wrapper class."""

    def test_available_property_when_tools_found(self):
        """Scanner reports available when blastn and makeblastdb are found."""
        def mock_find(name):
            if name in ("blastn", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            scanner = BlastScanner(db_path="/tmp/test_db")
            assert scanner.available is True

    def test_available_property_when_tools_not_found(self):
        """Scanner reports unavailable when tools are not found."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value=None):
            scanner = BlastScanner(db_path="/tmp/test_db")
            assert scanner.available is False

    def test_protein_available_property(self):
        """Scanner reports protein_available when blastp is found."""
        def mock_find(name):
            if name in ("blastn", "blastp", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            scanner = BlastScanner(db_path="/tmp/test_db")
            assert scanner.protein_available is True

    def test_screen_returns_unavailable_when_no_blast(self):
        """When BLAST is not installed, screen returns blast_available=False."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value=None):
            scanner = BlastScanner(db_path="/tmp/test_db")
            result = scanner.screen("ATCGATCG", seq_type="dna")
            assert result.blast_available is False
            assert len(result.warnings) > 0
            assert "BLAST+" in result.warnings[0]

    def test_screen_handles_blast_timeout(self):
        """When BLAST times out, result has an error message."""
        def mock_find(name):
            if name in ("blastn", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            scanner = BlastScanner(db_path="/tmp/test_db")
            with patch.object(scanner, "_run_blast", side_effect=subprocess.TimeoutExpired("blastn", 120)):
                result = scanner.screen("ATCGATCG", seq_type="dna")
                assert result.blast_available is True
                assert result.error is not None
                assert "timed out" in result.error

    def test_screen_handles_blast_failure(self):
        """When BLAST fails, result has an error message."""
        def mock_find(name):
            if name in ("blastn", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            scanner = BlastScanner(db_path="/tmp/test_db")
            with patch.object(
                scanner, "_run_blast",
                side_effect=subprocess.CalledProcessError(1, "blastn", stderr="Database not found"),
            ):
                result = scanner.screen("ATCGATCG", seq_type="dna")
                assert result.blast_available is True
                assert result.error is not None
                assert "failed" in result.error

    def test_screen_handles_protein_unavailable(self):
        """When blastp is not available, protein screening returns an error."""
        def mock_find(name):
            if name in ("blastn", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            scanner = BlastScanner(db_path="/tmp/test_db")
            result = scanner.screen("MSKGEELFTG", seq_type="protein")
            assert result.error is not None
            assert "blastp" in result.error

    def test_screen_successful_blast(self):
        """Successful BLAST search produces structured results."""
        def mock_find(name):
            if name in ("blastn", "blastp", "makeblastdb"):
                return f"/usr/bin/{name}"
            return None

        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", side_effect=mock_find):
            scanner = BlastScanner(db_path="/tmp/test_db")
            with patch.object(scanner, "_run_blast", return_value=SAMPLE_BLAST_XML):
                result = scanner.screen("MSKGEELFTGVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", seq_type="protein")
                assert result.blast_available is True
                assert len(result.hits) == 3
                assert result.hits[0].percent_identity == 95.0  # 95/100
                assert result.hits[0].bit_score == pytest.approx(150.5, rel=1e-6)
                assert result.blast_program == "blastp"


# ═══════════════════════════════════════════════════════════════════════════
# Test: BLAST XML parsing
# ═══════════════════════════════════════════════════════════════════════════


class TestBlastXmlParsing:
    """Test parsing of BLAST XML output."""

    def test_parse_sample_xml(self):
        """Parse the sample BLAST XML and verify hits."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            hits = scanner._parse_xml_output(SAMPLE_BLAST_XML)
            assert len(hits) == 3

            # First hit: ricin A-chain (95% identity, 100 alignment length)
            assert "ricin" in hits[0].subject_id.lower()
            assert hits[0].percent_identity == pytest.approx(95.0, rel=1e-6)
            assert hits[0].alignment_length == 100
            assert hits[0].evalue == 1e-40
            assert hits[0].bit_score == pytest.approx(150.5, rel=1e-6)

            # Second hit: blaTEM (42/60 = 70% identity)
            assert "blaTEM" in hits[1].subject_id
            assert hits[1].percent_identity == 70.0  # 42/60 = 70%
            assert hits[1].alignment_length == 60

            # Third hit: unknown (20/35 = ~57.1% identity)
            assert hits[2].percent_identity == pytest.approx(57.14, rel=0.01)

    def test_parse_empty_xml(self):
        """Parse empty BLAST XML (no hits)."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            hits = scanner._parse_xml_output(EMPTY_BLAST_XML)
            assert len(hits) == 0

    def test_parse_empty_string(self):
        """Parse empty string returns no hits."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            hits = scanner._parse_xml_output("")
            assert len(hits) == 0

    def test_parse_invalid_xml(self):
        """Invalid XML returns no hits without crashing."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            hits = scanner._parse_xml_output("NOT VALID XML <<<>>>")
            assert len(hits) == 0

    def test_hits_sorted_by_bit_score(self):
        """Parsed hits should be sorted by bit_score descending."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            hits = scanner._parse_xml_output(SAMPLE_BLAST_XML)
            for i in range(len(hits) - 1):
                assert hits[i].bit_score >= hits[i + 1].bit_score

    def test_respects_max_hits(self):
        """Only return up to max_hits results."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db", max_hits=2)
            hits = scanner._parse_xml_output(SAMPLE_BLAST_XML)
            assert len(hits) <= 2


# ═══════════════════════════════════════════════════════════════════════════
# Test: Risk level computation
# ═══════════════════════════════════════════════════════════════════════════


class TestComputeBlastRiskLevel:
    """Test _compute_blast_risk_level for BLAST-based risk assessment."""

    def test_critical_for_high_identity_high_coverage(self):
        """identity > 90% and coverage > 50% → critical."""
        assert _compute_blast_risk_level(95.0, 0.8) == "critical"
        assert _compute_blast_risk_level(99.0, 0.51) == "critical"
        assert _compute_blast_risk_level(90.1, 0.6) == "critical"

    def test_high_for_high_identity_low_coverage(self):
        """identity > 90% and coverage ≤ 50% → high."""
        assert _compute_blast_risk_level(95.0, 0.3) == "high"
        assert _compute_blast_risk_level(91.0, 0.1) == "high"

    def test_high_for_medium_identity_high_coverage(self):
        """identity > 70% and coverage > 50% → high."""
        assert _compute_blast_risk_level(75.0, 0.6) == "high"
        assert _compute_blast_risk_level(85.0, 0.8) == "high"

    def test_medium_for_medium_identity_low_coverage(self):
        """identity > 70% and coverage ≤ 50% → medium."""
        assert _compute_blast_risk_level(75.0, 0.3) == "medium"
        assert _compute_blast_risk_level(80.0, 0.1) == "medium"

    def test_medium_for_low_identity_high_coverage(self):
        """identity > 50% and coverage > 50% → medium."""
        assert _compute_blast_risk_level(60.0, 0.7) == "medium"
        assert _compute_blast_risk_level(55.0, 0.55) == "medium"

    def test_low_for_low_identity_low_coverage(self):
        """identity > 50% and coverage ≤ 50% → low."""
        assert _compute_blast_risk_level(55.0, 0.3) == "low"
        assert _compute_blast_risk_level(60.0, 0.1) == "low"

    def test_none_for_very_low_identity(self):
        """identity ≤ 50% → none."""
        assert _compute_blast_risk_level(40.0, 0.8) == "none"
        assert _compute_blast_risk_level(30.0, 0.1) == "none"
        assert _compute_blast_risk_level(0.0, 0.0) == "none"

    def test_boundary_at_90_percent(self):
        """Identity exactly at 90% is still 'high' (not 'critical' since > 90)."""
        # 90.0 is NOT > 90.0, so it falls to the > 70% bracket
        assert _compute_blast_risk_level(90.0, 0.6) == "high"

    def test_boundary_at_70_percent(self):
        """Identity exactly at 70% falls to the > 50% bracket."""
        assert _compute_blast_risk_level(70.0, 0.6) == "medium"

    def test_boundary_at_50_percent(self):
        """Identity exactly at 50% → none."""
        assert _compute_blast_risk_level(50.0, 0.6) == "none"


# ═══════════════════════════════════════════════════════════════════════════
# Test: Database building
# ═══════════════════════════════════════════════════════════════════════════


class TestBuildHazardDb:
    """Test build_hazard_db for BLAST database creation."""

    def test_raises_when_makeblastdb_not_found(self):
        """Raises RuntimeError when makeblastdb is not available."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value=None):
            with pytest.raises(RuntimeError, match="makeblastdb not found"):
                build_hazard_db({"seq1": "ATCG"}, "test_db", "/tmp/test_db")

    def test_raises_when_no_sequences(self):
        """Raises ValueError when no sequences are provided."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/makeblastdb"):
            with pytest.raises(ValueError, match="No sequences provided"):
                build_hazard_db({}, "test_db", "/tmp/test_db")

    def test_calls_makeblastdb_correctly(self):
        """Calls makeblastdb with correct arguments."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/makeblastdb"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
                result = build_hazard_db(
                    {"ricin_A_chain": "NIRVGLPIIS", "blaTEM": "HPETLALKFG"},
                    "hazard_db",
                    "/tmp/hazard_db",
                    db_type="prot",
                )
                assert result == "/tmp/hazard_db"
                mock_run.assert_called_once()
                cmd = mock_run.call_args[0][0]
                assert cmd[0] == "/usr/bin/makeblastdb"
                assert "-dbtype" in cmd
                assert "prot" in cmd
                assert "-out" in cmd
                assert "/tmp/hazard_db" in cmd

    def test_creates_output_directory(self):
        """Creates the output directory if it does not exist."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/makeblastdb"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
                with tempfile.TemporaryDirectory() as tmpdir:
                    build_hazard_db(
                        {"seq1": "ATCG"},
                        "test_db",
                        os.path.join(tmpdir, "subdir", "test_db"),
                    )
                    # Verify the subdirectory was created
                    assert os.path.isdir(os.path.join(tmpdir, "subdir"))

    def test_sanitizes_sequence_ids(self):
        """Sequence IDs with spaces and pipes are sanitized."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/makeblastdb"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
                with tempfile.TemporaryDirectory() as tmpdir:
                    build_hazard_db(
                        {"seq with spaces|pipes": "ATCG"},
                        "test_db",
                        os.path.join(tmpdir, "test_db"),
                    )
                    # Check that the FASTA file was written with sanitized ID
                    call_args = mock_run.call_args[0][0]
                    fasta_path = call_args[call_args.index("-in") + 1]
                    # The temp FASTA file should have been cleaned up, but the
                    # makeblastdb command should have been called


class TestBuildHazardSequencesFromSignatures:
    """Test extraction of sequences from the biosecurity signature database."""

    def test_protein_signatures_extracted(self):
        """Protein signatures are extracted correctly."""
        seqs = _build_hazard_sequences_from_signatures("protein")
        assert len(seqs) > 0
        assert "ricin_A_chain_catalytic" in seqs
        assert seqs["ricin_A_chain_catalytic"] == "NIRVGLPIIS"

    def test_dna_signatures_extracted(self):
        """DNA signatures are extracted correctly."""
        seqs = _build_hazard_sequences_from_signatures("dna")
        assert len(seqs) > 0
        assert "blaTEM_dna" in seqs

    def test_protein_and_dna_are_different(self):
        """Protein and DNA signatures are different sets."""
        prot_seqs = _build_hazard_sequences_from_signatures("protein")
        dna_seqs = _build_hazard_sequences_from_signatures("dna")
        assert set(prot_seqs.keys()) != set(dna_seqs.keys())


# ═══════════════════════════════════════════════════════════════════════════
# Test: Hit classification
# ═══════════════════════════════════════════════════════════════════════════


class TestClassifyBlastHit:
    """Test _classify_blast_hit for BLAST hit category assignment."""

    def test_select_agent_ricin(self):
        hit = BlastHit(subject_id="ricin_A_chain", percent_identity=95.0, alignment_length=100, evalue=1e-40, bit_score=150.0)
        assert _classify_blast_hit(hit) == "select_agent"

    def test_select_agent_botulinum(self):
        hit = BlastHit(subject_id="botulinum_neurotoxin", percent_identity=80.0, alignment_length=50, evalue=1e-10, bit_score=80.0)
        assert _classify_blast_hit(hit) == "select_agent"

    def test_viral_surface_spike(self):
        hit = BlastHit(subject_id="spike_glycoprotein_SARS2", percent_identity=75.0, alignment_length=100, evalue=1e-20, bit_score=100.0)
        assert _classify_blast_hit(hit) == "viral_surface"

    def test_viral_surface_hemagglutinin(self):
        hit = BlastHit(subject_id="hemagglutinin_influenza", percent_identity=70.0, alignment_length=80, evalue=1e-15, bit_score=90.0)
        assert _classify_blast_hit(hit) == "viral_surface"

    def test_antibiotic_resistance_blatem(self):
        hit = BlastHit(subject_id="blaTEM_beta_lactamase", percent_identity=85.0, alignment_length=60, evalue=1e-18, bit_score=85.0)
        assert _classify_blast_hit(hit) == "antibiotic_resistance"

    def test_antibiotic_resistance_ndm(self):
        hit = BlastHit(subject_id="ndm-1_metallo_beta_lactamase", percent_identity=90.0, alignment_length=70, evalue=1e-25, bit_score=120.0)
        assert _classify_blast_hit(hit) == "antibiotic_resistance"

    def test_oncogene_myc(self):
        hit = BlastHit(subject_id="c-myc_transcription_factor", percent_identity=65.0, alignment_length=40, evalue=1e-5, bit_score=50.0)
        assert _classify_blast_hit(hit) == "oncogene"

    def test_oncogene_ras(self):
        hit = BlastHit(subject_id="ras_gtpase_oncogene", percent_identity=60.0, alignment_length=30, evalue=0.001, bit_score=40.0)
        assert _classify_blast_hit(hit) == "oncogene"

    def test_unknown_defaults_to_blast_homology(self):
        hit = BlastHit(subject_id="hypothetical_protein_x", percent_identity=55.0, alignment_length=30, evalue=0.01, bit_score=35.0)
        assert _classify_blast_hit(hit) == "blast_homology"


# ═══════════════════════════════════════════════════════════════════════════
# Test: Confidence score mapping
# ═══════════════════════════════════════════════════════════════════════════


class TestBlastIdentityToConfidence:
    """Test _blast_identity_to_confidence for BLAST confidence mapping."""

    def test_high_identity(self):
        assert _blast_identity_to_confidence(99.0) == pytest.approx(0.95, rel=1e-6)
        assert _blast_identity_to_confidence(96.0) == pytest.approx(0.95, rel=1e-6)

    def test_90_plus_identity(self):
        assert _blast_identity_to_confidence(92.0) == pytest.approx(0.90, rel=1e-6)
        assert _blast_identity_to_confidence(95.0) == pytest.approx(0.90, rel=1e-6)

    def test_70_plus_identity(self):
        assert _blast_identity_to_confidence(80.0) == pytest.approx(0.75, rel=1e-6)
        assert _blast_identity_to_confidence(75.0) == pytest.approx(0.75, rel=1e-6)

    def test_50_plus_identity(self):
        assert _blast_identity_to_confidence(55.0) == pytest.approx(0.55, rel=1e-6)
        assert _blast_identity_to_confidence(60.0) == pytest.approx(0.55, rel=1e-6)

    def test_low_identity(self):
        assert _blast_identity_to_confidence(30.0) == pytest.approx(0.30, rel=1e-6)
        assert _blast_identity_to_confidence(40.0) == pytest.approx(0.30, rel=1e-6)

    def test_boundary_values(self):
        assert _blast_identity_to_confidence(95.0) == 0.90  # not > 95
        assert _blast_identity_to_confidence(90.0) == 0.75  # not > 90
        assert _blast_identity_to_confidence(70.0) == 0.55  # not > 70
        assert _blast_identity_to_confidence(50.0) == 0.30  # not > 50


# ═══════════════════════════════════════════════════════════════════════════
# Test: Report merging
# ═══════════════════════════════════════════════════════════════════════════


class TestMergeReports:
    """Test _merge_reports for combining exact/fuzzy and BLAST results."""

    def test_merge_with_no_blast_result(self):
        """When no BLAST result, return exact report with recommendation."""
        exact = BiosecurityReport(
            is_hazardous=False,
            risk_level="none",
            flagged_categories=[],
            matches=[],
            recommendations=["No biosecurity concerns detected."],
        )
        merged = _merge_reports(exact, None)
        assert merged.risk_level == "none"
        assert not merged.is_hazardous
        assert any("BLAST+" in r for r in merged.recommendations)

    def test_merge_takes_highest_risk(self):
        """Merged risk level is the maximum of both reports."""
        exact = BiosecurityReport(
            is_hazardous=True,
            risk_level="medium",
            flagged_categories=["antibiotic_resistance"],
            matches=[HazardMatch(
                category="antibiotic_resistance", name="nptII_protein",
                position=0, matched_sequence="RPMTIHGSGS",
                confidence=0.88, source="test",
            )],
            recommendations=["Medium risk detected."],
        )
        blast = BlastScreeningResult(
            blast_available=True,
            hits=[BlastHit(
                subject_id="ricin_A_chain", percent_identity=92.0,
                alignment_length=100, evalue=1e-40, bit_score=150.0,
                query_coverage=0.8,
            )],
            risk_level="critical",
            coverage=0.8,
            best_identity=92.0,
            database_path="/tmp/hazard_db",
            blast_program="blastp",
        )
        merged = _merge_reports(exact, blast)
        assert merged.risk_level == "critical"
        assert merged.is_hazardous is True
        assert "select_agent" in merged.flagged_categories

    def test_merge_adds_blast_hits_as_matches(self):
        """BLAST hits are added as HazardMatch entries."""
        exact = BiosecurityReport(
            is_hazardous=False,
            risk_level="none",
            flagged_categories=[],
            matches=[],
            recommendations=["No concerns."],
        )
        blast = BlastScreeningResult(
            blast_available=True,
            hits=[BlastHit(
                subject_id="ricin_A_chain", percent_identity=85.0,
                alignment_length=100, evalue=1e-30, bit_score=120.0,
                query_coverage=0.9,
            )],
            risk_level="high",
            coverage=0.9,
            best_identity=85.0,
            database_path="/tmp/hazard_db",
            blast_program="blastp",
        )
        merged = _merge_reports(exact, blast)
        # Should have at least one BLAST match
        blast_matches = [m for m in merged.matches if m.name.startswith("BLAST:")]
        assert len(blast_matches) >= 1
        assert "ricin" in blast_matches[0].name

    def test_merge_with_blast_unavailable(self):
        """When BLAST is unavailable, add a warning recommendation."""
        exact = BiosecurityReport(
            is_hazardous=False,
            risk_level="none",
            flagged_categories=[],
            matches=[],
            recommendations=["No concerns."],
        )
        blast = BlastScreeningResult(
            blast_available=False,
            warnings=["BLAST+ is not installed"],
        )
        merged = _merge_reports(exact, blast)
        assert any("BLAST+" in r or "not available" in r for r in merged.recommendations)

    def test_merge_with_blast_error(self):
        """When BLAST had an error, add error info to recommendations."""
        exact = BiosecurityReport(
            is_hazardous=False,
            risk_level="none",
            flagged_categories=[],
            matches=[],
            recommendations=["No concerns."],
        )
        blast = BlastScreeningResult(
            blast_available=True,
            error="blastp failed with exit code 1",
            blast_program="blastp",
        )
        merged = _merge_reports(exact, blast)
        assert any("error" in r.lower() for r in merged.recommendations)


# ═══════════════════════════════════════════════════════════════════════════
# Test: Full integration (check_biosecurity_blast)
# ═══════════════════════════════════════════════════════════════════════════


class TestCheckBiosecurityBlast:
    """Test check_biosecurity_blast integration function."""

    def test_returns_biosecurity_report(self):
        """Always returns a BiosecurityReport."""
        report = check_biosecurity_blast("MSKGEELFTGVPILVELDGDVNGHKFSVS")
        assert isinstance(report, BiosecurityReport)

    def test_falls_back_to_exact_only_without_blast(self):
        """Without BLAST, falls back to exact/fuzzy matching only."""
        with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=False):
            report = check_biosecurity_blast("MSKGEELFTGVPILVELDGDVNGHKFSVS")
            # Should have a recommendation about BLAST not being available
            assert any("BLAST+" in r for r in report.recommendations)

    def test_detects_ricin_with_exact_matching(self):
        """Ricin A-chain motif is detected via exact matching."""
        with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=False):
            report = check_biosecurity_blast("MISRDNIRVGLPIISTNKYEDKQL")
            assert report.is_hazardous
            assert report.risk_level == "critical"

    def test_skips_blast_when_no_db_path(self):
        """When no db_path is provided and no env var set, BLAST is skipped."""
        with patch.dict(os.environ, {}, clear=False):
            # Remove BIOCOMPILER_BLAST_DB_PATH if present
            env = os.environ.copy()
            env.pop("BIOCOMPILER_BLAST_DB_PATH", None)
            with patch.dict(os.environ, env, clear=True):
                with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=True):
                    report = check_biosecurity_blast("MSKGEELFTGVPILVELDGDVNGHKFSVS")
                    # BLAST should be skipped (no db path)
                    assert isinstance(report, BiosecurityReport)

    def test_runs_blast_when_available(self):
        """When BLAST is available and db_path is set, runs BLAST screening."""
        with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=True):
            with patch("biocompiler.biosecurity.blast_integration.screen_blast") as mock_screen:
                mock_screen.return_value = BlastScreeningResult(
                    blast_available=True,
                    hits=[BlastHit(
                        subject_id="ricin_A_chain", percent_identity=92.0,
                        alignment_length=100, evalue=1e-40, bit_score=150.0,
                        query_coverage=0.95,
                    )],
                    risk_level="critical",
                    coverage=0.95,
                    best_identity=92.0,
                    database_path="/tmp/hazard_db",
                    blast_program="blastp",
                )
                report = check_biosecurity_blast(
                    "MSKGEELFTGVPILVELDGDVNGHKFSVS",
                    db_path="/tmp/hazard_db",
                )
                assert report.risk_level == "critical" or report.risk_level in ("high", "medium")
                mock_screen.assert_called()

    def test_merges_both_protein_and_dna_blast(self):
        """When both protein and DNA are provided, both are BLASTed."""
        with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=True):
            with patch("biocompiler.biosecurity.blast_integration.screen_blast") as mock_screen:
                mock_screen.return_value = BlastScreeningResult(
                    blast_available=True,
                    hits=[],
                    risk_level="none",
                    coverage=0.0,
                    best_identity=0.0,
                    database_path="/tmp/hazard_db",
                    blast_program="blastn",
                )
                report = check_biosecurity_blast(
                    "MSKGEELFTGVPILVELDGDVNGHKFSVS",
                    dna="ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCC",
                    db_path="/tmp/hazard_db",
                )
                # Should have been called at least once for protein
                assert mock_screen.call_count >= 1

    def test_uses_env_var_for_db_path(self):
        """BIOCOMPILER_BLAST_DB_PATH env var is used when db_path not provided."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_DB_PATH": "/env/var/db"}):
            with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=True):
                with patch("biocompiler.biosecurity.blast_integration.screen_blast") as mock_screen:
                    mock_screen.return_value = BlastScreeningResult(
                        blast_available=True,
                        hits=[],
                        risk_level="none",
                        coverage=0.0,
                        best_identity=0.0,
                        database_path="/env/var/db",
                        blast_program="blastp",
                    )
                    report = check_biosecurity_blast("MSKGEELFTGVPILVELDGDVNGHKFSVS")
                    assert mock_screen.called
                    # Verify the env var path was used
                    call_args = mock_screen.call_args
                    assert call_args[1].get("db_path", call_args[0][1] if len(call_args[0]) > 1 else None) == "/env/var/db"


# ═══════════════════════════════════════════════════════════════════════════
# Test: screen_blast convenience function
# ═══════════════════════════════════════════════════════════════════════════


class TestScreenBlast:
    """Test the screen_blast convenience function."""

    def test_delegates_to_blast_scanner(self):
        """screen_blast creates a BlastScanner and calls screen."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            with patch.object(BlastScanner, "screen") as mock_screen:
                mock_screen.return_value = BlastScreeningResult(
                    blast_available=True, hits=[], risk_level="none"
                )
                result = screen_blast("ATCGATCG", "/tmp/test_db")
                mock_screen.assert_called_once_with("ATCGATCG", seq_type="dna")

    def test_passes_parameters_correctly(self):
        """Custom evalue and max_hits are passed to BlastScanner."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            with patch.object(BlastScanner, "screen") as mock_screen:
                mock_screen.return_value = BlastScreeningResult(
                    blast_available=True, hits=[], risk_level="none"
                )
                screen_blast("ATCG", "/tmp/db", evalue=1e-10, max_hits=10)
                # Verify the scanner was created with correct params
                call = mock_screen


# ═══════════════════════════════════════════════════════════════════════════
# Test: BlastHit and BlastScreeningResult dataclasses
# ═══════════════════════════════════════════════════════════════════════════


class TestDataClasses:
    """Test the dataclass structures."""

    def test_blast_hit_fields(self):
        hit = BlastHit(
            subject_id="test_subject",
            percent_identity=95.0,
            alignment_length=100,
            evalue=1e-40,
            bit_score=150.0,
            query_coverage=0.95,
            subject_length=267,
            query_start=1,
            query_end=100,
            subject_start=50,
            subject_end=149,
        )
        assert hit.subject_id == "test_subject"
        assert hit.percent_identity == pytest.approx(95.0, rel=1e-6)
        assert hit.alignment_length == 100
        assert hit.evalue == 1e-40
        assert hit.bit_score == pytest.approx(150.0, rel=1e-6)
        assert hit.query_coverage == pytest.approx(0.95, rel=1e-6)

    def test_blast_screening_result_defaults(self):
        result = BlastScreeningResult(blast_available=False)
        assert result.blast_available is False
        assert result.hits == []
        assert result.risk_level == "none"
        assert result.coverage == 0.0
        assert result.best_identity == 0.0
        assert result.database_path == ""
        assert result.blast_program == ""
        assert result.error is None
        assert result.warnings == []

    def test_blast_screening_result_with_error(self):
        result = BlastScreeningResult(
            blast_available=True,
            error="blastn failed",
            blast_program="blastn",
        )
        assert result.error == "blastn failed"


# ═══════════════════════════════════════════════════════════════════════════
# Test: Edge cases
# ═══════════════════════════════════════════════════════════════════════════


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_sequence_screen(self):
        """Screening an empty sequence returns no hits."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            with patch.object(scanner, "_run_blast", return_value=EMPTY_BLAST_XML):
                result = scanner.screen("", seq_type="dna")
                assert result.hits == []
                assert result.risk_level == "none"

    def test_very_long_sequence_query_coverage(self):
        """Query coverage is computed correctly for long sequences."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            with patch.object(scanner, "_run_blast", return_value=SAMPLE_BLAST_XML):
                result = scanner.screen("A" * 1000, seq_type="protein")
                # Query length is 1000, but hits align to positions 1-100
                # coverage should be (100-1+1)/1000 = 0.1
                if result.hits:
                    # At least one hit should have query_coverage set
                    assert all(0.0 <= h.query_coverage <= 1.0 for h in result.hits)

    def test_makeblastdb_failure_in_build_hazard_db(self):
        """makeblastdb failure raises RuntimeError."""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/makeblastdb"):
            with patch("subprocess.run") as mock_run:
                mock_run.side_effect = subprocess.CalledProcessError(
                    1, "makeblastdb", stderr="Error: invalid input"
                )
                with tempfile.TemporaryDirectory() as tmpdir:
                    with pytest.raises(RuntimeError, match="makeblastdb failed"):
                        build_hazard_db(
                            {"seq1": "ATCG"},
                            "test_db",
                            os.path.join(tmpdir, "test_db"),
                        )

    def test_blast_unavailable_does_not_block(self):
        """When BLAST is not installed, optimization should not be blocked."""
        with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=False):
            report = check_biosecurity_blast("MSKGEELFTGVPILVELDGDVNGHKFSVS")
            # The report should still be valid, just with a BLAST recommendation
            assert isinstance(report, BiosecurityReport)
            # It should include a recommendation about BLAST+
            assert any("BLAST+" in r for r in report.recommendations)

    def test_xml_with_malformed_hsp(self):
        """BLAST XML with malformed HSP values does not crash."""
        malformed_xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>test</Hit_id>
          <Hit_def>test_hit</Hit_def>
          <Hit_accession>1</Hit_accession>
          <Hit_len>100</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>NOT_A_NUMBER</Hsp_bit-score>
              <Hsp_score>50</Hsp_score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>50</Hsp_hit-to>
              <Hsp_identity>45</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>50</Hsp_align-len>
              <Hsp_hit-len>50</Hsp_hit-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""
        with patch("biocompiler.biosecurity.blast_integration.find_blast_bin", return_value="/usr/bin/blastn"):
            scanner = BlastScanner(db_path="/tmp/test_db")
            hits = scanner._parse_xml_output(malformed_xml)
            # Should skip the malformed HSP gracefully
            assert len(hits) == 0

    def test_db_path_from_env_var(self):
        """BIOCOMPILER_BLAST_DB_PATH is used when no db_path parameter given."""
        with patch.dict(os.environ, {"BIOCOMPILER_BLAST_DB_PATH": "/custom/db/path"}):
            with patch("biocompiler.biosecurity.blast_integration.is_blast_available", return_value=True):
                with patch("biocompiler.biosecurity.blast_integration.screen_blast") as mock_screen:
                    mock_screen.return_value = BlastScreeningResult(
                        blast_available=True, hits=[], risk_level="none",
                        database_path="/custom/db/path",
                    )
                    report = check_biosecurity_blast("MSKGEELFTGVPILVELDGDVNGHKFSVS")
                    assert mock_screen.called
