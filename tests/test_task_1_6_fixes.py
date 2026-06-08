"""
Tests for Task 1.6: Fix CLI check command, export validation, and API endpoint gaps.

Covers:
1. CLI check command uses evaluate_all_predicates() instead of hardcoded 8 calls
2. CLI check accepts --organism-domain, --source-organism, --therapeutic, --strict-mode, --json
3. CLI CSP strategy produces predicate results and certificate
4. Export endpoint DNA character validation (ExportFastaInput, ExportGenbankInput, ScanInput)
5. Organism validation on export endpoints
6. locus_name length validation for GenBank (max 16 chars)
7. Pydantic response models for /enzymes, /provenance/{id}, /provenance, /export/fasta,
   /export/genbank, /validate-datasets
8. SBOL3 export API endpoint POST /export/sbol3
9. --format sbol3 CLI option for optimize command
10. HTTP status code consistency (InvalidSequenceError → 400, provenance corruption → 422)
"""

from __future__ import annotations

import argparse
import os
import tempfile
from io import StringIO
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from biocompiler import __version__
from biocompiler.cli import (
    build_parser,
    cmd_check,
    cmd_optimize,
    main,
)
from biocompiler.api import (
    create_app,
    _rate_limiter,
    EnzymeListResponse,
    ProvenanceDetailResponse,
    ProvenanceListResponse,
    ProvenanceRecordSummary,
    ExportFastaResponse,
    ExportGenbankResponse,
    ExportSbol3Input,
    ExportSbol3Response,
    DatasetValidationResponse,
    DatasetValidationResult,
    ExportFastaInput,
    ExportGenbankInput,
    ScanInput,
)


# ─── Fixtures ────────────────────────────────────────────────────────

@pytest.fixture()
def app():
    """Create a fresh FastAPI app for each test."""
    import biocompiler.api as _api_module
    _rate_limiter.clear()
    original_mode = _api_module._AUTH_MODE
    original_keys = set(_api_module._CONFIGURED_API_KEYS)
    _api_module._AUTH_MODE = "disabled"
    _api_module._CONFIGURED_API_KEYS = set()
    try:
        yield create_app()
    finally:
        _api_module._AUTH_MODE = original_mode
        _api_module._CONFIGURED_API_KEYS = original_keys


@pytest.fixture()
def client(app):
    """Create a TestClient from the app."""
    return TestClient(app)


VALID_DNA = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
VALID_PROTEIN = "MGSWKRQPPA"
INVALID_DNA = "ATGGCCXYZ"


# ═══════════════════════════════════════════════════════════════════════
# 1. CLI check command uses evaluate_all_predicates()
# ═══════════════════════════════════════════════════════════════════════

class TestCLICheckUsesRegistry:
    """Verify that cmd_check uses evaluate_all_predicates() from the registry."""

    def test_check_produces_more_than_8_results(self, tmp_path):
        """evaluate_all_predicates returns 12 results, not just 8."""
        fasta = tmp_path / "gene.fasta"
        fasta.write_text(">test\nATGGCTAAGCTGGATCC\n")
        args = argparse.Namespace(
            input=str(fasta),
            organism=None,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            organism_domain="auto",
            source_organism=None,
            therapeutic=False,
            strict_mode=False,
            json=True,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        import json
        output = json.loads(mock_out.getvalue())
        assert "predicate_results" in output
        # evaluate_all_predicates returns 12 predicates, not 8
        assert len(output["predicate_results"]) > 8

    def test_check_json_output_format(self, tmp_path):
        """cmd_check --json should produce valid JSON with expected fields."""
        fasta = tmp_path / "gene.fasta"
        fasta.write_text(">test\nATGGCTAAGCTGGATCC\n")
        args = argparse.Namespace(
            input=str(fasta),
            organism=None,
            species="human",
            enzymes="",
            splice_low=3.0,
            splice_high=6.0,
            organism_domain="auto",
            source_organism=None,
            therapeutic=False,
            strict_mode=False,
            json=True,
        )
        with patch("sys.stdout", new_callable=StringIO) as mock_out:
            cmd_check(args)
        import json
        output = json.loads(mock_out.getvalue())
        assert "version" in output
        assert "organism" in output
        assert "organism_domain" in output
        assert "certificate_level" in output
        assert "predicate_results" in output
        assert "sequence_length" in output
        assert "gc_content" in output


# ═══════════════════════════════════════════════════════════════════════
# 2. CLI check new options
# ═══════════════════════════════════════════════════════════════════════

class TestCLICheckNewOptions:
    """Verify that check subparser accepts new CLI options."""

    def test_check_accepts_organism_domain(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--organism-domain", "prokaryote"])
        assert args.organism_domain == "prokaryote"

    def test_check_accepts_source_organism(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--source-organism", "ecoli"])
        assert args.source_organism == "ecoli"

    def test_check_accepts_therapeutic(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--therapeutic"])
        assert args.therapeutic is True

    def test_check_accepts_strict_mode(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--strict-mode"])
        assert args.strict_mode is True

    def test_check_accepts_json(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--json"])
        assert args.json is True

    def test_check_accepts_organism(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta", "--organism", "ecoli"])
        assert args.organism == "ecoli"

    def test_check_default_organism_domain(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta"])
        assert args.organism_domain == "auto"

    def test_check_default_therapeutic(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta"])
        assert args.therapeutic is False

    def test_check_default_strict_mode(self):
        parser = build_parser()
        args = parser.parse_args(["check", "gene.fasta"])
        assert args.strict_mode is False


# ═══════════════════════════════════════════════════════════════════════
# 3. CLI optimize --format sbol3 option
# ═══════════════════════════════════════════════════════════════════════

class TestCLIOptimizeFormatOption:
    """Verify that optimize accepts --format sbol3."""

    def test_optimize_accepts_format_fasta(self):
        parser = build_parser()
        args = parser.parse_args(["optimize", "MSKGEELFTG", "--format", "fasta"])
        assert args.format == "fasta"

    def test_optimize_accepts_format_genbank(self):
        parser = build_parser()
        args = parser.parse_args(["optimize", "MSKGEELFTG", "--format", "genbank"])
        assert args.format == "genbank"

    def test_optimize_accepts_format_sbol3(self):
        parser = build_parser()
        args = parser.parse_args(["optimize", "MSKGEELFTG", "--format", "sbol3"])
        assert args.format == "sbol3"

    def test_optimize_format_default(self):
        parser = build_parser()
        args = parser.parse_args(["optimize", "MSKGEELFTG"])
        assert args.format == "fasta"

    def test_optimize_format_invalid_rejected(self):
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["optimize", "MSKGEELFTG", "--format", "pdf"])


# ═══════════════════════════════════════════════════════════════════════
# 4. Export endpoint DNA character validation
# ═══════════════════════════════════════════════════════════════════════

class TestExportFastaInputValidation:
    """Test ExportFastaInput DNA character validation."""

    def test_valid_dna_sequence(self):
        m = ExportFastaInput(sequence="ATGCGTACGT")
        assert m.sequence == "ATGCGTACGT"

    def test_invalid_dna_characters_rejected(self):
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            ExportFastaInput(sequence="ATGXYZ")

    def test_empty_sequence_rejected(self):
        with pytest.raises(ValueError, match="must not be empty"):
            ExportFastaInput(sequence="")

    def test_lowercase_normalized(self):
        m = ExportFastaInput(sequence="atgcgt")
        assert m.sequence == "ATGCGT"

    def test_unsupported_organism_rejected(self):
        with pytest.raises(ValueError, match="Unsupported organism"):
            ExportFastaInput(sequence=VALID_DNA, organism="Alien_xenomorph")


class TestExportGenbankInputValidation:
    """Test ExportGenbankInput DNA character and locus_name validation."""

    def test_valid_dna_sequence(self):
        m = ExportGenbankInput(sequence="ATGCGTACGT")
        assert m.sequence == "ATGCGTACGT"

    def test_invalid_dna_characters_rejected(self):
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            ExportGenbankInput(sequence="ATGXYZ123")

    def test_empty_sequence_rejected(self):
        with pytest.raises(ValueError, match="must not be empty"):
            ExportGenbankInput(sequence="")

    def test_locus_name_max_16_chars(self):
        m = ExportGenbankInput(sequence="ATGCGT", locus_name="ABC123")
        assert m.locus_name == "ABC123"

    def test_locus_name_too_long_rejected(self):
        with pytest.raises(ValueError, match="locus_name must be at most 16 characters"):
            ExportGenbankInput(sequence="ATGCGT", locus_name="A" * 17)

    def test_locus_name_exactly_16_chars(self):
        m = ExportGenbankInput(sequence="ATGCGT", locus_name="A" * 16)
        assert len(m.locus_name) == 16

    def test_unsupported_organism_rejected(self):
        with pytest.raises(ValueError, match="Unsupported organism"):
            ExportGenbankInput(sequence=VALID_DNA, organism="Fictional_org")


class TestScanInputValidation:
    """Test ScanInput DNA character validation."""

    def test_valid_dna_sequence(self):
        m = ScanInput(sequence="ATGCGTACGT")
        assert m.sequence == "ATGCGTACGT"

    def test_invalid_dna_characters_rejected(self):
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            ScanInput(sequence="ATGXYZ")

    def test_empty_sequence_rejected(self):
        with pytest.raises(ValueError, match="must not be empty"):
            ScanInput(sequence="")


# ═══════════════════════════════════════════════════════════════════════
# 5. Pydantic response models
# ═══════════════════════════════════════════════════════════════════════

class TestNewResponseModels:
    """Test new Pydantic response models."""

    def test_enzyme_list_response(self):
        m = EnzymeListResponse(enzymes={"EcoRI": "GAATTC"})
        assert m.enzymes["EcoRI"] == "GAATTC"

    def test_provenance_detail_response(self):
        m = ProvenanceDetailResponse(id="test-id", trail={"key": "value"})
        assert m.id == "test-id"
        assert m.trail == {"key": "value"}

    def test_provenance_record_summary(self):
        m = ProvenanceRecordSummary(gene_name="insulin", organism="Homo_sapiens")
        assert m.gene_name == "insulin"

    def test_provenance_list_response(self):
        m = ProvenanceListResponse(count=0, records=[])
        assert m.count == 0

    def test_export_fasta_response(self):
        m = ExportFastaResponse(format="fasta", content=">test\nATG")
        assert m.format == "fasta"
        assert "ATG" in m.content

    def test_export_genbank_response(self):
        m = ExportGenbankResponse(format="genbank", content="LOCUS BIOCOMPILER")
        assert m.format == "genbank"

    def test_dataset_validation_result(self):
        m = DatasetValidationResult(dataset="human", test_name="gc_content", passed=True)
        assert m.passed is True

    def test_dataset_validation_response(self):
        m = DatasetValidationResponse(total_tests=5, passed=4, failed=1, results=[])
        assert m.total_tests == 5


# ═══════════════════════════════════════════════════════════════════════
# 6. API endpoint response models
# ═══════════════════════════════════════════════════════════════════════

class TestAPIEndpointResponseModels:
    """Verify endpoints use the new response models."""

    def test_enzymes_endpoint_returns_enzyme_list(self, client):
        resp = client.get("/enzymes")
        assert resp.status_code == 200
        data = resp.json()
        assert "enzymes" in data
        assert isinstance(data["enzymes"], dict)

    def test_export_fasta_returns_typed_response(self, client):
        resp = client.post("/export/fasta", json={"sequence": VALID_DNA})
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "fasta"
        assert "content" in data

    def test_export_genbank_returns_typed_response(self, client):
        resp = client.post("/export/genbank", json={"sequence": VALID_DNA})
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "genbank"
        assert "content" in data


# ═══════════════════════════════════════════════════════════════════════
# 7. SBOL3 export API endpoint
# ═══════════════════════════════════════════════════════════════════════

class TestSBOL3ExportEndpoint:
    """Test POST /export/sbol3 endpoint."""

    def test_sbol3_endpoint_exists(self, client):
        resp = client.post("/export/sbol3", json={"sequence": VALID_DNA})
        # Should not be 404
        assert resp.status_code != 404

    def test_sbol3_export_valid_sequence(self, client):
        resp = client.post("/export/sbol3", json={"sequence": VALID_DNA})
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "sbol3"
        assert "content" in data
        # SBOL3 XML should contain namespace references
        assert "sbols.org" in data["content"] or "rdf" in data["content"].lower()

    def test_sbol3_export_invalid_sequence(self, client):
        resp = client.post("/export/sbol3", json={"sequence": "ATGXYZ"})
        assert resp.status_code == 422

    def test_sbol3_export_empty_sequence(self, client):
        resp = client.post("/export/sbol3", json={"sequence": ""})
        assert resp.status_code == 422

    def test_sbol3_export_unsupported_organism(self, client):
        resp = client.post("/export/sbol3", json={
            "sequence": VALID_DNA,
            "organism": "Alien_xenomorph",
        })
        assert resp.status_code == 422

    def test_sbol3_export_json_format(self, client):
        resp = client.post("/export/sbol3", json={
            "sequence": VALID_DNA,
            "format": "sbol3json",
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "sbol3json"

    def test_sbol3_export_invalid_format(self, client):
        resp = client.post("/export/sbol3", json={
            "sequence": VALID_DNA,
            "format": "invalid",
        })
        assert resp.status_code == 422


class TestSBOL3InputModel:
    """Test ExportSbol3Input Pydantic model."""

    def test_valid_sbol3_input(self):
        m = ExportSbol3Input(sequence=VALID_DNA)
        assert m.sequence == VALID_DNA

    def test_invalid_nucleotides(self):
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            ExportSbol3Input(sequence="ATGXYZ")

    def test_empty_sequence(self):
        with pytest.raises(ValueError, match="must not be empty"):
            ExportSbol3Input(sequence="")

    def test_invalid_format(self):
        with pytest.raises(ValueError, match="Unsupported SBOL format"):
            ExportSbol3Input(sequence=VALID_DNA, format="invalid")

    def test_default_values(self):
        m = ExportSbol3Input(sequence=VALID_DNA)
        assert m.organism == "Homo_sapiens"
        assert m.gene_name == "optimized_gene"
        assert m.format == "sbol3"


# ═══════════════════════════════════════════════════════════════════════
# 8. HTTP status code consistency
# ═══════════════════════════════════════════════════════════════════════

class TestHTTPStatusCodes:
    """Test that HTTP status codes are consistent and correct."""

    def test_invalid_dna_export_fasta_returns_400_or_422(self, client):
        """InvalidSequenceError should return 400, not 500."""
        resp = client.post("/export/fasta", json={"sequence": INVALID_DNA})
        # Should be 422 from Pydantic validation (since we validate at model level)
        assert resp.status_code in (400, 422)

    def test_invalid_dna_export_genbank_returns_400_or_422(self, client):
        resp = client.post("/export/genbank", json={"sequence": INVALID_DNA})
        assert resp.status_code in (400, 422)

    def test_invalid_dna_scan_returns_400_or_422(self, client):
        resp = client.post("/scan", json={"sequence": INVALID_DNA})
        assert resp.status_code in (400, 422)

    def test_provenance_not_found_returns_404(self, client):
        resp = client.get("/provenance/nonexistent-id")
        # Provenance store may validate UUID format first (422) or not found (404)
        assert resp.status_code in (404, 422)

    def test_export_fasta_missing_sequence_returns_422(self, client):
        resp = client.post("/export/fasta", json={})
        assert resp.status_code == 422

    def test_export_genbank_long_locus_name_returns_422(self, client):
        resp = client.post("/export/genbank", json={
            "sequence": VALID_DNA,
            "locus_name": "A" * 17,
        })
        assert resp.status_code == 422

    def test_check_invalid_nucleotides_returns_422(self, client):
        resp = client.post("/check", json={"sequence": INVALID_DNA})
        assert resp.status_code == 422

    def test_scan_invalid_nucleotides_returns_422(self, client):
        resp = client.post("/scan", json={"sequence": INVALID_DNA})
        assert resp.status_code == 422


# ═══════════════════════════════════════════════════════════════════════
# 9. Integration: Export endpoints with validated inputs
# ═══════════════════════════════════════════════════════════════════════

class TestExportEndpointValidationIntegration:
    """Integration tests for export endpoint validation."""

    def test_export_fasta_with_valid_organism(self, client):
        resp = client.post("/export/fasta", json={
            "sequence": VALID_DNA,
            "organism": "Homo_sapiens",
        })
        assert resp.status_code == 200

    def test_export_fasta_with_invalid_organism(self, client):
        resp = client.post("/export/fasta", json={
            "sequence": VALID_DNA,
            "organism": "Fake_organism",
        })
        assert resp.status_code == 422

    def test_export_genbank_with_valid_organism(self, client):
        resp = client.post("/export/genbank", json={
            "sequence": VALID_DNA,
            "organism": "Escherichia_coli",
        })
        assert resp.status_code == 200

    def test_export_genbank_with_invalid_organism(self, client):
        resp = client.post("/export/genbank", json={
            "sequence": VALID_DNA,
            "organism": "Nonexistent_species",
        })
        assert resp.status_code == 422

    def test_export_genbank_locus_name_16_chars_ok(self, client):
        resp = client.post("/export/genbank", json={
            "sequence": VALID_DNA,
            "locus_name": "A" * 16,
        })
        assert resp.status_code == 200

    def test_export_genbank_locus_name_17_chars_fails(self, client):
        resp = client.post("/export/genbank", json={
            "sequence": VALID_DNA,
            "locus_name": "A" * 17,
        })
        assert resp.status_code == 422

    def test_sbol3_route_in_openapi(self, client):
        """SBOL3 endpoint should appear in the OpenAPI schema."""
        resp = client.get("/openapi.json")
        schema = resp.json()
        assert "/export/sbol3" in schema["paths"]
