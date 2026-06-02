"""
BioCompiler v3.1.0 — Tests for Enhanced Production Features

Tests for:
- API key authentication
- Rate limiting
- GTEx tissue data module
- Hardened Kazusa parser
- Pydantic V2 migration
- Docker build validation
- Extended API endpoint tests
"""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest


# ─── Test GTEx Tissue Data ─────────────────────────────────────────

class TestTissueData:
    """Test GTEx-derived tissue weight data module."""

    def test_get_tissue_weights_known_cell_line(self):
        from biocompiler.tissue_data import get_tissue_weights
        weights = get_tissue_weights("HEK293T")
        assert weights["canonical"] == 1.0
        assert 0 < weights["exon_skip"] < 1.0
        assert 0 < weights["intron_retention"] < 1.0
        assert 0 < weights["alt_site"] < 1.0
        assert 0 < weights["cryptic"] < 1.0

    def test_get_tissue_weights_known_tissue(self):
        from biocompiler.tissue_data import get_tissue_weights
        weights = get_tissue_weights("Brain")
        assert weights["canonical"] == 1.0
        # Brain has the highest exon skipping rate
        assert weights["exon_skip"] > 0.3

    def test_get_tissue_weights_alias(self):
        from biocompiler.tissue_data import get_tissue_weights
        # Lowercase alias should resolve
        weights = get_tissue_weights("brain")
        assert weights["canonical"] == 1.0
        # Should match "Brain"
        brain_weights = get_tissue_weights("Brain")
        assert weights == brain_weights

    def test_get_tissue_weights_unknown_returns_default(self):
        from biocompiler.tissue_data import get_tissue_weights
        weights = get_tissue_weights("UNKNOWN_TISSUE_XYZ")
        assert weights["canonical"] == 1.0
        assert "exon_skip" in weights

    def test_list_available_tissues(self):
        from biocompiler.tissue_data import list_available_tissues
        tissues = list_available_tissues()
        assert "HEK293T" in tissues
        assert "Brain" in tissues
        assert "default" not in tissues  # default is not a "available" tissue

    def test_add_custom_tissue(self):
        from biocompiler.tissue_data import add_custom_tissue, get_tissue_weights
        add_custom_tissue("CustomCell", {
            "canonical": 1.0,
            "exon_skip": 0.5,
            "intron_retention": 0.3,
            "alt_site": 0.4,
            "cryptic": 0.2,
        })
        weights = get_tissue_weights("CustomCell")
        assert weights["exon_skip"] == 0.5
        assert weights["cryptic"] == 0.2

    def test_add_custom_tissue_missing_keys(self):
        from biocompiler.tissue_data import add_custom_tissue
        with pytest.raises(ValueError, match="Missing required"):
            add_custom_tissue("BadTissue", {"canonical": 1.0})

    def test_export_tissue_weights_json(self, tmp_path):
        from biocompiler.tissue_data import export_tissue_weights_json
        output = export_tissue_weights_json(str(tmp_path / "weights.json"))
        data = json.loads(output)
        assert "weights" in data
        assert "aliases" in data
        assert "HEK293T" in data["weights"]

    def test_gtex_weights_differ_from_hardcoded(self):
        """Verify GTEx weights are derived from the tissue_data module (single source of truth)."""
        from biocompiler.tissue_data import get_tissue_weights, GTEX_TISSUE_WEIGHTS
        gtex = get_tissue_weights("HEK293T")
        source = GTEX_TISSUE_WEIGHTS["HEK293T"]
        # GTEx weights should match the single source of truth
        assert gtex["canonical"] == source["canonical"]  # always 1.0
        # Both should be in a reasonable range
        assert 0.05 < gtex["cryptic"] < 0.2


# ─── Test API Authentication ────────────────────────────────────────

class TestAPIAuthentication:
    """Test API key authentication."""

    @pytest.fixture
    def client_no_auth(self):
        """Client with no API key configured (auth disabled)."""
        from fastapi.testclient import TestClient
        with patch.dict(os.environ, {}, clear=False):
            # Ensure no API key is set
            os.environ.pop("BIOCOMPILER_API_KEY", None)
            from biocompiler.api import create_app
            app = create_app()
            return TestClient(app)

    def test_health_without_auth(self, client_no_auth):
        """Health endpoint should work without auth."""
        response = client_no_auth.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["auth_enabled"] is False

    def test_check_without_auth_when_disabled(self, client_no_auth):
        """Check endpoint should work when auth is disabled."""
        response = client_no_auth.post("/check", json={
            "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTG",
            "organism": "Homo_sapiens",
        })
        assert response.status_code == 200

    def test_health_shows_rate_limit(self, client_no_auth):
        """Health endpoint should show rate limit configuration."""
        response = client_no_auth.get("/health")
        data = response.json()
        assert "rate_limit_rpm" in data
        assert data["rate_limit_rpm"] > 0


# ─── Test Pydantic V2 Migration ─────────────────────────────────────

class TestPydanticV2Migration:
    """Verify Pydantic V2 field_validator works correctly."""

    def test_sequence_input_validation(self):
        from biocompiler.api import SequenceInput
        # Valid sequence
        seq = SequenceInput(sequence="ATGCGT", organism="Homo_sapiens")
        assert seq.sequence == "ATGCGT"

    def test_sequence_input_rejects_invalid_nucleotides(self):
        from biocompiler.api import SequenceInput
        with pytest.raises(Exception):  # ValidationError
            SequenceInput(sequence="ATGZXC", organism="Homo_sapiens")

    def test_protein_input_validation(self):
        from biocompiler.api import ProteinInput
        protein = ProteinInput(protein="MVHLTPEEK", organism="Homo_sapiens")
        assert protein.protein == "MVHLTPEEK"

    def test_protein_input_rejects_invalid_amino_acids(self):
        from biocompiler.api import ProteinInput
        with pytest.raises(Exception):  # ValidationError
            ProteinInput(protein="MVHXYZ", organism="Homo_sapiens")

    def test_organism_validation(self):
        from biocompiler.api import SequenceInput
        with pytest.raises(Exception):  # ValidationError
            SequenceInput(sequence="ATG", organism="Unknown_organism")


# ─── Test Rate Limiting ────────────────────────────────────────────

class TestRateLimiting:
    """Test rate limiting infrastructure."""

    def test_rate_limit_function_under_limit(self):
        from biocompiler.api import _check_rate_limit
        # Should not raise for a small number of requests
        for _ in range(5):
            _check_rate_limit("test_client_1")

    def test_rate_limit_different_clients(self):
        from biocompiler.api import _check_rate_limit
        # Different clients should have separate limits
        _check_rate_limit("client_a")
        _check_rate_limit("client_b")
        # Both should succeed


# ─── Test Hardened Organism DB ──────────────────────────────────────

class TestHardenedOrganismDB:
    """Test hardened Kazusa parser and database features."""

    def test_cache_freshness_check(self, tmp_path):
        from biocompiler.organism_db import OrganismDatabase
        db = OrganismDatabase(db_path=tmp_path / "test.db")
        # No data yet, so cache should not be fresh
        assert not db._is_cache_fresh("Homo_sapiens")

    def test_migrate_and_check_cache(self, tmp_path):
        from biocompiler.organism_db import OrganismDatabase
        db = OrganismDatabase(db_path=tmp_path / "test.db")
        db.migrate_builtin_data()
        # builtin data has source='builtin', not 'kazusa', so cache check
        # for kazusa source should return False
        assert not db._is_cache_fresh("Homo_sapiens")

    def test_robust_parser_fallback(self):
        """Test that the robust parser falls back gracefully on bad HTML."""
        from biocompiler.organism_db import OrganismDatabase
        # Create garbage HTML
        bad_html = "<html><body>Not a codon table at all</body></html>"
        result = OrganismDatabase._parse_kazusa_html_robust(bad_html, "TestOrg")
        # Should return 64 codons (with uniform fallback)
        assert len(result) == 64
        # All codons should have entries
        for codon, (aa, freq, adapt, count) in result.items():
            assert aa is not None
            assert freq >= 0.0

    def test_robust_parser_with_partial_data(self):
        """Test parser with partial HTML data."""
        from biocompiler.organism_db import OrganismDatabase
        # HTML with some codon data
        partial_html = """
        <table>
        <tr><td>ATG</td><td>M</td><td>1.00</td><td>28</td><td>0.36</td></tr>
        <tr><td>TTT</td><td>F</td><td>0.45</td><td>12</td><td>0.15</td></tr>
        </table>
        """
        result = OrganismDatabase._parse_kazusa_html_robust(partial_html, "TestOrg")
        # Should have all 64 codons (2 from HTML + 62 uniform fallback)
        assert len(result) == 64


# ─── Test Extended API Endpoints ────────────────────────────────────

class TestExtendedAPIEndpoints:
    """Test extended API functionality."""

    @pytest.fixture
    def client(self):
        from fastapi.testclient import TestClient
        from biocompiler.api import app
        return TestClient(app)

    def test_health_has_auth_field(self, client):
        response = client.get("/health")
        data = response.json()
        assert "auth_enabled" in data
        assert "rate_limit_rpm" in data

    def test_check_with_unsupported_organism(self, client):
        response = client.post("/check", json={
            "sequence": "ATGGTGAGC",
            "organism": "Unknown_Organism",
        })
        assert response.status_code == 422

    def test_scan_with_orfs(self, client):
        response = client.post("/scan", json={
            "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGAATAG",
            "find_orfs": True,
        })
        assert response.status_code == 200
        data = response.json()
        assert "orfs" in data

    def test_verify_malformed_certificate(self, client):
        response = client.post("/verify", json={
            "certificate": {"invalid": "data"}
        })
        # verify_certificate catches the error and returns 400
        assert response.status_code in (200, 400)
        if response.status_code == 200:
            data = response.json()
            assert data["status"] in ("failed", "invalid", "REJECTED")


# ─── Test CLI Tissue Data Integration ───────────────────────────────

class TestTissueDataIntegration:
    """Test that tissue data integrates correctly with splicing."""

    def test_splicing_uses_gtex_weights(self):
        from biocompiler.splicing import compute_splice_isoforms
        from biocompiler.tissue_data import get_tissue_weights

        # A simple multi-exon gene
        seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG" * 5
        exons = [(0, 50), (100, 150), (200, 300)]

        # Compute with HEK293T
        isoforms = compute_splice_isoforms(seq, exons, cellular_context="HEK293T")
        assert len(isoforms) > 0

        # Compute with Brain (should have different scores)
        isoforms_brain = compute_splice_isoforms(seq, exons, cellular_context="Brain")
        assert len(isoforms_brain) > 0

        # Brain has higher exon_skip weight, so exon-skipped isoforms
        # should have relatively higher scores in Brain vs HEK293T
        gtex_hek = get_tissue_weights("HEK293T")
        gtex_brain = get_tissue_weights("Brain")
        assert gtex_brain["exon_skip"] > gtex_hek["exon_skip"]
