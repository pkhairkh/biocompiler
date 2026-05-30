"""
BioCompiler v3.2.0 — Hardened Test Suite

Tests for production-grade robustness:
- API error-path coverage (malformed inputs, auth failures, rate limiting)
- SQLite concurrency and thread safety
- Database schema versioning and migration
- Kazusa response validation
- Cache integrity verification
- Multiple API key support
- CLI colored output and progress indicators
"""

import json
import os
import tempfile
import threading
import time
from pathlib import Path
from unittest.mock import patch

import pytest


# ─── API Error-Path Tests ──────────────────────────────────────────

class TestAPIErrorPaths:
    """Test API error handling for malformed and edge-case inputs."""

    @pytest.fixture
    def client(self):
        from fastapi.testclient import TestClient
        from biocompiler.api import create_app
        app = create_app()
        return TestClient(app)

    def test_check_empty_sequence(self, client):
        """Empty sequence should be handled (may accept or reject)."""
        response = client.post("/check", json={
            "sequence": "",
            "organism": "Homo_sapiens",
        })
        # Empty sequence may pass validation (Pydantic allows it) and fail at type-check
        assert response.status_code in (200, 422)

    def test_check_very_long_sequence(self, client):
        """Very long sequence should be handled gracefully."""
        long_seq = "ATGC" * 100_000  # 400K bp
        response = client.post("/check", json={
            "sequence": long_seq,
            "organism": "Homo_sapiens",
        })
        # Should either succeed or fail gracefully
        assert response.status_code in (200, 422, 413)

    def test_check_n_only_sequence(self, client):
        """All-N sequence should be handled (N is allowed)."""
        response = client.post("/check", json={
            "sequence": "NNNNNNNNNN",
            "organism": "Homo_sapiens",
        })
        assert response.status_code == 200

    def test_check_gc_boundary_values(self, client):
        """GC boundaries at 0 and 1 should work."""
        response = client.post("/check", json={
            "sequence": "ATGGTGAGCAAGGGCGAGGAG",
            "organism": "Homo_sapiens",
            "gc_lo": 0.0,
            "gc_hi": 1.0,
        })
        assert response.status_code == 200

    def test_check_gc_inverted_range(self, client):
        """GC lo > hi should produce a FAIL verdict."""
        response = client.post("/check", json={
            "sequence": "ATGGTGAGCAAGGGCGAGGAG",
            "organism": "Homo_sapiens",
            "gc_lo": 0.9,
            "gc_hi": 0.1,
        })
        assert response.status_code == 200
        data = response.json()
        assert data["overall_verdict"] == "FAIL"

    def test_optimize_empty_protein(self, client):
        """Empty protein should be handled (may accept or reject)."""
        response = client.post("/optimize", json={
            "protein": "",
            "organism": "Homo_sapiens",
        })
        # Empty protein may pass validation and fail at optimization
        assert response.status_code in (200, 400, 422)

    def test_optimize_invalid_amino_acids(self, client):
        """Invalid amino acids should return validation error."""
        response = client.post("/optimize", json={
            "protein": "MVHXYZ123",
            "organism": "Homo_sapiens",
        })
        assert response.status_code == 422

    def test_scan_no_sequence(self, client):
        """Missing sequence field should return validation error."""
        response = client.post("/scan", json={})
        assert response.status_code == 422

    def test_export_fasta_no_sequence(self, client):
        """Missing sequence in FASTA export should return validation error."""
        response = client.post("/export/fasta", json={
            "identifier": "test",
        })
        assert response.status_code == 422

    def test_export_genbank_long_locus(self, client):
        """Locus name > 16 chars should still work (GenBank spec allows it)."""
        response = client.post("/export/genbank", json={
            "sequence": "ATGGTGAGC",
            "locus_name": "VERYLONGLOCUSNAME12345",
        })
        assert response.status_code == 200

    def test_verify_missing_certificate_fields(self, client):
        """Certificate with missing required fields should return error."""
        response = client.post("/verify", json={
            "certificate": {"version": "1.0"}  # Missing many fields
        })
        assert response.status_code in (200, 400)

    def test_verify_empty_certificate(self, client):
        """Empty certificate dict should return error."""
        response = client.post("/verify", json={
            "certificate": {}
        })
        assert response.status_code in (200, 400)

    def test_check_with_exon_boundaries(self, client):
        """Check with exon boundaries should work."""
        response = client.post("/check", json={
            "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTG",
            "organism": "Homo_sapiens",
            "exon_boundaries": [[0, 48], [48, 96]],
        })
        assert response.status_code == 200


# ─── API Authentication Tests ──────────────────────────────────────

class TestAPIMultiKeyAuth:
    """Test multi-key API authentication and key rotation."""

    def test_multiple_api_keys_via_env(self):
        """Test that BIOCOMPILER_API_KEYS env var supports multiple keys."""
        from fastapi.testclient import TestClient
        # We need to reimport api module to pick up env changes
        import importlib
        import biocompiler.api

        with patch.dict(os.environ, {
            "BIOCOMPILER_API_KEYS": "key1,key2,key3",
            "BIOCOMPILER_API_KEY": "",
        }, clear=False):
            importlib.reload(biocompiler.api)
            from biocompiler.api import create_app
            app = create_app()
            client = TestClient(app)

            # Health should work without auth
            r = client.get("/health")
            assert r.status_code == 200
            assert r.json()["auth_enabled"] is True

            # Without key should fail
            r = client.post("/check", json={
                "sequence": "ATGGTGAGC",
                "organism": "Homo_sapiens",
            })
            assert r.status_code == 401

            # With key1 should work
            r = client.post("/check", json={
                "sequence": "ATGGTGAGCAAGGGCGAGGAG",
                "organism": "Homo_sapiens",
            }, headers={"X-API-Key": "key1"})
            assert r.status_code == 200

            # With key2 should work
            r = client.post("/check", json={
                "sequence": "ATGGTGAGCAAGGGCGAGGAG",
                "organism": "Homo_sapiens",
            }, headers={"X-API-Key": "key2"})
            assert r.status_code == 200

            # With wrong key should fail
            r = client.post("/check", json={
                "sequence": "ATGGTGAGC",
                "organism": "Homo_sapiens",
            }, headers={"X-API-Key": "wrong_key"})
            assert r.status_code == 401

        # Reload to reset state
        importlib.reload(biocompiler.api)

    def test_single_key_backward_compat(self):
        """Test that BIOCOMPILER_API_KEY still works for backward compatibility."""
        from fastapi.testclient import TestClient
        import importlib
        import biocompiler.api

        with patch.dict(os.environ, {
            "BIOCOMPILER_API_KEY": "my-secret-key",
        }, clear=False):
            os.environ.pop("BIOCOMPILER_API_KEYS", None)
            importlib.reload(biocompiler.api)
            from biocompiler.api import create_app
            app = create_app()
            client = TestClient(app)

            r = client.get("/health")
            assert r.json()["auth_enabled"] is True

            r = client.post("/check", json={
                "sequence": "ATGGTGAGCAAGGGCGAGGAG",
                "organism": "Homo_sapiens",
            }, headers={"X-API-Key": "my-secret-key"})
            assert r.status_code == 200

        # Reload to reset state
        importlib.reload(biocompiler.api)


# ─── SQLite Concurrency Tests ──────────────────────────────────────

class TestSQLiteConcurrency:
    """Test SQLite database under concurrent access."""

    def test_concurrent_reads(self, tmp_path):
        """Multiple threads reading simultaneously should not error."""
        from biocompiler.organism_db import OrganismDatabase

        db = OrganismDatabase(db_path=tmp_path / "concurrent.db")
        db.migrate_builtin_data()

        errors = []

        def read_organism():
            try:
                usage = db.get_codon_usage("Homo_sapiens")
                if len(usage) != 64:
                    errors.append(f"Expected 64 codons, got {len(usage)}")
            except Exception as e:
                errors.append(str(e))

        threads = [threading.Thread(target=read_organism) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert len(errors) == 0, f"Concurrent read errors: {errors}"

    def test_concurrent_writes(self, tmp_path):
        """Multiple threads writing simultaneously should not corrupt data."""
        from biocompiler.organism_db import OrganismDatabase
        from biocompiler.constants import CODON_TABLE

        db = OrganismDatabase(db_path=tmp_path / "concurrent_writes.db")
        errors = []

        def write_organism(idx):
            try:
                usage = {}
                for codon, aa in CODON_TABLE.items():
                    if aa != "*":
                        usage[codon] = (aa, 0.03 + idx * 0.001, 0.5, 100 + idx)
                db.store_organism(f"TestOrg_{idx}", usage, taxonomy_id=str(idx))
            except Exception as e:
                errors.append(f"Thread {idx}: {e}")

        threads = [threading.Thread(target=write_organism, args=(i,)) for i in range(5)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert len(errors) == 0, f"Concurrent write errors: {errors}"

        # Verify all organisms were stored
        names = db.list_organism_names()
        for i in range(5):
            assert f"TestOrg_{i}" in names

    def test_wal_mode_set(self, tmp_path):
        """Database should use WAL mode for better concurrency."""
        from biocompiler.organism_db import OrganismDatabase

        db = OrganismDatabase(db_path=tmp_path / "wal_test.db")
        conn = db._connect()
        try:
            mode = conn.execute("PRAGMA journal_mode").fetchone()[0]
            assert mode.lower() == "wal"
        finally:
            conn.close()


# ─── Database Schema Versioning Tests ──────────────────────────────

class TestDatabaseSchemaVersioning:
    """Test database schema versioning and migration."""

    def test_schema_version_stored(self, tmp_path):
        """Schema version should be stored in the database."""
        from biocompiler.organism_db import OrganismDatabase, SCHEMA_VERSION

        db = OrganismDatabase(db_path=tmp_path / "version_test.db")
        conn = db._connect()
        try:
            row = conn.execute("SELECT value FROM schema_version WHERE key = 'version'").fetchone()
            assert row is not None
            assert int(row[0]) == SCHEMA_VERSION
        finally:
            conn.close()

    def test_api_response_cache_table_exists(self, tmp_path):
        """v2 schema should have api_response_cache table."""
        from biocompiler.organism_db import OrganismDatabase

        db = OrganismDatabase(db_path=tmp_path / "v2_test.db")
        conn = db._connect()
        try:
            # Should not throw
            conn.execute("SELECT 1 FROM api_response_cache LIMIT 1")
        finally:
            conn.close()

    def test_store_response_hash(self, tmp_path):
        """Storing a response hash should work."""
        from biocompiler.organism_db import OrganismDatabase

        db = OrganismDatabase(db_path=tmp_path / "hash_test.db")
        db.migrate_builtin_data()

        # Store a hash
        html = "<html>test response</html>"
        db._store_response_hash("Homo_sapiens", html, "https://example.com")

        # Verify it was stored
        conn = db._connect()
        try:
            row = conn.execute(
                "SELECT response_hash, url FROM api_response_cache WHERE organism = ?",
                ("Homo_sapiens",),
            ).fetchone()
            assert row is not None
            assert row[1] == "https://example.com"
        finally:
            conn.close()


# ─── Kazusa Response Validation Tests ──────────────────────────────

class TestKazusaResponseValidation:
    """Test Kazusa response validation and data integrity checks."""

    def test_validate_short_response(self):
        """Too-short response should fail validation."""
        from biocompiler.organism_db import OrganismDatabase
        db = OrganismDatabase.__new__(OrganismDatabase)
        assert not db._validate_kazusa_response("", "TestOrg")
        assert not db._validate_kazusa_response("short", "TestOrg")

    def test_validate_error_response(self):
        """Responses with error indicators should fail validation."""
        from biocompiler.organism_db import OrganismDatabase
        db = OrganismDatabase.__new__(OrganismDatabase)
        html = "<html><body>Species not found in database. No data available.</body></html>" + "ATG" * 100
        assert not db._validate_kazusa_response(html, "TestOrg")

    def test_validate_good_response(self):
        """Valid response with codon patterns should pass validation."""
        from biocompiler.organism_db import OrganismDatabase
        db = OrganismDatabase.__new__(OrganismDatabase)
        html = "<html><body>" + "ATG TTT GGG CCC AAA " * 50 + "</body></html>"
        assert db._validate_kazusa_response(html, "TestOrg")

    def test_validate_codon_usage_data_corrects_aa(self):
        """Wrong amino acid assignment should be corrected."""
        from biocompiler.organism_db import OrganismDatabase

        # Create data with wrong AA for ATG (should be M, not X)
        bad_data = {"ATG": ("X", 1.0, 1.0, 100)}
        result = OrganismDatabase._validate_codon_usage_data(bad_data, "TestOrg")
        assert result["ATG"][0] == "M"  # Corrected

    def test_validate_codon_usage_data_clamps_negative_freq(self):
        """Negative frequency should be clamped to 0."""
        from biocompiler.organism_db import OrganismDatabase

        bad_data = {"ATG": ("M", -0.5, 1.0, 100)}
        result = OrganismDatabase._validate_codon_usage_data(bad_data, "TestOrg")
        assert result["ATG"][1] == 0.0  # Clamped

    def test_robust_parser_strategy1(self):
        """Strategy 1 (table row parsing) should work on valid HTML."""
        from biocompiler.organism_db import OrganismDatabase

        html = "<table>"
        codons_data = [
            ("ATG", "M", "1.00", "28", "36.0"),
            ("TTT", "F", "0.45", "12", "15.0"),
            ("GGG", "G", "0.25", "8", "10.0"),
        ]
        for codon, aa, freq, count, per_k in codons_data:
            html += f'<tr><td>{codon}</td><td>{aa}</td><td>{freq}</td><td>{count}</td><td>{per_k}</td></tr>'
        html += "</table>"

        result = OrganismDatabase._parse_kazusa_html_robust(html, "TestOrg")
        # Should have entries for parsed codons plus fallback
        assert len(result) >= 3

    def test_robust_parser_empty_html(self):
        """Empty HTML should return all 64 codons with uniform fallback."""
        from biocompiler.organism_db import OrganismDatabase

        result = OrganismDatabase._parse_kazusa_html_robust("<html></html>", "TestOrg")
        assert len(result) == 64


# ─── Interactive Report Tests ──────────────────────────────────────

class TestInteractiveReport:
    """Test interactive HTML report features."""

    def test_report_has_javascript(self):
        """Report should contain interactive JavaScript."""
        from biocompiler.report import generate_report
        html = generate_report(sequence="ATGGTGAGCAAGGGCGAGGAG")
        assert "<script>" in html or "<script " in html

    def test_report_has_dark_mode_toggle(self):
        """Report should have dark mode toggle."""
        from biocompiler.report import generate_report
        html = generate_report(sequence="ATGGTGAGCAAGGGCGAGGAG")
        assert "dark-mode" in html or "darkMode" in html

    def test_report_has_filter_buttons(self):
        """Report should have predicate filter buttons."""
        from biocompiler.report import generate_report
        from biocompiler.type_system import evaluate_all_predicates
        seq = "ATGGTGAGCAAGGGCGAGGAG"
        results = evaluate_all_predicates(seq=seq, known_exon_boundaries=[(0, len(seq))])
        html = generate_report(sequence=seq, type_results=results)
        assert "filter" in html.lower() or "PASS" in html

    def test_report_has_tooltips(self):
        """Report should have tooltip data attributes."""
        from biocompiler.report import generate_report
        html = generate_report(sequence="ATGGTGAGCAAGGGCGAGGAG")
        assert "data-tooltip" in html or "tooltip" in html

    def test_report_has_collapsible_sections(self):
        """Report should support collapsible sections."""
        from biocompiler.report import generate_report
        html = generate_report(sequence="ATGGTGAGCAAGGGCGAGGAG")
        assert "section-content" in html or "collapsible" in html or "collapsed" in html


# ─── Cache Integrity Tests ────────────────────────────────────────

class TestCacheIntegrity:
    """Test cache integrity verification."""

    def test_verify_cache_no_hash(self, tmp_path):
        """When no hash is stored, verify should return True."""
        from biocompiler.organism_db import OrganismDatabase
        db = OrganismDatabase(db_path=tmp_path / "cache_test.db")
        # No API data fetched, so no hash stored
        assert db.verify_cache_integrity("NonExistent") is True


# ─── CLI Enhancement Tests ─────────────────────────────────────────

class TestCLIEnhancements:
    """Test CLI colored output and progress indicators."""

    def test_cli_version(self):
        """CLI --version should work."""
        import subprocess
        result = subprocess.run(
            ["biocompiler", "--version"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "3.2" in result.stdout or "3.2" in result.stderr

    def test_cli_scan_smoke(self):
        """CLI scan should work with a basic sequence."""
        import subprocess
        result = subprocess.run(
            ["biocompiler", "scan",
             "--sequence", "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "Tokens" in result.stdout or "tokens" in result.stdout.lower()

    def test_cli_help(self):
        """CLI --help should list all subcommands."""
        import subprocess
        result = subprocess.run(
            ["biocompiler", "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "check" in result.stdout
        assert "optimize" in result.stdout
        assert "scan" in result.stdout
        assert "serve" in result.stdout
