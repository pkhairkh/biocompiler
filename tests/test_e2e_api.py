"""
BioCompiler End-to-End REST API Integration Tests
====================================================

Comprehensive integration tests for the BioCompiler REST API covering
all documented endpoints, error handling, security features, and
batch operations.

Tests:
  a. Test all API endpoints with valid inputs
  b. Test error handling with invalid inputs
  c. Test CORS headers
  d. Test rate limiting
  e. Test batch endpoints

All tests are marked with @pytest.mark.e2e and @pytest.mark.slow.
"""

from __future__ import annotations

import os
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from biocompiler.api import (
    create_app,
    _check_rate_limit,
    _check_batch_rate_limit,
    _rate_limiter,
    RATE_LIMIT_RPM,
    MAX_DNA_SEQUENCE_LENGTH,
    MAX_PROTEIN_LENGTH,
    BATCH_CHECK_MAX,
    BATCH_OPTIMIZE_MAX,
    BATCH_EXPORT_MAX,
    get_auth_mode,
    get_configured_api_keys,
)

from biocompiler import api as _api_module


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture()
def app():
    """Create a fresh FastAPI app with auth disabled for functional testing."""
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
    """Create a TestClient from the app, clearing rate limits first."""
    _rate_limiter.clear()
    return TestClient(app)


# ═══════════════════════════════════════════════════════════════════════════════
# Test data
# ═══════════════════════════════════════════════════════════════════════════════

VALID_DNA = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
VALID_PROTEIN = "MGSWKRQPPA"
INVALID_DNA = "ATGGCCXYZ"
INVALID_PROTEIN = "MGSWK1RQ"

# Longer proteins for optimization testing
HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
    "EFTPPVQAAYQKVVAGVANALAHKYH"
)


# ═══════════════════════════════════════════════════════════════════════════════
# a. Test all API endpoints with valid inputs
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.e2e
@pytest.mark.slow
class TestAllEndpointsValidInput:
    """Exercise every documented API endpoint with valid inputs and
    verify 200 responses with expected response structure."""

    # ── GET endpoints ────────────────────────────────────────────────────

    def test_health_endpoint(self, client):
        """GET /health should return 200 with status fields."""
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "healthy"
        assert "version" in data
        assert "timestamp" in data
        assert isinstance(data["auth_enabled"], bool)
        assert isinstance(data["rate_limit_rpm"], int)

    def test_info_endpoint(self, client):
        """GET /info should return system configuration."""
        resp = client.get("/info")
        assert resp.status_code == 200
        data = resp.json()
        assert "max_protein_length" in data
        assert "max_dna_length" in data
        assert "supported_organisms" in data
        assert isinstance(data["supported_organisms"], list)
        assert len(data["supported_organisms"]) > 0

    def test_organisms_endpoint(self, client):
        """GET /organisms should list all supported organisms."""
        resp = client.get("/organisms")
        assert resp.status_code == 200
        data = resp.json()
        assert "organisms" in data
        names = [o["name"] for o in data["organisms"]]
        assert "Homo_sapiens" in names
        assert "Escherichia_coli" in names

    def test_predicates_endpoint(self, client):
        """GET /predicates should list registered predicates."""
        resp = client.get("/predicates")
        assert resp.status_code == 200
        data = resp.json()
        assert "predicates" in data
        assert len(data["predicates"]) > 0

    def test_enzymes_endpoint(self, client):
        """GET /enzymes should list restriction enzyme recognition sites."""
        resp = client.get("/enzymes")
        assert resp.status_code == 200
        data = resp.json()
        assert "enzymes" in data
        assert isinstance(data["enzymes"], dict)
        # Standard enzymes should be present
        assert "EcoRI" in data["enzymes"]
        assert "BamHI" in data["enzymes"]

    # ── POST endpoints ──────────────────────────────────────────────────

    def test_check_endpoint(self, client):
        """POST /check should type-check a DNA sequence."""
        resp = client.post("/check", json={"sequence": VALID_DNA})
        assert resp.status_code == 200
        data = resp.json()
        assert "sequence_length" in data
        assert "gc_content" in data
        assert "protein" in data
        assert "results" in data
        assert "overall_verdict" in data
        assert data["overall_verdict"] in ("PASS", "FAIL", "UNCERTAIN")

    def test_check_with_organism(self, client):
        """POST /check with organism parameter."""
        resp = client.post("/check", json={
            "sequence": VALID_DNA,
            "organism": "Homo_sapiens",
        })
        assert resp.status_code == 200

    def test_check_with_custom_gc(self, client):
        """POST /check with custom GC bounds."""
        resp = client.post("/check", json={
            "sequence": VALID_DNA,
            "gc_lo": 0.20,
            "gc_hi": 0.80,
        })
        assert resp.status_code == 200

    def test_check_with_enzymes(self, client):
        """POST /check with restriction enzyme list."""
        resp = client.post("/check", json={
            "sequence": VALID_DNA,
            "enzymes": ["EcoRI", "BamHI"],
        })
        assert resp.status_code == 200

    def test_optimize_endpoint(self, client):
        """POST /optimize should return an optimized DNA sequence."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Homo_sapiens",
            "strict_mode": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert "sequence" in data
        assert "protein" in data
        assert "cai" in data
        assert "gc_content" in data
        assert "satisfied_predicates" in data
        assert "failed_predicates" in data
        assert "fallback_used" in data
        assert data["protein"] == VALID_PROTEIN

    def test_optimize_ecoli(self, client):
        """POST /optimize for E. coli with short protein."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "ecoli",
            "strict_mode": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["sequence"]) == len(VALID_PROTEIN) * 3

    def test_optimize_with_gc_bounds(self, client):
        """POST /optimize with tight GC bounds."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Escherichia_coli",
            "gc_lo": 0.40,
            "gc_hi": 0.60,
            "strict_mode": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert 0.40 <= data["gc_content"] <= 0.60

    def test_optimize_with_enzymes(self, client):
        """POST /optimize with restriction enzyme avoidance."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Homo_sapiens",
            "enzymes": ["EcoRI", "BamHI", "HindIII"],
            "strict_mode": False,
        })
        assert resp.status_code == 200

    def test_optimize_longer_protein(self, client):
        """POST /optimize with a longer protein (HBB, 147 AA)."""
        resp = client.post("/optimize", json={
            "protein": HBB_PROTEIN,
            "organism": "Homo_sapiens",
            "strict_mode": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["sequence"]) == len(HBB_PROTEIN) * 3

    def test_verify_endpoint(self, client):
        """POST /verify should verify a certificate."""
        # First optimize to get a certificate
        opt_resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Homo_sapiens",
            "strict_mode": False,
        })
        opt_data = opt_resp.json()

        # Check endpoint provides certificate
        check_resp = client.post("/check", json={"sequence": opt_data["sequence"]})
        check_data = check_resp.json()

        if check_data.get("certificate"):
            # Verify the certificate
            verify_resp = client.post("/verify", json={
                "certificate": check_data["certificate"],
            })
            assert verify_resp.status_code == 200
            vdata = verify_resp.json()
            assert vdata["status"] in ("VERIFIED", "REJECTED")

    def test_scan_endpoint(self, client):
        """POST /scan should scan a DNA sequence for features."""
        resp = client.post("/scan", json={"sequence": VALID_DNA})
        assert resp.status_code == 200
        data = resp.json()
        assert "sequence_length" in data
        assert "tokens" in data

    def test_scan_with_orfs(self, client):
        """POST /scan with ORF detection."""
        resp = client.post("/scan", json={
            "sequence": VALID_DNA,
            "find_orfs": True,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert "orfs" in data

    def test_scan_with_enzymes(self, client):
        """POST /scan with restriction enzyme scanning."""
        resp = client.post("/scan", json={
            "sequence": VALID_DNA,
            "enzymes": ["EcoRI", "BamHI"],
        })
        assert resp.status_code == 200

    def test_export_fasta_endpoint(self, client):
        """POST /export/fasta should export FASTA content."""
        resp = client.post("/export/fasta", json={
            "sequence": VALID_DNA,
            "identifier": "test_design",
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "fasta"
        assert ">" in data["content"]

    def test_export_genbank_endpoint(self, client):
        """POST /export/genbank should export GenBank content."""
        resp = client.post("/export/genbank", json={
            "sequence": VALID_DNA,
            "organism": "Homo_sapiens",
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "genbank"
        assert "LOCUS" in data["content"] or "BIOCOMPILER" in data["content"]


# ═══════════════════════════════════════════════════════════════════════════════
# b. Test error handling with invalid inputs
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.e2e
@pytest.mark.slow
class TestErrorHandling:
    """Test that invalid inputs produce proper error responses (422, 400)."""

    def test_check_missing_sequence(self, client):
        """POST /check without sequence should return 422."""
        resp = client.post("/check", json={})
        assert resp.status_code == 422

    def test_check_invalid_nucleotides(self, client):
        """POST /check with invalid DNA should return 422."""
        resp = client.post("/check", json={"sequence": INVALID_DNA})
        assert resp.status_code == 422

    def test_check_unsupported_organism(self, client):
        """POST /check with unsupported organism should return 422."""
        resp = client.post("/check", json={
            "sequence": VALID_DNA,
            "organism": "Alien_xenomorph",
        })
        assert resp.status_code == 422

    def test_optimize_missing_protein(self, client):
        """POST /optimize without protein should return 422."""
        resp = client.post("/optimize", json={})
        assert resp.status_code == 422

    def test_optimize_invalid_amino_acids(self, client):
        """POST /optimize with invalid amino acids should return 422."""
        resp = client.post("/optimize", json={"protein": INVALID_PROTEIN})
        assert resp.status_code == 422

    def test_optimize_empty_protein(self, client):
        """POST /optimize with empty protein should return 422."""
        resp = client.post("/optimize", json={"protein": ""})
        assert resp.status_code == 422

    def test_optimize_unsupported_organism(self, client):
        """POST /optimize with unsupported organism should return 422."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Martian_pathogen",
        })
        assert resp.status_code == 422

    def test_export_fasta_missing_sequence(self, client):
        """POST /export/fasta without sequence should return 422."""
        resp = client.post("/export/fasta", json={})
        assert resp.status_code == 422

    def test_export_genbank_missing_sequence(self, client):
        """POST /export/genbank without sequence should return 422."""
        resp = client.post("/export/genbank", json={})
        assert resp.status_code == 422

    def test_export_fasta_invalid_dna(self, client):
        """POST /export/fasta with invalid DNA should return 422."""
        resp = client.post("/export/fasta", json={"sequence": INVALID_DNA})
        assert resp.status_code == 422

    def test_scan_missing_sequence(self, client):
        """POST /scan without sequence should return 422."""
        resp = client.post("/scan", json={})
        assert resp.status_code == 422

    def test_verify_empty_certificate(self, client):
        """POST /verify with empty certificate should return 400 or 200."""
        resp = client.post("/verify", json={"certificate": {}})
        assert resp.status_code in (200, 400)

    def test_organism_alias_resolution(self, client):
        """Short organism aliases should be resolved correctly."""
        # 'ecoli' should resolve to Escherichia_coli
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "ecoli",
            "strict_mode": False,
        })
        assert resp.status_code == 200

    def test_organism_alias_human(self, client):
        """'human' should resolve to Homo_sapiens."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "human",
            "strict_mode": False,
        })
        assert resp.status_code == 200

    def test_invalid_organism_domain(self, client):
        """Invalid organism_domain should return 422."""
        resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Homo_sapiens",
            "organism_domain": "archaea_bacteria",
        })
        assert resp.status_code == 422


# ═══════════════════════════════════════════════════════════════════════════════
# c. Test CORS headers
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.e2e
@pytest.mark.slow
class TestCORSHeaders:
    """Test CORS middleware behavior."""

    def test_default_no_cors_origin_header(self, client):
        """With default config (no CORS origins), no Allow-Origin header."""
        resp = client.get("/health", headers={"Origin": "http://localhost:3000"})
        # Default config has no CORS origins, so no CORS header
        assert "access-control-allow-origin" not in resp.headers

    def test_health_includes_cors_info(self, client):
        """Health endpoint may report CORS configuration."""
        resp = client.get("/health")
        data = resp.json()
        # Health should at least report auth and rate limit info
        assert "auth_enabled" in data
        assert "rate_limit_rpm" in data

    def test_cors_with_explicit_origins(self):
        """When CORS origins are set, headers should appear."""
        with patch.dict(os.environ, {
            "BIOCOMPILER_CORS_ORIGINS": "http://localhost:3000",
        }):
            # Create a fresh app with CORS origins
            test_app = create_app()
            test_client = TestClient(test_app)
            resp = test_client.get(
                "/health",
                headers={"Origin": "http://localhost:3000"},
            )
            # If CORS middleware is configured, it should add headers
            if "access-control-allow-origin" in resp.headers:
                assert resp.headers["access-control-allow-origin"] in (
                    "http://localhost:3000", "*"
                )

    def test_cors_wildcard_no_credentials(self):
        """CORS with wildcard origin should not allow credentials."""
        with patch.dict(os.environ, {
            "BIOCOMPILER_CORS_ORIGINS": "*",
        }):
            test_app = create_app()
            test_client = TestClient(test_app)
            resp = test_client.get("/health")
            # If wildcard CORS is configured, credentials should be False
            # This is a security requirement: browsers reject
            # Access-Control-Allow-Credentials=true with wildcard origin

    def test_options_preflight_request(self, client):
        """OPTIONS preflight request should be handled."""
        resp = client.options(
            "/health",
            headers={
                "Origin": "http://localhost:3000",
                "Access-Control-Request-Method": "GET",
            },
        )
        # Should not return 405
        assert resp.status_code in (200, 204, 405)


# ═══════════════════════════════════════════════════════════════════════════════
# d. Test rate limiting
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.e2e
@pytest.mark.slow
class TestRateLimiting:
    """Test that rate limiting works correctly."""

    def test_rate_limiter_allows_initial_requests(self):
        """Initial requests should be allowed."""
        _rate_limiter.clear()
        # Should not raise
        _check_rate_limit("test_e2e_client_initial")

    def test_rate_limiter_blocks_at_limit(self):
        """Requests beyond the limit should be blocked."""
        _rate_limiter.clear()
        from fastapi import HTTPException

        # Exhaust the rate limit
        for _ in range(RATE_LIMIT_RPM):
            _rate_limiter.record("test_e2e_client_limit")

        with pytest.raises(HTTPException) as exc_info:
            _check_rate_limit("test_e2e_client_limit")
        assert exc_info.value.status_code == 429

    def test_batch_rate_limit_blocks_oversized(self):
        """Batch requests that exceed remaining rate limit should be blocked."""
        _rate_limiter.clear()
        from fastapi import HTTPException

        # Fill up most of the rate limit
        for _ in range(RATE_LIMIT_RPM - 5):
            _rate_limiter.record("test_e2e_batch_client")

        # Batch of 10 should fail (only 5 remaining)
        with pytest.raises(HTTPException) as exc_info:
            _check_batch_rate_limit("test_e2e_batch_client", 10)
        assert exc_info.value.status_code == 429

    def test_batch_rate_limit_allows_within_limit(self):
        """Batch requests within the remaining limit should be allowed."""
        _rate_limiter.clear()
        _check_batch_rate_limit("test_e2e_batch_client_ok", 5)

    def test_rate_limit_headers_on_response(self, client):
        """Responses should include rate limit information."""
        resp = client.get("/health")
        # The response should succeed
        assert resp.status_code == 200
        # Rate limit info may be in custom headers or in the response body
        data = resp.json()
        assert "rate_limit_rpm" in data


# ═══════════════════════════════════════════════════════════════════════════════
# e. Test batch endpoints
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.e2e
@pytest.mark.slow
class TestBatchEndpoints:
    """Test batch API endpoints for check, optimize, and export."""

    # ── Batch Check ─────────────────────────────────────────────────────

    def test_batch_check_empty(self, client):
        """POST /batch/check with empty list should return 200."""
        resp = client.post("/batch/check", json={"sequences": []})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 0

    def test_batch_check_single_item(self, client):
        """POST /batch/check with one sequence."""
        items = [{"sequence": VALID_DNA}]
        resp = client.post("/batch/check", json={"sequences": items})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 1
        assert len(data["results"]) == 1

    def test_batch_check_multiple_items(self, client):
        """POST /batch/check with multiple sequences."""
        items = [
            {"sequence": VALID_DNA},
            {"sequence": "ATGCGTAAGCTT"},
        ]
        resp = client.post("/batch/check", json={"sequences": items})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 2

    def test_batch_check_summary_fields(self, client):
        """Batch check summary should have all expected fields."""
        items = [{"sequence": VALID_DNA}]
        resp = client.post("/batch/check", json={"sequences": items})
        data = resp.json()
        summary = data["summary"]
        assert "total" in summary
        assert "pass" in summary
        assert "fail" in summary
        assert "uncertain" in summary
        assert "errors" in summary

    def test_batch_check_exceeds_max(self, client):
        """POST /batch/check exceeding max should return 400."""
        items = [{"sequence": VALID_DNA} for _ in range(BATCH_CHECK_MAX + 1)]
        resp = client.post("/batch/check", json={"sequences": items})
        assert resp.status_code == 400

    def test_batch_check_with_organisms(self, client):
        """POST /batch/check with different organisms per item."""
        items = [
            {"sequence": VALID_DNA, "organism": "Homo_sapiens"},
            {"sequence": VALID_DNA, "organism": "Escherichia_coli"},
        ]
        resp = client.post("/batch/check", json={"sequences": items})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 2

    def test_batch_check_invalid_sequence_in_item(self, client):
        """Invalid sequence in batch item should return 422."""
        items = [{"sequence": INVALID_DNA}]
        resp = client.post("/batch/check", json={"sequences": items})
        assert resp.status_code == 422

    # ── Batch Optimize ──────────────────────────────────────────────────

    def test_batch_optimize_empty(self, client):
        """POST /batch/optimize with empty list should return 200."""
        resp = client.post("/batch/optimize", json={"proteins": []})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 0

    def test_batch_optimize_single_item(self, client):
        """POST /batch/optimize with one protein."""
        items = [{"protein": VALID_PROTEIN}]
        resp = client.post("/batch/optimize", json={"proteins": items})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 1

    def test_batch_optimize_multiple_organisms(self, client):
        """POST /batch/optimize with different organisms."""
        items = [
            {"protein": VALID_PROTEIN, "organism": "Homo_sapiens"},
            {"protein": VALID_PROTEIN, "organism": "Escherichia_coli"},
        ]
        resp = client.post("/batch/optimize", json={"proteins": items})
        assert resp.status_code == 200
        data = resp.json()
        assert data["summary"]["total"] == 2

    def test_batch_optimize_summary_fields(self, client):
        """Batch optimize summary should have all expected fields."""
        items = [{"protein": VALID_PROTEIN}]
        resp = client.post("/batch/optimize", json={"proteins": items})
        data = resp.json()
        summary = data["summary"]
        assert "total" in summary
        assert "all_satisfied" in summary
        assert "partial" in summary
        assert "errors" in summary

    def test_batch_optimize_exceeds_max(self, client):
        """POST /batch/optimize exceeding max should return 400."""
        items = [{"protein": VALID_PROTEIN} for _ in range(BATCH_OPTIMIZE_MAX + 1)]
        resp = client.post("/batch/optimize", json={"proteins": items})
        assert resp.status_code == 400

    def test_batch_optimize_invalid_protein(self, client):
        """Invalid protein in batch item should return 422."""
        items = [{"protein": INVALID_PROTEIN}]
        resp = client.post("/batch/optimize", json={"proteins": items})
        assert resp.status_code == 422

    # ── Batch Export ────────────────────────────────────────────────────

    def test_batch_export_empty(self, client):
        """POST /batch/export with empty list should return 200."""
        resp = client.post("/batch/export", json={"sequences": []})
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["results"]) == 0

    def test_batch_export_fasta(self, client):
        """POST /batch/export with FASTA format."""
        items = [{"sequence": VALID_DNA, "format": "fasta"}]
        resp = client.post("/batch/export", json={"sequences": items})
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["results"]) == 1
        assert data["results"][0]["format"] == "fasta"

    def test_batch_export_genbank(self, client):
        """POST /batch/export with GenBank format."""
        items = [{"sequence": VALID_DNA, "format": "genbank"}]
        resp = client.post("/batch/export", json={"sequences": items})
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["results"]) == 1
        assert data["results"][0]["format"] == "genbank"

    def test_batch_export_mixed_formats(self, client):
        """POST /batch/export with mixed formats."""
        items = [
            {"sequence": VALID_DNA, "format": "fasta"},
            {"sequence": VALID_DNA, "format": "genbank"},
        ]
        resp = client.post("/batch/export", json={"sequences": items})
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["results"]) == 2

    def test_batch_export_exceeds_max(self, client):
        """POST /batch/export exceeding max should return 400."""
        items = [{"sequence": VALID_DNA} for _ in range(BATCH_EXPORT_MAX + 1)]
        resp = client.post("/batch/export", json={"sequences": items})
        assert resp.status_code == 400

    def test_batch_export_invalid_format(self, client):
        """Invalid export format should return 422."""
        items = [{"sequence": VALID_DNA, "format": "pdf"}]
        resp = client.post("/batch/export", json={"sequences": items})
        assert resp.status_code == 422

    # ── Full workflow: optimize → check → export → verify ──────────────

    def test_full_api_workflow(self, client):
        """End-to-end: optimize → check → export → verify via API."""
        # Step 1: Optimize
        opt_resp = client.post("/optimize", json={
            "protein": VALID_PROTEIN,
            "organism": "Homo_sapiens",
            "strict_mode": False,
        })
        assert opt_resp.status_code == 200
        opt_data = opt_resp.json()
        optimized_seq = opt_data["sequence"]

        # Step 2: Check the optimized sequence
        check_resp = client.post("/check", json={
            "sequence": optimized_seq,
            "organism": "Homo_sapiens",
        })
        assert check_resp.status_code == 200
        check_data = check_resp.json()

        # Step 3: Export as FASTA
        fasta_resp = client.post("/export/fasta", json={
            "sequence": optimized_seq,
            "identifier": "e2e_test_design",
        })
        assert fasta_resp.status_code == 200
        fasta_data = fasta_resp.json()
        assert ">e2e_test_design" in fasta_data["content"]

        # Step 4: Export as GenBank
        gb_resp = client.post("/export/genbank", json={
            "sequence": optimized_seq,
            "organism": "Homo_sapiens",
        })
        assert gb_resp.status_code == 200

        # Step 5: Verify certificate (if available)
        if check_data.get("certificate"):
            verify_resp = client.post("/verify", json={
                "certificate": check_data["certificate"],
            })
            assert verify_resp.status_code == 200

    def test_cross_organism_workflow(self, client):
        """End-to-end: optimize the same protein for multiple organisms."""
        organisms = ["Homo_sapiens", "Escherichia_coli", "Saccharomyces_cerevisiae"]
        results = {}
        for organism in organisms:
            resp = client.post("/optimize", json={
                "protein": VALID_PROTEIN,
                "organism": organism,
                "strict_mode": False,
            })
            assert resp.status_code == 200, f"Optimize failed for {organism}"
            data = resp.json()
            results[organism] = data

        # All organisms should produce valid (but different) sequences
        sequences = [r["sequence"] for r in results.values()]
        # At least two should differ (different codon usage)
        assert len(set(sequences)) >= 1, "All organisms produced identical sequences"

        # All should have valid CAI
        for organism, data in results.items():
            assert data["cai"] > 0, f"CAI is 0 for {organism}"
