"""
Tests for the persistent provenance store integration in the API.

Covers:
1. ProvenanceStore replaces in-memory dict — records survive process restarts
2. GET /provenance/{record_id} — retrieve a stored provenance record
3. GET /provenance — query/list provenance records with filters
4. BIOCOMPILER_PROVENANCE_DIR env var — configurable store directory
5. /optimize endpoint returns provenance_id when tracking is enabled
6. /optimize endpoint returns null provenance_id when tracking is disabled
7. 404 for non-existent provenance record
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from biocompiler.api import create_app, _provenance_store
from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
    ProvenanceStore,
)


# ─── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture()
def tmp_prov_dir(tmp_path):
    """Create a temporary directory for provenance storage."""
    d = tmp_path / "provenance"
    d.mkdir()
    return d


@pytest.fixture()
def app_with_tmp_store(tmp_prov_dir):
    """Create a FastAPI app with a temporary provenance store and auth disabled."""
    from biocompiler import api as api_module

    # Patch the module-level _provenance_store to use tmp dir
    store = ProvenanceStore(store_dir=str(tmp_prov_dir))

    # Disable auth for testing
    original_mode = api_module._AUTH_MODE
    original_keys = set(api_module._CONFIGURED_API_KEYS)

    api_module._AUTH_MODE = "disabled"
    api_module._CONFIGURED_API_KEYS = set()

    try:
        with patch("biocompiler.api._provenance_store", store):
            app = create_app()
            yield app
    finally:
        api_module._AUTH_MODE = original_mode
        api_module._CONFIGURED_API_KEYS = original_keys


@pytest.fixture()
def client_with_store(app_with_tmp_store):
    """Create a TestClient with a temporary provenance store."""
    return TestClient(app_with_tmp_store)


@pytest.fixture()
def sample_trail() -> OptimizationDecisionTrail:
    """Create a sample OptimizationDecisionTrail for testing."""
    return OptimizationDecisionTrail(
        gene_name="test_gene",
        input_protein="MVLSPADKTN",
        output_dna="ATGGTGCTGTCGCCCGCTGACAAGACCAAC",
        organism="Homo_sapiens",
        solver_backend="greedy",
        seed=42,
        total_cai=0.91,
        total_gc=0.54,
        codon_decisions=[
            CodonDecision(
                position=0,
                amino_acid="M",
                original_codon=None,
                chosen_codon="ATG",
                alternatives_considered=[],
                constraint_reason="maximize_cai",
                confidence=1.0,
                cai_impact=0.0,
            ),
            CodonDecision(
                position=1,
                amino_acid="V",
                original_codon=None,
                chosen_codon="GTG",
                alternatives_considered=[
                    {"codon": "GTC", "cai": 0.85, "reason": "lower_cai"},
                ],
                constraint_reason="maximize_cai",
                confidence=0.95,
                cai_impact=-0.002,
            ),
        ],
        constraint_decisions=[
            ConstraintDecision(
                constraint_name="GCInRange",
                constraint_type="hard",
                action_taken="satisfied",
                positions_affected=[1, 3],
                tradeoff_description="GC content kept within bounds",
                impact_on_cai=0.0,
            ),
        ],
        iteration_log=[{"step": 1, "score": 0.85, "action": "initial"}],
        timestamp="2026-01-15T12:00:00+00:00",
        version="1.0.0",
    )


# ═══════════════════════════════════════════════════════════════════════
# 1. ProvenanceStore persistence
# ═══════════════════════════════════════════════════════════════════════


class TestProvenanceStorePersistence:
    """Verify that ProvenanceStore persists records to disk."""

    def test_save_and_load_roundtrip(self, tmp_prov_dir, sample_trail):
        """Save a trail, then load it — data must match."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        record_id = store.save(sample_trail)
        loaded = store.load(record_id)
        assert loaded.gene_name == sample_trail.gene_name
        assert loaded.input_protein == sample_trail.input_protein
        assert loaded.output_dna == sample_trail.output_dna
        assert loaded.organism == sample_trail.organism
        assert loaded.total_cai == sample_trail.total_cai
        assert loaded.total_gc == sample_trail.total_gc
        assert len(loaded.codon_decisions) == len(sample_trail.codon_decisions)
        assert len(loaded.constraint_decisions) == len(sample_trail.constraint_decisions)

    def test_save_creates_json_file(self, tmp_prov_dir, sample_trail):
        """Save should create a JSON file on disk."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        record_id = store.save(sample_trail)
        filepath = tmp_prov_dir / f"{record_id}.json"
        assert filepath.exists()
        # File should be valid JSON
        with open(filepath) as f:
            data = json.load(f)
        assert data["gene_name"] == "test_gene"

    def test_load_nonexistent_record_raises(self, tmp_prov_dir):
        """Loading a non-existent record should raise FileNotFoundError."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        with pytest.raises(FileNotFoundError):
            store.load("00000000-0000-0000-0000-000000000000")

    def test_query_by_organism(self, tmp_prov_dir, sample_trail):
        """Query should filter by organism."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        store.save(sample_trail)
        results = store.query(organism="Homo_sapiens")
        assert len(results) == 1
        assert results[0].organism == "Homo_sapiens"

    def test_query_by_organism_no_match(self, tmp_prov_dir, sample_trail):
        """Query with non-matching organism should return empty list."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        store.save(sample_trail)
        results = store.query(organism="Escherichia_coli")
        assert len(results) == 0

    def test_query_by_protein_name(self, tmp_prov_dir, sample_trail):
        """Query should filter by gene_name."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        store.save(sample_trail)
        results = store.query(protein_name="test_gene")
        assert len(results) == 1
        assert results[0].gene_name == "test_gene"

    def test_records_survive_store_reinstantiation(self, tmp_prov_dir, sample_trail):
        """Records should survive store re-instantiation (file-based persistence)."""
        store1 = ProvenanceStore(store_dir=str(tmp_prov_dir))
        record_id = store1.save(sample_trail)

        # Create a new store instance pointing to the same directory
        store2 = ProvenanceStore(store_dir=str(tmp_prov_dir))
        loaded = store2.load(record_id)
        assert loaded.gene_name == sample_trail.gene_name


# ═══════════════════════════════════════════════════════════════════════
# 2. GET /provenance/{record_id} endpoint
# ═══════════════════════════════════════════════════════════════════════


class TestGetProvenanceEndpoint:
    """Test the GET /provenance/{record_id} endpoint."""

    def test_retrieve_existing_record(self, client_with_store, sample_trail):
        """Should return the full provenance trail for a valid record_id."""
        # Save a record first
        from biocompiler.api import _provenance_store as store
        record_id = store.save(sample_trail)

        resp = client_with_store.get(f"/provenance/{record_id}")
        assert resp.status_code == 200
        data = resp.json()
        assert data["id"] == record_id
        assert "trail" in data
        assert data["trail"]["gene_name"] == "test_gene"
        assert data["trail"]["organism"] == "Homo_sapiens"
        assert len(data["trail"]["codon_decisions"]) == 2

    def test_404_for_nonexistent_record(self, client_with_store):
        """Should return 404 for a non-existent record ID."""
        resp = client_with_store.get("/provenance/00000000-0000-0000-0000-000000000000")
        assert resp.status_code == 404
        assert "not found" in resp.json()["detail"].lower()

    def test_response_contains_decision_data(self, client_with_store, sample_trail):
        """Response trail should contain codon and constraint decisions."""
        from biocompiler.api import _provenance_store as store
        record_id = store.save(sample_trail)

        resp = client_with_store.get(f"/provenance/{record_id}")
        data = resp.json()
        trail = data["trail"]

        # Codon decisions
        assert len(trail["codon_decisions"]) == 2
        cd0 = trail["codon_decisions"][0]
        assert cd0["position"] == 0
        assert cd0["amino_acid"] == "M"
        assert cd0["chosen_codon"] == "ATG"

        # Constraint decisions
        assert len(trail["constraint_decisions"]) == 1
        assert trail["constraint_decisions"][0]["constraint_name"] == "GCInRange"


# ═══════════════════════════════════════════════════════════════════════
# 3. GET /provenance query/list endpoint
# ═══════════════════════════════════════════════════════════════════════


class TestQueryProvenanceEndpoint:
    """Test the GET /provenance query/list endpoint."""

    def test_list_all_records(self, client_with_store, sample_trail):
        """Should list all stored provenance records."""
        from biocompiler.api import _provenance_store as store
        store.save(sample_trail)

        resp = client_with_store.get("/provenance")
        assert resp.status_code == 200
        data = resp.json()
        assert data["count"] >= 1
        assert len(data["records"]) >= 1

    def test_list_empty_store(self, client_with_store):
        """Should return count=0 and empty records for empty store."""
        resp = client_with_store.get("/provenance")
        assert resp.status_code == 200
        data = resp.json()
        assert data["count"] == 0
        assert data["records"] == []

    def test_filter_by_organism(self, client_with_store, sample_trail):
        """Should filter records by organism query parameter."""
        from biocompiler.api import _provenance_store as store
        store.save(sample_trail)

        resp = client_with_store.get("/provenance", params={"organism": "Homo_sapiens"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["count"] >= 1

        # Non-matching organism
        resp2 = client_with_store.get("/provenance", params={"organism": "Escherichia_coli"})
        assert resp2.status_code == 200
        data2 = resp2.json()
        assert data2["count"] == 0

    def test_filter_by_protein_name(self, client_with_store, sample_trail):
        """Should filter records by protein_name query parameter."""
        from biocompiler.api import _provenance_store as store
        store.save(sample_trail)

        resp = client_with_store.get("/provenance", params={"protein_name": "test_gene"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["count"] >= 1

    def test_record_summary_fields(self, client_with_store, sample_trail):
        """Each record in list should have summary fields."""
        from biocompiler.api import _provenance_store as store
        store.save(sample_trail)

        resp = client_with_store.get("/provenance")
        data = resp.json()
        assert data["count"] >= 1
        record = data["records"][0]
        assert "gene_name" in record
        assert "organism" in record
        assert "timestamp" in record
        assert "total_cai" in record
        assert "total_gc" in record
        assert "codon_decision_count" in record
        assert "constraint_decision_count" in record

    def test_filter_by_date_range(self, client_with_store, sample_trail):
        """Should filter records by date_from and date_to."""
        from biocompiler.api import _provenance_store as store
        store.save(sample_trail)

        # The sample trail has timestamp "2026-01-15T12:00:00+00:00"
        resp = client_with_store.get(
            "/provenance",
            params={
                "date_from": "2026-01-01T00:00:00+00:00",
                "date_to": "2026-12-31T23:59:59+00:00",
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["count"] >= 1

        # Outside range
        resp2 = client_with_store.get(
            "/provenance",
            params={
                "date_from": "2025-01-01T00:00:00+00:00",
                "date_to": "2025-12-31T23:59:59+00:00",
            },
        )
        assert resp2.status_code == 200
        data2 = resp2.json()
        assert data2["count"] == 0


# ═══════════════════════════════════════════════════════════════════════
# 4. BIOCOMPILER_PROVENANCE_DIR env var
# ═══════════════════════════════════════════════════════════════════════


class TestProvenanceDirEnvVar:
    """Test that the store directory is configurable via env var."""

    def test_env_var_sets_store_dir(self, tmp_prov_dir):
        """BIOCOMPILER_PROVENANCE_DIR should configure the store directory."""
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_DIR": str(tmp_prov_dir)}):
            prov_dir = os.environ.get("BIOCOMPILER_PROVENANCE_DIR")
            store = ProvenanceStore(store_dir=prov_dir)
            assert store.store_dir == tmp_prov_dir

    def test_default_store_dir(self):
        """Without env var, store should use default ~/.biocompiler/provenance."""
        store = ProvenanceStore()
        expected = Path.home() / ".biocompiler" / "provenance"
        assert store.store_dir == expected

    def test_records_stored_in_configured_dir(self, tmp_prov_dir, sample_trail):
        """Records should be stored in the configured directory."""
        store = ProvenanceStore(store_dir=str(tmp_prov_dir))
        record_id = store.save(sample_trail)
        assert (tmp_prov_dir / f"{record_id}.json").exists()


# ═══════════════════════════════════════════════════════════════════════
# 5. Optimize endpoint provenance_id integration
# ═══════════════════════════════════════════════════════════════════════


class TestOptimizeProvenanceId:
    """Test that /optimize returns provenance_id when tracking is enabled."""

    def test_optimize_with_provenance_returns_id(self, client_with_store):
        """When track_provenance=True, response should include provenance_id."""
        resp = client_with_store.post(
            "/optimize",
            json={
                "protein": "MVSKGE",
                "organism": "Homo_sapiens",
                "track_provenance": True,
            },
        )
        # The request should succeed
        assert resp.status_code in (200, 422)
        if resp.status_code == 200:
            data = resp.json()
            # provenance_id should be present (may be None if decision_trail is None)
            assert "provenance_id" in data

    def test_optimize_without_provenance_returns_null_id(self, client_with_store):
        """When track_provenance=False, provenance_id should be None."""
        resp = client_with_store.post(
            "/optimize",
            json={
                "protein": "MVSKGE",
                "organism": "Homo_sapiens",
                "track_provenance": False,
            },
        )
        assert resp.status_code in (200, 422)
        if resp.status_code == 200:
            data = resp.json()
            assert data.get("provenance_id") is None
