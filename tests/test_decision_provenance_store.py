"""
Tests for decision_provenance.ProvenanceStore path traversal fix.

Verifies that ProvenanceStore.load() correctly validates UUID record IDs
and rejects path traversal attempts, ensuring no arbitrary file read is
possible through a malicious record_id parameter.

Covers:
1. Valid UUID accepted in load() and save/load roundtrip
2. Path traversal via "../etc/passwd" is rejected
3. Path traversal via "..\\windows\\system32" is rejected
4. UUIDs with path separators (/) are rejected
5. UUIDs with backslash path separators are rejected
6. Empty string, non-UUID strings, uppercase UUIDs are rejected
7. export_audit_trail also validates record_id (delegates to load)
8. Path containment check detects escape from store directory
"""

from __future__ import annotations

import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pytest

from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
    ProvenanceStore,
)


# ── Shared helpers ─────────────────────────────────────────────────────

def _make_trail() -> OptimizationDecisionTrail:
    """Create a minimal OptimizationDecisionTrail for testing."""
    return OptimizationDecisionTrail(
        gene_name="TestGene",
        input_protein="MV",
        output_dna="ATGGTG",
        organism="Homo_sapiens",
        solver_backend="greedy",
        seed=42,
        total_cai=0.91,
        total_gc=0.54,
        codon_decisions=[
            CodonDecision(
                position=0, amino_acid="M", original_codon=None,
                chosen_codon="ATG", alternatives_considered=[],
                constraint_reason="maximize_cai", confidence=1.0,
            ),
        ],
        constraint_decisions=[
            ConstraintDecision(
                constraint_name="GCInRange",
                constraint_type="hard",
                action_taken="satisfied",
                positions_affected=[0],
                tradeoff_description="GC in range",
                impact_on_cai=-0.005,
            ),
        ],
        iteration_log=[],
        timestamp="2025-01-15T12:00:00+00:00",
        version="10.0.0",
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 1. UUID validation — static method tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceStoreUUIDValidation:
    """Tests for ProvenanceStore._validate_uuid() preventing path traversal."""

    def test_valid_uuid_accepted(self):
        """A valid lowercase UUID should pass validation without raising."""
        ProvenanceStore._validate_uuid("a1b2c3d4-e5f6-7890-abcd-ef1234567890")
        # If we get here, the UUID was accepted (no exception raised)
        assert True

    def test_path_traversal_dotdot_slash_rejected(self):
        """../../../etc/passwd should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("../../../etc/passwd")

    def test_path_traversal_dotdot_backslash_rejected(self):
        """..\\windows\\system32 should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("..\\windows\\system32")

    def test_forward_slash_rejected(self):
        """IDs containing forward slash should be rejected."""
        with pytest.raises(ValueError, match="path separators"):
            ProvenanceStore._validate_uuid("a1b2c3d4/../../etc/passwd")

    def test_backslash_rejected(self):
        """IDs containing backslash should be rejected."""
        with pytest.raises(ValueError, match="path separators"):
            ProvenanceStore._validate_uuid("a1b2c3d4\\..\\etc")

    def test_double_dot_rejected(self):
        """IDs containing '..' should be rejected even without separators."""
        with pytest.raises(ValueError, match="\\.\\."):
            ProvenanceStore._validate_uuid("a1b2..d4-e5f6-7890-abcd-ef1234567890")

    def test_simple_non_uuid_rejected(self):
        """A simple non-UUID string should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("not-a-uuid")

    def test_uppercase_uuid_rejected(self):
        """Uppercase UUID should be rejected (only lowercase hex accepted)."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("A1B2C3D4-E5F6-7890-ABCD-EF1234567890")

    def test_uuid_with_braces_rejected(self):
        """UUID with curly braces should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("{a1b2c3d4-e5f6-7890-abcd-ef1234567890}")

    def test_uuid_no_hyphens_rejected(self):
        """UUID without hyphens should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("a1b2c3d4e5f67890abcdef1234567890")

    def test_empty_string_rejected(self):
        """Empty string should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("")

    def test_absolute_path_rejected(self):
        """An absolute path should be rejected."""
        with pytest.raises(ValueError, match="path separators"):
            ProvenanceStore._validate_uuid("/etc/passwd")

    def test_null_byte_in_uuid_rejected(self):
        """UUID with embedded null byte should be rejected (fails regex)."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("a1b2c3d4\x00-e5f6-7890-abcd-ef1234567890")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ProvenanceStore.load() integration tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceStoreLoadPathTraversal:
    """Integration tests for ProvenanceStore.load() with UUID validation."""

    def test_save_and_load_roundtrip(self):
        """Save then load should return equivalent trail."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)
            loaded = store.load(record_id)
            assert loaded.gene_name == trail.gene_name
            assert loaded.total_cai == trail.total_cai
            assert len(loaded.codon_decisions) == len(trail.codon_decisions)

    def test_load_rejects_dotdot_path_traversal(self):
        """load() should reject ../../../etc/passwd."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("../../../etc/passwd")

    def test_load_rejects_backslash_path_traversal(self):
        """load() should reject ..\\windows\\system32."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("..\\windows\\system32")

    def test_load_rejects_forward_slash_in_id(self):
        """load() should reject IDs with forward slashes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="path separators"):
                store.load("a1b2c3d4/e5f6-7890-abcd-ef1234567890")

    def test_load_rejects_backslash_in_id(self):
        """load() should reject IDs with backslashes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="path separators"):
                store.load("a1b2c3d4\\e5f6-7890-abcd-ef1234567890")

    def test_load_rejects_simple_string(self):
        """load() should reject non-UUID strings."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("my-record")

    def test_load_rejects_uppercase_uuid(self):
        """load() should reject uppercase UUIDs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("A1B2C3D4-E5F6-7890-ABCD-EF1234567890")

    def test_load_rejects_uuid_with_braces(self):
        """load() should reject UUIDs with curly braces."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("{a1b2c3d4-e5f6-7890-abcd-ef1234567890}")

    def test_load_nonexistent_uuid_raises_file_not_found(self):
        """load() with a valid UUID that doesn't exist should raise FileNotFoundError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(FileNotFoundError):
                store.load("00000000-0000-0000-0000-000000000000")

    def test_load_returns_correct_trail_data(self):
        """Loaded trail should have all original data intact."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)
            loaded = store.load(record_id)
            assert loaded.gene_name == "TestGene"
            assert loaded.organism == "Homo_sapiens"
            assert loaded.total_cai == 0.91
            assert len(loaded.codon_decisions) == 1
            assert loaded.codon_decisions[0].chosen_codon == "ATG"
            assert len(loaded.constraint_decisions) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# 3. export_audit_trail also validates record_id
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportAuditTrailPathTraversal:
    """Verify that export_audit_trail() also rejects path traversal IDs
    (since it delegates to load())."""

    def test_export_audit_trail_rejects_path_traversal(self):
        """export_audit_trail should reject path traversal IDs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.export_audit_trail("../../../etc/passwd")

    def test_export_audit_trail_rejects_backslash_traversal(self):
        """export_audit_trail should reject backslash traversal IDs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.export_audit_trail("..\\windows\\system32")

    def test_export_audit_trail_valid_uuid_works(self):
        """export_audit_trail should work with a valid saved record."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)
            result = store.export_audit_trail(record_id)
            assert isinstance(result, str)
            # JSON format by default
            parsed = json.loads(result)
            assert parsed["gene_name"] == "TestGene"


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Path containment check (defense-in-depth)
# ═══════════════════════════════════════════════════════════════════════════════

class TestPathContainmentCheck:
    """Tests for the _safe_record_path() containment check."""

    def test_safe_record_path_returns_path_for_valid_uuid(self):
        """_safe_record_path should return a valid path for a legitimate UUID."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            path = store._safe_record_path("a1b2c3d4-e5f6-7890-abcd-ef1234567890")
            assert str(path).startswith(str(store._store_dir))
            assert path.name == "a1b2c3d4-e5f6-7890-abcd-ef1234567890.json"

    def test_safe_record_path_rejects_traversal(self):
        """_safe_record_path should reject path traversal IDs before path check."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store._safe_record_path("../../../etc/passwd")
