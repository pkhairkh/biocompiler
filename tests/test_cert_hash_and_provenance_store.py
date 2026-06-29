"""
Tests for certificate hash v2 soundness fix and ProvenanceStore path traversal fix.

ISSUE 1: Certificate SHA-256 hash covers only sequence, not predicate results.
  - v1 hash: SHA-256(sequence) — soundness bug
  - v2 hash: SHA-256(sequence + sorted predicate results + key opt params)

ISSUE 2: ProvenanceStore.load() has path traversal risk.
  - Must validate UUID format before using it as a file path component
  - Must reject IDs like "../../../etc/passwd"
"""

import hashlib
import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pytest

from biocompiler.provenance.certificate import (
    _compute_certificate_hash,
    _CURRENT_HASH_VERSION,
    _HASH_ALGORITHM,
    _V2_HASH_PARAM_KEYS,
    generate_certificate,
    verify_certificate,
)
from biocompiler.provenance import ProvenanceStore, ProvenanceTracker, DecisionRecord
from biocompiler.shared.types import Certificate, TypeCheckResult, Verdict, SLOTMode


# ── Shared fixtures ─────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"  # 45 bp
SAMPLE_SEQ_2 = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGAAA"  # differs at end


def _make_type_results(
    verdicts: list[tuple[str, Verdict]] | None = None,
) -> list[TypeCheckResult]:
    """Build a list of TypeCheckResult objects for testing."""
    if verdicts is None:
        verdicts = [
            ("GCInRange", Verdict.PASS),
            ("NoStopCodons", Verdict.PASS),
        ]
    return [
        TypeCheckResult(predicate=p, verdict=v, violation=None if v == Verdict.PASS else "See details")
        for p, v in verdicts
    ]


# ══════════════════════════════════════════════════════════════════════════
# 1. Certificate hash v2 — soundness tests
# ══════════════════════════════════════════════════════════════════════════

class TestCertificateHashV2:
    """Tests for the v2 certificate hash covering sequence + predicates + params."""

    def test_v2_hash_differs_from_v1_for_same_sequence(self):
        """v1 and v2 hashes should differ even for the same sequence and types."""
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        v1 = _compute_certificate_hash(SAMPLE_SEQ, types_list, params, hash_version=1)
        v2 = _compute_certificate_hash(SAMPLE_SEQ, types_list, params, hash_version=2)
        assert v1 != v2, "v1 and v2 hashes must differ for same inputs"

    def test_v1_hash_is_sequence_only(self):
        """v1 hash should be SHA-256(sequence), ignoring types and params."""
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params = {"organism": "Homo_sapiens"}
        v1 = _compute_certificate_hash(SAMPLE_SEQ, types_list, params, hash_version=1)
        expected = hashlib.sha256(SAMPLE_SEQ.encode()).hexdigest()
        assert v1 == expected

    def test_v2_hash_depends_on_predicate_verdicts(self):
        """Two certificates with same sequence but different predicate outcomes
        should get different v2 hashes (the core soundness fix)."""
        types_pass = [{"predicate": "GCInRange", "verdict": "PASS"}]
        types_fail = [{"predicate": "GCInRange", "verdict": "FAIL"}]
        params = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        hash_pass = _compute_certificate_hash(SAMPLE_SEQ, types_pass, params, hash_version=2)
        hash_fail = _compute_certificate_hash(SAMPLE_SEQ, types_fail, params, hash_version=2)
        assert hash_pass != hash_fail, (
            "v2 hash must differ when predicate verdicts differ"
        )

    def test_v2_hash_depends_on_organism(self):
        """v2 hash must differ when organism parameter changes."""
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params_human = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        params_ecoli = {"organism": "Escherichia_coli", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        hash_human = _compute_certificate_hash(SAMPLE_SEQ, types_list, params_human, hash_version=2)
        hash_ecoli = _compute_certificate_hash(SAMPLE_SEQ, types_list, params_ecoli, hash_version=2)
        assert hash_human != hash_ecoli

    def test_v2_hash_depends_on_gc_bounds(self):
        """v2 hash must differ when GC bounds change."""
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params_wide = {"organism": "Homo_sapiens", "gc_lo": 0.30, "gc_hi": 0.70, "cai_threshold": 0.5}
        params_narrow = {"organism": "Homo_sapiens", "gc_lo": 0.40, "gc_hi": 0.60, "cai_threshold": 0.5}
        hash_wide = _compute_certificate_hash(SAMPLE_SEQ, types_list, params_wide, hash_version=2)
        hash_narrow = _compute_certificate_hash(SAMPLE_SEQ, types_list, params_narrow, hash_version=2)
        assert hash_wide != hash_narrow

    def test_v2_hash_depends_on_cai_threshold(self):
        """v2 hash must differ when CAI threshold changes."""
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params_low = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        params_high = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.8}
        hash_low = _compute_certificate_hash(SAMPLE_SEQ, types_list, params_low, hash_version=2)
        hash_high = _compute_certificate_hash(SAMPLE_SEQ, types_list, params_high, hash_version=2)
        assert hash_low != hash_high

    def test_v2_hash_depends_on_sequence(self):
        """v2 hash must differ when sequence changes."""
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        hash1 = _compute_certificate_hash(SAMPLE_SEQ, types_list, params, hash_version=2)
        hash2 = _compute_certificate_hash(SAMPLE_SEQ_2, types_list, params, hash_version=2)
        assert hash1 != hash2

    def test_v2_hash_order_independent_predicates(self):
        """v2 hash should be the same regardless of predicate order (sorted)."""
        types_a = [
            {"predicate": "GCInRange", "verdict": "PASS"},
            {"predicate": "NoStopCodons", "verdict": "FAIL"},
        ]
        types_b = [
            {"predicate": "NoStopCodons", "verdict": "FAIL"},
            {"predicate": "GCInRange", "verdict": "PASS"},
        ]
        params = {"organism": "Homo_sapiens", "gc_lo": 0.3, "gc_hi": 0.7, "cai_threshold": 0.5}
        hash_a = _compute_certificate_hash(SAMPLE_SEQ, types_a, params, hash_version=2)
        hash_b = _compute_certificate_hash(SAMPLE_SEQ, types_b, params, hash_version=2)
        assert hash_a == hash_b, "v2 hash must be order-independent (sorted predicates)"

    def test_generate_certificate_uses_current_hash(self):
        """generate_certificate should produce a v3 hash by default.

        Since C14 (prior fix) bumped ``_CURRENT_HASH_VERSION`` to 3, freshly
        generated certificates carry ``hash_version=3``.
        """
        type_results = _make_type_results()
        cert = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"gene": "eGFP", "organism": "Homo_sapiens"},
        )
        assert cert.hash_version == _CURRENT_HASH_VERSION
        assert cert.hash_version == 3
        # The design_id should NOT be just SHA-256(sequence)
        seq_only_hash = hashlib.sha256(SAMPLE_SEQ.encode()).hexdigest()
        assert cert.design_id != seq_only_hash, (
            "v3 design_id must differ from sequence-only hash"
        )

    def test_generate_certificate_current_hash_matches_compute(self):
        """generate_certificate design_id should match _compute_certificate_hash v3.

        Since C14 (prior fix), freshly generated certificates use the v3 hash
        (expanded param set). v2 hashes are still computable for legacy
        verification (see ``TestVerifyCertificateHashVersions``).
        """
        type_results = _make_type_results()
        input_params = {"gene": "eGFP", "organism": "Homo_sapiens"}
        cert = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params=input_params,
        )
        # Reconstruct the same params that generate_certificate uses internally
        complete_params = dict(input_params)
        complete_params.setdefault("organism", "Homo_sapiens")
        complete_params.setdefault("cell_type", "HEK293T")
        complete_params.setdefault("gc_lo", 0.30)
        complete_params.setdefault("gc_hi", 0.70)
        complete_params.setdefault("cai_threshold", 0.5)
        complete_params.setdefault("enzymes", ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"])
        complete_params.setdefault("cryptic_splice_threshold", 3.0)  # maxentscan.CRYPTIC_SPLICE_THRESHOLD (default 3.0)
        complete_params.setdefault("exon_boundaries", [(0, len(SAMPLE_SEQ))])

        expected_hash = _compute_certificate_hash(
            sequence=SAMPLE_SEQ,
            types_list=cert.types,
            params=complete_params,
            hash_version=3,
        )
        assert cert.design_id == expected_hash

    def test_different_predicates_different_design_id(self):
        """End-to-end: same sequence, different predicate outcomes → different design_id."""
        type_results_pass = _make_type_results([("GCInRange", Verdict.PASS)])
        type_results_fail = _make_type_results([("GCInRange", Verdict.FAIL)])
        cert_pass = generate_certificate(SAMPLE_SEQ, type_results_pass, {"organism": "Homo_sapiens"})
        cert_fail = generate_certificate(SAMPLE_SEQ, type_results_fail, {"organism": "Homo_sapiens"})
        assert cert_pass.design_id != cert_fail.design_id


# ══════════════════════════════════════════════════════════════════════════
# 2. Certificate hash_version field
# ══════════════════════════════════════════════════════════════════════════

class TestCertificateHashVersion:
    """Tests for the hash_version field on the Certificate dataclass."""

    def test_new_certificate_has_current_hash_version(self):
        """Newly generated certificates should have hash_version=_CURRENT_HASH_VERSION (3)."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"gene": "eGFP"})
        assert cert.hash_version == _CURRENT_HASH_VERSION
        assert cert.hash_version == 3

    def test_certificate_to_dict_includes_hash_version_for_current(self):
        """to_dict should include hash_version for current (v3) certificates."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"gene": "eGFP"})
        d = cert.to_dict()
        assert "hash_version" in d
        assert d["hash_version"] == _CURRENT_HASH_VERSION
        assert d["hash_version"] == 3

    def test_certificate_to_dict_omits_hash_version_for_v1(self):
        """to_dict should NOT include hash_version for v1 (legacy) certificates."""
        cert = Certificate(
            version="1.0.0",
            design_id="abc123",
            sequence="ATGC",
            types=[{"predicate": "gc", "verdict": "PASS"}],
            provenance={"tool": "BC"},
            hash_version=1,
        )
        d = cert.to_dict()
        assert "hash_version" not in d, "v1 certs should omit hash_version for backward compat"

    def test_from_dict_defaults_to_v1_when_missing(self):
        """from_dict should default to hash_version=1 when field is absent."""
        d = {
            "version": "1.0.0",
            "design_id": "abc123",
            "sequence": "ATGC",
            "types": [{"predicate": "gc", "verdict": "PASS"}],
            "provenance": {"tool": "BC"},
        }
        cert = Certificate.from_dict(d)
        assert cert.hash_version == 1

    def test_from_dict_preserves_hash_version(self):
        """from_dict should preserve hash_version when present."""
        d = {
            "version": "1.0.0",
            "design_id": "abc123",
            "sequence": "ATGC",
            "types": [{"predicate": "gc", "verdict": "PASS"}],
            "provenance": {"tool": "BC"},
            # C14 (prior fix): use v3 (current) to match _CURRENT_HASH_VERSION.
            "hash_version": 3,
        }
        cert = Certificate.from_dict(d)
        assert cert.hash_version == 3

    def test_roundtrip_preserves_hash_version(self):
        """to_dict/from_dict roundtrip should preserve hash_version."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"gene": "eGFP"})
        d = cert.to_dict()
        restored = Certificate.from_dict(d)
        # C14 (prior fix): freshly generated certs now use v3.
        assert restored.hash_version == cert.hash_version == 3


# ══════════════════════════════════════════════════════════════════════════
# 3. verify_certificate backward compatibility with v1 and v2 hashes
# ══════════════════════════════════════════════════════════════════════════

class TestVerifyCertificateHashVersions:
    """Tests that verify_certificate handles both v1 and v2 hashes."""

    def test_verify_v2_certificate_succeeds(self):
        """verify_certificate should accept a valid v2 certificate."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"gene": "eGFP"})
        cert_dict = cert.to_dict()
        # We cannot fully verify without a working registry, but the hash check
        # should at least not reject on hash mismatch
        status, failures = verify_certificate(cert_dict)
        # The hash check should pass (design_id matches v2 computation)
        hash_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(hash_failures) == 0, f"Hash check failed: {hash_failures}"

    def test_verify_v1_legacy_certificate_hash_check(self):
        """verify_certificate should accept a v1 certificate with sequence-only hash."""
        # Construct a v1-style certificate (no hash_version field)
        seq_hash = hashlib.sha256(SAMPLE_SEQ.encode()).hexdigest()
        cert_dict = {
            "version": "1.0.0",
            "design_id": seq_hash,
            "sequence": SAMPLE_SEQ,
            "types": [],
            "provenance": {
                "tool": "BioCompiler",
                "version": "1.0.0",
                "timestamp": "2025-01-01T00:00:00Z",
                "input_hash": seq_hash,
                "parameters": {"organism": "Homo_sapiens"},
            },
            # No hash_version field → v1
        }
        status, failures = verify_certificate(cert_dict)
        hash_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(hash_failures) == 0, f"v1 hash check should pass: {hash_failures}"

    def test_verify_detects_tampered_v2_hash(self):
        """verify_certificate should detect a v2 certificate with wrong design_id."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"gene": "eGFP"})
        cert_dict = cert.to_dict()
        # Tamper with the design_id
        cert_dict["design_id"] = "0" * 64
        status, failures = verify_certificate(cert_dict)
        hash_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(hash_failures) > 0, "Should detect tampered v2 hash"

    def test_verify_detects_tampered_v1_hash(self):
        """verify_certificate should detect a v1 certificate with wrong design_id."""
        cert_dict = {
            "version": "1.0.0",
            "design_id": "0" * 64,  # Wrong hash
            "sequence": SAMPLE_SEQ,
            "types": [],
            "provenance": {
                "tool": "BioCompiler",
                "version": "1.0.0",
                "timestamp": "2025-01-01T00:00:00Z",
                "input_hash": "abc",
                "parameters": {"organism": "Homo_sapiens"},
            },
        }
        status, failures = verify_certificate(cert_dict)
        hash_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(hash_failures) > 0, "Should detect tampered v1 hash"


# ══════════════════════════════════════════════════════════════════════════
# 4. ProvenanceStore UUID validation — path traversal prevention
# ══════════════════════════════════════════════════════════════════════════

class TestProvenanceStoreUUIDValidation:
    """Tests for ProvenanceStore UUID validation preventing path traversal."""

    def test_valid_uuid_accepted(self):
        """A valid lowercase UUID should pass validation."""
        ProvenanceStore._validate_uuid("a1b2c3d4-e5f6-7890-abcd-ef1234567890")

    def test_path_traversal_dotdot_rejected(self):
        """../../../etc/passwd should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("../../../etc/passwd")

    def test_path_traversal_simple_rejected(self):
        """A simple non-UUID string should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("not-a-uuid")

    def test_path_traversal_absolute_path_rejected(self):
        """An absolute path should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("/etc/passwd")

    def test_path_traversal_mixed_rejected(self):
        """A string with path separators mixed in should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("a1b2c3d4/../../etc/passwd")

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

    def test_uuid_with_null_byte_rejected(self):
        """UUID with embedded null byte should be rejected."""
        with pytest.raises(ValueError, match="Invalid provenance record ID"):
            ProvenanceStore._validate_uuid("a1b2c3d4\x00-e5f6-7890-abcd-ef1234567890")


class TestProvenanceStoreLoad:
    """Integration tests for ProvenanceStore.load() with UUID validation."""

    def _make_tracker(self, seed: int = 42) -> ProvenanceTracker:
        """Create a simple ProvenanceTracker for testing."""
        tracker = ProvenanceTracker(seed=seed)
        tracker.record_decision(DecisionRecord(
            timestamp=datetime.now(timezone.utc).isoformat(),
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=["ATG"],
            rationale="Start codon",
            constraint_context={"cai": 1.0},
        ))
        return tracker

    def test_save_and_load_roundtrip(self):
        """Save then load should return equivalent tracker."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            tracker = self._make_tracker()
            record_id = store.save(tracker)
            loaded = store.load(record_id)
            assert loaded.seed == tracker.seed
            assert len(loaded) == len(tracker)

    def test_load_rejects_path_traversal(self):
        """load() should reject path traversal IDs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("../../../etc/passwd")

    def test_load_rejects_simple_string(self):
        """load() should reject non-UUID strings."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.load("my-record")

    def test_load_nonexistent_uuid_raises_file_not_found(self):
        """load() with a valid UUID that does not exist should raise FileNotFoundError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            with pytest.raises(FileNotFoundError):
                store.load("00000000-0000-0000-0000-000000000000")

    def test_delete_rejects_path_traversal(self):
        """delete() should reject path traversal IDs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            with pytest.raises(ValueError, match="Invalid provenance record ID"):
                store.delete("../../etc/shadow")

    def test_save_and_delete(self):
        """Save then delete should work with valid UUID."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            tracker = self._make_tracker()
            record_id = store.save(tracker)
            # Should exist
            loaded = store.load(record_id)
            assert loaded.seed == 42
            # Delete it
            store.delete(record_id)
            # Should no longer exist
            with pytest.raises(FileNotFoundError):
                store.load(record_id)

    def test_default_base_dir(self):
        """ProvenanceStore with no base_dir should use default path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            import os
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                store = ProvenanceStore()
                assert store._base_dir == Path(".") / "provenance_store"
            finally:
                os.chdir(old_cwd)

    def test_load_returns_correct_tracker_data(self):
        """Loaded tracker should have all original data intact."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(base_dir=tmpdir)
            tracker = ProvenanceTracker(seed=99)
            tracker.record_decision(DecisionRecord(
                timestamp="2025-06-15T12:00:00+00:00",
                decision_type="codon_selected",
                position=6,
                chosen_value="GCT",
                alternatives_considered=["GCC", "GCA"],
                rationale="Good CAI",
                constraint_context={"gc": 0.55},
                cai_impact=0.02,
                codon_before="GCC",
                codon_after="GCT",
            ))
            record_id = store.save(tracker)
            loaded = store.load(record_id)
            assert loaded.seed == 99
            assert len(loaded) == 1
            decisions = loaded.get_decisions_for_position(6)
            assert len(decisions) == 1
            assert decisions[0].chosen_value == "GCT"
            assert decisions[0].cai_impact == 0.02
            assert decisions[0].codon_before == "GCC"
            assert decisions[0].codon_after == "GCT"


# ══════════════════════════════════════════════════════════════════════════
# 5. Cross-cutting: verify_certificate with ProvenanceStore integration
# ══════════════════════════════════════════════════════════════════════════

class TestCertificateProvenanceIntegration:
    """Tests that certificates can be stored and verified via ProvenanceStore."""

    def test_certificate_stored_and_retrieved(self):
        """A certificate dict should survive JSON storage roundtrip."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"gene": "eGFP"})
        cert_dict = cert.to_dict()
        # Verify it is JSON-serializable
        json_str = json.dumps(cert_dict)
        restored = json.loads(json_str)
        # C14 (prior fix): freshly generated certs now use v3.
        assert restored["hash_version"] == 3
        assert restored["design_id"] == cert.design_id
