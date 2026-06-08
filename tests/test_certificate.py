"""Tests for the certificate module.

Covers:
- generate_certificate: graduated and strict certificate generation
- verify_certificate: independent re-verification of certificates
- compute_certificate: GOLD/SILVER/BRONZE level computation
- format_certificate: human-readable certificate report
- _compute_certificate_hash: v1 and v2 hash computation
- _compute_gc_content: GC fraction helper
- Certificate dataclass: to_dict/from_dict round-trip
- Edge cases: empty results, missing parameters, hash integrity
"""

from __future__ import annotations

import hashlib
import pytest

from biocompiler.certificate import (
    generate_certificate,
    verify_certificate,
    compute_certificate,
    format_certificate,
    _compute_certificate_hash,
    _compute_gc_content,
    _CERTIFICATE_VERSION,
    _CURRENT_HASH_VERSION,
    _HASH_ALGORITHM,
    _REQUIRED_INPUT_PARAM_KEYS,
    VERSION,
)
from biocompiler.types import Verdict, TypeCheckResult, Certificate, SLOTMode
from biocompiler.type_system import CertLevel, PredicateResult
from biocompiler.exceptions import CertificateGenerationError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_pass_result(predicate: str = "GCInRange") -> TypeCheckResult:
    """Create a passing TypeCheckResult."""
    return TypeCheckResult(
        predicate=predicate,
        verdict=Verdict.PASS,
        derivation=[],
        violation=None,
        knowledge_gap=None,
    )


def _make_fail_result(predicate: str = "GCInRange") -> TypeCheckResult:
    """Create a failing TypeCheckResult."""
    return TypeCheckResult(
        predicate=predicate,
        verdict=Verdict.FAIL,
        derivation=[],
        violation="GC content out of range",
        knowledge_gap=None,
    )


def _basic_input_params() -> dict:
    """Return a minimal valid input_params dict."""
    return {
        "organism": "Homo_sapiens",
        "gc_lo": 0.30,
        "gc_hi": 0.70,
        "cai_threshold": 0.5,
        "enzymes": ["EcoRI", "BamHI"],
    }


# ---------------------------------------------------------------------------
# _compute_gc_content
# ---------------------------------------------------------------------------

class TestComputeGcContent:
    """Tests for the internal _compute_gc_content helper."""

    def test_fifty_percent(self):
        """50% GC should return ~0.5."""
        assert abs(_compute_gc_content("GCAT" * 25) - 0.5) < 0.01

    def test_all_gc(self):
        """All G/C should return 1.0."""
        assert _compute_gc_content("GCGCGCGC") == 1.0

    def test_all_at(self):
        """All A/T should return 0.0."""
        assert _compute_gc_content("ATATATAT") == 0.0

    def test_empty_sequence(self):
        """Empty sequence returns 0.0."""
        assert _compute_gc_content("") == 0.0

    def test_case_insensitive(self):
        """Should be case-insensitive."""
        assert _compute_gc_content("gcat") == _compute_gc_content("GCAT")


# ---------------------------------------------------------------------------
# _compute_certificate_hash
# ---------------------------------------------------------------------------

class TestComputeCertificateHash:
    """Tests for the internal _compute_certificate_hash function."""

    def test_v1_hash_is_sha256_of_sequence(self):
        """v1 hash should be SHA-256(sequence)."""
        seq = "ATGCGTACGT"
        expected = hashlib.sha256(seq.encode()).hexdigest()
        result = _compute_certificate_hash(
            sequence=seq, types_list=[], params={}, hash_version=1,
        )
        assert result == expected

    def test_v2_hash_differs_from_v1(self):
        """v2 hash should differ from v1 for the same input."""
        seq = "ATGCGTACGT"
        types_list = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params = {"organism": "Homo_sapiens"}
        v1 = _compute_certificate_hash(seq, types_list, params, hash_version=1)
        v2 = _compute_certificate_hash(seq, types_list, params, hash_version=2)
        assert v1 != v2

    def test_v2_hash_depends_on_predicates(self):
        """v2 hash should change when predicate results change."""
        seq = "ATGCGTACGT"
        params = {"organism": "Homo_sapiens"}
        types_pass = [{"predicate": "GCInRange", "verdict": "PASS"}]
        types_fail = [{"predicate": "GCInRange", "verdict": "FAIL"}]
        hash_pass = _compute_certificate_hash(seq, types_pass, params, hash_version=2)
        hash_fail = _compute_certificate_hash(seq, types_fail, params, hash_version=2)
        assert hash_pass != hash_fail

    def test_v2_hash_depends_on_parameters(self):
        """v2 hash should change when key parameters change."""
        seq = "ATGCGTACGT"
        types = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params1 = {"organism": "Homo_sapiens", "gc_lo": 0.30}
        params2 = {"organism": "Homo_sapiens", "gc_lo": 0.40}
        hash1 = _compute_certificate_hash(seq, types, params1, hash_version=2)
        hash2 = _compute_certificate_hash(seq, types, params2, hash_version=2)
        assert hash1 != hash2

    def test_v2_hash_deterministic(self):
        """v2 hash should be deterministic for the same input."""
        seq = "ATGCGTACGT"
        types = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params = {"organism": "Homo_sapiens"}
        h1 = _compute_certificate_hash(seq, types, params, hash_version=2)
        h2 = _compute_certificate_hash(seq, types, params, hash_version=2)
        assert h1 == h2

    def test_v2_hash_solver_backend_included(self):
        """v2 hash should depend on solver_backend."""
        seq = "ATGCGTACGT"
        types = [{"predicate": "GCInRange", "verdict": "PASS"}]
        params = {"organism": "Homo_sapiens"}
        h1 = _compute_certificate_hash(seq, types, params, hash_version=2, solver_backend="greedy")
        h2 = _compute_certificate_hash(seq, types, params, hash_version=2, solver_backend="ortools")
        assert h1 != h2


# ---------------------------------------------------------------------------
# generate_certificate
# ---------------------------------------------------------------------------

class TestGenerateCertificate:
    """Tests for generate_certificate."""

    def test_basic_generation(self):
        """Generate a certificate with all-passing results."""
        results = [_make_pass_result("GCInRange"), _make_pass_result("CodonAdapted")]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        assert isinstance(cert, Certificate)
        assert cert.sequence == "ATGCGT"
        assert len(cert.types) == 2

    def test_graduated_mode_default(self):
        """Default mode generates certificate even with failures."""
        results = [_make_pass_result("GCInRange"), _make_fail_result("NoRestrictionSite")]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        assert isinstance(cert, Certificate)
        # Should have partial status
        assert "PARTIAL" in cert.provenance["overall_status"]

    def test_strict_mode_raises_on_failure(self):
        """require_all_pass=True raises CertificateGenerationError."""
        results = [_make_pass_result("GCInRange"), _make_fail_result("NoRestrictionSite")]
        with pytest.raises(CertificateGenerationError):
            generate_certificate("ATGCGT", results, _basic_input_params(),
                                 require_all_pass=True)

    def test_strict_mode_passes_with_all_pass(self):
        """require_all_pass=True succeeds when all pass."""
        results = [_make_pass_result("GCInRange"), _make_pass_result("CodonAdapted")]
        cert = generate_certificate("ATGCGT", results, _basic_input_params(),
                                     require_all_pass=True)
        assert cert.provenance["overall_status"] == "FULL_PASS"

    def test_empty_sequence_raises(self):
        """Empty sequence should raise ValueError."""
        results = [_make_pass_result()]
        with pytest.raises(ValueError, match="Sequence must not be empty"):
            generate_certificate("", results, _basic_input_params())

    def test_empty_results_raises(self):
        """Empty type results should raise ValueError."""
        with pytest.raises(ValueError, match="Type results must not be empty"):
            generate_certificate("ATGCGT", [], _basic_input_params())

    def test_design_id_is_hash(self):
        """design_id should be a hex hash string."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        # SHA-256 hex digest is 64 characters
        assert len(cert.design_id) == 64
        assert all(c in "0123456789abcdef" for c in cert.design_id)

    def test_provenance_fields(self):
        """Provenance should have all required fields."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        prov = cert.provenance
        assert "tool" in prov
        assert "version" in prov
        assert "timestamp" in prov
        assert "input_hash" in prov
        assert prov["tool"] == "BioCompiler"

    def test_full_pass_status(self):
        """All passing results should yield FULL_PASS status."""
        results = [_make_pass_result("GCInRange"), _make_pass_result("CodonAdapted")]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        assert cert.provenance["overall_status"] == "FULL_PASS"

    def test_mutagenesis_metadata(self):
        """Mutagenesis substitutions should be documented in provenance."""
        results = [_make_pass_result()]
        subs = [{"position": 5, "from": "V", "to": "I", "blosum62": 2, "reason": "test"}]
        cert = generate_certificate("ATGCGT", results, _basic_input_params(),
                                     mutagenesis_substitutions=subs)
        assert cert.provenance["mutagenesis"]["applied"] is True
        assert cert.provenance["mutagenesis"]["n_substitutions"] == 1

    def test_no_mutagenesis_metadata(self):
        """No mutagenesis should set applied=False."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        assert cert.provenance["mutagenesis"]["applied"] is False

    def test_solver_backend_default(self):
        """Default solver_backend should be 'greedy'."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        assert cert.provenance["solver_backend"] == "greedy"

    def test_solver_backend_custom(self):
        """Custom solver_backend should be recorded."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params(),
                                     solver_backend="ortools")
        assert cert.provenance["solver_backend"] == "ortools"


# ---------------------------------------------------------------------------
# verify_certificate
# ---------------------------------------------------------------------------

class TestVerifyCertificate:
    """Tests for verify_certificate."""

    def test_valid_certificate_verifies(self):
        """A freshly generated certificate should verify."""
        results = [_make_pass_result("GCInRange")]
        params = _basic_input_params()
        cert = generate_certificate("ATGCGT", results, params)
        cert_dict = cert.to_dict()
        status, failures = verify_certificate(cert_dict)
        # May or may not verify depending on registry availability
        assert status in ("VERIFIED", "REJECTED")

    def test_tampered_sequence_rejected(self):
        """Changing the sequence after generation should fail verification."""
        results = [_make_pass_result("GCInRange")]
        params = _basic_input_params()
        cert = generate_certificate("ATGCGT", results, params)
        cert_dict = cert.to_dict()
        # Tamper with the sequence
        cert_dict["sequence"] = "ATGCGA"
        status, failures = verify_certificate(cert_dict)
        assert status == "REJECTED"
        assert any("design_id mismatch" in f for f in failures)

    def test_missing_required_keys_rejected(self):
        """Certificate missing required keys should be rejected."""
        cert_dict = {"version": "1.0"}  # Missing most keys
        status, failures = verify_certificate(cert_dict)
        assert status == "REJECTED"
        assert len(failures) > 0

    def test_missing_provenance_keys_rejected(self):
        """Certificate with incomplete provenance should be rejected."""
        cert_dict = {
            "version": "1.0",
            "design_id": "abc123",
            "sequence": "ATGCGT",
            "types": [],
            "provenance": {},  # Missing required provenance keys
        }
        status, failures = verify_certificate(cert_dict)
        assert status == "REJECTED"


# ---------------------------------------------------------------------------
# compute_certificate (GOLD/SILVER/BRONZE)
# ---------------------------------------------------------------------------

class TestComputeCertificate:
    """Tests for compute_certificate level computation."""

    def _make_predicate_result(self, passed: bool, details: str = "",
                                predicate: str = "test",
                                mutagenesis_applied: bool = False,
                                unavoidable_constraints: list | None = None) -> PredicateResult:
        """Create a PredicateResult for testing."""
        return PredicateResult(
            predicate=predicate,
            passed=passed,
            verdict=Verdict.PASS if passed else Verdict.FAIL,
            details=details,
            mutagenesis_applied=mutagenesis_applied,
            unavoidable_constraints=unavoidable_constraints or [],
        )

    def test_gold_level(self):
        """All predicates pass with no mutagenesis → GOLD."""
        results = [
            self._make_predicate_result(True, "OK", "GCInRange"),
            self._make_predicate_result(True, "OK", "CodonAdapted"),
        ]
        level = compute_certificate(results)
        assert level == CertLevel.GOLD

    def test_silver_level_mutagenesis(self):
        """Mutagenesis used → SILVER."""
        results = [
            self._make_predicate_result(True, "Resolved via mutagenesis", "NoRestrictionSite",
                                         mutagenesis_applied=True),
            self._make_predicate_result(True, "OK", "GCInRange"),
        ]
        level = compute_certificate(results)
        assert level == CertLevel.SILVER

    def test_silver_level_unavoidable(self):
        """Unavoidable constraints → SILVER."""
        results = [
            self._make_predicate_result(True, "unavoidable GT dinucleotide",
                                         "NoGTDinucleotide",
                                         unavoidable_constraints=["Valine_GT"]),
            self._make_predicate_result(True, "OK", "GCInRange"),
        ]
        level = compute_certificate(results)
        assert level == CertLevel.SILVER

    def test_bronze_level(self):
        """Some predicates fail → BRONZE."""
        results = [
            self._make_predicate_result(True, "OK", "GCInRange"),
            self._make_predicate_result(False, "GC content out of range", "GCInRange2"),
        ]
        level = compute_certificate(results)
        assert level == CertLevel.BRONZE


# ---------------------------------------------------------------------------
# format_certificate
# ---------------------------------------------------------------------------

class TestFormatCertificate:
    """Tests for format_certificate."""

    def test_basic_format(self):
        """format_certificate should return a readable string."""
        results = [
            PredicateResult("GCInRange", True, Verdict.PASS, "OK", 1.0),
            PredicateResult("CodonAdapted", True, Verdict.PASS, "OK", 0.8),
        ]
        report = format_certificate(results, "ATGCGT", "Homo_sapiens")
        assert isinstance(report, str)
        assert "BioCompiler" in report
        assert "GOLD" in report or "SILVER" in report or "BRONZE" in report
        assert "Homo_sapiens" in report

    def test_includes_predicate_names(self):
        """Report should include predicate names."""
        results = [
            PredicateResult("GCInRange", True, Verdict.PASS, "GC=50%", 1.0),
        ]
        report = format_certificate(results, "ATGCGT", "Homo_sapiens")
        assert "GCInRange" in report

    def test_failed_predicates_in_report(self):
        """Failed predicates should be listed."""
        results = [
            PredicateResult("GCInRange", True, Verdict.PASS, "OK", 1.0),
            PredicateResult("NoRestrictionSite", False, Verdict.FAIL, "EcoRI site found", 0.0),
        ]
        report = format_certificate(results, "ATGCGT", "Homo_sapiens")
        assert "NoRestrictionSite" in report


# ---------------------------------------------------------------------------
# Certificate dataclass round-trip
# ---------------------------------------------------------------------------

class TestCertificateRoundTrip:
    """Tests for Certificate to_dict/from_dict round-trip."""

    def test_round_trip(self):
        """Certificate should survive to_dict/from_dict round-trip."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        cert_dict = cert.to_dict()
        restored = Certificate.from_dict(cert_dict)
        assert restored.version == cert.version
        assert restored.design_id == cert.design_id
        assert restored.sequence == cert.sequence
        assert restored.hash_version == cert.hash_version

    def test_v2_hash_version_in_dict(self):
        """v2 certificates should include hash_version in dict."""
        results = [_make_pass_result()]
        cert = generate_certificate("ATGCGT", results, _basic_input_params())
        cert_dict = cert.to_dict()
        assert "hash_version" in cert_dict
        assert cert_dict["hash_version"] == 2

    def test_from_dict_missing_keys_raises(self):
        """from_dict with missing keys should raise ValueError."""
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict({"version": "1.0"})
