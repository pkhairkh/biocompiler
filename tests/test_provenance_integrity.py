"""
Tests for provenance integrity and certificate security (Task 1.9).

Covers:
1. HMAC-SHA256 signing of provenance records
2. Signature verification on load
3. Tampering detection
4. Mandatory provenance mode (BIOCOMPILER_PROVENANCE_MANDATORY)
5. Certificate hash covers predicate results + organism + params
6. certificate_version field in generated certificates
7. Certificate level determination via structured flags
8. Backward-compatible fallback to string matching
9. input_params validation and organism default warning
"""

import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import pytest

from biocompiler.provenance.certificate import (
    compute_certificate,
    generate_certificate,
    verify_certificate,
    _CERTIFICATE_VERSION,
    _REQUIRED_INPUT_PARAM_KEYS,
)
from biocompiler.provenance.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
    ProvenanceStore,
)
from biocompiler.provenance.crypto import (
    sign_record,
    verify_record,
    ProvenanceIntegrityError,
    get_provenance_secret,
    is_mandatory_provenance,
)
from biocompiler.type_system import CertLevel, PredicateResult, SpliceVerdict
from biocompiler.shared.types import Certificate, TypeCheckResult, Verdict, SLOTMode


# ─── Shared helpers ─────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"


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


def _make_predicate_results(
    all_pass: bool = True,
    mutagenesis_flag: bool = False,
    unavoidable_list: list[str] | None = None,
    mutagenesis_in_details: bool = False,
    unavoidable_in_details: bool = False,
) -> list[PredicateResult]:
    """Build a list of PredicateResult objects with structured flags."""
    results = [
        PredicateResult("NoStopCodons", True, verdict=Verdict.PASS, details="No internal stop codons"),
        PredicateResult(
            "NoCrypticSplice", True, verdict=SpliceVerdict.PASS,
            details="No GT dinucleotides found",
        ),
        PredicateResult("NoCpGIsland", True, details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
        PredicateResult("NoRestrictionSite", True, details="No restriction sites found"),
    ]

    # Build the GT predicate with structured flags
    gt_details = "No GT dinucleotides found"
    if unavoidable_in_details:
        gt_details = "All 2 GT dinucleotides are unavoidable"
    elif mutagenesis_in_details:
        gt_details = "No GT dinucleotides found mutagenesis applied: pos 3:V→I"

    results.append(PredicateResult(
        "NoGTDinucleotide", True,
        details=gt_details,
        mutagenesis_applied=mutagenesis_flag,
        unavoidable_constraints=unavoidable_list or [],
    ))

    if not all_pass:
        results.append(PredicateResult(
            "ValidCodingSeq", False,
            details="Sequence length not divisible by 3",
        ))
    else:
        results.append(PredicateResult("ValidCodingSeq", True, details="All codons valid"))

    results.append(PredicateResult(
        "ConservationScore", True, details="All AA conservation scores >= -1",
    ))
    results.append(PredicateResult(
        "CodonOptimality", True, details="Worst CAI: GCT=0.7244, min=0.0",
    ))
    return results


def _make_trail() -> OptimizationDecisionTrail:
    """Build a simple OptimizationDecisionTrail for testing."""
    return OptimizationDecisionTrail(
        gene_name="test_gene",
        input_protein="MVLSPADKTN",
        output_dna="ATGGTGCTGCCTGCTCCTGCTG",
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
                constraint_name="NoCrypticSplice", constraint_type="hard",
                action_taken="satisfied", positions_affected=[3, 7],
                tradeoff_description="Chose GTC over GTG to avoid cryptic donor",
                impact_on_cai=-0.003,
            ),
        ],
        iteration_log=[{"step": 1, "score": 0.85}],
        timestamp=datetime.now(timezone.utc).isoformat(),
        version="1.0.0",
    )


# ═══════════════════════════════════════════════════════════════════
# 1. HMAC-SHA256 Signing of Provenance Records
# ═══════════════════════════════════════════════════════════════════

class TestHMACSigning:
    """Tests for provenance_crypto sign_record and verify_record."""

    def test_sign_record_returns_hex_string(self):
        """sign_record should return a hex-encoded HMAC-SHA256 string."""
        data = {"seed": 42, "decisions": []}
        sig = sign_record(data)
        assert isinstance(sig, str)
        # HMAC-SHA256 produces 32 bytes = 64 hex chars
        assert len(sig) == 64
        assert all(c in "0123456789abcdef" for c in sig)

    def test_sign_record_deterministic(self):
        """Signing the same data with the same key should produce the same signature."""
        data = {"seed": 42, "decisions": [{"position": 0, "chosen": "ATG"}]}
        key = b"test_key_12345678"
        sig1 = sign_record(data, key=key)
        sig2 = sign_record(data, key=key)
        assert sig1 == sig2

    def test_sign_record_key_sensitive(self):
        """Different keys should produce different signatures."""
        data = {"seed": 42}
        sig1 = sign_record(data, key=b"key_a_1234567890")
        sig2 = sign_record(data, key=b"key_b_1234567890")
        assert sig1 != sig2

    def test_sign_record_data_sensitive(self):
        """Different data should produce different signatures."""
        key = b"test_key_12345678"
        sig1 = sign_record({"seed": 42}, key=key)
        sig2 = sign_record({"seed": 43}, key=key)
        assert sig1 != sig2

    def test_verify_record_valid_signature(self):
        """verify_record should return True for a valid signature."""
        data = {"seed": 42, "decisions": []}
        sig = sign_record(data)
        assert verify_record(data, sig) is True

    def test_verify_record_tampered_data(self):
        """verify_record should return False if data was tampered with."""
        data = {"seed": 42, "decisions": []}
        sig = sign_record(data)
        tampered = {"seed": 43, "decisions": []}
        assert verify_record(tampered, sig) is False

    def test_verify_record_wrong_signature(self):
        """verify_record should return False for a wrong signature."""
        data = {"seed": 42}
        assert verify_record(data, "0" * 64) is False

    def test_sign_record_order_independent(self):
        """Signing should be deterministic regardless of dict key order."""
        data1 = {"a": 1, "b": 2}
        data2 = {"b": 2, "a": 1}
        key = b"test_key_12345678"
        sig1 = sign_record(data1, key=key)
        sig2 = sign_record(data2, key=key)
        assert sig1 == sig2

    def test_provenance_integrity_error_is_exception(self):
        """ProvenanceIntegrityError should be a proper exception class."""
        assert issubclass(ProvenanceIntegrityError, Exception)
        err = ProvenanceIntegrityError("test error")
        assert str(err) == "test error"


class TestProvenanceSecretLoading:
    """Tests for get_provenance_secret key loading."""

    def test_env_var_takes_priority(self):
        """BIOCOMPILER_PROVENANCE_SECRET env var should take priority."""
        import biocompiler.provenance.crypto as crypto_mod
        # Reset cache
        crypto_mod._cached_secret = None
        test_key = b"\x01" * 32
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_SECRET": test_key.hex()}):
            key = get_provenance_secret()
            assert key == test_key
        # Clean up
        crypto_mod._cached_secret = None

    def test_auto_generated_key_is_bytes(self):
        """When no env var is set, a key should be auto-generated as bytes."""
        import biocompiler.provenance.crypto as crypto_mod
        crypto_mod._cached_secret = None
        with patch.dict(os.environ, {}, clear=False):
            # Remove the env var if it exists
            os.environ.pop("BIOCOMPILER_PROVENANCE_SECRET", None)
            key = get_provenance_secret()
            assert isinstance(key, bytes)
            assert len(key) == 32
        crypto_mod._cached_secret = None


# ═══════════════════════════════════════════════════════════════════
# 2. ProvenanceStore HMAC Integration
# ═══════════════════════════════════════════════════════════════════

class TestProvenanceStoreHMAC:
    """Tests for HMAC signing in ProvenanceStore save/load."""

    def test_save_includes_hmac_signature(self):
        """Saved records should include an _hmac_signature field."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)

            filepath = Path(tmpdir) / f"{record_id}.json"
            with open(filepath) as f:
                data = json.load(f)
            assert "_hmac_signature" in data
            assert len(data["_hmac_signature"]) == 64

    def test_load_verifies_signature(self):
        """Loading a valid record should succeed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)
            loaded = store.load(record_id)
            assert loaded.organism == "Homo_sapiens"
            assert len(loaded.codon_decisions) == 1

    def test_load_detects_tampering(self):
        """Loading a tampered record should raise ProvenanceIntegrityError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)

            # Tamper with the file
            filepath = Path(tmpdir) / f"{record_id}.json"
            with open(filepath) as f:
                data = json.load(f)
            data["organism"] = "TAMPERED_ORGANISM"
            with open(filepath, "w") as f:
                json.dump(data, f)

            with pytest.raises(ProvenanceIntegrityError, match="HMAC signature mismatch"):
                store.load(record_id)

    def test_load_missing_signature_raises(self):
        """Loading a record without an HMAC signature should raise ProvenanceIntegrityError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)

            # Remove the signature
            filepath = Path(tmpdir) / f"{record_id}.json"
            with open(filepath) as f:
                data = json.load(f)
            del data["_hmac_signature"]
            with open(filepath, "w") as f:
                json.dump(data, f)

            with pytest.raises(ProvenanceIntegrityError, match="no HMAC signature"):
                store.load(record_id)

    def test_save_load_roundtrip(self):
        """Save + load should produce an equivalent trail."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)
            loaded = store.load(record_id)

            assert loaded.gene_name == trail.gene_name
            assert loaded.input_protein == trail.input_protein
            assert loaded.output_dna == trail.output_dna
            assert loaded.organism == trail.organism
            assert loaded.solver_backend == trail.solver_backend
            assert loaded.seed == trail.seed
            assert len(loaded.codon_decisions) == len(trail.codon_decisions)
            assert len(loaded.constraint_decisions) == len(trail.constraint_decisions)

    def test_query_skips_tampered_records(self):
        """query() should skip records with invalid HMAC signatures."""
        with tempfile.TemporaryDirectory() as tmpdir:
            store = ProvenanceStore(store_dir=tmpdir)
            trail = _make_trail()
            record_id = store.save(trail)

            # Tamper with the file
            filepath = Path(tmpdir) / f"{record_id}.json"
            with open(filepath) as f:
                data = json.load(f)
            data["organism"] = "TAMPERED"
            with open(filepath, "w") as f:
                json.dump(data, f)

            results = store.query(organism="TAMPERED")
            # The tampered record should be skipped
            assert len(results) == 0


# ═══════════════════════════════════════════════════════════════════
# 3. Mandatory Provenance Mode
# ═══════════════════════════════════════════════════════════════════

class TestMandatoryProvenance:
    """Tests for BIOCOMPILER_PROVENANCE_MANDATORY env var."""

    def test_default_is_not_mandatory(self):
        """By default, provenance is not mandatory."""
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("BIOCOMPILER_PROVENANCE_MANDATORY", None)
            assert is_mandatory_provenance() is False

    def test_set_to_true_is_mandatory(self):
        """BIOCOMPILER_PROVENANCE_MANDATORY=true should make it mandatory."""
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_MANDATORY": "true"}):
            assert is_mandatory_provenance() is True

    def test_set_to_1_is_mandatory(self):
        """BIOCOMPILER_PROVENANCE_MANDATORY=1 should make it mandatory."""
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_MANDATORY": "1"}):
            assert is_mandatory_provenance() is True

    def test_set_to_yes_is_mandatory(self):
        """BIOCOMPILER_PROVENANCE_MANDATORY=yes should make it mandatory."""
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_MANDATORY": "yes"}):
            assert is_mandatory_provenance() is True

    def test_set_to_false_is_not_mandatory(self):
        """BIOCOMPILER_PROVENANCE_MANDATORY=false should not be mandatory."""
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_MANDATORY": "false"}):
            assert is_mandatory_provenance() is False

    def test_set_to_0_is_not_mandatory(self):
        """BIOCOMPILER_PROVENANCE_MANDATORY=0 should not be mandatory."""
        with patch.dict(os.environ, {"BIOCOMPILER_PROVENANCE_MANDATORY": "0"}):
            assert is_mandatory_provenance() is False


# ═══════════════════════════════════════════════════════════════════
# 4. Certificate Hash Covers Predicate Results + Params
# ═══════════════════════════════════════════════════════════════════

class TestCertificateHash:
    """Tests for the v2 certificate hash format."""

    def test_design_id_changes_with_different_predicates(self):
        """Two certificates with different predicate results should have different design_ids."""
        type_results_1 = _make_type_results([("GCInRange", Verdict.PASS)])
        type_results_2 = _make_type_results([("GCInRange", Verdict.FAIL)])

        cert1 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results_1,
            input_params={"organism": "Homo_sapiens"},
        )
        cert2 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results_2,
            input_params={"organism": "Homo_sapiens"},
        )
        assert cert1.design_id != cert2.design_id

    def test_design_id_changes_with_different_organism(self):
        """Two certificates with different organisms should have different design_ids."""
        type_results = _make_type_results()

        cert1 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
        )
        cert2 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"organism": "Escherichia_coli"},
        )
        assert cert1.design_id != cert2.design_id

    def test_design_id_changes_with_different_gc_params(self):
        """Two certificates with different GC bounds should have different design_ids."""
        type_results = _make_type_results()

        cert1 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens", "gc_lo": 0.30, "gc_hi": 0.70},
        )
        cert2 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens", "gc_lo": 0.40, "gc_hi": 0.60},
        )
        assert cert1.design_id != cert2.design_id

    def test_design_id_changes_with_different_solver(self):
        """Two certificates with different solver_backends should have different design_ids."""
        type_results = _make_type_results()

        cert1 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
            solver_backend="greedy",
        )
        cert2 = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"organism": "Homo_sapiens"},
            solver_backend="ortools",
        )
        assert cert1.design_id != cert2.design_id

    def test_same_input_same_hash(self):
        """Same sequence + predicates + params should produce the same hash."""
        type_results = _make_type_results()
        params = {"organism": "Homo_sapiens", "gc_lo": 0.30}

        cert1 = generate_certificate(SAMPLE_SEQ, type_results, params, solver_backend="greedy")
        cert2 = generate_certificate(SAMPLE_SEQ, type_results, params, solver_backend="greedy")
        assert cert1.design_id == cert2.design_id


# ═══════════════════════════════════════════════════════════════════
# 5. Certificate Version Field
# ═══════════════════════════════════════════════════════════════════

class TestCertificateVersion:
    """Tests for the certificate_version field."""

    def test_certificate_version_in_provenance(self):
        """Generated certificate provenance should include certificate_version."""
        type_results = _make_type_results()
        cert = generate_certificate(SAMPLE_SEQ, type_results, {"organism": "Homo_sapiens"})
        assert "certificate_version" in cert.provenance
        assert cert.provenance["certificate_version"] == _CERTIFICATE_VERSION

    def test_certificate_version_is_2(self):
        """The current certificate version should be 2."""
        assert _CERTIFICATE_VERSION == 2


# ═══════════════════════════════════════════════════════════════════
# 6. Certificate Level — Structured Flags
# ═══════════════════════════════════════════════════════════════════

class TestCertificateLevelStructuredFlags:
    """Tests for certificate level computation using structured flags."""

    def test_gold_with_structured_flags(self):
        """All pass, no mutagenesis_applied, no unavoidable_constraints → GOLD."""
        results = _make_predicate_results(all_pass=True)
        assert compute_certificate(results) == CertLevel.GOLD

    def test_silver_with_mutagenesis_flag(self):
        """All pass with mutagenesis_applied=True → SILVER."""
        results = _make_predicate_results(all_pass=True, mutagenesis_flag=True)
        assert compute_certificate(results) == CertLevel.SILVER

    def test_silver_with_unavoidable_flag(self):
        """All pass with unavoidable_constraints=['GT_dinucleotide'] → SILVER."""
        results = _make_predicate_results(
            all_pass=True,
            unavoidable_list=["GT_dinucleotide_Valine"],
        )
        assert compute_certificate(results) == CertLevel.SILVER

    def test_bronze_with_failed_predicate(self):
        """Any failed predicate → BRONZE regardless of flags."""
        results = _make_predicate_results(
            all_pass=False,
            mutagenesis_flag=True,
        )
        assert compute_certificate(results) == CertLevel.BRONZE

    def test_structured_flags_override_string_matching(self):
        """Structured flags should be checked even if details do not contain keywords."""
        results = [
            PredicateResult(
                "NoGTDinucleotide", True,
                details="No GT dinucleotides found",  # No "mutagenesis" or "unavoidable" in text
                mutagenesis_applied=True,  # But the flag is set
            ),
        ]
        assert compute_certificate(results) == CertLevel.SILVER

    def test_backward_compat_string_matching_mutagenesis(self):
        """If mutagenesis_applied is not set but details mentions 'mutagenesis', still SILVER."""
        results = [
            PredicateResult(
                "NoGTDinucleotide", True,
                details="No GT dinucleotides found mutagenesis applied: pos 3:V→I",
            ),
        ]
        assert compute_certificate(results) == CertLevel.SILVER

    def test_backward_compat_string_matching_unavoidable(self):
        """If unavoidable_constraints is not set but details mentions 'unavoidable', still SILVER."""
        results = [
            PredicateResult(
                "NoGTDinucleotide", True,
                details="All 2 GT dinucleotides are unavoidable",
            ),
        ]
        assert compute_certificate(results) == CertLevel.SILVER

    def test_no_false_positive_from_details(self):
        """Details without 'mutagenesis' or 'unavoidable' and no flags → GOLD."""
        results = [
            PredicateResult(
                "NoGTDinucleotide", True,
                details="No GT dinucleotides found",
                mutagenesis_applied=False,
                unavoidable_constraints=[],
            ),
        ]
        assert compute_certificate(results) == CertLevel.GOLD


# ═══════════════════════════════════════════════════════════════════
# 7. Certificate Default Parameters — Validation & Warning
# ═══════════════════════════════════════════════════════════════════

class TestCertificateDefaultParameters:
    """Tests for input_params validation and organism default behavior."""

    def test_missing_required_keys_warns(self, caplog):
        """Missing required input_params keys should trigger a warning."""
        type_results = _make_type_results()
        with caplog.at_level(logging.WARNING):
            cert = generate_certificate(
                SAMPLE_SEQ, type_results, input_params={}
            )
        # Should have warning about missing keys
        assert any("missing required keys" in msg.lower() for msg in caplog.messages)

    def test_all_required_keys_no_warning(self, caplog):
        """All required keys present should not trigger a warning about missing keys."""
        type_results = _make_type_results()
        params = {
            "organism": "Escherichia_coli",
            "gc_lo": 0.30,
            "gc_hi": 0.70,
            "cai_threshold": 0.5,
            "enzymes": ["EcoRI"],
        }
        with caplog.at_level(logging.WARNING):
            cert = generate_certificate(
                SAMPLE_SEQ, type_results, input_params=params
            )
        # No warning about missing keys
        assert not any("missing required keys" in msg.lower() for msg in caplog.messages)

    def test_missing_organism_warns_about_default(self, caplog):
        """Missing organism should trigger a specific warning about defaulting."""
        type_results = _make_type_results()
        with caplog.at_level(logging.WARNING):
            cert = generate_certificate(
                SAMPLE_SEQ, type_results, input_params={"gc_lo": 0.3}
            )
        assert any("organism" in msg.lower() and "default" in msg.lower()
                    for msg in caplog.messages)

    def test_organism_preserved_when_provided(self):
        """When organism is provided, it should be used, not defaulted."""
        type_results = _make_type_results()
        cert = generate_certificate(
            SAMPLE_SEQ, type_results,
            input_params={"organism": "Escherichia_coli"},
        )
        assert cert.provenance["parameters"]["organism"] == "Escherichia_coli"

    def test_default_organism_is_homo_sapiens(self):
        """When organism is not provided, it defaults to Homo_sapiens (with warning)."""
        type_results = _make_type_results()
        cert = generate_certificate(
            SAMPLE_SEQ, type_results, input_params={}
        )
        assert cert.provenance["parameters"]["organism"] == "Homo_sapiens"

    def test_design_id_reflects_actual_organism(self):
        """The design_id hash should reflect the actual organism used."""
        type_results = _make_type_results()

        cert_human = generate_certificate(
            SAMPLE_SEQ, type_results,
            input_params={"organism": "Homo_sapiens"},
        )
        cert_ecoli = generate_certificate(
            SAMPLE_SEQ, type_results,
            input_params={"organism": "Escherichia_coli"},
        )
        # Different organisms → different hashes
        assert cert_human.design_id != cert_ecoli.design_id


# ═══════════════════════════════════════════════════════════════════
# 8. PredicateResult Structured Fields
# ═══════════════════════════════════════════════════════════════════

class TestPredicateResultFields:
    """Tests for new PredicateResult fields: mutagenesis_applied, unavoidable_constraints."""

    def test_default_mutagenesis_applied_is_false(self):
        """mutagenesis_applied should default to False."""
        r = PredicateResult("Test", True, details="ok")
        assert r.mutagenesis_applied is False

    def test_default_unavoidable_constraints_is_empty(self):
        """unavoidable_constraints should default to empty list."""
        r = PredicateResult("Test", True, details="ok")
        assert r.unavoidable_constraints == []

    def test_mutagenesis_applied_can_be_set(self):
        """mutagenesis_applied can be explicitly set to True."""
        r = PredicateResult("Test", True, details="ok", mutagenesis_applied=True)
        assert r.mutagenesis_applied is True

    def test_unavoidable_constraints_can_be_set(self):
        """unavoidable_constraints can be set to a list of constraint names."""
        r = PredicateResult(
            "Test", True, details="ok",
            unavoidable_constraints=["GT_dinucleotide_Valine"],
        )
        assert r.unavoidable_constraints == ["GT_dinucleotide_Valine"]

    def test_fields_do_not_affect_basic_construction(self):
        """Old code that does not pass the new fields should still work."""
        r = PredicateResult("Test", True, verdict=Verdict.PASS, details="ok")
        assert r.predicate == "Test"
        assert r.passed is True
        assert r.mutagenesis_applied is False
        assert r.unavoidable_constraints == []


# ═══════════════════════════════════════════════════════════════════
# 9. Verify Certificate with v2 Hash
# ═══════════════════════════════════════════════════════════════════

class TestCertificateVerification:
    """Tests for certificate verification with the v2 hash format."""

    def test_v2_certificate_verifies_design_id(self):
        """A v2 certificate should verify with the correct design_id."""
        type_results = _make_type_results()
        cert = generate_certificate(
            SAMPLE_SEQ, type_results,
            input_params={"organism": "Homo_sapiens", "gc_lo": 0.30, "gc_hi": 0.70},
        )
        cert_dict = cert.to_dict()
        # v2 certificates have certificate_version in provenance
        assert cert_dict["provenance"].get("certificate_version") == 2
        # The design_id should match the computed v2 hash
        status, failures = verify_certificate(cert_dict)
        # Note: verification may fail due to registry re-evaluation,
        # but design_id mismatch should NOT be in the failures
        design_id_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(design_id_failures) == 0

    def test_tampered_sequence_detected(self):
        """Changing the sequence after signing should cause design_id mismatch."""
        type_results = _make_type_results()
        cert = generate_certificate(
            SAMPLE_SEQ, type_results,
            input_params={"organism": "Homo_sapiens"},
        )
        cert_dict = cert.to_dict()
        # Tamper with the sequence
        cert_dict["sequence"] = "ATGCATGCATGCATGC"
        status, failures = verify_certificate(cert_dict)
        design_id_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(design_id_failures) > 0

    def test_tampered_predicate_results_detected(self):
        """Changing predicate results after signing should cause design_id mismatch."""
        type_results = _make_type_results([("GCInRange", Verdict.PASS)])
        cert = generate_certificate(
            SAMPLE_SEQ, type_results,
            input_params={"organism": "Homo_sapiens"},
        )
        cert_dict = cert.to_dict()
        # Tamper with predicate results
        cert_dict["types"][0]["verdict"] = "FAIL"
        status, failures = verify_certificate(cert_dict)
        design_id_failures = [f for f in failures if "design_id mismatch" in f]
        assert len(design_id_failures) > 0


# Need logging import for caplog tests
import logging
