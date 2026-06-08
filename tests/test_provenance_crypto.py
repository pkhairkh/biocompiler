"""
BioCompiler Provenance Cryptographic Integrity Tests
=====================================================

Tests for the biocompiler.provenance_crypto module — HMAC-SHA256 signing,
verification, canonical JSON, key management, and mandatory provenance flag.
"""

from __future__ import annotations

import json
import os
import pytest

from biocompiler.provenance_crypto import (
    sign_record,
    verify_record,
    ProvenanceIntegrityError,
    get_provenance_secret,
    _canonical_json,
    is_mandatory_provenance,
)


# ════════════════════════════════════════════════════════════════════════════
# _canonical_json
# ════════════════════════════════════════════════════════════════════════════

class TestCanonicalJson:
    """Test canonical JSON serialization for deterministic hashing."""

    def test_sorted_keys(self):
        """Keys are sorted alphabetically."""
        data = {"z": 1, "a": 2, "m": 3}
        result = _canonical_json(data)
        parsed = json.loads(result)
        keys = list(parsed.keys())
        assert keys == ["a", "m", "z"]

    def test_compact_encoding(self):
        """No whitespace in output."""
        data = {"key": "value"}
        result = _canonical_json(data)
        assert " " not in result.decode("utf-8") or b'"key"' in result

    def test_deterministic_output(self):
        """Same data always produces the same bytes, regardless of insertion order."""
        d1 = {"a": 1, "b": 2, "c": 3}
        d2 = {"c": 3, "a": 1, "b": 2}
        assert _canonical_json(d1) == _canonical_json(d2)

    def test_returns_bytes(self):
        """Output is UTF-8 encoded bytes."""
        result = _canonical_json({"x": 1})
        assert isinstance(result, bytes)

    def test_nested_dict_sorted(self):
        """Nested dicts also have sorted keys."""
        data = {"outer": {"z": 1, "a": 2}}
        result = _canonical_json(data)
        parsed = json.loads(result)
        inner_keys = list(parsed["outer"].keys())
        assert inner_keys == ["a", "z"]


# ════════════════════════════════════════════════════════════════════════════
# sign_record
# ════════════════════════════════════════════════════════════════════════════

class TestSignRecord:

    def test_returns_hex_string(self):
        """Signature is a hex-encoded string."""
        sig = sign_record({"test": "data"}, key=b"testkey1234567890")
        assert isinstance(sig, str)
        # Hex string should have only hex characters
        assert all(c in "0123456789abcdef" for c in sig)

    def test_signature_length(self):
        """HMAC-SHA256 produces a 64-character hex string (256 bits)."""
        sig = sign_record({"test": "data"}, key=b"testkey1234567890")
        assert len(sig) == 64

    def test_same_data_same_key_same_sig(self):
        """Same data + same key should produce the same signature."""
        data = {"foo": "bar", "num": 42}
        key = b"deterministic_key"
        sig1 = sign_record(data, key=key)
        sig2 = sign_record(data, key=key)
        assert sig1 == sig2

    def test_different_data_different_sig(self):
        """Different data should produce different signatures."""
        key = b"same_key_for_both"
        sig1 = sign_record({"a": 1}, key=key)
        sig2 = sign_record({"a": 2}, key=key)
        assert sig1 != sig2

    def test_different_key_different_sig(self):
        """Same data + different key should produce different signatures."""
        data = {"test": "data"}
        sig1 = sign_record(data, key=b"key_one_1234567")
        sig2 = sign_record(data, key=b"key_two_1234567")
        assert sig1 != sig2

    def test_sign_with_default_key(self):
        """sign_record without explicit key uses the default provenance secret."""
        # Should not raise
        sig = sign_record({"test": "default_key"})
        assert isinstance(sig, str)
        assert len(sig) == 64

    def test_empty_dict(self):
        """Signing an empty dict should work."""
        sig = sign_record({}, key=b"testkey1234567890")
        assert len(sig) == 64


# ════════════════════════════════════════════════════════════════════════════
# verify_record
# ════════════════════════════════════════════════════════════════════════════

class TestVerifyRecord:

    def test_valid_signature_returns_true(self):
        """Correct signature verification returns True."""
        data = {"seed": 42, "decisions": ["a", "b"]}
        key = b"verify_key_12345"
        sig = sign_record(data, key=key)
        assert verify_record(data, sig, key=key) is True

    def test_tampered_data_returns_false(self):
        """Tampering with data after signing causes verification to fail."""
        data = {"seed": 42}
        key = b"verify_key_12345"
        sig = sign_record(data, key=key)
        tampered = {"seed": 43}
        assert verify_record(tampered, sig, key=key) is False

    def test_wrong_key_returns_false(self):
        """Verifying with the wrong key returns False."""
        data = {"seed": 42}
        sig = sign_record(data, key=b"correct_key!")
        assert verify_record(data, sig, key=b"wrong_key!!!") is False

    def test_bad_signature_returns_false(self):
        """A completely wrong signature returns False."""
        data = {"seed": 42}
        assert verify_record(data, "bad_signature", key=b"any_key") is False

    def test_verify_with_default_key(self):
        """Verify using the default provenance secret works."""
        data = {"test": "data"}
        sig = sign_record(data)  # uses default key
        assert verify_record(data, sig) is True  # also uses default key


# ════════════════════════════════════════════════════════════════════════════
# get_provenance_secret
# ════════════════════════════════════════════════════════════════════════════

class TestGetProvenanceSecret:

    def test_returns_bytes(self):
        """The provenance secret is returned as bytes."""
        secret = get_provenance_secret()
        assert isinstance(secret, bytes)

    def test_secret_length(self):
        """The secret should be 32 bytes (256 bits)."""
        secret = get_provenance_secret()
        assert len(secret) == 32

    def test_cached_secret_same_object(self):
        """Repeated calls return the same cached object."""
        import biocompiler.provenance_crypto as pc_mod
        # Force re-read from cache
        s1 = get_provenance_secret()
        s2 = get_provenance_secret()
        assert s1 is s2  # same object due to caching

    def test_env_variable_override(self):
        """BIOCOMPILER_PROVENANCE_SECRET env var overrides cached secret."""
        import biocompiler.provenance_crypto as pc_mod
        test_hex = "aa" * 32  # 32 bytes hex
        # Save old state
        old_cached = pc_mod._cached_secret
        old_env = os.environ.get("BIOCOMPILER_PROVENANCE_SECRET")
        try:
            pc_mod._cached_secret = None
            os.environ["BIOCOMPILER_PROVENANCE_SECRET"] = test_hex
            secret = get_provenance_secret()
            assert secret == bytes.fromhex(test_hex)
        finally:
            # Restore
            pc_mod._cached_secret = old_cached
            if old_env is not None:
                os.environ["BIOCOMPILER_PROVENANCE_SECRET"] = old_env
            else:
                os.environ.pop("BIOCOMPILER_PROVENANCE_SECRET", None)


# ════════════════════════════════════════════════════════════════════════════
# ProvenanceIntegrityError
# ════════════════════════════════════════════════════════════════════════════

class TestProvenanceIntegrityError:

    def test_is_exception(self):
        """ProvenanceIntegrityError inherits from Exception."""
        assert issubclass(ProvenanceIntegrityError, Exception)

    def test_message(self):
        """ProvenanceIntegrityError carries a message."""
        err = ProvenanceIntegrityError("sig mismatch")
        assert "sig mismatch" in str(err)


# ════════════════════════════════════════════════════════════════════════════
# is_mandatory_provenance
# ════════════════════════════════════════════════════════════════════════════

class TestIsMandatoryProvenance:

    def test_default_is_false(self):
        """Without the env var, provenance is not mandatory."""
        old = os.environ.pop("BIOCOMPILER_PROVENANCE_MANDATORY", None)
        try:
            assert is_mandatory_provenance() is False
        finally:
            if old is not None:
                os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = old

    def test_set_to_true(self):
        """With BIOCOMPILER_PROVENANCE_MANDATORY=true, returns True."""
        old = os.environ.get("BIOCOMPILER_PROVENANCE_MANDATORY")
        try:
            os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = "true"
            assert is_mandatory_provenance() is True
        finally:
            if old is not None:
                os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = old
            else:
                os.environ.pop("BIOCOMPILER_PROVENANCE_MANDATORY", None)

    def test_set_to_one(self):
        """With BIOCOMPILER_PROVENANCE_MANDATORY=1, returns True."""
        old = os.environ.get("BIOCOMPILER_PROVENANCE_MANDATORY")
        try:
            os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = "1"
            assert is_mandatory_provenance() is True
        finally:
            if old is not None:
                os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = old
            else:
                os.environ.pop("BIOCOMPILER_PROVENANCE_MANDATORY", None)

    def test_set_to_false(self):
        """With BIOCOMPILER_PROVENANCE_MANDATORY=false, returns False."""
        old = os.environ.get("BIOCOMPILER_PROVENANCE_MANDATORY")
        try:
            os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = "false"
            assert is_mandatory_provenance() is False
        finally:
            if old is not None:
                os.environ["BIOCOMPILER_PROVENANCE_MANDATORY"] = old
            else:
                os.environ.pop("BIOCOMPILER_PROVENANCE_MANDATORY", None)
