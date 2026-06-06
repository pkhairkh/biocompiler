"""
BioCompiler Provenance Cryptographic Integrity — HMAC-SHA256 Signing
====================================================================

Provides HMAC-SHA256 signing and verification for provenance records,
ensuring that saved provenance data cannot be tampered with without
detection.

The signing key is loaded from the ``BIOCOMPILER_PROVENANCE_SECRET``
environment variable. If not set, a key is auto-generated on first use
and persisted to ``~/.biocompiler/provenance_secret`` for reuse across
restarts.

Usage::

    from biocompiler.provenance_crypto import sign_record, verify_record

    # Sign a provenance dict before saving
    data = {"seed": 42, "decisions": [...]}
    signature = sign_record(data)
    data["_hmac_signature"] = signature

    # Verify on load
    stored_sig = data.pop("_hmac_signature")
    verify_record(data, stored_sig)  # raises ProvenanceIntegrityError on mismatch
"""

from __future__ import annotations

import hashlib
import hmac
import json
import logging
import os
import secrets
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "sign_record",
    "verify_record",
    "ProvenanceIntegrityError",
    "get_provenance_secret",
]

# Environment variable for the provenance signing key
_PROVENANCE_SECRET_ENV = "BIOCOMPILER_PROVENANCE_SECRET"

# Path where auto-generated key is persisted
_PROVENANCE_SECRET_FILE = Path.home() / ".biocompiler" / "provenance_secret"

# Cached signing key (loaded once, reused)
_cached_secret: bytes | None = None


class ProvenanceIntegrityError(Exception):
    """Raised when a provenance record's HMAC signature fails verification."""
    pass


def _canonical_json(data: dict[str, Any]) -> bytes:
    """Produce a canonical JSON byte representation for hashing.

    Sorts keys recursively and uses compact encoding (no whitespace)
    to ensure deterministic output regardless of dict insertion order.

    Args:
        data: The dict to canonicalize.

    Returns:
        UTF-8 encoded canonical JSON bytes.
    """
    return json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")


def get_provenance_secret() -> bytes:
    """Get or generate the HMAC signing key.

    Resolution order:
    1. ``BIOCOMPILER_PROVENANCE_SECRET`` environment variable (hex-encoded)
    2. Persisted key at ``~/.biocompiler/provenance_secret``
    3. Auto-generate a new 32-byte key and persist it

    Returns:
        The signing key as bytes.

    Note:
        The key is cached in memory after first load for performance.
    """
    global _cached_secret
    if _cached_secret is not None:
        return _cached_secret

    # 1. Check environment variable
    env_secret = os.environ.get(_PROVENANCE_SECRET_ENV, "")
    if env_secret:
        try:
            _cached_secret = bytes.fromhex(env_secret)
            logger.debug("Loaded provenance secret from environment variable")
            return _cached_secret
        except ValueError:
            logger.warning(
                "Invalid BIOCOMPILER_PROVENANCE_SECRET (not valid hex), "
                "falling back to file or auto-generation"
            )

    # 2. Check persisted file
    if _PROVENANCE_SECRET_FILE.exists():
        try:
            file_secret = _PROVENANCE_SECRET_FILE.read_text().strip()
            if file_secret:
                _cached_secret = bytes.fromhex(file_secret)
                logger.debug("Loaded provenance secret from %s", _PROVENANCE_SECRET_FILE)
                return _cached_secret
        except (OSError, ValueError) as exc:
            logger.warning(
                "Could not read provenance secret from %s: %s",
                _PROVENANCE_SECRET_FILE, exc,
            )

    # 3. Auto-generate and persist
    _cached_secret = secrets.token_bytes(32)
    try:
        _PROVENANCE_SECRET_FILE.parent.mkdir(parents=True, exist_ok=True)
        _PROVENANCE_SECRET_FILE.write_text(_cached_secret.hex())
        _PROVENANCE_SECRET_FILE.chmod(0o600)
        logger.info(
            "Generated new provenance secret and saved to %s",
            _PROVENANCE_SECRET_FILE,
        )
    except OSError as exc:
        logger.warning(
            "Could not persist provenance secret to %s: %s. "
            "Key will not survive restart.",
            _PROVENANCE_SECRET_FILE, exc,
        )

    return _cached_secret


def sign_record(data: dict[str, Any], key: bytes | None = None) -> str:
    """Compute an HMAC-SHA256 signature over a provenance record dict.

    The signature is computed over the canonical JSON representation
    of the data dict (sorted keys, compact encoding). This ensures
    deterministic signing regardless of dict insertion order.

    Args:
        data: The provenance record dict to sign. Must be JSON-serializable.
        key: The HMAC key. If None, uses the default provenance secret.

    Returns:
        Hex-encoded HMAC-SHA256 signature string.
    """
    if key is None:
        key = get_provenance_secret()
    canonical = _canonical_json(data)
    return hmac.new(key, canonical, hashlib.sha256).hexdigest()


def verify_record(
    data: dict[str, Any],
    expected_signature: str,
    key: bytes | None = None,
) -> bool:
    """Verify the HMAC-SHA256 signature of a provenance record.

    Uses constant-time comparison to prevent timing attacks.

    Args:
        data: The provenance record dict to verify.
        expected_signature: The hex-encoded signature to match.
        key: The HMAC key. If None, uses the default provenance secret.

    Returns:
        True if the signature matches, False otherwise.
    """
    if key is None:
        key = get_provenance_secret()
    computed = sign_record(data, key=key)
    return hmac.compare_digest(computed, expected_signature)


def is_mandatory_provenance() -> bool:
    """Check whether provenance tracking is mandatory.

    When ``BIOCOMPILER_PROVENANCE_MANDATORY`` is set (to any truthy value),
    per-request ``track_provenance=False`` is overridden to force provenance ON.

    Returns:
        True if provenance is mandatory, False otherwise.
    """
    return os.environ.get("BIOCOMPILER_PROVENANCE_MANDATORY", "").lower() in (
        "1", "true", "yes", "on",
    )
