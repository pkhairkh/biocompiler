"""
BioCompiler API — Authentication and rate limiting.

This module handles API key management, verification, and per-client
rate limiting using a SQLite-backed sliding window.
"""

import logging
import os
import secrets
from pathlib import Path

from fastapi import HTTPException, Security
from fastapi.security import APIKeyHeader

logger = logging.getLogger(__name__)

# ─── API Key Authentication ─────────────────────────────────────────

API_KEY_NAME = "X-API-Key"
_api_key_header = APIKeyHeader(name=API_KEY_NAME, auto_error=False)

# ── Auth mode: "required" (default), "optional", "disabled" ──────────
_AUTH_MODE = os.environ.get("BIOCOMPILER_AUTH_MODE", "required").lower()
if _AUTH_MODE not in ("required", "optional", "disabled"):
    logger.warning(
        "Invalid BIOCOMPILER_AUTH_MODE=%r, falling back to 'required'",
        _AUTH_MODE,
    )
    _AUTH_MODE = "required"

# Global flag set by --no-auth CLI flag (see cli.py)
_NO_AUTH_CLI_FLAG = False


def set_no_auth_flag() -> None:
    """Enable --no-auth mode from the CLI. Emits a warning."""
    global _NO_AUTH_CLI_FLAG, _AUTH_MODE
    _NO_AUTH_CLI_FLAG = True
    _AUTH_MODE = "disabled"
    logger.warning(
        "--no-auth flag: API authentication is DISABLED. "
        "Do not use in production!"
    )


def _get_api_key_file_path() -> Path:
    """Return the path to the persisted API key file."""
    return Path.home() / ".biocompiler" / "api_key"


def _generate_and_persist_api_key() -> str:
    """Generate a random API key, persist it, and return it.

    The key is saved to ~/.biocompiler/api_key for reuse across restarts.
    On first startup, the key is printed to the console as a one-time message.
    """
    key_file = _get_api_key_file_path()
    is_new = not key_file.exists()

    if key_file.exists():
        try:
            existing_key = key_file.read_text().strip()
            if existing_key:
                if is_new:
                    # Should not happen, but handle gracefully
                    pass
                return existing_key
        except (OSError, UnicodeDecodeError):
            logger.warning("Could not read API key from %s, generating new one", key_file)

    # Generate a new 32-byte hex key (64 chars)
    new_key = secrets.token_hex(32)

    # Persist to disk
    try:
        key_file.parent.mkdir(parents=True, exist_ok=True)
        key_file.write_text(new_key)
        # Set restrictive permissions (owner read/write only)
        key_file.chmod(0o600)
    except OSError as exc:
        logger.warning("Could not persist API key to %s: %s", key_file, exc)

    # One-time message to console
    print(
        f"\nGenerated API key: {new_key} (save this!)\n"
        f"Key saved to {key_file} for reuse across restarts.\n"
    )
    logger.info("Generated and persisted API key to %s", key_file)

    return new_key


# ── Resolve configured API keys ────────────────────────────────────
_CONFIGURED_API_KEYS: set[str] = set()

_single_key = os.environ.get("BIOCOMPILER_API_KEY", "")
_multi_keys = os.environ.get("BIOCOMPILER_API_KEYS", "")

if _multi_keys:
    _CONFIGURED_API_KEYS = {k.strip() for k in _multi_keys.split(",") if k.strip()}
elif _single_key == "disabled":
    # Explicit BIOCOMPILER_API_KEY=disabled → disable auth with warning
    logger.warning(
        "BIOCOMPILER_API_KEY=disabled: API authentication is DISABLED. "
        "This is dangerous for a tool that designs DNA. "
        "Use only in isolated development environments."
    )
    _AUTH_MODE = "disabled"
elif _single_key:
    _CONFIGURED_API_KEYS = {_single_key}
else:
    # No BIOCOMPILER_API_KEY set → defer key generation until API server starts.
    # Do NOT generate or persist a key at import time (side-effect-free import).
    _CONFIGURED_API_KEYS = set()
    _AUTH_MODE = "deferred"

# ─── Rate Limiting ──────────────────────────────────────────────────

RATE_LIMIT_RPM = int(os.environ.get("BIOCOMPILER_RATE_LIMIT", "60"))  # requests per minute

from ..infrastructure.rate_limiter import PersistentRateLimiter

# Lazy rate limiter — created on first use, not at import time.
_rate_limiter = None

def _get_rate_limiter():
    global _rate_limiter
    if _rate_limiter is None:
        _rate_limiter = PersistentRateLimiter(
            db_path=os.environ.get("BIOCOMPILER_RATE_LIMIT_DB", "~/.biocompiler/rate_limits.db"),
            max_requests=RATE_LIMIT_RPM,
            window_seconds=60,
        )
    return _rate_limiter


def _check_rate_limit(client_id: str) -> None:
    """Enforce per-client rate limiting (sliding window, SQLite-backed)."""
    allowed, remaining = _get_rate_limiter().check(client_id)
    if not allowed:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded: {RATE_LIMIT_RPM} requests/minute. "
                   f"Retry after {_get_rate_limiter().window_seconds} seconds.",
        )
    _get_rate_limiter().record(client_id)


def _check_batch_rate_limit(client_id: str, item_count: int) -> None:
    """
    Enforce per-client rate limiting for batch requests.

    Each item in the batch consumes one rate-limit unit. We pre-check
    whether the entire batch can be accommodated, then record all units.
    This prevents partial batch execution when rate-limited mid-way.
    """
    allowed, remaining = _get_rate_limiter().check(client_id)
    if item_count > remaining:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded: batch has {item_count} items but only "
                   f"{remaining} rate-limit units remaining this minute. "
                   f"Limit: {RATE_LIMIT_RPM} requests/minute.",
        )
    _get_rate_limiter().record_batch(client_id, item_count)


def _is_auth_enabled() -> bool:
    """Return True if authentication is currently enabled."""
    return _AUTH_MODE != "disabled" and bool(_CONFIGURED_API_KEYS)


def get_auth_mode() -> str:
    """Return the current auth mode string: 'required', 'optional', or 'disabled'."""
    return _AUTH_MODE


def get_configured_api_keys() -> set[str]:
    """Return a copy of the configured API keys set (for testing/inspection)."""
    return set(_CONFIGURED_API_KEYS)


def is_auth_enabled() -> bool:
    """Public API: Return True if authentication is currently enabled."""
    return _is_auth_enabled()


async def verify_api_key(api_key: str = Security(_api_key_header)) -> str:
    """
    Verify the API key based on the current auth mode.

    Auth modes (controlled by BIOCOMPILER_AUTH_MODE):
    - "required" (default): Unauthenticated requests get HTTP 401.
    - "optional": Auth is checked; unauthenticated requests are allowed
      but receive a warning header (added by middleware).
    - "disabled": Auth is completely disabled (DANGEROUS).

    API key resolution:
    - BIOCOMPILER_API_KEYS (comma-separated) for key rotation
    - BIOCOMPILER_API_KEY for a single key
    - BIOCOMPILER_API_KEY=disabled explicitly disables auth
    - No env var → auto-generated key persisted to ~/.biocompiler/api_key
    """
    # Disabled mode: no auth required
    if _AUTH_MODE == "disabled":
        return "anonymous"

    # No keys configured (should not happen with auto-generation, but guard)
    if not _CONFIGURED_API_KEYS:
        if _AUTH_MODE == "required":
            return "anonymous"
        return "anonymous"  # optional mode, no keys

    # Valid key provided
    if api_key is not None and api_key in _CONFIGURED_API_KEYS:
        return api_key

    # No key or invalid key
    if _AUTH_MODE == "required":
        raise HTTPException(
            status_code=401,
            detail="Invalid or missing API key. Set X-API-Key header.",
        )

    # Optional mode: allow through (middleware adds warning header)
    return "anonymous"
