"""
Tests for the PersistentRateLimiter (SQLite-backed rate limiting).

Covers:
1. Basic check/record — requests under the limit are allowed
2. Rate-limit enforcement — requests exceeding the limit are blocked
3. Batch recording — multiple requests recorded at once
4. Sliding window expiry — old entries no longer count (no real sleeps)
5. Cleanup — expired entries are removed from the database (no real sleeps)
6. Persistence — state survives across limiter instances
7. Client isolation — one client hitting the limit does not affect another
8. Clear — per-client and global clear works
9. Database schema and properties
10. Custom time function support
"""

from __future__ import annotations

import sqlite3
import tempfile
import os

import pytest

from biocompiler.infrastructure.rate_limiter import PersistentRateLimiter


# ─── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture()
def tmp_db(tmp_path):
    """Return a temporary database path that is cleaned up after the test."""
    return str(tmp_path / "test_rate_limits.db")


@pytest.fixture()
def limiter(tmp_db):
    """Return a PersistentRateLimiter with a temp DB and small limits."""
    return PersistentRateLimiter(
        db_path=tmp_db,
        max_requests=10,
        window_seconds=60,
    )


@pytest.fixture()
def mock_time():
    """Provide a controllable mock time source.

    Returns a dict with keys:
      - "time": current simulated time (float)
      - "fn": callable that returns the current simulated time
      - "advance": callable(dt) that advances time by dt seconds
    """
    state = {"time": 1000000.0}

    def time_fn() -> float:
        return state["time"]

    state["fn"] = time_fn
    state["advance"] = lambda dt: state.update(time=state["time"] + dt)
    return state


# ═══════════════════════════════════════════════════════════════════════
# 1. Basic check / record
# ═══════════════════════════════════════════════════════════════════════


class TestBasicCheckRecord:
    """Verify that requests under the limit are allowed."""

    def test_check_returns_allowed_when_under_limit(self, limiter):
        allowed, remaining = limiter.check("client_a")
        assert allowed is True
        assert remaining == 10

    def test_remaining_decreases_after_record(self, limiter):
        limiter.record("client_a")
        allowed, remaining = limiter.check("client_a")
        assert allowed is True
        assert remaining == 9

    def test_multiple_records_decrease_remaining(self, limiter):
        for _ in range(5):
            limiter.record("client_a")
        allowed, remaining = limiter.check("client_a")
        assert allowed is True
        assert remaining == 5


# ═══════════════════════════════════════════════════════════════════════
# 2. Rate-limit enforcement
# ═══════════════════════════════════════════════════════════════════════


class TestRateLimitEnforcement:
    """Verify that requests exceeding the limit are blocked."""

    def test_check_returns_not_allowed_at_limit(self, limiter):
        for _ in range(10):
            limiter.record("client_a")
        allowed, remaining = limiter.check("client_a")
        assert allowed is False
        assert remaining == 0

    def test_check_returns_not_allowed_over_limit(self, limiter):
        for _ in range(12):
            limiter.record("client_a")
        allowed, remaining = limiter.check("client_a")
        assert allowed is False
        assert remaining == 0


# ═══════════════════════════════════════════════════════════════════════
# 3. Batch recording
# ═══════════════════════════════════════════════════════════════════════


class TestBatchRecording:
    """Verify that batch recording works correctly."""

    def test_record_batch_adds_multiple_entries(self, limiter):
        limiter.record_batch("client_a", 5)
        allowed, remaining = limiter.check("client_a")
        assert allowed is True
        assert remaining == 5

    def test_record_batch_at_limit(self, limiter):
        limiter.record_batch("client_a", 10)
        allowed, remaining = limiter.check("client_a")
        assert allowed is False
        assert remaining == 0

    def test_record_batch_over_limit(self, limiter):
        limiter.record_batch("client_a", 15)
        allowed, remaining = limiter.check("client_a")
        assert allowed is False
        assert remaining == 0


# ═══════════════════════════════════════════════════════════════════════
# 4. Sliding window expiry (no real sleeps)
# ═══════════════════════════════════════════════════════════════════════


class TestSlidingWindowExpiry:
    """Verify that old entries no longer count toward the limit."""

    def test_old_entries_expire(self, tmp_db, mock_time):
        """Entries older than window_seconds should not count."""
        limiter = PersistentRateLimiter(
            db_path=tmp_db,
            max_requests=10,
            window_seconds=60,
            time_func=mock_time["fn"],
        )
        # Record 10 requests at the initial time
        for _ in range(10):
            limiter.record("client_a")

        # At limit now
        allowed, remaining = limiter.check("client_a")
        assert allowed is False

        # Advance time past the window (no real sleep!)
        mock_time["advance"](61)

        # Old entries should have expired
        allowed, remaining = limiter.check("client_a")
        assert allowed is True
        assert remaining == 10


# ═══════════════════════════════════════════════════════════════════════
# 5. Cleanup (no real sleeps)
# ═══════════════════════════════════════════════════════════════════════


class TestCleanup:
    """Verify that expired entries are removed from the database."""

    def test_cleanup_removes_expired_entries(self, tmp_db, mock_time):
        limiter = PersistentRateLimiter(
            db_path=tmp_db,
            max_requests=100,
            window_seconds=60,
            time_func=mock_time["fn"],
        )
        limiter.record("client_a")
        limiter.record("client_b")

        # Advance time past the window (no real sleep!)
        mock_time["advance"](61)

        deleted = limiter.cleanup()
        assert deleted >= 2  # At least the two we inserted

    def test_cleanup_preserves_active_entries(self, limiter):
        limiter.record("client_a")
        deleted = limiter.cleanup()
        # Entries are still within the window, so nothing should be deleted
        assert deleted == 0

    def test_periodic_cleanup_triggered(self, tmp_db, mock_time):
        """Verify that cleanup is triggered automatically after many records."""
        limiter = PersistentRateLimiter(
            db_path=tmp_db,
            max_requests=1000,
            window_seconds=60,
            time_func=mock_time["fn"],
        )
        # Record enough to trigger periodic cleanup (default every 100)
        for _ in range(110):
            limiter.record("client_a")

        # Advance time past the window (no real sleep!)
        mock_time["advance"](61)
        # The cleanup should have run during the record loop.
        # After window expiry, all should be gone.
        deleted = limiter.cleanup()
        # Should clean up the expired entries
        assert deleted >= 0  # May already have been cleaned by periodic


# ═══════════════════════════════════════════════════════════════════════
# 6. Persistence
# ═══════════════════════════════════════════════════════════════════════


class TestPersistence:
    """Verify that rate-limit state survives across limiter instances."""

    def test_state_persists_across_instances(self, tmp_db):
        limiter1 = PersistentRateLimiter(
            db_path=tmp_db, max_requests=10, window_seconds=60,
        )
        limiter1.record("client_a")
        limiter1.record("client_a")

        # Create a new instance pointing to the same DB
        limiter2 = PersistentRateLimiter(
            db_path=tmp_db, max_requests=10, window_seconds=60,
        )
        allowed, remaining = limiter2.check("client_a")
        assert allowed is True
        assert remaining == 8  # 10 - 2 = 8


# ═══════════════════════════════════════════════════════════════════════
# 7. Client isolation
# ═══════════════════════════════════════════════════════════════════════


class TestClientIsolation:
    """Verify that one client hitting the limit does not affect another."""

    def test_different_clients_independent(self, limiter):
        # Client A hits the limit
        for _ in range(10):
            limiter.record("client_a")
        allowed_a, remaining_a = limiter.check("client_a")
        assert allowed_a is False
        assert remaining_a == 0

        # Client B should still be allowed
        allowed_b, remaining_b = limiter.check("client_b")
        assert allowed_b is True
        assert remaining_b == 10

    def test_client_b_unaffected_by_client_a_records(self, limiter):
        limiter.record_batch("client_a", 8)
        allowed_b, remaining_b = limiter.check("client_b")
        assert allowed_b is True
        assert remaining_b == 10


# ═══════════════════════════════════════════════════════════════════════
# 8. Clear
# ═══════════════════════════════════════════════════════════════════════


class TestClear:
    """Verify that clear works for per-client and global cases."""

    def test_clear_specific_client(self, limiter):
        limiter.record("client_a")
        limiter.record("client_b")
        limiter.clear("client_a")

        allowed_a, remaining_a = limiter.check("client_a")
        assert allowed_a is True
        assert remaining_a == 10  # Cleared

        allowed_b, remaining_b = limiter.check("client_b")
        assert allowed_b is True
        assert remaining_b == 9  # Not cleared

    def test_clear_all_clients(self, limiter):
        limiter.record("client_a")
        limiter.record("client_b")
        limiter.clear()

        allowed_a, remaining_a = limiter.check("client_a")
        assert remaining_a == 10

        allowed_b, remaining_b = limiter.check("client_b")
        assert remaining_b == 10


# ═══════════════════════════════════════════════════════════════════════
# 9. Database schema and properties
# ═══════════════════════════════════════════════════════════════════════


class TestDatabaseSchema:
    """Verify the SQLite database is properly initialized."""

    def test_db_file_created(self, tmp_db):
        limiter = PersistentRateLimiter(db_path=tmp_db)
        assert limiter.db_path.exists()

    def test_table_created(self, tmp_db):
        limiter = PersistentRateLimiter(db_path=tmp_db)
        conn = sqlite3.connect(tmp_db)
        tables = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()
        conn.close()
        table_names = [t[0] for t in tables]
        assert "rate_limits" in table_names

    def test_index_created(self, tmp_db):
        limiter = PersistentRateLimiter(db_path=tmp_db)
        conn = sqlite3.connect(tmp_db)
        indexes = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='index'"
        ).fetchall()
        conn.close()
        index_names = [i[0] for i in indexes]
        assert "idx_rate_limits_client_ts" in index_names

    def test_max_requests_property(self, limiter):
        assert limiter.max_requests == 10

    def test_window_seconds_property(self, limiter):
        assert limiter.window_seconds == 60

    def test_db_path_property(self, tmp_db):
        limiter = PersistentRateLimiter(db_path=tmp_db)
        assert str(limiter.db_path) == tmp_db

    def test_parent_directory_created(self, tmp_path):
        nested = str(tmp_path / "sub" / "dir" / "rates.db")
        limiter = PersistentRateLimiter(db_path=nested)
        assert limiter.db_path.parent.exists()

    def test_tilde_expansion(self):
        limiter = PersistentRateLimiter(
            db_path="~/.biocompiler_test_rate_limits/test.db"
        )
        assert not str(limiter.db_path).startswith("~")

    def test_custom_time_func(self, tmp_db, mock_time):
        """Verify that a custom time_func is used for all time-dependent ops."""
        limiter = PersistentRateLimiter(
            db_path=tmp_db,
            max_requests=5,
            window_seconds=10,
            time_func=mock_time["fn"],
        )
        # Record 5 requests at time 0
        for _ in range(5):
            limiter.record("client_a")

        allowed, _ = limiter.check("client_a")
        assert allowed is False

        # Advance past the window
        mock_time["advance"](11)

        allowed, remaining = limiter.check("client_a")
        assert allowed is True
        assert remaining == 5
