"""
Persistent SQLite-backed rate limiter for BioCompiler API.

Unlike the previous in-memory ``defaultdict(list)`` approach, this
implementation stores request timestamps in SQLite so that:

* Rate-limit state survives process restarts.
* Multiple worker processes share the same limit (no per-worker gaps).
* Old entries are cleaned up automatically.

The schema uses a single table with ``(client_id, timestamp)`` rows.
A sliding-window algorithm counts rows within the last *window_seconds*
for a given *client_id*.
"""

from __future__ import annotations

import sqlite3
import threading
import time
import logging
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

__all__ = ["PersistentRateLimiter"]


class PersistentRateLimiter:
    """SQLite-backed rate limiter that persists across restarts.

    Parameters
    ----------
    db_path:
        Path to the SQLite database file.  The ``~`` prefix is expanded
        and parent directories are created automatically.
    max_requests:
        Maximum number of requests allowed within the sliding window.
    window_seconds:
        Width of the sliding window in seconds.
    time_func:
        Callable returning the current time as a float (seconds since epoch).
        Defaults to ``time.time``.  Override in tests to control time.
    """

    def __init__(
        self,
        db_path: str = "~/.biocompiler/rate_limits.db",
        max_requests: int = 100,
        window_seconds: int = 3600,
        time_func: Any | None = None,
    ) -> None:
        self._db_path = Path(db_path).expanduser()
        self._db_path.parent.mkdir(parents=True, exist_ok=True)
        self._max_requests = max_requests
        self._window_seconds = window_seconds
        self._time_func = time_func or time.time
        self._request_counter = 0  # for periodic cleanup
        self._cleanup_every = 100  # clean up after this many record() calls
        self._lock = threading.Lock()
        self._init_db()

    # ── Database initialisation ────────────────────────────────────

    def _init_db(self) -> None:
        """Create the rate limit table if it does not exist."""
        with self._connect() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS rate_limits (
                    client_id  TEXT    NOT NULL,
                    timestamp  REAL    NOT NULL
                )
                """
            )
            # Index for fast lookups by client within the window
            conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_rate_limits_client_ts
                ON rate_limits (client_id, timestamp)
                """
            )
            conn.commit()

    def _connect(self) -> sqlite3.Connection:
        """Return a new SQLite connection with WAL journal mode."""
        conn = sqlite3.connect(str(self._db_path), timeout=10)
        conn.execute("PRAGMA journal_mode=WAL")
        return conn

    # ── Public API ─────────────────────────────────────────────────

    def check(self, client_id: str) -> tuple[bool, int]:
        """Check if *client_id* is within the rate limit.

        Returns
        -------
        (allowed, remaining)
            *allowed* is ``True`` when the client has not exceeded
            *max_requests* within the current window.  *remaining* is
            the number of requests the client can still make.
        """
        now = self._time_func()
        cutoff = now - self._window_seconds
        with self._connect() as conn:
            count = conn.execute(
                """
                SELECT COUNT(*) FROM rate_limits
                WHERE client_id = ? AND timestamp > ?
                """,
                (client_id, cutoff),
            ).fetchone()[0]
        remaining = max(0, self._max_requests - count)
        allowed = count < self._max_requests
        return allowed, remaining

    def record(self, client_id: str) -> None:
        """Record a request from *client_id* at the current time.

        Also triggers periodic cleanup of expired entries.
        """
        now = self._time_func()
        with self._connect() as conn:
            conn.execute(
                """
                INSERT INTO rate_limits (client_id, timestamp)
                VALUES (?, ?)
                """,
                (client_id, now),
            )
            conn.commit()

        # Periodic cleanup
        with self._lock:
            self._request_counter += 1
            if self._request_counter >= self._cleanup_every:
                self._request_counter = 0
                self.cleanup()

    def record_batch(self, client_id: str, count: int) -> None:
        """Record *count* requests from *client_id* at the current time.

        Used for batch endpoints where each item consumes one rate-limit
        unit.
        """
        now = self._time_func()
        with self._connect() as conn:
            conn.executemany(
                """
                INSERT INTO rate_limits (client_id, timestamp)
                VALUES (?, ?)
                """,
                [(client_id, now)] * count,
            )
            conn.commit()

        # Periodic cleanup
        with self._lock:
            self._request_counter += count
            if self._request_counter >= self._cleanup_every:
                self._request_counter = 0
                self.cleanup()

    def cleanup(self) -> int:
        """Remove expired entries older than the window.

        Returns the number of rows deleted.
        """
        cutoff = self._time_func() - self._window_seconds
        with self._connect() as conn:
            cursor = conn.execute(
                """
                DELETE FROM rate_limits WHERE timestamp <= ?
                """,
                (cutoff,),
            )
            conn.commit()
            deleted = cursor.rowcount
        if deleted:
            logger.debug("Rate-limiter cleanup: removed %d expired entries", deleted)
        return deleted

    # ── Utility / testing helpers ──────────────────────────────────

    def clear(self, client_id: Optional[str] = None) -> None:
        """Clear rate-limit records.

        If *client_id* is given, only that client's records are removed.
        Otherwise all records are deleted.  Useful for testing.
        """
        with self._connect() as conn:
            if client_id is not None:
                conn.execute(
                    "DELETE FROM rate_limits WHERE client_id = ?",
                    (client_id,),
                )
            else:
                conn.execute("DELETE FROM rate_limits")
            conn.commit()

    @property
    def db_path(self) -> Path:
        """Return the resolved database path."""
        return self._db_path

    @property
    def max_requests(self) -> int:
        return self._max_requests

    @property
    def window_seconds(self) -> int:
        return self._window_seconds
