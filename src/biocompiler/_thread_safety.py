"""
Thread Safety Utilities for BioCompiler

Provides simple lock-based wrappers for global mutable state that must
be safe for concurrent access (e.g., when the API serves multiple
requests simultaneously via threading).

Design principle: simplicity over performance.  We use threading.Lock
(not RLock, not rwlock) because correctness is the priority and the
critical sections are tiny (dict lookups, list appends).

Usage patterns:

    1. ThreadSafeDict — wraps a dict with a lock for read/write access
    2. ThreadSafeDefaultDict — wraps a defaultdict with a lock
    3. ThreadSafeLazy — wraps a lazily-computed value with a lock

All wrappers preserve the public API of the underlying data structure
so that existing code continues to work without changes.
"""

from __future__ import annotations

import threading
from collections import defaultdict
from typing import Any, Callable, Generic, TypeVar

__all__ = [
    "ThreadSafeDict",
    "ThreadSafeDefaultDict",
    "ThreadSafeLazy",
]

T = TypeVar("T")
K = TypeVar("K")
V = TypeVar("V")


class ThreadSafeDict(Generic[K, V]):
    """A dict wrapper that serialises all access with a lock.

    Instead of inheriting from dict (which leaks unguarded methods),
    this class stores the underlying dict as a private attribute and
    exposes only the operations we actually use in BioCompiler.

    For backward compatibility, ``__getitem__``, ``__setitem__``,
    ``__contains__``, and ``get()`` behave identically to a plain dict.
    """

    def __init__(self, initial: dict[K, V] | None = None) -> None:
        self._data: dict[K, V] = dict(initial) if initial else {}
        self._lock = threading.Lock()

    # ── Read operations ────────────────────────────────────────────

    def __getitem__(self, key: K) -> V:
        with self._lock:
            return self._data[key]

    def get(self, key: K, default: V | None = None) -> V | None:
        with self._lock:
            return self._data.get(key, default)

    def __contains__(self, key: Any) -> bool:
        with self._lock:
            return key in self._data

    def __len__(self) -> int:
        with self._lock:
            return len(self._data)

    def keys(self):  # type: ignore[override]
        with self._lock:
            return list(self._data.keys())

    def values(self):  # type: ignore[override]
        with self._lock:
            return list(self._data.values())

    def items(self):  # type: ignore[override]
        with self._lock:
            return list(self._data.items())

    # ── Write operations ───────────────────────────────────────────

    def __setitem__(self, key: K, value: V) -> None:
        with self._lock:
            self._data[key] = value

    def __delitem__(self, key: K) -> None:
        with self._lock:
            del self._data[key]

    def pop(self, key: K, *args: Any) -> V | None:
        with self._lock:
            return self._data.pop(key, *args)

    def update(self, other: dict[K, V] | None = None, **kwargs: Any) -> None:
        with self._lock:
            if other:
                self._data.update(other)
            if kwargs:
                self._data.update(kwargs)

    def clear(self) -> None:
        with self._lock:
            self._data.clear()

    def setdefault(self, key: K, default: V | None = None) -> V | None:
        with self._lock:
            return self._data.setdefault(key, default)

    # ── Atomic compound operations ─────────────────────────────────
    # These perform a read-modify-write under a single lock acquisition,
    # which is the main reason this class exists — plain dict access
    # from multiple threads can lose updates.

    def get_and_set(self, key: K, value: V) -> V | None:
        """Atomically get the old value and set a new one."""
        with self._lock:
            old = self._data.get(key)
            self._data[key] = value
            return old

    def update_if_absent(self, key: K, value: V) -> V:
        """Set *key* = *value* only if *key* is absent; return the final value."""
        with self._lock:
            if key in self._data:
                return self._data[key]
            self._data[key] = value
            return value

    # ── Snapshot ───────────────────────────────────────────────────

    def snapshot(self) -> dict[K, V]:
        """Return a shallow copy of the current state (thread-safe)."""
        with self._lock:
            return dict(self._data)

    # ── Representation ─────────────────────────────────────────────

    def __repr__(self) -> str:
        with self._lock:
            return f"ThreadSafeDict({self._data!r})"

    def __bool__(self) -> bool:
        with self._lock:
            return bool(self._data)


class ThreadSafeDefaultDict(Generic[K, V]):
    """A defaultdict wrapper that serialises all access with a lock.

    This is needed because ``defaultdict`` creates entries on ``__getitem__``
    misses, which is a read-modify-write that must be atomic.
    """

    def __init__(self, default_factory: Callable[[], V]) -> None:
        self._data: defaultdict[K, V] = defaultdict(default_factory)
        self._lock = threading.Lock()

    def __getitem__(self, key: K) -> V:
        with self._lock:
            return self._data[key]

    def __setitem__(self, key: K, value: V) -> None:
        with self._lock:
            self._data[key] = value

    def __contains__(self, key: Any) -> bool:
        with self._lock:
            return key in self._data

    def __len__(self) -> int:
        with self._lock:
            return len(self._data)

    def __delitem__(self, key: K) -> None:
        with self._lock:
            del self._data[key]

    def get(self, key: K, default: V | None = None) -> V | None:
        """Get without creating a default entry."""
        with self._lock:
            return self._data.get(key, default)

    def keys(self):  # type: ignore[override]
        with self._lock:
            return list(self._data.keys())

    def values(self):  # type: ignore[override]
        with self._lock:
            return list(self._data.values())

    def items(self):  # type: ignore[override]
        with self._lock:
            return list(self._data.items())

    def clear(self) -> None:
        with self._lock:
            self._data.clear()

    def snapshot(self) -> dict[K, V]:
        """Return a shallow copy of the current state (thread-safe)."""
        with self._lock:
            return dict(self._data)

    def __repr__(self) -> str:
        with self._lock:
            return f"ThreadSafeDefaultDict({dict(self._data)!r})"

    def __bool__(self) -> bool:
        with self._lock:
            return bool(self._data)


class ThreadSafeLazy(Generic[T]):
    """A lazily-computed value protected by a lock.

    Many BioCompiler modules initialise a global dict to ``{}`` and
    populate it on first access (e.g. ``SHARP_LI_ADAPTIVENESS_TABLES``).
    If two threads trigger the lazy init simultaneously, they could
    both compute and assign, wasting work or, worse, partially
    overwrite each other's results.

    Usage::

        _TABLES = ThreadSafeLazy(lambda: _compute_tables())

        def get_tables():
            return _TABLES.get()

    The computation function is called at most once, even if multiple
    threads call ``get()`` concurrently.
    """

    def __init__(self, compute_fn: Callable[[], T]) -> None:
        self._compute_fn = compute_fn
        self._value: T | None = None
        self._computed = False
        self._lock = threading.Lock()

    def get(self) -> T:
        """Return the value, computing it on first call (thread-safe)."""
        if self._computed:
            # Fast path: already computed, no lock needed for read
            # (Python's GIL makes reading _computed atomic for bool)
            return self._value  # type: ignore[return-value]

        with self._lock:
            if self._computed:
                return self._value  # type: ignore[return-value]
            self._value = self._compute_fn()
            self._computed = True
            return self._value  # type: ignore[return-value]

    def reset(self) -> None:
        """Force recomputation on next ``get()`` (mainly for testing)."""
        with self._lock:
            self._value = None
            self._computed = False

    @property
    def is_computed(self) -> bool:
        return self._computed

    def __repr__(self) -> str:
        state = "computed" if self._computed else "pending"
        return f"ThreadSafeLazy({state})"
