"""Deprecated: use biocompiler.shared.thread_safety instead."""
import warnings

warnings.warn(
    "biocompiler._thread_safety is deprecated — use biocompiler.shared.thread_safety instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.thread_safety import *  # noqa: F401,F403

__all__ = [
    "ThreadSafeDict",
    "ThreadSafeDefaultDict",
    "ThreadSafeLazy",
]
