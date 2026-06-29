"""
BioCompiler Infrastructure Layer
=================================

Infrastructure concerns that should be separated from domain logic:

  - rate_limiter.py     : SQLite-backed persistent rate limiting
  - biopython_compat.py : BioPython SeqRecord interop + subprocess wrappers
  - dna_chisel_compat.py: DNA Chisel integration + comparative benchmarking
  - file_io.py          : Common file I/O operations (path resolution, file reading)
  - lims.py             : LIMS platform integration (Benchling, LabGuru)

This package follows DDD separation of concerns: domain logic should NOT
depend on infrastructure directly. Instead, domain modules define interfaces
(protocols/ABCs) and infrastructure provides concrete implementations.

All modules here are also accessible via backward-compat shims at their
original locations (e.g. ``biocompiler.shared.rate_limiter``).
"""

from biocompiler.infrastructure.rate_limiter import PersistentRateLimiter  # noqa: F401
from .file_io import resolve_input, looks_like_path  # noqa: F401

__all__ = [
    "PersistentRateLimiter",
    "resolve_input",
    "looks_like_path",
]
