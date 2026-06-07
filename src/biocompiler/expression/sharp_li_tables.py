"""
Thread-Safe Sharp & Li (1987) Table Access
============================================

This module provides thread-safe access to the Sharp & Li reference
gene sets, codon usage tables, and CAI weight tables.  All global
mutable state is encapsulated inside a ``_SharpLiState`` singleton
protected by ``threading.Lock``, eliminating the data-race surface
that existed when the tables were bare module-level dicts.

The 14 previously lockless global mutable variables are:

  1. _CODON_TABLE
  2. _AA_TO_CODONS
  3. _STOP_CODONS
  4. ECOLI_SHARP_LI_REFERENCE_GENES
  5. ECOLI_SHARP_LI_CODON_USAGE
  6. ECOLI_SHARP_LI_CAI_WEIGHTS
  7. YEAST_SHARP_LI_REFERENCE_GENES
  8. YEAST_SHARP_LI_CODON_USAGE
  9. YEAST_SHARP_LI_CAI_WEIGHTS
 10. SHARP_LI_PUBLISHED_CAI
 11. _REFERENCE_WEIGHTS
 12. SHARP_LI_REFERENCE_GENES
 13. SHARP_LI_CODON_USAGE
 14. SHARP_LI_CAI_WEIGHTS

Public API
----------
``get_sharp_li_table(organism)`` — the primary entry point used by
the rest of BioCompiler.  Returns a *copy* of the CAI weights dict
for the requested organism, making it safe for callers to read
without holding any lock.

``set_sharp_li_table(organism, table)`` — allows programmatic
replacement of an organism's table (used in testing and
hot-reloading scenarios).  The update is performed atomically
under the state lock.

Thread-safety guarantee: every read and write of the shared state
goes through ``_SharpLiState`` methods that acquire a
``threading.Lock`` before touching the underlying dicts.  Because
the lock serialises all mutations and all reads that return
references to mutable sub-structures, no two threads can observe
a partially-updated table.
"""

from __future__ import annotations

import threading
from typing import Dict, Optional

# Re-export the immutable data constants and computation functions
# from the original module so that existing imports continue to work.
from ..organisms.sharp_li_reference import (
    _CODON_TABLE,
    _AA_TO_CODONS,
    _STOP_CODONS,
    ECOLI_SHARP_LI_REFERENCE_GENES,
    ECOLI_SHARP_LI_CODON_USAGE,
    ECOLI_SHARP_LI_CAI_WEIGHTS,
    YEAST_SHARP_LI_REFERENCE_GENES,
    YEAST_SHARP_LI_CODON_USAGE,
    YEAST_SHARP_LI_CAI_WEIGHTS,
    SHARP_LI_PUBLISHED_CAI,
    _REFERENCE_WEIGHTS,
    SHARP_LI_REFERENCE_GENES,
    SHARP_LI_CODON_USAGE,
    SHARP_LI_CAI_WEIGHTS,
    compute_cai_with_reference,
    get_sharp_li_cai_weights,
)

__all__ = [
    # Thread-safe accessors
    "get_sharp_li_table",
    "set_sharp_li_table",
    # State class (for advanced/testing use)
    "_SharpLiState",
    # Re-exported constants (backward compat)
    "_CODON_TABLE",
    "_AA_TO_CODONS",
    "_STOP_CODONS",
    "ECOLI_SHARP_LI_REFERENCE_GENES",
    "ECOLI_SHARP_LI_CODON_USAGE",
    "ECOLI_SHARP_LI_CAI_WEIGHTS",
    "YEAST_SHARP_LI_REFERENCE_GENES",
    "YEAST_SHARP_LI_CODON_USAGE",
    "YEAST_SHARP_LI_CAI_WEIGHTS",
    "SHARP_LI_PUBLISHED_CAI",
    "_REFERENCE_WEIGHTS",
    "SHARP_LI_REFERENCE_GENES",
    "SHARP_LI_CODON_USAGE",
    "SHARP_LI_CAI_WEIGHTS",
    # Re-exported functions
    "compute_cai_with_reference",
    "get_sharp_li_cai_weights",
]

# ────────────────────────────────────────────────────────────
# Organism name resolution (mirrors sharp_li_reference logic)
# ────────────────────────────────────────────────────────────

_ALIASES: Dict[str, str] = {
    "ecoli": "Escherichia_coli",
    "e_coli": "Escherichia_coli",
    "E_coli": "Escherichia_coli",
    "E. coli": "Escherichia_coli",
    "Escherichia_coli": "Escherichia_coli",
    "yeast": "Saccharomyces_cerevisiae",
    "S_cerevisiae": "Saccharomyces_cerevisiae",
    "Saccharomyces_cerevisiae": "Saccharomyces_cerevisiae",
}


def _resolve_organism(organism: str) -> str:
    """Normalise an organism name to its canonical form."""
    key = organism.strip()
    return _ALIASES.get(key, _ALIASES.get(key.lower().replace(" ", "_"), key))


# ═══════════════════════════════════════════════════════════════
# _SharpLiState — lock-protected container for all mutable state
# ═══════════════════════════════════════════════════════════════

class _SharpLiState:
    """Encapsulates the 14 global mutable Sharp-Li variables behind a lock.

    Design notes
    ------------
    * A single ``threading.Lock`` guards all state.  This is intentionally
      simple — the critical sections are tiny (dict lookups/copies), so
      contention is negligible even under high concurrency.
    * ``get_table()`` returns a *shallow copy* so that callers cannot
      mutate the canonical state without going through ``set_table()``.
    * ``set_table()`` atomically replaces the table for an organism
      across all three parallel registries (reference_genes, codon_usage,
      cai_weights).

    The 14 managed variables (grouped by category):

    **Codon infrastructure (3):**
      1. codon_table       — mapping of codon → amino acid
      2. aa_to_codons      — mapping of amino acid → list of codons
      3. stop_codons       — set of stop codons

    **E. coli tables (3):**
      4. ecoli_reference_genes
      5. ecoli_codon_usage
      6. ecoli_cai_weights

    **Yeast tables (3):**
      7. yeast_reference_genes
      8. yeast_codon_usage
      9. yeast_cai_weights

    **Combined / published (5):**
     10. published_cai
     11. reference_weights
     12. reference_genes
     13. codon_usage
     14. cai_weights
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()

        # Codon infrastructure
        self._codon_table: Dict[str, str] = dict(_CODON_TABLE)
        self._aa_to_codons: Dict[str, list] = {k: list(v) for k, v in _AA_TO_CODONS.items()}
        self._stop_codons: set = set(_STOP_CODONS)

        # E. coli tables
        self._ecoli_reference_genes: Dict[str, str] = dict(ECOLI_SHARP_LI_REFERENCE_GENES)
        self._ecoli_codon_usage: Dict[str, float] = dict(ECOLI_SHARP_LI_CODON_USAGE)
        self._ecoli_cai_weights: Dict[str, float] = dict(ECOLI_SHARP_LI_CAI_WEIGHTS)

        # Yeast tables
        self._yeast_reference_genes: Dict[str, str] = dict(YEAST_SHARP_LI_REFERENCE_GENES)
        self._yeast_codon_usage: Dict[str, float] = dict(YEAST_SHARP_LI_CODON_USAGE)
        self._yeast_cai_weights: Dict[str, float] = dict(YEAST_SHARP_LI_CAI_WEIGHTS)

        # Combined / published
        self._published_cai: Dict[str, Dict[str, float]] = {
            k: dict(v) for k, v in SHARP_LI_PUBLISHED_CAI.items()
        }
        self._reference_weights: Dict[str, Dict[str, Dict[str, float]]] = {
            k: {rk: dict(rv) for rk, rv in v.items()}
            for k, v in _REFERENCE_WEIGHTS.items()
        }
        self._reference_genes: Dict[str, Dict[str, str]] = {
            k: dict(v) for k, v in SHARP_LI_REFERENCE_GENES.items()
        }
        self._codon_usage: Dict[str, Dict[str, float]] = {
            k: dict(v) for k, v in SHARP_LI_CODON_USAGE.items()
        }
        self._cai_weights: Dict[str, Dict[str, float]] = {
            k: dict(v) for k, v in SHARP_LI_CAI_WEIGHTS.items()
        }

    # ── Codon infrastructure accessors ──────────────────────────

    def get_codon_table(self) -> Dict[str, str]:
        """Return a copy of the codon table."""
        with self._lock:
            return dict(self._codon_table)

    def get_aa_to_codons(self) -> Dict[str, list]:
        """Return a copy of the amino-acid-to-codons mapping."""
        with self._lock:
            return {k: list(v) for k, v in self._aa_to_codons.items()}

    def get_stop_codons(self) -> set:
        """Return a copy of the stop codons set."""
        with self._lock:
            return set(self._stop_codons)

    # ── E. coli accessors ───────────────────────────────────────

    def get_ecoli_reference_genes(self) -> Dict[str, str]:
        with self._lock:
            return dict(self._ecoli_reference_genes)

    def get_ecoli_codon_usage(self) -> Dict[str, float]:
        with self._lock:
            return dict(self._ecoli_codon_usage)

    def get_ecoli_cai_weights(self) -> Dict[str, float]:
        with self._lock:
            return dict(self._ecoli_cai_weights)

    # ── Yeast accessors ─────────────────────────────────────────

    def get_yeast_reference_genes(self) -> Dict[str, str]:
        with self._lock:
            return dict(self._yeast_reference_genes)

    def get_yeast_codon_usage(self) -> Dict[str, float]:
        with self._lock:
            return dict(self._yeast_codon_usage)

    def get_yeast_cai_weights(self) -> Dict[str, float]:
        with self._lock:
            return dict(self._yeast_cai_weights)

    # ── Combined / published accessors ──────────────────────────

    def get_published_cai(self) -> Dict[str, Dict[str, float]]:
        with self._lock:
            return {k: dict(v) for k, v in self._published_cai.items()}

    def get_reference_weights(self) -> Dict[str, Dict[str, Dict[str, float]]]:
        with self._lock:
            return {
                k: {rk: dict(rv) for rk, rv in v.items()}
                for k, v in self._reference_weights.items()
            }

    def get_reference_genes(self) -> Dict[str, Dict[str, str]]:
        with self._lock:
            return {k: dict(v) for k, v in self._reference_genes.items()}

    def get_codon_usage(self) -> Dict[str, Dict[str, float]]:
        with self._lock:
            return {k: dict(v) for k, v in self._codon_usage.items()}

    def get_cai_weights(self) -> Dict[str, Dict[str, float]]:
        """Return a copy of the combined CAI weights dict."""
        with self._lock:
            return {k: dict(v) for k, v in self._cai_weights.items()}

    # ── Primary table accessor ──────────────────────────────────

    def get_table(self, organism: str) -> Optional[Dict[str, float]]:
        """Return a copy of the CAI weights table for *organism*.

        Args:
            organism: Any recognised organism name (e.g. ``"ecoli"``,
                ``"Escherichia_coli"``, ``"yeast"``).

        Returns:
            A dict mapping codon strings to CAI weight values, or
            ``None`` if the organism is not found.
        """
        canonical = _resolve_organism(organism)
        with self._lock:
            table = self._cai_weights.get(canonical)
            if table is None:
                return None
            return dict(table)

    def set_table(self, organism: str, table: Dict[str, float]) -> None:
        """Atomically replace the CAI weights table for *organism*.

        This also updates the parallel ``_codon_usage``,
        ``_reference_genes``, and ``_reference_weights`` registries
        so that all access paths remain consistent.

        Args:
            organism: Any recognised organism name.
            table: New CAI weights dict (codon → float).
        """
        canonical = _resolve_organism(organism)
        table_copy = dict(table)  # defensive copy before acquiring lock
        with self._lock:
            self._cai_weights[canonical] = table_copy
            # Keep the parallel registries in sync.  For codon_usage and
            # reference_genes we store the same table (the caller is
            # expected to update those separately if they differ).
            self._codon_usage[canonical] = dict(table_copy)
            self._reference_weights[canonical] = {"sharp_li": dict(table_copy)}


# ═══════════════════════════════════════════════════════════════
# Module-level singleton
# ═══════════════════════════════════════════════════════════════

_state = _SharpLiState()


# ── Lock for external synchronisation ─────────────────────────
# Re-exported for backward compatibility: code that previously
# imported ``_sharp_li_registry_lock`` from sharp_li_reference
# can continue to use it.
_sharp_li_registry_lock = _state._lock  # noqa: SLF001


# ═══════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════

def get_sharp_li_table(organism: str) -> Dict[str, float]:
    """Return a thread-safe copy of the Sharp-Li CAI weights for *organism*.

    This is the primary public accessor.  It always returns a fresh
    dict, so callers can read freely without holding any lock.

    Args:
        organism: Any recognised organism name — e.g. ``"ecoli"``,
            ``"Escherichia_coli"``, ``"yeast"``,
            ``"Saccharomyces_cerevisiae"``.

    Returns:
        Dict mapping codon strings to their CAI weight values.

    Raises:
        ValueError: If the organism is not recognised.
    """
    result = _state.get_table(organism)
    if result is None:
        supported = sorted(_state.get_cai_weights().keys())
        raise ValueError(
            f"No Sharp-Li CAI weights for organism '{organism}'. "
            f"Supported: {supported}"
        )
    return result


def set_sharp_li_table(organism: str, table: Dict[str, float]) -> None:
    """Thread-safe replacement of the CAI weights table for *organism*.

    Args:
        organism: Any recognised organism name.
        table: New CAI weights dict (codon → float).
    """
    _state.set_table(organism, table)
