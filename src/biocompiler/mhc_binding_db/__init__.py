"""Offline MHC binding prediction database for BioCompiler immunogenicity analysis.

This package provides precomputed MHC binding predictions for common HLA
alleles, generated using PSSM-based scoring from the immunogenicity module.
The data includes known epitopes from the IEDB and SYFPEITHI databases, plus
PSSM-derived predictions for ~100 binders and ~200 non-binders per allele.

The database is designed for offline use — no external API calls or network
access is required at runtime.  This makes it suitable for reproducible
pipelines and air-gapped environments.  Each allele's database contains
both experimentally verified epitopes (sourced from IEDB/SYFPEITHI) and
computationally predicted binders, enabling quick lookup during gene
optimization to assess immunogenicity risk.

Binding predictions are classified into three tiers: strong binder, weak
binder, and non-binder.  The classification is based on estimated IC50
values: strong binders have IC50 < 50 nM, weak binders 50–500 nM, and
non-binders ≥ 500 nM.  These thresholds follow the IEDB classification
convention.

Two database classes are provided:

- :class:`MHCBindingDatabase` — multi-allele in-memory database with
  ``(allele, peptide)`` key lookup, JSON serialisation, and batch queries.
  This is the primary database class for new code.
- :class:`PrecomputedAlleleDatabase` — legacy single-allele database
  used by the precomputed subpackage for backward compatibility.

Usage (new multi-allele database)::

    from biocompiler.mhc_binding_db import MHCBindingDatabase, MHCBindingRecord

    # Create and populate a database
    db = MHCBindingDatabase()
    db.add(MHCBindingRecord(
        allele="HLA-A*02:01", peptide="LLFGYPVYV",
        ic50_nm=12.0, rank=0.5, binding_class="strong_binder",
        source="mhcflurry_predicted", method="mhcflurry_class1",
        timestamp="2025-01-01T00:00:00Z",
    ))

    # Lookup
    result = db.lookup("HLA-A*02:01", "LLFGYPVYV")

    # Filter
    db.strong_binders("HLA-A*02:01")
    db.weak_binders("HLA-A*02:01")

    # Stats
    db.stats()

Usage (legacy single-allele database)::

    from biocompiler.mhc_binding_db import PrecomputedAlleleDatabase, get_database

    # Load the database for a specific HLA allele
    db = get_database("HLA-A*01:01")

    # Check if a peptide is a known/predicted binder
    if db.is_binder("YLDVSSNYI"):
        print("Binder detected!")

Accuracy
--------
All binding scores and IC50 values are derived from PSSM-based predictions
with expected AUC-ROC of 0.60–0.75. These are heuristic approximations and
should not replace experimental validation or NetMHCpan predictions where
available.

References:
  Nielsen, M. & Andreatta, M. (2016). "NetMHCpan-3.0; improved prediction
  of binding to MHC class I molecules." *Journal of Immunological Research*
  2016:641649. doi:10.1155/2016/6416495
  Rammensee, H.-G. et al. (1999). "SYFPEITHI: database for MHC ligands
  and peptide motifs." *Immunogenetics* 50:213–219.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

# Import new schema types — these are the primary public API
from .schema import (
    MHCBindingDatabase,
    MHCBindingRecord,
    generate_fallback_database,
    get_default_alleles_for_organism,
    is_mhc_class2,
)

__all__ = [
    # New schema (primary API)
    "MHCBindingRecord",
    "MHCBindingDatabase",
    "generate_fallback_database",
    "get_default_alleles_for_organism",
    "is_mhc_class2",
    # Legacy precomputed types
    "PrecomputedEntry",
    "PrecomputedAlleleDatabase",
    "get_database",
]


@dataclass(frozen=True)
class PrecomputedEntry:
    """A single precomputed MHC binding prediction entry (legacy).

    Attributes
    ----------
    peptide : str
        Amino acid sequence of the peptide (9-mer for MHC-I,
        13-25aa for MHC-II).
    binding_score : float
        Normalised binding score in [0, 1] (1 = strongest binder).
    ic50_nm : float
        Estimated IC50 in nM.
    binding_class : str
        One of ``"strong_binder"``, ``"moderate_binder"``,
        ``"weak_binder"``, ``"non_binder"``.
    anchor_residues : dict[int, str]
        Position index (0-based) -> amino acid at anchor positions.
    anchor_scores : dict[int, float]
        Position index (0-based) -> PSSM score contribution at anchor.
    source : str
        Data source: ``"known_epitope"`` for experimentally verified
        binders from IEDB/SYFPEITHI, or ``"pssm_predicted"`` for
        computationally derived predictions.
    peptide_length : int
        Length of the peptide in residues.  For MHC-I this is the
        core binding length (8-9 aa); for MHC-II this records the
        variable-length peptide (typically 13-25 aa).  A value of 0
        indicates that the field was not set (backward compat with
        older MHC-I entries).
    """

    peptide: str
    binding_score: float
    ic50_nm: float
    binding_class: str
    anchor_residues: dict[int, str]
    anchor_scores: dict[int, float]
    source: str = "pssm_predicted"
    peptide_length: int = 0

    @property
    def is_binder(self) -> bool:
        """True if classified as strong or moderate binder."""
        return self.binding_class in ("strong_binder", "moderate_binder")


@dataclass
class PrecomputedAlleleDatabase:
    """Precomputed MHC binding database for a single HLA allele (legacy).

    This is the legacy single-allele database class used by the
    :mod:`precomputed` subpackage.  For new code, prefer the
    multi-allele :class:`MHCBindingDatabase` from :mod:`schema`.

    Parameters
    ----------
    allele : str
        HLA allele name (e.g. ``"HLA-A*01:01"``).
    peptide_length : int
        Peptide length for this allele (default 9 for MHC-I).
    anchor_description : str
        Human-readable description of anchor residues.
    known_epitopes : list[str]
        Experimentally verified binding epitopes from IEDB/SYFPEITHI.
    entries : list[PrecomputedEntry]
        All precomputed binding entries (binders + non-binders).
    """

    allele: str
    peptide_length: int = 9
    anchor_description: str = ""
    known_epitopes: list[str] = field(default_factory=list)
    entries: list[PrecomputedEntry] = field(default_factory=list)

    # ── Lookup caches (built lazily) ──────────────────────────────────
    _peptide_index: dict[str, PrecomputedEntry] = field(
        default_factory=dict, repr=False,
    )

    def __post_init__(self) -> None:
        """Build lookup index from entries."""
        for entry in self.entries:
            self._peptide_index[entry.peptide] = entry

    # ── Public API ────────────────────────────────────────────────────

    def lookup(self, peptide: str) -> Optional[PrecomputedEntry]:
        """Look up a peptide in the precomputed database.

        Parameters
        ----------
        peptide : str
            Amino acid sequence of the peptide.

        Returns
        -------
        PrecomputedEntry or None
            The precomputed entry if found, else ``None``.
        """
        return self._peptide_index.get(peptide)

    def is_binder(self, peptide: str) -> bool:
        """Check if a peptide is a known/predicted binder.

        Returns ``True`` if the peptide is classified as a strong or
        moderate binder. Returns ``False`` for weak binders, non-binders,
        and peptides not in the database.
        """
        entry = self._peptide_index.get(peptide)
        if entry is None:
            return False
        return entry.is_binder

    def get_binders(self) -> list[PrecomputedEntry]:
        """Return all entries classified as binders (strong + moderate)."""
        return [e for e in self.entries if e.is_binder]

    def get_non_binders(self) -> list[PrecomputedEntry]:
        """Return all entries classified as non-binders (weak + non)."""
        return [e for e in self.entries if not e.is_binder]

    def get_strong_binders(self) -> list[PrecomputedEntry]:
        """Return all entries classified as strong binders."""
        return [e for e in self.entries if e.binding_class == "strong_binder"]

    def get_moderate_binders(self) -> list[PrecomputedEntry]:
        """Return all entries classified as moderate binders."""
        return [e for e in self.entries if e.binding_class == "moderate_binder"]

    def get_weak_binders(self) -> list[PrecomputedEntry]:
        """Return all entries classified as weak binders."""
        return [e for e in self.entries if e.binding_class == "weak_binder"]

    def get_known_epitope_entries(self) -> list[PrecomputedEntry]:
        """Return entries sourced from known experimental epitopes."""
        return [e for e in self.entries if e.source == "known_epitope"]

    @property
    def binder_count(self) -> int:
        """Number of strong + moderate binders."""
        return sum(1 for e in self.entries if e.is_binder)

    @property
    def non_binder_count(self) -> int:
        """Number of weak binders + non-binders."""
        return sum(1 for e in self.entries if not e.is_binder)

    @property
    def total_entries(self) -> int:
        """Total number of precomputed entries."""
        return len(self.entries)

    def __repr__(self) -> str:
        return (
            f"PrecomputedAlleleDatabase("
            f"allele={self.allele!r}, "
            f"entries={self.total_entries}, "
            f"binders={self.binder_count}, "
            f"non_binders={self.non_binder_count})"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Lazy database loading
# ═══════════════════════════════════════════════════════════════════════════

_database_cache: dict[str, PrecomputedAlleleDatabase] = {}


def get_database(allele: str) -> Optional[PrecomputedAlleleDatabase]:
    """Lazy-load a pre-computed binding database for the given allele.

    This function loads and caches the database on first access for each
    allele, so subsequent calls for the same allele return the cached
    instance.

    Parameters
    ----------
    allele : str
        HLA allele name (e.g. ``"HLA-A*02:01"``).

    Returns
    -------
    PrecomputedAlleleDatabase or None
        The pre-computed binding database if data exists for the allele,
        else ``None``.
    """
    if allele in _database_cache:
        return _database_cache[allele]

    try:
        from .precomputed import get_precomputed_database, AVAILABLE_ALLELES
        if allele in AVAILABLE_ALLELES:
            db = get_precomputed_database(allele)
            _database_cache[allele] = db
            return db
        return None
    except Exception as exc:
        import logging
        logging.getLogger(__name__).debug(
            "Failed to load pre-computed database for %s: %s", allele, exc
        )
        return None
