"""Precomputed MHC binding databases for common HLA and mouse alleles.

This subpackage provides lazily-loaded, precomputed MHC binding databases
for common HLA and mouse MHC alleles. Each database contains binders and
non-binders generated from PSSM-based scoring, plus known epitopes
from IEDB and SYFPEITHI databases.

Supported alleles
-----------------
MHC class I (human):
- HLA-A*02:01 — Anchor: position 2 (L/M/V/I), position 9 (V/L/I/A)
- HLA-A*01:01 — Anchor: position 2 (T/S), position 9 (Y)
- HLA-A*03:01 — Anchor: position 2 (V/L/I), position 9 (K/R)
- HLA-B*07:02 — Anchor: position 2 (P), position 9 (L)
- HLA-B*08:01 — Anchor: position 3 (K/R), position 9 (L/I)

MHC class II (human):
- HLA-DRB1*01:01 — Anchor: P1 (F/Y/W), P4 (L/I/V), P6 (S/T), P9 (V/L)
- HLA-DRB1*04:01 — Anchor: P1 (F/Y/W), P4 (D/E), P6 (A/S/G), P9 (L/I/V)

MHC class I (mouse):
- H-2Kb — Anchor: position 5 (F/Y), position 8 (V/L/I)
- H-2Db — Anchor: position 5 (N/M), position 9 (M/V/I/L)

Usage
-----
>>> from biocompiler.mhc_binding_db.precomputed import get_precomputed_database
>>> db = get_precomputed_database("HLA-A*01:01")
>>> db.is_binder("YLDVSSNYI")
True
>>> db.binder_count
100

>>> from biocompiler.mhc_binding_db.precomputed import get_all_precomputed_databases
>>> all_dbs = get_all_precomputed_databases()
>>> list(all_dbs.keys())
['HLA-A*02:01', 'HLA-A*01:01', 'HLA-A*03:01', 'HLA-B*07:02',
 'HLA-B*08:01', 'HLA-DRB1*01:01', 'HLA-DRB1*04:01', 'H-2Kb', 'H-2Db']
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .. import PrecomputedAlleleDatabase

__all__ = [
    "get_all_precomputed_databases",
    "get_precomputed_database",
    "AVAILABLE_ALLELES",
]

#: Alleles for which precomputed databases are available.
AVAILABLE_ALLELES: list[str] = [
    "HLA-A*02:01",
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-DRB1*01:01",
    "HLA-DRB1*04:01",
    "H-2Kb",
    "H-2Db",
]

# Lazy cache — populated on first access
_databases: dict[str, PrecomputedAlleleDatabase] = {}


def get_precomputed_database(allele: str) -> PrecomputedAlleleDatabase:
    """Return the precomputed binding database for a single allele.

    Parameters
    ----------
    allele : str
        HLA allele name (e.g. ``"HLA-A*01:01"``). Must be one of
        :data:`AVAILABLE_ALLELES`.

    Returns
    -------
    PrecomputedAlleleDatabase
        The precomputed binding database for the requested allele.

    Raises
    ------
    ValueError
        If *allele* is not among the available precomputed alleles.

    Examples
    --------
    >>> db = get_precomputed_database("HLA-A*01:01")
    >>> db.allele
    'HLA-A*01:01'
    >>> db.total_entries
    300
    """
    if allele in _databases:
        return _databases[allele]

    if allele == "HLA-A*02:01":
        from .hla_a0201 import get_database
        db = get_database()
    elif allele == "HLA-A*01:01":
        from .hla_a0101 import get_database
        db = get_database()
    elif allele == "HLA-A*03:01":
        from .hla_a0301 import get_database
        db = get_database()
    elif allele == "HLA-B*07:02":
        from .hla_b0702 import get_database
        db = get_database()
    elif allele == "HLA-B*08:01":
        from .hla_b0801 import get_database
        db = get_database()
    elif allele == "HLA-DRB1*01:01":
        from .hla_drb1_0101 import get_database
        db = get_database()
    elif allele == "HLA-DRB1*04:01":
        from .hla_drb1_0401 import get_database
        db = get_database()
    elif allele == "H-2Kb":
        from .mouse_h2kb import get_database
        db = get_database()
    elif allele == "H-2Db":
        from .mouse_h2db import get_database
        db = get_database()
    else:
        raise ValueError(
            f"No precomputed database for allele {allele!r}. "
            f"Available alleles: {AVAILABLE_ALLELES}"
        )

    _databases[allele] = db
    return db


def get_all_precomputed_databases() -> dict[str, PrecomputedAlleleDatabase]:
    """Return all precomputed binding databases.

    Returns
    -------
    dict[str, PrecomputedAlleleDatabase]
        Mapping from allele name to the corresponding
        :class:`~biocompiler.mhc_binding_db.PrecomputedAlleleDatabase`.

    Examples
    --------
    >>> all_dbs = get_all_precomputed_databases()
    >>> len(all_dbs)
    9
    >>> all_dbs["HLA-B*07:02"].binder_count
    100
    """
    result: dict[str, PrecomputedAlleleDatabase] = {}
    for allele in AVAILABLE_ALLELES:
        result[allele] = get_precomputed_database(allele)
    return result
