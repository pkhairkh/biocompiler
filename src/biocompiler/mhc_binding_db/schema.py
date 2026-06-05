"""MHC binding database schema and loader.

This module defines the core data types and database class for MHC binding
prediction storage, lookup, and offline fallback generation.

The :class:`MHCBindingRecord` dataclass represents a single binding prediction
or measurement.  The :class:`MHCBindingDatabase` class provides an in-memory
store keyed by ``(allele, peptide)`` for fast lookups across multiple alleles,
with JSON serialisation, batch queries, and statistical summaries.

For environments without MHCflurry or NetMHCpan, :func:`generate_fallback_database`
produces a PSSM-based approximation that requires no external dependencies.

Usage
-----
>>> from biocompiler.mhc_binding_db.schema import MHCBindingDatabase, MHCBindingRecord
>>> db = MHCBindingDatabase()
>>> db.add(MHCBindingRecord(
...     allele="HLA-A*02:01",
...     peptide="LLFGYPVYV",
...     ic50_nm=12.0,
...     rank=0.5,
...     binding_class="strong_binder",
...     source="mhcflurry_predicted",
...     method="mhcflurry_class1_presentation",
...     timestamp="2025-01-01T00:00:00Z",
... ))
>>> result = db.lookup("HLA-A*02:01", "LLFGYPVYV")
>>> result.ic50_nm
12.0
"""
from __future__ import annotations

import json
import math
import random
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from itertools import product
from typing import Any, Optional

__all__ = [
    "MHCBindingRecord",
    "MHCBindingDatabase",
    "generate_fallback_database",
    "get_default_alleles_for_organism",
    "is_mhc_class2",
]

# ═══════════════════════════════════════════════════════════════════════════
# Valid values for constrained string fields
# ═══════════════════════════════════════════════════════════════════════════

_VALID_BINDING_CLASSES: frozenset[str] = frozenset({
    "strong_binder",
    "weak_binder",
    "non_binder",
})

_VALID_SOURCES: frozenset[str] = frozenset({
    "mhcflurry_predicted",
    "iedb_verified",
    "netmhcpan_predicted",
    "pssm_fallback",
})

# ═══════════════════════════════════════════════════════════════════════════
# Standard amino acids used for fallback generation
# ═══════════════════════════════════════════════════════════════════════════

STANDARD_AMINO_ACIDS: str = "ACDEFGHIKLMNPQRSTVWY"

# ═══════════════════════════════════════════════════════════════════════════
# PSSM fallback data — simplified position-specific scoring matrices
# for common HLA alleles (9-mer core).  These are intentionally coarse
# approximations derived from published binding motif literature.
# ═══════════════════════════════════════════════════════════════════════════

_PSSM_FALLBACK: dict[str, list[dict[str, float]]] = {
    "HLA-A*02:01": [
        {"L": 1.2, "M": 1.2, "I": 1.2, "V": 1.2, "A": 1.1, "F": 1.1},  # pos 0 (defaults 1.0)
        {"L": 2.0, "M": 2.0, "I": 1.8, "V": 1.8},  # pos 1 — anchor
        {},  # pos 2
        {},  # pos 3
        {},  # pos 4
        {},  # pos 5
        {},  # pos 6
        {},  # pos 7
        {"V": 1.5, "L": 1.5, "I": 1.3, "A": 1.2},  # pos 8 — anchor
    ],
    "HLA-A*01:01": [
        {"A": 1.1, "S": 1.1},
        {"T": 1.8, "S": 1.6, "D": 1.5, "E": 1.5},  # pos 1 — anchor
        {},
        {},
        {},
        {},
        {},
        {},
        {"Y": 1.8, "F": 1.6},  # pos 8 — anchor
    ],
    "HLA-A*03:01": [
        {"A": 1.1, "S": 1.1},
        {"V": 1.8, "I": 1.8, "L": 1.6, "M": 1.6},  # pos 1 — anchor
        {},
        {},
        {},
        {},
        {},
        {},
        {"K": 2.0, "R": 1.8, "H": 1.4},  # pos 8 — anchor
    ],
    "HLA-A*24:02": [
        {"Y": 1.2, "F": 1.1},
        {"Y": 2.0, "F": 2.0, "W": 1.8},  # pos 1 — anchor
        {},
        {},
        {},
        {},
        {},
        {},
        {"F": 1.5, "L": 1.5, "I": 1.3},  # pos 8 — anchor
    ],
    "HLA-B*07:02": [
        {"A": 1.1, "P": 1.1},
        {"P": 2.0, "A": 1.8},  # pos 1 — anchor
        {},
        {},
        {},
        {},
        {},
        {},
        {"L": 1.5, "I": 1.5, "V": 1.3},  # pos 8 — anchor
    ],
    "HLA-B*08:01": [
        {},
        {"K": 1.8, "R": 1.8},  # pos 1 — anchor
        {},
        {},
        {},
        {},
        {},
        {},
        {"L": 1.5, "I": 1.3, "V": 1.3},  # pos 8 — anchor
    ],
    "HLA-DRB1*01:01": [
        {"F": 1.8, "Y": 1.7, "W": 1.6, "L": 1.5, "I": 1.4, "V": 1.4, "M": 1.3},
        {},
        {},
        {"A": 1.6, "S": 1.4, "T": 1.4, "N": 1.3, "G": 1.2},
        {},
        {"L": 1.5, "I": 1.4, "V": 1.4, "M": 1.3, "F": 1.3},
        {},
        {},
        {"K": 1.3, "R": 1.3, "N": 1.2, "Q": 1.2, "E": 1.1, "D": 1.1},
    ],
    "HLA-DRB1*04:01": [
        {"F": 1.8, "Y": 1.7, "W": 1.6, "L": 1.3},
        {},
        {},
        {"D": 1.8, "E": 1.6},
        {},
        {"A": 1.6, "S": 1.4, "G": 1.3, "N": 1.2},
        {},
        {},
        {"L": 1.4, "I": 1.3, "V": 1.3, "F": 1.3},
    ],
}

# Disfavored residues at anchor positions for fallback scoring
_PSSM_DISFAVORED: dict[str, list[dict[str, float]]] = {
    "HLA-A*02:01": [
        {"D": 0.5, "E": 0.5, "K": 0.5, "R": 0.5},
        {"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "P": 0.4},
        {}, {}, {}, {}, {}, {},
        {"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4, "P": 0.4},
    ],
    "HLA-A*01:01": [
        {"W": 0.5, "R": 0.5},
        {"L": 0.4, "I": 0.4, "V": 0.5, "F": 0.4},
        {}, {}, {}, {}, {}, {},
        {"K": 0.4, "R": 0.4, "D": 0.5, "E": 0.5},
    ],
    "HLA-A*03:01": [
        {},
        {"D": 0.4, "E": 0.4, "N": 0.5, "Q": 0.5},
        {}, {}, {}, {}, {}, {},
        {"D": 0.3, "E": 0.3, "S": 0.5, "T": 0.5},
    ],
    "HLA-A*24:02": [
        {},
        {"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "P": 0.4},
        {}, {}, {}, {}, {}, {},
        {"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
    ],
    "HLA-B*07:02": [
        {},
        {"D": 0.3, "E": 0.3, "K": 0.3, "R": 0.3, "W": 0.4},
        {}, {}, {}, {}, {}, {},
        {"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
    ],
    "HLA-B*08:01": [
        {},
        {"D": 0.3, "E": 0.3, "P": 0.4, "G": 0.5},
        {}, {}, {}, {}, {}, {},
        {"D": 0.4, "E": 0.4, "K": 0.5},
    ],
    "HLA-DRB1*01:01": [
        {"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
        {},
        {},
        {"W": 0.5, "F": 0.6, "Y": 0.6},
        {},
        {"D": 0.5, "E": 0.5, "K": 0.5},
        {},
        {},
        {},
    ],
    "HLA-DRB1*04:01": [
        {"D": 0.4, "E": 0.4, "K": 0.5, "R": 0.5},
        {},
        {},
        {"K": 0.4, "R": 0.4, "W": 0.5},
        {},
        {"W": 0.5, "F": 0.6, "Y": 0.6},
        {},
        {},
        {"D": 0.5, "E": 0.5, "K": 0.5},
    ],
}

# IC50 mapping constants (log-linear, consistent with immunogenicity module)
_IC50_LOG_INTERCEPT: float = 3.949
_IC50_LOG_SLOPE: float = 2.5
_PSSM_CONTRAST_POWER: float = 2.0
_PSSM_UNKNOWN_AA_SCORE: float = 0.3

# ═══════════════════════════════════════════════════════════════════════════
# MHC class classification helper
# ═══════════════════════════════════════════════════════════════════════════

# MHC-II allele prefixes — these alleles bind longer peptides (13–25 aa)
_MHC_II_PREFIXES: tuple[str, ...] = (
    "HLA-DRB1",
    "HLA-DRB3",
    "HLA-DRB4",
    "HLA-DRB5",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DPA1",
    "HLA-DPB1",
    "H2-IAb",
    "H2-IAe",
    "H2-IEb",
)


def is_mhc_class2(allele: str) -> bool:
    """Return True if *allele* is an MHC class II allele.

    MHC-II molecules (HLA-DR, -DQ, -DP in humans; I-A, I-E in mice)
    bind longer peptides (13–25 residues) compared to MHC-I (8–11).

    Parameters
    ----------
    allele : str
        MHC allele name (e.g. ``"HLA-DRB1*01:01"``).

    Returns
    -------
    bool
        ``True`` if the allele is MHC class II.
    """
    return any(allele.startswith(prefix) for prefix in _MHC_II_PREFIXES)


# Peptide length ranges by MHC class
_MHC_I_PEPTIDE_LENGTH_RANGE: tuple[int, int] = (8, 11)
_MHC_II_PEPTIDE_LENGTH_RANGE: tuple[int, int] = (13, 25)


# ═══════════════════════════════════════════════════════════════════════════
# Data classes
# ═══════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class MHCBindingRecord:
    """A single MHC binding prediction or measurement record.

    Attributes
    ----------
    allele : str
        HLA allele name (e.g. ``"HLA-A*02:01"``).
    peptide : str
        Amino acid sequence of the peptide.  For MHC-I alleles the
        length must be 8–11 residues; for MHC-II alleles 13–25 residues.
    ic50_nm : float
        Predicted or measured IC50 in nM.
    rank : float or None
        Percentile rank from the prediction tool, if available.
    binding_class : str
        One of ``"strong_binder"``, ``"weak_binder"``, ``"non_binder"``.
    source : str
        Data provenance: ``"mhcflurry_predicted"``, ``"iedb_verified"``,
        ``"netmhcpan_predicted"``, or ``"pssm_fallback"``.
    method : str
        Specific method or model version used.
    timestamp : str
        ISO-8601 timestamp when the record was created or loaded.
    peptide_length_range : tuple[int, int] or None
        Optional override for peptide length validation.  When ``None``
        (default), the range is inferred from the allele: MHC-I uses
        (8, 11), MHC-II uses (13, 25).  Set explicitly to bypass
        allele-based inference.
    """

    allele: str
    peptide: str
    ic50_nm: float
    rank: Optional[float]
    binding_class: str
    source: str
    method: str
    timestamp: str
    peptide_length_range: Optional[tuple[int, int]] = None

    def __post_init__(self) -> None:
        """Validate field values after construction."""
        # Determine allowed peptide length range
        if self.peptide_length_range is not None:
            min_len, max_len = self.peptide_length_range
        elif is_mhc_class2(self.allele):
            min_len, max_len = _MHC_II_PEPTIDE_LENGTH_RANGE
        else:
            min_len, max_len = _MHC_I_PEPTIDE_LENGTH_RANGE

        pep_len = len(self.peptide)
        if not (min_len <= pep_len <= max_len):
            raise ValueError(
                f"Peptide length must be {min_len}-{max_len} residues for "
                f"{'MHC-II' if is_mhc_class2(self.allele) else 'MHC-I'} "
                f"allele {self.allele!r}, got {pep_len} "
                f"for peptide {self.peptide!r}"
            )
        if self.binding_class not in _VALID_BINDING_CLASSES:
            raise ValueError(
                f"binding_class must be one of {sorted(_VALID_BINDING_CLASSES)}, "
                f"got {self.binding_class!r}"
            )
        if self.source not in _VALID_SOURCES:
            raise ValueError(
                f"source must be one of {sorted(_VALID_SOURCES)}, "
                f"got {self.source!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# Database class
# ═══════════════════════════════════════════════════════════════════════════


class MHCBindingDatabase:
    """In-memory MHC binding database with fast ``(allele, peptide)`` lookup.

    The database stores :class:`MHCBindingRecord` instances keyed by the
    tuple ``(allele, peptide)``.  This allows efficient queries across
    multiple alleles without needing a separate database per allele.

    Examples
    --------
    >>> db = MHCBindingDatabase()
    >>> db.add(MHCBindingRecord(
    ...     allele="HLA-A*02:01", peptide="LLFGYPVYV",
    ...     ic50_nm=12.0, rank=0.5, binding_class="strong_binder",
    ...     source="mhcflurry_predicted", method="mhcflurry_class1",
    ...     timestamp="2025-01-01T00:00:00Z",
    ... ))
    >>> db.lookup("HLA-A*02:01", "LLFGYPVYV")
    MHCBindingRecord(allele='HLA-A*02:01', peptide='LLFGYPVYV', ...)
    >>> db.strong_binders("HLA-A*02:01")
    ['LLFGYPVYV']
    """

    def __init__(self) -> None:
        self._records: dict[tuple[str, str], MHCBindingRecord] = {}

    # ── Mutation ────────────────────────────────────────────────────────

    def add(self, record: MHCBindingRecord) -> None:
        """Add or replace a binding record.

        If a record with the same ``(allele, peptide)`` key already exists
        it is silently replaced.
        """
        self._records[(record.allele, record.peptide)] = record

    def add_batch(self, records: list[MHCBindingRecord]) -> int:
        """Add multiple records.  Returns the number of records added."""
        for rec in records:
            self._records[(rec.allele, rec.peptide)] = rec
        return len(records)

    # ── Lookup ──────────────────────────────────────────────────────────

    def lookup(self, allele: str, peptide: str) -> Optional[MHCBindingRecord]:
        """Look up a single record by allele and peptide.

        Returns ``None`` if no matching record exists.
        """
        return self._records.get((allele, peptide))

    def lookup_batch(
        self,
        allele: str,
        peptides: list[str],
    ) -> list[Optional[MHCBindingRecord]]:
        """Look up multiple peptides for a given allele.

        Returns a list of the same length as *peptides* where each element
        is the :class:`MHCBindingRecord` or ``None`` if not found.
        """
        return [self._records.get((allele, p)) for p in peptides]

    # ── Filtering ───────────────────────────────────────────────────────

    def binding_peptides(
        self,
        allele: str,
        threshold_ic50: float = 500.0,
    ) -> list[str]:
        """Return all peptides for *allele* with IC50 below *threshold_ic50*.

        Parameters
        ----------
        allele : str
            HLA allele name.
        threshold_ic50 : float
            IC50 threshold in nM (default 500.0, the standard weak-binder
            cutoff).

        Returns
        -------
        list[str]
            Peptide sequences sorted alphabetically.
        """
        results: list[str] = []
        for (a, pep), rec in self._records.items():
            if a == allele and rec.ic50_nm < threshold_ic50:
                results.append(pep)
        return sorted(results)

    def strong_binders(self, allele: str) -> list[str]:
        """Return all strong binder peptides (IC50 < 50 nM) for *allele*.

        Returns
        -------
        list[str]
            Peptide sequences sorted alphabetically.
        """
        results: list[str] = []
        for (a, pep), rec in self._records.items():
            if a == allele and rec.ic50_nm < 50.0:
                results.append(pep)
        return sorted(results)

    def weak_binders(self, allele: str) -> list[str]:
        """Return all weak binder peptides (50 nM ≤ IC50 < 500 nM) for *allele*.

        Returns
        -------
        list[str]
            Peptide sequences sorted alphabetically.
        """
        results: list[str] = []
        for (a, pep), rec in self._records.items():
            if a == allele and 50.0 <= rec.ic50_nm < 500.0:
                results.append(pep)
        return sorted(results)

    # ── Statistics ──────────────────────────────────────────────────────

    def stats(self) -> dict[str, Any]:
        """Compute summary statistics across all records.

        Returns
        -------
        dict
            ``total_records``: total count.
            ``alleles``: dict mapping allele name to a sub-dict with keys
            ``total``, ``strong_binder``, ``weak_binder``, ``non_binder``.
        """
        allele_stats: dict[str, dict[str, int]] = {}
        for rec in self._records.values():
            a = rec.allele
            if a not in allele_stats:
                allele_stats[a] = {
                    "total": 0,
                    "strong_binder": 0,
                    "weak_binder": 0,
                    "non_binder": 0,
                }
            allele_stats[a]["total"] += 1
            if rec.binding_class in allele_stats[a]:
                allele_stats[a][rec.binding_class] += 1

        return {
            "total_records": len(self._records),
            "alleles": allele_stats,
        }

    # ── Serialisation ───────────────────────────────────────────────────

    def load_from_json(self, filepath: str) -> int:
        """Load records from a JSON file, merging into the current database.

        The JSON file should contain a list of objects, each matching the
        :class:`MHCBindingRecord` field names.

        Parameters
        ----------
        filepath : str
            Path to the JSON file.

        Returns
        -------
        int
            Number of records loaded.
        """
        with open(filepath, "r", encoding="utf-8") as fh:
            data = json.load(fh)

        records_data = data if isinstance(data, list) else data.get("records", [])
        count = 0
        for item in records_data:
            rec = MHCBindingRecord(
                allele=item["allele"],
                peptide=item["peptide"],
                ic50_nm=float(item["ic50_nm"]),
                rank=item.get("rank"),
                binding_class=item["binding_class"],
                source=item["source"],
                method=item["method"],
                timestamp=item["timestamp"],
            )
            self._records[(rec.allele, rec.peptide)] = rec
            count += 1
        return count

    def save_to_json(self, filepath: str) -> None:
        """Export all records to a JSON file.

        The output is a JSON object with a ``"records"`` key containing a
        list of record dicts.

        Parameters
        ----------
        filepath : str
            Destination file path.
        """
        records_list = [asdict(rec) for rec in self._records.values()]
        payload = {
            "version": "1.0",
            "exported_at": datetime.now(timezone.utc).isoformat(),
            "total_records": len(records_list),
            "records": records_list,
        }
        with open(filepath, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, default=str)

    # ── Container protocol ──────────────────────────────────────────────

    def __len__(self) -> int:
        return len(self._records)

    def __contains__(self, key: tuple[str, str]) -> bool:
        return key in self._records

    def __repr__(self) -> str:
        alleles = sorted({a for a, _ in self._records})
        return (
            f"MHCBindingDatabase("
            f"records={len(self._records)}, "
            f"alleles={alleles})"
        )

    @property
    def alleles(self) -> list[str]:
        """Return sorted list of alleles present in the database."""
        return sorted({a for a, _ in self._records})

    @property
    def records(self) -> list[MHCBindingRecord]:
        """Return all records as a list."""
        return list(self._records.values())


# ═══════════════════════════════════════════════════════════════════════════
# Fallback database generation
# ═══════════════════════════════════════════════════════════════════════════


def _score_peptide_pssm(
    peptide: str,
    preferred: list[dict[str, float]],
    disfavored: list[dict[str, float]],
    anchor_positions: Optional[set[int]] = None,
) -> float:
    """Score a peptide against a PSSM (preferred + disfavored rows).

    Returns a normalised binding score in [0, 1].

    Parameters
    ----------
    peptide : str
        Amino acid sequence.
    preferred : list[dict[str, float]]
        Per-position preferred residue scores.
    disfavored : list[dict[str, float]]
        Per-position disfavored residue scores.
    anchor_positions : set[int] or None
        0-indexed positions that are anchor residues.  These receive
        higher weight in the weighted geometric mean.  When ``None``,
        all positions are weighted equally (flat geometric mean).
    """
    n = len(peptide)
    scores: list[float] = []
    weights: list[float] = []

    # Anchor positions get double weight (2.0), non-anchors get 1.0
    anchor_weight = 2.0
    for i in range(n):
        aa = peptide[i]
        # Check preferred first
        if i < len(preferred) and aa in preferred[i]:
            scores.append(preferred[i][aa])
        elif i < len(disfavored) and aa in disfavored[i]:
            scores.append(disfavored[i][aa])
        else:
            scores.append(1.0)  # neutral

        if anchor_positions and i in anchor_positions:
            weights.append(anchor_weight)
        else:
            weights.append(1.0)

    # Weighted geometric mean
    total_weight = sum(weights)
    log_sum = sum(w * math.log(max(s, 1e-10)) for w, s in zip(weights, scores))
    geo_mean = math.exp(log_sum / total_weight)

    # Compute max and min possible weighted geometric means for normalisation
    max_scores: list[float] = []
    min_scores: list[float] = []
    for i in range(n):
        pos_vals = [1.0]  # default score
        if i < len(preferred):
            pos_vals.extend(preferred[i].values())
        if i < len(disfavored):
            pos_vals.extend(disfavored[i].values())
        max_scores.append(max(pos_vals))
        min_scores.append(min(pos_vals))

    max_log_sum = sum(w * math.log(max(s, 1e-10)) for w, s in zip(weights, max_scores))
    min_log_sum = sum(w * math.log(max(s, 1e-10)) for w, s in zip(weights, min_scores))
    max_geo = math.exp(max_log_sum / total_weight)
    min_geo = math.exp(min_log_sum / total_weight)

    if max_geo <= min_geo:
        return 0.0

    raw = (geo_mean - min_geo) / (max_geo - min_geo)
    raw = max(0.0, min(1.0, raw))
    normalised = raw ** _PSSM_CONTRAST_POWER
    return max(0.0, min(1.0, normalised))


def _binding_score_to_ic50(score: float) -> float:
    """Map a binding score to estimated IC50 (nM) via log-linear mapping."""
    clamped = max(0.0, min(1.0, score))
    return 10.0 ** (_IC50_LOG_INTERCEPT - _IC50_LOG_SLOPE * clamped)


def _classify_binding_fallback(ic50: float) -> str:
    """Classify binding by IC50 thresholds.

    Uses the 3-class system: strong_binder (<50), weak_binder (50-500),
    non_binder (>=500).
    """
    if ic50 < 50.0:
        return "strong_binder"
    elif ic50 < 500.0:
        return "weak_binder"
    else:
        return "non_binder"


# Anchor positions for PSSM scoring, keyed by allele.
# For MHC-I 9-mers: positions 1 and 8 (0-indexed) are the primary anchors.
# For MHC-II: positions 0, 3, 5, 8 (0-indexed within the 9-mer core) are
# the key pocket positions (P1, P4, P6, P9 in 1-based notation).
_PSSM_ANCHOR_POSITIONS: dict[str, set[int]] = {
    "HLA-A*02:01": {1, 8},
    "HLA-A*01:01": {1, 8},
    "HLA-A*03:01": {1, 8},
    "HLA-A*24:02": {1, 8},
    "HLA-B*07:02": {1, 8},
    "HLA-B*08:01": {1, 8},
    "HLA-DRB1*01:01": {0, 3, 5, 8},
    "HLA-DRB1*04:01": {0, 3, 5, 8},
}


def generate_fallback_database(
    alleles: list[str],
    peptide_lengths: Optional[list[int]] = None,
) -> MHCBindingDatabase:
    """Generate a fallback MHC binding database using PSSM scoring.

    For each allele with a known PSSM, this generates representative
    peptides by sampling amino acids with anchor-position bias,
    scores them, and stores the results.  For alleles without a known
    PSSM, a generic scoring model is applied.

    This provides OFFLINE capability when MHCflurry or NetMHCpan are
    unavailable, but is NOT as accurate as those tools (expected
    AUC-ROC 0.60–0.75).

    Parameters
    ----------
    alleles : list[str]
        HLA allele names to generate predictions for.
    peptide_lengths : list[int] or None
        Peptide lengths to consider.  When ``None`` (default), MHC-I
        alleles use ``[9]`` and MHC-II alleles use ``[15]`` (the most
        common MHC-II binding length).  PSSM scoring is applied for
        alleles with known PSSM data; other lengths use a generic
        heuristic.

    Returns
    -------
    MHCBindingDatabase
        Database populated with PSSM-scored records, all marked as
        ``source="pssm_fallback"``.
    """
    db = MHCBindingDatabase()
    now = datetime.now(timezone.utc).isoformat()
    rng = random.Random(42)  # deterministic seed for reproducibility

    for allele in alleles:
        preferred = _PSSM_FALLBACK.get(allele, [])
        disfavored = _PSSM_DISFAVORED.get(allele, [])
        anchor_positions = _PSSM_ANCHOR_POSITIONS.get(allele, set())
        is_class2 = is_mhc_class2(allele)

        # Determine peptide lengths for this allele
        if peptide_lengths is not None:
            lengths = peptide_lengths
        elif is_class2:
            lengths = [15]  # 15-mer is the most common MHC-II binding length
        else:
            lengths = [9]   # 9-mer is the standard MHC-I binding length

        # PSSM core length is always 9 for the scoring matrices
        pssm_core_len = len(preferred) if preferred else 0

        for pep_len in lengths:
            if pssm_core_len > 0:
                _generate_pssm_peptides(
                    db, allele, pep_len, pssm_core_len,
                    preferred, disfavored, anchor_positions,
                    rng, now, is_class2,
                )
            else:
                # No PSSM data for this allele — use generic
                _generate_generic_peptides(
                    db, allele, pep_len, rng, now,
                )

    return db


def _generate_pssm_peptides(
    db: MHCBindingDatabase,
    allele: str,
    pep_len: int,
    pssm_core_len: int,
    preferred: list[dict[str, float]],
    disfavored: list[dict[str, float]],
    anchor_positions: set[int],
    rng: random.Random,
    timestamp: str,
    is_class2: bool,
    num_samples: int = 200,
) -> None:
    """Generate peptides scored against a PSSM.

    For MHC-I (pep_len == pssm_core_len), the PSSM is applied directly.
    For MHC-II (pep_len > pssm_core_len), the longer peptide is scored
    by sliding the 9-mer PSSM core across all possible register positions
    and taking the best score.
    """
    seen_peptides: set[str] = set()

    # Determine which positions to bias during generation
    # For MHC-I: anchor positions within the peptide
    # For MHC-II: anchor positions within the 9-mer core
    anchor_aa_per_pos: dict[int, list[str]] = {}
    for pos in anchor_positions:
        if pos < len(preferred) and preferred[pos]:
            anchor_aa_per_pos[pos] = list(preferred[pos].keys())

    for _ in range(num_samples * 10):  # oversample to get enough unique
        if len(seen_peptides) >= num_samples:
            break

        if pep_len == pssm_core_len:
            # MHC-I: direct 9-mer generation
            residues: list[str] = []
            for pos in range(pep_len):
                if pos in anchor_aa_per_pos:
                    aa = rng.choice(anchor_aa_per_pos[pos] + list(STANDARD_AMINO_ACIDS))
                else:
                    aa = rng.choice(list(STANDARD_AMINO_ACIDS))
                residues.append(aa)
            peptide = "".join(residues)
        else:
            # MHC-II or non-standard length: generate full-length peptide,
            # biasing the core region (first 9 positions as representative core)
            residues = []
            for pos in range(pep_len):
                # Map to core position for bias (use first 9 positions)
                core_pos = pos if pos < pssm_core_len else -1
                if core_pos in anchor_aa_per_pos:
                    aa = rng.choice(anchor_aa_per_pos[core_pos] + list(STANDARD_AMINO_ACIDS))
                else:
                    aa = rng.choice(list(STANDARD_AMINO_ACIDS))
                residues.append(aa)
            peptide = "".join(residues)

        if peptide in seen_peptides:
            continue
        seen_peptides.add(peptide)

        # Score the peptide
        if pep_len == pssm_core_len:
            # Direct PSSM scoring for MHC-I 9-mers
            score = _score_peptide_pssm(
                peptide, preferred, disfavored, anchor_positions,
            )
        else:
            # Sliding-window PSSM scoring for MHC-II longer peptides:
            # slide the 9-mer core across all possible registers and take
            # the best (highest) binding score.
            best_score = 0.0
            for offset in range(pep_len - pssm_core_len + 1):
                core_peptide = peptide[offset:offset + pssm_core_len]
                s = _score_peptide_pssm(
                    core_peptide, preferred, disfavored, anchor_positions,
                )
                if s > best_score:
                    best_score = s
            score = best_score

        ic50 = _binding_score_to_ic50(score)
        binding_class = _classify_binding_fallback(ic50)

        db.add(MHCBindingRecord(
            allele=allele,
            peptide=peptide,
            ic50_nm=round(ic50, 1),
            rank=None,
            binding_class=binding_class,
            source="pssm_fallback",
            method="pssm_anchor_scoring",
            timestamp=timestamp,
        ))


def _generate_generic_peptides(
    db: MHCBindingDatabase,
    allele: str,
    pep_len: int,
    rng: random.Random,
    timestamp: str,
    num_samples: int = 100,
) -> None:
    """Generate peptides with a generic (non-PSSM) scoring model.

    Uses a simple hydrophobicity-based heuristic: peptides with
    hydrophobic residues at the termini score higher for MHC-I, while
    MHC-II uses hydrophobic/aromatic residues at pocket positions.
    """
    is_class2 = is_mhc_class2(allele)

    # Simple hydrophobicity scale (Kyte-Doolittle, normalised)
    _hydro: dict[str, float] = {
        "I": 1.0, "V": 0.93, "L": 0.89, "F": 0.88, "C": 0.61,
        "M": 0.56, "A": 0.42, "G": 0.0, "T": -0.17, "S": -0.27,
        "W": -0.33, "Y": -0.39, "P": -0.47, "H": -0.61, "E": -0.78,
        "Q": -0.72, "D": -0.83, "N": -0.64, "K": -1.0, "R": -1.0,
    }

    seen: set[str] = set()
    for _ in range(num_samples * 20):
        if len(seen) >= num_samples:
            break

        residues = [rng.choice(list(STANDARD_AMINO_ACIDS)) for _ in range(pep_len)]
        peptide = "".join(residues)
        if peptide in seen:
            continue
        seen.add(peptide)

        if is_class2:
            # MHC-II: score based on hydrophobicity at pocket positions
            # (P1, P4, P6, P9 in 1-based → 0, 3, 5, 8 in 0-based)
            pocket_positions = [0, 3, 5, min(8, pep_len - 1)]
            pocket_hydro = sum(
                _hydro.get(peptide[p], 0.0) for p in pocket_positions if p < pep_len
            ) / len([p for p in pocket_positions if p < pep_len])
            # Map from [-1, 1] to binding score [0, 1]
            score = (pocket_hydro + 1.0) / 2.0
        else:
            # MHC-I: average hydrophobicity at positions 1 and -1 (anchors)
            if pep_len >= 4:
                anchor_hydro = (
                    _hydro.get(peptide[1], 0.0) + _hydro.get(peptide[-1], 0.0)
                ) / 2.0
            else:
                anchor_hydro = 0.0
            # Map from [-1, 1] to binding score [0, 1]
            score = (anchor_hydro + 1.0) / 2.0

        # Apply contrast
        score = score ** _PSSM_CONTRAST_POWER

        ic50 = _binding_score_to_ic50(score)
        binding_class = _classify_binding_fallback(ic50)

        db.add(MHCBindingRecord(
            allele=allele,
            peptide=peptide,
            ic50_nm=round(ic50, 1),
            rank=None,
            binding_class=binding_class,
            source="pssm_fallback",
            method="pssm_generic_hydrophobicity",
            timestamp=timestamp,
        ))


# ═══════════════════════════════════════════════════════════════════════════
# Organism-specific default alleles
# ═══════════════════════════════════════════════════════════════════════════


def get_default_alleles_for_organism(organism: str) -> list[str]:
    """Return the default MHC alleles for a given organism.

    These are the most commonly studied alleles for immunogenicity
    screening in each organism, covering both MHC class I and class II.

    Parameters
    ----------
    organism : str
        Organism name (e.g. ``"Homo_sapiens"``, ``"Mus_musculus"``).

    Returns
    -------
    list[str]
        List of MHC allele names.  Returns an empty list for
        unrecognised organisms (with a warning).

    Examples
    --------
    >>> get_default_alleles_for_organism("Homo_sapiens")
    ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*03:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-DRB1*01:01', 'HLA-DRB1*04:01']
    >>> get_default_alleles_for_organism("Mus_musculus")
    ['H-2Kb', 'H-2Db']
    """
    import logging
    _logger = logging.getLogger(__name__)

    _DEFAULT_ALLELES: dict[str, list[str]] = {
        "Homo_sapiens": [
            "HLA-A*02:01",
            "HLA-A*01:01",
            "HLA-A*03:01",
            "HLA-B*07:02",
            "HLA-B*08:01",
            "HLA-DRB1*01:01",
            "HLA-DRB1*04:01",
        ],
        "Mus_musculus": [
            "H-2Kb",
            "H-2Db",
        ],
    }

    alleles = _DEFAULT_ALLELES.get(organism)
    if alleles is not None:
        return list(alleles)  # return a copy

    _logger.warning(
        "No default MHC alleles defined for organism %r. "
        "Returning empty list. Supported organisms: %s",
        organism,
        sorted(_DEFAULT_ALLELES.keys()),
    )
    return []
