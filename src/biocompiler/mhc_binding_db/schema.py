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
        Amino acid sequence of the peptide (8–11 residues for MHC-I).
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
    """

    allele: str
    peptide: str
    ic50_nm: float
    rank: Optional[float]
    binding_class: str
    source: str
    method: str
    timestamp: str

    def __post_init__(self) -> None:
        """Validate field values after construction."""
        if not (8 <= len(self.peptide) <= 11):
            raise ValueError(
                f"Peptide length must be 8-11 residues, got {len(self.peptide)} "
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
) -> float:
    """Score a peptide against a PSSM (preferred + disfavored rows).

    Returns a normalised binding score in [0, 1].
    """
    n = len(peptide)
    scores: list[float] = []
    for i in range(n):
        aa = peptide[i]
        # Check preferred first
        if i < len(preferred) and aa in preferred[i]:
            scores.append(preferred[i][aa])
        elif i < len(disfavored) and aa in disfavored[i]:
            scores.append(disfavored[i][aa])
        else:
            scores.append(1.0)  # neutral

    # Geometric mean
    log_sum = sum(math.log(max(s, 1e-10)) for s in scores)
    geo_mean = math.exp(log_sum / n)

    # Compute max and min possible geometric means for normalisation
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

    max_geo = math.exp(sum(math.log(max(s, 1e-10)) for s in max_scores) / n)
    min_geo = math.exp(sum(math.log(max(s, 1e-10)) for s in min_scores) / n)

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


def generate_fallback_database(
    alleles: list[str],
    peptide_lengths: list[int],
) -> MHCBindingDatabase:
    """Generate a fallback MHC binding database using PSSM scoring.

    For each allele with a known PSSM, this generates representative
    9-mer peptides by sampling amino acids with anchor-position bias,
    scores them, and stores the results.  For alleles without a known
    PSSM, a generic scoring model is applied.

    This provides OFFLINE capability when MHCflurry or NetMHCpan are
    unavailable, but is NOT as accurate as those tools (expected
    AUC-ROC 0.60–0.75).

    Parameters
    ----------
    alleles : list[str]
        HLA allele names to generate predictions for.
    peptide_lengths : list[int]
        Peptide lengths to consider.  Only 9-mers are scored with PSSMs;
        other lengths use a generic heuristic.

    Returns
    -------
    MHCBindingDatabase
        Database populated with PSSM-scored records, all marked as
        ``source="pssm_fallback"``.
    """
    db = MHCBindingDatabase()
    now = datetime.now(timezone.utc).isoformat()
    rng = random.Random(42)  # deterministic seed for reproducibility

    # Anchor positions for 9-mers: positions 1 and 8 (0-indexed)
    anchor_positions = {1, 8}

    for allele in alleles:
        preferred = _PSSM_FALLBACK.get(allele, [])
        disfavored = _PSSM_DISFAVORED.get(allele, [])

        for pep_len in peptide_lengths:
            if pep_len != 9:
                # For non-9-mer peptides, use a generic heuristic:
                # generate a small sample and score with a simple model
                _generate_generic_peptides(
                    db, allele, pep_len, rng, now,
                )
                continue

            if not preferred:
                # No PSSM data for this allele — use generic
                _generate_generic_peptides(
                    db, allele, pep_len, rng, now,
                )
                continue

            # Generate peptides with biased anchor positions
            # Sample a limited set for tractability
            num_samples = 200
            seen_peptides: set[str] = set()

            # Get anchor-preferred amino acids
            anchor_aa_pos1 = list(preferred[1].keys()) if len(preferred) > 1 else list(STANDARD_AMINO_ACIDS)
            anchor_aa_pos8 = list(preferred[8].keys()) if len(preferred) > 8 else list(STANDARD_AMINO_ACIDS)

            for _ in range(num_samples * 10):  # oversample to get enough unique
                if len(seen_peptides) >= num_samples:
                    break

                # Build peptide: bias anchor positions toward preferred AAs
                residues: list[str] = []
                for pos in range(9):
                    if pos == 1 and anchor_aa_pos1:
                        aa = rng.choice(anchor_aa_pos1 + list(STANDARD_AMINO_ACIDS))
                    elif pos == 8 and anchor_aa_pos8:
                        aa = rng.choice(anchor_aa_pos8 + list(STANDARD_AMINO_ACIDS))
                    else:
                        aa = rng.choice(list(STANDARD_AMINO_ACIDS))
                    residues.append(aa)

                peptide = "".join(residues)
                if peptide in seen_peptides:
                    continue
                seen_peptides.add(peptide)

                score = _score_peptide_pssm(peptide, preferred, disfavored)
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
                    timestamp=now,
                ))

    return db


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
    hydrophobic residues at the termini score higher, mimicking
    the general anchor-position preference of MHC-I molecules.
    """
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

        # Score: average hydrophobicity at positions 1 and -1 (anchors)
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
