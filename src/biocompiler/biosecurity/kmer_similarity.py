"""
K-mer extraction and similarity scoring, plus risk-level classification helpers.
"""

from __future__ import annotations

from biocompiler.biosecurity.types import RiskLevel

# ─────────────────────────────────────────────────────────────────────────────
# K-mer similarity helpers
# ─────────────────────────────────────────────────────────────────────────────

_KMER_SIZE: int = 5
_SIMILARITY_THRESHOLD: float = 0.6


def _extract_kmers(sequence: str, k: int = _KMER_SIZE) -> set[str]:
    """Extract all k-mers from a sequence.

    Parameters
    ----------
    sequence : str
        The amino acid or nucleotide sequence.
    k : int
        K-mer length (default: ``_KMER_SIZE``).

    Returns
    -------
    set[str]
        Set of all k-mers found in the sequence.  Returns an empty
        set if the sequence is shorter than *k*.
    """
    sequence = sequence.upper()
    if len(sequence) < k:
        return set()
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def _compute_kmer_similarity(query_kmers: set[str], pathogen_kmers: set[str]) -> float:
    """Compute the Jaccard-like k-mer similarity between two sets.

    Defined as ``|intersection| / |query_kmers|`` — the fraction of the
    query's k-mers that also appear in the pathogen signature.  Returns
    ``0.0`` if *query_kmers* is empty.

    Parameters
    ----------
    query_kmers : set[str]
        K-mers from the query sequence.
    pathogen_kmers : set[str]
        K-mers from a pathogen signature.

    Returns
    -------
    float
        Similarity score in [0, 1].
    """
    if not query_kmers:
        return 0.0
    return len(query_kmers & pathogen_kmers) / len(query_kmers)


# ─────────────────────────────────────────────────────────────────────────────
# Risk-level classification
# ─────────────────────────────────────────────────────────────────────────────

_RISK_ORDER = {"none": 0, "low": 1, "medium": 2, "high": 3, "critical": 4}


def _max_risk(*levels: str) -> RiskLevel:
    """Return the highest risk level among the given levels."""
    if not levels:
        return "none"
    best = max(levels, key=lambda l: _RISK_ORDER.get(l, 0))
    return best  # type: ignore[return-value]
