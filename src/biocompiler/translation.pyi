"""Type stubs for biocompiler.translation — public API surface."""

from __future__ import annotations

from typing import Any


# ────────────────────────────────────────────────────────────
# Public functions
# ────────────────────────────────────────────────────────────

def translate(
    sequence: str,
    table: int = ...,
    to_stop: bool = ...,
) -> str:
    """Translate a DNA sequence to a protein string.

    Args:
        sequence: DNA coding sequence.
        table: Translation table number (NCBI standard).
        to_stop: If True, stop at the first stop codon.

    Returns:
        Protein string (single-letter amino acid codes).
    """
    ...


def compute_cai(
    sequence: str,
    organism: str = ...,
    species: str | None = ...,
) -> float:
    """Compute Codon Adaptation Index (CAI) for a coding sequence.

    CAI = geometric mean of relative adaptiveness values of codons used.
    This is a DETERMINISTIC computation: same sequence → same CAI.

    Args:
        sequence: DNA coding sequence (starts with ATG).
        organism: Target organism name (canonical, short key, or alias).
            All forms are resolved via resolve_organism().
        species: Alias for ``organism`` (deprecated, emits DeprecationWarning).

    Returns:
        CAI value in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.

    Raises:
        UnsupportedOrganismError: If the organism is not supported.
    """
    ...


def find_orfs(
    sequence: str,
    min_length: int = ...,
) -> list[dict[str, Any]]:
    """Find open reading frames in a DNA sequence.

    Args:
        sequence: DNA sequence to search.
        min_length: Minimum ORF length in codons.

    Returns:
        List of ORF dictionaries with keys: start, end, length, protein.
    """
    ...
