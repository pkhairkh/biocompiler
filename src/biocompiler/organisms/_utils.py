"""
Internal utilities for computing derived codon data from raw usage tables.

Eliminates duplicated computation logic across organism modules.
"""

from __future__ import annotations

from typing import TypeAlias

# Type alias for the codon usage table format:
# {codon: (amino_acid, fraction, per_thousand, count)}
CodonUsageTable: TypeAlias = dict[str, tuple[str, float, float, int]]

# Stop codon indicator used in codon usage tables.
# Codons mapping to "*" are stop codons and excluded from adaptiveness/preferred computations.
STOP_CODON_AA: str = "*"


def compute_codon_adaptiveness(
    usage: CodonUsageTable,
) -> dict[str, float]:
    """Compute relative codon adaptiveness values from a codon usage table.

    For each amino acid, the most frequently used codon gets adaptiveness 1.0,
    and other synonymous codons are scaled proportionally.

    Args:
        usage: Codon usage table mapping codon strings to
            (amino_acid, fraction, per_thousand, count) tuples.

    Returns:
        Dict mapping codon strings to relative adaptiveness values (0.0–1.0).
        Stop codons are excluded.
    """
    # Find the maximum per-thousand frequency for each amino acid
    aa_max_freq: dict[str, float] = {}
    for _codon, (aa, _frac, freq, _count) in usage.items():
        if aa == STOP_CODON_AA:
            continue
        current = aa_max_freq.get(aa, 0.0)
        if freq > current:
            aa_max_freq[aa] = freq

    # Compute adaptiveness as freq / max_freq for same amino acid
    adaptiveness: dict[str, float] = {}
    for codon, (aa, _frac, freq, _count) in usage.items():
        if aa == STOP_CODON_AA:
            continue
        max_freq = aa_max_freq[aa]
        adaptiveness[codon] = freq / max_freq if max_freq > 0.0 else 0.0

    return adaptiveness


def compute_preferred_codons(
    usage: CodonUsageTable,
) -> dict[str, str]:
    """Compute the preferred (highest-frequency) codon for each amino acid.

    Args:
        usage: Codon usage table mapping codon strings to
            (amino_acid, fraction, per_thousand, count) tuples.

    Returns:
        Dict mapping single-letter amino acid codes to the preferred codon.
        Stop codons are excluded.
    """
    # Group codons by amino acid with their frequencies
    aa_codons: dict[str, list[tuple[str, float]]] = {}
    for codon, (aa, _frac, freq, _count) in usage.items():
        if aa == STOP_CODON_AA:
            continue
        aa_codons.setdefault(aa, []).append((codon, freq))

    # Select the highest-frequency codon per amino acid
    preferred: dict[str, str] = {}
    for aa, codons in aa_codons.items():
        preferred[aa] = max(codons, key=lambda x: x[1])[0]

    return preferred
