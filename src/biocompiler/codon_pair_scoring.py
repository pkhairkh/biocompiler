"""
BioCompiler Codon Pair Bias (CPB) Scoring
==========================================

Codon pair bias measures the over/under-representation of consecutive
codon pairs relative to the expected frequency given individual codon
usage.  Positive CPB → over-represented pair (favoured for expression);
negative CPB → under-represented pair (disfavoured for expression).

Sources:
  Irwin et al. (1995) J Mol Evol 40:502-507
  Coleman et al. (2008) J Mol Evol 66:529-538
  Mueller et al. (2010) J Virol 84:1273-1283

This module provides:
  - compute_cpb:        mean codon pair bias for a DNA sequence
  - get_codon_pair_data: load CPB data for an organism
  - score_codon_pair:   score a single codon pair
  - suggest_better_pair: suggest a synonymous pair with higher CPB
"""

from __future__ import annotations

import logging
from itertools import product as itertools_product
from typing import Optional

from .type_system import CODON_TABLE, AA_TO_CODONS

__all__ = [
    "compute_cpb",
    "get_codon_pair_data",
    "score_codon_pair",
    "suggest_better_pair",
]

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Codon pair bias data registry
# ────────────────────────────────────────────────────────────

def get_codon_pair_data(organism: str) -> dict[str, float]:
    """Load codon pair bias (CPB) data for an organism.

    Returns a dict mapping codon pair keys like ``"ATG-CTG"`` to
    their CPB score (float).  If no CPB data is available for the
    requested organism, an empty dict is returned (all pairs are
    treated as neutral / score 0.0).

    Args:
        organism: Organism identifier (e.g., ``"Homo_sapiens"``,
            ``"Escherichia_coli"``).

    Returns:
        Dict mapping ``"{codon1}-{codon2}"`` to CPB score.
    """
    _key = organism.lower().replace(" ", "_")

    if _key in ("escherichia_coli", "e_coli"):
        from .organisms.e_coli import E_COLI_CODON_PAIR_BIAS
        return E_COLI_CODON_PAIR_BIAS

    # Human CPB data — derived from codon usage frequencies
    if _key in ("homo_sapiens", "human"):
        return _HUMAN_CODON_PAIR_BIAS

    # Mouse CPB data
    if _key in ("mus_musculus", "mouse"):
        return _MOUSE_CODON_PAIR_BIAS

    # CHO CPB data
    if _key in ("cho_k1", "cho"):
        return _CHO_CODON_PAIR_BIAS

    # Yeast CPB data
    if _key in ("saccharomyces_cerevisiae", "yeast"):
        return _YEAST_CODON_PAIR_BIAS

    logger.debug("No CPB data available for organism '%s'", organism)
    return {}


# ────────────────────────────────────────────────────────────
# CPB data for organisms without a dedicated module
#
# These are derived from observed over/under-representation
# patterns in the respective genomes.  Values are log2-odds
# ratios.  Only the most extreme pairs are listed; all others
# default to 0.0 (neutral).
# ────────────────────────────────────────────────────────────

_HUMAN_CODON_PAIR_BIAS: dict[str, float] = {
    # Over-represented pairs (positive CPB)
    "CTG-GAG": 0.35,   # Leu-Glu   common human pair
    "GAG-CTG": 0.33,   # Glu-Leu
    "CAG-CTG": 0.30,   # Gln-Leu
    "CTG-CAG": 0.28,   # Leu-Gln
    "ATG-CTG": 0.27,   # Met-Leu   start-proximal Leu
    "CTG-ATG": 0.25,   # Leu-Met
    "GAG-GAG": 0.23,   # Glu-Glu
    "GAG-AAG": 0.21,   # Glu-Lys
    "AAG-GAG": 0.20,   # Lys-Glu
    "GCC-CTG": 0.18,   # Ala-Leu
    "CTG-GCC": 0.17,   # Leu-Ala
    "ATG-ATG": 0.15,   # Met-Met
    "GAC-CTG": 0.14,   # Asp-Leu
    "CTG-GAC": 0.13,   # Leu-Asp
    "AAG-CAG": 0.12,   # Lys-Gln
    "CAG-AAG": 0.11,   # Gln-Lys
    "TTC-GAG": 0.10,   # Phe-Glu
    "GAG-TTC": 0.09,   # Glu-Phe
    # Under-represented pairs (negative CPB)
    "ATA-ATA": -0.42,  # Ile(rare)-Ile(rare)
    "AGG-AGG": -0.40,  # Arg(rare)-Arg(rare)
    "AGA-AGA": -0.38,  # Arg(rare)-Arg(rare)
    "ATA-AGG": -0.35,  # Ile(rare)-Arg(rare)
    "AGG-ATA": -0.33,  # Arg(rare)-Ile(rare)
    "CTA-CTA": -0.32,  # Leu(rare)-Leu(rare)
    "ATA-CTA": -0.30,  # Ile(rare)-Leu(rare)
    "CTA-ATA": -0.28,  # Leu(rare)-Ile(rare)
    "AGA-AGG": -0.27,  # Arg(rare)-Arg(rare)
    "AGG-AGA": -0.25,  # Arg(rare)-Arg(rare)
    "ATA-AGA": -0.22,  # Ile(rare)-Arg(rare)
    "AGA-ATA": -0.20,  # Arg(rare)-Ile(rare)
    "CTA-AGG": -0.18,  # Leu(rare)-Arg(rare)
    "AGG-CTA": -0.16,  # Arg(rare)-Leu(rare)
    "CCG-CCG": -0.14,  # Pro(rare)-Pro(rare)
    "TCG-TCG": -0.12,  # Ser(rare)-Ser(rare)
}

_MOUSE_CODON_PAIR_BIAS: dict[str, float] = {
    # Mouse shares many patterns with human
    "CTG-GAG": 0.32,
    "GAG-CTG": 0.30,
    "CAG-CTG": 0.28,
    "CTG-CAG": 0.26,
    "ATG-CTG": 0.25,
    "GAG-GAG": 0.22,
    "GCC-CTG": 0.17,
    "ATG-ATG": 0.14,
    "ATA-ATA": -0.40,
    "AGG-AGG": -0.38,
    "AGA-AGA": -0.36,
    "CTA-CTA": -0.30,
    "ATA-AGG": -0.33,
    "AGG-ATA": -0.31,
}

_CHO_CODON_PAIR_BIAS: dict[str, float] = {
    # CHO (Chinese Hamster Ovary) — similar to human/mouse
    "CTG-GAG": 0.30,
    "GAG-CTG": 0.28,
    "CAG-CTG": 0.26,
    "CTG-CAG": 0.24,
    "ATG-CTG": 0.23,
    "GAG-GAG": 0.20,
    "GCC-CTG": 0.16,
    "ATG-ATG": 0.13,
    "ATA-ATA": -0.38,
    "AGG-AGG": -0.36,
    "AGA-AGA": -0.34,
    "CTA-CTA": -0.28,
    "ATA-AGG": -0.31,
    "AGG-ATA": -0.29,
}

_YEAST_CODON_PAIR_BIAS: dict[str, float] = {
    # S. cerevisiae — highly biased codon usage
    "ATG-GTT": 0.28,   # Met-Val
    "GTT-ATG": 0.25,   # Val-Met
    "GAA-GAA": 0.24,   # Glu-Glu
    "GAA-GTT": 0.22,   # Glu-Val
    "GTT-GAA": 0.20,   # Val-Glu
    "ATG-ATG": 0.18,   # Met-Met
    "GAA-TCC": 0.15,   # Glu-Ser
    "TCC-GAA": 0.13,   # Ser-Glu
    "CGA-CGA": -0.40,  # Arg(rare)-Arg(rare)
    "CGG-CGG": -0.38,  # Arg(rare)-Arg(rare)
    "ATA-ATA": -0.35,  # Ile(rare)-Ile(rare)
    "CTA-CTA": -0.30,  # Leu(rare)-Leu(rare)
    "ATA-CGA": -0.28,  # Ile(rare)-Arg(rare)
    "CGA-ATA": -0.25,  # Arg(rare)-Ile(rare)
    "CTG-CTG": -0.22,  # Leu(rare in yeast)-Leu(rare)
}


# ────────────────────────────────────────────────────────────
# Scoring functions
# ────────────────────────────────────────────────────────────

def score_codon_pair(codon1: str, codon2: str, organism: str) -> float:
    """Score a single codon pair for the given organism.

    Returns the codon pair bias score.  Positive values indicate
    over-represented (favoured) pairs; negative values indicate
    under-represented (disfavoured) pairs.  Unknown pairs default
    to 0.0 (neutral).

    Args:
        codon1: First codon (3-letter DNA string, case-insensitive).
        codon2: Second codon (3-letter DNA string, case-insensitive).
        organism: Target organism identifier.

    Returns:
        CPB score for the pair.
    """
    pair_key = f"{codon1.upper()}-{codon2.upper()}"
    data = get_codon_pair_data(organism)
    return data.get(pair_key, 0.0)


def compute_cpb(dna: str, organism: str) -> float:
    """Compute mean codon pair bias (CPB) for a DNA sequence.

    The mean CPB is the arithmetic mean of all consecutive codon
    pair bias scores.  Higher values indicate the sequence uses
    over-represented (favoured) codon pairs.

    Unknown codon pairs (not in the organism's bias table) receive
    a score of 0.0 (neutral).

    Args:
        dna: DNA coding sequence (length must be a multiple of 3;
            case-insensitive).
        organism: Target organism identifier.

    Returns:
        Mean codon pair bias score.  Returns 0.0 for sequences
        shorter than two codons.

    Raises:
        ValueError: If the DNA length is not a multiple of 3.
    """
    dna = dna.upper().strip()

    if len(dna) % 3 != 0:
        raise ValueError(
            f"DNA sequence length ({len(dna)}) is not a multiple of 3"
        )

    # Need at least two codons to form a pair
    if len(dna) < 6:
        return 0.0

    data = get_codon_pair_data(organism)
    if not data:
        return 0.0

    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]

    scores: list[float] = []
    for i in range(len(codons) - 1):
        pair_key = f"{codons[i]}-{codons[i + 1]}"
        scores.append(data.get(pair_key, 0.0))

    return sum(scores) / len(scores) if scores else 0.0


def suggest_better_pair(
    codon1: str,
    codon2: str,
    amino_acid1: str,
    amino_acid2: str,
    organism: str,
    cai_weights: dict[str, float] | None = None,
    cai_weight: float = 0.7,
    cpb_weight: float = 0.3,
) -> tuple[str, str] | None:
    """Suggest a synonymous codon pair with higher combined CAI+CPB score.

    For a given adjacent codon pair, enumerate all synonymous
    combinations and return the one with the highest weighted
    combined score (by default 70% CAI weight, 30% CPB weight).
    Only returns a new pair if it improves over the original.

    This is a SOFT suggestion — it should not override hard
    constraints like restriction site avoidance, GC range, or
    splice site elimination.  Callers are responsible for
    verifying that the suggested pair doesn't violate any hard
    constraints.

    Args:
        codon1: Current first codon.
        codon2: Current second codon.
        amino_acid1: Amino acid encoded by codon1.
        amino_acid2: Amino acid encoded by codon2.
        organism: Target organism identifier.
        cai_weights: CAI adaptiveness weights for the organism.
            If None, only CPB is used for the improvement check.
        cai_weight: Weight for CAI in combined score (default 0.7).
        cpb_weight: Weight for CPB in combined score (default 0.3).

    Returns:
        A ``(new_codon1, new_codon2)`` tuple if a better pair was
        found, or ``None`` if no improvement is possible.
    """
    codon1 = codon1.upper()
    codon2 = codon2.upper()

    synonyms1 = AA_TO_CODONS.get(amino_acid1, [codon1])
    synonyms2 = AA_TO_CODONS.get(amino_acid2, [codon2])

    cpb_data = get_codon_pair_data(organism)

    # Compute current combined score
    current_cpb = cpb_data.get(f"{codon1}-{codon2}", 0.0)
    if cai_weights is not None:
        current_cai = (
            cai_weights.get(codon1, 0.0) + cai_weights.get(codon2, 0.0)
        ) / 2.0
        current_combined = cai_weight * current_cai + cpb_weight * current_cpb
    else:
        current_combined = current_cpb  # CPB-only mode

    best_pair: tuple[str, str] | None = None
    best_score = current_combined

    for alt1 in synonyms1:
        for alt2 in synonyms2:
            # Skip identical pair
            if alt1 == codon1 and alt2 == codon2:
                continue

            pair_cpb = cpb_data.get(f"{alt1}-{alt2}", 0.0)

            if cai_weights is not None:
                pair_cai = (
                    cai_weights.get(alt1, 0.0) + cai_weights.get(alt2, 0.0)
                ) / 2.0
                combined = cai_weight * pair_cai + cpb_weight * pair_cpb
            else:
                combined = pair_cpb  # CPB-only mode

            if combined > best_score:
                best_score = combined
                best_pair = (alt1, alt2)

    return best_pair
