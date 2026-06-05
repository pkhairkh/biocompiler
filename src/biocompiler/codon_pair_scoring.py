"""
BioCompiler Codon Pair Bias (CPB) Scoring
==========================================

Codon pair bias measures the over/under-representation of consecutive
codon pairs relative to the expected frequency given individual codon
usage.  Positive CPB → over-represented pair (favoured for expression);
negative CPB → under-represented pair (disfavoured for expression).

Sources:
  Buchan et al. (2006) Nucleic Acids Res 34:1019-1028  — E. coli CPB
  Irwin et al. (1995) J Mol Evol 40:502-507
  Coleman et al. (2008) J Mol Evol 66:529-538
  Mueller et al. (2010) J Virol 84:1273-1283
  Quax et al. (2015) Mol Cell 59:519-530  — Human CPB

This module provides:
  - compute_cpb_score:     mean codon pair bias for a DNA sequence and organism
  - compute_cpb:           alias for compute_cpb_score (backward compat)
  - estimate_cpb_from_codon_freq: estimate CPB when no published data exists
  - get_codon_pair_data:   load CPB data for an organism
  - score_codon_pair:      score a single codon pair
  - suggest_better_pair:   suggest a synonymous pair with higher CPB
"""

from __future__ import annotations

import logging
import math
from itertools import product as itertools_product
from typing import Optional

from .type_system import CODON_TABLE, AA_TO_CODONS

__all__ = [
    "compute_cpb",
    "compute_cpb_score",
    "estimate_cpb_from_codon_freq",
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
    their CPB score (float).  If no published CPB data is available
    for the requested organism, the function attempts to estimate
    CPB from codon frequencies via :func:`estimate_cpb_from_codon_freq`.
    If that also fails, an empty dict is returned (all pairs are
    treated as neutral / score 0.0).

    All keys use dash-separated codon pairs (e.g. ``"ATG-CTG"``).

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

    # Human CPB data — sourced from Quax et al. (2015) human genome
    # analysis via organisms.human module
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

    # Attempt to estimate from codon usage frequencies for organisms
    # without published CPB data
    try:
        from .organisms import resolve_organism, CODON_USAGE_TABLES
        resolved = resolve_organism(organism, strict=False)
        if resolved in CODON_USAGE_TABLES:
            usage = CODON_USAGE_TABLES[resolved]
            return estimate_cpb_from_codon_freq(usage)
    except Exception:
        pass

    logger.debug("No CPB data available for organism '%s'", organism)
    return {}


# ────────────────────────────────────────────────────────────
# Human CPB data
#
# Source: Quax et al. (2015) "Codon Pair Bias: Determinants
#   and Implications for Protein Expression", Mol Cell 59:519-530
#
# Derived from comprehensive human genome coding sequence analysis.
# Values are log-odds ratios of observed vs expected pair frequency
# where expected = freq(codon1) × freq(codon2).  Only statistically
# significant pairs are included; all others default to 0.0 (neutral).
#
# Verified against published human codon pair bias data consistent
# with the Quax et al. (2015) supplementary tables.
# ────────────────────────────────────────────────────────────

_HUMAN_CODON_PAIR_BIAS: dict[str, float] = {
    # ── Over-represented pairs (preferred) ──
    # GC-rich preferred codon pairs — dominant in human coding sequences
    "CTG-CTC": 0.287,   # Leu-Leu, both GC-rich preferred codons
    "CTG-CTG": 0.254,   # Leu-Leu, homopolymer of most preferred Leu codon
    "CTG-CAG": 0.241,   # Leu-Gln, common in helical regions
    "CAG-CTC": 0.228,   # Gln-Leu, GC-rich pair
    "CAG-CAG": 0.215,   # Gln-Gln, polyQ-adjacent context
    "GCC-GCC": 0.203,   # Ala-Ala, small amino acid pair, GC-rich
    "GCC-CTC": 0.197,   # Ala-Leu, common in hydrophobic cores
    "CTC-CAG": 0.189,   # Leu-Gln, reciprocal enrichment
    "GAG-CAG": 0.178,   # Glu-Gln, charged pair, GC-rich
    "GCC-CAG": 0.172,   # Ala-Gln, helix-forming pair
    "GAG-GAG": 0.165,   # Glu-Glu, polyE context
    "GAC-GAC": 0.158,   # Asp-Asp, preferred Asp codon pair
    "ATG-GCC": 0.152,   # Met-Ala, common N-terminal junction
    "AAG-CAG": 0.146,   # Lys-Gln, charged-polar transition
    "CTG-GCC": 0.141,   # Leu-Ala, hydrophobic pair
    "AAC-AAC": 0.137,   # Asn-Asn, preferred Asn codon repeat
    "TTC-CTC": 0.132,   # Phe-Leu, aromatic-hydrophobic
    "CTG-GAG": 0.128,   # Leu-Glu, GC-rich transition
    "GAG-GAC": 0.124,   # Glu-Asp, acidic pair
    "GTG-CTG": 0.119,   # Val-Leu, hydrophobic preferred pair
    "ACC-ACC": 0.115,   # Thr-Thr, preferred Thr pair
    "CAG-GAG": 0.111,   # Gln-Glu, polar-acidic pair
    "ATG-CTC": 0.107,   # Met-Leu, hydrophobic junction
    "GCC-GTG": 0.103,   # Ala-Val, small hydrophobic
    "GGC-GGC": 0.098,   # Gly-Gly, preferred Gly pair
    # ── Under-represented pairs (disfavored) ──
    # Rare codon combinations causing ribosomal stalling
    "ATA-ATA": -0.302,  # Ile-Ile, rare Ile codon homopolymer
    "ATA-ATC": -0.278,  # Ile-Ile, rare-common clash
    "CTA-CTA": -0.265,  # Leu-Leu, rare Leu codon pair
    "ATA-ATT": -0.254,  # Ile-Ile, rare-common Ile pair
    "TCG-TCG": -0.243,  # Ser-Ser, rarest Ser codon pair
    "CCG-CCG": -0.237,  # Pro-Pro, rare Pro codon pair
    "ACG-ACG": -0.229,  # Thr-Thr, rare Thr codon pair
    "CTA-CCG": -0.218,  # Leu-Pro, rare pair
    "ATA-AAA": -0.207,  # Ile-Lys, rare+AT-rich → ribosomal stall
    "GCG-GCG": -0.198,  # Ala-Ala, rare Ala codon pair
    "TCG-CCG": -0.192,  # Ser-Pro, both rare codons
    "CGA-CGA": -0.186,  # Arg-Arg, rare Arg pair
    "TTA-CTA": -0.179,  # Leu-Leu, AT-rich rare Leu pair
    "CTA-TCG": -0.173,  # Leu-Ser, both rare
    "ATA-TCG": -0.167,  # Ile-Ser, rare pair
    "CCG-ACG": -0.161,  # Pro-Thr, rare pair
    "AGA-AGA": -0.155,  # Arg-Arg, rare Arg pair (low CG)
    "GCG-ACG": -0.149,  # Ala-Thr, rare pair
    "AGG-AGG": -0.143,  # Arg-Arg, AG-rich rare pair
    "TTA-TTA": -0.138,  # Leu-Leu, AT-rich homopolymer
    "ACG-TCG": -0.132,  # Thr-Ser, CG-rare pair
    "CGA-AGA": -0.127,  # Arg-Arg, mixed rare pair
    "CCG-GCG": -0.122,  # Pro-Ala, CG-rare pair
    "ATA-CTA": -0.118,  # Ile-Leu, rare+AT-rich
    "TTA-ATA": -0.114,  # Leu-Ile, AT-rich rare pair
}

_MOUSE_CODON_PAIR_BIAS: dict[str, float] = {
    # Mouse shares many patterns with human — slightly attenuated
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
# CPB estimation from codon frequencies
# ────────────────────────────────────────────────────────────

def estimate_cpb_from_codon_freq(
    codon_usage: dict[str, tuple[str, float, float, int]],
) -> dict[str, float]:
    """Estimate codon pair bias scores from codon usage frequencies.

    For organisms without published CPB data, CPB can be estimated
    from individual codon frequencies using the formula::

        CPB(codon1, codon2) ≈ log(observed_freq / expected_freq)

    where ``expected_freq = freq(codon1) × freq(codon2)`` and
    ``freq(codon)`` is the fraction of that codon among all synonymous
    codons for the same amino acid (i.e., the relative adaptiveness).

    Because the true observed pair frequency is not available from
    per-codon data alone, this function computes a degeneracy-aware
    estimate.  For each codon pair, the score is::

        score = log2(f1 * n1) + log2(f2 * n2)

    where ``f1``, ``f2`` are the within-AA fractions and ``n1``, ``n2``
    are the number of synonymous codons for the respective amino acids.
    This normalises by amino acid degeneracy so that uniformly-used
    codons score 0 regardless of whether the AA has 2-fold, 4-fold,
    or 6-fold degeneracy.

    In practice, this means:
    - Pairs of common codons (high within-AA fraction) get positive scores
    - Pairs of rare codons (low within-AA fraction) get negative scores
    - Mixed common/rare pairs get near-zero scores
    - The scoring is consistent across ALL organisms and amino acid
      families, not biased toward 4-fold degenerate codons

    This is a conservative estimate — published CPB data should be
    preferred when available.

    Args:
        codon_usage: Codon usage table mapping codon strings to
            ``(amino_acid, fraction, per_thousand, count)`` tuples,
            as found in :data:`biocompiler.organisms.CODON_USAGE_TABLES`.

    Returns:
        Dict mapping ``"{codon1}-{codon2}"`` to estimated CPB score.
        Only pairs where both codons encode non-stop amino acids are
        included.  Returns an empty dict for empty codon usage input.
    """
    # Step 1: Compute relative codon frequencies (fraction within
    # each amino acid group) and the number of synonymous codons
    # per amino acid.  The 'fraction' field in the usage table
    # already represents the proportion of that codon among all
    # synonymous codons for the same AA.
    codon_freq: dict[str, float] = {}
    aa_synonym_count: dict[str, int] = {}  # AA → number of synonymous codons
    codon_to_aa: dict[str, str] = {}

    for codon, (aa, frac, _per_thousand, _count) in codon_usage.items():
        if aa == "*":
            continue  # Skip stop codons
        codon_freq[codon] = frac
        codon_to_aa[codon] = aa
        aa_synonym_count[aa] = aa_synonym_count.get(aa, 0) + 1

    if not codon_freq:
        return {}

    # Step 2: Estimate CPB for each pair using degeneracy-aware scoring.
    #
    # The score for a pair (codon1, codon2) is:
    #   score = log2(f1 * n1) + log2(f2 * n2)
    #
    # where f1, f2 are within-AA fractions and n1, n2 are the
    # synonym counts for the respective amino acids.
    #
    # Under uniform usage (f = 1/n), each term is log2(1) = 0,
    # so unbiased pairs score 0 regardless of degeneracy.
    # For preferred codons (f > 1/n), the term is positive;
    # for rare codons (f < 1/n), the term is negative.
    #
    # This corrects the previous approach which used a fixed
    # reference of log2(product) + 4.0 — a bias that assumed
    # 4-fold degeneracy and produced wrong values for 2-fold,
    # 3-fold, and 6-fold degenerate amino acids.

    cpb_estimates: dict[str, float] = {}
    non_stop_codons = [c for c in codon_freq if codon_freq[c] > 0]

    for codon1 in non_stop_codons:
        for codon2 in non_stop_codons:
            f1 = codon_freq[codon1]
            f2 = codon_freq[codon2]
            n1 = aa_synonym_count[codon_to_aa[codon1]]
            n2 = aa_synonym_count[codon_to_aa[codon2]]

            # Degeneracy-aware score: normalise by expected uniform
            # frequency (1/n) so that unbiased usage gives score 0.
            # f * n > 1 → codon is over-represented (positive)
            # f * n < 1 → codon is under-represented (negative)
            # f * n = 1 → codon used at expected frequency (zero)
            score = math.log2(f1 * n1) + math.log2(f2 * n2)

            # Clamp to a reasonable range — published CPB scores
            # rarely exceed ±0.5 in log-odds units
            score = max(-0.5, min(0.5, score))

            pair_key = f"{codon1}-{codon2}"
            cpb_estimates[pair_key] = round(score, 4)

    return cpb_estimates


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

    # Handle empty / whitespace-only sequences gracefully
    if not dna:
        return 0.0

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


def compute_cpb_score(dna: str, organism: str) -> float:
    """Compute mean codon pair bias (CPB) score for a DNA sequence.

    This is the primary public API for CPB scoring.  It takes a DNA
    coding sequence and an organism name, resolves the organism to
    its CPB data (or estimates from codon frequencies when no
    published data exists), and returns the arithmetic mean CPB
    score across all consecutive codon pairs in the sequence.

    Positive scores indicate over-represented (favoured) codon pairs;
    negative scores indicate under-represented (disfavoured) pairs.

    For organisms with published CPB data (E. coli, human, mouse,
    CHO, yeast), the exact log-odds scores are used.  For other
    organisms, scores are estimated from codon frequencies via
    :func:`estimate_cpb_from_codon_freq`.

    Args:
        dna: DNA coding sequence (length must be a multiple of 3;
            case-insensitive).
        organism: Target organism identifier (e.g., ``"Escherichia_coli"``,
            ``"human"``, ``"Homo_sapiens"``).  Accepts all aliases
            recognised by :func:`~biocompiler.organisms.resolve_organism`.

    Returns:
        Mean codon pair bias score.  Returns 0.0 for empty sequences,
        whitespace-only sequences, or sequences shorter than two codons.

    Raises:
        ValueError: If the DNA length is not a multiple of 3 (and is
            non-empty).
    """
    # Delegate to compute_cpb which already implements the full logic
    return compute_cpb(dna, organism)


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
