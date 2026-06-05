"""
tRNA Adaptation Index (tAI)
============================

Calculate the tRNA Adaptation Index for DNA coding sequences.

tAI is a measure of codon usage bias that accounts for tRNA gene copy
numbers (a proxy for tRNA abundance) and wobble base pairing rules.
It was introduced by dos Reis et al. (2004) and is considered a more
biologically meaningful metric than CAI because:

1. **tRNA gene copy numbers**: tAI uses the actual number of tRNA genes
   in the genome (from GtRNAdb) as a proxy for tRNA abundance, rather
   than codon frequencies which can be influenced by mutational bias.
2. **Wobble base pairing**: tAI accounts for the fact that not all
   anticodon-codon pairings are equal — some tRNAs can read multiple
   codons through wobble base pairing, with different efficiencies.
3. **Expression prediction**: tAI better predicts expression levels
   for heterologous genes because it captures the translational
   capacity of the host organism's tRNA pool.

The tAI for a gene is defined as the geometric mean of the relative
adaptiveness values of its codons:

    tAI = (∏_{i=1}^{L} w_i)^{1/L}

where w_i is the relative adaptiveness of codon i, defined as:

    w_i = ∑_{j} T_{aa,j} × s(codon, anticodon_j)

where T_{aa,j} is the gene copy number of the j-th tRNA for amino acid
aa, and s(codon, anticodon_j) is the wobble efficiency factor for the
pairing between codon and anticodon j.

Data sources
------------
tRNA gene copy numbers are from the Genomic tRNA Database (GtRNAdb):
    Chan, P.P. & Lowe, T.M. (2016) Nucleic Acids Research 44:D184-D189

References
----------
dos Reis, M., Savva, R. & Wernisch, L. (2004). Solving the riddle of
codon usage preferences: a test for translational selection.
*Nucleic Acids Research*, 32(17), 5036-5044.
doi:10.1093/nar/gkh834
"""

from __future__ import annotations

import math
from typing import Optional

from .constants import CODON_TABLE
from .organisms import resolve_organism

__all__ = [
    "calculate_tai",
    "TRNA_GENE_COPIES",
    "WOBBLE_RULES",
    "WOBBLE_EFFICIENCY",
    "SUPPORTED_ORGANISMS_TAI",
    "compute_codon_weights",
]


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Wobble Base Pairing Efficiency Factors
# ═══════════════════════════════════════════════════════════════════════════════
#
# Following dos Reis et al. (2004), the efficiency of wobble base pairing
# between the first position of the anticodon (position 34) and the third
# position of the codon (position 36) is:
#
#   Watson-Crick pairs:     A:U, U:A, G:C, C:G → s = 1.0
#   Wobble G:U pair:        G:U → s = 0.5 (modified to ~0.6 in some models)
#   Modified wobble:        Inosine pairs → s varies
#   Other non-canonical:    s = 0.0
#
# The s values below follow the dos Reis et al. (2004) model.
# Key: (anticodon_first_base, codon_third_base) → efficiency

WOBBLE_EFFICIENCY: dict[tuple[str, str], float] = {
    # Watson-Crick pairs (perfect match)
    ("A", "U"): 1.0,  # anticodon A reads codon U (A:U)
    ("U", "A"): 1.0,  # anticodon U reads codon A (U:A)
    ("G", "C"): 1.0,  # anticodon G reads codon C (G:C)
    ("C", "G"): 1.0,  # anticodon C reads codon G (C:G)

    # Wobble G:U pair (G in anticodon, U in codon)
    ("G", "U"): 0.5,

    # Modified uridines in anticodon (wobble position)
    # U reads both A and G through wobble (in RNA, U can pair with A and G)
    ("U", "G"): 0.2,  # weak U:G wobble

    # Inosine (I) — modified adenosine in anticodon
    # I can pair with U, C, or A at the codon third position
    ("I", "U"): 0.35,
    ("I", "C"): 0.65,
    ("I", "A"): 0.15,

    # Modified Uridines with specific pairing properties
    # cmo5U (uridine-5-oxyacetic acid) — reads U, A, G
    ("cmo5U", "U"): 0.45,
    ("cmo5U", "A"): 0.65,
    ("cmo5U", "G"): 0.35,

    # xm5U (5-methyluridine derivatives) — reads U, A (limited wobble)
    ("xm5U", "U"): 0.65,
    ("xm5U", "A"): 0.45,

    # k2C (lysidine) — modified C that reads A specifically
    ("k2C", "A"): 1.0,

    # All other pairings: efficiency = 0
}

# Default efficiency for unlisted pairings
_DEFAULT_WOBBLE_EFFICIENCY: float = 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Wobble Rules
# ═══════════════════════════════════════════════════════════════════════════════
#
# For each codon, which anticodons can read it, and with what efficiency.
# This is derived from the wobble efficiency table above.
# Format: codon -> list of (anticodon, efficiency)

WOBBLE_RULES: dict[str, list[tuple[str, float]]] = {
    # ── Phenylalanine ──
    "UUU": [("GAA", 1.0)],       # GAA anticodon reads UUU (Watson-Crick)
    "UUC": [("GAA", 1.0)],       # GAA also reads UUC (wobble G:U → G reads U in codon)
    # ── Leucine ──
    "UUA": [("UAA", 1.0)],       # UAA anticodon reads UUA
    "UUG": [("UAA", 0.2), ("CAA", 1.0)],  # UAA wobble + CAA Watson-Crick
    "CUU": [("GAG", 1.0)],       # GAG reads CUU
    "CUC": [("GAG", 1.0)],       # GAG reads CUC
    "CUA": [("UAG", 1.0)],       # UAG reads CUA
    "CUG": [("UAG", 0.2), ("CAG", 1.0)],  # UAG wobble + CAG Watson-Crick
    # ── Isoleucine ──
    "AUU": [("UAU", 1.0), ("CAU", 1.0)],  # UAU + CAU (modified)
    "AUC": [("UAU", 1.0), ("CAU", 1.0)],  # UAU + CAU
    "AUA": [("UAU", 0.2), ("k2C", 1.0)],  # UAU wobble + lysidine reads AUA
    # ── Methionine ──
    "AUG": [("CAU", 1.0)],       # Only CAU anticodon
    # ── Valine ──
    "GUU": [("GAC", 1.0)],       # GAC reads GUU
    "GUC": [("GAC", 1.0)],       # GAC reads GUC
    "GUA": [("UAC", 1.0)],       # UAC reads GUA
    "GUG": [("UAC", 0.2), ("CAC", 1.0)],  # UAC wobble + CAC
    # ── Serine ──
    "UCU": [("GGA", 1.0), ("IGA", 0.35)],  # GGA + IGA (Inosine)
    "UCC": [("GGA", 1.0), ("IGA", 0.65)],
    "UCA": [("IGA", 0.15), ("UGA", 1.0)],
    "UCG": [("IGA", 0.15), ("CGA", 1.0)],
    "AGU": [("GCU", 1.0)],
    "AGC": [("GCU", 1.0)],
    # ── Proline ──
    "CCU": [("GGG", 1.0)],
    "CCC": [("GGG", 1.0)],
    "CCA": [("UGG", 1.0)],
    "CCG": [("UGG", 0.2), ("CGG", 1.0)],
    # ── Threonine ──
    "ACU": [("GGU", 1.0), ("IGU", 0.35)],
    "ACC": [("GGU", 1.0), ("IGU", 0.65)],
    "ACA": [("IGU", 0.15), ("UGU", 1.0)],
    "ACG": [("IGU", 0.15), ("CGU", 1.0)],
    # ── Alanine ──
    "GCU": [("GGC", 1.0), ("IGC", 0.35)],
    "GCC": [("GGC", 1.0), ("IGC", 0.65)],
    "GCA": [("IGC", 0.15), ("UGC", 1.0)],
    "GCG": [("IGC", 0.15), ("CGC", 1.0)],
    # ── Tyrosine ──
    "UAU": [("GUA", 1.0)],
    "UAC": [("GUA", 1.0)],
    # ── Histidine ──
    "CAU": [("GUG", 1.0)],
    "CAC": [("GUG", 1.0)],
    # ── Glutamine ──
    "CAA": [("UUG", 1.0)],
    "CAG": [("UUG", 0.2), ("CUG", 1.0)],
    # ── Asparagine ──
    "AAU": [("GUU", 1.0)],
    "AAC": [("GUU", 1.0)],
    # ── Lysine ──
    "AAA": [("UUU", 1.0)],
    "AAG": [("UUU", 0.2), ("CUU", 1.0)],
    # ── Aspartate ──
    "GAU": [("GUC", 1.0)],
    "GAC": [("GUC", 1.0)],
    # ── Glutamate ──
    "GAA": [("UUC", 1.0)],
    "GAG": [("UUC", 0.2), ("CUC", 1.0)],
    # ── Cysteine ──
    "UGU": [("GCA", 1.0)],
    "UGC": [("GCA", 1.0)],
    # ── Tryptophan ──
    "UGG": [("CCA", 1.0)],
    # ── Arginine ──
    "CGU": [("GCG", 1.0), ("ICG", 0.35)],
    "CGC": [("GCG", 1.0), ("ICG", 0.65)],
    "CGA": [("ICG", 0.15), ("UCG", 1.0)],
    "CGG": [("ICG", 0.15), ("CCG", 1.0)],
    "AGA": [("UCU", 1.0)],
    "AGG": [("UCU", 0.2), ("CCU", 1.0)],
    # ── Glycine ──
    "GGU": [("GCC", 1.0)],
    "GGC": [("GCC", 1.0)],
    "GGA": [("UCC", 1.0)],
    "GGG": [("UCC", 0.2), ("CCC", 1.0)],
    # ── Stop codons ──
    "UAA": [("UUA", 1.0)],
    "UAG": [("CUA", 1.0)],
    "UGA": [("UCA", 1.0)],
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. tRNA Gene Copy Numbers
# ═══════════════════════════════════════════════════════════════════════════════
#
# Source: GtRNAdb (Genomic tRNA Database)
# Chan, P.P. & Lowe, T.M. (2016) Nucleic Acids Research 44:D184-D189
#
# Format: organism -> {anticodon: gene_copy_number}
# The anticodon is written 5'→3' (matching the wobble position at position 34).

TRNA_GENE_COPIES: dict[str, dict[str, int]] = {
    # ── Escherichia coli K-12 MG1655 ─────────────────────────────────
    # Source: GtRNAdb, E. coli K-12 MG1655 genome
    # 86 tRNA genes total
    "e_coli": {
        # Phenylalanine
        "GAA": 3,
        # Leucine
        "UAA": 1, "CAA": 6, "UAG": 1, "CAG": 3, "GAG": 2,
        # Isoleucine
        "UAU": 3, "CAU": 1,
        # Methionine (initiator + elongator)
        "CAU": 5,
        # Valine
        "GAC": 2, "UAC": 1, "CAC": 3,
        # Serine
        "GGA": 1, "CGA": 1, "UGA": 4, "GCU": 2,
        # Proline
        "GGG": 2, "UGG": 2, "CGG": 3,
        # Threonine
        "GGU": 2, "CGU": 3, "UGU": 2,
        # Alanine
        "GGC": 3, "UGC": 2, "CGC": 4,
        # Tyrosine
        "GUA": 2,
        # Histidine
        "GUG": 1,
        # Glutamine
        "UUG": 2, "CUG": 4,
        # Asparagine
        "GUU": 3,
        # Lysine
        "UUU": 2, "CUU": 3,
        # Aspartate
        "GUC": 3,
        # Glutamate
        "UUC": 2, "CUC": 4,
        # Cysteine
        "GCA": 1,
        # Tryptophan
        "CCA": 1,
        # Arginine
        "GCG": 2, "CCG": 3, "UCG": 1, "UCU": 1, "CCU": 1,
        # Glycine
        "GCC": 3, "UCC": 1, "CCC": 3,
        # Stop (release factors; no tRNAs, but included for completeness)
    },

    # ── Homo sapiens ──────────────────────────────────────────────────
    # Source: GtRNAdb, human genome (GRCh38)
    # ~610 tRNA genes total (many are pseudogenes; only functional counted)
    "human": {
        # Phenylalanine
        "GAA": 10,
        # Leucine
        "UAA": 5, "CAA": 12, "UAG": 4, "CAG": 8, "GAG": 7,
        # Isoleucine
        "UAU": 8, "CAU": 4,
        # Methionine
        "CAU": 16,
        # Valine
        "GAC": 7, "UAC": 5, "CAC": 10,
        # Serine
        "GGA": 5, "CGA": 4, "UGA": 8, "GCU": 7,
        # Proline
        "GGG": 5, "UGG": 6, "CGG": 7,
        # Threonine
        "GGU": 7, "CGU": 5, "UGU": 6,
        # Alanine
        "GGC": 10, "UGC": 6, "CGC": 8,
        # Tyrosine
        "GUA": 7,
        # Histidine
        "GUG": 5,
        # Glutamine
        "UUG": 8, "CUG": 11,
        # Asparagine
        "GUU": 10,
        # Lysine
        "UUU": 8, "CUU": 9,
        # Aspartate
        "GUC": 8,
        # Glutamate
        "UUC": 12, "CUC": 13,
        # Cysteine
        "GCA": 7,
        # Tryptophan
        "CCA": 6,
        # Arginine
        "GCG": 5, "CCG": 7, "UCG": 4, "UCU": 5, "CCU": 4,
        # Glycine
        "GCC": 10, "UCC": 4, "CCC": 8,
    },

    # ── Saccharomyces cerevisiae S288C ────────────────────────────────
    # Source: GtRNAdb, S. cerevisiae S288C
    # 275 tRNA genes total
    "yeast": {
        # Phenylalanine
        "GAA": 5,
        # Leucine
        "UAA": 4, "CAA": 8, "UAG": 3, "CAG": 5, "GAG": 4,
        # Isoleucine
        "UAU": 5, "CAU": 2,
        # Methionine
        "CAU": 8,
        # Valine
        "GAC": 5, "UAC": 3, "CAC": 6,
        # Serine
        "GGA": 4, "CGA": 3, "UGA": 6, "GCU": 5,
        # Proline
        "GGG": 4, "UGG": 4, "CGG": 5,
        # Threonine
        "GGU": 5, "CGU": 4, "UGU": 4,
        # Alanine
        "GGC": 6, "UGC": 4, "CGC": 5,
        # Tyrosine
        "GUA": 4,
        # Histidine
        "GUG": 3,
        # Glutamine
        "UUG": 5, "CUG": 7,
        # Asparagine
        "GUU": 5,
        # Lysine
        "UUU": 5, "CUU": 6,
        # Aspartate
        "GUC": 5,
        # Glutamate
        "UUC": 7, "CUC": 8,
        # Cysteine
        "GCA": 4,
        # Tryptophan
        "CCA": 3,
        # Arginine
        "GCG": 4, "CCG": 5, "UCG": 3, "UCU": 4, "CCU": 3,
        # Glycine
        "GCC": 6, "UCC": 3, "CCC": 5,
    },
}

# Canonical organism name mapping for tRNA data
_TAI_ORGANISM_ALIASES: dict[str, str] = {
    "Escherichia_coli": "e_coli",
    "e_coli": "e_coli",
    "ecoli": "e_coli",
    "E. coli": "e_coli",
    "Homo_sapiens": "human",
    "human": "human",
    "H. sapiens": "human",
    "h_sapiens": "human",
    "Saccharomyces_cerevisiae": "yeast",
    "yeast": "yeast",
    "S. cerevisiae": "yeast",
    "s_cerevisiae": "yeast",
}

SUPPORTED_ORGANISMS_TAI: list[str] = list(TRNA_GENE_COPIES.keys())


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Core tAI computation
# ═══════════════════════════════════════════════════════════════════════════════

# Epsilon floor for zero-adaptiveness codons
_TAI_EPSILON: float = 1e-10

# Minimum copy number to avoid division by zero
_MIN_COPY_NUMBER: int = 1


def _resolve_tai_organism(organism: str) -> str:
    """Resolve an organism name to a tAI data key.

    Parameters
    ----------
    organism : str
        Any organism name or alias.

    Returns
    -------
    str
        Key in TRNA_GENE_COPIES.

    Raises
    ------
    ValueError
        If no tRNA data is available for the organism.
    """
    # Try direct alias resolution
    key = _TAI_ORGANISM_ALIASES.get(organism)
    if key is not None and key in TRNA_GENE_COPIES:
        return key

    # Try resolve_organism then alias
    resolved = resolve_organism(organism, strict=False)
    key = _TAI_ORGANISM_ALIASES.get(resolved)
    if key is not None and key in TRNA_GENE_COPIES:
        return key

    # Try lowercase
    key = _TAI_ORGANISM_ALIASES.get(organism.lower())
    if key is not None and key in TRNA_GENE_COPIES:
        return key

    # Direct lookup
    if organism in TRNA_GENE_COPIES:
        return organism

    raise ValueError(
        f"No tRNA gene copy data available for organism '{organism}'. "
        f"Available organisms: {list(TRNA_GENE_COPIES.keys())}"
    )


def _dna_to_rna(dna: str) -> str:
    """Convert DNA sequence to RNA (T → U)."""
    return dna.upper().replace("T", "U")


def _compute_codon_weight(
    rna_codon: str,
    trna_copies: dict[str, int],
) -> float:
    """Compute the raw weight for a single codon.

    The weight is the sum over all anticodons that can read this codon,
    of (tRNA_gene_copies × wobble_efficiency):

        W(codon) = Σ_j T(j) × s(codon, anticodon_j)

    Parameters
    ----------
    rna_codon : str
        Codon in RNA format (U instead of T).
    trna_copies : dict
        tRNA gene copy numbers {anticodon: count}.

    Returns
    -------
    float
        Raw weight for this codon.
    """
    wobble_list = WOBBLE_RULES.get(rna_codon, [])
    total_weight = 0.0

    for anticodon, efficiency in wobble_list:
        copies = trna_copies.get(anticodon, 0)
        total_weight += copies * efficiency

    return total_weight


def compute_codon_weights(
    organism: str,
) -> dict[str, float]:
    """Compute relative adaptiveness weights for all codons for an organism.

    For each codon, the relative adaptiveness is:

        w(codon) = W(codon) / W_max(amino_acid)

    where W(codon) is the raw weight and W_max is the maximum raw weight
    among all synonymous codons for the same amino acid.

    Parameters
    ----------
    organism : str
        Organism name.

    Returns
    -------
    dict[str, float]
        Mapping of RNA codon → relative adaptiveness value in [0, 1].
    """
    tai_key = _resolve_tai_organism(organism)
    trna_copies = TRNA_GENE_COPIES[tai_key]

    # Compute raw weights for all codons
    raw_weights: dict[str, float] = {}
    for codon in WOBBLE_RULES:
        raw_weights[codon] = _compute_codon_weight(codon, trna_copies)

    # Group codons by amino acid and normalize
    aa_to_codons: dict[str, list[str]] = {}
    for codon in WOBBLE_RULES:
        # Convert RNA codon to DNA for CODON_TABLE lookup
        dna_codon = codon.replace("U", "T")
        aa = CODON_TABLE.get(dna_codon)
        if aa is None or aa == "*":
            continue
        aa_to_codons.setdefault(aa, []).append(codon)

    relative_weights: dict[str, float] = {}
    for aa, codons in aa_to_codons.items():
        if aa == "M":  # Only one codon, skip
            continue
        weights = [raw_weights.get(c, 0.0) for c in codons]
        max_w = max(weights) if weights else 0.0
        for codon, w in zip(codons, weights):
            if max_w > 0:
                relative_weights[codon] = w / max_w
            else:
                relative_weights[codon] = 0.0

    return relative_weights


def calculate_tai(
    dna: str,
    organism: str,
    *,
    skip_stop: bool = True,
    skip_met: bool = True,
) -> float:
    """Calculate the tRNA Adaptation Index (tAI) for a DNA sequence.

    tAI is the geometric mean of relative adaptiveness values for all
    codons in the sequence (excluding Met and stop codons by default).

    The tAI differs from CAI in that it uses tRNA gene copy numbers
    (a proxy for tRNA abundance) and wobble pairing rules, rather than
    codon frequency data from highly expressed genes.

    Parameters
    ----------
    dna : str
        DNA coding sequence (length must be a multiple of 3).
    organism : str
        Target organism name. Must have tRNA data in TRNA_GENE_COPIES.
        Currently supports: ``"e_coli"``, ``"human"``, ``"yeast"``
        (and their canonical/alias names).
    skip_stop : bool
        Whether to exclude stop codons from the calculation (default True).
    skip_met : bool
        Whether to exclude the Met (ATG) codon from the calculation
        (default True, following the CAI convention).

    Returns
    -------
    float
        tAI value in [0, 1].  Returns 0.0 for empty or invalid sequences.

    Raises
    ------
    ValueError
        If the DNA length is not a multiple of 3.
    ValueError
        If no tRNA data is available for the organism.

    References
    ----------
    dos Reis, M., Savva, R. & Wernisch, L. (2004). Solving the riddle of
    codon usage preferences: a test for translational selection.
    *Nucleic Acids Research*, 32(17), 5036-5044.
    """
    if not dna or len(dna) < 3:
        return 0.0

    dna = dna.upper().strip()
    if len(dna) % 3 != 0:
        raise ValueError(
            f"DNA sequence length ({len(dna)}) is not a multiple of 3"
        )

    # Resolve organism and get tRNA data
    tai_key = _resolve_tai_organism(organism)
    trna_copies = TRNA_GENE_COPIES[tai_key]

    # Compute codon weights
    weights = compute_codon_weights(organism)

    # Compute tAI as geometric mean
    log_sum: float = 0.0
    count: int = 0

    for i in range(0, len(dna), 3):
        dna_codon = dna[i:i + 3]
        rna_codon = _dna_to_rna(dna_codon)

        aa = CODON_TABLE.get(dna_codon)
        if aa is None:
            continue
        if skip_stop and aa == "*":
            continue
        if skip_met and aa == "M":
            continue

        # Get relative adaptiveness
        w = weights.get(rna_codon, 0.0)
        if w <= 0:
            w = _TAI_EPSILON
        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0

    tai = math.exp(log_sum / count)
    return round(tai, 4)
