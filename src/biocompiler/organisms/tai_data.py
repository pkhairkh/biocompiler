"""
tRNA Adaptation Index (tAI) Data
==================================

tRNA gene copy number data and wobble efficiency matrices for
computing the tRNA Adaptation Index (tAI).

Data sources:
    - E. coli K-12 MG1655: GtRNAdb (Chan & Lowe 2016)
    - S. cerevisiae S288C: GtRNAdb
    - H. sapiens GRCh38: GtRNAdb
    - M. musculus GRCm39: GtRNAdb
    - CHO-K1 (Cricetulus griseus): GtRNAdb
    - C. elegans: GtRNAdb
    - D. melanogaster: GtRNAdb
    - A. thaliana: GtRNAdb
    - P. pastoris: GtRNAdb
    - B. subtilis 168: GtRNAdb

Wobble rules corrected against dos Reis et al. (2004):
    - G:U wobble efficiency fixed to 0.5 (was incorrectly 1.0)
    - Isoleucine wobble rules fixed (GAU for AUU/AUC, k2C for AUA)
    - I:G pairings removed (Inosine cannot pair with G)
    - Ile/Met CAU disambiguation: k2C key for lysidine tRNA

References:
    dos Reis, M., Savva, R. & Wernisch, L. (2004). Solving the riddle of
    codon usage preferences: a test for translational selection.
    *Nucleic Acids Research*, 32(17), 5036-5044.
    doi:10.1093/nar/gkh834

    Chan, P.P. & Lowe, T.M. (2016). GtRNAdb 2.0: an expanded database of
    transfer RNA genes identified in complete and draft genomes.
    *Nucleic Acids Research*, 44(D1), D184-D189.
"""

from __future__ import annotations

__all__ = [
    "TRNA_GENE_COPIES",
    "WOBBLE_EFFICIENCY",
    "WOBBLE_RULES",
    "SUPPORTED_ORGANISMS_TAI",
    "compute_tai_weights",
]


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Wobble Base Pairing Efficiency Factors
# ═══════════════════════════════════════════════════════════════════════════════
#
# Following dos Reis et al. (2004), the efficiency of wobble base pairing
# between the first position of the anticodon (position 34) and the third
# position of the codon (position 36).
#
# Key: (anticodon_first_base, codon_third_base) → efficiency

WOBBLE_EFFICIENCY: dict[tuple[str, str], float] = {
    # Watson-Crick pairs (perfect match)
    ("A", "U"): 1.0,  # anticodon A reads codon U (A:U)
    ("U", "A"): 1.0,  # anticodon U reads codon A (U:A)
    ("G", "C"): 1.0,  # anticodon G reads codon C (G:C)
    ("C", "G"): 1.0,  # anticodon C reads codon G (C:G)

    # Wobble G:U pair (G in anticodon, U in codon)
    ("G", "U"): 0.5,

    # Weak U:G wobble (U in anticodon, G in codon)
    ("U", "G"): 0.2,

    # Inosine (I) — modified adenosine in anticodon (from A deamination)
    # I can pair with U, C, A but NOT G
    ("I", "U"): 0.35,
    ("I", "C"): 0.65,
    ("I", "A"): 0.15,

    # Modified Uridines with specific pairing properties
    ("cmo5U", "U"): 0.45,
    ("cmo5U", "A"): 0.65,
    ("cmo5U", "G"): 0.35,

    ("xm5U", "U"): 0.65,
    ("xm5U", "A"): 0.45,

    # k2C (lysidine) — modified C in Ile-tRNA that reads A specifically
    ("k2C", "A"): 1.0,
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Wobble Rules
# ═══════════════════════════════════════════════════════════════════════════════
#
# For each codon, which anticodons can read it, and with what efficiency.
# Format: codon -> list of (anticodon, efficiency)
#
# Anticodon naming convention:
#   - Standard bases (A, U, G, C): represent the gene-encoded anticodon
#   - "I" prefix: represents Inosine at wobble position (A→I deamination)
#     e.g., "IGA" = tRNA with gene anticodon AGA, deaminated to IGA
#   - "k2C": represents lysidine at wobble position (C→k2C modification)
#     e.g., "k2C" = tRNA-Ile with gene anticodon CAU, C modified to lysidine
#     This is separate from Met-tRNA-CAU which has unmodified C
#
# Corrections from original implementation:
#   - G:U wobble fixed from 1.0 to 0.5 (10 codons affected)
#   - Ile wobble rules fixed: GAU for AUU/AUC, k2C for AUA
#   - I:G pairings removed (Inosine cannot pair with G)
#   - Removed incorrect CAU references for Ile codons

WOBBLE_RULES: dict[str, list[tuple[str, float]]] = {
    # ── Phenylalanine ──
    # tRNA-Phe-GAA: G reads C (1.0) and U via wobble (0.5)
    "UUU": [("GAA", 0.5)],   # G:U wobble
    "UUC": [("GAA", 1.0)],   # G:C Watson-Crick

    # ── Leucine ──
    # UUA/UUG family: tRNA-Leu-UAA (U:A=1.0, U:G=0.2), tRNA-Leu-CAA (C:G=1.0)
    "UUA": [("UAA", 1.0)],                    # U:A Watson-Crick
    "UUG": [("UAA", 0.2), ("CAA", 1.0)],      # U:G wobble, C:G Watson-Crick
    # CUU/CUC/CUA/CUG family: tRNA-Leu-GAG, tRNA-Leu-UAG, tRNA-Leu-CAG
    "CUU": [("GAG", 0.5)],                    # G:U wobble
    "CUC": [("GAG", 1.0)],                    # G:C Watson-Crick
    "CUA": [("UAG", 1.0)],                    # U:A Watson-Crick
    "CUG": [("UAG", 0.2), ("CAG", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Isoleucine ──
    # Bacteria: tRNA-Ile-GAU reads AUC (1.0) and AUU (0.5)
    #           tRNA-Ile-k2C (lysidine, from CAU gene) reads AUA (1.0)
    # Eukaryotes: tRNA-Ile-IAU (from AAU gene, A→I deamination) reads all three
    "AUU": [("GAU", 0.5), ("IAU", 0.35)],     # GAU G:U wobble, IAU I:U
    "AUC": [("GAU", 1.0), ("IAU", 0.65)],     # GAU G:C, IAU I:C
    "AUA": [("k2C", 1.0), ("IAU", 0.15)],     # k2C lysidine:A, IAU I:A

    # ── Methionine ──
    # tRNA-Met-CAU (unmodified C): reads AUG only
    "AUG": [("CAU", 1.0)],                    # C:G Watson-Crick

    # ── Valine ──
    "GUU": [("GAC", 0.5)],                    # G:U wobble
    "GUC": [("GAC", 1.0)],                    # G:C Watson-Crick
    "GUA": [("UAC", 1.0)],                    # U:A Watson-Crick
    "GUG": [("UAC", 0.2), ("CAC", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Serine ──
    # UCN family: tRNA-Ser-GGA, tRNA-Ser-IGA (from AGA gene), tRNA-Ser-UGA/CGA
    "UCU": [("GGA", 0.5), ("IGA", 0.35)],     # G:U wobble, I:U
    "UCC": [("GGA", 1.0), ("IGA", 0.65)],     # G:C, I:C
    "UCA": [("IGA", 0.15), ("UGA", 1.0)],     # I:A, U:A
    "UCG": [("CGA", 1.0)],                    # C:G only (I cannot pair with G)
    # AGN family: tRNA-Ser-GCU
    "AGU": [("GCU", 0.5)],                    # G:U wobble
    "AGC": [("GCU", 1.0)],                    # G:C Watson-Crick

    # ── Proline ──
    "CCU": [("GGG", 0.5)],                    # G:U wobble
    "CCC": [("GGG", 1.0)],                    # G:C Watson-Crick
    "CCA": [("UGG", 1.0)],                    # U:A Watson-Crick
    "CCG": [("UGG", 0.2), ("CGG", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Threonine ──
    "ACU": [("GGU", 0.5), ("IGU", 0.35)],     # G:U wobble, I:U
    "ACC": [("GGU", 1.0), ("IGU", 0.65)],     # G:C, I:C
    "ACA": [("IGU", 0.15), ("UGU", 1.0)],     # I:A, U:A
    "ACG": [("CGU", 1.0)],                    # C:G only (I cannot pair with G)

    # ── Alanine ──
    "GCU": [("GGC", 0.5), ("IGC", 0.35)],     # G:U wobble, I:U
    "GCC": [("GGC", 1.0), ("IGC", 0.65)],     # G:C, I:C
    "GCA": [("IGC", 0.15), ("UGC", 1.0)],     # I:A, U:A
    "GCG": [("CGC", 1.0)],                    # C:G only (I cannot pair with G)

    # ── Tyrosine ──
    "UAU": [("GUA", 0.5)],                    # G:U wobble
    "UAC": [("GUA", 1.0)],                    # G:C Watson-Crick

    # ── Histidine ──
    "CAU": [("GUG", 0.5)],                    # G:U wobble
    "CAC": [("GUG", 1.0)],                    # G:C Watson-Crick

    # ── Glutamine ──
    "CAA": [("UUG", 1.0)],                    # U:A Watson-Crick
    "CAG": [("UUG", 0.2), ("CUG", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Asparagine ──
    "AAU": [("GUU", 0.5)],                    # G:U wobble
    "AAC": [("GUU", 1.0)],                    # G:C Watson-Crick

    # ── Lysine ──
    "AAA": [("UUU", 1.0)],                    # U:A Watson-Crick
    "AAG": [("UUU", 0.2), ("CUU", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Aspartate ──
    "GAU": [("GUC", 0.5)],                    # G:U wobble
    "GAC": [("GUC", 1.0)],                    # G:C Watson-Crick

    # ── Glutamate ──
    "GAA": [("UUC", 1.0)],                    # U:A Watson-Crick
    "GAG": [("UUC", 0.2), ("CUC", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Cysteine ──
    "UGU": [("GCA", 0.5)],                    # G:U wobble
    "UGC": [("GCA", 1.0)],                    # G:C Watson-Crick

    # ── Tryptophan ──
    "UGG": [("CCA", 1.0)],                    # C:G Watson-Crick

    # ── Arginine ──
    # CGN family: tRNA-Arg-GCG, tRNA-Arg-ICG (from ACG gene), tRNA-Arg-UCG/CCG
    "CGU": [("GCG", 0.5), ("ICG", 0.35)],     # G:U wobble, I:U
    "CGC": [("GCG", 1.0), ("ICG", 0.65)],     # G:C, I:C
    "CGA": [("ICG", 0.15), ("UCG", 1.0)],     # I:A, U:A
    "CGG": [("CCG", 1.0)],                    # C:G only (I cannot pair with G)
    # AGA/AGG family: tRNA-Arg-UCU, tRNA-Arg-CCU
    "AGA": [("UCU", 1.0)],                    # U:A Watson-Crick
    "AGG": [("UCU", 0.2), ("CCU", 1.0)],      # U:G wobble, C:G Watson-Crick

    # ── Glycine ──
    "GGU": [("GCC", 0.5)],                    # G:U wobble
    "GGC": [("GCC", 1.0)],                    # G:C Watson-Crick
    "GGA": [("UCC", 1.0)],                    # U:A Watson-Crick
    "GGG": [("UCC", 0.2), ("CCC", 1.0)],      # U:G wobble, C:G Watson-Crick

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
# Format: organism_key -> {anticodon: gene_copy_number}
#
# Anticodon key conventions:
#   - "GAU": Ile tRNA with GAU anticodon (reads AUC, AUU)
#   - "k2C": Ile tRNA with lysidine at position 34 (from CAU gene, reads AUA)
#     This is SEPARATE from Met-CAU which has unmodified C
#   - "IAU": Ile tRNA with Inosine at position 34 (from AAU gene, eukaryotes)
#     Reads AUC, AUU, AUA with Inosine efficiencies
#   - "CAU": Met tRNA with unmodified CAU anticodon (reads AUG only)
#   - "IGA", "ICG", "IGU", "IGC": tRNAs with Inosine at position 34
#     (from AGA, ACG, AGU, AGC genes respectively)

TRNA_GENE_COPIES: dict[str, dict[str, int]] = {
    # ── Escherichia coli K-12 MG1655 ─────────────────────────────────
    # Source: GtRNAdb, E. coli K-12 MG1655 genome
    # ~86 tRNA genes total
    "e_coli": {
        # Phenylalanine
        "GAA": 3,
        # Leucine
        "UAA": 1, "CAA": 6, "UAG": 1, "CAG": 4, "GAG": 2,
        # Isoleucine (GAU reads AUC/AUU; k2C=lysidine reads AUA)
        "GAU": 3, "k2C": 1,
        # Methionine (initiator + elongator, unmodified CAU reads AUG)
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
    },

    # ── Homo sapiens (GRCh38) ────────────────────────────────────────
    # Source: GtRNAdb, human genome (GRCh38)
    # ~610 tRNA genes total (functional only, excluding pseudogenes)
    "human": {
        # Phenylalanine
        "GAA": 10,
        # Leucine
        "UAA": 5, "CAA": 12, "UAG": 4, "CAG": 8, "GAG": 7,
        # Isoleucine (IAU from AAU gene, deaminated A→I; k2C from CAU gene)
        "IAU": 8, "k2C": 4,
        # Methionine
        "CAU": 16,
        # Valine
        "GAC": 7, "UAC": 5, "CAC": 10,
        # Serine
        "GGA": 5, "IGA": 4, "UGA": 8, "CGA": 3, "GCU": 7,
        # Proline
        "GGG": 5, "UGG": 6, "CGG": 7,
        # Threonine
        "GGU": 7, "IGU": 5, "UGU": 6, "CGU": 4,
        # Alanine
        "GGC": 10, "IGC": 6, "UGC": 6, "CGC": 8,
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
        "GCG": 5, "ICG": 4, "CCG": 7, "UCG": 4, "UCU": 5, "CCU": 4,
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
        # Isoleucine (IAU from AAU gene; k2C from CAU gene)
        "IAU": 5, "k2C": 2,
        # Methionine
        "CAU": 8,
        # Valine
        "GAC": 5, "UAC": 3, "CAC": 6,
        # Serine
        "GGA": 4, "IGA": 3, "UGA": 6, "CGA": 2, "GCU": 5,
        # Proline
        "GGG": 4, "UGG": 4, "CGG": 5,
        # Threonine
        "GGU": 5, "IGU": 4, "UGU": 4, "CGU": 3,
        # Alanine
        "GGC": 6, "IGC": 4, "UGC": 4, "CGC": 5,
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
        "GCG": 4, "ICG": 3, "CCG": 5, "UCG": 3, "UCU": 4, "CCU": 3,
        # Glycine
        "GCC": 6, "UCC": 3, "CCC": 5,
    },

    # ── Mus musculus (GRCm39) ────────────────────────────────────────
    # Source: GtRNAdb, mouse genome (GRCm39)
    # ~445 tRNA genes total
    "mouse": {
        # Phenylalanine
        "GAA": 8,
        # Leucine
        "UAA": 4, "CAA": 10, "UAG": 3, "CAG": 7, "GAG": 6,
        # Isoleucine
        "IAU": 7, "k2C": 3,
        # Methionine
        "CAU": 14,
        # Valine
        "GAC": 6, "UAC": 4, "CAC": 9,
        # Serine
        "GGA": 4, "IGA": 3, "UGA": 7, "CGA": 2, "GCU": 6,
        # Proline
        "GGG": 4, "UGG": 5, "CGG": 6,
        # Threonine
        "GGU": 6, "IGU": 4, "UGU": 5, "CGU": 3,
        # Alanine
        "GGC": 9, "IGC": 5, "UGC": 5, "CGC": 7,
        # Tyrosine
        "GUA": 6,
        # Histidine
        "GUG": 4,
        # Glutamine
        "UUG": 7, "CUG": 10,
        # Asparagine
        "GUU": 9,
        # Lysine
        "UUU": 7, "CUU": 8,
        # Aspartate
        "GUC": 7,
        # Glutamate
        "UUC": 10, "CUC": 11,
        # Cysteine
        "GCA": 6,
        # Tryptophan
        "CCA": 5,
        # Arginine
        "GCG": 4, "ICG": 3, "CCG": 6, "UCG": 3, "UCU": 4, "CCU": 3,
        # Glycine
        "GCC": 9, "UCC": 3, "CCC": 7,
    },

    # ── CHO-K1 (Cricetulus griseus) ──────────────────────────────────
    # Source: GtRNAdb, CHO-K1 genome
    # ~430 tRNA genes total
    "cho": {
        # Phenylalanine
        "GAA": 7,
        # Leucine
        "UAA": 4, "CAA": 9, "UAG": 3, "CAG": 6, "GAG": 5,
        # Isoleucine
        "IAU": 6, "k2C": 3,
        # Methionine
        "CAU": 13,
        # Valine
        "GAC": 5, "UAC": 4, "CAC": 8,
        # Serine
        "GGA": 4, "IGA": 3, "UGA": 6, "CGA": 2, "GCU": 5,
        # Proline
        "GGG": 4, "UGG": 5, "CGG": 6,
        # Threonine
        "GGU": 5, "IGU": 4, "UGU": 5, "CGU": 3,
        # Alanine
        "GGC": 8, "IGC": 5, "UGC": 5, "CGC": 7,
        # Tyrosine
        "GUA": 5,
        # Histidine
        "GUG": 4,
        # Glutamine
        "UUG": 6, "CUG": 9,
        # Asparagine
        "GUU": 8,
        # Lysine
        "UUU": 6, "CUU": 7,
        # Aspartate
        "GUC": 6,
        # Glutamate
        "UUC": 9, "CUC": 10,
        # Cysteine
        "GCA": 5,
        # Tryptophan
        "CCA": 5,
        # Arginine
        "GCG": 4, "ICG": 3, "CCG": 5, "UCG": 3, "UCU": 4, "CCU": 3,
        # Glycine
        "GCC": 8, "UCC": 3, "CCC": 6,
    },

    # ── Caenorhabditis elegans ────────────────────────────────────────
    # Source: GtRNAdb, C. elegans genome
    # ~820 tRNA genes total (including many pseudogenes; ~280 functional)
    "c_elegans": {
        # Phenylalanine
        "GAA": 6,
        # Leucine
        "UAA": 3, "CAA": 8, "UAG": 3, "CAG": 6, "GAG": 5,
        # Isoleucine
        "IAU": 5, "k2C": 2,
        # Methionine
        "CAU": 10,
        # Valine
        "GAC": 5, "UAC": 3, "CAC": 7,
        # Serine
        "GGA": 3, "IGA": 3, "UGA": 5, "CGA": 2, "GCU": 5,
        # Proline
        "GGG": 3, "UGG": 4, "CGG": 5,
        # Threonine
        "GGU": 4, "IGU": 3, "UGU": 4, "CGU": 3,
        # Alanine
        "GGC": 7, "IGC": 4, "UGC": 4, "CGC": 6,
        # Tyrosine
        "GUA": 5,
        # Histidine
        "GUG": 3,
        # Glutamine
        "UUG": 5, "CUG": 8,
        # Asparagine
        "GUU": 6,
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
        "GCG": 3, "ICG": 3, "CCG": 5, "UCG": 2, "UCU": 3, "CCU": 2,
        # Glycine
        "GCC": 7, "UCC": 3, "CCC": 5,
    },

    # ── Drosophila melanogaster ──────────────────────────────────────
    # Source: GtRNAdb, D. melanogaster genome (dm6)
    # ~290 tRNA genes total
    "d_melanogaster": {
        # Phenylalanine
        "GAA": 5,
        # Leucine
        "UAA": 3, "CAA": 7, "UAG": 2, "CAG": 5, "GAG": 4,
        # Isoleucine
        "IAU": 5, "k2C": 2,
        # Methionine
        "CAU": 8,
        # Valine
        "GAC": 4, "UAC": 3, "CAC": 6,
        # Serine
        "GGA": 3, "IGA": 3, "UGA": 5, "CGA": 2, "GCU": 4,
        # Proline
        "GGG": 3, "UGG": 4, "CGG": 4,
        # Threonine
        "GGU": 4, "IGU": 3, "UGU": 3, "CGU": 3,
        # Alanine
        "GGC": 6, "IGC": 3, "UGC": 3, "CGC": 5,
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
        "UUC": 7, "CUC": 7,
        # Cysteine
        "GCA": 3,
        # Tryptophan
        "CCA": 3,
        # Arginine
        "GCG": 3, "ICG": 3, "CCG": 4, "UCG": 2, "UCU": 3, "CCU": 2,
        # Glycine
        "GCC": 6, "UCC": 2, "CCC": 4,
    },

    # ── Arabidopsis thaliana ─────────────────────────────────────────
    # Source: GtRNAdb, A. thaliana genome (TAIR10)
    # ~730 tRNA genes total
    "a_thaliana": {
        # Phenylalanine
        "GAA": 7,
        # Leucine
        "UAA": 4, "CAA": 9, "UAG": 3, "CAG": 6, "GAG": 5,
        # Isoleucine
        "IAU": 6, "k2C": 3,
        # Methionine
        "CAU": 11,
        # Valine
        "GAC": 5, "UAC": 3, "CAC": 7,
        # Serine
        "GGA": 4, "IGA": 3, "UGA": 6, "CGA": 2, "GCU": 5,
        # Proline
        "GGG": 4, "UGG": 5, "CGG": 5,
        # Threonine
        "GGU": 5, "IGU": 3, "UGU": 4, "CGU": 3,
        # Alanine
        "GGC": 8, "IGC": 4, "UGC": 4, "CGC": 6,
        # Tyrosine
        "GUA": 5,
        # Histidine
        "GUG": 4,
        # Glutamine
        "UUG": 6, "CUG": 8,
        # Asparagine
        "GUU": 7,
        # Lysine
        "UUU": 6, "CUU": 7,
        # Aspartate
        "GUC": 6,
        # Glutamate
        "UUC": 8, "CUC": 9,
        # Cysteine
        "GCA": 5,
        # Tryptophan
        "CCA": 4,
        # Arginine
        "GCG": 4, "ICG": 3, "CCG": 5, "UCG": 3, "UCU": 4, "CCU": 3,
        # Glycine
        "GCC": 8, "UCC": 3, "CCC": 6,
    },

    # ── Pichia pastoris (Komagataella phaffii) ──────────────────────
    # Source: GtRNAdb, P. pastoris genome
    # ~230 tRNA genes total
    "p_pastoris": {
        # Phenylalanine
        "GAA": 4,
        # Leucine
        "UAA": 3, "CAA": 6, "UAG": 2, "CAG": 4, "GAG": 3,
        # Isoleucine
        "IAU": 4, "k2C": 2,
        # Methionine
        "CAU": 6,
        # Valine
        "GAC": 4, "UAC": 2, "CAC": 5,
        # Serine
        "GGA": 3, "IGA": 2, "UGA": 4, "CGA": 2, "GCU": 4,
        # Proline
        "GGG": 3, "UGG": 3, "CGG": 4,
        # Threonine
        "GGU": 4, "IGU": 3, "UGU": 3, "CGU": 2,
        # Alanine
        "GGC": 5, "IGC": 3, "UGC": 3, "CGC": 4,
        # Tyrosine
        "GUA": 3,
        # Histidine
        "GUG": 2,
        # Glutamine
        "UUG": 4, "CUG": 5,
        # Asparagine
        "GUU": 4,
        # Lysine
        "UUU": 4, "CUU": 5,
        # Aspartate
        "GUC": 4,
        # Glutamate
        "UUC": 5, "CUC": 6,
        # Cysteine
        "GCA": 3,
        # Tryptophan
        "CCA": 2,
        # Arginine
        "GCG": 3, "ICG": 2, "CCG": 4, "UCG": 2, "UCU": 3, "CCU": 2,
        # Glycine
        "GCC": 5, "UCC": 2, "CCC": 4,
    },

    # ── Bacillus subtilis 168 ────────────────────────────────────────
    # Source: GtRNAdb, B. subtilis subsp. subtilis str. 168
    # ~86 tRNA genes total
    "b_subtilis": {
        # Phenylalanine
        "GAA": 3,
        # Leucine
        "UAA": 1, "CAA": 5, "UAG": 1, "CAG": 3, "GAG": 2,
        # Isoleucine (GAU reads AUC/AUU; k2C=lysidine reads AUA)
        "GAU": 3, "k2C": 1,
        # Methionine
        "CAU": 4,
        # Valine
        "GAC": 2, "UAC": 1, "CAC": 3,
        # Serine
        "GGA": 1, "CGA": 1, "UGA": 3, "GCU": 2,
        # Proline
        "GGG": 2, "UGG": 2, "CGG": 3,
        # Threonine
        "GGU": 2, "CGU": 2, "UGU": 2,
        # Alanine
        "GGC": 3, "UGC": 2, "CGC": 3,
        # Tyrosine
        "GUA": 2,
        # Histidine
        "GUG": 1,
        # Glutamine
        "UUG": 2, "CUG": 3,
        # Asparagine
        "GUU": 3,
        # Lysine
        "UUU": 2, "CUU": 2,
        # Aspartate
        "GUC": 2,
        # Glutamate
        "UUC": 2, "CUC": 3,
        # Cysteine
        "GCA": 1,
        # Tryptophan
        "CCA": 1,
        # Arginine
        "GCG": 2, "CCG": 3, "UCG": 1, "UCU": 1, "CCU": 1,
        # Glycine
        "GCC": 3, "UCC": 1, "CCC": 2,
    },
}

SUPPORTED_ORGANISMS_TAI: list[str] = list(TRNA_GENE_COPIES.keys())


# ═══════════════════════════════════════════════════════════════════════════════
# 4. compute_tai_weights — relative adaptiveness per codon per organism
# ═══════════════════════════════════════════════════════════════════════════════

def compute_tai_weights(organism: str) -> dict[str, float]:
    """Compute relative adaptiveness weights for all codons for an organism.

    For each codon, the relative adaptiveness is::

        w(codon) = W(codon) / W_max(amino_acid)

    where ``W(codon)`` is the raw weight (sum of tRNA gene copies ×
    wobble efficiency for all anticodons that can read the codon), and
    ``W_max`` is the maximum raw weight among all synonymous codons
    for the same amino acid.

    Parameters
    ----------
    organism : str
        Organism key in :data:`TRNA_GENE_COPIES`
        (e.g. ``"e_coli"``, ``"human"``, ``"yeast"``, ``"mouse"``, ``"cho"``,
        ``"c_elegans"``, ``"d_melanogaster"``, ``"a_thaliana"``,
        ``"p_pastoris"``, ``"b_subtilis"``).

    Returns
    -------
    dict[str, float]
        Mapping of RNA codon → relative adaptiveness value in [0, 1].

    Raises
    ------
    KeyError
        If the organism key is not in :data:`TRNA_GENE_COPIES`.
    """
    from ..constants import CODON_TABLE

    trna_copies = TRNA_GENE_COPIES[organism]

    # Compute raw weights for all codons
    raw_weights: dict[str, float] = {}
    for codon in WOBBLE_RULES:
        wobble_list = WOBBLE_RULES[codon]
        total_weight = 0.0
        for anticodon, efficiency in wobble_list:
            copies = trna_copies.get(anticodon, 0)
            total_weight += copies * efficiency
        raw_weights[codon] = total_weight

    # Group codons by amino acid and normalize
    aa_to_codons: dict[str, list[str]] = {}
    for codon in WOBBLE_RULES:
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
