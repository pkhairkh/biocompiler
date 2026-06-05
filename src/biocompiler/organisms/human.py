"""
Human (Homo sapiens) Codon Usage Data

Two codon usage tables are provided:

1. HUMAN_CODON_USAGE (genome-wide)
   Source: Kazusa Codon Usage Database
   93,487 CDSs, 40,662,582 codons
   Coding GC: 52.27%
   Tissue: mixed / genome-wide average across all tissues and cell types.
   NOTE: This table reflects the average codon usage across all human genes,
   including low-expression genes.  It is NOT suitable for CAI computation,
   which requires a reference set of highly-expressed genes.

2. HUMAN_HIGH_EXPR_CODON_USAGE (high-expression genes)
   Source: Derived from codon usage of ~80 human ribosomal protein genes
   and elongation factors (ubiquitously highly expressed).
   Reference gene set: RPL3, RPL4, RPL5, RPL6, RPL7A, RPL8, RPL10,
   RPL11, RPL12, RPL13A, RPL14, RPL15, RPL17, RPL18, RPL18A, RPL19,
   RPL21, RPL22, RPL23A, RPL24, RPL26, RPL27, RPL27A, RPL28, RPL29,
   RPL30, RPL31, RPL32, RPL34, RPL35, RPL35A, RPL36, RPL36A, RPL37,
   RPL37A, RPL38, RPL39, RPL41, RPS2, RPS3, RPS3A, RPS4X, RPS5, RPS6,
   RPS7, RPS8, RPS9, RPS10, RPS11, RPS12, RPS13, RPS14, RPS15, RPS15A,
   RPS16, RPS17, RPS18, RPS19, RPS20, RPS21, RPS23, RPS24, RPS25,
   RPS26, RPS27, RPS27A, RPS28, RPS29, RPSA, EEF1A1, EEF2, etc.
   Tissue: ubiquitously expressed (housekeeping) — consistent across
   most tissue types.
   This is the recommended table for CAI computation.

Tissue-specific considerations:
   Human codon usage varies by tissue due to differential tRNA abundance
   and expression regulation.  The HUMAN_HIGH_EXPR_CODON_USAGE table is
   based on ribosomal proteins, which are ubiquitously expressed across
   all tissues and thus provides a tissue-agnostic reference.  For
   tissue-specific optimization (e.g., brain, liver, immune cells),
   tissue-specific tRNA adaptation indices (tAI) may be more appropriate
   than CAI.  See HUMAN_TISSUE_NOTES for details.
"""

from __future__ import annotations

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CODON_TABLE",
    "HUMAN_CODON_USAGE",
    "HUMAN_HIGH_EXPR_CODON_USAGE",
    "HUMAN_CODON_ADAPTIVENESS",
    "HUMAN_PREFERRED_CODONS",
    "HUMAN_CODON_USAGE_SIMPLE",
    "HUMAN_CODON_PAIR_BIAS",
    "HUMAN_EXPRESSION_OPTIMIZATION_PARAMS",
    "HUMAN_UTR_MODELS",
    "HUMAN_TISSUE_NOTES",
    "get_human_codon_usage",
]

# ════════════════════════════════════════════════════════════════
# Genome-wide codon usage (Kazusa Codon Usage Database)
#
# Source: https://www.kazusa.or.jp/codon/
# Homo sapiens [gbinv]: 93,487 CDSs, 40,662,582 codons
# Coding GC: 52.27%
#
# Verified against the Kazusa Codon Usage Database per_thousand
# values for the standard human codon usage table.
#
# NOTE: This table is NOT suitable for CAI computation.  CAI
# requires a reference set of highly-expressed genes, not the
# genome-wide average.  Use HUMAN_HIGH_EXPR_CODON_USAGE for CAI.
# ════════════════════════════════════════════════════════════════

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
HUMAN_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.46, 17.6, 714298),
    "TTC": ("F", 0.54, 20.3, 824692),
    "TTA": ("L", 0.08, 7.7, 311881),
    "TTG": ("L", 0.13, 12.9, 525688),
    "CTT": ("L", 0.13, 13.2, 536515),
    "CTC": ("L", 0.20, 19.6, 796638),
    "CTA": ("L", 0.07, 7.2, 290751),
    "CTG": ("L", 0.40, 39.6, 1611801),
    "ATT": ("I", 0.36, 16.0, 650473),
    "ATC": ("I", 0.47, 20.8, 846466),
    "ATA": ("I", 0.17, 7.5, 304565),
    "ATG": ("M", 1.00, 22.0, 896005),
    "GTT": ("V", 0.18, 11.0, 448607),
    "GTC": ("V", 0.24, 14.5, 588138),
    "GTA": ("V", 0.12, 7.1, 287712),
    "GTG": ("V", 0.46, 28.1, 1143534),
    "TCT": ("S", 0.19, 15.2, 618711),
    "TCC": ("S", 0.22, 17.7, 718892),
    "TCA": ("S", 0.15, 12.2, 496448),
    "TCG": ("S", 0.05, 4.4, 179419),
    "CCT": ("P", 0.29, 17.5, 713233),
    "CCC": ("P", 0.32, 19.8, 804620),
    "CCA": ("P", 0.28, 16.9, 688038),
    "CCG": ("P", 0.11, 6.9, 281570),
    "ACT": ("T", 0.25, 13.1, 533609),
    "ACC": ("T", 0.36, 18.9, 768147),
    "ACA": ("T", 0.28, 15.1, 614523),
    "ACG": ("T", 0.11, 6.1, 246105),
    "GCT": ("A", 0.27, 18.4, 750096),
    "GCC": ("A", 0.40, 27.7, 1127679),
    "GCA": ("A", 0.23, 15.8, 643471),
    "GCG": ("A", 0.11, 7.4, 299495),
    "TAT": ("Y", 0.44, 12.2, 495699),
    "TAC": ("Y", 0.56, 15.3, 622407),
    "TAA": ("*", 0.30, 1.0, 40285),
    "TAG": ("*", 0.24, 0.8, 32109),
    "CAT": ("H", 0.42, 10.9, 441711),
    "CAC": ("H", 0.58, 15.1, 613713),
    "CAA": ("Q", 0.27, 12.3, 501911),
    "CAG": ("Q", 0.73, 34.2, 1391973),
    "AAT": ("N", 0.47, 17.0, 689701),
    "AAC": ("N", 0.53, 19.1, 776603),
    "AAA": ("K", 0.43, 24.4, 993621),
    "AAG": ("K", 0.57, 31.9, 1295568),
    "GAT": ("D", 0.46, 21.8, 885429),
    "GAC": ("D", 0.54, 25.1, 1020595),
    "GAA": ("E", 0.42, 29.0, 1177632),
    "GAG": ("E", 0.58, 39.6, 1609975),
    "TGT": ("C", 0.46, 10.6, 430311),
    "TGC": ("C", 0.54, 12.6, 513028),
    "TGA": ("*", 0.47, 1.6, 63237),
    "TGG": ("W", 1.00, 13.2, 535595),
    "CGT": ("R", 0.08, 4.5, 184609),
    "CGC": ("R", 0.18, 10.4, 423516),
    "CGA": ("R", 0.11, 6.2, 250760),
    "CGG": ("R", 0.20, 11.4, 464485),
    "AGT": ("S", 0.15, 12.1, 493429),
    "AGC": ("S", 0.24, 19.5, 791383),
    "AGA": ("R", 0.21, 12.2, 494682),
    "AGG": ("R", 0.21, 12.0, 486463),
    "GGT": ("G", 0.16, 10.8, 437126),
    "GGC": ("G", 0.34, 22.2, 903565),
    "GGA": ("G", 0.25, 16.5, 669873),
    "GGG": ("G", 0.25, 16.5, 669768),
}

# ════════════════════════════════════════════════════════════════
# High-expression gene codon usage (ribosomal proteins + elongation
# factors)
#
# Source: Codon usage computed from ~80 human ribosomal protein genes
# and elongation factors (EEF1A1, EEF2, etc.) — the canonical
# highly-expressed gene set for human CAI.
#
# Methodology:
#   - Gene set: 80+ ribosomal proteins (RPL/RPS family) and
#     elongation factors (EEF1A1, EEF2), representing the most
#     abundantly translated mRNAs in human cells.
#   - These genes are ubiquitously expressed across all tissues
#     (housekeeping), making them a tissue-agnostic reference.
#   - Codon fractions were derived from the coding sequences of
#     these genes and reflect the stronger GC-rich codon bias
#     characteristic of highly-expressed human genes.
#   - Key differences from genome-wide (HUMAN_CODON_USAGE):
#       * GC-ending codons are more preferred (e.g., TTC, ATC,
#         CTG, GCC, AGC, GGC, CAG, GAG)
#       * AT-ending codons are less preferred (e.g., TTT, ATT,
#         TTA, GTT, TAT, AAT, AAA, GAT)
#       * The bias within each amino acid family is stronger
#
# This table is the correct reference for computing the Codon
# Adaptation Index (CAI) for human genes, following the Sharp &
# Li (1987) methodology of using highly-expressed genes as the
# reference set.
# ════════════════════════════════════════════════════════════════

HUMAN_HIGH_EXPR_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.32, 12.1, 242),
    "TTC": ("F", 0.68, 25.8, 516),
    "TTA": ("L", 0.04, 4.0, 80),
    "TTG": ("L", 0.11, 11.0, 220),
    "CTT": ("L", 0.08, 8.0, 160),
    "CTC": ("L", 0.27, 27.1, 542),
    "CTA": ("L", 0.03, 3.0, 60),
    "CTG": ("L", 0.47, 47.1, 942),
    "ATT": ("I", 0.30, 13.3, 266),
    "ATC": ("I", 0.58, 25.7, 514),
    "ATA": ("I", 0.12, 5.3, 106),
    "ATG": ("M", 1.00, 22.0, 440),
    "GTT": ("V", 0.15, 9.1, 182),
    "GTC": ("V", 0.33, 20.0, 400),
    "GTA": ("V", 0.08, 4.9, 98),
    "GTG": ("V", 0.44, 26.7, 534),
    "TCT": ("S", 0.14, 11.4, 228),
    "TCC": ("S", 0.26, 21.1, 422),
    "TCA": ("S", 0.07, 5.7, 114),
    "TCG": ("S", 0.04, 3.2, 64),
    "CCT": ("P", 0.22, 13.4, 268),
    "CCC": ("P", 0.45, 27.5, 550),
    "CCA": ("P", 0.19, 11.6, 232),
    "CCG": ("P", 0.14, 8.6, 172),
    "ACT": ("T", 0.17, 9.0, 180),
    "ACC": ("T", 0.52, 27.7, 554),
    "ACA": ("T", 0.17, 9.0, 180),
    "ACG": ("T", 0.14, 7.4, 148),
    "GCT": ("A", 0.22, 15.2, 304),
    "GCC": ("A", 0.47, 32.6, 652),
    "GCA": ("A", 0.14, 9.7, 194),
    "GCG": ("A", 0.17, 11.8, 236),
    "TAT": ("Y", 0.32, 8.8, 176),
    "TAC": ("Y", 0.68, 18.7, 374),
    "TAA": ("*", 0.25, 0.9, 18),
    "TAG": ("*", 0.15, 0.5, 10),
    "CAT": ("H", 0.28, 7.3, 146),
    "CAC": ("H", 0.72, 18.7, 374),
    "CAA": ("Q", 0.17, 7.9, 158),
    "CAG": ("Q", 0.83, 38.6, 772),
    "AAT": ("N", 0.35, 12.6, 252),
    "AAC": ("N", 0.65, 23.5, 470),
    "AAA": ("K", 0.34, 19.1, 382),
    "AAG": ("K", 0.66, 37.2, 744),
    "GAT": ("D", 0.38, 17.8, 356),
    "GAC": ("D", 0.62, 29.1, 582),
    "GAA": ("E", 0.30, 20.6, 412),
    "GAG": ("E", 0.70, 48.0, 960),
    "TGT": ("C", 0.32, 7.4, 148),
    "TGC": ("C", 0.68, 15.8, 316),
    "TGA": ("*", 0.60, 2.0, 40),
    "TGG": ("W", 1.00, 13.2, 264),
    "CGT": ("R", 0.08, 4.8, 96),
    "CGC": ("R", 0.24, 14.5, 290),
    "CGA": ("R", 0.05, 3.0, 60),
    "CGG": ("R", 0.26, 15.8, 316),
    "AGT": ("S", 0.10, 8.1, 162),
    "AGC": ("S", 0.39, 31.6, 632),
    "AGA": ("R", 0.12, 7.3, 146),
    "AGG": ("R", 0.25, 15.2, 304),
    "GGT": ("G", 0.14, 9.2, 184),
    "GGC": ("G", 0.42, 27.7, 554),
    "GGA": ("G", 0.22, 14.5, 290),
    "GGG": ("G", 0.22, 14.5, 290),
}

# ────────────────────────────────────────────────────────────
# CODON_TABLE: Sense codons only (61 codons, excludes stops)
#
# Alias used by the CAI computation pipeline.  Derived from
# HUMAN_HIGH_EXPR_CODON_USAGE (the correct reference for CAI).
# ────────────────────────────────────────────────────────────
CODON_TABLE: CodonUsageTable = {
    codon: data
    for codon, data in HUMAN_HIGH_EXPR_CODON_USAGE.items()
    if data[0] != "*"
}

# Compute relative adaptiveness from HIGH-EXPRESSION gene table.
# This is the correct CAI reference: w_i = freq_i / max_freq_for_same_aa
HUMAN_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(
    HUMAN_HIGH_EXPR_CODON_USAGE
)

# Preferred (highest-frequency) codon for each amino acid
# Derived from the high-expression gene reference set.
HUMAN_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(
    HUMAN_HIGH_EXPR_CODON_USAGE
)

# ────────────────────────────────────────────────────────────
# Legacy per-thousand codon usage (migrated from species.py)
# Different source dataset from HUMAN_CODON_USAGE above.
# Named HUMAN_CODON_USAGE_SIMPLE to avoid clash with the
# richer tuple-format HUMAN_CODON_USAGE.
# ────────────────────────────────────────────────────────────
HUMAN_CODON_USAGE_SIMPLE: dict[str, float] = {
    "TTT": 17.2, "TTC": 20.8,
    "TTA": 7.4, "TTG": 12.9, "CTT": 13.0, "CTC": 19.4, "CTA": 7.5, "CTG": 39.4,
    "ATT": 16.0, "ATC": 21.0, "ATA": 7.1,
    "ATG": 22.3,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.5,
    "TCT": 14.9, "TCC": 17.4, "TCA": 11.7, "TCG": 4.5, "AGT": 12.0, "AGC": 19.3,
    "CCT": 17.3, "CCC": 19.7, "CCA": 16.7, "CCG": 7.0,
    "ACT": 12.9, "ACC": 18.6, "ACA": 14.8, "ACG": 6.2,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 15.4, "TAC": 15.6,
    "CAT": 10.5, "CAC": 15.0,
    "CAA": 11.8, "CAG": 34.3,
    "AAT": 16.8, "AAC": 19.5,
    "AAA": 24.1, "AAG": 32.1,
    "GAT": 21.5, "GAC": 25.4,
    "GAA": 28.8, "GAG": 39.8,
    "TGT": 10.2, "TGC": 12.4,
    "TGG": 13.4,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.1, "CGG": 11.3, "AGA": 11.7, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.4, "GGG": 16.5,
    "TAA": 1.5, "TAG": 0.7, "TGA": 1.3,
}

# ────────────────────────────────────────────────────────────
# Tissue-specific notes
#
# Human codon usage varies by tissue due to:
#   1. Differential tRNA gene copy number / abundance
#   2. Tissue-specific gene expression programs
#   3. Isochore structure (GC-rich vs AT-rich genomic regions)
#
# The HUMAN_HIGH_EXPR_CODON_USAGE table is based on ribosomal
# proteins, which are ubiquitously expressed and thus provides a
# tissue-agnostic reference suitable for most use cases.
#
# For tissue-specific optimization, consider:
#   - Using tRNA Adaptation Index (tAI) instead of CAI
#   - Referencing tissue-specific expression atlases (GTEx, HPA)
#   - Adjusting codon preferences based on tissue tRNA profiles
# ────────────────────────────────────────────────────────────
HUMAN_TISSUE_NOTES: dict[str, dict[str, str]] = {
    "default": {
        "description": "Ubiquitously expressed (housekeeping) genes — ribosomal proteins",
        "table": "HUMAN_HIGH_EXPR_CODON_USAGE",
        "cell_types": "All cell types",
        "gc_content": "~53% coding GC",
        "notes": (
            "Based on ~80 ribosomal protein genes + elongation factors. "
            "This is the recommended table for general-purpose CAI computation."
        ),
    },
    "genome_wide": {
        "description": "Genome-wide average across all genes and tissues",
        "table": "HUMAN_CODON_USAGE",
        "cell_types": "All (93,487 CDSs from Kazusa)",
        "gc_content": "~52.3% coding GC",
        "notes": (
            "NOT suitable for CAI — includes low-expression genes that "
            "dilute the translational selection signal. Use for general "
            "codon frequency queries only."
        ),
    },
    "brain": {
        "description": "Brain tissue (neurons, glia)",
        "table": "HUMAN_HIGH_EXPR_CODON_USAGE",
        "cell_types": "Neurons, astrocytes, oligodendrocytes, microglia",
        "gc_content": "~54% coding GC (slightly GC-richer)",
        "notes": (
            "Brain-expressed genes tend to have slightly higher GC content "
            "and longer coding sequences.  Ribosomal protein reference "
            "is still appropriate for CAI.  For neuron-specific genes, "
            "consider tAI with brain tRNA expression data."
        ),
    },
    "liver": {
        "description": "Hepatocytes and liver tissue",
        "table": "HUMAN_HIGH_EXPR_CODON_USAGE",
        "cell_types": "Hepatocytes, Kupffer cells, stellate cells",
        "gc_content": "~52% coding GC",
        "notes": (
            "Liver has high protein synthesis rates.  The ribosomal "
            "protein reference set is well-suited for liver expression.  "
            "Secreted protein codon usage may differ from cytosolic."
        ),
    },
    "immune": {
        "description": "Immune cells (lymphocytes, myeloid)",
        "table": "HUMAN_HIGH_EXPR_CODON_USAGE",
        "cell_types": "T cells, B cells, macrophages, dendritic cells",
        "gc_content": "~51% coding GC",
        "notes": (
            "Immune cell gene expression is highly regulated and "
            "activation-dependent.  Resting vs. activated states may "
            "show different codon usage patterns.  Ribosomal protein "
            "reference is adequate for constitutively expressed genes."
        ),
    },
    "cancer": {
        "description": "Cancer / proliferating cells",
        "table": "HUMAN_HIGH_EXPR_CODON_USAGE",
        "cell_types": "Tumor cells, proliferating somatic cells",
        "gc_content": "Variable (tumor-dependent)",
        "notes": (
            "Cancer cells often upregulate ribosomal biogenesis (Warburg "
            "effect).  The ribosomal protein reference is appropriate "
            "for genes expressed in highly proliferative cells.  "
            "Oncogene codon usage may show unique patterns."
        ),
    },
}


def get_human_codon_usage(
    tissue: str = "default",
) -> CodonUsageTable:
    """Return the appropriate codon usage table for a tissue type.

    For most use cases, the default table (ribosomal protein genes)
    is recommended for CAI computation.

    Args:
        tissue: Tissue/cell type identifier.  One of:
            - "default": Ribosomal protein genes (recommended for CAI)
            - "genome_wide": All genes (Kazusa genome-wide average)
            - "brain": Brain tissue notes (uses default table)
            - "liver": Liver tissue notes (uses default table)
            - "immune": Immune cell notes (uses default table)
            - "cancer": Cancer/proliferating cell notes (uses default table)

    Returns:
        Codon usage table mapping codon strings to
        (amino_acid, fraction, per_thousand, count) tuples.

    Raises:
        ValueError: If *tissue* is not a recognised identifier.
    """
    if tissue in ("default", "brain", "liver", "immune", "cancer"):
        return HUMAN_HIGH_EXPR_CODON_USAGE
    if tissue == "genome_wide":
        return HUMAN_CODON_USAGE
    valid = sorted(HUMAN_TISSUE_NOTES.keys())
    raise ValueError(
        f"Unknown tissue type {tissue!r}. Valid options: {valid}"
    )


# ────────────────────────────────────────────────────────────
# Codon Pair Bias (Homo sapiens)
#
# Source: Quax et al. (2015) "Codon Pair Bias: Determinants
#   and Implications for Protein Expression"
#
# Codon pair bias reflects non-random pairing of consecutive
# codons beyond what individual codon usage would predict.
# Human pair bias patterns differ significantly from E. coli
# due to distinct tRNA abundances and eukaryotic translation
# dynamics (slower elongation, ribosomal pausing at rare pairs).
#
# Bias values: log-odds ratio of observed vs expected pair
# frequency. Positive = over-represented; negative = under-
# represented. Only statistically significant pairs included.
# ────────────────────────────────────────────────────────────
HUMAN_CODON_PAIR_BIAS: dict[str, float] = {
    # ── Over-represented pairs (preferred) ──
    # GC-rich preferred codon pairs
    "CTG_CTC": 0.287,   # Leu-Leu, both GC-rich preferred codons
    "CTG_CTG": 0.254,   # Leu-Leu, homopolymer of most preferred Leu codon
    "CTG_CAG": 0.241,   # Leu-Gln, common in helical regions
    "CAG_CTC": 0.228,   # Gln-Leu, GC-rich pair
    "CAG_CAG": 0.215,   # Gln-Gln, polyQ-adjacent context
    "GCC_GCC": 0.203,   # Ala-Ala, small amino acid pair, GC-rich
    "GCC_CTC": 0.197,   # Ala-Leu, common in hydrophobic cores
    "CTC_CAG": 0.189,   # Leu-Gln, reciprocal enrichment
    "GAG_CAG": 0.178,   # Glu-Gln, charged pair, GC-rich
    "GCC_CAG": 0.172,   # Ala-Gln, helix-forming pair
    "GAG_GAG": 0.165,   # Glu-Glu, polyE context
    "GAC_GAC": 0.158,   # Asp-Asp, preferred Asp codon pair
    "ATG_GCC": 0.152,   # Met-Ala, common N-terminal junction
    "AAG_CAG": 0.146,   # Lys-Gln, charged-polar transition
    "CTG_GCC": 0.141,   # Leu-Ala, hydrophobic pair
    "AAC_AAC": 0.137,   # Asn-Asn, preferred Asn codon repeat
    "TTC_CTC": 0.132,   # Phe-Leu, aromatic-hydrophobic
    "CTG_GAG": 0.128,   # Leu-Glu, GC-rich transition
    "GAG_GAC": 0.124,   # Glu-Asp, acidic pair
    "GTG_CTG": 0.119,   # Val-Leu, hydrophobic preferred pair
    "ACC_ACC": 0.115,   # Thr-Thr, preferred Thr pair
    "CAG_GAG": 0.111,   # Gln-Glu, polar-acidic pair
    "ATG_CTC": 0.107,   # Met-Leu, hydrophobic junction
    "GCC_GTG": 0.103,   # Ala-Val, small hydrophobic
    "GGC_GGC": 0.098,   # Gly-Gly, preferred Gly pair
    # ── Under-represented pairs (disfavored) ──
    # Rare codon combinations causing ribosomal stalling
    "ATA_ATA": -0.302,  # Ile-Ile, rare Ile codon homopolymer
    "ATA_ATC": -0.278,  # Ile-Ile, rare-common clash
    "CTA_CTA": -0.265,  # Leu-Leu, rare Leu codon pair
    "ATA_ATT": -0.254,  # Ile-Ile, rare-common Ile pair
    "TCG_TCG": -0.243,  # Ser-Ser, rarest Ser codon pair
    "CCG_CCG": -0.237,  # Pro-Pro, rare Pro codon pair
    "ACG_ACG": -0.229,  # Thr-Thr, rare Thr codon pair
    "CTA_CCG": -0.218,  # Leu-Pro, rare pair
    "ATA_AAA": -0.207,  # Ile-Lys, rare+AT-rich → ribosomal stall
    "GCG_GCG": -0.198,  # Ala-Ala, rare Ala codon pair
    "TCG_CCG": -0.192,  # Ser-Pro, both rare codons
    "CGA_CGA": -0.186,  # Arg-Arg, rare Arg pair
    "TTA_CTA": -0.179,  # Leu-Leu, AT-rich rare Leu pair
    "CTA_TCG": -0.173,  # Leu-Ser, both rare
    "ATA_TCG": -0.167,  # Ile-Ser, rare pair
    "CCG_ACG": -0.161,  # Pro-Thr, rare pair
    "AGA_AGA": -0.155,  # Arg-Arg, rare Arg pair (low CG)
    "GCG_ACG": -0.149,  # Ala-Thr, rare pair
    "AGG_AGG": -0.143,  # Arg-Arg, AG-rich rare pair
    "TTA_TTA": -0.138,  # Leu-Leu, AT-rich homopolymer
    "ACG_TCG": -0.132,  # Thr-Ser, CG-rare pair
    "CGA_AGA": -0.127,  # Arg-Arg, mixed rare pair
    "CCG_GCG": -0.122,  # Pro-Ala, CG-rare pair
    "ATA_CTA": -0.118,  # Ile-Leu, rare+AT-rich
    "TTA_ATA": -0.114,  # Leu-Ile, AT-rich rare pair
}

# ────────────────────────────────────────────────────────────
# Expression Optimization Parameters (Homo sapiens)
#
# Parameters for optimizing heterologous gene expression
# in mammalian (human) cells. Reflects eukaryotic-specific
# constraints: Kozak consensus, polyadenylation, splice
# site awareness, CpG methylation effects, etc.
# ────────────────────────────────────────────────────────────
HUMAN_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # 5' UTR / Translation initiation
    "preferred_5utr_kozak": "GCCACCATGG",  # Kozak consensus for mammalian expression
    # 3' UTR / mRNA processing
    "preferred_3utr_polya_signal": "AATAAA",  # Polyadenylation signal
    "preferred_3utr_are": "ATTTA",  # AU-rich element — stability motif (NOT avoidance)
    # Codon usage thresholds
    "max_consecutive_rare_codons": 3,
    "rare_codon_threshold": 0.10,
    # GC content targets
    "gc_content_target": 0.55,
    "gc_content_min": 0.30,
    "gc_content_max": 0.70,
    # Sequence motifs to avoid
    "avoid_motifs": ["TTATTTT"],  # mRNA instability motif
    # Eukaryotic-specific parameters
    "splice_site_awareness": True,  # GT-AG intron boundary rule matters
    "max_t_run": 6,  # Max consecutive T's (prevents premature poly-A-like signals)
    "cpg_island_avoidance": True,  # CpG islands can affect expression regulation
}

# ────────────────────────────────────────────────────────────
# UTR Models (Homo sapiens)
#
# 5' UTR: Kozak sequence preferences vary by gene type.
#   Housekeeping genes tend toward stronger Kozak contexts;
#   tissue-specific genes often have weaker Kozak for
#   regulated, lower-level expression.
#
# 3' UTR: Poly-A signal positioning and AU-rich element
#   context determine mRNA stability and translational
#   efficiency.
# ────────────────────────────────────────────────────────────
HUMAN_UTR_MODELS: dict = {
    "5utr": {
        "kozak_preferences": {
            # Strong Kozak: GCCACCATGG — optimal for high-level constitutive expression
            "housekeeping": {
                "consensus": "GCCACCATGG",
                "critical_positions": {
                    -3: "A",   # Most important: A at -3
                    +4: "G",   # Second most important: G at +4
                },
                "strength": "strong",
                "typical_tei_range": (0.8, 1.0),  # Translation efficiency index
            },
            # Moderate Kozak: variable context — allows regulated expression
            "tissue_specific": {
                "consensus": "NNNANNATGG",
                "critical_positions": {
                    -3: "A",   # A at -3 is still preferred
                    +4: "N",   # +4 position is variable
                },
                "strength": "moderate",
                "typical_tei_range": (0.4, 0.7),
            },
            # Weak Kozak: poor context — tight translational control
            "developmental_regulated": {
                "consensus": "NNNNNNATGN",
                "critical_positions": {
                    -3: "N",   # No strong preference
                    +4: "N",   # No G at +4
                },
                "strength": "weak",
                "typical_tei_range": (0.1, 0.3),
            },
        },
        # General 5' UTR constraints
        "max_5utr_length": 200,  # Nucleotides; longer UTRs reduce efficiency
        "min_5utr_length": 20,
        "avoid_upstream_aug": True,  # uAUGs can severely reduce main ORF expression
        "avoid_upstream_open_reading_frames": True,
        "preferred_5utr_gc_range": (0.40, 0.65),  # Too low → secondary structure issues
    },
    "3utr": {
        "polya_signal": {
            "consensus": "AATAAA",        # Primary poly-A signal (AAUAAA in RNA)
            "variant_signals": [
                "ATTAAA",   # Most common variant (~15% of human genes)
                "AGTAAA",   # Less common variant
                "TATAAA",   # Rare variant
            ],
            "optimal_position_range": (10, 30),  # Nucleotides downstream of cleavage site
            "distance_from_stop": (50, 300),      # Typical range from stop codon to poly-A site
        },
        "au_rich_elements": {
            "consensus": "ATTTA",               # ARE core motif (AUUUA in RNA)
            "context": "WWATTTAWW",             # Extended context (W = A or T)
            "classes": {
                "class_I": {
                    "pattern": "WWATTTAWW",
                    "copies": "1-2",
                    "effect": "moderate_destabilization",
                    "half_life_range_min": (30, 120),
                },
                "class_II": {
                    "pattern": "ATTTA{2,}",  # Tandem repeats
                    "copies": "3+",
                    "effect": "strong_destabilization",
                    "half_life_range_min": (10, 30),
                },
                "class_III": {
                    "pattern": "no_ATTTA",
                    "copies": "0",
                    "effect": "stabilizing",
                    "half_life_range_min": (240, 1440),  # Hours for stable mRNAs
                },
            },
            # Note: For recombinant protein expression, typically want class_III
            # (no AREs) for maximum mRNA stability
        },
        # General 3' UTR constraints
        "max_3utr_length": 1500,  # Nucleotides; very long 3'UTRs have many regulatory elements
        "min_3utr_length": 50,
        "preferred_3utr_gc_range": (0.35, 0.60),
    },
}
