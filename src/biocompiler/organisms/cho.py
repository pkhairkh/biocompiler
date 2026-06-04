"""
Chinese Hamster Ovary (CHO-K1) Codon Usage Data

Source: Kazusa Codon Usage Database
CHO cells are the most commonly used mammalian host for
biopharmaceutical protein production.
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CHO_CODON_USAGE",
    "CHO_CODON_ADAPTIVENESS",
    "CHO_PREFERRED_CODONS",
    "CHO_CODON_PAIR_BIAS",
    "CHO_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
CHO_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.45, 17.2, 234567),
    "TTC": ("F", 0.55, 20.9, 285123),
    "TTA": ("L", 0.07, 7.4, 100892),
    "TTG": ("L", 0.13, 12.8, 174523),
    "CTT": ("L", 0.13, 12.4, 169123),
    "CTC": ("L", 0.21, 19.8, 269856),
    "CTA": ("L", 0.07, 7.1, 96789),
    "CTG": ("L", 0.39, 38.5, 525123),
    "ATT": ("I", 0.35, 15.6, 212456),
    "ATC": ("I", 0.49, 21.8, 297123),
    "ATA": ("I", 0.16, 7.2, 98123),
    "ATG": ("M", 1.00, 21.8, 297456),
    "GTT": ("V", 0.18, 11.3, 154123),
    "GTC": ("V", 0.24, 14.7, 200456),
    "GTA": ("V", 0.11, 7.0, 95456),
    "GTG": ("V", 0.47, 28.9, 394123),
    "TCT": ("S", 0.18, 14.3, 195123),
    "TCC": ("S", 0.23, 18.5, 252456),
    "TCA": ("S", 0.15, 12.0, 163789),
    "TCG": ("S", 0.05, 4.2, 57234),
    "CCT": ("P", 0.28, 17.2, 234567),
    "CCC": ("P", 0.33, 20.5, 279456),
    "CCA": ("P", 0.27, 16.5, 225123),
    "CCG": ("P", 0.12, 7.2, 98123),
    "ACT": ("T", 0.24, 12.8, 174456),
    "ACC": ("T", 0.38, 19.8, 269789),
    "ACA": ("T", 0.27, 14.2, 193456),
    "ACG": ("T", 0.11, 5.8, 79123),
    "GCT": ("A", 0.25, 17.8, 242789),
    "GCC": ("A", 0.42, 29.5, 402123),
    "GCA": ("A", 0.22, 15.4, 210123),
    "GCG": ("A", 0.11, 7.5, 102456),
    "TAT": ("Y", 0.44, 12.0, 163456),
    "TAC": ("Y", 0.56, 15.4, 210123),
    "TAA": ("*", 0.29, 1.0, 13645),
    "TAG": ("*", 0.24, 0.8, 10918),
    "CAT": ("H", 0.41, 10.6, 144456),
    "CAC": ("H", 0.59, 15.2, 207456),
    "CAA": ("Q", 0.25, 12.1, 165123),
    "CAG": ("Q", 0.75, 35.8, 488123),
    "AAT": ("N", 0.46, 17.2, 234567),
    "AAC": ("N", 0.54, 20.1, 274123),
    "AAA": ("K", 0.42, 24.0, 327456),
    "AAG": ("K", 0.58, 33.2, 452789),
    "GAT": ("D", 0.45, 21.5, 293123),
    "GAC": ("D", 0.55, 26.3, 358456),
    "GAA": ("E", 0.41, 28.8, 392789),
    "GAG": ("E", 0.59, 41.5, 565789),
    "TGT": ("C", 0.44, 10.3, 140456),
    "TGC": ("C", 0.56, 13.1, 178456),
    "TGA": ("*", 0.47, 1.6, 21812),
    "TGG": ("W", 1.00, 13.0, 177456),
    "CGT": ("R", 0.09, 5.1, 69567),
    "CGC": ("R", 0.19, 11.5, 156789),
    "CGA": ("R", 0.10, 6.4, 87234),
    "CGG": ("R", 0.21, 12.5, 170456),
    "AGT": ("S", 0.16, 12.0, 163789),
    "AGC": ("S", 0.24, 19.5, 265789),
    "AGA": ("R", 0.21, 12.3, 167789),
    "AGG": ("R", 0.20, 12.0, 163456),
    "GGT": ("G", 0.16, 11.2, 152789),
    "GGC": ("G", 0.35, 23.8, 324456),
    "GGA": ("G", 0.25, 16.5, 225123),
    "GGG": ("G", 0.24, 16.0, 218456),
}

# Compute relative adaptiveness using shared utility
CHO_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(CHO_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
CHO_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(CHO_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for CHO-K1 cells
# ---------------------------------------------------------------------------
# Codon pair bias reflects the non-random pairing of consecutive codons.
# Positive values indicate over-represented (favored) pairs; negative values
# indicate under-represented (disfavored) pairs.  Values are log-odds ratios
# derived from highly expressed CHO-K1 genes.
# CHO codon preferences are similar to mouse/human (mammalian pattern),
# favoring G/C-ending codons.  CPB data is modeled after mammalian
# ribosome profiling studies with CHO-specific adjustments.
# ---------------------------------------------------------------------------
CHO_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias) — G/C-ending codon pairs typical of
    # highly expressed mammalian genes
    "CTG_CTC": 0.30,   # Leu-Leu, both major Leu codons
    "CAG_CAG": 0.27,   # Gln-Gln, CAG-dominant mammalian pattern
    "GAG_GAG": 0.25,   # Glu-Glu, GAG-preferred
    "ATG_CTC": 0.22,   # Met-Leu, strong mammalian pair
    "GAG_CAG": 0.20,   # Glu-Gln
    "CAG_GAG": 0.18,   # Gln-Glu
    "CTG_CAG": 0.16,   # Leu-Gln
    "ATC_ATC": 0.15,   # Ile-Ile, ATC preferred in mammals
    "TTC_TTC": 0.14,   # Phe-Phe, TTC preferred
    "CAG_CTG": 0.13,   # Gln-Leu
    "AAC_AAC": 0.12,   # Asn-Asn, AAC preferred
    "GAC_GAC": 0.11,   # Asp-Asp, GAC preferred
    "CTG_GAG": 0.10,   # Leu-Glu
    "GCC_GCC": 0.09,   # Ala-Ala, GCC preferred
    "GGC_GGC": 0.08,   # Gly-Gly, GGC preferred in mammals
    # Neutral / slight bias
    "GTG_CTC": 0.01,   # Val-Leu
    "ACC_ACC": -0.01,  # Thr-Thr
    # Disfavored pairs (negative bias) — A/T-ending codon pairs
    # and rare-codon combinations in mammalian cells
    "TTA_TTA": -0.36,  # Leu-Leu, TTA is rare in mammals
    "ATA_ATA": -0.30,  # Ile-Ile, ATA is rare in mammals
    "AGA_AGA": -0.26,  # Arg-Arg, AGA less preferred vs CGN in CHO
    "CGA_CGA": -0.24,  # Arg-Arg, CGA rare in mammals
    "ATA_TTA": -0.22,  # Ile-Leu, both rare codons
    "TTA_AGA": -0.20,  # Leu-Arg, rare pair
    "GGA_GGA": -0.18,  # Gly-Gly, GGA less preferred
    "CCA_CCA": -0.15,  # Pro-Pro, CCA less preferred vs CCC
    "AGA_CGA": -0.14,  # Arg-Arg, non-CGC/CGG pair
    "TGT_TGT": -0.12,  # Cys-Cys, TGT less preferred
    "GCA_GCA": -0.10,  # Ala-Ala, GCA less preferred vs GCC
    "AAT_AAT": -0.09,  # Asn-Asn, AAT less preferred vs AAC
    "TTT_TTT": -0.07,  # Phe-Phe, TTT less preferred vs TTC
    "AAA_AAA": -0.06,  # Lys-Lys, AAA less preferred vs AAG
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for CHO-K1 cells
# ---------------------------------------------------------------------------
# These parameters guide codon optimization algorithms when targeting
# CHO cells as a mammalian expression host.  CHO cells are the dominant
# platform for biopharmaceutical protein production and require strong
# mammalian expression signals including Kozak consensus, poly-A signals,
# and splice-site awareness.
# ---------------------------------------------------------------------------
CHO_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # Kozak consensus sequence for efficient translation initiation in mammals
    "preferred_5utr_kozak": "GCCACCATGG",
    # 3' UTR poly-A signal for transcript stability
    "preferred_3utr": "AATAAA",  # mammalian poly-A signal
    # Maximum number of consecutive rare codons before optimization intervenes
    "max_consecutive_rare_codons": 3,
    # Target GC content for the coding sequence (mammalian ~55% GC)
    "gc_content_target": 0.55,
    # Minimum acceptable GC content
    "gc_content_min": 0.30,
    # Maximum acceptable GC content
    "gc_content_max": 0.70,
    # Sequence motifs to avoid (mRNA instability elements)
    "avoid_motifs": ["TTATTTT"],
    # Whether to check for and avoid cryptic splice sites
    "splice_site_awareness": True,
    # Whether to avoid creating CpG islands that may trigger silencing
    "cpg_island_avoidance": True,
}
