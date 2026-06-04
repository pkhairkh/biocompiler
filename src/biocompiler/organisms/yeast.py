"""
Baker's Yeast (Saccharomyces cerevisiae) Codon Usage Data

Source: Kazusa Codon Usage Database
Yeast is a common expression host for recombinant protein production.
High-expression genes used for CAI computation.
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "YEAST_CODON_USAGE",
    "YEAST_CODON_ADAPTIVENESS",
    "YEAST_PREFERRED_CODONS",
    "YEAST_CODON_PAIR_BIAS",
    "YEAST_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
YEAST_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.59, 18.3, 90456),
    "TTC": ("F", 0.41, 12.7, 62789),
    "TTA": ("L", 0.28, 13.2, 65234),
    "TTG": ("L", 0.11, 5.2, 25678),
    "CTT": ("L", 0.13, 6.1, 30123),
    "CTC": ("L", 0.05, 2.4, 11856),
    "CTA": ("L", 0.14, 6.7, 33123),
    "CTG": ("L", 0.29, 13.7, 67789),
    "ATT": ("I", 0.46, 20.6, 101892),
    "ATC": ("I", 0.30, 13.3, 65789),
    "ATA": ("I", 0.24, 10.7, 52893),
    "ATG": ("M", 1.00, 20.8, 102892),
    "GTT": ("V", 0.39, 18.2, 89956),
    "GTC": ("V", 0.21, 9.8, 48456),
    "GTA": ("V", 0.17, 8.0, 39567),
    "GTG": ("V", 0.22, 10.4, 51345),
    "TCT": ("S", 0.26, 14.7, 72678),
    "TCC": ("S", 0.16, 9.0, 44567),
    "TCA": ("S", 0.21, 11.9, 58789),
    "TCG": ("S", 0.10, 5.6, 27678),
    "CCT": ("P", 0.31, 13.2, 65234),
    "CCC": ("P", 0.15, 6.4, 31678),
    "CCA": ("P", 0.42, 17.9, 88456),
    "CCG": ("P", 0.12, 5.1, 25234),
    "ACT": ("T", 0.35, 16.1, 79567),
    "ACC": ("T", 0.22, 10.1, 49893),
    "ACA": ("T", 0.30, 13.8, 68234),
    "ACG": ("T", 0.13, 6.0, 29678),
    "GCT": ("A", 0.36, 18.4, 90956),
    "GCC": ("A", 0.22, 11.2, 55345),
    "GCA": ("A", 0.29, 14.8, 73123),
    "GCG": ("A", 0.13, 6.6, 32678),
    "TAT": ("Y", 0.56, 15.2, 75123),
    "TAC": ("Y", 0.44, 11.9, 58893),
    "TAA": ("*", 0.48, 1.7, 8412),
    "TAG": ("*", 0.22, 0.8, 3956),
    "CAT": ("H", 0.64, 13.2, 65234),
    "CAC": ("H", 0.36, 7.4, 36567),
    "CAA": ("Q", 0.69, 27.2, 134456),
    "CAG": ("Q", 0.31, 12.2, 60345),
    "AAT": ("N", 0.59, 17.7, 87567),
    "AAC": ("N", 0.41, 12.3, 60789),
    "AAA": ("K", 0.58, 30.3, 149789),
    "AAG": ("K", 0.42, 21.9, 108234),
    "GAT": ("D", 0.64, 33.4, 165123),
    "GAC": ("D", 0.36, 18.8, 92956),
    "GAA": ("E", 0.70, 45.3, 223789),
    "GAG": ("E", 0.30, 19.4, 95893),
    "TGT": ("C", 0.63, 7.7, 38078),
    "TGC": ("C", 0.37, 4.5, 22234),
    "TGA": ("*", 0.30, 1.1, 5432),
    "TGG": ("W", 1.00, 10.3, 50956),
    "CGT": ("R", 0.15, 6.8, 33678),
    "CGC": ("R", 0.06, 2.7, 13345),
    "CGA": ("R", 0.07, 3.2, 15823),
    "CGG": ("R", 0.04, 1.8, 8893),
    "AGT": ("S", 0.15, 8.5, 42078),
    "AGC": ("S", 0.12, 6.8, 33678),
    "AGA": ("R", 0.48, 21.5, 106234),
    "AGG": ("R", 0.20, 9.0, 44567),
    "GGT": ("G", 0.48, 18.3, 90456),
    "GGC": ("G", 0.19, 7.3, 36078),
    "GGA": ("G", 0.21, 8.0, 39567),
    "GGG": ("G", 0.12, 4.6, 22734),
}

# Compute relative adaptiveness using shared utility
YEAST_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(YEAST_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
YEAST_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(YEAST_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for S. cerevisiae
# ---------------------------------------------------------------------------
# Codon pair bias reflects the non-random pairing of consecutive codons.
# Positive values indicate over-represented (favored) pairs; negative values
# indicate under-represented (disfavored) pairs.  Values are log-odds ratios
# derived from highly expressed S. cerevisiae genes.
# Reference: rare-codon clustering and CPB effects on translational efficiency
# in yeast (based on Coleman et al. and follow-up studies).
# ---------------------------------------------------------------------------
YEAST_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias)
    "AAA_AAA": 0.32,   # Lys-Lys, A-rich, common in yeast genes
    "GAA_GAA": 0.28,   # Glu-Glu, highly expressed gene signature
    "GAA_GAT": 0.21,   # Glu-Asp
    "GAT_GAA": 0.19,   # Asp-Glu
    "AAA_GAA": 0.18,   # Lys-Glu
    "GAA_AAA": 0.17,   # Glu-Lys
    "AAT_GAT": 0.14,   # Asn-Asp
    "TTT_GAA": 0.12,   # Phe-Glu
    "GAT_AAA": 0.11,   # Asp-Lys
    "GCT_GCT": 0.09,   # Ala-Ala
    "ATT_GAA": 0.08,   # Ile-Glu
    "TCT_GCT": 0.07,   # Ser-Ala
    "GAA_AAT": 0.06,   # Glu-Asn
    "CCT_GAA": 0.05,   # Pro-Glu
    "AGA_AGA": 0.04,   # Arg-Arg (AGA is preferred Arg in yeast)
    # Neutral / slight bias
    "GTT_GCT": 0.00,   # Val-Ala
    "ACT_ACT": -0.01,  # Thr-Thr
    # Disfavored pairs (negative bias)
    "CGG_CGG": -0.35,  # Arg-Arg, CGG is rare in yeast
    "CTC_CTC": -0.28,  # Leu-Leu, CTC rare in yeast
    "CGC_CGC": -0.26,  # Arg-Arg, CGC rare in yeast
    "CCG_CCG": -0.24,  # Pro-Pro, CCG rare in yeast
    "GCG_GCG": -0.22,  # Ala-Ala, GCG rare in yeast
    "CGG_CGC": -0.20,  # Arg-Arg, rare codon pairing
    "CTC_CTA": -0.19,  # Leu-Leu, non-preferred pair
    "TCG_TCG": -0.17,  # Ser-Ser, TCG rare in yeast
    "ATA_ATA": -0.15,  # Ile-Ile, ATA is minor Ile codon
    "CGA_CGG": -0.14,  # Arg-Arg, rare pair
    "ACG_ACG": -0.12,  # Thr-Thr, ACG rare in yeast
    "GGG_GGG": -0.10,  # Gly-Gly, GGG rare in yeast
    "CAT_CGC": -0.09,  # His-Arg, includes rare Arg codon
    "CTA_CTC": -0.08,  # Leu-Leu, non-preferred pair
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for S. cerevisiae
# ---------------------------------------------------------------------------
# These parameters guide codon optimization algorithms when targeting
# yeast as an expression host.  S. cerevisiae has a notably AT-rich genome
# (~62% AT) which drives its codon preferences toward A/T-ending codons.
# ---------------------------------------------------------------------------
YEAST_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # 5' UTR context: consensus for efficient translation initiation in yeast
    "preferred_5utr": "AAATATCTTT",
    # Promoter elements: TATA-box and transcription start site context
    "preferred_promoter": "TATA-box + transcription start",
    # Maximum number of consecutive rare codons before optimization intervenes
    "max_consecutive_rare_codons": 2,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content for the coding sequence (S. cerevisiae ~38% GC)
    "gc_content_target": 0.38,
    # Minimum acceptable GC content
    "gc_content_min": 0.25,
    # Maximum acceptable GC content
    "gc_content_max": 0.55,
    # Sequence motifs to avoid (instability / degradation signals)
    "avoid_motifs": ["ATTTA"],
    # Maximum length of consecutive T runs (polypyrimidine tract effects)
    "max_t_run": 6,
    # Whether codon pair bias should be considered during optimization
    "codon_pair_awareness": True,
}
