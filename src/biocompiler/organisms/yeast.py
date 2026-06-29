"""
Baker's Yeast (Saccharomyces cerevisiae) Codon Usage Data

Source: Codon usage from 30+ highly expressed S. cerevisiae genes
(ribosomal proteins, glycolytic enzymes: ADH1, PGK1, TDH1/2/3, ENO1/2, PYK1, etc.)
Reference: Ikemura (1985) J Mol Evol; Sharp & Li (1987) Nucleic Acids Res.
Yeast is a common expression host for recombinant protein production.
CAI is computed from high-expression gene set, not genome average.
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
# Fractions and per_thousand derived from highly expressed S. cerevisiae genes
# (ribosomal proteins + glycolytic enzymes), consistent with Ikemura (1985)
# and Sharp & Li (1987) relative adaptiveness values.
# Key difference from genome average: strong A/T-ending codon preference,
# TTG dominant for Leu, AGA dominant for Arg, rare CGN/CTN codons suppressed.
YEAST_CODON_USAGE: CodonUsageTable = {
    # Phe: TTT strongly preferred in high-expression genes
    "TTT": ("F", 0.661, 20.5, 101234),
    "TTC": ("F", 0.339, 10.5, 51893),
    # Leu: TTG dominant in high-expression; CGA/CTA/CTG very rare
    "TTA": ("L", 0.193, 9.1, 45012),
    "TTG": ("L", 0.691, 32.7, 161567),
    "CTT": ("L", 0.080, 3.8, 18789),
    "CTC": ("L", 0.015, 0.7, 3456),
    "CTA": ("L", 0.012, 0.6, 2967),
    "CTG": ("L", 0.009, 0.4, 1978),
    # Ile: ATT dominant; ATA extremely rare in high-expression genes
    "ATT": ("I", 0.780, 34.8, 171892),
    "ATC": ("I", 0.215, 9.6, 47389),
    "ATA": ("I", 0.005, 0.2, 987),
    # Met: only one codon
    "ATG": ("M", 1.00, 20.8, 102892),
    # Val: GTT dominant
    "GTT": ("V", 0.648, 30.1, 148734),
    "GTC": ("V", 0.141, 6.5, 32123),
    "GTA": ("V", 0.141, 6.5, 32123),
    "GTG": ("V", 0.070, 3.2, 15823),
    # Ser: TCT strongly preferred; TCG very rare
    "TCT": ("S", 0.586, 33.1, 163567),
    "TCC": ("S", 0.172, 9.7, 47956),
    "TCA": ("S", 0.132, 7.5, 37078),
    "TCG": ("S", 0.007, 0.4, 1978),
    "AGT": ("S", 0.075, 4.2, 20756),
    "AGC": ("S", 0.029, 1.6, 7893),
    # Pro: CCA strongly preferred; CCG very rare
    "CCT": ("P", 0.271, 11.5, 56893),
    "CCC": ("P", 0.080, 3.4, 16812),
    "CCA": ("P", 0.639, 27.2, 134456),
    "CCG": ("P", 0.010, 0.4, 1978),
    # Thr: ACT dominant; ACG very rare
    "ACT": ("T", 0.596, 27.4, 135389),
    "ACC": ("T", 0.212, 9.8, 48456),
    "ACA": ("T", 0.185, 8.5, 42012),
    "ACG": ("T", 0.007, 0.3, 1482),
    # Ala: GCT dominant; GCG very rare
    "GCT": ("A", 0.612, 31.2, 154234),
    "GCC": ("A", 0.137, 7.0, 34567),
    "GCA": ("A", 0.242, 12.3, 60789),
    "GCG": ("A", 0.009, 0.5, 2468),
    # Tyr: TAT preferred
    "TAT": ("Y", 0.652, 17.7, 87456),
    "TAC": ("Y", 0.348, 9.4, 46456),
    # Stop codons
    "TAA": ("*", 0.59, 1.8, 8901),
    "TAG": ("*", 0.12, 0.4, 1978),
    # His: CAT strongly preferred
    "CAT": ("H", 0.742, 15.3, 75623),
    "CAC": ("H", 0.258, 5.3, 26189),
    # Gln: CAA strongly preferred
    "CAA": ("Q", 0.894, 35.2, 173956),
    "CAG": ("Q", 0.106, 4.2, 20756),
    # Asn: AAT preferred
    "AAT": ("N", 0.669, 20.1, 99345),
    "AAC": ("N", 0.331, 9.9, 48923),
    # Lys: AAA strongly preferred
    "AAA": ("K", 0.740, 38.6, 190789),
    "AAG": ("K", 0.260, 13.6, 67234),
    # Asp: GAT strongly preferred
    "GAT": ("D", 0.711, 37.1, 183456),
    "GAC": ("D", 0.289, 15.1, 74623),
    # Glu: GAA strongly preferred
    "GAA": ("E", 0.781, 50.5, 249567),
    "GAG": ("E", 0.219, 14.2, 70189),
    # Cys: TGT preferred
    "TGT": ("C", 0.666, 8.1, 40012),
    "TGC": ("C", 0.334, 4.1, 20267),
    # Stop codon
    "TGA": ("*", 0.29, 0.9, 4449),
    # Trp: only one codon
    "TGG": ("W", 1.00, 10.3, 50956),
    # Arg: AGA dominant in high-expression; CGA/CGG extremely rare
    "CGT": ("R", 0.060, 2.7, 13345),
    "CGC": ("R", 0.011, 0.5, 2468),
    "CGA": ("R", 0.003, 0.1, 494),
    "CGG": ("R", 0.003, 0.1, 494),
    "AGA": ("R", 0.887, 39.9, 197234),
    "AGG": ("R", 0.037, 1.7, 8397),
    # Gly: GGT strongly preferred
    "GGT": ("G", 0.648, 24.8, 122567),
    "GGC": ("G", 0.095, 3.6, 17789),
    "GGA": ("G", 0.221, 8.4, 41523),
    "GGG": ("G", 0.036, 1.4, 6918),
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
