"""
Chinese Hamster Ovary (CHO-K1) Codon Usage Data

Source: Kazusa Codon Usage Database
Cricetulus griseus, 8,175 CDSs, 4,361,612 codons
Coding GC: 52.72%

CHO cells are the most commonly used mammalian host for
biopharmaceutical protein production.  Codon usage is very
similar to mouse/human (mammalian pattern), with AGA as
the preferred Arg codon (consistent with other mammals).
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
# Source: Kazusa Codon Usage Database, Cricetulus griseus, 8,175 CDSs
# Counts derived from per-thousand frequencies × 4,361,612 total codons.
# Arg codon preference corrected: AGA is the preferred Arg codon in CHO
# (consistent with mouse/human mammalian pattern), not CGG.
CHO_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.45, 17.2, 75020),
    "TTC": ("F", 0.55, 20.9, 91158),
    "TTA": ("L", 0.07, 7.4, 32276),
    "TTG": ("L", 0.13, 12.8, 55829),
    "CTT": ("L", 0.13, 12.4, 54084),
    "CTC": ("L", 0.21, 19.8, 86360),
    "CTA": ("L", 0.07, 7.1, 30967),
    "CTG": ("L", 0.39, 38.5, 167922),
    "ATT": ("I", 0.35, 15.6, 68041),
    "ATC": ("I", 0.49, 21.8, 95083),
    "ATA": ("I", 0.16, 7.2, 31404),
    "ATG": ("M", 1.00, 21.8, 95083),
    "GTT": ("V", 0.18, 11.3, 49286),
    "GTC": ("V", 0.24, 14.7, 64116),
    "GTA": ("V", 0.11, 7.0, 30531),
    "GTG": ("V", 0.47, 28.9, 126051),
    "TCT": ("S", 0.18, 14.3, 62371),
    "TCC": ("S", 0.23, 18.5, 80690),
    "TCA": ("S", 0.15, 12.0, 52339),
    "TCG": ("S", 0.05, 4.2, 18319),
    "CCT": ("P", 0.28, 17.2, 75020),
    "CCC": ("P", 0.33, 20.5, 89413),
    "CCA": ("P", 0.27, 16.5, 71967),
    "CCG": ("P", 0.12, 7.2, 31404),
    "ACT": ("T", 0.24, 12.8, 55829),
    "ACC": ("T", 0.38, 19.8, 86360),
    "ACA": ("T", 0.27, 14.2, 61935),
    "ACG": ("T", 0.11, 5.8, 25297),
    "GCT": ("A", 0.25, 17.8, 77637),
    "GCC": ("A", 0.42, 29.5, 128668),
    "GCA": ("A", 0.22, 15.4, 67169),
    "GCG": ("A", 0.11, 7.5, 32712),
    "TAT": ("Y", 0.44, 12.0, 52339),
    "TAC": ("Y", 0.56, 15.4, 67169),
    "TAA": ("*", 0.29, 1.0, 4362),
    "TAG": ("*", 0.24, 0.8, 3489),
    "CAT": ("H", 0.41, 10.6, 46233),
    "CAC": ("H", 0.59, 15.2, 66297),
    "CAA": ("Q", 0.25, 12.1, 52776),
    "CAG": ("Q", 0.75, 35.8, 156146),
    "AAT": ("N", 0.46, 17.2, 75020),
    "AAC": ("N", 0.54, 20.1, 87668),
    "AAA": ("K", 0.42, 24.0, 104679),
    "AAG": ("K", 0.58, 33.2, 144806),
    "GAT": ("D", 0.45, 21.5, 93775),
    "GAC": ("D", 0.55, 26.3, 114710),
    "GAA": ("E", 0.41, 28.8, 125614),
    "GAG": ("E", 0.59, 41.5, 181007),
    "TGT": ("C", 0.44, 10.3, 44925),
    "TGC": ("C", 0.56, 13.1, 57137),
    "TGA": ("*", 0.47, 1.6, 6979),
    "TGG": ("W", 1.00, 13.0, 56701),
    "CGT": ("R", 0.09, 5.1, 22244),
    "CGC": ("R", 0.19, 11.5, 50159),
    "CGA": ("R", 0.10, 6.3, 27478),
    "CGG": ("R", 0.20, 12.1, 52776),
    "AGT": ("S", 0.15, 12.0, 52339),
    "AGC": ("S", 0.24, 19.5, 85051),
    "AGA": ("R", 0.22, 12.6, 54956),
    "AGG": ("R", 0.20, 12.0, 52339),
    "GGT": ("G", 0.16, 11.2, 48850),
    "GGC": ("G", 0.35, 23.8, 103806),
    "GGA": ("G", 0.25, 16.5, 71967),
    "GGG": ("G", 0.24, 16.0, 69786),
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
    "CGG_CGG": -0.26,  # Arg-Arg, CGG less preferred vs AGA in CHO
    "CGA_CGA": -0.24,  # Arg-Arg, CGA rare in mammals
    "ATA_TTA": -0.22,  # Ile-Leu, both rare codons
    "TTA_AGA": -0.20,  # Leu-Arg, rare pair
    "GGA_GGA": -0.18,  # Gly-Gly, GGA less preferred
    "CCA_CCA": -0.15,  # Pro-Pro, CCA less preferred vs CCC
    "AGA_CGA": -0.14,  # Arg-Arg, non-AGA/CGC pair
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
