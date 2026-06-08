"""
Oryza sativa (Rice) Codon Usage Data

Source: Kazusa Codon Usage Database
Oryza sativa [gbinv]: 42,091 CDSs, 21,741,635 codons
Coding GC: ~55.2%

O. sativa (rice) is the most important cereal crop and a model organism
for monocot molecular biology and biotechnology.  Its codon usage reflects
a GC-rich monocot pattern that differs significantly from dicots like
Arabidopsis:
  - GC-ending codons are STRONGLY preferred for most amino acids
  - CCG is much more common for Pro than in dicots
  - GCG is more common for Ala than in dicots
  - Coding GC (~55%) is substantially higher than Arabidopsis (~44%)
  - This GC bias is characteristic of grasses (Poaceae) and reflects
    the GC-heterogeneous isochore structure of the rice genome

Key features:
  - Strong preference for GC-ending codons across most amino acid families
  - CCG is a major Pro codon (unlike Arabidopsis where CCG is rare)
  - GCG is a significant Ala codon (unlike Arabidopsis where GCG is minor)
  - GCC is the preferred Ala codon (not GCT as in Arabidopsis)
  - AGA is the preferred Arg codon (consistent with plant pattern)
  - Important host for transgenic cereal crop development

References:
  - Liu et al. (2004) Plant Physiol 135:1014-1024
  - Wang & Hickey (2007) Mol Biol Evol 24:423-432
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 4530
  - Oryza sativa japonica group (Kazusa ID for japonica cultivar)
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "ORYZA_CODON_USAGE",
    "ORYZA_CODON_ADAPTIVENESS",
    "ORYZA_PREFERRED_CODONS",
    "ORYZA_CODON_PAIR_BIAS",
    "ORYZA_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# Source: Kazusa Codon Usage Database, Oryza sativa
# 42,091 CDSs, 21,741,635 codons, coding GC: ~55.2%
# Counts derived from per-thousand frequencies × total codons.
ORYZA_CODON_USAGE: CodonUsageTable = {
    # Phe: TTC strongly preferred (GC-ending, monocot pattern)
    "TTT": ("F", 0.36, 13.2, 286992),
    "TTC": ("F", 0.64, 23.4, 508754),
    # Leu: CTG strongly preferred; TTA very rare
    "TTA": ("L", 0.06, 5.6, 121753),
    "TTG": ("L", 0.11, 10.4, 226113),
    "CTT": ("L", 0.12, 11.1, 241332),
    "CTC": ("L", 0.24, 22.4, 487013),
    "CTA": ("L", 0.06, 5.8, 126101),
    "CTG": ("L", 0.41, 38.4, 834879),
    # Ile: ATC strongly preferred (GC-ending)
    "ATT": ("I", 0.31, 13.5, 293512),
    "ATC": ("I", 0.58, 25.2, 547890),
    "ATA": ("I", 0.11, 4.8, 104360),
    # Met: only one codon
    "ATG": ("M", 1.00, 22.8, 495710),
    # Val: GTG strongly preferred; GTA rare
    "GTT": ("V", 0.17, 11.0, 239158),
    "GTC": ("V", 0.24, 15.6, 339170),
    "GTA": ("V", 0.08, 5.3, 115231),
    "GTG": ("V", 0.51, 33.2, 721822),
    # Ser: TCT and TCC both common; TCG moderate in rice (GC-rich)
    "TCT": ("S", 0.20, 14.4, 313079),
    "TCC": ("S", 0.24, 17.2, 373956),
    "TCA": ("S", 0.12, 8.7, 189133),
    "TCG": ("S", 0.10, 7.3, 158720),
    "AGT": ("S", 0.12, 8.9, 193481),
    "AGC": ("S", 0.22, 16.0, 347866),
    # Pro: CCA and CCG both common (CCG much more common than in dicots)
    "CCT": ("P", 0.21, 13.2, 286992),
    "CCC": ("P", 0.21, 13.4, 291339),
    "CCA": ("P", 0.30, 19.2, 417439),
    "CCG": ("P", 0.28, 17.8, 387002),
    # Thr: ACC preferred; ACG moderate in rice (GC-rich)
    "ACT": ("T", 0.22, 13.7, 297861),
    "ACC": ("T", 0.39, 24.3, 528322),
    "ACA": ("T", 0.21, 13.0, 282641),
    "ACG": ("T", 0.18, 11.2, 243506),
    # Ala: GCC preferred; GCG significant in rice (monocot GC-rich pattern)
    "GCT": ("A", 0.22, 16.0, 347866),
    "GCC": ("A", 0.33, 24.0, 521800),
    "GCA": ("A", 0.20, 14.7, 319522),
    "GCG": ("A", 0.25, 18.4, 400006),
    # Tyr: TAC preferred (GC-ending)
    "TAT": ("Y", 0.38, 10.6, 230421),
    "TAC": ("Y", 0.62, 17.3, 376130),
    # Stop codons: TAA dominant in plants
    "TAA": ("*", 0.46, 1.2, 26090),
    "TAG": ("*", 0.22, 0.6, 13045),
    # His: CAC preferred
    "CAT": ("H", 0.36, 8.7, 189133),
    "CAC": ("H", 0.64, 15.5, 337095),
    # Gln: CAG strongly preferred
    "CAA": ("Q", 0.28, 16.1, 350001),
    "CAG": ("Q", 0.72, 41.6, 904452),
    # Asn: AAC preferred (GC-ending)
    "AAT": ("N", 0.37, 13.7, 297861),
    "AAC": ("N", 0.63, 23.3, 506561),
    # Lys: AAG strongly preferred (GC-ending)
    "AAA": ("K", 0.35, 21.0, 456574),
    "AAG": ("K", 0.65, 39.0, 847924),
    # Asp: GAC preferred (GC-ending)
    "GAT": ("D", 0.38, 19.5, 423962),
    "GAC": ("D", 0.62, 31.7, 689211),
    # Glu: GAG strongly preferred (GC-ending)
    "GAA": ("E", 0.37, 28.0, 608763),
    "GAG": ("E", 0.63, 47.9, 1041436),
    # Cys: TGC preferred (GC-ending)
    "TGT": ("C", 0.36, 7.0, 152191),
    "TGC": ("C", 0.64, 12.5, 271770),
    # Stop codon
    "TGA": ("*", 0.32, 0.8, 17393),
    # Trp: only one codon
    "TGG": ("W", 1.00, 11.5, 250028),
    # Arg: AGA preferred (consistent with plants); CGC and CGG moderate
    "CGT": ("R", 0.09, 5.6, 121753),
    "CGC": ("R", 0.13, 8.2, 178281),
    "CGA": ("R", 0.05, 3.1, 67399),
    "CGG": ("R", 0.08, 5.0, 108708),
    "AGA": ("R", 0.46, 28.7, 623986),
    "AGG": ("R", 0.19, 12.1, 263074),
    # Gly: GGC preferred (GC-ending)
    "GGT": ("G", 0.17, 12.2, 265248),
    "GGC": ("G", 0.37, 26.3, 571805),
    "GGA": ("G", 0.25, 17.8, 387002),
    "GGG": ("G", 0.21, 15.1, 328299),
}

# Compute relative adaptiveness using shared utility
ORYZA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(ORYZA_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
ORYZA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(ORYZA_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for Oryza sativa
# ---------------------------------------------------------------------------
# Rice codon pair bias reflects the strong GC-rich monocot pattern.
# Monocot CPB is less well-characterized than mammalian/E. coli but
# follows general eukaryotic patterns with grass-specific features:
#   - GC-ending codon pairs are strongly favored (reflects ~55% coding GC)
#   - CCG pairs are common (CCG is a major Pro codon in rice, unlike dicots)
#   - GCG pairs are more acceptable than in dicots
#   - AGA-rich pairs are common (AGA is the dominant Arg codon)
#   - AT-rich rare codon pairs (TTA, ATA, CGA, TCG) are strongly disfavored
# Reference: Quax et al. (2015); monocot-specific adjustments from
# Liu et al. (2004) and Wang & Hickey (2007).
# ---------------------------------------------------------------------------
ORYZA_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias) — GC-ending codon pairs, rice-specific
    "CAG_CAG": 0.25,   # Gln-Gln, CAG strongly preferred in rice
    "GAG_GAG": 0.23,   # Glu-Glu, GAG strongly preferred
    "CTG_CTC": 0.21,   # Leu-Leu, GC-rich preferred Leu codons
    "CTG_CTG": 0.20,   # Leu-Leu, CTG dominant
    "GAC_GAC": 0.18,   # Asp-Asp, GAC preferred
    "CAG_CTG": 0.16,   # Gln-Leu, GC-rich pair
    "CTG_CAG": 0.15,   # Leu-Gln
    "GAG_CAG": 0.14,   # Glu-Gln, charged pair
    "ATG_CTG": 0.13,   # Met-Leu
    "AGA_AGA": 0.12,   # Arg-Arg, AGA preferred in plants
    "GGC_GGC": 0.10,   # Gly-Gly, GGC preferred
    "GCC_GCC": 0.09,   # Ala-Ala, GCC preferred in rice
    "CCG_CCG": 0.07,   # Pro-Pro, CCG significant in rice (monocot-specific)
    "ATC_ATC": 0.06,   # Ile-Ile, ATC preferred
    "GCG_GCC": 0.05,   # Ala-Ala, GC-rich pair
    # Neutral / slight bias
    "GCT_GCT": 0.01,   # Ala-Ala
    "CCA_CCA": -0.01,  # Pro-Pro, CCA preferred
    # Disfavored pairs (negative bias) — AT-rich rare codon pairs
    "TTA_TTA": -0.36,  # Leu-Leu, TTA very rare in rice
    "ATA_ATA": -0.32,  # Ile-Ile, ATA rare
    "CGA_CGA": -0.27,  # Arg-Arg, CGA very rare (AGA preferred)
    "CGG_CGG": -0.24,  # Arg-Arg, CGG less common in rice
    "TCG_TCG": -0.22,  # Ser-Ser, TCG less common in rice (but more than dicots)
    "GCG_GCG": -0.15,  # Ala-Ala, GCG moderate in rice (less disfavored than dicots)
    "ACG_ACG": -0.14,  # Thr-Thr, ACG moderate in rice
    "TTA_CGA": -0.16,  # Leu-Arg, both rare
    "ATA_TTA": -0.13,  # Ile-Leu, both rare
    "CTA_CTA": -0.11,  # Leu-Leu, CTA less preferred
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for Oryza sativa
# ---------------------------------------------------------------------------
# O. sativa (rice) is the primary model for monocot molecular biology
# and biotechnology.  Key considerations:
#   - Plant cryptic splice sites must be avoided (GT-AG rule)
#   - CpG methylation can silence transgenes in rice (epigenetic silencing)
#   - Polyadenylation signals must be avoided in CDS
#   - Monocot-specific promoter elements (e.g., maize Ubiquitin-1 promoter)
#   - Rice has strong codon bias in highly expressed genes
#   - Target GC is higher than for dicots (~55% vs ~44%)
# ---------------------------------------------------------------------------
ORYZA_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # 5' UTR context: consensus for efficient translation initiation in monocots
    "preferred_5utr_kozak": "AACATGG",  # Plant consensus around ATG
    # Promoter: maize Ubiquitin-1 promoter is standard for rice transformation
    "preferred_promoter": "maize Ubi-1 / Actin1",
    # Maximum number of consecutive rare codons before optimization intervenes
    "max_consecutive_rare_codons": 3,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content for the coding sequence (O. sativa ~55% GC)
    "gc_content_target": 0.55,
    # Minimum acceptable GC content
    "gc_content_min": 0.35,
    # Maximum acceptable GC content
    "gc_content_max": 0.70,
    # Sequence motifs to avoid
    "avoid_motifs": [
        "ATTTA",           # mRNA instability element
        "AATAAA",          # Cryptic polyadenylation signal
        "AATAAT",          # Variant poly-A signal
        "GTAGGT",          # Cryptic splice donor
    ],
    # Plant-specific: CpG methylation silencing is a major concern in rice
    "cpg_island_avoidance": True,
    # Plant-specific: avoid cryptic splice sites
    "splice_site_awareness": True,
    # Maximum length of consecutive T runs (prevents premature poly-A-like signals)
    "max_t_run": 5,
    # Whether codon pair bias should be considered during optimization
    "codon_pair_awareness": True,
}
