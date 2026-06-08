"""
Arabidopsis thaliana (Thale Cress) Codon Usage Data

Source: Kazusa Codon Usage Database
Arabidopsis thaliana [gbinv]: 29,378 CDSs, 10,949,632 codons
Coding GC: 44.31%

A. thaliana is the primary model organism for plant molecular biology
and biotechnology.  Its codon usage reflects a moderately GC-rich
dicot pattern: GC-ending codons are preferred for most amino acids,
consistent with the genome-wide GC composition of ~44% coding GC.
Key features:
  - Strong preference for GCT (Ala), TCT (Ser), CCT (Pro), ACT (Thr)
    and other NNT codons, reflecting the AT/GC balance of dicot genes.
  - AGA is the preferred Arg codon (consistent with plant pattern).
  - GTG and GGC are preferred for Val and Gly respectively.
  - Stop codon: TAA most common, followed by TGA.

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 3702
  - Liu et al. (2004) Plant Physiol 135:1014-1024
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "ARABIDOPSIS_CODON_USAGE",
    "ARABIDOPSIS_CODON_ADAPTIVENESS",
    "ARABIDOPSIS_PREFERRED_CODONS",
    "ARABIDOPSIS_CODON_PAIR_BIAS",
    "ARABIDOPSIS_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# Source: Kazusa Codon Usage Database, Arabidopsis thaliana
# 29,378 CDSs, 10,949,632 codons, coding GC: 44.31%
# Counts derived from per-thousand frequencies × total codons.
ARABIDOPSIS_CODON_USAGE: CodonUsageTable = {
    # Phe: TTC preferred (GC-ending, consistent with dicot pattern)
    "TTT": ("F", 0.43, 15.4, 168625),
    "TTC": ("F", 0.57, 20.5, 224468),
    # Leu: CTG preferred; TTA very rare in highly expressed genes
    "TTA": ("L", 0.10, 9.6, 105117),
    "TTG": ("L", 0.15, 14.8, 162055),
    "CTT": ("L", 0.16, 15.9, 174100),
    "CTC": ("L", 0.18, 17.4, 190524),
    "CTA": ("L", 0.08, 7.8, 85407),
    "CTG": ("L", 0.33, 32.1, 351470),
    # Ile: ATC preferred
    "ATT": ("I", 0.39, 16.3, 178480),
    "ATC": ("I", 0.50, 21.1, 231018),
    "ATA": ("I", 0.11, 4.7, 51463),
    # Met: only one codon
    "ATG": ("M", 1.00, 22.4, 245272),
    # Val: GTG preferred
    "GTT": ("V", 0.23, 14.8, 162055),
    "GTC": ("V", 0.23, 14.9, 163150),
    "GTA": ("V", 0.12, 7.8, 85407),
    "GTG": ("V", 0.42, 27.3, 298925),
    # Ser: TCT strongly preferred; TCG rare
    "TCT": ("S", 0.27, 19.4, 212423),
    "TCC": ("S", 0.20, 14.3, 156580),
    "TCA": ("S", 0.17, 12.2, 133585),
    "TCG": ("S", 0.06, 4.3, 47083),
    "AGT": ("S", 0.15, 10.8, 118256),
    "AGC": ("S", 0.15, 10.7, 117161),
    # Pro: CCA preferred; CCG rare
    "CCT": ("P", 0.30, 16.2, 177384),
    "CCC": ("P", 0.20, 10.8, 118256),
    "CCA": ("P", 0.38, 20.7, 226658),
    "CCG": ("P", 0.12, 6.4, 70078),
    # Thr: ACT preferred; ACG rare
    "ACT": ("T", 0.33, 17.6, 192713),
    "ACC": ("T", 0.24, 13.0, 142345),
    "ACA": ("T", 0.31, 16.6, 181767),
    "ACG": ("T", 0.12, 6.2, 67888),
    # Ala: GCT preferred; GCG moderate
    "GCT": ("A", 0.35, 22.3, 244178),
    "GCC": ("A", 0.23, 14.9, 163150),
    "GCA": ("A", 0.27, 17.4, 190524),
    "GCG": ("A", 0.15, 9.8, 107307),
    # Tyr: TAC preferred
    "TAT": ("Y", 0.44, 11.7, 128112),
    "TAC": ("Y", 0.56, 14.9, 163150),
    # Stop codons: TAA dominant in plants
    "TAA": ("*", 0.48, 1.4, 15330),
    "TAG": ("*", 0.20, 0.6, 6570),
    # His: CAC preferred
    "CAT": ("H", 0.42, 9.7, 106213),
    "CAC": ("H", 0.58, 13.4, 146725),
    # Gln: CAG preferred
    "CAA": ("Q", 0.37, 16.5, 180673),
    "CAG": ("Q", 0.63, 28.1, 307685),
    # Asn: AAC preferred
    "AAT": ("N", 0.44, 15.2, 166436),
    "AAC": ("N", 0.56, 19.4, 212423),
    # Lys: AAG preferred
    "AAA": ("K", 0.42, 24.5, 268165),
    "AAG": ("K", 0.58, 33.9, 371173),
    # Asp: GAC preferred
    "GAT": ("D", 0.44, 21.7, 237590),
    "GAC": ("D", 0.56, 27.7, 303265),
    # Glu: GAG preferred
    "GAA": ("E", 0.42, 29.8, 326300),
    "GAG": ("E", 0.58, 41.1, 450231),
    # Cys: TGC preferred
    "TGT": ("C", 0.42, 7.4, 81028),
    "TGC": ("C", 0.58, 10.3, 112782),
    # Stop codon
    "TGA": ("*", 0.32, 0.9, 9855),
    # Trp: only one codon
    "TGG": ("W", 1.00, 11.6, 127016),
    # Arg: AGA strongly preferred in plants; CGG/CGA rare
    "CGT": ("R", 0.09, 5.1, 55843),
    "CGC": ("R", 0.07, 3.8, 41609),
    "CGA": ("R", 0.06, 3.6, 39419),
    "CGG": ("R", 0.06, 3.4, 37229),
    "AGA": ("R", 0.49, 28.1, 307685),
    "AGG": ("R", 0.23, 13.3, 145625),
    # Gly: GGC preferred
    "GGT": ("G", 0.22, 15.3, 167533),
    "GGC": ("G", 0.34, 23.6, 258412),
    "GGA": ("G", 0.27, 18.8, 205854),
    "GGG": ("G", 0.17, 11.7, 128112),
}

# Compute relative adaptiveness using shared utility
ARABIDOPSIS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(ARABIDOPSIS_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
ARABIDOPSIS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(ARABIDOPSIS_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for Arabidopsis thaliana
# ---------------------------------------------------------------------------
# Codon pair bias reflects the non-random pairing of consecutive codons.
# Plant codon pair bias is less well-characterized than mammalian/E. coli
# but follows general eukaryotic patterns with plant-specific features:
#   - GC-ending codon pairs are generally favored
#   - AGA-rich pairs are common (AGA is the dominant Arg codon in plants)
#   - AT-rich rare codon pairs are disfavored
# Reference: Quax et al. (2015); plant-specific adjustments from
# Duret & Galtier (2009) and follow-up studies.
# ---------------------------------------------------------------------------
ARABIDOPSIS_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias) — GC-ending codon pairs and AGA pairs
    "CAG_CAG": 0.22,   # Gln-Gln, CAG-dominant plant pattern
    "GAG_GAG": 0.20,   # Glu-Glu, GAG-preferred
    "CTG_CTC": 0.18,   # Leu-Leu, GC-rich preferred Leu codons
    "GAC_GAC": 0.16,   # Asp-Asp, GAC preferred
    "ATG_CTG": 0.14,   # Met-Leu
    "GAG_CAG": 0.13,   # Glu-Gln
    "CAG_CTG": 0.12,   # Gln-Leu
    "AGA_AGA": 0.11,   # Arg-Arg, AGA strongly preferred in plants
    "CTG_CAG": 0.10,   # Leu-Gln
    "GGC_GGC": 0.09,   # Gly-Gly, GGC preferred
    "ATC_ATC": 0.08,   # Ile-Ile, ATC preferred
    "TTC_CTC": 0.07,   # Phe-Leu
    "GCC_GCC": 0.06,   # Ala-Ala
    "AAC_AAC": 0.05,   # Asn-Asn, AAC preferred
    # Neutral / slight bias
    "GCT_GCT": 0.01,   # Ala-Ala, GCT most common Ala
    "ACT_ACT": -0.01,  # Thr-Thr, ACT most common Thr
    # Disfavored pairs (negative bias) — AT-rich and rare codon pairs
    "TTA_TTA": -0.32,  # Leu-Leu, TTA very rare in plants
    "ATA_ATA": -0.28,  # Ile-Ile, ATA rare in plants
    "CGA_CGA": -0.24,  # Arg-Arg, CGA very rare (AGA preferred)
    "CGG_CGG": -0.22,  # Arg-Arg, CGG rare in plants
    "TCG_TCG": -0.20,  # Ser-Ser, TCG rare in plants
    "CCG_CCG": -0.18,  # Pro-Pro, CCG rare in plants
    "GCG_GCG": -0.16,  # Ala-Ala, GCG less preferred
    "TTA_CGA": -0.15,  # Leu-Arg, both rare
    "ACG_ACG": -0.14,  # Thr-Thr, ACG rare in plants
    "ATA_TTA": -0.12,  # Ile-Leu, both rare
    "CTA_CTA": -0.10,  # Leu-Leu, CTA less preferred
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for Arabidopsis thaliana
# ---------------------------------------------------------------------------
# These parameters guide codon optimization algorithms when targeting
# A. thaliana as a plant expression host.  Key considerations:
#   - Plant cryptic splice sites must be avoided (GT-AG rule)
#   - CpG methylation can silence transgenes in plants
#   - Polyadenylation signals must be avoided in CDS
#   - Plant-specific mRNA instability elements (AU-rich elements)
# ---------------------------------------------------------------------------
ARABIDOPSIS_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # 5' UTR context: consensus for efficient translation initiation in plants
    # Kozak-like context (plants use a similar but distinct consensus)
    "preferred_5utr_kozak": "AACATGG",  # Plant consensus around ATG
    # Promoter: 35S CaMV promoter is standard for Arabidopsis
    "preferred_promoter": "35S CaMV",
    # Maximum number of consecutive rare codons before optimization intervenes
    "max_consecutive_rare_codons": 3,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content for the coding sequence (A. thaliana ~44% GC)
    "gc_content_target": 0.44,
    # Minimum acceptable GC content
    "gc_content_min": 0.30,
    # Maximum acceptable GC content
    "gc_content_max": 0.60,
    # Sequence motifs to avoid
    "avoid_motifs": [
        "ATTTA",           # mRNA instability element
        "AATAAA",          # Cryptic polyadenylation signal
        "AATAAT",          # Variant poly-A signal
        "GTAGGT",          # Cryptic splice donor
    ],
    # Plant-specific: avoid CpG dinucleotides (methylation silencing)
    "cpg_island_avoidance": True,
    # Plant-specific: avoid cryptic splice sites
    "splice_site_awareness": True,
    # Maximum length of consecutive T runs (prevents premature poly-A-like signals)
    "max_t_run": 5,
    # Whether codon pair bias should be considered during optimization
    "codon_pair_awareness": True,
}
