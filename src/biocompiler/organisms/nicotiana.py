"""
Nicotiana benthamiana (Tobacco) Codon Usage Data

Source: Codon usage from N. benthamiana genome sequencing project
and published codon usage tables.
NCBI TaxID: 4100
Coding GC: ~42.7%

N. benthamiana is the most widely used host for transient protein
expression in plants via Agrobacterium infiltration (agroinfiltration).
Its codon usage is very similar to A. thaliana (both are dicots) but
with slightly more AT-rich bias due to the Solanaceae genome
composition.

Key features:
  - Very similar codon preferences to A. thaliana (both dicots)
  - Slightly stronger preference for AT-ending codons than Arabidopsis
  - AGA is the strongly preferred Arg codon (consistent with plants)
  - GCT preferred for Ala, TCT for Ser, consistent with dicot pattern
  - Widely used for recombinant protein production via transient expression

References:
  - Goodin et al. (2008) Plant J 55:337-350
  - Nakasugi et al. (2014) BMC Genomics 15:31
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 4100
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "NICOTIANA_CODON_USAGE",
    "NICOTIANA_CODON_ADAPTIVENESS",
    "NICOTIANA_PREFERRED_CODONS",
    "NICOTIANA_CODON_PAIR_BIAS",
    "NICOTIANA_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# Source: Codon usage computed from N. benthamiana genome annotation
# (Nakasugi et al. 2014, Bombarde et al. updated annotation)
# Estimated from ~25,000 CDSs; coding GC ~42.7%
# Per-thousand values derived from published codon frequencies.
NICOTIANA_CODON_USAGE: CodonUsageTable = {
    # Phe: TTC moderately preferred (similar to Arabidopsis)
    "TTT": ("F", 0.45, 15.2, 42712),
    "TTC": ("F", 0.55, 18.6, 52279),
    # Leu: CTG preferred; TTA rare
    "TTA": ("L", 0.11, 9.4, 26416),
    "TTG": ("L", 0.15, 13.0, 36538),
    "CTT": ("L", 0.16, 14.2, 39911),
    "CTC": ("L", 0.17, 15.1, 42437),
    "CTA": ("L", 0.08, 7.3, 20514),
    "CTG": ("L", 0.33, 29.1, 81791),
    # Ile: ATC preferred
    "ATT": ("I", 0.40, 16.0, 44946),
    "ATC": ("I", 0.48, 19.4, 54513),
    "ATA": ("I", 0.12, 5.0, 14049),
    # Met: only one codon
    "ATG": ("M", 1.00, 21.6, 60692),
    # Val: GTG preferred
    "GTT": ("V", 0.24, 14.5, 40753),
    "GTC": ("V", 0.22, 13.5, 37944),
    "GTA": ("V", 0.13, 7.9, 22197),
    "GTG": ("V", 0.41, 25.2, 70836),
    # Ser: TCT preferred; TCG rare
    "TCT": ("S", 0.28, 18.8, 52843),
    "TCC": ("S", 0.19, 12.8, 35976),
    "TCA": ("S", 0.17, 11.5, 32322),
    "TCG": ("S", 0.05, 3.5, 9834),
    "AGT": ("S", 0.16, 10.9, 30632),
    "AGC": ("S", 0.15, 10.0, 28105),
    # Pro: CCA preferred; CCG rare
    "CCT": ("P", 0.29, 14.8, 41586),
    "CCC": ("P", 0.19, 9.7, 27267),
    "CCA": ("P", 0.39, 19.9, 55937),
    "CCG": ("P", 0.13, 6.5, 18279),
    # Thr: ACT preferred; ACG rare
    "ACT": ("T", 0.33, 16.2, 45527),
    "ACC": ("T", 0.23, 11.4, 32039),
    "ACA": ("T", 0.32, 15.7, 44121),
    "ACG": ("T", 0.12, 5.8, 16303),
    # Ala: GCT preferred; GCG moderate
    "GCT": ("A", 0.36, 21.0, 59036),
    "GCC": ("A", 0.22, 12.8, 35976),
    "GCA": ("A", 0.27, 16.0, 44946),
    "GCG": ("A", 0.15, 8.9, 25014),
    # Tyr: TAC preferred
    "TAT": ("Y", 0.45, 10.7, 30069),
    "TAC": ("Y", 0.55, 13.2, 37096),
    # Stop codons: TAA dominant in plants
    "TAA": ("*", 0.47, 1.3, 3654),
    "TAG": ("*", 0.22, 0.6, 1686),
    # His: CAC preferred
    "CAT": ("H", 0.43, 9.0, 25297),
    "CAC": ("H", 0.57, 11.9, 33444),
    # Gln: CAG preferred
    "CAA": ("Q", 0.38, 15.5, 43564),
    "CAG": ("Q", 0.62, 25.3, 71118),
    # Asn: AAC preferred
    "AAT": ("N", 0.45, 14.6, 41031),
    "AAC": ("N", 0.55, 17.9, 50300),
    # Lys: AAG preferred
    "AAA": ("K", 0.43, 23.5, 66052),
    "AAG": ("K", 0.57, 31.2, 87683),
    # Asp: GAC preferred
    "GAT": ("D", 0.44, 20.3, 57053),
    "GAC": ("D", 0.56, 26.0, 73072),
    # Glu: GAG preferred
    "GAA": ("E", 0.43, 28.3, 79525),
    "GAG": ("E", 0.57, 37.5, 105375),
    # Cys: TGC slightly preferred
    "TGT": ("C", 0.43, 7.2, 20237),
    "TGC": ("C", 0.57, 9.5, 26702),
    # Stop codon
    "TGA": ("*", 0.31, 0.9, 2530),
    # Trp: only one codon
    "TGG": ("W", 1.00, 10.8, 30346),
    # Arg: AGA strongly preferred in plants; CGN codons rare
    "CGT": ("R", 0.08, 4.5, 12649),
    "CGC": ("R", 0.06, 3.4, 9556),
    "CGA": ("R", 0.06, 3.2, 8994),
    "CGG": ("R", 0.05, 2.9, 8151),
    "AGA": ("R", 0.52, 29.2, 82064),
    "AGG": ("R", 0.23, 12.8, 35976),
    # Gly: GGC preferred
    "GGT": ("G", 0.23, 14.5, 40753),
    "GGC": ("G", 0.32, 20.4, 57334),
    "GGA": ("G", 0.28, 17.7, 49745),
    "GGG": ("G", 0.17, 10.8, 30346),
}

# Compute relative adaptiveness using shared utility
NICOTIANA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(NICOTIANA_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
NICOTIANA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(NICOTIANA_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for Nicotiana benthamiana
# ---------------------------------------------------------------------------
# Very similar to Arabidopsis (both dicots), with slightly more AT bias.
# Plant codon pair bias is less well-characterized than mammalian/E. coli
# but N. benthamiana pair bias can be reasonably approximated from
# Arabidopsis data with minor adjustments for Solanaceae-specific patterns.
# ---------------------------------------------------------------------------
NICOTIANA_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias)
    "CAG_CAG": 0.21,   # Gln-Gln, CAG-dominant plant pattern
    "GAG_GAG": 0.19,   # Glu-Glu, GAG-preferred
    "CTG_CTC": 0.17,   # Leu-Leu, GC-rich preferred Leu codons
    "GAC_GAC": 0.15,   # Asp-Asp, GAC preferred
    "AGA_AGA": 0.13,   # Arg-Arg, AGA strongly preferred in plants
    "ATG_CTG": 0.12,   # Met-Leu
    "GAG_CAG": 0.11,   # Glu-Gln
    "CAG_CTG": 0.10,   # Gln-Leu
    "GGC_GGC": 0.08,   # Gly-Gly, GGC preferred
    "ATC_ATC": 0.07,   # Ile-Ile, ATC preferred
    "TTC_CTC": 0.06,   # Phe-Leu
    "GCT_GCT": 0.05,   # Ala-Ala, GCT most common Ala in Nicotiana
    "ACT_ACT": 0.04,   # Thr-Thr, ACT most common Thr
    # Neutral / slight bias
    "GCC_GCC": 0.01,   # Ala-Ala
    "CCA_CCA": -0.01,  # Pro-Pro, CCA preferred
    # Disfavored pairs (negative bias)
    "TTA_TTA": -0.30,  # Leu-Leu, TTA very rare in plants
    "ATA_ATA": -0.27,  # Ile-Ile, ATA rare
    "CGA_CGA": -0.23,  # Arg-Arg, CGA very rare
    "CGG_CGG": -0.21,  # Arg-Arg, CGG rare
    "TCG_TCG": -0.19,  # Ser-Ser, TCG rare
    "CCG_CCG": -0.17,  # Pro-Pro, CCG rare
    "GCG_GCG": -0.15,  # Ala-Ala, GCG less preferred
    "TTA_CGA": -0.14,  # Leu-Arg, both rare
    "ACG_ACG": -0.13,  # Thr-Thr, ACG rare
    "ATA_TTA": -0.11,  # Ile-Leu, both rare
    "CTA_CTA": -0.09,  # Leu-Leu, CTA less preferred
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for Nicotiana benthamiana
# ---------------------------------------------------------------------------
# N. benthamiana is the primary host for transient protein expression
# in plants via agroinfiltration.  Key considerations:
#   - Same plant-specific constraints as Arabidopsis
#   - CpG methylation is a major concern for transgene silencing
#   - Cryptic splice sites must be avoided
#   - Polyadenylation signals in CDS must be eliminated
#   - Typical expression is driven by 35S CaMV or leaf-specific promoters
# ---------------------------------------------------------------------------
NICOTIANA_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # 5' UTR context: plant translation initiation consensus
    "preferred_5utr_kozak": "AACATGG",
    # Promoter: 35S CaMV or leaf-specific promoters for agroinfiltration
    "preferred_promoter": "35S CaMV / Leaf-specific",
    # Maximum number of consecutive rare codons
    "max_consecutive_rare_codons": 3,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content (N. benthamiana ~43% GC)
    "gc_content_target": 0.43,
    "gc_content_min": 0.30,
    "gc_content_max": 0.60,
    # Sequence motifs to avoid
    "avoid_motifs": [
        "ATTTA",           # mRNA instability element
        "AATAAA",          # Cryptic polyadenylation signal
        "AATAAT",          # Variant poly-A signal
        "GTAGGT",          # Cryptic splice donor
    ],
    # Plant-specific: CpG methylation silencing
    "cpg_island_avoidance": True,
    # Plant-specific: cryptic splice site avoidance
    "splice_site_awareness": True,
    # Maximum length of consecutive T runs
    "max_t_run": 5,
    # Codon pair bias awareness
    "codon_pair_awareness": True,
}
