"""
Drosophila melanogaster (Fruit Fly) Codon Usage Data

Source: Kazusa Codon Usage Database
Drosophila melanogaster [gbinv]: 23,056 CDSs, 13,740,762 codons
Coding GC: ~54.1%

D. melanogaster is the most widely used model organism in insect genetics
and developmental biology.  Its codon usage reflects a moderately GC-rich
dipteran pattern with strong translational selection in highly expressed
genes.

Key features:
  - GC-ending codons are generally preferred (consistent with ~54% coding GC)
  - AGA is the preferred Arg codon (consistent with insect pattern)
  - GCT preferred for Ala, TCT for Ser — insect consensus
  - CCG is rare for Pro (consistent with Drosophila genome)
  - Drosophila has one of the best-characterized codon bias patterns
    among insects due to decades of population genetics research
  - Highly expressed genes show very strong codon bias (high Fop)

References:
  - Drosophila codon usage: Sharp & Li (1986) Mol Biol Evol 3:125-131
  - Vicario et al. (2007) PLoS Genet 3(2):e16
  - Heger & Ponting (2007) Genome Res 17:538-549
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 7227
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "DROSOPHILA_CODON_USAGE",
    "DROSOPHILA_CODON_ADAPTIVENESS",
    "DROSOPHILA_PREFERRED_CODONS",
    "DROSOPHILA_CODON_PAIR_BIAS",
    "DROSOPHILA_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# Source: Kazusa Codon Usage Database, Drosophila melanogaster
# 23,056 CDSs, 13,740,762 codons, coding GC: ~54.1%
# Counts derived from per-thousand frequencies × total codons.
DROSOPHILA_CODON_USAGE: CodonUsageTable = {
    # Phe: TTC preferred (GC-ending, Drosophila pattern)
    "TTT": ("F", 0.38, 14.5, 199241),
    "TTC": ("F", 0.62, 23.7, 325636),
    # Leu: CTG strongly preferred; TTA very rare in highly expressed genes
    "TTA": ("L", 0.06, 6.0, 82445),
    "TTG": ("L", 0.11, 10.8, 148400),
    "CTT": ("L", 0.10, 9.8, 134661),
    "CTC": ("L", 0.22, 21.2, 291306),
    "CTA": ("L", 0.06, 5.9, 81071),
    "CTG": ("L", 0.45, 43.5, 597723),
    # Ile: ATC preferred; ATA rare
    "ATT": ("I", 0.33, 14.6, 200609),
    "ATC": ("I", 0.56, 25.1, 344873),
    "ATA": ("I", 0.11, 4.9, 67330),
    # Met: only one codon
    "ATG": ("M", 1.00, 23.3, 320160),
    # Val: GTG preferred
    "GTT": ("V", 0.18, 11.2, 153896),
    "GTC": ("V", 0.24, 15.1, 207486),
    "GTA": ("V", 0.09, 5.8, 79696),
    "GTG": ("V", 0.49, 30.4, 417719),
    # Ser: TCT preferred; TCG and AGT rare
    "TCT": ("S", 0.24, 17.0, 233593),
    "TCC": ("S", 0.24, 16.8, 230845),
    "TCA": ("S", 0.13, 9.1, 125041),
    "TCG": ("S", 0.05, 3.6, 49467),
    "AGT": ("S", 0.12, 8.3, 114074),
    "AGC": ("S", 0.22, 15.7, 215710),
    # Pro: CCA preferred; CCG rare
    "CCT": ("P", 0.24, 13.8, 189623),
    "CCC": ("P", 0.28, 16.2, 222600),
    "CCA": ("P", 0.34, 19.6, 269320),
    "CCG": ("P", 0.14, 8.1, 111301),
    # Thr: ACC preferred; ACG rare
    "ACT": ("T", 0.24, 13.2, 181358),
    "ACC": ("T", 0.40, 22.1, 303691),
    "ACA": ("T", 0.23, 12.5, 171759),
    "ACG": ("T", 0.13, 7.1, 97559),
    # Ala: GCT preferred; GCG moderate
    "GCT": ("A", 0.30, 20.5, 281685),
    "GCC": ("A", 0.30, 20.8, 285808),
    "GCA": ("A", 0.22, 15.4, 211608),
    "GCG": ("A", 0.18, 12.2, 167637),
    # Tyr: TAC preferred
    "TAT": ("Y", 0.38, 11.2, 153896),
    "TAC": ("Y", 0.62, 18.3, 251456),
    # Stop codons: TAA dominant in Drosophila
    "TAA": ("*", 0.50, 1.4, 19237),
    "TAG": ("*", 0.22, 0.6, 8245),
    # His: CAC preferred
    "CAT": ("H", 0.38, 9.8, 134661),
    "CAC": ("H", 0.62, 16.0, 219852),
    # Gln: CAG strongly preferred
    "CAA": ("Q", 0.28, 14.3, 196493),
    "CAG": ("Q", 0.72, 36.8, 505660),
    # Asn: AAC preferred
    "AAT": ("N", 0.38, 14.2, 195120),
    "AAC": ("N", 0.62, 23.3, 320160),
    # Lys: AAG preferred
    "AAA": ("K", 0.38, 23.1, 317414),
    "AAG": ("K", 0.62, 37.8, 519361),
    # Asp: GAC preferred
    "GAT": ("D", 0.40, 21.8, 299549),
    "GAC": ("D", 0.60, 32.9, 452111),
    # Glu: GAG strongly preferred
    "GAA": ("E", 0.38, 30.1, 413617),
    "GAG": ("E", 0.62, 49.1, 674631),
    # Cys: TGC preferred
    "TGT": ("C", 0.37, 7.7, 105803),
    "TGC": ("C", 0.63, 13.1, 180004),
    # Stop codon
    "TGA": ("*", 0.28, 0.8, 10993),
    # Trp: only one codon
    "TGG": ("W", 1.00, 11.8, 162141),
    # Arg: AGA strongly preferred in insects; CGA/CGG very rare
    "CGT": ("R", 0.10, 5.8, 79696),
    "CGC": ("R", 0.13, 7.6, 104429),
    "CGA": ("R", 0.04, 2.2, 30230),
    "CGG": ("R", 0.05, 3.1, 42597),
    "AGA": ("R", 0.49, 28.8, 395734),
    "AGG": ("R", 0.19, 11.0, 151148),
    # Gly: GGC preferred
    "GGT": ("G", 0.18, 12.1, 166263),
    "GGC": ("G", 0.37, 24.6, 338023),
    "GGA": ("G", 0.26, 17.2, 236321),
    "GGG": ("G", 0.19, 12.8, 175882),
}

# Compute relative adaptiveness using shared utility
DROSOPHILA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(DROSOPHILA_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
DROSOPHILA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(DROSOPHILA_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for Drosophila melanogaster
# ---------------------------------------------------------------------------
# Drosophila has one of the best-characterized codon pair biases among
# insects, owing to extensive population genetics and ribosome profiling
# studies.  Key features:
#   - GC-ending codon pairs are generally favored (reflects ~54% coding GC)
#   - AGA-rich pairs are common (AGA is the dominant Arg codon)
#   - Rare codon pairs (ATA, TTA, CGA, TCG, CCG, ACG) are disfavored
#   - Codon pair bias is correlated with translation elongation rates
# Reference: Quax et al. (2015); Drosophila-specific data from
# Heger & Ponting (2007) and vicario et al. (2007).
# ---------------------------------------------------------------------------
DROSOPHILA_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias) — GC-ending codon pairs
    "CTG_CTC": 0.24,   # Leu-Leu, GC-rich preferred codons
    "CTG_CTG": 0.22,   # Leu-Leu, CTG dominant
    "CAG_CAG": 0.21,   # Gln-Gln, CAG dominant in Drosophila
    "GAG_GAG": 0.20,   # Glu-Glu, GAG strongly preferred
    "CTG_CAG": 0.18,   # Leu-Gln, GC-rich pair
    "GAG_CAG": 0.16,   # Glu-Gln, charged pair
    "CAG_CTG": 0.15,   # Gln-Leu
    "GAC_GAC": 0.14,   # Asp-Asp, GAC preferred
    "ATG_CTG": 0.13,   # Met-Leu
    "AGA_AGA": 0.12,   # Arg-Arg, AGA strongly preferred in insects
    "GGC_GGC": 0.10,   # Gly-Gly, GGC preferred
    "TTC_CTC": 0.08,   # Phe-Leu
    "ATC_ATC": 0.07,   # Ile-Ile, ATC preferred
    "GCC_GCC": 0.06,   # Ala-Ala, GCC preferred
    "ACC_ACC": 0.05,   # Thr-Thr, ACC preferred
    # Neutral / slight bias
    "GCT_GCT": 0.01,   # Ala-Ala, GCT common
    "ACT_ACT": -0.01,  # Thr-Thr, ACT common
    # Disfavored pairs (negative bias) — rare codon pairs
    "TTA_TTA": -0.35,  # Leu-Leu, TTA very rare in Drosophila
    "ATA_ATA": -0.30,  # Ile-Ile, ATA rare
    "CGA_CGA": -0.26,  # Arg-Arg, CGA very rare (AGA preferred)
    "CGG_CGG": -0.24,  # Arg-Arg, CGG rare
    "TCG_TCG": -0.22,  # Ser-Ser, TCG rare
    "CCG_CCG": -0.20,  # Pro-Pro, CCG rare in Drosophila
    "GCG_GCG": -0.18,  # Ala-Ala, GCG less preferred
    "ACG_ACG": -0.16,  # Thr-Thr, ACG rare
    "TTA_CGA": -0.15,  # Leu-Arg, both rare
    "CTA_CTA": -0.13,  # Leu-Leu, CTA less preferred
    "ATA_TTA": -0.11,  # Ile-Leu, both rare
    "CTA_CCG": -0.09,  # Leu-Pro, both less preferred
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for Drosophila melanogaster
# ---------------------------------------------------------------------------
# D. melanogaster is the standard insect model organism.  When used as
# an expression host (typically S2 cell line), key considerations:
#   - Drosophila S2 cells are widely used for recombinant protein expression
#   - Less stringent GT avoidance than mammals (insect introns are AT-rich)
#   - CpG methylation is minimal in Drosophila (virtually absent)
#   - Gal4/UAS system is the standard expression driver
#   - Polyadenylation signals must be avoided in CDS
#   - Drosophila genes show very strong codon bias in highly expressed genes
# ---------------------------------------------------------------------------
DROSOPHILA_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # Expression system: Gal4/UAS (standard Drosophila expression system)
    "preferred_promoter": "Gal4/UAS / Actin5C",
    # Maximum number of consecutive rare codons
    "max_consecutive_rare_codons": 3,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content (Drosophila coding GC ~54%)
    "gc_content_target": 0.54,
    "gc_content_min": 0.35,
    "gc_content_max": 0.65,
    # Sequence motifs to avoid
    "avoid_motifs": [
        "ATTTA",           # mRNA instability element
        "AATAAA",          # Cryptic polyadenylation signal
        "AATAAT",          # Variant poly-A signal
    ],
    # Insect-specific: less stringent splice site avoidance than mammals
    # Drosophila introns are typically short (60-80 nt) and AT-rich
    "splice_site_awareness": True,
    # Drosophila has virtually no DNA methylation
    "cpg_island_avoidance": False,
    # Maximum length of consecutive T runs
    "max_t_run": 5,
    # Codon pair bias awareness
    "codon_pair_awareness": True,
}
