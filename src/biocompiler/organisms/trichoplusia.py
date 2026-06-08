"""
Trichoplusia ni Hi5 (Cabbage Looper) Codon Usage Data

Source: Codon usage from T. ni genome and Hi5 cell line sequencing.
NCBI TaxID: 7111
Coding GC: ~42.0%

Hi5 cells (derived from T. ni ovary, also known as Tn5 or High Five)
are the second major baculovirus expression host alongside Sf9.
Hi5 cells often produce higher recombinant protein yields than Sf9,
particularly for secreted and membrane proteins, making them preferred
for biopharmaceutical production.

Key features:
  - Very similar codon preferences to Sf9 (both lepidopteran insects)
  - Slightly more AT-rich than Sf9 (coding GC ~42% vs ~43.5%)
  - AGA strongly preferred for Arg (consistent with insects)
  - GCT, TCT, CCA preferred for Ala, Ser, Pro respectively
  - Important for secreted protein production due to superior
    secretory pathway capacity vs Sf9

References:
  - Davis et al. (1992) Biotechnol Bioeng 40:142-147
  - Granados et al. (2007) for T. ni genome data
  - NCBI TaxID: 7111
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "TRICHOPLUSIA_CODON_USAGE",
    "TRICHOPLUSIA_CODON_ADAPTIVENESS",
    "TRICHOPLUSIA_PREFERRED_CODONS",
    "TRICHOPLUSIA_CODON_PAIR_BIAS",
    "TRICHOPLUSIA_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# Source: Codon usage from T. ni genome annotation and Hi5 cell line
# expressed sequence data.
# ~13,800 CDSs; coding GC ~42.0%
# Per-thousand values derived from published lepidopteran codon frequencies.
TRICHOPLUSIA_CODON_USAGE: CodonUsageTable = {
    # Phe: TTC moderately preferred (slightly more AT-rich than Sf9)
    "TTT": ("F", 0.46, 15.4, 25840),
    "TTC": ("F", 0.54, 18.2, 30549),
    # Leu: CTG preferred; TTA and CTA rare
    "TTA": ("L", 0.10, 8.4, 14106),
    "TTG": ("L", 0.14, 12.1, 20325),
    "CTT": ("L", 0.16, 13.5, 22679),
    "CTC": ("L", 0.18, 15.3, 25703),
    "CTA": ("L", 0.07, 6.2, 10415),
    "CTG": ("L", 0.35, 29.6, 49718),
    # Ile: ATC preferred; ATA rare
    "ATT": ("I", 0.40, 16.6, 27875),
    "ATC": ("I", 0.48, 20.0, 33587),
    "ATA": ("I", 0.12, 5.0, 8397),
    # Met: only one codon
    "ATG": ("M", 1.00, 20.5, 34428),
    # Val: GTG preferred
    "GTT": ("V", 0.23, 13.6, 22837),
    "GTC": ("V", 0.21, 12.6, 21156),
    "GTA": ("V", 0.13, 7.5, 12601),
    "GTG": ("V", 0.43, 25.4, 42656),
    # Ser: TCT preferred; TCG rare
    "TCT": ("S", 0.27, 18.0, 30228),
    "TCC": ("S", 0.19, 12.8, 21497),
    "TCA": ("S", 0.17, 11.2, 18810),
    "TCG": ("S", 0.05, 3.4, 5712),
    "AGT": ("S", 0.17, 11.2, 18810),
    "AGC": ("S", 0.15, 10.2, 17130),
    # Pro: CCA preferred; CCG rare
    "CCT": ("P", 0.29, 14.6, 24530),
    "CCC": ("P", 0.19, 9.8, 16462),
    "CCA": ("P", 0.40, 20.1, 33754),
    "CCG": ("P", 0.12, 6.2, 10415),
    # Thr: ACT preferred; ACG rare
    "ACT": ("T", 0.34, 16.4, 27539),
    "ACC": ("T", 0.24, 11.8, 19817),
    "ACA": ("T", 0.30, 14.8, 24865),
    "ACG": ("T", 0.12, 5.8, 9742),
    # Ala: GCT preferred; GCG moderate
    "GCT": ("A", 0.36, 21.0, 35276),
    "GCC": ("A", 0.22, 12.8, 21497),
    "GCA": ("A", 0.27, 15.6, 26197),
    "GCG": ("A", 0.15, 8.8, 14782),
    # Tyr: TAC preferred
    "TAT": ("Y", 0.46, 10.8, 18137),
    "TAC": ("Y", 0.54, 12.8, 21497),
    # Stop codons: TAA common in insects
    "TAA": ("*", 0.45, 1.1, 1848),
    "TAG": ("*", 0.23, 0.6, 1008),
    # His: CAC preferred
    "CAT": ("H", 0.43, 9.2, 15454),
    "CAC": ("H", 0.57, 12.2, 20493),
    # Gln: CAG preferred
    "CAA": ("Q", 0.36, 14.8, 24865),
    "CAG": ("Q", 0.64, 26.4, 44337),
    # Asn: AAC preferred
    "AAT": ("N", 0.45, 14.4, 24192),
    "AAC": ("N", 0.55, 17.6, 29562),
    # Lys: AAG preferred
    "AAA": ("K", 0.42, 22.4, 37624),
    "AAG": ("K", 0.58, 30.8, 51723),
    # Asp: GAC preferred
    "GAT": ("D", 0.44, 19.8, 33251),
    "GAC": ("D", 0.56, 25.4, 42656),
    # Glu: GAG preferred
    "GAA": ("E", 0.43, 27.2, 45691),
    "GAG": ("E", 0.57, 36.2, 60812),
    # Cys: TGC slightly preferred
    "TGT": ("C", 0.44, 7.0, 11764),
    "TGC": ("C", 0.56, 8.9, 14952),
    # Stop codon
    "TGA": ("*", 0.32, 0.8, 1344),
    # Trp: only one codon
    "TGG": ("W", 1.00, 11.0, 18467),
    # Arg: AGA strongly preferred in insects; CGN codons rare
    "CGT": ("R", 0.09, 5.2, 8736),
    "CGC": ("R", 0.07, 4.0, 6720),
    "CGA": ("R", 0.06, 3.2, 5376),
    "CGG": ("R", 0.05, 3.0, 5040),
    "AGA": ("R", 0.50, 28.4, 47696),
    "AGG": ("R", 0.23, 13.0, 21840),
    # Gly: GGC preferred
    "GGT": ("G", 0.23, 13.6, 22837),
    "GGC": ("G", 0.33, 19.6, 32904),
    "GGA": ("G", 0.27, 16.2, 27206),
    "GGG": ("G", 0.17, 10.2, 17130),
}

# Compute relative adaptiveness using shared utility
TRICHOPLUSIA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(TRICHOPLUSIA_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
TRICHOPLUSIA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(TRICHOPLUSIA_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for Trichoplusia ni Hi5
# ---------------------------------------------------------------------------
# Very similar to Sf9 (both lepidopteran).  Slightly more AT-rich bias.
# Hi5 cells have been shown to have somewhat different translational
# efficiency profiles than Sf9 for secreted proteins, but codon pair
# bias patterns are largely conserved across lepidopterans.
# ---------------------------------------------------------------------------
TRICHOPLUSIA_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias)
    "CAG_CAG": 0.18,   # Gln-Gln
    "GAG_GAG": 0.16,   # Glu-Glu, GAG preferred
    "CTG_CTC": 0.15,   # Leu-Leu, GC-rich
    "GAC_GAC": 0.13,   # Asp-Asp, GAC preferred
    "ATG_CTG": 0.12,   # Met-Leu
    "AGA_AGA": 0.12,   # Arg-Arg, AGA preferred in insects
    "GAG_CAG": 0.10,   # Glu-Gln
    "CTG_CAG": 0.09,   # Leu-Gln
    "CAG_CTG": 0.08,   # Gln-Leu
    "GGC_GGC": 0.06,   # Gly-Gly
    "ATC_ATC": 0.05,   # Ile-Ile, ATC preferred
    "TTC_CTC": 0.04,   # Phe-Leu
    "GCT_GCT": 0.04,   # Ala-Ala, GCT dominant
    "ACT_ACT": 0.03,   # Thr-Thr, ACT dominant
    # Neutral / slight bias
    "CCA_CCA": 0.01,   # Pro-Pro, CCA dominant
    "GCC_GCC": -0.01,  # Ala-Ala
    # Disfavored pairs (negative bias)
    "TTA_TTA": -0.29,  # Leu-Leu, TTA rare
    "ATA_ATA": -0.26,  # Ile-Ile, ATA rare
    "CGA_CGA": -0.23,  # Arg-Arg, CGA rare
    "CGG_CGG": -0.21,  # Arg-Arg, CGG rare
    "TCG_TCG": -0.19,  # Ser-Ser, TCG rare
    "CCG_CCG": -0.17,  # Pro-Pro, CCG rare
    "GCG_GCG": -0.15,  # Ala-Ala, GCG less preferred
    "ACG_ACG": -0.13,  # Thr-Thr, ACG rare
    "TTA_CGA": -0.12,  # Leu-Arg, both rare
    "ATA_TTA": -0.10,  # Ile-Leu, both rare
    "CTA_CTA": -0.08,  # Leu-Leu, CTA less preferred
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for Trichoplusia ni Hi5
# ---------------------------------------------------------------------------
# Hi5 cells are used with baculovirus for recombinant protein production.
# Key considerations:
#   - Same insect-specific constraints as Sf9
#   - Hi5 cells are especially good for secreted/membrane proteins
#   - Less stringent GT avoidance than mammals
#   - CpG methylation less of a concern
#   - Secretory pathway capacity is superior to Sf9
# ---------------------------------------------------------------------------
TRICHOPLUSIA_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # Baculovirus expression: polyhedrin or p10 promoter-driven
    "preferred_promoter": "polyhedrin / p10",
    # Maximum number of consecutive rare codons
    "max_consecutive_rare_codons": 4,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content (Hi5 coding GC ~42%)
    "gc_content_target": 0.42,
    "gc_content_min": 0.30,
    "gc_content_max": 0.60,
    # Sequence motifs to avoid
    "avoid_motifs": [
        "ATTTA",           # mRNA instability element
        "AATAAA",          # Polyadenylation signal
    ],
    # Insect-specific: less stringent splice site avoidance than mammals
    "splice_site_awareness": True,
    # Insects: less CpG silencing concern than plants/mammals
    "cpg_island_avoidance": False,
    # Maximum length of consecutive T runs
    "max_t_run": 6,
    # Codon pair bias awareness
    "codon_pair_awareness": True,
}
