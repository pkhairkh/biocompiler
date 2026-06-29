"""
Spodoptera frugiperda Sf9 (Fall Armyworm) Codon Usage Data

Source: Codon usage from S. frugiperda genome and Sf9 cell line
sequencing projects.
NCBI TaxID: 7108
Coding GC: ~43.5%

Sf9 cells (derived from S. frugiperda ovary) are one of the two
primary hosts for baculovirus-mediated protein expression (the
Baculovirus Expression Vector System, BEVS).  Codon usage reflects
a lepidopteran insect pattern with moderate GC content.

Key features:
  - Codon preferences are intermediate between mammalian and Drosophila
  - G/C-ending codons moderately preferred (consistent with ~43% GC)
  - AGA is the preferred Arg codon (consistent with insects)
  - TCT preferred for Ser, consistent with insect pattern
  - Less stringent GT avoidance than mammals (insect introns are AT-rich)
  - Important host for recombinant protein production, especially for
    proteins requiring post-translational modifications

References:
  - Kawai & Oka (1999) based on Sf9 cell line codon usage
  - Posey et al. (2020) Genome Biol Evol 12:1664-1678
  - NCBI TaxID: 7108
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "SPODOPTERA_CODON_USAGE",
    "SPODOPTERA_CODON_ADAPTIVENESS",
    "SPODOPTERA_PREFERRED_CODONS",
    "SPODOPTERA_CODON_PAIR_BIAS",
    "SPODOPTERA_EXPRESSION_OPTIMIZATION_PARAMS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# Source: Codon usage from S. frugiperda genome annotation
# and Sf9 cell line expressed sequence data.
# ~14,500 CDSs; coding GC ~43.5%
# Per-thousand values derived from published lepidopteran codon frequencies.
SPODOPTERA_CODON_USAGE: CodonUsageTable = {
    # Phe: TTC moderately preferred
    "TTT": ("F", 0.44, 15.0, 26040),
    "TTC": ("F", 0.56, 19.2, 33337),
    # Leu: CTG preferred; TTA and CTA rare
    "TTA": ("L", 0.09, 7.8, 13541),
    "TTG": ("L", 0.14, 12.4, 21531),
    "CTT": ("L", 0.15, 13.2, 22921),
    "CTC": ("L", 0.19, 17.0, 29512),
    "CTA": ("L", 0.07, 6.3, 10933),
    "CTG": ("L", 0.36, 32.0, 55572),
    # Ile: ATC preferred; ATA rare
    "ATT": ("I", 0.38, 16.0, 27776),
    "ATC": ("I", 0.50, 21.3, 36988),
    "ATA": ("I", 0.12, 5.2, 9028),
    # Met: only one codon
    "ATG": ("M", 1.00, 21.0, 36469),
    # Val: GTG preferred
    "GTT": ("V", 0.22, 13.4, 23260),
    "GTC": ("V", 0.22, 13.6, 23607),
    "GTA": ("V", 0.12, 7.4, 12844),
    "GTG": ("V", 0.44, 27.0, 46885),
    # Ser: TCT preferred; TCG rare
    "TCT": ("S", 0.26, 18.0, 31251),
    "TCC": ("S", 0.20, 13.9, 24133),
    "TCA": ("S", 0.16, 11.2, 19444),
    "TCG": ("S", 0.05, 3.7, 6422),
    "AGT": ("S", 0.16, 11.0, 19097),
    "AGC": ("S", 0.17, 12.0, 20833),
    # Pro: CCA preferred; CCG rare
    "CCT": ("P", 0.28, 14.5, 25173),
    "CCC": ("P", 0.20, 10.4, 18056),
    "CCA": ("P", 0.39, 20.5, 35590),
    "CCG": ("P", 0.13, 6.6, 11458),
    # Thr: ACT preferred; ACG rare
    "ACT": ("T", 0.32, 16.0, 27776),
    "ACC": ("T", 0.25, 12.5, 21702),
    "ACA": ("T", 0.31, 15.6, 27081),
    "ACG": ("T", 0.12, 6.0, 10417),
    # Ala: GCT preferred; GCG moderate
    "GCT": ("A", 0.34, 20.8, 36110),
    "GCC": ("A", 0.23, 14.1, 24479),
    "GCA": ("A", 0.27, 16.5, 28646),
    "GCG": ("A", 0.16, 9.8, 17014),
    # Tyr: TAC preferred
    "TAT": ("Y", 0.44, 10.7, 18577),
    "TAC": ("Y", 0.56, 13.7, 23784),
    # Stop codons: TAA common in insects
    "TAA": ("*", 0.46, 1.2, 2083),
    "TAG": ("*", 0.22, 0.6, 1042),
    # His: CAC preferred
    "CAT": ("H", 0.42, 9.3, 16146),
    "CAC": ("H", 0.58, 12.9, 22396),
    # Gln: CAG preferred
    "CAA": ("Q", 0.35, 15.0, 26040),
    "CAG": ("Q", 0.65, 27.9, 48437),
    # Asn: AAC preferred
    "AAT": ("N", 0.44, 14.6, 25346),
    "AAC": ("N", 0.56, 18.6, 32290),
    # Lys: AAG preferred
    "AAA": ("K", 0.41, 22.8, 39583),
    "AAG": ("K", 0.59, 32.9, 57118),
    # Asp: GAC preferred
    "GAT": ("D", 0.43, 20.2, 35068),
    "GAC": ("D", 0.57, 26.7, 46353),
    # Glu: GAG preferred
    "GAA": ("E", 0.42, 28.1, 48783),
    "GAG": ("E", 0.58, 38.8, 67360),
    # Cys: TGC preferred
    "TGT": ("C", 0.43, 7.2, 12499),
    "TGC": ("C", 0.57, 9.5, 16493),
    # Stop codon
    "TGA": ("*", 0.32, 0.8, 1389),
    # Trp: only one codon
    "TGG": ("W", 1.00, 11.4, 19791),
    # Arg: AGA strongly preferred in insects; CGA/CGG rare
    "CGT": ("R", 0.10, 5.7, 9896),
    "CGC": ("R", 0.08, 4.5, 7813),
    "CGA": ("R", 0.06, 3.3, 5729),
    "CGG": ("R", 0.06, 3.5, 6076),
    "AGA": ("R", 0.48, 27.6, 47915),
    "AGG": ("R", 0.22, 12.4, 21528),
    # Gly: GGC preferred
    "GGT": ("G", 0.22, 13.8, 23957),
    "GGC": ("G", 0.34, 21.3, 36988),
    "GGA": ("G", 0.27, 17.0, 29512),
    "GGG": ("G", 0.17, 10.7, 18577),
}

# Compute relative adaptiveness using shared utility
SPODOPTERA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(SPODOPTERA_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
SPODOPTERA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(SPODOPTERA_CODON_USAGE)

# ---------------------------------------------------------------------------
# Codon Pair Bias (CPB) for Spodoptera frugiperda Sf9
# ---------------------------------------------------------------------------
# Insect codon pair bias is less well-studied than mammalian/E. coli
# but follows general eukaryotic patterns with insect-specific features:
#   - AGA-rich pairs are common (AGA is the dominant Arg codon)
#   - GC-ending codon pairs generally favored
#   - Less pronounced pair bias than mammals
#   - Based on lepidopteran ribosome profiling studies and
#     extrapolation from Drosophila CPB data with Sf9 adjustments.
# ---------------------------------------------------------------------------
SPODOPTERA_CODON_PAIR_BIAS: dict[str, float] = {
    # Favored pairs (positive bias)
    "CAG_CAG": 0.19,   # Gln-Gln
    "GAG_GAG": 0.17,   # Glu-Glu, GAG preferred
    "CTG_CTC": 0.16,   # Leu-Leu, GC-rich
    "GAC_GAC": 0.14,   # Asp-Asp, GAC preferred
    "ATG_CTG": 0.13,   # Met-Leu
    "AGA_AGA": 0.12,   # Arg-Arg, AGA strongly preferred in insects
    "GAG_CAG": 0.11,   # Glu-Gln
    "CTG_CAG": 0.10,   # Leu-Gln
    "CAG_CTG": 0.09,   # Gln-Leu
    "GGC_GGC": 0.07,   # Gly-Gly
    "ATC_ATC": 0.06,   # Ile-Ile, ATC preferred
    "TTC_CTC": 0.05,   # Phe-Leu
    "AAC_AAC": 0.04,   # Asn-Asn, AAC preferred
    # Neutral / slight bias
    "GCT_GCT": 0.01,   # Ala-Ala
    "ACT_ACT": -0.01,  # Thr-Thr
    # Disfavored pairs (negative bias)
    "TTA_TTA": -0.28,  # Leu-Leu, TTA rare in insects
    "ATA_ATA": -0.25,  # Ile-Ile, ATA rare
    "CGA_CGA": -0.22,  # Arg-Arg, CGA rare
    "CGG_CGG": -0.20,  # Arg-Arg, CGG rare
    "TCG_TCG": -0.18,  # Ser-Ser, TCG rare
    "CCG_CCG": -0.16,  # Pro-Pro, CCG rare
    "GCG_GCG": -0.14,  # Ala-Ala, GCG less preferred
    "ACG_ACG": -0.12,  # Thr-Thr, ACG rare
    "TTA_CGA": -0.11,  # Leu-Arg, both rare
    "ATA_TTA": -0.10,  # Ile-Leu, both rare
    "CTA_CTA": -0.08,  # Leu-Leu, CTA less preferred
}

# ---------------------------------------------------------------------------
# Expression Optimization Parameters for Spodoptera frugiperda Sf9
# ---------------------------------------------------------------------------
# Sf9 cells are used with the Baculovirus Expression Vector System (BEVS).
# Key considerations:
#   - Less stringent GT avoidance than mammals (insect introns are AT-rich
#     and use different splice signals)
#   - CpG methylation less of a concern than in plants/mammals
#   - Recombinant protein expression is baculovirus-driven (polyhedrin
#     or p10 promoters), not mammalian Kozak-dependent
#   - mRNA stability elements still relevant
#   - Post-translational modifications differ from mammalian cells
# ---------------------------------------------------------------------------
SPODOPTERA_EXPRESSION_OPTIMIZATION_PARAMS: dict = {
    # Baculovirus expression: polyhedrin or p10 promoter-driven
    "preferred_promoter": "polyhedrin / p10",
    # Maximum number of consecutive rare codons
    "max_consecutive_rare_codons": 4,
    # Fraction threshold below which a codon is considered "rare"
    "rare_codon_threshold": 0.10,
    # Target GC content (Sf9 coding GC ~43.5%)
    "gc_content_target": 0.44,
    "gc_content_min": 0.30,
    "gc_content_max": 0.60,
    # Sequence motifs to avoid
    "avoid_motifs": [
        "ATTTA",           # mRNA instability element
        "AATAAA",          # Polyadenylation signal (can cause premature poly-A)
    ],
    # Insect-specific: less stringent splice site avoidance than mammals
    # Insect introns are short and AT-rich, splice signals differ
    "splice_site_awareness": True,
    # Insects have DNA methylation but less CpG silencing concern
    "cpg_island_avoidance": False,
    # Maximum length of consecutive T runs
    "max_t_run": 6,
    # Codon pair bias awareness
    "codon_pair_awareness": True,
}
