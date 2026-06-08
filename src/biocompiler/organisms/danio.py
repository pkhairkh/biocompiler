"""
Danio rerio (Zebrafish) Codon Usage Data

Source: Kazusa Codon Usage Database
Danio rerio, ~26,000 CDSs, coding GC: ~53.5%

D. rerio (zebrafish) is one of the most important vertebrate model
organisms for developmental biology, genetics, and drug discovery.
Its codon usage reflects a moderately GC-rich vertebrate pattern
similar to other fish.

Key features:
  - GC-ending codons preferred for most amino acids (fish pattern)
  - AGA preferred for Arg (vertebrate pattern)
  - CTG preferred for Leu
  - Important model for transgenic studies

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 7955
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "DANIO_CODON_USAGE",
    "DANIO_CODON_ADAPTIVENESS",
    "DANIO_PREFERRED_CODONS",
]

DANIO_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.42, 15.8, 316000),
    "TTC": ("F", 0.58, 21.8, 436000),
    "TTA": ("L", 0.07, 7.5, 150000),
    "TTG": ("L", 0.12, 12.5, 250000),
    "CTT": ("L", 0.12, 12.8, 256000),
    "CTC": ("L", 0.21, 22.0, 440000),
    "CTA": ("L", 0.07, 7.2, 144000),
    "CTG": ("L", 0.41, 43.0, 860000),
    "ATT": ("I", 0.34, 15.8, 316000),
    "ATC": ("I", 0.52, 24.2, 484000),
    "ATA": ("I", 0.14, 6.5, 130000),
    "ATG": ("M", 1.00, 23.0, 460000),
    "GTT": ("V", 0.18, 12.0, 240000),
    "GTC": ("V", 0.24, 16.0, 320000),
    "GTA": ("V", 0.10, 6.8, 136000),
    "GTG": ("V", 0.48, 32.0, 640000),
    "TCT": ("S", 0.18, 15.0, 300000),
    "TCC": ("S", 0.22, 18.2, 364000),
    "TCA": ("S", 0.13, 10.8, 216000),
    "TCG": ("S", 0.05, 4.2, 84000),
    "CCT": ("P", 0.26, 15.0, 300000),
    "CCC": ("P", 0.30, 17.5, 350000),
    "CCA": ("P", 0.30, 17.5, 350000),
    "CCG": ("P", 0.14, 8.2, 164000),
    "ACT": ("T", 0.23, 13.5, 270000),
    "ACC": ("T", 0.40, 23.5, 470000),
    "ACA": ("T", 0.25, 14.8, 296000),
    "ACG": ("T", 0.12, 7.0, 140000),
    "GCT": ("A", 0.24, 17.5, 350000),
    "GCC": ("A", 0.42, 30.5, 610000),
    "GCA": ("A", 0.22, 16.0, 320000),
    "GCG": ("A", 0.12, 8.8, 176000),
    "TAT": ("Y", 0.42, 11.5, 230000),
    "TAC": ("Y", 0.58, 15.8, 316000),
    "TAA": ("*", 0.35, 1.0, 20000),
    "TAG": ("*", 0.20, 0.6, 12000),
    "CAT": ("H", 0.40, 10.2, 204000),
    "CAC": ("H", 0.60, 15.2, 304000),
    "CAA": ("Q", 0.25, 12.5, 250000),
    "CAG": ("Q", 0.75, 37.5, 750000),
    "AAT": ("N", 0.44, 17.0, 340000),
    "AAC": ("N", 0.56, 21.6, 432000),
    "AAA": ("K", 0.40, 25.0, 500000),
    "AAG": ("K", 0.60, 37.5, 750000),
    "GAT": ("D", 0.44, 23.0, 460000),
    "GAC": ("D", 0.56, 29.2, 584000),
    "GAA": ("E", 0.42, 32.0, 640000),
    "GAG": ("E", 0.58, 44.0, 880000),
    "TGT": ("C", 0.42, 9.8, 196000),
    "TGC": ("C", 0.58, 13.5, 270000),
    "TGA": ("*", 0.45, 1.4, 28000),
    "TGG": ("W", 1.00, 13.0, 260000),
    "CGT": ("R", 0.08, 5.0, 100000),
    "CGC": ("R", 0.18, 11.2, 224000),
    "CGA": ("R", 0.08, 5.2, 104000),
    "CGG": ("R", 0.18, 11.0, 220000),
    "AGA": ("R", 0.30, 18.5, 370000),
    "AGG": ("R", 0.18, 11.2, 224000),
    "AGT": ("S", 0.14, 11.5, 230000),
    "AGC": ("S", 0.28, 23.0, 460000),
    "GGT": ("G", 0.18, 13.0, 260000),
    "GGC": ("G", 0.36, 26.0, 520000),
    "GGA": ("G", 0.26, 18.8, 376000),
    "GGG": ("G", 0.20, 14.5, 290000),
}

DANIO_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(DANIO_CODON_USAGE)
DANIO_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(DANIO_CODON_USAGE)
