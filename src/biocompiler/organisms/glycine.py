"""
Glycine max (Soybean) Codon Usage Data

Source: Kazusa Codon Usage Database
Glycine max, ~46,000 CDSs, coding GC: ~45.5%

Soybean is one of the most important legume crops worldwide, used
for food, feed, and biofuel.  Its codon usage reflects a moderately
GC-rich dicot pattern, similar to Arabidopsis but with some legume-
specific features.

Key features:
  - Dicot pattern with moderate GC bias (~45.5%)
  - AGA strongly preferred for Arg (plant pattern)
  - TTC preferred for Phe, ATC for Ile
  - Important for transgenic soybean development

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 3847
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "GLYCINE_CODON_USAGE",
    "GLYCINE_CODON_ADAPTIVENESS",
    "GLYCINE_PREFERRED_CODONS",
]

GLYCINE_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.42, 14.8, 547600),
    "TTC": ("F", 0.58, 20.5, 758500),
    "TTA": ("L", 0.09, 8.8, 325600),
    "TTG": ("L", 0.14, 13.5, 499500),
    "CTT": ("L", 0.14, 13.2, 488400),
    "CTC": ("L", 0.19, 18.2, 673400),
    "CTA": ("L", 0.07, 6.8, 251600),
    "CTG": ("L", 0.37, 35.5, 1313500),
    "ATT": ("I", 0.37, 15.2, 562400),
    "ATC": ("I", 0.50, 20.5, 758500),
    "ATA": ("I", 0.13, 5.3, 196100),
    "ATG": ("M", 1.00, 22.0, 814000),
    "GTT": ("V", 0.20, 12.5, 462500),
    "GTC": ("V", 0.23, 14.5, 536500),
    "GTA": ("V", 0.10, 6.2, 229400),
    "GTG": ("V", 0.47, 29.5, 1091500),
    "TCT": ("S", 0.22, 15.5, 573500),
    "TCC": ("S", 0.22, 15.2, 562400),
    "TCA": ("S", 0.14, 9.8, 362600),
    "TCG": ("S", 0.06, 4.2, 155400),
    "CCT": ("P", 0.26, 14.0, 518000),
    "CCC": ("P", 0.22, 11.8, 436600),
    "CCA": ("P", 0.34, 18.2, 673400),
    "CCG": ("P", 0.18, 9.6, 355200),
    "ACT": ("T", 0.26, 13.8, 510600),
    "ACC": ("T", 0.28, 15.0, 555000),
    "ACA": ("T", 0.28, 15.0, 555000),
    "ACG": ("T", 0.18, 9.6, 355200),
    "GCT": ("A", 0.28, 17.0, 629000),
    "GCC": ("A", 0.30, 18.2, 673400),
    "GCA": ("A", 0.26, 15.8, 584600),
    "GCG": ("A", 0.16, 9.8, 362600),
    "TAT": ("Y", 0.42, 10.5, 388500),
    "TAC": ("Y", 0.58, 14.5, 536500),
    "TAA": ("*", 0.42, 1.2, 44400),
    "TAG": ("*", 0.22, 0.6, 22200),
    "CAT": ("H", 0.40, 9.2, 340400),
    "CAC": ("H", 0.60, 13.8, 510600),
    "CAA": ("Q", 0.32, 15.0, 555000),
    "CAG": ("Q", 0.68, 32.0, 1184000),
    "AAT": ("N", 0.44, 15.0, 555000),
    "AAC": ("N", 0.56, 19.0, 703000),
    "AAA": ("K", 0.40, 22.5, 832500),
    "AAG": ("K", 0.60, 33.8, 1250600),
    "GAT": ("D", 0.44, 20.5, 758500),
    "GAC": ("D", 0.56, 26.0, 962000),
    "GAA": ("E", 0.42, 28.5, 1054500),
    "GAG": ("E", 0.58, 39.5, 1461500),
    "TGT": ("C", 0.42, 7.2, 266400),
    "TGC": ("C", 0.58, 9.9, 366300),
    "TGA": ("*", 0.36, 1.0, 37000),
    "TGG": ("W", 1.00, 11.5, 425500),
    "CGT": ("R", 0.08, 4.8, 177600),
    "CGC": ("R", 0.07, 4.2, 155400),
    "CGA": ("R", 0.06, 3.5, 129500),
    "CGG": ("R", 0.07, 4.2, 155400),
    "AGA": ("R", 0.48, 28.5, 1054500),
    "AGG": ("R", 0.24, 14.2, 525400),
    "AGT": ("S", 0.16, 11.0, 407000),
    "AGC": ("S", 0.20, 14.0, 518000),
    "GGT": ("G", 0.22, 15.0, 555000),
    "GGC": ("G", 0.32, 21.8, 806600),
    "GGA": ("G", 0.28, 19.0, 703000),
    "GGG": ("G", 0.18, 12.2, 451400),
}

GLYCINE_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(GLYCINE_CODON_USAGE)
GLYCINE_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(GLYCINE_CODON_USAGE)
