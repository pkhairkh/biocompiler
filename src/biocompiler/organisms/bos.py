"""
Bos taurus (Cow) Codon Usage Data

Source: Kazusa Codon Usage Database
Bos taurus, ~22,000 CDSs, coding GC: ~52.2%

The cow is an important agricultural species and model for
veterinary biopharmaceuticals.  Its codon usage reflects a
moderately GC-rich mammalian pattern similar to human.

Key features:
  - GC-ending codons preferred (mammalian pattern)
  - Similar to human codon usage
  - AGA preferred for Arg
  - CTG preferred for Leu

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 9913
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "BOS_CODON_USAGE",
    "BOS_CODON_ADAPTIVENESS",
    "BOS_PREFERRED_CODONS",
]

BOS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.44, 16.5, 363000),
    "TTC": ("F", 0.56, 21.0, 462000),
    "TTA": ("L", 0.07, 7.2, 158400),
    "TTG": ("L", 0.12, 12.0, 264000),
    "CTT": ("L", 0.12, 12.2, 268400),
    "CTC": ("L", 0.21, 20.5, 451000),
    "CTA": ("L", 0.07, 6.8, 149600),
    "CTG": ("L", 0.41, 40.0, 880000),
    "ATT": ("I", 0.35, 15.8, 347600),
    "ATC": ("I", 0.49, 22.0, 484000),
    "ATA": ("I", 0.16, 7.2, 158400),
    "ATG": ("M", 1.00, 22.0, 484000),
    "GTT": ("V", 0.18, 11.2, 246400),
    "GTC": ("V", 0.24, 14.8, 325600),
    "GTA": ("V", 0.10, 6.2, 136400),
    "GTG": ("V", 0.48, 29.8, 655600),
    "TCT": ("S", 0.18, 14.5, 319000),
    "TCC": ("S", 0.22, 17.5, 385000),
    "TCA": ("S", 0.14, 11.2, 246400),
    "TCG": ("S", 0.05, 4.0, 88000),
    "CCT": ("P", 0.28, 16.5, 363000),
    "CCC": ("P", 0.32, 19.0, 418000),
    "CCA": ("P", 0.27, 16.0, 352000),
    "CCG": ("P", 0.13, 7.8, 171600),
    "ACT": ("T", 0.24, 13.0, 286000),
    "ACC": ("T", 0.38, 20.5, 451000),
    "ACA": ("T", 0.27, 14.5, 319000),
    "ACG": ("T", 0.11, 6.0, 132000),
    "GCT": ("A", 0.25, 17.5, 385000),
    "GCC": ("A", 0.42, 29.5, 649000),
    "GCA": ("A", 0.22, 15.5, 341000),
    "GCG": ("A", 0.11, 7.8, 171600),
    "TAT": ("Y", 0.43, 11.8, 259600),
    "TAC": ("Y", 0.57, 15.6, 343200),
    "TAA": ("*", 0.30, 1.0, 22000),
    "TAG": ("*", 0.20, 0.7, 15400),
    "CAT": ("H", 0.41, 10.4, 228800),
    "CAC": ("H", 0.59, 15.0, 330000),
    "CAA": ("Q", 0.25, 12.5, 275000),
    "CAG": ("Q", 0.75, 37.5, 825000),
    "AAT": ("N", 0.46, 17.0, 374000),
    "AAC": ("N", 0.54, 20.0, 440000),
    "AAA": ("K", 0.42, 24.5, 539000),
    "AAG": ("K", 0.58, 33.8, 743600),
    "GAT": ("D", 0.44, 21.5, 473000),
    "GAC": ("D", 0.56, 27.5, 605000),
    "GAA": ("E", 0.41, 29.5, 649000),
    "GAG": ("E", 0.59, 42.5, 935000),
    "TGT": ("C", 0.44, 10.0, 220000),
    "TGC": ("C", 0.56, 12.8, 281600),
    "TGA": ("*", 0.50, 1.5, 33000),
    "TGG": ("W", 1.00, 13.0, 286000),
    "CGT": ("R", 0.08, 5.0, 110000),
    "CGC": ("R", 0.18, 11.0, 242000),
    "CGA": ("R", 0.09, 5.6, 123200),
    "CGG": ("R", 0.20, 12.2, 268400),
    "AGA": ("R", 0.28, 17.0, 374000),
    "AGG": ("R", 0.17, 10.4, 228800),
    "AGT": ("S", 0.15, 11.2, 246400),
    "AGC": ("S", 0.26, 20.0, 440000),
    "GGT": ("G", 0.16, 11.2, 246400),
    "GGC": ("G", 0.35, 24.5, 539000),
    "GGA": ("G", 0.26, 18.2, 400400),
    "GGG": ("G", 0.23, 16.0, 352000),
}

BOS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(BOS_CODON_USAGE)
BOS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(BOS_CODON_USAGE)
