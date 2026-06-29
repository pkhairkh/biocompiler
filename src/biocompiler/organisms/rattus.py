"""
Rattus norvegicus (Rat) Codon Usage Data

Source: Kazusa Codon Usage Database
Rattus norvegicus, ~25,000 CDSs, coding GC: ~51.5%

R. norvegicus (rat) is a key model organism in pharmacology,
toxicology, and biomedical research.  Its codon usage is very
similar to mouse (both murine rodents) with GC-ending codons
preferred.

Key features:
  - Very similar to Mus musculus (murine pattern)
  - GC-ending codons preferred
  - AGA preferred for Arg
  - CTG preferred for Leu

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 10116
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "RATTUS_CODON_USAGE",
    "RATTUS_CODON_ADAPTIVENESS",
    "RATTUS_PREFERRED_CODONS",
]

RATTUS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.44, 16.5, 330000),
    "TTC": ("F", 0.56, 21.0, 420000),
    "TTA": ("L", 0.07, 7.0, 140000),
    "TTG": ("L", 0.12, 12.2, 244000),
    "CTT": ("L", 0.12, 12.0, 240000),
    "CTC": ("L", 0.21, 20.5, 410000),
    "CTA": ("L", 0.07, 6.8, 136000),
    "CTG": ("L", 0.41, 40.0, 800000),
    "ATT": ("I", 0.35, 15.5, 310000),
    "ATC": ("I", 0.49, 21.5, 430000),
    "ATA": ("I", 0.16, 7.0, 140000),
    "ATG": ("M", 1.00, 21.5, 430000),
    "GTT": ("V", 0.18, 11.2, 224000),
    "GTC": ("V", 0.24, 14.8, 296000),
    "GTA": ("V", 0.10, 6.5, 130000),
    "GTG": ("V", 0.48, 29.5, 590000),
    "TCT": ("S", 0.18, 14.2, 284000),
    "TCC": ("S", 0.22, 17.5, 350000),
    "TCA": ("S", 0.14, 11.0, 220000),
    "TCG": ("S", 0.05, 4.0, 80000),
    "CCT": ("P", 0.27, 16.0, 320000),
    "CCC": ("P", 0.33, 19.5, 390000),
    "CCA": ("P", 0.27, 16.0, 320000),
    "CCG": ("P", 0.13, 7.8, 156000),
    "ACT": ("T", 0.24, 12.8, 256000),
    "ACC": ("T", 0.38, 20.0, 400000),
    "ACA": ("T", 0.27, 14.2, 284000),
    "ACG": ("T", 0.11, 5.8, 116000),
    "GCT": ("A", 0.25, 17.5, 350000),
    "GCC": ("A", 0.42, 29.5, 590000),
    "GCA": ("A", 0.22, 15.5, 310000),
    "GCG": ("A", 0.11, 7.8, 156000),
    "TAT": ("Y", 0.43, 11.8, 236000),
    "TAC": ("Y", 0.57, 15.6, 312000),
    "TAA": ("*", 0.30, 1.0, 20000),
    "TAG": ("*", 0.20, 0.7, 14000),
    "CAT": ("H", 0.41, 10.2, 204000),
    "CAC": ("H", 0.59, 14.8, 296000),
    "CAA": ("Q", 0.25, 12.2, 244000),
    "CAG": ("Q", 0.75, 36.5, 730000),
    "AAT": ("N", 0.45, 16.8, 336000),
    "AAC": ("N", 0.55, 20.5, 410000),
    "AAA": ("K", 0.42, 24.5, 490000),
    "AAG": ("K", 0.58, 33.8, 676000),
    "GAT": ("D", 0.44, 21.5, 430000),
    "GAC": ("D", 0.56, 27.2, 544000),
    "GAA": ("E", 0.41, 29.0, 580000),
    "GAG": ("E", 0.59, 41.8, 836000),
    "TGT": ("C", 0.44, 10.0, 200000),
    "TGC": ("C", 0.56, 12.8, 256000),
    "TGA": ("*", 0.50, 1.5, 30000),
    "TGG": ("W", 1.00, 12.8, 256000),
    "CGT": ("R", 0.09, 5.2, 104000),
    "CGC": ("R", 0.18, 11.0, 220000),
    "CGA": ("R", 0.09, 5.6, 112000),
    "CGG": ("R", 0.20, 12.2, 244000),
    "AGA": ("R", 0.26, 16.0, 320000),
    "AGG": ("R", 0.18, 11.0, 220000),
    "AGT": ("S", 0.15, 11.2, 224000),
    "AGC": ("S", 0.26, 20.0, 400000),
    "GGT": ("G", 0.16, 11.0, 220000),
    "GGC": ("G", 0.34, 23.5, 470000),
    "GGA": ("G", 0.26, 17.8, 356000),
    "GGG": ("G", 0.24, 16.5, 330000),
}

RATTUS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(RATTUS_CODON_USAGE)
RATTUS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(RATTUS_CODON_USAGE)
