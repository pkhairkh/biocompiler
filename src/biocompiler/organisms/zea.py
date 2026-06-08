"""
Zea mays (Maize) Codon Usage Data

Source: Kazusa Codon Usage Database
Zea mays, ~30,000 CDSs, coding GC: ~55.8%

Maize is one of the world's most important cereal crops and a model
organism for monocot molecular biology.  Its codon usage reflects
a very GC-rich monocot pattern, similar to rice but with even
stronger GC bias.

Key features:
  - Very strong GC-ending codon preference (GC ~56%)
  - GCC strongly preferred for Ala
  - CCG common for Pro (monocot pattern)
  - GCG significant for Ala
  - Important for transgenic crop development

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 4577
  - Liu et al. (2004) Plant Physiol 135:1014-1024
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "ZEA_CODON_USAGE",
    "ZEA_CODON_ADAPTIVENESS",
    "ZEA_PREFERRED_CODONS",
]

ZEA_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.34, 12.5, 375000),
    "TTC": ("F", 0.66, 24.3, 729000),
    "TTA": ("L", 0.05, 5.0, 150000),
    "TTG": ("L", 0.11, 10.8, 324000),
    "CTT": ("L", 0.11, 10.5, 315000),
    "CTC": ("L", 0.24, 23.0, 690000),
    "CTA": ("L", 0.05, 5.2, 156000),
    "CTG": ("L", 0.44, 42.5, 1275000),
    "ATT": ("I", 0.29, 12.8, 384000),
    "ATC": ("I", 0.60, 26.5, 795000),
    "ATA": ("I", 0.11, 4.8, 144000),
    "ATG": ("M", 1.00, 23.5, 705000),
    "GTT": ("V", 0.16, 10.5, 315000),
    "GTC": ("V", 0.25, 16.2, 486000),
    "GTA": ("V", 0.07, 4.8, 144000),
    "GTG": ("V", 0.52, 33.8, 1014000),
    "TCT": ("S", 0.18, 13.5, 405000),
    "TCC": ("S", 0.24, 18.0, 540000),
    "TCA": ("S", 0.10, 7.5, 225000),
    "TCG": ("S", 0.12, 9.0, 270000),
    "CCT": ("P", 0.18, 11.5, 345000),
    "CCC": ("P", 0.20, 12.8, 384000),
    "CCA": ("P", 0.28, 17.8, 534000),
    "CCG": ("P", 0.34, 21.8, 654000),
    "ACT": ("T", 0.20, 12.0, 360000),
    "ACC": ("T", 0.42, 25.0, 750000),
    "ACA": ("T", 0.18, 10.8, 324000),
    "ACG": ("T", 0.20, 12.0, 360000),
    "GCT": ("A", 0.20, 15.0, 450000),
    "GCC": ("A", 0.38, 28.5, 855000),
    "GCA": ("A", 0.18, 13.5, 405000),
    "GCG": ("A", 0.24, 18.0, 540000),
    "TAT": ("Y", 0.36, 10.0, 300000),
    "TAC": ("Y", 0.64, 17.8, 534000),
    "TAA": ("*", 0.42, 1.0, 30000),
    "TAG": ("*", 0.22, 0.5, 15000),
    "CAT": ("H", 0.34, 8.2, 246000),
    "CAC": ("H", 0.66, 15.8, 474000),
    "CAA": ("Q", 0.25, 14.8, 444000),
    "CAG": ("Q", 0.75, 44.2, 1326000),
    "AAT": ("N", 0.35, 13.0, 390000),
    "AAC": ("N", 0.65, 24.0, 720000),
    "AAA": ("K", 0.32, 19.8, 594000),
    "AAG": ("K", 0.68, 42.0, 1260000),
    "GAT": ("D", 0.36, 18.5, 555000),
    "GAC": ("D", 0.64, 33.0, 990000),
    "GAA": ("E", 0.35, 26.8, 804000),
    "GAG": ("E", 0.65, 49.8, 1494000),
    "TGT": ("C", 0.34, 6.5, 195000),
    "TGC": ("C", 0.66, 12.6, 378000),
    "TGA": ("*", 0.36, 0.8, 24000),
    "TGG": ("W", 1.00, 11.8, 354000),
    "CGT": ("R", 0.10, 6.2, 186000),
    "CGC": ("R", 0.15, 9.2, 276000),
    "CGA": ("R", 0.05, 3.0, 90000),
    "CGG": ("R", 0.10, 6.2, 186000),
    "AGA": ("R", 0.42, 26.0, 780000),
    "AGG": ("R", 0.18, 11.2, 336000),
    "AGT": ("S", 0.12, 9.0, 270000),
    "AGC": ("S", 0.24, 18.0, 540000),
    "GGT": ("G", 0.18, 14.0, 420000),
    "GGC": ("G", 0.38, 29.5, 885000),
    "GGA": ("G", 0.25, 19.5, 585000),
    "GGG": ("G", 0.19, 14.8, 444000),
}

ZEA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(ZEA_CODON_USAGE)
ZEA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(ZEA_CODON_USAGE)
