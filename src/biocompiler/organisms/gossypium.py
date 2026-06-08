"""
Gossypium hirsutum (Cotton) Codon Usage Data

Source: Kazusa Codon Usage Database
Gossypium hirsutum, ~37,000 CDSs, coding GC: ~43.8%

Upland cotton is the most widely cultivated cotton species and an
important fiber crop.  Its codon usage reflects a moderately AT-rich
dicot pattern, similar to Arabidopsis.

Key features:
  - Dicot pattern with moderate AT bias (~44% GC)
  - AGA strongly preferred for Arg (plant pattern)
  - TCT preferred for Ser (dicot pattern)
  - Important for transgenic cotton (Bt cotton, herbicide resistance)

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 3635
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "GOSSYPIUM_CODON_USAGE",
    "GOSSYPIUM_CODON_ADAPTIVENESS",
    "GOSSYPIUM_PREFERRED_CODONS",
]

GOSSYPIUM_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.44, 15.5, 456500),
    "TTC": ("F", 0.56, 19.7, 580150),
    "TTA": ("L", 0.10, 9.5, 279750),
    "TTG": ("L", 0.15, 14.2, 417900),
    "CTT": ("L", 0.15, 14.0, 412000),
    "CTC": ("L", 0.18, 16.8, 494760),
    "CTA": ("L", 0.07, 6.5, 191350),
    "CTG": ("L", 0.35, 32.5, 956750),
    "ATT": ("I", 0.38, 16.0, 471200),
    "ATC": ("I", 0.49, 20.8, 612560),
    "ATA": ("I", 0.13, 5.5, 161900),
    "ATG": ("M", 1.00, 22.0, 647600),
    "GTT": ("V", 0.20, 12.8, 376960),
    "GTC": ("V", 0.22, 14.0, 412000),
    "GTA": ("V", 0.10, 6.5, 191350),
    "GTG": ("V", 0.48, 30.5, 898320),
    "TCT": ("S", 0.24, 17.0, 500480),
    "TCC": ("S", 0.20, 14.2, 417900),
    "TCA": ("S", 0.15, 10.5, 309150),
    "TCG": ("S", 0.06, 4.2, 123660),
    "CCT": ("P", 0.28, 14.5, 426770),
    "CCC": ("P", 0.20, 10.5, 309150),
    "CCA": ("P", 0.35, 18.2, 535960),
    "CCG": ("P", 0.17, 8.8, 259040),
    "ACT": ("T", 0.28, 14.2, 417900),
    "ACC": ("T", 0.28, 14.0, 412000),
    "ACA": ("T", 0.28, 14.2, 417900),
    "ACG": ("T", 0.16, 8.2, 241490),
    "GCT": ("A", 0.30, 16.5, 485700),
    "GCC": ("A", 0.28, 15.5, 456500),
    "GCA": ("A", 0.26, 14.5, 426770),
    "GCG": ("A", 0.16, 8.8, 259040),
    "TAT": ("Y", 0.44, 11.0, 323800),
    "TAC": ("Y", 0.56, 14.0, 412000),
    "TAA": ("*", 0.42, 1.2, 35340),
    "TAG": ("*", 0.22, 0.6, 17670),
    "CAT": ("H", 0.42, 9.5, 279750),
    "CAC": ("H", 0.58, 13.1, 385650),
    "CAA": ("Q", 0.35, 15.5, 456500),
    "CAG": ("Q", 0.65, 28.8, 847680),
    "AAT": ("N", 0.45, 16.0, 471200),
    "AAC": ("N", 0.55, 19.5, 574320),
    "AAA": ("K", 0.42, 24.0, 706800),
    "AAG": ("K", 0.58, 33.2, 977480),
    "GAT": ("D", 0.44, 20.5, 603760),
    "GAC": ("D", 0.56, 26.0, 765200),
    "GAA": ("E", 0.42, 29.0, 853540),
    "GAG": ("E", 0.58, 40.0, 1178000),
    "TGT": ("C", 0.44, 7.8, 229690),
    "TGC": ("C", 0.56, 9.9, 291560),
    "TGA": ("*", 0.36, 1.0, 29450),
    "TGG": ("W", 1.00, 11.2, 329720),
    "CGT": ("R", 0.08, 5.0, 147250),
    "CGC": ("R", 0.07, 4.2, 123660),
    "CGA": ("R", 0.06, 3.8, 111880),
    "CGG": ("R", 0.07, 4.2, 123660),
    "AGA": ("R", 0.48, 28.2, 830100),
    "AGG": ("R", 0.24, 14.0, 412000),
    "AGT": ("S", 0.16, 11.2, 329720),
    "AGC": ("S", 0.19, 13.5, 397560),
    "GGT": ("G", 0.22, 14.5, 426770),
    "GGC": ("G", 0.30, 19.8, 583080),
    "GGA": ("G", 0.28, 18.5, 544700),
    "GGG": ("G", 0.20, 13.2, 388560),
}

GOSSYPIUM_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(GOSSYPIUM_CODON_USAGE)
GOSSYPIUM_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(GOSSYPIUM_CODON_USAGE)
