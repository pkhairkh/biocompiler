"""
Gallus gallus (Chicken) Codon Usage Data

Source: Kazusa Codon Usage Database
Gallus gallus, ~18,000 CDSs, coding GC: ~47.5%

The chicken is an important agricultural species and model organism
for developmental biology.  Its codon usage reflects a moderately
GC-rich avian pattern.

Key features:
  - GC-ending codons preferred for most amino acids
  - Slightly more AT-rich than mammalian pattern
  - AGA preferred for Arg (consistent with vertebrate pattern)
  - Important for poultry vaccine and protein production

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 9031
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "GALLUS_CODON_USAGE",
    "GALLUS_CODON_ADAPTIVENESS",
    "GALLUS_PREFERRED_CODONS",
]

GALLUS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.45, 16.8, 252000),
    "TTC": ("F", 0.55, 20.5, 307500),
    "TTA": ("L", 0.08, 7.8, 117000),
    "TTG": ("L", 0.14, 13.5, 202500),
    "CTT": ("L", 0.13, 12.5, 187500),
    "CTC": ("L", 0.20, 19.2, 288000),
    "CTA": ("L", 0.07, 6.8, 102000),
    "CTG": ("L", 0.38, 36.5, 547500),
    "ATT": ("I", 0.36, 16.2, 243000),
    "ATC": ("I", 0.48, 21.5, 322500),
    "ATA": ("I", 0.16, 7.0, 105000),
    "ATG": ("M", 1.00, 21.5, 322500),
    "GTT": ("V", 0.18, 11.8, 177000),
    "GTC": ("V", 0.22, 14.2, 213000),
    "GTA": ("V", 0.11, 7.0, 105000),
    "GTG": ("V", 0.49, 31.5, 472500),
    "TCT": ("S", 0.19, 14.2, 213000),
    "TCC": ("S", 0.22, 16.5, 247500),
    "TCA": ("S", 0.14, 10.5, 157500),
    "TCG": ("S", 0.06, 4.5, 67500),
    "CCT": ("P", 0.27, 15.5, 232500),
    "CCC": ("P", 0.28, 16.0, 240000),
    "CCA": ("P", 0.29, 16.8, 252000),
    "CCG": ("P", 0.16, 9.2, 138000),
    "ACT": ("T", 0.24, 13.0, 195000),
    "ACC": ("T", 0.36, 19.5, 292500),
    "ACA": ("T", 0.26, 14.0, 210000),
    "ACG": ("T", 0.14, 7.5, 112500),
    "GCT": ("A", 0.26, 17.0, 255000),
    "GCC": ("A", 0.38, 25.0, 375000),
    "GCA": ("A", 0.24, 15.8, 237000),
    "GCG": ("A", 0.12, 7.8, 117000),
    "TAT": ("Y", 0.43, 11.2, 168000),
    "TAC": ("Y", 0.57, 14.8, 222000),
    "TAA": ("*", 0.38, 1.0, 15000),
    "TAG": ("*", 0.20, 0.5, 7500),
    "CAT": ("H", 0.42, 10.0, 150000),
    "CAC": ("H", 0.58, 13.8, 207000),
    "CAA": ("Q", 0.27, 13.0, 195000),
    "CAG": ("Q", 0.73, 35.0, 525000),
    "AAT": ("N", 0.46, 16.5, 247500),
    "AAC": ("N", 0.54, 19.4, 291000),
    "AAA": ("K", 0.44, 25.5, 382500),
    "AAG": ("K", 0.56, 32.4, 486000),
    "GAT": ("D", 0.44, 21.8, 327000),
    "GAC": ("D", 0.56, 27.8, 417000),
    "GAA": ("E", 0.42, 30.5, 457500),
    "GAG": ("E", 0.58, 42.0, 630000),
    "TGT": ("C", 0.44, 9.8, 147000),
    "TGC": ("C", 0.56, 12.5, 187500),
    "TGA": ("*", 0.42, 1.2, 18000),
    "TGG": ("W", 1.00, 12.8, 192000),
    "CGT": ("R", 0.08, 5.2, 78000),
    "CGC": ("R", 0.16, 10.2, 153000),
    "CGA": ("R", 0.08, 5.0, 75000),
    "CGG": ("R", 0.19, 12.0, 180000),
    "AGA": ("R", 0.32, 20.2, 303000),
    "AGG": ("R", 0.17, 10.8, 162000),
    "AGT": ("S", 0.16, 11.5, 172500),
    "AGC": ("S", 0.23, 17.5, 262500),
    "GGT": ("G", 0.18, 12.8, 192000),
    "GGC": ("G", 0.34, 24.0, 360000),
    "GGA": ("G", 0.27, 19.0, 285000),
    "GGG": ("G", 0.21, 14.8, 222000),
}

GALLUS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(GALLUS_CODON_USAGE)
GALLUS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(GALLUS_CODON_USAGE)
