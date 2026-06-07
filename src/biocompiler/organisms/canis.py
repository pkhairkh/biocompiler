"""
Canis familiaris (Dog) Codon Usage Data

Source: Kazusa Codon Usage Database
Canis lupus familiaris, ~19,000 CDSs, coding GC: ~43.5%

The domestic dog is an important model in veterinary medicine,
comparative genomics, and pharmacology.  Its codon usage reflects
a moderately AT-rich mammalian pattern.

Key features:
  - Slightly more AT-rich than human/mouse
  - GC-ending codons still preferred for most amino acids
  - AGA preferred for Arg
  - Important for veterinary biopharmaceuticals

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 9615
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CANIS_CODON_USAGE",
    "CANIS_CODON_ADAPTIVENESS",
    "CANIS_PREFERRED_CODONS",
]

CANIS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.47, 17.5, 262500),
    "TTC": ("F", 0.53, 19.8, 297000),
    "TTA": ("L", 0.09, 8.5, 127500),
    "TTG": ("L", 0.14, 13.2, 198000),
    "CTT": ("L", 0.14, 13.0, 195000),
    "CTC": ("L", 0.19, 17.8, 267000),
    "CTA": ("L", 0.07, 6.5, 97500),
    "CTG": ("L", 0.37, 34.8, 522000),
    "ATT": ("I", 0.37, 16.5, 247500),
    "ATC": ("I", 0.47, 21.0, 315000),
    "ATA": ("I", 0.16, 7.2, 108000),
    "ATG": ("M", 1.00, 21.5, 322500),
    "GTT": ("V", 0.19, 12.0, 180000),
    "GTC": ("V", 0.22, 13.8, 207000),
    "GTA": ("V", 0.11, 7.0, 105000),
    "GTG": ("V", 0.48, 30.0, 450000),
    "TCT": ("S", 0.20, 14.5, 217500),
    "TCC": ("S", 0.21, 15.2, 228000),
    "TCA": ("S", 0.15, 10.8, 162000),
    "TCG": ("S", 0.05, 3.8, 57000),
    "CCT": ("P", 0.28, 15.0, 225000),
    "CCC": ("P", 0.28, 15.0, 225000),
    "CCA": ("P", 0.28, 15.0, 225000),
    "CCG": ("P", 0.16, 8.6, 129000),
    "ACT": ("T", 0.25, 13.2, 198000),
    "ACC": ("T", 0.34, 18.0, 270000),
    "ACA": ("T", 0.26, 13.8, 207000),
    "ACG": ("T", 0.15, 8.0, 120000),
    "GCT": ("A", 0.26, 16.0, 240000),
    "GCC": ("A", 0.38, 23.5, 352500),
    "GCA": ("A", 0.24, 14.8, 222000),
    "GCG": ("A", 0.12, 7.4, 111000),
    "TAT": ("Y", 0.44, 11.5, 172500),
    "TAC": ("Y", 0.56, 14.6, 219000),
    "TAA": ("*", 0.35, 1.0, 15000),
    "TAG": ("*", 0.20, 0.6, 9000),
    "CAT": ("H", 0.42, 10.0, 150000),
    "CAC": ("H", 0.58, 13.8, 207000),
    "CAA": ("Q", 0.28, 13.5, 202500),
    "CAG": ("Q", 0.72, 34.5, 517500),
    "AAT": ("N", 0.47, 17.0, 255000),
    "AAC": ("N", 0.53, 19.2, 288000),
    "AAA": ("K", 0.45, 26.0, 390000),
    "AAG": ("K", 0.55, 31.8, 477000),
    "GAT": ("D", 0.45, 22.0, 330000),
    "GAC": ("D", 0.55, 27.0, 405000),
    "GAA": ("E", 0.43, 30.0, 450000),
    "GAG": ("E", 0.57, 39.8, 597000),
    "TGT": ("C", 0.45, 9.5, 142500),
    "TGC": ("C", 0.55, 11.6, 174000),
    "TGA": ("*", 0.45, 1.3, 19500),
    "TGG": ("W", 1.00, 12.5, 187500),
    "CGT": ("R", 0.08, 5.0, 75000),
    "CGC": ("R", 0.15, 9.2, 138000),
    "CGA": ("R", 0.08, 5.0, 75000),
    "CGG": ("R", 0.18, 11.0, 165000),
    "AGA": ("R", 0.32, 19.5, 292500),
    "AGG": ("R", 0.19, 11.6, 174000),
    "AGT": ("S", 0.16, 11.0, 165000),
    "AGC": ("S", 0.23, 16.8, 252000),
    "GGT": ("G", 0.18, 12.5, 187500),
    "GGC": ("G", 0.32, 22.0, 330000),
    "GGA": ("G", 0.28, 19.2, 288000),
    "GGG": ("G", 0.22, 15.2, 228000),
}

CANIS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(CANIS_CODON_USAGE)
CANIS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(CANIS_CODON_USAGE)
