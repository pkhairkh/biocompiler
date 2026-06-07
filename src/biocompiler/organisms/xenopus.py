"""
Xenopus laevis (African Clawed Frog) Codon Usage Data

Source: Kazusa Codon Usage Database
Xenopus laevis, ~15,000 CDSs, coding GC: ~45.2%

X. laevis is a key model organism for developmental biology, cell
biology, and toxicology.  Its codon usage reflects a moderately
AT-rich vertebrate pattern.

Key features:
  - Similar to other vertebrates but slightly more AT-rich
  - GC-ending codons still preferred for most amino acids
  - AGA preferred for Arg
  - Important for oocyte expression systems

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 8355
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "XENOPUS_CODON_USAGE",
    "XENOPUS_CODON_ADAPTIVENESS",
    "XENOPUS_PREFERRED_CODONS",
]

XENOPUS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.46, 16.0, 240000),
    "TTC": ("F", 0.54, 18.8, 282000),
    "TTA": ("L", 0.10, 9.5, 142500),
    "TTG": ("L", 0.15, 14.2, 213000),
    "CTT": ("L", 0.14, 13.5, 202500),
    "CTC": ("L", 0.19, 18.2, 273000),
    "CTA": ("L", 0.07, 6.8, 102000),
    "CTG": ("L", 0.35, 33.5, 502500),
    "ATT": ("I", 0.38, 17.0, 255000),
    "ATC": ("I", 0.47, 21.0, 315000),
    "ATA": ("I", 0.15, 6.8, 102000),
    "ATG": ("M", 1.00, 21.5, 322500),
    "GTT": ("V", 0.20, 12.5, 187500),
    "GTC": ("V", 0.22, 13.8, 207000),
    "GTA": ("V", 0.12, 7.5, 112500),
    "GTG": ("V", 0.46, 28.8, 432000),
    "TCT": ("S", 0.20, 14.5, 217500),
    "TCC": ("S", 0.20, 14.5, 217500),
    "TCA": ("S", 0.15, 10.8, 162000),
    "TCG": ("S", 0.06, 4.3, 64500),
    "CCT": ("P", 0.28, 14.5, 217500),
    "CCC": ("P", 0.25, 13.0, 195000),
    "CCA": ("P", 0.30, 15.5, 232500),
    "CCG": ("P", 0.17, 8.8, 132000),
    "ACT": ("T", 0.26, 12.8, 192000),
    "ACC": ("T", 0.32, 15.8, 237000),
    "ACA": ("T", 0.26, 12.8, 192000),
    "ACG": ("T", 0.16, 7.8, 117000),
    "GCT": ("A", 0.26, 16.0, 240000),
    "GCC": ("A", 0.36, 22.2, 333000),
    "GCA": ("A", 0.24, 14.8, 222000),
    "GCG": ("A", 0.14, 8.6, 129000),
    "TAT": ("Y", 0.44, 11.5, 172500),
    "TAC": ("Y", 0.56, 14.6, 219000),
    "TAA": ("*", 0.38, 1.0, 15000),
    "TAG": ("*", 0.22, 0.6, 9000),
    "CAT": ("H", 0.42, 10.0, 150000),
    "CAC": ("H", 0.58, 13.8, 207000),
    "CAA": ("Q", 0.30, 14.5, 217500),
    "CAG": ("Q", 0.70, 33.8, 507000),
    "AAT": ("N", 0.46, 17.0, 255000),
    "AAC": ("N", 0.54, 20.0, 300000),
    "AAA": ("K", 0.46, 26.5, 397500),
    "AAG": ("K", 0.54, 31.0, 465000),
    "GAT": ("D", 0.44, 22.0, 330000),
    "GAC": ("D", 0.56, 28.0, 420000),
    "GAA": ("E", 0.44, 30.0, 450000),
    "GAG": ("E", 0.56, 38.0, 570000),
    "TGT": ("C", 0.44, 9.5, 142500),
    "TGC": ("C", 0.56, 12.1, 181500),
    "TGA": ("*", 0.40, 1.2, 18000),
    "TGG": ("W", 1.00, 12.5, 187500),
    "CGT": ("R", 0.08, 5.0, 75000),
    "CGC": ("R", 0.16, 9.8, 147000),
    "CGA": ("R", 0.08, 5.0, 75000),
    "CGG": ("R", 0.18, 11.0, 165000),
    "AGA": ("R", 0.32, 19.5, 292500),
    "AGG": ("R", 0.18, 11.0, 165000),
    "AGT": ("S", 0.16, 11.5, 172500),
    "AGC": ("S", 0.23, 16.5, 247500),
    "GGT": ("G", 0.20, 14.0, 210000),
    "GGC": ("G", 0.32, 22.5, 337500),
    "GGA": ("G", 0.28, 19.5, 292500),
    "GGG": ("G", 0.20, 14.0, 210000),
}

XENOPUS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(XENOPUS_CODON_USAGE)
XENOPUS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(XENOPUS_CODON_USAGE)
