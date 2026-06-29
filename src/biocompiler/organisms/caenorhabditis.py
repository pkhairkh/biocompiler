"""
Caenorhabditis elegans Codon Usage Data

Source: Kazusa Codon Usage Database
Caenorhabditis elegans, ~20,000 CDSs, coding GC: ~43.0%

C. elegans is a nematode model organism widely used in genetics,
developmental biology, and neuroscience.  Its codon usage reflects
a moderately AT-rich genome with strong translational selection
in highly expressed genes.

Key features:
  - A/T-ending codons generally preferred (GC ~43%)
  - TTG preferred for Leu (nematode pattern)
  - AGA preferred for Arg
  - AAA preferred for Lys (AT-rich)

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 6239
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CAENORHABDITIS_CODON_USAGE",
    "CAENORHABDITIS_CODON_ADAPTIVENESS",
    "CAENORHABDITIS_PREFERRED_CODONS",
]

CAENORHABDITIS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.52, 16.2, 243000),
    "TTC": ("F", 0.48, 14.9, 223500),
    "TTA": ("L", 0.14, 10.5, 157500),
    "TTG": ("L", 0.24, 17.8, 267000),
    "CTT": ("L", 0.14, 10.4, 156000),
    "CTC": ("L", 0.12, 9.0, 135000),
    "CTA": ("L", 0.06, 4.5, 67500),
    "CTG": ("L", 0.30, 22.4, 336000),
    "ATT": ("I", 0.44, 17.5, 262500),
    "ATC": ("I", 0.42, 16.7, 250500),
    "ATA": ("I", 0.14, 5.6, 84000),
    "ATG": ("M", 1.00, 20.5, 307500),
    "GTT": ("V", 0.25, 13.8, 207000),
    "GTC": ("V", 0.20, 11.0, 165000),
    "GTA": ("V", 0.15, 8.3, 124500),
    "GTG": ("V", 0.40, 22.0, 330000),
    "TCT": ("S", 0.24, 16.5, 247500),
    "TCC": ("S", 0.18, 12.4, 186000),
    "TCA": ("S", 0.16, 11.0, 165000),
    "TCG": ("S", 0.07, 4.8, 72000),
    "CCT": ("P", 0.28, 12.5, 187500),
    "CCC": ("P", 0.18, 8.0, 120000),
    "CCA": ("P", 0.36, 16.0, 240000),
    "CCG": ("P", 0.18, 8.0, 120000),
    "ACT": ("T", 0.30, 13.8, 207000),
    "ACC": ("T", 0.25, 11.5, 172500),
    "ACA": ("T", 0.28, 12.8, 192000),
    "ACG": ("T", 0.17, 7.8, 117000),
    "GCT": ("A", 0.30, 17.5, 262500),
    "GCC": ("A", 0.24, 14.0, 210000),
    "GCA": ("A", 0.28, 16.2, 243000),
    "GCG": ("A", 0.18, 10.5, 157500),
    "TAT": ("Y", 0.48, 11.0, 165000),
    "TAC": ("Y", 0.52, 11.9, 178500),
    "TAA": ("*", 0.45, 1.2, 18000),
    "TAG": ("*", 0.20, 0.5, 7500),
    "CAT": ("H", 0.48, 9.5, 142500),
    "CAC": ("H", 0.52, 10.3, 154500),
    "CAA": ("Q", 0.38, 17.5, 262500),
    "CAG": ("Q", 0.62, 28.5, 427500),
    "AAT": ("N", 0.48, 16.5, 247500),
    "AAC": ("N", 0.52, 17.8, 267000),
    "AAA": ("K", 0.52, 28.0, 420000),
    "AAG": ("K", 0.48, 25.8, 387000),
    "GAT": ("D", 0.48, 24.0, 360000),
    "GAC": ("D", 0.52, 26.0, 390000),
    "GAA": ("E", 0.52, 32.0, 480000),
    "GAG": ("E", 0.48, 29.5, 442500),
    "TGT": ("C", 0.48, 8.5, 127500),
    "TGC": ("C", 0.52, 9.2, 138000),
    "TGA": ("*", 0.35, 0.9, 13500),
    "TGG": ("W", 1.00, 11.5, 172500),
    "CGT": ("R", 0.10, 5.8, 87000),
    "CGC": ("R", 0.08, 4.6, 69000),
    "CGA": ("R", 0.06, 3.5, 52500),
    "CGG": ("R", 0.08, 4.6, 69000),
    "AGA": ("R", 0.42, 24.2, 363000),
    "AGG": ("R", 0.26, 15.0, 225000),
    "AGT": ("S", 0.18, 12.2, 183000),
    "AGC": ("S", 0.17, 11.6, 174000),
    "GGT": ("G", 0.24, 14.0, 210000),
    "GGC": ("G", 0.24, 14.0, 210000),
    "GGA": ("G", 0.32, 18.5, 277500),
    "GGG": ("G", 0.20, 11.6, 174000),
}

CAENORHABDITIS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(CAENORHABDITIS_CODON_USAGE)
CAENORHABDITIS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(CAENORHABDITIS_CODON_USAGE)
