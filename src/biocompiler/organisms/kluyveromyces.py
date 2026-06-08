"""
Kluyveromyces lactis Codon Usage Data

Source: Kazusa Codon Usage Database
Kluyveromyces lactis NRRL Y-1140, ~5,300 CDSs, coding GC: ~38.6%

K. lactis is a non-conventional yeast used for protein expression,
particularly for secreted proteins and lactose-utilizing applications.
Its codon usage is similar to S. cerevisiae (AT-rich yeast pattern)
but with some important differences.

Key features:
  - A/T-ending codons generally preferred (yeast pattern)
  - AGA strongly preferred for Arg (like S. cerevisiae)
  - Slightly less extreme AT bias than S. cerevisiae
  - Important alternative to Pichia for some protein expressions

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 4959
  - Bolotin-Fukuhara et al. (2000) FEMS Microbiol Rev 24:1-12
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "KLUYVEROMYCES_CODON_USAGE",
    "KLUYVEROMYCES_CODON_ADAPTIVENESS",
    "KLUYVEROMYCES_PREFERRED_CODONS",
]

KLUYVEROMYCES_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.60, 18.5, 85000),
    "TTC": ("F", 0.40, 12.4, 57000),
    "TTA": ("L", 0.22, 10.0, 46000),
    "TTG": ("L", 0.35, 16.0, 73500),
    "CTT": ("L", 0.12, 5.5, 25250),
    "CTC": ("L", 0.07, 3.2, 14700),
    "CTA": ("L", 0.05, 2.3, 10570),
    "CTG": ("L", 0.19, 8.7, 40000),
    "ATT": ("I", 0.52, 19.8, 91000),
    "ATC": ("I", 0.42, 16.0, 73500),
    "ATA": ("I", 0.06, 2.3, 10570),
    "ATG": ("M", 1.00, 20.2, 92800),
    "GTT": ("V", 0.40, 17.5, 80500),
    "GTC": ("V", 0.20, 8.8, 40400),
    "GTA": ("V", 0.22, 9.6, 44100),
    "GTG": ("V", 0.18, 7.9, 36300),
    "TCT": ("S", 0.32, 16.5, 75800),
    "TCC": ("S", 0.20, 10.3, 47300),
    "TCA": ("S", 0.18, 9.2, 42300),
    "TCG": ("S", 0.06, 3.1, 14250),
    "CCT": ("P", 0.28, 11.5, 52800),
    "CCC": ("P", 0.15, 6.2, 28500),
    "CCA": ("P", 0.42, 17.2, 79000),
    "CCG": ("P", 0.15, 6.2, 28500),
    "ACT": ("T", 0.40, 16.0, 73500),
    "ACC": ("T", 0.25, 10.0, 46000),
    "ACA": ("T", 0.25, 10.0, 46000),
    "ACG": ("T", 0.10, 4.0, 18400),
    "GCT": ("A", 0.38, 16.5, 75800),
    "GCC": ("A", 0.22, 9.6, 44100),
    "GCA": ("A", 0.28, 12.2, 56000),
    "GCG": ("A", 0.12, 5.2, 23900),
    "TAT": ("Y", 0.52, 12.5, 57500),
    "TAC": ("Y", 0.48, 11.5, 52800),
    "TAA": ("*", 0.48, 1.5, 6900),
    "TAG": ("*", 0.15, 0.5, 2300),
    "CAT": ("H", 0.58, 10.5, 48200),
    "CAC": ("H", 0.42, 7.6, 34900),
    "CAA": ("Q", 0.68, 22.0, 101000),
    "CAG": ("Q", 0.32, 10.4, 47800),
    "AAT": ("N", 0.55, 16.0, 73500),
    "AAC": ("N", 0.45, 13.1, 60200),
    "AAA": ("K", 0.58, 28.5, 131000),
    "AAG": ("K", 0.42, 20.6, 94700),
    "GAT": ("D", 0.58, 25.0, 115000),
    "GAC": ("D", 0.42, 18.1, 83200),
    "GAA": ("E", 0.65, 38.0, 175000),
    "GAG": ("E", 0.35, 20.4, 93800),
    "TGT": ("C", 0.55, 7.5, 34500),
    "TGC": ("C", 0.45, 6.1, 28000),
    "TGA": ("*", 0.37, 1.0, 4600),
    "TGG": ("W", 1.00, 9.8, 45000),
    "CGT": ("R", 0.10, 4.5, 20700),
    "CGC": ("R", 0.04, 1.8, 8280),
    "CGA": ("R", 0.03, 1.4, 6430),
    "CGG": ("R", 0.04, 1.8, 8280),
    "AGA": ("R", 0.58, 26.2, 120000),
    "AGG": ("R", 0.21, 9.5, 43700),
    "AGT": ("S", 0.14, 7.2, 33100),
    "AGC": ("S", 0.10, 5.1, 23400),
    "GGT": ("G", 0.40, 14.0, 64300),
    "GGC": ("G", 0.18, 6.3, 28950),
    "GGA": ("G", 0.28, 9.8, 45000),
    "GGG": ("G", 0.14, 4.9, 22500),
}

KLUYVEROMYCES_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(KLUYVEROMYCES_CODON_USAGE)
KLUYVEROMYCES_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(KLUYVEROMYCES_CODON_USAGE)
