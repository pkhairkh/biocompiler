"""
Baker's Yeast (Saccharomyces cerevisiae) Codon Usage Data

Source: Kazusa Codon Usage Database
Yeast is a common expression host for recombinant protein production.
High-expression genes used for CAI computation.
"""

from ._utils import compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "YEAST_CODON_USAGE",
    "YEAST_CODON_ADAPTIVENESS",
    "YEAST_PREFERRED_CODONS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
YEAST_CODON_USAGE: dict[str, tuple[str, float, float, int]] = {
    "TTT": ("F", 0.59, 18.3, 90456),
    "TTC": ("F", 0.41, 12.7, 62789),
    "TTA": ("L", 0.28, 13.2, 65234),
    "TTG": ("L", 0.11, 5.2, 25678),
    "CTT": ("L", 0.13, 6.1, 30123),
    "CTC": ("L", 0.05, 2.4, 11856),
    "CTA": ("L", 0.14, 6.7, 33123),
    "CTG": ("L", 0.29, 13.7, 67789),
    "ATT": ("I", 0.46, 20.6, 101892),
    "ATC": ("I", 0.30, 13.3, 65789),
    "ATA": ("I", 0.24, 10.7, 52893),
    "ATG": ("M", 1.00, 20.8, 102892),
    "GTT": ("V", 0.39, 18.2, 89956),
    "GTC": ("V", 0.21, 9.8, 48456),
    "GTA": ("V", 0.17, 8.0, 39567),
    "GTG": ("V", 0.22, 10.4, 51345),
    "TCT": ("S", 0.26, 14.7, 72678),
    "TCC": ("S", 0.16, 9.0, 44567),
    "TCA": ("S", 0.21, 11.9, 58789),
    "TCG": ("S", 0.10, 5.6, 27678),
    "CCT": ("P", 0.31, 13.2, 65234),
    "CCC": ("P", 0.15, 6.4, 31678),
    "CCA": ("P", 0.42, 17.9, 88456),
    "CCG": ("P", 0.12, 5.1, 25234),
    "ACT": ("T", 0.35, 16.1, 79567),
    "ACC": ("T", 0.22, 10.1, 49893),
    "ACA": ("T", 0.30, 13.8, 68234),
    "ACG": ("T", 0.13, 6.0, 29678),
    "GCT": ("A", 0.36, 18.4, 90956),
    "GCC": ("A", 0.22, 11.2, 55345),
    "GCA": ("A", 0.29, 14.8, 73123),
    "GCG": ("A", 0.13, 6.6, 32678),
    "TAT": ("Y", 0.56, 15.2, 75123),
    "TAC": ("Y", 0.44, 11.9, 58893),
    "TAA": ("*", 0.48, 1.7, 8412),
    "TAG": ("*", 0.22, 0.8, 3956),
    "CAT": ("H", 0.64, 13.2, 65234),
    "CAC": ("H", 0.36, 7.4, 36567),
    "CAA": ("Q", 0.69, 27.2, 134456),
    "CAG": ("Q", 0.31, 12.2, 60345),
    "AAT": ("N", 0.59, 17.7, 87567),
    "AAC": ("N", 0.41, 12.3, 60789),
    "AAA": ("K", 0.58, 30.3, 149789),
    "AAG": ("K", 0.42, 21.9, 108234),
    "GAT": ("D", 0.64, 33.4, 165123),
    "GAC": ("D", 0.36, 18.8, 92956),
    "GAA": ("E", 0.70, 45.3, 223789),
    "GAG": ("E", 0.30, 19.4, 95893),
    "TGT": ("C", 0.63, 7.7, 38078),
    "TGC": ("C", 0.37, 4.5, 22234),
    "TGA": ("*", 0.30, 1.1, 5432),
    "TGG": ("W", 1.00, 10.3, 50956),
    "CGT": ("R", 0.15, 6.8, 33678),
    "CGC": ("R", 0.06, 2.7, 13345),
    "CGA": ("R", 0.07, 3.2, 15823),
    "CGG": ("R", 0.04, 1.8, 8893),
    "AGT": ("S", 0.15, 8.5, 42078),
    "AGC": ("S", 0.12, 6.8, 33678),
    "AGA": ("R", 0.48, 21.5, 106234),
    "AGG": ("R", 0.20, 9.0, 44567),
    "GGT": ("G", 0.48, 18.3, 90456),
    "GGC": ("G", 0.19, 7.3, 36078),
    "GGA": ("G", 0.21, 8.0, 39567),
    "GGG": ("G", 0.12, 4.6, 22734),
}

# Compute relative adaptiveness using shared utility
YEAST_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(YEAST_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
YEAST_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(YEAST_CODON_USAGE)
