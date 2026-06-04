"""
Chinese Hamster Ovary (CHO-K1) Codon Usage Data

Source: Kazusa Codon Usage Database
CHO cells are the most commonly used mammalian host for
biopharmaceutical protein production.
"""

from ._utils import compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CHO_CODON_USAGE",
    "CHO_CODON_ADAPTIVENESS",
    "CHO_PREFERRED_CODONS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
CHO_CODON_USAGE: dict[str, tuple[str, float, float, int]] = {
    "TTT": ("F", 0.45, 17.2, 234567),
    "TTC": ("F", 0.55, 20.9, 285123),
    "TTA": ("L", 0.07, 7.4, 100892),
    "TTG": ("L", 0.13, 12.8, 174523),
    "CTT": ("L", 0.13, 12.4, 169123),
    "CTC": ("L", 0.21, 19.8, 269856),
    "CTA": ("L", 0.07, 7.1, 96789),
    "CTG": ("L", 0.39, 38.5, 525123),
    "ATT": ("I", 0.35, 15.6, 212456),
    "ATC": ("I", 0.49, 21.8, 297123),
    "ATA": ("I", 0.16, 7.2, 98123),
    "ATG": ("M", 1.00, 21.8, 297456),
    "GTT": ("V", 0.18, 11.3, 154123),
    "GTC": ("V", 0.24, 14.7, 200456),
    "GTA": ("V", 0.11, 7.0, 95456),
    "GTG": ("V", 0.47, 28.9, 394123),
    "TCT": ("S", 0.18, 14.3, 195123),
    "TCC": ("S", 0.23, 18.5, 252456),
    "TCA": ("S", 0.15, 12.0, 163789),
    "TCG": ("S", 0.05, 4.2, 57234),
    "CCT": ("P", 0.28, 17.2, 234567),
    "CCC": ("P", 0.33, 20.5, 279456),
    "CCA": ("P", 0.27, 16.5, 225123),
    "CCG": ("P", 0.12, 7.2, 98123),
    "ACT": ("T", 0.24, 12.8, 174456),
    "ACC": ("T", 0.38, 19.8, 269789),
    "ACA": ("T", 0.27, 14.2, 193456),
    "ACG": ("T", 0.11, 5.8, 79123),
    "GCT": ("A", 0.25, 17.8, 242789),
    "GCC": ("A", 0.42, 29.5, 402123),
    "GCA": ("A", 0.22, 15.4, 210123),
    "GCG": ("A", 0.11, 7.5, 102456),
    "TAT": ("Y", 0.44, 12.0, 163456),
    "TAC": ("Y", 0.56, 15.4, 210123),
    "TAA": ("*", 0.29, 1.0, 13645),
    "TAG": ("*", 0.24, 0.8, 10918),
    "CAT": ("H", 0.41, 10.6, 144456),
    "CAC": ("H", 0.59, 15.2, 207456),
    "CAA": ("Q", 0.25, 12.1, 165123),
    "CAG": ("Q", 0.75, 35.8, 488123),
    "AAT": ("N", 0.46, 17.2, 234567),
    "AAC": ("N", 0.54, 20.1, 274123),
    "AAA": ("K", 0.42, 24.0, 327456),
    "AAG": ("K", 0.58, 33.2, 452789),
    "GAT": ("D", 0.45, 21.5, 293123),
    "GAC": ("D", 0.55, 26.3, 358456),
    "GAA": ("E", 0.41, 28.8, 392789),
    "GAG": ("E", 0.59, 41.5, 565789),
    "TGT": ("C", 0.44, 10.3, 140456),
    "TGC": ("C", 0.56, 13.1, 178456),
    "TGA": ("*", 0.47, 1.6, 21812),
    "TGG": ("W", 1.00, 13.0, 177456),
    "CGT": ("R", 0.09, 5.1, 69567),
    "CGC": ("R", 0.19, 11.5, 156789),
    "CGA": ("R", 0.10, 6.4, 87234),
    "CGG": ("R", 0.21, 12.5, 170456),
    "AGT": ("S", 0.16, 12.0, 163789),
    "AGC": ("S", 0.24, 19.5, 265789),
    "AGA": ("R", 0.21, 12.3, 167789),
    "AGG": ("R", 0.20, 12.0, 163456),
    "GGT": ("G", 0.16, 11.2, 152789),
    "GGC": ("G", 0.35, 23.8, 324456),
    "GGA": ("G", 0.25, 16.5, 225123),
    "GGG": ("G", 0.24, 16.0, 218456),
}

# Compute relative adaptiveness using shared utility
CHO_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(CHO_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
CHO_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(CHO_CODON_USAGE)
