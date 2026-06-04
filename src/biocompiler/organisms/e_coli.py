"""
Escherichia coli Codon Usage Data

Source: Kazusa Codon Usage Database
K-12 MG1655, high-expression genes
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "E_COLI_CODON_USAGE",
    "E_COLI_CODON_ADAPTIVENESS",
    "E_COLI_PREFERRED_CODONS",
    "ECOLI_CODON_USAGE",
]

E_COLI_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.35, 22.0, 142302),
    "TTC": ("F", 0.65, 17.2, 111196),
    "TTA": ("L", 0.14, 13.7, 88522),
    "TTG": ("L", 0.13, 13.0, 84208),
    "CTT": ("L", 0.12, 11.3, 73210),
    "CTC": ("L", 0.10, 10.2, 66066),
    "CTA": ("L", 0.04, 3.7, 24008),
    "CTG": ("L", 0.47, 51.9, 335604),
    "ATT": ("I", 0.31, 29.7, 192124),
    "ATC": ("I", 0.63, 25.8, 166840),
    "ATA": ("I", 0.06, 4.2, 27278),
    "ATG": ("M", 1.00, 27.0, 174662),
    "GTT": ("V", 0.36, 18.6, 120532),
    "GTC": ("V", 0.26, 14.5, 93830),
    "GTA": ("V", 0.17, 17.8, 115160),
    "GTG": ("V", 0.37, 27.1, 175262),
    "TCT": ("S", 0.17, 10.3, 66428),
    "TCC": ("S", 0.16, 10.1, 65376),
    "TCA": ("S", 0.13, 7.6, 48868),
    "TCG": ("S", 0.14, 8.6, 55740),
    "CCT": ("P", 0.17, 7.3, 47268),
    "CCC": ("P", 0.12, 5.6, 36010),
    "CCA": ("P", 0.20, 8.5, 55060),
    "CCG": ("P", 0.51, 21.5, 139120),
    "ACT": ("T", 0.18, 12.2, 78750),
    "ACC": ("T", 0.43, 22.6, 146198),
    "ACA": ("T", 0.14, 8.1, 52520),
    "ACG": ("T", 0.26, 14.0, 90460),
    "GCT": ("A", 0.18, 15.6, 100968),
    "GCC": ("A", 0.26, 25.0, 161630),
    "GCA": ("A", 0.22, 20.4, 132050),
    "GCG": ("A", 0.34, 31.3, 202322),
    "TAT": ("Y", 0.42, 16.0, 103610),
    "TAC": ("Y", 0.58, 12.4, 79978),
    "TAA": ("*", 0.30, 2.0, 12896),
    "TAG": ("*", 0.17, 0.3, 1764),
    "CAT": ("H", 0.44, 12.2, 78888),
    "CAC": ("H", 0.56, 9.4, 60712),
    "CAA": ("Q", 0.30, 15.2, 98296),
    "CAG": ("Q", 0.70, 28.6, 185102),
    "AAT": ("N", 0.39, 18.6, 120394),
    "AAC": ("N", 0.61, 19.2, 124280),
    "AAA": ("K", 0.74, 33.6, 217604),
    "AAG": ("K", 0.26, 12.0, 77590),
    "GAT": ("D", 0.60, 32.4, 209784),
    "GAC": ("D", 0.40, 19.2, 124256),
    "GAA": ("E", 0.70, 38.6, 249610),
    "GAG": ("E", 0.30, 17.4, 112590),
    "TGT": ("C", 0.44, 5.2, 33700),
    "TGC": ("C", 0.56, 6.2, 40194),
    "TGA": ("*", 0.53, 1.1, 7032),
    "TGG": ("W", 1.00, 11.0, 71272),
    "CGT": ("R", 0.37, 20.4, 131864),
    "CGC": ("R", 0.36, 21.0, 135900),
    "CGA": ("R", 0.07, 3.6, 23456),
    "CGG": ("R", 0.11, 5.2, 33700),
    "AGT": ("S", 0.15, 8.2, 53160),
    "AGC": ("S", 0.25, 16.4, 106064),
    "AGA": ("R", 0.07, 2.9, 18700),
    "AGG": ("R", 0.03, 1.2, 7938),
    "GGT": ("G", 0.35, 24.4, 157960),
    "GGC": ("G", 0.37, 29.4, 190240),
    "GGA": ("G", 0.10, 7.8, 50396),
    "GGG": ("G", 0.15, 11.6, 75120),
}

# Compute adaptiveness using shared utility
E_COLI_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(E_COLI_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
E_COLI_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(E_COLI_CODON_USAGE)

# ────────────────────────────────────────────────────────────
# Legacy per-thousand codon usage (migrated from species.py)
# Different source dataset from E_COLI_CODON_USAGE above.
# Kept for backward compatibility with the species.py API.
# ────────────────────────────────────────────────────────────
ECOLI_CODON_USAGE: dict[str, float] = {
    "TTT": 17.6, "TTC": 20.3,
    "TTA": 7.6, "TTG": 11.0, "CTT": 10.5, "CTC": 10.5, "CTA": 3.9, "CTG": 51.0,
    "ATT": 29.8, "ATC": 25.1, "ATA": 4.2,
    "ATG": 27.0,
    "GTT": 18.3, "GTC": 15.0, "GTA": 10.8, "GTG": 27.8,
    "TCT": 8.5, "TCC": 8.5, "TCA": 7.3, "TCG": 4.3, "AGT": 9.6, "AGC": 15.4,
    "CCT": 7.0, "CCC": 5.5, "CCA": 8.4, "CCG": 23.2,
    "ACT": 12.9, "ACC": 25.7, "ACA": 7.1, "ACG": 6.3,
    "GCT": 18.5, "GCC": 27.1, "GCA": 20.2, "GCG": 7.4,
    "TAT": 16.3, "TAC": 14.9,
    "CAT": 13.5, "CAC": 9.8,
    "CAA": 14.6, "CAG": 29.0,
    "AAT": 17.1, "AAC": 21.3,
    "AAA": 33.5, "AAG": 24.1,
    "GAT": 31.0, "GAC": 21.4,
    "GAA": 39.2, "GAG": 19.6,
    "TGT": 5.1, "TGC": 5.5,
    "TGG": 12.9,
    "CGT": 20.0, "CGC": 21.5, "CGA": 3.5, "CGG": 5.4, "AGA": 2.1, "AGG": 1.2,
    "GGT": 24.5, "GGC": 28.6, "GGA": 8.0, "GGG": 6.8,
    "TAA": 2.0, "TAG": 0.3, "TGA": 1.1,
}
