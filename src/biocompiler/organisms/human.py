"""
Human (Homo sapiens) Codon Usage Data

Source: Kazusa Codon Usage Database
93,487 CDSs, 40,662,582 codons
Coding GC: 52.27%
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "HUMAN_CODON_USAGE",
    "HUMAN_CODON_ADAPTIVENESS",
    "HUMAN_PREFERRED_CODONS",
    "HUMAN_CODON_USAGE_SIMPLE",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
HUMAN_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.46, 17.6, 714298),
    "TTC": ("F", 0.54, 20.3, 824692),
    "TTA": ("L", 0.08, 7.7, 311881),
    "TTG": ("L", 0.13, 12.9, 525688),
    "CTT": ("L", 0.13, 13.2, 536515),
    "CTC": ("L", 0.20, 19.6, 796638),
    "CTA": ("L", 0.07, 7.2, 290751),
    "CTG": ("L", 0.40, 39.6, 1611801),
    "ATT": ("I", 0.36, 16.0, 650473),
    "ATC": ("I", 0.47, 20.8, 846466),
    "ATA": ("I", 0.17, 7.5, 304565),
    "ATG": ("M", 1.00, 22.0, 896005),
    "GTT": ("V", 0.18, 11.0, 448607),
    "GTC": ("V", 0.24, 14.5, 588138),
    "GTA": ("V", 0.12, 7.1, 287712),
    "GTG": ("V", 0.46, 28.1, 1143534),
    "TCT": ("S", 0.19, 15.2, 618711),
    "TCC": ("S", 0.22, 17.7, 718892),
    "TCA": ("S", 0.15, 12.2, 496448),
    "TCG": ("S", 0.05, 4.4, 179419),
    "CCT": ("P", 0.29, 17.5, 713233),
    "CCC": ("P", 0.32, 19.8, 804620),
    "CCA": ("P", 0.28, 16.9, 688038),
    "CCG": ("P", 0.11, 6.9, 281570),
    "ACT": ("T", 0.25, 13.1, 533609),
    "ACC": ("T", 0.36, 18.9, 768147),
    "ACA": ("T", 0.28, 15.1, 614523),
    "ACG": ("T", 0.11, 6.1, 246105),
    "GCT": ("A", 0.27, 18.4, 750096),
    "GCC": ("A", 0.40, 27.7, 1127679),
    "GCA": ("A", 0.23, 15.8, 643471),
    "GCG": ("A", 0.11, 7.4, 299495),
    "TAT": ("Y", 0.44, 12.2, 495699),
    "TAC": ("Y", 0.56, 15.3, 622407),
    "TAA": ("*", 0.30, 1.0, 40285),
    "TAG": ("*", 0.24, 0.8, 32109),
    "CAT": ("H", 0.42, 10.9, 441711),
    "CAC": ("H", 0.58, 15.1, 613713),
    "CAA": ("Q", 0.27, 12.3, 501911),
    "CAG": ("Q", 0.73, 34.2, 1391973),
    "AAT": ("N", 0.47, 17.0, 689701),
    "AAC": ("N", 0.53, 19.1, 776603),
    "AAA": ("K", 0.43, 24.4, 993621),
    "AAG": ("K", 0.57, 31.9, 1295568),
    "GAT": ("D", 0.46, 21.8, 885429),
    "GAC": ("D", 0.54, 25.1, 1020595),
    "GAA": ("E", 0.42, 29.0, 1177632),
    "GAG": ("E", 0.58, 39.6, 1609975),
    "TGT": ("C", 0.46, 10.6, 430311),
    "TGC": ("C", 0.54, 12.6, 513028),
    "TGA": ("*", 0.47, 1.6, 63237),
    "TGG": ("W", 1.00, 13.2, 535595),
    "CGT": ("R", 0.08, 4.5, 184609),
    "CGC": ("R", 0.18, 10.4, 423516),
    "CGA": ("R", 0.11, 6.2, 250760),
    "CGG": ("R", 0.20, 11.4, 464485),
    "AGT": ("S", 0.15, 12.1, 493429),
    "AGC": ("S", 0.24, 19.5, 791383),
    "AGA": ("R", 0.21, 12.2, 494682),
    "AGG": ("R", 0.21, 12.0, 486463),
    "GGT": ("G", 0.16, 10.8, 437126),
    "GGC": ("G", 0.34, 22.2, 903565),
    "GGA": ("G", 0.25, 16.5, 669873),
    "GGG": ("G", 0.25, 16.5, 669768),
}

# Compute relative adaptiveness: w_i = freq_i / max_freq_for_same_aa
HUMAN_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(HUMAN_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
HUMAN_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(HUMAN_CODON_USAGE)

# ────────────────────────────────────────────────────────────
# Legacy per-thousand codon usage (migrated from species.py)
# Different source dataset from HUMAN_CODON_USAGE above.
# Named HUMAN_CODON_USAGE_SIMPLE to avoid clash with the
# richer tuple-format HUMAN_CODON_USAGE.
# ────────────────────────────────────────────────────────────
HUMAN_CODON_USAGE_SIMPLE: dict[str, float] = {
    "TTT": 17.2, "TTC": 20.8,
    "TTA": 7.4, "TTG": 12.9, "CTT": 13.0, "CTC": 19.4, "CTA": 7.5, "CTG": 39.4,
    "ATT": 16.0, "ATC": 21.0, "ATA": 7.1,
    "ATG": 22.3,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.5,
    "TCT": 14.9, "TCC": 17.4, "TCA": 11.7, "TCG": 4.5, "AGT": 12.0, "AGC": 19.3,
    "CCT": 17.3, "CCC": 19.7, "CCA": 16.7, "CCG": 7.0,
    "ACT": 12.9, "ACC": 18.6, "ACA": 14.8, "ACG": 6.2,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 15.4, "TAC": 15.6,
    "CAT": 10.5, "CAC": 15.0,
    "CAA": 11.8, "CAG": 34.3,
    "AAT": 16.8, "AAC": 19.5,
    "AAA": 24.1, "AAG": 32.1,
    "GAT": 21.5, "GAC": 25.4,
    "GAA": 28.8, "GAG": 39.8,
    "TGT": 10.2, "TGC": 12.4,
    "TGG": 13.4,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.1, "CGG": 11.3, "AGA": 11.7, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.4, "GGG": 16.5,
    "TAA": 1.5, "TAG": 0.7, "TGA": 1.3,
}
