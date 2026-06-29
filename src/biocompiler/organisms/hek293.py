"""
HEK293T Cell Line Codon Usage Data

Source: Derived from human high-expression codon usage with
HEK293T-specific adjustments based on ribosome profiling data.

HEK293T (Human Embryonic Kidney 293T) is one of the most widely used
mammalian cell lines for transient and stable recombinant protein
expression.  The 293T variant expresses SV40 large T antigen, enabling
episomal replication of SV40-origin plasmids.

Codon usage is similar to human high-expression genes but with
moderate adjustments reflecting the transformed cell line's
translation machinery preferences.

Key features:
  - Very similar to human high-expression codon usage
  - GC-ending codons preferred for most amino acids
  - AGA is the preferred Arg codon (consistent with human pattern)
  - CTG strongly preferred for Leu
  - CAG strongly preferred for Gln

References:
  - Thomas & Smart (2005) J Virol Methods 132:13-21
  - NCBI TaxID: 9606 (Homo sapiens — HEK293 is a human cell line)
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "HEK293_CODON_USAGE",
    "HEK293_CODON_ADAPTIVENESS",
    "HEK293_PREFERRED_CODONS",
]

HEK293_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.34, 12.8, 256),
    "TTC": ("F", 0.66, 24.9, 498),
    "TTA": ("L", 0.05, 4.2, 84),
    "TTG": ("L", 0.12, 11.4, 228),
    "CTT": ("L", 0.09, 8.4, 168),
    "CTC": ("L", 0.25, 26.2, 524),
    "CTA": ("L", 0.04, 3.6, 72),
    "CTG": ("L", 0.45, 47.0, 940),
    "ATT": ("I", 0.32, 14.0, 280),
    "ATC": ("I", 0.56, 24.8, 496),
    "ATA": ("I", 0.12, 5.4, 108),
    "ATG": ("M", 1.00, 21.5, 430),
    "GTT": ("V", 0.16, 9.4, 188),
    "GTC": ("V", 0.31, 19.2, 384),
    "GTA": ("V", 0.09, 5.4, 108),
    "GTG": ("V", 0.44, 27.2, 544),
    "TCT": ("S", 0.15, 12.0, 240),
    "TCC": ("S", 0.24, 19.6, 392),
    "TCA": ("S", 0.08, 6.4, 128),
    "TCG": ("S", 0.04, 3.2, 64),
    "CCT": ("P", 0.23, 13.8, 276),
    "CCC": ("P", 0.43, 26.0, 520),
    "CCA": ("P", 0.20, 12.0, 240),
    "CCG": ("P", 0.14, 8.4, 168),
    "ACT": ("T", 0.20, 10.2, 204),
    "ACC": ("T", 0.50, 25.6, 512),
    "ACA": ("T", 0.16, 8.2, 164),
    "ACG": ("T", 0.14, 7.2, 144),
    "GCT": ("A", 0.23, 15.4, 308),
    "GCC": ("A", 0.45, 30.4, 608),
    "GCA": ("A", 0.15, 10.2, 204),
    "GCG": ("A", 0.17, 11.4, 228),
    "TAT": ("Y", 0.34, 9.0, 180),
    "TAC": ("Y", 0.66, 17.4, 348),
    "TAA": ("*", 0.28, 0.9, 18),
    "TAG": ("*", 0.18, 0.6, 12),
    "CAT": ("H", 0.30, 7.6, 152),
    "CAC": ("H", 0.70, 17.8, 356),
    "CAA": ("Q", 0.19, 8.4, 168),
    "CAG": ("Q", 0.81, 35.8, 716),
    "AAT": ("N", 0.38, 13.2, 264),
    "AAC": ("N", 0.62, 21.6, 432),
    "AAA": ("K", 0.36, 20.0, 400),
    "AAG": ("K", 0.64, 35.6, 712),
    "GAT": ("D", 0.40, 18.4, 368),
    "GAC": ("D", 0.60, 27.6, 552),
    "GAA": ("E", 0.32, 21.4, 428),
    "GAG": ("E", 0.68, 45.4, 908),
    "TGT": ("C", 0.34, 7.6, 152),
    "TGC": ("C", 0.66, 14.8, 296),
    "TGA": ("*", 0.54, 1.8, 36),
    "TGG": ("W", 1.00, 12.8, 256),
    "CGT": ("R", 0.09, 5.0, 100),
    "CGC": ("R", 0.22, 12.8, 256),
    "CGA": ("R", 0.06, 3.6, 72),
    "CGG": ("R", 0.24, 14.2, 284),
    "AGA": ("R", 0.24, 14.2, 284),
    "AGG": ("R", 0.15, 8.8, 176),
    "AGT": ("S", 0.11, 8.6, 172),
    "AGC": ("S", 0.38, 30.4, 608),
    "GGT": ("G", 0.16, 10.0, 200),
    "GGC": ("G", 0.40, 25.2, 504),
    "GGA": ("G", 0.24, 15.0, 300),
    "GGG": ("G", 0.20, 12.6, 252),
}

HEK293_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(HEK293_CODON_USAGE)
HEK293_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(HEK293_CODON_USAGE)
