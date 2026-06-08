"""
Pseudomonas putida Codon Usage Data

Source: Kazusa Codon Usage Database
Pseudomonas putida KT2440, ~5,300 CDSs, coding GC: ~61.5%

P. putida is a Gram-negative soil bacterium increasingly used as an
industrial host for biocatalysis and metabolic engineering.  Its
codon usage reflects a very GC-rich genome, one of the highest
among common biotech hosts.

Key features:
  - Very strong GC-ending codon preference (GC ~61.5%)
  - CGC strongly preferred for Arg
  - GCC strongly preferred for Ala
  - AT-ending codons (TTA, ATA, GTA) very rare
  - CCG very common for Pro

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 160488
  - Nelson et al. (2002) Environ Microbiol 4:799-808
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "PSEUDOMONAS_CODON_USAGE",
    "PSEUDOMONAS_CODON_ADAPTIVENESS",
    "PSEUDOMONAS_PREFERRED_CODONS",
]

PSEUDOMONAS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.30, 10.8, 43200),
    "TTC": ("F", 0.70, 25.2, 100800),
    "TTA": ("L", 0.03, 2.8, 11200),
    "TTG": ("L", 0.12, 11.0, 44000),
    "CTT": ("L", 0.06, 5.8, 23200),
    "CTC": ("L", 0.20, 18.2, 72800),
    "CTA": ("L", 0.04, 3.5, 14000),
    "CTG": ("L", 0.55, 50.5, 202000),
    "ATT": ("I", 0.22, 10.5, 42000),
    "ATC": ("I", 0.68, 32.5, 130000),
    "ATA": ("I", 0.10, 4.8, 19200),
    "ATG": ("M", 1.00, 23.0, 92000),
    "GTT": ("V", 0.14, 10.0, 40000),
    "GTC": ("V", 0.28, 20.2, 80800),
    "GTA": ("V", 0.07, 5.0, 20000),
    "GTG": ("V", 0.51, 36.8, 147200),
    "TCT": ("S", 0.10, 8.2, 32800),
    "TCC": ("S", 0.26, 21.0, 84000),
    "TCA": ("S", 0.06, 4.8, 19200),
    "TCG": ("S", 0.20, 16.2, 64800),
    "CCT": ("P", 0.08, 5.2, 20800),
    "CCC": ("P", 0.18, 11.8, 47200),
    "CCA": ("P", 0.20, 13.0, 52000),
    "CCG": ("P", 0.54, 35.0, 140000),
    "ACT": ("T", 0.12, 7.5, 30000),
    "ACC": ("T", 0.45, 28.0, 112000),
    "ACA": ("T", 0.08, 5.0, 20000),
    "ACG": ("T", 0.35, 22.0, 88000),
    "GCT": ("A", 0.14, 14.0, 56000),
    "GCC": ("A", 0.40, 40.0, 160000),
    "GCA": ("A", 0.16, 16.0, 64000),
    "GCG": ("A", 0.30, 30.0, 120000),
    "TAT": ("Y", 0.30, 9.5, 38000),
    "TAC": ("Y", 0.70, 22.2, 88800),
    "TAA": ("*", 0.50, 1.2, 4800),
    "TAG": ("*", 0.15, 0.4, 1600),
    "CAT": ("H", 0.32, 8.2, 32800),
    "CAC": ("H", 0.68, 17.5, 70000),
    "CAA": ("Q", 0.18, 10.5, 42000),
    "CAG": ("Q", 0.82, 47.8, 191200),
    "AAT": ("N", 0.28, 13.0, 52000),
    "AAC": ("N", 0.72, 33.5, 134000),
    "AAA": ("K", 0.35, 21.0, 84000),
    "AAG": ("K", 0.65, 39.0, 156000),
    "GAT": ("D", 0.35, 22.0, 88000),
    "GAC": ("D", 0.65, 40.8, 163200),
    "GAA": ("E", 0.52, 38.0, 152000),
    "GAG": ("E", 0.48, 35.0, 140000),
    "TGT": ("C", 0.30, 5.5, 22000),
    "TGC": ("C", 0.70, 12.8, 51200),
    "TGA": ("*", 0.35, 0.9, 3600),
    "TGG": ("W", 1.00, 14.5, 58000),
    "CGT": ("R", 0.28, 18.0, 72000),
    "CGC": ("R", 0.44, 28.5, 114000),
    "CGA": ("R", 0.06, 3.8, 15200),
    "CGG": ("R", 0.12, 7.8, 31200),
    "AGA": ("R", 0.04, 2.8, 11200),
    "AGG": ("R", 0.06, 4.0, 16000),
    "AGT": ("S", 0.08, 6.2, 24800),
    "AGC": ("S", 0.30, 24.0, 96000),
    "GGT": ("G", 0.26, 22.0, 88000),
    "GGC": ("G", 0.45, 38.0, 152000),
    "GGA": ("G", 0.14, 12.0, 48000),
    "GGG": ("G", 0.15, 12.8, 51200),
}

PSEUDOMONAS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(PSEUDOMONAS_CODON_USAGE)
PSEUDOMONAS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(PSEUDOMONAS_CODON_USAGE)
