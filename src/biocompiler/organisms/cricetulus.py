"""
Cricetulus griseus (Chinese Hamster) Codon Usage Data

Source: Kazusa Codon Usage Database
Cricetulus griseus, ~8,000 CDSs, coding GC: ~52.7%

C. griseus is the parent species of CHO (Chinese Hamster Ovary) cells,
the dominant mammalian host for biopharmaceutical production.  This
module provides the wild-type Chinese hamster codon usage, which is
very similar to CHO-K1 but may differ slightly in some codon
frequencies due to the cell-line adaptation process.

Key features:
  - Nearly identical to CHO-K1 codon usage
  - GC-ending codons preferred (mammalian pattern)
  - AGA preferred for Arg
  - CTG strongly preferred for Leu

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 10029
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CRICETULUS_CODON_USAGE",
    "CRICETULUS_CODON_ADAPTIVENESS",
    "CRICETULUS_PREFERRED_CODONS",
]

CRICETULUS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.45, 17.0, 74800),
    "TTC": ("F", 0.55, 20.8, 91520),
    "TTA": ("L", 0.07, 7.2, 31680),
    "TTG": ("L", 0.13, 12.6, 55440),
    "CTT": ("L", 0.13, 12.2, 53680),
    "CTC": ("L", 0.21, 19.6, 86240),
    "CTA": ("L", 0.07, 7.0, 30800),
    "CTG": ("L", 0.39, 38.2, 168080),
    "ATT": ("I", 0.35, 15.4, 67760),
    "ATC": ("I", 0.49, 21.6, 95040),
    "ATA": ("I", 0.16, 7.0, 30800),
    "ATG": ("M", 1.00, 21.6, 95040),
    "GTT": ("V", 0.18, 11.2, 49280),
    "GTC": ("V", 0.24, 14.6, 64240),
    "GTA": ("V", 0.11, 6.8, 29920),
    "GTG": ("V", 0.47, 28.6, 125840),
    "TCT": ("S", 0.18, 14.2, 62480),
    "TCC": ("S", 0.23, 18.2, 80080),
    "TCA": ("S", 0.15, 11.8, 51920),
    "TCG": ("S", 0.05, 4.0, 17600),
    "CCT": ("P", 0.28, 17.0, 74800),
    "CCC": ("P", 0.33, 20.2, 88880),
    "CCA": ("P", 0.27, 16.4, 72160),
    "CCG": ("P", 0.12, 7.4, 32560),
    "ACT": ("T", 0.24, 12.6, 55440),
    "ACC": ("T", 0.38, 19.6, 86240),
    "ACA": ("T", 0.27, 14.0, 61600),
    "ACG": ("T", 0.11, 5.6, 24640),
    "GCT": ("A", 0.25, 17.6, 77440),
    "GCC": ("A", 0.42, 29.2, 128480),
    "GCA": ("A", 0.22, 15.2, 66880),
    "GCG": ("A", 0.11, 7.6, 33440),
    "TAT": ("Y", 0.44, 11.8, 51920),
    "TAC": ("Y", 0.56, 15.2, 66880),
    "TAA": ("*", 0.28, 1.0, 4400),
    "TAG": ("*", 0.22, 0.8, 3520),
    "CAT": ("H", 0.41, 10.4, 45760),
    "CAC": ("H", 0.59, 15.0, 66000),
    "CAA": ("Q", 0.25, 12.0, 52800),
    "CAG": ("Q", 0.75, 35.6, 156640),
    "AAT": ("N", 0.46, 17.0, 74800),
    "AAC": ("N", 0.54, 20.0, 88000),
    "AAA": ("K", 0.42, 24.0, 105600),
    "AAG": ("K", 0.58, 33.0, 145200),
    "GAT": ("D", 0.45, 21.4, 94160),
    "GAC": ("D", 0.55, 26.2, 115280),
    "GAA": ("E", 0.41, 28.6, 125840),
    "GAG": ("E", 0.59, 41.2, 181280),
    "TGT": ("C", 0.44, 10.2, 44880),
    "TGC": ("C", 0.56, 13.0, 57200),
    "TGA": ("*", 0.50, 1.6, 7040),
    "TGG": ("W", 1.00, 12.8, 56320),
    "CGT": ("R", 0.09, 5.2, 22880),
    "CGC": ("R", 0.19, 11.2, 49280),
    "CGA": ("R", 0.10, 6.2, 27280),
    "CGG": ("R", 0.20, 12.0, 52800),
    "AGA": ("R", 0.23, 14.2, 62480),
    "AGG": ("R", 0.19, 11.8, 51920),
    "AGT": ("S", 0.15, 11.8, 51920),
    "AGC": ("S", 0.24, 19.2, 84480),
    "GGT": ("G", 0.16, 11.0, 48400),
    "GGC": ("G", 0.35, 23.6, 103840),
    "GGA": ("G", 0.25, 16.4, 72160),
    "GGG": ("G", 0.24, 15.8, 69520),
}

CRICETULUS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(CRICETULUS_CODON_USAGE)
CRICETULUS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(CRICETULUS_CODON_USAGE)
