"""
Bacillus subtilis Codon Usage Data

Source: Kazusa Codon Usage Database
Bacillus subtilis subsp. subtilis str. 168, ~4,200 CDSs, coding GC: ~43.5%

B. subtilis is a Gram-positive bacterium widely used as an industrial
host for enzyme production and as a model organism for Gram-positive
bacteria.  Its codon usage reflects a moderately GC-rich genome with
strong translational selection.

Key features:
  - GC-ending codons preferred for most amino acids
  - CUG strongly preferred for Leu
  - AGA/AGG very rare for Arg (unlike eukaryotes)
  - CGC preferred for Arg
  - AT-rich rare codons (TTA, ATA) are very uncommon

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 224308
  - Moszer et al. (2002) Nature 390:249-256
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "BACILLUS_CODON_USAGE",
    "BACILLUS_CODON_ADAPTIVENESS",
    "BACILLUS_PREFERRED_CODONS",
]

BACILLUS_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.38, 13.8, 48300),
    "TTC": ("F", 0.62, 22.5, 78750),
    "TTA": ("L", 0.04, 4.2, 14700),
    "TTG": ("L", 0.14, 15.2, 53200),
    "CTT": ("L", 0.08, 8.8, 30800),
    "CTC": ("L", 0.16, 17.8, 62300),
    "CTA": ("L", 0.03, 3.4, 11900),
    "CTG": ("L", 0.55, 60.0, 210000),
    "ATT": ("I", 0.35, 15.6, 54600),
    "ATC": ("I", 0.55, 24.5, 85750),
    "ATA": ("I", 0.10, 4.4, 15400),
    "ATG": ("M", 1.00, 25.2, 88200),
    "GTT": ("V", 0.22, 14.5, 50750),
    "GTC": ("V", 0.22, 14.8, 51800),
    "GTA": ("V", 0.12, 8.2, 28700),
    "GTG": ("V", 0.44, 29.5, 103250),
    "TCT": ("S", 0.14, 10.2, 35700),
    "TCC": ("S", 0.22, 16.0, 56000),
    "TCA": ("S", 0.08, 5.8, 20300),
    "TCG": ("S", 0.16, 11.5, 40250),
    "CCT": ("P", 0.14, 6.8, 23800),
    "CCC": ("P", 0.14, 6.9, 24150),
    "CCA": ("P", 0.25, 12.0, 42000),
    "CCG": ("P", 0.47, 22.8, 79800),
    "ACT": ("T", 0.18, 10.2, 35700),
    "ACC": ("T", 0.40, 22.5, 78750),
    "ACA": ("T", 0.10, 5.6, 19600),
    "ACG": ("T", 0.32, 18.0, 63000),
    "GCT": ("A", 0.18, 16.5, 57750),
    "GCC": ("A", 0.30, 27.2, 95200),
    "GCA": ("A", 0.20, 18.0, 63000),
    "GCG": ("A", 0.32, 29.0, 101500),
    "TAT": ("Y", 0.38, 12.0, 42000),
    "TAC": ("Y", 0.62, 19.5, 68250),
    "TAA": ("*", 0.58, 1.5, 5250),
    "TAG": ("*", 0.12, 0.3, 1050),
    "CAT": ("H", 0.42, 9.2, 32200),
    "CAC": ("H", 0.58, 12.8, 44800),
    "CAA": ("Q", 0.22, 10.8, 37800),
    "CAG": ("Q", 0.78, 38.2, 133700),
    "AAT": ("N", 0.36, 14.0, 49000),
    "AAC": ("N", 0.64, 25.0, 87500),
    "AAA": ("K", 0.58, 32.5, 113750),
    "AAG": ("K", 0.42, 23.5, 82250),
    "GAT": ("D", 0.42, 22.0, 77000),
    "GAC": ("D", 0.58, 30.5, 106750),
    "GAA": ("E", 0.62, 40.0, 140000),
    "GAG": ("E", 0.38, 24.5, 85750),
    "TGT": ("C", 0.38, 4.5, 15750),
    "TGC": ("C", 0.62, 7.3, 25550),
    "TGA": ("*", 0.30, 0.8, 2800),
    "TGG": ("W", 1.00, 13.0, 45500),
    "CGT": ("R", 0.32, 18.5, 64750),
    "CGC": ("R", 0.38, 22.0, 77000),
    "CGA": ("R", 0.06, 3.5, 12250),
    "CGG": ("R", 0.12, 7.0, 24500),
    "AGA": ("R", 0.04, 2.5, 8750),
    "AGG": ("R", 0.08, 4.8, 16800),
    "AGT": ("S", 0.12, 8.5, 29750),
    "AGC": ("S", 0.28, 20.0, 70000),
    "GGT": ("G", 0.28, 21.5, 75250),
    "GGC": ("G", 0.40, 30.5, 106750),
    "GGA": ("G", 0.18, 13.8, 48300),
    "GGG": ("G", 0.14, 10.8, 37800),
}

BACILLUS_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(BACILLUS_CODON_USAGE)
BACILLUS_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(BACILLUS_CODON_USAGE)
