"""
Corynebacterium glutamicum Codon Usage Data

Source: Kazusa Codon Usage Database
Corynebacterium glutamicum ATCC 13032, ~3,000 CDSs, coding GC: ~53.8%

C. glutamicum is a Gram-positive bacterium of immense industrial
importance for amino acid production (glutamate, lysine, etc.) and
increasingly for synthetic biology applications.  Its codon usage
reflects a moderately GC-rich genome.

Key features:
  - GC-ending codons preferred (GC ~53.8%)
  - CGC preferred for Arg (GC-rich Gram-positive pattern)
  - CCG common for Pro
  - AT-ending codons relatively rare

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 196627
  - Kalinowski et al. (2003) J Biotechnol 104:5-25
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "CORYNEBACTERIUM_CODON_USAGE",
    "CORYNEBACTERIUM_CODON_ADAPTIVENESS",
    "CORYNEBACTERIUM_PREFERRED_CODONS",
]

CORYNEBACTERIUM_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.38, 12.5, 30000),
    "TTC": ("F", 0.62, 20.5, 49200),
    "TTA": ("L", 0.03, 3.0, 7200),
    "TTG": ("L", 0.12, 12.8, 30720),
    "CTT": ("L", 0.07, 7.5, 18000),
    "CTC": ("L", 0.19, 20.2, 48480),
    "CTA": ("L", 0.04, 4.2, 10080),
    "CTG": ("L", 0.55, 58.5, 140400),
    "ATT": ("I", 0.30, 13.5, 32400),
    "ATC": ("I", 0.58, 26.2, 62880),
    "ATA": ("I", 0.12, 5.4, 12960),
    "ATG": ("M", 1.00, 22.5, 54000),
    "GTT": ("V", 0.18, 12.2, 29280),
    "GTC": ("V", 0.26, 17.8, 42720),
    "GTA": ("V", 0.09, 6.2, 14880),
    "GTG": ("V", 0.47, 32.2, 77280),
    "TCT": ("S", 0.12, 9.5, 22800),
    "TCC": ("S", 0.24, 18.8, 45120),
    "TCA": ("S", 0.08, 6.2, 14880),
    "TCG": ("S", 0.18, 14.0, 33600),
    "CCT": ("P", 0.10, 6.8, 16320),
    "CCC": ("P", 0.20, 13.5, 32400),
    "CCA": ("P", 0.25, 16.8, 40320),
    "CCG": ("P", 0.45, 30.2, 72480),
    "ACT": ("T", 0.14, 8.5, 20400),
    "ACC": ("T", 0.42, 25.5, 61200),
    "ACA": ("T", 0.10, 6.0, 14400),
    "ACG": ("T", 0.34, 20.5, 49200),
    "GCT": ("A", 0.18, 15.8, 37920),
    "GCC": ("A", 0.38, 33.2, 79680),
    "GCA": ("A", 0.20, 17.5, 42000),
    "GCG": ("A", 0.24, 21.0, 50400),
    "TAT": ("Y", 0.35, 9.8, 23520),
    "TAC": ("Y", 0.65, 18.2, 43680),
    "TAA": ("*", 0.52, 1.2, 2880),
    "TAG": ("*", 0.12, 0.3, 720),
    "CAT": ("H", 0.38, 8.5, 20400),
    "CAC": ("H", 0.62, 13.8, 33120),
    "CAA": ("Q", 0.22, 12.0, 28800),
    "CAG": ("Q", 0.78, 42.5, 102000),
    "AAT": ("N", 0.32, 13.5, 32400),
    "AAC": ("N", 0.68, 28.8, 69120),
    "AAA": ("K", 0.48, 25.5, 61200),
    "AAG": ("K", 0.52, 27.8, 66720),
    "GAT": ("D", 0.38, 20.0, 48000),
    "GAC": ("D", 0.62, 32.5, 78000),
    "GAA": ("E", 0.55, 35.0, 84000),
    "GAG": ("E", 0.45, 28.5, 68400),
    "TGT": ("C", 0.32, 5.0, 12000),
    "TGC": ("C", 0.68, 10.6, 25440),
    "TGA": ("*", 0.36, 0.9, 2160),
    "TGG": ("W", 1.00, 12.5, 30000),
    "CGT": ("R", 0.25, 16.0, 38400),
    "CGC": ("R", 0.42, 26.8, 64320),
    "CGA": ("R", 0.05, 3.2, 7680),
    "CGG": ("R", 0.12, 7.8, 18720),
    "AGA": ("R", 0.06, 3.8, 9120),
    "AGG": ("R", 0.10, 6.2, 14880),
    "AGT": ("S", 0.08, 6.5, 15600),
    "AGC": ("S", 0.30, 23.5, 56400),
    "GGT": ("G", 0.26, 22.0, 52800),
    "GGC": ("G", 0.42, 35.5, 85200),
    "GGA": ("G", 0.18, 15.2, 36480),
    "GGG": ("G", 0.14, 11.8, 28320),
}

CORYNEBACTERIUM_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(CORYNEBACTERIUM_CODON_USAGE)
CORYNEBACTERIUM_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(CORYNEBACTERIUM_CODON_USAGE)
