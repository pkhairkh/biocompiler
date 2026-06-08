"""
Pichia pastoris (Komagataella phaffii) Codon Usage Data

Source: Kazusa Codon Usage Database
Komagataella phaffii (Pichia pastoris) GS115
~5,000 CDSs, coding GC: ~41.5%

P. pastoris is one of the most important yeast expression systems for
recombinant protein production, especially for secreted proteins and
biopharmaceuticals.  Its codon usage reflects a moderately AT-rich
yeast pattern, similar to S. cerevisiae but with some key differences.

Key features:
  - A/T-ending codons generally preferred (yeast pattern)
  - AGA strongly preferred for Arg (like S. cerevisiae)
  - GCT preferred for Ala, GTT for Val (yeast consensus)
  - Stronger codon bias than S. cerevisiae for some amino acids

References:
  - Codon Usage Database: https://www.kazusa.or.jp/codon/
  - NCBI TaxID: 644223 (K. phaffii GS115)
  - Ahn et al. (2007) J Microbiol Biotechnol 17:1228-1235
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "PICHIA_CODON_USAGE",
    "PICHIA_CODON_ADAPTIVENESS",
    "PICHIA_PREFERRED_CODONS",
]

PICHIA_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.58, 18.9, 89000),
    "TTC": ("F", 0.42, 13.7, 64000),
    "TTA": ("L", 0.25, 10.8, 51000),
    "TTG": ("L", 0.32, 13.9, 65000),
    "CTT": ("L", 0.12, 5.2, 24000),
    "CTC": ("L", 0.06, 2.6, 12000),
    "CTA": ("L", 0.04, 1.7, 8000),
    "CTG": ("L", 0.21, 9.1, 43000),
    "ATT": ("I", 0.55, 20.5, 96000),
    "ATC": ("I", 0.38, 14.2, 67000),
    "ATA": ("I", 0.07, 2.6, 12000),
    "ATG": ("M", 1.00, 20.2, 95000),
    "GTT": ("V", 0.42, 18.6, 87000),
    "GTC": ("V", 0.20, 8.9, 42000),
    "GTA": ("V", 0.17, 7.5, 35000),
    "GTG": ("V", 0.21, 9.3, 44000),
    "TCT": ("S", 0.32, 17.5, 82000),
    "TCC": ("S", 0.19, 10.4, 49000),
    "TCA": ("S", 0.17, 9.3, 44000),
    "TCG": ("S", 0.04, 2.2, 10000),
    "CCT": ("P", 0.25, 10.8, 51000),
    "CCC": ("P", 0.12, 5.2, 24000),
    "CCA": ("P", 0.48, 20.8, 98000),
    "CCG": ("P", 0.15, 6.5, 30000),
    "ACT": ("T", 0.42, 18.2, 85000),
    "ACC": ("T", 0.24, 10.4, 49000),
    "ACA": ("T", 0.24, 10.4, 49000),
    "ACG": ("T", 0.10, 4.3, 20000),
    "GCT": ("A", 0.40, 20.0, 94000),
    "GCC": ("A", 0.20, 10.0, 47000),
    "GCA": ("A", 0.26, 13.0, 61000),
    "GCG": ("A", 0.14, 7.0, 33000),
    "TAT": ("Y", 0.55, 13.0, 61000),
    "TAC": ("Y", 0.45, 10.6, 50000),
    "TAA": ("*", 0.48, 1.2, 6000),
    "TAG": ("*", 0.18, 0.4, 2000),
    "CAT": ("H", 0.60, 10.4, 49000),
    "CAC": ("H", 0.40, 6.9, 32000),
    "CAA": ("Q", 0.72, 24.8, 116000),
    "CAG": ("Q", 0.28, 9.6, 45000),
    "AAT": ("N", 0.58, 17.8, 83000),
    "AAC": ("N", 0.42, 12.9, 60000),
    "AAA": ("K", 0.62, 28.6, 134000),
    "AAG": ("K", 0.38, 17.5, 82000),
    "GAT": ("D", 0.60, 24.5, 115000),
    "GAC": ("D", 0.40, 16.3, 76000),
    "GAA": ("E", 0.68, 34.2, 160000),
    "GAG": ("E", 0.32, 16.1, 75000),
    "TGT": ("C", 0.58, 6.9, 32000),
    "TGC": ("C", 0.42, 5.0, 23000),
    "TGA": ("*", 0.34, 0.8, 4000),
    "TGG": ("W", 1.00, 9.6, 45000),
    "CGT": ("R", 0.12, 4.3, 20000),
    "CGC": ("R", 0.04, 1.4, 7000),
    "CGA": ("R", 0.03, 1.1, 5000),
    "CGG": ("R", 0.04, 1.4, 7000),
    "AGA": ("R", 0.58, 20.8, 97000),
    "AGG": ("R", 0.19, 6.8, 32000),
    "AGT": ("S", 0.17, 9.3, 43000),
    "AGC": ("S", 0.11, 6.0, 28000),
    "GGT": ("G", 0.42, 15.6, 73000),
    "GGC": ("G", 0.16, 5.9, 28000),
    "GGA": ("G", 0.28, 10.4, 49000),
    "GGG": ("G", 0.14, 5.2, 24000),
}

PICHIA_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(PICHIA_CODON_USAGE)
PICHIA_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(PICHIA_CODON_USAGE)
