"""
PER.C6 Human Cell Line Codon Usage Data

Source: Human codon usage with PER.C6-specific adjustments.
NCBI TaxID: 9606 (Homo sapiens — PER.C6 is a human cell line)

PER.C6 is a human retinoblast cell line (derived from human fetal retina)
used for recombinant protein production, particularly adenoviral vectors
and monoclonal antibodies.  Its codon usage follows the human pattern
with GC-ending codons preferred.

Key features:
  - Human codon usage pattern (GC-ending preferred)
  - Suitable for adenoviral vector production
  - CTG strongly preferred for Leu, CAG for Gln
  - AGA/CGG both significant for Arg

References:
  - Fallaux et al. (1998) Hum Gene Ther 9:1909-1917
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "PER_C6_CODON_USAGE",
    "PER_C6_CODON_ADAPTIVENESS",
    "PER_C6_PREFERRED_CODONS",
]

PER_C6_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.44, 16.8, 67200),
    "TTC": ("F", 0.56, 21.4, 85600),
    "TTA": ("L", 0.07, 7.2, 28800),
    "TTG": ("L", 0.12, 12.0, 48000),
    "CTT": ("L", 0.12, 12.4, 49600),
    "CTC": ("L", 0.21, 20.2, 80800),
    "CTA": ("L", 0.07, 7.0, 28000),
    "CTG": ("L", 0.41, 40.0, 160000),
    "ATT": ("I", 0.35, 15.6, 62400),
    "ATC": ("I", 0.49, 21.8, 87200),
    "ATA": ("I", 0.16, 7.2, 28800),
    "ATG": ("M", 1.00, 22.0, 88000),
    "GTT": ("V", 0.18, 11.2, 44800),
    "GTC": ("V", 0.24, 14.8, 59200),
    "GTA": ("V", 0.11, 6.8, 27200),
    "GTG": ("V", 0.47, 29.2, 116800),
    "TCT": ("S", 0.18, 14.6, 58400),
    "TCC": ("S", 0.23, 18.4, 73600),
    "TCA": ("S", 0.14, 11.4, 45600),
    "TCG": ("S", 0.05, 4.0, 16000),
    "CCT": ("P", 0.28, 17.0, 68000),
    "CCC": ("P", 0.34, 20.8, 83200),
    "CCA": ("P", 0.28, 17.0, 68000),
    "CCG": ("P", 0.10, 6.2, 24800),
    "ACT": ("T", 0.24, 13.0, 52000),
    "ACC": ("T", 0.38, 20.4, 81600),
    "ACA": ("T", 0.27, 14.4, 57600),
    "ACG": ("T", 0.11, 5.8, 23200),
    "GCT": ("A", 0.26, 18.0, 72000),
    "GCC": ("A", 0.42, 29.2, 116800),
    "GCA": ("A", 0.22, 15.2, 60800),
    "GCG": ("A", 0.10, 6.8, 27200),
    "TAT": ("Y", 0.43, 11.8, 47200),
    "TAC": ("Y", 0.57, 15.6, 62400),
    "TAA": ("*", 0.28, 1.0, 4000),
    "TAG": ("*", 0.22, 0.8, 3200),
    "CAT": ("H", 0.41, 10.6, 42400),
    "CAC": ("H", 0.59, 15.2, 60800),
    "CAA": ("Q", 0.25, 12.0, 48000),
    "CAG": ("Q", 0.75, 36.0, 144000),
    "AAT": ("N", 0.46, 17.2, 68800),
    "AAC": ("N", 0.54, 20.2, 80800),
    "AAA": ("K", 0.42, 24.2, 96800),
    "AAG": ("K", 0.58, 33.4, 133600),
    "GAT": ("D", 0.44, 21.4, 85600),
    "GAC": ("D", 0.56, 27.2, 108800),
    "GAA": ("E", 0.41, 29.4, 117600),
    "GAG": ("E", 0.59, 42.4, 169600),
    "TGT": ("C", 0.44, 10.2, 40800),
    "TGC": ("C", 0.56, 12.8, 51200),
    "TGA": ("*", 0.50, 1.6, 6400),
    "TGG": ("W", 1.00, 13.0, 52000),
    "CGT": ("R", 0.08, 4.8, 19200),
    "CGC": ("R", 0.18, 10.8, 43200),
    "CGA": ("R", 0.10, 6.2, 24800),
    "CGG": ("R", 0.21, 12.6, 50400),
    "AGA": ("R", 0.24, 14.4, 57600),
    "AGG": ("R", 0.19, 11.4, 45600),
    "AGT": ("S", 0.15, 11.6, 46400),
    "AGC": ("S", 0.25, 19.6, 78400),
    "GGT": ("G", 0.16, 11.2, 44800),
    "GGC": ("G", 0.35, 24.4, 97600),
    "GGA": ("G", 0.26, 18.0, 72000),
    "GGG": ("G", 0.23, 16.0, 64000),
}

PER_C6_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(PER_C6_CODON_USAGE)
PER_C6_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(PER_C6_CODON_USAGE)
