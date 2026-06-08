"""
NS0 Mouse Myeloma Cell Line Codon Usage Data

Source: Kazusa Codon Usage Database / Mus musculus with NS0-specific adjustments.
NCBI TaxID: 10090 (Mus musculus)

NS0 is a mouse myeloma cell line widely used for the production of
therapeutic monoclonal antibodies and other recombinant proteins.
Its codon usage is similar to mouse high-expression genes but with
some myeloma-specific biases.

Key features:
  - Similar to Mus musculus codon usage (mammalian pattern)
  - GC-ending codons preferred for most amino acids
  - AGA is the preferred Arg codon
  - CTG strongly preferred for Leu

References:
  - Barnes et al. (2000) Cytotechnology 33:13-23
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "NS0_CODON_USAGE",
    "NS0_CODON_ADAPTIVENESS",
    "NS0_PREFERRED_CODONS",
]

NS0_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.44, 16.8, 82320),
    "TTC": ("F", 0.56, 21.4, 104860),
    "TTA": ("L", 0.07, 7.2, 35280),
    "TTG": ("L", 0.13, 12.8, 62720),
    "CTT": ("L", 0.12, 12.0, 58800),
    "CTC": ("L", 0.21, 20.5, 100450),
    "CTA": ("L", 0.07, 6.8, 33320),
    "CTG": ("L", 0.40, 39.2, 192080),
    "ATT": ("I", 0.36, 16.0, 78400),
    "ATC": ("I", 0.48, 21.2, 103880),
    "ATA": ("I", 0.16, 7.0, 34300),
    "ATG": ("M", 1.00, 21.5, 105350),
    "GTT": ("V", 0.19, 11.5, 56350),
    "GTC": ("V", 0.24, 14.8, 72520),
    "GTA": ("V", 0.11, 6.8, 33320),
    "GTG": ("V", 0.46, 28.5, 139650),
    "TCT": ("S", 0.19, 14.8, 72520),
    "TCC": ("S", 0.22, 17.2, 84280),
    "TCA": ("S", 0.14, 11.0, 53900),
    "TCG": ("S", 0.05, 3.9, 19110),
    "CCT": ("P", 0.27, 16.2, 79380),
    "CCC": ("P", 0.34, 20.4, 99960),
    "CCA": ("P", 0.27, 16.2, 79380),
    "CCG": ("P", 0.12, 7.2, 35280),
    "ACT": ("T", 0.24, 12.5, 61250),
    "ACC": ("T", 0.38, 19.8, 97020),
    "ACA": ("T", 0.27, 14.0, 68600),
    "ACG": ("T", 0.11, 5.7, 27930),
    "GCT": ("A", 0.25, 17.2, 84280),
    "GCC": ("A", 0.42, 29.0, 142100),
    "GCA": ("A", 0.22, 15.2, 74480),
    "GCG": ("A", 0.11, 7.6, 37240),
    "TAT": ("Y", 0.43, 11.8, 57820),
    "TAC": ("Y", 0.57, 15.6, 76440),
    "TAA": ("*", 0.30, 1.0, 4900),
    "TAG": ("*", 0.22, 0.7, 3430),
    "CAT": ("H", 0.42, 10.4, 50960),
    "CAC": ("H", 0.58, 14.4, 70560),
    "CAA": ("Q", 0.26, 12.4, 60760),
    "CAG": ("Q", 0.74, 35.4, 173460),
    "AAT": ("N", 0.45, 17.0, 83300),
    "AAC": ("N", 0.55, 20.8, 101920),
    "AAA": ("K", 0.42, 24.2, 118580),
    "AAG": ("K", 0.58, 33.4, 163660),
    "GAT": ("D", 0.44, 21.2, 103880),
    "GAC": ("D", 0.56, 27.0, 132300),
    "GAA": ("E", 0.42, 29.2, 143080),
    "GAG": ("E", 0.58, 40.4, 197960),
    "TGT": ("C", 0.44, 9.8, 48020),
    "TGC": ("C", 0.56, 12.5, 61250),
    "TGA": ("*", 0.48, 1.5, 7350),
    "TGG": ("W", 1.00, 12.8, 62720),
    "CGT": ("R", 0.09, 5.3, 25970),
    "CGC": ("R", 0.18, 10.8, 52920),
    "CGA": ("R", 0.10, 6.0, 29400),
    "CGG": ("R", 0.19, 11.5, 56350),
    "AGA": ("R", 0.26, 15.8, 77420),
    "AGG": ("R", 0.18, 10.9, 53410),
    "AGT": ("S", 0.15, 11.5, 56350),
    "AGC": ("S", 0.25, 19.4, 95060),
    "GGT": ("G", 0.17, 11.5, 56350),
    "GGC": ("G", 0.34, 23.2, 113680),
    "GGA": ("G", 0.26, 17.6, 86240),
    "GGG": ("G", 0.23, 15.6, 76440),
}

NS0_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(NS0_CODON_USAGE)
NS0_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(NS0_CODON_USAGE)
