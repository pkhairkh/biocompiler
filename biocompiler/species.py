"""
BioCompiler Species Data v8.1.0
================================
Codon usage tables and CAI reference values for common expression hosts.

Supports: E. coli, Human (H. sapiens), S. cerevisiae (Yeast), CHO-K1 (C. griseus)
"""

from typing import Dict

# E. coli codon usage (per thousand, K-12 MG1655)
ECOLI_CODON_USAGE: Dict[str, float] = {
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


def compute_cai_weights(usage: Dict[str, float]) -> Dict[str, float]:
    """Compute CAI weights from codon usage. Most frequent codon per AA = 1.0."""
    from .type_system import AA_TO_CODONS
    weights: Dict[str, float] = {}
    for aa, codons in AA_TO_CODONS.items():
        if aa == "*":
            continue
        freqs = [usage.get(c, 0.1) for c in codons]
        max_freq = max(freqs) if freqs else 1.0
        for codon, freq in zip(codons, freqs):
            weights[codon] = freq / max_freq if max_freq > 0 else 0.0
    return weights


ECOLI_CAI: Dict[str, float] = compute_cai_weights(ECOLI_CODON_USAGE)

HUMAN_CODON_USAGE: Dict[str, float] = {
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

HUMAN_CAI: Dict[str, float] = compute_cai_weights(HUMAN_CODON_USAGE)

# S. cerevisiae (Yeast) codon usage (per thousand)
# Source: Nakamura et al., Codon usage tabulated from GenBank
YEAST_CODON_USAGE: Dict[str, float] = {
    "TTT": 18.2, "TTC": 21.3,
    "TTA": 13.6, "TTG": 27.5, "CTT": 12.1, "CTC": 5.4, "CTA": 7.0, "CTG": 10.7,
    "ATT": 30.1, "ATC": 17.0, "ATA": 17.8,
    "ATG": 20.9,
    "GTT": 22.0, "GTC": 11.8, "GTA": 10.4, "GTG": 10.8,
    "TCT": 23.5, "TCC": 14.2, "TCA": 18.4, "TCG": 8.5, "AGT": 14.7, "AGC": 9.7,
    "CCT": 15.3, "CCC": 6.8, "CCA": 31.8, "CCG": 5.3,
    "ACT": 20.1, "ACC": 12.8, "ACA": 17.8, "ACG": 8.0,
    "GCT": 35.8, "GCC": 22.4, "GCA": 30.9, "GCG": 6.1,
    "TAT": 18.8, "TAC": 14.9,
    "CAT": 13.6, "CAC": 7.9,
    "CAA": 27.1, "CAG": 12.1,
    "AAT": 35.7, "AAC": 24.8,
    "AAA": 41.9, "AAG": 30.7,
    "GAT": 37.0, "GAC": 20.2,
    "GAA": 45.6, "GAG": 19.2,
    "TGT": 6.6, "TGC": 4.7,
    "TGG": 10.2,
    "CGT": 6.4, "CGC": 2.6, "CGA": 3.0, "CGG": 1.7, "AGA": 21.3, "AGG": 9.2,
    "GGT": 23.8, "GGC": 9.6, "GGA": 29.2, "GGG": 6.0,
    "TAA": 1.0, "TAG": 0.5, "TGA": 0.7,
}

YEAST_CAI: Dict[str, float] = compute_cai_weights(YEAST_CODON_USAGE)

# CHO-K1 (Chinese Hamster Ovary) codon usage (per thousand)
# Source: Codon usage tabulated from Cricetulus griseus GenBank entries
CHO_CODON_USAGE: Dict[str, float] = {
    "TTT": 16.4, "TTC": 21.2,
    "TTA": 7.0, "TTG": 12.5, "CTT": 12.8, "CTC": 19.6, "CTA": 7.2, "CTG": 40.8,
    "ATT": 15.5, "ATC": 22.3, "ATA": 6.8,
    "ATG": 21.8,
    "GTT": 10.6, "GTC": 15.2, "GTA": 6.9, "GTG": 28.1,
    "TCT": 14.6, "TCC": 17.8, "TCA": 11.4, "TCG": 4.3, "AGT": 11.8, "AGC": 19.7,
    "CCT": 17.0, "CCC": 20.1, "CCA": 16.4, "CCG": 7.2,
    "ACT": 12.6, "ACC": 19.2, "ACA": 14.6, "ACG": 6.0,
    "GCT": 18.0, "GCC": 28.2, "GCA": 15.5, "GCG": 7.3,
    "TAT": 15.0, "TAC": 16.0,
    "CAT": 10.2, "CAC": 15.5,
    "CAA": 11.5, "CAG": 35.2,
    "AAT": 16.4, "AAC": 20.0,
    "AAA": 23.5, "AAG": 33.0,
    "GAT": 21.2, "GAC": 25.8,
    "GAA": 28.5, "GAG": 40.5,
    "TGT": 10.0, "TGC": 12.8,
    "TGG": 13.0,
    "CGT": 4.3, "CGC": 10.8, "CGA": 5.9, "CGG": 11.5, "AGA": 11.5, "AGG": 12.2,
    "GGT": 10.5, "GGC": 22.8, "GGA": 16.2, "GGG": 16.8,
    "TAA": 1.4, "TAG": 0.7, "TGA": 1.2,
}

CHO_CAI: Dict[str, float] = compute_cai_weights(CHO_CODON_USAGE)

SPECIES: Dict[str, Dict[str, float]] = {
    "ecoli": ECOLI_CAI,
    "human": HUMAN_CAI,
    "yeast": YEAST_CAI,
    "cho": CHO_CAI,
}
