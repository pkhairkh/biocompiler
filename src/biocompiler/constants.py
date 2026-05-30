"""
BioCompiler Constants

Single canonical source for all biological constants.
No more duplication across poc/data/mvp modules.
"""

# ==============================================================================
# 1. Genetic Code
# ==============================================================================

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

# Reverse lookup: amino acid -> list of codons
AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    if _aa != "*":
        AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# ==============================================================================
# 2. Scanner Consensus Sequences
# ==============================================================================

DONOR_CONSENSUS = "GT"
ACCEPTOR_CONSENSUS = "AG"
KOZAK_CONSENSUS = "GCCACC"
INSTABILITY_MOTIF = "ATTTA"
MIN_INTRON_LENGTH = 30  # Minimum intron length in nt
POLYPYRIMIDINE_WINDOW = 40  # Upstream window for acceptor scoring
POLYPYRIMIDINE_THRESHOLD = 0.5

# ==============================================================================
# 3. Restriction Enzyme Sites
# ==============================================================================

RESTRICTION_ENZYMES: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "XhoI": "CTCGAG",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "NdeI": "CATATG",
}

# ==============================================================================
# 4. Nucleotide encoding for CSP solver
# ==============================================================================

BASE_MAP = {0: "A", 1: "C", 2: "G", 3: "T"}
BASE_REV = {"A": 0, "C": 1, "G": 2, "T": 3}

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G",
              "a": "t", "t": "a", "g": "c", "c": "g",
              "N": "N", "n": "n"}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(COMPLEMENT[base] for base in reversed(seq))
