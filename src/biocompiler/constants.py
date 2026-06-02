"""
BioCompiler Constants

Single canonical source for all biological constants.
No more duplication across poc/data/mvp modules.

Extended with:
- More restriction enzymes from REBASE
- IUPAC ambiguity code support
- Safer reverse_complement with validation
- Organism-specific min intron lengths
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
KOZAK_CONSENSUS = "GCCACC"  # Exact match fallback; scanner uses PWM scoring
INSTABILITY_MOTIF = "ATTTA"

# Minimum intron lengths by organism (nt)
# Source: Nucleic Acids Research, various publications
MIN_INTRON_LENGTHS: dict[str, int] = {
    "Homo_sapiens": 30,
    "Mus_musculus": 30,
    "Drosophila_melanogaster": 30,
    "Saccharomyces_cerevisiae": 50,
    "Caenorhabditis_elegans": 30,
    "default": 30,
}

MIN_INTRON_LENGTH = MIN_INTRON_LENGTHS["default"]  # Backward compat

POLYPYRIMIDINE_WINDOW = 40  # Upstream window for acceptor scoring
POLYPYRIMIDINE_THRESHOLD = 0.5

# ==============================================================================
# 3. Restriction Enzyme Sites
# Extended from REBASE database
# ==============================================================================

RESTRICTION_ENZYMES: dict[str, str] = {
    # Common cloning enzymes
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
    # Additional common enzymes
    "NcoI": "CCATGG",
    "NheI": "GCTAGC",
    "KpnI": "GGTACC",
    "SmaI": "CCCGGG",
    "SacI": "GAGCTC",
    "SpeI": "ACTAGT",
    "ApaI": "GGGCCC",
    "ClaI": "ATCGAT",
    "EcoRV": "GATATC",
    "BglII": "AGATCT",
    "MluI": "ACGCGT",
    "AscI": "GGCGCGCC",
    "FseI": "GGCCGGCC",
    "PacI": "TTAATTAA",
    # Note: SfiI (GGCCNNNNNGGCC) uses IUPAC N wildcards and requires
    # regex-based matching. Use IUPAC_EXPAND for proper handling.
    # Simple string matching will NOT find SfiI sites correctly.
    "SfiI": "GGCCNNNNNGGCC",  # Requires wildcard matching
    "SbfI": "CCTGCAGG",
    "BsiWI": "CGTACG",
    "BsrGI": "TGTACA",
    "AgeI": "ACCGGT",
    "MfeI": "CAATTG",
}

# ==============================================================================
# 4. Nucleotide encoding for CSP solver
# ==============================================================================

BASE_MAP = {0: "A", 1: "C", 2: "G", 3: "T"}
BASE_REV = {"A": 0, "C": 1, "G": 2, "T": 3}

# IUPAC ambiguity codes
COMPLEMENT: dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "a": "t", "t": "a", "g": "c", "c": "g",
    # IUPAC ambiguity codes
    "R": "Y", "Y": "R",  # purine/pyrimidine
    "S": "S", "W": "W",  # strong/weak
    "K": "M", "M": "K",  # keto/amino
    "B": "V", "V": "B",  # not-A/not-T
    "D": "H", "H": "D",  # not-C/not-G
    "N": "N",             # any base
    # lowercase
    "r": "y", "y": "r",
    "s": "s", "w": "w",
    "k": "m", "m": "k",
    "b": "v", "v": "b",
    "d": "h", "h": "d",
    "n": "n",
}

# IUPAC base to set of concrete bases
IUPAC_EXPAND: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Supports IUPAC ambiguity codes. Raises ValueError for unknown characters.
    """
    try:
        return "".join(COMPLEMENT[base] for base in reversed(seq))
    except KeyError as e:
        raise ValueError(
            f"Unknown base '{e.args[0]}' in sequence. "
            f"Supported: A, C, G, T, N, and IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V)."
        ) from None


# ==============================================================================
# 5. Amino Acid Constants
# ==============================================================================

# Standard 20 amino acids (BLOSUM62 index order)
STANDARD_AAS: list[str] = list("ARNDCQEGHILKMFPSTWYV")

# ────────────────────────────────────────────────────────────
# BLOSUM62 Substitution Matrix (20×20, nested dict format)
# Source: Henikoff & Henikoff (1992) PNAS 89:10915
# Access: BLOSUM62[aa1][aa2] = score
# ────────────────────────────────────────────────────────────

_BLOSUM62_ROWS = [
    #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
    [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
    [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
    [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
    [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
    [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
    [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
    [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
    [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
    [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
    [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
    [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],  # K
    [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
    [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
    [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
    [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
    [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
    [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
    [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
    [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
]

_BLOSUM_INDEX = list("ARNDCQEGHILKMFPSTWYV")

BLOSUM62: dict[str, dict[str, int]] = {}
for _i, _a1 in enumerate(_BLOSUM_INDEX):
    BLOSUM62[_a1] = {}
    for _j, _a2 in enumerate(_BLOSUM_INDEX):
        BLOSUM62[_a1][_a2] = _BLOSUM62_ROWS[_i][_j]


# ────────────────────────────────────────────────────────────
# Kyte-Doolittle Hydropathy Scale
# Source: Kyte & Doolittle (1982) J Mol Biol 157:105
# ────────────────────────────────────────────────────────────

HYDROPATHY: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# Hydrophobic residues (Kyte-Doolittle > 1.0)
HYDROPHOBIC_AAS: set[str] = {"A", "I", "L", "M", "F", "W", "V"}
