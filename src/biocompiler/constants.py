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
