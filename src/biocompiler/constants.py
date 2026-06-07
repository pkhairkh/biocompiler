"""Deprecated: use biocompiler.shared.constants instead."""
import warnings

warnings.warn(
    "biocompiler.constants is deprecated — use biocompiler.shared.constants instead",
    DeprecationWarning,
    stacklevel=2,
)

from biocompiler.shared.constants import *  # noqa: F401,F403

__all__ = [
    # Genetic code
    "CODON_TABLE", "STOP_CODONS", "START_CODON", "AA_TO_CODONS",
    # Scanner consensus
    "DONOR_CONSENSUS", "ACCEPTOR_CONSENSUS", "KOZAK_CONSENSUS",
    "INSTABILITY_MOTIF", "MIN_INTRON_LENGTHS", "MIN_INTRON_LENGTH",
    "POLYPYRIMIDINE_WINDOW", "POLYPYRIMIDINE_THRESHOLD",
    # Restriction enzymes
    "RESTRICTION_ENZYMES",
    # Nucleotide encoding
    "BASE_MAP", "BASE_REV", "COMPLEMENT", "IUPAC_EXPAND",
    "VALID_IUPAC_BASES",
    "reverse_complement",
    # Amino acid constants
    "STANDARD_AAS", "STANDARD_AAS_BLOSUM_ORDER", "BLOSUM62",
    "HYDROPATHY", "HYDROPHOBIC_AAS",
    # Engine shared constants
    "DEFAULT_ENGINE_TIMEOUT", "DEFAULT_BATCH_SIZE",
    "DEFAULT_SOLUBILITY_WINDOW", "DEFAULT_SOLUBILITY_SMOOTHING",
    "DEFAULT_MHC_PEPTIDE_LENGTH",
]
