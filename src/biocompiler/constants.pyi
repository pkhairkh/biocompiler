"""Type stubs for biocompiler.constants — canonical biological constants."""

from __future__ import annotations

from typing import Any


# ────────────────────────────────────────────────────────────
# Genetic Code
# ────────────────────────────────────────────────────────────

CODON_TABLE: dict[str, str]
STOP_CODONS: set[str]
START_CODON: str
AA_TO_CODONS: dict[str, list[str]]


# ────────────────────────────────────────────────────────────
# Scanner Consensus Sequences
# ────────────────────────────────────────────────────────────

DONOR_CONSENSUS: str
ACCEPTOR_CONSENSUS: str
KOZAK_CONSENSUS: str
INSTABILITY_MOTIF: str
MIN_INTRON_LENGTHS: dict[str, int]
MIN_INTRON_LENGTH: int
POLYPYRIMIDINE_WINDOW: int
POLYPYRIMIDINE_THRESHOLD: float


# ────────────────────────────────────────────────────────────
# Restriction Enzymes
# ────────────────────────────────────────────────────────────

RESTRICTION_ENZYMES: dict[str, str]


# ────────────────────────────────────────────────────────────
# Nucleotide encoding
# ────────────────────────────────────────────────────────────

BASE_MAP: dict[int, str]
BASE_REV: dict[str, int]
COMPLEMENT: dict[str, str]
IUPAC_EXPAND: dict[str, str]
VALID_IUPAC_BASES: set[str]

def reverse_complement(seq: str) -> str: ...


# ────────────────────────────────────────────────────────────
# Amino Acid Constants
# ────────────────────────────────────────────────────────────

STANDARD_AAS: str
STANDARD_AAS_BLOSUM_ORDER: list[str]
BLOSUM62: dict[str, dict[str, int]]
HYDROPATHY: dict[str, float]
HYDROPHOBIC_AAS: set[str]


# ────────────────────────────────────────────────────────────
# Engine Shared Constants
# ────────────────────────────────────────────────────────────

DEFAULT_ENGINE_TIMEOUT: float
DEFAULT_BATCH_SIZE: int
DEFAULT_SOLUBILITY_WINDOW: int
DEFAULT_SOLUBILITY_SMOOTHING: int
DEFAULT_MHC_PEPTIDE_LENGTH: int
