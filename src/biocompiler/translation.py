"""
BioCompiler Translation Engine — Deterministic FST

FIXES from toy model:
- Multi-organism CAI support
- Handles partial codons at end of sequence
- Detailed translation metadata
"""

from __future__ import annotations

import logging
import math
from typing import TypedDict

from .constants import CODON_TABLE, START_CODON, reverse_complement
from .scanner import validate_dna_sequence
from .organisms import (
    CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS,
)
from .exceptions import UnsupportedOrganismError

__all__ = [
    "translate",
    "compute_cai",
    "find_orfs",
    "ORFResult",
    "DEFAULT_MIN_ORF_LENGTH_AA",
]

logger = logging.getLogger(__name__)

# --- Named constants (replacing magic numbers) ---

# Number of nucleotides per codon
_CODON_LENGTH: int = 3

# Default minimum ORF length in amino acids for find_orfs()
DEFAULT_MIN_ORF_LENGTH_AA: int = 30

# Decimal places for rounding CAI values
_CAI_ROUND_PRECISION: int = 4

# Epsilon floor for zero-adaptiveness codons in CAI computation
_ZERO_ADAPTIVENESS_EPSILON: float = 1e-10


class ORFResult(TypedDict):
    """Structured type for ORF results returned by find_orfs."""

    start: int
    end: int
    frame: int
    strand: str
    protein: str
    length: int


def translate(sequence: str, to_stop: bool = True) -> str:
    """
    Translate a DNA sequence to a protein sequence via the standard genetic code.

    This is a deterministic FST: codon → amino acid mapping.
    Translates from frame 0 until a stop codon or end of sequence.

    Args:
        sequence: DNA coding sequence
        to_stop: if True, stop at first stop codon; if False, include stops as '*'

    Returns:
        Protein sequence as single-letter amino acid codes.
    """
    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return ""

    # Warn about partial codons
    if len(sequence) % _CODON_LENGTH != 0:
        logger.warning(
            "Sequence length %d is not a multiple of %d; last %d base(s) will be ignored",
            len(sequence), _CODON_LENGTH, len(sequence) % _CODON_LENGTH,
        )

    protein: list[str] = []
    for i in range(0, len(sequence) - (_CODON_LENGTH - 1), _CODON_LENGTH):
        codon = sequence[i:i + _CODON_LENGTH]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            logger.warning("Unknown codon '%s' at position %d — mapping to 'X'", codon, i)
            aa = "X"
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def compute_cai(sequence: str, organism: str = "Homo_sapiens") -> float:
    """
    Compute Codon Adaptation Index (CAI) for a coding sequence.

    CAI = geometric mean of relative adaptiveness values of codons used.
    This is a DETERMINISTIC computation: same sequence → same CAI.

    Args:
        sequence: DNA coding sequence (starts with ATG)
        organism: organism name (must be in SUPPORTED_ORGANISMS)

    Returns:
        CAI value in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.

    Raises:
        UnsupportedOrganismError: if organism is not supported
    """
    sequence = validate_dna_sequence(sequence)
    if not sequence:
        return 0.0
    if organism not in SUPPORTED_ORGANISMS:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)

    adaptiveness = CODON_ADAPTIVENESS_TABLES[organism]
    ratios: list[float] = []

    for i in range(0, len(sequence) - (_CODON_LENGTH - 1), _CODON_LENGTH):
        codon = sequence[i:i + _CODON_LENGTH]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*" or aa == "M":
            continue
        w = adaptiveness.get(codon, 0.0)
        if w <= 0:
            w = _ZERO_ADAPTIVENESS_EPSILON  # Floor for zero-adaptiveness codons
        ratios.append(w)

    if not ratios:
        return 0.0

    # Geometric mean via log-based computation
    log_sum: float = sum(math.log(r) for r in ratios)
    cai = math.exp(log_sum / len(ratios))
    return round(cai, _CAI_ROUND_PRECISION)


def find_orfs(sequence: str, min_length_aa: int = DEFAULT_MIN_ORF_LENGTH_AA) -> list[ORFResult]:
    """
    Find all Open Reading Frames in all 6 frames (3 forward + 3 reverse complement).

    This is a production feature — the toy model didn't do ORF finding at all.

    Args:
        sequence: DNA sequence
        min_length_aa: minimum ORF length in amino acids

    Returns:
        List of ORF dicts with keys: start, end, frame, strand, protein, length
    """
    sequence = validate_dna_sequence(sequence)
    orfs: list[ORFResult] = []

    def _find_forward_orfs(seq: str, strand: str) -> list[ORFResult]:
        """Scan a single-strand sequence for ORFs in all 3 forward frames."""
        found: list[ORFResult] = []
        for frame in range(_CODON_LENGTH):
            i: int = frame
            while i < len(seq) - (_CODON_LENGTH - 1):
                codon = seq[i:i + _CODON_LENGTH]
                if codon == START_CODON:
                    # Found start codon — translate until stop
                    protein: list[str] = []
                    orf_start: int = i
                    j: int = i
                    while j < len(seq) - (_CODON_LENGTH - 1):
                        c = seq[j:j + _CODON_LENGTH]
                        aa = CODON_TABLE.get(c)
                        if aa is None:
                            logger.warning(
                                "Unknown codon '%s' at position %d in ORF scan — mapping to 'X'",
                                c, j,
                            )
                            aa = "X"
                        if aa == "*":
                            break
                        protein.append(aa)
                        j += _CODON_LENGTH
                    if len(protein) >= min_length_aa:
                        found.append({
                            "start": orf_start,
                            "end": j + _CODON_LENGTH,
                            "frame": frame,
                            "strand": strand,
                            "protein": "".join(protein),
                            "length": len(protein),
                        })
                    i = j + _CODON_LENGTH
                else:
                    i += _CODON_LENGTH
        return found

    # Forward strand
    orfs.extend(_find_forward_orfs(sequence, "+"))
    # Reverse complement strand
    rc_seq = reverse_complement(sequence)
    rc_orfs = _find_forward_orfs(rc_seq, "-")
    # Convert RC positions back to original coordinates
    seq_len = len(sequence)
    for orf in rc_orfs:
        orig_end = seq_len - orf["start"]
        orig_start = seq_len - orf["end"]
        orf["start"] = orig_start
        orf["end"] = orig_end
    orfs.extend(rc_orfs)

    orfs.sort(key=lambda o: o["start"])
    return orfs
