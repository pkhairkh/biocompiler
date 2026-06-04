"""
BioCompiler Translation Engine — Deterministic FST

FIXES from toy model:
- Multi-organism CAI support
- Handles partial codons at end of sequence
- Detailed translation metadata
"""

import logging
import math
from .constants import CODON_TABLE, START_CODON, reverse_complement
from .scanner import validate_dna_sequence
from .organisms import (
    CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS,
)
from .exceptions import UnsupportedOrganismError

logger = logging.getLogger(__name__)

# Epsilon floor for zero-adaptiveness codons in CAI computation
_ZERO_ADAPTIVENESS_EPSILON: float = 1e-10


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
    if len(sequence) % 3 != 0:
        logger.warning(
            "Sequence length %d is not a multiple of 3; last %d base(s) will be ignored",
            len(sequence), len(sequence) % 3,
        )

    protein: list[str] = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
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

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
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
    log_sum = sum(math.log(r) for r in ratios)
    cai = math.exp(log_sum / len(ratios))
    return round(cai, 4)


def find_orfs(sequence: str, min_length_aa: int = 30) -> list[dict]:
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
    orfs: list[dict] = []

    def _find_forward_orfs(seq: str, strand: str) -> list[dict]:
        """Scan a single-strand sequence for ORFs in all 3 forward frames."""
        found: list[dict] = []
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                codon = seq[i:i+3]
                if codon == START_CODON:
                    # Found start codon — translate until stop
                    protein = []
                    orf_start = i
                    j = i
                    while j < len(seq) - 2:
                        c = seq[j:j+3]
                        aa = CODON_TABLE.get(c, "X")
                        if aa == "*":
                            break
                        protein.append(aa)
                        j += 3
                    if len(protein) >= min_length_aa:
                        found.append({
                            "start": orf_start,
                            "end": j + 3,
                            "frame": frame,
                            "strand": strand,
                            "protein": "".join(protein),
                            "length": len(protein),
                        })
                    i = j + 3
                else:
                    i += 3
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
