"""
BioCompiler Scanner — Multi-DFA Motif Detection

FIXES from toy model:
- Start/stop codons scanned in ALL 3 reading frames (not just frame 0)
- Restriction sites checked on BOTH strands (forward + reverse complement)
- Proper token frame annotation
- Logging instead of print
"""

import logging
from .constants import (
    DONOR_CONSENSUS, ACCEPTOR_CONSENSUS, KOZAK_CONSENSUS, INSTABILITY_MOTIF,
    RESTRICTION_ENZYMES, POLYPYRIMIDINE_WINDOW, POLYPYRIMIDINE_THRESHOLD,
    reverse_complement,
)
from .types import Token
from .exceptions import InvalidSequenceError

logger = logging.getLogger(__name__)


def validate_dna_sequence(seq: str) -> str:
    """Validate and normalize a DNA sequence. Raises InvalidSequenceError for bad input."""
    seq = seq.upper()
    valid = set("ACGTN")
    invalid = set(seq) - valid
    if invalid:
        raise InvalidSequenceError(seq, invalid)
    return seq


def gc_content(seq: str) -> float:
    """Compute GC content as a fraction [0.0, 1.0]. Deterministic."""
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return round(gc / len(seq), 4)


def scan_sequence(seq: str, restriction_enzymes: list[str] | None = None, scan_all_frames: bool = True) -> list[Token]:
    """
    Scan a nucleotide sequence for biological motifs.

    This is a deterministic DFA-based scanner. Each position is examined
    for all motif types simultaneously.

    Args:
        seq: DNA sequence to scan
        restriction_enzymes: list of enzyme names to check for
        scan_all_frames: if True, scan start/stop codons in all 3 reading frames

    Returns:
        Ordered list of tokens sorted by (position, element_type).
    """
    seq = validate_dna_sequence(seq)
    if not seq:
        return []
    tokens: list[Token] = []

    # --- Splice donor sites (GT) ---
    for i in range(len(seq) - 1):
        if seq[i:i+2] == DONOR_CONSENSUS:
            tokens.append(Token(i, "splice_donor", seq[i:i+2], 5.0))

    # --- Splice acceptor sites (AG) with polypyrimidine tract scoring ---
    for i in range(len(seq) - 1):
        if seq[i:i+2] == ACCEPTOR_CONSENSUS:
            upstream = seq[max(0, i - POLYPYRIMIDINE_WINDOW):i]
            ct_count = upstream.count('C') + upstream.count('T')
            score = ct_count / max(len(upstream), 1)
            if score > POLYPYRIMIDINE_THRESHOLD:
                tokens.append(Token(i, "splice_acceptor", seq[i:i+2], score * 10.0))

    # --- Start codons in ALL reading frames ---
    if scan_all_frames:
        for frame in range(3):
            for i in range(frame, len(seq) - 2, 3):
                if seq[i:i+3] == "ATG":
                    tokens.append(Token(i, "start_codon", "ATG", 1.0, frame=frame))
    else:
        for i in range(0, len(seq) - 2, 3):
            if seq[i:i+3] == "ATG":
                tokens.append(Token(i, "start_codon", "ATG", 1.0, frame=0))

    # --- Stop codons in ALL reading frames ---
    if scan_all_frames:
        for frame in range(3):
            for i in range(frame, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if codon in ("TAA", "TAG", "TGA"):
                    tokens.append(Token(i, "stop_codon", codon, 1.0, frame=frame))
    else:
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in ("TAA", "TAG", "TGA"):
                tokens.append(Token(i, "stop_codon", codon, 1.0, frame=0))

    # --- Kozak consensus ---
    for i in range(len(seq) - len(KOZAK_CONSENSUS) + 1):
        if seq[i:i+len(KOZAK_CONSENSUS)] == KOZAK_CONSENSUS:
            tokens.append(Token(i, "kozak", KOZAK_CONSENSUS, 1.0))

    # --- Instability motifs (ATTTA) ---
    for i in range(len(seq) - 4):
        if seq[i:i+5] == INSTABILITY_MOTIF:
            tokens.append(Token(i, "instability_motif", seq[i:i+5], 1.0))

    # --- Restriction enzyme sites (BOTH strands) ---
    if restriction_enzymes:
        for enz_name in restriction_enzymes:
            if enz_name in RESTRICTION_ENZYMES:
                site = RESTRICTION_ENZYMES[enz_name]
                site_rc = reverse_complement(site)
                # Forward strand
                for i in range(len(seq) - len(site) + 1):
                    if seq[i:i+len(site)] == site:
                        tokens.append(Token(i, "restriction_site", site, 1.0, strand="+"))
                # Reverse complement strand
                if site_rc != site:  # Avoid double-counting palindromes
                    for i in range(len(seq) - len(site_rc) + 1):
                        if seq[i:i+len(site_rc)] == site_rc:
                            tokens.append(Token(i, "restriction_site", site_rc, 1.0, strand="-"))

    tokens.sort(key=lambda t: (t.position, t.element_type))
    logger.debug("Scanned %d nt sequence, found %d tokens", len(seq), len(tokens))
    return tokens
