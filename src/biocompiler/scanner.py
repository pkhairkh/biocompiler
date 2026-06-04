"""
BioCompiler Scanner — Multi-DFA Motif Detection

Production-grade scanner with:
- Start/stop codons scanned in ALL 3 reading frames (not just frame 0)
- Restriction sites checked on BOTH strands (forward + reverse complement)
- MaxEntScan-based splice site scoring (not constant 5.0 for donors)
- Kozak consensus scoring with position weights (not exact string match)
- Logging instead of print
"""

import logging
from .constants import (
    DONOR_CONSENSUS, ACCEPTOR_CONSENSUS, INSTABILITY_MOTIF,
    RESTRICTION_ENZYMES, POLYPYRIMIDINE_WINDOW,
    IUPAC_EXPAND,
    reverse_complement,
)
from .types import Token
from .exceptions import InvalidSequenceError
from .maxentscan import score_donor, score_acceptor

logger = logging.getLogger(__name__)

# Kozak consensus position weights for scoring
# Reference: Kozak M. (1987) "An analysis of 5'-noncoding sequences from 699 vertebrate mRNAs"
# Nucleic Acids Res. 15(20):8125-48
# Position -3 (A/G): weight 0.3
# Position -2 (C): weight 0.2
# Position -1 (C): weight 0.4 (most conserved non-ATG position)
# Position +4 (G): weight 0.1
# GCCACCATGG is the "optimal" Kozak consensus
KOZAK_POSITION_WEIGHTS: dict[int, dict[str, float]] = {
    -3: {"A": 1.0, "G": 0.8, "C": 0.2, "T": 0.1},  # A/G strongly preferred
    -2: {"C": 1.0, "A": 0.4, "G": 0.3, "T": 0.2},   # C preferred
    -1: {"C": 1.0, "A": 0.2, "G": 0.3, "T": 0.1},   # C most conserved
    +4: {"G": 1.0, "A": 0.3, "C": 0.2, "T": 0.2},   # G moderately preferred
}

# Minimum MaxEntScan score for a splice site to be considered functional
# Based on Yeo & Burge (2004): scores >3.0 are likely functional splice sites
SPLICE_DONOR_MIN_SCORE = 3.0
SPLICE_ACCEPTOR_MIN_SCORE = 3.0

# Fallback score assigned to splice donors when MaxEntScan is disabled
DONOR_FALLBACK_SCORE: float = 5.0

# Minimum C+T fraction in the upstream polypyrimidine tract for acceptor calling
POLYPYRIMIDINE_MIN_FRACTION: float = 0.5

# Multiplier that converts polypyrimidine fraction [0,1] to an acceptor score
ACCEPTOR_SCORE_MULTIPLIER: float = 10.0

# Minimum Kozak consensus score for a site to be reported as a kozak token
KOZAK_REPORT_THRESHOLD: float = 0.7


def _iupac_match(seq: str, pattern: str) -> bool:
    """
    Check if a DNA sequence matches an IUPAC pattern.
    
    Supports ambiguity codes: R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC,
    B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT.
    
    Args:
        seq: concrete DNA sequence (only ACGT)
        pattern: IUPAC pattern (may contain ambiguity codes)
    
    Returns:
        True if seq matches pattern
    """
    if len(seq) != len(pattern):
        return False
    for s_base, p_base in zip(seq, pattern):
        allowed = IUPAC_EXPAND.get(p_base.upper(), p_base.upper())
        if s_base not in allowed:
            return False
    return True


def validate_dna_sequence(seq: str) -> str:
    """
    Validate and normalize a DNA sequence.

    Converts to uppercase and checks that every character is in {A, C, G, T, N}.

    Args:
        seq: raw DNA sequence (case-insensitive)

    Returns:
        Uppercased, validated DNA string.

    Raises:
        InvalidSequenceError: if any character is not A/C/G/T/N.
    """
    seq = seq.upper()
    valid = set("ACGTN")
    invalid = set(seq) - valid
    if invalid:
        raise InvalidSequenceError(seq, invalid)
    return seq


def gc_content(seq: str) -> float:
    """
    Compute GC content as a fraction in [0.0, 1.0].

    Deterministic: depends only on the input sequence.

    Args:
        seq: DNA sequence (case-insensitive)

    Returns:
        GC fraction rounded to 4 decimal places.  Returns 0.0 for empty input.
    """
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return round(gc / len(seq), 4)


def _score_kozak(seq: str, atg_pos: int) -> float:
    """
    Score the Kozak consensus context around an ATG at the given position.

    Uses position-specific weights based on Kozak (1987) conservation data.
    Score range: 0.0 (no consensus) to 1.0 (perfect consensus GCCACCATGG).

    Args:
        seq: DNA sequence (already uppercased)
        atg_pos: position of the 'A' in ATG

    Returns:
        Kozak score in [0.0, 1.0]
    """
    score = 0.0
    total_weight = 0.0
    for offset, weights in KOZAK_POSITION_WEIGHTS.items():
        pos = atg_pos + offset
        if 0 <= pos < len(seq):
            base = seq[pos]
            weight = sum(weights.values())  # max possible for this position
            total_weight += weight
            score += weights.get(base, 0.0) * weight
    return round(score / total_weight, 4) if total_weight > 0 else 0.0


def scan_sequence(
    seq: str,
    restriction_enzymes: list[str] | None = None,
    scan_all_frames: bool = True,
    use_maxentscan: bool = True,
    donor_threshold: float = SPLICE_DONOR_MIN_SCORE,
    acceptor_threshold: float = SPLICE_ACCEPTOR_MIN_SCORE,
) -> list[Token]:
    """
    Scan a nucleotide sequence for biological motifs.

    This is a deterministic DFA-based scanner. Each position is examined
    for all motif types simultaneously.

    Args:
        seq: DNA sequence to scan
        restriction_enzymes: list of enzyme names to check for
        scan_all_frames: if True, scan start/stop codons in all 3 reading frames
        use_maxentscan: if True, use MaxEntScan scoring for splice sites;
                        if False, use simple consensus + polypyrimidine scoring
        donor_threshold: minimum MaxEntScan score for donor site reporting
        acceptor_threshold: minimum MaxEntScan score for acceptor site reporting

    Returns:
        Ordered list of tokens sorted by (position, element_type).
    """
    seq = validate_dna_sequence(seq)
    if not seq:
        return []
    tokens: list[Token] = []

    # --- Splice donor sites (GT) with MaxEntScan scoring ---
    for i in range(len(seq) - 1):
        if seq[i:i+2] == DONOR_CONSENSUS:
            if use_maxentscan:
                mes_score = score_donor(seq, i)
                if mes_score >= donor_threshold:
                    tokens.append(Token(i, "splice_donor", seq[i:i+2], mes_score))
            else:
                tokens.append(Token(i, "splice_donor", seq[i:i+2], DONOR_FALLBACK_SCORE))

    # --- Splice acceptor sites (AG) with MaxEntScan scoring ---
    for i in range(len(seq) - 1):
        if seq[i:i+2] == ACCEPTOR_CONSENSUS:
            if use_maxentscan:
                mes_score = score_acceptor(seq, i)
                if mes_score >= acceptor_threshold:
                    tokens.append(Token(i, "splice_acceptor", seq[i:i+2], mes_score))
            else:
                # Fallback: polypyrimidine tract scoring
                upstream = seq[max(0, i - POLYPYRIMIDINE_WINDOW):i]
                ct_count = upstream.count('C') + upstream.count('T')
                score = ct_count / max(len(upstream), 1)
                if score > POLYPYRIMIDINE_MIN_FRACTION:
                    tokens.append(Token(i, "splice_acceptor", seq[i:i+2], score * ACCEPTOR_SCORE_MULTIPLIER))

    # --- Start codons in ALL reading frames ---
    if scan_all_frames:
        for frame in range(3):
            for i in range(frame, len(seq) - 2, 3):
                if seq[i:i+3] == "ATG":
                    kozak_score = _score_kozak(seq, i)
                    tokens.append(Token(i, "start_codon", "ATG", kozak_score, frame=frame))
    else:
        for i in range(0, len(seq) - 2, 3):
            if seq[i:i+3] == "ATG":
                kozak_score = _score_kozak(seq, i)
                tokens.append(Token(i, "start_codon", "ATG", kozak_score, frame=0))

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

    # --- Kozak consensus (weighted scoring, not exact match) ---
    for i in range(len(seq) - 2):
        if seq[i:i+3] == "ATG":
            kozak_score = _score_kozak(seq, i)
            if kozak_score >= KOZAK_REPORT_THRESHOLD:  # Only report strong Kozak contexts
                tokens.append(Token(i, "kozak", seq[max(0,i-3):i+5], kozak_score))

    # --- Instability motifs (ATTTA) ---
    for i in range(len(seq) - 4):
        if seq[i:i+5] == INSTABILITY_MOTIF:
            tokens.append(Token(i, "instability_motif", seq[i:i+5], 1.0))

    # --- Restriction enzyme sites (BOTH strands, IUPAC-aware) ---
    if restriction_enzymes:
        for enz_name in restriction_enzymes:
            if enz_name in RESTRICTION_ENZYMES:
                site = RESTRICTION_ENZYMES[enz_name]
                has_iupac = any(b not in "ACGT" for b in site.upper())
                # Forward strand
                for i in range(len(seq) - len(site) + 1):
                    window = seq[i:i+len(site)]
                    if has_iupac:
                        if _iupac_match(window, site):
                            tokens.append(Token(i, "restriction_site", window, 1.0, strand="+"))
                    else:
                        if window == site:
                            tokens.append(Token(i, "restriction_site", site, 1.0, strand="+"))
                # Reverse complement strand (only for non-IUPAC sites; IUPAC RC is complex)
                if not has_iupac:
                    site_rc = reverse_complement(site)
                    if site_rc != site:  # Avoid double-counting palindromes
                        for i in range(len(seq) - len(site_rc) + 1):
                            if seq[i:i+len(site_rc)] == site_rc:
                                tokens.append(Token(i, "restriction_site", site_rc, 1.0, strand="-"))

    tokens.sort(key=lambda t: (t.position, t.element_type))
    logger.debug("Scanned %d nt sequence, found %d tokens", len(seq), len(tokens))
    return tokens
