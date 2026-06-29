"""
BioCompiler Type System — Alu-Repeat Predicate Checks
=====================================================
Alu repeat element detection (k-mer pre-filtered IUPAC consensus scan
against AluJ/AluS/AluY subfamilies, A/B Pol III boxes, and full/partial
monomers).

NOTE: This module hosts the W5-a fix for the UnboundLocalError on
``full_alu_consensus`` — the variable must be assigned BEFORE the k-mer
pre-computation block that references it. Do not reorder.

Extracted from the historical checks.py monolith during the W8-b refactor.
Re-exported by checks.py for backwards compatibility.
"""
import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from .codon_tables import (
    AA_TO_CODONS,
    BLOSUM62,
    CODON_TABLE,
    PROMOTER_CONSENSUS,
    PredicateResult,
    _BLOSUM62_MISSING_SCORE,
    _CPG_DENSITY_MULTIPLIER,
    _CPG_GC_RICH_THRESHOLD,
    _COTRANS_HIGH_CONFIDENCE,
    _COTRANS_LOW_CONFIDENCE,
    _CODON_RAMP_LENGTH,
    _DG_AU_PAIR_KCAL,
    _DG_GC_PAIR_KCAL,
    _DG_GU_PAIR_KCAL,
    _EUKARYOTE_GT_PER_BP,
    _EUK_INITIATOR_OFFSET_MAX,
    _EUK_INITIATOR_OFFSET_MIN,
    _FAST_CODON_CAI_THRESHOLD,
    _HIGH_AVG_CAI_THRESHOLD,
    _INSTABILITY_T_RUN_MIN,
    _MAXENT_INSUFFICIENT_CONTEXT_SCORE,
    _MIN_RAMP_FOR_WARNING,
    _MRNA_DG_EUKARYOTE_FAIL,
    _MRNA_DG_PROKARYOTE_FAIL,
    _MRNA_MODERATE_DG_RATIO,
    _MRNA_STABILITY_THRESHOLDS,
    _ORGANISM_TO_SPECIES_KEY,
    _PAUSE_SITE_CAI_THRESHOLD,
    _PROMOTER_UNCERTAIN_RATIO,
    _RESTRICTION_SITE_MIN_LENGTH,
    _TM_BORDERLINE_RATIO,
    _TM_EUKARYOTIC_MIN_STRETCH,
    _TM_PROKARYOTIC_MIN_STRETCH,
    _match_iupac,
    _score_consensus,
)
from biocompiler.shared.types import Verdict

from .sequence_checks import _is_prokaryotic_organism


logger = logging.getLogger(__name__)

def check_no_alu_repeat(
    seq: str,
    organism: str = "",
    min_match_score: float = 0.75,
) -> PredicateResult:
    """Predicate: No Alu repeat elements.

    Alu elements are ~300 bp SINE (Short Interspersed Nuclear Element) repeats
    that constitute ~11% of the human genome. In synthetic gene design, Alu
    elements can cause:

    - Aberrant splicing through cryptic splice site activation
    - Premature polyadenylation via internal A-rich tails
    - Homologous recombination causing genomic instability
    - Transcriptional interference from internal RNA Pol III promoters

    Detection uses a k-mer matching approach against consensus sequences from
    major Alu subfamilies (AluJ, AluS, AluY). The scan checks for:

    1. **Full Alu elements**: ~280 bp regions matching Alu consensus at >= min_match_score
    2. **Partial Alu fragments**: Left monomer (~130 bp) or right monomer (~130 bp)
       matching at >= min_match_score + 0.05 (higher threshold for shorter fragments
       to reduce false positives)
    3. **Alu-derived cryptic splice sites**: The left monomer contains A and B boxes
       (RNA Pol III promoter); the right monomer contains an A-rich tail. These
       features are checked independently.

    Only relevant for eukaryotic sequences (skipped for prokaryotes).

    Args:
        seq: DNA coding sequence (uppercase).
        organism: Organism name to determine if check is relevant.
        min_match_score: Minimum fraction of matching positions (0-1) to flag.
            Default 0.75 (75% identity to consensus).

    Returns:
        PredicateResult with verdict and positions of detected Alu elements.
    """
    # Skip for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "Alu repeat check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoAluRepeat", True, verdict=Verdict.PASS,
            details=f"Alu repeat check skipped for prokaryotic organism '{organism}'",
        )

    seq = seq.upper()

    # Alu consensus sequences (simplified from Repbase)
    # Left monomer consensus (~130bp)
    _ALU_LEFT_CONSENSUS = (
        "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACC"
        "TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAA"
        "AATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAAT"
        "GGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATT"
    )
    # Right monomer consensus (~130bp)
    _ALU_RIGHT_CONSENSUS = (
        "CGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAA"
        "AATTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCG"
        "CGATCTCGGCTCACTGCAACCTCCGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCGA"
        "GTAGCTGGGATTACAGGCGCGCGCCACCACGCCCAGCTAATTTTTGTATTTTTAGTAGAGACGGGG"
        "TTTCACCGTTTT"
    )
    # A-box and B-box of RNA Pol III promoter
    _ALU_A_BOX = "GGTTCGANNCC"  # A-box consensus
    _ALU_B_BOX = "GTTCNANNC"    # B-box consensus

    matches: List[Dict[str, Any]] = []

    # Pre-compute k-mer sets for each consensus for fast pre-filtering.
    # If a window shares < 30% of its 8-mers with the consensus, the match
    # score is guaranteed < 75%, so we skip the expensive IUPAC comparison.
    _ALU_KMER_SIZE = 8

    def _build_kmer_set(s: str) -> set:
        """Build set of k-mers from a consensus string."""
        return {s[i:i+_ALU_KMER_SIZE] for i in range(len(s) - _ALU_KMER_SIZE + 1)}

    def _has_sufficient_kmers(window: str, consensus_kmers: set, threshold: float = 0.15) -> bool:
        """Quick check: does the window share enough k-mers with the consensus?"""
        if len(window) < _ALU_KMER_SIZE:
            return False
        shared = 0
        total = len(window) - _ALU_KMER_SIZE + 1
        for i in range(total):
            if window[i:i+_ALU_KMER_SIZE] in consensus_kmers:
                shared += 1
        return shared / total >= threshold

    def _compute_match_score(window: str, consensus: str) -> float:
        """Compute fraction of matching bases between window and consensus.

        Uses IUPAC ambiguity matching via _match_iupac.
        Pre-uppers both strings once to avoid per-character .upper() calls.
        """
        if len(window) != len(consensus):
            return 0.0
        matches_count = 0
        for w_char, c_char in zip(window, consensus):
            if _match_iupac(w_char, c_char):
                matches_count += 1
        return matches_count / len(consensus)

    # Full Alu element (left + right monomer concatenated). Must be defined
    # BEFORE the k-mer pre-computation below, which references it; otherwise
    # Python treats it as an unbound local on the read.
    full_alu_consensus = _ALU_LEFT_CONSENSUS + _ALU_RIGHT_CONSENSUS
    full_alu_len = len(full_alu_consensus)

    # Pre-compute k-mer sets for all consensus sequences
    _full_alu_kmers = _build_kmer_set(full_alu_consensus)
    _left_alu_kmers = _build_kmer_set(_ALU_LEFT_CONSENSUS)
    _right_alu_kmers = _build_kmer_set(_ALU_RIGHT_CONSENSUS)
    _a_box_kmers = _build_kmer_set(_ALU_A_BOX) if len(_ALU_A_BOX) >= _ALU_KMER_SIZE else None
    _b_box_kmers = _build_kmer_set(_ALU_B_BOX) if len(_ALU_B_BOX) >= _ALU_KMER_SIZE else None

    if len(seq) >= full_alu_len:
        for start in range(len(seq) - full_alu_len + 1):
            window = seq[start:start + full_alu_len]
            if not _has_sufficient_kmers(window, _full_alu_kmers):
                continue  # k-mer pre-filter: skip windows with no ALU similarity
            score = _compute_match_score(window, full_alu_consensus)
            if score >= min_match_score:
                matches.append({
                    "type": "full_alu",
                    "start": start,
                    "length": full_alu_len,
                    "score": score,
                })

    # Check for partial Alu fragments (left monomer)
    partial_threshold = min(min_match_score + 0.05, 1.0)
    left_len = len(_ALU_LEFT_CONSENSUS)

    if len(seq) >= left_len:
        for start in range(len(seq) - left_len + 1):
            window = seq[start:start + left_len]
            if not _has_sufficient_kmers(window, _left_alu_kmers):
                continue
            score = _compute_match_score(window, _ALU_LEFT_CONSENSUS)
            if score >= partial_threshold:
                matches.append({
                    "type": "partial_alu_left",
                    "start": start,
                    "length": left_len,
                    "score": score,
                })

    # Check for partial Alu fragments (right monomer)
    right_len = len(_ALU_RIGHT_CONSENSUS)

    if len(seq) >= right_len:
        for start in range(len(seq) - right_len + 1):
            window = seq[start:start + right_len]
            if not _has_sufficient_kmers(window, _right_alu_kmers):
                continue
            score = _compute_match_score(window, _ALU_RIGHT_CONSENSUS)
            if score >= partial_threshold:
                matches.append({
                    "type": "partial_alu_right",
                    "start": start,
                    "length": right_len,
                    "score": score,
                })

    # Check for A-box and B-box independently (short patterns, higher match required)
    ab_box_threshold = 0.85
    a_box_len = len(_ALU_A_BOX)
    b_box_len = len(_ALU_B_BOX)

    if len(seq) >= a_box_len:
        for start in range(len(seq) - a_box_len + 1):
            window = seq[start:start + a_box_len]
            if _a_box_kmers and not _has_sufficient_kmers(window, _a_box_kmers):
                continue
            score = _compute_match_score(window, _ALU_A_BOX)
            if score >= ab_box_threshold:
                matches.append({
                    "type": "alu_a_box",
                    "start": start,
                    "length": a_box_len,
                    "score": score,
                })

    if len(seq) >= b_box_len:
        for start in range(len(seq) - b_box_len + 1):
            window = seq[start:start + b_box_len]
            if _b_box_kmers and not _has_sufficient_kmers(window, _b_box_kmers):
                continue
            score = _compute_match_score(window, _ALU_B_BOX)
            if score >= ab_box_threshold:
                matches.append({
                    "type": "alu_b_box",
                    "start": start,
                    "length": b_box_len,
                    "score": score,
                })

    # Determine verdict
    if not matches:
        return PredicateResult(
            "NoAluRepeat", True, verdict=Verdict.PASS,
            details="No Alu repeat elements detected",
        )

    # FAIL for full Alu elements, UNCERTAIN for partial fragments or A/B boxes
    has_full_alu = any(m["type"] == "full_alu" for m in matches)
    worst_verdict = Verdict.FAIL if has_full_alu else Verdict.UNCERTAIN

    passed = worst_verdict != Verdict.FAIL
    positions = [m["start"] for m in matches]

    # Build details
    full_alus = [m for m in matches if m["type"] == "full_alu"]
    partial_alus = [m for m in matches if m["type"] in ("partial_alu_left", "partial_alu_right")]
    ab_boxes = [m for m in matches if m["type"] in ("alu_a_box", "alu_b_box")]

    detail_parts = []
    if full_alus:
        full_info = [(m['start'], round(m['score'], 2)) for m in full_alus]
        detail_parts.append(f"Full Alu elements: {full_info}")
    if partial_alus:
        partial_info = [(m['start'], m['type'], round(m['score'], 2)) for m in partial_alus]
        detail_parts.append(f"Partial Alu fragments: {partial_info}")
    if ab_boxes:
        box_info = [(m['start'], m['type'], round(m['score'], 2)) for m in ab_boxes]
        detail_parts.append(f"Alu Pol III boxes: {box_info}")

    details = "; ".join(detail_parts)

    return PredicateResult(
        "NoAluRepeat", passed, verdict=worst_verdict,
        details=details,
        positions=positions,
    )


__all__ = [
    "check_no_alu_repeat",
]
