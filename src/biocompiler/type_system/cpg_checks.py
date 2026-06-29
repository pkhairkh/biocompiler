"""
BioCompiler Type System — CpG Island Predicate Checks
=====================================================
CpG island detection predicate check (organism-aware GC/density windowed
scan using the fast NUMBA dinucleotide counter).

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

from .sequence_checks import (
    _count_dinucs_fast,
    _is_prokaryotic_organism,
)


logger = logging.getLogger(__name__)

def check_no_cpg_island(seq: str, window: int = 200, threshold: float = 0.6, organism: str = "") -> PredicateResult:
    """Predicate 3: No CpG islands (Obs/Exp CG ratio > threshold in any window).

    CpG island avoidance is primarily relevant for mammalian expression
    systems where CpG methylation can lead to gene silencing.  For
    prokaryotic organisms (e.g. *E. coli*), CpG islands have no known
    regulatory significance, so the check is skipped when a prokaryotic
    organism is specified.

    Optimized: Uses sliding window with O(1) updates per step instead
    of O(W) full window scan, reducing total complexity from O(N*W) to O(N).

    Args:
        seq: DNA sequence to evaluate.
        window: Sliding window size in nucleotides (default 200).
        threshold: Maximum allowed Obs/Exp CG ratio (default 0.6).
        organism: Target organism name.  If prokaryotic, the check is
            skipped and PASS is returned immediately.  If empty (default),
            the check runs as before for backward compatibility.

    Returns:
        PredicateResult with PASS/FAIL verdict.
    """
    # Skip CpG checking for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "CpG island check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoCpGIsland", True, verdict=Verdict.PASS,
            details=f"CpG island check skipped for prokaryotic organism '{organism}'",
        )

    n = len(seq)
    if n < window:
        return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                               details=f"Sequence length {n} < window size {window}")

    # Fast short-circuit: if total CG count is 0, no CpG island is possible
    total_cg = _count_dinucs_fast(seq, "CG")[0]
    if total_cg == 0:
        return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                               details="No CG dinucleotides found in sequence")

    # Pre-compute CG positions for fast lookup
    # cg_at[i] = 1 if seq[i:i+2] == "CG", else 0
    cg_at = [0] * (n - 1)
    for i in range(n - 1):
        if seq[i] == 'C' and seq[i + 1] == 'G':
            cg_at[i] = 1

    # Initialize first window
    c_count = seq[:window].count("C")
    g_count = seq[:window].count("G")
    cg_count = sum(cg_at[:window - 1])

    worst_ratio = 0.0
    worst_start = -1

    # Check first window
    expected = (c_count * g_count) / window if window > 0 else 0
    obs_exp = cg_count / expected if expected > 0 else 0.0
    if obs_exp > worst_ratio:
        worst_ratio = obs_exp
        worst_start = 0

    # Slide window — O(1) per step
    for start in range(1, n - window + 1):
        # Remove outgoing base at start-1, add incoming base at start+window-1
        outgoing = seq[start - 1]
        incoming = seq[start + window - 1]

        if outgoing == 'C':
            c_count -= 1
        elif outgoing == 'G':
            g_count -= 1
        if incoming == 'C':
            c_count += 1
        elif incoming == 'G':
            g_count += 1

        # Update CG count: remove cg_at[start-1], add cg_at[start+window-2]
        cg_count -= cg_at[start - 1]
        if start + window - 2 < n - 1:
            cg_count += cg_at[start + window - 2]

        expected = (c_count * g_count) / window if window > 0 else 0
        obs_exp = cg_count / expected if expected > 0 else 0.0
        if obs_exp > worst_ratio:
            worst_ratio = obs_exp
            worst_start = start

    if worst_ratio > threshold:
        # GC-content-aware relaxation: for GC-rich sequences, CpG density
        # is naturally higher. Only FAIL if CpG density is >2x the expected
        # density for the sequence's GC content.
        seq_gc = (seq.count("G") + seq.count("C")) / len(seq) if len(seq) > 0 else 0
        if seq_gc > _CPG_GC_RICH_THRESHOLD:
            relaxed_threshold = threshold * _CPG_DENSITY_MULTIPLIER
            if worst_ratio > relaxed_threshold:
                return PredicateResult("NoCpGIsland", False, verdict=Verdict.FAIL,
                                       details=(f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} "
                                                f"> {relaxed_threshold:.3f} (GC-rich={seq_gc:.1%}, "
                                                f"relaxed threshold)"),
                                       positions=[worst_start])
            return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                                   details=(f"CpG Obs/Exp={worst_ratio:.3f} within relaxed threshold "
                                            f"{relaxed_threshold:.3f} for GC-rich sequence ({seq_gc:.1%})"))
        return PredicateResult("NoCpGIsland", False, verdict=Verdict.FAIL,
                               details=f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} > {threshold}",
                               positions=[worst_start])
    return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                           details=f"Worst CpG Obs/Exp ratio {worst_ratio:.3f} <= {threshold}")




__all__ = [
    "check_no_cpg_island",
]
