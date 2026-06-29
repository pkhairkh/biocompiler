"""
BioCompiler Type System — Splice-Site Predicate Checks
======================================================
Splice-site and GT dinucleotide predicate checks: cryptic splice sites
(MaxEntScan-based), GT dinucleotide (hard), avoidable GT (organism-aware),
and soft GT dinucleotide (organism-aware with certainty levels).

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
    _compute_max_gt_count,
    find_cross_codon_gt,
)


logger = logging.getLogger(__name__)

def check_no_cryptic_splice(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0, organism: str = "") -> PredicateResult:
    """Predicate 2: No cryptic splice sites (dual-threshold PASS/UNCERTAIN/FAIL).

    Uses the proper MaxEntScan log-odds scoring model (Yeo & Burge 2004) from
    the maxentscan module to evaluate both donor (GT) and acceptor (AG) splice
    sites.  Numeric scores are converted to SpliceVerdict thresholds:
      - score < low_thresh  -> PASS
      - low_thresh <= score < high_thresh -> UNCERTAIN
      - score >= high_thresh -> FAIL

    Organism-specific thresholds:
    - Prokaryotes: auto-PASS (no splicing in prokaryotes)
    - Eukaryotes: high_thresh=8.0 (stricter than default of 6.0 to reduce
      false positives from common coding-sequence GT/AG dinucleotides)

    Skipped for prokaryotic organisms (splice sites are a eukaryote-specific
    concern).
    """
    # Skip cryptic splice check for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "Cryptic splice check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoCrypticSplice", True, verdict=Verdict.PASS,
            details=f"Cryptic splice check skipped for prokaryotic organism '{organism}'",
        )

    # For eukaryotes, use stricter high_thresh of 8.0 to reduce false positives
    effective_high_thresh = high_thresh
    if organism and not _is_prokaryotic_organism(organism):
        effective_high_thresh = max(high_thresh, 8.0)

    from biocompiler.sequence.maxentscan import score_donor, score_acceptor

    seq = seq.upper()
    max_score = _MAXENT_INSUFFICIENT_CONTEXT_SCORE
    worst_pos = -1
    worst_verdict = Verdict.PASS

    # Scan donor sites (GT dinucleotides) — skip if no GT in sequence
    if "GT" not in seq:
        pass  # no donor sites possible
    else:
      for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            score = score_donor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0
            if score < low_thresh:
                v = Verdict.PASS
            elif score < effective_high_thresh:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.FAIL
            if score > max_score:
                max_score = score
                worst_pos = i
                worst_verdict = v

    # Scan acceptor sites (AG dinucleotides) — skip if no AG in sequence
    if "AG" not in seq:
        pass  # no acceptor sites possible
    else:
      for i in range(len(seq) - 1):
        if seq[i:i+2] == "AG":
            score = score_acceptor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0
            if score < low_thresh:
                v = Verdict.PASS
            elif score < effective_high_thresh:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.FAIL
            if score > max_score:
                max_score = score
                worst_pos = i
                worst_verdict = v

    if worst_pos < 0:
        return PredicateResult("NoCrypticSplice", True, verdict=Verdict.PASS,
                               details="No splice dinucleotides found")

    passed = worst_verdict != Verdict.FAIL
    return PredicateResult("NoCrypticSplice", passed, verdict=worst_verdict,
                           details=f"Worst splice score {max_score:.2f} at pos {worst_pos}",
                           positions=[worst_pos] if worst_pos >= 0 else [])


def check_no_gt_dinucleotide(seq: str) -> PredicateResult:
    """Predicate 5: No GT dinucleotides (5' splice donor mimic), including cross-codon.

    This is the STRICT version — any GT fails the predicate.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for fast count-based short-circuit; falls back to pure-Python otherwise.
    """
    # Fast short-circuit: if GT count is 0, no need to enumerate positions
    gt_count = _count_dinucs_fast(seq, "GT")[0]
    if gt_count == 0:
        return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS, details="No GT dinucleotides found")
    # Need positions for the result — enumerate only when count > 0
    positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
    return PredicateResult("NoGTDinucleotide", False, verdict=Verdict.FAIL,
                           details=f"GT dinucleotides at {positions}",
                           positions=positions)


def check_no_avoidable_gt(seq: str, organism: str = "") -> PredicateResult:
    """Predicate 5 (relaxed): No avoidable GT dinucleotides.

    A GT is "unavoidable" if ALL synonymous codons for that amino acid
    also contain GT or create a cross-codon GT.  This predicate PASSES
    if every remaining GT in the sequence is unavoidable — i.e., there
    is no synonymous substitution that could remove it.

    Specifically:
    - Within-codon GT: unavoidable if every synonymous codon for the AA
      also contains "GT" (e.g., Valine GTN where all 4 codons start with GT)
    - Cross-codon GT: unavoidable if no combination of synonymous codons
      for the two adjacent AAs eliminates the boundary GT

    Skipped for prokaryotic organisms (GT splice donor sites are a
    eukaryote-specific concern).
    """
    # Skip GT dinucleotide check for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "GT dinucleotide check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details=f"GT dinucleotide check skipped for prokaryotic organism '{organism}'",
        )

    # Fast short-circuit: if GT count is 0, no need to enumerate positions
    gt_count = _count_dinucs_fast(seq, "GT")[0]
    if gt_count == 0:
        return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS, details="No GT dinucleotides found")

    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]

    avoidable_positions = []
    unavoidable_positions = []

    for pos in gt_positions:
        codon_idx = pos // 3  # which codon does position 'pos' fall in?
        codon_start = codon_idx * 3
        next_codon_start = codon_start + 3

        # Determine whether this GT is within a single codon or crosses a boundary
        if pos + 1 < next_codon_start:
            # Within-codon GT (both bases in the same codon)
            codon = seq[codon_start:codon_start + 3]
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                unavoidable_positions.append(pos)
                continue

            # Check if any synonymous codon avoids GT
            has_avoidable = False
            for alt in AA_TO_CODONS.get(aa, []):
                if "GT" not in alt:
                    # Also check this alt does not create cross-codon GT
                    # with the previous codon's last base
                    if codon_start > 0:
                        prev_base = seq[codon_start - 1]
                        if prev_base + alt[0] == "GT":
                            continue  # would create cross-codon GT
                    # And check it does not create GT with the next codon's first base
                    if next_codon_start + 3 <= len(seq):
                        next_base = seq[next_codon_start]
                        if alt[-1] + next_base == "GT":
                            continue  # would create cross-codon GT
                    has_avoidable = True
                    break

            if has_avoidable:
                avoidable_positions.append(pos)
            else:
                unavoidable_positions.append(pos)
        else:
            # Cross-codon GT (pos is last base of one codon, pos+1 is first of next)
            prev_codon_start = codon_start  # codon containing 'pos'
            curr_codon_start = next_codon_start  # codon containing 'pos+1'

            if curr_codon_start + 3 > len(seq):
                unavoidable_positions.append(pos)
                continue

            prev_codon = seq[prev_codon_start:prev_codon_start + 3]
            curr_codon = seq[curr_codon_start:curr_codon_start + 3]
            prev_aa = CODON_TABLE.get(prev_codon)
            curr_aa = CODON_TABLE.get(curr_codon)

            if prev_aa is None or curr_aa is None:
                unavoidable_positions.append(pos)
                continue

            # If one side is a stop codon, we can only try changing the other side
            if prev_aa == "*" and curr_aa == "*":
                unavoidable_positions.append(pos)
                continue

            has_avoidable = False

            if prev_aa == "*":
                for c_alt in AA_TO_CODONS.get(curr_aa, [curr_codon]):
                    if prev_codon[-1] + c_alt[0] != "GT":
                        has_avoidable = True
                        break
            elif curr_aa == "*":
                for p_alt in AA_TO_CODONS.get(prev_aa, [prev_codon]):
                    if p_alt[-1] + curr_codon[0] != "GT":
                        has_avoidable = True
                        break
            else:
                prev_alts = AA_TO_CODONS.get(prev_aa, [prev_codon])
                curr_alts = AA_TO_CODONS.get(curr_aa, [curr_codon])

                for p_alt in prev_alts:
                    for c_alt in curr_alts:
                        if p_alt[-1] + c_alt[0] != "GT":
                            has_avoidable = True
                            break
                    if has_avoidable:
                        break

            if has_avoidable:
                avoidable_positions.append(pos)
            else:
                unavoidable_positions.append(pos)

    if avoidable_positions:
        return PredicateResult("NoGTDinucleotide", False, verdict=Verdict.FAIL,
                               details=(f"Avoidable GT dinucleotides at {avoidable_positions}; "
                                        f"unavoidable at {unavoidable_positions}"),
                               positions=avoidable_positions)
    return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS,
                           details=(f"All {len(unavoidable_positions)} GT dinucleotides are "
                                    f"unavoidable (no synonymous substitution can remove them)"),
                           positions=unavoidable_positions)


def check_no_gt_dinucleotide_soft(
    seq: str,
    organism: str = "",
    max_gt_count: int | None = None,
) -> PredicateResult:
    """Predicate 5 (soft): No GT dinucleotides with eukaryote-aware tolerance.

    This is the organism-aware soft-constraint version of the GT dinucleotide
    check, designed for eukaryotic gene optimization where destroying CAI to
    eliminate every GT is counter-productive.

    Evaluation semantics:

    - **Prokaryotes**: Hard constraint — any GT dinucleotide is a FAIL.

    - **Eukaryotes**: Soft constraint — GTs are reported but the predicate
      uses ``LIKELY_FAIL`` (not ``FAIL``) when GTs exceed ``max_gt_count``,
      indicating a soft violation that should not block the optimization.
      The predicate PASSES (``PASS``) if GT count <= ``max_gt_count``.

    - ``max_gt_count``: If not provided, auto-computed from sequence length
      using :func:`_compute_max_gt_count`.

    Args:
        seq: DNA sequence to evaluate.
        organism: Target organism name. If prokaryotic, any GT is FAIL.
        max_gt_count: Maximum GT count before triggering SOFT_FAIL for
            eukaryotes. If None, auto-computed from sequence length.

    Returns:
        PredicateResult with verdict.
    """
    seq = seq.upper()
    # Fast short-circuit: use NUMBA kernel for count, skip position enumeration if 0
    gt_count_fast = _count_dinucs_fast(seq, "GT")[0]
    if gt_count_fast == 0:
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details="No GT dinucleotides found",
        )

    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i + 2] == "GT"]

    if not gt_positions:
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details="No GT dinucleotides found",
        )

    # Compute max_gt_count if not provided
    if max_gt_count is None:
        max_gt_count = _compute_max_gt_count(len(seq), organism)

    # Count in-codon vs cross-codon GTs for reporting
    in_codon_gt = []
    cross_codon_gt = []
    for pos in gt_positions:
        codon_of_g = pos // 3
        codon_of_t = (pos + 1) // 3
        if codon_of_g == codon_of_t:
            in_codon_gt.append(pos)
        else:
            cross_codon_gt.append(pos)

    gt_count = len(gt_positions)

    # Prokaryotes: hard constraint (FAIL for any GT)
    if organism and _is_prokaryotic_organism(organism):
        return PredicateResult(
            "NoGTDinucleotide", False, verdict=Verdict.FAIL,
            details=(
                f"GT dinucleotides: {gt_count} "
                f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
                f"Hard constraint for prokaryotes: max_gt_count=0."
            ),
            positions=gt_positions,
        )

    # Eukaryotes: soft constraint
    if gt_count <= max_gt_count:
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details=(
                f"GT dinucleotides: {gt_count} <= max_gt_count={max_gt_count} "
                f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
                f"Acceptable for eukaryotes — in-codon GTs from optimal codons "
                f"are biologically common."
            ),
            positions=gt_positions,
        )

    # GT count exceeds tolerance — soft fail
    return PredicateResult(
        "NoGTDinucleotide", True, verdict=Verdict.LIKELY_FAIL,
        details=(
            f"GT dinucleotides: {gt_count} > max_gt_count={max_gt_count} "
            f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
            f"Soft constraint for eukaryotes: in-codon GTs from optimal codons "
            f"are acceptable (CAI > GT avoidance). Consider if these GTs form "
            f"strong cryptic splice donors (MaxEntScan score >= threshold)."
        ),
        positions=gt_positions,
    )




__all__ = [
    "check_no_cryptic_splice",
    "check_no_gt_dinucleotide",
    "check_no_avoidable_gt",
    "check_no_gt_dinucleotide_soft",
]
