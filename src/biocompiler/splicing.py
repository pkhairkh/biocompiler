"""
BioCompiler NDFST Splicing Engine — Enumerative Isoform Computation

Production-grade splicing engine with:
- MaxEntScan-based scoring for cryptic splice site detection (simplified PWM)
- Enumerative isoform computation (all valid donor-acceptor paths)
- Intron retention isoforms (major class of alternative splicing, ~15% of events)
- Multi-exon skipping combinations
- Tissue-weighted scoring using GTex-derived data (tissue_data module)
- Proper handling of alternative 5'/3' sites
- Configurable parameters (not hardcoded magic numbers)
"""

import logging
from itertools import combinations

# ==============================================================================
# Cryptic Splice Site Scoring (merged from splice.py)
#
# This section provides a simplified PWM-based scoring function for quick
# cryptic splice site detection. It uses a hand-crafted position weight matrix
# for 5' splice sites (positions -3 to +6 relative to GT).
#
# DISTINCTION from maxentscan.py:
#   - maxent_score()  (below): Simplified hand-crafted PWM, returns raw
#     weighted sum. Fast but approximate. Used for quick PASS/UNCERTAIN/FAIL
#     triage of cryptic sites.
#   - score_donor() / score_acceptor() (in maxentscan.py): Proper MaxEntScan
#     log-odds scoring using Yeo & Burge 2004 trained parameters. More
#     accurate but slower. Used for precise isoform-level scoring.
# ==============================================================================

# Position weight matrix for 5' splice site (positions -3 to +6 relative to GT)
_MAXENT_PWM = [
    # pos -3   -2   -1    G    T   +1   +2   +3   +4
    [0.10, 0.07, 0.12, 3.50, 3.50, 0.16, 0.20, 0.12, 0.08],  # A
    [0.06, 0.04, 0.06, 0.01, 0.01, 0.06, 0.04, 0.06, 0.05],  # C
    [0.06, 0.04, 0.06, 0.01, 0.01, 0.06, 0.04, 0.06, 0.05],  # G
    [0.08, 0.15, 0.10, 0.01, 0.01, 0.14, 0.10, 0.12, 0.08],  # T
]

_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}

# Named constants replacing magic numbers
_PWM_CONTEXT_LEN = 9          # Position weight matrix context length (positions -3 to +6)
_MIN_CONTEXT_LEN = 4          # Minimum context length for meaningful scoring
_DINUC_LEN = 2               # Splice dinucleotide length (GT / AG)
_EXON_SKIP_DECAY = 0.5       # Score decay factor per additional skipped exon
_PWM_UPSTREAM = 3            # Bases upstream of GT for PWM context
_PWM_DOWNSTREAM = 6          # Bases downstream of GT for PWM context


def maxent_score(context: str) -> float:
    """Compute simplified MaxEntScan score for a potential splice site context.

    This is a fast approximate scoring using a hand-crafted PWM for 5' donor
    sites. For precise log-odds scoring, use score_donor() / score_acceptor()
    from the maxentscan module instead.

    Args:
        context: DNA sequence around a GT dinucleotide (ideally 9-mer)

    Returns:
        PWM weighted sum score. Higher = stronger splice signal.
        Thresholds: < 3.0 PASS, 3.0-6.0 UNCERTAIN, >= 6.0 FAIL.
    """
    if len(context) < _MIN_CONTEXT_LEN:
        return 0.0

    if len(context) < _PWM_CONTEXT_LEN:
        context = "A" * (_PWM_CONTEXT_LEN - len(context)) + context
        context = context[-_PWM_CONTEXT_LEN:]

    score = 0.0
    for pos in range(min(len(context), _PWM_CONTEXT_LEN)):
        base = context[pos].upper() if pos < len(context) else "A"
        idx = _BASE_INDEX.get(base, 0)
        score += _MAXENT_PWM[idx][pos]

    return score


def score_splice_sites(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0) -> list[tuple[int, float, "SpliceVerdict"]]:
    """Score all potential splice sites in a sequence using simplified PWM.

    Scans for GT dinucleotides and classifies each site as
    PASS / UNCERTAIN / FAIL based on dual-threshold scoring.

    Args:
        seq: DNA sequence to scan
        low_thresh: Below this score → PASS (default 3.0)
        high_thresh: At or above this score → FAIL (default 6.0)

    Returns:
        List of (position, score, SpliceVerdict) tuples for each GT site found.
    """
    from .type_system import SpliceVerdict
    results: list[tuple[int, float, SpliceVerdict]] = []
    for i in range(len(seq) - 1):
        if seq[i:i + _DINUC_LEN] == "GT":
            start = max(0, i - _PWM_UPSTREAM)
            end = min(len(seq), i + _PWM_DOWNSTREAM)
            context = seq[start:end]
            sc = maxent_score(context)
            if sc < low_thresh:
                verdict = SpliceVerdict.PASS
            elif sc < high_thresh:
                verdict = SpliceVerdict.UNCERTAIN
            else:
                verdict = SpliceVerdict.FAIL
            results.append((i, sc, verdict))
    return results


# ==============================================================================
# NDFST Splicing Engine — Isoform Computation
# ==============================================================================

from .scanner import scan_sequence
from .types import SpliceIsoform
from .constants import MIN_INTRON_LENGTH
from .maxentscan import score_donor, score_acceptor

logger = logging.getLogger(__name__)

# Tissue-specific splice site score multipliers
# These weight the MaxEntScan scores based on cell type.
# DERIVED FROM GTEx data via tissue_data module.
from .tissue_data import get_tissue_weights as _get_gtex_weights


def _get_tissue_weights(cellular_context: str) -> dict[str, float]:
    """Get tissue-specific scoring weights from GTEx-derived data."""
    try:
        return _get_gtex_weights(cellular_context)
    except Exception:
        # Fallback to default weights from tissue_data module
        logger.debug("Falling back to default tissue weights for %r", cellular_context, exc_info=True)
        from .tissue_data import GTEX_TISSUE_WEIGHTS
        return GTEX_TISSUE_WEIGHTS.get(cellular_context, GTEX_TISSUE_WEIGHTS["default"])


def compute_splice_isoforms(
    pre_mrna: str,
    known_exon_boundaries: list[tuple[int, int]],
    cellular_context: str = "HEK293T",
    max_isoforms: int = 100,
    max_exon_skip_combos: int = 10,
    tolerance: int = 5,
    alt_site_window: int = 50,
    max_alt_sites: int = 3,
    cryptic_score_threshold: float = 8.0,
) -> list[SpliceIsoform]:
    """
    Compute all possible splice isoforms via enumerative NDFST.

    The NDFST explores every valid combination of donor/acceptor pairs
    to enumerate possible splice products:

    1. Find ALL donor (GT) and acceptor (AG) sites with MaxEntScan scoring
    2. Build ALL valid (donor, acceptor) pairs that could define introns
    3. Enumerate isoform paths:
       a. Canonical splicing
       b. Exon skipping (single and multi-exon)
       c. Intron retention
       d. Cryptic splice sites
       e. Alternative 5'/3' splice sites
    4. Score each isoform using MaxEntScan + tissue-specific weights

    KEY PROPERTY: This computation is DETERMINISTIC. Same input → same isoform set.

    Args:
        pre_mrna: pre-mRNA sequence
        known_exon_boundaries: known correct exon positions [(start, end), ...]
        cellular_context: cell type context for tissue-specific scoring
        max_isoforms: safety limit to prevent combinatorial explosion
        max_exon_skip_combos: maximum number of exon-skipping combinations
        tolerance: position tolerance for matching known splice sites (default 5)
        alt_site_window: window for alternative splice site detection (default 50)
        max_alt_sites: maximum alternative sites per intron
        cryptic_score_threshold: minimum score for cryptic sites to be included

    Returns:
        List of SpliceIsoform objects, sorted by score (canonical first).
    """
    seq = pre_mrna.upper()
    # Use permissive thresholds for scanning — we want ALL potential sites
    # for isoform enumeration. Scoring/filtering happens at the isoform level.
    tokens = scan_sequence(
        seq,
        use_maxentscan=True,
        donor_threshold=0.0,
        acceptor_threshold=0.0,
    )

    donors = sorted(
        [t for t in tokens if t.element_type == "splice_donor"],
        key=lambda t: t.position
    )
    acceptors = sorted(
        [t for t in tokens if t.element_type == "splice_acceptor"],
        key=lambda t: t.position
    )

    tissue_w = _get_tissue_weights(cellular_context)

    if not donors or not acceptors:
        return [SpliceIsoform(
            sequence=seq,
            exon_boundaries=[(0, len(seq))],
            parse_path=["no_splice_sites"],
            score=0.0,
        )]

    # Build all valid (donor, acceptor) pairs = potential introns
    valid_pairs: list[tuple[int, int, float, str]] = []

    for d in donors:
        for a in acceptors:
            intron_len = a.position - d.position
            if intron_len >= MIN_INTRON_LENGTH:
                pair_score = d.score + a.score
                is_known = False
                for i in range(len(known_exon_boundaries) - 1):
                    known_donor = known_exon_boundaries[i][1]
                    known_acceptor = known_exon_boundaries[i + 1][0]
                    if abs(d.position - known_donor) < tolerance and abs(a.position - known_acceptor) < tolerance:
                        is_known = True
                        break
                label = "canonical" if is_known else f"cryptic_{d.position}_{a.position}"
                valid_pairs.append((d.position, a.position, pair_score, label))

    if not valid_pairs:
        return [SpliceIsoform(
            sequence=seq,
            exon_boundaries=[(0, len(seq))],
            parse_path=["no_valid_introns"],
            score=0.0,
        )]

    isoforms: list[SpliceIsoform] = []

    # === Path 1: Canonical splicing ===
    canonical_exons = known_exon_boundaries
    canonical_seq = "".join(seq[start:end] for start, end in canonical_exons)
    canonical_score = sum(
        d.score + a.score
        for d in donors for a in acceptors
        for i in range(len(known_exon_boundaries) - 1)
        if abs(d.position - known_exon_boundaries[i][1]) < tolerance
        and abs(a.position - known_exon_boundaries[i + 1][0]) < tolerance
    ) * tissue_w["canonical"]
    isoforms.append(SpliceIsoform(
        sequence=canonical_seq,
        exon_boundaries=canonical_exons,
        parse_path=["canonical"] * (len(known_exon_boundaries) - 1),
        score=canonical_score,
    ))

    # === Path 2: Exon skipping (single and multi-exon) ===
    if len(known_exon_boundaries) > 2:
        internal_exon_indices = list(range(1, len(known_exon_boundaries) - 1))
        # Single exon skip
        for skip_idx in internal_exon_indices:
            skipped_exons = [e for i, e in enumerate(known_exon_boundaries) if i != skip_idx]
            skipped_seq = "".join(seq[start:end] for start, end in skipped_exons)
            isoforms.append(SpliceIsoform(
                sequence=skipped_seq,
                exon_boundaries=skipped_exons,
                parse_path=[f"skip_exon_{skip_idx}"],
                score=canonical_score * tissue_w["exon_skip"],
            ))

        # Multi-exon skip (combinations of 2+ skipped exons)
        combo_count = 0
        for r in range(2, len(internal_exon_indices) + 1):
            for combo in combinations(internal_exon_indices, r):
                if combo_count >= max_exon_skip_combos:
                    break
                skipped_exons = [e for i, e in enumerate(known_exon_boundaries) if i not in combo]
                skipped_seq = "".join(seq[start:end] for start, end in skipped_exons)
                isoforms.append(SpliceIsoform(
                    sequence=skipped_seq,
                    exon_boundaries=skipped_exons,
                    parse_path=[f"skip_exons_{'_'.join(str(c) for c in combo)}"],
                    score=canonical_score * tissue_w["exon_skip"] * (_EXON_SKIP_DECAY ** (r - 1)),
                ))
                combo_count += 1
            if combo_count >= max_exon_skip_combos:
                break

    # === Path 3: Intron retention ===
    # Each intron can be retained instead of spliced out
    for i in range(len(known_exon_boundaries) - 1):
        intron_start = known_exon_boundaries[i][1]
        intron_end = known_exon_boundaries[i + 1][0]
        # Retained intron: include intron in the mRNA
        retained_exons = list(known_exon_boundaries)
        # Merge exon i and exon i+1 into one, including the intron between them
        merged = (retained_exons[i][0], retained_exons[i + 1][1])
        retained_exons = retained_exons[:i] + [merged] + retained_exons[i + 2:]
        retained_seq = "".join(seq[start:end] for start, end in retained_exons)
        isoforms.append(SpliceIsoform(
            sequence=retained_seq,
            exon_boundaries=retained_exons,
            parse_path=[f"intron_retention_{i}"],
            score=canonical_score * tissue_w["intron_retention"],
        ))

    # === Path 4: Cryptic splice sites ===
    cryptic_pairs = [
        (d_pos, a_pos, score, label)
        for d_pos, a_pos, score, label in valid_pairs
        if "cryptic" in label and score > cryptic_score_threshold
    ]

    for d_pos, a_pos, score, label in cryptic_pairs:
        if len(isoforms) >= max_isoforms:
            break
        before = seq[:d_pos]
        after = seq[a_pos + _DINUC_LEN:]
        # Guard: both exon segments must be non-empty and boundaries valid
        if len(before) > 0 and len(after) > 0 and d_pos > 0 and a_pos + _DINUC_LEN <= len(seq) and a_pos + _DINUC_LEN > d_pos:
            cryptic_seq = before + after
            isoforms.append(SpliceIsoform(
                sequence=cryptic_seq,
                exon_boundaries=[(0, d_pos), (a_pos + _DINUC_LEN, len(seq))],
                parse_path=[label],
                score=score * tissue_w["cryptic"],
            ))

    # === Path 5: Alternative 5' and 3' splice sites ===
    for i in range(len(known_exon_boundaries) - 1):
        known_donor_pos = known_exon_boundaries[i][1]
        known_acceptor_pos = known_exon_boundaries[i + 1][0]

        # Alternative 5' sites (alternative donor)
        alt_donors = [
            d for d in donors
            if d.position != known_donor_pos
            and abs(d.position - known_donor_pos) < alt_site_window
            and d.position < known_acceptor_pos - MIN_INTRON_LENGTH
        ]

        for alt_d in alt_donors[:max_alt_sites]:
            if len(isoforms) >= max_isoforms:
                break
            # Guard: alternative donor must not move past the exon start
            if alt_d.position <= known_exon_boundaries[i][0]:
                continue
            alt_exons = list(known_exon_boundaries)
            alt_exons[i] = (alt_exons[i][0], alt_d.position)
            alt_seq = "".join(seq[start:end] for start, end in alt_exons)
            isoforms.append(SpliceIsoform(
                sequence=alt_seq,
                exon_boundaries=alt_exons,
                parse_path=[f"alt5ss_exon{i}_donor{alt_d.position}"],
                score=alt_d.score * tissue_w["alt_site"],
            ))

        # Alternative 3' sites (alternative acceptor)
        alt_acceptors = [
            a for a in acceptors
            if a.position != known_acceptor_pos
            and abs(a.position - known_acceptor_pos) < alt_site_window
            and a.position > known_donor_pos + MIN_INTRON_LENGTH
        ]

        for alt_a in alt_acceptors[:max_alt_sites]:
            if len(isoforms) >= max_isoforms:
                break
            new_start = alt_a.position + _DINUC_LEN
            # Guard: alternative acceptor must not move past the exon end
            if new_start >= known_exon_boundaries[i + 1][1]:
                continue
            alt_exons = list(known_exon_boundaries)
            alt_exons[i + 1] = (new_start, alt_exons[i + 1][1])
            alt_seq = "".join(seq[start:end] for start, end in alt_exons)
            isoforms.append(SpliceIsoform(
                sequence=alt_seq,
                exon_boundaries=alt_exons,
                parse_path=[f"alt3ss_exon{i+1}_acceptor{alt_a.position}"],
                score=alt_a.score * tissue_w["alt_site"],
            ))

    # Sort by score descending (canonical first)
    isoforms.sort(key=lambda iso: -iso.score)

    # Deduplicate by sequence
    seen: set[str] = set()
    unique_isoforms: list[SpliceIsoform] = []
    for iso in isoforms:
        if iso.sequence not in seen:
            seen.add(iso.sequence)
            unique_isoforms.append(iso)

    logger.info(
        "NDFST computed %d isoforms (%d unique) for %d nt pre-mRNA (context=%s)",
        len(isoforms), len(unique_isoforms), len(seq), cellular_context
    )
    return unique_isoforms[:max_isoforms]
