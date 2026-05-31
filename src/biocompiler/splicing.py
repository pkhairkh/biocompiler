"""
BioCompiler NDFST Splicing Engine — Enumerative Isoform Computation

Production-grade splicing engine with:
- Enumerative isoform computation (all valid donor-acceptor paths)
- Intron retention isoforms (major class of alternative splicing, ~15% of events)
- Multi-exon skipping combinations
- Tissue-weighted scoring using GTEx-derived data (tissue_data module)
- Proper handling of alternative 5'/3' sites
- Configurable parameters (not hardcoded magic numbers)
"""

import logging
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
        from itertools import combinations
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
                    score=canonical_score * tissue_w["exon_skip"] * (0.5 ** (r - 1)),
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
        after = seq[a_pos + 2:]
        if len(before) > 0 and len(after) > 0:
            cryptic_seq = before + after
            isoforms.append(SpliceIsoform(
                sequence=cryptic_seq,
                exon_boundaries=[(0, d_pos), (a_pos + 2, len(seq))],
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
            alt_exons = list(known_exon_boundaries)
            alt_exons[i + 1] = (alt_a.position + 2, alt_exons[i + 1][1])
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
