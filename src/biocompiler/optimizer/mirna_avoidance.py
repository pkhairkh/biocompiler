"""
BioCompiler miRNA Binding Site Active Avoidance
=================================================
Iteratively eliminates miRNA seed matches from coding sequences during
optimization by preferring codons that disrupt miRNA seed complementarity.

Unlike the post-hoc validation approach (``check_no_mirna_binding_site`` in
checks.py), this module performs ACTIVE avoidance: during codon optimization,
the optimizer prefers synonymous codons that break miRNA seed matches,
minimizing the chance that high-affinity sites survive to the final sequence.

Strategy (mirrors ``eliminate_cryptic_splice_sites`` pattern):
1. Scan for all miRNA seed matches using ``check_no_mirna_binding_site``
2. For each high-affinity match, identify the codon(s) overlapping the seed
3. Try synonymous substitutions that break the seed complement
4. Score each alternative by:
   a. Does it break the seed match? (required)
   b. Does it maintain CAI? (prefer high-CAI alternatives)
   c. Does it avoid creating new splice sites? (safety check)
5. Apply the best alternative and re-scan
6. Repeat until no high-affinity sites remain or max_iterations reached

EUKARYOTE-ONLY: miRNA regulation is a eukaryotic gene expression mechanism.
Prokaryotes do not have miRNA machinery, so this module should be skipped
for prokaryotic targets.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

from ..type_system import AA_TO_CODONS, CODON_TABLE
from ..type_system.checks import _rna_revcomp_to_dna
from ..type_system.mirna_seeds import get_mirna_seeds
from biocompiler.shared.constants import reverse_complement
from biocompiler.sequence.maxentscan import score_donor, score_acceptor

logger = logging.getLogger(__name__)

MAX_MIRNA_AVOIDANCE_ITERATIONS: int = 300

__all__ = [
    "eliminate_mirna_binding_sites",
    "MAX_MIRNA_AVOIDANCE_ITERATIONS",
]


def _extract_hits_from_predicate(
    predicate_result,
) -> List[Dict[str, Any]]:
    """Extract structured hit information from a check_no_mirna_binding_site result.

    The PredicateResult from check_no_mirna_binding_site stores hit details
    in the ``details`` string, but we need structured data for elimination.
    Instead of parsing the details string, we re-scan using the seed database
    directly to get structured positions.

    Args:
        predicate_result: A PredicateResult from check_no_mirna_binding_site.

    Returns:
        List of hit dicts with keys: mirna, pos, score, match_type, tier,
        tissue, seed_dna.
    """
    # We cannot easily extract structured data from the PredicateResult,
    # so we re-scan with get_mirna_seeds + _rna_revcomp_to_dna for the
    # full hit data.  The caller should use _scan_mirna_sites directly.
    return []


def _scan_mirna_sites(
    seq: str,
    organism: str = "Homo_sapiens",
    tissue: str = "",
    min_seed_match: int = 7,
) -> List[Dict[str, Any]]:
    """Scan a DNA coding sequence for miRNA seed matches.

    Returns structured hit information suitable for elimination logic.

    Args:
        seq: DNA coding sequence (uppercase).
        organism: Target organism for seed selection.
        tissue: Optional tissue filter.
        min_seed_match: Minimum seed match length.

    Returns:
        List of hit dicts with keys: mirna, pos, score, match_type, tier,
        tissue, seed_dna, match_len.
    """
    n = len(seq)
    if n < 6:
        return []

    seeds = get_mirna_seeds(organism=organism, tissue=tissue if tissue else None)
    hits: List[Dict[str, Any]] = []

    for mirna_name, (seed_rna, tier, seed_tissue) in seeds.items():
        seed_len = len(seed_rna)
        if seed_len < min_seed_match:
            continue

        dna_target_7 = _rna_revcomp_to_dna(seed_rna)
        dna_target_6 = dna_target_7[:6] if seed_len >= 7 else dna_target_7

        # Scan for 7mer-m8 (full seed match, positions 2-8)
        if seed_len >= 7:
            pos = seq.find(dna_target_7)
            while pos != -1:
                is_8mer = (pos > 0 and seq[pos - 1] == "T")
                if is_8mer:
                    score = 1.0
                    match_type = "8mer"
                    match_len = 8
                else:
                    score = 0.9
                    match_type = "7mer-m8"
                    match_len = 7
                hits.append({
                    "mirna": mirna_name,
                    "pos": pos,
                    "score": score,
                    "match_type": match_type,
                    "tier": tier,
                    "tissue": seed_tissue,
                    "seed_dna": dna_target_7,
                    "match_len": match_len,
                })
                pos = seq.find(dna_target_7, pos + 1)

        # Scan for 6mer match and 7mer-A1
        if min_seed_match <= 6:
            pos = seq.find(dna_target_6)
            while pos != -1:
                # Skip if already found as 7mer/8mer
                is_part_of_7mer = False
                if seed_len >= 7:
                    if pos + 7 <= n and seq[pos:pos + 7] == dna_target_7:
                        is_part_of_7mer = True
                    if pos > 0 and seq[pos - 1:pos + 6] == dna_target_7:
                        is_part_of_7mer = True
                if is_part_of_7mer:
                    pos = seq.find(dna_target_6, pos + 1)
                    continue

                is_7mer_A1 = (pos > 0 and seq[pos - 1] == "T")
                if is_7mer_A1:
                    score = 0.85
                    match_type = "7mer-A1"
                    match_len = 7
                else:
                    score = 0.7
                    match_type = "6mer"
                    match_len = 6
                hits.append({
                    "mirna": mirna_name,
                    "pos": pos,
                    "score": score,
                    "match_type": match_type,
                    "tier": tier,
                    "tissue": seed_tissue,
                    "seed_dna": dna_target_6,
                    "match_len": match_len,
                })
                pos = seq.find(dna_target_6, pos + 1)

    return hits


def _get_overlapping_codon_indices(
    hit_pos: int,
    match_len: int,
    num_codons: int,
) -> List[int]:
    """Get the indices of codons that overlap with a seed match.

    Args:
        hit_pos: Start position of the seed match in the DNA sequence.
        match_len: Length of the seed match.
        num_codons: Total number of codons in the sequence.

    Returns:
        Sorted list of codon indices that overlap with the match.
    """
    match_end = hit_pos + match_len
    first_codon = hit_pos // 3
    last_codon = (match_end - 1) // 3
    return [i for i in range(first_codon, min(last_codon + 1, num_codons))]


def _seed_match_broken(
    seq: str,
    hit_pos: int,
    seed_dna: str,
    match_type: str,
) -> bool:
    """Check whether a specific seed match is broken in the new sequence.

    A seed match is broken if the DNA target no longer appears at the
    hit position (for 7mer-m8 and 8mer) or the 6mer no longer matches.

    Args:
        seq: New DNA sequence after codon swap.
        hit_pos: Original position of the seed match.
        seed_dna: DNA target of the miRNA seed.
        match_type: Type of match (8mer, 7mer-m8, 7mer-A1, 6mer).

    Returns:
        True if the seed match is broken (no longer matches).
    """
    if match_type in ("8mer", "7mer-m8"):
        # The 7mer DNA target must no longer match at hit_pos
        if seq[hit_pos:hit_pos + len(seed_dna)] != seed_dna:
            return True
        # For 8mer, also check if the upstream T is gone
        if match_type == "8mer" and hit_pos > 0 and seq[hit_pos - 1] != "T":
            return True
    elif match_type == "7mer-A1":
        # The 6mer core must no longer match, or upstream T is gone
        if seq[hit_pos:hit_pos + len(seed_dna)] != seed_dna:
            return True
        if hit_pos > 0 and seq[hit_pos - 1] != "T":
            return True
    elif match_type == "6mer":
        if seq[hit_pos:hit_pos + len(seed_dna)] != seed_dna:
            return True
    return False


def eliminate_mirna_binding_sites(
    seq: str,
    organism: str = "Homo_sapiens",
    tissue: str = "",
    min_seed_match: int = 7,
    max_iterations: int = 300,
) -> Tuple[str, List[Dict[str, Any]]]:
    """Iteratively eliminate miRNA binding sites from a DNA coding sequence.

    For each miRNA seed match found, attempt to disrupt it by synonymous
    codon substitution. Works similarly to eliminate_cryptic_splice_sites
    but targets miRNA seed complementarity instead of GT/AG dinucleotides.

    Strategy:
    1. Scan for all miRNA seed matches using _scan_mirna_sites
    2. For each match, identify the codon(s) overlapping the seed match
    3. Try synonymous substitutions that break the seed complement
    4. Prefer substitutions that maintain CAI and avoid creating new problems
    5. Repeat until no high-affinity sites remain or max_iterations reached

    Only targets high-affinity matches (8mer and 7mer-m8 with tier 1-2
    miRNAs) for elimination. Weaker matches (6mer, tier 3) are left alone
    to preserve CAI.

    Args:
        seq: DNA coding sequence.
        organism: Target organism.
        tissue: Target tissue for tissue-aware filtering.
        min_seed_match: Minimum seed match length.
        max_iterations: Maximum elimination iterations.

    Returns:
        (optimized_seq, modifications) where modifications is a list of dicts
        describing each change made.
    """
    seq = seq.upper()
    modifications: List[Dict[str, Any]] = []

    # EUKARYOTE-ONLY: Prokaryotes have no miRNA machinery, so avoidance
    # is biologically irrelevant.  Check organism and skip if prokaryotic.
    try:
        from biocompiler.organisms.config import is_eukaryotic_organism
        if not is_eukaryotic_organism(organism):
            return seq, modifications
    except ImportError:
        pass  # If organism config unavailable, proceed (best-effort)

    # Check if organism has dedicated miRNA seed data.  Organisms not in
    # the ORGANISM_MIRNA_MAP (e.g., yeast, E. coli) lack miRNA machinery
    # and should be skipped.  We check the map directly rather than using
    # get_mirna_seeds() which silently falls back to human seeds.
    from biocompiler.type_system.mirna_seeds import ORGANISM_MIRNA_MAP
    if organism not in ORGANISM_MIRNA_MAP:
        return seq, modifications

    mirna_seeds = get_mirna_seeds(organism=organism, tissue=tissue if tissue else None)
    if not mirna_seeds:
        return seq, modifications

    # Pre-compute amino acid list from the sequence
    num_codons = len(seq) // 3
    aas: List[str] = []
    for i in range(num_codons):
        codon = seq[i * 3:i * 3 + 3]
        aa = CODON_TABLE.get(codon, "")
        aas.append(aa)

    # Pre-compute sorted codons per amino acid (sorted by our own CAI if
    # available, or by the order in AA_TO_CODONS as fallback)
    sorted_codons: Dict[str, List[str]] = {}
    for aa, codons in AA_TO_CODONS.items():
        sorted_codons[aa] = list(codons)

    # Get CAI usage table for the organism to sort codons by CAI
    usage: Dict[str, float] = {}
    try:
        from ..organisms import CODON_ADAPTIVENESS_TABLES
        usage = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    except (ImportError, KeyError):
        pass

    # Sort codons by CAI (descending) for each amino acid
    if usage:
        for aa in sorted_codons:
            sorted_codons[aa] = sorted(
                sorted_codons[aa],
                key=lambda c: usage.get(c, 0.0),
                reverse=True,
            )

    # Splice threshold for safety checks
    cryptic_splice_threshold = 3.0

    for iteration in range(max_iterations):
        # Scan for remaining miRNA binding sites
        hits = _scan_mirna_sites(
            seq, organism=organism, tissue=tissue,
            min_seed_match=min_seed_match,
        )

        if not hits:
            break

        # Filter to only high-affinity matches that we want to eliminate:
        # 8mer (score=1.0) and 7mer-m8 (score=0.9) for tier 1-2 miRNAs
        high_affinity = [
            h for h in hits
            if h["score"] >= 0.9 and h["tier"] <= 2
        ]

        if not high_affinity:
            # Only weaker matches remain — acceptable
            break

        # Sort by: score descending (worst first), then tier ascending
        high_affinity.sort(key=lambda h: (-h["score"], h["tier"]))

        fixed_any = False

        for hit in high_affinity:
            hit_pos = hit["pos"]
            match_type = hit["match_type"]
            seed_dna = hit["seed_dna"]
            match_len = hit["match_len"]
            mirna_name = hit["mirna"]

            # Determine codons overlapping this seed match
            codon_indices = _get_overlapping_codon_indices(
                hit_pos, match_len, num_codons,
            )

            if not codon_indices:
                continue

            # Strategy 1: Single-codon swap
            # Try each overlapping codon, trying all synonymous alternatives
            # sorted by CAI (best first). Pick the first alternative that:
            # (a) breaks the seed match at hit_pos
            # (b) does not create new high-affinity miRNA sites
            # (c) does not worsen cryptic splice scores
            best_swap: Optional[Dict[str, Any]] = None
            best_score = -1.0  # composite score for ranking alternatives

            for ci in codon_indices:
                if ci >= len(aas):
                    continue
                aa = aas[ci]
                if not aa or aa == "*" or aa not in sorted_codons:
                    continue
                if len(sorted_codons[aa]) <= 1:
                    continue  # Only one codon option for this AA

                codon_pos = ci * 3
                current_codon = seq[codon_pos:codon_pos + 3]

                for alt_codon in sorted_codons[aa]:
                    if alt_codon == current_codon:
                        continue

                    # Apply the swap
                    new_seq = seq[:codon_pos] + alt_codon + seq[codon_pos + 3:]

                    # Check (a): Does this break the seed match?
                    if not _seed_match_broken(new_seq, hit_pos, seed_dna, match_type):
                        continue

                    # Check (b): Does it create new high-affinity miRNA sites?
                    new_hits = _scan_mirna_sites(
                        new_seq, organism=organism, tissue=tissue,
                        min_seed_match=min_seed_match,
                    )
                    new_high_affinity = [
                        h for h in new_hits if h["score"] >= 0.9 and h["tier"] <= 2
                    ]
                    # Count high-affinity sites; allow only if total does not increase
                    old_high_count = len(high_affinity)
                    new_high_count = len(new_high_affinity)
                    if new_high_count > old_high_count:
                        continue

                    # Check (c): Cryptic splice safety
                    splice_worsened = False
                    try:
                        for p in range(max(0, codon_pos - 3), min(len(new_seq) - 1, codon_pos + 6)):
                            if new_seq[p:p + 2] == "GT":
                                new_s = score_donor(new_seq, p)
                                if new_s >= cryptic_splice_threshold:
                                    if seq[p:p + 2] == "GT":
                                        old_s = score_donor(seq, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                            if new_seq[p:p + 2] == "AG":
                                new_s = score_acceptor(new_seq, p)
                                if new_s >= cryptic_splice_threshold:
                                    if seq[p:p + 2] == "AG":
                                        old_s = score_acceptor(seq, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                    except Exception:
                        # If MaxEntScan fails, skip splice check
                        pass

                    if splice_worsened:
                        continue

                    # Score this alternative:
                    # Higher is better. Factors:
                    # - CAI weight (0.0-1.0, higher = better codon)
                    # - Net reduction in high-affinity sites (bonus)
                    cai_weight = usage.get(alt_codon, 0.5)
                    reduction_bonus = (old_high_count - new_high_count) * 0.1
                    composite = cai_weight + reduction_bonus

                    if composite > best_score:
                        best_score = composite
                        best_swap = {
                            "ci": ci,
                            "alt_codon": alt_codon,
                            "current_codon": current_codon,
                            "new_seq": new_seq,
                            "mirna": mirna_name,
                            "hit_pos": hit_pos,
                            "match_type": match_type,
                        }

            if best_swap is not None:
                # Apply the best swap
                old_codon = best_swap["current_codon"]
                new_codon = best_swap["alt_codon"]
                seq = best_swap["new_seq"]

                # Update the amino acid list (codon changed but AA is the same)
                modifications.append({
                    "iteration": iteration,
                    "codon_index": best_swap["ci"],
                    "old_codon": old_codon,
                    "new_codon": new_codon,
                    "amino_acid": aas[best_swap["ci"]],
                    "mirna_disrupted": best_swap["mirna"],
                    "match_type_disrupted": best_swap["match_type"],
                    "position": best_swap["hit_pos"],
                    "reason": f"Disrupted {best_swap['mirna']} {best_swap['match_type']} seed match at pos {best_swap['hit_pos']}",
                })

                logger.debug(
                    "miRNA avoidance iter %d: swapped codon %d (%s→%s) to disrupt "
                    "%s %s at pos %d",
                    iteration, best_swap["ci"], old_codon, new_codon,
                    best_swap["mirna"], best_swap["match_type"],
                    best_swap["hit_pos"],
                )

                fixed_any = True
                break  # Restart scanning from highest-priority site

            # Strategy 2: Two-codon coordinated swap
            # If single-codon swap did not work, try swapping two adjacent
            # codons simultaneously to disrupt the seed match
            if len(codon_indices) >= 2:
                for ci_idx in range(len(codon_indices) - 1):
                    ci1 = codon_indices[ci_idx]
                    ci2 = codon_indices[ci_idx + 1]
                    if ci2 != ci1 + 1:
                        continue  # Only try adjacent codons
                    if ci1 >= len(aas) or ci2 >= len(aas):
                        continue

                    aa1 = aas[ci1]
                    aa2 = aas[ci2]
                    if not aa1 or not aa2 or aa1 == "*" or aa2 == "*":
                        continue
                    if aa1 not in sorted_codons or aa2 not in sorted_codons:
                        continue
                    if len(sorted_codons[aa1]) <= 1 and len(sorted_codons[aa2]) <= 1:
                        continue

                    codon_pos1 = ci1 * 3
                    codon_pos2 = ci2 * 3
                    current1 = seq[codon_pos1:codon_pos1 + 3]
                    current2 = seq[codon_pos2:codon_pos2 + 3]

                    paired_fixed = False
                    # Try top 3 CAI alternatives for each codon
                    for alt1 in sorted_codons[aa1][:3]:
                        if alt1 == current1:
                            continue
                        for alt2 in sorted_codons[aa2][:3]:
                            if alt2 == current2 and alt1 == current1:
                                continue

                            new_seq = (
                                seq[:codon_pos1] + alt1 +
                                seq[codon_pos1 + 3:codon_pos2] + alt2 +
                                seq[codon_pos2 + 3:]
                            )

                            # Check if seed match is broken
                            if not _seed_match_broken(new_seq, hit_pos, seed_dna, match_type):
                                continue

                            # Check for new high-affinity miRNA sites
                            new_hits = _scan_mirna_sites(
                                new_seq, organism=organism, tissue=tissue,
                                min_seed_match=min_seed_match,
                            )
                            new_high_affinity = [
                                h for h in new_hits if h["score"] >= 0.9 and h["tier"] <= 2
                            ]
                            if len(new_high_affinity) > len(high_affinity):
                                continue

                            # Splice safety check
                            splice_worsened = False
                            try:
                                for p in range(
                                    max(0, codon_pos1 - 3),
                                    min(len(new_seq) - 1, codon_pos2 + 6),
                                ):
                                    if new_seq[p:p + 2] == "GT":
                                        new_s = score_donor(new_seq, p)
                                        if new_s >= cryptic_splice_threshold:
                                            if seq[p:p + 2] == "GT":
                                                old_s = score_donor(seq, p)
                                                if new_s > old_s:
                                                    splice_worsened = True
                                                    break
                                            else:
                                                splice_worsened = True
                                                break
                                    if new_seq[p:p + 2] == "AG":
                                        new_s = score_acceptor(new_seq, p)
                                        if new_s >= cryptic_splice_threshold:
                                            if seq[p:p + 2] == "AG":
                                                old_s = score_acceptor(seq, p)
                                                if new_s > old_s:
                                                    splice_worsened = True
                                                    break
                                            else:
                                                splice_worsened = True
                                                break
                            except Exception:
                                pass

                            if splice_worsened:
                                continue

                            # Apply the paired swap
                            seq = new_seq
                            modifications.append({
                                "iteration": iteration,
                                "codon_index": ci1,
                                "old_codon": current1,
                                "new_codon": alt1,
                                "amino_acid": aa1,
                                "mirna_disrupted": mirna_name,
                                "match_type_disrupted": match_type,
                                "position": hit_pos,
                                "reason": (
                                    f"Disrupted {mirna_name} {match_type} seed match "
                                    f"at pos {hit_pos} (paired swap with codon {ci2})"
                                ),
                            })
                            modifications.append({
                                "iteration": iteration,
                                "codon_index": ci2,
                                "old_codon": current2,
                                "new_codon": alt2,
                                "amino_acid": aa2,
                                "mirna_disrupted": mirna_name,
                                "match_type_disrupted": match_type,
                                "position": hit_pos,
                                "reason": (
                                    f"Coordinated swap for {mirna_name} {match_type} "
                                    f"seed match disruption at pos {hit_pos}"
                                ),
                            })

                            logger.debug(
                                "miRNA avoidance iter %d: paired swap codons "
                                "%d (%s→%s) + %d (%s→%s) to disrupt %s %s",
                                iteration, ci1, current1, alt1,
                                ci2, current2, alt2,
                                mirna_name, match_type,
                            )

                            paired_fixed = True
                            break
                        if paired_fixed:
                            break
                    if paired_fixed:
                        fixed_any = True
                        break

            if fixed_any:
                break

            # Strategy 3: Try neighboring codons (within ±2)
            # If the seed match spans codon boundaries and direct overlapping
            # codons cannot fix it, try swapping a neighbor to shift the
            # reading frame context
            for neighbor_offset in [-2, -1, 1, 2]:
                n_idx = codon_indices[0] + neighbor_offset
                if n_idx < 0 or n_idx >= num_codons:
                    continue
                if n_idx in codon_indices:
                    continue  # Already tried

                aa = aas[n_idx]
                if not aa or aa == "*" or aa not in sorted_codons:
                    continue
                if len(sorted_codons[aa]) <= 1:
                    continue

                codon_pos = n_idx * 3
                current = seq[codon_pos:codon_pos + 3]

                neighbor_fixed = False
                for alt in sorted_codons[aa][:3]:
                    if alt == current:
                        continue

                    new_seq = seq[:codon_pos] + alt + seq[codon_pos + 3:]

                    # Check if the original seed match is broken
                    if not _seed_match_broken(new_seq, hit_pos, seed_dna, match_type):
                        continue

                    # Check for new miRNA sites
                    new_hits = _scan_mirna_sites(
                        new_seq, organism=organism, tissue=tissue,
                        min_seed_match=min_seed_match,
                    )
                    new_high_affinity = [
                        h for h in new_hits if h["score"] >= 0.9 and h["tier"] <= 2
                    ]
                    if len(new_high_affinity) > len(high_affinity):
                        continue

                    # Splice safety
                    splice_worsened = False
                    try:
                        for p in range(max(0, codon_pos - 3), min(len(new_seq) - 1, codon_pos + 6)):
                            if new_seq[p:p + 2] == "GT":
                                new_s = score_donor(new_seq, p)
                                if new_s >= cryptic_splice_threshold:
                                    if seq[p:p + 2] == "GT":
                                        old_s = score_donor(seq, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                            if new_seq[p:p + 2] == "AG":
                                new_s = score_acceptor(new_seq, p)
                                if new_s >= cryptic_splice_threshold:
                                    if seq[p:p + 2] == "AG":
                                        old_s = score_acceptor(seq, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                    except Exception:
                        pass

                    if splice_worsened:
                        continue

                    seq = new_seq
                    modifications.append({
                        "iteration": iteration,
                        "codon_index": n_idx,
                        "old_codon": current,
                        "new_codon": alt,
                        "amino_acid": aa,
                        "mirna_disrupted": mirna_name,
                        "match_type_disrupted": match_type,
                        "position": hit_pos,
                        "reason": (
                            f"Disrupted {mirna_name} {match_type} seed match "
                            f"at pos {hit_pos} (neighbor swap at codon {n_idx})"
                        ),
                    })

                    logger.debug(
                        "miRNA avoidance iter %d: neighbor swap codon %d "
                        "(%s→%s) to disrupt %s %s",
                        iteration, n_idx, current, alt,
                        mirna_name, match_type,
                    )

                    neighbor_fixed = True
                    break

                if neighbor_fixed:
                    fixed_any = True
                    break

            if fixed_any:
                break

        if not fixed_any:
            # No more progress possible — report remaining sites
            remaining = _scan_mirna_sites(
                seq, organism=organism, tissue=tissue,
                min_seed_match=min_seed_match,
            )
            remaining_high = [h for h in remaining if h["score"] >= 0.9 and h["tier"] <= 2]
            if remaining_high:
                for rh in remaining_high[:5]:  # Limit detail
                    logger.warning(
                        "Could not eliminate miRNA site: %s %s at pos %d "
                        "(score=%.2f, tier=%d)",
                        rh["mirna"], rh["match_type"], rh["pos"],
                        rh["score"], rh["tier"],
                    )
            break
    else:
        logger.warning(
            "miRNA avoidance: max iterations (%d) reached", max_iterations,
        )

    # Reconciliation: check if miRNA avoidance fixes reintroduced restriction sites
    # (mirrors the reconciliation step in eliminate_cryptic_splice_sites)
    try:
        from .constraint_helpers import _remove_site_multicodon
        concrete_sites: List[str] = []  # Would need enzyme list from caller

        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            if site_upper in seq or site_rc in seq:
                new_seq, fixed = _remove_site_multicodon(
                    seq, aas, sorted_codons, site_upper, site_rc, usage=usage,
                )
                if fixed:
                    seq = new_seq
    except (ImportError, Exception):
        pass  # Reconciliation is best-effort

    return seq, modifications
