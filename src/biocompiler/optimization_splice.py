"""
BioCompiler Splice Site Elimination Logic
============================================
Cryptic splice donor/acceptor site elimination for the optimization pipeline.

Extracted from optimization.py for maintainability.
"""

import logging

from .type_system import AA_TO_CODONS
from .constants import reverse_complement
from .maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
from .optimization_helpers import (
    _find_gt_free_codons,
    _find_ag_free_codons,
    _is_unavoidable_gt_aa,
    _gt_free_cai_ratio,
)

logger = logging.getLogger(__name__)

MAX_SPLICE_ELIMINATION_ITERATIONS: int = 300
ELIMINATED_SITE_SCORE: float = -999.0
TOP_CAI_ALTERNATIVES: int = 3


def eliminate_cryptic_splice_sites(
    sequence: str,
    aas: list[str],
    sorted_codons: dict[str, list[str]],
    usage: dict[str, float],
    cryptic_splice_threshold: float,
    concrete_sites: list[str],
) -> tuple[str, list[str]]:
    """Eliminate cryptic splice donor and acceptor sites.

    EUKARYOTE-ONLY: Prokaryotes have no spliceosome, so cryptic splice
    sites are biologically irrelevant. The caller should skip this
    function for prokaryotic targets.

    Strategy:
    1. Find all strong donor (GT) and acceptor (AG) positions
    2. Sort by score descending for priority
    3. Try GT-free / AG-free codon swaps first (guaranteed elimination)
    4. If that fails, try single-codon swap with different codon
    5. If that fails, try 2-codon coordinated swap to disrupt context

    Args:
        sequence: Current DNA sequence (uppercase).
        aas: List of amino acid codes (one per codon).
        sorted_codons: AA → codons sorted by CAI (descending).
        usage: Codon → CAI adaptiveness weight dict.
        cryptic_splice_threshold: Score threshold for cryptic splice sites.
        concrete_sites: List of concrete restriction site sequences.

    Returns:
        Tuple of (adjusted_sequence, warnings_list).
    """
    warnings: list[str] = []

    for iteration in range(MAX_SPLICE_ELIMINATION_ITERATIONS):
        max_d = max_donor_score(sequence)
        max_a = max_acceptor_score(sequence)
        if max_d < cryptic_splice_threshold and max_a < cryptic_splice_threshold:
            break

        fixed_any = False

        # Try to eliminate strong donors (sorted by score descending for priority)
        if max_d >= cryptic_splice_threshold:
            # Collect all strong donor positions with scores, sort by score descending
            donor_sites = []
            for i in range(len(sequence) - 1):
                if sequence[i:i+2] == "GT":
                    s = score_donor(sequence, i)
                    if s >= cryptic_splice_threshold:
                        donor_sites.append((i, s))
            donor_sites.sort(key=lambda x: x[1], reverse=True)

            for gt_pos, gt_score in donor_sites:
                codon_idx = gt_pos // 3
                if codon_idx >= len(aas):
                    continue
                aa = aas[codon_idx]

                # Check if GT is unavoidable for this amino acid
                if _is_unavoidable_gt_aa(aa):
                    continue

                # Check if avoiding GT would sacrifice too much CAI.
                cai_ratio = _gt_free_cai_ratio(aa, usage)
                if cai_ratio < 0.5:
                    continue

                # Strategy 1: GT-free codon swap (guaranteed to eliminate GT)
                gt_free = _find_gt_free_codons(aa)
                # Sort by CAI (best first) so we prefer high-CAI alternatives
                gt_free_sorted = sorted(
                    gt_free,
                    key=lambda c: usage.get(c, 0.0),
                    reverse=True,
                )
                if gt_free_sorted:
                    for alt in gt_free_sorted:
                        seq_list = list(sequence)
                        seq_list[codon_idx*3:codon_idx*3+3] = list(alt)
                        test = "".join(seq_list)
                        # Verify the GT at this position is eliminated
                        if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                            # GT still present at this position (shouldn't happen with GT-free codon,
                            # but GT might straddle codon boundary)
                            new_s = score_donor(test, gt_pos)
                            if new_s < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                        else:
                            # GT eliminated — verify no new strong donors created globally
                            new_max_d = max_donor_score(test)
                            if new_max_d < max_d or new_max_d < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                            # Even if new donors appear, this position is fixed;
                            # they'll be handled in subsequent iterations
                            sequence = test
                            fixed_any = True
                            break

                if fixed_any:
                    break  # Restart scanning from highest-scoring site

                # Strategy 2: For Valine (no GT-free codons) or if Strategy 1 didn't work,
                # try 2-codon coordinated swap (GT codon + each neighbor) to disrupt
                # the 9-mer splice context.
                # Also try single-codon swap for V positions where different V codons
                # have different splice scores due to context.
                current = sequence[codon_idx*3:codon_idx*3+3]
                # First try single-codon swap (different V codon may give different score)
                for v_alt in sorted_codons[aa]:
                    if v_alt == current:
                        continue
                    seq_list = list(sequence)
                    seq_list[codon_idx*3:codon_idx*3+3] = list(v_alt)
                    test = "".join(seq_list)
                    if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                        new_s = score_donor(test, gt_pos)
                    else:
                        new_s = ELIMINATED_SITE_SCORE  # GT eliminated (cross-boundary)
                    if new_s < cryptic_splice_threshold:
                        sequence = test
                        fixed_any = True
                        break
                if fixed_any:
                    break  # Restart scanning

                # Then try 2-codon coordinated swap
                for neighbor_offset in [-2, -1, 1, 2]:
                    n_idx = codon_idx + neighbor_offset
                    if 0 <= n_idx < len(aas):
                        n_aa = aas[n_idx]
                        n_current = sequence[n_idx*3:n_idx*3+3]
                        # For the GT codon, try top 3 alternatives by CAI
                        for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                            # For the neighbor, try all alternatives
                            for n_alt in sorted_codons[n_aa]:
                                if n_alt == n_current and v_alt == current:
                                    continue
                                test = list(sequence)
                                test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                test[n_idx*3:n_idx*3+3] = list(n_alt)
                                test_str = "".join(test)
                                if gt_pos < len(test_str) - 1 and test_str[gt_pos:gt_pos+2] == "GT":
                                    new_s = score_donor(test_str, gt_pos)
                                else:
                                    new_s = ELIMINATED_SITE_SCORE  # GT eliminated
                                if new_s < cryptic_splice_threshold:
                                    sequence = test_str
                                    fixed_any = True
                                    break
                            if fixed_any:
                                break
                        if fixed_any:
                            break
                    if fixed_any:
                        break

                if fixed_any:
                    break  # Restart scanning

        # Try to eliminate strong acceptors (same strategy, using AG-free codons)
        if not fixed_any and max_a >= cryptic_splice_threshold:
            # Collect all strong acceptor positions with scores, sort by score descending
            acceptor_sites = []
            for i in range(len(sequence) - 1):
                if sequence[i:i+2] == "AG":
                    s = score_acceptor(sequence, i)
                    if s >= cryptic_splice_threshold:
                        acceptor_sites.append((i, s))
            acceptor_sites.sort(key=lambda x: x[1], reverse=True)

            for ag_pos, ag_score in acceptor_sites:
                codon_idx = ag_pos // 3
                if codon_idx >= len(aas):
                    continue
                aa = aas[codon_idx]

                # Strategy 1: AG-free codon swap (guaranteed to eliminate AG)
                ag_free = _find_ag_free_codons(aa)
                ag_free_sorted = sorted(
                    ag_free,
                    key=lambda c: usage.get(c, 0.0),
                    reverse=True,
                )
                if ag_free_sorted:
                    for alt in ag_free_sorted:
                        test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                        if ag_pos < len(test) - 1 and test[ag_pos:ag_pos+2] == "AG":
                            new_s = score_acceptor(test, ag_pos)
                            if new_s < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                        else:
                            # AG eliminated
                            new_max_a = max_acceptor_score(test)
                            if new_max_a < max_a or new_max_a < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                            sequence = test
                            fixed_any = True
                            break

                if fixed_any:
                    break

                # Strategy 2: Single-codon swap then 2-codon context disruption for AG positions
                current = sequence[codon_idx*3:codon_idx*3+3]
                # First try single-codon swap (different codon may give different score)
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                    if ag_pos < len(test) - 1 and test[ag_pos:ag_pos+2] == "AG":
                        new_s = score_acceptor(test, ag_pos)
                    else:
                        new_s = ELIMINATED_SITE_SCORE  # AG eliminated
                    if new_s < cryptic_splice_threshold:
                        sequence = test
                        fixed_any = True
                        break
                if fixed_any:
                    break

                # Then try 2-codon coordinated swap
                for neighbor_offset in [-2, -1, 1, 2]:
                    n_idx = codon_idx + neighbor_offset
                    if 0 <= n_idx < len(aas):
                        n_aa = aas[n_idx]
                        n_current = sequence[n_idx*3:n_idx*3+3]
                        for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                            for n_alt in sorted_codons[n_aa]:
                                if n_alt == n_current and v_alt == current:
                                    continue
                                test = list(sequence)
                                test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                test[n_idx*3:n_idx*3+3] = list(n_alt)
                                test_str = "".join(test)
                                if ag_pos < len(test_str) - 1 and test_str[ag_pos:ag_pos+2] == "AG":
                                    new_s = score_acceptor(test_str, ag_pos)
                                else:
                                    new_s = ELIMINATED_SITE_SCORE
                                if new_s < cryptic_splice_threshold:
                                    sequence = test_str
                                    fixed_any = True
                                    break
                            if fixed_any:
                                break
                        if fixed_any:
                            break
                    if fixed_any:
                        break

                if fixed_any:
                    break

        if not fixed_any:
            # No more progress — some sites are unrepairable by codon swaps
            max_d = max_donor_score(sequence)
            max_a = max_acceptor_score(sequence)
            if max_d >= cryptic_splice_threshold:
                warnings.append(
                    f"Cryptic splice donors remain: max_donor={max_d:.2f} "
                    f"(threshold={cryptic_splice_threshold}). "
                    f"Some positions may require amino acid substitution (e.g., V->I)"
                )
            if max_a >= cryptic_splice_threshold:
                warnings.append(
                    f"Cryptic splice acceptors remain: max_acceptor={max_a:.2f} "
                    f"(threshold={cryptic_splice_threshold})"
                )
            break
    else:
        warnings.append("Cryptic splice elimination: max iterations reached")

    # Reconciliation: check if cryptic splice fixes reintroduced restriction sites
    for site_upper in concrete_sites:
        site_rc = reverse_complement(site_upper)
        if site_upper in sequence or site_rc in sequence:
            from .optimization_helpers import _remove_site_multicodon
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
            )
            if fixed:
                sequence = new_seq

    return sequence, warnings
