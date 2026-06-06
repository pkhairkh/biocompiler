"""
BioCompiler CpG Disruption Logic
==================================
CpG dinucleotide disruption and reconciliation for the optimization pipeline.

Extracted from optimization.py for maintainability.
"""

import logging

from .constants import reverse_complement
from .maxentscan import score_donor, score_acceptor
from .optimization_helpers import _remove_site_multicodon

logger = logging.getLogger(__name__)

MAX_CPG_DISRUPTION_ITERATIONS: int = 200


def disrupt_cpg_dinucleotides(
    sequence: str,
    aas: list[str],
    sorted_codons: dict[str, list[str]],
    usage: dict[str, float],
    cryptic_splice_threshold: float,
) -> tuple[str, list[str]]:
    """Disrupt CpG dinucleotides to avoid CpG islands.

    EUKARYOTE-ONLY: CpG islands are a eukaryotic gene regulation concern.
    Prokaryotes don't methylate CpG dinucleotides, so avoidance is unnecessary.
    The caller should skip this function for prokaryotic targets.

    Strategy:
    1. Find all CG dinucleotide positions
    2. For each CG, try single-codon swap on the C-codon and/or G-codon
    3. If single swap fails, try coordinated 2-codon swap for cross-codon CG
    4. All swaps must not worsen cryptic splice scores

    Args:
        sequence: Current DNA sequence (uppercase).
        aas: List of amino acid codes (one per codon).
        sorted_codons: AA → codons sorted by CAI (descending).
        usage: Codon → CAI adaptiveness weight dict.
        cryptic_splice_threshold: Score threshold for splice site checking.

    Returns:
        Tuple of (adjusted_sequence, warnings_list).
    """
    warnings: list[str] = []

    for _cpg_iteration in range(MAX_CPG_DISRUPTION_ITERATIONS):
        cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
        if not cpg_positions:
            break
        fixed = False
        for pos in cpg_positions:
            left_ci = pos // 3          # codon containing the C
            right_ci = (pos + 1) // 3   # codon containing the G
            is_cross_codon = (left_ci != right_ci)

            # Strategy 1: Single-codon swap — try both the C-codon and the G-codon
            for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                if ci < 0 or ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci*3:ci*3+3]
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    test = sequence[:ci*3] + alt + sequence[ci*3+3:]
                    # CRITICAL: Don't worsen any cryptic splice scores.
                    codon_start = ci * 3
                    codon_end = ci * 3 + 3
                    splice_worsened = False
                    for p in range(max(0, codon_start - 3), min(len(test) - 1, codon_end + 6)):
                        if test[p:p+2] == "GT":
                            new_s = score_donor(test, p)
                            if new_s >= cryptic_splice_threshold:
                                if sequence[p:p+2] == "GT":
                                    old_s = score_donor(sequence, p)
                                    if new_s > old_s:
                                        splice_worsened = True
                                        break
                                else:
                                    splice_worsened = True
                                    break
                        if test[p:p+2] == "AG":
                            new_s = score_acceptor(test, p)
                            if new_s >= cryptic_splice_threshold:
                                if sequence[p:p+2] == "AG":
                                    old_s = score_acceptor(sequence, p)
                                    if new_s > old_s:
                                        splice_worsened = True
                                        break
                                else:
                                    splice_worsened = True
                                    break
                    if splice_worsened:
                        continue
                    # Check that the specific CG at pos is eliminated and net CpG decreases
                    if test[pos:pos+2] != "CG":
                        new_cpg_count = sum(1 for i in range(len(test) - 1) if test[i:i+2] == "CG")
                        old_cpg_count = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
                        if new_cpg_count < old_cpg_count:
                            sequence = test
                            fixed = True
                            break
                if fixed:
                    break

            # Strategy 2: Coordinated 2-codon swap for cross-codon CG
            if not fixed and is_cross_codon and 0 <= left_ci < len(aas) and 0 <= right_ci < len(aas):
                left_aa = aas[left_ci]
                right_aa = aas[right_ci]
                left_current = sequence[left_ci*3:left_ci*3+3]
                right_current = sequence[right_ci*3:right_ci*3+3]
                for left_alt in sorted_codons[left_aa]:
                    if left_alt == left_current:
                        continue
                    for right_alt in sorted_codons[right_aa]:
                        if right_alt == right_current and left_alt == left_current:
                            continue
                        # Skip if the combination still has CG at the boundary
                        if left_alt[-1] == "C" and right_alt[0] == "G":
                            continue
                        test = list(sequence)
                        test[left_ci*3:left_ci*3+3] = list(left_alt)
                        test[right_ci*3:right_ci*3+3] = list(right_alt)
                        test_str = "".join(test)
                        # Splice check across wider region (both codons)
                        splice_worsened = False
                        check_start = max(0, left_ci * 3 - 3)
                        check_end = min(len(test_str) - 1, right_ci * 3 + 9)
                        for p in range(check_start, check_end):
                            if test_str[p:p+2] == "GT":
                                new_s = score_donor(test_str, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "GT":
                                        old_s = score_donor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                            if test_str[p:p+2] == "AG":
                                new_s = score_acceptor(test_str, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "AG":
                                        old_s = score_acceptor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                        if splice_worsened:
                            continue
                        # Verify net CpG reduction
                        new_cpg_count = sum(1 for i in range(len(test_str) - 1) if test_str[i:i+2] == "CG")
                        old_cpg_count = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
                        if new_cpg_count < old_cpg_count:
                            sequence = test_str
                            fixed = True
                            break
                    if fixed:
                        break

            if fixed:
                break
        if not fixed:
            break

    return sequence, warnings


def reconcile_cpg_sites(
    sequence: str,
    aas: list[str],
    sorted_codons: dict[str, list[str]],
    usage: dict[str, float],
    cryptic_splice_threshold: float,
    concrete_sites: list[str],
) -> tuple[str, list[str]]:
    """CpG reconciliation after restriction site reconciliation.

    This is a second pass of CpG disruption that additionally ensures
    restriction sites are not reintroduced. It runs after the main
    restriction site reconciliation step.

    EUKARYOTE-ONLY: The caller should skip this function for prokaryotic targets.

    Args:
        sequence: Current DNA sequence (uppercase).
        aas: List of amino acid codes (one per codon).
        sorted_codons: AA → codons sorted by CAI (descending).
        usage: Codon → CAI adaptiveness weight dict.
        cryptic_splice_threshold: Score threshold for splice site checking.
        concrete_sites: List of concrete restriction site sequences.

    Returns:
        Tuple of (adjusted_sequence, warnings_list).
    """
    warnings: list[str] = []

    for _cpg_iter in range(MAX_CPG_DISRUPTION_ITERATIONS):
        cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
        if not cpg_positions:
            break
        fixed = False
        for pos in cpg_positions:
            left_ci = pos // 3
            right_ci = (pos + 1) // 3
            is_cross_codon = (left_ci != right_ci)

            # Strategy 1: Single-codon swap — try both the C-codon and the G-codon
            for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                if ci < 0 or ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci*3:ci*3+3]
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    test = sequence[:ci*3] + alt + sequence[ci*3+3:]
                    # Must not reintroduce restriction sites
                    site_ok = all(
                        s not in test and reverse_complement(s) not in test
                        for s in concrete_sites
                    )
                    if not site_ok:
                        continue
                    # Must not worsen cryptic splice scores
                    codon_start = ci * 3
                    codon_end = ci * 3 + 3
                    splice_worsened = False
                    for p in range(max(0, codon_start - 3), min(len(test) - 1, codon_end + 6)):
                        if test[p:p+2] == "GT":
                            new_s = score_donor(test, p)
                            if new_s >= cryptic_splice_threshold:
                                if sequence[p:p+2] == "GT":
                                    old_s = score_donor(sequence, p)
                                    if new_s > old_s:
                                        splice_worsened = True
                                        break
                                else:
                                    splice_worsened = True
                                    break
                        if test[p:p+2] == "AG":
                            new_s = score_acceptor(test, p)
                            if new_s >= cryptic_splice_threshold:
                                if sequence[p:p+2] == "AG":
                                    old_s = score_acceptor(sequence, p)
                                    if new_s > old_s:
                                        splice_worsened = True
                                        break
                                else:
                                    splice_worsened = True
                                    break
                    if splice_worsened:
                        continue
                    # Check that the specific CG at pos is eliminated and net CpG decreases
                    if test[pos:pos+2] != "CG":
                        new_cpg_count = sum(1 for i in range(len(test) - 1) if test[i:i+2] == "CG")
                        old_cpg_count = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
                        if new_cpg_count < old_cpg_count:
                            sequence = test
                            fixed = True
                            break
                if fixed:
                    break

            # Strategy 2: Coordinated 2-codon swap for cross-codon CG
            if not fixed and is_cross_codon and 0 <= left_ci < len(aas) and 0 <= right_ci < len(aas):
                left_aa = aas[left_ci]
                right_aa = aas[right_ci]
                left_current = sequence[left_ci*3:left_ci*3+3]
                right_current = sequence[right_ci*3:right_ci*3+3]
                for left_alt in sorted_codons[left_aa]:
                    if left_alt == left_current:
                        continue
                    for right_alt in sorted_codons[right_aa]:
                        if right_alt == right_current and left_alt == left_current:
                            continue
                        if left_alt[-1] == "C" and right_alt[0] == "G":
                            continue
                        test = list(sequence)
                        test[left_ci*3:left_ci*3+3] = list(left_alt)
                        test[right_ci*3:right_ci*3+3] = list(right_alt)
                        test_str = "".join(test)
                        # Must not reintroduce restriction sites
                        site_ok = all(
                            s not in test_str and reverse_complement(s) not in test_str
                            for s in concrete_sites
                        )
                        if not site_ok:
                            continue
                        # Splice check across wider region (both codons)
                        splice_worsened = False
                        check_start = max(0, left_ci * 3 - 3)
                        check_end = min(len(test_str) - 1, right_ci * 3 + 9)
                        for p in range(check_start, check_end):
                            if test_str[p:p+2] == "GT":
                                new_s = score_donor(test_str, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "GT":
                                        old_s = score_donor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                            if test_str[p:p+2] == "AG":
                                new_s = score_acceptor(test_str, p)
                                if new_s >= cryptic_splice_threshold:
                                    if sequence[p:p+2] == "AG":
                                        old_s = score_acceptor(sequence, p)
                                        if new_s > old_s:
                                            splice_worsened = True
                                            break
                                    else:
                                        splice_worsened = True
                                        break
                        if splice_worsened:
                            continue
                        # Verify net CpG reduction
                        new_cpg_count = sum(1 for i in range(len(test_str) - 1) if test_str[i:i+2] == "CG")
                        old_cpg_count = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
                        if new_cpg_count < old_cpg_count:
                            sequence = test_str
                            fixed = True
                            break
                    if fixed:
                        break

            if fixed:
                break
        if not fixed:
            break

    return sequence, warnings
