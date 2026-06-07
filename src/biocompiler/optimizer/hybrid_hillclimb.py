"""
CAI hill climbing methods for HybridOptimizer.

Extracted from hybrid_optimizer.py to decompose the monolith.
Contains _cai_hill_climb() and its validation helpers.
"""

from __future__ import annotations

import math
import logging
from typing import Any

from ..type_system import CODON_TABLE, AA_TO_CODONS
from ..incremental import IncrementalSequenceState
from .hybrid_types import Violation, SEVERITY_WEIGHTS

logger = logging.getLogger(__name__)


def _cai_hill_climb(
    optimizer, seq: str, protein: str, avoid_gt: bool
) -> tuple[str, float, int]:
    """Batch CAI hill climbing: upgrade codons while maintaining constraints.

    For each position where the current codon isn't optimal (w < 1.0):
    1. Find ALL upgradeable positions and their best upgrade codons
    2. Apply ALL non-conflicting upgrades at once (build new sequence)
    3. Validate the combined result against all constraints
    4. If valid, accept all upgrades in a single pass
    5. If not valid, fall back to one-at-a-time upgrades
    6. Only iterate when an upgrade creates a new violation

    Key improvement over one-at-a-time approach:
    - Previously: try position → apply → try next position → apply → ...
    - Now: find ALL upgradeable positions → apply all at once → validate
    This reduces the number of iterations from O(n) to O(1) in the
    common case where upgrades don't conflict.
    """
    improvements = 0
    aas = list(protein)
    n_codons = len(aas)

    for _iteration in range(optimizer.max_hill_climb_iterations):
        # ── Step 1: Find ALL upgradeable positions ──
        # For each position, determine the best upgrade codon
        upgrade_plan: dict[int, str] = {}  # ci -> best upgrade codon

        for ci in range(n_codons):
            aa = aas[ci]
            if aa == "*":
                continue
            current = seq[ci * 3:ci * 3 + 3]
            current_cai = optimizer.species_cai.get(current, 0.0)
            optimal = optimizer.optimal_codon.get(aa, current)
            optimal_cai = optimizer.species_cai.get(optimal, 0.0)

            if optimal_cai > current_cai and current != optimal:
                upgrade_plan[ci] = optimal

        if not upgrade_plan:
            break

        # ── Step 2: Apply ALL upgrades at once ──
        # Build the new sequence with all planned upgrades
        seq_list = list(seq)
        for ci, new_codon in upgrade_plan.items():
            start = ci * 3
            seq_list[start:start + 3] = list(new_codon)
        batch_seq = "".join(seq_list)

        # ── Step 3: Validate the combined result ──
        if _is_valid_batch_upgrade(optimizer, seq, batch_seq, upgrade_plan, avoid_gt):
            # All upgrades are valid — accept them all at once
            seq = batch_seq
            improvements += len(upgrade_plan)
            continue

        # ── Step 4: Batch failed — try non-conflicting groups ──
        # Group upgrades by proximity. Upgrades that are far apart
        # (>= 3 codons apart) have non-overlapping local windows and
        # can be applied simultaneously.
        any_improved = False

        # Sort by CAI improvement potential (biggest first)
        sorted_cis = sorted(
            upgrade_plan.keys(),
            key=lambda ci: (
                optimizer.species_cai.get(upgrade_plan[ci], 0.0)
                - optimizer.species_cai.get(seq[ci * 3:ci * 3 + 3], 0.0)
            ),
            reverse=True,
        )

        applied_positions: set[int] = set()

        for ci in sorted_cis:
            # Skip if a nearby position was already modified in this round
            # (positions within 2 codons share local windows)
            if any(abs(ci - applied) <= 2 for applied in applied_positions):
                continue

            new_codon = upgrade_plan[ci]
            current = seq[ci * 3:ci * 3 + 3]
            if current == new_codon:
                continue

            # Try the planned upgrade
            test_seq = seq[:ci * 3] + new_codon + seq[ci * 3 + 3:]
            if _is_valid_upgrade(optimizer, seq, test_seq, ci, avoid_gt):
                seq = test_seq
                improvements += 1
                any_improved = True
                applied_positions.add(ci)
                continue

            # Try other alternatives in CAI order
            aa = aas[ci]
            for alt in optimizer.sorted_codons.get(aa, []):
                if alt == current or alt == new_codon:
                    continue
                alt_cai = optimizer.species_cai.get(alt, 0.0)
                cur_cai = optimizer.species_cai.get(current, 0.0)
                if alt_cai <= cur_cai:
                    continue  # Only try upgrades

                test_seq = seq[:ci * 3] + alt + seq[ci * 3 + 3:]
                if _is_valid_upgrade(optimizer, seq, test_seq, ci, avoid_gt):
                    seq = test_seq
                    improvements += 1
                    any_improved = True
                    applied_positions.add(ci)
                    break

        if not any_improved:
            break

    cai = optimizer._compute_cai(seq)
    return seq, cai, improvements


def _is_valid_batch_upgrade(
    optimizer,
    old_seq: str,
    new_seq: str,
    upgrade_plan: dict[int, str],
    avoid_gt: bool,
) -> bool:
    """Check if a batch of codon upgrades is valid (no new violations).

    Instead of checking each upgrade individually, validates the entire
    new sequence against all constraints. This is more efficient than
    checking one at a time when many upgrades are applied simultaneously.

    Args:
        old_seq: The sequence before upgrades.
        new_seq: The sequence after all upgrades are applied.
        upgrade_plan: Mapping of codon indices to their new codons.
        avoid_gt: Whether to avoid GT dinucleotides.

    Returns:
        True if all upgrades are valid (no new violations introduced).
    """
    # 1. No new restriction sites
    if optimizer._rs_sites:
        for site, site_rc in optimizer._rs_sites:
            if site in new_seq or (site_rc and site_rc in new_seq):
                return False

    # 2. GC still in range
    gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
    if not (optimizer.gc_lo <= gc <= optimizer.gc_hi):
        return False

    # 3. No new ATTTA motifs
    old_attta = old_seq.count("ATTTA")
    new_attta = new_seq.count("ATTTA")
    if new_attta > old_attta:
        return False

    # 4. No long T runs (6+ consecutive T)
    i = 0
    while i < len(new_seq):
        if new_seq[i] == 'T':
            j = i
            while j < len(new_seq) and new_seq[j] == 'T':
                j += 1
            if j - i >= 6:
                return False
            i = j
        else:
            i += 1

    # 5. No new avoidable GTs that aren't CAI-critical (if avoiding GT)
    if avoid_gt:
        for pos in range(len(new_seq) - 1):
            if new_seq[pos:pos+2] == "GT":
                if not optimizer._is_unavoidable_gt(new_seq, pos):
                    return False

    # 6. No new premature stop codons at upgraded positions
    for ci, new_codon in upgrade_plan.items():
        if new_codon in ("TAA", "TAG", "TGA"):
            return False

    return True


def _is_valid_upgrade(
    optimizer,
    old_seq: str,
    new_seq: str,
    ci: int,
    avoid_gt: bool,
) -> bool:
    """Check if a codon swap is a valid upgrade (no new violations).

    Instead of checking the entire sequence, only check constraints
    in the neighborhood of the swapped codon for speed.
    """
    # 1. No new restriction sites (local window check around swapped codon)
    # Use Aho-Corasick scanner for O(window_len + matches) detection
    # instead of iterating through each enzyme separately.
    max_site_len = optimizer._max_rs_site_len
    if max_site_len > 0:
        check_start = max(0, ci * 3 - max_site_len + 1)
        check_end = min(len(new_seq), ci * 3 + 3 + max_site_len - 1)
        if optimizer._ac_scanner is not None:
            if optimizer._ac_scanner.has_any_match_in_region(
                new_seq, check_start, check_end
            ):
                return False
        else:
            region = new_seq[check_start:check_end]
            for site, site_rc in optimizer._rs_sites:
                if site in region or (site_rc and site_rc in region):
                    return False

    # 2. GC still in range
    gc = (new_seq.count("G") + new_seq.count("C")) / max(len(new_seq), 1)
    if not (optimizer.gc_lo <= gc <= optimizer.gc_hi):
        return False

    # 3. No new ATTTA motifs — only block if the swap INTRODUCES a new
    # ATTTA, not if ATTTA already existed before the swap.
    old_attta = old_seq.count("ATTTA")
    new_attta = new_seq.count("ATTTA")
    if new_attta > old_attta:
        return False

    # 4. No long T runs
    i = max(0, ci * 3 - 6)
    end = min(len(new_seq), ci * 3 + 9)
    window = new_seq[i:end]
    j = 0
    while j < len(window):
        if window[j] == 'T':
            k = j
            while k < len(window) and window[k] == 'T':
                k += 1
            if k - j >= 6:
                return False
            j = k
        else:
            j += 1

    # 5. No new avoidable GTs that aren't CAI-critical (if avoiding GT)
    if avoid_gt:
        # Check GT positions in neighborhood
        start = max(0, ci * 3 - 1)
        end_pos = min(len(new_seq) - 1, ci * 3 + 4)
        for pos in range(start, end_pos):
            if new_seq[pos:pos+2] == "GT":
                if not optimizer._is_unavoidable_gt(new_seq, pos):
                    return False

    # 6. No new premature stop codons
    new_codon = new_seq[ci * 3:ci * 3 + 3]
    if new_codon in ("TAA", "TAG", "TGA"):
        return False

    return True


