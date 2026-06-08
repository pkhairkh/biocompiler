"""
Post-processing passes for the BioOptimizer pipeline.

Contains GC adjustment, CpG elimination, CAI recovery, mRNA stability,
and related helper functions.  These are extracted from pipeline.py for
maintainability; the public API is unchanged (pipeline.py re-exports them
as BioOptimizer methods).

All functions that need optimizer configuration take a ``cfg: _OptConfig``
parameter instead of ``self``.
"""

from typing import Dict, List, Optional, Tuple, Any

import logging
import math

from ..type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, PredicateResult,
    check_no_cpg_island,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
)
from ..incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from ..sliding_gc import check_sliding_gc
from .cai import _count_dinucs_fast, _BatchSwapScorer
from .constraints import _count_gts, _is_unavoidable_gt
from .greedy import TOP_CAI_ALTERNATIVES
from .utils import _OptConfig, MAX_GC_ADJUSTMENT_ITERATIONS, MAX_ATTTA_MOTIF_ITERATIONS
from . import validation as _validation

logger = logging.getLogger(__name__)

def _check_no_restriction_sites(cfg: "_OptConfig", seq: str) -> bool:
    """Check that sequence doesn't contain any restriction enzyme sites."""
    from ..restriction_sites import get_recognition_site
    for enzyme in cfg.enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        if site in seq:
            return False
    return True


def _check_cpg_swap_gt_safe(old_list: list, new_list: list, codon_start: int) -> bool:
    """Check that a codon swap doesn't create new avoidable GT dinucleotides.

    We allow unavoidable GTs (e.g., Valine codons) but not new avoidable ones.
    This is less strict than requiring GT count to not increase at all.
    """
    old_seq = "".join(old_list)
    new_seq = "".join(new_list)
    codon_end = codon_start + 3

    # Check the local region around the swapped codon for new GTs
    check_start = max(0, codon_start - 1)
    check_end = min(len(new_seq) - 1, codon_end + 1)

    for i in range(check_start, check_end):
        if new_seq[i:i+2] == "GT" and old_seq[i:i+2] != "GT":
            # New GT created — check if it's unavoidable
            if not _is_unavoidable_gt(new_seq, i):
                return False  # Avoidable new GT — reject this swap

    return True


def fix_gc_after_cpg(cfg: "_OptConfig", seq: str) -> str:
    """Fix GC content after CpG reconciliation without reintroducing CG dinucleotides.

    Only adjusts GC if it's outside [0.30, 0.70].
    """
    gc_val = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
    if 0.30 <= gc_val <= 0.70:
        return seq

    seq_list = list(seq)
    gc_count = sum(1 for b in seq_list if b in "GC")
    n_bases = len(seq_list)
    target_gc = 0.50 if gc_val > 0.70 else 0.30

    for iteration in range(MAX_GC_ADJUSTMENT_ITERATIONS):
        gc_val = gc_count / n_bases
        if 0.30 <= gc_val <= 0.70:
            break

        best_alt = None
        best_ci = -1
        best_diff = abs(gc_val - target_gc)
        best_gc_delta = 0

        for ci in range(len(seq_list) // 3):
            codon_start = ci * 3
            if codon_start + 3 > len(seq_list):
                break
            codon = "".join(seq_list[codon_start:codon_start + 3])
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue
            current_gc = sum(1 for b in codon if b in "GC")
            for alt in AA_TO_CODONS.get(aa, []):
                if alt == codon:
                    continue
                alt_gc = sum(1 for b in alt if b in "GC")
                new_gc_count = gc_count - current_gc + alt_gc
                new_frac = new_gc_count / n_bases
                diff = abs(new_frac - target_gc)
                if diff < best_diff:
                    # Check this swap doesn't reintroduce CGs
                    test_list = seq_list[:]
                    for k, b in enumerate(alt):
                        test_list[codon_start + k] = b
                    test_seq = "".join(test_list)
                    # Must not increase CG count
                    new_cg = _count_dinucs_fast(test_seq, "CG")[0]
                    old_cg = _count_dinucs_fast(seq, "CG")[0]
                    if new_cg > old_cg:
                        continue
                    # Must not reintroduce restriction sites
                    if not _check_no_restriction_sites(cfg, test_seq):
                        continue
                    best_diff = diff
                    best_alt = alt
                    best_ci = ci
                    best_gc_delta = alt_gc - current_gc

        if best_alt is None:
            break
        codon_start = best_ci * 3
        for k, b in enumerate(best_alt):
            seq_list[codon_start + k] = b
        gc_count += best_gc_delta

    return "".join(seq_list)


def step_avoid_cpg_islands(cfg: "_OptConfig", seq: str) -> str:
    """Avoid CpG Islands step (optimized with incremental state).

    CpG island avoidance by synonymous substitution. Also handles
    cross-codon CG dinucleotides by two-codon coordination.

    Uses IncrementalSequenceState for O(1) GT/CG tracking and
    direct CG position lookup instead of window scanning.
    """
    state = IncrementalSequenceState(seq)
    codon_cache = CodonCache(cfg.species_cai)
    changed = True
    iterations = 0
    max_iterations = 80

    while changed and iterations < max_iterations:
        changed = False
        iterations += 1

        for start in range(0, len(state.sequence) - cfg.cpg_window + 1, 3):
            window_seq = state.sequence[start:start + cfg.cpg_window]
            c_count = window_seq.count("C")
            g_count = window_seq.count("G")
            cg_count = _count_dinucs_fast(window_seq, "CG")[0]
            expected = (c_count * g_count) / len(window_seq) if len(window_seq) > 0 else 0
            obs_exp = cg_count / expected if expected > 0 else 0.0

            if obs_exp <= cfg.cpg_threshold:
                continue

            # Find CG dinucleotides in this window and try to break them
            # Use incremental state's CG position data for efficiency
            cg_in_window = [p for p in state.cg_positions_list()
                           if start <= p < start + cfg.cpg_window - 1]

            for i in cg_in_window:
                left_ci = i // 3
                right_ci = (i + 1) // 3
                is_cross_codon = (left_ci != right_ci)

                # Strategy 1: Single-codon swap
                codon_fixed = False
                for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                    codon = state.get_codon(ci)
                    aa = state.get_aa(ci)
                    if aa is None or aa == "*":
                        continue

                    for alt in codon_cache.get_sorted_codons(aa):
                        if alt == codon:
                            continue
                        if not is_cross_codon and "CG" in alt:
                            continue
                        if cfg.avoid_gt and "GT" in alt:
                            continue
                        if codon_cache.get_cai(alt) < cfg.min_cai:
                            continue

                        # O(1) boundary check
                        left_gt, right_gt = state.boundary_creates_gt(ci, alt)
                        if left_gt or right_gt:
                            continue

                        # O(1) GT/CG change check
                        if cfg.avoid_gt and state.would_gt_increase(ci, alt):
                            continue

                        # Apply swap and check the specific CG position is gone
                        old_codon = state.swap_codon(ci, alt)
                        if state.sequence[i] != "C" or state.sequence[i+1] != "G":
                            changed = True
                            codon_fixed = True
                            break
                        else:
                            state.swap_codon(ci, old_codon)  # Rollback

                    if codon_fixed:
                        break

                if codon_fixed:
                    break

                # Strategy 2: Coordinated 2-codon swap for cross-codon CG
                if is_cross_codon and not codon_fixed:
                    left_codon = state.get_codon(left_ci)
                    right_codon = state.get_codon(right_ci)
                    left_aa = state.get_aa(left_ci)
                    right_aa = state.get_aa(right_ci)

                    if left_aa and right_aa and left_aa != "*" and right_aa != "*":
                        old_gt_count = state.gt_count
                        for left_alt in codon_cache.get_sorted_codons(left_aa):
                            if left_alt == left_codon:
                                continue
                            if cfg.avoid_gt and "GT" in left_alt:
                                continue
                            if codon_cache.get_cai(left_alt) < cfg.min_cai:
                                continue
                            for right_alt in codon_cache.get_sorted_codons(right_aa):
                                if left_alt[-1] == "C" and right_alt[0] == "G":
                                    continue
                                if cfg.avoid_gt and "GT" in right_alt:
                                    continue
                                if codon_cache.get_cai(right_alt) < cfg.min_cai:
                                    continue
                                # Check cross-codon GT effects
                                left_gt, _ = state.boundary_creates_gt(left_ci, left_alt)
                                if left_gt:
                                    continue
                                if left_alt[-1] + right_alt[0] == "GT":
                                    continue
                                _, right_gt = state.boundary_creates_gt(right_ci, right_alt)
                                if right_gt:
                                    continue

                                # Apply both swaps and check incrementally
                                old_left = state.swap_codon(left_ci, left_alt)
                                old_right = state.swap_codon(right_ci, right_alt)

                                # O(1) GT count check (state.gt_count is updated incrementally)
                                if state.gt_count <= old_gt_count:
                                    # Verify CG at position i is gone
                                    if state.sequence[i] != "C" or state.sequence[i+1] != "G":
                                        changed = True
                                        codon_fixed = True
                                        break

                                # Rollback
                                state.swap_codon(right_ci, old_right)
                                state.swap_codon(left_ci, old_left)

                            if codon_fixed:
                                break

                if changed:
                    break
            if changed:
                break

    return state.sequence

# Deprecated alias — use step_avoid_cpg_islands instead
_phase5_avoid_cpg_islands = step_avoid_cpg_islands


def step_remove_instability_motifs(cfg: "_OptConfig", seq: str) -> str:
    """Remove ATTTA mRNA instability motifs by synonymous codon substitution.

    ATTTA motifs are associated with mRNA instability in eukaryotic genes.
    This step identifies ATTTA pentamers and attempts to disrupt them by
    swapping one of the overlapping codons to a synonymous alternative that
    does not contain the motif, while preserving GT/CG/RS constraints where
    possible. Falls back to accepting minor constraint relaxation if needed
    to eliminate ATTTA (instability removal has higher priority).
    """
    from ..restriction_sites import get_recognition_site as _get_site

    for iteration in range(MAX_ATTTA_MOTIF_ITERATIONS):
        pos = seq.find("ATTTA")
        if pos == -1:
            break

        n_codons = len(seq) // 3
        best_swap = None  # (test_seq, score) where lower score = better
        # Score: 0 = perfect (ATTTA gone, no new issues), 1 = ATTTA gone + new RS, 2 = ATTTA gone + worse GT

        # Try swapping each overlapping codon
        for offset in range(5):
            ci = (pos + offset) // 3
            if ci < 0 or ci >= n_codons:
                continue
            codon = seq[ci*3:ci*3+3]
            aa = CODON_TABLE.get(codon, "")
            if not aa or aa == "*":
                continue
            alternatives = sorted(
                AA_TO_CODONS.get(aa, [codon]),
                key=lambda c: cfg.species_cai.get(c, 0.0),
                reverse=True,
            )
            for alt in alternatives:
                if alt == codon:
                    continue
                test = seq[:ci*3] + alt + seq[ci*3+3:]
                if "ATTTA" not in test:
                    # Check constraints
                    rs_ok = all(
                        not (_s := _get_site(enz)) or _s not in test
                        for enz in cfg.enzymes
                    )
                    from ..type_system import check_no_avoidable_gt as _check_gt
                    gt_result = _check_gt(test)
                    old_gt = _check_gt(seq)
                    gt_ok = gt_result.passed or not (not gt_result.passed and old_gt.passed)

                    if rs_ok and gt_ok:
                        # Perfect swap — use immediately
                        return step_remove_instability_motifs(cfg, test)
                    elif gt_ok and not rs_ok:
                        score = 1  # New RS but GT preserved (RS easier to fix later)
                    elif rs_ok and not gt_ok:
                        score = 2  # GT worsened but no new RS
                    else:
                        score = 3  # Both worsened

                    if best_swap is None or score < best_swap[1]:
                        best_swap = (test, score)

        if best_swap is not None:
            # Use the best available swap (even if not perfect)
            # After applying, run a reconciliation pass to fix any introduced issues
            seq = best_swap[0]
            # Re-check and fix restriction sites
            for enz in cfg.enzymes:
                site = _get_site(enz)
                if site and site in seq:
                    ci_rs = seq.find(site) // 3
                    if 0 <= ci_rs < n_codons:
                        aa = CODON_TABLE.get(seq[ci_rs*3:ci_rs*3+3], "")
                        current = seq[ci_rs*3:ci_rs*3+3]
                        for alt in AA_TO_CODONS.get(aa, [current]):
                            if alt != current:
                                test2 = seq[:ci_rs*3] + alt + seq[ci_rs*3+3:]
                                if site not in test2 and "ATTTA" not in test2:
                                    seq = test2
                                    break
            continue

        # No single-codon swap worked; try 2-codon coordinated swap
        fixed = False
        for offset1 in range(5):
            ci1 = (pos + offset1) // 3
            if ci1 < 0 or ci1 >= n_codons:
                continue
            aa1 = CODON_TABLE.get(seq[ci1*3:ci1*3+3], "")
            if not aa1 or aa1 == "*":
                continue
            for offset2 in range(offset1 + 1, 5):
                ci2 = (pos + offset2) // 3
                if ci2 < 0 or ci2 >= n_codons or ci2 == ci1:
                    continue
                aa2 = CODON_TABLE.get(seq[ci2*3:ci2*3+3], "")
                if not aa2 or aa2 == "*":
                    continue
                for alt1 in AA_TO_CODONS.get(aa1, []):
                    for alt2 in AA_TO_CODONS.get(aa2, []):
                        test_list = list(seq)
                        test_list[ci1*3:ci1*3+3] = list(alt1)
                        test_list[ci2*3:ci2*3+3] = list(alt2)
                        test = "".join(test_list)
                        if "ATTTA" not in test:
                            seq = test
                            fixed = True
                            break
                    if fixed:
                        break
                if fixed:
                    break
            if fixed:
                break

        if not fixed:
            break

    return seq


def step_cpg_reconciliation(cfg: "_OptConfig", seq: str) -> str:
    """Aggressive CpG reconciliation pass to eliminate CpG islands.

    This runs AFTER the CAI hill climb and reoptimize steps, because
    those steps can reintroduce CG dinucleotides by preferring high-CAI
    C-ending codons (e.g., GCC for Ala, CCC for Pro) which create
    cross-codon CG when followed by a G-starting codon.

    Strategy:
    1. Identify all windows that fail the CpG island check (Obs/Exp > threshold)
    2. For each failing window, identify all CG dinucleotides
    3. For cross-codon CG (codon ending with C + codon starting with G):
       - Try swapping the C-ending codon to a non-C-ending synonym (sorted by CAI)
       - If the G-starting codon has non-G-starting synonyms, try those too
    4. For within-codon CG:
       - Swap to a CG-free synonym (sorted by CAI)
    5. Accept CAI loss as necessary — CpG island avoidance is a hard constraint
    6. After CG elimination, re-check GC content and adjust if needed
    7. Iterate until all windows pass or no more progress can be made

    Uses IncrementalSequenceState for O(1) GT/CG tracking and
    CodonCache for pre-sorted codon lists.
    """
    state = IncrementalSequenceState(seq)
    codon_cache = CodonCache(cfg.species_cai)
    enzyme_cache = EnzymeSiteCache(cfg.enzymes) if cfg.enzymes else None
    max_iterations = 100

    def _check_gt_safe(ci: int, old_had_gt: dict) -> bool:
        """Check if swap at codon_idx created new avoidable GTs."""
        cs = ci * 3
        ce = cs + 3
        chk_s = max(0, cs - 1)
        chk_e = min(state._n - 1, ce + 1)
        for p in range(chk_s, chk_e):
            if state.has_gt_at(p) and not old_had_gt.get(p, False):
                if not _is_unavoidable_gt(state.sequence, p):
                    return False
        return True

    def _save_gt_region(ci: int) -> dict:
        """Save GT state in the region around a codon."""
        cs = ci * 3
        ce = cs + 3
        chk_s = max(0, cs - 1)
        chk_e = min(state._n - 1, ce + 1)
        return {p: state.has_gt_at(p) for p in range(chk_s, chk_e)}

    for iteration in range(max_iterations):
        current_seq = state.sequence
        old_cg_count = state.cg_count

        # Check if the CpG island predicate passes
        cpg_result = check_no_cpg_island(current_seq, cfg.cpg_window, cfg.cpg_threshold)
        if cpg_result.passed:
            break

        # Find the worst window
        worst_ratio = 0.0
        worst_start = -1
        for start in range(0, len(current_seq) - cfg.cpg_window + 1):
            window_seq = current_seq[start:start + cfg.cpg_window]
            c_count = window_seq.count("C")
            g_count = window_seq.count("G")
            cg_count = _count_dinucs_fast(window_seq, "CG")[0]
            expected = (c_count * g_count) / len(window_seq) if len(window_seq) > 0 else 0
            obs_exp = cg_count / expected if expected > 0 else 0.0
            if obs_exp > worst_ratio:
                worst_ratio = obs_exp
                worst_start = start

        if worst_start < 0 or worst_ratio <= cfg.cpg_threshold:
            break

        # Find CG dinucleotides in the worst window using incremental state
        window_end = min(worst_start + cfg.cpg_window, len(current_seq) - 1)
        fixed = False

        # Collect all CG positions in the worst window with their codon info
        cg_targets = []
        for cg_pos in state.cg_positions_list():
            if worst_start <= cg_pos < window_end:
                c_codon_idx = cg_pos // 3
                g_codon_idx = (cg_pos + 1) // 3
                cg_targets.append((cg_pos, c_codon_idx, g_codon_idx))

        # Sort targets: prioritize within-codon CG first (easier to fix)
        for cg_pos, c_codon_idx, g_codon_idx in cg_targets:
            if fixed:
                break

            if c_codon_idx == g_codon_idx:
                # Within-codon CG — swap the codon to a CG-free alternative
                aa = state.get_aa(c_codon_idx)
                if aa is None or aa == "*":
                    continue

                for alt in codon_cache.get_cg_free_codons(aa):
                    current_codon = state.get_codon(c_codon_idx)
                    if alt == current_codon:
                        continue
                    # Save GT state, apply swap, check constraints
                    old_had_gt = _save_gt_region(c_codon_idx)
                    old_codon = state.swap_codon(c_codon_idx, alt)
                    # Must not reintroduce restriction sites
                    if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                        state.swap_codon(c_codon_idx, old_codon)
                        continue
                    # Must not create new avoidable GTs
                    if not _check_gt_safe(c_codon_idx, old_had_gt):
                        state.swap_codon(c_codon_idx, old_codon)
                        continue
                    # Must reduce net CG count (O(1) check)
                    if state.cg_count < old_cg_count:
                        fixed = True
                        break
                    else:
                        state.swap_codon(c_codon_idx, old_codon)  # Rollback
            else:
                # Cross-codon CG: codon c_codon_idx ends with C, codon g_codon_idx starts with G
                # Strategy 1: Swap the C-ending codon to not end with C
                c_aa = state.get_aa(c_codon_idx)
                if c_aa is not None and c_aa != "*":
                    # Get non-C-ending alternatives from cache, then filter
                    non_c_end_alts = [
                        alt for alt in codon_cache.get_sorted_codons(c_aa)
                        if alt != state.get_codon(c_codon_idx) and alt[-1] != "C" and "CG" not in alt
                    ]
                    for alt in non_c_end_alts:
                        old_had_gt = _save_gt_region(c_codon_idx)
                        old_codon = state.swap_codon(c_codon_idx, alt)
                        # Must not reintroduce restriction sites
                        if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                            state.swap_codon(c_codon_idx, old_codon)
                            continue
                        # Must not create new avoidable GTs
                        if not _check_gt_safe(c_codon_idx, old_had_gt):
                            state.swap_codon(c_codon_idx, old_codon)
                            continue
                        # Must reduce net CG count (O(1) check)
                        if state.cg_count < old_cg_count:
                            fixed = True
                            break
                        else:
                            state.swap_codon(c_codon_idx, old_codon)  # Rollback

                if fixed:
                    break

                # Strategy 2: Swap the G-starting codon to not start with G (if possible)
                g_aa = state.get_aa(g_codon_idx)
                if g_aa is not None and g_aa != "*":
                    # Get non-G-starting alternatives from cache, then filter
                    non_g_start_alts = [
                        alt for alt in codon_cache.get_sorted_codons(g_aa)
                        if alt != state.get_codon(g_codon_idx) and alt[0] != "G" and "CG" not in alt
                    ]
                    for alt in non_g_start_alts:
                        old_had_gt = _save_gt_region(g_codon_idx)
                        old_codon = state.swap_codon(g_codon_idx, alt)
                        # Must not reintroduce restriction sites
                        if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                            state.swap_codon(g_codon_idx, old_codon)
                            continue
                        # Must not create new avoidable GTs
                        if not _check_gt_safe(g_codon_idx, old_had_gt):
                            state.swap_codon(g_codon_idx, old_codon)
                            continue
                        # Must reduce net CG count (O(1) check)
                        if state.cg_count < old_cg_count:
                            fixed = True
                            break
                        else:
                            state.swap_codon(g_codon_idx, old_codon)  # Rollback

                if fixed:
                    break

                # Strategy 3: Paired swap — change BOTH codons
                c_aa = state.get_aa(c_codon_idx)
                g_aa = state.get_aa(g_codon_idx)
                if c_aa is not None and g_aa is not None and c_aa != "*" and g_aa != "*":
                    c_alts = [
                        alt for alt in codon_cache.get_sorted_codons(c_aa)
                        if alt != state.get_codon(c_codon_idx) and alt[-1] != "C" and "CG" not in alt
                    ]
                    g_alts = [
                        alt for alt in codon_cache.get_sorted_codons(g_aa)
                        if alt != state.get_codon(g_codon_idx) and alt[0] != "G" and "CG" not in alt
                    ]
                    for c_alt in c_alts[:TOP_CAI_ALTERNATIVES]:
                        for g_alt in g_alts[:TOP_CAI_ALTERNATIVES]:
                            # Save GT state for both regions
                            old_had_gt_c = _save_gt_region(c_codon_idx)
                            old_had_gt_g = _save_gt_region(g_codon_idx)
                            # Apply both swaps
                            old_c = state.swap_codon(c_codon_idx, c_alt)
                            old_g = state.swap_codon(g_codon_idx, g_alt)
                            # Must not reintroduce restriction sites
                            if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                state.swap_codon(g_codon_idx, old_g)
                                state.swap_codon(c_codon_idx, old_c)
                                continue
                            # Must not create new avoidable GTs in either region
                            if not _check_gt_safe(c_codon_idx, old_had_gt_c) or not _check_gt_safe(g_codon_idx, old_had_gt_g):
                                state.swap_codon(g_codon_idx, old_g)
                                state.swap_codon(c_codon_idx, old_c)
                                continue
                            # Must reduce net CG count (O(1) check)
                            if state.cg_count < old_cg_count:
                                fixed = True
                                break
                            else:
                                state.swap_codon(g_codon_idx, old_g)  # Rollback
                                state.swap_codon(c_codon_idx, old_c)  # Rollback
                        if fixed:
                            break

        if not fixed:
            # No progress — try a more aggressive approach:
            # scan ALL codon positions and replace C-ending codons that
            # create cross-codon CG, regardless of window
            aggressive_fixed = False
            for ci in range(state.num_codons):
                codon = state.get_codon(ci)
                aa = state.get_aa(ci)
                if aa is None or aa == "*":
                    continue
                # Check if this codon ends with C and the next starts with G
                next_ci = ci + 1
                if next_ci < state.num_codons:
                    next_codon = state.get_codon(next_ci)
                    if codon[-1] == "C" and next_codon[0] == "G":
                        # This creates a cross-codon CG — try to fix
                        non_c_end_alts = [
                            alt for alt in codon_cache.get_sorted_codons(aa)
                            if alt != codon and alt[-1] != "C" and "CG" not in alt
                        ]
                        for alt in non_c_end_alts:
                            old_had_gt = _save_gt_region(ci)
                            old_codon = state.swap_codon(ci, alt)
                            if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                                state.swap_codon(ci, old_codon)
                                continue
                            if not _check_gt_safe(ci, old_had_gt):
                                state.swap_codon(ci, old_codon)
                                continue
                            # O(1) CG count check
                            if state.cg_count < old_cg_count:
                                aggressive_fixed = True
                                break
                            else:
                                state.swap_codon(ci, old_codon)  # Rollback

                if not aggressive_fixed and "CG" in codon:
                    # Within-codon CG — try to fix
                    for alt in codon_cache.get_cg_free_codons(aa):
                        if alt == codon:
                            continue
                        old_had_gt = _save_gt_region(ci)
                        old_codon = state.swap_codon(ci, alt)
                        if enzyme_cache is not None and enzyme_cache.check_any_site_present(state.sequence):
                            state.swap_codon(ci, old_codon)
                            continue
                        if not _check_gt_safe(ci, old_had_gt):
                            state.swap_codon(ci, old_codon)
                            continue
                        # O(1) CG count check
                        if state.cg_count < old_cg_count:
                            aggressive_fixed = True
                            break
                        else:
                            state.swap_codon(ci, old_codon)  # Rollback

                if aggressive_fixed:
                    break

            if not aggressive_fixed:
                break  # Truly no more progress possible

    # After CpG reconciliation, re-check GC content and adjust if needed
    result_seq = state.sequence
    gc_val = (result_seq.count("G") + result_seq.count("C")) / max(len(result_seq), 1)
    if not (0.30 <= gc_val <= 0.70):
        # GC drifted — adjust with single-codon swaps that don't reintroduce CGs
        result_seq = fix_gc_after_cpg(cfg, result_seq)

    return result_seq


def cai_recon_check_constraints(cfg: "_OptConfig", old_seq: str, test_seq: str, codon_pos: int) -> bool:
    """Check that a codon swap doesn't violate any hard constraints.

    Hard constraints:
    1. No new cryptic splice sites (MaxEnt score >= threshold)
    2. No restriction enzyme sites
    3. No new avoidable GT dinucleotides

    We allow unavoidable GTs (e.g., Valine codons) but not new avoidable ones.
    For new GTs, we check if they create a cryptic splice site.
    """
    # 1. Check for new avoidable GT dinucleotides and cryptic splice sites
    codon_end = codon_pos + 3
    check_start = max(0, codon_pos - 1)
    check_end = min(len(test_seq) - 1, codon_end + 1)

    for i in range(len(test_seq) - 1):
        if test_seq[i:i+2] == "GT" and old_seq[i:i+2] != "GT":
            # New GT created — is it unavoidable?
            if _is_unavoidable_gt(test_seq, i):
                continue  # Unavoidable GT is OK
            # Avoidable new GT — check if it creates a cryptic splice site
            from ..maxentscan import score_donor as _score_donor
            donor_score = _score_donor(test_seq, i)
            if donor_score >= cfg.splice_low:
                return False  # Would create cryptic splice site

    # 2. Check restriction sites
    if not _check_no_restriction_sites(cfg, test_seq):
        return False

    # 3. Check CpG islands
    cpg_result = check_no_cpg_island(test_seq, cfg.cpg_window, cfg.cpg_threshold)
    if not cpg_result.passed:
        # Check if the CpG failure is new (wasn't there before)
        old_cpg_result = check_no_cpg_island(old_seq, cfg.cpg_window, cfg.cpg_threshold)
        if not old_cpg_result.passed:
            pass  # Was already failing — don't make it worse but don't block
        else:
            # New CpG island created — reject
            return False

    return True


def cai_recon_check_global_constraints(cfg: "_OptConfig", test_seq: str) -> bool:
    """Check that the full sequence satisfies all hard constraints."""
    # 1. Check cryptic splice sites
    from ..maxentscan import max_donor_score as _max_donor, max_acceptor_score as _max_acceptor
    if _max_donor(test_seq) >= cfg.splice_low:
        return False
    if _max_acceptor(test_seq) >= cfg.splice_low:
        return False

    # 2. Check restriction sites
    if not _check_no_restriction_sites(cfg, test_seq):
        return False

    # 3. Check for avoidable GTs
    for i in range(len(test_seq) - 1):
        if test_seq[i:i+2] == "GT":
            if not _is_unavoidable_gt(test_seq, i):
                return False

    return True


def cai_recon_try_paired(cfg: "_OptConfig", seq_list: list, codon_pos: int, new_codon: str) -> Optional[list]:
    """Try a paired codon swap to enable a CAI upgrade.

    When upgrading codon at codon_pos to new_codon creates a cross-codon GT,
    try simultaneously adjusting the adjacent codon to avoid it.
    Returns the new seq_list if successful, None otherwise.
    """
    new_end = new_codon[-1]
    next_pos = codon_pos + 3

    # Check if new codon creates GT with the next codon
    if next_pos + 3 <= len(seq_list):
        next_base = seq_list[next_pos]
        if new_end + next_base == "GT":
            next_codon = "".join(seq_list[next_pos:next_pos + 3])
            next_aa = CODON_TABLE.get(next_codon)
            if next_aa is not None and next_aa != "*":
                for alt2 in sorted(
                    AA_TO_CODONS.get(next_aa, []),
                    key=lambda c: cfg.species_cai.get(c, 0.0),
                    reverse=True,
                ):
                    if new_end + alt2[0] == "GT":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(new_codon):
                        test_list[codon_pos + k] = b
                    for k, b in enumerate(alt2):
                        test_list[next_pos + k] = b
                    test_seq = "".join(test_list)
                    # Check no new avoidable GTs overall
                    old_gt_count = _count_gts("".join(seq_list))
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count <= old_gt_count:
                        if _check_no_restriction_sites(cfg, test_seq):
                            return test_list

    # Check if new codon creates GT with the previous codon
    if codon_pos >= 3:
        prev_end = seq_list[codon_pos - 1]
        if prev_end + new_codon[0] == "GT":
            prev_pos = codon_pos - 3
            prev_codon = "".join(seq_list[prev_pos:prev_pos + 3])
            prev_aa = CODON_TABLE.get(prev_codon)
            if prev_aa is not None and prev_aa != "*":
                for alt0 in sorted(
                    AA_TO_CODONS.get(prev_aa, []),
                    key=lambda c: cfg.species_cai.get(c, 0.0),
                    reverse=True,
                ):
                    if alt0[-1] + new_codon[0] == "GT":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(alt0):
                        test_list[prev_pos + k] = b
                    for k, b in enumerate(new_codon):
                        test_list[codon_pos + k] = b
                    test_seq = "".join(test_list)
                    old_gt_count = _count_gts("".join(seq_list))
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count <= old_gt_count:
                        if _check_no_restriction_sites(cfg, test_seq):
                            return test_list

    return None


def cai_recon_try_paired_state(
    state: IncrementalSequenceState, codon_idx: int,
    new_codon: str, codon_cache: CodonCache,
    enzyme_cache: Optional[EnzymeSiteCache],
) -> Tuple[bool, Optional[Tuple[int, str]]]:
    """Try a paired codon swap to enable a CAI upgrade (incremental version).

    The main codon at codon_idx is already swapped to new_codon in state.
    Try adjusting an adjacent codon to resolve GT conflicts.

    Returns (success, adj_rollback) where adj_rollback is (adj_idx, old_adj_codon)
    if a paired swap was applied and needs to be tracked for potential rollback.
    If success=False, the adjacent swap (if any) has already been rolled back;
    only the main swap needs rolling back by the caller.
    """
    new_end = new_codon[-1]
    next_idx = codon_idx + 1

    # Check if new codon creates GT with the next codon
    if next_idx < state.num_codons:
        next_codon = state.get_codon(next_idx)
        next_aa = state.get_aa(next_idx)
        if new_end + next_codon[0] == "GT" and next_aa is not None and next_aa != "*":
            old_gt_count = state.gt_count
            for alt2 in codon_cache.get_sorted_codons(next_aa):
                if new_end + alt2[0] == "GT":
                    continue
                # O(1) swap + O(1) GT count check + rollback
                old_adj = state.swap_codon(next_idx, alt2)
                if state.gt_count <= old_gt_count:
                    # Check restriction sites
                    rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                             if enzyme_cache else True)
                    if rs_ok:
                        return (True, (next_idx, old_adj))
                state.swap_codon(next_idx, old_adj)  # Rollback adjacent

    # Check if new codon creates GT with the previous codon
    prev_idx = codon_idx - 1
    if prev_idx >= 0:
        prev_codon = state.get_codon(prev_idx)
        prev_aa = state.get_aa(prev_idx)
        if prev_codon[-1] + new_codon[0] == "GT" and prev_aa is not None and prev_aa != "*":
            old_gt_count = state.gt_count
            for alt0 in codon_cache.get_sorted_codons(prev_aa):
                if alt0[-1] + new_codon[0] == "GT":
                    continue
                old_adj = state.swap_codon(prev_idx, alt0)
                if state.gt_count <= old_gt_count:
                    rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                             if enzyme_cache else True)
                    if rs_ok:
                        return (True, (prev_idx, old_adj))
                state.swap_codon(prev_idx, old_adj)  # Rollback adjacent

    return (False, None)


def cai_recon_rollback_paired(
    state: IncrementalSequenceState, codon_idx: int, old_codon: str,
) -> None:
    """Rollback a failed paired swap in _step_cai_reconciliation.

    This is called when a paired swap passes local constraints but fails
    global constraints. It rolls back the main codon swap.
    Note: the adjacent swap should have already been rolled back by
    _cai_recon_try_paired_state when it returned False, or should be
    handled separately when it returned True but global check failed.
    """
    state.swap_codon(codon_idx, old_codon)


def step_cai_reconciliation(cfg: "_OptConfig", seq: str) -> str:
    """CAI Reconciliation pass: upgrade low-CAI codons while maintaining hard constraints.

    This pass runs AFTER all hard constraints have been satisfied (cryptic splice,
    restriction sites, CpG islands). It greedily upgrades each low-CAI codon to the
    best synonymous alternative that doesn't violate any hard constraint.

    Algorithm:
    1. Identify all codon positions where CAI < min_cai (or below a working threshold)
    2. Sort positions by CAI improvement potential (best possible CAI minus current CAI)
    3. For each position, try upgrading to the highest-CAI synonymous codon that:
       a. Doesn't create a cryptic splice site (MaxEnt score < threshold)
       b. Doesn't create a restriction enzyme site
       c. Doesn't create a CpG island
    4. If direct upgrade creates a new avoidable GT, try paired codon swaps
       (upgrade the target + modify an adjacent codon to avoid the GT)
    5. Iterate until no more upgrades are possible

    Key insight: hard constraints (no cryptic splice, no restriction sites) are
    satisfied FIRST, then CAI is maximized SUBJECT to those constraints.

    Uses IncrementalSequenceState for O(1) GT/CG tracking and
    CodonCache for pre-sorted codon lists.
    """
    import math
    state = IncrementalSequenceState(seq)
    codon_cache = CodonCache(cfg.species_cai)
    enzyme_cache = EnzymeSiteCache(cfg.enzymes) if cfg.enzymes else None
    max_iterations = 30

    for iteration in range(max_iterations):
        current_seq = state.sequence
        any_upgrade = False

        # 1. Collect positions where CAI is low, with improvement potential
        low_cai_positions = []
        for codon_idx in range(state.num_codons):
            codon = state.get_codon(codon_idx)
            aa = state.get_aa(codon_idx)
            if aa is None or aa == "*":
                continue
            current_cai = codon_cache.get_cai(codon)
            if current_cai >= cfg.min_cai:
                continue  # Already above threshold
            # Compute best possible CAI for this AA
            best_possible = codon_cache.get_cai(codon_cache.get_best_codon(aa))
            improvement = best_possible - current_cai
            if improvement > 0:
                pos = codon_idx * 3
                low_cai_positions.append((pos, codon, aa, current_cai, best_possible, improvement))

        if not low_cai_positions:
            break

        # 2. Sort by improvement potential (descending) — fix positions with
        #    the most room for improvement first
        low_cai_positions.sort(key=lambda x: x[5], reverse=True)

        # 3. Try to upgrade each position
        for pos, codon, aa, current_cai, best_possible, improvement in low_cai_positions:
            codon_idx = pos // 3
            # Get synonymous codons sorted by CAI (highest first) from cache
            candidates = codon_cache.get_sorted_codons(aa)

            upgraded = False
            for alt in candidates:
                alt_cai = codon_cache.get_cai(alt)
                if alt_cai <= current_cai:
                    break  # No improvement possible (candidates are sorted)

                # O(1) swap + check + rollback pattern
                old_codon = state.swap_codon(codon_idx, alt)
                test_seq = state.sequence

                # Check hard constraints
                if not cai_recon_check_constraints(cfg, current_seq, test_seq, pos):
                    # Direct upgrade fails — try paired codon swap
                    # State already has the main swap applied; try adjacent swaps
                    paired_ok, adj_rollback = cai_recon_try_paired_state(
                        state, codon_idx, alt, codon_cache, enzyme_cache
                    )
                    if paired_ok:
                        # Verify the paired swap satisfies global constraints
                        if cai_recon_check_global_constraints(cfg, state.sequence):
                            any_upgrade = True
                            upgraded = True
                            break
                        else:
                            # Rollback the adjacent swap first, then main swap
                            if adj_rollback is not None:
                                state.swap_codon(adj_rollback[0], adj_rollback[1])
                    # Rollback the main swap
                    state.swap_codon(codon_idx, old_codon)
                    continue

                # Safe upgrade (already applied via swap_codon)
                any_upgrade = True
                upgraded = True
                break

            if upgraded:
                break  # Restart with updated sequence

        if not any_upgrade:
            break

    return state.sequence


def step_gt_reconciliation(cfg: "_OptConfig", seq: str) -> str:
    """Fix avoidable GT dinucleotides that may have been introduced by
    the ATTTA removal or CpG reconciliation steps.

    Strategy: For each avoidable GT, try synonymous codon swaps that
    eliminate the GT without reintroducing ATTTA, CpG islands, or
    restriction sites.

    Uses IncrementalSequenceState for O(1) GT tracking and
    CodonCache for pre-sorted codon lists.
    """
    state = IncrementalSequenceState(seq)
    codon_cache = CodonCache(cfg.species_cai)
    enzyme_cache = EnzymeSiteCache(cfg.enzymes) if cfg.enzymes else None
    initial_cg_count = state.cg_count  # Track initial CG count for cheap CpG check

    def _get_avoidable_gt_positions():
        """Get avoidable GT positions from incremental state (O(K) where K = num GTs)."""
        return [p for p in state.gt_positions_list()
                if not _is_unavoidable_gt(state.sequence, p)]

    for iteration in range(50):
        avoidable_positions = _get_avoidable_gt_positions()
        if not avoidable_positions:
            break

        n_avoidable = len(avoidable_positions)
        fixed = False
        for pos in avoidable_positions:
            n_codons = state.num_codons
            # Try codons overlapping the GT
            for offset in [0, 1]:
                ci = (pos + offset) // 3
                if ci < 0 or ci >= n_codons:
                    continue
                aa = state.get_aa(ci)
                if not aa or aa == "*":
                    continue
                current = state.get_codon(ci)
                # Try GT-free alternatives first, then all alternatives sorted by CAI
                gt_free = [c for c in codon_cache.get_gt_free_codons(aa) if c != current]
                other_alts = [
                    c for c in codon_cache.get_sorted_codons(aa)
                    if "GT" in c and c != current
                ]
                alternatives = gt_free + other_alts

                for alt in alternatives:
                    # O(1) predictive boundary check before swap
                    if cfg.avoid_gt:
                        left_gt, right_gt = state.boundary_creates_gt(ci, alt)
                        if left_gt or right_gt:
                            continue
                    # O(1) swap + check + rollback
                    old_codon = state.swap_codon(ci, alt)
                    # Quick check: did GT count actually decrease? O(1)
                    if state.gt_count >= n_avoidable + (n_avoidable - len(avoidable_positions)):
                        # More precise: count avoidable GTs
                        new_avoidable = _get_avoidable_gt_positions()
                        if len(new_avoidable) >= n_avoidable:
                            state.swap_codon(ci, old_codon)
                            continue
                    else:
                        new_avoidable = _get_avoidable_gt_positions()
                        if len(new_avoidable) >= n_avoidable:
                            state.swap_codon(ci, old_codon)
                            continue
                    # Check no ATTTA reintroduced
                    if "ATTTA" in state.sequence:
                        state.swap_codon(ci, old_codon)
                        continue
                    # Check no new restriction sites (O(E*L) with cache)
                    rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                             if enzyme_cache else True)
                    if not rs_ok:
                        state.swap_codon(ci, old_codon)
                        continue
                    # Quick CG check: only do expensive CpG island check if CG count increased
                    if state.cg_count > initial_cg_count:
                        # Full CpG island check only if CG actually increased
                        cpg_result = check_no_cpg_island(state.sequence, cfg.cpg_window, cfg.cpg_threshold)
                        if not cpg_result.passed:
                            state.swap_codon(ci, old_codon)
                            continue
                    fixed = True
                    break
                if fixed:
                    break
            if fixed:
                break

        if not fixed:
            # Try 2-codon coordinated swap for cross-codon GTs
            for pos in avoidable_positions:
                left_ci = pos // 3
                right_ci = (pos + 1) // 3
                if left_ci == right_ci or left_ci < 0 or right_ci < 0:
                    continue
                n_codons = state.num_codons
                if left_ci >= n_codons or right_ci >= n_codons:
                    continue

                left_aa = state.get_aa(left_ci)
                right_aa = state.get_aa(right_ci)
                if not left_aa or not right_aa or left_aa == "*" or right_aa == "*":
                    continue

                for left_alt in codon_cache.get_sorted_codons(left_aa):
                    for right_alt in codon_cache.get_sorted_codons(right_aa):
                        # Skip if GT still at boundary
                        if left_alt[-1] == "G" and right_alt[0] == "T":
                            continue
                        # O(1) paired swap + check + rollback
                        old_left = state.swap_codon(left_ci, left_alt)
                        old_right = state.swap_codon(right_ci, right_alt)

                        if "ATTTA" in state.sequence:
                            state.swap_codon(right_ci, old_right)
                            state.swap_codon(left_ci, old_left)
                            continue
                        # Check GT improved
                        new_avoidable = _get_avoidable_gt_positions()
                        if len(new_avoidable) >= n_avoidable:
                            state.swap_codon(right_ci, old_right)
                            state.swap_codon(left_ci, old_left)
                            continue
                        # Check RS and CpG
                        rs_ok = (not enzyme_cache.check_any_site_present(state.sequence)
                                 if enzyme_cache else True)
                        if not rs_ok:
                            state.swap_codon(right_ci, old_right)
                            state.swap_codon(left_ci, old_left)
                            continue
                        # Quick CG check: only expensive CpG island check if CG increased
                        if state.cg_count > initial_cg_count:
                            cpg_result = check_no_cpg_island(state.sequence, cfg.cpg_window, cfg.cpg_threshold)
                            if not cpg_result.passed:
                                state.swap_codon(right_ci, old_right)
                                state.swap_codon(left_ci, old_left)
                                continue
                        fixed = True
                        break
                    if fixed:
                        break
                if fixed:
                    break
            if not fixed:
                break

    return state.sequence


def try_paired_cai_upgrade(
    cfg: "_OptConfig", seq_list: list, codon_pos: int, new_codon: str, old_gt_count: int
) -> Optional[list]:
    """Try a paired codon swap to enable a CAI upgrade.

    When upgrading codon at codon_pos to new_codon creates a cross-codon GT,
    try simultaneously adjusting the adjacent codon to avoid it.
    Returns the new seq_list if successful, None otherwise.
    """
    new_end = new_codon[-1]  # Last base of the new codon
    next_pos = codon_pos + 3

    # Check if new codon creates GT with the next codon
    if next_pos + 3 <= len(seq_list):
        next_base = seq_list[next_pos]
        if new_end + next_base == "GT":
            # Try to fix by changing the next codon
            next_codon = "".join(seq_list[next_pos:next_pos + 3])
            next_aa = CODON_TABLE.get(next_codon)
            if next_aa is not None and next_aa != "*":
                for alt2 in sorted(
                    AA_TO_CODONS.get(next_aa, []),
                    key=lambda c: cfg.species_cai.get(c, 0.0),
                    reverse=True,
                ):
                    if new_end + alt2[0] == "GT":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(new_codon):
                        test_list[codon_pos + k] = b
                    for k, b in enumerate(alt2):
                        test_list[next_pos + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count <= old_gt_count:
                        # Check restriction sites
                        from ..restriction_sites import get_recognition_site
                        rs_ok = True
                        for enzyme in cfg.enzymes:
                            site = get_recognition_site(enzyme)
                            if site is None:
                                continue
                            if site in test_seq:
                                rs_ok = False
                                break
                        if rs_ok:
                            return test_list

    # Check if new codon creates GT with the previous codon
    if codon_pos >= 3:
        prev_end = seq_list[codon_pos - 1]
        if prev_end + new_codon[0] == "GT":
            # Try to fix by changing the previous codon
            prev_pos = codon_pos - 3
            prev_codon = "".join(seq_list[prev_pos:prev_pos + 3])
            prev_aa = CODON_TABLE.get(prev_codon)
            if prev_aa is not None and prev_aa != "*":
                for alt0 in sorted(
                    AA_TO_CODONS.get(prev_aa, []),
                    key=lambda c: cfg.species_cai.get(c, 0.0),
                    reverse=True,
                ):
                    if alt0[-1] + new_codon[0] == "GT":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(alt0):
                        test_list[prev_pos + k] = b
                    for k, b in enumerate(new_codon):
                        test_list[codon_pos + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count <= old_gt_count:
                        from ..restriction_sites import get_recognition_site
                        rs_ok = True
                        for enzyme in cfg.enzymes:
                            site = get_recognition_site(enzyme)
                            if site is None:
                                continue
                            if site in test_seq:
                                rs_ok = False
                                break
                        if rs_ok:
                            return test_list

    return None


def try_paired_cai_upgrade_incremental(
    state: IncrementalSequenceState, codon_idx: int,
    new_codon: str, codon_cache: CodonCache,
    enzyme_cache: Optional[EnzymeSiteCache] = None
) -> bool:
    """Try a paired codon swap using IncrementalSequenceState.

    When upgrading a codon would create a cross-codon GT, this tries
    simultaneously adjusting the adjacent codon to avoid it.
    Returns True if a paired swap was applied.
    """
    start = codon_idx * 3
    # The new GT must be at the right boundary: new_codon[-1] + next_base == GT
    # or at the left boundary: prev_base + new_codon[0] == GT
    # We need to find which adjacent codon to adjust

    next_start = start + 3

    # Try adjusting the next codon
    if next_start + 3 <= len(state.sequence):
        next_codon = state.get_codon(codon_idx + 1)
        next_aa = state.get_aa(codon_idx + 1)
        if next_aa is not None and next_aa != "*":
            for next_alt in codon_cache.get_sorted_codons(next_aa):
                if next_alt == next_codon:
                    continue
                # Check that new_codon[-1] + next_alt[0] != "GT"
                if new_codon[-1] + next_alt[0] == "GT":
                    continue
                # Try the paired swap
                old_codon = state.swap_codon(codon_idx, new_codon)
                old_next = state.swap_codon(codon_idx + 1, next_alt)

                # Check that GT count didn't increase beyond what we expected
                if not state.would_gt_increase(codon_idx, new_codon):
                    # Also check CG didn't increase
                    if not state.would_cg_increase(codon_idx, new_codon):
                        return True
                # Rollback
                state.swap_codon(codon_idx + 1, old_next)
                state.swap_codon(codon_idx, old_codon)

    # Try adjusting the previous codon
    if codon_idx > 0:
        prev_codon = state.get_codon(codon_idx - 1)
        prev_aa = state.get_aa(codon_idx - 1)
        if prev_aa is not None and prev_aa != "*":
            for prev_alt in codon_cache.get_sorted_codons(prev_aa):
                if prev_alt == prev_codon:
                    continue
                if prev_alt[-1] + new_codon[0] == "GT":
                    continue
                old_prev = state.swap_codon(codon_idx - 1, prev_alt)
                old_codon = state.swap_codon(codon_idx, new_codon)

                if not state.would_gt_increase(codon_idx, new_codon):
                    if not state.would_cg_increase(codon_idx, new_codon):
                        return True
                state.swap_codon(codon_idx, old_codon)
                state.swap_codon(codon_idx - 1, old_prev)

    return False


def step_cai_hill_climb(cfg: "_OptConfig", seq: str) -> str:
    """CAI Hill Climb step: CAI hill climbing (optimized with incremental state).

    For each codon position, try upgrading to a higher-CAI synonym
    if it doesn't introduce new constraint violations (GT, RS, etc.).
    This is the key step that recovers CAI lost during constraint fixing.

    Uses IncrementalSequenceState for O(1) GT/CG tracking instead of
    O(N) full-sequence rescans, and CodonCache for pre-sorted codon lists.

    When NUMBA is available, the ``batch_codon_swap_score`` kernel scores
    all candidate swaps at a position in a single vectorized pass, then
    constraint checks (GT, CG, RS) are applied only to improving candidates.
    Falls back to per-position Python evaluation when NUMBA is absent.
    """
    state = IncrementalSequenceState(seq)
    codon_cache = CodonCache(cfg.species_cai)
    enzyme_cache = EnzymeSiteCache(cfg.enzymes) if cfg.enzymes else None
    max_iterations = 10

    # ── Batch scorer (NUMBA-accelerated when available) ───────
    batch_scorer = _BatchSwapScorer(cfg.species_cai)

    # Initialize incremental CAI tracking for O(1) log_sum updates
    seq_codons_init = [state.get_codon(i) for i in range(state.num_codons)]
    batch_scorer.reset_incremental_state(seq_codons_init)

    for iteration in range(max_iterations):
        any_upgrade = False

        for codon_idx in range(state.num_codons):
            codon = state.get_codon(codon_idx)
            aa = state.get_aa(codon_idx)
            if aa is None or aa == "*":
                continue

            current_cai = codon_cache.get_cai(codon)

            # Get synonymous codons sorted by CAI (highest first) — from cache
            candidates = codon_cache.get_sorted_codons(aa)

            # ── Batch CAI scoring ──────────────────────────────
            # Use the NUMBA kernel (or its Python fallback) to score
            # all candidate swaps at once.  Then iterate over
            # improving candidates in CAI order for constraint checks.
            improving_candidates = [
                (alt, alt_cai)
                for alt in candidates
                if (alt_cai := codon_cache.get_cai(alt)) > current_cai
            ]

            if not improving_candidates:
                continue

            # Batch-score only the improving candidates
            improving_alts = [alt for alt, _ in improving_candidates]
            improving_alts_only = [alt for alt, _ in improving_candidates]

            # Build seq_codons for the batch scorer
            seq_codons = [
                state.get_codon(i) for i in range(state.num_codons)
            ]
            cai_scores = batch_scorer.score_candidates(
                seq_codons, codon_idx, improving_alts_only,
            )

            # Iterate over improving candidates (highest CAI first)
            for k, (alt, _) in enumerate(improving_candidates):
                # Quick O(1) boundary check: would this swap increase GT count?
                if cfg.avoid_gt and state.would_gt_increase(codon_idx, alt):
                    # Check if ALL new GTs are unavoidable
                    new_gt_positions = state.new_gt_positions_after_swap(codon_idx, alt)
                    if new_gt_positions:
                        # Temporarily swap to get the full sequence for _is_unavoidable_gt
                        old_codon = state.swap_codon(codon_idx, alt)
                        full_seq = state.sequence
                        all_unavoidable = all(
                            _is_unavoidable_gt(full_seq, pos) for pos in new_gt_positions
                        )
                        if not all_unavoidable:
                            # Try paired codon swap to fix the new avoidable cross-codon GT
                            if len(new_gt_positions) == 1:
                                paired = try_paired_cai_upgrade_incremental(
                                    state, codon_idx, alt, codon_cache, enzyme_cache
                                )
                                if paired:
                                    batch_scorer.update_incremental_state(codon, alt)
                                    any_upgrade = True
                                    break
                            # Rollback
                            state.swap_codon(codon_idx, old_codon)
                            continue
                        # All new GTs are unavoidable — accept the upgrade
                        batch_scorer.update_incremental_state(old_codon, alt)
                        any_upgrade = True
                        break
                    # No new GTs despite would_gt_increase returning True? Shouldn't happen, but accept

                # Quick O(1) check: would this swap increase CG count?
                if state.would_cg_increase(codon_idx, alt):
                    continue

                # Check restriction sites — only in the affected region
                if enzyme_cache is not None:
                    start = codon_idx * 3
                    old_codon = state.swap_codon(codon_idx, alt)
                    # Regional RS check (much faster than full-sequence scan)
                    rs_ok = not enzyme_cache.check_sites_in_region(
                        state.sequence, start, start + 3
                    )
                    if not rs_ok:
                        # Full-sequence fallback for safety
                        rs_ok = not enzyme_cache.check_any_site_present(state.sequence)
                    if not rs_ok:
                        state.swap_codon(codon_idx, old_codon)  # Rollback
                        continue
                    batch_scorer.update_incremental_state(old_codon, alt)
                    any_upgrade = True
                    break
                else:
                    # No enzymes to check — apply the upgrade
                    old_codon = state.swap_codon(codon_idx, alt)
                    batch_scorer.update_incremental_state(old_codon, alt)
                    any_upgrade = True
                    break

        if not any_upgrade:
            break

    return state.sequence


def step_reoptimize(cfg: "_OptConfig", seq: str) -> str:
    """Re-optimization step: Iterative re-optimization pass (optimized with incremental state).

    Repeats until no more improvements can be made:
    1. Per-codon CAI optimization with GT avoidance
    2. Cross-codon GT resolution
    3. Restriction site removal

    Uses IncrementalSequenceState for O(1) GT/CG tracking instead of
    O(N) full-sequence rescans.
    """
    state = IncrementalSequenceState(seq)
    codon_cache = CodonCache(cfg.species_cai)
    enzyme_cache = EnzymeSiteCache(cfg.enzymes) if cfg.enzymes else None
    max_iterations = 20

    for iteration in range(max_iterations):
        old_gt_count = state.gt_count
        improved = False

        # Step 1: Per-codon optimization - try to swap to GT-free codons
        for codon_idx in range(state.num_codons):
            codon = state.get_codon(codon_idx)
            aa = state.get_aa(codon_idx)
            if aa is None or aa == "*":
                continue

            candidates = codon_cache.get_sorted_codons(aa)
            if not candidates:
                continue

            # If current codon has GT, try to swap
            if "GT" in codon:
                for alt in candidates:
                    if "GT" in alt:
                        continue
                    # O(1) boundary check
                    left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
                    if left_gt or right_gt:
                        continue
                    # O(1) GT count change check
                    if not state.would_gt_increase(codon_idx, alt):
                        state.swap_codon(codon_idx, alt)
                        improved = True
                        break

        # Step 2: Cross-codon GT resolution
        current_seq = state.sequence
        for pos in find_cross_codon_gt(current_seq):
            codon_idx = pos // 3
            next_codon_idx = codon_idx + 1

            if next_codon_idx >= state.num_codons:
                continue

            aa1 = state.get_aa(codon_idx)
            aa2 = state.get_aa(next_codon_idx)
            if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                continue

            for c1 in codon_cache.get_sorted_codons(aa1):
                if "GT" in c1:
                    continue
                for c2 in codon_cache.get_sorted_codons(aa2):
                    if "GT" in c2:
                        continue
                    if c1[-1] + c2[0] == "GT":
                        continue
                    
                    # Try the paired swap with O(1) tracking
                    old1 = state.swap_codon(codon_idx, c1)
                    old2 = state.swap_codon(next_codon_idx, c2)
                    
                    if state.gt_count < old_gt_count:
                        improved = True
                        break
                    else:
                        # Rollback
                        state.swap_codon(next_codon_idx, old2)
                        state.swap_codon(codon_idx, old1)
            
            if improved:
                break

        if not improved and codon_idx > 0:
            prev_codon_idx = codon_idx - 1
            aa0 = state.get_aa(prev_codon_idx)
            aa1 = state.get_aa(codon_idx)
            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                for c0 in codon_cache.get_sorted_codons(aa0):
                    if "GT" in c0:
                        continue
                    for c1 in codon_cache.get_sorted_codons(aa1):
                        if "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        
                        old0 = state.swap_codon(prev_codon_idx, c0)
                        old1_codon = state.swap_codon(codon_idx, c1)
                        
                        if state.gt_count < old_gt_count:
                            improved = True
                            break
                        else:
                            state.swap_codon(codon_idx, old1_codon)
                            state.swap_codon(prev_codon_idx, old0)
                    
                    if improved:
                        break

        # Step 3: Restriction site removal
        if enzyme_cache is not None:
            current_seq = state.sequence
            for enzyme, site in enzyme_cache._sites.items():
                if site is None:
                    continue
                p = current_seq.find(site)
                while p != -1:
                    codon_starts = set()
                    for j in range(p, p + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(current_seq):
                            codon_starts.add(cs)
                    rs_resolved = False
                    for cs in sorted(codon_starts):
                        ci = cs // 3
                        codon = state.get_codon(ci)
                        aa = state.get_aa(ci)
                        if aa is None or aa == "*":
                            continue
                        for alt in codon_cache.get_sorted_codons(aa):
                            if alt == codon:
                                continue
                            old_codon = state.swap_codon(ci, alt)
                            if site not in state.sequence:
                                if not cfg.avoid_gt or state.gt_count <= old_gt_count:
                                    rs_resolved = True
                                    improved = True
                                    break
                            state.swap_codon(ci, old_codon)  # Rollback
                        if rs_resolved:
                            break
                    if rs_resolved:
                        current_seq = state.sequence
                        p = current_seq.find(site)
                    else:
                        p = current_seq.find(site, p + 1)

        if not improved:
            break

    return state.sequence

# Deprecated alias — use step_reoptimize instead
_phase7_reoptimize = step_reoptimize


def step_mrna_stability_improvement(cfg: "_OptConfig", seq: str) -> str:
    """Soft mRNA stability improvement pass.

    Identifies destabilizing motifs (ATTTA, extended AREs, long A/T runs)
    and applies synonymous codon changes to remove them, WITHOUT breaking
    any hard constraints (restriction sites, GC range, stop codons,
    cryptic splice sites).

    This step runs after all other optimization steps, so hard constraints
    are already satisfied.  We only apply changes that preserve all
    existing constraint satisfaction.

    Returns (seq, mrna_stability_score, motifs_removed, stability_improvement)
    for the caller to store on the optimizer instance.
    """
    if not cfg.optimize_mrna_stability:
        return (seq, None, 0, None)

    from ..mrna_stability import score_mrna_stability, suggest_mutations_for_stability
    from ..restriction_sites import get_recognition_site as _get_site
    from ..type_system import check_no_stop_codons as _check_stops

    # Map species key back to organism name for the stability module
    species_to_organism = {
        "human": "Homo_sapiens",
        "ecoli": "Escherichia_coli",
        "mouse": "Mus_musculus",
        "cho": "CHO_K1",
        "yeast": "Saccharomyces_cerevisiae",
    }
    organism = species_to_organism.get(cfg.species, "Homo_sapiens")

    initial_stability = score_mrna_stability(seq, organism)
    suggestions = suggest_mutations_for_stability(seq, organism)

    motifs_removed = 0
    for suggestion in suggestions:
        pos = suggestion['position']
        # The mrna_stability module returns 'position' as 0-based
        # nucleotide index of the codon start.  Use it directly.
        codon_start = pos
        new_codon = suggestion['suggested_codon']

        # Safety check: ensure codon start is in range
        if codon_start + 3 > len(seq):
            continue

        # Apply the change and verify no hard constraints are broken
        test_seq = seq[:codon_start] + new_codon + seq[codon_start + 3:]

        # Verify: same protein (synonymous)
        test_protein = _validation.translate(test_seq)
        orig_protein = _validation.translate(seq)
        if test_protein != orig_protein:
            continue

        # Verify: GC still in reasonable range (0.30–0.70)
        test_gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
        if not (0.30 <= test_gc <= 0.70):
            continue

        # Verify: no new restriction sites introduced
        rs_violated = False
        for enz in cfg.enzymes:
            site = _get_site(enz)
            if site and site in test_seq and site not in seq:
                rs_violated = True
                break
        if rs_violated:
            continue

        # Verify: no new stop codons
        stop_result = _check_stops(test_seq)
        if not stop_result.passed:
            continue

        # Verify: no worsening of cryptic splice scores
        try:
            from ..maxentscan import max_donor_score, max_acceptor_score
            old_max_d = max_donor_score(seq)
            old_max_a = max_acceptor_score(seq)
            new_max_d = max_donor_score(test_seq)
            new_max_a = max_acceptor_score(test_seq)
            if new_max_d > old_max_d + 0.1 or new_max_a > old_max_a + 0.1:
                continue
        except (ImportError, Exception):
            pass  # If MaxEntScan unavailable, skip splice check

        # Safe to apply
        seq = test_seq
        motifs_removed += 1
        logger.debug(
            "mRNA stability: replaced %s->%s at codon %d (removing %s)",
            suggestion.get('original_codon', '???'), new_codon,
            codon_start // 3,
            suggestion.get('motif_removed', suggestion.get('motif', '???')),
        )

    final_stability = score_mrna_stability(seq, organism)
    # score_mrna_stability returns an MRNAStabilityScore dataclass;
    # extract the float overall_score
    init_score = initial_stability.overall_score if hasattr(initial_stability, 'overall_score') else float(initial_stability)
    final_score = final_stability.overall_score if hasattr(final_stability, 'overall_score') else float(final_stability)
    mrna_stability_score = final_score
    # _destabilizing_motifs_removed was already incremented in the loop above
    if init_score > 0:
        stability_improvement = round(final_score - init_score, 6)
    else:
        stability_improvement = round(final_score, 6)

    logger.info(
        "mRNA stability: initial=%.4f final=%.4f improvement=%.4f motifs_removed=%d",
        init_score, final_score, stability_improvement,
        motifs_removed,
    )

    return (seq, mrna_stability_score, motifs_removed, stability_improvement)


# ── Deprecated aliases (for backward compatibility) ──────────────
_phase5_avoid_cpg_islands = step_avoid_cpg_islands
_phase6_cai_hill_climb = step_cai_hill_climb
_phase7_reoptimize = step_reoptimize
