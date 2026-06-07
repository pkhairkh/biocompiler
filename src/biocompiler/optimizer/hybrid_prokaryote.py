"""
Prokaryotic optimization path for HybridOptimizer.

Extracted from hybrid_optimizer.py to decompose the monolith.
Contains _optimize_prokaryote_fast() and its helpers.
"""

from __future__ import annotations

import math
import logging
from typing import Any

from ..type_system import CODON_TABLE, AA_TO_CODONS
from ..decision_provenance import CodonDecision, ConstraintDecision
from .hybrid_types import HybridResult, GT_CAI_COST_THRESHOLD as _GT_CAI_COST_THRESHOLD

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from ..numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        count_gc as _numba_count_gc,
        compute_cai_kernel as _numba_cai_kernel,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False

logger = logging.getLogger(__name__)


def _optimize_prokaryote_fast(optimizer, protein: str) -> HybridResult:
    """Ultra-fast optimization path for prokaryotic organisms.

    Prokaryotes (E. coli, B. subtilis, etc.) have no spliceosome,
    so MaxEntScan splice scoring is biologically irrelevant.  CpG
    island disruption and GT/AG avoidance are also meaningless.

    This method exploits those facts to skip ALL eukaryote-specific
    machinery and use a streamlined single-pass approach:

    1. Greedy init with highest-CAI codons + incremental GC/CAI tracking
    2. Fix restriction sites (single scan, local window checks)
    3. Fix GC content (incremental tracking, single-codon swaps)
    4. Fix ATTTA motifs (fast find + codon swap)
    5. Fix T-runs (single scan + codon swap)

    Key performance wins over the general path:
    - No IncrementalSequenceState object (saves ~0.3ms of GT/CG/AG setup)
    - No Violation objects or heapq operations
    - No Phase 2b (MaxEntScan) at all
    - No Phase 3 (CAI hill climbing) — start with max CAI
    - Incremental GC tracking during construction (no .count() scan)
    - Precomputed max adaptiveness for fast CAI computation
    - Local window checks instead of full-sequence scans after swaps
    - Single string scan for all pattern violations simultaneously

    Target: <0.5ms for GFP (714bp) in E. coli

    Args:
        protein: Amino acid sequence (single-letter codes, no stop).

    Returns:
        HybridResult with the optimized DNA sequence and metrics.
    """
    # NOTE: Prokaryotes skip GT avoidance entirely.  Prokaryotes have no
    # spliceosome, so GT dinucleotide avoidance is biologically irrelevant.
    # This is enforced by optimizer.avoid_gt=False (set in __init__ for
    # prokaryotes) and by the is_prokaryote flag.  The GT-related code
    # paths (gt_free, ag_free, _fix_avoidable_gt, etc.) are never
    # reached on this fast path.
    species_cai = optimizer.species_cai
    max_adapt = optimizer._max_adapt
    optimal_codon = optimizer.optimal_codon
    sorted_codons = optimizer.sorted_codons
    codon_gc = optimizer.codon_gc
    rs_sites = optimizer._rs_sites
    ac_scanner = optimizer._ac_scanner
    gc_lo = optimizer.gc_lo
    gc_hi = optimizer.gc_hi
    max_rs_len = optimizer._max_rs_site_len

    # ── Provenance collector reference ──
    _prov = optimizer.provenance_collector

    # ── Phase 1: Greedy init with incremental GC + CAI tracking ──
    # Build sequence and track GC + CAI log product simultaneously
    codon_list: list[str] = []
    gc_count = 0
    log_cai_sum = 0.0
    n_cai_codons = 0

    for pos, aa in enumerate(protein):
        if aa == "*":
            codon_list.append("TAA")
            continue
        best = optimal_codon.get(aa)
        if not best:
            sl = sorted_codons.get(aa, [])
            best = sl[0] if sl else "NNN"
        codon_list.append(best)
        gc_count += codon_gc.get(best, 0)
        # Incremental CAI
        adapt = species_cai.get(best, 0.0)
        max_a = max_adapt.get(aa, 0.0)
        if max_a > 0 and adapt > 0:
            log_cai_sum += math.log(adapt / max_a)
            n_cai_codons += 1

        # Record codon decision for provenance
        if _prov is not None:
            chosen_cai = adapt
            alternatives: list[dict[str, Any]] = []
            for codon in sorted_codons.get(aa, []):
                if codon == best:
                    continue
                cai_val = species_cai.get(codon, 0.0)
                gc_bases = sum(1 for b in codon if b in "GC")
                # Prokaryotes: skip eukaryotic-only constraints
                # (splice, CpG, GT dinucleotide) — these are
                # biologically irrelevant for organisms without
                # a spliceosome.
                violates: list[str] = []
                rejected_because: str | None = None
                if cai_val < chosen_cai:
                    rejected_because = "Lower CAI"
                else:
                    rejected_because = "Lower CAI"
                alternatives.append({
                    "codon": codon,
                    "cai_contribution": round(cai_val, 4),
                    "gc_contribution": round(gc_bases / 3.0, 2),
                    "violates_constraints": violates,
                    "rejected_because": rejected_because,
                })
            # Compute confidence
            sl = sorted_codons.get(aa, [])
            if len(sl) > 1:
                second_cai = species_cai.get(sl[1], 0.0)
                confidence = min(1.0, 0.5 + (chosen_cai - second_cai) * 5)
                confidence = max(0.0, confidence)
            else:
                confidence = 1.0
            # Compute CAI impact vs optimal
            cai_impact_val = 0.0
            opt = optimal_codon.get(aa, best)
            if opt != best:
                opt_cai = species_cai.get(opt, 0.0)
                cai_impact_val = chosen_cai - opt_cai  # negative = CAI lost
            _prov.record_codon_decision(CodonDecision(
                position=pos,
                amino_acid=aa,
                original_codon=None,
                chosen_codon=best,
                alternatives_considered=alternatives,
                constraint_reason="maximize_cai",
                confidence=round(confidence, 4),
                cai_impact=round(cai_impact_val, 6),
            ))

    # Build sequence as list for efficient slicing during swaps
    seq_chars = list("".join(codon_list))
    n_bases = len(seq_chars)
    n_codons = n_bases // 3

    # Compute initial CAI from incremental tracking
    phase1_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0

    # ── Helper: update CAI tracking for a single codon swap ──
    def _update_cai_for_swap(ci: int, old_codon: str, new_codon: str) -> None:
        nonlocal log_cai_sum, n_cai_codons
        # Remove old codon contribution
        aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
        if aa and aa != "*":
            old_adapt = species_cai.get(old_codon, 0.0)
            max_a = max_adapt.get(aa, 0.0)
            if max_a > 0 and old_adapt > 0:
                log_cai_sum -= math.log(old_adapt / max_a)
                n_cai_codons -= 1
            # Add new codon contribution
            new_adapt = species_cai.get(new_codon, 0.0)
            if max_a > 0 and new_adapt > 0:
                log_cai_sum += math.log(new_adapt / max_a)
                n_cai_codons += 1

    # ── Helper: check for restriction sites in local window ──
    def _has_rs_local(ci: int) -> bool:
        """Check if any restriction site overlaps with codon at ci.

        Uses Aho-Corasick scanner for O(window_len + matches) single-pass
        detection instead of iterating through each enzyme separately.
        Falls back to per-enzyme string search if no scanner is available.
        """
        if ac_scanner is not None:
            start = ci * 3
            check_start = max(0, start - max_rs_len + 1)
            check_end = min(n_bases, start + 3 + max_rs_len - 1)
            return ac_scanner.has_any_match_in_region(
                "".join(seq_chars), check_start, check_end
            )
        if not rs_sites:
            return False
        start = ci * 3
        check_start = max(0, start - max_rs_len + 1)
        check_end = min(n_bases, start + 3 + max_rs_len - 1)
        region = "".join(seq_chars[check_start:check_end])
        for site, site_rc in rs_sites:
            if site in region or (site_rc and site_rc in region):
                return True
        return False

    # ── Helper: check for restriction sites in full sequence ──
    def _has_rs_full() -> bool:
        """Check if any restriction site is present in the full sequence.

        Uses Aho-Corasick scanner for O(L+M) single-pass detection
        instead of iterating through each enzyme separately.
        """
        if ac_scanner is not None:
            return ac_scanner.has_any_match("".join(seq_chars))
        if not rs_sites:
            return False
        seq_str = "".join(seq_chars)
        for site, site_rc in rs_sites:
            if site in seq_str or (site_rc and site_rc in seq_str):
                return True
        return False

    # ── Helper: apply a codon swap ──
    def _apply_swap(ci: int, new_codon: str) -> str:
        nonlocal gc_count
        start = ci * 3
        old_codon = "".join(seq_chars[start:start + 3])
        old_gc = codon_gc.get(old_codon, 0)
        new_gc = codon_gc.get(new_codon, 0)
        gc_count += (new_gc - old_gc)
        seq_chars[start] = new_codon[0]
        seq_chars[start + 1] = new_codon[1]
        seq_chars[start + 2] = new_codon[2]
        _update_cai_for_swap(ci, old_codon, new_codon)
        return old_codon

    # ── Helper: rollback a codon swap ──
    def _rollback_swap(ci: int, old_codon: str) -> None:
        nonlocal gc_count
        start = ci * 3
        current_codon = "".join(seq_chars[start:start + 3])
        current_gc = codon_gc.get(current_codon, 0)
        old_gc = codon_gc.get(old_codon, 0)
        gc_count += (old_gc - current_gc)
        seq_chars[start] = old_codon[0]
        seq_chars[start + 1] = old_codon[1]
        seq_chars[start + 2] = old_codon[2]
        # Restore CAI tracking
        _update_cai_for_swap(ci, current_codon, old_codon)

    violations_fixed = 0
    warnings: list[str] = []
    rs_fixed_count = 0
    gc_fixed_count = 0
    attta_fixed_count = 0
    trun_fixed_count = 0

    # ── Phase 2a: Fix restriction sites (single scan) ──
    # Use Aho-Corasick scanner for O(L+M) multi-pattern detection
    # instead of per-enzyme O(N*L*site_len) scanning.
    if rs_sites or ac_scanner is not None:
        seq_str = "".join(seq_chars)
        for _iter in range(100):  # Max iterations total
            # Single-pass scan for ALL restriction sites simultaneously
            if ac_scanner is not None:
                matches = ac_scanner.scan(seq_str)
                if not matches:
                    break
                pos, site_match, _enzyme = matches[0]
                site_len = len(site_match)
            else:
                # Fallback: per-enzyme scan
                found = False
                for site, site_rc in rs_sites:
                    p = seq_str.find(site)
                    if p == -1 and site_rc:
                        p = seq_str.find(site_rc)
                    if p != -1:
                        pos = p
                        site_len = len(site)
                        site_match = site
                        found = True
                        break
                if not found:
                    break

            # Get overlapping codon indices
            first_ci = pos // 3
            last_ci = (pos + site_len - 1) // 3
            fixed = False

            # CAI-aware: collect ALL possible single-codon swaps across
            # overlapping positions, then pick the one with highest CAI.
            best_swap_info: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
            for ci in range(max(0, first_ci), min(n_codons, last_ci + 1)):
                aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    if not _has_rs_local(ci):
                        alt_cai = species_cai.get(alt, 0.0)
                        log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                        if best_swap_info is None or log_cai > best_swap_info[2]:
                            best_swap_info = (ci, alt, log_cai)
                        _rollback_swap(ci, old_codon)
                    else:
                        _rollback_swap(ci, old_codon)

            if best_swap_info is not None:
                ci, alt, _ = best_swap_info
                _apply_swap(ci, alt)
                seq_str = "".join(seq_chars)
                violations_fixed += 1
                rs_fixed_count += 1
                fixed = True

            if not fixed:
                # Try two-codon coordinated fix with CAI scoring
                # Evaluate ALL pair combinations and pick the one with
                # the highest combined CAI (minimises CAI loss).
                best_pair_info: tuple[int, str, str, float] | None = None  # (ci1, c1, c2, log_cai)
                alts1_full = sorted_codons.get(aa1, []) if aa1 else []
                alts2_full = sorted_codons.get(aa2, []) if aa2 else []

                for idx in range(max(0, first_ci), min(n_codons - 1, last_ci + 1)):
                    ci1, ci2 = idx, idx + 1
                    aa1 = protein[ci1] if ci1 < len(protein) else None
                    aa2 = protein[ci2] if ci2 < len(protein) else None
                    if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                        continue
                    current1 = "".join(seq_chars[ci1*3:ci1*3+3])
                    current2 = "".join(seq_chars[ci2*3:ci2*3+3])
                    alts1 = sorted_codons.get(aa1, [])
                    alts2 = sorted_codons.get(aa2, [])

                    for c1 in alts1:
                        for c2 in alts2:
                            if c1 == current1 and c2 == current2:
                                continue
                            old1 = _apply_swap(ci1, c1)
                            old2 = _apply_swap(ci2, c2)
                            if not _has_rs_local(ci1) and not _has_rs_local(ci2):
                                c1_cai = species_cai.get(c1, 0.0)
                                c2_cai = species_cai.get(c2, 0.0)
                                combined_log = (
                                    (math.log(c1_cai) if c1_cai > 0 else float('-inf'))
                                    + (math.log(c2_cai) if c2_cai > 0 else float('-inf'))
                                )
                                if best_pair_info is None or combined_log > best_pair_info[3]:
                                    best_pair_info = (ci1, c1, c2, combined_log)
                                _rollback_swap(ci2, old2)
                                _rollback_swap(ci1, old1)
                            else:
                                _rollback_swap(ci2, old2)
                                _rollback_swap(ci1, old1)

                if best_pair_info is not None:
                    ci1, c1, c2, _ = best_pair_info
                    _apply_swap(ci1, c1)
                    _apply_swap(ci2, c2)
                    seq_str = "".join(seq_chars)
                    violations_fixed += 1
                    rs_fixed_count += 1
                    fixed = True

                if not fixed:
                    warnings.append(
                        f"Cannot remove restriction site {site_match} at pos {pos}"
                    )
                    break

    # ── Phase 2b: Fix GC content (incremental) ──
    gc_val = gc_count / n_bases if n_bases > 0 else 0.0
    if not (gc_lo <= gc_val <= gc_hi):
        target = gc_lo if gc_val < gc_lo else gc_hi
        need_more_gc = gc_val < gc_lo

        for _iter in range(200):
            if gc_lo <= gc_val <= gc_hi:
                break

            best_swap = None
            best_ci = -1
            best_score = -1.0
            best_gc_delta = 0

            for ci in range(n_codons):
                aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])
                current_gc_val = codon_gc.get(current, 0)

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    alt_gc = codon_gc.get(alt, 0)
                    gc_delta = alt_gc - current_gc_val

                    if need_more_gc and gc_delta <= 0:
                        continue
                    if not need_more_gc and gc_delta >= 0:
                        continue

                    new_gc_count = gc_count + gc_delta
                    new_frac = new_gc_count / n_bases
                    diff = abs(new_frac - target)
                    alt_cai = species_cai.get(alt, 0.0)
                    score = (1.0 - diff) + alt_cai * 0.01

                    if score > best_score:
                        best_score = score
                        best_swap = alt
                        best_ci = ci
                        best_gc_delta = gc_delta

            if best_swap is None:
                break

            # Apply the swap and check RS
            old_codon = _apply_swap(best_ci, best_swap)
            if _has_rs_local(best_ci):
                _rollback_swap(best_ci, old_codon)
                break  # Can't fix without breaking RS

            gc_val = gc_count / n_bases
            violations_fixed += 1
            gc_fixed_count += 1

    # ── Phase 2c: Fix ATTTA motifs ──
    seq_str = "".join(seq_chars)
    for _iter in range(100):
        pos = seq_str.find("ATTTA")
        if pos == -1:
            break

        first_ci = max(0, (pos // 3) - 1)
        last_ci = min(n_codons, ((pos + 4) // 3) + 2)
        fixed = False

        for ci in range(first_ci, last_ci):
            aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
            if aa is None or aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])

            for alt in sorted_codons.get(aa, []):
                if alt == current:
                    continue
                old_codon = _apply_swap(ci, alt)
                # Check ATTTA eliminated and no new RS
                new_local = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                if "ATTTA" not in new_local and not _has_rs_local(ci):
                    seq_str = "".join(seq_chars)
                    violations_fixed += 1
                    attta_fixed_count += 1
                    fixed = True
                    break
                else:
                    _rollback_swap(ci, old_codon)

            if fixed:
                break

        if not fixed:
            warnings.append(f"Cannot remove ATTTA motif at pos {pos}")
            break

    # ── Phase 2d: Fix T-runs (6+ consecutive T) ──
    for _iter in range(100):
        # Find longest T-run
        max_run = 0
        max_pos = -1
        i = 0
        while i < n_bases:
            if seq_chars[i] == 'T':
                j = i
                while j < n_bases and seq_chars[j] == 'T':
                    j += 1
                if j - i > max_run:
                    max_run = j - i
                    max_pos = i
                i = j
            else:
                i += 1

        if max_run < 6:
            break

        # Fix at the center of the run
        ci = (max_pos + max_run // 2) // 3
        if ci >= n_codons:
            ci = n_codons - 1
        fixed = False

        aa = protein[ci] if ci < len(protein) else None  # type: ignore[assignment]
        if aa is not None and aa != "*":
            current = "".join(seq_chars[ci*3:ci*3+3])

            for alt in sorted_codons.get(aa, []):
                if alt == current:
                    continue
                old_codon = _apply_swap(ci, alt)
                # Quick T-run check in local window
                check_start = max(0, ci * 3 - 6)
                check_end = min(n_bases, ci * 3 + 9)
                has_long_run = False
                j = check_start
                while j < check_end:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_end and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            has_long_run = True
                            break
                        j = k
                    else:
                        j += 1

                if not has_long_run and not _has_rs_local(ci):
                    violations_fixed += 1
                    trun_fixed_count += 1
                    fixed = True
                    break
                else:
                    _rollback_swap(ci, old_codon)

        if not fixed:
            warnings.append(
                f"Cannot fix T-run of {max_run} at pos {max_pos}"
            )
            break

    # ── Compute final metrics ──
    seq = "".join(seq_chars)
    final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
    gc = gc_count / n_bases if n_bases > 0 else 0.0

    # ── Phase 2g: CAI hill climbing for prokaryotes ──
    # After fixing all constraints, try to upgrade any suboptimal codons
    # back to higher-CAI alternatives without breaking constraints.
    hill_climb_improvements = 0
    for _hc_iter in range(3):  # Limited passes to avoid excessive runtime
        improved = False
        for ci in range(n_codons):
            aa = protein[ci] if ci < len(protein) else None
            if aa is None or aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])
            current_cai = species_cai.get(current, 0.0)

            # Find the best alternative that's better than current
            for alt in sorted_codons.get(aa, []):
                if alt == current:
                    continue
                alt_cai = species_cai.get(alt, 0.0)
                if alt_cai <= current_cai:
                    continue
                # Check that this swap doesn't break constraints
                old_codon = _apply_swap(ci, alt)
                # Check no RS introduced locally
                if _has_rs_local(ci):
                    _rollback_swap(ci, old_codon)
                    continue
                # Check GC still in range
                new_gc_val = gc_count / n_bases
                if not (gc_lo <= new_gc_val <= gc_hi):
                    _rollback_swap(ci, old_codon)
                    continue
                # Check no ATTTA introduced locally
                local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                if "ATTTA" in local_window:
                    _rollback_swap(ci, old_codon)
                    continue
                # Check no long T-run introduced
                check_start = max(0, ci * 3 - 6)
                check_end = min(n_bases, ci * 3 + 9)
                has_long_run = False
                j = check_start
                while j < check_end:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_end and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            has_long_run = True
                            break
                        j = k
                    else:
                        j += 1
                if has_long_run:
                    _rollback_swap(ci, old_codon)
                    continue
                # This swap is safe and improves CAI
                hill_climb_improvements += 1
                improved = True
                break  # First valid alt in sorted order is the best

        if not improved:
            break

    # ── Phase 2h: CAI micro-optimization with constraint recheck ──
    # After all constraint phases, some codons may be suboptimal because
    # a restriction site existed earlier but was later broken by other
    # codon changes.  Re-check each non-optimal codon to see if the
    # optimal (or a higher-CAI) codon is now safe.  Also try coordinated
    # two-codon swaps: if upgrading codon A creates a RS, maybe changing
    # a neighboring codon B simultaneously resolves it with a net CAI gain.
    micro_opt_improvements = 0
    for _micro_iter in range(5):
        any_micro_improved = False

        for ci in range(n_codons):
            aa = protein[ci] if ci < len(protein) else None
            if aa is None or aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])
            current_cai = species_cai.get(current, 0.0)
            optimal = optimal_codon.get(aa, current)
            optimal_cai = species_cai.get(optimal, 0.0)

            # Already optimal? Skip
            if current == optimal or optimal_cai <= current_cai:
                continue

            # Try upgrading to the optimal codon
            old_codon = _apply_swap(ci, optimal)
            if not _has_rs_local(ci):
                # Also verify other constraints
                new_gc_val = gc_count / n_bases
                local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                has_attta = "ATTTA" in local_window
                # T-run check
                check_s = max(0, ci * 3 - 6)
                check_e = min(n_bases, ci * 3 + 9)
                has_long_run = False
                j = check_s
                while j < check_e:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_e and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            has_long_run = True
                            break
                        j = k
                    else:
                        j += 1
                if (gc_lo <= new_gc_val <= gc_hi) and not has_attta and not has_long_run:
                    micro_opt_improvements += 1
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    continue
            _rollback_swap(ci, old_codon)

            # Optimal codon creates a RS or violates another constraint.
            # Try coordinated two-codon swaps with immediate neighbors.
            # For each neighbor, try upgrading ci to optimal while
            # simultaneously changing the neighbor to a CAI-friendly
            # alternative that breaks the RS.
            for neighbor_ci in [ci - 1, ci + 1]:
                if neighbor_ci < 0 or neighbor_ci >= n_codons:
                    continue
                neighbor_aa = protein[neighbor_ci] if neighbor_ci < len(protein) else None
                if neighbor_aa is None or neighbor_aa == "*":
                    continue
                neighbor_current = "".join(seq_chars[neighbor_ci*3:neighbor_ci*3+3])
                neighbor_current_cai = species_cai.get(neighbor_current, 0.0)

                # Try each alternative for the neighbor
                best_net_log = 0.0  # net log-CAI gain must be positive
                best_neighbor_alt = None

                for neighbor_alt in sorted_codons.get(neighbor_aa, []):
                    if neighbor_alt == neighbor_current:
                        continue
                    # Apply neighbor swap first, then try optimal at ci
                    old_neighbor = _apply_swap(neighbor_ci, neighbor_alt)
                    old_ci = _apply_swap(ci, optimal)

                    if not _has_rs_local(ci) and not _has_rs_local(neighbor_ci):
                        # Check other constraints
                        new_gc_val = gc_count / n_bases
                        local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                        has_attta = "ATTTA" in local_window
                        check_s = max(0, ci * 3 - 6)
                        check_e = min(n_bases, ci * 3 + 9)
                        has_long_run = False
                        j = check_s
                        while j < check_e:
                            if seq_chars[j] == 'T':
                                k = j
                                while k < check_e and seq_chars[k] == 'T':
                                    k += 1
                                if k - j >= 6:
                                    has_long_run = True
                                    break
                                j = k
                            else:
                                j += 1
                        if (gc_lo <= new_gc_val <= gc_hi) and not has_attta and not has_long_run:
                            # Compute net CAI gain
                            neighbor_alt_cai = species_cai.get(neighbor_alt, 0.0)
                            optimal_cai_val = species_cai.get(optimal, 0.0)
                            net_log_gain = 0.0
                            if optimal_cai_val > 0 and current_cai > 0:
                                net_log_gain += math.log(optimal_cai_val) - math.log(current_cai)
                            if neighbor_alt_cai > 0 and neighbor_current_cai > 0:
                                net_log_gain += math.log(neighbor_alt_cai) - math.log(neighbor_current_cai)
                            if net_log_gain > best_net_log:
                                best_net_log = net_log_gain
                                best_neighbor_alt = neighbor_alt

                    _rollback_swap(ci, old_ci)
                    _rollback_swap(neighbor_ci, old_neighbor)

                if best_neighbor_alt is not None and best_net_log > 0:
                    # Apply the best coordinated swap
                    _apply_swap(neighbor_ci, best_neighbor_alt)
                    _apply_swap(ci, optimal)
                    micro_opt_improvements += 1
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    break  # Move to next codon

            # Also try non-optimal but higher-CAI alternatives for this codon
            # (without coordinated swaps — simple single upgrade)
            current_after = "".join(seq_chars[ci*3:ci*3+3])
            current_after_cai = species_cai.get(current_after, 0.0)
            for alt in sorted_codons.get(aa, []):
                if alt == current_after:
                    continue
                alt_cai = species_cai.get(alt, 0.0)
                if alt_cai <= current_after_cai:
                    continue
                old_codon2 = _apply_swap(ci, alt)
                if _has_rs_local(ci):
                    _rollback_swap(ci, old_codon2)
                    continue
                new_gc_val = gc_count / n_bases
                if not (gc_lo <= new_gc_val <= gc_hi):
                    _rollback_swap(ci, old_codon2)
                    continue
                local_window = "".join(seq_chars[max(0, ci*3-5):min(n_bases, ci*3+8)])
                if "ATTTA" in local_window:
                    _rollback_swap(ci, old_codon2)
                    continue
                check_s = max(0, ci * 3 - 6)
                check_e = min(n_bases, ci * 3 + 9)
                has_long_run = False
                j = check_s
                while j < check_e:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_e and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            has_long_run = True
                            break
                        j = k
                    else:
                        j += 1
                if has_long_run:
                    _rollback_swap(ci, old_codon2)
                    continue
                hill_climb_improvements += 1
                any_micro_improved = True
                break

        if not any_micro_improved:
            break

    # ── Phase 2i: CAI micro-optimization — recheck suboptimal codons ──
    # After all constraint phases (including Phase 2h coordinated swaps),
    # some codons that were suboptimally chosen (e.g., to avoid a
    # restriction site) may now be safe to upgrade because nearby codons
    # changed during later phases.  This final pass systematically tries
    # ALL higher-CAI alternatives, not just the optimal codon.
    for _micro_iter in range(3):
        any_micro_upgrade = False
        for ci in range(n_codons):
            aa = protein[ci] if ci < len(protein) else None
            if aa is None or aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])
            current_cai = species_cai.get(current, 0.0)

            # Try all higher-CAI alternatives (not just optimal)
            for alt in sorted_codons.get(aa, []):
                alt_cai = species_cai.get(alt, 0.0)
                if alt_cai <= current_cai:
                    break  # sorted_codons is CAI-descending, no point continuing
                if alt == current:
                    continue

                old_codon = _apply_swap(ci, alt)

                # Check ALL constraints (no GT/AG for prokaryotes)
                rs_ok = not _has_rs_local(ci)
                gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi

                # Check local region for ATTTA
                local_start = max(0, ci * 3 - 5)
                local_end = min(n_bases, ci * 3 + 8)
                local_str = "".join(seq_chars[local_start:local_end])
                attta_ok = "ATTTA" not in local_str

                # Check T-runs in local window
                check_s = max(0, ci * 3 - 6)
                check_e = min(n_bases, ci * 3 + 9)
                trun_ok = True
                j = check_s
                while j < check_e:
                    if seq_chars[j] == 'T':
                        k = j
                        while k < check_e and seq_chars[k] == 'T':
                            k += 1
                        if k - j >= 6:
                            trun_ok = False
                            break
                        j = k
                    else:
                        j += 1

                if rs_ok and gc_ok and attta_ok and trun_ok:
                    any_micro_upgrade = True
                    break  # Keep this upgrade, move to next codon
                else:
                    _rollback_swap(ci, old_codon)

        if not any_micro_upgrade:
            break

    # Recompute final CAI after hill climbing and micro-optimization
    final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
    gc = gc_count / n_bases if n_bases > 0 else 0.0
    seq = "".join(seq_chars)

    # ── Phase 4: Codon pair bias optimization (prokaryote fast path) ──
    cpb_improvements = 0
    mean_cpb = 0.0
    if optimizer.consider_codon_pair_bias:
        seq, final_cai, cpb_improvements, mean_cpb = (
            _codon_pair_bias_optimize_prokaryote(optimizer, 
                seq, protein, final_cai
            )
        )
        # Recompute gc from the sequence (may have changed slightly)
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)

    # ── Record constraint decisions in provenance ──
    # This records the _decision_trail for the prokaryote fast path,
    # capturing which constraints were active and what actions were taken.
    if _prov is not None:
        if rs_fixed_count > 0:
            _prov.record_constraint_decision(ConstraintDecision(
                constraint_name="NoRestrictionSite",
                constraint_type="hard",
                action_taken="satisfied",
                positions_affected=[],
                tradeoff_description=(
                    f"Fixed {rs_fixed_count} restriction site(s) "
                    f"by swapping to CAI-optimal alternative codons"
                ),
                impact_on_cai=0.0,  # CAI-aware fixes minimize CAI loss
            ))
        if gc_fixed_count > 0:
            _prov.record_constraint_decision(ConstraintDecision(
                constraint_name="GCInRange",
                constraint_type="hard",
                action_taken="satisfied",
                positions_affected=[],
                tradeoff_description=(
                    f"Fixed GC content ({gc_fixed_count} codon swap(s)) "
                    f"to bring GC fraction into [{gc_lo:.2f}, {gc_hi:.2f}]"
                ),
                impact_on_cai=0.0,
            ))
        if attta_fixed_count > 0:
            _prov.record_constraint_decision(ConstraintDecision(
                constraint_name="NoInstabilityMotif",
                constraint_type="soft",
                action_taken="satisfied",
                positions_affected=[],
                tradeoff_description=(
                    f"Removed {attta_fixed_count} ATTTA instability motif(s)"
                ),
                impact_on_cai=0.0,
            ))
        if trun_fixed_count > 0:
            _prov.record_constraint_decision(ConstraintDecision(
                constraint_name="NoLongTRun",
                constraint_type="soft",
                action_taken="satisfied",
                positions_affected=[],
                tradeoff_description=(
                    f"Fixed {trun_fixed_count} T-run(s) of 6+ consecutive T"
                ),
                impact_on_cai=0.0,
            ))
        if hill_climb_improvements > 0:
            _prov.record_iteration({
                "step": "cai_hill_climb",
                "action": "upgrade_codons",
                "improvements": hill_climb_improvements,
                "score": round(final_cai, 6),
            })
        _prov.record_iteration({
            "step": "prokaryote_fast_path_complete",
            "action": "finalize",
            "violations_fixed": violations_fixed,
            "hill_climb_improvements": hill_climb_improvements,
            "cai": round(final_cai, 6),
            "gc": round(gc, 4),
        })

    return HybridResult(
        sequence=seq,
        cai=final_cai,
        gc_content=round(gc, 4),
        violations_fixed=violations_fixed,
        hill_climb_improvements=hill_climb_improvements,
        iterations_used=0,
        phase1_cai=phase1_cai,
        phase2_cai=final_cai,
        phase3_cai=final_cai,
        phase4_cai=final_cai,
        cpb_improvements=cpb_improvements,
        mean_cpb=mean_cpb,
        warnings=warnings,
    )


def _codon_pair_bias_optimize_prokaryote(
    optimizer,
    seq: str,
    protein: str,
    current_cai: float,
    cai_threshold_fraction: float = 0.95,
    max_iterations: int = 3,
) -> tuple[str, float, int, float]:
    """Phase 4 for prokaryote fast path: codon pair bias optimization.

    Simplified version for prokaryotes — no GT/AG/splice checks needed.
    Only verifies restriction sites, GC range, ATTTA, and T-runs.

    Args:
        seq: Current optimized DNA sequence.
        protein: Amino acid sequence.
        current_cai: CAI after prokaryote fast path.
        cai_threshold_fraction: Minimum fraction of current_CAI to
            maintain (default 0.95).
        max_iterations: Maximum number of full-sequence passes.

    Returns:
        (sequence, cai, cpb_improvements, mean_cpb) tuple.
    """
    cai_floor = current_cai * cai_threshold_fraction
    cpb_improvements = 0
    aas = list(protein)
    n_codons = len(aas)

    for _iteration in range(max_iterations):
        any_improved = False

        for ci in range(n_codons - 1):
            aa1 = aas[ci]
            aa2 = aas[ci + 1]

            # Skip stop codons
            if aa1 == "*" or aa2 == "*":
                continue

            current_c1 = seq[ci * 3:ci * 3 + 3]
            current_c2 = seq[(ci + 1) * 3:(ci + 1) * 3 + 3]

            current_cpb = optimizer._score_cpb_pair(current_c1, current_c2)

            # Find the best synonymous pair by CPB score
            best_pair: tuple[str, str] | None = None
            best_cpb = current_cpb
            current_cai_sum = (
                optimizer.species_cai.get(current_c1, 0.0)
                + optimizer.species_cai.get(current_c2, 0.0)
            )

            synonyms1 = AA_TO_CODONS.get(aa1, [current_c1])
            synonyms2 = AA_TO_CODONS.get(aa2, [current_c2])

            for alt1 in synonyms1:
                for alt2 in synonyms2:
                    if alt1 == current_c1 and alt2 == current_c2:
                        continue

                    pair_cpb = optimizer._score_cpb_pair(alt1, alt2)

                    if pair_cpb <= best_cpb:
                        continue

                    alt_cai_sum = (
                        optimizer.species_cai.get(alt1, 0.0)
                        + optimizer.species_cai.get(alt2, 0.0)
                    )

                    if alt_cai_sum < current_cai_sum:
                        cai_drop = current_cai_sum - alt_cai_sum
                        max_allowed_drop = (1.0 - cai_threshold_fraction) * 2.0
                        if cai_drop > max_allowed_drop:
                            continue

                    if pair_cpb > best_cpb:
                        best_cpb = pair_cpb
                        best_pair = (alt1, alt2)

            if best_pair is None:
                continue

            alt1, alt2 = best_pair

            # Apply the swap
            test_seq = (
                seq[:ci * 3]
                + alt1
                + seq[ci * 3 + 3:(ci + 1) * 3]
                + alt2
                + seq[(ci + 1) * 3 + 3:]
            )

            # Verify CAI floor
            test_cai = optimizer._compute_cai(test_seq)
            if test_cai < cai_floor:
                continue

            # Prokaryote constraint checks (no GT/AG/splice)
            valid = True

            # 1. No new restriction sites
            if optimizer._rs_sites:
                for site, site_rc in optimizer._rs_sites:
                    if site in test_seq or (site_rc and site_rc in test_seq):
                        valid = False
                        break

            # 2. GC still in range
            if valid:
                gc = (test_seq.count("G") + test_seq.count("C")) / max(len(test_seq), 1)
                if not (optimizer.gc_lo <= gc <= optimizer.gc_hi):
                    valid = False

            # 3. No new ATTTA motifs
            if valid:
                if test_seq.count("ATTTA") > seq.count("ATTTA"):
                    valid = False

            # 4. No new premature stop codons at the swapped positions
            if valid:
                if alt1 in ("TAA", "TAG", "TGA") or alt2 in ("TAA", "TAG", "TGA"):
                    valid = False

            if valid:
                seq = test_seq
                cpb_improvements += 1
                any_improved = True

        if not any_improved:
            break

    final_cai = optimizer._compute_cai(seq)
    mean_cpb = optimizer._compute_mean_cpb(seq)
    return seq, final_cai, cpb_improvements, mean_cpb
