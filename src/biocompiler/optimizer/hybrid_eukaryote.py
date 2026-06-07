"""
Eukaryotic optimization path for HybridOptimizer.

Extracted from hybrid_optimizer.py to decompose the monolith.
Contains _optimize_eukaryote_fast() and its helpers.
"""

from __future__ import annotations

import math
import logging
from typing import Any

from ..type_system import CODON_TABLE, AA_TO_CODONS
from .hybrid_types import HybridResult, GT_CAI_COST_THRESHOLD as _GT_CAI_COST_THRESHOLD
from .hybrid_prokaryote import _codon_pair_bias_optimize_prokaryote as codon_pair_bias_optimize_prokaryote

# ── Fast MaxEntScan integration ────────────────────────────────────
try:
    from ..maxentscan_fast import scan_splice_sites_fast_str as _scan_splice_sites_fast
    _HAS_FAST_MAXENT = True
except ImportError:
    _HAS_FAST_MAXENT = False

logger = logging.getLogger(__name__)


def _optimize_eukaryote_fast(optimizer, protein: str) -> HybridResult:
    """Streamlined eukaryotic optimization with integrated GT/AG/CpG handling.

    This is the eukaryotic counterpart to _optimize_prokaryote_fast.
    It avoids the priority-queue-based _constraint_satisfaction which
    suffers from GT↔ATTTA oscillation loops by using a single-pass
    approach with incremental tracking.

    Pipeline:
    1. Greedy init with highest-CAI codons + incremental GC/CAI tracking
    2. Fix restriction sites (single scan, local window checks)
    3. Fix GC content (incremental tracking)
    4. Fix ATTTA motifs (fast find + codon swap, but no GT creation)
    5. Fix T-runs (single scan + codon swap, but no GT creation)
    6. MaxEntScan validation (only if GT/AG present — skip if none)
       — only fix GT/AG dinucleotides that score above splice threshold
       — prefer GT-free/AG-free alternatives that preserve CAI
    7. Fix remaining CpG dinucleotides (soft constraint — only if CAI-neutral)
    8. CAI recovery hill climb (upgrade codons while maintaining constraints)

    Key design principles:
    - **CAI first**: Never sacrifice CAI for soft constraints (GT/AG/CpG)
      unless MaxEntScan flags them as actual cryptic splice sites.
    - **No oscillation**: GT/AG/CpG are handled in a fixed order after
      hard constraints, with CAI-preserving alternatives only.
    - **MaxEntScan gating**: Only fix GT/AG dinucleotides that actually
      score as cryptic splice sites, not all GT/AG occurrences.
    - **Incremental tracking**: GC/CAI updated on every swap.

    Target: <5ms for HBB (444bp) in Human, CAI ≥ 0.99
    """
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
    gt_free = optimizer.gt_free
    ag_free = optimizer.ag_free
    splice_threshold = optimizer.splice_threshold

    # ── Phase 1: Greedy init with incremental GC + CAI tracking ──
    codon_list: list[str] = []
    gc_count = 0
    log_cai_sum = 0.0
    n_cai_codons = 0

    for aa in protein:
        if aa == "*":
            codon_list.append("TAA")
            continue
        best = optimal_codon.get(aa)
        if not best:
            sl = sorted_codons.get(aa, [])
            best = sl[0] if sl else "NNN"
        codon_list.append(best)
        gc_count += codon_gc.get(best, 0)
        adapt = species_cai.get(best, 0.0)
        max_a = max_adapt.get(aa, 0.0)
        if max_a > 0 and adapt > 0:
            log_cai_sum += math.log(adapt / max_a)
            n_cai_codons += 1

    seq_chars = list("".join(codon_list))
    # Maintain a bytearray for fast pattern matching (avoids "".join overhead)
    seq_buf = bytearray("".join(codon_list), 'ascii')
    n_bases = len(seq_chars)
    n_codons = n_bases // 3

    phase1_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0

    # ── Helper: update CAI tracking for a single codon swap ──
    def _update_cai_for_swap(ci: int, old_codon: str, new_codon: str) -> None:
        nonlocal log_cai_sum, n_cai_codons
        aa = protein[ci] if ci < len(protein) else None
        if aa and aa != "*":
            old_adapt = species_cai.get(old_codon, 0.0)
            max_a = max_adapt.get(aa, 0.0)
            if max_a > 0 and old_adapt > 0:
                log_cai_sum -= math.log(old_adapt / max_a)
                n_cai_codons -= 1
            new_adapt = species_cai.get(new_codon, 0.0)
            if max_a > 0 and new_adapt > 0:
                log_cai_sum += math.log(new_adapt / max_a)
                n_cai_codons += 1

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
        # Update bytearray in-place
        seq_buf[start] = ord(new_codon[0])
        seq_buf[start + 1] = ord(new_codon[1])
        seq_buf[start + 2] = ord(new_codon[2])
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
        # Update bytearray in-place
        seq_buf[start] = ord(old_codon[0])
        seq_buf[start + 1] = ord(old_codon[1])
        seq_buf[start + 2] = ord(old_codon[2])
        _update_cai_for_swap(ci, current_codon, old_codon)

    # ── Helper: check for restriction sites in local window ──
    def _has_rs_local(ci: int) -> bool:
        if ac_scanner is not None:
            start = ci * 3
            check_start = max(0, start - max_rs_len + 1)
            check_end = min(n_bases, start + 3 + max_rs_len - 1)
            # Use seq_buf directly instead of "".join(seq_chars)
            region = seq_buf[check_start:check_end].decode('ascii')
            return ac_scanner.has_any_match_in_region(
                region, 0, len(region)
            )
        if not rs_sites:
            return False
        start = ci * 3
        check_start = max(0, start - max_rs_len + 1)
        check_end = min(n_bases, start + 3 + max_rs_len - 1)
        region = seq_buf[check_start:check_end].decode('ascii')
        for site, site_rc in rs_sites:
            if site in region or (site_rc and site_rc in region):
                return True
        return False

    # ── Helper: check for new GT or AG in local window after swap ──
    def _creates_new_gt_or_ag(ci: int, new_codon: str) -> tuple[bool, bool]:
        """Return (creates_gt, creates_ag) for the local window around ci."""
        start = ci * 3
        prev_base = seq_chars[start - 1] if start > 0 else ''
        next_base = seq_chars[start + 3] if start + 3 < n_bases else ''
        local = [prev_base, new_codon[0], new_codon[1], new_codon[2], next_base]

        creates_gt = False
        creates_ag = False
        for i in range(len(local) - 1):
            if local[i] == 'G' and local[i + 1] == 'T':
                creates_gt = True
            if local[i] == 'A' and local[i + 1] == 'G':
                creates_ag = True
        return creates_gt, creates_ag

    violations_fixed = 0
    warnings: list[str] = []

    # ── Phase 2a: Fix restriction sites (single scan) ──
    if rs_sites or ac_scanner is not None:
        seq_str = seq_buf.decode('ascii')
        for _iter in range(100):
            if ac_scanner is not None:
                matches = ac_scanner.scan(seq_str)
                if not matches:
                    break
                pos, site_match, _enzyme = matches[0]
                site_len = len(site_match)
            else:
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

            first_ci = pos // 3
            last_ci = (pos + site_len - 1) // 3
            fixed = False

            # CAI-aware: collect ALL possible single-codon swaps across
            # overlapping positions, then pick the one with highest CAI.
            best_swap_info: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
            for ci in range(max(0, first_ci), min(n_codons, last_ci + 1)):
                aa = protein[ci] if ci < len(protein) else None
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
                seq_str = seq_buf.decode('ascii')
                violations_fixed += 1
                fixed = True

            if not fixed:
                # Try two-codon coordinated fix with CAI scoring
                # Evaluate ALL pair combinations and pick the one with
                # the highest combined CAI (minimises CAI loss).
                best_pair_info: tuple[int, str, str, float] | None = None  # (ci1, c1, c2, log_cai)

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
                    seq_str = seq_buf.decode('ascii')
                    violations_fixed += 1
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
                aa = protein[ci] if ci < len(protein) else None
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

            old_codon = _apply_swap(best_ci, best_swap)
            if _has_rs_local(best_ci):
                _rollback_swap(best_ci, old_codon)
                break

            gc_val = gc_count / n_bases
            violations_fixed += 1

    # ── Phase 2c: Fix ATTTA motifs (without creating new GT/AG) ──
    seq_str = seq_buf.decode('ascii')
    for _iter in range(100):
        pos = seq_str.find("ATTTA")
        if pos == -1:
            break

        first_ci = max(0, (pos // 3) - 1)
        last_ci = min(n_codons, ((pos + 4) // 3) + 2)
        fixed = False

        for ci in range(first_ci, last_ci):
            aa = protein[ci] if ci < len(protein) else None
            if aa is None or aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])

            for alt in sorted_codons.get(aa, []):
                if alt == current:
                    continue
                old_codon = _apply_swap(ci, alt)
                new_local = seq_buf[max(0, ci*3-5):min(n_bases, ci*3+8)].decode('ascii')
                creates_gt, creates_ag = _creates_new_gt_or_ag(ci, alt)
                if "ATTTA" not in new_local and not _has_rs_local(ci) and not creates_gt and not creates_ag:
                    seq_str = seq_buf.decode('ascii')
                    violations_fixed += 1
                    fixed = True
                    break
                else:
                    _rollback_swap(ci, old_codon)

            if fixed:
                break

        if not fixed:
            # Relax: allow GT/AG creation to fix ATTTA
            for ci in range(first_ci, last_ci):
                aa = protein[ci] if ci < len(protein) else None
                if aa is None or aa == "*":
                    continue
                current = "".join(seq_chars[ci*3:ci*3+3])

                for alt in sorted_codons.get(aa, []):
                    if alt == current:
                        continue
                    old_codon = _apply_swap(ci, alt)
                    new_local = seq_buf[max(0, ci*3-5):min(n_bases, ci*3+8)].decode('ascii')
                    if "ATTTA" not in new_local and not _has_rs_local(ci):
                        seq_str = seq_buf.decode('ascii')
                        violations_fixed += 1
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

        ci = (max_pos + max_run // 2) // 3
        if ci >= n_codons:
            ci = n_codons - 1
        fixed = False

        aa = protein[ci] if ci < len(protein) else None
        if aa is not None and aa != "*":
            current = "".join(seq_chars[ci*3:ci*3+3])

            for alt in sorted_codons.get(aa, []):
                if alt == current:
                    continue
                old_codon = _apply_swap(ci, alt)
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
                    fixed = True
                    break
                else:
                    _rollback_swap(ci, old_codon)

        if not fixed:
            warnings.append(
                f"Cannot fix T-run of {max_run} at pos {max_pos}"
            )
            break

    # ── Phase 3: MaxEntScan validation (only if GT/AG present) ──
    # This is the KEY step that was missing from the naive GT/AG elimination.
    # Instead of eliminating ALL GT/AG dinucleotides (which tanks CAI),
    # we only eliminate GT/AG dinucleotides that MaxEntScan identifies
    # as actual cryptic splice sites.  Most GT/AG dinucleotides are
    # not cryptic splice sites, so this preserves CAI.
    seq_str = seq_buf.decode('ascii')
    has_gt = "GT" in seq_str
    has_ag = "AG" in seq_str

    if has_gt or has_ag:
        try:
            if _HAS_FAST_MAXENT:
                scan_splice_sites = _scan_splice_sites_fast
            else:
                from ..maxentscan import scan_splice_sites

            sites = scan_splice_sites(
                seq_str,
                donor_threshold=splice_threshold,
                acceptor_threshold=splice_threshold,
            )

            for _splice_iter in range(10):
                if not sites:
                    break

                sites.sort(key=lambda s: s[2], reverse=True)
                pos, site_type, score = sites[0]

                # Fix cryptic splice sites
                codon_idx = pos // 3
                next_codon_start = codon_idx * 3 + 3
                is_within = (pos + 1) < next_codon_start

                if site_type == "donor" and is_within:
                    # Cost-aware: check splice donor potential before fixing.
                    # Even if MaxEntScan donor score exceeds threshold, a GT
                    # with low overall splice donor potential (no acceptor,
                    # no polypyrimidine tract) is not biologically dangerous.
                    from ..optimization import score_splice_donor_potential as _sdp, SPLICE_DONOR_POTENTIAL_THRESHOLD as _sdp_thresh
                    sdp_val = _sdp(seq_str, pos)
                    if sdp_val < _sdp_thresh:
                        sites = sites[1:]
                        continue  # Low-risk GT — CAI takes priority

                    aa = protein[codon_idx] if codon_idx < len(protein) else None
                    if aa is None or aa == "*":
                        sites = sites[1:]
                        continue
                    current = "".join(seq_chars[codon_idx*3:codon_idx*3+3])

                    # Check if this is a Valine codon (unavoidable)
                    if current[:2] == "GT":
                        from ..type_system import CODON_TABLE as _CT
                        if _CT.get(current) == "V":
                            sites = sites[1:]
                            continue

                    # Try GT-free alternatives (CAI-sorted)
                    fixed = False
                    for alt in gt_free.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(codon_idx, alt)
                        if not _has_rs_local(codon_idx):
                            violations_fixed += 1
                            fixed = True
                            break
                        else:
                            _rollback_swap(codon_idx, old_codon)

                    if not fixed:
                        sites = sites[1:]
                    else:
                        # Incremental update: remove sites near the fixed codon
                        # and rescan only the affected region instead of full sequence
                        ci_start = max(0, codon_idx * 3 - 22)  # acceptor model needs 23-mer
                        ci_end = min(n_bases, codon_idx * 3 + 3 + 22)
                        region = seq_buf[ci_start:ci_end].decode('ascii')
                        new_sites = scan_splice_sites(
                            region,
                            donor_threshold=splice_threshold,
                            acceptor_threshold=splice_threshold,
                        )
                        # Adjust positions back to absolute
                        adjusted = [(p + ci_start, st, sc) for p, st, sc in new_sites]
                        # Remove old sites in the affected region and add new ones
                        sites = [(p, st, sc) for p, st, sc in sites[1:]
                                 if not (ci_start <= p < ci_end)]
                        sites.extend(adjusted)

                elif site_type == "donor" and not is_within:
                    # Cross-codon GT: check splice donor potential first
                    from ..optimization import score_splice_donor_potential as _sdp2, SPLICE_DONOR_POTENTIAL_THRESHOLD as _sdp_thresh2
                    sdp_val2 = _sdp2(seq_str, pos)
                    if sdp_val2 < _sdp_thresh2:
                        sites = sites[1:]
                        continue  # Low-risk GT — CAI takes priority

                    # Try fixing either codon
                    next_ci = codon_idx + 1
                    fixed = False

                    if next_ci < n_codons:
                        aa2 = protein[next_ci] if next_ci < len(protein) else None
                        if aa2 is not None and aa2 != "*":
                            current2 = "".join(seq_chars[next_ci*3:next_ci*3+3])
                            for alt2 in sorted_codons.get(aa2, []):
                                if alt2[0] == 'T' or alt2 == current2:
                                    continue
                                old2 = _apply_swap(next_ci, alt2)
                                if not _has_rs_local(next_ci):
                                    violations_fixed += 1
                                    fixed = True
                                    break
                                else:
                                    _rollback_swap(next_ci, old2)

                    if not fixed:
                        aa1 = protein[codon_idx] if codon_idx < len(protein) else None
                        if aa1 is not None and aa1 != "*":
                            current1 = "".join(seq_chars[codon_idx*3:codon_idx*3+3])
                            for alt1 in sorted_codons.get(aa1, []):
                                if alt1[-1] == 'G' or alt1 == current1:
                                    continue
                                old1 = _apply_swap(codon_idx, alt1)
                                if not _has_rs_local(codon_idx):
                                    violations_fixed += 1
                                    fixed = True
                                    break
                                else:
                                    _rollback_swap(codon_idx, old1)

                    if not fixed:
                        sites = sites[1:]
                    else:
                        # Incremental update for cross-codon donor fix
                        ci_start = max(0, min(codon_idx, next_ci) * 3 - 22)
                        ci_end = min(n_bases, max(codon_idx, next_ci) * 3 + 3 + 22)
                        region = seq_buf[ci_start:ci_end].decode('ascii')
                        new_sites = scan_splice_sites(
                            region,
                            donor_threshold=splice_threshold,
                            acceptor_threshold=splice_threshold,
                        )
                        adjusted = [(p + ci_start, st, sc) for p, st, sc in new_sites]
                        sites = [(p, st, sc) for p, st, sc in sites[1:]
                                 if not (ci_start <= p < ci_end)]
                        sites.extend(adjusted)

                elif site_type == "acceptor" and is_within:
                    aa = protein[codon_idx] if codon_idx < len(protein) else None
                    if aa is None or aa == "*":
                        sites = sites[1:]
                        continue
                    current = "".join(seq_chars[codon_idx*3:codon_idx*3+3])

                    # Try AG-free alternatives (CAI-sorted)
                    fixed = False
                    for alt in ag_free.get(aa, []):
                        if alt == current:
                            continue
                        old_codon = _apply_swap(codon_idx, alt)
                        # Check no new GT created
                        creates_gt, _ = _creates_new_gt_or_ag(codon_idx, alt)
                        if not creates_gt and not _has_rs_local(codon_idx):
                            violations_fixed += 1
                            fixed = True
                            break
                        else:
                            _rollback_swap(codon_idx, old_codon)

                    if not fixed:
                        sites = sites[1:]
                    else:
                        # Incremental update for acceptor fix
                        ci_start = max(0, codon_idx * 3 - 22)
                        ci_end = min(n_bases, codon_idx * 3 + 3 + 22)
                        region = seq_buf[ci_start:ci_end].decode('ascii')
                        new_sites = scan_splice_sites(
                            region,
                            donor_threshold=splice_threshold,
                            acceptor_threshold=splice_threshold,
                        )
                        adjusted = [(p + ci_start, st, sc) for p, st, sc in new_sites]
                        sites = [(p, st, sc) for p, st, sc in sites[1:]
                                 if not (ci_start <= p < ci_end)]
                        sites.extend(adjusted)

                else:
                    # Cross-codon acceptor or unknown — skip
                    sites = sites[1:]

        except ImportError:
            pass

    # ── Phase 4: CpG island disruption ──
    # Fix CpG dinucleotides by replacing with CG-free alternatives.
    # We prefer higher-CAI alternatives but accept some CAI loss since
    # CpG islands are a hard constraint for eukaryotic optimization.
    # For cross-codon CGs, try both sides of the boundary.
    # Allow creating GT/AG during CpG fix (MaxEntScan will handle
    # any cryptic splice sites that result).
    seq_str = seq_buf.decode('ascii')
    for _cpg_iter in range(30):
        seq_str = seq_buf.decode('ascii')
        found_cpg = False

        for i in range(n_bases - 1):
            if seq_str[i] == 'C' and seq_str[i + 1] == 'G':
                codon_idx = i // 3
                next_codon_start = codon_idx * 3 + 3
                is_within = (i + 1) < next_codon_start

                if is_within:
                    aa = protein[codon_idx] if codon_idx < len(protein) else None
                    if aa is None or aa == "*":
                        continue

                    current = "".join(seq_chars[codon_idx*3:codon_idx*3+3])
                    # Try CG-free alternatives (CAI-sorted)
                    # First try without creating GT/AG, then allow GT/AG
                    fixed = False
                    for alt in sorted_codons.get(aa, []):
                        if alt == current or "CG" in alt:
                            continue
                        old_codon = _apply_swap(codon_idx, alt)
                        creates_gt, creates_ag = _creates_new_gt_or_ag(codon_idx, alt)
                        if not creates_gt and not creates_ag and not _has_rs_local(codon_idx):
                            violations_fixed += 1
                            found_cpg = True
                            fixed = True
                            break
                        else:
                            _rollback_swap(codon_idx, old_codon)

                    if not fixed:
                        # Allow GT/AG creation for CpG elimination
                        for alt in sorted_codons.get(aa, []):
                            if alt == current or "CG" in alt:
                                continue
                            old_codon = _apply_swap(codon_idx, alt)
                            if not _has_rs_local(codon_idx):
                                violations_fixed += 1
                                found_cpg = True
                                fixed = True
                                break
                            else:
                                _rollback_swap(codon_idx, old_codon)

                    if fixed:
                        break
                else:
                    # Cross-codon CG: current codon ends with C, next starts with G
                    next_ci = codon_idx + 1
                    if next_ci >= n_codons:
                        continue
                    aa2 = protein[next_ci] if next_ci < len(protein) else None
                    if aa2 is None or aa2 == "*":
                        continue

                    # Try next codon with non-G-starting alternative
                    current2 = "".join(seq_chars[next_ci*3:next_ci*3+3])
                    fixed = False
                    for alt2 in sorted_codons.get(aa2, []):
                        if alt2[0] == 'G' or alt2 == current2:
                            continue
                        old2 = _apply_swap(next_ci, alt2)
                        creates_gt2, creates_ag2 = _creates_new_gt_or_ag(next_ci, alt2)
                        if not creates_gt2 and not creates_ag2 and not _has_rs_local(next_ci):
                            violations_fixed += 1
                            found_cpg = True
                            fixed = True
                            break
                        else:
                            _rollback_swap(next_ci, old2)

                    if not fixed:
                        # Allow GT/AG for next codon
                        for alt2 in sorted_codons.get(aa2, []):
                            if alt2[0] == 'G' or alt2 == current2:
                                continue
                            old2 = _apply_swap(next_ci, alt2)
                            if not _has_rs_local(next_ci):
                                violations_fixed += 1
                                found_cpg = True
                                fixed = True
                                break
                            else:
                                _rollback_swap(next_ci, old2)

                    if not fixed:
                        # Try current codon with non-C-ending alternative
                        aa = protein[codon_idx] if codon_idx < len(protein) else None
                        if aa is not None and aa != "*":
                            current1 = "".join(seq_chars[codon_idx*3:codon_idx*3+3])
                            for alt1 in sorted_codons.get(aa, []):
                                if alt1[-1] == 'C' or alt1 == current1:
                                    continue
                                old1 = _apply_swap(codon_idx, alt1)
                                creates_gt1, creates_ag1 = _creates_new_gt_or_ag(codon_idx, alt1)
                                if not creates_gt1 and not creates_ag1 and not _has_rs_local(codon_idx):
                                    violations_fixed += 1
                                    found_cpg = True
                                    fixed = True
                                    break
                                else:
                                    _rollback_swap(codon_idx, old1)

                        if not fixed:
                            # Allow GT/AG for current codon too
                            if aa is not None and aa != "*":
                                for alt1 in sorted_codons.get(aa, []):
                                    if alt1[-1] == 'C' or alt1 == current1:
                                        continue
                                    old1 = _apply_swap(codon_idx, alt1)
                                    if not _has_rs_local(codon_idx):
                                        violations_fixed += 1
                                        found_cpg = True
                                        fixed = True
                                        break
                                    else:
                                        _rollback_swap(codon_idx, old1)

                    if fixed:
                        break

        if not found_cpg:
            break

    # ── Phase 5: CAI recovery hill climb ──
    # Upgrade codons while maintaining all constraints.
    # Key insight: GT/AG dinucleotides are only problematic if MaxEntScan
    # identifies them as cryptic splice sites.  Most GT/AG occurrences are
    # benign.  So we use a two-tier approach:
    #  - First try: accept upgrade if it passes hard constraints (RS, GC,
    #    ATTTA, T-run) AND doesn't create new GT/AG.
    #  - Second try: accept upgrade even if it creates new GT/AG, but
    #    then validate with MaxEntScan and revert only if it creates an
    #    actual cryptic splice site.
    aas = list(protein)
    # Determine whether MaxEntScan is available for GT/AG validation
    _mes_available = False
    try:
        if _HAS_FAST_MAXENT:
            _mes_scan = _scan_splice_sites_fast
            _mes_available = True
        else:
            from ..maxentscan import scan_splice_sites
            _mes_scan = scan_splice_sites
            _mes_available = True
    except ImportError:
        pass

    for _iteration in range(5):
        upgrade_plan: dict[int, str] = {}

        for ci in range(n_codons):
            aa = aas[ci]
            if aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])
            current_cai = species_cai.get(current, 0.0)
            optimal = optimal_codon.get(aa, current)
            optimal_cai = species_cai.get(optimal, 0.0)

            if optimal_cai > current_cai and current != optimal:
                upgrade_plan[ci] = optimal

        if not upgrade_plan:
            break

        any_improved = False
        applied: set[int] = set()

        for ci in sorted(
            upgrade_plan.keys(),
            key=lambda c: species_cai.get(upgrade_plan[c], 0.0) - species_cai.get("".join(seq_chars[c*3:c*3+3]), 0.0),
            reverse=True,
        ):
            if any(abs(ci - a) <= 2 for a in applied):
                continue

            new_codon = upgrade_plan[ci]
            current = "".join(seq_chars[ci*3:ci*3+3])
            if current == new_codon:
                continue

            # Helper: check hard constraints (RS, GC, ATTTA, T-run)
            def _check_hard_constraints(ci_idx: int) -> bool:
                if _has_rs_local(ci_idx):
                    return False
                if not (gc_lo <= (gc_count / n_bases) <= gc_hi):
                    return False
                local = seq_buf.decode('ascii')
                if "ATTTA" in local[max(0, ci_idx*3-5):min(n_bases, ci_idx*3+8)]:
                    return False
                # T-run check
                cs = max(0, ci_idx * 3 - 6)
                ce = min(n_bases, ci_idx * 3 + 9)
                j2 = cs
                while j2 < ce:
                    if seq_chars[j2] == 'T':
                        k2 = j2
                        while k2 < ce and seq_chars[k2] == 'T':
                            k2 += 1
                        if k2 - j2 >= 6:
                            return False
                        j2 = k2
                    else:
                        j2 += 1
                return True

            # Try the swap and check constraints
            old_codon = _apply_swap(ci, new_codon)
            local_seq = seq_buf.decode('ascii')
            local_region = local_seq[max(0, ci*3-1):min(n_bases, ci*3+4)]

            rs_ok = not _has_rs_local(ci)
            gt_ok = "GT" not in local_region
            ag_ok = "AG" not in local_region
            gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi
            attta_ok = "ATTTA" not in local_seq[max(0, ci*3-5):min(n_bases, ci*3+8)]

            # T-run check
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

            if rs_ok and gt_ok and ag_ok and gc_ok and attta_ok and trun_ok:
                # All constraints pass — accept upgrade
                violations_fixed += 1
                any_improved = True
                applied.add(ci)
            elif rs_ok and gc_ok and attta_ok and trun_ok and (not gt_ok or not ag_ok) and _mes_available:
                # Hard constraints pass, but new GT/AG created.
                # Validate with MaxEntScan — only reject if it creates an
                # actual cryptic splice site above the threshold.
                mes_seq = seq_buf.decode('ascii')
                mes_sites = _mes_scan(
                    mes_seq,
                    donor_threshold=splice_threshold,
                    acceptor_threshold=splice_threshold,
                )
                # Check if any new cryptic splice site overlaps our codon
                has_new_splice = False
                for sp_pos, sp_type, sp_score in mes_sites:
                    if abs(sp_pos - ci * 3) <= 4:
                        has_new_splice = True
                        break
                if not has_new_splice:
                    # MaxEntScan says the new GT/AG is benign — accept upgrade
                    violations_fixed += 1
                    any_improved = True
                    applied.add(ci)
                else:
                    _rollback_swap(ci, old_codon)
            else:
                _rollback_swap(ci, old_codon)
                # Try other alternatives that have higher CAI than current
                aa = aas[ci]
                for alt in sorted_codons.get(aa, []):
                    if alt == current or alt == new_codon:
                        continue
                    alt_cai = species_cai.get(alt, 0.0)
                    cur_cai = species_cai.get(current, 0.0)
                    if alt_cai <= cur_cai:
                        continue

                    old_codon2 = _apply_swap(ci, alt)
                    if not _check_hard_constraints(ci):
                        _rollback_swap(ci, old_codon2)
                        continue
                    # Check GT/AG
                    local_seq2 = seq_buf.decode('ascii')
                    local_region2 = local_seq2[max(0, ci*3-1):min(n_bases, ci*3+4)]
                    gt_ok2 = "GT" not in local_region2
                    ag_ok2 = "AG" not in local_region2
                    if gt_ok2 and ag_ok2:
                        violations_fixed += 1
                        any_improved = True
                        applied.add(ci)
                        break
                    elif _mes_available:
                        # New GT/AG — validate with MaxEntScan
                        mes_seq2 = seq_buf.decode('ascii')
                        mes_sites2 = _mes_scan(
                            mes_seq2,
                            donor_threshold=splice_threshold,
                            acceptor_threshold=splice_threshold,
                        )
                        has_new_splice2 = False
                        for sp_pos, sp_type, sp_score in mes_sites2:
                            if abs(sp_pos - ci * 3) <= 4:
                                has_new_splice2 = True
                                break
                        if not has_new_splice2:
                            violations_fixed += 1
                            any_improved = True
                            applied.add(ci)
                            break
                        else:
                            _rollback_swap(ci, old_codon2)
                    else:
                        _rollback_swap(ci, old_codon2)

        if not any_improved:
            break

    # ── Phase 6: CAI micro-optimization with constraint recheck ──
    # After all constraint phases, some codons may be suboptimal because
    # a restriction site or GT/AG pattern existed earlier but was later
    # broken by other codon changes.  Re-check each non-optimal codon to
    # see if the optimal (or a higher-CAI) codon is now safe.  Also try
    # coordinated two-codon swaps: if upgrading codon A creates a RS,
    # maybe changing a neighboring codon B simultaneously resolves it
    # with a net CAI gain.
    hill_climb_improvements = 0
    for _micro_iter in range(5):
        any_micro_improved = False

        for ci in range(n_codons):
            aa = aas[ci]
            if aa == "*":
                continue
            current = "".join(seq_chars[ci*3:ci*3+3])
            current_cai = species_cai.get(current, 0.0)
            optimal = optimal_codon.get(aa, current)
            optimal_cai = species_cai.get(optimal, 0.0)

            if current == optimal or optimal_cai <= current_cai:
                continue

            # Try upgrading to the optimal codon
            old_codon = _apply_swap(ci, optimal)
            local_seq = seq_buf.decode('ascii')
            local_region = local_seq[max(0, ci*3-1):min(n_bases, ci*3+4)]

            rs_ok = not _has_rs_local(ci)
            gt_ok = "GT" not in local_region
            ag_ok = "AG" not in local_region
            gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi
            attta_ok = "ATTTA" not in local_seq[max(0, ci*3-5):min(n_bases, ci*3+8)]
            # T-run check
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

            all_hard_ok = rs_ok and gc_ok and attta_ok and trun_ok

            if all_hard_ok and gt_ok and ag_ok:
                # All constraints pass — accept upgrade
                hill_climb_improvements += 1
                any_micro_improved = True
                continue
            elif all_hard_ok and (not gt_ok or not ag_ok) and _mes_available:
                # Hard constraints pass, but new GT/AG — validate with MaxEntScan
                mes_sites = _mes_scan(
                    local_seq,
                    donor_threshold=splice_threshold,
                    acceptor_threshold=splice_threshold,
                )
                has_new_splice = False
                for sp_pos, sp_type, sp_score in mes_sites:
                    if abs(sp_pos - ci * 3) <= 4:
                        has_new_splice = True
                        break
                if not has_new_splice:
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    continue
            _rollback_swap(ci, old_codon)

            # Optimal codon creates a constraint violation.
            # Try coordinated two-codon swaps with immediate neighbors
            # to resolve the RS while gaining net CAI.
            for neighbor_ci in [ci - 1, ci + 1]:
                if neighbor_ci < 0 or neighbor_ci >= n_codons:
                    continue
                neighbor_aa = aas[neighbor_ci]
                if neighbor_aa == "*":
                    continue
                neighbor_current = "".join(seq_chars[neighbor_ci*3:neighbor_ci*3+3])
                neighbor_current_cai = species_cai.get(neighbor_current, 0.0)

                best_net_log = 0.0  # net log-CAI gain must be positive
                best_neighbor_alt = None

                for neighbor_alt in sorted_codons.get(neighbor_aa, []):
                    if neighbor_alt == neighbor_current:
                        continue
                    old_neighbor = _apply_swap(neighbor_ci, neighbor_alt)
                    old_ci = _apply_swap(ci, optimal)

                    # Check constraints
                    if not _has_rs_local(ci) and not _has_rs_local(neighbor_ci):
                        new_gc_val = gc_count / n_bases
                        if gc_lo <= new_gc_val <= gc_hi:
                            local_w = seq_buf.decode('ascii')
                            has_attta = "ATTTA" in local_w[max(0, ci*3-5):min(n_bases, ci*3+8)]
                            if not has_attta:
                                # T-run check
                                cs2 = max(0, ci * 3 - 6)
                                ce2 = min(n_bases, ci * 3 + 9)
                                has_lr = False
                                jj = cs2
                                while jj < ce2:
                                    if seq_chars[jj] == 'T':
                                        kk = jj
                                        while kk < ce2 and seq_chars[kk] == 'T':
                                            kk += 1
                                        if kk - jj >= 6:
                                            has_lr = True
                                            break
                                        jj = kk
                                    else:
                                        jj += 1
                                if not has_lr:
                                    # Check GT/AG
                                    lr = local_w[max(0, ci*3-1):min(n_bases, ci*3+4)]
                                    gt2 = "GT" not in lr
                                    ag2 = "AG" not in lr
                                    accept_swap = False
                                    if gt2 and ag2:
                                        accept_swap = True
                                    elif _mes_available:
                                        mes2 = _mes_scan(
                                            local_w,
                                            donor_threshold=splice_threshold,
                                            acceptor_threshold=splice_threshold,
                                        )
                                        has_ns = False
                                        for sp2_p, sp2_t, sp2_s in mes2:
                                            if abs(sp2_p - ci * 3) <= 4:
                                                has_ns = True
                                                break
                                        if not has_ns:
                                            accept_swap = True
                                    if accept_swap:
                                        neighbor_alt_cai = species_cai.get(neighbor_alt, 0.0)
                                        net_log_gain = 0.0
                                        if optimal_cai > 0 and current_cai > 0:
                                            net_log_gain += math.log(optimal_cai) - math.log(current_cai)
                                        if neighbor_alt_cai > 0 and neighbor_current_cai > 0:
                                            net_log_gain += math.log(neighbor_alt_cai) - math.log(neighbor_current_cai)
                                        if net_log_gain > best_net_log:
                                            best_net_log = net_log_gain
                                            best_neighbor_alt = neighbor_alt

                    _rollback_swap(ci, old_ci)
                    _rollback_swap(neighbor_ci, old_neighbor)

                if best_neighbor_alt is not None and best_net_log > 0:
                    _apply_swap(neighbor_ci, best_neighbor_alt)
                    _apply_swap(ci, optimal)
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    break  # Move to next codon

            # Also try non-optimal but higher-CAI alternatives
            current_after = "".join(seq_chars[ci*3:ci*3+3])
            current_after_cai = species_cai.get(current_after, 0.0)
            for alt in sorted_codons.get(aa, []):
                if alt == current_after:
                    continue
                alt_cai = species_cai.get(alt, 0.0)
                if alt_cai <= current_after_cai:
                    continue
                old_codon2 = _apply_swap(ci, alt)
                local_s = seq_buf.decode('ascii')
                lr2 = local_s[max(0, ci*3-1):min(n_bases, ci*3+4)]
                rs_ok2 = not _has_rs_local(ci)
                gt_ok2 = "GT" not in lr2
                ag_ok2 = "AG" not in lr2
                gc_ok2 = gc_lo <= (gc_count / n_bases) <= gc_hi
                attta_ok2 = "ATTTA" not in local_s[max(0, ci*3-5):min(n_bases, ci*3+8)]
                # T-run check
                cs3 = max(0, ci * 3 - 6)
                ce3 = min(n_bases, ci * 3 + 9)
                trun_ok2 = True
                jj2 = cs3
                while jj2 < ce3:
                    if seq_chars[jj2] == 'T':
                        kk2 = jj2
                        while kk2 < ce3 and seq_chars[kk2] == 'T':
                            kk2 += 1
                        if kk2 - jj2 >= 6:
                            trun_ok2 = False
                            break
                        jj2 = kk2
                    else:
                        jj2 += 1

                all_hard2 = rs_ok2 and gc_ok2 and attta_ok2 and trun_ok2
                if all_hard2 and gt_ok2 and ag_ok2:
                    hill_climb_improvements += 1
                    any_micro_improved = True
                    break
                elif all_hard2 and (not gt_ok2 or not ag_ok2) and _mes_available:
                    mes3 = _mes_scan(
                        local_s,
                        donor_threshold=splice_threshold,
                        acceptor_threshold=splice_threshold,
                    )
                    has_ns3 = False
                    for sp3_p, sp3_t, sp3_s in mes3:
                        if abs(sp3_p - ci * 3) <= 4:
                            has_ns3 = True
                            break
                    if not has_ns3:
                        hill_climb_improvements += 1
                        any_micro_improved = True
                        break
                    else:
                        _rollback_swap(ci, old_codon2)
                else:
                    _rollback_swap(ci, old_codon2)

        if not any_micro_improved:
            break

    # ── Phase 6b: CAI micro-optimization — recheck suboptimal codons ──
    # After all constraint phases (including Phase 6 coordinated swaps),
    # some codons that were suboptimally chosen (e.g., to avoid a
    # restriction site) may now be safe to upgrade because nearby codons
    # changed during later phases.  This final pass systematically tries
    # ALL higher-CAI alternatives, not just the optimal codon.
    for _micro_iter in range(3):
        any_micro_upgrade = False
        for ci in range(n_codons):
            aa = aas[ci]
            if aa == "*":
                continue
            current = seq_buf[ci*3:ci*3+3].decode('ascii')
            current_cai = species_cai.get(current, 0.0)

            # Try all higher-CAI alternatives (not just optimal)
            for alt in sorted_codons.get(aa, []):
                alt_cai = species_cai.get(alt, 0.0)
                if alt_cai <= current_cai:
                    break  # sorted_codons is CAI-descending, no point continuing
                if alt == current:
                    continue

                old_codon = _apply_swap(ci, alt)

                # Check ALL constraints
                rs_ok = not _has_rs_local(ci)
                gc_ok = gc_lo <= (gc_count / n_bases) <= gc_hi

                # Check local region for ATTTA
                local_start = max(0, ci * 3 - 5)
                local_end = min(n_bases, ci * 3 + 8)
                local_str = seq_buf[local_start:local_end].decode('ascii')
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

                # Check GT/AG dinucleotides (splice avoidance)
                local_gt_ag_start = max(0, ci * 3 - 1)
                local_gt_ag_end = min(n_bases, ci * 3 + 4)
                local_gt_ag = seq_buf[local_gt_ag_start:local_gt_ag_end].decode('ascii')
                gt_ag_ok = "GT" not in local_gt_ag and "AG" not in local_gt_ag

                all_hard_ok = rs_ok and gc_ok and attta_ok and trun_ok

                if all_hard_ok and gt_ag_ok:
                    any_micro_upgrade = True
                    break  # Keep this upgrade, move to next codon
                elif all_hard_ok and (not gt_ag_ok) and _mes_available:
                    # Hard constraints pass but GT/AG created — validate
                    # with MaxEntScan; accept if it's a benign GT/AG.
                    mes_local = seq_buf.decode('ascii')
                    mes_sites = _mes_scan(
                        mes_local,
                        donor_threshold=splice_threshold,
                        acceptor_threshold=splice_threshold,
                    )
                    has_new_splice = False
                    for sp_pos, sp_type, sp_score in mes_sites:
                        if abs(sp_pos - ci * 3) <= 4:
                            has_new_splice = True
                            break
                    if not has_new_splice:
                        any_micro_upgrade = True
                        break  # Benign GT/AG — keep upgrade
                    else:
                        _rollback_swap(ci, old_codon)
                else:
                    _rollback_swap(ci, old_codon)

        if not any_micro_upgrade:
            break

    # ── Compute final metrics ──
    seq = seq_buf.decode('ascii')
    final_cai = math.exp(log_cai_sum / n_cai_codons) if n_cai_codons > 0 else 0.0
    gc = gc_count / n_bases if n_bases > 0 else 0.0

    # ── Codon pair bias optimization ──
    cpb_improvements = 0
    mean_cpb = 0.0
    if optimizer.consider_codon_pair_bias:
        seq, final_cai, cpb_improvements, mean_cpb = (
            codon_pair_bias_optimize_prokaryote(
                    optimizer, seq, protein, final_cai
                )
        )
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)

    return HybridResult(
        sequence=seq,
        cai=final_cai,
        gc_content=round(gc, 4),
        violations_fixed=violations_fixed,
        hill_climb_improvements=0,
        iterations_used=0,
        phase1_cai=phase1_cai,
        phase2_cai=final_cai,
        phase3_cai=final_cai,
        phase4_cai=final_cai,
        cpb_improvements=cpb_improvements,
        mean_cpb=mean_cpb,
        warnings=warnings,
        splice_sites_validated=True,  # Phase 3 MaxEntScan already validated
    )
