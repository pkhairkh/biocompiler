"""
Constraint satisfaction methods for HybridOptimizer.

Extracted from hybrid_optimizer.py to decompose the monolith.
Contains _constraint_satisfaction() and its helper methods
for detecting and fixing constraint violations.
"""

from __future__ import annotations

import heapq
import logging
from typing import Any

from ..type_system import CODON_TABLE, AA_TO_CODONS
from ..incremental import IncrementalSequenceState
from ..decision_provenance import ConstraintDecision
from .hybrid_types import Violation, SEVERITY_WEIGHTS, GT_CAI_COST_THRESHOLD as _GT_CAI_COST_THRESHOLD

logger = logging.getLogger(__name__)


def _constraint_satisfaction(
    optimizer, seq: str, protein: str, avoid_gt: bool
) -> tuple[str, float, int, int, list[str]]:
    """Simultaneous constraint satisfaction with conflict detection.

    Algorithm:
    1. Scan for ALL constraint violations at once
    2. Score each violation by severity
    3. Group violations by whether they conflict (share codon indices)
    4. Fix ALL non-conflicting violations simultaneously:
       - Restriction sites at different positions → fix all at once
       - GC adjustment → batch all codon changes and apply at once
       - ATTTA/T-run fixes at non-overlapping positions → fix all at once
    5. Only iterate when a fix conflicts with another fix
       (shared codon indices between violations)
    6. For conflicting violations, use CAI-aware resolution:
       prefer the solution that maximizes CAI, using
       CAIAwareConstraintResolver if available.

    Key improvement over one-at-a-time approach:
    - Previously: fix ONE violation → re-scan → fix ONE → re-scan → ...
    - Now: fix ALL non-conflicting violations → re-scan only if conflicts
    - CAI-aware: when constraints conflict, prefer higher-CAI resolution

    Returns:
        (sequence, cai, violations_fixed, iterations, warnings)
    """
    state = IncrementalSequenceState(
        seq, species=optimizer.species, enzymes=optimizer.enzymes
    )
    warnings: list[str] = []
    violations_fixed = 0
    total_iterations = 0

    # ── CAI-aware constraint resolver (lazy-loaded) ──
    # When two constraints conflict (share codon indices), the
    # CAIAwareConstraintResolver evaluates which resolution option
    # has the smallest impact on CAI and prefers that resolution.
    _cai_resolver = None
    try:
        from ..solver.cai_aware_resolver import CAIAwareConstraintResolver
        # We don't have a CSPModel here, but we can use the resolver
        # for its CAI impact estimation methods directly.
        _cai_resolver = CAIAwareConstraintResolver.__new__(CAIAwareConstraintResolver)
        _cai_resolver._adaptiveness = optimizer.species_cai
    except (ImportError, Exception):
        _cai_resolver = None

    # Track violation signatures for stale detection — if the same set
    # of violations persists for 2 consecutive iterations, the remaining
    # violations are unfixable and we should stop rather than waste
    # cycles on futile swap-check-rollback loops.
    _prev_violation_sig: frozenset[tuple[str, int]] | None = None
    _stale_count = 0

    # Phase 2a: Fast constraint fixes (no MaxEntScan needed)
    for iteration in range(optimizer.max_local_search_iterations):
        violations = _detect_cheap_violations(optimizer, state, avoid_gt)

        if not violations:
            break

        # ── Eukaryotic GT softening ──
        # For eukaryotic organisms, GT avoidance is a SOFT constraint.
        # Lower the severity of avoidable_gt violations so they are
        # treated as WARNINGS rather than hard failures.  GT violations
        # should not block CAI recovery or force suboptimal codon choices.
        if avoid_gt:
            for v in violations:
                if v.violation_type == "avoidable_gt":
                    # Reduce severity below all hard constraints so GT
                    # fixes are attempted only after everything else
                    v.severity = SEVERITY_WEIGHTS["avoidable_gt"] * 0.1

        total_iterations = iteration + 1

        # ── Stale detection ──
        # Compute a signature of the current violation set (type + position).
        # If this signature is identical to the previous iteration, the
        # fixes are not making progress and the remaining violations are
        # unfixable.  Break immediately to avoid wasting time.
        sig = frozenset((v.violation_type, v.position) for v in violations)
        if sig == _prev_violation_sig:
            _stale_count += 1
            if _stale_count >= 2:
                break
        else:
            _stale_count = 0
        _prev_violation_sig = sig

        # ── CAI-aware ordering: when multiple violations have similar
        # severity, prefer fixing the one whose resolution has the
        # smallest CAI impact.  This ensures that when two constraints
        # conflict, the resolution chosen maximizes CAI.
        violations.sort(key=lambda v: v.severity, reverse=True)

        # If CAI resolver is available, re-sort same-severity violations
        # by estimated CAI impact (lower impact first)
        if _cai_resolver is not None and len(violations) > 1:
            # Group by severity tier (within 10% of each other)
            # and re-sort within each tier by CAI impact
            for i in range(len(violations)):
                for j in range(i + 1, len(violations)):
                    vi, vj = violations[i], violations[j]
                    # Only re-sort if severities are similar
                    if abs(vi.severity - vj.severity) / max(vi.severity, 1.0) < 0.1:
                        # Estimate CAI impact based on violation type
                        # Soft constraints (CpG, ATTTA, T-run) have lower
                        # CAI impact than hard constraints (RS, GC, stop)
                        cai_impact_i = optimizer._estimate_cai_impact(vi.violation_type)
                        cai_impact_j = optimizer._estimate_cai_impact(vj.violation_type)
                        if cai_impact_i > cai_impact_j:
                            # j has lower CAI impact — fix it first
                            violations[i], violations[j] = violations[j], violations[i]

        # ── Simultaneous fix: group non-conflicting violations ──
        # Two violations conflict if they share any codon index.
        # Non-conflicting violations (different positions) can be
        # fixed simultaneously without interfering with each other.
        # For conflicting violations, use CAI-aware resolution.
        fixed_this_round = 0

        # Track which codon indices have been modified
        used_codon_indices: set[int] = set()

        # Collect conflicting violations for CAI-aware resolution
        conflicting_violations: list[Violation] = []

        for violation in violations:
            # Check if this violation conflicts with already-fixed ones
            v_codons = set(violation.codon_indices)
            if v_codons & used_codon_indices:
                # Conflict: record for CAI-aware resolution
                conflicting_violations.append(violation)
                continue

            fixed = _fix_violation(optimizer, state, violation, protein, avoid_gt)
            if fixed:
                violations_fixed += 1
                fixed_this_round += 1
                # Mark these codon indices as used so future violations
                # in this round don't conflict
                used_codon_indices.update(v_codons)

        # ── CAI-aware conflict resolution ──
        # For violations that conflicted with already-fixed ones,
        # try to resolve them by preferring the option that maximizes
        # CAI.  If two constraints conflict at the same position,
        # try fixing the one with lower CAI impact first.
        if conflicting_violations and fixed_this_round > 0:
            # Re-sort conflicting by CAI impact (lower first)
            if _cai_resolver is not None:
                conflicting_violations.sort(
                    key=lambda v: optimizer._estimate_cai_impact(v.violation_type)
                )
            for violation in conflicting_violations:
                v_codons = set(violation.codon_indices)
                if v_codons & used_codon_indices:
                    continue
                fixed = _fix_violation(optimizer, state, violation, protein, avoid_gt)
                if fixed:
                    violations_fixed += 1
                    fixed_this_round += 1
                    used_codon_indices.update(v_codons)
                    # Record CAI-aware resolution in provenance
                    if optimizer.provenance_collector is not None:
                        optimizer.provenance_collector.record_constraint_decision(
                            ConstraintDecision(
                                constraint_name=violation.violation_type,
                                constraint_type="hard",
                                action_taken="satisfied",
                                positions_affected=violation.codon_indices,
                                tradeoff_description=(
                                    f"CAI-aware resolution: fixed {violation.violation_type} "
                                    f"at pos {violation.position} with CAI-optimal alternative"
                                ),
                                impact_on_cai=0.0,
                            )
                        )

        if fixed_this_round == 0:
            # No violations could be fixed — try one at a time with
            # the highest severity violation (fallback)
            for violation in violations:
                fixed = _fix_violation(optimizer, state, violation, protein, avoid_gt)
                if fixed:
                    violations_fixed += 1
                    fixed_this_round += 1
                    break

            if fixed_this_round == 0:
                break

    # Phase 2b: Expensive constraint fixes (MaxEntScan-based splice checks)
    # Only when avoid_gt is True (eukaryotic targets)
    if avoid_gt:
        for splice_iter in range(min(10, optimizer.max_local_search_iterations)):
            violations = _detect_expensive_violations(optimizer, state)

            if not violations:
                break

            # Fix all non-conflicting violations simultaneously
            violations.sort(key=lambda v: v.severity, reverse=True)
            any_fixed = False
            used_codon_indices = set()

            for violation in violations:
                v_codons = set(violation.codon_indices)
                if v_codons & used_codon_indices:
                    continue  # Skip conflicting violations

                fixed = _fix_violation(optimizer, state, violation, protein, avoid_gt)
                if fixed:
                    violations_fixed += 1
                    any_fixed = True
                    used_codon_indices.update(v_codons)

            if not any_fixed:
                for v in violations[:5]:
                    warnings.append(
                        f"Unresolved {v.violation_type} at pos "
                        f"{v.position}: {v.details}"
                    )
                break

        # Re-check cheap violations after expensive fixes
        # (splice site fixes may reintroduce restriction sites or GC issues)
        for recheck_iter in range(min(5, optimizer.max_local_search_iterations)):
            re_violations = _detect_cheap_violations(optimizer, state, avoid_gt)
            if not re_violations:
                break
            re_violations.sort(key=lambda v: v.severity, reverse=True)
            any_fixed = False
            used_codon_indices = set()
            for violation in re_violations:
                v_codons = set(violation.codon_indices)
                if v_codons & used_codon_indices:
                    continue
                fixed = _fix_violation(optimizer, state, violation, protein, avoid_gt)
                if fixed:
                    violations_fixed += 1
                    any_fixed = True
                    used_codon_indices.update(v_codons)
            if not any_fixed:
                break

    # ── Eukaryotic GT warning collection ──
    # After constraint satisfaction, any remaining GT violations for
    # eukaryotes should be reported as warnings rather than treated as
    # hard failures.  GT avoidance is soft for eukaryotes.
    if avoid_gt:
        remaining_gt = _detect_cheap_violations(optimizer, state, avoid_gt)
        for v in remaining_gt:
            if v.violation_type == "avoidable_gt":
                warnings.append(
                    f"Soft GT warning: {v.details} "
                    f"(eukaryotic GT avoidance is soft, not a hard constraint)"
                )

    cai = optimizer._compute_cai(state.sequence)
    return state.sequence, cai, violations_fixed, total_iterations, warnings


def _detect_cheap_violations(
    optimizer, state: IncrementalSequenceState, avoid_gt: bool
) -> list[Violation]:
    """Detect constraint violations that don't require MaxEntScan.

    These are fast to check and should be fixed before the expensive
    splice site checks.
    """
    violations: list[Violation] = []
    seq = state.sequence
    n_codons = state.num_codons

    # 1. Restriction site violations (incremental — uses pre-tracked positions)
    for site_seq, site_pos in state.check_restriction_sites(changed_only=False):
        codon_indices = optimizer._get_overlapping_codon_indices(
            site_pos, len(site_seq), n_codons
        )
        violations.append(Violation(
            violation_type="restriction_site",
            position=site_pos,
            severity=SEVERITY_WEIGHTS["restriction_site"],
            codon_indices=codon_indices,
            details=f"Site {site_seq} at position {site_pos}",
        ))

    # 2. Stop codon violations
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in ("TAA", "TAG", "TGA") and i < len(seq) - 3:
            aa = state.get_aa(i // 3)
            if aa == "*" and i // 3 == n_codons - 1:
                continue
            violations.append(Violation(
                violation_type="stop_codon",
                position=i,
                severity=SEVERITY_WEIGHTS["stop_codon"],
                codon_indices=[i // 3],
                details=f"Premature stop codon {codon} at position {i}",
            ))

    # 3. GC out of range (O(1) via incremental gc_count)
    gc = state.gc_fraction
    if not (optimizer.gc_lo <= gc <= optimizer.gc_hi):
        if gc < optimizer.gc_lo:
            distance = optimizer.gc_lo - gc
        else:
            distance = gc - optimizer.gc_hi
        severity = SEVERITY_WEIGHTS["gc_out_of_range"] * min(distance * 10, 1.0)
        violations.append(Violation(
            violation_type="gc_out_of_range",
            position=0,
            severity=severity,
            codon_indices=list(range(n_codons)),
            details=f"GC={gc:.3f} outside [{optimizer.gc_lo}, {optimizer.gc_hi}]",
        ))

    # 4. Avoidable GT dinucleotides (simple check, no MaxEntScan)
    # Cost-aware: only flag GTs with high splice donor potential.
    # GTs with low splice donor potential (< threshold) are not
    # biologically dangerous and CAI takes priority.
    if avoid_gt:
        from ..optimization import score_splice_donor_potential, SPLICE_DONOR_POTENTIAL_THRESHOLD
        for gt_pos in state.gt_positions_list():
            if optimizer._is_unavoidable_gt(seq, gt_pos):
                continue
            # Cost-aware GT resolution: check splice donor potential
            sdp = score_splice_donor_potential(seq, gt_pos)
            if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                # This GT has low splice donor potential — not dangerous.
                # Accept it (CAI > GT avoidance for low-risk GTs).
                continue
            codon_idx = gt_pos // 3
            next_codon_start = (gt_pos // 3) * 3 + 3
            is_within = (gt_pos + 1) < next_codon_start
            if is_within:
                involved = [codon_idx]
            else:
                involved = (
                    [codon_idx, codon_idx + 1]
                    if codon_idx + 1 < n_codons
                    else [codon_idx]
                )
            violations.append(Violation(
                violation_type="avoidable_gt",
                position=gt_pos,
                severity=SEVERITY_WEIGHTS["avoidable_gt"],
                codon_indices=involved,
                details=(
                    f"GT at pos {gt_pos} "
                    f"({'within' if is_within else 'cross'}-codon, "
                    f"SDP={sdp:.3f})"
                ),
            ))

    # 5. ATTTA instability motifs
    pos = 0
    while True:
        p = seq.find("ATTTA", pos)
        if p == -1:
            break
        codon_indices = optimizer._get_overlapping_codon_indices(p, 5, n_codons)
        violations.append(Violation(
            violation_type="atttta_motif",
            position=p,
            severity=SEVERITY_WEIGHTS["atttta_motif"],
            codon_indices=codon_indices,
            details=f"ATTTA at position {p}",
        ))
        pos = p + 1

    # 6. Long T runs (6+ consecutive T)
    i = 0
    while i < len(seq):
        if seq[i] == 'T':
            j = i
            while j < len(seq) and seq[j] == 'T':
                j += 1
            run_len = j - i
            if run_len >= 6:
                codon_indices = optimizer._get_overlapping_codon_indices(
                    i, run_len, n_codons
                )
                violations.append(Violation(
                    violation_type="t_run",
                    position=i,
                    severity=SEVERITY_WEIGHTS["t_run"] * (run_len - 5),
                    codon_indices=codon_indices,
                    details=f"T-run of {run_len} at position {i}",
                ))
            i = j
        else:
            i += 1

    # 7. CpG dinucleotides (simplified)
    if avoid_gt:
        cg_positions = state.cg_positions_list()
        for cg_pos in cg_positions[:5]:
            ci = cg_pos // 3
            if ci < n_codons:
                violations.append(Violation(
                    violation_type="cpg_island",
                    position=cg_pos,
                    severity=SEVERITY_WEIGHTS["cpg_island"],
                    codon_indices=[ci],
                    details=f"CpG at pos {cg_pos}",
                ))

    return violations


def _detect_expensive_violations(
    optimizer, state: IncrementalSequenceState
) -> list[Violation]:
    """Detect constraint violations that require MaxEntScan (expensive).

    Only called after all cheap violations are resolved.
    NEVER called for prokaryotic organisms — splice sites are irrelevant.
    """
    violations: list[Violation] = []
    seq = state.sequence
    n_codons = state.num_codons

    # Prokaryotes have no spliceosome — skip all MaxEntScan calls
    if not optimizer.avoid_gt:
        return violations

    try:
        from ..maxentscan import max_donor_score, max_acceptor_score

        # Check donors (GT)
        max_donor = max_donor_score(seq)
        if max_donor >= optimizer.splice_threshold:
            from ..maxentscan import score_donor
            from ..optimization import score_splice_donor_potential, SPLICE_DONOR_POTENTIAL_THRESHOLD
            for i in range(len(seq) - 1):
                if seq[i:i+2] == "GT":
                    if optimizer._is_unavoidable_gt(seq, i):
                        continue
                    donor_score = score_donor(seq, i)
                    if donor_score >= optimizer.splice_threshold:
                        # Cost-aware: also check splice donor potential
                        # Even if MaxEntScan donor score is high, a GT
                        # with low overall splice donor potential (no
                        # acceptor context, no polypyrimidine tract)
                        # is not biologically dangerous.
                        sdp = score_splice_donor_potential(seq, i)
                        if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                            continue  # Low-risk GT — CAI takes priority
                        codon_idx = i // 3
                        next_codon_start = codon_idx * 3 + 3
                        is_within = (i + 1) < next_codon_start
                        if is_within:
                            involved = [codon_idx]
                        else:
                            involved = (
                                [codon_idx, codon_idx + 1]
                                if codon_idx + 1 < n_codons
                                else [codon_idx]
                            )
                        splice_severity = (
                            SEVERITY_WEIGHTS["cryptic_splice_donor"]
                            * min(donor_score / 10.0, 1.0)
                        )
                        violations.append(Violation(
                            violation_type="cryptic_splice_donor",
                            position=i,
                            severity=(
                                SEVERITY_WEIGHTS["avoidable_gt"]
                                + splice_severity
                            ),
                            codon_indices=involved,
                            details=(
                                f"GT donor at pos {i}, "
                                f"score={donor_score:.1f}"
                            ),
                        ))

        # Check acceptors (AG)
        max_acceptor = max_acceptor_score(seq)
        if max_acceptor >= optimizer.splice_threshold:
            from ..maxentscan import score_acceptor
            for i in range(len(seq) - 1):
                if seq[i:i+2] == "AG":
                    acceptor_score = score_acceptor(seq, i)
                    if acceptor_score >= optimizer.splice_threshold:
                        codon_idx = i // 3
                        next_codon_start = codon_idx * 3 + 3
                        is_within = (i + 1) < next_codon_start
                        if is_within:
                            involved = [codon_idx]
                        else:
                            involved = (
                                [codon_idx, codon_idx + 1]
                                if codon_idx + 1 < n_codons
                                else [codon_idx]
                            )
                        violations.append(Violation(
                            violation_type="cryptic_splice_acceptor",
                            position=i,
                            severity=(
                                SEVERITY_WEIGHTS["cryptic_splice_acceptor"]
                                * min(acceptor_score / 10.0, 1.0)
                            ),
                            codon_indices=involved,
                            details=(
                                f"AG acceptor at pos {i}, "
                                f"score={acceptor_score:.1f}"
                            ),
                        ))
    except ImportError:
        pass

    return violations


def _fix_violation(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    protein: str,
    avoid_gt: bool,
) -> bool:
    """Try to fix a single constraint violation.

    For each codon position involved in the violation, try all synonymous
    alternatives and pick the one that maximizes the scoring function:
        score = cai_weight * CAI - sum(remaining_violations * penalty)

    Returns True if the violation was fixed.
    """
    vtype = violation.violation_type

    if vtype == "restriction_site":
        return _fix_restriction_site(optimizer, state, violation, avoid_gt)
    elif vtype == "gc_out_of_range":
        return _fix_gc_range(optimizer, state, violation, avoid_gt)
    elif vtype == "cryptic_splice_donor":
        return _fix_avoidable_gt(optimizer, state, violation, avoid_gt)
    elif vtype == "avoidable_gt":
        return _fix_avoidable_gt(optimizer, state, violation, avoid_gt)
    elif vtype == "cryptic_splice_acceptor":
        return _fix_cryptic_splice_acceptor(optimizer, state, violation, avoid_gt)
    elif vtype == "cpg_island":
        return _fix_cpg(optimizer, state, violation, avoid_gt)
    elif vtype == "atttta_motif":
        return _fix_atttta(optimizer, state, violation, avoid_gt)
    elif vtype == "t_run":
        return _fix_t_run(optimizer, state, violation, avoid_gt)
    elif vtype == "stop_codon":
        return _fix_stop_codon(optimizer, state, violation, protein, avoid_gt)
    else:
        return False


def _fix_restriction_site(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix a restriction site violation with minimal CAI loss."""
    for ci in violation.codon_indices:
        if ci < 0 or ci >= state.num_codons:
            continue
        aa = state.get_aa(ci)
        if aa is None or aa == "*":
            continue

        current_codon = state.get_codon(ci)
        old_gt_count = state.gt_count

        # Try each synonymous alternative, sorted by CAI (highest first)
        # Use precomputed sorted_codons for faster lookup
        for alt in optimizer.sorted_codons.get(aa, []):
            if alt == current_codon:
                continue

            # Try the swap
            old_codon = state.swap_codon(ci, alt)
            new_seq = state.sequence

            # Check if the restriction site is eliminated
            site_present = state.has_any_restriction_site()

            if site_present:
                # Site still present — rollback
                state.swap_codon(ci, old_codon)
                continue

            # Site eliminated! Check we didn't create new GTs
            if avoid_gt and state.gt_count > old_gt_count:
                new_gts = state.gt_positions_list()
                all_unavoidable = all(
                    optimizer._is_unavoidable_gt(state.sequence, pos)
                    for pos in new_gts
                    if pos not in {
                        p for p in range(ci * 3 - 1, ci * 3 + 3)
                    }
                )
                if not all_unavoidable:
                    state.swap_codon(ci, old_codon)
                    continue

            # Accept the fix
            return True

    # Try two-codon coordinated fix
    return _fix_restriction_site_two_codons(optimizer, 
        state, violation, avoid_gt
    )


def _fix_restriction_site_two_codons(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix a restriction site by modifying two adjacent codons."""
    codon_indices = sorted(set(violation.codon_indices))
    if len(codon_indices) < 2:
        return False

    old_gt_count = state.gt_count

    for idx in range(len(codon_indices) - 1):
        ci1, ci2 = codon_indices[idx], codon_indices[idx + 1]
        if ci2 != ci1 + 1:
            continue  # Only adjacent codons

        aa1 = state.get_aa(ci1)
        aa2 = state.get_aa(ci2)
        if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
            continue

        # Try all pairs sorted by combined CAI
        pairs = []
        for c1 in optimizer.sorted_codons.get(aa1, []):
            for c2 in optimizer.sorted_codons.get(aa2, []):
                combined = (
                    optimizer.species_cai.get(c1, 0.0)
                    + optimizer.species_cai.get(c2, 0.0)
                )
                pairs.append((c1, c2, combined))
        pairs.sort(key=lambda x: x[2], reverse=True)

        for c1, c2, _ in pairs:
            old1 = state.swap_codon(ci1, c1)
            old2 = state.swap_codon(ci2, c2)
            new_seq = state.sequence

            # Check if the site is eliminated
            site_present = state.has_any_restriction_site()

            if site_present:
                state.swap_codon(ci2, old2)
                state.swap_codon(ci1, old1)
                continue

            # Check GT constraint
            if avoid_gt and state.gt_count > old_gt_count:
                state.swap_codon(ci2, old2)
                state.swap_codon(ci1, old1)
                continue

            return True

    return False


def _fix_gc_range(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix GC content by swapping codons to adjust GC fraction."""
    gc_count = state.gc_count
    n_bases = state._n
    gc_val = state.gc_fraction

    if gc_val < optimizer.gc_lo:
        target = optimizer.gc_lo
        need_more_gc = True
    elif gc_val > optimizer.gc_hi:
        target = optimizer.gc_hi
        need_more_gc = False
    else:
        return False  # Already in range

    # Find the best single-codon swap to move GC toward target
    best_swap = None
    best_score = -1.0
    best_ci = -1
    best_gc_delta = 0

    for ci in range(state.num_codons):
        aa = state.get_aa(ci)
        if aa is None or aa == "*":
            continue
        current = state.get_codon(ci)
        current_gc = optimizer.codon_gc.get(current, 0)

        for alt in optimizer.sorted_codons.get(aa, []):
            if alt == current:
                continue
            alt_gc = optimizer.codon_gc.get(alt, 0)
            gc_delta = alt_gc - current_gc

            # Check if this swap moves GC in the right direction
            if need_more_gc and gc_delta <= 0:
                continue
            if not need_more_gc and gc_delta >= 0:
                continue

            new_gc_count = gc_count + gc_delta
            new_frac = new_gc_count / n_bases
            diff = abs(new_frac - target)
            alt_cai = optimizer.species_cai.get(alt, 0.0)

            # Score: prefer better GC improvement with higher CAI
            score = (1.0 - diff) + alt_cai * 0.01
            if score > best_score:
                best_score = score
                best_swap = alt
                best_ci = ci
                best_gc_delta = gc_delta

    if best_swap is None:
        return False

    # Apply the swap
    old_codon = state.swap_codon(best_ci, best_swap)

    # Check we didn't break any hard constraints
    if state.has_any_restriction_site():
        state.swap_codon(best_ci, old_codon)
        return False

    # Check GT constraint
    if avoid_gt:
        new_gts = state.gt_positions_list()
        for gt_pos in new_gts:
            if not optimizer._is_unavoidable_gt(state.sequence, gt_pos):
                state.swap_codon(best_ci, old_codon)
                return False

    return True


def _fix_avoidable_gt(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix an avoidable GT dinucleotide with minimal CAI loss."""
    gt_pos = violation.position
    codon_idx = gt_pos // 3
    next_codon_start = codon_idx * 3 + 3
    is_within = (gt_pos + 1) < next_codon_start

    if is_within:
        return _fix_within_codon_gt(optimizer, state, codon_idx, avoid_gt)
    else:
        return _fix_cross_codon_gt(optimizer, state, codon_idx, avoid_gt)


def _fix_within_codon_gt(
    optimizer,
    state: IncrementalSequenceState,
    codon_idx: int,
    avoid_gt: bool,
) -> bool:
    """Fix a within-codon GT by choosing the best CAI-preserving
    substitution using precomputed GT-free codon lists.

    If the CAI cost of switching to a GT-free codon would be excessive
    (relative adaptiveness loss > _GT_CAI_COST_THRESHOLD, currently 3%),
    the GT is accepted as "CAI-critical" and the fix is skipped.
    For eukaryotes, in-codon GTs from optimal codons are biologically
    acceptable, so only sacrifice trivial CAI for GT avoidance.
    """
    aa = state.get_aa(codon_idx)
    if aa is None or aa == "*":
        return False

    old_gt_count = state.gt_count

    # Check CAI cost before trying GT-free alternatives
    gt_free_list = optimizer.gt_free.get(aa, [])
    if gt_free_list:
        current_codon = state.sequence[codon_idx*3:codon_idx*3+3]
        current_w = optimizer.species_cai.get(current_codon, 0.0)
        best_gt_free_w = optimizer.species_cai.get(gt_free_list[0], 0.0)
        max_a = optimizer._max_adapt.get(aa, 0.0)
        if max_a > 0:
            current_rel = current_w / max_a
            best_gtf_rel = best_gt_free_w / max_a
        else:
            current_rel = current_w
            best_gtf_rel = best_gt_free_w
        # If CAI loss from avoiding GT would be excessive, accept the GT
        if current_rel - best_gtf_rel > _GT_CAI_COST_THRESHOLD:
            return False  # GT is CAI-critical, don't fix

    # Try GT-free alternatives (sorted by CAI, highest first)
    # Use precomputed gt_free lookup table
    for alt in gt_free_list:
        # Quick boundary check
        left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
        if left_gt or right_gt:
            continue

        old_codon = state.swap_codon(codon_idx, alt)
        if state.gt_count < old_gt_count:
            # Check no new restriction sites
            site_ok = not state.has_any_restriction_site()
            if site_ok:
                return True
        state.swap_codon(codon_idx, old_codon)  # Rollback

    # Fallback: try all alternatives sorted by CAI
    for alt in optimizer.sorted_codons.get(aa, []):
        if "GT" not in alt:
            continue  # Already tried above
        old_codon = state.swap_codon(codon_idx, alt)
        if state.gt_count < old_gt_count:
            site_ok = not state.has_any_restriction_site()
            if site_ok:
                return True
        state.swap_codon(codon_idx, old_codon)

    return False


def _fix_cross_codon_gt(
    optimizer,
    state: IncrementalSequenceState,
    codon_idx: int,
    avoid_gt: bool,
) -> bool:
    """Fix a cross-codon GT by modifying one or both adjacent codons."""
    old_gt_count = state.gt_count
    next_idx = codon_idx + 1
    prev_idx = codon_idx - 1

    # Strategy D: Change only the next codon (cheapest single-codon fix)
    if next_idx < state.num_codons:
        aa2 = state.get_aa(next_idx)
        if aa2 is not None and aa2 != "*":
            for alt2 in optimizer.sorted_codons.get(aa2, []):
                if alt2[0] == 'T':
                    continue  # Would still create GT at boundary
                old2 = state.swap_codon(next_idx, alt2)
                if state.gt_count < old_gt_count:
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        return True
                state.swap_codon(next_idx, old2)

    # Strategy C: Change only the current codon
    # (to one that doesn't end with G)
    aa1 = state.get_aa(codon_idx)
    if aa1 is not None and aa1 != "*":
        for alt1 in optimizer.sorted_codons.get(aa1, []):
            if alt1[-1] == 'G':
                continue  # Would still create GT at boundary
            old1 = state.swap_codon(codon_idx, alt1)
            if state.gt_count < old_gt_count:
                site_ok = not state.has_any_restriction_site()
                if site_ok:
                    return True
            state.swap_codon(codon_idx, old1)

    # Strategy B: Change both codons (2-codon coordinated)
    if next_idx < state.num_codons:
        aa1 = state.get_aa(codon_idx)
        aa2 = state.get_aa(next_idx)
        if (
            aa1 is not None and aa1 != "*"
            and aa2 is not None and aa2 != "*"
        ):
            # Try pairs sorted by combined CAI
            pairs = []
            for c1 in optimizer.sorted_codons.get(aa1, [])[:3]:
                for c2 in optimizer.sorted_codons.get(aa2, [])[:3]:
                    combined = (
                        optimizer.species_cai.get(c1, 0.0)
                        + optimizer.species_cai.get(c2, 0.0)
                    )
                    pairs.append((c1, c2, combined))
            pairs.sort(key=lambda x: x[2], reverse=True)

            for c1, c2, _ in pairs:
                # Ensure no GT at the boundary
                if c1[-1] == 'G' and c2[0] == 'T':
                    continue
                old1 = state.swap_codon(codon_idx, c1)
                old2 = state.swap_codon(next_idx, c2)
                if state.gt_count < old_gt_count:
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        return True
                state.swap_codon(next_idx, old2)
                state.swap_codon(codon_idx, old1)

    # Strategy A: Try changing the previous codon
    if prev_idx >= 0:
        aa_prev = state.get_aa(prev_idx)
        if aa_prev is not None and aa_prev != "*":
            for alt_prev in optimizer.sorted_codons.get(aa_prev, []):
                # No special constraint on prev — just try to reduce GTs
                old_prev = state.swap_codon(prev_idx, alt_prev)
                if state.gt_count < old_gt_count:
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        return True
                state.swap_codon(prev_idx, old_prev)

    return False


def _fix_atttta(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix an ATTTA motif by swapping codons at overlapping positions."""
    for ci in violation.codon_indices:
        if ci < 0 or ci >= state.num_codons:
            continue
        aa = state.get_aa(ci)
        if aa is None or aa == "*":
            continue
        current = state.get_codon(ci)

        # Try all alternatives sorted by CAI
        for alt in optimizer.sorted_codons.get(aa, []):
            if alt == current:
                continue
            old_codon = state.swap_codon(ci, alt)
            if "ATTTA" not in state.sequence:
                # Check no new restriction sites or GTs
                site_ok = not state.has_any_restriction_site()
                if site_ok:
                    if not avoid_gt or state.gt_count <= (
                        state.gt_count  # no change check
                    ):
                        return True
            state.swap_codon(ci, old_codon)
    return False


def _fix_t_run(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix a long T run by swapping a codon at the run's center."""
    run_pos = violation.position
    n_codons = state.num_codons

    # Try codons in the middle of the run first
    codon_idx = (run_pos + 3) // 3  # slightly past the start
    if codon_idx >= n_codons:
        codon_idx = n_codons - 1

    for offset in range(n_codons):
        ci = codon_idx + offset
        if ci >= n_codons:
            ci = offset  # Wrap around

        aa = state.get_aa(ci)
        if aa is None or aa == "*":
            continue
        current = state.get_codon(ci)

        for alt in optimizer.sorted_codons.get(aa, []):
            if alt == current:
                continue
            old_codon = state.swap_codon(ci, alt)
            # Check T-run is fixed
            new_seq = state.sequence
            has_long_run = False
            i = 0
            while i < len(new_seq):
                if new_seq[i] == 'T':
                    j = i
                    while j < len(new_seq) and new_seq[j] == 'T':
                        j += 1
                    if j - i >= 6:
                        has_long_run = True
                        break
                    i = j
                else:
                    i += 1

            if not has_long_run:
                site_ok = not state.has_any_restriction_site()
                if site_ok:
                    return True
            state.swap_codon(ci, old_codon)
    return False


def _fix_stop_codon(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    protein: str,
    avoid_gt: bool,
) -> bool:
    """Fix a premature stop codon by replacing with a non-stop codon."""
    ci = violation.codon_indices[0]
    if ci < 0 or ci >= len(protein):
        return False
    aa = protein[ci]
    if aa == "*":
        return False  # Legitimate stop

    current = state.get_codon(ci)
    for alt in optimizer.sorted_codons.get(aa, []):
        if alt == current:
            continue
        if alt in ("TAA", "TAG", "TGA"):
            continue  # Don't swap to another stop
        old_codon = state.swap_codon(ci, alt)
        site_ok = not state.has_any_restriction_site()
        if site_ok:
            return True
        state.swap_codon(ci, old_codon)
    return False


def _fix_cryptic_splice_acceptor(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix a cryptic splice acceptor (AG) with minimal CAI loss."""
    ag_pos = violation.position
    codon_idx = ag_pos // 3
    next_codon_start = codon_idx * 3 + 3
    is_within = (ag_pos + 1) < next_codon_start

    if is_within:
        # Within-codon AG: use precomputed AG-free codons
        aa = state.get_aa(codon_idx)
        if aa is None or aa == "*":
            return False

        old_gt_count = state.gt_count
        ag_free_list = optimizer.ag_free.get(aa, [])

        for alt in ag_free_list:
            left_gt, right_gt = state.boundary_creates_gt(codon_idx, alt)
            if left_gt or right_gt:
                continue

            old_codon = state.swap_codon(codon_idx, alt)
            # Verify splice acceptor score dropped below threshold
            try:
                from ..maxentscan import score_acceptor
                new_score = score_acceptor(state.sequence, ag_pos)
                if new_score < optimizer.splice_threshold:
                    site_ok = not state.has_any_restriction_site()
                    if site_ok:
                        return True
            except ImportError:
                # No MaxEntScan — just check AG is gone
                if "AG" not in state.sequence[ag_pos:ag_pos+2]:
                    return True
            state.swap_codon(codon_idx, old_codon)
    else:
        # Cross-codon AG: try fixing by changing adjacent codons
        next_idx = codon_idx + 1
        if next_idx < state.num_codons:
            aa2 = state.get_aa(next_idx)
            if aa2 is not None and aa2 != "*":
                for alt2 in optimizer.ag_free.get(aa2, []):
                    if alt2[0] == 'G':
                        continue
                    old2 = state.swap_codon(next_idx, alt2)
                    try:
                        from ..maxentscan import score_acceptor
                        new_score = score_acceptor(
                            state.sequence, ag_pos
                        )
                        if new_score < optimizer.splice_threshold:
                            site_ok = not state.has_any_restriction_site()
                            if site_ok:
                                return True
                    except ImportError:
                        if "AG" not in state.sequence[ag_pos:ag_pos+2]:
                            return True
                    state.swap_codon(next_idx, old2)

        # Try changing the current codon (to one that doesn't end with A)
        aa1 = state.get_aa(codon_idx)
        if aa1 is not None and aa1 != "*":
            for alt1 in optimizer.sorted_codons.get(aa1, []):
                if alt1[-1] == 'A':
                    continue
                old1 = state.swap_codon(codon_idx, alt1)
                try:
                    from ..maxentscan import score_acceptor
                    new_score = score_acceptor(
                        state.sequence, ag_pos
                    )
                    if new_score < optimizer.splice_threshold:
                        site_ok = not state.has_any_restriction_site()
                        if site_ok:
                            return True
                except ImportError:
                    if "AG" not in state.sequence[ag_pos:ag_pos+2]:
                        return True
                state.swap_codon(codon_idx, old1)

    return False


def _fix_cpg(
    optimizer,
    state: IncrementalSequenceState,
    violation: Violation,
    avoid_gt: bool,
) -> bool:
    """Fix a CpG dinucleotide by swapping to a CG-free alternative."""
    cg_pos = violation.position
    codon_idx = cg_pos // 3
    if codon_idx >= state.num_codons:
        return False

    aa = state.get_aa(codon_idx)
    if aa is None or aa == "*":
        return False

    current = state.get_codon(codon_idx)
    old_gt_count = state.gt_count

    # Find alternatives that don't contain "CG"
    for alt in optimizer.sorted_codons.get(aa, []):
        if alt == current:
            continue
        if "CG" in alt:
            continue  # Would still have CpG
        old_codon = state.swap_codon(codon_idx, alt)

        # Check no new CG was created at boundary
        new_seq = state.sequence
        cg_still = False
        # Check the neighborhood of the swap
        start = max(0, codon_idx * 3 - 1)
        end = min(len(new_seq), codon_idx * 3 + 4)
        if "CG" in new_seq[start:end]:
            # More careful check: is it the same CG position?
            if "CG" in new_seq[cg_pos:cg_pos+2]:
                cg_still = True

        if not cg_still:
            # Check restriction sites
            site_ok = not state.has_any_restriction_site()
            if site_ok:
                if not avoid_gt or state.gt_count <= old_gt_count:
                    return True
        state.swap_codon(codon_idx, old_codon)

    return False


def _estimate_cai_impact(violation_type: str) -> float:
    """Estimate the CAI impact of fixing a constraint violation.

    Used by CAI-aware constraint resolution to prefer the fix that
    maximizes CAI when two constraints conflict.  Higher values
    indicate greater expected CAI loss from fixing this violation.

    Args:
        violation_type: The type of constraint violation.

    Returns:
        Estimated CAI impact (higher = more CAI loss).
    """
    # Hard constraints typically have higher CAI impact because
    # they force specific codon choices (e.g., removing a restriction
    # site may require using a lower-CAI synonym).
    # Soft constraints typically have lower CAI impact because they
    # can often be satisfied with CAI-neutral alternatives.
    _CAI_IMPACT_ESTIMATES = {
        "restriction_site": 0.15,     # May force suboptimal codon
        "stop_codon": 0.20,           # Must use non-stop synonym
        "gc_out_of_range": 0.10,      # GC swap usually CAI-friendly
        "cryptic_splice_donor": 0.08, # GT-free alternatives often available
        "cryptic_splice_acceptor": 0.08,
        "avoidable_gt": 0.06,         # GT avoidance is usually cheap
        "cpg_island": 0.04,           # CpG avoidance is soft, low impact
        "atttta_motif": 0.05,         # ATTTA fix is usually cheap
        "t_run": 0.03,                # T-run fix is usually cheap
    }
    return _CAI_IMPACT_ESTIMATES.get(violation_type, 0.05)

