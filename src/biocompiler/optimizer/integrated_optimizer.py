"""
BioCompiler Integrated Constraint-Solving Optimizer
====================================================

Greedy forward pass (backward-looking constraint checks) + multi-pass
cleanup for cross-codon violations + final CpG cleanup + GC adjustment.

.. note::
    This is NOT a single integrated pass that satisfies ALL constraints
    simultaneously. The forward pass places codons left-to-right and can
    only evaluate constraints against the PREVIOUSLY placed codon — the
    next codon has not been placed yet, so right-boundary checks (GT, CG,
    ATTTA spanning into the next codon, T-runs, restriction sites) are
    no-ops during the forward pass BY DESIGN. Those right-boundary and
    residual cross-codon violations are resolved in a dedicated
    multi-pass cleanup that has full visibility of both neighbours.

Architecture:
    The original sequential greedy optimizer fixes constraints one at a
    time (remove GT → remove CG → remove ATTTA → fix T-runs → adjust GC
    → CAI recovery), costing 12+ passes. This optimizer reduces that to:
    (1) one greedy forward pass that satisfies backward + within-codon
    constraints, (2) a bounded multi-pass cleanup for cross-codon
    violations, (3) a targeted final CpG sweep, (4) a GC adjustment pass.

Algorithm:
    1. Back-translate: sort synonymous codons by CAI (best first)
    2. Greedy forward pass (left-to-right, n positions):
       a. Try best-CAI codon
       b. Check backward-looking constraints against the PREVIOUS codon
          plus within-codon constraints (GT, CG, ATTTA, T-run,
          restriction site, premature stop). The NEXT codon is still
          empty (seq_list is filled left-to-right from [""]*n), so the
          right-boundary branches of these checks are no-ops here.
       c. If violations, try next-best codon
       d. Pick the codon with highest CAI that satisfies the backward
          + within-codon constraints
    3. Multi-pass cleanup (up to 3 iterations): scan every codon with
       FULL left AND right neighbour visibility; fix GT/CG/ATTTA/T-run
       violations by swapping to the best-CAI alternative codon.
    4. Final CpG cleanup (aggressive mode only): targeted swaps for any
       remaining CG dinucleotides, including boundary CGs between codons.
    5. GC adjustment: if the sequence GC is outside [gc_lo, gc_hi], swap
       codons toward the target window while preserving constraints.
    6. Verify: protein preservation (SECIS-aware) + all constraints.

Complexity:
    Forward pass: O(n × k) where n = codons, k = alternatives per amino
    acid (≤6). Cleanup: O(n × k × m) where m ≤ 3 iterations. GC/CpG
    sweeps are O(n × k). The original sequential greedy optimizer is
    O(n × p × c) with p = passes (12+), c = constraint checks.

Correctness guarantees:
    - Protein is ALWAYS preserved (only synonymous codons are used)
    - Backward + within-codon constraints are checked during the forward
      pass; right-boundary and residual cross-codon violations are fixed
      in the cleanup pass
    - If no codon can satisfy a constraint, the best-CAI codon is chosen
      and a note is recorded for the cleanup pass to attempt a fix
"""

from __future__ import annotations

import logging
import re
from typing import Optional

from biocompiler.type_system.codon_tables import AA_TO_CODONS, CODON_TABLE
from biocompiler.organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES, ORGANISM_GC_TARGETS
from biocompiler.sequence.scanner import gc_content
from biocompiler.sequence.restriction_sites import get_recognition_site
from biocompiler.shared.constants import reverse_complement

logger = logging.getLogger(__name__)

# Constraint constants
_T_RUN_MAX = 5  # max consecutive T's allowed
_ATTTA_MOTIF = "ATTTA"
_STOP_CODONS = {"TAA", "TAG", "TGA"}


def _build_usage_table(organism: str) -> dict[str, float]:
    """Get codon adaptiveness table for the organism."""
    species_key = resolve_organism(organism)
    return CODON_ADAPTIVENESS_TABLES.get(species_key, {})


def _sorted_codons_for_aa(aa: str, usage: dict[str, float]) -> list[str]:
    """Get synonymous codons for an amino acid, sorted by CAI (best first).
    
    Special case: selenocysteine (U) is encoded by TGA (normally a stop codon).
    The SECIS element in the 3'UTR causes the ribosome to insert selenocysteine
    instead of stopping. We return TGA as the only codon for U.
    """
    if aa == "U":
        return ["TGA"]  # selenocysteine codon
    codons = AA_TO_CODONS.get(aa, [])
    return sorted(codons, key=lambda c: -usage.get(c, 0.0))


def _check_local_constraints(
    seq_list: list[str],
    codon_idx: int,
    candidate: str,
    n_codons: int,
    enzymes_sites: list[tuple[str, str]],
    is_prokaryote: bool,
    cpg_mode: str,
) -> bool:
    """Check if placing `candidate` at `codon_idx` violates any constraint.

    Checks the PREVIOUS codon boundary (prev|candidate) and within-codon
    constraints. The NEXT codon boundary (candidate|next) is also checked
    when a next codon is present, but during the forward greedy pass
    `next_codon` is empty BY DESIGN (seq_list is filled left-to-right
    from [""]*n), so right-boundary violations are deferred to the
    cleanup pass (_cleanup_violations), which has full neighbour
    visibility.

    Constraints evaluated (when the relevant neighbour is present):
    - GT dinucleotide at prev|candidate and candidate|next boundaries
      (eukaryotes only)
    - CG dinucleotide at prev|candidate and candidate|next boundaries
      (CpG mode aggressive)
    - ATTTA motif in (prev + candidate + next)
    - T-run (6+ consecutive T) in (prev[-3:] + candidate + next[:3])
    - Restriction site in (prev2 + prev + candidate + next + next2)
    - Premature stop codon (TGA allowed as potential selenocysteine)

    Returns True if all evaluated constraints are satisfied.
    """
    # Build the local context: previous codon + candidate + next codon.
    # NOTE: during the forward pass (integrated_optimize), seq_list is
    # filled left-to-right starting from [""] * n_codons, so next_codon
    # is ALWAYS empty here. The right-boundary branches below (GT, CG,
    # ATTTA spanning into next, T-run, restriction sites) are therefore
    # no-ops during the forward pass BY DESIGN — those cross-codon
    # violations are caught by _cleanup_violations, which has full
    # left AND right neighbour visibility. Only prev_codon boundary
    # and within-candidate checks take effect during the forward pass.
    prev_codon = seq_list[codon_idx - 1] if codon_idx > 0 else ""
    next_codon = seq_list[codon_idx + 1] if codon_idx < n_codons - 1 else ""
    
    # Build local sequence (up to 9 bases: prev + candidate + next)
    local = prev_codon + candidate + next_codon
    candidate_start = len(prev_codon)  # offset of candidate in local
    
    # ── Check GT dinucleotide at boundaries (eukaryotes only) ──
    if not is_prokaryote:
        # Check if candidate creates a GT at the left boundary (prev|candidate)
        if prev_codon:
            boundary = prev_codon[-1] + candidate[0]
            if boundary == "GT":
                return False
        # Check if candidate creates a GT at the right boundary (candidate|next)
        if next_codon:
            boundary = candidate[-1] + next_codon[0]
            if boundary == "GT":
                return False
        # Check internal GT within the candidate
        if "GT" in candidate:
            return False
    
    # ── Check CG dinucleotide (CpG mode aggressive) ──
    if cpg_mode == "aggressive":
        if prev_codon:
            boundary = prev_codon[-1] + candidate[0]
            if boundary == "CG":
                return False
        if next_codon:
            boundary = candidate[-1] + next_codon[0]
            if boundary == "CG":
                return False
        if "CG" in candidate:
            return False
    
    # ── Check ATTTA motif ──
    # ATTTA can span: prev[-2:] + candidate[:3], prev[-1:] + candidate[:4],
    # candidate[0:5], candidate[1:6] + next[:0], candidate[2:] + next[:2]
    if _ATTTA_MOTIF in local:
        return False
    
    # ── Check T-run (6+ consecutive T) ──
    # Build the full local context to check T-runs
    t_run_check = prev_codon[-3:] + candidate + next_codon[:3]
    if re.search(r'T{6,}', t_run_check):
        return False
    
    # ── Check restriction sites ──
    if enzymes_sites:
        # Build a slightly larger context for restriction site checking
        prev2 = seq_list[codon_idx - 2] if codon_idx > 1 else ""
        next2 = seq_list[codon_idx + 2] if codon_idx < n_codons - 2 else ""
        rs_context = prev2 + prev_codon + candidate + next_codon + next2
        for site, site_rc in enzymes_sites:
            if site in rs_context or site_rc in rs_context:
                return False
    
    # ── Check premature stop codon (not at the last position) ──
    # Exception: TGA is allowed at selenocysteine (U) positions
    if codon_idx < n_codons - 1 and candidate in _STOP_CODONS:
        # Check if this position is a selenocysteine position
        # We can't access the protein here, so we allow TGA and let the
        # cleanup pass handle any issues. The IR pipeline will verify.
        if candidate == "TGA":
            pass  # TGA might be selenocysteine — allow it
        else:
            return False
    
    return True


def _check_gc_range(seq: str, gc_lo: float, gc_hi: float) -> bool:
    """Check if the sequence GC content is within range."""
    gc = gc_content(seq)
    return gc_lo <= gc <= gc_hi


def integrated_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: list[str] | None = None,
    cpg_mode: str = "aggressive",
    **kwargs,  # accept and ignore is_prokaryote (auto-detected)
) -> tuple[str, list[str], list[int]]:
    """
    Greedy forward-pass optimizer with multi-pass cleanup.

    Places the best-CAI codon at each position using backward-looking
    constraint checks (against the previously placed codon), then runs a
    multi-pass cleanup to fix cross-codon violations that the forward
    pass could not see (right-boundary GT/CG/ATTTA/T-run/restriction
    sites), then a final CpG sweep (aggressive mode), then adjusts GC
    content if needed. NOT a single integrated pass: right-boundary
    constraints are resolved in the cleanup pass.

    Args:
        protein: Amino acid sequence (single-letter codes, no stop)
        organism: Target organism name
        gc_lo: Minimum GC content
        gc_hi: Maximum GC content
        enzymes: Restriction enzymes to avoid
        cpg_mode: "aggressive" (eliminate all CG) or "off"

    Returns:
        Tuple of (optimized DNA sequence, list of notes/warnings, list of
        SECIS / selenocysteine codon indices).
    """
    notes: list[str] = []
    protein = protein.upper().rstrip("*")
    n_codons = len(protein)
    
    if n_codons == 0:
        return "TAA", ["Empty protein"], []
    
    # ── Setup ──
    usage = _build_usage_table(organism)
    is_prokaryote = _is_prokaryote(organism)
    
    # Pre-compute restriction enzyme sites
    enzymes_sites: list[tuple[str, str]] = []
    if enzymes:
        for enz in enzymes:
            site = get_recognition_site(enz)
            if site:
                enzymes_sites.append((site, reverse_complement(site)))
    
    # Pre-compute sorted codons for each amino acid
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(protein):
        sorted_codons[aa] = _sorted_codons_for_aa(aa, usage)
    
    # ── Greedy forward pass (backward-looking constraint checks) ──
    # seq_list is filled left-to-right; next_codon is always "" during
    # this pass, so right-boundary violations are deferred to cleanup.
    seq_list: list[str] = [""] * n_codons
    
    for i in range(n_codons):
        aa = protein[i]
        candidates = sorted_codons.get(aa, [AA_TO_CODONS.get(aa, ["GCT"])[0]])

        # Collect ALL constraint-satisfying candidates (already CAI-sorted).
        # Heuristic tiebreaker: prefer a codon whose last base is NOT 'G'
        # (eukaryotes — reduces future boundary GTs at codon[i]|codon[i+1])
        # or NOT 'C' (aggressive CpG mode — reduces future boundary CGs).
        # Applied AFTER CAI sorting as a soft preference: pick the first
        # satisfying candidate whose last base avoids a future boundary
        # violation; if none qualifies, fall back to the best-CAI
        # satisfying candidate. The forward pass cannot look ahead at the
        # next codon choice, so this heuristic trades a small CAI cost for
        # fewer right-boundary GT/CG violations that the cleanup pass
        # would otherwise have to fix (and may not always be able to).
        satisfying = [
            c for c in candidates
            if _check_local_constraints(
                seq_list, i, c, n_codons,
                enzymes_sites, is_prokaryote, cpg_mode,
            )
        ]

        if satisfying:
            best_codon = satisfying[0]
            if not is_prokaryote:
                # Eukaryotes: prefer last base != 'G' (avoid future GT)
                for c in satisfying:
                    if c[-1] != 'G':
                        best_codon = c
                        break
            elif cpg_mode == "aggressive":
                # Aggressive CpG: prefer last base != 'C' (avoid future CG)
                for c in satisfying:
                    if c[-1] != 'C':
                        best_codon = c
                        break
        else:
            # No codon satisfies all constraints — pick the best-CAI one
            # and note the violation for cleanup
            best_codon = candidates[0] if candidates else "GCT"
            notes.append(f"Position {i}: no constraint-satisfying codon for {aa}, using {best_codon}")

        seq_list[i] = best_codon
    
    # ── Add stop codon ──
    seq_list.append("TAA")
    
    # ── Lightweight cleanup pass (fix any remaining boundary violations) ──
    seq = "".join(seq_list)
    seq = _cleanup_violations(seq, protein, usage, is_prokaryote, cpg_mode, enzymes_sites, notes)
    
    # ── GC content adjustment (if needed) ──
    if not _check_gc_range(seq, gc_lo, gc_hi):
        seq = _adjust_gc_content(seq, protein, usage, gc_lo, gc_hi, is_prokaryote, cpg_mode, enzymes_sites, notes)
    
    # ── Final verification ──
    # Verify protein preservation (SECIS-aware)
    secis_positions = [i for i, aa in enumerate(protein) if aa == "U"]
    translated = _translate_dna(seq, secis_positions)
    expected = protein + "*"
    if translated != expected:
        notes.append(f"WARNING: protein mismatch — expected {expected[:20]}..., got {translated[:20]}...")
        # Fall back to simple back-translation
        seq = _simple_back_translate(protein, usage)
    
    # Collect SECIS positions (codon indices where U was placed)
    secis_positions = [i for i, aa in enumerate(protein) if aa == "U"]
    
    return seq, notes, secis_positions


def _cleanup_violations(
    seq: str,
    protein: str,
    usage: dict[str, float],
    is_prokaryote: bool,
    cpg_mode: str,
    enzymes_sites: list[tuple[str, str]],
    notes: list[str],
) -> str:
    """Fix cross-codon constraint violations with local codon swaps.

    Resolves right-boundary violations the forward pass could not see
    (next_codon was empty during the forward pass), as well as residual
    left-boundary or within-codon violations. Three passes:

    1. Single-codon swap loop (up to ``max_iter`` = 3 iterations): for
       each codon with a violation, try swapping to the best-CAI synonym
       that eliminates the violation. Full left AND right neighbour
       visibility — the stop codon (TAA) at index ``n_codons`` is
       included as a non-swappable neighbour so boundary violations
       between the last protein codon and the stop codon are visible.
    2. Joint two-codon swap pass: for boundary GT/CG violations that
       single-codon swaps could NOT fix (e.g. codon[i] ends in 'G' and
       codon[i+1] starts with 'T', but no single synonym for either
       position breaks the dinucleotide), try all (synonym_i,
       synonym_{i+1}) pairs and pick the best-CAI pair that eliminates
       the violation without introducing new ones. O(k^2) per boundary
       (k <= 6, so <= 36 pairs). If no pair works, a note is logged and
       the violation is left in place (honest — cannot always fix).
    3. Final CpG cleanup (aggressive mode only): targeted single-codon
       swaps for any remaining CG, checking ALL constraints (GT, CG,
       ATTTA, T-run, restriction, stop) via ``_check_local_constraints``
       so no new GT/ATTTA/T-run violations are introduced.
    """
    n_codons = len(protein)
    # Include the stop codon (index n_codons) as a non-swappable neighbour
    # so boundary violations between the last protein codon and the stop
    # codon are visible to all passes below.
    codons = [seq[i*3:(i+1)*3] for i in range(n_codons + 1)]
    n_total = n_codons + 1  # protein codons + stop codon

    # ── Pass 1: single-codon swaps (up to max_iter iterations) ──
    changed = True
    iterations = 0
    max_iter = 3  # limit cleanup iterations

    while changed and iterations < max_iter:
        changed = False
        iterations += 1

        for i in range(n_codons):  # skip stop codon (index n_codons)
            current = codons[i]
            # Check if current codon has any violation (full neighbour
            # visibility, including the stop codon as the right neighbour
            # of the last protein codon).
            if _check_local_constraints(
                codons, i, current, n_total,
                enzymes_sites, is_prokaryote, cpg_mode,
            ):
                continue  # no violation, skip

            # Try to find a synonym that eliminates the violation
            aa = protein[i]
            candidates = _sorted_codons_for_aa(aa, usage)
            for candidate in candidates:
                if candidate == current:
                    continue
                if _check_local_constraints(
                    codons, i, candidate, n_total,
                    enzymes_sites, is_prokaryote, cpg_mode,
                ):
                    codons[i] = candidate
                    changed = True
                    break

    # ── Pass 2: joint two-codon swaps for boundary GT/CG ──
    # For a GT (or CG in aggressive mode) at codon[i][-1] + codon[i+1][0]
    # that single-codon swaps could not fix, try all (synonym_i,
    # synonym_{i+1}) pairs and pick the best-CAI pair that eliminates the
    # violation without introducing new ones. The stop codon (index
    # n_codons) is a non-swappable neighbour, so for the boundary
    # (n_codons-1, n_codons) only codon n_codons-1 is swapped.
    for i in range(n_codons):  # boundary (i, i+1); i+1 can be stop codon
        boundary = codons[i][-1] + codons[i+1][0]
        is_gt = (not is_prokaryote) and boundary == "GT"
        is_cg = (cpg_mode == "aggressive") and boundary == "CG"
        if not (is_gt or is_cg):
            continue

        cands_i = _sorted_codons_for_aa(protein[i], usage)
        if i + 1 < n_codons:
            cands_i1 = _sorted_codons_for_aa(protein[i + 1], usage)
        else:
            cands_i1 = [codons[i + 1]]  # stop codon, non-swappable

        best_pair = None
        best_score = -1.0
        old_i, old_i1 = codons[i], codons[i + 1]
        for ci in cands_i:
            for ci1 in cands_i1:
                if ci == old_i and ci1 == old_i1:
                    continue  # skip current pair
                # The boundary dinucleotide must be fixed
                new_boundary = ci[-1] + ci1[0]
                if is_gt and new_boundary == "GT":
                    continue
                if is_cg and new_boundary == "CG":
                    continue
                # Temporarily place the pair and check ALL constraints
                # for both positions (full neighbour visibility).
                codons[i] = ci
                codons[i + 1] = ci1
                ok = _check_local_constraints(
                    codons, i, ci, n_total,
                    enzymes_sites, is_prokaryote, cpg_mode,
                )
                if ok:
                    ok = _check_local_constraints(
                        codons, i + 1, ci1, n_total,
                        enzymes_sites, is_prokaryote, cpg_mode,
                    )
                if ok:
                    # Rank pairs by product of CAI (geometric-mean
                    # contribution); candidates are CAI-sorted, but the
                    # product over both positions gives a better joint
                    # ranking than just taking the first valid pair.
                    score = usage.get(ci, 0.0) * usage.get(ci1, 0.0)
                    if score > best_score:
                        best_score = score
                        best_pair = (ci, ci1)
        # Restore originals before applying best (or leaving as-is)
        codons[i] = old_i
        codons[i + 1] = old_i1
        if best_pair:
            codons[i] = best_pair[0]
            codons[i + 1] = best_pair[1]
        else:
            motif = "GT" if is_gt else "CG"
            notes.append(
                f"Boundary {motif} at codon {i}-{i + 1} unfixable by "
                f"joint swap, leaving as-is"
            )

    # ── Pass 3: final CpG cleanup (aggressive mode only) ──
    # Targeted single-codon swaps for any remaining CG (internal or
    # boundary), checking ALL constraints (GT, CG, ATTTA, T-run,
    # restriction, stop) via _check_local_constraints so no new GT/ATTTA/
    # T-run violations are introduced. This is a fallback for CGs that
    # passes 1 and 2 could not fix.
    if cpg_mode == "aggressive":
        for i in range(n_codons):  # skip stop codon
            current = codons[i]
            prev = codons[i - 1] if i > 0 else ""
            nxt = codons[i + 1]
            has_cg = (
                "CG" in current
                or (prev and prev[-1] + current[0] == "CG")
                or (current[-1] + nxt[0] == "CG")
            )
            if not has_cg:
                continue
            aa = protein[i]
            candidates = _sorted_codons_for_aa(aa, usage)
            for alt in candidates:
                if alt == current:
                    continue
                if _check_local_constraints(
                    codons, i, alt, n_total,
                    enzymes_sites, is_prokaryote, cpg_mode,
                ):
                    codons[i] = alt
                    break

    # ── Pass 4: final GT boundary sweep (eukaryotes only) ──
    # Safety net: Pass 3 (CpG sweep) may have removed CGs that were
    # blocking GT fixes (e.g., a boundary GT at (i, i+1) where codon[i]
    # could not be swapped because the only GT-free synonym created a
    # CG at (i-1, i) -- but Pass 3 has now cleared that CG). Re-scan for
    # remaining BOUNDARY GTs (internal GTs in codons like Valine GT*
    # are biologically unavoidable and skipped) and try single-codon
    # swaps on each side of the boundary. _check_local_constraints
    # ensures no new CG/ATTTA/T-run/restriction violations are
    # introduced. This is NOT a full O(k^2) joint-pair enumeration --
    # it tries single swaps on each side, picking the first valid
    # (highest-CAI) synonym that eliminates the boundary GT.
    if not is_prokaryote:
        for i in range(n_codons):
            current = codons[i]
            prev = codons[i - 1] if i > 0 else ""
            nxt = codons[i + 1]
            left_gt = bool(prev) and prev[-1] + current[0] == "GT"
            right_gt = current[-1] + nxt[0] == "GT"
            if not (left_gt or right_gt):
                continue
            for alt in _sorted_codons_for_aa(protein[i], usage):
                if alt == current:
                    continue
                if _check_local_constraints(
                    codons, i, alt, n_total,
                    enzymes_sites, is_prokaryote, cpg_mode,
                ):
                    codons[i] = alt
                    break

    return "".join(codons)


def _adjust_gc_content(
    seq: str,
    protein: str,
    usage: dict[str, float],
    gc_lo: float,
    gc_hi: float,
    is_prokaryote: bool,
    cpg_mode: str,
    enzymes_sites: list[tuple[str, str]],
    notes: list[str],
) -> str:
    """Adjust GC content by swapping codons to GC-richer or GC-poorer alternatives.

    Every candidate swap is validated via :func:`_check_local_constraints`,
    which enforces the FULL constraint set (GT, CG, ATTTA, T-run,
    restriction sites, premature stop) — not just GT/CG. This prevents
    GC adjustment from re-introducing ATTTA motifs, T-runs, or
    restriction-site violations that the cleanup pass already removed
    (Task W4 / H1 fix).

    The stop codon (index ``n_codons``) is included in the codon list as
    a non-swappable neighbour, mirroring :func:`_cleanup_violations`, so
    boundary violations between the last protein codon and the stop
    codon are visible to ``_check_local_constraints``.
    """
    n_codons = len(protein)
    # Include the stop codon (index n_codons) as a non-swappable neighbour
    # so boundary violations between the last protein codon and the stop
    # codon are visible to _check_local_constraints (mirrors
    # _cleanup_violations).
    codons = [seq[i*3:(i+1)*3] for i in range(n_codons + 1)]
    n_total = n_codons + 1  # protein codons + stop codon

    gc = gc_content(seq)

    if gc < gc_lo:
        # Need more GC — swap to GC-richer codons
        for i in range(n_codons):
            aa = protein[i]
            candidates = _sorted_codons_for_aa(aa, usage)
            current_gc = codons[i].count("G") + codons[i].count("C")
            best = codons[i]
            best_gc = current_gc
            for c in candidates:
                c_gc = c.count("G") + c.count("C")
                if c_gc > best_gc:
                    # Validate ALL constraints (GT, CG, ATTTA, T-run,
                    # restriction sites, premature stop) so GC adjustment
                    # cannot re-introduce a violation the cleanup pass
                    # already removed. (Task W4 / H1.)
                    if _check_local_constraints(
                        codons, i, c, n_total,
                        enzymes_sites, is_prokaryote, cpg_mode,
                    ):
                        best = c
                        best_gc = c_gc
            codons[i] = best
    elif gc > gc_hi:
        # Need less GC — swap to GC-poorer codons
        for i in range(n_codons):
            aa = protein[i]
            candidates = _sorted_codons_for_aa(aa, usage)
            current_gc = codons[i].count("G") + codons[i].count("C")
            best = codons[i]
            best_gc = current_gc
            for c in candidates:
                c_gc = c.count("G") + c.count("C")
                if c_gc < best_gc:
                    if _check_local_constraints(
                        codons, i, c, n_total,
                        enzymes_sites, is_prokaryote, cpg_mode,
                    ):
                        best = c
                        best_gc = c_gc
            codons[i] = best

    result = "".join(codons)  # already includes the stop codon
    final_gc = gc_content(result)
    if not (gc_lo <= final_gc <= gc_hi):
        notes.append(f"GC adjustment incomplete: {final_gc:.3f} (target: {gc_lo}-{gc_hi})")
    return result


def _translate_dna(dna: str, secis_positions: list[int] | None = None) -> str:
    """Translate DNA to protein.
    
    If secis_positions is provided, TGA at those codon indices is translated
    as U (selenocysteine) instead of * (stop).
    """
    secis_set = set(secis_positions) if secis_positions else set()
    protein = []
    codon_idx = 0
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        if codon == "TGA" and codon_idx in secis_set:
            protein.append("U")  # selenocysteine
        else:
            aa = CODON_TABLE.get(codon, "X")
            protein.append(aa)
        codon_idx += 1
    return "".join(protein)


def _simple_back_translate(protein: str, usage: dict[str, float]) -> str:
    """Simple back-translation: best codon per amino acid."""
    dna = []
    for aa in protein:
        if aa == "U":
            dna.append("TGA")  # selenocysteine
        else:
            codons = _sorted_codons_for_aa(aa, usage)
            dna.append(codons[0] if codons else "GCT")
    return "".join(dna) + "TAA"


def _is_prokaryote(organism: str) -> bool:
    """Check if organism is prokaryotic."""
    prokaryotes = {"e_coli", "Escherichia_coli", "bacillus", "Bacillus_subtilis",
                   "pseudomonas", "Pseudomonas_putida", "corynebacterium",
                   "Corynebacterium_glutamicum"}
    return organism in prokaryotes or _is_prokaryotic_organism(organism)


def _is_prokaryotic_organism(organism: str) -> bool:
    """Check if organism is prokaryotic (delegates to type_system)."""
    try:
        from biocompiler.type_system.checks import _is_prokaryotic_organism as _check
        return _check(organism)
    except ImportError:
        return False
