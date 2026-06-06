"""
BioCompiler GC Adjustment Logic
=================================
GC content adjustment and reconciliation for the optimization pipeline.

Extracted from optimization.py for maintainability.
"""

import logging
from typing import Any

from .constants import reverse_complement
from .organisms import ORGANISM_GC_TARGETS
from .aho_corasick import AhoCorasickScanner  # type: ignore[import-untyped]

# ── NUMBA integration ──────────────────────────────────────────────
try:
    from .numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        count_gc as _numba_count_gc,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False

logger = logging.getLogger(__name__)

MAX_GC_ADJUSTMENT_ITERATIONS: int = 200


def _compute_gc_count(sequence: str) -> int:
    """Compute GC base count, using NUMBA if available."""
    if _HAS_NUMBA:
        try:
            return _numba_count_gc(_seq_to_bytes(sequence))
        except Exception:
            return sum(1 for b in sequence if b in "GC")
    return sum(1 for b in sequence if b in "GC")


def adjust_gc_content(
    sequence: str,
    aas: list[str],
    sorted_codons: dict[str, list[str]],
    usage: dict[str, float],
    gc_lo: float,
    gc_hi: float,
    organism: str,
    concrete_sites: list[str],
    concrete_scanner: Any | None,
    remove_site_multicodon_fn: Any,
) -> tuple[str, list[str]]:
    """Adjust GC content to fall within [gc_lo, gc_hi].

    This implements the GC adjustment step and the subsequent reconciliation
    step that checks if GC adjustment reintroduced restriction sites.

    Strategy: GC must be in [gc_lo, gc_hi] (hard constraint).
    If in range, we gently nudge toward organism target but NEVER at the
    cost of significant CAI reduction. The organism GC target is aspirational,
    not mandatory — a sequence with CAI=0.99 and GC=0.61 (slightly above
    human's 0.41 target) is better than CAI=0.82 and GC=0.46.

    Args:
        sequence: Current DNA sequence (uppercase).
        aas: List of amino acid codes (one per codon).
        sorted_codons: AA → codons sorted by CAI (descending).
        usage: Codon → CAI adaptiveness weight dict.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        organism: Canonical organism name (for GC target lookup).
        concrete_sites: List of concrete restriction site sequences.
        concrete_scanner: AhoCorasickScanner for concrete sites (or None).
        remove_site_multicodon_fn: Reference to _remove_site_multicodon function.

    Returns:
        Tuple of (adjusted_sequence, warnings_list).
    """
    warnings: list[str] = []
    n_bases = len(sequence)
    gc_count = _compute_gc_count(sequence)
    gc_val = gc_count / n_bases
    organism_gc_range = ORGANISM_GC_TARGETS.get(organism, (gc_lo, gc_hi))
    organism_gc = (organism_gc_range[0] + organism_gc_range[1]) / 2.0
    target_gc = max(gc_lo, min(gc_hi, organism_gc))

    gc_out_of_range = not (gc_lo <= gc_val <= gc_hi)

    if gc_out_of_range:
        # Hard constraint: MUST get GC into range
        # Target the nearest bound
        if gc_val < gc_lo:
            phase_target = gc_lo
        else:
            phase_target = gc_hi

        for iteration in range(MAX_GC_ADJUSTMENT_ITERATIONS):
            if gc_lo <= gc_val <= gc_hi:
                break
            best_alt = None
            best_ci = -1
            best_diff = abs(gc_val - phase_target)
            best_gc_delta = 0
            best_cai = -1.0  # CAI tiebreaker: prefer higher CAI among equal GC improvement
            for ci in range(len(aas)):
                aa = aas[ci]
                current = sequence[ci * 3:ci * 3 + 3]
                current_gc = sum(1 for b in current if b in "GC")
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    alt_gc = sum(1 for b in alt if b in "GC")
                    new_gc_count = gc_count - current_gc + alt_gc
                    new_frac = new_gc_count / n_bases
                    diff = abs(new_frac - phase_target)
                    alt_cai = usage.get(alt, 0.0)
                    # Prefer better GC improvement; among equal GC improvement, prefer higher CAI
                    if diff < best_diff or (diff == best_diff and alt_cai > best_cai):
                        best_diff = diff
                        best_alt = alt
                        best_ci = ci
                        best_gc_delta = alt_gc - current_gc
                        best_cai = alt_cai
            if best_alt is None:
                break
            # PERF (Optimization D): Use list mutation for codon swap
            seq_list = list(sequence)
            seq_list[best_ci * 3: best_ci * 3 + 3] = list(best_alt)
            sequence = "".join(seq_list)
            gc_count += best_gc_delta
            gc_val = gc_count / n_bases
        else:
            warnings.append(f"GC adjustment: max iterations reached, current GC={gc_val:.3f}")

    # Step: Reconciliation — check if GC adjustment reintroduced restriction sites
    # Use Aho-Corasick scanner for fast multi-site detection if available
    if concrete_scanner is not None:
        remaining_matches = concrete_scanner.scan(sequence)
        for pos, site_match, _label in remaining_matches:
            site_rc = reverse_complement(site_match)
            # Try one more round of multi-codon removal
            new_seq, fixed = remove_site_multicodon_fn(
                sequence, aas, sorted_codons, site_match, site_rc, usage=usage
            )
            if fixed:
                sequence = new_seq
                # Re-check GC
                gc_count = _compute_gc_count(sequence)
                gc_val = gc_count / n_bases
                if not (gc_lo <= gc_val <= gc_hi):
                    # GC drifted — try to fix with single-codon swaps that don't reintroduce sites
                    for ci in range(len(aas)):
                        aa = aas[ci]
                        current = sequence[ci * 3:ci * 3 + 3]
                        current_gc = sum(1 for b in current if b in "GC")
                        for alt in sorted_codons[aa]:
                            if alt == current:
                                continue
                            alt_gc = sum(1 for b in alt if b in "GC")
                            new_gc_count = gc_count - current_gc + alt_gc
                            new_frac = new_gc_count / n_bases
                            # Check this swap doesn't reintroduce any site
                            seq_list = list(sequence)
                            seq_list[ci * 3: ci * 3 + 3] = list(alt)
                            test = "".join(seq_list)
                            site_ok = not concrete_scanner.has_any_match(test)
                            if site_ok and abs(new_frac - target_gc) < abs(gc_val - target_gc):
                                sequence = test
                                gc_count = new_gc_count
                                gc_val = gc_count / n_bases
                                break
    else:
        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            if site_upper in sequence or site_rc in sequence:
                # Try one more round of multi-codon removal
                new_seq, fixed = remove_site_multicodon_fn(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    # Re-check GC
                    gc_count = _compute_gc_count(sequence)
                    gc_val = gc_count / n_bases
                    if not (gc_lo <= gc_val <= gc_hi):
                        # GC drifted — try to fix with single-codon swaps that don't reintroduce sites
                        for ci in range(len(aas)):
                            aa = aas[ci]
                            current = sequence[ci * 3:ci * 3 + 3]
                            current_gc = sum(1 for b in current if b in "GC")
                            for alt in sorted_codons[aa]:
                                if alt == current:
                                    continue
                                alt_gc = sum(1 for b in alt if b in "GC")
                                new_gc_count = gc_count - current_gc + alt_gc
                                new_frac = new_gc_count / n_bases
                                # Check this swap doesn't reintroduce any site
                                seq_list = list(sequence)
                                seq_list[ci * 3: ci * 3 + 3] = list(alt)
                                test = "".join(seq_list)
                                site_ok = all(
                                    s not in test and reverse_complement(s) not in test
                                    for s in concrete_sites
                                )
                                if site_ok and abs(new_frac - target_gc) < abs(gc_val - target_gc):
                                    sequence = test
                                    gc_count = new_gc_count
                                    gc_val = gc_count / n_bases
                                    break
            else:
                # Could not remove — already warned in Remove Restriction Sites step
                pass

    return sequence, warnings
