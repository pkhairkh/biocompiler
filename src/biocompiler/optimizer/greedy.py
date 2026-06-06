"""
Greedy optimizer: multi-step greedy codon optimization.

Contains the main greedy optimization loop, splice donor potential scoring,
GT-aware codon selection, eukaryote-specific CAI recovery, CpG elimination,
IUPAC site expansion, and predicate checking.
"""

from typing import Callable, List, Dict, Optional, Tuple, Set, Any

import logging
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import product as itertools_product

from ..type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, PredicateResult,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_avoidable_gt,
    check_no_gt_dinucleotide_soft, _compute_max_gt_count,
    check_valid_coding_seq,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
)
from ..organisms import CODON_ADAPTIVENESS_TABLES, ORGANISM_GC_TARGETS, resolve_organism, ORGANISM_ALIASES, SPECIES_SHORT_NAMES, SUPPORTED_ORGANISMS
from ..constants import reverse_complement, RESTRICTION_ENZYMES, IUPAC_EXPAND, VALID_IUPAC_BASES

# Optimization thresholds and sentinel values
# (Originally defined inline in optimization.py; moved here during decomposition.)
IUPAC_EXPANSION_CAP: int = 4096
ELIMINATED_SITE_SCORE: float = -999.0
TOP_CAI_ALTERNATIVES: int = 3
T_RUN_LENGTH_THRESHOLD: int = 6
SPLICE_DONOR_POTENTIAL_THRESHOLD: float = 0.5
_MAX_ACCEPTOR_SEARCH_DIST: int = 200
EUKARYOTE_CAI_GT_COST_THRESHOLD: float = 0.02
GT_BOUNDARY_CAI_TOLERANCE: float = 0.03
GT_CAI_LOG_ADAPTIVENESS_COST: float = 0.03
from ..scanner import gc_content
from ..sliding_gc import check_sliding_gc, fix_sliding_gc_violations
from ..maxentscan import score_donor, score_acceptor, max_donor_score, max_acceptor_score
from ..mutagenesis import propose_mutagenesis, MutagenesisReport, MutagenesisProposal
from ..incremental import IncrementalSequenceState, CodonCache, EnzymeSiteCache
from ..certificate import format_certificate
from ..exceptions import InvalidProteinError, UnsupportedOrganismError, OptimizationConstraintError, BiosecurityError
from ..decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
)
from ..aho_corasick import AhoCorasickScanner, build_scanner_from_enzymes, build_scanner_from_sites  # type: ignore[import-untyped]
from ..objectives import resolve_objective as _resolve_objective, cai_objective as _cai_objective

# Internal optimizer sub-module imports
from .cai import HAS_NUMBA, _HAS_NUMBA, _compute_cai_fast, _count_dinucs_fast, _BatchSwapScorer
from .constraints import (
    _find_site_in_sequence, _get_overlapping_codons, _remove_site_multicodon,
    _find_gt_free_codons, _find_ag_free_codons,
    _organism_to_species_key, _species_key_to_organism,
    _back_translate_protein,
    _count_gts, _is_unavoidable_gt, _has_gt, _codon_creates_boundary_gt,
)
from .utils import (
    protein_to_aa_list, ConvergenceTracker, OptimizationResult, FullConstructResult,
    MAX_RESTRICTION_SITE_ITERATIONS, MAX_IUPAC_SITE_ITERATIONS,
    MAX_ATTTA_MOTIF_ITERATIONS, MAX_T_RUN_ITERATIONS,
    MAX_GC_ADJUSTMENT_ITERATIONS, MAX_SPLICE_ELIMINATION_ITERATIONS,
    MAX_CPG_DISRUPTION_ITERATIONS,
    DEFAULT_MAX_ITERATIONS,
    CONVERGENCE_IMPROVEMENT_THRESHOLD,
    CONVERGENCE_PATIENCE,
    OSCILLATION_WINDOW,
)

logger = logging.getLogger(__name__)

def score_splice_donor_potential(dna: str, position: int) -> float:
    """Score how likely a GT dinucleotide is to function as a splice donor.

    Not all GT dinucleotides are splice donors.  A true splice donor needs:
    - GT at the 5' end of an intron
    - AG at the 3' end (acceptor) typically 30-200 nt downstream
    - A polypyrimidine tract upstream of the AG
    - A branch point (YNYRAY) ~18-40 nt upstream of the AG

    This function returns a score from 0.0 (unlikely splice donor) to 1.0
    (strong, functional splice donor) based on:

    1. **MaxEntScan donor score** (primary, ~70% weight): How well the
       9-mer context matches the splice donor PWM.
    2. **Downstream acceptor context** (~20% weight): Whether there's a
       downstream AG with a reasonable acceptor score within intron-length
       distance.
    3. **Polypyrimidine tract quality** (~10% weight): Whether there's a
       C/T-rich region between the GT and the best downstream AG.

    Score interpretation:
    - < 0.3: Low risk — GT very unlikely to function as splice donor.
      In-codon GTs from optimal codons almost always score here.
    - 0.3–0.5: Moderate risk — possible cryptic splice site.
    - > 0.5: High risk — likely functional splice donor; should be avoided.

    Pre-conditions:
    - dna is a valid uppercase DNA string
    - 0 <= position < len(dna) - 1
    - dna[position:position+2] == "GT"

    Post-conditions:
    - returns a float in [0.0, 1.0]
    - deterministic: same input always produces same output
    """
    # ── Component 1: MaxEntScan donor score ──
    donor_score = score_donor(dna, position)

    # Normalize donor score to [0, 1].
    # Typical ranges:
    #   -5 to 0:  non-donor
    #   0 to 3:   very weak
    #   3 to 8:   moderate cryptic
    #   8+:       strong canonical
    if donor_score < 0:
        donor_component = 0.0
    elif donor_score < 3.0:
        donor_component = 0.05 + 0.10 * (donor_score / 3.0)   # 0.05–0.15
    elif donor_score < 8.0:
        donor_component = 0.15 + 0.35 * ((donor_score - 3.0) / 5.0)  # 0.15–0.50
    else:
        donor_component = 0.50 + min(0.20, 0.20 * ((donor_score - 8.0) / 5.0))  # 0.50–0.70

    # ── Component 2: Downstream acceptor context ──
    # Search for the strongest AG acceptor downstream within intron distance.
    acceptor_component = 0.0
    best_ag_score = -999.0
    best_ag_pos = -1
    search_end = min(len(dna) - 1, position + _MAX_ACCEPTOR_SEARCH_DIST)

    for i in range(position + 4, search_end):  # +4: skip the GT itself
        if dna[i] == "A" and dna[i + 1] == "G":
            ag_score = score_acceptor(dna, i)
            if ag_score > best_ag_score:
                best_ag_score = ag_score
                best_ag_pos = i

    if best_ag_pos >= 0:
        dist = best_ag_pos - position
        # Typical intron: 30-5000 nt; short introns in yeast ~80-500 nt
        # Only count acceptors at reasonable intron-like distances
        if 20 <= dist <= _MAX_ACCEPTOR_SEARCH_DIST:
            if best_ag_score >= 8.0:
                acceptor_component = 0.20  # Strong acceptor
            elif best_ag_score >= 3.0:
                acceptor_component = 0.12  # Moderate acceptor
            else:
                acceptor_component = 0.05  # Weak acceptor

    # ── Component 3: Polypyrimidine tract quality ──
    # Check for C/T-rich region upstream of the best AG (between GT and AG).
    py_component = 0.0
    if best_ag_pos >= 0:
        tract_start = max(position + 4, best_ag_pos - 30)
        tract_end = max(position + 4, best_ag_pos - 5)
        if tract_end > tract_start:
            tract = dna[tract_start:tract_end]
            py_fraction = (tract.count("C") + tract.count("T")) / len(tract)
            if py_fraction >= 0.7:
                py_component = 0.10  # Good polypyrimidine tract
            elif py_fraction >= 0.5:
                py_component = 0.05  # Moderate

    # Combined score
    return min(1.0, max(0.0, donor_component + acceptor_component + py_component))


# ==============================================================================
# Eukaryotic CAI Recovery
# ==============================================================================
# NOTE (Task 1.8): EUKARYOTE_CAI_GT_COST_THRESHOLD, GT_BOUNDARY_CAI_TOLERANCE,
# and GT_CAI_LOG_ADAPTIVENESS_COST are now imported from .constants via aliased
# imports at the top of this file.


def _gt_aware_select_codon(
    aa: str,
    next_aa: str | None,
    sorted_codons: dict[str, list[str]],
    usage: dict[str, float],
    cai_tolerance: float = GT_BOUNDARY_CAI_TOLERANCE,
) -> str:
    """Select the best codon for an amino acid, preferring GT-free boundary if
    the CAI cost is within tolerance.

    This is a "best effort" approach that maximizes CAI while minimizing GT
    dinucleotides at codon boundaries. For eukaryotes, if the optimal codon
    would create a GT at the boundary with the next codon (i.e., the optimal
    codon ends with 'G' and the next codon's optimal starts with 'T'), we
    check if there's an alternative codon with similar CAI (within
    ``cai_tolerance`` relative) that avoids the boundary GT.

    Decision logic:
    1. Select the optimal (highest-CAI) codon for this amino acid.
    2. If no next amino acid, or the optimal codon doesn't create a boundary
       GT, use it.
    3. If the optimal codon WOULD create a boundary GT with the next codon's
       optimal start:
       a. Check each alternative codon (in CAI order) to see if it avoids
          the boundary GT and has CAI within tolerance.
       b. If such an alternative exists, use it (best CAI among valid alts).
       c. If no alternative within CAI tolerance, use the optimal codon and
          accept the GT (CAI > GT avoidance).

    For prokaryotes, this function is never called (GT avoidance is skipped).

    Args:
        aa: Amino acid at the current position.
        next_aa: Amino acid at the next position (None if last position or
            next is a stop codon).
        sorted_codons: Codons per AA sorted by CAI descending.
        usage: Codon adaptiveness table (codon → w value).
        cai_tolerance: Maximum relative CAI loss acceptable for avoiding a
            boundary GT. Default: 0.10 (10%).

    Returns:
        The selected codon string.
    """
    candidates = sorted_codons.get(aa, [])
    if not candidates:
        return AA_TO_CODONS.get(aa, [""])[0]

    optimal = candidates[0]
    optimal_cai = usage.get(optimal, 0.0)

    # If no next AA, just use the optimal codon
    if next_aa is None:
        return optimal

    # Determine what the next position's first base would be if we use
    # the optimal codon for next_aa.  We need to look ahead to see if
    # our codon's last base + next codon's first base = "GT".
    next_candidates = sorted_codons.get(next_aa, [])
    if not next_candidates:
        return optimal

    next_optimal = next_candidates[0]
    # Check if optimal + next_optimal creates a boundary GT
    if optimal[-1] + next_optimal[0] != "GT":
        return optimal  # No boundary GT — use optimal

    # Boundary GT detected.  Look for alternatives within CAI tolerance.
    min_acceptable_cai = optimal_cai * (1.0 - cai_tolerance)
    best_alt = None
    best_alt_cai = -1.0

    for alt in candidates:
        if alt == optimal:
            continue
        alt_cai = usage.get(alt, 0.0)
        # Check if this alternative avoids the boundary GT
        if alt[-1] + next_optimal[0] == "GT":
            continue  # Still creates boundary GT
        # Check CAI within tolerance
        if alt_cai >= min_acceptable_cai and alt_cai > best_alt_cai:
            best_alt = alt
            best_alt_cai = alt_cai

    if best_alt is not None:
        return best_alt  # Alternative found within CAI tolerance

    # No alternative within tolerance — accept the boundary GT
    return optimal


def _is_in_codon_gt(seq: str, pos: int) -> bool:
    """Check whether a GT dinucleotide at position *pos* is entirely within
    a single codon (in-codon) rather than spanning a codon boundary (cross-codon).

    Pre-conditions:
    - seq is a valid uppercase DNA string
    - 0 <= pos < len(seq) - 1
    - seq[pos:pos+2] == "GT"

    Post-conditions:
    - Returns True if both the G and T are in the same codon
    - Returns False if the GT spans a codon boundary (G at end of one codon,
      T at start of the next)
    """
    codon_of_g = pos // 3
    codon_of_t = (pos + 1) // 3
    return codon_of_g == codon_of_t


def _eukaryote_cai_recovery(
    sequence: str,
    protein: str,
    usage: dict[str, float],
    enzymes: list[str] | None = None,
    cai_cost_threshold: float = EUKARYOTE_CAI_GT_COST_THRESHOLD,
) -> tuple[str, int]:
    """Recover CAI by swapping suboptimal codons to optimal ones for eukaryotes.

    For eukaryotes, GT avoidance should be a SOFT preference, not a hard
    constraint.  When the optimizer has replaced an optimal codon with a
    suboptimal alternative (typically to avoid GT dinucleotides), this
    function swaps back using **cost-aware GT resolution** with splice
    donor scoring:

    - If the optimal codon creates a GT with **low** splice donor potential
      (score < SPLICE_DONOR_POTENTIAL_THRESHOLD), always swap to optimal —
      the GT is unlikely to be a real splice donor, so CAI wins.
    - If the optimal codon creates a GT with **high** splice donor potential
      (score >= SPLICE_DONOR_POTENTIAL_THRESHOLD), only swap if the CAI cost
      exceeds ``cai_cost_threshold`` — the GT might be dangerous, so keep
      the suboptimal codon unless the CAI loss is too large.
    - If the optimal codon does NOT create any new GT dinucleotides, always
      swap (no tradeoff needed).

    Pre-conditions:
    - sequence translates to protein
    - len(sequence) == len(protein) * 3
    - usage maps codon strings to adaptiveness values (0.0–1.0)

    Post-conditions:
    - returned sequence translates to the same protein
    - CAI of returned sequence >= CAI of input sequence
    - no new restriction sites are introduced

    Args:
        sequence: Current optimized DNA sequence.
        protein: Amino acid sequence (no stop).
        usage: Codon adaptiveness table (codon → w value).
        enzymes: List of restriction enzyme names to avoid creating sites for.
        cai_cost_threshold: Maximum acceptable CAI cost for GT avoidance.
            If using a GT-free codon would drop CAI by more than this,
            keep the GT-containing optimal codon.

    Returns:
        Tuple of (recovered sequence, number of codons upgraded).
    """
    # Precompute optimal codons per amino acid
    optimal_codons: dict[str, str] = {}
    for aa in set(protein):
        if aa == "*":
            continue
        codons = AA_TO_CODONS.get(aa, [])
        if not codons:
            continue
        best = max(codons, key=lambda c: usage.get(c, 0.0))
        optimal_codons[aa] = best

    # Precompute restriction site sequences to check
    rs_sites: list[tuple[str, str]] = []
    max_site_len = 0
    if enzymes:
        from ..restriction_sites import get_recognition_site as _get_site
        for enz in enzymes:
            site = _get_site(enz)
            if site is not None:
                site_rc = reverse_complement(site)
                rs_sites.append((site, site_rc))
                max_site_len = max(max_site_len, len(site))

    seq_list = list(sequence)
    n_codons = len(protein)
    upgrades = 0

    # Local-region buffer sizes for splice donor context checking.
    # score_splice_donor_potential() needs up to 3 bases upstream (donor
    # 9-mer context) and _MAX_ACCEPTOR_SEARCH_DIST (80) bases downstream
    # (acceptor search).  Using a local region avoids O(n) full-sequence
    # joins inside the per-codon loop (Fix A: O(n²) → O(n)).
    _LOCAL_GT_UPSTREAM = 4    # bases upstream of codon for donor context
    _LOCAL_GT_DOWNSTREAM = _MAX_ACCEPTOR_SEARCH_DIST + 3  # downstream for acceptor

    def _should_swap_to_optimal(ci: int, current: str, optimal: str,
                                 cai_cost: float) -> bool:
        """Decide whether to swap current → optimal at codon position ci.

        Uses cost-aware GT resolution with splice donor scoring:
        - If optimal doesn't introduce new GTs: always swap (pure CAI gain).
        - If optimal introduces GT(s) with low splice donor potential: swap
          (GT is unlikely to be functional).
        - If optimal introduces GT(s) with high splice donor potential: only
          swap if CAI cost exceeds the threshold (tradeoff).

        PERF (Fix A): Uses local-region string construction instead of
        rebuilding the entire sequence.  Only the ~90-char window around
        the swapped codon is joined, avoiding O(n) per-codon overhead.
        """
        # Build a local region around the swapped codon instead of the
        # full sequence.  This is sufficient for GT detection (which
        # only spans the codon ±1 base) and splice donor scoring (which
        # needs ~3 bases upstream and ~80 bases downstream).
        region_start = max(0, ci * 3 - _LOCAL_GT_UPSTREAM)
        region_end = min(len(seq_list), ci * 3 + _LOCAL_GT_DOWNSTREAM)

        # Current local region (before swap)
        current_region = "".join(seq_list[region_start:region_end])

        # Test local region (after swap) — copy the slice and modify
        test_region_chars = list(seq_list[region_start:region_end])
        local_codon_start = ci * 3 - region_start
        test_region_chars[local_codon_start] = optimal[0]
        test_region_chars[local_codon_start + 1] = optimal[1]
        test_region_chars[local_codon_start + 2] = optimal[2]
        test_region = "".join(test_region_chars)

        # Find all GT positions in the swapped region that would be new
        # (present in test_region but not in current_region)
        check_start_local = max(0, local_codon_start - 1)
        check_end_local = min(len(test_region), local_codon_start + 4)

        new_gt_positions = []
        for p in range(check_start_local, check_end_local - 1):
            if test_region[p:p+2] == "GT" and current_region[p:p+2] != "GT":
                new_gt_positions.append(p)

        if not new_gt_positions:
            # No new GTs introduced — pure CAI gain, always swap
            return True

        # Check splice donor potential of each new GT using the local region.
        # The region is sized to contain the full donor+acceptor context.
        max_splice_score = 0.0
        for gt_pos in new_gt_positions:
            sdp = score_splice_donor_potential(test_region, gt_pos)
            if sdp > max_splice_score:
                max_splice_score = sdp

        if max_splice_score < SPLICE_DONOR_POTENTIAL_THRESHOLD:
            # New GTs have low splice donor potential — accept them for CAI
            return True

        # New GTs have high splice donor potential — only swap if CAI
        # cost is significant enough to justify the risk
        return cai_cost > cai_cost_threshold

    # ── Pass 1: Single-codon swaps ──
    # For each position where the current codon is suboptimal, decide
    # whether to swap using cost-aware GT resolution.
    blocked_positions: list[int] = []  # Positions blocked by restriction sites

    for ci in range(n_codons):
        aa = protein[ci]
        if aa == "*":
            continue
        current = "".join(seq_list[ci * 3:ci * 3 + 3])
        optimal = optimal_codons.get(aa)
        if optimal is None or current == optimal:
            continue

        # Check CAI cost of using the suboptimal codon
        current_w = usage.get(current, 0.0)
        optimal_w = usage.get(optimal, 0.0)
        cai_cost = optimal_w - current_w

        # Cost-aware decision: should we swap?
        if not _should_swap_to_optimal(ci, current, optimal, cai_cost):
            continue

        # Swap to optimal codon
        old_codon = current
        seq_list[ci * 3] = optimal[0]
        seq_list[ci * 3 + 1] = optimal[1]
        seq_list[ci * 3 + 2] = optimal[2]

        # Check for new restriction sites in local region
        # PERF (Fix A): Build only the local region instead of the full sequence
        if rs_sites and max_site_len > 0:
            check_start = max(0, ci * 3 - max_site_len + 1)
            check_end = min(len(seq_list), ci * 3 + 3 + max_site_len - 1)
            local_region = "".join(seq_list[check_start:check_end])
            site_found = False
            for site, site_rc in rs_sites:
                if site in local_region or (site_rc and site_rc in local_region):
                    site_found = True
                    break
            if site_found:
                # Undo swap — restriction site would be created
                seq_list[ci * 3] = old_codon[0]
                seq_list[ci * 3 + 1] = old_codon[1]
                seq_list[ci * 3 + 2] = old_codon[2]
                blocked_positions.append(ci)
                continue

        upgrades += 1

    # ── Pass 2: Coordinated swaps for positions blocked by restriction sites ──
    # When a single-codon swap is blocked because it would create a restriction
    # site that spans a codon boundary, we can sometimes resolve it by also
    # changing an adjacent codon to break the restriction site.  This is only
    # done when the combined CAI improvement outweighs the CAI cost of the
    # adjacent codon change.
    if blocked_positions and rs_sites:
        # Precompute sorted codons per amino acid (by CAI descending)
        sorted_codons_map: dict[str, list[str]] = {}
        for aa in set(protein):
            if aa == "*":
                continue
            codons = AA_TO_CODONS.get(aa, [])
            sorted_codons_map[aa] = sorted(
                codons, key=lambda c: usage.get(c, 0.0), reverse=True
            )

        for ci in blocked_positions:
            aa = protein[ci]
            current = "".join(seq_list[ci * 3:ci * 3 + 3])
            optimal = optimal_codons.get(aa)
            if optimal is None or current == optimal:
                continue
            current_w = usage.get(current, 0.0)
            optimal_w = usage.get(optimal, 0.0)
            benefit = optimal_w - current_w
            if benefit <= cai_cost_threshold:
                continue

            # Check if the swap is justified by splice donor scoring
            if not _should_swap_to_optimal(ci, current, optimal, benefit):
                continue

            # Try adjacent codons to break the restriction site
            best_combo: tuple[int, str, float] | None = None  # (adj_ci, adj_codon, net_cai)
            for adj_offset in [-2, -1, 1, 2]:
                adj_ci = ci + adj_offset
                if adj_ci < 0 or adj_ci >= n_codons:
                    continue
                adj_aa = protein[adj_ci]
                if adj_aa == "*":
                    continue
                adj_current = "".join(seq_list[adj_ci * 3:adj_ci * 3 + 3])
                adj_current_w = usage.get(adj_current, 0.0)
                for adj_alt in sorted_codons_map.get(adj_aa, []):
                    if adj_alt == adj_current:
                        continue
                    # Apply both swaps and check for restriction sites
                    adj_old = adj_current
                    seq_list[ci * 3] = optimal[0]
                    seq_list[ci * 3 + 1] = optimal[1]
                    seq_list[ci * 3 + 2] = optimal[2]
                    seq_list[adj_ci * 3] = adj_alt[0]
                    seq_list[adj_ci * 3 + 1] = adj_alt[1]
                    seq_list[adj_ci * 3 + 2] = adj_alt[2]

                    # PERF (Fix A): Build only the local region instead of full sequence
                    check_start = max(0, min(ci, adj_ci) * 3 - max_site_len + 1)
                    check_end = min(len(seq_list), max(ci, adj_ci) * 3 + 3 + max_site_len - 1)
                    local_region = "".join(seq_list[check_start:check_end])
                    site_found = False
                    for site, site_rc in rs_sites:
                        if site in local_region or (site_rc and site_rc in local_region):
                            site_found = True
                            break

                    # Undo both swaps
                    seq_list[ci * 3] = current[0]
                    seq_list[ci * 3 + 1] = current[1]
                    seq_list[ci * 3 + 2] = current[2]
                    seq_list[adj_ci * 3] = adj_old[0]
                    seq_list[adj_ci * 3 + 1] = adj_old[1]
                    seq_list[adj_ci * 3 + 2] = adj_old[2]

                    if not site_found:
                        # Check net CAI: benefit of optimal at ci minus cost of adj_alt
                        adj_alt_w = usage.get(adj_alt, 0.0)
                        adj_cost = adj_current_w - adj_alt_w  # positive = CAI loss
                        net_cai = benefit - adj_cost
                        # Only accept if net CAI is positive AND the coordinated
                        # swap is better than keeping the status quo
                        if net_cai > 0 and (best_combo is None or net_cai > best_combo[2]):
                            best_combo = (adj_ci, adj_alt, net_cai)

            if best_combo is not None:
                adj_ci, adj_codon, _net = best_combo
                # Apply both swaps
                seq_list[ci * 3] = optimal[0]
                seq_list[ci * 3 + 1] = optimal[1]
                seq_list[ci * 3 + 2] = optimal[2]
                seq_list[adj_ci * 3] = adj_codon[0]
                seq_list[adj_ci * 3 + 1] = adj_codon[1]
                seq_list[adj_ci * 3 + 2] = adj_codon[2]
                upgrades += 2  # Count both the optimal swap and the adjacent change

    return "".join(seq_list), upgrades


# ==============================================================================
# Systematic CpG Dinucleotide Elimination
# ==============================================================================

def _eliminate_cpg_dinucleotides(
    sequence: str,
    protein: str,
    usage: dict[str, float],
    enzymes: list[str] | None = None,
    max_iterations: int = 50,
    cpg_window: int = 200,
    cpg_threshold: float = 0.6,
    organism: str = "",
    gc_lo: float = 0.0,
    gc_hi: float = 1.0,
    max_cai_cost: float = 1.0,
) -> tuple[str, list[str]]:
    """Systematically eliminate CpG dinucleotides to avoid CpG islands.

    This is a post-CAI-maximization pass that scans for ALL CG dinucleotides
    (both within codons and at codon boundaries) and attempts to replace them
    with synonymous codons that:
    1. Eliminate the specific CG dinucleotide
    2. Minimize CAI loss (prefer highest-CAI alternative)
    3. Do not create restriction sites
    4. Do not create new CG dinucleotides at the same or adjacent positions
    5. Keep GC content within the target range [gc_lo, gc_hi]

    Unlike the previous "best-effort" approach which broke on the first
    unfixable CpG position or returned early when the CpG island check
    passed, this function:
    - Always attempts to eliminate ALL CG dinucleotides, even if the
      sequence already passes the CpG island Obs/Exp ratio check
    - Continues trying ALL CpG positions even if some fail
    - Makes multiple passes to handle cascading fixes
    - For within-codon CpG: replaces with the highest-CAI synonymous codon
      that does not contain "CG"
    - For boundary CpG: tries changing the downstream codon first (avoid
      starting with G), then the upstream codon (avoid ending with C),
      then coordinated 2-codon swap — always preferring minimal CAI loss
    - Tracks and reports which CpGs were eliminated and which remain
    - Only stops when no more CG dinucleotides can be eliminated or
      max_iterations is reached

    Pre-conditions:
    - sequence translates to protein
    - len(sequence) == len(protein) * 3
    - usage maps codon strings to adaptiveness values (0.0–1.0)

    Post-conditions:
    - returned sequence translates to the same protein
    - CG dinucleotide count in returned sequence <= CG count in input sequence
    - warnings list describes any remaining CpG positions

    Args:
        sequence: Current optimized DNA sequence.
        protein: Amino acid sequence (no stop).
        usage: Codon adaptiveness table (codon → w value).
        enzymes: List of restriction enzyme names to avoid creating sites for.
        max_iterations: Maximum number of elimination passes.
        cpg_window: Window size for CpG island detection.
        cpg_threshold: Obs/Exp ratio threshold for CpG islands.
        organism: Target organism (for prokaryote skip).
        gc_lo: Minimum GC content fraction.
        gc_hi: Maximum GC content fraction.

    Returns:
        Tuple of (optimized sequence, list of warning strings about remaining CpGs).
    """
    # Skip for prokaryotic organisms — CpG islands are a eukaryotic concern
    if organism:
        from ..organism_config import is_eukaryotic_organism
        if not is_eukaryotic_organism(organism):
            return sequence, []

    n_codons = len(protein)
    if n_codons == 0:
        return sequence, []

    # Precompute sorted codons per amino acid (by CAI descending)
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(protein):
        if aa == "*":
            continue
        codons = AA_TO_CODONS.get(aa, [])
        sorted_codons[aa] = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)

    # Precompute restriction site sequences to avoid
    rs_sites: list[tuple[str, str]] = []
    max_site_len = 0
    if enzymes:
        from ..restriction_sites import get_recognition_site as _get_site
        for enz in enzymes:
            site = _get_site(enz)
            if site is not None:
                site_rc = reverse_complement(site)
                rs_sites.append((site, site_rc))
                max_site_len = max(max_site_len, len(site))

    # ── Cost-aware CpG elimination ──
    # For eukaryotes, CpG avoidance is a SOFT preference, not a hard constraint.
    # Only eliminate CG dinucleotides when they contribute to a CpG island
    # (Obs/Exp ratio > cpg_threshold).  If the sequence already passes the
    # CpG island check, do NOT eliminate individual CGs — the CAI cost is
    # too high.  Individual CG dinucleotides in a CDS are common in
    # high-expression genes and are not biologically problematic unless they
    # cluster into CpG islands.
    from ..type_system import check_no_cpg_island

    # Check if the sequence already passes the CpG island check
    cpg_result = check_no_cpg_island(sequence, cpg_window, cpg_threshold)
    
    # Count total CG dinucleotides — even if the sequence passes the CpG island
    # ratio check, we should still eliminate individual CGs if they're present.
    # Short sequences (< cpg_window) may pass the island check despite having
    # many CG dinucleotides because the windowed scan doesn't apply.
    total_cg = sum(1 for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG")
    
    if cpg_result.passed and total_cg == 0:
        # No CpG island and no individual CGs — nothing to do
        return sequence, []

    # ── Adaptive CAI cost threshold ──
    # When the sequence FAILS the CpG island check, we need to eliminate CGs
    # more aggressively (higher max_cai_cost) to break the CpG island.
    # When it PASSES, CpG avoidance is a soft preference and we use the
    # user-specified (or default low) max_cai_cost to protect CAI.
    # For short sequences (< cpg_window), the island check trivially passes,
    # but a high CG density is still biologically relevant.  In this case,
    # use a more permissive threshold proportional to the CG density.
    effective_max_cai_cost = max_cai_cost
    if not cpg_result.passed:
        # CpG island detected — allow more CAI sacrifice to fix it.
        # Use 0.20 as the adaptive ceiling (enough to swap CGT→AGA, etc.)
        # but still protect against the worst swaps.
        effective_max_cai_cost = max(max_cai_cost, 0.20)
        logger.debug(
            "CpG island detected (Obs/Exp > %.2f): using adaptive "
            "max_cai_cost=%.4f (up from %.4f) to break the island",
            cpg_threshold, effective_max_cai_cost, max_cai_cost,
        )
    elif cpg_result.passed and total_cg > 0 and len(sequence) < cpg_window:
        # Short sequence: CpG island check doesn't apply but CGs are present.
        # Use a density-adaptive threshold: if >30% of dinucleotides are CG,
        # treat it as aggressively as a CpG island.
        cg_density = total_cg / max(len(sequence) - 1, 1)
        if cg_density > 0.30:
            effective_max_cai_cost = max(max_cai_cost, 0.20)
            logger.debug(
                "Short sequence (%d bp < window %d) with high CG density "
                "(%.1f%%): using adaptive max_cai_cost=%.4f",
                len(sequence), cpg_window, cg_density * 100,
                effective_max_cai_cost,
            )

    if cpg_result.passed and total_cg > 0:
        # Sequence passes island check but has individual CGs — still attempt
        # elimination as a best-effort pass (soft preference)
        logger.debug(
            "Sequence passes CpG island check but has %d CG dinucleotides — "
            "attempting elimination as soft preference (max_cai_cost=%.4f)",
            total_cg, effective_max_cai_cost,
        )

    warnings: list[str] = []
    seq = sequence
    eliminated_count = 0
    skipped_cai_cost = 0  # Track CpG fixes skipped due to CAI cost
    initial_cg_count = _count_dinucs_fast(seq, "CG")[0]
    # Track the initial GC content to allow swaps that move GC toward the
    # target range even if the current GC is outside the range.
    initial_gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)

    # PERF (Fix B): Track GC count incrementally instead of recounting the
    # entire sequence on every candidate swap.  This avoids O(n) GC counting
    # in the inner loops of CpG elimination.
    seq_gc_count = seq.count("G") + seq.count("C")
    seq_len = len(seq)

    def _gc_ok_incremental(current_gc_count: int, alt: str, current_codon: str) -> bool:
        """Check if a swap is acceptable from a GC perspective (incremental).

        Computes the new GC count by adjusting for the codon swap instead of
        recounting the entire sequence.

        Args:
            current_gc_count: Current G+C count in the sequence.
            alt: The replacement codon.
            current_codon: The codon being replaced.
        """
        if gc_lo <= 0.0 and gc_hi >= 1.0:
            return True
        alt_gc = sum(1 for b in alt if b in "GC")
        cur_gc = sum(1 for b in current_codon if b in "GC")
        new_gc_count = current_gc_count - cur_gc + alt_gc
        test_gc = new_gc_count / max(seq_len, 1)
        # Case 1: Test GC is within range — always OK
        if gc_lo <= test_gc <= gc_hi:
            return True
        # Case 2: Current GC is already outside range — accept if test GC
        # is closer to the target range (moves in the right direction)
        current_gc = current_gc_count / max(seq_len, 1)
        if not (gc_lo <= current_gc <= gc_hi):
            def _dist_to_range(gc_val: float) -> float:
                if gc_val < gc_lo:
                    return gc_lo - gc_val
                elif gc_val > gc_hi:
                    return gc_val - gc_hi
                return 0.0
            return _dist_to_range(test_gc) <= _dist_to_range(current_gc)
        # Case 3: Current GC is within range but test GC is not — reject
        return False

    def _gc_ok_two_codon_swap(current_gc_count: int,
                               left_alt: str, left_current: str,
                               right_alt: str, right_current: str) -> bool:
        """Check if a two-codon swap is acceptable from a GC perspective."""
        if gc_lo <= 0.0 and gc_hi >= 1.0:
            return True
        gc_delta = ((sum(1 for b in left_alt if b in "GC") - sum(1 for b in left_current if b in "GC")) +
                    (sum(1 for b in right_alt if b in "GC") - sum(1 for b in right_current if b in "GC")))
        new_gc_count = current_gc_count + gc_delta
        test_gc = new_gc_count / max(seq_len, 1)
        if gc_lo <= test_gc <= gc_hi:
            return True
        current_gc = current_gc_count / max(seq_len, 1)
        if not (gc_lo <= current_gc <= gc_hi):
            def _dist_to_range(gc_val: float) -> float:
                if gc_val < gc_lo:
                    return gc_lo - gc_val
                elif gc_val > gc_hi:
                    return gc_val - gc_hi
                return 0.0
            return _dist_to_range(test_gc) <= _dist_to_range(current_gc)
        return False

    def _creates_boundary_cg(test_seq: str, codon_idx: int) -> bool:
        """Check if modifying codon at codon_idx creates a new boundary CG.

        PERF (Fix B): Uses only local bases around the codon boundary.
        Accepts either a full sequence or a local region (when
        local_offset is provided).
        """
        codon_start = codon_idx * 3
        # Check boundary with previous codon
        if codon_idx > 0 and codon_start > 0:
            if test_seq[codon_start - 1:codon_start + 1] == "CG":
                return True
        # Check boundary with next codon
        next_start = codon_start + 3
        if next_start < len(test_seq):
            if test_seq[next_start - 1:next_start + 1] == "CG":
                return True
        return False

    def _creates_boundary_cg_local(local_region: str, local_offset: int,
                                    codon_idx: int, region_start: int) -> bool:
        """Check boundary CG using a local region instead of full sequence.

        Args:
            local_region: The local substring of the sequence.
            local_offset: Position of codon_idx*3 relative to region_start
                          (i.e., codon_idx*3 - region_start).
            codon_idx: The codon index being modified.
            region_start: The start position of local_region in the full sequence.
        """
        codon_start = local_offset  # position of codon in local_region
        # Check boundary with previous codon
        if codon_idx > 0 and codon_start > 0:
            if local_region[codon_start - 1:codon_start + 1] == "CG":
                return True
        # Check boundary with next codon
        next_start = codon_start + 3
        if next_start + 1 <= len(local_region):
            if local_region[next_start - 1:next_start + 1] == "CG":
                return True
        return False

    def _check_rs(test_seq: str, codon_start: int) -> bool:
        """Check if test_seq has restriction sites near codon_start.
        Returns True if OK (no sites), False if a site was created."""
        if not rs_sites:
            return True
        check_start = max(0, codon_start - max_site_len + 1)
        check_end = min(len(test_seq), codon_start + 3 + max_site_len - 1)
        local_region = test_seq[check_start:check_end]
        for site, site_rc in rs_sites:
            if site in local_region or (site_rc and site_rc in local_region):
                return False
        return True

    for iteration in range(max_iterations):
        # Find all CG dinucleotide positions
        cpg_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "CG"]
        if not cpg_positions:
            break  # All CGs eliminated

        any_fixed = False

        for pos in cpg_positions:
            # Re-check: this position may have been fixed by an earlier swap
            if seq[pos:pos+2] != "CG":
                continue

            left_ci = pos // 3          # codon containing the C
            right_ci = (pos + 1) // 3   # codon containing the G
            is_cross_codon = (left_ci != right_ci)

            fixed = False

            if not is_cross_codon:
                # ── Within-codon CpG ──
                # The CG is entirely within one codon. Replace with a
                # synonymous codon that doesn't contain "CG".
                if left_ci >= n_codons:
                    continue
                aa = protein[left_ci]
                if aa == "*":
                    continue
                current = seq[left_ci*3:left_ci*3+3]

                # PERF (Fix B): Build only the local region for restriction
                # site and boundary-CG checks instead of the full sequence.
                # A CpG within a codon only affects that codon ±1 neighbor.
                _buf = max(max_site_len, 4)  # enough for RS and boundary CG
                _region_start = max(0, left_ci * 3 - _buf + 1)
                _region_end = min(seq_len, left_ci * 3 + 3 + _buf - 1)
                _current_region = seq[_region_start:_region_end]
                _local_codon_start = left_ci * 3 - _region_start

                # Try alternatives sorted by CAI (best first)
                best_alt: str | None = None
                best_alt_cai = -1.0
                for alt in sorted_codons.get(aa, []):
                    if alt == current or "CG" in alt:
                        continue

                    # Build local test region by splicing alt into current region
                    _test_region = (_current_region[:_local_codon_start] +
                                    alt +
                                    _current_region[_local_codon_start + 3:])

                    # Check restriction sites in local region
                    if rs_sites:
                        site_found = False
                        for site, site_rc in rs_sites:
                            if site in _test_region or (site_rc and site_rc in _test_region):
                                site_found = True
                                break
                        if site_found:
                            continue

                    # Check GC content (incremental)
                    if not _gc_ok_incremental(seq_gc_count, alt, current):
                        continue

                    # Check that we don't create new boundary CGs (local)
                    if _creates_boundary_cg_local(_test_region, _local_codon_start, left_ci, _region_start):
                        continue

                    alt_cai = usage.get(alt, 0.0)
                    if alt_cai > best_alt_cai:
                        best_alt = alt
                        best_alt_cai = alt_cai

                if best_alt is not None:
                    # CAI cost gate: if the best CG-free alternative has
                    # too much CAI loss, keep the CG (CAI > CpG avoidance).
                    current_cai = usage.get(current, 0.0)
                    cai_loss = current_cai - best_alt_cai
                    if cai_loss > effective_max_cai_cost:
                        skipped_cai_cost += 1
                        logger.debug(
                            "CpG elimination skipped at codon %d (%s→%s): "
                            "CAI cost %.4f > max_cai_cost %.4f",
                            left_ci, current, best_alt,
                            cai_loss, effective_max_cai_cost,
                        )
                    else:
                        # Update GC count incrementally
                        _old_gc = sum(1 for b in current if b in "GC")
                        _new_gc = sum(1 for b in best_alt if b in "GC")
                        seq_gc_count += _new_gc - _old_gc
                        seq = seq[:left_ci*3] + best_alt + seq[left_ci*3+3:]
                        fixed = True
                        any_fixed = True
                        eliminated_count += 1
            else:
                # ── Cross-codon CpG ──
                # The C is the last base of codon left_ci, the G is the
                # first base of codon right_ci. Try:
                # 1. Change right codon to not start with G
                # 2. Change left codon to not end with C
                # 3. Coordinated 2-codon swap
                # For each strategy, pick the highest-CAI alternative.

                # Strategy 1: Change right codon (to not start with G)
                if 0 <= right_ci < n_codons:
                    right_aa = protein[right_ci]
                    if right_aa != "*":
                        right_current = seq[right_ci*3:right_ci*3+3]

                        # PERF (Fix B): Local region for right codon checks
                        _buf = max(max_site_len, 4)
                        _r_region_start = max(0, right_ci * 3 - _buf + 1)
                        _r_region_end = min(seq_len, right_ci * 3 + 3 + _buf - 1)
                        _r_current_region = seq[_r_region_start:_r_region_end]
                        _r_local_start = right_ci * 3 - _r_region_start

                        best_right_alt: str | None = None
                        best_right_cai = -1.0
                        for alt in sorted_codons.get(right_aa, []):
                            if alt == right_current or alt[0] == 'G':
                                continue

                            # Build local test region
                            _test_region = (_r_current_region[:_r_local_start] +
                                            alt +
                                            _r_current_region[_r_local_start + 3:])

                            # Check restriction sites in local region
                            if rs_sites:
                                site_found = False
                                for site, site_rc in rs_sites:
                                    if site in _test_region or (site_rc and site_rc in _test_region):
                                        site_found = True
                                        break
                                if site_found:
                                    continue

                            # Check GC content (incremental)
                            if not _gc_ok_incremental(seq_gc_count, alt, right_current):
                                continue

                            # Check no new boundary CGs (local)
                            if _creates_boundary_cg_local(_test_region, _r_local_start, right_ci, _r_region_start):
                                continue

                            # Make sure the specific CG at pos is gone.
                            # For cross-codon CG at pos: C is last base of
                            # left_ci, G is first base of right_ci.
                            # With alt swapped in at right_ci, the base at
                            # pos+1 is alt[0] which is not G (already checked).
                            # So the CG at pos is definitely gone.
                            # But verify alt doesn't create CG at pos-1:pos+1
                            # (the new right codon's first base with left
                            # codon's last base).
                            if pos > 0 and _r_local_start > 0:
                                _new_boundary = seq[pos] + alt[0]
                                if _new_boundary == "CG":
                                    continue

                            alt_cai = usage.get(alt, 0.0)
                            if alt_cai > best_right_cai:
                                best_right_alt = alt
                                best_right_cai = alt_cai

                        if best_right_alt is not None:
                            # CAI cost gate for cross-codon CpG strategy 1
                            right_current_cai = usage.get(right_current, 0.0)
                            right_cai_loss = right_current_cai - best_right_cai
                            if right_cai_loss > effective_max_cai_cost:
                                skipped_cai_cost += 1
                                logger.debug(
                                    "CpG elimination skipped (strategy 1, right codon %d: %s→%s): "
                                    "CAI cost %.4f > max_cai_cost %.4f",
                                    right_ci, right_current, best_right_alt,
                                    right_cai_loss, effective_max_cai_cost,
                                )
                            else:
                                _old_gc = sum(1 for b in right_current if b in "GC")
                                _new_gc = sum(1 for b in best_right_alt if b in "GC")
                                seq_gc_count += _new_gc - _old_gc
                                seq = seq[:right_ci*3] + best_right_alt + seq[right_ci*3+3:]
                                fixed = True
                                any_fixed = True
                                eliminated_count += 1

                # Strategy 2: Change left codon (to not end with C)
                if not fixed and 0 <= left_ci < n_codons:
                    left_aa = protein[left_ci]
                    if left_aa != "*":
                        left_current = seq[left_ci*3:left_ci*3+3]

                        # PERF (Fix B): Local region for left codon checks
                        _buf = max(max_site_len, 4)
                        _l_region_start = max(0, left_ci * 3 - _buf + 1)
                        _l_region_end = min(seq_len, left_ci * 3 + 3 + _buf - 1)
                        _l_current_region = seq[_l_region_start:_l_region_end]
                        _l_local_start = left_ci * 3 - _l_region_start

                        best_left_alt: str | None = None
                        best_left_cai = -1.0
                        for alt in sorted_codons.get(left_aa, []):
                            if alt == left_current or alt[-1] == 'C':
                                continue

                            # Build local test region
                            _test_region = (_l_current_region[:_l_local_start] +
                                            alt +
                                            _l_current_region[_l_local_start + 3:])

                            # Check restriction sites in local region
                            if rs_sites:
                                site_found = False
                                for site, site_rc in rs_sites:
                                    if site in _test_region or (site_rc and site_rc in _test_region):
                                        site_found = True
                                        break
                                if site_found:
                                    continue

                            # Check GC content (incremental)
                            if not _gc_ok_incremental(seq_gc_count, alt, left_current):
                                continue

                            # Check no new boundary CGs (local)
                            if _creates_boundary_cg_local(_test_region, _l_local_start, left_ci, _l_region_start):
                                continue

                            # Make sure the specific CG at pos is gone.
                            # For cross-codon CG at pos: C is last base of
                            # left_ci (=alt[-1], which is not C, already checked),
                            # G is first base of right_ci. So CG at pos is gone.
                            # But verify alt doesn't create CG at pos:pos+2
                            # (alt's last base with right codon's first base).
                            if pos + 2 < seq_len:
                                _new_boundary = alt[-1] + seq[pos + 1]
                                if _new_boundary == "CG":
                                    continue

                            alt_cai = usage.get(alt, 0.0)
                            if alt_cai > best_left_cai:
                                best_left_alt = alt
                                best_left_cai = alt_cai

                        if best_left_alt is not None:
                            # CAI cost gate for cross-codon CpG strategy 2
                            left_current_cai = usage.get(left_current, 0.0)
                            left_cai_loss = left_current_cai - best_left_cai
                            if left_cai_loss > effective_max_cai_cost:
                                skipped_cai_cost += 1
                                logger.debug(
                                    "CpG elimination skipped (strategy 2, left codon %d: %s→%s): "
                                    "CAI cost %.4f > max_cai_cost %.4f",
                                    left_ci, left_current, best_left_alt,
                                    left_cai_loss, effective_max_cai_cost,
                                )
                            else:
                                _old_gc = sum(1 for b in left_current if b in "GC")
                                _new_gc = sum(1 for b in best_left_alt if b in "GC")
                                seq_gc_count += _new_gc - _old_gc
                                seq = seq[:left_ci*3] + best_left_alt + seq[left_ci*3+3:]
                                fixed = True
                                any_fixed = True
                                eliminated_count += 1

                # Strategy 3: Coordinated 2-codon swap
                # PERF (Fix B): The original code built full test sequences
                # via "".join(test_list) inside the nested product loop over
                # left×right codon alternatives — O(n) per combination.
                # Fix: Only check the local region around the CpG site.
                # A CpG dinucleotide spans at most 2 adjacent codons.
                # Rebuild only those codons + a small buffer for restriction
                # site checking, not the entire sequence.
                if (not fixed and 0 <= left_ci < n_codons
                        and 0 <= right_ci < n_codons):
                    left_aa = protein[left_ci]
                    right_aa = protein[right_ci]
                    if left_aa != "*" and right_aa != "*":
                        left_current = seq[left_ci*3:left_ci*3+3]
                        right_current = seq[right_ci*3:right_ci*3+3]

                        # Precompute local region bounds for both codons
                        _buf = max(max_site_len, 4)
                        _s3_region_start = max(0, left_ci * 3 - _buf + 1)
                        _s3_region_end = min(seq_len, right_ci * 3 + 3 + _buf - 1)
                        _s3_current_region = seq[_s3_region_start:_s3_region_end]
                        _s3_local_left = left_ci * 3 - _s3_region_start
                        _s3_local_right = right_ci * 3 - _s3_region_start
                        _s3_local_pos = pos - _s3_region_start

                        # Try all combinations of left+right alternatives
                        best_swap: tuple[str, str, float] | None = None  # (left_alt, right_alt, cai_sum)
                        for left_alt in sorted_codons.get(left_aa, []):
                            if left_alt == left_current:
                                continue
                            for right_alt in sorted_codons.get(right_aa, []):
                                if right_alt == right_current and left_alt == left_current:
                                    continue
                                # Skip if boundary still has CG
                                if left_alt[-1] == 'C' and right_alt[0] == 'G':
                                    continue
                                # Skip if within-codon CG created
                                if "CG" in left_alt or "CG" in right_alt:
                                    continue

                                # Build local test region by splicing both alts
                                if left_ci == right_ci - 1:
                                    # Adjacent codons — single splice
                                    _test_region = (_s3_current_region[:_s3_local_left] +
                                                    left_alt + right_alt +
                                                    _s3_current_region[_s3_local_right + 3:])
                                else:
                                    # Non-adjacent — two separate splices
                                    _test_region = (_s3_current_region[:_s3_local_left] +
                                                    left_alt +
                                                    _s3_current_region[_s3_local_left + 3:_s3_local_right] +
                                                    right_alt +
                                                    _s3_current_region[_s3_local_right + 3:])

                                # Check restriction sites in local region
                                site_ok = True
                                for site, site_rc in rs_sites:
                                    if site in _test_region or (site_rc and site_rc in _test_region):
                                        site_ok = False
                                        break
                                if not site_ok:
                                    continue

                                # Check GC content (incremental for two-codon swap)
                                if not _gc_ok_two_codon_swap(seq_gc_count, left_alt, left_current, right_alt, right_current):
                                    continue

                                # Check no new boundary CGs (local)
                                if _creates_boundary_cg_local(_test_region, _s3_local_left, left_ci, _s3_region_start):
                                    continue
                                if _creates_boundary_cg_local(_test_region, _s3_local_right, right_ci, _s3_region_start):
                                    continue

                                # Verify the specific CG at pos is gone
                                if 0 <= _s3_local_pos < len(_test_region) - 1:
                                    if _test_region[_s3_local_pos:_s3_local_pos + 2] == "CG":
                                        continue
                                else:
                                    # Position out of local region — check directly
                                    # Since we changed both codons around the CG,
                                    # and we already checked left_alt[-1]!='C' or
                                    # right_alt[0]!='G', the CG should be gone.
                                    pass

                                cai_sum = (usage.get(left_alt, 0.0) +
                                           usage.get(right_alt, 0.0))
                                if best_swap is None or cai_sum > best_swap[2]:
                                    best_swap = (left_alt, right_alt, cai_sum)

                        if best_swap is not None:
                            left_alt, right_alt, _ = best_swap
                            # CAI cost gate for cross-codon CpG strategy 3
                            left_current_cai_s3 = usage.get(left_current, 0.0)
                            right_current_cai_s3 = usage.get(right_current, 0.0)
                            left_alt_cai_s3 = usage.get(left_alt, 0.0)
                            right_alt_cai_s3 = usage.get(right_alt, 0.0)
                            cai_loss_l = left_current_cai_s3 - left_alt_cai_s3
                            cai_loss_r = right_current_cai_s3 - right_alt_cai_s3
                            if cai_loss_l > effective_max_cai_cost or cai_loss_r > effective_max_cai_cost:
                                skipped_cai_cost += 1
                                logger.debug(
                                    "CpG elimination skipped (strategy 3, codons %d+%d: %s→%s, %s→%s): "
                                    "CAI cost L=%.4f R=%.4f > max_cai_cost %.4f",
                                    left_ci, right_ci,
                                    left_current, left_alt,
                                    right_current, right_alt,
                                    cai_loss_l, cai_loss_r, effective_max_cai_cost,
                                )
                            else:
                                # Update GC count incrementally
                                _old_gc_l = sum(1 for b in left_current if b in "GC")
                                _new_gc_l = sum(1 for b in left_alt if b in "GC")
                                _old_gc_r = sum(1 for b in right_current if b in "GC")
                                _new_gc_r = sum(1 for b in right_alt if b in "GC")
                                seq_gc_count += (_new_gc_l - _old_gc_l) + (_new_gc_r - _old_gc_r)
                                seq = (seq[:left_ci*3] + left_alt +
                                       seq[left_ci*3+3:right_ci*3] + right_alt +
                                       seq[right_ci*3+3:])
                                fixed = True
                                any_fixed = True
                                eliminated_count += 1

        # If no CG was fixed in this pass, stop trying
        if not any_fixed:
            break

        # Do NOT break early when CpG island check passes — we want to
        # eliminate ALL CG dinucleotides, not just enough to pass the
        # Obs/Exp ratio check.

    # Report remaining CpGs
    remaining_cpgs = [i for i in range(len(seq) - 1) if seq[i:i+2] == "CG"]
    if remaining_cpgs:
        # Only warn if the sequence still fails the CpG island check
        final_cpg_result = check_no_cpg_island(seq, cpg_window, cpg_threshold, organism=organism)
        if not final_cpg_result.passed:
            warnings.append(
                f"CpG island avoidance: {len(remaining_cpgs)} CG dinucleotide(s) remain. "
                f"No synonymous substitution could eliminate them without creating "
                f"restriction sites or other violations. "
                f"{final_cpg_result.details}"
            )

    # Log summary including CAI cost skips
    if skipped_cai_cost > 0:
        logger.warning(
            "CpG elimination: %d position(s) skipped because CAI cost exceeded "
            "max_cai_cost=%.4f. CAI is preserved at the cost of retaining CG "
            "dinucleotides. For eukaryotes, CpG avoidance is a soft preference — "
            "individual CGs in CDS are common and not biologically problematic "
            "unless they form CpG islands.",
            skipped_cai_cost, effective_max_cai_cost,
        )

    logger.debug(
        "CpG elimination: %d/%d CG dinucleotides eliminated in %d passes "
        "(%d skipped due to CAI cost)",
        eliminated_count, initial_cg_count, iteration + 1 if max_iterations > 0 else 0,
        skipped_cai_cost,
    )

    return seq, warnings


# ==============================================================================
# Greedy Optimizer
# ==============================================================================



def _greedy_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    cryptic_splice_threshold: float = 3.0,
    seed: int | None = None,
    provenance_collector: DecisionProvenanceCollector | None = None,
    is_prokaryote: bool = False,
) -> tuple[str, list[str]]:
    """
    Greedy multi-objective codon optimization with coordinated constraint solving.

    Step ordering prioritizes hard constraints (restriction sites) over soft constraints (CAI).
    Reconciliation pass ensures earlier steps aren't undone by later ones.

    Steps:
    1. Best codon per position (maximize CAI)
    2. Remove restriction sites (multi-codon coordinated)
    3. Remove ATTTA instability motifs
    4. Fix 6+ consecutive T runs
    5. Adjust GC content (hard constraint, organism target aspiration)
    6. Reconciliation — restriction sites vs GC
    7. Eliminate cryptic splice donor/acceptor sites (SKIPPED for prokaryotes)
    7.5. Disrupt CpG dinucleotides to avoid CpG islands (SKIPPED for prokaryotes)
    8. Reconciliation — restriction sites after splice/CpG fixes (SKIPPED for prokaryotes)
    8.5. CpG reconciliation (SKIPPED for prokaryotes)

    Pre-conditions:
    - protein is a valid amino acid sequence (no invalid codes)
    - organism is in SUPPORTED_ORGANISMS
    - 0.0 <= gc_lo < gc_hi <= 1.0
    - cryptic_splice_threshold > 0

    Post-conditions:
    - returned sequence translates to the input protein
    - len(returned sequence) == len(protein) * 3
    - all codons in sequence are valid for their amino acid

    Args:
        is_prokaryote: When True, skip eukaryote-specific constraint steps
            (cryptic splice elimination, CpG disruption, and their
            reconciliation passes). Prokaryotes have no spliceosome, so
            GT/AG avoidance and CpG island disruption are biologically
            irrelevant and unnecessarily lower CAI.

    Note: The ``seed`` parameter is currently unused because the greedy
    optimizer is fully deterministic. It is reserved for future
    randomized optimization steps.
    """
    # Set deterministic seed if provided (reserved for future randomized steps).
    # Use a per-call random.Random instance instead of random.seed() to avoid
    # polluting global random state (thread-safety / reproducibility concern).
    _rng = None
    if seed is not None:
        import random as _random_mod
        _rng = _random_mod.Random(seed)

    # Validate pre-conditions
    if not (0.0 <= gc_lo < gc_hi <= 1.0):
        raise ValueError(f"Invalid GC range: gc_lo={gc_lo}, gc_hi={gc_hi}")
    if cryptic_splice_threshold <= 0:
        raise ValueError(f"Threshold must be positive, got {cryptic_splice_threshold}")

    usage = CODON_ADAPTIVENESS_TABLES.get(organism)
    if usage is None:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)
    aas = protein_to_aa_list(protein)
    restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.values())
    warnings: list[str] = []

    # sorted_codons: codons for each AA sorted by CAI (descending) — used in ALL
    # constraint resolution steps so that the highest-CAI alternative that
    # fixes the constraint is always preferred.
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(aas):
        codons = AA_TO_CODONS[aa]
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted

    # Step: Maximize CAI — Best codon per position (maximize CAI)
    # For eukaryotes, use GT-aware codon selection that avoids boundary GTs
    # when there's an alternative within 10% relative CAI.
    if provenance_collector is not None:
        # Expanded loop with per-codon provenance tracking
        chosen_codons: list[str] = []
        for pos, aa in enumerate(aas):
            # For eukaryotes, use GT-aware codon selection to avoid
            # boundary GTs when possible without destroying CAI
            if not is_prokaryote:
                next_aa = aas[pos + 1] if pos + 1 < len(aas) else None
                chosen = _gt_aware_select_codon(aa, next_aa, sorted_codons, usage)
            else:
                candidates = sorted_codons[aa]
                chosen = candidates[0]

            # Build alternatives list for provenance (excluding chosen codon)
            alternatives: list[dict[str, Any]] = []
            chosen_cai = usage.get(chosen, 0.0)
            for codon in candidates:
                if codon == chosen:
                    continue  # chosen codon is already in chosen_codon field
                cai_val = usage.get(codon, 0.0)
                gc_bases = sum(1 for b in codon if b in "GC")
                gc_contribution = gc_bases / 3.0
                # Check constraint violations for this codon
                violates: list[str] = []
                if "GT" in codon:
                    violates.append("cryptic_splice_donor")
                if "AG" in codon:
                    violates.append("cryptic_splice_acceptor")
                gc_bases_total = sum(1 for b in codon if b in "GC")
                if gc_bases_total == 0 and gc_lo > 0.5:
                    violates.append("gc_too_low")
                elif gc_bases_total == 3 and gc_hi < 0.5:
                    violates.append("gc_too_high")
                # Determine rejection reason
                rejected_because: str | None = None
                if violates:
                    rejected_because = f"Violates: {', '.join(violates)}"
                elif cai_val < chosen_cai:
                    rejected_because = "Lower CAI"
                else:
                    rejected_because = "Lower CAI"
                alternatives.append({
                    "codon": codon,
                    "cai_contribution": round(cai_val, 4),
                    "gc_contribution": round(gc_contribution, 2),
                    "violates_constraints": violates,
                    "rejected_because": rejected_because,
                })

            # Compute confidence: 1.0 if the best codon is clearly better,
            # lower if alternatives are close in CAI
            if len(candidates) > 1:
                second_best_cai = usage.get(candidates[1], 0.0)
                confidence = min(1.0, 0.5 + (chosen_cai - second_best_cai) * 5)
                confidence = max(0.0, confidence)
            else:
                confidence = 1.0

            provenance_collector.record_codon_decision(CodonDecision(
                position=pos,
                amino_acid=aa,
                original_codon=None,
                chosen_codon=chosen,
                alternatives_considered=alternatives,
                constraint_reason="Maximize CAI while maintaining GC in range",
                confidence=round(confidence, 4),
            ))
            chosen_codons.append(chosen)
        sequence = "".join(chosen_codons)
    else:
        # Original fast path (no provenance overhead)
        # For eukaryotes, use GT-aware codon selection
        if not is_prokaryote:
            chosen_codons_fast: list[str] = []
            for pos, aa in enumerate(aas):
                next_aa = aas[pos + 1] if pos + 1 < len(aas) else None
                chosen_codons_fast.append(
                    _gt_aware_select_codon(aa, next_aa, sorted_codons, usage)
                )
            sequence = "".join(chosen_codons_fast)
        else:
            sequence = "".join(sorted_codons[aa][0] for aa in aas)
    if len(sequence) != len(aas) * 3:
        raise ValueError("Maximize CAI step: sequence length mismatch")

    # Step: Remove Restriction Sites (HIGHEST PRIORITY — multi-codon coordinated)
    # Process concrete sites first, then IUPAC sites
    concrete_sites = []
    iupac_sites = []
    for site in restriction_sites:
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            iupac_sites.append(site_upper)
        else:
            concrete_sites.append(site_upper)

    # Build Aho-Corasick scanner for O(L+M) multi-pattern detection of all
    # concrete sites + their reverse complements simultaneously.
    # This replaces per-site O(N*L*site_len) scanning with a single O(L+M) pass.
    concrete_scanner = build_scanner_from_sites(concrete_sites) if concrete_sites else None

    # Remove concrete sites using Aho-Corasick for fast detection
    if concrete_scanner is not None:
        # Fast path: scan for ALL sites at once, fix one at a time
        for iteration in range(MAX_RESTRICTION_SITE_ITERATIONS * len(concrete_sites)):
            matches = concrete_scanner.scan(sequence)
            if not matches:
                break

            # Fix the first match found
            pos, site_match, _label = matches[0]
            site_len = len(site_match)
            site_rc = reverse_complement(site_match)

            # Try multi-codon coordinated removal (CAI-aware)
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_match, site_rc, usage=usage
            )
            if fixed:
                sequence = new_seq
                continue

            # Fallback: try single-codon swap — CAI-aware: find the best
            # swap across ALL overlapping codon positions, not just the first
            overlapping = _get_overlapping_codons(pos, site_len, len(aas))
            best_single_swap: tuple[int, str, float] | None = None  # (ci, alt, log_cai)
            for ci in overlapping:
                if ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci * 3: ci * 3 + 3]
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    seq_list = list(sequence)
                    seq_list[ci * 3: ci * 3 + 3] = list(alt)
                    test = "".join(seq_list)
                    if site_match not in test and site_rc not in test:
                        alt_cai = usage.get(alt, 0.0)
                        log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                        if best_single_swap is None or log_cai > best_single_swap[2]:
                            best_single_swap = (ci, alt, log_cai)
            if best_single_swap is not None:
                ci, alt, _ = best_single_swap
                seq_list = list(sequence)
                seq_list[ci * 3: ci * 3 + 3] = list(alt)
                sequence = "".join(seq_list)
                continue

            # Try neighboring codons
            overlapping = _get_overlapping_codons(pos, site_len, len(aas))
            neighbor_fixed = False
            for ci in overlapping:
                if ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci * 3: ci * 3 + 3]
                for alt in sorted_codons[aa]:
                    if alt != current:
                        seq_list = list(sequence)
                        seq_list[ci * 3: ci * 3 + 3] = list(alt)
                        test = "".join(seq_list)
                        if site_match not in test and site_rc not in test:
                            sequence = test
                            neighbor_fixed = True
                            break
                if neighbor_fixed:
                    break
            if neighbor_fixed:
                continue

            warnings.append(f"Cannot remove {site_match} at iteration {iteration}")
            break
        else:
            # Check if any site is still present
            remaining = concrete_scanner.scan(sequence)
            if remaining:
                remaining_sites = [s for _, s, _ in remaining]
                logger.warning(
                    "Restriction site elimination did not converge after %d iterations. "
                    "%d restriction sites remain: %s",
                    MAX_RESTRICTION_SITE_ITERATIONS, len(remaining_sites),
                    remaining_sites[:5],
                )
                for _, site_still, _ in remaining[:3]:
                    warnings.append(
                        f"Restriction site elimination capped at {MAX_RESTRICTION_SITE_ITERATIONS} iterations. "
                        f"Site {site_still} could not be eliminated."
                    )
    else:
        # Fallback: per-site scan (no Aho-Corasick scanner available)
        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            for iteration in range(MAX_RESTRICTION_SITE_ITERATIONS):
                positions = _find_site_in_sequence(sequence, site_upper, site_rc)
                if not positions:
                    break

                # Try multi-codon coordinated removal (CAI-aware)
                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    continue

                # Fallback: CAI-aware single-codon swap across all overlapping
                # positions — find the swap with minimal CAI loss
                pos = positions[0]
                overlapping = _get_overlapping_codons(pos, len(site_upper), len(aas))
                best_single_swap: tuple[int, str, float] | None = None
                for ci in overlapping:
                    if ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * 3: ci * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        seq_list = list(sequence)
                        seq_list[ci * 3: ci * 3 + 3] = list(alt)
                        test = "".join(seq_list)
                        if site_upper not in test and site_rc not in test:
                            alt_cai = usage.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_single_swap is None or log_cai > best_single_swap[2]:
                                best_single_swap = (ci, alt, log_cai)
                if best_single_swap is not None:
                    ci, alt, _ = best_single_swap
                    seq_list = list(sequence)
                    seq_list[ci * 3: ci * 3 + 3] = list(alt)
                    sequence = "".join(seq_list)
                    continue

                # Try neighboring codons
                overlapping = _get_overlapping_codons(pos, len(site_upper), len(aas))
                neighbor_fixed = False
                for ci in overlapping:
                    if ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * 3: ci * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt != current:
                            seq_list = list(sequence)
                            seq_list[ci * 3: ci * 3 + 3] = list(alt)
                            test = "".join(seq_list)
                            if site_upper not in test and site_rc not in test:
                                sequence = test
                                neighbor_fixed = True
                                break
                    if neighbor_fixed:
                        break
                if neighbor_fixed:
                    continue

                warnings.append(f"Cannot remove {site_upper} at iteration {iteration}")
                break
            else:
                # Check if site is still present
                if site_upper in sequence or site_rc in sequence:
                    logger.warning(
                        "Restriction site elimination did not converge after %d iterations. "
                        "Site %s remains.",
                        MAX_RESTRICTION_SITE_ITERATIONS, site_upper,
                    )
                    warnings.append(
                        f"Restriction site elimination capped at {MAX_RESTRICTION_SITE_ITERATIONS} iterations. "
                        f"Site {site_upper} could not be eliminated."
                    )

    # Remove IUPAC sites (expand to concrete variants, check each)
    for site_upper in iupac_sites:
        concrete_variants = _expand_iupac_site(site_upper)
        if not concrete_variants:
            continue
        for variant in concrete_variants:
            variant_rc = reverse_complement(variant)
            for iteration in range(MAX_IUPAC_SITE_ITERATIONS):
                positions = _find_site_in_sequence(sequence, variant, variant_rc)
                if not positions:
                    break

                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, variant, variant_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    continue

                # CAI-aware single-codon fallback for IUPAC sites
                pos = positions[0]
                overlapping = _get_overlapping_codons(pos, len(variant), len(aas))
                best_iupac_swap: tuple[int, str, float] | None = None
                for ci in overlapping:
                    if ci >= len(aas):
                        continue
                    aa = aas[ci]
                    current = sequence[ci * 3: ci * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt == current:
                            continue
                        test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
                        if variant not in test and variant_rc not in test:
                            alt_cai = usage.get(alt, 0.0)
                            log_cai = math.log(alt_cai) if alt_cai > 0 else float('-inf')
                            if best_iupac_swap is None or log_cai > best_iupac_swap[2]:
                                best_iupac_swap = (ci, alt, log_cai)
                if best_iupac_swap is not None:
                    ci, alt, _ = best_iupac_swap
                    sequence = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
                else:
                    warnings.append(f"Cannot remove IUPAC {site_upper} variant {variant} at iteration {iteration}")
                    break
            else:
                if variant in sequence or variant_rc in sequence:
                    logger.warning(
                        "IUPAC site elimination did not converge after %d iterations. "
                        "Site %s variant %s remains.",
                        MAX_IUPAC_SITE_ITERATIONS, site_upper, variant,
                    )
                    warnings.append(
                        f"IUPAC site elimination capped at {MAX_IUPAC_SITE_ITERATIONS} iterations. "
                        f"Variant {variant} of {site_upper} could not be eliminated."
                    )

    # Step: Remove ATTTA instability motifs
    # PERF (Optimization D): Use list mutation for codon swaps
    for iteration in range(MAX_ATTTA_MOTIF_ITERATIONS):
        pos = sequence.find("ATTTA")
        if pos == -1:
            break
        codon_idx = pos // 3
        fixed = False
        for ci in range(max(0, codon_idx - 1), min(len(aas), codon_idx + 2)):
            aa = aas[ci]
            current = sequence[ci * 3:ci * 3 + 3]
            for alt in sorted_codons[aa]:
                if alt != current:
                    seq_list = list(sequence)
                    seq_list[ci * 3:ci * 3 + 3] = list(alt)
                    test = "".join(seq_list)
                    if "ATTTA" not in test:
                        sequence = test
                        fixed = True
                        break
            if fixed:
                break
        if not fixed:
            warnings.append(f"ATTTA motif: cannot remove at iteration {iteration}")
            break
    else:
        _remaining_attta = sequence.count("ATTTA")
        if _remaining_attta > 0:
            logger.warning(
                "ATTTA motif elimination did not converge after %d iterations. "
                "%d ATTTA motifs remain.",
                MAX_ATTTA_MOTIF_ITERATIONS, _remaining_attta,
            )
            warnings.append(
                f"ATTTA motif elimination capped at {MAX_ATTTA_MOTIF_ITERATIONS} iterations. "
                f"{_remaining_attta} motifs could not be eliminated."
            )

    # Step: Fix 6+ consecutive T runs
    for iteration in range(MAX_T_RUN_ITERATIONS):
        max_run, max_pos = 0, -1
        i = 0
        while i < len(sequence):
            if sequence[i] == "T":
                j = i
                while j < len(sequence) and sequence[j] == "T":
                    j += 1
                if j - i > max_run:
                    max_run, max_pos = j - i, i
                i = j
            else:
                i += 1
        if max_run < T_RUN_LENGTH_THRESHOLD:
            break
        codon_idx = (max_pos + max_run // 2) // 3
        if codon_idx < len(aas):
            aa = aas[codon_idx]
            current = sequence[codon_idx * 3:codon_idx * 3 + 3]
            fixed = False
            for alt in sorted_codons[aa]:
                if alt != current:
                    seq_list = list(sequence)
                    seq_list[codon_idx * 3:codon_idx * 3 + 3] = list(alt)
                    test = "".join(seq_list)
                    if not any(test[i:i + T_RUN_LENGTH_THRESHOLD] == "T" * T_RUN_LENGTH_THRESHOLD for i in range(len(test) - 5)):
                        sequence = test
                        fixed = True
                        break
            if not fixed:
                warnings.append(f"Consecutive T run: cannot fix at iteration {iteration}")
                break
    else:
        _remaining_t_runs = sum(
            1 for i in range(len(sequence) - T_RUN_LENGTH_THRESHOLD + 1)
            if sequence[i:i + T_RUN_LENGTH_THRESHOLD] == "T" * T_RUN_LENGTH_THRESHOLD
        )
        if _remaining_t_runs > 0:
            logger.warning(
                "T-run elimination did not converge after %d iterations. "
                "%d runs of %d+ T remain.",
                MAX_T_RUN_ITERATIONS, _remaining_t_runs, T_RUN_LENGTH_THRESHOLD,
            )
            warnings.append(
                f"T-run elimination capped at {MAX_T_RUN_ITERATIONS} iterations. "
                f"{_remaining_t_runs} runs of {T_RUN_LENGTH_THRESHOLD}+ T could not be eliminated."
            )

    # Step: Adjust GC content
    # Strategy: GC must be in [gc_lo, gc_hi] (hard constraint).
    # If in range, we gently nudge toward organism target but NEVER at the
    # cost of significant CAI reduction. The organism GC target is aspirational,
    # not mandatory — a sequence with CAI=0.99 and GC=0.61 (slightly above
    # human's 0.41 target) is better than CAI=0.82 and GC=0.46.
    # PERF (Optimization F): Cache GC count for incremental updates
    n_bases = len(sequence)
    if _HAS_NUMBA:
        try:
            gc_count = _numba_count_gc(_seq_to_bytes(sequence))
        except Exception:
            logger.warning("NUMBA GC count failed, falling back to pure-Python", exc_info=True)
            gc_count = sum(1 for b in sequence if b in "GC")
    else:
        gc_count = sum(1 for b in sequence if b in "GC")
    gc_val = gc_count / n_bases
    organism_gc_range = ORGANISM_GC_TARGETS.get(organism, (gc_lo, gc_hi))
    organism_gc = (organism_gc_range[0] + organism_gc_range[1]) / 2.0
    target_gc = max(gc_lo, min(gc_hi, organism_gc))

    gc_out_of_range = not (gc_lo <= gc_val <= gc_hi)

    if gc_out_of_range:
        # Hard constraint: MUST get GC into range
        # (gc_count and n_bases already computed above)
        # Target the nearest bound
        if gc_val < gc_lo:
            phase_target = gc_lo
        else:
            phase_target = gc_hi

        # PERF (Fix C): Use a priority queue (heapq) of codon swap
        # candidates sorted by GC impact instead of scanning ALL codons
        # every iteration.  With up to 200 iterations, the O(200 × n × k)
        # scan becomes O(k × n × log(n×k)) initial build + O(I × k × log(n×k))
        # amortized updates, where I = iterations and k = avg alternatives/AA.
        import heapq

        # Build initial heap of all swap candidates.
        # Each entry: (diff_from_target, -cai_for_tiebreak, ci,
        #              current_codon_at_insert, alt_codon, gc_delta)
        # We use -cai so that higher CAI sorts first (heapq is min-heap).
        # current_codon_at_insert enables lazy staleness detection.
        _gc_heap: list[tuple[float, float, int, str, str, int]] = []
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
                heapq.heappush(_gc_heap, (diff, -alt_cai, ci, current, alt, alt_gc - current_gc))

        for iteration in range(MAX_GC_ADJUSTMENT_ITERATIONS):
            if gc_lo <= gc_val <= gc_hi:
                break

            best_alt = None
            best_ci = -1
            best_gc_delta = 0

            # Pop the best candidate, skipping stale entries
            while _gc_heap:
                diff, neg_cai, ci, codon_at_insert, alt, gc_delta = heapq.heappop(_gc_heap)
                # Check if this entry is stale (the codon at ci has changed)
                current_at_ci = sequence[ci * 3:ci * 3 + 3]
                if current_at_ci != codon_at_insert:
                    # Stale entry — re-add with current codon
                    current_gc_ci = sum(1 for b in current_at_ci if b in "GC")
                    alt_gc = sum(1 for b in alt if b in "GC")
                    new_gc_count = gc_count - current_gc_ci + alt_gc
                    new_frac = new_gc_count / n_bases
                    new_diff = abs(new_frac - phase_target)
                    alt_cai = usage.get(alt, 0.0)
                    heapq.heappush(_gc_heap, (new_diff, -alt_cai, ci, current_at_ci, alt, alt_gc - current_gc_ci))
                    continue
                # Valid entry — use it
                best_alt = alt
                best_ci = ci
                best_gc_delta = gc_delta
                break

            if best_alt is None:
                break
            # Apply the swap
            seq_list = list(sequence)
            seq_list[best_ci * 3: best_ci * 3 + 3] = list(best_alt)
            sequence = "".join(seq_list)
            gc_count += best_gc_delta
            gc_val = gc_count / n_bases

            # Re-add candidates for the swapped position with the new codon
            aa = aas[best_ci]
            new_current = sequence[best_ci * 3:best_ci * 3 + 3]
            new_current_gc = sum(1 for b in new_current if b in "GC")
            for alt in sorted_codons[aa]:
                if alt == new_current:
                    continue
                alt_gc = sum(1 for b in alt if b in "GC")
                new_gc_count = gc_count - new_current_gc + alt_gc
                new_frac = new_gc_count / n_bases
                diff = abs(new_frac - phase_target)
                alt_cai = usage.get(alt, 0.0)
                heapq.heappush(_gc_heap, (diff, -alt_cai, best_ci, new_current, alt, alt_gc - new_current_gc))
        else:
            logger.warning(
                "GC adjustment did not converge after %d iterations. "
                "Current GC=%.3f (target range [%.3f, %.3f])",
                MAX_GC_ADJUSTMENT_ITERATIONS, gc_val, gc_lo, gc_hi,
            )
            warnings.append(
                f"GC adjustment capped at {MAX_GC_ADJUSTMENT_ITERATIONS} iterations. "
                f"Current GC={gc_val:.3f} (target range [{gc_lo:.3f}, {gc_hi:.3f}])."
            )

    # Step: Reconciliation — check if GC adjustment reintroduced restriction sites
    # Use Aho-Corasick scanner for fast multi-site detection if available
    if concrete_scanner is not None:
        remaining_matches = concrete_scanner.scan(sequence)
        for pos, site_match, _label in remaining_matches:
            site_rc = reverse_complement(site_match)
            # Try one more round of multi-codon removal
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_match, site_rc, usage=usage
            )
            if fixed:
                sequence = new_seq
                # Re-check GC
                if _HAS_NUMBA:
                    try:
                        gc_count = _numba_count_gc(_seq_to_bytes(sequence))
                    except Exception:
                        gc_count = sum(1 for b in sequence if b in "GC")
                else:
                    gc_count = sum(1 for b in sequence if b in "GC")
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
                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq
                    # Re-check GC
                    # PERF (Optimization F): Update cached GC count
                    if _HAS_NUMBA:
                        try:
                            gc_count = _numba_count_gc(_seq_to_bytes(sequence))
                        except Exception:
                            gc_count = sum(1 for b in sequence if b in "GC")
                    else:
                        gc_count = sum(1 for b in sequence if b in "GC")
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
                                # PERF (Optimization D): Use list mutation
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

    # Step: Eliminate cryptic splice donor/acceptor sites
    # EUKARYOTE-ONLY: Prokaryotes have no spliceosome, so cryptic splice
    # sites are biologically irrelevant. Skipping this step recovers
    # significant CAI on prokaryotic targets.
    #
    # IMPORTANT for eukaryotes: GT avoidance is a SOFT preference, not a
    # hard constraint. In-codon GTs from optimal codons (GGT, TGT, GTT,
    # etc.) are acceptable because they are unavoidable for high CAI.
    # Cross-codon GTs are only eliminated when the CAI cost of doing so
    # is < EUKARYOTE_CAI_GT_COST_THRESHOLD (default 0.05). This ensures
    # that CAI is prioritized over GT avoidance for eukaryotes.
    if not is_prokaryote:
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

                    # ── Splice donor potential check ──
                    # Not all GT dinucleotides are equally dangerous.  A GT
                    # with low splice donor potential (< SPLICE_DONOR_POTENTIAL_THRESHOLD)
                    # is unlikely to function as a cryptic splice donor, even if
                    # its MaxEntScan donor score exceeds the threshold.  This can
                    # happen when the surrounding context doesn't support splicing
                    # (e.g., no downstream AG acceptor, no polypyrimidine tract).
                    # For such GTs, CAI should always win.
                    sdp = score_splice_donor_potential(sequence, gt_pos)
                    if sdp < SPLICE_DONOR_POTENTIAL_THRESHOLD:
                        # This GT has low splice donor potential — not dangerous.
                        # Accept it and move on (CAI > GT avoidance here).
                        continue

                    # ── Eukaryotic GT-vs-CAI tradeoff ──
                    # For eukaryotes, in-codon GTs from optimal codons are
                    # acceptable (biologically common in high-expression genes).
                    # Only eliminate GT if the CAI cost is < threshold.
                    is_in_codon = _is_in_codon_gt(sequence, gt_pos)
                    current_codon = sequence[codon_idx*3:codon_idx*3+3]
                    optimal_codon = sorted_codons[aa][0]

                    if is_in_codon:
                        # In-codon GT: acceptable if the current codon is optimal
                        # or if swapping to a GT-free codon would cost too much CAI
                        if current_codon == optimal_codon:
                            # In-codon GT from optimal codon (e.g., GGT for Gly,
                            # TGT for Cys, GTT for Val) — acceptable for eukaryotes
                            continue
                        # Check CAI cost of best GT-free alternative
                        current_w = usage.get(current_codon, 0.0)
                        gt_free = _find_gt_free_codons(aa)
                        if gt_free:
                            best_gt_free_w = max(usage.get(c, 0.0) for c in gt_free)
                            if current_w - best_gt_free_w > EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                # CAI cost too high — keep the GT-containing codon
                                continue
                        else:
                            # No GT-free alternative (e.g., Valine) — must accept
                            continue
                    else:
                        # Cross-codon GT: for eukaryotes, only eliminate if the
                        # CAI cost of the best fix is < threshold. Check the CAI
                        # cost of changing the codon at gt_pos to a non-T-starting
                        # alternative (for the T-side) or non-G-ending alternative
                        # (for the G-side).
                        # The G is at the end of codon_idx, the T is at the start
                        # of codon_idx+1.
                        next_codon_idx = (gt_pos + 1) // 3
                        # Check CAI cost of changing either codon
                        _cross_cai_cost_ok = False
                        # Try changing the G-ending codon
                        if codon_idx < len(aas):
                            _g_aa = aas[codon_idx]
                            _g_current = sequence[codon_idx*3:codon_idx*3+3]
                            _g_current_w = usage.get(_g_current, 0.0)
                            _g_non_g_end = [c for c in sorted_codons[_g_aa] if c[-1] != "G"]
                            if _g_non_g_end:
                                _best_non_g_end_w = usage.get(_g_non_g_end[0], 0.0)
                                if _g_current_w - _best_non_g_end_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                    _cross_cai_cost_ok = True
                        # Try changing the T-starting codon
                        if not _cross_cai_cost_ok and next_codon_idx < len(aas):
                            _t_aa = aas[next_codon_idx]
                            _t_current = sequence[next_codon_idx*3:next_codon_idx*3+3]
                            _t_current_w = usage.get(_t_current, 0.0)
                            _t_non_t_start = [c for c in sorted_codons[_t_aa] if c[0] != "T"]
                            if _t_non_t_start:
                                _best_non_t_start_w = usage.get(_t_non_t_start[0], 0.0)
                                if _t_current_w - _best_non_t_start_w < EUKARYOTE_CAI_GT_COST_THRESHOLD:
                                    _cross_cai_cost_ok = True
                        if not _cross_cai_cost_ok:
                            # CAI cost too high for both possible fixes — accept
                            # the cross-codon GT for eukaryotes
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

                    # Strategy 3 (Issue 2): Deep backtracking — 3-codon coordinated swap
                    # When 2-codon swap can't eliminate GT, try modifying the GT codon
                    # plus TWO neighboring codons simultaneously. This is especially
                    # effective for Valine GTs where all codons contain GT but the
                    # splice score can be reduced by changing the surrounding context.
                    if not fixed_any:
                        # Try the adjacent codon pair (codon_idx-1, codon_idx, codon_idx+1)
                        for n1_offset in [-1, 1]:
                            n1_idx = codon_idx + n1_offset
                            if not (0 <= n1_idx < len(aas)):
                                continue
                            n1_aa = aas[n1_idx]
                            n1_current = sequence[n1_idx*3:n1_idx*3+3]
                            # Second neighbor — try both sides for maximum context coverage
                            for n2_offset in [-2, 2, -3, 3]:
                                n2_idx = codon_idx + n2_offset
                                if not (0 <= n2_idx < len(aas)) or n2_idx == n1_idx:
                                    continue
                                n2_aa = aas[n2_idx]
                                n2_current = sequence[n2_idx*3:n2_idx*3+3]
                                # For Valine GTs: try ALL V codons (all contain GT but
                                # give different splice scores due to 3rd base context).
                                # For neighbors: try all alternatives to maximize the
                                # chance of finding a low-scoring 9-mer context.
                                v_limit = len(sorted_codons[aa])  # Try ALL V alternatives
                                n1_limit = len(sorted_codons[n1_aa])  # Try ALL neighbor alts
                                n2_limit = min(len(sorted_codons[n2_aa]), 6)  # Cap for performance
                                for v_alt in sorted_codons[aa][:v_limit]:
                                    for n1_alt in sorted_codons[n1_aa][:n1_limit]:
                                        if n1_alt == n1_current and v_alt == current:
                                            continue
                                        for n2_alt in sorted_codons[n2_aa][:n2_limit]:
                                            if n2_alt == n2_current and n1_alt == n1_current and v_alt == current:
                                                continue
                                            test = list(sequence)
                                            test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                            test[n1_idx*3:n1_idx*3+3] = list(n1_alt)
                                            test[n2_idx*3:n2_idx*3+3] = list(n2_alt)
                                            test_str = "".join(test)
                                            if gt_pos < len(test_str) - 1 and test_str[gt_pos:gt_pos+2] == "GT":
                                                new_s = score_donor(test_str, gt_pos)
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

                    # Strategy 4 (Issue 2): Frame-shift approach — when GT is at a
                    # codon boundary and can't be eliminated within the current codon,
                    # adjust the PREVIOUS codon to shift the reading frame. This is
                    # especially effective for HBB where GT/AG dinucleotides persist
                    # because the preceding codon ends with G and the next starts with T.
                    if not fixed_any and codon_idx > 0:
                        prev_idx = codon_idx - 1
                        prev_aa = aas[prev_idx]
                        prev_current = sequence[prev_idx*3:prev_idx*3+3]
                        for prev_alt in sorted_codons[prev_aa]:
                            if prev_alt == prev_current:
                                continue
                            # Check if this previous codon swap changes the boundary
                            # Check GT at boundary between prev_alt and current codon
                            if prev_alt[-1] == "G" and sequence[codon_idx*3] == "T":
                                # This would create a new cross-codon GT — skip
                                continue
                            test = sequence[:prev_idx*3] + prev_alt + sequence[prev_idx*3+3:]
                            # Now check if the original GT at gt_pos is still there
                            if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                                new_s = score_donor(test, gt_pos)
                            else:
                                new_s = ELIMINATED_SITE_SCORE
                            if new_s < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                            # Also try combining previous codon swap with GT codon swap
                            for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                                if v_alt == current:
                                    continue
                                test2 = list(test)  # test already has prev_alt applied
                                test2[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                test2_str = "".join(test2)
                                if gt_pos < len(test2_str) - 1 and test2_str[gt_pos:gt_pos+2] == "GT":
                                    new_s = score_donor(test2_str, gt_pos)
                                else:
                                    new_s = ELIMINATED_SITE_SCORE
                                if new_s < cryptic_splice_threshold:
                                    sequence = test2_str
                                    fixed_any = True
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

                    # Strategy 3 (Issue 2): Deep backtracking for AG — 3-codon coordinated swap
                    if not fixed_any:
                        for n1_offset in [-1, 1]:
                            n1_idx = codon_idx + n1_offset
                            if not (0 <= n1_idx < len(aas)):
                                continue
                            n1_aa = aas[n1_idx]
                            n1_current = sequence[n1_idx*3:n1_idx*3+3]
                            for n2_offset in [-2, 2, -3, 3]:
                                n2_idx = codon_idx + n2_offset
                                if not (0 <= n2_idx < len(aas)) or n2_idx == n1_idx:
                                    continue
                                n2_aa = aas[n2_idx]
                                n2_current = sequence[n2_idx*3:n2_idx*3+3]
                                # Same expanded search as donor Strategy 3
                                v_limit = len(sorted_codons[aa])
                                n1_limit = len(sorted_codons[n1_aa])
                                n2_limit = min(len(sorted_codons[n2_aa]), 6)
                                for v_alt in sorted_codons[aa][:v_limit]:
                                    for n1_alt in sorted_codons[n1_aa][:n1_limit]:
                                        if n1_alt == n1_current and v_alt == current:
                                            continue
                                        for n2_alt in sorted_codons[n2_aa][:n2_limit]:
                                            if n2_alt == n2_current and n1_alt == n1_current and v_alt == current:
                                                continue
                                            test = list(sequence)
                                            test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                            test[n1_idx*3:n1_idx*3+3] = list(n1_alt)
                                            test[n2_idx*3:n2_idx*3+3] = list(n2_alt)
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

                    # Strategy 4 (Issue 2): Frame-shift approach for AG
                    if not fixed_any and codon_idx > 0:
                        prev_idx = codon_idx - 1
                        prev_aa = aas[prev_idx]
                        prev_current = sequence[prev_idx*3:prev_idx*3+3]
                        for prev_alt in sorted_codons[prev_aa]:
                            if prev_alt == prev_current:
                                continue
                            if prev_alt[-1] == "A" and sequence[codon_idx*3] == "G":
                                continue  # Would create new cross-codon AG
                            test = sequence[:prev_idx*3] + prev_alt + sequence[prev_idx*3+3:]
                            if ag_pos < len(test) - 1 and test[ag_pos:ag_pos+2] == "AG":
                                new_s = score_acceptor(test, ag_pos)
                            else:
                                new_s = ELIMINATED_SITE_SCORE
                            if new_s < cryptic_splice_threshold:
                                sequence = test
                                fixed_any = True
                                break
                            for v_alt in sorted_codons[aa][:TOP_CAI_ALTERNATIVES]:
                                if v_alt == current:
                                    continue
                                test2 = list(test)
                                test2[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                test2_str = "".join(test2)
                                if ag_pos < len(test2_str) - 1 and test2_str[ag_pos:ag_pos+2] == "AG":
                                    new_s = score_acceptor(test2_str, ag_pos)
                                else:
                                    new_s = ELIMINATED_SITE_SCORE
                                if new_s < cryptic_splice_threshold:
                                    sequence = test2_str
                                    fixed_any = True
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
            _remaining_d = max_donor_score(sequence)
            _remaining_a = max_acceptor_score(sequence)
            logger.warning(
                "Cryptic splice elimination did not converge after %d iterations. "
                "max_donor=%.2f, max_acceptor=%.2f (threshold=%.2f)",
                MAX_SPLICE_ELIMINATION_ITERATIONS, _remaining_d, _remaining_a,
                cryptic_splice_threshold,
            )
            warnings.append(
                f"Cryptic splice elimination capped at {MAX_SPLICE_ELIMINATION_ITERATIONS} iterations. "
                f"max_donor={_remaining_d:.2f}, max_acceptor={_remaining_a:.2f} "
                f"(threshold={cryptic_splice_threshold})."
            )

    # Step: Disrupt CpG dinucleotides to avoid CpG islands
    # EUKARYOTE-ONLY: CpG islands are a eukaryotic gene regulation concern.
    # Prokaryotes don't methylate CpG dinucleotides, so avoidance is unnecessary.
    #
    # Key improvement: Don't break on the first unfixed CpG position.
    # Instead, continue trying all positions. Accept swaps that eliminate
    # a specific CG at pos even if a new CG is created elsewhere (the new
    # one might be fixable in a subsequent iteration or might not contribute
    # to a CpG island). Only require that the specific CG at pos is eliminated,
    # not that the global CpG count decreases.
    if not is_prokaryote:
        for _cpg_iteration in range(MAX_CPG_DISRUPTION_ITERATIONS):
            cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if not cpg_positions:
                break
            any_fixed_this_iter = False
            for pos in cpg_positions:
                # Re-check: a previous fix may have eliminated this CG
                if sequence[pos:pos+2] != "CG":
                    continue
                left_ci = pos // 3          # codon containing the C
                right_ci = (pos + 1) // 3   # codon containing the G
                is_cross_codon = (left_ci != right_ci)

                # Strategy 1: Single-codon swap — try both the C-codon and the G-codon
                fixed = False
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
                        # Accept if the specific CG at pos is eliminated
                        # (relaxed: don't require global CpG decrease)
                        if test[pos:pos+2] != "CG":
                            sequence = test
                            fixed = True
                            any_fixed_this_iter = True
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
                            # Accept if the specific CG at pos is eliminated
                            # (relaxed: don't require global CpG decrease)
                            if test_str[pos:pos+2] != "CG":
                                sequence = test_str
                                fixed = True
                                any_fixed_this_iter = True
                                break
                        if fixed:
                            break

            if not any_fixed_this_iter:
                break
        else:
            _remaining_cpg = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if _remaining_cpg:
                logger.warning(
                    "CpG disruption did not converge after %d iterations. "
                    "%d CpG dinucleotides remain.",
                    MAX_CPG_DISRUPTION_ITERATIONS, len(_remaining_cpg),
                )
                warnings.append(
                    f"CpG disruption capped at {MAX_CPG_DISRUPTION_ITERATIONS} iterations. "
                    f"{len(_remaining_cpg)} CpG dinucleotides could not be eliminated."
                )

    # Step: Reconciliation after cryptic splice elimination
    # EUKARYOTE-ONLY: Only needed if splice/CpG steps ran
    if not is_prokaryote:
        # Check if cryptic splice fixes reintroduced restriction sites
        for site_upper in concrete_sites:
            site_rc = reverse_complement(site_upper)
            if site_upper in sequence or site_rc in sequence:
                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, site_upper, site_rc, usage=usage
                )
                if fixed:
                    sequence = new_seq

    # Step: CpG reconciliation after restriction site reconciliation
    # EUKARYOTE-ONLY: Only needed if CpG avoidance is active
    # Same improvement as the CpG disruption step: don't break on first
    # unfixed position, and accept swaps that eliminate a specific CG
    # even if new CGs are created elsewhere.
    if not is_prokaryote:
        for _cpg_iter in range(MAX_CPG_DISRUPTION_ITERATIONS):
            cpg_positions = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if not cpg_positions:
                break
            any_fixed_this_iter = False
            for pos in cpg_positions:
                # Re-check: a previous fix may have eliminated this CG
                if sequence[pos:pos+2] != "CG":
                    continue
                left_ci = pos // 3
                right_ci = (pos + 1) // 3
                is_cross_codon = (left_ci != right_ci)

                # Strategy 1: Single-codon swap — try both the C-codon and the G-codon
                fixed = False
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
                        # Accept if the specific CG at pos is eliminated
                        # (relaxed: don't require global CpG decrease)
                        if test[pos:pos+2] != "CG":
                            sequence = test
                            fixed = True
                            any_fixed_this_iter = True
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
                            # Accept if the specific CG at pos is eliminated
                            # (relaxed: don't require global CpG decrease)
                            if test_str[pos:pos+2] != "CG":
                                sequence = test_str
                                fixed = True
                                any_fixed_this_iter = True
                                break
                        if fixed:
                            break

            if not any_fixed_this_iter:
                break
        else:
            _remaining_cpg2 = [i for i in range(len(sequence) - 1) if sequence[i:i+2] == "CG"]
            if _remaining_cpg2:
                logger.warning(
                    "CpG reconciliation did not converge after %d iterations. "
                    "%d CpG dinucleotides remain.",
                    MAX_CPG_DISRUPTION_ITERATIONS, len(_remaining_cpg2),
                )
                warnings.append(
                    f"CpG reconciliation capped at {MAX_CPG_DISRUPTION_ITERATIONS} iterations. "
                    f"{len(_remaining_cpg2)} CpG dinucleotides could not be eliminated."
                )

    # Step: CAI Hill Climb — upgrade codons to higher-CAI alternatives
    # without violating any constraint. This recovers CAI lost during
    # constraint resolution by swapping to higher-CAI synonymous codons
    # that don't reintroduce any forbidden pattern.
    #
    # Speed optimizations (Issue 3):
    # - Incremental GC tracking: avoid O(N) gc_content() calls per position
    # - Localized RS check: only check restriction sites near the changed codon
    # - Batch CAI hill-climb: collect all improvements, apply non-conflicting ones
    _MAX_HILL_CLIMB_ITERATIONS = 10

    # Pre-compute site lengths and RCs for localized checking
    _local_rs_check_radius = max((len(s) for s in concrete_sites), default=0) + 3

    # Incremental GC tracking for the hill climb
    _hc_gc_count = sum(1 for b in sequence if b in "GC")
    _hc_n_bases = len(sequence)

    # Pre-build Aho-Corasick scanner for fast local RS check if available
    _hc_scanner = concrete_scanner  # reuse the scanner from Step 2

    for _hc_iter in range(_MAX_HILL_CLIMB_ITERATIONS):
        # Batch mode: collect all possible CAI upgrades, then apply
        # non-conflicting ones together (Issue 3: speed)
        _hc_upgrades: list[tuple[int, str, float, str]] = []  # (ci, alt, alt_cai, current)

        for ci in range(len(aas)):
            aa = aas[ci]
            current = sequence[ci * 3:ci * 3 + 3]
            current_cai = usage.get(current, 0.0)
            for alt in sorted_codons[aa]:
                alt_cai = usage.get(alt, 0.0)
                if alt_cai <= current_cai:
                    break  # sorted_codons is CAI-descending; no improvement possible

                # Quick incremental GC check (O(1) instead of O(N))
                current_gc_bases = sum(1 for b in current if b in "GC")
                alt_gc_bases = sum(1 for b in alt if b in "GC")
                new_gc_count = _hc_gc_count - current_gc_bases + alt_gc_bases
                test_gc = new_gc_count / _hc_n_bases
                if not (gc_lo <= test_gc <= gc_hi):
                    continue

                test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]

                # Localized RS check: only scan the region around the changed codon
                # (Issue 3: avoid O(N) full-sequence scan per position)
                site_ok = True
                if _hc_scanner is not None:
                    # Scan only the region that could contain new sites
                    scan_start = max(0, ci * 3 - _local_rs_check_radius)
                    scan_end = min(len(test), ci * 3 + 3 + _local_rs_check_radius)
                    local_region = test[scan_start:scan_end]
                    local_matches = _hc_scanner.scan(local_region)
                    # Check if any match corresponds to a genuinely new site
                    for _m_pos, _m_site, _ in local_matches:
                        # Map back to absolute position
                        abs_pos = scan_start + _m_pos
                        # Check if this site was already in the old sequence
                        old_local = sequence[scan_start:scan_end]
                        if _m_site not in old_local or abs_pos != sequence.find(_m_site, max(0, abs_pos - len(_m_site))):
                            # More precise: check if the site exists in the old sequence
                            # at any position overlapping the changed region
                            if abs_pos + len(_m_site) > ci * 3 and abs_pos < ci * 3 + 3:
                                # Site overlaps the changed region — check if it's new
                                if _m_site not in sequence or reverse_complement(_m_site) not in sequence:
                                    site_ok = False
                                    break
                                else:
                                    # Site exists in old sequence but may have moved;
                                    # do a full check to be safe
                                    if _m_site in test and reverse_complement(_m_site) in test:
                                        pass  # Both existed before, OK
                                    elif _m_site in test and _m_site not in sequence:
                                        site_ok = False
                                        break
                    # Fallback: full scan for safety on the first iteration
                    if _hc_iter == 0 and site_ok:
                        full_matches = _hc_scanner.scan(test)
                        if full_matches:
                            # Verify no new sites introduced
                            old_matches = _hc_scanner.scan(sequence)
                            new_match_sites = set((p, s) for p, s, _ in full_matches) - set((p, s) for p, s, _ in old_matches)
                            if new_match_sites:
                                site_ok = False
                else:
                    # No scanner: check concrete sites directly (localized)
                    for site_upper in concrete_sites:
                        site_rc = reverse_complement(site_upper)
                        # Only check in the region that could be affected
                        region_start = max(0, ci * 3 - len(site_upper) + 1)
                        region_end = min(len(test), ci * 3 + 3 + len(site_upper) - 1)
                        if site_upper in test[region_start:region_end] or site_rc in test[region_start:region_end]:
                            # Check if this is a genuinely new site
                            if site_upper not in sequence[region_start:region_end] and site_rc not in sequence[region_start:region_end]:
                                site_ok = False
                                break
                if not site_ok:
                    continue

                # No ATTTA motif increase
                if "ATTTA" in test:
                    old_attta = sequence.count("ATTTA")
                    new_attta = test.count("ATTTA")
                    if new_attta > old_attta:
                        continue

                # No 6+ T runs — only check near the changed codon (localized)
                _t_run_region_start = max(0, ci * 3 - T_RUN_LENGTH_THRESHOLD)
                _t_run_region_end = min(len(test), ci * 3 + 3 + T_RUN_LENGTH_THRESHOLD)
                _region_has_long_t = False
                j = _t_run_region_start
                while j < _t_run_region_end:
                    if test[j] == "T":
                        k = j
                        while k < len(test) and test[k] == "T":
                            k += 1
                        if k - j >= T_RUN_LENGTH_THRESHOLD:
                            _region_has_long_t = True
                            break
                        j = k
                    else:
                        j += 1
                if _region_has_long_t:
                    continue

                # No worsening of cryptic splice scores (eukaryotes only)
                # Prokaryotes have no spliceosome, so splice score checks
                # are irrelevant and unnecessarily block CAI upgrades.
                if not is_prokaryote:
                    if max_donor_score(test) > max_donor_score(sequence) + 0.5:
                        continue
                    if max_acceptor_score(test) > max_acceptor_score(sequence) + 0.5:
                        continue

                # Record this upgrade candidate
                _hc_upgrades.append((ci, alt, alt_cai, current))
                break  # Take the best (highest-CAI) alt for this position

        if not _hc_upgrades:
            break

        # Apply upgrades in batch (non-conflicting = different codon positions)
        _applied_any = False
        _applied_positions: set[int] = set()
        for ci, alt, alt_cai, current in _hc_upgrades:
            if ci in _applied_positions:
                continue
            # Re-validate: the sequence may have changed from prior batch swaps
            if sequence[ci * 3:ci * 3 + 3] != current:
                continue  # Position was already modified

            test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]

            # Quick re-validation of key constraints
            current_gc_bases = sum(1 for b in current if b in "GC")
            alt_gc_bases = sum(1 for b in alt if b in "GC")
            new_gc = _hc_gc_count - current_gc_bases + alt_gc_bases
            test_gc = new_gc / _hc_n_bases
            if not (gc_lo <= test_gc <= gc_hi):
                continue

            # Localized RS re-check
            rs_ok = True
            for site_upper in concrete_sites:
                site_rc = reverse_complement(site_upper)
                region_start = max(0, ci * 3 - len(site_upper) + 1)
                region_end = min(len(test), ci * 3 + 3 + len(site_upper) - 1)
                if site_upper in test[region_start:region_end] or site_rc in test[region_start:region_end]:
                    if site_upper not in sequence[region_start:region_end] and site_rc not in sequence[region_start:region_end]:
                        rs_ok = False
                        break
            if not rs_ok:
                continue

            # Apply the upgrade
            sequence = test
            _hc_gc_count = new_gc  # Incremental GC update
            _applied_any = True
            _applied_positions.add(ci)

        if not _applied_any:
            break

    # Step: CAI Recovery Pass (Issue 1)
    # After the hill climb, some positions may still have suboptimal codons
    # because the hill climb only tries one position at a time. This pass
    # systematically checks EVERY position and ALWAYS picks the highest-CAI
    # synonymous codon that doesn't violate any constraint.
    # For prokaryotes, this should close the 0.997→1.0 gap since there are
    # no splice constraints to block upgrades.
    _CAI_RECOVERY_MAX_ITERS = 3
    _rec_gc_count = _hc_gc_count  # Carry forward incremental GC tracking

    for _rec_iter in range(_CAI_RECOVERY_MAX_ITERS):
        _any_recovery = False
        for ci in range(len(aas)):
            aa = aas[ci]
            if aa == "*" or aa == "M":
                continue  # Skip stop and Met (only one codon)
            current = sequence[ci * 3:ci * 3 + 3]
            current_w = usage.get(current, 0.0)
            best_codon = sorted_codons[aa][0]  # Highest-CAI codon for this AA
            best_w = usage.get(best_codon, 0.0)

            if best_w <= current_w or best_codon == current:
                continue  # Already optimal

            # Try the best codon first, then fall back to next-best
            for alt in sorted_codons[aa]:
                alt_w = usage.get(alt, 0.0)
                if alt_w <= current_w:
                    break  # No improvement possible

                # Incremental GC check (O(1))
                cur_gc = sum(1 for b in current if b in "GC")
                alt_gc = sum(1 for b in alt if b in "GC")
                new_gc = _rec_gc_count - cur_gc + alt_gc
                if not (gc_lo <= new_gc / _hc_n_bases <= gc_hi):
                    continue

                test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]

                # Localized RS check
                rs_ok = True
                for site_upper in concrete_sites:
                    site_rc = reverse_complement(site_upper)
                    region_start = max(0, ci * 3 - len(site_upper) + 1)
                    region_end = min(len(test), ci * 3 + 3 + len(site_upper) - 1)
                    if site_upper in test[region_start:region_end] or site_rc in test[region_start:region_end]:
                        if site_upper not in sequence[region_start:region_end] and site_rc not in sequence[region_start:region_end]:
                            rs_ok = False
                            break
                if not rs_ok:
                    continue

                # No ATTTA increase
                if "ATTTA" in test and test.count("ATTTA") > sequence.count("ATTTA"):
                    continue

                # No 6+ T runs (localized check)
                _t_ok = True
                _t_start = max(0, ci * 3 - T_RUN_LENGTH_THRESHOLD)
                _t_end = min(len(test), ci * 3 + 3 + T_RUN_LENGTH_THRESHOLD)
                j = _t_start
                while j < _t_end:
                    if test[j] == "T":
                        k = j
                        while k < len(test) and test[k] == "T":
                            k += 1
                        if k - j >= T_RUN_LENGTH_THRESHOLD:
                            _t_ok = False
                            break
                        j = k
                    else:
                        j += 1
                if not _t_ok:
                    continue

                # No worsening of splice scores (eukaryotes only)
                if not is_prokaryote:
                    if max_donor_score(test) > max_donor_score(sequence) + 0.5:
                        continue
                    if max_acceptor_score(test) > max_acceptor_score(sequence) + 0.5:
                        continue

                # All checks passed — accept the upgrade
                sequence = test
                _rec_gc_count = new_gc
                _any_recovery = True
                logger.debug(
                    "CAI recovery: upgraded codon %d from %s to %s (w %.4f→%.4f)",
                    ci, current, alt, current_w, alt_w,
                )
                break  # Move to next position

            # Paired CAI recovery: if single-codon upgrade was blocked (likely
            # by a restriction site), try a paired swap — upgrade the target codon
            # AND adjust an adjacent codon to eliminate the new restriction site.
            # Net CAI must still improve.
            if not _any_recovery:
                import math as _rec_math
                for _adj_offset in (1, -1):
                    _adj_ci = ci + _adj_offset
                    if _adj_ci < 0 or _adj_ci >= len(aas):
                        continue
                    _adj_aa = aas[_adj_ci]
                    if _adj_aa == "*" or _adj_aa == "M":
                        continue
                    _adj_current = sequence[_adj_ci * 3:_adj_ci * 3 + 3]
                    _adj_current_w = usage.get(_adj_current, 0.0)

                    for alt in sorted_codons[aa]:
                        alt_w = usage.get(alt, 0.0)
                        if alt_w <= current_w:
                            break

                        _adj_sorted = sorted(
                            AA_TO_CODONS.get(_adj_aa, []),
                            key=lambda c: usage.get(c, 0.0),
                            reverse=True,
                        )
                        for _adj_alt in _adj_sorted:
                            if _adj_alt == _adj_current:
                                continue
                            _adj_alt_w = usage.get(_adj_alt, 0.0)
                            # Net log-CAI change must be positive
                            _old_log = _rec_math.log(max(current_w, 1e-10)) + _rec_math.log(max(_adj_current_w, 1e-10))
                            _new_log = _rec_math.log(max(alt_w, 1e-10)) + _rec_math.log(max(_adj_alt_w, 1e-10))
                            if _new_log <= _old_log:
                                continue

                            # Build test sequence with both swaps
                            _lo_ci = min(ci, _adj_ci)
                            _hi_ci = max(ci, _adj_ci)
                            _test_seq = (
                                sequence[:_lo_ci * 3]
                                + (alt if _lo_ci == ci else _adj_alt)
                                + sequence[_lo_ci * 3 + 3:_hi_ci * 3]
                                + (alt if _hi_ci == ci else _adj_alt)
                                + sequence[_hi_ci * 3 + 3:]
                            )

                            # Check: no new restriction sites
                            _rs_ok = True
                            for _site_upper in concrete_sites:
                                _site_rc = reverse_complement(_site_upper)
                                if _site_upper in _test_seq or (_site_rc and _site_rc in _test_seq):
                                    # Check if this is genuinely new
                                    if _site_upper not in sequence or _site_rc not in sequence:
                                        _rs_ok = False
                                        break
                                    else:
                                        # Count occurrences — no net increase
                                        if _test_seq.count(_site_upper) + _test_seq.count(_site_rc) > sequence.count(_site_upper) + sequence.count(_site_rc):
                                            _rs_ok = False
                                            break
                            if not _rs_ok:
                                continue

                            # Check: GC in range
                            _cur_gc_both = sum(1 for b in current if b in "GC") + sum(1 for b in _adj_current if b in "GC")
                            _alt_gc_both = sum(1 for b in alt if b in "GC") + sum(1 for b in _adj_alt if b in "GC")
                            _new_gc = _rec_gc_count - _cur_gc_both + _alt_gc_both
                            if not (gc_lo <= _new_gc / _hc_n_bases <= gc_hi):
                                continue

                            # Check: no new ATTTA
                            if _test_seq.count("ATTTA") > sequence.count("ATTTA"):
                                continue

                            # Check: no 6+ T runs (localized)
                            _t_ok = True
                            _t_start = max(0, _lo_ci * 3 - T_RUN_LENGTH_THRESHOLD)
                            _t_end = min(len(_test_seq), _hi_ci * 3 + 3 + T_RUN_LENGTH_THRESHOLD)
                            _j = _t_start
                            while _j < _t_end:
                                if _test_seq[_j] == "T":
                                    _k = _j
                                    while _k < len(_test_seq) and _test_seq[_k] == "T":
                                        _k += 1
                                    if _k - _j >= T_RUN_LENGTH_THRESHOLD:
                                        _t_ok = False
                                        break
                                    _j = _k
                                else:
                                    _j += 1
                            if not _t_ok:
                                continue

                            # No worsening of splice scores (eukaryotes only)
                            if not is_prokaryote:
                                if max_donor_score(_test_seq) > max_donor_score(sequence) + 0.5:
                                    continue
                                if max_acceptor_score(_test_seq) > max_acceptor_score(sequence) + 0.5:
                                    continue

                            # Accept the paired upgrade
                            sequence = _test_seq
                            _rec_gc_count = _new_gc
                            _any_recovery = True
                            logger.debug(
                                "CAI recovery: paired upgrade codon %d (%s→%s) + "
                                "codon %d (%s→%s)",
                                ci, current, alt,
                                _adj_ci, _adj_current, _adj_alt,
                            )
                            break
                        if _any_recovery:
                            break
                    if _any_recovery:
                        break

        if not _any_recovery:
            break

    # Post-condition: verify sequence still encodes the same protein
    from ..translation import translate
    translated = translate(sequence)
    if translated != protein:
        raise ValueError(
            f"Post-condition violation: optimizer changed the protein. "
            f"Expected '{protein[:20]}...', got '{translated[:20]}...'"
        )
    if len(sequence) != len(aas) * 3:
        raise ValueError(
            f"Post-condition violation: sequence length {len(sequence)} "
            f"!= expected {len(aas) * 3}"
        )

    for w in warnings:
        logger.warning(w)

    return sequence, warnings


# ==============================================================================
# IUPAC Expansion
# ==============================================================================

def _expand_iupac_site(pattern: str) -> list[str]:
    """Expand an IUPAC restriction site pattern into all concrete ACGT sequences.

    E.g., GGCCNNNNNGGCC expands into 4^5 = 1024 concrete sequences.
    For very large expansions, we cap at 4096 to avoid combinatorial explosion.

    Pre-conditions:
    - pattern is a non-empty string containing IUPAC codes

    Post-conditions:
    - all returned strings contain only ACGT characters
    - len(result[0]) == len(pattern) for all results
    """
    if len(pattern) == 0:
        raise ValueError("Pattern must not be empty")

    if not any(b not in "ACGT" for b in pattern):
        return [pattern]

    total_combos = 1
    for b in pattern:
        if b not in "ACGT":
            total_combos *= len(IUPAC_EXPAND.get(b, "A"))

    if total_combos > IUPAC_EXPANSION_CAP:
        logger.warning(
            "IUPAC site %s expands to %d variants (>%d), skipping",
            pattern, total_combos, IUPAC_EXPANSION_CAP,
        )
        return []

    results = [""]
    for b in pattern:
        bases = IUPAC_EXPAND.get(b, b)
        results = [r + x for r in results for x in bases]
    return results


# ==============================================================================
# Predicate Checking (Delegates to Type System — SOC)
# ==============================================================================

def _check_predicates_via_type_system(
    sequence: str,
    gc_lo: float,
    gc_hi: float,
    restriction_sites: list[str],
    cai_threshold: float,
    organism: str,
    cryptic_splice_threshold: float = 3.0,
) -> tuple[list[str], list[str]]:
    """Check all type predicates against the optimized sequence.

    DELEGATES to the type system's evaluate_all_predicates rather than
    re-implementing predicate logic here. This is the single source of truth.

    Pre-conditions:
    - sequence is a valid DNA string
    - organism is in SUPPORTED_ORGANISMS
    - 0.0 <= gc_lo < gc_hi <= 1.0
    - cai_threshold > 0

    Post-conditions:
    - satisfied + failed covers all checked predicates
    - satisfied and failed are disjoint
    """
    from ..type_system import evaluate_all_predicates
    from ..types import Verdict

    if not (0.0 <= gc_lo < gc_hi <= 1.0):
        raise ValueError(f"Invalid GC bounds: gc_lo={gc_lo}, gc_hi={gc_hi}")
    if cai_threshold <= 0:
        raise ValueError(f"CAI threshold must be positive, got {cai_threshold}")

    # Build exon boundaries for a coding sequence (single exon)
    exon_boundaries = [(0, len(sequence))]

    # Get enzyme names from sequences
    enzyme_names = []
    for site in restriction_sites:
        found = False
        for name, seq in RESTRICTION_ENZYMES.items():
            if seq.upper() == site.upper():
                enzyme_names.append(name)
                found = True
                break
        if not found:
            enzyme_names.append(site)  # Use raw sequence as name

    results = evaluate_all_predicates(
        seq=sequence,
        known_exon_boundaries=exon_boundaries,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        enzymes=enzyme_names,
        cryptic_splice_threshold=cryptic_splice_threshold,
    )

    satisfied = []
    failed = []
    for r in results:
        predicate_name = r.predicate
        if r.verdict in (Verdict.PASS, Verdict.LIKELY_PASS):
            satisfied.append(predicate_name)
        elif r.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
            failed.append(predicate_name)
        else:
            # UNCERTAIN: treat as failed for optimization purposes
            failed.append(predicate_name)

    # Verify disjoint
    overlap = set(satisfied) & set(failed)
    if overlap:
        raise ValueError(
            f"Predicates cannot be both satisfied and failed: {overlap}"
        )

    return satisfied, failed


# ==============================================================================
# Main Optimization Entry Point
# ==============================================================================

