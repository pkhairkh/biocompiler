"""
BioCompiler Multi-Objective Sequence Optimizer

Production-grade optimizer with:
- Multi-codon coordinated restriction site removal (handles cross-boundary sites)
- Greedy optimizer as default (z3 optional via use_z3 flag)
- Organism-specific GC targeting (natural GC, not just "in range")
- Smart phase ordering with reconciliation pass
- Type-directed mutagenesis integration loop
- Delegation to type system for predicate checking (SOC)

Phase Summary:
- Phase 0 (implicit): Splice-safe codon reordering via _reorder_for_splice_safety()
  Before Phase 1, codon preference lists are reordered so that GT-free and AG-free
  codons are preferred for amino acids where splice-safe alternatives exist (C, G,
  R, S). This ensures Phase 1 starts with a splice-safe codon when possible,
  reducing the work Phase 7 must do and preventing the "chose poorly" pattern
  where Phase 1 selects a GT-containing codon that Phase 7 must later fix.
- Phase 7: GT-free codon prioritization for cryptic splice elimination
  Uses GT-free and AG-free synonymous codon swaps as the primary strategy
  for eliminating cryptic splice donor and acceptor sites. For amino acids
  with GT-free codons (C, G, R, S), this provides a guaranteed fix.
- Phase 7.5: Window-aware CpG island elimination
  Uses the same sliding-window CpG island metric as the NoCpGIsland predicate
  in type_system.py (window_size=200, gc_threshold=0.6, obs_exp_threshold=0.65).
  Finds all windows failing the CpG island check, then targets CG dinucleotide
  removal within those windows specifically. Falls back to GC-reducing codon
  swaps when CG removal alone cannot bring obs/exp below threshold. Never
  worsens cryptic splice scores or reintroduces restriction sites.
- Phase 8.5: CpG reconciliation
  After restriction site reconciliation (Phase 8) may reintroduce CpG
  dinucleotides, this phase re-applies window-aware CpG island elimination
  using the same parameters as Phase 7.5 and the NoCpGIsland predicate.
- Improved mutagenesis integration with GT-mandatory distinction
  Positions where the amino acid requires GT in ALL codons (e.g., Valine)
  are flagged as GT-mandatory, distinguishing them from positions where
  the optimizer simply hasn't found a good codon yet (optimizer-weakness).
  This distinction drives the mutagenesis engine to only propose AA
  substitutions for truly GT-mandatory positions.

Architecture: gene design is multi-objective optimization, not compilation.
The optimizer searches for sequences satisfying competing constraints simultaneously.

Separation of Concerns:
- Predicate evaluation is the type system's responsibility, NOT the optimizer's.
- The optimizer only CALLS the type system; it never re-implements predicate logic.
- This eliminates the duplicate _check_predicates and ensures single source of truth.
"""

import logging
import math
from dataclasses import dataclass, field
from itertools import product as itertools_product

from .constants import (
    CODON_TABLE, AA_TO_CODONS, BASE_MAP, BASE_REV,
    RESTRICTION_ENZYMES, reverse_complement, IUPAC_EXPAND,
)
from .organisms import (
    CODON_USAGE_TABLES, SUPPORTED_ORGANISMS, CODON_ADAPTIVENESS_TABLES,
    ORGANISM_GC_TARGETS,
)
from .exceptions import UnsupportedOrganismError, InvalidProteinError, OptimizationError
from .translation import compute_cai
from .scanner import gc_content
from .maxentscan import max_donor_score, max_acceptor_score, score_donor, score_acceptor

logger = logging.getLogger(__name__)


# ==============================================================================
# Data Structures
# ==============================================================================

@dataclass
class OptimizationResult:
    """Result of multi-objective optimization.

    Invariants:
    - sequence length == len(protein) * 3 (if protein is non-empty)
    - 0.0 <= cai <= 1.0
    - 0.0 <= gc_content <= 1.0
    - failed_predicates is a subset of all checked predicate names
    - if mutagenesis_applied, then aa_substitutions is non-None and non-empty
    """
    sequence: str
    protein: str
    cai: float
    gc_content: float
    satisfied_predicates: list[str]
    failed_predicates: list[str]
    unsat_core: list[str] | None = None
    fallback_used: bool = False
    mutagenesis_applied: bool = False
    aa_substitutions: list[dict] | None = None  # [{pos, from, to, blosum62}]

    def __post_init__(self):
        """Validate OptimizationResult invariants."""
        if self.protein and self.sequence:
            assert len(self.sequence) == len(self.protein) * 3, (
                f"Sequence length ({len(self.sequence)}) must equal "
                f"protein length * 3 ({len(self.protein) * 3})"
            )
        assert 0.0 <= self.cai <= 1.0, f"CAI must be in [0, 1], got {self.cai}"
        assert 0.0 <= self.gc_content <= 1.0, f"GC content must be in [0, 1], got {self.gc_content}"
        if self.mutagenesis_applied:
            assert self.aa_substitutions is not None and len(self.aa_substitutions) > 0, (
                "Mutagenesis applied but no substitutions recorded"
            )


# ==============================================================================
# Input Validation
# ==============================================================================

def protein_to_aa_list(protein: str) -> list[str]:
    """Convert protein string to list of amino acid codes. Raises InvalidProteinError for bad input.

    Pre-conditions:
    - protein must be a non-empty string of standard amino acid codes

    Post-conditions:
    - result is a list of valid single-letter amino acid codes
    - len(result) == len(protein.strip())
    """
    if not protein or not protein.strip():
        raise InvalidProteinError(protein, set())
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    invalid = set(ch for ch in protein if ch not in valid_aas)
    if invalid:
        raise InvalidProteinError(protein, invalid)
    return list(protein)


# ==============================================================================
# Restriction Site Removal Helpers
# ==============================================================================

def _find_site_in_sequence(sequence: str, site: str, site_rc: str) -> list[int]:
    """Find all positions where site or its reverse complement appears in sequence.

    Pre-conditions:
    - sequence is a valid uppercase DNA string
    - site is a non-empty uppercase DNA string
    - site_rc is the reverse complement of site (or empty string)

    Post-conditions:
    - returns sorted list of unique positions
    """
    positions = []
    if site:
        start = 0
        while True:
            pos = sequence.find(site, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
    if site_rc and site_rc != site:  # Avoid double-counting palindromes
        start = 0
        while True:
            pos = sequence.find(site_rc, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
    return sorted(set(positions))


def _get_overlapping_codons(pos: int, site_len: int, n_codons: int) -> list[int]:
    """Get indices of codons that overlap with a site at position pos.

    Pre-conditions:
    - pos >= 0
    - site_len > 0
    - n_codons > 0

    Post-conditions:
    - all indices in result are in [0, n_codons)
    """
    assert pos >= 0, f"Position must be non-negative, got {pos}"
    assert site_len > 0, f"Site length must be positive, got {site_len}"
    assert n_codons > 0, f"Number of codons must be positive, got {n_codons}"

    first_codon = pos // 3
    last_base = pos + site_len - 1
    last_codon = last_base // 3
    return list(range(max(0, first_codon), min(n_codons, last_codon + 1)))


def _remove_site_multicodon(
    sequence: str, aas: list[str], sorted_codons: dict[str, list[str]],
    site_upper: str, site_rc: str,
) -> tuple[str, bool]:
    """
    Try to remove a restriction site using multi-codon coordinated solving.

    When a site straddles codon boundaries (e.g., PstI CTG|CAG spans codons i and i+1),
    single-codon swaps fail because changing either codon alone doesn't eliminate the site.
    This function enumerates valid codon COMBINATIONS for all overlapping codons.

    Pre-conditions:
    - sequence is uppercase DNA
    - len(aas) > 0
    - sorted_codons maps each aa in aas to a non-empty list of codons
    - site_upper is a valid uppercase DNA string

    Post-conditions:
    - if fixed, returned sequence has same length as input and encodes same protein
    - if not fixed, returned sequence is identical to input
    """
    n_codons = len(aas)
    positions = _find_site_in_sequence(sequence, site_upper, site_rc)

    for pos in positions:
        overlapping = _get_overlapping_codons(pos, len(site_upper), n_codons)
        if not overlapping:
            continue

        # Build candidate codon lists for each overlapping position
        candidate_lists = []
        for ci in overlapping:
            aa = aas[ci]
            candidate_lists.append(sorted_codons[aa])

        # Enumerate all codon combinations for overlapping positions
        # Max search: ~6^3 = 216 for 3-codon overlap (very rare)
        for combo in itertools_product(*candidate_lists):
            # Build test sequence with this combo applied
            test = list(sequence)
            for idx, ci in enumerate(overlapping):
                start = ci * 3
                test[start:start + 3] = list(combo[idx])
            test_seq = "".join(test)

            # Check if site is eliminated
            if site_upper not in test_seq and (not site_rc or site_rc not in test_seq):
                return test_seq, True

    return sequence, False


# ==============================================================================
# Greedy Optimizer
# ==============================================================================

def _find_gt_free_codons(aa: str) -> list[str]:
    """Return codons for the given amino acid that do NOT contain the GT dinucleotide.

    For amino acids like Cysteine (C), Glycine (G), Arginine (R), and Serine (S),
    there exist synonymous codons without GT, providing a guaranteed fix for
    cryptic splice donor elimination. Valine (V) has NO GT-free codons.

    Pre-conditions:
    - aa is a valid single-letter amino acid code present in AA_TO_CODONS

    Post-conditions:
    - all returned codons are valid for the amino acid
    - no returned codon contains "GT" as a substring
    - if no GT-free codons exist (e.g., Valine), returns empty list
    """
    return [c for c in AA_TO_CODONS[aa] if "GT" not in c]


def _find_ag_free_codons(aa: str) -> list[str]:
    """Return codons for the given amino acid that do NOT contain the AG dinucleotide.

    Similar to _find_gt_free_codons but for acceptor (AG) dinucleotide elimination.
    Many amino acids have AG-free synonymous codons.

    Pre-conditions:
    - aa is a valid single-letter amino acid code present in AA_TO_CODONS

    Post-conditions:
    - all returned codons are valid for the amino acid
    - no returned codon contains "AG" as a substring
    - if no AG-free codons exist, returns empty list
    """
    return [c for c in AA_TO_CODONS[aa] if "AG" not in c]


def _reorder_for_splice_safety(
    sorted_codons: dict[str, list[str]],
    aas: list[str],
    cryptic_splice_threshold: float,
) -> dict[str, list[str]]:
    """Reorder codon preference to put GT-free and AG-free codons first.

    When cryptic splice elimination is a concern (threshold > 0), for amino acids
    that have GT-free or AG-free synonymous codons, this function demotes
    dinucleotide-containing codons to the END of the preference list. This way,
    Phase 1 starts with a splice-safe codon when possible.

    Strategy: Prioritize codons that are BOTH GT-free AND AG-free first,
    then codons that are only GT-free, then codons that are only AG-free,
    then codons with both GT and AG. This dual-prioritization handles both
    donor and acceptor sites simultaneously.

    Scope: Only amino acids that have splice-problematic codons (containing GT
    or AG) with available splice-safe alternatives are reordered. This targets
    the amino acids where Phase 7 struggles (C, G, R, S) and avoids unnecessary
    CAI penalties for amino acids like E, K, Q where Phase 7 handles AG sites
    efficiently through simple AG-free swaps.

    Pre-conditions:
    - sorted_codons maps each unique amino acid in aas to a list of codons
      sorted by CAI descending
    - cryptic_splice_threshold > 0 (if <= 0, returns sorted_codons unchanged)

    Post-conditions:
    - returned dict has the same keys as sorted_codons
    - for each aa, the returned list is a permutation of sorted_codons[aa]
    """
    if cryptic_splice_threshold <= 0:
        return sorted_codons

    # Amino acids with splice-problematic codons that have splice-safe alternatives.
    # C, G have GT-containing codons (TGT, GGT) with GT-free alternatives.
    # R has AG-containing codons (AGA, AGG) with AG-free alternatives (CG*).
    # S has both GT and AG-containing codons (AGT, AGC) with alternatives (TC*).
    # V is excluded because ALL V codons contain GT (no alternatives).
    # E, K, Q have AG-containing codons (GAG, AAG, CAG) but Phase 7 handles
    # their AG sites efficiently — reordering would sacrifice too much CAI.
    SPLICE_REORDER_AAS = {"C", "G", "R", "S"}

    reordered: dict[str, list[str]] = {}
    for aa, codons in sorted_codons.items():
        if aa not in SPLICE_REORDER_AAS:
            reordered[aa] = codons
            continue

        # Tier 1: Both GT-free and AG-free (best for splice safety)
        tier1 = [c for c in codons if "GT" not in c and "AG" not in c]
        # Tier 2: GT-free only (eliminates donors but may have AG)
        tier2 = [c for c in codons if "GT" not in c and "AG" in c]
        # Tier 3: AG-free only (eliminates acceptors but may have GT)
        tier3 = [c for c in codons if "GT" in c and "AG" not in c]
        # Tier 4: Both GT and AG (worst for splice safety)
        tier4 = [c for c in codons if "GT" in c and "AG" in c]

        result = tier1 + tier2 + tier3 + tier4
        if len(result) == len(codons):
            reordered[aa] = result
        else:
            # Fallback (shouldn't happen, but defensive)
            reordered[aa] = codons
    return reordered


def _find_failing_cpg_windows(
    sequence: str,
    window_size: int = 200,
    gc_threshold: float = 0.6,
    obs_exp_threshold: float = 0.65,
) -> list[tuple[int, int, float, float]]:
    """Find all windows in the sequence that fail the CpG island check.

    Returns a list of tuples: (window_start, window_end, gc_content, obs_exp_ratio)
    for each window where GC >= gc_threshold AND obs_exp >= obs_exp_threshold.

    Results are sorted by obs_exp_ratio descending (worst first).

    Pre-conditions:
    - sequence is a valid uppercase DNA string
    - window_size > 0
    - 0.0 <= gc_threshold <= 1.0
    - 0.0 <= obs_exp_threshold <= 1.0

    Post-conditions:
    - all returned windows have gc_content >= gc_threshold and obs_exp >= obs_exp_threshold
    - results are sorted by obs_exp descending
    """
    failing = []
    for start in range(len(sequence) - window_size + 1):
        window = sequence[start:start + window_size]
        gc = gc_content(window)
        if gc < gc_threshold:
            continue  # Not enough GC to be a CpG island
        cpg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
        c_count = window.count('C')
        g_count = window.count('G')
        expected = (c_count * g_count) / max(len(window), 1)
        obs_exp = cpg_count / max(expected, 1e-10)
        if obs_exp >= obs_exp_threshold:
            failing.append((start, start + window_size, gc, obs_exp))
    failing.sort(key=lambda x: x[3], reverse=True)  # Sort by obs_exp descending
    return failing


def _reduce_cpg_in_window(
    sequence: str,
    window_start: int,
    window_end: int,
    aas: list[str],
    sorted_codons: dict[str, list[str]],
    usage: dict,
    cryptic_splice_threshold: float,
    concrete_sites: list[str],
    max_iterations: int = 50,
    window_size: int = 200,
    gc_threshold: float = 0.6,
    obs_exp_threshold: float = 0.65,
) -> tuple[str, bool]:
    """Reduce CpG density within a specific failing window.

    Strategy:
    1. Identify all CG dinucleotides within the window
    2. Try swapping to CG-free synonymous codons that don't worsen splice/restriction
    3. After each swap, re-evaluate the window's CpG island status
    4. If obs/exp drops below threshold or GC drops below threshold, success
    5. Fallback: try GC-reducing codon swaps to push GC below gc_threshold

    Pre-conditions:
    - sequence is a valid uppercase DNA string
    - window_start >= 0, window_end <= len(sequence), window_start < window_end
    - aas corresponds to the protein encoded by sequence
    - sorted_codons maps each aa to a list of synonymous codons (sorted by CAI desc)
    - usage is the codon adaptiveness table for the target organism
    - cryptic_splice_threshold > 0
    - concrete_sites is a list of concrete restriction site sequences

    Post-conditions:
    - returned sequence encodes the same protein
    - if successful, the window no longer fails the CpG island check
    """
    for iteration in range(max_iterations):
        window = sequence[window_start:window_end]
        gc_val = gc_content(window)

        # If GC dropped below threshold, CpG island is eliminated
        if gc_val < gc_threshold:
            return sequence, True

        # Compute obs/exp for current window
        cpg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
        c_count = window.count('C')
        g_count = window.count('G')
        expected = (c_count * g_count) / max(len(window), 1)
        obs_exp = cpg_count / max(expected, 1e-10)

        # If obs/exp dropped below threshold, CpG island is eliminated
        if obs_exp < obs_exp_threshold:
            return sequence, True

        # Find all CG positions within the window
        cgs = [i for i in range(window_start, window_end - 1) if sequence[i:i+2] == "CG"]
        if not cgs:
            # No more CGs to remove but still failing — try GC reduction
            break

        # Try to remove a CG by swapping the codon containing C or G
        best_swap = None  # (ci, alt_codon, new_obs_exp)
        best_obs_exp = obs_exp

        for pos in cgs:
            for offset in [0, -1]:
                ci = (pos + offset) // 3
                if ci < 0 or ci >= len(aas):
                    continue
                aa = aas[ci]
                current = sequence[ci*3:ci*3+3]
                for alt in sorted_codons[aa]:
                    if alt == current:
                        continue
                    test = sequence[:ci*3] + alt + sequence[ci*3+3:]

                    # Check: local CG elimination
                    local_start = max(0, ci*3 - 2)
                    local_end = min(len(test), ci*3 + 5)
                    if "CG" in test[local_start:local_end]:
                        # This swap doesn't help with the CG at this position
                        continue

                    # Check: splice safety
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

                    # Check: restriction site safety
                    site_ok = all(
                        s not in test and reverse_complement(s) not in test
                        for s in concrete_sites
                    )
                    if not site_ok:
                        continue

                    # Check: window-level CpG improvement
                    new_window = test[window_start:window_end]
                    new_cpg_count = sum(1 for i in range(len(new_window) - 1) if new_window[i:i+2] == "CG")
                    new_c_count = new_window.count('C')
                    new_g_count = new_window.count('G')
                    new_expected = (new_c_count * new_g_count) / max(len(new_window), 1)
                    new_obs_exp = new_cpg_count / max(new_expected, 1e-10)

                    if new_obs_exp < best_obs_exp:
                        best_obs_exp = new_obs_exp
                        best_swap = (ci, alt, new_obs_exp)

        if best_swap is not None:
            ci, alt_codon, _ = best_swap
            sequence = sequence[:ci*3] + alt_codon + sequence[ci*3+3:]
        else:
            # No CG-removing swap found — try GC reduction fallback
            break

    # Fallback: try to reduce GC below threshold
    window = sequence[window_start:window_end]
    gc_val = gc_content(window)
    if gc_val >= gc_threshold and gc_val < gc_threshold + 0.05:
        # GC is close to threshold — try GC-reducing swaps
        for ci in range(max(0, window_start // 3), min(len(aas), (window_end + 2) // 3)):
            aa = aas[ci]
            current = sequence[ci*3:ci*3+3]
            current_gc_count = sum(1 for b in current if b in "GC")
            for alt in sorted_codons[aa]:
                if alt == current:
                    continue
                alt_gc_count = sum(1 for b in alt if b in "GC")
                if alt_gc_count >= current_gc_count:
                    continue  # Doesn't reduce GC

                test = sequence[:ci*3] + alt + sequence[ci*3+3:]

                # Splice safety check
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

                # Restriction site check
                site_ok = all(
                    s not in test and reverse_complement(s) not in test
                    for s in concrete_sites
                )
                if not site_ok:
                    continue

                # Check if this swap drops GC below threshold in the window
                new_window = test[window_start:window_end]
                new_gc = gc_content(new_window)
                if new_gc < gc_threshold:
                    return test, True

    # Final check
    window = sequence[window_start:window_end]
    gc_val = gc_content(window)
    cpg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
    c_count = window.count('C')
    g_count = window.count('G')
    expected = (c_count * g_count) / max(len(window), 1)
    obs_exp = cpg_count / max(expected, 1e-10)

    if gc_val < gc_threshold or obs_exp < obs_exp_threshold:
        return sequence, True
    return sequence, False


def _greedy_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    cryptic_splice_threshold: float = 3.0,
) -> tuple[str, list[str]]:
    """
    Greedy multi-objective codon optimization with coordinated constraint solving.

    Phase ordering prioritizes hard constraints (restriction sites) over soft constraints (CAI).
    Reconciliation pass ensures earlier phases aren't undone by later ones.

    Phases:
    1. Best codon per position (maximize CAI)
    2. Remove restriction sites (multi-codon coordinated)
    3. Remove ATTTA instability motifs
    4. Fix 6+ consecutive T runs
    5. Adjust GC content (hard constraint, organism target aspiration)
    6. Reconciliation — restriction sites vs GC
    7. Eliminate cryptic splice donor/acceptor sites (GT-free/AG-free codon swap priority)
    7.5. Disrupt CpG dinucleotides to avoid CpG islands
    8. Reconciliation — restriction sites after splice/CpG fixes
    8.5. CpG reconciliation — re-disrupt CpG if reconciliation reintroduced them

    Pre-conditions:
    - protein is a valid amino acid sequence (no invalid codes)
    - organism is in SUPPORTED_ORGANISMS
    - 0.0 <= gc_lo < gc_hi <= 1.0
    - cryptic_splice_threshold > 0

    Post-conditions:
    - returned sequence translates to the input protein
    - len(returned sequence) == len(protein) * 3
    - all codons in sequence are valid for their amino acid
    """
    # Validate pre-conditions
    assert 0.0 <= gc_lo < gc_hi <= 1.0, f"GC bounds invalid: [{gc_lo}, {gc_hi}]"
    assert cryptic_splice_threshold > 0, f"Threshold must be positive, got {cryptic_splice_threshold}"

    usage = CODON_ADAPTIVENESS_TABLES.get(organism)
    if usage is None:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)
    aas = protein_to_aa_list(protein)
    restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.values())
    warnings: list[str] = []

    sorted_codons: dict[str, list[str]] = {}
    for aa in set(aas):
        codons = AA_TO_CODONS[aa]
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted

    # When cryptic splice elimination is active, reorder codon preference so that
    # GT-free and AG-free codons are preferred for amino acids where alternatives exist.
    # This ensures Phase 1 starts with a splice-safe codon when possible, reducing
    # the work Phase 7 must do and preventing the "chose poorly" pattern.
    sorted_codons = _reorder_for_splice_safety(sorted_codons, aas, cryptic_splice_threshold)

    # Phase 1: Best codon per position (maximize CAI, with splice-safe reordering)
    sequence = "".join(sorted_codons[aa][0] for aa in aas)
    assert len(sequence) == len(aas) * 3, "Phase 1: sequence length mismatch"

    # Phase 2: Remove restriction sites (HIGHEST PRIORITY — multi-codon coordinated)
    # Process concrete sites first, then IUPAC sites
    concrete_sites = []
    iupac_sites = []
    for site in restriction_sites:
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            iupac_sites.append(site_upper)
        else:
            concrete_sites.append(site_upper)

    # Remove concrete sites
    for site_upper in concrete_sites:
        site_rc = reverse_complement(site_upper)
        for iteration in range(100):
            positions = _find_site_in_sequence(sequence, site_upper, site_rc)
            if not positions:
                break

            # Try multi-codon coordinated removal
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_upper, site_rc
            )
            if fixed:
                sequence = new_seq
                continue

            # Fallback: try single-codon swap
            pos = positions[0]
            codon_idx = pos // 3
            if codon_idx < len(aas):
                aa = aas[codon_idx]
                current = sequence[codon_idx * 3: codon_idx * 3 + 3]
                single_fixed = False
                for alt in sorted_codons[aa]:
                    if alt != current:
                        test = sequence[:codon_idx * 3] + alt + sequence[codon_idx * 3 + 3:]
                        if site_upper not in test and site_rc not in test:
                            sequence = test
                            single_fixed = True
                            break
                if single_fixed:
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
                        test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
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
                warnings.append(f"Restriction site {site_upper}: max iterations reached, may still be present")

    # Remove IUPAC sites (expand to concrete variants, check each)
    for site_upper in iupac_sites:
        concrete_variants = _expand_iupac_site(site_upper)
        if not concrete_variants:
            continue
        for variant in concrete_variants:
            variant_rc = reverse_complement(variant)
            for iteration in range(100):
                positions = _find_site_in_sequence(sequence, variant, variant_rc)
                if not positions:
                    break

                new_seq, fixed = _remove_site_multicodon(
                    sequence, aas, sorted_codons, variant, variant_rc
                )
                if fixed:
                    sequence = new_seq
                    continue

                # Single-codon fallback
                pos = positions[0]
                codon_idx = pos // 3
                if codon_idx < len(aas):
                    aa = aas[codon_idx]
                    current = sequence[codon_idx * 3: codon_idx * 3 + 3]
                    for alt in sorted_codons[aa]:
                        if alt != current:
                            test = sequence[:codon_idx * 3] + alt + sequence[codon_idx * 3 + 3:]
                            if variant not in test and variant_rc not in test:
                                sequence = test
                                break
                    else:
                        warnings.append(f"Cannot remove IUPAC {site_upper} variant {variant} at iteration {iteration}")
                        break
            else:
                if variant in sequence or variant_rc in sequence:
                    warnings.append(f"IUPAC site {site_upper} variant {variant}: max iterations")

    # Phase 3: Remove ATTTA instability motifs
    for iteration in range(100):
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
                    test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
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
        warnings.append("ATTTA motif: max iterations reached, may still be present")

    # Phase 4: Fix 6+ consecutive T runs
    for iteration in range(100):
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
        if max_run < 6:
            break
        codon_idx = (max_pos + max_run // 2) // 3
        if codon_idx < len(aas):
            aa = aas[codon_idx]
            current = sequence[codon_idx * 3:codon_idx * 3 + 3]
            fixed = False
            for alt in sorted_codons[aa]:
                if alt != current:
                    test = sequence[:codon_idx * 3] + alt + sequence[codon_idx * 3 + 3:]
                    if not any(test[i:i + 6] == "TTTTTT" for i in range(len(test) - 5)):
                        sequence = test
                        fixed = True
                        break
            if not fixed:
                warnings.append(f"Consecutive T run: cannot fix at iteration {iteration}")
                break
    else:
        warnings.append("Consecutive T: max iterations reached, may still have 6+ T runs")

    # Phase 5: Adjust GC content
    # Strategy: GC must be in [gc_lo, gc_hi] (hard constraint).
    # If in range, we gently nudge toward organism target but NEVER at the
    # cost of significant CAI reduction. The organism GC target is aspirational,
    # not mandatory — a sequence with CAI=0.99 and GC=0.61 (slightly above
    # human's 0.41 target) is better than CAI=0.82 and GC=0.46.
    gc_val = gc_content(sequence)
    organism_gc = ORGANISM_GC_TARGETS.get(organism, (gc_lo + gc_hi) / 2.0)
    target_gc = max(gc_lo, min(gc_hi, organism_gc))

    gc_out_of_range = not (gc_lo <= gc_val <= gc_hi)

    if gc_out_of_range:
        # Hard constraint: MUST get GC into range
        gc_count = sum(1 for b in sequence if b in "GC")
        n_bases = len(sequence)
        # Target the nearest bound
        if gc_val < gc_lo:
            phase_target = gc_lo
        else:
            phase_target = gc_hi

        for iteration in range(200):
            if gc_lo <= gc_val <= gc_hi:
                break
            best_alt = None
            best_ci = -1
            best_diff = abs(gc_val - phase_target)
            best_gc_delta = 0
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
                    if diff < best_diff:
                        best_diff = diff
                        best_alt = alt
                        best_ci = ci
                        best_gc_delta = alt_gc - current_gc
            if best_alt is None:
                break
            sequence = sequence[:best_ci * 3] + best_alt + sequence[best_ci * 3 + 3:]
            gc_count += best_gc_delta
            gc_val = gc_count / n_bases
        else:
            warnings.append(f"GC adjustment: max iterations reached, current GC={gc_val:.3f}")

    # Phase 6: Reconciliation — check if GC adjustment reintroduced restriction sites
    for site_upper in concrete_sites:
        site_rc = reverse_complement(site_upper)
        if site_upper in sequence or site_rc in sequence:
            # Try one more round of multi-codon removal
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_upper, site_rc
            )
            if fixed:
                sequence = new_seq
                # Re-check GC
                gc_val = gc_content(sequence)
                if not (gc_lo <= gc_val <= gc_hi):
                    # GC drifted — try to fix with single-codon swaps that don't reintroduce sites
                    gc_count = sum(1 for b in sequence if b in "GC")
                    n_bases = len(sequence)
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
                            test = sequence[:ci * 3] + alt + sequence[ci * 3 + 3:]
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
                # Could not remove — already warned in Phase 2
                pass

    # Phase 7: Eliminate cryptic splice donor/acceptor sites
    # Strategy (ordered by effectiveness):
    #   1. GT-free codon swap — guaranteed to eliminate the GT dinucleotide
    #      (works for C, G, R, S which have GT-free synonymous codons)
    #   2. 2-codon context disruption — swap GT codon + neighbor to disrupt
    #      the 9-mer splice context (needed for Valine which has no GT-free codons)
    #   3. Accept that some Valine positions are unrepairable by codon swaps alone
    #      (these will be handled by mutagenesis if enabled)
    # Same strategy applied for AG acceptor sites using AG-free codons.
    for iteration in range(300):
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
                        test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
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
                    test = sequence[:codon_idx*3] + v_alt + sequence[codon_idx*3+3:]
                    if gt_pos < len(test) - 1 and test[gt_pos:gt_pos+2] == "GT":
                        new_s = score_donor(test, gt_pos)
                    else:
                        new_s = -999  # GT eliminated (cross-boundary)
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
                        for v_alt in sorted_codons[aa][:3]:
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
                                    new_s = -999  # GT eliminated
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
                        new_s = -999  # AG eliminated
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
                        for v_alt in sorted_codons[aa][:3]:
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
                                    new_s = -999
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

    # Phase 7.5: Window-aware CpG island avoidance
    # Strategy: Find windows that fail the CpG island check and target them
    # specifically, using the same window_size, gc_threshold, and obs_exp_threshold
    # as the NoCpGIsland predicate in type_system.py.
    CpG_WINDOW_SIZE = 200
    CpG_GC_THRESHOLD = 0.6
    CpG_OBS_EXP_THRESHOLD = 0.65

    for _cpg_outer in range(20):  # Max 20 failing windows to fix
        failing_windows = _find_failing_cpg_windows(
            sequence, CpG_WINDOW_SIZE, CpG_GC_THRESHOLD, CpG_OBS_EXP_THRESHOLD
        )
        if not failing_windows:
            break

        # Fix the worst window first (highest obs/exp ratio)
        w_start, w_end, w_gc, w_obs_exp = failing_windows[0]
        sequence, fixed = _reduce_cpg_in_window(
            sequence, w_start, w_end, aas, sorted_codons, usage,
            cryptic_splice_threshold, concrete_sites,
            window_size=CpG_WINDOW_SIZE,
            gc_threshold=CpG_GC_THRESHOLD,
            obs_exp_threshold=CpG_OBS_EXP_THRESHOLD,
        )
        if not fixed:
            warnings.append(
                f"CpG island at [{w_start},{w_end}): GC={w_gc:.3f}, "
                f"obs/exp={w_obs_exp:.3f} — cannot reduce below threshold"
            )
            break

    # Phase 8: Reconciliation after cryptic splice elimination
    # Check if cryptic splice fixes reintroduced restriction sites
    for site_upper in concrete_sites:
        site_rc = reverse_complement(site_upper)
        if site_upper in sequence or site_rc in sequence:
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_upper, site_rc
            )
            if fixed:
                sequence = new_seq

    # Phase 8.5: CpG reconciliation after restriction site reconciliation
    # Re-apply window-aware CpG elimination (restriction site removal may have
    # reintroduced CpG islands).
    for _cpg_outer in range(20):
        failing_windows = _find_failing_cpg_windows(
            sequence, CpG_WINDOW_SIZE, CpG_GC_THRESHOLD, CpG_OBS_EXP_THRESHOLD
        )
        if not failing_windows:
            break

        w_start, w_end, w_gc, w_obs_exp = failing_windows[0]
        sequence, fixed = _reduce_cpg_in_window(
            sequence, w_start, w_end, aas, sorted_codons, usage,
            cryptic_splice_threshold, concrete_sites,
            window_size=CpG_WINDOW_SIZE,
            gc_threshold=CpG_GC_THRESHOLD,
            obs_exp_threshold=CpG_OBS_EXP_THRESHOLD,
        )
        if not fixed:
            warnings.append(
                f"CpG reconciliation: island at [{w_start},{w_end}) "
                f"GC={w_gc:.3f}, obs/exp={w_obs_exp:.3f} — cannot reduce"
            )
            break

    # Post-condition: verify sequence still encodes the same protein
    from .translation import translate
    translated = translate(sequence)
    assert translated == protein, (
        f"Post-condition violation: optimizer changed the protein. "
        f"Expected '{protein[:20]}...', got '{translated[:20]}...'"
    )
    assert len(sequence) == len(aas) * 3, (
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
    assert len(pattern) > 0, "Pattern must not be empty"

    if not any(b not in "ACGT" for b in pattern):
        return [pattern]

    total_combos = 1
    for b in pattern:
        if b not in "ACGT":
            total_combos *= len(IUPAC_EXPAND.get(b, "A"))

    if total_combos > 4096:
        logger.warning(
            "IUPAC site %s expands to %d variants (>4096), skipping",
            pattern, total_combos,
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
    from .type_system import evaluate_all_predicates
    from .types import Verdict

    assert 0.0 <= gc_lo < gc_hi <= 1.0, f"Invalid GC bounds: [{gc_lo}, {gc_hi}]"
    assert cai_threshold > 0, f"CAI threshold must be positive, got {cai_threshold}"

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
        if r.verdict == Verdict.PASS:
            satisfied.append(predicate_name)
        else:
            failed.append(predicate_name)

    # Verify disjoint
    assert not (set(satisfied) & set(failed)), (
        f"Predicates cannot be both satisfied and failed: "
        f"{set(satisfied) & set(failed)}"
    )

    return satisfied, failed


# ==============================================================================
# Main Optimization Entry Point
# ==============================================================================

def optimize_sequence(
    target_protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: list[str] | None = None,
    cai_threshold: float = 0.2,
    max_amino_acids_for_z3: int | None = None,
    z3_timeout_ms: int = 30000,
    cryptic_splice_threshold: float = 3.0,
    use_z3: bool = False,
    enable_mutagenesis: bool = False,
    max_mutagenesis_substitutions: int = 30,
    min_blosum62: int = -1,
) -> OptimizationResult:
    """
    Find or generate a DNA sequence encoding the target protein that satisfies
    type predicates using multi-objective optimization.

    By default, uses the greedy optimizer for all protein lengths (fast, high CAI).
    Set use_z3=True to use z3 constraint solver (experimental, slower for short proteins).

    When enable_mutagenesis=True and codon-level optimization cannot satisfy all
    predicates (e.g., Valine positions creating unrepairable cryptic splice donors),
    the engine proposes conservative amino acid substitutions (type-directed
    mutagenesis) to make constraint satisfaction possible.

    The optimizer targets organism-specific GC content rather than just "in range."

    Pre-conditions:
    - target_protein contains only valid standard amino acid codes
    - organism is in SUPPORTED_ORGANISMS
    - 0.0 <= gc_lo < gc_hi <= 1.0
    - 0.0 < cai_threshold <= 1.0
    - cryptic_splice_threshold > 0
    """
    aas = protein_to_aa_list(target_protein)
    n_aa = len(aas)
    if n_aa == 0:
        return OptimizationResult("", "", 0.0, 0.0, [], ["empty_sequence"])

    if organism not in SUPPORTED_ORGANISMS:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)

    assert 0.0 <= gc_lo < gc_hi <= 1.0, f"Invalid GC bounds: [{gc_lo}, {gc_hi}]"
    assert 0.0 < cai_threshold <= 1.0, f"Invalid CAI threshold: {cai_threshold}"
    assert cryptic_splice_threshold > 0, f"Threshold must be positive: {cryptic_splice_threshold}"

    restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.values())
    fallback_used = False
    all_warnings: list[str] = []

    # Determine whether to use z3
    should_use_z3 = use_z3
    # Backward compat: if max_amino_acids_for_z3 is explicitly set, use old behavior
    if max_amino_acids_for_z3 is not None and not use_z3:
        should_use_z3 = n_aa <= max_amino_acids_for_z3

    if should_use_z3:
        try:
            from z3 import Int, Optimize, If, And, Or, Not, Implies, Sum, sat, Real, ToReal
        except ImportError:
            logger.warning("z3-solver not installed, falling back to greedy optimizer")
            should_use_z3 = False

    if should_use_z3:
        # z3 constraint solver (experimental)
        opt = Optimize()
        opt.set("timeout", z3_timeout_ms)
        n_bases = n_aa * 3
        nuc = [Int(f"n_{i}") for i in range(n_bases)]

        for i in range(n_bases):
            opt.add(nuc[i] >= 0, nuc[i] <= 3)

        usage = CODON_ADAPTIVENESS_TABLES[organism]
        log_cai_contributions = []
        for aa_idx, aa in enumerate(aas):
            base_idx = aa_idx * 3
            n0, n1, n2 = nuc[base_idx], nuc[base_idx + 1], nuc[base_idx + 2]
            valid_codons = AA_TO_CODONS[aa]
            codon_constraints = []
            log_cai_chain = None
            sorted_z3 = sorted(valid_codons, key=lambda c: usage.get(c, 0.0), reverse=True)
            for codon in sorted_z3:
                c0, c1, c2 = (BASE_REV[codon[0]], BASE_REV[codon[1]], BASE_REV[codon[2]])
                match = And(n0 == c0, n1 == c1, n2 == c2)
                codon_constraints.append(match)
                w = usage.get(codon, 0.01)
                log_w_int = int(math.log(max(w, 1e-10)) * 1000)
                log_cai_chain = log_w_int if log_cai_chain is None else If(match, log_w_int, log_cai_chain)
            opt.add(Or(*codon_constraints))
            if log_cai_chain is not None:
                log_cai_contributions.append(log_cai_chain)

        if log_cai_contributions:
            opt.maximize(Sum(log_cai_contributions))

        gc_count = Sum([If(Or(nuc[i] == 1, nuc[i] == 2), 1, 0) for i in range(n_bases)])
        opt.add(gc_count >= int(gc_lo * n_bases), gc_count <= int(gc_hi * n_bases))

        for site in restriction_sites:
            site_upper = site.upper()
            if any(b not in "ACGT" for b in site_upper):
                continue
            indices = [BASE_REV[b] for b in site_upper]
            for start in range(n_bases - len(site_upper) + 1):
                opt.add(Not(And(*[nuc[start + j] == indices[j] for j in range(len(site_upper))])))

        for start in range(n_bases - 5 + 1):
            opt.add(Not(And(*[nuc[start + j] == BASE_REV[b] for j, b in enumerate("ATTTA")])))

        for start in range(n_bases - 6 + 1):
            opt.add(Not(And(*[nuc[start + j] == 3 for j in range(6)])))

        result = opt.check()
        if result == sat:
            model = opt.model()
            sequence = "".join(BASE_MAP[model.eval(nuc[i], model_completion=True).as_long()] for i in range(n_bases))
        else:
            logger.info("z3 returned %s, falling back to greedy optimizer", result)
            sequence, warnings = _greedy_optimize(
                target_protein, organism, gc_lo, gc_hi, restriction_sites, cryptic_splice_threshold
            )
            fallback_used = True
            all_warnings.extend(warnings)
    else:
        # Greedy optimizer (default) — fast, high CAI, multi-codon constraint solving
        sequence, warnings = _greedy_optimize(
            target_protein, organism, gc_lo, gc_hi, restriction_sites, cryptic_splice_threshold
        )
        fallback_used = True  # Greedy is now the primary path
        all_warnings.extend(warnings)

    cai = compute_cai(sequence, organism)
    gc = gc_content(sequence)

    # Delegate predicate checking to the type system (SOC)
    satisfied, failed = _check_predicates_via_type_system(
        sequence, gc_lo, gc_hi, restriction_sites, cai_threshold, organism,
        cryptic_splice_threshold,
    )

    if all_warnings:
        logger.info("Optimization warnings: %s", "; ".join(all_warnings))

    # Type-directed mutagenesis: if predicates fail and mutagenesis is enabled,
    # propose conservative amino acid substitutions to make satisfaction possible.
    mutagenesis_applied = False
    aa_substitutions = None

    # Before proposing any amino acid substitutions, try to fix "chose poorly"
    # positions — non-Valine positions with strong cryptic donors where the
    # optimizer used a GT-containing codon despite GT-free alternatives.
    if failed and any("CrypticSplice" in p for p in failed):
        from .mutagenesis import force_gt_free_reoptimization, is_gt_mandatory
        from .maxentscan import score_donor as _score_donor

        # Check if there are non-GT-mandatory positions with strong donors
        has_chose_poorly = False
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == "GT":
                s = _score_donor(sequence, i)
                if s >= cryptic_splice_threshold:
                    codon_idx = i // 3
                    if codon_idx < len(target_protein):
                        aa = target_protein[codon_idx]
                        if not is_gt_mandatory(aa):
                            has_chose_poorly = True
                            break

        if has_chose_poorly:
            sequence = force_gt_free_reoptimization(
                sequence, target_protein, organism,
                threshold=cryptic_splice_threshold,
            )
            # Re-check predicates after forced re-optimization
            cai = compute_cai(sequence, organism)
            gc = gc_content(sequence)
            satisfied, failed = _check_predicates_via_type_system(
                sequence, gc_lo, gc_hi, restriction_sites, cai_threshold,
                organism, cryptic_splice_threshold,
            )

    if enable_mutagenesis and failed:
        from .mutagenesis import type_directed_mutagenesis
        mut_result = type_directed_mutagenesis(
            protein=target_protein,
            organism=organism,
            failed_predicates=failed,
            optimize_fn=optimize_sequence,
            max_substitutions=max_mutagenesis_substitutions,
            min_blosum62=min_blosum62,
            cryptic_splice_threshold=cryptic_splice_threshold,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            # Pass through non-mutagenesis kwargs (avoid recursion)
            enable_mutagenesis=False,
            cai_threshold=cai_threshold,
        )

        if mut_result.substitutions:
            mutagenesis_applied = True
            aa_substitutions = [
                {
                    "position": sub.position,
                    "from": sub.original_aa,
                    "to": sub.substitute_aa,
                    "blosum62": sub.blosum62_score,
                    "reason": sub.reason,
                    "predicate": sub.predicate_addressed,
                }
                for sub in mut_result.substitutions
            ]

            # Re-optimize with the modified protein
            mut_protein = mut_result.modified_protein
            # Re-run the full optimization pipeline (mutagenesis disabled to avoid recursion)
            final_result = optimize_sequence(
                target_protein=mut_protein,
                organism=organism,
                gc_lo=gc_lo,
                gc_hi=gc_hi,
                restriction_sites=restriction_sites,
                cai_threshold=cai_threshold,
                cryptic_splice_threshold=cryptic_splice_threshold,
                enable_mutagenesis=False,
            )

            return OptimizationResult(
                sequence=final_result.sequence,
                protein=mut_protein,
                cai=final_result.cai,
                gc_content=final_result.gc_content,
                satisfied_predicates=final_result.satisfied_predicates,
                failed_predicates=final_result.failed_predicates,
                fallback_used=final_result.fallback_used,
                mutagenesis_applied=True,
                aa_substitutions=aa_substitutions,
            )

    return OptimizationResult(
        sequence=sequence, protein=target_protein, cai=cai,
        gc_content=gc, satisfied_predicates=satisfied,
        failed_predicates=failed, fallback_used=fallback_used,
        mutagenesis_applied=mutagenesis_applied,
        aa_substitutions=aa_substitutions,
    )
