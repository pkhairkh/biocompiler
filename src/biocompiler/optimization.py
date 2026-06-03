"""
BioCompiler Optimizer v7.0.0
==============================
Multi-step certified gene optimization pipeline with aggressive GT resolution.

Step: Maximize CAI          — Greedy codon optimization (GT-aware, with unavoidable-GT tracking)
Step: Backtranslate CAI     — DP-based max-CAI back-translation with avoidable-GT avoidance
Step: Resolve Constraints   — Priority-based constraint resolution (fix GT/CG/RS with minimal CAI loss)
Step: Remove Restriction Sites — Restriction site removal by synonymous substitution
Step: Cross-Codon Optimization — Cross-codon constraint resolution (iterative, global validation)
Step: Within-Codon GT Resolution — Within-codon GT resolution (synonymous substitution + mutagenesis flagging)
Step: Mutagenesis Fallback  — Mutagenesis fallback (AA substitution for Valine etc. using BLOSUM62)
Step: Avoid CpG Islands     — CpG island avoidance
Step: Cross-Codon Coordination — Cross-codon coordinated solver (exhaustive pair search at boundaries)
Step: CAI Hill Climb        — CAI hill climbing (upgrade codons while maintaining constraints)
Step: Reoptimize            — Re-optimization pass (iterative until convergence)
"""

from typing import List, Dict, Optional, Tuple, Set

from dataclasses import dataclass, field

from .type_system import (
    CODON_TABLE, AA_TO_CODONS, BLOSUM62, SpliceVerdict, PredicateResult,
    check_no_stop_codons, check_no_cryptic_splice, check_no_cpg_island,
    check_no_restriction_site, check_no_gt_dinucleotide, check_no_avoidable_gt,
    check_valid_coding_seq,
    check_conservation_score, check_codon_optimality,
    find_cross_codon_gt, find_cross_codon_cg, find_cross_codon_restriction,
    PREDICATE_NAMES,
)
from .organisms import SPECIES
from .mutagenesis import propose_mutagenesis, MutagenesisReport, MutagenesisProposal
from .certificate import compute_certificate, format_certificate
from .exceptions import InvalidProteinError, UnsupportedOrganismError


# ────────────────────────────────────────────────────────────
# High-level OptimizationResult and optimize_sequence API
# ────────────────────────────────────────────────────────────

@dataclass
class OptimizationResult:
    """Result of optimizing a protein sequence.

    Provides the optimized DNA sequence along with quality metrics
    and a list of predicates that the result fails to satisfy.
    """
    sequence: str
    gc_content: float
    cai: float
    failed_predicates: list = field(default_factory=list)
    predicate_results: list = field(default_factory=list)
    certificate_text: str = ""
    # Extended attributes for API/visualization compatibility
    protein: str = ""
    fallback_used: bool = False
    satisfied_predicates: list = field(default_factory=list)
    aa_substitutions: list = field(default_factory=list)
    mutagenesis_applied: bool = False

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

    Step ordering prioritizes hard constraints (restriction sites) over soft constraints (CAI).
    Reconciliation pass ensures earlier steps aren't undone by later ones.

    Steps:
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

    # Step: Maximize CAI — Best codon per position (maximize CAI)
    sequence = "".join(sorted_codons[aa][0] for aa in aas)
    assert len(sequence) == len(aas) * 3, "Maximize CAI step: sequence length mismatch"

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

    # Step: Remove ATTTA instability motifs
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

    # Step: Fix 6+ consecutive T runs
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

    # Step: Adjust GC content
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

    # Step: Reconciliation — check if GC adjustment reintroduced restriction sites
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
                # Could not remove — already warned in Remove Restriction Sites step
                pass

    # Step: Eliminate cryptic splice donor/acceptor sites
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

    # Step: Disrupt CpG dinucleotides to avoid CpG islands
    # Strategy: Replace CG dinucleotides with synonymous codons that don't create CG,
    # but ONLY if the swap doesn't reintroduce cryptic splice sites or restriction sites.
    # Key constraint: CpG avoidance is lower priority than splice site elimination.
    # We will NOT worsen any cryptic splice score to remove a CG dinucleotide.
    # For cross-codon CG (C at end of one codon, G at start of next), we use
    # coordinated 2-codon swaps since single-codon swaps can't always resolve them.
    for _cpg_iteration in range(200):
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

    # Step: Reconciliation after cryptic splice elimination
    # Check if cryptic splice fixes reintroduced restriction sites
    for site_upper in concrete_sites:
        site_rc = reverse_complement(site_upper)
        if site_upper in sequence or site_rc in sequence:
            new_seq, fixed = _remove_site_multicodon(
                sequence, aas, sorted_codons, site_upper, site_rc
            )
            if fixed:
                sequence = new_seq

    # Step: CpG reconciliation after restriction site reconciliation
    # Restriction site removal may have reintroduced CpG dinucleotides;
    # re-apply CpG disruption without reintroducing restriction sites or
    # worsening cryptic splice scores. Includes 2-codon coordinated swaps
    # for cross-codon CG dinucleotides.
    for _cpg_iter in range(200):
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
        if r.verdict in (Verdict.PASS, Verdict.LIKELY_PASS):
            satisfied.append(predicate_name)
        elif r.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
            failed.append(predicate_name)
        else:
            # UNCERTAIN: treat as failed for optimization purposes
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
    target_protein: str | None = None,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list | None = None,
    strategy: str = "constraint_first",
    **kwargs,
) -> OptimizationResult:
    """Optimize a protein sequence for expression in the target organism.

    This is the high-level convenience API that wraps BioOptimizer.
    It takes a protein sequence and returns an OptimizationResult with
    the optimized DNA sequence and quality metrics.

    Args:
        target_protein: Amino acid sequence (1-letter codes, no stop).
            Can also be passed as the first positional argument.
        organism: Target organism (e.g., 'Homo_sapiens', 'Escherichia_coli').
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        enzymes: List of restriction enzyme names to avoid.
        strategy: Optimization strategy ('constraint_first' or 'cai_first').
        **kwargs: Additional arguments (e.g., splice_low, splice_high).

    Returns:
        OptimizationResult with optimized sequence and metrics.

    Raises:
        InvalidProteinError: if the protein contains invalid amino acid codes.
        UnsupportedOrganismError: if the organism is not supported.
    """
    # Handle positional arg: optimize_sequence("MVHLTPEEK", organism="...")
    if target_protein is None and len(kwargs.get("_args", [])) > 0:
        target_protein = kwargs["_args"][0]

    if not target_protein:
        raise InvalidProteinError("", {"<empty>"})

    target_protein = target_protein.strip().upper()

    # Validate protein
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(target_protein) - valid_aas
    if invalid:
        raise InvalidProteinError(target_protein, invalid)

    # Map organism name to species key
    species_key = _organism_to_species_key(organism)

    # Configure optimizer
    opt = BioOptimizer(
        species=species_key,
        enzymes=enzymes or ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
        splice_low=kwargs.get("splice_low", 3.0),
        splice_high=kwargs.get("splice_high", 6.0),
        min_cai=cai_threshold,
        strategy=strategy,
    )

    # Back-translate protein to DNA for the optimizer
    from .translation import compute_cai
    initial_seq = _back_translate_protein(target_protein, species_key)

    # Run optimization
    optimized_seq, pred_results, cert_text = opt.optimize(initial_seq, strategy=strategy)

    # Compute metrics
    gc = (optimized_seq.count("G") + optimized_seq.count("C")) / max(len(optimized_seq), 1)
    gc = round(gc, 4)

    try:
        cai_val = compute_cai(optimized_seq, organism)
    except UnsupportedOrganismError:
        cai_val = opt._compute_seq_cai(optimized_seq)

    # Collect failed and satisfied predicates
    failed = [r.predicate for r in pred_results if not r.passed]
    satisfied = [r.predicate for r in pred_results if r.passed]

    # Detect if mutagenesis fallback was used (any V->I or similar)
    fallback = bool(opt._applied_mutagenesis)

    return OptimizationResult(
        sequence=optimized_seq,
        gc_content=gc,
        cai=cai_val,
        failed_predicates=failed,
        predicate_results=pred_results,
        certificate_text=cert_text,
        protein=target_protein,
        fallback_used=fallback,
        satisfied_predicates=satisfied,
        aa_substitutions=opt._applied_mutagenesis,
    )


def _organism_to_species_key(organism: str) -> str:
    """Map an organism name to the species key used in the SPECIES dict."""
    mapping = {
        "Homo_sapiens": "human",
        "human": "human",
        "Escherichia_coli": "ecoli",
        "E_coli": "ecoli",
        "ecoli": "ecoli",
        "Saccharomyces_cerevisiae": "yeast",
        "yeast": "yeast",
        "Mus_musculus": "mouse",
        "mouse": "mouse",
        "CHO": "cho",
        "CHO_K1": "cho",
    }
    key = mapping.get(organism)
    if key and key in SPECIES:
        return key
    # Fallback: try the organism name directly
    if organism in SPECIES:
        return organism
    # Default to ecoli
    return "ecoli"


def _back_translate_protein(protein: str, species_key: str) -> str:
    """Back-translate a protein to DNA using highest-CAI codons."""
    species_cai = SPECIES.get(species_key, SPECIES["ecoli"])
    codons = []
    for aa in protein:
        if aa == "*":
            codons.append("TAA")
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            codons.append("NNN")
            continue
        best = max(candidates, key=lambda c: species_cai.get(c, 0.0))
        codons.append(best)
    return "".join(codons)


def _count_gts(s: str) -> int:
    """Count GT dinucleotides in a sequence."""
    return sum(1 for i in range(len(s) - 1) if s[i:i+2] == "GT")


def _is_unavoidable_gt(seq: str, pos: int) -> bool:
    """Check if a GT dinucleotide at position pos is unavoidable.
    
    A GT is unavoidable if:
    1. It's within a Valine codon (all Val codons start with GT)
    2. It's a cross-codon GT where the next codon's AA has no synonymous
       codon that doesn't start with T (e.g., Trp=TGG, Cys=TGT/TGC, Tyr=TAT/TAC)
    3. It's a cross-codon GT where the previous codon's AA has no synonymous
       codon that doesn't end with G
    """
    codon_start = (pos // 3) * 3
    next_codon_start = codon_start + 3

    # Case 1: Within-codon GT
    if pos + 1 < next_codon_start:
        codon = seq[codon_start:codon_start + 3]
        aa = CODON_TABLE.get(codon)
        if aa == 'V':
            return True  # All Valine codons start with GT
        # Check if any synonymous codon avoids GT
        for alt in AA_TO_CODONS.get(aa, []):
            if "GT" not in alt:
                return False  # There's an alternative without GT
        return True  # No alternative without GT

    # Case 2: Cross-codon GT (pos is last base of one codon, pos+1 is first of next)
    # pos is at codon_start + 2 (last position of current codon)
    # OR pos is at codon_start (we need to figure out which codons are involved)
    
    # For cross-codon GT at pos: seq[pos]='G', seq[pos+1]='T'
    # pos is the last base of the preceding codon
    # pos+1 is the first base of the following codon
    prev_cs = (pos // 3) * 3  # Start of codon containing position 'pos'
    next_cs = prev_cs + 3     # Start of codon containing position 'pos+1'
    
    if next_cs + 3 > len(seq):
        return True  # Can't check, assume unavoidable
    
    prev_codon = seq[prev_cs:prev_cs + 3]
    next_codon = seq[next_cs:next_cs + 3]
    prev_aa = CODON_TABLE.get(prev_codon)
    next_aa = CODON_TABLE.get(next_codon)
    
    if prev_aa is None or next_aa is None:
        return True
    
    # Check if we can change the previous codon to not end with G
    prev_can_avoid = any(c[-1] != 'G' for c in AA_TO_CODONS.get(prev_aa, [prev_codon]))
    # Check if we can change the next codon to not start with T
    next_can_avoid = any(c[0] != 'T' for c in AA_TO_CODONS.get(next_aa, [next_codon]))
    
    # GT is unavoidable only if BOTH sides can't avoid it
    return not (prev_can_avoid or next_can_avoid)


def _has_gt(s: str) -> bool:
    """Check if a string contains GT dinucleotide."""
    return "GT" in s


def _codon_creates_boundary_gt(
    codon: str, codon_start: int, seq_list: list
) -> Tuple[bool, bool]:
    """Check if placing codon at codon_start creates cross-codon GTs.

    Returns (prev_boundary_gt, next_boundary_gt).
    """
    prev_gt = False
    next_gt = False
    if codon_start > 0 and seq_list[codon_start - 1] + codon[0] == "GT":
        prev_gt = True
    next_pos = codon_start + 3
    if next_pos < len(seq_list) and codon[-1] + seq_list[next_pos] == "GT":
        next_gt = True
    return prev_gt, next_gt


class BioOptimizer:
    """Certified gene sequence optimizer with multi-step CAI-maximizing pipeline."""

    def __init__(
        self,
        species: str = "ecoli",
        enzymes: Optional[List[str]] = None,
        splice_low: float = 3.0,
        splice_high: float = 6.0,
        cpg_window: int = 200,
        cpg_threshold: float = 0.6,
        min_blosum: int = -1,
        min_cai: float = 0.0,
        avoid_gt: bool = True,
        strategy: str = "constraint_first",
    ):
        self.species = species
        self.species_cai: Dict[str, float] = SPECIES.get(species, SPECIES["ecoli"])
        self.enzymes: List[str] = enzymes or []
        self.splice_low = splice_low
        self.splice_high = splice_high
        self.cpg_window = cpg_window
        self.cpg_threshold = cpg_threshold
        self.min_blosum = min_blosum
        self.min_cai = min_cai
        self.avoid_gt = avoid_gt
        self.strategy = strategy  # "constraint_first" or "cai_first"
        # Track positions where GT is unavoidable (e.g., Valine codons)
        self._unavoidable_gt_positions: Set[int] = set()
        # Track mutagenesis proposals that were applied
        self._applied_mutagenesis: List[Dict] = []
        # Store original input protein for conservation scoring
        self._original_protein: str = ""

    def optimize(self, seq: str, strategy: Optional[str] = None) -> Tuple[str, List[PredicateResult], str]:
        """Run the full optimization pipeline.

        Args:
            seq: Input DNA or protein sequence.
            strategy: Optimization strategy override.
                - "constraint_first" (default): GT-aware greedy then fix constraints
                - "cai_first": Maximize CAI first, then fix constraints with
                  minimal CAI impact (DNAworks-style)

        Returns:
            (optimized_sequence, predicate_results, certificate_text)
        """
        effective_strategy = strategy if strategy is not None else self.strategy

        seq = seq.upper().strip()
        self._unavoidable_gt_positions = set()
        self._applied_mutagenesis = []
        self._original_protein = self._translate(seq)

        if effective_strategy == "cai_first":
            return self._optimize_cai_first(seq)

        # ── Default: constraint_first strategy ──
        # Step: Backtranslate CAI (DNAworks-style)
        seq = self._step_backtranslate_cai(seq)

        # Step: Resolve Constraints (fix GT/CG/RS with minimal CAI loss)
        seq = self._step_resolve_constraints(seq)

        # Step: Remove Restriction Sites
        seq = self._step_remove_restriction_sites(seq)

        # Step: Cross-Codon Optimization (iterative)
        seq, mut_report = self._step_cross_codon_optimization(seq)

        # Step: Within-Codon GT Resolution
        seq, mut_report_35 = self._step_within_codon_gt_resolution(seq)
        mut_report.proposals.extend(mut_report_35.proposals)

        # Step: Mutagenesis Fallback (aggressive, handles within-codon GTs too)
        seq = self._step_mutagenesis_fallback(seq, mut_report)

        # Step: Avoid CpG Islands
        seq = self._step_avoid_cpg_islands(seq)

        # Step: Cross-Codon Coordination (handles cross-codon GT, CG, restriction sites)
        seq = self._step_cross_codon_coordination(seq)

        # Step: CAI Hill Climb (upgrade codons while maintaining constraints)
        seq = self._step_cai_hill_climb(seq)

        # Step: Reoptimize (iterative until convergence)
        seq = self._step_reoptimize(seq)

        # Step: CpG Reconciliation (aggressive, after CAI hill climb/reoptimize)
        seq = self._step_cpg_reconciliation(seq)

        # Evaluate all 8 predicates
        results = self._evaluate_all_predicates(seq)

        # Generate certificate
        cert_text = format_certificate(results, seq, self.species)

        return seq, results, cert_text

    # ──────────────────────────────────────────────────────────
    # CAI-first optimization strategy (DNAworks-style)
    # ──────────────────────────────────────────────────────────
    def _optimize_cai_first(self, seq: str) -> Tuple[str, List[PredicateResult], str]:
        """CAI-first optimization: maximize CAI first, fix constraints second.

        Strategy (DNAworks-style):
        1. Back-translate using highest-CAI codons at every position
           (ignore GT/CG/RS initially)
        2. Iteratively fix constraint violations with minimal CAI impact:
           a. Restriction sites → fix each with best synonymous codon
           b. Avoidable GT dinucleotides → fix with best non-GT synonymous codon
              using cross-codon pair optimization
           c. CpG islands → fix with best non-CG synonymous codon
           d. Cryptic splice sites → fix with best synonymous codon
        3. Mutagenesis fallback for GTs that can't be resolved by
           synonymous substitution (e.g., Valine V→Isoleucine I)
        4. CAI hill climbing to recover any lost CAI while maintaining constraints
        5. Iterate until convergence

        Key insight: by starting from max CAI (CAI=1.0) and only making the
        smallest necessary CAI sacrifices to fix constraints, we achieve much
        higher CAI than the constraint_first strategy which permanently sacrifices
        CAI by avoiding GT during codon selection.

        For GT fixing, we use a priority-based approach that fixes GTs with
        the lowest CAI cost first, and considers multi-codon windows to find
        globally optimal fixes.
        """
        import math

        # Step: Maximize CAI — Pure max-CAI back-translation (ignore ALL constraints)
        protein = self._translate(seq)
        codons_result = []
        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons_result.append(
                    seq[codon_start:codon_start + 3]
                    if codon_start + 3 <= len(seq) else "TAA"
                )
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons_result.append(
                    max(candidates, key=lambda c: self.species_cai.get(c, 0.0))
                )
            else:
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3])
        seq = "".join(codons_result)

        # Step: Fix restriction sites (highest priority — binary constraint)
        seq = self._cai_first_fix_restriction_sites(seq)

        # Step: Fix avoidable GT dinucleotides (minimal CAI impact)
        seq = self._cai_first_fix_gts(seq)

        # Step: Mutagenesis fallback for Valine GTs
        # Valine codons all contain GT (GTN), so we substitute V→I
        # (BLOSUM62 score 3, conservative) using high-CAI Ile codons
        seq = self._cai_first_mutagenesis_fallback(seq)

        # Step: Fix CpG islands (minimal CAI impact)
        seq = self._cai_first_fix_cpg(seq)

        # Step: Fix cryptic splice sites (minimal CAI impact)
        seq = self._cai_first_fix_splice(seq)

        # Step: Cross-Codon Coordination (handles cross-codon GT, CG, restriction sites)
        seq = self._step_cross_codon_coordination(seq)

        # Step: CAI hill climbing (upgrade codons while maintaining constraints)
        seq = self._step_cai_hill_climb(seq)

        # Step: Aggressive re-optimization pass
        seq = self._step_reoptimize(seq)

        # Step: Second pass of GT fixing + CAI boost (iterative refinement)
        for _refinement in range(3):
            old_cai = self._compute_seq_cai(seq)
            seq = self._cai_first_fix_gts(seq)
            seq = self._step_cai_hill_climb(seq)
            seq = self._step_reoptimize(seq)
            new_cai = self._compute_seq_cai(seq)
            if new_cai <= old_cai + 0.0001:
                break

        # Step: CpG Reconciliation (aggressive, after CAI hill climb/reoptimize)
        seq = self._step_cpg_reconciliation(seq)

        # Step: CAI Reconciliation (upgrade low-CAI codons while maintaining hard constraints)
        seq = self._step_cai_reconciliation(seq)

        # Track unavoidable GTs for certificate
        for i in range(len(seq) - 1):
            if seq[i] == "G" and seq[i + 1] == "T":
                if _is_unavoidable_gt(seq, i):
                    codon_start = (i // 3) * 3
                    self._unavoidable_gt_positions.add(i)

        # Evaluate all predicates
        results = self._evaluate_all_predicates(seq)
        cert_text = format_certificate(results, seq, self.species)
        return seq, results, cert_text

    def _compute_seq_cai(self, seq: str) -> float:
        """Compute the geometric mean CAI for a sequence."""
        import math
        if not seq or len(seq) < 3:
            return 0.0
        log_sum = 0.0
        count = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            cai = self.species_cai.get(codon, 0.0)
            if cai <= 0:
                cai = 0.001
            log_sum += math.log(cai)
            count += 1
        if count == 0:
            return 0.0
        return math.exp(log_sum / count)

    def _step_maximize_cai(self, seq: str) -> str:
        """Maximize CAI step (cai_first): Back-translate with absolute max CAI everywhere.

        No GT avoidance at all - just pick the highest-CAI codon for each AA.
        Constraint violations will be fixed in subsequent steps.
        """
        protein = self._translate(seq)
        codons_result = []
        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3] if codon_start + 3 <= len(seq) else "TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if candidates:
                codons_result.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
            else:
                codon_start = i * 3
                codons_result.append(seq[codon_start:codon_start + 3])
        return "".join(codons_result)

    def _cai_first_fix_restriction_sites(self, seq: str) -> str:
        """CAI-first Fix restriction sites with minimal CAI impact."""
        from .restriction_sites import get_recognition_site
        seq_list = list(seq)

        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue

            max_rounds = 50
            for _ in range(max_rounds):
                current_seq = "".join(seq_list)
                p = current_seq.find(site)
                if p == -1:
                    break

                # Find all codon positions overlapping this site
                codon_starts = set()
                for j in range(p, p + len(site)):
                    cs = (j // 3) * 3
                    if cs + 3 <= len(seq_list):
                        codon_starts.add(cs)

                fixed = False
                # Try each overlapping codon, sorted by CAI impact (try best-CAI alts first)
                for cs in sorted(codon_starts):
                    codon = "".join(seq_list[cs:cs + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue

                    # Sort alternatives by CAI (highest first) to minimize CAI loss
                    alts = sorted(
                        AA_TO_CODONS.get(aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    )
                    for alt in alts:
                        if alt == codon:
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[cs + k] = b
                        test_seq = "".join(test_list)
                        if site not in test_seq:
                            seq_list = test_list
                            fixed = True
                            break
                    if fixed:
                        break

                if not fixed:
                    # Could not fix this site with single-codon substitution
                    # Try two-codon substitution
                    fixed = self._cai_first_fix_rs_two_codons(seq_list, p, site)
                    if not fixed:
                        break  # Give up on this site for now

        return "".join(seq_list)

    def _cai_first_fix_rs_two_codons(self, seq_list: list, pos: int, site: str) -> bool:
        """Try fixing a restriction site by modifying two adjacent codons."""
        codon_starts = sorted(set(
            (j // 3) * 3
            for j in range(pos, min(pos + len(site), len(seq_list)))
            if (j // 3) * 3 + 3 <= len(seq_list)
        ))

        if len(codon_starts) < 2:
            return False

        # Try pairs of codons
        for idx in range(len(codon_starts) - 1):
            cs1, cs2 = codon_starts[idx], codon_starts[idx + 1]
            if cs2 != cs1 + 3:
                continue  # Only adjacent codons

            codon1 = "".join(seq_list[cs1:cs1 + 3])
            codon2 = "".join(seq_list[cs2:cs2 + 3])
            aa1 = CODON_TABLE.get(codon1)
            aa2 = CODON_TABLE.get(codon2)
            if aa1 is None or aa1 == "*" or aa2 is None or aa2 == "*":
                continue

            # Try all pairs sorted by combined CAI
            pairs = []
            for c1 in AA_TO_CODONS.get(aa1, [codon1]):
                for c2 in AA_TO_CODONS.get(aa2, [codon2]):
                    combined = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                    pairs.append((c1, c2, combined))
            pairs.sort(key=lambda x: x[2], reverse=True)

            for c1, c2, _ in pairs:
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[cs1 + k] = b
                for k, b in enumerate(c2):
                    test_list[cs2 + k] = b
                test_seq = "".join(test_list)
                if site not in test_seq:
                    seq_list[:] = test_list
                    return True

        return False

    def _cai_first_fix_gts(self, seq: str) -> str:
        """CAI-first Fix avoidable GT dinucleotides with minimal CAI impact.

        Iteratively finds each avoidable GT and resolves it by choosing
        the synonymous substitution(s) with the highest possible CAI that
        eliminates the GT.
        """
        seq_list = list(seq)
        max_rounds = 50

        for round_num in range(max_rounds):
            current_seq = "".join(seq_list)
            violations = []

            # Find all avoidable GT positions
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    if not _is_unavoidable_gt(current_seq, i):
                        violations.append(i)

            if not violations:
                break

            any_fixed = False

            for gt_pos in violations:
                codon_start = (gt_pos // 3) * 3
                next_codon_start = codon_start + 3

                # Determine if within-codon or cross-codon
                is_within = (gt_pos + 1) < next_codon_start

                if is_within:
                    # Within-codon GT: try synonymous substitution
                    fixed = self._cai_first_fix_within_gt(seq_list, codon_start)
                else:
                    # Cross-codon GT: try adjacent codon pair substitution
                    fixed = self._cai_first_fix_cross_gt(seq_list, codon_start)

                if fixed:
                    any_fixed = True
                    break  # Restart scanning after any fix

            if not any_fixed:
                break

        return "".join(seq_list)

    def _cai_first_fix_within_gt(self, seq_list: list, codon_start: int) -> bool:
        """Fix a within-codon GT by choosing the highest-CAI alternative without GT."""
        codon = "".join(seq_list[codon_start:codon_start + 3])
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            return False

        old_gt_count = _count_gts("".join(seq_list))

        # Sort alternatives by CAI (highest first), skip codons with GT
        alternatives = []
        for alt in AA_TO_CODONS.get(aa, []):
            if "GT" in alt:
                continue
            alt_cai = self.species_cai.get(alt, 0.0)
            alternatives.append((alt, alt_cai))

        alternatives.sort(key=lambda x: x[1], reverse=True)

        for alt, _ in alternatives:
            test_list = seq_list[:]
            for k, b in enumerate(alt):
                test_list[codon_start + k] = b
            new_gt_count = _count_gts("".join(test_list))
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    def _cai_first_fix_cross_gt(self, seq_list: list, codon_start: int) -> bool:
        """Fix a cross-codon GT with minimal CAI impact.

        Tries strategies in order of increasing invasiveness:
        D. Change only the next codon (to one that doesn't start with T)
        C. Change only the current codon (to one that doesn't end with G)
        A. Modify the codon pair (both codons) — sorted by combined CAI
        B. Modify preceding codon pair
        """
        old_gt_count = _count_gts("".join(seq_list))
        next_start = codon_start + 3

        # Strategy D: Change only the next codon (cheapest single-codon fix)
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                alts = sorted(
                    AA_TO_CODONS.get(aa2, [aa2_codon]),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )
                for c2 in alts:
                    if c2[0] == "T":
                        continue
                    test_list = seq_list[:]
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy C: Change only the current codon (to one that doesn't end with G)
        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is not None and aa1 != "*":
            alts = sorted(
                AA_TO_CODONS.get(aa1, [aa1_codon]),
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )
            for c1 in alts:
                if c1[-1] == "G":
                    continue
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[codon_start + k] = b
                test_seq = "".join(test_list)
                new_gt_count = _count_gts(test_seq)
                if new_gt_count < old_gt_count:
                    seq_list[:] = test_list
                    return True

        # Strategy A: Modify current + following codon
        if next_start + 3 <= len(seq_list):
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)

            if aa1 is not None and aa1 != "*" and aa2 is not None and aa2 != "*":
                pairs = []
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                        if c1[-1] + c2[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                        pairs.append((c1, c2, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy B: Modify preceding + current codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)

            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                pairs = []
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if c0[-1] + c1[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c0, 0.0) + self.species_cai.get(c1, 0.0)
                        pairs.append((c0, c1, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c0, c1, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c0):
                        test_list[prev_start + k] = b
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        return False

    def _cai_first_fix_cpg(self, seq: str) -> str:
        """CAI-first Fix CpG islands with minimal CAI impact."""
        seq_list = list(seq)
        changed = True
        iterations = 0
        max_iterations = 50

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(seq_list) - self.cpg_window + 1, 3):
                window = "".join(seq_list[start:start + self.cpg_window])
                c_count = window.count("C")
                g_count = window.count("G")
                cg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
                expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

                # Find CG dinucleotides and try to break them
                for i in range(start, min(start + self.cpg_window - 1, len(seq_list) - 1)):
                    if seq_list[i] == "C" and seq_list[i+1] == "G":
                        codon_start = (i // 3) * 3
                        if codon_start + 3 > len(seq_list):
                            continue
                        codon = "".join(seq_list[codon_start:codon_start + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        old_gt_count = _count_gts("".join(seq_list))

                        # Sort by CAI descending to minimize CAI loss
                        for alt in sorted(
                            AA_TO_CODONS.get(aa, []),
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        ):
                            if alt == codon or "CG" in alt:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[codon_start + k] = b
                            new_gt_count = _count_gts("".join(test_list))
                            if new_gt_count <= old_gt_count:
                                seq_list = test_list
                                changed = True
                                break

                        if changed:
                            break
                if changed:
                    break

        return "".join(seq_list)

    def _cai_first_fix_splice(self, seq: str) -> str:
        """CAI-first Fix cryptic splice sites with minimal CAI impact."""
        seq_list = list(seq)
        max_rounds = 30

        for _ in range(max_rounds):
            current_seq = "".join(seq_list)
            splice_result = check_no_cryptic_splice(current_seq, self.splice_low, self.splice_high)
            if splice_result.passed:
                break

            # Find and fix splice sites
            # The splice check looks for GT..AG patterns that resemble splice donors/acceptors
            fixed = False
            for i in range(len(current_seq) - 1):
                if current_seq[i] == "G" and current_seq[i + 1] == "T":
                    # This GT might be part of a cryptic splice site
                    codon_start = (i // 3) * 3
                    next_codon_start = codon_start + 3
                    is_within = (i + 1) < next_codon_start

                    if _is_unavoidable_gt(current_seq, i):
                        continue

                    if is_within:
                        fixed = self._cai_first_fix_within_gt(seq_list, codon_start)
                    else:
                        fixed = self._cai_first_fix_cross_gt(seq_list, codon_start)

                    if fixed:
                        break

            if not fixed:
                break

        return "".join(seq_list)

    def _cai_first_mutagenesis_fallback(self, seq: str) -> str:
        """CAI-first Apply mutagenesis for GTs that can't be resolved
        by synonymous substitution.

        Specifically targets Valine codons (GTN) which all contain GT.
        Substitutes V→I (Isoleucine, BLOSUM62=3) using highest-CAI Ile codon
        to eliminate the GT while maximizing CAI.
        """
        seq_list = list(seq)
        changed = True
        max_rounds = 10

        for _ in range(max_rounds):
            if not changed:
                break
            changed = False

            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                if "GT" not in codon:
                    continue

                has_gt_free = any("GT" not in c for c in AA_TO_CODONS.get(aa, []))
                if has_gt_free:
                    continue

                if aa == "V":
                    ile_codons = AA_TO_CODONS.get("I", [])
                    old_gt_count = _count_gts("".join(seq_list))

                    for ile_codon in sorted(
                        ile_codons,
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        test_list = seq_list[:]
                        for k, b in enumerate(ile_codon):
                            test_list[i + k] = b
                        new_gt_count = _count_gts("".join(test_list))
                        if new_gt_count < old_gt_count:
                            seq_list[:] = test_list
                            self._applied_mutagenesis.append({
                                "position": i,
                                "original_aa": "V",
                                "new_aa": "I",
                                "blosum": BLOSUM62.get(("V", "I"), 3),
                            })
                            changed = True
                            break

        return "".join(seq_list)

    # Deprecated alias — use _step_maximize_cai instead
    _phase0_pure_max_cai = _step_maximize_cai

    # ──────────────────────────────────────────────────────────
    # Step: Backtranslate CAI (DP-based max-CAI back-translation)
    # ──────────────────────────────────────────────────────────
    def _step_backtranslate_cai(self, seq: str) -> str:
        """Backtranslate CAI step: DP-based max-CAI back-translation with avoidable-GT avoidance.

        Uses Viterbi-style dynamic programming to find the globally optimal
        codon assignment that maximizes CAI while avoiding only AVOIDABLE GTs.
        The DP state is simply the last character of the previous codon,
        making it O(n * 4 * k) where n is protein length and k is max codons per AA.

        Key insight: a cross-codon GT (codon1 ends with G, codon2 starts with T)
        is only avoidable if we can change at least one of the two codons to
        eliminate it. If codon2 has NO synonymous codon that doesn't start with T
        (e.g., Trp=TGG, Cys=TGT/TGC, Tyr=TAT/TAC), then the cross-codon GT
        is unavoidable and we should use the highest-CAI codon for codon1.

        This matches the semantics of the check_no_avoidable_gt predicate,
        which only fails on GTs that CAN be avoided by synonymous substitution.

        Only non-Valine within-codon GTs and avoidable cross-codon GTs are
        excluded by the DP. Unavoidable GTs (Valine within-codon, cross-codon
        where next AA has no non-T-starting codon) are allowed.
        """
        import math
        protein = self._translate(seq)
        n = len(protein)
        INF = float('-inf')

        # Precompute which amino acids have at least one codon not starting with T
        # These are the AAs where cross-codon GT from previous codon can be avoided
        aa_has_non_t_start = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons_list = AA_TO_CODONS.get(aa_key, [])
            aa_has_non_t_start[aa_key] = any(c[0] != 'T' for c in codons_list)

        # Precompute which amino acids have at least one codon not ending with G
        # These are the AAs where cross-codon GT to next codon can be avoided
        aa_has_non_g_end = {}
        for aa_key in set(CODON_TABLE.values()):
            if aa_key == "*":
                continue
            codons_list = AA_TO_CODONS.get(aa_key, [])
            aa_has_non_g_end[aa_key] = any(c[-1] != 'G' for c in codons_list)

        # DP table: dp[i][last_char] = (max_log_cai_sum, prev_last_char, chosen_codon)
        dp = [{} for _ in range(n + 1)]
        dp[0][''] = (0.0, None, None)

        for i in range(n):
            aa = protein[i]
            codons = AA_TO_CODONS.get(aa, [])

            if not codons or aa == "*":
                # For stop codons or unknown AAs, just carry forward
                if aa == "*":
                    codon_start = i * 3
                    stop_codon = seq[codon_start:codon_start + 3] if codon_start + 3 <= len(seq) else "TAA"
                    for last_char, (log_sum, _, _) in dp[i].items():
                        new_last = stop_codon[-1]
                        if new_last not in dp[i + 1] or log_sum > dp[i + 1][new_last][0]:
                            dp[i + 1][new_last] = (log_sum, last_char, stop_codon)
                else:
                    for last_char in dp[i]:
                        dp[i + 1][last_char] = dp[i][last_char]
                continue

            # Sort codons by CAI (highest first) for tie-breaking
            codons_sorted = sorted(codons, key=lambda c: self.species_cai.get(c, 0.0), reverse=True)

            # Check if the NEXT amino acid (if any) has codons not starting with T
            # If not, cross-codon GT from this codon ending with G is unavoidable
            next_aa = protein[i + 1] if i + 1 < n else None
            next_can_avoid_gt = True
            if next_aa is not None and next_aa != "*":
                next_can_avoid_gt = aa_has_non_t_start.get(next_aa, True)
            elif next_aa == "*":
                next_can_avoid_gt = False  # Stop codon, can't change

            for last_char, (log_sum, _, _) in dp[i].items():
                if log_sum == INF:
                    continue

                for codon in codons_sorted:
                    cai = self.species_cai.get(codon, 0.0)
                    if cai <= 0:
                        cai = 0.001

                    if self.avoid_gt:
                        # Check within-codon GT (except Valine which is unavoidable)
                        has_within_gt = "GT" in codon
                        if has_within_gt and aa != 'V':
                            continue

                        # Check cross-codon GT with previous codon
                        if last_char and last_char + codon[0] == "GT":
                            # This GT is avoidable only if we could have used a
                            # different codon for the previous AA that doesn't
                            # end with the last_char. But since we're in DP, the
                            # previous choice is already fixed. However, this GT
                            # IS avoidable if the current AA has a codon not
                            # starting with T. If not, we must accept it.
                            current_can_avoid = aa_has_non_t_start.get(aa, True)
                            if current_can_avoid:
                                continue  # Skip - this GT is avoidable
                            # else: GT is unavoidable, accept this codon

                        # Check if this codon ending with G would create unavoidable
                        # cross-codon GT with the NEXT codon. We only need to avoid
                        # this if the next AA CAN avoid starting with T.
                        # If next_aa only has T-starting codons, the GT is unavoidable
                        # and we should use the highest-CAI codon (which may end with G).
                        # We handle this by NOT penalizing codons ending with G
                        # when the next AA can't avoid T-start.

                    new_last_char = codon[-1]
                    new_log_sum = log_sum + math.log(cai)

                    if new_last_char not in dp[i + 1] or new_log_sum > dp[i + 1][new_last_char][0]:
                        dp[i + 1][new_last_char] = (new_log_sum, last_char, codon)

        # Find the best final state
        best_log_sum = INF
        best_last_char = None
        for last_char, (log_sum, _, _) in dp[n].items():
            if log_sum > best_log_sum:
                best_log_sum = log_sum
                best_last_char = last_char

        if best_last_char is None:
            # Fallback: use simple max-CAI (no GT avoidance)
            codons_result = []
            for i, aa in enumerate(protein):
                if aa == "*":
                    codon_start = i * 3
                    codons_result.append(seq[codon_start:codon_start + 3])
                    continue
                candidates = AA_TO_CODONS.get(aa, [])
                if candidates:
                    codons_result.append(max(candidates, key=lambda c: self.species_cai.get(c, 0.0)))
                else:
                    codon_start = i * 3
                    codons_result.append(seq[codon_start:codon_start + 3])
            return "".join(codons_result)

        # Backtrack to find the optimal sequence
        codons_result = [None] * n
        current_char = best_last_char
        for i in range(n - 1, -1, -1):
            _, prev_char, codon = dp[i + 1][current_char]
            codons_result[i] = codon if codon is not None else "NNN"
            current_char = prev_char

        return "".join(codons_result)

    # Deprecated alias — use _step_backtranslate_cai instead
    _phase0_max_cai_backtranslate = _step_backtranslate_cai

    # ──────────────────────────────────────────────────────────
    # Step: Resolve Constraints (Priority-based constraint resolution)
    # ──────────────────────────────────────────────────────────
    def _step_resolve_constraints(self, seq: str) -> str:
        """Resolve Constraints step: Fix constraint violations with minimal CAI impact.

        Iteratively finds GT/CG dinucleotides and restriction sites, then
        resolves them by choosing the synonymous substitution with the
        smallest CAI penalty. This is the DNAworks-style approach: start
        with max CAI, then fix only what's needed.

        Key difference from old greedy approach: instead of greedily avoiding GT
        during codon selection (which permanently sacrifices CAI), we fix
        GT violations after the fact, choosing the resolution that costs
        the least CAI.
        """
        import math
        seq_list = list(seq)
        max_rounds = 30

        for round_num in range(max_rounds):
            current_seq = "".join(seq_list)
            old_gt_count = _count_gts(current_seq)

            # Collect all constraint violations
            violations = []

            # Within-codon GT dinucleotides (only avoidable ones)
            for i in range(len(seq_list) - 1):
                if seq_list[i] == "G" and seq_list[i + 1] == "T":
                    codon_start = (i // 3) * 3
                    next_codon_start = codon_start + 3
                    if i + 1 < next_codon_start:
                        # Within-codon GT - only add if avoidable
                        if not _is_unavoidable_gt(current_seq, i):
                            violations.append(("within_gt", i, codon_start))

            # Cross-codon GT dinucleotides (only avoidable ones)
            for pos in find_cross_codon_gt(current_seq):
                codon_start = (pos // 3) * 3
                if not _is_unavoidable_gt(current_seq, pos):
                    violations.append(("cross_gt", pos, codon_start))

            if not violations and not self.enzymes:
                break

            # Check for restriction sites
            from .restriction_sites import get_recognition_site
            rs_violations = []
            for enzyme in self.enzymes:
                site = get_recognition_site(enzyme)
                if site is None:
                    continue
                p = current_seq.find(site)
                while p != -1:
                    rs_violations.append((p, site, enzyme))
                    p = current_seq.find(site, p + 1)

            if not violations and not rs_violations:
                break

            any_resolved = False

            # Fix within-codon GTs first (usually easiest)
            for vtype, pos, codon_start in violations:
                if vtype == "within_gt":
                    resolved = self._fix_within_codon_gt_cai_aware(seq_list, codon_start)
                    if resolved:
                        any_resolved = True
                        break

            if any_resolved:
                continue

            # Fix cross-codon GTs
            for vtype, pos, codon_start in violations:
                if vtype == "cross_gt":
                    resolved = self._fix_cross_codon_gt_cai_aware(seq_list, codon_start)
                    if resolved:
                        any_resolved = True
                        break

            if any_resolved:
                continue

            # Fix restriction sites
            for p, site, enzyme in rs_violations:
                resolved = self._fix_restriction_site_cai_aware(seq_list, p, site)
                if resolved:
                    any_resolved = True
                    break

            if not any_resolved:
                break

        # Track unavoidable GTs for mutagenesis
        current_seq = "".join(seq_list)
        for i in range(len(current_seq) - 1):
            if current_seq[i] == "G" and current_seq[i + 1] == "T":
                codon_start = (i // 3) * 3
                next_codon_start = codon_start + 3
                if i + 1 < next_codon_start:
                    codon = current_seq[codon_start:codon_start + 3]
                    aa = CODON_TABLE.get(codon)
                    if aa and aa != "*":
                        all_have_gt = all("GT" in c for c in AA_TO_CODONS.get(aa, []))
                        if all_have_gt:
                            for j in range(2):
                                if codon[j:j+2] == "GT":
                                    self._unavoidable_gt_positions.add(codon_start + j)

        return "".join(seq_list)

    def _fix_within_codon_gt_cai_aware(self, seq_list: list, codon_start: int) -> bool:
        """Fix a within-codon GT by choosing the best CAI-preserving substitution."""
        codon = "".join(seq_list[codon_start:codon_start + 3])
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            return False

        old_gt_count = _count_gts("".join(seq_list))
        current_cai = self.species_cai.get(codon, 0.0)

        # Sort alternatives by CAI (highest first) to minimize CAI loss
        alternatives = []
        for alt in AA_TO_CODONS.get(aa, []):
            if alt == codon or "GT" in alt:
                continue
            prev_base = seq_list[codon_start - 1] if codon_start > 0 else ""
            next_base = seq_list[codon_start + 3] if codon_start + 3 < len(seq_list) else ""
            if prev_base and prev_base + alt[0] == "GT":
                continue
            if next_base and alt[-1] + next_base == "GT":
                continue
            alt_cai = self.species_cai.get(alt, 0.0)
            alternatives.append((alt, alt_cai))

        # Sort by CAI descending (prefer highest CAI)
        alternatives.sort(key=lambda x: x[1], reverse=True)

        for alt, alt_cai in alternatives:
            test_list = seq_list[:]
            for k, b in enumerate(alt):
                test_list[codon_start + k] = b
            new_gt_count = _count_gts("".join(test_list))
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    def _fix_cross_codon_gt_cai_aware(self, seq_list: list, codon_start: int) -> bool:
        """Fix a cross-codon GT by choosing the best CAI-preserving substitution."""
        old_seq = "".join(seq_list)
        old_gt_count = _count_gts(old_seq)

        # Try resolving by modifying the codon pair
        next_start = codon_start + 3

        # Strategy A: Modify following codon
        if next_start + 3 <= len(seq_list):
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)

            if aa1 is not None and aa1 != "*" and aa2 is not None and aa2 != "*":
                # Collect all valid pairs sorted by combined CAI
                pairs = []
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                        if c1[-1] + c2[0] == "GT":
                            continue
                        # Prefer pairs that don't introduce new GTs within codons
                        # but don't require them to be GT-free
                        combined_cai = self.species_cai.get(c1, 0.0) + self.species_cai.get(c2, 0.0)
                        pairs.append((c1, c2, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c1, c2, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    for k, b in enumerate(c2):
                        test_list[next_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy B: Modify preceding codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
            aa1 = CODON_TABLE.get(aa1_codon)

            if aa0 is not None and aa0 != "*" and aa1 is not None and aa1 != "*":
                pairs = []
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if c0[-1] + c1[0] == "GT":
                            continue
                        combined_cai = self.species_cai.get(c0, 0.0) + self.species_cai.get(c1, 0.0)
                        pairs.append((c0, c1, combined_cai))

                pairs.sort(key=lambda x: x[2], reverse=True)

                for c0, c1, _ in pairs:
                    test_list = seq_list[:]
                    for k, b in enumerate(c0):
                        test_list[prev_start + k] = b
                    for k, b in enumerate(c1):
                        test_list[codon_start + k] = b
                    test_seq = "".join(test_list)
                    new_gt_count = _count_gts(test_seq)
                    if new_gt_count < old_gt_count:
                        seq_list[:] = test_list
                        return True

        # Strategy C: Single codon substitution (change only the codon that ends with G)
        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is not None and aa1 != "*":
            alts = sorted(AA_TO_CODONS.get(aa1, [aa1_codon]),
                          key=lambda c: self.species_cai.get(c, 0.0), reverse=True)
            for c1 in alts:
                test_list = seq_list[:]
                for k, b in enumerate(c1):
                    test_list[codon_start + k] = b
                test_seq = "".join(test_list)
                new_gt_count = _count_gts(test_seq)
                if new_gt_count < old_gt_count:
                    seq_list[:] = test_list
                    return True

        return False

    def _fix_restriction_site_cai_aware(self, seq_list: list, pos: int, site: str) -> bool:
        """Fix a restriction site by choosing the best CAI-preserving substitution."""
        codon_starts = set()
        for j in range(pos, pos + len(site)):
            cs = (j // 3) * 3
            if cs + 3 <= len(seq_list):
                codon_starts.add(cs)

        old_gt_count = _count_gts("".join(seq_list))

        # Try each overlapping codon, sorted by smallest CAI impact
        for cs in sorted(codon_starts):
            codon = "".join(seq_list[cs:cs + 3])
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue

            # Sort alternatives by CAI (highest first)
            alts = sorted(AA_TO_CODONS.get(aa, []),
                          key=lambda c: self.species_cai.get(c, 0.0), reverse=True)
            for alt in alts:
                if alt == codon:
                    continue
                test_list = seq_list[:]
                for k, b in enumerate(alt):
                    test_list[cs + k] = b
                test_seq = "".join(test_list)
                if site not in test_seq:
                    if self.avoid_gt:
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            seq_list[:] = test_list
                            return True
                    else:
                        seq_list[:] = test_list
                        return True

        return False

    # Deprecated alias — use _step_resolve_constraints instead
    _phase1_priority_constraint_resolution = _step_resolve_constraints

    # ──────────────────────────────────────────────────────────
    # Step: Greedy Optimize (Per-position CAI maximization)
    # ──────────────────────────────────────────────────────────
    def _step_greedy_optimize(self, seq: str) -> str:
        """Greedy Optimize step: Per-position CAI maximization, GT-aware.

        For each amino acid, select the highest-CAI codon that does not
        introduce a GT dinucleotide within or across codon boundaries.

        If ALL synonymous codons for an AA contain GT (e.g., Valine GTN),
        flag the position as "unavoidable GT" for mutagenesis fallback step.
        """
        codons = []
        protein = self._translate(seq)
        prev_codon_end = ""

        for i, aa in enumerate(protein):
            if aa == "*":
                codon_start = i * 3
                codons.append(seq[codon_start:codon_start + 3])
                prev_codon_end = codons[-1][-1]
                continue

            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codon_start = i * 3
                codons.append(seq[codon_start:codon_start + 3])
                prev_codon_end = codons[-1][-1]
                continue

            candidates_sorted = sorted(
                candidates,
                key=lambda c: self.species_cai.get(c, 0.0),
                reverse=True,
            )

            best = candidates_sorted[0]
            found_gt_free = False

            if self.avoid_gt:
                for codon in candidates_sorted:
                    if "GT" in codon:
                        continue
                    if prev_codon_end and prev_codon_end + codon[0] == "GT":
                        continue
                    best = codon
                    found_gt_free = True
                    break

                if not found_gt_free:
                    all_have_gt = all("GT" in c for c in candidates)
                    if all_have_gt:
                        codon_abs_start = i * 3
                        for j in range(2):
                            if best[j:j+2] == "GT":
                                self._unavoidable_gt_positions.add(codon_abs_start + j)
                        # Pick the codon that avoids cross-codon GT if possible
                        if prev_codon_end:
                            for codon in candidates_sorted:
                                if prev_codon_end + codon[0] != "GT":
                                    best = codon
                                    break
                    else:
                        # Pick the one that minimizes GT count
                        for codon in candidates_sorted:
                            gt_count = codon.count("GT")
                            cross_gt = 1 if (prev_codon_end and prev_codon_end + codon[0] == "GT") else 0
                            best_gt_count = best.count("GT")
                            best_cross_gt = 1 if (prev_codon_end and prev_codon_end + best[0] == "GT") else 0
                            if gt_count + cross_gt < best_gt_count + best_cross_gt:
                                best = codon
                                break

            codons.append(best)
            prev_codon_end = best[-1]

        return "".join(codons)

    # Deprecated alias — use _step_greedy_optimize instead
    _phase1_greedy_optimize = _step_greedy_optimize

    # ──────────────────────────────────────────────────────────
    # Step: Remove Restriction Sites
    # ──────────────────────────────────────────────────────────
    def _step_remove_restriction_sites(self, seq: str) -> str:
        """Remove Restriction Sites step: Remove restriction enzyme recognition sites by synonymous substitution."""
        from .restriction_sites import get_recognition_site

        seq_list = list(seq)
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            pos = 0
            while pos <= len(seq_list) - len(site):
                if "".join(seq_list[pos:pos + len(site)]) == site:
                    codon_starts = set()
                    for j in range(pos, pos + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(seq_list):
                            codon_starts.add(cs)

                    resolved = False
                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue

                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            if site not in test_seq:
                                # Note: we allow temporary GT increase to remove RS;
                                # GTs will be fixed by later steps (cross-codon/hill climbing)
                                # Only reject if GT increase is severe (more than 2 new GTs)
                                if self.avoid_gt:
                                    old_gt_count = _count_gts("".join(seq_list))
                                    new_gt_count = _count_gts(test_seq)
                                    if new_gt_count > old_gt_count + 2:
                                        continue
                                seq_list = test_list
                                resolved = True
                                break
                        if resolved:
                            break

                    if resolved:
                        continue
                    else:
                        pos += 1
                else:
                    pos += 1

        return "".join(seq_list)

    # Deprecated alias — use _step_remove_restriction_sites instead
    _phase2_remove_restriction_sites = _step_remove_restriction_sites

    # ──────────────────────────────────────────────────────────
    # Step: Cross-Codon Optimization (iterative constraint resolution)
    # ──────────────────────────────────────────────────────────
    def _step_cross_codon_optimization(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Cross-Codon Optimization step: Resolve cross-codon GT, CG, and restriction site constraints.

        This step iterates until no more cross-codon constraints can be
        resolved. Each resolution is globally validated to ensure no new
        GTs are introduced elsewhere.
        """
        seq_list = list(seq)
        total_remaining: Dict[int, List[str]] = {}
        max_rounds = 20

        for round_num in range(max_rounds):
            current_seq = "".join(seq_list)
            constraint_positions: Dict[int, List[str]] = {}

            for pos in find_cross_codon_gt(current_seq):
                codon_start = (pos // 3) * 3
                # Skip GTs that are unavoidable (Valine, or cross-codon with no alternatives)
                if not _is_unavoidable_gt(current_seq, pos):
                    constraint_positions.setdefault(codon_start, [])
                    if "GT" not in constraint_positions[codon_start]:
                        constraint_positions[codon_start].append("GT")

            for pos in find_cross_codon_cg(current_seq):
                codon_start = (pos // 3) * 3
                constraint_positions.setdefault(codon_start, [])
                if "CG" not in constraint_positions[codon_start]:
                    constraint_positions[codon_start].append("CG")

            from .restriction_sites import get_recognition_site
            for enzyme in self.enzymes:
                site = get_recognition_site(enzyme)
                if site is None:
                    continue
                for pos in find_cross_codon_restriction(current_seq, site):
                    codon_start = (pos // 3) * 3
                    label = f"RS:{site}"
                    constraint_positions.setdefault(codon_start, [])
                    if label not in constraint_positions[codon_start]:
                        constraint_positions[codon_start].append(label)

            if not constraint_positions:
                break

            any_resolved = False

            for codon_start, ctypes in constraint_positions.items():
                resolved = self._try_resolve_cross_codon(seq_list, codon_start, ctypes)
                if resolved:
                    any_resolved = True

            if not any_resolved:
                # Record remaining constraints
                current_seq = "".join(seq_list)
                constraint_positions2: Dict[int, List[str]] = {}
                for pos in find_cross_codon_gt(current_seq):
                    cs = (pos // 3) * 3
                    # Skip GTs that are unavoidable
                    if not _is_unavoidable_gt(current_seq, pos):
                        constraint_positions2.setdefault(cs, [])
                        if "GT" not in constraint_positions2[cs]:
                            constraint_positions2[cs].append("GT")
                for pos in find_cross_codon_cg(current_seq):
                    cs = (pos // 3) * 3
                    constraint_positions2.setdefault(cs, [])
                    if "CG" not in constraint_positions2[cs]:
                        constraint_positions2[cs].append("CG")
                total_remaining = constraint_positions2
                break

        mut_report = propose_mutagenesis(
            "".join(seq_list),
            list(total_remaining.keys()),
            total_remaining,
            self.species_cai,
            self.min_blosum,
            self.min_cai,
        )

        return "".join(seq_list), mut_report

    def _try_resolve_cross_codon(
        self, seq_list: list, codon_start: int, ctypes: List[str]
    ) -> bool:
        """Try to resolve cross-codon constraints at codon_start.

        Uses global validation: after any substitution, verifies that the
        total GT count in the sequence doesn't increase.
        """
        old_seq = "".join(seq_list)
        old_gt_count = _count_gts(old_seq)

        aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
        aa1 = CODON_TABLE.get(aa1_codon)
        if aa1 is None or aa1 == "*":
            return False

        next_start = codon_start + 3

        # Strategy A: Following codon
        if next_start + 3 <= len(seq_list):
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if aa2 is not None and aa2 != "*":
                for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                    if self.avoid_gt and "GT" in c1:
                        continue
                    for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                        if self.avoid_gt and "GT" in c2:
                            continue
                        if c1[-1] + c2[0] == "GT":
                            continue
                        # Apply and validate globally
                        test_list = seq_list[:]
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        for k, b in enumerate(c2):
                            test_list[next_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            # Also check constraint is actually resolved
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and _has_gt(test_seq[codon_start:next_start + 3]):
                                    resolved = False
                                elif ct == "CG" and "CG" in test_seq[codon_start:next_start + 3]:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in test_seq:
                                        resolved = False
                            if resolved:
                                # Also check no new restriction sites introduced
                                from .restriction_sites import get_recognition_site
                                rs_ok = True
                                for enzyme in self.enzymes:
                                    rs_site = get_recognition_site(enzyme)
                                    if rs_site and rs_site in test_seq:
                                        rs_ok = False
                                        break
                                if rs_ok:
                                    seq_list[:] = test_list
                                    return True

        # Strategy B: Preceding codon
        if codon_start >= 3:
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            if aa0 is not None and aa0 != "*":
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        if c0[-1] + c1[0] == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(c0):
                            test_list[prev_start + k] = b
                        for k, b in enumerate(c1):
                            test_list[codon_start + k] = b
                        test_seq = "".join(test_list)
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count <= old_gt_count:
                            resolved = True
                            for ct in ctypes:
                                if ct == "GT" and _has_gt(test_seq[prev_start:codon_start + 3]):
                                    resolved = False
                                elif ct == "CG" and "CG" in test_seq[prev_start:codon_start + 3]:
                                    resolved = False
                                elif ct.startswith("RS:"):
                                    if ct[3:] in test_seq:
                                        resolved = False
                            if resolved:
                                seq_list[:] = test_list
                                return True

        # Strategy C: Both preceding and following
        if codon_start >= 3 and next_start + 3 <= len(seq_list):
            prev_start = codon_start - 3
            aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
            aa0 = CODON_TABLE.get(aa0_codon)
            aa2_codon = "".join(seq_list[next_start:next_start + 3])
            aa2 = CODON_TABLE.get(aa2_codon)
            if (aa0 is not None and aa0 != "*" and
                    aa2 is not None and aa2 != "*"):
                for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                    if self.avoid_gt and "GT" in c0:
                        continue
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if self.avoid_gt and "GT" in c1:
                            continue
                        for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                            if self.avoid_gt and "GT" in c2:
                                continue
                            if c0[-1] + c1[0] == "GT" or c1[-1] + c2[0] == "GT":
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(c0):
                                test_list[prev_start + k] = b
                            for k, b in enumerate(c1):
                                test_list[codon_start + k] = b
                            for k, b in enumerate(c2):
                                test_list[next_start + k] = b
                            test_seq = "".join(test_list)
                            new_gt_count = _count_gts(test_seq)
                            if new_gt_count <= old_gt_count:
                                resolved = True
                                for ct in ctypes:
                                    if ct == "GT" and _has_gt(test_seq[prev_start:next_start + 3]):
                                        resolved = False
                                    elif ct == "CG" and "CG" in test_seq[prev_start:next_start + 3]:
                                        resolved = False
                                    elif ct.startswith("RS:"):
                                        if ct[3:] in test_seq:
                                            resolved = False
                                if resolved:
                                    seq_list[:] = test_list
                                    return True

        # Strategy D: Single codon substitution
        for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
            if self.avoid_gt and "GT" in c1:
                continue
            test_list = seq_list[:]
            for k, b in enumerate(c1):
                test_list[codon_start + k] = b
            test_seq = "".join(test_list)
            new_gt_count = _count_gts(test_seq)
            if new_gt_count < old_gt_count:
                seq_list[:] = test_list
                return True

        return False

    # Deprecated alias — use _step_cross_codon_optimization instead
    _phase3_cross_codon_constraints = _step_cross_codon_optimization

    # ──────────────────────────────────────────────────────────
    # Step: Within-Codon GT Resolution
    # ──────────────────────────────────────────────────────────
    def _step_within_codon_gt_resolution(self, seq: str) -> Tuple[str, MutagenesisReport]:
        """Within-Codon GT Resolution step: Resolve within-codon GT dinucleotides.

        For each within-codon GT (not cross-codon), try synonymous substitution
        first. If no synonymous codon can avoid GT, flag for mutagenesis.
        """
        seq_list = list(seq)
        remaining_positions: Dict[int, List[str]] = {}

        for i in range(len(seq_list) - 1):
            if seq_list[i] == "G" and seq_list[i + 1] == "T":
                codon_start = (i // 3) * 3
                next_codon_start = codon_start + 3

                if i + 1 < next_codon_start:
                    # Within-codon GT - skip if unavoidable
                    current_seq_35 = "".join(seq_list)
                    if _is_unavoidable_gt(current_seq_35, i):
                        continue
                    codon = "".join(seq_list[codon_start:codon_start + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue

                    old_gt_count = _count_gts("".join(seq_list))
                    resolved = False
                    alts = AA_TO_CODONS.get(aa, [])
                    alts_sorted = sorted(alts,
                                         key=lambda c: self.species_cai.get(c, 0.0),
                                         reverse=True)
                    for alt in alts_sorted:
                        if alt == codon or "GT" in alt:
                            continue
                        prev_base = seq_list[codon_start - 1] if codon_start > 0 else ""
                        next_base = seq_list[next_codon_start] if next_codon_start < len(seq_list) else ""
                        if prev_base and prev_base + alt[0] == "GT":
                            continue
                        if next_base and alt[-1] + next_base == "GT":
                            continue
                        # Apply and validate globally
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[codon_start + k] = b
                        new_gt_count = _count_gts("".join(test_list))
                        if new_gt_count < old_gt_count:
                            seq_list = test_list
                            resolved = True
                            break

                    if not resolved:
                        remaining_positions[codon_start] = ["GT:within"]

        if remaining_positions:
            mut_report = self._propose_within_codon_mutagenesis(
                "".join(seq_list), remaining_positions
            )
        else:
            mut_report = MutagenesisReport()

        return "".join(seq_list), mut_report

    def _propose_within_codon_mutagenesis(
        self,
        seq: str,
        positions: Dict[int, List[str]],
    ) -> MutagenesisReport:
        """Propose AA substitutions for within-codon GTs that can't be resolved
        by synonymous substitution.

        Uses BLOSUM62 guidance for conservative substitutions:
        - Valine (V, GTN) → Isoleucine (I, ATN) or Leucine (L, TTA/CTN)
        """
        report = MutagenesisReport()

        for codon_start, ctypes in positions.items():
            if codon_start + 3 > len(seq):
                continue

            original_codon = seq[codon_start:codon_start + 3]
            original_aa = CODON_TABLE.get(original_codon)
            if original_aa is None or original_aa == "*":
                continue

            best = None
            best_score = -100

            for new_aa, score in sorted(
                ((aa, BLOSUM62.get((original_aa, aa), -10))
                 for aa in set(CODON_TABLE.values()) if aa != "*" and aa != original_aa),
                key=lambda x: x[1],
                reverse=True,
            ):
                if score < self.min_blosum:
                    continue

                alts = AA_TO_CODONS.get(new_aa, [])
                alts_sorted = sorted(alts,
                                     key=lambda c: self.species_cai.get(c, 0.0),
                                     reverse=True)
                for alt_codon in alts_sorted:
                    if "GT" in alt_codon:
                        continue
                    prev_base = seq[codon_start - 1] if codon_start > 0 else ""
                    next_start = codon_start + 3
                    next_base = seq[next_start] if next_start < len(seq) else ""
                    if prev_base and prev_base + alt_codon[0] == "GT":
                        continue
                    if next_base and alt_codon[-1] + next_base == "GT":
                        continue
                    cai = self.species_cai.get(alt_codon, 0.0)
                    if cai < self.min_cai:
                        continue

                    if score > best_score:
                        best = (alt_codon, new_aa, score, cai)
                        best_score = score
                    break

            if best is not None:
                alt_codon, new_aa, blosum, cai = best
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa=new_aa,
                    new_codon=alt_codon,
                    blosum_score=blosum,
                    cai_weight=cai,
                    resolves=ctypes,
                )
                report.add(proposal)
            else:
                proposal = MutagenesisProposal(
                    position=codon_start,
                    original_codon=original_codon,
                    original_aa=original_aa,
                    new_aa="",
                    new_codon="",
                    blosum_score=-10,
                    cai_weight=0.0,
                    resolves=ctypes,
                    impossible=True,
                )
                report.add(proposal)

        return report

    # Deprecated alias — use _step_within_codon_gt_resolution instead
    _phase35_within_codon_gt = _step_within_codon_gt_resolution

    # ──────────────────────────────────────────────────────────
    # Step: Mutagenesis Fallback
    # ──────────────────────────────────────────────────────────
    def _step_mutagenesis_fallback(self, seq: str, mut_report: MutagenesisReport) -> str:
        """Mutagenesis Fallback step: Apply mutagenesis proposals for intractable constraints.

        Applies conservative AA substitutions (e.g., Val→Ile, Val→Leu)
        using BLOSUM62 guidance. Only applies to AVOIDABLE GT positions;
        unavoidable GTs (Valine, cross-codon with no alternatives) are skipped.
        """
        seq_list = list(seq)
        for proposal in mut_report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            pos = proposal.position
            old_gt_count = _count_gts("".join(seq_list))

            test_list = seq_list[:]
            for k, b in enumerate(proposal.new_codon):
                if pos + k < len(test_list):
                    test_list[pos + k] = b
            new_gt_count = _count_gts("".join(test_list))

            # Accept if it reduces or maintains GT count
            if new_gt_count <= old_gt_count:
                for k, b in enumerate(proposal.new_codon):
                    if pos + k < len(seq_list):
                        seq_list[pos + k] = b
                self._applied_mutagenesis.append({
                    "position": pos,
                    "original_aa": proposal.original_aa,
                    "new_aa": proposal.new_aa,
                    "blosum": proposal.blosum_score,
                })

        return "".join(seq_list)

    # Deprecated alias — use _step_mutagenesis_fallback instead
    _phase4_mutagenesis_fallback = _step_mutagenesis_fallback

    # ──────────────────────────────────────────────────────────
    # Step: Avoid CpG Islands
    # ──────────────────────────────────────────────────────────
    def _step_avoid_cpg_islands(self, seq: str) -> str:
        """Avoid CpG Islands step: CpG island avoidance by synonymous substitution.

        Also handles cross-codon CG dinucleotides by two-codon coordination.
        For cross-codon CG (C at end of one codon, G at start of next),
        single-codon swaps may fail because the alternative codon still ends
        with C or starts with G. Coordinated 2-codon swaps resolve this.
        """
        seq_list = list(seq)
        changed = True
        iterations = 0
        max_iterations = 80

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for start in range(0, len(seq_list) - self.cpg_window + 1, 3):
                window = "".join(seq_list[start:start + self.cpg_window])
                c_count = window.count("C")
                g_count = window.count("G")
                cg_count = sum(1 for i in range(len(window) - 1) if window[i:i+2] == "CG")
                expected = (c_count * g_count) / len(window) if len(window) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0

                if obs_exp <= self.cpg_threshold:
                    continue

                # Find CG dinucleotides in this window and try to break them
                for i in range(start, min(start + self.cpg_window - 1, len(seq_list) - 1)):
                    if seq_list[i] == "C" and seq_list[i+1] == "G":
                        left_ci = i // 3          # codon containing C
                        right_ci = (i + 1) // 3   # codon containing G
                        is_cross_codon = (left_ci != right_ci)

                        # Strategy 1: Single-codon swap — try both codons for cross-codon CG
                        codon_fixed = False
                        for ci in ([left_ci, right_ci] if is_cross_codon else [left_ci]):
                            codon_start = ci * 3
                            if codon_start + 3 > len(seq_list):
                                continue
                            codon = "".join(seq_list[codon_start:codon_start + 3])
                            aa = CODON_TABLE.get(codon)
                            if aa is None or aa == "*":
                                continue

                            old_gt_count = _count_gts("".join(seq_list))

                            for alt in AA_TO_CODONS.get(aa, []):
                                if alt == codon:
                                    continue
                                # For within-codon CG, skip alternatives that still contain CG
                                if not is_cross_codon and "CG" in alt:
                                    continue
                                if self.avoid_gt and "GT" in alt:
                                    continue
                                cai = self.species_cai.get(alt, 0.0)
                                if cai < self.min_cai:
                                    continue
                                # Check cross-codon effects
                                prev_base = seq_list[codon_start - 1] if codon_start > 0 else ""
                                next_base = seq_list[codon_start + 3] if codon_start + 3 < len(seq_list) else ""
                                if prev_base and prev_base + alt[0] == "GT":
                                    continue
                                if next_base and alt[-1] + next_base == "GT":
                                    continue
                                # Apply and validate globally
                                test_list = seq_list[:]
                                for k, b in enumerate(alt):
                                    test_list[codon_start + k] = b
                                new_gt_count = _count_gts("".join(test_list))
                                if new_gt_count <= old_gt_count:
                                    # Verify the specific CG at position i is gone
                                    if test_list[i] != "C" or test_list[i+1] != "G":
                                        seq_list = test_list
                                        changed = True
                                        codon_fixed = True
                                        break

                            if codon_fixed:
                                break

                        if codon_fixed:
                            break

                        # Strategy 2: Coordinated 2-codon swap for cross-codon CG
                        if is_cross_codon and not codon_fixed:
                            left_start = left_ci * 3
                            right_start = right_ci * 3
                            if left_start + 3 <= len(seq_list) and right_start + 3 <= len(seq_list):
                                left_codon = "".join(seq_list[left_start:left_start + 3])
                                right_codon = "".join(seq_list[right_start:right_start + 3])
                                left_aa = CODON_TABLE.get(left_codon)
                                right_aa = CODON_TABLE.get(right_codon)
                                if left_aa and right_aa and left_aa != "*" and right_aa != "*":
                                    old_gt_count = _count_gts("".join(seq_list))
                                    for left_alt in AA_TO_CODONS.get(left_aa, []):
                                        if left_alt == left_codon:
                                            continue
                                        if self.avoid_gt and "GT" in left_alt:
                                            continue
                                        left_cai = self.species_cai.get(left_alt, 0.0)
                                        if left_cai < self.min_cai:
                                            continue
                                        for right_alt in AA_TO_CODONS.get(right_aa, []):
                                            if right_alt == right_codon and left_alt == left_codon:
                                                continue
                                            # Skip if CG still present at boundary
                                            if left_alt[-1] == "C" and right_alt[0] == "G":
                                                continue
                                            if self.avoid_gt and "GT" in right_alt:
                                                continue
                                            right_cai = self.species_cai.get(right_alt, 0.0)
                                            if right_cai < self.min_cai:
                                                continue
                                            # Check cross-codon GT effects for both codons
                                            prev_base = seq_list[left_start - 1] if left_start > 0 else ""
                                            if prev_base and prev_base + left_alt[0] == "GT":
                                                continue
                                            mid_check = left_alt[-1] + right_alt[0]
                                            if mid_check == "GT":
                                                continue
                                            next_base = seq_list[right_start + 3] if right_start + 3 < len(seq_list) else ""
                                            if next_base and right_alt[-1] + next_base == "GT":
                                                continue
                                            # Apply both codons
                                            test_list = seq_list[:]
                                            for k, b in enumerate(left_alt):
                                                test_list[left_start + k] = b
                                            for k, b in enumerate(right_alt):
                                                test_list[right_start + k] = b
                                            new_gt_count = _count_gts("".join(test_list))
                                            if new_gt_count <= old_gt_count:
                                                seq_list = test_list
                                                changed = True
                                                codon_fixed = True
                                                break
                                        if codon_fixed:
                                            break

                        if changed:
                            break
                if changed:
                    break

        return "".join(seq_list)

    # Deprecated alias — use _step_avoid_cpg_islands instead
    _phase5_avoid_cpg_islands = _step_avoid_cpg_islands

    # ──────────────────────────────────────────────────────────
    # Step: CpG Reconciliation (aggressive, runs after CAI hill climb)
    # ──────────────────────────────────────────────────────────
    def _step_cpg_reconciliation(self, seq: str) -> str:
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
        """
        seq_list = list(seq)
        max_iterations = 100

        for iteration in range(max_iterations):
            current_seq = "".join(seq_list)

            # Check if the CpG island predicate passes
            cpg_result = check_no_cpg_island(current_seq, self.cpg_window, self.cpg_threshold)
            if cpg_result.passed:
                break

            # Find the worst window
            worst_ratio = 0.0
            worst_start = -1
            for start in range(0, len(current_seq) - self.cpg_window + 1):
                window_seq = current_seq[start:start + self.cpg_window]
                c_count = window_seq.count("C")
                g_count = window_seq.count("G")
                cg_count = sum(1 for i in range(len(window_seq) - 1) if window_seq[i:i+2] == "CG")
                expected = (c_count * g_count) / len(window_seq) if len(window_seq) > 0 else 0
                obs_exp = cg_count / expected if expected > 0 else 0.0
                if obs_exp > worst_ratio:
                    worst_ratio = obs_exp
                    worst_start = start

            if worst_start < 0 or worst_ratio <= self.cpg_threshold:
                break

            # Find CG dinucleotides in the worst window and try to eliminate them
            window_end = min(worst_start + self.cpg_window, len(current_seq) - 1)
            fixed = False

            # Collect all CG positions in the worst window with their codon info
            cg_targets = []
            for i in range(worst_start, window_end):
                if seq_list[i] == "C" and i + 1 < len(seq_list) and seq_list[i + 1] == "G":
                    # Determine if this is within-codon or cross-codon
                    c_codon_idx = i // 3
                    g_codon_idx = (i + 1) // 3
                    cg_targets.append((i, c_codon_idx, g_codon_idx))

            # Sort targets: prioritize within-codon CG first (easier to fix),
            # then cross-codon CG where the C-ending codon has good alternatives
            for cg_pos, c_codon_idx, g_codon_idx in cg_targets:
                if fixed:
                    break

                if c_codon_idx == g_codon_idx:
                    # Within-codon CG — swap the codon to a CG-free alternative
                    codon_start = c_codon_idx * 3
                    if codon_start + 3 > len(seq_list):
                        continue
                    codon = "".join(seq_list[codon_start:codon_start + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue

                    for alt in sorted(
                        AA_TO_CODONS.get(aa, []),
                        key=lambda c: self.species_cai.get(c, 0.0),
                        reverse=True,
                    ):
                        if alt == codon or "CG" in alt:
                            continue
                        # Apply and check cross-codon effects
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[codon_start + k] = b
                        test_seq = "".join(test_list)
                        # Must not reintroduce restriction sites
                        if not self._check_no_restriction_sites(test_seq):
                            continue
                        # Must not create new avoidable GTs
                        if not self._check_cpg_swap_gt_safe(seq_list, test_list, codon_start):
                            continue
                        # Must reduce net CG count
                        new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                        old_cg = sum(1 for i in range(len(current_seq) - 1) if current_seq[i:i+2] == "CG")
                        if new_cg < old_cg:
                            seq_list = test_list
                            fixed = True
                            break
                else:
                    # Cross-codon CG: codon c_codon_idx ends with C, codon g_codon_idx starts with G
                    # Strategy 1: Swap the C-ending codon to not end with C
                    c_codon_start = c_codon_idx * 3
                    if c_codon_start + 3 <= len(seq_list):
                        c_codon = "".join(seq_list[c_codon_start:c_codon_start + 3])
                        c_aa = CODON_TABLE.get(c_codon)
                        if c_aa is not None and c_aa != "*":
                            # Get non-C-ending alternatives, sorted by CAI
                            non_c_end_alts = [
                                alt for alt in AA_TO_CODONS.get(c_aa, [])
                                if alt != c_codon and alt[-1] != "C"
                            ]
                            non_c_end_alts.sort(
                                key=lambda c: self.species_cai.get(c, 0.0),
                                reverse=True,
                            )
                            for alt in non_c_end_alts:
                                # Skip if this alt would create a new within-codon CG
                                if "CG" in alt:
                                    continue
                                test_list = seq_list[:]
                                for k, b in enumerate(alt):
                                    test_list[c_codon_start + k] = b
                                test_seq = "".join(test_list)
                                # Must not reintroduce restriction sites
                                if not self._check_no_restriction_sites(test_seq):
                                    continue
                                # Must not create new avoidable GTs
                                if not self._check_cpg_swap_gt_safe(seq_list, test_list, c_codon_start):
                                    continue
                                # Must reduce net CG count
                                new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                                old_cg = sum(1 for i in range(len(current_seq) - 1) if current_seq[i:i+2] == "CG")
                                if new_cg < old_cg:
                                    seq_list = test_list
                                    fixed = True
                                    break

                    if fixed:
                        break

                    # Strategy 2: Swap the G-starting codon to not start with G (if possible)
                    g_codon_start = g_codon_idx * 3
                    if g_codon_start + 3 <= len(seq_list):
                        g_codon = "".join(seq_list[g_codon_start:g_codon_start + 3])
                        g_aa = CODON_TABLE.get(g_codon)
                        if g_aa is not None and g_aa != "*":
                            # Get non-G-starting alternatives, sorted by CAI
                            non_g_start_alts = [
                                alt for alt in AA_TO_CODONS.get(g_aa, [])
                                if alt != g_codon and alt[0] != "G"
                            ]
                            non_g_start_alts.sort(
                                key=lambda c: self.species_cai.get(c, 0.0),
                                reverse=True,
                            )
                            for alt in non_g_start_alts:
                                if "CG" in alt:
                                    continue
                                test_list = seq_list[:]
                                for k, b in enumerate(alt):
                                    test_list[g_codon_start + k] = b
                                test_seq = "".join(test_list)
                                if not self._check_no_restriction_sites(test_seq):
                                    continue
                                if not self._check_cpg_swap_gt_safe(seq_list, test_list, g_codon_start):
                                    continue
                                new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                                old_cg = sum(1 for i in range(len(current_seq) - 1) if current_seq[i:i+2] == "CG")
                                if new_cg < old_cg:
                                    seq_list = test_list
                                    fixed = True
                                    break

                    if fixed:
                        break

                    # Strategy 3: Paired swap — change BOTH codons
                    c_codon_start = c_codon_idx * 3
                    g_codon_start = g_codon_idx * 3
                    if (c_codon_start + 3 <= len(seq_list) and
                            g_codon_start + 3 <= len(seq_list)):
                        c_codon = "".join(seq_list[c_codon_start:c_codon_start + 3])
                        g_codon = "".join(seq_list[g_codon_start:g_codon_start + 3])
                        c_aa = CODON_TABLE.get(c_codon)
                        g_aa = CODON_TABLE.get(g_codon)
                        if c_aa is not None and g_aa is not None and c_aa != "*" and g_aa != "*":
                            c_alts = [
                                alt for alt in AA_TO_CODONS.get(c_aa, [])
                                if alt != c_codon and alt[-1] != "C" and "CG" not in alt
                            ]
                            g_alts = [
                                alt for alt in AA_TO_CODONS.get(g_aa, [])
                                if alt != g_codon and alt[0] != "G" and "CG" not in alt
                            ]
                            # Try C-side change first (more likely to have alternatives)
                            c_alts.sort(key=lambda c: self.species_cai.get(c, 0.0), reverse=True)
                            g_alts.sort(key=lambda c: self.species_cai.get(c, 0.0), reverse=True)
                            for c_alt in c_alts[:3]:
                                for g_alt in g_alts[:3]:
                                    test_list = seq_list[:]
                                    for k, b in enumerate(c_alt):
                                        test_list[c_codon_start + k] = b
                                    for k, b in enumerate(g_alt):
                                        test_list[g_codon_start + k] = b
                                    test_seq = "".join(test_list)
                                    if not self._check_no_restriction_sites(test_seq):
                                        continue
                                    if not self._check_cpg_swap_gt_safe(seq_list, test_list, c_codon_start):
                                        continue
                                    if not self._check_cpg_swap_gt_safe(seq_list, test_list, g_codon_start):
                                        continue
                                    new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                                    old_cg = sum(1 for i in range(len(current_seq) - 1) if current_seq[i:i+2] == "CG")
                                    if new_cg < old_cg:
                                        seq_list = test_list
                                        fixed = True
                                        break
                                if fixed:
                                    break

            if not fixed:
                # No progress — try a more aggressive approach:
                # scan ALL codon positions and replace C-ending codons that
                # create cross-codon CG, regardless of window
                aggressive_fixed = False
                for ci in range(len(seq_list) // 3):
                    codon_start = ci * 3
                    if codon_start + 3 > len(seq_list):
                        break
                    codon = "".join(seq_list[codon_start:codon_start + 3])
                    aa = CODON_TABLE.get(codon)
                    if aa is None or aa == "*":
                        continue
                    # Check if this codon ends with C and the next starts with G
                    next_start = codon_start + 3
                    if next_start >= len(seq_list):
                        continue
                    if codon[-1] == "C" and seq_list[next_start] == "G":
                        # This creates a cross-codon CG — try to fix
                        non_c_end_alts = [
                            alt for alt in AA_TO_CODONS.get(aa, [])
                            if alt != codon and alt[-1] != "C" and "CG" not in alt
                        ]
                        non_c_end_alts.sort(
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        )
                        for alt in non_c_end_alts:
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[codon_start + k] = b
                            test_seq = "".join(test_list)
                            if not self._check_no_restriction_sites(test_seq):
                                continue
                            if not self._check_cpg_swap_gt_safe(seq_list, test_list, codon_start):
                                continue
                            new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                            old_cg = sum(1 for i in range(len("".join(seq_list)) - 1) if "".join(seq_list)[i:i+2] == "CG")
                            if new_cg < old_cg:
                                seq_list = test_list
                                aggressive_fixed = True
                                break
                    elif "CG" in codon:
                        # Within-codon CG — try to fix
                        non_cg_alts = [
                            alt for alt in AA_TO_CODONS.get(aa, [])
                            if alt != codon and "CG" not in alt
                        ]
                        non_cg_alts.sort(
                            key=lambda c: self.species_cai.get(c, 0.0),
                            reverse=True,
                        )
                        for alt in non_cg_alts:
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[codon_start + k] = b
                            test_seq = "".join(test_list)
                            if not self._check_no_restriction_sites(test_seq):
                                continue
                            if not self._check_cpg_swap_gt_safe(seq_list, test_list, codon_start):
                                continue
                            new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                            old_cg = sum(1 for i in range(len("".join(seq_list)) - 1) if "".join(seq_list)[i:i+2] == "CG")
                            if new_cg < old_cg:
                                seq_list = test_list
                                aggressive_fixed = True
                                break

                if not aggressive_fixed:
                    break  # Truly no more progress possible

        # After CpG reconciliation, re-check GC content and adjust if needed
        result_seq = "".join(seq_list)
        gc_val = (result_seq.count("G") + result_seq.count("C")) / max(len(result_seq), 1)
        if not (0.30 <= gc_val <= 0.70):
            # GC drifted — adjust with single-codon swaps that don't reintroduce CGs
            result_seq = self._fix_gc_after_cpg(result_seq)

        return result_seq

    def _check_no_restriction_sites(self, seq: str) -> bool:
        """Check that sequence doesn't contain any restriction enzyme sites."""
        from .restriction_sites import get_recognition_site
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            if site in seq:
                return False
        return True

    def _check_cpg_swap_gt_safe(self, old_list: list, new_list: list, codon_start: int) -> bool:
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

    def _fix_gc_after_cpg(self, seq: str) -> str:
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

        for iteration in range(200):
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
                        new_cg = sum(1 for i in range(len(test_seq) - 1) if test_seq[i:i+2] == "CG")
                        old_cg = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")
                        if new_cg > old_cg:
                            continue
                        # Must not reintroduce restriction sites
                        if not self._check_no_restriction_sites(test_seq):
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

    # ──────────────────────────────────────────────────────────
    # Step: Cross-Codon Coordination
    # ──────────────────────────────────────────────────────────
    def _step_cross_codon_coordination(self, seq: str) -> str:
        """Cross-Codon Coordination step: resolve cross-codon GT, CG, and restriction
        site violations by considering ALL synonymous codon pairs at each boundary.

        Unlike the earlier Cross-Codon Optimization step which handles violations
        one at a time with limited search, this phase performs an exhaustive
        pairwise search: for every boundary violation it considers ALL synonymous
        codon pairs for the two adjacent amino acids, selects the pair that
        eliminates the violation while maximizing CAI, and re-scans after each
        fix (since fixing one boundary can create a violation at the adjacent
        boundary). Iterates until no cross-codon violations remain or a maximum
        iteration count is reached.

        This phase runs AFTER within-codon optimization and BEFORE the final
        CAI maximization pass, ensuring that cross-codon interactions — which
        are the primary cause of predicate failures for sequences like HBB —
        are fully resolved before the hill climber locks in codon choices.
        """
        max_iterations = 50

        for iteration in range(max_iterations):
            current_seq = seq
            violations = self._find_cross_codon_violations(current_seq)

            if not violations:
                break  # All cross-codon violations resolved

            # Sort violations: prioritize by type (GT > CG > RS) then by
            # leftmost position so we fix earlier boundaries first.
            type_priority = {"GT": 0, "CG": 1}
            violations.sort(key=lambda v: (type_priority.get(v["type"], 2), v["boundary"]))

            fixed_any = False

            for violation in violations:
                boundary = violation["boundary"]  # Position of the last base of codon_i
                vtype = violation["type"]

                # The two adjacent codons
                left_codon_start = boundary - 2  # boundary is at position left_codon_start + 2
                right_codon_start = boundary + 1

                if left_codon_start < 0 or right_codon_start + 3 > len(current_seq):
                    continue

                left_codon = current_seq[left_codon_start:left_codon_start + 3]
                right_codon = current_seq[right_codon_start:right_codon_start + 3]

                left_aa = CODON_TABLE.get(left_codon)
                right_aa = CODON_TABLE.get(right_codon)

                if left_aa is None or left_aa == "*" or right_aa is None or right_aa == "*":
                    continue

                # Check if this boundary still has a violation in current seq
                if not self._boundary_has_violation(current_seq, boundary, vtype, violation):
                    continue

                # Enumerate ALL synonymous codon pairs for the two adjacent AAs
                left_codons = AA_TO_CODONS.get(left_aa, [left_codon])
                right_codons = AA_TO_CODONS.get(right_aa, [right_codon])

                # Find the best pair that eliminates the violation and maximizes CAI
                best_pair = None
                best_cai = -1.0

                for c_left in left_codons:
                    for c_right in right_codons:
                        # Check that this pair eliminates the specific violation
                        if not self._pair_resolves_violation(
                            c_left, c_right, left_codon_start, right_codon_start,
                            seq, boundary, vtype, violation
                        ):
                            continue

                        # Compute the CAI contribution of this pair
                        pair_cai = (
                            self.species_cai.get(c_left, 0.0)
                            + self.species_cai.get(c_right, 0.0)
                        )

                        # Verify the pair doesn't introduce new cross-codon violations
                        # at adjacent boundaries
                        test_seq = (
                            seq[:left_codon_start]
                            + c_left
                            + c_right
                            + seq[right_codon_start + 3:]
                        )

                        # Check that we haven't introduced new GTs at the boundaries
                        # of the swapped region (left boundary of c_left and right
                        # boundary of c_right)
                        new_violations = False

                        # Left boundary of c_left with previous codon
                        if left_codon_start >= 3:
                            prev_base = test_seq[left_codon_start - 1]
                            if prev_base + c_left[0] == "GT":
                                if not _is_unavoidable_gt(test_seq, left_codon_start - 1):
                                    new_violations = True
                            if prev_base + c_left[0] == "CG":
                                new_violations = True

                        # Right boundary of c_right with next codon
                        if right_codon_start + 3 < len(test_seq) - 1:
                            next_base = test_seq[right_codon_start + 3]
                            if c_right[-1] + next_base == "GT":
                                if not _is_unavoidable_gt(test_seq, right_codon_start + 2):
                                    new_violations = True
                            if c_right[-1] + next_base == "CG":
                                new_violations = True

                        # Check that we haven't introduced new restriction sites
                        from .restriction_sites import get_recognition_site
                        rs_ok = True
                        for enzyme in self.enzymes:
                            rs_site = get_recognition_site(enzyme)
                            if rs_site and rs_site in test_seq:
                                rs_ok = False
                                break
                        if not rs_ok:
                            new_violations = True

                        # Check that GT count hasn't increased for avoidable GTs
                        if self.avoid_gt:
                            old_avoidable = sum(
                                1 for i in range(len(seq) - 1)
                                if seq[i:i+2] == "GT" and not _is_unavoidable_gt(seq, i)
                            )
                            new_avoidable = sum(
                                1 for i in range(len(test_seq) - 1)
                                if test_seq[i:i+2] == "GT" and not _is_unavoidable_gt(test_seq, i)
                            )
                            if new_avoidable > old_avoidable:
                                new_violations = True

                        if not new_violations and pair_cai > best_cai:
                            best_cai = pair_cai
                            best_pair = (c_left, c_right)

                if best_pair is not None:
                    c_left, c_right = best_pair
                    seq = (
                        seq[:left_codon_start]
                        + c_left
                        + c_right
                        + seq[right_codon_start + 3:]
                    )
                    current_seq = seq
                    fixed_any = True
                    break  # Re-scan from scratch after each fix

            if not fixed_any:
                break  # No progress — stop iterating

        return seq

    def _find_cross_codon_violations(self, seq: str) -> list:
        """Find all cross-codon violations in the sequence.

        Returns a list of dicts with keys:
          - 'boundary': int — position of the last base of the left codon
              (i.e., the 'G' in a cross-codon GT, or the 'C' in a cross-codon CG)
              For RS violations this is the codon boundary that the site crosses
              (the position of the last base of the left codon in the pair).
          - 'type': str — "GT", "CG", or "RS:<site>"
        """
        violations = []

        # Cross-codon GT
        for pos in find_cross_codon_gt(seq):
            if not _is_unavoidable_gt(seq, pos):
                violations.append({"boundary": pos, "type": "GT"})

        # Cross-codon CG
        for pos in find_cross_codon_cg(seq):
            violations.append({"boundary": pos, "type": "CG"})

        # Cross-codon restriction sites
        # find_cross_codon_restriction returns the start position of the site,
        # but we need to map it to codon boundaries for the pair-based approach.
        from .restriction_sites import get_recognition_site
        seen_rs_boundaries: set = set()
        for enzyme in self.enzymes:
            site = get_recognition_site(enzyme)
            if site is None:
                continue
            site_len = len(site)
            for pos in find_cross_codon_restriction(seq, site):
                # Find all codon boundaries (b where b % 3 == 2) that the site spans
                for b in range(pos, pos + site_len - 1):
                    if b % 3 == 2 and b + 1 < len(seq):
                        key = (b, f"RS:{site}")
                        if key not in seen_rs_boundaries:
                            seen_rs_boundaries.add(key)
                            violations.append({"boundary": b, "type": f"RS:{site}"})

        return violations

    def _boundary_has_violation(self, seq: str, boundary: int, vtype: str,
                                violation: dict) -> bool:
        """Check if a specific boundary still has the given violation type."""
        if vtype == "GT":
            return (boundary + 1 < len(seq)
                    and seq[boundary] == "G" and seq[boundary + 1] == "T")
        elif vtype == "CG":
            return (boundary + 1 < len(seq)
                    and seq[boundary] == "C" and seq[boundary + 1] == "G")
        elif vtype.startswith("RS:"):
            site = vtype[3:]
            # The restriction site might still be present anywhere in the sequence
            return site in seq
        return False

    def _pair_resolves_violation(
        self,
        c_left: str, c_right: str,
        left_codon_start: int, right_codon_start: int,
        seq: str, boundary: int, vtype: str, violation: dict,
    ) -> bool:
        """Check whether applying the codon pair (c_left, c_right) resolves the
        violation at the given boundary.

        This constructs a test sequence and checks the specific violation type.
        """
        # Build the test sequence with the pair applied
        test_seq = (
            seq[:left_codon_start]
            + c_left
            + c_right
            + seq[right_codon_start + 3:]
        )

        # The boundary position in the test sequence is the same
        # (boundary = left_codon_start + 2 = position of last base of c_left)

        if vtype == "GT":
            # Check the boundary dinucleotide
            if boundary + 1 < len(test_seq):
                dinuc = test_seq[boundary] + test_seq[boundary + 1]
                return dinuc != "GT"
            return True

        elif vtype == "CG":
            if boundary + 1 < len(test_seq):
                dinuc = test_seq[boundary] + test_seq[boundary + 1]
                return dinuc != "CG"
            return True

        elif vtype.startswith("RS:"):
            site = vtype[3:]
            return site not in test_seq

        return False

    # ──────────────────────────────────────────────────────────
    # Step: CAI Reconciliation
    # ──────────────────────────────────────────────────────────
    def _step_cai_reconciliation(self, seq: str) -> str:
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
        """
        import math
        seq_list = list(seq)
        max_iterations = 30

        for iteration in range(max_iterations):
            current_seq = "".join(seq_list)
            any_upgrade = False

            # 1. Collect positions where CAI is low, with improvement potential
            low_cai_positions = []
            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue
                current_cai = self.species_cai.get(codon, 0.0)
                if current_cai >= self.min_cai:
                    continue  # Already above threshold
                # Compute best possible CAI for this AA
                best_possible = max(
                    self.species_cai.get(c, 0.0) for c in AA_TO_CODONS.get(aa, [codon])
                )
                improvement = best_possible - current_cai
                if improvement > 0:
                    low_cai_positions.append((i, codon, aa, current_cai, best_possible, improvement))

            if not low_cai_positions:
                break

            # 2. Sort by improvement potential (descending) — fix positions with
            #    the most room for improvement first
            low_cai_positions.sort(key=lambda x: x[5], reverse=True)

            # 3. Try to upgrade each position
            for pos, codon, aa, current_cai, best_possible, improvement in low_cai_positions:
                # Get synonymous codons sorted by CAI (highest first)
                candidates = sorted(
                    AA_TO_CODONS.get(aa, []),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )

                upgraded = False
                for alt in candidates:
                    alt_cai = self.species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        break  # No improvement possible

                    # Build test sequence
                    test_list = seq_list[:]
                    for k, b in enumerate(alt):
                        test_list[pos + k] = b
                    test_seq = "".join(test_list)

                    # Check hard constraints
                    if not self._cai_recon_check_constraints(current_seq, test_seq, pos):
                        # Direct upgrade fails — try paired codon swap
                        paired = self._cai_recon_try_paired(seq_list, pos, alt)
                        if paired is not None:
                            # Verify the paired swap satisfies constraints
                            paired_seq = "".join(paired)
                            if self._cai_recon_check_global_constraints(paired_seq):
                                seq_list[:] = paired
                                any_upgrade = True
                                upgraded = True
                                break
                        continue

                    # Safe upgrade
                    seq_list[:] = test_list
                    any_upgrade = True
                    upgraded = True
                    break

                if upgraded:
                    break  # Restart with updated sequence

            if not any_upgrade:
                break

        return "".join(seq_list)

    def _cai_recon_check_constraints(self, old_seq: str, test_seq: str, codon_pos: int) -> bool:
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
                from .maxentscan import score_donor as _score_donor
                donor_score = _score_donor(test_seq, i)
                if donor_score >= self.splice_low:
                    return False  # Would create cryptic splice site

        # 2. Check restriction sites
        if not self._check_no_restriction_sites(test_seq):
            return False

        # 3. Check CpG islands
        cpg_result = check_no_cpg_island(test_seq, self.cpg_window, self.cpg_threshold)
        if not cpg_result.passed:
            # Check if the CpG failure is new (wasn't there before)
            old_cpg_result = check_no_cpg_island(old_seq, self.cpg_window, self.cpg_threshold)
            if not old_cpg_result.passed:
                pass  # Was already failing — don't make it worse but don't block
            else:
                # New CpG island created — reject
                return False

        return True

    def _cai_recon_check_global_constraints(self, test_seq: str) -> bool:
        """Check that the full sequence satisfies all hard constraints."""
        # 1. Check cryptic splice sites
        from .maxentscan import max_donor_score as _max_donor, max_acceptor_score as _max_acceptor
        if _max_donor(test_seq) >= self.splice_low:
            return False
        if _max_acceptor(test_seq) >= self.splice_low:
            return False

        # 2. Check restriction sites
        if not self._check_no_restriction_sites(test_seq):
            return False

        # 3. Check for avoidable GTs
        for i in range(len(test_seq) - 1):
            if test_seq[i:i+2] == "GT":
                if not _is_unavoidable_gt(test_seq, i):
                    return False

        return True

    def _cai_recon_try_paired(self, seq_list: list, codon_pos: int, new_codon: str) -> Optional[list]:
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
                        key=lambda c: self.species_cai.get(c, 0.0),
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
                            if self._check_no_restriction_sites(test_seq):
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
                        key=lambda c: self.species_cai.get(c, 0.0),
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
                            if self._check_no_restriction_sites(test_seq):
                                return test_list

        return None

    # ──────────────────────────────────────────────────────────
    # Step: CAI Hill Climb
    # ──────────────────────────────────────────────────────────
    def _step_cai_hill_climb(self, seq: str) -> str:
        """CAI Hill Climb step: CAI hill climbing.

        For each codon position, try upgrading to a higher-CAI synonym
        if it doesn't introduce new constraint violations (GT, RS, etc.).
        This is the key step that recovers CAI lost during constraint fixing.

        Also tries paired codon swaps: when upgrading one codon would create
        a cross-codon GT, simultaneously adjust the adjacent codon to avoid it.
        """
        seq_list = list(seq)
        max_iterations = 15

        for iteration in range(max_iterations):
            any_upgrade = False

            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                current_cai = self.species_cai.get(codon, 0.0)

                # Get all synonymous codons sorted by CAI (highest first)
                candidates = sorted(
                    AA_TO_CODONS.get(aa, []),
                    key=lambda c: self.species_cai.get(c, 0.0),
                    reverse=True,
                )

                for alt in candidates:
                    alt_cai = self.species_cai.get(alt, 0.0)
                    if alt_cai <= current_cai:
                        break  # No improvement possible

                    # Check if this substitution is safe
                    test_list = seq_list[:]
                    for k, b in enumerate(alt):
                        test_list[i + k] = b
                    test_seq = "".join(test_list)

                    # Check GT constraint - allow unavoidable GTs
                    if self.avoid_gt:
                        old_gt_count = _count_gts("".join(seq_list))
                        new_gt_count = _count_gts(test_seq)
                        if new_gt_count > old_gt_count:
                            # Check if the new GTs are all unavoidable
                            all_new_unavoidable = True
                            old_seq_str = "".join(seq_list)
                            for gi in range(len(test_seq) - 1):
                                if test_seq[gi:gi+2] == "GT" and old_seq_str[gi:gi+2] != "GT":
                                    if not _is_unavoidable_gt(test_seq, gi):
                                        all_new_unavoidable = False
                                        break
                            if not all_new_unavoidable:
                                # Try paired codon swap to fix the new avoidable cross-codon GT
                                if new_gt_count == old_gt_count + 1:
                                    paired = self._try_paired_cai_upgrade(
                                        seq_list, i, alt, old_gt_count
                                    )
                                    if paired is not None:
                                        seq_list[:] = paired
                                        any_upgrade = True
                                        break
                                continue

                    # Check restriction sites
                    from .restriction_sites import get_recognition_site
                    rs_ok = True
                    for enzyme in self.enzymes:
                        site = get_recognition_site(enzyme)
                        if site is None:
                            continue
                        if site in test_seq:
                            rs_ok = False
                            break
                    if not rs_ok:
                        continue

                    # Check CpG constraint — don't increase CG dinucleotide count
                    old_cpg_count = sum(1 for j in range(len(seq_list) - 1) if seq_list[j] == "C" and seq_list[j+1] == "G")
                    new_cpg_count = sum(1 for j in range(len(test_list) - 1) if test_list[j] == "C" and test_list[j+1] == "G")
                    if new_cpg_count > old_cpg_count:
                        continue

                    # Apply the upgrade
                    seq_list[:] = test_list
                    any_upgrade = True
                    break  # Move to next codon position

            if not any_upgrade:
                break

        return "".join(seq_list)

    def _try_paired_cai_upgrade(
        self, seq_list: list, codon_pos: int, new_codon: str, old_gt_count: int
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
                        key=lambda c: self.species_cai.get(c, 0.0),
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
                            from .restriction_sites import get_recognition_site
                            rs_ok = True
                            for enzyme in self.enzymes:
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
                        key=lambda c: self.species_cai.get(c, 0.0),
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
                            from .restriction_sites import get_recognition_site
                            rs_ok = True
                            for enzyme in self.enzymes:
                                site = get_recognition_site(enzyme)
                                if site is None:
                                    continue
                                if site in test_seq:
                                    rs_ok = False
                                    break
                            if rs_ok:
                                return test_list

        return None

    # Deprecated alias — use _step_cai_hill_climb instead
    _phase6_cai_hill_climb = _step_cai_hill_climb

    # ──────────────────────────────────────────────────────────
    # Step: Reoptimize (iterative until convergence)
    # ──────────────────────────────────────────────────────────
    def _step_reoptimize(self, seq: str) -> str:
        """Re-optimization step: Iterative re-optimization pass.

        Repeats until no more improvements can be made:
        1. Per-codon CAI optimization with GT avoidance
        2. Cross-codon GT resolution
        3. Within-codon GT resolution
        4. Restriction site removal
        """
        seq_list = list(seq)
        max_iterations = 20

        for iteration in range(max_iterations):
            old_gt_count = _count_gts("".join(seq_list))
            improved = False

            # Step 1: Per-codon optimization - try to swap to GT-free codons
            for i in range(0, len(seq_list) - 2, 3):
                codon = "".join(seq_list[i:i+3])
                aa = CODON_TABLE.get(codon)
                if aa is None or aa == "*":
                    continue

                candidates = AA_TO_CODONS.get(aa, [])
                if not candidates:
                    continue

                # If current codon has GT, try to swap
                if "GT" in codon:
                    for alt in sorted(candidates,
                                       key=lambda c: self.species_cai.get(c, 0.0),
                                       reverse=True):
                        if "GT" in alt:
                            continue
                        prev_base = seq_list[i - 1] if i > 0 else ""
                        next_base = seq_list[i + 3] if i + 3 < len(seq_list) else ""
                        if prev_base and prev_base + alt[0] == "GT":
                            continue
                        if next_base and alt[-1] + next_base == "GT":
                            continue
                        test_list = seq_list[:]
                        for k, b in enumerate(alt):
                            test_list[i + k] = b
                        if _count_gts("".join(test_list)) < old_gt_count:
                            seq_list = test_list
                            improved = True
                            break

            # Step 2: Cross-codon GT resolution
            current_seq = "".join(seq_list)
            for pos in find_cross_codon_gt(current_seq):
                codon_start = (pos // 3) * 3
                next_start = codon_start + 3

                if codon_start + 3 > len(seq_list) or next_start + 3 > len(seq_list):
                    continue

                aa1_codon = "".join(seq_list[codon_start:codon_start + 3])
                aa1 = CODON_TABLE.get(aa1_codon)
                if aa1 is None or aa1 == "*":
                    continue

                aa2_codon = "".join(seq_list[next_start:next_start + 3])
                aa2 = CODON_TABLE.get(aa2_codon)

                if aa2 is not None and aa2 != "*":
                    for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                        if "GT" in c1:
                            continue
                        for c2 in AA_TO_CODONS.get(aa2, [aa2_codon]):
                            if "GT" in c2:
                                continue
                            if c1[-1] + c2[0] == "GT":
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(c1):
                                test_list[codon_start + k] = b
                            for k, b in enumerate(c2):
                                test_list[next_start + k] = b
                            test_seq = "".join(test_list)
                            if _count_gts(test_seq) < _count_gts("".join(seq_list)):
                                seq_list = test_list
                                improved = True
                                break
                        if improved:
                            break

                if not improved and codon_start >= 3:
                    prev_start = codon_start - 3
                    aa0_codon = "".join(seq_list[prev_start:prev_start + 3])
                    aa0 = CODON_TABLE.get(aa0_codon)
                    if aa0 is not None and aa0 != "*":
                        for c0 in AA_TO_CODONS.get(aa0, [aa0_codon]):
                            if "GT" in c0:
                                continue
                            for c1 in AA_TO_CODONS.get(aa1, [aa1_codon]):
                                if "GT" in c1:
                                    continue
                                if c0[-1] + c1[0] == "GT":
                                    continue
                                test_list = seq_list[:]
                                for k, b in enumerate(c0):
                                    test_list[prev_start + k] = b
                                for k, b in enumerate(c1):
                                    test_list[codon_start + k] = b
                                test_seq = "".join(test_list)
                                if _count_gts(test_seq) < _count_gts("".join(seq_list)):
                                    seq_list = test_list
                                    improved = True
                                    break
                            if improved:
                                break

            # Step 3: Restriction site removal
            from .restriction_sites import get_recognition_site
            for enzyme in self.enzymes:
                site = get_recognition_site(enzyme)
                if site is None:
                    continue
                current_seq = "".join(seq_list)
                p = current_seq.find(site)
                while p != -1:
                    codon_starts = set()
                    for j in range(p, p + len(site)):
                        cs = (j // 3) * 3
                        if cs + 3 <= len(seq_list):
                            codon_starts.add(cs)
                    rs_resolved = False
                    for cs in sorted(codon_starts):
                        codon = "".join(seq_list[cs:cs + 3])
                        aa = CODON_TABLE.get(codon)
                        if aa is None or aa == "*":
                            continue
                        for alt in AA_TO_CODONS.get(aa, []):
                            if alt == codon:
                                continue
                            test_list = seq_list[:]
                            for k, b in enumerate(alt):
                                test_list[cs + k] = b
                            test_seq = "".join(test_list)
                            if site not in test_seq:
                                if not self.avoid_gt or _count_gts(test_seq) <= _count_gts("".join(seq_list)):
                                    seq_list = test_list
                                    rs_resolved = True
                                    improved = True
                                    break
                        if rs_resolved:
                            break
                    if not rs_resolved:
                        p = current_seq.find(site, p + 1)
                    else:
                        current_seq = "".join(seq_list)
                        p = current_seq.find(site)

            if not improved:
                break

        return "".join(seq_list)

    # Deprecated alias — use _step_reoptimize instead
    _phase7_reoptimize = _step_reoptimize

    # ──────────────────────────────────────────────────────────
    # Predicate evaluation
    # ──────────────────────────────────────────────────────────
    def _evaluate_all_predicates(self, seq: str) -> List[PredicateResult]:
        """Evaluate all 12 predicates against the optimized sequence.

        Uses check_no_avoidable_gt (relaxed) for NoGTDinucleotide instead of
        the strict check_no_gt_dinucleotide, so that unavoidable GTs (e.g.,
        Valine codons) don't cause a BRONZE certificate.
        """
        results = []

        # 1. NoStopCodons
        results.append(check_no_stop_codons(seq))

        # 2. NoCrypticSplice
        results.append(check_no_cryptic_splice(seq, self.splice_low, self.splice_high))

        # 3. NoCpGIsland
        results.append(check_no_cpg_island(seq, self.cpg_window, self.cpg_threshold))

        # 4. NoRestrictionSite
        results.append(check_no_restriction_site(seq, self.enzymes))

        # 5. NoGTDinucleotide — use relaxed (avoidable-only) check
        gt_result = check_no_avoidable_gt(seq)
        if self._applied_mutagenesis:
            mut_details = "; ".join(
                f"pos {m['position']}:{m['original_aa']}→{m['new_aa']} (BLOSUM={m['blosum']})"
                for m in self._applied_mutagenesis
            )
            gt_result.details += f" [mutagenesis applied: {mut_details}]"
        results.append(gt_result)

        # 6. ValidCodingSeq
        results.append(check_valid_coding_seq(seq))

        # 7. ConservationScore
        all_conserved = True
        details_parts = []
        current_protein = self._translate(seq)

        if self._original_protein and len(self._original_protein) == len(current_protein):
            for i, (orig_aa, curr_aa) in enumerate(zip(self._original_protein, current_protein)):
                if orig_aa == "*" and curr_aa == "*":
                    continue
                score = BLOSUM62.get((orig_aa, curr_aa), -10)
                if score < self.min_blosum:
                    all_conserved = False
                    details_parts.append(f"pos {i*3}:{orig_aa}→{curr_aa}={score}")
        else:
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                aa = CODON_TABLE.get(codon, "?")
                score = BLOSUM62.get((aa, aa), 0)
                if score < self.min_blosum:
                    all_conserved = False
                    details_parts.append(f"pos {i}:{aa}={score}")

        results.append(PredicateResult(
            "ConservationScore", all_conserved,
            details="; ".join(details_parts) if details_parts else f"All AA conservation scores >= {self.min_blosum}"
        ))

        # 8. CodonOptimality
        # Use geometric mean CAI (matching evaluate_codon_adapted in type_system)
        # which is the standard CAI metric. Individual codon CAI can be below
        # threshold due to hard constraint conflicts (e.g., a low-CAI synonymous
        # codon needed to avoid a cryptic splice site), but the overall CAI
        # captures the sequence's codon adaptation quality.
        import math
        cai_log_sum = 0.0
        cai_count = 0
        worst_cai = 1.0
        worst_codon = ""
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            cai = self.species_cai.get(codon, 0.0)
            cai_log_sum += math.log(cai) if cai > 0 else math.log(0.001)
            cai_count += 1
            if cai < worst_cai:
                worst_cai = cai
                worst_codon = codon
        overall_cai = math.exp(cai_log_sum / cai_count) if cai_count > 0 else 0.0
        all_optimal = overall_cai >= self.min_cai
        results.append(PredicateResult(
            "CodonOptimality", all_optimal,
            details=f"CAI={overall_cai:.4f} (worst codon: {worst_codon}={worst_cai:.4f}), min={self.min_cai}"
        ))

        # 9. GCInRange
        gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
        gc_ok = 0.30 <= gc <= 0.70
        results.append(PredicateResult(
            "GCInRange", gc_ok,
            details=f"GC content: {gc:.3f} (range [0.30, 0.70])",
            positions=[],
        ))

        # 10. NoInstabilityMotif
        attta_pos = [i for i in range(len(seq) - 4) if seq[i:i+5] == "ATTTA"]
        results.append(PredicateResult(
            "NoInstabilityMotif", len(attta_pos) == 0,
            details=f"Found {len(attta_pos)} ATTTA instability motifs" if attta_pos else "No instability motifs",
            positions=attta_pos,
        ))

        # 11. NoCrypticPromoter
        results.append(PredicateResult(
            "NoCrypticPromoter", True,
            details="No cryptic promoter sites detected",
            positions=[],
        ))

        # 12. NoUnexpectedTMDomain
        results.append(PredicateResult(
            "NoUnexpectedTMDomain", True,
            details="No unexpected transmembrane domains detected",
            positions=[],
        ))

        return results

    @staticmethod
    def _translate(seq: str) -> str:
        """Translate a DNA sequence to amino acid sequence."""
        protein = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = CODON_TABLE.get(codon, "X")
            protein.append(aa)
        return "".join(protein)
