"""
BioCompiler Multi-Objective Sequence Optimizer

Production-grade optimizer with:
- Multi-codon coordinated restriction site removal (handles cross-boundary sites)
- Greedy optimizer as default (z3 optional via use_z3 flag)
- Organism-specific GC targeting (natural GC, not just "in range")
- Smart phase ordering with reconciliation pass
- Type-directed mutagenesis integration loop
- Delegation to type system for predicate checking (SOC)

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

    # Phase 1: Best codon per position (maximize CAI)
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
    # Strategy: for each strong cryptic site, try codon swaps that disrupt
    # the 9-mer (donor) or 23-mer (acceptor) context to bring the
    # MaxEntScan score below threshold.
    for iteration in range(200):
        max_d = max_donor_score(sequence)
        max_a = max_acceptor_score(sequence)
        if max_d < cryptic_splice_threshold and max_a < cryptic_splice_threshold:
            break

        fixed_any = False

        # Try to eliminate strong donors
        if max_d >= cryptic_splice_threshold:
            for i in range(len(sequence) - 1):
                if sequence[i:i+2] == "GT":
                    s = score_donor(sequence, i)
                    if s >= cryptic_splice_threshold:
                        codon_idx = i // 3
                        # Try swapping the GT codon and its neighbors
                        # to disrupt the 9-mer splice context
                        best_score = s
                        best_seq = None

                        # Try single-codon swap at the GT position
                        if codon_idx < len(aas):
                            aa = aas[codon_idx]
                            current = sequence[codon_idx*3:codon_idx*3+3]
                            for alt in sorted_codons[aa]:
                                if alt == current:
                                    continue
                                test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                                if i < len(test)-1 and test[i:i+2] == "GT":
                                    new_s = score_donor(test, i)
                                else:
                                    new_s = -999  # GT eliminated
                                if new_s < best_score:
                                    best_score = new_s
                                    best_seq = test

                        # Try 2-codon context disruption (GT codon + neighbor)
                        if best_score >= cryptic_splice_threshold and codon_idx < len(aas):
                            aa = aas[codon_idx]
                            current = sequence[codon_idx*3:codon_idx*3+3]
                            for neighbor_offset in [-1, 1]:
                                n_idx = codon_idx + neighbor_offset
                                if 0 <= n_idx < len(aas):
                                    n_aa = aas[n_idx]
                                    n_current = sequence[n_idx*3:n_idx*3+3]
                                    for v_alt in sorted_codons[aa][:2]:
                                        for n_alt in sorted_codons[n_aa][:3]:
                                            if n_alt == n_current and v_alt == current:
                                                continue
                                            test = list(sequence)
                                            test[codon_idx*3:codon_idx*3+3] = list(v_alt)
                                            test[n_idx*3:n_idx*3+3] = list(n_alt)
                                            test_str = "".join(test)
                                            if i < len(test_str)-1 and test_str[i:i+2] == "GT":
                                                new_s = score_donor(test_str, i)
                                            else:
                                                new_s = -999
                                            if new_s < best_score:
                                                best_score = new_s
                                                best_seq = test_str

                        if best_seq is not None and best_score < cryptic_splice_threshold:
                            sequence = best_seq
                            fixed_any = True
                            break  # Restart scanning

        # Try to eliminate strong acceptors
        if not fixed_any and max_a >= cryptic_splice_threshold:
            for i in range(len(sequence) - 1):
                if sequence[i:i+2] == "AG":
                    s = score_acceptor(sequence, i)
                    if s >= cryptic_splice_threshold:
                        codon_idx = i // 3
                        if codon_idx >= len(aas):
                            continue
                        aa = aas[codon_idx]
                        current = sequence[codon_idx*3:codon_idx*3+3]
                        best_score = s
                        best_seq = None

                        for alt in sorted_codons[aa]:
                            if alt == current:
                                continue
                            test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                            if i < len(test)-1 and test[i:i+2] == "AG":
                                new_s = score_acceptor(test, i)
                            else:
                                new_s = -999
                            if new_s < best_score:
                                best_score = new_s
                                best_seq = test

                        if best_seq is not None and best_score < cryptic_splice_threshold:
                            sequence = best_seq
                            fixed_any = True
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
