"""
BioCompiler Optimization Helpers
=================================
Restriction site removal helpers, GT/AG-free codon utilities,
and IUPAC expansion for the optimization pipeline.

Extracted from optimization.py for maintainability.
"""

import logging
import math
from itertools import product as itertools_product

from ..type_system import AA_TO_CODONS
from ..constants import reverse_complement, IUPAC_EXPAND

__all__ = [
    "MAX_RESTRICTION_SITE_ITERATIONS",
    "MAX_IUPAC_SITE_ITERATIONS",
    "IUPAC_EXPANSION_CAP",
]

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Named Constants (duplicated from optimization.py for
# self-containment — these are small integers, not tables)
# ────────────────────────────────────────────────────────────
MAX_RESTRICTION_SITE_ITERATIONS: int = 100
MAX_IUPAC_SITE_ITERATIONS: int = 100
IUPAC_EXPANSION_CAP: int = 4096


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
    usage: dict[str, float] | None = None,
) -> tuple[str, bool]:
    """
    Try to remove a restriction site using multi-codon coordinated solving.

    When a site straddles codon boundaries (e.g., PstI CTG|CAG spans codons i and i+1),
    single-codon swaps fail because changing either codon alone doesn't eliminate the site.
    This function enumerates valid codon COMBINATIONS for all overlapping codons.

    When ``usage`` (a codon→CAI adaptiveness dict) is provided, the function
    ranks all viable codon combinations by CAI impact and picks the one that
    minimises CAI loss.  Without ``usage``, it returns the first combo that
    eliminates the site (legacy behaviour).

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

        if usage is not None:
            # CAI-aware: enumerate ALL combos, rank by CAI, pick the best
            best_seq: str | None = None
            best_cai_sum: float = float('-inf')

            for combo in itertools_product(*candidate_lists):
                # Build test sequence with this combo applied
                test = list(sequence)
                for idx, ci in enumerate(overlapping):
                    start = ci * 3
                    test[start:start + 3] = list(combo[idx])
                test_seq = "".join(test)

                # Check if site is eliminated
                if site_upper not in test_seq and (not site_rc or site_rc not in test_seq):
                    # Compute total CAI contribution of the changed codons
                    cai_sum = sum(
                        math.log(usage.get(combo[idx], 1e-10))
                        for idx, ci in enumerate(overlapping)
                    )
                    if cai_sum > best_cai_sum:
                        best_cai_sum = cai_sum
                        best_seq = test_seq

            if best_seq is not None:
                return best_seq, True
        else:
            # Legacy: return first combo that eliminates the site
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
# GT/AG-Free Codon Utilities
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


def _is_unavoidable_gt_aa(aa: str) -> bool:
    """Check if ALL synonymous codons for an amino acid contain GT.

    For amino acids like Valine (V), every codon starts with GT,
    making GT avoidance impossible at these positions. Such positions
    are "unavoidable GT" and should not count against optimization quality.

    Pre-conditions:
    - aa is a valid single-letter amino acid code

    Post-conditions:
    - returns True if all codons for aa contain "GT"
    - returns False if at least one GT-free alternative exists
    """
    codons = AA_TO_CODONS.get(aa, [])
    if not codons:
        return True  # Unknown AA, assume unavoidable
    return all("GT" in c for c in codons)


def _gt_free_cai_ratio(aa: str, usage: dict[str, float]) -> float:
    """Compute the CAI ratio of best GT-free codon vs best GT-containing codon.

    Returns the ratio of the best GT-free codon's CAI weight to the best
    GT-containing codon's CAI weight. A ratio < 0.5 means avoiding GT
    would sacrifice more than half the CAI at this position, which is
    considered excessive.

    For amino acids where ALL codons contain GT (like Valine), returns 0.0.
    For amino acids where NO codons contain GT, returns 1.0 (no restriction needed).

    Pre-conditions:
    - aa is a valid amino acid code
    - usage is a codon->CAI weight dict

    Post-conditions:
    - returns a float in [0.0, 1.0]
    """
    codons = AA_TO_CODONS.get(aa, [])
    if not codons:
        return 0.0

    gt_free = [c for c in codons if "GT" not in c]
    gt_containing = [c for c in codons if "GT" in c]

    if not gt_containing:
        return 1.0  # No GT-containing codons, no restriction needed
    if not gt_free:
        return 0.0  # All codons contain GT, unavoidable

    best_gt_free_cai = max(usage.get(c, 0.01) for c in gt_free)
    best_gt_containing_cai = max(usage.get(c, 0.01) for c in gt_containing)

    if best_gt_containing_cai <= 0:
        return 1.0
    return best_gt_free_cai / best_gt_containing_cai


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
