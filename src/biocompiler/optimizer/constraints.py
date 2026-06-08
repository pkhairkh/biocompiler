"""
Constraint handling helpers for the optimizer.

Contains restriction site removal, GT/AG-free codon finding,
organism name resolution, back-translation, and cross-codon
constraint detection utilities.
"""

from typing import Tuple

import math
from itertools import product as itertools_product

from ..type_system import CODON_TABLE, AA_TO_CODONS
from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism, SPECIES_SHORT_NAMES
from ..constants import reverse_complement, RESTRICTION_ENZYMES
from .cai import _count_dinucs_fast


# ────────────────────────────────────────────────────────────
# Restriction Site Removal Helpers
# ────────────────────────────────────────────────────────────

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
    if pos < 0:
        raise ValueError(f"Position must be non-negative, got {pos}")
    if site_len <= 0:
        raise ValueError(f"Site length must be positive, got {site_len}")
    if n_codons <= 0:
        raise ValueError(f"Number of codons must be positive, got {n_codons}")

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


# ────────────────────────────────────────────────────────────
# GT/AG-Free Codon Helpers
# ────────────────────────────────────────────────────────────

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


# ────────────────────────────────────────────────────────────
# Organism Name Resolution Helpers
# ────────────────────────────────────────────────────────────

def _organism_to_species_key(organism: str) -> str:
    """Map an organism name to the species key used in SPECIES_SHORT_NAMES.

    Delegates to :func:`~biocompiler.organisms.resolve_organism` for
    name resolution and then looks up the short key in
    :data:`~biocompiler.organisms.SPECIES_SHORT_NAMES`.

    .. deprecated::
        Use :func:`~biocompiler.organisms.resolve_organism` and
        CODON_ADAPTIVENESS_TABLES directly instead of mapping to
        SPECIES dict keys.  Retained for backward compatibility.
    """
    canonical = resolve_organism(organism, strict=False)
    key = SPECIES_SHORT_NAMES.get(canonical)
    if key:
        return key
    # Fallback: try the organism name directly in SPECIES_SHORT_NAMES values
    if organism in SPECIES_SHORT_NAMES.values():
        return organism
    # Default to ecoli
    return "ecoli"


def _species_key_to_organism(species_key: str) -> str:
    """Map a species key or organism alias to the canonical organism name
    used in CODON_ADAPTIVENESS_TABLES.

    Delegates to :func:`~biocompiler.organisms.resolve_organism` for
    name resolution, which accepts both short aliases (e.g. 'ecoli')
    and full canonical names (e.g. 'Escherichia_coli'), as well as
    display names ('E. coli') and abbreviated binomials ('e_coli').

    .. deprecated::
        Use :func:`~biocompiler.organisms.resolve_organism` directly
        instead of this wrapper.  It is retained for internal use only.
    """
    return resolve_organism(species_key, default="Homo_sapiens")


# ────────────────────────────────────────────────────────────
# Back-Translation Helpers
# ────────────────────────────────────────────────────────────

def _back_translate_protein(
    protein: str,
    species_key: str,
    strategy: str = "greedy",
    enzymes: list[str] | None = None,
    is_eukaryote: bool = True,
) -> str:
    """Back-translate a protein to DNA using highest-CAI codons.

    Uses CODON_ADAPTIVENESS_TABLES (the same table used by compute_cai) so
    the initial sequence starts with optimal codons that will be reflected
    in the final CAI score.

    When ``strategy='dp'`` (or automatically for sequences < 200 aa), a
    dynamic-programming approach is used that considers the top 3 codons
    per position and finds the globally optimal sequence that:
      - Maximizes CAI
      - Avoids common restriction enzyme recognition sites
      - Minimizes GT/AG dinucleotides (for eukaryotes)

    Args:
        protein: Amino acid sequence (1-letter codes).
        species_key: Short species key (e.g. 'ecoli', 'human').
        strategy: 'greedy' for per-position best-CAI; 'dp' for global
            optimisation with cross-codon effects.
        enzymes: Restriction enzyme names to avoid (only used for DP).
        is_eukaryote: If True, penalise GT/AG dinucleotides in DP.

    Returns:
        DNA sequence string.
    """
    organism = _species_key_to_organism(species_key)
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
    if adaptiveness is None:
        # Fallback to Homo_sapiens if organism not found
        adaptiveness = CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]

    # Decide whether to use DP
    # Explicit 'dp' strategy always uses DP; 'greedy' always skips DP.
    # For other strategies (e.g. 'constraint_first'), use DP for short
    # sequences where it's fast and produces a better starting point.
    if strategy == "dp":
        use_dp = True
    elif strategy == "greedy":
        use_dp = False
    else:
        use_dp = len(protein) < 200
    if use_dp:
        return _back_translate_protein_dp(
            protein, adaptiveness, enzymes=enzymes, is_eukaryote=is_eukaryote,
        )

    # Greedy: simply pick the highest-CAI codon per position
    codons = []
    for aa in protein:
        if aa == "*":
            codons.append("TAA")
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            codons.append("NNN")
            continue
        best = max(candidates, key=lambda c: adaptiveness.get(c, 0.0))
        codons.append(best)
    return "".join(codons)


# ────────────────────────────────────────────────────────────
# DP-based back-translation
# ────────────────────────────────────────────────────────────

# Maximum number of top codons to consider per amino-acid position in DP
_DP_TOP_K: int = 3
# Penalty multiplier for GT/AG dinucleotides (eukaryotes)
_GT_AG_PENALTY: float = 0.02
# Penalty for creating a restriction site (effectively -inf)
_RESTRICTION_PENALTY: float = 0.5


def _build_restriction_site_set(
    enzymes: list[str] | None,
) -> set[str]:
    """Build the set of restriction site sequences (fwd + rc) to avoid."""
    sites: set[str] = set()
    if enzymes is None:
        enzymes = list(RESTRICTION_ENZYMES.keys())
    for name in enzymes:
        seq = RESTRICTION_ENZYMES.get(name)
        if seq is None:
            continue
        # Skip IUPAC wildcard sites (contain N or other ambiguity codes)
        if any(b not in "ATCG" for b in seq.upper()):
            continue
        seq_upper = seq.upper()
        sites.add(seq_upper)
        rc = reverse_complement(seq_upper)
        if rc != seq_upper:
            sites.add(rc)
    return sites


def _contains_restriction_site(
    seq: str, sites: set[str], window: int
) -> bool:
    """Check whether *seq* contains any of the restriction *sites*.

    Only checks the last *window* characters (enough to catch any new site
    introduced by the most recent codon).
    """
    start = max(0, len(seq) - window)
    tail = seq[start:]
    for site in sites:
        if site in tail:
            return True
    return False


def _count_dinucleotides(seq: str, di: str) -> int:
    """Count occurrences of dinucleotide *di* in *seq*.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for single-pass counting; falls back to pure-Python str.find
    otherwise.
    """
    return _count_dinucs_fast(seq, di)[0]


def _back_translate_protein_dp(
    protein: str,
    adaptiveness: dict[str, float],
    enzymes: list[str] | None = None,
    is_eukaryote: bool = True,
) -> str:
    """DP-based back-translation that considers cross-codon effects.

    For each amino-acid position, the top K codons by CAI are considered.
    The DP state tracks the last two codons (6 nt) to evaluate:
      - GT / AG dinucleotide penalties (eukaryotes)
      - Restriction site avoidance

    The objective is to maximise the sum of log-adaptiveness values
    (equivalent to maximising the geometric mean = CAI) while penalising
    undesirable features.

    Complexity: O(N * K^3) where N = protein length, K = _DP_TOP_K.
    For N < 200 and K = 3 this is negligible (< 1 ms).
    """
    import math as _math

    K = _DP_TOP_K
    restriction_sites = _build_restriction_site_set(enzymes)
    # Maximum restriction site length — only need to check this many
    # trailing characters for newly introduced sites
    max_site_len = max((len(s) for s in restriction_sites), default=0)
    # Window of trailing sequence to check for restriction sites:
    # last codon (3) + overlap from previous codon (max_site_len - 1)
    rest_check_window = max_site_len + 2  # conservative

    # Pre-compute top-K codons per amino acid
    top_codons_per_aa: list[list[tuple[str, float]]] = []
    for aa in protein:
        if aa == "*":
            # Stop codon — only one choice needed
            top_codons_per_aa.append([("TAA", 1.0)])
            continue
        candidates = AA_TO_CODONS.get(aa, [])
        if not candidates:
            top_codons_per_aa.append([("NNN", 0.0)])
            continue
        # Sort by adaptiveness descending, take top K
        sorted_cands = sorted(
            candidates,
            key=lambda c: adaptiveness.get(c, 0.0),
            reverse=True,
        )
        top = [(c, adaptiveness.get(c, 0.0)) for c in sorted_cands[:K]]
        top_codons_per_aa.append(top)

    n = len(protein)

    # DP state: (prev_codon_idx, curr_codon_idx) -> best log-CAI score
    # We track the last two codon choices to evaluate cross-codon effects.
    # prev_codon_idx refers to position i-2, curr_codon_idx to position i-1.
    # At position i, we extend to a new codon choice.

    # For position 0: no previous context
    dp: dict[tuple[int, int], tuple[float, list[int]]] = {}
    for ci, (codon_i, adapt_i) in enumerate(top_codons_per_aa[0]):
        score = _math.log(adapt_i) if adapt_i > 0 else -20.0
        dp[(0, ci)] = (score, [ci])

    if n == 1:
        # Trivial case
        best_state = max(dp, key=lambda k: dp[k][0])
        best_codon = top_codons_per_aa[0][dp[best_state][1][0]][0]
        return best_codon

    # Position 1: extend from position 0
    new_dp: dict[tuple[int, int], tuple[float, list[int]]] = {}
    for (prev_ci_key, curr_ci_key), (score, path) in dp.items():
        curr_ci = path[0]  # index into top_codons_per_aa[0]
        curr_codon = top_codons_per_aa[0][curr_ci][0]
        for ci1, (codon1, adapt1) in enumerate(top_codons_per_aa[1]):
            s = score + (_math.log(adapt1) if adapt1 > 0 else -20.0)
            # Check cross-codon effects between position 0 and 1
            junction = curr_codon + codon1  # 6 nt
            if is_eukaryote:
                # Penalise GT and AG at the cross-codon boundary
                if curr_codon[-1] + codon1[0] == "GT":
                    s -= _GT_AG_PENALTY * 50  # scale penalty to log space
                if curr_codon[-1] + codon1[0] == "AG":
                    s -= _GT_AG_PENALTY * 50
            # Check restriction sites in the junction
            if restriction_sites:
                for site in restriction_sites:
                    if len(site) <= len(junction) and site in junction:
                        s -= _RESTRICTION_PENALTY
            state = (curr_ci, ci1)
            if state not in new_dp or s > new_dp[state][0]:
                new_dp[state] = (s, path + [ci1])
    dp = new_dp

    # Positions 2..n-1: full DP with two-codon lookback
    for pos in range(2, n):
        new_dp = {}
        for (prev_ci, curr_ci), (score, path) in dp.items():
            prev_codon = top_codons_per_aa[pos - 2][prev_ci][0] if pos >= 2 else ""
            curr_codon = top_codons_per_aa[pos - 1][curr_ci][0]
            for ci, (codon, adapt) in enumerate(top_codons_per_aa[pos]):
                s = score + (_math.log(adapt) if adapt > 0 else -20.0)

                # Cross-codon GT/AG penalty (curr_codon | codon boundary)
                if is_eukaryote:
                    boundary = curr_codon[-1] + codon[0]
                    if boundary == "GT":
                        s -= _GT_AG_PENALTY * 50
                    elif boundary == "AG":
                        s -= _GT_AG_PENALTY * 50

                # Restriction site check: only need to check the region
                # that could contain a new site introduced by *codon*.
                if restriction_sites:
                    check_seq = curr_codon + codon
                    if prev_codon:
                        check_seq = prev_codon[-1] + check_seq
                    found_restriction = False
                    for site in restriction_sites:
                        if site in check_seq:
                            found_restriction = True
                            break
                    if found_restriction:
                        s -= _RESTRICTION_PENALTY

                state = (curr_ci, ci)
                if state not in new_dp or s > new_dp[state][0]:
                    new_dp[state] = (s, path + [ci])
        dp = new_dp

    # Extract best path
    if not dp:
        # Fallback: greedy
        codons = []
        for aa in protein:
            if aa == "*":
                codons.append("TAA")
                continue
            candidates = AA_TO_CODONS.get(aa, [])
            if not candidates:
                codons.append("NNN")
                continue
            best = max(candidates, key=lambda c: adaptiveness.get(c, 0.0))
            codons.append(best)
        return "".join(codons)

    best_state = max(dp, key=lambda k: dp[k][0])
    best_path = dp[best_state][1]

    # Reconstruct sequence
    result_codons = []
    for pos, idx in enumerate(best_path):
        result_codons.append(top_codons_per_aa[pos][idx][0])
    return "".join(result_codons)


# ────────────────────────────────────────────────────────────
# GT Dinucleotide Helpers
# ────────────────────────────────────────────────────────────

def _count_gts(s: str) -> int:
    """Count GT dinucleotides in a sequence.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for single-pass counting; falls back to pure-Python otherwise.
    """
    return _count_dinucs_fast(s, "GT")[0]


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


__all__ = [
    "_find_site_in_sequence",
    "_get_overlapping_codons",
    "_remove_site_multicodon",
    "_find_gt_free_codons",
    "_find_ag_free_codons",
    "_organism_to_species_key",
    "_species_key_to_organism",
    "_back_translate_protein",
    "_back_translate_protein_dp",
    "_build_restriction_site_set",
    "_contains_restriction_site",
    "_count_dinucleotides",
    "_count_gts",
    "_is_unavoidable_gt",
    "_has_gt",
    "_codon_creates_boundary_gt",
    "_DP_TOP_K",
    "_GT_AG_PENALTY",
    "_RESTRICTION_PENALTY",
]
