"""
CSP Sequence Optimizer using z3 Constraint Solver

Uses the z3 theorem prover to FIND or GENERATE nucleotide sequences that
pass all type predicates. Each nucleotide position is encoded as a z3 Int
variable, and constraints enforce biological correctness while maximizing
Codon Adaptation Index (CAI).

Key design:
- Nucleotide encoding: 0=A, 1=C, 2=G, 3=T
- Hard constraints: correct amino acid encoding, InFrame, GCInRange,
  NoRestrictionSite, NoInstabilityMotif
- Objective: maximize CAI (sum of log relative adaptiveness)
- Post-hoc verification: NoCrypticSplice (checked with MaxEntScan)
"""

from typing import List, Tuple, Dict, Optional, Set
from dataclasses import dataclass, field

try:
    from z3 import (
        Int, Optimize, If, And, Or, Not, Implies, Sum, ArithRef,
        sat, unsat, unknown, BoolRef, ArithRef as Z3ArithRef,
        Solver, set_param,
    )
except ImportError:
    raise ImportError("z3-solver is required. Install with: pip install z3-solver")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BASE_MAP = {0: "A", 1: "C", 2: "G", 3: "T"}
BASE_REV = {"A": 0, "C": 1, "G": 2, "T": 3}

# Standard genetic code: codon -> amino acid (single-letter code)
CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

# ---------------------------------------------------------------------------
# Codon Usage Tables (relative adaptiveness values)
# Relative adaptiveness = freq(codon) / max_freq(codons for same AA)
# Values from Kazusa human codon usage database (simplified)
# ---------------------------------------------------------------------------
HUMAN_CODON_USAGE: Dict[str, float] = {
    # Phe
    "TTT": 0.46, "TTC": 1.00,
    # Leu
    "TTA": 0.08, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 1.00,
    # Ile
    "ATT": 0.36, "ATC": 1.00, "ATA": 0.14,
    # Met
    "ATG": 1.00,
    # Val
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 1.00,
    # Ser
    "TCT": 0.18, "TCC": 0.23, "TCA": 0.14, "TCG": 0.05,
    "AGT": 0.15, "AGC": 1.00,
    # Pro
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.42, "CCG": 1.00,
    # Thr
    "ACT": 0.25, "ACC": 1.00, "ACA": 0.30, "ACG": 0.12,
    # Ala
    "GCT": 0.27, "GCC": 1.00, "GCA": 0.27, "GCG": 0.10,
    # Tyr
    "TAT": 0.43, "TAC": 1.00,
    # His
    "CAT": 0.42, "CAC": 1.00,
    # Gln
    "CAA": 0.26, "CAG": 1.00,
    # Asn
    "AAT": 0.47, "AAC": 1.00,
    # Lys
    "AAA": 0.43, "AAG": 1.00,
    # Asp
    "GAT": 0.46, "GAC": 1.00,
    # Glu
    "GAA": 0.42, "GAG": 1.00,
    # Cys
    "TGT": 0.46, "TGC": 1.00,
    # Trp
    "TGG": 1.00,
    # Arg
    "CGT": 0.09, "CGC": 0.21, "CGA": 0.07, "CGG": 0.12,
    "AGA": 0.21, "AGG": 1.00,
    # Gly
    "GGT": 0.16, "GGC": 1.00, "GGA": 0.25, "GGG": 0.25,
}

E_COLI_CODON_USAGE: Dict[str, float] = {
    # Phe
    "TTT": 0.35, "TTC": 1.00,
    # Leu
    "TTA": 0.14, "TTG": 0.13, "CTT": 0.12, "CTC": 0.10, "CTA": 0.04, "CTG": 1.00,
    # Ile
    "ATT": 0.31, "ATC": 1.00, "ATA": 0.06,
    # Met
    "ATG": 1.00,
    # Val
    "GTT": 0.36, "GTC": 0.26, "GTA": 0.17, "GTG": 1.00,
    # Ser
    "TCT": 0.17, "TCC": 0.16, "TCA": 0.13, "TCG": 0.14,
    "AGT": 0.15, "AGC": 1.00,
    # Pro
    "CCT": 0.17, "CCC": 0.12, "CCA": 0.20, "CCG": 1.00,
    # Thr
    "ACT": 0.18, "ACC": 1.00, "ACA": 0.14, "ACG": 0.26,
    # Ala
    "GCT": 0.18, "GCC": 0.26, "GCA": 0.22, "GCG": 1.00,
    # Tyr
    "TAT": 0.42, "TAC": 1.00,
    # His
    "CAT": 0.44, "CAC": 1.00,
    # Gln
    "CAA": 0.30, "CAG": 1.00,
    # Asn
    "AAT": 0.39, "AAC": 1.00,
    # Lys
    "AAA": 0.74, "AAG": 1.00,
    # Asp
    "GAT": 0.60, "GAC": 1.00,
    # Glu
    "GAA": 0.70, "GAG": 1.00,
    # Cys
    "TGT": 0.44, "TGC": 1.00,
    # Trp
    "TGG": 1.00,
    # Arg
    "CGT": 0.37, "CGC": 0.36, "CGA": 0.07, "CGG": 0.11,
    "AGA": 0.07, "AGG": 1.00,
    # Gly
    "GGT": 0.35, "GGC": 1.00, "GGA": 0.10, "GGG": 0.15,
}

CODON_USAGE_TABLES = {
    "human": HUMAN_CODON_USAGE,
    "e_coli": E_COLI_CODON_USAGE,
}

# Reverse lookup: amino acid -> list of codons
AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    if _aa != "*":
        AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# Common restriction enzyme sites (5'->3')
COMMON_RESTRICTION_SITES: Dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "NdeI": "CATATG",
}

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------
@dataclass
class TypePredicate:
    """A type predicate that a DNA sequence must satisfy."""
    name: str
    description: str


@dataclass
class OptimizationResult:
    """Result of CSP optimization."""
    sequence: str
    protein: str
    cai: float
    gc_content: float
    satisfied_predicates: List[str]
    failed_predicates: List[str]
    unsat_core: Optional[List[str]] = None
    fallback_used: bool = False


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def protein_to_aa_list(protein: str) -> List[str]:
    """Convert protein string to list of amino acid single-letter codes."""
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    result = []
    for i, ch in enumerate(protein):
        if ch not in valid_aas:
            raise ValueError(f"Invalid amino acid '{ch}' at position {i}")
        result.append(ch)
    return result


def codon_to_indices(codon: str) -> Tuple[int, int, int]:
    """Convert a 3-letter codon string to (n0, n1, n2) indices."""
    return (BASE_REV[codon[0]], BASE_REV[codon[1]], BASE_REV[codon[2]])


def indices_to_codon(n0: int, n1: int, n2: int) -> str:
    """Convert (n0, n1, n2) indices to a 3-letter codon string."""
    return BASE_MAP[n0] + BASE_MAP[n1] + BASE_MAP[n2]


def compute_cai(sequence: str, organism: str = "human") -> float:
    """
    Compute the Codon Adaptation Index for a coding sequence.

    CAI = geometric mean of relative adaptiveness values for all codons.
    """
    usage = CODON_USAGE_TABLES.get(organism, HUMAN_CODON_USAGE)
    seq = sequence.upper()

    if len(seq) < 3 or len(seq) % 3 != 0:
        return 0.0

    log_sum = 0.0
    n_codons = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in usage:
            w = usage[codon]
            if w > 0:
                import math
                log_sum += math.log(w)
                n_codons += 1
            else:
                # Rare codon with 0 usage - give small penalty
                log_sum += math.log(0.01)
                n_codons += 1
        else:
            return 0.0

    if n_codons == 0:
        return 0.0

    import math
    cai = math.exp(log_sum / n_codons)
    return round(cai, 4)


def compute_gc_content(sequence: str) -> float:
    """Compute GC content as a fraction [0, 1]."""
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(gc / len(seq), 4)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq.upper()))


# ---------------------------------------------------------------------------
# Greedy solver (fallback for long sequences)
# ---------------------------------------------------------------------------
def _greedy_optimize(
    protein: str,
    organism: str = "human",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: Optional[List[str]] = None,
    cai_threshold: float = 0.2,
) -> str:
    """
    Greedy codon optimization: choose the highest-CAI codon at each position,
    then apply local fixes for constraint violations.

    Used as a fallback when z3 optimization is too slow for long sequences.
    """
    usage = CODON_USAGE_TABLES.get(organism, HUMAN_CODON_USAGE)
    aas = protein_to_aa_list(protein)
    restriction_sites = restriction_sites or []

    # Sort codons for each AA by usage (descending)
    sorted_codons: Dict[str, List[str]] = {}
    for aa in set(aas):
        codons = AA_TO_CODONS[aa]
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted

    # Phase 1: Choose best codon for each position
    seq_parts: List[str] = []
    for aa in aas:
        seq_parts.append(sorted_codons[aa][0])  # highest usage codon
    sequence = "".join(seq_parts)

    # Phase 2: Fix restriction sites by replacing codons
    for site in restriction_sites:
        site_upper = site.upper()
        site_rc = reverse_complement(site_upper)
        max_iterations = 100
        for _ in range(max_iterations):
            pos = sequence.find(site_upper)
            pos_rc = sequence.find(site_rc)
            if pos == -1 and pos_rc == -1:
                break
            # Fix the first occurrence
            target_pos = pos if pos != -1 else pos_rc
            codon_idx = target_pos // 3
            offset_in_codon = target_pos % 3
            # Try replacing the codon that contains the start of the site
            if codon_idx < len(aas):
                aa = aas[codon_idx]
                current_codon = sequence[codon_idx*3 : codon_idx*3+3]
                # Try alternative codons
                for alt_codon in sorted_codons[aa]:
                    if alt_codon != current_codon:
                        test_seq = sequence[:codon_idx*3] + alt_codon + sequence[codon_idx*3+3:]
                        if site_upper not in test_seq and site_rc not in test_seq:
                            sequence = test_seq
                            break
                else:
                    # Try replacing adjacent codon
                    if codon_idx + 1 < len(aas):
                        next_aa = aas[codon_idx + 1]
                        next_codon = sequence[(codon_idx+1)*3 : (codon_idx+1)*3+3]
                        for alt_codon in sorted_codons[next_aa]:
                            if alt_codon != next_codon:
                                test_seq = sequence[:(codon_idx+1)*3] + alt_codon + sequence[(codon_idx+1)*3+3:]
                                if site_upper not in test_seq and site_rc not in test_seq:
                                    sequence = test_seq
                                    break

    # Phase 3: Fix ATTTA instability motifs
    max_iterations = 100
    for _ in range(max_iterations):
        pos = sequence.find("ATTTA")
        if pos == -1:
            break
        codon_idx = pos // 3
        for ci in range(max(0, codon_idx - 1), min(len(aas), codon_idx + 2)):
            aa = aas[ci]
            current_codon = sequence[ci*3 : ci*3+3]
            for alt_codon in sorted_codons[aa]:
                if alt_codon != current_codon:
                    test_seq = sequence[:ci*3] + alt_codon + sequence[ci*3+3:]
                    if "ATTTA" not in test_seq:
                        sequence = test_seq
                        break
            else:
                continue
            break

    # Phase 4: Fix 6+ consecutive T
    max_iterations = 100
    for _ in range(max_iterations):
        # Find longest T run
        max_run = 0
        max_run_pos = -1
        i = 0
        while i < len(sequence):
            if sequence[i] == "T":
                j = i
                while j < len(sequence) and sequence[j] == "T":
                    j += 1
                run_len = j - i
                if run_len > max_run:
                    max_run = run_len
                    max_run_pos = i
                i = j
            else:
                i += 1

        if max_run < 6:
            break

        # Break the T run by replacing a codon
        codon_idx = (max_run_pos + max_run // 2) // 3
        if codon_idx < len(aas):
            aa = aas[codon_idx]
            current_codon = sequence[codon_idx*3 : codon_idx*3+3]
            for alt_codon in sorted_codons[aa]:
                if alt_codon != current_codon:
                    test_seq = sequence[:codon_idx*3] + alt_codon + sequence[codon_idx*3+3:]
                    # Check no 6+ T run
                    if not any(test_seq[i:i+6] == "TTTTTT" for i in range(len(test_seq) - 5)):
                        sequence = test_seq
                        break

    # Phase 5: Adjust GC content if needed
    gc = compute_gc_content(sequence)
    if gc < gc_lo or gc > gc_hi:
        # Try swapping codons to adjust GC
        for iteration in range(200):
            gc = compute_gc_content(sequence)
            if gc_lo <= gc <= gc_hi:
                break

            for ci in range(len(aas)):
                aa = aas[ci]
                current_codon = sequence[ci*3 : ci*3+3]
                current_gc = sum(1 for b in current_codon if b in "GC")

                best_alt = None
                best_gc_diff = abs(gc - (gc_lo + gc_hi) / 2)

                for alt_codon in sorted_codons[aa]:
                    if alt_codon == current_codon:
                        continue
                    alt_gc = sum(1 for b in alt_codon if b in "GC")
                    # Would this move GC in the right direction?
                    new_total_gc = sum(1 for b in sequence if b in "GC") - current_gc + alt_gc
                    new_gc_frac = new_total_gc / len(sequence)
                    diff = abs(new_gc_frac - (gc_lo + gc_hi) / 2)
                    if diff < best_gc_diff:
                        best_gc_diff = diff
                        best_alt = alt_codon

                if best_alt is not None:
                    sequence = sequence[:ci*3] + best_alt + sequence[ci*3+3:]

    return sequence


# ---------------------------------------------------------------------------
# z3-based CSP optimizer
# ---------------------------------------------------------------------------
def optimize_sequence(
    target_protein: str,
    organism: str = "human",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: Optional[List[str]] = None,
    cai_threshold: float = 0.2,
    max_amino_acids_for_z3: int = 80,
    z3_timeout_ms: int = 30000,
) -> OptimizationResult:
    """
    Find or generate a DNA sequence encoding the target protein that satisfies
    all type predicates, using z3 constraint solver with greedy fallback.

    Args:
        target_protein: amino acid sequence (single-letter codes)
        organism: 'human' or 'e_coli' for codon usage table
        gc_lo: minimum GC content fraction
        gc_hi: maximum GC content fraction
        restriction_sites: list of restriction enzyme site sequences to avoid
        cai_threshold: minimum relative adaptiveness per codon (0-1)
        max_amino_acids_for_z3: use z3 for sequences up to this length,
            fallback to greedy for longer sequences
        z3_timeout_ms: z3 solver timeout in milliseconds

    Returns:
        OptimizationResult with the optimized sequence and metadata
    """
    aas = protein_to_aa_list(target_protein)
    n_aa = len(aas)
    n_bases = n_aa * 3

    if n_aa == 0:
        return OptimizationResult(
            sequence="", protein="", cai=0.0, gc_content=0.0,
            satisfied_predicates=[], failed_predicates=["empty_sequence"],
        )

    usage = CODON_USAGE_TABLES.get(organism, HUMAN_CODON_USAGE)
    restriction_sites = restriction_sites or list(COMMON_RESTRICTION_SITES.values())

    # Choose solver based on sequence length
    fallback_used = False
    if n_aa > max_amino_acids_for_z3:
        sequence = _greedy_optimize(
            target_protein, organism, gc_lo, gc_hi, restriction_sites, cai_threshold
        )
        fallback_used = True
    else:
        # ---- z3 optimization ----
        opt = Optimize()
        set_param("timeout", z3_timeout_ms)

        # Create nucleotide variables: nuc[i] for i in 0..n_bases-1
        nuc = [Int(f"n_{i}") for i in range(n_bases)]

        # Constraint: each nucleotide is A(0), C(1), G(2), or T(3)
        for i in range(n_bases):
            opt.add(nuc[i] >= 0, nuc[i] <= 3)

        # Constraint: each codon encodes the correct amino acid + CAI objective
        cai_contributions = []

        for aa_idx, aa in enumerate(aas):
            base_idx = aa_idx * 3
            n0, n1, n2 = nuc[base_idx], nuc[base_idx + 1], nuc[base_idx + 2]

            valid_codons = AA_TO_CODONS[aa]
            codon_constraints = []
            cai_if_chain = None

            sorted_codons_z3 = sorted(valid_codons, key=lambda c: usage.get(c, 0.0), reverse=True)

            for codon in sorted_codons_z3:
                c0, c1, c2 = codon_to_indices(codon)
                codon_match = And(n0 == c0, n1 == c1, n2 == c2)
                codon_constraints.append(codon_match)

                w = usage.get(codon, 0.01)
                if cai_if_chain is None:
                    cai_if_chain = If(codon_match, w, 0.0)
                else:
                    cai_if_chain = If(codon_match, w, cai_if_chain)

            opt.add(Or(*codon_constraints))

            if cai_if_chain is not None:
                cai_contributions.append(cai_if_chain)

        # Objective: maximize sum of CAI contributions
        total_cai_sum = Sum(cai_contributions)
        opt.maximize(total_cai_sum)

        # GCInRange
        gc_count = Sum([If(Or(nuc[i] == 1, nuc[i] == 2), 1, 0) for i in range(n_bases)])
        gc_count_lo = int(gc_lo * n_bases)
        gc_count_hi = int(gc_hi * n_bases)
        opt.add(gc_count >= gc_count_lo, gc_count <= gc_count_hi)

        # NoRestrictionSite
        for site in restriction_sites:
            site_upper = site.upper()
            site_len = len(site_upper)
            site_indices = [BASE_REV[b] for b in site_upper]
            for start in range(n_bases - site_len + 1):
                window_match = And(
                    *[nuc[start + j] == site_indices[j] for j in range(site_len)]
                )
                opt.add(Not(window_match))

        # NoInstabilityMotif: no ATTTA
        attta_indices = [BASE_REV[b] for b in "ATTTA"]
        for start in range(n_bases - 5 + 1):
            window_match = And(
                *[nuc[start + j] == attta_indices[j] for j in range(5)]
            )
            opt.add(Not(window_match))

        # No 6+ consecutive T
        for start in range(n_bases - 6 + 1):
            all_t = And(*[nuc[start + j] == 3 for j in range(6)])
            opt.add(Not(all_t))

        # NoCrypticSplice (approximate): weaken GT donor context
        for i in range(n_bases - 1):
            is_gt = And(nuc[i] == 2, nuc[i + 1] == 3)
            # Force position +3 to be C or T (break donor consensus)
            if i + 3 < n_bases:
                opt.add(Implies(is_gt, Or(nuc[i + 2] == 1, nuc[i + 2] == 3)))
            # Force position -1 to not be C
            if i >= 1:
                opt.add(Implies(is_gt, Or(nuc[i - 1] == 0, nuc[i - 1] == 2, nuc[i - 1] == 3)))

        # Weaken AG acceptor context
        for i in range(n_bases - 1):
            is_ag = And(nuc[i] == 0, nuc[i + 1] == 2)
            if i >= 1:
                opt.add(Implies(is_ag, Or(nuc[i - 1] == 0, nuc[i - 1] == 2)))

        # Solve
        result = opt.check()

        if result == sat:
            model = opt.model()
            seq_bases = []
            for i in range(n_bases):
                val = model.eval(nuc[i], model_completion=True)
                int_val = val.as_long()
                seq_bases.append(BASE_MAP[int_val])
            sequence = "".join(seq_bases)
        else:
            # z3 UNSAT or timeout — fall back to greedy
            sequence = _greedy_optimize(
                target_protein, organism, gc_lo, gc_hi, restriction_sites, cai_threshold
            )
            fallback_used = True

    # ---- Post-processing: minimize cryptic splice sites (always) ----
    sequence = _minimize_cryptic_splice(sequence, target_protein, organism, gc_lo, gc_hi, restriction_sites)

    cai = compute_cai(sequence, organism)
    gc = compute_gc_content(sequence)
    satisfied, failed = _check_predicates(
        sequence, gc_lo, gc_hi, restriction_sites, cai_threshold, organism
    )
    return OptimizationResult(
        sequence=sequence, protein=target_protein, cai=cai,
        gc_content=gc, satisfied_predicates=satisfied,
        failed_predicates=failed, fallback_used=fallback_used,
    )


def _minimize_cryptic_splice(
    sequence: str,
    protein: str,
    organism: str,
    gc_lo: float,
    gc_hi: float,
    restriction_sites: List[str],
) -> str:
    """
    Post-processing: iteratively swap codons to reduce MaxEntScan scores
    below the cryptic splice site threshold (3.0).
    
    This is applied AFTER the initial optimization (z3 or greedy) to
    ensure NoCrypticSplice passes.
    """
    from maxentscan import max_donor_score, max_acceptor_score, scan_splice_sites
    
    aas = protein_to_aa_list(protein)
    usage = CODON_USAGE_TABLES.get(organism, HUMAN_CODON_USAGE)
    
    # Sort codons for each AA by usage (descending)
    sorted_codons: Dict[str, List[str]] = {}
    for aa in set(aas):
        codons = AA_TO_CODONS[aa]
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted
    
    for iteration in range(500):
        max_d = max_donor_score(sequence)
        max_a = max_acceptor_score(sequence)
        if max_d < 3.0 and max_a < 3.0:
            break
        
        current_worst = max(max_d, max_a)
        
        # Find high-scoring sites and only try codons near them
        sites = scan_splice_sites(sequence, 2.0, 2.0)  # lower threshold to find all relevant sites
        if not sites:
            break
        
        # Get codon indices near high-scoring sites (within 9 bases = 3 codons)
        codons_to_try = set()
        for pos, stype, score in sites:
            if score < 2.5:
                continue  # skip very low scores
            for offset in range(-3, 4):
                ci = (pos + offset) // 3
                if 0 <= ci < len(aas):
                    codons_to_try.add(ci)
        
        best_improvement = 0.0
        best_swap = None
        
        for ci in codons_to_try:
            aa = aas[ci]
            current_codon = sequence[ci*3 : ci*3+3]
            
            for alt_codon in sorted_codons[aa]:
                if alt_codon == current_codon:
                    continue
                
                test_seq = sequence[:ci*3] + alt_codon + sequence[ci*3+3:]
                
                # Must not violate other constraints
                has_restriction = False
                for site in restriction_sites:
                    if site.upper() in test_seq or reverse_complement(site.upper()) in test_seq:
                        has_restriction = True
                        break
                if has_restriction:
                    continue
                if "ATTTA" in test_seq:
                    continue
                if any(test_seq[i:i+6] == "TTTTTT" for i in range(len(test_seq) - 5)):
                    continue
                new_gc = compute_gc_content(test_seq)
                if not (gc_lo <= new_gc <= gc_hi):
                    continue
                
                new_max_d = max_donor_score(test_seq)
                new_max_a = max_acceptor_score(test_seq)
                new_worst = max(new_max_d, new_max_a)
                
                improvement = current_worst - new_worst
                if improvement > best_improvement:
                    best_improvement = improvement
                    best_swap = (ci, alt_codon)
        
        if best_swap is None or best_improvement <= 0:
            break
        
        ci, alt_codon = best_swap
        sequence = sequence[:ci*3] + alt_codon + sequence[ci*3+3:]
    
    return sequence


def _check_predicates(
    sequence: str,
    gc_lo: float,
    gc_hi: float,
    restriction_sites: List[str],
    cai_threshold: float,
    organism: str,
) -> Tuple[List[str], List[str]]:
    """Check all type predicates against the optimized sequence."""
    satisfied: List[str] = []
    failed: List[str] = []

    # InFrame: always satisfied by codon-based construction
    satisfied.append("InFrame")

    # GCInRange
    gc = compute_gc_content(sequence)
    if gc_lo <= gc <= gc_hi:
        satisfied.append("GCInRange")
    else:
        failed.append(f"GCInRange(gc={gc:.3f}, range=[{gc_lo},{gc_hi}])")

    # NoRestrictionSite
    has_restriction = False
    for site in restriction_sites:
        site_upper = site.upper()
        site_rc = reverse_complement(site_upper)
        if site_upper in sequence or site_rc in sequence:
            has_restriction = True
            break
    if not has_restriction:
        satisfied.append("NoRestrictionSite")
    else:
        failed.append("NoRestrictionSite")

    # NoInstabilityMotif
    has_atta = "ATTTA" in sequence
    has_t6 = any(sequence[i:i+6] == "TTTTTT" for i in range(len(sequence) - 5))
    if not has_atta and not has_t6:
        satisfied.append("NoInstabilityMotif")
    else:
        reasons = []
        if has_atta:
            reasons.append("ATTTA found")
        if has_t6:
            reasons.append("6+ consecutive T")
        failed.append(f"NoInstabilityMotif({', '.join(reasons)})")

    # CodonAdapted
    cai = compute_cai(sequence, organism)
    if cai >= cai_threshold:
        satisfied.append("CodonAdapted")
    else:
        failed.append(f"CodonAdapted(cai={cai:.4f}, threshold={cai_threshold})")

    # NoCrypticSplice (checked with MaxEntScan — import here to avoid circular)
    from maxentscan import max_donor_score, max_acceptor_score
    max_donor = max_donor_score(sequence)
    max_acceptor = max_acceptor_score(sequence)
    if max_donor < 3.0 and max_acceptor < 3.0:
        satisfied.append("NoCrypticSplice")
    else:
        failed.append(
            f"NoCrypticSplice(max_donor={max_donor:.2f}, max_acceptor={max_acceptor:.2f})"
        )

    return satisfied, failed


if __name__ == "__main__":
    # Quick test with a small protein
    test_protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    print(f"Testing with Human Insulin ({len(test_protein)} aa):")
    result = optimize_sequence(test_protein, organism="human")
    print(f"  Sequence length: {len(result.sequence)} bp")
    print(f"  CAI: {result.cai}")
    print(f"  GC content: {result.gc_content:.1%}")
    print(f"  Satisfied: {result.satisfied_predicates}")
    print(f"  Failed: {result.failed_predicates}")
    print(f"  Fallback used: {result.fallback_used}")
