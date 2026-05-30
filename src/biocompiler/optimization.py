"""
BioCompiler CSP Sequence Optimizer — z3 + Greedy Fallback

Production-grade optimizer with:
- Thread-safe z3 usage (local solver timeout, not global set_param)
- Convergence checking for greedy optimizer loops
- Geometric mean CAI objective for z3 (sum-of-logs, not sum)
- Configurable cryptic splice threshold
- Clear failure reporting
"""

import logging
import math
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field

from .constants import (
    CODON_TABLE, AA_TO_CODONS, BASE_MAP, BASE_REV,
    RESTRICTION_ENZYMES, reverse_complement,
)
from .organisms import CODON_USAGE_TABLES, SUPPORTED_ORGANISMS, CODON_ADAPTIVENESS_TABLES
from .exceptions import UnsupportedOrganismError, InvalidProteinError, OptimizationError
from .translation import compute_cai
from .scanner import gc_content
from .maxentscan import max_donor_score, max_acceptor_score

logger = logging.getLogger(__name__)


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


def protein_to_aa_list(protein: str) -> List[str]:
    """Convert protein string to list of amino acid codes. Raises InvalidProteinError for bad input."""
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    invalid = set(ch for ch in protein if ch not in valid_aas)
    if invalid:
        raise InvalidProteinError(protein, invalid)
    return list(protein)


def _greedy_optimize(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: Optional[List[str]] = None,
    cryptic_splice_threshold: float = 3.0,
) -> Tuple[str, List[str]]:
    """
    Greedy codon optimization: choose highest-CAI codon, then fix violations.

    Returns:
        Tuple of (optimized_sequence, list_of_warnings)
    """
    usage = CODON_ADAPTIVENESS_TABLES.get(organism)
    if usage is None:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)
    aas = protein_to_aa_list(protein)
    restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.values())
    warnings: List[str] = []

    sorted_codons: Dict[str, List[str]] = {}
    for aa in set(aas):
        codons = AA_TO_CODONS[aa]
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted

    # Phase 1: Best codon per position
    sequence = "".join(sorted_codons[aa][0] for aa in aas)

    # Phase 2: Fix restriction sites
    for site in restriction_sites:
        site_upper = site.upper()
        site_rc = reverse_complement(site_upper)
        for iteration in range(100):
            pos = sequence.find(site_upper)
            pos_rc = sequence.find(site_rc)
            if pos == -1 and pos_rc == -1:
                break
            target_pos = pos if pos != -1 else pos_rc
            codon_idx = target_pos // 3
            if codon_idx < len(aas):
                aa = aas[codon_idx]
                current = sequence[codon_idx*3: codon_idx*3+3]
                fixed = False
                for alt in sorted_codons[aa]:
                    if alt != current:
                        test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                        if site_upper not in test and site_rc not in test:
                            sequence = test
                            fixed = True
                            break
                if not fixed:
                    warnings.append(f"Cannot remove {site} at iteration {iteration}")
                    break
        else:
            warnings.append(f"Restriction site {site}: max iterations reached, may still be present")

    # Phase 3: Fix ATTTA
    for iteration in range(100):
        pos = sequence.find("ATTTA")
        if pos == -1:
            break
        codon_idx = pos // 3
        fixed = False
        for ci in range(max(0, codon_idx-1), min(len(aas), codon_idx+2)):
            aa = aas[ci]
            current = sequence[ci*3:ci*3+3]
            for alt in sorted_codons[aa]:
                if alt != current:
                    test = sequence[:ci*3] + alt + sequence[ci*3+3:]
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

    # Phase 4: Fix 6+ consecutive T
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
            current = sequence[codon_idx*3:codon_idx*3+3]
            fixed = False
            for alt in sorted_codons[aa]:
                if alt != current:
                    test = sequence[:codon_idx*3] + alt + sequence[codon_idx*3+3:]
                    if not any(test[i:i+6] == "TTTTTT" for i in range(len(test)-5)):
                        sequence = test
                        fixed = True
                        break
            if not fixed:
                warnings.append(f"Consecutive T run: cannot fix at iteration {iteration}")
                break
    else:
        warnings.append("Consecutive T: max iterations reached, may still have 6+ T runs")

    # Phase 5: Adjust GC — use running count instead of recomputing each iteration
    gc_val = gc_content(sequence)
    target_gc = (gc_lo + gc_hi) / 2.0
    gc_count = sum(1 for b in sequence if b in "GC")
    n_bases = len(sequence)
    for iteration in range(200):
        if gc_lo <= gc_val <= gc_hi:
            break
        best_alt = None
        best_diff = abs(gc_val - target_gc)
        best_gc_delta = 0
        for ci in range(len(aas)):
            aa = aas[ci]
            current = sequence[ci*3:ci*3+3]
            current_gc = sum(1 for b in current if b in "GC")
            for alt in sorted_codons[aa]:
                if alt == current:
                    continue
                alt_gc = sum(1 for b in alt if b in "GC")
                new_gc_count = gc_count - current_gc + alt_gc
                new_frac = new_gc_count / n_bases
                diff = abs(new_frac - target_gc)
                if diff < best_diff:
                    best_diff = diff
                    best_alt = alt
                    best_gc_delta = alt_gc - current_gc
            if best_alt:
                # Apply immediately and update running count
                old_codon = sequence[ci*3:ci*3+3]
                sequence = sequence[:ci*3] + best_alt + sequence[ci*3+3:]
                gc_count += best_gc_delta
                gc_val = gc_count / n_bases
                best_alt = None  # Reset for next outer iteration
                break
    else:
        warnings.append(f"GC adjustment: max iterations reached, current GC={gc_val:.3f}")

    # Phase 6: Check cryptic splice sites
    max_d = max_donor_score(sequence)
    max_a = max_acceptor_score(sequence)
    if max_d >= cryptic_splice_threshold or max_a >= cryptic_splice_threshold:
        warnings.append(
            f"Cryptic splice sites remain: max_donor={max_d:.2f}, max_acceptor={max_a:.2f} "
            f"(threshold={cryptic_splice_threshold})"
        )

    for w in warnings:
        logger.warning(w)

    return sequence, warnings


def optimize_sequence(
    target_protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_sites: Optional[List[str]] = None,
    cai_threshold: float = 0.2,
    max_amino_acids_for_z3: int = 80,
    z3_timeout_ms: int = 30000,
    cryptic_splice_threshold: float = 3.0,
) -> OptimizationResult:
    """
    Find or generate a DNA sequence encoding the target protein that satisfies
    all type predicates, using z3 constraint solver with greedy fallback.

    The z3 solver uses a LOCAL solver timeout (not global set_param) for thread safety,
    and optimizes sum-of-logs (geometric mean proxy) for CAI.
    """
    aas = protein_to_aa_list(target_protein)
    n_aa = len(aas)
    if n_aa == 0:
        return OptimizationResult("", "", 0.0, 0.0, [], ["empty_sequence"])

    if organism not in SUPPORTED_ORGANISMS:
        raise UnsupportedOrganismError(organism, SUPPORTED_ORGANISMS)

    restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.values())
    fallback_used = False
    all_warnings: List[str] = []

    if n_aa > max_amino_acids_for_z3:
        sequence, warnings = _greedy_optimize(
            target_protein, organism, gc_lo, gc_hi, restriction_sites, cryptic_splice_threshold
        )
        fallback_used = True
        all_warnings.extend(warnings)
    else:
        try:
            from z3 import Int, Optimize, If, And, Or, Not, Implies, Sum, sat, Real, ToReal
        except ImportError:
            logger.warning("z3-solver not installed, falling back to greedy optimizer")
            sequence, warnings = _greedy_optimize(
                target_protein, organism, gc_lo, gc_hi, restriction_sites, cryptic_splice_threshold
            )
            fallback_used = True
            all_warnings.extend(warnings)
        else:
            # Use LOCAL solver timeout for thread safety — NOT global set_param
            opt = Optimize()
            opt.set("timeout", z3_timeout_ms)
            n_bases = n_aa * 3
            nuc = [Int(f"n_{i}") for i in range(n_bases)]

            for i in range(n_bases):
                opt.add(nuc[i] >= 0, nuc[i] <= 3)

            usage = CODON_ADAPTIVENESS_TABLES[organism]
            # Use sum-of-logs for CAI geometric mean proxy
            log_cai_contributions = []
            for aa_idx, aa in enumerate(aas):
                base_idx = aa_idx * 3
                n0, n1, n2 = nuc[base_idx], nuc[base_idx+1], nuc[base_idx+2]
                valid_codons = AA_TO_CODONS[aa]
                codon_constraints = []
                log_cai_chain = None
                sorted_z3 = sorted(valid_codons, key=lambda c: usage.get(c, 0.0), reverse=True)
                for codon in sorted_z3:
                    c0, c1, c2 = (BASE_REV[codon[0]], BASE_REV[codon[1]], BASE_REV[codon[2]])
                    match = And(n0 == c0, n1 == c1, n2 == c2)
                    codon_constraints.append(match)
                    w = usage.get(codon, 0.01)
                    # Use log(w) scaled by 1000 to avoid floating point in z3
                    # z3 works better with integers, so we approximate
                    log_w_int = int(math.log(max(w, 1e-10)) * 1000)
                    log_cai_chain = log_w_int if log_cai_chain is None else If(match, log_w_int, log_cai_chain)
                opt.add(Or(*codon_constraints))
                if log_cai_chain is not None:
                    log_cai_contributions.append(log_cai_chain)

            # Maximize sum of log(CAI) ≈ maximize geometric mean CAI
            if log_cai_contributions:
                opt.maximize(Sum(log_cai_contributions))

            # GCInRange
            gc_count = Sum([If(Or(nuc[i] == 1, nuc[i] == 2), 1, 0) for i in range(n_bases)])
            opt.add(gc_count >= int(gc_lo * n_bases), gc_count <= int(gc_hi * n_bases))

            # NoRestrictionSite
            for site in restriction_sites:
                site_upper = site.upper()
                indices = [BASE_REV[b] for b in site_upper]
                for start in range(n_bases - len(site_upper) + 1):
                    opt.add(Not(And(*[nuc[start+j] == indices[j] for j in range(len(site_upper))])))

            # NoInstabilityMotif
            for start in range(n_bases - 5 + 1):
                opt.add(Not(And(*[nuc[start+j] == BASE_REV[b] for j, b in enumerate("ATTTA")])))

            # No 6+ consecutive T
            for start in range(n_bases - 6 + 1):
                opt.add(Not(And(*[nuc[start+j] == 3 for j in range(6)])))

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

    cai = compute_cai(sequence, organism)
    gc = gc_content(sequence)
    satisfied, failed = _check_predicates(
        sequence, gc_lo, gc_hi, restriction_sites, cai_threshold, organism,
        cryptic_splice_threshold,
    )

    if all_warnings:
        logger.info("Optimization warnings: %s", "; ".join(all_warnings))

    return OptimizationResult(
        sequence=sequence, protein=target_protein, cai=cai,
        gc_content=gc, satisfied_predicates=satisfied,
        failed_predicates=failed, fallback_used=fallback_used,
    )


def _check_predicates(
    sequence: str, gc_lo: float, gc_hi: float,
    restriction_sites: List[str], cai_threshold: float, organism: str,
    cryptic_splice_threshold: float = 3.0,
) -> Tuple[List[str], List[str]]:
    """Check all type predicates against the optimized sequence."""
    satisfied, failed = [], []

    satisfied.append("InFrame")
    gc = gc_content(sequence)
    gc_ok = gc_lo <= gc <= gc_hi
    (satisfied if gc_ok else failed).append(
        "GCInRange" if gc_ok else f"GCInRange(gc={gc:.3f})"
    )

    has_restriction = any(
        site.upper() in sequence or reverse_complement(site.upper()) in sequence
        for site in restriction_sites
    )
    (satisfied if not has_restriction else failed).append("NoRestrictionSite")

    has_atta = "ATTTA" in sequence
    has_t6 = any(sequence[i:i+6] == "TTTTTT" for i in range(len(sequence)-5))
    inst_ok = not (has_atta or has_t6)
    (satisfied if inst_ok else failed).append("NoInstabilityMotif")

    cai = compute_cai(sequence, organism)
    cai_ok = cai >= cai_threshold
    (satisfied if cai_ok else failed).append(
        "CodonAdapted" if cai_ok else f"CodonAdapted(cai={cai:.4f})"
    )

    max_d = max_donor_score(sequence)
    max_a = max_acceptor_score(sequence)
    cryptic_ok = max_d < cryptic_splice_threshold and max_a < cryptic_splice_threshold
    (satisfied if cryptic_ok else failed).append(
        "NoCrypticSplice" if cryptic_ok
        else f"NoCrypticSplice(donor={max_d:.2f},acceptor={max_a:.2f})"
    )

    return satisfied, failed
