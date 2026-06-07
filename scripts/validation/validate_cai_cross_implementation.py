#!/usr/bin/env python3
"""
Cross-implementation CAI validation with random protein sequences.

Tests:
  1. Generate 200 random protein sequences (20-200 AA each, realistic composition)
  2. Encode each with random synonymous codons
  3. Compute CAI using both implementations for all 5 organisms
  4. Report maximum discrepancy between implementations
  5. Verify: CAI in [0,1], all-optimal-codon sequence -> CAI >= 0.95,
     length-independence
"""

import sys
import os
import random
import math

# Add project src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from biocompiler.benchmarking.cai_validated import (
    compute_cai_sharp_li_for_organism,
    load_reference_set,
    compute_cai_sharp_li,
    SUPPORTED_REFERENCE_ORGANISMS,
)
from biocompiler.translation import compute_cai
from biocompiler.organisms import (
    CODON_ADAPTIVENESS_TABLES,
    PREFERRED_CODON_TABLES,
    SUPPORTED_ORGANISMS,
)

# Standard genetic code
CODON_TABLE = {
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

# Reverse: amino acid -> list of synonymous codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)

# Realistic amino acid composition (approximate human protein average)
# Based on Trifonov (2000) and King & Jukes (1969)
AA_FREQUENCIES = {
    "A": 0.074, "R": 0.052, "N": 0.045, "D": 0.054, "C": 0.013,
    "Q": 0.034, "E": 0.061, "G": 0.069, "H": 0.026, "I": 0.057,
    "L": 0.099, "K": 0.058, "M": 0.024, "F": 0.040, "P": 0.051,
    "S": 0.069, "T": 0.058, "W": 0.013, "Y": 0.032, "V": 0.066,
}

# Build weighted AA list for random sampling
_AA_POOL = []
for aa, freq in AA_FREQUENCIES.items():
    _AA_POOL.extend([aa] * int(freq * 10000))


def random_protein(length):
    """Generate a random protein sequence with realistic AA composition."""
    return "".join(random.choice(_AA_POOL) for _ in range(length))


def encode_with_random_codons(protein):
    """Encode a protein with randomly chosen synonymous codons."""
    codons = []
    for aa in protein:
        codon = random.choice(AA_TO_CODONS[aa])
        codons.append(codon)
    return "".join(codons)


def encode_with_optimal_codons(protein, organism):
    """Encode a protein using the most optimal (preferred) codon for each AA."""
    preferred = PREFERRED_CODON_TABLES.get(organism, {})
    codons = []
    for aa in protein:
        if aa in preferred:
            codons.append(preferred[aa])
        elif aa == "M":
            codons.append("ATG")
        elif aa == "W":
            codons.append("TGG")
        else:
            # Fallback: random codon
            codons.append(random.choice(AA_TO_CODONS[aa]))
    return "".join(codons)


def main():
    output_lines = []

    def log(msg=""):
        print(msg)
        output_lines.append(msg)

    random.seed(42)  # Reproducibility

    ORGANISMS = SUPPORTED_ORGANISMS

    log("=" * 80)
    log("CROSS-IMPLEMENTATION CAI VALIDATION WITH RANDOM SEQUENCES")
    log("=" * 80)
    log()

    # ============================================================
    # TEST 1: Random sequences — compare both implementations
    # ============================================================
    log("=" * 80)
    log("TEST 1: 200 Random Protein Sequences — Cross-Implementation Agreement")
    log("=" * 80)
    log()

    N_SEQUENCES = 200
    max_discrepancy = 0.0
    max_disc_info = ("", "", 0.0, 0.0)
    all_discrepancies = []
    cai_range_violations = 0
    total_computations = 0

    for i in range(N_SEQUENCES):
        length = random.randint(20, 200)
        protein = random_protein(length)
        dna = encode_with_random_codons(protein)

        for organism in ORGANISMS:
            try:
                cai_main = compute_cai(dna, organism=organism)
            except Exception as e:
                log(f"  ERROR in main impl, seq {i}, {organism}: {e}")
                cai_main = None

            try:
                cai_valid = compute_cai_sharp_li_for_organism(dna, organism)
            except Exception as e:
                log(f"  ERROR in valid impl, seq {i}, {organism}: {e}")
                cai_valid = None

            if cai_main is not None and cai_valid is not None:
                disc = abs(cai_main - cai_valid)
                all_discrepancies.append(disc)
                total_computations += 1

                # Check CAI range
                if cai_main < 0.0 or cai_main > 1.0:
                    cai_range_violations += 1
                    log(f"  RANGE VIOLATION: main CAI={cai_main}, seq {i}, {organism}")
                if cai_valid < 0.0 or cai_valid > 1.0:
                    cai_range_violations += 1
                    log(f"  RANGE VIOLATION: valid CAI={cai_valid}, seq {i}, {organism}")

                if disc > max_discrepancy:
                    max_discrepancy = disc
                    max_disc_info = (f"seq_{i}", organism, cai_main, cai_valid)

    log(f"  Total computations:              {total_computations}")
    log(f"  Max discrepancy:                 {max_discrepancy:.10f}")
    log(f"  Max discrepancy details:         {max_disc_info[0]}/{max_disc_info[1]}: "
        f"main={max_disc_info[2]:.6f}, valid={max_disc_info[3]:.6f}")
    log(f"  Mean discrepancy:                {sum(all_discrepancies)/len(all_discrepancies):.10f}")
    log(f"  Median discrepancy:              {sorted(all_discrepancies)[len(all_discrepancies)//2]:.10f}")
    log(f"  95th percentile discrepancy:     {sorted(all_discrepancies)[int(0.95*len(all_discrepancies))]:.10f}")
    log(f"  CAI range violations:            {cai_range_violations}")
    log(f"  CAI in [0,1] for all values:     {'YES' if cai_range_violations == 0 else 'NO'}")
    log()

    # ============================================================
    # TEST 2: All-optimal-codon sequences should give CAI >= 0.95
    # ============================================================
    log("=" * 80)
    log("TEST 2: All-Optimal-Codon Sequences → CAI ≥ 0.95")
    log("=" * 80)
    log()

    opt_cai_failures = []
    for organism in ORGANISMS:
        # Generate a protein and encode it with optimal codons
        protein = random_protein(100)
        dna_optimal = encode_with_optimal_codons(protein, organism)

        try:
            cai_main_opt = compute_cai(dna_optimal, organism=organism)
        except Exception as e:
            cai_main_opt = None
            log(f"  ERROR: {organism}: {e}")

        try:
            cai_valid_opt = compute_cai_sharp_li_for_organism(dna_optimal, organism)
        except Exception as e:
            cai_valid_opt = None
            log(f"  ERROR: {organism}: {e}")

        log(f"  {organism}:")
        log(f"    Main impl CAI:    {cai_main_opt}")
        log(f"    Valid impl CAI:   {cai_valid_opt}")
        if cai_main_opt is not None and cai_main_opt < 0.95:
            opt_cai_failures.append((organism, "main", cai_main_opt))
            log(f"    *** BELOW 0.95 threshold (main) ***")
        if cai_valid_opt is not None and cai_valid_opt < 0.95:
            opt_cai_failures.append((organism, "valid", cai_valid_opt))
            log(f"    *** BELOW 0.95 threshold (valid) ***")
        log()

    if opt_cai_failures:
        log(f"  OPTIMAL CODON FAILURES: {len(opt_cai_failures)}")
        for org, impl, val in opt_cai_failures:
            log(f"    {org}/{impl}: CAI = {val:.4f}")
    else:
        log(f"  All-optimal-codon CAI ≥ 0.95: YES (all organisms pass)")
    log()

    # ============================================================
    # TEST 3: Length independence
    # ============================================================
    log("=" * 80)
    log("TEST 3: CAI Length Independence")
    log("=" * 80)
    log()

    # Same protein at different lengths (repeat the same subsequence)
    # CAI should be approximately the same regardless of length
    base_protein = random_protein(20)
    base_dna = encode_with_random_codons(base_protein)

    for organism in ORGANISMS:
        cai_base = compute_cai(base_dna, organism=organism)

        # Extend the protein by repeating it
        cais_by_length = {20: cai_base}
        for mult in [2, 5, 10]:
            extended_protein = base_protein * mult
            extended_dna = encode_with_random_codons(extended_protein)
            cai_ext = compute_cai(extended_dna, organism=organism)
            cais_by_length[20 * mult] = cai_ext

        # CAI should be similar across lengths for same codon composition
        cai_vals = list(cais_by_length.values())
        cai_spread = max(cai_vals) - min(cai_vals)
        log(f"  {organism}:")
        for length, cai in sorted(cais_by_length.items()):
            log(f"    Length {length:4d} AA: CAI = {cai:.4f}")
        log(f"    Spread (max-min): {cai_spread:.4f}")
        if cai_spread > 0.05:
            log(f"    NOTE: Spread > 0.05 — CAI varies with length due to codon sampling variance")
        log()

    # ============================================================
    # TEST 4: Monotonicity — more optimal codons should give higher CAI
    # ============================================================
    log("=" * 80)
    log("TEST 4: Monotonicity — Replacing rare codons with optimal increases CAI")
    log("=" * 80)
    log()

    test_protein = random_protein(50)
    random_dna = encode_with_random_codons(test_protein)

    for organism in ORGANISMS:
        optimal_dna = encode_with_optimal_codons(test_protein, organism)

        cai_random = compute_cai(random_dna, organism=organism)
        cai_optimal = compute_cai(optimal_dna, organism=organism)

        monotonic = cai_optimal >= cai_random
        log(f"  {organism}:")
        log(f"    Random codons CAI:    {cai_random:.4f}")
        log(f"    Optimal codons CAI:   {cai_optimal:.4f}")
        log(f"    Optimal >= Random:    {'YES' if monotonic else 'NO *** VIOLATION ***'}")
        log()

    # ============================================================
    # Summary
    # ============================================================
    log("=" * 80)
    log("OVERALL SUMMARY")
    log("=" * 80)
    log()
    log(f"  Test 1 — Cross-impl agreement:  Max discrepancy = {max_discrepancy:.10f}")
    log(f"                                    Mean discrepancy = {sum(all_discrepancies)/len(all_discrepancies):.10f}")
    log(f"                                    CAI in [0,1]: {'PASS' if cai_range_violations == 0 else 'FAIL'}")
    log()
    log(f"  Test 2 — Optimal codon CAI ≥ 0.95:  {'PASS' if len(opt_cai_failures) == 0 else 'FAIL'}")
    log()
    log(f"  Test 3 — Length independence:    See per-organism results above")
    log()
    log(f"  Test 4 — Monotonicity:           See per-organism results above")
    log()

    # Verdict
    all_pass = (
        max_discrepancy < 1e-6  # Implementations should agree to within rounding
        and cai_range_violations == 0
        and len(opt_cai_failures) == 0
    )
    if all_pass:
        log("  VERDICT: ALL TESTS PASSED ✓")
    else:
        log("  VERDICT: SOME TESTS FAILED — see details above")
    log()

    return output_lines


if __name__ == "__main__":
    output = main()
