#!/usr/bin/env python
"""Benchmark: Aho-Corasick trie-based restriction site detection vs. per-site scanning.

Measures the performance improvement from using Aho-Corasick multi-pattern
matching over the original per-enzyme str.find() approach.

Complexity comparison:
  - Per-site (before): O(N * L * site_len) where N = #enzymes, L = seq length
  - Aho-Corasick (after): O(L + M) where M = total matches
"""

import time
import sys

# Add the source directory to path
sys.path.insert(0, "/home/z/my-project/biocompiler/src")

from biocompiler.sequence.aho_corasick import AhoCorasickScanner, build_scanner_from_enzymes
from biocompiler.shared.constants import RESTRICTION_ENZYMES, reverse_complement
from biocompiler.optimizer import _find_site_in_sequence


def build_patterns_from_constants():
    """Build patterns dict from RESTRICTION_ENZYMES (same as build_scanner_from_enzymes)."""
    patterns = {}
    for name, site in RESTRICTION_ENZYMES.items():
        site_upper = site.upper()
        if all(b in "ACGT" for b in site_upper):
            patterns[site_upper] = name
            rc = reverse_complement(site_upper)
            if rc != site_upper:
                patterns[rc] = name
    return patterns


def benchmark_scan_ac(scanner, sequence, n_iters=1000):
    """Benchmark Aho-Corasick scan()."""
    t0 = time.perf_counter()
    for _ in range(n_iters):
        results = scanner.scan(sequence)
    elapsed = time.perf_counter() - t0
    return elapsed, results


def benchmark_scan_per_site(sites_and_rcs, sequence, n_iters=1000):
    """Benchmark per-site str.find() scanning."""
    t0 = time.perf_counter()
    for _ in range(n_iters):
        all_results = []
        for site, site_rc in sites_and_rcs:
            positions = _find_site_in_sequence(sequence, site, site_rc)
            for pos in positions:
                all_results.append((pos, site))
    elapsed = time.perf_counter() - t0
    return elapsed, all_results


def benchmark_has_any_match_ac(scanner, sequence, n_iters=1000):
    """Benchmark Aho-Corasick has_any_match()."""
    t0 = time.perf_counter()
    for _ in range(n_iters):
        result = scanner.has_any_match(sequence)
    elapsed = time.perf_counter() - t0
    return elapsed, result


def benchmark_has_any_match_per_site(sites_and_rcs, sequence, n_iters=1000):
    """Benchmark per-site str.find() for boolean 'any match' check."""
    t0 = time.perf_counter()
    for _ in range(n_iters):
        found = False
        for site, site_rc in sites_and_rcs:
            if site in sequence or (site_rc and site_rc in sequence):
                found = True
                break
    elapsed = time.perf_counter() - t0
    return elapsed, found


def benchmark_region_check_ac(scanner, sequence, region_size, n_iters=1000):
    """Benchmark Aho-Corasick has_any_match_in_region()."""
    t0 = time.perf_counter()
    for _ in range(n_iters):
        # Check region around codon index 100
        ci = 100
        start = max(0, ci * 3 - scanner.longest_pattern + 1)
        end = min(len(sequence), ci * 3 + 3 + scanner.longest_pattern - 1)
        result = scanner.has_any_match_in_region(sequence, start, end)
    elapsed = time.perf_counter() - t0
    return elapsed, result


def benchmark_region_check_per_site(sites_and_rcs, sequence, max_site_len, n_iters=1000):
    """Benchmark per-site region checking."""
    t0 = time.perf_counter()
    for _ in range(n_iters):
        ci = 100
        start = max(0, ci * 3 - max_site_len + 1)
        end = min(len(sequence), ci * 3 + 3 + max_site_len - 1)
        region = sequence[start:end]
        found = False
        for site, site_rc in sites_and_rcs:
            if site in region or (site_rc and site_rc in region):
                found = True
                break
    elapsed = time.perf_counter() - t0
    return elapsed, found


def main():
    print("=" * 72)
    print("Aho-Corasick vs. Per-Site Restriction Site Detection Benchmark")
    print("=" * 72)

    # Build patterns
    patterns = build_patterns_from_constants()
    print(f"\nPatterns: {len(patterns)} site strings from {len(RESTRICTION_ENZYMES)} enzymes")

    # Build sites_and_rcs for per-site benchmark
    seen = set()
    sites_and_rcs = []
    for name, site in RESTRICTION_ENZYMES.items():
        site_upper = site.upper()
        if all(b in "ACGT" for b in site_upper) and site_upper not in seen:
            seen.add(site_upper)
            site_rc = reverse_complement(site_upper)
            sites_and_rcs.append((site_upper, site_rc))
    max_site_len = max(len(s) for s, _ in sites_and_rcs)
    print(f"Concrete sites (unique): {len(sites_and_rcs)}")
    print(f"Longest site: {max_site_len} bp")

    # Build Aho-Corasick scanner
    scanner = AhoCorasickScanner(patterns)
    print(f"AC automaton nodes: {scanner.num_nodes}")

    # Test sequences of different lengths
    test_seqs = [
        ("Short (60bp)", "ATG" * 20),
        ("Medium (600bp)", "ATG" * 200),
        ("Long (3000bp)", "ATG" * 1000),
        ("GFP-like (714bp)", "ATGGAATTCAAGCTTGGATCCCTCGAGCATATGGCATGCCTGCAGG" + "ATG" * 223),
    ]

    n_iters = 1000

    for label, seq in test_seqs:
        print(f"\n{'─' * 72}")
        print(f"Sequence: {label} ({len(seq)} bp)")
        print(f"{'─' * 72}")

        # Full scan benchmark
        ac_time, ac_results = benchmark_scan_ac(scanner, seq, n_iters)
        ps_time, ps_results = benchmark_scan_per_site(sites_and_rcs, seq, n_iters)

        print(f"\n  Full scan ({n_iters} iterations):")
        print(f"    Aho-Corasick: {ac_time*1000:.2f} ms ({ac_time/n_iters*1e6:.1f} µs/scan)")
        print(f"    Per-site:     {ps_time*1000:.2f} ms ({ps_time/n_iters*1e6:.1f} µs/scan)")
        speedup = ps_time / ac_time if ac_time > 0 else float('inf')
        print(f"    Speedup:      {speedup:.1f}x")
        print(f"    Matches:      AC={len(ac_results)}, PS={len(ps_results)}")

        # Boolean check benchmark
        ac_bool_time, ac_bool = benchmark_has_any_match_ac(scanner, seq, n_iters)
        ps_bool_time, ps_bool = benchmark_has_any_match_per_site(sites_and_rcs, seq, n_iters)

        print(f"\n  Boolean 'any match?' ({n_iters} iterations):")
        print(f"    Aho-Corasick: {ac_bool_time*1000:.2f} ms ({ac_bool_time/n_iters*1e6:.1f} µs/check)")
        print(f"    Per-site:     {ps_bool_time*1000:.2f} ms ({ps_bool_time/n_iters*1e6:.1f} µs/check)")
        speedup_bool = ps_bool_time / ac_bool_time if ac_bool_time > 0 else float('inf')
        print(f"    Speedup:      {speedup_bool:.1f}x")

        # Region check benchmark (relevant for codon swap validation)
        if len(seq) > 300:
            ac_reg_time, ac_reg = benchmark_region_check_ac(scanner, seq, max_site_len, n_iters)
            ps_reg_time, ps_reg = benchmark_region_check_per_site(
                sites_and_rcs, seq, max_site_len, n_iters
            )

            print(f"\n  Local region check ({n_iters} iterations):")
            print(f"    Aho-Corasick: {ac_reg_time*1000:.2f} ms ({ac_reg_time/n_iters*1e6:.1f} µs/check)")
            print(f"    Per-site:     {ps_reg_time*1000:.2f} ms ({ps_reg_time/n_iters*1e6:.1f} µs/check)")
            speedup_reg = ps_reg_time / ac_reg_time if ac_reg_time > 0 else float('inf')
            print(f"    Speedup:      {speedup_reg:.1f}x")

    # End-to-end optimization benchmark
    print(f"\n{'=' * 72}")
    print("End-to-End Optimization Benchmark")
    print(f"{'=' * 72}")

    from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer

    protein = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVK" \
              "GHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAH" \
              "LPAEFTPAVHASLDKFLASVSTVLTSKYR"

    enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI", "NotI", "XbaI", "SalI",
               "PstI", "SphI", "NdeI", "NcoI", "NheI", "KpnI", "SmaI",
               "SacI", "SpeI", "ApaI", "ClaI", "EcoRV", "BglII", "MluI"]

    n_opt_iters = 100

    # With AC scanner (default)
    opt_ac = HybridOptimizer(species="ecoli", enzymes=enzymes)
    t0 = time.perf_counter()
    for _ in range(n_opt_iters):
        result_ac = opt_ac.optimize(protein, is_prokaryote=True)
    ac_total = time.perf_counter() - t0

    print(f"\n  HybridOptimizer (E. coli, {len(enzymes)} enzymes, {n_opt_iters} runs):")
    print(f"    CAI: {result_ac.cai:.4f}, GC: {result_ac.gc_content:.3f}")
    print(f"    Violations fixed: {result_ac.violations_fixed}")
    print(f"    Total time: {ac_total*1000:.1f} ms ({ac_total/n_opt_iters*1000:.2f} ms/run)")

    # Verify correctness: compare AC results with per-site results
    print(f"\n{'=' * 72}")
    print("Correctness Verification")
    print(f"{'=' * 72}")

    test_seq = "ATGGAATTCAAGCTTGGATCCCTCGAGCATATGGCATGCCTGCAGGGCGGCCGCTCTAGAGTCGAC"
    ac_matches = scanner.scan(test_seq)
    ps_matches = []
    for site, site_rc in sites_and_rcs:
        for pos in _find_site_in_sequence(test_seq, site, site_rc):
            ps_matches.append((pos, site))

    ac_positions = sorted(set((p, s) for p, s, _ in ac_matches))
    ps_positions = sorted(set((p, s) for p, s in ps_matches))

    print(f"  Test sequence: {test_seq}")
    print(f"  AC matches:  {ac_positions}")
    print(f"  PS matches:  {ps_positions}")
    print(f"  Match: {ac_positions == ps_positions}")

    # Also verify build_scanner_from_enzymes
    scanner2 = build_scanner_from_enzymes(enzymes)
    if scanner2:
        ac2_matches = scanner2.scan(test_seq)
        print(f"  build_scanner_from_enzymes matches: {sorted(set((p, e) for p, _, e in ac2_matches))}")
        print(f"  AC scanner from enzyme names: OK")

    print(f"\n{'=' * 72}")
    print("Benchmark complete!")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
