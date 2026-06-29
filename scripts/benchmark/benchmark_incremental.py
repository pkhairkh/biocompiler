#!/usr/bin/env python3
"""Benchmark: Old full-scan RS checking vs new incremental RS tracking.

Measures the speedup from using state.has_any_restriction_site() instead
of 'site in state.sequence' after each codon swap in the optimizer loop.

This is the key optimization: in the optimizer, every codon swap is followed
by a check for restriction sites. The old approach scans the entire sequence
(O(N*S)), while the new approach uses pre-tracked positions (O(S)).

Usage:
    cd /home/z/my-project/biocompiler
    python benchmark_incremental.py
"""
import time
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from biocompiler.optimizer.incremental import IncrementalSequenceState, EnzymeSiteCache
from biocompiler.sequence.restriction_sites import get_recognition_site
from biocompiler.shared.constants import reverse_complement


def get_rs_sites(enzymes):
    """Build the (site, site_rc) list like HybridOptimizer does."""
    sites = []
    for enz in enzymes:
        site = get_recognition_site(enz)
        if site is not None:
            site_rc = reverse_complement(site)
            sites.append((site, site_rc))
    return sites


def benchmark_rs_checking():
    """Benchmark the key optimization: RS checking after codon swaps."""
    enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'XbaI',
               'SalI', 'PstI', 'SphI', 'KpnI', 'SacI']
    rs_sites = get_rs_sites(enzymes)
    
    print("=" * 70)
    print("RS Check Benchmark: Full-Scan vs Incremental After Codon Swap")
    print("=" * 70)
    
    for seq_len in [600, 1800, 6000]:
        seq = 'ATG' * (seq_len // 3)
        n_codons = seq_len // 3
        
        print(f"\nSequence: {seq_len}bp ({n_codons} codons)")
        print("-" * 50)
        
        # ---- OLD APPROACH: full-scan RS check after each swap ----
        state_old = IncrementalSequenceState(seq, species='ecoli', enzymes=enzymes)
        
        # Simulate the optimizer's inner loop: swap codon, check RS
        t0 = time.perf_counter()
        for ci in range(0, min(n_codons, 200), 2):
            old = state_old.swap_codon(ci, 'GGC')
            # OLD: check RS by scanning entire sequence
            site_present = False
            for site, site_rc in rs_sites:
                if site in state_old.sequence or (site_rc and site_rc in state_old.sequence):
                    site_present = True
                    break
            state_old.swap_codon(ci, old)
        old_time = (time.perf_counter() - t0) / min(n_codons // 2, 100)
        
        # ---- NEW APPROACH: incremental RS check after each swap ----
        state_new = IncrementalSequenceState(seq, species='ecoli', enzymes=enzymes)
        
        t0 = time.perf_counter()
        for ci in range(0, min(n_codons, 200), 2):
            old = state_new.swap_codon(ci, 'GGC')
            # NEW: check RS using pre-tracked positions (O(S))
            site_present = state_new.has_any_restriction_site()
            state_new.swap_codon(ci, old)
        new_time = (time.perf_counter() - t0) / min(n_codons // 2, 100)
        
        speedup = old_time / new_time if new_time > 0 else float('inf')
        print(f"  Old (full scan):    {old_time*1000:.4f}ms per swap+check")
        print(f"  New (incremental):  {new_time*1000:.4f}ms per swap+check")
        print(f"  Speedup:            {speedup:.1f}x")
    
    # ---- GC checking benchmark ----
    print("\n" + "=" * 70)
    print("GC Check Benchmark: O(N) recount vs O(1) incremental")
    print("=" * 70)
    
    for seq_len in [600, 1800, 6000]:
        seq = 'ATG' * (seq_len // 3)
        state = IncrementalSequenceState(seq, species='ecoli')
        
        # OLD: sum over all bases
        t0 = time.perf_counter()
        for _ in range(10000):
            gc = sum(1 for b in state.sequence if b in 'GC')
            gc_frac = gc / len(state.sequence)
        old_gc_time = (time.perf_counter() - t0) / 10000
        
        # NEW: incremental gc_fraction
        t0 = time.perf_counter()
        for _ in range(10000):
            gc_frac = state.gc_fraction
        new_gc_time = (time.perf_counter() - t0) / 10000
        
        print(f"\n  {seq_len}bp:")
        print(f"    Old (O(N) sum):   {old_gc_time*1000000:.1f}µs")
        print(f"    New (O(1) incr):  {new_gc_time*1000000:.1f}µs")
        print(f"    Speedup:          {old_gc_time/new_gc_time:.0f}x")
    
    # ---- Splice site checking benchmark ----
    print("\n" + "=" * 70)
    print("Splice Site Benchmark: Full scan vs incremental (human species)")
    print("=" * 70)
    
    for seq_len in [600, 1800]:
        # Use a sequence with some GT/AG dinucleotides for splice scoring
        seq = 'ATG' * (seq_len // 3)
        state = IncrementalSequenceState(seq, species='human', enzymes=[])
        
        # Full splice check
        t0 = time.perf_counter()
        for _ in range(10):
            result = state.check_splice_sites(changed_only=False, threshold=3.0)
        full_splice_time = (time.perf_counter() - t0) / 10
        
        # Incremental splice check (after changing a codon)
        state.update_codon(50, 'GGT')  # GGT introduces a GT dinucleotide
        t0 = time.perf_counter()
        for _ in range(10):
            result = state.check_splice_sites(changed_only=True, threshold=3.0)
        incr_splice_time = (time.perf_counter() - t0) / 10
        
        speedup = full_splice_time / incr_splice_time if incr_splice_time > 0 else float('inf')
        print(f"\n  {seq_len}bp:")
        print(f"    Full scan:        {full_splice_time*1000:.3f}ms")
        print(f"    Incremental:      {incr_splice_time*1000:.4f}ms")
        print(f"    Speedup:          {speedup:.0f}x")


if __name__ == '__main__':
    benchmark_rs_checking()
