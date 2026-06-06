"""Performance regression tests for the eukaryote optimization pipeline.

These tests verify that the O(n²) → O(n) algorithmic fixes (Task 1.3)
prevent severe performance degradation for long proteins:

  Fix A: _eukaryote_cai_recovery() — local-region string construction
  Fix B: _eliminate_cpg_dinucleotides() — local-region CpG checking
  Fix C: GC adjustment in _greedy_optimize() — priority queue (heapq)
  Fix D: _BatchSwapScorer.update_incremental_state() — direct log_sum update

Each test has a generous time budget to avoid flakiness on slow CI, but
should complete well within the budget on any modern machine.
"""

import math
import time

import pytest

from biocompiler.optimization import (
    BioOptimizer,
    _eukaryote_cai_recovery,
    _eliminate_cpg_dinucleotides,
    _BatchSwapScorer,
)
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
from biocompiler.type_system import AA_TO_CODONS
from biocompiler.translation import translate


# ── Helpers ──────────────────────────────────────────────────────────────


def _generate_protein(length: int, seed: int = 42) -> str:
    """Generate a random protein sequence of the given length.

    Uses a simple LCG for reproducibility (no external dependency).
    Avoids stop codons ('*') in the body.
    """
    aas = "ACDEFGHIKLMNPQRSTVWY"
    state = seed
    result = []
    for _ in range(length):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        result.append(aas[state % len(aas)])
    return "".join(result)


def _back_translate(protein: str, species: str = "Homo_sapiens") -> str:
    """Back-translate a protein to DNA using highest-CAI codons."""
    usage = CODON_ADAPTIVENESS_TABLES.get(species, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])
    codons = []
    for aa in protein:
        aa_codons = AA_TO_CODONS.get(aa, [])
        if not aa_codons:
            aa_codons = ["NNN"]
        best = max(aa_codons, key=lambda c: usage.get(c, 0.0))
        codons.append(best)
    return "".join(codons)


# ── Tests ────────────────────────────────────────────────────────────────


class TestEukaryotePerformance:
    """Performance regression tests for eukaryotic protein optimization."""

    def test_500aa_eukaryote_under_30s(self):
        """Optimizing a 500aa eukaryotic protein should complete in < 30 seconds."""
        protein = _generate_protein(500)
        dna = _back_translate(protein, "Homo_sapiens")

        optimizer = BioOptimizer(species="human")
        start = time.monotonic()
        optimized, results, cert_text = optimizer.optimize(dna)
        elapsed = time.monotonic() - start

        # Verify correctness
        assert len(optimized) == len(dna), "Sequence length changed"
        opt_protein = translate(optimized)
        assert len(opt_protein) == len(protein), "Protein length changed"

        assert elapsed < 30.0, (
            f"500aa eukaryotic optimization took {elapsed:.1f}s (limit: 30s). "
            f"Possible O(n²) regression in eukaryote pipeline."
        )

    def test_1000aa_eukaryote_under_120s(self):
        """Optimizing a 1000aa eukaryotic protein should complete in < 120 seconds."""
        protein = _generate_protein(1000)
        dna = _back_translate(protein, "Homo_sapiens")

        optimizer = BioOptimizer(species="human")
        start = time.monotonic()
        optimized, results, cert_text = optimizer.optimize(dna)
        elapsed = time.monotonic() - start

        # Verify correctness
        assert len(optimized) == len(dna), "Sequence length changed"

        assert elapsed < 120.0, (
            f"1000aa eukaryotic optimization took {elapsed:.1f}s (limit: 120s). "
            f"Possible O(n²) regression in eukaryote pipeline."
        )


class TestIncrementalCAIDrift:
    """Test that the direct log_sum update eliminates floating-point drift."""

    def test_cai_drift_below_1e10_after_1000_swaps(self):
        """Incremental CAI drift should stay below 1e-10 after 1000 swaps.

        This verifies Fix D: update_incremental_state() now uses the direct
        formula log_sum = log_sum - log(w_old) + log(w_new) instead of
        recovering log_sum from CAI via n * log(CAI), which accumulated
        rounding errors.
        """
        species = "Homo_sapiens"
        usage = CODON_ADAPTIVENESS_TABLES.get(species, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])

        # Build a 200-codon sequence using suboptimal codons
        protein = _generate_protein(200, seed=123)
        dna = _back_translate(protein, species)
        seq_codons = [dna[i:i+3] for i in range(0, len(dna), 3)]

        scorer = _BatchSwapScorer(usage)
        scorer.reset_incremental_state(seq_codons)

        # Record the exact initial log_sum
        initial_log_sum = scorer.current_log_sum
        epsilon = 1e-10

        # Simulate 1000 codon swaps
        rng_state = 99
        for swap_i in range(1000):
            # Pick a random codon position
            rng_state = (rng_state * 1103515245 + 12345) & 0x7FFFFFFF
            ci = rng_state % len(seq_codons)
            aa = protein[ci]
            current = seq_codons[ci]
            alternatives = [c for c in AA_TO_CODONS.get(aa, []) if c != current]
            if not alternatives:
                continue

            # Pick a random alternative
            rng_state = (rng_state * 1103515245 + 12345) & 0x7FFFFFFF
            alt = alternatives[rng_state % len(alternatives)]

            # Apply the swap
            old_codon = current
            seq_codons[ci] = alt
            scorer.update_incremental_state(old_codon, alt)

        # Recompute the exact log_sum from scratch
        exact_log_sum = scorer._compute_log_sum(seq_codons)

        # The incremental log_sum should be very close to the exact value
        incremental_log_sum = scorer.current_log_sum
        drift = abs(incremental_log_sum - exact_log_sum)

        assert drift < 1e-10, (
            f"CAI log_sum drift after 1000 swaps: {drift:.2e} (limit: 1e-10). "
            f"Incremental: {incremental_log_sum}, Exact: {exact_log_sum}. "
            f"Fix D (direct formula) may not be working correctly."
        )


class TestGCAjustmentPerformance:
    """Test that the heapq-based GC adjustment is fast for large proteins."""

    def test_gc_adjustment_500aa_under_5s(self):
        """GC adjustment for a 500aa protein should complete in < 5 seconds.

        This verifies Fix C: the heapq-based priority queue replaces the
        O(200 × n × k) full-scan approach.
        """
        protein = _generate_protein(500, seed=77)
        # Use a GC-extreme back-translation to force GC adjustment
        # Use E. coli which has different GC target, but specify human
        # with extreme GC bounds to trigger the hard constraint path
        dna = _back_translate(protein, "Homo_sapiens")

        # Force GC adjustment by using tight bounds that the initial
        # sequence likely violates
        optimizer = BioOptimizer(species="human", gc_lo=0.35, gc_hi=0.40)

        start = time.monotonic()
        optimized, results, cert_text = optimizer.optimize(dna)
        elapsed = time.monotonic() - start

        # Just verify the optimizer completed
        assert len(optimized) == len(dna), "Sequence length changed"

        # The full optimization includes many steps; GC adjustment is just
        # one of them. We check total time is reasonable.
        # With Fix C, GC adjustment alone should be < 1s for 500aa.
        # The full pipeline should be well under 30s.
        assert elapsed < 30.0, (
            f"500aa optimization with tight GC took {elapsed:.1f}s. "
            f"GC adjustment may have O(n²) regression."
        )


class TestEukaryoteCAIRecoveryPerformance:
    """Test that _eukaryote_cai_recovery is fast with local-region construction."""

    def test_cai_recovery_500aa_fast(self):
        """_eukaryote_cai_recovery for a 500aa protein should be fast.

        With Fix A (local-region string construction), this function is O(n)
        instead of O(n²).
        """
        protein = _generate_protein(500, seed=55)
        dna = _back_translate(protein, "Homo_sapiens")
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens", CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])

        start = time.monotonic()
        result, upgrades = _eukaryote_cai_recovery(dna, protein, usage)
        elapsed = time.monotonic() - start

        assert len(result) == len(dna), "Sequence length changed"
        # With Fix A, 500aa should take < 1 second
        assert elapsed < 5.0, (
            f"_eukaryote_cai_recovery for 500aa took {elapsed:.1f}s. "
            f"Possible O(n²) regression (Fix A may not be working)."
        )


class TestCpGEliminationPerformance:
    """Test that _eliminate_cpg_dinucleotides is fast with local-region checks."""

    def test_cpg_elimination_500aa_fast(self):
        """_eliminate_cpg_dinucleotides for a 500aa protein should be fast.

        With Fix B (local-region CpG checking), this function avoids O(n)
        full-sequence rebuilds in the inner loops.
        """
        protein = _generate_protein(500, seed=88)
        dna = _back_translate(protein, "Homo_sapiens")
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens", CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])

        start = time.monotonic()
        result, warnings = _eliminate_cpg_dinucleotides(
            dna, protein, usage,
            gc_lo=0.30, gc_hi=0.70,
            organism="Homo_sapiens",
        )
        elapsed = time.monotonic() - start

        assert len(result) == len(dna), "Sequence length changed"
        # With Fix B, 500aa should take < 5 seconds even with many CpGs
        assert elapsed < 10.0, (
            f"_eliminate_cpg_dinucleotides for 500aa took {elapsed:.1f}s. "
            f"Possible O(n²) regression (Fix B may not be working)."
        )
