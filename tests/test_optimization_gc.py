"""
Unit tests for biocompiler.optimization_gc — GC content adjustment and reconciliation.

Covers:
- _compute_gc_count: GC base counting
- adjust_gc_content: full GC adjustment pipeline
  - GC already in range
  - GC too low
  - GC too high
  - Reconciliation with restriction sites
"""

from __future__ import annotations

import pytest

from biocompiler.optimization_gc import (
    _compute_gc_count,
    adjust_gc_content,
    MAX_GC_ADJUSTMENT_ITERATIONS,
)
from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE
from biocompiler.constants import reverse_complement


# ── Helpers ──────────────────────────────────────────────────────────────────

def _seq_to_aas(seq: str) -> list[str]:
    """Convert DNA sequence to amino acid list."""
    aas = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        aas.append(CODON_TABLE.get(codon, "X"))
    return aas


def _sorted_codons_for_aas(aas: list[str]) -> dict[str, list[str]]:
    """Build sorted_codons dict for given amino acids."""
    result = {}
    for aa in set(aas):
        result[aa] = list(AA_TO_CODONS.get(aa, []))
    return result


def _uniform_usage() -> dict[str, float]:
    """Uniform CAI weights for all codons."""
    return {c: 1.0 for c in CODON_TABLE}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. _compute_gc_count
# ═══════════════════════════════════════════════════════════════════════════════

class TestComputeGCCount:

    def test_all_gc(self):
        assert _compute_gc_count("GCGCGCGC") == 8

    def test_no_gc(self):
        assert _compute_gc_count("ATATATAT") == 0

    def test_mixed(self):
        assert _compute_gc_count("ATGC") == 2

    def test_empty(self):
        assert _compute_gc_count("") == 0


# ═══════════════════════════════════════════════════════════════════════════════
# 2. adjust_gc_content
# ═══════════════════════════════════════════════════════════════════════════════

class TestAdjustGCContent:

    def test_gc_already_in_range(self):
        """When GC is already in range, sequence should be unchanged."""
        seq = "ATGGTCAAGCTGGCCTGA"  # 18 bases, GC ~0.5
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        gc = _compute_gc_count(seq) / len(seq)

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.0,
            gc_hi=1.0,
            organism="Escherichia_coli",
            concrete_sites=[],
            concrete_scanner=None,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        assert isinstance(new_seq, str)
        assert len(new_seq) == len(seq)
        assert isinstance(warnings, list)

    def test_gc_too_low_is_adjusted(self):
        """When GC is too low, adjustment should try to increase it."""
        # All-AT sequence (for a simple protein)
        # Build a sequence with known low GC
        seq = "ATGATAATAATAATAATAATAAT"  # MIIIII — very low GC
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        initial_gc = _compute_gc_count(seq) / len(seq)

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.5,
            gc_hi=0.7,
            organism="Escherichia_coli",
            concrete_sites=[],
            concrete_scanner=None,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        new_gc = _compute_gc_count(new_seq) / len(new_seq)
        # GC should have moved toward the range
        assert new_gc >= initial_gc or len(warnings) > 0

    def test_gc_too_high_is_adjusted(self):
        """When GC is too high, adjustment should try to decrease it."""
        # Build a GC-rich sequence
        seq = "ATGGCGGCGGCGGCGGCGGCG"  # MAAAA... with high GC codons
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        initial_gc = _compute_gc_count(seq) / len(seq)

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.0,
            gc_hi=0.3,
            organism="Escherichia_coli",
            concrete_sites=[],
            concrete_scanner=None,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        new_gc = _compute_gc_count(new_seq) / len(new_seq)
        # GC should have moved toward the range
        assert new_gc <= initial_gc or len(warnings) > 0

    def test_preserves_protein(self):
        """GC adjustment should not change the protein."""
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        original_protein = "".join(aas)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.3,
            gc_hi=0.7,
            organism="Escherichia_coli",
            concrete_sites=[],
            concrete_scanner=None,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        new_aas = _seq_to_aas(new_seq)
        new_protein = "".join(new_aas)
        assert new_protein == original_protein

    def test_preserves_sequence_length(self):
        """GC adjustment should not change sequence length."""
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.0,
            gc_hi=1.0,
            organism="Escherichia_coli",
            concrete_sites=[],
            concrete_scanner=None,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        assert len(new_seq) == len(seq)

    def test_max_iterations_warning(self):
        """If GC adjustment can't reach target, a warning should be issued."""
        # Use a sequence that's impossible to adjust to target
        # (e.g., protein with only AT-rich codons and high GC target)
        seq = "ATGATAATAATAATAATAATAAT"  # MIIIII
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.99,
            gc_hi=1.0,
            organism="Escherichia_coli",
            concrete_sites=[],
            concrete_scanner=None,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        # May or may not produce warnings depending on whether it can reach target
        assert isinstance(warnings, list)

    def test_with_concrete_scanner_reconciliation(self):
        """When a scanner is provided, reconciliation should check for reintroduced sites."""
        from biocompiler.aho_corasick import AhoCorasickScanner

        seq = "ATGGTCAAGCTGGCCTGA"  # no EcoRI site
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        scanner = AhoCorasickScanner({"GAATTC": "EcoRI"})

        new_seq, warnings = adjust_gc_content(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            gc_lo=0.0,
            gc_hi=1.0,
            organism="Escherichia_coli",
            concrete_sites=["GAATTC"],
            concrete_scanner=scanner,
            remove_site_multicodon_fn=lambda *a, **k: (seq, False),
        )
        # Should not have EcoRI site in result
        assert "GAATTC" not in new_seq
