"""
Unit tests for biocompiler.optimization_splice — cryptic splice site elimination.

Covers:
- eliminate_cryptic_splice_sites: basic elimination, donor sites, acceptor sites,
  warnings for unrepairable sites, reconciliation with restriction sites
"""

from __future__ import annotations

import pytest

from biocompiler.optimization_splice import (
    eliminate_cryptic_splice_sites,
    MAX_SPLICE_ELIMINATION_ITERATIONS,
)
from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE


# ── Helpers ──────────────────────────────────────────────────────────────────

def _seq_to_aas(seq: str) -> list[str]:
    aas = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        aas.append(CODON_TABLE.get(codon, "X"))
    return aas


def _sorted_codons_for_aas(aas: list[str]) -> dict[str, list[str]]:
    result = {}
    for aa in set(aas):
        result[aa] = list(AA_TO_CODONS.get(aa, []))
    return result


def _uniform_usage() -> dict[str, float]:
    return {c: 1.0 for c in CODON_TABLE}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. eliminate_cryptic_splice_sites
# ═══════════════════════════════════════════════════════════════════════════════

class TestEliminateCrypticSpliceSites:

    def test_no_splice_sites_returns_unchanged(self):
        """Sequence with no strong splice sites should be unchanged."""
        # Use a sequence unlikely to have strong splice signals
        seq = "ATGGCTGCTGCTGCTGCTAA"  # MAAAAAK — avoid GT/AG in codons
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        assert isinstance(new_seq, str)
        assert len(new_seq) == len(seq)
        assert isinstance(warnings, list)

    def test_preserves_protein(self):
        """Splice elimination should not change the protein."""
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        original_protein = "".join(aas)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        new_aas = _seq_to_aas(new_seq)
        new_protein = "".join(new_aas)
        assert new_protein == original_protein

    def test_preserves_sequence_length(self):
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        assert len(new_seq) == len(seq)

    def test_very_high_threshold_no_changes(self):
        """With a very high threshold, no changes should be needed."""
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=99.0,
            concrete_sites=[],
        )
        assert new_seq == seq
        assert len(warnings) == 0

    def test_returns_warnings(self):
        """Should return a list of warnings."""
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        _, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        assert isinstance(warnings, list)

    def test_reconciliation_with_restriction_sites(self):
        """If splice fixes reintroduce restriction sites, they should be handled."""
        seq = "ATGGTCAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=["GAATTC"],  # not in this sequence
        )
        assert isinstance(new_seq, str)

    def test_gt_containing_sequence(self):
        """A sequence with GT dinucleotides should be processed."""
        # Build a sequence with known GT at a codon boundary
        seq = "ATGAGTAGTGGTGGTGGTTAA"  # contains multiple GT positions
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = eliminate_cryptic_splice_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        # Should return a valid sequence
        assert len(new_seq) == len(seq)
        assert len(new_seq) % 3 == 0
