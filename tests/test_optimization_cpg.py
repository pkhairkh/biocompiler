"""
Unit tests for biocompiler.optimization_cpg — CpG dinucleotide disruption and reconciliation.

Covers:
- disrupt_cpg_dinucleotides: CpG disruption with splice score constraints
- reconcile_cpg_sites: CpG disruption with restriction site constraints
"""

from __future__ import annotations

import pytest

from biocompiler.optimization_cpg import (
    disrupt_cpg_dinucleotides,
    reconcile_cpg_sites,
    MAX_CPG_DISRUPTION_ITERATIONS,
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


def _count_cpg(seq: str) -> int:
    return sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. disrupt_cpg_dinucleotides
# ═══════════════════════════════════════════════════════════════════════════════

class TestDisruptCpgDinucleotides:

    def test_no_cpg_returns_unchanged(self):
        """Sequence with no CG dinucleotides should be unchanged."""
        seq = "ATGATTAATTAATTAA"  # no CG
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = disrupt_cpg_dinucleotides(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
        )
        assert new_seq == seq
        assert len(warnings) == 0

    def test_preserves_protein(self):
        """CpG disruption should not change the protein."""
        seq = "ATGGCGAAGCTGGCCTGA"  # contains CG at positions 3, 9, 12
        aas = _seq_to_aas(seq)
        original_protein = "".join(aas)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = disrupt_cpg_dinucleotides(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
        )
        new_aas = _seq_to_aas(new_seq)
        new_protein = "".join(new_aas)
        assert new_protein == original_protein

    def test_preserves_sequence_length(self):
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = disrupt_cpg_dinucleotides(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
        )
        assert len(new_seq) == len(seq)

    def test_cpg_count_decreases_or_same(self):
        """After disruption, CpG count should not increase."""
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        initial_cpg = _count_cpg(seq)

        new_seq, warnings = disrupt_cpg_dinucleotides(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
        )
        new_cpg = _count_cpg(new_seq)
        assert new_cpg <= initial_cpg

    def test_returns_warnings_list(self):
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        _, warnings = disrupt_cpg_dinucleotides(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
        )
        assert isinstance(warnings, list)

    def test_high_splice_threshold_allows_more_disruption(self):
        """With a very high splice threshold, CpG disruption should be easier."""
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        initial_cpg = _count_cpg(seq)

        new_seq, warnings = disrupt_cpg_dinucleotides(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=99.0,  # very high — no splice constraint
        )
        new_cpg = _count_cpg(new_seq)
        assert new_cpg <= initial_cpg


# ═══════════════════════════════════════════════════════════════════════════════
# 2. reconcile_cpg_sites
# ═══════════════════════════════════════════════════════════════════════════════

class TestReconcileCpgSites:

    def test_no_cpg_returns_unchanged(self):
        seq = "ATGATTAATTAATTAA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = reconcile_cpg_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        assert new_seq == seq

    def test_preserves_protein(self):
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        original_protein = "".join(aas)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = reconcile_cpg_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        new_aas = _seq_to_aas(new_seq)
        assert "".join(new_aas) == original_protein

    def test_preserves_sequence_length(self):
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = reconcile_cpg_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        assert len(new_seq) == len(seq)

    def test_does_not_reintroduce_restriction_sites(self):
        """Reconciliation should not reintroduce restriction sites."""
        seq = "ATGGCGAAGCTGGCCTGA"  # no EcoRI site
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()

        new_seq, warnings = reconcile_cpg_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=99.0,
            concrete_sites=["GAATTC"],
        )
        assert "GAATTC" not in new_seq

    def test_cpg_count_does_not_increase(self):
        seq = "ATGGCGAAGCTGGCCTGA"
        aas = _seq_to_aas(seq)
        sorted_codons = _sorted_codons_for_aas(aas)
        usage = _uniform_usage()
        initial_cpg = _count_cpg(seq)

        new_seq, warnings = reconcile_cpg_sites(
            sequence=seq,
            aas=aas,
            sorted_codons=sorted_codons,
            usage=usage,
            cryptic_splice_threshold=3.0,
            concrete_sites=[],
        )
        assert _count_cpg(new_seq) <= initial_cpg
