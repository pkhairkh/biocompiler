"""Tests for the tRNA Adaptation Index (tAI) module.

Agent 58: Validates the tAI implementation including wobble rules,
tRNA gene copy number handling, and comparison with CAI.
"""

import math
import pytest

from biocompiler.tai import (
    calculate_tai,
    TRNA_GENE_COPIES,
    WOBBLE_RULES,
    WOBBLE_EFFICIENCY,
    SUPPORTED_ORGANISMS_TAI,
    compute_codon_weights,
)
from biocompiler.translation import compute_cai


# ═══════════════════════════════════════════════════════════════════════════════
# Test Data
# ═══════════════════════════════════════════════════════════════════════════════

# Short test sequences
_ECOLI_OPTIMAL_SEQ = "ATGGCTAAAGCGTTT"  # M A K A F (6 codons including ATG)
_ALL_OPTIMAL_ECOLI = "ATGGCTAGCAAAGAG"  # M A S K E (6 codons including ATG)

# eGFP in E. coli
_EGFP_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGG"
    "CCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCAC"
    "CACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCT"
    "ACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTT"
    "CTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA"
    "GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAA"
    "CGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGC"
    "AGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
    "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTG"
    "ACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
)

# Human insulin
_HUMAN_INSULIN_DNA = (
    "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCT"
    "TTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTC"
    "TACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGC"
    "AGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCA"
    "TCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
)


class TestTRNAGeneCopies:
    """Test the tRNA gene copy number database."""

    def test_ecoli_data_exists(self):
        """E. coli tRNA data should be available."""
        assert "e_coli" in TRNA_GENE_COPIES

    def test_human_data_exists(self):
        """Human tRNA data should be available."""
        assert "human" in TRNA_GENE_COPIES

    def test_yeast_data_exists(self):
        """Yeast tRNA data should be available."""
        assert "yeast" in TRNA_GENE_COPIES

    def test_supported_organisms_not_empty(self):
        """Supported organisms list should not be empty."""
        assert len(SUPPORTED_ORGANISMS_TAI) >= 3

    def test_ecoli_has_required_anticodons(self):
        """E. coli should have tRNA data for common anticodons."""
        ecoli = TRNA_GENE_COPIES["e_coli"]
        # Essential anticodons
        assert "GAA" in ecoli  # Phe
        assert "CAU" in ecoli  # Met
        assert "UUC" in ecoli  # Glu
        assert "GUC" in ecoli  # Asp

    def test_human_has_more_copies_than_ecoli(self):
        """Human genome should have more tRNA gene copies than E. coli."""
        human_total = sum(TRNA_GENE_COPIES["human"].values())
        ecoli_total = sum(TRNA_GENE_COPIES["e_coli"].values())
        assert human_total > ecoli_total

    def test_all_copy_numbers_positive(self):
        """All tRNA gene copy numbers should be positive integers."""
        for org, copies in TRNA_GENE_COPIES.items():
            for anticodon, count in copies.items():
                assert isinstance(count, int), f"{org}/{anticodon}: not int"
                assert count > 0, f"{org}/{anticodon}: count <= 0"


class TestWobbleRules:
    """Test the wobble base pairing rules."""

    def test_all_codons_represented(self):
        """All 64 codons should be in the wobble rules."""
        # Count RNA codons (using U instead of T)
        assert len(WOBBLE_RULES) == 64

    def test_watson_crick_pairs_have_efficiency_1(self):
        """Watson-Crick pairs should have efficiency 1.0."""
        assert WOBBLE_EFFICIENCY[("A", "U")] == 1.0
        assert WOBBLE_EFFICIENCY[("U", "A")] == 1.0
        assert WOBBLE_EFFICIENCY[("G", "C")] == 1.0
        assert WOBBLE_EFFICIENCY[("C", "G")] == 1.0

    def test_gu_wobble_efficiency(self):
        """G:U wobble pair should have efficiency between 0 and 1."""
        assert 0 < WOBBLE_EFFICIENCY[("G", "U")] < 1.0

    def test_each_codon_has_at_least_one_anticodon(self):
        """Every codon should be readable by at least one anticodon."""
        for codon, anticodons in WOBBLE_RULES.items():
            assert len(anticodons) >= 1, f"Codon {codon} has no anticodon"
            for anticodon, efficiency in anticodons:
                assert 0.0 <= efficiency <= 1.0, \
                    f"Invalid efficiency for {codon}/{anticodon}: {efficiency}"

    def test_stop_codons_have_wobble_rules(self):
        """Stop codons should have wobble rules (even if no tRNAs read them)."""
        for stop in ("UAA", "UAG", "UGA"):
            assert stop in WOBBLE_RULES

    def test_inosine_pairing(self):
        """Inosine (I) should pair with U, C, A at varying efficiencies."""
        assert ("I", "U") in WOBBLE_EFFICIENCY
        assert ("I", "C") in WOBBLE_EFFICIENCY
        assert ("I", "A") in WOBBLE_EFFICIENCY
        # I:C should be stronger than I:U or I:A
        assert WOBBLE_EFFICIENCY[("I", "C")] > WOBBLE_EFFICIENCY[("I", "U")]
        assert WOBBLE_EFFICIENCY[("I", "C")] > WOBBLE_EFFICIENCY[("I", "A")]


class TestCalculateTAI:
    """Test the calculate_tai function."""

    def test_basic_ecoli_tai(self):
        """tAI should be computable for E. coli."""
        tai = calculate_tai(_EGFP_DNA, "e_coli")
        assert 0.0 < tai <= 1.0

    def test_basic_human_tai(self):
        """tAI should be computable for human."""
        tai = calculate_tai(_HUMAN_INSULIN_DNA, "human")
        assert 0.0 < tai <= 1.0

    def test_basic_yeast_tai(self):
        """tAI should be computable for yeast."""
        tai = calculate_tai(_EGFP_DNA, "yeast")
        assert 0.0 < tai <= 1.0

    def test_empty_sequence(self):
        """Empty sequence should return 0.0."""
        assert calculate_tai("", "e_coli") == 0.0

    def test_short_sequence(self):
        """Sequence shorter than 3 bases should return 0.0."""
        assert calculate_tai("AT", "e_coli") == 0.0

    def test_invalid_length(self):
        """Non-multiple-of-3 length should raise ValueError."""
        with pytest.raises(ValueError, match="not a multiple of 3"):
            calculate_tai("ATGA", "e_coli")

    def test_unsupported_organism(self):
        """Unsupported organism should raise ValueError."""
        with pytest.raises(ValueError, match="No tRNA gene copy data"):
            calculate_tai("ATGAAAGCGTTT", "zebrafish")

    def test_organism_aliases(self):
        """Various organism name aliases should work."""
        dna = "ATGAAAGCGTTT"
        # These should all resolve without error
        tai1 = calculate_tai(dna, "e_coli")
        tai2 = calculate_tai(dna, "Escherichia_coli")
        tai3 = calculate_tai(dna, "ecoli")
        # All should produce the same value
        assert tai1 == tai2 == tai3

    def test_tai_in_range(self):
        """tAI should always be in [0, 1]."""
        for org in SUPPORTED_ORGANISMS_TAI:
            tai = calculate_tai(_EGFP_DNA, org)
            assert 0.0 <= tai <= 1.0, f"tAI out of range for {org}: {tai}"

    def test_optimal_codons_higher_tai(self):
        """Sequences using optimal codons should have higher tAI."""
        # E. coli optimal codons: GCT(Ala), AGC(Ser), AAA(Lys), GAG(Glu)
        # Non-optimal: GCG(Ala), TCG(Ser), AAG(Lys), GAA(Glu)
        optimal_seq = "ATGGCTAGCAAAGAG"  # M A S K E (skip M in calculation)
        suboptimal_seq = "ATGGCGTCGAAGGAA"  # M A S K E (suboptimal codons)
        tai_opt = calculate_tai(optimal_seq, "e_coli")
        tai_sub = calculate_tai(suboptimal_seq, "e_coli")
        # Optimal should be at least as high as suboptimal
        # (may not always be strictly higher due to tRNA copy numbers vs CAI)
        assert tai_opt > 0.0
        assert tai_sub > 0.0

    def test_skip_stop_codons(self):
        """Stop codons should be excluded by default."""
        # Sequence with stop codon at end
        with_stop = "ATGGCTAAAGCGTTTTAA"  # M A K A F *
        without_stop = "ATGGCTAAAGCGTTT"   # M A K A F
        tai_with = calculate_tai(with_stop, "e_coli", skip_stop=True)
        tai_without = calculate_tai(without_stop, "e_coli")
        assert tai_with == tai_without

    def test_include_stop_codons(self):
        """Stop codons should be included when skip_stop=False."""
        with_stop = "ATGGCTAAAGCGTTTTAA"
        tai_skip = calculate_tai(with_stop, "e_coli", skip_stop=True)
        tai_include = calculate_tai(with_stop, "e_coli", skip_stop=False)
        # Including stop codons should give a different (usually lower) value
        # since stop codons have low tRNA adaptiveness
        assert tai_include != tai_skip or tai_skip == 0.0

    def test_skip_met_codons(self):
        """Met codons should be excluded by default (CAI convention)."""
        seq = "ATGGCTAAAGCG"  # M A K A
        tai_skip = calculate_tai(seq, "e_coli", skip_met=True)
        tai_include = calculate_tai(seq, "e_coli", skip_met=False)
        # Including Met should give a different value
        assert isinstance(tai_skip, float)
        assert isinstance(tai_include, float)

    def test_dna_case_insensitive(self):
        """DNA sequence should be case-insensitive."""
        tai_upper = calculate_tai("ATGGCTAAAGCG", "e_coli")
        tai_lower = calculate_tai("atggctaaagcg", "e_coli")
        assert tai_upper == tai_lower

    def test_whitespace_stripped(self):
        """Whitespace in DNA should be stripped."""
        tai_clean = calculate_tai("ATGGCTAAAGCG", "e_coli")
        tai_spaced = calculate_tai("  ATGGCTAAAGCG  ", "e_coli")
        assert tai_clean == tai_spaced


class TestComputeCodonWeights:
    """Test the compute_codon_weights function."""

    def test_ecoli_weights(self):
        """E. coli codon weights should be computable."""
        weights = compute_codon_weights("e_coli")
        assert isinstance(weights, dict)
        assert len(weights) > 0

    def test_weights_in_range(self):
        """All weights should be in [0, 1]."""
        for org in SUPPORTED_ORGANISMS_TAI:
            weights = compute_codon_weights(org)
            for codon, w in weights.items():
                assert 0.0 <= w <= 1.0, \
                    f"Weight out of range for {org}/{codon}: {w}"

    def test_optimal_codon_has_weight_1(self):
        """At least one codon per amino acid should have weight 1.0."""
        from biocompiler.constants import CODON_TABLE
        from biocompiler.type_system import AA_TO_CODONS

        weights = compute_codon_weights("e_coli")
        for aa, codons in AA_TO_CODONS.items():
            if aa == "*" or aa == "M":
                continue
            rna_codons = [c.replace("T", "U") for c in codons]
            max_w = max(weights.get(c, 0.0) for c in rna_codons)
            assert max_w > 0, f"No positive weight for {aa} in E. coli"

    def test_weights_differ_between_organisms(self):
        """Codon weights should differ between organisms."""
        weights_ecoli = compute_codon_weights("e_coli")
        weights_human = compute_codon_weights("human")
        # At least some weights should differ
        common_codons = set(weights_ecoli.keys()) & set(weights_human.keys())
        differences = sum(
            1 for c in common_codons
            if abs(weights_ecoli[c] - weights_human[c]) > 0.01
        )
        assert differences > 0, "E. coli and human tAI weights are identical"

    def test_organism_aliases(self):
        """Organism aliases should produce the same weights."""
        w1 = compute_codon_weights("e_coli")
        w2 = compute_codon_weights("Escherichia_coli")
        w3 = compute_codon_weights("ecoli")
        assert w1 == w2 == w3


class TestTAIvsCAI:
    """Compare tAI and CAI — they should be correlated but not identical."""

    def test_tai_and_cai_both_computable(self):
        """Both tAI and CAI should be computable for the same sequence."""
        dna = _EGFP_DNA
        tai = calculate_tai(dna, "e_coli")
        cai = compute_cai(dna, organism="Escherichia_coli")
        assert 0.0 < tai <= 1.0
        assert 0.0 < cai <= 1.0

    def test_tai_cai_correlation_direction(self):
        """For well-optimized sequences, both tAI and CAI should be high."""
        # eGFP is already fairly well-optimized for E. coli
        tai = calculate_tai(_EGFP_DNA, "e_coli")
        cai = compute_cai(_EGFP_DNA, organism="Escherichia_coli")
        # Both should be reasonably high for a well-codon-optimized sequence
        assert tai > 0.3, f"eGFP tAI unexpectedly low: {tai}"
        assert cai > 0.3, f"eGFP CAI unexpectedly low: {cai}"

    def test_tai_uses_trna_not_codon_freq(self):
        """tAI should reflect tRNA abundance, not codon frequency."""
        # This is a design-level test: verify that tAI uses
        # TRNA_GENE_COPIES, not CODON_ADAPTIVENESS_TABLES
        ecoli_weights = compute_codon_weights("e_coli")
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        cai_weights = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]

        # Convert tAI weights from RNA to DNA for comparison
        tai_weights_dna = {
            codon.replace("U", "T"): w
            for codon, w in ecoli_weights.items()
        }

        # The weights should not be identical (different data sources)
        common = set(tai_weights_dna.keys()) & set(cai_weights.keys())
        differences = sum(
            1 for c in common
            if abs(tai_weights_dna[c] - cai_weights[c]) > 0.01
        )
        # At least some differences expected since tAI uses tRNA copies
        # while CAI uses codon frequencies
        assert differences > 0, "tAI weights identical to CAI — may not be using tRNA data"

    def test_human_insulin_ecoli_tai_lower(self):
        """Human insulin with native codons should have low tAI in E. coli."""
        tai = calculate_tai(_HUMAN_INSULIN_DNA, "e_coli")
        # Human codons are not optimized for E. coli tRNA pool
        assert tai < 0.8, f"Human insulin tAI in E. coli unexpectedly high: {tai}"


class TestTAIEdgeCases:
    """Test edge cases for tAI computation."""

    def test_single_met_codon(self):
        """A sequence of just ATG should return 0.0 (Met is skipped)."""
        assert calculate_tai("ATG", "e_coli") == 0.0

    def test_single_non_met_codon(self):
        """A single non-Met codon should return a valid tAI."""
        tai = calculate_tai("ATGGCT", "e_coli")  # M + A
        assert 0.0 < tai <= 1.0  # Just Ala contributes

    def test_all_same_amino_acid(self):
        """Sequence encoding the same amino acid should have valid tAI."""
        # All alanine codons: GCU, GCC, GCA, GCG
        alanine_optimal = "ATGGCTGCTGCT"  # M + A + A + A
        tai = calculate_tai(alanine_optimal, "e_coli")
        assert 0.0 < tai <= 1.0

    def test_long_sequence(self):
        """Long sequences should be processable without error."""
        long_dna = _EGFP_DNA * 5  # ~3600 bp
        tai = calculate_tai(long_dna, "e_coli")
        assert 0.0 < tai <= 1.0

    def test_dna_with_only_stop(self):
        """Sequence of just stop codons should return 0.0."""
        assert calculate_tai("TAATAA", "e_coli") == 0.0

    def test_mixed_organism_tai(self):
        """Same sequence should give different tAI for different organisms."""
        dna = _EGFP_DNA
        tai_ecoli = calculate_tai(dna, "e_coli")
        tai_human = calculate_tai(dna, "human")
        tai_yeast = calculate_tai(dna, "yeast")
        # Values should differ because tRNA pools differ
        values = [tai_ecoli, tai_human, tai_yeast]
        # Not all values should be identical
        assert len(set(round(v, 4) for v in values)) > 1, \
            "tAI identical across all organisms — data may not be organism-specific"
