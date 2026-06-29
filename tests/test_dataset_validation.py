"""
BioCompiler Dataset Validation Tests

Tests the optimizer and core algorithms against well-known biological
sequences from common datasets. These tests validate that:

1. Translation fidelity is preserved during optimization
2. GC content falls within organism-specific expected ranges
3. CAI scores are within expected bounds for optimized sequences
4. Cross-organism optimization produces consistent results
5. Optimization actually improves over random codon assignment
6. Short peptides and long proteins are both handled correctly

Dataset sources:
- Human: TP53, BRCA1, CFTR, VEGFA, MYC, HBB (UniProt reference sequences)
- E. coli: LacZ, GFP, bla/ampR, recA, rpoB (standard cloning/molecular markers)
- Yeast: GAL4, ADH1, PGK1 (highly expressed S. cerevisiae genes)
- Synthetic: BH3 domain, WW domain, zinc finger, insulin B chain, IL-2, EGFP
"""

import os
import sys
import pytest

# Ensure src is on path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from biocompiler.validation.dataset_validation import (
    HUMAN_REFERENCE_GENES,
    ECOLI_REFERENCE_GENES,
    YEAST_REFERENCE_GENES,
    SYNTHETIC_BENCHMARKS,
    ALL_DATASETS,
    validate_translation_fidelity,
    validate_gc_content,
    validate_cai_bounds,
    validate_cross_organism_consistency,
    validate_protein_length,
    validate_optimization_improvement,
    validate_no_cpg_island,
    run_dataset_validation,
)
from biocompiler.expression.translation import translate, compute_cai
from biocompiler.sequence.scanner import gc_content
from biocompiler.optimizer import optimize_sequence
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS


# ============================================================================
# Dataset Completeness Tests — ensure data is not garbage
# ============================================================================

class TestDatasetCompleteness:
    """Validate that the reference datasets are complete and well-formed."""

    def test_human_genes_present(self):
        """Human reference set should have at least 5 genes."""
        assert len(HUMAN_REFERENCE_GENES) >= 5

    def test_ecoli_genes_present(self):
        """E. coli reference set should have at least 4 genes."""
        assert len(ECOLI_REFERENCE_GENES) >= 4

    def test_yeast_genes_present(self):
        """Yeast reference set should have at least 2 genes."""
        assert len(YEAST_REFERENCE_GENES) >= 2

    def test_synthetic_benchmarks_present(self):
        """Synthetic benchmark set should have at least 3 entries."""
        assert len(SYNTHETIC_BENCHMARKS) >= 3

    def test_all_genes_have_required_fields(self):
        """Every gene entry should have required fields."""
        required = {"description", "organism", "protein", "expected_gc_range", "protein_length"}
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                missing = required - set(gene_data.keys())
                assert not missing, (
                    f"Gene {ds_name}/{gene_name} missing fields: {missing}"
                )

    def test_all_proteins_are_valid_amino_acids(self):
        """Every protein sequence should contain only valid amino acid codes."""
        valid_aas = set(AA_TO_CODONS.keys())
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                protein = gene_data["protein"]
                invalid = set(protein) - valid_aas
                assert not invalid, (
                    f"Gene {ds_name}/{gene_name} has invalid amino acids: {invalid}"
                )

    def test_all_protein_lengths_match(self):
        """The protein_length field should match the actual protein length."""
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                expected = gene_data["protein_length"]
                actual = len(gene_data["protein"])
                assert actual == expected, (
                    f"Gene {ds_name}/{gene_name}: expected length {expected}, got {actual}"
                )

    def test_gc_ranges_are_valid(self):
        """GC ranges should be valid (lo < hi, both in [0, 1])."""
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                lo, hi = gene_data["expected_gc_range"]
                assert 0.0 <= lo < hi <= 1.0, (
                    f"Gene {ds_name}/{gene_name}: invalid GC range [{lo}, {hi}]"
                )

    def test_organisms_are_supported(self):
        """All organisms in the datasets should be supported by BioCompiler."""
        from biocompiler.organisms import SUPPORTED_ORGANISMS
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                organism = gene_data["organism"]
                assert organism in SUPPORTED_ORGANISMS, (
                    f"Gene {ds_name}/{gene_name}: unsupported organism '{organism}'"
                )


# ============================================================================
# Human Gene Validation Tests
# ============================================================================

class TestHumanGeneOptimization:
    """Validate optimization against human reference genes."""

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_translation_fidelity(self, gene_name):
        """Optimized human gene sequences must encode the correct protein."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        result = validate_translation_fidelity(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="human",
        )
        assert result.passed, (
            f"Translation fidelity failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_gc_content(self, gene_name):
        """Optimized human gene sequences should have GC in expected range."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        result = validate_gc_content(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="human",
            expected_gc_range=gene["expected_gc_range"],
        )
        assert result.passed, (
            f"GC content failed for {gene_name}: {result.actual} (expected {result.expected})"
        )

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_cai_bounds(self, gene_name):
        """Optimized human gene sequences should have CAI in expected range."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        cai_range = gene.get("expected_cai_human", (0.2, 1.0))
        result = validate_cai_bounds(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="human",
            expected_cai_range=cai_range,
        )
        assert result.passed, (
            f"CAI bounds failed for {gene_name}: {result.actual} (expected {result.expected})"
        )

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_protein_length(self, gene_name):
        """Optimized human gene sequences should have the correct number of codons."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        result = validate_protein_length(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="human",
            expected_length=gene["protein_length"],
        )
        assert result.passed, (
            f"Protein length failed for {gene_name}: {result.actual} (expected {result.expected})"
        )


# ============================================================================
# E. coli Gene Validation Tests
# ============================================================================

class TestEcoliGeneOptimization:
    """Validate optimization against E. coli reference genes."""

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_translation_fidelity(self, gene_name):
        """Optimized E. coli gene sequences must encode the correct protein."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        result = validate_translation_fidelity(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="ecoli",
        )
        assert result.passed, (
            f"Translation fidelity failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_gc_content(self, gene_name):
        """Optimized E. coli gene sequences should have GC in expected range."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        result = validate_gc_content(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="ecoli",
            expected_gc_range=gene["expected_gc_range"],
        )
        assert result.passed, (
            f"GC content failed for {gene_name}: {result.actual} (expected {result.expected})"
        )

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_cai_bounds(self, gene_name):
        """Optimized E. coli gene sequences should have CAI in expected range."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        cai_range = gene.get("expected_cai_ecoli", (0.2, 1.0))
        result = validate_cai_bounds(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="ecoli",
            expected_cai_range=cai_range,
        )
        assert result.passed, (
            f"CAI bounds failed for {gene_name}: {result.actual} (expected {result.expected})"
        )

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_protein_length(self, gene_name):
        """Optimized E. coli gene sequences should have correct codon count."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        result = validate_protein_length(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="ecoli",
            expected_length=gene["protein_length"],
        )
        assert result.passed, (
            f"Protein length failed for {gene_name}: {result.actual} (expected {result.expected})"
        )


# ============================================================================
# Yeast Gene Validation Tests
# ============================================================================

class TestYeastGeneOptimization:
    """Validate optimization against S. cerevisiae reference genes."""

    @pytest.mark.parametrize("gene_name", list(YEAST_REFERENCE_GENES.keys()))
    def test_yeast_translation_fidelity(self, gene_name):
        """Optimized yeast gene sequences must encode the correct protein."""
        gene = YEAST_REFERENCE_GENES[gene_name]
        result = validate_translation_fidelity(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="yeast",
        )
        assert result.passed, (
            f"Translation fidelity failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(YEAST_REFERENCE_GENES.keys()))
    def test_yeast_cai_bounds(self, gene_name):
        """Optimized yeast gene sequences should have CAI in expected range."""
        gene = YEAST_REFERENCE_GENES[gene_name]
        cai_range = gene.get("expected_cai_yeast", (0.2, 1.0))
        result = validate_cai_bounds(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="yeast",
            expected_cai_range=cai_range,
        )
        assert result.passed, (
            f"CAI bounds failed for {gene_name}: {result.actual} (expected {result.expected})"
        )


# ============================================================================
# Synthetic Benchmark Validation Tests
# ============================================================================

class TestSyntheticBenchmarks:
    """Validate optimization against synthetic benchmark proteins."""

    @pytest.mark.parametrize("gene_name", list(SYNTHETIC_BENCHMARKS.keys()))
    def test_synthetic_translation_fidelity(self, gene_name):
        """Optimized synthetic proteins must encode the correct sequence."""
        gene = SYNTHETIC_BENCHMARKS[gene_name]
        result = validate_translation_fidelity(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="synthetic",
        )
        assert result.passed, (
            f"Translation fidelity failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(SYNTHETIC_BENCHMARKS.keys()))
    def test_synthetic_protein_length(self, gene_name):
        """Optimized synthetic proteins should have correct codon count."""
        gene = SYNTHETIC_BENCHMARKS[gene_name]
        result = validate_protein_length(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="synthetic",
            expected_length=gene["protein_length"],
        )
        assert result.passed, (
            f"Protein length failed for {gene_name}: {result.actual} (expected {result.expected})"
        )


# ============================================================================
# Cross-Organism Consistency Tests
# ============================================================================

class TestCrossOrganismConsistency:
    """Validate that organism-specific optimization actually works."""

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_gene_cross_organism(self, gene_name):
        """Human genes optimized for E. coli vs human should show CAI differences."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        result = validate_cross_organism_consistency(
            protein=gene["protein"],
            gene_name=gene_name,
            dataset_name="human",
        )
        assert result.passed, (
            f"Cross-organism consistency failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_gene_cross_organism(self, gene_name):
        """E. coli genes optimized for different organisms should show CAI differences."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        result = validate_cross_organism_consistency(
            protein=gene["protein"],
            gene_name=gene_name,
            dataset_name="ecoli",
        )
        assert result.passed, (
            f"Cross-organism consistency failed for {gene_name}: {result.actual}"
        )


# ============================================================================
# Optimization Improvement Tests
# ============================================================================

class TestOptimizationImprovement:
    """Validate that optimization actually improves over random codon assignment."""

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_optimization_improvement(self, gene_name):
        """Optimized human sequences should have higher CAI than random."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        result = validate_optimization_improvement(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="human",
        )
        assert result.passed, (
            f"Optimization improvement failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_optimization_improvement(self, gene_name):
        """Optimized E. coli sequences should have higher CAI than random."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        result = validate_optimization_improvement(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="ecoli",
        )
        assert result.passed, (
            f"Optimization improvement failed for {gene_name}: {result.actual}"
        )

    @pytest.mark.parametrize("gene_name", list(SYNTHETIC_BENCHMARKS.keys()))
    def test_synthetic_optimization_improvement(self, gene_name):
        """Optimized synthetic sequences should have higher CAI than random."""
        gene = SYNTHETIC_BENCHMARKS[gene_name]
        result = validate_optimization_improvement(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="synthetic",
        )
        assert result.passed, (
            f"Optimization improvement failed for {gene_name}: {result.actual}"
        )


# ============================================================================
# CAI Computation Correctness Tests
# ============================================================================

class TestCAIComputation:
    """Validate CAI computation against known properties."""

    def test_all_preferred_codons_gives_cai_1(self):
        """A sequence using only preferred codons should have CAI close to 1.0."""
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES, PREFERRED_CODON_TABLES
        for organism in ["Homo_sapiens", "Escherichia_coli", "Saccharomyces_cerevisiae"]:
            preferred = PREFERRED_CODON_TABLES[organism]
            adaptiveness = CODON_ADAPTIVENESS_TABLES[organism]
            # Build a sequence using only preferred codons
            seq_parts = []
            # Use a few amino acids
            for aa in ["A", "G", "L", "K", "E", "V"]:
                codon = preferred[aa]
                seq_parts.append(codon * 5)  # 5 copies of each
            seq = "ATG" + "".join(seq_parts) + "TAA"
            cai = compute_cai(seq, organism)
            assert cai >= 0.95, (
                f"All-preferred-codon sequence for {organism} should have CAI >= 0.95, got {cai:.4f}"
            )

    def test_cai_is_organism_dependent(self):
        """Same sequence should have different CAI values for different organisms."""
        # Use a diverse protein with many codon choices
        # ATG (Met), AAA (Lys), TTT (Phe), GGG (Gly), CCC (Pro) — mix of preferred and non-preferred
        seq = "ATGAAATTTGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTAA"
        cai_human = compute_cai(seq, "Homo_sapiens")
        cai_ecoli = compute_cai(seq, "Escherichia_coli")
        cai_yeast = compute_cai(seq, "Saccharomyces_cerevisiae")
        # At least one pair should differ
        max_diff = max(abs(cai_human - cai_ecoli),
                       abs(cai_human - cai_yeast),
                       abs(cai_ecoli - cai_yeast))
        assert max_diff > 0.01, (
            f"CAI should differ across organisms: human={cai_human}, ecoli={cai_ecoli}, yeast={cai_yeast}"
        )

    def test_cai_deterministic(self):
        """Same sequence + same organism should always give the same CAI."""
        seq = "ATGAAATTTGGGCCCTAA"
        cai1 = compute_cai(seq, "Homo_sapiens")
        cai2 = compute_cai(seq, "Homo_sapiens")
        assert cai1 == cai2, f"CAI should be deterministic: {cai1} != {cai2}"

    def test_random_codon_assignment_has_lower_cai(self):
        """Randomly assigned codons should have lower CAI than optimized."""
        import random
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        # Random codon assignment
        random.seed(42)  # Reproducible
        random_seq_parts = []
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, ["ATG"])
            random_seq_parts.append(random.choice(codons))
        random_seq = "".join(random_seq_parts)
        random_cai = compute_cai(random_seq, "Homo_sapiens")

        # Optimized
        result = optimize_sequence(
            target_protein=protein,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            cai_threshold=0.2,
            strict_mode=False,
        )

        assert result.cai > random_cai, (
            f"Optimized CAI ({result.cai:.4f}) should exceed random CAI ({random_cai:.4f})"
        )


# ============================================================================
# Specific Gene Correctness Tests
# ============================================================================

class TestSpecificGeneCorrectness:
    """Test specific well-known genes against expected properties."""

    def test_hbb_translation_correct(self):
        """HBB (hemoglobin beta) should translate to the known protein sequence."""
        hbb_cds = (
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAG"
            "TTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTG"
            "GGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGG"
            "TGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCAC"
            "TGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCC"
            "CATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTA"
            "ATGCCCTGGCCCACAAGTATCACTAA"
        )
        protein = translate(hbb_cds)
        assert protein.startswith("MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"), (
            f"HBB translation incorrect: got {protein[:35]}"
        )
        assert len(protein.rstrip("*")) == 147, (
            f"HBB protein length should be 147, got {len(protein.rstrip('*'))}"
        )

    def test_gfp_translation_correct(self):
        """GFP should translate to the known protein sequence."""
        gfp_cds = (
            "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
            "CACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGC"
            "ACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCG"
            "AGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGA"
            "ACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTT"
            "AATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAAC"
            "TATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGA"
            "CACAACATTGAAGATGGTTCTTTAAATCAAG"
        )
        protein = translate(gfp_cds)
        assert protein.startswith("MSKGEELFT"), (
            f"GFP translation incorrect: got {protein[:10]}"
        )

    def test_tp53_known_protein_sequence(self):
        """TP53 protein from dataset should match UniProt reference."""
        tp53 = HUMAN_REFERENCE_GENES["TP53"]["protein"]
        # First few amino acids of TP53 from UniProt P04637
        assert tp53.startswith("MEEPQSDPSV"), (
            f"TP53 protein sequence does not match UniProt: starts with {tp53[:10]}"
        )

    def test_ecoli_recA_high_cai(self):
        """RecA is a highly expressed E. coli gene — optimization should produce high CAI."""
        recA = ECOLI_REFERENCE_GENES["recA"]
        result = optimize_sequence(
            target_protein=recA["protein"],
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            cai_threshold=0.2,
            strict_mode=False,
        )
        # High-expression genes should optimize to high CAI (optimizer picks preferred codons)
        assert result.cai >= 0.5, (
            f"RecA (high-expression E. coli gene) should optimize to CAI >= 0.5, got {result.cai:.4f}"
        )

    def test_yeast_adh1_high_cai(self):
        """ADH1 is one of the most highly expressed yeast genes."""
        adh1 = YEAST_REFERENCE_GENES["ADH1"]
        result = optimize_sequence(
            target_protein=adh1["protein"],
            organism="Saccharomyces_cerevisiae",
            gc_lo=0.20,
            gc_hi=0.70,
            cai_threshold=0.2,
            strict_mode=False,
        )
        assert result.cai >= 0.5, (
            f"ADH1 (high-expression yeast gene) should optimize to CAI >= 0.5, got {result.cai:.4f}"
        )


# ============================================================================
# Edge Cases in Dataset Context
# ============================================================================

class TestDatasetEdgeCases:
    """Test edge cases that are specific to real biological data."""

    def test_short_peptide_optimizes(self):
        """Very short peptides (BH3 domain, 23 aa) should optimize correctly."""
        bh3 = SYNTHETIC_BENCHMARKS["BH3_domain"]
        result = optimize_sequence(
            target_protein=bh3["protein"],
            organism=bh3["organism"],
            gc_lo=0.20,
            gc_hi=0.80,
            cai_threshold=0.1,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == bh3["protein"], (
            f"BH3 domain: translation mismatch after optimization"
        )

    def test_long_protein_optimizes(self):
        """Long proteins (TP53, 263 aa) should optimize successfully."""
        tp53 = HUMAN_REFERENCE_GENES["TP53"]
        result = optimize_sequence(
            target_protein=tp53["protein"],
            organism=tp53["organism"],
            gc_lo=0.30,
            gc_hi=0.70,
            cai_threshold=0.2,
            strict_mode=False,
        )
        # Should produce a valid optimized sequence
        assert len(result.sequence) > 0, "TP53 optimization produced empty sequence"
        assert len(result.sequence) % 3 == 0, "TP53 sequence not divisible by 3"
        assert len(result.sequence) == len(tp53["protein"]) * 3, (
            "TP53 sequence length should match protein"
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == tp53["protein"], (
            f"TP53: translation mismatch after optimization"
        )

    def test_met_start_preserved(self):
        """Optimized sequences for proteins starting with Met should start with ATG."""
        # Only test proteins that start with M (Met) — not all synthetic benchmarks do
        for gene_name, gene in SYNTHETIC_BENCHMARKS.items():
            if not gene["protein"].startswith("M"):
                continue
            result = optimize_sequence(
                target_protein=gene["protein"],
                organism=gene["organism"],
                gc_lo=0.30,
                gc_hi=0.70,
                strict_mode=False,
            )
            assert result.sequence.startswith("ATG"), (
                f"{gene_name}: optimized sequence for Met-starting protein should start with ATG, "
                f"starts with {result.sequence[:3]}"
            )

    def test_protein_with_many_leucines(self):
        """Proteins rich in Leucine (6 codons) stress-test codon selection."""
        # Construct a protein that is 50% leucine
        protein = "M" + "L" * 40 + "K" * 10
        result = optimize_sequence(
            target_protein=protein,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            cai_threshold=0.2,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein, (
            f"Leucine-rich protein: translation mismatch"
        )

    def test_protein_with_rare_amino_acids(self):
        """Proteins with Trp (1 codon) and Met (1 codon) should work trivially."""
        protein = "MWWWMWWWK"  # Trp and Met each have only 1 codon
        result = optimize_sequence(
            target_protein=protein,
            organism="Homo_sapiens",
            gc_lo=0.20,
            gc_hi=0.80,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein


# ============================================================================
# NoCpGIsland Validation Tests
# ============================================================================

class TestNoCpGIsland:
    """Validate that optimized sequences avoid CpG islands where possible.

    CpG island avoidance is now a systematic optimization pass (not just
    best-effort). For prokaryotes (E. coli), CpG islands are biologically
    irrelevant so the check is automatically skipped. For eukaryotes,
    the optimizer systematically eliminates CG dinucleotides that
    contribute to CpG islands. GC-rich genes may still contain CG
    dinucleotides if no synonymous substitution can eliminate them
    without creating restriction sites or other violations.
    """

    @pytest.mark.parametrize("gene_name", list(HUMAN_REFERENCE_GENES.keys()))
    def test_human_no_cpg_island(self, gene_name):
        """Optimized human genes should avoid CpG islands where possible."""
        gene = HUMAN_REFERENCE_GENES[gene_name]
        result = validate_no_cpg_island(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="human",
        )
        assert result.passed, f"CpG island found in {gene_name}: {result.actual}"

    @pytest.mark.parametrize("gene_name", list(ECOLI_REFERENCE_GENES.keys()))
    def test_ecoli_no_cpg_island(self, gene_name):
        """Optimized E. coli genes pass CpG check (automatically skipped for prokaryotes)."""
        gene = ECOLI_REFERENCE_GENES[gene_name]
        result = validate_no_cpg_island(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="ecoli",
        )
        # CpG islands are biologically irrelevant for prokaryotes,
        # so the check is automatically skipped and always passes.
        assert result.passed, f"CpG island found in {gene_name}: {result.actual}"

    @pytest.mark.parametrize("gene_name", list(SYNTHETIC_BENCHMARKS.keys()))
    def test_synthetic_no_cpg_island(self, gene_name):
        """Optimized synthetic proteins should avoid CpG islands where possible."""
        gene = SYNTHETIC_BENCHMARKS[gene_name]
        result = validate_no_cpg_island(
            protein=gene["protein"],
            organism=gene["organism"],
            gene_name=gene_name,
            dataset_name="synthetic",
        )
        assert result.passed, f"CpG island found in {gene_name}: {result.actual}"

    def test_no_cpg_island_aggregate_pass_rate(self):
        """At least 80% of genes should avoid CpG islands."""
        report = run_dataset_validation(
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=True,
        )
        cpg_results = [r for r in report.results if r.test_type == "no_cpg_island"]
        assert len(cpg_results) > 0, "No CpG island tests were run"
        cpg_passed = sum(1 for r in cpg_results if r.passed)
        cpg_rate = cpg_passed / len(cpg_results)
        assert cpg_rate >= 0.80, (
            f"CpG island avoidance rate {cpg_rate:.1%} is below 80% threshold. "
            f"Passed: {cpg_passed}/{len(cpg_results)}"
        )


# ============================================================================
# Full Dataset Validation Runner
# ============================================================================

class TestFullDatasetValidation:
    """Run the complete dataset validation suite."""

    def test_full_validation_pass_rate(self):
        """Full validation should achieve at least 85% pass rate.

        NoCpGIsland results are excluded from the pass rate since they are
        informational / best-effort and should not gate the overall threshold.
        """
        report = run_dataset_validation(
            include_cross_organism=True,
            include_optimization_improvement=True,
            include_no_cpg_island=True,
        )
        # Exclude informational NoCpGIsland tests from the overall pass rate
        non_cpg = [r for r in report.results if r.test_type != "no_cpg_island"]
        non_cpg_passed = sum(1 for r in non_cpg if r.passed)
        non_cpg_rate = non_cpg_passed / max(len(non_cpg), 1)
        assert non_cpg_rate >= 0.85, (
            f"Full validation pass rate {non_cpg_rate:.1%} is below 85% threshold "
            f"(excluding NoCpGIsland). "
            f"Failed tests:\n" +
            "\n".join(
                f"  {r.dataset_name}/{r.gene_name}/{r.test_type}: {r.actual}"
                for r in non_cpg if not r.passed
            )
        )

    def test_all_translation_fidelity_passes(self):
        """ALL translation fidelity tests must pass — zero tolerance for wrong protein."""
        report = run_dataset_validation(
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        fidelity_results = [r for r in report.results if r.test_type == "translation_fidelity"]
        assert len(fidelity_results) > 0, "No translation fidelity tests were run"
        all_pass = all(r.passed for r in fidelity_results)
        assert all_pass, (
            f"Translation fidelity failures:\n" +
            "\n".join(
                f"  {r.gene_name}: {r.actual}"
                for r in fidelity_results if not r.passed
            )
        )

    def test_all_protein_lengths_correct(self):
        """ALL protein length tests must pass — zero tolerance for wrong length."""
        report = run_dataset_validation(
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        length_results = [r for r in report.results if r.test_type == "protein_length"]
        assert len(length_results) > 0, "No protein length tests were run"
        all_pass = all(r.passed for r in length_results)
        assert all_pass, (
            f"Protein length failures:\n" +
            "\n".join(
                f"  {r.gene_name}: {r.actual} (expected {r.expected})"
                for r in length_results if not r.passed
            )
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
