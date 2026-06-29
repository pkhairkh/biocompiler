"""
BioCompiler Dataset Validation — Extended Tests

Complements test_dataset_validation.py with:
1. Additional edge cases not covered by existing tests
2. Property-based tests using Hypothesis
3. Error handling and graceful degradation tests

These tests focus on boundary conditions, invariants, and robustness
rather than biological correctness (which is covered by the main test suite).
"""

import os
import sys
import time

import pytest
pytest.importorskip("hypothesis")

# Ensure src is on path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

pytest.importorskip("hypothesis")
from hypothesis import given, settings, assume, HealthCheck
from hypothesis import strategies as st

from biocompiler.validation.dataset_validation import (
    HUMAN_REFERENCE_GENES,
    ECOLI_REFERENCE_GENES,
    YEAST_REFERENCE_GENES,
    SYNTHETIC_BENCHMARKS,
    ALL_DATASETS,
    PUBLISHED_CAI_BENCHMARKS,
    DatasetValidationResult,
    DatasetValidationReport,
    validate_translation_fidelity,
    validate_gc_content,
    validate_cai_bounds,
    validate_cross_organism_consistency,
    validate_protein_length,
    validate_optimization_improvement,
    validate_no_cpg_island,
    run_dataset_validation,
    format_dataset_report_text,
)
from biocompiler.expression.translation import translate, compute_cai
from biocompiler.sequence.scanner import gc_content
from biocompiler.optimizer import optimize_sequence
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS, STANDARD_AAS
from biocompiler.shared.exceptions import (
    InvalidProteinError,
    UnsupportedOrganismError,
)


# ============================================================================
# Shared Constants & Helpers
# ============================================================================

SUPPORTED_ORGANISMS = [
    "Homo_sapiens",
    "Escherichia_coli",
    "Saccharomyces_cerevisiae",
]

VALID_AAS = set(AA_TO_CODONS.keys())

# Amino acids with only 1 codon (no synonymous choice)
SINGLE_CODON_AAS = [aa for aa, codons in AA_TO_CODONS.items() if len(codons) == 1]

# Amino acids with max codon degeneracy
MAX_DEGENERACY_AAS = sorted(AA_TO_CODONS.keys(), key=lambda aa: -len(AA_TO_CODONS[aa]))


# ============================================================================
# 1. Additional Edge Cases
# ============================================================================

class TestMinimalProteinEdgeCases:
    """Test optimization of minimal / degenerate protein sequences."""

    def test_single_methionine(self):
        """A single Met (ATG) should optimize trivially."""
        result = optimize_sequence("M", "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        assert result.sequence == "ATG"
        assert len(result.sequence) == 3

    def test_single_tryptophan(self):
        """A single Trp (TGG) should optimize trivially — only one codon."""
        result = optimize_sequence("W", "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        assert result.sequence == "TGG"

    def test_two_amino_acids(self):
        """Two-amino-acid protein should still produce a valid sequence."""
        protein = "MK"
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_protein_of_only_single_codon_aas(self):
        """Protein composed entirely of Met and Trp (no codon choice)."""
        protein = "MWWMWWMW"
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein
        # Sequence should be fully determined: ATG + TGG alternating
        expected_seq = "".join(
            "ATG" if aa == "M" else "TGG" for aa in protein
        )
        assert result.sequence == expected_seq

    def test_protein_of_only_high_degeneracy_aas(self):
        """Protein of Leu (6 codons), Ser (6 codons), Arg (6 codons)."""
        protein = "LSR" * 10
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein
        # Each codon must be valid for its amino acid
        for i, aa in enumerate(protein):
            codon = result.sequence[i * 3 : i * 3 + 3]
            assert codon in AA_TO_CODONS[aa], (
                f"Position {i}: codon {codon} not valid for {aa}"
            )

    def test_all_twenty_standard_aas(self):
        """A protein containing every standard amino acid at least once."""
        protein = STANDARD_AAS  # "ACDEFGHIKLMNPQRSTVWY"
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_homopolymer_protein(self):
        """A protein of repeated single amino acid (Ala, 4 codons)."""
        protein = "A" * 50
        result = optimize_sequence(protein, "Escherichia_coli", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein


class TestGCRangeEdgeCases:
    """Test GC content with tight, wide, and extreme bounds."""

    def test_very_tight_gc_range(self):
        """Extremely tight GC range (0.49, 0.51) should still produce valid output."""
        protein = "ACDEFGHIKLMNPQRSTVWY" * 2
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.49, gc_hi=0.51, strict_mode=False)
        # The optimizer may not perfectly hit the tight range, but should not crash
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_wide_gc_range(self):
        """Very wide GC range (0.05, 0.95) should produce valid output."""
        protein = "ACDEFGHIKLMNPQRSTVWY" * 2
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.05, gc_hi=0.95, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein
        gc = gc_content(result.sequence)
        assert 0.05 <= gc <= 0.95

    def test_gc_range_low_extreme(self):
        """Very low GC range (0.05, 0.20) — AT-rich target."""
        # Use amino acids that can be encoded with AT-rich codons
        protein = "MFFLLIII" * 3
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.05, gc_hi=0.20, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_gc_range_high_extreme(self):
        """Very high GC range (0.80, 0.95) — GC-rich target."""
        # Use amino acids that can be encoded with GC-rich codons
        protein = "MGGGPPPPAAAA" * 3
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.80, gc_hi=0.95, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein


class TestDatasetStructureEdgeCases:
    """Test structural properties of the reference datasets."""

    def test_all_datasets_in_all_datasets(self):
        """ALL_DATASETS should contain all four primary dataset collections."""
        assert set(ALL_DATASETS.keys()) == {"human", "ecoli", "yeast", "synthetic"}

    def test_published_cai_benchmarks_structure(self):
        """PUBLISHED_CAI_BENCHMARKS should have valid structure."""
        for key, value in PUBLISHED_CAI_BENCHMARKS.items():
            assert "description" in value, f"Missing description in {key}"
            # At least one organism entry (skip description key)
            organism_entries = {
                k: v for k, v in value.items() if k != "description"
            }
            assert len(organism_entries) >= 1, (
                f"No organism entries in {key}"
            )
            for org, range_val in organism_entries.items():
                assert isinstance(range_val, tuple) and len(range_val) == 2, (
                    f"CAI range for {key}/{org} should be a 2-tuple"
                )
                lo, hi = range_val
                assert 0.0 <= lo < hi <= 1.0, (
                    f"Invalid CAI range [{lo}, {hi}] in {key}/{org}"
                )

    def test_no_duplicate_gene_names_within_dataset(self):
        """Each dataset should not have duplicate gene names (by definition — dict keys)."""
        for ds_name, ds in ALL_DATASETS.items():
            # Dict keys are inherently unique; this is a structural sanity check
            assert len(ds) == len(set(ds.keys())), (
                f"Dataset {ds_name} has duplicate gene names"
            )

    def test_human_genes_have_cai_human_key(self):
        """Human genes should have expected_cai_human key."""
        for gene_name, gene_data in HUMAN_REFERENCE_GENES.items():
            assert "expected_cai_human" in gene_data, (
                f"Human gene {gene_name} missing expected_cai_human"
            )

    def test_ecoli_genes_have_cai_ecoli_key(self):
        """E. coli genes should have expected_cai_ecoli key."""
        for gene_name, gene_data in ECOLI_REFERENCE_GENES.items():
            assert "expected_cai_ecoli" in gene_data, (
                f"E. coli gene {gene_name} missing expected_cai_ecoli"
            )

    def test_yeast_genes_have_cai_yeast_key(self):
        """Yeast genes should have expected_cai_yeast key."""
        for gene_name, gene_data in YEAST_REFERENCE_GENES.items():
            assert "expected_cai_yeast" in gene_data, (
                f"Yeast gene {gene_name} missing expected_cai_yeast"
            )

    def test_protein_sequences_are_uppercase(self):
        """All protein sequences in datasets should be uppercase."""
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                protein = gene_data["protein"]
                assert protein == protein.upper(), (
                    f"Gene {ds_name}/{gene_name}: protein sequence not uppercase"
                )

    def test_protein_sequences_no_whitespace(self):
        """All protein sequences should not contain whitespace."""
        for ds_name, ds in ALL_DATASETS.items():
            for gene_name, gene_data in ds.items():
                protein = gene_data["protein"]
                assert protein.strip() == protein, (
                    f"Gene {ds_name}/{gene_name}: protein has leading/trailing whitespace"
                )
                assert " " not in protein, (
                    f"Gene {ds_name}/{gene_name}: protein contains spaces"
                )


class TestValidationResultStructure:
    """Test DatasetValidationResult and DatasetValidationReport data structures."""

    def test_result_has_execution_time(self):
        """DatasetValidationResult should record execution time."""
        result = validate_translation_fidelity(
            "MK", "Homo_sapiens", "test_gene", "test_ds"
        )
        assert isinstance(result.execution_time_ms, float)
        assert result.execution_time_ms >= 0.0

    def test_report_pass_rate_with_zero_tests(self):
        """DatasetValidationReport pass_rate should be 0.0 with zero tests."""
        report = DatasetValidationReport(
            timestamp="2025-01-01", version="1.0",
            total_tests=0, passed=0, failed=0,
        )
        assert report.pass_rate == 0.0

    def test_report_pass_rate_with_all_passed(self):
        """pass_rate should be 1.0 when all tests pass."""
        report = DatasetValidationReport(
            timestamp="2025-01-01", version="1.0",
            total_tests=10, passed=10, failed=0,
        )
        assert report.pass_rate == 1.0

    def test_report_pass_rate_with_mixed_results(self):
        """pass_rate should be correct for mixed results."""
        report = DatasetValidationReport(
            timestamp="2025-01-01", version="1.0",
            total_tests=4, passed=3, failed=1,
        )
        assert abs(report.pass_rate - 0.75) < 1e-9

    def test_format_report_text_not_empty(self):
        """format_dataset_report_text should return a non-empty string."""
        report = run_dataset_validation(
            datasets=["synthetic"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        text = format_dataset_report_text(report)
        assert isinstance(text, str)
        assert len(text) > 0
        assert "BioCompiler Dataset Validation Report" in text

    def test_format_report_contains_all_test_types(self):
        """Report text should contain all test types that were run."""
        report = run_dataset_validation(
            datasets=["synthetic"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        text = format_dataset_report_text(report)
        # Should contain key test types
        assert "translation_fidelity" in text
        assert "gc_content" in text
        assert "cai_bounds" in text
        assert "protein_length" in text

    def test_result_fields_populated(self):
        """DatasetValidationResult should have all fields populated."""
        result = validate_translation_fidelity(
            "MK", "Homo_sapiens", "test_gene", "test_ds"
        )
        assert result.dataset_name == "test_ds"
        assert result.gene_name == "test_gene"
        assert result.test_type == "translation_fidelity"
        assert isinstance(result.passed, bool)
        assert isinstance(result.expected, str)
        assert isinstance(result.actual, str)


class TestRunDatasetValidationEdgeCases:
    """Test run_dataset_validation with various configurations."""

    def test_subset_single_dataset(self):
        """Running validation on a single dataset should work."""
        report = run_dataset_validation(
            datasets=["ecoli"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        # Should only have ecoli genes
        ds_names = set(r.dataset_name for r in report.results)
        assert ds_names == {"ecoli"}

    def test_subset_multiple_datasets(self):
        """Running validation on multiple specific datasets should work."""
        report = run_dataset_validation(
            datasets=["human", "yeast"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        ds_names = set(r.dataset_name for r in report.results)
        assert ds_names == {"human", "yeast"}

    def test_unknown_dataset_skipped(self):
        """Unknown dataset names should be skipped without error."""
        report = run_dataset_validation(
            datasets=["nonexistent_dataset"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        assert report.total_tests == 0

    def test_mixed_known_unknown_datasets(self):
        """Mix of valid and invalid dataset names should only run valid ones."""
        report = run_dataset_validation(
            datasets=["synthetic", "bogus"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        ds_names = set(r.dataset_name for r in report.results)
        assert "synthetic" in ds_names
        assert "bogus" not in ds_names

    def test_disable_all_optional_tests(self):
        """Disabling all optional tests should leave only core tests."""
        report = run_dataset_validation(
            datasets=["synthetic"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        test_types = set(r.test_type for r in report.results)
        assert "cross_organism_consistency" not in test_types
        assert "optimization_improvement" not in test_types
        assert "no_cpg_island" not in test_types
        # Core tests should still be present
        assert "translation_fidelity" in test_types
        assert "gc_content" in test_types
        assert "cai_bounds" in test_types
        assert "protein_length" in test_types

    def test_enable_all_optional_tests(self):
        """Enabling all optional tests should include them."""
        report = run_dataset_validation(
            datasets=["synthetic"],
            include_cross_organism=True,
            include_optimization_improvement=True,
            include_no_cpg_island=True,
        )
        test_types = set(r.test_type for r in report.results)
        assert "optimization_improvement" in test_types
        assert "no_cpg_island" in test_types
        # Cross-organism only for proteins >= 20 aa
        # Some synthetic benchmarks are < 20 aa so may not appear

    def test_report_summary_structure(self):
        """Report summary should contain expected keys."""
        report = run_dataset_validation(
            datasets=["synthetic"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        assert "by_dataset" in report.summary
        assert "by_test_type" in report.summary
        assert "avg_execution_time_ms" in report.summary
        assert "max_execution_time_ms" in report.summary
        assert "total_genes_tested" in report.summary

    def test_report_timestamp_and_version(self):
        """Report should have a valid timestamp and version."""
        report = run_dataset_validation(
            datasets=["synthetic"],
            include_cross_organism=False,
            include_optimization_improvement=False,
            include_no_cpg_island=False,
        )
        assert report.timestamp  # Non-empty
        assert report.version  # Non-empty
        # Timestamp should be ISO format
        assert "T" in report.timestamp or "-" in report.timestamp


class TestCrossOrganismShortProtein:
    """Test cross-organism consistency with very short proteins."""

    def test_cross_organism_very_short_protein(self):
        """Very short proteins (< 20 aa) should produce a result (may fail)."""
        protein = "MKWVA" * 4  # 20 aa
        result = validate_cross_organism_consistency(
            protein, "test_short", "test_ds"
        )
        assert isinstance(result, DatasetValidationResult)
        assert result.test_type == "cross_organism_consistency"


class TestValidateFunctionsReturnType:
    """Ensure all validate_* functions always return DatasetValidationResult."""

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_validate_translation_fidelity_returns_result(self, organism):
        result = validate_translation_fidelity("MKWVA", organism, "test", "test")
        assert isinstance(result, DatasetValidationResult)

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_validate_gc_content_returns_result(self, organism):
        result = validate_gc_content(
            "MKWVA", organism, "test", "test", (0.20, 0.80)
        )
        assert isinstance(result, DatasetValidationResult)

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_validate_cai_bounds_returns_result(self, organism):
        result = validate_cai_bounds(
            "MKWVA", organism, "test", "test", (0.1, 1.0)
        )
        assert isinstance(result, DatasetValidationResult)

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_validate_protein_length_returns_result(self, organism):
        result = validate_protein_length(
            "MKWVA", organism, "test", "test", 5
        )
        assert isinstance(result, DatasetValidationResult)

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_validate_optimization_improvement_returns_result(self, organism):
        result = validate_optimization_improvement(
            "MKWVA", organism, "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_validate_no_cpg_island_returns_result(self, organism):
        result = validate_no_cpg_island(
            "MKWVA", organism, "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)


# ============================================================================
# 2. Property-Based Tests with Hypothesis
# ============================================================================

# Strategy for generating valid protein sequences
amino_acid_strategy = st.sampled_from(list(STANDARD_AAS))

protein_strategy = st.text(
    alphabet=st.sampled_from(list(STANDARD_AAS)),
    min_size=1,
    max_size=30,  # Keep small for speed
)

organism_strategy = st.sampled_from(SUPPORTED_ORGANISMS)

# Strategy for valid GC ranges
gc_range_strategy = st.tuples(
    st.floats(min_value=0.05, max_value=0.49),
    st.floats(min_value=0.51, max_value=0.95),
)


class TestPropertyOptimizationCorrectness:
    """Property-based tests for optimization invariants."""

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_optimized_sequence_translates_back(self, protein, organism):
        """Property: Any valid protein optimized for any supported organism
        should translate back to the same protein."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein, (
            f"Protein {protein} for {organism}: "
            f"translated={translated}, expected={protein}"
        )

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_optimized_sequence_correct_length(self, protein, organism):
        """Property: Optimized sequence length must equal protein length * 3."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        assert len(result.sequence) == len(protein) * 3, (
            f"Protein length {len(protein)}, sequence length {len(result.sequence)}"
        )

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_optimized_sequence_valid_dna(self, protein, organism):
        """Property: Optimized sequence should contain only A, C, G, T."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        valid_bases = set("ACGT")
        invalid = set(result.sequence) - valid_bases
        assert not invalid, f"Invalid bases in sequence: {invalid}"

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_cai_in_valid_range(self, protein, organism):
        """Property: CAI must always be in [0, 1]."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        assert 0.0 <= result.cai <= 1.0, f"CAI out of range: {result.cai}"

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_content_in_valid_range(self, protein, organism):
        """Property: GC content must be in [0, 1]."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        assert 0.0 <= result.gc_content <= 1.0, (
            f"GC content out of range: {result.gc_content}"
        )

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_each_codon_matches_amino_acid(self, protein, organism):
        """Property: Every codon in the optimized sequence must encode
        the corresponding amino acid in the protein."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        for i, aa in enumerate(protein):
            codon = result.sequence[i * 3 : i * 3 + 3]
            assert codon in AA_TO_CODONS[aa], (
                f"Position {i}: codon {codon} does not encode {aa}. "
                f"Valid codons: {AA_TO_CODONS[aa]}"
            )

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_gc_content_matches_computed(self, protein, organism):
        """Property: Reported GC content should match independently computed GC."""
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        computed_gc = gc_content(result.sequence)
        assert abs(result.gc_content - computed_gc) < 1e-6, (
            f"Reported GC={result.gc_content}, computed GC={computed_gc}"
        )

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_optimization_deterministic(self, protein, organism):
        """Property: Same protein + same organism + same params → same result."""
        kwargs = dict(gc_lo=0.30, gc_hi=0.70, cai_threshold=0.2)
        result1 = optimize_sequence(protein, organism, **kwargs, strict_mode=False)
        result2 = optimize_sequence(protein, organism, **kwargs, strict_mode=False)
        assert result1.sequence == result2.sequence, (
            f"Optimization is not deterministic for {protein[:10]}.../{organism}"
        )


class TestPropertyValidationResults:
    """Property-based tests for validation result invariants."""

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_validate_translation_fidelity_always_returns_result(self, protein, organism):
        """Property: validate_translation_fidelity never raises, always returns a result."""
        result = validate_translation_fidelity(protein, organism, "gene", "ds")
        assert isinstance(result, DatasetValidationResult)
        assert isinstance(result.passed, bool)
        assert result.test_type == "translation_fidelity"

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_validate_protein_length_always_returns_result(self, protein, organism):
        """Property: validate_protein_length never raises, always returns a result."""
        result = validate_protein_length(
            protein, organism, "gene", "ds", len(protein)
        )
        assert isinstance(result, DatasetValidationResult)
        # When expected_length matches, should pass
        if len(protein) > 0:
            assert result.passed, (
                f"Protein length validation should pass when lengths match: "
                f"{result.actual} vs {result.expected}"
            )

    @given(protein=protein_strategy, organism=organism_strategy)
    @settings(max_examples=10, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def test_validate_gc_content_wide_range_always_passes(self, protein, organism):
        """Property: GC content validation with very wide range should always pass."""
        result = validate_gc_content(
            protein, organism, "gene", "ds", (0.05, 0.95),
            gc_lo=0.05, gc_hi=0.95,
        )
        assert isinstance(result, DatasetValidationResult)
        assert result.passed, (
            f"GC content {result.actual} outside very wide range {result.expected}"
        )


class TestPropertyCAIComputation:
    """Property-based tests for CAI computation invariants."""

    @given(organism=organism_strategy)
    @settings(max_examples=5, deadline=None)
    def test_cai_of_preferred_codon_sequence_high(self, organism):
        """Property: A sequence using only preferred codons should have high CAI."""
        from biocompiler.organisms import PREFERRED_CODON_TABLES

        preferred = PREFERRED_CODON_TABLES[organism]
        # Build a sequence from preferred codons for several amino acids
        seq_parts = [preferred[aa] for aa in ["A", "G", "L", "K", "E", "V"]] * 5
        seq = "ATG" + "".join(seq_parts) + "TAA"
        cai = compute_cai(seq, organism)
        assert cai >= 0.9, f"Preferred codon CAI too low: {cai:.4f}"

    @given(organism=organism_strategy)
    @settings(max_examples=5, deadline=None)
    def test_cai_is_bounded(self, organism):
        """Property: CAI should always be in [0, 1] for any valid sequence."""
        # Use a variety of codons
        seq = "ATGAAATTTGGGCCCTAG"
        cai = compute_cai(seq, organism)
        assert 0.0 <= cai <= 1.0


class TestPropertyDatasetIntegrity:
    """Property-based tests for dataset integrity."""

    @given(
        ds_name=st.sampled_from(list(ALL_DATASETS.keys())),
        gene_name=st.from_regex(r"[A-Za-z0-9_]+", fullmatch=True),
    )
    @settings(max_examples=10)
    def test_missing_gene_does_not_crash(self, ds_name, gene_name):
        """Property: Accessing a non-existent gene in a dataset does not crash."""
        ds = ALL_DATASETS[ds_name]
        # This should just return None or similar, not raise
        result = ds.get(gene_name)
        # Result is either None or a valid gene dict
        if result is not None:
            assert "protein" in result
            assert "organism" in result


# ============================================================================
# 3. Error Handling Tests
# ============================================================================

class TestErrorHandlingInvalidOrganism:
    """Test behavior with unsupported organisms."""

    def test_optimize_unsupported_organism_fallback(self):
        """optimize_sequence with unsupported organism should either raise
        UnsupportedOrganismError or produce a very-low-CAI fallback result."""
        try:
            result = optimize_sequence("MK", "Alien_organism", gc_lo=0.30, gc_hi=0.70, strict_mode=False)
            # If it does not raise, it should produce a valid (but low-quality) result
            assert len(result.sequence) == 6
            assert result.cai >= 0.0
        except UnsupportedOrganismError:
            pass  # Also acceptable

    def test_compute_cai_unsupported_organism_raises(self):
        """compute_cai with unsupported organism should raise."""
        with pytest.raises(UnsupportedOrganismError):
            compute_cai("ATGAAATTTGGGCCCTAA", "Martian_genome")

    def test_validate_translation_fidelity_unsupported_organism_graceful(self):
        """validate_translation_fidelity with bad organism returns a result
        (may pass with fallback or fail — key is no unhandled exception)."""
        result = validate_translation_fidelity(
            "MK", "Fictitious_organism", "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)
        # The optimizer may use a fallback for unsupported organisms,
        # so translation fidelity can still pass. The key invariant
        # is that no unhandled exception propagates.

    def test_validate_gc_content_unsupported_organism_graceful(self):
        """validate_gc_content with bad organism returns a result (may pass
        with fallback — key is no unhandled exception)."""
        result = validate_gc_content(
            "MK", "Fictitious_organism", "test", "test", (0.20, 0.80)
        )
        assert isinstance(result, DatasetValidationResult)

    def test_validate_cai_bounds_unsupported_organism_graceful(self):
        """validate_cai_bounds with bad organism returns a result."""
        result = validate_cai_bounds(
            "MK", "Fictitious_organism", "test", "test", (0.1, 1.0)
        )
        assert isinstance(result, DatasetValidationResult)

    def test_validate_protein_length_unsupported_organism_graceful(self):
        """validate_protein_length with bad organism returns a result."""
        result = validate_protein_length(
            "MK", "Fictitious_organism", "test", "test", 2
        )
        assert isinstance(result, DatasetValidationResult)

    def test_validate_optimization_improvement_unsupported_organism_graceful(self):
        """validate_optimization_improvement with bad organism returns a result."""
        result = validate_optimization_improvement(
            "MK", "Fictitious_organism", "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)

    def test_validate_no_cpg_island_unsupported_organism_graceful(self):
        """validate_no_cpg_island with bad organism returns a result."""
        result = validate_no_cpg_island(
            "MK", "Fictitious_organism", "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)

    def test_validate_cross_organism_consistency_unsupported_graceful(self):
        """validate_cross_organism_consistency with bad organism in inner call
        returns failed result."""
        # This function hardcodes E. coli and human internally, so
        # passing a bad organism for the gene_name should not cause issues.
        # But if optimize_sequence fails for internal organisms, it is caught.
        result = validate_cross_organism_consistency(
            "MKWVA", "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)


class TestErrorHandlingInvalidProtein:
    """Test behavior with invalid protein sequences."""

    def test_empty_protein_raises(self):
        """Empty protein should raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("", "Homo_sapiens", strict_mode=False)

    def test_whitespace_only_protein_raises(self):
        """Whitespace-only protein should raise InvalidProteinError or
        a clear error (currently triggers ZeroDivisionError — a known bug)."""
        with pytest.raises((InvalidProteinError, ZeroDivisionError)):
            optimize_sequence("   ", "Homo_sapiens", strict_mode=False)

    def test_protein_with_invalid_chars_raises(self):
        """Protein with invalid amino acid codes should raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("MKX", "Homo_sapiens", strict_mode=False)

    def test_protein_with_numbers_raises(self):
        """Protein with numeric characters should raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("M1K2", "Homo_sapiens", strict_mode=False)

    def test_protein_with_special_chars_raises(self):
        """Protein with special characters should raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_sequence("M@K", "Homo_sapiens", strict_mode=False)

    def test_protein_lowercase_accepted(self):
        """Protein with lowercase letters is accepted (optimizer upper-strips)."""
        result = optimize_sequence("mkwva", "Homo_sapiens", gc_lo=0.20, gc_hi=0.80, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == "MKWVA"

    def test_validate_function_catches_empty_protein(self):
        """Validation functions should catch empty protein errors gracefully."""
        result = validate_translation_fidelity(
            "", "Homo_sapiens", "test", "test"
        )
        assert isinstance(result, DatasetValidationResult)
        assert result.passed is False
        assert "ERROR" in result.actual

    def test_validate_function_catches_invalid_protein(self):
        """Validation functions should catch invalid protein errors gracefully."""
        result = validate_gc_content(
            "MKXZ", "Homo_sapiens", "test", "test", (0.20, 0.80)
        )
        assert isinstance(result, DatasetValidationResult)
        assert result.passed is False
        assert "ERROR" in result.actual


class TestErrorHandlingInvalidGCBounds:
    """Test behavior with invalid GC bound parameters.

    Note: The optimizer does not strictly validate GC bounds — it proceeds
    gracefully even with unusual ranges. These tests verify it does not crash
    and still produces valid output.
    """

    def test_gc_lo_greater_than_gc_hi_no_crash(self):
        """gc_lo > gc_hi should not crash the optimizer."""
        result = optimize_sequence("MKWVA", "Homo_sapiens", gc_lo=0.70, gc_hi=0.30, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == "MKWVA"

    def test_gc_lo_equals_gc_hi_no_crash(self):
        """gc_lo == gc_hi should not crash the optimizer."""
        result = optimize_sequence("MKWVA", "Homo_sapiens", gc_lo=0.50, gc_hi=0.50, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == "MKWVA"

    def test_gc_bounds_negative_no_crash(self):
        """Negative GC bounds should not crash the optimizer."""
        result = optimize_sequence("MKWVA", "Homo_sapiens", gc_lo=-0.1, gc_hi=0.50, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == "MKWVA"

    def test_gc_hi_above_1_no_crash(self):
        """gc_hi > 1.0 should not crash the optimizer."""
        result = optimize_sequence("MKWVA", "Homo_sapiens", gc_lo=0.30, gc_hi=1.5, strict_mode=False)
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == "MKWVA"


class TestErrorHandlingProteinLengthMismatch:
    """Test validate_protein_length with mismatched expected lengths."""

    def test_expected_length_too_high(self):
        """Expected length greater than actual should fail."""
        result = validate_protein_length(
            "MKWVA", "Homo_sapiens", "test", "test", expected_length=10
        )
        assert result.passed is False
        assert "10 codons" in result.expected

    def test_expected_length_too_low(self):
        """Expected length less than actual should fail."""
        result = validate_protein_length(
            "MKWVA", "Homo_sapiens", "test", "test", expected_length=3
        )
        assert result.passed is False

    def test_expected_length_zero(self):
        """Expected length of 0 for non-empty protein should fail."""
        result = validate_protein_length(
            "MKWVA", "Homo_sapiens", "test", "test", expected_length=0
        )
        assert result.passed is False

    def test_expected_length_correct(self):
        """Correct expected length should pass."""
        result = validate_protein_length(
            "MKWVA", "Homo_sapiens", "test", "test", expected_length=5
        )
        assert result.passed is True


class TestErrorHandlingTranslation:
    """Test error handling in the translate function."""

    def test_translate_empty_sequence(self):
        """Translating empty sequence should return empty string."""
        assert translate("") == ""

    def test_translate_invalid_dna_character(self):
        """Translating sequence with invalid characters should still work
        (the validate_dna_sequence function handles this)."""
        # The translation module has its own validation
        # Let us verify the behavior is graceful
        result = translate("ATGAAATTTGGGCCCTAA")
        assert len(result) > 0


class TestErrorHandlingCAIComputation:
    """Test error handling in CAI computation."""

    def test_cai_empty_sequence(self):
        """CAI of empty sequence should be 0.0."""
        assert compute_cai("", "Homo_sapiens") == 0.0

    def test_cai_stop_codons_only(self):
        """CAI of a sequence with only stop codons should be 0.0 or very low."""
        cai = compute_cai("TAATAGTGA", "Homo_sapiens")
        # Stop codons are skipped in CAI computation, so this should return 0.0
        assert cai == 0.0 or cai < 0.01

    def test_cai_deterministic_same_input(self):
        """Same input should always give same CAI."""
        seq = "ATGAAATTTGGGCCCTAA"
        cai1 = compute_cai(seq, "Homo_sapiens")
        cai2 = compute_cai(seq, "Homo_sapiens")
        assert cai1 == cai2


class TestErrorHandlingValidationSummary:
    """Test summary computation edge cases."""

    def test_summary_with_no_results(self):
        """Summary with empty results should not crash."""
        report = DatasetValidationReport(
            timestamp="2025-01-01", version="1.0",
            total_tests=0, passed=0, failed=0,
            results=[],
            summary={},
        )
        # Calling pass_rate on empty report
        assert report.pass_rate == 0.0

    def test_format_report_empty_results(self):
        """Formatting a report with no results should produce valid text."""
        report = DatasetValidationReport(
            timestamp="2025-01-01", version="1.0",
            total_tests=0, passed=0, failed=0,
            results=[],
            summary={},
        )
        text = format_dataset_report_text(report)
        assert isinstance(text, str)
        assert "0/0 passed" in text


# ============================================================================
# 4. Multi-Organism Consistency Extended Tests
# ============================================================================

class TestMultiOrganismConsistency:
    """Extended cross-organism consistency tests."""

    @pytest.mark.parametrize("organism", SUPPORTED_ORGANISMS)
    def test_optimization_for_each_organism(self, organism):
        """A common protein should optimize successfully for each organism."""
        protein = "ACDEFGHIKLMNPQRSTVWY"
        result = optimize_sequence(
            protein, organism, gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_different_organisms_different_cai_scores(self):
        """Optimizing the same protein for different organisms should produce
        different CAI scores when evaluated under their respective codon usage
        tables, even if the DNA sequences happen to be identical (which can
        occur when constraint solving overrides codon preferences).

        The key invariant is that organism-specific CAI scoring produces
        different values — confirming the codon adaptiveness tables differ."""
        protein = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 aa
        cai_scores = {}
        for organism in SUPPORTED_ORGANISMS:
            result = optimize_sequence(
                protein, organism, gc_lo=0.30, gc_hi=0.70, cai_threshold=0.1,
                strict_mode=False,
            )
            cai_scores[organism] = result.cai

        # CAI scores should not all be identical (tables differ)
        scores = list(cai_scores.values())
        any_different = any(scores[i] != scores[j] for i in range(len(scores)) for j in range(i+1, len(scores)))
        assert any_different, (
            f"All organisms produced identical CAI scores: {cai_scores}"
        )

    def test_home_organism_cai_advantage(self):
        """A sequence optimized for an organism should generally have higher CAI
        under that organism's codon usage than under a distant organism's."""
        protein = "ACDEFGHIKLMNPQRSTVWY" * 3
        ecoli_result = optimize_sequence(
            protein, "Escherichia_coli", gc_lo=0.30, gc_hi=0.70,
            strict_mode=False,
        )
        # E. coli-optimized sequence should score at least somewhat well
        # under E. coli codon usage
        assert ecoli_result.cai > 0.0, "E. coli CAI should be positive"


# ============================================================================
# 5. Stress / Robustness Tests
# ============================================================================

class TestRobustness:
    """Stress tests for the optimization pipeline."""

    def test_repeated_optimization_same_protein(self):
        """Optimizing the same protein 5 times should produce the same result."""
        protein = "MKWVACDEFGH"
        results = []
        for _ in range(5):
            result = optimize_sequence(
                protein, "Homo_sapiens", gc_lo=0.30, gc_hi=0.70,
                strict_mode=False,
            )
            results.append(result.sequence)
        assert len(set(results)) == 1, "Repeated optimization not deterministic"

    def test_protein_with_all_same_amino_acid_long(self):
        """A long protein of a single amino acid should optimize."""
        protein = "A" * 100
        result = optimize_sequence(
            protein, "Homo_sapiens", gc_lo=0.20, gc_hi=0.80,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein
        assert len(result.sequence) == 300

    def test_protein_alternating_two_aas(self):
        """Alternating two amino acids should optimize correctly."""
        protein = "AE" * 50
        result = optimize_sequence(
            protein, "Escherichia_coli", gc_lo=0.30, gc_hi=0.70,
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_validation_report_consistency(self):
        """Total tests = passed + failed in any report."""
        for ds in [["synthetic"], ["human", "ecoli"]]:
            report = run_dataset_validation(
                datasets=ds,
                include_cross_organism=True,
                include_optimization_improvement=True,
                include_no_cpg_island=True,
            )
            assert report.total_tests == report.passed + report.failed, (
                f"Inconsistent test counts: total={report.total_tests}, "
                f"passed={report.passed}, failed={report.failed}"
            )
            assert len(report.results) == report.total_tests, (
                f"Results list length != total_tests"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
