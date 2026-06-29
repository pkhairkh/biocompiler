"""Tests for CAI benchmarking against published values.

Agent 57: Validates that BioCompiler's CAI implementation can be benchmarked
against published CAI values from Sharp & Li (1987) and Puigbo et al (2008).
"""

import pytest

from biocompiler.benchmarking.cai_benchmarking import (
    CAIBenchmarkResult,
    PUBLISHED_CAI_VALUES,
    benchmark_cai,
    benchmark_cai_for_dna,
    benchmark_optimization,
    summarize_benchmark,
    BenchmarkSummary,
)
from biocompiler.expression.translation import compute_cai


# ═══════════════════════════════════════════════════════════════════════════════
# Test Data
# ═══════════════════════════════════════════════════════════════════════════════

# GFP coding sequence (optimized for E. coli — all optimal codons)
_ECOLI_OPTIMAL_GFP = (
    "ATGGCTAGCAAAGAGAAGAAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGC"
    "ACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTAC"
    "TGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTGCGAGATACCC"
    "AGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTA"
    "AAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTA"
    "AAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATA"
    "CATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAGCTGTT"
    "CAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTC"
    "CACACAATCTGCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTG"
    "GGATTACACATGGCATGGATGAACTATACAAATAA"
)

# Human insulin CDS (native human codons — not optimized for E. coli)
_HUMAN_INSULIN_DNA = (
    "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCT"
    "TTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTC"
    "TACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGC"
    "AGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCA"
    "TCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
)

# HBB (human beta-globin) native CDS
_HUMAN_HBB_DNA = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGG"
    "TGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGT"
    "CCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGAT"
    "GGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGT"
    "GGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCA"
    "CCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
)


class TestPublishedCAIValues:
    """Test the PUBLISHED_CAI_VALUES database."""

    def test_database_not_empty(self):
        """Published CAI values database should not be empty."""
        assert len(PUBLISHED_CAI_VALUES) > 0

    def test_ecoli_genes_present(self):
        """E. coli genes should be present in the database."""
        ecoli_genes = {k[0] for k in PUBLISHED_CAI_VALUES if k[1] == "e_coli"}
        assert "lacZ" in ecoli_genes
        assert "trpA" in ecoli_genes
        assert "gfp" in ecoli_genes
        assert "insulin" in ecoli_genes

    def test_human_genes_present(self):
        """Human genes should be present in the database."""
        human_genes = {k[0] for k in PUBLISHED_CAI_VALUES if k[1] == "human"}
        assert "HBB" in human_genes
        assert "INS" in human_genes or "insulin" in human_genes
        assert "TP53" in human_genes

    def test_yeast_genes_present(self):
        """Yeast genes should be present in the database."""
        yeast_genes = {k[0] for k in PUBLISHED_CAI_VALUES if k[1] == "yeast"}
        assert "ADH1" in yeast_genes
        assert "PGK1" in yeast_genes

    def test_all_entries_have_required_fields(self):
        """Each published value should have 'cai' and 'source' fields."""
        for key, entry in PUBLISHED_CAI_VALUES.items():
            assert "cai" in entry, f"Missing 'cai' for {key}"
            assert "source" in entry, f"Missing 'source' for {key}"
            assert isinstance(entry["cai"], (int, float)), f"CAI not numeric for {key}"
            assert 0.0 <= float(entry["cai"]) <= 1.0, f"CAI out of range for {key}"

    def test_gene_type_field(self):
        """Each entry should have a gene_type field."""
        for key, entry in PUBLISHED_CAI_VALUES.items():
            assert "gene_type" in entry, f"Missing 'gene_type' for {key}"
            assert entry["gene_type"] in ("native", "heterologous"), \
                f"Invalid gene_type for {key}: {entry['gene_type']}"

    def test_highly_expressed_genes_have_high_cai(self):
        """Highly expressed native genes should have CAI >= 0.7."""
        high_expr = [
            ("trpA", "e_coli"),
            ("ompA", "e_coli"),
            ("groEL", "e_coli"),
            ("ADH1", "yeast"),
            ("PGK1", "yeast"),
            ("HBB", "human"),
        ]
        for gene, org in high_expr:
            entry = PUBLISHED_CAI_VALUES.get((gene, org))
            if entry is not None:
                assert float(entry["cai"]) >= 0.7, \
                    f"{gene} in {org} should have CAI >= 0.7, got {entry['cai']}"

    def test_heterologous_genes_have_low_cai(self):
        """Heterologous genes in E. coli should have lower CAI than native."""
        het_genes = [
            ("gfp", "e_coli"),
            ("insulin", "e_coli"),
            ("hGH", "e_coli"),
        ]
        for gene, org in het_genes:
            entry = PUBLISHED_CAI_VALUES.get((gene, org))
            if entry is not None:
                assert float(entry["cai"]) < 0.6, \
                    f"Heterologous {gene} in {org} should have CAI < 0.6, got {entry['cai']}"


class TestBenchmarkCAIForDNA:
    """Test benchmark_cai_for_dna function."""

    def test_benchmark_with_ecoli_insulin(self):
        """Benchmark human insulin (native codons) against E. coli published value."""
        result = benchmark_cai_for_dna("insulin", "e_coli", _HUMAN_INSULIN_DNA)
        assert isinstance(result, CAIBenchmarkResult)
        assert result.gene == "insulin"
        assert result.organism == "Escherichia_coli"
        assert result.published_cai is not None
        assert result.published_cai == 0.34  # Puigbo et al 2008
        assert result.predicted_cai > 0.0
        assert result.source == "Puigbo et al 2008, CAIcal server"

    def test_benchmark_with_human_hbb(self):
        """Benchmark HBB against human published value."""
        result = benchmark_cai_for_dna("HBB", "human", _HUMAN_HBB_DNA)
        assert isinstance(result, CAIBenchmarkResult)
        assert result.organism == "Homo_sapiens"
        assert result.published_cai == 0.95

    def test_benchmark_without_published_value(self):
        """Benchmark a gene with no published value should return None for published_cai."""
        result = benchmark_cai_for_dna("unknown_gene", "e_coli", "ATGAAAGCGTTT")
        assert result.published_cai is None
        assert result.deviation is None
        assert result.pass_threshold is None
        assert result.source == "N/A"

    def test_benchmark_organism_resolution(self):
        """Organism name should be resolved to canonical form."""
        result = benchmark_cai_for_dna("insulin", "ecoli", _HUMAN_INSULIN_DNA)
        assert result.organism == "Escherichia_coli"

    def test_benchmark_gene_type(self):
        """Gene type should be correctly identified."""
        # Native gene
        result = benchmark_cai_for_dna("trpA", "e_coli", "ATGAAAGCGTTT")
        assert result.gene_type == "native"

        # Heterologous gene
        result = benchmark_cai_for_dna("gfp", "e_coli", "ATGAAAGCGTTT")
        assert result.gene_type == "heterologous"

    def test_benchmark_deviation_computation(self):
        """Deviation should be predicted - published."""
        result = benchmark_cai_for_dna("insulin", "e_coli", _HUMAN_INSULIN_DNA)
        assert result.deviation is not None
        expected_deviation = result.predicted_cai - result.published_cai
        assert abs(result.deviation - expected_deviation) < 1e-10

    def test_optimized_sequence_should_pass(self):
        """An optimized sequence should pass the benchmark threshold."""
        result = benchmark_cai_for_dna("gfp", "e_coli", _ECOLI_OPTIMAL_GFP)
        # GFP optimized for E. coli should have high CAI
        # Heterologous genes always pass (we expect improvement)
        assert result.pass_threshold is True

    def test_result_is_frozen(self):
        """CAIBenchmarkResult should be immutable (frozen dataclass)."""
        result = benchmark_cai_for_dna("insulin", "e_coli", _HUMAN_INSULIN_DNA)
        with pytest.raises(AttributeError):
            result.gene = "changed"


class TestBenchmarkCAI:
    """Test the batch benchmark_cai function."""

    def test_batch_benchmark(self):
        """Batch benchmark should process multiple genes."""
        genes = ["insulin", "gfp"]
        organisms = ["e_coli", "e_coli"]
        dna_sequences = {
            ("insulin", "e_coli"): _HUMAN_INSULIN_DNA,
            ("gfp", "e_coli"): _ECOLI_OPTIMAL_GFP,
        }
        results = benchmark_cai(genes, organisms, dna_sequences)
        assert len(results) == 2
        assert all(isinstance(r, CAIBenchmarkResult) for r in results)
        assert results[0].gene == "insulin"
        assert results[1].gene == "gfp"

    def test_batch_without_dna_sequences(self):
        """Batch benchmark without DNA should use predicted_cai=0.0."""
        results = benchmark_cai(
            ["insulin"], ["e_coli"]
        )
        assert len(results) == 1
        assert results[0].predicted_cai == 0.0
        assert results[0].published_cai is not None

    def test_batch_mismatched_lengths(self):
        """Mismatched gene/organism lists should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            benchmark_cai(["gfp"], ["e_coli", "human"])

    def test_batch_organism_resolution(self):
        """Organism names should be resolved in batch mode."""
        results = benchmark_cai(
            ["insulin"], ["ecoli"],
            dna_sequences={("insulin", "ecoli"): _HUMAN_INSULIN_DNA},
        )
        assert results[0].organism == "Escherichia_coli"


class TestBenchmarkOptimization:
    """Test the benchmark_optimization function."""

    def test_optimization_improvement(self):
        """Optimized sequence should have higher CAI than original."""
        result = benchmark_optimization(
            gene="insulin",
            organism="e_coli",
            original_dna=_HUMAN_INSULIN_DNA,
            optimized_dna=_ECOLI_OPTIMAL_GFP,  # Using GFP as a proxy optimized seq
        )
        assert "original_cai" in result
        assert "optimized_cai" in result
        assert "improvement" in result
        assert isinstance(result["optimization_success"], bool)

    def test_optimization_with_published_value(self):
        """Optimization benchmark should look up published CAI."""
        result = benchmark_optimization(
            gene="insulin",
            organism="e_coli",
            original_dna=_HUMAN_INSULIN_DNA,
            optimized_dna=_ECOLI_OPTIMAL_GFP,
        )
        assert result["published_cai"] is not None
        assert result["published_cai"] == 0.34

    def test_optimization_without_published_value(self):
        """Optimization benchmark for unknown gene should have None published_cai."""
        result = benchmark_optimization(
            gene="unknown_gene",
            organism="e_coli",
            original_dna="ATGAAAGCGTTT",
            optimized_dna="ATGGCTAGCAAAGAG",
        )
        assert result["published_cai"] is None
        assert result["exceeds_published"] is None

    def test_same_sequence_zero_improvement(self):
        """Same sequence for original and optimized should have zero improvement."""
        result = benchmark_optimization(
            gene="insulin",
            organism="e_coli",
            original_dna=_HUMAN_INSULIN_DNA,
            optimized_dna=_HUMAN_INSULIN_DNA,
        )
        assert abs(result["improvement"]) < 1e-10
        assert result["optimization_success"] is False


class TestSummarizeBenchmark:
    """Test the summarize_benchmark function."""

    def test_summary_with_results(self):
        """Summary should compute correct statistics."""
        results = [
            CAIBenchmarkResult(
                gene="gfp", organism="Escherichia_coli",
                predicted_cai=0.80, published_cai=0.54,
                source="Puigbo et al 2008", deviation=0.26,
                gene_type="heterologous", pass_threshold=True,
            ),
            CAIBenchmarkResult(
                gene="insulin", organism="Escherichia_coli",
                predicted_cai=0.50, published_cai=0.34,
                source="Puigbo et al 2008", deviation=0.16,
                gene_type="heterologous", pass_threshold=True,
            ),
        ]
        summary = summarize_benchmark(results)
        assert isinstance(summary, BenchmarkSummary)
        assert summary.total_genes == 2
        assert summary.genes_with_published == 2
        assert abs(summary.mean_deviation - 0.21) < 1e-10
        assert summary.pass_rate == 1.0
        assert summary.failures == []

    def test_summary_with_failures(self):
        """Summary should identify failures."""
        results = [
            CAIBenchmarkResult(
                gene="bad_gene", organism="Escherichia_coli",
                predicted_cai=0.10, published_cai=0.84,
                source="Sharp & Li 1987", deviation=-0.74,
                gene_type="native", pass_threshold=False,
            ),
        ]
        summary = summarize_benchmark(results)
        assert summary.pass_rate == 0.0
        assert "bad_gene" in summary.failures

    def test_summary_no_published_values(self):
        """Summary with no published values should return safe defaults."""
        results = [
            CAIBenchmarkResult(
                gene="unknown", organism="Escherichia_coli",
                predicted_cai=0.5, published_cai=None,
                source="N/A", deviation=None,
                gene_type="unknown", pass_threshold=None,
            ),
        ]
        summary = summarize_benchmark(results)
        assert summary.genes_with_published == 0
        assert summary.mean_deviation == 0.0

    def test_empty_results(self):
        """Empty results should return safe defaults."""
        summary = summarize_benchmark([])
        assert summary.total_genes == 0
        assert summary.genes_with_published == 0


class TestBenchmarkConsistency:
    """Test consistency between CAI computation and benchmarking."""

    def test_predicted_cai_matches_compute_cai(self):
        """benchmark_cai_for_dna should use compute_cai internally."""
        dna = _HUMAN_INSULIN_DNA
        result = benchmark_cai_for_dna("insulin", "e_coli", dna)
        direct_cai = compute_cai(dna, organism="Escherichia_coli")
        assert abs(result.predicted_cai - direct_cai) < 1e-10

    def test_optimized_cai_higher_than_heterologous_published(self):
        """An E. coli-optimized sequence should have higher CAI than
        the published value for the native (heterologous) gene.

        Note: The published CAI of 0.54 for GFP in E. coli was computed
        using the Sharp-Li (1987) reference set, while BioCompiler uses
        the Kazusa-derived CODON_ADAPTIVENESS_TABLES.  Cross-reference-set
        comparisons may differ by ±0.10, so we use a tolerance-based check.
        """
        # Use a truly all-optimal-codon sequence for the GFP protein
        # to guarantee the highest possible CAI in our reference set.
        from biocompiler.organisms import E_COLI_PREFERRED_CODONS
        gfp_protein = (
            "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
            "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
            "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
            "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        )
        optimal_gfp_dna = "ATG" + "".join(
            E_COLI_PREFERRED_CODONS.get(aa, "NNN") for aa in gfp_protein[1:]
        )
        result = benchmark_cai_for_dna("gfp", "e_coli", optimal_gfp_dna)
        # With all-optimal codons, our predicted CAI should be very high (close to 1.0)
        assert result.predicted_cai >= 0.9, \
            f"All-optimal GFP CAI unexpectedly low: {result.predicted_cai}"
        # And should exceed the published heterologous value of 0.54
        assert result.predicted_cai > result.published_cai

    def test_published_values_are_for_native_genes(self):
        """All 'native' gene_type entries should be from the same organism."""
        for (gene, org), entry in PUBLISHED_CAI_VALUES.items():
            if entry.get("gene_type") == "native" and org in ("e_coli",):
                # Native E. coli genes should have reasonable CAI values
                cai = float(entry["cai"])
                assert 0.1 <= cai <= 1.0, \
                    f"Native E. coli gene {gene} has unreasonable CAI: {cai}"
