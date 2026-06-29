"""Regression tests for wet-lab validation benchmarks.

These tests run the optimizer on known protein sequences and verify
that the CAI scores are within expected ranges. Any regression in
CAI scores (e.g., due to GT-avoidance-over-CAI bugs) will be caught.

Thresholds are set based on current optimizer behavior with a safety
margin to catch regressions. For eukaryotic organisms, GT avoidance
and CpG elimination depress CAI relative to the theoretical maximum,
so thresholds are adjusted accordingly.

Run with:
    pytest tests/test_wetlab_benchmark_regression.py -v
    pytest tests/test_wetlab_benchmark_regression.py -v -m benchmark  # slow benchmarks
    pytest tests/test_wetlab_benchmark_regression.py -v -m "not slow"  # fast only
"""
import pytest


@pytest.mark.slow
@pytest.mark.benchmark
class TestWetLabBenchmarkRegression:
    """Regression tests against published experimental data."""

    def test_ecoli_gfp_cai_above_threshold(self):
        """E. coli GFP should achieve CAI >= 0.90 after optimization."""
        from biocompiler.optimizer.pipeline import optimize_sequence

        gfp = (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"
            "LVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDT"
            "LVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQL"
            "ADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        )
        result = optimize_sequence(gfp, organism="E. coli")
        assert result.cai >= 0.90, (
            f"E. coli GFP CAI={result.cai:.4f}, expected >= 0.90"
        )

    def test_human_insulin_cai_above_threshold(self):
        """Human insulin in human should achieve CAI >= 0.65.

        Note: Eukaryotic GT avoidance and CpG elimination significantly
        depress CAI for insulin (which contains many Val/Ala residues
        encoded by GTN/GCN codons). The threshold reflects the current
        optimizer's CAI after all eukaryotic constraints are applied.
        A regression below 0.65 would indicate a bug in the CAI
        recovery pass or an over-aggressive constraint.
        """
        from biocompiler.optimizer.pipeline import optimize_sequence

        insulin = (
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
            "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        )
        result = optimize_sequence(insulin, organism="Homo_sapiens")
        assert result.cai >= 0.65, (
            f"Human insulin CAI={result.cai:.4f}, expected >= 0.65 — "
            f"GT avoidance or CpG elimination may be too aggressive"
        )

    def test_yeast_insulin_cai_above_threshold(self):
        """Yeast insulin should achieve CAI >= 0.80 (GT avoidance should not tank CAI below this).

        Note: Yeast GT avoidance depresses CAI for insulin (many Val/Ala
        residues with GTN/GCN codons). The threshold reflects current
        optimizer behavior after eukaryotic constraints.
        """
        from biocompiler.optimizer.pipeline import optimize_sequence

        insulin = (
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
            "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        )
        result = optimize_sequence(insulin, organism="Saccharomyces_cerevisiae")
        assert result.cai >= 0.80, (
            f"Yeast insulin CAI={result.cai:.4f}, expected >= 0.80 "
            f"— GT avoidance may be too aggressive"
        )

    def test_ecoli_hbb_cai_above_threshold(self):
        """E. coli HBB should achieve CAI >= 0.90."""
        from biocompiler.optimizer.pipeline import optimize_sequence

        hbb = (
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
            "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
            "KEFTPPVQAAYQKVVAGVANALAHKYH"
        )
        result = optimize_sequence(hbb, organism="E. coli")
        assert result.cai >= 0.90, (
            f"E. coli HBB CAI={result.cai:.4f}, expected >= 0.90"
        )

    def test_yeast_gfp_cai_not_depressed_by_gt_avoidance(self):
        """Yeast GFP should have CAI >= 0.85 — GT avoidance should not tank CAI.

        Uses strict_mode=False because yeast optimization of GFP may
        produce GCInRange/SlidingGC warnings that do not affect CAI.
        """
        from biocompiler.optimizer.pipeline import optimize_sequence

        gfp = (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"
            "LVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDT"
            "LVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQL"
            "ADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        )
        result = optimize_sequence(gfp, organism="Saccharomyces_cerevisiae", strict_mode=False)
        assert result.cai >= 0.85, (
            f"Yeast GFP CAI={result.cai:.4f}, expected >= 0.85"
        )

    def test_ecoli_insulin_cai_above_threshold(self):
        """E. coli insulin should achieve CAI >= 0.90."""
        from biocompiler.optimizer.pipeline import optimize_sequence

        insulin = (
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
            "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        )
        result = optimize_sequence(insulin, organism="E. coli")
        assert result.cai >= 0.90, (
            f"E. coli insulin CAI={result.cai:.4f}, expected >= 0.90"
        )

    def test_human_hbb_cai_above_threshold(self):
        """Human HBB should achieve CAI >= 0.80 in human.

        Note: Eukaryotic GT avoidance and CpG elimination depress CAI
        for HBB. The threshold reflects current optimizer behavior.
        """
        from biocompiler.optimizer.pipeline import optimize_sequence

        hbb = (
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
            "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
            "KEFTPPVQAAYQKVVAGVANALAHKYH"
        )
        result = optimize_sequence(hbb, organism="Homo_sapiens")
        assert result.cai >= 0.80, (
            f"Human HBB CAI={result.cai:.4f}, expected >= 0.80 — "
            f"GT avoidance or CpG elimination may be too aggressive"
        )


@pytest.mark.benchmark
class TestBenchmarkRunner:
    """Tests for the benchmark runner infrastructure itself."""

    def test_published_expression_data_imports(self):
        """Published expression data module should import cleanly."""
        from biocompiler.validation.published_expression_data import (
            PublishedExpressionResult,
            KUDLA_2009_GFP_DATASET,
            WELCH_2009_DATASET,
            PUIGBO_2008_DATASET,
            SHARP_LI_1987_DATASET,
            ALL_PUBLISHED_DATASETS,
        )
        assert len(ALL_PUBLISHED_DATASETS) >= 4
        assert len(KUDLA_2009_GFP_DATASET) > 0
        assert len(WELCH_2009_DATASET) > 0
        assert len(PUIGBO_2008_DATASET) > 0
        assert len(SHARP_LI_1987_DATASET) > 0

    def test_benchmark_runner_imports(self):
        """Benchmark runner module should import cleanly."""
        from biocompiler.validation.benchmark_runner import (
            BenchmarkResult,
            BenchmarkReport,
            run_single_benchmark,
            run_benchmark_suite,
            format_report_text,
        )

    def test_published_data_has_valid_entries(self):
        """All published expression data entries should have valid fields."""
        from biocompiler.validation.published_expression_data import (
            ALL_PUBLISHED_DATASETS,
        )

        for dataset_name, entries in ALL_PUBLISHED_DATASETS.items():
            assert len(entries) > 0, f"Dataset {dataset_name} is empty"
            for entry in entries:
                assert entry.gene_name, f"Empty gene_name in {dataset_name}"
                assert entry.organism, f"Empty organism in {dataset_name}"
                assert entry.protein_sequence, f"Empty protein_sequence in {dataset_name}"
                assert 0.0 <= entry.cai_predicted <= 1.0, (
                    f"Invalid cai_predicted={entry.cai_predicted} in {dataset_name}"
                )
                assert entry.source, f"Empty source in {dataset_name}"
                assert entry.doi, f"Empty doi in {dataset_name}"

    def test_format_report_text(self):
        """format_report_text should produce a non-empty string."""
        from biocompiler.validation.benchmark_runner import (
            BenchmarkReport,
            BenchmarkResult,
            format_report_text,
        )

        report = BenchmarkReport(
            timestamp="2025-01-01T00:00:00Z",
            biocompiler_version="1.0.0",
            results=[
                BenchmarkResult(
                    gene_name="GFP",
                    organism="Escherichia_coli",
                    measured_expression=1.0,
                    predicted_cai=0.95,
                    expected_cai=0.90,
                    cai_error=0.05,
                    cai_relative_error=0.056,
                    optimization_time_seconds=0.01,
                    passed=True,
                ),
            ],
            summary={
                "total_benchmarks": 1,
                "passed": 1,
                "failed": 0,
                "pass_rate": 1.0,
                "mean_cai_error": 0.05,
                "max_cai_error": 0.05,
                "cai_tolerance": 0.05,
            },
        )
        text = format_report_text(report)
        assert "BIOCOMPILER WET-LAB VALIDATION BENCHMARK REPORT" in text
        assert "GFP" in text
        assert "SUMMARY" in text

    @pytest.mark.slow
    def test_run_single_benchmark_ecoli_gfp(self):
        """run_single_benchmark should return a BenchmarkResult for E. coli GFP."""
        from biocompiler.validation.benchmark_runner import run_single_benchmark

        result = run_single_benchmark(
            protein="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
            organism="Escherichia_coli",
            expected_cai=0.95,
            gene_name="GFP",
            cai_tolerance=0.10,
        )
        assert result.gene_name == "GFP"
        assert result.organism == "Escherichia_coli"
        assert 0.0 < result.predicted_cai <= 1.0
        assert result.optimization_time_seconds > 0
