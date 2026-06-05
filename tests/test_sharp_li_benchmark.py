"""
Tests for the Sharp-Li vs Kazusa CAI Benchmark
=================================================

Validates that:
1. benchmark_sharp_li_cai() returns valid results for E. coli genes
2. Sharp-Li CAI values are closer to published values than Kazusa (on average)
3. The benchmark report can be printed without errors
"""

from __future__ import annotations

import math

import pytest

from biocompiler.benchmarking.sharp_li_benchmark import (
    SHARP_LI_ECOLI_REFERENCE,
    benchmark_sharp_li_cai,
    print_benchmark_report,
)
from biocompiler.benchmarking.cai_validated import (
    compute_cai_sharp_li,
    _compute_adaptiveness_table,
    load_reference_set,
)
from biocompiler.benchmarking.cai_published_values import (
    PUBLISHED_CAI_VALUES,
    VALIDATION_SEQUENCES,
)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Benchmark returns valid results for E. coli genes
# ═══════════════════════════════════════════════════════════════════════════════

class TestBenchmarkResults:
    """Validate that benchmark_sharp_li_cai returns well-structured results."""

    @pytest.fixture(scope="class")
    def benchmark_results(self):
        """Run the benchmark once and share results across tests."""
        return benchmark_sharp_li_cai()

    def test_returns_dict(self, benchmark_results):
        """The benchmark returns a dict."""
        assert isinstance(benchmark_results, dict)

    def test_has_required_keys(self, benchmark_results):
        """The benchmark dict has all required keys."""
        required_keys = {
            "per_gene_results",
            "mean_kazusa_error",
            "mean_sharp_li_error",
            "sharp_li_is_closer",
        }
        assert required_keys.issubset(benchmark_results.keys())

    def test_per_gene_results_is_list(self, benchmark_results):
        """per_gene_results is a list."""
        assert isinstance(benchmark_results["per_gene_results"], list)

    def test_at_least_one_ecoli_gene(self, benchmark_results):
        """At least one E. coli gene with a DNA sequence is benchmarked."""
        assert len(benchmark_results["per_gene_results"]) >= 1

    def test_per_gene_result_structure(self, benchmark_results):
        """Each per-gene result has the required keys."""
        required_keys = {
            "gene", "organism", "published_cai",
            "kazusa_cai", "sharp_li_cai",
            "kazusa_error", "sharp_li_error",
        }
        for result in benchmark_results["per_gene_results"]:
            assert required_keys.issubset(result.keys()), (
                f"Missing keys in result for {result.get('gene', 'unknown')}: "
                f"{required_keys - result.keys()}"
            )

    def test_all_genes_are_ecoli(self, benchmark_results):
        """All benchmarked genes are E. coli genes."""
        for result in benchmark_results["per_gene_results"]:
            assert result["organism"] == "Escherichia_coli", (
                f"Unexpected organism {result['organism']} for gene {result['gene']}"
            )

    def test_cai_values_in_valid_range(self, benchmark_results):
        """All computed CAI values are in [0, 1]."""
        for result in benchmark_results["per_gene_results"]:
            kaz = result["kazusa_cai"]
            sl = result["sharp_li_cai"]
            if not math.isnan(kaz):
                assert 0.0 <= kaz <= 1.0, (
                    f"Kazusa CAI for {result['gene']} out of range: {kaz}"
                )
            if not math.isnan(sl):
                assert 0.0 <= sl <= 1.0, (
                    f"Sharp-Li CAI for {result['gene']} out of range: {sl}"
                )

    def test_published_cai_values_match_database(self, benchmark_results):
        """Published CAI values in the results match PUBLISHED_CAI_VALUES."""
        for result in benchmark_results["per_gene_results"]:
            key = (result["gene"], result["organism"])
            assert key in PUBLISHED_CAI_VALUES, (
                f"Gene {result['gene']}/{result['organism']} not in PUBLISHED_CAI_VALUES"
            )
            expected = PUBLISHED_CAI_VALUES[key]["expected_cai"]
            assert result["published_cai"] == expected, (
                f"Published CAI mismatch for {result['gene']}: "
                f"{result['published_cai']} != {expected}"
            )

    def test_errors_are_non_negative(self, benchmark_results):
        """All absolute errors are non-negative."""
        for result in benchmark_results["per_gene_results"]:
            kaz_err = result["kazusa_error"]
            sl_err = result["sharp_li_error"]
            if not math.isnan(kaz_err):
                assert kaz_err >= 0.0, (
                    f"Negative Kazusa error for {result['gene']}: {kaz_err}"
                )
            if not math.isnan(sl_err):
                assert sl_err >= 0.0, (
                    f"Negative Sharp-Li error for {result['gene']}: {sl_err}"
                )

    def test_mean_errors_are_non_negative(self, benchmark_results):
        """Mean errors are non-negative or NaN."""
        mean_kaz = benchmark_results["mean_kazusa_error"]
        mean_sl = benchmark_results["mean_sharp_li_error"]
        if not math.isnan(mean_kaz):
            assert mean_kaz >= 0.0
        if not math.isnan(mean_sl):
            assert mean_sl >= 0.0

    def test_sharp_li_is_closer_is_bool(self, benchmark_results):
        """sharp_li_is_closer is a boolean."""
        assert isinstance(benchmark_results["sharp_li_is_closer"], bool)

    def test_known_ecoli_genes_are_present(self, benchmark_results):
        """Key E. coli genes with DNA sequences are included in the benchmark."""
        genes_benchmarked = {r["gene"] for r in benchmark_results["per_gene_results"]}

        # These genes have DNA sequences in VALIDATION_SEQUENCES
        expected_genes = {"trpA", "recA", "ompA", "groEL"}
        missing = expected_genes - genes_benchmarked
        assert not missing, (
            f"Expected E. coli genes missing from benchmark: {missing}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Sharp-Li CAI values are closer to published values than Kazusa
# ═══════════════════════════════════════════════════════════════════════════════

class TestSharpLiCloserToPublished:
    """Validate that Sharp-Li reference produces CAI closer to published values."""

    @pytest.fixture(scope="class")
    def benchmark_results(self):
        return benchmark_sharp_li_cai()

    def test_sharp_li_closer_on_average(self, benchmark_results):
        """Sharp-Li reference set produces smaller average error than Kazusa.

        The Sharp-Li reference set is reconstructed with modifications that
        approximate the original 24-gene reference set from Sharp & Li (1987).
        These modifications (Ala preference swap + stronger Leu/Ser/Arg bias)
        bring the computed CAI values closer to the published values on average.
        """
        assert benchmark_results["sharp_li_is_closer"], (
            f"Expected Sharp-Li to be closer: "
            f"mean_kazusa_error={benchmark_results['mean_kazusa_error']:.4f}, "
            f"mean_sharp_li_error={benchmark_results['mean_sharp_li_error']:.4f}"
        )

    def test_sharp_li_closer_for_more_genes(self, benchmark_results):
        """Sharp-Li is closer for more individual genes than Kazusa."""
        sl_wins = 0
        kaz_wins = 0
        for r in benchmark_results["per_gene_results"]:
            if math.isnan(r["kazusa_error"]) or math.isnan(r["sharp_li_error"]):
                continue
            if r["sharp_li_error"] < r["kazusa_error"]:
                sl_wins += 1
            elif r["kazusa_error"] < r["sharp_li_error"]:
                kaz_wins += 1

        assert sl_wins >= kaz_wins, (
            f"Expected Sharp-Li to win on at least as many genes as Kazusa: "
            f"Sharp-Li wins={sl_wins}, Kazusa wins={kaz_wins}"
        )

    def test_sharp_li_closer_for_lacz(self, benchmark_results):
        """Sharp-Li should produce a lower CAI for lacZ than Kazusa.

        lacZ has a published CAI of 0.27 (low expression). The Sharp-Li
        reference set's stronger bias should penalise non-optimal codons
        more heavily, producing a CAI closer to the published value.
        """
        lacz_result = None
        for result in benchmark_results["per_gene_results"]:
            if result["gene"] == "lacZ":
                lacz_result = result
                break

        if lacz_result is None:
            pytest.skip("lacZ DNA sequence not available for benchmarking")

        assert lacz_result["sharp_li_error"] < lacz_result["kazusa_error"], (
            f"Expected Sharp-Li error ({lacz_result['sharp_li_error']:.4f}) < "
            f"Kazusa error ({lacz_result['kazusa_error']:.4f}) for lacZ"
        )

    def test_sharp_li_gives_lower_cai_for_lacz(self, benchmark_results):
        """Sharp-Li reference produces a lower CAI for lacZ than Kazusa.

        This validates the key property of the Sharp-Li reference: its
        stronger codon bias correctly reduces the CAI for lowly expressed
        genes that use non-optimal codons.
        """
        lacz_result = None
        for result in benchmark_results["per_gene_results"]:
            if result["gene"] == "lacZ":
                lacz_result = result
                break

        if lacz_result is None:
            pytest.skip("lacZ DNA sequence not available for benchmarking")

        assert lacz_result["sharp_li_cai"] < lacz_result["kazusa_cai"], (
            f"Expected Sharp-Li CAI ({lacz_result['sharp_li_cai']:.4f}) < "
            f"Kazusa CAI ({lacz_result['kazusa_cai']:.4f}) for lacZ"
        )

    def test_sharp_li_reference_has_gcg_preferred_for_ala(self):
        """Sharp-Li reference should have GCG as the preferred Ala codon.

        This is a key modification: in ribosomal protein genes, GCG is the
        most frequent Ala codon, unlike in the broader Kazusa set where
        GCC dominates.
        """
        ala_freqs = SHARP_LI_ECOLI_REFERENCE["A"]
        preferred_ala = max(ala_freqs, key=ala_freqs.get)
        assert preferred_ala == "GCG", (
            f"Expected GCG as preferred Ala codon, got {preferred_ala}"
        )

    def test_kazusa_reference_has_gcc_preferred_for_ala(self):
        """Kazusa reference should have GCC as the preferred Ala codon.

        This confirms that the two reference sets differ in their Ala
        codon preference, which is the primary modification.
        """
        kazusa_ref = load_reference_set("Escherichia_coli")
        ala_freqs = kazusa_ref["A"]
        preferred_ala = max(ala_freqs, key=ala_freqs.get)
        assert preferred_ala == "GCC", (
            f"Expected GCC as preferred Ala codon in Kazusa, got {preferred_ala}"
        )

    def test_sharp_li_has_stronger_bias_for_leu(self):
        """Sharp-Li reference should have stronger Leu codon bias.

        Non-preferred Leu codons should have lower frequencies in the
        Sharp-Li set compared to the Kazusa set.
        """
        kazusa_ref = load_reference_set("Escherichia_coli")
        kazusa_leu = kazusa_ref["L"]
        sharp_li_leu = SHARP_LI_ECOLI_REFERENCE["L"]

        kazusa_max = max(kazusa_leu.values())
        sharp_li_max = max(sharp_li_leu.values())

        # The preferred codon (CTG) should have the same frequency
        assert kazusa_leu["CTG"] == kazusa_max
        assert sharp_li_leu["CTG"] == sharp_li_max

        # Non-preferred codons should have lower frequencies in Sharp-Li
        for codon in kazusa_leu:
            if codon != "CTG":
                assert sharp_li_leu.get(codon, 0) <= kazusa_leu[codon], (
                    f"Leu codon {codon}: Sharp-Li freq ({sharp_li_leu.get(codon, 0)}) "
                    f"should be <= Kazusa freq ({kazusa_leu[codon]})"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. The benchmark report can be printed without errors
# ═══════════════════════════════════════════════════════════════════════════════

class TestBenchmarkReport:
    """Validate that the benchmark report prints without errors."""

    def test_print_report_with_results(self, capsys):
        """print_benchmark_report runs without errors on normal results."""
        results = benchmark_sharp_li_cai()
        print_benchmark_report(results)
        captured = capsys.readouterr()
        assert len(captured.out) > 0, "Report produced no output"
        assert "Sharp-Li vs Kazusa" in captured.out
        assert "Summary" in captured.out

    def test_print_report_with_empty_results(self, capsys):
        """print_benchmark_report handles empty results gracefully."""
        empty_results = {
            "per_gene_results": [],
            "mean_kazusa_error": float("nan"),
            "mean_sharp_li_error": float("nan"),
            "sharp_li_is_closer": False,
        }
        # Should not raise
        print_benchmark_report(empty_results)
        captured = capsys.readouterr()
        assert "No E. coli genes" in captured.out

    def test_print_report_contains_gene_names(self, capsys):
        """The report output contains the gene names from the benchmark."""
        results = benchmark_sharp_li_cai()
        print_benchmark_report(results)
        captured = capsys.readouterr()

        for result in results["per_gene_results"]:
            assert result["gene"] in captured.out, (
                f"Gene {result['gene']} not found in report output"
            )

    def test_print_report_shows_comparison(self, capsys):
        """The report indicates which reference set is closer overall."""
        results = benchmark_sharp_li_cai()
        print_benchmark_report(results)
        captured = capsys.readouterr()

        # Should mention one of the reference sets as closer
        assert "closer" in captured.out.lower() or "CLOSER" in captured.out


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Direct CAI computation tests with both reference sets
# ═══════════════════════════════════════════════════════════════════════════════

class TestDirectCAIComputation:
    """Direct CAI computation tests using both reference sets on known sequences."""

    def test_both_references_produce_valid_cai_for_trpa(self):
        """Both reference sets produce valid CAI for trpA."""
        seq_data = VALIDATION_SEQUENCES.get(("trpA", "Escherichia_coli"))
        if seq_data is None:
            pytest.skip("trpA DNA sequence not available")

        dna = seq_data.get("dna_sequence_full") or seq_data["dna_sequence"]
        kazusa_ref = load_reference_set("Escherichia_coli")

        kazusa_cai = compute_cai_sharp_li(
            dna, kazusa_ref, skip_met=True, min_adaptiveness=0.01,
        )
        sharp_li_cai = compute_cai_sharp_li(
            dna, SHARP_LI_ECOLI_REFERENCE, skip_met=True, min_adaptiveness=0.01,
        )

        assert 0.0 <= kazusa_cai <= 1.0, f"Kazusa CAI out of range: {kazusa_cai}"
        assert 0.0 <= sharp_li_cai <= 1.0, f"Sharp-Li CAI out of range: {sharp_li_cai}"

    def test_both_references_produce_valid_cai_for_groel(self):
        """Both reference sets produce valid CAI for groEL."""
        seq_data = VALIDATION_SEQUENCES.get(("groEL", "Escherichia_coli"))
        if seq_data is None:
            pytest.skip("groEL DNA sequence not available")

        dna = seq_data.get("dna_sequence_full") or seq_data["dna_sequence"]
        kazusa_ref = load_reference_set("Escherichia_coli")

        kazusa_cai = compute_cai_sharp_li(
            dna, kazusa_ref, skip_met=True, min_adaptiveness=0.01,
        )
        sharp_li_cai = compute_cai_sharp_li(
            dna, SHARP_LI_ECOLI_REFERENCE, skip_met=True, min_adaptiveness=0.01,
        )

        assert 0.0 <= kazusa_cai <= 1.0, f"Kazusa CAI out of range: {kazusa_cai}"
        assert 0.0 <= sharp_li_cai <= 1.0, f"Sharp-Li CAI out of range: {sharp_li_cai}"

    def test_sharp_li_gives_lower_cai_for_lacz_direct(self):
        """Sharp-Li gives lower CAI for lacZ than Kazusa (direct computation)."""
        seq_data = VALIDATION_SEQUENCES.get(("lacZ", "Escherichia_coli"))
        if seq_data is None:
            pytest.skip("lacZ DNA sequence not available")

        dna = seq_data.get("dna_sequence_full") or seq_data["dna_sequence"]
        kazusa_ref = load_reference_set("Escherichia_coli")

        kazusa_cai = compute_cai_sharp_li(
            dna, kazusa_ref, skip_met=True, min_adaptiveness=0.01,
        )
        sharp_li_cai = compute_cai_sharp_li(
            dna, SHARP_LI_ECOLI_REFERENCE, skip_met=True, min_adaptiveness=0.01,
        )

        assert sharp_li_cai < kazusa_cai, (
            f"Expected Sharp-Li CAI < Kazusa CAI for lacZ: "
            f"Sharp-Li={sharp_li_cai:.4f}, Kazusa={kazusa_cai:.4f}"
        )

    def test_reference_sets_differ_in_alanine_preference(self):
        """The two reference sets have different preferred Ala codons."""
        kazusa_ref = load_reference_set("Escherichia_coli")

        kazusa_ala_pref = max(kazusa_ref["A"], key=kazusa_ref["A"].get)
        sharp_li_ala_pref = max(SHARP_LI_ECOLI_REFERENCE["A"],
                                key=SHARP_LI_ECOLI_REFERENCE["A"].get)

        assert kazusa_ala_pref != sharp_li_ala_pref, (
            f"Expected different preferred Ala codons: "
            f"Kazusa={kazusa_ala_pref}, Sharp-Li={sharp_li_ala_pref}"
        )
