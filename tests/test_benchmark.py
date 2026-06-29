"""Comprehensive tests for BioCompiler benchmark module.

Covers:
  1. _compute_cai: geometric-mean CAI with known sequences
  2. _count_gt: GT dinucleotide counting
  3. _count_cpg_ratio: CpG Obs/Exp ratio
  4. _count_restriction_sites: known enzyme site detection
  5. _build_best_codon_sequence: best-CAI codon selection
  6. BenchmarkResult dataclass: construction & field access
  7. BenchmarkReport: construction, pass_rate property
  8. format_benchmark_report_json: valid JSON output
  9. format_benchmark_report_text: text output format
 10. run_structured_benchmarks: returns BenchmarkReport
 11. REFERENCE_GENES: has expected keys (HBB, INS, EGFP)
 12. GENE_PANEL: has expected number of genes
"""

import json
import math

import pytest

from biocompiler.benchmarking.core import (
    GENE_PANEL,
    REFERENCE_GENES,
    BenchmarkReport,
    BenchmarkResult,
    _build_best_codon_sequence,
    _compute_cai,
    _count_cpg_ratio,
    _count_gt,
    _count_restriction_sites,
    format_benchmark_report_json,
    format_benchmark_report_text,
    run_structured_benchmarks,
)
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES


# ────────────────────────────────────────────────────────────────
# 1. _compute_cai
# ────────────────────────────────────────────────────────────────

class TestComputeCai:
    """Tests for _compute_cai (geometric mean CAI)."""

    def test_empty_sequence_returns_zero(self):
        """Empty sequence → CAI 0.0."""
        assert _compute_cai("", CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]) == 0.0

    def test_short_sequence_returns_zero(self):
        """Sequence shorter than 3 nt → CAI 0.0."""
        assert _compute_cai("AT", CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]) == 0.0

    def test_single_codon_known_cai(self):
        """Single codon with known CAI weight returns that weight."""
        cai_table = {"ATG": 1.0}
        result = _compute_cai("ATG", cai_table)
        assert result == pytest.approx(1.0, abs=1e-6)

    def test_single_codon_partial_cai(self):
        """Single codon with partial CAI weight."""
        cai_table = {"ATG": 0.5}
        result = _compute_cai("ATG", cai_table)
        assert result == pytest.approx(0.5, abs=1e-6)

    def test_two_codons_geometric_mean(self):
        """Two codons → geometric mean of their CAI weights."""
        cai_table = {"ATG": 1.0, "GCT": 0.25}
        result = _compute_cai("ATGGCT", cai_table)
        expected = math.exp((math.log(1.0) + math.log(0.25)) / 2)
        assert result == pytest.approx(expected, abs=1e-6)

    def test_unknown_codon_uses_floor(self):
        """Unknown codons get floor value of 0.001 to avoid log(0)."""
        cai_table = {"ATG": 1.0}  # GCT not in table
        result = _compute_cai("ATGGCT", cai_table)
        expected = math.exp((math.log(1.0) + math.log(0.001)) / 2)
        assert result == pytest.approx(expected, abs=1e-4)

    def test_zero_cai_codon_uses_floor(self):
        """Codon with CAI=0 gets floor value of 0.001."""
        cai_table = {"ATG": 1.0, "GCT": 0.0}
        result = _compute_cai("ATGGCT", cai_table)
        expected = math.exp((math.log(1.0) + math.log(0.001)) / 2)
        assert result == pytest.approx(expected, abs=1e-4)

    def test_ecoli_species_table(self):
        """Using real E. coli CAI table gives a value in [0, 1]."""
        seq = "ATGAAAGCTTAA"
        result = _compute_cai(seq, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])
        assert 0.0 < result <= 1.0

    def test_human_species_table(self):
        """Using real human CAI table gives a value in [0, 1]."""
        seq = "ATGAAAGCTTAA"
        result = _compute_cai(seq, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])
        assert 0.0 < result <= 1.0

    def test_all_high_cai_codons(self):
        """Sequence of all highest-CAI codons should approach 1.0."""
        # Build a table where all codons have CAI = 1.0
        cai_table = {"ATG": 1.0, "GCT": 1.0, "TAA": 1.0}
        result = _compute_cai("ATGGCTTAA", cai_table)
        assert result == pytest.approx(1.0, abs=1e-6)


# ────────────────────────────────────────────────────────────────
# 2. _count_gt
# ────────────────────────────────────────────────────────────────

class TestCountGt:
    """Tests for _count_gt (GT dinucleotide counting)."""

    def test_no_gt_returns_zero(self):
        assert _count_gt("ACACAC") == 0

    def test_single_gt(self):
        assert _count_gt("ACGTAC") == 1

    def test_two_gt(self):
        assert _count_gt("ACGTGTAC") == 2

    def test_overlapping_gt_not_counted(self):
        """GTG should count as 1 GT, not 2 (non-overlapping by position)."""
        # Positions: 0-1: GT, 1-2: TG → only 1 GT
        assert _count_gt("GTG") == 1

    def test_empty_sequence(self):
        assert _count_gt("") == 0

    def test_single_base(self):
        assert _count_gt("G") == 0

    def test_all_gt(self):
        """GTTGGT has GT at positions 0 and 3."""
        assert _count_gt("GTTGGT") == 2

    def test_case_sensitive(self):
        """GT counting is case-sensitive — lowercase 'gt' does NOT match."""
        assert _count_gt("acgt") == 0


# ────────────────────────────────────────────────────────────────
# 3. _count_cpg_ratio
# ────────────────────────────────────────────────────────────────

class TestCountCpgRatio:
    """Tests for _count_cpg_ratio (CpG Obs/Exp ratio)."""

    def test_no_cpg_returns_zero(self):
        """No CG dinucleotides → ratio 0."""
        assert _count_cpg_ratio("ATATAT") == 0.0

    def test_empty_sequence_returns_zero(self):
        assert _count_cpg_ratio("") == 0.0

    def test_single_cg(self):
        """Single CG in ACGT: C=1, G=1, CG=1, expected=1*1/4=0.25, ratio=4.0."""
        result = _count_cpg_ratio("ACGT")
        expected = 1.0 / ((1 * 1) / 4)
        assert result == pytest.approx(expected, abs=1e-6)

    def test_known_cpg_ratio(self):
        """Manual computation for 'CGCG': C=2, G=2, CG=2 (pos 0-1, 2-3)."""
        seq = "CGCG"
        c = seq.count("C")  # 2
        g = seq.count("G")  # 2
        cg = 2  # positions 0-1 and 2-3
        expected = cg / ((c * g) / len(seq))
        result = _count_cpg_ratio(seq)
        assert result == pytest.approx(expected, abs=1e-6)

    def test_no_c_returns_zero(self):
        """No C bases → expected = 0 → ratio = 0."""
        assert _count_cpg_ratio("AGTAGTA") == 0.0

    def test_no_g_returns_zero(self):
        """No G bases → expected = 0 → ratio = 0."""
        assert _count_cpg_ratio("ACTACTA") == 0.0


# ────────────────────────────────────────────────────────────────
# 4. _count_restriction_sites
# ────────────────────────────────────────────────────────────────

class TestCountRestrictionSites:
    """Tests for _count_restriction_sites (forward strand enzyme site counting)."""

    def test_no_sites(self):
        """Sequence without any enzyme recognition sites."""
        seq = "ACACACACAC"
        enzymes = ["EcoRI", "BamHI"]
        assert _count_restriction_sites(seq, enzymes) == 0

    def test_single_ecori_site(self):
        """EcoRI site (GAATTC) present once."""
        seq = "ACGGAATTCACG"
        assert _count_restriction_sites(seq, ["EcoRI"]) == 1

    def test_two_ecori_sites(self):
        """EcoRI site present twice."""
        seq = "GAATTCACGAATTC"
        assert _count_restriction_sites(seq, ["EcoRI"]) == 2

    def test_bamhi_site(self):
        """BamHI site (GGATCC) detected."""
        seq = "ACGGATCCTT"
        assert _count_restriction_sites(seq, ["BamHI"]) == 1

    def test_hindiii_site(self):
        """HindIII site (AAGCTT) detected."""
        seq = "ACAAGCTTTT"
        assert _count_restriction_sites(seq, ["HindIII"]) == 1

    def test_xhoi_site(self):
        """XhoI site (CTCGAG) detected."""
        seq = "ACCTCGAGTT"
        assert _count_restriction_sites(seq, ["XhoI"]) == 1

    def test_multiple_enzymes(self):
        """Multiple enzyme sites in one sequence."""
        seq = "GAATTCGGATCC"  # EcoRI + BamHI
        assert _count_restriction_sites(seq, ["EcoRI", "BamHI"]) == 2

    def test_unknown_enzyme_skipped(self):
        """Unknown enzyme name is skipped without error."""
        seq = "GAATTC"
        assert _count_restriction_sites(seq, ["UnknownEnzyme"]) == 0

    def test_empty_sequence(self):
        assert _count_restriction_sites("", ["EcoRI"]) == 0

    def test_empty_enzyme_list(self):
        assert _count_restriction_sites("GAATTC", []) == 0


# ────────────────────────────────────────────────────────────────
# 5. _build_best_codon_sequence
# ────────────────────────────────────────────────────────────────

class TestBuildBestCodonSequence:
    """Tests for _build_best_codon_sequence (highest-CAI codon selection)."""

    def test_single_amino_acid_human(self):
        """M (Met) → ATG (only codon)."""
        result = _build_best_codon_sequence("M", "human")
        assert result == "ATG"

    def test_single_amino_acid_ecoli(self):
        """M (Met) → ATG (only codon)."""
        result = _build_best_codon_sequence("M", "ecoli")
        assert result == "ATG"

    def test_length_matches_protein(self):
        """Output DNA length = 3 × protein length."""
        protein = "MSKGEELFTG"
        result = _build_best_codon_sequence(protein, "human")
        assert len(result) == len(protein) * 3

    def test_all_standard_amino_acids(self):
        """All 20 standard AAs should produce a valid DNA sequence."""
        protein = "ACDEFGHIKLMNPQRSTVWY"
        result = _build_best_codon_sequence(protein, "human")
        # 20 amino acids → 60 nucleotides
        assert len(result) == 20 * 3

    def test_unknown_amino_acid_produces_nnn(self):
        """Unknown amino acid (e.g. 'X') → 'NNN' placeholder."""
        result = _build_best_codon_sequence("X", "human")
        assert result == "NNN"

    def test_default_species_is_human(self):
        """Default species parameter is 'human'."""
        result_human = _build_best_codon_sequence("A", "human")
        result_default = _build_best_codon_sequence("A")
        assert result_human == result_default

    def test_result_only_contains_valid_codons(self):
        """Each triplet should be a valid codon from AA_TO_CODONS or 'NNN'."""
        from biocompiler.shared.constants import AA_TO_CODONS
        all_codons = set()
        for codons in AA_TO_CODONS.values():
            all_codons.update(codons)
        protein = "ACDEFGHIKLMNPQRSTVWY"
        result = _build_best_codon_sequence(protein, "ecoli")
        for i in range(0, len(result), 3):
            codon = result[i:i+3]
            assert codon in all_codons or codon == "NNN", f"Invalid codon: {codon}"


# ────────────────────────────────────────────────────────────────
# 6. BenchmarkResult dataclass
# ────────────────────────────────────────────────────────────────

class TestBenchmarkResult:
    """Tests for BenchmarkResult dataclass construction and field access."""

    def test_construction_required_fields(self):
        """BenchmarkResult with required fields only."""
        r = BenchmarkResult(
            gene_name="HBB",
            test_name="translation_length",
            passed=True,
            expected="protein_length=147",
            actual="protein_length=147",
        )
        assert r.gene_name == "HBB"
        assert r.test_name == "translation_length"
        assert r.passed is True
        assert r.expected == "protein_length=147"
        assert r.actual == "protein_length=147"

    def test_default_details_is_none(self):
        """Details field defaults to None."""
        r = BenchmarkResult(
            gene_name="INS", test_name="gc", passed=True,
            expected="", actual="",
        )
        assert r.details is None

    def test_default_execution_time_is_zero(self):
        """Execution time defaults to 0.0."""
        r = BenchmarkResult(
            gene_name="EGFP", test_name="cai", passed=False,
            expected="", actual="",
        )
        assert r.execution_time_ms == 0.0

    def test_optional_fields_set(self):
        """BenchmarkResult with all optional fields."""
        r = BenchmarkResult(
            gene_name="HBB",
            test_name="type_predicates",
            passed=True,
            expected=">=4 predicates PASS",
            actual="PASS=6, total=8",
            details="NoStopCodons=PASS; NoCpGIsland=PASS",
            execution_time_ms=12.5,
        )
        assert r.details == "NoStopCodons=PASS; NoCpGIsland=PASS"
        assert r.execution_time_ms == 12.5

    def test_field_access_after_construction(self):
        """All fields are accessible via dot notation."""
        r = BenchmarkResult(
            gene_name="HBB", test_name="gc", passed=True,
            expected="GC in [0.35, 0.55]", actual="GC = 0.42",
            details="ok", execution_time_ms=3.2,
        )
        assert r.gene_name == "HBB"
        assert r.test_name == "gc"
        assert r.passed is True
        assert r.expected == "GC in [0.35, 0.55]"
        assert r.actual == "GC = 0.42"
        assert r.details == "ok"
        assert r.execution_time_ms == 3.2


# ────────────────────────────────────────────────────────────────
# 7. BenchmarkReport
# ────────────────────────────────────────────────────────────────

class TestBenchmarkReport:
    """Tests for BenchmarkReport construction and pass_rate property."""

    def test_default_construction(self):
        """Default BenchmarkReport has empty results and zero counts."""
        report = BenchmarkReport()
        assert report.results == []
        assert report.passed == 0
        assert report.failed == 0
        assert report.total_tests == 0
        assert report.pass_rate == 0.0

    def test_construction_with_results(self):
        """BenchmarkReport with explicit results counts."""
        results = [
            BenchmarkResult("HBB", "gc", True, "", ""),
            BenchmarkResult("HBB", "cai", True, "", ""),
            BenchmarkResult("INS", "gc", False, "", ""),
        ]
        report = BenchmarkReport(
            timestamp="2025-01-01T00:00:00Z",
            version="9.0.0",
            total_tests=3,
            passed=2,
            failed=1,
            results=results,
        )
        assert report.passed == 2
        assert report.failed == 1
        assert report.total_tests == 3
        assert len(report.results) == 3

    def test_pass_rate_all_passed(self):
        """pass_rate = 1.0 when all tests pass."""
        results = [
            BenchmarkResult("HBB", "gc", True, "", ""),
            BenchmarkResult("INS", "gc", True, "", ""),
        ]
        report = BenchmarkReport(passed=2, failed=0, results=results)
        assert report.pass_rate == 1.0

    def test_pass_rate_half_passed(self):
        """pass_rate = 0.5 when half tests pass."""
        results = [
            BenchmarkResult("HBB", "gc", True, "", ""),
            BenchmarkResult("INS", "gc", False, "", ""),
        ]
        report = BenchmarkReport(passed=1, failed=1, results=results)
        assert report.pass_rate == pytest.approx(0.5, abs=1e-6)

    def test_pass_rate_no_results(self):
        """pass_rate = 0.0 with no results (avoids division by zero)."""
        report = BenchmarkReport()
        assert report.pass_rate == 0.0

    def test_timestamp_and_version(self):
        """Timestamp and version are stored correctly."""
        report = BenchmarkReport(timestamp="2025-06-01T12:00:00Z", version="1.0.0")
        assert report.timestamp == "2025-06-01T12:00:00Z"
        assert report.version == "1.0.0"

    def test_summary_property(self):
        """Summary dict can be get and set."""
        report = BenchmarkReport(summary={"by_gene": {"HBB": 2}})
        assert report.summary == {"by_gene": {"HBB": 2}}
        report.summary = {"by_gene": {"INS": 1}}
        assert report.summary == {"by_gene": {"INS": 1}}

    def test_total_tests_returns_max(self):
        """total_tests returns max of stored value and len(results)."""
        results = [
            BenchmarkResult("HBB", "gc", True, "", ""),
            BenchmarkResult("INS", "gc", True, "", ""),
        ]
        # When stored total_tests < len(results), total_tests = len(results)
        report = BenchmarkReport(total_tests=0, results=results, passed=2, failed=0)
        assert report.total_tests == 2

    def test_passed_aliases_successful(self):
        """BenchmarkReport.passed is an alias for .successful."""
        results = [
            BenchmarkResult("HBB", "gc", True, "", ""),
        ]
        report = BenchmarkReport(passed=1, failed=0, results=results)
        assert report.passed == report.successful


# ────────────────────────────────────────────────────────────────
# 8. format_benchmark_report_json
# ────────────────────────────────────────────────────────────────

class TestFormatBenchmarkReportJson:
    """Tests for format_benchmark_report_json."""

    def test_valid_json_output(self):
        """Output is valid JSON."""
        results = [
            BenchmarkResult("HBB", "gc", True, "GC in [0.35, 0.55]", "GC = 0.42"),
        ]
        report = BenchmarkReport(
            timestamp="2025-01-01T00:00:00Z",
            version="9.0.0",
            total_tests=1,
            passed=1,
            failed=0,
            results=results,
        )
        output = format_benchmark_report_json(report)
        parsed = json.loads(output)
        assert isinstance(parsed, dict)

    def test_json_contains_required_fields(self):
        """JSON output contains all required top-level keys."""
        results = [
            BenchmarkResult("HBB", "gc", True, "exp", "act"),
        ]
        report = BenchmarkReport(
            timestamp="2025-01-01T00:00:00Z",
            version="9.0.0",
            total_tests=1,
            passed=1,
            failed=0,
            results=results,
        )
        parsed = json.loads(format_benchmark_report_json(report))
        assert "timestamp" in parsed
        assert "version" in parsed
        assert "total_tests" in parsed
        assert "passed" in parsed
        assert "failed" in parsed
        assert "pass_rate" in parsed
        assert "summary" in parsed
        assert "results" in parsed

    def test_json_results_structure(self):
        """Each result in JSON has gene, test, passed, expected, actual, details, time_ms."""
        results = [
            BenchmarkResult(
                "EGFP", "cai", True, "CAI in [0.6, 1.0]", "CAI = 0.85",
                details="ok", execution_time_ms=5.3,
            ),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=1,
            passed=1, failed=0, results=results,
        )
        parsed = json.loads(format_benchmark_report_json(report))
        r = parsed["results"][0]
        assert r["gene"] == "EGFP"
        assert r["test"] == "cai"
        assert r["passed"] is True
        assert r["expected"] == "CAI in [0.6, 1.0]"
        assert r["actual"] == "CAI = 0.85"
        assert r["details"] == "ok"
        assert r["time_ms"] == 5.3

    def test_json_empty_report(self):
        """Empty report produces valid JSON with empty results."""
        report = BenchmarkReport()
        parsed = json.loads(format_benchmark_report_json(report))
        assert parsed["results"] == []
        assert parsed["total_tests"] == 0


# ────────────────────────────────────────────────────────────────
# 9. format_benchmark_report_text
# ────────────────────────────────────────────────────────────────

class TestFormatBenchmarkReportText:
    """Tests for format_benchmark_report_text."""

    def test_output_is_string(self):
        """Text output is a string."""
        report = BenchmarkReport(
            timestamp="2025-01-01", version="9.0", results=[],
        )
        output = format_benchmark_report_text(report)
        assert isinstance(output, str)

    def test_contains_header(self):
        """Text output contains header lines."""
        report = BenchmarkReport(
            timestamp="2025-01-01", version="9.0", results=[],
        )
        output = format_benchmark_report_text(report)
        assert "BioCompiler Benchmark Report" in output
        assert "Version: 9.0" in output
        assert "Timestamp: 2025-01-01" in output

    def test_contains_results_summary(self):
        """Text output contains pass/fail summary."""
        results = [
            BenchmarkResult("HBB", "gc", True, "exp", "act"),
            BenchmarkResult("INS", "cai", False, "exp2", "act2"),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=2,
            passed=1, failed=1, results=results,
        )
        output = format_benchmark_report_text(report)
        assert "1/2 passed" in output

    def test_pass_indicator(self):
        """Passed tests show [PASS] and failed show [FAIL]."""
        results = [
            BenchmarkResult("HBB", "gc", True, "exp", "act"),
            BenchmarkResult("INS", "cai", False, "exp", "act"),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=2,
            passed=1, failed=1, results=results,
        )
        output = format_benchmark_report_text(report)
        assert "[PASS]" in output
        assert "[FAIL]" in output

    def test_gene_and_test_name_in_output(self):
        """Gene name and test name appear in the output."""
        results = [
            BenchmarkResult("HBB", "gc_content_range", True, "exp", "act"),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=1,
            passed=1, failed=0, results=results,
        )
        output = format_benchmark_report_text(report)
        assert "HBB/gc_content_range" in output

    def test_details_included_when_present(self):
        """Details line is included when BenchmarkResult.details is set."""
        results = [
            BenchmarkResult("HBB", "gc", True, "exp", "act", details="some detail"),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=1,
            passed=1, failed=0, results=results,
        )
        output = format_benchmark_report_text(report)
        assert "some detail" in output
        assert "Details:" in output

    def test_details_omitted_when_none(self):
        """Details line is omitted when BenchmarkResult.details is None."""
        results = [
            BenchmarkResult("HBB", "gc", True, "exp", "act", details=None),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=1,
            passed=1, failed=0, results=results,
        )
        output = format_benchmark_report_text(report)
        assert "Details:" not in output

    def test_expected_and_actual_in_output(self):
        """Expected and actual values appear in the output."""
        results = [
            BenchmarkResult("HBB", "gc", True, "GC in [0.35, 0.55]", "GC = 0.42"),
        ]
        report = BenchmarkReport(
            timestamp="", version="1.0", total_tests=1,
            passed=1, failed=0, results=results,
        )
        output = format_benchmark_report_text(report)
        assert "GC in [0.35, 0.55]" in output
        assert "GC = 0.42" in output


# ────────────────────────────────────────────────────────────────
# 10. run_structured_benchmarks
# ────────────────────────────────────────────────────────────────

class TestRunStructuredBenchmarks:
    """Tests for run_structured_benchmarks."""

    def test_returns_benchmark_report(self):
        """run_structured_benchmarks returns a BenchmarkReport instance."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        assert isinstance(report, BenchmarkReport)

    def test_report_has_results(self):
        """Report contains at least one result per gene."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        assert len(report.results) >= 1

    def test_report_has_timestamp(self):
        """Report has a non-empty timestamp."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        assert report.timestamp != ""

    def test_report_has_version(self):
        """Report has a non-empty version."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        assert report.version != ""

    def test_report_has_summary(self):
        """Report has a summary dict."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        assert isinstance(report.summary, dict)

    def test_single_gene(self):
        """Running with a single gene produces correct number of tests."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        # At least 3 tests: translation_length, gc_content_range, cai_range
        assert len(report.results) >= 3

    def test_all_default_genes(self):
        """Running without gene_names uses all REFERENCE_GENES."""
        report = run_structured_benchmarks()
        # 3 genes × (at least 3 tests each) = at least 9 results
        assert len(report.results) >= 9

    def test_invalid_gene_name_skipped(self):
        """Invalid gene name is silently skipped."""
        report = run_structured_benchmarks(gene_names=["NONEXISTENT"])
        assert len(report.results) == 0

    def test_results_have_correct_gene_name(self):
        """All results have the correct gene_name."""
        report = run_structured_benchmarks(gene_names=["INS"])
        for r in report.results:
            assert r.gene_name == "INS"

    def test_pass_rate_is_valid(self):
        """pass_rate is between 0 and 1."""
        report = run_structured_benchmarks(gene_names=["HBB"])
        assert 0.0 <= report.pass_rate <= 1.0


# ────────────────────────────────────────────────────────────────
# 11. REFERENCE_GENES
# ────────────────────────────────────────────────────────────────

class TestReferenceGenes:
    """Tests for REFERENCE_GENES constant."""

    def test_has_hbb_key(self):
        assert "HBB" in REFERENCE_GENES

    def test_has_ins_key(self):
        assert "INS" in REFERENCE_GENES

    def test_has_egfp_key(self):
        assert "EGFP" in REFERENCE_GENES

    def test_exactly_three_keys(self):
        assert len(REFERENCE_GENES) == 3

    def test_hbb_has_required_fields(self):
        hbb = REFERENCE_GENES["HBB"]
        assert "description" in hbb
        assert "organism" in hbb
        assert "exon_boundaries" in hbb
        assert "known_protein_length" in hbb
        assert "expected_gc_range" in hbb
        assert "expected_cai_range" in hbb
        assert "known_splice_events" in hbb
        assert "pre_mrna" in hbb

    def test_hbb_description(self):
        assert "Beta-Globin" in REFERENCE_GENES["HBB"]["description"]

    def test_ins_description(self):
        assert "Insulin" in REFERENCE_GENES["INS"]["description"]

    def test_egfp_description(self):
        assert "EGFP" in REFERENCE_GENES["EGFP"]["description"] or "Green" in REFERENCE_GENES["EGFP"]["description"]

    def test_exon_boundaries_are_tuples(self):
        for gene_name, gene_data in REFERENCE_GENES.items():
            boundaries = gene_data["exon_boundaries"]
            assert isinstance(boundaries, list)
            for start, end in boundaries:
                assert isinstance(start, int)
                assert isinstance(end, int)
                assert end > start

    def test_gc_range_is_tuple(self):
        for gene_name, gene_data in REFERENCE_GENES.items():
            lo, hi = gene_data["expected_gc_range"]
            assert 0.0 <= lo <= 1.0
            assert 0.0 <= hi <= 1.0
            assert lo < hi

    def test_pre_mrna_not_empty(self):
        for gene_name, gene_data in REFERENCE_GENES.items():
            assert len(gene_data["pre_mrna"]) > 0


# ────────────────────────────────────────────────────────────────
# 12. GENE_PANEL
# ────────────────────────────────────────────────────────────────

class TestGenePanel:
    """Tests for GENE_PANEL constant."""

    def test_has_expected_number_of_genes(self):
        """GENE_PANEL should have 12 genes."""
        assert len(GENE_PANEL) == 12

    def test_contains_hbb(self):
        assert "HBB" in GENE_PANEL

    def test_contains_ins(self):
        assert "INS" in GENE_PANEL

    def test_contains_egfp(self):
        assert "EGFP" in GENE_PANEL

    def test_contains_mcherry(self):
        assert "mCherry" in GENE_PANEL

    def test_contains_cas9_frag(self):
        assert "Cas9_frag" in GENE_PANEL

    def test_contains_lacz_frag(self):
        assert "LacZ_frag" in GENE_PANEL

    def test_each_entry_is_tuple(self):
        """Each GENE_PANEL entry is (protein_sequence, organism) tuple."""
        for gene_name, (protein, organism) in GENE_PANEL.items():
            assert isinstance(protein, str)
            assert isinstance(organism, str)
            assert len(protein) > 0

    def test_protein_sequences_are_valid(self):
        """All protein sequences contain only standard amino acids."""
        standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for gene_name, (protein, organism) in GENE_PANEL.items():
            invalid = set(protein) - standard_aa
            assert len(invalid) == 0, f"{gene_name} has invalid AAs: {invalid}"

    def test_organisms_are_known(self):
        """All organisms in GENE_PANEL are in the ORGANISM_TO_SPECIES mapping."""
        from biocompiler.benchmarking.core import ORGANISM_TO_SPECIES
        for gene_name, (protein, organism) in GENE_PANEL.items():
            assert organism in ORGANISM_TO_SPECIES, f"{gene_name} has unknown organism: {organism}"
