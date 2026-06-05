"""
Property-based tests for biocompiler.benchmark using Hypothesis.

Covers three core properties:
  1. _compute_cai returns value in [0, 1] for valid codon sequences
  2. _count_gt count matches manual count for any DNA string
  3. BenchmarkResult fields are consistent (types, invariants)
"""

import math
from dataclasses import fields as dataclass_fields

import pytest
from hypothesis import given, settings, assume
from hypothesis import strategies as st

from biocompiler.benchmark import (
    BenchmarkResult,
    BenchmarkReport,
    _compute_cai,
    _count_gt,
)
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
from biocompiler.type_system import CODON_TABLE


# ────────────────────────────────────────────────────────────
# Shared Strategies
# ────────────────────────────────────────────────────────────

# Single DNA base
dna_base = st.sampled_from("ACGT")

# Arbitrary DNA string (may contain partial codons, lower-case, etc.)
dna_string = st.text(alphabet="ACGT", min_size=0, max_size=300)

# Valid codon sampled from the standard genetic code
valid_codon = st.sampled_from(list(CODON_TABLE.keys()))

# Non-empty valid coding DNA (length always divisible by 3)
valid_coding_seq = st.integers(min_value=1, max_value=50).flatmap(
    lambda n: st.lists(valid_codon, min_size=n, max_size=n).map(
        lambda codons: "".join(codons)
    )
)

# Species key used by _compute_cai via CODON_ADAPTIVENESS_TABLES
species_key = st.sampled_from(list(CODON_ADAPTIVENESS_TABLES.keys()))

# Arbitrary CAI table: codon -> float in (0, 1]
cai_table_entry = st.tuples(valid_codon, st.floats(min_value=0.001, max_value=1.0))
cai_table = st.lists(cai_table_entry, min_size=1, max_size=20).map(
    lambda pairs: dict(pairs)
)

# BenchmarkResult required fields
gene_name_strat = st.sampled_from(["HBB", "INS", "EGFP", "mCherry", "TNF_alpha"])
test_name_strat = st.sampled_from([
    "translation_length", "gc_content_range", "cai_range",
    "type_predicates", "restriction_sites",
])
bool_strat = st.booleans()
text_strat = st.text(min_size=0, max_size=100)


# ────────────────────────────────────────────────────────────
# Property 1: _compute_cai returns value in [0, 1] for valid
#             codon sequences
# ────────────────────────────────────────────────────────────

class TestComputeCaiRange:
    """Property: _compute_cai always returns a value in [0, 1]."""

    @given(seq=valid_coding_seq, species=species_key)
    @settings(max_examples=50, deadline=5000)
    def test_cai_in_unit_interval_real_species(self, seq, species):
        """With a real species CAI table, CAI is in [0, 1]."""
        species_cai = CODON_ADAPTIVENESS_TABLES[species]
        result = _compute_cai(seq, species_cai)
        assert 0.0 <= result <= 1.0, (
            f"CAI {result} out of [0,1] for species={species}, seq_len={len(seq)}"
        )

    @given(seq=valid_coding_seq, table=cai_table)
    @settings(max_examples=50, deadline=5000)
    def test_cai_in_unit_interval_arbitrary_table(self, seq, table):
        """With an arbitrary table where all values are in (0, 1], CAI is in (0, 1]."""
        result = _compute_cai(seq, table)
        assert 0.0 < result <= 1.0, (
            f"CAI {result} out of (0,1] for arbitrary table with {len(table)} entries"
        )

    @given(seq=valid_coding_seq, species=species_key)
    @settings(max_examples=30, deadline=5000)
    def test_cai_geometric_mean_property(self, seq, species):
        """CAI equals the geometric mean of per-codon adaptiveness values."""
        species_cai = CODON_ADAPTIVENESS_TABLES[species]
        result = _compute_cai(seq, species_cai)

        # Manually compute geometric mean
        log_sum = 0.0
        count = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            cai = species_cai.get(codon, 0.0)
            if cai <= 0:
                cai = 0.001
            log_sum += math.log(cai)
            count += 1

        if count == 0:
            expected = 0.0
        else:
            expected = math.exp(log_sum / count)

        assert result == pytest.approx(expected, abs=1e-10), (
            f"CAI {result} != geometric mean {expected}"
        )

    @given(n=st.integers(min_value=1, max_value=10))
    @settings(max_examples=20, deadline=3000)
    def test_all_perfect_cai_gives_one(self, n):
        """A CAI table where every codon maps to 1.0 yields CAI = 1.0."""
        perfect_table = {codon: 1.0 for codon in CODON_TABLE}
        seq = "ATG" * n  # n methionines, all CAI=1.0
        result = _compute_cai(seq, perfect_table)
        assert result == pytest.approx(1.0, abs=1e-10)

    def test_empty_sequence_returns_zero(self):
        """Empty sequence → CAI 0.0 (not in (0,1], but 0 is valid)."""
        assert _compute_cai("", CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]) == 0.0

    def test_short_sequence_returns_zero(self):
        """Sequence shorter than one codon → CAI 0.0."""
        assert _compute_cai("AT", CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]) == 0.0


# ────────────────────────────────────────────────────────────
# Property 2: _count_gt count matches manual count
# ────────────────────────────────────────────────────────────

class TestCountGtManualMatch:
    """Property: _count_gt matches a manual sliding-window count for any string."""

    @staticmethod
    def _manual_gt_count(seq: str) -> int:
        """Naive GT dinucleotide counter — reference implementation."""
        return sum(1 for i in range(len(seq) - 1) if seq[i] == "G" and seq[i + 1] == "T")

    @given(seq=dna_string)
    @settings(max_examples=200, deadline=3000)
    def test_count_matches_manual(self, seq):
        """_count_gt matches the reference implementation for any DNA string."""
        assert _count_gt(seq) == self._manual_gt_count(seq), (
            f"Mismatch for {seq!r}: got {_count_gt(seq)}, expected {self._manual_gt_count(seq)}"
        )

    @given(seq=dna_string)
    @settings(max_examples=100, deadline=3000)
    def test_count_is_non_negative(self, seq):
        """GT count is always non-negative."""
        assert _count_gt(seq) >= 0

    @given(seq=dna_string)
    @settings(max_examples=100, deadline=3000)
    def test_count_upper_bound(self, seq):
        """GT count cannot exceed len(seq) - 1."""
        if len(seq) <= 1:
            assert _count_gt(seq) == 0
        else:
            assert _count_gt(seq) <= len(seq) - 1

    @given(n=st.integers(min_value=0, max_value=20))
    @settings(max_examples=30, deadline=3000)
    def test_repeated_gt_pattern(self, n):
        """Sequence of n 'GT' repeats has exactly n GT dinucleotides."""
        seq = "GT" * n
        assert _count_gt(seq) == n

    @given(n=st.integers(min_value=0, max_value=20))
    @settings(max_examples=30, deadline=3000)
    def test_repeated_gtg_pattern(self, n):
        """Sequence of n 'GTG' repeats: GTG overlaps, each has 1 GT at pos 0-1,
        plus between adjacent repeats the last G + next G don't form GT.
        So: n GTs within codons, plus (n-1) inter-codon boundary checks.
        Between 'GTG' and 'GTG': ...G|G... → no GT. Total = n."""
        seq = "GTG" * n
        assert _count_gt(seq) == self._manual_gt_count(seq)

    def test_empty_string_zero(self):
        """Empty string has 0 GT dinucleotides."""
        assert _count_gt("") == 0

    def test_single_char_zero(self):
        """Single character string has 0 GT dinucleotides."""
        assert _count_gt("G") == 0
        assert _count_gt("T") == 0

    @given(seq=st.text(alphabet="AC", min_size=0, max_size=100))
    @settings(max_examples=50, deadline=3000)
    def test_no_g_or_t_means_zero(self, seq):
        """A sequence without G or without T has zero GT dinucleotides."""
        # seq only has A and C — no G at all, so no GT
        assert _count_gt(seq) == 0

    @given(seq=valid_coding_seq)
    @settings(max_examples=50, deadline=5000)
    def test_valid_coding_seq_count_matches(self, seq):
        """GT count matches manual count for valid coding sequences."""
        assert _count_gt(seq) == self._manual_gt_count(seq)


# ────────────────────────────────────────────────────────────
# Property 3: BenchmarkResult fields are consistent
# ────────────────────────────────────────────────────────────

class TestBenchmarkResultFieldConsistency:
    """Property: BenchmarkResult fields maintain type and value invariants."""

    EXPECTED_FIELDS = {
        "gene_name": str,
        "test_name": str,
        "passed": bool,
        "expected": str,
        "actual": str,
        "details": (str, type(None)),
        "execution_time_ms": (float, int),
    }

    @given(
        gene_name=gene_name_strat,
        test_name=test_name_strat,
        passed=bool_strat,
        expected=text_strat,
        actual=text_strat,
    )
    @settings(max_examples=50, deadline=3000)
    def test_required_fields_types(self, gene_name, test_name, passed, expected, actual):
        """All required fields have the expected types."""
        r = BenchmarkResult(
            gene_name=gene_name,
            test_name=test_name,
            passed=passed,
            expected=expected,
            actual=actual,
        )
        assert isinstance(r.gene_name, str)
        assert isinstance(r.test_name, str)
        assert isinstance(r.passed, bool)
        assert isinstance(r.expected, str)
        assert isinstance(r.actual, str)

    @given(
        gene_name=gene_name_strat,
        test_name=test_name_strat,
        passed=bool_strat,
        expected=text_strat,
        actual=text_strat,
    )
    @settings(max_examples=50, deadline=3000)
    def test_details_defaults_to_none(self, gene_name, test_name, passed, expected, actual):
        """details defaults to None when not provided."""
        r = BenchmarkResult(
            gene_name=gene_name,
            test_name=test_name,
            passed=passed,
            expected=expected,
            actual=actual,
        )
        assert r.details is None

    @given(
        gene_name=gene_name_strat,
        test_name=test_name_strat,
        passed=bool_strat,
        expected=text_strat,
        actual=text_strat,
    )
    @settings(max_examples=50, deadline=3000)
    def test_execution_time_defaults_to_zero(self, gene_name, test_name, passed, expected, actual):
        """execution_time_ms defaults to 0.0 when not provided."""
        r = BenchmarkResult(
            gene_name=gene_name,
            test_name=test_name,
            passed=passed,
            expected=expected,
            actual=actual,
        )
        assert r.execution_time_ms == 0.0

    @given(
        gene_name=gene_name_strat,
        test_name=test_name_strat,
        passed=bool_strat,
        expected=text_strat,
        actual=text_strat,
        details=st.one_of(st.none(), st.text(min_size=0, max_size=200)),
        execution_time_ms=st.floats(min_value=0.0, max_value=10000.0, allow_nan=False, allow_infinity=False),
    )
    @settings(max_examples=80, deadline=3000)
    def test_all_fields_stored_correctly(
        self, gene_name, test_name, passed, expected, actual, details, execution_time_ms
    ):
        """All fields (including optional) are stored exactly as given."""
        r = BenchmarkResult(
            gene_name=gene_name,
            test_name=test_name,
            passed=passed,
            expected=expected,
            actual=actual,
            details=details,
            execution_time_ms=execution_time_ms,
        )
        assert r.gene_name == gene_name
        assert r.test_name == test_name
        assert r.passed == passed
        assert r.expected == expected
        assert r.actual == actual
        assert r.details == details
        assert r.execution_time_ms == execution_time_ms

    @given(
        details=st.one_of(st.none(), st.text(min_size=0, max_size=50)),
    )
    @settings(max_examples=30, deadline=3000)
    def test_details_type_is_correct(self, details):
        """details is always str or None."""
        r = BenchmarkResult(
            gene_name="TEST", test_name="test", passed=True,
            expected="", actual="", details=details,
        )
        assert r.details is None or isinstance(r.details, str)

    @given(
        execution_time_ms=st.floats(min_value=0.0, max_value=100000.0, allow_nan=False, allow_infinity=False),
    )
    @settings(max_examples=30, deadline=3000)
    def test_execution_time_is_non_negative(self, execution_time_ms):
        """execution_time_ms is always non-negative."""
        r = BenchmarkResult(
            gene_name="TEST", test_name="test", passed=True,
            expected="", actual="", execution_time_ms=execution_time_ms,
        )
        assert r.execution_time_ms >= 0.0

    def test_dataclass_has_seven_fields(self):
        """BenchmarkResult has exactly 7 fields."""
        field_names = [f.name for f in dataclass_fields(BenchmarkResult)]
        assert len(field_names) == 7, f"Expected 7 fields, got {len(field_names)}: {field_names}"

    def test_expected_field_names(self):
        """BenchmarkResult has the expected field names."""
        field_names = {f.name for f in dataclass_fields(BenchmarkResult)}
        expected = {"gene_name", "test_name", "passed", "expected", "actual", "details", "execution_time_ms"}
        assert field_names == expected, f"Expected {expected}, got {field_names}"


class TestBenchmarkResultPassFieldConsistency:
    """Property: The 'passed' field is consistently a bool and matches
    semantic expectations."""

    @given(passed=bool_strat)
    @settings(max_examples=20, deadline=3000)
    def test_passed_is_bool(self, passed):
        """passed field is always a bool."""
        r = BenchmarkResult(
            gene_name="X", test_name="t", passed=passed,
            expected="", actual="",
        )
        assert isinstance(r.passed, bool)

    @given(passed=bool_strat)
    @settings(max_examples=20, deadline=3000)
    def test_passed_value_preserved(self, passed):
        """passed field stores the exact value given."""
        r = BenchmarkResult(
            gene_name="X", test_name="t", passed=passed,
            expected="", actual="",
        )
        assert r.passed is passed


class TestBenchmarkReportConsistency:
    """Property: BenchmarkReport fields maintain consistency invariants."""

    @given(
        n_passed=st.integers(min_value=0, max_value=20),
        n_failed=st.integers(min_value=0, max_value=20),
    )
    @settings(max_examples=50, deadline=3000)
    def test_pass_rate_is_in_unit_interval(self, n_passed, n_failed):
        """pass_rate is always in [0.0, 1.0]."""
        results = [
            BenchmarkResult("G", "t", True, "", "")
            for _ in range(n_passed)
        ] + [
            BenchmarkResult("G", "t", False, "", "")
            for _ in range(n_failed)
        ]
        total = n_passed + n_failed
        report = BenchmarkReport(
            passed=n_passed,
            failed=n_failed,
            results=results,
            total_tests=total,
        )
        assert 0.0 <= report.pass_rate <= 1.0

    @given(
        n_passed=st.integers(min_value=0, max_value=20),
        n_failed=st.integers(min_value=0, max_value=20),
    )
    @settings(max_examples=50, deadline=3000)
    def test_pass_rate_equals_passed_over_total(self, n_passed, n_failed):
        """pass_rate = passed / total (or 0 when total=0)."""
        results = [
            BenchmarkResult("G", "t", True, "", "")
            for _ in range(n_passed)
        ] + [
            BenchmarkResult("G", "t", False, "", "")
            for _ in range(n_failed)
        ]
        total = n_passed + n_failed
        report = BenchmarkReport(
            passed=n_passed,
            failed=n_failed,
            results=results,
            total_tests=total,
        )
        if total == 0:
            assert report.pass_rate == 0.0
        else:
            expected_rate = n_passed / total
            assert report.pass_rate == pytest.approx(expected_rate, abs=1e-10)

    @given(
        n_passed=st.integers(min_value=0, max_value=10),
        n_failed=st.integers(min_value=0, max_value=10),
    )
    @settings(max_examples=30, deadline=3000)
    def test_passed_equals_successful(self, n_passed, n_failed):
        """BenchmarkReport.passed is an alias for .successful (from BatchResult)."""
        results = [
            BenchmarkResult("G", "t", True, "", "")
            for _ in range(n_passed)
        ] + [
            BenchmarkResult("G", "t", False, "", "")
            for _ in range(n_failed)
        ]
        report = BenchmarkReport(
            passed=n_passed,
            failed=n_failed,
            results=results,
        )
        assert report.passed == report.successful

    @given(
        n_passed=st.integers(min_value=0, max_value=10),
        n_failed=st.integers(min_value=0, max_value=10),
    )
    @settings(max_examples=30, deadline=3000)
    def test_total_tests_at_least_len_results(self, n_passed, n_failed):
        """total_tests >= len(results)."""
        results = [
            BenchmarkResult("G", "t", True, "", "")
            for _ in range(n_passed)
        ] + [
            BenchmarkResult("G", "t", False, "", "")
            for _ in range(n_failed)
        ]
        report = BenchmarkReport(
            passed=n_passed,
            failed=n_failed,
            results=results,
        )
        assert report.total_tests >= len(report.results)

    def test_default_report_has_zero_pass_rate(self):
        """Default BenchmarkReport has pass_rate 0.0."""
        report = BenchmarkReport()
        assert report.pass_rate == 0.0

    def test_default_report_has_no_results(self):
        """Default BenchmarkReport has no results."""
        report = BenchmarkReport()
        assert report.results == []

    @given(
        summary=st.dictionaries(
            st.text(min_size=1, max_size=10, alphabet="abcdefghijklmnopqrstuvwxyz"),
            st.integers(min_value=0, max_value=100),
            max_size=5,
        )
    )
    @settings(max_examples=30, deadline=3000)
    def test_summary_stored_and_retrieved(self, summary):
        """Summary dict is stored and can be retrieved."""
        report = BenchmarkReport(summary=summary)
        assert report.summary == summary

    @given(
        summary1=st.dictionaries(
            st.text(min_size=1, max_size=5, alphabet="abc"),
            st.integers(min_value=0, max_value=10),
            max_size=3,
        ),
        summary2=st.dictionaries(
            st.text(min_size=1, max_size=5, alphabet="xyz"),
            st.integers(min_value=0, max_value=10),
            max_size=3,
        ),
    )
    @settings(max_examples=20, deadline=3000)
    def test_summary_setter_overwrites(self, summary1, summary2):
        """Setting summary overwrites the previous value."""
        assume(summary1 != summary2)
        report = BenchmarkReport(summary=summary1)
        assert report.summary == summary1
        report.summary = summary2
        assert report.summary == summary2
