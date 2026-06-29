"""
BioCompiler Splice Scoring Validation Tests
=============================================

Tests for the splice scoring validation module that demonstrates the
anti-correlation between the deprecated splicing.maxent_score() and
the proper maxentscan.score_donor()/score_acceptor() implementations.
"""

from __future__ import annotations

import io
import sys
import warnings

import pytest

from biocompiler.benchmarking.splice_scoring_validation import (
    SpliceValidationResult,
    validate_splice_scoring,
    print_splice_validation_report,
    _pearson_correlation,
)
from biocompiler.sequence.maxentscan import score_donor, score_acceptor


# ==============================================================================
# 1. validate_splice_scoring() runs without errors
# ==============================================================================


class TestValidateSpliceScoringRuns:
    """Verify that validate_splice_scoring() executes without errors."""

    def test_validate_returns_result(self):
        """validate_splice_scoring() returns a SpliceValidationResult."""
        result = validate_splice_scoring()
        assert isinstance(result, SpliceValidationResult)

    def test_validate_returns_test_sites(self):
        """Result contains a non-empty list of test sites."""
        result = validate_splice_scoring()
        assert isinstance(result.test_sites, list)
        assert len(result.test_sites) > 0

    def test_each_test_site_has_required_keys(self):
        """Each test site dict has all required keys."""
        required_keys = {
            "sequence",
            "site_type",
            "expected_strength",
            "deprecated_score",
            "correct_score",
        }
        result = validate_splice_scoring()
        for site in result.test_sites:
            assert required_keys.issubset(site.keys()), (
                f"Missing keys: {required_keys - site.keys()} in site {site}"
            )

    def test_each_test_site_has_valid_type(self):
        """Each test site has a valid site_type."""
        valid_types = {"donor", "acceptor"}
        result = validate_splice_scoring()
        for site in result.test_sites:
            assert site["site_type"] in valid_types, (
                f"Invalid site_type: {site['site_type']}"
            )

    def test_each_test_site_has_valid_strength(self):
        """Each test site has a valid expected_strength."""
        valid_strengths = {"strong", "moderate", "weak"}
        result = validate_splice_scoring()
        for site in result.test_sites:
            assert site["expected_strength"] in valid_strengths, (
                f"Invalid expected_strength: {site['expected_strength']}"
            )

    def test_scores_are_floats(self):
        """All scores are float type."""
        result = validate_splice_scoring()
        for site in result.test_sites:
            assert isinstance(site["deprecated_score"], float), (
                f"deprecated_score is not float: {type(site['deprecated_score'])}"
            )
            assert isinstance(site["correct_score"], float), (
                f"correct_score is not float: {type(site['correct_score'])}"
            )

    def test_result_has_required_fields(self):
        """SpliceValidationResult has all required fields."""
        result = validate_splice_scoring()
        assert isinstance(result.correlation, float)
        assert isinstance(result.is_anti_correlated, bool)
        assert isinstance(result.correct_scores_rank_order, bool)
        assert isinstance(result.deprecated_scores_rank_order, bool)

    def test_has_both_donor_and_acceptor_sites(self):
        """Test sites include both donor and acceptor types."""
        result = validate_splice_scoring()
        site_types = {s["site_type"] for s in result.test_sites}
        assert "donor" in site_types
        assert "acceptor" in site_types

    def test_has_both_strong_and_weak_sites(self):
        """Test sites include both strong and weak strengths."""
        result = validate_splice_scoring()
        strengths = {s["expected_strength"] for s in result.test_sites}
        assert "strong" in strengths
        assert "weak" in strengths


# ==============================================================================
# 2. maxentscan scores are positively correlated with known site strength
# ==============================================================================


class TestMaxentScanCorrelatedWithStrength:
    """Verify that proper maxentscan scores correlate positively with
    known splice site strength (strong > moderate > weak)."""

    def test_strong_donors_score_higher_than_weak(self):
        """Strong donor sites score higher than weak donor sites in maxentscan."""
        result = validate_splice_scoring()
        strong_donors = [
            s["correct_score"]
            for s in result.test_sites
            if s["site_type"] == "donor" and s["expected_strength"] == "strong"
        ]
        weak_donors = [
            s["correct_score"]
            for s in result.test_sites
            if s["site_type"] == "donor" and s["expected_strength"] == "weak"
        ]
        # Exclude impossible scores
        strong_valid = [s for s in strong_donors if s > -50.0]
        weak_valid = [s for s in weak_donors if s > -50.0]

        assert len(strong_valid) > 0, "No valid strong donor scores"
        assert len(weak_valid) > 0, "No valid weak donor scores"
        assert min(strong_valid) > max(weak_valid), (
            f"Strong donors ({strong_valid}) should outscore weak ({weak_valid})"
        )

    def test_strong_acceptors_score_higher_than_weak(self):
        """Strong acceptor sites score higher than weak acceptor sites."""
        result = validate_splice_scoring()
        strong_acc = [
            s["correct_score"]
            for s in result.test_sites
            if s["site_type"] == "acceptor" and s["expected_strength"] == "strong"
        ]
        weak_acc = [
            s["correct_score"]
            for s in result.test_sites
            if s["site_type"] == "acceptor" and s["expected_strength"] == "weak"
        ]
        strong_valid = [s for s in strong_acc if s > -50.0]
        weak_valid = [s for s in weak_acc if s > -50.0]

        assert len(strong_valid) > 0, "No valid strong acceptor scores"
        assert len(weak_valid) > 0, "No valid weak acceptor scores"
        assert min(strong_valid) > max(weak_valid), (
            f"Strong acceptors ({strong_valid}) should outscore weak ({weak_valid})"
        )

    def test_correct_rank_order_passes(self):
        """The correct_scores_rank_order field is True (strong > weak)."""
        result = validate_splice_scoring()
        assert result.correct_scores_rank_order, (
            "Proper maxentscan should rank strong sites above weak sites"
        )

    def test_strong_donor_scores_positive(self):
        """Strong canonical donor sites score positively (> 0) with maxentscan."""
        result = validate_splice_scoring()
        strong_donors = [
            s["correct_score"]
            for s in result.test_sites
            if s["site_type"] == "donor" and s["expected_strength"] == "strong"
        ]
        for score in strong_donors:
            assert score > 0.0, f"Strong donor should score > 0, got {score}"

    def test_strong_acceptor_scores_positive(self):
        """Strong canonical acceptor sites score positively (> 0) with maxentscan."""
        result = validate_splice_scoring()
        strong_acc = [
            s["correct_score"]
            for s in result.test_sites
            if s["site_type"] == "acceptor" and s["expected_strength"] == "strong"
        ]
        for score in strong_acc:
            if score > -50.0:  # Skip impossible scores
                assert score > 0.0, (
                    f"Strong acceptor should score > 0, got {score}"
                )


# ==============================================================================
# 3. Correlation between deprecated and correct methods is computed
# ==============================================================================


class TestCorrelationComputation:
    """Verify that the correlation between methods is properly computed."""

    def test_correlation_is_float(self):
        """Correlation is a float."""
        result = validate_splice_scoring()
        assert isinstance(result.correlation, float)

    def test_correlation_is_in_valid_range(self):
        """Correlation is between -1.0 and 1.0."""
        result = validate_splice_scoring()
        assert -1.0 <= result.correlation <= 1.0, (
            f"Correlation {result.correlation} outside [-1, 1]"
        )

    def test_anti_correlated_flag_matches_correlation(self):
        """is_anti_correlated is True iff correlation < 0."""
        result = validate_splice_scoring()
        if result.correlation < 0:
            assert result.is_anti_correlated
        else:
            assert not result.is_anti_correlated

    def test_deprecated_rank_order_is_wrong(self):
        """The deprecated method fails to rank strong > weak (as expected)."""
        result = validate_splice_scoring()
        # The deprecated PWM gives higher scores to sequences with A at GT
        # positions, which are NOT strong splice sites. So it should NOT
        # correctly rank strong > weak.
        assert not result.deprecated_scores_rank_order, (
            "Deprecated maxent_score should NOT correctly rank "
            "strong > weak splice sites (it is anti-correlated)"
        )

    def test_pearson_correlation_helper(self):
        """_pearson_correlation computes correct values for known inputs."""
        # Perfect positive correlation
        assert _pearson_correlation([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]) == pytest.approx(
            1.0, abs=1e-10
        )
        # Perfect negative correlation
        assert _pearson_correlation([1.0, 2.0, 3.0], [6.0, 4.0, 2.0]) == pytest.approx(
            -1.0, abs=1e-10
        )
        # No correlation (constant x)
        assert _pearson_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == 0.0
        # No correlation (constant y)
        assert _pearson_correlation([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]) == 0.0
        # Too few points
        assert _pearson_correlation([1.0], [2.0]) == 0.0
        assert _pearson_correlation([], []) == 0.0

    def test_anti_correlation_detected(self):
        """The validation detects anti-correlation between methods."""
        result = validate_splice_scoring()
        # The key finding: deprecated maxent_score is anti-correlated
        # with proper MaxEntScan scoring
        assert result.is_anti_correlated, (
            f"Anti-correlation expected but correlation = {result.correlation:.4f}"
        )


# ==============================================================================
# 4. Report printing function
# ==============================================================================


class TestReportPrinting:
    """Verify that print_splice_validation_report produces formatted output."""

    def test_report_prints_without_error(self):
        """print_splice_validation_report runs without raising."""
        result = validate_splice_scoring()
        # Capture stdout to avoid cluttering test output
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        assert len(output) > 0, "Report should produce output"

    def test_report_contains_table(self):
        """Report includes a comparison table with header and data."""
        result = validate_splice_scoring()
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        # Check for table header elements
        assert "Sequence" in output
        assert "Type" in output
        assert "Strength" in output
        assert "Deprecated" in output
        assert "Correct" in output

    def test_report_contains_correlation(self):
        """Report includes the Pearson correlation value."""
        result = validate_splice_scoring()
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        assert "Pearson correlation" in output
        assert "anti-correlated" in output.lower()

    def test_report_contains_rank_order(self):
        """Report includes rank order check results."""
        result = validate_splice_scoring()
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        assert "rank order" in output.lower()

    def test_report_contains_recommendation(self):
        """Report includes a recommendation to use maxentscan."""
        result = validate_splice_scoring()
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        assert "RECOMMENDATION" in output
        assert "maxentscan" in output

    def test_report_shows_anti_correlation_when_detected(self):
        """Report highlights anti-correlation finding when detected."""
        result = validate_splice_scoring()
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        if result.is_anti_correlated:
            assert "ANTI-CORRELATION DETECTED" in output
            assert "Root cause" in output

    def test_report_includes_test_sequences(self):
        """Report shows the actual test sequences in the table."""
        result = validate_splice_scoring()
        captured = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured
            print_splice_validation_report(result)
        finally:
            sys.stdout = old_stdout

        output = captured.getvalue()
        # The canonical donor sequence should appear in the table
        assert "CAGGTGAGT" in output


# ==============================================================================
# 5. Integration: end-to-end scoring consistency
# ==============================================================================


class TestScoringIntegration:
    """Integration tests verifying that the validation module correctly
    wraps the underlying scoring functions."""

    def test_canonical_donor_correct_score_matches_direct_call(self):
        """The correct_score for the canonical donor matches a direct
        call to score_donor()."""
        result = validate_splice_scoring()
        # Find the canonical donor site
        canonical = [
            s
            for s in result.test_sites
            if s["sequence"] == "CAGGTGAGT" and s["site_type"] == "donor"
        ]
        assert len(canonical) == 1, "Expected exactly one canonical donor site"
        site = canonical[0]

        # Direct call to score_donor
        direct_score = score_donor("CAGGTGAGT", 3)
        assert site["correct_score"] == pytest.approx(
            direct_score, abs=0.1
        ), (
            f"Validation score {site['correct_score']} != "
            f"direct score_donor {direct_score}"
        )

    def test_deprecated_score_matches_direct_call(self):
        """The deprecated_score matches a direct call to maxent_score()
        (with warning suppression)."""
        from biocompiler.sequence.splicing import maxent_score

        result = validate_splice_scoring()
        for site in result.test_sites:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                direct = maxent_score(site["sequence"])
            assert site["deprecated_score"] == pytest.approx(
                direct, abs=1e-10
            ), (
                f"Deprecated score mismatch for {site['sequence']!r}: "
                f"{site['deprecated_score']} != {direct}"
            )

    def test_validation_result_is_deterministic(self):
        """Running validate_splice_scoring() twice produces identical results."""
        result1 = validate_splice_scoring()
        result2 = validate_splice_scoring()

        assert len(result1.test_sites) == len(result2.test_sites)
        for s1, s2 in zip(result1.test_sites, result2.test_sites):
            assert s1["deprecated_score"] == pytest.approx(
                s2["deprecated_score"], abs=1e-10
            )
            assert s1["correct_score"] == pytest.approx(
                s2["correct_score"], abs=1e-10
            )

        assert result1.correlation == pytest.approx(
            result2.correlation, abs=1e-10
        )
        assert result1.is_anti_correlated == result2.is_anti_correlated
