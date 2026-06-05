"""
BioCompiler Deprecation Tests — splicing.maxent_score() and score_splice_sites()
=================================================================================

Task ID: F5.1

Tests that:
1. Calling maxent_score() emits a DeprecationWarning
2. Calling score_splice_sites() emits a DeprecationWarning
3. maxent_score_v2() produces scores correlated with maxentscan.score_donor()
4. maxent_score_v2() does NOT emit deprecation warnings
"""

import warnings

import pytest

from biocompiler.splicing import maxent_score, maxent_score_v2, score_splice_sites
from biocompiler.maxentscan import score_donor


# ==============================================================================
# 1. maxent_score() emits DeprecationWarning
# ==============================================================================


class TestMaxentScoreDeprecation:
    """Verify that calling maxent_score() emits a DeprecationWarning."""

    def test_maxent_score_emits_deprecation_warning(self):
        """Calling maxent_score() emits a DeprecationWarning."""
        with pytest.warns(DeprecationWarning, match=r"splicing\.maxent_score\(\) is deprecated"):
            maxent_score("CAGGTAAGT")

    def test_maxent_score_warning_message_content(self):
        """DeprecationWarning mentions anti-correlated PWM and maxentscan alternative."""
        with pytest.warns(DeprecationWarning) as record:
            maxent_score("CAGGTAAGT")
        assert len(record) >= 1
        msg = str(record[0].message)
        assert "anti-correlated" in msg
        assert "maxentscan.score_donor" in msg
        assert "v10.0" in msg

    def test_maxent_score_emits_warning_every_call(self):
        """Each call to maxent_score() produces its own DeprecationWarning."""
        with pytest.warns(DeprecationWarning) as record:
            maxent_score("CAGGTAAGT")
            maxent_score("ATGGTCATC")
            maxent_score("AAGT")
        # Each call should produce at least one warning
        assert len(record) >= 3

    def test_maxent_score_still_returns_float(self):
        """Despite deprecation, maxent_score() still returns a valid float."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = maxent_score("CAGGTAAGT")
        assert isinstance(result, float)
        assert result >= 0.0  # simplified PWM always returns non-negative

    def test_maxent_score_short_context_emits_warning(self):
        """Even for short contexts, maxent_score() emits the deprecation warning."""
        with pytest.warns(DeprecationWarning):
            maxent_score("")

    def test_maxent_score_empty_context_emits_warning(self):
        """Empty string context also triggers the deprecation warning."""
        with pytest.warns(DeprecationWarning):
            maxent_score("")


# ==============================================================================
# 2. score_splice_sites() emits DeprecationWarning
# ==============================================================================


class TestScoreSpliceSitesDeprecation:
    """Verify that calling score_splice_sites() emits a DeprecationWarning."""

    def test_score_splice_sites_emits_deprecation_warning(self):
        """Calling score_splice_sites() emits a DeprecationWarning."""
        with pytest.warns(DeprecationWarning, match=r"splicing\.score_splice_sites\(\) is deprecated"):
            score_splice_sites("AAACAGGTAAGTAAAA")

    def test_score_splice_sites_warning_mentions_alternatives(self):
        """DeprecationWarning mentions maxentscan.scan_splice_sites alternative."""
        with pytest.warns(DeprecationWarning) as record:
            score_splice_sites("AAACAGGTAAGTAAAA")
        assert len(record) >= 1
        msg = str(record[0].message)
        assert "maxentscan.scan_splice_sites" in msg
        assert "v10.0" in msg

    def test_score_splice_sites_no_gt_still_emits_warning(self):
        """Even with no GT dinucleotides, the warning is still emitted."""
        with pytest.warns(DeprecationWarning):
            score_splice_sites("AAAAAAAAAA")

    def test_score_splice_sites_still_returns_list(self):
        """Despite deprecation, score_splice_sites() still returns a valid list."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = score_splice_sites("AAACAGGTAAGTAAAA")
        assert isinstance(result, list)


# ==============================================================================
# 3. maxent_score_v2() produces scores correlated with maxentscan.score_donor()
# ==============================================================================


class TestMaxentScoreV2:
    """Verify that maxent_score_v2() correctly delegates to maxentscan."""

    def test_maxent_score_v2_strong_donor_positive(self):
        """Strong consensus donor (CAGGTAAGT) produces a positive score."""
        score = maxent_score_v2("CAGGTAAGT")
        assert score > 0.0, f"Strong donor should score > 0, got {score}"

    def test_maxent_score_v2_correlated_with_score_donor(self):
        """maxent_score_v2() produces the same score as maxentscan.score_donor()."""
        context = "CAGGTAAGT"
        v2_score = maxent_score_v2(context)
        # Build the same sequence that maxent_score_v2 builds internally
        # For a 9-mer with GT at position 3, no padding is needed
        mes_score = score_donor(context, 3)
        assert abs(v2_score - mes_score) < 0.001, (
            f"maxent_score_v2({context!r}) = {v2_score} != "
            f"score_donor({context!r}, 3) = {mes_score}"
        )

    def test_maxent_score_v2_multiple_contexts_correlated(self):
        """Multiple donor contexts show positive correlation with score_donor()."""
        contexts = [
            "CAGGTAAGT",  # strong consensus
            "ATGGTCATC",  # weak
            "TTTGTAAGT",  # T-rich upstream
            "CCGGTAAGT",  # C-rich upstream
        ]
        v2_scores = [maxent_score_v2(ctx) for ctx in contexts]
        mes_scores = [score_donor(ctx, ctx.find("GT")) for ctx in contexts]

        # Correlation: same ranking order
        for i in range(len(contexts)):
            for j in range(i + 1, len(contexts)):
                v2_order = v2_scores[i] > v2_scores[j]
                mes_order = mes_scores[i] > mes_scores[j]
                assert v2_order == mes_order, (
                    f"Ranking mismatch for {contexts[i]!r} vs {contexts[j]!r}: "
                    f"v2={v2_scores[i]:.2f} vs {v2_scores[j]:.2f}, "
                    f"mes={mes_scores[i]:.2f} vs {mes_scores[j]:.2f}"
                )

    def test_maxent_score_v2_strong_outscores_weak(self):
        """Strong consensus donor scores higher than weak donor."""
        strong = maxent_score_v2("CAGGTAAGT")
        weak = maxent_score_v2("ATGGTCATC")
        assert strong > weak, (
            f"Strong donor ({strong:.2f}) should outscore weak ({weak:.2f})"
        )

    def test_maxent_score_v2_no_gt_returns_impossible(self):
        """Context without GT dinucleotide returns -50.0."""
        score = maxent_score_v2("AAAAAAAAA")
        assert score == -50.0, f"Non-GT context should return -50.0, got {score}"

    def test_maxent_score_v2_short_context_returns_impossible(self):
        """Context shorter than 9-mer returns -50.0."""
        score = maxent_score_v2("AAGT")
        assert score == -50.0, f"Short context should return -50.0, got {score}"

    def test_maxent_score_v2_gt_at_different_position(self):
        """GT not at canonical position 3 still gets scored correctly."""
        # GT at position 1 — needs left-padding
        context = "AGTAAGTCC"
        score = maxent_score_v2(context)
        # Should produce a valid score (not -50.0) since GT is present and
        # maxent_score_v2 pads as needed
        # The GT is at position 1, needs 2 chars of left padding
        # Padded: "AAAGTAAGTCC" with GT at position 3
        expected = score_donor("AA" + context, 3)
        assert abs(score - expected) < 0.001, (
            f"maxent_score_v2({context!r}) = {score}, "
            f"expected {expected}"
        )

    def test_maxent_score_v2_returns_float(self):
        """maxent_score_v2 returns a float."""
        score = maxent_score_v2("CAGGTAAGT")
        assert isinstance(score, float)


# ==============================================================================
# 4. maxent_score_v2() does NOT emit deprecation warnings
# ==============================================================================


class TestMaxentScoreV2NoDeprecation:
    """Verify that maxent_score_v2() does NOT emit DeprecationWarning."""

    def test_maxent_score_v2_no_deprecation_warning(self):
        """Calling maxent_score_v2() does NOT emit any DeprecationWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")  # Catch all warnings
            maxent_score_v2("CAGGTAAGT")
        deprecation_warnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        assert len(deprecation_warnings) == 0, (
            f"maxent_score_v2() should not emit DeprecationWarning, "
            f"but got: {[str(w.message) for w in deprecation_warnings]}"
        )

    def test_maxent_score_v2_no_warnings_at_all(self):
        """Calling maxent_score_v2() does not emit any warnings."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            maxent_score_v2("CAGGTAAGT")
            maxent_score_v2("ATGGTCATC")
        assert len(caught) == 0, (
            f"maxent_score_v2() should not emit any warnings, "
            f"but got: {[(str(w.category.__name__), str(w.message)) for w in caught]}"
        )

    def test_maxent_score_v2_no_gt_no_warning(self):
        """Even when returning -50.0 (no GT), no warning is emitted."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            maxent_score_v2("AAAAAAAAA")
        deprecation_warnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        assert len(deprecation_warnings) == 0
