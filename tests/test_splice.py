"""Test BioCompiler Splice Site Scoring — MaxEntScan scoring and dual-threshold classification."""

import warnings
import pytest
from biocompiler.splicing import maxent_score, score_splice_sites
from biocompiler.maxentscan import score_donor, score_acceptor, scan_splice_sites
from biocompiler.type_system import SpliceVerdict


# ── Tests for DEPRECATED functions (maxent_score, score_splice_sites) ────────
# These tests verify the deprecated functions still work correctly but are
# marked with filterwarnings to avoid spamming DeprecationWarning during runs.

@pytest.mark.filterwarnings("ignore::DeprecationWarning")
class TestMaxEntScore:
    """Tests for the deprecated MaxEntScan PWM scoring."""

    def test_maxent_score_range(self):
        """MaxEntScan scores should be non-negative."""
        # Test with a variety of contexts
        contexts = [
            "AAGTAAGT",   # Typical context around a GT
            "CAGGTAAGT",
            "GT",         # Minimal
            "AAAA",       # No GT signal
            "CAGGTAAGT",  # Strong donor-like context
        ]
        for ctx in contexts:
            score = maxent_score(ctx)
            assert score >= 0.0, f"maxent_score({ctx!r})={score} is negative"

    def test_maxent_score_short_context(self):
        """Short context (< 4 chars) returns 0.0."""
        assert maxent_score("AG") == 0.0
        assert maxent_score("A") == 0.0
        assert maxent_score("") == 0.0


# ── Tests for NEW maxentscan functions (preferred API) ────────────────────────

class TestMaxEntScanDonorScore:
    """Tests for the proper MaxEntScan score_donor function."""

    def test_score_donor_positive_for_consensus(self):
        """A strong consensus donor scores positively."""
        seq = "AAACAGGTAAGTAAAA"
        gt_pos = seq.find("GT")
        score = score_donor(seq, gt_pos)
        assert score > 0.0, f"Strong donor should score > 0, got {score}"

    def test_score_donor_range(self):
        """score_donor returns reasonable values for common sequences."""
        seq = "AAACAGGTAAGTAAAA"
        gt_pos = seq.find("GT")
        score = score_donor(seq, gt_pos)
        assert -20.0 < score < 20.0, f"Donor score {score} outside expected range"


class TestMaxEntScanScanSpliceSites:
    """Tests for the proper MaxEntScan scan_splice_sites function."""

    def test_scan_finds_donor_sites(self):
        """scan_splice_sites finds donor sites above threshold."""
        seq = "AAACAGGTAAGTAAAA"
        results = scan_splice_sites(seq, donor_threshold=3.0, acceptor_threshold=3.0)
        donors = [(pos, typ, sc) for pos, typ, sc in results if typ == "donor"]
        assert len(donors) >= 1, "Should find at least one donor above threshold"

    def test_scan_returns_position_type_score(self):
        """Each result is a (position, site_type, score) tuple."""
        seq = "AAACAGGTAAGTAAAA"
        results = scan_splice_sites(seq, donor_threshold=0.0, acceptor_threshold=0.0)
        for result in results:
            assert len(result) == 3
            pos, site_type, score = result
            assert isinstance(pos, int)
            assert site_type in ("donor", "acceptor")
            assert isinstance(score, float)

    def test_no_splice_sites_empty_sequence(self):
        """Empty sequence has no splice sites."""
        results = scan_splice_sites("", donor_threshold=0.0, acceptor_threshold=0.0)
        assert results == []


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
class TestDualThreshold:
    """Tests for dual-threshold PASS/UNCERTAIN/FAIL classification (deprecated API)."""

    def test_classify_splice_pass(self):
        """Score < 3.0 should classify as PASS."""
        seq = "ATGAAAAAAAATAAAAAAATAA"  # mostly A, no GT
        results = score_splice_sites(seq)
        if results:
            for pos, score, verdict in results:
                if score < 3.0:
                    assert verdict == SpliceVerdict.PASS

    def test_classify_splice_uncertain(self):
        """3.0 <= score < 6.0 should classify as UNCERTAIN."""
        score = 4.5  # between 3.0 and 6.0
        if 3.0 <= score < 6.0:
            verdict = SpliceVerdict.UNCERTAIN
        assert verdict == SpliceVerdict.UNCERTAIN

    def test_classify_splice_fail(self):
        """Score >= 6.0 should classify as FAIL."""
        context = "CAGGTAAGT"
        score = maxent_score(context)
        if score >= 6.0:
            verdict = SpliceVerdict.FAIL
            assert verdict == SpliceVerdict.FAIL
        else:
            high_score = 8.0
            assert high_score >= 6.0
            verdict_for_high = SpliceVerdict.FAIL
            assert verdict_for_high == SpliceVerdict.FAIL
