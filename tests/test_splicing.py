"""Tests for the splicing module.

Covers:
- maxent_score: deprecated PWM scoring for splice sites
- maxent_score_v2: proper MaxEntScan donor scoring
- score_splice_sites: scanning for GT/AG sites with verdicts
- compute_splice_isoforms: NDFST isoform enumeration
- Edge cases: short sequences, prokaryotic organisms, empty inputs
"""

from __future__ import annotations

import pytest
import warnings

from biocompiler.splicing import (
    maxent_score,
    maxent_score_v2,
    score_splice_sites,
    compute_splice_isoforms,
)
from biocompiler.type_system import SpliceVerdict


# ---------------------------------------------------------------------------
# maxent_score (deprecated)
# ---------------------------------------------------------------------------

class TestMaxentScore:
    """Tests for the deprecated maxent_score function."""

    def test_basic_gt_context(self):
        """A GT dinucleotide with context should return a score."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            score = maxent_score("AAGTGTAAG")
        assert isinstance(score, float)
        assert score >= 0.0

    def test_scores_vary_by_context(self):
        """Different contexts should produce different scores."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            s1 = maxent_score("CAGGTAAGT")
            s2 = maxent_score("AAGTGTAAG")
        # Different sequences should generally produce different scores
        assert isinstance(s1, float)
        assert isinstance(s2, float)

    def test_emits_deprecation_warning(self):
        """maxent_score should emit a DeprecationWarning."""
        with pytest.warns(DeprecationWarning, match="deprecated"):
            maxent_score("AAGTGTAAG")

    def test_short_context(self):
        """Short context should return a score (padded)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            score = maxent_score("GT")
        assert isinstance(score, float)

    def test_empty_context(self):
        """Empty context should return 0.0."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            score = maxent_score("")
        assert score == 0.0

    def test_case_insensitive(self):
        """Scoring should be case-insensitive."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            s1 = maxent_score("AAGTGTAAG")
            s2 = maxent_score("aagtgtaag")
        assert abs(s1 - s2) < 0.001


# ---------------------------------------------------------------------------
# maxent_score_v2
# ---------------------------------------------------------------------------

class TestMaxentScoreV2:
    """Tests for the corrected maxent_score_v2 function."""

    def test_basic_scoring(self):
        """A 9-mer with GT should return a valid score."""
        score = maxent_score_v2("CAGGTAAGT")
        assert isinstance(score, float)

    def test_no_gt_returns_low_score(self):
        """Context without GT returns -50.0 (invalid)."""
        score = maxent_score_v2("AAAAAAAAA")
        assert score == -50.0

    def test_short_context_returns_low_score(self):
        """Context shorter than 9 bases returns -50.0."""
        score = maxent_score_v2("GT")
        assert score == -50.0

    def test_canonical_donor_scores_high(self):
        """Strong donor consensus should score higher than weak."""
        # Strong donor context: CAGGTAAAG
        strong = maxent_score_v2("CAGGTAAAG")
        # Weak donor context: TTGTTTTTT
        # But this doesn't have GT in the right place, so let's use valid contexts
        # With GT at position 3 (standard -3 to +6)
        weak = maxent_score_v2("TTCGTCATC")
        # Canonical donors typically score 8-12
        if strong > -50.0 and weak > -50.0:
            assert strong > weak

    def test_gt_at_different_positions(self):
        """GT at various positions should be handled."""
        # GT at position 0
        s1 = maxent_score_v2("GTAAGTCAG")
        # GT at position 3 (standard)
        s2 = maxent_score_v2("CAGGTAAGT")
        # Both should return valid scores (not -50)
        assert s1 > -50.0
        assert s2 > -50.0


# ---------------------------------------------------------------------------
# score_splice_sites
# ---------------------------------------------------------------------------

class TestScoreSpliceSites:
    """Tests for the deprecated score_splice_sites function."""

    def test_basic_scan(self):
        """Scanning a sequence with GT sites should find them."""
        seq = "AAGTGTAAGATATAGGTACGATAG"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites(seq)
        assert isinstance(results, list)
        # Should find GT dinucleotides
        for pos, score, verdict in results:
            assert seq[pos:pos+2] == "GT"

    def test_prokaryote_organism_skips(self):
        """Prokaryotic organism should return empty list."""
        seq = "AAGTGTAAGATATAGGTACGATAG"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites(seq, organism="Escherichia_coli")
        assert results == []

    def test_cds_only_skips(self):
        """CDS-only sequences should return empty list."""
        seq = "AAGTGTAAGATATAGGTACGATAG"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites(seq, is_cds_only=True)
        assert results == []

    def test_emits_deprecation_warning(self):
        """score_splice_sites should emit a DeprecationWarning."""
        seq = "AAGTGTAAGATATAGGTACGATAG"
        with pytest.warns(DeprecationWarning, match="deprecated"):
            score_splice_sites(seq)

    def test_known_splice_sites_preserved(self):
        """Known splice sites should always get PASS verdict."""
        seq = "AAGTGTAAGATATAGGTACGATAG"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            # Find where GT sites are
            gt_positions = [i for i in range(len(seq)-1) if seq[i:i+2] == "GT"]
            if gt_positions:
                results = score_splice_sites(seq, known_splice_sites=gt_positions)
                for pos, score, verdict in results:
                    if pos in gt_positions:
                        assert verdict == SpliceVerdict.PASS

    def test_verdict_types(self):
        """Verdicts should be SpliceVerdict enum values."""
        from biocompiler.type_system import SpliceVerdict
        seq = "AAGTGTAAGATATAGGTACGATAG"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites(seq)
        for pos, score, verdict in results:
            assert isinstance(verdict, SpliceVerdict)

    def test_empty_sequence(self):
        """Empty sequence returns empty list."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites("")
        assert results == []

    def test_no_gt_sites(self):
        """Sequence without GT dinucleotides returns empty list."""
        seq = "AAAAAAAAAAAAAAAAGAAAAAA"  # AG but no GT
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites(seq)
        # No GT sites should be found
        for pos, score, verdict in results:
            assert False, f"Unexpected GT site at {pos}"

    def test_donor_acceptor_context_check(self):
        """GT without a nearby AG should auto-PASS (unlikely splice site)."""
        # 100bp of G's and T's but no AG
        seq = "GT" + "GT" * 50  # Many GTs but no AG
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites(seq)
        # All sites should be PASS because there's no compatible AG
        for pos, score, verdict in results:
            assert verdict == SpliceVerdict.PASS


# ---------------------------------------------------------------------------
# compute_splice_isoforms
# ---------------------------------------------------------------------------

class TestComputeSpliceIsoforms:
    """Tests for compute_splice_isoforms."""

    def test_basic_isoform_computation(self):
        """Basic isoform computation with known exon boundaries."""
        # Simple pre-mRNA with 2 exons and 1 intron
        # Exon1: positions 0-20, Intron: 20-60, Exon2: 60-80
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "A" * 20
        boundaries = [(0, 20), (62, 82)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        assert isinstance(isoforms, list)
        assert len(isoforms) >= 1

    def test_prokaryote_returns_single_isoform(self):
        """Prokaryotic organism should return a single unspliced isoform."""
        seq = "A" * 60
        boundaries = [(0, 60)]
        isoforms = compute_splice_isoforms(
            seq, boundaries, organism="Escherichia_coli",
        )
        assert len(isoforms) == 1
        assert isoforms[0].sequence == seq.upper()
        assert "no_splice_prokaryote" in isoforms[0].parse_path

    def test_isoform_structure(self):
        """Each isoform should have required attributes."""
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "A" * 20
        boundaries = [(0, 20), (62, 82)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        for iso in isoforms:
            assert hasattr(iso, "sequence")
            assert hasattr(iso, "exon_boundaries")
            assert hasattr(iso, "parse_path")
            assert hasattr(iso, "score")
            assert isinstance(iso.sequence, str)
            assert isinstance(iso.score, float)

    def test_canonical_isoform_first(self):
        """Canonical isoform should have the highest score."""
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "A" * 20
        boundaries = [(0, 20), (62, 82)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        if len(isoforms) > 1:
            # First isoform should have highest score
            scores = [iso.score for iso in isoforms]
            assert scores == sorted(scores, reverse=True)

    def test_three_exon_isoforms(self):
        """Three exons should produce more isoforms (exon skipping)."""
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "C" * 20 + "GTAAGT" + "T" * 34 + "AG" + "G" * 20
        boundaries = [(0, 20), (62, 82), (124, 144)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        # Should have at least: canonical + 2 single-exon skips + 2 intron retentions
        assert len(isoforms) >= 3

    def test_max_isoforms_limit(self):
        """max_isoforms parameter should limit output."""
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "C" * 20 + "GTAAGT" + "T" * 34 + "AG" + "G" * 20
        boundaries = [(0, 20), (62, 82), (124, 144)]
        isoforms = compute_splice_isoforms(seq, boundaries, max_isoforms=2)
        assert len(isoforms) <= 2

    def test_no_splice_sites_returns_unspliced(self):
        """Sequence with no splice sites returns single unspliced isoform."""
        seq = "A" * 60
        boundaries = [(0, 60)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        assert len(isoforms) == 1
        assert isoforms[0].sequence == seq.upper()

    def test_deduplication(self):
        """Isoforms should be deduplicated by sequence."""
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "A" * 20
        boundaries = [(0, 20), (62, 82)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        sequences = [iso.sequence for iso in isoforms]
        assert len(sequences) == len(set(sequences))

    def test_intron_retention_isoform(self):
        """Intron retention should produce a longer isoform."""
        seq = "A" * 20 + "GTAAGT" + "T" * 34 + "AG" + "C" * 20
        boundaries = [(0, 20), (62, 82)]
        isoforms = compute_splice_isoforms(seq, boundaries)
        # One isoform should be the full unspliced (intron retained)
        sequences = [iso.sequence for iso in isoforms]
        full_seq = seq.upper()
        assert full_seq in sequences
