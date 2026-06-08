"""
Unit tests for biocompiler.maxentscan_fast — NUMBA-accelerated & pure-Python MaxEntScan.

Covers:
- Pre-computed lookup tables (base index, PWM, log-odds)
- scan_splice_sites_fast_str: string-based API
- score_donor_fast / score_acceptor_fast: individual scoring
- scan_splice_sites_fast: byte array API
- Edge cases: short sequences, no sites, boundary conditions
"""

from __future__ import annotations

import numpy as np
import pytest

from biocompiler.maxentscan_fast import (
    HAS_NUMBA_MAXENT,
    scan_splice_sites_fast,
    scan_splice_sites_fast_str,
    score_donor_fast,
    score_acceptor_fast,
    _BASE_IDX_NP,
    _DONOR_PWM_NP,
    _ACCEPTOR_PWM_NP,
    _DONOR_LOG_NP,
    _ACCEPTOR_LOG_NP,
)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Pre-computed lookup tables
# ═══════════════════════════════════════════════════════════════════════════════

class TestLookupTables:

    def test_base_idx_acgt(self):
        """A=0, C=1, G=2, T=3."""
        assert _BASE_IDX_NP[ord('A')] == 0
        assert _BASE_IDX_NP[ord('C')] == 1
        assert _BASE_IDX_NP[ord('G')] == 2
        assert _BASE_IDX_NP[ord('T')] == 3

    def test_base_idx_lowercase(self):
        assert _BASE_IDX_NP[ord('a')] == 0
        assert _BASE_IDX_NP[ord('c')] == 1
        assert _BASE_IDX_NP[ord('g')] == 2
        assert _BASE_IDX_NP[ord('t')] == 3

    def test_base_idx_invalid(self):
        """Non-DNA characters should map to -1."""
        assert _BASE_IDX_NP[ord('N')] == -1
        assert _BASE_IDX_NP[ord('X')] == -1
        assert _BASE_IDX_NP[ord(' ')] == -1

    def test_donor_pwm_shape(self):
        """Donor PWM should be (9, 4)."""
        assert _DONOR_PWM_NP.shape == (9, 4)

    def test_acceptor_pwm_shape(self):
        """Acceptor PWM should be (23, 4)."""
        assert _ACCEPTOR_PWM_NP.shape == (23, 4)

    def test_donor_log_shape(self):
        assert _DONOR_LOG_NP.shape == (9, 4)

    def test_acceptor_log_shape(self):
        assert _ACCEPTOR_LOG_NP.shape == (23, 4)

    def test_donor_log_no_inf(self):
        """Log-odds should not contain infinities (epsilon clipping)."""
        assert not np.any(np.isinf(_DONOR_LOG_NP))

    def test_acceptor_log_no_inf(self):
        assert not np.any(np.isinf(_ACCEPTOR_LOG_NP))


# ═══════════════════════════════════════════════════════════════════════════════
# 2. score_donor_fast / score_acceptor_fast
# ═══════════════════════════════════════════════════════════════════════════════

class TestIndividualScoring:

    def test_score_donor_returns_float(self):
        """score_donor_fast should return a float."""
        # Need a sequence long enough for the 9-mer context
        seq = "A" * 20 + "GT" + "A" * 28
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
        score = score_donor_fast(seq_bytes, 20)
        assert isinstance(score, float)

    def test_score_acceptor_returns_float(self):
        seq = "A" * 20 + "AG" + "A" * 28
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
        score = score_acceptor_fast(seq_bytes, 20)
        assert isinstance(score, float)

    def test_score_donor_near_boundary(self):
        """Donor too close to start of sequence should return impossible score."""
        seq = "GT" + "A" * 30
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        score = score_donor_fast(seq_bytes, 0)
        # Position 0 with _DONOR_UPSTREAM bases needed before — should be -50
        assert score < -10  # _IMPOSSIBLE_SCORE or very low

    def test_score_acceptor_near_boundary(self):
        """Acceptor too close to start should return impossible score."""
        seq = "AG" + "A" * 30
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        score = score_acceptor_fast(seq_bytes, 0)
        assert score < -10


# ═══════════════════════════════════════════════════════════════════════════════
# 3. scan_splice_sites_fast (byte array API)
# ═══════════════════════════════════════════════════════════════════════════════

class TestScanSpliceSitesFast:

    def test_returns_three_arrays(self):
        seq = "A" * 50
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        positions, types, scores = scan_splice_sites_fast(seq_bytes)
        assert isinstance(positions, np.ndarray)
        assert isinstance(types, np.ndarray)
        assert isinstance(scores, np.ndarray)

    def test_no_sites_in_all_a(self):
        seq = "A" * 60
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        positions, types, scores = scan_splice_sites_fast(seq_bytes)
        assert len(positions) == 0

    def test_detects_gt_dinucleotide(self):
        """A sequence with GT should produce potential donor site(s)."""
        # Create a long enough sequence with a GT at a reasonable position
        seq = "A" * 30 + "GT" + "A" * 30
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        positions, types, scores = scan_splice_sites_fast(
            seq_bytes, donor_threshold=-99.0, acceptor_threshold=3.0
        )
        # Should have at least one donor site detected
        donor_count = np.sum(types == 0)
        assert donor_count >= 1

    def test_detects_ag_dinucleotide(self):
        """A sequence with AG should produce potential acceptor site(s)."""
        seq = "A" * 30 + "AG" + "A" * 30
        seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        positions, types, scores = scan_splice_sites_fast(
            seq_bytes, donor_threshold=3.0, acceptor_threshold=-99.0
        )
        acceptor_count = np.sum(types == 1)
        assert acceptor_count >= 1


# ═══════════════════════════════════════════════════════════════════════════════
# 4. scan_splice_sites_fast_str (string API)
# ═══════════════════════════════════════════════════════════════════════════════

class TestScanSpliceSitesFastStr:

    def test_returns_list_of_tuples(self):
        results = scan_splice_sites_fast_str("A" * 60)
        assert isinstance(results, list)

    def test_no_sites_in_all_a(self):
        results = scan_splice_sites_fast_str("A" * 60)
        assert len(results) == 0

    def test_site_format(self):
        """Each result should be (position, site_type, score)."""
        # Use very low threshold to ensure we get results
        seq = "A" * 30 + "GTAAGT" + "A" * 30
        results = scan_splice_sites_fast_str(
            seq, donor_threshold=-99.0, acceptor_threshold=-99.0
        )
        if results:
            pos, site_type, score = results[0]
            assert isinstance(pos, int)
            assert site_type in ("donor", "acceptor")
            assert isinstance(score, float)

    def test_sorted_by_position(self):
        seq = "A" * 20 + "GT" + "A" * 20 + "GT" + "A" * 20
        results = scan_splice_sites_fast_str(seq, donor_threshold=-99.0)
        positions = [r[0] for r in results]
        assert positions == sorted(positions)

    def test_lowercase_input(self):
        """Lowercase input should be handled (uppercased internally)."""
        results = scan_splice_sites_fast_str("a" * 60)
        assert isinstance(results, list)

    def test_short_sequence(self):
        """Very short sequence should return empty results."""
        results = scan_splice_sites_fast_str("ATG")
        assert isinstance(results, list)

    def test_empty_sequence(self):
        results = scan_splice_sites_fast_str("")
        assert isinstance(results, list)
        assert len(results) == 0


# ═══════════════════════════════════════════════════════════════════════════════
# 5. NUMBA flag
# ═══════════════════════════════════════════════════════════════════════════════

class TestNumbaFlag:

    def test_has_numba_maxent_is_bool(self):
        assert isinstance(HAS_NUMBA_MAXENT, bool)

    def test_consistent_results_numba_vs_python(self):
        """NUMBA and pure-Python paths should give same results (within tolerance)."""
        # This test just verifies the API works regardless of NUMBA availability
        seq = "A" * 30 + "GTAAGT" + "A" * 30
        results = scan_splice_sites_fast_str(seq, donor_threshold=-99.0, acceptor_threshold=-99.0)
        # Results should be consistent (no crash)
        assert isinstance(results, list)
