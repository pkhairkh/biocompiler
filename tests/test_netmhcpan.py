"""Test BioCompiler NetMHCpan integration — API client, parsing, caching, fallback.

All tests use mocking to avoid calling the real NetMHCpan API.
"""
import pytest
from unittest.mock import patch, MagicMock
import urllib.error

from biocompiler.netmhcpan import (
    NetMHCpanError,
    MHCBindingResult,
    NetMHCpanCache,
    NetMHCpanClient,
    STRONG_BINDER_RANK_THRESHOLD,
    WEAK_BINDER_RANK_THRESHOLD,
    clear_cache,
    is_netmhcpan_available,
    classify_binding_rank,
    parse_netmhcpan_output,
    _rank_to_binding_score,
    _is_mhc_ii_allele,
    _looks_like_netmhcpan_output,
    _extract_result_url,
    _extract_error_message,
)


# ---------------------------------------------------------------------------
# Sample NetMHCpan output for testing
# ---------------------------------------------------------------------------

SAMPLE_NETMHCPAN_OUTPUT = """\
# NetMHCpan version 4.1
# Input is in FASTA format
# Peptide length 9
# Allele HLA-A*02:01
----------------------------------------------------------------------------------------------------
  Pos  HLA         Peptide  Core Of Gp Gl Ip Il Ic  Identity  1-log50k(IC50)  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
    1  HLA-A*02:01  SIINFEKLW     SIINFEKLW  0  0  0  0  0  0   Protein1            0.025           <0.5%  0.12
    2  HLA-A*02:01  IINFEKLWA     IINFEKLWA  0  0  0  0  0  0   Protein1            0.450           WB     1.50
    3  HLA-A*02:01  INFEKLWAA     INFEKLWAA  0  0  0  0  0  0   Protein1            0.015                  8.00
----------------------------------------------------------------------------------------------------
"""

SAMPLE_NETMHCII_OUTPUT = """\
# NetMHCIIpan version 4.0
# Input is in FASTA format
# Allele HLA-DRB1*01:01
----------------------------------------------------------------------------------------------------
  Pos  HLA               Peptide  Identity  Score  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
    1  HLA-DRB1*01:01  AAAAAAAAAAAAAAA   Protein1   0.890  SB     0.25
    2  HLA-DRB1*01:01  AAAAAAAAAAAAAAB   Protein1   0.120         15.0
----------------------------------------------------------------------------------------------------
"""

SAMPLE_ERROR_HTML = """\
<html>
<body>
<h1>Error</h1>
<p>Unrecognized allele: HLA-INVALID*99:99</p>
</body>
</html>
"""

SAMPLE_JOB_SUBMIT_HTML = """\
<html>
<head>
<meta http-equiv="refresh" content="0;url=https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=abc123">
</head>
<body>
<p>Your job has been submitted. Job ID: abc123</p>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# TestNetMHCpanBindingResult
# ---------------------------------------------------------------------------

class TestNetMHCpanBindingResult:
    """Tests for the MHCBindingResult data class."""

    def test_create_result(self):
        """MHCBindingResult can be created with all fields."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=8,
            binding_score=0.85,
            ic50_nm=42.0,
            binding_class="strong_binder",
            rank=0.12,
            method="netmhcpan",
        )
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "SIINFEKL"
        assert result.binding_score == pytest.approx(0.85)
        assert result.ic50_nm == pytest.approx(42.0)
        assert result.binding_class == "strong_binder"
        assert result.rank == pytest.approx(0.12)
        assert result.method == "netmhcpan"

    def test_result_defaults(self):
        """MHCBindingResult has sensible defaults."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="AAAAAAAAA",
            start_position=0,
            end_position=8,
        )
        assert result.binding_score == 0.0
        assert result.ic50_nm is None
        assert result.binding_class == "non_binder"
        assert result.rank is None
        assert result.method == "netmhcpan"
        assert result.anchor_residues == {}
        assert result.anchor_scores == {}


# ---------------------------------------------------------------------------
# TestClassifyBindingRank
# ---------------------------------------------------------------------------

class TestClassifyBindingRank:
    """Tests for the classify_binding_rank function."""

    def test_strong_binder(self):
        """Rank < 0.5% is classified as strong_binder."""
        assert classify_binding_rank(0.01) == "strong_binder"
        assert classify_binding_rank(0.49) == "strong_binder"
        assert classify_binding_rank(0.1) == "strong_binder"

    def test_weak_binder(self):
        """0.5% <= Rank < 2% is classified as weak_binder."""
        assert classify_binding_rank(0.5) == "weak_binder"
        assert classify_binding_rank(1.0) == "weak_binder"
        assert classify_binding_rank(1.99) == "weak_binder"

    def test_non_binder(self):
        """Rank >= 2% is classified as non_binder."""
        assert classify_binding_rank(2.0) == "non_binder"
        assert classify_binding_rank(5.0) == "non_binder"
        assert classify_binding_rank(50.0) == "non_binder"

    def test_threshold_constants(self):
        """Verify the threshold constants match expected values."""
        assert STRONG_BINDER_RANK_THRESHOLD == 0.5
        assert WEAK_BINDER_RANK_THRESHOLD == 2.0


# ---------------------------------------------------------------------------
# TestRankToBindingScore
# ---------------------------------------------------------------------------

class TestRankToBindingScore:
    """Tests for the _rank_to_binding_score conversion function."""

    def test_zero_rank(self):
        """Rank 0 → score 1.0 (strongest possible)."""
        assert _rank_to_binding_score(0) == 1.0

    def test_very_low_rank(self):
        """Very low rank → high binding score."""
        score = _rank_to_binding_score(0.01)
        assert score > 0.9

    def test_moderate_rank(self):
        """Rank ~2% → moderate binding score."""
        score = _rank_to_binding_score(2.0)
        assert 0.3 < score < 1.0

    def test_high_rank(self):
        """High rank → low binding score."""
        score = _rank_to_binding_score(50.0)
        assert score < 0.7  # Monotonically decreasing from 1.0

    def test_score_range(self):
        """Binding score is always in [0, 1]."""
        for rank in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]:
            score = _rank_to_binding_score(rank)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for rank {rank}"

    def test_monotonicity(self):
        """Higher rank → lower binding score (monotonically decreasing)."""
        ranks = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
        scores = [_rank_to_binding_score(r) for r in ranks]
        for i in range(len(scores) - 1):
            assert scores[i] >= scores[i + 1], (
                f"Non-monotonic: rank={ranks[i]} → score={scores[i]}, "
                f"rank={ranks[i+1]} → score={scores[i+1]}"
            )


# ---------------------------------------------------------------------------
# TestIsMHCIIAllele
# ---------------------------------------------------------------------------

class TestIsMHCIIAllele:
    """Tests for the _is_mhc_ii_allele helper function."""

    def test_mhc_ii_alleles(self):
        """Common MHC-II alleles are correctly identified."""
        assert _is_mhc_ii_allele("HLA-DRB1*01:01") is True
        assert _is_mhc_ii_allele("HLA-DQB1*03:01") is True
        assert _is_mhc_ii_allele("HLA-DPB1*04:01") is True
        assert _is_mhc_ii_allele("H2-IAb") is True

    def test_mhc_i_alleles(self):
        """MHC-I alleles are not identified as MHC-II."""
        assert _is_mhc_ii_allele("HLA-A*02:01") is False
        assert _is_mhc_ii_allele("HLA-B*07:02") is False
        assert _is_mhc_ii_allele("H2-Kb") is False
        assert _is_mhc_ii_allele("H2-Db") is False


# ---------------------------------------------------------------------------
# TestParseNetMHCpanOutput
# ---------------------------------------------------------------------------

class TestParseNetMHCpanOutput:
    """Tests for the parse_netmhcpan_output function."""

    def test_parse_mhc_i_output(self):
        """Parse MHC-I output with three binders of different strengths."""
        results = parse_netmhcpan_output(SAMPLE_NETMHCPAN_OUTPUT, "HLA-A*02:01")
        assert len(results) == 3

        # First: strong binder (rank 0.12)
        r0 = results[0]
        assert r0["peptide"] == "SIINFEKLW"
        assert r0["binding_class"] == "strong_binder"
        assert r0["rank"] == pytest.approx(0.12, abs=0.01)

        # Second: weak binder (rank 1.50)
        r1 = results[1]
        assert r1["peptide"] == "IINFEKLWA"
        assert r1["binding_class"] == "weak_binder"
        assert r1["rank"] == pytest.approx(1.50, abs=0.01)

        # Third: non-binder (rank 8.00)
        r2 = results[2]
        assert r2["peptide"] == "INFEKLWAA"
        assert r2["binding_class"] == "non_binder"
        assert r2["rank"] == pytest.approx(8.00, abs=0.01)

    def test_parse_mhc_ii_output(self):
        """Parse MHC-II output."""
        results = parse_netmhcpan_output(SAMPLE_NETMHCII_OUTPUT, "HLA-DRB1*01:01")
        assert len(results) == 2

        r0 = results[0]
        assert r0["binding_class"] == "strong_binder"

    def test_parse_empty_output(self):
        """Empty output returns no results."""
        results = parse_netmhcpan_output("", "HLA-A*02:01")
        assert results == []

    def test_parse_header_only(self):
        """Output with only headers returns no results."""
        header_only = """\
# NetMHCpan version 4.1
# Input is in FASTA format
----------------------------------------------------------------------------------------------------
  Pos  HLA  Peptide  Core  Identity  Score  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
"""
        results = parse_netmhcpan_output(header_only, "HLA-A*02:01")
        assert results == []

    def test_parse_preserves_allele(self):
        """Parsed results contain the correct allele."""
        results = parse_netmhcpan_output(SAMPLE_NETMHCPAN_OUTPUT, "HLA-A*02:01")
        for r in results:
            assert "HLA-A*02:01" in r["allele"] or "A*02:01" in r["allele"]

    def test_position_is_zero_based(self):
        """Positions in parsed output are 0-based (NetMHCpan output is 1-based)."""
        results = parse_netmhcpan_output(SAMPLE_NETMHCPAN_OUTPUT, "HLA-A*02:01")
        if results:
            # First position in NetMHCpan output is 1 → should be 0 in parsed
            assert results[0]["position"] == 0

    def test_ic50_conversion(self):
        """IC50 is computed from the 1-log50k score."""
        results = parse_netmhcpan_output(SAMPLE_NETMHCPAN_OUTPUT, "HLA-A*02:01")
        if results:
            # IC50 should be a positive number
            for r in results:
                assert r["ic50_nm"] > 0


# ---------------------------------------------------------------------------
# TestLooksLikeNetMHCpanOutput
# ---------------------------------------------------------------------------

class TestLooksLikeNetMHCpanOutput:
    """Tests for the _looks_like_netmhcpan_output helper."""

    def test_recognizes_netmhcpan_output(self):
        """Recognizes actual NetMHCpan output."""
        assert _looks_like_netmhcpan_output(SAMPLE_NETMHCPAN_OUTPUT) is True

    def test_recognizes_netmhciipan_output(self):
        """Recognizes NetMHCIIpan output."""
        assert _looks_like_netmhcpan_output(SAMPLE_NETMHCII_OUTPUT) is True

    def test_rejects_empty_string(self):
        """Rejects empty string."""
        assert _looks_like_netmhcpan_output("") is False

    def test_rejects_random_html(self):
        """Rejects generic HTML."""
        assert _looks_like_netmhcpan_output("<html><body>Hello</body></html>") is False


# ---------------------------------------------------------------------------
# TestExtractResultUrl
# ---------------------------------------------------------------------------

class TestExtractResultUrl:
    """Tests for the _extract_result_url helper."""

    def test_extract_from_meta_refresh(self):
        """Extracts URL from meta-refresh tag."""
        url = _extract_result_url(SAMPLE_JOB_SUBMIT_HTML)
        assert url is not None
        assert "jobid=abc123" in url

    def test_extract_from_jobid(self):
        """Extracts job ID and constructs URL."""
        html = '<p>Your job ID is: jobid=xyz789</p>'
        url = _extract_result_url(html)
        assert url is not None
        assert "xyz789" in url

    def test_no_url_in_html(self):
        """Returns None when no result URL is found."""
        html = "<html><body>No redirect here</body></html>"
        url = _extract_result_url(html)
        assert url is None


# ---------------------------------------------------------------------------
# TestExtractErrorMessage
# ---------------------------------------------------------------------------

class TestExtractErrorMessage:
    """Tests for the _extract_error_message helper."""

    def test_extract_allele_error(self):
        """Extracts 'Unrecognized allele' error."""
        msg = _extract_error_message(SAMPLE_ERROR_HTML)
        assert msg is not None
        assert "Unrecognized allele" in msg

    def test_extract_generic_error(self):
        """Extracts a generic error message."""
        html = "<p>Error: Invalid sequence format</p>"
        msg = _extract_error_message(html)
        assert msg is not None
        assert "Invalid sequence format" in msg

    def test_no_error(self):
        """Returns None when no error is found."""
        html = "<html><body>Success</body></html>"
        msg = _extract_error_message(html)
        assert msg is None


# ---------------------------------------------------------------------------
# TestNetMHCpanCache
# ---------------------------------------------------------------------------

class TestNetMHCpanCache:
    """Tests for the NetMHCpanCache class."""

    def test_cache_create(self):
        """Cache can be created with default settings."""
        cache = NetMHCpanCache()
        assert cache.size == 0
        assert cache.hits == 0
        assert cache.misses == 0

    def test_cache_put_get(self):
        """Put a result into cache and get it back."""
        cache = NetMHCpanCache()
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            start_position=0,
            end_position=8,
            binding_score=0.9,
            binding_class="strong_binder",
            rank=0.12,
        )
        cache.put("HLA-A*02:01", "SIINFEKL", result)
        retrieved = cache.get("HLA-A*02:01", "SIINFEKL")
        assert retrieved is not None
        assert retrieved.allele == "HLA-A*02:01"
        assert retrieved.peptide == "SIINFEKL"
        assert retrieved.binding_score == pytest.approx(0.9)

    def test_cache_miss(self):
        """Getting a non-existent key returns None."""
        cache = NetMHCpanCache()
        assert cache.get("HLA-A*02:01", "SIINFEKL") is None
        assert cache.misses == 1

    def test_cache_hit_rate(self):
        """Hit rate tracks cache performance."""
        cache = NetMHCpanCache()
        result = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=8,
        )
        cache.put("HLA-A*02:01", "SIINFEKL", result)

        # Hit
        cache.get("HLA-A*02:01", "SIINFEKL")
        # Miss
        cache.get("HLA-A*02:01", "AAAAAAAAA")

        assert cache.hits == 1
        assert cache.misses == 1
        assert cache.hit_rate == pytest.approx(0.5, abs=0.01)

    def test_cache_eviction(self):
        """When max_size is exceeded, oldest entries are evicted."""
        cache = NetMHCpanCache(max_size=3)

        for i in range(5):
            result = MHCBindingResult(
                allele="HLA-A*02:01", peptide=f"PEPTIDE{i}",
                start_position=0, end_position=8,
            )
            cache.put("HLA-A*02:01", f"PEPTIDE{i}", result)

        # Only the last 3 should remain (FIFO eviction)
        assert cache.get("HLA-A*02:01", "PEPTIDE0") is None
        assert cache.get("HLA-A*02:01", "PEPTIDE1") is None
        assert cache.get("HLA-A*02:01", "PEPTIDE2") is not None
        assert cache.get("HLA-A*02:01", "PEPTIDE3") is not None
        assert cache.get("HLA-A*02:01", "PEPTIDE4") is not None

    def test_cache_clear(self):
        """Clear resets the cache and statistics."""
        cache = NetMHCpanCache()
        result = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=8,
        )
        cache.put("HLA-A*02:01", "SIINFEKL", result)
        cache.get("HLA-A*02:01", "SIINFEKL")
        assert cache.size == 1

        cache.clear()
        assert cache.size == 0
        assert cache.hits == 0
        assert cache.misses == 0

    def test_cache_key_includes_epitope_length(self):
        """Cache keys differentiate by epitope length."""
        cache = NetMHCpanCache()
        r9 = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKLA",
            start_position=0, end_position=8, binding_score=0.5,
        )
        r10 = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKLA",
            start_position=0, end_position=9, binding_score=0.8,
        )
        cache.put("HLA-A*02:01", "SIINFEKLA", r9, epitope_length=9)
        cache.put("HLA-A*02:01", "SIINFEKLA", r10, epitope_length=10)

        got9 = cache.get("HLA-A*02:01", "SIINFEKLA", epitope_length=9)
        got10 = cache.get("HLA-A*02:01", "SIINFEKLA", epitope_length=10)
        assert got9 is not None
        assert got10 is not None
        assert got9.binding_score != got10.binding_score

    def test_module_clear_cache(self):
        """Module-level clear_cache() resets the default cache."""
        cache = _get_cache()
        cache.put("HLA-A*02:01", "TEST", MHCBindingResult(
            allele="HLA-A*02:01", peptide="TEST",
            start_position=0, end_position=3,
        ))
        assert cache.size >= 1

        clear_cache()
        assert cache.size == 0

    def test_put_batch(self):
        """put_batch stores multiple results at once."""
        cache = NetMHCpanCache()
        results = [
            MHCBindingResult(allele="HLA-A*02:01", peptide="PEP1",
                           start_position=0, end_position=2),
            MHCBindingResult(allele="HLA-A*02:01", peptide="PEP2",
                           start_position=3, end_position=5),
        ]
        cache.put_batch(results)
        assert cache.size == 2
        assert cache.get("HLA-A*02:01", "PEP1") is not None
        assert cache.get("HLA-A*02:01", "PEP2") is not None


def _get_cache():
    """Helper to get the module-level default cache."""
    from biocompiler.netmhcpan import _get_default_cache
    return _get_default_cache()


# ---------------------------------------------------------------------------
# TestNetMHCpanClient
# ---------------------------------------------------------------------------

class TestNetMHCpanClient:
    """Tests for the NetMHCpanClient class (with mocked API calls)."""

    def _make_client(self, **kwargs):
        """Create a client with default settings for testing."""
        clear_cache()  # Ensure fresh cache for each client
        defaults = dict(max_retries=1, timeout=5.0, use_cache=True)
        defaults.update(kwargs)
        return NetMHCpanClient(**defaults)

    def test_client_creation(self):
        """Client can be created with default settings."""
        client = self._make_client()
        assert client.api_url is not None
        assert client.timeout == 5.0
        assert client.max_retries == 1
        assert client.use_cache is True

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_success(self, mock_urlopen):
        """predict_mhc_i_binding returns result when API succeeds."""
        # Mock the API response
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCPAN_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        result = client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")

        assert result is not None
        assert result.allele == "HLA-A*02:01"
        assert result.method == "netmhcpan"

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_caches_result(self, mock_urlopen):
        """predict_mhc_i_binding caches results on second call."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCPAN_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        result1 = client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")
        result2 = client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")

        # API should only be called once (second hit is from cache)
        assert mock_urlopen.call_count == 1
        assert result1.allele == result2.allele

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_no_cache(self, mock_urlopen):
        """predict_mhc_i_binding skips cache when use_cache=False."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCPAN_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client(use_cache=False)
        result1 = client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")
        result2 = client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")

        # API should be called twice (no caching)
        assert mock_urlopen.call_count == 2

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_network_error(self, mock_urlopen):
        """predict_mhc_i_binding raises NetMHCpanError on network failure."""
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="API call failed"):
            client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_http_400(self, mock_urlopen):
        """predict_mhc_i_binding raises NetMHCpanError on HTTP 400 error."""
        mock_urlopen.side_effect = urllib.error.HTTPError(
            url="", code=400, msg="Bad Request",
            hdrs={}, fp=MagicMock(read=MagicMock(return_value=b"Invalid allele")),
        )

        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="API call failed"):
            client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_rate_limit_retry(self, mock_urlopen):
        """predict_mhc_i_binding retries on HTTP 429 rate limit."""
        # First call: rate limited
        rate_limit_resp = urllib.error.HTTPError(
            url="", code=429, msg="Too Many Requests",
            hdrs={}, fp=MagicMock(read=MagicMock(return_value=b"rate limited")),
        )
        # Second call: success
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCPAN_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)

        mock_urlopen.side_effect = [rate_limit_resp, mock_resp]

        client = self._make_client(max_retries=3)
        result = client.predict_mhc_i_binding("AAAAAAAV", "HLA-A*02:01")

        assert result is not None
        assert mock_urlopen.call_count == 2

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_empty_peptide(self, mock_urlopen):
        """predict_mhc_i_binding raises NetMHCpanError on empty peptide."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="must not be empty"):
            client.predict_mhc_i_binding("", "HLA-A*02:01")

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_binding_no_results(self, mock_urlopen):
        """predict_mhc_i_binding returns non-binder when API returns no data."""
        # Return output that looks like NetMHCpan but has no data lines
        no_data_output = """\
# NetMHCpan version 4.1
# Input is in FASTA format
# Peptide length 9
# Allele HLA-A*02:01
----------------------------------------------------------------------------------------------------
  Pos  HLA  Peptide  Core Of Gp Gl Ip Il Ic  Identity  1-log50k(IC50)  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
"""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = no_data_output.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        result = client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")

        assert result.binding_class == "non_binder"
        assert result.binding_score == 0.0

    @patch("urllib.request.urlopen")
    def test_predict_mhc_ii_binding_success(self, mock_urlopen):
        """predict_mhc_ii_binding returns result when API succeeds."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCII_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        result = client.predict_mhc_ii_binding(
            "AAAAAAAAAAAAAAA", "HLA-DRB1*01:01",
        )

        assert result is not None
        assert result.allele == "HLA-DRB1*01:01"
        assert result.method == "netmhcpan"

    @patch("urllib.request.urlopen")
    def test_batch_predict(self, mock_urlopen):
        """batch_predict returns results for multiple alleles."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCPAN_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        protein = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 29 AAs
        results = client.batch_predict(
            protein, ["HLA-A*02:01"], epitope_lengths=[9],
        )

        assert isinstance(results, list)
        assert len(results) > 0

    @patch("urllib.request.urlopen")
    def test_batch_predict_mixed_alleles(self, mock_urlopen):
        """batch_predict handles both MHC-I and MHC-II alleles."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_NETMHCPAN_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        protein = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 39 AAs
        results = client.batch_predict(
            protein, ["HLA-A*02:01", "HLA-DRB1*01:01"],
            epitope_lengths=[9],
        )

        assert isinstance(results, list)
        assert len(results) > 0

    @patch("urllib.request.urlopen")
    def test_batch_predict_invalid_protein(self, mock_urlopen):
        """batch_predict raises NetMHCpanError on invalid protein."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="invalid"):
            client.batch_predict("MXG", ["HLA-A*02:01"])


# ---------------------------------------------------------------------------
# TestFallbackToPSSM
# ---------------------------------------------------------------------------

class TestFallbackToPSSM:
    """Tests for fallback from NetMHCpan to PSSM in immunogenicity module."""

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_falls_back_on_api_failure(self, mock_urlopen):
        """predict_mhc_i_binding falls back to PSSM when use_netmhcpan=True but API fails."""
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        from biocompiler.immunogenicity import predict_mhc_i_binding, clear_cache
        clear_cache()

        # Should fall back to PSSM and return results
        protein = "MAGPKWVTFIS"
        results = predict_mhc_i_binding(protein, use_netmhcpan=True)
        assert isinstance(results, list)
        # PSSM should still produce results for known alleles
        # (even if some peptides don't match PSSM length)

    @patch("urllib.request.urlopen")
    def test_predict_mhc_i_pssm_default(self, mock_urlopen):
        """predict_mhc_i_binding uses PSSM by default (use_netmhcpan=False)."""
        # Should not call the API at all
        from biocompiler.immunogenicity import predict_mhc_i_binding, clear_cache
        clear_cache()

        protein = "MAGPKWVTFIS"
        results = predict_mhc_i_binding(protein, use_netmhcpan=False)
        assert isinstance(results, list)
        # API should not have been called
        mock_urlopen.assert_not_called()

    @patch("urllib.request.urlopen")
    def test_predict_mhc_ii_falls_back_on_api_failure(self, mock_urlopen):
        """predict_mhc_ii_binding falls back to PSSM when API fails."""
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        from biocompiler.immunogenicity import predict_mhc_ii_binding, clear_cache
        clear_cache()

        protein = "MAGPKWVTFISLLFLFSSAYS"
        results = predict_mhc_ii_binding(protein, use_netmhcpan=True)
        assert isinstance(results, list)

    @patch("urllib.request.urlopen")
    def test_predict_all_falls_back_on_api_failure(self, mock_urlopen):
        """predict_all falls back to PSSM when API fails."""
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        from biocompiler.immunogenicity import predict_all, clear_cache
        clear_cache()

        protein = "MAGPKWVTFISLLFLFSSAYS"
        result = predict_all(protein, use_netmhcpan=True)
        assert result is not None
        assert isinstance(result.mhc_i_results, list)


# ---------------------------------------------------------------------------
# TestIsNetMHCpanAvailable
# ---------------------------------------------------------------------------

class TestIsNetMHCpanAvailable:
    """Tests for the is_netmhcpan_available function."""

    @patch("urllib.request.urlopen")
    def test_available_when_server_responds(self, mock_urlopen):
        """Returns True when server responds with < 500 status."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        assert is_netmhcpan_available(timeout=5.0) is True

    @patch("urllib.request.urlopen")
    def test_unavailable_when_server_down(self, mock_urlopen):
        """Returns False when server is unreachable."""
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")
        assert is_netmhcpan_available(timeout=1.0) is False

    @patch("urllib.request.urlopen")
    def test_unavailable_on_timeout(self, mock_urlopen):
        """Returns False on timeout."""
        mock_urlopen.side_effect = TimeoutError()
        assert is_netmhcpan_available(timeout=1.0) is False


# ---------------------------------------------------------------------------
# TestNetMHCpanError
# ---------------------------------------------------------------------------

class TestNetMHCpanError:
    """Tests for the NetMHCpanError exception class."""

    def test_error_creation(self):
        """NetMHCpanError can be created with reason and allele."""
        error = NetMHCpanError("API failed", allele="HLA-A*02:01")
        assert "API failed" in str(error)
        assert "HLA-A*02:01" in str(error)

    def test_error_is_immunogenicity_error(self):
        """NetMHCpanError is a subclass of ImmunogenicityError."""
        from biocompiler.exceptions import ImmunogenicityError
        assert issubclass(NetMHCpanError, ImmunogenicityError)

    def test_error_without_allele(self):
        """NetMHCpanError works without allele."""
        error = NetMHCpanError("Network timeout")
        assert "Network timeout" in str(error)

    def test_error_catchable_as_engine_error(self):
        """NetMHCpanError can be caught as EngineError."""
        from biocompiler.exceptions import EngineError
        with pytest.raises(EngineError):
            raise NetMHCpanError("test")
