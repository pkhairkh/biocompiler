"""Unit tests for biocompiler.netmhcpan — pure-function tests that do NOT require NetMHCpan.

These tests exercise the module's public API, input validation, fallback
behaviour, cache management, and pure functions directly, without requiring
the NetMHCpan web API or any network connectivity.

Test categories
---------------
1. Named constants — correct values, types, and invariants
2. Function existence and signatures — all public symbols exist and are callable
3. Input validation — empty peptides, invalid proteins, length mismatches
4. Fallback behaviour — non-binder result when API returns no data
5. MHCBindingResult dataclass — construction, defaults, field types
6. NetMHCpanCache — put/get/eviction/hit_rate/clear/put_batch/key isolation
7. classify_binding_rank — boundary values for all three classification bands
8. _rank_to_binding_score — range, monotonicity, edge cases
9. parse_netmhcpan_output — header skipping, data parsing, position conversion
10. Helper functions — _is_mhc_ii_allele, _looks_like_netmhcpan_output,
    _extract_result_url, _extract_error_message
11. NetMHCpanError — exception hierarchy and attributes
12. Module-level helpers — clear_cache, _get_default_cache, _get_default_client
"""
from __future__ import annotations

import inspect
from unittest.mock import MagicMock, patch

import pytest

from biocompiler.netmhcpan import (
    DEFAULT_API_URL,
    DEFAULT_MHC_I_EPITOPE_LENGTHS,
    DEFAULT_MHCII_API_URL,
    DEFAULT_TIMEOUT,
    MAX_POLL_ATTEMPTS,
    MAX_RETRIES,
    MHC_II_EPITOPE_LENGTH,
    MHCBindingResult,
    NetMHCpanCache,
    NetMHCpanClient,
    NetMHCpanError,
    POLL_INTERVAL,
    RETRY_BASE_DELAY,
    STRONG_BINDER_RANK_THRESHOLD,
    WEAK_BINDER_RANK_THRESHOLD,
    _extract_error_message,
    _extract_result_url,
    _get_default_cache,
    _get_default_client,
    _is_mhc_ii_allele,
    _looks_like_netmhcpan_output,
    _parse_data_line,
    _rank_to_binding_score,
    batch_predict,
    classify_binding_rank,
    clear_cache,
    is_netmhcpan_available,
    parse_netmhcpan_output,
    predict_mhc_i_binding,
    predict_mhc_ii_binding,
)
from biocompiler.exceptions import EngineError, ImmunogenicityError

# ---------------------------------------------------------------------------
# Sample NetMHCpan output fixtures
# ---------------------------------------------------------------------------

SAMPLE_MHC_I_OUTPUT = """\
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

SAMPLE_MHC_II_OUTPUT = """\
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

SAMPLE_ERROR_HTML = """\
<html>
<body>
<h1>Error</h1>
<p>Unrecognized allele: HLA-INVALID*99:99</p>
</body>
</html>
"""


# ═══════════════════════════════════════════════════════════════════════════
# 1. Named constants
# ═══════════════════════════════════════════════════════════════════════════

class TestNamedConstants:
    """Verify all exported constants have expected values and types."""

    def test_default_api_url_is_dtu(self) -> None:
        """DEFAULT_API_URL points to the DTU webface2 CGI endpoint."""
        assert isinstance(DEFAULT_API_URL, str)
        assert "services.healthtech.dtu.dk" in DEFAULT_API_URL
        assert "webface2" in DEFAULT_API_URL

    def test_default_mhcii_api_url_is_dtu(self) -> None:
        """DEFAULT_MHCII_API_URL points to the DTU webface2 CGI endpoint."""
        assert isinstance(DEFAULT_MHCII_API_URL, str)
        assert "services.healthtech.dtu.dk" in DEFAULT_MHCII_API_URL

    def test_default_timeout_positive(self) -> None:
        """DEFAULT_TIMEOUT must be a positive float."""
        assert isinstance(DEFAULT_TIMEOUT, float)
        assert DEFAULT_TIMEOUT > 0.0

    def test_max_retries_positive_int(self) -> None:
        """MAX_RETRIES must be a positive integer."""
        assert isinstance(MAX_RETRIES, int)
        assert MAX_RETRIES >= 1

    def test_retry_base_delay_positive(self) -> None:
        """RETRY_BASE_DELAY must be positive."""
        assert isinstance(RETRY_BASE_DELAY, float)
        assert RETRY_BASE_DELAY > 0.0

    def test_poll_interval_positive(self) -> None:
        """POLL_INTERVAL must be positive."""
        assert isinstance(POLL_INTERVAL, float)
        assert POLL_INTERVAL > 0.0

    def test_max_poll_attempts_positive(self) -> None:
        """MAX_POLL_ATTEMPTS must be a positive integer."""
        assert isinstance(MAX_POLL_ATTEMPTS, int)
        assert MAX_POLL_ATTEMPTS >= 1

    def test_strong_binder_threshold(self) -> None:
        """STRONG_BINDER_RANK_THRESHOLD is 0.5%."""
        assert isinstance(STRONG_BINDER_RANK_THRESHOLD, float)
        assert STRONG_BINDER_RANK_THRESHOLD == pytest.approx(0.5)

    def test_weak_binder_threshold(self) -> None:
        """WEAK_BINDER_RANK_THRESHOLD is 2.0%."""
        assert isinstance(WEAK_BINDER_RANK_THRESHOLD, float)
        assert WEAK_BINDER_RANK_THRESHOLD == pytest.approx(2.0)

    def test_strong_less_than_weak(self) -> None:
        """Strong binder threshold must be strictly less than weak."""
        assert STRONG_BINDER_RANK_THRESHOLD < WEAK_BINDER_RANK_THRESHOLD

    def test_mhc_ii_epitope_length(self) -> None:
        """MHC_II_EPITOPE_LENGTH is 15 (standard 15-mer)."""
        assert isinstance(MHC_II_EPITOPE_LENGTH, int)
        assert MHC_II_EPITOPE_LENGTH == 15

    def test_default_mhc_i_epitope_lengths(self) -> None:
        """DEFAULT_MHC_I_EPITOPE_LENGTHS contains [8, 9, 10, 11]."""
        assert isinstance(DEFAULT_MHC_I_EPITOPE_LENGTHS, list)
        assert DEFAULT_MHC_I_EPITOPE_LENGTHS == [8, 9, 10, 11]

    def test_max_retries_is_3(self) -> None:
        """MAX_RETRIES is 3 (documented default)."""
        assert MAX_RETRIES == 3

    def test_retry_base_delay_is_2(self) -> None:
        """RETRY_BASE_DELAY is 2.0 seconds (documented default)."""
        assert RETRY_BASE_DELAY == pytest.approx(2.0)


# ═══════════════════════════════════════════════════════════════════════════
# 2. Function existence and signatures
# ═══════════════════════════════════════════════════════════════════════════

class TestFunctionExistenceAndSignatures:
    """Verify all public symbols exist and are callable with documented params."""

    # --- Core prediction functions ---

    def test_predict_mhc_i_binding_is_callable(self) -> None:
        """Module-level predict_mhc_i_binding is callable."""
        assert callable(predict_mhc_i_binding)

    def test_predict_mhc_i_binding_signature(self) -> None:
        """predict_mhc_i_binding accepts (peptide_sequence, allele, epitope_length)."""
        sig = inspect.signature(predict_mhc_i_binding)
        params = list(sig.parameters.keys())
        assert "peptide_sequence" in params
        assert "allele" in params
        assert "epitope_length" in params

    def test_predict_mhc_ii_binding_is_callable(self) -> None:
        """Module-level predict_mhc_ii_binding is callable."""
        assert callable(predict_mhc_ii_binding)

    def test_predict_mhc_ii_binding_signature(self) -> None:
        """predict_mhc_ii_binding accepts (peptide_sequence, allele)."""
        sig = inspect.signature(predict_mhc_ii_binding)
        params = list(sig.parameters.keys())
        assert "peptide_sequence" in params
        assert "allele" in params

    def test_batch_predict_is_callable(self) -> None:
        """Module-level batch_predict is callable."""
        assert callable(batch_predict)

    def test_batch_predict_signature(self) -> None:
        """batch_predict accepts (protein_sequence, alleles, epitope_lengths)."""
        sig = inspect.signature(batch_predict)
        params = list(sig.parameters.keys())
        assert "protein_sequence" in params
        assert "alleles" in params
        assert "epitope_lengths" in params

    # --- Utility functions ---

    def test_classify_binding_rank_is_callable(self) -> None:
        """classify_binding_rank is callable."""
        assert callable(classify_binding_rank)

    def test_classify_binding_rank_signature(self) -> None:
        """classify_binding_rank accepts a single float parameter."""
        sig = inspect.signature(classify_binding_rank)
        assert list(sig.parameters.keys()) == ["rank"]

    def test_parse_netmhcpan_output_is_callable(self) -> None:
        """parse_netmhcpan_output is callable."""
        assert callable(parse_netmhcpan_output)

    def test_parse_netmhcpan_output_signature(self) -> None:
        """parse_netmhcpan_output accepts (output_text, allele)."""
        sig = inspect.signature(parse_netmhcpan_output)
        params = list(sig.parameters.keys())
        assert "output_text" in params
        assert "allele" in params

    def test_is_netmhcpan_available_is_callable(self) -> None:
        """is_netmhcpan_available is callable."""
        assert callable(is_netmhcpan_available)

    def test_is_netmhcpan_available_signature(self) -> None:
        """is_netmhcpan_available accepts (timeout)."""
        sig = inspect.signature(is_netmhcpan_available)
        assert "timeout" in sig.parameters

    def test_clear_cache_is_callable(self) -> None:
        """clear_cache is callable."""
        assert callable(clear_cache)

    def test_clear_cache_no_params(self) -> None:
        """clear_cache takes no parameters."""
        sig = inspect.signature(clear_cache)
        assert len(sig.parameters) == 0

    # --- Client class ---

    def test_netmhcpan_client_is_class(self) -> None:
        """NetMHCpanClient is a class."""
        assert isinstance(NetMHCpanClient, type)

    def test_client_predict_mhc_i_signature(self) -> None:
        """NetMHCpanClient.predict_mhc_i_binding has expected params."""
        sig = inspect.signature(NetMHCpanClient.predict_mhc_i_binding)
        params = list(sig.parameters.keys())
        assert "peptide_sequence" in params
        assert "allele" in params
        assert "epitope_length" in params

    def test_client_predict_mhc_ii_signature(self) -> None:
        """NetMHCpanClient.predict_mhc_ii_binding has expected params."""
        sig = inspect.signature(NetMHCpanClient.predict_mhc_ii_binding)
        params = list(sig.parameters.keys())
        assert "peptide_sequence" in params
        assert "allele" in params

    def test_client_batch_predict_signature(self) -> None:
        """NetMHCpanClient.batch_predict has expected params."""
        sig = inspect.signature(NetMHCpanClient.batch_predict)
        params = list(sig.parameters.keys())
        assert "protein_sequence" in params
        assert "alleles" in params
        assert "epitope_lengths" in params


# ═══════════════════════════════════════════════════════════════════════════
# 3. Input validation
# ═══════════════════════════════════════════════════════════════════════════

class TestInputValidation:
    """Test input validation in NetMHCpanClient — empty peptides, invalid proteins."""

    def _make_client(self, **kwargs):
        """Create a client with default settings for testing."""
        clear_cache()
        defaults = dict(max_retries=1, timeout=5.0, use_cache=True)
        defaults.update(kwargs)
        return NetMHCpanClient(**defaults)

    def test_empty_peptide_raises_mhc_i(self) -> None:
        """predict_mhc_i_binding raises NetMHCpanError on empty peptide."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="must not be empty"):
            client.predict_mhc_i_binding("", "HLA-A*02:01")

    def test_whitespace_peptide_raises_mhc_i(self) -> None:
        """predict_mhc_i_binding raises NetMHCpanError on whitespace-only peptide."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="must not be empty"):
            client.predict_mhc_i_binding("   ", "HLA-A*02:01")

    def test_empty_peptide_raises_mhc_ii(self) -> None:
        """predict_mhc_ii_binding raises NetMHCpanError on empty peptide."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="must not be empty"):
            client.predict_mhc_ii_binding("", "HLA-DRB1*01:01")

    def test_whitespace_peptide_raises_mhc_ii(self) -> None:
        """predict_mhc_ii_binding raises NetMHCpanError on whitespace-only peptide."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="must not be empty"):
            client.predict_mhc_ii_binding("   ", "HLA-DRB1*01:01")

    def test_batch_predict_invalid_protein_raises(self) -> None:
        """batch_predict raises NetMHCpanError on invalid protein (non-standard AA)."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="invalid"):
            client.batch_predict("MXG", ["HLA-A*02:01"])

    def test_batch_predict_empty_alleles_returns_empty(self) -> None:
        """batch_predict with empty alleles list returns empty results."""
        client = self._make_client()
        results = client.batch_predict("MAGPKWVTFIS", [])
        assert results == []

    def test_peptide_uppercased_before_validation(self) -> None:
        """Lowercase peptide is uppercased before empty-check (so '  ' still fails)."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError):
            client.predict_mhc_i_binding("  ", "HLA-A*02:01")

    def test_validation_error_includes_allele(self) -> None:
        """NetMHCpanError from empty peptide includes allele information."""
        client = self._make_client()
        with pytest.raises(NetMHCpanError) as exc_info:
            client.predict_mhc_i_binding("", "HLA-A*02:01")
        assert exc_info.value.allele == "HLA-A*02:01"


# ═══════════════════════════════════════════════════════════════════════════
# 4. Fallback behaviour
# ═══════════════════════════════════════════════════════════════════════════

class TestFallbackBehaviour:
    """Test fallback when NetMHCpan API returns no data or fails."""

    def _make_client(self, **kwargs):
        """Create a client with default settings for testing."""
        clear_cache()
        defaults = dict(max_retries=1, timeout=5.0, use_cache=True)
        defaults.update(kwargs)
        return NetMHCpanClient(**defaults)

    @patch("urllib.request.urlopen")
    def test_mhc_i_no_parsed_data_returns_non_binder(self, mock_urlopen) -> None:
        """predict_mhc_i_binding returns non_binder when API returns header-only output."""
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
        assert result.method == "netmhcpan"

    @patch("urllib.request.urlopen")
    def test_mhc_ii_no_parsed_data_returns_non_binder(self, mock_urlopen) -> None:
        """predict_mhc_ii_binding returns non_binder when API returns header-only output."""
        # Need output that _looks_like_netmhcpan_output returns True for,
        # but parse_netmhcpan_output returns empty list.
        header_only_output = """\
# NetMHCIIpan version 4.0
# Input is in FASTA format
----------------------------------------------------------------------------------------------------
  Pos  HLA  Peptide  Identity  Score  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
"""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = header_only_output.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = self._make_client()
        result = client.predict_mhc_ii_binding("AAAAAAAAAAAAAAA", "HLA-DRB1*01:01")
        assert result.binding_class == "non_binder"
        assert result.binding_score == 0.0

    @patch("urllib.request.urlopen")
    def test_mhc_i_network_error_raises(self, mock_urlopen) -> None:
        """predict_mhc_i_binding raises NetMHCpanError on network failure."""
        import urllib.error
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="API call failed"):
            client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")

    @patch("urllib.request.urlopen")
    def test_batch_predict_api_failure_falls_back_to_non_binder(self, mock_urlopen) -> None:
        """batch_predict appends non_binder with method='netmhcpan_failed' on API error."""
        import urllib.error
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        client = self._make_client()
        results = client.batch_predict("AAAAAAAAA", ["HLA-A*02:01"], epitope_lengths=[9])
        # Protein length 9, 1 allele, 1 epitope length → 1 result
        assert len(results) == 1
        assert results[0].binding_class == "non_binder"
        assert results[0].method == "netmhcpan_failed"
        assert results[0].binding_score == 0.0

    @patch("urllib.request.urlopen")
    def test_mhc_i_http_400_raises(self, mock_urlopen) -> None:
        """predict_mhc_i_binding raises NetMHCpanError on HTTP 400 (no retry)."""
        import urllib.error
        mock_urlopen.side_effect = urllib.error.HTTPError(
            url="", code=400, msg="Bad Request",
            hdrs={}, fp=MagicMock(read=MagicMock(return_value=b"Invalid allele")),
        )
        client = self._make_client()
        with pytest.raises(NetMHCpanError, match="API call failed"):
            client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")


# ═══════════════════════════════════════════════════════════════════════════
# 5. MHCBindingResult dataclass
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCBindingResult:
    """Test MHCBindingResult construction, defaults, and field types."""

    def test_create_result_with_all_fields(self) -> None:
        """MHCBindingResult can be created with all fields."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="SIINFEKLW",
            start_position=0,
            end_position=8,
            binding_score=0.85,
            ic50_nm=42.0,
            binding_class="strong_binder",
            rank=0.12,
            method="netmhcpan",
        )
        assert result.allele == "HLA-A*02:01"
        assert result.peptide == "SIINFEKLW"
        assert result.binding_score == pytest.approx(0.85)
        assert result.ic50_nm == pytest.approx(42.0)
        assert result.binding_class == "strong_binder"
        assert result.rank == pytest.approx(0.12)
        assert result.method == "netmhcpan"

    def test_result_defaults(self) -> None:
        """MHCBindingResult has sensible defaults for optional fields."""
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

    def test_binding_class_values(self) -> None:
        """Binding class can be set to all three valid values."""
        for bc in ("strong_binder", "weak_binder", "non_binder"):
            result = MHCBindingResult(
                allele="X", peptide="Y", start_position=0, end_position=0,
                binding_class=bc,
            )
            assert result.binding_class == bc

    def test_method_field_values(self) -> None:
        """Method field can be set to netmhcpan, pssm_fallback, or netmhcpan_failed."""
        for method in ("netmhcpan", "pssm_fallback", "netmhcpan_failed"):
            result = MHCBindingResult(
                allele="X", peptide="Y", start_position=0, end_position=0,
                method=method,
            )
            assert result.method == method

    def test_anchor_dicts_are_independent(self) -> None:
        """Each MHCBindingResult gets independent anchor dicts (no shared mutable default)."""
        r1 = MHCBindingResult(allele="A", peptide="P", start_position=0, end_position=0)
        r2 = MHCBindingResult(allele="A", peptide="P", start_position=0, end_position=0)
        r1.anchor_residues[0] = "K"
        assert 0 not in r2.anchor_residues

    def test_ic50_nm_can_be_none(self) -> None:
        """ic50_nm is None by default (no prediction data)."""
        result = MHCBindingResult(allele="X", peptide="Y", start_position=0, end_position=0)
        assert result.ic50_nm is None

    def test_rank_can_be_none(self) -> None:
        """rank is None by default (no rank data)."""
        result = MHCBindingResult(allele="X", peptide="Y", start_position=0, end_position=0)
        assert result.rank is None


# ═══════════════════════════════════════════════════════════════════════════
# 6. NetMHCpanCache
# ═══════════════════════════════════════════════════════════════════════════

class TestNetMHCpanCache:
    """Test NetMHCpanCache — put/get/eviction/hit_rate/clear/put_batch/key isolation."""

    def _make_result(self, allele: str = "HLA-A*02:01", peptide: str = "SIINFEKL",
                     score: float = 0.9) -> MHCBindingResult:
        """Helper to create a minimal MHCBindingResult for cache tests."""
        return MHCBindingResult(
            allele=allele, peptide=peptide,
            start_position=0, end_position=8,
            binding_score=score,
        )

    def test_create_cache(self) -> None:
        """Cache can be created with default settings."""
        cache = NetMHCpanCache()
        assert cache.size == 0
        assert cache.hits == 0
        assert cache.misses == 0

    def test_put_and_get(self) -> None:
        """Put a result and retrieve it."""
        cache = NetMHCpanCache()
        result = self._make_result()
        cache.put("HLA-A*02:01", "SIINFEKL", result)
        retrieved = cache.get("HLA-A*02:01", "SIINFEKL")
        assert retrieved is not None
        assert retrieved.allele == "HLA-A*02:01"
        assert retrieved.peptide == "SIINFEKL"
        assert retrieved.binding_score == pytest.approx(0.9)

    def test_get_miss_returns_none(self) -> None:
        """Getting a non-existent key returns None and increments misses."""
        cache = NetMHCpanCache()
        assert cache.get("HLA-A*02:01", "SIINFEKL") is None
        assert cache.misses == 1

    def test_hit_rate_tracking(self) -> None:
        """hit_rate reflects hits / (hits + misses)."""
        cache = NetMHCpanCache()
        cache.put("HLA-A*02:01", "SIINFEKL", self._make_result())
        cache.get("HLA-A*02:01", "SIINFEKL")  # hit
        cache.get("HLA-A*02:01", "MISS")       # miss
        assert cache.hits == 1
        assert cache.misses == 1
        assert cache.hit_rate == pytest.approx(0.5, abs=0.01)

    def test_hit_rate_no_accesses(self) -> None:
        """hit_rate is 0.0 when there have been no accesses."""
        cache = NetMHCpanCache()
        assert cache.hit_rate == 0.0

    def test_eviction_when_full(self) -> None:
        """FIFO eviction when max_size is exceeded."""
        cache = NetMHCpanCache(max_size=3)
        for i in range(5):
            cache.put("HLA-A*02:01", f"PEPTIDE{i}", self._make_result(peptide=f"PEPTIDE{i}"))
        # Only last 3 should remain
        assert cache.get("HLA-A*02:01", "PEPTIDE0") is None
        assert cache.get("HLA-A*02:01", "PEPTIDE1") is None
        assert cache.get("HLA-A*02:01", "PEPTIDE2") is not None
        assert cache.get("HLA-A*02:01", "PEPTIDE3") is not None
        assert cache.get("HLA-A*02:01", "PEPTIDE4") is not None

    def test_clear_resets_cache_and_stats(self) -> None:
        """clear() resets size, hits, and misses."""
        cache = NetMHCpanCache()
        cache.put("HLA-A*02:01", "SIINFEKL", self._make_result())
        cache.get("HLA-A*02:01", "SIINFEKL")
        cache.clear()
        assert cache.size == 0
        assert cache.hits == 0
        assert cache.misses == 0

    def test_cache_key_includes_epitope_length(self) -> None:
        """Cache keys differentiate by epitope length."""
        cache = NetMHCpanCache()
        r9 = MHCBindingResult(allele="HLA-A*02:01", peptide="SIINFEKLA",
                              start_position=0, end_position=8, binding_score=0.5)
        r10 = MHCBindingResult(allele="HLA-A*02:01", peptide="SIINFEKLA",
                               start_position=0, end_position=9, binding_score=0.8)
        cache.put("HLA-A*02:01", "SIINFEKLA", r9, epitope_length=9)
        cache.put("HLA-A*02:01", "SIINFEKLA", r10, epitope_length=10)
        got9 = cache.get("HLA-A*02:01", "SIINFEKLA", epitope_length=9)
        got10 = cache.get("HLA-A*02:01", "SIINFEKLA", epitope_length=10)
        assert got9 is not None
        assert got10 is not None
        assert got9.binding_score != got10.binding_score

    def test_put_batch(self) -> None:
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

    def test_cache_key_deterministic(self) -> None:
        """Same allele + peptide + length always produces the same cache key."""
        key1 = NetMHCpanCache._key("HLA-A*02:01", "SIINFEKL", 9)
        key2 = NetMHCpanCache._key("HLA-A*02:01", "SIINFEKL", 9)
        assert key1 == key2

    def test_cache_key_length(self) -> None:
        """Cache key is 16 hex characters (from SHA-256 digest)."""
        key = NetMHCpanCache._key("HLA-A*02:01", "SIINFEKL", 9)
        assert len(key) == 16
        assert all(c in "0123456789abcdef" for c in key)

    def test_different_peptides_different_keys(self) -> None:
        """Different peptides produce different cache keys."""
        key1 = NetMHCpanCache._key("HLA-A*02:01", "SIINFEKL", 9)
        key2 = NetMHCpanCache._key("HLA-A*02:01", "AAAAAAAAA", 9)
        assert key1 != key2

    def test_different_alleles_different_keys(self) -> None:
        """Different alleles produce different cache keys."""
        key1 = NetMHCpanCache._key("HLA-A*02:01", "SIINFEKL", 9)
        key2 = NetMHCpanCache._key("HLA-B*07:02", "SIINFEKL", 9)
        assert key1 != key2


# ═══════════════════════════════════════════════════════════════════════════
# 7. classify_binding_rank
# ═══════════════════════════════════════════════════════════════════════════

class TestClassifyBindingRank:
    """Test classify_binding_rank boundary values for all three bands."""

    def test_strong_binder(self) -> None:
        """Rank < 0.5% is classified as strong_binder."""
        assert classify_binding_rank(0.01) == "strong_binder"
        assert classify_binding_rank(0.49) == "strong_binder"
        assert classify_binding_rank(0.1) == "strong_binder"

    def test_weak_binder(self) -> None:
        """0.5% <= Rank < 2% is classified as weak_binder."""
        assert classify_binding_rank(0.5) == "weak_binder"
        assert classify_binding_rank(1.0) == "weak_binder"
        assert classify_binding_rank(1.99) == "weak_binder"

    def test_non_binder(self) -> None:
        """Rank >= 2% is classified as non_binder."""
        assert classify_binding_rank(2.0) == "non_binder"
        assert classify_binding_rank(5.0) == "non_binder"
        assert classify_binding_rank(50.0) == "non_binder"
        assert classify_binding_rank(100.0) == "non_binder"

    def test_returns_string(self) -> None:
        """classify_binding_rank always returns a string."""
        for rank in [-1.0, 0.0, 0.01, 0.5, 1.0, 2.0, 5.0, 100.0]:
            assert isinstance(classify_binding_rank(rank), str)

    def test_boundary_at_strong_threshold(self) -> None:
        """Exactly 0.5 is weak_binder (not strictly less than 0.5)."""
        assert classify_binding_rank(0.5) == "weak_binder"

    def test_boundary_at_weak_threshold(self) -> None:
        """Exactly 2.0 is non_binder (not strictly less than 2.0)."""
        assert classify_binding_rank(2.0) == "non_binder"


# ═══════════════════════════════════════════════════════════════════════════
# 8. _rank_to_binding_score
# ═══════════════════════════════════════════════════════════════════════════

class TestRankToBindingScore:
    """Test _rank_to_binding_score range, monotonicity, and edge cases."""

    def test_zero_rank_gives_max_score(self) -> None:
        """Rank 0 → score 1.0 (strongest possible)."""
        assert _rank_to_binding_score(0) == 1.0

    def test_negative_rank_gives_max_score(self) -> None:
        """Negative rank → score 1.0."""
        assert _rank_to_binding_score(-1.0) == 1.0

    def test_very_low_rank_high_score(self) -> None:
        """Very low rank → high binding score."""
        score = _rank_to_binding_score(0.01)
        assert score > 0.9

    def test_score_range(self) -> None:
        """Binding score is always in [0, 1]."""
        for rank in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]:
            score = _rank_to_binding_score(rank)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for rank {rank}"

    def test_monotonicity(self) -> None:
        """Higher rank → lower binding score (monotonically decreasing)."""
        ranks = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
        scores = [_rank_to_binding_score(r) for r in ranks]
        for i in range(len(scores) - 1):
            assert scores[i] >= scores[i + 1], (
                f"Non-monotonic: rank={ranks[i]} → score={scores[i]}, "
                f"rank={ranks[i+1]} → score={scores[i+1]}"
            )

    def test_high_rank_gives_low_score(self) -> None:
        """Rank 100 → score ≤ 0.5 (log10(100)/4 = 0.5, so score = 1 - 0.5)."""
        score = _rank_to_binding_score(100.0)
        assert score <= 0.5
        assert score >= 0.0

    def test_result_is_rounded(self) -> None:
        """Score is rounded to 6 decimal places."""
        score = _rank_to_binding_score(1.5)
        # Should have at most 6 decimal digits
        assert score == round(score, 6)


# ═══════════════════════════════════════════════════════════════════════════
# 9. parse_netmhcpan_output
# ═══════════════════════════════════════════════════════════════════════════

class TestParseNetMHCpanOutput:
    """Test parse_netmhcpan_output and _parse_data_line — parsing, header skipping, positions."""

    def test_parse_mhc_i_output(self) -> None:
        """Parse MHC-I output with three binders of different strengths."""
        results = parse_netmhcpan_output(SAMPLE_MHC_I_OUTPUT, "HLA-A*02:01")
        assert len(results) == 3

        r0 = results[0]
        assert r0["peptide"] == "SIINFEKLW"
        assert r0["binding_class"] == "strong_binder"
        assert r0["rank"] == pytest.approx(0.12, abs=0.01)

        r1 = results[1]
        assert r1["peptide"] == "IINFEKLWA"
        assert r1["binding_class"] == "weak_binder"
        assert r1["rank"] == pytest.approx(1.50, abs=0.01)

        r2 = results[2]
        assert r2["peptide"] == "INFEKLWAA"
        assert r2["binding_class"] == "non_binder"
        assert r2["rank"] == pytest.approx(8.00, abs=0.01)

    def test_parse_mhc_ii_output(self) -> None:
        """Parse MHC-II output."""
        results = parse_netmhcpan_output(SAMPLE_MHC_II_OUTPUT, "HLA-DRB1*01:01")
        assert len(results) == 2
        r0 = results[0]
        assert r0["binding_class"] == "strong_binder"

    def test_parse_empty_output(self) -> None:
        """Empty output returns no results."""
        results = parse_netmhcpan_output("", "HLA-A*02:01")
        assert results == []

    def test_parse_header_only(self) -> None:
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

    def test_position_is_zero_based(self) -> None:
        """Positions in parsed output are 0-based (NetMHCpan output is 1-based)."""
        results = parse_netmhcpan_output(SAMPLE_MHC_I_OUTPUT, "HLA-A*02:01")
        if results:
            # First position in NetMHCpan output is 1 → should be 0 in parsed
            assert results[0]["position"] == 0

    def test_ic50_positive(self) -> None:
        """IC50 values are positive numbers."""
        results = parse_netmhcpan_output(SAMPLE_MHC_I_OUTPUT, "HLA-A*02:01")
        for r in results:
            assert r["ic50_nm"] > 0

    def test_parse_skips_comment_lines(self) -> None:
        """Lines starting with # are skipped."""
        text = "# This is a comment\n1  HLA-A*02:01  SIINFEKL  SIINFEKL  0 0 0 0 0 0  Protein1  0.500  SB  0.30\n"
        results = parse_netmhcpan_output(text, "HLA-A*02:01")
        assert len(results) == 1

    def test_parse_skips_separator_lines(self) -> None:
        """Lines starting with - are skipped."""
        text = "----separator----\n1  HLA-A*02:01  SIINFEKL  SIINFEKL  0 0 0 0 0 0  Protein1  0.500  SB  0.30\n"
        results = parse_netmhcpan_output(text, "HLA-A*02:01")
        assert len(results) == 1

    def test_parse_data_line_returns_none_for_short_line(self) -> None:
        """_parse_data_line returns None for lines with too few tokens."""
        assert _parse_data_line("1  HLA") is None

    def test_parse_data_line_returns_none_for_no_peptide(self) -> None:
        """_parse_data_line returns None if no peptide token is found."""
        assert _parse_data_line("1  HLA-A*02:01  12345  Protein1  0.500  0.30") is None

    def test_parsed_dict_keys(self) -> None:
        """Each parsed result dict has all expected keys."""
        results = parse_netmhcpan_output(SAMPLE_MHC_I_OUTPUT, "HLA-A*02:01")
        expected_keys = {"position", "allele", "peptide", "score", "ic50_nm", "rank", "binding_class"}
        for r in results:
            assert expected_keys.issubset(r.keys())


# ═══════════════════════════════════════════════════════════════════════════
# 10. Helper functions
# ═══════════════════════════════════════════════════════════════════════════

class TestIsMHCIIAllele:
    """Test _is_mhc_ii_allele for correct classification."""

    def test_mhc_ii_dr_allele(self) -> None:
        """HLA-DR alleles are MHC-II."""
        assert _is_mhc_ii_allele("HLA-DRB1*01:01") is True

    def test_mhc_ii_dq_allele(self) -> None:
        """HLA-DQ alleles are MHC-II."""
        assert _is_mhc_ii_allele("HLA-DQB1*03:01") is True

    def test_mhc_ii_dp_allele(self) -> None:
        """HLA-DP alleles are MHC-II."""
        assert _is_mhc_ii_allele("HLA-DPB1*04:01") is True

    def test_mhc_ii_mouse_allele(self) -> None:
        """Mouse H2-I alleles are MHC-II."""
        assert _is_mhc_ii_allele("H2-IAb") is True

    def test_mhc_i_allele_not_ii(self) -> None:
        """HLA-A alleles are NOT MHC-II."""
        assert _is_mhc_ii_allele("HLA-A*02:01") is False

    def test_mhc_i_b_allele_not_ii(self) -> None:
        """HLA-B alleles are NOT MHC-II."""
        assert _is_mhc_ii_allele("HLA-B*07:02") is False

    def test_mhc_i_mouse_allele_not_ii(self) -> None:
        """Mouse H2-Kb alleles are NOT MHC-II."""
        assert _is_mhc_ii_allele("H2-Kb") is False

    def test_case_insensitive(self) -> None:
        """Allele classification is case-insensitive."""
        assert _is_mhc_ii_allele("hla-drb1*01:01") is True
        assert _is_mhc_ii_allele("HLA-DRB1*01:01") is True


class TestLooksLikeNetMHCpanOutput:
    """Test _looks_like_netmhcpan_output detection."""

    def test_recognizes_mhc_i_output(self) -> None:
        """Recognizes MHC-I output."""
        assert _looks_like_netmhcpan_output(SAMPLE_MHC_I_OUTPUT) is True

    def test_recognizes_mhc_ii_output(self) -> None:
        """Recognizes MHC-II output."""
        assert _looks_like_netmhcpan_output(SAMPLE_MHC_II_OUTPUT) is True

    def test_rejects_empty_string(self) -> None:
        """Rejects empty string."""
        assert _looks_like_netmhcpan_output("") is False

    def test_rejects_random_html(self) -> None:
        """Rejects generic HTML."""
        assert _looks_like_netmhcpan_output("<html><body>Hello</body></html>") is False

    def test_requires_at_least_two_markers(self) -> None:
        """Requires at least 2 characteristic markers to match."""
        # Only one marker → not enough
        assert _looks_like_netmhcpan_output("# Rank") is False
        # Two markers → enough
        assert _looks_like_netmhcpan_output("# NetMHCpan\n%Rank") is True


class TestExtractResultUrl:
    """Test _extract_result_url for job URL extraction."""

    def test_extract_from_meta_refresh(self) -> None:
        """Extracts URL from meta-refresh tag."""
        url = _extract_result_url(SAMPLE_JOB_SUBMIT_HTML)
        assert url is not None
        assert "jobid=abc123" in url

    def test_extract_from_jobid(self) -> None:
        """Extracts job ID and constructs URL."""
        html = '<p>Your job ID is: jobid=xyz789</p>'
        url = _extract_result_url(html)
        assert url is not None
        assert "xyz789" in url

    def test_no_url_in_html(self) -> None:
        """Returns None when no result URL is found."""
        html = "<html><body>No redirect here</body></html>"
        url = _extract_result_url(html)
        assert url is None


class TestExtractErrorMessage:
    """Test _extract_error_message for error detection."""

    def test_extract_allele_error(self) -> None:
        """Extracts 'Unrecognized allele' error."""
        msg = _extract_error_message(SAMPLE_ERROR_HTML)
        assert msg is not None
        assert "Unrecognized allele" in msg

    def test_extract_generic_error(self) -> None:
        """Extracts a generic error message."""
        html = "<p>Error: Invalid sequence format</p>"
        msg = _extract_error_message(html)
        assert msg is not None
        assert "Invalid sequence format" in msg

    def test_no_error(self) -> None:
        """Returns None when no error is found."""
        html = "<html><body>Success</body></html>"
        msg = _extract_error_message(html)
        assert msg is None

    def test_unrecognized_allele_lower(self) -> None:
        """Detects lowercase 'unrecognized allele' pattern."""
        html = "<p>unrecognized allele: XYZ</p>"
        msg = _extract_error_message(html)
        assert msg is not None

    def test_invalid_sequence_pattern(self) -> None:
        """Detects 'invalid sequence' pattern."""
        html = "<p>invalid sequence provided</p>"
        msg = _extract_error_message(html)
        assert msg is not None


# ═══════════════════════════════════════════════════════════════════════════
# 11. NetMHCpanError
# ═══════════════════════════════════════════════════════════════════════════

class TestNetMHCpanError:
    """Test NetMHCpanError exception class hierarchy and attributes."""

    def test_error_creation(self) -> None:
        """NetMHCpanError can be created with reason and allele."""
        error = NetMHCpanError("API failed", allele="HLA-A*02:01")
        assert "API failed" in str(error)
        assert "HLA-A*02:01" in str(error)

    def test_error_is_immunogenicity_error(self) -> None:
        """NetMHCpanError is a subclass of ImmunogenicityError."""
        assert issubclass(NetMHCpanError, ImmunogenicityError)

    def test_error_is_engine_error(self) -> None:
        """NetMHCpanError is a subclass of EngineError."""
        assert issubclass(NetMHCpanError, EngineError)

    def test_error_without_allele(self) -> None:
        """NetMHCpanError works without allele."""
        error = NetMHCpanError("Network timeout")
        assert "Network timeout" in str(error)
        assert error.allele is None

    def test_error_catchable_as_engine_error(self) -> None:
        """NetMHCpanError can be caught as EngineError."""
        with pytest.raises(EngineError):
            raise NetMHCpanError("test")

    def test_error_catchable_as_immunogenicity_error(self) -> None:
        """NetMHCpanError can be caught as ImmunogenicityError."""
        with pytest.raises(ImmunogenicityError):
            raise NetMHCpanError("test")

    def test_error_reason_attribute(self) -> None:
        """NetMHCpanError stores reason attribute."""
        error = NetMHCpanError("some reason")
        # Note: ImmunogenicityError.__init__ may overwrite self.reason with the
        # formatted message, so we check the reason is captured in the string.
        assert "some reason" in str(error)

    def test_error_allele_attribute(self) -> None:
        """NetMHCpanError stores allele attribute."""
        error = NetMHCpanError("test", allele="HLA-B*07:02")
        assert error.allele == "HLA-B*07:02"


# ═══════════════════════════════════════════════════════════════════════════
# 12. Module-level helpers
# ═══════════════════════════════════════════════════════════════════════════

class TestModuleHelpers:
    """Test module-level helper functions and default instances."""

    def test_clear_cache_no_error(self) -> None:
        """clear_cache() runs without error."""
        clear_cache()  # should not raise

    def test_get_default_cache_returns_cache(self) -> None:
        """_get_default_cache returns a NetMHCpanCache instance."""
        cache = _get_default_cache()
        assert isinstance(cache, NetMHCpanCache)

    def test_get_default_client_returns_client(self) -> None:
        """_get_default_client returns a NetMHCpanClient instance."""
        client = _get_default_client()
        assert isinstance(client, NetMHCpanClient)

    def test_clear_cache_resets_default(self) -> None:
        """clear_cache() resets the module-level default cache."""
        cache = _get_default_cache()
        cache.put("HLA-A*02:01", "TEST", MHCBindingResult(
            allele="HLA-A*02:01", peptide="TEST",
            start_position=0, end_position=3,
        ))
        assert cache.size >= 1
        clear_cache()
        assert cache.size == 0

    @patch("urllib.request.urlopen")
    def test_is_netmhcpan_available_returns_bool(self, mock_urlopen) -> None:
        """is_netmhcpan_available returns a boolean."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp
        result = is_netmhcpan_available(timeout=1.0)
        assert isinstance(result, bool)

    @patch("urllib.request.urlopen")
    def test_is_netmhcpan_available_true_when_reachable(self, mock_urlopen) -> None:
        """is_netmhcpan_available returns True when server responds."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp
        assert is_netmhcpan_available(timeout=1.0) is True

    @patch("urllib.request.urlopen")
    def test_is_netmhcpan_available_false_when_down(self, mock_urlopen) -> None:
        """is_netmhcpan_available returns False when server is unreachable."""
        import urllib.error
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")
        assert is_netmhcpan_available(timeout=1.0) is False


# ═══════════════════════════════════════════════════════════════════════════
# 13. NetMHCpanClient construction
# ═══════════════════════════════════════════════════════════════════════════

class TestNetMHCpanClientConstruction:
    """Test NetMHCpanClient instantiation and parameter handling."""

    def test_default_construction(self) -> None:
        """Client can be created with all defaults."""
        client = NetMHCpanClient()
        assert client.api_url is not None
        assert client.timeout > 0
        assert client.max_retries >= 1
        assert client.use_cache is True

    def test_custom_construction(self) -> None:
        """Client accepts custom parameters."""
        cache = NetMHCpanCache()
        client = NetMHCpanClient(
            api_url="https://example.com",
            timeout=30.0,
            max_retries=5,
            cache=cache,
            use_cache=False,
        )
        assert client.api_url == "https://example.com"
        assert client.timeout == 30.0
        assert client.max_retries == 5
        assert client._cache is cache
        assert client.use_cache is False

    def test_client_uses_default_cache_when_none(self) -> None:
        """Client uses the module-level default cache when cache=None."""
        client = NetMHCpanClient(cache=None)
        assert client._cache is _get_default_cache()

    def test_client_uses_custom_cache(self) -> None:
        """Client uses a provided custom cache."""
        custom_cache = NetMHCpanCache()
        client = NetMHCpanClient(cache=custom_cache)
        assert client._cache is custom_cache

    @patch("urllib.request.urlopen")
    def test_client_cache_hit_skips_api(self, mock_urlopen) -> None:
        """Client returns cached result on second call without hitting API."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_MHC_I_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        clear_cache()
        client = NetMHCpanClient(max_retries=1, timeout=5.0)
        result1 = client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")
        result2 = client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")
        # API should only be called once (second hit is from cache)
        assert mock_urlopen.call_count == 1
        assert result1.allele == result2.allele

    @patch("urllib.request.urlopen")
    def test_client_no_cache_calls_api_each_time(self, mock_urlopen) -> None:
        """Client with use_cache=False calls API on every invocation."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.read.return_value = SAMPLE_MHC_I_OUTPUT.encode("utf-8")
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        client = NetMHCpanClient(max_retries=1, timeout=5.0, use_cache=False)
        client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")
        client.predict_mhc_i_binding("SIINFEKLW", "HLA-A*02:01")
        # API should be called twice (no caching)
        assert mock_urlopen.call_count == 2
