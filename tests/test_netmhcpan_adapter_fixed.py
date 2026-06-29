"""Tests for the fixed NetMHCpan adapter — fallback chain integration.

These tests verify that the adapter functions added in F4.7 work correctly
as part of the MHC prediction fallback chain:

  1. MHCflurry (if installed)
  2. NetMHCpan (if installed)
  3. Precomputed database
  4. PSSM fallback

Key properties under test:
  - is_netmhcpan_installed() returns a bool
  - is_netmhcpan_available() returns a bool
  - predict_binding_netmhcpan() returns None when NetMHCpan is not available
  - batch_predict_binding_netmhcpan() returns None when not available
  - The adapter never crashes, regardless of environment
"""
from __future__ import annotations

import inspect
from unittest.mock import MagicMock, patch

import pytest

from biocompiler.immunogenicity.netmhcpan import (
    MHCBindingResult,
    NetMHCpanError,
    _predict_binding_local,
    batch_predict_binding_netmhcpan,
    clear_cache,
    is_netmhcpan_available,
    is_netmhcpan_installed,
    predict_binding_netmhcpan,
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. is_netmhcpan_installed
# ═══════════════════════════════════════════════════════════════════════════

class TestIsNetMHCpanInstalled:
    """Test is_netmhcpan_installed() — local binary detection."""

    def setup_method(self):
        """Reset the installed cache before each test."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

    def test_returns_bool(self):
        """is_netmhcpan_installed() always returns a bool."""
        result = is_netmhcpan_installed()
        assert isinstance(result, bool)

    @patch("shutil.which", return_value=None)
    def test_returns_false_when_not_on_path(self, mock_which):
        """Returns False when no NetMHCpan binary is found on PATH."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None
        result = is_netmhcpan_installed()
        assert result is False

    @patch("shutil.which", return_value="/usr/local/bin/netMHCpan")
    def test_returns_true_when_on_path(self, mock_which):
        """Returns True when a NetMHCpan binary is found on PATH."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None
        result = is_netmhcpan_installed()
        assert result is True

    @patch("shutil.which")
    def test_checks_all_binary_names(self, mock_which):
        """Checks all names in binary list until one is found."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

        # First two names not found, third found
        mock_which.side_effect = [None, None, "/usr/bin/NetMHCpan"]
        result = is_netmhcpan_installed()
        assert result is True
        assert mock_which.call_count == 3

    @patch("shutil.which", return_value="/usr/bin/netMHCpan")
    def test_result_is_cached(self, mock_which):
        """Result is cached after first call (no repeated PATH lookups)."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

        is_netmhcpan_installed()  # First call — does lookup
        is_netmhcpan_installed()  # Second call — uses cache
        assert mock_which.call_count == 1  # Only called once


# ═══════════════════════════════════════════════════════════════════════════
# 2. is_netmhcpan_available
# ═══════════════════════════════════════════════════════════════════════════

class TestIsNetMHCpanAvailable:
    """Test is_netmhcpan_available() — combined local + web check."""

    def setup_method(self):
        """Reset caches before each test."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

    def test_returns_bool(self):
        """is_netmhcpan_available() always returns a bool."""
        result = is_netmhcpan_available(timeout=0.01)
        assert isinstance(result, bool)

    @patch("shutil.which", return_value="/usr/bin/netMHCpan")
    def test_available_when_installed(self, mock_which):
        """Returns True when local binary is installed (no web check needed)."""
        result = is_netmhcpan_available(timeout=0.01)
        assert result is True

    @patch("shutil.which", return_value=None)
    @patch("urllib.request.urlopen")
    def test_available_when_web_api_reachable(self, mock_urlopen, mock_which):
        """Returns True when web API is reachable but local binary absent."""
        mock_resp = MagicMock()
        mock_resp.status = 200
        mock_resp.__enter__ = MagicMock(return_value=mock_resp)
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        result = is_netmhcpan_available(timeout=0.01)
        assert result is True

    @patch("shutil.which", return_value=None)
    @patch("urllib.request.urlopen")
    def test_unavailable_when_neither_installed_nor_reachable(self, mock_urlopen, mock_which):
        """Returns False when local binary is absent and web API unreachable."""
        import urllib.error
        mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

        result = is_netmhcpan_available(timeout=0.01)
        assert result is False

    @patch("shutil.which", return_value=None)
    @patch("urllib.request.urlopen")
    def test_unavailable_on_timeout(self, mock_urlopen, mock_which):
        """Returns False on timeout."""
        mock_urlopen.side_effect = TimeoutError()
        result = is_netmhcpan_available(timeout=0.01)
        assert result is False


# ═══════════════════════════════════════════════════════════════════════════
# 3. predict_binding_netmhcpan
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictBindingNetMHCpan:
    """Test predict_binding_netmhcpan() — adapter for the fallback chain."""

    def setup_method(self):
        """Reset caches before each test."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_returns_none_when_not_available(self, mock_avail):
        """Returns None when NetMHCpan is not available."""
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_does_not_crash_when_unavailable(self, mock_avail):
        """The adapter never crashes, even with unusual inputs, when unavailable."""
        # Empty-ish but valid-looking args — should just return None
        result = predict_binding_netmhcpan("HLA-A*02:01", "")
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=False)
    @patch("biocompiler.netmhcpan.predict_mhc_i_binding")
    def test_delegates_to_web_api_when_not_installed(self, mock_predict, mock_installed, mock_avail):
        """Delegates to predict_mhc_i_binding (web API) when not installed locally."""
        mock_predict.return_value = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=8,
            binding_score=0.9, method="netmhcpan",
        )
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is not None
        assert result.allele == "HLA-A*02:01"
        mock_predict.assert_called_once()

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=True)
    @patch("biocompiler.netmhcpan._predict_binding_local")
    def test_delegates_to_local_when_installed(self, mock_local, mock_installed, mock_avail):
        """Delegates to _predict_binding_local when installed locally."""
        mock_local.return_value = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=8,
            binding_score=0.8, method="netmhcpan_local",
        )
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is not None
        assert result.method == "netmhcpan_local"
        mock_local.assert_called_once()

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=False)
    @patch("biocompiler.netmhcpan.predict_mhc_i_binding")
    def test_returns_none_on_prediction_error(self, mock_predict, mock_installed, mock_avail):
        """Returns None (not an exception) when prediction fails."""
        mock_predict.side_effect = NetMHCpanError("API call failed")
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=False)
    @patch("biocompiler.netmhcpan.predict_mhc_i_binding")
    def test_returns_none_on_unexpected_error(self, mock_predict, mock_installed, mock_avail):
        """Returns None even on unexpected errors (not just NetMHCpanError)."""
        mock_predict.side_effect = RuntimeError("Unexpected error")
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=True)
    @patch("biocompiler.netmhcpan._predict_binding_local")
    def test_returns_none_on_local_error(self, mock_local, mock_installed, mock_avail):
        """Returns None when local prediction fails."""
        mock_local.side_effect = NetMHCpanError("Local binary error")
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is None

    def test_function_signature(self):
        """predict_binding_netmhcpan has the expected signature."""
        sig = inspect.signature(predict_binding_netmhcpan)
        params = list(sig.parameters.keys())
        assert "allele" in params
        assert "peptide" in params
        assert "epitope_length" in params

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=False)
    @patch("biocompiler.netmhcpan.predict_mhc_i_binding")
    def test_result_is_mhc_binding_result(self, mock_predict, mock_installed, mock_avail):
        """Successful result is an MHCBindingResult instance."""
        mock_predict.return_value = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=8,
        )
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert isinstance(result, MHCBindingResult)


# ═══════════════════════════════════════════════════════════════════════════
# 4. batch_predict_binding_netmhcpan
# ═══════════════════════════════════════════════════════════════════════════

class TestBatchPredictBindingNetMHCpan:
    """Test batch_predict_binding_netmhcpan() — batch adapter for fallback chain."""

    def setup_method(self):
        """Reset caches before each test."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_returns_none_when_not_available(self, mock_avail):
        """Returns None when NetMHCpan is not available."""
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_does_not_crash_with_empty_alleles(self, mock_avail):
        """Does not crash even with empty allele list when unavailable."""
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", [])
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.batch_predict")
    def test_delegates_to_batch_predict(self, mock_batch, mock_avail):
        """Delegates to batch_predict when available."""
        mock_batch.return_value = [
            MHCBindingResult(
                allele="HLA-A*02:01", peptide="MAGPKWVTF",
                start_position=0, end_position=8,
            )
        ]
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
        assert result is not None
        assert len(result) == 1
        mock_batch.assert_called_once()

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.batch_predict")
    def test_returns_none_on_error(self, mock_batch, mock_avail):
        """Returns None when batch_predict raises an exception."""
        mock_batch.side_effect = NetMHCpanError("API call failed")
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.batch_predict")
    def test_returns_none_on_unexpected_error(self, mock_batch, mock_avail):
        """Returns None even on unexpected errors."""
        mock_batch.side_effect = RuntimeError("Unexpected")
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
        assert result is None

    def test_function_signature(self):
        """batch_predict_binding_netmhcpan has the expected signature."""
        sig = inspect.signature(batch_predict_binding_netmhcpan)
        params = list(sig.parameters.keys())
        assert "protein_sequence" in params
        assert "alleles" in params
        assert "epitope_lengths" in params


# ═══════════════════════════════════════════════════════════════════════════
# 5. _predict_binding_local
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictBindingLocal:
    """Test _predict_binding_local() — local binary execution."""

    def test_raises_on_empty_peptide(self):
        """Raises NetMHCpanError on empty peptide."""
        with pytest.raises(NetMHCpanError, match="must not be empty"):
            _predict_binding_local("HLA-A*02:01", "")

    @patch("shutil.which", return_value=None)
    def test_raises_when_binary_not_found(self, mock_which):
        """Raises NetMHCpanError when no binary is found on PATH."""
        with pytest.raises(NetMHCpanError, match="not found"):
            _predict_binding_local("HLA-A*02:01", "SIINFEKL")

    @patch("shutil.which", return_value="/usr/bin/netMHCpan")
    @patch("subprocess.run")
    def test_parses_successful_output(self, mock_run, mock_which):
        """Parses output from local NetMHCpan binary successfully."""
        sample_output = """\
# NetMHCpan version 4.1
----------------------------------------------------------------------------------------------------
  Pos  HLA         Peptide  Core Of Gp Gl Ip Il Ic  Identity  1-log50k(IC50)  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
    1  HLA-A*02:01  SIINFEKLW     SIINFEKLW  0  0  0  0  0  0   Protein1            0.025           <0.5%  0.12
----------------------------------------------------------------------------------------------------
"""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = sample_output
        mock_result.stderr = ""
        mock_run.return_value = mock_result

        result = _predict_binding_local("HLA-A*02:01", "SIINFEKLW", epitope_length=9)
        assert isinstance(result, MHCBindingResult)
        assert result.method == "netmhcpan_local"
        assert result.allele == "HLA-A*02:01"

    @patch("shutil.which", return_value="/usr/bin/netMHCpan")
    @patch("subprocess.run")
    def test_raises_on_nonzero_exit(self, mock_run, mock_which):
        """Raises NetMHCpanError when local binary exits with non-zero code."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = ""
        mock_result.stderr = "Error: invalid allele"
        mock_run.return_value = mock_result

        with pytest.raises(NetMHCpanError, match="exited with code 1"):
            _predict_binding_local("HLA-INVALID", "SIINFEKL")

    @patch("shutil.which", return_value="/usr/bin/netMHCpan")
    @patch("subprocess.run")
    def test_raises_on_timeout(self, mock_run, mock_which):
        """Raises NetMHCpanError when local binary times out."""
        import subprocess
        mock_run.side_effect = subprocess.TimeoutExpired(cmd="netMHCpan", timeout=120)

        with pytest.raises(NetMHCpanError, match="timed out"):
            _predict_binding_local("HLA-A*02:01", "SIINFEKL")

    @patch("shutil.which", return_value="/usr/bin/netMHCpan")
    @patch("subprocess.run")
    def test_returns_non_binder_on_no_parsed_data(self, mock_run, mock_which):
        """Returns non_binder result when output has no parseable data lines."""
        no_data_output = """\
# NetMHCpan version 4.1
----------------------------------------------------------------------------------------------------
  Pos  HLA  Peptide  Core Of Gp Gl Ip Il Ic  Identity  1-log50k(IC50)  BindLevel  %Rank
----------------------------------------------------------------------------------------------------
"""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = no_data_output
        mock_result.stderr = ""
        mock_run.return_value = mock_result

        result = _predict_binding_local("HLA-A*02:01", "SIINFEKLW", epitope_length=9)
        assert result.binding_class == "non_binder"
        assert result.method == "netmhcpan_local"


# ═══════════════════════════════════════════════════════════════════════════
# 6. Adapter robustness — never crash
# ═══════════════════════════════════════════════════════════════════════════

class TestAdapterRobustness:
    """Verify that the adapter never crashes in any scenario."""

    def setup_method(self):
        """Reset caches before each test."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_predict_adapter_no_crash_unavailable(self, mock_avail):
        """predict_binding_netmhcpan never crashes when unavailable."""
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_batch_adapter_no_crash_unavailable(self, mock_avail):
        """batch_predict_binding_netmhcpan never crashes when unavailable."""
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
        assert result is None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=False)
    @patch("biocompiler.netmhcpan.predict_mhc_i_binding")
    def test_predict_adapter_no_crash_on_any_exception(self, mock_predict, mock_installed, mock_avail):
        """predict_binding_netmhcpan catches any exception and returns None."""
        for exc in [
            NetMHCpanError("test"),
            ValueError("bad value"),
            RuntimeError("runtime error"),
            ConnectionError("connection lost"),
            OSError("os error"),
        ]:
            mock_predict.side_effect = exc
            result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
            assert result is None, f"Expected None for {type(exc).__name__}"

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.batch_predict")
    def test_batch_adapter_no_crash_on_any_exception(self, mock_batch, mock_avail):
        """batch_predict_binding_netmhcpan catches any exception and returns None."""
        for exc in [
            NetMHCpanError("test"),
            ValueError("bad value"),
            RuntimeError("runtime error"),
            ConnectionError("connection lost"),
        ]:
            mock_batch.side_effect = exc
            result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
            assert result is None, f"Expected None for {type(exc).__name__}"


# ═══════════════════════════════════════════════════════════════════════════
# 7. Fallback chain integration
# ═══════════════════════════════════════════════════════════════════════════

class TestFallbackChainIntegration:
    """Verify the adapter fits correctly into the fallback chain pattern."""

    def setup_method(self):
        """Reset caches before each test."""
        import biocompiler.immunogenicity.netmhcpan as mod
        mod._netmhcpan_installed_cache = None

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_fallback_chain_skips_netmhcpan_when_unavailable(self, mock_avail):
        """Simulated fallback chain correctly skips NetMHCpan when unavailable."""
        # Simulate the fallback chain pattern
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")

        # In a real fallback chain, the caller would do:
        # if result is None:
        #     result = next_predictor(allele, peptide)
        assert result is None, "Chain should skip to next predictor"

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=True)
    @patch("biocompiler.netmhcpan.is_netmhcpan_installed", return_value=False)
    @patch("biocompiler.netmhcpan.predict_mhc_i_binding")
    def test_fallback_chain_uses_netmhcpan_when_available(self, mock_predict, mock_installed, mock_avail):
        """Simulated fallback chain uses NetMHCpan when available."""
        mock_predict.return_value = MHCBindingResult(
            allele="HLA-A*02:01", peptide="SIINFEKL",
            start_position=0, end_position=8,
            binding_score=0.9, binding_class="strong_binder",
            rank=0.12, method="netmhcpan",
        )
        result = predict_binding_netmhcpan("HLA-A*02:01", "SIINFEKL")
        assert result is not None
        assert isinstance(result, MHCBindingResult)

    @patch("biocompiler.netmhcpan.is_netmhcpan_available", return_value=False)
    def test_batch_fallback_returns_none_for_chain(self, mock_avail):
        """Batch fallback chain correctly receives None."""
        result = batch_predict_binding_netmhcpan("MAGPKWVTFIS", ["HLA-A*02:01"])
        # In a real fallback chain:
        # if result is None:
        #     result = next_batch_predictor(protein, alleles)
        assert result is None


# ═══════════════════════════════════════════════════════════════════════════
# 8. Exports and public API
# ═══════════════════════════════════════════════════════════════════════════

class TestPublicAPI:
    """Verify new functions are exported in __all__."""

    def test_is_netmhcpan_installed_in_all(self):
        """is_netmhcpan_installed is exported in __all__."""
        from biocompiler.immunogenicity.netmhcpan import __all__
        assert "is_netmhcpan_installed" in __all__

    def test_is_netmhcpan_available_in_all(self):
        """is_netmhcpan_available is exported in __all__."""
        from biocompiler.immunogenicity.netmhcpan import __all__
        assert "is_netmhcpan_available" in __all__

    def test_predict_binding_netmhcpan_in_all(self):
        """predict_binding_netmhcpan is exported in __all__."""
        from biocompiler.immunogenicity.netmhcpan import __all__
        assert "predict_binding_netmhcpan" in __all__

    def test_batch_predict_binding_netmhcpan_in_all(self):
        """batch_predict_binding_netmhcpan is exported in __all__."""
        from biocompiler.immunogenicity.netmhcpan import __all__
        assert "batch_predict_binding_netmhcpan" in __all__

    def test_all_new_functions_importable(self):
        """All new public functions can be imported from the module."""
        from biocompiler.immunogenicity.netmhcpan import (
            is_netmhcpan_installed,
            is_netmhcpan_available,
            predict_binding_netmhcpan,
            batch_predict_binding_netmhcpan,
        )
        assert callable(is_netmhcpan_installed)
        assert callable(is_netmhcpan_available)
        assert callable(predict_binding_netmhcpan)
        assert callable(batch_predict_binding_netmhcpan)
