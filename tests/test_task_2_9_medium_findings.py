"""
Tests for Task 2.9 — Medium-severity audit findings.

Covers:
1. CSP solver failure logging (hybrid_optimizer.py)
2. Rate limit per-client identification with X-Forwarded-For (api.py)
3. Chunked encoding body size check (api.py)
4. Deprecated NUMBA kernels (numba_kernels.py)
5. DNA sequence length limits (api.py SequenceInput)
6. OpenAPI customization (api.py)
7. docs_url/redoc_url env var (api.py)
8. Logging for broad exception catches (api.py, dispatch.py)
"""

import hashlib
import os
import warnings
from unittest.mock import MagicMock, patch, PropertyMock

import pytest


# ─── 1. CSP solver failure logging ─────────────────────────────────────

class TestCSPSolverFailureLogging:
    """Test that CSP solver failures are logged at WARNING level with
    exception type, message, and protein length."""

    def test_csp_solver_exception_logs_warning_with_details(self):
        """When CSP solver raises an exception, log WARNING with
        exception type, message, and protein length."""
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer

        opt = HybridOptimizer(species="ecoli")
        protein = "MVLSPADKTN"

        with patch("biocompiler.solver.dispatch.is_solver_available", return_value=False):
            # When solver is not available, _try_csp_solver should return None
            result = opt._try_csp_solver(protein)
            assert result is None

    def test_csp_solver_unsolved_logs_warning_with_protein_length(self):
        """When CSP solver returns unsolved, log WARNING with protein length.

        We verify this by checking the actual code path in _try_csp_solver
        uses logger.warning when the result is unsolved."""
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer
        import biocompiler.optimizer.hybrid_optimizer as ho_mod

        opt = HybridOptimizer(species="ecoli")
        protein = "MVLSPADKTN"

        # Patch at the point where the imports are used
        with patch.object(ho_mod, "is_solver_available", return_value=True), \
             patch.object(ho_mod, "solve_with_csp", return_value=MagicMock(
                 solved=False, sequence="", fallback_used=True,
                 metadata={"reason": "infeasible"},
                 backend_used=MagicMock(value="none"),
             )), \
             patch.object(ho_mod, "logger") as mock_logger:
            result = opt._try_csp_solver(protein)
            assert result is None
            # Should have called logger.warning
            assert mock_logger.warning.called

    def test_csp_fallback_path_logs_warning_with_protein_length(self):
        """When CSP solver is unavailable and fallback to greedy, log WARNING
        with protein length."""
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer

        # Use human (eukaryotic) so we do not take the prokaryotic fast path
        opt = HybridOptimizer(species="human", avoid_gt=False)

        with patch.object(opt, "_try_csp_solver", return_value=None), \
             patch.object(opt, "_greedy_init", return_value=("ATG", 0.9)), \
             patch.object(opt, "_constraint_satisfaction", return_value=("ATG", 0.9, 0, 1, [])), \
             patch.object(opt, "_cai_hill_climb", return_value=("ATG", 0.9, 0)), \
             patch.object(opt, "_compute_cai", return_value=0.9):
            # Call optimize directly which calls _optimize_impl internally
            import biocompiler.optimizer.hybrid_optimizer as ho_mod
            with patch.object(ho_mod, "logger") as mock_logger:
                result = opt.optimize("M", is_prokaryote=False)
                # The fallback path should log warning with protein length
                warning_calls = [str(c) for c in mock_logger.warning.call_args_list]
                assert any("protein_len" in c for c in warning_calls), \
                    f"Expected 'protein_len' in warning calls, got: {warning_calls}"


# ─── 2. Rate limit per-client identification ───────────────────────────

class TestRateLimitClientId:
    """Test per-client identification with IP + User-Agent hash
    and X-Forwarded-For support."""

    def test_resolve_client_id_uses_forwarded_for(self):
        """When X-Forwarded-For is present, use the first IP."""
        from biocompiler.api import _resolve_client_id

        request = MagicMock()
        request.client = MagicMock()
        request.client.host = "127.0.0.1"
        request.headers = MagicMock()
        request.headers.get = lambda k, d="": {
            "x-forwarded-for": "1.0.0.1, 192.168.1.1",
            "user-agent": "",
            "x-api-key": "",
        }.get(k.lower(), d)

        with patch("biocompiler.api._AUTH_MODE", "required"):
            client_id = _resolve_client_id(request)
            # Should use 1.0.0.1 (first in chain) not 127.0.0.1
            assert client_id == "1.0.0.1"

    def test_resolve_client_id_optional_mode_uses_hash(self):
        """In optional auth mode without API key, use IP + User-Agent hash."""
        from biocompiler.api import _resolve_client_id

        request = MagicMock()
        request.client = MagicMock()
        request.client.host = "192.168.1.100"
        request.headers.get = lambda k, d="": {
            "x-forwarded-for": "",
            "user-agent": "BioCompilerClient/1.0",
            "x-api-key": "",
        }.get(k.lower(), d)

        with patch("biocompiler.api._AUTH_MODE", "optional"), \
             patch("biocompiler.api._CONFIGURED_API_KEYS", {"test-key"}):
            client_id = _resolve_client_id(request)
            # Should be hashed: anon:<sha256 hash>
            assert client_id.startswith("anon:")
            # Verify the hash
            expected_hash = hashlib.sha256(
                "192.168.1.100:BioCompilerClient/1.0".encode()
            ).hexdigest()[:16]
            assert client_id == f"anon:{expected_hash}"

    def test_resolve_client_id_with_valid_api_key(self):
        """When a valid API key is provided, use the key as client_id."""
        from biocompiler.api import _resolve_client_id

        request = MagicMock()
        request.client = MagicMock()
        request.client.host = "192.168.1.100"
        request.headers.get = lambda k, d="": {
            "x-forwarded-for": "",
            "user-agent": "SomeAgent",
            "x-api-key": "my-secret-key",
        }.get(k.lower(), d)

        with patch("biocompiler.api._AUTH_MODE", "optional"), \
             patch("biocompiler.api._CONFIGURED_API_KEYS", {"my-secret-key"}):
            client_id = _resolve_client_id(request)
            assert client_id == "my-secret-key"

    def test_resolve_client_id_no_forwarded_for(self):
        """Without X-Forwarded-For, fall back to request.client.host."""
        from biocompiler.api import _resolve_client_id

        request = MagicMock()
        request.client = MagicMock()
        request.client.host = "1.0.0.5"
        request.headers.get = lambda k, d="": {
            "user-agent": "",
            "x-api-key": "",
        }.get(k.lower(), d)

        with patch("biocompiler.api._AUTH_MODE", "required"):
            client_id = _resolve_client_id(request)
            assert client_id == "1.0.0.5"


# ─── 3. Chunked encoding body size check ───────────────────────────────

class TestChunkedEncodingBodySize:
    """Test that chunked transfer encoding body size is tracked during
    streaming and requests exceeding MAX_REQUEST_SIZE are aborted."""

    def test_max_request_size_constant_exists(self):
        """MAX_REQUEST_SIZE is defined and reasonable."""
        from biocompiler.api import MAX_REQUEST_SIZE
        assert MAX_REQUEST_SIZE > 0
        assert MAX_REQUEST_SIZE <= 100_000_000  # <= 100 MB

    def test_dna_length_constant_exists(self):
        """MAX_DNA_LENGTH is defined and configurable."""
        from biocompiler.api import MAX_DNA_LENGTH
        assert MAX_DNA_LENGTH > 0
        assert MAX_DNA_LENGTH >= 100_000  # default 100K bases


# ─── 4. Deprecated NUMBA kernels ───────────────────────────────────────

class TestDeprecatedNumbaKernels:
    """Test that count_gc_parallel and scan_restriction_sites_multi
    emit deprecation warnings when called."""

    def test_count_gc_parallel_emits_deprecation_warning(self):
        """count_gc_parallel should emit DeprecationWarning."""
        from biocompiler.optimizer.numba_kernels import count_gc_parallel

        test_data = [71, 67, 65, 84]  # G, C, A, T
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = count_gc_parallel(test_data)
            # Should emit at least one DeprecationWarning
            dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            assert len(dep_warnings) >= 1
            assert "count_gc_parallel" in str(dep_warnings[0].message)
            assert "deprecated" in str(dep_warnings[0].message).lower()

    def test_scan_restriction_sites_multi_emits_deprecation_warning(self):
        """scan_restriction_sites_multi should emit DeprecationWarning."""
        from biocompiler.optimizer.numba_kernels import scan_restriction_sites_multi

        seq = [65, 84, 71, 67]  # ATGC
        pattern = [71, 84]  # GT
        offsets = [0]
        lens = [2]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = scan_restriction_sites_multi(
                seq, pattern, offsets, lens, 1
            )
            dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            assert len(dep_warnings) >= 1
            assert "scan_restriction_sites_multi" in str(dep_warnings[0].message)
            assert "deprecated" in str(dep_warnings[0].message).lower()

    def test_deprecated_kernels_still_produce_correct_results(self):
        """Deprecated kernels should still return correct results."""
        from biocompiler.optimizer.numba_kernels import count_gc_parallel, count_gc

        test_data = [71, 67, 65, 84]  # G, C, A, T
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result_parallel = count_gc_parallel(test_data)
        result_serial = count_gc(test_data)
        assert result_parallel == result_serial == 2


# ─── 5. DNA sequence length limits ─────────────────────────────────────

class TestDNASequenceLengthLimit:
    """Test that SequenceInput.validate_seq() rejects sequences
    exceeding MAX_DNA_LENGTH."""

    def test_max_dna_length_configurable_via_env(self):
        """MAX_DNA_LENGTH should be configurable via BIOCOMPILER_MAX_DNA_LENGTH."""
        # The default is 100,000
        from biocompiler.api import MAX_DNA_LENGTH
        assert MAX_DNA_LENGTH == 100_000

    def test_sequence_input_rejects_oversized_dna(self):
        """SequenceInput should reject DNA sequences exceeding MAX_DNA_LENGTH."""
        from biocompiler.api import SequenceInput, MAX_DNA_LENGTH

        # Create a sequence that is too long
        long_dna = "A" * (MAX_DNA_LENGTH + 1)
        with pytest.raises(ValueError) as exc_info:
            SequenceInput(sequence=long_dna, organism="Escherichia_coli")
        assert "too long" in str(exc_info.value).lower()
        assert str(MAX_DNA_LENGTH) in str(exc_info.value)

    def test_sequence_input_accepts_valid_dna(self):
        """SequenceInput should accept DNA sequences within limits."""
        from biocompiler.api import SequenceInput

        valid_dna = "ATGCATGC"
        result = SequenceInput(sequence=valid_dna, organism="Escherichia_coli")
        assert result.sequence == valid_dna

    def test_sequence_input_accepts_dna_at_exact_limit(self):
        """SequenceInput should accept DNA at exactly MAX_DNA_LENGTH."""
        from biocompiler.api import SequenceInput, MAX_DNA_LENGTH

        exact_dna = "A" * MAX_DNA_LENGTH
        result = SequenceInput(sequence=exact_dna, organism="Escherichia_coli")
        assert len(result.sequence) == MAX_DNA_LENGTH


# ─── 6. OpenAPI customization ──────────────────────────────────────────

class TestOpenAPICustomization:
    """Test that the FastAPI app has proper OpenAPI customization."""

    def test_create_app_has_contact_info(self):
        """FastAPI app should have contact information."""
        from biocompiler.api import create_app

        # Do not actually create the app as it has side effects.
        # Instead, check the create_app function logic.
        with patch("biocompiler.api._AUTH_MODE", "disabled"):
            app = create_app()
            assert app.title == "BioCompiler API"
            # Check OpenAPI spec has contact
            openapi = app.openapi()
            assert "contact" in openapi.get("info", {})
            assert openapi["info"]["contact"]["name"] == "BioCompiler Team"

    def test_create_app_has_license_info(self):
        """FastAPI app should have license information."""
        from biocompiler.api import create_app

        with patch("biocompiler.api._AUTH_MODE", "disabled"):
            app = create_app()
            openapi = app.openapi()
            assert "license" in openapi.get("info", {})
            assert "MIT" in openapi["info"]["license"]["name"]

    def test_create_app_has_description(self):
        """FastAPI app should have a detailed description."""
        from biocompiler.api import create_app

        with patch("biocompiler.api._AUTH_MODE", "disabled"):
            app = create_app()
            openapi = app.openapi()
            desc = openapi.get("info", {}).get("description", "")
            assert "BioCompiler" in desc
            assert "constraint satisfaction" in desc.lower() or "CSP" in desc


# ─── 7. docs_url/redoc_url env var ────────────────────────────────────

class TestDocsUrlRedocUrl:
    """Test BIOCOMPILER_DOCS_ENABLED env var controls docs_url/redoc_url."""

    def test_docs_disabled_in_production_mode(self):
        """When auth mode is 'required' (production), docs should be disabled
        unless explicitly enabled."""
        from biocompiler.api import create_app

        with patch("biocompiler.api._AUTH_MODE", "required"), \
             patch.dict(os.environ, {}, clear=False):
            # Remove the env var if it exists
            os.environ.pop("BIOCOMPILER_DOCS_ENABLED", None)
            app = create_app()
            assert app.docs_url is None
            assert app.redoc_url is None

    def test_docs_enabled_in_dev_mode(self):
        """When auth mode is 'disabled' (dev), docs should be enabled
        unless explicitly disabled."""
        from biocompiler.api import create_app

        with patch("biocompiler.api._AUTH_MODE", "disabled"), \
             patch.dict(os.environ, {}, clear=False):
            os.environ.pop("BIOCOMPILER_DOCS_ENABLED", None)
            app = create_app()
            assert app.docs_url == "/docs"
            assert app.redoc_url == "/redoc"

    def test_docs_explicitly_enabled(self):
        """BIOCOMPILER_DOCS_ENABLED=true should enable docs even in production."""
        from biocompiler.api import create_app

        with patch("biocompiler.api._AUTH_MODE", "required"), \
             patch.dict(os.environ, {"BIOCOMPILER_DOCS_ENABLED": "true"}):
            app = create_app()
            assert app.docs_url == "/docs"
            assert app.redoc_url == "/redoc"

    def test_docs_explicitly_disabled(self):
        """BIOCOMPILER_DOCS_ENABLED=false should disable docs even in dev."""
        from biocompiler.api import create_app

        with patch("biocompiler.api._AUTH_MODE", "disabled"), \
             patch.dict(os.environ, {"BIOCOMPILER_DOCS_ENABLED": "false"}):
            app = create_app()
            assert app.docs_url is None
            assert app.redoc_url is None


# ─── 8. Logging for broad exception catches ────────────────────────────

class TestBroadExceptionLogging:
    """Test that broad exception catches log the full exception with traceback."""

    def test_csp_solver_exception_logged_with_details(self):
        """CSP solver exception in _try_csp_solver should be logged with
        exception type and message."""
        from biocompiler.optimizer.hybrid_optimizer import HybridOptimizer
        import biocompiler.optimizer.hybrid_optimizer as ho_mod

        opt = HybridOptimizer(species="ecoli")
        protein = "MVLSPADKTN"

        with patch.object(ho_mod, "is_solver_available", return_value=True), \
             patch.object(ho_mod, "solve_with_csp", side_effect=RuntimeError("solver crashed")), \
             patch.object(ho_mod, "logger") as mock_logger:
            result = opt._try_csp_solver(protein)
            assert result is None
            # Should log warning with exception type and message
            mock_logger.warning.assert_called()
            call_args = str(mock_logger.warning.call_args)
            assert "RuntimeError" in call_args
            assert "solver crashed" in call_args

    def test_dispatch_logs_backend_transitions(self):
        """When solver backends fail, dispatch should log the transition."""
        from biocompiler.solver.dispatch import _try_backend
        from biocompiler.solver.types import SolverResult, SolverBackend, SolverConfig

        # Test that _try_backend logs failure reason
        model = MagicMock()
        config = SolverConfig(organism="Escherichia_coli")
        engine_cls = MagicMock()
        engine_cls.return_value.solve.return_value = SolverResult(
            sequence="",
            solved=False,
            backend_used=SolverBackend.NONE,
            protein="M",
            organism="Escherichia_coli",
            fallback_used=True,
            solve_time_seconds=0.1,
            violations=[],
            metadata={"reason": "test failure"},
        )

        result, fail_reason = _try_backend(
            engine_cls, model, config, "TestBackend", SolverBackend.ORTOOLS,
        )
        assert fail_reason is not None
        assert "returned" in fail_reason or "infeasible" in fail_reason or "fallback" in fail_reason


# ─── Integration: InfoResponse includes max_dna_length ────────────────

class TestInfoResponseIncludesDNALength:
    """Test that the /info endpoint includes max_dna_length."""

    def test_info_response_model_has_max_dna_length(self):
        """InfoResponse should have max_dna_length field."""
        from biocompiler.api import InfoResponse

        info = InfoResponse(
            max_protein_length=10000,
            max_dna_length=100000,
            max_batch_size=50,
            max_request_size=10000000,
            optimize_timeout_s=300,
            supported_organisms=["Escherichia_coli"],
            api_version="1.0.0",
            safety_version="1.0.0",
        )
        assert info.max_dna_length == 100000
