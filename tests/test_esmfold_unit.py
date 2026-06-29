"""Unit tests for biocompiler.esmfold — pure-function tests that do NOT require ESMFold.

These tests exercise the module's public API, input validation, fallback
behaviour, and named constants directly, without requiring the ESMFold
library, the ESM Atlas API, or a local ``esm`` package installation.

Test categories
---------------
1. Named constants — correct values, types, and invariants
2. Function existence and signatures — all public symbols exist and are callable
3. Input validation — _validate_protein rejects invalid/empty sequences
4. Fallback behaviour — predict_structure degrades gracefully offline
5. ESMFoldResult dataclass — construction, aliases, confidence_level
6. ESMFoldCache — put/get/eviction/hit_rate/clear
7. classify_plddt — boundary values for all four confidence bands
8. parse_pdb — minimal PDB string parsing
9. Batch helpers — validate_batch_input, estimate_batch_time
"""
from __future__ import annotations

import inspect
from unittest.mock import patch

import pytest

from biocompiler.engines.esmfold import (
    ESMFoldCache,
    ESMFoldError,
    ESMFoldResult,
    ESMFOLD_PLDDT_CORRELATION,
    MAX_BATCH_SIZE,
    MAX_PROTEIN_LENGTH,
    MAX_RETRIES,
    RETRY_BASE_DELAY,
    DEFAULT_API_URL,
    DEFAULT_TIMEOUT,
    STANDARD_AMINO_ACIDS,
    BatchStructureRequest,
    BatchStructureResult,
    classify_plddt,
    clear_cache,
    estimate_batch_time,
    estimate_contact_map,
    is_esmfold_available,
    parse_pdb,
    predict_structure,
    predict_structure_batch,
    analyze_structure,
    compute_backbone_dihedrals,
    validate_batch_input,
    format_batch_report,
    _validate_protein,
    _build_result_from_pdb,
    _get_default_cache,
)
from biocompiler.shared.exceptions import EngineError

# Minimal PDB fixture used by parse_pdb / contact map / dihedral tests
MINI_PDB = (
    "ATOM      1  CA  MET A   1       1.000   2.000   3.000  1.00 85.00           C\n"
    "ATOM      2  CA  ALA A   2       4.000   2.000   3.000  1.00 90.00           C\n"
    "ATOM      3  CA  GLY A   3       7.000   2.000   3.000  1.00 78.00           C\n"
    "END\n"
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Named constants
# ═══════════════════════════════════════════════════════════════════════════

class TestNamedConstants:
    """Verify all exported constants have expected values and types."""

    def test_standard_amino_acids_is_20(self) -> None:
        """STANDARD_AMINO_ACIDS contains exactly the 20 canonical single-letter codes."""
        assert isinstance(STANDARD_AMINO_ACIDS, set)
        assert len(STANDARD_AMINO_ACIDS) == 20
        expected = set("ACDEFGHIKLMNPQRSTVWY")
        assert STANDARD_AMINO_ACIDS == expected

    def test_esmfold_plddt_correlation_value(self) -> None:
        """ESMFOLD_PLDDT_CORRELATION should be 0.8 (Lin et al., Science 2023)."""
        assert isinstance(ESMFOLD_PLDDT_CORRELATION, float)
        assert ESMFOLD_PLDDT_CORRELATION == pytest.approx(0.8)

    def test_esmfold_plddt_correlation_in_unit_interval(self) -> None:
        """Pearson correlation must be in [0, 1]."""
        assert 0.0 <= ESMFOLD_PLDDT_CORRELATION <= 1.0

    def test_default_api_url(self) -> None:
        """DEFAULT_API_URL is the ESM Atlas endpoint."""
        assert isinstance(DEFAULT_API_URL, str)
        assert DEFAULT_API_URL == "https://api.esmatlas.com/fetchPredictedStructure"

    def test_default_timeout_positive(self) -> None:
        """DEFAULT_TIMEOUT must be a positive float."""
        assert isinstance(DEFAULT_TIMEOUT, float)
        assert DEFAULT_TIMEOUT > 0.0

    def test_max_retries_positive_int(self) -> None:
        """MAX_RETRIES must be a positive integer."""
        assert isinstance(MAX_RETRIES, int)
        assert MAX_RETRIES >= 1

    def test_retry_base_delay_positive(self) -> None:
        """RETRY_BASE_DELAY must be a positive number of seconds."""
        assert isinstance(RETRY_BASE_DELAY, float)
        assert RETRY_BASE_DELAY > 0.0

    def test_max_batch_size_positive(self) -> None:
        """MAX_BATCH_SIZE must be a positive integer."""
        assert isinstance(MAX_BATCH_SIZE, int)
        assert MAX_BATCH_SIZE >= 1

    def test_max_protein_length_positive(self) -> None:
        """MAX_PROTEIN_LENGTH must be a positive integer."""
        assert isinstance(MAX_PROTEIN_LENGTH, int)
        assert MAX_PROTEIN_LENGTH >= 1

    def test_max_batch_size_is_50(self) -> None:
        """MAX_BATCH_SIZE is 50 (documented default)."""
        assert MAX_BATCH_SIZE == 50

    def test_max_protein_length_is_1000(self) -> None:
        """MAX_PROTEIN_LENGTH is 1000 (documented default)."""
        assert MAX_PROTEIN_LENGTH == 1000

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

    def test_predict_structure_is_callable(self) -> None:
        """predict_structure is a callable function."""
        assert callable(predict_structure)

    def test_predict_structure_signature(self) -> None:
        """predict_structure accepts (protein, organism, use_api, api_url, timeout)."""
        sig = inspect.signature(predict_structure)
        params = list(sig.parameters.keys())
        assert "protein" in params
        assert "organism" in params
        assert "use_api" in params
        assert "api_url" in params
        assert "timeout" in params

    def test_predict_structure_batch_is_callable(self) -> None:
        """predict_structure_batch is a callable function."""
        assert callable(predict_structure_batch)

    def test_predict_structure_batch_signature(self) -> None:
        """predict_structure_batch accepts (sequences, max_concurrent, ...)."""
        sig = inspect.signature(predict_structure_batch)
        params = list(sig.parameters.keys())
        assert "sequences" in params
        assert "max_concurrent" in params
        assert "use_api" in params

    def test_analyze_structure_is_callable(self) -> None:
        """analyze_structure is a callable function."""
        assert callable(analyze_structure)

    def test_analyze_structure_signature(self) -> None:
        """analyze_structure accepts (protein, use_api, api_url, timeout)."""
        sig = inspect.signature(analyze_structure)
        params = list(sig.parameters.keys())
        assert "protein" in params
        assert "use_api" in params

    # --- Utility functions ---

    def test_is_esmfold_available_is_callable(self) -> None:
        """is_esmfold_available is a callable function."""
        assert callable(is_esmfold_available)

    def test_parse_pdb_is_callable(self) -> None:
        """parse_pdb is a callable function."""
        assert callable(parse_pdb)

    def test_compute_backbone_dihedrals_is_callable(self) -> None:
        """compute_backbone_dihedrals is a callable function."""
        assert callable(compute_backbone_dihedrals)

    def test_classify_plddt_is_callable(self) -> None:
        """classify_plddt is a callable function."""
        assert callable(classify_plddt)

    def test_estimate_contact_map_is_callable(self) -> None:
        """estimate_contact_map is a callable function."""
        assert callable(estimate_contact_map)

    def test_validate_batch_input_is_callable(self) -> None:
        """validate_batch_input is a callable function."""
        assert callable(validate_batch_input)

    def test_estimate_batch_time_is_callable(self) -> None:
        """estimate_batch_time is a callable function."""
        assert callable(estimate_batch_time)

    def test_format_batch_report_is_callable(self) -> None:
        """format_batch_report is a callable function."""
        assert callable(format_batch_report)

    def test_clear_cache_is_callable(self) -> None:
        """clear_cache is a callable function."""
        assert callable(clear_cache)

    # --- Signature parameter counts ---

    def test_classify_plddt_signature(self) -> None:
        """classify_plddt accepts a single float parameter."""
        sig = inspect.signature(classify_plddt)
        assert list(sig.parameters.keys()) == ["mean_plddt"]

    def test_parse_pdb_signature(self) -> None:
        """parse_pdb accepts a single string parameter."""
        sig = inspect.signature(parse_pdb)
        assert list(sig.parameters.keys()) == ["pdb_string"]

    def test_estimate_contact_map_signature(self) -> None:
        """estimate_contact_map accepts pdb_string and distance_threshold."""
        sig = inspect.signature(estimate_contact_map)
        params = list(sig.parameters.keys())
        assert "pdb_string" in params
        assert "distance_threshold" in params

    def test_validate_batch_input_signature(self) -> None:
        """validate_batch_input accepts a list of proteins."""
        sig = inspect.signature(validate_batch_input)
        assert list(sig.parameters.keys()) == ["proteins"]

    def test_estimate_batch_time_signature(self) -> None:
        """estimate_batch_time accepts (num_proteins, avg_length, concurrent)."""
        sig = inspect.signature(estimate_batch_time)
        params = list(sig.parameters.keys())
        assert "num_proteins" in params
        assert "avg_length" in params
        assert "concurrent" in params


# ═══════════════════════════════════════════════════════════════════════════
# 3. Input validation (_validate_protein)
# ═══════════════════════════════════════════════════════════════════════════

class TestValidateProtein:
    """Test _validate_protein — valid/invalid/empty sequences."""

    def test_valid_protein_passes(self) -> None:
        """Valid protein sequence does not raise."""
        _validate_protein("ACDEFGHIKLMNPQRSTVWY")  # should not raise
        # Implicit assertion: no exception raised

    def test_short_valid_protein(self) -> None:
        """Short but valid protein sequence does not raise."""
        _validate_protein("MAG")
        # Implicit assertion: no exception raised

    def test_single_residue_passes(self) -> None:
        """Single standard amino acid does not raise."""
        _validate_protein("M")
        # Implicit assertion: no exception raised

    def test_empty_protein_raises(self) -> None:
        """Empty protein string raises ESMFoldError."""
        with pytest.raises(ESMFoldError):
            _validate_protein("")

    def test_whitespace_only_raises(self) -> None:
        """Whitespace-only protein string raises ESMFoldError."""
        with pytest.raises(ESMFoldError):
            _validate_protein("   ")

    def test_non_standard_aa_raises(self) -> None:
        """Non-standard amino acid (X) raises ESMFoldError."""
        with pytest.raises(ESMFoldError):
            _validate_protein("MXG")

    def test_digit_in_sequence_raises(self) -> None:
        """Digit in sequence raises ESMFoldError."""
        with pytest.raises(ESMFoldError):
            _validate_protein("MA1G")

    def test_lowercase_protein_raises(self) -> None:
        """Lowercase amino acids are NOT accepted by _validate_protein.

        The underlying engine_base.validate_protein_sequence uppercases
        and then checks, so lowercase is accepted after normalization.
        This test verifies the behaviour (it should pass after uppercasing).
        """
        # The implementation uppercases then checks, so lowercase is accepted
        _validate_protein("mag")  # should not raise
        # Implicit assertion: no exception raised

    def test_esmfold_error_is_engine_error(self) -> None:
        """ESMFoldError should be a subclass of EngineError."""
        assert issubclass(ESMFoldError, EngineError)

    def test_esmfold_error_has_reason_and_protein(self) -> None:
        """ESMFoldError stores reason and protein attributes.

        Note: EngineError.__init__ overwrites self.reason with the
        formatted message, so err.reason is the full message string.
        """
        err = ESMFoldError("test reason", protein="MXG")
        assert "test reason" in str(err)
        assert err.protein == "MXG"

    def test_esmfold_error_without_protein(self) -> None:
        """ESMFoldError can be created without a protein."""
        err = ESMFoldError("test reason")
        assert "test reason" in str(err)
        assert err.protein is None

    def test_validation_error_chained_from_value_error(self) -> None:
        """ESMFoldError from _validate_protein chains from ValueError."""
        with pytest.raises(ESMFoldError) as exc_info:
            _validate_protein("MXG")
        assert exc_info.value.__cause__ is not None
        assert isinstance(exc_info.value.__cause__, ValueError)


# ═══════════════════════════════════════════════════════════════════════════
# 4. Fallback behaviour (predict_structure offline)
# ═══════════════════════════════════════════════════════════════════════════

class TestFallbackBehaviour:
    """Test predict_structure fallback when API and local esm are unavailable."""

    @patch("biocompiler.engines.esmfold._predict_via_api", return_value=None)
    @patch("biocompiler.engines.esmfold._predict_via_local_esm", return_value=None)
    def test_offline_returns_heuristic_fallback(self, mock_local, mock_api) -> None:
        """When API and local esm are unavailable, heuristic fallback is used."""
        result = predict_structure("MAG", use_api=False)
        # The heuristic fallback should be invoked
        assert result.success is True
        assert result.method == "heuristic_fallback"
        assert result.sequence == "MAG"
        assert len(result.plddt_scores) == 3
        assert result.mean_plddt > 0.0

    @patch("biocompiler.engines.esmfold._predict_via_api", return_value=None)
    @patch("biocompiler.engines.esmfold._predict_via_local_esm", return_value=None)
    def test_heuristic_fallback_no_pdb_string(self, mock_local, mock_api) -> None:
        """Heuristic fallback results have no PDB string."""
        result = predict_structure("MAG", use_api=False)
        assert result.pdb_string == ""

    @patch("biocompiler.engines.esmfold._predict_via_api", return_value=None)
    @patch("biocompiler.engines.esmfold._predict_via_local_esm", return_value=None)
    def test_heuristic_fallback_confidence_level_very_low(self, mock_local, mock_api) -> None:
        """Heuristic fallback results have confidence_level='very_low'."""
        result = predict_structure("MAG", use_api=False)
        assert result.confidence_level == "very_low"

    @patch("biocompiler.engines.esmfold._predict_via_api", return_value=None)
    @patch("biocompiler.engines.esmfold._predict_via_local_esm", return_value=None)
    def test_invalid_protein_raises_even_offline(self, mock_local, mock_api) -> None:
        """Invalid protein raises ESMFoldError even in offline mode (no retry)."""
        with pytest.raises(ESMFoldError):
            predict_structure("MXG", use_api=False)

    @patch("biocompiler.engines.esmfold._predict_via_api", return_value=None)
    @patch("biocompiler.engines.esmfold._predict_via_local_esm", return_value=None)
    def test_empty_protein_raises_even_offline(self, mock_local, mock_api) -> None:
        """Empty protein raises ESMFoldError even in offline mode."""
        with pytest.raises(ESMFoldError):
            predict_structure("", use_api=False)

    def test_predict_structure_accepts_use_api_false(self) -> None:
        """predict_structure can be called with use_api=False to skip API."""
        # Just verify it does not crash — it will use heuristic fallback
        result = predict_structure("MAG", use_api=False)
        assert isinstance(result, ESMFoldResult)
        assert result.sequence == "MAG"

    @patch("biocompiler.engines.esmfold._predict_via_api", return_value=None)
    @patch("biocompiler.engines.esmfold._predict_via_local_esm", return_value=None)
    def test_predict_structure_accepts_organism(self, mock_local, mock_api) -> None:
        """predict_structure accepts an organism parameter without error."""
        result = predict_structure("MAG", organism="Homo_sapiens", use_api=False)
        assert result.sequence == "MAG"


# ═══════════════════════════════════════════════════════════════════════════
# 5. ESMFoldResult dataclass
# ═══════════════════════════════════════════════════════════════════════════

class TestESMFoldResult:
    """Test ESMFoldResult construction, defaults, aliases, and confidence_level."""

    def test_default_construction(self) -> None:
        """ESMFoldResult can be constructed with all defaults."""
        result = ESMFoldResult()
        assert result.sequence == ""
        assert result.primary_score == 0.0
        assert result.success is True
        assert result.error is None
        assert result.execution_time_s == 0.0
        assert result.engine_name == "esmfold"
        assert result.primary_score_label == "pLDDT"
        assert result.plddt_scores == []
        assert result.pae_matrix is None
        assert result.pdb_string == ""
        assert result.model_name == "esmfold_v1"
        assert result.method == "esmfold_api"

    def test_construction_with_protein_alias(self) -> None:
        """protein alias maps to sequence field."""
        result = ESMFoldResult(protein="MAG")
        assert result.sequence == "MAG"
        assert result.protein == "MAG"

    def test_construction_with_mean_plddt_alias(self) -> None:
        """mean_plddt alias maps to primary_score field."""
        result = ESMFoldResult(mean_plddt=85.0)
        assert result.primary_score == pytest.approx(85.0)
        assert result.mean_plddt == pytest.approx(85.0)

    def test_construction_with_plddt_alias(self) -> None:
        """plddt alias maps to primary_score field."""
        result = ESMFoldResult(plddt=90.0)
        assert result.primary_score == pytest.approx(90.0)
        assert result.plddt == pytest.approx(90.0)

    def test_construction_with_confidence_class_alias(self) -> None:
        """confidence_class alias maps to classification field."""
        result = ESMFoldResult(confidence_class="Confident")
        assert result.classification == "Confident"
        assert result.confidence_class == "Confident"

    def test_auto_classification_on_construction(self) -> None:
        """When classification is not given but primary_score > 0, auto-compute."""
        result = ESMFoldResult(mean_plddt=85.0)
        assert result.classification != ""
        # 85.0 should classify as "Confident" (between 70+eps and 90+eps)
        assert result.classification == "Confident"

    def test_protein_property_setter(self) -> None:
        """protein property setter updates sequence."""
        result = ESMFoldResult()
        result.protein = "KLV"
        assert result.sequence == "KLV"
        assert result.protein == "KLV"

    def test_mean_plddt_property_setter(self) -> None:
        """mean_plddt property setter updates primary_score."""
        result = ESMFoldResult()
        result.mean_plddt = 75.0
        assert result.primary_score == pytest.approx(75.0)

    def test_plddt_property_setter(self) -> None:
        """plddt property setter updates primary_score."""
        result = ESMFoldResult()
        result.plddt = 80.0
        assert result.primary_score == pytest.approx(80.0)

    def test_confidence_class_property_setter(self) -> None:
        """confidence_class property setter updates classification."""
        result = ESMFoldResult()
        result.confidence_class = "Very high (experimental)"
        assert result.classification == "Very high (experimental)"

    def test_pae_alias(self) -> None:
        """pae property is alias for pae_matrix."""
        matrix = [[1.0, 2.0], [3.0, 4.0]]
        result = ESMFoldResult(pae_matrix=matrix)
        assert result.pae == matrix
        # Setter
        result.pae = [[5.0, 6.0]]
        assert result.pae_matrix == [[5.0, 6.0]]

    # --- confidence_level property ---

    def test_confidence_level_failed(self) -> None:
        """Failed prediction returns confidence_level='none'."""
        result = ESMFoldResult(success=False)
        assert result.confidence_level == "none"

    def test_confidence_level_heuristic_fallback(self) -> None:
        """Heuristic fallback always returns confidence_level='very_low'."""
        result = ESMFoldResult(mean_plddt=50.0, method="heuristic_fallback")
        assert result.confidence_level == "very_low"

    def test_confidence_level_high(self) -> None:
        """Mean pLDDT >= 70 with API method returns confidence_level='high'."""
        result = ESMFoldResult(mean_plddt=85.0, method="esmfold_api")
        assert result.confidence_level == "high"

    def test_confidence_level_medium(self) -> None:
        """Mean pLDDT 50-70 with API method returns confidence_level='medium'."""
        result = ESMFoldResult(mean_plddt=60.0, method="esmfold_api")
        assert result.confidence_level == "medium"

    def test_confidence_level_low(self) -> None:
        """Mean pLDDT < 50 with API method returns confidence_level='low'."""
        result = ESMFoldResult(mean_plddt=30.0, method="esmfold_api")
        assert result.confidence_level == "low"

    def test_confidence_level_local_method(self) -> None:
        """Local esm method uses same thresholds as API."""
        result = ESMFoldResult(mean_plddt=85.0, method="esmfold_local")
        assert result.confidence_level == "high"

    def test_plddt_scores_default_empty_list(self) -> None:
        """plddt_scores defaults to empty list, not None."""
        result = ESMFoldResult()
        assert result.plddt_scores == []

    def test_plddt_scores_explicit(self) -> None:
        """plddt_scores can be set explicitly."""
        scores = [85.0, 90.0, 78.0]
        result = ESMFoldResult(plddt_scores=scores)
        assert result.plddt_scores == scores

    def test_passed_property(self) -> None:
        """passed property reflects success field."""
        assert ESMFoldResult(success=True).passed is True
        assert ESMFoldResult(success=False).passed is False


# ═══════════════════════════════════════════════════════════════════════════
# 6. ESMFoldCache
# ═══════════════════════════════════════════════════════════════════════════

class TestESMFoldCache:
    """Test ESMFoldCache — put/get/eviction/hit_rate/clear."""

    def _make_result(self, protein: str = "MAG", mean_plddt: float = 84.33) -> ESMFoldResult:
        """Helper to create a minimal ESMFoldResult for cache tests."""
        return ESMFoldResult(
            protein=protein,
            mean_plddt=mean_plddt,
            plddt_scores=[85.0, 90.0, 78.0],
            pdb_string=MINI_PDB,
            success=True,
        )

    def test_create_memory_only(self) -> None:
        """ESMFoldCache can be created without a cache directory."""
        cache = ESMFoldCache()
        assert cache.size == 0

    def test_put_and_get(self) -> None:
        """Put a result and retrieve it."""
        cache = ESMFoldCache()
        result = self._make_result()
        cache.put("MAG", result)
        retrieved = cache.get("MAG")
        assert retrieved is not None
        assert retrieved.sequence == "MAG"
        assert retrieved.mean_plddt == pytest.approx(84.33, abs=0.01)

    def test_get_miss_returns_none(self) -> None:
        """Getting a non-existent key returns None."""
        cache = ESMFoldCache()
        assert cache.get("NONEXISTENT") is None

    def test_eviction_when_full(self) -> None:
        """FIFO eviction when max_size is exceeded."""
        cache = ESMFoldCache(max_size=3)
        for i in range(5):
            cache.put(f"SEQ{i}", self._make_result(protein=f"SEQ{i}"))
        # Only last 3 should remain
        assert cache.get("SEQ0") is None
        assert cache.get("SEQ1") is None
        assert cache.get("SEQ2") is not None
        assert cache.get("SEQ3") is not None
        assert cache.get("SEQ4") is not None

    def test_hit_rate_tracking(self) -> None:
        """hit_rate reflects hits / (hits + misses)."""
        cache = ESMFoldCache()
        cache.put("MAG", self._make_result())
        cache.get("MAG")   # hit
        cache.get("MISS")  # miss
        assert cache.hits == 1
        assert cache.misses == 1
        assert cache.hit_rate == pytest.approx(0.5, abs=0.01)

    def test_hit_rate_no_accesses(self) -> None:
        """hit_rate is 0.0 when there have been no accesses."""
        cache = ESMFoldCache()
        assert cache.hit_rate == 0.0

    def test_clear_resets_cache_and_stats(self) -> None:
        """clear() resets size, hits, and misses."""
        cache = ESMFoldCache()
        cache.put("MAG", self._make_result())
        cache.get("MAG")
        cache.clear()
        assert cache.size == 0
        assert cache.hits == 0
        assert cache.misses == 0

    def test_cache_key_deterministic(self) -> None:
        """Same protein sequence always produces the same cache key."""
        key1 = ESMFoldCache._key("MAG")
        key2 = ESMFoldCache._key("MAG")
        assert key1 == key2

    def test_cache_key_length(self) -> None:
        """Cache key is 16 hex characters (from SHA-256 digest)."""
        key = ESMFoldCache._key("MAG")
        assert len(key) == 16
        assert all(c in "0123456789abcdef" for c in key)

    def test_different_sequences_different_keys(self) -> None:
        """Different protein sequences produce different cache keys."""
        key1 = ESMFoldCache._key("MAG")
        key2 = ESMFoldCache._key("KLV")
        assert key1 != key2


# ═══════════════════════════════════════════════════════════════════════════
# 7. classify_plddt
# ═══════════════════════════════════════════════════════════════════════════

class TestClassifyPlddt:
    """Test classify_plddt boundary values for all four confidence bands."""

    def test_very_high_band(self) -> None:
        """pLDDT >= 90.1 → 'Very high (experimental)'."""
        assert classify_plddt(95.0) == "Very high (experimental)"
        assert classify_plddt(90.1) == "Very high (experimental)"

    def test_confident_band(self) -> None:
        """pLDDT 70.1–90.0 → 'Confident'."""
        assert classify_plddt(80.0) == "Confident"
        assert classify_plddt(70.1) == "Confident"
        assert classify_plddt(90.0) == "Confident"

    def test_low_confidence_band(self) -> None:
        """pLDDT 50.1–70.0 → 'Low confidence'."""
        assert classify_plddt(60.0) == "Low confidence"
        assert classify_plddt(50.1) == "Low confidence"
        assert classify_plddt(70.0) == "Low confidence"

    def test_very_low_band(self) -> None:
        """pLDDT < 50.1 → 'Very low'."""
        assert classify_plddt(40.0) == "Very low"
        assert classify_plddt(50.0) == "Very low"
        assert classify_plddt(0.0) == "Very low"

    def test_boundary_90(self) -> None:
        """Exact 90.0 is 'Confident' (not > 90)."""
        assert classify_plddt(90.0) == "Confident"

    def test_boundary_70(self) -> None:
        """Exact 70.0 is 'Low confidence' (not > 70)."""
        assert classify_plddt(70.0) == "Low confidence"

    def test_boundary_50(self) -> None:
        """Exact 50.0 is 'Very low' (not > 50)."""
        assert classify_plddt(50.0) == "Very low"

    def test_returns_string(self) -> None:
        """classify_plddt always returns a string."""
        for score in [-1.0, 0.0, 25.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0]:
            assert isinstance(classify_plddt(score), str)


# ═══════════════════════════════════════════════════════════════════════════
# 8. parse_pdb
# ═══════════════════════════════════════════════════════════════════════════

class TestParsePdb:
    """Test parse_pdb with a minimal PDB string."""

    def test_returns_dict_with_required_keys(self) -> None:
        """parse_pdb returns dict with atoms, residues, chains, plddt_scores."""
        parsed = parse_pdb(MINI_PDB)
        assert isinstance(parsed, dict)
        assert "atoms" in parsed
        assert "residues" in parsed
        assert "chains" in parsed
        assert "plddt_scores" in parsed

    def test_atoms_count(self) -> None:
        """3 CA atoms in MINI_PDB."""
        parsed = parse_pdb(MINI_PDB)
        assert len(parsed["atoms"]) == 3

    def test_residues_count(self) -> None:
        """3 residues in MINI_PDB."""
        parsed = parse_pdb(MINI_PDB)
        assert len(parsed["residues"]) == 3

    def test_chains(self) -> None:
        """Chain A is present."""
        parsed = parse_pdb(MINI_PDB)
        assert "A" in parsed["chains"]

    def test_plddt_scores(self) -> None:
        """pLDDT scores extracted from B-factors."""
        parsed = parse_pdb(MINI_PDB)
        assert len(parsed["plddt_scores"]) == 3
        assert parsed["plddt_scores"] == pytest.approx([85.0, 90.0, 78.0], abs=0.5)

    def test_empty_pdb(self) -> None:
        """Empty PDB string returns empty structures."""
        parsed = parse_pdb("")
        assert parsed["atoms"] == []
        assert parsed["residues"] == []
        assert parsed["chains"] == []
        assert parsed["plddt_scores"] == []

    def test_non_atom_lines_skipped(self) -> None:
        """Non-ATOM lines (HEADER, TER, END, etc.) are skipped."""
        pdb = "HEADER    TEST\nEND\n"
        parsed = parse_pdb(pdb)
        assert parsed["atoms"] == []


# ═══════════════════════════════════════════════════════════════════════════
# 9. Batch helpers
# ═══════════════════════════════════════════════════════════════════════════

class TestValidateBatchInput:
    """Test validate_batch_input for valid and invalid protein lists."""

    def test_valid_batch_no_errors(self) -> None:
        """Valid list of proteins returns empty error list."""
        errors = validate_batch_input(["MKWVTFISLLFLFSSAYS", "MAG"])
        assert isinstance(errors, list)
        assert len(errors) == 0

    def test_batch_too_large(self) -> None:
        """More than MAX_BATCH_SIZE proteins produces an error."""
        errors = validate_batch_input(["MAG"] * (MAX_BATCH_SIZE + 1))
        assert len(errors) > 0

    def test_protein_too_long(self) -> None:
        """Protein exceeding MAX_PROTEIN_LENGTH produces an error."""
        errors = validate_batch_input(["A" * (MAX_PROTEIN_LENGTH + 1)])
        assert len(errors) > 0

    def test_non_standard_aa(self) -> None:
        """Protein with non-standard AA produces an error."""
        errors = validate_batch_input(["MXG"])
        assert len(errors) > 0

    def test_empty_batch_valid(self) -> None:
        """Empty batch list returns no errors (trivially valid)."""
        errors = validate_batch_input([])
        assert errors == []

    def test_single_valid_protein(self) -> None:
        """Single valid protein returns no errors."""
        errors = validate_batch_input(["MAG"])
        assert errors == []


class TestEstimateBatchTime:
    """Test estimate_batch_time for correct estimation logic."""

    def test_positive_inputs(self) -> None:
        """Positive inputs return a positive float."""
        t = estimate_batch_time(5, 100)
        assert isinstance(t, float)
        assert t > 0

    def test_zero_proteins(self) -> None:
        """Zero proteins returns 0.0."""
        assert estimate_batch_time(0, 100) == 0.0

    def test_zero_avg_length(self) -> None:
        """Zero avg_length returns 0.0."""
        assert estimate_batch_time(5, 0) == 0.0

    def test_negative_proteins(self) -> None:
        """Negative num_proteins returns 0.0."""
        assert estimate_batch_time(-5, 100) == 0.0

    def test_more_concurrency_reduces_time(self) -> None:
        """Higher concurrency should reduce estimated time."""
        t1 = estimate_batch_time(10, 100, concurrent=1)
        t2 = estimate_batch_time(10, 100, concurrent=5)
        assert t2 < t1


class TestBatchStructureRequest:
    """Test BatchStructureRequest dataclass."""

    def test_default_values(self) -> None:
        """BatchStructureRequest defaults: names=None, use_cache=True, max_concurrent=3."""
        req = BatchStructureRequest(proteins=["MAG"])
        assert req.names is None
        assert req.use_cache is True
        assert req.max_concurrent == 3
        assert req.timeout_per_protein == 120.0
        assert req.stop_on_failure is False

    def test_custom_values(self) -> None:
        """BatchStructureRequest accepts custom values."""
        req = BatchStructureRequest(
            proteins=["MAG", "KLV"],
            names=["p1", "p2"],
            use_cache=False,
            max_concurrent=5,
            timeout_per_protein=60.0,
            stop_on_failure=True,
        )
        assert req.proteins == ["MAG", "KLV"]
        assert req.names == ["p1", "p2"]
        assert req.use_cache is False
        assert req.max_concurrent == 5
        assert req.timeout_per_protein == 60.0
        assert req.stop_on_failure is True


class TestBatchStructureResult:
    """Test BatchStructureResult dataclass."""

    def test_creation(self) -> None:
        """BatchStructureResult can be created with all fields."""
        result = BatchStructureResult(
            results=[{"name": "p1", "status": "success"}],
            names=["p1"],
            total=1,
            successful=1,
            failed=0,
            from_cache=0,
            total_time_s=1.5,
            summary={"total": 1, "successful": 1, "failed": 0},
        )
        assert result.total == 1
        assert result.successful == 1
        assert result.failed == 0


# ═══════════════════════════════════════════════════════════════════════════
# 10. Module-level helpers
# ═══════════════════════════════════════════════════════════════════════════

class TestModuleHelpers:
    """Test module-level helper functions."""

    def test_clear_cache_no_error(self) -> None:
        """clear_cache() runs without error."""
        clear_cache()  # should not raise
        # Implicit assertion: no exception raised

    def test_get_default_cache_returns_cache(self) -> None:
        """_get_default_cache returns an ESMFoldCache instance."""
        cache = _get_default_cache()
        assert isinstance(cache, ESMFoldCache)

    def test_is_esmfold_available_returns_bool(self) -> None:
        """is_esmfold_available returns a boolean."""
        result = is_esmfold_available()
        assert isinstance(result, bool)

    def test_build_result_from_pdb(self) -> None:
        """_build_result_from_pdb constructs ESMFoldResult from PDB string."""
        result = _build_result_from_pdb("MAG", MINI_PDB, "esmfold_v1")
        assert result.success is True
        assert result.sequence == "MAG"
        assert result.model_name == "esmfold_v1"
        assert result.pdb_string == MINI_PDB
        assert len(result.plddt_scores) == 3
        # mean of [85, 90, 78] ≈ 84.33
        assert result.mean_plddt == pytest.approx(84.33, abs=0.01)

    def test_build_result_from_pdb_method_kwarg(self) -> None:
        """_build_result_from_pdb method kwarg is propagated."""
        result = _build_result_from_pdb("MAG", MINI_PDB, "esmfold_v1", method="esmfold_local")
        assert result.method == "esmfold_local"


# ═══════════════════════════════════════════════════════════════════════════
# 11. Contact map and dihedral edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestEstimateContactMap:
    """Test estimate_contact_map symmetry and basic properties."""

    def test_symmetric_matrix(self) -> None:
        """Contact map is symmetric."""
        cmap = estimate_contact_map(MINI_PDB, distance_threshold=5.0)
        for i in range(len(cmap)):
            for j in range(len(cmap)):
                assert cmap[i][j] == cmap[j][i]

    def test_diagonal_zero(self) -> None:
        """Contact map diagonal is all zeros."""
        cmap = estimate_contact_map(MINI_PDB)
        for i in range(len(cmap)):
            assert cmap[i][i] == 0

    def test_empty_pdb_empty_contact_map(self) -> None:
        """Empty PDB returns empty contact map."""
        cmap = estimate_contact_map("")
        assert cmap == []

    def test_default_threshold_is_8(self) -> None:
        """Default distance_threshold is 8.0 Angstroms."""
        sig = inspect.signature(estimate_contact_map)
        assert sig.parameters["distance_threshold"].default == 8.0


class TestComputeBackboneDihedrals:
    """Test compute_backbone_dihedrals returns correct structure."""

    def test_returns_phi_and_psi(self) -> None:
        """Returns dict with 'phi' and 'psi' keys."""
        dihedrals = compute_backbone_dihedrals(MINI_PDB)
        assert isinstance(dihedrals, dict)
        assert "phi" in dihedrals
        assert "psi" in dihedrals

    def test_length_matches_residues(self) -> None:
        """phi and psi lists have length equal to number of residues."""
        dihedrals = compute_backbone_dihedrals(MINI_PDB)
        # 3 residues → 3 phi and 3 psi values
        assert len(dihedrals["phi"]) == 3
        assert len(dihedrals["psi"]) == 3

    def test_values_are_float_or_none(self) -> None:
        """All values are float or None."""
        dihedrals = compute_backbone_dihedrals(MINI_PDB)
        for val in dihedrals["phi"]:
            assert val is None or isinstance(val, float)
        for val in dihedrals["psi"]:
            assert val is None or isinstance(val, float)

    def test_empty_pdb_returns_empty(self) -> None:
        """Empty PDB returns empty phi/psi lists."""
        dihedrals = compute_backbone_dihedrals("")
        assert dihedrals["phi"] == []
        assert dihedrals["psi"] == []
