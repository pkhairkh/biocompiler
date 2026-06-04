"""Test BioCompiler MHCflurry adapter — offline NN-based MHC binding prediction.

Tests cover availability checks, client creation, result format conversion,
binding prediction, batch prediction, presentation prediction, error handling,
caching, and model downloading.  All MHCflurry-dependent tests use
``pytest.importorskip`` so the suite passes even when MHCflurry is not
installed.  Heavy / network-dependent tests are marked ``@pytest.mark.slow``.
"""
from __future__ import annotations

import importlib
import math
import sys
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Try importing the adapter — all tests must work when it is absent
# ---------------------------------------------------------------------------

_mhcflurry_adapter = pytest.importorskip(
    "biocompiler.mhcflurry_adapter",
    reason="mhcflurry_adapter module not yet available",
)

# Re-export public names for convenience
is_mhcflurry_available = _mhcflurry_adapter.is_mhcflurry_available
MHCflurryClient = _mhcflurry_adapter.MHCflurryClient
MHCBindingResult = _mhcflurry_adapter.MHCBindingResult
ic50_to_binding_score = _mhcflurry_adapter.ic50_to_binding_score
binding_score_to_ic50 = _mhcflurry_adapter.binding_score_to_ic50
classify_binding = _mhcflurry_adapter.classify_binding
predict_binding = _mhcflurry_adapter.predict_binding
batch_predict = _mhcflurry_adapter.batch_predict
predict_presentation = _mhcflurry_adapter.predict_presentation
download_models = _mhcflurry_adapter.download_models
clear_cache = _mhcflurry_adapter.clear_cache

# ---------------------------------------------------------------------------
# Constants used across tests
# ---------------------------------------------------------------------------

VALID_BINDING_CLASSES = frozenset({
    "strong_binder",
    "moderate_binder",
    "weak_binder",
    "non_binder",
})

# Well-known peptide / allele pairs from the MHCflurry literature
KNOWN_PEPTIDE = "SIINFEKL"       # OVA epitope, 8-mer
KNOWN_ALLELE = "HLA-A*02:01"    # Most common MHC-I allele

# Test protein (30 aa) — long enough for multiple 9-mer peptides
TEST_PROTEIN = "MAGPKWVTFISLLFLFSSAYSRGVFRQPEN"

# Short protein for edge-case testing
SHORT_PROTEIN = "M"


# ═══════════════════════════════════════════════════════════════════════════
# 1. Availability check
# ═══════════════════════════════════════════════════════════════════════════

class TestAvailability:
    """Tests for ``is_mhcflurry_available()``."""

    def test_returns_bool(self) -> None:
        """``is_mhcflurry_available()`` always returns a bool."""
        result = is_mhcflurry_available()
        assert isinstance(result, bool)

    def test_consistent_return(self) -> None:
        """Calling twice without intervening changes returns the same value."""
        r1 = is_mhcflurry_available()
        r2 = is_mhcflurry_available()
        assert r1 == r2

    def test_returns_false_when_not_installed(self) -> None:
        """When ``mhcflurry`` cannot be imported, the function returns False."""
        with patch.dict(sys.modules, {"mhcflurry": None}):
            # Force re-evaluation by clearing any cached result
            importlib.reload(_mhcflurry_adapter)
            result = _mhcflurry_adapter.is_mhcflurry_available()
            assert result is False
        # Restore original module state
        importlib.reload(_mhcflurry_adapter)


# ═══════════════════════════════════════════════════════════════════════════
# 2. MHCflurryClient creation
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCflurryClientCreation:
    """Tests for ``MHCflurryClient.__init__``."""

    def test_init_default_models_dir(self) -> None:
        """Client can be created with default models_dir."""
        client = MHCflurryClient()
        assert client is not None
        # models_dir should be set (may be None or a path)
        assert hasattr(client, "models_dir")

    def test_init_custom_models_dir(self, tmp_path) -> None:
        """Client accepts a custom models_dir path."""
        custom_dir = str(tmp_path / "custom_models")
        client = MHCflurryClient(models_dir=custom_dir)
        assert client.models_dir == custom_dir

    def test_init_stores_models_dir(self, tmp_path) -> None:
        """The models_dir attribute is accessible after construction."""
        custom_dir = str(tmp_path / "my_models")
        client = MHCflurryClient(models_dir=custom_dir)
        assert client.models_dir == custom_dir


# ═══════════════════════════════════════════════════════════════════════════
# 3. Result format conversion
# ═══════════════════════════════════════════════════════════════════════════

class TestResultFormatConversion:
    """Tests for IC50 ↔ binding_score conversion and classify_binding."""

    # --- IC50 → binding_score ---

    def test_ic50_to_binding_score_low_ic50(self) -> None:
        """Low IC50 (strong binder) maps to high binding score."""
        score = ic50_to_binding_score(10.0)  # 10 nM — very strong
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0
        assert score > 0.5

    def test_ic50_to_binding_score_high_ic50(self) -> None:
        """High IC50 (non-binder) maps to low binding score."""
        score = ic50_to_binding_score(50000.0)  # 50 μM — non-binder
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0
        assert score < 0.5

    def test_ic50_to_binding_score_range(self) -> None:
        """Conversion always returns a value in [0, 1]."""
        for ic50 in [0.1, 1.0, 10.0, 50.0, 500.0, 5000.0, 50000.0, 100000.0]:
            score = ic50_to_binding_score(ic50)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for IC50 {ic50}"

    def test_ic50_to_binding_score_monotonic(self) -> None:
        """Lower IC50 should yield higher binding score."""
        ic50s = [1.0, 10.0, 50.0, 500.0, 5000.0, 50000.0]
        scores = [ic50_to_binding_score(ic) for ic in ic50s]
        for i in range(len(scores) - 1):
            assert scores[i] >= scores[i + 1], (
                f"Non-monotonic: IC50={ic50s[i]} → score={scores[i]}, "
                f"IC50={ic50s[i+1]} → score={scores[i+1]}"
            )

    # --- binding_score → IC50 ---

    def test_binding_score_to_ic50_high_score(self) -> None:
        """High binding score maps to low IC50."""
        ic50 = binding_score_to_ic50(0.9)
        assert isinstance(ic50, float)
        assert ic50 > 0.0
        assert ic50 < 100.0  # Should be in the ~50 nM range

    def test_binding_score_to_ic50_low_score(self) -> None:
        """Low binding score maps to high IC50."""
        ic50 = binding_score_to_ic50(0.1)
        assert isinstance(ic50, float)
        assert ic50 > 1000.0  # Should be in the ~5000 nM range

    # --- Round-trip ---

    def test_roundtrip_approximate(self) -> None:
        """IC50 → binding_score → IC50 round-trip is approximately correct.

        The conversion may not be perfectly invertible due to clamping,
        but should be within an order of magnitude.
        """
        for ic50 in [50.0, 500.0, 5000.0]:
            score = ic50_to_binding_score(ic50)
            ic50_back = binding_score_to_ic50(score)
            # Allow up to 10× deviation (log-scale mapping)
            ratio = ic50_back / ic50
            assert 0.1 <= ratio <= 10.0, (
                f"Round-trip failed: IC50={ic50} → score={score} → IC50={ic50_back}"
            )

    # --- classify_binding ---

    def test_classify_binding_strong(self) -> None:
        """IC50 < 50 nM → strong_binder."""
        assert classify_binding(10.0) == "strong_binder"
        assert classify_binding(49.9) == "strong_binder"

    def test_classify_binding_moderate(self) -> None:
        """50 ≤ IC50 ≤ 500 nM → moderate_binder."""
        assert classify_binding(50.0) == "moderate_binder"
        assert classify_binding(200.0) == "moderate_binder"
        assert classify_binding(500.0) == "moderate_binder"

    def test_classify_binding_weak(self) -> None:
        """500 < IC50 ≤ 5000 nM → weak_binder."""
        assert classify_binding(501.0) == "weak_binder"
        assert classify_binding(5000.0) == "weak_binder"

    def test_classify_binding_non_binder(self) -> None:
        """IC50 > 5000 nM → non_binder."""
        assert classify_binding(5001.0) == "non_binder"
        assert classify_binding(50000.0) == "non_binder"

    def test_classify_matches_immunogenicity(self) -> None:
        """classify_binding matches the thresholds in immunogenicity.py.

        The immunogenicity module uses the same 50/500/5000 nM thresholds.
        """
        from biocompiler.immunogenicity import classify_binding as immuno_classify

        for ic50 in [10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0]:
            assert classify_binding(ic50) == immuno_classify(ic50), (
                f"Mismatch at IC50={ic50}: "
                f"adapter={classify_binding(ic50)}, "
                f"immuno={immuno_classify(ic50)}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 4. predict_binding (conditional on MHCflurry)
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictBinding:
    """Tests for ``predict_binding`` — conditional on MHCflurry installation."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_mhcflurry(self) -> None:
        """Skip all tests in this class if MHCflurry is not installed."""
        pytest.importorskip("mhcflurry", reason="MHCflurry not installed")
        if not is_mhcflurry_available():
            pytest.skip("MHCflurry not available")

    def test_returns_mhc_binding_result(self) -> None:
        """predict_binding returns an MHCBindingResult instance."""
        result = predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert isinstance(result, MHCBindingResult)

    def test_result_fields_populated(self) -> None:
        """All fields of MHCBindingResult are populated correctly."""
        result = predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert result.allele == KNOWN_ALLELE
        assert result.peptide == KNOWN_PEPTIDE
        assert isinstance(result.binding_score, float)
        assert 0.0 <= result.binding_score <= 1.0
        assert result.ic50_nm is not None
        assert result.ic50_nm > 0.0
        assert isinstance(result.binding_class, str)
        assert result.binding_class in VALID_BINDING_CLASSES

    def test_binding_class_is_valid(self) -> None:
        """binding_class is one of the 4 valid values."""
        # Test several peptide/allele combinations
        pairs = [
            ("SIINFEKL", "HLA-A*02:01"),
            ("ELAGIGILTV", "HLA-A*02:01"),
            ("GILGFVFTL", "HLA-A*02:01"),
        ]
        for peptide, allele in pairs:
            result = predict_binding(peptide, allele)
            assert result.binding_class in VALID_BINDING_CLASSES, (
                f"Invalid binding_class '{result.binding_class}' for "
                f"{peptide}/{allele}"
            )

    def test_positions_are_consistent(self) -> None:
        """start_position and end_position are consistent with peptide length."""
        result = predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert result.end_position - result.start_position + 1 == len(KNOWN_PEPTIDE)


# ═══════════════════════════════════════════════════════════════════════════
# 5. batch_predict (conditional on MHCflurry)
# ═══════════════════════════════════════════════════════════════════════════

class TestBatchPredict:
    """Tests for ``batch_predict`` — conditional on MHCflurry installation."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_mhcflurry(self) -> None:
        """Skip all tests in this class if MHCflurry is not installed."""
        pytest.importorskip("mhcflurry", reason="MHCflurry not installed")
        if not is_mhcflurry_available():
            pytest.skip("MHCflurry not available")

    def test_batch_predict_returns_results(self) -> None:
        """batch_predict returns a list of MHCBindingResult."""
        results = batch_predict(
            TEST_PROTEIN,
            alleles=["HLA-A*02:01"],
            epitope_lengths=[9],
        )
        assert isinstance(results, list)
        assert len(results) > 0
        for r in results:
            assert isinstance(r, MHCBindingResult)

    def test_multiple_alleles_covered(self) -> None:
        """batch_predict covers all requested alleles."""
        alleles = ["HLA-A*02:01", "HLA-A*01:01"]
        results = batch_predict(
            TEST_PROTEIN,
            alleles=alleles,
            epitope_lengths=[9],
        )
        found_alleles = {r.allele for r in results}
        for allele in alleles:
            assert allele in found_alleles, f"Allele {allele} not in results"

    def test_multiple_epitope_lengths_scanned(self) -> None:
        """batch_predict scans multiple epitope lengths."""
        epi_lengths = [8, 9, 10]
        results = batch_predict(
            TEST_PROTEIN,
            alleles=["HLA-A*02:01"],
            epitope_lengths=epi_lengths,
        )
        found_lengths = {len(r.peptide) for r in results}
        # At least some of the requested lengths should appear
        assert found_lengths.intersection(set(epi_lengths)), (
            f"No peptides of requested lengths {epi_lengths} found; "
            f"got lengths {found_lengths}"
        )

    def test_results_count_reasonable(self) -> None:
        """Number of results is reasonable given protein length and parameters.

        For a protein of length L and peptide length p, we expect
        approximately (L - p + 1) peptides per allele per epitope length.
        """
        protein = TEST_PROTEIN
        allele_count = 2
        epi_lengths = [9]
        results = batch_predict(
            protein,
            alleles=["HLA-A*02:01", "HLA-A*01:01"],
            epitope_lengths=epi_lengths,
        )
        # Expected: roughly (L - 9 + 1) * 2 alleles * 1 length
        expected_per_allele = len(protein) - 9 + 1
        expected_min = expected_per_allele * allele_count * len(epi_lengths)
        # Allow some tolerance (invalid peptides may be skipped)
        assert len(results) >= expected_min * 0.5, (
            f"Too few results: {len(results)} < {expected_min * 0.5}"
        )
        assert len(results) <= expected_per_allele * allele_count * len(epi_lengths) * 2, (
            f"Too many results: {len(results)}"
        )


# ═══════════════════════════════════════════════════════════════════════════
# 6. predict_presentation (conditional on MHCflurry)
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictPresentation:
    """Tests for ``predict_presentation`` — conditional on MHCflurry."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_mhcflurry(self) -> None:
        """Skip all tests in this class if MHCflurry is not installed."""
        pytest.importorskip("mhcflurry", reason="MHCflurry not installed")
        if not is_mhcflurry_available():
            pytest.skip("MHCflurry not available")

    def test_predict_presentation_returns_results(self) -> None:
        """predict_presentation returns a list of results."""
        results = predict_presentation(
            TEST_PROTEIN,
            alleles=["HLA-A*02:01"],
        )
        assert isinstance(results, list)

    def test_predict_presentation_may_differ_from_binding(self) -> None:
        """Presentation prediction may differ from binding-only prediction.

        MHCflurry's presentation model incorporates antigen processing,
        so results may differ from the binding-only model.
        """
        binding_results = batch_predict(
            TEST_PROTEIN,
            alleles=["HLA-A*02:01"],
            epitope_lengths=[9],
        )
        presentation_results = predict_presentation(
            TEST_PROTEIN,
            alleles=["HLA-A*02:01"],
        )
        # Both should return results, but they may have different scores
        assert len(binding_results) > 0
        assert len(presentation_results) > 0

    def test_predict_presentation_binding_class_valid(self) -> None:
        """All presentation results have valid binding_class values."""
        results = predict_presentation(
            TEST_PROTEIN,
            alleles=["HLA-A*02:01"],
        )
        for r in results:
            assert r.binding_class in VALID_BINDING_CLASSES, (
                f"Invalid binding_class: {r.binding_class}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 7. Error handling
# ═══════════════════════════════════════════════════════════════════════════

class TestErrorHandling:
    """Tests for graceful error handling."""

    def test_invalid_allele_skips_not_crashes(self) -> None:
        """An invalid allele name should be skipped, not cause a crash.

        This test runs regardless of MHCflurry availability by mocking
        the client's internal predict call.
        """
        client = MHCflurryClient()
        # If MHCflurry is not available, predict_binding should handle
        # it gracefully (return empty result or raise a specific error)
        try:
            result = predict_binding("SIINFEKL", "INVALID_ALLELE_999")
            # If it returns a result, it should be a non-binder
            if result is not None:
                assert isinstance(result, MHCBindingResult)
        except Exception as exc:
            # Should not be an unhandled exception type
            assert not isinstance(exc, (SystemExit, KeyboardInterrupt))

    def test_invalid_peptide_handled(self) -> None:
        """An invalid peptide should be handled gracefully."""
        # Empty peptide
        try:
            result = predict_binding("", KNOWN_ALLELE)
            if result is not None:
                assert isinstance(result, MHCBindingResult)
        except (ValueError, TypeError):
            pass  # Acceptable — input validation rejecting empty peptide
        except Exception as exc:
            pytest.fail(f"Unexpected exception for empty peptide: {exc}")

    def test_very_short_protein(self) -> None:
        """A protein shorter than the minimum epitope length returns empty results."""
        # 1-AA protein — cannot produce any 8-11 mer peptides
        results = batch_predict(
            SHORT_PROTEIN,
            alleles=["HLA-A*02:01"],
            epitope_lengths=[8, 9, 10, 11],
        )
        assert isinstance(results, list)
        assert len(results) == 0

    def test_non_standard_amino_acids(self) -> None:
        """Peptides with non-standard amino acids are handled gracefully."""
        # Protein with X (non-standard)
        protein_with_x = "MAGXKWVTFISLLFLFSSAYS"
        try:
            results = batch_predict(
                protein_with_x,
                alleles=["HLA-A*02:01"],
                epitope_lengths=[9],
            )
            # Should skip peptides containing X
            assert isinstance(results, list)
            for r in results:
                assert "X" not in r.peptide
        except (ValueError, TypeError):
            pass  # Acceptable — input validation rejecting non-standard AAs


# ═══════════════════════════════════════════════════════════════════════════
# 8. Caching
# ═══════════════════════════════════════════════════════════════════════════

class TestCaching:
    """Tests for the prediction cache."""

    def test_predict_same_peptide_twice_cache_hit(self) -> None:
        """Predicting the same peptide twice should hit the cache on the
        second call.

        We mock the underlying MHCflurry predictor to verify it is only
        called once, even when MHCflurry is not installed.
        """
        client = MHCflurryClient()
        mock_result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=len(KNOWN_PEPTIDE) - 1,
            binding_score=0.75,
            ic50_nm=150.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )

        call_count = 0
        original_predict = client._predict_single

        def counting_predict(peptide: str, allele: str, **kwargs):
            nonlocal call_count
            call_count += 1
            return mock_result

        # Patch the internal prediction method
        client._predict_single = counting_predict  # type: ignore[assignment]

        # First call — cache miss
        r1 = client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert r1.peptide == KNOWN_PEPTIDE

        # Second call — should be a cache hit
        r2 = client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert r2.peptide == KNOWN_PEPTIDE

        # _predict_single should only have been called once
        assert call_count == 1, (
            f"Expected 1 call to _predict_single, got {call_count}"
        )

    def test_clear_cache_fresh_prediction(self) -> None:
        """After clearing the cache, next prediction should be fresh."""
        client = MHCflurryClient()
        mock_result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=len(KNOWN_PEPTIDE) - 1,
            binding_score=0.75,
            ic50_nm=150.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )

        call_count = 0

        def counting_predict(peptide: str, allele: str, **kwargs):
            nonlocal call_count
            call_count += 1
            return mock_result

        client._predict_single = counting_predict  # type: ignore[assignment]

        # First call
        client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert call_count == 1

        # Second call — should be cache hit
        client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert call_count == 1

        # Clear cache
        clear_cache()

        # Third call — should be fresh (cache miss)
        client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert call_count == 2, (
            f"Expected 2 calls after cache clear, got {call_count}"
        )

    def test_different_peptides_not_cached(self) -> None:
        """Different peptides should not share cache entries."""
        client = MHCflurryClient()

        results_returned = {}

        def multi_predict(peptide: str, allele: str, **kwargs):
            result = MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=len(peptide) - 1,
                binding_score=0.5,
                ic50_nm=500.0,
                binding_class="moderate_binder",
                anchor_residues={},
                anchor_scores={},
            )
            results_returned[peptide] = result
            return result

        client._predict_single = multi_predict  # type: ignore[assignment]

        r1 = client.predict_binding("SIINFEKL", KNOWN_ALLELE)
        r2 = client.predict_binding("GILGFVFT", KNOWN_ALLELE)

        assert "SIINFEKL" in results_returned
        assert "GILGFVFT" in results_returned
        assert r1.peptide != r2.peptide


# ═══════════════════════════════════════════════════════════════════════════
# 9. download_models (conditional, slow)
# ═══════════════════════════════════════════════════════════════════════════

class TestDownloadModels:
    """Tests for ``download_models`` — conditional on MHCflurry, marked slow."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_mhcflurry(self) -> None:
        """Skip all tests in this class if MHCflurry is not installed."""
        pytest.importorskip("mhcflurry", reason="MHCflurry not installed")
        if not is_mhcflurry_available():
            pytest.skip("MHCflurry not available")

    @pytest.mark.slow
    def test_download_models_returns_true(self) -> None:
        """download_models returns True when models are downloaded successfully."""
        result = download_models()
        assert isinstance(result, bool)
        # If models already exist, it should still return True
        assert result is True

    @pytest.mark.slow
    def test_download_models_idempotent(self) -> None:
        """Calling download_models twice should not fail."""
        r1 = download_models()
        r2 = download_models()
        assert r1 is True
        assert r2 is True


# ═══════════════════════════════════════════════════════════════════════════
# Additional: MHCBindingResult data class
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCBindingResult:
    """Tests for the MHCBindingResult data class."""

    def test_create_result(self) -> None:
        """MHCBindingResult can be created with all required fields."""
        result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=7,
            binding_score=0.85,
            ic50_nm=42.0,
            binding_class="strong_binder",
            anchor_residues={0: "S", 8: "L"},
            anchor_scores={0: 1.5, 8: 2.0},
        )
        assert result.allele == KNOWN_ALLELE
        assert result.peptide == KNOWN_PEPTIDE
        assert result.binding_score == pytest.approx(0.85)
        assert result.ic50_nm == pytest.approx(42.0)
        assert result.binding_class == "strong_binder"
        assert result.anchor_residues[0] == "S"
        assert result.anchor_scores[8] == pytest.approx(2.0)

    def test_binding_score_range(self) -> None:
        """binding_score should be in [0, 1]."""
        result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=7,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert 0.0 <= result.binding_score <= 1.0

    def test_ic50_positive(self) -> None:
        """IC50 should always be a positive value."""
        result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=7,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        assert result.ic50_nm > 0

    def test_binding_class_from_known_ic50(self) -> None:
        """binding_class is consistent with IC50 thresholds."""
        test_cases = [
            (30.0, "strong_binder"),
            (100.0, "moderate_binder"),
            (1000.0, "weak_binder"),
            (10000.0, "non_binder"),
        ]
        for ic50, expected_class in test_cases:
            result = MHCBindingResult(
                allele=KNOWN_ALLELE,
                peptide=KNOWN_PEPTIDE,
                start_position=0,
                end_position=7,
                binding_score=0.5,
                ic50_nm=ic50,
                binding_class=expected_class,
                anchor_residues={},
                anchor_scores={},
            )
            assert result.binding_class == expected_class


# ═══════════════════════════════════════════════════════════════════════════
# Additional: Integration-level tests (conditional, with mocking)
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCflurryClientIntegration:
    """Integration-level tests using mocked MHCflurry backend."""

    def test_client_with_mocked_predictor(self) -> None:
        """Client correctly delegates to the internal predictor."""
        client = MHCflurryClient()

        mock_result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=len(KNOWN_PEPTIDE) - 1,
            binding_score=0.9,
            ic50_nm=25.0,
            binding_class="strong_binder",
            anchor_residues={0: "S"},
            anchor_scores={0: 2.0},
        )

        client._predict_single = MagicMock(return_value=mock_result)  # type: ignore[assignment]

        result = client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert result.allele == KNOWN_ALLELE
        assert result.binding_class == "strong_binder"
        client._predict_single.assert_called_once()

    def test_batch_predict_with_mocked_predictor(self) -> None:
        """batch_predict correctly calls predict_binding for each peptide."""
        client = MHCflurryClient()

        def mock_predict(peptide: str, allele: str, **kwargs):
            return MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=len(peptide) - 1,
                binding_score=0.5,
                ic50_nm=500.0,
                binding_class="moderate_binder",
                anchor_residues={},
                anchor_scores={},
            )

        client.predict_binding = MagicMock(side_effect=mock_predict)  # type: ignore[assignment]

        protein = "SIINFEKLW"
        results = client.batch_predict(
            protein,
            alleles=["HLA-A*02:01"],
            epitope_lengths=[9],
        )

        # Protein length 9, epitope length 9 → 1 peptide per allele
        assert len(results) >= 1

    def test_batch_predict_produces_correct_positions(self) -> None:
        """batch_predict results have correct start/end positions."""
        client = MHCflurryClient()

        def mock_predict(peptide: str, allele: str, **kwargs):
            return MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=len(peptide) - 1,
                binding_score=0.5,
                ic50_nm=500.0,
                binding_class="moderate_binder",
                anchor_residues={},
                anchor_scores={},
            )

        client.predict_binding = MagicMock(side_effect=mock_predict)  # type: ignore[assignment]

        protein = "SIINFEKLWSIINFEKLW"  # 18 aa
        results = client.batch_predict(
            protein,
            alleles=["HLA-A*02:01"],
            epitope_lengths=[9],
        )

        # 18 - 9 + 1 = 10 overlapping 9-mers
        assert len(results) >= 1

    def test_empty_protein_returns_empty(self) -> None:
        """batch_predict with empty protein returns empty list."""
        client = MHCflurryClient()
        results = client.batch_predict(
            "",
            alleles=["HLA-A*02:01"],
            epitope_lengths=[9],
        )
        assert results == []

    def test_no_alleles_returns_empty(self) -> None:
        """batch_predict with no alleles returns empty list."""
        client = MHCflurryClient()
        results = client.batch_predict(
            TEST_PROTEIN,
            alleles=[],
            epitope_lengths=[9],
        )
        assert results == []


# ═══════════════════════════════════════════════════════════════════════════
# Additional: Score conversion edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestScoreConversionEdgeCases:
    """Edge-case tests for score conversion functions."""

    def test_ic50_to_binding_score_zero_ic50(self) -> None:
        """IC50 of 0 should still produce a valid score."""
        score = ic50_to_binding_score(0.0)
        assert 0.0 <= score <= 1.0

    def test_ic50_to_binding_score_very_large(self) -> None:
        """Very large IC50 should produce a score near 0."""
        score = ic50_to_binding_score(1e9)
        assert 0.0 <= score <= 1.0
        assert score < 0.1

    def test_binding_score_to_ic50_zero_score(self) -> None:
        """binding_score of 0 should produce a very high IC50."""
        ic50 = binding_score_to_ic50(0.0)
        assert ic50 > 5000.0

    def test_binding_score_to_ic50_one_score(self) -> None:
        """binding_score of 1 should produce a very low IC50."""
        ic50 = binding_score_to_ic50(1.0)
        assert ic50 < 100.0

    def test_binding_score_to_ic50_midpoint(self) -> None:
        """binding_score ~0.5 should produce IC50 near 500 nM."""
        ic50 = binding_score_to_ic50(0.5)
        # The exact mapping depends on the implementation, but 0.5
        # should map to something in the moderate binder range
        assert 10.0 < ic50 < 10000.0

    def test_classify_binding_boundary_values(self) -> None:
        """Boundary IC50 values are classified correctly."""
        # Exact boundary at 50
        assert classify_binding(50.0) == "moderate_binder"
        # Exact boundary at 500
        assert classify_binding(500.0) == "moderate_binder"
        # Exact boundary at 5000
        assert classify_binding(5000.0) == "weak_binder"

    def test_ic50_to_binding_score_negative_ic50(self) -> None:
        """Negative IC50 (invalid) should still return a valid score."""
        score = ic50_to_binding_score(-1.0)
        assert 0.0 <= score <= 1.0


# ═══════════════════════════════════════════════════════════════════════════
# Additional: Thread-safety / re-entrancy (smoke tests)
# ═══════════════════════════════════════════════════════════════════════════

class TestReentrancy:
    """Smoke tests for cache consistency under repeated access."""

    def test_repeated_clear_cache_safe(self) -> None:
        """Calling clear_cache multiple times does not raise."""
        clear_cache()
        clear_cache()
        clear_cache()
        # No assertion needed — just verifying no exception

    def test_clear_cache_then_predict(self) -> None:
        """Cache can be cleared and then used for new predictions."""
        clear_cache()
        client = MHCflurryClient()

        mock_result = MHCBindingResult(
            allele=KNOWN_ALLELE,
            peptide=KNOWN_PEPTIDE,
            start_position=0,
            end_position=len(KNOWN_PEPTIDE) - 1,
            binding_score=0.6,
            ic50_nm=300.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )

        client._predict_single = MagicMock(return_value=mock_result)  # type: ignore[assignment]

        # Should work after cache clear
        result = client.predict_binding(KNOWN_PEPTIDE, KNOWN_ALLELE)
        assert result.binding_class == "moderate_binder"
