"""Integration tests for BioCompiler ViennaRNA module.

Covers:
1. ``is_viennarna_available`` returns a bool
2. ``compute_5prime_dg`` with a short sequence (skip if ViennaRNA unavailable)
3. ``check_mrna_structure_viennarna`` returns a typed result (skip if unavailable)
4. Fallback behaviour when ViennaRNA is not installed
5. Named constants (DEFAULT_5PRIME_WINDOW, etc.) exist and are reasonable

Uses ``pytest.importorskip`` / ``skipIf`` for ViennaRNA-dependent tests.
Always-importable module guarantee: importing ``biocompiler.viennarna`` must
never raise even when the ViennaRNA library is absent.

Code targets Python 3.10+.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

# ---------------------------------------------------------------------------
# Module import — guaranteed to succeed regardless of ViennaRNA installation
# ---------------------------------------------------------------------------

import biocompiler.viennarna as vr
from biocompiler.viennarna import (
    DEFAULT_5PRIME_WINDOW,
    DEFAULT_DG_THRESHOLD,
    DEFAULT_FULL_LENGTH_CUTOFF,
    DEFAULT_FOLD_TIMEOUT_SECONDS,
    DEFAULT_OVERLAP_THRESHOLD,
    DEFAULT_STEP,
    DEFAULT_WINDOW_SIZE,
    EXPECTED_VIENNARNA_VERSION,
    MAX_FOLD_TIMEOUT_SECONDS,
    NEAREST_NEIGHBOR_AU,
    NEAREST_NEIGHBOR_GC,
    NEAREST_NEIGHBOR_GU,
    REGION_5UTR,
    REGION_FULL,
    REGION_START_CODON,
    AccessibilityResult,
    MFEResult,
    MRNAStructureResult,
    StemLoop,
    check_mrna_structure_viennarna,
    compute_5prime_dg,
    is_viennarna_available,
    predict_mfe,
)

# ---------------------------------------------------------------------------
# Test sequences
# ---------------------------------------------------------------------------

# Short DNA sequence suitable for 5'-region folding
SEQ_SHORT_50 = "ATGGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

# GC-rich sequence — strongly structured
SEQ_GC_RICH = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"

# All-A sequence — essentially no structure
SEQ_ALL_A = "A" * 60

# Empty sequence
SEQ_EMPTY = ""


# ═══════════════════════════════════════════════════════════════════════════
# 1. is_viennarna_available returns bool
# ═══════════════════════════════════════════════════════════════════════════


class TestIsViennaRNAAvailable:
    """Verify that ``is_viennarna_available()`` returns a bool."""

    def test_returns_bool(self) -> None:
        """The function must always return a bool, never raise."""
        result = is_viennarna_available()
        assert isinstance(result, bool), (
            f"Expected bool, got {type(result).__name__}"
        )

    def test_consistent_on_repeated_calls(self) -> None:
        """Repeated calls return the same value (environment is stable)."""
        first = is_viennarna_available()
        second = is_viennarna_available()
        assert first == second, (
            f"Inconsistent availability: first={first}, second={second}"
        )

    def test_import_succeeds_regardless(self) -> None:
        """Importing the module must never raise, even without ViennaRNA."""
        import importlib

        mod = importlib.import_module("biocompiler.viennarna")
        assert hasattr(mod, "is_viennarna_available")


# ═══════════════════════════════════════════════════════════════════════════
# 2. compute_5prime_dg with short sequence
# ═══════════════════════════════════════════════════════════════════════════


class TestCompute5PrimeDG:
    """Tests for ``compute_5prime_dg()`` — skip if ViennaRNA unavailable."""

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_returns_float(self) -> None:
        """compute_5prime_dg returns a numeric value."""
        result = compute_5prime_dg(SEQ_SHORT_50)
        assert isinstance(result, (int, float)), (
            f"Expected numeric, got {type(result).__name__}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_gc_rich_negative_dg(self) -> None:
        """GC-rich 5' region should have negative ΔG (stable structure)."""
        result = compute_5prime_dg(SEQ_GC_RICH)
        assert result < 0, (
            f"GC-rich 5' region should have negative ΔG, got {result}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_all_a_near_zero_dg(self) -> None:
        """All-A 5' region should have ΔG close to zero (no structure)."""
        result = compute_5prime_dg(SEQ_ALL_A)
        assert result > -5.0, (
            f"All-A 5' region shouldn't have very negative ΔG, got {result}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_custom_window_parameter(self) -> None:
        """Passing a custom window size works without error."""
        result = compute_5prime_dg(SEQ_SHORT_50, window=30)
        assert isinstance(result, (int, float))

    def test_empty_sequence_returns_zero(self) -> None:
        """Empty sequence returns 0.0 regardless of ViennaRNA availability."""
        result = compute_5prime_dg(SEQ_EMPTY)
        assert result == 0.0

    def test_very_short_sequence_returns_zero(self) -> None:
        """Sequence shorter than 4 nt returns 0.0 (trivial fold)."""
        result = compute_5prime_dg("ATG")
        assert result == 0.0


# ═══════════════════════════════════════════════════════════════════════════
# 3. check_mrna_structure_viennarna returns typed result
# ═══════════════════════════════════════════════════════════════════════════


class TestCheckMRNAStructureViennaRNA:
    """Tests for ``check_mrna_structure_viennarna()`` — skip if unavailable."""

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_returns_dict_with_required_keys(self) -> None:
        """Result contains all required keys of MRNAStructureResult."""
        result = check_mrna_structure_viennarna(SEQ_SHORT_50)
        assert isinstance(result, dict)
        for key in ("dg", "method", "structure", "viennarna_used"):
            assert key in result, f"Missing key: {key!r}"

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_dg_is_numeric(self) -> None:
        """dg field is a number."""
        result = check_mrna_structure_viennarna(SEQ_SHORT_50)
        assert isinstance(result["dg"], (int, float)), (
            f"dg should be numeric, got {type(result['dg']).__name__}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_method_is_string(self) -> None:
        """method field is a string."""
        result = check_mrna_structure_viennarna(SEQ_SHORT_50)
        assert isinstance(result["method"], str)

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_structure_is_string(self) -> None:
        """structure field is a string (dot-bracket or empty)."""
        result = check_mrna_structure_viennarna(SEQ_SHORT_50)
        assert isinstance(result["structure"], str)

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_viennarna_used_is_bool(self) -> None:
        """viennarna_used field is a bool."""
        result = check_mrna_structure_viennarna(SEQ_SHORT_50)
        assert isinstance(result["viennarna_used"], bool)

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_method_is_viennarna_when_available(self) -> None:
        """When ViennaRNA is available, method should be a ViennaRNA backend."""
        result = check_mrna_structure_viennarna(SEQ_SHORT_50)
        assert result["viennarna_used"] is True, (
            "ViennaRNA is available but viennarna_used is False"
        )
        assert result["method"] in ("viennarna_python", "viennarna_cli"), (
            f"Expected ViennaRNA method, got {result['method']!r}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_gc_rich_negative_dg(self) -> None:
        """GC-rich sequence should produce a negative ΔG."""
        result = check_mrna_structure_viennarna(SEQ_GC_RICH)
        assert result["dg"] < 0, (
            f"GC-rich should have negative ΔG, got {result['dg']}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_structure_length_matches_window(self) -> None:
        """Structure length should equal the window size."""
        window_end = 50
        result = check_mrna_structure_viennarna(
            SEQ_SHORT_50, window_start=0, window_end=window_end,
        )
        assert len(result["structure"]) == window_end, (
            f"Structure length {len(result['structure'])} != window {window_end}"
        )

    @pytest.mark.skipif(
        not is_viennarna_available(),
        reason="ViennaRNA not available",
    )
    def test_custom_window_parameters(self) -> None:
        """Custom window_start and window_end are respected."""
        result = check_mrna_structure_viennarna(
            SEQ_SHORT_50, window_start=0, window_end=30,
        )
        assert isinstance(result["dg"], (int, float))
        assert len(result["structure"]) == 30

    def test_very_short_sequence_trivial(self) -> None:
        """Sequences shorter than 4 nt use the trivial path (always works)."""
        result = check_mrna_structure_viennarna("ATG")
        assert isinstance(result, dict)
        assert result["method"] == "trivial"
        assert result["viennarna_used"] is False
        assert result["dg"] == 0.0

    def test_mrna_structure_result_typed_dict_exists(self) -> None:
        """The MRNAStructureResult TypedDict is importable and has annotations."""
        assert hasattr(MRNAStructureResult, "__annotations__")
        for key in ("dg", "method", "structure", "viennarna_used"):
            assert key in MRNAStructureResult.__annotations__, (
                f"MRNAStructureResult missing annotation for {key!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 4. Fallback behaviour when ViennaRNA is not installed
# ═══════════════════════════════════════════════════════════════════════════


class TestFallbackBehaviour:
    """Verify graceful degradation when ViennaRNA is unavailable.

    Uses ``unittest.mock.patch`` to simulate a missing ViennaRNA installation
    so these tests always run.
    """

    def test_predict_mfe_returns_unavailable_result(self) -> None:
        """predict_mfe returns MFEResult(success=False) when ViennaRNA missing."""
        with patch.object(vr, "_fold_mfe") as mock_fold:
            mock_fold.return_value = MFEResult(
                success=False,
                method="unavailable",
                error="Simulated: ViennaRNA not available",
            )
            result = predict_mfe(SEQ_SHORT_50)
            assert isinstance(result, MFEResult)
            assert result.success is False
            assert result.method == "unavailable"

    def test_check_mrna_falls_back_to_toy_model(self) -> None:
        """check_mrna_structure_viennarna falls back to toy_hairpin_fallback."""
        with patch.object(vr, "_fold_mfe") as mock_fold:
            mock_fold.return_value = MFEResult(
                success=False,
                method="unavailable",
                error="Simulated: ViennaRNA not available",
            )
            result = check_mrna_structure_viennarna(SEQ_GC_RICH)
            assert isinstance(result, dict)
            assert result["method"] == "toy_hairpin_fallback"
            assert result["viennarna_used"] is False
            # Toy model should still compute a ΔG
            assert isinstance(result["dg"], (int, float))

    def test_check_mrna_toy_model_gc_negative(self) -> None:
        """Toy model fallback produces negative ΔG for GC-rich sequences."""
        with patch.object(vr, "_fold_mfe") as mock_fold:
            mock_fold.return_value = MFEResult(
                success=False,
                method="unavailable",
                error="Simulated: ViennaRNA not available",
            )
            result = check_mrna_structure_viennarna(SEQ_GC_RICH)
            # GC pairs contribute -1.5 each via NEAREST_NEIGHBOR_GC
            assert result["dg"] < 0, (
                f"Toy model should detect GC pairs, got ΔG={result['dg']}"
            )

    def test_check_mrna_toy_model_all_a_near_zero(self) -> None:
        """Toy model fallback gives near-zero ΔG for all-A sequences."""
        with patch.object(vr, "_fold_mfe") as mock_fold:
            mock_fold.return_value = MFEResult(
                success=False,
                method="unavailable",
                error="Simulated: ViennaRNA not available",
            )
            result = check_mrna_structure_viennarna(SEQ_ALL_A)
            # All-A has no complementary pairs; toy model should give ~0
            assert result["dg"] == pytest.approx(0.0, abs=0.01), (
                f"All-A toy model should give ~0 ΔG, got {result['dg']}"
            )

    def test_compute_5prime_dg_returns_zero_on_failure(self) -> None:
        """compute_5prime_dg returns 0.0 when folding fails."""
        with patch.object(vr, "_fold_mfe") as mock_fold:
            mock_fold.return_value = MFEResult(
                success=False,
                method="unavailable",
                error="Simulated: ViennaRNA not available",
            )
            result = compute_5prime_dg(SEQ_SHORT_50)
            assert result == 0.0

    def test_predict_mfe_empty_sequence(self) -> None:
        """predict_mfe with empty sequence returns failure result."""
        result = predict_mfe(SEQ_EMPTY)
        assert isinstance(result, MFEResult)
        assert result.success is False

    def test_predict_mfe_trivial_short_sequence(self) -> None:
        """predict_mfe with <4 nt sequence uses trivial path (success)."""
        # len < 4 triggers the trivial path; "ATCG" has len=4 which does NOT
        # trigger it (condition is strict <).  Use a 3-nt sequence instead.
        result = predict_mfe("ATG")
        assert isinstance(result, MFEResult)
        assert result.success is True
        assert result.method == "trivial"

    def test_unavailable_mfe_result_fields(self) -> None:
        """MFEResult from unavailable backend has expected field values."""
        with patch.object(vr, "_fold_mfe") as mock_fold:
            mock_fold.return_value = MFEResult(
                success=False,
                method="unavailable",
                error="Test error",
                sequence="AUGCGAU",
            )
            result = predict_mfe("ATGCGAT")  # longer than 4 to reach _fold_mfe
            # Note: predict_mfe extracts region & converts to RNA first,
            # then calls _fold_mfe with the RNA. The returned MFEResult
            # is passed through, so we check the mocked fields.
            assert result.success is False
            assert result.method == "unavailable"
            assert result.error is not None


# ═══════════════════════════════════════════════════════════════════════════
# 5. Named constants exist and are reasonable
# ═══════════════════════════════════════════════════════════════════════════


class TestNamedConstants:
    """Verify all named constants exist and have reasonable values."""

    # -- Window / threshold constants --

    def test_default_5prime_window(self) -> None:
        """DEFAULT_5PRIME_WINDOW should be a positive integer (biologically ~50 nt)."""
        assert isinstance(DEFAULT_5PRIME_WINDOW, int)
        assert DEFAULT_5PRIME_WINDOW > 0
        assert 20 <= DEFAULT_5PRIME_WINDOW <= 200, (
            f"DEFAULT_5PRIME_WINDOW={DEFAULT_5PRIME_WINDOW} seems unreasonable"
        )

    def test_default_dg_threshold(self) -> None:
        """DEFAULT_DG_THRESHOLD should be negative (stable structures have ΔG < 0)."""
        assert isinstance(DEFAULT_DG_THRESHOLD, float)
        assert DEFAULT_DG_THRESHOLD < 0, (
            f"DEFAULT_DG_THRESHOLD should be negative, got {DEFAULT_DG_THRESHOLD}"
        )

    def test_default_window_size(self) -> None:
        """DEFAULT_WINDOW_SIZE should be a positive integer (sliding window)."""
        assert isinstance(DEFAULT_WINDOW_SIZE, int)
        assert DEFAULT_WINDOW_SIZE > 0

    def test_default_step(self) -> None:
        """DEFAULT_STEP should be a positive integer ≤ DEFAULT_WINDOW_SIZE."""
        assert isinstance(DEFAULT_STEP, int)
        assert DEFAULT_STEP > 0
        assert DEFAULT_STEP <= DEFAULT_WINDOW_SIZE, (
            f"DEFAULT_STEP ({DEFAULT_STEP}) should be ≤ DEFAULT_WINDOW_SIZE ({DEFAULT_WINDOW_SIZE})"
        )

    def test_default_full_length_cutoff(self) -> None:
        """DEFAULT_FULL_LENGTH_CUTOFF should be a positive integer."""
        assert isinstance(DEFAULT_FULL_LENGTH_CUTOFF, int)
        assert DEFAULT_FULL_LENGTH_CUTOFF > 0

    # -- Timeout constants --

    def test_default_fold_timeout_seconds(self) -> None:
        """DEFAULT_FOLD_TIMEOUT_SECONDS should be a positive int (reasonable wait)."""
        assert isinstance(DEFAULT_FOLD_TIMEOUT_SECONDS, int)
        assert DEFAULT_FOLD_TIMEOUT_SECONDS > 0
        assert DEFAULT_FOLD_TIMEOUT_SECONDS <= 60, (
            f"DEFAULT_FOLD_TIMEOUT_SECONDS={DEFAULT_FOLD_TIMEOUT_SECONDS} too large"
        )

    def test_max_fold_timeout_seconds(self) -> None:
        """MAX_FOLD_TIMEOUT_SECONDS should be ≥ DEFAULT_FOLD_TIMEOUT_SECONDS."""
        assert isinstance(MAX_FOLD_TIMEOUT_SECONDS, int)
        assert MAX_FOLD_TIMEOUT_SECONDS >= DEFAULT_FOLD_TIMEOUT_SECONDS, (
            f"MAX ({MAX_FOLD_TIMEOUT_SECONDS}) < DEFAULT ({DEFAULT_FOLD_TIMEOUT_SECONDS})"
        )

    # -- Overlap threshold --

    def test_default_overlap_threshold(self) -> None:
        """DEFAULT_OVERLAP_THRESHOLD should be in (0, 1]."""
        assert isinstance(DEFAULT_OVERLAP_THRESHOLD, float)
        assert 0.0 < DEFAULT_OVERLAP_THRESHOLD <= 1.0

    # -- Nearest-neighbor energies --

    def test_nearest_neighbor_gc(self) -> None:
        """NEAREST_NEIGHBOR_GC should be negative (GC pairs are stabilizing)."""
        assert isinstance(NEAREST_NEIGHBOR_GC, float)
        assert NEAREST_NEIGHBOR_GC < 0

    def test_nearest_neighbor_au(self) -> None:
        """NEAREST_NEIGHBOR_AU should be negative but weaker than GC."""
        assert isinstance(NEAREST_NEIGHBOR_AU, float)
        assert NEAREST_NEIGHBOR_AU < 0
        assert NEAREST_NEIGHBOR_AU > NEAREST_NEIGHBOR_GC, (
            "AU pairs should be weaker (less negative) than GC pairs"
        )

    def test_nearest_neighbor_gu(self) -> None:
        """NEAREST_NEIGHBOR_GU should be negative but weaker than GC and AU."""
        assert isinstance(NEAREST_NEIGHBOR_GU, float)
        assert NEAREST_NEIGHBOR_GU < 0
        assert NEAREST_NEIGHBOR_GU > NEAREST_NEIGHBOR_GC, (
            "GU wobble should be weaker than GC"
        )

    def test_gc_strongest_pair(self) -> None:
        """GC pair energy should be the most negative (strongest)."""
        assert NEAREST_NEIGHBOR_GC <= NEAREST_NEIGHBOR_AU
        assert NEAREST_NEIGHBOR_GC <= NEAREST_NEIGHBOR_GU

    # -- Region labels --

    def test_region_labels_are_strings(self) -> None:
        """Region label constants are non-empty strings."""
        for label in (REGION_FULL, REGION_5UTR, REGION_START_CODON):
            assert isinstance(label, str)
            assert len(label) > 0

    def test_region_labels_distinct(self) -> None:
        """Region labels should be distinct."""
        labels = {REGION_FULL, REGION_5UTR, REGION_START_CODON}
        assert len(labels) == 3, "Region labels must be unique"

    # -- Version constant --

    def test_expected_viennarna_version_is_tuple(self) -> None:
        """EXPECTED_VIENNARNA_VERSION should be a 3-element tuple of ints."""
        assert isinstance(EXPECTED_VIENNARNA_VERSION, tuple)
        assert len(EXPECTED_VIENNARNA_VERSION) == 3
        for component in EXPECTED_VIENNARNA_VERSION:
            assert isinstance(component, int)
            assert component >= 0

    # -- Consistency checks --

    def test_5prime_window_lte_full_cutoff(self) -> None:
        """5' window should be ≤ full-length cutoff (5'UTR is a sub-region)."""
        assert DEFAULT_5PRIME_WINDOW <= DEFAULT_FULL_LENGTH_CUTOFF, (
            f"5PRIME_WINDOW ({DEFAULT_5PRIME_WINDOW}) > FULL_CUTOFF ({DEFAULT_FULL_LENGTH_CUTOFF})"
        )
