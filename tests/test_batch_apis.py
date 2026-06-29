"""Tests for batch APIs across all engines.

Verifies that every batch function returns BatchResult[T] with the correct
element type, and that BatchResult properties (results, successful, failed,
total, success_count, failure_count, total_time_s) work correctly.

Engines tested:
  - esmfold:  predict_structure_batch  -> BatchResult[ESMFoldResult]
  - camsol:   compute_solubility_batch -> BatchResult[CamSolResult]
  - immunogenicity: compute_immunogenicity_batch -> BatchResult[ImmunogenicityResult]
  - foldx:    run_stability_batch      -> BatchResult[FoldXResult]
"""

from __future__ import annotations

import sys
import time

import pytest

# ---------------------------------------------------------------------------
# Path setup — import from src/
# ---------------------------------------------------------------------------
sys.path.insert(0, "src")

from biocompiler.engines.base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
)

# Engine-specific imports
from biocompiler.engines.esmfold import (
    ESMFoldResult,
    predict_structure_batch,
)
from biocompiler.engines.camsol import (
    CamSolResult,
    compute_solubility_batch,
)
from biocompiler.immunogenicity.core import (
    ImmunogenicityResult,
    compute_immunogenicity_batch,
)
from biocompiler.engines.foldx import (
    FoldXResult,
    run_stability_batch,
)

# ---------------------------------------------------------------------------
# Test sequences — short valid protein sequences
# ---------------------------------------------------------------------------
VALID_SEQS = [
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAF",
    "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQD",
]

SHORT_SEQS = [
    "MVLSPADKTNVKAAWGKVGA",
    "MKWVTFISLLFLFSSAYSRG",
    "MADQLTEEQIAEFKEAFSLF",
]

# Invalid sequences that should produce failed results
INVALID_SEQS = [
    "12345",       # non-amino-acid characters
    "",            # empty string
    "   ",         # whitespace only
]


# =====================================================================
# 1. ESMFold batch API
# =====================================================================

class TestPredictStructureBatch:
    """Test predict_structure_batch returns BatchResult[ESMFoldResult]."""

    def test_returns_batch_result(self):
        """predict_structure_batch returns a BatchResult instance."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert isinstance(result, BatchResult)

    def test_results_are_esmfold_results(self):
        """Each item in results is an ESMFoldResult instance."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert isinstance(result.results, list)
        for r in result.results:
            assert isinstance(r, ESMFoldResult)

    def test_results_length_matches_input(self):
        """Number of results equals number of input sequences."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert len(result.results) == len(SHORT_SEQS)

    def test_successful_plus_failed_equals_total(self):
        """successful + failed == total."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert result.successful + result.failed == result.total

    def test_total_time_s_positive(self):
        """total_time_s should be >= 0 (may be 0 for very fast offline)."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert result.total_time_s >= 0

    def test_success_count_alias(self):
        """success_count is an alias for successful."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert result.success_count == result.successful

    def test_failure_count_alias(self):
        """failure_count is an alias for failed."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert result.failure_count == result.failed

    def test_results_is_list(self):
        """results attribute is a list."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        assert isinstance(result.results, list)

    def test_esmfold_result_inherits_base(self):
        """ESMFoldResult is instance of BaseEngineResult."""
        result = predict_structure_batch(SHORT_SEQS, max_concurrent=1)
        if result.results:
            assert isinstance(result.results[0], BaseEngineResult)


# =====================================================================
# 2. CamSol batch API
# =====================================================================

class TestComputeSolubilityBatch:
    """Test compute_solubility_batch returns BatchResult[CamSolResult]."""

    def test_returns_batch_result(self):
        """compute_solubility_batch returns a BatchResult instance."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert isinstance(result, BatchResult)

    def test_results_are_camsol_results(self):
        """Each item in results is a CamSolResult instance."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert isinstance(result.results, list)
        for r in result.results:
            assert isinstance(r, CamSolResult)

    def test_results_length_matches_input(self):
        """Number of results equals number of input sequences."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert len(result.results) == len(SHORT_SEQS)

    def test_all_successful(self):
        """Valid sequences should all produce successful results."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert result.successful == len(SHORT_SEQS)
        assert result.failed == 0

    def test_successful_plus_failed_equals_total(self):
        """successful + failed == total."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert result.successful + result.failed == result.total

    def test_total_time_s_positive(self):
        """total_time_s > 0 for real computation."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert result.total_time_s > 0

    def test_success_count_alias(self):
        """success_count is an alias for successful."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert result.success_count == result.successful

    def test_failure_count_alias(self):
        """failure_count is an alias for failed."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert result.failure_count == result.failed

    def test_results_is_list(self):
        """results attribute is a list."""
        result = compute_solubility_batch(SHORT_SEQS)
        assert isinstance(result.results, list)

    def test_camsol_result_inherits_base(self):
        """CamSolResult is instance of BaseEngineResult."""
        result = compute_solubility_batch(SHORT_SEQS)
        if result.results:
            assert isinstance(result.results[0], BaseEngineResult)

    def test_camsol_primary_score_label(self):
        """CamSolResult has primary_score_label = 'solubility'."""
        result = compute_solubility_batch(SHORT_SEQS)
        if result.results:
            assert result.results[0].primary_score_label == "solubility"

    def test_camsol_engine_name(self):
        """CamSolResult has engine_name = 'camsol'."""
        result = compute_solubility_batch(SHORT_SEQS)
        if result.results:
            assert result.results[0].engine_name == "camsol"

    def test_empty_input_returns_empty_batch(self):
        """Empty input list returns empty BatchResult."""
        result = compute_solubility_batch([])
        assert isinstance(result, BatchResult)
        assert result.total == 0
        assert result.results == []

    def test_invalid_sequence_produces_failed_result(self):
        """Invalid sequence should produce a failed CamSolResult in batch."""
        result = compute_solubility_batch(["123INVALID"])
        assert isinstance(result, BatchResult)
        assert len(result.results) == 1
        # The result should exist; may be successful=False
        r = result.results[0]
        assert isinstance(r, CamSolResult)
        assert r.success is False


# =====================================================================
# 3. Immunogenicity batch API
# =====================================================================

class TestComputeImmunogenicityBatch:
    """Test compute_immunogenicity_batch returns BatchResult[ImmunogenicityResult]."""

    def test_returns_batch_result(self):
        """compute_immunogenicity_batch returns a BatchResult instance."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert isinstance(result, BatchResult)

    def test_results_are_immunogenicity_results(self):
        """Each item in results is an ImmunogenicityResult instance."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert isinstance(result.results, list)
        for r in result.results:
            assert isinstance(r, ImmunogenicityResult)

    def test_results_length_matches_input(self):
        """Number of results equals number of input sequences."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert len(result.results) == len(SHORT_SEQS)

    def test_successful_plus_failed_equals_total(self):
        """successful + failed == total."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert result.successful + result.failed == result.total

    def test_total_time_s_positive(self):
        """total_time_s > 0 for real computation."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert result.total_time_s > 0

    def test_success_count_alias(self):
        """success_count is an alias for successful."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert result.success_count == result.successful

    def test_failure_count_alias(self):
        """failure_count is an alias for failed."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert result.failure_count == result.failed

    def test_results_is_list(self):
        """results attribute is a list."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        assert isinstance(result.results, list)

    def test_immunogenicity_result_inherits_base(self):
        """ImmunogenicityResult is instance of BaseEngineResult."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        if result.results:
            assert isinstance(result.results[0], BaseEngineResult)

    def test_immunogenicity_engine_name(self):
        """ImmunogenicityResult has engine_name = 'immunogenicity'."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        if result.results:
            assert result.results[0].engine_name == "immunogenicity"

    def test_immunogenicity_primary_score_label(self):
        """ImmunogenicityResult has primary_score_label = 'immunogenicity'."""
        result = compute_immunogenicity_batch(SHORT_SEQS)
        if result.results:
            assert result.results[0].primary_score_label == "immunogenicity"


# =====================================================================
# 4. FoldX batch API
# =====================================================================

class TestRunStabilityBatch:
    """Test run_stability_batch returns BatchResult[FoldXResult]."""

    def test_returns_batch_result(self):
        """run_stability_batch returns a BatchResult instance."""
        result = run_stability_batch(SHORT_SEQS)
        assert isinstance(result, BatchResult)

    def test_results_are_foldx_results(self):
        """Each item in results is a FoldXResult instance."""
        result = run_stability_batch(SHORT_SEQS)
        assert isinstance(result.results, list)
        for r in result.results:
            assert isinstance(r, FoldXResult)

    def test_results_length_matches_input(self):
        """Number of results equals number of input sequences."""
        result = run_stability_batch(SHORT_SEQS)
        assert len(result.results) == len(SHORT_SEQS)

    def test_all_successful(self):
        """Valid sequences should all produce successful results."""
        result = run_stability_batch(SHORT_SEQS)
        assert result.successful == len(SHORT_SEQS)
        assert result.failed == 0

    def test_successful_plus_failed_equals_total(self):
        """successful + failed == total."""
        result = run_stability_batch(SHORT_SEQS)
        assert result.successful + result.failed == result.total

    def test_total_time_s_non_negative(self):
        """total_time_s >= 0 for real computation (may be 0 for very fast empirical)."""
        result = run_stability_batch(SHORT_SEQS)
        assert result.total_time_s >= 0

    def test_success_count_alias(self):
        """success_count is an alias for successful."""
        result = run_stability_batch(SHORT_SEQS)
        assert result.success_count == result.successful

    def test_failure_count_alias(self):
        """failure_count is an alias for failed."""
        result = run_stability_batch(SHORT_SEQS)
        assert result.failure_count == result.failed

    def test_results_is_list(self):
        """results attribute is a list."""
        result = run_stability_batch(SHORT_SEQS)
        assert isinstance(result.results, list)

    def test_foldx_result_inherits_base(self):
        """FoldXResult is instance of BaseEngineResult."""
        result = run_stability_batch(SHORT_SEQS)
        if result.results:
            assert isinstance(result.results[0], BaseEngineResult)

    def test_foldx_engine_name(self):
        """FoldXResult has engine_name = 'foldx'."""
        result = run_stability_batch(SHORT_SEQS)
        if result.results:
            assert result.results[0].engine_name == "foldx"

    def test_foldx_primary_score_label(self):
        """FoldXResult has primary_score_label with delta-delta G."""
        result = run_stability_batch(SHORT_SEQS)
        if result.results:
            assert result.results[0].primary_score_label == "ΔΔG"


# =====================================================================
# 5. BatchResult property tests (using actual batch outputs)
# =====================================================================

class TestBatchResultProperties:
    """Test BatchResult properties using real engine batch outputs."""

    @pytest.fixture()
    def camsol_batch(self):
        """Return a real CamSol batch result for property testing."""
        return compute_solubility_batch(SHORT_SEQS)

    @pytest.fixture()
    def foldx_batch(self):
        """Return a real FoldX batch result for property testing."""
        return run_stability_batch(SHORT_SEQS)

    def test_results_is_list_camsol(self, camsol_batch):
        """results is a list for CamSol batch."""
        assert isinstance(camsol_batch.results, list)

    def test_successful_plus_failed_equals_total_camsol(self, camsol_batch):
        """successful + failed == total for CamSol batch."""
        assert camsol_batch.successful + camsol_batch.failed == camsol_batch.total

    def test_total_time_s_positive_camsol(self, camsol_batch):
        """total_time_s > 0 for CamSol batch."""
        assert camsol_batch.total_time_s > 0

    def test_success_count_alias_camsol(self, camsol_batch):
        """success_count == successful for CamSol batch."""
        assert camsol_batch.success_count == camsol_batch.successful

    def test_failure_count_alias_camsol(self, camsol_batch):
        """failure_count == failed for CamSol batch."""
        assert camsol_batch.failure_count == camsol_batch.failed

    def test_results_is_list_foldx(self, foldx_batch):
        """results is a list for FoldX batch."""
        assert isinstance(foldx_batch.results, list)

    def test_successful_plus_failed_equals_total_foldx(self, foldx_batch):
        """successful + failed == total for FoldX batch."""
        assert foldx_batch.successful + foldx_batch.failed == foldx_batch.total

    def test_total_time_s_non_negative_foldx(self, foldx_batch):
        """total_time_s >= 0 for FoldX batch (may round to 0 for very fast)."""
        assert foldx_batch.total_time_s >= 0

    def test_success_count_alias_foldx(self, foldx_batch):
        """success_count == successful for FoldX batch."""
        assert foldx_batch.success_count == foldx_batch.successful

    def test_failure_count_alias_foldx(self, foldx_batch):
        """failure_count == failed for FoldX batch."""
        assert foldx_batch.failure_count == foldx_batch.failed


# =====================================================================
# 6. BatchResult constructed directly (unit-level property tests)
# =====================================================================

class TestBatchResultDirectConstruction:
    """Test BatchResult properties via direct construction (no engine)."""

    def test_empty_batch_defaults(self):
        """Empty BatchResult has total=0, successful=0, failed=0."""
        br = BatchResult()
        assert br.results == []
        assert br.total == 0
        assert br.successful == 0
        assert br.failed == 0

    def test_auto_counting_from_results(self):
        """BatchResult auto-counts successful/failed from results."""
        r1 = CamSolResult(
            sequence="ACDE", primary_score=1.0, classification="soluble",
            success=True,
        )
        r2 = CamSolResult(
            sequence="FGHI", primary_score=0.0, classification="insoluble",
            success=False, error="test error",
        )
        br = BatchResult(results=[r1, r2])
        assert br.total == 2
        assert br.successful == 1
        assert br.failed == 1

    def test_success_count_alias(self):
        """success_count property aliases successful."""
        r1 = CamSolResult(
            sequence="ACDE", primary_score=1.0, classification="soluble",
            success=True,
        )
        br = BatchResult(results=[r1])
        assert br.success_count == br.successful == 1

    def test_failure_count_alias(self):
        """failure_count property aliases failed."""
        r1 = CamSolResult(
            sequence="ACDE", primary_score=0.0, classification="insoluble",
            success=False, error="fail",
        )
        br = BatchResult(results=[r1])
        assert br.failure_count == br.failed == 1

    def test_total_time_s_stored(self):
        """total_time_s is stored correctly."""
        br = BatchResult(total_time_s=1.234)
        assert br.total_time_s == 1.234

    def test_total_time_s_positive_requirement(self):
        """total_time_s should be > 0 for real computation."""
        # Simulate a real batch timing
        with EngineTimer() as timer:
            time.sleep(0.001)
        br = BatchResult(total_time_s=round(timer.elapsed, 4))
        assert br.total_time_s > 0

    def test_results_is_list(self):
        """results attribute is a list."""
        br = BatchResult(results=[])
        assert isinstance(br.results, list)

    def test_successful_plus_failed_equals_total_with_mixed(self):
        """successful + failed == total with mixed success/failure."""
        ok = CamSolResult(
            sequence="A", primary_score=1.0, classification="soluble", success=True,
        )
        fail = CamSolResult(
            sequence="B", primary_score=0.0, classification="insoluble",
            success=False, error="err",
        )
        br = BatchResult(results=[ok, fail, ok])
        assert br.successful + br.failed == br.total == 3


# =====================================================================
# 7. Cross-engine batch consistency
# =====================================================================

class TestCrossEngineBatchConsistency:
    """Verify all four engines return BatchResult with consistent properties."""

    @pytest.fixture(scope="class")
    def all_batch_results(self):
        """Run all four batch APIs and return their results."""
        seqs = SHORT_SEQS
        return {
            "esmfold": predict_structure_batch(seqs, max_concurrent=1),
            "camsol": compute_solubility_batch(seqs),
            "immunogenicity": compute_immunogenicity_batch(seqs),
            "foldx": run_stability_batch(seqs),
        }

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_returns_batch_result(self, all_batch_results, engine):
        """Each engine returns a BatchResult instance."""
        assert isinstance(all_batch_results[engine], BatchResult)

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_results_is_list(self, all_batch_results, engine):
        """Each engine's results is a list."""
        assert isinstance(all_batch_results[engine].results, list)

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_results_length_matches_input(self, all_batch_results, engine):
        """Each engine produces same number of results as inputs."""
        assert len(all_batch_results[engine].results) == len(SHORT_SEQS)

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_successful_plus_failed_equals_total(self, all_batch_results, engine):
        """successful + failed == total for every engine."""
        br = all_batch_results[engine]
        assert br.successful + br.failed == br.total

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_total_time_s_non_negative(self, all_batch_results, engine):
        """total_time_s >= 0 for every engine."""
        assert all_batch_results[engine].total_time_s >= 0

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_success_count_alias(self, all_batch_results, engine):
        """success_count aliases successful for every engine."""
        br = all_batch_results[engine]
        assert br.success_count == br.successful

    @pytest.mark.parametrize("engine", ["esmfold", "camsol", "immunogenicity", "foldx"])
    def test_failure_count_alias(self, all_batch_results, engine):
        """failure_count aliases failed for every engine."""
        br = all_batch_results[engine]
        assert br.failure_count == br.failed

    def test_esmfold_results_are_esmfold_result(self, all_batch_results):
        """ESMFold batch produces ESMFoldResult items."""
        for r in all_batch_results["esmfold"].results:
            assert isinstance(r, ESMFoldResult)

    def test_camsol_results_are_camsol_result(self, all_batch_results):
        """CamSol batch produces CamSolResult items."""
        for r in all_batch_results["camsol"].results:
            assert isinstance(r, CamSolResult)

    def test_immunogenicity_results_are_immunogenicity_result(self, all_batch_results):
        """Immunogenicity batch produces ImmunogenicityResult items."""
        for r in all_batch_results["immunogenicity"].results:
            assert isinstance(r, ImmunogenicityResult)

    def test_foldx_results_are_foldx_result(self, all_batch_results):
        """FoldX batch produces FoldXResult items."""
        for r in all_batch_results["foldx"].results:
            assert isinstance(r, FoldXResult)

    def test_all_results_inherit_base_engine_result(self, all_batch_results):
        """All results from all engines inherit from BaseEngineResult."""
        for engine_name, br in all_batch_results.items():
            for r in br.results:
                assert isinstance(r, BaseEngineResult), (
                    f"{engine_name} result {r} is not BaseEngineResult"
                )
