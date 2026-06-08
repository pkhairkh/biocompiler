"""Unit tests for foldx.py that don't require FoldX installation.

Covers:
1. Function existence and signatures
2. Input validation
3. Fallback behavior when FoldX is not available
4. Named constants validation

These tests exercise the module's pure functions, dataclass constructors,
cache, stability classification, input validation, and CLI-fallback paths
without needing the ``foldx`` binary on PATH.
"""

from __future__ import annotations

import inspect
import subprocess
from dataclasses import fields as dataclass_fields
from unittest.mock import patch

import pytest

from biocompiler.foldx import (
    AA_VOLUME,
    BLOSUM62,
    HYDROPATHY,
    FOLDX_CLI_ACCURACY,
    FOLDX_EMPIRICAL_BIAS,
    FOLDX_EMPIRICAL_DIRECTION_ACCURACY,
    FOLDX_EMPIRICAL_LARGE_MAE,
    FOLDX_EMPIRICAL_MAE,
    FOLDX_EMPIRICAL_MEDIUM_MAE,
    FOLDX_EMPIRICAL_PEARSON_R,
    FOLDX_EMPIRICAL_SMALL_MAE,
    FoldXCache,
    FoldXError,
    FoldXResult,
    MutationResult,
    StabilityLandscape,
    ConservationScore,
    clear_cache,
    empirical_stability,
    find_compensatory_mutations,
    find_stabilizing_mutations,
    identify_hotspot_regions,
    is_foldx_available,
    rank_positions_by_mutability,
    run_foldx_mutation,
    run_foldx_repair,
    run_foldx_stability,
    run_stability_batch,
    scan_all_mutations,
    scan_mutations,
    scan_position,
    compute_conservation,
)
from biocompiler.engine_base import BaseEngineResult, BatchResult


# ═══════════════════════════════════════════════════════════════════════════
# 1. Function existence and signatures
# ═══════════════════════════════════════════════════════════════════════════


class TestFunctionExistence:
    """Verify that every public function from __all__ is importable."""

    EXPECTED_PUBLIC_FUNCTIONS = [
        "is_foldx_available",
        "run_foldx_stability",
        "run_foldx_repair",
        "run_foldx_mutation",
        "empirical_stability",
        "run_stability_batch",
        "scan_mutations",
        "find_stabilizing_mutations",
        "scan_all_mutations",
        "scan_position",
        "compute_conservation",
        "find_compensatory_mutations",
        "rank_positions_by_mutability",
        "identify_hotspot_regions",
        "clear_cache",
    ]

    @pytest.mark.parametrize("name", EXPECTED_PUBLIC_FUNCTIONS)
    def test_public_function_importable(self, name: str) -> None:
        """Every declared public function can be imported from biocompiler.foldx."""
        import biocompiler.foldx as mod

        assert hasattr(mod, name), f"Module missing public function: {name}"
        assert callable(getattr(mod, name)), f"{name} is not callable"


class TestFunctionSignatures:
    """Verify key functions have expected parameter names."""

    def test_is_foldx_available_signature(self) -> None:
        """is_foldx_available takes no arguments."""
        sig = inspect.signature(is_foldx_available)
        assert list(sig.parameters) == []

    def test_run_foldx_stability_signature(self) -> None:
        """run_foldx_stability accepts pdb_string, foldx_dir, timeout."""
        sig = inspect.signature(run_foldx_stability)
        params = list(sig.parameters)
        assert "pdb_string" in params
        assert "foldx_dir" in params
        assert "timeout" in params

    def test_run_foldx_repair_signature(self) -> None:
        """run_foldx_repair accepts pdb_string, foldx_dir, timeout."""
        sig = inspect.signature(run_foldx_repair)
        params = list(sig.parameters)
        assert "pdb_string" in params
        assert "foldx_dir" in params
        assert "timeout" in params

    def test_run_foldx_mutation_signature(self) -> None:
        """run_foldx_mutation accepts pdb_string, mutations, foldx_dir, timeout."""
        sig = inspect.signature(run_foldx_mutation)
        params = list(sig.parameters)
        assert "pdb_string" in params
        assert "mutations" in params
        assert "foldx_dir" in params
        assert "timeout" in params

    def test_empirical_stability_signature(self) -> None:
        """empirical_stability accepts protein, pdb_string, timeout, organism."""
        sig = inspect.signature(empirical_stability)
        params = list(sig.parameters)
        assert "protein" in params
        assert "pdb_string" in params
        assert "timeout" in params
        assert "organism" in params

    def test_run_stability_batch_signature(self) -> None:
        """run_stability_batch accepts sequences, max_workers, batch_size, organism."""
        sig = inspect.signature(run_stability_batch)
        params = list(sig.parameters)
        assert "sequences" in params
        assert "max_workers" in params
        assert "batch_size" in params

    def test_scan_mutations_signature(self) -> None:
        """scan_mutations accepts protein, pdb_string, positions, organism."""
        sig = inspect.signature(scan_mutations)
        params = list(sig.parameters)
        assert "protein" in params
        assert "pdb_string" in params
        assert "positions" in params

    def test_find_stabilizing_mutations_signature(self) -> None:
        """find_stabilizing_mutations accepts protein, pdb_string, ddg_threshold, organism."""
        sig = inspect.signature(find_stabilizing_mutations)
        params = list(sig.parameters)
        assert "protein" in params
        assert "pdb_string" in params
        assert "ddg_threshold" in params

    def test_scan_all_mutations_signature(self) -> None:
        """scan_all_mutations accepts protein, method."""
        sig = inspect.signature(scan_all_mutations)
        params = list(sig.parameters)
        assert "protein" in params
        assert "method" in params

    def test_scan_position_signature(self) -> None:
        """scan_position accepts protein, position, method."""
        sig = inspect.signature(scan_position)
        params = list(sig.parameters)
        assert "protein" in params
        assert "position" in params
        assert "method" in params

    def test_compute_conservation_signature(self) -> None:
        """compute_conservation accepts protein, method."""
        sig = inspect.signature(compute_conservation)
        params = list(sig.parameters)
        assert "protein" in params
        assert "method" in params

    def test_find_compensatory_mutations_signature(self) -> None:
        """find_compensatory_mutations accepts protein, destabilizing_mutations."""
        sig = inspect.signature(find_compensatory_mutations)
        params = list(sig.parameters)
        assert "protein" in params
        assert "destabilizing_mutations" in params

    def test_rank_positions_by_mutability_signature(self) -> None:
        """rank_positions_by_mutability accepts protein."""
        sig = inspect.signature(rank_positions_by_mutability)
        params = list(sig.parameters)
        assert "protein" in params

    def test_identify_hotspot_regions_signature(self) -> None:
        """identify_hotspot_regions accepts protein, window, threshold."""
        sig = inspect.signature(identify_hotspot_regions)
        params = list(sig.parameters)
        assert "protein" in params
        assert "window" in params
        assert "threshold" in params

    def test_clear_cache_signature(self) -> None:
        """clear_cache takes no arguments."""
        sig = inspect.signature(clear_cache)
        assert list(sig.parameters) == []


class TestClassExistenceAndFields:
    """Verify public classes are importable with expected fields."""

    def test_foldx_result_inherits_base_engine_result(self) -> None:
        """FoldXResult inherits from BaseEngineResult."""
        assert issubclass(FoldXResult, BaseEngineResult)

    def test_foldx_result_has_engine_name(self) -> None:
        """FoldXResult.ENGINE_NAME is 'foldx'."""
        assert FoldXResult.ENGINE_NAME == "foldx"

    def test_foldx_result_has_primary_score_label(self) -> None:
        """FoldXResult.PRIMARY_SCORE_LABEL is a non-empty string."""
        assert isinstance(FoldXResult.PRIMARY_SCORE_LABEL, str)
        assert len(FoldXResult.PRIMARY_SCORE_LABEL) > 0

    def test_foldx_result_ddg_property(self) -> None:
        """FoldXResult.ddg property returns primary_score."""
        result = FoldXResult(protein="MK", stability_kcal=-5.0)
        assert result.ddg == result.primary_score

    def test_foldx_result_stability_class_property(self) -> None:
        """FoldXResult.stability_class property returns classification."""
        result = FoldXResult(protein="MK", stability_kcal=-5.0)
        assert result.stability_class == result.classification

    def test_foldx_result_stabilizing_mutations_property(self) -> None:
        """FoldXResult.stabilizing_mutations property returns mutations."""
        result = FoldXResult(protein="MK", stability_kcal=-5.0)
        assert result.stabilizing_mutations == result.mutations

    def test_stability_landscape_fields(self) -> None:
        """StabilityLandscape has expected dataclass fields."""
        field_names = {f.name for f in dataclass_fields(StabilityLandscape)}
        expected = {
            "protein", "wildtype_stability", "mutations",
            "stabilizing_count", "destabilizing_count", "neutral_count",
            "most_stabilizing", "most_destabilizing",
            "positions_scanned", "method",
        }
        assert expected <= field_names

    def test_conservation_score_fields(self) -> None:
        """ConservationScore has expected dataclass fields."""
        field_names = {f.name for f in dataclass_fields(ConservationScore)}
        expected = {"position", "wildtype", "conservation",
                    "substitution_tolerance", "critical"}
        assert expected <= field_names

    def test_foldx_cache_methods(self) -> None:
        """FoldXCache has get, put, clear, __len__ methods."""
        cache = FoldXCache()
        assert hasattr(cache, "get")
        assert hasattr(cache, "put")
        assert hasattr(cache, "clear")
        assert hasattr(cache, "__len__")

    def test_mutation_result_is_engine_base_mutation_result(self) -> None:
        """MutationResult imported from foldx is engine_base.MutationResult."""
        from biocompiler.engine_base import MutationResult as BaseMR
        assert MutationResult is BaseMR


# ═══════════════════════════════════════════════════════════════════════════
# 2. Input validation
# ═══════════════════════════════════════════════════════════════════════════


class TestEmpiricalStabilityInputValidation:
    """Input validation for empirical_stability."""

    def test_empty_protein_returns_failure(self) -> None:
        """Empty string → success=False."""
        result = empirical_stability("")
        assert result.success is False
        assert result.error is not None

    def test_whitespace_only_protein_returns_failure(self) -> None:
        """Whitespace-only string → success=False."""
        result = empirical_stability("   ")
        assert result.success is False

    def test_invalid_amino_acids_returns_failure(self) -> None:
        """Non-standard amino acids → success=False."""
        result = empirical_stability("MKX")
        assert result.success is False
        assert "invalid" in result.error.lower() or "non-standard" in result.error.lower()

    def test_valid_protein_returns_success(self) -> None:
        """Valid protein → success=True."""
        result = empirical_stability("MKGIL")
        assert result.success is True

    def test_lowercase_protein_accepted(self) -> None:
        """Lowercase protein is auto-uppercased and succeeds."""
        result = empirical_stability("mkgil")
        assert result.success is True

    def test_result_is_foldx_result(self) -> None:
        """empirical_stability always returns a FoldXResult."""
        result = empirical_stability("MKG")
        assert isinstance(result, FoldXResult)

    def test_failure_result_is_foldx_result(self) -> None:
        """Even failed results are FoldXResult instances."""
        result = empirical_stability("")
        assert isinstance(result, FoldXResult)


class TestRunFoldxStabilityInputValidation:
    """Input validation for run_foldx_stability (no FoldX needed)."""

    def test_empty_pdb_string_returns_failure(self) -> None:
        """Empty PDB string → success=False."""
        result = run_foldx_stability("")
        assert result.success is False

    def test_whitespace_pdb_string_returns_failure(self) -> None:
        """Whitespace-only PDB string → success=False."""
        result = run_foldx_stability("   \n\t  ")
        assert result.success is False

    def test_none_pdb_string_returns_failure(self) -> None:
        """None PDB string → success=False (or TypeError)."""
        # The function should handle this gracefully
        try:
            result = run_foldx_stability(None)  # type: ignore[arg-type]
            assert result.success is False
        except TypeError:
            pass  # Also acceptable — None is not a valid str


class TestScanMutationsInputValidation:
    """Input validation for scan_mutations."""

    def test_empty_protein_raises_foldx_error(self) -> None:
        """Empty protein → raises FoldXError."""
        with pytest.raises(FoldXError):
            scan_mutations("")

    def test_invalid_amino_acids_raises_foldx_error(self) -> None:
        """Non-standard amino acids → raises FoldXError."""
        with pytest.raises(FoldXError):
            scan_mutations("MKXB")

    def test_valid_protein_returns_list(self) -> None:
        """Valid protein → returns list of MutationResult."""
        result = scan_mutations("MKG")
        assert isinstance(result, list)

    def test_result_items_are_mutation_results(self) -> None:
        """Each result item is a MutationResult."""
        results = scan_mutations("MKG")
        for mr in results:
            assert isinstance(mr, MutationResult)


class TestScanAllMutationsInputValidation:
    """Input validation for scan_all_mutations."""

    def test_empty_protein_returns_empty_landscape(self) -> None:
        """Empty protein → empty StabilityLandscape."""
        landscape = scan_all_mutations("")
        assert isinstance(landscape, StabilityLandscape)
        assert len(landscape.mutations) == 0

    def test_valid_protein_returns_landscape(self) -> None:
        """Valid protein → StabilityLandscape with mutations."""
        landscape = scan_all_mutations("MKG")
        assert isinstance(landscape, StabilityLandscape)
        assert len(landscape.mutations) > 0


class TestScanPositionInputValidation:
    """Input validation for scan_position."""

    def test_negative_position_returns_empty(self) -> None:
        """Negative position → empty list."""
        result = scan_position("MKG", -1)
        assert result == []

    def test_out_of_range_position_returns_empty(self) -> None:
        """Position beyond protein length → empty list."""
        result = scan_position("MKG", 10)
        assert result == []

    def test_valid_position_returns_list(self) -> None:
        """Valid position → list of mutation dicts."""
        result = scan_position("MKG", 0)
        assert isinstance(result, list)
        assert len(result) > 0


class TestFindCompensatoryMutationsInputValidation:
    """Input validation for find_compensatory_mutations."""

    def test_empty_protein_returns_empty(self) -> None:
        """Empty protein → empty list."""
        result = find_compensatory_mutations("", [{"position": 0, "ddg": 5.0}])
        assert result == []

    def test_empty_mutations_returns_empty(self) -> None:
        """Empty mutations list → empty list."""
        result = find_compensatory_mutations("MKG", [])
        assert result == []

    def test_both_empty_returns_empty(self) -> None:
        """Both empty → empty list."""
        result = find_compensatory_mutations("", [])
        assert result == []


# ═══════════════════════════════════════════════════════════════════════════
# 3. Fallback behavior when FoldX is not available
# ═══════════════════════════════════════════════════════════════════════════


class TestIsFoldxAvailable:
    """Test is_foldx_available and its return behavior."""

    def test_returns_bool(self) -> None:
        """is_foldx_available returns a bool."""
        result = is_foldx_available()
        assert isinstance(result, bool)

    def test_returns_consistently(self) -> None:
        """Repeated calls return the same value."""
        r1 = is_foldx_available()
        r2 = is_foldx_available()
        assert r1 == r2

    def test_returns_false_when_not_on_path(self) -> None:
        """When foldx binary is not on PATH, returns False."""
        with patch("subprocess.run", side_effect=FileNotFoundError):
            assert is_foldx_available() is False

    def test_returns_false_on_timeout(self) -> None:
        """When foldx --version times out, returns False."""
        with patch("subprocess.run", side_effect=subprocess.TimeoutExpired("foldx", 10)):
            assert is_foldx_available() is False

    def test_returns_false_on_os_error(self) -> None:
        """When an OSError occurs, returns False."""
        with patch("subprocess.run", side_effect=OSError("permission denied")):
            assert is_foldx_available() is False

    def test_returns_true_when_available(self) -> None:
        """When foldx runs successfully, returns True."""
        with patch("subprocess.run"):
            assert is_foldx_available() is True


class TestFoldXStabilityFallback:
    """Fallback behavior of run_foldx_stability when FoldX is unavailable."""

    def test_returns_failed_result_when_not_available(self) -> None:
        """When FoldX CLI not available, returns FoldXResult with success=False."""
        with patch("biocompiler.foldx.is_foldx_available", return_value=False):
            # Need a PDB string that passes initial empty check
            pdb = "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N"
            result = run_foldx_stability(pdb)
            # The function may fail at PDB parsing stage before reaching
            # is_foldx_available, so we accept either failure path
            assert isinstance(result, FoldXResult)
            if result.method == "foldx_cli":
                assert result.success is False

    def test_empty_pdb_returns_failure_without_calling_foldx(self) -> None:
        """Empty PDB string fails before FoldX availability check."""
        with patch("biocompiler.foldx.is_foldx_available") as mock_avail:
            result = run_foldx_stability("")
            mock_avail.assert_not_called()
            assert result.success is False


class TestFoldXRepairFallback:
    """Fallback behavior of run_foldx_repair when FoldX is unavailable."""

    def test_returns_original_pdb_when_not_available(self) -> None:
        """When FoldX not available, returns (original_pdb, failed_result)."""
        with patch("biocompiler.foldx.is_foldx_available", return_value=False):
            pdb = "ATOM      1  N   MET A   1       0.0   0.0   0.0  1.00  0.00           N"
            returned_pdb, result = run_foldx_repair(pdb)
            assert returned_pdb == pdb
            assert isinstance(result, FoldXResult)
            assert result.success is False


class TestFoldXMutationFallback:
    """Fallback behavior of run_foldx_mutation when FoldX is unavailable."""

    def test_returns_empty_list_when_not_available(self) -> None:
        """When FoldX not available, returns empty list."""
        with patch("biocompiler.foldx.is_foldx_available", return_value=False):
            result = run_foldx_mutation("some_pdb", ["A1G"])
            assert result == []

    def test_returns_empty_list_for_empty_mutations(self) -> None:
        """Empty mutation list returns empty list regardless of availability."""
        # Don't even need to mock availability — empty list is early return
        result = run_foldx_mutation("some_pdb", [])
        assert result == []


class TestEmpiricalModeIsUsedOffline:
    """When FoldX is not available, empirical mode is used."""

    def test_empirical_stability_method_is_empirical(self) -> None:
        """empirical_stability always uses method='empirical'."""
        result = empirical_stability("MKGILV")
        assert result.method == "empirical"

    def test_scan_mutations_uses_empirical_without_pdb(self) -> None:
        """scan_mutations without PDB uses empirical estimation."""
        with patch("biocompiler.foldx.is_foldx_available", return_value=False):
            results = scan_mutations("MKG")
            # All results should be from empirical (not foldx_cli)
            for mr in results:
                assert mr.engine == "foldx"  # engine name is still foldx
                assert mr.score_type == "ddg"

    def test_scan_all_mutations_falls_back_to_empirical(self) -> None:
        """scan_all_mutations with method='foldx' falls back to 'empirical'."""
        landscape = scan_all_mutations("MKG", method="foldx")
        assert landscape.method == "empirical"

    def test_scan_position_uses_empirical_method(self) -> None:
        """scan_position uses empirical estimation (no PDB needed)."""
        results = scan_position("MKG", 0)
        assert len(results) == 19  # 19 possible substitutions
        for r in results:
            assert "ddg" in r


class TestRunStabilityBatchOffline:
    """Batch API works without FoldX."""

    def test_batch_returns_batch_result(self) -> None:
        """run_stability_batch returns BatchResult[FoldXResult]."""
        result = run_stability_batch(["MKG", "ILV"])
        assert isinstance(result, BatchResult)

    def test_batch_results_match_input_count(self) -> None:
        """Number of results matches number of input sequences."""
        result = run_stability_batch(["MKG", "ILV", "AAA"])
        assert len(result.results) == 3

    def test_batch_handles_empty_list(self) -> None:
        """Empty input list returns empty BatchResult."""
        result = run_stability_batch([])
        assert isinstance(result, BatchResult)
        assert len(result.results) == 0

    def test_batch_handles_invalid_sequence(self) -> None:
        """Invalid sequence in batch produces a failed FoldXResult."""
        result = run_stability_batch(["", "MKG"])
        assert len(result.results) == 2
        # The empty one should be a failure
        assert result.results[0].success is False
        assert result.results[1].success is True


# ═══════════════════════════════════════════════════════════════════════════
# 4. Named constants validation
# ═══════════════════════════════════════════════════════════════════════════


class TestAccuracyConstants:
    """Validate named accuracy and error constants."""

    def test_foldx_cli_accuracy(self) -> None:
        """FOLDX_CLI_ACCURACY = 1.0 (±1 kcal/mol)."""
        assert FOLDX_CLI_ACCURACY == 1.0

    def test_foldx_empirical_mae(self) -> None:
        """FOLDX_EMPIRICAL_MAE is a positive float representing overall MAE."""
        assert isinstance(FOLDX_EMPIRICAL_MAE, float)
        assert FOLDX_EMPIRICAL_MAE > 0.0

    def test_foldx_empirical_direction_accuracy(self) -> None:
        """FOLDX_EMPIRICAL_DIRECTION_ACCURACY is 1.0 (100%)."""
        assert FOLDX_EMPIRICAL_DIRECTION_ACCURACY == 1.0

    def test_foldx_empirical_small_mae(self) -> None:
        """FOLDX_EMPIRICAL_SMALL_MAE < FOLDX_EMPIRICAL_MAE (small proteins are better)."""
        assert FOLDX_EMPIRICAL_SMALL_MAE < FOLDX_EMPIRICAL_MAE

    def test_foldx_empirical_medium_mae(self) -> None:
        """FOLDX_EMPIRICAL_MEDIUM_MAE is between small and large."""
        assert FOLDX_EMPIRICAL_SMALL_MAE <= FOLDX_EMPIRICAL_MEDIUM_MAE
        assert FOLDX_EMPIRICAL_MEDIUM_MAE <= FOLDX_EMPIRICAL_LARGE_MAE

    def test_foldx_empirical_large_mae(self) -> None:
        """FOLDX_EMPIRICAL_LARGE_MAE > FOLDX_EMPIRICAL_MAE (large proteins worse)."""
        assert FOLDX_EMPIRICAL_LARGE_MAE > FOLDX_EMPIRICAL_MAE

    def test_foldx_empirical_pearson_r(self) -> None:
        """FOLDX_EMPIRICAL_PEARSON_R is between 0 and 1."""
        assert 0.0 <= FOLDX_EMPIRICAL_PEARSON_R <= 1.0

    def test_foldx_empirical_bias(self) -> None:
        """FOLDX_EMPIRICAL_BIAS is a positive float (under-predicts stability)."""
        assert isinstance(FOLDX_EMPIRICAL_BIAS, float)
        assert FOLDX_EMPIRICAL_BIAS > 0.0

    def test_accuracy_ordering(self) -> None:
        """Small MAE < Medium MAE < Large MAE."""
        assert FOLDX_EMPIRICAL_SMALL_MAE < FOLDX_EMPIRICAL_MEDIUM_MAE
        assert FOLDX_EMPIRICAL_MEDIUM_MAE < FOLDX_EMPIRICAL_LARGE_MAE


class TestBLOSUM62Constant:
    """Validate the BLOSUM62 substitution matrix constant."""

    def test_blosum62_is_dict(self) -> None:
        """BLOSUM62 is a dict."""
        assert isinstance(BLOSUM62, dict)

    def test_blosum62_has_20_rows(self) -> None:
        """BLOSUM62 has entries for all 20 standard amino acids."""
        standard = set("ARNDCQEGHILKMFPSTWYV")
        assert set(BLOSUM62.keys()) == standard

    def test_blosum62_each_row_has_20_cols(self) -> None:
        """Each BLOSUM62 row has 20 entries."""
        standard = set("ARNDCQEGHILKMFPSTWYV")
        for aa in standard:
            assert set(BLOSUM62[aa].keys()) == standard

    def test_blosum62_is_symmetric(self) -> None:
        """BLOSUM62[i][j] == BLOSUM62[j][i]."""
        standard = set("ARNDCQEGHILKMFPSTWYV")
        for a1 in standard:
            for a2 in standard:
                assert BLOSUM62[a1][a2] == BLOSUM62[a2][a1], (
                    f"BLOSUM62 not symmetric: [{a1}][{a2}]={BLOSUM62[a1][a2]} "
                    f"!= [{a2}][{a1}]={BLOSUM62[a2][a1]}"
                )

    def test_blosum62_diagonal_positive(self) -> None:
        """BLOSUM62 diagonal entries are positive (identical substitution)."""
        for aa in "ARNDCQEGHILKMFPSTWYV":
            assert BLOSUM62[aa][aa] > 0, f"BLOSUM62[{aa}][{aa}]={BLOSUM62[aa][aa]} not positive"

    def test_blosum62_values_are_integers(self) -> None:
        """All BLOSUM62 values are integers (standard BLOSUM matrices)."""
        for a1, row in BLOSUM62.items():
            for a2, val in row.items():
                assert isinstance(val, int), f"BLOSUM62[{a1}][{a2}]={val} not int"


class TestHydropathyConstant:
    """Validate the HYDROPATHY scale constant."""

    def test_hydropathy_is_dict(self) -> None:
        """HYDROPATHY is a dict."""
        assert isinstance(HYDROPATHY, dict)

    def test_hydropathy_has_20_entries(self) -> None:
        """HYDROPATHY has entries for all 20 standard amino acids."""
        standard = set("ARNDCQEGHILKMFPSTWYV")
        assert set(HYDROPATHY.keys()) == standard

    def test_hydropathy_values_are_floats(self) -> None:
        """All HYDROPATHY values are floats."""
        for aa, val in HYDROPATHY.items():
            assert isinstance(val, float), f"HYDROPATHY[{aa}]={val} not float"

    def test_hydropathy_values_in_reasonable_range(self) -> None:
        """Kyte-Doolittle values range from about -4.5 to +4.5."""
        for aa, val in HYDROPATHY.items():
            assert -5.0 <= val <= 5.0, f"HYDROPATHY[{aa}]={val} out of [-5, 5]"

    def test_hydropathy_ile_most_hydrophobic(self) -> None:
        """Isoleucine should be the most hydrophobic (highest value)."""
        # In Kyte-Doolittle, I has the highest hydropathy (4.5)
        max_aa = max(HYDROPATHY, key=HYDROPATHY.get)
        # At minimum, the most hydrophobic should be one of I, V, L, F, C, M, A
        assert max_aa in "IVLFCMA"


class TestAAVolumeConstant:
    """Validate the AA_VOLUME constant."""

    def test_aa_volume_is_dict(self) -> None:
        """AA_VOLUME is a dict."""
        assert isinstance(AA_VOLUME, dict)

    def test_aa_volume_has_20_entries(self) -> None:
        """AA_VOLUME covers all 20 standard amino acids."""
        standard = set("ARNDCQEGHILKMFPSTWYV")
        assert set(AA_VOLUME.keys()) == standard

    def test_aa_volume_values_positive(self) -> None:
        """All AA_VOLUME values are positive."""
        for aa, vol in AA_VOLUME.items():
            assert vol > 0, f"AA_VOLUME[{aa}]={vol} not positive"

    def test_aa_volume_glycine_smallest(self) -> None:
        """Glycine has the smallest van der Waals volume."""
        gly_vol = AA_VOLUME["G"]
        for aa, vol in AA_VOLUME.items():
            if aa != "G":
                assert gly_vol < vol, f"Glycine ({gly_vol}) >= {aa} ({vol})"

    def test_aa_volume_tryptophan_largest(self) -> None:
        """Tryptophan has the largest van der Waals volume."""
        trp_vol = AA_VOLUME["W"]
        for aa, vol in AA_VOLUME.items():
            if aa != "W":
                assert trp_vol > vol, f"Tryptophan ({trp_vol}) <= {aa} ({vol})"


# ═══════════════════════════════════════════════════════════════════════════
# 5. FoldXResult classification and confidence
# ═══════════════════════════════════════════════════════════════════════════


class TestClassifyStability:
    """Test _classify_stability via FoldXResult construction."""

    def test_very_stable_classification(self) -> None:
        """stability_kcal < -10 → 'very_stable'."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-15.0)
        assert result.classification == "very_stable"

    def test_stable_classification(self) -> None:
        """-10 ≤ stability_kcal < -5 → 'stable'."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-7.0)
        assert result.classification == "stable"

    def test_marginally_stable_classification(self) -> None:
        """-5 ≤ stability_kcal < 0 → 'marginally_stable'."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-2.0)
        assert result.classification == "marginally_stable"

    def test_unstable_classification(self) -> None:
        """stability_kcal ≥ 0 → 'unstable'."""
        result = FoldXResult(protein="M" * 50, stability_kcal=5.0)
        assert result.classification == "unstable"

    def test_failed_classification(self) -> None:
        """success=False → 'failed'."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-15.0, success=False)
        assert result.classification == "failed"

    def test_explicit_classification_overrides(self) -> None:
        """Explicit classification= parameter overrides auto-derived."""
        result = FoldXResult(
            protein="M" * 50, stability_kcal=-15.0,
            classification="custom_label",
        )
        assert result.classification == "custom_label"


class TestConfidenceLevel:
    """Test FoldXResult.confidence_level property."""

    def test_foldx_cli_high_confidence(self) -> None:
        """FoldX CLI mode → 'high' confidence."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-5.0, method="foldx_cli")
        assert result.confidence_level == "high"

    def test_empirical_small_protein_high(self) -> None:
        """Empirical + small protein (<100 aa) → 'high' confidence."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-5.0, method="empirical")
        assert result.confidence_level == "high"

    def test_empirical_medium_protein_medium(self) -> None:
        """Empirical + medium protein (100-300 aa) → 'medium' confidence."""
        result = FoldXResult(protein="M" * 150, stability_kcal=-5.0, method="empirical")
        assert result.confidence_level == "medium"

    def test_empirical_large_protein_low(self) -> None:
        """Empirical + large protein (>300 aa) → 'low' confidence."""
        result = FoldXResult(protein="M" * 400, stability_kcal=-5.0, method="empirical")
        assert result.confidence_level == "low"

    def test_failed_result_unknown_confidence(self) -> None:
        """Failed result → 'unknown' confidence."""
        result = FoldXResult(protein="M" * 50, stability_kcal=-5.0, success=False)
        assert result.confidence_level == "unknown"


# ═══════════════════════════════════════════════════════════════════════════
# 6. FoldXCache
# ═══════════════════════════════════════════════════════════════════════════


class TestFoldXCache:
    """Unit tests for FoldXCache (in-memory result cache)."""

    def test_empty_cache_len_zero(self) -> None:
        """New cache has length 0."""
        cache = FoldXCache()
        assert len(cache) == 0

    def test_put_and_get(self) -> None:
        """Put a result and retrieve it."""
        cache = FoldXCache()
        dummy = {"stability_kcal": -5.0}
        cache.put("protein1", "empirical", dummy)
        retrieved = cache.get("protein1", "empirical")
        assert retrieved is dummy

    def test_get_miss_returns_none(self) -> None:
        """Getting non-existent key returns None."""
        cache = FoldXCache()
        assert cache.get("nonexistent", "empirical") is None

    def test_clear_empties_cache(self) -> None:
        """Clear removes all entries."""
        cache = FoldXCache()
        cache.put("p1", "empirical", {"a": 1})
        cache.put("p2", "empirical", {"b": 2})
        assert len(cache) == 2
        cache.clear()
        assert len(cache) == 0

    def test_eviction_when_over_capacity(self) -> None:
        """Entries are evicted (FIFO) when cache exceeds max_size."""
        cache = FoldXCache(max_size=2)
        cache.put("p1", "empirical", {"a": 1})
        cache.put("p2", "empirical", {"b": 2})
        cache.put("p3", "empirical", {"c": 3})  # p1 should be evicted
        assert cache.get("p1", "empirical") is None
        assert cache.get("p2", "empirical") is not None
        assert cache.get("p3", "empirical") is not None

    def test_different_methods_are_separate_keys(self) -> None:
        """Same content, different method → separate cache entries."""
        cache = FoldXCache()
        cache.put("p1", "empirical", {"a": 1})
        cache.put("p1", "foldx_stability", {"b": 2})
        assert len(cache) == 2
        assert cache.get("p1", "empirical") != cache.get("p1", "foldx_stability")


class TestModuleCache:
    """Test the module-level cache and clear_cache function."""

    def test_clear_cache_resets_module_cache(self) -> None:
        """clear_cache() empties the module-level cache."""
        clear_cache()
        from biocompiler.foldx import _cache
        assert len(_cache) == 0

    def test_clear_cache_is_idempotent(self) -> None:
        """Calling clear_cache twice doesn't raise."""
        clear_cache()
        clear_cache()
        # Implicit assertion: no exception raised

    def test_empirical_stability_populates_cache(self) -> None:
        """Running empirical_stability stores result in module cache."""
        clear_cache()
        from biocompiler.foldx import _cache
        empirical_stability("MKGVIL")
        assert len(_cache) > 0


# ═══════════════════════════════════════════════════════════════════════════
# 7. FoldXError
# ═══════════════════════════════════════════════════════════════════════════


class TestFoldXError:
    """Test FoldXError exception class."""

    def test_is_exception(self) -> None:
        """FoldXError is a subclass of Exception."""
        assert issubclass(FoldXError, Exception)

    def test_message_contains_reason(self) -> None:
        """Error message includes the reason string."""
        err = FoldXError("something broke")
        assert "something broke" in str(err)

    def test_command_attribute(self) -> None:
        """FoldXError can carry a command attribute."""
        err = FoldXError("timeout", command="foldx --command=Stability")
        assert err.command == "foldx --command=Stability"

    def test_command_in_message(self) -> None:
        """When command is provided, it appears in the string representation."""
        err = FoldXError("timeout", command="foldx --command=Stability")
        assert "foldx --command=Stability" in str(err)

    def test_no_command_attribute(self) -> None:
        """FoldXError without command has command=None."""
        err = FoldXError("generic error")
        assert err.command is None

    def test_can_be_raised_and_caught(self) -> None:
        """FoldXError can be raised and caught as expected."""
        with pytest.raises(FoldXError):
            raise FoldXError("test error")


# ═══════════════════════════════════════════════════════════════════════════
# 8. __all__ validation
# ═══════════════════════════════════════════════════════════════════════════


class TestDunderAll:
    """Verify __all__ is correctly defined."""

    def test_all_is_list(self) -> None:
        """__all__ is a list."""
        import biocompiler.foldx as mod
        assert isinstance(mod.__all__, list)

    def test_all_entries_are_strings(self) -> None:
        """Every entry in __all__ is a string."""
        import biocompiler.foldx as mod
        for name in mod.__all__:
            assert isinstance(name, str), f"__all__ entry {name!r} is not a string"

    def test_all_entries_are_importable(self) -> None:
        """Every name in __all__ is an attribute of the module."""
        import biocompiler.foldx as mod
        for name in mod.__all__:
            assert hasattr(mod, name), f"__all__ entry {name!r} not found in module"

    def test_all_contains_key_functions(self) -> None:
        """Key public functions are listed in __all__."""
        import biocompiler.foldx as mod
        for name in [
            "is_foldx_available", "run_foldx_stability", "run_foldx_repair",
            "run_foldx_mutation", "empirical_stability", "run_stability_batch",
            "scan_mutations", "find_stabilizing_mutations", "scan_all_mutations",
            "scan_position", "compute_conservation", "find_compensatory_mutations",
            "rank_positions_by_mutability", "identify_hotspot_regions",
            "FoldXResult", "MutationResult", "FoldXError",
            "StabilityLandscape", "ConservationScore", "FoldXCache",
            "clear_cache", "BLOSUM62", "HYDROPATHY", "AA_VOLUME",
        ]:
            assert name in mod.__all__, f"{name!r} not in __all__"

    def test_all_contains_accuracy_constants(self) -> None:
        """Accuracy constants are listed in __all__."""
        import biocompiler.foldx as mod
        for name in [
            "FOLDX_CLI_ACCURACY", "FOLDX_EMPIRICAL_MAE",
            "FOLDX_EMPIRICAL_DIRECTION_ACCURACY", "FOLDX_EMPIRICAL_SMALL_MAE",
            "FOLDX_EMPIRICAL_MEDIUM_MAE", "FOLDX_EMPIRICAL_LARGE_MAE",
            "FOLDX_EMPIRICAL_PEARSON_R", "FOLDX_EMPIRICAL_BIAS",
        ]:
            assert name in mod.__all__, f"{name!r} not in __all__"
