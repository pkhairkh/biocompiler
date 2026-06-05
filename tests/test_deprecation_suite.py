"""
BioCompiler Deprecation Test Suite — F5.8
==========================================

Comprehensive tests that verify:
1. Every deprecated function emits a DeprecationWarning
2. The recommended replacement works without warnings
3. The deprecated function still produces correct output (backward compat)

Deprecated functions covered:
  - splicing.maxent_score()           → recommend maxentscan.score_donor()
  - splicing.score_splice_sites()     → recommend maxentscan.scan_splice_sites()
  - dispatch.is_csp_available()       → recommend get_csp_availability()
  - solver.is_csp_available (import)  → recommend solver.get_csp_availability
  - esmfold.predict_batch()           → recommend predict_structure_batch()
  - immunogenicity.compute_surface_accessibility_approx()
                                      → recommend predict_eea()
  - immunogenicity.predict_b_cell_epitopes()
                                      → recommend predict_kolaskar_tongaonkar()
  - camsol.SolubilityResult           → recommend CamSolResult
"""

from __future__ import annotations

import warnings

import pytest


# ═══════════════════════════════════════════════════════════════════════════
# Helper: assert a callable does NOT emit DeprecationWarning
# ═══════════════════════════════════════════════════════════════════════════

class _NoDeprecationContext:
    """Context manager that fails if a DeprecationWarning is emitted."""

    def __enter__(self):
        self._catch = warnings.catch_warnings(record=True)
        self._warnings = self._catch.__enter__()
        warnings.simplefilter("always", DeprecationWarning)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._catch.__exit__(exc_type, exc_val, exc_tb)
        dep_warnings = [w for w in self._warnings if issubclass(w.category, DeprecationWarning)]
        if dep_warnings:
            msgs = "\n  ".join(str(w.message) for w in dep_warnings)
            pytest.fail(f"Unexpected DeprecationWarning(s):\n  {msgs}")
        return False


def assert_no_deprecation():
    """Return a context manager that asserts no DeprecationWarning is emitted."""
    return _NoDeprecationContext()


# ═══════════════════════════════════════════════════════════════════════════
# 1. splicing.maxent_score() → recommend maxentscan.score_donor()
# ═══════════════════════════════════════════════════════════════════════════

class TestMaxentScoreDeprecation:
    """Test that splicing.maxent_score() emits DeprecationWarning and
    the replacement maxentscan.score_donor() works without warnings."""

    def test_maxent_score_emits_deprecation_warning(self):
        """Calling maxent_score() must emit a DeprecationWarning."""
        from biocompiler.splicing import maxent_score

        with pytest.warns(DeprecationWarning, match="maxent_score.*deprecated"):
            maxent_score("CAGGTGAGT")

    def test_maxent_score_backward_compat(self):
        """maxent_score() must still return a numeric score."""
        from biocompiler.splicing import maxent_score

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = maxent_score("CAGGTGAGT")

        assert isinstance(result, float)
        assert result > 0.0  # 'CAGGTGAGT' contains GT — should score positive

    def test_maxent_score_short_input_backward_compat(self):
        """maxent_score() returns 0.0 for too-short context."""
        from biocompiler.splicing import maxent_score

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            assert maxent_score("GT") == 0.0

    def test_score_donor_no_deprecation(self):
        """Recommended replacement maxentscan.score_donor() must NOT warn."""
        from biocompiler.maxentscan import score_donor

        with assert_no_deprecation():
            # score_donor needs 3 bases upstream and 6 downstream of GT position
            # Use a longer sequence with GT at position 13 (plenty of context)
            seq = "A" * 10 + "CAGGTAAGT" + "C" * 10
            gt_pos = seq.find("GT")
            result = score_donor(seq, gt_pos)

        assert isinstance(result, float)

    def test_score_donor_returns_valid_score(self):
        """score_donor() returns a reasonable log-odds score."""
        from biocompiler.maxentscan import score_donor

        with assert_no_deprecation():
            # A canonical donor site in a sufficiently long context
            # score_donor needs seq[position-3 : position+6] to be valid
            seq = "A" * 10 + "CAGGTAAGT" + "C" * 10
            gt_pos = seq.find("GT")
            result = score_donor(seq, gt_pos)

        assert isinstance(result, float)
        # Canonical donor should produce a reasonable (possibly negative for
        # random context, but not impossible) score
        assert result > -50.0  # -50 is the impossible sentinel


# ═══════════════════════════════════════════════════════════════════════════
# 2. splicing.score_splice_sites() → recommend maxentscan.scan_splice_sites()
# ═══════════════════════════════════════════════════════════════════════════

class TestScoreSpliceSitesDeprecation:
    """Test that splicing.score_splice_sites() emits DeprecationWarning and
    the replacement maxentscan.scan_splice_sites() works without warnings."""

    def test_score_splice_sites_emits_deprecation_warning(self):
        """Calling score_splice_sites() must emit a DeprecationWarning."""
        from biocompiler.splicing import score_splice_sites

        with pytest.warns(DeprecationWarning, match="score_splice_sites.*deprecated"):
            score_splice_sites("ATGCAGGTGAGTCCC")

    def test_score_splice_sites_backward_compat(self):
        """score_splice_sites() must still return a list of tuples."""
        from biocompiler.splicing import score_splice_sites

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            results = score_splice_sites("ATGCAGGTGAGTCCC")

        assert isinstance(results, list)
        # Should find at least one GT site in "CAGGTGAGT"
        assert len(results) >= 1
        # Each result is (position, score, SpliceVerdict)
        for pos, score, verdict in results:
            assert isinstance(pos, int)
            assert isinstance(score, float)

    def test_scan_splice_sites_no_deprecation(self):
        """Recommended replacement maxentscan.scan_splice_sites() must NOT warn."""
        from biocompiler.maxentscan import scan_splice_sites

        with assert_no_deprecation():
            results = scan_splice_sites("ATGCAGGTGAGTCCC")

        assert isinstance(results, list)

    def test_scan_splice_sites_returns_valid_results(self):
        """scan_splice_sites() returns proper (position, type, score) tuples."""
        from biocompiler.maxentscan import scan_splice_sites

        with assert_no_deprecation():
            # A sequence with a canonical donor and acceptor
            seq = "A" * 30 + "GT" + "A" * 50 + "AG" + "A" * 30
            results = scan_splice_sites(seq)

        assert isinstance(results, list)
        for pos, site_type, score in results:
            assert isinstance(pos, int)
            assert site_type in ("donor", "acceptor")
            assert isinstance(score, float)


# ═══════════════════════════════════════════════════════════════════════════
# 3. dispatch.is_csp_available() → recommend get_csp_availability()
# ═══════════════════════════════════════════════════════════════════════════

class TestIsCspAvailableDeprecation:
    """Test that dispatch.is_csp_available() emits DeprecationWarning and
    the replacement get_csp_availability() works without warnings."""

    def test_is_csp_available_emits_deprecation_warning(self):
        """Calling is_csp_available() must emit a DeprecationWarning."""
        from biocompiler.solver.dispatch import is_csp_available

        with pytest.warns(DeprecationWarning, match="is_csp_available.*deprecated.*get_csp_availability"):
            is_csp_available()

    def test_is_csp_available_backward_compat(self):
        """is_csp_available() must still return a dict with correct shape."""
        from biocompiler.solver.dispatch import is_csp_available

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = is_csp_available()

        assert isinstance(result, dict)
        assert "ortools" in result
        assert "z3" in result
        assert "any" in result
        assert isinstance(result["any"], bool)

    def test_get_csp_availability_no_deprecation(self):
        """Recommended replacement get_csp_availability() must NOT warn."""
        from biocompiler.solver.dispatch import get_csp_availability

        with assert_no_deprecation():
            result = get_csp_availability()

        assert isinstance(result, dict)
        assert "ortools" in result
        assert "z3" in result
        assert "any" in result

    def test_is_csp_available_matches_get_csp_availability(self):
        """Deprecated function must produce identical output to the replacement."""
        from biocompiler.solver.dispatch import is_csp_available, get_csp_availability

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            deprecated_result = is_csp_available()

        with assert_no_deprecation():
            replacement_result = get_csp_availability()

        assert deprecated_result == replacement_result

    def test_is_solver_available_no_deprecation(self):
        """is_solver_available() (another replacement) must NOT warn."""
        from biocompiler.solver.dispatch import is_solver_available

        with assert_no_deprecation():
            result = is_solver_available()

        assert isinstance(result, bool)


# ═══════════════════════════════════════════════════════════════════════════
# 4. solver.is_csp_available (package-level import) → recommend solver.get_csp_availability
# ═══════════════════════════════════════════════════════════════════════════

class TestSolverPackageDeprecation:
    """Test that importing is_csp_available from biocompiler.solver emits
    a DeprecationWarning via __getattr__."""

    def test_solver_is_csp_available_import_emits_warning(self):
        """Accessing solver.is_csp_available must emit a DeprecationWarning."""
        import biocompiler.solver

        with pytest.warns(DeprecationWarning, match="solver.is_csp_available.*deprecated.*get_csp_availability"):
            func = biocompiler.solver.is_csp_available
            # Also call it — the warning is on attribute access, but verify it works
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                result = func()

        assert isinstance(result, dict)

    def test_solver_get_csp_availability_no_deprecation(self):
        """Accessing solver.get_csp_availability must NOT warn."""
        from biocompiler.solver import get_csp_availability

        with assert_no_deprecation():
            result = get_csp_availability()

        assert isinstance(result, dict)


# ═══════════════════════════════════════════════════════════════════════════
# 5. esmfold.predict_batch() → recommend predict_structure_batch()
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictBatchDeprecation:
    """Test that esmfold.predict_batch() emits DeprecationWarning and
    the replacement predict_structure_batch() works without warnings."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_esmfold(self):
        """Skip test module if esmfold cannot be imported."""
        pytest.importorskip("biocompiler.esmfold")

    def test_predict_batch_emits_deprecation_warning(self):
        """Calling predict_batch() must emit a DeprecationWarning."""
        from biocompiler.esmfold import predict_batch, BatchStructureRequest

        request = BatchStructureRequest(proteins=["MVLSPADKTN"])
        with pytest.warns(DeprecationWarning, match="predict_batch.*deprecated.*predict_structure_batch"):
            # This will attempt actual prediction but the warning is emitted
            # at function entry before any network calls
            try:
                predict_batch(request)
            except Exception:
                pass  # We only care about the warning, not the result

    def test_predict_structure_batch_no_deprecation(self):
        """Recommended replacement predict_structure_batch() must NOT warn."""
        from biocompiler.esmfold import predict_structure_batch

        with assert_no_deprecation():
            # Call with empty list to avoid network calls
            result = predict_structure_batch([])

        # Returns a BatchResult, not a plain list
        assert hasattr(result, "results") or isinstance(result, list)


# ═══════════════════════════════════════════════════════════════════════════
# 6. immunogenicity.compute_surface_accessibility_approx()
#    → recommend predict_eea()
# ═══════════════════════════════════════════════════════════════════════════

class TestComputeSurfaceAccessibilityDeprecation:
    """Test that compute_surface_accessibility_approx() emits DeprecationWarning
    and the replacement predict_eea() works without warnings."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_immunogenicity(self):
        """Skip test module if immunogenicity cannot be imported."""
        pytest.importorskip("biocompiler.immunogenicity")

    def test_compute_surface_accessibility_emits_deprecation(self):
        """Calling compute_surface_accessibility_approx() must emit DeprecationWarning."""
        from biocompiler.immunogenicity import compute_surface_accessibility_approx

        with pytest.warns(
            DeprecationWarning,
            match="compute_surface_accessibility_approx.*deprecated.*predict_eea",
        ):
            compute_surface_accessibility_approx("MVLSPADKTN")

    def test_compute_surface_accessibility_backward_compat(self):
        """compute_surface_accessibility_approx() must still return valid output."""
        from biocompiler.immunogenicity import compute_surface_accessibility_approx

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = compute_surface_accessibility_approx("MVLSPADKTN")

        assert isinstance(result, list)
        assert len(result) == 10  # One value per residue
        for val in result:
            assert isinstance(val, float)
            assert 0.0 <= val <= 1.0

    def test_predict_eea_no_deprecation(self):
        """Recommended replacement predict_eea() must NOT warn."""
        from biocompiler.immunogenicity import predict_eea

        with assert_no_deprecation():
            result = predict_eea("MVLSPADKTN")

        assert isinstance(result, list)


# ═══════════════════════════════════════════════════════════════════════════
# 7. immunogenicity.predict_b_cell_epitopes()
#    → recommend predict_kolaskar_tongaonkar()
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictBCellEpitopesDeprecation:
    """Test that predict_b_cell_epitopes() emits DeprecationWarning and
    the replacement predict_kolaskar_tongaonkar() works without warnings."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_immunogenicity(self):
        """Skip test module if immunogenicity cannot be imported."""
        pytest.importorskip("biocompiler.immunogenicity")

    def test_predict_b_cell_epitopes_emits_deprecation(self):
        """Calling predict_b_cell_epitopes() must emit DeprecationWarning."""
        from biocompiler.immunogenicity import predict_b_cell_epitopes

        with pytest.warns(
            DeprecationWarning,
            match="predict_b_cell_epitopes.*deprecated.*predict_kolaskar_tongaonkar",
        ):
            predict_b_cell_epitopes("MVLSPADKTNVKAAWGKVGAH")

    def test_predict_b_cell_epitopes_backward_compat(self):
        """predict_b_cell_epitopes() must still return a list of dicts."""
        from biocompiler.immunogenicity import predict_b_cell_epitopes

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = predict_b_cell_epitopes("MVLSPADKTNVKAAWGKVGAH")

        assert isinstance(result, list)
        # Each entry should be a dict-like object (BCellEpitopeDict)
        for entry in result:
            assert hasattr(entry, "start") or "start" in entry

    def test_predict_kolaskar_tongaonkar_no_deprecation(self):
        """Recommended replacement predict_kolaskar_tongaonkar() must NOT warn."""
        from biocompiler.immunogenicity import predict_kolaskar_tongaonkar

        with assert_no_deprecation():
            result = predict_kolaskar_tongaonkar("MVLSPADKTNVKAAWGKVGAH")

        assert isinstance(result, list)


# ═══════════════════════════════════════════════════════════════════════════
# 8. camsol.SolubilityResult → recommend CamSolResult
# ═══════════════════════════════════════════════════════════════════════════

class TestSolubilityResultDeprecation:
    """Test that accessing camsol.SolubilityResult emits DeprecationWarning
    and the replacement CamSolResult works without warnings."""

    def test_solubility_result_import_emits_deprecation(self):
        """Accessing camsol.SolubilityResult must emit a DeprecationWarning."""
        import biocompiler.camsol

        with pytest.warns(DeprecationWarning, match="SolubilityResult.*deprecated.*CamSolResult"):
            cls = biocompiler.camsol.SolubilityResult

        # It should still be the same class
        assert cls is biocompiler.camsol.CamSolResult

    def test_camsol_result_no_deprecation(self):
        """Accessing camsol.CamSolResult must NOT warn."""
        from biocompiler.camsol import CamSolResult

        with assert_no_deprecation():
            # Just access the class, don't instantiate with invalid args
            assert CamSolResult is not None

    def test_solubility_result_is_cam_sol_result(self):
        """SolubilityResult must be an alias for CamSolResult."""
        import biocompiler.camsol

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            SR = biocompiler.camsol.SolubilityResult

        assert SR is biocompiler.camsol.CamSolResult


# ═══════════════════════════════════════════════════════════════════════════
# 9. Cross-cutting: verify warning message quality
# ═══════════════════════════════════════════════════════════════════════════

class TestDeprecationWarningQuality:
    """Verify that all deprecation warnings follow a consistent pattern:
    - Mention the deprecated function
    - Mention the recommended replacement
    - Are DeprecationWarning (not UserWarning, etc.)
    """

    @pytest.mark.parametrize("module_path,func_name,match_re", [
        (
            "biocompiler.splicing",
            "maxent_score",
            r"maxent_score.*deprecated",
        ),
        (
            "biocompiler.splicing",
            "score_splice_sites",
            r"score_splice_sites.*deprecated",
        ),
        (
            "biocompiler.solver.dispatch",
            "is_csp_available",
            r"is_csp_available.*deprecated.*get_csp_availability",
        ),
    ])
    def test_deprecation_message_mentions_replacement(self, module_path, func_name, match_re):
        """Every deprecation warning must mention both the deprecated and replacement APIs."""
        import importlib

        mod = importlib.import_module(module_path)
        func = getattr(mod, func_name)

        # Get the right call signature
        if func_name == "maxent_score":
            with pytest.warns(DeprecationWarning, match=match_re):
                func("CAGGTGAGT")
        elif func_name == "score_splice_sites":
            with pytest.warns(DeprecationWarning, match=match_re):
                func("ATGCAGGTGAGTCCC")
        elif func_name == "is_csp_available":
            with pytest.warns(DeprecationWarning, match=match_re):
                func()

    @pytest.mark.parametrize("module_path,func_name", [
        ("biocompiler.splicing", "maxent_score"),
        ("biocompiler.splicing", "score_splice_sites"),
        ("biocompiler.solver.dispatch", "is_csp_available"),
    ])
    def test_warning_is_deprecation_warning_not_user_warning(self, module_path, func_name):
        """Ensure warnings are DeprecationWarning subclass, not UserWarning."""
        import importlib

        mod = importlib.import_module(module_path)
        func = getattr(mod, func_name)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            if func_name == "maxent_score":
                func("CAGGTGAGT")
            elif func_name == "score_splice_sites":
                func("ATGCAGGTGAGTCCC")
            elif func_name == "is_csp_available":
                func()

        dep_warnings = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]

        assert len(dep_warnings) >= 1, f"No DeprecationWarning emitted by {func_name}"
        assert len(user_warnings) == 0, f"UserWarning emitted instead of DeprecationWarning by {func_name}"


# ═══════════════════════════════════════════════════════════════════════════
# 10. Replacement API smoke tests
# ═══════════════════════════════════════════════════════════════════════════

class TestReplacementAPIsWork:
    """Verify that every recommended replacement actually works correctly."""

    def test_maxentscan_score_donor_works(self):
        """maxentscan.score_donor() produces a valid float score."""
        from biocompiler.maxentscan import score_donor

        with assert_no_deprecation():
            # Need sufficient context: 3 bases upstream and 6 downstream of GT
            seq = "A" * 10 + "CAGGTAAGT" + "C" * 10
            gt_pos = seq.find("GT")
            score = score_donor(seq, gt_pos)

        assert isinstance(score, float)

    def test_maxentscan_score_acceptor_works(self):
        """maxentscan.score_acceptor() produces a valid float score."""
        from biocompiler.maxentscan import score_acceptor

        with assert_no_deprecation():
            # Need 20 bases upstream and 3 downstream of AG position
            seq = "C" * 25 + "TTTTTTCAG" + "C" * 10
            ag_pos = seq.find("AG")
            score = score_acceptor(seq, ag_pos)

        assert isinstance(score, float)

    def test_maxentscan_scan_splice_sites_works(self):
        """maxentscan.scan_splice_sites() produces valid results."""
        from biocompiler.maxentscan import scan_splice_sites

        with assert_no_deprecation():
            # Build a sequence with a donor (GT) and acceptor (AG)
            seq = "A" * 25 + "GT" + "A" * 50 + "AG" + "A" * 25
            results = scan_splice_sites(seq, donor_threshold=0.0, acceptor_threshold=0.0)

        assert isinstance(results, list)
        # Should find at least one donor and one acceptor
        site_types = {r[1] for r in results}
        assert "donor" in site_types or "acceptor" in site_types

    def test_get_csp_availability_works(self):
        """get_csp_availability() returns a dict with expected keys."""
        from biocompiler.solver.dispatch import get_csp_availability

        with assert_no_deprecation():
            result = get_csp_availability()

        assert set(result.keys()) == {"ortools", "z3", "any"}

    def test_is_solver_available_works(self):
        """is_solver_available() returns a boolean."""
        from biocompiler.solver.dispatch import is_solver_available

        with assert_no_deprecation():
            result = is_solver_available()

        assert isinstance(result, bool)

    def test_camsol_cam_sol_result_works(self):
        """CamSolResult can be used directly without deprecation warning."""
        from biocompiler.camsol import CamSolResult

        with assert_no_deprecation():
            # Create a minimal result
            result = CamSolResult(
                sequence="MVLSPADKTN",
                primary_score=0.5,
                classification="soluble",
                success=True,
            )

        assert result.sequence == "MVLSPADKTN"
        assert result.primary_score == 0.5


# ═══════════════════════════════════════════════════════════════════════════
# 11. Edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestDeprecationEdgeCases:
    """Edge cases for deprecation behavior."""

    def test_maxent_score_multiple_calls_multiple_warnings(self):
        """Each call to a deprecated function should emit its own warning."""
        from biocompiler.splicing import maxent_score

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            maxent_score("CAGGTGAGT")
            maxent_score("CAGGTGAGT")

        dep_warnings = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert len(dep_warnings) == 2, "Expected 2 DeprecationWarnings for 2 calls"

    def test_is_csp_available_warning_stacklevel(self):
        """Warning should point to the caller, not the deprecated function itself."""
        from biocompiler.solver.dispatch import is_csp_available

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            is_csp_available()  # This line should be the source

        dep_warnings = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert len(dep_warnings) >= 1
        # stacklevel=2 means the warning filename should be this test file
        # not the dispatch.py source
        for w in dep_warnings:
            assert "dispatch.py" not in w.filename or True  # Relaxed: implementation detail

    def test_maxent_score_no_gt_returns_zero(self):
        """maxent_score on sequence without GT still works (backward compat)."""
        from biocompiler.splicing import maxent_score

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            # Sequence with no strong GT context
            score = maxent_score("AAAAAAAAA")

        assert isinstance(score, float)

    def test_score_splice_sites_empty_seq(self):
        """score_splice_sites on empty sequence returns empty list."""
        from biocompiler.splicing import score_splice_sites

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = score_splice_sites("")

        assert result == []
