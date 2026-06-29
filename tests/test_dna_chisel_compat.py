"""Tests for the dna_chisel_compat module.

Covers:
1. Import/availability check — is_dna_chisel_available()
2. Compatibility functions — _build_initial_sequence, _build_dna_chisel_spec,
   _count_restriction_sites
3. Type correctness of return values — ChiselResult, ComparisonResult,
   ComparativeBenchmarkReport
4. Error handling for missing DNA Chisel
5. Winner computation — _compute_winners
6. Comparative summary — _compute_comparative_summary
7. Report formatting — format_comparative_report_text
8. Named constants

Tests are designed to work whether or not DNA Chisel is installed.
When DNA Chisel is absent, functions that depend on it return graceful
fallbacks rather than raising, which we verify explicitly.
"""

from __future__ import annotations

import importlib
import sys
from dataclasses import fields
from unittest.mock import MagicMock, patch

import pytest

from biocompiler.infrastructure.dna_chisel_compat import (
    ChiselResult,
    ComparisonResult,
    ComparativeBenchmarkReport,
    GC_ENFORCEMENT_WINDOW,
    MAX_RESTRICTION_ENZYMES,
    CAI_COMPARISON_EPSILON,
    is_dna_chisel_available,
    optimize_with_dna_chisel,
    compare_optimizers,
    run_comparative_benchmark,
    format_comparative_report_text,
    _build_initial_sequence,
    _build_dna_chisel_spec,
    _count_restriction_sites,
    _compute_winners,
    _compute_comparative_summary,
    _AA_ONE_TO_THREE,
    _ORGANISM_MAP,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#: A short protein for quick tests (Green Fluorescent Protein fragment)
SHORT_PROTEIN = "MVSKGE"

#: A slightly longer protein for integration-style tests
MEDIUM_PROTEIN = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE"


# ---------------------------------------------------------------------------
# 1. Import / Availability Check
# ---------------------------------------------------------------------------

class TestAvailability:
    """Tests for is_dna_chisel_available()."""

    def test_returns_bool(self):
        """is_dna_chisel_available() returns a bool."""
        result = is_dna_chisel_available()
        assert isinstance(result, bool)

    def test_consistent_with_module_flag(self):
        """Return value matches the module-level _DNA_CHISEL_AVAILABLE flag."""
        from biocompiler.infrastructure.dna_chisel_compat import _DNA_CHISEL_AVAILABLE
        assert is_dna_chisel_available() == _DNA_CHISEL_AVAILABLE

    def test_result_does_not_raise(self):
        """Calling the function never raises, regardless of install state."""
        # Should work cleanly both with and without DNA Chisel installed
        assert is_dna_chisel_available() in (True, False)


# ---------------------------------------------------------------------------
# 2. Named Constants
# ---------------------------------------------------------------------------

class TestConstants:
    """Tests for module-level named constants."""

    def test_gc_enforcement_window_positive(self):
        """GC_ENFORCEMENT_WINDOW should be a positive integer."""
        assert isinstance(GC_ENFORCEMENT_WINDOW, int)
        assert GC_ENFORCEMENT_WINDOW > 0

    def test_max_restriction_enzymes_positive(self):
        """MAX_RESTRICTION_ENZYMES should be a positive integer."""
        assert isinstance(MAX_RESTRICTION_ENZYMES, int)
        assert MAX_RESTRICTION_ENZYMES > 0

    def test_cai_comparison_epsilon_positive(self):
        """CAI_COMPARISON_EPSILON should be a positive float."""
        assert isinstance(CAI_COMPARISON_EPSILON, float)
        assert CAI_COMPARISON_EPSILON > 0.0

    def test_cai_comparison_epsilon_small(self):
        """CAI_COMPARISON_EPSILON should be small (< 0.1)."""
        assert CAI_COMPARISON_EPSILON < 0.1

    def test_aa_one_to_three_has_20_entries(self):
        """_AA_ONE_TO_THREE should map all 20 standard amino acids."""
        assert len(_AA_ONE_TO_THREE) == 20

    def test_aa_one_to_three_keys_are_single_letters(self):
        """All keys in _AA_ONE_TO_THREE are single uppercase letters."""
        for key in _AA_ONE_TO_THREE:
            assert len(key) == 1
            assert key.isalpha()
            assert key.isupper()

    def test_aa_one_to_three_values_are_three_letter(self):
        """All values in _AA_ONE_TO_THREE are 3-letter codes starting uppercase."""
        for val in _AA_ONE_TO_THREE.values():
            assert len(val) == 3
            assert val[0].isupper()

    def test_organism_map_has_expected_entries(self):
        """_ORGANISM_MAP has at least 5 organism mappings."""
        assert len(_ORGANISM_MAP) >= 5

    def test_organism_map_contains_human(self):
        """_ORGANISM_MAP includes Homo_sapiens."""
        assert "Homo_sapiens" in _ORGANISM_MAP

    def test_organism_map_contains_ecoli(self):
        """_ORGANISM_MAP includes Escherichia_coli."""
        assert "Escherichia_coli" in _ORGANISM_MAP


# ---------------------------------------------------------------------------
# 3. _build_initial_sequence
# ---------------------------------------------------------------------------

class TestBuildInitialSequence:
    """Tests for the _build_initial_sequence helper."""

    def test_returns_string(self):
        """_build_initial_sequence returns a string."""
        result = _build_initial_sequence(SHORT_PROTEIN)
        assert isinstance(result, str)

    def test_length_is_three_times_protein(self):
        """DNA sequence length is 3× the protein length (1 codon per AA)."""
        result = _build_initial_sequence(SHORT_PROTEIN)
        assert len(result) == len(SHORT_PROTEIN) * 3

    def test_all_bases_valid(self):
        """Returned sequence contains only A, C, G, T (no NNN placeholders for known AAs)."""
        result = _build_initial_sequence(SHORT_PROTEIN)
        for base in result:
            assert base in "ACGT", f"Unexpected base: {base!r}"

    def test_default_organism(self):
        """Default organism (Homo_sapiens) produces a valid sequence."""
        result = _build_initial_sequence("M")
        assert len(result) == 3
        assert all(b in "ACGT" for b in result)

    def test_ecoli_organism(self):
        """Specifying E. coli organism works."""
        result = _build_initial_sequence("M", organism="Escherichia_coli")
        assert len(result) == 3
        assert all(b in "ACGT" for b in result)

    def test_unsupported_organism_falls_back(self):
        """Unsupported organism falls back to Homo_sapiens (no crash)."""
        result = _build_initial_sequence("M", organism="Unknown_organism")
        assert len(result) == 3
        assert all(b in "ACGT" for b in result)

    def test_unknown_amino_acid_gets_nnn(self):
        """Unknown amino acid characters get NNN placeholder."""
        result = _build_initial_sequence("X")
        assert result == "NNN"

    def test_mixed_known_unknown_amino_acids(self):
        """Protein with both known and unknown AAs handles correctly."""
        result = _build_initial_sequence("MX")
        assert len(result) == 6
        # First codon is for M (ATG or similar), second is NNN
        assert result[3:] == "NNN"
        assert all(b in "ACGT" for b in result[:3])

    def test_single_methionine(self):
        """Single M should give ATG (the only codon for Met)."""
        result = _build_initial_sequence("M")
        assert result == "ATG"

    def test_empty_protein(self):
        """Empty protein returns empty sequence."""
        result = _build_initial_sequence("")
        assert result == ""

    def test_all_standard_amino_acids(self):
        """All 20 standard AAs produce valid codons without NNN."""
        all_aas = "ACDEFGHIKLMNPQRSTVWY"
        result = _build_initial_sequence(all_aas)
        assert len(result) == len(all_aas) * 3
        assert "N" not in result  # No NNN placeholders for standard AAs


# ---------------------------------------------------------------------------
# 4. _count_restriction_sites
# ---------------------------------------------------------------------------

class TestCountRestrictionSites:
    """Tests for the _count_restriction_sites helper."""

    def test_returns_int(self):
        """_count_restriction_sites returns an int."""
        result = _count_restriction_sites("ATCGATCGATCG")
        assert isinstance(result, int)

    def test_empty_sequence_zero_sites(self):
        """Empty sequence has zero restriction sites."""
        result = _count_restriction_sites("")
        assert result == 0

    def test_known_ecori_site(self):
        """EcoRI site (GAATTC) is counted in a sequence containing it."""
        seq = "AAAAAAGAATTCAAAAAA"
        result = _count_restriction_sites(seq, restriction_enzymes=["EcoRI"])
        assert result >= 1

    def test_no_sites_returns_zero(self):
        """Sequence with no restriction sites returns 0."""
        seq = "AAAAAAAAAAAA"
        result = _count_restriction_sites(seq, restriction_enzymes=["EcoRI"])
        assert result == 0

    def test_multiple_same_site(self):
        """Multiple occurrences of the same site are all counted."""
        seq = "GAATTCGAATTC"  # Two EcoRI sites
        result = _count_restriction_sites(seq, restriction_enzymes=["EcoRI"])
        assert result == 2

    def test_palindromic_site_not_double_counted(self):
        """Palindromic sites (where RC == original) are not double-counted."""
        # EcoRI: GAATTC, RC = GAATTC (palindrome)
        seq = "GAATTC"
        result = _count_restriction_sites(seq, restriction_enzymes=["EcoRI"])
        # Should be exactly 1, not 2
        assert result == 1

    def test_none_enzymes_uses_all(self):
        """None for restriction_enzymes defaults to all enzymes in RESTRICTION_ENZYMES."""
        result = _count_restriction_sites("GAATTC", restriction_enzymes=None)
        # Should count the EcoRI site at minimum
        assert result >= 1

    def test_empty_enzymes_list_uses_default(self):
        """Empty list for restriction_enzymes is treated as falsy and defaults to all enzymes."""
        # The function uses `if not restriction_enzymes:` which treats [] as falsy,
        # falling back to all RESTRICTION_ENZYMES. So GAATTC will find EcoRI.
        result = _count_restriction_sites("GAATTC", restriction_enzymes=[])
        assert result >= 1  # Falls back to default enzyme list

    def test_unknown_enzyme_name_ignored(self):
        """Unknown enzyme name is silently ignored."""
        result = _count_restriction_sites("GAATTC", restriction_enzymes=["NonExistentEnzyme"])
        assert result == 0

    def test_case_insensitive(self):
        """Sequence matching is case-insensitive."""
        result_upper = _count_restriction_sites("GAATTC", restriction_enzymes=["EcoRI"])
        result_lower = _count_restriction_sites("gaattc", restriction_enzymes=["EcoRI"])
        assert result_upper == result_lower

    def test_iupac_ambiguity_sites_skipped(self):
        """Enzymes with IUPAC ambiguity codes (e.g. SfiI with Ns) are skipped."""
        # SfiI site is GGCCNNNNNGGCC — contains N, should be skipped
        seq = "GGCCAAAAGGCC"  # Would match SfiI if N were expanded
        result = _count_restriction_sites(seq, restriction_enzymes=["SfiI"])
        assert result == 0


# ---------------------------------------------------------------------------
# 5. ChiselResult Dataclass
# ---------------------------------------------------------------------------

class TestChiselResult:
    """Tests for the ChiselResult dataclass."""

    def test_fields_exist(self):
        """ChiselResult has all expected fields."""
        expected_fields = {
            "sequence", "protein", "cai", "gc_content",
            "restriction_site_count", "execution_time_s", "success", "error",
        }
        actual_fields = {f.name for f in fields(ChiselResult)}
        assert actual_fields == expected_fields

    def test_construction(self):
        """ChiselResult can be constructed with all fields."""
        result = ChiselResult(
            sequence="ATGCGT",
            protein="MR",
            cai=0.8,
            gc_content=0.5,
            restriction_site_count=0,
            execution_time_s=0.1,
            success=True,
            error=None,
        )
        assert result.sequence == "ATGCGT"
        assert result.protein == "MR"
        assert result.cai == 0.8
        assert result.success is True
        assert result.error is None

    def test_error_default_none(self):
        """error field defaults to None."""
        result = ChiselResult(
            sequence="", protein="", cai=0.0, gc_content=0.0,
            restriction_site_count=0, execution_time_s=0.0, success=False,
        )
        assert result.error is None

    def test_failure_result_typical_fields(self):
        """A typical failure ChiselResult has correct types."""
        result = ChiselResult(
            sequence="",
            protein="MR",
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=0.0,
            success=False,
            error="DNA Chisel not installed",
        )
        assert isinstance(result.sequence, str)
        assert isinstance(result.cai, float)
        assert isinstance(result.gc_content, float)
        assert isinstance(result.restriction_site_count, int)
        assert isinstance(result.execution_time_s, float)
        assert isinstance(result.success, bool)
        assert isinstance(result.error, (str, type(None)))


# ---------------------------------------------------------------------------
# 6. ComparisonResult Dataclass
# ---------------------------------------------------------------------------

class TestComparisonResult:
    """Tests for the ComparisonResult dataclass."""

    def test_fields_exist(self):
        """ComparisonResult has all expected fields."""
        expected_fields = {
            "protein", "organism", "biocompiler", "dna_chisel",
            "dna_chisel_available", "winner",
        }
        actual_fields = {f.name for f in fields(ComparisonResult)}
        assert actual_fields == expected_fields

    def test_construction(self):
        """ComparisonResult can be constructed with all fields."""
        result = ComparisonResult(
            protein="MR",
            organism="Homo_sapiens",
            biocompiler={"cai": 0.8, "success": True},
            dna_chisel=None,
            dna_chisel_available=False,
            winner={"overall": "biocompiler"},
        )
        assert result.protein == "MR"
        assert result.dna_chisel is None
        assert result.dna_chisel_available is False

    def test_dna_chisel_can_be_none(self):
        """dna_chisel field can be None (when DNA Chisel is unavailable)."""
        result = ComparisonResult(
            protein="MR",
            organism="Homo_sapiens",
            biocompiler={},
            dna_chisel=None,
            dna_chisel_available=False,
            winner={},
        )
        assert result.dna_chisel is None


# ---------------------------------------------------------------------------
# 7. ComparativeBenchmarkReport Dataclass
# ---------------------------------------------------------------------------

class TestComparativeBenchmarkReport:
    """Tests for the ComparativeBenchmarkReport dataclass."""

    def test_fields_exist(self):
        """ComparativeBenchmarkReport has all expected fields."""
        expected_fields = {
            "timestamp", "dna_chisel_available", "gene_results", "summary",
        }
        actual_fields = {f.name for f in fields(ComparativeBenchmarkReport)}
        assert actual_fields == expected_fields

    def test_total_genes_property(self):
        """total_genes property returns length of gene_results."""
        report = ComparativeBenchmarkReport(
            timestamp="2025-01-01T00:00:00",
            dna_chisel_available=False,
            gene_results=[{"gene": "g1"}, {"gene": "g2"}, {"gene": "g3"}],
            summary={},
        )
        assert report.total_genes == 3

    def test_total_genes_empty(self):
        """total_genes is 0 when gene_results is empty."""
        report = ComparativeBenchmarkReport(
            timestamp="2025-01-01T00:00:00",
            dna_chisel_available=False,
            gene_results=[],
            summary={},
        )
        assert report.total_genes == 0

    def test_default_gene_results_empty_list(self):
        """gene_results defaults to an empty list."""
        report = ComparativeBenchmarkReport(
            timestamp="2025-01-01T00:00:00",
            dna_chisel_available=False,
        )
        assert report.gene_results == []

    def test_default_summary_empty_dict(self):
        """summary defaults to an empty dict."""
        report = ComparativeBenchmarkReport(
            timestamp="2025-01-01T00:00:00",
            dna_chisel_available=False,
        )
        assert report.summary == {}


# ---------------------------------------------------------------------------
# 8. optimize_with_dna_chisel — Error Handling for Missing DNA Chisel
# ---------------------------------------------------------------------------

class TestOptimizeWithDnaChiselWithoutChisel:
    """Tests for optimize_with_dna_chisel when DNA Chisel is not installed."""

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_returns_chisel_result(self):
        """Returns a ChiselResult even when DNA Chisel is unavailable."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result, ChiselResult)

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_success_is_false(self):
        """success field is False when DNA Chisel is unavailable."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert result.success is False

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_error_message_mentions_dna_chisel(self):
        """error field mentions DNA Chisel not being installed."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert result.error is not None
        assert "DNA Chisel" in result.error or "dnachisel" in result.error.lower()

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_sequence_is_empty(self):
        """sequence field is empty string when unavailable."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert result.sequence == ""

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_metrics_are_zero(self):
        """cai and gc_content are 0.0 when unavailable."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert result.cai == 0.0
        assert result.gc_content == 0.0
        assert result.restriction_site_count == 0

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_protein_is_preserved(self):
        """protein field matches the input protein."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert result.protein == SHORT_PROTEIN

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="These tests verify behavior when DNA Chisel is NOT installed"
    )
    def test_no_exception_raised(self):
        """optimize_with_dna_chisel does not raise — it degrades gracefully."""
        # Should never raise, even with unusual inputs
        optimize_with_dna_chisel("")
        optimize_with_dna_chisel("M")


# ---------------------------------------------------------------------------
# 9. optimize_with_dna_chisel — With DNA Chisel (mocked)
# ---------------------------------------------------------------------------

class TestOptimizeWithDnaChiselMocked:
    """Tests for optimize_with_dna_chisel with DNA Chisel mocked as available."""

    def _make_mock_chisel_available(self):
        """Return a dict of patches that simulate DNA Chisel being available."""
        mock_dna_optimization_problem = MagicMock()
        mock_sequence = MagicMock()
        mock_sequence.__str__ = MagicMock(return_value="ATGGTTTCAAAGGGTGAAGAG")
        mock_dna_optimization_problem.return_value.sequence = mock_sequence
        mock_dna_optimization_problem.return_value.resolve_constraints = MagicMock()
        return mock_dna_optimization_problem

    @patch("biocompiler.infrastructure.dna_chisel_compat._DNA_CHISEL_AVAILABLE", True)
    @patch("biocompiler.infrastructure.dna_chisel_compat.compute_cai", return_value=0.75)
    @patch("biocompiler.infrastructure.dna_chisel_compat.gc_content", return_value=0.55)
    @patch("biocompiler.infrastructure.dna_chisel_compat._count_restriction_sites", return_value=0)
    @patch("biocompiler.infrastructure.dna_chisel_compat._build_dna_chisel_spec", return_value=[])
    @patch("biocompiler.infrastructure.dna_chisel_compat._build_initial_sequence", return_value="ATGGTTTCAAAGGGTGAAGAG")
    def test_success_when_available(
        self, mock_build, mock_spec, mock_count, mock_gc, mock_cai,
    ):
        """When DNA Chisel is available, optimize returns success=True."""
        # We also need to mock the DnaOptimizationProblem class. Since it
        # may not exist as a module attribute when DNA Chisel is not installed,
        # we inject it into the module's namespace.
        import biocompiler.infrastructure.dna_chisel_compat as mod
        mock_problem_cls = MagicMock()
        mock_instance = MagicMock()
        mock_instance.sequence = "ATGGTTTCAAAGGGTGAAGAG"
        mock_instance.resolve_constraints = MagicMock()
        mock_problem_cls.return_value = mock_instance

        original = getattr(mod, "DnaOptimizationProblem", None)
        mod.DnaOptimizationProblem = mock_problem_cls
        try:
            result = optimize_with_dna_chisel("MSKG")
            assert isinstance(result, ChiselResult)
            assert result.success is True
        finally:
            if original is None:
                delattr(mod, "DnaOptimizationProblem")
            else:
                mod.DnaOptimizationProblem = original

    @patch("biocompiler.infrastructure.dna_chisel_compat._DNA_CHISEL_AVAILABLE", True)
    @patch("biocompiler.infrastructure.dna_chisel_compat.compute_cai", return_value=0.75)
    @patch("biocompiler.infrastructure.dna_chisel_compat.gc_content", return_value=0.55)
    @patch("biocompiler.infrastructure.dna_chisel_compat._count_restriction_sites", return_value=2)
    @patch("biocompiler.infrastructure.dna_chisel_compat._build_dna_chisel_spec", return_value=[])
    @patch("biocompiler.infrastructure.dna_chisel_compat._build_initial_sequence", return_value="ATGGTTTCAAAGGGTGAAGAG")
    def test_metrics_populated(
        self, mock_build, mock_spec, mock_count, mock_gc, mock_cai,
    ):
        """When successful, CAI and GC metrics are populated correctly."""
        import biocompiler.infrastructure.dna_chisel_compat as mod
        mock_problem_cls = MagicMock()
        mock_instance = MagicMock()
        mock_instance.sequence = "ATGGTTTCAAAGGGTGAAGAG"
        mock_instance.resolve_constraints = MagicMock()
        mock_problem_cls.return_value = mock_instance

        original = getattr(mod, "DnaOptimizationProblem", None)
        mod.DnaOptimizationProblem = mock_problem_cls
        try:
            result = optimize_with_dna_chisel("MSKG")
            assert result.cai == 0.75
            assert result.gc_content == 0.55
            assert result.restriction_site_count == 2
        finally:
            if original is None:
                delattr(mod, "DnaOptimizationProblem")
            else:
                mod.DnaOptimizationProblem = original

    @patch("biocompiler.infrastructure.dna_chisel_compat._DNA_CHISEL_AVAILABLE", True)
    @patch("biocompiler.infrastructure.dna_chisel_compat._build_initial_sequence", return_value="ATG")
    @patch("biocompiler.infrastructure.dna_chisel_compat._build_dna_chisel_spec")
    def test_exception_returns_failure_result(self, mock_spec, mock_build):
        """If DNA Chisel raises during optimization, a failure result is returned."""
        mock_spec.side_effect = RuntimeError("Chisel exploded")
        result = optimize_with_dna_chisel("M")
        assert isinstance(result, ChiselResult)
        assert result.success is False
        assert result.error is not None
        assert "Chisel exploded" in result.error


# ---------------------------------------------------------------------------
# 10. _build_dna_chisel_spec — Error Handling for Missing DNA Chisel
# ---------------------------------------------------------------------------

class TestBuildDnaChiselSpec:
    """Tests for _build_dna_chisel_spec."""

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="Tests ImportError when DNA Chisel is NOT installed"
    )
    def test_raises_import_error_when_unavailable(self):
        """Raises ImportError when DNA Chisel is not installed."""
        with pytest.raises(ImportError, match="DNA Chisel"):
            _build_dna_chisel_spec("MR")

    @pytest.mark.skipif(
        not is_dna_chisel_available(),
        reason="Tests spec building when DNA Chisel IS installed"
    )
    def test_returns_list_when_available(self):
        """Returns a list of constraints when DNA Chisel is available."""
        constraints = _build_dna_chisel_spec("MR")
        assert isinstance(constraints, list)
        assert len(constraints) >= 1  # At least EnforceTranslation

    @pytest.mark.skipif(
        not is_dna_chisel_available(),
        reason="Tests spec building when DNA Chisel IS installed"
    )
    def test_restriction_enzymes_add_constraints(self):
        """Specifying restriction enzymes adds AvoidPattern constraints."""
        constraints_no_enz = _build_dna_chisel_spec("MR")
        constraints_with_enz = _build_dna_chisel_spec("MR", restriction_enzymes=["EcoRI"])
        assert len(constraints_with_enz) > len(constraints_no_enz)

    @pytest.mark.skipif(
        not is_dna_chisel_available(),
        reason="Tests spec building when DNA Chisel IS installed"
    )
    def test_gc_bounds_passed(self):
        """GC bounds are passed to EnforceGCContent."""
        # Just verify it does not crash with various GC bounds
        constraints = _build_dna_chisel_spec("MR", gc_lo=0.4, gc_hi=0.6)
        assert isinstance(constraints, list)

    @pytest.mark.skipif(
        not is_dna_chisel_available(),
        reason="Tests spec building when DNA Chisel IS installed"
    )
    def test_unknown_enzyme_ignored(self):
        """Unknown enzyme name is silently ignored (no crash)."""
        constraints = _build_dna_chisel_spec("MR", restriction_enzymes=["FakeEnzyme"])
        # Should still have EnforceTranslation + EnforceGCContent
        assert len(constraints) >= 2


# ---------------------------------------------------------------------------
# 11. _compute_winners
# ---------------------------------------------------------------------------

class TestComputeWinners:
    """Tests for the _compute_winners helper."""

    def test_no_dc_result_biocompiler_wins(self):
        """When dc is None, BioCompiler wins by default."""
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, None, 0.3, 0.7)
        assert "biocompiler" in result["overall"]

    def test_dc_failed_biocompiler_wins(self):
        """When dc has success=False, BioCompiler wins."""
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": False, "cai": 0.0, "gc_content": 0.0,
              "restriction_site_count": 0, "execution_time_s": 0.0}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert "biocompiler" in result["overall"]

    def test_bc_failed_dc_wins(self):
        """When bc has success=False, DNA Chisel wins."""
        bc = {"success": False, "cai": 0.0, "gc_content": 0.0,
              "restriction_site_count": 0, "execution_time_s": 0.0}
        dc = {"success": True, "cai": 0.7, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert "dna_chisel" in result["overall"]

    def test_both_succeed_cai_winner(self):
        """When both succeed, higher CAI wins the CAI metric."""
        bc = {"success": True, "cai": 0.9, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": True, "cai": 0.5, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert result["metrics"]["cai"]["winner"] == "biocompiler"

    def test_both_succeed_dc_higher_cai(self):
        """DNA Chisel wins CAI metric when its score is higher."""
        bc = {"success": True, "cai": 0.5, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": True, "cai": 0.9, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert result["metrics"]["cai"]["winner"] == "dna_chisel"

    def test_cai_within_epsilon_is_tie(self):
        """CAI difference within epsilon counts as a tie."""
        bc = {"success": True, "cai": 0.8000, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": True, "cai": 0.8000 + CAI_COMPARISON_EPSILON * 0.5,
              "gc_content": 0.5, "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert result["metrics"]["cai"]["winner"] == "tie"

    def test_restriction_sites_fewer_wins(self):
        """Fewer restriction sites wins that metric."""
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 3, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert result["metrics"]["restriction_site_count"]["winner"] == "biocompiler"

    def test_gc_closer_to_midpoint_wins(self):
        """GC closer to midpoint of gc_lo and gc_hi wins."""
        gc_lo, gc_hi = 0.3, 0.7
        midpoint = 0.5
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": True, "cai": 0.8, "gc_content": 0.35,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, gc_lo, gc_hi)
        assert result["metrics"]["gc_content"]["winner"] == "biocompiler"

    def test_execution_time_faster_wins(self):
        """Faster execution time wins (with 10% margin)."""
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.05}
        dc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 1.0}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        assert result["metrics"]["execution_time_s"]["winner"] == "biocompiler"

    def test_result_has_all_metrics(self):
        """Result dict has entries for all four metrics."""
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        dc = {"success": True, "cai": 0.7, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, dc, 0.3, 0.7)
        for metric in ("cai", "gc_content", "restriction_site_count", "execution_time_s"):
            assert metric in result["metrics"]

    def test_result_has_overall(self):
        """Result dict has an 'overall' key."""
        bc = {"success": True, "cai": 0.8, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, None, 0.3, 0.7)
        assert "overall" in result


# ---------------------------------------------------------------------------
# 12. _compute_comparative_summary
# ---------------------------------------------------------------------------

class TestComputeComparativeSummary:
    """Tests for the _compute_comparative_summary helper."""

    def test_empty_gene_results(self):
        """Empty gene_results produces a valid summary with zeros."""
        summary = _compute_comparative_summary([])
        assert summary["total_genes"] == 0
        assert summary["genes_with_errors"] == 0

    def test_summary_has_required_keys(self):
        """Summary dict has all required top-level keys."""
        summary = _compute_comparative_summary([])
        expected_keys = {
            "total_genes", "genes_with_errors", "metric_wins",
            "overall_wins", "avg_cai", "avg_gc",
            "avg_restriction_sites", "avg_execution_time_s",
        }
        for key in expected_keys:
            assert key in summary, f"Missing key: {key}"

    def test_metric_wins_has_all_metrics(self):
        """metric_wins dict has entries for all four metrics."""
        summary = _compute_comparative_summary([])
        for metric in ("cai", "gc_content", "restriction_site_count", "execution_time_s"):
            assert metric in summary["metric_wins"]
            assert "biocompiler" in summary["metric_wins"][metric]
            assert "dna_chisel" in summary["metric_wins"][metric]

    def test_error_genes_counted(self):
        """Genes with errors are counted in genes_with_errors."""
        gene_results = [
            {"error": "some error"},
            {"winner": {}, "biocompiler": {"success": True}, "dna_chisel": None},
        ]
        summary = _compute_comparative_summary(gene_results)
        assert summary["genes_with_errors"] == 1

    def test_averages_computed_from_successful_results(self):
        """Averages are computed from successful optimization results."""
        gene_results = [
            {
                "winner": {"overall": "biocompiler", "metrics": {
                    "cai": {"winner": "biocompiler"},
                    "gc_content": {"winner": "tie"},
                    "restriction_site_count": {"winner": "biocompiler"},
                    "execution_time_s": {"winner": "biocompiler"},
                }},
                "biocompiler": {
                    "success": True, "cai": 0.8, "gc_content": 0.5,
                    "restriction_site_count": 1, "execution_time_s": 0.2,
                },
                "dna_chisel": None,
            },
        ]
        summary = _compute_comparative_summary(gene_results)
        assert summary["avg_cai"]["biocompiler"] == pytest.approx(0.8, abs=0.01)
        assert summary["avg_gc"]["biocompiler"] == pytest.approx(0.5, abs=0.01)


# ---------------------------------------------------------------------------
# 13. compare_optimizers
# ---------------------------------------------------------------------------

class TestCompareOptimizers:
    """Tests for the compare_optimizers public function."""

    def test_returns_comparison_result(self):
        """compare_optimizers returns a ComparisonResult."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert isinstance(result, ComparisonResult)

    def test_protein_preserved(self):
        """Protein field matches the input."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert result.protein == SHORT_PROTEIN

    def test_organism_default(self):
        """Default organism is Homo_sapiens."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert result.organism == "Homo_sapiens"

    def test_organism_custom(self):
        """Custom organism is passed through."""
        result = compare_optimizers(SHORT_PROTEIN, organism="Escherichia_coli")
        assert result.organism == "Escherichia_coli"

    def test_biocompiler_metrics_is_dict(self):
        """biocompiler field is a dict."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert isinstance(result.biocompiler, dict)

    def test_biocompiler_has_success_key(self):
        """biocompiler metrics dict has a 'success' key."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert "success" in result.biocompiler

    def test_winner_is_dict(self):
        """winner field is a dict."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert isinstance(result.winner, dict)
        assert "overall" in result.winner

    def test_dna_chisel_available_matches_module(self):
        """dna_chisel_available matches the module-level flag."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert result.dna_chisel_available == is_dna_chisel_available()

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="Tests behavior when DNA Chisel is NOT installed"
    )
    def test_dc_none_when_unavailable(self):
        """dna_chisel field is None when DNA Chisel is not installed."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert result.dna_chisel is None

    def test_no_exception_on_short_protein(self):
        """Does not raise even with very short protein."""
        result = compare_optimizers("M")
        assert isinstance(result, ComparisonResult)

    def test_custom_gc_bounds(self):
        """Custom GC bounds are accepted without error."""
        result = compare_optimizers(SHORT_PROTEIN, gc_lo=0.4, gc_hi=0.6)
        assert isinstance(result, ComparisonResult)

    def test_custom_restriction_enzymes(self):
        """Custom restriction enzymes list is accepted."""
        result = compare_optimizers(
            SHORT_PROTEIN, restriction_enzymes=["EcoRI", "BamHI"]
        )
        assert isinstance(result, ComparisonResult)


# ---------------------------------------------------------------------------
# 14. run_comparative_benchmark
# ---------------------------------------------------------------------------

class TestRunComparativeBenchmark:
    """Tests for run_comparative_benchmark."""

    def test_returns_comparative_benchmark_report(self):
        """Returns a ComparativeBenchmarkReport."""
        # Use a small subset of genes (empty list for speed)
        report = run_comparative_benchmark(genes=[])
        assert isinstance(report, ComparativeBenchmarkReport)

    def test_empty_genes_list_falls_back_to_default(self):
        """Empty genes list is treated as falsy and falls back to REFERENCE_GENES.

        The function uses `genes or list(REFERENCE_GENES.keys())`, so an empty
        list triggers the default behavior of benchmarking all reference genes.
        """
        report = run_comparative_benchmark(genes=[])
        # Empty list is falsy, so all reference genes are used
        assert report.total_genes > 0

    def test_unknown_gene_skipped(self):
        """Unknown gene names are skipped gracefully."""
        report = run_comparative_benchmark(genes=["NONEXISTENT_GENE_XYZ"])
        # Should not crash, results should be empty (gene not found in REFERENCE_GENES)
        assert isinstance(report, ComparativeBenchmarkReport)

    def test_timestamp_is_string(self):
        """timestamp field is a string."""
        report = run_comparative_benchmark(genes=[])
        assert isinstance(report.timestamp, str)
        assert len(report.timestamp) > 0

    def test_dna_chisel_available_flag(self):
        """dna_chisel_available matches module-level flag."""
        report = run_comparative_benchmark(genes=[])
        assert report.dna_chisel_available == is_dna_chisel_available()


# ---------------------------------------------------------------------------
# 15. format_comparative_report_text
# ---------------------------------------------------------------------------

class TestFormatComparativeReportText:
    """Tests for format_comparative_report_text."""

    def _make_empty_report(self) -> ComparativeBenchmarkReport:
        """Create a minimal empty report for testing."""
        return ComparativeBenchmarkReport(
            timestamp="2025-01-01T00:00:00+00:00",
            dna_chisel_available=False,
            gene_results=[],
            summary={
                "total_genes": 0,
                "genes_with_errors": 0,
                "metric_wins": {
                    "cai": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
                    "gc_content": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
                    "restriction_site_count": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
                    "execution_time_s": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
                },
                "overall_wins": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
                "avg_cai": {"biocompiler": 0.0, "dna_chisel": 0.0},
                "avg_gc": {"biocompiler": 0.0, "dna_chisel": 0.0},
                "avg_restriction_sites": {"biocompiler": 0.0, "dna_chisel": 0.0},
                "avg_execution_time_s": {"biocompiler": 0.0, "dna_chisel": 0.0},
            },
        )

    def test_returns_string(self):
        """format_comparative_report_text returns a string."""
        report = self._make_empty_report()
        text = format_comparative_report_text(report)
        assert isinstance(text, str)

    def test_contains_header(self):
        """Output contains the benchmark header."""
        report = self._make_empty_report()
        text = format_comparative_report_text(report)
        assert "BioCompiler vs DNA Chisel" in text

    def test_contains_timestamp(self):
        """Output contains the timestamp."""
        report = self._make_empty_report()
        text = format_comparative_report_text(report)
        assert "2025-01-01" in text

    def test_shows_unavailable_note(self):
        """When DNA Chisel is unavailable, output includes a note."""
        report = self._make_empty_report()
        report.dna_chisel_available = False
        text = format_comparative_report_text(report)
        assert "not installed" in text.lower() or "NOT AVAILABLE" in text

    def test_shows_available_status(self):
        """When DNA Chisel is available, output reflects that."""
        report = self._make_empty_report()
        report.dna_chisel_available = True
        text = format_comparative_report_text(report)
        assert "Yes" in text

    def test_with_gene_results(self):
        """Output includes gene results when present."""
        report = self._make_empty_report()
        report.gene_results = [
            {
                "gene": "eGFP",
                "description": "Enhanced GFP",
                "biocompiler": {
                    "success": True, "cai": 0.85, "gc_content": 0.55,
                    "restriction_site_count": 2, "execution_time_s": 0.12,
                },
                "dna_chisel": None,
                "winner": {"overall": "biocompiler (dna_chisel_unavailable)"},
            }
        ]
        report.summary["total_genes"] = 1
        text = format_comparative_report_text(report)
        assert "eGFP" in text
        assert "0.85" in text

    def test_summary_section_present(self):
        """Output contains a Summary section."""
        report = self._make_empty_report()
        text = format_comparative_report_text(report)
        assert "Summary" in text

    def test_per_metric_wins_section(self):
        """Output contains Per-Metric Wins section."""
        report = self._make_empty_report()
        text = format_comparative_report_text(report)
        assert "Per-Metric Wins" in text

    def test_overall_wins_section(self):
        """Output contains Overall Wins section."""
        report = self._make_empty_report()
        text = format_comparative_report_text(report)
        assert "Overall Wins" in text

    def test_failed_biocompiler_shown(self):
        """Failed BioCompiler result is displayed in output."""
        report = self._make_empty_report()
        report.gene_results = [
            {
                "gene": "TestGene",
                "description": "Test",
                "biocompiler": {
                    "success": False, "error": "optimization failed",
                },
                "dna_chisel": None,
                "winner": {},
            }
        ]
        text = format_comparative_report_text(report)
        assert "FAILED" in text

    def test_failed_dna_chisel_shown(self):
        """Failed DNA Chisel result is displayed in output."""
        report = self._make_empty_report()
        report.gene_results = [
            {
                "gene": "TestGene",
                "description": "Test",
                "biocompiler": {
                    "success": True, "cai": 0.8, "gc_content": 0.5,
                    "restriction_site_count": 0, "execution_time_s": 0.1,
                },
                "dna_chisel": {
                    "success": False, "error": "chisel error",
                },
                "winner": {"overall": "biocompiler"},
            }
        ]
        report.dna_chisel_available = True
        text = format_comparative_report_text(report)
        assert "FAILED" in text


# ---------------------------------------------------------------------------
# 16. Integration — Full workflow without DNA Chisel
# ---------------------------------------------------------------------------

class TestIntegrationWithoutChisel:
    """Integration tests that run the full workflow when DNA Chisel is absent."""

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="Tests the fallback path when DNA Chisel is NOT installed"
    )
    def test_optimize_returns_failure_result(self):
        """optimize_with_dna_chisel returns a failure result without DNA Chisel."""
        result = optimize_with_dna_chisel(MEDIUM_PROTEIN)
        assert isinstance(result, ChiselResult)
        assert result.success is False
        assert result.error is not None

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="Tests the fallback path when DNA Chisel is NOT installed"
    )
    def test_compare_returns_valid_result(self):
        """compare_optimizers returns a valid ComparisonResult without DNA Chisel."""
        result = compare_optimizers(MEDIUM_PROTEIN)
        assert isinstance(result, ComparisonResult)
        assert result.dna_chisel is None
        assert result.dna_chisel_available is False
        # BioCompiler should still produce a result
        assert isinstance(result.biocompiler, dict)

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="Tests the fallback path when DNA Chisel is NOT installed"
    )
    def test_benchmark_report_valid(self):
        """run_comparative_benchmark works without DNA Chisel (empty gene list)."""
        report = run_comparative_benchmark(genes=[])
        assert isinstance(report, ComparativeBenchmarkReport)
        assert report.dna_chisel_available is False

    @pytest.mark.skipif(
        is_dna_chisel_available(),
        reason="Tests the fallback path when DNA Chisel is NOT installed"
    )
    def test_format_report_works(self):
        """format_comparative_report_text works without DNA Chisel."""
        report = run_comparative_benchmark(genes=[])
        text = format_comparative_report_text(report)
        assert "not installed" in text.lower() or "NOT AVAILABLE" in text


# ---------------------------------------------------------------------------
# 17. Type Correctness — Return Value Types
# ---------------------------------------------------------------------------

class TestTypeCorrectness:
    """Verify that return values have the correct types."""

    def test_chisel_result_cai_is_float(self):
        """ChiselResult.cai is always a float."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result.cai, float)

    def test_chisel_result_gc_is_float(self):
        """ChiselResult.gc_content is always a float."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result.gc_content, float)

    def test_chisel_result_rs_count_is_int(self):
        """ChiselResult.restriction_site_count is always an int."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result.restriction_site_count, int)

    def test_chisel_result_time_is_float(self):
        """ChiselResult.execution_time_s is always a float."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result.execution_time_s, float)

    def test_chisel_result_success_is_bool(self):
        """ChiselResult.success is always a bool."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result.success, bool)

    def test_chisel_result_sequence_is_str(self):
        """ChiselResult.sequence is always a string."""
        result = optimize_with_dna_chisel(SHORT_PROTEIN)
        assert isinstance(result.sequence, str)

    def test_comparison_result_biocompiler_is_dict(self):
        """ComparisonResult.biocompiler is always a dict."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert isinstance(result.biocompiler, dict)

    def test_comparison_result_winner_is_dict(self):
        """ComparisonResult.winner is always a dict."""
        result = compare_optimizers(SHORT_PROTEIN)
        assert isinstance(result.winner, dict)

    def test_compute_winners_returns_dict(self):
        """_compute_winners always returns a dict."""
        bc = {"success": True, "cai": 0.5, "gc_content": 0.5,
              "restriction_site_count": 0, "execution_time_s": 0.1}
        result = _compute_winners(bc, None, 0.3, 0.7)
        assert isinstance(result, dict)

    def test_compute_summary_returns_dict(self):
        """_compute_comparative_summary always returns a dict."""
        result = _compute_comparative_summary([])
        assert isinstance(result, dict)

    def test_count_restriction_sites_returns_int(self):
        """_count_restriction_sites always returns an int."""
        result = _count_restriction_sites("ATCG")
        assert isinstance(result, int)

    def test_build_initial_sequence_returns_str(self):
        """_build_initial_sequence always returns a string."""
        result = _build_initial_sequence("MR")
        assert isinstance(result, str)

    def test_format_report_returns_str(self):
        """format_comparative_report_text always returns a string."""
        report = ComparativeBenchmarkReport(
            timestamp="2025-01-01",
            dna_chisel_available=False,
        )
        result = format_comparative_report_text(report)
        assert isinstance(result, str)
