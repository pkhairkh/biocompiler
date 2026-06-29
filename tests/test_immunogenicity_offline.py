"""Test offline immunogenicity prediction — predict_immunogenicity, scan_peptides.

Covers:
  - predict_immunogenicity() works offline (precomputed + PSSM fallback)
  - scan_peptides() generates correct overlapping peptide windows
  - PSSM fallback produces reasonable IC50 values and binding classifications
  - MHC-II alleles with 15-mer peptides are handled correctly
"""

import pytest

from biocompiler.immunogenicity.core import (
    predict_immunogenicity,
    scan_peptides,
    PeptideResult,
    ImmunogenicityPrediction,
    PRECOMPUTED_BINDERS,
    classify_binding,
    binding_score_to_ic50,
    score_peptide_pssm,
    MHC_I_PSSM,
    MHC_II_PSSM,
    MHC_II_CORE_LENGTH,
    DEFAULT_MHC_I_ALLELES,
    DEFAULT_MHC_II_ALLELES,
)

from biocompiler.immunogenicity.predicates import (
    evaluate_low_immunogenicity,
    evaluate_no_strong_t_cell_epitope,
    evaluate_no_dominant_b_cell_epitope,
    evaluate_population_coverage_safe,
)

from biocompiler.shared.types import Verdict, TypeCheckResult


# ── Test proteins (all standard amino acids only) ────────────

HBB_HUMAN = (
    "MVHLTPEEKSAVTALWGKVNVADIVGHALSDLHAKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)

# Known HLA-A*02:01 binder from influenza M1 protein
KNOWN_BINDER = "GILGFVFTL"

# A short protein for window tests — standard AAs only
SHORT_PROTEIN = "ACDEFGHIKLMNPQRSTVWY"


# ═══════════════════════════════════════════════════════════════
# TestPredictImmunogenicity
# ═══════════════════════════════════════════════════════════════

class TestPredictImmunogenicity:
    """Tests for predict_immunogenicity() offline API."""

    def test_precomputed_lookup_returns_high_confidence(self):
        """Known binder in PRECOMPUTED_BINDERS should use precomputed_lookup method."""
        result = predict_immunogenicity("HLA-A*02:01", "GILGFVFTL")
        assert isinstance(result, ImmunogenicityPrediction)
        assert result.method == "precomputed_lookup"
        assert result.confidence == "high"
        assert result.ic50_nm < 50.0  # Known strong binder
        assert result.binding_class == "strong_binder"

    def test_precomputed_lookup_correct_binding_class(self):
        """All precomputed entries should have correct binding classification."""
        for (allele, peptide), ic50 in PRECOMPUTED_BINDERS.items():
            result = predict_immunogenicity(allele, peptide)
            expected_class = classify_binding(ic50)
            assert result.binding_class == expected_class, (
                f"({allele}, {peptide}): ic50={ic50}, "
                f"expected {expected_class}, got {result.binding_class}"
            )

    def test_pssm_fallback_for_mhc_i(self):
        """Unknown peptide for a known MHC-I allele uses PSSM fallback."""
        result = predict_immunogenicity("HLA-A*02:01", "AAAAAAAAA")
        assert isinstance(result, ImmunogenicityPrediction)
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert result.confidence in ("low", "medium")
        assert result.ic50_nm > 0

    def test_pssm_fallback_for_mhc_ii(self):
        """MHC-II allele with a 15-mer uses PSSM core scanning."""
        result = predict_immunogenicity("HLA-DRB1*01:01", "PKYVKQNTLKLATLV")
        assert isinstance(result, ImmunogenicityPrediction)
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert result.confidence in ("low", "medium")
        assert result.ic50_nm > 0

    def test_precomputed_lookup_for_mhc_ii(self):
        """Precomputed MHC-II binder should be found via lookup."""
        result = predict_immunogenicity("HLA-DRB1*01:01", "PKYVKQNTLKLAT")
        assert result.method == "precomputed_lookup"
        assert result.confidence == "high"

    def test_unknown_allele_returns_non_binder(self):
        """Allele with no PSSM should return non_binder."""
        result = predict_immunogenicity("HLA-XXX*99:99", "GILGFVFTL")
        assert result.binding_class == "non_binder"
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert result.confidence == "low"

    def test_wrong_length_mhc_i_returns_non_binder(self):
        """MHC-I peptide with wrong length should return non_binder."""
        result = predict_immunogenicity("HLA-A*02:01", "GILGFVFT")  # 8-mer instead of 9
        assert result.binding_class == "non_binder"

    def test_result_fields_present(self):
        """ImmunogenicityPrediction has all required fields."""
        result = predict_immunogenicity("HLA-A*02:01", "GILGFVFTL")
        assert hasattr(result, "allele")
        assert hasattr(result, "peptide")
        assert hasattr(result, "ic50_nm")
        assert hasattr(result, "binding_class")
        assert hasattr(result, "method")
        assert hasattr(result, "confidence")

    def test_never_requires_external_tools(self):
        """predict_immunogenicity should never import mhcflurry or netmhcpan."""
        for allele in DEFAULT_MHC_I_ALLELES + DEFAULT_MHC_II_ALLELES:
            result = predict_immunogenicity(allele, "GILGFVFTL")
            assert result.method in ("precomputed_lookup", "pssm_fallback")


# ═══════════════════════════════════════════════════════════════
# TestScanPeptides
# ═══════════════════════════════════════════════════════════════

class TestScanPeptides:
    """Tests for scan_peptides() overlapping peptide generation."""

    def test_correct_number_of_peptides(self):
        """scan_peptides should generate len(protein) - peptide_length + 1 peptides."""
        protein = SHORT_PROTEIN  # All standard AAs
        results = scan_peptides(protein, "HLA-A*02:01", peptide_length=9)
        expected_count = len(protein) - 9 + 1
        assert len(results) == expected_count, (
            f"Expected {expected_count} peptides, got {len(results)}"
        )

    def test_peptide_windows_are_overlapping(self):
        """Peptide windows should overlap by peptide_length - 1 residues."""
        protein = SHORT_PROTEIN
        results = scan_peptides(protein, "HLA-A*02:01", peptide_length=9)
        # Sort by position to check overlap (results are sorted by IC50 by default)
        by_pos = sorted(results, key=lambda r: r.position)
        for i in range(len(by_pos) - 1):
            assert by_pos[i].position + 1 == by_pos[i + 1].position, (
                f"Non-consecutive positions: {by_pos[i].position} -> {by_pos[i + 1].position}"
            )

    def test_peptide_content_matches_protein(self):
        """Each peptide should be a substring of the protein at its position."""
        protein = SHORT_PROTEIN
        results = scan_peptides(protein, "HLA-A*02:01", peptide_length=9)
        for r in results:
            expected = protein[r.position : r.position + 9]
            assert r.peptide == expected, (
                f"Position {r.position}: expected '{expected}', got '{r.peptide}'"
            )

    def test_results_sorted_by_ic50(self):
        """Results should be sorted by IC50 (strongest binder first)."""
        results = scan_peptides(HBB_HUMAN, "HLA-A*02:01", peptide_length=9)
        if len(results) < 2:
            pytest.skip("Need at least 2 results to test sorting")
        for i in range(len(results) - 1):
            assert results[i].ic50_nm <= results[i + 1].ic50_nm, (
                f"Results not sorted by IC50: {results[i].ic50_nm} > {results[i + 1].ic50_nm}"
            )

    def test_empty_protein_returns_empty(self):
        """Empty protein should return empty list."""
        assert scan_peptides("", "HLA-A*02:01") == []
        assert scan_peptides("", "HLA-A*02:01", peptide_length=15) == []

    def test_short_protein_returns_empty(self):
        """Protein shorter than peptide_length should return empty list."""
        assert scan_peptides("ACDE", "HLA-A*02:01", peptide_length=9) == []

    def test_result_type_is_peptide_result(self):
        """Each result should be a PeptideResult instance."""
        results = scan_peptides(HBB_HUMAN, "HLA-A*02:01", peptide_length=9)
        for r in results:
            assert isinstance(r, PeptideResult)

    def test_peptide_result_fields(self):
        """PeptideResult should have position, peptide, ic50_nm, binding_class."""
        results = scan_peptides(HBB_HUMAN, "HLA-A*02:01", peptide_length=9)
        if not results:
            pytest.skip("No results")
        r = results[0]
        assert hasattr(r, "position")
        assert hasattr(r, "peptide")
        assert hasattr(r, "ic50_nm")
        assert hasattr(r, "binding_class")
        assert isinstance(r.position, int)
        assert isinstance(r.peptide, str)
        assert isinstance(r.ic50_nm, float)
        assert r.binding_class in (
            "strong_binder", "moderate_binder", "weak_binder", "non_binder",
        )

    def test_single_residue_protein(self):
        """Protein of exactly peptide_length should yield one peptide."""
        results = scan_peptides("ACDEFGHIK", "HLA-A*02:01", peptide_length=9)
        assert len(results) == 1
        assert results[0].peptide == "ACDEFGHIK"


# ═══════════════════════════════════════════════════════════════
# TestPSSMFallback
# ═══════════════════════════════════════════════════════════════

class TestPSSMFallback:
    """Tests that PSSM fallback produces reasonable results."""

    def test_known_binder_scores_high(self):
        """Known HLA-A*02:01 binder GILGFVFTL should score well via PSSM."""
        score = score_peptide_pssm("GILGFVFTL", "HLA-A*02:01")
        assert score > 0.0, "Known binder should have positive PSSM score"
        ic50 = binding_score_to_ic50(score)
        assert ic50 < 5000.0, f"Known binder IC50 should be < 5000 nM, got {ic50:.1f}"

    def test_random_peptide_ic50_reasonable(self):
        """Random peptide should produce an IC50 in reasonable range."""
        result = predict_immunogenicity("HLA-A*02:01", "MTKRRSGNM")
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert 1.0 <= result.ic50_nm <= 100000.0, (
            f"IC50 {result.ic50_nm:.1f} out of reasonable range"
        )

    def test_binding_class_is_valid(self):
        """Binding classification should always be a valid value."""
        valid_classes = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        for peptide in ["GILGFVFTL", "MTKRRSGNM", "AAAAAAAAA", "LLLLLLLLL"]:
            result = predict_immunogenicity("HLA-A*02:01", peptide)
            assert result.binding_class in valid_classes, (
                f"Invalid binding class: {result.binding_class}"
            )

    def test_pssm_ic50_monotonic_with_score(self):
        """Higher PSSM scores should correspond to lower IC50 values."""
        peptides = ["GILGFVFTL", "AAAAAAAAA", "LLLLLLLLL", "DEKRDEKRD"]
        scores_ic50s = []
        for pep in peptides:
            score = score_peptide_pssm(pep, "HLA-A*02:01")
            ic50 = binding_score_to_ic50(score)
            scores_ic50s.append((score, ic50))
        # Sort by score descending; IC50 should also be sorted (ascending)
        scores_ic50s.sort(key=lambda x: -x[0])
        for i in range(len(scores_ic50s) - 1):
            s1, ic50_1 = scores_ic50s[i]
            s2, ic50_2 = scores_ic50s[i + 1]
            assert s1 >= s2, "Scores not in descending order"
            if s1 > s2:
                assert ic50_1 <= ic50_2, (
                    f"Score {s1:.4f} > {s2:.4f} but IC50 {ic50_1:.1f} > {ic50_2:.1f}"
                )

    def test_confidence_levels_valid(self):
        """Confidence should be one of: high, medium, low."""
        valid_confidences = {"high", "medium", "low"}
        for peptide in ["GILGFVFTL", "MTKRRSGNM"]:
            result = predict_immunogenicity("HLA-A*02:01", peptide)
            assert result.confidence in valid_confidences, (
                f"Invalid confidence: {result.confidence}"
            )


# ═══════════════════════════════════════════════════════════════
# TestMHCII
# ═══════════════════════════════════════════════════════════════

class TestMHCII:
    """Tests for MHC-II alleles with 15-mer peptides."""

    def test_mhc_ii_15mer_scan(self):
        """scan_peptides with 15-mer window for MHC-II allele."""
        results = scan_peptides(HBB_HUMAN, "HLA-DRB1*01:01", peptide_length=15)
        assert isinstance(results, list)
        assert len(results) > 0, "Should produce at least one 15-mer peptide"
        expected_count = len(HBB_HUMAN) - 15 + 1
        assert len(results) == expected_count, (
            f"Expected {expected_count} 15-mers, got {len(results)}"
        )

    def test_mhc_ii_peptide_length(self):
        """Each MHC-II peptide should be 15 residues long."""
        results = scan_peptides(HBB_HUMAN, "HLA-DRB1*01:01", peptide_length=15)
        for r in results:
            assert len(r.peptide) == 15, (
                f"Expected 15-mer, got {len(r.peptide)}-mer: {r.peptide}"
            )

    def test_mhc_ii_core_scanning(self):
        """MHC-II predict_immunogenicity should scan 9-mer cores within the peptide."""
        result = predict_immunogenicity("HLA-DRB1*01:01", "PKYVKQNTLKLATLV")
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert result.ic50_nm < 100000.0

    def test_mhc_ii_binding_classes(self):
        """MHC-II binding classes should be valid."""
        valid_classes = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        results = scan_peptides(HBB_HUMAN, "HLA-DRB1*01:01", peptide_length=15)
        for r in results:
            assert r.binding_class in valid_classes, (
                f"Invalid binding class: {r.binding_class}"
            )

    def test_mhc_ii_all_default_alleles(self):
        """All default MHC-II alleles should work with scan_peptides."""
        for allele in DEFAULT_MHC_II_ALLELES:
            results = scan_peptides(HBB_HUMAN, allele, peptide_length=15)
            assert isinstance(results, list), f"scan_peptides failed for {allele}"
            assert len(results) > 0, f"No results for {allele}"

    def test_mhc_ii_short_peptide_returns_non_binder(self):
        """Peptide shorter than 9 (MHC-II core) should return non_binder."""
        result = predict_immunogenicity("HLA-DRB1*01:01", "PKYVKQN")
        assert result.binding_class == "non_binder"

    def test_mhc_ii_exact_core_length(self):
        """A 9-mer peptide given to MHC-II should be scored as a single core."""
        result = predict_immunogenicity("HLA-DRB1*01:01", "PKYVKQNTL")
        assert result.method in ("pssm_fallback", "precomputed_lookup")
        assert result.ic50_nm > 0

    def test_mhc_ii_results_sorted_by_ic50(self):
        """scan_peptides for MHC-II should be sorted by IC50."""
        results = scan_peptides(HBB_HUMAN, "HLA-DRB1*01:01", peptide_length=15)
        if len(results) < 2:
            pytest.skip("Need at least 2 results")
        for i in range(len(results) - 1):
            assert results[i].ic50_nm <= results[i + 1].ic50_nm


# ═══════════════════════════════════════════════════════════════
# TestImmunogenicityOfflinePredicates
# ═══════════════════════════════════════════════════════════════

class TestImmunogenicityOfflinePredicates:
    """Tests for immuno_predicates working with the offline API."""

    def test_evaluate_low_immunogenicity_returns_type_check_result(self):
        """evaluate_low_immunogenicity should return a TypeCheckResult."""
        result = evaluate_low_immunogenicity(HBB_HUMAN)
        assert isinstance(result, TypeCheckResult)
        assert result.verdict in (
            Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
            Verdict.LIKELY_FAIL, Verdict.FAIL,
        )

    def test_evaluate_no_strong_t_cell_returns_type_check_result(self):
        """evaluate_no_strong_t_cell_epitope should return a TypeCheckResult."""
        result = evaluate_no_strong_t_cell_epitope(HBB_HUMAN)
        assert isinstance(result, TypeCheckResult)

    def test_evaluate_no_dominant_b_cell_returns_type_check_result(self):
        """evaluate_no_dominant_b_cell_epitope should return a TypeCheckResult."""
        result = evaluate_no_dominant_b_cell_epitope(HBB_HUMAN)
        assert isinstance(result, TypeCheckResult)

    def test_evaluate_population_coverage_returns_type_check_result(self):
        """evaluate_population_coverage_safe should return a TypeCheckResult."""
        result = evaluate_population_coverage_safe(HBB_HUMAN)
        assert isinstance(result, TypeCheckResult)

    def test_all_predicates_handle_empty_protein(self):
        """All predicates should handle empty protein gracefully."""
        for eval_fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = eval_fn("")
            assert isinstance(result, TypeCheckResult)
            assert result.verdict == Verdict.UNCERTAIN

    def test_predicates_never_require_external_tools(self):
        """All predicates should work without any external tools."""
        for protein in [HBB_HUMAN, "GILGFVFTL"]:
            evaluate_low_immunogenicity(protein)
            evaluate_no_strong_t_cell_epitope(protein)
            evaluate_no_dominant_b_cell_epitope(protein)
            evaluate_population_coverage_safe(protein)
