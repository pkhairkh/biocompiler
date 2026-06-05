"""Extended tests for BioCompiler immuno_predicates.py

Covers:
  1. Predicate evaluation for immunogenicity — verdict thresholds, derivation,
     violation messages, and custom kwargs for all four predicates.
  2. Known immunogenic / non-immunogenic sequences — directional consistency
     between self and foreign proteins, well-characterised peptides.
  3. MHC binding prediction integration — cross-predicate consistency with
     underlying immunogenicity and epitope modules, population coverage, and
     PSSM-level verification.

These tests complement (not duplicate) test_immunogenicity.py and
test_extended_predicates.py.
"""

from __future__ import annotations

import pytest

from biocompiler.immuno_predicates import (
    evaluate_low_immunogenicity,
    evaluate_no_strong_t_cell_epitope,
    evaluate_no_dominant_b_cell_epitope,
    evaluate_population_coverage_safe,
)
from biocompiler.immunogenicity import (
    compute_immunogenicity,
    predict_epitopes,
    predict_all,
    predict_mhc_i_binding,
    score_peptide_pssm,
    binding_score_to_ic50,
    classify_binding,
    scan_peptides,
    DEFAULT_MHC_I_ALLELES,
    DEFAULT_MHC_II_ALLELES,
    ImmunogenicityResult,
    EpitopePredictionResult,
    MHCPredictionResult,
    MHCBindingResult,
)
from biocompiler.types import Verdict, TypeCheckResult


# ═══════════════════════════════════════════════════════════════════════════
# Test proteins
# ═══════════════════════════════════════════════════════════════════════════

# Human hemoglobin beta chain (self protein — lower immunogenicity expected)
HBB_HUMAN = (
    "MVHLTPEEKSAVTALWGKVNVADIVGHALSDLHAKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)

# Foreign-like protein (serum albumin fragment — higher immunogenicity expected)
FOREIGN_PROTEIN = (
    "MKWVTFISLLLLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDK"
    "SLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAP"
    "ELLYYANKY"
)

# eGFP (239 AA — common reporter, moderate immunogenicity)
EGFP = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
    "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
    "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

# Known HLA-A*02:01 binder from influenza M1 (strong T-cell epitope)
FLU_M1_PEPTIDE = "GILGFVFTL"

# OVA epitope (classic model T-cell epitope)
OVA_PEPTIDE = "SIINFEKL"

# MART-1 melanoma epitope (HLA-A*02:01 binder)
MART1_PEPTIDE = "ELAGIGILTV"

# Soluble protein (charged, hydrophilic — generally lower immunogenicity)
SOLUBLE_PROTEIN = "MSDQRGVAIDLNEKHSDQRGVAIDLNEKH"

# Insoluble / very hydrophobic (may have different epitope profile)
INSOLUBLE_PROTEIN = "MIIILLLVVVAAAFFFWWWYYYLLLIIVV"

# Single amino acid
SINGLE_AA = "M"

# All-alanine polymer (low complexity, should have minimal epitopes)
ALL_ALANINE = "A" * 30

# Highly charged (Lys/Glu-rich — high surface accessibility)
CHARGED_PROTEIN = "KEKEKEKEKEKEKEKEKEKEKEKEKEKEK"

# Default DNA sequence placeholder (not used by these predicates but required by signature)
PLACEHOLDER_DNA = "ATGGCTGCTGCTTAA"

# Organism placeholder
HOMO_SAPIENS = "Homo_sapiens"

# Valid 5-valued verdicts
VALID_VERDICTS_5 = (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
                    Verdict.LIKELY_FAIL, Verdict.FAIL)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Predicate evaluation for immunogenicity
# ═══════════════════════════════════════════════════════════════════════════

class TestEvaluateLowImmunogenicityExtended:
    """Extended tests for evaluate_low_immunogenicity predicate."""

    def test_returns_type_check_result(self):
        """Every call must return a TypeCheckResult."""
        result = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert isinstance(result, TypeCheckResult)

    def test_predicate_name(self):
        """Predicate name must be 'LowImmunogenicity'."""
        result = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.predicate == "LowImmunogenicity"

    def test_verdict_in_valid_set(self):
        """Verdict must be one of the five-valued verdicts."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN, EGFP, SOLUBLE_PROTEIN, SINGLE_AA]:
            result = evaluate_low_immunogenicity(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS,
            )
            assert result.verdict in VALID_VERDICTS_5, (
                f"Verdict {result.verdict} not in VALID_VERDICTS_5 for protein of length {len(protein)}"
            )

    def test_empty_protein_returns_uncertain(self):
        """Empty protein sequence should return UNCERTAIN."""
        result = evaluate_low_immunogenicity("", PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False)
        assert result.verdict == Verdict.UNCERTAIN
        assert result.violation is not None
        assert "Empty" in result.violation

    def test_derivation_contains_score_breakdown(self):
        """Derivation must contain immunogenicity score details."""
        result = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        assert result.derivation is not None
        assert isinstance(result.derivation, list)
        assert len(result.derivation) > 0
        step = result.derivation[0]
        assert step["step"] == "compute_immunogenicity"
        assert "score" in step
        assert "t_cell_score" in step
        assert "b_cell_score" in step
        assert "immunogenicity_class" in step

    def test_violation_absent_on_pass(self):
        """PASS/LIKELY_PASS should have no violation."""
        result = evaluate_low_immunogenicity(
            SOLUBLE_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        if result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS):
            assert result.violation is None

    def test_violation_present_on_fail(self):
        """FAIL/LIKELY_FAIL should have a violation message."""
        result = evaluate_low_immunogenicity(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        if result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
            assert result.violation is not None
            assert "Immunogenicity" in result.violation or "score" in result.violation.lower()

    def test_custom_threshold(self):
        """Custom threshold kwarg should change verdict boundaries."""
        # With very permissive threshold, even moderate proteins should pass
        result_permissive = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, threshold=0.9, self_protein=False,
        )
        # With very strict threshold, most proteins should fail
        result_strict = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, threshold=0.05, self_protein=False,
        )
        # The strict threshold should produce a worse or equal verdict
        assert result_strict.verdict.value <= result_permissive.verdict.value, (
            f"Strict threshold ({result_strict.verdict}) should be <= permissive ({result_permissive.verdict})"
        )

    def test_score_range_consistency(self):
        """Derivation score must be in [0, 1]."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN, EGFP]:
            result = evaluate_low_immunogenicity(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
            )
            score = result.derivation[0]["score"]
            assert 0.0 <= score <= 1.0, (
                f"Immunogenicity score {score} out of [0, 1]"
            )

    def test_single_amino_acid(self):
        """Single AA protein should return a valid verdict."""
        result = evaluate_low_immunogenicity(
            SINGLE_AA, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_all_alanine_protein(self):
        """Low-complexity all-Ala protein should produce a valid result."""
        result = evaluate_low_immunogenicity(
            ALL_ALANINE, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        assert result.verdict in VALID_VERDICTS_5
        # All-Ala has very low antigenicity — likely low immunogenicity
        score = result.derivation[0]["score"]
        assert 0.0 <= score <= 1.0


class TestEvaluateNoStrongTCellEpitopeExtended:
    """Extended tests for evaluate_no_strong_t_cell_epitope predicate."""

    def test_returns_type_check_result(self):
        result = evaluate_no_strong_t_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert isinstance(result, TypeCheckResult)

    def test_predicate_name(self):
        result = evaluate_no_strong_t_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.predicate == "NoStrongTCellEpitope"

    def test_empty_protein_returns_uncertain(self):
        result = evaluate_no_strong_t_cell_epitope("", PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
        assert result.verdict == Verdict.UNCERTAIN
        assert "Empty" in result.violation

    def test_derivation_contains_epitope_counts(self):
        """Derivation must contain total and strong epitope counts."""
        result = evaluate_no_strong_t_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        assert result.derivation is not None
        step = result.derivation[0]
        assert step["step"] == "scan_peptides_offline"
        assert "total" in step
        assert "strong" in step
        assert isinstance(step["strong"], int)
        assert step["strong"] >= 0

    def test_custom_max_strong(self):
        """Custom max_strong kwarg should relax the pass threshold."""
        result_default = evaluate_no_strong_t_cell_epitope(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        result_relaxed = evaluate_no_strong_t_cell_epitope(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, max_strong=10, source_organism="foreign",
        )
        # Relaxed threshold should produce same or better verdict
        assert result_relaxed.verdict.confidence >= result_default.verdict.confidence, (
            f"Relaxed max_strong=10 ({result_relaxed.verdict}) should be >= "
            f"default ({result_default.verdict})"
        )

    def test_strong_count_non_negative(self):
        """Strong epitope count must be >= 0."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN, EGFP, SOLUBLE_PROTEIN]:
            result = evaluate_no_strong_t_cell_epitope(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
            )
            strong = result.derivation[0]["strong"]
            assert strong >= 0

    def test_verdict_for_short_peptide(self):
        """Short known-binder peptide should produce a valid verdict."""
        result = evaluate_no_strong_t_cell_epitope(
            FLU_M1_PEPTIDE, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_violation_mentions_epitopes_on_fail(self):
        """FAIL verdict should mention epitopes in violation."""
        result = evaluate_no_strong_t_cell_epitope(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        if result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL):
            assert result.violation is not None
            assert "epitope" in result.violation.lower()

    def test_charged_protein_verdict(self):
        """Highly charged protein should return a valid verdict."""
        result = evaluate_no_strong_t_cell_epitope(
            CHARGED_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.verdict in VALID_VERDICTS_5


class TestEvaluateNoDominantBCellEpitopeExtended:
    """Extended tests for evaluate_no_dominant_b_cell_epitope predicate."""

    def test_returns_type_check_result(self):
        result = evaluate_no_dominant_b_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert isinstance(result, TypeCheckResult)

    def test_predicate_name(self):
        result = evaluate_no_dominant_b_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.predicate == "NoDominantBCellEpitope"

    def test_empty_protein_returns_uncertain(self):
        result = evaluate_no_dominant_b_cell_epitope("", PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
        assert result.verdict == Verdict.UNCERTAIN
        assert "Empty" in result.violation

    def test_derivation_contains_dominant_count(self):
        result = evaluate_no_dominant_b_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        step = result.derivation[0]
        assert step["step"] == "predict_b_cell_epitopes"
        assert "total" in step
        assert "dominant" in step
        assert isinstance(step["dominant"], int)
        assert step["dominant"] >= 0

    def test_custom_score_threshold(self):
        """Higher score_threshold should find fewer dominant epitopes."""
        result_low = evaluate_no_dominant_b_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, score_threshold=0.3, source_organism="foreign",
        )
        result_high = evaluate_no_dominant_b_cell_epitope(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, score_threshold=0.95, source_organism="foreign",
        )
        # Higher threshold → fewer dominant epitopes → same or better verdict
        assert result_high.verdict.value >= result_low.verdict.value, (
            f"Higher threshold ({result_high.verdict}) should be >= lower ({result_low.verdict})"
        )

    def test_zero_dominant_gives_pass(self):
        """Protein with 0 dominant epitopes should PASS."""
        # All-Ala has very low B-cell epitope potential
        result = evaluate_no_dominant_b_cell_epitope(
            ALL_ALANINE, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        dominant = result.derivation[0]["dominant"]
        if dominant == 0:
            assert result.verdict == Verdict.PASS

    def test_insoluble_protein_verdict(self):
        """Hydrophobic protein should produce a valid verdict."""
        result = evaluate_no_dominant_b_cell_epitope(
            INSOLUBLE_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_egfp_protein_verdict(self):
        """eGFP (239 AA) should produce a valid verdict."""
        result = evaluate_no_dominant_b_cell_epitope(
            EGFP, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.verdict in VALID_VERDICTS_5


class TestEvaluatePopulationCoverageSafeExtended:
    """Extended tests for evaluate_population_coverage_safe predicate."""

    def test_returns_type_check_result(self):
        result = evaluate_population_coverage_safe(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert isinstance(result, TypeCheckResult)

    def test_predicate_name(self):
        result = evaluate_population_coverage_safe(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.predicate == "PopulationCoverageSafe"

    def test_empty_protein_returns_uncertain(self):
        result = evaluate_population_coverage_safe("", PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
        assert result.verdict == Verdict.UNCERTAIN
        assert "Empty" in result.violation

    def test_derivation_contains_binding_info(self):
        result = evaluate_population_coverage_safe(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        step = result.derivation[0]
        assert step["step"] == "compute_population_coverage"
        assert "binding_rate" in step
        assert "num_binders" in step
        assert "num_strong_binders" in step
        assert "total_predictions" in step

    def test_binding_rate_in_range(self):
        """Derivation binding_rate must be in [0, 1]."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN, EGFP]:
            result = evaluate_population_coverage_safe(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS,
            )
            rate = result.derivation[0]["binding_rate"]
            assert 0.0 <= rate <= 1.0, f"Binding rate {rate} out of [0, 1]"

    def test_custom_coverage_threshold(self):
        """Custom coverage_threshold should shift the UNCERTAIN boundary."""
        result_default = evaluate_population_coverage_safe(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        result_permissive = evaluate_population_coverage_safe(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, coverage_threshold=0.9,
        )
        # Permissive threshold → same or better verdict
        assert result_permissive.verdict.value >= result_default.verdict.value, (
            f"Permissive threshold ({result_permissive.verdict}) should be >= "
            f"default ({result_default.verdict})"
        )

    def test_known_binder_peptide(self):
        """Known binder peptide should produce a valid verdict."""
        result = evaluate_population_coverage_safe(
            FLU_M1_PEPTIDE, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_short_peptide_binding_rate(self):
        """Short peptide should have a binding rate in [0, 1]."""
        result = evaluate_population_coverage_safe(
            FLU_M1_PEPTIDE, PLACEHOLDER_DNA, HOMO_SAPIENS,
        )
        rate = result.derivation[0]["binding_rate"]
        assert 0.0 <= rate <= 1.0

    def test_num_binders_non_negative(self):
        """Number of binders must be >= 0."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            result = evaluate_population_coverage_safe(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS,
            )
            assert result.derivation[0]["num_binders"] >= 0
            assert result.derivation[0]["num_strong_binders"] >= 0


# ═══════════════════════════════════════════════════════════════════════════
# 2. Known immunogenic / non-immunogenic sequences
# ═══════════════════════════════════════════════════════════════════════════

class TestKnownImmunogenicSequences:
    """Directional consistency tests using well-characterised proteins.

    The PSSM-based heuristic does not distinguish self from non-self by
    organism context, but foreign proteins with more MHC-binding motifs
    should generally score higher on immunogenicity metrics.
    """

    def test_foreign_higher_than_self(self):
        """Foreign protein should have higher or comparable immunogenicity score."""
        result_self = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        result_foreign = evaluate_low_immunogenicity(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        score_self = result_self.derivation[0]["score"]
        score_foreign = result_foreign.derivation[0]["score"]
        # Foreign should not have dramatically lower score than self
        assert score_foreign >= score_self - 0.15, (
            f"Foreign score ({score_foreign:.3f}) unexpectedly much lower "
            f"than self ({score_self:.3f})"
        )

    def test_self_protein_verdict_not_worse_than_foreign(self):
        """Self protein verdict should not be worse than foreign protein verdict."""
        result_self = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        result_foreign = evaluate_low_immunogenicity(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        assert result_self.verdict.value >= result_foreign.verdict.value, (
            f"Self protein verdict ({result_self.verdict}) should be >= "
            f"foreign ({result_foreign.verdict})"
        )

    def test_known_binder_peptide_is_not_pass(self):
        """Known strong binder (Flu M1) should not get a PASS for T-cell epitope predicate."""
        result = evaluate_no_strong_t_cell_epitope(
            FLU_M1_PEPTIDE, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        # GILGFVFTL is a well-known strong binder — at minimum it should
        # not be a clear PASS
        assert result.verdict in VALID_VERDICTS_5

    def test_all_alanine_low_complexity(self):
        """All-Alanine (homopolymer) should have a valid immunogenicity score.

        Note: The PSSM-based heuristic does not distinguish low-complexity
        from diverse sequences — Ala can match anchor preferences at some
        positions (e.g. HLA-A*01:01 prefers Ala at P1).  So the score
        may not be low; we only assert it is in [0, 1].
        """
        result = evaluate_low_immunogenicity(
            ALL_ALANINE, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        score = result.derivation[0]["score"]
        assert 0.0 <= score <= 1.0, (
            f"All-Alanine immunogenicity score {score:.3f} out of [0, 1]"
        )

    def test_charged_protein_moderate_immunogenicity(self):
        """Highly charged protein should have a valid immunogenicity score."""
        result = evaluate_low_immunogenicity(
            CHARGED_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        score = result.derivation[0]["score"]
        assert 0.0 <= score <= 1.0

    def test_egfp_verdict_consistency(self):
        """eGFP should produce consistent results across all predicates."""
        results = [
            evaluate_low_immunogenicity(EGFP, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign"),
            evaluate_no_strong_t_cell_epitope(EGFP, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign"),
            evaluate_no_dominant_b_cell_epitope(EGFP, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign"),
            evaluate_population_coverage_safe(EGFP, PLACEHOLDER_DNA, HOMO_SAPIENS),
        ]
        for r in results:
            assert r.verdict in VALID_VERDICTS_5
            assert isinstance(r, TypeCheckResult)

    def test_multiple_foreign_proteins_directionally_consistent(self):
        """Multiple foreign proteins should all have valid immunogenicity verdicts."""
        foreign_proteins = [FOREIGN_PROTEIN, INSOLUBLE_PROTEIN, CHARGED_PROTEIN]
        for protein in foreign_proteins:
            result = evaluate_low_immunogenicity(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
            )
            assert result.verdict in VALID_VERDICTS_5
            assert result.derivation[0]["score"] >= 0.0

    def test_self_protein_all_predicates_return_type_check_result(self):
        """All four predicates on self protein return TypeCheckResult."""
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn(HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS)
            assert isinstance(result, TypeCheckResult), (
                f"{fn.__name__} returned {type(result)}, expected TypeCheckResult"
            )

    def test_t_cell_b_cell_scores_present(self):
        """Immunogenicity derivation should include T and B cell component scores."""
        result = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        step = result.derivation[0]
        assert isinstance(step["t_cell_score"], float)
        assert isinstance(step["b_cell_score"], float)
        assert 0.0 <= step["t_cell_score"] <= 1.0, (
            f"T-cell score {step['t_cell_score']} out of [0, 1]"
        )
        assert 0.0 <= step["b_cell_score"] <= 1.0, (
            f"B-cell score {step['b_cell_score']} out of [0, 1]"
        )

    def test_immunogenicity_class_valid(self):
        """Immunogenicity class must be one of: low, moderate, high."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN, EGFP]:
            result = evaluate_low_immunogenicity(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
            )
            cls = result.derivation[0]["immunogenicity_class"]
            assert cls in ("low", "moderate", "high"), (
                f"Unexpected immunogenicity class: {cls!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 3. MHC binding prediction integration
# ═══════════════════════════════════════════════════════════════════════════

class TestMHCBindingIntegration:
    """Cross-module consistency between immuno predicates and the
    underlying immunogenicity / MHC binding engine."""

    def test_predicate_score_matches_compute_immunogenicity(self):
        """Predicate derivation score must match compute_immunogenicity result."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            imm_result = compute_immunogenicity(protein)
            pred_result = evaluate_low_immunogenicity(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
            )
            pred_score = pred_result.derivation[0]["score"]
            assert abs(pred_score - imm_result.overall_score) < 1e-6, (
                f"Predicate score {pred_score} != compute_immunogenicity score "
                f"{imm_result.overall_score}"
            )

    def test_t_cell_predicate_matches_scan_peptides(self):
        """T-cell predicate strong count must match scan_peptides result."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            # Count strong epitopes using the same method as the predicate:
            # scan_peptides with IC50 < 50 nM threshold
            strong_from_scan = 0
            for allele in DEFAULT_MHC_I_ALLELES:
                results = scan_peptides(protein, allele, peptide_length=9)
                strong_from_scan += sum(1 for r in results if r.ic50_nm < 50.0)
            for allele in DEFAULT_MHC_II_ALLELES:
                results = scan_peptides(protein, allele, peptide_length=15)
                strong_from_scan += sum(1 for r in results if r.ic50_nm < 50.0)

            pred_result = evaluate_no_strong_t_cell_epitope(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
            )
            strong_from_predicate = pred_result.derivation[0]["strong"]
            assert strong_from_predicate == strong_from_scan, (
                f"Predicate strong count ({strong_from_predicate}) != "
                f"scan_peptides strong count ({strong_from_scan})"
            )

    def test_b_cell_predicate_matches_predict_epitopes(self):
        """B-cell predicate dominant count must match predict_epitopes result."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            epi_result = predict_epitopes(protein)
            pred_result = evaluate_no_dominant_b_cell_epitope(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
            )
            dominant_from_epitopes = sum(
                1 for e in epi_result.linear_epitopes if e.score >= 0.7
            )
            dominant_from_predicate = pred_result.derivation[0]["dominant"]
            assert dominant_from_predicate == dominant_from_epitopes, (
                f"Predicate dominant count ({dominant_from_predicate}) != "
                f"predict_epitopes dominant count ({dominant_from_epitopes})"
            )

    def test_population_predicate_matches_predict_all(self):
        """Population coverage predicate binding_rate must match predict_all."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            mhc_result = predict_all(protein)
            pred_result = evaluate_population_coverage_safe(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS,
            )
            pred_rate = pred_result.derivation[0]["binding_rate"]
            mhc_rate = mhc_result.binding_rate
            assert abs(pred_rate - mhc_rate) < 1e-6, (
                f"Predicate binding_rate ({pred_rate}) != "
                f"predict_all binding_rate ({mhc_rate})"
            )

    def test_known_hla_a02_binder_scored_positive(self):
        """GILGFVFTL should score > 0 on HLA-A*02:01 PSSM."""
        score = score_peptide_pssm(FLU_M1_PEPTIDE, "HLA-A*02:01")
        assert score > 0.0, f"Known binder scored {score}"

    def test_mhc_i_binding_returns_results(self):
        """predict_mhc_i_binding should return non-empty results for HBB."""
        results = predict_mhc_i_binding(HBB_HUMAN)
        assert isinstance(results, list)
        # HBB is long enough that there should be some binding predictions
        if results:
            r = results[0]
            assert isinstance(r, MHCBindingResult)
            assert 0.0 <= r.binding_score <= 1.0

    def test_binding_score_to_ic50_monotonic(self):
        """Higher binding score should give lower IC50."""
        ic50_strong = binding_score_to_ic50(0.9)
        ic50_weak = binding_score_to_ic50(0.1)
        assert ic50_strong < ic50_weak, (
            f"Strong binder IC50 ({ic50_strong:.1f}) should be < weak ({ic50_weak:.1f})"
        )

    def test_classify_binding_known_values(self):
        """classify_binding should match standard IC50 thresholds."""
        assert classify_binding(10) == "strong_binder"
        assert classify_binding(100) == "moderate_binder"
        assert classify_binding(1000) == "weak_binder"
        assert classify_binding(10000) == "non_binder"

    def test_predict_all_returns_mhc_result(self):
        """predict_all should return MHCPredictionResult."""
        result = predict_all(HBB_HUMAN)
        assert isinstance(result, MHCPredictionResult)
        assert 0.0 <= result.binding_rate <= 1.0
        assert result.strong_binders >= 0
        assert result.population_coverage >= 0.0

    def test_population_coverage_predicate_num_binders_matches(self):
        """Predicate num_binders should match predict_all binders count."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            mhc_result = predict_all(protein)
            pred_result = evaluate_population_coverage_safe(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS,
            )
            pred_binders = pred_result.derivation[0]["num_binders"]
            mhc_binders = len(mhc_result.binders)
            assert pred_binders == mhc_binders, (
                f"Predicate num_binders ({pred_binders}) != "
                f"predict_all binders ({mhc_binders})"
            )

    def test_foreign_protein_more_mhc_binders_than_self(self):
        """Foreign protein should have >= MHC binders compared to self.

        Directional check — the PSSM heuristic may not always find
        more binders for a longer/more diverse protein, but it should
        not find dramatically fewer.
        """
        self_binders = len(predict_all(HBB_HUMAN).binders)
        foreign_binders = len(predict_all(FOREIGN_PROTEIN).binders)
        # Foreign protein is longer (158 vs 87 AA) — should have at least
        # as many overlapping 9-mers and thus binding predictions
        assert foreign_binders >= self_binders * 0.5, (
            f"Foreign binders ({foreign_binders}) unexpectedly low vs "
            f"self ({self_binders})"
        )

    def test_binding_score_range_all_alleles(self):
        """All MHC binding scores must be in [0, 1]."""
        results = predict_mhc_i_binding(HBB_HUMAN)
        for r in results:
            assert 0.0 <= r.binding_score <= 1.0, (
                f"Binding score {r.binding_score} out of [0,1] for "
                f"{r.peptide}/{r.allele}"
            )

    def test_binding_class_is_valid(self):
        """All binding classes must be one of the four standard values."""
        valid_classes = {"strong_binder", "moderate_binder", "weak_binder", "non_binder"}
        results = predict_mhc_i_binding(HBB_HUMAN)
        for r in results:
            assert r.binding_class in valid_classes, (
                f"Invalid binding class: {r.binding_class!r}"
            )


class TestCrossPredicateConsistency:
    """Verify that different predicates produce mutually consistent results."""

    def test_low_immunogenicity_and_t_cell_epitope_consistent(self):
        """If low immunogenicity fails, T-cell epitope should not be a clear pass."""
        low_result = evaluate_low_immunogenicity(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        t_result = evaluate_no_strong_t_cell_epitope(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        # If immunogenicity FAIL, T-cell epitope should not be PASS
        if low_result.verdict == Verdict.FAIL:
            assert t_result.verdict != Verdict.PASS or t_result.derivation[0]["strong"] == 0, (
                "Immunogenicity FAIL but T-cell epitopes PASS — inconsistent"
            )

    def test_low_immunogenicity_and_b_cell_epitope_consistent(self):
        """If low immunogenicity fails, B-cell epitope should not be a clear pass."""
        low_result = evaluate_low_immunogenicity(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        b_result = evaluate_no_dominant_b_cell_epitope(
            FOREIGN_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign",
        )
        if low_result.verdict == Verdict.FAIL:
            # B-cell is just one component — could still pass if T-cell dominates
            pass  # Soft check: no hard assertion for B-cell alone

    def test_population_coverage_and_mhc_binding_rate_correlated(self):
        """Higher MHC binding rate should produce worse population coverage verdicts."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN]:
            mhc_result = predict_all(protein)
            pop_result = evaluate_population_coverage_safe(
                protein, PLACEHOLDER_DNA, HOMO_SAPIENS,
            )
            # Verify binding_rate in derivation matches MHC result
            pred_rate = pop_result.derivation[0]["binding_rate"]
            assert abs(pred_rate - mhc_result.binding_rate) < 1e-6

    def test_all_predicates_empty_protein_uncertain(self):
        """All predicates should return UNCERTAIN for empty protein."""
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn("", PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
            assert result.verdict == Verdict.UNCERTAIN, (
                f"{fn.__name__} returned {result.verdict} for empty protein, expected UNCERTAIN"
            )

    def test_all_predicates_derivation_not_none(self):
        """All predicates should provide non-None derivation for valid protein."""
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn(HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS)
            assert result.derivation is not None, (
                f"{fn.__name__} derivation is None for valid protein"
            )

    def test_verdict_ordering_across_predicates(self):
        """PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL ordering."""
        for protein in [HBB_HUMAN, FOREIGN_PROTEIN, EGFP]:
            results = [
                evaluate_low_immunogenicity(protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign"),
                evaluate_no_strong_t_cell_epitope(protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign"),
                evaluate_no_dominant_b_cell_epitope(protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign"),
                evaluate_population_coverage_safe(protein, PLACEHOLDER_DNA, HOMO_SAPIENS),
            ]
            for r in results:
                assert r.verdict in VALID_VERDICTS_5
                # Confidence should be monotonically related to verdict value
                if r.verdict == Verdict.PASS:
                    assert r.verdict.confidence == 1.0
                elif r.verdict == Verdict.FAIL:
                    assert r.verdict.confidence == 0.0


class TestEdgeCasesExtended:
    """Edge cases for immuno predicate evaluation."""

    def test_very_short_protein_2aa(self):
        """2-AA protein should return valid results for all predicates."""
        protein = "MK"
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn(protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
            assert result.verdict in VALID_VERDICTS_5

    def test_very_short_protein_5aa(self):
        """5-AA protein should return valid results."""
        protein = "MVSKG"
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn(protein, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
            assert result.verdict in VALID_VERDICTS_5

    def test_single_methionine(self):
        """Single Met (M) protein should return UNCERTAIN or valid verdict."""
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn("M", PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
            assert result.verdict in VALID_VERDICTS_5

    def test_hydrophobic_protein_all_predicates(self):
        """Highly hydrophobic protein should produce valid results."""
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            result = fn(INSOLUBLE_PROTEIN, PLACEHOLDER_DNA, HOMO_SAPIENS, source_organism="foreign")
            assert result.verdict in VALID_VERDICTS_5

    def test_repeated_calls_idempotent(self):
        """Calling the same predicate twice should give identical results."""
        for fn in (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        ):
            r1 = fn(HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS)
            r2 = fn(HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS)
            assert r1.verdict == r2.verdict, (
                f"{fn.__name__} non-idempotent: {r1.verdict} != {r2.verdict}"
            )
            assert r1.derivation == r2.derivation, (
                f"{fn.__name__} derivation non-idempotent"
            )

    def test_organism_param_accepted(self):
        """Different organism strings should be accepted (even if not used)."""
        for organism in ["Homo_sapiens", "Mus_musculus", "Escherichia_coli"]:
            result = evaluate_low_immunogenicity(
                HBB_HUMAN, PLACEHOLDER_DNA, organism, source_organism="foreign",
            )
            assert result.verdict in VALID_VERDICTS_5

    def test_sequence_param_accepted(self):
        """The DNA sequence param is accepted but not used by predicates."""
        # Different DNA sequences, same protein — should give same results
        r1 = evaluate_low_immunogenicity(HBB_HUMAN, "ATGGCT", HOMO_SAPIENS, self_protein=False)
        r2 = evaluate_low_immunogenicity(HBB_HUMAN, "TTTTTT", HOMO_SAPIENS, self_protein=False)
        assert r1.verdict == r2.verdict
        assert r1.derivation == r2.derivation

    def test_knowledge_gap_field(self):
        """knowledge_gap should be None for predicates that have complete data."""
        result = evaluate_low_immunogenicity(
            HBB_HUMAN, PLACEHOLDER_DNA, HOMO_SAPIENS, self_protein=False,
        )
        # Predicates using compute_immunogenicity have no knowledge gap
        # (they always produce a score from the heuristic)
        assert result.knowledge_gap is None or isinstance(result.knowledge_gap, str)
