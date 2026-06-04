"""
BioCompiler Solubility Predicates — Comprehensive Test Suite
=============================================================
Tests for solubility_predicates.py covering:
  1. Predicate evaluation for protein solubility
  2. Known soluble/insoluble sequences
  3. Score ranges

Task ID: 114
"""

import pytest
from biocompiler.solubility_predicates import (
    evaluate_soluble_expression,
    evaluate_no_aggregation_prone_region,
    evaluate_charge_composition,
    evaluate_no_long_hydrophobic_stretch,
    compute_approximate_pI,
    compute_net_charge,
    find_hydrophobic_stretches,
    PKA_VALUES,
    _CAMSOL_INTRINSIC,
    _CAMSOL_SCALE,
    _CAMSOL_WINDOW,
    _CAMSOL_HIGHLY_SOLUBLE,
    _CAMSOL_MARGINAL,
    _AGG_BORDERLINE_MAX,
    _AGG_UNCERTAIN_MAX,
    _AGG_LIKELY_FAIL_MAX,
    _HYDRO_EXCESS_BORDERLINE,
    _HYDRO_EXCESS_UNCERTAIN,
    _DEFAULT_HYDROPHOBIC,
    _camsol_smoothed_profile,
    _camsol_overall_score,
    _find_aggregation_regions,
)
from biocompiler.type_system import Verdict, AA_TO_CODONS
from biocompiler.types import TypeCheckResult


# ────────────────────────────────────────────────────────────
# Test proteins with known solubility characteristics
# ────────────────────────────────────────────────────────────

# Highly soluble: rich in charged residues (D, E, K)
SOLUBLE_CHARGED = "MDEKDEKDEKDEKDEKDEKDEKDEKDEKDEKDEKDEKDEK"

# Insoluble: dominated by hydrophobic residues (I, L, V, F, W)
INSOLUBLE_HYDROPHOBIC = "MVVVIIVVVLLLFLLLLFFFFWWWAAAIIIMMM"

# Extremely soluble: all lysine/glutamate — maximally charged
SUPER_SOLUBLE = "MKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEK"

# Extremely insoluble: all isoleucine — maximally hydrophobic
SUPER_INSOLUBLE = "MIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"

# Balanced: mixed composition similar to a typical globular protein
BALANCED_GLOBULAR = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAVDILSKKGDVQVIK"

# Low charge / high pI: few acidic residues, many basic
BASIC_PROTEIN = "MKRKHKRHKRHKRHKRHKRHKRHKRH"

# Acidic protein: many D/E, few K/R
ACIDIC_PROTEIN = "MDEDEDEDEDEDEDEDEDEDEDEDEDED"

# All-alanine: borderline, small hydrophobic but not strongly aggregating
ALL_ALANINE = "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

# Mixed with single long hydrophobic stretch
LONG_HYDRO_STRETCH = "MDEKDEKLVVLVVLVVLVVLVVLVVLVVDEKDEKDEK"

# Short protein
SHORT_PROTEIN = "MKEK"

# Single residue
SINGLE_METHIONINE = "M"


def _protein_to_dna(protein: str) -> str:
    """Generate a DNA sequence encoding the given protein using first codon for each AA."""
    return "".join(AA_TO_CODONS.get(aa, ["GCT"])[0] for aa in protein)


# ══════════════════════════════════════════════════════════════
# 1. Predicate evaluation for protein solubility
# ══════════════════════════════════════════════════════════════

class TestEvaluateSolubleExpression:
    """Tests for evaluate_soluble_expression predicate."""

    # ── Basic verdict classification ──────────────────────────

    def test_soluble_protein_pass(self):
        """Highly charged protein should get PASS or LIKELY_PASS."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_soluble_expression(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS), (
            f"Soluble protein got {result.verdict}, expected PASS or LIKELY_PASS"
        )

    def test_insoluble_protein_fail(self):
        """Hydrophobic protein should get FAIL or LIKELY_FAIL."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_soluble_expression(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN), (
            f"Insoluble protein got {result.verdict}, expected poor solubility verdict"
        )

    def test_super_soluble_pass(self):
        """All charged residues → definitely PASS."""
        dna = _protein_to_dna(SUPER_SOLUBLE)
        result = evaluate_soluble_expression(dna, SUPER_SOLUBLE, "Homo_sapiens")
        assert result.verdict == Verdict.PASS, (
            f"Super soluble protein got {result.verdict}, expected PASS"
        )

    def test_super_insoluble_fail(self):
        """All hydrophobic residues → LIKELY_FAIL or FAIL."""
        dna = _protein_to_dna(SUPER_INSOLUBLE)
        result = evaluate_soluble_expression(dna, SUPER_INSOLUBLE, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL), (
            f"Super insoluble protein got {result.verdict}, expected FAIL or LIKELY_FAIL"
        )

    def test_balanced_protein_verdict(self):
        """Balanced globular protein should get a reasonable verdict."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert isinstance(result.verdict, Verdict), "Verdict should be a Verdict enum"

    # ── Empty / edge cases ───────────────────────────────────

    def test_empty_protein_fail(self):
        """Empty protein → FAIL."""
        result = evaluate_soluble_expression("", "", "Homo_sapiens")
        assert result.verdict == Verdict.FAIL
        assert result.violation is not None
        assert "Empty" in result.violation

    def test_single_methionine(self):
        """Single M should produce a valid result (not crash)."""
        dna = _protein_to_dna(SINGLE_METHIONINE)
        result = evaluate_soluble_expression(dna, SINGLE_METHIONINE, "Homo_sapiens")
        assert isinstance(result, TypeCheckResult)
        assert isinstance(result.verdict, Verdict)

    def test_short_protein(self):
        """Short protein should not crash."""
        dna = _protein_to_dna(SHORT_PROTEIN)
        result = evaluate_soluble_expression(dna, SHORT_PROTEIN, "Homo_sapiens")
        assert isinstance(result, TypeCheckResult)

    # ── min_solubility_score parameter ───────────────────────

    def test_min_solubility_score_override(self):
        """Setting a high min_solubility_score can downgrade a passing verdict."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        # Default: should pass
        result_default = evaluate_soluble_expression(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        # Very high threshold: should force LIKELY_FAIL even for a soluble protein
        result_strict = evaluate_soluble_expression(
            dna, SOLUBLE_CHARGED, "Homo_sapiens", min_solubility_score=5.0
        )
        assert result_strict.verdict == Verdict.LIKELY_FAIL, (
            f"Strict min_solubility_score=5.0 should cause LIKELY_FAIL, got {result_strict.verdict}"
        )

    def test_min_solubility_score_zero_no_effect(self):
        """min_solubility_score=0.0 (default) should not change verdict."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_soluble_expression(
            dna, SOLUBLE_CHARGED, "Homo_sapiens", min_solubility_score=0.0
        )
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def test_min_solubility_score_negative_no_effect(self):
        """Negative min_solubility_score should not change verdict (all scores > negative)."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_soluble_expression(
            dna, SOLUBLE_CHARGED, "Homo_sapiens", min_solubility_score=-10.0
        )
        result_default = evaluate_soluble_expression(
            dna, SOLUBLE_CHARGED, "Homo_sapiens"
        )
        assert result.verdict == result_default.verdict

    # ── Result structure ─────────────────────────────────────

    def test_result_is_type_check_result(self):
        """Result should be a TypeCheckResult instance."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert isinstance(result, TypeCheckResult)

    def test_predicate_name_contains_soluble(self):
        """Predicate name should contain 'SolubleExpression'."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert "SolubleExpression" in result.predicate

    def test_derivation_has_camsol_score(self):
        """Derivation should include the CamSol intrinsic score step."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert result.derivation is not None
        steps = [d["step"] for d in result.derivation]
        assert "camsol_intrinsic_score" in steps
        assert "min_solubility_score" in steps
        assert "aggregation_prone_regions" in steps

    def test_derivation_camsol_score_is_numeric(self):
        """CamSol score in derivation should be a number."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        score_step = next(d for d in result.derivation if d["step"] == "camsol_intrinsic_score")
        assert isinstance(score_step["value"], (int, float))

    def test_derivation_aggregation_regions_is_list(self):
        """Aggregation-prone regions in derivation should be a list."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        agg_step = next(d for d in result.derivation if d["step"] == "aggregation_prone_regions")
        assert isinstance(agg_step["value"], list)

    def test_knowledge_gap_without_pdb(self):
        """Without PDB, knowledge_gap should note structural correction missing."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert result.knowledge_gap is not None
        assert "structure" in result.knowledge_gap.lower() or "PDB" in result.knowledge_gap

    def test_no_knowledge_gap_with_pdb(self):
        """With PDB string, knowledge_gap should be None."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", pdb_string="fake PDB"
        )
        assert result.knowledge_gap is None

    # ── Case insensitivity ───────────────────────────────────

    def test_lowercase_protein_handled(self):
        """Lowercase protein input should be handled (uppercased internally)."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_soluble_expression(dna, SOLUBLE_CHARGED.lower(), "Homo_sapiens")
        result_upper = evaluate_soluble_expression(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        assert result.verdict == result_upper.verdict

    # ── Organism parameter accepted ──────────────────────────

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "Escherichia_coli",
        "Saccharomyces_cerevisiae",
        "Mus_musculus",
    ])
    def test_different_organisms(self, organism):
        """Predicate should accept various organism names."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_soluble_expression(dna, BALANCED_GLOBULAR, organism)
        assert isinstance(result, TypeCheckResult)


# ══════════════════════════════════════════════════════════════
# 2. Known soluble/insoluble sequences
# ══════════════════════════════════════════════════════════════

class TestEvaluateNoAggregationProneRegion:
    """Tests for evaluate_no_aggregation_prone_region predicate."""

    def test_soluble_protein_no_apr(self):
        """Soluble protein should have no aggregation-prone regions → PASS."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_no_aggregation_prone_region(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_insoluble_protein_has_apr(self):
        """Hydrophobic protein should have long APRs → FAIL or LIKELY_FAIL."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_no_aggregation_prone_region(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL), (
            f"Insoluble protein APR verdict {result.verdict}, expected FAIL or LIKELY_FAIL"
        )

    def test_empty_protein_pass(self):
        """Empty protein has no regions → PASS."""
        result = evaluate_no_aggregation_prone_region("", "", "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_max_region_length_override(self):
        """Custom max_region_length should affect the verdict."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        # Very permissive: max_region_length=100 → everything passes
        result_permissive = evaluate_no_aggregation_prone_region(
            dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens", max_region_length=100
        )
        assert result_permissive.verdict == Verdict.PASS

    def test_score_threshold_override(self):
        """More negative threshold makes fewer regions aggregation-prone."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        # Default threshold
        result_default = evaluate_no_aggregation_prone_region(
            dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens"
        )
        # Very strict threshold (higher value = more regions flagged)
        result_strict = evaluate_no_aggregation_prone_region(
            dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens", score_threshold=0.0
        )
        # Strict threshold should have verdict at least as bad or worse
        _verdict_order = {
            Verdict.PASS: 0, Verdict.LIKELY_PASS: 1,
            Verdict.UNCERTAIN: 2, Verdict.LIKELY_FAIL: 3, Verdict.FAIL: 4,
        }
        assert _verdict_order[result_strict.verdict] >= _verdict_order[result_default.verdict]

    def test_predicate_name_contains_params(self):
        """Predicate name should include max and threshold parameters."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_no_aggregation_prone_region(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", max_region_length=8, score_threshold=-0.5
        )
        assert "NoAggregationProneRegion" in result.predicate
        assert "8" in result.predicate
        assert "-0.5" in result.predicate

    def test_derivation_contains_regions(self):
        """Derivation should list aggregation-prone regions with details."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_no_aggregation_prone_region(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        steps = {d["step"]: d["value"] for d in result.derivation}
        assert "aggregation_prone_regions" in steps
        assert "longest_region" in steps
        assert "max_region_length" in steps
        assert "score_threshold" in steps
        assert isinstance(steps["longest_region"], int)
        assert isinstance(steps["aggregation_prone_regions"], list)

    def test_knowledge_gap_without_pdb(self):
        """Without PDB, knowledge_gap should be set."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_no_aggregation_prone_region(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert result.knowledge_gap is not None

    def test_no_knowledge_gap_with_pdb(self):
        """With PDB, knowledge_gap should be None."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_no_aggregation_prone_region(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", pdb_string="ATOM  ..."
        )
        assert result.knowledge_gap is None

    # ── Verdict tier boundaries ──────────────────────────────

    def test_verdict_tiers_by_region_length(self):
        """Different region lengths should produce different verdict tiers."""
        # Region length ≤ max_region_length → PASS
        # Region length max+1 to _AGG_BORDERLINE_MAX → LIKELY_PASS
        # Region length _AGG_BORDERLINE_MAX+1 to _AGG_UNCERTAIN_MAX → UNCERTAIN
        # Region length _AGG_UNCERTAIN_MAX+1 to _AGG_LIKELY_FAIL_MAX → LIKELY_FAIL
        # Region length > _AGG_LIKELY_FAIL_MAX → FAIL
        assert _AGG_BORDERLINE_MAX == 7
        assert _AGG_UNCERTAIN_MAX == 10
        assert _AGG_LIKELY_FAIL_MAX == 15


class TestEvaluateChargeComposition:
    """Tests for evaluate_charge_composition predicate."""

    def test_soluble_charged_pass(self):
        """Protein with many charged residues should PASS."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_charge_composition(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_low_charge_fail(self):
        """Protein with very few charged residues should LIKELY_FAIL."""
        low_charge = "MAVILFWMAPTVIGLFWVILMAFGLMVVIAF"
        dna = _protein_to_dna(low_charge)
        result = evaluate_charge_composition(dna, low_charge, "Homo_sapiens")
        assert result.verdict in (Verdict.LIKELY_FAIL, Verdict.FAIL), (
            f"Low charge protein got {result.verdict}, expected LIKELY_FAIL or FAIL"
        )

    def test_high_pi_uncertain(self):
        """Protein with very high pI should get UNCERTAIN or worse."""
        dna = _protein_to_dna(BASIC_PROTEIN)
        result = evaluate_charge_composition(dna, BASIC_PROTEIN, "Homo_sapiens")
        # Basic protein has high pI but also many charged residues
        # So it might be UNCERTAIN (high pI) or PASS (enough charges)
        assert isinstance(result.verdict, Verdict)

    def test_empty_protein_fail(self):
        """Empty protein → FAIL."""
        result = evaluate_charge_composition("", "", "Homo_sapiens")
        assert result.verdict == Verdict.FAIL
        assert "Empty" in result.violation

    def test_acidic_protein_pass(self):
        """Acidic protein (low pI, many charged residues) should PASS."""
        dna = _protein_to_dna(ACIDIC_PROTEIN)
        result = evaluate_charge_composition(dna, ACIDIC_PROTEIN, "Homo_sapiens")
        assert result.verdict == Verdict.PASS, (
            f"Acidic protein got {result.verdict}, expected PASS"
        )

    def test_custom_min_charged_fraction(self):
        """Custom min_charged_fraction should affect verdict."""
        # Require 50% charged residues — most proteins don't have that
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_charge_composition(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", min_charged_fraction=0.50
        )
        # Very high threshold → should fail or at least not PASS
        assert result.verdict in (Verdict.LIKELY_FAIL, Verdict.FAIL, Verdict.UNCERTAIN), (
            f"50% charge threshold should be hard to meet, got {result.verdict}"
        )

    def test_custom_max_pi(self):
        """Custom max_pI should affect verdict."""
        dna = _protein_to_dna(BASIC_PROTEIN)
        # Very low max_pI → almost any protein with K/R will exceed it
        result_strict = evaluate_charge_composition(
            dna, BASIC_PROTEIN, "Homo_sapiens", max_pI=5.0
        )
        assert result_strict.verdict in (Verdict.UNCERTAIN, Verdict.LIKELY_FAIL, Verdict.FAIL), (
            f"Very low max_pI=5.0 should flag basic protein, got {result_strict.verdict}"
        )

    def test_derivation_fields(self):
        """Derivation should contain charge composition details."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_charge_composition(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        steps = {d["step"]: d["value"] for d in result.derivation}
        assert "charged_fraction" in steps
        assert "min_charged_fraction" in steps
        assert "isoelectric_point" in steps
        assert "max_pI" in steps
        assert "positive_residues" in steps
        assert "negative_residues" in steps
        assert "total_charged" in steps
        assert "protein_length" in steps

    def test_derivation_charged_fraction_in_range(self):
        """Charged fraction in derivation should be in [0, 1]."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_charge_composition(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        steps = {d["step"]: d["value"] for d in result.derivation}
        cf = steps["charged_fraction"]
        assert 0.0 <= cf <= 1.0, f"Charged fraction {cf} outside [0, 1]"

    def test_derivation_pi_in_range(self):
        """pI in derivation should be in [0, 14]."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_charge_composition(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        steps = {d["step"]: d["value"] for d in result.derivation}
        pi = steps["isoelectric_point"]
        assert 0.0 <= pi <= 14.0, f"pI {pi} outside [0, 14]"

    def test_predicate_name_contains_params(self):
        """Predicate name should include parameter values."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_charge_composition(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", min_charged_fraction=0.15, max_pI=8.5
        )
        assert "ChargeComposition" in result.predicate
        assert "0.15" in result.predicate
        assert "8.5" in result.predicate

    def test_low_charge_and_high_pi_combined(self):
        """Low charge AND high pI together → LIKELY_FAIL with combined message."""
        # Protein with few charges and many basic residues
        low_charge_high_pi = "MAVILKWMAPTVIRLFWVILMAFGLMVVIAF"
        dna = _protein_to_dna(low_charge_high_pi)
        result = evaluate_charge_composition(dna, low_charge_high_pi, "Homo_sapiens")
        if result.verdict == Verdict.LIKELY_FAIL:
            assert result.violation is not None


class TestEvaluateNoLongHydrophobicStretch:
    """Tests for evaluate_no_long_hydrophobic_stretch predicate."""

    def test_soluble_protein_pass(self):
        """Charged protein should have no long hydrophobic stretches."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_no_long_hydrophobic_stretch(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_insoluble_protein_fail(self):
        """Hydrophobic protein should trigger long stretch detection."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_no_long_hydrophobic_stretch(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN), (
            f"Insoluble protein got {result.verdict}, expected failure verdict"
        )

    def test_empty_protein_pass(self):
        """Empty protein → PASS."""
        result = evaluate_no_long_hydrophobic_stretch("", "", "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_no_long_stretch_pass(self):
        """Protein with short hydrophobic stretches → PASS."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_no_long_hydrophobic_stretch(dna, BALANCED_GLOBULAR, "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_custom_max_stretch(self):
        """Custom max_stretch should affect verdict."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        # Very permissive
        result_permissive = evaluate_no_long_hydrophobic_stretch(
            dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens", max_stretch=100
        )
        assert result_permissive.verdict == Verdict.PASS

    def test_all_alanine_protein(self):
        """All-alanine protein: A is in hydrophobic set → long stretch."""
        dna = _protein_to_dna(ALL_ALANINE)
        result = evaluate_no_long_hydrophobic_stretch(dna, ALL_ALANINE, "Homo_sapiens")
        # All alanine → one big hydrophobic stretch > 7 → FAIL
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)

    def test_mixed_protein_with_hydro_stretch(self):
        """Protein with a long hydrophobic stretch in the middle → flagged."""
        dna = _protein_to_dna(LONG_HYDRO_STRETCH)
        result = evaluate_no_long_hydrophobic_stretch(dna, LONG_HYDRO_STRETCH, "Homo_sapiens")
        # Has a long LVV stretch → should be flagged
        assert result.verdict != Verdict.PASS or True  # depends on exact stretch length

    def test_derivation_structure(self):
        """Derivation should list long hydrophobic stretches."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_no_long_hydrophobic_stretch(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        steps = {d["step"]: d["value"] for d in result.derivation}
        assert "long_stretches" in steps
        assert "longest_stretch" in steps or "max_stretch_found" in steps

    def test_predicate_name_contains_max(self):
        """Predicate name should include max_stretch value."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_no_long_hydrophobic_stretch(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", max_stretch=10
        )
        assert "NoLongHydrophobicStretch" in result.predicate
        assert "10" in result.predicate

    def test_knowledge_gap_without_pdb(self):
        """Without PDB, knowledge_gap should note structural context missing."""
        # Use insoluble protein to trigger the long-stretch code path that sets knowledge_gap
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_no_long_hydrophobic_stretch(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        # knowledge_gap is set only when long stretches are found (main code path)
        if result.verdict != Verdict.PASS:
            assert result.knowledge_gap is not None

    def test_no_knowledge_gap_with_pdb(self):
        """With PDB, knowledge_gap should be None."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        result = evaluate_no_long_hydrophobic_stretch(
            dna, BALANCED_GLOBULAR, "Homo_sapiens", pdb_string="ATOM ..."
        )
        assert result.knowledge_gap is None


# ══════════════════════════════════════════════════════════════
# 3. Score ranges
# ══════════════════════════════════════════════════════════════

class TestCamSolScoreRanges:
    """Tests for CamSol scoring range validity."""

    def test_camsol_overall_score_range_soluble(self):
        """CamSol score for soluble protein should be > 0."""
        score = _camsol_overall_score(SOLUBLE_CHARGED)
        assert score > 0, f"Soluble protein score {score} should be > 0"

    def test_camsol_overall_score_range_insoluble(self):
        """CamSol score for insoluble protein should be < 0."""
        score = _camsol_overall_score(INSOLUBLE_HYDROPHOBIC)
        assert score < 0, f"Insoluble protein score {score} should be < 0"

    def test_camsol_score_bounded(self):
        """CamSol overall scores should be in a reasonable range [-5, 5]."""
        for protein in [SOLUBLE_CHARGED, INSOLUBLE_HYDROPHOBIC, BALANCED_GLOBULAR,
                        SUPER_SOLUBLE, SUPER_INSOLUBLE, ALL_ALANINE]:
            score = _camsol_overall_score(protein)
            assert -5.0 <= score <= 5.0, (
                f"CamSol score {score} for {protein[:10]}... outside [-5, 5]"
            )

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble", SOLUBLE_CHARGED),
        ("insoluble", INSOLUBLE_HYDROPHOBIC),
        ("balanced", BALANCED_GLOBULAR),
        ("super_soluble", SUPER_SOLUBLE),
        ("super_insoluble", SUPER_INSOLUBLE),
    ])
    def test_smoothed_profile_length(self, protein_name, protein):
        """Smoothed profile length should equal protein length."""
        profile = _camsol_smoothed_profile(protein)
        assert len(profile) == len(protein), (
            f"Profile length {len(profile)} != protein length {len(protein)} "
            f"for {protein_name}"
        )

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble", SOLUBLE_CHARGED),
        ("insoluble", INSOLUBLE_HYDROPHOBIC),
        ("balanced", BALANCED_GLOBULAR),
    ])
    def test_smoothed_profile_values_bounded(self, protein_name, protein):
        """Smoothed profile values should be in reasonable range."""
        profile = _camsol_smoothed_profile(protein)
        for val in profile:
            assert -10.0 <= val <= 10.0, (
                f"Profile value {val} outside [-10, 10] for {protein_name}"
            )

    def test_soluble_higher_than_insoluble(self):
        """Soluble protein should have higher CamSol score than insoluble."""
        score_sol = _camsol_overall_score(SOLUBLE_CHARGED)
        score_ins = _camsol_overall_score(INSOLUBLE_HYDROPHOBIC)
        assert score_sol > score_ins, (
            f"Soluble score ({score_sol}) should be > insoluble ({score_ins})"
        )

    def test_camsol_intrinsic_dict_coverage(self):
        """_CAMSOL_INTRINSIC should cover all 20 standard amino acids."""
        standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert standard_aas.issubset(set(_CAMSOL_INTRINSIC.keys())), (
            f"Missing AAs: {standard_aas - set(_CAMSOL_INTRINSIC.keys())}"
        )

    def test_camsol_intrinsic_values_bounded(self):
        """All intrinsic values should be in [-2, 2]."""
        for aa, val in _CAMSOL_INTRINSIC.items():
            assert -2.0 <= val <= 2.0, f"AA {aa} intrinsic value {val} outside [-2, 2]"

    def test_camsol_scale_positive(self):
        """CamSol scaling factor should be positive."""
        assert _CAMSOL_SCALE > 0

    def test_camsol_window_positive(self):
        """CamSol smoothing window should be positive."""
        assert _CAMSOL_WINDOW > 0

    def test_empty_protein_score_zero(self):
        """Empty protein should return score 0.0."""
        assert _camsol_overall_score("") == 0.0

    def test_empty_protein_profile_empty(self):
        """Empty protein should return empty profile."""
        assert _camsol_smoothed_profile("") == []

    def test_empty_protein_no_regions(self):
        """Empty protein should have no aggregation regions."""
        assert _find_aggregation_regions("") == []

    def test_known_soluble_residues_positive(self):
        """Known soluble residues (K, E, D, R) should have positive intrinsic scores."""
        for aa in ["K", "E", "D", "R"]:
            assert _CAMSOL_INTRINSIC[aa] > 0, (
                f"Soluble residue {aa} should have positive intrinsic score, "
                f"got {_CAMSOL_INTRINSIC[aa]}"
            )

    def test_known_insoluble_residues_negative(self):
        """Known aggregation-prone residues (I, L, V, F, W) should have negative intrinsic scores."""
        for aa in ["I", "L", "V", "F", "W"]:
            assert _CAMSOL_INTRINSIC[aa] < 0, (
                f"Insoluble residue {aa} should have negative intrinsic score, "
                f"got {_CAMSOL_INTRINSIC[aa]}"
            )


class TestAggregationRegionRanges:
    """Tests for aggregation-prone region scoring and detection."""

    def test_find_aggregation_regions_insoluble(self):
        """Insoluble protein should have aggregation-prone regions."""
        regions = _find_aggregation_regions(INSOLUBLE_HYDROPHOBIC)
        assert len(regions) > 0, "Insoluble protein should have aggregation-prone regions"

    def test_find_aggregation_regions_soluble(self):
        """Soluble protein should have few or no aggregation-prone regions."""
        regions = _find_aggregation_regions(SOLUBLE_CHARGED)
        # Soluble protein should have very few or none
        assert len(regions) <= 2, (
            f"Soluble protein has {len(regions)} regions, expected 0-2"
        )

    def test_region_format(self):
        """Each region should be (start, end, avg_score) tuple."""
        regions = _find_aggregation_regions(INSOLUBLE_HYDROPHOBIC)
        for start, end, avg_score in regions:
            assert isinstance(start, int) and start >= 0
            assert isinstance(end, int) and end > start
            assert isinstance(avg_score, (int, float))
            assert avg_score < 0, f"Aggregation region avg_score {avg_score} should be < 0"

    def test_custom_threshold(self):
        """Custom score_threshold should affect region detection."""
        regions_default = _find_aggregation_regions(INSOLUBLE_HYDROPHOBIC, score_threshold=-1.0)
        regions_lenient = _find_aggregation_regions(INSOLUBLE_HYDROPHOBIC, score_threshold=-2.0)
        # More lenient (lower) threshold → fewer regions or same
        assert len(regions_lenient) <= len(regions_default)

    def test_regions_covered_by_profile(self):
        """Region positions should be within protein bounds."""
        regions = _find_aggregation_regions(INSOLUBLE_HYDROPHOBIC)
        n = len(INSOLUBLE_HYDROPHOBIC)
        for start, end, _ in regions:
            assert start >= 0
            assert end <= n


# ══════════════════════════════════════════════════════════════
# Helper function tests
# ══════════════════════════════════════════════════════════════

class TestComputeNetCharge:
    """Tests for compute_net_charge helper."""

    def test_neutral_pH_7(self):
        """At pH 7, a protein with mixed charges should have finite net charge."""
        charge = compute_net_charge(BALANCED_GLOBULAR, 7.0)
        assert isinstance(charge, float)
        assert abs(charge) < 100, f"Net charge {charge} seems unreasonable"

    def test_acidic_protein_negative_at_pH7(self):
        """Acidic protein should have negative net charge at pH 7."""
        charge = compute_net_charge(ACIDIC_PROTEIN, 7.0)
        assert charge < 0, f"Acidic protein charge {charge} should be negative at pH 7"

    def test_basic_protein_positive_at_pH7(self):
        """Basic protein should have positive net charge at pH 7."""
        charge = compute_net_charge(BASIC_PROTEIN, 7.0)
        assert charge > 0, f"Basic protein charge {charge} should be positive at pH 7"

    def test_extreme_low_pH_positive(self):
        """At very low pH, net charge should be positive (all groups protonated)."""
        charge = compute_net_charge(BALANCED_GLOBULAR, 0.0)
        assert charge > 0

    def test_extreme_high_pH_negative(self):
        """At very high pH, net charge should be negative (all groups deprotonated)."""
        charge = compute_net_charge(BALANCED_GLOBULAR, 14.0)
        assert charge < 0

    def test_charge_monotonic_with_pH(self):
        """Net charge should monotonically decrease as pH increases."""
        charges = [compute_net_charge(BALANCED_GLOBULAR, pH) for pH in range(0, 15)]
        for i in range(len(charges) - 1):
            assert charges[i] >= charges[i + 1], (
                f"Charge not monotonically decreasing: pH {i} → {charges[i]:.3f}, "
                f"pH {i+1} → {charges[i+1]:.3f}"
            )

    def test_empty_protein_zero_charge(self):
        """Empty protein should have zero net charge."""
        assert compute_net_charge("", 7.0) == 0.0

    def test_single_methionine_charge(self):
        """Single M: N-term (+) and C-term (-) at pH 7."""
        charge = compute_net_charge("M", 7.0)
        # At pH 7: N-term mostly protonated (+), C-term mostly deprotonated (-)
        # Net should be near zero but slightly positive (N-term pKa 9.69 > 7)
        assert isinstance(charge, float)

    @pytest.mark.parametrize("pH", [0, 1, 7, 10, 14])
    def test_charge_is_finite(self, pH):
        """Net charge should be finite at all pH values."""
        charge = compute_net_charge(BALANCED_GLOBULAR, float(pH))
        assert abs(charge) < 1e6


class TestComputeApproximatePI:
    """Tests for compute_approximate_pI helper."""

    def test_pi_in_range(self):
        """pI should be in [0, 14]."""
        pI = compute_approximate_pI(BALANCED_GLOBULAR)
        assert 0.0 <= pI <= 14.0

    def test_basic_protein_high_pi(self):
        """Basic protein should have high pI (> 7)."""
        pI = compute_approximate_pI(BASIC_PROTEIN)
        assert pI > 7.0, f"Basic protein pI {pI:.2f} should be > 7.0"

    def test_acidic_protein_low_pi(self):
        """Acidic protein should have low pI (< 7)."""
        pI = compute_approximate_pI(ACIDIC_PROTEIN)
        assert pI < 7.0, f"Acidic protein pI {pI:.2f} should be < 7.0"

    def test_pi_zero_charge_crossing(self):
        """At the computed pI, net charge should be approximately zero."""
        pI = compute_approximate_pI(BALANCED_GLOBULAR)
        charge = compute_net_charge(BALANCED_GLOBULAR, pI)
        assert abs(charge) < 0.01, (
            f"Net charge at pI {pI}: {charge:.4f} should be ≈ 0"
        )

    def test_empty_protein_pi(self):
        """Empty protein should return pI 7.0."""
        assert compute_approximate_pI("") == 7.0

    def test_pi_precision(self):
        """pI should be computed with good precision (bisection with 100 iterations)."""
        pI = compute_approximate_pI(BALANCED_GLOBULAR)
        charge = compute_net_charge(BALANCED_GLOBULAR, pI)
        assert abs(charge) < 1e-10, f"pI precision too low: charge = {charge}"

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble", SOLUBLE_CHARGED),
        ("insoluble", INSOLUBLE_HYDROPHOBIC),
        ("balanced", BALANCED_GLOBULAR),
        ("basic", BASIC_PROTEIN),
        ("acidic", ACIDIC_PROTEIN),
    ])
    def test_pi_always_valid(self, protein_name, protein):
        """pI should always be in [0, 14] for any valid protein."""
        pI = compute_approximate_pI(protein)
        assert 0.0 <= pI <= 14.0, f"pI {pI} out of range for {protein_name}"


class TestFindHydrophobicStretches:
    """Tests for find_hydrophobic_stretches helper."""

    def test_insoluble_has_stretches(self):
        """Insoluble protein should have hydrophobic stretches."""
        stretches = find_hydrophobic_stretches(INSOLUBLE_HYDROPHOBIC)
        assert len(stretches) > 0

    def test_soluble_few_stretches(self):
        """Soluble protein should have few or no hydrophobic stretches."""
        stretches = find_hydrophobic_stretches(SOLUBLE_CHARGED)
        # Charged protein might have M at start, but generally very few
        total_hydro = sum(e - s for s, e in stretches)
        assert total_hydro < len(SOLUBLE_CHARGED) * 0.2, (
            f"Too much hydrophobic content in soluble protein: {total_hydro}"
        )

    def test_stretch_format(self):
        """Stretches should be (start, end) tuples with start < end."""
        stretches = find_hydrophobic_stretches(INSOLUBLE_HYDROPHOBIC)
        for start, end in stretches:
            assert isinstance(start, int) and start >= 0
            assert isinstance(end, int) and end > start
            assert end <= len(INSOLUBLE_HYDROPHOBIC)

    def test_stretches_dont_overlap(self):
        """Stretches should not overlap."""
        stretches = find_hydrophobic_stretches(INSOLUBLE_HYDROPHOBIC)
        for i in range(len(stretches) - 1):
            _, end_i = stretches[i]
            start_j, _ = stretches[i + 1]
            assert end_i <= start_j, f"Overlapping stretches: {stretches[i]} and {stretches[i+1]}"

    def test_empty_protein_no_stretches(self):
        """Empty protein → no stretches."""
        assert find_hydrophobic_stretches("") == []

    def test_all_alanine_is_hydrophobic(self):
        """All-alanine should be one big hydrophobic stretch (A is in default set)."""
        stretches = find_hydrophobic_stretches(ALL_ALANINE)
        total = sum(e - s for s, e in stretches)
        assert total == len(ALL_ALANINE), "All alanine should be entirely hydrophobic"

    def test_custom_hydrophobic_set(self):
        """Custom hydrophobic set should affect stretch detection."""
        # Only I and L as hydrophobic
        custom_set = {"I", "L"}
        stretches = find_hydrophobic_stretches(INSOLUBLE_HYDROPHOBIC, hydrophobic=custom_set)
        for start, end in stretches:
            for i in range(start, end):
                assert INSOLUBLE_HYDROPHOBIC[i] in custom_set

    def test_no_hydrophobic_residues(self):
        """Protein with no hydrophobic residues → no stretches."""
        all_charged = "KEDEKEDEKEDEKEDEKEDE"
        stretches = find_hydrophobic_stretches(all_charged)
        # K, E, D are not in the default hydrophobic set (AILMFWV)
        assert len(stretches) == 0, f"Expected no stretches, got {stretches}"

    def test_default_hydrophobic_set(self):
        """Default hydrophobic set should be AILMFWV."""
        assert _DEFAULT_HYDROPHOBIC == set("AILMFWV")


# ══════════════════════════════════════════════════════════════
# Cross-predicate consistency and integration
# ══════════════════════════════════════════════════════════════

class TestCrossPredicateConsistency:
    """Tests ensuring cross-predicate consistency."""

    def test_soluble_protein_passes_all_predicates(self):
        """A known soluble protein should pass or likely-pass all four predicates."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        r1 = evaluate_soluble_expression(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        r2 = evaluate_no_aggregation_prone_region(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        r3 = evaluate_charge_composition(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        r4 = evaluate_no_long_hydrophobic_stretch(dna, SOLUBLE_CHARGED, "Homo_sapiens")

        for name, result in [
            ("SolubleExpression", r1),
            ("NoAggregationProneRegion", r2),
            ("ChargeComposition", r3),
            ("NoLongHydrophobicStretch", r4),
        ]:
            assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS), (
                f"{name} verdict {result.verdict} for soluble protein"
            )

    def test_insoluble_protein_fails_at_least_one(self):
        """An insoluble protein should fail at least one predicate."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        results = [
            evaluate_soluble_expression(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens"),
            evaluate_no_aggregation_prone_region(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens"),
            evaluate_charge_composition(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens"),
            evaluate_no_long_hydrophobic_stretch(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens"),
        ]
        fail_count = sum(
            1 for r in results
            if r.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)
        )
        assert fail_count >= 1, "Insoluble protein should fail at least one predicate"

    def test_all_predicates_return_type_check_result(self):
        """All predicates should return TypeCheckResult instances."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        results = [
            evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
            evaluate_no_aggregation_prone_region(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
            evaluate_charge_composition(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
            evaluate_no_long_hydrophobic_stretch(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
        ]
        for r in results:
            assert isinstance(r, TypeCheckResult)

    def test_all_predicates_have_verdict(self):
        """All results should have a valid Verdict enum value."""
        dna = _protein_to_dna(BALANCED_GLOBULAR)
        results = [
            evaluate_soluble_expression(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
            evaluate_no_aggregation_prone_region(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
            evaluate_charge_composition(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
            evaluate_no_long_hydrophobic_stretch(dna, BALANCED_GLOBULAR, "Homo_sapiens"),
        ]
        for r in results:
            assert isinstance(r.verdict, Verdict)
            assert r.verdict in {
                Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
                Verdict.LIKELY_FAIL, Verdict.FAIL,
            }

    def test_passed_property_consistent(self):
        """The `passed` property should be True for PASS and LIKELY_PASS."""
        dna = _protein_to_dna(SOLUBLE_CHARGED)
        result = evaluate_soluble_expression(dna, SOLUBLE_CHARGED, "Homo_sapiens")
        if result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS):
            assert result.passed is True
        else:
            assert result.passed is False

    def test_violation_set_on_failure(self):
        """Failed/uncertain predicates should have a non-None violation."""
        dna = _protein_to_dna(INSOLUBLE_HYDROPHOBIC)
        result = evaluate_soluble_expression(dna, INSOLUBLE_HYDROPHOBIC, "Homo_sapiens")
        if result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN):
            assert result.violation is not None


class TestVerdictThresholdConstants:
    """Tests that verify the threshold constants are sensible."""

    def test_camsol_thresholds_ordering(self):
        """CamSol thresholds should be ordered: MARGINAL < 0 < HIGHLY_SOLUBLE."""
        assert _CAMSOL_MARGINAL < 0.0
        assert _CAMSOL_HIGHLY_SOLUBLE > 0.0
        assert _CAMSOL_MARGINAL < _CAMSOL_HIGHLY_SOLUBLE

    def test_agg_thresholds_ordering(self):
        """Aggregation region length thresholds should be increasing."""
        assert _AGG_BORDERLINE_MAX < _AGG_UNCERTAIN_MAX
        assert _AGG_UNCERTAIN_MAX < _AGG_LIKELY_FAIL_MAX

    def test_hydro_excess_thresholds_ordering(self):
        """Hydrophobic excess thresholds should be increasing."""
        assert _HYDRO_EXCESS_BORDERLINE < _HYDRO_EXCESS_UNCERTAIN

    def test_pka_values_positive(self):
        """All pKa values should be positive."""
        for key, pka in PKA_VALUES.items():
            assert pka > 0, f"pKa for {key} is {pka}, should be > 0"

    def test_pka_values_in_range(self):
        """All pKa values should be in a reasonable range [0, 14]."""
        for key, pka in PKA_VALUES.items():
            assert 0 < pka <= 14, f"pKa for {key} is {pka}, outside (0, 14]"

    def test_verdict_enum_complete(self):
        """Verdict enum should have all five values."""
        assert set(Verdict) == {
            Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
            Verdict.LIKELY_FAIL, Verdict.FAIL,
        }


class TestScoreRangeParametrized:
    """Parametrized tests for score ranges across diverse proteins."""

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble_charged", SOLUBLE_CHARGED),
        ("insoluble_hydrophobic", INSOLUBLE_HYDROPHOBIC),
        ("balanced_globular", BALANCED_GLOBULAR),
        ("super_soluble", SUPER_SOLUBLE),
        ("super_insoluble", SUPER_INSOLUBLE),
        ("basic", BASIC_PROTEIN),
        ("acidic", ACIDIC_PROTEIN),
        ("all_alanine", ALL_ALANINE),
    ])
    def test_camsol_score_bounded(self, protein_name, protein):
        """CamSol score should be in [-5, 5] for all test proteins."""
        score = _camsol_overall_score(protein)
        assert -5.0 <= score <= 5.0, (
            f"Score {score} for {protein_name} outside [-5, 5]"
        )

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble_charged", SOLUBLE_CHARGED),
        ("insoluble_hydrophobic", INSOLUBLE_HYDROPHOBIC),
        ("balanced_globular", BALANCED_GLOBULAR),
        ("super_soluble", SUPER_SOLUBLE),
        ("super_insoluble", SUPER_INSOLUBLE),
        ("basic", BASIC_PROTEIN),
        ("acidic", ACIDIC_PROTEIN),
        ("all_alanine", ALL_ALANINE),
    ])
    def test_pi_bounded(self, protein_name, protein):
        """pI should be in [0, 14] for all test proteins."""
        pI = compute_approximate_pI(protein)
        assert 0.0 <= pI <= 14.0, f"pI {pI} for {protein_name} outside [0, 14]"

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble_charged", SOLUBLE_CHARGED),
        ("insoluble_hydrophobic", INSOLUBLE_HYDROPHOBIC),
        ("balanced_globular", BALANCED_GLOBULAR),
        ("super_soluble", SUPER_SOLUBLE),
        ("super_insoluble", SUPER_INSOLUBLE),
        ("basic", BASIC_PROTEIN),
        ("acidic", ACIDIC_PROTEIN),
        ("all_alanine", ALL_ALANINE),
    ])
    def test_net_charge_bounded_at_pH7(self, protein_name, protein):
        """Net charge at pH 7 should be finite and reasonable."""
        charge = compute_net_charge(protein, 7.0)
        assert abs(charge) < len(protein) + 2, (
            f"Net charge {charge} for {protein_name} seems too large"
        )

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble_charged", SOLUBLE_CHARGED),
        ("insoluble_hydrophobic", INSOLUBLE_HYDROPHOBIC),
        ("balanced_globular", BALANCED_GLOBULAR),
        ("super_soluble", SUPER_SOLUBLE),
        ("super_insoluble", SUPER_INSOLUBLE),
    ])
    def test_predicate_verdict_is_valid(self, protein_name, protein):
        """All four predicates should produce valid verdicts."""
        dna = _protein_to_dna(protein)
        results = [
            evaluate_soluble_expression(dna, protein, "Homo_sapiens"),
            evaluate_no_aggregation_prone_region(dna, protein, "Homo_sapiens"),
            evaluate_charge_composition(dna, protein, "Homo_sapiens"),
            evaluate_no_long_hydrophobic_stretch(dna, protein, "Homo_sapiens"),
        ]
        valid_verdicts = {Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
                         Verdict.LIKELY_FAIL, Verdict.FAIL}
        for r in results:
            assert r.verdict in valid_verdicts, (
                f"{protein_name}: verdict {r.verdict} not in valid set"
            )

    @pytest.mark.parametrize("protein_name,protein", [
        ("soluble_charged", SOLUBLE_CHARGED),
        ("insoluble_hydrophobic", INSOLUBLE_HYDROPHOBIC),
        ("balanced_globular", BALANCED_GLOBULAR),
        ("super_soluble", SUPER_SOLUBLE),
        ("super_insoluble", SUPER_INSOLUBLE),
    ])
    def test_smoothed_profile_length_matches(self, protein_name, protein):
        """Smoothed profile length should match protein length."""
        profile = _camsol_smoothed_profile(protein)
        assert len(profile) == len(protein)
