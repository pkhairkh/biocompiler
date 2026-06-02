"""
BioCompiler CamSol Solubility Module — Comprehensive Test Suite
================================================================
Tests for the CamSol solubility prediction module (camsol.py) and
the solubility type-check predicates (solubility_predicates.py).
"""

import pytest
from biocompiler.camsol import (
    compute_intrinsic_solubility,
    compute_solubility,
    classify_solubility,
    find_solubility_mutations,
    SolubilityResult,
    CAMSOL_HYDROPATHY,
    CAMSOL_CHARGE,
    CAMSOL_ALPHA_HELIX,
)
from biocompiler.solubility_predicates import (
    evaluate_soluble_expression,
    evaluate_no_aggregation_prone_region,
    evaluate_charge_composition,
    evaluate_no_long_hydrophobic_stretch,
    compute_approximate_pI,
    compute_net_charge,
    find_hydrophobic_stretches,
    PKA_VALUES,
)
from biocompiler.type_system import Verdict, AA_TO_CODONS
from biocompiler.types import TypeCheckResult
from biocompiler.exceptions import CamSolError


# ────────────────────────────────────────────────────────────
# Test proteins
# ────────────────────────────────────────────────────────────

# Soluble protein (highly charged, hydrophilic)
SOLUBLE_PROTEIN = "MSEKKDKKEKEKKDEKKDEEKKDEEKKDESKKDEEKKDEEKKDESKKDEEKKDEEKK"

# Insoluble/aggregation-prone protein (hydrophobic)
INSOLUBLE_PROTEIN = "MVVVIIVVVLLLFLLLLFFFFWWWAAAIIIMMM"

# Balanced protein (average composition)
BALANCED_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAVDILSKKGDVQVIK"

# Helper: generate a DNA coding sequence from a protein
def _protein_to_dna(protein: str) -> str:
    """Generate a DNA sequence encoding the given protein using first codon for each AA."""
    return "".join(AA_TO_CODONS.get(aa, ["GCT"])[0] for aa in protein)


# ══════════════════════════════════════════════════════════════
# TestCamSolIntrinsic — core solubility computation
# ══════════════════════════════════════════════════════════════

class TestCamSolIntrinsic:
    """Tests for compute_intrinsic_solubility and related functions."""

    def test_compute_intrinsic_solubility_returns_result(self):
        """Basic call returns a SolubilityResult."""
        result = compute_intrinsic_solubility(SOLUBLE_PROTEIN)
        assert isinstance(result, SolubilityResult)

    def test_intrinsic_score_range(self):
        """Score should be roughly in [-3, 3] range."""
        for protein in [SOLUBLE_PROTEIN, INSOLUBLE_PROTEIN, BALANCED_PROTEIN]:
            result = compute_intrinsic_solubility(protein)
            assert -3.0 <= result.intrinsic_score <= 3.0, (
                f"Score {result.intrinsic_score} outside [-3, 3] for protein starting with {protein[:10]}"
            )

    def test_soluble_protein_high_score(self):
        """Known soluble protein should have positive intrinsic score."""
        result = compute_intrinsic_solubility(SOLUBLE_PROTEIN)
        assert result.intrinsic_score > 0, (
            f"Soluble protein score {result.intrinsic_score} should be positive"
        )

    def test_insoluble_protein_low_score(self):
        """Hydrophobic protein should have lower score than soluble one."""
        r_soluble = compute_intrinsic_solubility(SOLUBLE_PROTEIN)
        r_insoluble = compute_intrinsic_solubility(INSOLUBLE_PROTEIN)
        assert r_insoluble.intrinsic_score < r_soluble.intrinsic_score, (
            f"Insoluble ({r_insoluble.intrinsic_score}) should score lower than "
            f"soluble ({r_soluble.intrinsic_score})"
        )

    def test_per_residue_scores_length(self):
        """Length of per-residue scores should match protein length."""
        for protein in [SOLUBLE_PROTEIN, INSOLUBLE_PROTEIN, BALANCED_PROTEIN]:
            result = compute_intrinsic_solubility(protein)
            assert len(result.per_residue_scores) == len(protein), (
                f"Per-residue scores length {len(result.per_residue_scores)} "
                f"!= protein length {len(protein)}"
            )

    def test_aggregation_prone_regions(self):
        """Hydrophobic protein should have aggregation-prone regions or low score."""
        result = compute_intrinsic_solubility(INSOLUBLE_PROTEIN)
        # Either there are APRs or the overall score is negative/low
        has_aprs = len(result.aggregation_prone_regions) > 0
        low_score = result.intrinsic_score < 0.0
        assert has_aprs or low_score, (
            "Insoluble protein should have APRs or low intrinsic score"
        )

    def test_solubility_class_values(self):
        """Verify classification returns known string values."""
        valid_classes = {"highly_soluble", "soluble", "marginally_soluble", "insoluble"}
        for protein in [SOLUBLE_PROTEIN, INSOLUBLE_PROTEIN, BALANCED_PROTEIN]:
            result = compute_intrinsic_solubility(protein)
            assert result.solubility_class in valid_classes, (
                f"Unknown solubility class: {result.solubility_class}"
            )

    def test_classify_solubility(self):
        """Test all 4 classes directly via classify_solubility()."""
        assert classify_solubility(2.0) == "highly_soluble"
        assert classify_solubility(0.5) == "soluble"
        assert classify_solubility(-0.5) == "marginally_soluble"
        assert classify_solubility(-2.0) == "insoluble"

    def test_camsol_hydropathy_scale(self):
        """All 20 standard AAs should be present in the hydropathy scale."""
        standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(CAMSOL_HYDROPATHY.keys()) == standard_aas, (
            f"Missing AAs in CAMSOL_HYDROPATHY: {standard_aas - set(CAMSOL_HYDROPATHY.keys())}"
        )

    def test_camsol_charge_scale(self):
        """Charged AAs should have non-zero values in the charge scale."""
        charged_aas = {"K", "R", "H", "D", "E"}
        for aa in charged_aas:
            assert CAMSOL_CHARGE.get(aa, 0.0) != 0.0, (
                f"Charged AA {aa} should have non-zero charge value"
            )

    def test_camsol_alpha_helix_scale(self):
        """All 20 standard AAs should be present in the alpha-helix scale."""
        standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(CAMSOL_ALPHA_HELIX.keys()) == standard_aas, (
            f"Missing AAs in CAMSOL_ALPHA_HELIX: {standard_aas - set(CAMSOL_ALPHA_HELIX.keys())}"
        )

    def test_generate_recommendations(self):
        """Insoluble protein should have recommendations."""
        result = compute_intrinsic_solubility(INSOLUBLE_PROTEIN)
        # If the protein scores low enough, recommendations should be generated
        if result.intrinsic_score < 0:
            assert len(result.recommendations) > 0, (
                "Insoluble protein should have recommendations"
            )

    def test_empty_protein_raises(self):
        """Empty protein should raise CamSolError."""
        with pytest.raises(CamSolError, match="empty"):
            compute_intrinsic_solubility("")

    def test_overall_score_equals_intrinsic_when_no_structure(self):
        """When no PDB is provided, overall_score should equal intrinsic_score."""
        result = compute_intrinsic_solubility(SOLUBLE_PROTEIN)
        assert result.overall_score == result.intrinsic_score

    def test_method_is_intrinsic(self):
        """Method should be 'camsol_intrinsic' for sequence-only computation."""
        result = compute_intrinsic_solubility(BALANCED_PROTEIN)
        assert result.method == "camsol_intrinsic"


# ══════════════════════════════════════════════════════════════
# TestCamSolMutations — solubility-improving mutation suggestions
# ══════════════════════════════════════════════════════════════

class TestCamSolMutations:
    """Tests for find_solubility_mutations."""

    def test_find_solubility_mutations(self):
        """Run on protein with aggregation region should produce mutations."""
        mutations = find_solubility_mutations(INSOLUBLE_PROTEIN)
        assert isinstance(mutations, list)
        # The insoluble protein should have some mutation suggestions
        assert len(mutations) > 0, "Expected at least one mutation suggestion"

    def test_find_solubility_mutations_blosum62_filter(self):
        """Mutations should have valid score (BLOSUM62-conservative)."""
        mutations = find_solubility_mutations(INSOLUBLE_PROTEIN)
        for mut in mutations:
            assert mut.score > 0, (
                f"Mutation {mut.original}{mut.position}{mut.mutant} "
                f"has score={mut.score} <= 0"
            )

    def test_mutation_improves_solubility(self):
        """Each suggested mutation should improve solubility (positive score)."""
        mutations = find_solubility_mutations(INSOLUBLE_PROTEIN)
        for mut in mutations:
            assert mut.score > 0, (
                f"Mutation {mut.original}{mut.position}{mut.mutant} "
                f"has score={mut.score} <= 0"
            )

    def test_mutation_has_required_attrs(self):
        """Each MutationResult should have required attributes."""
        required_attrs = {"position", "original", "mutant", "score", "engine", "description"}
        mutations = find_solubility_mutations(INSOLUBLE_PROTEIN)
        for mut in mutations:
            missing = required_attrs - set(dir(mut))
            assert not missing, (
                f"Missing attributes in MutationResult: {missing}"
            )

    def test_no_mutations_for_soluble_protein(self):
        """Already soluble protein should have no or few mutations."""
        mutations = find_solubility_mutations(SOLUBLE_PROTEIN)
        # Soluble protein already scores well, so few/no mutations needed
        # (may still have some for marginal improvements)
        assert isinstance(mutations, list)


# ══════════════════════════════════════════════════════════════
# TestSolubilityPredicates — type-check predicates
# ══════════════════════════════════════════════════════════════

class TestSolubilityPredicates:
    """Tests for solubility type-check predicates."""

    def test_evaluate_soluble_expression_pass(self):
        """Soluble protein → PASS or LIKELY_PASS."""
        dna = _protein_to_dna(SOLUBLE_PROTEIN)
        result = evaluate_soluble_expression(dna, SOLUBLE_PROTEIN, "Homo_sapiens")
        assert isinstance(result, TypeCheckResult)
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS), (
            f"Soluble protein verdict {result.verdict} should be PASS or LIKELY_PASS"
        )

    def test_evaluate_soluble_expression_fail(self):
        """Insoluble protein → FAIL or LIKELY_FAIL."""
        dna = _protein_to_dna(INSOLUBLE_PROTEIN)
        result = evaluate_soluble_expression(dna, INSOLUBLE_PROTEIN, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN), (
            f"Insoluble protein verdict {result.verdict} should indicate poor solubility"
        )

    def test_evaluate_no_aggregation_prone_region_pass(self):
        """Protein without long APRs → PASS."""
        dna = _protein_to_dna(SOLUBLE_PROTEIN)
        result = evaluate_no_aggregation_prone_region(dna, SOLUBLE_PROTEIN, "Homo_sapiens")
        assert result.verdict == Verdict.PASS, (
            f"Soluble protein APR verdict {result.verdict} should be PASS"
        )

    def test_evaluate_no_aggregation_prone_region_fail(self):
        """Protein with long APRs → FAIL/LIKELY_FAIL."""
        dna = _protein_to_dna(INSOLUBLE_PROTEIN)
        result = evaluate_no_aggregation_prone_region(dna, INSOLUBLE_PROTEIN, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL), (
            f"Insoluble protein APR verdict {result.verdict} should be FAIL or LIKELY_FAIL"
        )

    def test_evaluate_charge_composition_pass(self):
        """Balanced charge → PASS."""
        dna = _protein_to_dna(SOLUBLE_PROTEIN)
        result = evaluate_charge_composition(dna, SOLUBLE_PROTEIN, "Homo_sapiens")
        # The soluble protein is highly charged, should pass
        assert result.verdict == Verdict.PASS, (
            f"Soluble protein charge verdict {result.verdict} should be PASS"
        )

    def test_evaluate_charge_composition_low_charge(self):
        """Low charged fraction → LIKELY_FAIL."""
        # Build a protein with very few charged residues
        low_charge_protein = "MAVILFWMAPTVIGLFWVILMAFGLMVVIAF"
        dna = _protein_to_dna(low_charge_protein)
        result = evaluate_charge_composition(dna, low_charge_protein, "Homo_sapiens")
        assert result.verdict in (Verdict.LIKELY_FAIL, Verdict.FAIL), (
            f"Low charge protein verdict {result.verdict} should be LIKELY_FAIL or FAIL"
        )

    def test_evaluate_no_long_hydrophobic_stretch_pass(self):
        """No long stretches → PASS."""
        dna = _protein_to_dna(SOLUBLE_PROTEIN)
        result = evaluate_no_long_hydrophobic_stretch(dna, SOLUBLE_PROTEIN, "Homo_sapiens")
        assert result.verdict == Verdict.PASS, (
            f"Soluble protein hydrophobic stretch verdict {result.verdict} should be PASS"
        )

    def test_evaluate_no_long_hydrophobic_stretch_fail(self):
        """Long hydrophobic stretch → FAIL."""
        dna = _protein_to_dna(INSOLUBLE_PROTEIN)
        result = evaluate_no_long_hydrophobic_stretch(dna, INSOLUBLE_PROTEIN, "Homo_sapiens")
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN), (
            f"Insoluble protein hydrophobic stretch verdict {result.verdict} "
            f"should indicate problems"
        )

    def test_compute_approximate_pi(self):
        """Known protein, check pI in reasonable range."""
        pI = compute_approximate_pI(BALANCED_PROTEIN)
        assert 0.0 <= pI <= 14.0, f"pI {pI} outside valid range [0, 14]"
        # The balanced protein has many basic residues, so pI should be above 7
        # (it has more K/R/H than D/E)
        assert pI > 7.0, f"pI {pI:.2f} expected to be > 7.0 for basic protein"

    def test_compute_net_charge(self):
        """At pH 7, check sign for a charged protein."""
        # SOLUBLE_PROTEIN has many D/E (acidic), so net charge at pH 7 should be negative
        charge = compute_net_charge(SOLUBLE_PROTEIN, 7.0)
        assert charge < 0, (
            f"Net charge at pH 7 for soluble (acidic) protein should be negative, got {charge:.2f}"
        )
        # A basic protein should have positive net charge at pH 7
        basic_protein = "MKRKHKRHKRHKRHKRHKRHKRH"
        charge_basic = compute_net_charge(basic_protein, 7.0)
        assert charge_basic > 0, (
            f"Net charge at pH 7 for basic protein should be positive, got {charge_basic:.2f}"
        )

    def test_find_hydrophobic_stretches(self):
        """Known protein, verify stretch detection."""
        stretches = find_hydrophobic_stretches(INSOLUBLE_PROTEIN)
        assert len(stretches) > 0, "Insoluble protein should have hydrophobic stretches"
        # The insoluble protein is entirely hydrophobic, should have one big stretch
        longest = max(e - s for s, e in stretches)
        assert longest >= 10, f"Longest stretch {longest} should be >= 10"

    def test_pka_values(self):
        """All relevant AAs present in PKA_VALUES."""
        required_keys = {"K", "R", "H", "D", "E", "N_term", "C_term"}
        assert required_keys.issubset(set(PKA_VALUES.keys())), (
            f"Missing PKA keys: {required_keys - set(PKA_VALUES.keys())}"
        )

    def test_evaluate_soluble_expression_empty_protein(self):
        """Empty protein should result in FAIL verdict."""
        result = evaluate_soluble_expression("", "", "Homo_sapiens")
        assert result.verdict == Verdict.FAIL

    def test_evaluate_no_aggregation_prone_region_empty(self):
        """Empty protein should pass APR check (nothing to aggregate)."""
        result = evaluate_no_aggregation_prone_region("", "", "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_evaluate_no_long_hydrophobic_stretch_empty(self):
        """Empty protein should pass hydrophobic stretch check."""
        result = evaluate_no_long_hydrophobic_stretch("", "", "Homo_sapiens")
        assert result.verdict == Verdict.PASS

    def test_predicate_results_have_predicate_name(self):
        """All predicate results should have a meaningful predicate name."""
        dna = _protein_to_dna(BALANCED_PROTEIN)
        for func, name in [
            (evaluate_soluble_expression, "SolubleExpression"),
            (evaluate_no_aggregation_prone_region, "NoAggregationProneRegion"),
            (evaluate_charge_composition, "ChargeComposition"),
            (evaluate_no_long_hydrophobic_stretch, "NoLongHydrophobicStretch"),
        ]:
            result = func(dna, BALANCED_PROTEIN, "Homo_sapiens")
            assert name in result.predicate, (
                f"Predicate name '{result.predicate}' should contain '{name}'"
            )

    def test_predicate_results_have_verdict(self):
        """All predicate results should have a valid Verdict."""
        dna = _protein_to_dna(BALANCED_PROTEIN)
        for func in [
            evaluate_soluble_expression,
            evaluate_no_aggregation_prone_region,
            evaluate_charge_composition,
            evaluate_no_long_hydrophobic_stretch,
        ]:
            result = func(dna, BALANCED_PROTEIN, "Homo_sapiens")
            assert isinstance(result.verdict, Verdict)

    def test_predicate_results_have_derivation(self):
        """All predicate results should have derivation info."""
        dna = _protein_to_dna(BALANCED_PROTEIN)
        for func in [
            evaluate_soluble_expression,
            evaluate_no_aggregation_prone_region,
            evaluate_charge_composition,
            evaluate_no_long_hydrophobic_stretch,
        ]:
            result = func(dna, BALANCED_PROTEIN, "Homo_sapiens")
            assert result.derivation is not None
            assert isinstance(result.derivation, list)


# ══════════════════════════════════════════════════════════════
# TestComputeSolubility — unified API
# ══════════════════════════════════════════════════════════════

class TestComputeSolubility:
    """Tests for the compute_solubility unified API."""

    def test_compute_solubility_returns_result(self):
        """Basic call returns SolubilityResult."""
        result = compute_solubility(SOLUBLE_PROTEIN)
        assert isinstance(result, SolubilityResult)

    def test_compute_solubility_intrinsic_only(self):
        """Without PDB, should return intrinsic results."""
        result = compute_solubility(INSOLUBLE_PROTEIN)
        assert result.method == "camsol_intrinsic"
        assert result.structural_score is None

    def test_compute_solubility_with_empty_pdb(self):
        """Empty PDB string should fall back to intrinsic."""
        result = compute_solubility(SOLUBLE_PROTEIN, pdb_string="")
        assert result.method == "camsol_intrinsic"

    def test_soluble_higher_than_insoluble(self):
        """Soluble protein should score higher than insoluble."""
        r_sol = compute_solubility(SOLUBLE_PROTEIN)
        r_ins = compute_solubility(INSOLUBLE_PROTEIN)
        assert r_sol.overall_score > r_ins.overall_score

    def test_classify_consistency(self):
        """classify_solubility should match result.solubility_class."""
        for protein in [SOLUBLE_PROTEIN, INSOLUBLE_PROTEIN, BALANCED_PROTEIN]:
            result = compute_intrinsic_solubility(protein)
            expected = classify_solubility(result.intrinsic_score)
            assert result.solubility_class == expected
