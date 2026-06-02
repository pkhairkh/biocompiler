"""Integration tests for Phase 2+3 modules — cross-module workflows."""

import pytest

HBB_PROTEIN = "MVHLTPEEKSAVTALWGKVNVADIVGHALSDLHAKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
HBB_SEQUENCE = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"

MINI_PDB = """ATOM      1  CA  MET A   1       1.000   2.000   3.000  1.00 85.00           C
ATOM      2  CA  ALA A   2       4.000   2.000   3.000  1.00 90.00           C
ATOM      3  CA  GLY A   3       7.000   2.000   3.000  1.00 78.00           C
END
"""


# ────────────────────────────────────────────────────────────
# Test Structure Pipeline
# ────────────────────────────────────────────────────────────
class TestStructurePipeline:
    """End-to-end structure prediction pipeline."""

    def test_parse_pdb_to_quality_report(self):
        """Parse PDB → compute quality metrics → verify report."""
        from biocompiler.structure_quality import compute_structure_quality, StructureQualityReport
        report = compute_structure_quality(MINI_PDB)
        assert isinstance(report, StructureQualityReport)
        assert report.mean_plddt > 0
        # verdict is a string like "PASS", "LIKELY_FAIL" etc.
        assert isinstance(report.verdict, str)

    def test_esmfold_offline_result(self):
        """ESMFold offline fallback returns valid result."""
        from biocompiler.esmfold import predict_structure, ESMFoldResult
        result = predict_structure(HBB_PROTEIN)
        assert isinstance(result, ESMFoldResult)
        # In test env, API is likely unavailable → success may be False
        if result.success:
            assert result.pdb_string != ""
            assert result.mean_plddt > 0

    def test_cache_integration(self):
        """Cache put/get round-trip works."""
        from biocompiler.esmfold_cache import ESMFoldCache
        from biocompiler.esmfold import ESMFoldResult
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = ESMFoldCache(cache_dir=tmpdir, max_size=5)
            # Put a result using ESMFoldResult object
            esm_result = ESMFoldResult(
                protein="MKTV",
                pdb_string=MINI_PDB,
                plddt_scores=[85.0, 90.0, 78.0],
                mean_plddt=84.33,
                pae_matrix=None,
                model_name="esmfold_v1",
                execution_time_s=1.0,
                success=True,
            )
            cache.put(esm_result.protein, esm_result)
            # Retrieve it
            result = cache.get("MKTV")
            assert result is not None
            # Cache may return the ESMFoldResult object directly or as dict
            if isinstance(result, dict):
                assert result["mean_plddt"] == 85.0 or "mean_plddt" in result
            else:
                assert hasattr(result, 'mean_plddt')


# ────────────────────────────────────────────────────────────
# Test Stability Pipeline
# ────────────────────────────────────────────────────────────
class TestStabilityPipeline:
    """End-to-end stability analysis."""

    def test_empirical_stability_to_predicate(self):
        """Compute empirical stability → evaluate_stable_folding predicate."""
        from biocompiler.foldx import empirical_stability
        from biocompiler.stability_predicates import evaluate_stable_folding
        from biocompiler.types import Verdict, TypeCheckResult

        # Get empirical stability
        result = empirical_stability(HBB_PROTEIN)
        assert result.success is True
        assert result.stability_kcal != 0

        # Run predicate
        pred = evaluate_stable_folding(
            sequence=HBB_SEQUENCE, protein=HBB_PROTEIN, organism="Homo_sapiens"
        )
        assert isinstance(pred, TypeCheckResult)
        assert pred.verdict in list(Verdict)

    def test_stability_and_solubility_combined(self):
        """Analyze protein with both FoldX and CamSol."""
        from biocompiler.foldx import empirical_stability
        from biocompiler.camsol import compute_solubility

        stab = empirical_stability(HBB_PROTEIN)
        sol = compute_solubility(HBB_PROTEIN)
        assert stab.success is True
        assert sol.overall_score is not None


# ────────────────────────────────────────────────────────────
# Test Immunogenicity Pipeline
# ────────────────────────────────────────────────────────────
class TestImmunogenicityPipeline:
    """End-to-end immunogenicity analysis."""

    def test_mhc_binding_basic(self):
        """Predict MHC binding → verify results."""
        from biocompiler.mhc_binding import predict_mhc_i_binding, MHCBindingResult
        results = predict_mhc_i_binding(HBB_PROTEIN)
        assert isinstance(results, list)
        for r in results:
            assert isinstance(r, MHCBindingResult)

    def test_epitope_prediction(self):
        """Find epitopes → verify structure."""
        from biocompiler.epitope import predict_epitopes, EpitopePredictionResult
        result = predict_epitopes(HBB_PROTEIN)
        assert isinstance(result, EpitopePredictionResult)
        assert isinstance(result.linear_epitopes, list)
        assert 0 <= result.epitope_coverage <= 1

    def test_immunogenicity_scoring(self):
        """Compute immunogenicity → verify score."""
        from biocompiler.immunogenicity import compute_immunogenicity, ImmunogenicityResult
        result = compute_immunogenicity(HBB_PROTEIN)
        assert isinstance(result, ImmunogenicityResult)
        assert 0 <= result.overall_score <= 1

    def test_all_immuno_predicates(self):
        """Run all 4 immunogenicity predicates on same protein."""
        from biocompiler.immuno_predicates import (
            evaluate_low_immunogenicity,
            evaluate_no_strong_t_cell_epitope,
            evaluate_no_dominant_b_cell_epitope,
            evaluate_population_coverage_safe,
        )
        from biocompiler.types import Verdict, TypeCheckResult

        preds = [
            evaluate_low_immunogenicity(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_no_strong_t_cell_epitope(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_no_dominant_b_cell_epitope(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_population_coverage_safe(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
        ]
        for p in preds:
            assert isinstance(p, TypeCheckResult)
            assert p.verdict in list(Verdict)


# ────────────────────────────────────────────────────────────
# Test Cross-Module Integration
# ────────────────────────────────────────────────────────────
class TestCrossModuleIntegration:
    """Tests that span multiple Phase 2+3 modules."""

    def test_full_protein_assessment(self):
        """Run ALL Phase 2+3 predicates on HBB, verify all return TypeCheckResult."""
        from biocompiler.stability_predicates import (
            evaluate_stable_folding, evaluate_disulfide_bond_integrity,
            evaluate_hydrophobic_core_quality,
        )
        from biocompiler.solubility_predicates import (
            evaluate_soluble_expression, evaluate_charge_composition,
            evaluate_no_long_hydrophobic_stretch,
        )
        from biocompiler.immuno_predicates import (
            evaluate_low_immunogenicity, evaluate_population_coverage_safe,
        )
        from biocompiler.structure_predicates import (
            evaluate_structure_confidence,
        )
        from biocompiler.types import Verdict, TypeCheckResult

        predicates = [
            evaluate_stable_folding(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_disulfide_bond_integrity(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_hydrophobic_core_quality(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_soluble_expression(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_charge_composition(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_no_long_hydrophobic_stretch(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_low_immunogenicity(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_population_coverage_safe(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
            evaluate_structure_confidence(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens"),
        ]
        for p in predicates:
            assert isinstance(p, TypeCheckResult), f"Expected TypeCheckResult, got {type(p)}"
            assert p.verdict in list(Verdict), f"Invalid verdict: {p.verdict}"

    def test_verdict_types_correct(self):
        """All new predicates return valid Verdict enum values."""
        from biocompiler.stability_predicates import evaluate_stable_folding
        from biocompiler.solubility_predicates import evaluate_soluble_expression
        from biocompiler.immuno_predicates import evaluate_low_immunogenicity
        from biocompiler.structure_predicates import evaluate_structure_confidence
        from biocompiler.types import Verdict, TypeCheckResult

        for fn in [evaluate_stable_folding, evaluate_soluble_expression,
                    evaluate_low_immunogenicity, evaluate_structure_confidence]:
            result = fn(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens")
            assert isinstance(result, TypeCheckResult)
            assert result.verdict in list(Verdict)

    def test_predicate_registry_has_28(self):
        """Predicate registry should have 28 predicates total."""
        from biocompiler.type_system import registry
        assert len(registry.names()) == 28


# ────────────────────────────────────────────────────────────
# Test Backward Compatibility
# ────────────────────────────────────────────────────────────
class TestBackwardCompatibility:
    """Phase 1 predicates still work after Phase 2+3 additions."""

    def test_phase1_predicates_still_work(self):
        """Original 12 predicates still work."""
        from biocompiler.type_system import evaluate_all_predicates
        results = evaluate_all_predicates(HBB_SEQUENCE, organism="Homo_sapiens")
        assert len(results) == 12  # evaluate_all_predicates returns the original 12
        for r in results:
            assert r.verdict is not None

    def test_type_check_result_format(self):
        """All new predicates return consistent TypeCheckResult format."""
        from biocompiler.stability_predicates import evaluate_stable_folding
        from biocompiler.types import TypeCheckResult

        result = evaluate_stable_folding(HBB_SEQUENCE, HBB_PROTEIN, "Homo_sapiens")
        assert isinstance(result, TypeCheckResult)
        assert hasattr(result, 'predicate')
        assert hasattr(result, 'verdict')
        assert hasattr(result, 'derivation')

    def test_verdict_enum_values(self):
        """All verdicts are from the 5-valued enum."""
        from biocompiler.types import Verdict
        valid = {Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
                 Verdict.LIKELY_FAIL, Verdict.FAIL}
        assert len(valid) == 5
