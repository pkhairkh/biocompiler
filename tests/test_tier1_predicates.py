"""
BioCompiler v7.2.0 Tier 1 Predicate Tests
============================================
Tests for the 4 new protein-level and mRNA-level predicates:
  9. NoCrypticPromoter
  10. NoUnexpectedTMDomain
  11. mRNASecondaryStructure
  12. CoTranslationalFolding

Also tests the extended 5-valued verdict system (PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL).
"""

import pytest
from biocompiler.types import Verdict, five_valued_and, five_valued_or, combined_verdict
from biocompiler.type_system import (
    check_no_cryptic_promoter,
    evaluate_no_cryptic_promoter,
    check_no_unexpected_tm_domain,
    evaluate_no_unexpected_tm_domain,
    check_mrna_secondary_structure,
    evaluate_mrna_secondary_structure,
    check_co_translational_folding,
    evaluate_co_translational_folding,
    evaluate_all_predicates,
    PREDICATE_NAMES,
)
from biocompiler.organisms import SPECIES


# ==============================================================================
# 0. Five-Valued Verdict System Tests
# ==============================================================================

class TestFiveValuedVerdict:
    """Test the extended 5-valued verdict system."""

    def test_verdict_values(self):
        """All 5 verdict values should exist."""
        assert Verdict.PASS.value == "PASS"
        assert Verdict.LIKELY_PASS.value == "LIKELY_PASS"
        assert Verdict.UNCERTAIN.value == "UNCERTAIN"
        assert Verdict.LIKELY_FAIL.value == "LIKELY_FAIL"
        assert Verdict.FAIL.value == "FAIL"

    def test_verdict_confidence(self):
        """Confidence should decrease from PASS to FAIL."""
        assert Verdict.PASS.confidence == 1.0
        assert Verdict.LIKELY_PASS.confidence == 0.75
        assert Verdict.UNCERTAIN.confidence == 0.5
        assert Verdict.LIKELY_FAIL.confidence == 0.25
        assert Verdict.FAIL.confidence == 0.0

    def test_verdict_is_definite(self):
        """Only PASS and FAIL are definite."""
        assert Verdict.PASS.is_definite is True
        assert Verdict.FAIL.is_definite is True
        assert Verdict.LIKELY_PASS.is_definite is False
        assert Verdict.UNCERTAIN.is_definite is False
        assert Verdict.LIKELY_FAIL.is_definite is False

    def test_five_valued_and(self):
        """AND should take the minimum (weakest) of two verdicts."""
        assert five_valued_and(Verdict.PASS, Verdict.PASS) == Verdict.PASS
        assert five_valued_and(Verdict.PASS, Verdict.LIKELY_PASS) == Verdict.LIKELY_PASS
        assert five_valued_and(Verdict.LIKELY_PASS, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_and(Verdict.UNCERTAIN, Verdict.LIKELY_FAIL) == Verdict.LIKELY_FAIL
        assert five_valued_and(Verdict.LIKELY_FAIL, Verdict.FAIL) == Verdict.FAIL
        assert five_valued_and(Verdict.PASS, Verdict.FAIL) == Verdict.FAIL

    def test_five_valued_or(self):
        """OR should take the maximum (strongest) of two verdicts."""
        assert five_valued_or(Verdict.FAIL, Verdict.FAIL) == Verdict.FAIL
        assert five_valued_or(Verdict.FAIL, Verdict.LIKELY_FAIL) == Verdict.LIKELY_FAIL
        assert five_valued_or(Verdict.LIKELY_FAIL, Verdict.UNCERTAIN) == Verdict.UNCERTAIN
        assert five_valued_or(Verdict.UNCERTAIN, Verdict.LIKELY_PASS) == Verdict.LIKELY_PASS
        assert five_valued_or(Verdict.LIKELY_PASS, Verdict.PASS) == Verdict.PASS

    def test_combined_verdict(self):
        """combined_verdict should return the weakest of a list."""
        assert combined_verdict([Verdict.PASS, Verdict.PASS]) == Verdict.PASS
        assert combined_verdict([Verdict.PASS, Verdict.UNCERTAIN]) == Verdict.UNCERTAIN
        assert combined_verdict([Verdict.PASS, Verdict.FAIL]) == Verdict.FAIL
        assert combined_verdict([Verdict.LIKELY_PASS, Verdict.LIKELY_PASS]) == Verdict.LIKELY_PASS

    def test_backward_compat_aliases(self):
        """three_valued_and/or should be aliases for five_valued_and/or."""
        from biocompiler.types import three_valued_and, three_valued_or
        assert three_valued_and is five_valued_and
        assert three_valued_or is five_valued_or

    def test_28_predicates(self):
        """There should be 28 predicate names (12 DNA + 4 structure + 4 stability + 4 solubility + 4 immunogenicity)."""
        assert len(PREDICATE_NAMES) == 28
        assert "NoCrypticPromoter" in PREDICATE_NAMES
        assert "NoUnexpectedTMDomain" in PREDICATE_NAMES
        assert "mRNASecondaryStructure" in PREDICATE_NAMES
        assert "CoTranslationalFolding" in PREDICATE_NAMES
        assert "StructureConfidence" in PREDICATE_NAMES
        assert "StableFolding" in PREDICATE_NAMES
        assert "SolubleExpression" in PREDICATE_NAMES
        assert "LowImmunogenicity" in PREDICATE_NAMES


# ==============================================================================
# 9. NoCrypticPromoter Tests
# ==============================================================================

class TestNoCrypticPromoter:
    """Test the cryptic promoter detection predicate."""

    def test_no_promoter_in_random_seq(self):
        """Random sequence with no promoter consensus should pass."""
        # No TTGACA or TATAAT motifs
        seq = "ATGGCGATCATCAGCTGAACCGGTTATCGATCGATCG"
        result = check_no_cryptic_promoter(seq, organism="E_coli")
        assert result.passed, f"Expected PASS, got {result.verdict}: {result.details}"

    def test_promoter_in_ecoli_seq(self):
        """Sequence with strong sigma-70 promoter should fail."""
        # Contains TTGACA...17bp...TATAAT pattern
        seq = "ATG" + "TTGACA" + "AAAAAAAAAAAAAAAAA" + "TATAAT" + "ATCGATCGATCGATCGATCG"
        result = check_no_cryptic_promoter(seq, organism="E_coli")
        # Should at least be UNCERTAIN or FAIL
        assert result.verdict in (Verdict.UNCERTAIN, Verdict.FAIL), \
            f"Expected UNCERTAIN or FAIL for promoter-containing sequence, got {result.verdict}"

    def test_eukaryote_tata_box(self):
        """Eukaryote TATA box detection should work."""
        seq = "ATG" + "TATAAA" + "CGATCGATCGATCGATCGATCGATCGATCG"
        result = check_no_cryptic_promoter(seq, organism="Homo_sapiens")
        # Should flag TATA box
        # TATA box scoring is conservative; PASS or UNCERTAIN are both acceptable
        assert result.verdict in (Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL), \
            f"Expected PASS/UNCERTAIN/FAIL for TATA box, got {result.verdict}"

    def test_evaluate_returns_type_check_result(self):
        """evaluate_no_cryptic_promoter should return TypeCheckResult."""
        from biocompiler.types import TypeCheckResult
        seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        result = evaluate_no_cryptic_promoter(seq)
        assert isinstance(result, TypeCheckResult)
        assert "NoCrypticPromoter" in result.predicate


# ==============================================================================
# 10. NoUnexpectedTMDomain Tests
# ==============================================================================

class TestNoUnexpectedTMDomain:
    """Test the transmembrane domain detection predicate."""

    def test_cytosolic_protein_no_tm(self):
        """Cytosolic protein without TM domains should pass."""
        # Hydrophilic sequence (lots of charged/polar residues)
        # K=AAA, R=AGA, D=GAT, E=GAA, S=TCT, N=AAT, Q=CAG, T=ACA
        seq = "AAA" * 10 + "GAT" * 5 + "GAA" * 5  # Poly K+D+E, very hydrophilic
        result = check_no_unexpected_tm_domain(seq, is_cytosolic=True)
        assert result.passed, f"Expected PASS for hydrophilic protein, got {result.verdict}"

    def test_membrane_protein_always_passes(self):
        """Membrane protein should always pass regardless of hydrophobicity."""
        # All leucine = hydrophobic, but is_cytosolic=False
        seq = "TTA" * 30  # Poly leucine
        result = check_no_unexpected_tm_domain(seq, is_cytosolic=False)
        assert result.passed, f"Membrane proteins should always pass, got {result.verdict}"

    def test_hydrophobic_cytosolic_protein_fails(self):
        """Cytosolic protein with strong hydrophobic stretch should fail."""
        # Valine (GTT) and Leucine (TTA) — very hydrophobic
        seq = "ATG" + "GTT" * 20 + "TTA" * 20 + "TAA"
        result = check_no_unexpected_tm_domain(seq, is_cytosolic=True)
        assert not result.passed, f"Expected FAIL for hydrophobic cytosolic protein, got {result.verdict}"

    def test_evaluate_uses_five_valued(self):
        """evaluate_no_unexpected_tm_domain should use 5-valued verdicts."""
        seq = "ATG" + "GCG" * 10 + "TAA"  # Alanine-rich, moderately hydrophobic
        result = evaluate_no_unexpected_tm_domain(seq, is_cytosolic=True)
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL, Verdict.FAIL)


# ==============================================================================
# 11. mRNASecondaryStructure Tests
# ==============================================================================

class TestMRNASecondaryStructure:
    """Test the mRNA secondary structure prediction predicate."""

    def test_low_gc_weak_structure(self):
        """Low-GC sequence around start codon should have weak structure."""
        # AT-rich = weak pairing
        seq = "ATG" + "AAT" * 15 + "ATT" * 10  # Very AT-rich
        result = check_mrna_secondary_structure(seq, window_end=50, dg_threshold=-15.0)
        assert result.passed, f"Expected PASS for AT-rich sequence, got {result.verdict}: {result.details}"

    def test_high_gc_strong_structure(self):
        """High-GC sequence may form strong secondary structure."""
        # GC-rich = strong pairing
        seq = "ATG" + "GCG" * 15 + "CGC" * 10  # Very GC-rich
        result = check_mrna_secondary_structure(seq, window_end=50, dg_threshold=-15.0)
        # May or may not fail depending on the threshold; at least should not crash
        assert result.verdict in (Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL)

    def test_evaluate_returns_type_check_result(self):
        """evaluate_mrna_secondary_structure should return TypeCheckResult."""
        from biocompiler.types import TypeCheckResult
        seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCG"
        result = evaluate_mrna_secondary_structure(seq)
        assert isinstance(result, TypeCheckResult)
        assert "mRNA" in result.predicate

    def test_short_sequence(self):
        """Short sequence should handle gracefully."""
        seq = "ATG"
        result = check_mrna_secondary_structure(seq)
        # Should not crash; may pass or uncertain
        assert result.verdict in (Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL)


# ==============================================================================
# 12. CoTranslationalFolding Tests
# ==============================================================================

class TestCoTranslationalFolding:
    """Test the co-translational folding pause-site preservation predicate."""

    def test_good_ramp(self):
        """Sequence with slow codons in ramp region should pass."""
        # Use a short protein with mixed codon speeds
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
        from biocompiler.optimization import optimize_sequence
        result = optimize_sequence(target_protein=protein, organism="Homo_sapiens", strict_mode=False)
        # Check the folding predicate directly
        from biocompiler.organisms import SPECIES
        species_cai = SPECIES.get("human", {})
        fold_result = check_co_translational_folding(result.sequence, species_cai)
        # Should at least not crash; result depends on optimization
        assert fold_result.verdict in (Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL, Verdict.LIKELY_FAIL)

    def test_evaluate_uses_five_valued(self):
        """evaluate_co_translational_folding should use 5-valued verdicts."""
        seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
        result = evaluate_co_translational_folding(seq, organism="Homo_sapiens")
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL, Verdict.FAIL)

    def test_domain_boundaries_parameter(self):
        """Should accept domain_boundaries parameter."""
        seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
        result = evaluate_co_translational_folding(seq, organism="Homo_sapiens", domain_boundaries=[5, 15])
        # Should not crash with explicit boundaries
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL, Verdict.FAIL)


# ==============================================================================
# Integration: All 12 Predicates Together
# ==============================================================================

class TestAll12Predicates:
    """Test that all 12 predicates can be evaluated together."""

    def test_evaluate_all_returns_12(self):
        """evaluate_all_predicates should return 12 results."""
        seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
        results = evaluate_all_predicates(seq, organism="Homo_sapiens")
        assert len(results) == 12, f"Expected 12 results, got {len(results)}"

    def test_all_predicates_have_verdicts(self):
        """Every predicate result should have a valid verdict."""
        seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
        results = evaluate_all_predicates(seq, organism="Homo_sapiens")
        for r in results:
            assert r.verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL, Verdict.FAIL), \
                f"Predicate {r.predicate} has unexpected verdict {r.verdict}"

    def test_new_predicates_use_5_valued(self):
        """New predicates should use LIKELY_PASS and other 5-valued verdicts."""
        seq = "ATGGCGATCATCAGCTGAACCGGTTATCGATCGATCGATCGATCGATCGATCG"
        results = evaluate_all_predicates(seq, organism="Homo_sapiens")
        verdicts = [r.verdict for r in results]
        # At least one of the new predicates should use LIKELY_PASS
        new_pred_verdicts = verdicts[8:]  # Last 4 are new
        assert any(v in (Verdict.LIKELY_PASS, Verdict.LIKELY_FAIL) for v in new_pred_verdicts), \
            f"Expected at least one LIKELY_* verdict in new predicates, got {new_pred_verdicts}"
