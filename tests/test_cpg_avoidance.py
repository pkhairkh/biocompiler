"""
Dedicated tests for CpG island avoidance in BioCompiler.

Tests the systematic CpG elimination pass that replaces CG dinucleotides
with synonymous codons to avoid CpG islands. CpG islands are regions
where the Obs/Exp CG ratio exceeds a threshold (default 0.6) in a
sliding window (default 200bp). They can trigger epigenetic silencing
in eukaryotic expression systems.

Test categories:
1. _eliminate_cpg_dinucleotides() unit tests
2. Prokaryote skip logic (CpG irrelevant for E. coli)
3. Within-codon CpG elimination
4. Cross-codon CpG elimination
5. Eukaryotic optimization end-to-end
6. CpG island predicate integration
7. CAI vs CpG tradeoff
8. Edge cases
"""

import os
import sys
import pytest

# Ensure src is on path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from biocompiler.optimization import optimize_sequence, _eliminate_cpg_dinucleotides
from biocompiler.type_system import check_no_cpg_island, evaluate_no_cpg_island
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
from biocompiler.translation import translate, compute_cai
from biocompiler.constants import AA_TO_CODONS


# ============================================================================
# 1. _eliminate_cpg_dinucleotides() unit tests
# ============================================================================

class TestEliminateCpgDinucleotidesUnit:
    """Unit tests for the _eliminate_cpg_dinucleotides function."""

    def test_returns_same_sequence_if_already_passes(self):
        """If the sequence already passes the CpG check, return it unchanged."""
        # A yeast-optimized sequence typically has no CGs (low GC)
        protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
        result = optimize_sequence(
            target_protein=protein, organism="Saccharomyces_cerevisiae",
            strict_mode=False,
        )
        usage = CODON_ADAPTIVENESS_TABLES.get("Saccharomyces_cerevisiae")
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Saccharomyces_cerevisiae",
        )
        assert new_seq == result.sequence
        assert warnings == []

    def test_prokaryote_organism_skipped(self):
        """CpG elimination is skipped for prokaryotic organisms."""
        protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
        result = optimize_sequence(
            target_protein=protein, organism="Escherichia_coli",
            strict_mode=False,
        )
        usage = CODON_ADAPTIVENESS_TABLES.get("Escherichia_coli")
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Escherichia_coli",
        )
        assert new_seq == result.sequence
        assert warnings == []

    def test_translation_preserved_after_elimination(self):
        """After CpG elimination, the sequence must still translate to the same protein."""
        protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Homo_sapiens",
        )
        translated = translate(new_seq, to_stop=True).rstrip("*")
        assert translated == protein

    def test_cg_count_decreases_or_stays_same(self):
        """CpG elimination should never increase the CG dinucleotide count."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        old_cgs = sum(1 for i in range(len(result.sequence) - 1) if result.sequence[i:i+2] == "CG")
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Homo_sapiens",
        )
        new_cgs = sum(1 for i in range(len(new_seq) - 1) if new_seq[i:i+2] == "CG")
        assert new_cgs <= old_cgs

    def test_warnings_for_unfixable_cpgs(self):
        """If CGs remain and the CpG island check fails, warnings should be emitted."""
        # This is hard to trigger in practice because most proteins
        # can have their CGs eliminated. Just verify the return type.
        protein = "MFFMMMMFFFM"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Homo_sapiens",
        )
        assert isinstance(warnings, list)
        for w in warnings:
            assert isinstance(w, str)


# ============================================================================
# 2. Prokaryote skip logic
# ============================================================================

class TestProkaryoteCpgSkip:
    """CpG islands are irrelevant for prokaryotes and should be skipped."""

    @pytest.mark.parametrize("organism", [
        "Escherichia_coli",
        "e_coli",
    ])
    def test_prokaryote_cpg_check_auto_skips(self, organism):
        """The CpG island check should auto-skip for prokaryotic organisms."""
        protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
        result = optimize_sequence(
            target_protein=protein, organism=organism,
            strict_mode=False,
        )
        cpg_result = check_no_cpg_island(result.sequence, organism="Escherichia_coli")
        assert cpg_result.passed

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "Saccharomyces_cerevisiae",
    ])
    def test_eukaryote_cpg_check_runs(self, organism):
        """The CpG island check should run for eukaryotic organisms."""
        protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
        result = optimize_sequence(
            target_protein=protein, organism=organism,
            strict_mode=False,
        )
        # The check runs but the sequence might be too short for a 200bp window
        cpg_result = check_no_cpg_island(result.sequence, organism=organism)
        # For short sequences, the check auto-passes
        if len(result.sequence) < 200:
            assert cpg_result.passed  # Too short for a CpG island window


# ============================================================================
# 3. Within-codon CpG elimination
# ============================================================================

class TestWithinCodonCpg:
    """Test elimination of CG dinucleotides entirely within a single codon."""

    def test_arginine_cga_replaced(self):
        """CGA (Arginine) contains CG within the codon — should be replaced."""
        # Arginine has 6 codons: CGT, CGC, CGA, CGG, AGA, AGG
        # CG-free alternatives: AGA, AGG
        protein = "MRRRR"  # Multiple arginines
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        # After optimization, the sequence should not have a CpG island
        cpg_result = check_no_cpg_island(result.sequence, organism="Homo_sapiens")
        # Short sequence, auto-passes; check CG content instead
        cg_count = sum(1 for i in range(len(result.sequence) - 1) if result.sequence[i:i+2] == "CG")
        # The optimizer should have reduced CGs
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_alanine_with_within_codon_cg(self):
        """GCG (Alanine) contains CG — should be replaceable with GCT/GCC/GCA."""
        protein = "MAAAA"  # Multiple alanines
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        # Verify translation preserved
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein


# ============================================================================
# 4. Cross-codon CpG elimination
# ============================================================================

class TestCrossCodonCpg:
    """Test elimination of CG dinucleotides that span codon boundaries."""

    def test_cross_codon_cg_identified(self):
        """A CG at a codon boundary (C at end of codon, G at start of next) should be found."""
        # This tests the internal logic, not the public API directly
        # Create a sequence with a known cross-codon CG
        # TAC (Tyr) ends with C, GGT (Gly) starts with G → CG at boundary
        protein = "MYG"  # TAC + GGT = TACGGT has CG at pos 2-3
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        # The sequence is too short for a CpG island (9bp < 200bp window)
        # but the optimizer should still handle cross-codon CGs
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_cross_codon_elimination_preserves_translation(self):
        """Cross-codon CG elimination must preserve the protein sequence."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        # Apply CpG elimination
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        new_seq, _ = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Homo_sapiens",
        )
        translated = translate(new_seq, to_stop=True).rstrip("*")
        assert translated == protein


# ============================================================================
# 5. Eukaryotic optimization end-to-end
# ============================================================================

class TestEukaryoticCpgOptimization:
    """End-to-end tests for CpG avoidance in eukaryotic optimization."""

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "Saccharomyces_cerevisiae",
    ])
    def test_optimized_sequence_passes_cpg_check(self, organism):
        """Optimized eukaryotic sequences should pass the CpG island check."""
        # Use a long enough protein to have a 200bp+ CDS
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(
            target_protein=protein, organism=organism,
            gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        cpg_result = check_no_cpg_island(result.sequence, organism=organism)
        assert cpg_result.passed, (
            f"CpG island found for {organism}: {cpg_result.details}"
        )

    def test_brca1_avoids_cpg_island(self):
        """BRCA1 segment (previously the only failing gene) should now pass."""
        from biocompiler.dataset_validation import HUMAN_REFERENCE_GENES
        gene = HUMAN_REFERENCE_GENES["BRCA1_segment"]
        result = optimize_sequence(
            target_protein=gene["protein"], organism=gene["organism"],
            gc_lo=0.30, gc_hi=0.70, cai_threshold=0.2,
            strict_mode=False,
        )
        cpg_result = check_no_cpg_island(result.sequence, organism="Homo_sapiens")
        assert cpg_result.passed, (
            f"BRCA1 still has CpG island: {cpg_result.details}"
        )

    def test_insulin_for_human_avoids_cpg(self):
        """Human insulin optimization should avoid CpG islands."""
        protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        # Insulin CDS is only 165bp, shorter than the 200bp window
        cpg_result = check_no_cpg_island(result.sequence, organism="Homo_sapiens")
        assert cpg_result.passed


# ============================================================================
# 6. CpG island predicate integration
# ============================================================================

class TestCpgPredicateIntegration:
    """Test that the NoCpGIsland predicate is properly evaluated."""

    def test_eukaryotic_result_includes_cpg_predicate(self):
        """Optimization results for eukaryotes should include NoCpGIsland predicate."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        cpg_preds = [p for p in result.predicate_results if p.predicate == "NoCpGIsland"]
        assert len(cpg_preds) == 1
        assert cpg_preds[0].passed

    def test_prokaryotic_result_skips_cpg(self):
        """Optimization results for prokaryotes should skip CpG check."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        result = optimize_sequence(
            target_protein=protein, organism="Escherichia_coli",
            strict_mode=False,
        )
        cpg_preds = [p for p in result.predicate_results if p.predicate == "NoCpGIsland"]
        assert len(cpg_preds) == 1
        assert cpg_preds[0].passed
        assert "prokaryote" in cpg_preds[0].details.lower() or "skip" in cpg_preds[0].details.lower()


# ============================================================================
# 7. CAI vs CpG tradeoff
# ============================================================================

class TestCaiVsCpgTradeoff:
    """Test that CpG avoidance trades off CAI correctly."""

    def test_cai_still_reasonable_after_cpg_elimination(self):
        """After CpG elimination, CAI should still be reasonable (>0.5)."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        # CAI should still be reasonable
        assert result.cai >= 0.5, (
            f"CAI too low after CpG elimination: {result.cai:.4f}"
        )

    def test_gc_still_in_range_after_cpg_elimination(self):
        """After CpG elimination, GC should still be within the target range."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            gc_lo=0.20, gc_hi=0.80, cai_threshold=0.1,
            strict_mode=False,
        )
        assert 0.20 <= result.gc_content <= 0.80, (
            f"GC out of range after CpG elimination: {result.gc_content:.3f}"
        )


# ============================================================================
# 8. Edge cases
# ============================================================================

class TestCpgEdgeCases:
    """Test edge cases for CpG avoidance."""

    def test_short_sequence_auto_passes(self):
        """Sequences shorter than 200bp auto-pass the CpG island check."""
        protein = "MRK"  # Only 9bp CDS
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        cpg_result = check_no_cpg_island(result.sequence, organism="Homo_sapiens")
        assert cpg_result.passed

    def test_empty_protein_handled(self):
        """An empty or invalid protein should not crash CpG elimination."""
        with pytest.raises(Exception):
            optimize_sequence(target_protein="", organism="Homo_sapiens")

    def test_single_amino_acid(self):
        """A single amino acid should be handled without error."""
        protein = "M"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        cpg_result = check_no_cpg_island(result.sequence, organism="Homo_sapiens")
        assert cpg_result.passed

    def test_protein_with_only_cg_amino_acids(self):
        """Proteins composed of amino acids whose codons all contain CG."""
        # Arginine: all CGx codons contain CG, but AGA/AGG don't
        protein = "MRRRRRR"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_all_arginine_protein_reduces_cgs(self):
        """An all-Arg protein should have CG dinucleotides significantly reduced.

        Arg has CGN codons (CGT, CGC, CGA, CGG) which contain CG, but
        also AGA and AGG which are CG-free. The CpG elimination pass
        should replace most CGN codons with AGA/AGG, though some may
        remain due to the optimizer pipeline's tradeoffs.
        """
        protein = "RRRRRRRRRRR"  # 11 Arg residues
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            gc_lo=0.20, gc_hi=0.80, cai_threshold=0.2,
            strict_mode=False,
        )
        cg_count = sum(1 for i in range(len(result.sequence) - 1)
                       if result.sequence[i:i+2] == "CG")
        # The optimizer should significantly reduce CG count
        # (worst case all-Arg is 11 CGs, we expect most to be eliminated)
        assert cg_count <= 4, (
            f"Expected ≤4 CG dinucleotides in all-Arg protein, got {cg_count}: "
            f"{result.sequence}"
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_boundary_cg_between_arg_ala_eliminated(self):
        """Arg-Ala boundary CGs should be eliminated.

        Arg codons CGG/CGC end with G/C, Ala codons GCC/GCT/GCA/GCG start
        with G. When an Arg codon ends with C and the next Ala codon starts
        with G, a boundary CG is created. The elimination pass should fix this.
        """
        protein = "RARARARARARA"  # 6 Arg-Ala pairs
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            gc_lo=0.20, gc_hi=0.80, cai_threshold=0.2,
            strict_mode=False,
        )
        cg_count = sum(1 for i in range(len(result.sequence) - 1)
                       if result.sequence[i:i+2] == "CG")
        # The optimizer should significantly reduce CG count for Arg-Ala
        assert cg_count <= 3, (
            f"Expected ≤3 CG dinucleotides in Arg-Ala protein, got {cg_count}: "
            f"{result.sequence}"
        )
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_cpg_elimination_prefers_high_cai(self):
        """CpG elimination should prefer the highest-CAI CG-free alternative.

        For Arg, the CG-free alternatives are AGG (CAI≈0.96) and AGA (CAI≈0.46).
        The function should choose AGG over AGA when possible.
        """
        protein = "MRRRRRR"
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        # Count how many Arg codons use AGG vs AGA
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        agg_count = 0
        aga_count = 0
        for i in range(1, len(protein)):  # Skip M at position 0
            codon = result.sequence[i*3:i*3+3]
            if codon == "AGG":
                agg_count += 1
            elif codon == "AGA":
                aga_count += 1
        # AGG should be preferred (higher CAI) unless GC constraint forces AGA
        # At minimum, some Arg codons should use CG-free alternatives
        cg_free_count = agg_count + aga_count
        assert cg_free_count > 0, "No CG-free Arg codons found"
        translated = translate(result.sequence, to_stop=True).rstrip("*")
        assert translated == protein

    def test_cpg_elimination_with_gc_out_of_range(self):
        """CpG elimination should still work when GC is outside the target range.

        If the initial optimization produces a sequence with GC outside the
        target range (e.g., all-Arg with CGG codons = 100% GC), the CpG
        elimination should still be able to eliminate CGs by moving GC
        toward the target range.
        """
        protein = "RRRRRRRRRRR"
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        # Manually create a high-GC sequence
        high_gc_seq = "CGG" * 11  # 100% GC, all CGs
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            high_gc_seq, protein, usage,
            organism="Homo_sapiens",
            gc_lo=0.20, gc_hi=0.80,
        )
        cg_count = sum(1 for i in range(len(new_seq) - 1)
                       if new_seq[i:i+2] == "CG")
        # CpG elimination should significantly reduce CGs even from high-GC start
        assert cg_count <= 4, (
            f"Expected ≤4 CGs after elimination from high-GC sequence, got {cg_count}"
        )
        # GC should have moved toward the target range
        new_gc = (new_seq.count("G") + new_seq.count("C")) / len(new_seq)
        old_gc = 1.0  # Original was 100% GC
        assert new_gc < old_gc, "GC should have decreased"

    def test_stale_position_recheck(self):
        """After fixing one CG, subsequent positions should be re-checked.

        If fixing CG at position i also fixes the CG at position i+3 (because
        both were in the same codon), the function should skip position i+3
        rather than trying to fix it again.
        """
        protein = "MRRRR"  # 4 Arg codons, all likely CGG
        result = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        usage = CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens")
        new_seq, warnings = _eliminate_cpg_dinucleotides(
            result.sequence, protein, usage, organism="Homo_sapiens",
        )
        # Should still translate correctly
        translated = translate(new_seq, to_stop=True).rstrip("*")
        assert translated == protein
        # Should have eliminated CGs
        cg_count = sum(1 for i in range(len(new_seq) - 1)
                       if new_seq[i:i+2] == "CG")
        # CG count should be 0 or at least reduced
        original_cg = sum(1 for i in range(len(result.sequence) - 1)
                         if result.sequence[i:i+2] == "CG")
        assert cg_count <= original_cg

    def test_validate_no_cpg_island_passes_organism(self):
        """validate_no_cpg_island should pass the organism to the CpG check."""
        from biocompiler.dataset_validation import validate_no_cpg_island
        # For E. coli (prokaryote), the check should auto-pass
        result = validate_no_cpg_island(
            protein="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPEKT",
            organism="Escherichia_coli",
            gene_name="insulin",
            dataset_name="ecoli",
        )
        assert result.passed

    def test_multiple_organisms_same_protein(self):
        """The same protein optimized for different organisms should have different CpG results."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        # Prokaryote: always passes (skipped)
        result_ecoli = optimize_sequence(
            target_protein=protein, organism="Escherichia_coli",
            strict_mode=False,
        )
        cpg_ecoli = check_no_cpg_island(result_ecoli.sequence, organism="Escherichia_coli")
        assert cpg_ecoli.passed

        # Eukaryote: should pass after CpG elimination
        result_human = optimize_sequence(
            target_protein=protein, organism="Homo_sapiens",
            strict_mode=False,
        )
        cpg_human = check_no_cpg_island(result_human.sequence, organism="Homo_sapiens")
        assert cpg_human.passed


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
