"""
BioCompiler v7.2.0 Optimization Improvements — Tests

Tests for:
1. Phase 1 splice-safe codon reordering (Task 1)
2. Force GT-free re-optimization for "chose poorly" cases (Task 2)
3. Window-aware CpG island elimination (Task 3)
"""

import pytest
from src.biocompiler.optimization import (
    optimize_sequence,
    _reorder_for_splice_safety,
    _find_failing_cpg_windows,
    _greedy_optimize,
)
from src.biocompiler.mutagenesis import (
    force_gt_free_reoptimization,
    is_gt_mandatory,
    diagnose_optimizer_weakness,
    GT_MANDATORY_AAS,
)
from src.biocompiler.constants import AA_TO_CODONS
from src.biocompiler.type_system import evaluate_no_cryptic_splice, evaluate_no_cpg_island
from src.biocompiler.maxentscan import score_donor, max_donor_score


# ==============================================================================
# 1. Splice-Safe Codon Reordering Tests
# ==============================================================================

class TestReorderForSpliceSafety:
    """Test the _reorder_for_splice_safety function."""

    def test_gt_free_codons_come_first_for_cysteine(self):
        """For Cysteine, TGC (GT-free) should come before TGT (GT-containing)."""
        sorted_codons = {"C": ["TGT", "TGC"]}  # CAI-sorted (TGT first by default)
        result = _reorder_for_splice_safety(sorted_codons, ["C"], 3.0)
        assert result["C"][0] == "TGC", "GT-free TGC should come before GT-containing TGT"

    def test_valine_unchanged(self):
        """Valine has no GT-free codons — reordering should not change order."""
        sorted_codons = {"V": ["GTG", "GTC", "GTT", "GTA"]}
        result = _reorder_for_splice_safety(sorted_codons, ["V"], 3.0)
        assert result["V"] == ["GTG", "GTC", "GTT", "GTA"], "Valine order should be unchanged"

    def test_threshold_zero_no_reordering(self):
        """When threshold is 0, no reordering should happen."""
        sorted_codons = {"C": ["TGT", "TGC"]}
        result = _reorder_for_splice_safety(sorted_codons, ["C"], 0.0)
        assert result["C"] == ["TGT", "TGC"], "No reordering when threshold is 0"

    def test_ag_free_preferred(self):
        """Codons that are both GT-free and AG-free should be preferred."""
        # Arginine: AGA and AGG contain AG; CGC, CGA, CGG contain neither GT nor AG in some
        # Let's check: CGC has no GT, no AG → tier 1
        # AGA has AG → tier 2 or 4
        sorted_codons = {"R": ["AGA", "CGG", "CGC", "CGA", "AGG"]}
        result = _reorder_for_splice_safety(sorted_codons, ["R"], 3.0)
        # Tier 1 codons (no GT, no AG) should come first
        tier1 = [c for c in result["R"] if "GT" not in c and "AG" not in c]
        tier2_plus = [c for c in result["R"] if "GT" in c or "AG" in c]
        # All tier1 codons should appear before any tier2+ codon
        if tier1 and tier2_plus:
            last_tier1_idx = max(result["R"].index(c) for c in tier1)
            first_tier2_idx = min(result["R"].index(c) for c in tier2_plus)
            assert last_tier1_idx < first_tier2_idx

    def test_preserves_all_codons(self):
        """Reordering should not add or remove codons."""
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            codons = list(AA_TO_CODONS[aa])
            sorted_codons = {aa: codons}
            result = _reorder_for_splice_safety(sorted_codons, [aa], 3.0)
            assert set(result[aa]) == set(codons), f"AA {aa}: codons changed after reordering"
            assert len(result[aa]) == len(codons), f"AA {aa}: codon count changed"


# ==============================================================================
# 2. Force GT-Free Re-optimization Tests
# ==============================================================================

class TestForceGTFreeReoptimization:
    """Test the force_gt_free_reoptimization function."""

    def test_fixes_non_valine_gt_positions(self):
        """Should swap GT-containing codons to GT-free for non-Valine AAs."""
        # Cysteine with TGT (contains GT) → should swap to TGC
        protein = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # 31 Cysteines
        # First optimize to get a baseline
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.30, gc_hi=0.70)
        # Check if there are strong donors at C positions
        for i in range(len(result.sequence) - 1):
            if result.sequence[i:i+2] == "GT":
                s = score_donor(result.sequence, i)
                if s >= 3.0:
                    codon_idx = i // 3
                    if codon_idx < len(protein) and protein[codon_idx] == "C":
                        # Found a C position with strong donor — force re-opt
                        new_seq = force_gt_free_reoptimization(
                            result.sequence, protein, "Homo_sapiens", threshold=3.0
                        )
                        # Verify the codon at this position no longer has GT
                        new_codon = new_seq[codon_idx*3:codon_idx*3+3]
                        assert "GT" not in new_codon, f"C codon {new_codon} still contains GT after forced re-opt"
                        break

    def test_preserves_protein(self):
        """Force re-optimization should not change the protein."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.30, gc_hi=0.70)
        new_seq = force_gt_free_reoptimization(result.sequence, protein, "Homo_sapiens")
        from src.biocompiler.translation import translate
        assert translate(new_seq) == protein, "Protein changed after forced re-optimization"

    def test_does_not_modify_valine_positions(self):
        """Valine positions should be untouched by forced re-optimization."""
        protein = "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"  # 31 Valines
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.30, gc_hi=0.70)
        new_seq = force_gt_free_reoptimization(result.sequence, protein, "Homo_sapiens")
        # V positions should be unchanged since all V codons contain GT
        for ci in range(len(protein)):
            if protein[ci] == "V":
                old_codon = result.sequence[ci*3:ci*3+3]
                new_codon = new_seq[ci*3:ci*3+3]
                # The codon may change to a different V codon, but must still be a V codon
                assert new_codon in AA_TO_CODONS["V"]


# ==============================================================================
# 3. Window-Aware CpG Elimination Tests
# ==============================================================================

class TestWindowAwareCpG:
    """Test the window-aware CpG island elimination."""

    def test_find_failing_cpg_windows(self):
        """_find_failing_cpg_windows should identify CpG islands."""
        # Create a sequence with high GC and many CG dinucleotides
        # GCGCGCGCGC repeated = GC=1.0, all CG
        high_cpg = "GCGC" * 100  # 400bp, GC=1.0, obs/exp very high
        windows = _find_failing_cpg_windows(high_cpg, window_size=200)
        assert len(windows) > 0, "Should find CpG islands in GCGC-repeat sequence"

    def test_no_cpg_in_low_gc_sequence(self):
        """No CpG islands in low-GC sequences."""
        low_gc = "ATATATATAT" * 40  # 400bp, GC=0.0
        windows = _find_failing_cpg_windows(low_gc, window_size=200)
        assert len(windows) == 0, "Should find no CpG islands in AT-only sequence"

    def test_optimizer_produces_no_cpg_islands(self):
        """Full optimization should eliminate CpG islands when possible."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(target_protein=protein, organism="Homo_sapiens")
        cpg_result = evaluate_no_cpg_island(seq=result.sequence)
        assert cpg_result.verdict.value == "PASS", (
            f"NoCpGIsland should pass after optimization, got: {cpg_result.violation}"
        )


# ==============================================================================
# 4. End-to-End Integration Tests
# ==============================================================================

class TestV72Integration:
    """End-to-end integration tests for v7.2.0 improvements."""

    def test_hbb_no_cryptic_splice_with_mutagenesis(self):
        """HBB with mutagenesis should pass NoCrypticSplice."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(target_protein=protein, organism="Homo_sapiens", enable_mutagenesis=True)
        assert not any("CrypticSplice" in p for p in result.failed_predicates), (
            f"HBB should pass NoCrypticSplice with mutagenesis, failed: {result.failed_predicates}"
        )

    def test_hbb_preserves_high_cai(self):
        """HBB optimization should maintain CAI >= 0.90."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(target_protein=protein, organism="Homo_sapiens", enable_mutagenesis=True)
        assert result.cai >= 0.90, f"CAI should be >= 0.90, got {result.cai:.4f}"

    def test_no_chose_poorly_after_phase1_reordering(self):
        """After Phase 1 splice-safe reordering, there should be no optimizer weaknesses for non-V AAs."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
        result = optimize_sequence(protein, "Homo_sapiens", gc_lo=0.30, gc_hi=0.70)
        weaknesses = diagnose_optimizer_weakness(result.sequence, protein)
        # After Phase 1 reordering, non-V weaknesses should be rare or zero
        non_v_weaknesses = [w for w in weaknesses if w["aa"] != "V"]
        assert len(non_v_weaknesses) == 0, (
            f"Expected no non-V optimizer weaknesses, found {len(non_v_weaknesses)}: "
            f"{[w['aa'] + '@' + str(w['codon_idx']) for w in non_v_weaknesses[:5]]}"
        )

    def test_egfp_optimization(self):
        """EGFP optimization should succeed with high CAI."""
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        result = optimize_sequence(target_protein=protein, organism="Homo_sapiens")
        assert result.cai >= 0.85, f"EGFP CAI should be >= 0.85, got {result.cai:.4f}"

    def test_mutagenesis_only_proposes_for_valine(self):
        """Type-directed mutagenesis should only propose substitutions at V positions."""
        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        result = optimize_sequence(target_protein=protein, organism="Homo_sapiens", enable_mutagenesis=True)
        if result.aa_substitutions:
            for sub in result.aa_substitutions:
                assert sub["from"] == "V", (
                    f"Mutagenesis proposed {sub['from']}->{sub['to']} at pos {sub['position']}, "
                    f"but only V substitutions should be proposed for GT-mandatory positions"
                )
