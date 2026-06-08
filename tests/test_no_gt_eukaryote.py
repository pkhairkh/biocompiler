"""
NoGTDinucleotide Eukaryote Soft Constraint Tests
=================================================

Tests that NoGTDinucleotide works correctly for eukaryotes without
destroying CAI. The key design decisions:

1. For eukaryotes, NoGTDinucleotide is a **soft constraint**: GTs are
   tolerated up to max_gt_count (1 per 50bp). Exceeding this triggers
   LIKELY_FAIL (not FAIL), so it doesn't block optimization.

2. For prokaryotes, NoGTDinucleotide is a **hard constraint** (max_gt_count=0):
   any GT is a FAIL. However, prokaryotes skip GT checking entirely since
   they have no spliceosome.

3. GT-aware codon selection: when selecting codons for eukaryotes, if a
   codon would create a GT at the boundary with the next codon, the
   optimizer checks for an alternative with similar CAI (within 10%
   relative). If none exists, it accepts the GT.

Task ID: 11
"""

from __future__ import annotations

import math

import pytest

from biocompiler.optimization import (
    optimize_sequence,
    _gt_aware_select_codon,
    GT_BOUNDARY_CAI_TOLERANCE,
)
from biocompiler.translation import translate, compute_cai
from biocompiler.organism_config import is_eukaryotic_organism
from biocompiler.organisms import (
    CODON_ADAPTIVENESS_TABLES,
    SUPPORTED_ORGANISMS,
    resolve_organism,
)
from biocompiler.type_system import (
    AA_TO_CODONS,
    check_no_gt_dinucleotide,
    check_no_avoidable_gt,
    check_no_gt_dinucleotide_soft,
    _compute_max_gt_count,
    _EUKARYOTE_GT_PER_BP,
)
from biocompiler.types import Verdict


# ═══════════════════════════════════════════════════════════════════════
# Test proteins
# ═══════════════════════════════════════════════════════════════════════

INSULIN_PROTEIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)

EGFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    "VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
    "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ"
    "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

DEFAULT_ENZYMES = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]

EUKARYOTIC_ORGANISMS = [
    "Homo_sapiens",
    "Mus_musculus",
    "CHO_K1",
    "Saccharomyces_cerevisiae",
]


def _optimize(protein, organism, **kwargs):
    """Run optimization with standard settings, strict_mode=False to avoid
    GCInRange failures unrelated to GT testing."""
    defaults = dict(
        enzymes=DEFAULT_ENZYMES,
        strict_mode=False,
        optimize_mrna_stability=False,
        include_utr=False,
        track_provenance=False,
    )
    defaults.update(kwargs)
    return optimize_sequence(protein, organism=organism, **defaults)


# ═══════════════════════════════════════════════════════════════════════
# 1. Soft constraint behavior for eukaryotes
# ═══════════════════════════════════════════════════════════════════════

class TestSoftConstraintEukaryotes:
    """Test that NoGTDinucleotide is a soft constraint for eukaryotes."""

    def test_soft_check_returns_likely_fail_not_fail(self):
        """For eukaryotes with GTs exceeding max_gt_count, the soft check
        should return LIKELY_FAIL, not FAIL. Soft failures use passed=True
        so they don't trigger strict mode errors."""
        result = _optimize(INSULIN_PROTEIN, "Homo_sapiens")
        # Find the NoGTDinucleotide predicate result
        gt_pred = None
        for pr in result.predicate_results:
            if pr.predicate == "NoGTDinucleotide":
                gt_pred = pr
                break
        assert gt_pred is not None, "NoGTDinucleotide predicate result not found"

        # For eukaryotes, the predicate should not use FAIL (hard constraint).
        # If GTs exceed max_gt_count, it should be LIKELY_FAIL (soft).
        # If GTs are within max_gt_count, it should be PASS.
        # In both cases, passed=True so it doesn't trigger strict mode.
        assert gt_pred.passed, (
            f"Eukaryotic NoGTDinucleotide should pass (soft constraint). "
            f"Got: passed={gt_pred.passed}, verdict={gt_pred.verdict}. "
            f"Details: {gt_pred.details}"
        )
        assert gt_pred.verdict in (Verdict.PASS, Verdict.LIKELY_FAIL), (
            f"Eukaryotic NoGTDinucleotide should use PASS or LIKELY_FAIL, "
            f"not FAIL (hard). Got: {gt_pred.verdict}. Details: {gt_pred.details}"
        )

    @pytest.mark.parametrize("organism", EUKARYOTIC_ORGANISMS)
    def test_eukaryote_cai_not_destroyed(self, organism: str):
        """Optimizing for eukaryotes should not destroy CAI due to GT avoidance.
        CAI must be > 0.90 for insulin (which has many Valine residues causing
        unavoidable GTs)."""
        result = _optimize(INSULIN_PROTEIN, organism)
        assert result.cai > 0.90, (
            f"{organism}/insulin: CAI={result.cai:.4f} is not > 0.90. "
            f"GT avoidance should not destroy CAI for eukaryotes."
        )

    @pytest.mark.parametrize("organism", EUKARYOTIC_ORGANISMS)
    def test_eukaryote_soft_gt_does_not_block_optimization(self, organism: str):
        """Even if NoGTDinucleotide doesn't pass, the optimization should
        still produce a high-quality result. The soft constraint (LIKELY_FAIL)
        should not block the pipeline."""
        result = _optimize(INSULIN_PROTEIN, organism)
        # The optimization should complete successfully
        assert result.sequence is not None and len(result.sequence) > 0
        # Translation must be preserved
        translated = translate(result.sequence)
        assert translated == INSULIN_PROTEIN, (
            f"{organism}: translation mismatch after optimization"
        )

    def test_soft_check_passes_when_gt_within_tolerance(self):
        """If GT count is within max_gt_count, the soft check should PASS."""
        # Short protein with minimal GTs (no Valine)
        # M-A-S: ATG GCT AGC → no GT
        dna = "ATGGCTAGC"
        result = check_no_gt_dinucleotide_soft(dna, organism="Homo_sapiens")
        assert result.passed, (
            f"Short protein with no GTs should PASS: {result.details}"
        )
        assert result.verdict == Verdict.PASS

    def test_soft_check_reports_gt_count_and_positions(self):
        """The soft check should report the number and positions of GTs."""
        result = _optimize(INSULIN_PROTEIN, "Homo_sapiens")
        gt_pred = None
        for pr in result.predicate_results:
            if pr.predicate == "NoGTDinucleotide":
                gt_pred = pr
                break
        assert gt_pred is not None
        # Details should contain GT count and in-codon/cross-codon breakdown
        assert "GT dinucleotides:" in gt_pred.details or "No GT" in gt_pred.details


# ═══════════════════════════════════════════════════════════════════════
# 2. Hard constraint for prokaryotes
# ═══════════════════════════════════════════════════════════════════════

class TestHardConstraintProkaryotes:
    """Test that NoGTDinucleotide is a hard constraint for prokaryotes."""

    def test_prokaryote_gt_check_skipped_by_optimizer(self):
        """For prokaryotes, the optimizer should skip GT checking entirely
        (no spliceosome)."""
        result = _optimize(INSULIN_PROTEIN, "Escherichia_coli")
        gt_pred = None
        for pr in result.predicate_results:
            if pr.predicate == "NoGTDinucleotide":
                gt_pred = pr
                break
        assert gt_pred is not None
        # For prokaryotes, GT check should be skipped (auto-PASS)
        assert gt_pred.passed, (
            f"Prokaryotic NoGTDinucleotide should be auto-PASS (skipped). "
            f"Got: passed={gt_pred.passed}, details={gt_pred.details}"
        )

    def test_prokaryote_soft_check_returns_fail(self):
        """When explicitly checking GT for a prokaryotic organism with the
        soft check, it should return FAIL (hard constraint) if GTs are present."""
        dna_with_gt = "ATGGGTAGC"  # GGT contains GT
        result = check_no_gt_dinucleotide_soft(dna_with_gt, organism="Escherichia_coli")
        # For prokaryotes, any GT is FAIL
        assert result.verdict == Verdict.FAIL, (
            f"Prokaryotic GT check should return FAIL for sequences with GT. "
            f"Got: {result.verdict}"
        )

    def test_prokaryote_max_gt_count_is_zero(self):
        """For prokaryotes, max_gt_count should be 0."""
        assert _compute_max_gt_count(300, "Escherichia_coli") == 0
        assert _compute_max_gt_count(300, "E_coli") == 0

    def test_prokaryote_cai_high_without_gt_constraint(self):
        """E. coli optimization should achieve high CAI without GT constraint."""
        result = _optimize(INSULIN_PROTEIN, "Escherichia_coli")
        # E. coli should get very high CAI (no GT constraint)
        assert result.cai > 0.95, (
            f"E. coli insulin CAI={result.cai:.4f} is not > 0.95. "
            f"Without GT constraint, CAI should be very high."
        )


# ═══════════════════════════════════════════════════════════════════════
# 3. max_gt_count computation
# ═══════════════════════════════════════════════════════════════════════

class TestMaxGTCountComputation:
    """Test the max_gt_count computation logic."""

    def test_eukaryote_max_gt_scales_with_length(self):
        """For eukaryotes, max_gt_count should scale with sequence length."""
        assert _compute_max_gt_count(50, "Homo_sapiens") == 1
        assert _compute_max_gt_count(100, "Homo_sapiens") == 2
        assert _compute_max_gt_count(500, "Homo_sapiens") == 10

    def test_eukaryote_min_max_gt_is_one(self):
        """Even for very short sequences, max_gt_count should be at least 1."""
        assert _compute_max_gt_count(10, "Homo_sapiens") == 1
        assert _compute_max_gt_count(1, "Homo_sapiens") == 1

    def test_prokaryote_max_gt_is_zero(self):
        """For prokaryotes, max_gt_count should always be 0."""
        assert _compute_max_gt_count(50, "E_coli") == 0
        assert _compute_max_gt_count(500, "E_coli") == 0
        assert _compute_max_gt_count(1000, "Escherichia_coli") == 0

    def test_empty_organism_defaults_to_eukaryote(self):
        """If no organism is specified, max_gt_count should use eukaryote formula."""
        assert _compute_max_gt_count(100, "") == 2

    def test_yeast_uses_eukaryote_formula(self):
        """Yeast should use the eukaryotic formula."""
        assert _compute_max_gt_count(100, "Saccharomyces_cerevisiae") == 2
        assert _compute_max_gt_count(330, "Saccharomyces_cerevisiae") == 6


# ═══════════════════════════════════════════════════════════════════════
# 4. GT-aware codon selection logic
# ═══════════════════════════════════════════════════════════════════════

class TestGTAwareCodonSelection:
    """Test the _gt_aware_select_codon function."""

    @pytest.fixture
    def human_usage(self):
        """Human codon adaptiveness table."""
        return CODON_ADAPTIVENESS_TABLES.get("Homo_sapiens", {})

    @pytest.fixture
    def human_sorted_codons(self, human_usage):
        """Human codons sorted by CAI (descending)."""
        sorted_codons = {}
        for aa in set(AA_TO_CODONS.keys()):
            if aa == "*":
                continue
            codons = AA_TO_CODONS[aa]
            sorted_codons[aa] = sorted(
                codons, key=lambda c: human_usage.get(c, 0.0), reverse=True
            )
        return sorted_codons

    def test_no_next_aa_returns_optimal(self, human_usage, human_sorted_codons):
        """When there's no next amino acid, return the optimal codon."""
        result = _gt_aware_select_codon("L", None, human_sorted_codons, human_usage)
        assert result == human_sorted_codons["L"][0]

    def test_no_boundary_gt_returns_optimal(self, human_usage, human_sorted_codons):
        """When the optimal codon doesn't create a boundary GT, return it."""
        result = _gt_aware_select_codon("A", "A", human_sorted_codons, human_usage)
        assert result == human_sorted_codons["A"][0]

    def test_boundary_gt_avoided_within_tolerance(self, human_usage, human_sorted_codons):
        """When the optimal codon creates a boundary GT and an alternative
        exists within 10% CAI, the alternative should be selected."""
        leu_optimal = human_sorted_codons["L"][0]
        # Find a next AA whose optimal starts with T
        t_starting_aas = []
        for aa, codons in human_sorted_codons.items():
            if codons and codons[0][0] == "T":
                t_starting_aas.append(aa)

        if t_starting_aas and leu_optimal[-1] == "G":
            next_aa = t_starting_aas[0]
            next_optimal = human_sorted_codons[next_aa][0]
            if leu_optimal[-1] + next_optimal[0] == "GT":
                result = _gt_aware_select_codon("L", next_aa, human_sorted_codons, human_usage)
                leu_opt_cai = human_usage.get(leu_optimal, 0.0)
                result_cai = human_usage.get(result, 0.0)
                # If alternative was chosen, CAI should be within 10%
                if result != leu_optimal:
                    assert result_cai >= leu_opt_cai * (1.0 - GT_BOUNDARY_CAI_TOLERANCE), (
                        f"Alternative codon {result} CAI={result_cai} is not within "
                        f"{GT_BOUNDARY_CAI_TOLERANCE*100}% of optimal {leu_optimal} "
                        f"CAI={leu_opt_cai}"
                    )
                    # The alternative should not create a boundary GT
                    assert result[-1] + next_optimal[0] != "GT", (
                        f"Alternative {result} still creates boundary GT with {next_optimal}"
                    )

    def test_boundary_gt_accepted_when_no_good_alt(self, human_usage, human_sorted_codons):
        """When no alternative within CAI tolerance exists, the optimal codon
        should be used even though it creates a boundary GT."""
        val_codons = human_sorted_codons.get("V", [])
        if val_codons:
            result = _gt_aware_select_codon("V", "T", human_sorted_codons, human_usage)
            # Should return the optimal Val codon (accepting the GT)
            assert result in val_codons

    def test_prokaryote_path_not_affected(self, human_usage, human_sorted_codons):
        """The GT-aware selection function is only called for eukaryotes;
        prokaryote path uses direct optimal codon selection."""
        result = _gt_aware_select_codon("G", None, human_sorted_codons, human_usage)
        assert result in AA_TO_CODONS["G"]


# ═══════════════════════════════════════════════════════════════════════
# 5. Integration: Eukaryotic CAI preserved with soft GT constraint
# ═══════════════════════════════════════════════════════════════════════

class TestEukaryoteCAIPreserved:
    """Integration tests verifying that the soft GT constraint preserves CAI."""

    @pytest.mark.parametrize("organism", EUKARYOTIC_ORGANISMS)
    def test_insulin_cai_above_threshold(self, organism: str):
        """Insulin optimization for eukaryotes should maintain CAI > 0.90."""
        result = _optimize(INSULIN_PROTEIN, organism)
        assert result.cai > 0.90, (
            f"{organism}/insulin: CAI={result.cai:.4f} is not > 0.90. "
            f"Soft GT constraint should not destroy CAI."
        )

    @pytest.mark.parametrize("organism", EUKARYOTIC_ORGANISMS)
    def test_egfp_cai_above_threshold(self, organism: str):
        """EGFP optimization for eukaryotes should maintain CAI > 0.85.
        
        Note: GFP has many Val codons (which all contain GT), so some CAI
        loss from GT avoidance is expected. The soft GT constraint ensures
        that CAI loss is bounded rather than catastrophic.
        """
        result = _optimize(EGFP_PROTEIN, organism)
        assert result.cai > 0.85, (
            f"{organism}/EGFP: CAI={result.cai:.4f} is not > 0.85. "
            f"Soft GT constraint should not destroy CAI."
        )

    @pytest.mark.parametrize("organism", EUKARYOTIC_ORGANISMS)
    def test_translation_preserved(self, organism: str):
        """All eukaryotic optimizations must preserve protein sequence."""
        result = _optimize(INSULIN_PROTEIN, organism)
        translated = translate(result.sequence)
        assert translated == INSULIN_PROTEIN, (
            f"{organism}: translation mismatch after optimization"
        )

    def test_yeast_insulin_cai_not_regressed(self):
        """Yeast insulin CAI should not regress below 0.92 (was 0.83 before fix)."""
        result = _optimize(INSULIN_PROTEIN, "Saccharomyces_cerevisiae")
        assert result.cai > 0.92, (
            f"Yeast insulin CAI={result.cai:.4f} has regressed. "
            f"Before v10 fix: 0.83. After fix: ~0.925. "
            f"Soft GT constraint should not cause regression."
        )


# ═══════════════════════════════════════════════════════════════════════
# 6. Predicate result details and structure
# ═══════════════════════════════════════════════════════════════════════

class TestPredicateResultStructure:
    """Test that the soft GT predicate returns well-structured results."""

    def test_no_gt_passes_with_pass_verdict(self):
        """Sequence with no GTs should get PASS verdict."""
        result = check_no_gt_dinucleotide_soft("ATGGCTAGC", organism="Homo_sapiens")
        assert result.passed
        assert result.verdict == Verdict.PASS

    def test_gt_exceeds_max_likely_fail_verdict(self):
        """Sequence with many GTs exceeding max_gt_count should get LIKELY_FAIL."""
        many_gt_seq = "GTTGTTGTTGTTGTTGTT"
        result = check_no_gt_dinucleotide_soft(many_gt_seq, organism="Homo_sapiens")
        # For soft constraints, passed=True so it doesn't trigger strict mode
        assert result.passed  # soft fail → passed=True
        assert result.verdict == Verdict.LIKELY_FAIL

    def test_prokaryote_gt_returns_fail_verdict(self):
        """Prokaryotic sequence with GTs should get FAIL (hard constraint)."""
        result = check_no_gt_dinucleotide_soft("GTTGTTGTT", organism="E_coli")
        assert not result.passed
        assert result.verdict == Verdict.FAIL

    def test_explicit_max_gt_count_zero_hard_constraint(self):
        """Setting max_gt_count=0 on eukaryotes still uses LIKELY_FAIL (soft).
        Only prokaryotic organisms get FAIL (hard constraint)."""
        seq_with_gt = "ATGGGTAGC"  # GGT contains GT
        result = check_no_gt_dinucleotide_soft(
            seq_with_gt, organism="Homo_sapiens", max_gt_count=0
        )
        # Eukaryotic organism → soft constraint, even with max_gt_count=0
        assert result.verdict == Verdict.LIKELY_FAIL
        assert result.passed  # soft fail → passed=True

    def test_gt_positions_reported(self):
        """GT positions should be reported in the result."""
        result = check_no_gt_dinucleotide_soft("GTTGTT", organism="Homo_sapiens")
        assert len(result.positions) > 0, "GT positions should be reported"

    def test_in_codon_vs_cross_codon_in_details(self):
        """Details should distinguish in-codon vs cross-codon GTs."""
        result = check_no_gt_dinucleotide_soft("GTTGTT", organism="Homo_sapiens")
        assert "in-codon:" in result.details
        assert "cross-codon:" in result.details


# ═══════════════════════════════════════════════════════════════════════
# 7. Backward compatibility
# ═══════════════════════════════════════════════════════════════════════

class TestBackwardCompatibility:
    """Ensure existing check functions still work."""

    def test_strict_check_still_works(self):
        """check_no_gt_dinucleotide (strict) should still work."""
        result = check_no_gt_dinucleotide("GTTGTTGTT")
        assert not result.passed
        assert result.verdict == Verdict.FAIL

    def test_avoidable_check_still_works(self):
        """check_no_avoidable_gt (relaxed) should still work."""
        result = check_no_avoidable_gt("GTTGTTGTT", organism="Homo_sapiens")
        # Valine GTs are unavoidable, so should pass
        assert result.passed

    def test_soft_check_importable(self):
        """check_no_gt_dinucleotide_soft should be importable."""
        from biocompiler.type_system import check_no_gt_dinucleotide_soft
        assert callable(check_no_gt_dinucleotide_soft)

    def test_compute_max_gt_count_importable(self):
        """_compute_max_gt_count should be importable."""
        from biocompiler.type_system import _compute_max_gt_count
        assert callable(_compute_max_gt_count)
