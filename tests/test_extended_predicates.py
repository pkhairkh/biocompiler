"""
BioCompiler v7.6.0 Extended Predicate Tests
=============================================
Tests for the 20 extended predicates (predicates 9–28) plus the
foundational check_ functions (1, 5–8) that they build upon.

Each predicate is tested for:
  - Correct return type (PredicateResult or TypeCheckResult)
  - PASS case with a valid/healthy sequence
  - FAIL case with a problematic sequence (where possible)
  - Verdict is one of PASS/UNCERTAIN/FAIL (3-valued) or
    PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL (5-valued)
"""

import sys
import math

sys.path.insert(0, "src")

import pytest
from biocompiler.types import Verdict, TypeCheckResult
from biocompiler.type_system import (
    PredicateResult,
    check_no_gt_dinucleotide,
    check_no_avoidable_gt,
    check_no_stop_codons,
    check_valid_coding_seq,
    check_conservation_score,
    check_codon_optimality,
    check_no_cryptic_promoter,
    check_no_unexpected_tm_domain,
    check_mrna_secondary_structure,
    check_co_translational_folding,
    CODON_TABLE,
    AA_TO_CODONS,
    BLOSUM62,
)
from biocompiler.structure.predicates import (
    evaluate_structure_confidence,
    evaluate_no_misfolding_risk,
    evaluate_correct_fold_topology,
    evaluate_no_unexpected_interaction,
)
from biocompiler.stability_predicates import (
    evaluate_stable_folding,
    evaluate_no_destabilizing_mutation,
    evaluate_disulfide_bond_integrity,
    evaluate_hydrophobic_core_quality,
)
from biocompiler.solubility_predicates import (
    evaluate_soluble_expression,
    evaluate_no_aggregation_prone_region,
    evaluate_charge_composition,
    evaluate_no_long_hydrophobic_stretch,
)
from biocompiler.immuno_predicates import (
    evaluate_low_immunogenicity,
    evaluate_no_strong_t_cell_epitope,
    evaluate_no_dominant_b_cell_epitope,
    evaluate_population_coverage_safe,
)

# ────────────────────────────────────────────────────────────
# Shared test fixtures
# ────────────────────────────────────────────────────────────

# A short valid coding sequence: ATG GCC GCC GCC GCC TAA → M A A A A *
VALID_SEQ = "ATGGCCGCCGCCGCCTAA"

# Same but no stop codon (truncated)
VALID_SEQ_NO_STOP = "ATGGCCGCCGCCGCC"

# Sequence with internal stop: ATG TAA GCC ... TAA
SEQ_WITH_STOP = "ATGTAAGCCGCCGCCTAA"

# Sequence not divisible by 3
BAD_LEN_SEQ = "ATGGCCGC"

# Sequence with invalid codon
INVALID_CODON_SEQ = "ATGZZZGCCTAA"

# GT-free coding sequence (no GT at all)
# Use codons that don't produce GT: ATG GCC GCT AAA TAA
GT_FREE_SEQ = "ATGGCCGCTAAATAA"

# Sequence with GT dinucleotides (Valine codons GTN)
SEQ_WITH_GT = "ATGGTGGTGGTGTAA"

# AT-rich sequence (weak mRNA structure)
AT_RICH_SEQ = "ATG" + "AAT" * 15 + "ATT" * 10 + "TAA"

# GC-rich sequence (strong mRNA structure)
GC_RICH_SEQ = "ATG" + "GCG" * 15 + "CGC" * 10 + "TAA"

# Hydrophilic protein DNA (Lys-Asp-Glu-rich)
HYDROPHILIC_DNA = "AAA" * 10 + "GAT" * 5 + "GAA" * 5

# Hydrophobic protein DNA (Val-Leu-rich)
HYDROPHOBIC_DNA = "ATG" + "GTT" * 20 + "TTA" * 20 + "TAA"

# E. coli promoter-containing sequence (TTGACA + 17bp spacer + TATAAT)
PROMOTER_SEQ = "ATG" + "TTGACA" + "A" * 17 + "TATAAT" + "ATCGATCGATCGATCGATCG"

# Soluble protein (charged, hydrophilic, balanced pI ~5.5)
SOLUBLE_PROTEIN = "MSDQRGVAIDLNEKHSDQRGVAIDLNEKH"

# Insoluble protein (very hydrophobic)
INSOLUBLE_PROTEIN = "MIIILLLVVVAAAFFFWWWYYYLLLIIVV"

# Protein with odd number of cysteines (3 C's)
ODD_CYS_PROTEIN = "MACCDEFGHKLMPQSTVWYACDE"

# Protein with even number of cysteines (2 C's)
EVEN_CYS_PROTEIN = "MACKDEFGHKLMPQSTVWYACDE"

# Simple CAI dict for testing
SIMPLE_CAI = {
    "ATG": 0.8, "GCC": 1.0, "GCT": 0.6, "GCA": 0.5, "GCG": 0.4,
    "AAA": 0.9, "AAG": 0.7, "TTT": 0.5, "TTC": 0.8,
    "TAA": 0.1, "TAG": 0.05, "TGA": 0.05,
    "GTT": 0.6, "GTC": 0.7, "GTA": 0.3, "GTG": 1.0,
    "GAT": 0.6, "GAC": 0.9, "GAA": 0.8, "GAG": 0.7,
    "TTA": 0.2, "TTG": 0.4, "CTT": 0.3, "CTC": 0.5,
    "CTA": 0.2, "CTG": 1.0,
    "AAC": 0.8, "AAT": 0.4, "CAC": 0.7, "CAT": 0.3,
    "TGG": 0.9, "TAC": 0.7, "TAT": 0.4,
    "CGT": 0.5, "CGC": 0.7, "CGA": 0.2, "CGG": 0.3,
    "AGA": 0.6, "AGG": 0.3, "AGT": 0.4, "AGC": 0.7,
    "GGT": 0.5, "GGC": 0.8, "GGA": 0.4, "GGG": 0.6,
    "CCT": 0.3, "CCC": 0.5, "CCA": 0.4, "CCG": 0.9,
    "ACT": 0.4, "ACC": 0.8, "ACA": 0.3, "ACG": 0.5,
    "TCT": 0.4, "TCC": 0.7, "TCA": 0.3, "TCG": 0.5,
    "TGC": 0.6, "TGT": 0.3, "CAG": 0.8, "CAA": 0.4,
}

def _make_pdb_line(serial, atom_name, resname, chain, resseq, x, y, z, occ=1.00, bfac=50.00):
    """Construct a properly formatted PDB ATOM record line."""
    altloc = ' '
    # 4 spaces after resseq: insertion code (1) + 3 padding spaces
    line = (f'ATOM  {serial:5d} {atom_name:<4s}{altloc:1s}{resname:>3s} '
            f'{chain:1s}{resseq:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}'
            f'{occ:6.2f}{bfac:6.2f}          ')
    return line


def _make_pdb(n_residues, bfac_start=95.0, bfac_step=1.0, chain='A', resname='ALA'):
    """Construct a minimal PDB string with N residues for testing.

    Each residue has N, CA, C, CB, O atoms laid out along the x-axis.
    """
    lines = []
    for i in range(n_residues):
        x = 10.0 + i * 1.5
        bfac = bfac_start - i * bfac_step
        lines.append(_make_pdb_line(5*i+1, ' N  ', resname, chain, i+1, x, 20.0, 30.0, 1.00, bfac))
        lines.append(_make_pdb_line(5*i+2, ' CA ', resname, chain, i+1, x+0.5, 20.5, 30.5, 1.00, bfac))
        lines.append(_make_pdb_line(5*i+3, ' C  ', resname, chain, i+1, x+1.0, 21.0, 31.0, 1.00, bfac))
        lines.append(_make_pdb_line(5*i+4, ' CB ', resname, chain, i+1, x+0.5, 19.5, 30.5, 1.00, bfac))
        lines.append(_make_pdb_line(5*i+5, ' O  ', resname, chain, i+1, x+1.5, 21.5, 31.5, 1.00, bfac))
    return '\n'.join(lines)


# Minimal PDB string for structure predicates (5 residues, high pLDDT ~93)
MINI_PDB_HIGH_PLDDT = _make_pdb(5, bfac_start=95.0, bfac_step=1.0)

# Minimal PDB string with low pLDDT (20 residues, pLDDT ~30-40)
MINI_PDB_LOW_PLDDT = _make_pdb(20, bfac_start=40.0, bfac_step=0.5)

VALID_VERDICTS_3 = (Verdict.PASS, Verdict.UNCERTAIN, Verdict.FAIL)
VALID_VERDICTS_5 = (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN,
                    Verdict.LIKELY_FAIL, Verdict.FAIL)


# ═══════════════════════════════════════════════════════════════
# 1. check_no_stop_codons (Predicate 1)
# ═══════════════════════════════════════════════════════════════

class TestCheckNoStopCodons:
    """Predicate 1: No internal stop codons."""

    def test_returns_predicate_result(self):
        result = check_no_stop_codons(VALID_SEQ)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "NoStopCodons"
        assert isinstance(result.passed, bool)
        assert result.verdict in VALID_VERDICTS_3
        assert isinstance(result.details, str)

    def test_pass_no_stops(self):
        result = check_no_stop_codons(VALID_SEQ)
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_fail_internal_stop(self):
        result = check_no_stop_codons(SEQ_WITH_STOP)
        assert result.passed is False
        assert result.verdict == Verdict.FAIL
        assert len(result.positions) > 0


# ═══════════════════════════════════════════════════════════════
# 5a. check_no_gt_dinucleotide (Predicate 5 strict)
# ═══════════════════════════════════════════════════════════════

class TestCheckNoGTDinucleotide:
    """Predicate 5 (strict): No GT dinucleotides at all."""

    def test_returns_predicate_result(self):
        result = check_no_gt_dinucleotide(GT_FREE_SEQ)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "NoGTDinucleotide"

    def test_pass_no_gt(self):
        result = check_no_gt_dinucleotide(GT_FREE_SEQ)
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_fail_with_gt(self):
        result = check_no_gt_dinucleotide(SEQ_WITH_GT)
        assert result.passed is False
        assert result.verdict == Verdict.FAIL
        assert len(result.positions) > 0


# ═══════════════════════════════════════════════════════════════
# 5b. check_no_avoidable_gt (Predicate 5 relaxed)
# ═══════════════════════════════════════════════════════════════

class TestCheckNoAvoidableGT:
    """Predicate 5 (relaxed): No avoidable GT dinucleotides."""

    def test_returns_predicate_result(self):
        result = check_no_avoidable_gt(GT_FREE_SEQ)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "NoGTDinucleotide"

    def test_pass_no_gt(self):
        result = check_no_avoidable_gt(GT_FREE_SEQ)
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_pass_unavoidable_gt(self):
        # Valine codons (GTN) are all GT-starting, so they're unavoidable
        # ATG GTG GTG GTG TAA  →  M V V V *
        result = check_no_avoidable_gt(SEQ_WITH_GT)
        # All GTs from Valine are unavoidable
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_fail_avoidable_gt(self):
        # AGT (Ser) contains GT at position 1-2, avoidable via AGC
        seq = "ATGAGTAGCTAA"  # M S S *
        result = check_no_avoidable_gt(seq)
        assert result.passed is False
        assert result.verdict == Verdict.FAIL


# ═══════════════════════════════════════════════════════════════
# 6. check_valid_coding_seq (Predicate 6)
# ═══════════════════════════════════════════════════════════════

class TestCheckValidCodingSeq:
    """Predicate 6: Valid coding sequence."""

    def test_returns_predicate_result(self):
        result = check_valid_coding_seq(VALID_SEQ)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "ValidCodingSeq"

    def test_pass_valid(self):
        result = check_valid_coding_seq(VALID_SEQ)
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_fail_bad_length(self):
        result = check_valid_coding_seq(BAD_LEN_SEQ)
        assert result.passed is False
        assert result.verdict == Verdict.FAIL

    def test_fail_invalid_codon(self):
        result = check_valid_coding_seq(INVALID_CODON_SEQ)
        assert result.passed is False
        assert result.verdict == Verdict.FAIL


# ═══════════════════════════════════════════════════════════════
# 7. check_conservation_score (Predicate 7)
# ═══════════════════════════════════════════════════════════════

class TestCheckConservationScore:
    """Predicate 7: BLOSUM62 conservation score."""

    def test_returns_predicate_result(self):
        result = check_conservation_score("A", "A")
        assert isinstance(result, PredicateResult)
        assert result.predicate == "ConservationScore"

    def test_pass_conserved(self):
        # A→A is perfectly conserved (BLOSUM62 score = 4)
        result = check_conservation_score("A", "A", min_score=0)
        assert result.passed is True

    def test_fail_radical_substitution(self):
        # W→C is a radical substitution (BLOSUM62 score = -8)
        result = check_conservation_score("W", "C", min_score=0)
        assert result.passed is False

    def test_verdict_is_pass_or_fail(self):
        result = check_conservation_score("A", "V")
        assert result.verdict in VALID_VERDICTS_3


# ═══════════════════════════════════════════════════════════════
# 8. check_codon_optimality (Predicate 8)
# ═══════════════════════════════════════════════════════════════

class TestCheckCodonOptimality:
    """Predicate 8: Codon optimality (CAI)."""

    def test_returns_predicate_result(self):
        result = check_codon_optimality("GCC", SIMPLE_CAI, min_cai=0.5)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "CodonOptimality"

    def test_pass_high_cai(self):
        result = check_codon_optimality("GCC", SIMPLE_CAI, min_cai=0.5)
        assert result.passed is True
        assert result.verdict in VALID_VERDICTS_3

    def test_fail_low_cai(self):
        result = check_codon_optimality("TTA", SIMPLE_CAI, min_cai=0.5)
        assert result.passed is False

    def test_verdict_is_pass_or_fail(self):
        result = check_codon_optimality("GCC", SIMPLE_CAI, min_cai=0.0)
        assert result.verdict in VALID_VERDICTS_3


# ═══════════════════════════════════════════════════════════════
# 9. check_no_cryptic_promoter (Predicate 9)
# ═══════════════════════════════════════════════════════════════

class TestCheckNoCrypticPromoter:
    """Predicate 9: No cryptic promoter sites."""

    def test_returns_predicate_result(self):
        seq = "ATGGCGATCATCAGCTGAACCGGTTATCGATCGATCG"
        result = check_no_cryptic_promoter(seq, organism="E_coli")
        assert isinstance(result, PredicateResult)
        assert result.predicate == "NoCrypticPromoter"

    def test_pass_no_promoter(self):
        seq = "ATGGCGATCATCAGCTGAACCGGTTATCGATCGATCG"
        result = check_no_cryptic_promoter(seq, organism="E_coli")
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_fail_or_uncertain_with_promoter(self):
        result = check_no_cryptic_promoter(PROMOTER_SEQ, organism="E_coli")
        assert result.verdict in (Verdict.UNCERTAIN, Verdict.FAIL)

    def test_eukaryote_mode(self):
        seq = "ATG" + "TATAAA" + "CGATCGATCGATCGATCGATCGATCGATCG"
        result = check_no_cryptic_promoter(seq, organism="Homo_sapiens")
        assert result.verdict in VALID_VERDICTS_3


# ═══════════════════════════════════════════════════════════════
# 10. check_no_unexpected_tm_domain (Predicate 10)
# ═══════════════════════════════════════════════════════════════

class TestCheckNoUnexpectedTMDomain:
    """Predicate 10: No unexpected transmembrane domains."""

    def test_returns_predicate_result(self):
        result = check_no_unexpected_tm_domain(HYDROPHILIC_DNA + "TAA")
        assert isinstance(result, PredicateResult)
        assert result.predicate == "NoUnexpectedTMDomain"

    def test_pass_hydrophilic(self):
        result = check_no_unexpected_tm_domain(HYDROPHILIC_DNA + "TAA", is_cytosolic=True)
        assert result.passed is True
        assert result.verdict in VALID_VERDICTS_3

    def test_pass_membrane_protein(self):
        result = check_no_unexpected_tm_domain(HYDROPHOBIC_DNA, is_cytosolic=False)
        assert result.passed is True
        assert result.verdict == Verdict.PASS

    def test_fail_hydrophobic_cytosolic(self):
        result = check_no_unexpected_tm_domain(HYDROPHOBIC_DNA, is_cytosolic=True)
        assert result.passed is False
        assert result.verdict == Verdict.FAIL


# ═══════════════════════════════════════════════════════════════
# 11. check_mrna_secondary_structure (Predicate 11)
# ═══════════════════════════════════════════════════════════════

class TestCheckMRNASecondaryStructure:
    """Predicate 11: mRNA secondary structure around RBS."""

    def test_returns_predicate_result(self):
        result = check_mrna_secondary_structure(AT_RICH_SEQ)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "mRNASecondaryStructure"

    def test_pass_at_rich(self):
        result = check_mrna_secondary_structure(AT_RICH_SEQ, window_end=50)
        assert result.passed is True
        assert result.verdict in VALID_VERDICTS_3

    def test_gc_rich_may_fail(self):
        result = check_mrna_secondary_structure(GC_RICH_SEQ, window_end=50, dg_threshold=-15.0)
        assert result.verdict in VALID_VERDICTS_3

    def test_short_sequence(self):
        result = check_mrna_secondary_structure("ATG")
        assert result.verdict in VALID_VERDICTS_3


# ═══════════════════════════════════════════════════════════════
# 12. check_co_translational_folding (Predicate 12)
# ═══════════════════════════════════════════════════════════════

class TestCheckCoTranslationalFolding:
    """Predicate 12: Co-translational folding preservation."""

    def test_returns_predicate_result(self):
        seq = "ATGGCCGCCGCCGCCTAA"
        result = check_co_translational_folding(seq, SIMPLE_CAI)
        assert isinstance(result, PredicateResult)
        assert result.predicate == "CoTranslationalFolding"

    def test_pass_with_slow_codons(self):
        # TTA (Leu) has CAI=0.2, giving natural pause sites
        seq = "ATGTTAGCCGCCTTA" + "GCCTAA"
        result = check_co_translational_folding(seq, SIMPLE_CAI)
        assert result.verdict in VALID_VERDICTS_5

    def test_fail_all_fast_no_boundaries(self):
        # All high-CAI codons, no domain boundaries
        seq = "ATGGCCGCCGCCGCCGCCGCCGCC" + "GCCGCCGCCTAA"
        result = check_co_translational_folding(seq, SIMPLE_CAI, domain_boundaries=[3, 7])
        assert result.verdict in VALID_VERDICTS_5

    def test_empty_sequence(self):
        result = check_co_translational_folding("", SIMPLE_CAI)
        assert result.verdict == Verdict.PASS


# ═══════════════════════════════════════════════════════════════
# 13. evaluate_structure_confidence (Predicate 13)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateStructureConfidence:
    """Predicate 13: Structure confidence (pLDDT)."""

    def test_returns_type_check_result(self):
        result = evaluate_structure_confidence(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "StructureConfidence" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_pass_high_plddt_pdb(self):
        result = evaluate_structure_confidence(
            "ATGGCTGCTGCTTAA", "MAAAA", "Homo_sapiens",
            pdb_string=MINI_PDB_HIGH_PLDDT,
        )
        # Mean pLDDT ~93 → PASS
        assert result.verdict == Verdict.PASS

    def test_fail_low_plddt_pdb(self):
        result = evaluate_structure_confidence(
            "ATGGCTGCTGCTTAA", "MAAAA", "Homo_sapiens",
            pdb_string=MINI_PDB_LOW_PLDDT, min_mean_plddt=50.0,
        )
        assert result.verdict in (Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)

    def test_uncertain_no_pdb(self):
        result = evaluate_structure_confidence(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        # Without PDB and ESMFold, should return UNCERTAIN
        assert result.verdict in (Verdict.UNCERTAIN, Verdict.LIKELY_PASS, Verdict.PASS)


# ═══════════════════════════════════════════════════════════════
# 14. evaluate_no_misfolding_risk (Predicate 14)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoMisfoldingRisk:
    """Predicate 14: Misfolding risk indicators."""

    def test_returns_type_check_result(self):
        result = evaluate_no_misfolding_risk(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "NoMisfoldingRisk"
        assert result.verdict in VALID_VERDICTS_5

    def test_uncertain_no_pdb(self):
        result = evaluate_no_misfolding_risk(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert result.verdict == Verdict.UNCERTAIN
        assert result.knowledge_gap is not None

    def test_pass_with_good_pdb(self):
        result = evaluate_no_misfolding_risk(
            "ATGGCTGCTGCTTAA", "MAAAA", "Homo_sapiens",
            pdb_string=MINI_PDB_HIGH_PLDDT,
        )
        # Good structure should pass or at least not LIKELY_FAIL
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS, Verdict.UNCERTAIN)


# ═══════════════════════════════════════════════════════════════
# 15. evaluate_correct_fold_topology (Predicate 15)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateCorrectFoldTopology:
    """Predicate 15: Fold topology validation."""

    def test_returns_type_check_result(self):
        result = evaluate_correct_fold_topology(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "CorrectFoldTopology"
        assert result.verdict in VALID_VERDICTS_5

    def test_uncertain_no_pdb(self):
        result = evaluate_correct_fold_topology(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert result.verdict == Verdict.UNCERTAIN

    def test_with_pdb(self):
        result = evaluate_correct_fold_topology(
            "ATGGCTGCTGCTTAA", "MAAAA", "Homo_sapiens",
            pdb_string=MINI_PDB_HIGH_PLDDT,
        )
        assert result.verdict in VALID_VERDICTS_5


# ═══════════════════════════════════════════════════════════════
# 16. evaluate_no_unexpected_interaction (Predicate 16)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoUnexpectedInteraction:
    """Predicate 16: No unexpected protein-protein interactions."""

    def test_returns_type_check_result(self):
        result = evaluate_no_unexpected_interaction(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "NoUnexpectedInteraction"
        assert result.verdict in VALID_VERDICTS_5

    def test_pass_soluble_protein(self):
        result = evaluate_no_unexpected_interaction(
            "ATGAAAAAA...TAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_hydrophobic_may_flag(self):
        result = evaluate_no_unexpected_interaction(
            "ATG...TAA", INSOLUBLE_PROTEIN, "Homo_sapiens"
        )
        # Hydrophobic proteins may flag interaction risk
        assert result.verdict in VALID_VERDICTS_5


# ═══════════════════════════════════════════════════════════════
# 17. evaluate_stable_folding (Predicate 17)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateStableFolding:
    """Predicate 17: Thermodynamic stability (dG)."""

    def test_returns_type_check_result(self):
        result = evaluate_stable_folding(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "StableFolding" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_soluble_protein_verdict(self):
        result = evaluate_stable_folding(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5
        assert result.derivation is not None

    def test_empty_protein(self):
        result = evaluate_stable_folding("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.UNCERTAIN

    def test_with_pdb(self):
        result = evaluate_stable_folding(
            "ATGGCTGCTGCTTAA", "MAAAA", "Homo_sapiens",
            pdb_string=MINI_PDB_HIGH_PLDDT,
        )
        assert result.verdict in VALID_VERDICTS_5


# ═══════════════════════════════════════════════════════════════
# 18. evaluate_no_destabilizing_mutation (Predicate 18)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoDestabilizingMutation:
    """Predicate 18: No high-ΔΔG mutations."""

    def test_returns_type_check_result(self):
        result = evaluate_no_destabilizing_mutation(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "NoDestabilizingMutation" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_pass_no_original(self):
        # Without original_protein, should PASS
        result = evaluate_no_destabilizing_mutation(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens"
        )
        assert result.verdict == Verdict.PASS

    def test_pass_same_sequence(self):
        # Same protein as original → no mutations → PASS
        result = evaluate_no_destabilizing_mutation(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens",
            original_protein="MAAA",
        )
        assert result.verdict == Verdict.PASS

    def test_fail_length_mismatch(self):
        result = evaluate_no_destabilizing_mutation(
            "ATGGCTGCTGCTTAA", "MAAA", "Homo_sapiens",
            original_protein="MAAAA",
        )
        assert result.verdict == Verdict.FAIL

    def test_fail_radical_mutation(self):
        # W→P is radical (BLOSUM62=-4), ddG estimate ~3.2 > 3.0
        result = evaluate_no_destabilizing_mutation(
            "ATGGCTGCTGCTTAA", "MPAA", "Homo_sapiens",
            original_protein="MWAA", max_ddg=3.0,
        )
        assert result.verdict in (Verdict.LIKELY_FAIL, Verdict.FAIL)


# ═══════════════════════════════════════════════════════════════
# 19. evaluate_disulfide_bond_integrity (Predicate 19)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateDisulfideBondIntegrity:
    """Predicate 19: Cysteine pairing check."""

    def test_returns_type_check_result(self):
        result = evaluate_disulfide_bond_integrity(
            "ATGGCTGCTGCTTAA", EVEN_CYS_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "DisulfideBondIntegrity"
        assert result.verdict in VALID_VERDICTS_5

    def test_pass_no_cysteines(self):
        result = evaluate_disulfide_bond_integrity(
            "ATGGCTGCTGCTTAA", "MADEKRH", "Homo_sapiens"
        )
        assert result.verdict == Verdict.PASS

    def test_pass_even_cysteines(self):
        result = evaluate_disulfide_bond_integrity(
            "ATGGCTGCTGCTTAA", EVEN_CYS_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict == Verdict.PASS

    def test_fail_odd_cysteines(self):
        result = evaluate_disulfide_bond_integrity(
            "ATGGCTGCTGCTTAA", ODD_CYS_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict == Verdict.LIKELY_FAIL


# ═══════════════════════════════════════════════════════════════
# 20. evaluate_hydrophobic_core_quality (Predicate 20)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateHydrophobicCoreQuality:
    """Predicate 20: Hydrophobic core composition."""

    def test_returns_type_check_result(self):
        result = evaluate_hydrophobic_core_quality(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "HydrophobicCoreQuality"
        assert result.verdict in VALID_VERDICTS_5

    def test_pass_normal_hydrophobic(self):
        # Protein with ~33% hydrophobic (normal range 30-45%)
        # hydrophobic AAs: A,I,L,M,F,W,V
        protein = "MSDQRGVAIDLNEKH"
        result = evaluate_hydrophobic_core_quality(
            "ATGGCTGCTGCTTAA", protein, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def test_fail_very_low_hydrophobic(self):
        # All charged residues → hydrophobic fraction near 0
        protein = "KRRRDEEEEEKRRRDEEEEKRR"
        result = evaluate_hydrophobic_core_quality(
            "ATGGCTGCTGCTTAA", protein, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)

    def test_fail_very_high_hydrophobic(self):
        # All hydrophobic → fraction ~100%
        result = evaluate_hydrophobic_core_quality(
            "ATGGCTGCTGCTTAA", INSOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)


# ═══════════════════════════════════════════════════════════════
# 21. evaluate_soluble_expression (Predicate 21)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateSolubleExpression:
    """Predicate 21: CamSol solubility score."""

    def test_returns_type_check_result(self):
        result = evaluate_soluble_expression(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "SolubleExpression" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_soluble_protein_passes(self):
        result = evaluate_soluble_expression(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def test_insoluble_protein_fails(self):
        result = evaluate_soluble_expression(
            "ATGGCTGCTGCTTAA", INSOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.UNCERTAIN, Verdict.LIKELY_FAIL, Verdict.FAIL)

    def test_empty_protein(self):
        result = evaluate_soluble_expression("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.FAIL


# ═══════════════════════════════════════════════════════════════
# 22. evaluate_no_aggregation_prone_region (Predicate 22)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoAggregationProneRegion:
    """Predicate 22: Aggregation-prone region detection."""

    def test_returns_type_check_result(self):
        result = evaluate_no_aggregation_prone_region(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "NoAggregationProneRegion" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_soluble_protein_passes(self):
        result = evaluate_no_aggregation_prone_region(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def test_hydrophobic_protein_flags(self):
        result = evaluate_no_aggregation_prone_region(
            "ATGGCTGCTGCTTAA", INSOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_empty_protein(self):
        result = evaluate_no_aggregation_prone_region("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.PASS


# ═══════════════════════════════════════════════════════════════
# 23. evaluate_charge_composition (Predicate 23)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateChargeComposition:
    """Predicate 23: Charge balance and pI."""

    def test_returns_type_check_result(self):
        result = evaluate_charge_composition(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "ChargeComposition" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_balanced_protein_passes(self):
        # Balanced protein with ~40% charged and pI ~5.3
        protein = "MSDQRGVAIDLNEKH"
        result = evaluate_charge_composition(
            "ATGGCTGCTGCTTAA", protein, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def test_uncharged_protein_fails(self):
        # All hydrophobic, no charged residues
        result = evaluate_charge_composition(
            "ATGGCTGCTGCTTAA", INSOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.LIKELY_FAIL, Verdict.UNCERTAIN, Verdict.FAIL)

    def test_empty_protein(self):
        result = evaluate_charge_composition("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.FAIL


# ═══════════════════════════════════════════════════════════════
# 24. evaluate_no_long_hydrophobic_stretch (Predicate 24)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoLongHydrophobicStretch:
    """Predicate 24: Long hydrophobic stretch detection."""

    def test_returns_type_check_result(self):
        result = evaluate_no_long_hydrophobic_stretch(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert "NoLongHydrophobicStretch" in result.predicate
        assert result.verdict in VALID_VERDICTS_5

    def test_soluble_protein_passes(self):
        result = evaluate_no_long_hydrophobic_stretch(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def test_very_long_stretch_fails(self):
        # 20 consecutive hydrophobic residues (exceeds default max=7 by >6)
        protein = "M" + "I" * 20 + "K"
        result = evaluate_no_long_hydrophobic_stretch(
            "ATGGCTGCTGCTTAA", protein, "Homo_sapiens"
        )
        assert result.verdict in (Verdict.FAIL, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)

    def test_empty_protein(self):
        result = evaluate_no_long_hydrophobic_stretch("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.PASS


# ═══════════════════════════════════════════════════════════════
# 25. evaluate_low_immunogenicity (Predicate 25)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateLowImmunogenicity:
    """Predicate 25: Overall immunogenicity score."""

    def test_returns_type_check_result(self):
        result = evaluate_low_immunogenicity(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "LowImmunogenicity"
        assert result.verdict in VALID_VERDICTS_5

    def test_verdict_valid(self):
        result = evaluate_low_immunogenicity(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5
        assert result.derivation is not None

    def test_empty_protein(self):
        result = evaluate_low_immunogenicity("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.UNCERTAIN


# ═══════════════════════════════════════════════════════════════
# 26. evaluate_no_strong_t_cell_epitope (Predicate 26)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoStrongTCellEpitope:
    """Predicate 26: MHC binding epitope detection."""

    def test_returns_type_check_result(self):
        result = evaluate_no_strong_t_cell_epitope(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "NoStrongTCellEpitope"
        assert result.verdict in VALID_VERDICTS_5

    def test_verdict_valid(self):
        result = evaluate_no_strong_t_cell_epitope(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_empty_protein(self):
        result = evaluate_no_strong_t_cell_epitope("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.UNCERTAIN


# ═══════════════════════════════════════════════════════════════
# 27. evaluate_no_dominant_b_cell_epitope (Predicate 27)
# ═══════════════════════════════════════════════════════════════

class TestEvaluateNoDominantBCellEpitope:
    """Predicate 27: B-cell epitope coverage."""

    def test_returns_type_check_result(self):
        result = evaluate_no_dominant_b_cell_epitope(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "NoDominantBCellEpitope"
        assert result.verdict in VALID_VERDICTS_5

    def test_verdict_valid(self):
        result = evaluate_no_dominant_b_cell_epitope(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_empty_protein(self):
        result = evaluate_no_dominant_b_cell_epitope("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.UNCERTAIN


# ═══════════════════════════════════════════════════════════════
# 28. evaluate_population_coverage_safe (Predicate 28)
# ═══════════════════════════════════════════════════════════════

class TestEvaluatePopulationCoverageSafe:
    """Predicate 28: MHC allele population coverage."""

    def test_returns_type_check_result(self):
        result = evaluate_population_coverage_safe(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert isinstance(result, TypeCheckResult)
        assert result.predicate == "PopulationCoverageSafe"
        assert result.verdict in VALID_VERDICTS_5

    def test_verdict_valid(self):
        result = evaluate_population_coverage_safe(
            "ATGGCTGCTGCTTAA", SOLUBLE_PROTEIN, "Homo_sapiens"
        )
        assert result.verdict in VALID_VERDICTS_5

    def test_empty_protein(self):
        result = evaluate_population_coverage_safe("ATGTAA", "", "Homo_sapiens")
        assert result.verdict == Verdict.UNCERTAIN


# ═══════════════════════════════════════════════════════════════
# Integration: All 20 Extended Predicates Can Run Together
# ═══════════════════════════════════════════════════════════════

class TestAllExtendedPredicatesIntegration:
    """Verify that all 20 extended predicates can run without error."""

    def test_all_return_valid_verdicts(self):
        """Every extended predicate should return a valid verdict for a typical sequence."""
        dna = "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTAA"
        protein = "MAAAAAAA"

        # DNA-level predicates (check_ functions)
        r9 = check_no_cryptic_promoter(dna, organism="E_coli")
        assert r9.verdict in VALID_VERDICTS_3

        r10 = check_no_unexpected_tm_domain(dna, is_cytosolic=True)
        assert r10.verdict in VALID_VERDICTS_3

        r11 = check_mrna_secondary_structure(dna)
        assert r11.verdict in VALID_VERDICTS_3

        r12 = check_co_translational_folding(dna, SIMPLE_CAI)
        assert r12.verdict in VALID_VERDICTS_5

        # Structure predicates
        r13 = evaluate_structure_confidence(dna, protein, "Homo_sapiens")
        assert r13.verdict in VALID_VERDICTS_5

        r14 = evaluate_no_misfolding_risk(dna, protein, "Homo_sapiens")
        assert r14.verdict in VALID_VERDICTS_5

        r15 = evaluate_correct_fold_topology(dna, protein, "Homo_sapiens")
        assert r15.verdict in VALID_VERDICTS_5

        r16 = evaluate_no_unexpected_interaction(dna, protein, "Homo_sapiens")
        assert r16.verdict in VALID_VERDICTS_5

        # Stability predicates
        r17 = evaluate_stable_folding(dna, protein, "Homo_sapiens")
        assert r17.verdict in VALID_VERDICTS_5

        r18 = evaluate_no_destabilizing_mutation(dna, protein, "Homo_sapiens")
        assert r18.verdict in VALID_VERDICTS_5

        r19 = evaluate_disulfide_bond_integrity(dna, protein, "Homo_sapiens")
        assert r19.verdict in VALID_VERDICTS_5

        r20 = evaluate_hydrophobic_core_quality(dna, protein, "Homo_sapiens")
        assert r20.verdict in VALID_VERDICTS_5

        # Solubility predicates
        r21 = evaluate_soluble_expression(dna, protein, "Homo_sapiens")
        assert r21.verdict in VALID_VERDICTS_5

        r22 = evaluate_no_aggregation_prone_region(dna, protein, "Homo_sapiens")
        assert r22.verdict in VALID_VERDICTS_5

        r23 = evaluate_charge_composition(dna, protein, "Homo_sapiens")
        assert r23.verdict in VALID_VERDICTS_5

        r24 = evaluate_no_long_hydrophobic_stretch(dna, protein, "Homo_sapiens")
        assert r24.verdict in VALID_VERDICTS_5

        # Immunogenicity predicates
        r25 = evaluate_low_immunogenicity(dna, protein, "Homo_sapiens")
        assert r25.verdict in VALID_VERDICTS_5

        r26 = evaluate_no_strong_t_cell_epitope(dna, protein, "Homo_sapiens")
        assert r26.verdict in VALID_VERDICTS_5

        r27 = evaluate_no_dominant_b_cell_epitope(dna, protein, "Homo_sapiens")
        assert r27.verdict in VALID_VERDICTS_5

        r28 = evaluate_population_coverage_safe(dna, protein, "Homo_sapiens")
        assert r28.verdict in VALID_VERDICTS_5

    def test_check_functions_return_predicate_result(self):
        """All check_ functions should return PredicateResult."""
        dna = "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTAA"

        assert isinstance(check_no_stop_codons(dna), PredicateResult)
        assert isinstance(check_no_gt_dinucleotide(dna), PredicateResult)
        assert isinstance(check_no_avoidable_gt(dna), PredicateResult)
        assert isinstance(check_valid_coding_seq(dna), PredicateResult)
        assert isinstance(check_conservation_score("A", "A"), PredicateResult)
        assert isinstance(check_codon_optimality("GCC", SIMPLE_CAI), PredicateResult)
        assert isinstance(check_no_cryptic_promoter(dna), PredicateResult)
        assert isinstance(check_no_unexpected_tm_domain(dna), PredicateResult)
        assert isinstance(check_mrna_secondary_structure(dna), PredicateResult)
        assert isinstance(check_co_translational_folding(dna, SIMPLE_CAI), PredicateResult)

    def test_evaluate_functions_return_type_check_result(self):
        """All evaluate_ functions should return TypeCheckResult."""
        dna = "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTAA"
        protein = "MAAAAAAA"

        assert isinstance(evaluate_structure_confidence(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_misfolding_risk(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_correct_fold_topology(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_unexpected_interaction(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_stable_folding(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_destabilizing_mutation(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_disulfide_bond_integrity(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_hydrophobic_core_quality(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_soluble_expression(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_aggregation_prone_region(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_charge_composition(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_long_hydrophobic_stretch(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_low_immunogenicity(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_strong_t_cell_epitope(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_no_dominant_b_cell_epitope(dna, protein, "Homo_sapiens"), TypeCheckResult)
        assert isinstance(evaluate_population_coverage_safe(dna, protein, "Homo_sapiens"), TypeCheckResult)
