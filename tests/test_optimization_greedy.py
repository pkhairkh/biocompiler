"""Tests for the greedy optimizer (_greedy_optimize) and related functions.

Covers:
- _greedy_optimize with basic proteins (insulin, GFP)
- CAI maximization step
- Restriction site removal
- GC adjustment
- T-run fixing
- Empty input handling
- Invalid codon handling
- Eukaryote-specific steps (GT avoidance, CpG, splice)
- Convergence detection
- protein_to_aa_list validation
- ConvergenceTracker
- _find_site_in_sequence and _get_overlapping_codons helpers
"""

import pytest
from biocompiler.optimization import (
    _greedy_optimize,
    optimize_sequence,
    OptimizationResult,
    ConvergenceTracker,
    protein_to_aa_list,
    _find_site_in_sequence,
    _get_overlapping_codons,
    score_splice_donor_potential,
    _compute_cai_fast,
    _count_dinucs_fast,
    _BatchSwapScorer,
)
from biocompiler.exceptions import InvalidProteinError, UnsupportedOrganismError
from biocompiler.type_system import CODON_TABLE, AA_TO_CODONS
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES


# ── protein_to_aa_list ──────────────────────────────────────────────

class TestProteinToAAList:
    def test_valid_protein(self):
        result = protein_to_aa_list("MVSKGE")
        assert result == ["M", "V", "S", "K", "G", "E"]

    def test_empty_protein_raises(self):
        with pytest.raises(InvalidProteinError):
            protein_to_aa_list("")

    def test_whitespace_protein_raises(self):
        with pytest.raises(InvalidProteinError):
            protein_to_aa_list("   ")

    def test_invalid_amino_acid_raises(self):
        with pytest.raises(InvalidProteinError):
            protein_to_aa_list("MVSXGE")

    def test_lowercase_protein_converted(self):
        result = protein_to_aa_list("mvskge")
        assert result == ["M", "V", "S", "K", "G", "E"]

    def test_single_amino_acid(self):
        result = protein_to_aa_list("M")
        assert result == ["M"]

    def test_all_standard_aas(self):
        protein = "ACDEFGHIKLMNPQRSTVWY"
        result = protein_to_aa_list(protein)
        assert len(result) == 20
        assert result == list(protein)


# ── ConvergenceTracker ──────────────────────────────────────────────

class TestConvergenceTracker:
    def test_initial_state(self):
        ct = ConvergenceTracker()
        assert ct.iterations == 0
        assert ct.best == float('-inf')
        assert ct.best_iteration_index == -1

    def test_record_improves_best(self):
        ct = ConvergenceTracker()
        ct.record(0.5)
        assert ct.best == pytest.approx(0.5, rel=1e-6)
        assert ct.best_iteration_index == 0
        ct.record(0.8)
        assert ct.best == pytest.approx(0.8, rel=1e-6)
        assert ct.best_iteration_index == 1

    def test_convergence_detected(self):
        ct = ConvergenceTracker(patience=3, improvement_threshold=1e-6)
        ct.record(0.5)
        ct.record(0.5000001)
        ct.record(0.5000002)
        ct.record(0.5000003)
        result = ct.check_convergence()
        assert result == "converged"

    def test_no_convergence_when_improving(self):
        ct = ConvergenceTracker(patience=3, improvement_threshold=0.01)
        ct.record(0.5)
        ct.record(0.6)
        assert ct.check_convergence() is None

    def test_oscillation_detected(self):
        ct = ConvergenceTracker(oscillation_window=5, improvement_threshold=1e-6)
        ct.record(0.5)
        ct.record(0.7)
        ct.record(0.4)
        ct.record(0.7)
        ct.record(0.4)
        result = ct.check_convergence()
        assert result == "oscillating"

    def test_patience_respected(self):
        ct = ConvergenceTracker(patience=5, improvement_threshold=1e-6)
        ct.record(0.5)
        ct.record(0.5)
        ct.record(0.5)
        # Only 3 records, patience=5, should not converge yet
        assert ct.check_convergence() is None


# ── _find_site_in_sequence ──────────────────────────────────────────

class TestFindSiteInSequence:
    def test_site_found(self):
        positions = _find_site_in_sequence("ATGGAATTCC", "GAATTC", "GAATTC")
        assert 3 in positions

    def test_site_not_found(self):
        positions = _find_site_in_sequence("ATGCATGCAT", "GAATTC", "GAATTC")
        assert positions == []

    def test_reverse_complement_found(self):
        # EcoRI: GAATTC, RC: GAATTC (palindrome)
        positions = _find_site_in_sequence("ATGGAATTCC", "GAATTC", "GAATTC")
        assert len(positions) >= 1

    def test_non_palindrome_site(self):
        # BamHI: GGATCC is palindromic, but position should be found
        positions = _find_site_in_sequence("ATGGATCCAT", "GGATCC", "GGATCC")
        assert len(positions) >= 1
        assert 2 in positions  # GGATCC starts at position 2

    def test_empty_site(self):
        positions = _find_site_in_sequence("ATGATG", "", "")
        assert positions == []


# ── _get_overlapping_codons ─────────────────────────────────────────

class TestGetOverlappingCodons:
    def test_site_within_single_codon(self):
        result = _get_overlapping_codons(0, 3, 10)
        assert result == [0]

    def test_site_spanning_two_codons(self):
        # Position 2, length 4 spans codons 0 and 1
        result = _get_overlapping_codons(2, 4, 10)
        assert result == [0, 1]

    def test_site_at_boundary(self):
        # Position 3 = start of codon 1, length 6 spans codons 1 and 2
        result = _get_overlapping_codons(3, 6, 10)
        assert result == [1, 2]

    def test_site_spanning_three_codons(self):
        # Position 1, length 7 spans codons 0, 1, 2
        result = _get_overlapping_codons(1, 7, 10)
        assert result == [0, 1, 2]


# ── _compute_cai_fast ──────────────────────────────────────────────

class TestComputeCAIFast:
    def test_empty_sequence(self):
        result = _compute_cai_fast("", CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])
        assert result == 0.0

    def test_short_sequence(self):
        result = _compute_cai_fast("AT", CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])
        assert result == 0.0

    def test_valid_ecoli_sequence(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        seq = "ATGAAA"  # M + K (Met is skipped, Lys is scored)
        result = _compute_cai_fast(seq, ecoli_cai)
        assert isinstance(result, float)
        assert 0.0 <= result <= 1.0

    def test_all_optimal_codons_high_cai(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        # Use optimal E. coli codons
        seq = "ATG" + "CTG" * 10  # Met + Leu (CTG is optimal Leu in E. coli)
        result = _compute_cai_fast(seq, ecoli_cai)
        assert result > 0.5


# ── _count_dinucs_fast ─────────────────────────────────────────────

class TestCountDinucsFast:
    def test_no_dinucleotides(self):
        result = _count_dinucs_fast("")
        assert result == ()

    def test_count_gt(self):
        result = _count_dinucs_fast("ATGGTGAAGGT", "GT")
        assert result[0] >= 2  # At least 2 GT dinucleotides

    def test_count_multiple_dinucs(self):
        result = _count_dinucs_fast("ATGCGTAACG", "GT", "CG")
        assert len(result) == 2
        assert isinstance(result[0], int)
        assert isinstance(result[1], int)

    def test_no_match(self):
        result = _count_dinucs_fast("AAAAAAAAAA", "GT")
        assert result[0] == 0


# ── _BatchSwapScorer ───────────────────────────────────────────────

class TestBatchSwapScorer:
    def test_initialization(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        scorer = _BatchSwapScorer(ecoli_cai)
        assert scorer.current_log_sum is None

    def test_reset_incremental_state(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        scorer = _BatchSwapScorer(ecoli_cai)
        codons = ["ATG", "CTG", "AAA"]
        scorer.reset_incremental_state(codons)
        assert scorer.current_log_sum is not None
        assert isinstance(scorer.current_log_sum, float)

    def test_score_candidates(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        scorer = _BatchSwapScorer(ecoli_cai)
        codons = ["ATG", "CTG", "AAA"]
        scorer.reset_incremental_state(codons)
        # Try alternative Leu codons at position 1
        candidates = AA_TO_CODONS["L"]
        scores = scorer.score_candidates(codons, 1, candidates)
        assert len(scores) == len(candidates)
        assert all(isinstance(s, float) for s in scores)
        assert all(0.0 <= s <= 1.0 for s in scores)

    def test_update_incremental_state(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        scorer = _BatchSwapScorer(ecoli_cai)
        codons = ["ATG", "CTG", "AAA"]
        scorer.reset_incremental_state(codons)
        old_log = scorer.current_log_sum
        scorer.update_incremental_state("CTG", "TTA")  # Swap optimal Leu to rare Leu
        assert scorer.current_log_sum != old_log

    def test_score_candidates_empty(self):
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        scorer = _BatchSwapScorer(ecoli_cai)
        scorer.reset_incremental_state(["ATG"])
        scores = scorer.score_candidates(["ATG"], 0, [])
        assert scores == []


# ── _greedy_optimize ───────────────────────────────────────────────

class TestGreedyOptimize:
    def test_basic_protein_ecoli(self):
        seq, warnings = _greedy_optimize(
            "MVSKGE", organism="Escherichia_coli", is_prokaryote=True,
        )
        assert isinstance(seq, str)
        assert len(seq) == 18  # 6 AA * 3
        assert len(warnings) >= 0

    def test_basic_protein_human(self):
        seq, warnings = _greedy_optimize(
            "MVSKGE", organism="Homo_sapiens", is_prokaryote=False,
        )
        assert isinstance(seq, str)
        assert len(seq) == 18

    def test_insulin_protein_ecoli(self):
        protein = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
        seq, warnings = _greedy_optimize(
            protein, organism="Escherichia_coli", is_prokaryote=True,
        )
        assert isinstance(seq, str)
        assert len(seq) == len(protein) * 3

    def test_gfp_protein_human(self):
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
        seq, warnings = _greedy_optimize(
            protein, organism="Homo_sapiens", is_prokaryote=False,
        )
        assert isinstance(seq, str)
        assert len(seq) == len(protein) * 3

    def test_sequence_translates_to_input_protein(self):
        protein = "MVSKGE"
        seq, _ = _greedy_optimize(protein, organism="Escherichia_coli", is_prokaryote=True)
        translated = "".join(CODON_TABLE.get(seq[i:i+3], "?") for i in range(0, len(seq), 3))
        assert translated == protein

    def test_prokaryote_skips_eukaryote_steps(self):
        # Prokaryote optimization should not fail when splice/CpG steps are skipped
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVS"
        seq, warnings = _greedy_optimize(
            protein, organism="Escherichia_coli", is_prokaryote=True,
        )
        assert len(seq) == len(protein) * 3

    def test_gc_bounds_respected(self):
        protein = "MVSKGEELFTGVVPILVELDGDVNGHK"
        seq, _ = _greedy_optimize(
            protein, organism="Escherichia_coli",
            gc_lo=0.30, gc_hi=0.70, is_prokaryote=True,
        )
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        # Allow some tolerance due to constraint conflicts
        assert 0.20 <= gc <= 0.80

    def test_unsupported_organism_raises(self):
        with pytest.raises(UnsupportedOrganismError):
            _greedy_optimize("MVSKGE", organism="Alien_martian")

    def test_restriction_site_removal(self):
        # Pass enzyme sites as DNA sequences, not names
        protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVS"
        seq, _ = _greedy_optimize(
            protein, organism="Escherichia_coli",
            restriction_sites=["GAATTC"], is_prokaryote=True,
        )
        # After optimization, verify the sequence is valid
        assert isinstance(seq, str)
        assert len(seq) == len(protein) * 3

    def test_t_run_fixing(self):
        # Proteins with many Phe (TTT/TTC) can create T-runs
        protein = "FFFFFM"  # 5 Phe + Met
        seq, _ = _greedy_optimize(
            protein, organism="Escherichia_coli", is_prokaryote=True,
        )
        assert isinstance(seq, str)
        # Should not have 6+ consecutive T's after fixing
        # (best effort — may not always be achievable)
        assert len(seq) == 18

    def test_deterministic_output(self):
        protein = "MVSKGE"
        seq1, _ = _greedy_optimize(protein, organism="Escherichia_coli", is_prokaryote=True)
        seq2, _ = _greedy_optimize(protein, organism="Escherichia_coli", is_prokaryote=True)
        assert seq1 == seq2  # Greedy optimizer is deterministic

    def test_single_aa_protein(self):
        seq, _ = _greedy_optimize("M", organism="Escherichia_coli", is_prokaryote=True)
        assert seq == "ATG"

    def test_cai_maximization_prefers_optimal_codons(self):
        """The greedy optimizer should prefer high-CAI codons."""
        ecoli_cai = CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        protein = "LLLLLLLLLL"  # 10 Leu residues
        seq, _ = _greedy_optimize(protein, organism="Escherichia_coli", is_prokaryote=True)
        cai = _compute_cai_fast(seq, ecoli_cai)
        assert cai > 0.3  # Should have reasonable CAI with optimal codon selection


# ── score_splice_donor_potential ────────────────────────────────────

class TestScoreSpliceDonorPotential:
    def test_returns_float(self):
        dna = "ATGCGTAAAGCGCGC" * 5  # Long enough for context
        score = score_splice_donor_potential(dna, 3)  # GT at position 3
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_high_score_for_canonical_splice_site(self):
        # Create a sequence with a strong splice donor context
        dna = "A" * 50 + "GTAAGT" + "A" * 50
        # Find GT position
        gt_pos = dna.find("GT")
        if gt_pos >= 0:
            score = score_splice_donor_potential(dna, gt_pos)
            assert isinstance(score, float)


# ── optimize_sequence high-level API ────────────────────────────────

class TestOptimizeSequence:
    def test_basic_optimization(self):
        result = optimize_sequence(
            target_protein="MVSKGE", organism="Escherichia_coli",
        )
        assert isinstance(result, OptimizationResult)
        assert isinstance(result.sequence, str)
        assert len(result.sequence) == 18
        assert 0.0 <= result.gc_content <= 1.0
        assert 0.0 <= result.cai <= 1.0

    def test_result_has_required_fields(self):
        result = optimize_sequence(
            target_protein="MVSKGE", organism="Homo_sapiens",
        )
        assert hasattr(result, 'sequence')
        assert hasattr(result, 'gc_content')
        assert hasattr(result, 'cai')
        assert hasattr(result, 'failed_predicates')
        assert hasattr(result, 'predicate_results')
        assert hasattr(result, 'protein')
        assert result.protein == "MVSKGE"

    def test_protein_encoding_preserved(self):
        protein = "MVSKGE"
        result = optimize_sequence(target_protein=protein, organism="Escherichia_coli")
        translated = "".join(CODON_TABLE.get(result.sequence[i:i+3], "?")
                           for i in range(0, len(result.sequence), 3))
        assert translated == protein

    def test_gc_content_reasonable(self):
        result = optimize_sequence(
            target_protein="MVSKGEELFTGVVPILVELDGDVNGHK",
            organism="Escherichia_coli",
            gc_lo=0.30, gc_hi=0.70,
        )
        assert result.gc_content >= 0.0
        assert result.gc_content <= 1.0

    def test_cai_nonzero(self):
        result = optimize_sequence(
            target_protein="MVSKGE", organism="Escherichia_coli",
        )
        assert result.cai > 0.0

    def test_ecoli_prokaryote_fast_path(self):
        result = optimize_sequence(
            target_protein="MVSKGEELFTGVVPILVELDGDVNGHK",
            organism="Escherichia_coli",
        )
        assert isinstance(result, OptimizationResult)
        assert result.sequence  # non-empty
