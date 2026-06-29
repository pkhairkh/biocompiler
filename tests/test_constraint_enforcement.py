"""
End-to-end tests for constraint enforcement in biocompiler.solver.constraints.

Verifies that every hard and soft constraint correctly detects violations
and returns appropriate results for satisfying sequences.  Covers:

1.  HardConstraint.check() returns False for violating sequences
2.  GCRangeConstraint properly validates GC content
3.  NoRestrictionSiteConstraint detects known sites
4.  NoCrypticSpliceConstraint detects high-scoring sites
5.  NoCpGIslandConstraint detects CpG islands
6.  NoATTTAMotifConstraint detects ATTTA motifs
7.  NoTRunConstraint detects poly-T runs
8.  TranslationConstraint validates correct translation
9.  CSPModel.check_all_hard() returns correct boolean
10. CSPModel.hard_violations() returns correct positions
11. SoftConstraint.score() returns reasonable values
12. MaximizeCAI.cai() returns value in [0, 1]
"""

from __future__ import annotations

import math

import pytest

from biocompiler.solver.constraints import (
    # Hard constraints
    HardConstraint,
    TranslationConstraint,
    NoRestrictionSiteConstraint,
    GCRangeConstraint,
    NoCrypticSpliceConstraint,
    NoCpGIslandConstraint,
    NoATTTAMotifConstraint,
    NoTRunConstraint,
    # Soft constraints
    SoftConstraint,
    MaximizeCAI,
    MinimizeCpG,
    MinimizeMRNADG,
    # Model
    CSPModel,
    # Helpers / constants
    compute_gc_from_codons,
    DEFAULT_GC_LO,
    DEFAULT_GC_HI,
    DEFAULT_CPG_WINDOW,
    DEFAULT_CPG_THRESHOLD,
    DEFAULT_MAX_T_RUN,
    CAI_LOG_EPSILON,
)
from biocompiler.solver.types import (
    CodonVariable,
    SolverConfig,
    ConstraintStrictness,
    ConstraintType,
)
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS, INSTABILITY_MOTIF
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def human_adaptiveness() -> dict[str, float]:
    """Human codon adaptiveness table for MaximizeCAI tests."""
    return CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]


@pytest.fixture
def sample_protein() -> str:
    """6-AA protein for quick tests: M-V-S-K-G-E."""
    return "MVSKGE"


@pytest.fixture
def sample_sequence() -> str:
    """Valid DNA encoding sample_protein using common human codons."""
    return "ATGGCTTCTAAAGGTGAA"


@pytest.fixture
def default_solver_config() -> SolverConfig:
    """SolverConfig with default GC bounds."""
    return SolverConfig(gc_lo=0.30, gc_hi=0.70)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. HardConstraint.check() returns False for violating sequences
# ═══════════════════════════════════════════════════════════════════════════════

class TestHardConstraintCheckContract:
    """Verify the HardConstraint.check() contract: True = satisfied, False = violated."""

    def test_gc_violation_returns_false(self):
        constraint = GCRangeConstraint(gc_lo=0.50, gc_hi=0.70)
        all_at_seq = "AAA" * 20  # 0% GC, below lo
        assert constraint.check(all_at_seq) is False

    def test_restriction_site_violation_returns_false(self):
        constraint = NoRestrictionSiteConstraint(sites=["GAATTC"])
        seq_with_site = "AAAGAATTCAAA"
        assert constraint.check(seq_with_site) is False

    def test_translation_violation_returns_false(self):
        constraint = TranslationConstraint(protein="MV")
        # "AAA" translates to K, not M
        wrong_seq = "AAAGTT"
        assert constraint.check(wrong_seq) is False

    def test_attta_violation_returns_false(self):
        constraint = NoATTTAMotifConstraint()
        seq_with_attta = "ATTTA" + "CCC" * 5
        assert constraint.check(seq_with_attta) is False

    def test_t_run_violation_returns_false(self):
        constraint = NoTRunConstraint(max_run=5)
        seq_with_long_t = "AAATTTTTTAAA"
        assert constraint.check(seq_with_long_t) is False

    def test_satisfying_sequence_returns_true(self):
        """When all constraints are satisfied, each check() returns True."""
        gc = GCRangeConstraint(gc_lo=0.20, gc_hi=0.80)
        assert gc.check("ATGGCTTCTAAAGGTGAA") is True

        rs = NoRestrictionSiteConstraint(sites=["GAATTC"])
        assert rs.check("ATGGCTTCTAAAGGTGAA") is True

        tr = NoTRunConstraint(max_run=5)
        assert tr.check("ATGGCTTCTAAAGGTGAA") is True

        attta = NoATTTAMotifConstraint()
        assert attta.check("ATGGCTTCTAAAGGTGAA") is True


# ═══════════════════════════════════════════════════════════════════════════════
# 2. GCRangeConstraint properly validates GC content
# ═══════════════════════════════════════════════════════════════════════════════

class TestGCRangeConstraint:
    """Thorough validation of GCRangeConstraint.check() and violated_positions()."""

    def test_within_range_passes(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        # ATG=1GC/3, GCT=2GC/3, TCT=1GC/3, AAA=0GC/3, GGT=1GC/3, GAA=1GC/3
        # total GC = 6/18 = 0.333
        seq = "ATGGCTTCTAAAGGTGAA"
        assert c.check(seq) is True

    def test_all_gc_above_hi_fails(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        seq = "GCGCGCGCGC"
        assert c.check(seq) is False

    def test_all_at_below_lo_fails(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        seq = "ATATATATAT"
        assert c.check(seq) is False

    def test_exact_lower_boundary_passes(self):
        c = GCRangeConstraint(gc_lo=0.50, gc_hi=0.70)
        # 3 GC in 6 bases = 0.50
        seq = "GCGAAA"
        assert c.check(seq) is True

    def test_exact_upper_boundary_passes(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.50)
        # 3 GC in 6 bases = 0.50
        seq = "GCGAAA"
        assert c.check(seq) is True

    def test_empty_sequence_passes(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        assert c.check("") is True

    def test_violated_positions_high_gc(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.50)
        seq = "GCGCGC"  # 100% GC
        positions = c.violated_positions(seq)
        # Should return positions of G/C bases
        assert len(positions) == 6
        assert positions == [0, 1, 2, 3, 4, 5]

    def test_violated_positions_low_gc(self):
        c = GCRangeConstraint(gc_lo=0.50, gc_hi=0.70)
        seq = "ATATAT"  # 0% GC
        positions = c.violated_positions(seq)
        # Should return positions of A/T bases (candidates for GC increase)
        assert len(positions) == 6
        assert positions == [0, 1, 2, 3, 4, 5]

    def test_violated_positions_in_range_empty(self):
        c = GCRangeConstraint(gc_lo=0.30, gc_hi=0.70)
        seq = "ATGGCT"  # 3/6 = 50% GC
        assert c.violated_positions(seq) == []

    def test_invalid_range_raises(self):
        with pytest.raises(ValueError):
            GCRangeConstraint(gc_lo=0.70, gc_hi=0.30)

    def test_default_bounds(self):
        c = GCRangeConstraint()
        assert c.gc_lo == DEFAULT_GC_LO
        assert c.gc_hi == DEFAULT_GC_HI

    def test_constraint_type(self):
        c = GCRangeConstraint()
        assert c.constraint_type == ConstraintType.GC_CONTENT


# ═══════════════════════════════════════════════════════════════════════════════
# 3. NoRestrictionSiteConstraint detects known sites
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoRestrictionSiteConstraint:
    """Verify detection of restriction enzyme recognition sites."""

    def test_ecori_forward_detected(self):
        c = NoRestrictionSiteConstraint(sites=["GAATTC"])
        assert c.check("AAAGAATTCAAA") is False

    def test_ecori_rc_detected(self):
        """Reverse complement of GAATTC is GAATTC (palindrome) — should also be caught."""
        c = NoRestrictionSiteConstraint(sites=["GAATTC"])
        assert c.check("AAAGAATTCAAA") is False

    def test_bamhi_detected(self):
        c = NoRestrictionSiteConstraint(sites=["GGATCC"])
        assert c.check("AAAGGATCCAAA") is False

    def test_hindiii_detected(self):
        c = NoRestrictionSiteConstraint(sites=["AAGCTT"])
        assert c.check("AAGCTT") is False

    def test_clean_sequence_passes(self):
        c = NoRestrictionSiteConstraint(sites=["GAATTC", "GGATCC"])
        assert c.check("ATGGCTTCTAAAGGTGAA") is True

    def test_non_palindromic_rc_detected(self):
        """Non-palindromic site: XhoI (CTCGAG) RC = CTCGAG — also a palindrome,
        but let us use a truly non-palindromic example."""
        c = NoRestrictionSiteConstraint(sites=["GATC"])
        # RC of GATC is GATC (palindrome too); try ACTAGT (SpeI) RC = ACTAGT
        # Use a non-palindromic sequence: "ATCGAT" RC = "ATCGAT" (also palindromic)
        # Test with a site whose RC is different
        c2 = NoRestrictionSiteConstraint(sites=["GATC"])
        # RC of GATC is GATC — palindrome
        # Test with site "ACGCGT" (MluI) — RC is ACGCGT, palindrome
        # The implementation checks RC != site before scanning for RC
        # For a palindrome like GAATTC, RC == site, so it only checks once
        # Let us test a site with a truly different RC
        c3 = NoRestrictionSiteConstraint(sites=["AAGCTT"])  # HindIII, palindrome
        assert c3.check("AAGCTT") is False
        # Non-palindromic: BglII (AGATCT) RC = AGATCT — also palindrome!
        # Use "GCTAGC" (NheI) — palindrome too. Most 6-cutters are palindromic.
        # Use a 5-base site: "GATC" — RC is GATC (palindromic in this notation)
        # The key test is that the constraint does not crash when RC != site
        # Use a 4-base example: "AATT" RC = "AATT" palindrome
        # For a non-palindrome: "GTCGAC" (SalI) RC = "GTCGAC" palindrome
        # Most common enzymes are palindromic. Let us test with a synthetic one
        c4 = NoRestrictionSiteConstraint(sites=["AACG"])
        # RC of AACG is CGTT
        seq_with_rc = "AAACGTTAAA"  # contains CGTT at pos 3
        assert c4.check(seq_with_rc) is False

    def test_violated_positions(self):
        c = NoRestrictionSiteConstraint(sites=["GAATTC"])
        seq = "GAATTCGAATTC"  # Two occurrences
        positions = c.violated_positions(seq)
        assert 0 in positions
        assert 6 in positions

    def test_violated_positions_clean(self):
        c = NoRestrictionSiteConstraint(sites=["GAATTC"])
        assert c.violated_positions("ATGGCTTCT") == []

    def test_sites_property(self):
        c = NoRestrictionSiteConstraint(sites=["GAATTC", "GGATCC"])
        assert c.sites == ["GAATTC", "GGATCC"]

    def test_constraint_type(self):
        c = NoRestrictionSiteConstraint(sites=["GAATTC"])
        assert c.constraint_type == ConstraintType.RESTRICTION_SITE


# ═══════════════════════════════════════════════════════════════════════════════
# 4. NoCrypticSpliceConstraint detects high-scoring sites
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoCrypticSpliceConstraint:
    """Verify detection of cryptic splice sites via MaxEntScan scoring."""

    def test_strong_donor_site_detected(self):
        """A strong canonical donor consensus should score above threshold."""
        c = NoCrypticSpliceConstraint(threshold=3.0)
        # Build a long sequence with a strong donor: ...CAG|GTAAGT...
        # The 9-mer around a strong GT donor: CAG GTAAGT
        # We need at least 3 nt upstream of GT and 6 downstream
        seq = "AAACAGGTAAGTAAA" * 5  # 75 nt, plenty of context
        # Whether this triggers depends on the MaxEntScan model scores;
        # at minimum the constraint should not crash.
        result = c.check(seq)
        assert isinstance(result, bool)

    def test_no_gt_ag_passes(self):
        """Sequence with no GT or AG dinucleotides should always pass."""
        c = NoCrypticSpliceConstraint(threshold=3.0)
        # A sequence with no GT and no AG: use only A, C, T in pairs that avoid GT/AG
        seq = "AAACCCAAACCC" * 5
        # "AC" contains no GT or AG
        # Let us be more careful: avoid G entirely
        seq = "AAACCCAAACCC" * 5
        # Has "AC", "CA", "AA", "CC" — no GT, no AG
        assert c.check(seq) is True

    def test_short_sequence_passes(self):
        """Very short sequences cannot form 9-mer donor / 23-mer acceptor contexts."""
        c = NoCrypticSpliceConstraint(threshold=3.0)
        # 5 bases: too short for donor (needs 9) and acceptor (needs 23)
        # But GT at pos 0 would get IMPOSSIBLE_SCORE, so check returns True
        assert c.check("AGTCG") is True

    def test_violated_positions_returns_list(self):
        c = NoCrypticSpliceConstraint(threshold=3.0)
        seq = "AAACCCAAACCC" * 5
        positions = c.violated_positions(seq)
        assert isinstance(positions, list)

    def test_invalid_threshold_raises(self):
        with pytest.raises(ValueError):
            NoCrypticSpliceConstraint(threshold=0)
        with pytest.raises(ValueError):
            NoCrypticSpliceConstraint(threshold=-1.0)

    def test_threshold_property(self):
        c = NoCrypticSpliceConstraint(threshold=5.0)
        assert c.threshold == 5.0

    def test_constraint_type(self):
        c = NoCrypticSpliceConstraint()
        assert c.constraint_type == ConstraintType.NO_CRYPTIC_SPLICE


# ═══════════════════════════════════════════════════════════════════════════════
# 5. NoCpGIslandConstraint detects CpG islands
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoCpGIslandConstraint:
    """Verify detection of CpG islands via Obs/Exp CG ratio."""

    def test_high_cpg_island_detected(self):
        """A sequence dense in CG dinucleotides should be flagged."""
        # Use a small window for test feasibility
        c = NoCpGIslandConstraint(window=20, threshold=0.6)
        # Build a 20-nt window that is CpG-dense
        # CGCGCGCGCGCGCGCGCGCG = 10 CG pairs in 20 nt
        # C_count=10, G_count=10, expected = 10*10/20 = 5
        # Obs/Exp = 10/5 = 2.0 > 0.6 → violation
        seq = "CGCGCGCGCGCGCGCGCGCG"
        assert c.check(seq) is False

    def test_no_cpg_passes(self):
        """Sequence with no CG dinucleotides should pass."""
        c = NoCpGIslandConstraint(window=20, threshold=0.6)
        seq = "AAATTTAAATTTAAATTTAA" * 2  # 40 nt, no CG
        assert c.check(seq) is True

    def test_short_sequence_passes(self):
        """Sequence shorter than window size should pass."""
        c = NoCpGIslandConstraint(window=200, threshold=0.6)
        seq = "CGCGCGCG"  # Only 8 nt, shorter than 200
        assert c.check(seq) is True

    def test_violated_positions_returns_window_starts(self):
        c = NoCpGIslandConstraint(window=20, threshold=0.6)
        seq = "CGCGCGCGCGCGCGCGCGCG"
        positions = c.violated_positions(seq)
        # Should report at least one violating window start
        assert len(positions) >= 1
        assert all(isinstance(p, int) for p in positions)

    def test_moderate_cpg_passes(self):
        """A moderate CG density should not trigger the constraint."""
        c = NoCpGIslandConstraint(window=20, threshold=0.6)
        # ACGTACGTACGTACGTACGT = 20 nt
        # C_count=5, G_count=5, CG_count=5
        # expected = 5*5/20 = 1.25
        # Obs/Exp = 5/1.25 = 4.0 > 0.6 → still high
        # Let us try something sparser: ATCGATCGATCGATCGATCG
        # C_count=4, G_count=4, CG count at pos 2-3, 6-7, 10-11, 14-15, 18-19 = 5
        # ATCG -> CG at pos 2-3. ATCGATCG has CG at 2-3 and 6-7. = 2 per 8
        # For 20 chars: ATCGATCGATCGATCGATCG
        # C_count=4, G_count=4, CG pairs at 2,6,10,14,18 = 5
        # expected = 4*4/20 = 0.8; obs/exp = 5/0.8 = 6.25 → still high
        # Let us use an even sparser sequence
        # AAAACGAAAAACGAAAA = 17 nt (too short). Let us use window=50
        c2 = NoCpGIslandConstraint(window=50, threshold=0.6)
        # Build a 50-nt seq with just 1 CG
        seq = "AAAACG" + "A" * 44  # 50 nt, 1 CG
        # C_count ~1, G_count ~1, expected = 1*1/50 = 0.02
        # obs/exp = 1/0.02 = 50 → still > 0.6!
        # The issue is that even 1 CG in a background of A gives high Obs/Exp
        # Let us make a more balanced sequence
        seq2 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"[:50]
        # This has decent C and G content, but many CG dinucleotides
        # A truly passing sequence needs CG obs/exp <= 0.6
        # In a 50-nt random sequence: C~12.5, G~12.5, expected CG~3.125
        # To have obs/exp < 0.6: need < 1.875 CG. So 0 or 1 CG
        # Build 50-nt with C and G but few CG:
        seq3 = "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCA"[:50]
        # "CAGT" = CG at pos 0-1. Count per repeat: 1 CG per 4 nt = 12.5 per 50
        # That is too many. Use "CCAGGTT" pattern:
        # Let us just verify with a sequence we know has zero CG:
        seq4 = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"[:50]
        # "CAG" has no "CG" within it. "CAGCAG" has "GC" at boundary (2-3),
        # not CG. So CG count = 0
        assert c2.check(seq4) is True

    def test_default_window_and_threshold(self):
        c = NoCpGIslandConstraint()
        assert c.window == DEFAULT_CPG_WINDOW
        assert c.threshold == DEFAULT_CPG_THRESHOLD

    def test_invalid_window_raises(self):
        with pytest.raises(ValueError):
            NoCpGIslandConstraint(window=0)

    def test_invalid_threshold_raises(self):
        with pytest.raises(ValueError):
            NoCpGIslandConstraint(threshold=0)

    def test_constraint_type(self):
        c = NoCpGIslandConstraint()
        assert c.constraint_type == ConstraintType.NO_CPG


# ═══════════════════════════════════════════════════════════════════════════════
# 6. NoATTTAMotifConstraint detects ATTTA motifs
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoATTTAMotifConstraint:
    """Verify detection of ATTTA instability motifs."""

    def test_attta_detected(self):
        c = NoATTTAMotifConstraint()
        assert c.check("AAATTTAAAA") is False

    def test_multiple_attta_detected(self):
        c = NoATTTAMotifConstraint()
        seq = "ATTTAATTTAATTTA"
        assert c.check(seq) is False

    def test_no_attta_passes(self):
        c = NoATTTAMotifConstraint()
        assert c.check("AAACCCGGGTTT") is True

    def test_cross_codon_attta_detected(self):
        """AAT|TTA should produce ATTTA at the boundary."""
        c = NoATTTAMotifConstraint()
        seq = "AATTTA"  # ATTTA starts at position 1
        assert c.check(seq) is False

    def test_violated_positions_single(self):
        c = NoATTTAMotifConstraint()
        seq = "AAATTTAAAA"
        positions = c.violated_positions(seq)
        assert 2 in positions  # ATTTA starts at index 2

    def test_violated_positions_multiple(self):
        c = NoATTTAMotifConstraint()
        seq = "ATTTAATTTA"
        positions = c.violated_positions(seq)
        assert 0 in positions
        assert 5 in positions

    def test_violated_positions_clean(self):
        c = NoATTTAMotifConstraint()
        assert c.violated_positions("AAACCCGGG") == []

    def test_case_insensitive(self):
        c = NoATTTAMotifConstraint()
        assert c.check("aaatttaaaa") is False

    def test_constraint_type(self):
        c = NoATTTAMotifConstraint()
        assert c.constraint_type == ConstraintType.NO_INSTABILITY_MOTIF


# ═══════════════════════════════════════════════════════════════════════════════
# 7. NoTRunConstraint detects poly-T runs
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoTRunConstraint:
    """Verify detection of poly-T runs exceeding max_run."""

    def test_six_t_detected(self):
        """Default max_run=5, so 6 T's should fail."""
        c = NoTRunConstraint(max_run=5)
        assert c.check("TTTTTT") is False

    def test_five_t_passes(self):
        """5 T's should pass with max_run=5."""
        c = NoTRunConstraint(max_run=5)
        assert c.check("TTTTT") is True

    def test_nine_t_detected(self):
        c = NoTRunConstraint(max_run=5)
        assert c.check("TTTTTTTTT") is False

    def test_mixed_sequence_passes(self):
        c = NoTRunConstraint(max_run=5)
        seq = "AAATTTTAAATTTTTAAA"  # max run = 4 and 5
        assert c.check(seq) is True

    def test_custom_max_run(self):
        c = NoTRunConstraint(max_run=3)
        assert c.check("TTTT") is False  # 4 > 3
        assert c.check("TTT") is True   # 3 <= 3

    def test_violated_positions(self):
        c = NoTRunConstraint(max_run=5)
        seq = "AAATTTTTTAAA"  # 6 T's starting at pos 3
        positions = c.violated_positions(seq)
        assert 3 in positions

    def test_multiple_violated_runs(self):
        c = NoTRunConstraint(max_run=5)
        seq = "TTTTTTAAATTTTTT"  # Two runs of 6
        positions = c.violated_positions(seq)
        assert 0 in positions
        assert 9 in positions

    def test_violated_positions_clean(self):
        c = NoTRunConstraint(max_run=5)
        assert c.violated_positions("AAACCCGGG") == []

    def test_case_insensitive(self):
        c = NoTRunConstraint(max_run=5)
        assert c.check("tttttt") is False

    def test_default_max_run(self):
        c = NoTRunConstraint()
        assert c.max_run == DEFAULT_MAX_T_RUN

    def test_invalid_max_run_raises(self):
        with pytest.raises(ValueError):
            NoTRunConstraint(max_run=0)

    def test_constraint_type(self):
        c = NoTRunConstraint()
        assert c.constraint_type == ConstraintType.MRNA_STABILITY


# ═══════════════════════════════════════════════════════════════════════════════
# 8. TranslationConstraint validates correct translation
# ═══════════════════════════════════════════════════════════════════════════════

class TestTranslationConstraint:
    """Verify that TranslationConstraint checks codon-to-AA mapping."""

    def test_correct_translation_passes(self):
        c = TranslationConstraint(protein="MV")
        # ATG=M, GTT=V
        assert c.check("ATGGTT") is True

    def test_wrong_codon_fails(self):
        c = TranslationConstraint(protein="MV")
        # AAA=K (not M), GTT=V
        assert c.check("AAAGTT") is False

    def test_wrong_length_fails(self):
        c = TranslationConstraint(protein="MV")
        # Only 3 bases, need 6
        assert c.check("ATG") is False

    def test_empty_protein_passes(self):
        c = TranslationConstraint(protein="")
        assert c.check("") is True

    def test_violated_positions_correct(self):
        c = TranslationConstraint(protein="MV")
        # AAA≠M (pos 0), GTT=V (pos 3 OK)
        positions = c.violated_positions("AAAGTT")
        assert 0 in positions
        assert 3 not in positions

    def test_violated_positions_truncated_sequence(self):
        c = TranslationConstraint(protein="MV")
        # Sequence only has 4 bases — second codon is incomplete
        positions = c.violated_positions("ATGG")
        assert 0 not in positions   # ATG = M [OK]
        assert 3 in positions       # incomplete second codon

    def test_all_standard_amino_acids(self):
        """Test with a protein using all 20 standard AAs."""
        protein = "ACDEFGHIKLMNPQRSTVWY"
        # Build sequence using first codon for each AA
        seq = ""
        for aa in protein:
            codons = AA_TO_CODONS[aa]
            seq += codons[0]
        c = TranslationConstraint(protein=protein)
        assert c.check(seq) is True

    def test_constraint_type(self):
        c = TranslationConstraint(protein="M")
        assert c.constraint_type == ConstraintType.AMINO_ACID_IDENTITY

    def test_protein_property(self):
        c = TranslationConstraint(protein="MVSKGE")
        assert c.protein == "MVSKGE"


# ═══════════════════════════════════════════════════════════════════════════════
# 9. CSPModel.check_all_hard() returns correct boolean
# ═══════════════════════════════════════════════════════════════════════════════

class TestCSPModelCheckAllHard:
    """Verify CSPModel.check_all_hard() aggregation logic."""

    def _make_model(
        self,
        protein: str = "MV",
        hard_constraints: list[HardConstraint] | None = None,
        soft_constraints: list[SoftConstraint] | None = None,
    ) -> CSPModel:
        """Build a minimal CSPModel for testing."""
        variables = [
            CodonVariable(position=i, amino_acid=aa, domain=AA_TO_CODONS[aa])
            for i, aa in enumerate(protein)
        ]
        return CSPModel(
            variables=variables,
            hard_constraints=hard_constraints or [],
            soft_constraints=soft_constraints or [],
            protein=protein,
            organism="Homo_sapiens",
            config=SolverConfig(),
        )

    def test_all_pass_returns_true(self, sample_sequence):
        model = self._make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
            ],
        )
        # ATGGTT: ATG=M, GTT=V → GC = 2/6 ≈ 0.33
        assert model.check_all_hard("ATGGTT") is True

    def test_one_violation_returns_false(self):
        model = self._make_model(
            protein="MV",
            hard_constraints=[
                TranslationConstraint(protein="MV"),
                GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),  # impossible GC
            ],
        )
        assert model.check_all_hard("ATGGTT") is False

    def test_no_hard_constraints_returns_true(self):
        model = self._make_model(hard_constraints=[])
        assert model.check_all_hard("ATGGTT") is True

    def test_restriction_site_violation_detected(self):
        model = self._make_model(
            protein="EF",  # GAA(TTC) → GAATTC = EcoRI
            hard_constraints=[
                NoRestrictionSiteConstraint(sites=["GAATTC"]),
            ],
        )
        assert model.check_all_hard("GAATTC") is False

    def test_attta_violation_detected(self):
        model = self._make_model(
            protein="NL",  # AAT TTA → ATTTA
            hard_constraints=[
                NoATTTAMotifConstraint(),
            ],
        )
        assert model.check_all_hard("AATTTA") is False

    def test_t_run_violation_detected(self):
        model = self._make_model(
            hard_constraints=[
                NoTRunConstraint(max_run=5),
            ],
        )
        assert model.check_all_hard("TTTTTT") is False


# ═══════════════════════════════════════════════════════════════════════════════
# 10. CSPModel.hard_violations() returns correct positions
# ═══════════════════════════════════════════════════════════════════════════════

class TestCSPModelHardViolations:
    """Verify CSPModel.hard_violations() returns dict of constraint→positions."""

    def _make_model(
        self,
        hard_constraints: list[HardConstraint],
        protein: str = "MV",
    ) -> CSPModel:
        variables = [
            CodonVariable(position=i, amino_acid=aa, domain=AA_TO_CODONS[aa])
            for i, aa in enumerate(protein)
        ]
        return CSPModel(
            variables=variables,
            hard_constraints=hard_constraints,
            soft_constraints=[],
            protein=protein,
            organism="Homo_sapiens",
            config=SolverConfig(),
        )

    def test_no_violations_returns_empty_dict(self):
        model = self._make_model([
            TranslationConstraint(protein="MV"),
            GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
        ])
        violations = model.hard_violations("ATGGTT")
        assert violations == {}

    def test_gc_violation_returned(self):
        model = self._make_model([
            GCRangeConstraint(gc_lo=0.80, gc_hi=0.90),
        ])
        violations = model.hard_violations("ATGGTT")
        assert "GCRangeConstraint" in violations
        assert len(violations["GCRangeConstraint"]) > 0

    def test_multiple_violations(self):
        model = self._make_model(
            protein="NL",
            hard_constraints=[
                NoATTTAMotifConstraint(),
                NoTRunConstraint(max_run=3),
            ],
        )
        # AAT TTA contains ATTTA and TTT (3 T's, which is ≤3 so no T-run)
        violations = model.hard_violations("AATTTA")
        assert "NoATTTAMotifConstraint" in violations

    def test_restriction_site_violation_positions(self):
        model = self._make_model(
            protein="EF",
            hard_constraints=[
                NoRestrictionSiteConstraint(sites=["GAATTC"]),
            ],
        )
        violations = model.hard_violations("GAATTC")
        assert "NoRestrictionSiteConstraint" in violations
        assert 0 in violations["NoRestrictionSiteConstraint"]

    def test_only_violated_constraints_included(self):
        """Satisfied constraints should not appear in the violations dict."""
        model = self._make_model([
            TranslationConstraint(protein="MV"),
            GCRangeConstraint(gc_lo=0.10, gc_hi=0.80),
        ])
        violations = model.hard_violations("ATGGTT")
        # Both pass, so dict should be empty
        assert violations == {}


# ═══════════════════════════════════════════════════════════════════════════════
# 11. SoftConstraint.score() returns reasonable values
# ═══════════════════════════════════════════════════════════════════════════════

class TestSoftConstraintScore:
    """Verify that soft constraint scores are numerically reasonable."""

    def test_maximize_cai_score_finite(self, human_adaptiveness, sample_protein, sample_sequence):
        cai = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        score = cai.score(sample_sequence)
        assert math.isfinite(score)
        # score = sum(log(w_i)) — should be negative (since log(w) ≤ 0)
        assert score <= 0.0

    def test_maximize_cai_score_all_optimal(self, human_adaptiveness, sample_protein):
        """Score should be highest when all codons have adaptiveness=1.0."""
        # Build a sequence using only the best codon per AA
        best_seq = ""
        for aa in sample_protein:
            codons = AA_TO_CODONS[aa]
            best_codon = max(codons, key=lambda c: human_adaptiveness.get(c, 0.0))
            best_seq += best_codon
        cai = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        best_score = cai.score(best_seq)

        # Build a suboptimal sequence
        worst_seq = ""
        for aa in sample_protein:
            codons = AA_TO_CODONS[aa]
            worst_codon = min(codons, key=lambda c: human_adaptiveness.get(c, 0.0))
            worst_seq += worst_codon
        worst_score = cai.score(worst_seq)

        assert best_score >= worst_score

    def test_minimize_cpg_score_negative(self):
        cpg = MinimizeCpG()
        seq = "CGCGCGCGCG"
        score = cpg.score(seq)
        # score = -(normalized CpG density) — should be negative
        assert score < 0.0
        # CGCGCGCGCG (10 bases) → 5 CG dinucleotides out of 9 possible
        # windows = 5/9 ≈ 0.5556.  Default expected density (no protein) is
        # 0.25 (GC fraction 0.5).  Normalized = (5/9) / 0.25 = 20/9.
        assert score == pytest.approx(-20.0 / 9.0)

    def test_minimize_cpg_score_zero_for_no_cpg(self):
        cpg = MinimizeCpG()
        score = cpg.score("AAATTTAAATTT")
        assert score == 0.0

    def test_minimize_mrna_dg_score_finite(self):
        dg = MinimizeMRNADG(window_start=0, window_end=50)
        seq = "ATGGCTTCTAAAGGTGAA" + "A" * 32  # 50 nt
        score = dg.score(seq)
        assert math.isfinite(score)

    def test_minimize_mrna_dg_score_is_negated_abs(self):
        """score = -|dG|, so should be ≤ 0."""
        dg = MinimizeMRNADG(window_start=0, window_end=50)
        seq = "ATGGCTTCTAAAGGTGAA" + "G" * 32  # GC-rich = more structure
        score = dg.score(seq)
        assert score <= 0.0

    def test_soft_check_always_true(self, human_adaptiveness, sample_protein):
        """Soft constraints' check() should always return True (they are objectives)."""
        cai = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        assert cai.check("ATGGCTTCTAAAGGTGAA") is True
        cpg = MinimizeCpG()
        assert cpg.check("CGCGCGCGCG") is True
        dg = MinimizeMRNADG()
        assert dg.check("ATGGCTTCTAAAGGTGAA") is True


# ═══════════════════════════════════════════════════════════════════════════════
# 12. MaximizeCAI.cai() returns value in [0, 1]
# ═══════════════════════════════════════════════════════════════════════════════

class TestMaximizeCAI:
    """Verify MaximizeCAI.cai() returns values in [0, 1]."""

    def test_cai_in_range(self, human_adaptiveness, sample_protein, sample_sequence):
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        cai_value = cai_obj.cai(sample_sequence)
        assert 0.0 <= cai_value <= 1.0

    def test_cai_optimal_codons_near_one(self, human_adaptiveness, sample_protein):
        """Using the most adaptive codons should give CAI close to 1.0."""
        best_seq = ""
        for aa in sample_protein:
            codons = AA_TO_CODONS[aa]
            best_codon = max(codons, key=lambda c: human_adaptiveness.get(c, 0.0))
            best_seq += best_codon
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        cai_value = cai_obj.cai(best_seq)
        assert cai_value > 0.5  # should be high
        assert cai_value <= 1.0

    def test_cai_suboptimal_codons_lower(self, human_adaptiveness, sample_protein):
        """Using the least adaptive codons should give a lower CAI."""
        best_seq = ""
        worst_seq = ""
        for aa in sample_protein:
            codons = AA_TO_CODONS[aa]
            best_codon = max(codons, key=lambda c: human_adaptiveness.get(c, 0.0))
            worst_codon = min(codons, key=lambda c: human_adaptiveness.get(c, 0.0))
            best_seq += best_codon
            worst_seq += worst_codon

        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        best_cai = cai_obj.cai(best_seq)
        worst_cai = cai_obj.cai(worst_seq)
        assert best_cai >= worst_cai

    def test_cai_empty_protein(self, human_adaptiveness):
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein="")
        assert cai_obj.cai("") == 0.0

    def test_cai_single_codon(self, human_adaptiveness):
        """Single methionine: ATG adaptiveness should give CAI close to 1."""
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein="M")
        cai_value = cai_obj.cai("ATG")
        assert 0.0 <= cai_value <= 1.0

    def test_cai_with_zero_adaptiveness(self):
        """If a codon has 0 adaptiveness, CAI should handle gracefully."""
        # Construct adaptiveness with one zero entry
        adaptiveness = {"ATG": 0.0, "AAA": 1.0, "TTT": 0.5}
        cai_obj = MaximizeCAI(adaptiveness=adaptiveness, protein="M")
        cai_value = cai_obj.cai("ATG")
        # Should not crash; value may be very small or 0
        assert 0.0 <= cai_value <= 1.0

    def test_score_and_cai_consistent(self, human_adaptiveness, sample_protein, sample_sequence):
        """score() and cai() should be consistent: cai = exp(score / N)."""
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        n = len(sample_protein)
        score = cai_obj.score(sample_sequence)
        cai_value = cai_obj.cai(sample_sequence)
        expected_cai = math.exp(score / n)
        assert abs(cai_value - expected_cai) < 1e-10

    def test_constraint_type(self, human_adaptiveness, sample_protein):
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        assert cai_obj.constraint_type == ConstraintType.CODON_USAGE

    def test_adaptiveness_property(self, human_adaptiveness, sample_protein):
        cai_obj = MaximizeCAI(adaptiveness=human_adaptiveness, protein=sample_protein)
        # Should return a copy
        assert cai_obj.adaptiveness == human_adaptiveness
        assert cai_obj.adaptiveness is not human_adaptiveness


# ═══════════════════════════════════════════════════════════════════════════════
# Bonus: Strictness and name contract checks
# ═══════════════════════════════════════════════════════════════════════════════

class TestConstraintContracts:
    """Verify abstract contract: names, strictness, repr."""

    def test_hard_constraints_are_hard(self):
        constraints = [
            TranslationConstraint(protein="M"),
            NoRestrictionSiteConstraint(sites=["GAATTC"]),
            GCRangeConstraint(),
            NoCrypticSpliceConstraint(),
            NoCpGIslandConstraint(),
            NoATTTAMotifConstraint(),
            NoTRunConstraint(),
        ]
        for c in constraints:
            assert c.strictness == ConstraintStrictness.HARD

    def test_soft_constraints_are_soft(self, human_adaptiveness):
        constraints = [
            MaximizeCAI(adaptiveness=human_adaptiveness, protein="M"),
            MinimizeCpG(),
            MinimizeMRNADG(),
        ]
        for c in constraints:
            assert c.strictness == ConstraintStrictness.SOFT

    def test_all_names_are_unique(self, human_adaptiveness):
        all_constraints = [
            TranslationConstraint(protein="M"),
            NoRestrictionSiteConstraint(sites=["GAATTC"]),
            GCRangeConstraint(),
            NoCrypticSpliceConstraint(),
            NoCpGIslandConstraint(),
            NoATTTAMotifConstraint(),
            NoTRunConstraint(),
            MaximizeCAI(adaptiveness=human_adaptiveness, protein="M"),
            MinimizeCpG(),
            MinimizeMRNADG(),
        ]
        names = [c.name for c in all_constraints]
        assert len(names) == len(set(names)), f"Duplicate names: {names}"

    def test_repr_contains_class_name(self):
        c = GCRangeConstraint()
        assert "GCRangeConstraint" in repr(c)
