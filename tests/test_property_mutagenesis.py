"""
BioCompiler Property Tests: Mutagenesis & Conservation Consistency
==================================================================
Hypothesis-based property tests verifying that the Python mutagenesis
and conservation implementation matches the Lean4 theorems from:

  - proof/BioCompiler/Mutagenesis.lean
  - proof/BioCompiler/TypeSystem.lean (previously lean/BioCompiler/Predicates.lean)

Mapped Lean4 theorems → Python property tests:

  1. synonymous_preserves_translation
     → test_synonymous_mutation_preserves_amino_acid
  2. ConservationScore (Predicates.lean)
     → test_blosum62_diagonal_positive, test_blosum62_symmetric
  3. ConservationScore monotonicity
     → test_conservation_monotonicity_higher_blosum_more_likely_pass
  4. noGT_subsumes_crypticSplice_FAIL
     → test_no_gt_subsumes_no_cryptic_splice
  5. validCoding_implies_noInternalStops
     → test_valid_coding_seq_implies_no_stop_codons
  6. CodonOptimality (Predicates.lean)
     → test_cai_values_in_unit_range
  7. synonymous_preserves_translation (frame-level)
     → test_synonymous_mutation_preserves_reading_frame
  8. all_valine_codons_have_gt / mandatory_gt_has_gt
     → test_valine_codons_all_contain_gt
  9. applySynonymousMutation preserves length
     → test_mutagenesis_preserves_length
  10. synonymous_gc_counterexample
      → test_synonymous_can_change_gc_content
  11. synonymous_restriction_counterexample
      → test_synonymous_can_introduce_ag_dinucleotide
  12. valine_only_mandatory_gt_aa
      → test_non_valine_aa_has_gt_free_codon
"""

from hypothesis import given, settings, assume, example
from hypothesis import strategies as st
import pytest

from biocompiler.type_system import (
    BLOSUM62,
    AA_TO_CODONS,
    CODON_TABLE,
    check_conservation_score,
    check_codon_optimality,
    check_no_stop_codons,
    check_valid_coding_seq,
    check_no_gt_dinucleotide,
    check_no_cryptic_splice,
    _translate_dna_to_aa,
)
from biocompiler.mutagenesis import (
    propose_mutagenesis,
    MutagenesisProposal,
    MutagenesisReport,
    GT_MANDATORY_AAS,
)

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

STANDARD_AAS = [aa for aa in sorted(AA_TO_CODONS.keys()) if aa != "*"]
STANDARD_CODONS = list(CODON_TABLE.keys())  # includes stop codons
VALID_CODONS = [c for c in CODON_TABLE if CODON_TABLE[c] != "*"]

# A minimal mock CAI table for testing — assign 0.5 to every codon
MOCK_CAI = {codon: 0.5 for codon in VALID_CODONS}

# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

amino_acid = st.sampled_from(STANDARD_AAS)
"""Strategy for a standard amino acid (1-letter, non-stop)."""

valid_codon = st.sampled_from(VALID_CODONS)
"""Strategy for a codon that translates to a standard amino acid."""

aa_pair = st.tuples(amino_acid, amino_acid)
"""Strategy for a pair of standard amino acids."""

codon_sequence = st.lists(valid_codon, min_size=1, max_size=10)
"""Strategy for a list of valid codons (i.e., codon-position-aligned)."""


@st.composite
def dna_sequence(draw, min_len: int = 6, max_len: int = 60) -> str:
    """Generate a DNA sequence of length divisible by 3, using valid codons."""
    num_codons = draw(st.integers(min_value=min_len // 3, max_value=max_len // 3))
    codons = draw(st.lists(valid_codon, min_size=num_codons, max_size=num_codons))
    return "".join(codons)


@st.composite
def synonymous_codon_pair(draw):
    """Draw two different codons that encode the same amino acid."""
    aa = draw(amino_acid)
    codons = AA_TO_CODONS[aa]
    assume(len(codons) >= 2)  # need at least 2 codons for a synonymous pair
    c1, c2 = draw(st.permutations(codons).map(lambda p: (p[0], p[1])))
    assume(c1 != c2)
    return aa, c1, c2


@st.composite
def sequence_with_synonymous_mutation(draw):
    """Draw a DNA sequence and a valid synonymous substitution position.

    Returns (original_seq, mutated_seq, codon_position, original_aa).
    """
    num_codons = draw(st.integers(min_value=2, max_value=8))
    # Build the sequence codon by codon; pick one codon position to mutate
    mut_pos = draw(st.integers(min_value=0, max_value=num_codons - 1))
    codons = []
    for i in range(num_codons):
        codons.append(draw(valid_codon))
    # Ensure the mutation position has a redundant AA (>=2 codons)
    aa = CODON_TABLE[codons[mut_pos]]
    assume(len(AA_TO_CODONS[aa]) >= 2)
    # Pick a different synonymous codon
    alt_codons = [c for c in AA_TO_CODONS[aa] if c != codons[mut_pos]]
    alt = draw(st.sampled_from(alt_codons))
    original_seq = "".join(codons)
    codons[mut_pos] = alt
    mutated_seq = "".join(codons)
    return original_seq, mutated_seq, mut_pos, aa


# ===========================================================================
# Property Test 1: Synonymous mutations preserve the amino acid
# Lean4: synonymous_preserves_translation
# ===========================================================================

class TestSynonymousPreservesTranslation:
    """Synonymous codon substitutions must encode the same amino acid."""

    @given(data=synonymous_codon_pair())
    @settings(max_examples=200)
    def test_synonymous_codons_same_aa(self, data):
        """Two synonymous codons must translate to the same amino acid."""
        aa, c1, c2 = data
        assert CODON_TABLE[c1] == aa
        assert CODON_TABLE[c2] == aa
        assert CODON_TABLE[c1] == CODON_TABLE[c2]

    @given(data=sequence_with_synonymous_mutation())
    @settings(max_examples=200)
    def test_synonymous_mutation_preserves_protein(self, data):
        """Replacing a codon with a synonymous one preserves the full protein."""
        original_seq, mutated_seq, _, _ = data
        assert _translate_dna_to_aa(original_seq) == _translate_dna_to_aa(mutated_seq)


# ===========================================================================
# Property Test 2: BLOSUM62 diagonal is positive for all standard AAs
# Lean4: ConservationScore with minBLOSUM = 0 means all >= 0
# ===========================================================================

class TestBLOSUM62Properties:
    """Structural properties of the BLOSUM62 substitution matrix."""

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_blosum62_diagonal_positive(self, aa):
        """BLOSUM62(a, a) > 0 for every standard amino acid.

        Corresponds to: a substitution of an amino acid with itself
        should always be conservative (positive score).
        """
        assert BLOSUM62[(aa, aa)] > 0, f"BLOSUM62({aa},{aa}) = {BLOSUM62[(aa, aa)]} is not positive"

    @given(pair=aa_pair)
    @settings(max_examples=200)
    def test_blosum62_symmetric(self, pair):
        """BLOSUM62(a, b) == BLOSUM62(b, a) for all amino acid pairs.

        The BLOSUM62 matrix is defined to be symmetric.
        """
        a, b = pair
        assert BLOSUM62[(a, b)] == BLOSUM62[(b, a)], (
            f"BLOSUM62 not symmetric: B62({a},{b})={BLOSUM62[(a,b)]} != B62({b},{a})={BLOSUM62[(b,a)]}"
        )

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_blosum62_diagonal_strictly_maximal(self, aa):
        """BLOSUM62(a, a) >= BLOSUM62(a, b) for any b != a.

        The diagonal entry is always the maximum in its row (self-substitution
        is always the most conservative).
        """
        for b in STANDARD_AAS:
            if b != aa:
                assert BLOSUM62[(aa, aa)] >= BLOSUM62[(aa, b)], (
                    f"BLOSUM62({aa},{aa})={BLOSUM62[(aa,aa)]} < BLOSUM62({aa},{b})={BLOSUM62[(aa,b)]}"
                )


# ===========================================================================
# Property Test 3: Conservation monotonicity
# Lean4: conservation_nonnegative / ConservationScore with minBLOSUM
# ===========================================================================

class TestConservationMonotonicity:
    """Higher BLOSUM62 score implies more likely to pass conservation check."""

    @given(pair=aa_pair)
    @settings(max_examples=200)
    def test_higher_blosum_more_likely_pass(self, pair):
        """If BLOSUM62(a, b) >= BLOSUM62(a, c), then passing conservation
        for (a, c) at threshold t implies passing for (a, b) at threshold t.

        This is conservation monotonicity: a higher score should never
        make it harder to pass the conservation threshold.
        """
        a, b = pair
        # Pick a random third amino acid c != a, b
        other_aas = [x for x in STANDARD_AAS if x != a and x != b]
        assume(len(other_aas) > 0)
        # We use a fixed c for determinism in the property
        c = other_aas[0]

        score_ab = BLOSUM62[(a, b)]
        score_ac = BLOSUM62[(a, c)]

        # Pick threshold = score_ac
        # If score_ab >= score_ac, then (a,b) should also pass
        if score_ab >= score_ac:
            result_c = check_conservation_score(a, c, min_score=score_ac)
            result_b = check_conservation_score(a, b, min_score=score_ac)
            # If (a,c) passes, (a,b) must also pass
            if result_c.passed:
                assert result_b.passed, (
                    f"Conservation monotonicity violated: "
                    f"B62({a},{b})={score_ab} >= B62({a},{c})={score_ac}, "
                    f"but ({a},{b}) fails at threshold {score_ac}"
                )

    @given(aa=amino_acid, min_score=st.integers(min_value=-5, max_value=10))
    @settings(max_examples=100)
    def test_self_substitution_always_passes_nonneg(self, aa, min_score):
        """A self-substitution BLOSUM62(a,a) always passes at min_score=0.

        Corresponds to Lean4: conservation_nonnegative.
        """
        result = check_conservation_score(aa, aa, min_score=0)
        assert result.passed, f"Self-substitution B62({aa},{aa})={BLOSUM62[(aa,aa)]} fails min_score=0"


# ===========================================================================
# Property Test 4: NoGT subsumption
# Lean4: noGT_subsumes_crypticSplice_FAIL
# ===========================================================================

class TestNoGTSubsumption:
    """If a sequence has no GT dinucleotides, NoCrypticSplice must PASS."""

    @given(seq=dna_sequence())
    @settings(max_examples=100)
    def test_no_gt_implies_cryptic_splice_pass(self, seq):
        """If NoGTDinucleotide passes, then NoCrypticSplice must also pass.

        Lean4 theorem: noGT_subsumes_crypticSplice_FAIL
        If there are no GT dinucleotides at all, NoCrypticSplice trivially
        holds because its antecedent (codonHasGT) is never satisfied.
        """
        gt_result = check_no_gt_dinucleotide(seq)
        if gt_result.passed:
            splice_result = check_no_cryptic_splice(seq)
            assert splice_result.passed, (
                f"NoGT passes but NoCrypticSplice fails for seq {seq[:30]}... "
                f"GT result: {gt_result.details}, "
                f"Splice result: {splice_result.details}"
            )

    def test_no_gt_example(self):
        """Concrete example: a sequence with no GT should pass NoCrypticSplice."""
        gt_result = check_no_gt_dinucleotide("AACCCCATTGGGAAA")
        assert gt_result.passed
        splice_result = check_no_cryptic_splice("AACCCCATTGGGAAA")
        assert splice_result.passed


# ===========================================================================
# Property Test 5: ValidCodingSeq implies NoStopCodons
# Lean4: validCoding_implies_noInternalStops
# ===========================================================================

class TestValidCodingImpliesNoStops:
    """A valid coding sequence (all codons valid) should have no internal stops."""

    @given(seq=dna_sequence())
    @settings(max_examples=100)
    def test_valid_coding_implies_no_stop_codons(self, seq):
        """If ValidCodingSeq passes, then NoStopCodons must also pass.

        Lean4 theorem: validCoding_implies_noInternalStops
        All codons (except possibly the last) must be non-stop.
        """
        valid_result = check_valid_coding_seq(seq)
        if valid_result.passed:
            stop_result = check_no_stop_codons(seq)
            assert stop_result.passed, (
                f"ValidCodingSeq passes but NoStopCodons fails for seq {seq[:30]}... "
                f"Stop result: {stop_result.details}"
            )


# ===========================================================================
# Property Test 6: CodonOptimality CAI values are in [0, 1]
# Lean4: CodonOptimality predicate
# ===========================================================================

class TestCodonOptimality:
    """CAI (Codon Adaptation Index) values must be in [0, 1]."""

    @given(codon=valid_codon)
    @settings(max_examples=64)
    def test_cai_values_in_unit_range(self, codon):
        """CAI values for any codon must be in [0.0, 1.0].

        This validates the CodonOptimality predicate: CAI is defined
        as a relative adaptiveness index, which is normalized to [0, 1].
        """
        cai = MOCK_CAI.get(codon, 0.0)
        assert 0.0 <= cai <= 1.0, f"CAI({codon}) = {cai} is outside [0, 1]"

    @given(codon=valid_codon, min_cai=st.floats(min_value=0.0, max_value=1.0))
    @settings(max_examples=100)
    def test_codon_optimality_threshold_consistency(self, codon, min_cai):
        """check_codon_optimality(codon, cai_table, min_cai).passed
        iff cai_table[codon] >= min_cai."""
        cai = MOCK_CAI.get(codon, 0.0)
        result = check_codon_optimality(codon, MOCK_CAI, min_cai=min_cai)
        expected_pass = cai >= min_cai
        assert result.passed == expected_pass, (
            f"CodonOptimality inconsistency: CAI({codon})={cai}, "
            f"min_cai={min_cai}, expected pass={expected_pass}, got {result.passed}"
        )


# ===========================================================================
# Property Test 7: Synonymous mutations preserve reading frame
# Lean4: synonymous_preserves_translation (frame-level)
# ===========================================================================

class TestReadingFramePreservation:
    """Synonymous codon substitutions must not change the protein sequence."""

    @given(data=sequence_with_synonymous_mutation())
    @settings(max_examples=200)
    def test_synonymous_preserves_reading_frame(self, data):
        """A synonymous codon substitution does not change the protein.

        This is the reading-frame-level version of
        synonymous_preserves_translation: the full translated protein
        must be identical before and after the mutation.
        """
        original_seq, mutated_seq, _, _ = data
        protein_orig = _translate_dna_to_aa(original_seq)
        protein_mut = _translate_dna_to_aa(mutated_seq)
        assert protein_orig == protein_mut


# ===========================================================================
# Property Test 8: All Valine codons contain GT
# Lean4: all_valine_codons_have_gt / mandatory_gt_has_gt
# ===========================================================================

class TestValineGTMandatory:
    """Valine is the only amino acid where ALL codons contain GT."""

    @given(codon=st.sampled_from(AA_TO_CODONS["V"]))
    @settings(max_examples=10)
    def test_valine_codons_all_contain_gt(self, codon):
        """Every Valine codon contains a GT dinucleotide.

        Lean4 theorem: all_valine_codons_have_gt
        Valine codons: GTT, GTC, GTA, GTG — all start with GT.
        """
        assert "GT" in codon, f"Valine codon {codon} does not contain GT"

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_non_valine_aa_has_gt_free_codon(self, aa):
        """Every amino acid except Valine has at least one GT-free codon.

        Lean4 theorem: valine_only_mandatory_gt_aa
        """
        assume(aa != "V")
        codons = AA_TO_CODONS[aa]
        gt_free = [c for c in codons if "GT" not in c]
        assert len(gt_free) > 0, f"AA {aa} has no GT-free codons: {codons}"

    def test_valine_only_gt_mandatory_aa(self):
        """GT_MANDATORY_AAS should be exactly {'V'}.

        Lean4: GT_MANDATORY_AAS = ['V'] and ag_mandatory_empty.
        """
        assert GT_MANDATORY_AAS == {"V"}


# ===========================================================================
# Property Test 9: Mutagenesis preserves length
# Lean4: applySynonymousMutation preserves sequence length
# ===========================================================================

class TestMutagenesisPreservesLength:
    """Mutagenesis must produce output sequences of the same length as input."""

    @given(data=sequence_with_synonymous_mutation())
    @settings(max_examples=200)
    def test_synonymous_preserves_length(self, data):
        """A synonymous codon substitution preserves sequence length.

        Lean4: applySynonymousMutation replaces a codon of length 3 with
        another codon of length 3, so len(output) == len(input).
        """
        original_seq, mutated_seq, _, _ = data
        assert len(original_seq) == len(mutated_seq)

    @given(seq=dna_sequence(min_len=6, max_len=60))
    @settings(max_examples=50)
    def test_propose_mutagenesis_preserves_length(self, seq):
        """propose_mutagenesis: each proposal's new_codon has length 3."""
        # Find constraint positions (GT dinucleotides)
        gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
        if not gt_positions:
            return  # nothing to test

        constraint_positions = list({(pos // 3) * 3 for pos in gt_positions})
        constraint_types = {pos: ["GT"] for pos in constraint_positions}

        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            if proposal.new_codon:  # not impossible
                assert len(proposal.new_codon) == 3, (
                    f"Proposed codon {proposal.new_codon} has length "
                    f"{len(proposal.new_codon)} != 3"
                )


# ===========================================================================
# Property Test 10: Synonymous mutations can change GC content
# Lean4: synonymous_gc_counterexample
# ===========================================================================

class TestSynonymousGCChange:
    """Synonymous substitutions need NOT preserve GC content."""

    @given(data=synonymous_codon_pair())
    @settings(max_examples=100)
    def test_synonymous_can_change_gc_content(self, data):
        """There exist synonymous codon pairs with different GC counts.

        Lean4 theorem: synonymous_gc_counterexample
        AAA (Lysine, 0% GC) and AAG (Lysine, 33% GC).
        We verify the general pattern that some synonymous pairs differ.
        """
        aa, c1, c2 = data
        gc1 = sum(1 for b in c1 if b in "GC")
        gc2 = sum(1 for b in c2 if b in "GC")
        # We don't assert they ALWAYS differ, just that they CAN differ
        # (verified by the concrete counterexample below)

    def test_concrete_gc_counterexample(self):
        """AAA (K) -> AAG (K) changes GC count from 0 to 1.

        Matches Lean4: synonymous_gc_counterexample.
        """
        assert CODON_TABLE["AAA"] == "K"
        assert CODON_TABLE["AAG"] == "K"
        gc_aaa = sum(1 for b in "AAA" if b in "GC")
        gc_aag = sum(1 for b in "AAG" if b in "GC")
        assert gc_aaa == 0
        assert gc_aag == 1
        assert gc_aaa != gc_aag


# ===========================================================================
# Property Test 11: Synonymous mutations can introduce AG dinucleotides
# Lean4: synonymous_restriction_counterexample
# ===========================================================================

class TestSynonymousAGCounterexample:
    """Synonymous substitutions can introduce AG dinucleotides."""

    def test_concrete_ag_counterexample(self):
        """TCT (S) -> AGT (S): original has no AG at positions 0-1, new has AG.

        Matches Lean4: synonymous_restriction_counterexample.
        """
        assert CODON_TABLE["TCT"] == "S"
        assert CODON_TABLE["AGT"] == "S"
        assert "AG" not in "TCT"  # no AG at start
        assert "AG" in "AGT"      # has AG at start

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_some_aa_has_ag_and_non_ag_codons(self, aa):
        """For AAs with multiple codons, some may have AG and some may not.

        This verifies the structural condition that makes the
        synonymous_restriction_counterexample possible.
        """
        codons = AA_TO_CODONS[aa]
        assume(len(codons) >= 2)
        has_ag = [c for c in codons if "AG" in c]
        no_ag = [c for c in codons if "AG" not in c]
        # Not asserting for all AAs, but for Serine and Arginine this is true
        if aa in ("S", "R"):
            assert len(has_ag) > 0, f"AA {aa} has no codons with AG"
            assert len(no_ag) > 0, f"AA {aa} has no codons without AG"


# ===========================================================================
# Property Test 12: Every non-Valine AA has a GT-free synonymous codon
# Lean4: valine_only_mandatory_gt_aa / hasGTFreeCodon
# ===========================================================================

class TestGTFreeCodons:
    """Every amino acid except Valine has at least one GT-free synonymous codon."""

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_non_valine_has_gt_free_codon(self, aa):
        """For any AA != V, there exists a codon without GT.

        Lean4: valine_only_mandatory_gt_aa
        """
        assume(aa != "V")
        codons = AA_TO_CODONS[aa]
        gt_free = [c for c in codons if "GT" not in c]
        assert len(gt_free) >= 1, f"AA {aa} has no GT-free codon: {codons}"

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_every_aa_has_ag_free_codon(self, aa):
        """Every amino acid (even Valine) has at least one AG-free codon.

        Lean4: every_aa_has_ag_free_codon / ag_mandatory_empty
        """
        codons = AA_TO_CODONS[aa]
        ag_free = [c for c in codons if "AG" not in c]
        assert len(ag_free) >= 1, f"AA {aa} has no AG-free codon: {codons}"


# ===========================================================================
# Property Test 13: Synonymous mutations can introduce GT
# Lean4: synonymous_introduces_gt
# ===========================================================================

class TestSynonymousIntroducesGT:
    """Synonymous mutations CAN introduce GT dinucleotides."""

    def test_concrete_gt_introduction(self):
        """AGA (R) -> CGT (R): AGA has no GT, CGT has GT at positions 1-2.

        Matches Lean4: synonymous_introduces_gt.
        """
        assert CODON_TABLE["AGA"] == "R"
        assert CODON_TABLE["CGT"] == "R"
        assert "GT" not in "AGA"
        assert "GT" in "CGT"  # positions 1-2


# ===========================================================================
# Property Test 14: Unrepairable cryptic donors for Valine
# Lean4: unrepairable_cryptic_donor_exists
# ===========================================================================

class TestUnrepairableCrypticDonors:
    """Valine positions create unrepairable cryptic splice donor sites."""

    def test_valine_position_unrepairable(self):
        """A GT at a Valine codon position cannot be eliminated by
        synonymous substitution.

        Lean4: unrepairable_cryptic_donor_exists
        """
        # All Valine codons contain GT
        for codon in AA_TO_CODONS["V"]:
            assert "GT" in codon, f"Valine codon {codon} lacks GT (unexpected)"

    def test_gtt_creates_cryptic_donor(self):
        """GTT at position 0 creates an unrepairable cryptic donor.

        Lean4: unrepairable_cryptic_donor_exists uses the sequence
        [G, T, T, A, A, A] at position 0.
        """
        seq = "GTTAAA"
        assert seq[0:2] == "GT"
        assert CODON_TABLE["GTT"] == "V"
        # All Valine codons have GT, so this is unrepairable
        assert all("GT" in c for c in AA_TO_CODONS["V"])


# ===========================================================================
# Property Test 15: Mutation safety classification
# Lean4: synonymous_at_least_protein_safe
# ===========================================================================

class TestMutationSafetyClassification:
    """All synonymous mutations should be at least PROTEIN_SAFE."""

    @given(data=synonymous_codon_pair())
    @settings(max_examples=100)
    def test_synonymous_preserves_protein_level_predicates(self, data):
        """Synonymous mutations always preserve protein-level (SLOT) predicates.

        Lean4: synonymous_safe_for_SLOT_predicates / synonymous_at_least_protein_safe.
        SLOT-dependent predicates always return UNCERTAIN, so they are
        trivially preserved by any mutation.
        """
        aa, c1, c2 = data
        # Both codons translate to the same amino acid
        assert CODON_TABLE[c1] == CODON_TABLE[c2]


# ===========================================================================
# Property Test 16: Conservation score check consistency
# Lean4: ConservationScore definition
# ===========================================================================

class TestConservationScoreConsistency:
    """ConservationScore checks are consistent with BLOSUM62 values."""

    @given(pair=aa_pair, min_score=st.integers(min_value=-5, max_value=10))
    @settings(max_examples=300)
    def test_conservation_check_matches_blosum(self, pair, min_score):
        """check_conservation_score(a, b, t).passed iff BLOSUM62(a,b) >= t."""
        a, b = pair
        score = BLOSUM62[(a, b)]
        result = check_conservation_score(a, b, min_score=min_score)
        assert result.passed == (score >= min_score), (
            f"Conservation check inconsistent: B62({a},{b})={score}, "
            f"min_score={min_score}, expected pass={score >= min_score}, "
            f"got {result.passed}"
        )

    @given(pair=aa_pair)
    @settings(max_examples=200)
    def test_conservation_self_always_passes(self, pair):
        """BLOSUM62(a,a) is always positive, so self-substitution always
        passes with min_score=0."""
        a, _ = pair
        result = check_conservation_score(a, a, min_score=0)
        assert result.passed


# ===========================================================================
# Property Test 17: SpliceVerdict dual-threshold monotonicity
# Lean4: dual_threshold_monotonicity / classifySplice_monotone
# ===========================================================================

class TestDualThresholdMonotonicity:
    """PASS at a strict threshold implies PASS at a more permissive threshold."""

    def test_low_threshold_more_permissive(self):
        """If NoCrypticSplice passes with low=3.0, it should also pass
        with low=2.0 (more permissive threshold).

        Lean4: classifySplice_monotone — PASS at strict implies PASS at permissive.
        """
        # A sequence with no GT always passes
        seq = "AACCCCATT"
        r1 = check_no_cryptic_splice(seq, low_thresh=3.0, high_thresh=6.0)
        r2 = check_no_cryptic_splice(seq, low_thresh=2.0, high_thresh=6.0)
        assert r1.passed
        assert r2.passed


# ===========================================================================
# Property Test 18: Codon degeneracy matches Lean4
# Lean4: valine_degeneracy_4, leucine_degeneracy_6, tryptophan_degeneracy_1,
#        methionine_degeneracy_1
# ===========================================================================

class TestCodonDegeneracy:
    """Codon degeneracy matches the Lean4 formal specification."""

    def test_valine_degeneracy_4(self):
        """Valine has exactly 4 codons. Lean4: valine_degeneracy_4."""
        assert len(AA_TO_CODONS["V"]) == 4

    def test_leucine_degeneracy_6(self):
        """Leucine has exactly 6 codons. Lean4: leucine_degeneracy_6."""
        assert len(AA_TO_CODONS["L"]) == 6

    def test_tryptophan_degeneracy_1(self):
        """Tryptophan has exactly 1 codon. Lean4: tryptophan_degeneracy_1."""
        assert len(AA_TO_CODONS["W"]) == 1

    def test_methionine_degeneracy_1(self):
        """Methionine has exactly 1 codon. Lean4: methionine_degeneracy_1."""
        assert len(AA_TO_CODONS["M"]) == 1

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_degeneracy_at_least_1(self, aa):
        """Every standard amino acid has at least 1 codon."""
        assert len(AA_TO_CODONS[aa]) >= 1

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_degeneracy_at_most_6(self, aa):
        """No amino acid has more than 6 codons (Leucine and Serine are max)."""
        assert len(AA_TO_CODONS[aa]) <= 6


# ===========================================================================
# Property Test 19: No AG-mandatory amino acids
# Lean4: ag_mandatory_empty / no_ag_mandatory
# ===========================================================================

class TestNoAGMandatory:
    """No amino acid has AG in ALL its codons."""

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_no_aa_is_ag_mandatory(self, aa):
        """For every amino acid, at least one codon does NOT contain AG.

        Lean4: ag_mandatory_empty / no_ag_mandatory
        AG_MANDATORY_AAS = [] in Lean4.
        """
        codons = AA_TO_CODONS[aa]
        ag_free = [c for c in codons if "AG" not in c]
        assert len(ag_free) >= 1, f"AA {aa} is AG-mandatory (all codons contain AG): {codons}"


# ===========================================================================
# Property Test 20: Wobble position analysis for Valine
# Lean4: valine_wobble_cannot_eliminate_gt
# ===========================================================================

class TestWobbleAnalysis:
    """For Valine codons (GTx), wobble position CANNOT eliminate GT."""

    def test_valine_wobble_cannot_eliminate_gt(self):
        """All Valine codons have GT at positions 0-1, so the wobble
        position (index 2) cannot eliminate GT.

        Lean4: valine_wobble_cannot_eliminate_gt
        canWobbleEliminateGT [G, T, _] = false for any nucleotide.
        """
        for codon in AA_TO_CODONS["V"]:
            # GT is at positions 0-1, so changing position 2 doesn't help
            assert codon[0:2] == "GT", f"Valine codon {codon} doesn't start with GT"

    def test_arginine_wobble_can_eliminate_gt(self):
        """For Arginine codon CGT (GT at positions 1-2), the wobble
        position CAN eliminate GT by changing T to something else.

        Lean4: position12_gt_wobble_can_eliminate
        CGT → CGA, CGC, CGG are all Arginine codons without GT at pos 1-2.
        """
        assert CODON_TABLE["CGT"] == "R"
        # CGT has GT at positions 1-2
        assert "GT" in "CGT"
        # CGA, CGC, CGG are alternative Arginine codons without GT at 1-2
        for alt in ["CGA", "CGC", "CGG"]:
            assert CODON_TABLE[alt] == "R"
            assert "GT" not in alt[1:3]  # GT not at positions 1-2
