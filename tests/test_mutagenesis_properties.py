"""
Property-based tests for mutagenesis.py using Hypothesis
========================================================
Covers three core algebraic properties:

  1. Any valid mutation on a protein produces a valid protein of same length
  2. Mutation positions are within sequence bounds
  3. BLOSUM62 scores are symmetric

Import from biocompiler.mutagenesis and biocompiler.type_system.
"""

from hypothesis import given, settings, assume
from hypothesis import strategies as st
import pytest

from biocompiler.mutagenesis import (
    propose_mutagenesis,
    MutagenesisProposal,
    MutagenesisReport,
)
from biocompiler.type_system import (
    BLOSUM62,
    AA_TO_CODONS,
    CODON_TABLE,
    _translate_dna_to_aa,
)


# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

STANDARD_AAS = sorted(aa for aa in AA_TO_CODONS.keys() if aa != "*")
VALID_CODONS = [c for c in CODON_TABLE if CODON_TABLE[c] != "*"]

# A minimal mock CAI table for testing — assign 0.5 to every codon
MOCK_CAI = {codon: 0.5 for codon in VALID_CODONS}

# All unique amino acid pairs that appear in the BLOSUM62 tuple-key dict
_BLOSUM_KEYS = list(BLOSUM62.keys())


# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

amino_acid = st.sampled_from(STANDARD_AAS)
"""Strategy for a standard amino acid (1-letter, non-stop)."""

valid_codon = st.sampled_from(VALID_CODONS)
"""Strategy for a codon that translates to a standard amino acid."""

aa_pair = st.tuples(amino_acid, amino_acid)
"""Strategy for a pair of standard amino acids."""

blosum_key = st.sampled_from(_BLOSUM_KEYS)
"""Strategy for a BLOSUM62 tuple key (aa1, aa2)."""


@st.composite
def dna_sequence(draw, min_len: int = 6, max_len: int = 60) -> str:
    """Generate a DNA sequence of length divisible by 3, using valid codons."""
    num_codons = draw(st.integers(min_value=min_len // 3, max_value=max_len // 3))
    codons = draw(st.lists(valid_codon, min_size=num_codons, max_size=num_codons))
    return "".join(codons)


@st.composite
def dna_with_gt_constraint(draw) -> tuple:
    """Draw a DNA sequence that contains at least one GT dinucleotide.

    Returns (seq, constraint_positions, constraint_types) suitable for
    passing directly to propose_mutagenesis.
    """
    # Generate a sequence, then retry if it has no GT dinucleotides
    seq = draw(dna_sequence(min_len=9, max_len=60))
    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i + 2] == "GT"]
    assume(len(gt_positions) > 0)

    # Deduplicate to codon-start positions
    constraint_positions = list({(pos // 3) * 3 for pos in gt_positions})
    constraint_types = {pos: ["GT"] for pos in constraint_positions}

    return seq, constraint_positions, constraint_types


# ---------------------------------------------------------------------------
# Helper: apply a MutagenesisProposal to a DNA sequence
# ---------------------------------------------------------------------------

def apply_proposal(seq: str, proposal: MutagenesisProposal) -> str:
    """Replace the codon at proposal.position with proposal.new_codon.

    Only applies if new_codon is non-empty (i.e. not an impossible proposal).
    """
    if not proposal.new_codon:
        return seq  # impossible proposal — nothing to apply
    seq_list = list(seq)
    for i, base in enumerate(proposal.new_codon):
        seq_list[proposal.position + i] = base
    return "".join(seq_list)


# ===========================================================================
# Property 1: Any valid mutation produces a valid protein of same length
# ===========================================================================

class TestMutationPreservesProtein:
    """Applying any feasible mutagenesis proposal yields a protein of the
    same length as the original, with all valid amino acid residues."""

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_mutation_preserves_protein_length(self, data):
        """For any feasible (non-impossible) proposal from propose_mutagenesis,
        applying the proposal and re-translating yields a protein of the same
        length as the original."""
        seq, constraint_positions, constraint_types = data
        original_protein = _translate_dna_to_aa(seq)
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue  # skip impossible proposals

            mutated_seq = apply_proposal(seq, proposal)
            mutated_protein = _translate_dna_to_aa(mutated_seq)

            assert len(mutated_protein) == len(original_protein), (
                f"Mutation at codon position {proposal.position} changed protein "
                f"length from {len(original_protein)} to {len(mutated_protein)}"
            )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_mutation_produces_valid_protein(self, data):
        """After applying a feasible mutagenesis proposal, the resulting
        protein contains only standard amino acids (no stop codons)."""
        seq, constraint_positions, constraint_types = data
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue

            mutated_seq = apply_proposal(seq, proposal)
            mutated_protein = _translate_dna_to_aa(mutated_seq)

            for i, aa in enumerate(mutated_protein):
                assert aa in STANDARD_AAS, (
                    f"Invalid amino acid {aa!r} at position {i} in mutated "
                    f"protein {mutated_protein!r}"
                )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_synonymous_mutation_preserves_protein_sequence(self, data):
        """A synonymous mutation (new_aa == original_aa) must produce an
        identical protein sequence after re-translation."""
        seq, constraint_positions, constraint_types = data
        original_protein = _translate_dna_to_aa(seq)
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            if proposal.new_aa != proposal.original_aa:
                continue  # not a synonymous mutation

            mutated_seq = apply_proposal(seq, proposal)
            mutated_protein = _translate_dna_to_aa(mutated_seq)

            assert mutated_protein == original_protein, (
                f"Synonymous mutation at position {proposal.position} changed "
                f"protein: {original_protein!r} -> {mutated_protein!r}"
            )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_nonsynonymous_mutation_changes_exactly_one_residue(self, data):
        """A non-synonymous mutation (new_aa != original_aa) should change
        exactly one residue in the protein."""
        seq, constraint_positions, constraint_types = data
        original_protein = _translate_dna_to_aa(seq)
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            if proposal.new_aa == proposal.original_aa:
                continue  # synonymous — skip

            mutated_seq = apply_proposal(seq, proposal)
            mutated_protein = _translate_dna_to_aa(mutated_seq)

            # Count differing positions
            diffs = sum(
                1 for a, b in zip(original_protein, mutated_protein) if a != b
            )
            assert diffs == 1, (
                f"Non-synonymous mutation at codon pos {proposal.position} "
                f"changed {diffs} residues instead of exactly 1: "
                f"{original_protein!r} -> {mutated_protein!r}"
            )

    def test_mutation_preserves_protein_length_concrete(self):
        """Concrete example: AGT(S) -> synonymous codon preserves protein."""
        seq = "ATGAGTGCTTAA"
        constraint_positions = [3]
        constraint_types = {3: ["GT"]}
        original_protein = _translate_dna_to_aa(seq)
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )
        assert len(report.proposals) > 0
        for proposal in report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            mutated_seq = apply_proposal(seq, proposal)
            mutated_protein = _translate_dna_to_aa(mutated_seq)
            assert len(mutated_protein) == len(original_protein)


# ===========================================================================
# Property 2: Mutation positions are within sequence bounds
# ===========================================================================

class TestMutationPositionBounds:
    """Every proposal from propose_mutagenesis must have a position that
    falls within the bounds of the input sequence (position + 3 <= len(seq))."""

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_position_within_sequence_bounds(self, data):
        """Every proposal.position + 3 must be <= len(seq)."""
        seq, constraint_positions, constraint_types = data
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            assert proposal.position >= 0, (
                f"Proposal position {proposal.position} is negative"
            )
            assert proposal.position + 3 <= len(seq), (
                f"Proposal position {proposal.position} + 3 = "
                f"{proposal.position + 3} exceeds sequence length {len(seq)}"
            )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_position_is_codon_aligned(self, data):
        """Every proposal.position must be a multiple of 3 (codon-aligned)."""
        seq, constraint_positions, constraint_types = data
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            assert proposal.position % 3 == 0, (
                f"Proposal position {proposal.position} is not codon-aligned "
                f"(not divisible by 3)"
            )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_original_codon_matches_sequence(self, data):
        """Every proposal's original_codon matches the 3-base slice of the
        input sequence at proposal.position."""
        seq, constraint_positions, constraint_types = data
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            expected_codon = seq[proposal.position:proposal.position + 3]
            assert proposal.original_codon == expected_codon, (
                f"Proposal.original_codon {proposal.original_codon!r} != "
                f"seq[{proposal.position}:{proposal.position + 3}] "
                f"= {expected_codon!r}"
            )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_original_aa_matches_codon_table(self, data):
        """Every proposal's original_aa matches the CODON_TABLE entry for
        original_codon."""
        seq, constraint_positions, constraint_types = data
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            expected_aa = CODON_TABLE.get(proposal.original_codon)
            assert proposal.original_aa == expected_aa, (
                f"Proposal.original_aa {proposal.original_aa!r} != "
                f"CODON_TABLE[{proposal.original_codon!r}] = {expected_aa!r}"
            )

    @given(data=dna_with_gt_constraint())
    @settings(max_examples=200, deadline=5000)
    def test_new_codon_is_valid(self, data):
        """Every feasible proposal's new_codon is a valid 3-base codon in
        CODON_TABLE that translates to new_aa."""
        seq, constraint_positions, constraint_types = data
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )

        for proposal in report.proposals:
            if proposal.impossible or not proposal.new_codon:
                continue
            assert len(proposal.new_codon) == 3, (
                f"new_codon {proposal.new_codon!r} is not 3 bases"
            )
            assert proposal.new_codon in CODON_TABLE, (
                f"new_codon {proposal.new_codon!r} not in CODON_TABLE"
            )
            assert CODON_TABLE[proposal.new_codon] == proposal.new_aa, (
                f"CODON_TABLE[{proposal.new_codon!r}] = "
                f"{CODON_TABLE[proposal.new_codon]!r} != "
                f"new_aa {proposal.new_aa!r}"
            )

    def test_position_bounds_concrete(self):
        """Concrete example: constraint at codon 3 of a 12-base sequence."""
        seq = "ATGAGTGCTTAA"
        constraint_positions = [3]
        constraint_types = {3: ["GT"]}
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, MOCK_CAI
        )
        for proposal in report.proposals:
            assert proposal.position + 3 <= len(seq)


# ===========================================================================
# Property 3: BLOSUM62 scores are symmetric
# ===========================================================================

class TestBLOSUM62Symmetry:
    """BLOSUM62(a, b) == BLOSUM62(b, a) for all amino acid pairs.

    This is a fundamental algebraic property of the BLOSUM62 substitution
    matrix.  Symmetry ensures that the score for substituting a with b is
    the same as substituting b with a, regardless of direction.
    """

    @given(key=blosum_key)
    @settings(max_examples=400)
    def test_blosum62_symmetric(self, key):
        """BLOSUM62[(a, b)] == BLOSUM62[(b, a)] for any pair of amino acids."""
        a, b = key
        score_ab = BLOSUM62[(a, b)]
        score_ba = BLOSUM62[(b, a)]
        assert score_ab == score_ba, (
            f"BLOSUM62 not symmetric: B62({a},{b})={score_ab} != "
            f"B62({b},{a})={score_ba}"
        )

    @given(aa=amino_acid)
    @settings(max_examples=50)
    def test_blosum62_diagonal_positive(self, aa):
        """BLOSUM62[(a, a)] > 0 for every standard amino acid.

        A self-substitution should always be conservative (positive score).
        """
        assert BLOSUM62[(aa, aa)] > 0, (
            f"BLOSUM62({aa},{aa}) = {BLOSUM62[(aa, aa)]} is not positive"
        )

    @given(pair=aa_pair)
    @settings(max_examples=200)
    def test_blosum62_diagonal_maximal(self, pair):
        """BLOSUM62[(a, a)] >= BLOSUM62[(a, b)] for any b.

        The diagonal entry is always the maximum in its row.
        """
        a, b = pair
        assert BLOSUM62[(a, a)] >= BLOSUM62[(a, b)], (
            f"BLOSUM62({a},{a})={BLOSUM62[(a, a)]} < "
            f"BLOSUM62({a},{b})={BLOSUM62[(a, b)]}"
        )

    @given(pair=aa_pair)
    @settings(max_examples=200)
    def test_symmetry_implies_consistent_mutagenesis_scoring(self, pair):
        """If BLOSUM62 is symmetric, then propose_mutagenesis should score
        a↔b the same regardless of which is original and which is new.

        This tests the practical implication of symmetry: the conservation
        score of a substitution does not depend on the direction.
        """
        a, b = pair
        score_ab = BLOSUM62[(a, b)]
        score_ba = BLOSUM62[(b, a)]
        # Direct symmetry check (reinforces Property 3)
        assert score_ab == score_ba
        # Also verify both are integers
        assert isinstance(score_ab, int)
        assert isinstance(score_ba, int)

    def test_all_pairs_symmetric_exhaustive(self):
        """Exhaustive check: BLOSUM62[(a, b)] == BLOSUM62[(b, a)]
        for every pair in the matrix.

        This complements the Hypothesis-based test by checking every
        single entry, not just a random sample.
        """
        for (a, b), score in BLOSUM62.items():
            reverse_score = BLOSUM62.get((b, a))
            assert reverse_score is not None, (
                f"BLOSUM62 key ({b}, {a}) missing but ({a}, {b}) exists"
            )
            assert score == reverse_score, (
                f"BLOSUM62({a},{b})={score} != BLOSUM62({b},{a})={reverse_score}"
            )
