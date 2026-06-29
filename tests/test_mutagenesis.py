"""Test BioCompiler Mutagenesis Engine — proposals, chose_poorly/impossible detection."""

import pytest
from biocompiler.optimizer.mutagenesis import (
    propose_mutagenesis, MutagenesisProposal, MutagenesisReport,
)
from biocompiler.type_system import AA_TO_CODONS, CODON_TABLE
from biocompiler.organisms import ECOLI_CAI


class TestMutagenesisReport:
    """Tests for MutagenesisReport data structure."""

    def test_empty_report(self):
        """Empty report has no chose_poorly or impossible."""
        report = MutagenesisReport()
        assert report.has_chose_poorly is False
        assert report.has_impossible is False
        assert len(report.proposals) == 0

    def test_report_with_impossible(self):
        """Report with an impossible proposal detects it."""
        proposal = MutagenesisProposal(
            position=0,
            original_codon="GTT",
            original_aa="V",
            new_aa="",
            new_codon="",
            blosum_score=-10,
            cai_weight=0.0,
            resolves=["GT"],
            impossible=True,
        )
        report = MutagenesisReport()
        report.add(proposal)
        assert report.has_impossible is True
        assert report.has_chose_poorly is False


class TestProposeMutagenesis:
    """Tests for the propose_mutagenesis function."""

    def test_propose_mutagenesis_synonymous(self):
        """For a within-codon GT in a non-Valine AA, propose synonymous codon first.

        Serine has codons TCT, TCC, TCA, TCG, AGT, AGC.
        TCG contains GT... wait, no. AGT contains GT at positions 1-2.
        Let us use a Serine codon AGT which has GT.
        """
        # Build a short sequence: ATG (M) AGT (S, has GT) GCT (A) TAA (stop)
        # Constraint at codon position 3 (AGT, Serine)
        seq = "ATGAGTGCTTAA"
        # The AGT codon starts at position 3
        constraint_positions = [3]
        constraint_types = {3: ["GT"]}

        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, ECOLI_CAI
        )

        # There should be at least one proposal
        assert len(report.proposals) > 0

        # The proposal should suggest a synonymous codon (same AA = Serine)
        proposal = report.proposals[0]
        # Serine has TCT, TCC, TCA, TCG (no GT), AGT, AGC (no GT)
        # So there should be a synonymous codon without GT
        assert proposal.original_aa == "S"
        assert proposal.new_aa == "S"  # synonymous
        assert "GT" not in proposal.new_codon

    def test_propose_mutagenesis_impossible(self):
        """When no substitution can resolve the constraint, reports impossible.

        We create a very constrained scenario where no amino acid substitution
        can resolve the GT constraint because all alternatives also contain GT
        or are blocked.
        """
        # For an impossible case, we need a scenario where:
        # - No synonymous codon resolves the constraint
        # - No conservative AA substitution resolves the constraint
        # This is hard to achieve naturally, so we test with very high min_blosum
        # and min_cai thresholds to force impossible.

        # Use a Valine codon (all codons have GT) with very restrictive thresholds
        seq = "ATGGTTGCTTAA"  # GTT = Valine (all Val codons have GT)
        constraint_positions = [3]
        constraint_types = {3: ["GT"]}

        # With very high min_blosum, no AA substitution will pass
        report = propose_mutagenesis(
            seq, constraint_positions, constraint_types, ECOLI_CAI,
            min_blosum=100,  # Impossibly high - no substitution can achieve this
            min_cai=0.99,    # Also very restrictive
        )

        # With such restrictive thresholds, should report impossible
        assert len(report.proposals) > 0
        # The proposal should be flagged as impossible
        assert report.has_impossible is True
