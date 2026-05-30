#!/usr/bin/env python3
"""
Unit tests for BioCompiler biological data module.
"""

import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from biocompiler_data import (
    HBB_PREMRNA, HBB_MRNA, HBB_CDS, HBB_PROTEIN,
    HBB_EXONS, HBB_INTRONS, HBB_SPLICE_SITES,
    EGFP_CDS, EGFP_PROTEIN,
    CODON_USAGE, CODON_ADAPTIVENESS, PREFERRED_CODONS,
    splice_premrna, get_intron_sequence, translate_dna,
)

import pytest


class TestHBBPreMRNALength:
    def test_premrna_length(self):
        """HBB pre-mRNA should be 1608 bp."""
        assert len(HBB_PREMRNA) == 1608

    def test_mrna_length(self):
        """HBB mRNA should be 628 bp."""
        assert len(HBB_MRNA) == 628

    def test_cds_length(self):
        """HBB CDS should be 444 bp (147 aa * 3)."""
        assert len(HBB_CDS) == 444

    def test_protein_length(self):
        """HBB protein should be 147 aa."""
        assert len(HBB_PROTEIN) == 147


class TestSplicingProducesMRNA:
    def test_splice_produces_mrna(self):
        """Splicing pre-mRNA with exon boundaries should produce mRNA."""
        mrna = splice_premrna(HBB_PREMRNA, HBB_EXONS)
        assert mrna == HBB_MRNA

    def test_exon_sizes_sum_to_mrna(self):
        """Exon sizes should sum to mRNA length."""
        total_exon = sum(length for _, _, length in HBB_EXONS)
        assert total_exon == len(HBB_MRNA)


class TestCodonTableCompleteness:
    def test_codon_usage_count(self):
        """Codon usage table should have 64 entries."""
        assert len(CODON_USAGE) == 64

    def test_codon_adaptiveness_count(self):
        """Codon adaptiveness should have 61 sense codons."""
        assert len(CODON_ADAPTIVENESS) == 61

    def test_preferred_codons_per_aa(self):
        """Preferred codons should be defined for standard amino acids."""
        standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
        for aa in standard_aas:
            assert aa in PREFERRED_CODONS, f"Missing preferred codon for {aa}"

    def test_cds_starts_with_atg(self):
        """CDS should start with ATG."""
        assert HBB_CDS[:3] == "ATG"

    def test_cds_translates_to_protein(self):
        """HBB CDS should translate to expected protein."""
        protein = translate_dna(HBB_CDS)
        assert protein == HBB_PROTEIN

    def test_egfp_length(self):
        """EGFP CDS should be 720 bp."""
        assert len(EGFP_CDS) == 720

    def test_egfp_translates(self):
        """EGFP CDS should translate correctly."""
        protein = translate_dna(EGFP_CDS)
        assert protein == EGFP_PROTEIN

    def test_intron_gt_ag(self):
        """Introns should follow GT...AG rule."""
        intron1 = get_intron_sequence(HBB_PREMRNA, 143, 272)
        intron2 = get_intron_sequence(HBB_PREMRNA, 496, 1345)
        assert intron1[:2] == "GT" and intron1[-2:] == "AG"
        assert intron2[:2] == "GT" and intron2[-2:] == "AG"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
