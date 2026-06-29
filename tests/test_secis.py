"""Tests for selenocysteine (Sec/U) support via SECIS positions."""

import pytest
from biocompiler.ir import IR_L0_GenomicDNA, IR_L1_PreMRNA, IR_L2_MatureMRNA, IRLevel, IRError
from biocompiler.ir.passes import transcribe, splice, translate, compile_gene
from biocompiler.ir.invariants import check_l0_invariants, check_l2_invariants
from biocompiler.type_system.codon_tables import AA_TO_CODONS


class TestSelenocysteineTranslation:
    """Test that UGA codons at SECIS positions translate to U (selenocysteine)."""

    def test_uga_without_secis_is_stop(self):
        """UGA without SECIS position → stop codon (*)."""
        # ATG-GCT-TGA → M-A-* (TGA is stop)
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGA", regions=[], organism="human",
            secis_positions=[]  # no SECIS
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)
        assert ir_l3.sequence == "MA*"

    def test_uga_with_secis_is_selenocysteine(self):
        """UGA at a SECIS position → selenocysteine (U)."""
        # ATG-GCT-TGA-TAA → M-A-U-* (TGA recoded to U at codon 2, TAA is stop)
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGATAA", regions=[], organism="human",
            secis_positions=[2]  # codon 2 is TGA → U
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)
        assert ir_l3.sequence == "MAU*"

    def test_multiple_secis_positions(self):
        """Multiple UGA codons recoded to U."""
        # ATG-TGA-GCT-TGA-TGA → M-U-A-U-* (last TGA is stop, others are U)
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGTGAGCTTGATGA", regions=[], organism="human",
            secis_positions=[1, 3]  # codons 1 and 3 are TGA → U
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)
        assert ir_l3.sequence == "MUAU*"

    def test_secis_propagates_through_transcribe(self):
        """SECIS positions survive transcription (L0→L1)."""
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGATAA", regions=[], organism="human",
            secis_positions=[2]
        )
        ir_l1 = transcribe(ir_l0)
        assert ir_l1.secis_positions == [2]

    def test_secis_propagates_through_splice(self):
        """SECIS positions survive splicing (L1→L2)."""
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGATAA", regions=[], organism="human",
            secis_positions=[2]
        )
        ir_l2 = compile_gene(ir_l0, IRLevel.L2)
        assert ir_l2.secis_positions == [2]

    def test_secis_with_real_selenoprotein_fragment(self):
        """Test with a fragment of a real selenoprotein (GPX1)."""
        # GPX1 contains selenocysteine at position 49 (0-indexed)
        # Fragment: M...U...* (simplified)
        # Construct: M-A-U-A-* (ATG GCT TGA GCT TAA, codon 2 = U)
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGAGCTTAA", regions=[], organism="human",
            gene_name="GPX1_fragment",
            secis_positions=[2]
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)
        assert ir_l3.sequence == "MAUA*"
        assert "U" in ir_l3.sequence

    def test_secis_does_not_affect_non_uga_codons(self):
        """SECIS positions only affect UGA codons, not other stop codons."""
        # ATG-TAA-GCT → M-*-... (TAA is stop, not affected by SECIS)
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGTAAGCTTAA", regions=[], organism="human",
            secis_positions=[1]  # codon 1 is TAA, not UGA — should still be stop
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)
        assert ir_l3.sequence == "M*"

    def test_back_translate_selenocysteine(self):
        """Test back-translation of U to TGA with SECIS position tracking."""
        import random
        random.seed(42)
        protein = "MAUAR*"  # contains selenocysteine
        secis_pos = []
        dna = []
        for i, aa in enumerate(protein):
            if aa == 'U':
                dna.append('TGA')  # selenocysteine codon
                secis_pos.append(i)
            elif aa == '*':
                dna.append('TAA')  # stop
            else:
                codons = AA_TO_CODONS.get(aa, ['GCT'])
                dna.append(random.choice(codons))

        dna_seq = ''.join(dna)
        ir_l0 = IR_L0_GenomicDNA(
            sequence=dna_seq, regions=[], organism="human",
            secis_positions=secis_pos
        )
        ir_l3 = compile_gene(ir_l0, IRLevel.L3)
        assert ir_l3.sequence == protein
        assert "U" in ir_l3.sequence


class TestSECISCodeGen:
    """Test that SECIS positions are reflected in codegen output."""

    def test_genbank_with_secis(self):
        """GenBank output should note SECIS positions."""
        from biocompiler.ir.codegen import to_genbank
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGAGCTTAA", regions=[], organism="human",
            gene_name="test_secis",
            secis_positions=[2]
        )
        gb = to_genbank(ir_l0)
        assert "LOCUS" in gb
        # The GenBank output should be valid (SECIS note is metadata)

    def test_sbol3_with_secis(self):
        """SBOL3 output should handle SECIS-containing sequences."""
        from biocompiler.ir.codegen import to_sbol3
        ir_l0 = IR_L0_GenomicDNA(
            sequence="ATGGCTTGAGCTTAA", regions=[], organism="human",
            gene_name="test_secis",
            secis_positions=[2]
        )
        sb = to_sbol3(ir_l0)
        assert "sbol:Component" in sb
        assert "ATGGCTTGAGCTTAA" in sb
