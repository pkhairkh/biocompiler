"""Test that optimize_sequence() now runs the full IR pipeline.

These tests verify the TIGHTEN-1 change: the main ``optimize_sequence()``
API has been wired into the BioCompiler IR module so that every call
runs the optimized DNA through the L0→L1→L2→L3 lowering passes, verifies
that the IR-translated protein matches the input protein, and attaches
the IR objects + codegen outputs (GenBank / FASTA / SBOL3) to the
``OptimizationResult``.

Before TIGHTEN-1 the IR module was real but isolated — users had to
call ``compile_gene(ir_l0, IRLevel.L3)`` manually.  After TIGHTEN-1
the "compiler" claim is true for every ``optimize_sequence()`` call.
"""
import pytest

from biocompiler.optimizer.pipeline_core import optimize_sequence


class TestIRIntegration:
    """Verify optimize_sequence() runs the IR pipeline and attaches
    IR objects + codegen outputs to the OptimizationResult."""

    def test_optimize_returns_ir_objects(self):
        """Result must carry non-None ir_l0 and ir_l3 IR objects."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.ir_l0 is not None
        assert result.ir_l3 is not None
        assert result.ir_verified is True

    def test_optimize_produces_genbank(self):
        """GenBank codegen output must be populated and well-formed."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.genbank is not None
        assert "LOCUS" in result.genbank

    def test_optimize_produces_fasta(self):
        """FASTA codegen output must be populated and well-formed."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.fasta is not None
        assert ">" in result.fasta

    def test_optimize_produces_sbol3(self):
        """SBOL3 codegen output must be populated and well-formed."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.sbol3 is not None
        assert "sbol:" in result.sbol3

    def test_ir_verifies_protein_preservation(self):
        """IR-translated protein must start with the input protein
        and ir_verified must be True for a valid optimization."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.ir_verified is True
        assert result.ir_l3.sequence.startswith("MVHLTPEEK")

    def test_hbb_full_pipeline(self):
        """Full HBB N-terminus must round-trip through the IR.

        The IR-translated polypeptide must equal the input protein
        plus a trailing '*' (the stop codon that the integration code
        appends before passing the sequence to compile_gene).
        """
        hbb = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
        result = optimize_sequence(hbb, organism="human", strict_mode=False)
        assert result.ir_verified is True
        assert result.ir_l3.sequence == hbb + "*"
        assert result.genbank is not None
        assert result.fasta is not None

    def test_ir_l0_sequence_matches_optimized_dna(self):
        """The IR-L0 object's sequence must be the optimized DNA plus
        an appended stop codon (when the optimizer did not include one)."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        # ir_l0.sequence should end with a stop codon
        assert result.ir_l0.sequence[-3:] in {"TAA", "TAG", "TGA"}
        # The first len(protein)*3 bases must match result.sequence
        prefix_len = len(result.protein) * 3
        assert result.ir_l0.sequence[:prefix_len] == result.sequence

    def test_ir_l3_level_is_polypeptide(self):
        """The IR-L3 object's level must be IRLevel.L3."""
        from biocompiler.ir import IRLevel
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.ir_l3.level == IRLevel.L3

    def test_ir_l0_level_is_genomic_dna(self):
        """The IR-L0 object's level must be IRLevel.L0."""
        from biocompiler.ir import IRLevel
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert result.ir_l0.level == IRLevel.L0

    def test_ir_provenance_trail_is_stamped(self):
        """IR objects must carry the pass/source_level provenance trail
        produced by the lowering passes (transcribe, splice, translate)."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        # The L3 object's metadata should record that it was produced
        # by the 'translate' pass from 'L2'.
        meta = result.ir_l3.metadata
        assert meta.get("pass") == "translate"
        assert meta.get("source_level") == "L2"

    def test_genbank_contains_organism(self):
        """GenBank output must reference the target organism."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        # The codegen _organism_display() resolves 'human' → 'Homo sapiens'
        assert "Homo sapiens" in result.genbank

    def test_fasta_header_contains_level(self):
        """FASTA header must identify the IR level (L3 = protein)."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert "|L3|" in result.fasta

    def test_sbol3_contains_dna_elements(self):
        """SBOL3 output must include the DNA sequence in sbol:elements."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        assert "sbol:elements" in result.sbol3
        # The DNA sequence (uppercase) must appear in the SBOL3 document
        assert result.sequence[:9] in result.sbol3  # first codon (ATG) is enough

    def test_prokaryote_also_gets_ir(self):
        """Prokaryotic optimization must also produce IR objects."""
        result = optimize_sequence("MVHLTPEEK", organism="ecoli", strict_mode=False)
        assert result.ir_l0 is not None
        assert result.ir_l3 is not None
        assert result.ir_verified is True
        assert result.ir_l3.sequence == "MVHLTPEEK*"

    def test_ir_does_not_modify_result_sequence(self):
        """The IR pipeline must not modify result.sequence — the
        appended stop codon is for IR/codegen purposes only."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        # result.sequence length must equal len(protein) * 3 (no stop)
        assert len(result.sequence) == len(result.protein) * 3
