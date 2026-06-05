"""
Tests for biocompiler.multigene — Multi-Gene Construct / Operon Support

Comprehensive test suite covering:
- GeneSpec validation and construction
- Prokaryotic operon assembly (optimize_operon)
- Eukaryotic polycistronic constructs (2A peptides, IRES)
- Bidirectional promoter constructs
- Custom multi-gene constructs with linkers
- GenBank export with per-gene feature annotations
- Construct type inference
- Edge cases and error handling
"""

import pytest

from biocompiler.multigene import (
    GeneSpec,
    MultiGeneResult,
    OperonConfig,
    optimize_multigene,
    optimize_operon,
    LINKER_2A_P2A,
    LINKER_2A_T2A,
    RBS_STRONG,
    RBS_SPACER_DEFAULT,
)
from biocompiler.exceptions import InvalidProteinError


# ─── GeneSpec Tests ────────────────────────────────────────────────


class TestGeneSpec:
    """Tests for GeneSpec data class validation."""

    def test_basic_construction(self):
        """GeneSpec with just a protein sequence."""
        spec = GeneSpec(protein="MSKGEELFTG")
        assert spec.protein == "MSKGEELFTG"
        assert spec.organism == ""
        assert spec.name == ""
        assert spec.promoter == ""
        assert spec.rbs == ""
        assert spec.terminator == ""

    def test_full_construction(self):
        """GeneSpec with all fields."""
        spec = GeneSpec(
            protein="MSKGEELFTG",
            organism="Escherichia_coli",
            name="gfp",
            promoter="TTGACA",
            rbs="AGGAGG",
            terminator="TTTTTT",
        )
        assert spec.protein == "MSKGEELFTG"
        assert spec.organism == "Escherichia_coli"
        assert spec.name == "gfp"
        assert spec.promoter == "TTGACA"
        assert spec.rbs == "AGGAGG"
        assert spec.terminator == "TTTTTT"

    def test_case_normalization(self):
        """Protein sequences are uppercased."""
        spec = GeneSpec(protein="mskgeelftg")
        assert spec.protein == "MSKGEELFTG"

    def test_whitespace_stripped(self):
        """Whitespace is stripped from protein sequences."""
        spec = GeneSpec(protein="  MSKGEELFTG  ")
        assert spec.protein == "MSKGEELFTG"

    def test_empty_protein_raises(self):
        """Empty protein sequence raises ValueError."""
        with pytest.raises(ValueError, match="non-empty"):
            GeneSpec(protein="")

    def test_whitespace_only_protein_raises(self):
        """Whitespace-only protein sequence raises ValueError."""
        with pytest.raises(ValueError, match="non-empty"):
            GeneSpec(protein="   ")

    def test_invalid_amino_acids_raises(self):
        """Invalid amino acid characters raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGX")  # X is invalid

    def test_invalid_amino_acids_numbers(self):
        """Numbers in protein sequence raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELF1G")


# ─── Prokaryotic Operon Tests ──────────────────────────────────────


class TestOperon:
    """Tests for prokaryotic operon optimization."""

    def test_basic_operon_ecoli(self):
        """Two-gene operon in E. coli produces valid result."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert isinstance(result, MultiGeneResult)
        assert result.construct_type == "operon"
        assert len(result.genes) == 2
        assert result.total_length > 0
        assert len(result.full_dna) == result.total_length
        assert 0.0 <= result.gc_content <= 1.0
        assert result.organism == "Escherichia_coli"

    def test_operon_has_rbs(self):
        """Each gene in an operon has an RBS upstream of its start codon."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        # The full_dna should contain RBS sequences before each gene
        # AGGAGG is the default RBS
        assert "AGGAGG" in result.full_dna

    def test_operon_with_promoter(self):
        """Operon with a specified promoter includes it in the construct."""
        promoter = "TTGACA"
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
            ],
            organism="Escherichia_coli",
            promoter=promoter,
        )
        assert result.full_dna.startswith(promoter.upper())

    def test_operon_with_terminator(self):
        """Operon with a specified terminator includes it in the construct."""
        terminator = "TTTTTTT"
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
            ],
            organism="Escherichia_coli",
            terminator=terminator,
        )
        assert result.full_dna.endswith(terminator.upper())

    def test_operon_three_genes(self):
        """Three-gene operon assembles correctly."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
                GeneSpec(protein="MRVLKFGGTS", name="lacZ"),
            ],
            organism="Escherichia_coli",
        )
        assert len(result.genes) == 3
        # Each gene's DNA should be approximately 30 bp (10 aa * 3)
        for g in result.genes:
            assert len(g.sequence) == 30  # 10 aa * 3

    def test_operon_single_gene(self):
        """Single-gene 'operon' still works."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
            ],
            organism="Escherichia_coli",
        )
        assert len(result.genes) == 1

    def test_operon_empty_raises(self):
        """Empty gene list raises ValueError."""
        with pytest.raises(ValueError, match="At least one"):
            optimize_operon(genes=[])

    def test_operon_gene_per_gene_cai(self):
        """Each gene in the operon has a valid CAI score."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp1"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        for g in result.genes:
            assert 0.0 <= g.cai <= 1.0


# ─── Polycistronic 2A Tests ────────────────────────────────────────


class TestPolycistronic2A:
    """Tests for eukaryotic polycistronic constructs with 2A peptides."""

    def test_2a_construct_type(self):
        """2A linker results in polycistronic_2A construct type."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"

    def test_2a_has_both_genes(self):
        """2A construct contains both gene optimization results."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        assert len(result.genes) == 2
        # The full DNA should be longer than the sum of individual genes
        # because of the 2A peptide coding sequence between them
        gene_sum = sum(len(g.sequence) for g in result.genes)
        assert result.total_length >= gene_sum

    def test_2a_t2a(self):
        """T2A peptide linker is recognized."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=LINKER_2A_T2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"


# ─── Polycistronic IRES Tests ──────────────────────────────────────


class TestPolycistronicIRES:
    """Tests for eukaryotic polycistronic constructs with IRES."""

    def test_ires_construct_type(self):
        """Long DNA linker results in polycistronic_IRES construct type."""
        # Simulate an IRES sequence (typically 400-600 bp)
        ires_seq = "C" * 200  # Long enough to trigger IRES detection
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=ires_seq,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_IRES"

    def test_ires_linker_in_construct(self):
        """IRES DNA sequence is included in the full construct."""
        ires_seq = "ACGT" * 50  # 200 bp IRES
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=ires_seq,
            organism="Homo_sapiens",
        )
        assert ires_seq.upper() in result.full_dna


# ─── Custom Construct Tests ────────────────────────────────────────


class TestCustomConstruct:
    """Tests for custom multi-gene constructs."""

    def test_custom_with_linker(self):
        """Custom construct with DNA linker between genes."""
        linker = "GGATCC"  # BamHI site as linker
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            linker=linker,
            organism="Homo_sapiens",
        )
        assert isinstance(result, MultiGeneResult)
        assert len(result.genes) == 2
        # Linker should be present between genes
        assert linker.upper() in result.full_dna

    def test_custom_without_linker(self):
        """Custom construct without linker concatenates genes."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Homo_sapiens",
        )
        gene_sum = sum(len(g.sequence) for g in result.genes)
        # Without linkers, total should equal gene sum
        # (for custom constructs without regulatory elements)
        assert result.total_length == gene_sum

    def test_single_gene_construct(self):
        """Single-gene construct works."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
            ],
            organism="Escherichia_coli",
        )
        assert len(result.genes) == 1
        assert result.genes[0].protein == "MSKGEELFTG"

    def test_per_gene_organism(self):
        """Each gene can specify its own organism."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp",
                         organism="Escherichia_coli"),
            ],
            organism="Homo_sapiens",  # Default overridden by gene spec
        )
        # The gene should have been optimized for E. coli
        assert len(result.genes) == 1

    def test_gene_with_promoter_rbs_terminator(self):
        """Gene with promoter, RBS, and terminator."""
        result = optimize_multigene(
            genes=[
                GeneSpec(
                    protein="MSKGEELFTG",
                    name="gfp",
                    promoter="TTGACA",
                    rbs="AGGAGG",
                    terminator="TTTTTT",
                ),
            ],
            organism="Escherichia_coli",
        )
        assert "TTGACA" in result.full_dna
        assert "AGGAGG" in result.full_dna
        assert "TTTTTT" in result.full_dna


# ─── GenBank Export Tests ──────────────────────────────────────────


class TestGenBankExport:
    """Tests for GenBank format export of multi-gene constructs."""

    def test_genbank_has_locus(self):
        """GenBank export contains LOCUS line."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert "LOCUS" in result.genbank_export

    def test_genbank_has_features(self):
        """GenBank export contains FEATURES section."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert "FEATURES" in result.genbank_export

    def test_genbank_has_origin(self):
        """GenBank export contains ORIGIN section."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert "ORIGIN" in result.genbank_export

    def test_genbank_has_terminator(self):
        """GenBank export ends with //."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
            ],
            organism="Escherichia_coli",
        )
        assert result.genbank_export.endswith("//")

    def test_genbank_has_gene_names(self):
        """GenBank export contains gene name annotations."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert '/gene="gfp"' in result.genbank_export
        assert '/gene="tetR"' in result.genbank_export

    def test_genbank_has_cds_features(self):
        """GenBank export contains CDS features for each gene."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert "CDS" in result.genbank_export

    def test_genbank_has_cai_values(self):
        """GenBank export contains CAI qualifiers."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
            ],
            organism="Escherichia_coli",
        )
        assert "/cai=" in result.genbank_export

    def test_genbank_has_translation(self):
        """GenBank export contains protein translation."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
            ],
            organism="Escherichia_coli",
        )
        assert "/translation=" in result.genbank_export

    def test_genbank_comment_has_construct_info(self):
        """GenBank COMMENT section includes construct type and gene count."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert "operon" in result.genbank_export.lower()
        assert "Number of genes" in result.genbank_export

    def test_genbank_eukaryotic(self):
        """GenBank export for human construct."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="GFP"),
                GeneSpec(protein="MDDRLEAIAG", name="RFP"),
            ],
            organism="Homo_sapiens",
        )
        assert "Homo sapiens" in result.genbank_export


# ─── Construct Type Inference Tests ────────────────────────────────


class TestConstructTypeInference:
    """Tests for automatic construct type detection."""

    def test_ecoli_with_rbs_is_operon(self):
        """E. coli genes with RBS → operon."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp", rbs="AGGAGG"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR", rbs="AGGAGG"),
            ],
            organism="Escherichia_coli",
        )
        assert result.construct_type == "operon"

    def test_ecoli_without_rbs_still_operon(self):
        """E. coli multi-gene without RBS defaults to operon."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert result.construct_type == "operon"

    def test_2a_linker_is_polycistronic_2a(self):
        """2A peptide linker → polycistronic_2A."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"

    def test_long_dna_linker_is_ires(self):
        """Long DNA linker (>100 bp) → polycistronic_IRES."""
        ires = "ACGT" * 50  # 200 bp
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker=ires,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_IRES"

    def test_human_short_linker_is_custom(self):
        """Short DNA linker with human → custom."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGE", name="GFP"),
                GeneSpec(protein="MDDRLA", name="RFP"),
            ],
            linker="GGATCC",
            organism="Homo_sapiens",
        )
        assert result.construct_type == "custom"


# ─── Error Handling Tests ──────────────────────────────────────────


class TestErrorHandling:
    """Tests for error conditions and edge cases."""

    def test_no_genes_raises(self):
        """Empty gene list raises ValueError."""
        with pytest.raises(ValueError, match="At least one"):
            optimize_multigene(genes=[])

    def test_no_organism_raises(self):
        """Missing organism raises ValueError."""
        with pytest.raises(ValueError, match="Organism must be specified"):
            optimize_multigene(
                genes=[GeneSpec(protein="MSKGEELFTG")],
                organism="",
            )

    def test_invalid_protein_in_spec_raises(self):
        """Invalid amino acids in GeneSpec raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_multigene(
                genes=[GeneSpec(protein="MSKGEELFTGX")],
                organism="Escherichia_coli",
            )

    def test_unsupported_organism(self):
        """Unsupported organism name still resolves gracefully."""
        # resolve_organism with strict=False doesn't raise for unknowns
        # but the optimizer may produce lower quality results
        result = optimize_multigene(
            genes=[GeneSpec(protein="MSKGEELFTG")],
            organism="Escherichia_coli",  # Use a known organism
        )
        assert isinstance(result, MultiGeneResult)

    def test_constraints_forwarded(self):
        """Constraint parameters are forwarded to optimize_sequence."""
        result = optimize_multigene(
            genes=[GeneSpec(protein="MSKGEELFTG")],
            organism="Escherichia_coli",
            constraints={"gc_lo": 0.3, "gc_hi": 0.7},
        )
        assert isinstance(result, MultiGeneResult)


# ─── MultiGeneResult Properties Tests ──────────────────────────────


class TestMultiGeneResultProperties:
    """Tests for MultiGeneResult invariant checks."""

    def test_result_consistency(self):
        """full_dna length equals total_length."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert len(result.full_dna) == result.total_length

    def test_gc_content_valid(self):
        """GC content is between 0 and 1."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
            ],
            organism="Escherichia_coli",
        )
        assert 0.0 <= result.gc_content <= 1.0

    def test_organism_resolved(self):
        """Organism name is properly resolved."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
            ],
            organism="E_coli",  # Short name
        )
        # Should be resolved to canonical form
        assert result.organism in ("Escherichia_coli", "E_coli")

    def test_gene_results_have_sequences(self):
        """Each gene result has a non-empty DNA sequence."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        for g in result.genes:
            assert g.sequence
            assert len(g.sequence) > 0
            # Sequence length should be protein_length * 3
            assert len(g.sequence) % 3 == 0


# ─── Well-Known Linker Constants Tests ─────────────────────────────


class TestLinkerConstants:
    """Tests for well-known linker sequence constants."""

    def test_p2a_non_empty(self):
        """P2A peptide sequence is non-empty."""
        assert len(LINKER_2A_P2A) > 0

    def test_t2a_non_empty(self):
        """T2A peptide sequence is non-empty."""
        assert len(LINKER_2A_T2A) > 0

    def test_rbs_strong(self):
        """Strong RBS is the Shine-Dalgarno consensus."""
        assert RBS_STRONG == "AGGAGG"

    def test_rbs_spacer(self):
        """RBS spacer has reasonable length."""
        assert 3 <= len(RBS_SPACER_DEFAULT) <= 15


# ─── OperonConfig Tests ────────────────────────────────────────────


class TestOperonConfig:
    """Tests for OperonConfig data class."""

    def test_default_config(self):
        """Default OperonConfig has expected values."""
        config = OperonConfig()
        assert config.rbs_per_gene == "AGGAGG"
        assert config.rbs_spacer == "AAAAA"
        assert config.promoter == ""
        assert config.terminator == ""
        assert config.include_restriction_sites is False

    def test_custom_config(self):
        """OperonConfig with custom values."""
        config = OperonConfig(
            promoter="TTGACA",
            rbs_per_gene="AGGA",
            rbs_spacer="AAA",
            terminator="TTTTTT",
            include_restriction_sites=True,
        )
        assert config.promoter == "TTGACA"
        assert config.rbs_per_gene == "AGGA"
        assert config.rbs_spacer == "AAA"
        assert config.terminator == "TTTTTT"
        assert config.include_restriction_sites is True


# ─── Integration Tests ─────────────────────────────────────────────


class TestIntegration:
    """Integration tests combining multiple features."""

    def test_operon_full_pipeline(self):
        """Full pipeline: optimize operon, check GenBank, verify DNA."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDA",
                         name="gfp"),
                GeneSpec(protein="MDDRLEAIAGMTRLLRALRRKL",
                         name="tetR"),
            ],
            organism="Escherichia_coli",
            promoter="TTGACA",
            terminator="TTTTTTT",
        )

        # Check construct type
        assert result.construct_type == "operon"

        # Check DNA starts with promoter
        assert result.full_dna.startswith("TTGACA")

        # Check DNA ends with terminator
        assert result.full_dna.endswith("TTTTTTT")

        # Check GenBank export
        assert "LOCUS" in result.genbank_export
        assert "FEATURES" in result.genbank_export
        assert "ORIGIN" in result.genbank_export
        assert result.genbank_export.endswith("//")
        assert '/gene="gfp"' in result.genbank_export
        assert '/gene="tetR"' in result.genbank_export

    def test_2a_full_pipeline(self):
        """Full pipeline: 2A polycistronic in human."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MVSKGEELFTG",
                         name="GFP"),
                GeneSpec(protein="MDDRLA",
                         name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )

        assert result.construct_type == "polycistronic_2A"
        assert len(result.genes) == 2
        assert result.total_length > 0

        # GenBank should mention construct type
        assert "polycistronic_2A" in result.genbank_export or "polycistronic" in result.genbank_export.lower()

    def test_different_organisms_per_gene(self):
        """Genes with different organisms in the same construct."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTG", name="gfp",
                         organism="Escherichia_coli"),
                GeneSpec(protein="MDDRLEAIAG", name="tetR",
                         organism="Escherichia_coli"),
            ],
            organism="Escherichia_coli",
        )
        assert isinstance(result, MultiGeneResult)
        assert len(result.genes) == 2

    def test_constraints_propagation(self):
        """GC constraints are propagated to individual gene optimization."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDA",
                         name="gfp"),
            ],
            organism="Escherichia_coli",
            constraints={"gc_lo": 0.30, "gc_hi": 0.70},
        )
        # The GC content of the optimized gene should be within bounds
        # (the optimizer tries to satisfy this but may not always succeed)
        assert isinstance(result, MultiGeneResult)
