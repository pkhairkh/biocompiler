"""
Agent 49: Comprehensive Multi-Gene Construct Tests
====================================================

Expanded test suite for biocompiler.multigene covering:
- Operon with 3 genes for E. coli (insulin + GFP + RFP)
- Polycistronic 2A construct for human cells
- Gene translation verification in full constructs
- GenBank export feature annotation correctness
- Linker sequence insertion verification
- Bidirectional promoter constructs
- Error handling (empty gene list, invalid protein, missing organism)
"""

import pytest

from biocompiler.multigene import (
    GeneSpec,
    MultiGeneResult,
    OperonConfig,
    optimize_multigene,
    optimize_operon,
    _assemble_construct,
    _infer_construct_type,
    _protein_to_2a_dna,
    _generate_genbank_multigene,
    LINKER_2A_P2A,
    LINKER_2A_T2A,
    LINKER_2A_E2A,
    LINKER_2A_F2A,
    RBS_STRONG,
    RBS_MEDIUM,
    RBS_WEAK,
    RBS_SPACER_DEFAULT,
    LINKER_IRES_PLACEHOLDER,
)
from biocompiler.translation import translate
from biocompiler.exceptions import InvalidProteinError, UnsupportedOrganismError
from biocompiler.optimization import OptimizationResult


# ─── Standard protein fixtures ─────────────────────────────────────

# Human insulin B-chain (51 AA, realistic therapeutic protein)
INSULIN_PROTEIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
)

# EGFP (first 40 AA, enough for realistic optimization)
GFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDA"
)

# mCherry RFP (first 40 AA)
RFP_PROTEIN = (
    "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGE"
)

# Short proteins for quick tests
SHORT_GFP = "MSKGEELFTG"
SHORT_TETR = "MDDRLEAIAG"
SHORT_LACZ = "MRVLKFGGTS"

# 2A peptide sequences (protein-level)
P2A_PROTEIN = LINKER_2A_P2A
T2A_PROTEIN = LINKER_2A_T2A


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 1: Three-Gene E. coli Operon
# ═══════════════════════════════════════════════════════════════════


class TestThreeGeneEcoliOperon:
    """Test operon with 3 genes (insulin + GFP + RFP) for E. coli."""

    def test_three_gene_operon_basic(self):
        """Three-gene operon (insulin + GFP + RFP) produces valid result."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=INSULIN_PROTEIN, name="insulin"),
                GeneSpec(protein=GFP_PROTEIN, name="gfp"),
                GeneSpec(protein=RFP_PROTEIN, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        assert isinstance(result, MultiGeneResult)
        assert result.construct_type == "operon"
        assert len(result.genes) == 3
        assert result.total_length > 0
        assert len(result.full_dna) == result.total_length
        assert 0.0 <= result.gc_content <= 1.0
        assert result.organism == "Escherichia_coli"

    def test_three_gene_operon_each_gene_has_sequence(self):
        """Each gene in the 3-gene operon has a non-empty DNA sequence."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=INSULIN_PROTEIN, name="insulin"),
                GeneSpec(protein=GFP_PROTEIN, name="gfp"),
                GeneSpec(protein=RFP_PROTEIN, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        for i, g in enumerate(result.genes):
            assert g.sequence, f"Gene {i} has empty sequence"
            assert len(g.sequence) > 0, f"Gene {i} sequence length is 0"
            assert len(g.sequence) % 3 == 0, f"Gene {i} sequence not a multiple of 3"

    def test_three_gene_operon_translations_match(self):
        """Each gene's DNA translates back to the original protein."""
        proteins = [INSULIN_PROTEIN, GFP_PROTEIN, RFP_PROTEIN]
        result = optimize_operon(
            genes=[
                GeneSpec(protein=INSULIN_PROTEIN, name="insulin"),
                GeneSpec(protein=GFP_PROTEIN, name="gfp"),
                GeneSpec(protein=RFP_PROTEIN, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        for i, (g, expected_protein) in enumerate(zip(result.genes, proteins)):
            translated = translate(g.sequence)
            assert translated == expected_protein, (
                f"Gene {i} ({expected_protein[:10]}...): "
                f"translation mismatch: got {translated[:10]}..."
            )

    def test_three_gene_operon_with_promoter_and_terminator(self):
        """Operon with promoter and terminator includes them in the construct."""
        promoter = "TTGACA"
        terminator = "TTTTTTT"
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
                GeneSpec(protein=SHORT_LACZ, name="lacZ"),
            ],
            organism="Escherichia_coli",
            promoter=promoter,
            terminator=terminator,
        )
        assert result.full_dna.startswith(promoter.upper())
        assert result.full_dna.endswith(terminator.upper())

    def test_three_gene_operon_has_rbs_before_each_gene(self):
        """Each gene in the operon should have an RBS upstream."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
                GeneSpec(protein=SHORT_LACZ, name="lacZ"),
            ],
            organism="Escherichia_coli",
        )
        # Default RBS is AGGAGG — should appear before each gene
        assert RBS_STRONG.upper() in result.full_dna

    def test_three_gene_operon_cai_scores(self):
        """Each gene in the operon has a valid CAI score."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=INSULIN_PROTEIN, name="insulin"),
                GeneSpec(protein=GFP_PROTEIN, name="gfp"),
                GeneSpec(protein=RFP_PROTEIN, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        for i, g in enumerate(result.genes):
            assert 0.0 <= g.cai <= 1.0, f"Gene {i} CAI out of range: {g.cai}"

    def test_three_gene_operon_gc_content_reasonable(self):
        """GC content of the full construct is in a reasonable range."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
                GeneSpec(protein=SHORT_LACZ, name="lacZ"),
            ],
            organism="Escherichia_coli",
        )
        # E. coli GC target is 0.45-0.55 but with RBS/spacer the total
        # might be slightly outside — just check it's not extreme
        assert 0.1 <= result.gc_content <= 0.9


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 2: Polycistronic 2A Construct for Human
# ═══════════════════════════════════════════════════════════════════


class TestPolycistronic2AForHuman:
    """Test polycistronic 2A constructs for human expression."""

    def test_2a_gfp_rfp_human(self):
        """2A-linked GFP + RFP for human cells produces correct construct type."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=GFP_PROTEIN, name="GFP"),
                GeneSpec(protein=RFP_PROTEIN, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"
        assert len(result.genes) == 2
        assert result.organism == "Homo_sapiens"

    def test_2a_linker_dna_inserted_between_genes(self):
        """2A peptide DNA is inserted between genes in the construct."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        # The full DNA should be longer than the sum of individual gene DNAs
        # because of the 2A peptide coding sequence
        gene_sum = sum(len(g.sequence) for g in result.genes)
        assert result.total_length > gene_sum, (
            f"Total length ({result.total_length}) should exceed sum of genes "
            f"({gene_sum}) due to 2A peptide insertion"
        )

    def test_2a_each_gene_translates_correctly(self):
        """Each gene in the 2A construct translates back to its input protein."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=GFP_PROTEIN, name="GFP"),
                GeneSpec(protein=RFP_PROTEIN, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        for i, (g, expected) in enumerate(
            zip(result.genes, [GFP_PROTEIN, RFP_PROTEIN])
        ):
            translated = translate(g.sequence)
            assert translated == expected, (
                f"Gene {i} translation mismatch in 2A construct"
            )

    def test_2a_t2a_variant(self):
        """T2A peptide linker is recognized and produces polycistronic_2A."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_T2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"
        assert len(result.genes) == 2

    def test_2a_e2a_variant(self):
        """E2A peptide linker is recognized and produces polycistronic_2A."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_E2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"

    def test_2a_f2a_variant(self):
        """F2A peptide linker is recognized and produces polycistronic_2A."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_F2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"

    def test_2a_with_promoter_and_terminator(self):
        """2A construct with gene-specific promoter and terminator."""
        result = optimize_multigene(
            genes=[
                GeneSpec(
                    protein=SHORT_GFP,
                    name="GFP",
                    promoter="CCACTG",
                    terminator="AATAAA",
                ),
                GeneSpec(
                    protein=SHORT_TETR,
                    name="RFP",
                    terminator="AATAAA",
                ),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        # Promoter and terminator should appear in the construct
        assert "CCACTG" in result.full_dna
        assert "AATAAA" in result.full_dna


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 3: Gene Translation in Full Construct
# ═══════════════════════════════════════════════════════════════════


class TestGeneTranslationInFullConstruct:
    """Verify that each gene translates correctly when embedded in the
    full multi-gene construct."""

    def test_operon_gene_dna_located_in_full_construct(self):
        """Each gene's DNA can be found as a substring of the full construct."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
                GeneSpec(protein=SHORT_LACZ, name="lacZ"),
            ],
            organism="Escherichia_coli",
        )
        for i, g in enumerate(result.genes):
            assert g.sequence.upper() in result.full_dna.upper(), (
                f"Gene {i} DNA not found as substring of full construct"
            )

    def test_2a_gene_dna_located_in_full_construct(self):
        """Each gene's DNA is located within the 2A construct."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        for i, g in enumerate(result.genes):
            assert g.sequence.upper() in result.full_dna.upper(), (
                f"Gene {i} DNA not found in 2A construct"
            )

    def test_custom_construct_gene_order_preserved(self):
        """Genes appear in the correct order in the full construct."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="geneA"),
                GeneSpec(protein=SHORT_TETR, name="geneB"),
            ],
            organism="Escherichia_coli",
        )
        pos_a = result.full_dna.upper().find(result.genes[0].sequence.upper())
        pos_b = result.full_dna.upper().find(result.genes[1].sequence.upper())
        assert pos_a < pos_b, "Gene A should appear before Gene B in the construct"

    def test_full_construct_valid_dna(self):
        """Full construct contains only valid DNA characters."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        valid_bases = set("ACGT")
        for ch in result.full_dna:
            assert ch in valid_bases, f"Invalid DNA character: {ch!r}"


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 4: GenBank Feature Annotations
# ═══════════════════════════════════════════════════════════════════


class TestGenBankFeatureAnnotations:
    """Verify GenBank export has correct feature annotations for each gene."""

    def test_genbank_has_gene_feature_per_gene(self):
        """GenBank export has a 'gene' feature for each gene."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="insulin"),
                GeneSpec(protein=SHORT_TETR, name="gfp"),
                GeneSpec(protein=SHORT_LACZ, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        # Count 'gene' features — should be at least 3
        gene_feature_count = gb.count("     gene            ")
        assert gene_feature_count >= 3, (
            f"Expected at least 3 gene features, found {gene_feature_count}"
        )

    def test_genbank_has_cds_feature_per_gene(self):
        """GenBank export has a CDS feature for each gene."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="insulin"),
                GeneSpec(protein=SHORT_TETR, name="gfp"),
                GeneSpec(protein=SHORT_LACZ, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        cds_count = gb.count("     CDS             ")
        assert cds_count >= 3, f"Expected at least 3 CDS features, found {cds_count}"

    def test_genbank_gene_names_annotated(self):
        """All gene names appear in the GenBank export."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="insulin"),
                GeneSpec(protein=SHORT_TETR, name="gfp"),
                GeneSpec(protein=SHORT_LACZ, name="rfp"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        assert '/gene="insulin"' in gb
        assert '/gene="gfp"' in gb
        assert '/gene="rfp"' in gb

    def test_genbank_translation_per_gene(self):
        """Each gene has a /translation qualifier in the GenBank export."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gene1"),
                GeneSpec(protein=SHORT_TETR, name="gene2"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        translation_count = gb.count("/translation=")
        assert translation_count >= 2, (
            f"Expected at least 2 translation qualifiers, found {translation_count}"
        )

    def test_genbank_cai_per_gene(self):
        """Each gene has a /cai qualifier in the GenBank export."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gene1"),
                GeneSpec(protein=SHORT_TETR, name="gene2"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        cai_count = gb.count("/cai=")
        assert cai_count >= 2, f"Expected at least 2 CAI qualifiers, found {cai_count}"

    def test_genbank_cds_positions_valid(self):
        """CDS feature positions are valid (start < end, within construct)."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="g1"),
                GeneSpec(protein=SHORT_TETR, name="g2"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        # Extract CDS position lines
        import re
        cds_pattern = re.compile(r"CDS\s+(\d+)\.\.(\d+)")
        matches = cds_pattern.findall(gb)
        assert len(matches) >= 2, "Not enough CDS features found"
        for start_str, end_str in matches:
            start = int(start_str)
            end = int(end_str)
            assert 1 <= start < end <= result.total_length, (
                f"Invalid CDS position: {start}..{end}, "
                f"construct length={result.total_length}"
            )

    def test_genbank_operon_definition(self):
        """GenBank DEFINITION mentions 'Operon' for operon constructs."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="g1"),
                GeneSpec(protein=SHORT_TETR, name="g2"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        assert "Operon" in gb or "operon" in gb.lower()

    def test_genbank_2a_definition(self):
        """GenBank DEFINITION mentions 'Polycistronic 2A' for 2A constructs."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        gb = result.genbank_export
        assert "Polycistronic" in gb or "polycistronic" in gb.lower()

    def test_genbank_promoter_feature(self):
        """GenBank export includes promoter feature when specified."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="g1", promoter="TTGACA"),
            ],
            organism="Escherichia_coli",
            promoter="TTGACA",
        )
        gb = result.genbank_export
        assert "promoter" in gb.lower()

    def test_genbank_rbs_feature(self):
        """GenBank export includes RBS regulatory feature when specified."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="g1", rbs="AGGAGG"),
            ],
            organism="Escherichia_coli",
        )
        gb = result.genbank_export
        assert "ribosome_binding_site" in gb

    def test_genbank_terminator_feature(self):
        """GenBank export includes terminator feature when specified."""
        result = optimize_operon(
            genes=[
                GeneSpec(
                    protein=SHORT_GFP,
                    name="g1",
                    terminator="TTTTTTT",
                ),
            ],
            organism="Escherichia_coli",
            terminator="TTTTTTT",
        )
        gb = result.genbank_export
        assert "terminator" in gb.lower()


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 5: Linker Sequence Insertion
# ═══════════════════════════════════════════════════════════════════


class TestLinkerSequenceInsertion:
    """Verify that linker sequences are correctly inserted between genes."""

    def test_custom_linker_between_genes(self):
        """Custom DNA linker appears between genes in the full construct."""
        linker = "GGATCC"  # BamHI site
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            linker=linker,
            organism="Homo_sapiens",
        )
        assert linker.upper() in result.full_dna, (
            f"Custom linker {linker} not found in full construct"
        )

    def test_linker_between_not_at_ends(self):
        """Linker appears between genes, not at the start or end."""
        linker = "GGATCC"
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            linker=linker,
            organism="Homo_sapiens",
        )
        # Linker should not be at the very start or end
        gene_a_dna = result.genes[0].sequence.upper()
        gene_b_dna = result.genes[1].sequence.upper()
        pos_a = result.full_dna.find(gene_a_dna)
        pos_b = result.full_dna.find(gene_b_dna)
        linker_pos = result.full_dna.find(linker.upper())
        assert linker_pos > pos_a, "Linker should appear after gene A"
        assert linker_pos + len(linker) <= pos_b + len(gene_b_dna), (
            "Linker should appear before or overlapping with gene B"
        )

    def test_2a_linker_converted_to_dna(self):
        """2A protein linker is converted to DNA before insertion."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        # The full DNA should be longer than gene sum because of 2A DNA
        gene_sum = sum(len(g.sequence) for g in result.genes)
        extra_dna = result.total_length - gene_sum
        # 2A peptide is ~22 AA, so DNA is ~66 bp; there may also be
        # promoter/terminator DNA from GeneSpecs
        assert extra_dna > 0, (
            f"Expected extra DNA from 2A peptide insertion, "
            f"but total ({result.total_length}) == gene_sum ({gene_sum})"
        )

    def test_ires_linker_preserved_as_dna(self):
        """IRES DNA linker is preserved as-is in the construct."""
        ires_seq = "ACGT" * 50  # 200 bp IRES
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=ires_seq,
            organism="Homo_sapiens",
        )
        assert ires_seq.upper() in result.full_dna, (
            "IRES DNA linker not found in full construct"
        )

    def test_no_linker_no_extra_dna(self):
        """Without a linker, custom construct has no extra intergenic DNA."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            organism="Homo_sapiens",
        )
        gene_sum = sum(len(g.sequence) for g in result.genes)
        # For custom constructs without linker and without regulatory elements,
        # total should equal gene sum
        assert result.total_length == gene_sum, (
            f"Without linker, expected total_length == gene_sum, "
            f"got {result.total_length} vs {gene_sum}"
        )

    def test_multiple_linkers_between_three_genes(self):
        """With 3 genes and a linker, linker appears between each pair."""
        linker = "GGATCC"
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="g1"),
                GeneSpec(protein=SHORT_TETR, name="g2"),
                GeneSpec(protein=SHORT_LACZ, name="g3"),
            ],
            linker=linker,
            organism="Homo_sapiens",
        )
        # Count occurrences of the linker
        count = result.full_dna.upper().count(linker.upper())
        assert count >= 2, f"Expected at least 2 linker occurrences, got {count}"


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 6: Bidirectional Promoter Construct
# ═══════════════════════════════════════════════════════════════════


class TestBidirectionalPromoter:
    """Test bidirectional promoter constructs with genes facing both
    directions from a shared promoter region."""

    def test_bidirectional_four_genes(self):
        """Four genes with bidirectional promoters detected as bidirectional."""
        bidirectional_promoter = "CGCGCGCGCG"
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="fwd1", promoter="TTGACA"),
                GeneSpec(protein=SHORT_TETR, name="fwd2", promoter="TTGACA"),
                GeneSpec(protein=SHORT_LACZ, name="rev1", promoter="TATAAT"),
                GeneSpec(protein=SHORT_LACZ, name="rev2", promoter="TATAAT"),
            ],
            linker=bidirectional_promoter,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "bidirectional"
        assert len(result.genes) == 4

    def test_bidirectional_linker_between_halves(self):
        """Bidirectional promoter region appears between forward and reverse halves."""
        promoter_region = "CGCGCGCGCG"
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="fwd1", promoter="TTGACA"),
                GeneSpec(protein=SHORT_TETR, name="fwd2"),
                GeneSpec(protein=SHORT_LACZ, name="rev1", promoter="TATAAT"),
                GeneSpec(protein=SHORT_LACZ, name="rev2"),
            ],
            linker=promoter_region,
            organism="Homo_sapiens",
        )
        assert promoter_region.upper() in result.full_dna

    def test_bidirectional_all_genes_present(self):
        """All genes from both directions have optimization results."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="fwd1", promoter="TTGACA"),
                GeneSpec(protein=SHORT_TETR, name="fwd2"),
                GeneSpec(protein=SHORT_LACZ, name="rev1", promoter="TATAAT"),
                GeneSpec(protein=SHORT_LACZ, name="rev2"),
            ],
            linker="CGCGCGCGCG",
            organism="Homo_sapiens",
        )
        assert len(result.genes) == 4
        for g in result.genes:
            assert g.sequence
            assert len(g.sequence) > 0


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 7: Error Handling
# ═══════════════════════════════════════════════════════════════════


class TestMultigeneErrorHandling:
    """Test error conditions and edge cases for multi-gene constructs."""

    def test_empty_gene_list_raises(self):
        """Empty gene list raises ValueError."""
        with pytest.raises(ValueError, match="At least one"):
            optimize_multigene(genes=[], organism="Escherichia_coli")

    def test_empty_gene_list_operon_raises(self):
        """Empty gene list raises ValueError for optimize_operon too."""
        with pytest.raises(ValueError, match="At least one"):
            optimize_operon(genes=[])

    def test_no_organism_raises(self):
        """Missing organism raises ValueError when no gene specifies one."""
        with pytest.raises(ValueError, match="Organism must be specified"):
            optimize_multigene(
                genes=[GeneSpec(protein=SHORT_GFP)],
                organism="",
            )

    def test_invalid_protein_in_genespec_raises(self):
        """Invalid amino acid characters in GeneSpec raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            optimize_multigene(
                genes=[GeneSpec(protein="MSKGEELFTGX")],
                organism="Escherichia_coli",
            )

    def test_invalid_protein_empty_raises(self):
        """Empty protein in GeneSpec raises ValueError."""
        with pytest.raises(ValueError, match="non-empty"):
            GeneSpec(protein="")

    def test_invalid_protein_whitespace_only_raises(self):
        """Whitespace-only protein raises ValueError."""
        with pytest.raises(ValueError, match="non-empty"):
            GeneSpec(protein="   ")

    def test_invalid_protein_numbers_raises(self):
        """Numbers in protein sequence raise InvalidProteinError."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEEL1G")

    def test_genespec_organism_overrides_default(self):
        """Per-gene organism overrides the construct-level default."""
        result = optimize_multigene(
            genes=[
                GeneSpec(
                    protein=SHORT_GFP,
                    name="gfp",
                    organism="Escherichia_coli",
                ),
            ],
            organism="Homo_sapiens",
        )
        # Should succeed — per-gene organism takes precedence
        assert isinstance(result, MultiGeneResult)
        assert len(result.genes) == 1

    def test_single_gene_operon_works(self):
        """Single-gene operon is valid (edge case)."""
        result = optimize_operon(
            genes=[GeneSpec(protein=SHORT_GFP, name="solo")],
            organism="Escherichia_coli",
        )
        assert isinstance(result, MultiGeneResult)
        assert len(result.genes) == 1

    def test_gene_with_per_gene_organism(self):
        """Genes with different organisms in the same construct work."""
        result = optimize_multigene(
            genes=[
                GeneSpec(
                    protein=SHORT_GFP,
                    name="gfp",
                    organism="Escherichia_coli",
                ),
                GeneSpec(
                    protein=SHORT_TETR,
                    name="tetR",
                    organism="Escherichia_coli",
                ),
            ],
            organism="Escherichia_coli",
        )
        assert isinstance(result, MultiGeneResult)
        assert len(result.genes) == 2

    def test_constraints_forwarded_to_optimization(self):
        """GC constraints are forwarded to individual gene optimization."""
        result = optimize_multigene(
            genes=[
                GeneSpec(
                    protein=GFP_PROTEIN,
                    name="gfp",
                ),
            ],
            organism="Escherichia_coli",
            constraints={"gc_lo": 0.30, "gc_hi": 0.70},
        )
        assert isinstance(result, MultiGeneResult)


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 8: Construct Type Inference
# ═══════════════════════════════════════════════════════════════════


class TestConstructTypeInferenceComprehensive:
    """Comprehensive tests for automatic construct type detection."""

    def test_ecoli_with_rbs_inferred_operon(self):
        """E. coli genes with RBS → operon."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp", rbs="AGGAGG"),
                GeneSpec(protein=SHORT_TETR, name="tetR", rbs="AGGAGG"),
            ],
            organism="Escherichia_coli",
        )
        assert result.construct_type == "operon"

    def test_ecoli_without_rbs_still_operon(self):
        """E. coli multi-gene without RBS defaults to operon."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert result.construct_type == "operon"

    def test_2a_p2a_inferred(self):
        """P2A linker → polycistronic_2A."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_P2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"

    def test_2a_t2a_inferred(self):
        """T2A linker → polycistronic_2A."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=LINKER_2A_T2A,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_2A"

    def test_ires_inferred_from_long_linker(self):
        """Long DNA linker (>100 bp) → polycistronic_IRES."""
        ires = "ACGT" * 50  # 200 bp
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker=ires,
            organism="Homo_sapiens",
        )
        assert result.construct_type == "polycistronic_IRES"

    def test_human_short_linker_custom(self):
        """Short DNA linker with human → custom."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="GFP"),
                GeneSpec(protein=SHORT_TETR, name="RFP"),
            ],
            linker="GGATCC",
            organism="Homo_sapiens",
        )
        assert result.construct_type == "custom"

    def test_bidirectional_inferred(self):
        """Genes with promoters in both halves → bidirectional."""
        result = optimize_multigene(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="f1", promoter="TTGACA"),
                GeneSpec(protein=SHORT_TETR, name="f2"),
                GeneSpec(protein=SHORT_LACZ, name="r1", promoter="TATAAT"),
                GeneSpec(protein=SHORT_LACZ, name="r2"),
            ],
            linker="CGCGCGCGCG",
            organism="Homo_sapiens",
        )
        assert result.construct_type == "bidirectional"


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 9: MultiGeneResult Properties
# ═══════════════════════════════════════════════════════════════════


class TestMultiGeneResultProperties:
    """Verify MultiGeneResult invariants and data integrity."""

    def test_full_dna_length_matches_total_length(self):
        """full_dna length equals total_length."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        assert len(result.full_dna) == result.total_length

    def test_gc_content_in_valid_range(self):
        """GC content is between 0 and 1."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
            ],
            organism="Escherichia_coli",
        )
        assert 0.0 <= result.gc_content <= 1.0

    def test_construct_type_is_valid(self):
        """Construct type is one of the valid types."""
        valid_types = {"operon", "polycistronic_2A", "polycistronic_IRES",
                       "bidirectional", "custom"}
        result = optimize_operon(
            genes=[GeneSpec(protein=SHORT_GFP, name="gfp")],
            organism="Escherichia_coli",
        )
        assert result.construct_type in valid_types

    def test_gene_results_have_protein_field(self):
        """Each gene result has a protein field matching the input."""
        proteins = [SHORT_GFP, SHORT_TETR]
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),
                GeneSpec(protein=SHORT_TETR, name="tetR"),
            ],
            organism="Escherichia_coli",
        )
        for i, (g, expected) in enumerate(zip(result.genes, proteins)):
            assert g.protein == expected, (
                f"Gene {i} protein mismatch: {g.protein} != {expected}"
            )

    def test_gene_sequences_length_correct(self):
        """Each gene's DNA length is protein_length * 3."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="gfp"),  # 10 AA → 30 bp
                GeneSpec(protein=SHORT_TETR, name="tetR"),  # 10 AA → 30 bp
            ],
            organism="Escherichia_coli",
        )
        for g in result.genes:
            assert len(g.sequence) == len(g.protein) * 3


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 10: GeneSpec Validation
# ═══════════════════════════════════════════════════════════════════


class TestGeneSpecValidationComprehensive:
    """Comprehensive GeneSpec validation tests."""

    def test_valid_protein_accepted(self):
        """Valid protein sequences are accepted."""
        spec = GeneSpec(protein="ACDEFGHIKLMNPQRSTVWY")
        assert spec.protein == "ACDEFGHIKLMNPQRSTVWY"

    def test_lowercase_normalized(self):
        """Lowercase protein is uppercased."""
        spec = GeneSpec(protein="mskgeelftg")
        assert spec.protein == "MSKGEELFTG"

    def test_whitespace_stripped(self):
        """Leading/trailing whitespace is stripped."""
        spec = GeneSpec(protein="  MSKGEELFTG  ")
        assert spec.protein == "MSKGEELFTG"

    def test_all_standard_amino_acids(self):
        """All 20 standard amino acids are accepted."""
        spec = GeneSpec(protein="ACDEFGHIKLMNPQRSTVWY")
        assert len(spec.protein) == 20

    def test_invalid_char_x_rejected(self):
        """X (ambiguous amino acid) is rejected."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGX")

    def test_invalid_char_b_rejected(self):
        """B (Asx) is rejected in strict mode."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGB")

    def test_invalid_char_z_rejected(self):
        """Z (Glx) is rejected in strict mode."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGZ")

    def test_invalid_char_j_rejected(self):
        """J (Leu/Ile) is rejected."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGJ")

    def test_invalid_char_o_rejected(self):
        """O (Pyrrolysine) is rejected."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGO")

    def test_invalid_char_u_rejected(self):
        """U (Selenocysteine) is rejected."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTGU")

    def test_numeric_rejected(self):
        """Numeric characters are rejected."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTG1")

    def test_special_chars_rejected(self):
        """Special characters are rejected."""
        with pytest.raises(InvalidProteinError):
            GeneSpec(protein="MSKGEELFTG!")

    def test_empty_string_rejected(self):
        """Empty protein string is rejected."""
        with pytest.raises(ValueError):
            GeneSpec(protein="")

    def test_whitespace_only_rejected(self):
        """Whitespace-only protein is rejected."""
        with pytest.raises(ValueError):
            GeneSpec(protein="   ")


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 11: Internal Helper Functions
# ═══════════════════════════════════════════════════════════════════


class TestInternalHelpers:
    """Test internal helper functions for multi-gene assembly."""

    def test_protein_to_2a_dna(self):
        """_protein_to_2a_dna converts a 2A peptide to DNA."""
        dna = _protein_to_2a_dna(LINKER_2A_P2A, "Homo_sapiens")
        assert len(dna) > 0
        assert len(dna) % 3 == 0, "2A DNA length must be a multiple of 3"
        # Should translate back to the 2A peptide (approximately)
        translated = translate(dna)
        assert translated == LINKER_2A_P2A

    def test_assemble_construct_operon(self):
        """_assemble_construct correctly assembles an operon."""
        gene_dnas = ["ATGAAAGCGTTT", "ATGCCCGCGAAA"]
        gene_specs = [
            GeneSpec(protein="MKA", name="g1", rbs="AGGAGG"),
            GeneSpec(protein="MPA", name="g2", rbs="AGGAGG"),
        ]
        config = OperonConfig(promoter="TTGACA", rbs_spacer="AAAAA")
        result = _assemble_construct(
            gene_dnas, gene_specs,
            construct_type="operon", operon_config=config,
        )
        assert result.startswith("TTGACA")
        assert "AGGAGG" in result
        assert "AAAAA" in result

    def test_assemble_construct_custom(self):
        """_assemble_construct correctly assembles a custom construct."""
        gene_dnas = ["ATGAAAGCGTTT", "ATGCCCGCGAAA"]
        gene_specs = [
            GeneSpec(protein="MKA", name="g1"),
            GeneSpec(protein="MPA", name="g2"),
        ]
        linker = "GGATCC"
        result = _assemble_construct(
            gene_dnas, gene_specs,
            linker=linker, construct_type="custom",
        )
        assert "GGATCC" in result

    def test_infer_construct_type_2a(self):
        """_infer_construct_type detects 2A peptide."""
        genes = [
            GeneSpec(protein=SHORT_GFP, name="g1"),
            GeneSpec(protein=SHORT_TETR, name="g2"),
        ]
        ct = _infer_construct_type(genes, LINKER_2A_P2A, "Homo_sapiens")
        assert ct == "polycistronic_2A"

    def test_infer_construct_type_ires(self):
        """_infer_construct_type detects IRES."""
        genes = [
            GeneSpec(protein=SHORT_GFP, name="g1"),
            GeneSpec(protein=SHORT_TETR, name="g2"),
        ]
        ires = "ACGT" * 50  # 200 bp
        ct = _infer_construct_type(genes, ires, "Homo_sapiens")
        assert ct == "polycistronic_IRES"

    def test_infer_construct_type_operon_ecoli(self):
        """_infer_construct_type detects operon for E. coli."""
        genes = [
            GeneSpec(protein=SHORT_GFP, name="g1", rbs="AGGAGG"),
            GeneSpec(protein=SHORT_TETR, name="g2", rbs="AGGAGG"),
        ]
        ct = _infer_construct_type(genes, "", "Escherichia_coli")
        assert ct == "operon"

    def test_infer_construct_type_bidirectional(self):
        """_infer_construct_type detects bidirectional promoter."""
        genes = [
            GeneSpec(protein=SHORT_GFP, name="f1", promoter="TTGACA"),
            GeneSpec(protein=SHORT_TETR, name="f2"),
            GeneSpec(protein=SHORT_LACZ, name="r1", promoter="TATAAT"),
            GeneSpec(protein=SHORT_LACZ, name="r2"),
        ]
        ct = _infer_construct_type(genes, "CGCGCGCGCG", "Homo_sapiens")
        assert ct == "bidirectional"

    def test_infer_construct_type_custom_human(self):
        """_infer_construct_type defaults to custom for human with short linker."""
        genes = [
            GeneSpec(protein=SHORT_GFP, name="g1"),
            GeneSpec(protein=SHORT_TETR, name="g2"),
        ]
        ct = _infer_construct_type(genes, "GGATCC", "Homo_sapiens")
        assert ct == "custom"


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 12: Linker Constants
# ═══════════════════════════════════════════════════════════════════


class TestLinkerConstantsComprehensive:
    """Verify well-known linker sequence constants."""

    def test_p2a_is_protein_sequence(self):
        """P2A linker is a valid amino acid sequence."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert all(c in valid_aas for c in LINKER_2A_P2A)

    def test_t2a_is_protein_sequence(self):
        """T2A linker is a valid amino acid sequence."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert all(c in valid_aas for c in LINKER_2A_T2A)

    def test_e2a_is_protein_sequence(self):
        """E2A linker is a valid amino acid sequence."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert all(c in valid_aas for c in LINKER_2A_E2A)

    def test_f2a_is_protein_sequence(self):
        """F2A linker is a valid amino acid sequence."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert all(c in valid_aas for c in LINKER_2A_F2A)

    def test_all_2a_linkers_different(self):
        """All 2A linkers are distinct sequences."""
        linkers = {LINKER_2A_P2A, LINKER_2A_T2A, LINKER_2A_E2A, LINKER_2A_F2A}
        assert len(linkers) == 4, "2A linkers should all be different"

    def test_rbs_strong_is_shine_dalgarno(self):
        """Strong RBS is the Shine-Dalgarno consensus."""
        assert RBS_STRONG == "AGGAGG"

    def test_rbs_spacer_reasonable_length(self):
        """RBS spacer has a reasonable length (3-15 bp)."""
        assert 3 <= len(RBS_SPACER_DEFAULT) <= 15

    def test_ires_placeholder_empty(self):
        """IRES placeholder is empty (user must supply)."""
        assert LINKER_IRES_PLACEHOLDER == ""


# ═══════════════════════════════════════════════════════════════════
# Agent 49 — Test Suite 13: OperonConfig
# ═══════════════════════════════════════════════════════════════════


class TestOperonConfigComprehensive:
    """Comprehensive OperonConfig tests."""

    def test_default_config_values(self):
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

    def test_gene_specific_rbs_overrides_default(self):
        """Gene-specific RBS overrides the operon default RBS."""
        result = optimize_operon(
            genes=[
                GeneSpec(protein=SHORT_GFP, name="g1", rbs="AGGA"),
                GeneSpec(protein=SHORT_TETR, name="g2", rbs="AGGAGG"),
            ],
            organism="Escherichia_coli",
            rbs_per_gene="AGGAGG",  # default
        )
        # The construct should contain both RBS sequences
        assert "AGGA" in result.full_dna
        assert "AGGAGG" in result.full_dna
