"""Tests for BioCompiler IR codegen — GenBank / FASTA backend output.

Verifies that:

* :func:`to_genbank` produces well-formed GenBank records (LOCUS,
  FEATURES, ORIGIN, ``//`` terminator).
* :func:`to_fasta_dna` / :func:`to_fasta_rna` / :func:`to_fasta_protein`
  produce FASTA records with the correct ``>gene|organism|level`` header
  and 60-char line wrapping.
* The HBB gene (human hemoglobin beta, N-terminal 31 aa + stop) can be
  compiled L0 → L3 and the protein FASTA matches the canonical
  ``MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*``.
* Round-trip: ``to_genbank()`` → ``parse_genbank_sequence()`` recovers
  the original DNA sequence byte-for-byte.
* Eukaryotic genes with exon / intron / UTR regions produce correct
  ``join()`` CDS features and per-exon annotations.
"""

from __future__ import annotations

import sys

import pytest

sys.path.insert(0, "src")

from biocompiler.ir import (
    IR_L0_GenomicDNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IRLevel,
    GeneRegion,
    compile_gene,
)
from biocompiler.ir.codegen import (
    to_genbank,
    to_fasta_dna,
    to_fasta_rna,
    to_fasta_protein,
    to_fasta,
    to_sbol3,
    parse_genbank_sequence,
    parse_genbank_features,
    FASTA_LINE_WIDTH,
)


# ────────────────────────────────────────────────────────────────────
# Fixtures: HBB gene (N-terminal 31 aa + stop, the canonical fragment
# used in the IR pipeline tests).
# ────────────────────────────────────────────────────────────────────

HBB_DNA = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG"
    "GTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
    "GTAA"  # AGG=R completes codon 31, then TAA=*
)
HBB_PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"


@pytest.fixture
def hbb_l0() -> IR_L0_GenomicDNA:
    return IR_L0_GenomicDNA(
        sequence=HBB_DNA, regions=[], organism="Homo_sapiens", gene_name="HBB"
    )


@pytest.fixture
def hbb_l3() -> IR_L3_Polypeptide:
    return compile_gene(IR_L0_GenomicDNA(
        sequence=HBB_DNA, regions=[], organism="Homo_sapiens", gene_name="HBB"
    ), IRLevel.L3)


@pytest.fixture
def simple_l0() -> IR_L0_GenomicDNA:
    """A minimal prokaryotic gene: ATG GCT TAA → MAK*."""
    return IR_L0_GenomicDNA(
        sequence="ATGGCTTAA", regions=[],
        organism="e_coli", gene_name="simple",
    )


@pytest.fixture
def euk_l0() -> IR_L0_GenomicDNA:
    """A eukaryotic-style gene with 5'UTR + 2 exons + intron + 3'UTR."""
    # 5'UTR=GGG, exon1=ATGGCT, intron=TTTTTT, exon2=AAGTAA, 3'UTR=CCC
    seq = "GGG" + "ATGGCT" + "TTTTTT" + "AAGTAA" + "CCC"
    regions = [
        GeneRegion(0, 3, "5_utr"),
        GeneRegion(3, 9, "exon"),
        GeneRegion(9, 15, "intron"),
        GeneRegion(15, 21, "exon"),
        GeneRegion(21, 24, "3_utr"),
    ]
    return IR_L0_GenomicDNA(
        sequence=seq, regions=regions,
        organism="human", gene_name="testEuk",
    )


# ════════════════════════════════════════════════════════════════════
# GenBank format — structural validity
# ════════════════════════════════════════════════════════════════════
class TestGenBankStructure:
    """A valid GenBank record must have LOCUS, FEATURES, ORIGIN, //."""

    def test_starts_with_LOCUS(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert gb.startswith("LOCUS"), f"GenBank must start with LOCUS: {gb[:50]!r}"

    def test_ends_with_double_slash(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert gb.rstrip().endswith("//"), "GenBank must end with //"

    def test_has_FEATURES_section(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nFEATURES" in gb, "missing FEATURES section"

    def test_has_ORIGIN_section(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nORIGIN" in gb, "missing ORIGIN section"

    def test_has_DEFINITION(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nDEFINITION" in gb

    def test_has_ACCESSION(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nACCESSION" in gb
        # Accession should be BC-prefixed.
        for line in gb.splitlines():
            if line.startswith("ACCESSION"):
                acc = line.split()[1]
                assert acc.startswith("BC"), f"accession {acc!r} should start with BC"
                break

    def test_has_VERSION(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nVERSION" in gb

    def test_has_SOURCE_and_ORGANISM(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nSOURCE" in gb
        assert "ORGANISM" in gb

    def test_has_KEYWORDS(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "\nKEYWORDS" in gb

    def test_LOCUS_contains_length_and_DNA(self, simple_l0):
        gb = to_genbank(simple_l0)
        locus_line = gb.splitlines()[0]
        assert "9 bp" in locus_line
        assert "DNA" in locus_line
        assert "linear" in locus_line

    def test_LOCUS_contains_gene_name(self, simple_l0):
        gb = to_genbank(simple_l0)
        locus_line = gb.splitlines()[0]
        assert "simple" in locus_line


# ════════════════════════════════════════════════════════════════════
# GenBank FEATURES table
# ════════════════════════════════════════════════════════════════════
class TestGenBankFeatures:
    def test_prokaryotic_has_gene_and_CDS(self, simple_l0):
        gb = to_genbank(simple_l0)
        features = parse_genbank_features(gb)
        keys = [k for k, _ in features]
        assert "gene" in keys
        assert "CDS" in keys

    def test_prokaryotic_CDS_spans_whole_sequence(self, simple_l0):
        gb = to_genbank(simple_l0)
        features = dict(parse_genbank_features(gb))
        # 9 bp → 1..9
        assert features["CDS"] == "1..9"
        assert features["gene"] == "1..9"

    def test_eukaryotic_CDS_uses_join(self, euk_l0):
        gb = to_genbank(euk_l0)
        features = dict(parse_genbank_features(gb))
        # Exon1 = [3,9) → 4..9, exon2 = [15,21) → 16..21
        assert features["CDS"] == "join(4..9,16..21)"

    def test_eukaryotic_has_two_exon_features(self, euk_l0):
        gb = to_genbank(euk_l0)
        features = parse_genbank_features(gb)
        exons = [(k, v) for k, v in features if k == "exon"]
        assert len(exons) == 2
        assert exons[0] == ("exon", "4..9")
        assert exons[1] == ("exon", "16..21")

    def test_eukaryotic_has_UTR_features(self, euk_l0):
        gb = to_genbank(euk_l0)
        features = dict(parse_genbank_features(gb))
        # 5'UTR [0,3) → 1..3, 3'UTR [21,24) → 22..24
        assert features.get("5'UTR") == "1..3"
        assert features.get("3'UTR") == "22..24"

    def test_eukaryotic_has_intron_feature(self, euk_l0):
        gb = to_genbank(euk_l0)
        features = dict(parse_genbank_features(gb))
        assert features.get("intron") == "10..15"

    def test_CDS_has_codon_start_and_transl_table(self, simple_l0):
        gb = to_genbank(simple_l0)
        assert "/codon_start=1" in gb
        assert "/transl_table=1" in gb

    def test_features_carry_gene_qualifier(self, simple_l0):
        gb = to_genbank(simple_l0)
        # Should have /gene="simple" qualifier lines
        assert '/gene="simple"' in gb

    def test_features_carry_organism_qualifier(self, simple_l0):
        gb = to_genbank(simple_l0)
        # e_coli → "e coli" (underscore → space)
        assert '/organism="e coli"' in gb


# ════════════════════════════════════════════════════════════════════
# GenBank ORIGIN sequence section
# ════════════════════════════════════════════════════════════════════
class TestGenBankOrigin:
    def test_origin_has_sequence(self, simple_l0):
        gb = to_genbank(simple_l0)
        seq = parse_genbank_sequence(gb)
        assert seq == "ATGGCTTAA"

    def test_sequence_is_uppercase(self, simple_l0):
        ir = IR_L0_GenomicDNA(
            sequence="atggcttaa", regions=[],
            organism="e_coli", gene_name="lc",
        )
        gb = to_genbank(ir)
        seq = parse_genbank_sequence(gb)
        assert seq == "ATGGCTTAA"

    def test_sequence_lines_max_60_chars(self):
        # 120 bp → 2 lines of 60 each.
        seq = "ATG" + "GCT" * 39 + "TAA"  # 3 + 117 + 3 = 123 bp
        ir = IR_L0_GenomicDNA(
            sequence=seq, regions=[],
            organism="e_coli", gene_name="long",
        )
        gb = to_genbank(ir)
        in_origin = False
        seq_line_lengths = []
        for line in gb.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if not in_origin:
                continue
            if line.startswith("//"):
                break
            # Each data line: 9-char number, space, then groups of 10.
            # The "bases portion" is everything after the leading number.
            bases_part = line.split(None, 1)[1] if " " in line else ""
            # Without spaces, bases count should be ≤ 60.
            base_count = len(bases_part.replace(" ", ""))
            seq_line_lengths.append(base_count)
        assert all(n <= 60 for n in seq_line_lengths), seq_line_lengths

    def test_sequence_groups_of_10(self, simple_l0):
        gb = to_genbank(simple_l0)
        # Find an ORIGIN data line and verify it's grouped in 10s.
        in_origin = False
        for line in gb.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if not in_origin or line.startswith("//"):
                break
            # First data line for 9-bp seq: "        1 ATGGCTTAA"
            # Single group, no spaces within (only 9 bases).
            assert "ATGGCTTAA" in line.replace(" ", "")
            break


# ════════════════════════════════════════════════════════════════════
# Round-trip: GenBank → parse → sequence recovery
# ════════════════════════════════════════════════════════════════════
class TestGenBankRoundTrip:
    def test_roundtrip_recovers_sequence(self, simple_l0):
        gb = to_genbank(simple_l0)
        recovered = parse_genbank_sequence(gb)
        assert recovered == simple_l0.sequence.upper()

    def test_roundtrip_hbb(self, hbb_l0):
        gb = to_genbank(hbb_l0)
        recovered = parse_genbank_sequence(gb)
        assert recovered == HBB_DNA

    def test_roundtrip_eukaryotic(self, euk_l0):
        gb = to_genbank(euk_l0)
        recovered = parse_genbank_sequence(gb)
        assert recovered == euk_l0.sequence.upper()

    def test_roundtrip_preserves_length(self, hbb_l0):
        gb = to_genbank(hbb_l0)
        recovered = parse_genbank_sequence(gb)
        assert len(recovered) == len(HBB_DNA)

    def test_parse_features_returns_list(self, simple_l0):
        gb = to_genbank(simple_l0)
        features = parse_genbank_features(gb)
        assert isinstance(features, list)
        assert len(features) >= 2  # at least gene + CDS

    def test_reimported_sequence_can_rebuild_l0(self, hbb_l0):
        """Reconstruct an IR_L0 from a parsed GenBank record."""
        gb = to_genbank(hbb_l0)
        recovered_seq = parse_genbank_sequence(gb)
        rebuilt_l0 = IR_L0_GenomicDNA(
            sequence=recovered_seq, regions=[],
            organism=hbb_l0.organism, gene_name=hbb_l0.gene_name,
        )
        assert rebuilt_l0.sequence == hbb_l0.sequence.upper()


# ════════════════════════════════════════════════════════════════════
# FASTA — header format
# ════════════════════════════════════════════════════════════════════
class TestFastaHeader:
    def test_dna_header_format(self, simple_l0):
        fa = to_fasta_dna(simple_l0)
        header = fa.splitlines()[0]
        assert header == ">simple|e_coli|L0|len=9"

    def test_protein_header_format(self):
        ir = IR_L3_Polypeptide(sequence="MAK*", organism="e_coli", gene_name="t")
        fa = to_fasta_protein(ir)
        assert fa.splitlines()[0] == ">t|e_coli|L3|len=4"

    def test_rna_header_format(self):
        ir = IR_L2_MatureMRNA(
            sequence="AUGGCUUAA", five_utr="", cds="AUGGCUUAA",
            three_utr="", organism="e_coli", gene_name="t",
        )
        fa = to_fasta_rna(ir)
        assert fa.splitlines()[0] == ">t|e_coli|L2|len=9"

    def test_header_starts_with_greater_than(self, simple_l0):
        fa = to_fasta_dna(simple_l0)
        assert fa.startswith(">")

    def test_header_contains_gene_name(self, simple_l0):
        fa = to_fasta_dna(simple_l0)
        assert "simple" in fa.splitlines()[0]

    def test_header_contains_organism(self, simple_l0):
        fa = to_fasta_dna(simple_l0)
        assert "e_coli" in fa.splitlines()[0]

    def test_header_contains_level(self, simple_l0):
        fa = to_fasta_dna(simple_l0)
        assert "|L0|" in fa.splitlines()[0]


# ════════════════════════════════════════════════════════════════════
# FASTA — body / wrapping
# ════════════════════════════════════════════════════════════════════
class TestFastaBody:
    def test_dna_body_uppercase(self):
        ir = IR_L0_GenomicDNA(
            sequence="atggcttaa", regions=[],
            organism="e_coli", gene_name="lc",
        )
        fa = to_fasta_dna(ir)
        body = fa.splitlines()[1]
        assert body == "ATGGCTTAA"

    def test_protein_body_uppercase(self):
        ir = IR_L3_Polypeptide(sequence="mak*", organism="e_coli", gene_name="lc")
        fa = to_fasta_protein(ir)
        body = fa.splitlines()[1]
        assert body == "MAK*"

    def test_fasta_wraps_at_60_chars(self):
        # 120-char protein (no stop) → exactly 2 lines of 60 each.
        long_protein = "M" * 120
        ir = IR_L3_Polypeptide(
            sequence=long_protein, organism="e_coli", gene_name="long",
        )
        fa = to_fasta_protein(ir)
        body_lines = fa.splitlines()[1:]
        assert all(len(line) <= 60 for line in body_lines)
        # 120 chars → 2 lines of 60 each.
        assert len(body_lines) == 2
        assert len(body_lines[0]) == 60
        assert len(body_lines[1]) == 60

    def test_fasta_wraps_overflow_into_third_line(self):
        # 121-char protein → 2 full lines of 60 + 1 line of 1.
        long_protein = "M" * 120 + "*"
        ir = IR_L3_Polypeptide(
            sequence=long_protein, organism="e_coli", gene_name="long2",
        )
        fa = to_fasta_protein(ir)
        body_lines = fa.splitlines()[1:]
        assert len(body_lines) == 3
        assert len(body_lines[0]) == 60
        assert len(body_lines[1]) == 60
        assert body_lines[2] == "*"

    def test_short_sequence_single_line(self):
        ir = IR_L3_Polypeptide(sequence="MAK*", organism="e_coli", gene_name="s")
        fa = to_fasta_protein(ir)
        # Header + 1 body line
        assert len(fa.splitlines()) == 2

    def test_fasta_ends_with_newline(self, simple_l0):
        fa = to_fasta_dna(simple_l0)
        assert fa.endswith("\n")


# ════════════════════════════════════════════════════════════════════
# FASTA — dispatcher
# ════════════════════════════════════════════════════════════════════
class TestFastaDispatcher:
    def test_dispatcher_l0(self, simple_l0):
        fa = to_fasta(simple_l0)
        assert "|L0|" in fa

    def test_dispatcher_l2(self):
        ir = IR_L2_MatureMRNA(
            sequence="AUGGCUUAA", five_utr="", cds="AUGGCUUAA",
            three_utr="", organism="e_coli", gene_name="t",
        )
        fa = to_fasta(ir)
        assert "|L2|" in fa

    def test_dispatcher_l3(self):
        ir = IR_L3_Polypeptide(sequence="MAK*", organism="e_coli", gene_name="t")
        fa = to_fasta(ir)
        assert "|L3|" in fa

    def test_dispatcher_rejects_l4(self):
        # L4 isn't supported by FASTA codegen.
        from biocompiler.ir import IR_L4_FoldedProtein
        ir = IR_L4_FoldedProtein(sequence="MAK*", organism="e_coli")
        with pytest.raises(TypeError):
            to_fasta(ir)


# ════════════════════════════════════════════════════════════════════
# HBB end-to-end demo: L0 → GenBank, L3 → protein FASTA
# ════════════════════════════════════════════════════════════════════
class TestHBBDemo:
    def test_hbb_l0_to_genbank(self, hbb_l0):
        gb = to_genbank(hbb_l0)
        # Structural checks.
        assert gb.startswith("LOCUS")
        assert gb.rstrip().endswith("//")
        assert "FEATURES" in gb
        assert "ORIGIN" in gb
        # Sequence matches.
        assert parse_genbank_sequence(gb) == HBB_DNA
        # HBB name appears in LOCUS.
        assert "HBB" in gb.splitlines()[0]

    def test_hbb_l3_protein_matches_canonical(self, hbb_l3):
        """HBB protein from the IR pipeline matches the canonical sequence."""
        assert hbb_l3.sequence == HBB_PROTEIN

    def test_hbb_l3_to_fasta_protein(self, hbb_l3):
        fa = to_fasta_protein(hbb_l3)
        header = fa.splitlines()[0]
        body = "".join(fa.splitlines()[1:])
        assert header == ">HBB|Homo_sapiens|L3|len=32"
        assert body == HBB_PROTEIN

    def test_hbb_pipeline_to_fasta_protein(self, hbb_l0):
        """Full pipeline: L0 → L3 → protein FASTA."""
        l3 = compile_gene(hbb_l0, IRLevel.L3)
        fa = to_fasta_protein(l3)
        body = "".join(fa.splitlines()[1:])
        assert body == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"

    def test_hbb_genbank_has_hbb_features(self, hbb_l0):
        gb = to_genbank(hbb_l0)
        features = parse_genbank_features(gb)
        keys = [k for k, _ in features]
        assert "gene" in keys
        assert "CDS" in keys
        # HBB is 96 bp, no introns → CDS spans 1..96.
        fdict = dict(features)
        assert fdict["CDS"] == "1..96"


# ════════════════════════════════════════════════════════════════════
# Provenance / metadata propagation in codegen
# ════════════════════════════════════════════════════════════════════
class TestCodegenProvenance:
    def test_genbank_carries_organism_in_SOURCE(self):
        ir = IR_L0_GenomicDNA(
            sequence="ATGGCTTAA", regions=[],
            organism="Mus_musculus", gene_name="m",
        )
        gb = to_genbank(ir)
        assert "Mus musculus" in gb  # underscores → spaces

    def test_fasta_protein_carries_organism(self):
        ir = IR_L3_Polypeptide(
            sequence="MAK*", organism="Escherichia_coli", gene_name="g",
        )
        fa = to_fasta_protein(ir)
        assert "Escherichia_coli" in fa  # FASTA preserves underscore form

    def test_genbank_without_gene_name_uses_default(self):
        ir = IR_L0_GenomicDNA(
            sequence="ATGGCTTAA", regions=[], organism="e_coli",
        )
        gb = to_genbank(ir)
        # Should not crash and should still be valid.
        assert gb.startswith("LOCUS")
        # The fallback gene name "gene" should appear in /gene= qualifier.
        assert '/gene="gene"' in gb

    def test_fasta_without_gene_name_uses_BC_prefix(self):
        ir = IR_L3_Polypeptide(sequence="MAK*", organism="e_coli")
        fa = to_fasta_protein(ir)
        header = fa.splitlines()[0]
        # Falls back to BC_<hex>.
        assert header.startswith(">BC_")


# ════════════════════════════════════════════════════════════════════
# SBOL3 (RDF/Turtle) — synthetic biology interchange standard
# ════════════════════════════════════════════════════════════════════
class TestSBOL3Codegen:
    """SBOL3 codegen emits valid Turtle RDF with Component, Sequence,
    and SequenceFeature objects.

    These tests verify structural validity (prefixes, top-level object
    types), the displayId, the Sequence elements round-trip, and that
    CDS regions appear as SequenceFeature sub-objects with Range
    locations carrying the correct 1-based inclusive coordinates.
    """

    # ── Structural validity ────────────────────────────────────────
    def test_has_prefix_declarations(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "@prefix sbol:" in doc
        assert "@prefix rdf:" in doc
        assert "@prefix rdfs:" in doc

    def test_has_sbol_Component(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "a sbol:Component" in doc

    def test_has_sbol_Sequence(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "a sbol:Sequence" in doc

    def test_has_biopax_DnaRegion_type(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "biopax:DnaRegion" in doc

    def test_has_dna_encoding(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "sbol:encoding" in doc
        assert "enc:DNA" in doc

    def test_has_elements_with_sequence(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert 'sbol:elements "ATGGCTTAA"' in doc

    def test_document_ends_with_newline(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert doc.endswith("\n")

    # ── Component displayId ────────────────────────────────────────
    def test_component_displayId_uses_gene_name(self, simple_l0):
        doc = to_sbol3(simple_l0)
        # simple_l0 has gene_name="simple"
        assert 'sbol:displayId "simple"' in doc

    def test_component_name_uses_gene_name(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert 'sbol:name "simple"' in doc

    def test_component_has_description(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert 'sbol:description "Gene design compiled by BioCompiler"' in doc

    def test_component_has_hasSequence(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "sbol:hasSequence" in doc

    def test_component_has_hasFeature(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "sbol:hasFeature" in doc

    def test_unique_uri_based_on_gene_name(self, simple_l0):
        doc = to_sbol3(simple_l0)
        # The component URI should embed the gene_name "simple".
        assert "http://biocompiler.org/gene/simple/1" in doc

    # ── Sequence elements ─────────────────────────────────────────
    def test_sequence_elements_match_ir(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert f'sbol:elements "{simple_l0.sequence.upper()}"' in doc

    def test_sequence_uppercases_lowercase_input(self):
        ir = IR_L0_GenomicDNA(
            sequence="atggcttaa", regions=[],
            organism="e_coli", gene_name="lc",
        )
        doc = to_sbol3(ir)
        assert 'sbol:elements "ATGGCTTAA"' in doc

    def test_sequence_displayId_suffixed(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert 'sbol:displayId "simple_sequence"' in doc

    # ── Features: SequenceFeature + Range locations ───────────────
    def test_prokaryotic_has_gene_feature(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "a sbol:SequenceFeature" in doc
        # gene role: SO:0000704
        assert "sbol:role so:0000704" in doc

    def test_prokaryotic_has_cds_feature_whole_sequence(self, simple_l0):
        doc = to_sbol3(simple_l0)
        # CDS role: SO:0000316
        assert "sbol:role so:0000316" in doc
        # Whole sequence 9 bp → Range start=1, end=9
        assert "sbol:start 1" in doc
        assert "sbol:end 9" in doc

    def test_features_have_Range_locations(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "a sbol:Range" in doc

    def test_features_have_hasLocation(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "sbol:hasLocation" in doc

    def test_feature_locations_reference_sequence(self, simple_l0):
        doc = to_sbol3(simple_l0)
        # Range objects point back to the Sequence via sbol:sequence.
        assert "sbol:sequence <" in doc

    def test_prokaryotic_emits_two_features(self, simple_l0):
        """Prokaryotic gene: 1 gene feature + 1 CDS feature = 2 features."""
        doc = to_sbol3(simple_l0)
        # Count "a sbol:SequenceFeature" occurrences.
        n = doc.count("a sbol:SequenceFeature")
        assert n == 2, f"expected 2 SequenceFeatures, got {n}"

    # ── Eukaryotic: multiple regions produce multiple features ────
    def test_eukaryotic_emits_feature_per_region_plus_gene(self, euk_l0):
        """euk_l0 has 5 regions (5'UTR, exon, intron, exon, 3'UTR).
        Expected features: 1 gene + 2 CDS (from exons) + 5'UTR + intron + 3'UTR = 6.
        """
        doc = to_sbol3(euk_l0)
        n = doc.count("a sbol:SequenceFeature")
        assert n == 6, f"expected 6 SequenceFeatures, got {n}"

    def test_eukaryotic_has_two_cds_features(self, euk_l0):
        doc = to_sbol3(euk_l0)
        # Two CDS features → two Range blocks with role so:0000316
        n_cds = doc.count("sbol:role so:0000316")
        assert n_cds == 2, f"expected 2 CDS features, got {n_cds}"

    def test_eukaryotic_cds_coordinates_correct(self, euk_l0):
        """Exon1 [3,9) → Range start=4, end=9; Exon2 [15,21) → start=16, end=21."""
        doc = to_sbol3(euk_l0)
        assert "sbol:start 4" in doc
        assert "sbol:end 9" in doc
        assert "sbol:start 16" in doc
        assert "sbol:end 21" in doc

    def test_eukaryotic_has_utr_features(self, euk_l0):
        doc = to_sbol3(euk_l0)
        # 5'UTR role: SO:0000204, 3'UTR role: SO:0000205
        assert "sbol:role so:0000204" in doc
        assert "sbol:role so:0000205" in doc

    def test_eukaryotic_has_intron_feature(self, euk_l0):
        doc = to_sbol3(euk_l0)
        # Intron role: SO:0000188
        assert "sbol:role so:0000188" in doc

    def test_eukaryotic_has_six_Range_locations(self, euk_l0):
        """6 SequenceFeatures → 6 Range locations."""
        doc = to_sbol3(euk_l0)
        n = doc.count("a sbol:Range")
        assert n == 6, f"expected 6 Range locations, got {n}"

    def test_eukaryotic_gene_feature_spans_whole_sequence(self, euk_l0):
        doc = to_sbol3(euk_l0)
        # euk_l0 sequence is 24 bp; gene feature should span 1..24.
        assert "sbol:start 1" in doc
        assert "sbol:end 24" in doc

    # ── HBB end-to-end demo: L0 → SBOL3 ───────────────────────────
    def test_hbb_l0_to_sbol3_structure(self, hbb_l0):
        """HBB (96 bp, prokaryotic) → valid SBOL3 document."""
        doc = to_sbol3(hbb_l0)
        # Structural checks.
        assert "@prefix sbol:" in doc
        assert "a sbol:Component" in doc
        assert "a sbol:Sequence" in doc
        assert "a sbol:SequenceFeature" in doc
        assert "a sbol:Range" in doc
        # displayId and sequence.
        assert 'sbol:displayId "HBB"' in doc
        assert 'sbol:elements "' + HBB_DNA + '"' in doc

    def test_hbb_sbol3_uri_uses_gene_name(self, hbb_l0):
        doc = to_sbol3(hbb_l0)
        assert "http://biocompiler.org/gene/HBB/1" in doc

    def test_hbb_sbol3_has_gene_and_cds_features(self, hbb_l0):
        doc = to_sbol3(hbb_l0)
        # 1 gene feature + 1 CDS feature (whole sequence, no exons).
        n = doc.count("a sbol:SequenceFeature")
        assert n == 2
        # CDS spans the whole 96 bp sequence.
        assert "sbol:start 1" in doc
        assert "sbol:end 96" in doc

    def test_hbb_sbol3_sequence_length_correct(self, hbb_l0):
        doc = to_sbol3(hbb_l0)
        # HBB_DNA is 96 bp.
        assert len(HBB_DNA) == 96
        # The Range end of the gene feature = 96.
        assert "sbol:end 96" in doc

    # ── Edge cases & provenance ───────────────────────────────────
    def test_sbol3_falls_back_to_gene_when_no_name(self):
        ir = IR_L0_GenomicDNA(
            sequence="ATGGCTTAA", regions=[], organism="e_coli",
        )
        doc = to_sbol3(ir)
        assert 'sbol:displayId "gene"' in doc
        assert 'sbol:name "gene"' in doc

    def test_sbol3_sanitises_special_chars_in_name(self):
        ir = IR_L0_GenomicDNA(
            sequence="ATGGCTTAA", regions=[],
            organism="e_coli", gene_name="my gene!",
        )
        doc = to_sbol3(ir)
        # displayId must match [A-Za-z_][A-Za-z0-9_]* — spaces and !
        # get replaced with underscores.
        assert 'sbol:displayId "my_gene_"' in doc

    def test_sbol3_prefixes_for_rdf_validity(self, simple_l0):
        doc = to_sbol3(simple_l0)
        # Every @prefix line ends with a period.
        for line in doc.splitlines():
            if line.startswith("@prefix"):
                assert line.endswith(" ."), f"bad prefix line: {line!r}"

    def test_sbol3_no_duplicate_blank_lines_at_end(self, simple_l0):
        doc = to_sbol3(simple_l0)
        # Strip the trailing newline and verify no run of >1 blank line
        # at the very end of the document.
        assert not doc.endswith("\n\n\n"), "trailing blank lines"

    def test_sbol3_dna_region_type_uri(self, simple_l0):
        doc = to_sbol3(simple_l0)
        assert "biopax:DnaRegion" in doc

    def test_sbol3_constraints_present_for_each_feature(self, simple_l0):
        """Each SequenceFeature must have a Range location (hasLocation)."""
        doc = to_sbol3(simple_l0)
        n_features = doc.count("a sbol:SequenceFeature")
        n_locations = doc.count("sbol:hasLocation")
        assert n_locations == n_features, (
            f"{n_features} features but only {n_locations} hasLocation refs"
        )
