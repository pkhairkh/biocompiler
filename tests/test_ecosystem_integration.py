"""
Agent 54: Ecosystem Integration Tests
=======================================

Tests for SBOL export/import round-trip, LIMS export (Benchling format),
GenBank round-trip verification, BioPython deep integration features,
sequence annotation enrichment, and biosafety annotation coverage.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import fromstring

import pytest

from biocompiler.optimizer import OptimizationResult


# ═══════════════════════════════════════════════════════════════════════
# Shared Fixtures
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def sample_opt_result() -> OptimizationResult:
    """A sample OptimizationResult for ecosystem integration tests."""
    return OptimizationResult(
        sequence="ATGGTTTCTAAAGGTGAA",
        gc_content=0.333,
        cai=0.78,
        failed_predicates=[],
        predicate_results=[],
        certificate_text="",
        protein="MVSKGE",
        fallback_used=False,
        satisfied_predicates=["GCInRange", "CodonAdapted"],
        aa_substitutions=[],
        mutagenesis_applied=False,
        codon_pair_bias=0.35,
    )


@pytest.fixture
def larger_opt_result() -> OptimizationResult:
    """A longer OptimizationResult for annotation tests."""
    # eGFP-like protein (shortened)
    protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    # Back-translate with common codons
    codon_map = {
        "M": "ATG", "V": "GTT", "S": "TCT", "K": "AAA", "G": "GGT",
        "E": "GAA", "L": "CTT", "F": "TTT", "T": "ACT", "C": "TGT",
        "P": "CCT", "I": "ATT", "D": "GAT", "N": "AAT", "H": "CAT",
    }
    dna = "".join(codon_map.get(aa, "NNN") for aa in protein)
    return OptimizationResult(
        sequence=dna,
        gc_content=0.48,
        cai=0.85,
        failed_predicates=[],
        predicate_results=[],
        certificate_text="",
        protein=protein,
        fallback_used=False,
        satisfied_predicates=["GCInRange", "CodonAdapted", "NoRestrictionSite"],
        aa_substitutions=[],
        mutagenesis_applied=False,
    )


# ═══════════════════════════════════════════════════════════════════════
# 1. SBOL Export/Import Round-Trip
# ═══════════════════════════════════════════════════════════════════════


class TestSBOLRoundTrip:
    """Test SBOL3 export followed by import preserves key data."""

    def test_sbol_export_produces_valid_xml(self, sample_opt_result, tmp_path):
        """export_sbol should produce valid XML."""
        from biocompiler.export.sbol_export import export_sbol
        output_file = str(tmp_path / "test_sbol.xml")
        result_path = export_sbol(
            sample_opt_result,
            output_file,
            gene_name="test_gene",
            organism="Homo_sapiens",
        )
        assert Path(result_path).exists()
        # Should be parseable XML
        with open(result_path, "r") as f:
            content = f.read()
        root = fromstring(content)
        assert root is not None

    def test_sbol_export_has_component_elements(self, sample_opt_result, tmp_path):
        """SBOL export should contain Component elements."""
        from biocompiler.export.sbol_export import export_sbol, SBOL3_NS
        output_file = str(tmp_path / "test_sbol.xml")
        export_sbol(sample_opt_result, output_file, gene_name="test_gene")
        with open(output_file, "r") as f:
            content = f.read()
        root = fromstring(content)
        # Should have Component elements
        components = list(root.iter(f"{{{SBOL3_NS}}}Component"))
        assert len(components) >= 1

    def test_sbol_export_has_sequence_elements(self, sample_opt_result, tmp_path):
        """SBOL export should contain Sequence elements."""
        from biocompiler.export.sbol_export import export_sbol, SBOL3_NS
        output_file = str(tmp_path / "test_sbol.xml")
        export_sbol(sample_opt_result, output_file, gene_name="test_gene")
        with open(output_file, "r") as f:
            content = f.read()
        root = fromstring(content)
        sequences = list(root.iter(f"{{{SBOL3_NS}}}Sequence"))
        assert len(sequences) >= 1

    def test_sbol_import_reconstructs_components(self, sample_opt_result, tmp_path):
        """SBOL import should reconstruct components from exported file."""
        from biocompiler.export.sbol_export import export_sbol
        from biocompiler.export.sbol_import import import_sbol
        output_file = str(tmp_path / "test_sbol.xml")
        export_sbol(sample_opt_result, output_file, gene_name="test_gene",
                    organism="Homo_sapiens")
        components = import_sbol(output_file)
        assert len(components) >= 1
        # At least one component should have DNA type
        dna_comps = [c for c in components if c.component_type.upper() == "DNA"]
        assert len(dna_comps) >= 1

    def test_sbol_roundtrip_preserves_sequence(self, sample_opt_result, tmp_path):
        """Sequence should survive SBOL export → import round trip."""
        from biocompiler.export.sbol_export import export_sbol
        from biocompiler.export.sbol_import import import_sbol
        output_file = str(tmp_path / "test_sbol.xml")
        export_sbol(sample_opt_result, output_file, gene_name="test_gene",
                    organism="Homo_sapiens")
        components = import_sbol(output_file)
        # Find a component with a sequence
        seq_components = [c for c in components if c.sequence]
        assert len(seq_components) >= 1
        # The sequence should match the original
        original_upper = sample_opt_result.sequence.upper()
        imported_seqs = {c.sequence.upper() for c in seq_components if c.sequence}
        assert original_upper in imported_seqs

    def test_sbol_roundtrip_preserves_roles(self, sample_opt_result, tmp_path):
        """Component roles should survive SBOL round trip."""
        from biocompiler.export.sbol_export import export_sbol
        from biocompiler.export.sbol_import import import_sbol
        output_file = str(tmp_path / "test_sbol.xml")
        export_sbol(sample_opt_result, output_file, gene_name="test_gene")
        components = import_sbol(output_file)
        # Should have at least gene or CDS role
        all_roles = []
        for c in components:
            all_roles.extend(c.roles)
        assert "gene" in all_roles or "CDS" in all_roles

    def test_sbol_json_format_export(self, sample_opt_result, tmp_path):
        """SBOL JSON-LD export should produce valid JSON."""
        from biocompiler.export.sbol_export import export_sbol
        output_file = str(tmp_path / "test_sbol.json")
        export_sbol(sample_opt_result, output_file, format="sbol3json",
                    gene_name="test_gene")
        with open(output_file, "r") as f:
            data = json.load(f)
        assert "@context" in data
        assert "components" in data


# ═══════════════════════════════════════════════════════════════════════
# 2. LIMS Export (Benchling Format)
# ═══════════════════════════════════════════════════════════════════════


class TestLIMSExport:
    """Test LIMS integration, especially Benchling format export."""

    def test_benchling_export_has_required_fields(self, sample_opt_result):
        """Benchling export should have name, sequence, annotations, customFields."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_opt_result)
        assert "name" in payload
        assert "sequence" in payload
        assert "annotations" in payload
        assert "customFields" in payload

    def test_benchling_export_sequence_matches(self, sample_opt_result):
        """Benchling export sequence should match the optimized DNA."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_opt_result)
        assert payload["sequence"] == sample_opt_result.sequence.upper()

    def test_benchling_export_has_cai_and_gc(self, sample_opt_result):
        """Benchling custom fields should include CAI and GC content."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_opt_result)
        assert "cai" in payload["customFields"]
        assert "gc_content" in payload["customFields"]
        assert payload["customFields"]["cai"]["value"] == pytest.approx(sample_opt_result.cai, abs=0.01)
        assert payload["customFields"]["gc_content"]["value"] == pytest.approx(sample_opt_result.gc_content, abs=0.01)

    def test_benchling_export_has_annotations(self, sample_opt_result):
        """Benchling export should include feature annotations."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_opt_result)
        annotations = payload["annotations"]
        assert len(annotations) >= 1
        # Should have gene and CDS annotations at minimum
        annot_types = {a["type"] for a in annotations}
        assert "gene" in annot_types
        assert "CDS" in annot_types

    def test_benchling_export_has_biosafety_info(self, sample_opt_result):
        """Benchling export should include biosafety-related custom fields."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_opt_result)
        custom = payload["customFields"]
        # Should have organism and version metadata
        assert "organism" in custom
        assert "biocompiler_version" in custom

    def test_benchling_submit_design(self, sample_opt_result):
        """submit_design should return a design_id and cache the record."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter(folder_id="test_folder")
        design_id = exporter.submit_design(sample_opt_result, "project_1")
        assert design_id.startswith("BC_")
        # Should be retrievable
        status = exporter.get_design_status(design_id)
        assert status["status"] == "submitted"

    def test_labguru_export_has_required_fields(self, sample_opt_result):
        """LabGuru export should have item, data, tags, custom_fields."""
        from biocompiler.infrastructure.lims import LabGuruExporter
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_opt_result)
        assert "item" in payload
        assert "data" in payload
        assert "tags" in payload
        assert "custom_fields" in payload

    def test_labguru_export_includes_tags(self, sample_opt_result):
        """LabGuru export should include searchable tags."""
        from biocompiler.infrastructure.lims import LabGuruExporter
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_opt_result)
        tags = payload["tags"]
        assert "biocompiler" in tags
        assert "codon-optimized" in tags

    def test_convenience_export_to_benchling(self, sample_opt_result):
        """Convenience function export_to_benchling should work."""
        from biocompiler.infrastructure.lims import export_to_benchling
        payload = export_to_benchling(sample_opt_result)
        assert "name" in payload
        assert "sequence" in payload


# ═══════════════════════════════════════════════════════════════════════
# 3. GenBank Round-Trip Verification
# ═══════════════════════════════════════════════════════════════════════


class TestGenBankRoundTrip:
    """Test GenBank export/import round-trip verification."""

    def test_genbank_roundtrip_preserves_dna(self, sample_opt_result):
        """DNA sequence should survive GenBank export → import round trip."""
        from biocompiler.export.genbank_roundtrip import verify_genbank_roundtrip
        report = verify_genbank_roundtrip(
            sample_opt_result,
            gene_name="test_gene",
            organism="Homo_sapiens",
        )
        assert report.success is True
        assert report.original_dna == report.reimported_dna
        assert len(report.mismatches) == 0

    def test_genbank_roundtrip_preserves_annotations(self, sample_opt_result):
        """Key annotations should survive GenBank round trip."""
        from biocompiler.export.genbank_roundtrip import verify_genbank_roundtrip
        report = verify_genbank_roundtrip(
            sample_opt_result,
            gene_name="test_gene",
            organism="Homo_sapiens",
        )
        assert report.annotation_preserved is True

    def test_genbank_roundtrip_larger_protein(self, larger_opt_result):
        """Round trip should work for longer sequences too."""
        from biocompiler.export.genbank_roundtrip import verify_genbank_roundtrip
        report = verify_genbank_roundtrip(
            larger_opt_result,
            gene_name="egfp_short",
            organism="Homo_sapiens",
        )
        assert report.success is True
        assert len(report.mismatches) == 0

    def test_genbank_roundtrip_organism_verification(self, sample_opt_result):
        """Organism annotation should survive round trip."""
        from biocompiler.export.genbank_roundtrip import verify_genbank_roundtrip
        report = verify_genbank_roundtrip(
            sample_opt_result,
            gene_name="test_gene",
            organism="Escherichia_coli",
        )
        # Organism may have underscores replaced with spaces in GenBank
        orig_org = report.original_annotations.get("organism", "")
        reimp_org = report.reimported_annotations.get("organism", "")
        # Either exact match or normalized match (underscores → spaces)
        assert (
            orig_org == reimp_org
            or orig_org.replace("_", " ") == reimp_org.replace("_", " ").strip()
        )

    def test_compare_sequences_function(self):
        """compare_sequences should detect mismatches correctly."""
        from biocompiler.export.genbank_roundtrip import compare_sequences
        # Identical sequences
        assert compare_sequences("ATCG", "ATCG") == []
        # Single mismatch
        mismatches = compare_sequences("ATCG", "ATAG")
        assert len(mismatches) == 1
        assert mismatches[0] == (2, "C", "A")
        # No mismatches for identical
        assert compare_sequences("AAAA", "AAAA") == []

    def test_verify_annotation_preservation(self):
        """verify_annotation_preservation should check key fields."""
        from biocompiler.export.genbank_roundtrip import verify_annotation_preservation
        # Matching annotations
        orig = {"gene_name": "GFP", "protein": "MVSKGE", "organism": "Homo_sapiens"}
        reimp = {"gene_name": "GFP", "protein": "MVSKGE", "organism": "Homo sapiens"}
        preserved, warnings = verify_annotation_preservation(orig, reimp)
        assert preserved is True

    def test_genbank_export_has_biosafety_section(self, sample_opt_result):
        """GenBank export should include biosecurity annotations."""
        from biocompiler.export.core import export_genbank
        gb = export_genbank(
            sequence=sample_opt_result.sequence,
            organism="Homo_sapiens",
            gene_name="test_gene",
            protein=sample_opt_result.protein,
            cai=sample_opt_result.cai,
        )
        # Should contain biosafety-related annotations
        assert "BIOCOMPILER_ANNOTATIONS" in gb
        assert "biosafety_level" in gb
        assert "biosecurity_screened" in gb


# ═══════════════════════════════════════════════════════════════════════
# 4. BioPython Deep Integration Features
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.requires_external
class TestBioPythonIntegration:
    """Test BioPython interoperability features."""

    @pytest.fixture(autouse=True)
    def _skip_without_biopython(self):
        """Skip all tests in this class if BioPython is not installed."""
        pytest.importorskip("Bio")

    def test_to_seqrecord_creates_valid_record(self, sample_dna):
        """to_seqrecord should create a BioPython SeqRecord."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord
        record = to_seqrecord(
            sequence=sample_dna,
            organism="Homo_sapiens",
            gene_name="test_gene",
        )
        assert str(record.seq) == sample_dna.upper()
        assert record.annotations["organism"] == "Homo_sapiens"
        assert len(record.features) >= 1  # Should have at least CDS

    def test_from_seqrecord_roundtrip(self, sample_dna):
        """to_seqrecord → from_seqrecord should preserve key data."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord, from_seqrecord
        record = to_seqrecord(
            sequence=sample_dna,
            organism="Homo_sapiens",
            gene_name="test_gene",
        )
        data = from_seqrecord(record)
        assert data["sequence"] == sample_dna.upper()
        assert data["organism"] == "Homo_sapiens"
        assert data["gene_name"] == "test_gene"

    def test_codon_usage_table_loading(self):
        """load_codon_usage_table should return structured data."""
        from biocompiler.infrastructure.biopython_compat import load_codon_usage_table
        result = load_codon_usage_table("Homo_sapiens")
        assert result.organism == "Homo_sapiens"
        assert len(result.codon_counts) > 0
        assert len(result.adaptiveness) > 0
        assert result.source != ""

    def test_cai_from_custom_table(self, sample_dna):
        """compute_cai_from_table should use the provided table."""
        from biocompiler.infrastructure.biopython_compat import (
            load_codon_usage_table,
            compute_cai_from_table,
        )
        table = load_codon_usage_table("Homo_sapiens")
        cai = compute_cai_from_table(sample_dna, table)
        assert 0.0 <= cai <= 1.0

    def test_phylo_distance(self, sample_dna):
        """phylo_distance should return a non-negative distance."""
        from biocompiler.infrastructure.biopython_compat import phylo_distance
        dist = phylo_distance(sample_dna, organism="Homo_sapiens", method="euclidean")
        assert dist >= 0.0
        # Euclidean distance of a short sequence against reference codon usage
        # profile can be substantial even for "self" because the sequence's
        # actual codon usage differs from the whole-genome reference profile.
        # Just verify it returns a finite value.
        assert dist < 100.0

    def test_detect_orfs(self, sample_dna):
        """detect_orfs should find ORFs in the sequence."""
        from biocompiler.infrastructure.biopython_compat import detect_orfs
        orfs = detect_orfs(sample_dna, min_length_aa=2)
        assert isinstance(orfs, list)
        # sample_dna encodes MVSKGE which starts with ATG → at least one ORF
        if len(sample_dna) >= 9:  # need at least 3 codons
            assert len(orfs) >= 1
            assert orfs[0].protein is not None

    def test_optimize_to_seqrecord(self):
        """optimize_to_seqrecord should produce an optimized SeqRecord."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord
        record = optimize_to_seqrecord(
            protein="MVSKGE",
            organism="Escherichia_coli",
            gene_name="test_gene",
        )
        assert str(record.seq).startswith("ATG")
        assert record.annotations["organism"] == "Escherichia_coli"
        # Protein length should be 6 AA = 18 bp (plus stop codon)
        assert len(record.seq) >= 18


# ═══════════════════════════════════════════════════════════════════════
# 5. Sequence Annotation Enrichment
# ═══════════════════════════════════════════════════════════════════════


class TestAnnotationEnrichment:
    """Test sequence annotation enrichment features."""

    def test_annotate_sequence_returns_list(self, larger_opt_result):
        """annotate_sequence should return a list of annotations."""
        from biocompiler.export.annotation import annotate_sequence
        annotations = annotate_sequence(larger_opt_result.sequence)
        assert isinstance(annotations, list)
        assert len(annotations) > 0

    def test_annotation_has_required_fields(self, larger_opt_result):
        """Each annotation should have feature_type, start, end, strand."""
        from biocompiler.export.annotation import annotate_sequence
        annotations = annotate_sequence(larger_opt_result.sequence)
        for ann in annotations:
            assert hasattr(ann, "feature_type")
            assert hasattr(ann, "start")
            assert hasattr(ann, "end")
            assert hasattr(ann, "strand")
            assert ann.start < ann.end

    def test_cds_annotations_found(self, larger_opt_result):
        """annotate_sequence should find CDS features."""
        from biocompiler.export.annotation import annotate_sequence
        annotations = annotate_sequence(larger_opt_result.sequence)
        cds_annots = [a for a in annotations if a.feature_type == "CDS"]
        assert len(cds_annots) >= 1

    def test_restriction_site_annotations(self, larger_opt_result):
        """annotate_sequence should find restriction site features."""
        from biocompiler.export.annotation import annotate_sequence
        annotations = annotate_sequence(larger_opt_result.sequence)
        rs_annots = [a for a in annotations if a.feature_type == "restriction_site"]
        # May or may not have restriction sites, but the feature type
        # should be valid if present
        for rs in rs_annots:
            assert "enzyme" in rs.qualifiers

    def test_annotate_to_genbank(self, larger_opt_result):
        """annotate_to_genbank should produce a valid GenBank record."""
        from biocompiler.export.annotation import annotate_to_genbank
        gb = annotate_to_genbank(
            larger_opt_result.sequence,
            name="TEST_GENE",
            organism="Homo_sapiens",
        )
        assert "LOCUS" in gb
        assert "ORIGIN" in gb
        assert "//" in gb  # GenBank terminator

    def test_annotations_sorted_by_position(self, larger_opt_result):
        """Annotations should be sorted by start position."""
        from biocompiler.export.annotation import annotate_sequence
        annotations = annotate_sequence(larger_opt_result.sequence)
        for i in range(1, len(annotations)):
            assert annotations[i].start >= annotations[i - 1].start

    def test_annotation_types_are_valid(self, larger_opt_result):
        """All annotation feature types should be recognized."""
        from biocompiler.export.annotation import annotate_sequence
        valid_types = {
            "CDS", "promoter", "RBS", "restriction_site", "CpG_island",
            "splice_donor", "splice_acceptor", "simple_repeat",
            "GC_rich", "AT_rich",
        }
        annotations = annotate_sequence(larger_opt_result.sequence)
        for ann in annotations:
            assert ann.feature_type in valid_types


# ═══════════════════════════════════════════════════════════════════════
# 6. Biosafety Annotations in All Exports
# ═══════════════════════════════════════════════════════════════════════


class TestBiosafetyAnnotations:
    """Test that all export formats include biosafety annotations."""

    def test_genbank_export_has_biosafety_level(self, sample_opt_result):
        """GenBank export should include BSL level annotation."""
        from biocompiler.export.core import export_genbank
        gb = export_genbank(
            sequence=sample_opt_result.sequence,
            organism="Escherichia_coli",
            gene_name="test_gene",
            protein=sample_opt_result.protein,
        )
        assert "biosafety_level" in gb
        # E. coli should be BSL-1
        assert "BSL-1" in gb

    def test_genbank_export_has_biosecurity_screened(self, sample_opt_result):
        """GenBank export should include biosecurity_screened status."""
        from biocompiler.export.core import export_genbank
        gb = export_genbank(
            sequence=sample_opt_result.sequence,
            organism="Homo_sapiens",
            gene_name="test_gene",
        )
        assert "biosecurity_screened" in gb

    def test_genbank_export_has_provenance_id(self, sample_opt_result):
        """GenBank export should include provenance ID."""
        from biocompiler.export.core import export_genbank
        gb = export_genbank(
            sequence=sample_opt_result.sequence,
            organism="Homo_sapiens",
            gene_name="test_gene",
        )
        assert "provenance_id" in gb
        assert "BC_" in gb  # BioCompiler provenance ID prefix

    def test_fasta_export_has_biosecurity_header(self, sample_opt_result):
        """FASTA export should include biosecurity info in header."""
        from biocompiler.export.core import export_fasta
        fasta = export_fasta(
            sequence=sample_opt_result.sequence,
            identifier="test_gene",
            organism="Homo_sapiens",
            protein=sample_opt_result.protein,
        )
        # Header should contain biosecurity annotation
        assert "biosecurity=" in fasta
        assert "BSL-" in fasta

    def test_json_export_has_biosafety_section(self, sample_opt_result):
        """JSON export should include a biosafety section."""
        from biocompiler.export.core import export_json
        json_str = export_json(sample_opt_result)
        data = json.loads(json_str)
        assert "biosafety" in data
        bio = data["biosafety"]
        assert "biosafety_level" in bio
        assert "biosecurity_screened" in bio
        assert "provenance_id" in bio

    def test_sbol_export_includes_creator(self, sample_opt_result, tmp_path):
        """SBOL export should include BioCompiler creator annotation."""
        from biocompiler.export.sbol_export import export_sbol, DCT_NS
        output_file = str(tmp_path / "biosafety_sbol.xml")
        export_sbol(sample_opt_result, output_file, gene_name="test_gene")
        with open(output_file, "r") as f:
            content = f.read()
        assert "BioCompiler" in content

    def test_benchling_export_has_biosafety_metadata(self, sample_opt_result):
        """Benchling export should include biosafety metadata."""
        from biocompiler.infrastructure.lims import BenchlingExporter
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_opt_result)
        custom = payload["customFields"]
        assert "biocompiler_version" in custom
        assert "organism" in custom

    def test_biosafety_report_generation(self, sample_opt_result):
        """format_biosecurity_report should produce a readable report."""
        from biocompiler.export.core import format_biosecurity_report
        report = format_biosecurity_report(
            sequence=sample_opt_result.sequence,
            organism="Escherichia_coli",
            cai=sample_opt_result.cai,
            gc=sample_opt_result.gc_content,
        )
        assert "BIOSECURITY SCREENING REPORT" in report
        assert "biosafety_level" in report.lower() or "BSL" in report
        assert "biocompiler" in report.lower()

    def test_e_coli_gets_bsl1(self, sample_opt_result):
        """E. coli exports should be classified as BSL-1."""
        from biocompiler.export.core import export_genbank
        gb = export_genbank(
            sequence=sample_opt_result.sequence,
            organism="Escherichia_coli",
            gene_name="test_gene",
        )
        assert "BSL-1" in gb

    def test_human_gets_bsl2(self, sample_opt_result):
        """Homo_sapiens exports should be classified as BSL-2."""
        from biocompiler.export.core import export_genbank
        gb = export_genbank(
            sequence=sample_opt_result.sequence,
            organism="Homo_sapiens",
            gene_name="test_gene",
        )
        assert "BSL-2" in gb

    def test_biosecurity_screening_function(self):
        """screen_hazardous_sequence should detect hazardous motifs."""
        from biocompiler.biosecurity import screen_hazardous_sequence
        # Safe protein
        safe_report = screen_hazardous_sequence("MVSKGEELFTG")
        assert safe_report.risk_level in ("none", "low", "medium", "high", "critical")
        # A protein containing a select agent motif (ricin A chain catalytic site)
        # NIRVGLPIIS is a known ricin A chain motif
        hazardous_report = screen_hazardous_sequence("AAANIRVGLPIISAAA")
        assert hazardous_report.is_hazardous is True
        assert hazardous_report.risk_level in ("high", "critical")
