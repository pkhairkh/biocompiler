"""
Tests for the biocompiler.export_sbol module.

Covers:
- export_sbol3: basic output, valid XML structure, correct SBOL3 components
- RDF/XML structural validity (namespaces, Component, Sequence)
- Metadata annotations (organism, protein, CAI, GC)
- Certificate annotation embedding
- Input validation (empty sequence, empty protein, empty organism)
- SBOL3 identity uniqueness
- Sequence content preservation
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import pytest

from biocompiler.export_sbol import export_sbol3, _build_rdf_xml, _generate_sbol_identity


# ─── Test fixtures ────────────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"  # 45 bp
SAMPLE_PROTEIN = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSG"
SAMPLE_ORGANISM = "Escherichia_coli"

SAMPLE_CERTIFICATE = {
    "version": "1.0",
    "design_id": "abc123def456",
    "sequence": SAMPLE_SEQ,
    "types": [
        {"predicate": "gc_content", "verdict": "PASS"},
        {"predicate": "no_stop_codons", "verdict": "PASS"},
    ],
    "provenance": {
        "tool": "BioCompiler",
        "version": "10.0.0",
        "timestamp": "2025-01-15T12:00:00Z",
        "overall_status": "FULL_PASS",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Basic output validity
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportSbol3Basic:

    def test_returns_string(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert isinstance(result, str)

    def test_starts_with_xml_declaration(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert result.startswith("<?xml")

    def test_contains_rdf_root_element(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "rdf:RDF" in result

    def test_valid_xml(self):
        """The output should be parseable as XML."""
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        # Should not raise
        ET.fromstring(result)
        root = ET.fromstring(result)
        assert root is not None

    def test_output_is_not_empty(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert len(result) > 100


# ═══════════════════════════════════════════════════════════════════════════════
# 2. SBOL3 Component structure
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Component:

    def test_has_component_element(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        # Find Component elements in the SBOL3 namespace
        ns = {"sbol": "http://sbols.org/v3#"}
        components = root.findall(".//sbol:Component", ns)
        assert len(components) >= 1

    def test_component_has_display_id(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        display_id = comp.find("sbol:displayId", ns)
        assert display_id is not None
        assert display_id.text is not None

    def test_component_has_name(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        name = comp.find("sbol:name", ns)
        assert name is not None

    def test_component_has_description(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        desc = comp.find("sbol:description", ns)
        assert desc is not None
        assert "BioCompiler" in desc.text

    def test_component_has_sequence_reference(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        has_seq = comp.find("sbol:hasSequence", ns)
        assert has_seq is not None

    def test_component_has_type(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        types = comp.findall("sbol:type", ns)
        assert len(types) >= 1
        # Should reference DnaRegion
        type_refs = [t.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource") for t in types]
        assert any("DnaRegion" in str(ref) for ref in type_refs)

    def test_component_has_roles(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        roles = comp.findall("sbol:role", ns)
        assert len(roles) >= 1


# ═══════════════════════════════════════════════════════════════════════════════
# 3. SBOL3 Sequence structure
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Sequence:

    def test_has_sequence_element(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        sequences = root.findall(".//sbol:Sequence", ns)
        assert len(sequences) >= 1

    def test_sequence_has_elements(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        assert seq is not None
        elements = seq.find("sbol:elements", ns)
        assert elements is not None
        assert elements.text is not None
        assert len(elements.text) > 0

    def test_sequence_content_matches_input(self):
        """The DNA sequence in the SBOL3 document should match the input."""
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert elements.text == SAMPLE_SEQ.upper()

    def test_sequence_has_encoding(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        assert seq is not None
        encoding = seq.find("sbol:encoding", ns)
        assert encoding is not None
        resource = encoding.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource")
        assert "iupacNucleicAcid" in resource

    def test_sequence_has_display_id(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        assert seq is not None
        display_id = seq.find("sbol:displayId", ns)
        assert display_id is not None


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Namespace validity
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Namespaces:

    def test_has_rdf_namespace(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "http://www.w3.org/1999/02/22-rdf-syntax-ns#" in result

    def test_has_sbol3_namespace(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "http://sbols.org/v3#" in result

    def test_has_rdfs_namespace(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "http://www.w3.org/2000/01/rdf-schema#" in result

    def test_has_biocompiler_namespace(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "https://biocompiler.dev/vocab/" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Metadata annotations
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Metadata:

    def test_organism_in_output(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert SAMPLE_ORGANISM in result

    def test_protein_in_output(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert SAMPLE_PROTEIN in result

    def test_gc_content_in_output(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "gcContent" in result

    def test_biocompiler_version_in_output(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        from biocompiler import __version__
        assert __version__ in result

    def test_optimization_date_in_output(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "optimizationDate" in result

    def test_custom_metadata_in_output(self):
        metadata = {"geneName": "GFP", "project": "test_project"}
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            metadata=metadata,
        )
        assert "geneName" in result
        assert "GFP" in result
        assert "project" in result
        assert "test_project" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Certificate annotation
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Certificate:

    def test_no_certificate_by_default(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "hasCertificate" not in result

    def test_certificate_annotation_present(self):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "hasCertificate" in result
        assert "true" in result

    def test_certificate_design_id_in_output(self):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "certificateDesignId" in result
        assert SAMPLE_CERTIFICATE["design_id"] in result

    def test_certificate_status_in_output(self):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "certificateStatus" in result
        assert "FULL_PASS" in result

    def test_certificate_data_embedded(self):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "certificateData" in result

    def test_certificate_xml_escaped(self):
        """Certificate data should be properly XML-escaped."""
        cert_with_special = dict(SAMPLE_CERTIFICATE)
        cert_with_special["design_id"] = "test<>&\"id"
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=cert_with_special,
        )
        # Should not raise XML parse error
        ET.fromstring(result)


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Input validation
# ═══════════════════════════════════════════════════════════════════════════════
        root = ET.fromstring(result)
        assert root is not None

class TestExportSbol3Validation:

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError, match="Sequence must not be empty"):
            export_sbol3("", SAMPLE_PROTEIN, SAMPLE_ORGANISM)

    def test_empty_protein_raises(self):
        with pytest.raises(ValueError, match="Protein must not be empty"):
            export_sbol3(SAMPLE_SEQ, "", SAMPLE_ORGANISM)

    def test_empty_organism_raises(self):
        with pytest.raises(ValueError, match="Organism must not be empty"):
            export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, "")

    def test_whitespace_sequence_normalized(self):
        """Spaces in sequences should be removed."""
        result = export_sbol3("ATG GCA TGC", "MAM", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert " " not in elements.text
        assert elements.text == "ATGGCATGC"

    def test_lowercase_sequence_uppercased(self):
        """Lowercase sequence should be uppercased."""
        result = export_sbol3("atggcatgc", "MAM", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert elements.text == elements.text.upper()


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Identity uniqueness
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Identity:

    def test_unique_identities_across_calls(self):
        """Each call should generate unique SBOL3 identities."""
        result1 = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        result2 = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        # Extract identity URIs
        root1 = ET.fromstring(result1)
        root2 = ET.fromstring(result2)
        ns = {"sbol": "http://sbols.org/v3#"}
        comp1_id = root1.find(".//sbol:Component", ns).get(
            "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
        )
        comp2_id = root2.find(".//sbol:Component", ns).get(
            "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
        )
        assert comp1_id != comp2_id

    def test_generate_sbol_identity_format(self):
        identity = _generate_sbol_identity("https://example.com", "myComponent")
        assert identity.startswith("https://example.com/myComponent/")
        # Should have a UUID suffix
        parts = identity.split("/")
        assert len(parts) >= 3


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Manual RDF/XML builder
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuildRdfXml:

    def test_produces_valid_xml(self):
        result = _build_rdf_xml(
            component_identity="https://example.com/comp/1",
            component_display_id="comp1",
            component_name="Test Component",
            component_description="A test",
            sequence_identity="https://example.com/seq/1",
            sequence_display_id="seq1",
            sequence_elements="ATGCATGC",
            sequence_length=8,
            organism="Escherichia_coli",
            protein="MAM",
        )
        ET.fromstring(result)  # Should not raise
        root = ET.fromstring(result)
        assert root is not None

    def test_includes_all_metadata(self):
        result = _build_rdf_xml(
            component_identity="https://example.com/comp/1",
            component_display_id="comp1",
            component_name="Test",
            component_description="Test desc",
            sequence_identity="https://example.com/seq/1",
            sequence_display_id="seq1",
            sequence_elements="ATGCATGC",
            sequence_length=8,
            organism="Escherichia_coli",
            protein="MAM",
            metadata={"cai": "0.95", "gcContent": "0.5000"},
        )
        assert "cai" in result
        assert "0.95" in result
        assert "gcContent" in result

    def test_certificate_embedding(self):
        result = _build_rdf_xml(
            component_identity="https://example.com/comp/1",
            component_display_id="comp1",
            component_name="Test",
            component_description="Test desc",
            sequence_identity="https://example.com/seq/1",
            sequence_display_id="seq1",
            sequence_elements="ATGCATGC",
            sequence_length=8,
            organism="Escherichia_coli",
            protein="MAM",
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "hasCertificate" in result
        assert "certificateDesignId" in result
        assert SAMPLE_CERTIFICATE["design_id"] in result


# ═══════════════════════════════════════════════════════════════════════════════
# 10. Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportSbol3EdgeCases:

    def test_very_short_sequence(self):
        result = export_sbol3("ATG", "M", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": "http://sbols.org/v3#"}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert elements.text == "ATG"

    def test_organism_with_special_chars(self):
        """Organism names with underscores should be preserved."""
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, "Homo_sapiens")
        assert "Homo_sapiens" in result
        ET.fromstring(result)  # Should be valid XML

    def test_none_metadata(self):
        """None metadata should not cause errors."""
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, metadata=None)
        ET.fromstring(result)
        root = ET.fromstring(result)
        assert root is not None

    def test_none_certificate(self):
        """None certificate should not cause errors."""
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, certificate=None)
        ET.fromstring(result)
        assert "hasCertificate" not in result

    def test_empty_metadata_dict(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, metadata={})
        ET.fromstring(result)
        root = ET.fromstring(result)
        assert root is not None

    def test_metadata_with_none_value(self):
        """Metadata keys with None values should be skipped."""
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            metadata={"geneName": "GFP", "nullField": None},
        )
        assert "geneName" in result
        assert "GFP" in result
        ET.fromstring(result)
