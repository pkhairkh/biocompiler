"""
Tests for biocompiler.sbol3_export — SBOL3 RDF/XML and JSON-LD generation.

Covers:
- generate_sbol3_xml: valid XML output, SBOL3 namespaces, Component, Sequence,
  Constraint objects
- generate_sbol3_json: valid JSON-LD output, @graph structure
- export_sbol3: primary API, format dispatch, OptimizationResult integration
- Input validation (empty sequence, empty protein name, empty organism)
- Metadata annotations (CAI, GC content, predicates, fallback, certificate)
- SBOL3 identity uniqueness
- Sequence content preservation
- Constraint element structure
- Namespace correctness (http://sbols.org/v3#)
"""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from typing import Any

import pytest

from biocompiler.sbol3_export import (
    generate_sbol3_xml,
    generate_sbol3_json,
    export_sbol3,
    _generate_identity,
    SBOL3_NS,
    RDF_NS,
    DCT_NS,
    BIOPAX_NS,
    BC_NS,
    DNAREGION_TYPE,
    IUPAC_DNA_ENCODING,
    SO_CDS,
    SO_GENE,
)


# ─── Test fixtures ────────────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"
SAMPLE_PROTEIN = "GFP"
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
# 1. XML output basic validity
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenerateSbol3XmlBasic:

    def test_returns_string(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert isinstance(result, str)

    def test_starts_with_xml_declaration(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert result.startswith("<?xml")

    def test_contains_rdf_root_element(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "rdf:RDF" in result

    def test_valid_xml_parseable(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        assert root is not None

    def test_output_is_not_empty(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert len(result) > 100


# ═══════════════════════════════════════════════════════════════════════════════
# 2. SBOL3 Component structure
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Component:

    def test_has_component_element(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        components = root.findall(".//sbol:Component", ns)
        assert len(components) >= 1

    def test_component_has_display_id(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        display_id = comp.find("sbol:displayId", ns)
        assert display_id is not None
        assert display_id.text == SAMPLE_PROTEIN

    def test_component_has_name(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        name = comp.find("sbol:name", ns)
        assert name is not None
        assert name.text == SAMPLE_PROTEIN

    def test_component_has_description(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        desc = comp.find("sbol:description", ns)
        assert desc is not None
        assert "BioCompiler" in desc.text

    def test_component_has_dna_region_type(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        types = comp.findall("sbol:type", ns)
        assert len(types) >= 1
        type_refs = [t.get(f"{{{RDF_NS}}}resource") for t in types]
        assert DNAREGION_TYPE in type_refs

    def test_component_has_roles(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        roles = comp.findall("sbol:role", ns)
        assert len(roles) >= 2
        role_refs = [r.get(f"{{{RDF_NS}}}resource") for r in roles]
        assert SO_GENE in role_refs
        assert SO_CDS in role_refs

    def test_component_has_sequence_reference(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        has_seq = comp.find("sbol:hasSequence", ns)
        assert has_seq is not None
        assert has_seq.get(f"{{{RDF_NS}}}resource") is not None

    def test_component_has_constraint_reference(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        assert comp is not None
        has_constraint = comp.find("sbol:hasConstraint", ns)
        assert has_constraint is not None
        assert has_constraint.get(f"{{{RDF_NS}}}resource") is not None


# ═══════════════════════════════════════════════════════════════════════════════
# 3. SBOL3 Sequence structure
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Sequence:

    def test_has_sequence_element(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        sequences = root.findall(".//sbol:Sequence", ns)
        assert len(sequences) >= 1

    def test_sequence_has_elements(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        assert seq is not None
        elements = seq.find("sbol:elements", ns)
        assert elements is not None
        assert elements.text == SAMPLE_SEQ.upper()

    def test_sequence_has_encoding(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        assert seq is not None
        encoding = seq.find("sbol:encoding", ns)
        assert encoding is not None
        resource = encoding.get(f"{{{RDF_NS}}}resource")
        assert resource == IUPAC_DNA_ENCODING

    def test_sequence_has_display_id(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        assert seq is not None
        display_id = seq.find("sbol:displayId", ns)
        assert display_id is not None
        assert SAMPLE_PROTEIN in display_id.text


# ═══════════════════════════════════════════════════════════════════════════════
# 4. SBOL3 Constraint structure
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Constraint:

    def test_has_constraint_element(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        constraints = root.findall(".//sbol:Constraint", ns)
        assert len(constraints) >= 1

    def test_constraint_has_subject(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        constraint = root.find(".//sbol:Constraint", ns)
        assert constraint is not None
        subject = constraint.find("sbol:subject", ns)
        assert subject is not None
        resource = subject.get(f"{{{RDF_NS}}}resource")
        assert resource is not None

    def test_constraint_has_object(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        constraint = root.find(".//sbol:Constraint", ns)
        assert constraint is not None
        obj = constraint.find("sbol:object", ns)
        assert obj is not None
        resource = obj.get(f"{{{RDF_NS}}}resource")
        assert resource is not None

    def test_constraint_has_restriction(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        constraint = root.find(".//sbol:Constraint", ns)
        assert constraint is not None
        restriction = constraint.find("sbol:restriction", ns)
        assert restriction is not None
        resource = restriction.get(f"{{{RDF_NS}}}resource")
        assert resource is not None
        assert "sbols.org/v3#" in resource

    def test_constraint_links_component_to_sequence(self):
        """The constraint subject should reference the component,
        and the object should reference the sequence."""
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}

        comp = root.find(".//sbol:Component", ns)
        constraint = root.find(".//sbol:Constraint", ns)
        seq = root.find(".//sbol:Sequence", ns)

        comp_id = comp.get(f"{{{RDF_NS}}}about")
        seq_id = seq.get(f"{{{RDF_NS}}}about")

        subject_ref = constraint.find("sbol:subject", ns).get(f"{{{RDF_NS}}}resource")
        object_ref = constraint.find("sbol:object", ns).get(f"{{{RDF_NS}}}resource")

        assert subject_ref == comp_id
        assert object_ref == seq_id


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Namespace validity
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Namespaces:

    def test_has_rdf_namespace(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert RDF_NS in result

    def test_has_sbol3_namespace(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert SBOL3_NS in result

    def test_has_dcterms_namespace(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert DCT_NS in result

    def test_has_biopax_namespace(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert BIOPAX_NS in result

    def test_has_biocompiler_namespace(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert BC_NS in result

    def test_sbol3_namespace_is_v3(self):
        """The SBOL namespace must be http://sbols.org/v3#."""
        assert SBOL3_NS == "http://sbols.org/v3#"


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Metadata annotations
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Metadata:

    def test_organism_in_output(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert SAMPLE_ORGANISM in result

    def test_protein_name_in_output(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert SAMPLE_PROTEIN in result

    def test_cai_in_output(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, cai=0.95
        )
        assert "bc:cai" in result
        assert "0.9500" in result

    def test_gc_content_in_output(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, gc_content=0.52
        )
        assert "bc:gcContent" in result
        assert "0.5200" in result

    def test_satisfied_predicates_in_output(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            satisfied_predicates=["GCInRange", "CodonAdapted"],
        )
        assert "bc:satisfiedPredicates" in result
        assert "GCInRange" in result

    def test_failed_predicates_in_output(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            failed_predicates=["HasStopCodon"],
        )
        assert "bc:failedPredicates" in result
        assert "HasStopCodon" in result

    def test_fallback_used_in_output(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, fallback_used=True
        )
        assert "bc:fallbackUsed" in result
        assert "true" in result

    def test_no_fallback_by_default(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "bc:fallbackUsed" not in result

    def test_biocompiler_version_in_output(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        from biocompiler import __version__
        assert __version__ in result

    def test_optimization_date_in_output(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "bc:optimizationDate" in result

    def test_sequence_length_in_output(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "bc:sequenceLength" in result
        assert str(len(SAMPLE_SEQ)) in result

    def test_creator_annotation(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "dcterms:creator" in result
        assert "BioCompiler" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Certificate annotation
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Certificate:

    def test_no_certificate_by_default(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert "bc:hasCertificate" not in result

    def test_certificate_annotation_present(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "bc:hasCertificate" in result
        assert "true" in result

    def test_certificate_design_id_in_output(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "bc:certificateDesignId" in result
        assert SAMPLE_CERTIFICATE["design_id"] in result

    def test_certificate_data_embedded(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=SAMPLE_CERTIFICATE,
        )
        assert "bc:certificateData" in result

    def test_certificate_xml_escaped(self):
        """Certificate data with special chars should be properly XML-escaped."""
        cert_with_special = dict(SAMPLE_CERTIFICATE)
        cert_with_special["design_id"] = "test<>&\"id"
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            certificate=cert_with_special,
        )
        # Should not raise XML parse error
        root = ET.fromstring(result)
        assert root is not None


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Input validation
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportSbol3Validation:

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError, match="Sequence must not be empty"):
            generate_sbol3_xml("", SAMPLE_PROTEIN, SAMPLE_ORGANISM)

    def test_empty_protein_name_raises(self):
        with pytest.raises(ValueError, match="Protein name must not be empty"):
            generate_sbol3_xml(SAMPLE_SEQ, "", SAMPLE_ORGANISM)

    def test_empty_organism_raises(self):
        with pytest.raises(ValueError, match="Organism must not be empty"):
            generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, "")

    def test_whitespace_sequence_normalized(self):
        """Spaces in sequences should be removed in the output."""
        result = generate_sbol3_xml("ATG GCA TGC", "MAM", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert " " not in elements.text
        assert elements.text == "ATGGCATGC"

    def test_lowercase_sequence_uppercased(self):
        result = generate_sbol3_xml("atggcatgc", "MAM", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert elements.text == "ATGGCATGC"


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Identity uniqueness
# ═══════════════════════════════════════════════════════════════════════════════

class TestSbol3Identity:

    def test_unique_identities_across_calls(self):
        """Each call should generate unique SBOL3 identities."""
        result1 = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        result2 = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root1 = ET.fromstring(result1)
        root2 = ET.fromstring(result2)
        ns = {"sbol": SBOL3_NS}
        comp1_id = root1.find(".//sbol:Component", ns).get(f"{{{RDF_NS}}}about")
        comp2_id = root2.find(".//sbol:Component", ns).get(f"{{{RDF_NS}}}about")
        assert comp1_id != comp2_id

    def test_generate_identity_format(self):
        identity = _generate_identity("https://example.com", "myComponent")
        assert identity.startswith("https://example.com/myComponent/")
        parts = identity.split("/")
        assert len(parts) >= 3


# ═══════════════════════════════════════════════════════════════════════════════
# 10. JSON-LD output
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenerateSbol3Json:

    def test_returns_valid_json(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        assert isinstance(data, dict)

    def test_has_context(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        assert "@context" in data
        assert "sbol" in data["@context"]
        assert data["@context"]["sbol"] == SBOL3_NS

    def test_has_graph(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        assert "@graph" in data
        assert len(data["@graph"]) >= 3  # Component, Sequence, Constraint

    def test_component_in_graph(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        components = [g for g in data["@graph"] if g.get("@type") == "sbol:Component"]
        assert len(components) == 1
        assert components[0]["sbol:displayId"] == SAMPLE_PROTEIN

    def test_sequence_in_graph(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        sequences = [g for g in data["@graph"] if g.get("@type") == "sbol:Sequence"]
        assert len(sequences) == 1
        assert sequences[0]["sbol:elements"] == SAMPLE_SEQ.upper()

    def test_constraint_in_graph(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        constraints = [g for g in data["@graph"] if g.get("@type") == "sbol:Constraint"]
        assert len(constraints) == 1
        assert "sbol:subject" in constraints[0]
        assert "sbol:object" in constraints[0]
        assert "sbol:restriction" in constraints[0]

    def test_cai_in_json(self):
        result = generate_sbol3_json(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, cai=0.95
        )
        data = json.loads(result)
        comp = [g for g in data["@graph"] if g.get("@type") == "sbol:Component"][0]
        assert comp["bc:cai"] == 0.95

    def test_gc_content_in_json(self):
        result = generate_sbol3_json(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, gc_content=0.52
        )
        data = json.loads(result)
        comp = [g for g in data["@graph"] if g.get("@type") == "sbol:Component"][0]
        assert comp["bc:gcContent"] == 0.52

    def test_predicates_in_json(self):
        result = generate_sbol3_json(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            satisfied_predicates=["GCInRange"],
            failed_predicates=["HasStopCodon"],
        )
        data = json.loads(result)
        comp = [g for g in data["@graph"] if g.get("@type") == "sbol:Component"][0]
        assert comp["bc:satisfiedPredicates"] == ["GCInRange"]
        assert comp["bc:failedPredicates"] == ["HasStopCodon"]

    def test_organism_in_json(self):
        result = generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        data = json.loads(result)
        comp = [g for g in data["@graph"] if g.get("@type") == "sbol:Component"][0]
        assert comp["bc:targetOrganism"] == SAMPLE_ORGANISM

    def test_json_validation_empty_sequence(self):
        with pytest.raises(ValueError, match="Sequence must not be empty"):
            generate_sbol3_json("", SAMPLE_PROTEIN, SAMPLE_ORGANISM)

    def test_json_validation_empty_protein_name(self):
        with pytest.raises(ValueError, match="Protein name must not be empty"):
            generate_sbol3_json(SAMPLE_SEQ, "", SAMPLE_ORGANISM)

    def test_json_validation_empty_organism(self):
        with pytest.raises(ValueError, match="Organism must not be empty"):
            generate_sbol3_json(SAMPLE_SEQ, SAMPLE_PROTEIN, "")


# ═══════════════════════════════════════════════════════════════════════════════
# 11. export_sbol3 — primary API
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportSbol3:

    def test_default_format_is_xml(self):
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        assert result.startswith("<?xml")

    def test_explicit_xml_format(self):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, format="sbol3"
        )
        assert result.startswith("<?xml")

    def test_json_format(self):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, format="sbol3json"
        )
        data = json.loads(result)
        assert "@graph" in data

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Unsupported SBOL format"):
            export_sbol3(
                SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, format="sbol2"
            )

    def test_gc_content_auto_computed(self):
        """When gc_content is not provided, it should be auto-computed from sequence."""
        result = export_sbol3(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        # Should contain a gcContent annotation (auto-computed)
        assert "bc:gcContent" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 12. OptimizationResult integration
# ═══════════════════════════════════════════════════════════════════════════════

class TestOptimizationResultIntegration:

    @pytest.fixture
    def mock_optimization_result(self):
        """Create a mock OptimizationResult for testing."""
        from biocompiler.optimization import OptimizationResult
        return OptimizationResult(
            sequence="ATGAGCAAAGGTGAAGAACTGTTTACCGGT",
            gc_content=0.433,
            cai=0.876,
            protein="MSKGEELFTG",
            failed_predicates=[],
            satisfied_predicates=["GCInRange", "CodonAdapted", "NoStopCodons"],
            fallback_used=False,
        )

    def test_with_optimization_result(self, mock_optimization_result):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            optimization_result=mock_optimization_result,
        )
        # Should contain CAI from optimization result
        assert "bc:cai" in result
        assert "0.8760" in result

    def test_optimization_result_extracts_gc(self, mock_optimization_result):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            optimization_result=mock_optimization_result,
        )
        assert "bc:gcContent" in result

    def test_optimization_result_extracts_predicates(self, mock_optimization_result):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            optimization_result=mock_optimization_result,
        )
        assert "GCInRange" in result
        assert "CodonAdapted" in result

    def test_optimization_result_extracts_fallback(self, mock_optimization_result):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            optimization_result=mock_optimization_result,
        )
        # fallback_used is False, so should not appear
        assert "bc:fallbackUsed" not in result

    def test_explicit_params_override_optimization_result(self, mock_optimization_result):
        """Explicitly provided CAI should override the optimization result's CAI."""
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            optimization_result=mock_optimization_result,
            cai=0.999,  # Override
        )
        assert "0.9990" in result

    def test_invalid_optimization_result_type_raises(self):
        with pytest.raises(TypeError, match="Expected OptimizationResult"):
            export_sbol3(
                SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
                optimization_result="not_a_result",
            )

    def test_optimization_result_json_format(self, mock_optimization_result):
        result = export_sbol3(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            optimization_result=mock_optimization_result,
            format="sbol3json",
        )
        data = json.loads(result)
        comp = [g for g in data["@graph"] if g.get("@type") == "sbol:Component"][0]
        assert "bc:cai" in comp


# ═══════════════════════════════════════════════════════════════════════════════
# 13. Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_very_short_sequence(self):
        result = generate_sbol3_xml("ATG", "M", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert elements.text == "ATG"

    def test_organism_with_special_chars(self):
        """Organism names with underscores should be preserved."""
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, "Homo_sapiens")
        assert "Homo_sapiens" in result
        ET.fromstring(result)

    def test_protein_name_with_spaces(self):
        """Spaces in protein name should be replaced with underscores in displayId."""
        result = generate_sbol3_xml(SAMPLE_SEQ, "Green Fluorescent Protein", SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        display_id = comp.find("sbol:displayId", ns)
        assert " " not in display_id.text
        assert display_id.text == "Green_Fluorescent_Protein"

    def test_xml_escape_special_chars_in_protein_name(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, "GFP<test>&more", SAMPLE_ORGANISM)
        root = ET.fromstring(result)  # Should be valid XML
        assert root is not None

    def test_none_cai_not_in_output(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, cai=None)
        assert "bc:cai" not in result

    def test_none_gc_explicit(self):
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM, gc_content=None)
        assert "bc:gcContent" not in result

    def test_custom_base_uri(self):
        result = generate_sbol3_xml(
            SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM,
            base_uri="https://myorg.org/sbol3",
        )
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        comp = root.find(".//sbol:Component", ns)
        comp_id = comp.get(f"{{{RDF_NS}}}about")
        assert comp_id.startswith("https://myorg.org/sbol3/")

    def test_sequence_preserved_in_xml(self):
        """The DNA sequence in the SBOL3 document should match the input exactly."""
        result = generate_sbol3_xml(SAMPLE_SEQ, SAMPLE_PROTEIN, SAMPLE_ORGANISM)
        root = ET.fromstring(result)
        ns = {"sbol": SBOL3_NS}
        seq = root.find(".//sbol:Sequence", ns)
        elements = seq.find("sbol:elements", ns)
        assert elements.text == SAMPLE_SEQ.upper()
