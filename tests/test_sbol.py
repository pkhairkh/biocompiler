"""
Tests for BioCompiler SBOL3 export and import.

Covers:
- SBOLComponent dataclass creation and validation
- export_sbol: single OptimizationResult → SBOL3 XML
- export_sbol: JSON-LD format
- export_sbol_collection: multiple results → SBOL3 Collection
- import_sbol: parse SBOL3 XML → list[SBOLComponent]
- import_sbol: parse SBOL3 JSON → list[SBOLComponent]
- Round-trip: export → import preserves key data
- sbol_to_genespecs: SBOL components → GeneSpec objects
- Measure extraction (CAI, GC)
- Role mapping (SO URIs ↔ role names)
- Edge cases: empty results, invalid format, missing sequences
"""

import json
import os
import tempfile
from pathlib import Path
from xml.etree.ElementTree import parse as xml_parse

import pytest

from biocompiler.sbol_export import (
    SBOLComponent,
    export_sbol,
    export_sbol_collection,
    SBOL3_NS,
    RDF_NS,
    SO_NS,
    SO_CDS,
    SO_PROMOTER,
    SO_TERMINATOR,
)
from biocompiler.sbol_import import (
    import_sbol,
    sbol_to_genespecs,
)


# ─── Fixtures ──────────────────────────────────────────────────────


@pytest.fixture
def sample_optimization_result():
    """Create a sample OptimizationResult for testing."""
    from biocompiler.optimization import OptimizationResult
    # Use a valid short protein: M (Met) → ATG, K (Lys) → AAA
    # DNA: ATG AAA TGA (stop codon) → but we need CDS only (no stop in protein)
    # Let's use M+K → ATG AAA (6 bp)
    return OptimizationResult(
        sequence="ATGAAAAAA",  # MKK → but let's make it longer for CAI validation
        gc_content=0.10,
        cai=0.95,
        protein="MKK",  # 3 amino acids = 9 bp... let me fix
        failed_predicates=[],
        satisfied_predicates=["GCInRange", "CodonAdapted"],
    )


@pytest.fixture
def valid_optimization_result():
    """Create a valid OptimizationResult with matching protein/sequence lengths."""
    from biocompiler.optimization import OptimizationResult
    # Protein "MSKGEELFTG" → 10 amino acids = 30 bp
    # Simple codon assignment: M=ATG, S=AGC, K=AAA, G=GGT, E=GAA, E=GAA,
    #                          L=CTG, F=TTC, T=ACC, G=GGT
    seq = "ATGAGCAAAGGTGAAGAACTGTTTACCGGT"
    protein = "MSKGEELFTG"
    assert len(seq) == len(protein) * 3, f"Seq length {len(seq)} != {len(protein)*3}"
    return OptimizationResult(
        sequence=seq,
        gc_content=0.433,
        cai=0.876,
        protein=protein,
        failed_predicates=[],
        satisfied_predicates=["GCInRange", "CodonAdapted", "NoStopCodons"],
    )


@pytest.fixture
def second_optimization_result():
    """A second OptimizationResult for collection tests."""
    from biocompiler.optimization import OptimizationResult
    # Protein "MVSKGEELFTG" → 11 amino acids = 33 bp
    seq = "ATGGTTAGCAAAGGTGAAGAACTGTTTACCGGT"
    protein = "MVSKGEELFTG"
    assert len(seq) == len(protein) * 3
    return OptimizationResult(
        sequence=seq,
        gc_content=0.424,
        cai=0.891,
        protein=protein,
        failed_predicates=[],
        satisfied_predicates=["GCInRange", "CodonAdapted"],
    )


@pytest.fixture
def tmp_dir():
    """Provide a temporary directory for SBOL output files."""
    with tempfile.TemporaryDirectory() as d:
        yield d


# ─── SBOLComponent Tests ──────────────────────────────────────────


class TestSBOLComponent:
    """Tests for the SBOLComponent dataclass."""

    def test_create_basic_component(self):
        comp = SBOLComponent(
            identity="https://example.org/gfp",
            display_id="gfp",
            component_type="DNA",
            sequence="ATGGCT",
            roles=["CDS"],
        )
        assert comp.display_id == "gfp"
        assert comp.component_type == "DNA"
        assert comp.sequence == "ATGGCT"
        assert comp.roles == ["CDS"]

    def test_display_id_normalisation(self):
        comp = SBOLComponent(
            identity="https://example.org/gfp gene",
            display_id="gfp gene",
            component_type="DNA",
            sequence="ATG",
            roles=["gene"],
        )
        assert comp.display_id == "gfp_gene"  # spaces → underscores

    def test_empty_display_id_gets_generated(self):
        comp = SBOLComponent(
            identity="https://example.org/x",
            display_id="",
            component_type="DNA",
            sequence="ATG",
            roles=["CDS"],
        )
        assert comp.display_id.startswith("component_")

    def test_invalid_component_type_raises(self):
        with pytest.raises(ValueError, match="component_type must be"):
            SBOLComponent(
                identity="https://example.org/x",
                display_id="x",
                component_type="RNA",
                sequence="AUG",
                roles=["CDS"],
            )

    def test_protein_component(self):
        comp = SBOLComponent(
            identity="https://example.org/prot",
            display_id="protein1",
            component_type="Protein",
            sequence="MSKGEELFTG",
            roles=["CDS"],
        )
        assert comp.component_type == "PROTEIN"  # normalised to upper


# ─── Export Tests ──────────────────────────────────────────────────


class TestExportSbol:
    """Tests for export_sbol function."""

    def test_export_xml_basic(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "test_sbol3.xml")
        result_path = export_sbol(
            valid_optimization_result,
            output,
            gene_name="gfp",
            organism="Escherichia_coli",
        )
        assert os.path.exists(result_path)

        # Parse and verify
        tree = xml_parse(result_path)
        root = tree.getroot()

        # Should have at least one Component
        components = list(root.iter(f"{{{SBOL3_NS}}}Component"))
        assert len(components) >= 1

        # Should have Sequence elements
        sequences = list(root.iter(f"{{{SBOL3_NS}}}Sequence"))
        assert len(sequences) >= 1

    def test_export_xml_contains_measures(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "test_measures.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        tree = xml_parse(output)
        root = tree.getroot()

        # Find Measure elements
        measures = list(root.iter(f"{{{SBOL3_NS}}}Measure"))
        assert len(measures) >= 2  # CAI and GC

        # Check that measure values exist
        measure_values = []
        for m in measures:
            val_elem = m.find(f"{{{SBOL3_NS}}}value")
            if val_elem is not None:
                measure_values.append(float(val_elem.text))
        assert len(measure_values) >= 2

    def test_export_xml_contains_roles(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "test_roles.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        tree = xml_parse(output)
        root = tree.getroot()

        # Find role elements
        roles = []
        for role_elem in root.iter(f"{{{SBOL3_NS}}}role"):
            role_uri = role_elem.get(f"{{{RDF_NS}}}resource", "")
            roles.append(role_uri)

        # Should have gene and CDS roles
        assert any("gene" in r.lower() or "0704" in r for r in roles)
        assert any("CDS" in r or "0316" in r for r in roles)

    def test_export_xml_has_provenance(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "test_prov.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        tree = xml_parse(output)
        root = tree.getroot()

        # Should have Activity element (provenance)
        activities = list(root.iter(f"{{{SBOL3_NS}}}Activity"))
        # Note: Activity uses PROV namespace, check that too
        prov_activities = list(root.iter(
            "{http://www.w3.org/ns/prov#}Activity"
        ))
        assert len(activities) + len(prov_activities) >= 1

        # Should have creator annotation
        creators = []
        for comp_elem in root.iter(f"{{{SBOL3_NS}}}Component"):
            for child in comp_elem:
                if "creator" in child.tag:
                    if child.text:
                        creators.append(child.text)
        assert any("BioCompiler" in c for c in creators)

    def test_export_json_format(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "test_sbol3.json")
        result_path = export_sbol(
            valid_optimization_result,
            output,
            format="sbol3json",
            gene_name="gfp",
        )
        assert os.path.exists(result_path)

        with open(result_path) as f:
            data = json.load(f)

        assert "components" in data
        assert len(data["components"]) >= 1
        assert "sequences" in data

    def test_export_invalid_format_raises(self, valid_optimization_result, tmp_dir):
        with pytest.raises(ValueError, match="Unsupported SBOL format"):
            export_sbol(
                valid_optimization_result,
                os.path.join(tmp_dir, "bad.xml"),
                format="sbol2",
            )

    def test_export_wrong_type_raises(self, tmp_dir):
        with pytest.raises(TypeError, match="Expected OptimizationResult"):
            export_sbol(
                "not_a_result",
                os.path.join(tmp_dir, "bad.xml"),
            )

    def test_export_creates_parent_dirs(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "sub", "dir", "test.xml")
        result_path = export_sbol(valid_optimization_result, output, gene_name="gfp")
        assert os.path.exists(result_path)

    def test_export_returns_absolute_path(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "test.xml")
        result_path = export_sbol(valid_optimization_result, output, gene_name="gfp")
        assert os.path.isabs(result_path)


class TestExportSbolCollection:
    """Tests for export_sbol_collection function."""

    def test_export_collection_xml(
        self, valid_optimization_result, second_optimization_result, tmp_dir
    ):
        output = os.path.join(tmp_dir, "collection.xml")
        result_path = export_sbol_collection(
            [valid_optimization_result, second_optimization_result],
            output,
            collection_name="test_collection",
            organism="Escherichia_coli",
        )
        assert os.path.exists(result_path)

        tree = xml_parse(result_path)
        root = tree.getroot()

        # Should have Collection element
        collections = list(root.iter(f"{{{SBOL3_NS}}}Collection"))
        assert len(collections) == 1

        # Collection should have members
        members = []
        for coll in collections:
            for child in coll:
                if "member" in child.tag:
                    members.append(child.get(f"{{{RDF_NS}}}resource", ""))
        assert len(members) >= 2  # Two gene-level components

    def test_export_collection_empty_results(self, tmp_dir):
        output = os.path.join(tmp_dir, "empty_collection.xml")
        result_path = export_sbol_collection([], output, collection_name="empty")
        assert os.path.exists(result_path)

    def test_export_collection_single_result(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "single_collection.xml")
        result_path = export_sbol_collection(
            [valid_optimization_result],
            output,
            collection_name="single",
        )
        assert os.path.exists(result_path)


# ─── Import Tests ──────────────────────────────────────────────────


class TestImportSbol:
    """Tests for import_sbol function."""

    def test_import_xml_round_trip(self, valid_optimization_result, tmp_dir):
        """Export then import should preserve key data."""
        output = os.path.join(tmp_dir, "roundtrip.xml")
        export_sbol(
            valid_optimization_result,
            output,
            gene_name="gfp",
            organism="Escherichia_coli",
        )

        components = import_sbol(output)
        assert len(components) >= 1

        # Find the gene component
        gene_comps = [c for c in components if "gene" in c.roles]
        assert len(gene_comps) >= 1

        gene = gene_comps[0]
        assert gene.display_id == "gfp"
        assert gene.component_type == "DNA"
        assert len(gene.sequence) > 0
        assert "gene" in gene.roles

    def test_import_json_round_trip(self, valid_optimization_result, tmp_dir):
        """Export then import JSON format should preserve key data."""
        output = os.path.join(tmp_dir, "roundtrip.json")
        export_sbol(
            valid_optimization_result,
            output,
            format="sbol3json",
            gene_name="gfp",
        )

        components = import_sbol(output)
        assert len(components) >= 1

    def test_import_preserves_sequence(self, valid_optimization_result, tmp_dir):
        """Sequence data should survive round-trip."""
        output = os.path.join(tmp_dir, "seq_test.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        components = import_sbol(output)
        # At least one component should have the original sequence
        seqs = [c.sequence for c in components if c.sequence]
        assert any(s.upper() == valid_optimization_result.sequence.upper() for s in seqs)

    def test_import_preserves_roles(self, valid_optimization_result, tmp_dir):
        """SO role annotations should survive round-trip."""
        output = os.path.join(tmp_dir, "roles_test.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        components = import_sbol(output)
        all_roles = []
        for c in components:
            all_roles.extend(c.roles)

        # Should have gene and CDS roles
        assert "gene" in all_roles
        assert "CDS" in all_roles

    def test_import_raw_xml_text(self, valid_optimization_result, tmp_dir):
        """import_sbol should accept raw XML text, not just file paths."""
        output = os.path.join(tmp_dir, "raw_text.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        with open(output) as f:
            xml_text = f.read()

        components = import_sbol(xml_text)
        assert len(components) >= 1

    def test_import_invalid_format_raises(self):
        """Non-XML, non-JSON text should raise FileFormatError."""
        from biocompiler.exceptions import FileFormatError
        with pytest.raises(FileFormatError):
            import_sbol("this is not valid SBOL at all")


# ─── sbol_to_genespecs Tests ──────────────────────────────────────


class TestSbolToGenespecs:
    """Tests for sbol_to_genespecs function."""

    def test_convert_cds_component(self):
        comp = SBOLComponent(
            identity="https://example.org/gfp",
            display_id="gfp",
            component_type="DNA",
            sequence="ATGAGCAAAGGTGAAGAACTGTTTACCGGT",  # MSKGEELFTG
            roles=["CDS"],
        )
        specs = sbol_to_genespecs([comp])
        assert len(specs) == 1
        assert specs[0].name == "gfp"
        assert specs[0].protein == "MSKGEELFTG"

    def test_convert_gene_component(self):
        comp = SBOLComponent(
            identity="https://example.org/insulin",
            display_id="insulin",
            component_type="DNA",
            sequence="ATGAGCAAAGGTGAAGAACTG",  # MSKGEEL = 7aa = 21bp
            roles=["gene"],
        )
        specs = sbol_to_genespecs([comp])
        assert len(specs) == 1
        assert specs[0].name == "insulin"

    def test_skip_protein_component(self):
        comp = SBOLComponent(
            identity="https://example.org/prot",
            display_id="protein1",
            component_type="Protein",
            sequence="MSKGEELFTG",
            roles=["CDS"],
        )
        specs = sbol_to_genespecs([comp])
        assert len(specs) == 0

    def test_skip_promoter_only(self):
        comp = SBOLComponent(
            identity="https://example.org/prom",
            display_id="promoter1",
            component_type="DNA",
            sequence="TTGACA",
            roles=["promoter"],
        )
        specs = sbol_to_genespecs([comp])
        assert len(specs) == 0  # promoter without CDS/gene role

    def test_skip_non_coding_length(self):
        comp = SBOLComponent(
            identity="https://example.org/short",
            display_id="short",
            component_type="DNA",
            sequence="ATGA",  # 4 bp, not divisible by 3
            roles=["CDS"],
        )
        specs = sbol_to_genespecs([comp])
        assert len(specs) == 0

    def test_skip_empty_sequence(self):
        comp = SBOLComponent(
            identity="https://example.org/empty",
            display_id="empty",
            component_type="DNA",
            sequence="",
            roles=["CDS"],
        )
        specs = sbol_to_genespecs([comp])
        assert len(specs) == 0

    def test_multiple_components(self):
        comps = [
            SBOLComponent(
                identity="https://example.org/g1",
                display_id="gene1",
                component_type="DNA",
                sequence="ATGAGCAAAGGTGAAGAACTG",  # 21bp = 7aa
                roles=["gene"],
            ),
            SBOLComponent(
                identity="https://example.org/g2",
                display_id="gene2",
                component_type="DNA",
                sequence="ATGGTTAGCAAAGGTGAAGAACTGTTTACC",  # 30bp = 10aa
                roles=["CDS"],
            ),
        ]
        specs = sbol_to_genespecs(comps)
        assert len(specs) == 2
        assert specs[0].name == "gene1"
        assert specs[1].name == "gene2"


# ─── Round-Trip Integration Tests ──────────────────────────────────


class TestRoundTrip:
    """End-to-end round-trip tests: export → import → genespecs."""

    def test_full_round_trip(self, valid_optimization_result, tmp_dir):
        """Export → Import → sbol_to_genespecs should produce valid GeneSpecs."""
        output = os.path.join(tmp_dir, "full_roundtrip.xml")
        export_sbol(
            valid_optimization_result,
            output,
            gene_name="gfp",
            organism="Escherichia_coli",
        )

        components = import_sbol(output)
        specs = sbol_to_genespecs(components)

        assert len(specs) >= 1
        # The first GeneSpec should have a protein that matches
        spec = specs[0]
        assert spec.name == "gfp"
        assert len(spec.protein) > 0

    def test_collection_round_trip(
        self, valid_optimization_result, second_optimization_result, tmp_dir
    ):
        """Collection export → import should find all gene components."""
        output = os.path.join(tmp_dir, "coll_roundtrip.xml")
        export_sbol_collection(
            [valid_optimization_result, second_optimization_result],
            output,
            collection_name="test_lib",
            organism="Escherichia_coli",
        )

        components = import_sbol(output)
        # Should find gene-level components for both results
        gene_comps = [c for c in components if "gene" in c.roles]
        assert len(gene_comps) >= 2

    def test_json_round_trip(self, valid_optimization_result, tmp_dir):
        """JSON format round-trip."""
        output = os.path.join(tmp_dir, "json_roundtrip.json")
        export_sbol(
            valid_optimization_result,
            output,
            format="sbol3json",
            gene_name="gfp",
        )

        components = import_sbol(output)
        assert len(components) >= 1

        # Convert to GeneSpecs
        specs = sbol_to_genespecs(components)
        assert len(specs) >= 1


# ─── Measure Extraction Tests ──────────────────────────────────────


class TestMeasureExtraction:
    """Tests for CAI and GC measure extraction from SBOL documents."""

    def test_measures_in_xml_export(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "measures.xml")
        export_sbol(valid_optimization_result, output, gene_name="gfp")

        tree = xml_parse(output)
        root = tree.getroot()

        # Extract measures
        from biocompiler.sbol_import import _extract_measures_xml
        measures = _extract_measures_xml(root)

        # At least one component should have measures
        assert len(measures) >= 1

        # Find the gene component measures
        for comp_id, comp_measures in measures.items():
            if "CAI" in comp_measures:
                assert 0 <= comp_measures["CAI"] <= 1.0
            if "GC_content" in comp_measures:
                assert 0 <= comp_measures["GC_content"] <= 1.0

    def test_measures_in_json_export(self, valid_optimization_result, tmp_dir):
        output = os.path.join(tmp_dir, "measures.json")
        export_sbol(
            valid_optimization_result,
            output,
            format="sbol3json",
            gene_name="gfp",
        )

        with open(output) as f:
            data = json.load(f)

        from biocompiler.sbol_import import _extract_measures_json
        measures = _extract_measures_json(data)

        assert len(measures) >= 1


# ─── Role Mapping Tests ───────────────────────────────────────────


class TestRoleMapping:
    """Tests for SO role URI ↔ role name mapping."""

    def test_cds_role_maps(self):
        from biocompiler.sbol_export import _resolve_role_uri
        uri = _resolve_role_uri("CDS")
        assert "0316" in uri  # SO:0000316

    def test_promoter_role_maps(self):
        from biocompiler.sbol_export import _resolve_role_uri
        uri = _resolve_role_uri("promoter")
        assert "0167" in uri  # SO:0000167

    def test_terminator_role_maps(self):
        from biocompiler.sbol_export import _resolve_role_uri
        uri = _resolve_role_uri("terminator")
        assert "0141" in uri  # SO:0000141

    def test_rbs_role_maps(self):
        from biocompiler.sbol_export import _resolve_role_uri
        uri = _resolve_role_uri("RBS")
        assert "0139" in uri  # SO:0000139

    def test_unknown_role_passes_through(self):
        from biocompiler.sbol_export import _resolve_role_uri
        uri = _resolve_role_uri("custom_role")
        assert "0000000" in uri  # Unknown → SO root

    def test_uri_passes_through(self):
        from biocompiler.sbol_export import _resolve_role_uri
        custom_uri = "http://example.org/custom_role"
        assert _resolve_role_uri(custom_uri) == custom_uri

    def test_reverse_mapping(self):
        from biocompiler.sbol_import import _so_uri_to_role
        assert _so_uri_to_role(SO_CDS) == "CDS"
        assert _so_uri_to_role(SO_PROMOTER) == "promoter"
        assert _so_uri_to_role(SO_TERMINATOR) == "terminator"
        assert _so_uri_to_role("http://unknown.org/") == "unknown"
