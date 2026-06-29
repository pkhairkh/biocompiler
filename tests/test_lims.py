"""Tests for BioCompiler LIMS Integration Hooks.

Covers:
- LIMSIntegration base class interface
- BenchlingExporter: payload format, annotations, custom fields
- LabGuruExporter: payload format, tags, custom fields
- LIMSSubmissionRecord: dataclass behavior
- Convenience functions: export_to_benchling, export_to_labguru
- Cross-export delegation (BenchlingExporter.export_to_labguru and vice versa)
"""

import pytest
from biocompiler.optimizer import OptimizationResult
from biocompiler.infrastructure.lims import (
    LIMSIntegration,
    BenchlingExporter,
    LabGuruExporter,
    LIMSSubmissionRecord,
    export_to_benchling,
    export_to_labguru,
)


# ────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture
def sample_result() -> OptimizationResult:
    """A sample OptimizationResult for LIMS tests."""
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
    )


@pytest.fixture
def result_with_mutagenesis() -> OptimizationResult:
    """OptimizationResult with mutagenesis applied."""
    return OptimizationResult(
        sequence="ATGGTTTCTAAAGGTGAA",
        gc_content=0.333,
        cai=0.75,
        failed_predicates=["NoCrypticSplice"],
        predicate_results=[],
        certificate_text="",
        protein="MVSKGE",
        fallback_used=True,
        satisfied_predicates=["GCInRange"],
        aa_substitutions=[
            {"position": 3, "original": "K", "substitution": "R"},
        ],
        mutagenesis_applied=True,
        codon_pair_bias=0.42,
    )


# ────────────────────────────────────────────────────────────
# LIMSSubmissionRecord
# ────────────────────────────────────────────────────────────

class TestLIMSSubmissionRecord:
    """Tests for the LIMSSubmissionRecord dataclass."""

    def test_default_submitted_at(self):
        """submitted_at should be auto-populated with ISO timestamp."""
        record = LIMSSubmissionRecord(
            design_id="BC_test_12345678",
            project_id="proj_1",
        )
        assert record.submitted_at  # Not empty
        assert "T" in record.submitted_at  # ISO format

    def test_custom_submitted_at(self):
        """submitted_at can be explicitly set."""
        record = LIMSSubmissionRecord(
            design_id="BC_test_12345678",
            project_id="proj_1",
            submitted_at="2024-01-15T12:00:00Z",
        )
        assert record.submitted_at == "2024-01-15T12:00:00Z"

    def test_default_status(self):
        """Default status should be 'submitted'."""
        record = LIMSSubmissionRecord(
            design_id="BC_test_12345678",
            project_id="proj_1",
        )
        assert record.status == "submitted"


# ────────────────────────────────────────────────────────────
# LIMSIntegration base class
# ────────────────────────────────────────────────────────────

class TestLIMSIntegrationBase:
    """Tests for the LIMSIntegration abstract base class."""

    def test_cannot_instantiate_directly(self):
        """LIMSIntegration is abstract and cannot be instantiated."""
        with pytest.raises(TypeError):
            LIMSIntegration()  # type: ignore[abstract]

    def test_subclass_must_implement_methods(self):
        """Subclass without all abstract methods cannot be instantiated."""
        class IncompleteLIMS(LIMSIntegration):
            def submit_design(self, result, project_id):
                return "id"
            # Missing: get_design_status, export_to_benchling, export_to_labguru

        with pytest.raises(TypeError):
            IncompleteLIMS()  # type: ignore[abstract]

    def test_complete_subclass_can_instantiate(self):
        """A fully implemented subclass can be instantiated."""
        class CompleteLIMS(LIMSIntegration):
            def submit_design(self, result, project_id):
                return "design_123"
            def get_design_status(self, design_id):
                return {"status": "ok"}
            def export_to_benchling(self, result):
                return {}
            def export_to_labguru(self, result):
                return {}

        lims = CompleteLIMS()
        assert lims.submit_design(None, "proj") == "design_123"

    def test_generate_design_id(self):
        """_generate_design_id should produce unique IDs."""
        class CompleteLIMS(LIMSIntegration):
            def submit_design(self, result, project_id):
                return "id"
            def get_design_status(self, design_id):
                return {}
            def export_to_benchling(self, result):
                return {}
            def export_to_labguru(self, result):
                return {}

        lims = CompleteLIMS()
        id1 = lims._generate_design_id("proj1")
        id2 = lims._generate_design_id("proj1")
        assert id1.startswith("BC_proj1_")
        assert id2.startswith("BC_proj1_")
        assert id1 != id2  # UUID component makes them unique

    def test_base_url_trailing_slash_stripped(self):
        """Trailing slash is stripped from base_url."""
        class CompleteLIMS(LIMSIntegration):
            def submit_design(self, result, project_id):
                return "id"
            def get_design_status(self, design_id):
                return {}
            def export_to_benchling(self, result):
                return {}
            def export_to_labguru(self, result):
                return {}

        lims = CompleteLIMS(base_url="https://example.com/")
        assert lims.base_url == "https://example.com"


# ────────────────────────────────────────────────────────────
# BenchlingExporter
# ────────────────────────────────────────────────────────────

class TestBenchlingExporter:
    """Tests for the BenchlingExporter class."""

    def test_export_basic_fields(self, sample_result):
        """Export should produce name, sequence, annotations, customFields."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        assert "name" in payload
        assert "sequence" in payload
        assert "annotations" in payload
        assert "customFields" in payload
        assert payload["sequence"] == "ATGGTTTCTAAAGGTGAA"

    def test_export_name_format(self, sample_result):
        """Name should include organism and length."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        assert "18bp" in payload["name"]
        # Default organism is 'unknown' since sample_result has no organism_name
        assert "unknown" in payload["name"] or "BioCompiler" in payload["name"]

    def test_export_sequence_uppercase(self, sample_result):
        """Sequence in export should be uppercase."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)
        assert payload["sequence"] == payload["sequence"].upper()

    def test_export_annotations(self, sample_result):
        """Annotations should include gene, CDS, and regulatory features."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        annotations = payload["annotations"]
        assert len(annotations) >= 3  # gene, CDS, regulatory

        ann_types = [a["type"] for a in annotations]
        assert "gene" in ann_types
        assert "CDS" in ann_types
        assert "regulatory" in ann_types

    def test_export_annotation_bounds(self, sample_result):
        """Annotations should cover the full sequence."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        for ann in payload["annotations"]:
            assert ann["start"] == 0
            assert ann["end"] == len(sample_result.sequence)
            assert ann["strand"] in (1, -1)

    def test_export_custom_fields(self, sample_result):
        """Custom fields should include CAI, GC, organism, version."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        cf = payload["customFields"]
        assert "cai" in cf
        assert "gc_content" in cf
        assert "organism" in cf
        assert "biocompiler_version" in cf
        assert "sequence_length" in cf
        assert "optimization_date" in cf

        # Check values
        assert cf["cai"]["value"] == pytest.approx(0.78, abs=0.01)
        assert cf["gc_content"]["value"] == pytest.approx(0.333, abs=0.01)
        assert cf["sequence_length"]["value"] == 18

    def test_export_protein_length_custom_field(self, sample_result):
        """Protein length custom field should be present when protein is provided."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        cf = payload["customFields"]
        assert "protein_length" in cf
        assert cf["protein_length"]["value"] == 6

    def test_export_satisfied_predicates(self, sample_result):
        """Satisfied predicates should be in custom fields."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        cf = payload["customFields"]
        assert "satisfied_predicates" in cf
        assert "GCInRange" in cf["satisfied_predicates"]["value"]

    def test_export_no_protein_no_length(self):
        """Without protein, protein_length should be absent."""
        result = OptimizationResult(
            sequence="ATGGTTTCTAAAGGTGAA",
            gc_content=0.333,
            cai=0.78,
        )
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(result)

        cf = payload["customFields"]
        assert "protein_length" not in cf

    def test_export_mutagenesis_custom_fields(self, result_with_mutagenesis):
        """Mutagenesis metadata should appear in custom fields."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(result_with_mutagenesis)

        cf = payload["customFields"]
        assert "mutagenesis_applied" in cf
        assert cf["mutagenesis_applied"]["value"] is True
        assert "aa_substitutions" in cf

    def test_export_failed_predicates(self, result_with_mutagenesis):
        """Failed predicates should be in custom fields."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(result_with_mutagenesis)

        cf = payload["customFields"]
        assert "failed_predicates" in cf
        assert "NoCrypticSplice" in cf["failed_predicates"]["value"]

    def test_export_codon_pair_bias(self, result_with_mutagenesis):
        """Codon pair bias should be in custom fields when present."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(result_with_mutagenesis)

        cf = payload["customFields"]
        assert "codon_pair_bias" in cf
        assert cf["codon_pair_bias"]["value"] == pytest.approx(0.42, abs=0.01)

    def test_export_is_circular(self, sample_result):
        """Export should default to linear (isCircular=False)."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)
        assert payload["isCircular"] is False

    def test_folder_id(self, sample_result):
        """folder_id should be included in payload when set."""
        exporter = BenchlingExporter(folder_id="folder_abc123")
        payload = exporter.export_to_benchling(sample_result)
        assert payload["folderId"] == "folder_abc123"

    def test_submit_design(self, sample_result):
        """submit_design should return a design ID and cache the record."""
        exporter = BenchlingExporter()
        design_id = exporter.submit_design(sample_result, "proj_1")

        assert design_id.startswith("BC_proj_1_")
        # Check status
        status = exporter.get_design_status(design_id)
        assert status["status"] == "submitted"
        assert status["lims_system"] == "benchling"

    def test_get_design_status_not_found(self):
        """get_design_status should return 'not_found' for unknown IDs."""
        exporter = BenchlingExporter()
        status = exporter.get_design_status("nonexistent_id")
        assert status["status"] == "not_found"


# ────────────────────────────────────────────────────────────
# LabGuruExporter
# ────────────────────────────────────────────────────────────

class TestLabGuruExporter:
    """Tests for the LabGuruExporter class."""

    def test_export_basic_structure(self, sample_result):
        """Export should produce item, data, tags, custom_fields."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        assert "item" in payload
        assert "data" in payload
        assert "tags" in payload
        assert "custom_fields" in payload

    def test_export_item_fields(self, sample_result):
        """Item should have name, description, type."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        item = payload["item"]
        assert "name" in item
        assert "description" in item
        assert item["type"] == "DNA Sequence"
        assert "18bp" in item["name"]

    def test_export_data_is_sequence(self, sample_result):
        """data field should be the raw DNA sequence."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        assert payload["data"] == "ATGGTTTCTAAAGGTGAA"

    def test_export_tags(self, sample_result):
        """Tags should include biocompiler, organism, cai, gc."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        tags = payload["tags"]
        assert "biocompiler" in tags
        assert "codon-optimized" in tags
        assert any("cai:" in t for t in tags)
        assert any("gc:" in t for t in tags)

    def test_export_all_predicates_pass_tag(self, sample_result):
        """All-predicates-pass tag should be present when no predicates fail."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        tags = payload["tags"]
        assert "all-predicates-pass" in tags

    def test_export_mutagenesis_tag(self, result_with_mutagenesis):
        """Mutagenesis tag should be present when mutagenesis was applied."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(result_with_mutagenesis)

        tags = payload["tags"]
        assert "mutagenesis" in tags

    def test_export_fallback_tag(self, result_with_mutagenesis):
        """Fallback-used tag should be present when fallback was used."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(result_with_mutagenesis)

        tags = payload["tags"]
        assert "fallback-used" in tags

    def test_export_custom_fields(self, sample_result):
        """Custom fields should include CAI, GC, organism, version."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        cf = payload["custom_fields"]
        assert "cai" in cf
        assert "gc_content" in cf
        assert "organism" in cf
        assert "biocompiler_version" in cf
        assert "sequence_length" in cf
        assert cf["cai"] == pytest.approx(0.78, abs=0.01)
        assert cf["gc_content"] == pytest.approx(0.333, abs=0.01)

    def test_export_protein_custom_fields(self, sample_result):
        """Protein metadata should be in custom fields when available."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        cf = payload["custom_fields"]
        assert "protein_length" in cf
        assert cf["protein_length"] == 6
        assert "protein_sequence" in cf
        assert cf["protein_sequence"] == "MVSKGE"

    def test_export_codon_pair_bias(self, result_with_mutagenesis):
        """Codon pair bias should be in custom fields when present."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(result_with_mutagenesis)

        cf = payload["custom_fields"]
        assert "codon_pair_bias" in cf
        assert cf["codon_pair_bias"] == pytest.approx(0.42, abs=0.01)

    def test_export_description_content(self, sample_result):
        """Description should contain key optimization info."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        desc = payload["item"]["description"]
        assert "BioCompiler" in desc
        assert "CAI" in desc
        assert "GC" in desc

    def test_submit_design(self, sample_result):
        """submit_design should return a design ID and cache the record."""
        exporter = LabGuruExporter()
        design_id = exporter.submit_design(sample_result, "proj_1")

        assert design_id.startswith("BC_proj_1_")
        status = exporter.get_design_status(design_id)
        assert status["status"] == "submitted"
        assert status["lims_system"] == "labguru"

    def test_get_design_status_not_found(self):
        """get_design_status should return 'not_found' for unknown IDs."""
        exporter = LabGuruExporter()
        status = exporter.get_design_status("nonexistent_id")
        assert status["status"] == "not_found"

    def test_project_id_in_item(self, sample_result):
        """project_id should be set in item when configured."""
        exporter = LabGuruExporter(project_id=42)
        payload = exporter.export_to_labguru(sample_result)

        assert payload["item"]["project_id"] == 42

    def test_submit_with_project_id(self, sample_result):
        """submit_design should set project_id in payload."""
        exporter = LabGuruExporter(project_id=99)
        design_id = exporter.submit_design(sample_result, "proj_1")

        record = exporter._submission_cache[design_id]
        assert record.payload["item"]["project_id"] == 99


# ────────────────────────────────────────────────────────────
# Cross-export delegation
# ────────────────────────────────────────────────────────────

class TestCrossExport:
    """Test that BenchlingExporter can export to LabGuru and vice versa."""

    def test_benchling_exports_to_labguru(self, sample_result):
        """BenchlingExporter.export_to_labguru should produce LabGuru format."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_labguru(sample_result)

        assert "item" in payload
        assert "data" in payload
        assert "tags" in payload

    def test_labguru_exports_to_benchling(self, sample_result):
        """LabGuruExporter.export_to_benchling should produce Benchling format."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_benchling(sample_result)

        assert "name" in payload
        assert "sequence" in payload
        assert "annotations" in payload
        assert "customFields" in payload


# ────────────────────────────────────────────────────────────
# Convenience functions
# ────────────────────────────────────────────────────────────

class TestConvenienceFunctions:
    """Tests for module-level convenience functions."""

    def test_export_to_benchling(self, sample_result):
        """export_to_benchling should produce valid Benchling payload."""
        payload = export_to_benchling(sample_result)

        assert "name" in payload
        assert "sequence" in payload
        assert "annotations" in payload
        assert "customFields" in payload

    def test_export_to_labguru(self, sample_result):
        """export_to_labguru should produce valid LabGuru payload."""
        payload = export_to_labguru(sample_result)

        assert "item" in payload
        assert "data" in payload
        assert "tags" in payload
        assert "custom_fields" in payload


# ────────────────────────────────────────────────────────────
# Edge cases
# ────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_empty_protein(self):
        """Export with empty protein should not crash."""
        result = OptimizationResult(
            sequence="ATGGTTTCTAAAGGTGAA",
            gc_content=0.333,
            cai=0.78,
            protein="",
        )
        payload = export_to_benchling(result)
        assert payload["sequence"] == "ATGGTTTCTAAAGGTGAA"

        payload_lg = export_to_labguru(result)
        assert payload_lg["data"] == "ATGGTTTCTAAAGGTGAA"

    def test_long_sequence(self):
        """Export with a long sequence should work."""
        long_seq = "ATG" + "GCT" * 1000 + "TAA"  # ~3009 bp
        result = OptimizationResult(
            sequence=long_seq,
            gc_content=0.5,
            cai=0.9,
            protein="MA" + "A" * 999 + "*",  # approximate
        )
        # Just check it does not crash
        payload = export_to_benchling(result)
        assert len(payload["sequence"]) == len(long_seq)

    def test_multiple_submissions_unique_ids(self, sample_result):
        """Multiple submissions should produce unique design IDs."""
        exporter = BenchlingExporter()
        ids = set()
        for _ in range(10):
            design_id = exporter.submit_design(sample_result, "proj_1")
            ids.add(design_id)
        assert len(ids) == 10  # All unique

    def test_organism_name_underscores(self, sample_result):
        """Organism names with underscores should be handled properly."""
        exporter = LabGuruExporter()
        payload = exporter.export_to_labguru(sample_result)

        # LabGuru description should have spaces for display
        desc = payload["item"]["description"]
        assert "BioCompiler" in desc

    def test_benchling_annotation_colors(self, sample_result):
        """Annotations should have valid hex color strings."""
        exporter = BenchlingExporter()
        payload = exporter.export_to_benchling(sample_result)

        for ann in payload["annotations"]:
            assert "color" in ann
            assert ann["color"].startswith("#")
