"""Tests for the parts module — Part dataclass and PartLibrary class."""

from __future__ import annotations

import json
import os
import tempfile

import pytest

from biocompiler.parts import Part, PartLibrary, DEFAULT_PARTS


# ---------------------------------------------------------------------------
# Part dataclass
# ---------------------------------------------------------------------------

class TestPart:
    """Tests for the Part dataclass."""

    def test_construction(self):
        """Can construct a valid Part."""
        p = Part(name="test_promoter", part_type="promoter", sequence="ATGC")
        assert p.name == "test_promoter"
        assert p.part_type == "promoter"
        assert p.sequence == "ATGC"

    def test_description_default(self):
        """Description defaults to empty string."""
        p = Part(name="x", part_type="rbs", sequence="ATGC")
        assert p.description == ""

    def test_metadata_default(self):
        """Metadata defaults to empty dict."""
        p = Part(name="x", part_type="rbs", sequence="ATGC")
        assert p.metadata == {}

    def test_sequence_uppercased(self):
        """Sequence is normalized to uppercase."""
        p = Part(name="x", part_type="rbs", sequence="atgc")
        assert p.sequence == "ATGC"

    def test_invalid_part_type_raises(self):
        """Invalid part_type raises ValueError."""
        with pytest.raises(ValueError, match="Invalid part_type"):
            Part(name="x", part_type="invalid_type", sequence="ATGC")

    def test_empty_name_raises(self):
        """Empty name raises ValueError."""
        with pytest.raises(ValueError, match="name"):
            Part(name="", part_type="promoter", sequence="ATGC")

    def test_empty_sequence_raises(self):
        """Empty sequence raises ValueError."""
        with pytest.raises(ValueError, match="sequence"):
            Part(name="x", part_type="promoter", sequence="")

    def test_valid_types(self):
        """All valid part types are accepted."""
        for pt in ("promoter", "cds", "terminator", "rbs", "linker"):
            p = Part(name=f"test_{pt}", part_type=pt, sequence="ATGC")
            assert p.part_type == pt

    def test_with_metadata(self):
        """Metadata is stored correctly."""
        p = Part(
            name="B0034",
            part_type="rbs",
            sequence="AAAGAGGAGATATACAT",
            description="Strong RBS",
            metadata={"organism": "E_coli", "strength": "strong"},
        )
        assert p.metadata["organism"] == "E_coli"
        assert p.metadata["strength"] == "strong"


# ---------------------------------------------------------------------------
# Default parts
# ---------------------------------------------------------------------------

class TestDefaultParts:
    """Tests for the built-in DEFAULT_PARTS list."""

    def test_not_empty(self):
        """DEFAULT_PARTS is not empty."""
        assert len(DEFAULT_PARTS) > 0

    def test_all_are_part_instances(self):
        """All default parts are Part instances."""
        for p in DEFAULT_PARTS:
            assert isinstance(p, Part)

    def test_contains_t7_promoter(self):
        """T7 promoter is in the default parts."""
        names = [p.name for p in DEFAULT_PARTS]
        assert "T7_promoter" in names

    def test_contains_lac_promoter(self):
        """Lac promoter is in the default parts."""
        names = [p.name for p in DEFAULT_PARTS]
        assert "lac_promoter" in names

    def test_contains_t7_terminator(self):
        """T7 terminator is in the default parts."""
        names = [p.name for p in DEFAULT_PARTS]
        assert "T7_terminator" in names

    def test_contains_rbs(self):
        """At least one RBS is in the default parts."""
        rbs_parts = [p for p in DEFAULT_PARTS if p.part_type == "rbs"]
        assert len(rbs_parts) >= 1

    def test_all_sequences_valid_dna(self):
        """All default part sequences contain only ACGT."""
        for p in DEFAULT_PARTS:
            for base in p.sequence:
                assert base in "ACGT", f"Invalid base {base!r} in {p.name}"

    def test_unique_names(self):
        """All default part names are unique."""
        names = [p.name for p in DEFAULT_PARTS]
        assert len(names) == len(set(names))

    def test_covers_multiple_types(self):
        """Default parts cover at least 3 different types."""
        types = {p.part_type for p in DEFAULT_PARTS}
        assert len(types) >= 3


# ---------------------------------------------------------------------------
# PartLibrary
# ---------------------------------------------------------------------------

class TestPartLibrary:
    """Tests for the PartLibrary class."""

    def test_default_construction(self):
        """Default construction loads built-in parts."""
        lib = PartLibrary()
        assert len(lib) > 0

    def test_get_existing_part(self):
        """get() returns a known part."""
        lib = PartLibrary()
        part = lib.get("T7_promoter")
        assert part.name == "T7_promoter"
        assert part.part_type == "promoter"

    def test_get_nonexistent_raises(self):
        """get() raises KeyError for unknown parts."""
        lib = PartLibrary()
        with pytest.raises(KeyError, match="not found"):
            lib.get("NONEXISTENT_PART")

    def test_search_by_type(self):
        """search() filters by part_type."""
        lib = PartLibrary()
        promoters = lib.search("promoter")
        assert len(promoters) >= 1
        for p in promoters:
            assert p.part_type == "promoter"

    def test_search_by_type_and_organism(self):
        """search() filters by part_type and organism."""
        lib = PartLibrary()
        ecoli_promoters = lib.search("promoter", organism="E_coli")
        assert len(ecoli_promoters) >= 1
        for p in ecoli_promoters:
            assert p.part_type == "promoter"
            assert p.metadata.get("organism", "").lower() == "e_coli"

    def test_search_no_match(self):
        """search() returns empty list when no matches."""
        lib = PartLibrary()
        results = lib.search("cds", organism="NONEXISTENT")
        assert results == []

    def test_add_part(self):
        """add() adds a new part to the library."""
        lib = PartLibrary()
        initial_count = len(lib)
        lib.add(Part(name="custom_cds", part_type="cds", sequence="ATGCGT"))
        assert len(lib) == initial_count + 1
        assert lib.get("custom_cds").sequence == "ATGCGT"

    def test_add_replaces_existing(self):
        """add() replaces a part with the same name."""
        lib = PartLibrary()
        lib.add(Part(name="T7_promoter", part_type="promoter", sequence="AAA"))
        assert lib.get("T7_promoter").sequence == "AAA"

    def test_list_parts(self):
        """list_parts() returns sorted list of part names."""
        lib = PartLibrary()
        names = lib.list_parts()
        assert names == sorted(names)
        assert "T7_promoter" in names

    def test_contains(self):
        """__contains__ checks for part existence."""
        lib = PartLibrary()
        assert "T7_promoter" in lib
        assert "NONEXISTENT" not in lib

    def test_repr(self):
        """repr includes part count."""
        lib = PartLibrary()
        r = repr(lib)
        assert "parts" in r

    def test_len(self):
        """len() returns number of parts."""
        lib = PartLibrary()
        assert len(lib) == len(lib._parts)


# ---------------------------------------------------------------------------
# PartLibrary — file loading
# ---------------------------------------------------------------------------

class TestPartLibraryFileLoading:
    """Tests for loading parts from JSON and YAML files."""

    def test_load_json_file(self):
        """Can load parts from a JSON file."""
        parts_data = [
            {
                "name": "test_json_prom",
                "part_type": "promoter",
                "sequence": "ATGCATGC",
                "description": "Test promoter from JSON",
                "metadata": {"organism": "E_coli"},
            }
        ]
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump(parts_data, f)
            tmp_path = f.name

        try:
            lib = PartLibrary(library_path=tmp_path)
            part = lib.get("test_json_prom")
            assert part.part_type == "promoter"
            assert part.sequence == "ATGCATGC"
        finally:
            os.unlink(tmp_path)

    def test_load_json_overlays_defaults(self):
        """JSON file parts overlay but don't replace all defaults."""
        parts_data = [
            {
                "name": "my_custom",
                "part_type": "linker",
                "sequence": "GGGGGG",
            }
        ]
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump(parts_data, f)
            tmp_path = f.name

        try:
            lib = PartLibrary(library_path=tmp_path)
            # Custom part is loaded
            assert "my_custom" in lib
            # Default parts are still there
            assert "T7_promoter" in lib
        finally:
            os.unlink(tmp_path)

    def test_file_not_found_raises(self):
        """Non-existent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            PartLibrary(library_path="/nonexistent/path/parts.json")

    def test_unsupported_format_raises(self):
        """Unsupported file format raises ValueError."""
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            f.write(b"not valid")
            tmp_path = f.name
        try:
            with pytest.raises(ValueError, match="Unsupported"):
                PartLibrary(library_path=tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_invalid_json_raises(self):
        """Malformed JSON raises appropriate error."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            f.write("not valid json {{{")
            tmp_path = f.name
        try:
            with pytest.raises(json.JSONDecodeError):
                PartLibrary(library_path=tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_missing_required_field_raises(self):
        """Part entry missing required fields raises ValueError."""
        parts_data = [{"name": "no_type_or_seq"}]  # Missing part_type and sequence
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump(parts_data, f)
            tmp_path = f.name
        try:
            with pytest.raises(ValueError, match="name.*part_type.*sequence"):
                PartLibrary(library_path=tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_non_list_data_raises(self):
        """Top-level data that isn't a list raises ValueError."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump({"not": "a list"}, f)
            tmp_path = f.name
        try:
            with pytest.raises(ValueError, match="list"):
                PartLibrary(library_path=tmp_path)
        finally:
            os.unlink(tmp_path)
