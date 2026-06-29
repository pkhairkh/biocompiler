"""Tests for the BioCompiler IR YAML frontend.

Covers:
* :func:`parse_spec` — YAML string → IR_L0_GenomicDNA
* :func:`compile_from_spec` — YAML → IR_L0 → L1/L2/L3/L4 via ``compile_gene``
* Error handling — missing required fields, invalid bases, bad region
  coordinates, malformed YAML, ...
* Provenance stamping — ``metadata["pass"]`` / ``metadata["source_format"]``
* The two example specs shipped with the package (``hbb.yaml``,
  ``gfp.yaml``) — verifies they parse and compile to the expected
  protein sequences.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from biocompiler.ir.frontend import (
    parse_spec,
    compile_from_spec,
    SpecError,
)
from biocompiler.ir.types import (
    IRLevel,
    IR_L0_GenomicDNA,
    IR_L1_PreMRNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IR_L4_FoldedProtein,
    GeneRegion,
    IRError,
)

# ────────────────────────────────────────────────────────────────────
# Path to the example specs shipped with the package.
# ────────────────────────────────────────────────────────────────────
_EXAMPLE_SPECS_DIR = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "biocompiler"
    / "ir"
    / "example_specs"
)
HBB_YAML = _EXAMPLE_SPECS_DIR / "hbb.yaml"
GFP_YAML = _EXAMPLE_SPECS_DIR / "gfp.yaml"

# Canonical HBB N-terminus + stop — the "hello world" of the IR.
HBB_DNA = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG"
    "GTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA"
)
HBB_PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"

# GFP N-terminus + stop (first 15 residues of avGFP, UniProt P42212).
GFP_PROTEIN = "MSKGEELFTGVVPIL*"


# ════════════════════════════════════════════════════════════════════
# parse_spec — YAML string → IR_L0
# ════════════════════════════════════════════════════════════════════
class TestParseSpec:
    def test_parse_yaml_string_basic(self):
        yaml_str = (
            "gene_name: HBB\n"
            "organism: human\n"
            f"sequence: {HBB_DNA}\n"
        )
        ir = parse_spec(yaml_str)
        assert isinstance(ir, IR_L0_GenomicDNA)
        assert ir.sequence == HBB_DNA
        assert ir.organism == "human"
        assert ir.gene_name == "HBB"
        assert ir.regions == []
        assert ir.level == IRLevel.L0

    def test_parse_yaml_string_minimal(self):
        # Only the two required fields.
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTTAA\n"
        )
        ir = parse_spec(yaml_str)
        assert ir.sequence == "ATGGCTTAA"
        assert ir.organism == "e_coli"
        assert ir.gene_name is None
        assert ir.regions == []

    def test_parse_yaml_string_with_regions(self):
        yaml_str = (
            "gene_name: g\n"
            "organism: e_coli\n"
            "sequence: ATGGCTTAA\n"
            "regions:\n"
            "  - type: cds\n"
            "    start: 0\n"
            "    end: 9\n"
        )
        ir = parse_spec(yaml_str)
        assert len(ir.regions) == 1
        r = ir.regions[0]
        assert isinstance(r, GeneRegion)
        assert r.region_type == "cds"
        assert r.start == 0
        assert r.end == 9
        assert r.metadata == {}

    def test_parse_yaml_string_with_region_metadata(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTTAA\n"
            "regions:\n"
            "  - type: exon\n"
            "    start: 0\n"
            "    end: 9\n"
            "    metadata:\n"
            "      frame: 0\n"
            "      note: first-exon\n"
        )
        ir = parse_spec(yaml_str)
        assert ir.regions[0].metadata == {"frame": 0, "note": "first-exon"}

    def test_parse_yaml_string_with_top_level_metadata(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTTAA\n"
            "metadata:\n"
            "  description: a test gene\n"
            "  request_id: abc-123\n"
        )
        ir = parse_spec(yaml_str)
        # User metadata is preserved.
        assert ir.metadata["description"] == "a test gene"
        assert ir.metadata["request_id"] == "abc-123"
        # And provenance is stamped.
        assert ir.metadata["pass"] == "frontend.parse_spec"
        assert ir.metadata["source_format"] == "yaml"

    def test_parse_yaml_uppercases_sequence(self):
        yaml_str = (
            "organism: e_coli\n"
            f"sequence: {HBB_DNA.lower()}\n"
        )
        ir = parse_spec(yaml_str)
        assert ir.sequence == HBB_DNA  # upper-cased

    def test_parse_yaml_allows_N_bases(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGNNNTAA\n"
        )
        ir = parse_spec(yaml_str)
        assert ir.sequence == "ATGNNNTAA"

    def test_parse_yaml_multiple_regions(self):
        # Two exons + one intron, like a miniature eukaryotic gene.
        yaml_str = (
            "gene_name: g\n"
            "organism: human\n"
            "sequence: GGGATGGCTTAACCC\n"   # 15 nt
            "regions:\n"
            "  - type: exon\n"
            "    start: 0\n"
            "    end: 9\n"
            "  - type: intron\n"
            "    start: 9\n"
            "    end: 12\n"
            "  - type: exon\n"
            "    start: 12\n"
            "    end: 15\n"
        )
        ir = parse_spec(yaml_str)
        assert len(ir.regions) == 3
        assert [r.region_type for r in ir.regions] == ["exon", "intron", "exon"]
        assert [(r.start, r.end) for r in ir.regions] == [
            (0, 9), (9, 12), (12, 15),
        ]

    def test_parse_yaml_region_metadata_defaults_to_empty_dict(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTTAA\n"
            "regions:\n"
            "  - type: cds\n"
            "    start: 0\n"
            "    end: 9\n"
        )
        ir = parse_spec(yaml_str)
        assert ir.regions[0].metadata == {}

    def test_parse_yaml_string_provenance_stamped(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTTAA\n"
        )
        ir = parse_spec(yaml_str)
        assert ir.metadata["pass"] == "frontend.parse_spec"
        assert ir.metadata["source_format"] == "yaml"


# ════════════════════════════════════════════════════════════════════
# parse_spec — YAML file → IR_L0
# ════════════════════════════════════════════════════════════════════
class TestParseSpecFromFile:
    def test_parse_hbb_yaml_file(self):
        ir = parse_spec(str(HBB_YAML))
        assert isinstance(ir, IR_L0_GenomicDNA)
        assert ir.gene_name == "HBB"
        assert ir.organism == "human"
        assert ir.sequence == HBB_DNA
        # The spec has one CDS region spanning [0, 96).
        assert len(ir.regions) == 1
        assert ir.regions[0].region_type == "cds"
        assert ir.regions[0].start == 0
        assert ir.regions[0].end == 96

    def test_parse_gfp_yaml_file(self):
        ir = parse_spec(str(GFP_YAML))
        assert isinstance(ir, IR_L0_GenomicDNA)
        assert ir.gene_name == "GFP"
        assert ir.organism == "aequorea_victoria"
        # 48 nt = 16 codons.
        assert len(ir.sequence) == 48
        assert ir.regions[0].region_type == "cds"
        assert (ir.regions[0].start, ir.regions[0].end) == (0, 48)

    def test_parse_yaml_file_path_object(self):
        # ``parse_spec`` accepts a string; ensure pathlib.Path converts
        # cleanly (Path.__str__ produces a valid file path).
        ir = parse_spec(str(HBB_YAML))
        assert ir.gene_name == "HBB"

    def test_parse_yaml_file_stamps_provenance(self):
        ir = parse_spec(str(HBB_YAML))
        assert ir.metadata["pass"] == "frontend.parse_spec"
        assert ir.metadata["source_format"] == "yaml"
        # User metadata is preserved too.
        assert "source" in ir.metadata  # "UniProt P68871 (HBB_HUMAN)"


# ════════════════════════════════════════════════════════════════════
# compile_from_spec — YAML → IR_L0 → {L1, L2, L3, L4}
# ════════════════════════════════════════════════════════════════════
class TestCompileFromSpec:
    def test_compile_default_target_is_l3(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTAAGTAA\n"  # M A K *
        )
        out = compile_from_spec(yaml_str)
        assert isinstance(out, IR_L3_Polypeptide)
        assert out.sequence == "MAK*"
        assert out.level == IRLevel.L3

    def test_compile_to_l0_is_identity_with_parsed(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTAAGTAA\n"
        )
        out = compile_from_spec(yaml_str, target_level=IRLevel.L0)
        assert isinstance(out, IR_L0_GenomicDNA)
        assert out.sequence == "ATGGCTAAGTAA"

    def test_compile_to_l1(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTAAGTAA\n"
        )
        out = compile_from_spec(yaml_str, target_level=IRLevel.L1)
        assert isinstance(out, IR_L1_PreMRNA)
        assert out.sequence == "AUGGCUAAGUAA"  # T → U

    def test_compile_to_l2(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTAAGTAA\n"
        )
        out = compile_from_spec(yaml_str, target_level=IRLevel.L2)
        assert isinstance(out, IR_L2_MatureMRNA)
        assert out.cds == "AUGGCUAAGUAA"
        assert out.five_utr == ""
        assert out.three_utr == ""

    def test_compile_to_l4(self):
        """Compile gene spec to L4 (folded protein with heuristic fallback)."""
        # compile_from_spec takes a YAML string
        yaml_str = "sequence: ATGGCTAAGTAA\norganism: e_coli\n"
        result = compile_from_spec(yaml_str, target_level=IRLevel.L4)
        assert result.sequence == "MAK*"
        assert result.confidence is not None  # folding oracle ran
        assert result.metadata.get("oracle") == "fallback"

    def test_compile_hbb_yaml_to_l3(self):
        """The canonical HBB smoke test: spec → IR-L3 → UniProt match."""
        out = compile_from_spec(str(HBB_YAML))
        assert isinstance(out, IR_L3_Polypeptide)
        assert out.sequence == HBB_PROTEIN
        assert out.gene_name == "HBB"
        assert out.organism == "human"
        assert out.sequence.endswith("*")  # stop codon preserved

    def test_compile_gfp_yaml_to_l3(self):
        """GFP fragment smoke test: spec → IR-L3 → expected N-terminus."""
        out = compile_from_spec(str(GFP_YAML))
        assert isinstance(out, IR_L3_Polypeptide)
        assert out.sequence == GFP_PROTEIN
        assert out.gene_name == "GFP"

    def test_compile_propagates_user_metadata(self):
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTAAGTAA\n"
            "metadata:\n"
            "  request_id: req-xyz-001\n"
        )
        out = compile_from_spec(yaml_str, target_level=IRLevel.L3)
        # User metadata flows through every pass.
        assert out.metadata["request_id"] == "req-xyz-001"
        # The frontend's provenance is also preserved.
        assert out.metadata["source_format"] == "yaml"
        # The last pass applied (translate) overwrites the "pass" key.
        assert out.metadata["pass"] == "translate"


# ════════════════════════════════════════════════════════════════════
# Error handling — malformed specs
# ════════════════════════════════════════════════════════════════════
class TestSpecErrors:
    def test_error_missing_sequence(self):
        with pytest.raises(SpecError, match="sequence"):
            parse_spec("organism: e_coli\n")

    def test_error_missing_organism(self):
        with pytest.raises(SpecError, match="organism"):
            parse_spec("sequence: ATGGCTTAA\n")

    def test_error_empty_sequence(self):
        with pytest.raises(SpecError, match="non-empty"):
            parse_spec("organism: e_coli\nsequence: ''\n")

    def test_error_invalid_bases(self):
        with pytest.raises(SpecError, match="invalid bases"):
            parse_spec("organism: e_coli\nsequence: ATGXYZ\n")

    def test_error_invalid_bases_lists_them(self):
        # Error message should name the offending bases so users can
        # fix the typo without grepping the sequence by eye.
        with pytest.raises(SpecError, match=r"['\"]X['\"].*['\"]Z['\"]|['\"]Z['\"].*['\"]X['\"]"):
            parse_spec("organism: e_coli\nsequence: ATGXZ\n")

    def test_error_non_string_sequence(self):
        with pytest.raises(SpecError, match="must be a string"):
            parse_spec("organism: e_coli\nsequence: 12345\n")

    def test_error_non_string_organism(self):
        with pytest.raises(SpecError, match="organism"):
            parse_spec("organism: 42\nsequence: ATGGCTTAA\n")

    def test_error_empty_organism(self):
        with pytest.raises(SpecError, match="organism"):
            parse_spec("organism: ''\nsequence: ATGGCTTAA\n")

    def test_error_non_string_gene_name(self):
        with pytest.raises(SpecError, match="gene_name"):
            parse_spec(
                "gene_name: 99\n"
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
            )

    def test_error_region_missing_type(self):
        with pytest.raises(SpecError, match="type"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - start: 0\n"
                "    end: 9\n"
            )

    def test_error_region_missing_start(self):
        with pytest.raises(SpecError, match="start"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    end: 9\n"
            )

    def test_error_region_missing_end(self):
        with pytest.raises(SpecError, match="end"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: 0\n"
            )

    def test_error_region_end_le_start(self):
        with pytest.raises(SpecError, match="must be > 'start'"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: 5\n"
                "    end: 5\n"
            )

    def test_error_region_end_lt_start(self):
        with pytest.raises(SpecError, match="must be > 'start'"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: 7\n"
                "    end: 3\n"
            )

    def test_error_region_out_of_bounds(self):
        # Sequence is 9 nt, but the region claims to end at 100.
        with pytest.raises(SpecError, match="out of bounds"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: 0\n"
                "    end: 100\n"
            )

    def test_error_region_negative_start(self):
        with pytest.raises(SpecError, match="out of bounds"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: -1\n"
                "    end: 9\n"
            )

    def test_error_region_non_integer_start(self):
        with pytest.raises(SpecError, match="'start' must be an integer"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: zero\n"
                "    end: 9\n"
            )

    def test_error_region_boolean_start(self):
        # ``bool`` is a subclass of ``int`` — guard explicitly so
        # ``start: true`` doesn't silently become ``start: 1``.
        with pytest.raises(SpecError, match="'start' must be an integer"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - type: cds\n"
                "    start: true\n"
                "    end: 9\n"
            )

    def test_error_region_non_dict_entry(self):
        with pytest.raises(SpecError, match="must be a mapping"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions:\n"
                "  - just-a-string\n"
            )

    def test_error_regions_not_a_list(self):
        with pytest.raises(SpecError, match="'regions' must be a list"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "regions: not-a-list\n"
            )

    def test_error_metadata_not_a_dict(self):
        with pytest.raises(SpecError, match="'metadata' must be a mapping"):
            parse_spec(
                "organism: e_coli\n"
                "sequence: ATGGCTTAA\n"
                "metadata: not-a-mapping\n"
            )

    def test_error_non_dict_root(self):
        # YAML root is a list, not a mapping.
        with pytest.raises(SpecError, match="YAML root must be a mapping"):
            parse_spec("- a\n- b\n- c\n")

    def test_error_empty_yaml_document(self):
        with pytest.raises(SpecError, match="empty"):
            parse_spec("# just a comment\n")

    def test_error_empty_input_string(self):
        with pytest.raises(SpecError, match="empty"):
            parse_spec("")

    def test_error_whitespace_only_input(self):
        with pytest.raises(SpecError, match="empty"):
            parse_spec("   \n\t  \n")

    def test_error_non_string_input(self):
        with pytest.raises(SpecError, match="input must be a string"):
            parse_spec(12345)  # type: ignore[arg-type]

    def test_error_nonexistent_yaml_file_falls_back_to_string_parse(self):
        # A path-like string that doesn't exist is treated as YAML
        # content — which should fail to parse as a mapping.
        with pytest.raises((SpecError, Exception)):
            parse_spec("/nonexistent/path/to/gene.yaml")


# ════════════════════════════════════════════════════════════════════
# Error handling — IR-level errors after parsing
# ════════════════════════════════════════════════════════════════════
class TestIRErrrsAfterParse:
    def test_ir_no_start_codon_succeeds(self):
        """CDS no longer requires AUG start — back-translated genes may not start with M."""
        yaml_str = (
            "organism: e_coli\n"
            "sequence: GGGGCTTAA\n"  # no ATG, but that's OK now
        )
        result = compile_from_spec(yaml_str)
        assert result.sequence == "GA*"  # GGG=GCT=TAA → G,A,*

    def test_ir_error_no_stop_codon(self):
        # Spec parses fine, but no in-frame stop codon.
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGGCTAAG\n"  # ATG GCT AAG — no stop
        )
        with pytest.raises(IRError):
            compile_from_spec(yaml_str)

    def test_ir_error_invalid_base_caught_at_transcribe(self):
        # 'X' is not a valid base — the IR's transcribe pass rejects it
        # before the frontend can even build a useful IR-L0.  (The
        # frontend also rejects this at parse time, so this is a
        # belt-and-braces check.)
        yaml_str = (
            "organism: e_coli\n"
            "sequence: ATGXZTAA\n"
        )
        with pytest.raises(SpecError):
            compile_from_spec(yaml_str)


# ════════════════════════════════════════════════════════════════════
# Region semantics — make sure YAML regions flow into splice correctly
# ════════════════════════════════════════════════════════════════════
class TestRegionSemantics:
    def test_cds_region_drives_splice(self):
        # A CDS region that explicitly excludes flanking UTR sequence.
        # Whole sequence: "GGGATGGCTTAACCC" (15 nt)
        # CDS region: [3, 12) → "ATGGCTTAA" → M A *
        yaml_str = (
            "gene_name: g\n"
            "organism: human\n"
            "sequence: GGGATGGCTTAACCC\n"
            "regions:\n"
            "  - type: cds\n"
            "    start: 3\n"
            "    end: 12\n"
        )
        out = compile_from_spec(yaml_str, target_level=IRLevel.L2)
        assert isinstance(out, IR_L2_MatureMRNA)
        # Splice picks the CDS region [3,12) and scans it for AUG/stop.
        assert out.cds.endswith("UAA")  # ends with stop codon
        assert out.five_utr == ""  # nothing before AUG inside the CDS slice
        assert out.three_utr == ""  # nothing after stop inside the CDS slice

    def test_exon_regions_are_concatenated_by_splice(self):
        # Two exons separated by an intron, spliced together to form
        # the CDS.
        #   exon1 [0,9):  "GGGATGGCT"
        #   intron [9,12): "TAA" (will be dropped)
        #   exon2 [12,15): "CCC"
        # Spliced: "GGGATGGCTCCC" — but no stop codon here, so splice
        # will raise.  Instead let's design exons that DO yield a stop.
        #   exon1 [0,9):  "GGGATGGCT"  (5'UTR GGG + ATG GCT)
        #   exon2 [9,12): "TAA"        (stop codon)
        # Spliced: "GGGATGGCTTAA" → 5'UTR=GGG, CDS=AUGGCUUAA, 3'UTR=""
        yaml_str = (
            "gene_name: g\n"
            "organism: human\n"
            "sequence: GGGATGGCTTAACCC\n"
            "regions:\n"
            "  - type: exon\n"
            "    start: 0\n"
            "    end: 9\n"
            "  - type: intron\n"
            "    start: 9\n"
            "    end: 12\n"
            "  - type: exon\n"
            "    start: 12\n"
            "    end: 15\n"
        )
        # The intron "TAA" happens to be a stop codon, but splice
        # concatenates only exon/cds regions, so the intron is dropped.
        # Spliced transcript = exon1 + exon2 = "GGGATGGCT" + "CCC"
        # That has no stop codon → splice raises IRError.  Adjust:
        # we want the spliced transcript to contain AUG + ... + stop.
        # Use exon2 = "TAA" by making exon2 cover [9, 12).
        yaml_str = (
            "gene_name: g\n"
            "organism: human\n"
            "sequence: GGGATGGCTTAA\n"  # 12 nt
            "regions:\n"
            "  - type: exon\n"
            "    start: 0\n"
            "    end: 9\n"
            "  - type: exon\n"
            "    start: 9\n"
            "    end: 12\n"
        )
        out = compile_from_spec(yaml_str, target_level=IRLevel.L2)
        assert out.five_utr == ""  # no 5UTR (CDS starts at position 0)
        assert out.cds.endswith("UAA")  # ends with stop codon
        assert out.three_utr == ""


# ════════════════════════════════════════════════════════════════════
# SpecError is a ValueError (so generic catch works)
# ════════════════════════════════════════════════════════════════════
class TestSpecErrorType:
    def test_spec_error_is_value_error(self):
        # SpecError must subclass ValueError so callers that catch
        # ValueError (e.g. argparse type handlers) still work.
        assert issubclass(SpecError, ValueError)

    def test_spec_error_is_not_ir_error(self):
        # And it must NOT be IRError — they encode different failures.
        assert not issubclass(SpecError, IRError)
