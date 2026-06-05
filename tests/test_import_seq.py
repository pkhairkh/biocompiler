"""
Tests for biocompiler.import_seq -- FASTA, GenBank, and auto-detect import.

Covers:
  1. FASTA import (single record, multi-FASTA, edge cases)
  2. GenBank import (full record, features, exon boundaries, protein)
  3. Auto-detect format (import_sequence)
  4. Invalid input handling (bad DNA, missing headers, bad format)
  5. Format detection heuristics (_looks_like_path, _resolve_input)
  6. Helper functions (_parse_exon_boundaries, _clean_qualifier_value)
"""

from __future__ import annotations

import os
import textwrap

import pytest

from biocompiler.import_seq import (
    import_fasta,
    import_genbank,
    import_sequence,
    _parse_exon_boundaries,
    _clean_qualifier_value,
    _looks_like_path,
    _resolve_input,
)
from biocompiler.exceptions import FileFormatError, InvalidSequenceError


# -- Test data constants ---------------------------------------------------

SAMPLE_FASTA_SINGLE = ">seq1 test sequence\nATGCGTACGT\nATGCGTACGT\n"

SAMPLE_FASTA_MULTI = (
    ">seq1 first sequence\n"
    "ATGCGTACGT\n"
    ">seq2 second sequence|organism=Homo_sapiens\n"
    "GCTAGCTAGC\n"
    ">seq3 third one\n"
    "AATTCCGG\n"
)

SAMPLE_FASTA_WINDOWS = ">seq1\r\nATGCGTACGT\r\n"

SAMPLE_FASTA_EMPTY_LINES = ">seq1 test\n\nATGCGTACGT\n\nATGCGTACGT\n"

SAMPLE_GENBANK = (
    "LOCUS       HBB_HUMAN      626 bp    DNA     linear   PRI 15-JAN-2024\n"
    "DEFINITION  Homo sapiens hemoglobin subunit beta (HBB), mRNA.\n"
    "\n"
    "ACCESSION   U01317\n"
    "VERSION     U01317.1\n"
    "SOURCE      Homo sapiens (human)\n"
    "  ORGANISM  Homo sapiens\n"
    "            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata;\n"
    "            Mammalia; Eutheria; Euarchontoglires; Primates;\n"
    "            Catarrhini; Hominidae; Homo.\n"
    "FEATURES             Location/Qualifiers\n"
    "     source          1..626\n"
    "                     /organism=\"Homo sapiens\"\n"
    "                     /mol_type=\"mRNA\"\n"
    "     gene            51..494\n"
    "                     /gene=\"HBB\"\n"
    "     CDS             join(51..92,273..494)\n"
    "                     /gene=\"HBB\"\n"
    "                     /note=\"hemoglobin beta\"\n"
    '                     /translation="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"\n'
    "ORIGIN\n"
    "        1 atggtgcatc tgactcctga ggagaagtct gcggtaccct cttctgcatc tttcatacgg\n"
    "       61 agaagagcca tggtgcatct gactcctgag gagaagtctg ccgtactgcc ctgtggggca\n"
    "      121 aggtgaacgt ggatgaagtt ggtggtgagg ccctgggcag gctgctggtg gtctaccctt\n"
    "      181 ggacccagag gttctttgat ccaacctggc cgcgcttggt ggtctaccct tggacccaga\n"
    "      241 ggttctttga accaacttgg ccgcgcttgg tggtctaccc ttggacccag agggttcttt\n"
    "      301 gaccgacctg gccgcgcttg gtggtctacc cttggaccca gagggttctt tgaccgacct\n"
    "      361 ggccgcgctt ggtggtctac ccttggaccc agagggttct ttgaccgacc tggccgcgct\n"
    "      421 tggtggtcta cccttggacc cagagggttc tttgaccgac ctggccgcgc ttggtggtct\n"
    "      481 acccttggac ccagagggtt ctttgaccga cctggccgcg cttggtggtc tacccttgga\n"
    "      541 cccagagggt tctttgaccg acctggccgc gcttggtggt ctacccttgg acccagaggg\n"
    "      601 ttctttgacc gacctggccg cgcttggtgg\n"
    "//\n"
)

SAMPLE_GENBANK_NO_FEATURES = (
    "LOCUS       SIMPLE        18 bp    DNA     linear   PRI 15-JAN-2024\n"
    "DEFINITION  Simple test sequence.\n"
    "\n"
    "SOURCE      Homo sapiens\n"
    "  ORGANISM  Homo sapiens\n"
    "ORIGIN\n"
    "        1 atgcatgcat gcatgcatg\n"
    "//\n"
)

SAMPLE_GENBANK_COMPLEMENT = (
    "LOCUS       TEST_COMPL    100 bp    DNA     linear   PRI 15-JAN-2024\n"
    "DEFINITION  Complement test.\n"
    "\n"
    "SOURCE      Escherichia coli\n"
    "  ORGANISM  Escherichia coli\n"
    "FEATURES             Location/Qualifiers\n"
    "     source          1..100\n"
    "                     /organism=\"Escherichia coli\"\n"
    "     CDS             complement(join(10..50,70..95))\n"
    "                     /gene=\"testGene\"\n"
    '                     /translation="MTESTPROTEIN"\n'
    "ORIGIN\n"
    "        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat\n"
    "       51 gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc\n"
    "//\n"
)


# ======================================================================
# 1. FASTA Import Tests
# ======================================================================

class TestImportFasta:
    """Tests for import_fasta()."""

    def test_single_record(self):
        """Parse a single FASTA record from text."""
        records = import_fasta(SAMPLE_FASTA_SINGLE)
        assert len(records) == 1
        rec = records[0]
        assert rec["id"] == "seq1"
        assert rec["description"] == "test sequence"
        assert rec["sequence"] == "ATGCGTACGTATGCGTACGT"
        assert rec["gc_content"] > 0.0

    def test_multi_fasta(self):
        """Parse multiple FASTA records from a multi-FASTA string."""
        records = import_fasta(SAMPLE_FASTA_MULTI)
        assert len(records) == 3
        assert records[0]["id"] == "seq1"
        assert records[1]["id"] == "seq2"
        assert records[2]["id"] == "seq3"

    def test_organism_extraction(self):
        """Extract organism from pipe-delimited FASTA header."""
        records = import_fasta(SAMPLE_FASTA_MULTI)
        rec2 = records[1]
        assert rec2["organism"] == "Homo_sapiens"

    def test_no_organism_in_header(self):
        """Organism is empty string when not present in header."""
        records = import_fasta(SAMPLE_FASTA_SINGLE)
        assert records[0]["organism"] == ""

    def test_gc_content_computation(self):
        """GC content is computed correctly on import."""
        # GCTAGCTAGC: G=3, C=3, total=10, GC=0.6
        records = import_fasta(SAMPLE_FASTA_MULTI)
        rec2 = records[1]
        assert rec2["sequence"] == "GCTAGCTAGC"
        assert abs(rec2["gc_content"] - 0.6) < 0.01

    def test_windows_line_endings(self):
        """FASTA with Windows \\r\\n line endings is parsed correctly."""
        records = import_fasta(SAMPLE_FASTA_WINDOWS)
        assert len(records) == 1
        assert records[0]["id"] == "seq1"
        assert records[0]["sequence"] == "ATGCGTACGT"

    def test_empty_lines_between_blocks(self):
        """Empty lines between sequence blocks are handled gracefully."""
        records = import_fasta(SAMPLE_FASTA_EMPTY_LINES)
        assert len(records) == 1
        assert records[0]["sequence"] == "ATGCGTACGTATGCGTACGT"

    def test_from_file(self, tmp_path):
        """Parse FASTA from an actual file on disk."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(SAMPLE_FASTA_SINGLE)
        records = import_fasta(str(fasta_file))
        assert len(records) == 1
        assert records[0]["id"] == "seq1"

    def test_lowercase_input_uppered(self):
        """Lowercase DNA characters in FASTA are converted to uppercase."""
        text = ">seq1 lower\natgc\n"
        records = import_fasta(text)
        assert records[0]["sequence"] == "ATGC"

    def test_fasta_with_n_bases(self):
        """FASTA with N bases is accepted (N is valid in DNA)."""
        text = ">seq1 with N\nATGCNATGCN\n"
        records = import_fasta(text)
        assert "N" in records[0]["sequence"]

    def test_header_only_id(self):
        """FASTA header with only an ID and no description."""
        text = ">seq1\nATGC\n"
        records = import_fasta(text)
        assert records[0]["id"] == "seq1"
        assert records[0]["description"] == ""

    def test_pipe_delimited_header(self):
        """FASTA header with pipe-delimited metadata (e.g. NCBI style)."""
        text = ">gi|12345|gb|U01317.1| HBB mRNA\nATGC\n"
        records = import_fasta(text)
        assert records[0]["id"] == "gi|12345|gb|U01317.1|"


class TestImportFastaInvalid:
    """Tests for invalid FASTA input handling."""

    def test_no_header_raises(self):
        """Text without a '>' header raises FileFormatError."""
        with pytest.raises(FileFormatError, match="FASTA"):
            import_fasta("ATGCATGC")

    def test_sequence_before_header_raises(self):
        """Sequence data before any header raises FileFormatError."""
        with pytest.raises(FileFormatError, match="Sequence data found before any header"):
            import_fasta("ATGC\n>seq1\nATGC")

    def test_empty_input_raises(self):
        """Empty string raises FileFormatError (no records found)."""
        with pytest.raises(FileFormatError, match="No FASTA records"):
            import_fasta("")

    def test_invalid_dna_in_fasta_raises(self):
        """FASTA with non-DNA characters raises FileFormatError."""
        text = ">seq1 bad\nATGCXYZ\n"
        with pytest.raises(FileFormatError, match="Invalid DNA"):
            import_fasta(text)

    def test_only_header_no_sequence(self):
        """Header with no sequence data produces record with empty sequence."""
        text = ">seq1 empty\n"
        records = import_fasta(text)
        assert len(records) == 1
        assert records[0]["sequence"] == ""


# ======================================================================
# 2. GenBank Import Tests
# ======================================================================

class TestImportGenbank:
    """Tests for import_genbank()."""

    def test_basic_genbank_parsing(self):
        """Parse a full GenBank record with LOCUS, DEFINITION, FEATURES, ORIGIN."""
        result = import_genbank(SAMPLE_GENBANK)
        assert result["locus"] == "HBB_HUMAN"
        assert "hemoglobin" in result["definition"].lower() or "HBB" in result["definition"]
        assert result["organism"] == "Homo sapiens"
        assert len(result["sequence"]) > 0
        assert result["length"] == len(result["sequence"])

    def test_gc_content(self):
        """GC content is computed for the GenBank ORIGIN sequence."""
        result = import_genbank(SAMPLE_GENBANK)
        assert 0.0 <= result["gc_content"] <= 1.0

    def test_features_extracted(self):
        """Features are extracted from the FEATURES section."""
        result = import_genbank(SAMPLE_GENBANK)
        assert len(result["features"]) > 0
        feature_types = [f["type"] for f in result["features"]]
        assert "source" in feature_types
        assert "gene" in feature_types
        assert "CDS" in feature_types

    def test_gene_name_extraction(self):
        """Gene name is extracted from /gene qualifier."""
        result = import_genbank(SAMPLE_GENBANK)
        assert result["gene_name"] == "HBB"

    def test_exon_boundaries_from_join(self):
        """Exon boundaries are parsed from CDS join() locations."""
        result = import_genbank(SAMPLE_GENBANK)
        boundaries = result["exon_boundaries"]
        assert len(boundaries) == 2
        # join(51..92,273..494) -> [(50, 92), (272, 494)]
        assert boundaries[0] == (50, 92)
        assert boundaries[1] == (272, 494)

    def test_protein_from_translation(self):
        """Protein translation is extracted from /translation qualifier."""
        result = import_genbank(SAMPLE_GENBANK)
        assert result["protein"] != ""
        assert result["protein"].startswith("MVH")

    def test_cds_qualifiers(self):
        """Qualifiers like /gene, /note, /translation are parsed."""
        result = import_genbank(SAMPLE_GENBANK)
        cds = None
        for f in result["features"]:
            if f["type"] == "CDS":
                cds = f
                break
        assert cds is not None
        assert cds["qualifiers"]["gene"] == "HBB"
        assert "hemoglobin" in cds["qualifiers"].get("note", "")

    def test_feature_strand_forward(self):
        """Forward-strand features have strand '+'."""
        result = import_genbank(SAMPLE_GENBANK)
        cds = next(f for f in result["features"] if f["type"] == "CDS")
        assert cds["strand"] == "+"

    def test_feature_strand_complement(self):
        """Complement features have strand '-'."""
        result = import_genbank(SAMPLE_GENBANK_COMPLEMENT)
        cds = next(f for f in result["features"] if f["type"] == "CDS")
        assert cds["strand"] == "-"

    def test_genbank_from_file(self, tmp_path):
        """Parse GenBank from an actual file on disk."""
        gb_file = tmp_path / "test.gb"
        gb_file.write_text(SAMPLE_GENBANK)
        result = import_genbank(str(gb_file))
        assert result["locus"] == "HBB_HUMAN"

    def test_genbank_no_features(self):
        """GenBank record with no FEATURES section parses cleanly."""
        result = import_genbank(SAMPLE_GENBANK_NO_FEATURES)
        assert result["locus"] == "SIMPLE"
        assert result["sequence"] != ""
        assert result["features"] == []
        assert result["gene_name"] == ""
        assert result["exon_boundaries"] == []
        assert result["protein"] == ""

    def test_complement_join_exon_boundaries(self):
        """Exon boundaries from complement(join(...)) are extracted correctly."""
        result = import_genbank(SAMPLE_GENBANK_COMPLEMENT)
        boundaries = result["exon_boundaries"]
        assert len(boundaries) == 2
        # complement(join(10..50,70..95)) -> [(9, 50), (69, 95)]
        assert boundaries[0] == (9, 50)
        assert boundaries[1] == (69, 95)

    def test_complement_protein_extraction(self):
        """Protein is extracted from complement CDS."""
        result = import_genbank(SAMPLE_GENBANK_COMPLEMENT)
        assert result["protein"] == "MTESTPROTEIN"

    def test_windows_line_endings_genbank(self):
        """GenBank with Windows line endings is parsed correctly."""
        text = SAMPLE_GENBANK_NO_FEATURES.replace("\n", "\r\n")
        result = import_genbank(text)
        assert result["locus"] == "SIMPLE"


class TestImportGenbankInvalid:
    """Tests for invalid GenBank input handling."""

    def test_non_acgtn_chars_stripped_in_origin(self):
        """Non-ACGTN characters in ORIGIN section are stripped during extraction.

        The _extract_origin_sequence function removes all non-DNA characters
        before validation, so the remaining valid bases pass validation.
        This is by design: GenBank files may contain spacing and numbering
        that get filtered out.
        """
        gb = (
            "LOCUS       BAD_SEQ       10 bp    DNA     linear   PRI 15-JAN-2024\n"
            "DEFINITION  Bad sequence.\n"
            "\n"
            "SOURCE      Test\n"
            "  ORGANISM  Test\n"
            "ORIGIN\n"
            "        1 atgcxyzqwe\n"
            "//\n"
        )
        result = import_genbank(gb)
        # x,y,z,q,w,e are stripped; only ATGC remains
        assert result["sequence"] == "ATGC"

    def test_missing_origin_returns_empty_sequence(self):
        """GenBank without an ORIGIN section returns empty sequence."""
        no_origin = (
            "LOCUS       NO_ORIGIN     0 bp    DNA     linear   PRI 15-JAN-2024\n"
            "DEFINITION  No origin.\n"
            "\n"
            "SOURCE      Test\n"
            "  ORGANISM  Test\n"
            "//\n"
        )
        result = import_genbank(no_origin)
        assert result["sequence"] == ""
        assert result["gc_content"] == 0.0


# ======================================================================
# 3. Auto-detect Format Tests (import_sequence)
# ======================================================================

class TestImportSequenceAutoDetect:
    """Tests for import_sequence() auto-detection."""

    def test_detect_fasta(self):
        """Auto-detect FASTA format (starts with >)."""
        result = import_sequence(SAMPLE_FASTA_SINGLE)
        assert result["format"] == "fasta"
        assert result["sequence"] != ""
        assert "records" in result
        assert result["gc_content"] > 0.0

    def test_detect_genbank(self):
        """Auto-detect GenBank format (starts with LOCUS)."""
        result = import_sequence(SAMPLE_GENBANK_NO_FEATURES)
        assert result["format"] == "genbank"
        assert result["sequence"] != ""

    def test_detect_plain_text(self):
        """Auto-detect plain DNA text (no FASTA or GenBank markers)."""
        result = import_sequence("ATGCGTACGTATGCGTACGT")
        assert result["format"] == "plain"
        assert result["sequence"] == "ATGCGTACGTATGCGTACGT"

    def test_plain_text_gc_content(self):
        """GC content is computed for plain text input."""
        # GGCC: all GC, gc_content=1.0
        result = import_sequence("GGCC")
        assert result["format"] == "plain"
        assert abs(result["gc_content"] - 1.0) < 0.01

    def test_fasta_result_keys(self):
        """FASTA auto-detect result has all expected keys."""
        result = import_sequence(SAMPLE_FASTA_SINGLE)
        expected_keys = {
            "format", "records", "sequence", "gc_content",
            "organism", "gene_name", "exon_boundaries", "protein",
            "locus", "definition", "features",
        }
        assert set(result.keys()) == expected_keys

    def test_genbank_result_keys(self):
        """GenBank auto-detect result has all expected keys."""
        result = import_sequence(SAMPLE_GENBANK)
        expected_keys = {
            "format", "records", "sequence", "gc_content",
            "organism", "gene_name", "exon_boundaries", "protein",
            "locus", "definition", "features",
        }
        assert set(result.keys()) == expected_keys

    def test_plain_result_keys(self):
        """Plain text auto-detect result has all expected keys."""
        result = import_sequence("ATGC")
        expected_keys = {
            "format", "records", "sequence", "gc_content",
            "organism", "gene_name", "exon_boundaries", "protein",
            "locus", "definition", "features",
        }
        assert set(result.keys()) == expected_keys

    def test_genbank_exon_boundaries_propagated(self):
        """GenBank exon boundaries are propagated through import_sequence."""
        result = import_sequence(SAMPLE_GENBANK)
        assert len(result["exon_boundaries"]) == 2

    def test_genbank_protein_propagated(self):
        """GenBank protein is propagated through import_sequence."""
        result = import_sequence(SAMPLE_GENBANK)
        assert result["protein"] != ""

    def test_genbank_gene_name_propagated(self):
        """GenBank gene_name is propagated through import_sequence."""
        result = import_sequence(SAMPLE_GENBANK)
        assert result["gene_name"] == "HBB"

    def test_fasta_locus_from_id(self):
        """FASTA auto-detect sets locus from the sequence ID."""
        result = import_sequence(SAMPLE_FASTA_SINGLE)
        assert result["locus"] == "seq1"

    def test_fasta_definition_from_header(self):
        """FASTA auto-detect sets definition from the header description."""
        result = import_sequence(SAMPLE_FASTA_SINGLE)
        assert result["definition"] == "test sequence"


class TestImportSequenceInvalid:
    """Invalid input handling for import_sequence()."""

    def test_invalid_plain_dna_raises(self):
        """Plain text with invalid DNA characters raises FileFormatError."""
        with pytest.raises(FileFormatError, match="Invalid DNA"):
            import_sequence("ATGCXYZQQQ")

    def test_empty_input_raises(self):
        """Empty string cannot be auto-detected and raises FileFormatError."""
        with pytest.raises(FileFormatError):
            import_sequence("")

    def test_whitespace_only_raises(self):
        """Whitespace-only input raises FileFormatError."""
        with pytest.raises(FileFormatError):
            import_sequence("   \n  \t  ")


# ======================================================================
# 4. Helper Function Tests
# ======================================================================

class TestParseExonBoundaries:
    """Tests for _parse_exon_boundaries()."""

    def test_simple_range(self):
        """Simple range location: 1..495 -> [(0, 495)]."""
        result = _parse_exon_boundaries("1..495")
        assert result == [(0, 495)]

    def test_join_two_exons(self):
        """join(1..92,273..495) -> [(0, 92), (272, 495)]."""
        result = _parse_exon_boundaries("join(1..92,273..495)")
        assert result == [(0, 92), (272, 495)]

    def test_join_three_exons(self):
        """join(1..92,273..495,600..900) -> 3 boundaries."""
        result = _parse_exon_boundaries("join(1..92,273..495,600..900)")
        assert len(result) == 3
        assert result[0] == (0, 92)
        assert result[1] == (272, 495)
        assert result[2] == (599, 900)

    def test_complement_join(self):
        """complement(join(1..92,273..495)) -> boundaries extracted."""
        result = _parse_exon_boundaries("complement(join(1..92,273..495))")
        assert result == [(0, 92), (272, 495)]

    def test_single_position(self):
        """Single position (no range) -> single base boundary."""
        result = _parse_exon_boundaries("42")
        assert result == [(41, 42)]

    def test_empty_string(self):
        """Empty location string returns empty list."""
        result = _parse_exon_boundaries("")
        assert result == []

    def test_one_based_to_zero_based_conversion(self):
        """Verify 1-based inclusive -> 0-based exclusive conversion."""
        # 1..10 in GenBank = bases 1-10 inclusive = positions 0-9 in 0-based
        # 0-based exclusive end = 10
        result = _parse_exon_boundaries("1..10")
        assert result == [(0, 10)]


class TestCleanQualifierValue:
    """Tests for _clean_qualifier_value()."""

    def test_quoted_value(self):
        """Quoted value is unquoted."""
        assert _clean_qualifier_value('"Homo sapiens"') == "Homo sapiens"

    def test_unquoted_value(self):
        """Unquoted value is left as-is."""
        assert _clean_qualifier_value("some_value") == "some_value"

    def test_multiline_quoted(self):
        """Multi-line quoted value is joined and quotes removed."""
        raw = '"MGKL"\n"VTVL"'
        result = _clean_qualifier_value(raw)
        assert "MGKL" in result
        assert "VTVL" in result
        # Internal quote pairing should be removed
        assert '""' not in result

    def test_whitespace_cleanup(self):
        """Leading/trailing whitespace is stripped."""
        assert _clean_qualifier_value('  "value"  ') == "value"


class TestLooksLikePath:
    """Tests for _looks_like_path() heuristic."""

    def test_absolute_unix_path(self):
        """Absolute Unix path is detected."""
        assert _looks_like_path("/home/user/seq.fasta") is True

    def test_relative_path_with_slash(self):
        """Relative path with ./ is detected."""
        assert _looks_like_path("./data/seq.gb") is True

    def test_parent_path(self):
        """Path starting with .. is detected."""
        assert _looks_like_path("../data/seq.fasta") is True

    def test_file_extension(self):
        """Filename with .fasta extension is detected."""
        assert _looks_like_path("sequence.fasta") is True

    def test_genbank_extension(self):
        """Filename with .gb extension is detected."""
        assert _looks_like_path("sequence.gb") is True

    def test_fasta_text_not_path(self):
        """FASTA text content is not detected as a path."""
        assert _looks_like_path(">seq1\nATGC\n") is False

    def test_dna_sequence_not_path(self):
        """DNA sequence text is not detected as a path."""
        assert _looks_like_path("ATGCGTACGT") is False

    def test_empty_string_not_path(self):
        """Empty string is not a path."""
        assert _looks_like_path("") is False

    def test_windows_path(self):
        """Windows path is detected."""
        assert _looks_like_path("C:\\Users\\seq.fasta") is True


class TestResolveInput:
    """Tests for _resolve_input()."""

    def test_existing_file_is_read(self, tmp_path):
        """When input is a path to an existing file, its content is returned."""
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        result = _resolve_input(str(f), "test")
        assert result == "hello world"

    def test_nonexistent_path_treated_as_text(self):
        """When input looks like a path but file doesn't exist, it's treated as text."""
        # A path that doesn't exist -- _resolve_input returns the text as-is
        result = _resolve_input("/nonexistent/path.fasta", "test")
        # It returns the original string since the file doesn't exist
        assert result == "/nonexistent/path.fasta"

    def test_raw_text_returned_as_is(self):
        """Raw text content (not a path) is returned as-is."""
        result = _resolve_input(">seq1\nATGC\n", "test")
        assert result == ">seq1\nATGC\n"


# ======================================================================
# 5. Integration-style tests
# ======================================================================

class TestImportSequenceFileIO:
    """Tests that import from actual files on disk works end-to-end."""

    def test_fasta_file_roundtrip(self, tmp_path):
        """Write and re-read a FASTA file."""
        f = tmp_path / "test.fasta"
        f.write_text(SAMPLE_FASTA_SINGLE)
        result = import_sequence(str(f))
        assert result["format"] == "fasta"
        assert result["sequence"] == "ATGCGTACGTATGCGTACGT"

    def test_genbank_file_roundtrip(self, tmp_path):
        """Write and re-read a GenBank file."""
        f = tmp_path / "test.gb"
        f.write_text(SAMPLE_GENBANK)
        result = import_sequence(str(f))
        assert result["format"] == "genbank"
        assert result["locus"] == "HBB_HUMAN"

    def test_plain_text_file(self, tmp_path):
        """A .txt file with plain DNA sequence."""
        f = tmp_path / "seq.txt"
        f.write_text("ATGCGTACGT")
        result = import_sequence(str(f))
        # .txt extension triggers _looks_like_path -> reads the file
        # Content doesn't start with > or LOCUS, so detected as plain
        assert result["format"] == "plain"
        assert result["sequence"] == "ATGCGTACGT"

    def test_unreadable_file_raises(self, tmp_path):
        """A path that exists but cannot be read raises FileFormatError."""
        f = tmp_path / "noperm.fasta"
        f.write_text(SAMPLE_FASTA_SINGLE)
        # Remove read permission
        os.chmod(str(f), 0o000)
        try:
            with pytest.raises(FileFormatError, match="Cannot read file"):
                import_fasta(str(f))
        finally:
            # Restore permissions so tmp_path cleanup works
            os.chmod(str(f), 0o644)
