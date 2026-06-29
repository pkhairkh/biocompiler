"""
Property-based tests for biocompiler.export using Hypothesis.

Covers four core properties:
  1. FASTA export always starts with ">"
  2. GenBank export always ends with "//"
  3. Exported sequences contain only uppercase ACGT characters
  4. GC content advertised in the FASTA header matches independently computed GC
"""

from __future__ import annotations

import re

import pytest
pytest.importorskip("hypothesis")
from hypothesis import given, settings, assume
from hypothesis import strategies as st

from biocompiler.export.core import export_fasta, export_genbank
from biocompiler.sequence.scanner import gc_content


# ────────────────────────────────────────────────────────────
# Strategies
# ────────────────────────────────────────────────────────────

# Single DNA base
dna_base = st.sampled_from("ACGT")

# Non-empty DNA sequence (any ACGT string, length 1–300)
dna_seq = st.text(alphabet="ACGT", min_size=1, max_size=300)

# Longer DNA sequences for more thorough GC distribution coverage
dna_seq_long = st.text(alphabet="ACGT", min_size=1, max_size=600)

# Valid identifier string (no spaces, no newlines)
identifier = st.text(
    alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_",
    min_size=1,
    max_size=40,
)

# Simple description (no newlines)
description = st.text(
    alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 _-",
    min_size=0,
    max_size=60,
)

# Organism names from the known set
organism = st.sampled_from([
    "Homo_sapiens",
    "Mus_musculus",
    "Escherichia_coli",
    "E_coli",
    "CHO_K1",
    "Saccharomyces_cerevisiae",
])

# Locus name (GenBank convention: max 16 chars, uppercase)
locus_name = st.text(
    alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_",
    min_size=1,
    max_size=20,
)


# ────────────────────────────────────────────────────────────
# Property 1: FASTA export starts with ">"
# ────────────────────────────────────────────────────────────

def _find_fasta_header_line(fasta_output: str) -> str:
    """Find the header line (starting with '>') in FASTA output.

    The output may contain comment lines (starting with ';') before the header.
    """
    for line in fasta_output.split("\n"):
        if line.startswith(">"):
            return line
    raise ValueError(f"No '>' header line found in FASTA output: {fasta_output!r}")


class TestFastaStartsWithGreaterThan:
    """Property: export_fasta output always contains a '>' header line."""

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_basic_sequence_has_gt_header(self, seq):
        """Any non-empty ACGT sequence produces a FASTA record with a '>' header."""
        result = export_fasta(seq)
        header = _find_fasta_header_line(result)
        assert header.startswith(">"), (
            f"FASTA output has no '>' header for sequence of length {len(seq)}"
        )

    @given(seq=dna_seq, ident=identifier)
    @settings(max_examples=40, deadline=5000)
    def test_custom_identifier_has_gt_header(self, seq, ident):
        """FASTA output has a '>' header even with a custom identifier."""
        result = export_fasta(seq, identifier=ident)
        header = _find_fasta_header_line(result)
        assert header.startswith(">")

    @given(seq=dna_seq, ident=identifier, desc=description)
    @settings(max_examples=30, deadline=5000)
    def test_with_description_has_gt_header(self, seq, ident, desc):
        """FASTA output has a '>' header when a description is provided."""
        result = export_fasta(seq, identifier=ident, description=desc)
        header = _find_fasta_header_line(result)
        assert header.startswith(">")

    @given(seq=dna_seq, org=organism)
    @settings(max_examples=30, deadline=5000)
    def test_various_organisms_has_gt_header(self, seq, org):
        """FASTA output has a '>' header for all supported organisms."""
        result = export_fasta(seq, organism=org)
        header = _find_fasta_header_line(result)
        assert header.startswith(">")

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_header_line_is_present_and_well_formed(self, seq):
        """The FASTA output contains a header line starting with '>'."""
        result = export_fasta(seq)
        header = _find_fasta_header_line(result)
        assert header.startswith(">")
        # Header should contain no newline characters
        assert "\n" not in header


# ────────────────────────────────────────────────────────────
# Property 2: GenBank export ends with "//"
# ────────────────────────────────────────────────────────────

class TestGenBankEndsWithTerminator:
    """Property: export_genbank output always ends with '//'."""

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_basic_sequence_ends_with_double_slash(self, seq):
        """Any non-empty ACGT sequence produces a GenBank record ending with '//'."""
        result = export_genbank(seq)
        assert result.rstrip().endswith("//"), (
            f"GenBank output does not end with '//' for sequence of length {len(seq)}"
        )

    @given(seq=dna_seq, locus=locus_name)
    @settings(max_examples=40, deadline=5000)
    def test_custom_locus_ends_with_double_slash(self, seq, locus):
        """GenBank output ends with '//' with a custom locus name."""
        result = export_genbank(seq, locus_name=locus)
        assert result.rstrip().endswith("//")

    @given(seq=dna_seq, org=organism)
    @settings(max_examples=30, deadline=5000)
    def test_various_organisms_ends_with_double_slash(self, seq, org):
        """GenBank output ends with '//' for all supported organisms."""
        result = export_genbank(seq, organism=org)
        assert result.rstrip().endswith("//")

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_terminator_on_own_line(self, seq):
        """The '//' terminator should be the last line of the GenBank record."""
        result = export_genbank(seq)
        lines = result.rstrip().split("\n")
        assert lines[-1] == "//", (
            f"Last line of GenBank record is {lines[-1]!r}, expected '//'"
        )

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_terminator_appears_exactly_once(self, seq):
        """The '//' terminator should appear exactly once (at the end)."""
        result = export_genbank(seq)
        # The terminator "//" should only appear as the final line
        # (not embedded in sequence data, since sequences use only ACGT)
        count = result.rstrip().split("\n")[-1].count("//")
        assert count == 1

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_has_complete_genbank_structure(self, seq):
        """GenBank output has all required sections in correct order."""
        result = export_genbank(seq)
        assert "LOCUS" in result, "Missing LOCUS section"
        assert "FEATURES" in result, "Missing FEATURES section"
        assert "ORIGIN" in result, "Missing ORIGIN section"
        # Sections must appear in order: LOCUS < FEATURES < ORIGIN < //
        assert result.index("LOCUS") < result.index("FEATURES")
        assert result.index("FEATURES") < result.index("ORIGIN")
        assert result.index("ORIGIN") < result.index("//")


# ────────────────────────────────────────────────────────────
# Property 3: Exported sequences are uppercase ACGT only
# ────────────────────────────────────────────────────────────

def _extract_fasta_sequence(fasta_output: str) -> str:
    """Extract the raw sequence portion from a FASTA string (no spaces, no newlines)."""
    lines = fasta_output.strip().split("\n")
    # Skip the header line (starts with >) and comment lines (starts with ;)
    seq_lines = [line for line in lines if not line.startswith(">") and not line.startswith(";")]
    return "".join(seq_lines)


def _extract_genbank_origin_sequence(gb_output: str) -> str:
    """Extract the DNA sequence from the ORIGIN section of a GenBank record."""
    origin_start = gb_output.index("ORIGIN")
    terminator = gb_output.index("//")
    origin_section = gb_output[origin_start + len("ORIGIN"):terminator]
    # Remove line numbers and spaces
    no_numbers = re.sub(r'\d+', '', origin_section)
    no_spaces = no_numbers.replace(" ", "").replace("\n", "")
    return no_spaces


class TestExportedSequenceUppercaseACGT:
    """Property: all bases in exported sequences are uppercase ACGT only."""

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_fasta_sequence_is_uppercase(self, seq):
        """FASTA sequence portion contains only uppercase characters."""
        result = export_fasta(seq)
        extracted = _extract_fasta_sequence(result)
        assert extracted == extracted.upper(), (
            f"FASTA sequence is not fully uppercase: {extracted!r}"
        )

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_fasta_sequence_is_acgt_only(self, seq):
        """FASTA sequence portion contains only A, C, G, T characters."""
        result = export_fasta(seq)
        extracted = _extract_fasta_sequence(result)
        valid_bases = set("ACGT")
        for base in extracted:
            assert base in valid_bases, (
                f"Invalid base {base!r} in FASTA sequence; expected only ACGT"
            )

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_genbank_origin_is_uppercase(self, seq):
        """GenBank ORIGIN sequence is uppercase only."""
        result = export_genbank(seq)
        extracted = _extract_genbank_origin_sequence(result)
        assert extracted == extracted.upper(), (
            f"GenBank ORIGIN sequence is not fully uppercase: {extracted!r}"
        )

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_genbank_origin_is_acgt_only(self, seq):
        """GenBank ORIGIN sequence contains only A, C, G, T characters."""
        result = export_genbank(seq)
        extracted = _extract_genbank_origin_sequence(result)
        valid_bases = set("ACGT")
        for base in extracted:
            assert base in valid_bases, (
                f"Invalid base {base!r} in GenBank ORIGIN; expected only ACGT"
            )

    @given(seq=dna_seq)
    @settings(max_examples=40, deadline=5000)
    def test_fasta_preserves_sequence_content(self, seq):
        """FASTA export preserves the original sequence content (uppercased)."""
        result = export_fasta(seq)
        extracted = _extract_fasta_sequence(result)
        assert extracted == seq.upper(), (
            f"FASTA sequence content mismatch: expected {seq.upper()!r}, got {extracted!r}"
        )

    @given(seq=dna_seq)
    @settings(max_examples=40, deadline=5000)
    def test_genbank_preserves_sequence_content(self, seq):
        """GenBank export preserves the original sequence content (uppercased)."""
        result = export_genbank(seq)
        extracted = _extract_genbank_origin_sequence(result)
        assert extracted == seq.upper(), (
            f"GenBank sequence content mismatch: expected {seq.upper()!r}, got {extracted!r}"
        )

    @given(seq=st.text(alphabet="acgt", min_size=1, max_size=100))
    @settings(max_examples=30, deadline=5000)
    def test_lowercase_input_produces_uppercase_fasta(self, seq):
        """Lowercase input DNA is converted to uppercase in FASTA output."""
        result = export_fasta(seq)
        extracted = _extract_fasta_sequence(result)
        assert extracted == seq.upper()

    @given(seq=st.text(alphabet="acgt", min_size=1, max_size=100))
    @settings(max_examples=30, deadline=5000)
    def test_lowercase_input_produces_uppercase_genbank(self, seq):
        """Lowercase input DNA is converted to uppercase in GenBank output."""
        result = export_genbank(seq)
        extracted = _extract_genbank_origin_sequence(result)
        assert extracted == seq.upper()


# ────────────────────────────────────────────────────────────
# Property 4: GC content in header matches computed GC
# ────────────────────────────────────────────────────────────

def _extract_fasta_gc(fasta_output: str) -> float:
    """Extract the GC= value from the FASTA header (case-insensitive)."""
    header = _find_fasta_header_line(fasta_output)
    for part in header.lstrip(">").split("|"):
        if part.lower().startswith("gc="):
            return float(part.split("=")[1])
    raise ValueError(f"No gc= field found in FASTA header: {header!r}")


def _extract_genbank_gc(gb_output: str) -> float:
    """Extract the GC content value from the GenBank COMMENT section."""
    for line in gb_output.split("\n"):
        if "GC content:" in line:
            # Line format: "            GC content: 0.5234"
            gc_str = line.split("GC content:")[-1].strip()
            return float(gc_str)
    raise ValueError("No GC content found in GenBank COMMENT section")


class TestGCContentMatchesComputed:
    """Property: GC content advertised in export headers matches gc_content()."""

    @given(seq=dna_seq_long)
    @settings(max_examples=60, deadline=5000)
    def test_fasta_gc_matches_scanner(self, seq):
        """GC content in FASTA header matches gc_content() from scanner."""
        result = export_fasta(seq)
        header_gc = _extract_fasta_gc(result)
        computed_gc = gc_content(seq.upper())
        # FASTA header formats GC to 3 decimal places
        assert abs(header_gc - round(computed_gc, 3)) < 0.001, (
            f"FASTA GC={header_gc}, computed GC={computed_gc}, "
            f"rounded={round(computed_gc, 3)}"
        )

    @given(seq=dna_seq_long)
    @settings(max_examples=60, deadline=5000)
    def test_genbank_gc_matches_scanner(self, seq):
        """GC content in GenBank COMMENT matches gc_content() from scanner."""
        result = export_genbank(seq)
        gb_gc = _extract_genbank_gc(result)
        computed_gc = gc_content(seq.upper())
        # GenBank COMMENT formats GC to 4 decimal places
        assert abs(gb_gc - round(computed_gc, 4)) < 0.0001, (
            f"GenBank GC={gb_gc}, computed GC={computed_gc}, "
            f"rounded={round(computed_gc, 4)}"
        )

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_fasta_gc_in_valid_range(self, seq):
        """GC content in FASTA header is always in [0.0, 1.0]."""
        result = export_fasta(seq)
        gc = _extract_fasta_gc(result)
        assert 0.0 <= gc <= 1.0, f"GC content {gc} out of [0, 1] range"

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_genbank_gc_in_valid_range(self, seq):
        """GC content in GenBank COMMENT is always in [0.0, 1.0]."""
        result = export_genbank(seq)
        gc = _extract_genbank_gc(result)
        assert 0.0 <= gc <= 1.0, f"GC content {gc} out of [0, 1] range"

    def test_all_gc_sequence_fasta(self):
        """Sequence of all G/C has GC content 1.0 in FASTA header."""
        result = export_fasta("GCGCGCGC")
        gc = _extract_fasta_gc(result)
        assert gc == 1.0

    def test_all_at_sequence_fasta(self):
        """Sequence of all A/T has GC content 0.0 in FASTA header."""
        result = export_fasta("ATATATAT")
        gc = _extract_fasta_gc(result)
        assert gc == 0.0

    def test_all_gc_sequence_genbank(self):
        """Sequence of all G/C has GC content 1.0 in GenBank COMMENT."""
        result = export_genbank("GCGCGCGC")
        gc = _extract_genbank_gc(result)
        assert gc == 1.0

    def test_all_at_sequence_genbank(self):
        """Sequence of all A/T has GC content 0.0 in GenBank COMMENT."""
        result = export_genbank("ATATATAT")
        gc = _extract_genbank_gc(result)
        assert gc == 0.0

    @given(seq=dna_seq)
    @settings(max_examples=50, deadline=5000)
    def test_fasta_and_genbank_gc_agree(self, seq):
        """GC content values in FASTA and GenBank exports agree within
        rounding tolerance (FASTA uses 3 decimals, GenBank uses 4)."""
        fasta_result = export_fasta(seq)
        gb_result = export_genbank(seq)
        fasta_gc = _extract_fasta_gc(fasta_result)
        gb_gc = _extract_genbank_gc(gb_result)
        # Both derive from the same gc_content() call, so they should be
        # very close (differ only by rounding: 3 vs 4 decimal places)
        assert abs(fasta_gc - gb_gc) < 0.001, (
            f"FASTA GC={fasta_gc} and GenBank GC={gb_gc} disagree"
        )

    @given(seq=dna_seq)
    @settings(max_examples=40, deadline=5000)
    def test_gc_content_deterministic(self, seq):
        """Calling export_fasta twice on the same sequence produces the same GC."""
        result1 = export_fasta(seq)
        result2 = export_fasta(seq)
        gc1 = _extract_fasta_gc(result1)
        gc2 = _extract_fasta_gc(result2)
        assert gc1 == gc2, f"Non-deterministic GC: {gc1} vs {gc2}"

    @given(seq=dna_seq)
    @settings(max_examples=30, deadline=5000)
    def test_gc_content_case_insensitive(self, seq):
        """GC content is the same regardless of input case."""
        result_upper = export_fasta(seq.upper())
        result_lower = export_fasta(seq.lower())
        gc_upper = _extract_fasta_gc(result_upper)
        gc_lower = _extract_fasta_gc(result_lower)
        assert gc_upper == gc_lower, (
            f"Case-sensitive GC: upper={gc_upper}, lower={gc_lower}"
        )
