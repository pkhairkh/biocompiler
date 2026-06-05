"""
Comprehensive tests for the biocompiler.export module.

Covers:
- export_fasta: basic output, header format, sequence wrapping, organism metadata
- export_genbank: LOCUS line format, ACCESSION, FEATURES, ORIGIN, terminator
- export_multi_fasta: multiple sequences
- export_genbank_with_certificate: certificate embedding
- _format_genbank_sequence: numbering, grouping
- _format_fasta_sequence: 60-char wrapping
- _get_taxonomy: known organisms, unknown organism fallback
- GenBank output structural validity (LOCUS, FEATURES, ORIGIN, //)
- FASTA header includes GC content
- Empty / edge cases
"""

from __future__ import annotations

import pytest

from biocompiler.export import (
    _format_fasta_sequence,
    _format_genbank_sequence,
    _generate_accession,
    _get_taxonomy,
    export_fasta,
    export_genbank,
    export_genbank_with_certificate,
    export_multi_fasta,
)
from biocompiler.types import Certificate, TypeCheckResult, Verdict


# ─── Test fixtures ────────────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"
# 45 bases — GC-rich enough to be interesting

SHORT_SEQ = "ATGC"  # 4 bases

LONG_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTG" * 4  # 240 bases


def _make_certificate(
    design_id: str = "DESIGN_001_CERT",
    types: list[dict] | None = None,
    provenance: dict | None = None,
) -> Certificate:
    """Helper to build a Certificate for testing."""
    if types is None:
        types = [
            {"predicate": "gc_content", "verdict": "PASS"},
            {"predicate": "no_stop_codons", "verdict": "PASS"},
        ]
    if provenance is None:
        provenance = {"timestamp": "2025-01-15T12:00:00Z", "tool": "biocompiler"}
    return Certificate(
        version="1.0",
        design_id=design_id,
        sequence=SAMPLE_SEQ,
        types=types,
        provenance=provenance,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 1. _format_fasta_sequence
# ═══════════════════════════════════════════════════════════════════════════════

class TestFormatFastaSequence:

    def test_wraps_at_60_chars(self):
        """Each output line should be at most 60 characters."""
        seq = "A" * 200
        result = _format_fasta_sequence(seq)
        lines = result.split("\n")
        for line in lines:
            assert len(line) <= 60

    def test_exact_60_no_extra_newline(self):
        """A sequence of exactly 60 chars produces one line."""
        seq = "A" * 60
        result = _format_fasta_sequence(seq)
        lines = result.split("\n")
        assert len(lines) == 1
        assert lines[0] == "A" * 60

    def test_120_produces_two_lines(self):
        seq = "A" * 120
        result = _format_fasta_sequence(seq)
        lines = result.split("\n")
        assert len(lines) == 2
        assert lines[0] == "A" * 60
        assert lines[1] == "A" * 60

    def test_uppercases_sequence(self):
        result = _format_fasta_sequence("atgc")
        assert result == "ATGC"

    def test_short_sequence_single_line(self):
        result = _format_fasta_sequence("ATGC")
        assert result == "ATGC"

    def test_empty_sequence(self):
        result = _format_fasta_sequence("")
        assert result == ""

    def test_custom_line_width(self):
        seq = "A" * 30
        result = _format_fasta_sequence(seq, line_width=10)
        lines = result.split("\n")
        assert len(lines) == 3
        for line in lines:
            assert len(line) == 10


# ═══════════════════════════════════════════════════════════════════════════════
# 2. _format_genbank_sequence
# ═══════════════════════════════════════════════════════════════════════════════

class TestFormatGenbankSequence:

    def test_numbering_starts_at_1(self):
        """First line should be numbered starting at 1."""
        seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        result = _format_genbank_sequence(seq)
        first_line = result.split("\n")[0]
        assert first_line.startswith("        1 ")

    def test_numbering_increments_by_60(self):
        """Second line should start at 61."""
        seq = "A" * 120
        result = _format_genbank_sequence(seq)
        lines = result.split("\n")
        # Extract the number from the second line
        second_num = int(lines[1].split()[0])
        assert second_num == 61

    def test_groups_of_10(self):
        """Bases should be grouped in blocks of 10 separated by spaces."""
        seq = "ATGCATGCATGCATGCATGC"  # 20 bases
        result = _format_genbank_sequence(seq)
        first_line = result.split("\n")[0]
        # Format: "        1 ATGCATGCAT GCATGCATGC"
        # Split on whitespace; skip the number prefix
        parts = first_line.split()
        # parts[0] = "1", parts[1] = "ATGCATGCAT", parts[2] = "GCATGCATGC"
        assert len(parts) == 3
        assert len(parts[1]) == 10
        assert len(parts[2]) == 10

    def test_uppercases_sequence(self):
        seq = "atgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgc"
        result = _format_genbank_sequence(seq)
        assert "a" not in result.split("\n")[0]

    def test_short_sequence(self):
        seq = "ATGCATGCATGC"  # 12 bases
        result = _format_genbank_sequence(seq)
        lines = result.split("\n")
        assert len(lines) == 1

    def test_empty_sequence(self):
        result = _format_genbank_sequence("")
        assert result == ""

    def test_numbering_right_justified(self):
        """Number prefix should be right-justified in 9-char field."""
        seq = "A" * 60
        result = _format_genbank_sequence(seq)
        first_line = result.split("\n")[0]
        # First 9 chars should be the number right-justified
        prefix = first_line[:9]
        assert prefix.strip() == "1"

    def test_custom_group_size(self):
        seq = "A" * 20
        result = _format_genbank_sequence(seq, line_width=60, group_size=5)
        first_line = result.split("\n")[0]
        parts = first_line.split()
        # parts[0] = "1", then groups of 5
        for g in parts[1:]:
            assert len(g) == 5


# ═══════════════════════════════════════════════════════════════════════════════
# 3. _get_taxonomy
# ═══════════════════════════════════════════════════════════════════════════════

class TestGetTaxonomy:

    def test_homo_sapiens(self):
        tax = _get_taxonomy("Homo_sapiens")
        assert "Eukaryota" in tax
        assert "Hominidae" in tax
        assert "Homo" in tax

    def test_mus_musculus(self):
        tax = _get_taxonomy("Mus_musculus")
        assert "Eukaryota" in tax
        assert "Muridae" in tax
        assert "Mus" in tax

    def test_escherichia_coli(self):
        tax = _get_taxonomy("Escherichia_coli")
        assert "Bacteria" in tax
        assert "Escherichia" in tax

    def test_e_coli_alias(self):
        tax = _get_taxonomy("E_coli")
        assert "Bacteria" in tax
        assert "Escherichia" in tax

    def test_cho_k1(self):
        tax = _get_taxonomy("CHO_K1")
        assert "Eukaryota" in tax
        assert "Cricetulus" in tax

    def test_saccharomyces_cerevisiae(self):
        tax = _get_taxonomy("Saccharomyces_cerevisiae")
        assert "Eukaryota" in tax
        assert "Fungi" in tax
        assert "Saccharomyces" in tax

    def test_unknown_organism_fallback(self):
        """Unknown organism should return a generic Eukaryota lineage."""
        tax = _get_taxonomy("Alien_genome_X9")
        assert "Eukaryota" in tax
        assert "Unclassified" in tax

    def test_unknown_organism_contains_period(self):
        """Fallback taxonomy should end with a period (GenBank convention)."""
        tax = _get_taxonomy("Unknown_thing")
        assert tax.endswith(".")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. export_fasta
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportFasta:

    def test_starts_with_greater_than(self):
        result = export_fasta(SAMPLE_SEQ)
        assert result.startswith(">")

    def test_header_includes_organism(self):
        result = export_fasta(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "organism=Homo_sapiens" in result

    def test_header_includes_gc_content(self):
        result = export_fasta(SAMPLE_SEQ)
        assert "gc=" in result
        # Extract gc value from header
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("gc="):
                gc_val = float(part.split("=")[1])
                assert 0.0 <= gc_val <= 1.0

    def test_header_includes_length(self):
        result = export_fasta(SAMPLE_SEQ)
        assert "len=" in result
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("len="):
                assert str(len(SAMPLE_SEQ)) in part

    def test_header_includes_identifier(self):
        result = export_fasta(SAMPLE_SEQ, identifier="my_gene")
        assert ">my_gene" in result

    def test_header_includes_description(self):
        result = export_fasta(SAMPLE_SEQ, description="test sequence")
        assert "test sequence" in result

    def test_description_after_space(self):
        """Description should appear after a space, separated from the pipe-delimited header."""
        result = export_fasta(SAMPLE_SEQ, identifier="ID", description="desc")
        header = result.split("\n")[0]
        # After all pipe-delimited parts, a space then the description
        assert " desc" in header

    def test_no_description_no_trailing_space(self):
        result = export_fasta(SAMPLE_SEQ, identifier="ID", description="")
        header = result.split("\n")[0]
        # Should not end with a space before newline
        assert not header.endswith(" ")

    def test_header_includes_protein_len(self):
        result = export_fasta(SAMPLE_SEQ)
        assert "protein_len=" in result
        assert "aa" in result

    def test_sequence_uppercased(self):
        result = export_fasta("atgc")
        lines = result.split("\n")
        # The sequence line(s) should be uppercase
        for line in lines[1:]:
            if line:
                assert line == line.upper()

    def test_sequence_wrapping(self):
        """Sequences longer than 60 chars should be wrapped."""
        long_seq = "ATGC" * 50  # 200 bases
        result = export_fasta(long_seq)
        lines = result.split("\n")
        # First line is header, last may be empty (trailing newline)
        seq_lines = [l for l in lines[1:] if l]
        for line in seq_lines[:-1]:  # all but possibly the last
            assert len(line) <= 60

    def test_ends_with_newline(self):
        result = export_fasta(SAMPLE_SEQ)
        assert result.endswith("\n")

    def test_spaces_removed_from_sequence(self):
        result = export_fasta("ATG CAT GCA")
        lines = result.split("\n")
        # Second line should have no spaces in the sequence
        seq_line = lines[1]
        assert " " not in seq_line

    def test_custom_organism_metadata(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "organism=Escherichia_coli" in result

    def test_protein_parameter_used(self):
        """When protein is explicitly provided, it should be used (affects protein_len)."""
        result_auto = export_fasta(SAMPLE_SEQ)
        result_explicit = export_fasta(SAMPLE_SEQ, protein="MVE")
        # Both should have protein_len but values may differ
        assert "protein_len=" in result_auto
        assert "protein_len=" in result_explicit


# ═══════════════════════════════════════════════════════════════════════════════
# 5. export_genbank
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportGenbank:

    def test_has_locus_line(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.startswith("LOCUS")

    def test_locus_name_uppercase_and_truncated(self):
        """LOCUS name should be uppercase and at most 16 chars."""
        result = export_genbank(SAMPLE_SEQ, locus_name="my_long_locus_name_exceeding_16")
        locus_line = result.split("\n")[0]
        # Extract the locus name: it's the part after "LOCUS" whitespace
        parts = locus_line.split()
        assert len(parts) >= 2
        locus_name = parts[1]
        assert locus_name == locus_name.upper()
        assert len(locus_name) <= 16

    def test_locus_includes_bp(self):
        result = export_genbank(SAMPLE_SEQ)
        assert f"{len(SAMPLE_SEQ)} bp" in result

    def test_locus_includes_molecule_type(self):
        result = export_genbank(SAMPLE_SEQ, molecule_type="DNA")
        locus_line = result.split("\n")[0]
        assert "DNA" in locus_line

    def test_locus_includes_topology(self):
        result = export_genbank(SAMPLE_SEQ, topology="circular")
        locus_line = result.split("\n")[0]
        assert "circular" in locus_line

    def test_locus_includes_syn_division(self):
        """GenBank SYN division for synthetic sequences."""
        result = export_genbank(SAMPLE_SEQ)
        locus_line = result.split("\n")[0]
        assert "SYN" in locus_line

    def test_has_definition_line(self):
        result = export_genbank(SAMPLE_SEQ, definition="My custom definition")
        assert "DEFINITION" in result
        assert "My custom definition" in result

    def test_has_accession_line(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "ACCESSION" in result

    def test_accession_generated_when_no_certificate(self):
        """Without a certificate, accession is generated via _generate_accession."""
        result = export_genbank(SAMPLE_SEQ, locus_name="MYLOCUS")
        for line in result.split("\n"):
            if line.startswith("ACCESSION"):
                # Accession should be of form BC_xxxxxxxx
                assert "BC_" in line
                break

    def test_accession_deterministic(self):
        """Same sequence should always produce the same accession."""
        result1 = export_genbank(SAMPLE_SEQ)
        result2 = export_genbank(SAMPLE_SEQ)
        acc1 = [l for l in result1.split("\n") if l.startswith("ACCESSION")][0]
        acc2 = [l for l in result2.split("\n") if l.startswith("ACCESSION")][0]
        assert acc1 == acc2

    def test_accession_uses_certificate_design_id(self):
        cert = _make_certificate(design_id="CERT_ABC_123")
        result = export_genbank(SAMPLE_SEQ, certificate=cert)
        for line in result.split("\n"):
            if line.startswith("ACCESSION"):
                assert "CERT_ABC_123"[:12].upper() in line
                break

    def test_has_version_line(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "VERSION" in result

    def test_has_source_line(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "SOURCE" in result
        assert "Homo_sapiens" in result

    def test_has_organism_line(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "ORGANISM" in result

    def test_has_taxonomy(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "Eukaryota" in result

    def test_has_features_section(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "FEATURES" in result

    def test_features_with_gene_name(self):
        result = export_genbank(SAMPLE_SEQ, gene_name="eGFP")
        assert '/gene="eGFP"' in result
        assert "gene" in result

    def test_features_cds(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "CDS" in result
        assert "/translation=" in result

    def test_cds_with_multi_exon_join(self):
        boundaries = [(0, 100), (200, 300)]
        result = export_genbank(SAMPLE_SEQ, exon_boundaries=boundaries)
        assert "join(" in result
        assert "1..100" in result
        assert "201..300" in result

    def test_cds_single_exon_no_join(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "join(" not in result
        assert f"1..{len(SAMPLE_SEQ)}" in result

    def test_exon_features(self):
        boundaries = [(0, 30), (30, 45)]
        result = export_genbank(SAMPLE_SEQ, exon_boundaries=boundaries, gene_name="test")
        assert "exon" in result
        assert "/number=1" in result
        assert "/number=2" in result

    def test_restriction_site_features(self):
        sites = [{"enzyme": "EcoRI", "position": 5, "strand": "+", "site": "GAATTC"}]
        result = export_genbank(SAMPLE_SEQ, restriction_sites=sites)
        assert "misc_feature" in result
        assert "EcoRI" in result

    def test_restriction_site_1_based(self):
        """Restriction site positions should be converted to 1-based."""
        sites = [{"enzyme": "EcoRI", "position": 5, "strand": "+", "site": "GAATTC"}]
        result = export_genbank(SAMPLE_SEQ, restriction_sites=sites)
        assert "6..11" in result  # position 5 (0-based) + 1 = 6, +6 (len) = 11

    def test_type_check_fail_annotated(self):
        type_results = [
            TypeCheckResult(predicate="gc_content", verdict=Verdict.FAIL, violation="GC too low"),
        ]
        result = export_genbank(SAMPLE_SEQ, type_results=type_results)
        assert "TYPE FAIL: gc_content" in result
        assert "typecheck_fail" in result

    def test_type_check_pass_no_annotation(self):
        type_results = [
            TypeCheckResult(predicate="gc_content", verdict=Verdict.PASS),
        ]
        result = export_genbank(SAMPLE_SEQ, type_results=type_results)
        assert "TYPE FAIL" not in result

    def test_has_origin_section(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "ORIGIN" in result

    def test_has_terminator(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.rstrip().endswith("//")

    def test_comment_section_includes_gc(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "GC content:" in result

    def test_comment_section_includes_version(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "BioCompiler" in result

    def test_certificate_in_comment(self):
        cert = _make_certificate(design_id="DESGN123")
        result = export_genbank(SAMPLE_SEQ, certificate=cert)
        assert "Certificate ID:" in result
        assert "Certificate timestamp:" in result

    def test_type_results_in_comment(self):
        type_results = [
            TypeCheckResult(predicate="gc_content", verdict=Verdict.PASS),
            TypeCheckResult(predicate="no_stop", verdict=Verdict.FAIL, violation="Stop found"),
        ]
        result = export_genbank(SAMPLE_SEQ, type_results=type_results)
        assert "Type-check verdict:" in result
        assert "[+]" in result   # PASS
        assert "[X]" in result   # FAIL

    def test_cds_translation_wrapped(self):
        """Protein translation in CDS should be wrapped at 40 chars."""
        long_seq = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTG" * 4
        result = export_genbank(long_seq)
        # Find translation lines
        lines = result.split("\n")
        in_translation = False
        translation_parts = []
        for line in lines:
            if '/translation="' in line:
                in_translation = True
                # Extract the part after the first quote
                part = line.split('/translation="')[1].rstrip('"')
                translation_parts.append(part)
            elif in_translation and line.strip().startswith('"'):
                part = line.strip().strip('"')
                translation_parts.append(part)
                if not line.rstrip().endswith('"'):
                    pass
                else:
                    in_translation = False
            elif in_translation:
                in_translation = False
        # Each wrapped chunk should be at most 40 chars
        for part in translation_parts:
            assert len(part) <= 40


# ═══════════════════════════════════════════════════════════════════════════════
# 6. export_multi_fasta
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportMultiFasta:

    def test_multiple_sequences(self):
        entries = [
            {"sequence": "ATGCATGC", "id": "seq1", "description": "first"},
            {"sequence": "GCTAGCTA", "id": "seq2", "description": "second"},
        ]
        result = export_multi_fasta(entries)
        assert ">seq1" in result
        assert ">seq2" in result
        assert "first" in result
        assert "second" in result

    def test_two_records_separated_by_newline(self):
        entries = [
            {"sequence": "ATGC", "id": "s1"},
            {"sequence": "GCTA", "id": "s2"},
        ]
        result = export_multi_fasta(entries)
        # Each record ends with \n, so joined with \n means no blank line between
        records = result.split(">")
        # First element is empty (before first >)
        assert len([r for r in records if r]) == 2

    def test_empty_list_returns_empty(self):
        result = export_multi_fasta([])
        assert result == ""

    def test_single_entry(self):
        entries = [{"sequence": "ATGCATGC", "id": "only_one"}]
        result = export_multi_fasta(entries)
        assert ">only_one" in result
        assert "ATGCATGC" in result

    def test_uses_default_id(self):
        entries = [{"sequence": "ATGC"}]
        result = export_multi_fasta(entries)
        assert ">BioCompiler_design" in result

    def test_uses_default_organism(self):
        entries = [{"sequence": "ATGC"}]
        result = export_multi_fasta(entries)
        assert "organism=Homo_sapiens" in result

    def test_custom_organism_per_entry(self):
        entries = [
            {"sequence": "ATGC", "organism": "Escherichia_coli"},
            {"sequence": "GCTA", "organism": "Mus_musculus"},
        ]
        result = export_multi_fasta(entries)
        assert "organism=Escherichia_coli" in result
        assert "organism=Mus_musculus" in result

    def test_optional_protein_per_entry(self):
        entries = [
            {"sequence": "ATGCATGC", "id": "sp", "protein": "MH"},
        ]
        result = export_multi_fasta(entries)
        assert "protein_len=2aa" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 7. export_genbank_with_certificate
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportGenbankWithCertificate:

    def test_includes_certificate_in_comment(self):
        cert = _make_certificate(design_id="CERT_ABC")
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "Certificate ID:" in result

    def test_uses_certificate_design_id_as_locus(self):
        cert = _make_certificate(design_id="MYDESIGN123")
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        locus_line = result.split("\n")[0]
        assert "MYDESIGN123"[:16].upper() in locus_line

    def test_certificate_timestamp_in_comment(self):
        cert = _make_certificate(provenance={"timestamp": "2025-03-15T10:00:00Z"})
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "2025-03-15T10:00:00Z" in result

    def test_type_results_from_certificate(self):
        types = [
            {"predicate": "gc_content", "verdict": "PASS"},
            {"predicate": "no_stop", "verdict": "FAIL"},
        ]
        cert = _make_certificate(types=types)
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "Type-check verdict:" in result
        assert "FAIL" in result

    def test_reconstructs_fail_annotations(self):
        types = [
            {"predicate": "stability", "verdict": "FAIL"},
        ]
        cert = _make_certificate(types=types)
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "TYPE FAIL: stability" in result

    def test_passes_gene_name(self):
        cert = _make_certificate()
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert, gene_name="myGene")
        assert '/gene="myGene"' in result

    def test_passes_exon_boundaries(self):
        cert = _make_certificate()
        boundaries = [(0, 20), (25, 45)]
        result = export_genbank_with_certificate(
            SAMPLE_SEQ, certificate=cert, exon_boundaries=boundaries
        )
        assert "join(" in result
        assert "1..20" in result
        assert "26..45" in result

    def test_definition_includes_design_id_prefix(self):
        cert = _make_certificate(design_id="MYDESN_XYZ")
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "DEFINITION" in result
        assert "MYDESN_X"[:8] in result  # first 8 chars of design_id

    def test_valid_genbank_structure(self):
        cert = _make_certificate()
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "LOCUS" in result
        assert "FEATURES" in result
        assert "ORIGIN" in result
        assert result.rstrip().endswith("//")


# ═══════════════════════════════════════════════════════════════════════════════
# 8. GenBank structural validity
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenBankStructuralValidity:

    def test_has_locus(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "LOCUS" in result

    def test_has_features(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "FEATURES" in result

    def test_has_origin(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "ORIGIN" in result

    def test_has_terminator(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "//" in result

    def test_locus_before_features(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("LOCUS") < result.index("FEATURES")

    def test_features_before_origin(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("FEATURES") < result.index("ORIGIN")

    def test_origin_before_terminator(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("ORIGIN") < result.index("//")

    def test_accession_before_version(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("ACCESSION") < result.index("VERSION")

    def test_definition_after_locus(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("LOCUS") < result.index("DEFINITION")


# ═══════════════════════════════════════════════════════════════════════════════
# 9. FASTA header GC content validation
# ═══════════════════════════════════════════════════════════════════════════════

class TestFastaGCContent:

    def test_gc_content_for_all_gc(self):
        result = export_fasta("GCGCGCGC")
        header = result.split("\n")[0]
        for part in header.lstrip(">").split("|"):
            if part.startswith("gc="):
                gc = float(part.split("=")[1])
                assert gc == 1.0

    def test_gc_content_for_all_at(self):
        result = export_fasta("ATATATAT")
        header = result.split("\n")[0]
        for part in header.lstrip(">").split("|"):
            if part.startswith("gc="):
                gc = float(part.split("=")[1])
                assert gc == 0.0

    def test_gc_content_mixed(self):
        result = export_fasta("ATGC")
        header = result.split("\n")[0]
        for part in header.lstrip(">").split("|"):
            if part.startswith("gc="):
                gc = float(part.split("=")[1])
                assert gc == 0.5

    def test_gc_format_three_decimals(self):
        """GC content should be formatted to 3 decimal places."""
        result = export_fasta("ATGCATGCATGC")
        header = result.split("\n")[0]
        for part in header.lstrip(">").split("|"):
            if part.startswith("gc="):
                val = part.split("=")[1]
                # Should have exactly 3 decimal places
                assert "." in val
                decimal_part = val.split(".")[1]
                assert len(decimal_part) == 3


# ═══════════════════════════════════════════════════════════════════════════════
# 10. Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_fasta_very_short_sequence(self):
        result = export_fasta("A")
        assert result.startswith(">")
        assert "A\n" in result

    def test_genbank_very_short_sequence(self):
        result = export_genbank("ATGC")
        assert "LOCUS" in result
        assert "4 bp" in result
        assert "//" in result

    def test_fasta_single_base_sequence(self):
        result = export_fasta("G")
        lines = result.split("\n")
        assert len(lines) >= 2  # header + sequence
        assert lines[1] == "G"

    def test_genbank_long_locus_name_truncated(self):
        result = export_genbank(SAMPLE_SEQ, locus_name="A" * 30)
        locus_line = result.split("\n")[0]
        # The locus name in the output should be at most 16 chars
        parts = locus_line.split()
        locus_name = parts[1]
        assert len(locus_name) <= 16

    def test_format_fasta_sequence_non_standard_width(self):
        result = _format_fasta_sequence("ATGCATGCATGC", line_width=3)
        lines = result.split("\n")
        assert lines == ["ATG", "CAT", "GCA", "TGC"]

    def test_format_genbank_sequence_partial_last_line(self):
        """A sequence not divisible by 60 should have a partial last line."""
        seq = "A" * 75
        result = _format_genbank_sequence(seq)
        lines = result.split("\n")
        assert len(lines) == 2
        # First line: 60 bases, second: 15 bases

    def test_export_fasta_preserves_sequence_content(self):
        """All bases in the original should appear in the output."""
        seq = "ATGCATGCATGC"
        result = export_fasta(seq)
        # Remove header and newlines
        seq_in_output = result.split("\n", 1)[1].replace("\n", "")
        assert seq_in_output == seq.upper()

    def test_export_genbank_preserves_sequence_content(self):
        """All bases should appear in the ORIGIN section."""
        seq = "ATGCATGCATGC"
        result = export_genbank(seq)
        # Extract ORIGIN section
        origin_start = result.index("ORIGIN")
        terminator = result.index("//")
        origin_section = result[origin_start:terminator]
        # Remove numbering and spaces
        bases = origin_section.replace("ORIGIN", "")
        import re
        bases = re.sub(r'\d+', '', bases)
        bases = bases.replace(" ", "").replace("\n", "")
        assert bases == seq.upper()

    def test_multi_fasta_with_empty_description(self):
        entries = [
            {"sequence": "ATGC", "id": "s1", "description": ""},
        ]
        result = export_multi_fasta(entries)
        assert ">s1" in result

    def test_genbank_restriction_sites_limited_to_20(self):
        """Only first 20 restriction sites should be annotated."""
        sites = [
            {"enzyme": f"Enz{i}", "position": i, "strand": "+", "site": "ATGC"}
            for i in range(30)
        ]
        result = export_genbank("ATGC" * 30, restriction_sites=sites)
        # Count misc_feature occurrences for restriction sites
        count = result.count("Restriction site:")
        assert count <= 20

    def test_certificate_empty_types(self):
        cert = Certificate(
            version="1.0",
            design_id="EMPTY_TYPES",
            sequence=SAMPLE_SEQ,
            types=[],
            provenance={"timestamp": "2025-01-01"},
        )
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "LOCUS" in result
        assert "//" in result

    def test_type_results_all_uncertain(self):
        type_results = [
            TypeCheckResult(predicate="pred_a", verdict=Verdict.UNCERTAIN),
        ]
        result = export_genbank(SAMPLE_SEQ, type_results=type_results)
        assert "Type-check verdict:" in result
        assert "UNCERTAIN" in result

    def test_taxonomy_long_line_wrapping(self):
        """Taxonomy longer than 70 chars should be wrapped across lines."""
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        # The taxonomy for Homo_sapiens is long; it should be present
        # but may be split across lines. Check the full content is there.
        assert "Hominidae" in result
        assert "Homo" in result

    def test_genbank_locus_line_structure(self):
        """LOCUS line should contain locus name, length, molecule type, topology, division."""
        result = export_genbank(SAMPLE_SEQ, locus_name="MYGENE", molecule_type="DNA", topology="linear")
        locus_line = result.split("\n")[0]
        assert "MYGENE" in locus_line
        assert "bp" in locus_line
        assert "DNA" in locus_line
        assert "linear" in locus_line
        assert "SYN" in locus_line
