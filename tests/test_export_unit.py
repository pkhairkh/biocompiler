"""
Additional unit tests for biocompiler.export — covering uncovered branches.

Focuses on:
- export_batch_fasta
- export_json
- export_full_construct
- _format_full_construct_features
- _serialize_for_json
- _format_sequence_numbered
- _wrap_text
- Edge cases for export_genbank (mRNA molecule type, circular topology, etc.)
"""

from __future__ import annotations

import json
import pytest

from biocompiler.export.core import (
    export_fasta,
    export_genbank,
    export_batch_fasta,
    export_full_construct,
    export_json,
    _format_sequence_numbered,
    _wrap_text,
    _serialize_for_json,
    _generate_accession,
    FastaSequenceEntry,
    RestrictionSiteInfo,
)
from biocompiler.shared.types import Certificate, TypeCheckResult, Verdict


# ── Helpers ──────────────────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"


def _make_certificate(**kwargs):
    defaults = dict(
        version="1.0",
        design_id="DESIGN_001_CERT",
        sequence=SAMPLE_SEQ,
        types=[{"predicate": "gc_content", "verdict": "PASS"}],
        provenance={"timestamp": "2025-01-15T12:00:00Z"},
    )
    defaults.update(kwargs)
    return Certificate(**defaults)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. _format_sequence_numbered
# ═══════════════════════════════════════════════════════════════════════════════

class TestFormatSequenceNumbered:

    def test_basic_format(self):
        result = _format_sequence_numbered("ATGCATGCATGCATGCATGC")
        assert "1" in result
        assert "ATGCATGCAT" in result

    def test_uppercases(self):
        result = _format_sequence_numbered("atgcatgc")
        assert "ATGC" in result

    def test_custom_group_size(self):
        result = _format_sequence_numbered("ATGCATGC", line_width=8, group_size=4)
        parts = result.split()
        assert len(parts) >= 2  # number + groups

    def test_long_sequence_multiple_lines(self):
        result = _format_sequence_numbered("A" * 120)
        lines = result.split("\n")
        assert len(lines) >= 2

    def test_empty_sequence(self):
        result = _format_sequence_numbered("")
        assert result == ""


# ═══════════════════════════════════════════════════════════════════════════════
# 2. _wrap_text
# ═══════════════════════════════════════════════════════════════════════════════

class TestWrapText:

    def test_basic_wrap(self):
        result = _wrap_text("Hello World", width=5)
        assert isinstance(result, str)

    def test_with_indent(self):
        result = _wrap_text("Hello World", width=20, indent=4)
        lines = result.split("\n")
        for line in lines:
            assert line.startswith("    ")  # 4-space indent

    def test_short_text_no_wrap(self):
        result = _wrap_text("Hi", width=80)
        assert result.strip() == "Hi"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. _generate_accession
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenerateAccession:

    def test_format(self):
        acc = _generate_accession()
        assert acc.startswith("BC_")
        assert len(acc) == 11  # BC_ + 8 hex chars

    def test_uniqueness(self):
        acc1 = _generate_accession()
        acc2 = _generate_accession()
        assert acc1 != acc2


# ═══════════════════════════════════════════════════════════════════════════════
# 4. export_batch_fasta
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportBatchFasta:

    def test_basic_batch(self):
        results = [
            {"sequence": "ATGCATGC", "identifier": "gene1"},
            {"sequence": "GCTAGCTA", "identifier": "gene2"},
        ]
        fasta = export_batch_fasta(results, organism="Escherichia_coli")
        assert ">gene1" in fasta
        assert ">gene2" in fasta
        assert "organism=Escherichia_coli" in fasta

    def test_empty_list(self):
        fasta = export_batch_fasta([])
        assert fasta == "\n" or fasta == ""

    def test_skip_empty_sequences(self):
        results = [
            {"sequence": "ATGCATGC", "identifier": "gene1"},
            {"sequence": "", "identifier": "gene2"},
        ]
        fasta = export_batch_fasta(results)
        assert ">gene1" in fasta
        assert ">gene2" not in fasta

    def test_custom_cai_per_entry(self):
        results = [
            {"sequence": "ATGCATGC", "identifier": "gene1", "cai": 0.95},
        ]
        fasta = export_batch_fasta(results)
        assert "cai=0.9500" in fasta

    def test_default_organism(self):
        results = [{"sequence": "ATGC"}]
        fasta = export_batch_fasta(results)
        assert "organism=Homo_sapiens" in fasta

    def test_per_entry_organism(self):
        results = [
            {"sequence": "ATGC", "organism": "Escherichia_coli"},
        ]
        fasta = export_batch_fasta(results, organism="Homo_sapiens")
        assert "organism=Escherichia_coli" in fasta


# ═══════════════════════════════════════════════════════════════════════════════
# 5. export_full_construct
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportFullConstruct:

    def test_basic_construct(self):
        gb = export_full_construct(
            utr5="AAAAAA",
            cds="ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC",
            utr3="TTTTTT",
            organism="Homo_sapiens",
            gene_name="eGFP",
        )
        assert "LOCUS" in gb
        assert "FEATURES" in gb
        assert "ORIGIN" in gb
        assert "//" in gb

    def test_empty_utrs(self):
        gb = export_full_construct(
            utr5="",
            cds=SAMPLE_SEQ,
            utr3="",
            organism="Homo_sapiens",
        )
        assert "LOCUS" in gb
        assert "//" in gb

    def test_empty_cds_raises(self):
        with pytest.raises(ValueError, match="CDS is required"):
            export_full_construct(
                utr5="",
                cds="",
                utr3="TTTTTT",
            )

    def test_all_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            export_full_construct(
                utr5="",
                cds="",
                utr3="",
            )

    def test_includes_5utr_feature(self):
        gb = export_full_construct(
            utr5="AAAAAA",
            cds=SAMPLE_SEQ,
            utr3="TTTTTT",
            gene_name="test",
        )
        assert "5'UTR" in gb

    def test_includes_3utr_feature(self):
        gb = export_full_construct(
            utr5="AAAAAA",
            cds=SAMPLE_SEQ,
            utr3="TTTTTT",
            gene_name="test",
        )
        assert "3'UTR" in gb

    def test_includes_mrna_feature(self):
        gb = export_full_construct(
            utr5="AAAAAA",
            cds=SAMPLE_SEQ,
            utr3="TTTTTT",
        )
        assert "mRNA" in gb

    def test_cai_qualifier(self):
        gb = export_full_construct(
            utr5="",
            cds=SAMPLE_SEQ,
            utr3="",
            cai=0.95,
        )
        assert 'cai="0.9500"' in gb

    def test_protein_translation(self):
        gb = export_full_construct(
            utr5="",
            cds=SAMPLE_SEQ,
            utr3="",
            gene_name="test",
        )
        assert '/translation="' in gb

    def test_sequence_length_correct(self):
        utr5 = "AAAAAA"
        cds = SAMPLE_SEQ
        utr3 = "TTTTTT"
        gb = export_full_construct(utr5=utr5, cds=cds, utr3=utr3)
        total_len = len(utr5) + len(cds) + len(utr3)
        assert f"{total_len} bp" in gb


# ═══════════════════════════════════════════════════════════════════════════════
# 6. export_json
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportJson:

    def test_rejects_non_optimization_result(self):
        with pytest.raises(TypeError, match="Expected OptimizationResult"):
            export_json("not a result")

    def test_rejects_dict(self):
        with pytest.raises(TypeError, match="Expected OptimizationResult"):
            export_json({"sequence": "ATG"})


# ═══════════════════════════════════════════════════════════════════════════════
# 7. _serialize_for_json
# ═══════════════════════════════════════════════════════════════════════════════

class TestSerializeForJson:

    def test_none(self):
        assert _serialize_for_json(None) is None

    def test_string(self):
        assert _serialize_for_json("hello") == "hello"

    def test_int(self):
        assert _serialize_for_json(42) == 42

    def test_float(self):
        assert _serialize_for_json(3.14) == 3.14

    def test_bool(self):
        assert _serialize_for_json(True) is True

    def test_list(self):
        result = _serialize_for_json([1, 2, 3])
        assert result == [1, 2, 3]

    def test_dict(self):
        result = _serialize_for_json({"a": 1})
        assert result == {"a": 1}

    def test_set_sorted(self):
        result = _serialize_for_json({3, 1, 2})
        assert result == [1, 2, 3]

    def test_verdict(self):
        result = _serialize_for_json(Verdict.PASS)
        assert result == "PASS"

    def test_dataclass(self):
        cert = _make_certificate()
        result = _serialize_for_json(cert)
        assert isinstance(result, dict)
        assert "version" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 8. export_genbank edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportGenbankEdgeCases:

    def test_mrna_molecule_type(self):
        gb = export_genbank(SAMPLE_SEQ, molecule_type="mRNA")
        assert "mRNA" in gb

    def test_rna_molecule_type(self):
        gb = export_genbank(SAMPLE_SEQ, molecule_type="RNA")
        assert "RNA" in gb

    def test_circular_topology(self):
        gb = export_genbank(SAMPLE_SEQ, topology="circular")
        assert "circular" in gb

    def test_invalid_topology_defaults_to_linear(self):
        gb = export_genbank(SAMPLE_SEQ, topology="weird")
        assert "linear" in gb

    def test_with_cai_qualifier(self):
        gb = export_genbank(SAMPLE_SEQ, cai=0.95)
        assert 'cai="0.9500"' in gb

    def test_organism_with_underscores(self):
        """Organism name should have underscores replaced with spaces in display."""
        gb = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "Escherichia coli" in gb

    def test_no_comments_mode(self):
        """FASTA without comments."""
        fasta = export_fasta(SAMPLE_SEQ, include_comments=False)
        assert not fasta.startswith(";")

    def test_fasta_with_cai(self):
        fasta = export_fasta(SAMPLE_SEQ, cai=0.99)
        assert "cai=0.9900" in fasta

    def test_restriction_site_strand_minus(self):
        sites = [{"enzyme": "EcoRI", "position": 5, "strand": "-", "site": "GAATTC"}]
        gb = export_genbank(SAMPLE_SEQ, restriction_sites=sites)
        assert "EcoRI" in gb
        assert "-" in gb  # strand annotation

    def test_organism_names_in_features(self):
        gb = export_genbank(SAMPLE_SEQ, organism="Mus_musculus", gene_name="Hbb")
        assert 'organism="Mus musculus"' in gb

    def test_regulatory_feature_for_coding_sequence(self):
        """Codon-optimized sequences should have a regulatory feature."""
        gb = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "regulatory" in gb
        assert "codon_optimization" in gb

    def test_cds_with_exon_join(self):
        boundaries = [(0, 30), (30, 45)]
        gb = export_genbank(SAMPLE_SEQ, exon_boundaries=boundaries)
        assert "join(" in gb
