"""Tests for BioPython round-trip compatibility — v12.1.0 additions.

Covered areas:
1. optimize_seqrecord() — end-to-end SeqRecord in, SeqRecord out
2. from_seqio() — import from Bio.SeqIO parsed objects
3. to_genbank_string() — export SeqRecord to GenBank via Bio.SeqIO
4. to_fasta_string() — export SeqRecord to FASTA via Bio.SeqIO
5. to_seqrecord(result=...) — accept OptimizationResult
6. Feature/annotation preservation during optimization
7. Graceful degradation when BioPython is not installed

All BioPython-dependent tests use ``pytest.importorskip("Bio")`` so the
entire suite is skipped cleanly when BioPython is absent.
"""

from __future__ import annotations

import os
import sys
import tempfile
from unittest import mock

import pytest

# ---------------------------------------------------------------------------
# 0. BioPython availability
# ---------------------------------------------------------------------------

Bio = pytest.importorskip("Bio", reason="BioPython not installed — skipping round-trip tests")

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation


# ---------------------------------------------------------------------------
# 1. optimize_seqrecord — end-to-end BioPython workflow
# ---------------------------------------------------------------------------


class TestOptimizeSeqrecord:
    """Tests for ``optimize_seqrecord`` — the primary round-trip entry point."""

    @pytest.fixture(autouse=True)
    def _check_deps(self):
        """Skip if biocompiler optimization pipeline is unavailable."""
        pytest.importorskip("biocompiler.optimization")
        pytest.importorskip("biocompiler.type_system")

    def _make_record_with_cds(self, dna: str = "ATGGTTTCTAAAGGTGAA",
                              gene_name: str = "eGFP",
                              organism: str = "Homo_sapiens") -> SeqRecord:
        """Helper: build a SeqRecord with a CDS feature and organism annotation."""
        from biocompiler.biopython_compat import to_seqrecord
        return to_seqrecord(dna, organism=organism, gene_name=gene_name)

    def test_returns_seqrecord(self, sample_protein):
        """optimize_seqrecord returns a BioPython SeqRecord."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds()
        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        assert isinstance(result, SeqRecord)

    def test_organism_in_annotations(self):
        """Organism is set in the output annotations."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds()
        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        assert result.annotations["organism"] == "Escherichia_coli"

    def test_organism_from_record_annotation(self):
        """If organism is not specified, it is read from the record annotation."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds(organism="Mus_musculus")
        result = optimize_seqrecord(rec)  # no explicit organism
        assert result.annotations["organism"] == "Mus_musculus"

    def test_organism_space_to_underscore_conversion(self):
        """Organism names with spaces are converted to underscore form."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds(organism="Homo sapiens")
        result = optimize_seqrecord(rec)
        assert result.annotations["organism"] == "Homo_sapiens"

    def test_has_cds_feature(self):
        """Resulting SeqRecord has a CDS feature."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds()
        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        cds = [f for f in result.features if f.type == "CDS"]
        assert len(cds) >= 1

    def test_has_gene_feature(self):
        """Gene name from the input record is preserved."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds(gene_name="testGene")
        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        gene_feats = [f for f in result.features if f.type == "gene"]
        assert len(gene_feats) >= 1
        assert gene_feats[0].qualifiers["gene"] == ["testGene"]

    def test_type_error_on_non_seqrecord(self):
        """optimize_seqrecord raises TypeError for non-SeqRecord input."""
        from biocompiler.biopython_compat import optimize_seqrecord

        with pytest.raises(TypeError, match="SeqRecord"):
            optimize_seqrecord("not_a_record", organism="Escherichia_coli")

    def test_value_error_on_empty_sequence(self):
        """optimize_seqrecord raises error if the record can't be translated."""
        from biocompiler.biopython_compat import optimize_seqrecord
        from biocompiler.exceptions import InvalidProteinError

        rec = SeqRecord(Seq("NNNNNN"), id="bad")
        with pytest.raises((ValueError, InvalidProteinError)):
            optimize_seqrecord(rec, organism="Escherichia_coli")

    def test_original_id_preserved_in_annotations(self):
        """Original record ID is stored for traceability."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds()
        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        assert "biocompiler_original_id" in result.annotations
        assert result.annotations["biocompiler_original_id"] == rec.id

    def test_molecule_type_dna(self):
        """Output record has molecule_type DNA."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = self._make_record_with_cds()
        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        assert result.annotations["molecule_type"] == "DNA"


# ---------------------------------------------------------------------------
# 2. Feature and annotation preservation
# ---------------------------------------------------------------------------


class TestFeaturePreservation:
    """Tests for preserving non-BioCompiler features during optimization."""

    @pytest.fixture(autouse=True)
    def _check_deps(self):
        pytest.importorskip("biocompiler.optimization")

    def test_promoter_feature_preserved(self):
        """Non-BioCompiler features (e.g. promoter) survive optimization."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"

        # Add CDS
        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"translation": ["MVSKGE"]},
        )
        rec.features.append(cds)

        # Add promoter feature (should be preserved)
        promoter = SeqFeature(
            FeatureLocation(0, 5), type="promoter",
            qualifiers={"note": ["T7 promoter"]},
        )
        rec.features.append(promoter)

        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        promoter_feats = [f for f in result.features if f.type == "promoter"]
        assert len(promoter_feats) == 1
        assert "T7 promoter" in promoter_feats[0].qualifiers.get("note", [])

    def test_terminator_feature_preserved(self):
        """Terminator features are preserved."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"

        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"translation": ["MVSKGE"]},
        )
        rec.features.append(cds)

        terminator = SeqFeature(
            FeatureLocation(13, 18), type="terminator",
            qualifiers={"note": ["T7 terminator"]},
        )
        rec.features.append(terminator)

        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        term_feats = [f for f in result.features if f.type == "terminator"]
        assert len(term_feats) == 1

    def test_custom_annotation_preserved(self):
        """Custom annotations (not BioCompiler-managed) are carried over."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["custom_field"] = "preserved_value"

        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"translation": ["MVSKGE"]},
        )
        rec.features.append(cds)

        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        assert result.annotations.get("custom_field") == "preserved_value"

    def test_biocompiler_features_regenerated(self):
        """BioCompiler features (gene, CDS, exon) are regenerated, not duplicated."""
        from biocompiler.biopython_compat import optimize_seqrecord

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"

        # Original gene + CDS features
        gene = SeqFeature(
            FeatureLocation(0, 18), type="gene",
            qualifiers={"gene": ["oldGene"]},
        )
        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"gene": ["oldGene"], "translation": ["MVSKGE"]},
        )
        rec.features.extend([gene, cds])

        result = optimize_seqrecord(rec, organism="Escherichia_coli")
        # Should have exactly one gene and one CDS (regenerated)
        gene_feats = [f for f in result.features if f.type == "gene"]
        cds_feats = [f for f in result.features if f.type == "CDS"]
        assert len(gene_feats) == 1
        assert len(cds_feats) == 1


# ---------------------------------------------------------------------------
# 3. from_seqio — import from Bio.SeqIO parsed objects
# ---------------------------------------------------------------------------


class TestFromSeqio:
    """Tests for ``from_seqio`` — importing from Bio.SeqIO sources."""

    def test_from_fasta_file(self, tmp_path):
        """from_seqio can parse a FASTA file and return list of dicts."""
        from biocompiler.biopython_compat import from_seqio

        # Create a temporary FASTA file
        fasta_content = ">gene1\nATGGTTTCTAAAGGTGAA\n>gene2\nATGGCCAAAGGGTTT\n"
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(fasta_content)

        results = from_seqio(str(fasta_path), format="fasta")
        assert isinstance(results, list)
        assert len(results) == 2
        assert results[0]["sequence"] == "ATGGTTTCTAAAGGTGAA"
        assert results[1]["sequence"] == "ATGGCCAAAGGGTTT"

    def test_from_genbank_file(self, tmp_path):
        """from_seqio can parse a GenBank file."""
        from biocompiler.biopython_compat import from_seqio, to_seqrecord, to_genbank_string

        rec = to_seqrecord("ATGGTTTCTAAAGGTGAA", organism="Escherichia_coli", gene_name="test")
        gb_str = to_genbank_string(rec)

        gb_path = tmp_path / "test.gbk"
        gb_path.write_text(gb_str)

        results = from_seqio(str(gb_path), format="genbank")
        assert len(results) >= 1
        assert results[0]["sequence"] == "ATGGTTTCTAAAGGTGAA"

    def test_file_not_found_raises_error(self):
        """from_seqio raises FileNotFoundError for missing files."""
        from biocompiler.biopython_compat import from_seqio

        with pytest.raises(FileNotFoundError):
            from_seqio("/nonexistent/path/file.fasta", format="fasta")

    def test_returns_list_of_dicts(self, tmp_path):
        """Each result from from_seqio is a dict with expected keys."""
        from biocompiler.biopython_compat import from_seqio

        fasta_content = ">gene1\nATGGTTTCTAAAGGTGAA\n"
        fasta_path = tmp_path / "single.fasta"
        fasta_path.write_text(fasta_content)

        results = from_seqio(str(fasta_path), format="fasta")
        assert len(results) == 1
        expected_keys = {"sequence", "organism", "gene_name", "exon_boundaries",
                         "protein", "certificate", "features", "gc_content"}
        assert expected_keys.issubset(results[0].keys())

    def test_from_seqrecord_iterator(self):
        """from_seqio accepts an iterator of SeqRecord objects."""
        from biocompiler.biopython_compat import from_seqio

        records = [
            SeqRecord(Seq("ATGGTTTCT"), id="r1"),
            SeqRecord(Seq("ATGGCCAAA"), id="r2"),
        ]
        results = from_seqio(iter(records))
        assert len(results) == 2
        assert results[0]["sequence"] == "ATGGTTTCT"
        assert results[1]["sequence"] == "ATGGCCAAA"

    def test_from_single_seqrecord(self):
        """from_seqio accepts a single SeqRecord object."""
        from biocompiler.biopython_compat import from_seqio

        rec = SeqRecord(Seq("ATGGTTTCT"), id="single")
        results = from_seqio(rec)
        assert len(results) == 1
        assert results[0]["sequence"] == "ATGGTTTCT"


# ---------------------------------------------------------------------------
# 4. to_genbank_string — export to GenBank via Bio.SeqIO
# ---------------------------------------------------------------------------


class TestToGenbankString:
    """Tests for ``to_genbank_string`` — GenBank export via Bio.SeqIO."""

    def test_returns_string(self, sample_dna):
        """to_genbank_string returns a string."""
        from biocompiler.biopython_compat import to_seqrecord, to_genbank_string

        rec = to_seqrecord(sample_dna, gene_name="testGene")
        gb = to_genbank_string(rec)
        assert isinstance(gb, str)

    def test_genbank_format_markers(self, sample_dna):
        """Output contains GenBank format markers (LOCUS, ORIGIN, //)."""
        from biocompiler.biopython_compat import to_seqrecord, to_genbank_string

        rec = to_seqrecord(sample_dna, gene_name="testGene")
        gb = to_genbank_string(rec)
        assert "LOCUS" in gb
        assert "ORIGIN" in gb
        assert gb.strip().endswith("//")

    def test_sequence_in_output(self, sample_dna):
        """The DNA sequence appears in the GenBank output (lowercased in ORIGIN)."""
        from biocompiler.biopython_compat import to_seqrecord, to_genbank_string

        rec = to_seqrecord(sample_dna, gene_name="testGene")
        gb = to_genbank_string(rec)
        # GenBank ORIGIN section lowercases the sequence, so check lowercase
        gb_no_spaces = gb.replace(" ", "").replace("\n", "").lower()
        assert sample_dna.upper().lower() in gb_no_spaces

    def test_organism_in_output(self, sample_dna):
        """Organism name appears in the GenBank output."""
        from biocompiler.biopython_compat import to_seqrecord, to_genbank_string

        rec = to_seqrecord(sample_dna, organism="Escherichia_coli", gene_name="test")
        gb = to_genbank_string(rec)
        assert "Escherichia coli" in gb or "Escherichia_coli" in gb

    def test_molecule_type_auto_added(self, sample_dna):
        """If molecule_type is missing, it is automatically set to DNA."""
        from biocompiler.biopython_compat import to_genbank_string

        rec = SeqRecord(Seq(sample_dna.upper()), id="test")
        rec.annotations["organism"] = "Homo_sapiens"
        # No molecule_type set — to_genbank_string should add it
        gb = to_genbank_string(rec)
        assert "DNA" in gb


# ---------------------------------------------------------------------------
# 5. to_fasta_string — export to FASTA via Bio.SeqIO
# ---------------------------------------------------------------------------


class TestToFastaString:
    """Tests for ``to_fasta_string`` — FASTA export via Bio.SeqIO."""

    def test_returns_string(self, sample_dna):
        """to_fasta_string returns a string."""
        from biocompiler.biopython_compat import to_seqrecord, to_fasta_string

        rec = to_seqrecord(sample_dna, gene_name="testGene")
        fasta = to_fasta_string(rec)
        assert isinstance(fasta, str)

    def test_fasta_header(self, sample_dna):
        """Output starts with a FASTA header line (>)."""
        from biocompiler.biopython_compat import to_seqrecord, to_fasta_string

        rec = to_seqrecord(sample_dna, gene_name="testGene")
        fasta = to_fasta_string(rec)
        assert fasta.startswith(">")
        assert "testGene" in fasta

    def test_sequence_in_output(self, sample_dna):
        """The DNA sequence appears in the FASTA output."""
        from biocompiler.biopython_compat import to_seqrecord, to_fasta_string

        rec = to_seqrecord(sample_dna, gene_name="testGene")
        fasta = to_fasta_string(rec)
        # FASTA wraps at 60 chars, so remove newlines
        fasta_no_newlines = fasta.replace("\n", "")
        assert sample_dna.upper() in fasta_no_newlines


# ---------------------------------------------------------------------------
# 6. to_seqrecord(result=...) — accept OptimizationResult
# ---------------------------------------------------------------------------


class TestToSeqrecordWithResult:
    """Tests for to_seqrecord accepting an OptimizationResult object."""

    @pytest.fixture(autouse=True)
    def _check_deps(self):
        pytest.importorskip("biocompiler.optimization")

    def test_with_optimization_result(self, sample_protein):
        """to_seqrecord(result=opt_result) produces a valid SeqRecord."""
        from biocompiler.biopython_compat import to_seqrecord
        from biocompiler.optimization import optimize_sequence

        opt = optimize_sequence(sample_protein, organism="Escherichia_coli")
        rec = to_seqrecord(result=opt)
        assert isinstance(rec, SeqRecord)
        assert len(rec.seq) >= 3 * len(sample_protein)

    def test_organism_from_result(self, sample_protein):
        """Organism is extracted from the OptimizationResult when available."""
        from biocompiler.biopython_compat import to_seqrecord
        from biocompiler.optimization import optimize_sequence

        opt = optimize_sequence(sample_protein, organism="Escherichia_coli")
        rec = to_seqrecord(result=opt, organism="Escherichia_coli")
        assert rec.annotations["organism"] == "Escherichia_coli"

    def test_type_error_without_sequence_or_result(self):
        """to_seqrecord raises TypeError when neither sequence nor result is given."""
        from biocompiler.biopython_compat import to_seqrecord

        with pytest.raises(TypeError, match="sequence.*result"):
            to_seqrecord()

    def test_backward_compat_positional_sequence(self, sample_dna):
        """to_seqrecord('ATG...') still works (backward compatibility)."""
        from biocompiler.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert isinstance(rec, SeqRecord)
        assert str(rec.seq) == sample_dna.upper()


# ---------------------------------------------------------------------------
# 7. Full round-trip: SeqRecord → optimize → export → re-import
# ---------------------------------------------------------------------------


class TestFullRoundTrip:
    """End-to-end round-trip tests: import → optimize → export → re-import."""

    @pytest.fixture(autouse=True)
    def _check_deps(self):
        pytest.importorskip("biocompiler.optimization")
        pytest.importorskip("biocompiler.type_system")

    def test_seqrecord_roundtrip_preserves_protein(self, sample_protein):
        """Optimizing and re-importing preserves the protein sequence."""
        from biocompiler.biopython_compat import (
            optimize_seqrecord, from_seqrecord, to_fasta_string,
        )

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"

        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"translation": [sample_protein]},
        )
        rec.features.append(cds)

        optimized = optimize_seqrecord(rec, organism="Escherichia_coli")

        # Re-import the optimized record
        reimported = from_seqrecord(optimized)

        # The protein should match the original
        assert reimported["protein"] == sample_protein

    def test_fasta_roundtrip(self, sample_protein):
        """Optimized SeqRecord → FASTA string → parse back → sequence matches."""
        from biocompiler.biopython_compat import (
            optimize_seqrecord, to_fasta_string, from_seqio,
        )

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"

        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"translation": [sample_protein]},
        )
        rec.features.append(cds)

        optimized = optimize_seqrecord(rec, organism="Escherichia_coli")
        fasta_str = to_fasta_string(optimized)

        # Write to temp file and re-import
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.fasta', delete=False
        ) as f:
            f.write(fasta_str)
            tmp_path = f.name

        try:
            reimported = from_seqio(tmp_path, format="fasta")
            assert len(reimported) >= 1
            # The sequence should be the same as the optimized one
            assert reimported[0]["sequence"] == str(optimized.seq).upper()
        finally:
            os.unlink(tmp_path)

    def test_genbank_roundtrip(self, sample_protein):
        """Optimized SeqRecord → GenBank string → parse back → sequence matches."""
        from biocompiler.biopython_compat import (
            optimize_seqrecord, to_genbank_string, from_seqio,
        )

        rec = SeqRecord(Seq("ATGGTTTCTAAAGGTGAA"), id="test")
        rec.annotations["organism"] = "Escherichia_coli"
        rec.annotations["molecule_type"] = "DNA"

        cds = SeqFeature(
            FeatureLocation(0, 18), type="CDS",
            qualifiers={"translation": [sample_protein]},
        )
        rec.features.append(cds)

        optimized = optimize_seqrecord(rec, organism="Escherichia_coli")
        gb_str = to_genbank_string(optimized)

        # Write to temp file and re-import
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.gbk', delete=False
        ) as f:
            f.write(gb_str)
            tmp_path = f.name

        try:
            reimported = from_seqio(tmp_path, format="genbank")
            assert len(reimported) >= 1
            # The sequence should match
            assert reimported[0]["sequence"] == str(optimized.seq).upper()
        finally:
            os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# 8. Graceful degradation when BioPython is not installed
# ---------------------------------------------------------------------------


class TestRoundTripGracefulDegradation:
    """All new round-trip functions raise ImportError without BioPython."""

    @staticmethod
    def _make_biopython_unavailable():
        """Return a context manager that makes ``import Bio`` fail."""
        import builtins
        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "Bio" or name.startswith("Bio."):
                raise ImportError(f"No module named '{name}'")
            return real_import(name, *args, **kwargs)

        return mock.patch("builtins.__import__", _fake_import)

    def test_optimize_seqrecord_raises_importerror(self):
        """optimize_seqrecord raises ImportError without BioPython."""
        import biocompiler.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                rec = SeqRecord(Seq("ATG"), id="test")
                mod.optimize_seqrecord(rec)

    def test_from_seqio_raises_importerror(self):
        """from_seqio raises ImportError without BioPython."""
        import biocompiler.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.from_seqio("file.fasta")

    def test_to_genbank_string_raises_importerror(self):
        """to_genbank_string raises ImportError without BioPython."""
        import biocompiler.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                rec = SeqRecord(Seq("ATG"), id="test")
                mod.to_genbank_string(rec)

    def test_to_fasta_string_raises_importerror(self):
        """to_fasta_string raises ImportError without BioPython."""
        import biocompiler.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                rec = SeqRecord(Seq("ATG"), id="test")
                mod.to_fasta_string(rec)


# ---------------------------------------------------------------------------
# 9. Module-level checks (no BioPython required)
# ---------------------------------------------------------------------------


class TestRoundTripModuleChecks:
    """Verify new public names are exposed by the module (no BioPython needed)."""

    def test_optimize_seqrecord_in_all(self):
        """optimize_seqrecord is in __all__."""
        import biocompiler.biopython_compat as mod
        assert "optimize_seqrecord" in mod.__all__

    def test_from_seqio_in_all(self):
        """from_seqio is in __all__."""
        import biocompiler.biopython_compat as mod
        assert "from_seqio" in mod.__all__

    def test_to_genbank_string_in_all(self):
        """to_genbank_string is in __all__."""
        import biocompiler.biopython_compat as mod
        assert "to_genbank_string" in mod.__all__

    def test_to_fasta_string_in_all(self):
        """to_fasta_string is in __all__."""
        import biocompiler.biopython_compat as mod
        assert "to_fasta_string" in mod.__all__

    def test_all_new_functions_callable(self):
        """All new functions are callable."""
        import biocompiler.biopython_compat as mod
        assert callable(mod.optimize_seqrecord)
        assert callable(mod.from_seqio)
        assert callable(mod.to_genbank_string)
        assert callable(mod.to_fasta_string)
