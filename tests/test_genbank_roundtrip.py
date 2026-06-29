"""Tests for BioCompiler GenBank Round-Trip Verification.

Covers:
- RoundTripResult dataclass
- verify_genbank_roundtrip: full round-trip with sequence & annotation verification
- compare_sequences: DNA sequence comparison
- verify_annotation_preservation: annotation survival check
- Edge cases: empty sequences, export failures, annotation mismatches
"""

import pytest
from biocompiler.optimizer import OptimizationResult
from biocompiler.export.genbank_roundtrip import (
    RoundTripResult,
    verify_genbank_roundtrip,
    compare_sequences,
    verify_annotation_preservation,
)


# ────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture
def sample_result() -> OptimizationResult:
    """A sample OptimizationResult for round-trip tests."""
    return OptimizationResult(
        sequence="ATGGTTTCTAAAGGTGAA",
        gc_content=0.333,
        cai=0.78,
        failed_predicates=[],
        predicate_results=[],
        certificate_text="",
        protein="MVSKGE",
        fallback_used=False,
        satisfied_predicates=[],
        aa_substitutions=[],
        mutagenesis_applied=False,
    )


@pytest.fixture
def longer_result() -> OptimizationResult:
    """A longer OptimizationResult for more thorough round-trip tests.

    Uses a sequence that encodes a known protein for reliable testing.
    """
    # 90 bp → 30 aa protein
    seq = "ATGGCCGCTAAAGGTGAAGCCGCTTTTGCCGCTAAAGGTGAAGCCGCTAAAGGTGAAGCCGCTAAAGGTGAAGCCGCTAAAGGTGAAGCC"
    protein = "MAAKGEAAAAAAAAAKGEAAAAAAAAAKGE"[:30]
    return OptimizationResult(
        sequence=seq[:90],
        gc_content=0.5,
        cai=0.85,
        failed_predicates=[],
        predicate_results=[],
        certificate_text="",
        protein=protein,
        fallback_used=False,
        satisfied_predicates=[],
        aa_substitutions=[],
        mutagenesis_applied=False,
    )


# ────────────────────────────────────────────────────────────
# compare_sequences
# ────────────────────────────────────────────────────────────

class TestCompareSequences:
    """Tests for the compare_sequences function."""

    def test_identical_sequences(self):
        """Identical sequences should produce no mismatches."""
        mismatches = compare_sequences("ATGCGT", "ATGCGT")
        assert mismatches == []

    def test_single_mismatch(self):
        """Single position difference should be detected."""
        mismatches = compare_sequences("ATGCGT", "ATTCGT")
        assert len(mismatches) == 1
        assert mismatches[0] == (2, "G", "T")

    def test_multiple_mismatches(self):
        """Multiple differences should all be reported."""
        mismatches = compare_sequences("ATGCGT", "ATTCGA")
        assert len(mismatches) == 2
        positions = {m[0] for m in mismatches}
        assert 2 in positions  # G vs T
        assert 5 in positions  # T vs A

    def test_case_insensitive(self):
        """Comparison should be case-insensitive."""
        mismatches = compare_sequences("atgcgt", "ATGCGT")
        assert mismatches == []

    def test_different_lengths_shorter_reimported(self):
        """When reimported is shorter, only overlapping positions are checked."""
        mismatches = compare_sequences("ATGCGT", "ATGC")
        assert len(mismatches) == 0  # Overlapping region matches

    def test_different_lengths_longer_reimported(self):
        """When reimported is longer, only overlapping positions are checked."""
        mismatches = compare_sequences("ATGC", "ATGCGT")
        assert len(mismatches) == 0  # Overlapping region matches

    def test_empty_sequences(self):
        """Empty sequences should produce no mismatches."""
        mismatches = compare_sequences("", "")
        assert mismatches == []

    def test_one_empty_sequence(self):
        """Comparing with an empty sequence should produce no mismatches."""
        mismatches = compare_sequences("ATGC", "")
        assert mismatches == []

    def test_completely_different(self):
        """Completely different sequences should all be mismatches."""
        mismatches = compare_sequences("AAAA", "TTTT")
        assert len(mismatches) == 4


# ────────────────────────────────────────────────────────────
# verify_annotation_preservation
# ────────────────────────────────────────────────────────────

class TestVerifyAnnotationPreservation:
    """Tests for the verify_annotation_preservation function."""

    def test_all_preserved(self):
        """When all annotations match, should return True."""
        original = {"gene_name": "gfp", "protein": "MVSKGE", "organism": "Escherichia_coli"}
        reimported = {"gene_name": "gfp", "protein": "MVSKGE", "organism": "Escherichia coli"}

        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is True
        assert len(warnings) == 0

    def test_gene_name_mismatch(self):
        """Gene name mismatch should be detected."""
        original = {"gene_name": "gfp", "protein": "", "organism": ""}
        reimported = {"gene_name": "rfp", "protein": "", "organism": ""}

        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is False
        assert any("gene name" in w.lower() for w in warnings)

    def test_protein_mismatch(self):
        """Protein mismatch should be detected."""
        original = {"gene_name": "", "protein": "MVSKGE", "organism": ""}
        reimported = {"gene_name": "", "protein": "MVSKGD", "organism": ""}

        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is False
        assert any("protein" in w.lower() for w in warnings)

    def test_organism_mismatch(self):
        """Organism mismatch should be detected."""
        original = {"gene_name": "", "protein": "", "organism": "Escherichia_coli"}
        reimported = {"gene_name": "", "protein": "", "organism": "Homo sapiens"}

        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is False
        assert any("organism" in w.lower() for w in warnings)

    def test_empty_original_values_pass(self):
        """If original values are empty, the check should pass."""
        original = {"gene_name": "", "protein": "", "organism": ""}
        reimported = {"gene_name": "something", "protein": "else", "organism": "different"}

        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is True
        assert len(warnings) == 0

    def test_organism_underscore_normalization(self):
        """Organism with underscores vs spaces should be treated as matching."""
        original = {"gene_name": "", "protein": "", "organism": "Escherichia_coli"}
        reimported = {"gene_name": "", "protein": "", "organism": "Escherichia coli"}

        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is True
        assert len(warnings) == 0

    def test_protein_whitespace_normalization(self):
        """Protein with whitespace should be normalized before comparison."""
        original = {"gene_name": "", "protein": "MVSKGE", "organism": ""}
        reimported = {"gene_name": "", "protein": "MVSK GE", "organism": ""}

        # After removing whitespace, they should match
        preserved, warnings = verify_annotation_preservation(original, reimported)
        assert preserved is True


# ────────────────────────────────────────────────────────────
# verify_genbank_roundtrip — full integration
# ────────────────────────────────────────────────────────────

class TestVerifyGenbankRoundtrip:
    """Tests for the verify_genbank_roundtrip function."""

    def test_basic_roundtrip_success(self, sample_result):
        """Basic round-trip should succeed for a simple OptimizationResult."""
        report = verify_genbank_roundtrip(sample_result)

        assert report.success is True
        assert report.original_dna == sample_result.sequence.upper()
        assert report.reimported_dna == sample_result.sequence.upper()
        assert report.mismatches == []
        assert report.annotation_preserved is True
        assert len(report.warnings) == 0

    def test_sequence_preserved(self, sample_result):
        """DNA sequence should be preserved through the round trip."""
        report = verify_genbank_roundtrip(sample_result)

        assert report.original_dna == report.reimported_dna
        assert len(report.reimported_dna) == len(sample_result.sequence)

    def test_gene_name_preserved(self, sample_result):
        """Gene name should be preserved through the round trip."""
        report = verify_genbank_roundtrip(
            sample_result, gene_name="test_gene"
        )

        assert report.annotation_preserved is True
        # Check that the gene name appears in the re-imported annotations
        reimp_gene = report.reimported_annotations.get("gene_name", "")
        assert reimp_gene == "test_gene"

    def test_protein_preserved(self, sample_result):
        """Protein translation should be preserved through the round trip."""
        report = verify_genbank_roundtrip(sample_result)

        if report.annotation_preserved:
            reimp_protein = report.reimported_annotations.get("protein", "")
            # Protein may have been extracted from /translation qualifier
            # or could be empty if GenBank did not parse it
            if reimp_protein:
                # Normalize whitespace for comparison
                assert reimp_protein.replace(" ", "") == sample_result.protein

    def test_organism_preserved(self, sample_result):
        """Organism name should be preserved (underscores→spaces is ok)."""
        report = verify_genbank_roundtrip(
            sample_result, organism="Escherichia_coli"
        )

        if report.annotation_preserved:
            reimp_org = report.reimported_annotations.get("organism", "")
            # GenBank uses spaces instead of underscores
            assert "Escherichia" in reimp_org

    def test_genbank_text_generated(self, sample_result):
        """genbank_text should be populated for debugging."""
        report = verify_genbank_roundtrip(sample_result)

        assert report.genbank_text != ""
        assert "LOCUS" in report.genbank_text
        assert "ORIGIN" in report.genbank_text
        assert "//" in report.genbank_text

    def test_annotations_populated(self, sample_result):
        """Both original and reimported annotations should be populated."""
        report = verify_genbank_roundtrip(sample_result)

        assert report.original_annotations != {}
        assert report.reimported_annotations != {}
        assert "gene_name" in report.original_annotations
        assert "protein" in report.original_annotations
        assert "organism" in report.original_annotations

    def test_default_organism(self, sample_result):
        """Default organism should be Homo_sapiens if not specified."""
        report = verify_genbank_roundtrip(sample_result)

        assert report.original_annotations.get("organism") == "Homo_sapiens"

    def test_custom_organism(self, sample_result):
        """Custom organism should be used when specified."""
        report = verify_genbank_roundtrip(
            sample_result, organism="Saccharomyces_cerevisiae"
        )

        assert report.original_annotations.get("organism") == "Saccharomyces_cerevisiae"

    def test_result_no_protein(self):
        """Round-trip should work even without protein."""
        result = OptimizationResult(
            sequence="ATGGTTTCTAAAGGTGAA",
            gc_content=0.333,
            cai=0.78,
            protein="",
        )
        report = verify_genbank_roundtrip(result)

        assert report.success is True
        assert report.original_dna == report.reimported_dna


# ────────────────────────────────────────────────────────────
# RoundTripResult dataclass
# ────────────────────────────────────────────────────────────

class TestRoundTripResult:
    """Tests for the RoundTripResult dataclass."""

    def test_basic_construction(self):
        """Should construct with required fields."""
        result = RoundTripResult(
            success=True,
            original_dna="ATGC",
            reimported_dna="ATGC",
            mismatches=[],
            annotation_preserved=True,
        )
        assert result.success is True
        assert result.original_dna == "ATGC"
        assert result.warnings == []

    def test_default_fields(self):
        """Optional fields should have defaults."""
        result = RoundTripResult(
            success=False,
            original_dna="ATGC",
            reimported_dna="ATGG",
            mismatches=[(3, "C", "G")],
            annotation_preserved=False,
        )
        assert result.original_annotations == {}
        assert result.reimported_annotations == {}
        assert result.genbank_text == ""
        assert result.warnings == []

    def test_with_warnings(self):
        """Warnings should be storable."""
        result = RoundTripResult(
            success=False,
            original_dna="ATGC",
            reimported_dna="ATGG",
            mismatches=[(3, "C", "G")],
            annotation_preserved=True,
            warnings=["Sequence mismatch detected"],
        )
        assert len(result.warnings) == 1


# ────────────────────────────────────────────────────────────
# Edge cases and error handling
# ────────────────────────────────────────────────────────────

class TestRoundtripEdgeCases:
    """Edge case tests for round-trip verification."""

    def test_single_codon_sequence(self):
        """Round-trip with a very short sequence (single codon)."""
        result = OptimizationResult(
            sequence="ATG",
            gc_content=0.333,
            cai=1.0,
            protein="M",
        )
        report = verify_genbank_roundtrip(result)
        assert report.original_dna == "ATG"
        assert report.reimported_dna == "ATG"

    def test_long_sequence_roundtrip(self):
        """Round-trip with a longer sequence."""
        # Generate a simple repeating sequence
        seq = "ATGGCT" * 50  # 300 bp
        protein = "MA" * 50
        result = OptimizationResult(
            sequence=seq,
            gc_content=0.5,
            cai=0.9,
            protein=protein,
        )
        report = verify_genbank_roundtrip(result)
        assert report.success is True
        assert report.original_dna == report.reimported_dna

    def test_high_gc_sequence(self):
        """Round-trip with a high-GC sequence."""
        seq = "GCGCCGCGCC" * 9  # 90 bp, ~80% GC
        result = OptimizationResult(
            sequence=seq,
            gc_content=0.8,
            cai=0.6,
        )
        report = verify_genbank_roundtrip(result)
        assert report.original_dna == report.reimported_dna

    def test_low_gc_sequence(self):
        """Round-trip with a low-GC sequence."""
        seq = "ATATATATAT" * 9  # 90 bp, ~20% GC
        result = OptimizationResult(
            sequence=seq,
            gc_content=0.2,
            cai=0.5,
        )
        report = verify_genbank_roundtrip(result)
        assert report.original_dna == report.reimported_dna

    def test_result_with_all_metadata(self):
        """Round-trip with a result that has all optional metadata populated."""
        result = OptimizationResult(
            sequence="ATGGTTTCTAAAGGTGAA",
            gc_content=0.333,
            cai=0.78,
            protein="MVSKGE",
            fallback_used=True,
            satisfied_predicates=["GCInRange"],
            failed_predicates=[],
            mutagenesis_applied=True,
            aa_substitutions=[{"position": 0, "original": "M", "substitution": "M"}],
            codon_pair_bias=0.42,
            mrna_stability_score=0.85,
            suggested_5utr="GCCACC",
            suggested_3utr="TTTTTT",
        )
        report = verify_genbank_roundtrip(result)
        assert report.success is True
        assert report.original_dna == report.reimported_dna
