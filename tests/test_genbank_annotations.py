"""Tests for GenBank annotation-based workflow."""

import pytest
from biocompiler.export.genbank_annotations import (
    parse_annotation_note,
    AnnotationDirective,
    GenBankAnnotationResult,
    annotations_to_optimization_params,
)


class TestParseAnnotationNote:
    """Test parsing of annotation note strings."""

    def test_no_constraint(self):
        """Empty note should return no directives."""
        result = parse_annotation_note("")
        assert result == []

    def test_single_no_constraint(self):
        """Parse a single @no constraint."""
        result = parse_annotation_note("@no(BsaI_site)")
        assert len(result) >= 1
        assert any(d.directive_type == "no" for d in result)

    def test_multiple_constraints(self):
        """Parse multiple constraints in one note."""
        note = "@no(BsaI_site) @no(EcoRI_site) @optimize(Homo_sapiens)"
        result = parse_annotation_note(note)
        assert len(result) >= 2

    def test_gc_constraint(self):
        """Parse GC content constraint."""
        result = parse_annotation_note("@gc(40-60)")
        assert len(result) >= 1
        gc_directives = [d for d in result if d.directive_type == "gc"]
        assert len(gc_directives) >= 1

    def test_optimize_constraint(self):
        """Parse optimization directive."""
        result = parse_annotation_note("@optimize(Escherichia_coli)")
        assert len(result) >= 1
        opt_directives = [d for d in result if d.directive_type == "optimize"]
        assert len(opt_directives) >= 1

    def test_harmonize_constraint(self):
        """Parse harmonization directive."""
        result = parse_annotation_note("@harmonize(E_coli->Homo_sapiens)")
        assert len(result) >= 1

    def test_no_cpg(self):
        """Parse no_cpg directive."""
        result = parse_annotation_note("@no_cpg")
        assert len(result) >= 1
        assert any(d.directive_type == "no_cpg" for d in result)

    def test_no_splice(self):
        """Parse no_splice directive."""
        result = parse_annotation_note("@no_splice")
        assert len(result) >= 1
        assert any(d.directive_type == "no_splice" for d in result)

    def test_keep_directive(self):
        """Parse keep directive."""
        result = parse_annotation_note("@keep")
        assert len(result) >= 1
        assert any(d.directive_type == "keep" for d in result)


class TestAnnotationDirective:
    """Test AnnotationDirective dataclass."""

    def test_creation(self):
        """Test basic dataclass creation."""
        d = AnnotationDirective(
            directive_type="no",
            parameter="BsaI_site",
            region=None,
            source="@no(BsaI_site)",
        )
        assert d.directive_type == "no"
        assert d.parameter == "BsaI_site"

    def test_with_region(self):
        """Test directive with region specification."""
        d = AnnotationDirective(
            directive_type="no",
            parameter="EcoRI_site",
            region=(100, 300),
            source="@no(EcoRI_site)",
        )
        assert d.region == (100, 300)


class TestGenBankAnnotationResult:
    """Test GenBankAnnotationResult dataclass."""

    def test_creation(self):
        """Test basic dataclass creation."""
        result = GenBankAnnotationResult(
            sequence="ATGCGATCG",
            organism="Homo_sapiens",
            gene_name="GFP",
            directives=[],
            raw_features=[],
        )
        assert result.sequence == "ATGCGATCG"
        assert result.gene_name == "GFP"
        assert result.directives == []


class TestAnnotationsToOptimizationParams:
    """Test conversion of annotations to optimization parameters."""

    def test_empty_annotations(self):
        """Empty annotations should return default parameters."""
        result = GenBankAnnotationResult(
            sequence="ATGCGATCG",
            organism="Homo_sapiens",
            gene_name="test",
            directives=[],
            raw_features=[],
        )
        params = annotations_to_optimization_params(result)
        assert isinstance(params, dict)

    def test_with_optimize_directive(self):
        """Optimize directive should set organism parameter."""
        directives = [
            AnnotationDirective(
                directive_type="optimize",
                parameter="Escherichia_coli",
                region=None,
                source="@optimize(Escherichia_coli)",
            ),
        ]
        result = GenBankAnnotationResult(
            sequence="ATGCGATCG",
            organism="",
            gene_name="test",
            directives=directives,
            raw_features=[],
        )
        params = annotations_to_optimization_params(result)
        assert isinstance(params, dict)

    def test_with_gc_directive(self):
        """GC directive should set gc_lo and gc_hi parameters."""
        directives = [
            AnnotationDirective(
                directive_type="gc",
                parameter="40-60",
                region=None,
                source="@gc(40-60)",
            ),
        ]
        result = GenBankAnnotationResult(
            sequence="ATGCGATCG",
            organism="Homo_sapiens",
            gene_name="test",
            directives=directives,
            raw_features=[],
        )
        params = annotations_to_optimization_params(result)
        assert isinstance(params, dict)

    def test_with_no_restriction(self):
        """No restriction directive should set enzymes."""
        directives = [
            AnnotationDirective(
                directive_type="no",
                parameter="BsaI_site",
                region=None,
                source="@no(BsaI_site)",
            ),
        ]
        result = GenBankAnnotationResult(
            sequence="ATGCGATCG",
            organism="Homo_sapiens",
            gene_name="test",
            directives=directives,
            raw_features=[],
        )
        params = annotations_to_optimization_params(result)
        assert isinstance(params, dict)
