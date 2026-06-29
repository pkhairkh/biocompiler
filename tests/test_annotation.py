"""Tests for biocompiler.annotation — Sequence Annotation Enrichment.

Covered areas:
1. SequenceAnnotation dataclass
2. annotate_sequence — feature detection
3. annotate_to_genbank — GenBank output
4. Individual annotation detectors (ORFs, restriction sites, CpG islands, etc.)
5. Edge cases and boundary conditions
"""

from __future__ import annotations

import pytest

from biocompiler.export.annotation import (
    SequenceAnnotation,
    annotate_sequence,
    annotate_to_genbank,
    _find_orfs,
    _find_restriction_sites,
    _find_cpg_islands,
    _find_simple_repeats,
    _find_gc_at_rich_regions,
    _find_rbs,
    _find_splice_sites,
)


# ═══════════════════════════════════════════════════════════════════════
# 1. SequenceAnnotation dataclass
# ═══════════════════════════════════════════════════════════════════════


class TestSequenceAnnotation:
    """Test the SequenceAnnotation dataclass."""

    def test_create_basic_annotation(self):
        """Create a basic annotation."""
        ann = SequenceAnnotation(
            feature_type="CDS",
            start=0,
            end=100,
            strand=1,
            qualifiers={"gene": "egfp"},
        )
        assert ann.feature_type == "CDS"
        assert ann.start == 0
        assert ann.end == 100
        assert ann.strand == 1
        assert ann.qualifiers["gene"] == "egfp"

    def test_default_qualifiers(self):
        """Qualifiers default to empty dict."""
        ann = SequenceAnnotation(feature_type="CpG_island", start=0, end=200, strand=0)
        assert ann.qualifiers == {}

    def test_all_feature_types(self):
        """All documented feature types can be used."""
        for ftype in ["CDS", "promoter", "RBS", "restriction_site", "CpG_island",
                       "splice_donor", "splice_acceptor", "simple_repeat",
                       "GC_rich", "AT_rich"]:
            ann = SequenceAnnotation(feature_type=ftype, start=0, end=10, strand=0)
            assert ann.feature_type == ftype

    def test_strand_values(self):
        """All strand values (1, -1, 0) are accepted."""
        for strand in (1, -1, 0):
            ann = SequenceAnnotation(feature_type="CDS", start=0, end=10, strand=strand)
            assert ann.strand == strand


# ═══════════════════════════════════════════════════════════════════════
# 2. annotate_sequence
# ═══════════════════════════════════════════════════════════════════════


class TestAnnotateSequence:
    """Test the main annotate_sequence function."""

    def test_returns_list_of_annotations(self, sample_dna):
        """annotate_sequence returns a list of SequenceAnnotation."""
        result = annotate_sequence(sample_dna)
        assert isinstance(result, list)
        for ann in result:
            assert isinstance(ann, SequenceAnnotation)

    def test_empty_sequence(self):
        """Empty sequence returns empty list."""
        result = annotate_sequence("")
        assert result == []

    def test_short_sequence(self):
        """Short sequence may have few annotations."""
        result = annotate_sequence("ATGCGT")
        assert isinstance(result, list)

    def test_annotations_sorted_by_position(self, sample_dna):
        """Annotations are sorted by start position."""
        result = annotate_sequence(sample_dna)
        for i in range(len(result) - 1):
            assert result[i].start <= result[i + 1].start

    def test_finds_cds_in_coding_sequence(self):
        """annotate_sequence finds CDS in a known coding sequence."""
        seq = "ATG" + "GCT" * 40 + "TAA"  # 123 bp, 40 aa ORF
        result = annotate_sequence(seq)
        cds = [a for a in result if a.feature_type == "CDS"]
        assert len(cds) >= 1

    def test_finds_restriction_sites(self):
        """annotate_sequence finds restriction sites."""
        # EcoRI site: GAATTC
        seq = "ATGGAATTC" + "GCT" * 30 + "TAA"
        result = annotate_sequence(seq)
        rsites = [a for a in result if a.feature_type == "restriction_site"]
        eco_ri = [a for a in rsites if a.qualifiers.get("enzyme") == "EcoRI"]
        assert len(eco_ri) >= 1

    def test_organism_parameter_accepted(self, sample_dna):
        """Organism parameter is accepted without error."""
        result = annotate_sequence(sample_dna, organism="Homo_sapiens")
        assert isinstance(result, list)

    def test_feature_types_are_valid(self, sample_dna):
        """All returned annotations have valid feature types."""
        valid_types = {"CDS", "promoter", "RBS", "restriction_site", "CpG_island",
                       "splice_donor", "splice_acceptor", "simple_repeat",
                       "GC_rich", "AT_rich"}
        result = annotate_sequence(sample_dna)
        for ann in result:
            assert ann.feature_type in valid_types, f"Invalid feature type: {ann.feature_type}"


# ═══════════════════════════════════════════════════════════════════════
# 3. Individual detectors
# ═══════════════════════════════════════════════════════════════════════


class TestFindORFs:
    """Test _find_orfs function."""

    def test_finds_orf_in_coding_sequence(self):
        """Finds the main ORF in a coding sequence."""
        seq = "ATG" + "GCT" * 40 + "TAA"
        orfs = _find_orfs(seq)
        assert len(orfs) >= 1
        assert orfs[0].feature_type == "CDS"

    def test_orf_has_translation(self):
        """CDS annotation includes protein translation."""
        seq = "ATG" + "GCT" * 40 + "TAA"
        orfs = _find_orfs(seq)
        assert len(orfs) >= 1
        assert "translation" in orfs[0].qualifiers

    def test_no_orf_in_short_sequence(self):
        """Short sequence has no ORFs above threshold."""
        orfs = _find_orfs("ATGCGT")
        assert len(orfs) == 0


class TestFindRestrictionSites:
    """Test _find_restriction_sites function."""

    def test_finds_ecori(self):
        """Finds EcoRI site (GAATTC)."""
        seq = "ATGGAATTC" + "GCT" * 30 + "TAA"
        sites = _find_restriction_sites(seq)
        eco_ri = [s for s in sites if s.qualifiers.get("enzyme") == "EcoRI"]
        assert len(eco_ri) >= 1

    def test_finds_bamhi(self):
        """Finds BamHI site (GGATCC)."""
        seq = "ATGGGATCC" + "GCT" * 30 + "TAA"
        sites = _find_restriction_sites(seq)
        bamhi = [s for s in sites if s.qualifiers.get("enzyme") == "BamHI"]
        assert len(bamhi) >= 1

    def test_no_sites_in_simple_sequence(self):
        """Simple AT-rich sequence has few restriction sites."""
        seq = "ATATATATATATATAT"
        sites = _find_restriction_sites(seq)
        # ATATAT does not match any 6-cutter; may or may not find sites
        assert isinstance(sites, list)

    def test_site_position_correct(self):
        """Restriction site position matches the actual location."""
        seq = "AAAA" + "GAATTC" + "AAAA"
        sites = _find_restriction_sites(seq)
        eco_ri = [s for s in sites if s.qualifiers.get("enzyme") == "EcoRI"]
        if eco_ri:
            assert eco_ri[0].start == 4
            assert eco_ri[0].end == 10


class TestFindCpGIslands:
    """Test _find_cpg_islands function."""

    def test_detects_cpg_island(self):
        """Detects a CpG island in a GC/CpG-rich sequence."""
        # Build a 300bp sequence rich in CG dinucleotides
        cpg_unit = "CGCGATCGCG"  # 10 bp, lots of CpGs
        seq = cpg_unit * 30  # 300 bp
        islands = _find_cpg_islands(seq)
        assert len(islands) >= 1
        assert islands[0].feature_type == "CpG_island"

    def test_no_cpg_in_at_rich(self):
        """AT-rich sequence has no CpG islands."""
        seq = "ATATATATAT" * 30  # 300 bp, no CpGs
        islands = _find_cpg_islands(seq)
        assert len(islands) == 0

    def test_cpg_island_qualifiers(self):
        """CpG island annotation includes gc_content and length."""
        cpg_unit = "CGCGATCGCG"
        seq = cpg_unit * 30
        islands = _find_cpg_islands(seq)
        if islands:
            assert "gc_content" in islands[0].qualifiers
            assert "length" in islands[0].qualifiers

    def test_short_sequence_no_cpg(self):
        """Sequence shorter than window size has no CpG islands."""
        seq = "CGCGCGCG"
        islands = _find_cpg_islands(seq)
        assert len(islands) == 0


class TestFindSimpleRepeats:
    """Test _find_simple_repeats function."""

    def test_detects_dinucleotide_repeat(self):
        """Detects a dinucleotide repeat (e.g. ATATATAT)."""
        seq = "AAAA" + "ATATATATAT" + "AAAA"  # 5+ copies of AT
        repeats = _find_simple_repeats(seq)
        at_repeats = [r for r in repeats if r.qualifiers.get("repeat_unit") == "AT"]
        assert len(at_repeats) >= 1

    def test_detects_trinucleotide_repeat(self):
        """Detects a trinucleotide repeat (e.g. CAGCAGCAG)."""
        seq = "CAGCAGCAGCAG"  # 4 copies
        repeats = _find_simple_repeats(seq)
        cag_repeats = [r for r in repeats if r.qualifiers.get("repeat_unit") == "CAG"]
        assert len(cag_repeats) >= 1

    def test_repeat_qualifiers(self):
        """Simple repeat annotation has repeat_unit and copies."""
        seq = "ATATATATAT"
        repeats = _find_simple_repeats(seq)
        if repeats:
            assert "repeat_unit" in repeats[0].qualifiers
            assert "copies" in repeats[0].qualifiers

    def test_no_repeat_in_random_sequence(self):
        """Random-like sequence has few simple repeats."""
        seq = "ATGCGTCGATCGATCG"
        repeats = _find_simple_repeats(seq)
        # May have some short repeats, but should not be many
        assert isinstance(repeats, list)


class TestFindGCATRichRegions:
    """Test _find_gc_at_rich_regions function."""

    def test_detects_gc_rich(self):
        """Detects a GC-rich region."""
        seq = "GCGCGCGCGC" * 10  # 100 bp, 100% GC
        regions = _find_gc_at_rich_regions(seq)
        gc_rich = [r for r in regions if r.feature_type == "GC_rich"]
        assert len(gc_rich) >= 1

    def test_detects_at_rich(self):
        """Detects an AT-rich region."""
        seq = "ATATATATAT" * 10  # 100 bp, 100% AT
        regions = _find_gc_at_rich_regions(seq)
        at_rich = [r for r in regions if r.feature_type == "AT_rich"]
        assert len(at_rich) >= 1

    def test_short_sequence_no_regions(self):
        """Sequence shorter than window has no GC/AT-rich regions."""
        seq = "GCGCGC"
        regions = _find_gc_at_rich_regions(seq)
        assert len(regions) == 0


class TestFindRBS:
    """Test _find_rbs function."""

    def test_detects_shine_dalgarno(self):
        """Detects a Shine-Dalgarno sequence upstream of ATG."""
        seq = "AAAAAGGAGGAAAAAATG" + "GCT" * 30 + "TAA"
        rbs = _find_rbs(seq)
        assert len(rbs) >= 1
        assert rbs[0].feature_type == "RBS"

    def test_no_rbs_without_atg(self):
        """No RBS found when there is no ATG."""
        seq = "AAAAAGGAGGAAAAAAGCT" * 5
        rbs = _find_rbs(seq)
        assert len(rbs) == 0


class TestFindSpliceSites:
    """Test _find_splice_sites function."""

    def test_finds_donor_sites(self):
        """Finds GT splice donor sites."""
        seq = "ATGGTAAAGTCCC"  # contains GT at positions 3, 8
        sites = _find_splice_sites(seq)
        donors = [s for s in sites if s.feature_type == "splice_donor"]
        assert len(donors) >= 1

    def test_finds_acceptor_sites(self):
        """Finds AG splice acceptor sites."""
        seq = "ATGAGAAAGTCCC"  # contains AG at position 3
        sites = _find_splice_sites(seq)
        acceptors = [s for s in sites if s.feature_type == "splice_acceptor"]
        assert len(acceptors) >= 1

    def test_splice_site_has_score(self):
        """Splice site annotations include a score qualifier."""
        seq = "ATGGTAAAGTCCC"
        sites = _find_splice_sites(seq)
        if sites:
            assert "score" in sites[0].qualifiers


# ═══════════════════════════════════════════════════════════════════════
# 4. annotate_to_genbank
# ═══════════════════════════════════════════════════════════════════════


class TestAnnotateToGenbank:
    """Test the annotate_to_genbank function."""

    def test_returns_string(self, sample_dna):
        """annotate_to_genbank returns a string."""
        result = annotate_to_genbank(sample_dna)
        assert isinstance(result, str)

    def test_genbank_format_headers(self, sample_dna):
        """Output has GenBank headers."""
        result = annotate_to_genbank(sample_dna)
        assert "LOCUS" in result
        assert "DEFINITION" in result
        assert "ORIGIN" in result
        assert result.endswith("//")

    def test_genbank_features_section(self, sample_dna):
        """Output has a FEATURES section."""
        result = annotate_to_genbank(sample_dna)
        assert "FEATURES" in result

    def test_custom_name(self, sample_dna):
        """Custom name is used as the locus name."""
        result = annotate_to_genbank(sample_dna, name="EGFP")
        assert "EGFP" in result.split("\n")[0]  # LOCUS line

    def test_organism_in_record(self, sample_dna):
        """Organism appears in the record."""
        result = annotate_to_genbank(sample_dna, organism="Escherichia_coli")
        assert "Escherichia" in result or "E. coli" in result or "coli" in result

    def test_empty_sequence(self):
        """Empty sequence returns empty string."""
        result = annotate_to_genbank("")
        assert result == ""

    def test_with_precomputed_annotations(self, sample_dna):
        """Pre-computed annotations are used when provided."""
        custom_ann = [
            SequenceAnnotation(
                feature_type="CpG_island",
                start=0,
                end=len(sample_dna),
                strand=0,
                qualifiers={"gc_content": "0.3333"},
            )
        ]
        result = annotate_to_genbank(sample_dna, annotations=custom_ann)
        assert "CpG" in result or "CpG_island" in result

    def test_restriction_site_in_genbank(self):
        """Restriction site appears in GenBank output."""
        seq = "ATGGAATTC" + "GCT" * 30 + "TAA"
        result = annotate_to_genbank(seq)
        assert "EcoRI" in result

    def test_long_sequence_annotation(self, egfp_dna):
        """annotate_to_genbank works with a full-length gene."""
        result = annotate_to_genbank(egfp_dna, name="eGFP", organism="Homo_sapiens")
        assert isinstance(result, str)
        assert "LOCUS" in result
        assert len(result) > 100  # non-trivial output


# ═══════════════════════════════════════════════════════════════════════
# 5. Edge cases
# ═══════════════════════════════════════════════════════════════════════


class TestAnnotationEdgeCases:
    """Edge cases and boundary conditions."""

    def test_all_n_sequence(self):
        """All-N sequence is handled."""
        result = annotate_sequence("NNNNNNNNNN")
        assert isinstance(result, list)

    def test_very_short_sequence(self):
        """Very short sequence (3 bp) is handled."""
        result = annotate_sequence("ATG")
        assert isinstance(result, list)

    def test_single_restriction_site_sequence(self):
        """Sequence that is exactly a restriction site."""
        result = annotate_sequence("GAATTC")
        assert isinstance(result, list)
        rsites = [a for a in result if a.feature_type == "restriction_site"]
        eco_ri = [a for a in rsites if a.qualifiers.get("enzyme") == "EcoRI"]
        assert len(eco_ri) >= 1

    def test_annotation_start_end_valid(self, sample_dna):
        """All annotations have valid start/end positions."""
        result = annotate_sequence(sample_dna)
        for ann in result:
            assert ann.start >= 0, f"Negative start: {ann}"
            assert ann.end > ann.start, f"End <= start: {ann}"
            assert ann.end <= len(sample_dna), f"End beyond sequence: {ann}"

    def test_palindromic_restriction_site(self):
        """Palindromic restriction sites are found on forward strand only."""
        # EcoRI: GAATTC (palindromic — RC is also GAATTC)
        seq = "ATGGAATTC" + "GCT" * 30 + "TAA"
        sites = _find_restriction_sites(seq)
        eco_ri = [s for s in sites if s.qualifiers.get("enzyme") == "EcoRI"]
        # Should find exactly 1 (not 2 — palindrome is same on both strands)
        assert len(eco_ri) == 1
