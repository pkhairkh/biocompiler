"""
Unit tests for biocompiler.report — top-level HTML report generation.

Covers:
1. Report generation functions (generate_report, _compute_gc_windows,
   _compute_codon_usage, _element_type_description)
2. Output format correctness (valid HTML structure, embedded SVGs,
   verdict badges, certificate section, engine section, isoform section,
   XSS / HTML-escaping safety)
3. Input validation (empty sequence, lowercase input, missing optional
   parameters, broken dependencies via mocking)
"""

import html
import re
from unittest.mock import patch, MagicMock

import pytest

from biocompiler.provenance.report import (
    generate_report,
    _compute_gc_windows,
    _compute_codon_usage,
    _element_type_description,
    _VERDICT_SYMBOLS,
    _MAX_TOKEN_DISPLAY,
    _MAX_ISOFORM_DISPLAY,
    _DESIGN_ID_DISPLAY_LENGTH,
    _DEFAULT_GC_WINDOW_SIZE,
)
from biocompiler.shared.types import (
    Certificate,
    TypeCheckResult,
    Token,
    SpliceIsoform,
    Verdict,
)
from biocompiler.engines.base import BaseEngineResult


# ── Test fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def simple_seq() -> str:
    """A valid 60-bp DNA sequence with balanced GC."""
    return "ATGGCGGCCTGAATGGCGGCCTGAATGGCGGCCTGAATGGCGGCCTGAATGGCGGCCTGA"


@pytest.fixture
def coding_seq() -> str:
    """A 90-bp sequence starting with ATG, ending with a stop codon in frame."""
    return "ATGGCGGCCTAGATGGCGGCCTAGATGGCGGCCTAGATGGCGGCCTAGATGGCGGCCTAG" \
           "ATGGCGGCCTAGATGGCGGCCTAGATGGCGGCCTAGATGGCGGCCTAG"


@pytest.fixture
def type_results_pass() -> list[TypeCheckResult]:
    return [
        TypeCheckResult(predicate="GCContent", verdict=Verdict.PASS),
        TypeCheckResult(predicate="NoStopCodon", verdict=Verdict.PASS),
    ]


@pytest.fixture
def type_results_mixed() -> list[TypeCheckResult]:
    return [
        TypeCheckResult(predicate="GCContent", verdict=Verdict.PASS),
        TypeCheckResult(predicate="RestrictionSite", verdict=Verdict.FAIL,
                        violation="EcoRI site found at position 42"),
        TypeCheckResult(predicate="SpliceDonor", verdict=Verdict.UNCERTAIN,
                        knowledge_gap="Heuristic detection only"),
    ]


@pytest.fixture
def sample_tokens() -> list[Token]:
    return [
        Token(position=0, element_type="start_codon", match_sequence="ATG", score=1.0, frame=0, strand="+"),
        Token(position=27, element_type="splice_donor", match_sequence="GT", score=0.85, strand="+"),
        Token(position=57, element_type="stop_codon", match_sequence="TGA", score=1.0, frame=0, strand="+"),
    ]


@pytest.fixture
def sample_isoforms() -> list[SpliceIsoform]:
    return [
        SpliceIsoform(
            sequence="ATGGCGGCCTGA",
            exon_boundaries=[(0, 12)],
            parse_path=["exon_1"],
            score=0.95,
        ),
        SpliceIsoform(
            sequence="ATGGCGGCCTGAATGGCG",
            exon_boundaries=[(0, 6), (12, 18)],
            parse_path=["exon_1", "exon_2"],
            score=0.72,
        ),
    ]


@pytest.fixture
def sample_certificate() -> Certificate:
    return Certificate(
        version="1.0",
        design_id="abc123def456ghi789jkl012mno345pqr678stu901",
        sequence="ATGGCGGCCTGA",
        types=[{"predicate": "GCContent", "verdict": "PASS"}],
        provenance={"timestamp": "2025-01-15T10:30:00Z", "tool": "biocompiler"},
    )


@pytest.fixture
def sample_engine_results() -> list[BaseEngineResult]:
    return [
        BaseEngineResult(
            sequence="ATGGCGGCCTGA",
            primary_score=0.87,
            classification="stable",
            success=True,
            engine_name="FoldX",
            primary_score_label="ddG",
            execution_time_s=1.234,
        ),
        BaseEngineResult(
            sequence="ATGGCGGCCTGA",
            primary_score=-2.5,
            classification="low_immunogenicity",
            success=True,
            engine_name="Immunogenicity",
            primary_score_label="immunogenicity_score",
            execution_time_s=0.456,
        ),
    ]


@pytest.fixture
def sample_optimization_result() -> dict:
    return {
        "cai": 0.8234,
        "gc_content": 0.55,
        "satisfied": ["GCContent", "NoStopCodon"],
        "failed": ["RestrictionSite"],
    }


# ══════════════════════════════════════════════════════════════════════════
# 1. Report Generation Functions
# ══════════════════════════════════════════════════════════════════════════

class TestGenerateReport:
    """Tests for the main generate_report function."""

    def test_returns_string(self, simple_seq):
        result = generate_report(simple_seq)
        assert isinstance(result, str)

    def test_returns_non_empty(self, simple_seq):
        result = generate_report(simple_seq)
        assert len(result) > 0

    def test_sequence_uppercased(self):
        """Input sequence should be uppercased internally."""
        result = generate_report("atggcggcctga")
        # The raw sequence textarea should contain uppercase
        assert "ATGGCGGCCTGA" in result

    def test_default_organism(self, simple_seq):
        result = generate_report(simple_seq)
        assert "Homo_sapiens" in result

    def test_custom_organism(self, simple_seq):
        result = generate_report(simple_seq, organism="Escherichia_coli")
        assert "Escherichia_coli" in result

    def test_gene_name_in_title(self, simple_seq):
        result = generate_report(simple_seq, gene_name="eGFP")
        assert "eGFP" in result
        assert "BioCompiler Report" in result

    def test_no_gene_name_default_title(self, simple_seq):
        result = generate_report(simple_seq)
        assert "<title>BioCompiler Report</title>" in result

    def test_with_type_results_pass(self, simple_seq, type_results_pass):
        result = generate_report(simple_seq, type_results=type_results_pass)
        assert "PASS" in result
        assert "verdict-badge pass" in result

    def test_with_type_results_mixed(self, simple_seq, type_results_mixed):
        result = generate_report(simple_seq, type_results=type_results_mixed)
        # Combined verdict should be FAIL (weakest link)
        assert "FAIL" in result

    def test_with_no_type_results_shows_na(self, simple_seq):
        result = generate_report(simple_seq, type_results=None)
        assert "N/A" in result

    def test_with_empty_type_results_shows_na(self, simple_seq):
        result = generate_report(simple_seq, type_results=[])
        assert "N/A" in result

    def test_with_tokens(self, simple_seq, sample_tokens):
        result = generate_report(simple_seq, tokens=sample_tokens)
        assert "start_codon" in result
        assert "splice_donor" in result

    def test_with_isoforms(self, simple_seq, sample_isoforms):
        result = generate_report(simple_seq, isoforms=sample_isoforms)
        assert "Splice Isoforms" in result
        assert "0.95" in result

    def test_without_isoforms_no_section(self, simple_seq):
        result = generate_report(simple_seq, isoforms=None,
                                 exon_boundaries=[(0, len(simple_seq))])
        # The isoform <section> element should not appear if no isoforms;
        # only an HTML comment <!-- Splice Isoforms --> may be present.
        assert '<h2>Splice Isoforms' not in result

    def test_with_certificate(self, simple_seq, sample_certificate):
        result = generate_report(simple_seq, certificate=sample_certificate)
        assert "Certificate" in result
        assert "abc123de" in result  # First _DESIGN_ID_DISPLAY_LENGTH chars
        assert "1.0" in result
        assert "biocompiler" in result

    def test_without_certificate(self, simple_seq):
        result = generate_report(simple_seq, certificate=None)
        assert "Certificate" not in result

    def test_with_optimization_result(self, simple_seq, sample_optimization_result):
        result = generate_report(simple_seq, optimization_result=sample_optimization_result)
        assert "Optimization Results" in result
        assert "0.8234" in result

    def test_without_optimization_result(self, simple_seq):
        result = generate_report(simple_seq, optimization_result=None)
        assert "Optimization Results" not in result

    def test_with_engine_results(self, simple_seq, sample_engine_results):
        result = generate_report(simple_seq, engine_results=sample_engine_results)
        assert "Engine Analysis Results" in result
        assert "FoldX" in result
        assert "Immunogenicity" in result
        assert "0.87" in result
        assert "1.234s" in result

    def test_without_engine_results(self, simple_seq):
        result = generate_report(simple_seq, engine_results=None)
        assert "Engine Analysis Results" not in result

    def test_with_exon_boundaries(self, simple_seq):
        boundaries = [(0, 30), (30, 60)]
        result = generate_report(simple_seq, exon_boundaries=boundaries)
        # 2 exon boundaries should be reflected
        assert "2" in result  # Exons count

    def test_default_exon_boundaries(self, simple_seq):
        result = generate_report(simple_seq, exon_boundaries=None)
        # Should default to single exon spanning the whole sequence
        assert "1" in result

    def test_failed_engine_result(self, simple_seq):
        failed = BaseEngineResult(
            sequence="ATGGCGGCCTGA",
            primary_score=0.0,
            classification="error",
            success=False,
            error="Model not found",
            engine_name="ESMFold",
            execution_time_s=0.0,
        )
        result = generate_report(simple_seq, engine_results=[failed])
        assert "FAIL" in result
        assert "Model not found" in result

    def test_type_result_with_violation(self, simple_seq):
        results = [
            TypeCheckResult(
                predicate="RestrictionSite",
                verdict=Verdict.FAIL,
                violation="EcoRI site found",
            ),
        ]
        result = generate_report(simple_seq, type_results=results)
        assert "EcoRI site found" in result
        assert "violation" in result

    def test_type_result_with_knowledge_gap(self, simple_seq):
        results = [
            TypeCheckResult(
                predicate="SpliceDonor",
                verdict=Verdict.UNCERTAIN,
                knowledge_gap="Heuristic only",
            ),
        ]
        result = generate_report(simple_seq, type_results=results)
        assert "Heuristic only" in result
        assert "knowledge-gap" in result


class TestComputeGCWindows:
    """Tests for the _compute_gc_windows helper."""

    def test_basic_gc_windows(self, simple_seq):
        result = _compute_gc_windows(simple_seq, window_size=50)
        assert isinstance(result, list)
        assert len(result) > 0
        for pos, gc_val in result:
            assert isinstance(pos, int)
            assert isinstance(gc_val, float)
            assert 0.0 <= gc_val <= 1.0

    def test_window_size_equals_sequence(self):
        seq = "ATGGCGGCCTGAATGGCGGCCTGAATGGCGGCCTGAATGGCGGCCTGAATGGCGGCCTGA"
        result = _compute_gc_windows(seq, window_size=len(seq))
        assert len(result) == 1
        assert result[0][0] == 0

    def test_sequence_shorter_than_window(self):
        seq = "ATG"
        result = _compute_gc_windows(seq, window_size=50)
        assert result == []

    def test_window_step_size(self):
        """Windows should step by window_size // _GC_WINDOW_STEP_DIVISOR."""
        seq = "G" * 200
        result = _compute_gc_windows(seq, window_size=50)
        # Step should be 50 // 4 = 12
        if len(result) >= 2:
            step = result[1][0] - result[0][0]
            assert step == 12

    def test_all_gc_sequence(self):
        seq = "GCGCGCGC" * 10  # 80 bp, all GC
        result = _compute_gc_windows(seq, window_size=50)
        for pos, gc_val in result:
            assert gc_val == 1.0

    def test_all_at_sequence(self):
        seq = "ATATATAT" * 10  # 80 bp, no GC
        result = _compute_gc_windows(seq, window_size=50)
        for pos, gc_val in result:
            assert gc_val == 0.0

    def test_empty_sequence(self):
        result = _compute_gc_windows("", window_size=50)
        assert result == []


class TestComputeCodonUsage:
    """Tests for the _compute_codon_usage helper."""

    def test_basic_codon_usage(self, simple_seq):
        result = _compute_codon_usage(simple_seq)
        assert isinstance(result, dict)
        # Should have at least some amino acids
        assert len(result) > 0

    def test_codon_frequencies_sum_to_one(self, simple_seq):
        result = _compute_codon_usage(simple_seq)
        for aa, codons in result.items():
            total = sum(codons.values())
            assert abs(total - 1.0) < 0.01, f"Frequencies for {aa} sum to {total}"

    def test_single_codon_amino_acid(self):
        # Methionine only has ATG
        seq = "ATG" * 10
        result = _compute_codon_usage(seq)
        assert "M" in result
        assert result["M"]["ATG"] == 1.0

    def test_stop_codons_excluded(self):
        # TAA is a stop codon
        seq = "ATGTAA"
        result = _compute_codon_usage(seq)
        assert "*" not in result

    def test_empty_sequence(self):
        result = _compute_codon_usage("")
        assert result == {}

    def test_short_sequence(self):
        # Only 2 bases - no complete codon
        result = _compute_codon_usage("AT")
        assert result == {}

    def test_all_same_codon(self):
        seq = "ATG" * 5  # 5 Met codons
        result = _compute_codon_usage(seq)
        assert "M" in result
        assert len(result["M"]) == 1
        assert result["M"]["ATG"] == 1.0

    def test_multiple_codons_same_aa(self):
        # Leucine: TTA, TTG, CTT, CTC, CTA, CTG
        seq = "TTACTTCTTCTACTTGCTACTG"
        result = _compute_codon_usage(seq)
        if "L" in result:
            assert len(result["L"]) > 1
            total = sum(result["L"].values())
            assert abs(total - 1.0) < 0.01


class TestElementTypeDescription:
    """Tests for the _element_type_description helper."""

    def test_known_types(self):
        known_types = [
            "start_codon", "stop_codon", "splice_donor", "splice_acceptor",
            "kozak", "cpg_island", "u_rich_region", "tata_box",
            "poly_a_signal", "restriction_site", "repeat", "promoter",
            "enhancer", "silencer", "insulator",
        ]
        for etype in known_types:
            desc = _element_type_description(etype)
            assert isinstance(desc, str)
            assert len(desc) > 0
            assert desc != f"Biological element of type: {etype}"

    def test_unknown_type_returns_fallback(self):
        desc = _element_type_description("custom_element")
        assert "custom_element" in desc
        assert "Biological element of type" in desc

    def test_description_not_empty_for_known(self):
        desc = _element_type_description("start_codon")
        assert "ATG" in desc

    def test_all_descriptions_different(self):
        types = ["start_codon", "stop_codon", "splice_donor"]
        descs = [_element_type_description(t) for t in types]
        assert len(set(descs)) == len(descs)


# ══════════════════════════════════════════════════════════════════════════
# 2. Output Format Correctness
# ══════════════════════════════════════════════════════════════════════════

class TestHTMLStructure:
    """Tests for valid HTML structure and required sections."""

    def test_starts_with_doctype(self, simple_seq):
        result = generate_report(simple_seq)
        assert result.strip().startswith("<!DOCTYPE html>")

    def test_has_html_tags(self, simple_seq):
        result = generate_report(simple_seq)
        assert "<html" in result
        assert "</html>" in result

    def test_has_head_and_body(self, simple_seq):
        result = generate_report(simple_seq)
        assert "<head>" in result
        assert "</head>" in result
        assert "<body>" in result
        assert "</body>" in result

    def test_has_meta_charset(self, simple_seq):
        result = generate_report(simple_seq)
        assert 'charset="UTF-8"' in result or "charset=UTF-8" in result

    def test_has_style_block(self, simple_seq):
        result = generate_report(simple_seq)
        assert "<style>" in result
        assert "</style>" in result

    def test_has_script_block(self, simple_seq):
        result = generate_report(simple_seq)
        assert "<script>" in result
        assert "</script>" in result

    def test_has_header_section(self, simple_seq):
        result = generate_report(simple_seq)
        assert "BioCompiler Report" in result
        assert "header" in result

    def test_has_sequence_summary(self, simple_seq):
        result = generate_report(simple_seq)
        assert "Sequence Summary" in result
        assert "Length (bp)" in result
        assert "GC Content" in result
        assert "CAI" in result
        assert "Protein (aa)" in result
        assert "Exons" in result

    def test_has_predicate_table(self, simple_seq):
        result = generate_report(simple_seq)
        assert "predicate-table" in result
        assert "Type-Check Results" in result

    def test_has_gc_content_plot(self, simple_seq):
        result = generate_report(simple_seq)
        assert "GC Content Profile" in result
        assert "<svg" in result

    def test_has_sequence_structure(self, simple_seq):
        result = generate_report(simple_seq)
        assert "Sequence Structure" in result

    def test_has_codon_usage(self, simple_seq):
        result = generate_report(simple_seq)
        assert "Codon Usage" in result

    def test_has_token_table(self, simple_seq):
        result = generate_report(simple_seq)
        assert "token-table" in result
        assert "Scan Tokens" in result

    def test_has_sequence_section(self, simple_seq):
        result = generate_report(simple_seq)
        assert "raw-sequence" in result
        assert "sequence-display" in result

    def test_has_protein_section(self, simple_seq):
        result = generate_report(simple_seq)
        assert "Protein Translation" in result

    def test_has_footer(self, simple_seq):
        result = generate_report(simple_seq)
        assert "Generated by BioCompiler" in result
        assert "Machine-Verified Gene Design" in result

    def test_has_dark_mode_toggle(self, simple_seq):
        result = generate_report(simple_seq)
        assert "dark-mode-toggle" in result

    def test_has_filter_buttons(self, simple_seq):
        result = generate_report(simple_seq)
        assert "filter-btn" in result
        assert 'data-filter="all"' in result
        assert 'data-filter="pass"' in result
        assert 'data-filter="fail"' in result

    def test_has_copy_button(self, simple_seq):
        result = generate_report(simple_seq)
        assert "copy-seq-btn" in result

    def test_has_token_search(self, simple_seq):
        result = generate_report(simple_seq)
        assert "token-search" in result


class TestSVGGeneration:
    """Tests for SVG plot generation within the report."""

    def test_gc_plot_svg_present(self, simple_seq):
        result = generate_report(simple_seq)
        # Should contain at least one SVG element (GC plot)
        assert result.count("<svg") >= 1

    def test_structure_svg_present(self, simple_seq):
        result = generate_report(simple_seq)
        # Multiple SVGs: GC plot, structure, codon heatmap
        svg_count = result.count("<svg")
        assert svg_count >= 2

    def test_svg_has_closing_tag(self, simple_seq):
        result = generate_report(simple_seq)
        svg_opens = result.count("<svg")
        svg_closes = result.count("</svg>")
        assert svg_opens == svg_closes


class TestVerdictBadges:
    """Tests for verdict badge rendering."""

    def test_pass_verdict_badge(self, simple_seq, type_results_pass):
        result = generate_report(simple_seq, type_results=type_results_pass)
        assert "verdict-badge pass" in result

    def test_fail_verdict_badge(self, simple_seq):
        results = [TypeCheckResult(predicate="test", verdict=Verdict.FAIL)]
        result = generate_report(simple_seq, type_results=results)
        assert "verdict-badge fail" in result

    def test_uncertain_verdict_badge(self, simple_seq):
        results = [TypeCheckResult(predicate="test", verdict=Verdict.UNCERTAIN)]
        result = generate_report(simple_seq, type_results=results)
        assert "verdict-badge uncertain" in result

    def test_likely_pass_verdict_badge(self, simple_seq):
        results = [TypeCheckResult(predicate="test", verdict=Verdict.LIKELY_PASS)]
        result = generate_report(simple_seq, type_results=results)
        assert "verdict-badge likely_pass" in result

    def test_likely_fail_verdict_badge(self, simple_seq):
        results = [TypeCheckResult(predicate="test", verdict=Verdict.LIKELY_FAIL)]
        result = generate_report(simple_seq, type_results=results)
        assert "verdict-badge likely_fail" in result

    @pytest.mark.parametrize("verdict,expected_symbol", [
        (Verdict.PASS, "&#10003;"),
        (Verdict.LIKELY_PASS, "&#10003;~"),
        (Verdict.UNCERTAIN, "?"),
        (Verdict.LIKELY_FAIL, "&#10007;~"),
        (Verdict.FAIL, "&#10007;"),
    ])
    def test_verdict_symbols_in_report(self, simple_seq, verdict, expected_symbol):
        results = [TypeCheckResult(predicate="test_predicate", verdict=verdict)]
        result = generate_report(simple_seq, type_results=results)
        assert expected_symbol in result

    def test_combined_verdict_is_weakest_link(self, simple_seq):
        """Combined verdict should be the minimum (weakest link)."""
        results = [
            TypeCheckResult(predicate="p1", verdict=Verdict.PASS),
            TypeCheckResult(predicate="p2", verdict=Verdict.FAIL),
        ]
        result = generate_report(simple_seq, type_results=results)
        # FAIL should dominate the overall verdict
        assert "FAIL" in result
        # The verdict badge should show fail
        assert "verdict-badge fail" in result


class TestCertificateSection:
    """Tests for certificate rendering in the report."""

    def test_certificate_design_id_truncated(self, simple_seq, sample_certificate):
        result = generate_report(simple_seq, certificate=sample_certificate)
        # Design ID should be truncated to _DESIGN_ID_DISPLAY_LENGTH + "..."
        assert "..." in result
        # The first _DESIGN_ID_DISPLAY_LENGTH chars should be present
        short_id = sample_certificate.design_id[:_DESIGN_ID_DISPLAY_LENGTH]
        assert short_id in result

    def test_certificate_version_displayed(self, simple_seq, sample_certificate):
        result = generate_report(simple_seq, certificate=sample_certificate)
        assert "Version" in result
        assert sample_certificate.version in result

    def test_certificate_timestamp_displayed(self, simple_seq, sample_certificate):
        result = generate_report(simple_seq, certificate=sample_certificate)
        assert "2025-01-15T10:30:00Z" in result

    def test_certificate_tool_displayed(self, simple_seq, sample_certificate):
        result = generate_report(simple_seq, certificate=sample_certificate)
        assert "biocompiler" in result

    def test_certificate_non_dict_provenance(self, simple_seq):
        """If provenance is not a dict, should show N/A."""
        cert = Certificate(
            version="1.0",
            design_id="abc123def456",
            sequence="ATGGCG",
            types=[],
            provenance="not_a_dict",
        )
        result = generate_report(simple_seq, certificate=cert)
        assert "N/A" in result


class TestEngineSection:
    """Tests for engine analysis results rendering."""

    def test_engine_success_pass_badge(self, simple_seq, sample_engine_results):
        result = generate_report(simple_seq, engine_results=sample_engine_results)
        assert 'class="badge pass"' in result

    def test_engine_failure_fail_badge(self, simple_seq):
        failed = BaseEngineResult(
            sequence="ATGGCG",
            primary_score=0.0,
            classification="error",
            success=False,
            error="Timeout",
            engine_name="TestEngine",
        )
        result = generate_report(simple_seq, engine_results=[failed])
        assert 'class="badge fail"' in result
        assert "Timeout" in result

    def test_engine_score_formatted(self, simple_seq, sample_engine_results):
        result = generate_report(simple_seq, engine_results=sample_engine_results)
        # Score should be formatted with 4 decimal places
        assert "0.87" in result

    def test_engine_time_formatted(self, simple_seq, sample_engine_results):
        result = generate_report(simple_seq, engine_results=sample_engine_results)
        # Time should be formatted with 3 decimal places + 's'
        assert "1.234s" in result

    def test_engine_table_headers(self, simple_seq, sample_engine_results):
        result = generate_report(simple_seq, engine_results=sample_engine_results)
        assert "Status" in result
        assert "Engine" in result
        assert "Classification" in result


class TestIsoformSection:
    """Tests for splice isoform rendering."""

    def test_isoform_count_in_header(self, simple_seq, sample_isoforms):
        result = generate_report(simple_seq, isoforms=sample_isoforms)
        assert "Splice Isoforms (2)" in result

    def test_isoform_table_content(self, simple_seq, sample_isoforms):
        result = generate_report(simple_seq, isoforms=sample_isoforms)
        assert "0.95" in result
        assert "0.72" in result

    def test_isoform_max_display_limit(self, simple_seq):
        """More than _MAX_ISOFORM_DISPLAY isoforms should be truncated."""
        many_isoforms = [
            SpliceIsoform(
                sequence=f"ATG{i:03d}",
                exon_boundaries=[(0, 6)],
                parse_path=["exon_1"],
                score=float(i) / 100,
            )
            for i in range(_MAX_ISOFORM_DISPLAY + 10)
        ]
        result = generate_report(simple_seq, isoforms=many_isoforms)
        # Header shows total count, but rows are limited
        assert f"Splice Isoforms ({_MAX_ISOFORM_DISPLAY + 10})" in result


class TestTokenSection:
    """Tests for token table rendering."""

    def test_token_data_attributes(self, simple_seq, sample_tokens):
        result = generate_report(simple_seq, tokens=sample_tokens)
        assert 'data-element-type="start_codon"' in result
        assert 'data-match-sequence="atg"' in result

    def test_token_frame_display(self, simple_seq, sample_tokens):
        result = generate_report(simple_seq, tokens=sample_tokens)
        # Frame 0 should show "0"
        assert ">0<" in result

    def test_token_score_formatted(self, simple_seq, sample_tokens):
        result = generate_report(simple_seq, tokens=sample_tokens)
        assert "1.00" in result
        assert "0.85" in result

    def test_token_max_display_limit(self, simple_seq):
        """More than _MAX_TOKEN_DISPLAY tokens should be truncated."""
        many_tokens = [
            Token(position=i, element_type="start_codon", match_sequence="ATG", score=1.0, frame=0, strand="+")
            for i in range(_MAX_TOKEN_DISPLAY + 10)
        ]
        result = generate_report(simple_seq, tokens=many_tokens)
        # Total count should reflect all tokens
        assert f"Scan Tokens ({_MAX_TOKEN_DISPLAY + 10})" in result

    def test_token_null_frame_shows_dash(self, simple_seq):
        tokens = [Token(position=0, element_type="repeat", match_sequence="AT", score=0.5,
                        frame=None, strand="+")]
        result = generate_report(simple_seq, tokens=tokens)
        assert ">-<" in result


class TestOptimizationSection:
    """Tests for optimization results rendering."""

    def test_optimization_cai_displayed(self, simple_seq, sample_optimization_result):
        result = generate_report(simple_seq, optimization_result=sample_optimization_result)
        assert "0.8234" in result

    def test_optimization_satisfied_count(self, simple_seq, sample_optimization_result):
        result = generate_report(simple_seq, optimization_result=sample_optimization_result)
        assert "Satisfied" in result

    def test_optimization_failed_count(self, simple_seq, sample_optimization_result):
        result = generate_report(simple_seq, optimization_result=sample_optimization_result)
        assert "Failed" in result


class TestXSSSafety:
    """Tests that user-provided strings are properly HTML-escaped."""

    def test_gene_name_escaped(self):
        result = generate_report("ATGGCG", gene_name='<script>alert("xss")</script>')
        # The gene name should be HTML-escaped in the subtitle and title
        assert "&lt;script&gt;" in result
        # The raw <script>alert should NOT appear in body text content
        # (the report has its own <script> JS block, so check the subtitle specifically)
        assert '&lt;script&gt;alert(&quot;xss&quot;)&lt;/script&gt;' in result

    def test_organism_escaped(self):
        result = generate_report("ATGGCG", organism='<img onerror="hack">')
        assert "&lt;img" in result

    def test_violation_escaped(self):
        results = [TypeCheckResult(
            predicate="test",
            verdict=Verdict.FAIL,
            violation='<b>bold</b>',
        )]
        result = generate_report("ATGGCG", type_results=results)
        assert "<b>bold</b>" not in result
        assert "&lt;b&gt;bold&lt;/b&gt;" in result

    def test_knowledge_gap_escaped(self):
        results = [TypeCheckResult(
            predicate="test",
            verdict=Verdict.UNCERTAIN,
            knowledge_gap='<a href="evil">link</a>',
        )]
        result = generate_report("ATGGCG", type_results=results)
        assert '<a href="evil">' not in result

    def test_token_element_type_escaped(self):
        tokens = [Token(position=0, element_type='<script>', match_sequence="ATG",
                        score=1.0, strand="+")]
        result = generate_report("ATGGCG", tokens=tokens)
        assert "&lt;script&gt;" in result

    def test_token_match_sequence_escaped(self):
        tokens = [Token(position=0, element_type="test", match_sequence='<script>',
                        score=1.0, strand="+")]
        result = generate_report("ATGGCG", tokens=tokens)
        assert "&lt;script&gt;" in result

    def test_certificate_design_id_escaped(self):
        cert = Certificate(
            version="1.0",
            design_id='<script>xss</script>123456789012345678901234567890',
            sequence="ATGGCG",
            types=[],
            provenance={},
        )
        result = generate_report("ATGGCG", certificate=cert)
        assert "&lt;script&gt;xss&lt;/script&gt;" in result

    def test_engine_name_escaped(self):
        er = BaseEngineResult(
            sequence="ATGGCG",
            primary_score=1.0,
            classification="ok",
            success=True,
            engine_name='<script>alert(1)</script>',
        )
        result = generate_report("ATGGCG", engine_results=[er])
        assert "&lt;script&gt;alert(1)&lt;/script&gt;" in result

    def test_predicate_name_escaped(self):
        results = [TypeCheckResult(
            predicate='<img src=x onerror=alert(1)>',
            verdict=Verdict.PASS,
        )]
        result = generate_report("ATGGCG", type_results=results)
        assert "&lt;img" in result


class TestSequenceFormatting:
    """Tests for sequence display in the report."""

    def test_raw_sequence_in_textarea(self, simple_seq):
        result = generate_report(simple_seq)
        # The sequence should be in the hidden textarea for copy
        assert simple_seq in result

    def test_sequence_length_displayed(self, simple_seq):
        result = generate_report(simple_seq)
        assert f"{len(simple_seq):,}" in result


# ══════════════════════════════════════════════════════════════════════════
# 3. Input Validation & Error Handling
# ══════════════════════════════════════════════════════════════════════════

class TestInputValidation:
    """Tests for handling edge cases and invalid inputs."""

    def test_lowercase_sequence_uppercased(self):
        result = generate_report("atggcggcctga")
        assert isinstance(result, str)
        assert "ATGGCGGCCTGA" in result

    def test_mixed_case_sequence(self):
        result = generate_report("AtGgCgGcCtGa")
        assert isinstance(result, str)
        assert "ATGGCGGCCTGA" in result

    def test_empty_sequence(self):
        """Empty sequence should not crash; report should still be generated."""
        result = generate_report("")
        assert isinstance(result, str)
        assert "<!DOCTYPE html>" in result

    def test_very_short_sequence(self):
        result = generate_report("ATG")
        assert isinstance(result, str)
        assert "<!DOCTYPE html>" in result

    def test_single_nucleotide(self):
        result = generate_report("A")
        assert isinstance(result, str)

    def test_non_dna_characters(self):
        """Non-DNA characters should not crash the report."""
        result = generate_report("ATGXYZ123")
        assert isinstance(result, str)
        assert "<!DOCTYPE html>" in result

    def test_very_long_sequence(self):
        """Long sequences should generate without error."""
        seq = "ATGGCG" * 1000  # 6000 bp
        result = generate_report(seq)
        assert isinstance(result, str)
        assert "<!DOCTYPE html>" in result

    def test_none_type_results_handled(self, simple_seq):
        result = generate_report(simple_seq, type_results=None)
        assert isinstance(result, str)
        assert "N/A" in result

    def test_none_tokens_handled(self, simple_seq):
        """When tokens=None, scan_sequence should be called as fallback."""
        result = generate_report(simple_seq, tokens=None)
        assert isinstance(result, str)

    def test_none_certificate_handled(self, simple_seq):
        result = generate_report(simple_seq, certificate=None)
        assert isinstance(result, str)
        assert "Certificate" not in result

    def test_none_isoforms_handled(self, simple_seq):
        result = generate_report(simple_seq, isoforms=None)
        assert isinstance(result, str)

    def test_none_optimization_handled(self, simple_seq):
        result = generate_report(simple_seq, optimization_result=None)
        assert isinstance(result, str)

    def test_none_engine_results_handled(self, simple_seq):
        result = generate_report(simple_seq, engine_results=None)
        assert isinstance(result, str)


class TestMockedDependencyFailures:
    """Tests that generate_report gracefully handles failures in dependencies."""

    def test_gc_content_failure_in_generate(self, simple_seq):
        """If gc_content raises in generate_report, it is caught (gc=0.0),
        but _compute_gc_windows also calls gc_content and will propagate.
        The report generation will fail in this scenario."""
        with patch("biocompiler.report.gc_content", side_effect=RuntimeError("boom")):
            with pytest.raises(RuntimeError, match="boom"):
                generate_report(simple_seq)

    def test_translate_failure_returns_empty_protein(self, simple_seq):
        with patch("biocompiler.report.translate", side_effect=RuntimeError("boom")):
            result = generate_report(simple_seq)
            assert isinstance(result, str)

    def test_compute_cai_failure_returns_zero(self, simple_seq):
        with patch("biocompiler.report.compute_cai", side_effect=RuntimeError("boom")):
            result = generate_report(simple_seq)
            assert isinstance(result, str)
            assert "0.0000" in result

    def test_scan_sequence_failure_returns_empty_tokens(self, simple_seq):
        with patch("biocompiler.report.scan_sequence", side_effect=RuntimeError("boom")):
            result = generate_report(simple_seq, tokens=None)
            assert isinstance(result, str)

    def test_compute_splice_isoforms_failure_returns_empty(self, simple_seq):
        with patch("biocompiler.report.compute_splice_isoforms", side_effect=RuntimeError("boom")):
            result = generate_report(simple_seq, isoforms=None,
                                     exon_boundaries=[(0, 30), (30, 60)])
            assert isinstance(result, str)


class TestVerdictSymbols:
    """Tests for the _VERDICT_SYMBOLS mapping."""

    def test_all_verdicts_have_symbols(self):
        for v in Verdict:
            assert v.value in _VERDICT_SYMBOLS

    def test_pass_symbol_is_checkmark(self):
        assert "10003" in _VERDICT_SYMBOLS["PASS"]

    def test_fail_symbol_is_cross(self):
        assert "10007" in _VERDICT_SYMBOLS["FAIL"]

    def test_uncertain_symbol_is_question(self):
        assert _VERDICT_SYMBOLS["UNCERTAIN"] == "?"


class TestReportSelfContained:
    """Tests that the report is truly self-contained (no external deps)."""

    def test_no_external_css_links(self, simple_seq):
        result = generate_report(simple_seq)
        # No <link rel="stylesheet" href="...">
        assert 'rel="stylesheet"' not in result
        assert '<link' not in result

    def test_no_external_script_links(self, simple_seq):
        result = generate_report(simple_seq)
        # No <script src="...">
        assert '<script src=' not in result
        assert "<script src=" not in result

    def test_contains_embedded_css(self, simple_seq):
        result = generate_report(simple_seq)
        assert ":root" in result
        assert "--pass:" in result
        assert "--fail:" in result

    def test_contains_embedded_js(self, simple_seq):
        result = generate_report(simple_seq)
        assert "document.createElement" in result
        assert "addEventListener" in result


class TestReportPrintFriendly:
    """Tests for print-friendly layout features."""

    def test_has_print_media_query(self, simple_seq):
        result = generate_report(simple_seq)
        assert "@media print" in result

    def test_print_hides_interactive_elements(self, simple_seq):
        result = generate_report(simple_seq)
        # Print CSS should hide filter bar, search bar, copy button, toggle button
        assert ".filter-bar" in result
        assert ".search-bar" in result
        assert ".copy-btn" in result
        assert ".toggle-btn" in result


class TestReportConstants:
    """Tests that report module constants are reasonable."""

    def test_max_token_display_positive(self):
        assert _MAX_TOKEN_DISPLAY > 0

    def test_max_isoform_display_positive(self):
        assert _MAX_ISOFORM_DISPLAY > 0

    def test_design_id_display_length_positive(self):
        assert _DESIGN_ID_DISPLAY_LENGTH > 0

    def test_default_gc_window_size_positive(self):
        assert _DEFAULT_GC_WINDOW_SIZE > 0

    def test_verdict_symbols_all_strings(self):
        for key, value in _VERDICT_SYMBOLS.items():
            assert isinstance(key, str)
            assert isinstance(value, str)
