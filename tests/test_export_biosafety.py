"""
Comprehensive tests for biosafety annotations in biocompiler.export.

Covers:
- GenBank output contains BIOCOMPILER_ANNOTATIONS section
- GenBank output contains BIOSECURITY NOTICE
- FASTA header has biosecurity metadata
- FASTA type_results parameter for BSL escalation
- WARNING is included when predicates fail
- Biosecurity notice formatting varies by screening result
- BSL-1 for prokaryotes, BSL-2 for mammalian/failed predicates
- export_with_annotations() wraps GenBank and FASTA exports
- export_with_annotations() passes type_results to FASTA
- format_biosecurity_report() produces valid report
- _assess_biosafety_level() classification logic
- _is_biosecurity_screened() detection
- JSON export includes biosafety annotations
- Provenance ID format and uniqueness
"""

from __future__ import annotations

import pytest

from biocompiler.export import (
    _assess_biosafety_level,
    _is_biosecurity_screened,
    export_fasta,
    export_genbank,
    export_with_annotations,
    format_biosecurity_report,
    export_json,
)
from biocompiler.types import TypeCheckResult, Verdict


# ─── Test fixtures ────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"

PASS_RESULTS = [
    TypeCheckResult(predicate="gc_content", verdict=Verdict.PASS),
    TypeCheckResult(predicate="no_stop_codons", verdict=Verdict.PASS),
]

FAIL_RESULTS = [
    TypeCheckResult(predicate="gc_content", verdict=Verdict.PASS),
    TypeCheckResult(predicate="no_stop_codons", verdict=Verdict.FAIL, violation="Stop codon found"),
]

MIXED_RESULTS = [
    TypeCheckResult(predicate="gc_content", verdict=Verdict.PASS),
    TypeCheckResult(predicate="no_stop_codons", verdict=Verdict.PASS),
    TypeCheckResult(predicate="cpg_island", verdict=Verdict.UNCERTAIN),
]


# ═══════════════════════════════════════════════════════════════════════
# 1. _assess_biosafety_level
# ═══════════════════════════════════════════════════════════════════════

class TestAssessBiosafetyLevel:

    def test_prokaryote_bsl1(self):
        assert _assess_biosafety_level("Escherichia_coli") == "BSL-1"

    def test_e_coli_alias_bsl1(self):
        assert _assess_biosafety_level("E_coli") == "BSL-1"

    def test_bacillus_bsl1(self):
        assert _assess_biosafety_level("Bacillus_subtilis") == "BSL-1"

    def test_human_bsl2(self):
        assert _assess_biosafety_level("Homo_sapiens") == "BSL-2"

    def test_mouse_bsl2(self):
        assert _assess_biosafety_level("Mus_musculus") == "BSL-2"

    def test_cho_bsl2(self):
        assert _assess_biosafety_level("CHO_K1") == "BSL-2"

    def test_yeast_unknown(self):
        """Yeast is not in our explicit BSL-1/BSL-2 sets."""
        assert _assess_biosafety_level("Saccharomyces_cerevisiae") == "unknown"

    def test_prokaryote_escalated_with_fail(self):
        """Prokaryote with failed predicates should escalate to BSL-2."""
        assert _assess_biosafety_level("Escherichia_coli", FAIL_RESULTS) == "BSL-2"

    def test_prokaryote_stays_bsl1_with_pass(self):
        assert _assess_biosafety_level("Escherichia_coli", PASS_RESULTS) == "BSL-1"

    def test_human_stays_bsl2_with_pass(self):
        assert _assess_biosafety_level("Homo_sapiens", PASS_RESULTS) == "BSL-2"

    def test_human_stays_bsl2_with_fail(self):
        assert _assess_biosafety_level("Homo_sapiens", FAIL_RESULTS) == "BSL-2"

    def test_no_type_results(self):
        assert _assess_biosafety_level("Escherichia_coli", None) == "BSL-1"


# ═══════════════════════════════════════════════════════════════════════
# 2. _is_biosecurity_screened
# ═══════════════════════════════════════════════════════════════════════

class TestIsBiosecurityScreened:

    def test_all_pass(self):
        assert _is_biosecurity_screened(PASS_RESULTS) is True

    def test_has_fail(self):
        assert _is_biosecurity_screened(FAIL_RESULTS) is False

    def test_has_uncertain_only(self):
        """UNCERTAIN verdicts should not count as failed."""
        assert _is_biosecurity_screened(MIXED_RESULTS) is True

    def test_none_results(self):
        assert _is_biosecurity_screened(None) is False

    def test_empty_results(self):
        # Empty list means no type-check was performed, so not screened
        assert _is_biosecurity_screened([]) is False


# ═══════════════════════════════════════════════════════════════════════
# 3. GenBank biosafety annotations
# ═══════════════════════════════════════════════════════════════════════

class TestGenBankBiosafetyAnnotations:

    def test_has_biocompiler_annotations_section(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "BIOCOMPILER_ANNOTATIONS:" in result

    def test_has_optimized_by(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "optimized_by: biocompiler v" in result

    def test_has_organism(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "organism: Escherichia_coli" in result

    def test_has_cai_score(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli", cai=0.95)
        assert "cai_score: 0.9500" in result

    def test_has_cai_score_when_auto_computed(self):
        # Short sequences may auto-compute CAI, so check for cai_score line
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli", cai=None)
        # Either cai_score has a value or N/A
        assert "cai_score:" in result

    def test_has_gc_content(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "gc_content:" in result

    def test_has_passed_predicates(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        assert "passed_predicates:" in result
        assert "gc_content" in result

    def test_has_failed_predicates(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "failed_predicates:" in result
        assert "no_stop_codons" in result

    def test_has_biosecurity_screened_true(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        assert "biosecurity_screened: True" in result

    def test_has_biosecurity_screened_false_with_fail(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "biosecurity_screened: False" in result

    def test_has_biosafety_level_bsl1(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "biosafety_level: BSL-1" in result

    def test_has_biosafety_level_bsl2(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "biosafety_level: BSL-2" in result

    def test_has_provenance_id(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "provenance_id: BC_" in result

    def test_has_biosecurity_notice(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "BIOSECURITY NOTICE:" in result

    def test_notice_safe_when_passed(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        assert "safe for synthesis" in result

    def test_notice_warning_when_failed(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "FAILED one or more biosecurity predicates" in result

    def test_notice_bsl1_risk(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        assert "BSL-1" in result
        assert "Standard laboratory practices" in result

    def test_notice_bsl2_risk(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "BSL-2" in result
        assert "BSL-2 containment procedures" in result


# ═══════════════════════════════════════════════════════════════════════
# 4. GenBank WARNING when predicates fail
# ═══════════════════════════════════════════════════════════════════════

class TestGenBankWarningOnFailure:

    def test_warning_present_when_predicates_fail(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "WARNING:" in result
        assert "predicate(s) failed" in result

    def test_warning_not_present_when_all_pass(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        # "WARNING:" should not appear in the output when all pass
        assert "WARNING: 1 predicate(s) failed" not in result

    def test_warning_includes_predicate_names(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "no_stop_codons" in result

    def test_warning_includes_review_advice(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "Review before submitting for gene synthesis" in result


# ═══════════════════════════════════════════════════════════════════════
# 5. FASTA header biosafety metadata
# ═══════════════════════════════════════════════════════════════════════

class TestFastaBiosafetyMetadata:

    def test_header_has_biosecurity_level(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "biosecurity=BSL-1" in result

    def test_header_has_biosecurity_bsl2(self):
        result = export_fasta(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "biosecurity=BSL-2" in result

    def test_header_has_biocompiler_version(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "biocompiler_v" in result

    def test_header_has_organism(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "organism=Escherichia_coli" in result

    def test_header_has_cai(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli", cai=0.95)
        assert "CAI=0.9500" in result

    def test_header_has_gc(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "GC=" in result

    def test_header_format_matches_spec(self):
        """Header should follow: >name |organism={org}|CAI={cai}|GC={gc}|biosecurity={risk}|biocompiler_v11.1"""
        result = export_fasta(SAMPLE_SEQ, identifier="test_gene",
                             organism="Escherichia_coli", cai=0.95)
        # Find the header line (starts with >)
        header_line = None
        for line in result.split("\n"):
            if line.startswith(">"):
                header_line = line
                break
        assert header_line is not None
        header_content = header_line.lstrip(">")
        parts = header_content.split("|")
        assert parts[0] == "test_gene"
        assert any(p.startswith("organism=") for p in parts)
        assert any(p.startswith("CAI=") for p in parts)
        assert any(p.startswith("GC=") for p in parts)
        assert any(p.startswith("biosecurity=") for p in parts)
        assert any(p.startswith("biocompiler_v") for p in parts)


# ═══════════════════════════════════════════════════════════════════════
# 6. export_with_annotations
# ═══════════════════════════════════════════════════════════════════════

class TestExportWithAnnotations:

    def test_genbank_has_annotations(self):
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            cai=0.95,
            type_results=PASS_RESULTS,
            format="genbank",
        )
        assert "BIOCOMPILER_ANNOTATIONS:" in result
        assert "biosafety_level: BSL-1" in result

    def test_fasta_has_biosecurity(self):
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            cai=0.95,
            format="fasta",
        )
        assert "biosecurity=BSL-1" in result

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Unsupported format"):
            export_with_annotations(
                SAMPLE_SEQ,
                organism="Escherichia_coli",
                format="xyz",
            )

    def test_genbank_with_failed_predicates(self):
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            type_results=FAIL_RESULTS,
            format="genbank",
        )
        assert "WARNING:" in result
        assert "biosecurity_screened: False" in result

    def test_genbank_passes_kwargs(self):
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            format="genbank",
            gene_name="testGene",
        )
        assert 'gene="testGene"' in result

    def test_fasta_passes_kwargs(self):
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            format="fasta",
            identifier="my_seq",
        )
        assert ">my_seq" in result

    def test_auto_computes_gc(self):
        """GC should be auto-computed when not provided."""
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            gc=None,
            format="genbank",
        )
        assert "gc_content:" in result


# ═══════════════════════════════════════════════════════════════════════
# 7. format_biosecurity_report
# ═══════════════════════════════════════════════════════════════════════

class TestFormatBiosecurityReport:

    def test_report_has_header(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "BIOSECURITY SCREENING REPORT" in report

    def test_report_has_biosafety_level(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "Biosecurity level: BSL-1" in report

    def test_report_has_organism(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "Organism" in report
        assert "Escherichia coli" in report  # display format

    def test_report_has_screened_status(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=PASS_RESULTS)
        assert "Screened" in report

    def test_report_pass_when_screened(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=PASS_RESULTS)
        assert "PASS" in report

    def test_report_fail_when_not_screened(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=FAIL_RESULTS)
        assert "FAIL" in report

    def test_report_has_passed_predicates(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=PASS_RESULTS)
        assert "Passed predicates" in report
        assert "gc_content" in report

    def test_report_has_failed_predicates(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=FAIL_RESULTS)
        assert "Failed predicates" in report
        assert "no_stop_codons" in report

    def test_report_has_risk_assessment(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=PASS_RESULTS)
        assert "Risk Assessment Summary" in report
        assert "safe for synthesis" in report

    def test_report_warning_when_failed(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=FAIL_RESULTS)
        assert "WARNING" in report

    def test_report_bsl1_recommendation(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          type_results=PASS_RESULTS)
        assert "BSL-1" in report
        assert "Standard laboratory practices" in report

    def test_report_bsl2_recommendation(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Homo_sapiens",
                                          type_results=PASS_RESULTS)
        assert "BSL-2" in report
        assert "BSL-2 containment procedures" in report

    def test_report_unknown_organism(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Alien_genome")
        assert "unknown" in report
        assert "biosafety officer" in report

    def test_report_has_provenance_id(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "Provenance ID" in report
        assert "BC_" in report

    def test_report_has_cai_score(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          cai=0.95)
        assert "CAI score" in report
        assert "0.9500" in report

    def test_report_cai_na_when_missing(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli",
                                          cai=None)
        assert "CAI score        : N/A" in report

    def test_report_has_gc_content(self):
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "GC content" in report


# ═══════════════════════════════════════════════════════════════════════
# 8. Biosecurity notice formatting
# ═══════════════════════════════════════════════════════════════════════

class TestBiosecurityNoticeFormatting:

    def test_notice_pass_format(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        # Find the BIOSECURITY NOTICE section
        assert "BIOSECURITY NOTICE:" in result
        assert "passed all biosecurity screening predicates" in result

    def test_notice_fail_format(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=FAIL_RESULTS)
        assert "BIOSECURITY NOTICE:" in result
        assert "FAILED one or more biosecurity predicates" in result
        assert "biosafety officer" in result

    def test_notice_bsl1_risk_level(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        assert "Risk level: BSL-1" in result
        assert "Standard laboratory practices are sufficient" in result

    def test_notice_bsl2_risk_level(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens",
                               type_results=PASS_RESULTS)
        assert "Risk level: BSL-2" in result
        assert "Institutional BSL-2 containment procedures apply" in result

    def test_notice_unknown_risk_level(self):
        result = export_genbank(SAMPLE_SEQ, organism="Saccharomyces_cerevisiae")
        assert "Risk level: Unknown" in result
        assert "Consult institutional biosafety officer" in result

    def test_notice_appears_after_annotations(self):
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli",
                               type_results=PASS_RESULTS)
        annotations_pos = result.index("BIOCOMPILER_ANNOTATIONS:")
        notice_pos = result.index("BIOSECURITY NOTICE:")
        assert notice_pos > annotations_pos


# ═══════════════════════════════════════════════════════════════════════
# 9. Integration: GenBank structural validity with annotations
# ═══════════════════════════════════════════════════════════════════════

class TestGenBankStructuralValidityWithAnnotations:

    def test_has_locus(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert result.startswith("LOCUS")

    def test_has_features(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert "FEATURES" in result

    def test_has_origin(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert "ORIGIN" in result

    def test_has_terminator(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert result.rstrip().endswith("//")

    def test_locus_before_features(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert result.index("LOCUS") < result.index("FEATURES")

    def test_features_before_origin(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert result.index("FEATURES") < result.index("ORIGIN")

    def test_comment_section_contains_annotations(self):
        result = export_genbank(SAMPLE_SEQ, type_results=PASS_RESULTS)
        assert "COMMENT" in result
        assert "BIOCOMPILER_ANNOTATIONS:" in result


# ═══════════════════════════════════════════════════════════════════════
# 10. FASTA type_results parameter for BSL escalation
# ═══════════════════════════════════════════════════════════════════════

class TestFastaTypeResultsEscalation:

    def test_fasta_without_type_results_bsl1(self):
        """Without type_results, E. coli should be BSL-1."""
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "biosecurity=BSL-1" in result

    def test_fasta_with_pass_results_bsl1(self):
        """With all-pass type_results, E. coli should stay BSL-1."""
        result = export_fasta(
            SAMPLE_SEQ, organism="Escherichia_coli",
            type_results=PASS_RESULTS,
        )
        assert "biosecurity=BSL-1" in result

    def test_fasta_with_fail_results_escalates_bsl2(self):
        """Failed predicates should escalate E. coli to BSL-2 in FASTA header."""
        result = export_fasta(
            SAMPLE_SEQ, organism="Escherichia_coli",
            type_results=FAIL_RESULTS,
        )
        assert "biosecurity=BSL-2" in result

    def test_fasta_with_none_type_results_bsl1(self):
        """None type_results should not change organism-based BSL."""
        result = export_fasta(
            SAMPLE_SEQ, organism="Escherichia_coli",
            type_results=None,
        )
        assert "biosecurity=BSL-1" in result

    def test_fasta_human_stays_bsl2_with_fail(self):
        """Homo sapiens stays BSL-2 even with fail results."""
        result = export_fasta(
            SAMPLE_SEQ, organism="Homo_sapiens",
            type_results=FAIL_RESULTS,
        )
        assert "biosecurity=BSL-2" in result

    def test_fasta_header_spec_format(self):
        """FASTA header should match spec: >name|organism={org}|CAI={cai}|GC={gc}|biosecurity={risk}|biocompiler_v{ver}"""
        result = export_fasta(
            SAMPLE_SEQ, identifier="myGene",
            organism="Escherichia_coli", cai=0.95,
            type_results=PASS_RESULTS,
        )
        header_line = None
        for line in result.split("\n"):
            if line.startswith(">"):
                header_line = line
                break
        assert header_line is not None
        content = header_line.lstrip(">")
        parts = content.split("|")
        # Validate pipe-delimited structure
        assert parts[0] == "myGene"
        assert any(p == "organism=Escherichia_coli" for p in parts)
        assert any(p == "CAI=0.9500" for p in parts)
        assert any(p.startswith("GC=") for p in parts)
        assert any(p.startswith("biosecurity=BSL-1") for p in parts)
        assert any(p.startswith("biocompiler_v") for p in parts)


# ═══════════════════════════════════════════════════════════════════════
# 11. export_with_annotations passes type_results to FASTA
# ═══════════════════════════════════════════════════════════════════════

class TestExportWithAnnotationsTypeResultsPassthrough:

    def test_fasta_with_type_results_escalates_bsl(self):
        """export_with_annotations should pass type_results to FASTA for BSL escalation."""
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            type_results=FAIL_RESULTS,
            format="fasta",
        )
        assert "biosecurity=BSL-2" in result

    def test_fasta_with_pass_type_results_stays_bsl1(self):
        """export_with_annotations with PASS results should keep BSL-1 for E. coli."""
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            type_results=PASS_RESULTS,
            format="fasta",
        )
        assert "biosecurity=BSL-1" in result

    def test_fasta_no_type_results_default_bsl(self):
        """export_with_annotations without type_results uses organism-based BSL."""
        result = export_with_annotations(
            SAMPLE_SEQ,
            organism="Escherichia_coli",
            format="fasta",
        )
        assert "biosecurity=BSL-1" in result


# ═══════════════════════════════════════════════════════════════════════
# 12. Provenance ID format and uniqueness
# ═══════════════════════════════════════════════════════════════════════

class TestProvenanceId:

    def test_provenance_id_format(self):
        """Provenance ID should be BC_ followed by 12 hex characters."""
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        # Find provenance_id in the output
        for line in result.split("\n"):
            if "provenance_id: BC_" in line:
                # Extract the ID
                idx = line.index("provenance_id: BC_") + len("provenance_id: ")
                prov_id = line[idx:].strip()
                assert prov_id.startswith("BC_")
                hex_part = prov_id[3:]
                assert len(hex_part) == 12
                assert all(c in "0123456789ABCDEF" for c in hex_part)
                break
        else:
            pytest.fail("provenance_id not found in GenBank output")

    def test_provenance_id_unique_across_exports(self):
        """Each export should get a unique provenance_id."""
        result1 = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        result2 = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        # Extract provenance IDs
        ids = []
        for result in [result1, result2]:
            for line in result.split("\n"):
                if "provenance_id: BC_" in line:
                    idx = line.index("provenance_id: BC_") + len("provenance_id: ")
                    ids.append(line[idx:].strip())
                    break
        assert len(ids) == 2
        assert ids[0] != ids[1]

    def test_report_provenance_id_format(self):
        """format_biosecurity_report provenance ID should be BC_ + 12 hex chars."""
        report = format_biosecurity_report(SAMPLE_SEQ, organism="Escherichia_coli")
        for line in report.split("\n"):
            if "Provenance ID" in line:
                idx = line.index("BC_")
                prov_id = line[idx:].strip()
                assert prov_id.startswith("BC_")
                assert len(prov_id[3:]) == 12
                break


# ═══════════════════════════════════════════════════════════════════════
# 13. JSON export biosafety annotations
# ═══════════════════════════════════════════════════════════════════════

class TestJsonBiosafetyAnnotations:

    @staticmethod
    def _make_result(organism_name="Escherichia_coli", satisfied=None, failed=None):
        """Helper to create an OptimizationResult for JSON biosafety tests."""
        from biocompiler.optimization import OptimizationResult
        # SAMPLE_SEQ is 42 bases, so protein must be 14 aa
        result = OptimizationResult(
            sequence=SAMPLE_SEQ,
            protein="MVEKGGCGAGGAGC",  # 14 aa = 42/3
            cai=0.95,
            gc_content=0.65,
            satisfied_predicates=satisfied or ["gc_content", "no_stop_codons"],
            failed_predicates=failed or [],
            codon_pair_bias=None,
            mutagenesis_applied=False,
            fallback_used=False,
            suggested_5utr=None,
            suggested_3utr=None,
            utr_score_5=None,
            utr_score_3=None,
            aa_substitutions=[],
            mrna_stability_score=None,
            destabilizing_motifs_removed=0,
            stability_improvement=None,
            certificate_text=None,
            provenance=None,
            decision_trail=None,
        )
        # organism_name is set by the optimizer, not the constructor
        result.organism_name = organism_name
        return result

    def test_json_has_biosafety_section(self):
        """JSON export should include a biosafety section."""
        from biocompiler.export import export_json
        result = self._make_result()
        json_str = export_json(result)
        assert '"biosafety"' in json_str

    def test_json_biosafety_has_required_fields(self):
        """JSON biosafety section should contain all annotation fields."""
        import json
        from biocompiler.export import export_json
        result = self._make_result(
            satisfied=["gc_content"],
            failed=["no_stop_codons"],
        )
        json_str = export_json(result)
        data = json.loads(json_str)
        assert "biosafety" in data
        bs = data["biosafety"]
        assert "optimized_by" in bs
        assert "organism" in bs
        assert "cai_score" in bs
        assert "gc_content" in bs
        assert "passed_predicates" in bs
        assert "failed_predicates" in bs
        assert "biosecurity_screened" in bs
        assert "biosafety_level" in bs
        assert "provenance_id" in bs

    def test_json_biosafety_screened_false_with_failures(self):
        """biosecurity_screened should be False when there are failed predicates."""
        import json
        from biocompiler.export import export_json
        result = self._make_result(
            satisfied=["gc_content"],
            failed=["no_stop_codons"],
        )
        json_str = export_json(result)
        data = json.loads(json_str)
        assert data["biosafety"]["biosecurity_screened"] is False
        assert data["biosafety"]["failed_predicates"] == ["no_stop_codons"]

    def test_json_biosafety_screened_true_when_all_pass(self):
        """biosecurity_screened should be True when all predicates pass."""
        import json
        from biocompiler.export import export_json
        result = self._make_result(
            satisfied=["gc_content", "no_stop_codons"],
            failed=[],
        )
        json_str = export_json(result)
        data = json.loads(json_str)
        assert data["biosafety"]["biosecurity_screened"] is True
        assert data["biosafety"]["provenance_id"].startswith("BC_")

    def test_json_biosafety_level_ecoli(self):
        """JSON biosafety_level should be BSL-1 for E. coli."""
        import json
        from biocompiler.export import export_json
        result = self._make_result(organism_name="Escherichia_coli")
        json_str = export_json(result)
        data = json.loads(json_str)
        assert data["biosafety"]["biosafety_level"] == "BSL-1"

    def test_json_biosafety_level_human(self):
        """JSON biosafety_level should be BSL-2 for Homo sapiens."""
        import json
        from biocompiler.export import export_json
        result = self._make_result(organism_name="Homo_sapiens")
        json_str = export_json(result)
        data = json.loads(json_str)
        assert data["biosafety"]["biosafety_level"] == "BSL-2"


# ═══════════════════════════════════════════════════════════════════════
# 14. BIOCOMPILER_ANNOTATIONS field completeness
# ═══════════════════════════════════════════════════════════════════════

class TestBiocompilerAnnotationsFieldCompleteness:

    def test_all_8_annotation_fields_present(self):
        """BIOCOMPILER_ANNOTATIONS should contain all 8 required fields."""
        result = export_genbank(
            SAMPLE_SEQ, organism="Escherichia_coli",
            cai=0.95, type_results=PASS_RESULTS,
        )
        # Extract the BIOCOMPILER_ANNOTATIONS section
        in_section = False
        found_fields = []
        for line in result.split("\n"):
            if "BIOCOMPILER_ANNOTATIONS:" in line:
                in_section = True
                continue
            if in_section:
                stripped = line.strip()
                if stripped.startswith("WARNING") or stripped.startswith("BIOSECURITY"):
                    break
                if ":" in stripped:
                    field_name = stripped.split(":")[0].strip()
                    found_fields.append(field_name)
        required_fields = [
            "optimized_by", "organism", "cai_score", "gc_content",
            "passed_predicates", "failed_predicates",
            "biosecurity_screened", "biosafety_level", "provenance_id",
        ]
        for field in required_fields:
            assert field in found_fields, f"Missing field: {field}"

    def test_annotations_section_ordering(self):
        """BIOCOMPILER_ANNOTATIONS should appear before BIOSECURITY NOTICE."""
        result = export_genbank(SAMPLE_SEQ, organism="Escherichia_coli")
        annotations_pos = result.index("BIOCOMPILER_ANNOTATIONS:")
        notice_pos = result.index("BIOSECURITY NOTICE:")
        assert annotations_pos < notice_pos

    def test_warning_between_annotations_and_notice(self):
        """WARNING should appear between BIOCOMPILER_ANNOTATIONS and BIOSECURITY NOTICE."""
        result = export_genbank(
            SAMPLE_SEQ, organism="Escherichia_coli",
            type_results=FAIL_RESULTS,
        )
        annotations_pos = result.index("BIOCOMPILER_ANNOTATIONS:")
        warning_pos = result.index("WARNING:")
        notice_pos = result.index("BIOSECURITY NOTICE:")
        assert annotations_pos < warning_pos < notice_pos
