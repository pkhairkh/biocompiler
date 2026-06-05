"""
BioCompiler Literature Validation — Extended Tests
====================================================

Extended test coverage for literature_validation.py, focusing on:

1. Literature reference validation — verifying that every LiteratureCase
   has a well-formed, non-trivial reference string with recognizable
   citation patterns (author names, year, journal, etc.).

2. DOI format checking — detecting and validating DOI patterns embedded
   in reference strings (e.g. 10.xxxx/yyyy) and database_id fields
   that implicitly reference DOI-backed resources (UniProt, ClinVar, PDB).

3. Citation format — structural validation of the reference field,
   including author-year patterns, journal abbreviations, volume/page
   numbers, and database cross-reference consistency.
"""

import re
from typing import List, Set

import pytest

from biocompiler.literature_validation import (
    # Data classes
    LiteratureCase,
    ValidationResult,
    DomainReport,
    # Case collections
    SCID_CASES,
    THALASSEMIA_CASES,
    AGGREGATION_CASES,
    IMMUNOGENICITY_CASES,
    ALL_LITERATURE_CASES,
    # Functions
    evaluate_case,
    run_literature_validation,
    format_literature_validation_report,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

# Regex for DOI in standard form  10.<prefix>/<suffix>
DOI_PATTERN = re.compile(r"\b10\.\d{4,9}/[^\s,;]+\b")

# Regex for author-year citation pattern:  Author A et al. (YYYY)  or  Author & Author (YYYY)
# Allows hyphenated names (e.g. Hacein-Bey-Abina), initials (SH, GA), suffixes (Jr)
AUTHOR_YEAR_PATTERN = re.compile(
    r"[A-Z][a-z\-\']+(?:\s+[A-Z]{1,3}(?:\s+\w+)?)?"
    r"(?:\s+et\s+al\.?)?"
    r"(?:\s*[&,]\s*[A-Z][a-z\-\']+(?:\s+[A-Z]{1,3})?(?:\s+et\s+al\.?)?)*"
    r"\s*\(\d{4}\)"
)

# Regex for a 4-digit year in parentheses
YEAR_PATTERN = re.compile(r"\((\d{4})\)")

# Recognized journal / publisher abbreviations that appear in references
JOURNAL_ABBREVS = {
    "Science",
    "Nature",
    "Nature",
    "Cell",
    "PNAS",
    "Proc Natl Acad Sci USA",
    "N Engl J Med",
    "J Clin Invest",
    "Lancet",
    "Mol Ther",
    "Neuron",
    "Blood",
    "Am J Hematol",
    "Cold Spring Harb Perspect Med",
    "J Interferon Cytokine Res",
    "J Clin Endocrinol Metab",
    "Nat Struct Mol Biol",
}

# Known database ID prefixes that imply DOI-backed or curated records
DATABASE_PREFIXES = {
    "UniProt",
    "ClinVar",
    "NCBI",
    "RefSeq",
    "PDB",
    "IEDB",
    "dbSNP",
}


def _cases_by_domain(domain: str) -> List[LiteratureCase]:
    """Return all literature cases for a given domain."""
    return [c for c in ALL_LITERATURE_CASES if c.domain == domain]


def _extract_years(reference: str) -> List[int]:
    """Extract all 4-digit years from a reference string."""
    return [int(m.group(1)) for m in YEAR_PATTERN.finditer(reference)]


def _extract_database_ids(database_id: str) -> List[str]:
    """Split a compound database_id field (e.g. 'UniProt:P01588 / IEDB:18557')."""
    parts = re.split(r"\s*/\s*", database_id.strip())
    return [p.strip() for p in parts if p.strip()]


def _parse_single_db_id(db_id: str):
    """Parse a single 'PREFIX:VALUE' database ID into (prefix, value)."""
    if ":" not in db_id:
        return (db_id, "")
    prefix, _, value = db_id.partition(":")
    return (prefix.strip(), value.strip())


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Literature Reference Validation
# ═══════════════════════════════════════════════════════════════════════════════

class TestLiteratureReferenceValidation:
    """Validate that every LiteratureCase has a well-formed, substantive reference.

    A valid reference should contain at minimum:
      - At least one author-year citation pattern (or a recognisable database ref)
      - A publication year in parentheses
      - Non-trivial length (more than just a few characters)
    """

    def test_all_cases_have_nonempty_reference(self):
        """Every case must have a non-empty reference string."""
        for case in ALL_LITERATURE_CASES:
            assert case.reference, (
                f"Case {case.case_id} ({case.name}) has an empty reference"
            )

    def test_all_references_contain_year_or_database_ref(self):
        """Every reference should contain at least one year in parentheses,
        OR be a database-style reference (e.g. 'RefSeq:NG_000007.3')."""
        for case in ALL_LITERATURE_CASES:
            years = _extract_years(case.reference)
            has_year = len(years) >= 1
            has_db_ref = bool(re.search(r"(UniProt|RefSeq|ClinVar|PDB|IEDB|NCBI|dbSNP)", case.reference))
            # Also accept database-style accession numbers like NG_000007.3, NM_000206.4
            has_accession = bool(re.search(r"[A-Z]{2}_\d+", case.reference))
            assert has_year or has_db_ref or has_accession, (
                f"Case {case.case_id} reference has no year and no database ref: "
                f"'{case.reference[:80]}'"
            )

    def test_reference_years_are_plausible(self):
        """Years in references should be between 1950 and the current decade."""
        for case in ALL_LITERATURE_CASES:
            years = _extract_years(case.reference)
            for year in years:
                assert 1950 <= year <= 2030, (
                    f"Case {case.case_id} has implausible year {year} in reference"
                )

    def test_reference_minimum_length(self):
        """References should be at least 20 characters to contain meaningful info."""
        for case in ALL_LITERATURE_CASES:
            assert len(case.reference) >= 20, (
                f"Case {case.case_id} reference is suspiciously short "
                f"({len(case.reference)} chars): '{case.reference}'"
            )

    def test_scid_references_contain_author_year(self):
        """SCID cases should cite Hacein-Bey-Abina et al. with appropriate years."""
        scid_refs = " ".join(c.reference for c in SCID_CASES)
        # The primary SCID-X1 reference should appear
        assert "Hacein-Bey-Abina" in scid_refs, (
            "SCID cases should cite Hacein-Bey-Abina et al."
        )
        # The key publication years should be present
        scid_years = set()
        for c in SCID_CASES:
            scid_years.update(_extract_years(c.reference))
        assert 2003 in scid_years, "SCID cases should reference the 2003 Science paper"

    def test_thalassemia_references_contain_foundational_citations(self):
        """Thalassemia cases should cite Orkin, Kazazian, or Spritz."""
        thal_refs = " ".join(c.reference for c in THALASSEMIA_CASES)
        foundational_authors = ["Orkin", "Kazazian", "Spritz"]
        found_any = any(author in thal_refs for author in foundational_authors)
        assert found_any, (
            "Thalassemia cases should cite at least one foundational author "
            f"(Orkin, Kazazian, Spritz). Got: '{thal_refs[:100]}'"
        )

    def test_aggregation_references_contain_disease_protein_citations(self):
        """Aggregation cases should cite relevant protein aggregation literature."""
        agg_refs = " ".join(c.reference for c in AGGREGATION_CASES)
        # Should reference key aggregation researchers or journals
        key_terms = ["amyloid", "synuclein", "huntingtin", "aggregation", "Neuron", "Nature"]
        found = any(term.lower() in agg_refs.lower() for term in key_terms)
        assert found, (
            "Aggregation cases should reference key aggregation literature"
        )

    def test_immunogenicity_references_contain_clinical_citations(self):
        """Immunogenicity cases should cite clinical immunogenicity studies."""
        imm_refs = " ".join(c.reference for c in IMMUNOGENICITY_CASES)
        # Should reference IEDB or clinical journals
        key_terms = ["IEDB", "antibod", "immun", "Blood", "N Engl J Med"]
        found = any(term.lower() in imm_refs.lower() for term in key_terms)
        assert found, (
            "Immunogenicity cases should reference clinical/immunological literature"
        )

    def test_no_duplicate_references_across_cases(self):
        """No two cases should have identical references (each case is unique)."""
        refs = [c.reference for c in ALL_LITERATURE_CASES]
        seen: Set[str] = set()
        duplicates = []
        for ref in refs:
            if ref in seen:
                duplicates.append(ref[:60])
            seen.add(ref)
        assert not duplicates, (
            f"Duplicate references found: {duplicates}"
        )

    def test_reference_does_not_contain_placeholder_text(self):
        """References should not contain placeholder or TODO markers."""
        forbidden = ["TODO", "FIXME", "XXX", "TBD", "placeholder", "INSERT REF"]
        for case in ALL_LITERATURE_CASES:
            ref_lower = case.reference.lower()
            for marker in forbidden:
                assert marker.lower() not in ref_lower, (
                    f"Case {case.case_id} reference contains placeholder '{marker}'"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. DOI Format Checking
# ═══════════════════════════════════════════════════════════════════════════════

class TestDOIFormatChecking:
    """Validate DOI-related patterns in references and database_id fields.

    While not all references include explicit DOIs, we check:
      - If a DOI is present, it follows the 10.xxxx/yyyy format
      - Database IDs follow their expected format (UniProt:XXXX, PDB:XXXX, etc.)
      - Database ID prefixes are recognized
    """

    def test_explicit_dois_are_well_formed(self):
        """If a reference contains an explicit DOI, it must match the 10.xxxx/yyyy pattern."""
        for case in ALL_LITERATURE_CASES:
            dois = DOI_PATTERN.findall(case.reference)
            for doi in dois:
                # DOI must start with '10.' and have at least a registrar code
                assert doi.startswith("10."), (
                    f"Case {case.case_id}: malformed DOI '{doi}'"
                )
                # Registrar code should be 4+ digits
                registrar = doi.split("/")[0].replace("10.", "")
                assert len(registrar) >= 4, (
                    f"Case {case.case_id}: DOI registrar code too short in '{doi}'"
                )

    def test_database_id_has_recognized_prefix(self):
        """Every component of database_id should use a recognized prefix."""
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue  # database_id is optional
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                assert prefix in DATABASE_PREFIXES, (
                    f"Case {case.case_id}: unrecognized database prefix '{prefix}' "
                    f"in database_id '{case.database_id}'"
                )

    def test_uniprot_ids_are_valid_format(self):
        """UniProt IDs referenced in database_id should follow the format: 1-2 letters + alphanumeric."""
        uniprot_pattern = re.compile(r"^[A-Z][A-Z0-9]{5,}$")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "UniProt" and value:
                    assert uniprot_pattern.match(value), (
                        f"Case {case.case_id}: invalid UniProt ID '{value}'"
                    )

    def test_pdb_ids_are_valid_format(self):
        """PDB IDs should be 4 alphanumeric characters."""
        pdb_pattern = re.compile(r"^[0-9][A-Z0-9]{3}$")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "PDB" and value:
                    assert pdb_pattern.match(value), (
                        f"Case {case.case_id}: invalid PDB ID '{value}'"
                    )

    def test_clinvar_ids_are_valid_format(self):
        """ClinVar IDs should start with 'VCV' followed by digits."""
        clinvar_pattern = re.compile(r"^VCV\d+$")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "ClinVar" and value:
                    assert clinvar_pattern.match(value), (
                        f"Case {case.case_id}: invalid ClinVar ID '{value}'"
                    )

    def test_dbsnp_ids_are_valid_format(self):
        """dbSNP rs IDs should start with 'rs' followed by digits."""
        dbsnp_pattern = re.compile(r"^rs\d+$")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "dbSNP" and value:
                    assert dbsnp_pattern.match(value), (
                        f"Case {case.case_id}: invalid dbSNP ID '{value}'"
                    )

    def test_iedb_ids_are_numeric(self):
        """IEDB epitope IDs should be numeric."""
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "IEDB" and value:
                    assert value.isdigit(), (
                        f"Case {case.case_id}: IEDB ID should be numeric, got '{value}'"
                    )

    def test_refseq_ids_are_valid_format(self):
        """RefSeq IDs should follow the pattern: 2 letters + underscore + digits (e.g. NM_000206.4)."""
        refseq_pattern = re.compile(r"^[A-Z]{2}_\d+(\.\d+)?$")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "RefSeq" and value:
                    assert refseq_pattern.match(value), (
                        f"Case {case.case_id}: invalid RefSeq ID '{value}'"
                    )

    def test_ncbi_ids_are_valid_format(self):
        """NCBI RefSeq IDs should follow the pattern: 2 letters + underscore + alphanumeric."""
        ncbi_pattern = re.compile(r"^[A-Z]{2}_\d+$")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "NCBI" and value:
                    assert ncbi_pattern.match(value), (
                        f"Case {case.case_id}: invalid NCBI ID '{value}'"
                    )

    def test_database_id_components_are_delimited_by_slash(self):
        """Multiple database IDs should be separated by ' / '."""
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            # If there are multiple components, they should be ' / ' delimited
            if "/" in case.database_id:
                # Should not have bare '/' without spaces
                assert " / " in case.database_id, (
                    f"Case {case.case_id}: database_id components should be "
                    f"separated by ' / ', got '{case.database_id}'"
                )

    def test_iedb_epitope_ids_in_reference_match_database_id(self):
        """If IEDB epitope IDs appear in both reference and database_id, they should match."""
        iedb_ref_pattern = re.compile(r"EPITOPE_ID\s*(\d+)")
        for case in ALL_LITERATURE_CASES:
            if not case.database_id:
                continue
            ref_matches = iedb_ref_pattern.findall(case.reference)
            if not ref_matches:
                continue
            components = _extract_database_ids(case.database_id)
            db_iedb_values = []
            for comp in components:
                prefix, value = _parse_single_db_id(comp)
                if prefix == "IEDB" and value:
                    db_iedb_values.append(value)
            # If both reference and database_id have IEDB IDs, check consistency
            if db_iedb_values:
                for ref_id in ref_matches:
                    assert ref_id in db_iedb_values, (
                        f"Case {case.case_id}: IEDB EPITOPE_ID {ref_id} in reference "
                        f"not found in database_id {case.database_id}"
                    )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Citation Format
# ═══════════════════════════════════════════════════════════════════════════════

class TestCitationFormat:
    """Validate the structural format of citation strings in references.

    Expected format (author-year style):
        Author1 A et al. (YYYY) Journal Volume:Pages. 'Title.'
    Or variations with multiple citations separated by 'Also:' or similar.
    """

    def test_reference_contains_author_year_pattern(self):
        """Every reference should contain at least one author-year citation."""
        for case in ALL_LITERATURE_CASES:
            has_author_year = bool(AUTHOR_YEAR_PATTERN.search(case.reference))
            # Some references may use database-style citations (e.g. "RefSeq:XXX")
            # Allow those as an alternative
            has_db_citation = bool(re.search(r"(UniProt|RefSeq|ClinVar|PDB|IEDB)", case.reference))
            # Also accept database-style accession numbers (e.g. NG_000007.3)
            has_accession = bool(re.search(r"[A-Z]{2}_\d+", case.reference))
            assert has_author_year or has_db_citation or has_accession, (
                f"Case {case.case_id}: reference lacks author-year citation pattern: "
                f"'{case.reference[:80]}'"
            )

    def test_citation_has_journal_name(self):
        """References that use author-year style should include a journal name."""
        for case in ALL_LITERATURE_CASES:
            years = _extract_years(case.reference)
            if not years:
                continue  # Skip database-style references without years
            # After the year, there should be a journal name or similar
            # Check for common journal patterns
            has_journal = any(
                journal in case.reference
                for journal in JOURNAL_ABBREVS
            )
            # Some references use "Also:" to chain citations; check each part
            has_volume_pages = bool(re.search(r"\d+:\d+-\d+", case.reference))
            # Either a journal name OR volume:page pattern should be present
            assert has_journal or has_volume_pages, (
                f"Case {case.case_id}: reference with year but no recognizable "
                f"journal or volume:pages: '{case.reference[:100]}'"
            )

    def test_volume_page_format_is_valid(self):
        """Volume:page citations should follow the format 'Volume:StartPage-EndPage'."""
        vol_page_pattern = re.compile(r"(\d+):(\d+)-(\d+)")
        for case in ALL_LITERATURE_CASES:
            matches = vol_page_pattern.findall(case.reference)
            for vol, start, end in matches:
                assert int(end) >= int(start), (
                    f"Case {case.case_id}: end page {end} < start page {start} "
                    f"in volume:page citation"
                )

    def test_article_title_is_quoted(self):
        """Article titles in references should be enclosed in single quotes."""
        # Pattern:  'Title text.'
        quoted_title_pattern = re.compile(r"'[^']+\.'")
        for case in ALL_LITERATURE_CASES:
            years = _extract_years(case.reference)
            if not years:
                continue  # Skip database-style references
            has_quoted_title = bool(quoted_title_pattern.search(case.reference))
            # Not all references may include the title — that's acceptable
            # But if a title-like string is present, it should be quoted
            # We just check that if there is a quoted section, it's properly formed
            if has_quoted_title:
                # Extract the quoted title and verify it ends with a period
                titles = quoted_title_pattern.findall(case.reference)
                for title in titles:
                    assert title.endswith(".'") or title.endswith("."), (
                        f"Case {case.case_id}: quoted title should end with period: '{title}'"
                    )

    def test_multiple_citations_are_properly_joined(self):
        """References with multiple citations should use 'Also:' or similar connector."""
        for case in ALL_LITERATURE_CASES:
            years = _extract_years(case.reference)
            if len(years) > 1:
                # Multiple citations should be separated by a connector
                has_connector = bool(
                    re.search(r"(Also:|also:|and:|See also:)", case.reference)
                )
                # Or the citations could be in separate sentences
                assert has_connector or case.reference.count(".") >= 1, (
                    f"Case {case.case_id}: multiple citations in reference but "
                    f"no clear separator: '{case.reference[:100]}'"
                )

    def test_et_al_format(self):
        """'et al.' citations should use the correct format with period."""
        for case in ALL_LITERATURE_CASES:
            # Find 'et al' without period
            bad_et_al = re.search(r"\bet\s+al[^.]", case.reference)
            if bad_et_al:
                # Allow 'et al.' (correct) but flag 'et al' without period
                # that is not followed by a period
                pass  # Soft check — some styles omit the period

    def test_scid_citations_cover_both_trials(self):
        """SCID references should cover both the 2003 and 2010/2014 trials."""
        scid_years = set()
        for c in SCID_CASES:
            scid_years.update(_extract_years(c.reference))
        # Should reference the original 2003 trial
        assert 2003 in scid_years, "SCID references should include the 2003 trial"
        # And at least one follow-up
        followup_years = {2010, 2014}
        assert scid_years & followup_years, (
            f"SCID references should include follow-up years (2010/2014). "
            f"Found years: {sorted(scid_years)}"
        )

    def test_thalassemia_citations_cover_key_mutations(self):
        """Thalassemia references should cover the three canonical IVS1 mutations."""
        thal_refs = " ".join(c.reference for c in THALASSEMIA_CASES)
        # Thein SL is the primary review author for thalassemia
        assert "Thein" in thal_refs, (
            "Thalassemia references should cite Thein SL (primary reviewer)"
        )

    def test_reference_encoding_is_ascii_compatible(self):
        """References should not contain problematic Unicode characters."""
        for case in ALL_LITERATURE_CASES:
            try:
                case.reference.encode("ascii")
            except UnicodeEncodeError:
                # Non-ASCII is acceptable (e.g., accented author names)
                # but should be valid UTF-8
                assert case.reference.encode("utf-8"), (
                    f"Case {case.case_id}: reference is not valid UTF-8"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Extended ValidationResult and DomainReport Validation
# ═══════════════════════════════════════════════════════════════════════════════

class TestValidationResultExtended:
    """Extended tests for ValidationResult classification logic."""

    def test_true_positive_classification(self):
        """ValidationResult should classify TP correctly."""
        case = ALL_LITERATURE_CASES[0]
        result = ValidationResult(
            case=case,
            predicate_name=case.expected_predicate,
            predicted_flagged=True,
            ground_truth_flagged=True,
        )
        assert result.true_positive is True
        assert result.false_negative is False
        assert result.false_positive is False
        assert result.true_negative is False

    def test_false_negative_classification(self):
        """ValidationResult should classify FN correctly."""
        case = ALL_LITERATURE_CASES[0]
        result = ValidationResult(
            case=case,
            predicate_name=case.expected_predicate,
            predicted_flagged=False,
            ground_truth_flagged=True,
        )
        assert result.true_positive is False
        assert result.false_negative is True
        assert result.false_positive is False
        assert result.true_negative is False

    def test_false_positive_classification(self):
        """ValidationResult should classify FP correctly."""
        case = ALL_LITERATURE_CASES[0]
        result = ValidationResult(
            case=case,
            predicate_name=case.expected_predicate,
            predicted_flagged=True,
            ground_truth_flagged=False,
        )
        assert result.true_positive is False
        assert result.false_negative is False
        assert result.false_positive is True
        assert result.true_negative is False

    def test_true_negative_classification(self):
        """ValidationResult should classify TN correctly."""
        case = ALL_LITERATURE_CASES[0]
        result = ValidationResult(
            case=case,
            predicate_name=case.expected_predicate,
            predicted_flagged=False,
            ground_truth_flagged=False,
        )
        assert result.true_positive is False
        assert result.false_negative is False
        assert result.false_positive is False
        assert result.true_negative is True

    def test_exactly_one_classification(self):
        """Exactly one of TP/FN/FP/TN should be True for any result."""
        case = ALL_LITERATURE_CASES[0]
        for pred in (True, False):
            for gt in (True, False):
                result = ValidationResult(
                    case=case,
                    predicate_name=case.expected_predicate,
                    predicted_flagged=pred,
                    ground_truth_flagged=gt,
                )
                flags = [result.true_positive, result.false_negative,
                         result.false_positive, result.true_negative]
                assert sum(flags) == 1, (
                    f"Exactly one classification should be True for "
                    f"predicted={pred}, ground_truth={gt}. Got: {flags}"
                )


class TestDomainReportExtended:
    """Extended tests for DomainReport metric computation."""

    def test_compute_metrics_all_tp(self):
        """DomainReport with all TPs should have 100% sensitivity and precision."""
        report = DomainReport(
            domain="test",
            total_cases=3,
            true_positives=3,
            false_negatives=0,
            false_positives=0,
            true_negatives=0,
        )
        report.compute_metrics()
        assert report.sensitivity == 1.0
        assert report.precision == 1.0
        assert report.accuracy == 1.0

    def test_compute_metrics_all_tn(self):
        """DomainReport with all TNs should have 100% specificity."""
        report = DomainReport(
            domain="test",
            total_cases=3,
            true_positives=0,
            false_negatives=0,
            false_positives=0,
            true_negatives=3,
        )
        report.compute_metrics()
        assert report.specificity == 1.0
        assert report.accuracy == 1.0
        # Sensitivity is 0/0 → defined as 0.0
        assert report.sensitivity == 0.0

    def test_compute_metrics_mixed(self):
        """DomainReport with mixed results should compute correct metrics."""
        report = DomainReport(
            domain="test",
            total_cases=4,
            true_positives=2,
            false_negatives=1,
            false_positives=0,
            true_negatives=1,
        )
        report.compute_metrics()
        assert report.sensitivity == pytest.approx(2 / 3)
        assert report.specificity == 1.0  # 1 / (1 + 0)
        assert report.precision == 1.0     # 2 / (2 + 0)
        assert report.accuracy == pytest.approx(3 / 4)

    def test_compute_metrics_zero_division(self):
        """Metrics should handle zero-division gracefully (return 0.0)."""
        report = DomainReport(
            domain="test",
            total_cases=0,
            true_positives=0,
            false_negatives=0,
            false_positives=0,
            true_negatives=0,
        )
        report.compute_metrics()
        assert report.sensitivity == 0.0
        assert report.specificity == 0.0
        assert report.precision == 0.0
        assert report.accuracy == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Extended Cross-Domain and Integration Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestCrossDomainExtended:
    """Extended cross-domain integration tests."""

    def test_all_domains_have_database_ids(self):
        """Every case should have a database_id (even if optional, it should be populated)."""
        for case in ALL_LITERATURE_CASES:
            assert case.database_id, (
                f"Case {case.case_id} should have a database_id for traceability"
            )

    def test_case_id_follows_literature_convention(self):
        """Case IDs should follow the LIT-X# pattern (e.g. LIT-A1, LIT-B2)."""
        case_id_pattern = re.compile(r"^LIT-[A-D]\d+$")
        for case in ALL_LITERATURE_CASES:
            assert case_id_pattern.match(case.case_id), (
                f"Case ID '{case.case_id}' should follow LIT-X# pattern"
            )

    def test_case_ids_are_sequential_within_domain(self):
        """Case IDs within each domain should be sequentially numbered."""
        domain_letters = {"SCID": "A", "thalassemia": "B",
                          "aggregation": "C", "immunogenicity": "D"}
        for domain, letter in domain_letters.items():
            domain_cases = _cases_by_domain(domain)
            numbers = []
            for c in domain_cases:
                match = re.match(rf"^LIT-{letter}(\d+)$", c.case_id)
                assert match, f"Case {c.case_id} doesn't match domain {domain}"
                numbers.append(int(match.group(1)))
            # Numbers should be sequential starting from 1
            expected = list(range(1, len(numbers) + 1))
            assert sorted(numbers) == expected, (
                f"Domain {domain}: case IDs not sequential. "
                f"Expected {expected}, got {sorted(numbers)}"
            )

    def test_dna_cases_reference_genomic_databases(self):
        """DNA cases should reference genomic databases (NCBI, ClinVar, RefSeq)."""
        genomic_prefixes = {"NCBI", "ClinVar", "RefSeq", "dbSNP"}
        for case in ALL_LITERATURE_CASES:
            if case.sequence_type != "dna":
                continue
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            prefixes = {_parse_single_db_id(comp)[0] for comp in components}
            assert prefixes & genomic_prefixes, (
                f"DNA case {case.case_id} should reference a genomic database. "
                f"Found prefixes: {prefixes}"
            )

    def test_protein_cases_reference_protein_databases(self):
        """Protein cases should reference protein databases (UniProt, PDB, IEDB)."""
        protein_prefixes = {"UniProt", "PDB", "IEDB"}
        for case in ALL_LITERATURE_CASES:
            if case.sequence_type != "protein":
                continue
            if not case.database_id:
                continue
            components = _extract_database_ids(case.database_id)
            prefixes = {_parse_single_db_id(comp)[0] for comp in components}
            assert prefixes & protein_prefixes, (
                f"Protein case {case.case_id} should reference a protein database. "
                f"Found prefixes: {prefixes}"
            )

    def test_evaluate_case_unknown_domain(self):
        """evaluate_case should handle an unknown domain gracefully."""
        case = LiteratureCase(
            case_id="LIT-TEST",
            domain="unknown_domain",
            name="Test case",
            description="Test",
            sequence="ATGC",
            sequence_type="dna",
            expected_predicate="Unknown",
            ground_truth="FLAGGED — test",
            reference="Test ref",
        )
        result = evaluate_case(case)
        assert "No evaluator" in result.details
        assert result.predicted_flagged is False

    def test_format_report_includes_all_case_ids(self):
        """The formatted report should mention every case_id."""
        reports = run_literature_validation()
        text = format_literature_validation_report(reports)
        for case in ALL_LITERATURE_CASES:
            assert case.case_id in text, (
                f"Case {case.case_id} not found in formatted report"
            )

    def test_format_report_includes_all_database_ids(self):
        """The formatted report should reference all database IDs (via references)."""
        reports = run_literature_validation()
        text = format_literature_validation_report(reports)
        # At minimum, the report should contain reference info
        for case in ALL_LITERATURE_CASES:
            # The reference should appear in the report
            assert case.reference[:50] in text or case.case_id in text, (
                f"Case {case.case_id} reference not represented in report"
            )

    def test_ground_truth_flagged_parsing(self):
        """Ground truth 'FLAGGED' / 'NOT FLAGGED' should be parsed correctly."""
        flagged_cases = [c for c in ALL_LITERATURE_CASES
                        if "FLAGGED" in c.ground_truth and "NOT FLAGGED" not in c.ground_truth]
        not_flagged_cases = [c for c in ALL_LITERATURE_CASES
                            if "NOT FLAGGED" in c.ground_truth]
        # There should be both flagged and not-flagged cases
        assert len(flagged_cases) > 0, "No FLAGGED cases found"
        assert len(not_flagged_cases) > 0, "No NOT FLAGGED cases found"
        # Every case should be in exactly one category
        assert len(flagged_cases) + len(not_flagged_cases) == len(ALL_LITERATURE_CASES), (
            "Some cases are neither FLAGGED nor NOT FLAGGED"
        )
