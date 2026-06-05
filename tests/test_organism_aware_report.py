"""
Tests for organism-aware report generation (Task F1.11).

Verifies that:
- Reports include organism domain information
- Skipped constraints are listed with reasons
- CAI comparison section is generated when data is available
- Provenance reports include organism domain
- Provenance reports note when constraints were not applied due to organism type
"""

import json
import pytest

from biocompiler.report import (
    generate_report,
    _build_organism_aware_constraints,
    _build_cai_comparison_section,
    _EUKARYOTE_ONLY_CONSTRAINTS,
)
from biocompiler.provenance_reporting import (
    ProvenanceReport,
    ProvenanceQuery,
)
from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_sequence() -> str:
    """A short valid DNA sequence (starts with ATG)."""
    return "ATGGCGGCAGCAGCTGCAGCTGG"  # 24 bp = 8 codons


@pytest.fixture
def eukaryotic_org() -> str:
    return "Homo_sapiens"


@pytest.fixture
def prokaryotic_org() -> str:
    return "E_coli_K12"


@pytest.fixture
def minimal_trail() -> OptimizationDecisionTrail:
    """Create a minimal OptimizationDecisionTrail for testing."""
    return OptimizationDecisionTrail(
        gene_name="test_gene",
        organism="E_coli_K12",
        solver_backend="greedy",
        seed=42,
        timestamp="2025-01-01T00:00:00Z",
        version="1.0",
        input_protein="MAAA",
        output_dna="ATGGCGGCAGCGCT",
        codon_decisions=[
            CodonDecision(
                position=0,
                amino_acid="M",
                chosen_codon="ATG",
                original_codon="ATG",
                constraint_reason="maximize_cai",
                confidence=1.0,
                alternatives_considered=[],
            ),
            CodonDecision(
                position=1,
                amino_acid="A",
                chosen_codon="GCG",
                original_codon="GCA",
                constraint_reason="maximize_cai",
                confidence=0.9,
                alternatives_considered=[{"codon": "GCA", "cai": 0.5}],
            ),
            CodonDecision(
                position=2,
                amino_acid="A",
                chosen_codon="GCC",
                original_codon="GCT",
                constraint_reason="gc_content",
                confidence=0.8,
                alternatives_considered=[{"codon": "GCT", "cai": 0.6}],
            ),
            CodonDecision(
                position=3,
                amino_acid="A",
                chosen_codon="GCT",
                original_codon="GCT",
                constraint_reason="maximize_cai",
                confidence=1.0,
                alternatives_considered=[],
            ),
        ],
        constraint_decisions=[
            ConstraintDecision(
                constraint_name="GCRangeConstraint",
                constraint_type="gc_content",
                action_taken="satisfied",
                impact_on_cai=-0.02,
                positions_affected=[2],
                tradeoff_description="Chose GCC over GCT to meet GC target",
            ),
        ],
        iteration_log=[],
        total_cai=0.85,
        total_gc=0.65,
    )


# ---------------------------------------------------------------------------
# Tests for _build_organism_aware_constraints
# ---------------------------------------------------------------------------

class TestBuildOrganismAwareConstraints:
    """Tests for the _build_organism_aware_constraints helper."""

    def test_eukaryotic_organism_has_no_skipped_constraints(self, eukaryotic_org):
        result = _build_organism_aware_constraints(eukaryotic_org)
        assert result["is_eukaryote"] is True
        assert result["domain"] == "eukaryote"
        assert len(result["skipped_constraints"]) == 0
        assert len(result["active_constraints"]) > 0

    def test_prokaryotic_organism_has_skipped_constraints(self, prokaryotic_org):
        result = _build_organism_aware_constraints(prokaryotic_org)
        assert result["is_eukaryote"] is False
        assert result["domain"] == "prokaryote"
        assert len(result["skipped_constraints"]) > 0
        skipped_names = [name for name, _ in result["skipped_constraints"]]
        for euk_constraint in _EUKARYOTE_ONLY_CONSTRAINTS:
            assert euk_constraint in skipped_names

    def test_prokaryotic_skipped_constraints_have_reasons(self, prokaryotic_org):
        result = _build_organism_aware_constraints(prokaryotic_org)
        for name, reason in result["skipped_constraints"]:
            assert "skipped" in reason.lower()
            # Reason should mention prokaryotic or the constraint type
            assert name in reason or "prokaryote" in reason.lower()

    def test_active_constraints_have_descriptions(self, eukaryotic_org):
        result = _build_organism_aware_constraints(eukaryotic_org)
        for name, desc in result["active_constraints"]:
            assert isinstance(name, str)
            assert isinstance(desc, str)
            assert len(desc) > 0

    def test_organism_config_name_populated(self, eukaryotic_org, prokaryotic_org):
        result_euk = _build_organism_aware_constraints(eukaryotic_org)
        assert result_euk["config_name"] == "Homo sapiens"

        result_prok = _build_organism_aware_constraints(prokaryotic_org)
        assert "Escherichia coli" in result_prok["config_name"]

    def test_solver_provided_constraints_take_precedence(self, prokaryotic_org):
        opt_result = {
            "applied_constraints": ["TranslationConstraint", "MaximizeCAI"],
            "skipped_constraints": ["NoCrypticSpliceConstraint"],
        }
        result = _build_organism_aware_constraints(prokaryotic_org, opt_result)
        assert len(result["active_constraints"]) == 2
        assert len(result["skipped_constraints"]) == 1

    def test_unknown_organism_defaults_to_eukaryote(self):
        result = _build_organism_aware_constraints("Unknown_organism_xyz")
        assert result["is_eukaryote"] is True
        assert result["domain"] == "eukaryote"
        assert len(result["skipped_constraints"]) == 0


# ---------------------------------------------------------------------------
# Tests for _build_cai_comparison_section
# ---------------------------------------------------------------------------

class TestBuildCAIComparisonSection:
    """Tests for the _build_cai_comparison_section helper."""

    def test_returns_none_when_no_optimization_result(self, eukaryotic_org):
        result = _build_cai_comparison_section(eukaryotic_org, None)
        assert result is None

    def test_returns_none_when_no_cai_data(self, eukaryotic_org):
        result = _build_cai_comparison_section(eukaryotic_org, {})
        assert result is None

    def test_returns_comparison_when_both_cai_values_present(self, eukaryotic_org):
        opt = {
            "cai_all_constraints": 0.78,
            "cai_organism_aware_constraints": 0.85,
        }
        result = _build_cai_comparison_section(eukaryotic_org, opt)
        assert result is not None
        assert abs(result["cai_all_constraints"] - 0.78) < 1e-6
        assert abs(result["cai_organism_aware"] - 0.85) < 1e-6
        assert abs(result["recovery"] - 0.07) < 1e-6

    def test_recovery_is_positive_for_organism_aware(self, prokaryotic_org):
        opt = {
            "cai_all_constraints": 0.80,
            "cai_organism_aware_constraints": 0.88,
        }
        result = _build_cai_comparison_section(prokaryotic_org, opt)
        assert result is not None
        assert result["recovery"] > 0

    def test_fallback_to_main_cai(self, eukaryotic_org):
        opt = {
            "cai_all_constraints": 0.75,
            "cai": 0.75,
        }
        result = _build_cai_comparison_section(eukaryotic_org, opt)
        # When cai_organism_aware is not provided, it falls back to "cai"
        assert result is not None
        assert abs(result["cai_organism_aware"] - 0.75) < 1e-6

    def test_invalid_cai_values_return_none(self, eukaryotic_org):
        opt = {
            "cai_all_constraints": "not_a_number",
            "cai_organism_aware_constraints": 0.85,
        }
        result = _build_cai_comparison_section(eukaryotic_org, opt)
        assert result is None


# ---------------------------------------------------------------------------
# Tests for generate_report (HTML output)
# ---------------------------------------------------------------------------

class TestGenerateReportOrganismAware:
    """Tests that the HTML report includes organism-aware information."""

    def test_report_includes_organism_domain_in_header(
        self, sample_sequence, eukaryotic_org
    ):
        html = generate_report(sample_sequence, organism=eukaryotic_org)
        assert "eukaryote" in html
        assert eukaryotic_org in html

    def test_report_includes_organism_aware_constraints_section(
        self, sample_sequence, eukaryotic_org
    ):
        html = generate_report(sample_sequence, organism=eukaryotic_org)
        assert "Organism-Aware Constraints" in html
        assert "Constraints Applied" in html

    def test_report_shows_active_constraints(
        self, sample_sequence, eukaryotic_org
    ):
        html = generate_report(sample_sequence, organism=eukaryotic_org)
        assert "Active" in html
        # Key constraints should appear in the report
        assert "TranslationConstraint" in html

    def test_report_for_prokaryote_shows_skipped_constraints(
        self, sample_sequence, prokaryotic_org
    ):
        html = generate_report(sample_sequence, organism=prokaryotic_org)
        assert "Skipped" in html
        assert "NoCrypticSpliceConstraint" in html
        assert "prokaryote" in html.lower()

    def test_report_includes_cai_comparison_when_data_available(
        self, sample_sequence, prokaryotic_org
    ):
        opt = {
            "cai": 0.85,
            "gc_content": 0.55,
            "satisfied": [],
            "failed": [],
            "cai_all_constraints": 0.78,
            "cai_organism_aware_constraints": 0.85,
        }
        html = generate_report(
            sample_sequence,
            organism=prokaryotic_org,
            optimization_result=opt,
        )
        assert "CAI Comparison" in html
        assert "CAI with all constraints" in html
        assert "CAI with organism-aware constraints" in html
        assert "Recovery" in html

    def test_report_skipped_constraint_shows_reason(
        self, sample_sequence, prokaryotic_org
    ):
        html = generate_report(sample_sequence, organism=prokaryotic_org)
        # The reason should mention eukaryotic organism
        assert "eukaryotic organism" in html.lower() or "prokaryotic" in html.lower()

    def test_report_no_cai_comparison_without_data(
        self, sample_sequence, eukaryotic_org
    ):
        html = generate_report(sample_sequence, organism=eukaryotic_org)
        # Without optimization_result providing both CAI values,
        # the CAI comparison section should not appear
        assert "CAI Comparison" not in html

    def test_report_domain_in_subtitle(
        self, sample_sequence, prokaryotic_org
    ):
        html = generate_report(sample_sequence, organism=prokaryotic_org)
        # The subtitle should show domain
        assert "prokaryote" in html


# ---------------------------------------------------------------------------
# Tests for ProvenanceReport (Markdown)
# ---------------------------------------------------------------------------

class TestProvenanceReportOrganismAware:
    """Tests that provenance reports include organism domain information."""

    def test_markdown_report_includes_organism_domain(self, minimal_trail):
        md = ProvenanceReport.generate_markdown(minimal_trail)
        assert "Organism domain" in md
        assert "prokaryote" in md

    def test_markdown_report_notes_skipped_constraints_for_prokaryote(
        self, minimal_trail
    ):
        md = ProvenanceReport.generate_markdown(minimal_trail)
        assert "Organism-Aware Constraint Filtering" in md
        assert "NoCrypticSpliceConstraint" in md
        assert "NoCpGIslandConstraint" in md
        assert "NoTRunConstraint" in md
        assert "skipped" in md.lower()

    def test_markdown_report_no_skipping_for_eukaryote(self, minimal_trail):
        minimal_trail.organism = "Homo_sapiens"
        md = ProvenanceReport.generate_markdown(minimal_trail)
        assert "Organism domain" in md
        assert "eukaryote" in md
        # Eukaryotic organisms should NOT have the filtering section
        assert "Organism-Aware Constraint Filtering" not in md

    def test_html_report_includes_domain(self, minimal_trail):
        html = ProvenanceReport.generate_html(minimal_trail)
        assert "Domain" in html
        assert "prokaryote" in html

    def test_html_report_notes_skipped_for_prokaryote(self, minimal_trail):
        html = ProvenanceReport.generate_html(minimal_trail)
        assert "Organism-Aware Constraint Filtering" in html
        assert "NoCrypticSpliceConstraint" in html

    def test_html_report_no_filtering_for_eukaryote(self, minimal_trail):
        minimal_trail.organism = "Homo_sapiens"
        html = ProvenanceReport.generate_html(minimal_trail)
        assert "Organism-Aware Constraint Filtering" not in html

    def test_json_report_includes_domain(self, minimal_trail):
        json_str = ProvenanceReport.generate_json(minimal_trail)
        data = json.loads(json_str)
        assert "organism_domain" in data
        assert data["organism_domain"] == "prokaryote"
        assert "organism_is_eukaryote" in data
        assert data["organism_is_eukaryote"] is False

    def test_json_report_includes_skipped_constraints(self, minimal_trail):
        json_str = ProvenanceReport.generate_json(minimal_trail)
        data = json.loads(json_str)
        assert "skipped_constraints_for_organism" in data
        skipped = data["skipped_constraints_for_organism"]
        assert "NoCrypticSpliceConstraint" in skipped
        assert "NoCpGIslandConstraint" in skipped
        assert "NoTRunConstraint" in skipped

    def test_json_report_eukaryote_empty_skipped(self, minimal_trail):
        minimal_trail.organism = "Homo_sapiens"
        json_str = ProvenanceReport.generate_json(minimal_trail)
        data = json.loads(json_str)
        assert data["organism_domain"] == "eukaryote"
        assert data["organism_is_eukaryote"] is True
        assert data["skipped_constraints_for_organism"] == []
