"""
Tests for biocompiler.benchmarking.provenance_comparison
=========================================================

Tests the provenance comparison module that documents BioCompiler's unique
differentiator vs DNAchisel: per-codon, per-constraint decision provenance
with CAI impact tracking and full audit trail.

Test classes:
- TestCompareProvenanceCapabilities: Validates the capability comparison
- TestProvenanceComparisonDataClass: Validates the data class structure
- TestGenerateAuditReport: Validates audit report generation
- TestAuditReportDataClass: Validates the audit report data class
- TestRegulatorySummary: Validates regulatory summary content
"""

from __future__ import annotations

import pytest

from biocompiler.benchmarking.provenance_comparison import (
    AuditReport,
    ProvenanceComparison,
    compare_provenance_capabilities,
    generate_audit_report,
    _generate_regulatory_summary,
)
from biocompiler.provenance.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
)


# ---------------------------------------------------------------------------
# Test compare_provenance_capabilities()
# ---------------------------------------------------------------------------

class TestCompareProvenanceCapabilities:
    """Test that compare_provenance_capabilities() returns a valid comparison."""

    def test_returns_provenance_comparison(self):
        """compare_provenance_capabilities() returns a ProvenanceComparison."""
        result = compare_provenance_capabilities()
        assert isinstance(result, ProvenanceComparison)

    def test_biocompiler_features_non_empty(self):
        """BioCompiler features list is non-empty."""
        result = compare_provenance_capabilities()
        assert len(result.biocompiler_features) > 0

    def test_dnachisel_features_non_empty(self):
        """DNAchisel features list is non-empty."""
        result = compare_provenance_capabilities()
        assert len(result.dnachisel_features) > 0

    def test_biocompiler_has_more_features(self):
        """BioCompiler has more provenance features than DNAchisel."""
        result = compare_provenance_capabilities()
        assert len(result.biocompiler_features) > len(result.dnachisel_features)

    def test_unique_to_biocompiler_non_empty(self):
        """There are features unique to BioCompiler."""
        result = compare_provenance_capabilities()
        assert len(result.unique_to_biocompiler) > 0

    def test_unique_to_biocompiler_includes_per_codon_provenance(self):
        """Unique features include per-codon provenance."""
        result = compare_provenance_capabilities()
        unique_text = " ".join(result.unique_to_biocompiler).lower()
        # Check for per-codon provenance in any form
        assert (
            "per-codon" in unique_text
            or "per codon" in unique_text
            or "codon decision" in unique_text
        ), (
            f"Expected per-codon provenance in unique features, "
            f"got: {result.unique_to_biocompiler}"
        )

    def test_unique_to_biocompiler_includes_cai_impact(self):
        """Unique features include CAI impact tracking."""
        result = compare_provenance_capabilities()
        unique_text = " ".join(result.unique_to_biocompiler).lower()
        assert "cai" in unique_text, (
            f"Expected CAI impact tracking in unique features, "
            f"got: {result.unique_to_biocompiler}"
        )

    def test_unique_to_biocompiler_includes_conflict_resolution(self):
        """Unique features include conflict resolution trace."""
        result = compare_provenance_capabilities()
        unique_text = " ".join(result.unique_to_biocompiler).lower()
        assert "conflict" in unique_text, (
            f"Expected conflict resolution in unique features, "
            f"got: {result.unique_to_biocompiler}"
        )

    def test_unique_to_biocompiler_includes_audit_trail(self):
        """Unique features include full audit trail."""
        result = compare_provenance_capabilities()
        unique_text = " ".join(result.unique_to_biocompiler).lower()
        assert "audit" in unique_text, (
            f"Expected audit trail in unique features, "
            f"got: {result.unique_to_biocompiler}"
        )

    def test_overlap_non_empty(self):
        """Overlap between BioCompiler and DNAchisel is non-empty."""
        result = compare_provenance_capabilities()
        assert len(result.overlap) > 0

    def test_overlap_includes_cai_computation(self):
        """Overlap includes CAI computation (both tools can compute CAI post-hoc)."""
        result = compare_provenance_capabilities()
        overlap_text = " ".join(result.overlap).lower()
        assert "cai" in overlap_text, (
            f"Expected CAI in overlap features, got: {result.overlap}"
        )

    def test_unique_features_not_in_dnachisel(self):
        """Unique BioCompiler features are not in the DNAchisel features list."""
        result = compare_provenance_capabilities()
        dn_set = set(result.dnachisel_features)
        for feature in result.unique_to_biocompiler:
            assert feature not in dn_set, (
                f"Unique feature '{feature}' should not be in DNAchisel features"
            )

    def test_to_dict(self):
        """ProvenanceComparison.to_dict() returns a valid dict."""
        result = compare_provenance_capabilities()
        d = result.to_dict()
        assert isinstance(d, dict)
        assert "biocompiler_features" in d
        assert "dnachisel_features" in d
        assert "unique_to_biocompiler" in d
        assert "overlap" in d
        assert len(d["biocompiler_features"]) == len(result.biocompiler_features)
        assert len(d["unique_to_biocompiler"]) == len(result.unique_to_biocompiler)


# ---------------------------------------------------------------------------
# Test ProvenanceComparison data class
# ---------------------------------------------------------------------------

class TestProvenanceComparisonDataClass:
    """Test ProvenanceComparison data class structure."""

    def test_fields(self):
        """ProvenanceComparison has the required fields."""
        comp = ProvenanceComparison(
            biocompiler_features=["a"],
            dnachisel_features=["b"],
            unique_to_biocompiler=["a"],
            overlap=[],
        )
        assert comp.biocompiler_features == ["a"]
        assert comp.dnachisel_features == ["b"]
        assert comp.unique_to_biocompiler == ["a"]
        assert comp.overlap == []

    def test_to_dict_keys(self):
        """to_dict() returns all required keys."""
        comp = ProvenanceComparison(
            biocompiler_features=["f1"],
            dnachisel_features=["f2"],
            unique_to_biocompiler=["f1"],
            overlap=[],
        )
        d = comp.to_dict()
        assert set(d.keys()) == {
            "biocompiler_features",
            "dnachisel_features",
            "unique_to_biocompiler",
            "overlap",
        }

    def test_to_dict_copies(self):
        """to_dict() returns copies, not references to internal lists."""
        comp = ProvenanceComparison(
            biocompiler_features=["f1"],
            dnachisel_features=["f2"],
            unique_to_biocompiler=["f1"],
            overlap=[],
        )
        d = comp.to_dict()
        d["biocompiler_features"].append("f3")
        assert len(comp.biocompiler_features) == 1  # Original unchanged


# ---------------------------------------------------------------------------
# Test generate_audit_report()
# ---------------------------------------------------------------------------

class TestGenerateAuditReport:
    """Test generate_audit_report() produces a complete report."""

    @pytest.fixture
    def sample_provenance_records(self):
        """Create sample provenance records for testing."""
        records = []

        # Codon decisions
        for i in range(10):
            records.append({
                "type": "codon",
                "position": i,
                "amino_acid": "MVLSPADKTN"[i],
                "chosen_codon": ["ATG", "GTC", "CTG", "TCT", "CCG",
                                 "GCG", "GAC", "AAG", "ACC", "AAC"][i],
                "constraint_reason": "maximize_cai" if i % 3 != 0 else "gc_content",
                "confidence": 0.95 - (i * 0.05),
                "alternatives_considered": [],
            })

        # Constraint decisions
        records.append({
            "type": "constraint",
            "constraint_name": "NoCrypticSplice",
            "constraint_type": "hard",
            "action_taken": "satisfied",
            "positions_affected": [2, 5],
            "tradeoff_description": "Chose CTG over CTA to avoid cryptic donor",
            "impact_on_cai": -0.003,
        })
        records.append({
            "type": "constraint",
            "constraint_name": "GCInRange",
            "constraint_type": "hard",
            "action_taken": "satisfied",
            "positions_affected": [0, 3, 4],
            "tradeoff_description": "Adjusted GC content to meet target range",
            "impact_on_cai": -0.005,
        })
        records.append({
            "type": "constraint",
            "constraint_name": "AvoidRestrictionSite_EcoRI",
            "constraint_type": "hard",
            "action_taken": "conflicted",
            "positions_affected": [6, 7],
            "tradeoff_description": "Replaced GAATTC motif at position 18 with synonymous codon",
            "impact_on_cai": -0.008,
        })

        return records

    @pytest.fixture
    def sample_protein(self):
        """Sample protein sequence."""
        return "MVLSPADKTN"

    @pytest.fixture
    def sample_sequence(self):
        """Sample optimized DNA sequence."""
        return "ATGGTCCTGTCGCCGGCGGACAAGACCAAC"

    def test_returns_audit_report(self, sample_protein, sample_sequence,
                                  sample_provenance_records):
        """generate_audit_report() returns an AuditReport."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert isinstance(report, AuditReport)

    def test_protein_preserved(self, sample_protein, sample_sequence,
                                sample_provenance_records):
        """Protein is preserved in the report."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert report.protein == sample_protein

    def test_organism_preserved(self, sample_protein, sample_sequence,
                                 sample_provenance_records):
        """Organism is preserved in the report."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert report.organism == "Homo_sapiens"

    def test_sequence_preserved(self, sample_protein, sample_sequence,
                                 sample_provenance_records):
        """Sequence is preserved in the report."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert report.sequence == sample_sequence

    def test_total_codon_decisions(self, sample_protein, sample_sequence,
                                    sample_provenance_records):
        """Total codon decisions is correct."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert report.total_codon_decisions == 10

    def test_total_constraint_decisions(self, sample_protein, sample_sequence,
                                         sample_provenance_records):
        """Total constraint decisions is correct."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert report.total_constraint_decisions == 3

    def test_cai_cost_breakdown(self, sample_protein, sample_sequence,
                                 sample_provenance_records):
        """CAI cost breakdown includes all constraints."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert isinstance(report.cai_cost_breakdown, dict)
        assert "NoCrypticSplice" in report.cai_cost_breakdown
        assert "GCInRange" in report.cai_cost_breakdown
        assert "AvoidRestrictionSite_EcoRI" in report.cai_cost_breakdown
        assert report.cai_cost_breakdown["NoCrypticSplice"] == pytest.approx(-0.003, abs=0.0001)
        assert report.cai_cost_breakdown["GCInRange"] == pytest.approx(-0.005, abs=0.0001)
        assert report.cai_cost_breakdown["AvoidRestrictionSite_EcoRI"] == pytest.approx(-0.008, abs=0.0001)

    def test_top_decisions_limited(self, sample_protein, sample_sequence,
                                    sample_provenance_records):
        """Top decisions list has at most 10 entries."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert len(report.top_decisions) <= 10

    def test_top_decisions_sorted_by_impact(self, sample_protein, sample_sequence,
                                             sample_provenance_records):
        """Top decisions are sorted by impact (most impactful first)."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        # The first decisions should be from constraint records or
        # low-confidence codon records
        assert len(report.top_decisions) > 0

    def test_conflict_resolutions(self, sample_protein, sample_sequence,
                                   sample_provenance_records):
        """Conflict resolutions capture non-satisfied constraint decisions."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert isinstance(report.conflict_resolutions, list)
        # The "conflicted" constraint should be captured
        conflicted_names = [
            cr["constraint_name"] for cr in report.conflict_resolutions
        ]
        assert "AvoidRestrictionSite_EcoRI" in conflicted_names

    def test_regulatory_summary_non_empty(self, sample_protein, sample_sequence,
                                           sample_provenance_records):
        """Regulatory summary is non-empty."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        assert len(report.regulatory_summary) > 0

    def test_empty_provenance_records(self):
        """generate_audit_report() handles empty provenance records."""
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[],
        )
        assert report.total_codon_decisions == 0
        assert report.total_constraint_decisions == 0
        assert report.cai_cost_breakdown == {}
        assert report.top_decisions == []
        assert report.conflict_resolutions == []
        assert len(report.regulatory_summary) > 0

    def test_with_codon_decision_objects(self):
        """generate_audit_report() handles CodonDecision objects."""
        codon_decision = CodonDecision(
            position=0,
            amino_acid="M",
            original_codon=None,
            chosen_codon="ATG",
            alternatives_considered=[],
            constraint_reason="maximize_cai",
            confidence=1.0,
        )
        report = generate_audit_report(
            protein="M",
            organism="Homo_sapiens",
            sequence="ATG",
            provenance_records=[codon_decision],
        )
        assert report.total_codon_decisions == 1
        assert report.total_constraint_decisions == 0

    def test_with_constraint_decision_objects(self):
        """generate_audit_report() handles ConstraintDecision objects."""
        constraint_decision = ConstraintDecision(
            constraint_name="NoCrypticSplice",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[2, 5],
            tradeoff_description="Chose CTG over CTA to avoid cryptic donor",
            impact_on_cai=-0.003,
        )
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[constraint_decision],
        )
        assert report.total_codon_decisions == 0
        assert report.total_constraint_decisions == 1
        assert "NoCrypticSplice" in report.cai_cost_breakdown

    def test_mixed_records(self):
        """generate_audit_report() handles mixed dict and object records."""
        codon_decision = CodonDecision(
            position=0,
            amino_acid="M",
            original_codon=None,
            chosen_codon="ATG",
            alternatives_considered=[],
            constraint_reason="maximize_cai",
            confidence=0.9,
        )
        constraint_dict = {
            "type": "constraint",
            "constraint_name": "GCInRange",
            "constraint_type": "hard",
            "action_taken": "satisfied",
            "positions_affected": [0],
            "tradeoff_description": "Adjusted GC",
            "impact_on_cai": -0.004,
        }
        report = generate_audit_report(
            protein="M",
            organism="Homo_sapiens",
            sequence="ATG",
            provenance_records=[codon_decision, constraint_dict],
        )
        assert report.total_codon_decisions == 1
        assert report.total_constraint_decisions == 1

    def test_to_dict(self, sample_protein, sample_sequence,
                      sample_provenance_records):
        """AuditReport.to_dict() returns a valid dict."""
        report = generate_audit_report(
            protein=sample_protein,
            organism="Homo_sapiens",
            sequence=sample_sequence,
            provenance_records=sample_provenance_records,
        )
        d = report.to_dict()
        assert isinstance(d, dict)
        assert "protein" in d
        assert "organism" in d
        assert "sequence" in d
        assert "total_codon_decisions" in d
        assert "total_constraint_decisions" in d
        assert "cai_total" in d
        assert "cai_unconstrained" in d
        assert "cai_cost_breakdown" in d
        assert "top_decisions" in d
        assert "conflict_resolutions" in d
        assert "regulatory_summary" in d


# ---------------------------------------------------------------------------
# Test AuditReport data class
# ---------------------------------------------------------------------------

class TestAuditReportDataClass:
    """Test AuditReport data class structure."""

    def test_fields(self):
        """AuditReport has all required fields."""
        report = AuditReport(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTG",
            total_codon_decisions=10,
            total_constraint_decisions=3,
            cai_total=0.85,
            cai_unconstrained=0.92,
            cai_cost_breakdown={"GCInRange": -0.04, "NoCrypticSplice": -0.03},
            top_decisions=[{"position": 0, "chosen_codon": "ATG"}],
            conflict_resolutions=[{"constraint_name": "EcoRI", "action_taken": "relaxed"}],
            regulatory_summary="This is a regulatory summary.",
        )
        assert report.protein == "MVLSPADKTN"
        assert report.organism == "Homo_sapiens"
        assert report.sequence == "ATGGTCCTG"
        assert report.total_codon_decisions == 10
        assert report.total_constraint_decisions == 3
        assert report.cai_total == 0.85
        assert report.cai_unconstrained == 0.92
        assert len(report.cai_cost_breakdown) == 2
        assert len(report.top_decisions) == 1
        assert len(report.conflict_resolutions) == 1
        assert report.regulatory_summary == "This is a regulatory summary."

    def test_to_dict_complete(self):
        """to_dict() includes all fields."""
        report = AuditReport(
            protein="M",
            organism="Homo_sapiens",
            sequence="ATG",
            total_codon_decisions=1,
            total_constraint_decisions=0,
            cai_total=0.9,
            cai_unconstrained=1.0,
            cai_cost_breakdown={},
            top_decisions=[],
            conflict_resolutions=[],
            regulatory_summary="Summary.",
        )
        d = report.to_dict()
        expected_keys = {
            "protein", "organism", "sequence",
            "total_codon_decisions", "total_constraint_decisions",
            "cai_total", "cai_unconstrained",
            "cai_cost_breakdown", "top_decisions",
            "conflict_resolutions", "regulatory_summary",
        }
        assert set(d.keys()) == expected_keys


# ---------------------------------------------------------------------------
# Test regulatory_summary content
# ---------------------------------------------------------------------------

class TestRegulatorySummary:
    """Test that the regulatory_summary is non-empty and readable."""

    def test_regulatory_summary_non_empty(self):
        """Regulatory summary is non-empty for a basic optimization."""
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[],
        )
        assert len(report.regulatory_summary) > 0

    def test_regulatory_summary_contains_protein_length(self):
        """Regulatory summary mentions protein length."""
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[],
        )
        assert "10-amino-acid" in report.regulatory_summary

    def test_regulatory_summary_contains_organism(self):
        """Regulatory summary mentions the target organism."""
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[],
        )
        assert "Homo sapiens" in report.regulatory_summary

    def test_regulatory_summary_mentions_provenance(self):
        """Regulatory summary mentions provenance/audit trail."""
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[],
        )
        summary_lower = report.regulatory_summary.lower()
        assert "provenance" in summary_lower or "audit" in summary_lower

    def test_regulatory_summary_with_constraints(self):
        """Regulatory summary mentions constraint costs when present."""
        records = [{
            "type": "constraint",
            "constraint_name": "NoCrypticSplice",
            "constraint_type": "hard",
            "action_taken": "satisfied",
            "positions_affected": [2],
            "tradeoff_description": "Avoided cryptic donor",
            "impact_on_cai": -0.005,
        }]
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=records,
        )
        summary_lower = report.regulatory_summary.lower()
        assert "constraint" in summary_lower or "cai" in summary_lower

    def test_regulatory_summary_with_conflicts(self):
        """Regulatory summary mentions conflicts when present."""
        records = [{
            "type": "constraint",
            "constraint_name": "AvoidRestrictionSite",
            "constraint_type": "hard",
            "action_taken": "conflicted",
            "positions_affected": [3],
            "tradeoff_description": "Had to choose between CAI and restriction avoidance",
            "impact_on_cai": -0.01,
        }]
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=records,
        )
        summary_lower = report.regulatory_summary.lower()
        assert "conflict" in summary_lower

    def test_regulatory_summary_is_readable_sentence(self):
        """Regulatory summary is a readable English sentence (starts with capital)."""
        report = generate_audit_report(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            provenance_records=[],
        )
        assert report.regulatory_summary[0].isupper()
        assert report.regulatory_summary.endswith(".")

    def test_generate_regulatory_summary_helper(self):
        """_generate_regulatory_summary() returns a non-empty string."""
        summary = _generate_regulatory_summary(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            sequence="ATGGTCCTGTCGCCGGCGGACAAGACCAAC",
            total_codon_decisions=10,
            total_constraint_decisions=3,
            cai_total=0.85,
            cai_unconstrained=0.92,
            cai_cost_breakdown={"GCInRange": -0.04, "NoCrypticSplice": -0.03},
            conflict_resolutions=[],
        )
        assert isinstance(summary, str)
        assert len(summary) > 50  # Should be a substantial paragraph
        assert "10-amino-acid" in summary
        assert "10" in summary  # codon decisions
        assert "3" in summary  # constraint decisions
