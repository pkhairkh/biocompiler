"""
Tests for CAI-aware provenance (Task F2.3).

Covers:
1. DecisionRecord captures cai_impact, codon_before, codon_after
2. DecisionRecord serialization/deserialization with CAI fields
3. CodonDecision captures cai_impact
4. CodonDecision serialization/deserialization with cai_impact
5. ConflictProvenance captures codon_changes
6. record_relaxation_provenance accepts cai_impact_override and codon_changes
7. record_violation_provenance accepts cai_impact_override and codon_changes
8. record_cai_aware_provenance populates codon_changes
9. Total CAI cost is computed correctly from codon decisions
10. Total CAI cost is computed correctly from constraint decisions
11. ProvenanceReport includes CAI impact summary
12. ProvenanceReport includes top 5 most CAI-expensive fixes
13. ProvenanceReport JSON export includes cai_impact_summary
14. ProvenanceTracker records decisions with CAI impact
15. ProvenanceTracker serialization preserves CAI fields
"""

import json
from datetime import datetime, timezone

import pytest

from biocompiler.provenance import (
    DecisionRecord,
    ProvenanceTracker,
    OptimizationProvenance,
    OptimizationRecord,
    generate_provenance_report,
)
from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
    DecisionProvenanceCollector,
)
from biocompiler.provenance_reporting import (
    ProvenanceQuery,
    ProvenanceReport,
    explain_position,
)
from biocompiler.solver.conflict_provenance import (
    ConflictProvenance,
    ConflictResolverWithProvenance,
)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. DecisionRecord CAI fields
# ═══════════════════════════════════════════════════════════════════════════════

class TestDecisionRecordCAIFields:
    """Test that DecisionRecord captures cai_impact, codon_before, codon_after."""

    def test_default_cai_impact_is_zero(self):
        """Default cai_impact should be 0.0."""
        record = DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Test",
            constraint_context={},
        )
        assert record.cai_impact == 0.0

    def test_default_codon_before_is_empty(self):
        """Default codon_before should be empty string."""
        record = DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Test",
            constraint_context={},
        )
        assert record.codon_before == ""

    def test_default_codon_after_is_empty(self):
        """Default codon_after should be empty string."""
        record = DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Test",
            constraint_context={},
        )
        assert record.codon_after == ""

    def test_cai_impact_negative_means_loss(self):
        """Negative cai_impact should represent CAI sacrificed."""
        record = DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="constraint_relaxed",
            position=5,
            chosen_value="GTC",
            alternatives_considered=["GTG"],
            rationale="Chose GTC to avoid restriction site",
            constraint_context={"cai": 0.85},
            cai_impact=-0.003,
            codon_before="GTG",
            codon_after="GTC",
        )
        assert record.cai_impact == -0.003
        assert record.codon_before == "GTG"
        assert record.codon_after == "GTC"

    def test_cai_impact_positive_means_gain(self):
        """Positive cai_impact should represent CAI gained."""
        record = DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="codon_selected",
            position=2,
            chosen_value="CTG",
            alternatives_considered=["TTA"],
            rationale="High-CAI codon for Leu",
            constraint_context={"cai": 0.95},
            cai_impact=0.02,
            codon_before="TTA",
            codon_after="CTG",
        )
        assert record.cai_impact == 0.02


class TestDecisionRecordCAISerialization:
    """Test that DecisionRecord serialization includes CAI fields."""

    def _make_record(self, **kwargs):
        defaults = dict(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Test",
            constraint_context={},
        )
        defaults.update(kwargs)
        return DecisionRecord(**defaults)

    def test_to_dict_includes_cai_fields(self):
        """to_dict should include cai_impact, codon_before, codon_after."""
        record = self._make_record(
            cai_impact=-0.005,
            codon_before="GTG",
            codon_after="GTC",
        )
        d = record.to_dict()
        assert "cai_impact" in d
        assert d["cai_impact"] == -0.005
        assert "codon_before" in d
        assert d["codon_before"] == "GTG"
        assert "codon_after" in d
        assert d["codon_after"] == "GTC"

    def test_from_dict_with_cai_fields(self):
        """from_dict should restore cai_impact, codon_before, codon_after."""
        d = {
            "timestamp": "2026-01-01T00:00:00Z",
            "decision_type": "codon_selected",
            "position": 3,
            "chosen_value": "CTG",
            "alternatives_considered": ["TTA"],
            "rationale": "CAI optimization",
            "constraint_context": {"cai": 0.90},
            "cai_impact": -0.012,
            "codon_before": "TTA",
            "codon_after": "CTG",
        }
        record = DecisionRecord.from_dict(d)
        assert record.cai_impact == -0.012
        assert record.codon_before == "TTA"
        assert record.codon_after == "CTG"

    def test_from_dict_without_cai_fields_uses_defaults(self):
        """from_dict should use defaults when CAI fields are missing (backward compat)."""
        d = {
            "timestamp": "2026-01-01T00:00:00Z",
            "decision_type": "codon_selected",
            "position": 0,
            "chosen_value": "ATG",
            "alternatives_considered": [],
            "rationale": "Test",
            "constraint_context": {},
        }
        record = DecisionRecord.from_dict(d)
        assert record.cai_impact == 0.0
        assert record.codon_before == ""
        assert record.codon_after == ""

    def test_round_trip_serialization(self):
        """to_dict → from_dict should preserve all CAI fields."""
        original = self._make_record(
            cai_impact=-0.007,
            codon_before="GTA",
            codon_after="GTC",
        )
        restored = DecisionRecord.from_dict(original.to_dict())
        assert restored.cai_impact == original.cai_impact
        assert restored.codon_before == original.codon_before
        assert restored.codon_after == original.codon_after


# ═══════════════════════════════════════════════════════════════════════════════
# 2. CodonDecision cai_impact
# ═══════════════════════════════════════════════════════════════════════════════

class TestCodonDecisionCAIImpact:
    """Test that CodonDecision captures cai_impact."""

    def test_default_cai_impact_is_zero(self):
        """Default cai_impact should be 0.0."""
        cd = CodonDecision(
            position=0,
            amino_acid="M",
            original_codon=None,
            chosen_codon="ATG",
            alternatives_considered=[],
            constraint_reason="maximize_cai",
            confidence=1.0,
        )
        assert cd.cai_impact == 0.0

    def test_cai_impact_for_constrained_choice(self):
        """cai_impact should be set when a codon is chosen due to a constraint."""
        cd = CodonDecision(
            position=5,
            amino_acid="V",
            original_codon=None,
            chosen_codon="GTC",
            alternatives_considered=[
                {"codon": "GTG", "cai_contribution": 1.0, "gc_contribution": 1.0},
            ],
            constraint_reason="gc_content",
            confidence=0.8,
            cai_impact=-0.015,
        )
        assert cd.cai_impact == -0.015

    def test_to_dict_includes_cai_impact(self):
        """to_dict should include cai_impact."""
        cd = CodonDecision(
            position=3,
            amino_acid="L",
            original_codon="TTA",
            chosen_codon="CTG",
            alternatives_considered=[],
            constraint_reason="maximize_cai",
            confidence=1.0,
            cai_impact=0.0,
        )
        d = cd.to_dict()
        assert "cai_impact" in d
        assert d["cai_impact"] == 0.0

    def test_from_dict_restores_cai_impact(self):
        """from_dict should restore cai_impact."""
        d = {
            "position": 5,
            "amino_acid": "V",
            "original_codon": None,
            "chosen_codon": "GTC",
            "alternatives_considered": [],
            "constraint_reason": "gc_content",
            "confidence": 0.8,
            "cai_impact": -0.025,
        }
        cd = CodonDecision.from_dict(d)
        assert cd.cai_impact == -0.025

    def test_from_dict_without_cai_impact_defaults_zero(self):
        """from_dict should default cai_impact to 0.0 if missing."""
        d = {
            "position": 0,
            "amino_acid": "M",
            "original_codon": None,
            "chosen_codon": "ATG",
            "alternatives_considered": [],
            "constraint_reason": "maximize_cai",
            "confidence": 1.0,
        }
        cd = CodonDecision.from_dict(d)
        assert cd.cai_impact == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 3. ConflictProvenance codon_changes
# ═══════════════════════════════════════════════════════════════════════════════

class TestConflictProvenanceCodonChange:
    """Test that ConflictProvenance captures codon_changes."""

    def test_default_codon_changes_is_empty_list(self):
        """Default codon_changes should be an empty list."""
        cp = ConflictProvenance(
            conflicting_constraints=["A", "B"],
            resolution_method="priority_based",
            winner="A",
            loser="B",
            impact="A won",
            positions_affected=[1, 2],
        )
        assert cp.codon_changes == []

    def test_codon_changes_with_single_change(self):
        """codon_changes should capture a single codon substitution."""
        cp = ConflictProvenance(
            conflicting_constraints=["NoRestrictionSite", "MaximizeCAI"],
            resolution_method="cai_aware",
            winner="NoRestrictionSite",
            loser="MaximizeCAI",
            impact="Fixed",
            positions_affected=[5],
            cai_impact=-0.003,
            cai_delta=-0.003,
            codon_changes=[(5, "GTG", "GTC", -0.003)],
        )
        assert len(cp.codon_changes) == 1
        assert cp.codon_changes[0] == (5, "GTG", "GTC", -0.003)

    def test_codon_changes_with_multiple_changes(self):
        """codon_changes should capture multiple codon substitutions."""
        cp = ConflictProvenance(
            conflicting_constraints=["GCInRange", "MaximizeCAI"],
            resolution_method="cai_aware",
            winner="GCInRange",
            loser="MaximizeCAI",
            impact="Fixed GC",
            positions_affected=[3, 7],
            cai_impact=-0.01,
            cai_delta=-0.01,
            codon_changes=[
                (3, "CTC", "CTG", -0.004),
                (7, "GCG", "GCT", -0.006),
            ],
        )
        assert len(cp.codon_changes) == 2
        assert cp.codon_changes[0] == (3, "CTC", "CTG", -0.004)
        assert cp.codon_changes[1] == (7, "GCG", "GCT", -0.006)

    def test_repr_includes_codon_changes_count(self):
        """repr should include codon_changes count when non-empty."""
        cp = ConflictProvenance(
            conflicting_constraints=["A", "B"],
            resolution_method="cai_aware",
            winner="A",
            loser="B",
            impact="Fixed",
            positions_affected=[1],
            codon_changes=[(1, "AAA", "AAG", -0.001)],
        )
        r = repr(cp)
        assert "codon_changes=1" in r

    def test_repr_excludes_codon_changes_when_empty(self):
        """repr should not include codon_changes when empty."""
        cp = ConflictProvenance(
            conflicting_constraints=["A", "B"],
            resolution_method="priority_based",
            winner="A",
            loser="B",
            impact="A won",
            positions_affected=[1],
        )
        r = repr(cp)
        assert "codon_changes" not in r


# ═══════════════════════════════════════════════════════════════════════════════
# 4. ConflictResolverWithProvenance CAI methods
# ═══════════════════════════════════════════════════════════════════════════════

class TestConflictResolverCAIMethods:
    """Test that record_relaxation_provenance and record_violation_provenance
    accept and store CAI impact data."""

    def test_record_cai_aware_provenance_includes_codon_changes(self):
        """record_cai_aware_provenance should populate codon_changes."""
        resolver = ConflictResolverWithProvenance(
            track_provenance=True, organism="Homo_sapiens"
        )
        record = resolver.record_cai_aware_provenance(
            constraint_name="NoRestrictionSite",
            codon_position=5,
            old_codon="GTG",
            new_codon="GTC",
            cai_delta=-0.003,
        )
        assert len(record.codon_changes) == 1
        assert record.codon_changes[0] == (5, "GTG", "GTC", -0.003)
        assert record.cai_delta == -0.003
        assert record.resolution_method == "cai_aware"

    def test_record_relaxation_with_cai_override(self):
        """record_relaxation_provenance should use cai_impact_override when provided."""
        resolver = ConflictResolverWithProvenance(
            track_provenance=True, organism="Homo_sapiens"
        )
        record = resolver.record_relaxation_provenance(
            relaxed_constraint_name="MaximizeCAI",
            kept_constraint_name="GCInRange",
            positions_affected=[3, 7],
            sequence="ATGGTGCTG",
            cai_impact_override=-0.008,
            codon_changes=[
                (3, "CTC", "CTG", -0.004),
                (7, "GCG", "GCT", -0.004),
            ],
        )
        assert record.cai_impact == -0.008
        assert record.cai_delta == -0.008
        assert len(record.codon_changes) == 2

    def test_record_relaxation_without_cai_override_uses_heuristic(self):
        """record_relaxation_provenance should use heuristic when no override."""
        resolver = ConflictResolverWithProvenance(
            track_provenance=True, organism="Homo_sapiens"
        )
        record = resolver.record_relaxation_provenance(
            relaxed_constraint_name="GCInRange",
            kept_constraint_name="NoRestrictionSite",
            positions_affected=[5],
            sequence="ATGGTG",
        )
        # Should use heuristic estimate (GC constraint has -0.05)
        assert record.cai_impact == -0.05
        assert record.cai_delta is None
        assert record.codon_changes == []

    def test_record_violation_with_cai_override(self):
        """record_violation_provenance should use cai_impact_override when provided."""
        resolver = ConflictResolverWithProvenance(
            track_provenance=True, organism="Homo_sapiens"
        )

        # Create a minimal mock violation object
        class MockViolation:
            constraint_name = "GCInRange"
            severity = 0.5
            positions = [3, 7]
            class priority:
                name = "MEDIUM"

        violations = [MockViolation()]
        records = resolver.record_violation_provenance(
            violations=violations,
            sequence="ATGGTGCTG",
            cai_impact_override=-0.015,
            codon_changes=[(3, "TTA", "CTG", -0.008), (7, "GCG", "GCT", -0.007)],
        )
        assert len(records) == 1
        assert records[0].cai_impact == -0.015
        assert records[0].cai_delta == -0.015
        assert len(records[0].codon_changes) == 2

    def test_record_violation_without_cai_override(self):
        """record_violation_provenance should use heuristic when no override."""
        resolver = ConflictResolverWithProvenance(
            track_provenance=True, organism="Homo_sapiens"
        )

        class MockViolation:
            constraint_name = "NoCrypticSplice"
            severity = 0.3
            positions = [12]
            class priority:
                name = "HIGH"

        violations = [MockViolation()]
        records = resolver.record_violation_provenance(
            violations=violations,
            sequence="ATGGTGCTGAAAGGG",
        )
        assert len(records) == 1
        # Heuristic for splice constraint is 0.01
        assert records[0].cai_impact == 0.01
        assert records[0].codon_changes == []


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Total CAI cost computation
# ═══════════════════════════════════════════════════════════════════════════════

class TestTotalCAICostComputation:
    """Test that total CAI cost is computed correctly from decisions."""

    def test_cai_cost_from_codon_decisions(self):
        """Total CAI cost should sum all negative cai_impact from CodonDecisions."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MVV",
            organism="Homo_sapiens",
            constraints=["GCInRange"],
        )
        collector.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
            cai_impact=0.0,
        ))
        collector.record_codon_decision(CodonDecision(
            position=1, amino_acid="V", original_codon=None,
            chosen_codon="GTC", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.8,
            cai_impact=-0.005,
        ))
        collector.record_codon_decision(CodonDecision(
            position=2, amino_acid="V", original_codon=None,
            chosen_codon="GTG", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.7,
            cai_impact=-0.008,
        ))
        trail = collector.finalize(output_dna="ATGGTCGTG", cai=0.88, gc=0.56)

        total_cai_cost = sum(
            cd.cai_impact for cd in trail.codon_decisions if cd.cai_impact < 0
        )
        assert total_cai_cost == pytest.approx(-0.013)

    def test_cai_cost_from_constraint_decisions(self):
        """Total CAI cost should sum impact_on_cai from ConstraintDecisions."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MVL",
            organism="Homo_sapiens",
            constraints=["GCInRange", "NoCrypticSplice"],
        )
        collector.record_constraint_decision(ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[1, 2],
            tradeoff_description="GC kept in range",
            impact_on_cai=-0.005,
        ))
        collector.record_constraint_decision(ConstraintDecision(
            constraint_name="NoCrypticSplice",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[3],
            tradeoff_description="Avoided cryptic donor",
            impact_on_cai=-0.002,
        ))
        trail = collector.finalize(output_dna="ATGGTGCTG", cai=0.90, gc=0.55)

        total_cai_cost = sum(
            cd.impact_on_cai for cd in trail.constraint_decisions
        )
        assert total_cai_cost == pytest.approx(-0.007)

    def test_no_cai_cost_when_unconstrained(self):
        """Total CAI cost should be 0 when all cai_impact values are 0."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="M",
            organism="Escherichia_coli",
            constraints=[],
        )
        collector.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
            cai_impact=0.0,
        ))
        trail = collector.finalize(output_dna="ATG", cai=1.0, gc=0.0)

        total_cai_cost = sum(
            cd.cai_impact for cd in trail.codon_decisions if cd.cai_impact < 0
        )
        assert total_cai_cost == 0.0

    def test_provenance_query_cai_lost_to_constraints(self):
        """ProvenanceQuery.cai_lost_to_constraints should return total CAI loss."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MV",
            organism="Homo_sapiens",
            constraints=["GCInRange"],
        )
        collector.record_constraint_decision(ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[1],
            tradeoff_description="GC fixed",
            impact_on_cai=-0.01,
        ))
        collector.record_constraint_decision(ConstraintDecision(
            constraint_name="NoRestrictionSite_EcoRI",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[2],
            tradeoff_description="Avoided EcoRI",
            impact_on_cai=-0.003,
        ))
        trail = collector.finalize(output_dna="ATGGTC", cai=0.88, gc=0.50)
        query = ProvenanceQuery(trail)

        cai_lost = query.cai_lost_to_constraints()
        assert cai_lost == pytest.approx(-0.013)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Provenance reporting includes CAI impact summary
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceReportCAIImpact:
    """Test that ProvenanceReport includes CAI impact summary."""

    def _make_trail_with_cai_cost(self):
        """Create a trail with some CAI-costly codon decisions."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MVVLL",
            organism="Homo_sapiens",
            constraints=["GCInRange", "NoCrypticSplice"],
        )
        # M - no cost
        collector.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
            cai_impact=0.0,
        ))
        # V - small cost
        collector.record_codon_decision(CodonDecision(
            position=1, amino_acid="V", original_codon=None,
            chosen_codon="GTC", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.8,
            cai_impact=-0.003,
        ))
        # V - larger cost
        collector.record_codon_decision(CodonDecision(
            position=2, amino_acid="V", original_codon=None,
            chosen_codon="GTA", alternatives_considered=[],
            constraint_reason="avoid_cpg", confidence=0.6,
            cai_impact=-0.012,
        ))
        # L - medium cost
        collector.record_codon_decision(CodonDecision(
            position=3, amino_acid="L", original_codon=None,
            chosen_codon="CTC", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.7,
            cai_impact=-0.006,
        ))
        # L - tiny cost
        collector.record_codon_decision(CodonDecision(
            position=4, amino_acid="L", original_codon=None,
            chosen_codon="CTT", alternatives_considered=[],
            constraint_reason="no_cryptic_splice", confidence=0.9,
            cai_impact=-0.001,
        ))
        collector.record_constraint_decision(ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[1, 3],
            tradeoff_description="GC kept in range",
            impact_on_cai=-0.009,
        ))
        return collector.finalize(
            output_dna="ATGGTCGTACTCCTT", cai=0.85, gc=0.50
        )

    def test_markdown_report_includes_cai_cost_summary(self):
        """Markdown report should include 'Total CAI cost of constraint satisfaction'."""
        trail = self._make_trail_with_cai_cost()
        report = ProvenanceReport.generate_markdown(trail)
        assert "Total CAI cost of constraint satisfaction" in report

    def test_markdown_report_includes_top5_expensive_fixes(self):
        """Markdown report should include 'Top 5 Most CAI-Expensive Constraint Fixes'."""
        trail = self._make_trail_with_cai_cost()
        report = ProvenanceReport.generate_markdown(trail)
        assert "Top 5 Most CAI-Expensive Constraint Fixes" in report

    def test_markdown_report_top5_lists_positions(self):
        """Top 5 section should list position numbers."""
        trail = self._make_trail_with_cai_cost()
        report = ProvenanceReport.generate_markdown(trail)
        # The most expensive position is 2 (cai_impact=-0.012)
        assert "2" in report

    def test_html_report_includes_cai_cost(self):
        """HTML report should include CAI cost of constraint satisfaction."""
        trail = self._make_trail_with_cai_cost()
        report = ProvenanceReport.generate_html(trail)
        assert "CAI Cost of Constraint Satisfaction" in report

    def test_html_report_includes_top5_expensive_fixes(self):
        """HTML report should include Top 5 most CAI-expensive fixes."""
        trail = self._make_trail_with_cai_cost()
        report = ProvenanceReport.generate_html(trail)
        assert "Top 5 Most CAI-Expensive Constraint Fixes" in report

    def test_json_report_includes_cai_impact_summary(self):
        """JSON report should include cai_impact_summary in analysis."""
        trail = self._make_trail_with_cai_cost()
        report_str = ProvenanceReport.generate_json(trail)
        report_data = json.loads(report_str)

        assert "analysis" in report_data
        assert "cai_impact_summary" in report_data["analysis"]
        summary = report_data["analysis"]["cai_impact_summary"]

        assert "total_cai_cost_of_constraint_satisfaction" in summary
        assert "positions_with_cai_cost" in summary
        assert "top_5_most_expensive_fixes" in summary

        # Total cost should be sum of negative cai_impacts: -0.003 + -0.012 + -0.006 + -0.001
        assert summary["total_cai_cost_of_constraint_satisfaction"] == pytest.approx(-0.022)
        assert summary["positions_with_cai_cost"] == 4
        assert len(summary["top_5_most_expensive_fixes"]) == 4  # only 4 costly positions

    def test_json_report_top5_ordered_by_cost(self):
        """Top 5 most expensive fixes should be ordered by CAI impact (most negative first)."""
        trail = self._make_trail_with_cai_cost()
        report_str = ProvenanceReport.generate_json(trail)
        report_data = json.loads(report_str)
        fixes = report_data["analysis"]["cai_impact_summary"]["top_5_most_expensive_fixes"]

        # Most expensive first (most negative)
        assert fixes[0]["position"] == 2
        assert fixes[0]["cai_impact"] == pytest.approx(-0.012)
        assert fixes[1]["position"] == 3
        assert fixes[1]["cai_impact"] == pytest.approx(-0.006)

    def test_report_with_no_cai_cost(self):
        """Report should handle case where no positions incur CAI cost."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="M",
            organism="Escherichia_coli",
            constraints=[],
        )
        collector.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
            cai_impact=0.0,
        ))
        trail = collector.finalize(output_dna="ATG", cai=1.0, gc=0.0)

        md_report = ProvenanceReport.generate_markdown(trail)
        assert "Total CAI cost of constraint satisfaction" in md_report
        assert "0.0000" in md_report

        html_report = ProvenanceReport.generate_html(trail)
        assert "CAI Cost of Constraint Satisfaction" in html_report

        json_str = ProvenanceReport.generate_json(trail)
        json_data = json.loads(json_str)
        summary = json_data["analysis"]["cai_impact_summary"]
        assert summary["total_cai_cost_of_constraint_satisfaction"] == 0.0
        assert summary["positions_with_cai_cost"] == 0
        assert summary["top_5_most_expensive_fixes"] == []


# ═══════════════════════════════════════════════════════════════════════════════
# 7. ProvenanceTracker with CAI impact
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceTrackerCAIImpact:
    """Test that ProvenanceTracker records decisions with CAI impact and
    serialization preserves CAI fields."""

    def test_tracker_records_cai_impact(self):
        """ProvenanceTracker should store decisions with cai_impact."""
        tracker = ProvenanceTracker(seed=42)
        tracker.record_decision(DecisionRecord(
            timestamp=datetime.now(timezone.utc).isoformat(),
            decision_type="codon_selected",
            position=5,
            chosen_value="GTC",
            alternatives_considered=["GTG"],
            rationale="Chose GTC to avoid EcoRI site",
            constraint_context={"cai": 0.85},
            cai_impact=-0.003,
            codon_before="GTG",
            codon_after="GTC",
        ))

        decisions = tracker.get_decisions_for_position(5)
        assert len(decisions) == 1
        assert decisions[0].cai_impact == -0.003
        assert decisions[0].codon_before == "GTG"
        assert decisions[0].codon_after == "GTC"

    def test_tracker_serialization_preserves_cai_fields(self):
        """ProvenanceTracker.to_dict / from_dict should preserve CAI fields."""
        tracker = ProvenanceTracker(seed=42)
        tracker.record_decision(DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="constraint_relaxed",
            position=3,
            chosen_value="CTG",
            alternatives_considered=["CTC"],
            rationale="GC constraint",
            constraint_context={"gc": 0.55},
            cai_impact=-0.007,
            codon_before="CTC",
            codon_after="CTG",
        ))

        # Serialize
        data = tracker.to_dict()
        # Check CAI fields are in serialized decision
        assert data["decisions"][0]["cai_impact"] == -0.007
        assert data["decisions"][0]["codon_before"] == "CTC"
        assert data["decisions"][0]["codon_after"] == "CTG"

        # Deserialize
        restored = ProvenanceTracker.from_dict(data)
        decisions = restored.get_decisions_for_position(3)
        assert len(decisions) == 1
        assert decisions[0].cai_impact == -0.007
        assert decisions[0].codon_before == "CTC"
        assert decisions[0].codon_after == "CTG"

    def test_tracker_json_round_trip(self):
        """ProvenanceTracker.to_json / from_json should preserve CAI fields."""
        tracker = ProvenanceTracker(seed=42)
        tracker.record_decision(DecisionRecord(
            timestamp="2026-01-01T00:00:00Z",
            decision_type="codon_selected",
            position=7,
            chosen_value="GCT",
            alternatives_considered=["GCC"],
            rationale="Avoid CpG",
            constraint_context={"cai": 0.80},
            cai_impact=-0.004,
            codon_before="GCC",
            codon_after="GCT",
        ))

        json_str = tracker.to_json()
        restored = ProvenanceTracker.from_json(json_str)
        decisions = restored.get_decisions_for_position(7)
        assert decisions[0].cai_impact == -0.004
        assert decisions[0].codon_before == "GCC"
        assert decisions[0].codon_after == "GCT"


# ═══════════════════════════════════════════════════════════════════════════════
# 8. DecisionProvenanceCollector summary with CAI cost
# ═══════════════════════════════════════════════════════════════════════════════

class TestCollectorSummaryWithCAICost:
    """Test that DecisionProvenanceCollector.summary() includes CAI cost info."""

    def test_summary_includes_total_cai_cost(self):
        """summary() should include total CAI cost from codon decisions."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MV",
            organism="Homo_sapiens",
            constraints=["GCInRange"],
        )
        collector.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
            cai_impact=0.0,
        ))
        collector.record_codon_decision(CodonDecision(
            position=1, amino_acid="V", original_codon=None,
            chosen_codon="GTC", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.8,
            cai_impact=-0.005,
        ))
        summary = collector.summary()
        assert "Total CAI cost (codon decisions)" in summary
        assert "-0.0050" in summary

    def test_summary_includes_constrained_position_count(self):
        """summary() should include count of constrained codon positions."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="MVV",
            organism="Homo_sapiens",
            constraints=["GCInRange"],
        )
        collector.record_codon_decision(CodonDecision(
            position=0, amino_acid="M", original_codon=None,
            chosen_codon="ATG", alternatives_considered=[],
            constraint_reason="maximize_cai", confidence=1.0,
            cai_impact=0.0,
        ))
        collector.record_codon_decision(CodonDecision(
            position=1, amino_acid="V", original_codon=None,
            chosen_codon="GTC", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.8,
            cai_impact=-0.005,
        ))
        collector.record_codon_decision(CodonDecision(
            position=2, amino_acid="V", original_codon=None,
            chosen_codon="GTG", alternatives_considered=[],
            constraint_reason="gc_content", confidence=0.7,
            cai_impact=-0.008,
        ))
        summary = collector.summary()
        assert "Constrained codon positions" in summary
        assert "2" in summary
