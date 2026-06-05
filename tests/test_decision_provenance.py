"""
Integration tests for decision-level provenance pipeline.

Covers:
1. Basic codon decision recording with provenance tracking ON
2. Alternatives_considered quality (Leucine, Methionine)
3. Constraint decision recording with tight GC constraint
4. What-if analysis (GC relaxation)
5. Provenance query (constraints_that_reduced_cai, explain_position)
6. Backward compatibility (track_provenance=False)
7. Report generation (markdown report from decision trail)
"""

import pytest

from biocompiler.optimization import optimize_sequence
from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
    DecisionProvenanceCollector,
)
from biocompiler.provenance_reporting import ProvenanceQuery, ProvenanceReport, explain_position
from biocompiler.whatif_analysis import WhatIfAnalyzer, WhatIfScenario


# ────────────────────────────────────────────────────────────
# Shared test data
# ────────────────────────────────────────────────────────────

# Short 20aa protein containing Leu (L) and Met (M)
SHORT_PROTEIN = "MALWMRLLPLLALLALWGP"

# Insulin B chain (proinsulin-derived, 30aa, for what-if test)
INSULIN_PROTEIN = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Basic codon decision recording
# ═══════════════════════════════════════════════════════════════════════════════

class TestBasicCodonDecisionRecording:
    """Optimize a short protein with provenance tracking ON and verify the
    OptimizationResult has a populated decision_trail."""

    @pytest.fixture
    def result(self):
        """Optimize a 20aa protein with provenance tracking enabled."""
        return optimize_sequence(
            target_protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            track_provenance=True,
        )

    def test_result_has_decision_trail(self, result):
        """OptimizationResult should have decision_trail when tracking is ON."""
        assert result.decision_trail is not None, (
            "decision_trail should be populated when track_provenance=True"
        )

    def test_decision_trail_is_optimization_decision_trail(self, result):
        """decision_trail should be an OptimizationDecisionTrail instance."""
        assert isinstance(result.decision_trail, OptimizationDecisionTrail)

    def test_codon_decisions_present(self, result):
        """decision_trail should have codon_decisions."""
        trail = result.decision_trail
        assert len(trail.codon_decisions) > 0, (
            "decision_trail should contain codon decisions"
        )

    def test_codon_decisions_cover_protein_positions(self, result):
        """Codon decisions should cover all positions in the protein."""
        trail = result.decision_trail
        # Get the set of unique positions from codon decisions
        positions = set(cd.position for cd in trail.codon_decisions)
        expected_positions = set(range(len(SHORT_PROTEIN)))
        assert positions == expected_positions, (
            f"Expected positions {expected_positions}, got {positions}. "
            f"Missing: {expected_positions - positions}"
        )

    def test_each_codon_decision_has_required_fields(self, result):
        """Each CodonDecision should have chosen_codon, alternatives_considered,
        and constraint_reason."""
        trail = result.decision_trail
        for cd in trail.codon_decisions:
            assert isinstance(cd, CodonDecision)
            assert cd.chosen_codon, "chosen_codon should not be empty"
            assert isinstance(cd.alternatives_considered, list), (
                "alternatives_considered should be a list"
            )
            assert isinstance(cd.constraint_reason, str), (
                "constraint_reason should be a string"
            )
            assert cd.constraint_reason, "constraint_reason should not be empty"

    def test_amino_acids_match_input(self, result):
        """Each CodonDecision.amino_acid should match the input protein."""
        trail = result.decision_trail
        # Build a map of last decision per position
        last_per_pos = {}
        for cd in trail.codon_decisions:
            last_per_pos[cd.position] = cd
        for pos, cd in last_per_pos.items():
            if pos < len(SHORT_PROTEIN):
                assert cd.amino_acid == SHORT_PROTEIN[pos], (
                    f"Position {pos}: expected amino acid {SHORT_PROTEIN[pos]}, "
                    f"got {cd.amino_acid}"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Alternatives_considered quality
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlternativesConsideredQuality:
    """Verify the quality of alternatives_considered for specific amino acids."""

    @pytest.fixture
    def result(self):
        return optimize_sequence(
            target_protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            track_provenance=True,
        )

    def test_leucine_has_at_least_3_alternatives(self, result):
        """Leucine (L) has 6 codons; at least 3 alternatives should be considered."""
        trail = result.decision_trail
        query = ProvenanceQuery(trail)
        leucine_decisions = query.decisions_for_amino_acid("L")
        assert len(leucine_decisions) > 0, "Protein should contain Leucine"
        for cd in leucine_decisions:
            # alternatives_considered includes the chosen codon itself,
            # so count alternatives that are NOT the chosen one
            other_alts = [
                a for a in cd.alternatives_considered
                if isinstance(a, dict) and a.get("codon") != cd.chosen_codon
            ]
            assert len(other_alts) >= 3, (
                f"Leucine at position {cd.position} should have at least 3 "
                f"alternatives (excluding chosen), got {len(other_alts)}. "
                f"Total alternatives_considered: {len(cd.alternatives_considered)}"
            )

    def test_methionine_has_no_other_alternatives(self, result):
        """Methionine (M) has only ATG; no other alternatives should be listed."""
        trail = result.decision_trail
        query = ProvenanceQuery(trail)
        met_decisions = query.decisions_for_amino_acid("M")
        assert len(met_decisions) > 0, "Protein should contain Methionine"
        for cd in met_decisions:
            # Count alternatives that are different from ATG
            other_alts = [
                a for a in cd.alternatives_considered
                if isinstance(a, dict) and a.get("codon") != "ATG"
            ]
            assert len(other_alts) == 0, (
                f"Methionine at position {cd.position} should have no "
                f"alternatives besides ATG, got {len(other_alts)} others"
            )
            assert cd.chosen_codon == "ATG", (
                f"Methionine codon should be ATG, got {cd.chosen_codon}"
            )

    def test_alternatives_have_cai_contribution(self, result):
        """Each alternative dict should have cai_contribution and gc_contribution."""
        trail = result.decision_trail
        for cd in trail.codon_decisions:
            for alt in cd.alternatives_considered:
                if isinstance(alt, dict):
                    assert "cai_contribution" in alt, (
                        f"Alternative for position {cd.position} missing cai_contribution. "
                        f"Keys: {list(alt.keys())}"
                    )
                    assert "gc_contribution" in alt, (
                        f"Alternative for position {cd.position} missing gc_contribution. "
                        f"Keys: {list(alt.keys())}"
                    )
                    # cai_contribution should be a number
                    assert isinstance(alt["cai_contribution"], (int, float)), (
                        f"cai_contribution should be numeric, got {type(alt['cai_contribution'])}"
                    )
                    assert isinstance(alt["gc_contribution"], (int, float)), (
                        f"gc_contribution should be numeric, got {type(alt['gc_contribution'])}"
                    )

    def test_alternatives_have_codon_field(self, result):
        """Each alternative dict should have a 'codon' field identifying it."""
        trail = result.decision_trail
        for cd in trail.codon_decisions:
            for alt in cd.alternatives_considered:
                if isinstance(alt, dict):
                    assert "codon" in alt, (
                        f"Alternative for position {cd.position} missing 'codon' field. "
                        f"Keys: {list(alt.keys())}"
                    )
                    assert isinstance(alt["codon"], str), (
                        "codon field should be a string"
                    )
                    assert len(alt["codon"]) == 3, (
                        f"Codon should be 3 bases, got '{alt['codon']}'"
                    )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Constraint decision recording
# ═══════════════════════════════════════════════════════════════════════════════

class TestConstraintDecisionRecording:
    """Optimize with a tight GC constraint and verify constraint_decisions."""

    @pytest.fixture
    def result(self):
        """Optimize with a tight GC range that may force sub-optimal codon choices."""
        return optimize_sequence(
            target_protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.45,
            gc_hi=0.55,
            track_provenance=True,
        )

    def test_constraint_decisions_non_empty(self, result):
        """constraint_decisions should be non-empty after optimization."""
        trail = result.decision_trail
        assert trail is not None
        assert len(trail.constraint_decisions) > 0, (
            "constraint_decisions should not be empty after optimization"
        )

    def test_each_constraint_decision_has_required_fields(self, result):
        """Each ConstraintDecision should have action_taken and tradeoff_description."""
        trail = result.decision_trail
        for cd in trail.constraint_decisions:
            assert isinstance(cd, ConstraintDecision)
            assert cd.action_taken in ("satisfied", "relaxed", "conflicted", "overridden"), (
                f"Invalid action_taken: {cd.action_taken}"
            )
            assert isinstance(cd.tradeoff_description, str), (
                "tradeoff_description should be a string"
            )
            assert cd.tradeoff_description, (
                "tradeoff_description should not be empty"
            )

    def test_constraint_decisions_have_names(self, result):
        """Each ConstraintDecision should have a meaningful constraint_name."""
        trail = result.decision_trail
        for cd in trail.constraint_decisions:
            assert cd.constraint_name, (
                "constraint_name should not be empty"
            )

    def test_gc_constraint_is_recorded(self, result):
        """When GC constraint is active, GCInRange should appear in constraints."""
        trail = result.decision_trail
        constraint_names = [cd.constraint_name for cd in trail.constraint_decisions]
        assert "GCInRange" in constraint_names, (
            f"GCInRange should be in constraint decisions, got {constraint_names}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 4. What-if analysis
# ═══════════════════════════════════════════════════════════════════════════════

class TestWhatIfAnalysis:
    """Test what-if analysis: relaxing GC constraint for insulin in E. coli."""

    @pytest.fixture
    def baseline_result(self):
        """Optimize insulin for E. coli with GC [0.30, 0.70]."""
        return optimize_sequence(
            target_protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            track_provenance=True,
        )

    def test_gc_relaxation_scenario(self, baseline_result):
        """Relax GC to [0.30, 0.80] and verify scenario is produced."""
        analyzer = WhatIfAnalyzer(
            protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        scenario = analyzer.analyze_gc_relaxation(
            dna=baseline_result.sequence,
            current_gc_hi=0.70,
            alternative_gc_hi=0.80,
        )
        assert isinstance(scenario, WhatIfScenario)

    def test_gc_relaxation_is_feasible(self, baseline_result):
        """Relaxing a constraint should not make the problem infeasible."""
        analyzer = WhatIfAnalyzer(
            protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        scenario = analyzer.analyze_gc_relaxation(
            dna=baseline_result.sequence,
            current_gc_hi=0.70,
            alternative_gc_hi=0.80,
        )
        assert scenario.feasibility == "feasible", (
            f"Relaxing GC constraint should be feasible, got {scenario.feasibility}"
        )

    def test_gc_relaxation_predicts_cai_change(self, baseline_result):
        """The scenario should predict or measure a CAI change."""
        analyzer = WhatIfAnalyzer(
            protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        scenario = analyzer.analyze_gc_relaxation(
            dna=baseline_result.sequence,
            current_gc_hi=0.70,
            alternative_gc_hi=0.80,
        )
        if scenario.feasibility == "feasible":
            assert scenario.predicted_cai is not None, (
                "Feasible scenario should have predicted_cai"
            )
            assert 0.0 <= scenario.predicted_cai <= 1.0, (
                f"predicted_cai should be in [0,1], got {scenario.predicted_cai}"
            )
            # cai_delta should be computed
            assert scenario.cai_delta is not None, (
                "cai_delta should be computed when predicted_cai is available"
            )

    def test_scenario_has_baseline_cai(self, baseline_result):
        """The scenario should record the baseline CAI."""
        analyzer = WhatIfAnalyzer(
            protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        scenario = analyzer.analyze_gc_relaxation(
            dna=baseline_result.sequence,
            current_gc_hi=0.70,
            alternative_gc_hi=0.80,
        )
        assert scenario.baseline_cai >= 0.0, (
            "baseline_cai should be non-negative"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Provenance query
# ═══════════════════════════════════════════════════════════════════════════════

class TestProvenanceQuery:
    """Test ProvenanceQuery for constraints_that_reduced_cai and explain_position."""

    @pytest.fixture
    def trail(self):
        """Get a decision trail from an optimization run."""
        result = optimize_sequence(
            target_protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            track_provenance=True,
        )
        assert result.decision_trail is not None
        return result.decision_trail

    @pytest.fixture
    def query(self, trail):
        return ProvenanceQuery(trail)

    def test_constraints_that_reduced_cai_returns_list(self, query):
        """constraints_that_reduced_cai should return a list."""
        reducers = query.constraints_that_reduced_cai()
        assert isinstance(reducers, list)

    def test_constraints_that_reduced_cai_items_are_constraint_decisions(self, query):
        """Each item should be a ConstraintDecision."""
        reducers = query.constraints_that_reduced_cai()
        for cd in reducers:
            assert isinstance(cd, ConstraintDecision)

    def test_constraints_that_reduced_cai_negative_impact(self, query):
        """Each returned constraint should have negative impact_on_cai."""
        reducers = query.constraints_that_reduced_cai()
        for cd in reducers:
            assert cd.impact_on_cai < 0, (
                f"Constraint {cd.constraint_name} should have negative CAI impact, "
                f"got {cd.impact_on_cai}"
            )

    def test_explain_position_returns_string(self, trail):
        """explain_position should return a human-readable string."""
        explanation = explain_position(trail, position=0)
        assert isinstance(explanation, str)

    def test_explain_position_mentions_codon_choice(self, trail):
        """Explanation for position 0 should mention the codon choice."""
        explanation = explain_position(trail, position=0)
        # Position 0 is Methionine (M) in our protein
        assert "0" in explanation, "Should mention position number"
        # Should mention the chosen codon or amino acid
        query = ProvenanceQuery(trail)
        first_decision = query.decisions_at_position(0)
        assert (
            first_decision.chosen_codon in explanation
            or first_decision.amino_acid in explanation
            or "Methionine" in explanation
        ), (
            f"Explanation should mention codon choice ({first_decision.chosen_codon}), "
            f"amino acid ({first_decision.amino_acid}), or 'Methionine'. "
            f"Got: {explanation[:200]}"
        )

    def test_explain_position_mentions_reasoning(self, trail):
        """Explanation should mention the reasoning behind the codon choice."""
        explanation = explain_position(trail, position=0)
        # Should mention the reason or context (CAI, constraint, etc.)
        query = ProvenanceQuery(trail)
        first_decision = query.decisions_at_position(0)
        assert (
            "CAI" in explanation
            or "chosen" in explanation.lower()
            or first_decision.constraint_reason.lower().replace("_", " ") in explanation.lower()
        ), (
            "Explanation should mention the reasoning for the codon choice. "
            f"Got: {explanation[:200]}"
        )

    def test_explain_position_out_of_range_raises(self, trail):
        """Explaining a position that doesn't exist should raise KeyError."""
        with pytest.raises(KeyError):
            explain_position(trail, position=9999)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Backward compatibility
# ═══════════════════════════════════════════════════════════════════════════════

class TestBackwardCompatibility:
    """Optimize with provenance tracking OFF and verify it still works."""

    @pytest.fixture
    def result(self):
        return optimize_sequence(
            target_protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            track_provenance=False,
        )

    def test_decision_trail_is_none(self, result):
        """decision_trail should be None when tracking is OFF."""
        assert result.decision_trail is None, (
            "decision_trail should be None when track_provenance=False"
        )

    def test_optimization_still_works(self, result):
        """Optimization should produce valid results even without provenance."""
        assert result.sequence, "Should produce a DNA sequence"
        assert len(result.sequence) == len(SHORT_PROTEIN) * 3, (
            "Sequence length should be protein_length * 3"
        )
        assert 0.0 <= result.cai <= 1.0, (
            f"CAI should be in [0, 1], got {result.cai}"
        )
        assert 0.0 <= result.gc_content <= 1.0, (
            f"GC content should be in [0, 1], got {result.gc_content}"
        )

    def test_protein_preserved(self, result):
        """The optimized sequence should still encode the correct protein."""
        from biocompiler.translation import translate
        protein = translate(result.sequence)
        assert protein == SHORT_PROTEIN, (
            f"Protein should be preserved, expected {SHORT_PROTEIN}, got {protein}"
        )

    def test_provenance_record_still_populated(self, result):
        """The provenance field (OptimizationRecord) should still be populated
        for backward compatibility, even when decision_trail is None."""
        assert result.provenance is not None, (
            "provenance (OptimizationRecord) should still be populated"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Report generation
# ═══════════════════════════════════════════════════════════════════════════════

class TestReportGeneration:
    """Generate markdown report from a decision trail and verify content."""

    @pytest.fixture
    def trail(self):
        result = optimize_sequence(
            target_protein=SHORT_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            track_provenance=True,
        )
        assert result.decision_trail is not None
        return result.decision_trail

    def test_generate_markdown_returns_string(self, trail):
        """ProvenanceReport.generate_markdown should return a string."""
        report = ProvenanceReport.generate_markdown(trail)
        assert isinstance(report, str)

    def test_report_has_summary_section(self, trail):
        """Report should contain an Executive Summary section."""
        report = ProvenanceReport.generate_markdown(trail)
        # The report uses "Executive Summary" or "Summary"
        assert "Summary" in report, "Report should have a Summary section"

    def test_report_has_constraint_impact_section(self, trail):
        """Report should contain a Constraint Impact section."""
        report = ProvenanceReport.generate_markdown(trail)
        assert "Constraint" in report, "Report should have constraint information"

    def test_report_has_codon_usage_section(self, trail):
        """Report should contain a Codon Usage section."""
        report = ProvenanceReport.generate_markdown(trail)
        assert "Codon" in report, "Report should have codon usage information"

    def test_report_contains_cai(self, trail):
        """Report should mention the CAI value."""
        report = ProvenanceReport.generate_markdown(trail)
        assert "CAI" in report, "Report should mention CAI"

    def test_report_contains_gc(self, trail):
        """Report should mention the GC content."""
        report = ProvenanceReport.generate_markdown(trail)
        assert "GC" in report, "Report should mention GC"

    def test_report_contains_organism(self, trail):
        """Report should mention the organism."""
        report = ProvenanceReport.generate_markdown(trail)
        assert "Escherichia_coli" in report, "Report should mention the organism"

    def test_generate_json_returns_string(self, trail):
        """ProvenanceReport.generate_json should return a string."""
        report = ProvenanceReport.generate_json(trail)
        assert isinstance(report, str)

    def test_generate_html_returns_string(self, trail):
        """ProvenanceReport.generate_html should return a string."""
        report = ProvenanceReport.generate_html(trail)
        assert isinstance(report, str)

    def test_report_does_not_crash_on_minimal_trail(self):
        """Report generation should not crash on any valid input."""
        # Create a minimal trail
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="M",
            organism="Escherichia_coli",
            constraints=[],
        )
        collector.record_codon_decision(CodonDecision(
            position=0,
            amino_acid="M",
            original_codon=None,
            chosen_codon="ATG",
            alternatives_considered=[],
            constraint_reason="maximize_cai",
            confidence=1.0,
        ))
        trail = collector.finalize(output_dna="ATG", cai=1.0, gc=1.0)

        # Should not crash
        md_report = ProvenanceReport.generate_markdown(trail)
        assert isinstance(md_report, str)

        json_report = ProvenanceReport.generate_json(trail)
        assert isinstance(json_report, str)

        html_report = ProvenanceReport.generate_html(trail)
        assert isinstance(html_report, str)

    def test_report_does_not_crash_on_trail_with_many_alternatives(self):
        """Report generation should not crash with complex alternatives."""
        collector = DecisionProvenanceCollector()
        collector.start_optimization(
            protein="LL",
            organism="Escherichia_coli",
            constraints=["GCInRange"],
        )
        collector.record_codon_decision(CodonDecision(
            position=0,
            amino_acid="L",
            original_codon=None,
            chosen_codon="CTG",
            alternatives_considered=[
                {"codon": "TTA", "cai_contribution": 0.1, "gc_contribution": 0.0,
                 "violates_constraints": [], "rejected_because": "Lower CAI"},
                {"codon": "TTG", "cai_contribution": 0.2, "gc_contribution": 0.33,
                 "violates_constraints": [], "rejected_because": "Lower CAI"},
                {"codon": "CTC", "cai_contribution": 0.8, "gc_contribution": 1.0,
                 "violates_constraints": [], "rejected_because": "Lower CAI"},
                {"codon": "CTA", "cai_contribution": 0.3, "gc_contribution": 0.33,
                 "violates_constraints": [], "rejected_because": "Lower CAI"},
                {"codon": "CTT", "cai_contribution": 0.15, "gc_contribution": 0.0,
                 "violates_constraints": [], "rejected_because": "Lower CAI"},
                {"codon": "CTG", "cai_contribution": 1.0, "gc_contribution": 0.67,
                 "violates_constraints": [], "rejected_because": None},
            ],
            constraint_reason="maximize_cai",
            confidence=1.0,
        ))
        collector.record_codon_decision(CodonDecision(
            position=1,
            amino_acid="L",
            original_codon=None,
            chosen_codon="CTG",
            alternatives_considered=[
                {"codon": "TTA", "cai_contribution": 0.1, "gc_contribution": 0.0,
                 "violates_constraints": [], "rejected_because": "Lower CAI"},
                {"codon": "CTC", "cai_contribution": 0.8, "gc_contribution": 1.0,
                 "violates_constraints": ["gc_content"], "rejected_because": "Violates: gc_content"},
                {"codon": "CTG", "cai_contribution": 1.0, "gc_contribution": 0.67,
                 "violates_constraints": [], "rejected_because": None},
            ],
            constraint_reason="gc_content",
            confidence=0.7,
        ))
        collector.record_constraint_decision(ConstraintDecision(
            constraint_name="GCInRange",
            constraint_type="hard",
            action_taken="satisfied",
            positions_affected=[0, 1],
            tradeoff_description="GC kept in [0.30, 0.70] by choosing CTG over CTC",
            impact_on_cai=-0.005,
        ))
        trail = collector.finalize(output_dna="CTGCTG", cai=0.88, gc=0.67)

        md_report = ProvenanceReport.generate_markdown(trail)
        assert isinstance(md_report, str)
        assert "GCInRange" in md_report

        json_report = ProvenanceReport.generate_json(trail)
        assert isinstance(json_report, str)
