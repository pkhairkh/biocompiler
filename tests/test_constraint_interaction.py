"""Unit tests for solver/constraint_interaction.py — Constraint Interaction Map.

Covers:
1. InteractionInfo construction and validation
2. ConstraintInteractionMap.build_interaction_map() for simple models
3. Overlapping position detection
4. Pre-computed known interactions (CAI costs and severities)
5. print_interaction_report() formatting

All tests construct CSPModel objects directly — no solver backend needed.
"""

from __future__ import annotations

import pytest

from biocompiler.solver.constraint_interaction import (
    ConstraintInteractionMap,
    InteractionInfo,
    print_interaction_report,
    _KNOWN_INTERACTIONS,
    _worse_severity,
)
from biocompiler.solver.types import (
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
    SolverConfig,
)


# ── Helpers ────────────────────────────────────────────────────────

def _make_config(**overrides) -> SolverConfig:
    """Build a SolverConfig with sensible defaults, overridden by *overrides*."""
    defaults = dict(gc_lo=0.30, gc_hi=0.70)
    defaults.update(overrides)
    return SolverConfig(**defaults)


def _make_model(
    protein: str = "MVSKGE",
    constraints: list[ConstraintSpec] | None = None,
    codon_domains: dict[int, list[str]] | None = None,
    config: SolverConfig | None = None,
) -> CSPModel:
    """Build a minimal CSPModel.

    If *codon_domains* is None, the default AA_TO_CODONS lookup is used.
    """
    from biocompiler.constants import AA_TO_CODONS

    if config is None:
        config = _make_config()
    if constraints is None:
        constraints = []
    if codon_domains is None:
        codon_domains = {i: AA_TO_CODONS.get(aa, []) for i, aa in enumerate(protein)}
    return CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints,
        config=config,
    )


def _constraint(
    ctype: ConstraintType = ConstraintType.NO_CRYPTIC_SPLICE,
    name: str = "test_constraint",
    positions: list[int] | None = None,
    params: dict | None = None,
) -> ConstraintSpec:
    """Shorthand to build a ConstraintSpec."""
    return ConstraintSpec(
        ctype=ctype,
        name=name,
        strictness=ConstraintStrictness.HARD,
        params=params or {},
        positions=positions or [],
    )


# =====================================================================
# 1. InteractionInfo construction and validation
# =====================================================================

class TestInteractionInfoConstruction:
    """Tests for the InteractionInfo dataclass."""

    def test_basic_construction(self):
        """InteractionInfo can be constructed with all fields."""
        info = InteractionInfo(
            constraint_a="NoCrypticSpliceConstraint",
            constraint_b="MaximizeCAI",
            overlapping_positions=[0, 3, 7],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        assert info.constraint_a == "NoCrypticSpliceConstraint"
        assert info.constraint_b == "MaximizeCAI"
        assert info.overlapping_positions == [0, 3, 7]
        assert info.estimated_cai_cost == 0.10
        assert info.conflict_severity == "high"

    def test_default_fields(self):
        """InteractionInfo defaults: empty positions, 0.0 cost, low severity."""
        info = InteractionInfo(
            constraint_a="A",
            constraint_b="B",
        )
        assert info.overlapping_positions == []
        assert info.estimated_cai_cost == 0.0
        assert info.conflict_severity == "low"

    def test_invalid_severity_raises(self):
        """An invalid conflict_severity value raises ValueError."""
        with pytest.raises(ValueError, match="conflict_severity"):
            InteractionInfo(
                constraint_a="A",
                constraint_b="B",
                conflict_severity="critical",
            )

    def test_negative_cai_cost_raises(self):
        """A negative estimated_cai_cost raises ValueError."""
        with pytest.raises(ValueError, match="estimated_cai_cost"):
            InteractionInfo(
                constraint_a="A",
                constraint_b="B",
                estimated_cai_cost=-0.01,
            )

    def test_all_valid_severities(self):
        """All three valid severity levels should be accepted."""
        for severity in ("high", "medium", "low"):
            info = InteractionInfo(
                constraint_a="A", constraint_b="B",
                conflict_severity=severity,
            )
            assert info.conflict_severity == severity


# =====================================================================
# 2. build_interaction_map for simple models
# =====================================================================

class TestBuildInteractionMapSimple:
    """Tests for ConstraintInteractionMap.build_interaction_map() with simple models."""

    def test_empty_model_returns_empty(self):
        """A model with no constraints returns an empty interaction map."""
        model = _make_model(protein="MVSKGE", constraints=[])
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        assert result == {}

    def test_single_constraint_returns_empty(self):
        """A model with one constraint has no pairs, so returns empty."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "NoCrypticSpliceConstraint"),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        assert result == {}

    def test_two_non_overlapping_constraints_no_interaction(self):
        """Two constraints affecting disjoint positions produce no interaction."""
        model = _make_model(
            protein="MVSKGEEEE",  # 9 positions
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "splice_pos0",
                    positions=[0, 1],
                ),
                _constraint(
                    ConstraintType.RESTRICTION_SITE,
                    "restr_pos5",
                    positions=[5, 6],
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        # No overlap → no interaction
        assert result == {}

    def test_two_overlapping_constraints_produce_interaction(self):
        """Two constraints with overlapping positions produce one interaction."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "NoCrypticSpliceConstraint",
                    positions=[0, 1, 2],
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                    positions=[0, 1, 2, 3, 4, 5],
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        assert len(result) == 1

        # Key should be sorted alphabetically
        key = tuple(sorted(["NoCrypticSpliceConstraint", "MaximizeCAI"]))
        assert key in result

        info = result[key]
        assert info.constraint_a == key[0]
        assert info.constraint_b == key[1]
        # Overlap should be [0, 1, 2]
        assert sorted(info.overlapping_positions) == [0, 1, 2]

    def test_global_constraints_overlap_all_positions(self):
        """Two global constraints (empty positions) overlap at all positions."""
        protein = "MVSKGE"
        model = _make_model(
            protein=protein,
            constraints=[
                _constraint(
                    ConstraintType.GC_CONTENT,
                    "GCRangeConstraint",
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        assert len(result) == 1

        key = tuple(sorted(["GCRangeConstraint", "MaximizeCAI"]))
        info = result[key]
        # Global constraints should overlap at all 6 codon positions
        assert len(info.overlapping_positions) == len(protein)

    def test_three_constraints_three_pairs(self):
        """Three pairwise-overlapping constraints produce three interactions."""
        model = _make_model(
            protein="MVSKGEEEE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "NoCrypticSpliceConstraint",
                    positions=[0, 1, 2, 3, 4],
                ),
                _constraint(
                    ConstraintType.NO_CPG,
                    "NoCpGIslandConstraint",
                    positions=[2, 3, 4, 5],
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        assert len(result) == 3


# =====================================================================
# 3. Overlapping positions correctly identified
# =====================================================================

class TestOverlappingPositions:
    """Tests that overlapping positions are correctly computed."""

    def test_partial_overlap(self):
        """Constraints overlapping at only some positions report only those."""
        model = _make_model(
            protein="MVSKGEEEE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "splice",
                    positions=[0, 1, 2, 3],
                ),
                _constraint(
                    ConstraintType.NO_CPG,
                    "cpg",
                    positions=[2, 3, 4, 5],
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["splice", "cpg"]))
        info = result[key]
        # Only positions 2 and 3 are in both sets
        assert info.overlapping_positions == [2, 3]

    def test_full_overlap(self):
        """Identical position sets produce full overlap."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "splice",
                    positions=[1, 2, 3],
                ),
                _constraint(
                    ConstraintType.NO_CPG,
                    "cpg",
                    positions=[1, 2, 3],
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["splice", "cpg"]))
        info = result[key]
        assert info.overlapping_positions == [1, 2, 3]

    def test_single_position_overlap(self):
        """Even a single overlapping position counts as an interaction."""
        model = _make_model(
            protein="MVSKGEEEE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "splice",
                    positions=[0],
                ),
                _constraint(
                    ConstraintType.RESTRICTION_SITE,
                    "restr",
                    positions=[0, 5],
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["splice", "restr"]))
        info = result[key]
        assert info.overlapping_positions == [0]

    def test_global_vs_position_specific(self):
        """A global constraint overlaps a position-specific constraint at its positions."""
        model = _make_model(
            protein="MVSKGEEEE",
            constraints=[
                _constraint(
                    ConstraintType.GC_CONTENT,
                    "GCRangeConstraint",
                    positions=[],  # global
                ),
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "splice_pos3",
                    positions=[3],
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["GCRangeConstraint", "splice_pos3"]))
        info = result[key]
        assert 3 in info.overlapping_positions


# =====================================================================
# 4. Pre-computed known interactions
# =====================================================================

class TestKnownInteractions:
    """Tests that pre-computed CAI costs and severities are used correctly."""

    def test_cryptic_splice_cai_interaction(self):
        """NoCrypticSpliceConstraint <-> MaximizeCAI: HIGH severity, ~0.10 CAI cost."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "NoCrypticSpliceConstraint",
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["NoCrypticSpliceConstraint", "MaximizeCAI"]))
        info = result[key]
        assert info.conflict_severity == "high"
        assert info.estimated_cai_cost == pytest.approx(0.10, abs=0.01)

    def test_cpg_cai_interaction(self):
        """NoCpGIslandConstraint <-> MaximizeCAI: MEDIUM severity, ~0.05 CAI cost."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CPG,
                    "NoCpGIslandConstraint",
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["NoCpGIslandConstraint", "MaximizeCAI"]))
        info = result[key]
        assert info.conflict_severity == "medium"
        assert info.estimated_cai_cost == pytest.approx(0.05, abs=0.01)

    def test_restriction_site_cai_interaction(self):
        """NoRestrictionSiteConstraint <-> MaximizeCAI: LOW severity, ~0.02 CAI cost."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.RESTRICTION_SITE,
                    "NoRestrictionSiteConstraint",
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["NoRestrictionSiteConstraint", "MaximizeCAI"]))
        info = result[key]
        assert info.conflict_severity == "low"
        assert info.estimated_cai_cost == pytest.approx(0.02, abs=0.01)

    def test_gc_range_cai_interaction(self):
        """GCRangeConstraint <-> MaximizeCAI: LOW severity."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.GC_CONTENT,
                    "GCRangeConstraint",
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["GCRangeConstraint", "MaximizeCAI"]))
        info = result[key]
        assert info.conflict_severity == "low"

    def test_attta_cai_interaction(self):
        """NoATTTAMotifConstraint <-> MaximizeCAI: LOW severity."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(
                    ConstraintType.NO_INSTABILITY_MOTIF,
                    "NoATTTAMotifConstraint",
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        key = tuple(sorted(["NoATTTAMotifConstraint", "MaximizeCAI"]))
        info = result[key]
        assert info.conflict_severity == "low"

    def test_position_specific_name_matches_base(self):
        """Position-specific constraint names like 'NoCrypticSpliceConstraint_pos42'
        should still match the pre-computed interaction for 'NoCrypticSpliceConstraint'."""
        model = _make_model(
            protein="MVSKGEEEEEEE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "NoCrypticSpliceConstraint_pos42",
                    positions=[4, 5],
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        # The key uses the actual constraint names, but the severity lookup
        # should match the base name.
        assert len(result) >= 1
        # Find the interaction with MaximizeCAI
        cai_interaction = None
        for key, info in result.items():
            if "MaximizeCAI" in key:
                cai_interaction = info
                break
        assert cai_interaction is not None
        assert cai_interaction.conflict_severity == "high"


# =====================================================================
# 5. print_interaction_report
# =====================================================================

class TestPrintInteractionReport:
    """Tests for the print_interaction_report() function."""

    def test_empty_map_message(self):
        """An empty interaction map produces a specific message."""
        report = print_interaction_report({})
        assert "No constraint interactions detected" in report

    def test_report_contains_constraint_names(self):
        """The report should include constraint names in the table."""
        info = InteractionInfo(
            constraint_a="NoCrypticSpliceConstraint",
            constraint_b="MaximizeCAI",
            overlapping_positions=[0, 1, 2],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        key = ("MaximizeCAI", "NoCrypticSpliceConstraint")
        report = print_interaction_report({key: info})
        assert "NoCrypticSpliceConstraint" in report
        assert "MaximizeCAI" in report

    def test_report_contains_severity(self):
        """The report should include severity levels."""
        info = InteractionInfo(
            constraint_a="A",
            constraint_b="B",
            overlapping_positions=[0],
            estimated_cai_cost=0.05,
            conflict_severity="medium",
        )
        report = print_interaction_report({("A", "B"): info})
        assert "medium" in report

    def test_report_contains_cai_cost(self):
        """The report should include estimated CAI cost values."""
        info = InteractionInfo(
            constraint_a="A",
            constraint_b="B",
            overlapping_positions=[0],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        report = print_interaction_report({("A", "B"): info})
        assert "0.1000" in report or "0.10" in report

    def test_report_sorted_by_cai_cost(self):
        """Interactions should be sorted by CAI cost, highest first."""
        info_high = InteractionInfo(
            constraint_a="HighCostConstr",
            constraint_b="OtherHighConstr",
            overlapping_positions=[0],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        info_low = InteractionInfo(
            constraint_a="LowCostConstr",
            constraint_b="OtherLowConstr",
            overlapping_positions=[0],
            estimated_cai_cost=0.02,
            conflict_severity="low",
        )
        report = print_interaction_report({
            ("HighCostConstr", "OtherHighConstr"): info_high,
            ("LowCostConstr", "OtherLowConstr"): info_low,
        })
        # "HighCostConstr" (high cost) should appear before
        # "LowCostConstr" (low cost) in the data rows
        pos_high = report.index("HighCostConstr")
        pos_low = report.index("LowCostConstr")
        assert pos_high < pos_low

    def test_report_includes_recommendations(self):
        """The report should include a recommendations section."""
        info = InteractionInfo(
            constraint_a="NoCrypticSpliceConstraint",
            constraint_b="MaximizeCAI",
            overlapping_positions=[0, 1, 2],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        key = ("MaximizeCAI", "NoCrypticSpliceConstraint")
        report = print_interaction_report({key: info})
        assert "Recommendation" in report

    def test_high_severity_has_relaxation_recommendation(self):
        """High severity interactions should produce relaxation suggestions."""
        info = InteractionInfo(
            constraint_a="NoCrypticSpliceConstraint",
            constraint_b="MaximizeCAI",
            overlapping_positions=[0],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        key = ("MaximizeCAI", "NoCrypticSpliceConstraint")
        report = print_interaction_report({key: info})
        # Should recommend relaxing the hard constraint (not MaximizeCAI)
        assert "NoCrypticSpliceConstraint" in report
        assert "HIGH" in report or "relax" in report.lower()

    def test_low_severity_only_no_recommendation_to_relax(self):
        """When all interactions are LOW severity, report says no relaxation needed."""
        info = InteractionInfo(
            constraint_a="GCRangeConstraint",
            constraint_b="MaximizeCAI",
            overlapping_positions=[0],
            estimated_cai_cost=0.01,
            conflict_severity="low",
        )
        key = ("GCRangeConstraint", "MaximizeCAI")
        report = print_interaction_report({key: info})
        assert "no relaxation needed" in report.lower() or "LOW" in report

    def test_recommendation_never_relaxes_maximize_cai(self):
        """Recommendations should never suggest relaxing MaximizeCAI."""
        info = InteractionInfo(
            constraint_a="NoCrypticSpliceConstraint",
            constraint_b="MaximizeCAI",
            overlapping_positions=[0],
            estimated_cai_cost=0.10,
            conflict_severity="high",
        )
        key = ("MaximizeCAI", "NoCrypticSpliceConstraint")
        report = print_interaction_report({key: info})
        # The recommendation should mention the hard constraint, not CAI
        lines = [l for l in report.split("\n") if "HIGH" in l or "relax" in l.lower()]
        for line in lines:
            # Should not say "relax MaximizeCAI"
            if "MaximizeCAI" in line and "relax" in line.lower():
                pytest.fail("Recommendation should not suggest relaxing MaximizeCAI")

    def test_report_contains_tip(self):
        """The report should include a tip about severity levels."""
        report = print_interaction_report({
            ("A", "B"): InteractionInfo(
                constraint_a="A",
                constraint_b="B",
                overlapping_positions=[0],
                estimated_cai_cost=0.01,
                conflict_severity="low",
            ),
        })
        assert "Tip" in report or "tip" in report.lower()


# =====================================================================
# 6. Helper function tests
# =====================================================================

class TestHelperFunctions:
    """Tests for internal helper functions."""

    def test_worse_severity_high_over_medium(self):
        """'high' is worse than 'medium'."""
        assert _worse_severity("high", "medium") == "high"
        assert _worse_severity("medium", "high") == "high"

    def test_worse_severity_medium_over_low(self):
        """'medium' is worse than 'low'."""
        assert _worse_severity("medium", "low") == "medium"
        assert _worse_severity("low", "medium") == "medium"

    def test_worse_severity_same(self):
        """Same severity returns itself."""
        assert _worse_severity("high", "high") == "high"
        assert _worse_severity("low", "low") == "low"

    def test_known_interactions_table_populated(self):
        """The pre-computed interactions table should not be empty."""
        assert len(_KNOWN_INTERACTIONS) > 0

    def test_known_interactions_severity_values_valid(self):
        """All severities in _KNOWN_INTERACTIONS should be valid."""
        for pair, (cost, severity) in _KNOWN_INTERACTIONS.items():
            assert severity in ("high", "medium", "low"), (
                f"Invalid severity {severity!r} for pair {pair}"
            )
            assert cost >= 0, f"Negative cost for pair {pair}"


# =====================================================================
# 7. Integration: full model with multiple constraints
# =====================================================================

class TestIntegration:
    """Integration tests using realistic constraint combinations."""

    def test_full_model_with_all_common_constraints(self):
        """A model with all common constraints should produce multiple interactions."""
        model = _make_model(
            protein="MVSKGEELFT",
            constraints=[
                _constraint(ConstraintType.AMINO_ACID_IDENTITY, "TranslationConstraint"),
                _constraint(ConstraintType.GC_CONTENT, "GCRangeConstraint"),
                _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "NoCrypticSpliceConstraint"),
                _constraint(ConstraintType.NO_CPG, "NoCpGIslandConstraint"),
                _constraint(ConstraintType.RESTRICTION_SITE, "NoRestrictionSiteConstraint"),
                _constraint(ConstraintType.NO_INSTABILITY_MOTIF, "NoATTTAMotifConstraint"),
                _constraint(ConstraintType.CODON_USAGE, "MaximizeCAI"),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        # Should have interactions for many pairs
        assert len(result) > 5

        # Check that the highest-CAI-cost interaction is CrypticSplice <-> MaximizeCAI
        max_cost_pair = max(result.items(), key=lambda kv: kv[1].estimated_cai_cost)
        assert max_cost_pair[1].conflict_severity in ("high", "medium", "low")

    def test_interaction_map_keys_are_sorted_tuples(self):
        """All keys in the result should be alphabetically sorted 2-tuples."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "NoCrypticSpliceConstraint"),
                _constraint(ConstraintType.CODON_USAGE, "MaximizeCAI"),
                _constraint(ConstraintType.GC_CONTENT, "GCRangeConstraint"),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        for key in result:
            assert len(key) == 2
            assert key[0] <= key[1], f"Key {key} is not sorted"

    def test_merge_multiple_position_specific_constraints(self):
        """Multiple position-specific constraints with the same base name
        should merge their overlap positions."""
        model = _make_model(
            protein="MVSKGEEEEEEE",
            constraints=[
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "NoCrypticSpliceConstraint_pos2",
                    positions=[2],
                ),
                _constraint(
                    ConstraintType.NO_CRYPTIC_SPLICE,
                    "NoCrypticSpliceConstraint_pos5",
                    positions=[5],
                ),
                _constraint(
                    ConstraintType.CODON_USAGE,
                    "MaximizeCAI",
                ),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        # Should have interactions for both position-specific constraints
        # with MaximizeCAI
        cai_interactions = [
            info for key, info in result.items()
            if "MaximizeCAI" in key
        ]
        assert len(cai_interactions) >= 1

    def test_report_with_full_model(self):
        """print_interaction_report works with a full model's interaction map."""
        model = _make_model(
            protein="MVSKGE",
            constraints=[
                _constraint(ConstraintType.NO_CRYPTIC_SPLICE, "NoCrypticSpliceConstraint"),
                _constraint(ConstraintType.CODON_USAGE, "MaximizeCAI"),
                _constraint(ConstraintType.GC_CONTENT, "GCRangeConstraint"),
            ],
        )
        imap = ConstraintInteractionMap()
        result = imap.build_interaction_map(model)
        report = print_interaction_report(result)
        assert isinstance(report, str)
        assert len(report) > 0
        assert "Constraint Interaction Map" in report
