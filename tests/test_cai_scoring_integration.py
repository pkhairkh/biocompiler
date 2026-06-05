"""Tests for CAI-aware scoring integration (Task F2.2).

Validates that:
1. SoftConstraintScorer.score_cai_impact() correctly computes CAI deltas
   for alternative codons at various positions.
2. ConstraintEnforcer.enforce_with_cai_awareness() picks the least-CAI-damaging
   fix when resolving constraint violations.
3. ConflictProvenance records include CAI impact data (both heuristic
   cai_impact and measured cai_delta).
"""

from __future__ import annotations

import math
import pytest

from biocompiler.solver.constraints import (
    CSPModel,
    CodonVariable,
    MaximizeCAI,
    MinimizeCpG,
    MinimizeMRNADG,
    NoRestrictionSiteConstraint,
    NoTRunConstraint,
    TranslationConstraint,
    CAI_LOG_EPSILON,
)
from biocompiler.solver.types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintViolation,
    SolverConfig,
)
from biocompiler.solver.scoring import SoftConstraintScorer, CAIImpactResult
from biocompiler.solver.enforcement import (
    ConstraintEnforcer,
    CAIAwareFix,
    CAIAwareEnforcementResult,
)
from biocompiler.solver.conflict_provenance import (
    ConflictProvenance,
    ConflictResolverWithProvenance,
)
from biocompiler.constants import AA_TO_CODONS, CODON_TABLE


# ==============================================================================
# Test fixtures
# ==============================================================================

# Simple adaptiveness table: some codons have high w, others low w
# For amino acid 'A' (Ala): GCT=0.9, GCC=0.8, GCA=0.3, GCG=0.1
# For amino acid 'V' (Val): GTG=0.9, GTC=0.7, GTT=0.4, GTA=0.1
# For amino acid 'L' (Leu): CTT=0.8, CTG=0.7, CTC=0.6, CTA=0.2,
#                           TTA=0.1, TTG=0.5
# For amino acid 'M' (Met): ATG=1.0 (only one codon)
SIMPLE_ADAPTIVENESS: dict[str, float] = {
    # Ala
    "GCT": 0.9, "GCC": 0.8, "GCA": 0.3, "GCG": 0.1,
    # Val
    "GTG": 0.9, "GTC": 0.7, "GTT": 0.4, "GTA": 0.1,
    # Leu
    "CTT": 0.8, "CTG": 0.7, "CTC": 0.6, "CTA": 0.2,
    "TTA": 0.1, "TTG": 0.5,
    # Met
    "ATG": 1.0,
    # Fill in all other codons with moderate values
    "TTT": 0.5, "TTC": 0.6,
    "ATT": 0.5, "ATC": 0.6, "ATA": 0.3,
    "TCT": 0.5, "TCC": 0.6, "TCA": 0.3, "TCG": 0.2,
    "CCT": 0.5, "CCC": 0.6, "CCA": 0.3, "CCG": 0.2,
    "ACT": 0.5, "ACC": 0.6, "ACA": 0.3, "ACG": 0.2,
    "TAT": 0.5, "TAC": 0.6,
    "CAT": 0.5, "CAC": 0.6, "CAA": 0.3, "CAG": 0.7,
    "AAT": 0.5, "AAC": 0.6, "AAA": 0.3, "AAG": 0.7,
    "GAT": 0.5, "GAC": 0.6, "GAA": 0.3, "GAG": 0.7,
    "TGT": 0.5, "TGC": 0.6, "TGG": 1.0,
    "CGT": 0.5, "CGC": 0.6, "CGA": 0.2, "CGG": 0.3,
    "AGT": 0.5, "AGC": 0.6, "AGA": 0.2, "AGG": 0.3,
    "GGT": 0.5, "GGC": 0.6, "GGA": 0.3, "GGG": 0.4,
}


def _make_model(
    protein: str = "MVA",
    adaptiveness: dict[str, float] | None = None,
    config: SolverConfig | None = None,
    hard_constraints: list | None = None,
) -> CSPModel:
    """Build a minimal CSPModel for testing."""
    if adaptiveness is None:
        adaptiveness = SIMPLE_ADAPTIVENESS
    if config is None:
        config = SolverConfig()
    if hard_constraints is None:
        hard_constraints = [TranslationConstraint(protein)]

    # Build codon variables
    variables = []
    for i, aa in enumerate(protein):
        domain = AA_TO_CODONS.get(aa, ["ATG"])
        variables.append(CodonVariable(position=i, amino_acid=aa, domain=domain))

    # Build soft constraints
    soft_constraints = [
        MaximizeCAI(adaptiveness, protein),
        MinimizeCpG(),
        MinimizeMRNADG(),
    ]

    return CSPModel(
        variables=variables,
        hard_constraints=hard_constraints,
        soft_constraints=soft_constraints,
        protein=protein,
        organism="Homo_sapiens",
        config=config,
    )


def _model_sequence(model: CSPModel, codons: list[str]) -> str:
    """Build a DNA sequence from a list of codons matching the model's protein."""
    return "".join(codons)


# ==============================================================================
# 1. Tests for SoftConstraintScorer.score_cai_impact()
# ==============================================================================

class TestScoreCAIImpact:
    """Test score_cai_impact() for various positions and codons."""

    def test_high_adaptiveness_codon_gives_positive_delta(self):
        """Replacing a low-adaptiveness codon with a high one should give positive delta."""
        model = _make_model("MVA")
        # M=ATG(w=1.0), V=GTG(w=0.9), A=GCA(w=0.3)
        sequence = "ATGGTGGCA"
        scorer = SoftConstraintScorer(model.config)

        # At position 2 (A), GCA(w=0.3) -> GCT(w=0.9)
        deltas = scorer.score_cai_impact(model, sequence, 2, ["GCT", "GCC", "GCG"])

        n = len(model.protein)  # 3
        # delta = (log(0.9) - log(0.3)) / 3 for GCT
        expected_gct = (math.log(0.9) - math.log(0.3)) / n
        assert deltas["GCT"] == pytest.approx(expected_gct, rel=1e-6)
        assert deltas["GCT"] > 0  # improvement

    def test_low_adaptiveness_codon_gives_negative_delta(self):
        """Replacing a high-adaptiveness codon with a low one should give negative delta."""
        model = _make_model("MVA")
        # A position 2: GCT(w=0.9) -> GCG(w=0.1)
        sequence = "ATGGTGGCT"
        scorer = SoftConstraintScorer(model.config)

        deltas = scorer.score_cai_impact(model, sequence, 2, ["GCA", "GCG"])

        n = len(model.protein)
        expected_gcg = (math.log(0.1) - math.log(0.9)) / n
        assert deltas["GCG"] == pytest.approx(expected_gcg, rel=1e-6)
        assert deltas["GCG"] < 0  # loss

    def test_same_adaptiveness_codon_gives_zero_delta(self):
        """Replacing with same adaptiveness should give approximately zero delta."""
        model = _make_model("MVA")
        # Position 2: GCC(w=0.8) -> GCA(w=0.3) and GCC(w=0.8)
        sequence = "ATGGTGGCC"
        scorer = SoftConstraintScorer(model.config)

        # Compare GCC with itself (same codon)
        deltas = scorer.score_cai_impact(model, sequence, 2, ["GCC"])

        assert deltas["GCC"] == pytest.approx(0.0, abs=1e-10)

    def test_single_amino_acid_model(self):
        """Test with a single-codon protein (Met only)."""
        model = _make_model("M")
        sequence = "ATG"
        scorer = SoftConstraintScorer(model.config)

        # ATG is the only codon for Met — no alternatives
        deltas = scorer.score_cai_impact(model, sequence, 0, [])
        assert deltas == {}

    def test_out_of_range_position_raises(self):
        """Position out of range should raise ValueError."""
        model = _make_model("MVA")
        sequence = "ATGGTGGCA"
        scorer = SoftConstraintScorer(model.config)

        with pytest.raises(ValueError, match="out of range"):
            scorer.score_cai_impact(model, sequence, 5, ["GCT"])

        with pytest.raises(ValueError, match="out of range"):
            scorer.score_cai_impact(model, sequence, -1, ["GCT"])

    def test_no_maximize_cai_constraint_returns_zero(self):
        """When model has no MaximizeCAI, all deltas should be 0."""
        config = SolverConfig()
        protein = "MVA"
        variables = []
        for i, aa in enumerate(protein):
            domain = AA_TO_CODONS.get(aa, ["ATG"])
            variables.append(CodonVariable(position=i, amino_acid=aa, domain=domain))

        model = CSPModel(
            variables=variables,
            hard_constraints=[TranslationConstraint(protein)],
            soft_constraints=[MinimizeCpG(), MinimizeMRNADG()],  # No MaximizeCAI
            protein=protein,
            organism="Homo_sapiens",
            config=config,
        )

        sequence = "ATGGTGGCA"
        scorer = SoftConstraintScorer(model.config)
        deltas = scorer.score_cai_impact(model, sequence, 2, ["GCT", "GCC"])

        assert all(v == 0.0 for v in deltas.values())

    def test_longer_protein_smaller_deltas(self):
        """CAI deltas should be smaller for longer proteins (divided by N)."""
        model_short = _make_model("MVA")  # N=3
        model_long = _make_model("MVAMVAMVA")  # N=9

        # Both have the same A at position 2 (short) and position 2 (long)
        sequence_short = "ATGGTGGCA"
        # For the long one, build a matching sequence
        sequence_long = "ATGGTGGCA" + "ATGGTGGCA" + "ATGGTGGCA"

        scorer_short = SoftConstraintScorer(model_short.config)
        scorer_long = SoftConstraintScorer(model_long.config)

        deltas_short = scorer_short.score_cai_impact(
            model_short, sequence_short, 2, ["GCT"]
        )
        deltas_long = scorer_long.score_cai_impact(
            model_long, sequence_long, 2, ["GCT"]
        )

        # The delta for the longer protein should be 1/3 of the short one
        ratio = deltas_long["GCT"] / deltas_short["GCT"]
        assert ratio == pytest.approx(3.0 / 9.0, rel=1e-6)

    def test_multiple_alternatives_sorted_by_delta(self):
        """Verify that alternatives can be sorted by CAI delta."""
        model = _make_model("MVA")
        sequence = "ATGGTGGCA"  # A at pos 2 = GCA(w=0.3)
        scorer = SoftConstraintScorer(model.config)

        deltas = scorer.score_cai_impact(
            model, sequence, 2, ["GCT", "GCC", "GCA", "GCG"]
        )

        # Sort by delta descending (best first)
        sorted_alts = sorted(deltas.items(), key=lambda kv: kv[1], reverse=True)

        # GCT(w=0.9) should be best, then GCC(w=0.8), then GCA(same=0), then GCG(w=0.1)
        assert sorted_alts[0][0] == "GCT"
        assert sorted_alts[1][0] == "GCC"
        assert sorted_alts[-1][0] == "GCG"


# ==============================================================================
# 2. Tests for ConstraintEnforcer.enforce_with_cai_awareness()
# ==============================================================================

class TestEnforceWithCAIAwareness:
    """Test that enforce_with_cai_awareness() picks the least-CAI-damaging fix."""

    def test_no_violations_returns_unchanged(self):
        """When all constraints are satisfied, no fixes are applied."""
        model = _make_model("MVA")
        sequence = "ATGGTGGCA"  # Valid sequence
        enforcer = ConstraintEnforcer()

        result = enforcer.enforce_with_cai_awareness(model, sequence)

        assert result.all_hard_satisfied is True
        assert result.fixes == []
        assert result.total_cai_delta == 0.0
        assert result.remaining_violations == []
        assert result.sequence == sequence

    def test_fix_prefers_least_cai_damaging(self):
        """When fixing a restriction site, the fix should prefer the least-CAI-damaging codon."""
        # Create a model with a restriction site constraint for "GCT"
        # Protein: "A" -> codons: GCT, GCC, GCA, GCG
        # GCT contains "GC" which is part of many restriction sites
        # We'll use a restriction site that matches a substring in the sequence

        # Build a sequence where a restriction site exists and can be fixed
        # Protein = "AA" -> we need a restriction site that spans the sequence
        # Let's use a site like "GCTGC" which would be formed by GCT at pos 0 + GC from pos 1
        # Simpler: just avoid "GCT" itself as a "restriction site" (4-letter site)

        # Actually let's use a real test: force EcoRI (GAATTC) into the sequence
        # That's hard to construct. Let's instead test with NoTRunConstraint
        # which is easier to violate and fix.

        # Protein: "LLL" (Leu has 6 codons)
        # Use a sequence that creates a T-run
        # Leu codons: CTT, CTC, CTA, CTG, TTA, TTG
        # "CTTCTTTTA" = CTT CTT TTA — no T-run > 5
        # "CTTCTTTTG" = CTT CTT TTG — no T-run > 5
        # Let's create: "TTTTTATTTTT" — but that's not a valid protein
        # Better: protein "LL" with sequence "CTTTTA" has CTT + TTA
        # The cross-codon "TTTTA" has TTTT (4 Ts) — not > 5

        # Let me use a simpler approach: NoRestrictionSiteConstraint
        # with a site that appears in a specific codon choice.

        # Protein = "AA", sequence = "GCTGCT" (no site violated by default)
        # Add restriction site "GCTG" — that spans from codon 0 into codon 1
        # But restriction site checking is on the full sequence string

        # Simplest: use a site that's literally a codon. Say "GCT" as a "site".
        # Sequence: "GCTGCC" for protein "AA"
        # If we avoid "GCT", we can use GCC, GCA, or GCG instead.

        protein = "AA"
        hard_constraints = [
            TranslationConstraint(protein),
            NoRestrictionSiteConstraint(["GCTG"]),  # Site spans pos 0-3
        ]

        model = _make_model(
            protein=protein,
            hard_constraints=hard_constraints,
        )

        # Sequence "GCTGCC" contains "GCTG" — violates restriction site
        # GCC(w=0.8) at pos 0 would give "GCCGCC" which doesn't contain "GCTG"
        # GCA(w=0.3) at pos 0 would give "GCAGCC" which doesn't contain "GCTG"
        # GCG(w=0.1) at pos 0 would give "GCGGCC" — no "GCTG" either
        # The best CAI fix is GCC (w=0.8, least CAI damage from GCT w=0.9)

        sequence = "GCTGCC"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        # Should have applied a fix
        assert len(result.fixes) > 0

        # The fix should have chosen the least-CAI-damaging alternative
        # GCT(w=0.9) -> GCC(w=0.8) is the best alternative that resolves it
        fix = result.fixes[0]
        assert fix.old_codon == "GCT"
        # GCC should be preferred over GCA and GCG since it has highest
        # adaptiveness among the alternatives
        assert fix.new_codon == "GCC"
        assert fix.cai_delta < 0  # Some CAI was lost (0.9 -> 0.8)

        # The restriction site should now be resolved
        assert result.all_hard_satisfied is True

    def test_cai_aware_fix_records_cai_delta(self):
        """Each fix should record its CAI delta."""
        protein = "AA"
        hard_constraints = [
            TranslationConstraint(protein),
            NoRestrictionSiteConstraint(["GCTG"]),
        ]
        model = _make_model(protein=protein, hard_constraints=hard_constraints)

        sequence = "GCTGCC"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        for fix in result.fixes:
            assert isinstance(fix.cai_delta, float)
            assert isinstance(fix.codon_position, int)
            assert isinstance(fix.old_codon, str)
            assert isinstance(fix.new_codon, str)
            assert isinstance(fix.constraint_name, str)

    def test_total_cai_delta_is_sum_of_fixes(self):
        """total_cai_delta should be the sum of all fix deltas."""
        protein = "AA"
        hard_constraints = [
            TranslationConstraint(protein),
            NoRestrictionSiteConstraint(["GCTG"]),
        ]
        model = _make_model(protein=protein, hard_constraints=hard_constraints)

        sequence = "GCTGCC"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        expected_total = sum(f.cai_delta for f in result.fixes)
        assert result.total_cai_delta == pytest.approx(expected_total, rel=1e-10)

    def test_result_type_is_cai_aware_enforcement_result(self):
        """Result should be a CAIAwareEnforcementResult."""
        model = _make_model("MVA")
        sequence = "ATGGTGGCA"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        assert isinstance(result, CAIAwareEnforcementResult)

    def test_fix_type_is_cai_aware_fix(self):
        """Each fix should be a CAIAwareFix."""
        protein = "AA"
        hard_constraints = [
            TranslationConstraint(protein),
            NoRestrictionSiteConstraint(["GCTG"]),
        ]
        model = _make_model(protein=protein, hard_constraints=hard_constraints)

        sequence = "GCTGCC"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        for fix in result.fixes:
            assert isinstance(fix, CAIAwareFix)

    def test_no_valid_alternative_leaves_remaining_violation(self):
        """When no alternative codon can fix the violation, it remains."""
        # Protein "M" — Met only has one codon (ATG). If we add a
        # restriction site that includes "ATG", there's no alternative.
        protein = "M"
        hard_constraints = [
            TranslationConstraint(protein),
            NoRestrictionSiteConstraint(["ATG"]),  # Impossible to avoid!
        ]
        model = _make_model(protein=protein, hard_constraints=hard_constraints)

        sequence = "ATG"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        # Should have remaining violations since ATG is the only Met codon
        assert len(result.remaining_violations) > 0
        assert result.all_hard_satisfied is False


# ==============================================================================
# 3. Tests for ConflictProvenance CAI impact data
# ==============================================================================

class TestProvenanceCAIImpact:
    """Test that provenance records include CAI impact data."""

    def test_conflict_provenance_has_cai_delta_field(self):
        """ConflictProvenance should have a cai_delta field."""
        record = ConflictProvenance(
            conflicting_constraints=["GCRangeConstraint", "MaximizeCAI"],
            resolution_method="priority_based",
            winner="GCRangeConstraint",
            loser="MaximizeCAI",
            impact="Test impact",
            positions_affected=[0, 1, 2],
            cai_impact=-0.05,
            cai_delta=-0.003,
        )

        assert record.cai_delta == -0.003
        assert record.cai_impact == -0.05

    def test_conflict_provenance_cai_delta_defaults_to_none(self):
        """cai_delta should default to None when not provided."""
        record = ConflictProvenance(
            conflicting_constraints=["A", "B"],
            resolution_method="priority_based",
            winner="A",
            loser="B",
            impact="Test",
            positions_affected=[],
        )

        assert record.cai_delta is None

    def test_cai_aware_resolution_method_is_valid(self):
        """The 'cai_aware' resolution method should be accepted."""
        record = ConflictProvenance(
            conflicting_constraints=["NoRestrictionSiteConstraint", "MaximizeCAI"],
            resolution_method="cai_aware",
            winner="NoRestrictionSiteConstraint",
            loser="MaximizeCAI",
            impact="Fixed with minimal CAI loss",
            positions_affected=[3, 4, 5],
            cai_impact=0.01,
            cai_delta=-0.002,
        )

        assert record.resolution_method == "cai_aware"

    def test_invalid_resolution_method_rejected(self):
        """Invalid resolution methods should be rejected."""
        with pytest.raises(ValueError, match="Invalid resolution_method"):
            ConflictProvenance(
                conflicting_constraints=["A", "B"],
                resolution_method="invalid_method",
                winner="A",
                loser="B",
                impact="Test",
                positions_affected=[],
            )

    def test_record_cai_aware_provenance(self):
        """ConflictResolverWithProvenance.record_cai_aware_provenance should work."""
        resolver = ConflictResolverWithProvenance(track_provenance=True)

        record = resolver.record_cai_aware_provenance(
            constraint_name="NoRestrictionSiteConstraint",
            codon_position=2,
            old_codon="GCT",
            new_codon="GCC",
            cai_delta=-0.0015,
        )

        assert record.resolution_method == "cai_aware"
        assert record.winner == "NoRestrictionSiteConstraint"
        assert record.loser == "MaximizeCAI"
        assert record.cai_delta == -0.0015
        assert "GCT" in record.impact
        assert "GCC" in record.impact
        assert "NoRestrictionSiteConstraint" in record.impact
        # positions_affected defaults to codon_position * 3 + [0,1,2]
        assert record.positions_affected == [6, 7, 8]

    def test_record_cai_aware_provenance_custom_positions(self):
        """Custom positions_affected should be used when provided."""
        resolver = ConflictResolverWithProvenance(track_provenance=True)

        record = resolver.record_cai_aware_provenance(
            constraint_name="GCRangeConstraint",
            codon_position=5,
            old_codon="GCG",
            new_codon="GCA",
            cai_delta=-0.01,
            positions_affected=[15, 16],
        )

        assert record.positions_affected == [15, 16]

    def test_record_cai_aware_provenance_positive_delta(self):
        """When cai_delta is positive, impact should note CAI improvement."""
        resolver = ConflictResolverWithProvenance(track_provenance=True)

        record = resolver.record_cai_aware_provenance(
            constraint_name="NoTRunConstraint",
            codon_position=3,
            old_codon="CTT",
            new_codon="CTG",
            cai_delta=0.005,
        )

        assert "improved" in record.impact.lower() or "maintained" in record.impact.lower()

    def test_record_cai_aware_provenance_negative_delta(self):
        """When cai_delta is negative, impact should note CAI sacrifice."""
        resolver = ConflictResolverWithProvenance(track_provenance=True)

        record = resolver.record_cai_aware_provenance(
            constraint_name="NoRestrictionSiteConstraint",
            codon_position=0,
            old_codon="GCT",
            new_codon="GCA",
            cai_delta=-0.02,
        )

        assert "sacrificed" in record.impact.lower() or "least-damaging" in record.impact.lower()

    def test_record_relaxation_provenance_has_cai_delta_none(self):
        """Legacy record_relaxation_provenance should have cai_delta=None."""
        resolver = ConflictResolverWithProvenance(track_provenance=True)

        record = resolver.record_relaxation_provenance(
            relaxed_constraint_name="GCRangeConstraint",
            kept_constraint_name="NoRestrictionSiteConstraint",
            positions_affected=[10, 11, 12],
            sequence="ATGGTGGCA",
        )

        assert record.cai_delta is None
        assert isinstance(record.cai_impact, float)

    def test_provenance_repr_includes_cai_delta_when_set(self):
        """repr should include cai_delta when it's set."""
        record = ConflictProvenance(
            conflicting_constraints=["A", "B"],
            resolution_method="cai_aware",
            winner="A",
            loser="B",
            impact="Test",
            positions_affected=[0],
            cai_impact=-0.05,
            cai_delta=-0.003,
        )

        r = repr(record)
        assert "cai_delta=" in r
        assert "-0.003000" in r

    def test_provenance_repr_omits_cai_delta_when_none(self):
        """repr should omit cai_delta when it's None."""
        record = ConflictProvenance(
            conflicting_constraints=["A", "B"],
            resolution_method="priority_based",
            winner="A",
            loser="B",
            impact="Test",
            positions_affected=[0],
            cai_impact=-0.05,
        )

        r = repr(record)
        assert "cai_delta=" not in r


# ==============================================================================
# 4. Integration test: end-to-end CAI-aware enforcement with provenance
# ==============================================================================

class TestCAIAwareScoringIntegration:
    """End-to-end integration tests for the CAI-aware scoring pipeline."""

    def test_full_pipeline_enforcement_to_provenance(self):
        """Test the full pipeline: enforce -> record provenance -> verify records."""
        protein = "AA"
        hard_constraints = [
            TranslationConstraint(protein),
            NoRestrictionSiteConstraint(["GCTG"]),
        ]
        model = _make_model(protein=protein, hard_constraints=hard_constraints)

        sequence = "GCTGCC"
        enforcer = ConstraintEnforcer()
        result = enforcer.enforce_with_cai_awareness(model, sequence)

        # Should have fixed the violation
        assert result.all_hard_satisfied is True
        assert len(result.fixes) > 0

        # Record provenance for each fix
        resolver = ConflictResolverWithProvenance(track_provenance=True)
        provenance_records = []
        for fix in result.fixes:
            record = resolver.record_cai_aware_provenance(
                constraint_name=fix.constraint_name,
                codon_position=fix.codon_position,
                old_codon=fix.old_codon,
                new_codon=fix.new_codon,
                cai_delta=fix.cai_delta,
            )
            provenance_records.append(record)

        # Verify provenance records
        assert len(provenance_records) == len(result.fixes)
        for record in provenance_records:
            assert record.resolution_method == "cai_aware"
            assert record.cai_delta is not None
            assert isinstance(record.cai_delta, float)
            assert record.loser == "MaximizeCAI"

    def test_cai_impact_result_dataclass(self):
        """Test the CAIImpactResult dataclass."""
        result = CAIImpactResult(
            position=2,
            current_codon="GCA",
            cai_deltas={"GCT": 0.05, "GCC": 0.03, "GCG": -0.08},
        )

        assert result.position == 2
        assert result.current_codon == "GCA"
        assert result.best_codon == "GCT"  # highest delta
        assert result.best_delta == pytest.approx(0.05, rel=1e-6)

    def test_cai_impact_result_empty_deltas(self):
        """CAIImpactResult with empty deltas should have no best_codon."""
        result = CAIImpactResult(
            position=0,
            current_codon="ATG",
            cai_deltas={},
        )

        assert result.best_codon is None
        assert result.best_delta == 0.0

    def test_cai_aware_result_repr(self):
        """Test repr of CAIAwareEnforcementResult."""
        result = CAIAwareEnforcementResult(
            sequence="ATGGCC",
            fixes=[
                CAIAwareFix(
                    codon_position=0,
                    old_codon="GCT",
                    new_codon="GCC",
                    cai_delta=-0.01,
                    constraint_name="NoRestrictionSiteConstraint",
                )
            ],
            total_cai_delta=-0.01,
            remaining_violations=[],
            all_hard_satisfied=True,
        )

        r = repr(result)
        assert "fixes=1" in r
        assert "hard_satisfied=True" in r

    def test_cai_aware_fix_repr(self):
        """Test repr of CAIAwareFix."""
        fix = CAIAwareFix(
            codon_position=2,
            old_codon="GCT",
            new_codon="GCC",
            cai_delta=-0.005,
            constraint_name="NoRestrictionSiteConstraint",
        )

        r = repr(fix)
        assert "GCT" in r
        assert "GCC" in r
        assert "NoRestrictionSiteConstraint" in r
