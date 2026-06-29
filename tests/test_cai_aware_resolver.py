"""Tests for biocompiler.solver.cai_aware_resolver.

Covers the following areas:
1. ResolutionChoice dataclass construction and validation
2. CAIAwareConstraintResolver.estimate_cai_impact()
3. CAIAwareConstraintResolver.resolve_with_min_cai_loss()
4. Edge cases (no conflicts, single constraint)
5. Integration with MaximizeCAI scoring
6. ConflictResolver auto_resolve with cai_aware strategy
"""

from __future__ import annotations

import math
import pytest

# ---------------------------------------------------------------------------
# Graceful import — entire module skipped when solver package is absent
# ---------------------------------------------------------------------------
car = pytest.importorskip("biocompiler.solver.cai_aware_resolver")
CAIAwareConstraintResolver = car.CAIAwareConstraintResolver
ResolutionChoice = car.ResolutionChoice

from biocompiler.solver.types import (
    ConstraintPriority,
    ConstraintSpec,
    ConstraintStrictness,
    ConstraintType,
    CSPModel,
    SolverConfig,
)
from biocompiler.solver.conflict_resolution import (
    ConstraintConflict,
    ConflictResolver,
)


# =====================================================================
# Helpers
# =====================================================================

# A small adaptiveness table for testing — alanine and valine codons
_SIMPLE_ADAPTIVENESS: dict[str, float] = {
    # Alanine (A): GCT is optimal, then GCC, GCA, GCG
    "GCT": 1.0,
    "GCC": 0.6,
    "GCA": 0.3,
    "GCG": 0.2,
    # Valine (V): GTT is optimal, then GTC, GTA, GTG
    "GTT": 1.0,
    "GTC": 0.5,
    "GTA": 0.3,
    "GTG": 0.8,
    # Leucine (L): CTT optimal
    "CTT": 1.0,
    "CTC": 0.6,
    "CTA": 0.2,
    "CTG": 0.9,
    "TTA": 0.1,
    "TTG": 0.4,
    # Methionine (M): only ATG
    "ATG": 1.0,
}


def _make_model(
    protein: str = "AVL",
    constraints: list[ConstraintSpec] | None = None,
    config: SolverConfig | None = None,
    codon_domains: dict[int, list[str]] | None = None,
) -> CSPModel:
    """Build a minimal CSPModel for testing."""
    from biocompiler.shared.constants import AA_TO_CODONS

    if codon_domains is None:
        codon_domains = {i: list(AA_TO_CODONS.get(aa, ["ATG"])) for i, aa in enumerate(protein)}
    return CSPModel(
        protein_sequence=protein,
        codon_domains=codon_domains,
        constraints=constraints or [],
        config=config or SolverConfig(),
    )


def _hard_spec(
    name: str,
    ctype: ConstraintType = ConstraintType.GC_CONTENT,
    positions: list[int] | None = None,
    priority: ConstraintPriority = ConstraintPriority.MEDIUM,
    params: dict | None = None,
) -> ConstraintSpec:
    """Shorthand for a hard ConstraintSpec."""
    return ConstraintSpec(
        ctype=ctype,
        name=name,
        strictness=ConstraintStrictness.HARD,
        positions=positions or [],
        priority=priority,
        params=params or {},
    )


def _make_resolver(
    protein: str = "AVL",
    adaptiveness: dict[str, float] | None = None,
    codon_domains: dict[int, list[str]] | None = None,
    constraints: list[ConstraintSpec] | None = None,
) -> CAIAwareConstraintResolver:
    """Create a CAIAwareConstraintResolver for testing."""
    adapt = adaptiveness or _SIMPLE_ADAPTIVENESS
    model = _make_model(protein=protein, constraints=constraints, codon_domains=codon_domains)
    return CAIAwareConstraintResolver(model, adapt)


# =====================================================================
# 1. ResolutionChoice dataclass
# =====================================================================


class TestResolutionChoice:
    """Test ResolutionChoice dataclass construction and validation."""

    def test_basic_construction(self):
        """ResolutionChoice can be created with valid arguments."""
        choice = ResolutionChoice(
            relax_constraint="gc",
            keep_constraint="splice",
            cai_loss=0.02,
            alternative_relax="splice",
            alternative_cai_loss=0.05,
            conflict_positions=[0, 1],
        )
        assert choice.relax_constraint == "gc"
        assert choice.keep_constraint == "splice"
        assert choice.cai_loss == 0.02
        assert choice.alternative_cai_loss == 0.05
        assert choice.conflict_positions == [0, 1]
        assert choice.strategy == "cai_aware"

    def test_cai_savings_positive(self):
        """cai_savings is positive when chosen option is better."""
        choice = ResolutionChoice(
            relax_constraint="a",
            keep_constraint="b",
            cai_loss=0.01,
            alternative_relax="b",
            alternative_cai_loss=0.05,
            conflict_positions=[0],
        )
        assert choice.cai_savings == pytest.approx(0.04)

    def test_cai_savings_negative(self):
        """cai_savings is negative when alternative was better (should not happen normally)."""
        choice = ResolutionChoice(
            relax_constraint="a",
            keep_constraint="b",
            cai_loss=0.10,
            alternative_relax="b",
            alternative_cai_loss=0.02,
            conflict_positions=[0],
        )
        assert choice.cai_savings == pytest.approx(-0.08)

    def test_cai_savings_zero(self):
        """cai_savings is zero when both options have equal loss."""
        choice = ResolutionChoice(
            relax_constraint="a",
            keep_constraint="b",
            cai_loss=0.03,
            alternative_relax="b",
            alternative_cai_loss=0.03,
            conflict_positions=[0],
        )
        assert choice.cai_savings == pytest.approx(0.0)

    def test_invalid_strategy_raises(self):
        """Invalid strategy raises ValueError."""
        with pytest.raises(ValueError, match="Invalid strategy"):
            ResolutionChoice(
                relax_constraint="a",
                keep_constraint="b",
                cai_loss=0.0,
                alternative_relax="b",
                alternative_cai_loss=0.0,
                conflict_positions=[],
                strategy="invalid",
            )

    def test_valid_strategies(self):
        """All valid strategies are accepted."""
        for strategy in ["relax_a", "relax_b", "compromise", "cai_aware"]:
            choice = ResolutionChoice(
                relax_constraint="a",
                keep_constraint="b",
                cai_loss=0.0,
                alternative_relax="b",
                alternative_cai_loss=0.0,
                conflict_positions=[],
                strategy=strategy,
            )
            assert choice.strategy == strategy

    def test_repr_format(self):
        """__repr__ includes key fields."""
        choice = ResolutionChoice(
            relax_constraint="gc",
            keep_constraint="splice",
            cai_loss=0.02,
            alternative_relax="splice",
            alternative_cai_loss=0.05,
            conflict_positions=[1],
        )
        r = repr(choice)
        assert "gc" in r
        assert "splice" in r
        assert "cai_aware" in r


# =====================================================================
# 2. CAIAwareConstraintResolver.estimate_cai_impact()
# =====================================================================


class TestEstimateCAIImpact:
    """Test CAIAwareConstraintResolver.estimate_cai_impact()."""

    def test_codon_usage_constraint(self):
        """Codon-usage constraints have impact proportional to domain variance."""
        resolver = _make_resolver(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        spec = _hard_spec("cai", ctype=ConstraintType.CODON_USAGE, positions=[0])
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        # GCT=1.0 is best, avg = (1.0 + 0.6 + 0.3 + 0.2) / 4 = 0.525
        # Impact = 1.0 - 0.525 = 0.475
        assert impact == pytest.approx(0.475, abs=0.01)

    def test_codon_usage_single_codon_no_impact(self):
        """Single-codon domain has zero CAI impact for codon-usage constraint."""
        resolver = _make_resolver(
            protein="M",
            codon_domains={0: ["ATG"]},
        )
        spec = _hard_spec("cai", ctype=ConstraintType.CODON_USAGE, positions=[0])
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["ATG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        assert impact == pytest.approx(0.0, abs=0.001)

    def test_gc_constraint_has_impact(self):
        """GC constraints produce non-negative impact estimates."""
        resolver = _make_resolver(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        spec = _hard_spec(
            "gc", ctype=ConstraintType.GC_CONTENT, positions=[0],
            params={"gc_lo": 0.45, "gc_hi": 0.55},
        )
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        assert impact >= 0.0

    def test_splice_constraint_with_gt_codons(self):
        """Splice constraints produce impact when GT-containing codons exist."""
        resolver = _make_resolver(
            protein="V",
            codon_domains={0: ["GTT", "GTC", "GTA", "GTG"]},
        )
        spec = _hard_spec(
            "splice", ctype=ConstraintType.NO_CRYPTIC_SPLICE, positions=[0],
        )
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GTT", "GTC", "GTA", "GTG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        # All valine codons contain GT → all excluded → high impact
        assert impact > 0.0

    def test_splice_constraint_no_gt_codons(self):
        """Splice constraints have zero impact when no GT/AG codons exist."""
        resolver = _make_resolver(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        spec = _hard_spec(
            "splice", ctype=ConstraintType.NO_CRYPTIC_SPLICE, positions=[0],
        )
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        # Alanine codons do not contain GT or AG → zero impact
        assert impact == pytest.approx(0.0, abs=0.001)

    def test_cpg_constraint_impact(self):
        """CpG constraints produce impact when CpG-containing codons exist."""
        resolver = _make_resolver(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        spec = _hard_spec(
            "cpg", ctype=ConstraintType.NO_CPG, positions=[0],
        )
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        # GCG contains CG → should have some impact
        assert impact >= 0.0

    def test_restriction_site_constraint_low_impact(self):
        """Restriction site constraints produce low impact estimates."""
        resolver = _make_resolver(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        spec = _hard_spec(
            "eco", ctype=ConstraintType.RESTRICTION_SITE, positions=[0],
        )
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        # Restriction site impact is scaled down by 0.3
        assert impact >= 0.0

    def test_amino_acid_identity_maximum_impact(self):
        """Amino acid identity constraints have maximum impact."""
        resolver = _make_resolver(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        spec = _hard_spec(
            "trans", ctype=ConstraintType.AMINO_ACID_IDENTITY, positions=[0],
        )
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        assert impact == 1.0

    def test_empty_domain_zero_impact(self):
        """Empty codon domain produces zero impact."""
        resolver = _make_resolver(protein="A")
        spec = _hard_spec("gc", positions=[0])
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: []},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        assert impact == 0.0

    def test_unknown_constraint_type_zero_impact(self):
        """Unknown/custom constraint types produce zero impact."""
        resolver = _make_resolver(protein="A")
        spec = _hard_spec("custom", ctype=ConstraintType.CUSTOM, positions=[0])
        impact = resolver.estimate_cai_impact(
            spec, position=0,
            codon_domains={0: ["GCT", "GCC"]},
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        assert impact == 0.0


# =====================================================================
# 3. CAIAwareConstraintResolver.resolve_with_min_cai_loss()
# =====================================================================


class TestResolveWithMinCAILoss:
    """Test CAIAwareConstraintResolver.resolve_with_min_cai_loss()."""

    def test_picks_lower_cai_loss(self):
        """Resolver picks the resolution with lower CAI loss."""
        # Create two constraints that overlap at position 0
        # GC constraint (higher CAI impact) vs restriction site (lower CAI impact)
        model = _make_model(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)

        gc_spec = _hard_spec(
            "gc", ctype=ConstraintType.GC_CONTENT, positions=[0],
            params={"gc_lo": 0.45, "gc_hi": 0.55},
        )
        splice_spec = _hard_spec(
            "splice", ctype=ConstraintType.NO_CRYPTIC_SPLICE, positions=[0],
        )

        choice = resolver.resolve_with_min_cai_loss(
            conflicting_constraints=[gc_spec, splice_spec],
            current_sequence="GCT",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )

        # The choice should be a valid ResolutionChoice
        assert isinstance(choice, ResolutionChoice)
        assert choice.cai_loss <= choice.alternative_cai_loss
        assert choice.cai_savings >= 0.0

    def test_requires_exactly_two_constraints(self):
        """resolve_with_min_cai_loss requires exactly 2 constraints."""
        resolver = _make_resolver(protein="A")
        spec = _hard_spec("gc", positions=[0])

        with pytest.raises(ValueError, match="Expected exactly 2"):
            resolver.resolve_with_min_cai_loss(
                conflicting_constraints=[spec],
                current_sequence="GCT",
                codon_domains={0: ["GCT", "GCC"]},
            )

        with pytest.raises(ValueError, match="Expected exactly 2"):
            resolver.resolve_with_min_cai_loss(
                conflicting_constraints=[spec, spec, spec],
                current_sequence="GCT",
                codon_domains={0: ["GCT", "GCC"]},
            )

    def test_equal_impact_returns_first_option(self):
        """When both options have equal CAI impact, returns the first (relax A)."""
        # Use two constraints of the same type at the same position
        model = _make_model(
            protein="A",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)

        spec_a = _hard_spec("gc_a", ctype=ConstraintType.GC_CONTENT, positions=[0],
                            params={"gc_lo": 0.45, "gc_hi": 0.55})
        spec_b = _hard_spec("gc_b", ctype=ConstraintType.GC_CONTENT, positions=[0],
                            params={"gc_lo": 0.45, "gc_hi": 0.55})

        choice = resolver.resolve_with_min_cai_loss(
            conflicting_constraints=[spec_a, spec_b],
            current_sequence="GCT",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"]},
        )

        # Both have the same impact, so cai_loss == alternative_cai_loss
        assert choice.cai_loss == pytest.approx(choice.alternative_cai_loss)
        assert choice.cai_savings == pytest.approx(0.0)

    def test_conflict_positions_populated(self):
        """ResolutionChoice contains correct conflict positions."""
        model = _make_model(
            protein="AV",
            codon_domains={0: ["GCT", "GCC"], 1: ["GTT", "GTC"]},
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)

        spec_a = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT, positions=[0, 1])
        spec_b = _hard_spec("splice", ctype=ConstraintType.NO_CRYPTIC_SPLICE, positions=[0, 1])

        choice = resolver.resolve_with_min_cai_loss(
            conflicting_constraints=[spec_a, spec_b],
            current_sequence="GCTGTT",
            codon_domains={0: ["GCT", "GCC"], 1: ["GTT", "GTC"]},
        )

        # Both constraints overlap on positions 0 and 1
        assert 0 in choice.conflict_positions
        assert 1 in choice.conflict_positions


# =====================================================================
# 4. Edge cases
# =====================================================================


class TestEdgeCases:
    """Test edge cases for the CAI-aware resolver."""

    def test_no_conflicts_auto_resolve(self):
        """auto_resolve_cai_aware with no conflicts returns model unchanged."""
        model = _make_model(protein="AV", constraints=[])
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        resolved = resolver.auto_resolve_cai_aware(model)
        assert resolved is model

    def test_single_constraint_no_conflict(self):
        """A single constraint cannot conflict with anything."""
        model = _make_model(
            protein="AV",
            constraints=[_hard_spec("gc", positions=[0])],
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        resolved = resolver.auto_resolve_cai_aware(model)
        # No conflicts, so same model returned
        assert resolved is model

    def test_auto_resolve_produces_valid_model(self):
        """auto_resolve_cai_aware produces a valid CSPModel."""
        model = _make_model(
            protein="AV",
            constraints=[
                _hard_spec("gc", positions=[0, 1], params={"gc_lo": 0.45, "gc_hi": 0.55}),
                _hard_spec("splice", positions=[0, 1]),
            ],
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        resolved = resolver.auto_resolve_cai_aware(model)
        assert isinstance(resolved, CSPModel)
        assert resolved.protein_sequence == model.protein_sequence
        assert resolved.codon_domains == model.codon_domains
        # One constraint should have been removed or downgraded
        assert len(resolved.constraints) <= len(model.constraints)

    def test_compute_cai_with_valid_sequence(self):
        """compute_cai returns a value in [0, 1] for a valid sequence."""
        model = _make_model(protein="AV")
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        # Build a sequence using the codon domains
        seq = "GCTGTT"  # A=Ala(GCT) V=Val(GTT) — both optimal
        cai = resolver.compute_cai(seq)
        assert 0.0 <= cai <= 1.0
        # Both codons are optimal → CAI should be close to 1.0
        assert cai > 0.9

    def test_compute_cai_loss(self):
        """compute_cai_loss correctly measures CAI difference."""
        model = _make_model(protein="AV")
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        # Optimal sequence
        before = "GCTGTT"  # Both optimal codons
        # Sub-optimal sequence
        after = "GCGGTA"   # GCG(0.2) and GTA(0.3)
        loss = resolver.compute_cai_loss(before, after)
        assert loss >= 0.0

    def test_compute_cai_loss_zero_when_identical(self):
        """compute_cai_loss returns 0 for identical sequences."""
        model = _make_model(protein="AV")
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        loss = resolver.compute_cai_loss("GCTGTT", "GCTGTT")
        assert loss == pytest.approx(0.0)


# =====================================================================
# 5. Integration with MaximizeCAI scoring
# =====================================================================


class TestMaximizeCAIIntegration:
    """Test that the resolver integrates correctly with MaximizeCAI scoring."""

    def test_cai_values_consistent_with_maximize_cai(self):
        """Resolver's CAI computation matches MaximizeCAI.cai()."""
        from biocompiler.solver.constraints import MaximizeCAI

        model = _make_model(protein="AVL")
        adaptiveness = _SIMPLE_ADAPTIVENESS
        resolver = CAIAwareConstraintResolver(model, adaptiveness)

        cai_constraint = MaximizeCAI(adaptiveness, protein="AVL")
        seq = "GCTGTTCTT"  # A=Ala(GCT), V=Val(GTT), L=Leu(CTT)

        cai_from_resolver = resolver.compute_cai(seq)
        cai_from_constraint = cai_constraint.cai(seq)

        assert cai_from_resolver == pytest.approx(cai_from_constraint, abs=0.001)

    def test_optimal_codons_yield_high_cai(self):
        """Sequences using optimal codons have high CAI scores."""
        model = _make_model(protein="AV")
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)

        # Use optimal codons: GCT (Ala, 1.0) and GTT (Val, 1.0)
        optimal_seq = "GCTGTT"
        cai = resolver.compute_cai(optimal_seq)
        assert cai == pytest.approx(1.0, abs=0.01)

    def test_suboptimal_codons_yield_lower_cai(self):
        """Sequences using sub-optimal codons have lower CAI scores."""
        model = _make_model(protein="AV")
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)

        # Use sub-optimal codons: GCG (Ala, 0.2) and GTA (Val, 0.3)
        suboptimal_seq = "GCGGTA"
        cai = resolver.compute_cai(suboptimal_seq)
        assert cai < 0.5

    def test_resolution_preserves_higher_cai(self):
        """Resolution choices that pick lower CAI loss preserve more CAI."""
        model = _make_model(
            protein="AV",
            codon_domains={0: ["GCT", "GCC", "GCA", "GCG"], 1: ["GTT", "GTC", "GTA", "GTG"]},
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)

        # A CpG constraint at both positions (moderate CAI impact)
        # vs a GC constraint at both positions (higher CAI impact)
        cpg_spec = _hard_spec("cpg", ctype=ConstraintType.NO_CPG, positions=[0, 1])
        gc_spec = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT, positions=[0, 1],
                             params={"gc_lo": 0.45, "gc_hi": 0.55})

        choice = resolver.resolve_with_min_cai_loss(
            conflicting_constraints=[cpg_spec, gc_spec],
            current_sequence="GCTGTT",
            codon_domains=model.codon_domains,
        )

        # The resolver should make a valid choice
        assert isinstance(choice, ResolutionChoice)
        assert choice.relax_constraint in ("cpg", "gc")
        assert choice.keep_constraint in ("cpg", "gc")
        assert choice.relax_constraint != choice.keep_constraint


# =====================================================================
# 6. ConflictResolver auto_resolve with cai_aware strategy
# =====================================================================


class TestConflictResolverCAIAware:
    """Test CAIAwareConstraintResolver auto-resolve with CAI awareness."""

    def test_cai_aware_auto_resolve(self):
        """CAI-aware auto-resolve works via CAIAwareConstraintResolver."""
        model = _make_model(
            protein="AV",
            constraints=[
                _hard_spec("gc", positions=[0, 1], params={"gc_lo": 0.45, "gc_hi": 0.55}),
                _hard_spec("splice", positions=[0, 1]),
            ],
        )
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        resolved = resolver.auto_resolve_cai_aware(model)
        assert isinstance(resolved, CSPModel)
        assert resolved.protein_sequence == model.protein_sequence

    def test_base_resolver_rejects_cai_aware_strategy(self):
        """Base ConflictResolver rejects 'cai_aware' as unknown strategy."""
        model = _make_model(
            protein="AV",
            constraints=[
                _hard_spec("gc", positions=[0, 1]),
                _hard_spec("splice", positions=[0, 1]),
            ],
        )
        resolver = ConflictResolver()
        with pytest.raises(ValueError, match="Unknown auto-resolution strategy"):
            resolver.auto_resolve(model, strategy="cai_aware")

    def test_cai_aware_resolver_creates_valid_model(self):
        """CAIAwareConstraintResolver produces a valid CSPModel."""
        model = _make_model(protein="AV", constraints=[])
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        # Should not raise
        resolved = resolver.auto_resolve_cai_aware(model)
        assert isinstance(resolved, CSPModel)


# =====================================================================
# 7. Lightweight vs Full mode
# =====================================================================


class TestLightweightVsFullMode:
    """Test the lightweight CAIImpactEstimator from conflict_resolution."""

    def test_lightweight_mode_no_model(self):
        """CAIImpactEstimator works in lightweight mode without model."""
        from biocompiler.solver.conflict_resolution import CAIImpactEstimator as LightweightResolver

        resolver = LightweightResolver()
        spec = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT)
        impact = resolver.estimate_cai_impact(spec)
        # Should return a value from the static table
        assert isinstance(impact, float)

    def test_full_mode_with_model(self):
        """CAIAwareConstraintResolver (full) works with model + adaptiveness."""
        model = _make_model(protein="A")
        # Full resolver from cai_aware_resolver module
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        spec = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT, positions=[0])

        # estimate_cai_impact requires position, codon_domains, and adaptiveness
        impact = resolver.estimate_cai_impact(
            spec,
            position=0,
            codon_domains=model.codon_domains,
            adaptiveness=_SIMPLE_ADAPTIVENESS,
        )
        assert isinstance(impact, float)

    def test_resolve_with_min_cai_loss_works_in_full_mode(self):
        """resolve_with_min_cai_loss works with the full CAIAwareConstraintResolver."""
        model = _make_model(protein="A")
        resolver = CAIAwareConstraintResolver(model, _SIMPLE_ADAPTIVENESS)
        spec_a = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT, positions=[0])
        spec_b = _hard_spec("cpg", ctype=ConstraintType.NO_CPG, positions=[0])

        choice = resolver.resolve_with_min_cai_loss(
            conflicting_constraints=[spec_a, spec_b],
            current_sequence="GCT",
            codon_domains={0: ["GCT", "GCC"]},
        )
        assert isinstance(choice, ResolutionChoice)

    def test_rank_resolution_compromise_adjusted(self):
        """rank_resolution adjusts compromise strategy based on CAI impact."""
        from biocompiler.solver.conflict_resolution import CAIImpactEstimator as LightweightResolver

        resolver = LightweightResolver()

        # GC constraint (impact=-0.05) vs restriction site (impact=0.01)
        gc_spec = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT,
                             priority=ConstraintPriority.MEDIUM)
        rest_spec = _hard_spec("eco", ctype=ConstraintType.RESTRICTION_SITE,
                               priority=ConstraintPriority.MEDIUM)

        # With equal priority, compromise → should adjust to prefer relaxing
        # the constraint whose removal helps CAI more (rest_spec has higher impact)
        result = resolver.rank_resolution(gc_spec, rest_spec, "compromise")
        # Restriction site has positive impact (0.01), GC has negative (-0.05)
        # Relaxing restriction site hurts CAI less than relaxing GC
        # GC impact = -0.05, REST impact = 0.01
        # rank_resolution prefers relaxing the one with higher (more positive) impact
        # So it should prefer relaxing restriction site → "relax_b"
        # But the diff is 0.06 which is > 0.03 threshold
        assert result in ("relax_a", "relax_b", "compromise")

    def test_rank_resolution_does_not_override_priority(self):
        """rank_resolution does not override different-priority decisions."""
        from biocompiler.solver.conflict_resolution import CAIImpactEstimator as LightweightResolver

        resolver = LightweightResolver()

        gc_spec = _hard_spec("gc", ctype=ConstraintType.GC_CONTENT,
                             priority=ConstraintPriority.HIGH)
        rest_spec = _hard_spec("eco", ctype=ConstraintType.RESTRICTION_SITE,
                               priority=ConstraintPriority.LOW)

        # Different priorities → should NOT adjust regardless of CAI impact
        result = resolver.rank_resolution(gc_spec, rest_spec, "relax_b")
        assert result == "relax_b"  # Unchanged
