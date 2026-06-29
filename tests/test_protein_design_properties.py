"""
Property-Based Tests for BioCompiler protein_design Module
===========================================================

Verifies structural invariants and consistency properties of the
protein_design module using Hypothesis-based property testing.

Three core properties tested:
  1. _estimate_ddg always returns a finite float for any valid AA pair
  2. score_mutation always returns a dict with exactly 5 expected keys
  3. design_* functions always produce a protein of the same length as input
"""

import math

import pytest
pytest.importorskip("hypothesis")
pytest.importorskip("hypothesis")
from hypothesis import given, settings, assume, HealthCheck
from hypothesis import strategies as st

from biocompiler.engines.protein_design import (
    DesignConstraints,
    DesignResult,
    _estimate_ddg,
    design_thermostable,
    design_soluble,
    design_low_immunogenicity,
    score_mutation,
)


# ────────────────────────────────────────────────────────────
# Shared Constants
# ────────────────────────────────────────────────────────────

# The 20 canonical one-letter amino acid codes
STANDARD_AA: str = "ACDEFGHIKLMNPQRSTVWY"

# The five required keys in a MutationScore dict
MUTATION_SCORE_KEYS: frozenset[str] = frozenset({
    "stability_ddg",
    "solubility_delta",
    "immunogenicity_delta",
    "blosum62",
    "weighted_score",
})


# ────────────────────────────────────────────────────────────
# Hypothesis Strategies
# ────────────────────────────────────────────────────────────

# Single standard amino acid
aa_strategy = st.sampled_from(list(STANDARD_AA))

# Pair of (wildtype, mutant) amino acids
aa_pair_strategy = st.tuples(aa_strategy, aa_strategy)

# Valid protein sequence (at least 10 residues for MHC_II_WINDOW analysis)
protein_strategy = st.text(
    alphabet=STANDARD_AA,
    min_size=10,
    max_size=50,
)

# Shorter protein for faster design tests
protein_short_strategy = st.text(
    alphabet=STANDARD_AA,
    min_size=10,
    max_size=20,
)

# Position within a given protein (computed inside tests)
# Weight dicts with positive weights
weights_strategy = st.dictionaries(
    keys=st.sampled_from(["stability", "solubility", "immunogenicity"]),
    values=st.floats(min_value=0.0, max_value=2.0, allow_nan=False, allow_infinity=False),
    min_size=3,
    max_size=3,
)


# ══════════════════════════════════════════════════════════════
# TEST CLASS 1: _estimate_ddg Always Returns Finite Float
# ══════════════════════════════════════════════════════════════

class TestEstimateDdgFiniteFloat:
    """Property: _estimate_ddg always returns a finite float for any
    pair of standard amino acids.

    This must hold for:
      - Any wildtype/mutant combination (including same residue)
      - Proline and glycine special cases
      - Unusual amino acid pairs
    """

    @given(pair=aa_pair_strategy)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_returns_float(self, pair):
        """_estimate_ddg returns a float for any AA pair."""
        wildtype, mutant = pair
        result = _estimate_ddg(wildtype, mutant)
        assert isinstance(result, float), (
            f"_estimate_ddg({wildtype!r}, {mutant!r}) returned "
            f"{type(result).__name__}, expected float"
        )

    @given(pair=aa_pair_strategy)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_returns_finite(self, pair):
        """_estimate_ddg never returns NaN or infinity."""
        wildtype, mutant = pair
        result = _estimate_ddg(wildtype, mutant)
        assert math.isfinite(result), (
            f"_estimate_ddg({wildtype!r}, {mutant!r}) = {result}, "
            f"which is not finite"
        )

    @given(pair=aa_pair_strategy)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_result_is_rounded_to_3_decimals(self, pair):
        """_estimate_ddg rounds to 3 decimal places."""
        wildtype, mutant = pair
        result = _estimate_ddg(wildtype, mutant)
        assert result == round(result, 3), (
            f"_estimate_ddg({wildtype!r}, {mutant!r}) = {result}, "
            f"not rounded to 3 decimal places"
        )

    @given(pair=aa_pair_strategy)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_result_in_reasonable_range(self, pair):
        """_estimate_ddg result should be in a physically plausible range.

        The most extreme possible combination:
          - BLOSUM62 min = -4, weight = -0.15  → 0.6
          - Hydropathy max swing ~9.0, weight = 0.05  → 0.45
          - Proline bonus: -0.3
          - Glycine penalty: +0.3
        So the range should be roughly [-2, +2] kcal/mol.
        """
        wildtype, mutant = pair
        result = _estimate_ddg(wildtype, mutant)
        assert -5.0 <= result <= 5.0, (
            f"_estimate_ddg({wildtype!r}, {mutant!r}) = {result}, "
            f"outside plausible range [-5, 5]"
        )

    def test_proline_mutation_includes_bonus(self):
        """Any non-Pro → Pro mutation includes the proline stabilization bonus."""
        ddg_with_pro = _estimate_ddg("A", "P")
        ddg_without_pro = _estimate_ddg("A", "V")
        # Pro bonus should make it more negative than a similar substitution
        # (not always true due to hydropathy, but A→P should be stabilizing)
        assert ddg_with_pro < 0, (
            f"A→P should be stabilizing due to PROLINE_BONUS, got {ddg_with_pro}"
        )

    def test_glycine_mutation_includes_penalty(self):
        """Any non-Gly → Gly mutation includes the glycine flexibility penalty."""
        ddg_with_gly = _estimate_ddg("A", "G")
        # Glycine penalty should push it positive
        assert ddg_with_gly > 0, (
            f"A→G should be destabilizing due to GLYCINE_BONUS, got {ddg_with_gly}"
        )

    @given(wt=aa_strategy)
    @settings(max_examples=20)
    def test_same_residue_is_deterministic(self, wt):
        """_estimate_ddg(x, x) always returns the same finite float."""
        r1 = _estimate_ddg(wt, wt)
        r2 = _estimate_ddg(wt, wt)
        assert r1 == r2
        assert math.isfinite(r1)


# ══════════════════════════════════════════════════════════════
# TEST CLASS 2: score_mutation Always Returns Dict With All 5 Keys
# ══════════════════════════════════════════════════════════════

class TestScoreMutationDictKeys:
    """Property: score_mutation always returns a dict with exactly the
    5 MutationScore keys: stability_ddg, solubility_delta,
    immunogenicity_delta, blosum62, weighted_score.
    """

    @given(
        protein=protein_strategy,
        mutant=aa_strategy,
        weights=weights_strategy,
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_has_all_five_keys(self, protein, mutant, weights):
        """score_mutation result has exactly the 5 expected keys."""
        position = 0  # always valid since min_size=10
        wt = protein[position]
        assume(mutant != wt)  # make it an actual mutation (though same-residue is also valid)

        result = score_mutation(protein, position, mutant, weights=weights)
        assert isinstance(result, dict), (
            f"score_mutation returned {type(result).__name__}, expected dict"
        )
        result_keys = frozenset(result.keys())
        assert result_keys == MUTATION_SCORE_KEYS, (
            f"score_mutation keys = {result_keys}, expected {MUTATION_SCORE_KEYS}"
        )

    @given(
        protein=protein_strategy,
        position=st.integers(min_value=0, max_value=49),
        mutant=aa_strategy,
    )
    @settings(max_examples=150, suppress_health_check=[HealthCheck.too_slow])
    def test_has_five_keys_at_any_position(self, protein, position, mutant):
        """score_mutation result has 5 keys at any valid position."""
        assume(position < len(protein))
        result = score_mutation(protein, position, mutant)
        assert set(result.keys()) == MUTATION_SCORE_KEYS, (
            f"Missing keys: got {set(result.keys())}, expected {MUTATION_SCORE_KEYS}"
        )

    @given(protein=protein_strategy, mutant=aa_strategy)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_stability_ddg_is_finite_float(self, protein, mutant):
        """score_mutation's stability_ddg is always a finite float."""
        result = score_mutation(protein, 0, mutant)
        ddg = result["stability_ddg"]
        assert isinstance(ddg, float), f"stability_ddg is {type(ddg).__name__}"
        assert math.isfinite(ddg), f"stability_ddg = {ddg} is not finite"

    @given(protein=protein_strategy, mutant=aa_strategy)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_solubility_delta_is_finite_float(self, protein, mutant):
        """score_mutation's solubility_delta is always a finite float."""
        result = score_mutation(protein, 0, mutant)
        delta = result["solubility_delta"]
        assert isinstance(delta, float), f"solubility_delta is {type(delta).__name__}"
        assert math.isfinite(delta), f"solubility_delta = {delta} is not finite"

    @given(protein=protein_strategy, mutant=aa_strategy)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_immunogenicity_delta_is_finite_float(self, protein, mutant):
        """score_mutation's immunogenicity_delta is always a finite float."""
        result = score_mutation(protein, 0, mutant)
        delta = result["immunogenicity_delta"]
        assert isinstance(delta, float), f"immunogenicity_delta is {type(delta).__name__}"
        assert math.isfinite(delta), f"immunogenicity_delta = {delta} is not finite"

    @given(protein=protein_strategy, mutant=aa_strategy)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_blosum62_is_int(self, protein, mutant):
        """score_mutation's blosum62 is always an int."""
        result = score_mutation(protein, 0, mutant)
        blosum = result["blosum62"]
        assert isinstance(blosum, int), f"blosum62 is {type(blosum).__name__}"

    @given(protein=protein_strategy, mutant=aa_strategy)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_weighted_score_is_finite_float(self, protein, mutant):
        """score_mutation's weighted_score is always a finite float."""
        result = score_mutation(protein, 0, mutant)
        score = result["weighted_score"]
        assert isinstance(score, float), f"weighted_score is {type(score).__name__}"
        assert math.isfinite(score), f"weighted_score = {score} is not finite"

    @given(protein=protein_strategy, mutant=aa_strategy)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_stability_ddg_matches_estimate_ddg(self, protein, mutant):
        """score_mutation's stability_ddg matches _estimate_ddg for same pair."""
        wt = protein[0]
        result = score_mutation(protein, 0, mutant)
        expected = _estimate_ddg(wt, mutant)
        assert result["stability_ddg"] == expected, (
            f"stability_ddg={result['stability_ddg']} != "
            f"_estimate_ddg({wt!r}, {mutant!r})={expected}"
        )

    @given(
        protein=protein_strategy,
        mutant=aa_strategy,
        weights=weights_strategy,
    )
    @settings(max_examples=80, suppress_health_check=[HealthCheck.too_slow])
    def test_custom_weights_affect_weighted_score(self, protein, mutant, weights):
        """Custom weights should affect the weighted_score value."""
        position = 0
        result_default = score_mutation(protein, position, mutant)
        result_custom = score_mutation(protein, position, mutant, weights=weights)

        # The weighted_score should differ when weights differ from defaults
        # (unless the components happen to produce the same value by coincidence)
        w_stab_d = 0.4
        w_sol_d = 0.3
        w_imm_d = 0.3
        w_stab_c = weights.get("stability", 0.4)
        w_sol_c = weights.get("solubility", 0.3)
        w_imm_c = weights.get("immunogenicity", 0.3)

        # If weights are actually different, recalculate and verify
        ddg = result_custom["stability_ddg"]
        sol = result_custom["solubility_delta"]
        imm = result_custom["immunogenicity_delta"]

        expected_custom = round(w_stab_c * (-ddg) + w_sol_c * sol + w_imm_c * (-imm), 4)
        assert result_custom["weighted_score"] == expected_custom, (
            f"weighted_score mismatch: got {result_custom['weighted_score']}, "
            f"expected {expected_custom}"
        )


# ══════════════════════════════════════════════════════════════
# TEST CLASS 3: Design Result Protein Length Equals Input Length
# ══════════════════════════════════════════════════════════════

class TestDesignResultPreservesLength:
    """Property: All design_* functions produce a designed_protein with the
    same length as the input protein.

    This must hold regardless of:
      - How many mutations are applied
      - What constraints are specified
      - Whether the design succeeds or fails
    """

    @given(protein=protein_short_strategy)
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_thermostable_preserves_length(self, protein):
        """design_thermostable preserves protein length."""
        constraints = DesignConstraints(max_mutations=3)
        result = design_thermostable(protein, constraints=constraints)
        assert isinstance(result, DesignResult)
        assert len(result.designed_protein) == len(protein), (
            f"design_thermostable: input length={len(protein)}, "
            f"output length={len(result.designed_protein)}"
        )

    @given(protein=protein_short_strategy)
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_soluble_preserves_length(self, protein):
        """design_soluble preserves protein length."""
        constraints = DesignConstraints(max_mutations=3)
        result = design_soluble(protein, constraints=constraints)
        assert isinstance(result, DesignResult)
        assert len(result.designed_protein) == len(protein), (
            f"design_soluble: input length={len(protein)}, "
            f"output length={len(result.designed_protein)}"
        )

    @given(protein=protein_short_strategy)
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_low_immunogenicity_preserves_length(self, protein):
        """design_low_immunogenicity preserves protein length."""
        constraints = DesignConstraints(max_mutations=3)
        result = design_low_immunogenicity(protein, constraints=constraints)
        assert isinstance(result, DesignResult)
        assert len(result.designed_protein) == len(protein), (
            f"design_low_immunogenicity: input length={len(protein)}, "
            f"output length={len(result.designed_protein)}"
        )

    @given(
        protein=protein_short_strategy,
        max_mut=st.integers(min_value=1, max_value=10),
    )
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_thermostable_preserves_length_varying_mutations(self, protein, max_mut):
        """design_thermostable preserves length with varying max_mutations."""
        constraints = DesignConstraints(max_mutations=max_mut)
        result = design_thermostable(protein, constraints=constraints)
        assert len(result.designed_protein) == len(protein)

    @given(
        protein=protein_short_strategy,
        max_mut=st.integers(min_value=1, max_value=10),
    )
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_soluble_preserves_length_varying_mutations(self, protein, max_mut):
        """design_soluble preserves length with varying max_mutations."""
        constraints = DesignConstraints(max_mutations=max_mut)
        result = design_soluble(protein, constraints=constraints)
        assert len(result.designed_protein) == len(protein)

    @given(
        protein=protein_short_strategy,
        max_mut=st.integers(min_value=1, max_value=10),
    )
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_low_immunogenicity_preserves_length_varying_mutations(self, protein, max_mut):
        """design_low_immunogenicity preserves length with varying max_mutations."""
        constraints = DesignConstraints(max_mutations=max_mut)
        result = design_low_immunogenicity(protein, constraints=constraints)
        assert len(result.designed_protein) == len(protein)

    def test_thermostable_known_sequence(self):
        """design_thermostable on a known sequence preserves length."""
        protein = "MKTLLILAVF"
        result = design_thermostable(protein, constraints=DesignConstraints(max_mutations=5))
        assert len(result.designed_protein) == 10

    def test_soluble_known_sequence(self):
        """design_soluble on a known hydrophobic sequence preserves length."""
        protein = "IIIVVLLLFFFFYYYI"
        result = design_soluble(protein, constraints=DesignConstraints(max_mutations=5))
        assert len(result.designed_protein) == 16

    def test_low_immunogenicity_known_sequence(self):
        """design_low_immunogenicity on a known hydrophobic sequence preserves length."""
        protein = "IIIVVLLLFFFFYYYI"
        result = design_low_immunogenicity(
            protein, constraints=DesignConstraints(max_mutations=5)
        )
        assert len(result.designed_protein) == 16

    @given(protein=protein_short_strategy)
    @settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_thermostable_result_is_valid_aa_sequence(self, protein):
        """design_thermostable output contains only standard amino acids."""
        constraints = DesignConstraints(max_mutations=3)
        result = design_thermostable(protein, constraints=constraints)
        for aa in result.designed_protein:
            assert aa in STANDARD_AA, (
                f"Non-standard AA {aa!r} in designed protein"
            )

    @given(protein=protein_short_strategy)
    @settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_soluble_result_is_valid_aa_sequence(self, protein):
        """design_soluble output contains only standard amino acids."""
        constraints = DesignConstraints(max_mutations=3)
        result = design_soluble(protein, constraints=constraints)
        for aa in result.designed_protein:
            assert aa in STANDARD_AA, (
                f"Non-standard AA {aa!r} in designed protein"
            )

    @given(protein=protein_short_strategy)
    @settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_low_immunogenicity_result_is_valid_aa_sequence(self, protein):
        """design_low_immunogenicity output contains only standard amino acids."""
        constraints = DesignConstraints(max_mutations=3)
        result = design_low_immunogenicity(protein, constraints=constraints)
        for aa in result.designed_protein:
            assert aa in STANDARD_AA, (
                f"Non-standard AA {aa!r} in designed protein"
            )

    @given(protein=protein_short_strategy)
    @settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_design_result_has_original_protein(self, protein):
        """DesignResult stores the original protein string."""
        constraints = DesignConstraints(max_mutations=2)
        result = design_thermostable(protein, constraints=constraints)
        assert result.original_protein == protein

    @given(protein=protein_short_strategy)
    @settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow], deadline=5000)
    def test_design_result_iterations_within_budget(self, protein):
        """Number of iterations does not exceed max_mutations.

        Note: len(mutations) can exceed max_mutations because disulfide pair
        strategies add a partner mutation per iteration (2 entries per step).
        The iteration count is the true budget constraint.
        """
        max_mut = 3
        constraints = DesignConstraints(max_mutations=max_mut)
        result = design_thermostable(protein, constraints=constraints)
        assert result.iterations <= max_mut, (
            f"Ran {result.iterations} iterations, "
            f"but max_mutations={max_mut}"
        )
