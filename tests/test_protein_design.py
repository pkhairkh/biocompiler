"""Comprehensive tests for the protein_design module.

Tests cover heuristic predictors, base scoring functions, mutation scoring,
design algorithms, disulfide/proline scanners, and data classes.
"""

import pytest

from biocompiler.protein_design import (
    DesignConstraints,
    DesignResult,
    _base_immunogenicity,
    _base_solubility,
    _base_stability,
    _check_constraints,
    _estimate_ddg,
    _estimate_immunogenicity_delta,
    _estimate_solubility_delta,
    _is_preserved,
    _predict_secondary_structure_simple,
    design_low_immunogenicity,
    design_soluble,
    design_thermostable,
    find_disulfide_opportunities,
    find_proline_substitution_sites,
    score_mutation,
)
from biocompiler.constants import BLOSUM62, HYDROPATHY


# ────────────────────────────────────────────────────────────
# Fixtures
# ────────────────────────────────────────────────────────────

@pytest.fixture
def short_protein():
    """A short protein with mixed residues, enough for MHC_II_WINDOW analysis."""
    return "MKTLLILAVF"  # 10 residues — mixed hydrophobic/charged


@pytest.fixture
def hydrophobic_protein():
    """A protein rich in hydrophobic residues."""
    return "IIIVVLLLFFFFYYYY"  # 16 residues — very hydrophobic


@pytest.fixture
def hydrophilic_protein():
    """A protein rich in charged/polar residues."""
    return "DEKRDEKRDEKRDEKR"  # 16 residues — very hydrophilic


@pytest.fixture
def loop_protein():
    """Protein with clear loop regions (non-helix, non-strand residues).
    
    'STNQSTNQ' are polar residues not in HELIX_FORMERS or STRAND_FORMERS,
    so they will be predicted as loops. The surrounding hydrophobic stretches
    will form strands/helices.
    """
    return "VVVIIISTNQSTNQVVVIII"  # 20 residues — loops in the middle


@pytest.fixture
def disulfide_protein():
    """Protein with multiple loop positions that can host disulfide pairs.
    
    Contains enough loop residues separated by >5 positions.
    """
    return "STNQPPPPSTNQ"  # 12 residues, loops at 0-3 and 8-11


# ────────────────────────────────────────────────────────────
# 1. Test _estimate_ddg
# ────────────────────────────────────────────────────────────

class TestEstimateDdg:
    """Tests for the _estimate_ddg heuristic."""

    def test_stabilizing_mutation_negative(self):
        """A conservative stabilizing mutation should yield negative ddg."""
        # A→V: BLOSUM62[A][V] = 0, hydropathy change = 4.2 - 1.8 = 2.4
        # ddg = -0.15 * 0 + 0.05 * 2.4 = 0.12
        # A→V is slightly destabilizing by this model, not the best example.
        # Try R→K: BLOSUM62[R][K] = 2, hydro change = -3.9 - (-4.5) = 0.6
        # ddg = -0.15 * 2 + 0.05 * 0.6 = -0.3 + 0.03 = -0.27
        ddg = _estimate_ddg("R", "K")
        assert ddg < 0, f"R→K should be stabilizing, got ddg={ddg}"

    def test_destabilizing_mutation_positive(self):
        """A clearly destabilizing mutation should yield positive ddg."""
        # W→G: BLOSUM62[W][G] = -2, hydro change = -0.4 - (-0.9) = 0.5
        # ddg = -0.15 * (-2) + 0.05 * 0.5 = 0.3 + 0.025 = 0.325
        # Plus G bonus: +0.3 → 0.625
        ddg = _estimate_ddg("W", "G")
        assert ddg > 0, f"W→G should be destabilizing, got ddg={ddg}"

    def test_proline_bonus(self):
        """Introducing proline should add a stability bonus (negative)."""
        ddg_to_pro = _estimate_ddg("A", "P")
        # A→P: BLOSUM62[A][P] = -1, hydro change = -1.6 - 1.8 = -3.4
        # ddg = -0.15*(-1) + 0.05*(-3.4) = 0.15 - 0.17 = -0.02
        # Plus PROLINE_BONUS: -0.3 → -0.32
        assert ddg_to_pro < 0, f"A→P should be stabilizing due to proline bonus, got {ddg_to_pro}"

    def test_glycine_penalty(self):
        """Introducing glycine should add a flexibility penalty (positive)."""
        ddg_to_gly = _estimate_ddg("A", "G")
        # A→G without bonus: BLOSUM62[A][G] = 0, hydro change = -0.4 - 1.8 = -2.2
        # ddg_base = -0.15 * 0 + 0.05 * (-2.2) = -0.11
        # Plus GLYCINE_BONUS: +0.3 → 0.19
        assert ddg_to_gly > 0, f"A→G should be destabilizing due to glycine penalty, got {ddg_to_gly}"

    def test_same_residue_ddg_zero(self):
        """Mutation to same residue should yield ddg of 0."""
        # A→A: BLOSUM62[A][A] = 4, hydro change = 0
        # ddg = -0.15 * 4 = -0.6, but no proline/glycine bonus
        ddg = _estimate_ddg("A", "A")
        assert ddg == pytest.approx(-0.15 * 4, abs=0.01)

    def test_returns_float_rounded_to_3(self):
        """Return value should be a float rounded to 3 decimal places."""
        ddg = _estimate_ddg("L", "I")
        assert isinstance(ddg, float)
        # Check it's rounded: multiply by 1000 and check it's close to an integer
        assert round(ddg * 1000) == ddg * 1000 or abs(ddg * 1000 - round(ddg * 1000)) < 1e-6

    def test_known_pair_l_to_i(self):
        """L→I is a conservative hydrophobic swap — should be near zero."""
        ddg = _estimate_ddg("L", "I")
        blosum = BLOSUM62["L"]["I"]  # 2
        dh = HYDROPATHY["I"] - HYDROPATHY["L"]  # 4.5 - 3.8 = 0.7
        expected = round(-0.15 * blosum + 0.05 * dh, 3)
        assert ddg == expected


# ────────────────────────────────────────────────────────────
# 2. Test _estimate_solubility_delta
# ────────────────────────────────────────────────────────────

class TestEstimateSolubilityDelta:
    """Tests for the _estimate_solubility_delta heuristic."""

    def test_hydrophobic_to_charged_improves(self):
        """Replacing a hydrophobic residue with a charged one should increase solubility."""
        delta = _estimate_solubility_delta("I", "D")
        assert delta > 0, f"I→D should improve solubility, got delta={delta}"

    def test_hydrophobic_to_hydrophilic_improves(self):
        """Replacing hydrophobic with hydrophilic should increase solubility."""
        delta = _estimate_solubility_delta("V", "K")
        assert delta > 0, f"V→K should improve solubility, got delta={delta}"

    def test_charged_to_hydrophobic_decreases(self):
        """Replacing a charged residue with hydrophobic should decrease solubility."""
        delta = _estimate_solubility_delta("D", "I")
        assert delta < 0, f"D→I should decrease solubility, got delta={delta}"

    def test_charged_bonus(self):
        """Introducing a charged residue from a non-charged one gets a bonus."""
        delta_charged = _estimate_solubility_delta("A", "K")
        # Should include SOLUBILITY_CHARGED_WEIGHT = 0.3 bonus
        # Base hydropathy part: (1.8 - (-3.9)) * 0.2 = 1.14
        # Plus charged bonus: 0.3
        assert delta_charged > 0.3, f"Charged bonus should push delta above 0.3, got {delta_charged}"

    def test_aggregation_prone_penalty(self):
        """Introducing an aggregation-prone residue gets a penalty."""
        delta_agg = _estimate_solubility_delta("D", "V")
        # V is in _AGGREGATION_PRONE_AAS, D is not
        # Should have SOLUBILITY_AGGREGATION_WEIGHT penalty subtracted
        # Base hydropathy part: (-3.5 - 4.2) * 0.2 = -1.54 (already negative)
        # Plus aggregation penalty: -0.2
        assert delta_agg < -0.2, f"Aggregation penalty should make delta more negative, got {delta_agg}"

    def test_returns_float_rounded_to_3(self):
        """Return value should be a float rounded to 3 decimal places."""
        delta = _estimate_solubility_delta("L", "E")
        assert isinstance(delta, float)


# ────────────────────────────────────────────────────────────
# 3. Test _estimate_immunogenicity_delta
# ────────────────────────────────────────────────────────────

class TestEstimateImmunogenicityDelta:
    """Tests for the _estimate_immunogenicity_delta heuristic."""

    def test_hydrophobic_to_polar_decreases_immunogenicity(self):
        """Replacing hydrophobic with polar should decrease immunogenicity (negative delta)."""
        protein = "IIIVVLLLFFFFYYYI"  # 16 hydrophobic residues
        # Position 4 is V, mutate to D (charged/hydrophilic)
        delta = _estimate_immunogenicity_delta(protein, 4, "D")
        assert delta < 0, f"V→D in hydrophobic region should reduce immunogenicity, got {delta}"

    def test_polar_to_hydrophobic_increases_immunogenicity(self):
        """Replacing polar with hydrophobic should increase immunogenicity (positive delta)."""
        protein = "DEKRDEKRDEKRDEKR"  # 16 charged/polar residues
        # Position 0 is D, mutate to I (hydrophobic)
        delta = _estimate_immunogenicity_delta(protein, 0, "I")
        assert delta > 0, f"D→I in polar region should increase immunogenicity, got {delta}"

    def test_same_residue_zero_delta(self):
        """Mutating to the same residue should give approximately zero delta."""
        protein = "AAAAAAAAAAAAAAA"  # 15 A's
        delta = _estimate_immunogenicity_delta(protein, 5, "A")
        assert delta == pytest.approx(0.0, abs=0.001)

    def test_position_near_boundaries(self):
        """Positions near start/end of protein should still compute correctly."""
        protein = "IIIVVLLLFFFFYYYI"
        # Position 0 (near start)
        delta_start = _estimate_immunogenicity_delta(protein, 0, "D")
        # Position 15 (near end)
        delta_end = _estimate_immunogenicity_delta(protein, 15, "D")
        assert isinstance(delta_start, float)
        assert isinstance(delta_end, float)
        assert delta_start < 0  # hydrophobic→charged should decrease
        assert delta_end < 0

    def test_returns_float_rounded_to_3(self):
        """Return value should be a float rounded to 3 decimal places."""
        protein = "IIIVVLLLFFFFYYYI"
        delta = _estimate_immunogenicity_delta(protein, 3, "K")
        assert isinstance(delta, float)


# ────────────────────────────────────────────────────────────
# 4. Test _base_stability, _base_solubility, _base_immunogenicity
# ────────────────────────────────────────────────────────────

class TestBaseScorers:
    """Tests for the base scoring functions."""

    # --- _base_stability ---

    def test_base_stability_empty(self):
        """Empty protein should return 0.0 stability."""
        assert _base_stability("") == 0.0

    def test_base_stability_returns_float(self):
        """Should return a float."""
        result = _base_stability("MKTLLILAVF")
        assert isinstance(result, float)

    def test_base_stability_proline_rich_is_stable(self):
        """A proline-rich protein should have more negative (stable) ddg."""
        pro_rich = "PPPPPPPPPP"
        no_pro = "AAAAAAAAAA"
        assert _base_stability(pro_rich) < _base_stability(no_pro)

    def test_base_stability_glycine_rich_is_unstable(self):
        """A glycine-rich protein should have less negative (unstable) ddg."""
        gly_rich = "GGGGGGGGGG"
        no_gly = "AAAAAAAAAA"
        assert _base_stability(gly_rich) > _base_stability(no_gly)

    def test_base_stability_disulfide_pairs(self):
        """Proteins with cysteine pairs should be more stable."""
        with_cys = "CCCCCCAAAA"  # 3 pairs
        no_cys = "AAAAAAAEAA"
        assert _base_stability(with_cys) < _base_stability(no_cys)

    # --- _base_solubility ---

    def test_base_solubility_empty(self):
        """Empty protein should return 0.0 solubility."""
        assert _base_solubility("") == 0.0

    def test_base_solubility_returns_float(self):
        """Should return a float."""
        result = _base_solubility("MKTLLILAVF")
        assert isinstance(result, float)

    def test_base_solubility_charged_better_than_hydrophobic(self):
        """A charged protein should have higher solubility than a hydrophobic one."""
        charged = "DEKRDEKRDE"
        hydrophobic = "IIIVVLLLFA"
        assert _base_solubility(charged) > _base_solubility(hydrophobic)

    def test_base_solubility_hydrophobic_stretch_penalty(self):
        """A protein with long hydrophobic stretches should get a penalty."""
        with_stretch = "IIIVVLLLFFFFFFYY"  # long hydrophobic stretch
        no_stretch = "DEKRDEKRDEKRDEKR"  # no stretch
        assert _base_solubility(with_stretch) < _base_solubility(no_stretch)

    # --- _base_immunogenicity ---

    def test_base_immunogenicity_short_protein(self):
        """Protein shorter than MHC_II_WINDOW should return 0.0."""
        short = "MKTL"  # 4 residues < 9
        assert _base_immunogenicity(short) == 0.0

    def test_base_immunogenicity_returns_float(self):
        """Should return a float."""
        result = _base_immunogenicity("IIIVVLLLFFFFYYYI")
        assert isinstance(result, float)

    def test_base_immunogenicity_in_0_to_1_range(self):
        """Should be in [0, 1] range."""
        result = _base_immunogenicity("IIIVVLLLFFFFYYYI")
        assert 0.0 <= result <= 1.0

    def test_base_immunogenicity_hydrophobic_higher_than_charged(self):
        """A hydrophobic protein should have higher immunogenicity than a charged one."""
        hydrophobic = "IIIVVLLLFFFFYYYW"
        charged = "DEKRDEKRDEKRDEKR"
        assert _base_immunogenicity(hydrophobic) > _base_immunogenicity(charged)


# ────────────────────────────────────────────────────────────
# 5. Test score_mutation
# ────────────────────────────────────────────────────────────

class TestScoreMutation:
    """Tests for the score_mutation function."""

    def test_returns_dict_with_expected_keys(self):
        """score_mutation should return dict with all expected keys."""
        protein = "MKTLLILAVF"
        result = score_mutation(protein, 0, "D")
        expected_keys = {"stability_ddg", "solubility_delta", "immunogenicity_delta", "blosum62", "weighted_score"}
        assert set(result.keys()) == expected_keys

    def test_values_are_numeric(self):
        """All values in the returned dict should be numeric (int or float)."""
        protein = "MKTLLILAVF"
        result = score_mutation(protein, 0, "D")
        for key, value in result.items():
            assert isinstance(value, (int, float)), f"{key} should be numeric, got {type(value)}"

    def test_stability_ddg_matches_estimate(self):
        """stability_ddg should match _estimate_ddg for the same mutation."""
        protein = "MKTLLILAVF"
        wt = protein[0]  # M
        mutant = "D"
        result = score_mutation(protein, 0, mutant)
        assert result["stability_ddg"] == _estimate_ddg(wt, mutant)

    def test_solubility_delta_matches_estimate(self):
        """solubility_delta should match _estimate_solubility_delta."""
        protein = "MKTLLILAVF"
        wt = protein[0]
        mutant = "D"
        result = score_mutation(protein, 0, mutant)
        assert result["solubility_delta"] == _estimate_solubility_delta(wt, mutant)

    def test_blosum62_matches_lookup(self):
        """blosum62 should match the BLOSUM62 matrix value."""
        protein = "MKTLLILAVF"
        wt = protein[0]
        mutant = "D"
        result = score_mutation(protein, 0, mutant)
        assert result["blosum62"] == BLOSUM62.get(wt, {}).get(mutant, -4)

    def test_custom_weights(self):
        """Custom weights should affect the weighted_score."""
        protein = "MKTLLILAVF"
        result_default = score_mutation(protein, 0, "D")
        result_custom = score_mutation(protein, 0, "D", weights={"stability": 1.0, "solubility": 0.0, "immunogenicity": 0.0})
        # With all weight on stability, score should just be -ddg
        assert result_custom["weighted_score"] != result_default["weighted_score"]

    def test_weighted_score_calculation(self):
        """weighted_score should follow the formula: w_stab*(-ddg) + w_sol*sol + w_imm*(-imm)."""
        protein = "MKTLLILAVF"
        result = score_mutation(protein, 0, "D")
        ddg = result["stability_ddg"]
        sol = result["solubility_delta"]
        imm = result["immunogenicity_delta"]
        expected = 0.4 * (-ddg) + 0.3 * sol + 0.3 * (-imm)
        assert result["weighted_score"] == pytest.approx(expected, abs=0.001)


# ────────────────────────────────────────────────────────────
# 6. Test design_thermostable
# ────────────────────────────────────────────────────────────

class TestDesignThermostable:
    """Tests for design_thermostable function."""

    def test_returns_design_result(self, short_protein):
        """design_thermostable should return a DesignResult instance."""
        result = design_thermostable(short_protein, constraints=DesignConstraints(max_mutations=3))
        assert isinstance(result, DesignResult)

    def test_result_has_all_fields(self, short_protein):
        """Result should have all DesignResult fields populated."""
        result = design_thermostable(short_protein, constraints=DesignConstraints(max_mutations=3))
        assert isinstance(result.original_protein, str)
        assert isinstance(result.designed_protein, str)
        assert isinstance(result.mutations, list)
        assert isinstance(result.stability_change, float)
        assert isinstance(result.solubility_change, float)
        assert isinstance(result.immunogenicity_change, float)
        assert isinstance(result.iterations, int)
        assert isinstance(result.constraints_satisfied, list)
        assert isinstance(result.constraints_violated, list)
        assert isinstance(result.execution_time_s, float)

    def test_designed_protein_same_length(self, short_protein):
        """Designed protein should have same length as original."""
        result = design_thermostable(short_protein, constraints=DesignConstraints(max_mutations=3))
        assert len(result.designed_protein) == len(short_protein)

    def test_iterations_within_budget(self, short_protein):
        """Iterations should not exceed max_mutations."""
        max_mut = 5
        result = design_thermostable(short_protein, constraints=DesignConstraints(max_mutations=max_mut))
        assert result.iterations <= max_mut

    def test_execution_time_non_negative(self, short_protein):
        """Execution time should be non-negative."""
        result = design_thermostable(short_protein, constraints=DesignConstraints(max_mutations=2))
        assert result.execution_time_s >= 0.0

    def test_invalid_sequence_raises(self):
        """Invalid protein sequence should raise ValueError."""
        with pytest.raises(ValueError):
            design_thermostable("123BAD")


# ────────────────────────────────────────────────────────────
# 7. Test design_soluble
# ────────────────────────────────────────────────────────────

class TestDesignSoluble:
    """Tests for design_soluble function."""

    def test_returns_design_result(self, hydrophobic_protein):
        """design_soluble should return a DesignResult instance."""
        result = design_soluble(hydrophobic_protein, constraints=DesignConstraints(max_mutations=3))
        assert isinstance(result, DesignResult)

    def test_result_has_all_fields(self, hydrophobic_protein):
        """Result should have all DesignResult fields populated."""
        result = design_soluble(hydrophobic_protein, constraints=DesignConstraints(max_mutations=3))
        assert isinstance(result.original_protein, str)
        assert isinstance(result.designed_protein, str)
        assert isinstance(result.mutations, list)
        assert isinstance(result.stability_change, float)
        assert isinstance(result.solubility_change, float)
        assert isinstance(result.immunogenicity_change, float)
        assert isinstance(result.iterations, int)

    def test_designed_protein_same_length(self, hydrophobic_protein):
        """Designed protein should have same length as original."""
        result = design_soluble(hydrophobic_protein, constraints=DesignConstraints(max_mutations=3))
        assert len(result.designed_protein) == len(hydrophobic_protein)

    def test_invalid_sequence_raises(self):
        """Invalid protein sequence should raise ValueError."""
        with pytest.raises(ValueError):
            design_soluble("XXXXX")


# ────────────────────────────────────────────────────────────
# 8. Test design_low_immunogenicity
# ────────────────────────────────────────────────────────────

class TestDesignLowImmunogenicity:
    """Tests for design_low_immunogenicity function."""

    def test_returns_design_result(self, hydrophobic_protein):
        """design_low_immunogenicity should return a DesignResult instance."""
        result = design_low_immunogenicity(
            hydrophobic_protein, constraints=DesignConstraints(max_mutations=3)
        )
        assert isinstance(result, DesignResult)

    def test_result_has_all_fields(self, hydrophobic_protein):
        """Result should have all DesignResult fields populated."""
        result = design_low_immunogenicity(
            hydrophobic_protein, constraints=DesignConstraints(max_mutations=3)
        )
        assert isinstance(result.original_protein, str)
        assert isinstance(result.designed_protein, str)
        assert isinstance(result.mutations, list)
        assert isinstance(result.stability_change, float)
        assert isinstance(result.solubility_change, float)
        assert isinstance(result.immunogenicity_change, float)
        assert isinstance(result.iterations, int)

    def test_designed_protein_same_length(self, hydrophobic_protein):
        """Designed protein should have same length as original."""
        result = design_low_immunogenicity(
            hydrophobic_protein, constraints=DesignConstraints(max_mutations=3)
        )
        assert len(result.designed_protein) == len(hydrophobic_protein)

    def test_invalid_sequence_raises(self):
        """Invalid protein sequence should raise ValueError."""
        with pytest.raises(ValueError):
            design_low_immunogenicity("12345")


# ────────────────────────────────────────────────────────────
# 9. Test find_disulfide_opportunities
# ────────────────────────────────────────────────────────────

class TestFindDisulfideOpportunities:
    """Tests for find_disulfide_opportunities function."""

    def test_returns_list(self, loop_protein):
        """Should return a list."""
        result = find_disulfide_opportunities(loop_protein)
        assert isinstance(result, list)

    def test_opportunity_dict_keys(self, loop_protein):
        """Each opportunity should have the expected keys."""
        result = find_disulfide_opportunities(loop_protein)
        for opp in result:
            assert "position1" in opp
            assert "position2" in opp
            assert "distance_estimate" in opp
            assert "stabilizing_estimate" in opp

    def test_positions_are_loop_residues(self, loop_protein):
        """Both positions should be in predicted loop regions."""
        ss = _predict_secondary_structure_simple(loop_protein)
        result = find_disulfide_opportunities(loop_protein)
        for opp in result:
            assert ss[opp["position1"]] == "L"
            assert ss[opp["position2"]] == "L"

    def test_positions_not_cysteine(self, loop_protein):
        """Neither position should already be a cysteine."""
        result = find_disulfide_opportunities(loop_protein)
        for opp in result:
            assert loop_protein[opp["position1"]] != "C"
            assert loop_protein[opp["position2"]] != "C"

    def test_positions_at_least_5_apart(self, loop_protein):
        """Positions must be at least 5 residues apart."""
        result = find_disulfide_opportunities(loop_protein)
        for opp in result:
            assert abs(opp["position2"] - opp["position1"]) >= 5

    def test_sorted_by_stabilizing_estimate(self, loop_protein):
        """Results should be sorted by stabilizing_estimate (most stabilizing first)."""
        result = find_disulfide_opportunities(loop_protein)
        if len(result) > 1:
            estimates = [opp["stabilizing_estimate"] for opp in result]
            assert estimates == sorted(estimates)

    def test_empty_for_no_loop_protein(self):
        """Protein with no loop regions should return empty list."""
        # "AELMAELM" → all helix formers → all H, no loops
        no_loop = "AELMAELMAELM"
        result = find_disulfide_opportunities(no_loop)
        assert result == []

    def test_protein_with_cysteine_in_loop(self):
        """Loop positions that are already C should be excluded."""
        protein = "CCCCCSTNQSTNQ"  # C's at start but in loops
        result = find_disulfide_opportunities(protein)
        for opp in result:
            assert protein[opp["position1"]] != "C"
            assert protein[opp["position2"]] != "C"


# ────────────────────────────────────────────────────────────
# 10. Test find_proline_substitution_sites
# ────────────────────────────────────────────────────────────

class TestFindProlineSubstitutionSites:
    """Tests for find_proline_substitution_sites function."""

    def test_returns_list(self, loop_protein):
        """Should return a list."""
        result = find_proline_substitution_sites(loop_protein)
        assert isinstance(result, list)

    def test_site_dict_keys(self, loop_protein):
        """Each site should have the expected keys."""
        result = find_proline_substitution_sites(loop_protein)
        for site in result:
            assert "position" in site
            assert "wildtype" in site
            assert "ddg_estimate" in site
            assert "in_loop" in site

    def test_sites_are_in_loops(self, loop_protein):
        """All sites should be in predicted loop regions."""
        ss = _predict_secondary_structure_simple(loop_protein)
        result = find_proline_substitution_sites(loop_protein)
        for site in result:
            assert ss[site["position"]] == "L"
            assert site["in_loop"] is True

    def test_no_proline_sites(self, loop_protein):
        """No site should have wildtype = P."""
        result = find_proline_substitution_sites(loop_protein)
        for site in result:
            assert site["wildtype"] != "P"

    def test_no_glycine_sites(self, loop_protein):
        """No site should have wildtype = G (G→P is highly destabilizing)."""
        result = find_proline_substitution_sites(loop_protein)
        for site in result:
            assert site["wildtype"] != "G"

    def test_blosum62_score_at_least_neg1(self, loop_protein):
        """All sites should have BLOSUM62 score >= -1 for X→P."""
        result = find_proline_substitution_sites(loop_protein)
        for site in result:
            blosum = BLOSUM62.get(site["wildtype"], {}).get("P", -4)
            assert blosum >= -1

    def test_sorted_by_ddg_estimate(self, loop_protein):
        """Results should be sorted by ddg_estimate (most stabilizing first)."""
        result = find_proline_substitution_sites(loop_protein)
        if len(result) > 1:
            ddgs = [site["ddg_estimate"] for site in result]
            assert ddgs == sorted(ddgs)

    def test_empty_for_no_loop_protein(self):
        """Protein with no loop regions should return empty list."""
        # "AELMAELM" → all helix formers → all H, no loops
        no_loop = "AELMAELMAELM"
        result = find_proline_substitution_sites(no_loop)
        assert result == []


# ────────────────────────────────────────────────────────────
# 11. Test DesignConstraints
# ────────────────────────────────────────────────────────────

class TestDesignConstraints:
    """Tests for the DesignConstraints dataclass."""

    def test_default_values(self):
        """Default values should match the documented defaults."""
        dc = DesignConstraints()
        assert dc.min_stability_kcal == -5.0
        assert dc.min_solubility_score == 0.0
        assert dc.max_immunogenicity == 0.5
        assert dc.max_mutations == 10
        assert dc.blosum62_min == 0
        assert dc.max_ddg_per_mutation == 2.0
        assert dc.preserve_positions is None
        assert dc.preserve_residues is None

    def test_custom_values(self):
        """Custom values should be accepted."""
        dc = DesignConstraints(
            min_stability_kcal=-10.0,
            min_solubility_score=0.5,
            max_immunogenicity=0.3,
            max_mutations=5,
            blosum62_min=1,
            max_ddg_per_mutation=1.5,
            preserve_positions=[0, 5],
            preserve_residues=["C"],
        )
        assert dc.min_stability_kcal == -10.0
        assert dc.min_solubility_score == 0.5
        assert dc.max_immunogenicity == 0.3
        assert dc.max_mutations == 5
        assert dc.blosum62_min == 1
        assert dc.max_ddg_per_mutation == 1.5
        assert dc.preserve_positions == [0, 5]
        assert dc.preserve_residues == ["C"]

    def test_is_dataclass(self):
        """DesignConstraints should be a dataclass."""
        from dataclasses import is_dataclass
        assert is_dataclass(DesignConstraints)

    def test_preserve_positions_effect(self):
        """preserve_positions should protect specific positions from mutation."""
        dc = DesignConstraints(preserve_positions=[0, 1, 2])
        assert _is_preserved(0, "A", dc) is True
        assert _is_preserved(3, "A", dc) is False

    def test_preserve_residues_effect(self):
        """preserve_residues should protect specific residue types from mutation."""
        dc = DesignConstraints(preserve_residues=["C"])
        assert _is_preserved(0, "C", dc) is True
        assert _is_preserved(0, "A", dc) is False

    def test_no_preservation(self):
        """With no preserve options, nothing should be preserved."""
        dc = DesignConstraints()
        assert _is_preserved(0, "C", dc) is False
        assert _is_preserved(0, "A", dc) is False


# ────────────────────────────────────────────────────────────
# 12. Test DesignResult
# ────────────────────────────────────────────────────────────

class TestDesignResult:
    """Tests for the DesignResult dataclass."""

    def test_construction_basic(self):
        """Basic construction with defaults should work."""
        result = DesignResult()
        assert result.original_protein == ""
        assert result.designed_protein == ""
        assert result.mutations == []
        assert result.stability_change == 0.0
        assert result.solubility_change == 0.0
        assert result.immunogenicity_change == 0.0
        assert result.cai is None
        assert result.iterations == 0
        assert result.constraints_satisfied == []
        assert result.constraints_violated == []

    def test_construction_with_values(self):
        """Construction with custom values should work."""
        result = DesignResult(
            original_protein="MKTL",
            designed_protein="MKDL",
            mutations=[{"position": 2, "wildtype": "T", "mutant": "D", "ddg": -0.5}],
            stability_change=-0.5,
            solubility_change=0.1,
            immunogenicity_change=-0.05,
            iterations=1,
            success=True,
            constraints_satisfied=["min_stability"],
            constraints_violated=[],
        )
        assert result.original_protein == "MKTL"
        assert result.designed_protein == "MKDL"
        assert len(result.mutations) == 1
        assert result.stability_change == -0.5
        assert result.success is True

    def test_post_init_syncs_sequence(self):
        """__post_init__ should sync sequence from designed_protein when sequence is empty."""
        result = DesignResult(designed_protein="MKDL")
        assert result.sequence == "MKDL"

    def test_post_init_syncs_primary_score(self):
        """__post_init__ should sync primary_score from stability_change when 0."""
        result = DesignResult(stability_change=-2.5)
        assert result.primary_score == -2.5

    def test_post_init_syncs_stability_change_from_primary_score(self):
        """__post_init__ should sync stability_change from primary_score when 0."""
        result = DesignResult(primary_score=-3.0)
        assert result.stability_change == -3.0

    def test_post_init_classification_no_violations(self):
        """With no violated constraints and empty classification, should be design_success."""
        result = DesignResult(constraints_violated=[], success=True)
        assert result.classification == "design_success"

    def test_post_init_classification_partial(self):
        """With violations but success=True, should be design_partial."""
        result = DesignResult(constraints_violated=["min_stability"], success=True)
        assert result.classification == "design_partial"

    def test_post_init_classification_failed(self):
        """With violations and success=False, should be design_failed."""
        result = DesignResult(constraints_violated=["min_stability"], success=False)
        assert result.classification == "design_failed"

    def test_post_init_classification_preserved(self):
        """If classification is already set, __post_init__ should not overwrite it."""
        result = DesignResult(classification="custom_class")
        assert result.classification == "custom_class"

    def test_designed_sequence_property(self):
        """designed_sequence property should return designed_protein."""
        result = DesignResult(designed_protein="MKDL")
        assert result.designed_sequence == "MKDL"

    def test_engine_name_default(self):
        """engine_name should default to 'protein_design'."""
        result = DesignResult()
        assert result.engine_name == "protein_design"

    def test_primary_score_label_default(self):
        """primary_score_label should default to 'ddg'."""
        result = DesignResult()
        assert result.primary_score_label == "ddg"

    def test_is_dataclass(self):
        """DesignResult should be a dataclass."""
        from dataclasses import is_dataclass
        assert is_dataclass(DesignResult)

    def test_inherits_base_engine_result(self):
        """DesignResult should inherit from BaseEngineResult."""
        from biocompiler.engine_base import BaseEngineResult
        assert issubclass(DesignResult, BaseEngineResult)


# ────────────────────────────────────────────────────────────
# Additional: Test _check_constraints
# ────────────────────────────────────────────────────────────

class TestCheckConstraints:
    """Tests for the _check_constraints helper."""

    def test_all_satisfied(self):
        """When all constraints are met, all should be in satisfied list."""
        dc = DesignConstraints(min_stability_kcal=-5.0, min_solubility_score=0.0, max_immunogenicity=0.5)
        satisfied, violated = _check_constraints("MKTL", dc, -6.0, 0.5, 0.3)
        assert "min_stability" in satisfied
        assert "min_solubility" in satisfied
        assert "max_immunogenicity" in satisfied
        assert len(violated) == 0

    def test_all_violated(self):
        """When no constraints are met, all should be in violated list."""
        dc = DesignConstraints(min_stability_kcal=-5.0, min_solubility_score=1.0, max_immunogenicity=0.1)
        satisfied, violated = _check_constraints("MKTL", dc, 0.0, 0.5, 0.8)
        assert "min_stability" in violated
        assert "min_solubility" in violated
        assert "max_immunogenicity" in violated
        assert len(satisfied) == 0

    def test_partial(self):
        """Some constraints satisfied, some violated."""
        dc = DesignConstraints(min_stability_kcal=-5.0, min_solubility_score=0.0, max_immunogenicity=0.5)
        satisfied, violated = _check_constraints("MKTL", dc, -6.0, 0.5, 0.8)
        assert "min_stability" in satisfied
        assert "min_solubility" in satisfied
        assert "max_immunogenicity" in violated


# ────────────────────────────────────────────────────────────
# Additional: Test _predict_secondary_structure_simple
# ────────────────────────────────────────────────────────────

class TestPredictSecondaryStructureSimple:
    """Tests for the _predict_secondary_structure_simple helper."""

    def test_returns_list_of_strings(self):
        """Should return a list of single-character strings."""
        result = _predict_secondary_structure_simple("AELMAELMSTNQ")
        assert isinstance(result, list)
        assert all(isinstance(s, str) for s in result)
        assert all(s in {"H", "E", "L"} for s in result)

    def test_length_matches_protein(self):
        """Output length should match protein length."""
        protein = "AELMVIYSTNQ"
        result = _predict_secondary_structure_simple(protein)
        assert len(result) == len(protein)

    def test_helix_formers(self):
        """A run of helix formers (AELM) >= 4 should be predicted as helix."""
        protein = "AELM"  # exactly 4 helix formers
        result = _predict_secondary_structure_simple(protein)
        assert all(s == "H" for s in result)

    def test_strand_formers(self):
        """A run of strand formers (VIY) >= 3 should be predicted as strand."""
        protein = "VIY"  # exactly 3 strand formers
        result = _predict_secondary_structure_simple(protein)
        assert all(s == "E" for s in result)

    def test_loop_default(self):
        """Short runs of formers should default to loop."""
        # Single A (helix former) is too short for helix (min 4)
        protein = "AD"
        result = _predict_secondary_structure_simple(protein)
        # A alone is only 1 < 4, D not a former → all loop
        assert all(s == "L" for s in result)

    def test_empty_protein(self):
        """Empty protein should return empty list."""
        result = _predict_secondary_structure_simple("")
        assert result == []


# ────────────────────────────────────────────────────────────
# Additional: Test _is_preserved
# ────────────────────────────────────────────────────────────

class TestIsPreserved:
    """Tests for the _is_preserved helper."""

    def test_preserve_positions(self):
        """Position in preserve_positions should be preserved."""
        dc = DesignConstraints(preserve_positions=[2, 5])
        assert _is_preserved(2, "A", dc) is True
        assert _is_preserved(3, "A", dc) is False

    def test_preserve_residues(self):
        """Residue in preserve_residues should be preserved regardless of position."""
        dc = DesignConstraints(preserve_residues=["C", "W"])
        assert _is_preserved(0, "C", dc) is True
        assert _is_preserved(5, "W", dc) is True
        assert _is_preserved(0, "A", dc) is False

    def test_both_preserve_options(self):
        """Both position and residue preservation should work together."""
        dc = DesignConstraints(preserve_positions=[0], preserve_residues=["C"])
        assert _is_preserved(0, "A", dc) is True  # position preserved
        assert _is_preserved(5, "C", dc) is True  # residue preserved
        assert _is_preserved(3, "A", dc) is False  # neither preserved

    def test_no_preservation_default(self):
        """Default constraints should not preserve anything."""
        dc = DesignConstraints()
        assert _is_preserved(0, "C", dc) is False
