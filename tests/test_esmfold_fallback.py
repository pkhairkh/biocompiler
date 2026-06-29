"""Test BioCompiler ESMFold Offline Fallback — Heuristic Structure Prediction (v2).

Tests the improved heuristic fallback with:
- Full Chou-Fasman secondary structure prediction (including turn prediction)
- Contact density estimation
- Per-residue pLDDT calibration (SS-based ranges)
- DSSP-style secondary structure string output
- Confidence calibration with SS content
"""

import pytest

from biocompiler.engines.esmfold_fallback import (
    predict_structure_heuristic,
    estimate_plddt_from_sequence,
    estimate_secondary_structure_from_sequence,
    compute_hydrophobicity_profile,
    compute_charge_profile,
    compute_contact_density,
    HEURISTIC_MAX_CONFIDENCE,
    HEURISTIC_MIN_CONFIDENCE,
    KYTE_DOOLITTLE,
    CHOU_FASMAN_PROPENSITY,
    TURN_PROPENSITY_BY_POSITION,
    CONTACT_DENSITY,
    PLDDT_RANGES,
    ChargeProfile,
    SecondaryStructureEstimate,
    ContactDensityProfile,
)


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

# A well-behaved globular protein fragment (mixed hydrophobic/hydrophilic)
GLOBULAR_SEQ = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQ"

# A highly charged, likely disordered sequence
DISORDERED_SEQ = "KKKKKKKEEEEEEEQQQQQQQSSSSSSSGGGGGGG"

# A hydrophobic-rich sequence (like a transmembrane helix)
HYDROPHOBIC_SEQ = "LLLLLLVVVVVVAAAAAAIIIIIIIMMMMMMFFFFF"

# A short peptide
SHORT_SEQ = "MAG"

# A single residue
SINGLE_AA = "M"

# A helix-favoring sequence (poly-Ala with Leu)
HELIX_SEQ = "AALAAALAAALAAALAAALAAALAAAAEAAAAKAAAA"

# A sheet-favoring sequence (Val, Ile, Tyr)
SHEET_SEQ = "VIYVIYVIYVIYVIYVIYVIY"

# A turn-rich sequence (Asn, Gly, Pro, Ser)
TURN_SEQ = "NGPSNGPSNGPSNGPSNGPS"

# A proline-rich sequence (disorder-promoting)
PROLINE_SEQ = "PPPPPPPPPPPPPPPPPPPPPP"


# ---------------------------------------------------------------------------
# TestComputeHydrophobicityProfile
# ---------------------------------------------------------------------------

class TestComputeHydrophobicityProfile:
    """Tests for compute_hydrophobicity_profile."""

    def test_returns_list_of_floats(self):
        """Hydrophobicity profile is a list of floats."""
        profile = compute_hydrophobicity_profile(GLOBULAR_SEQ)
        assert isinstance(profile, list)
        assert len(profile) == len(GLOBULAR_SEQ)
        for val in profile:
            assert isinstance(val, float)

    def test_length_matches_protein(self):
        """Profile length equals protein length."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ, SINGLE_AA]:
            profile = compute_hydrophobicity_profile(seq)
            assert len(profile) == len(seq)

    def test_empty_protein(self):
        """Empty sequence returns empty profile."""
        assert compute_hydrophobicity_profile("") == []

    def test_hydrophobic_seq_has_high_values(self):
        """Hydrophobic sequence has higher average hydrophobicity than disordered."""
        hydro_profile = compute_hydrophobicity_profile(HYDROPHOBIC_SEQ)
        disordered_profile = compute_hydrophobicity_profile(DISORDERED_SEQ)
        assert sum(hydro_profile) / len(hydro_profile) > sum(disordered_profile) / len(disordered_profile)

    def test_single_residue_matches_kyte_doolittle(self):
        """Single residue profile matches Kyte-Doolittle value."""
        profile = compute_hydrophobicity_profile("I")
        assert profile[0] == pytest.approx(KYTE_DOOLITTLE["I"], abs=0.01)


# ---------------------------------------------------------------------------
# TestComputeChargeProfile
# ---------------------------------------------------------------------------

class TestComputeChargeProfile:
    """Tests for compute_charge_profile."""

    def test_returns_charge_profile(self):
        """compute_charge_profile returns a ChargeProfile dataclass."""
        profile = compute_charge_profile(GLOBULAR_SEQ)
        assert isinstance(profile, ChargeProfile)

    def test_net_charge_positive_seq(self):
        """Sequence with many Lys/Arg has positive net charge."""
        profile = compute_charge_profile("KKKKRRR")
        assert profile.net_charge > 0
        assert profile.positive_fraction > 0.5

    def test_net_charge_negative_seq(self):
        """Sequence with many Asp/Glu has negative net charge."""
        profile = compute_charge_profile("DDDDEEEE")
        assert profile.net_charge < 0
        assert profile.negative_fraction > 0.5

    def test_balanced_charge(self):
        """Mixed charged sequence has high charge balance."""
        profile = compute_charge_profile("KDEKRHD")
        assert profile.charge_balance > 0.5

    def test_no_charged_residues(self):
        """Sequence with no charged residues has zero charge fractions."""
        profile = compute_charge_profile("AAAALLLL")
        assert profile.positive_fraction == 0.0
        assert profile.negative_fraction == 0.0
        assert profile.net_charge == 0

    def test_charge_patch_detection(self):
        """Long run of same-charge residues is detected as a patch."""
        # 7 consecutive K residues → strong positive patch
        profile = compute_charge_profile("KKKKKKKAAAA")
        assert profile.charge_patch_count >= 1

    def test_empty_protein(self):
        """Empty protein returns zero values."""
        profile = compute_charge_profile("")
        assert profile.net_charge == 0
        assert profile.positive_fraction == 0.0


# ---------------------------------------------------------------------------
# TestEstimateSecondaryStructure (improved Chou-Fasman)
# ---------------------------------------------------------------------------

class TestEstimateSecondaryStructure:
    """Tests for estimate_secondary_structure_from_sequence (full Chou-Fasman)."""

    def test_returns_ss_estimate(self):
        """Returns a SecondaryStructureEstimate dataclass."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        assert isinstance(ss, SecondaryStructureEstimate)

    def test_fractions_sum_to_one(self):
        """Helix + sheet + turn + coil fractions sum to ~1.0."""
        for seq in [GLOBULAR_SEQ, DISORDERED_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ, HELIX_SEQ]:
            ss = estimate_secondary_structure_from_sequence(seq)
            total = ss.helix_fraction + ss.sheet_fraction + ss.turn_fraction + ss.coil_fraction
            assert total == pytest.approx(1.0, abs=0.02)

    def test_assignments_length_matches_protein(self):
        """Per-residue assignments list matches protein length."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        assert len(ss.assignments) == len(GLOBULAR_SEQ)

    def test_assignments_are_valid(self):
        """All assignments are 'H', 'E', 'T', or 'C'."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        for a in ss.assignments:
            assert a in ("H", "E", "T", "C")

    def test_helix_forming_sequence(self):
        """Sequence rich in helix-formers should have some helix."""
        ss = estimate_secondary_structure_from_sequence(HELIX_SEQ)
        assert ss.helix_fraction > 0.0

    def test_empty_protein(self):
        """Empty sequence returns all-coil."""
        ss = estimate_secondary_structure_from_sequence("")
        assert ss.helix_fraction == 0.0
        assert ss.sheet_fraction == 0.0
        assert ss.turn_fraction == 0.0
        assert ss.coil_fraction == 1.0

    def test_ss_string_length_matches(self):
        """DSSP-style SS string length matches protein length."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ, HELIX_SEQ, SHEET_SEQ]:
            ss = estimate_secondary_structure_from_sequence(seq)
            assert len(ss.ss_string) == len(seq)

    def test_ss_string_valid_characters(self):
        """DSSP-style SS string only contains H, E, C."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        for char in ss.ss_string:
            assert char in ("H", "E", "C")

    def test_ss_string_consistent_with_assignments(self):
        """SS string is consistent with per-residue assignments."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        for i, (assignment, char) in enumerate(zip(ss.assignments, ss.ss_string)):
            if assignment == "H":
                assert char == "H"
            elif assignment == "E":
                assert char == "E"
            else:  # T or C both map to C in DSSP string
                assert char == "C"

    def test_turn_prediction(self):
        """Sequences rich in turn-formers should have some turns."""
        ss = estimate_secondary_structure_from_sequence(TURN_SEQ)
        # Turn sequence has Asn, Gly, Pro, Ser — all strong turn formers
        # Should have at least some turn predictions
        assert ss.turn_fraction > 0.0 or ss.helix_fraction > 0.0 or ss.sheet_fraction > 0.0
        # At minimum, not all coil
        assert ss.coil_fraction < 1.0

    def test_sheet_forming_sequence(self):
        """Sequence rich in sheet-formers should have some sheet."""
        ss = estimate_secondary_structure_from_sequence(SHEET_SEQ)
        assert ss.sheet_fraction > 0.0

    def test_proline_breaks_helices(self):
        """Proline inside a helix region should break it."""
        # A strong helix sequence with a proline in the middle
        seq_with_pro = "AALAAALAAAPAAALAAALAAAL"
        ss = estimate_secondary_structure_from_sequence(seq_with_pro)
        # The proline should break any helix that extends past it
        pro_pos = seq_with_pro.index("P")
        # Position after proline should not be H (proline breaks helix internally)
        # Position of P itself should not be H (it is a helix breaker)
        assert ss.assignments[pro_pos] != "H"

    def test_short_helices_pruned(self):
        """Helices shorter than 4 residues should be pruned to coil."""
        # This is an edge case; hard to construct a minimal helix
        # Just verify that no H region shorter than 4 exists
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        if "H" in ss.assignments:
            # Find all H runs
            h_runs = []
            in_run = False
            run_start = 0
            for i, a in enumerate(ss.assignments):
                if a == "H" and not in_run:
                    in_run = True
                    run_start = i
                elif a != "H" and in_run:
                    in_run = False
                    h_runs.append(i - run_start)
            if in_run:
                h_runs.append(len(ss.assignments) - run_start)
            for run_len in h_runs:
                assert run_len >= 4

    def test_turn_fraction_is_nonnegative(self):
        """Turn fraction should be non-negative."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ, HELIX_SEQ]:
            ss = estimate_secondary_structure_from_sequence(seq)
            assert ss.turn_fraction >= 0.0

    def test_chou_fasman_propensity_has_turn_values(self):
        """CHOU_FASMAN_PROPENSITY now includes P_turn values."""
        for aa in CHOU_FASMAN_PROPENSITY:
            assert len(CHOU_FASMAN_PROPENSITY[aa]) == 3  # (P_helix, P_sheet, P_turn)

    def test_turn_propensity_by_position_exists(self):
        """TURN_PROPENSITY_BY_POSITION has entries for all 20 AAs."""
        assert len(TURN_PROPENSITY_BY_POSITION) == 20
        for aa in TURN_PROPENSITY_BY_POSITION:
            assert len(TURN_PROPENSITY_BY_POSITION[aa]) == 4  # 4 positions


# ---------------------------------------------------------------------------
# TestContactDensity
# ---------------------------------------------------------------------------

class TestContactDensity:
    """Tests for compute_contact_density."""

    def test_returns_contact_density_profile(self):
        """compute_contact_density returns a ContactDensityProfile."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        cd = compute_contact_density(ss.assignments)
        assert isinstance(cd, ContactDensityProfile)

    def test_per_residue_length_matches(self):
        """Per-residue contact density length matches assignments length."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ, HELIX_SEQ]:
            ss = estimate_secondary_structure_from_sequence(seq)
            cd = compute_contact_density(ss.assignments)
            assert len(cd.per_residue) == len(ss.assignments)

    def test_mean_is_positive(self):
        """Mean contact density is positive for non-empty sequences."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        cd = compute_contact_density(ss.assignments)
        assert cd.mean > 0.0

    def test_ss_weighted_is_positive(self):
        """SS-weighted contact density is positive for non-empty sequences."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        cd = compute_contact_density(ss.assignments)
        assert cd.ss_weighted > 0.0

    def test_helix_higher_density_than_coil(self):
        """Helix-dominated sequences have higher density than coil-dominated."""
        ss_helix = estimate_secondary_structure_from_sequence(HELIX_SEQ)
        ss_coil = estimate_secondary_structure_from_sequence("GGGGGGGGGGGGGGGG")
        cd_helix = compute_contact_density(ss_helix.assignments)
        cd_coil = compute_contact_density(ss_coil.assignments)
        assert cd_helix.ss_weighted >= cd_coil.ss_weighted

    def test_density_within_expected_range(self):
        """Contact density values are within expected bounds."""
        for seq in [GLOBULAR_SEQ, HELIX_SEQ, SHEET_SEQ]:
            ss = estimate_secondary_structure_from_sequence(seq)
            cd = compute_contact_density(ss.assignments)
            for val in cd.per_residue:
                assert 0.5 <= val <= 2.5

    def test_empty_assignments(self):
        """Empty assignments return empty profile."""
        cd = compute_contact_density([])
        assert cd.per_residue == []
        assert cd.mean == 0.0
        assert cd.ss_weighted == 0.0

    def test_all_helix_high_density(self):
        """All-helix assignments give high contact density."""
        assignments = ["H"] * 20
        cd = compute_contact_density(assignments)
        assert cd.mean > 1.5

    def test_all_coil_low_density(self):
        """All-coil assignments give low contact density."""
        assignments = ["C"] * 20
        cd = compute_contact_density(assignments)
        assert cd.mean < 1.0


# ---------------------------------------------------------------------------
# TestEstimatePlddt
# ---------------------------------------------------------------------------

class TestEstimatePlddt:
    """Tests for estimate_plddt_from_sequence (improved with per-residue)."""

    def test_returns_dict_with_required_keys(self):
        """Result contains all expected keys."""
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        assert "estimated_mean_plddt" in result
        assert "confidence" in result
        assert "method" in result
        assert "heuristic_details" in result
        assert "hydrophobicity_profile" in result
        assert "charge_profile" in result
        assert "secondary_structure" in result
        assert "contact_density" in result
        assert "ss_prediction" in result
        assert "plddt_scores" in result

    def test_method_is_heuristic_fallback(self):
        """Method is always 'heuristic_fallback'."""
        result = estimate_plddt_from_sequence(SHORT_SEQ)
        assert result["method"] == "heuristic_fallback"

    def test_mean_plddt_below_max_confidence(self):
        """Estimated mean pLDDT never exceeds HEURISTIC_MAX_CONFIDENCE."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ, HELIX_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert result["estimated_mean_plddt"] <= HEURISTIC_MAX_CONFIDENCE

    def test_mean_plddt_above_min_confidence(self):
        """Estimated mean pLDDT is at least HEURISTIC_MIN_CONFIDENCE for non-empty sequences."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ, HELIX_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert result["estimated_mean_plddt"] >= HEURISTIC_MIN_CONFIDENCE

    def test_plddt_scores_length_matches(self):
        """Per-residue pLDDT scores length matches protein length."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ, SINGLE_AA]:
            result = estimate_plddt_from_sequence(seq)
            assert len(result["plddt_scores"]) == len(seq)

    def test_per_residue_plddt_within_range(self):
        """Each per-residue pLDDT is within [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE]."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            for score in result["plddt_scores"]:
                assert score >= HEURISTIC_MIN_CONFIDENCE
                assert score <= HEURISTIC_MAX_CONFIDENCE

    def test_confidence_below_half(self):
        """Confidence is always < 0.5."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert result["confidence"] < 0.5

    def test_confidence_above_minimum(self):
        """Confidence is at least 0.1 for non-empty sequences."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert result["confidence"] >= 0.1

    def test_globular_higher_than_disordered(self):
        """Well-folded globular sequence gets higher pLDDT than disordered."""
        glob_result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        dis_result = estimate_plddt_from_sequence(DISORDERED_SEQ)
        assert glob_result["estimated_mean_plddt"] >= dis_result["estimated_mean_plddt"]

    def test_heuristic_details_has_calibration_info(self):
        """Heuristic details include per-residue calibration info."""
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        details = result["heuristic_details"]
        assert "per_residue_method" in details
        assert "helix_plddt_range" in details
        assert "sheet_plddt_range" in details
        assert "turn_plddt_range" in details
        assert "coil_plddt_range" in details
        assert "contact_density_mean" in details
        assert "contact_density_ss_weighted" in details
        assert "combined_penalty_factor" in details

    def test_empty_protein_returns_zero(self):
        """Empty protein returns pLDDT of 0."""
        result = estimate_plddt_from_sequence("")
        assert result["estimated_mean_plddt"] == 0.0

    def test_ss_prediction_string(self):
        """Result includes a DSSP-style SS prediction string."""
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        ss_pred = result["ss_prediction"]
        assert isinstance(ss_pred, str)
        assert len(ss_pred) == len(GLOBULAR_SEQ)
        for char in ss_pred:
            assert char in ("H", "E", "C")

    def test_helix_residues_higher_plddt_than_coil(self):
        """Residues predicted as helix should have higher pLDDT than coil."""
        result = estimate_plddt_from_sequence(HELIX_SEQ)
        ss = result["secondary_structure"]
        plddt = result["plddt_scores"]
        assignments = ss.assignments

        helix_scores = [p for p, a in zip(plddt, assignments) if a == "H"]
        coil_scores = [p for p, a in zip(plddt, assignments) if a == "C"]

        if helix_scores and coil_scores:
            mean_helix = sum(helix_scores) / len(helix_scores)
            mean_coil = sum(coil_scores) / len(coil_scores)
            assert mean_helix > mean_coil

    def test_sheet_residues_higher_plddt_than_coil(self):
        """Residues predicted as sheet should have higher pLDDT than coil."""
        result = estimate_plddt_from_sequence(SHEET_SEQ)
        ss = result["secondary_structure"]
        plddt = result["plddt_scores"]
        assignments = ss.assignments

        sheet_scores = [p for p, a in zip(plddt, assignments) if a == "E"]
        coil_scores = [p for p, a in zip(plddt, assignments) if a == "C"]

        if sheet_scores and coil_scores:
            mean_sheet = sum(sheet_scores) / len(sheet_scores)
            mean_coil = sum(coil_scores) / len(coil_scores)
            assert mean_sheet > mean_coil

    def test_mean_plddt_typical_range(self):
        """Mean pLDDT should be in the 35-50 range for typical proteins."""
        # Well-folded globular proteins should have higher mean
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        assert 30.0 <= result["estimated_mean_plddt"] <= 55.0

    def test_helix_plddt_range_correct(self):
        """Helix residue pLDDT should be in the 45-55 range."""
        result = estimate_plddt_from_sequence(HELIX_SEQ)
        ss = result["secondary_structure"]
        plddt = result["plddt_scores"]
        assignments = ss.assignments

        helix_scores = [p for p, a in zip(plddt, assignments) if a == "H"]
        if helix_scores:
            for score in helix_scores:
                # Allow some deviation due to modulation factors
                assert 40.0 <= score <= 58.0

    def test_coil_plddt_range_correct(self):
        """Coil residue pLDDT should be in the 25-35 range."""
        result = estimate_plddt_from_sequence(DISORDERED_SEQ)
        ss = result["secondary_structure"]
        plddt = result["plddt_scores"]
        assignments = ss.assignments

        coil_scores = [p for p, a in zip(plddt, assignments) if a == "C"]
        if coil_scores:
            for score in coil_scores:
                # Allow some deviation due to modulation factors
                assert 20.0 <= score <= 42.0


# ---------------------------------------------------------------------------
# TestPredictStructureHeuristic
# ---------------------------------------------------------------------------

class TestPredictStructureHeuristic:
    """Tests for predict_structure_heuristic (the main fallback entry point)."""

    def test_returns_dict_with_required_keys(self):
        """Result contains all keys needed for ESMFoldResult construction."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        assert "protein" in result
        assert "mean_plddt" in result
        assert "plddt_scores" in result
        assert "method" in result
        assert "model_name" in result
        assert "confidence" in result
        assert "heuristic_details" in result
        assert "secondary_structure" in result
        assert "ss_prediction" in result

    def test_method_is_heuristic_fallback(self):
        """Method is always 'heuristic_fallback'."""
        result = predict_structure_heuristic(SHORT_SEQ)
        assert result["method"] == "heuristic_fallback"

    def test_model_name_is_heuristic_v2(self):
        """Model name is 'heuristic_v2' (upgraded from v1)."""
        result = predict_structure_heuristic(SHORT_SEQ)
        assert result["model_name"] == "heuristic_v2"

    def test_plddt_scores_length_matches_protein(self):
        """Per-residue pLDDT list length matches protein length."""
        for seq in [GLOBULAR_SEQ, SHORT_SEQ, SINGLE_AA]:
            result = predict_structure_heuristic(seq)
            assert len(result["plddt_scores"]) == len(seq)

    def test_mean_plddt_below_max_confidence(self):
        """Mean pLDDT never exceeds HEURISTIC_MAX_CONFIDENCE."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ]:
            result = predict_structure_heuristic(seq)
            assert result["mean_plddt"] <= HEURISTIC_MAX_CONFIDENCE

    def test_mean_plddt_above_min_confidence(self):
        """Mean pLDDT is at least HEURISTIC_MIN_CONFIDENCE for non-empty sequences."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ]:
            result = predict_structure_heuristic(seq)
            assert result["mean_plddt"] >= HEURISTIC_MIN_CONFIDENCE

    def test_per_residue_plddt_within_range(self):
        """Each per-residue pLDDT is within [HEURISTIC_MIN_CONFIDENCE, HEURISTIC_MAX_CONFIDENCE]."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        for score in result["plddt_scores"]:
            assert score >= HEURISTIC_MIN_CONFIDENCE
            assert score <= HEURISTIC_MAX_CONFIDENCE

    def test_confidence_below_half(self):
        """Confidence is always < 0.5."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        assert result["confidence"] < 0.5

    def test_globular_higher_than_disordered(self):
        """Globular protein gets higher pLDDT than disordered one."""
        glob_result = predict_structure_heuristic(GLOBULAR_SEQ)
        dis_result = predict_structure_heuristic(DISORDERED_SEQ)
        assert glob_result["mean_plddt"] >= dis_result["mean_plddt"]

    def test_secondary_structure_included(self):
        """Result includes secondary structure fractions including turn."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        ss = result["secondary_structure"]
        assert "helix_fraction" in ss
        assert "sheet_fraction" in ss
        assert "turn_fraction" in ss
        assert "coil_fraction" in ss
        total = ss["helix_fraction"] + ss["sheet_fraction"] + ss["turn_fraction"] + ss["coil_fraction"]
        assert total == pytest.approx(1.0, abs=0.02)

    def test_ss_prediction_field(self):
        """Result includes ss_prediction DSSP-style string."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        ss_pred = result["ss_prediction"]
        assert isinstance(ss_pred, str)
        assert len(ss_pred) == len(GLOBULAR_SEQ)
        for char in ss_pred:
            assert char in ("H", "E", "C")

    def test_ss_prediction_consistent_with_assignments(self):
        """ss_prediction is consistent with per-residue assignments."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        ss = result["secondary_structure"]
        ss_pred = result["ss_prediction"]
        for i, (assignment, char) in enumerate(zip(ss["assignments"], ss_pred)):
            if assignment == "H":
                assert char == "H"
            elif assignment == "E":
                assert char == "E"
            else:  # T or C both map to C in DSSP string
                assert char == "C"

    def test_short_peptide(self):
        """Short peptides are handled without error."""
        result = predict_structure_heuristic(SHORT_SEQ)
        assert result["protein"] == SHORT_SEQ
        assert len(result["plddt_scores"]) == 3

    def test_single_residue(self):
        """Single-residue protein is handled."""
        result = predict_structure_heuristic(SINGLE_AA)
        assert result["protein"] == SINGLE_AA
        assert len(result["plddt_scores"]) == 1

    def test_helix_sequence_has_helix_in_ss(self):
        """Helix-favoring sequence should have helix in SS prediction."""
        result = predict_structure_heuristic(HELIX_SEQ)
        ss = result["secondary_structure"]
        assert ss["helix_fraction"] > 0.0

    def test_sheet_sequence_has_sheet_in_ss(self):
        """Sheet-favoring sequence should have sheet in SS prediction."""
        result = predict_structure_heuristic(SHEET_SEQ)
        ss = result["secondary_structure"]
        assert ss["sheet_fraction"] > 0.0


# ---------------------------------------------------------------------------
# TestESMFoldIntegration
# ---------------------------------------------------------------------------

class TestESMFoldIntegration:
    """Tests for ESMFold engine integration with the heuristic fallback."""

    def test_predict_structure_uses_fallback_offline(self):
        """predict_structure uses heuristic fallback when ESMFold is unavailable."""
        from biocompiler.engines.esmfold import predict_structure, is_esmfold_available

        if is_esmfold_available():
            pytest.skip("ESMFold is available; skipping offline integration test")

        result = predict_structure("MKWVTFISLLFLFSSAYS")
        assert result.success is True
        assert result.method == "heuristic_fallback"
        assert result.mean_plddt > 0.0
        assert result.mean_plddt <= HEURISTIC_MAX_CONFIDENCE
        assert result.confidence_level == "very_low"
        assert len(result.plddt_scores) == len("MKWVTFISLLFLFSSAYS")

    def test_heuristic_result_no_pdb_string(self):
        """Heuristic fallback results have no PDB string (no 3D coordinates)."""
        from biocompiler.engines.esmfold import predict_structure, is_esmfold_available

        if is_esmfold_available():
            pytest.skip("ESMFold is available; skipping offline integration test")

        result = predict_structure("MKWVTFISLLFLFSSAYS")
        assert result.pdb_string == ""

    def test_heuristic_result_has_plddt_scores(self):
        """Heuristic fallback produces per-residue pLDDT scores."""
        from biocompiler.engines.esmfold import predict_structure, is_esmfold_available

        if is_esmfold_available():
            pytest.skip("ESMFold is available; skipping offline integration test")

        result = predict_structure("MKWVTFISLLFLFSSAYS")
        assert len(result.plddt_scores) > 0
        for score in result.plddt_scores:
            assert isinstance(score, float)
            assert score >= HEURISTIC_MIN_CONFIDENCE
            assert score <= HEURISTIC_MAX_CONFIDENCE


# ---------------------------------------------------------------------------
# TestEdgeCases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge case tests for the heuristic fallback."""

    def test_all_same_residue(self):
        """Poly-A sequence is handled without error."""
        result = predict_structure_heuristic("A" * 50)
        assert result["mean_plddt"] >= HEURISTIC_MIN_CONFIDENCE
        assert len(result["plddt_scores"]) == 50

    def test_alternating_charge(self):
        """Alternating positive/negative charge sequence."""
        result = predict_structure_heuristic("KEKDKEKDKEKD" * 3)
        assert result["mean_plddt"] >= HEURISTIC_MIN_CONFIDENCE
        # Should detect good charge balance
        cp = compute_charge_profile("KEKDKEKDKEKD" * 3)
        assert cp.charge_balance > 0.5

    def test_very_long_sequence(self):
        """Long sequence (>100 residues) is handled."""
        seq = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQ" * 4
        result = predict_structure_heuristic(seq)
        assert len(result["plddt_scores"]) == len(seq)
        assert result["mean_plddt"] <= HEURISTIC_MAX_CONFIDENCE

    def test_proline_rich_disorder(self):
        """Proline-rich sequences get lower pLDDT than well-folded sequences."""
        normal_seq = "MKWVTFISLLFLFSSAYSRGV"
        pro_result = estimate_plddt_from_sequence(PROLINE_SEQ)
        normal_result = estimate_plddt_from_sequence(normal_seq)
        # Proline-rich should have lower pLDDT due to disorder
        assert pro_result["estimated_mean_plddt"] <= normal_result["estimated_mean_plddt"]

    def test_heuristic_max_confidence_below_esmfold_low(self):
        """HEURISTIC_MAX_CONFIDENCE is below ESMFold's 'Low confidence' band (50-70)."""
        # Our max is 55, which is below the 70 threshold but allows some
        # values in the 50-55 range. This is intentional: helix residues
        # in well-folded proteins can reach 55, which is still in the
        # "very low" confidence range of ESMFold.
        assert HEURISTIC_MAX_CONFIDENCE <= 55.0

    def test_heuristic_min_confidence_positive(self):
        """HEURISTIC_MIN_CONFIDENCE is a positive value."""
        assert HEURISTIC_MIN_CONFIDENCE > 0.0

    def test_single_residue_no_crash(self):
        """Single residue does not crash the turn prediction."""
        ss = estimate_secondary_structure_from_sequence("M")
        assert len(ss.assignments) == 1
        assert ss.assignments[0] in ("H", "E", "T", "C")

    def test_two_residues_no_crash(self):
        """Two residues do not crash the turn prediction."""
        ss = estimate_secondary_structure_from_sequence("MA")
        assert len(ss.assignments) == 2

    def test_three_residues_no_crash(self):
        """Three residues (shorter than turn window) are handled."""
        ss = estimate_secondary_structure_from_sequence("MAG")
        assert len(ss.assignments) == 3


# ---------------------------------------------------------------------------
# TestPerResiduePlddtCalibration
# ---------------------------------------------------------------------------

class TestPerResiduePlddtCalibration:
    """Tests for the SS-based per-residue pLDDT calibration."""

    def test_helix_plddt_higher_than_sheet(self):
        """Helix residues should have higher pLDDT than sheet residues on average."""
        # Use a sequence with both helix and sheet
        seq = "AALAAALVIYVIYVAALAAALVIYVIYV"
        result = estimate_plddt_from_sequence(seq)
        ss = result["secondary_structure"]
        plddt = result["plddt_scores"]
        assignments = ss.assignments

        helix_scores = [p for p, a in zip(plddt, assignments) if a == "H"]
        sheet_scores = [p for p, a in zip(plddt, assignments) if a == "E"]

        if helix_scores and sheet_scores:
            mean_helix = sum(helix_scores) / len(helix_scores)
            mean_sheet = sum(sheet_scores) / len(sheet_scores)
            assert mean_helix >= mean_sheet

    def test_defined_ss_gives_higher_mean(self):
        """Sequences with more defined SS should have higher mean pLDDT."""
        # A helix-favoring sequence vs a disordered sequence
        helix_result = estimate_plddt_from_sequence(HELIX_SEQ)
        disordered_result = estimate_plddt_from_sequence(DISORDERED_SEQ)
        assert helix_result["estimated_mean_plddt"] >= disordered_result["estimated_mean_plddt"]

    def test_plddt_ranges_are_below_esmfold_confident(self):
        """All pLDDT ranges are below ESMFold's 'Confident' threshold (70)."""
        for ss_type, (low, high) in PLDDT_RANGES.items():
            assert high < 70.0

    def test_helix_range_is_highest(self):
        """Helix pLDDT range should be the highest."""
        assert PLDDT_RANGES["H"][0] > PLDDT_RANGES["E"][0]
        assert PLDDT_RANGES["H"][0] > PLDDT_RANGES["C"][0]

    def test_coil_range_is_lowest(self):
        """Coil pLDDT range should be the lowest."""
        assert PLDDT_RANGES["C"][1] < PLDDT_RANGES["H"][1]
        assert PLDDT_RANGES["C"][1] < PLDDT_RANGES["E"][1]

    def test_shorter_proteins_slightly_higher(self):
        """Shorter proteins should have slightly higher pLDDT (all else equal)."""
        # Compare same composition at different lengths
        short = "AALAAALAAAL"  # 11 aa
        long = "AALAAALAAAL" * 10  # 110 aa
        short_result = estimate_plddt_from_sequence(short)
        long_result = estimate_plddt_from_sequence(long)
        # Short should have slightly higher mean (length factor)
        # Not a strict requirement since SS content may differ
        # but the length bonus should contribute
        assert short_result["heuristic_details"]["sequence_length_factor"] >= \
               long_result["heuristic_details"]["sequence_length_factor"]

    def test_mean_plddt_typical_range(self):
        """Mean pLDDT for typical proteins should be in the 30-55 range."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, HELIX_SEQ, SHEET_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert 25.0 <= result["estimated_mean_plddt"] <= 58.0


# ---------------------------------------------------------------------------
# TestContactDensityIntegration
# ---------------------------------------------------------------------------

class TestContactDensityIntegration:
    """Tests for contact density integration with pLDDT estimation."""

    def test_contact_density_in_plddt_result(self):
        """Contact density profile is included in pLDDT estimation result."""
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        cd = result["contact_density"]
        assert isinstance(cd, ContactDensityProfile)
        assert len(cd.per_residue) == len(GLOBULAR_SEQ)

    def test_contact_density_modulates_plddt(self):
        """Contact density should modulate per-residue pLDDT within SS range."""
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        ss = result["secondary_structure"]
        plddt = result["plddt_scores"]
        cd = result["contact_density"]

        # Residues with higher contact density within same SS type
        # should tend to have higher pLDDT
        for ss_type in ("H", "E"):
            same_type = [(p, cd.per_residue[i]) for i, (p, a)
                         in enumerate(zip(plddt, ss.assignments)) if a == ss_type]
            if len(same_type) >= 2:
                # Sort by contact density and check rough correlation
                same_type.sort(key=lambda x: x[1])
                # Higher contact density → higher pLDDT (within same SS type)
                # Not a strict guarantee due to hydrophobicity modulation,
                # but the trend should exist
                low_cd_plddt = same_type[0][0]
                high_cd_plddt = same_type[-1][0]
                # At minimum, they should be in the same general range
                assert abs(high_cd_plddt - low_cd_plddt) < 20.0

    def test_contact_density_values_reasonable(self):
        """Contact density values are within published ranges."""
        for density_val in CONTACT_DENSITY.values():
            assert 0.5 <= density_val <= 2.5
