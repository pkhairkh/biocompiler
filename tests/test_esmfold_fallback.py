"""Test BioCompiler ESMFold Offline Fallback — Heuristic Structure Prediction."""

import pytest

from biocompiler.engines.esmfold_fallback import (
    predict_structure_heuristic,
    estimate_plddt_from_sequence,
    estimate_secondary_structure_from_sequence,
    compute_hydrophobicity_profile,
    compute_charge_profile,
    HEURISTIC_MAX_CONFIDENCE,
    KYTE_DOOLITTLE,
    CHOU_FASMAN_PROPENSITY,
    ChargeProfile,
    SecondaryStructureEstimate,
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
# TestEstimateSecondaryStructure
# ---------------------------------------------------------------------------

class TestEstimateSecondaryStructure:
    """Tests for estimate_secondary_structure_from_sequence."""

    def test_returns_ss_estimate(self):
        """Returns a SecondaryStructureEstimate dataclass."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        assert isinstance(ss, SecondaryStructureEstimate)

    def test_fractions_sum_to_one(self):
        """Helix + sheet + coil fractions sum to ~1.0."""
        for seq in [GLOBULAR_SEQ, DISORDERED_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ]:
            ss = estimate_secondary_structure_from_sequence(seq)
            total = ss.helix_fraction + ss.sheet_fraction + ss.coil_fraction
            assert total == pytest.approx(1.0, abs=0.01)

    def test_assignments_length_matches_protein(self):
        """Per-residue assignments list matches protein length."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        assert len(ss.assignments) == len(GLOBULAR_SEQ)

    def test_assignments_are_valid(self):
        """All assignments are 'H', 'E', or 'C'."""
        ss = estimate_secondary_structure_from_sequence(GLOBULAR_SEQ)
        for a in ss.assignments:
            assert a in ("H", "E", "C")

    def test_helix_forming_sequence(self):
        """Sequence rich in helix-formers should have some helix."""
        # Poly-Ala with some Leu — strong helix propensity
        helix_seq = "AALAAALAAALAAALAAALAAAL"
        ss = estimate_secondary_structure_from_sequence(helix_seq)
        assert ss.helix_fraction > 0.0

    def test_empty_protein(self):
        """Empty sequence returns all-coil."""
        ss = estimate_secondary_structure_from_sequence("")
        assert ss.helix_fraction == 0.0
        assert ss.sheet_fraction == 0.0
        assert ss.coil_fraction == 1.0


# ---------------------------------------------------------------------------
# TestEstimatePlddt
# ---------------------------------------------------------------------------

class TestEstimatePlddt:
    """Tests for estimate_plddt_from_sequence."""

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

    def test_method_is_heuristic_fallback(self):
        """Method is always 'heuristic_fallback'."""
        result = estimate_plddt_from_sequence(SHORT_SEQ)
        assert result["method"] == "heuristic_fallback"

    def test_plddt_below_max_confidence(self):
        """Estimated pLDDT never exceeds HEURISTIC_MAX_CONFIDENCE."""
        for seq in [GLOBULAR_SEQ, HYDROPHOBIC_SEQ, SHORT_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert result["estimated_mean_plddt"] <= HEURISTIC_MAX_CONFIDENCE

    def test_plddt_non_negative(self):
        """Estimated pLDDT is never negative."""
        for seq in [GLOBULAR_SEQ, DISORDERED_SEQ, SHORT_SEQ]:
            result = estimate_plddt_from_sequence(seq)
            assert result["estimated_mean_plddt"] >= 0.0

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

    def test_heuristic_details_has_scoring_breakdown(self):
        """Heuristic details include scoring components."""
        result = estimate_plddt_from_sequence(GLOBULAR_SEQ)
        details = result["heuristic_details"]
        assert "base" in details
        assert "hydrophobicity_bonus" in details
        assert "secondary_structure_bonus" in details
        assert "charge_balance_bonus" in details
        assert "disorder_penalty" in details
        assert "charge_clustering_penalty" in details

    def test_empty_protein_returns_zero(self):
        """Empty protein returns pLDDT of 0."""
        result = estimate_plddt_from_sequence("")
        assert result["estimated_mean_plddt"] == 0.0


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

    def test_method_is_heuristic_fallback(self):
        """Method is always 'heuristic_fallback'."""
        result = predict_structure_heuristic(SHORT_SEQ)
        assert result["method"] == "heuristic_fallback"

    def test_model_name_is_heuristic_v1(self):
        """Model name is 'heuristic_v1'."""
        result = predict_structure_heuristic(SHORT_SEQ)
        assert result["model_name"] == "heuristic_v1"

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

    def test_per_residue_plddt_below_max_confidence(self):
        """Each per-residue pLDDT is at most HEURISTIC_MAX_CONFIDENCE."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        for score in result["plddt_scores"]:
            assert score <= HEURISTIC_MAX_CONFIDENCE

    def test_per_residue_plddt_non_negative(self):
        """Each per-residue pLDDT is non-negative."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        for score in result["plddt_scores"]:
            assert score >= 0.0

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
        """Result includes secondary structure fractions."""
        result = predict_structure_heuristic(GLOBULAR_SEQ)
        ss = result["secondary_structure"]
        assert "helix_fraction" in ss
        assert "sheet_fraction" in ss
        assert "coil_fraction" in ss
        total = ss["helix_fraction"] + ss["sheet_fraction"] + ss["coil_fraction"]
        assert total == pytest.approx(1.0, abs=0.01)

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


# ---------------------------------------------------------------------------
# TestESMFoldIntegration
# ---------------------------------------------------------------------------

class TestESMFoldIntegration:
    """Tests for ESMFold engine integration with the heuristic fallback."""

    def test_predict_structure_uses_fallback_offline(self):
        """predict_structure uses heuristic fallback when ESMFold is unavailable."""
        from biocompiler.esmfold import predict_structure, is_esmfold_available

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
        from biocompiler.esmfold import predict_structure, is_esmfold_available

        if is_esmfold_available():
            pytest.skip("ESMFold is available; skipping offline integration test")

        result = predict_structure("MKWVTFISLLFLFSSAYS")
        assert result.pdb_string == ""

    def test_heuristic_result_has_plddt_scores(self):
        """Heuristic fallback produces per-residue pLDDT scores."""
        from biocompiler.esmfold import predict_structure, is_esmfold_available

        if is_esmfold_available():
            pytest.skip("ESMFold is available; skipping offline integration test")

        result = predict_structure("MKWVTFISLLFLFSSAYS")
        assert len(result.plddt_scores) > 0
        for score in result.plddt_scores:
            assert isinstance(score, float)
            assert score >= 0.0
            assert score <= HEURISTIC_MAX_CONFIDENCE


# ---------------------------------------------------------------------------
# TestEdgeCases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge case tests for the heuristic fallback."""

    def test_all_same_residue(self):
        """Poly-A sequence is handled without error."""
        result = predict_structure_heuristic("A" * 50)
        assert result["mean_plddt"] >= 0.0
        assert len(result["plddt_scores"]) == 50

    def test_alternating_charge(self):
        """Alternating positive/negative charge sequence."""
        result = predict_structure_heuristic("KEKDKEKDKEKD" * 3)
        assert result["mean_plddt"] >= 0.0
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
        """Proline-rich sequences get disorder penalty."""
        proline_seq = "PPPPPPPPPPPPPPPPPPPPPP"
        normal_seq = "MKWVTFISLLFLFSSAYSRGV"
        pro_result = estimate_plddt_from_sequence(proline_seq)
        normal_result = estimate_plddt_from_sequence(normal_seq)
        # Proline-rich should have higher disorder penalty
        assert pro_result["heuristic_details"]["disorder_penalty"] > normal_result["heuristic_details"]["disorder_penalty"]

    def test_heuristic_max_confidence_constant(self):
        """HEURISTIC_MAX_CONFIDENCE is below ESMFold's 'Low confidence' threshold."""
        # ESMFold pLDDT 50-70 is "Low confidence"; heuristic should be below that
        assert HEURISTIC_MAX_CONFIDENCE < 50.0
