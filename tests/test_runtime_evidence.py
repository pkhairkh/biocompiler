"""
Tests for ``biocompiler.provenance.runtime_evidence``.

Each of the 34 narrowed axioms in ``proof/BioCompiler/SLOTVerification.lean``
gets a positive test (valid input → check passes) and a negative test
(invalid input → check fails).  Total: 34 × 2 = 68 tests, plus a handful
of bonus tests for ``run_all_evidence_checks``.

Reference: ``proof/BioCompiler/SLOTVerification.lean`` (lines 610-998).
"""

from __future__ import annotations

import math

import pytest

from biocompiler.provenance.runtime_evidence import (
    EvidenceCheckResult,
    run_all_evidence_checks,
    ALL_AXIOM_NAMES,
    TM_DOMAIN_WINDOW_SIZE,
    MRNA_STRUCTURE_WINDOW_SIZE,
    COTRANS_RAMP_CODONS,
    TM_DOMAIN_THRESHOLD,
    # TMHMM (3)
    check_tmhmm_window_size_contract,
    check_tmhmm_cytosolic_only_contract,
    check_tmhmm_threshold_sound_contract,
    # ViennaRNA (2)
    check_viennarna_window_size_contract,
    check_viennarna_deltaG_sound_contract,
    # AlphaFold (3)
    check_alphafold_ramp_window_contract,
    check_alphafold_cotrans_threshold_contract,
    check_alphafold_adaptation_index_sound_contract,
    # FoldX stable folding (3)
    check_foldx_stability_meaning_contract,
    check_foldx_estimated_deltaG_proxy_contract,
    check_foldx_stable_folding_sound_contract,
    # FoldX stability margin (2)
    check_foldx_ddg_threshold_meaningful_contract,
    check_foldx_stability_margin_sound_contract,
    # FoldX destabilizing mutation (2)
    check_foldx_max_ddg_meaningful_contract,
    check_foldx_destabilizing_mutation_sound_contract,
    # FoldX hydrophobic core (2)
    check_foldx_core_window_contract,
    check_foldx_core_quality_sound_contract,
    # ProteinSol (3)
    check_proteinsol_score_range_contract,
    check_proteinsol_gc_proxy_contract,
    check_proteinsol_solubility_sound_contract,
    # Aggrescan (3)
    check_aggrescan_window_size_contract,
    check_aggrescan_threshold_value_contract,
    check_aggrescan_no_aggregation_sound_contract,
    # ExPASy (3)
    check_expasy_pi_range_contract,
    check_expasy_gc_proxy_contract,
    check_expasy_charge_composition_sound_contract,
    # NetMHC (2)
    check_netmhc_score_nonneg_contract,
    check_netmhc_threshold_nonneg_contract,
    # NetMHCpan (2)
    check_netmhcpan_ic50_positive_contract,
    check_netmhcpan_threshold_positive_contract,
    # BepiPred (2)
    check_bepipred_score_nonneg_contract,
    check_bepipred_threshold_nonneg_contract,
    # IEDB (2)
    check_iedb_coverage_range_contract,
    check_iedb_threshold_range_contract,
)


# ─────────────────────────────────────────────────────────────────────────────
# Sanity: module-level invariants
# ─────────────────────────────────────────────────────────────────────────────

def test_all_axiom_names_count_is_34():
    """W1-A5 worklog specifies exactly 34 narrowed axioms."""
    assert len(ALL_AXIOM_NAMES) == 34
    # And each is unique
    assert len(set(ALL_AXIOM_NAMES)) == 34


def test_evidence_check_result_bool_convention():
    """EvidenceCheckResult.__bool__ should equal .passed."""
    r_pass = EvidenceCheckResult("axiom", True, "ok")
    r_fail = EvidenceCheckResult("axiom", False, "bad")
    assert bool(r_pass) is True
    assert bool(r_fail) is False


# ─────────────────────────────────────────────────────────────────────────────
# TMHMM (3 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestTMHMMChecks:
    """TMHMM tool-soundness evidence checks (3 axioms × 2 = 6 tests)."""

    def test_window_size_positive(self):
        r = check_tmhmm_window_size_contract(TM_DOMAIN_WINDOW_SIZE)
        assert r.axiom_name == "tmhmm_window_size_contract"
        assert r.passed is True

    def test_window_size_negative_wrong_value(self):
        r = check_tmhmm_window_size_contract(50)  # off-by-one
        assert r.passed is False

    def test_cytosolic_only_positive_true(self):
        r = check_tmhmm_cytosolic_only_contract(True)
        assert r.axiom_name == "tmhmm_cytosolic_only_contract"
        assert r.passed is True

    def test_cytosolic_only_negative_not_bool(self):
        r = check_tmhmm_cytosolic_only_contract(None)
        assert r.passed is False

    def test_threshold_sound_positive(self):
        # max hydrophobic fraction 0.4 < threshold 0.68 ⇒ passes
        r = check_tmhmm_threshold_sound_contract(0.4, 0.68, is_cytosolic=True)
        assert r.axiom_name == "tmhmm_threshold_sound_contract"
        assert r.passed is True

    def test_threshold_sound_negative(self):
        # max hydrophobic fraction 0.9 ≥ threshold 0.68 ⇒ fails
        r = check_tmhmm_threshold_sound_contract(0.9, 0.68, is_cytosolic=True)
        assert r.passed is False

    def test_threshold_sound_vacuous_when_not_cytosolic(self):
        """Bonus: when isCytosolic=False, the check is vacuously true."""
        r = check_tmhmm_threshold_sound_contract(0.99, 0.68, is_cytosolic=False)
        assert r.passed is True
        assert "Vacuously" in r.message


# ─────────────────────────────────────────────────────────────────────────────
# ViennaRNA (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestViennaRNAChecks:
    """ViennaRNA tool-soundness evidence checks (2 axioms × 2 = 4 tests)."""

    def test_window_size_positive(self):
        r = check_viennarna_window_size_contract(MRNA_STRUCTURE_WINDOW_SIZE)
        assert r.axiom_name == "viennarna_window_size_contract"
        assert r.passed is True

    def test_window_size_negative_wrong_value(self):
        r = check_viennarna_window_size_contract(31)
        assert r.passed is False

    def test_deltaG_sound_positive(self):
        # min estimated ΔG -5.0 > -10.0 threshold ⇒ passes
        r = check_viennarna_deltaG_sound_contract(-5.0, -10.0)
        assert r.axiom_name == "viennarna_deltaG_sound_contract"
        assert r.passed is True

    def test_deltaG_sound_negative(self):
        # min estimated ΔG -15.0 ≤ -10.0 threshold ⇒ fails
        r = check_viennarna_deltaG_sound_contract(-15.0, -10.0)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# AlphaFold co-translational (3 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestAlphaFoldChecks:
    """AlphaFold co-translational evidence checks (3 × 2 = 6 tests)."""

    def test_ramp_window_positive(self):
        r = check_alphafold_ramp_window_contract(COTRANS_RAMP_CODONS)
        assert r.axiom_name == "alphafold_ramp_window_contract"
        assert r.passed is True

    def test_ramp_window_negative_wrong_value(self):
        r = check_alphafold_ramp_window_contract(29)
        assert r.passed is False

    def test_cotrans_threshold_positive(self):
        r = check_alphafold_cotrans_threshold_contract(0.30)
        assert r.axiom_name == "alphafold_cotrans_threshold_contract"
        assert r.passed is True

    def test_cotrans_threshold_negative_out_of_range(self):
        r = check_alphafold_cotrans_threshold_contract(1.5)
        assert r.passed is False

    def test_adaptation_index_sound_positive(self):
        # rampAdaptationIndex 0.6 > threshold 0.30 ⇒ passes
        r = check_alphafold_adaptation_index_sound_contract(0.6, 0.30)
        assert r.axiom_name == "alphafold_adaptation_index_sound_contract"
        assert r.passed is True

    def test_adaptation_index_sound_negative(self):
        # rampAdaptationIndex 0.1 ≤ threshold 0.30 ⇒ fails
        r = check_alphafold_adaptation_index_sound_contract(0.1, 0.30)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# FoldX stable folding (3 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestFoldXStableFoldingChecks:
    """FoldX stable-folding evidence checks (3 × 2 = 6 tests)."""

    def test_stability_meaning_positive_default(self):
        """Default (no threshold arg) mirrors the Lean tautology 0 < 1."""
        r = check_foldx_stability_meaning_contract()
        assert r.axiom_name == "foldx_stability_meaning_contract"
        assert r.passed is True

    def test_stability_meaning_negative_non_finite_threshold(self):
        r = check_foldx_stability_meaning_contract(float("nan"))
        assert r.passed is False

    def test_estimated_deltaG_proxy_positive(self):
        r = check_foldx_estimated_deltaG_proxy_contract(-12.5)
        assert r.axiom_name == "foldx_estimated_deltaG_proxy_contract"
        assert r.passed is True

    def test_estimated_deltaG_proxy_negative_nan(self):
        r = check_foldx_estimated_deltaG_proxy_contract(float("nan"))
        assert r.passed is False

    def test_stable_folding_sound_positive(self):
        # estimated ΔG -8.0 < 0 ⇒ stable ⇒ passes
        r = check_foldx_stable_folding_sound_contract(-8.0)
        assert r.axiom_name == "foldx_stable_folding_sound_contract"
        assert r.passed is True

    def test_stable_folding_sound_negative(self):
        # estimated ΔG +3.0 ≥ 0 ⇒ not stable ⇒ fails
        r = check_foldx_stable_folding_sound_contract(3.0)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# FoldX stability margin (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestFoldXStabilityMarginChecks:
    """FoldX stability-margin evidence checks (2 × 2 = 4 tests)."""

    def test_ddg_threshold_meaningful_positive(self):
        r = check_foldx_ddg_threshold_meaningful_contract(5.0)
        assert r.axiom_name == "foldx_ddg_threshold_meaningful_contract"
        assert r.passed is True

    def test_ddg_threshold_meaningful_negative(self):
        r = check_foldx_ddg_threshold_meaningful_contract(-1.0)
        assert r.passed is False

    def test_stability_margin_sound_positive(self):
        # estimated ΔG -10.0 ≤ -ddgThreshold -5.0 ⇒ passes
        r = check_foldx_stability_margin_sound_contract(-10.0, 5.0)
        assert r.axiom_name == "foldx_stability_margin_sound_contract"
        assert r.passed is True

    def test_stability_margin_sound_negative(self):
        # estimated ΔG -3.0 > -ddgThreshold -5.0 ⇒ fails
        r = check_foldx_stability_margin_sound_contract(-3.0, 5.0)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# FoldX destabilizing mutation (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestFoldXDestabilizingMutationChecks:
    """FoldX destabilizing-mutation evidence checks (2 × 2 = 4 tests)."""

    def test_max_ddg_meaningful_positive(self):
        r = check_foldx_max_ddg_meaningful_contract(2.0)
        assert r.axiom_name == "foldx_max_ddg_meaningful_contract"
        assert r.passed is True

    def test_max_ddg_meaningful_negative(self):
        r = check_foldx_max_ddg_meaningful_contract(-0.5)
        assert r.passed is False

    def test_destabilizing_mutation_sound_positive(self):
        # estimated ΔG -7.0 ≤ -maxDDG -3.0 ⇒ passes
        r = check_foldx_destabilizing_mutation_sound_contract(-7.0, 3.0)
        assert r.axiom_name == "foldx_destabilizing_mutation_sound_contract"
        assert r.passed is True

    def test_destabilizing_mutation_sound_negative(self):
        # estimated ΔG -1.0 > -maxDDG -3.0 ⇒ fails
        r = check_foldx_destabilizing_mutation_sound_contract(-1.0, 3.0)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# FoldX hydrophobic core (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestFoldXHydrophobicCoreChecks:
    """FoldX hydrophobic-core evidence checks (2 × 2 = 4 tests)."""

    def test_core_window_positive(self):
        r = check_foldx_core_window_contract(TM_DOMAIN_WINDOW_SIZE)
        assert r.axiom_name == "foldx_core_window_contract"
        assert r.passed is True

    def test_core_window_negative_wrong_value(self):
        r = check_foldx_core_window_contract(100)
        assert r.passed is False

    def test_core_quality_sound_positive(self):
        # max hydrophobic fraction 0.75 ≥ threshold 0.68 ⇒ passes
        r = check_foldx_core_quality_sound_contract(0.75, 0.68)
        assert r.axiom_name == "foldx_core_quality_sound_contract"
        assert r.passed is True

    def test_core_quality_sound_negative(self):
        # max hydrophobic fraction 0.5 < threshold 0.68 ⇒ fails
        r = check_foldx_core_quality_sound_contract(0.5, 0.68)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# ProteinSol (3 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestProteinSolChecks:
    """ProteinSol evidence checks (3 × 2 = 6 tests)."""

    def test_score_range_positive(self):
        r = check_proteinsol_score_range_contract(0.45)
        assert r.axiom_name == "proteinsol_score_range_contract"
        assert r.passed is True

    def test_score_range_negative(self):
        r = check_proteinsol_score_range_contract(1.5)
        assert r.passed is False

    def test_gc_proxy_positive(self):
        r = check_proteinsol_gc_proxy_contract(0.5)
        assert r.axiom_name == "proteinsol_gc_proxy_contract"
        assert r.passed is True

    def test_gc_proxy_negative(self):
        r = check_proteinsol_gc_proxy_contract(-0.1)
        assert r.passed is False

    def test_solubility_sound_positive(self):
        # gcContent 0.55 ≥ minScore 0.40 ⇒ passes
        r = check_proteinsol_solubility_sound_contract(0.55, 0.40)
        assert r.axiom_name == "proteinsol_solubility_sound_contract"
        assert r.passed is True

    def test_solubility_sound_negative(self):
        # gcContent 0.30 < minScore 0.40 ⇒ fails
        r = check_proteinsol_solubility_sound_contract(0.30, 0.40)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# Aggrescan (3 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestAggrescanChecks:
    """Aggrescan evidence checks (3 × 2 = 6 tests)."""

    def test_window_size_positive(self):
        r = check_aggrescan_window_size_contract(TM_DOMAIN_WINDOW_SIZE)
        assert r.axiom_name == "aggrescan_window_size_contract"
        assert r.passed is True

    def test_window_size_negative_wrong_value(self):
        r = check_aggrescan_window_size_contract(20)
        assert r.passed is False

    def test_threshold_value_positive(self):
        r = check_aggrescan_threshold_value_contract(TM_DOMAIN_THRESHOLD)
        assert r.axiom_name == "aggrescan_threshold_value_contract"
        assert r.passed is True

    def test_threshold_value_negative_wrong_value(self):
        r = check_aggrescan_threshold_value_contract(0.50)
        assert r.passed is False

    def test_no_aggregation_sound_positive(self):
        # max hydrophobic fraction 0.4 < tmDomainThreshold 0.68 ⇒ passes
        r = check_aggrescan_no_aggregation_sound_contract(0.4)
        assert r.axiom_name == "aggrescan_no_aggregation_sound_contract"
        assert r.passed is True

    def test_no_aggregation_sound_negative(self):
        # max hydrophobic fraction 0.7 ≥ tmDomainThreshold 0.68 ⇒ fails
        r = check_aggrescan_no_aggregation_sound_contract(0.7)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# ExPASy (3 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestExPASyChecks:
    """ExPASy evidence checks (3 × 2 = 6 tests)."""

    def test_pi_range_positive(self):
        r = check_expasy_pi_range_contract(7.4)
        assert r.axiom_name == "expasy_pi_range_contract"
        assert r.passed is True

    def test_pi_range_negative(self):
        r = check_expasy_pi_range_contract(20.0)
        assert r.passed is False

    def test_gc_proxy_positive(self):
        r = check_expasy_gc_proxy_contract(0.6)
        assert r.axiom_name == "expasy_gc_proxy_contract"
        assert r.passed is True

    def test_gc_proxy_negative(self):
        r = check_expasy_gc_proxy_contract(1.2)
        assert r.passed is False

    def test_charge_composition_sound_positive(self):
        # gcContent 0.5 in [0.4, 0.6] ⇒ passes
        r = check_expasy_charge_composition_sound_contract(0.5, 0.4, 0.6)
        assert r.axiom_name == "expasy_charge_composition_sound_contract"
        assert r.passed is True

    def test_charge_composition_sound_negative(self):
        # gcContent 0.7 outside [0.4, 0.6] ⇒ fails
        r = check_expasy_charge_composition_sound_contract(0.7, 0.4, 0.6)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# NetMHC (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestNetMHCChecks:
    """NetMHC immunogenicity evidence checks (2 × 2 = 4 tests)."""

    def test_score_nonneg_positive(self):
        r = check_netmhc_score_nonneg_contract(0.5)
        assert r.axiom_name == "netmhc_score_nonneg_contract"
        assert r.passed is True

    def test_score_nonneg_negative(self):
        r = check_netmhc_score_nonneg_contract(-0.1)
        assert r.passed is False

    def test_threshold_nonneg_positive(self):
        r = check_netmhc_threshold_nonneg_contract(0.7)
        assert r.axiom_name == "netmhc_threshold_nonneg_contract"
        assert r.passed is True

    def test_threshold_nonneg_negative(self):
        r = check_netmhc_threshold_nonneg_contract(-0.5)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# NetMHCpan (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestNetMHCpanChecks:
    """NetMHCpan T-cell epitope evidence checks (2 × 2 = 4 tests)."""

    def test_ic50_positive_contract_positive(self):
        r = check_netmhcpan_ic50_positive_contract(500.0)
        assert r.axiom_name == "netmhcpan_ic50_positive_contract"
        assert r.passed is True

    def test_ic50_positive_contract_negative(self):
        r = check_netmhcpan_ic50_positive_contract(0.0)
        assert r.passed is False

    def test_threshold_positive_contract_positive(self):
        r = check_netmhcpan_threshold_positive_contract(500.0)
        assert r.axiom_name == "netmhcpan_threshold_positive_contract"
        assert r.passed is True

    def test_threshold_positive_contract_negative(self):
        r = check_netmhcpan_threshold_positive_contract(0.0)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# BepiPred (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestBepiPredChecks:
    """BepiPred B-cell epitope evidence checks (2 × 2 = 4 tests)."""

    def test_score_nonneg_positive(self):
        r = check_bepipred_score_nonneg_contract(0.35)
        assert r.axiom_name == "bepipred_score_nonneg_contract"
        assert r.passed is True

    def test_score_nonneg_negative(self):
        r = check_bepipred_score_nonneg_contract(-0.2)
        assert r.passed is False

    def test_threshold_nonneg_positive(self):
        r = check_bepipred_threshold_nonneg_contract(0.5)
        assert r.axiom_name == "bepipred_threshold_nonneg_contract"
        assert r.passed is True

    def test_threshold_nonneg_negative(self):
        r = check_bepipred_threshold_nonneg_contract(-0.1)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# IEDB (2 narrowed axioms)
# ─────────────────────────────────────────────────────────────────────────────

class TestIEDBChecks:
    """IEDB population-coverage evidence checks (2 × 2 = 4 tests)."""

    def test_coverage_range_positive(self):
        r = check_iedb_coverage_range_contract(0.85)
        assert r.axiom_name == "iedb_coverage_range_contract"
        assert r.passed is True

    def test_coverage_range_negative(self):
        r = check_iedb_coverage_range_contract(1.5)
        assert r.passed is False

    def test_threshold_range_positive(self):
        r = check_iedb_threshold_range_contract(0.9)
        assert r.axiom_name == "iedb_threshold_range_contract"
        assert r.passed is True

    def test_threshold_range_negative(self):
        r = check_iedb_threshold_range_contract(-0.1)
        assert r.passed is False


# ─────────────────────────────────────────────────────────────────────────────
# Aggregator: run_all_evidence_checks
# ─────────────────────────────────────────────────────────────────────────────

class TestRunAllEvidenceChecks:
    """Tests for the ``run_all_evidence_checks`` aggregator.

    Bonus tests (not counted in the 68) — verify that the aggregator
    dispatches to the right checks based on the dict keys present.
    """

    def test_empty_dict_returns_empty_list(self):
        assert run_all_evidence_checks({}) == []

    def test_non_dict_input_returns_failure(self):
        results = run_all_evidence_checks(None)  # type: ignore[arg-type]
        assert len(results) == 1
        assert results[0].passed is False

    def test_dispatches_tmhmm_checks(self):
        out = {
            "tmhmm_window_size": 51,
            "tmhmm_is_cytosolic": True,
            "tmhmm_max_hydrophobic_fraction": 0.4,
            "tmhmm_threshold": 0.68,
        }
        results = run_all_evidence_checks(out)
        names = {r.axiom_name for r in results}
        assert "tmhmm_window_size_contract" in names
        assert "tmhmm_cytosolic_only_contract" in names
        assert "tmhmm_threshold_sound_contract" in names
        assert all(r.passed for r in results)

    def test_dispatches_full_set_all_passing(self):
        """A complete tool-output dict with all valid values triggers
        all 34 checks (or close to it; some checks share keys)."""
        out = {
            # TMHMM
            "tmhmm_window_size": 51,
            "tmhmm_is_cytosolic": True,
            "tmhmm_max_hydrophobic_fraction": 0.4,
            "tmhmm_threshold": 0.68,
            # ViennaRNA
            "viennarna_window_size": 30,
            "viennarna_min_estimated_delta_g": -5.0,
            "viennarna_dg_threshold": -10.0,
            # AlphaFold
            "alphafold_ramp_codons": 30,
            "alphafold_cotrans_threshold": 0.30,
            "alphafold_ramp_adaptation_index": 0.6,
            # FoldX stable folding
            "foldx_stability_criterion_present": True,
            "foldx_estimated_delta_g": -8.0,
            "foldx_stable_folding_passed": True,
            # FoldX stability margin
            "foldx_ddg_threshold": 5.0,
            # FoldX destabilizing mutation
            "foldx_max_ddg": 3.0,
            # FoldX hydrophobic core
            "foldx_core_window_size": 51,
            "foldx_core_max_hydrophobic_fraction": 0.75,
            "foldx_core_threshold": 0.68,
            # ProteinSol
            "proteinsol_score": 0.45,
            "proteinsol_gc_content": 0.55,
            "proteinsol_min_score": 0.40,
            # Aggrescan
            "aggrescan_window_size": 51,
            "aggrescan_threshold": 0.68,
            "aggrescan_max_hydrophobic_fraction": 0.4,
            # ExPASy
            "expasy_pi": 7.4,
            "expasy_gc_content": 0.5,
            "expasy_pi_lo": 0.4,
            "expasy_pi_hi": 0.6,
            # NetMHC
            "netmhc_score": 0.5,
            "netmhc_max_score": 0.7,
            # NetMHCpan
            "netmhcpan_ic50": 500.0,
            "netmhcpan_ic50_threshold": 500.0,
            # BepiPred
            "bepipred_score": 0.35,
            "bepipred_score_threshold": 0.5,
            # IEDB
            "iedb_coverage": 0.85,
            "iedb_max_coverage": 0.9,
        }
        results = run_all_evidence_checks(out)
        # Every result must be a passing EvidenceCheckResult
        assert all(isinstance(r, EvidenceCheckResult) for r in results)
        assert all(r.passed for r in results), (
            "Failed: " + "; ".join(r.message for r in results if not r.passed)
        )
        # We should have run a substantial fraction of the 34 checks.
        # Some checks share keys (e.g., foldx_estimated_delta_g appears in
        # multiple checks), so the total run count can exceed 34.  We
        # require at least 34 — the size of ALL_AXIOM_NAMES.
        assert len(results) >= 34, (
            f"Expected ≥34 results, got {len(results)}: "
            f"{sorted(r.axiom_name for r in results)}"
        )

    def test_dispatches_with_failures(self):
        """If some checks fail, the aggregator should still return all results."""
        out = {
            "tmhmm_window_size": 50,  # wrong — fails
            "tmhmm_is_cytosolic": True,
            "tmhmm_max_hydrophobic_fraction": 0.9,  # too high — fails
            "tmhmm_threshold": 0.68,
        }
        results = run_all_evidence_checks(out)
        assert len(results) == 3
        # The first check should fail (wrong window size)
        ws = next(r for r in results if r.axiom_name == "tmhmm_window_size_contract")
        assert ws.passed is False
        # The third check should fail (fraction too high)
        ts = next(r for r in results if r.axiom_name == "tmhmm_threshold_sound_contract")
        assert ts.passed is False
