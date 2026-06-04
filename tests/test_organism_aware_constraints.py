"""Tests for organism-aware constraint selection in build_csp_model (F1.2).

Validates that the CSP solver pipeline skips eukaryote-specific constraints
(splice sites, CpG islands) when the target organism is prokaryotic, while
keeping organism-agnostic constraints (ATTTA motifs, T-runs) for all organisms.

This is the key fix that recovers ~0.27 CAI on prokaryotic targets by
removing constraints that unnecessarily restrict the codon search space.
"""

from __future__ import annotations

import pytest

from biocompiler.solver.constraints import (
    NoCrypticSpliceConstraint,
    NoCpGIslandConstraint,
    NoATTTAMotifConstraint,
    NoTRunConstraint,
    MinimizeCpG,
    MaximizeCAI,
    MinimizeMRNADG,
    TranslationConstraint,
    GCRangeConstraint,
    NoRestrictionSiteConstraint,
    build_csp_model,
)
from biocompiler.solver.types import SolverConfig


# ─── Helpers ─────────────────────────────────────────────────────────


def _hard_constraint_names(model) -> set[str]:
    """Return the set of hard constraint names in the model."""
    return {c.name for c in model.hard_constraints}


def _soft_constraint_names(model) -> set[str]:
    """Return the set of soft constraint names in the model."""
    return {c.name for c in model.soft_constraints}


# Use a short protein for fast test execution
_PROTEIN = "MVSKGE"

# Prokaryote organism name as it appears in SUPPORTED_ORGANISMS
_PROKARYOTE = "Escherichia_coli"


# ─── Prokaryotic organisms: eukaryote-only constraints are skipped ──


class TestProkaryoticConstraintSkips:
    """Prokaryotic organisms should NOT include eukaryote-only constraints."""

    def test_no_cryptic_splice_for_prokaryote(self) -> None:
        """NoCrypticSpliceConstraint should be skipped for prokaryotes."""
        model = build_csp_model(_PROTEIN, _PROKARYOTE)
        assert "NoCrypticSpliceConstraint" not in _hard_constraint_names(model), (
            f"NoCrypticSpliceConstraint should NOT be present for {_PROKARYOTE}"
        )

    def test_no_cpg_island_for_prokaryote(self) -> None:
        """NoCpGIslandConstraint should be skipped for prokaryotes."""
        model = build_csp_model(_PROTEIN, _PROKARYOTE)
        assert "NoCpGIslandConstraint" not in _hard_constraint_names(model), (
            f"NoCpGIslandConstraint should NOT be present for {_PROKARYOTE}"
        )

    def test_no_minimize_cpg_for_prokaryote(self) -> None:
        """MinimizeCpG soft constraint should be skipped for prokaryotes."""
        model = build_csp_model(_PROTEIN, _PROKARYOTE)
        assert "MinimizeCpG" not in _soft_constraint_names(model), (
            f"MinimizeCpG should NOT be present for {_PROKARYOTE}"
        )


# ─── Eukaryotic organisms: eukaryote-only constraints are included ──


class TestEukaryoticConstraintInclusion:
    """Eukaryotic organisms SHOULD include eukaryote-only constraints."""

    @pytest.mark.parametrize(
        "organism",
        ["Homo_sapiens", "Saccharomyces_cerevisiae", "Mus_musculus", "CHO_K1"],
    )
    def test_cryptic_splice_for_eukaryote(self, organism: str) -> None:
        """NoCrypticSpliceConstraint should be present for eukaryotes."""
        model = build_csp_model(_PROTEIN, organism)
        assert "NoCrypticSpliceConstraint" in _hard_constraint_names(model), (
            f"NoCrypticSpliceConstraint should be present for {organism}"
        )

    @pytest.mark.parametrize(
        "organism",
        ["Homo_sapiens", "Saccharomyces_cerevisiae", "Mus_musculus", "CHO_K1"],
    )
    def test_cpg_island_for_eukaryote(self, organism: str) -> None:
        """NoCpGIslandConstraint should be present for eukaryotes."""
        model = build_csp_model(_PROTEIN, organism)
        assert "NoCpGIslandConstraint" in _hard_constraint_names(model), (
            f"NoCpGIslandConstraint should be present for {organism}"
        )

    @pytest.mark.parametrize(
        "organism",
        ["Homo_sapiens", "Saccharomyces_cerevisiae", "Mus_musculus", "CHO_K1"],
    )
    def test_minimize_cpg_for_eukaryote(self, organism: str) -> None:
        """MinimizeCpG soft constraint should be present for eukaryotes."""
        model = build_csp_model(_PROTEIN, organism)
        assert "MinimizeCpG" in _soft_constraint_names(model), (
            f"MinimizeCpG should be present for {organism}"
        )


# ─── Universal constraints: present for ALL organisms ────────────────


class TestUniversalConstraints:
    """Constraints relevant to both domains should be present for ALL organisms."""

    @pytest.mark.parametrize(
        "organism",
        ["Escherichia_coli", "Homo_sapiens", "Saccharomyces_cerevisiae"],
    )
    def test_attta_motif_for_all_organisms(self, organism: str) -> None:
        """NoATTTAMotifConstraint should be present for all organisms."""
        model = build_csp_model(_PROTEIN, organism)
        assert "NoATTTAMotifConstraint" in _hard_constraint_names(model), (
            f"NoATTTAMotifConstraint should be present for {organism}"
        )

    @pytest.mark.parametrize(
        "organism",
        ["Escherichia_coli", "Homo_sapiens", "Saccharomyces_cerevisiae"],
    )
    def test_t_run_for_all_organisms(self, organism: str) -> None:
        """NoTRunConstraint should be present for all organisms."""
        model = build_csp_model(_PROTEIN, organism)
        assert "NoTRunConstraint" in _hard_constraint_names(model), (
            f"NoTRunConstraint should be present for {organism}"
        )

    @pytest.mark.parametrize(
        "organism",
        ["Escherichia_coli", "Homo_sapiens"],
    )
    def test_translation_for_all_organisms(self, organism: str) -> None:
        """TranslationConstraint should be present for all organisms."""
        model = build_csp_model(_PROTEIN, organism)
        assert "TranslationConstraint" in _hard_constraint_names(model), (
            f"TranslationConstraint should be present for {organism}"
        )

    @pytest.mark.parametrize(
        "organism",
        ["Escherichia_coli", "Homo_sapiens"],
    )
    def test_gc_range_for_all_organisms(self, organism: str) -> None:
        """GCRangeConstraint should be present for all organisms."""
        model = build_csp_model(_PROTEIN, organism)
        assert "GCRangeConstraint" in _hard_constraint_names(model), (
            f"GCRangeConstraint should be present for {organism}"
        )

    @pytest.mark.parametrize(
        "organism",
        ["Escherichia_coli", "Homo_sapiens"],
    )
    def test_maximize_cai_for_all_organisms(self, organism: str) -> None:
        """MaximizeCAI should be present for all organisms."""
        model = build_csp_model(_PROTEIN, organism)
        assert "MaximizeCAI" in _soft_constraint_names(model), (
            f"MaximizeCAI should be present for {organism}"
        )


# ─── auto_detect_organism_domain flag ───────────────────────────────


class TestAutoDetectOrganismDomain:
    """The auto_detect_organism_domain flag controls constraint filtering."""

    def test_flag_true_skips_eukaryote_constraints_for_prokaryote(self) -> None:
        """When True (default), eukaryote constraints are skipped for prokaryotes."""
        config = SolverConfig(auto_detect_organism_domain=True)
        model = build_csp_model(_PROTEIN, _PROKARYOTE, config)
        assert "NoCrypticSpliceConstraint" not in _hard_constraint_names(model)
        assert "NoCpGIslandConstraint" not in _hard_constraint_names(model)
        assert "MinimizeCpG" not in _soft_constraint_names(model)

    def test_flag_false_includes_eukaryote_constraints_for_prokaryote(self) -> None:
        """When False, all constraints are applied regardless of organism (backward compat)."""
        config = SolverConfig(auto_detect_organism_domain=False)
        model = build_csp_model(_PROTEIN, _PROKARYOTE, config)
        # With auto_detect OFF, even prokaryotes get eukaryotic constraints
        assert "NoCrypticSpliceConstraint" in _hard_constraint_names(model), (
            "With auto_detect_organism_domain=False, NoCrypticSpliceConstraint "
            "should be present even for prokaryotes (backward compat)"
        )
        assert "NoCpGIslandConstraint" in _hard_constraint_names(model), (
            "With auto_detect_organism_domain=False, NoCpGIslandConstraint "
            "should be present even for prokaryotes (backward compat)"
        )
        assert "MinimizeCpG" in _soft_constraint_names(model), (
            "With auto_detect_organism_domain=False, MinimizeCpG "
            "should be present even for prokaryotes (backward compat)"
        )

    def test_flag_true_includes_eukaryote_constraints_for_eukaryote(self) -> None:
        """When True, eukaryote constraints ARE included for eukaryotes."""
        config = SolverConfig(auto_detect_organism_domain=True)
        model = build_csp_model(_PROTEIN, "Homo_sapiens", config)
        assert "NoCrypticSpliceConstraint" in _hard_constraint_names(model)
        assert "NoCpGIslandConstraint" in _hard_constraint_names(model)
        assert "MinimizeCpG" in _soft_constraint_names(model)

    def test_default_flag_value_is_true(self) -> None:
        """SolverConfig.auto_detect_organism_domain defaults to True."""
        config = SolverConfig()
        assert config.auto_detect_organism_domain is True

    def test_universal_constraints_regardless_of_flag(self) -> None:
        """Universal constraints are present regardless of the flag value."""
        for flag in [True, False]:
            config = SolverConfig(auto_detect_organism_domain=flag)
            model = build_csp_model(_PROTEIN, _PROKARYOTE, config)
            assert "NoATTTAMotifConstraint" in _hard_constraint_names(model), (
                f"NoATTTAMotifConstraint should be present with flag={flag}"
            )
            assert "NoTRunConstraint" in _hard_constraint_names(model), (
                f"NoTRunConstraint should be present with flag={flag}"
            )


# ─── Constraint count verification ──────────────────────────────────


class TestConstraintCounts:
    """Verify that prokaryotic models have fewer constraints than eukaryotic."""

    def test_prokaryote_has_fewer_hard_constraints(self) -> None:
        """Prokaryotic model should have fewer hard constraints than eukaryotic."""
        model_ecoli = build_csp_model(_PROTEIN, _PROKARYOTE)
        model_human = build_csp_model(_PROTEIN, "Homo_sapiens")
        assert model_ecoli.num_hard_constraints < model_human.num_hard_constraints, (
            f"E. coli hard constraints ({model_ecoli.num_hard_constraints}) "
            f"should be < human ({model_human.num_hard_constraints})"
        )

    def test_prokaryote_has_fewer_soft_constraints(self) -> None:
        """Prokaryotic model should have fewer soft constraints than eukaryotic."""
        model_ecoli = build_csp_model(_PROTEIN, _PROKARYOTE)
        model_human = build_csp_model(_PROTEIN, "Homo_sapiens")
        assert model_ecoli.num_soft_constraints < model_human.num_soft_constraints, (
            f"E. coli soft constraints ({model_ecoli.num_soft_constraints}) "
            f"should be < human ({model_human.num_soft_constraints})"
        )

    def test_exact_constraint_difference(self) -> None:
        """The difference should be exactly the eukaryote-only constraints.

        Eukaryotes have 3 additional constraints:
        - NoCrypticSpliceConstraint (hard)
        - NoCpGIslandConstraint (hard)
        - MinimizeCpG (soft)
        """
        model_ecoli = build_csp_model(_PROTEIN, _PROKARYOTE)
        model_human = build_csp_model(_PROTEIN, "Homo_sapiens")
        hard_diff = model_human.num_hard_constraints - model_ecoli.num_hard_constraints
        soft_diff = model_human.num_soft_constraints - model_ecoli.num_soft_constraints
        assert hard_diff == 2, (
            f"Expected 2 fewer hard constraints for E. coli, got {hard_diff}"
        )
        assert soft_diff == 1, (
            f"Expected 1 fewer soft constraint for E. coli, got {soft_diff}"
        )


# ─── is_eukaryotic_organism helper integration ──────────────────────


class TestIsEukaryoticOrganismIntegration:
    """Verify that is_eukaryotic_organism works with solver organism names."""

    def test_escherichia_coli_is_prokaryote(self) -> None:
        """Escherichia_coli (the solver's name) should resolve as prokaryote."""
        from biocompiler.organism_config import is_eukaryotic_organism
        assert is_eukaryotic_organism("Escherichia_coli") is False

    def test_homo_sapiens_is_eukaryote(self) -> None:
        """Homo_sapiens should resolve as eukaryote."""
        from biocompiler.organism_config import is_eukaryotic_organism
        assert is_eukaryotic_organism("Homo_sapiens") is True

    def test_ecoli_alias_is_prokaryote(self) -> None:
        """The 'ecoli' alias should resolve as prokaryote."""
        from biocompiler.organism_config import is_eukaryotic_organism
        assert is_eukaryotic_organism("ecoli") is False

    def test_e_coli_alias_is_prokaryote(self) -> None:
        """The 'E_coli' alias should resolve as prokaryote."""
        from biocompiler.organism_config import is_eukaryotic_organism
        assert is_eukaryotic_organism("E_coli") is False

    def test_unknown_organism_defaults_eukaryote(self) -> None:
        """Unknown organisms should default to eukaryote (safe fallback)."""
        from biocompiler.organism_config import is_eukaryotic_organism
        assert is_eukaryotic_organism("Unknown_organism_xyz") is True
