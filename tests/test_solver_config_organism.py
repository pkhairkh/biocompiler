"""Tests for SolverConfig.organism field and is_eukaryotic property.

Task F1.9: Organism field on SolverConfig so the solver always knows which
organism is being targeted, enabling organism-aware constraint selection
throughout the pipeline.
"""

from __future__ import annotations

import pytest

from biocompiler.solver.types import SolverConfig
from biocompiler.solver.constraints import build_csp_model


# ─── SolverConfig.organism default ────────────────────────────────────


class TestSolverConfigOrganismDefault:
    """SolverConfig.organism defaults to 'Homo_sapiens'."""

    def test_default_organism(self) -> None:
        cfg = SolverConfig()
        assert cfg.organism == "Homo_sapiens"

    def test_explicit_organism(self) -> None:
        cfg = SolverConfig(organism="E_coli_K12")
        assert cfg.organism == "E_coli_K12"

    def test_organism_preserved_with_other_fields(self) -> None:
        cfg = SolverConfig(organism="Saccharomyces_cerevisiae", gc_lo=0.35)
        assert cfg.organism == "Saccharomyces_cerevisiae"
        assert cfg.gc_lo == 0.35


# ─── SolverConfig.is_eukaryotic property ──────────────────────────────


class TestSolverConfigIsEukaryotic:
    """SolverConfig.is_eukaryotic correctly detects organism domain."""

    def test_eukaryotic_organism(self) -> None:
        cfg = SolverConfig(organism="Homo_sapiens")
        assert cfg.is_eukaryotic is True

    def test_prokaryotic_organism(self) -> None:
        cfg = SolverConfig(organism="E_coli_K12")
        assert cfg.is_eukaryotic is False

    def test_yeast_is_eukaryotic(self) -> None:
        cfg = SolverConfig(organism="Saccharomyces_cerevisiae")
        assert cfg.is_eukaryotic is True

    def test_mouse_is_eukaryotic(self) -> None:
        cfg = SolverConfig(organism="Mus_musculus")
        assert cfg.is_eukaryotic is True

    def test_cho_is_eukaryotic(self) -> None:
        cfg = SolverConfig(organism="CHO_K1")
        assert cfg.is_eukaryotic is True

    def test_default_organism_is_eukaryotic(self) -> None:
        cfg = SolverConfig()
        assert cfg.is_eukaryotic is True  # Homo_sapiens default


# ─── SolverConfig.auto_detect_organism_domain ─────────────────────────


class TestAutoDetectOrganismDomain:
    """auto_detect_organism_domain flag controls is_eukaryotic behaviour."""

    def test_auto_detect_on_prokaryote(self) -> None:
        """With auto_detect=True, prokaryotic organism returns False."""
        cfg = SolverConfig(organism="E_coli_K12", auto_detect_organism_domain=True)
        assert cfg.is_eukaryotic is False

    def test_auto_detect_off_prokaryote_assumes_eukaryote(self) -> None:
        """With auto_detect=False, even prokaryotic organism returns True
        (conservative default)."""
        cfg = SolverConfig(organism="E_coli_K12", auto_detect_organism_domain=False)
        assert cfg.is_eukaryotic is True

    def test_auto_detect_on_eukaryote(self) -> None:
        """With auto_detect=True, eukaryotic organism returns True."""
        cfg = SolverConfig(organism="Homo_sapiens", auto_detect_organism_domain=True)
        assert cfg.is_eukaryotic is True

    def test_auto_detect_off_eukaryote(self) -> None:
        """With auto_detect=False, eukaryotic organism still returns True."""
        cfg = SolverConfig(organism="Homo_sapiens", auto_detect_organism_domain=False)
        assert cfg.is_eukaryotic is True

    def test_default_auto_detect_is_true(self) -> None:
        cfg = SolverConfig()
        assert cfg.auto_detect_organism_domain is True


# ─── build_csp_model uses config.organism ─────────────────────────────


class TestBuildCspModelUsesConfigOrganism:
    """build_csp_model falls back to config.organism when organism is None."""

    def test_organism_from_config(self) -> None:
        """When organism param is None, config.organism is used."""
        cfg = SolverConfig(organism="Homo_sapiens")
        model = build_csp_model(protein="MVLSPADKTN", organism=None, config=cfg)
        assert model.organism == "Homo_sapiens"

    def test_organism_param_overrides_config(self) -> None:
        """Explicit organism param takes precedence over config.organism."""
        cfg = SolverConfig(organism="Homo_sapiens")
        model = build_csp_model(protein="MVLSPADKTN", organism="Escherichia_coli", config=cfg)
        assert model.organism == "Escherichia_coli"

    def test_prokaryotic_config_skips_splice_constraints(self) -> None:
        """Prokaryotic config.organism causes splice constraints to be skipped."""
        cfg = SolverConfig(organism="Escherichia_coli")
        model = build_csp_model(protein="MVLSPADKTN", organism=None, config=cfg)
        # NoCrypticSpliceConstraint should be absent for prokaryotes
        constraint_names = [c.name for c in model.hard_constraints]
        assert "NoCrypticSpliceConstraint" not in constraint_names

    def test_eukaryotic_config_includes_splice_constraints(self) -> None:
        """Eukaryotic config.organism causes splice constraints to be included."""
        cfg = SolverConfig(organism="Homo_sapiens")
        model = build_csp_model(protein="MVLSPADKTN", organism=None, config=cfg)
        constraint_names = [c.name for c in model.hard_constraints]
        assert "NoCrypticSpliceConstraint" in constraint_names

    def test_default_config_uses_homo_sapiens(self) -> None:
        """Default SolverConfig uses Homo_sapiens as organism."""
        model = build_csp_model(protein="MVLSPADKTN", organism=None)
        assert model.organism == "Homo_sapiens"
