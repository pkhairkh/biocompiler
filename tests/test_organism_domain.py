"""Tests for organism domain classification (F1.1).

Validates that each organism in the registry is correctly classified
as eukaryotic or prokaryotic, and that the ``is_eukaryotic_organism``
helper returns the expected values — including the safe default for
unknown organisms.
"""

from __future__ import annotations

import pytest

from biocompiler.organism_config import (
    ORGANISM_CONFIGS,
    OrganismConfig,
    get_organism_config,
    is_eukaryotic_organism,
)


# ─── Registry-level domain classification ───────────────────────────


class TestRegistryDomainClassification:
    """Each entry in ORGANISM_CONFIGS should have the correct domain."""

    # Prokaryotes
    @pytest.mark.parametrize("key", ["E_coli_K12", "E_coli_BL21"])
    def test_ecoli_is_prokaryote(self, key: str) -> None:
        cfg = ORGANISM_CONFIGS[key]
        assert cfg.domain == "prokaryote", (
            f"{key} should be prokaryote, got {cfg.domain!r}"
        )
        assert not cfg.is_eukaryote, (
            f"{key}.is_eukaryote should be False"
        )

    # Eukaryotes
    @pytest.mark.parametrize(
        "key",
        ["Homo_sapiens", "Saccharomyces_cerevisiae", "Mus_musculus", "CHO_K1"],
    )
    def test_eukaryotic_organisms(self, key: str) -> None:
        cfg = ORGANISM_CONFIGS[key]
        assert cfg.domain == "eukaryote", (
            f"{key} should be eukaryote, got {cfg.domain!r}"
        )
        assert cfg.is_eukaryote, (
            f"{key}.is_eukaryote should be True"
        )


# ─── Fallback config ───────────────────────────────────────────────


class TestFallbackConfig:
    """The fallback config should default to eukaryote (safe default)."""

    def test_fallback_is_eukaryote(self) -> None:
        cfg = get_organism_config("Totally_unknown_organism_xyz")
        assert cfg.domain == "eukaryote"
        assert cfg.is_eukaryote is True

    def test_fallback_name(self) -> None:
        cfg = get_organism_config("does_not_exist")
        assert "fallback" in cfg.name.lower()


# ─── is_eukaryotic_organism helper ─────────────────────────────────


class TestIsEukaryoticOrganism:
    """The module-level helper should work for all known organisms."""

    @pytest.mark.parametrize(
        "key",
        ["Homo_sapiens", "Saccharomyces_cerevisiae", "Mus_musculus", "CHO_K1"],
    )
    def test_eukaryotes_return_true(self, key: str) -> None:
        assert is_eukaryotic_organism(key) is True

    @pytest.mark.parametrize("key", ["E_coli_K12", "E_coli_BL21"])
    def test_prokaryotes_return_false(self, key: str) -> None:
        assert is_eukaryotic_organism(key) is False

    def test_unknown_organism_returns_true(self) -> None:
        """Unknown organisms fall back to eukaryote — the safe default."""
        assert is_eukaryotic_organism("Unknown_organism_42") is True

    @pytest.mark.parametrize(
        "alias,expected",
        [
            ("human", True),
            ("mouse", True),
            ("yeast", True),
            ("cho", True),
            ("ecoli", False),
            ("E_coli", False),
            ("Escherichia_coli", False),
        ],
    )
    def test_legacy_aliases(self, alias: str, expected: bool) -> None:
        """Legacy aliases should resolve correctly through get_organism_config."""
        assert is_eukaryotic_organism(alias) is expected


# ─── OrganismConfig.is_eukaryote property ──────────────────────────


class TestIsEukaryoteProperty:
    """Direct tests on the OrganismConfig.is_eukaryote property."""

    def test_eukaryote_property(self) -> None:
        cfg = OrganismConfig(
            name="Test eukaryote",
            gc_target_lo=0.4,
            gc_target_hi=0.6,
            codon_usage_validated=False,
            rbs_calculator_available=False,
            domain="eukaryote",
        )
        assert cfg.is_eukaryote is True

    def test_prokaryote_property(self) -> None:
        cfg = OrganismConfig(
            name="Test prokaryote",
            gc_target_lo=0.5,
            gc_target_hi=0.5,
            codon_usage_validated=False,
            rbs_calculator_available=False,
            domain="prokaryote",
        )
        assert cfg.is_eukaryote is False

    def test_archaea_property(self) -> None:
        cfg = OrganismConfig(
            name="Test archaeon",
            gc_target_lo=0.5,
            gc_target_hi=0.5,
            codon_usage_validated=False,
            rbs_calculator_available=False,
            domain="archaea",
        )
        assert cfg.is_eukaryote is False

    def test_default_domain_is_eukaryote(self) -> None:
        """Omitting domain should default to 'eukaryote'."""
        cfg = OrganismConfig(
            name="Default domain test",
            gc_target_lo=0.4,
            gc_target_hi=0.6,
            codon_usage_validated=False,
            rbs_calculator_available=False,
        )
        assert cfg.domain == "eukaryote"
        assert cfg.is_eukaryote is True
