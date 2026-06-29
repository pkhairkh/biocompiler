"""Tests for biocompiler.organisms.config module.

Covers:
1. OrganismConfig construction, frozen immutability
2. ORGANISM_CONFIGS has all expected entries
3. get_organism_config() direct/alias/fallback lookup
4. GC targets are reasonable
5. preferred_codons populated correctly
"""

from __future__ import annotations

import dataclasses

import pytest

from biocompiler.organisms.config import (
    OrganismConfig,
    ORGANISM_CONFIGS,
    get_organism_config,
)


# ─── Constants ────────────────────────────────────────────────────────

EXPECTED_CONFIG_KEYS = {
    "E_coli_K12",
    "E_coli_BL21",
    "Homo_sapiens",
    "Saccharomyces_cerevisiae",
    "Mus_musculus",
    "CHO_K1",
}

# Standard 20 amino acid single-letter codes (no stop)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Known alias -> canonical key mapping (mirrors _ALIASES in organism_config.py)
KNOWN_ALIASES = {
    "Escherichia_coli": "E_coli_K12",
    "ecoli": "E_coli_K12",
    "E_coli": "E_coli_K12",
    "human": "Homo_sapiens",
    "mouse": "Mus_musculus",
    "cho": "CHO_K1",
    "yeast": "Saccharomyces_cerevisiae",
}


# ===========================================================================
# 1. OrganismConfig construction and frozen immutability
# ===========================================================================


class TestOrganismConfigConstruction:
    """Tests for OrganismConfig dataclass construction and immutability."""

    def test_construct_with_required_fields(self):
        """OrganismConfig can be constructed with all required fields."""
        cfg = OrganismConfig(
            name="Test Organism",
            gc_target_lo=0.40,
            gc_target_hi=0.60,
            codon_usage_validated=True,
            rbs_calculator_available=False,
        )
        assert cfg.name == "Test Organism"
        assert cfg.gc_target_lo == 0.40
        assert cfg.gc_target_hi == 0.60
        assert cfg.codon_usage_validated is True
        assert cfg.rbs_calculator_available is False

    def test_default_values(self):
        """Optional fields have correct defaults."""
        cfg = OrganismConfig(
            name="X",
            gc_target_lo=0.30,
            gc_target_hi=0.70,
            codon_usage_validated=False,
            rbs_calculator_available=False,
        )
        assert cfg.preferred_codons == {}
        assert cfg.avoided_motifs == []
        assert cfg.max_homopolymer_run == 6
        assert cfg.mrna_degradation_model == "none"

    def test_custom_optional_fields(self):
        """Optional fields can be overridden at construction."""
        cfg = OrganismConfig(
            name="X",
            gc_target_lo=0.30,
            gc_target_hi=0.70,
            codon_usage_validated=True,
            rbs_calculator_available=True,
            preferred_codons={"L": "CTG"},
            avoided_motifs=["GAATTC"],
            max_homopolymer_run=4,
            mrna_degradation_model="detailed",
        )
        assert cfg.preferred_codons == {"L": "CTG"}
        assert cfg.avoided_motifs == ["GAATTC"]
        assert cfg.max_homopolymer_run == 4
        assert cfg.mrna_degradation_model == "detailed"

    def test_frozen_immutability_name(self):
        """Cannot set attributes on a frozen dataclass instance."""
        cfg = OrganismConfig(
            name="X",
            gc_target_lo=0.30,
            gc_target_hi=0.70,
            codon_usage_validated=False,
            rbs_calculator_available=False,
        )
        with pytest.raises(dataclasses.FrozenInstanceError):
            cfg.name = "Y"

    def test_frozen_immutability_gc_lo(self):
        """Cannot mutate gc_target_lo on a frozen instance."""
        cfg = OrganismConfig(
            name="X",
            gc_target_lo=0.30,
            gc_target_hi=0.70,
            codon_usage_validated=False,
            rbs_calculator_available=False,
        )
        with pytest.raises(dataclasses.FrozenInstanceError):
            cfg.gc_target_lo = 0.50

    def test_frozen_immutability_preferred_codons(self):
        """Cannot mutate preferred_codons on a frozen instance."""
        cfg = OrganismConfig(
            name="X",
            gc_target_lo=0.30,
            gc_target_hi=0.70,
            codon_usage_validated=False,
            rbs_calculator_available=False,
            preferred_codons={"M": "ATG"},
        )
        with pytest.raises(dataclasses.FrozenInstanceError):
            cfg.preferred_codons = {"L": "CTG"}

    def test_is_dataclass_instance(self):
        """OrganismConfig instances are proper dataclass instances."""
        cfg = OrganismConfig(
            name="X",
            gc_target_lo=0.30,
            gc_target_hi=0.70,
            codon_usage_validated=False,
            rbs_calculator_available=False,
        )
        assert dataclasses.is_dataclass(cfg)
        assert dataclasses.is_dataclass(OrganismConfig)


# ===========================================================================
# 2. ORGANISM_CONFIGS registry has all expected entries
# ===========================================================================


class TestOrganismConfigsRegistry:
    """Tests for the ORGANISM_CONFIGS built-in registry."""

    def test_has_expected_keys(self):
        """ORGANISM_CONFIGS must contain all six expected organism keys."""
        assert EXPECTED_CONFIG_KEYS == set(ORGANISM_CONFIGS.keys())

    def test_all_values_are_organism_config(self):
        """Every value in ORGANISM_CONFIGS must be an OrganismConfig instance."""
        for key, cfg in ORGANISM_CONFIGS.items():
            assert isinstance(cfg, OrganismConfig), (
                f"ORGANISM_CONFIGS[{key!r}] is not an OrganismConfig instance"
            )

    def test_e_coli_k12_fields(self):
        """E_coli_K12 config has expected field values."""
        cfg = ORGANISM_CONFIGS["E_coli_K12"]
        assert cfg.name == "Escherichia coli K-12"
        assert cfg.gc_target_lo == pytest.approx(0.45)
        assert cfg.gc_target_hi == pytest.approx(0.55)
        assert cfg.codon_usage_validated is True
        assert cfg.rbs_calculator_available is False
        assert cfg.max_homopolymer_run == 5
        assert cfg.mrna_degradation_model == "simple"

    def test_e_coli_bl21_fields(self):
        """E_coli_BL21 config has expected field values."""
        cfg = ORGANISM_CONFIGS["E_coli_BL21"]
        assert cfg.name == "Escherichia coli BL21(DE3)"
        assert cfg.gc_target_lo == pytest.approx(0.45)
        assert cfg.gc_target_hi == pytest.approx(0.55)
        assert cfg.codon_usage_validated is True
        assert cfg.rbs_calculator_available is True  # BL21 has RBS calculator
        assert cfg.max_homopolymer_run == 5

    def test_homo_sapiens_fields(self):
        """Homo_sapiens config has expected field values."""
        cfg = ORGANISM_CONFIGS["Homo_sapiens"]
        assert cfg.name == "Homo sapiens"
        assert cfg.gc_target_lo == pytest.approx(0.40)
        assert cfg.gc_target_hi == pytest.approx(0.60)
        assert cfg.codon_usage_validated is True
        assert cfg.rbs_calculator_available is False
        assert cfg.mrna_degradation_model == "detailed"

    def test_saccharomyces_cerevisiae_fields(self):
        """Saccharomyces_cerevisiae config has expected field values."""
        cfg = ORGANISM_CONFIGS["Saccharomyces_cerevisiae"]
        assert cfg.name == "Saccharomyces cerevisiae"
        assert cfg.gc_target_lo == pytest.approx(0.35)
        assert cfg.gc_target_hi == pytest.approx(0.45)
        assert cfg.codon_usage_validated is True
        assert cfg.mrna_degradation_model == "simple"

    def test_mus_musculus_fields(self):
        """Mus_musculus config has expected field values."""
        cfg = ORGANISM_CONFIGS["Mus_musculus"]
        assert cfg.name == "Mus musculus"
        assert cfg.gc_target_lo == pytest.approx(0.40)
        assert cfg.gc_target_hi == pytest.approx(0.55)
        assert cfg.codon_usage_validated is False
        assert cfg.rbs_calculator_available is False

    def test_cho_k1_fields(self):
        """CHO_K1 config has expected field values."""
        cfg = ORGANISM_CONFIGS["CHO_K1"]
        assert cfg.name == "CHO-K1 (Cricetulus griseus)"
        assert cfg.gc_target_lo == pytest.approx(0.40)
        assert cfg.gc_target_hi == pytest.approx(0.60)
        assert cfg.codon_usage_validated is False

    def test_avoided_motifs_non_empty(self):
        """Each config should have at least one avoided motif."""
        for key, cfg in ORGANISM_CONFIGS.items():
            assert len(cfg.avoided_motifs) > 0, (
                f"ORGANISM_CONFIGS[{key!r}] has empty avoided_motifs"
            )

    def test_e_coli_avoided_motifs_include_ecori(self):
        """E. coli configs must include the EcoRI restriction site."""
        for key in ("E_coli_K12", "E_coli_BL21"):
            cfg = ORGANISM_CONFIGS[key]
            assert "GAATTC" in cfg.avoided_motifs, (
                f"{key} avoided_motifs missing EcoRI site 'GAATTC'"
            )

    def test_human_avoided_motifs_include_attta(self):
        """Homo_sapiens must include the ATTTA instability motif."""
        cfg = ORGANISM_CONFIGS["Homo_sapiens"]
        assert "ATTTA" in cfg.avoided_motifs

    def test_organism_names_non_empty(self):
        """Every config name must be a non-empty string."""
        for key, cfg in ORGANISM_CONFIGS.items():
            assert isinstance(cfg.name, str) and len(cfg.name) > 0, (
                f"ORGANISM_CONFIGS[{key!r}] has empty or non-string name"
            )


# ===========================================================================
# 3. get_organism_config() direct/alias/fallback lookup
# ===========================================================================


class TestGetOrganismConfig:
    """Tests for the get_organism_config() lookup helper."""

    # ── Direct hits ──────────────────────────────────────────────

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_direct_lookup_returns_config(self, key: str):
        """Direct key lookup returns the matching OrganismConfig."""
        cfg = get_organism_config(key)
        assert isinstance(cfg, OrganismConfig)
        assert cfg is ORGANISM_CONFIGS[key]

    def test_direct_lookup_e_coli_k12(self):
        """Direct lookup for E_coli_K12 returns correct config."""
        cfg = get_organism_config("E_coli_K12")
        assert cfg.name == "Escherichia coli K-12"
        assert cfg.gc_target_lo == pytest.approx(0.45)

    def test_direct_lookup_homo_sapiens(self):
        """Direct lookup for Homo_sapiens returns correct config."""
        cfg = get_organism_config("Homo_sapiens")
        assert cfg.name == "Homo sapiens"
        assert cfg.mrna_degradation_model == "detailed"

    # ── Alias resolution ─────────────────────────────────────────

    @pytest.mark.parametrize("alias,canonical", sorted(KNOWN_ALIASES.items()))
    def test_alias_lookup_returns_canonical(self, alias: str, canonical: str):
        """Alias keys resolve to their canonical OrganismConfig."""
        cfg = get_organism_config(alias)
        assert cfg is ORGANISM_CONFIGS[canonical]

    def test_ecoli_alias_resolves_to_k12(self):
        """'ecoli' alias resolves to E_coli_K12."""
        cfg = get_organism_config("ecoli")
        assert cfg.name == "Escherichia coli K-12"

    def test_human_alias_resolves(self):
        """'human' alias resolves to Homo_sapiens."""
        cfg = get_organism_config("human")
        assert cfg.name == "Homo sapiens"

    def test_mouse_alias_resolves(self):
        """'mouse' alias resolves to Mus_musculus."""
        cfg = get_organism_config("mouse")
        assert cfg.name == "Mus musculus"

    def test_cho_alias_resolves(self):
        """'cho' alias resolves to CHO_K1."""
        cfg = get_organism_config("cho")
        assert cfg.name == "CHO-K1 (Cricetulus griseus)"

    def test_yeast_alias_resolves(self):
        """'yeast' alias resolves to Saccharomyces_cerevisiae."""
        cfg = get_organism_config("yeast")
        assert cfg.name == "Saccharomyces cerevisiae"

    def test_escherichia_coli_alias(self):
        """'Escherichia_coli' alias resolves to E_coli_K12."""
        cfg = get_organism_config("Escherichia_coli")
        assert cfg.name == "Escherichia coli K-12"

    def test_e_coli_alias(self):
        """'E_coli' alias resolves to E_coli_K12."""
        cfg = get_organism_config("E_coli")
        assert cfg.name == "Escherichia coli K-12"

    # ── Fallback for unknown organisms ───────────────────────────

    def test_unknown_organism_returns_fallback(self):
        """Unknown organism key returns domain-appropriate fallback config."""
        cfg = get_organism_config("Nonexistent_organism")
        assert isinstance(cfg, OrganismConfig)
        # The fallback config name includes the domain
        assert "Unknown" in cfg.name and "fallback" in cfg.name

    def test_fallback_gc_targets_permissive(self):
        """Fallback config has permissive GC bounds [0.30, 0.70]."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.gc_target_lo == pytest.approx(0.30)
        assert cfg.gc_target_hi == pytest.approx(0.70)

    def test_fallback_codon_usage_not_validated(self):
        """Fallback config has codon_usage_validated=False."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.codon_usage_validated is False

    def test_fallback_no_rbs_calculator(self):
        """Fallback config has rbs_calculator_available=False."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.rbs_calculator_available is False

    def test_fallback_empty_preferred_codons(self):
        """Fallback config has empty preferred_codons."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.preferred_codons == {}

    def test_fallback_empty_avoided_motifs(self):
        """Fallback config has empty avoided_motifs."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.avoided_motifs == []

    def test_fallback_none_degradation_model(self):
        """Fallback config has mrna_degradation_model='none'."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.mrna_degradation_model == "none"

    def test_fallback_max_homopolymer_run(self):
        """Fallback config has max_homopolymer_run=6."""
        cfg = get_organism_config("Nonexistent_organism")
        assert cfg.max_homopolymer_run == 6

    # ── Repeated lookups are consistent ──────────────────────────

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_repeated_lookup_returns_same_object(self, key: str):
        """Repeated lookups return the exact same object (identity)."""
        cfg1 = get_organism_config(key)
        cfg2 = get_organism_config(key)
        assert cfg1 is cfg2

    def test_repeated_unknown_returns_same_fallback(self):
        """Repeated unknown lookups return the same domain-appropriate fallback object."""
        # Both are unknown eukaryotes (no bacterial indicators)
        cfg1 = get_organism_config("Xyz_organism")
        cfg2 = get_organism_config("Xyz_organism")
        assert cfg1 is cfg2


# ===========================================================================
# 4. GC targets are reasonable
# ===========================================================================


class TestGCTargetsReasonable:
    """Tests that GC target values are biologically plausible."""

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_gc_lo_in_valid_range(self, key: str):
        """GC low target must be >= 0.1 for all organisms."""
        cfg = ORGANISM_CONFIGS[key]
        assert cfg.gc_target_lo >= 0.1, (
            f"{key} gc_target_lo={cfg.gc_target_lo} is below 0.1"
        )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_gc_hi_in_valid_range(self, key: str):
        """GC high target must be <= 0.9 for all organisms."""
        cfg = ORGANISM_CONFIGS[key]
        assert cfg.gc_target_hi <= 0.9, (
            f"{key} gc_target_hi={cfg.gc_target_hi} is above 0.9"
        )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_gc_lo_strictly_below_hi(self, key: str):
        """GC low must be strictly below GC high."""
        cfg = ORGANISM_CONFIGS[key]
        assert cfg.gc_target_lo < cfg.gc_target_hi, (
            f"{key} gc_target_lo ({cfg.gc_target_lo}) not < gc_target_hi ({cfg.gc_target_hi})"
        )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_gc_window_reasonable_width(self, key: str):
        """GC target window width should be at least 0.05 (5%)."""
        cfg = ORGANISM_CONFIGS[key]
        width = cfg.gc_target_hi - cfg.gc_target_lo
        assert width >= 0.05, (
            f"{key} GC window width {width:.3f} is too narrow (< 0.05)"
        )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_gc_window_not_too_wide(self, key: str):
        """GC target window width should not exceed 0.40 (40%)."""
        cfg = ORGANISM_CONFIGS[key]
        width = cfg.gc_target_hi - cfg.gc_target_lo
        assert width <= 0.40, (
            f"{key} GC window width {width:.3f} is too wide (> 0.40)"
        )

    def test_gc_targets_vary_by_organism(self):
        """Not all organisms should have identical GC targets."""
        targets = set(
            (cfg.gc_target_lo, cfg.gc_target_hi)
            for cfg in ORGANISM_CONFIGS.values()
        )
        assert len(targets) > 1, "All organisms have identical GC targets"

    def test_yeast_has_lowest_gc_target(self):
        """Saccharomyces_cerevisiae should have the lowest GC target (AT-rich genome)."""
        yeast_lo = ORGANISM_CONFIGS["Saccharomyces_cerevisiae"].gc_target_lo
        for key, cfg in ORGANISM_CONFIGS.items():
            assert yeast_lo <= cfg.gc_target_lo, (
                f"Saccharomyces_cerevisiae gc_target_lo ({yeast_lo}) > "
                f"{key} gc_target_lo ({cfg.gc_target_lo})"
            )

    def test_e_coli_gc_targets_centered_near_50_percent(self):
        """E. coli GC targets should be centered near 50% (genome GC ~50.8%)."""
        for key in ("E_coli_K12", "E_coli_BL21"):
            cfg = ORGANISM_CONFIGS[key]
            midpoint = (cfg.gc_target_lo + cfg.gc_target_hi) / 2.0
            assert 0.45 <= midpoint <= 0.55, (
                f"{key} GC midpoint {midpoint:.3f} is not near 0.50"
            )


# ===========================================================================
# 5. preferred_codons populated correctly
# ===========================================================================


class TestPreferredCodons:
    """Tests for preferred_codons population across organism configs."""

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_preferred_codons_is_dict(self, key: str):
        """preferred_codons must be a dict for every organism."""
        cfg = ORGANISM_CONFIGS[key]
        assert isinstance(cfg.preferred_codons, dict), (
            f"{key} preferred_codons is not a dict"
        )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_preferred_codons_keys_are_standard_aa(self, key: str):
        """All keys in preferred_codons must be standard single-letter AA codes."""
        cfg = ORGANISM_CONFIGS[key]
        for aa in cfg.preferred_codons:
            assert aa in STANDARD_AA, (
                f"{key} preferred_codons has non-standard AA key {aa!r}"
            )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_preferred_codons_values_are_valid_codons(self, key: str):
        """All values in preferred_codons must be valid 3-letter DNA codons."""
        from biocompiler.type_system import CODON_TABLE

        cfg = ORGANISM_CONFIGS[key]
        for aa, codon in cfg.preferred_codons.items():
            assert len(codon) == 3, (
                f"{key} preferred_codons[{aa!r}] = {codon!r} is not 3 letters"
            )
            assert codon in CODON_TABLE, (
                f"{key} preferred_codons[{aa!r}] = {codon!r} is not a valid codon"
            )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_preferred_codons_match_their_amino_acid(self, key: str):
        """Each preferred codon must translate to its corresponding amino acid."""
        from biocompiler.type_system import CODON_TABLE

        cfg = ORGANISM_CONFIGS[key]
        for aa, codon in cfg.preferred_codons.items():
            translated = CODON_TABLE[codon]
            assert translated == aa, (
                f"{key} preferred_codons[{aa!r}] = {codon!r} translates to "
                f"{translated!r}, not {aa!r}"
            )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_preferred_codons_covers_most_amino_acids(self, key: str):
        """preferred_codons should cover at least 18 of the 20 standard AAs.

        Some organisms may not have preferred codons for amino acids with
        a single codon (M=ATG, W=TGG), but most should be covered.
        """
        cfg = ORGANISM_CONFIGS[key]
        covered = set(cfg.preferred_codons.keys())
        # Methionine (M) and Tryptophan (W) each have only one codon,
        # so some tables may omit them as trivially preferred.
        required = STANDARD_AA - {"M", "W"}
        missing = required - covered
        assert len(missing) <= 2, (
            f"{key} preferred_codons is missing too many AAs: {sorted(missing)}"
        )

    @pytest.mark.parametrize("key", sorted(EXPECTED_CONFIG_KEYS))
    def test_preferred_codons_non_empty(self, key: str):
        """preferred_codons must be non-empty for all organisms."""
        cfg = ORGANISM_CONFIGS[key]
        assert len(cfg.preferred_codons) > 0, (
            f"{key} preferred_codons is empty"
        )

    def test_e_coli_preferred_leucine_is_ctg(self):
        """E. coli preferred codon for Leucine should be CTG (highest frequency)."""
        cfg = ORGANISM_CONFIGS["E_coli_K12"]
        if "L" in cfg.preferred_codons:
            assert cfg.preferred_codons["L"] == "CTG", (
                f"E_coli_K12 preferred L codon is {cfg.preferred_codons['L']!r}, expected 'CTG'"
            )

    def test_human_preferred_leucine_is_ctg(self):
        """Human preferred codon for Leucine should be CTG (highest frequency)."""
        cfg = ORGANISM_CONFIGS["Homo_sapiens"]
        if "L" in cfg.preferred_codons:
            assert cfg.preferred_codons["L"] == "CTG"

    def test_e_coli_bl21_shares_preferred_codons_with_k12(self):
        """E_coli_BL21 and E_coli_K12 share the same preferred_codons."""
        k12 = ORGANISM_CONFIGS["E_coli_K12"]
        bl21 = ORGANISM_CONFIGS["E_coli_BL21"]
        assert k12.preferred_codons == bl21.preferred_codons, (
            "E_coli_K12 and E_coli_BL21 preferred_codons differ"
        )

    def test_no_stop_codon_in_preferred_codons(self):
        """No preferred_codons entry should map to a stop codon."""
        from biocompiler.type_system import CODON_TABLE

        for key, cfg in ORGANISM_CONFIGS.items():
            for aa, codon in cfg.preferred_codons.items():
                translated = CODON_TABLE[codon]
                assert translated != "*", (
                    f"{key} preferred_codons[{aa!r}] = {codon!r} is a stop codon"
                )
