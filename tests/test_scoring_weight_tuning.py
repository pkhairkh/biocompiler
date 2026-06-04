"""Tests for F2.8: Scoring weight tuning.

Validates that SolverConfig.for_organism() returns properly tuned configs
for each organism and that organism-specific weight adjustments work correctly.

Key biological rationale:
- Prokaryotes lack DNA methylation → CpG weight should be 0
- Eukaryotes with CpG methylation → non-zero CpG weight
- mRNA dG weight varies by degradation model fidelity
"""

from __future__ import annotations

import importlib
import sys
import pytest


# ── Direct module loading to avoid circular import issues ────────────

def _load_module_directly(module_name: str, file_path: str):
    """Load a module directly from file, bypassing package __init__.py."""
    # Set up biocompiler package if not already present
    if "biocompiler" not in sys.modules or not hasattr(sys.modules.get("biocompiler"), "__path__"):
        _pkg = type(sys)("biocompiler")
        _pkg.__path__ = ["src/biocompiler"]
        _pkg.__package__ = "biocompiler"
        sys.modules["biocompiler"] = _pkg

    if module_name in sys.modules:
        return sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(module_name, file_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod


_types_mod = _load_module_directly(
    "biocompiler.solver.types",
    "src/biocompiler/solver/types.py",
)
_organism_config_mod = _load_module_directly(
    "biocompiler.organism_config",
    "src/biocompiler/organism_config.py",
)

SolverConfig = _types_mod.SolverConfig
OrganismConfig = _organism_config_mod.OrganismConfig
ORGANISM_CONFIGS = _organism_config_mod.ORGANISM_CONFIGS
get_organism_config = _organism_config_mod.get_organism_config
is_eukaryotic_organism = _organism_config_mod.is_eukaryotic_organism


# ════════════════════════════════════════════════════════════════════
# 1. SolverConfig.for_organism() tests
# ════════════════════════════════════════════════════════════════════

class TestSolverConfigForOrganism:
    """Tests for SolverConfig.for_organism() class method."""

    # ── E. coli (prokaryote) ────────────────────────────────────────

    def test_ecoli_k12_cai_weight(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        assert cfg.cai_weight == 1.0

    def test_ecoli_k12_cpg_weight_zero(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        assert cfg.cpg_weight == 0.0

    def test_ecoli_k12_mrna_dg_weight(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        assert cfg.mrna_dg_weight == 0.2

    def test_ecoli_k12_avoid_cpg_false(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        assert cfg.avoid_cpg is False

    def test_ecoli_k12_gc_range(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        assert cfg.gc_lo == 0.45
        assert cfg.gc_hi == 0.55

    def test_ecoli_k12_splice_threshold_zero(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        assert cfg.cryptic_splice_threshold == 0.0

    def test_ecoli_bl21_matches_k12_weights(self):
        cfg_k12 = SolverConfig.for_organism("E_coli_K12")
        cfg_bl21 = SolverConfig.for_organism("E_coli_BL21")
        assert cfg_bl21.cai_weight == cfg_k12.cai_weight
        assert cfg_bl21.cpg_weight == cfg_k12.cpg_weight
        assert cfg_bl21.mrna_dg_weight == cfg_k12.mrna_dg_weight

    # ── Human (eukaryote) ──────────────────────────────────────────

    def test_human_cai_weight(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.cai_weight == 1.0

    def test_human_cpg_weight_nonzero(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.cpg_weight == 0.5

    def test_human_mrna_dg_weight(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.mrna_dg_weight == 0.3

    def test_human_avoid_cpg_true(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.avoid_cpg is True

    def test_human_gc_range(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.gc_lo == 0.40
        assert cfg.gc_hi == 0.60

    def test_human_splice_threshold(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.cryptic_splice_threshold == 3.0

    # ── Mouse (eukaryote) ──────────────────────────────────────────

    def test_mouse_cpg_weight(self):
        cfg = SolverConfig.for_organism("Mus_musculus")
        assert cfg.cpg_weight == 0.4

    def test_mouse_mrna_dg_weight(self):
        cfg = SolverConfig.for_organism("Mus_musculus")
        assert cfg.mrna_dg_weight == 0.2

    def test_mouse_gc_range(self):
        cfg = SolverConfig.for_organism("Mus_musculus")
        assert cfg.gc_lo == 0.40
        assert cfg.gc_hi == 0.55

    # ── CHO (eukaryote) ────────────────────────────────────────────

    def test_cho_cpg_weight(self):
        cfg = SolverConfig.for_organism("CHO_K1")
        assert cfg.cpg_weight == 0.5

    def test_cho_mrna_dg_weight(self):
        cfg = SolverConfig.for_organism("CHO_K1")
        assert cfg.mrna_dg_weight == 0.3

    # ── Yeast (eukaryote with low CpG concern) ─────────────────────

    def test_yeast_cpg_weight_low(self):
        cfg = SolverConfig.for_organism("Saccharomyces_cerevisiae")
        assert cfg.cpg_weight == 0.1

    def test_yeast_mrna_dg_weight(self):
        cfg = SolverConfig.for_organism("Saccharomyces_cerevisiae")
        assert cfg.mrna_dg_weight == 0.2

    def test_yeast_avoid_cpg_false(self):
        cfg = SolverConfig.for_organism("Saccharomyces_cerevisiae")
        assert cfg.avoid_cpg is False

    # ── Alias resolution ───────────────────────────────────────────

    def test_alias_ecoli(self):
        cfg = SolverConfig.for_organism("ecoli")
        cfg_direct = SolverConfig.for_organism("E_coli_K12")
        assert cfg.cpg_weight == cfg_direct.cpg_weight

    def test_alias_human(self):
        cfg = SolverConfig.for_organism("human")
        cfg_direct = SolverConfig.for_organism("Homo_sapiens")
        assert cfg.cpg_weight == cfg_direct.cpg_weight

    def test_alias_mouse(self):
        cfg = SolverConfig.for_organism("mouse")
        cfg_direct = SolverConfig.for_organism("Mus_musculus")
        assert cfg.cpg_weight == cfg_direct.cpg_weight

    # ── Unknown organism fallback ──────────────────────────────────

    def test_unknown_organism_returns_default(self):
        cfg = SolverConfig.for_organism("Unknown_organism_xyz")
        default = SolverConfig()
        assert cfg.cai_weight == default.cai_weight
        assert cfg.cpg_weight == default.cpg_weight
        assert cfg.mrna_dg_weight == default.mrna_dg_weight

    def test_unknown_organism_is_valid_config(self):
        cfg = SolverConfig.for_organism("does_not_exist")
        assert isinstance(cfg, SolverConfig)


# ════════════════════════════════════════════════════════════════════
# 2. Prokaryotic vs Eukaryotic CpG weight classification
# ════════════════════════════════════════════════════════════════════

class TestProkaryoticCpGWeight:
    """Prokaryotic configs should always have zero CpG weight."""

    @pytest.mark.parametrize("organism", [
        "E_coli_K12",
        "E_coli_BL21",
    ])
    def test_prokaryotic_cpg_weight_is_zero(self, organism):
        cfg = SolverConfig.for_organism(organism)
        assert cfg.cpg_weight == 0.0

    @pytest.mark.parametrize("organism", [
        "E_coli_K12",
        "E_coli_BL21",
    ])
    def test_prokaryotic_avoid_cpg_is_false(self, organism):
        cfg = SolverConfig.for_organism(organism)
        assert cfg.avoid_cpg is False


class TestEukaryoticCpGWeight:
    """Eukaryotic configs should have non-zero CpG weight."""

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "Mus_musculus",
        "CHO_K1",
    ])
    def test_eukaryotic_cpg_weight_nonzero(self, organism):
        cfg = SolverConfig.for_organism(organism)
        assert cfg.cpg_weight > 0.0

    @pytest.mark.parametrize("organism", [
        "Homo_sapiens",
        "Mus_musculus",
        "CHO_K1",
    ])
    def test_eukaryotic_avoid_cpg_is_true(self, organism):
        cfg = SolverConfig.for_organism(organism)
        assert cfg.avoid_cpg is True

    def test_yeast_eukaryote_but_low_cpg(self):
        cfg = SolverConfig.for_organism("Saccharomyces_cerevisiae")
        assert cfg.cpg_weight <= 0.1
        assert cfg.avoid_cpg is False


# ════════════════════════════════════════════════════════════════════
# 3. SoftConstraintScorer.adjust_weights_for_organism() tests
# ════════════════════════════════════════════════════════════════════

class TestSoftConstraintScorerAdjustWeights:
    """Tests for SoftConstraintScorer.adjust_weights_for_organism().

    These tests verify the weight adjustment logic by directly calling
    the method on a SoftConstraintScorer instance. The scorer's config
    is modified in-place.
    """

    @pytest.fixture
    def scorer(self):
        """Create a SoftConstraintScorer with default weights."""
        # Import SoftConstraintScorer lazily to avoid circular imports
        from biocompiler.solver.scoring import SoftConstraintScorer
        config = SolverConfig(cai_weight=1.0, cpg_weight=0.5, mrna_dg_weight=0.3)
        return SoftConstraintScorer(config)

    def test_prokaryotic_adjustment_zeros_cpg(self, scorer):
        scorer.adjust_weights_for_organism("E_coli_K12")
        assert scorer.config.cpg_weight == 0.0

    def test_prokaryotic_adjustment_disables_avoid_cpg(self, scorer):
        scorer.adjust_weights_for_organism("E_coli_K12")
        assert scorer.config.avoid_cpg is False

    def test_prokaryotic_cai_weight_unchanged(self, scorer):
        scorer.adjust_weights_for_organism("E_coli_K12")
        assert scorer.config.cai_weight == 1.0

    def test_eukaryotic_cpg_weight_preserved(self, scorer):
        scorer.adjust_weights_for_organism("Homo_sapiens")
        assert scorer.config.cpg_weight > 0.0

    def test_eukaryotic_avoid_cpg_not_disabled(self, scorer):
        scorer.adjust_weights_for_organism("Homo_sapiens")
        assert scorer.config.avoid_cpg is True

    def test_ecoli_alias_adjustment(self, scorer):
        scorer.adjust_weights_for_organism("ecoli")
        assert scorer.config.cpg_weight == 0.0

    def test_gc_bounds_updated_from_organism(self, scorer):
        scorer.adjust_weights_for_organism("E_coli_K12")
        org_cfg = get_organism_config("E_coli_K12")
        assert scorer.config.gc_lo == org_cfg.gc_target_lo
        assert scorer.config.gc_hi == org_cfg.gc_target_hi

    def test_mrna_dg_weight_capped_for_simple_model(self, scorer):
        """E. coli uses simple mRNA model; weight should be capped at 0.2."""
        scorer.config.mrna_dg_weight = 0.5
        scorer.adjust_weights_for_organism("E_coli_K12")
        assert scorer.config.mrna_dg_weight <= 0.2

    def test_mrna_dg_weight_preserved_for_detailed_model(self, scorer):
        """Human uses detailed mRNA model; weight should be preserved."""
        scorer.adjust_weights_for_organism("Homo_sapiens")
        assert scorer.config.mrna_dg_weight == 0.3

    def test_adjustment_is_in_place(self, scorer):
        config_id_before = id(scorer.config)
        scorer.adjust_weights_for_organism("E_coli_K12")
        config_id_after = id(scorer.config)
        assert config_id_before == config_id_after

    def test_unknown_organism_adjustment(self, scorer):
        """Unknown organism falls back to eukaryote; cpg_weight not zeroed."""
        scorer.adjust_weights_for_organism("totally_unknown_species")
        assert scorer.config.cpg_weight == 0.5


# ════════════════════════════════════════════════════════════════════
# 4. Integration: for_organism + organism_config consistency
# ════════════════════════════════════════════════════════════════════

class TestForOrganismOrganismConfigConsistency:
    """Integration tests ensuring SolverConfig.for_organism() aligns with
    OrganismConfig parameters."""

    def test_ecoli_gc_range_matches_organism_config(self):
        cfg = SolverConfig.for_organism("E_coli_K12")
        org = get_organism_config("E_coli_K12")
        assert cfg.gc_lo == org.gc_target_lo
        assert cfg.gc_hi == org.gc_target_hi

    def test_human_gc_range_matches_organism_config(self):
        cfg = SolverConfig.for_organism("Homo_sapiens")
        org = get_organism_config("Homo_sapiens")
        assert cfg.gc_lo == org.gc_target_lo
        assert cfg.gc_hi == org.gc_target_hi

    def test_prokaryotic_organism_has_prokaryote_domain(self):
        assert not is_eukaryotic_organism("E_coli_K12")
        assert not is_eukaryotic_organism("E_coli_BL21")

    def test_eukaryotic_organisms_have_eukaryote_domain(self):
        for org in ["Homo_sapiens", "Mus_musculus", "CHO_K1", "Saccharomyces_cerevisiae"]:
            assert is_eukaryotic_organism(org), f"{org} should be eukaryotic"

    def test_all_eukaryotes_have_positive_cpg_in_presets(self):
        eukaryotes = ["Homo_sapiens", "Mus_musculus", "CHO_K1"]
        for org in eukaryotes:
            cfg = SolverConfig.for_organism(org)
            assert cfg.cpg_weight > 0.0, f"{org} should have positive cpg_weight"

    def test_all_prokaryotes_have_zero_cpg_in_presets(self):
        prokaryotes = ["E_coli_K12", "E_coli_BL21"]
        for org in prokaryotes:
            cfg = SolverConfig.for_organism(org)
            assert cfg.cpg_weight == 0.0, f"{org} should have zero cpg_weight"
