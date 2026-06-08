"""Tests for Task 1.8: constants centralization, validation fixes, and exception logging.

Covers:
1. Constants module: all extracted constants are importable and have expected values
2. Constants module: YAML config override via BIOCOMPILER_CONFIG
3. Validation fixes: assert → ValueError in OptimizationResult, FullConstructResult,
   sliding_gc, multigene, etc.
4. Species-to-organism deduplication: BioOptimizer and HybridOptimizer use resolve_organism
5. Exception swallowing: logging added to NUMBA fallbacks and maxentscan checks
"""

from __future__ import annotations

import os
import tempfile
import pytest

# ── 1. Constants module importability ──────────────────────────────────────

class TestDeimmunizationConstants:
    """Verify deimmunization constants are centralized in constants.py."""

    def test_mhc_i_peptide_length(self):
        from biocompiler.constants import DEIMMUNIZATION_MHC_I_PEPTIDE_LENGTH
        assert DEIMMUNIZATION_MHC_I_PEPTIDE_LENGTH == 9

    def test_blosum62_range_divisor(self):
        from biocompiler.constants import DEIMMUNIZATION_BLOSUM62_RANGE_DIVISOR
        assert DEIMMUNIZATION_BLOSUM62_RANGE_DIVISOR == pytest.approx(4.0, rel=1e-6)

    def test_ddg_scaling_divisor(self):
        from biocompiler.constants import DEIMMUNIZATION_DDG_SCALING_DIVISOR
        assert DEIMMUNIZATION_DDG_SCALING_DIVISOR == pytest.approx(5.0, rel=1e-6)

    def test_default_blosum62_score(self):
        from biocompiler.constants import DEIMMUNIZATION_DEFAULT_BLOSUM62_SCORE
        assert DEIMMUNIZATION_DEFAULT_BLOSUM62_SCORE == -4

    def test_top_epitopes_for_scoring(self):
        from biocompiler.constants import DEIMMUNIZATION_TOP_EPITOPES_FOR_SCORING
        assert DEIMMUNIZATION_TOP_EPITOPES_FOR_SCORING == 5

    def test_cumulative_ddg_multiplier(self):
        from biocompiler.constants import DEIMMUNIZATION_CUMULATIVE_DDG_MULTIPLIER
        assert DEIMMUNIZATION_CUMULATIVE_DDG_MULTIPLIER == 2

    def test_max_cumulative_ddg_kcal(self):
        from biocompiler.constants import DEIMMUNIZATION_MAX_CUMULATIVE_DDG_KCAL
        assert DEIMMUNIZATION_MAX_CUMULATIVE_DDG_KCAL == pytest.approx(5.0, rel=1e-6)

    def test_anchor_positions(self):
        from biocompiler.constants import DEIMMUNIZATION_MHC_I_ANCHOR_POSITIONS
        assert DEIMMUNIZATION_MHC_I_ANCHOR_POSITIONS == {1, 8}

    def test_helix_favoring(self):
        from biocompiler.constants import DEIMMUNIZATION_HELIX_FAVORING
        assert "A" in DEIMMUNIZATION_HELIX_FAVORING
        assert "E" in DEIMMUNIZATION_HELIX_FAVORING

    def test_strong_binder_class(self):
        from biocompiler.constants import DEIMMUNIZATION_STRONG_BINDER_CLASS
        assert DEIMMUNIZATION_STRONG_BINDER_CLASS == "strong_binder"

    def test_deimmunization_imports_from_constants(self):
        """Verify deimmunization.py uses constants from constants.py."""
        from biocompiler import deimmunization
        # The module should have _MHC_I_PEPTIDE_LENGTH available (imported from constants)
        assert hasattr(deimmunization, '_MHC_I_PEPTIDE_LENGTH')
        assert deimmunization._MHC_I_PEPTIDE_LENGTH == 9


class TestCamSolConstants:
    """Verify CamSol constants are centralized in constants.py."""

    def test_weights(self):
        from biocompiler.constants import (
            CAMSOL_WEIGHT_HYDROPATHY, CAMSOL_WEIGHT_CHARGE,
            CAMSOL_WEIGHT_ALPHA_HELIX, CAMSOL_WEIGHT_BETA_STRAND,
            CAMSOL_WEIGHT_PROGLY,
        )
        assert CAMSOL_WEIGHT_HYDROPATHY == pytest.approx(0.35)
        assert CAMSOL_WEIGHT_CHARGE == pytest.approx(0.25)
        assert CAMSOL_WEIGHT_ALPHA_HELIX == pytest.approx(0.15)
        assert CAMSOL_WEIGHT_BETA_STRAND == pytest.approx(0.15)
        assert CAMSOL_WEIGHT_PROGLY == pytest.approx(0.10)
        # Weights should sum to 1.0
        total = (CAMSOL_WEIGHT_HYDROPATHY + CAMSOL_WEIGHT_CHARGE +
                 CAMSOL_WEIGHT_ALPHA_HELIX + CAMSOL_WEIGHT_BETA_STRAND +
                 CAMSOL_WEIGHT_PROGLY)
        assert total == pytest.approx(1.0)

    def test_aggregation_threshold(self):
        from biocompiler.constants import CAMSOL_AGGREGATION_THRESHOLD
        assert CAMSOL_AGGREGATION_THRESHOLD == -0.5

    def test_score_clamp(self):
        from biocompiler.constants import CAMSOL_SCORE_CLAMP_MIN, CAMSOL_SCORE_CLAMP_MAX
        assert CAMSOL_SCORE_CLAMP_MIN == -3.0
        assert CAMSOL_SCORE_CLAMP_MAX == pytest.approx(3.0, rel=1e-6)

    def test_idp_thresholds(self):
        from biocompiler.constants import (
            CAMSOL_IDP_DISORDER_FRACTION_THRESHOLD,
            CAMSOL_IDP_ORDER_FRACTION_THRESHOLD,
        )
        assert CAMSOL_IDP_DISORDER_FRACTION_THRESHOLD == pytest.approx(0.40)
        assert CAMSOL_IDP_ORDER_FRACTION_THRESHOLD == pytest.approx(0.31)

    def test_sasa_thresholds(self):
        from biocompiler.constants import CAMSOL_SASA_BURIED_THRESHOLD, CAMSOL_SASA_EXPOSED_THRESHOLD
        assert CAMSOL_SASA_BURIED_THRESHOLD < CAMSOL_SASA_EXPOSED_THRESHOLD

    def test_camsol_imports_from_constants(self):
        """Verify camsol.py uses constants from constants.py."""
        from biocompiler import camsol
        assert hasattr(camsol, '_WEIGHT_HYDROPATHY')
        assert camsol._WEIGHT_HYDROPATHY == pytest.approx(0.35)


class TestStabilityConstants:
    """Verify stability predicates constants are centralized in constants.py."""

    def test_hydro_frac_range(self):
        from biocompiler.constants import STABILITY_HYDRO_FRAC_LO, STABILITY_HYDRO_FRAC_HI
        assert STABILITY_HYDRO_FRAC_LO == pytest.approx(0.30)
        assert STABILITY_HYDRO_FRAC_HI == pytest.approx(0.45)
        assert STABILITY_HYDRO_FRAC_LO < STABILITY_HYDRO_FRAC_HI

    def test_disulfide_threshold(self):
        from biocompiler.constants import STABILITY_DISULFIDE_CB_DIST_THRESHOLD
        assert STABILITY_DISULFIDE_CB_DIST_THRESHOLD == pytest.approx(6.5)

    def test_prokaryotic_keywords(self):
        from biocompiler.constants import STABILITY_PROKARYOTIC_KEYWORDS
        assert "ecoli" in STABILITY_PROKARYOTIC_KEYWORDS
        assert "escherichia" in STABILITY_PROKARYOTIC_KEYWORDS

    def test_stability_imports_from_constants(self):
        """Verify stability_predicates.py uses constants from constants.py."""
        from biocompiler import stability_predicates
        assert hasattr(stability_predicates, '_HYDRO_FRAC_LO')
        assert stability_predicates._HYDRO_FRAC_LO == pytest.approx(0.30)


class TestSplicingConstants:
    """Verify splicing constants are centralized in constants.py."""

    def test_context_lengths(self):
        from biocompiler.constants import SPLICING_PWM_CONTEXT_LEN, SPLICING_MIN_CONTEXT_LEN
        assert SPLICING_PWM_CONTEXT_LEN == 9
        assert SPLICING_MIN_CONTEXT_LEN == 4

    def test_default_thresholds(self):
        from biocompiler.constants import SPLICING_DEFAULT_LOW_THRESH, SPLICING_DEFAULT_HIGH_THRESH
        assert SPLICING_DEFAULT_LOW_THRESH == pytest.approx(3.0, rel=1e-6)
        assert SPLICING_DEFAULT_HIGH_THRESH == pytest.approx(8.0, rel=1e-6)
        assert SPLICING_DEFAULT_LOW_THRESH < SPLICING_DEFAULT_HIGH_THRESH

    def test_max_isoforms(self):
        from biocompiler.constants import SPLICING_DEFAULT_MAX_ISOFORMS
        assert SPLICING_DEFAULT_MAX_ISOFORMS == 100

    def test_donor_acceptor_max_distance(self):
        from biocompiler.constants import SPLICING_DONOR_ACCEPTOR_MAX_DISTANCE
        assert SPLICING_DONOR_ACCEPTOR_MAX_DISTANCE == 500

    def test_splicing_imports_from_constants(self):
        """Verify splicing.py uses constants from constants.py."""
        from biocompiler import splicing
        assert hasattr(splicing, '_DEFAULT_LOW_THRESH')
        assert splicing._DEFAULT_LOW_THRESH == pytest.approx(3.0, rel=1e-6)


class TestOptimizationConstants:
    """Verify optimization constants are centralized in constants.py."""

    def test_t_run_threshold(self):
        from biocompiler.constants import OPT_T_RUN_LENGTH_THRESHOLD
        assert OPT_T_RUN_LENGTH_THRESHOLD == 6

    def test_iupac_expansion_cap(self):
        from biocompiler.constants import OPT_IUPAC_EXPANSION_CAP
        assert OPT_IUPAC_EXPANSION_CAP == 4096

    def test_splice_donor_threshold(self):
        from biocompiler.constants import OPT_SPLICE_DONOR_POTENTIAL_THRESHOLD
        assert OPT_SPLICE_DONOR_POTENTIAL_THRESHOLD == pytest.approx(0.5)

    def test_eukaryote_cai_gt_cost(self):
        from biocompiler.constants import OPT_EUKARYOTE_CAI_GT_COST_THRESHOLD
        assert OPT_EUKARYOTE_CAI_GT_COST_THRESHOLD == pytest.approx(0.10)


# ── 2. YAML config override ───────────────────────────────────────────────

class TestYAMLConfigOverride:
    """Test BIOCOMPILER_CONFIG env var override support."""

    def test_load_config_overrides_function_exists(self):
        from biocompiler.constants import load_config_overrides
        assert callable(load_config_overrides)

    def test_yaml_override_with_valid_file(self, tmp_path):
        """A YAML file can override numeric constants."""
        import importlib
        import biocompiler.constants as const

        yaml_content = "OPT_T_RUN_LENGTH_THRESHOLD: 8\n"
        yaml_file = tmp_path / "test_config.yaml"
        yaml_file.write_text(yaml_content)

        # Save original
        original = const.OPT_T_RUN_LENGTH_THRESHOLD

        # Set env and reload overrides
        os.environ["BIOCOMPILER_CONFIG"] = str(yaml_file)
        try:
            const.load_config_overrides()
            assert const.OPT_T_RUN_LENGTH_THRESHOLD == 8

            # Restore original
            const.OPT_T_RUN_LENGTH_THRESHOLD = original
        finally:
            del os.environ["BIOCOMPILER_CONFIG"]

    def test_yaml_override_missing_key_ignored(self, tmp_path):
        """Unknown keys in YAML are ignored with a debug log."""
        from biocompiler.constants import load_config_overrides
        import biocompiler.constants as const

        yaml_content = "NONEXISTENT_CONSTANT: 999\n"
        yaml_file = tmp_path / "test_config2.yaml"
        yaml_file.write_text(yaml_content)

        os.environ["BIOCOMPILER_CONFIG"] = str(yaml_file)
        try:
            # Should not raise
            load_config_overrides()
        finally:
            del os.environ["BIOCOMPILER_CONFIG"]

    def test_yaml_override_nonexistent_file(self, tmp_path):
        """A non-existent config file logs a warning and is ignored."""
        from biocompiler.constants import load_config_overrides

        os.environ["BIOCOMPILER_CONFIG"] = str(tmp_path / "nonexistent.yaml")
        try:
            # Should not raise
            load_config_overrides()
        finally:
            del os.environ["BIOCOMPILER_CONFIG"]


# ── 3. Validation fixes: assert → ValueError ──────────────────────────────

class TestOptimizationResultValidation:
    """Verify OptimizationResult raises ValueError instead of AssertionError."""

    def test_invalid_cai_raises_valueerror(self):
        from biocompiler.optimization import OptimizationResult
        with pytest.raises(ValueError, match="CAI must be in"):
            OptimizationResult(
                sequence="ATGAAATTT", protein="MKF",
                cai=1.5, gc_content=0.5,
            )

    def test_invalid_gc_raises_valueerror(self):
        from biocompiler.optimization import OptimizationResult
        with pytest.raises(ValueError, match="GC content must be in"):
            OptimizationResult(
                sequence="ATGAAATTT", protein="MKF",
                cai=0.5, gc_content=2.0,
            )

    def test_negative_cai_raises_valueerror(self):
        from biocompiler.optimization import OptimizationResult
        with pytest.raises(ValueError, match="CAI must be in"):
            OptimizationResult(
                sequence="ATGAAATTT", protein="MKF",
                cai=-0.1, gc_content=0.5,
            )

    def test_mutagenesis_without_substitutions_raises_valueerror(self):
        from biocompiler.optimization import OptimizationResult
        with pytest.raises(ValueError, match="Mutagenesis applied"):
            OptimizationResult(
                sequence="ATGAAATTT", protein="MKF",
                cai=0.5, gc_content=0.5,
                mutagenesis_applied=True,
            )

    def test_invalid_utr_score_raises_valueerror(self):
        from biocompiler.optimization import OptimizationResult
        with pytest.raises(ValueError, match="UTR 5' score"):
            OptimizationResult(
                sequence="ATGAAATTT", protein="MKF",
                cai=0.5, gc_content=0.5,
                utr_score_5=1.5,
            )

    def test_valid_result_no_exception(self):
        from biocompiler.optimization import OptimizationResult
        result = OptimizationResult(
            sequence="ATGAAATTT", protein="MKF",
            cai=0.5, gc_content=0.5,
        )
        assert result.cai == pytest.approx(0.5, rel=1e-6)


class TestFullConstructResultValidation:
    """Verify FullConstructResult raises ValueError instead of AssertionError."""

    def test_construct_mismatch_raises_valueerror(self):
        from biocompiler.optimization import FullConstructResult
        with pytest.raises(ValueError, match="full_construct must equal"):
            FullConstructResult(
                utr5="AAAA", cds="ATGAAATTT", utr3="CCCC",
                full_construct="WRONG", organism="ecoli",
                gc_content=0.5, cai=0.5,
            )

    def test_invalid_gc_raises_valueerror(self):
        from biocompiler.optimization import FullConstructResult
        with pytest.raises(ValueError, match="GC content must be in"):
            FullConstructResult(
                utr5="", cds="ATGAAATTT", utr3="",
                full_construct="ATGAAATTT", organism="ecoli",
                gc_content=2.0, cai=0.5,
            )


class TestSlidingGCValidation:
    """Verify sliding_gc raises ValueError instead of AssertionError."""

    def test_invalid_direction_raises_valueerror(self):
        from biocompiler.sliding_gc import WindowViolation
        with pytest.raises(ValueError, match="direction must be"):
            WindowViolation(start=0, end=10, gc_content=0.5, direction="invalid")

    def test_invalid_gc_content_raises_valueerror(self):
        from biocompiler.sliding_gc import WindowViolation
        with pytest.raises(ValueError, match="gc_content must be in"):
            WindowViolation(start=0, end=10, gc_content=2.0, direction="too_high")

    def test_negative_start_raises_valueerror(self):
        from biocompiler.sliding_gc import WindowViolation
        with pytest.raises(ValueError, match="start must be >= 0"):
            WindowViolation(start=-1, end=10, gc_content=0.5, direction="too_low")

    def test_end_before_start_raises_valueerror(self):
        from biocompiler.sliding_gc import WindowViolation
        with pytest.raises(ValueError, match="end must be > start"):
            WindowViolation(start=10, end=5, gc_content=0.5, direction="too_low")

    def test_check_sliding_gc_invalid_window_size(self):
        from biocompiler.sliding_gc import check_sliding_gc
        with pytest.raises(ValueError, match="window_size must be > 0"):
            check_sliding_gc("ATGCAA", window_size=0, gc_min=0.3, gc_max=0.7)

    def test_check_sliding_gc_invalid_gc_range(self):
        from biocompiler.sliding_gc import check_sliding_gc
        with pytest.raises(ValueError, match="gc_min must be in"):
            check_sliding_gc("ATGCAA", window_size=3, gc_min=0.8, gc_max=0.2)


class TestMultigeneValidation:
    """Verify multigene raises ValueError instead of AssertionError."""

    def test_invalid_gc_content(self):
        from biocompiler.multigene import MultiGeneResult
        with pytest.raises(ValueError, match="gc_content must be in"):
            MultiGeneResult(
                genes=[], full_dna="ATGAAATTT",
                total_length=9, genbank_export="",
                construct_type="operon", organism="ecoli",
                gc_content=2.0,
            )

    def test_invalid_construct_type(self):
        from biocompiler.multigene import MultiGeneResult
        with pytest.raises(ValueError, match="Invalid construct_type"):
            MultiGeneResult(
                genes=[], full_dna="ATGAAATTT",
                total_length=9, genbank_export="",
                construct_type="invalid_type", organism="ecoli",
                gc_content=0.5,
            )


class TestOptimizationFunctionValidation:
    """Verify optimization functions raise ValueError for bad inputs."""

    def test_invalid_gc_bounds_internal(self):
        """_greedy_optimize should raise ValueError for invalid GC bounds."""
        from biocompiler.optimization import _greedy_optimize
        with pytest.raises(ValueError, match="Invalid GC range"):
            _greedy_optimize(
                protein="MKF", organism="Escherichia_coli",
                gc_lo=0.8, gc_hi=0.2,
                cryptic_splice_threshold=3.0,
                is_prokaryote=True,
            )

    def test_negative_threshold_internal(self):
        """_greedy_optimize should raise ValueError for negative threshold."""
        from biocompiler.optimization import _greedy_optimize
        with pytest.raises(ValueError, match="Threshold must be positive"):
            _greedy_optimize(
                protein="MKF", organism="Escherichia_coli",
                gc_lo=0.3, gc_hi=0.7,
                cryptic_splice_threshold=-1.0,
                is_prokaryote=True,
            )


# ── 4. Species-to-organism deduplication ──────────────────────────────────

class TestSpeciesToOrganismDeduplication:
    """Verify BioOptimizer and HybridOptimizer use resolve_organism."""

    def test_biooptimizer_no_species_to_organism_dict(self):
        """BioOptimizer should no longer have _SPECIES_TO_ORGANISM."""
        from biocompiler.optimization import BioOptimizer
        assert not hasattr(BioOptimizer, '_SPECIES_TO_ORGANISM')

    def test_hybrid_optimizer_no_species_to_organism_dict(self):
        """HybridOptimizer should no longer have _SPECIES_TO_ORGANISM."""
        from biocompiler.hybrid_optimizer import HybridOptimizer
        assert not hasattr(HybridOptimizer, '_SPECIES_TO_ORGANISM')

    def test_biooptimizer_resolves_ecoli(self):
        """BioOptimizer resolves 'ecoli' to canonical 'Escherichia_coli'."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli")
        assert opt.organism_name == "Escherichia_coli"

    def test_biooptimizer_resolves_human(self):
        """BioOptimizer resolves 'human' to canonical 'Homo_sapiens'."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="human")
        assert opt.organism_name == "Homo_sapiens"

    def test_biooptimizer_resolves_alias(self):
        """BioOptimizer resolves extended aliases like 'E. coli'."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="E. coli")
        assert opt.organism_name == "Escherichia_coli"

    def test_biooptimizer_organism_name_kwarg(self):
        """BioOptimizer resolves organism_name kwarg."""
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(species="ecoli", organism_name="Homo_sapiens")
        assert opt.organism_name == "Homo_sapiens"

    def test_hybrid_optimizer_resolves_species(self):
        """HybridOptimizer resolves species via resolve_organism."""
        from biocompiler.hybrid_optimizer import HybridOptimizer
        # HybridOptimizer needs protein sequence; use a simple one
        opt = HybridOptimizer(species="ecoli")
        assert opt.organism == "Escherichia_coli"

    def test_resolve_organism_has_more_aliases(self):
        """resolve_organism should support more aliases than the old 10-entry dict."""
        from biocompiler.organisms import resolve_organism, ORGANISM_ALIASES
        # The old _SPECIES_TO_ORGANISM had 10 entries.
        # ORGANISM_ALIASES should have significantly more.
        assert len(ORGANISM_ALIASES) > 20
        # Test some aliases that the old dict didn't have
        assert resolve_organism("E. coli") == "Escherichia_coli"
        assert resolve_organism("h_sapiens") == "Homo_sapiens"
        assert resolve_organism("M. musculus") == "Mus_musculus"


# ── 5. Exception logging ─────────────────────────────────────────────────

class TestExceptionLogging:
    """Verify that broad exception catches now log warnings."""

    def test_numba_kernels_has_logger(self):
        """numba_kernels.py should have a logger for warmup failure."""
        from biocompiler import numba_kernels
        assert hasattr(numba_kernels, 'logger')

    def test_optimization_module_logger_exists(self):
        """optimization.py should have a logger for NUMBA fallbacks."""
        from biocompiler.optimization import logger
        assert logger is not None
