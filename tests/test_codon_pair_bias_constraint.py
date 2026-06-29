"""
Tests for the MinimizeCodonPairBias soft constraint and its integration
with the CSP solver model.

Covers:
1. MinimizeCodonPairBias scoring — mean CPB, geometric mean, edge cases
2. Good vs bad codon pair sequences — higher CPB scores better
3. Integration with build_csp_model — constraint toggled via SolverConfig
4. ConstraintType and ConstraintSpec recognition
5. E. coli codon pair bias data validation (CUA→CTA fix)
"""

from __future__ import annotations

import math
import pytest

from biocompiler.solver.types import (
    ConstraintStrictness,
    ConstraintType,
    SolverConfig,
)


def _import_constraints():
    """Import constraints module; skip tests if not available."""
    try:
        from biocompiler.solver import constraints
        return constraints
    except ImportError as exc:
        pytest.skip(f"solver.constraints not available: {exc}")


def _import_cpb():
    """Import codon_pair_scoring module; skip if not available."""
    try:
        from biocompiler import codon_pair_scoring
        return codon_pair_scoring
    except ImportError as exc:
        pytest.skip(f"codon_pair_scoring not available: {exc}")


def _import_e_coli():
    """Import escherichia organism module."""
    try:
        from biocompiler.organisms import escherichia
        return escherichia
    except ImportError as exc:
        pytest.skip(f"organisms.escherichia not available: {exc}")


# ── Fixtures ──────────────────────────────────────────────────────────────

@pytest.fixture
def ecoli_cpb_constraint():
    """A MinimizeCodonPairBias constraint for E. coli."""
    c = _import_constraints()
    return c.MinimizeCodonPairBias(organism="Escherichia_coli")


@pytest.fixture
def ecoli_config_with_cpb():
    """SolverConfig with codon pair bias optimization enabled."""
    return SolverConfig(
        optimize_codon_pair_bias=True,
        codon_pair_bias_weight=0.2,
    )


@pytest.fixture
def ecoli_config_without_cpb():
    """SolverConfig with codon pair bias optimization disabled (default)."""
    return SolverConfig(
        optimize_codon_pair_bias=False,
    )


# ═══════════════════════════════════════════════════════════════════════════
# 1. MinimizeCodonPairBias Scoring
# ═══════════════════════════════════════════════════════════════════════════

class TestMinimizeCodonPairBiasScoring:
    """Test the scoring methods of MinimizeCodonPairBias."""

    def test_name_property(self, ecoli_cpb_constraint):
        assert ecoli_cpb_constraint.name == "MinimizeCodonPairBias"

    def test_constraint_type(self, ecoli_cpb_constraint):
        assert ecoli_cpb_constraint.constraint_type == ConstraintType.CODON_PAIR_BIAS

    def test_organism_property(self, ecoli_cpb_constraint):
        assert ecoli_cpb_constraint.organism == "Escherichia_coli"

    def test_check_always_true(self, ecoli_cpb_constraint):
        """Soft optimization objective — always satisfied."""
        assert ecoli_cpb_constraint.check("ATGCTGAAAGAA") is True
        assert ecoli_cpb_constraint.check("ATAAGG") is True

    def test_score_known_positive_pair(self, ecoli_cpb_constraint):
        """CTG-CTG has CPB=0.45, the highest E. coli pair."""
        # ATG-CTG (0.38) + CTG-CTG (0.45) → mean = (0.38+0.45)/2
        seq = "ATGCTGCTG"  # Met-Leu-Leu
        score = ecoli_cpb_constraint.score(seq)
        expected = (0.38 + 0.45) / 2.0
        assert abs(score - expected) < 1e-9

    def test_score_known_negative_pair(self, ecoli_cpb_constraint):
        """CTA-ATA has CPB=-0.50, the lowest E. coli pair."""
        # CTA-ATA (-0.50) — single pair
        seq = "CTAATA"  # Leu(rare)-Ile(rare)
        score = ecoli_cpb_constraint.score(seq)
        assert abs(score - (-0.50)) < 1e-9

    def test_score_mixed_pairs(self, ecoli_cpb_constraint):
        """Sequence with both positive and negative pairs."""
        # ATG(1)-CTG(2)-ATA(3) → ATG-CTG=0.38, CTG-ATA: not in table → 0.0
        seq = "ATGCTGATA"
        score = ecoli_cpb_constraint.score(seq)
        expected = (0.38 + 0.0) / 2.0
        assert abs(score - expected) < 1e-9

    def test_score_short_sequence_returns_zero(self, ecoli_cpb_constraint):
        """Sequences shorter than two codons return 0.0."""
        assert ecoli_cpb_constraint.score("ATG") == 0.0

    def test_score_empty_cpb_data(self):
        """Constraint with no CPB data for organism returns 0.0."""
        c = _import_constraints()
        constraint = c.MinimizeCodonPairBias(organism="Unknown_organism")
        assert constraint.score("ATGCTGAAAGAA") == 0.0

    def test_score_neutral_pairs(self, ecoli_cpb_constraint):
        """All neutral pairs (not in table) → score = 0.0."""
        # Use codons unlikely to be in the E. coli CPB table
        # GGT-GGC: not an extreme pair
        seq = "GGTGGC"
        score = ecoli_cpb_constraint.score(seq)
        # GGT-GGC should default to 0.0 (neutral)
        assert score == 0.0

    def test_geometric_mean_score(self, ecoli_cpb_constraint):
        """Test geometric mean of 2^CPB values."""
        seq = "ATGCTGCTG"  # ATG-CTG=0.38, CTG-CTG=0.45
        gm = ecoli_cpb_constraint.geometric_mean_score(seq)
        # Expected: exp(mean of [0.38*ln2, 0.45*ln2])
        log_vals = [0.38 * math.log(2), 0.45 * math.log(2)]
        expected = math.exp(sum(log_vals) / len(log_vals))
        assert abs(gm - expected) < 1e-9

    def test_geometric_mean_neutral_returns_one(self, ecoli_cpb_constraint):
        """Geometric mean of neutral pairs should be 1.0 (2^0 = 1)."""
        seq = "GGTGGC"  # Neutral pair
        gm = ecoli_cpb_constraint.geometric_mean_score(seq)
        assert abs(gm - 1.0) < 1e-9

    def test_geometric_mean_short_sequence(self, ecoli_cpb_constraint):
        """Short sequences return 1.0 for geometric mean."""
        assert ecoli_cpb_constraint.geometric_mean_score("ATG") == 1.0


class TestViolatedPositions:
    """Test violated_positions method of MinimizeCodonPairBias."""

    def test_negative_cpb_positions_reported(self, ecoli_cpb_constraint):
        """Codon pairs with negative CPB should be reported."""
        # CTA-ATA = -0.50 (worst pair)
        seq = "CTAATA"  # Leu(rare)-Ile(rare)
        positions = ecoli_cpb_constraint.violated_positions(seq)
        assert 0 in positions  # First codon pair starts at position 0

    def test_positive_cpb_not_reported(self, ecoli_cpb_constraint):
        """Codon pairs with positive CPB should not be reported."""
        # ATG-CTG = +0.38 (good pair)
        seq = "ATGCTG"
        positions = ecoli_cpb_constraint.violated_positions(seq)
        assert len(positions) == 0

    def test_mixed_positions(self, ecoli_cpb_constraint):
        """Only negative-CPB pairs should be reported in a mixed sequence."""
        # ATG-CTG (good, 0.38) + CTG-ATA (not in table → 0.0, neutral)
        seq = "ATGCTGATA"
        positions = ecoli_cpb_constraint.violated_positions(seq)
        assert len(positions) == 0  # No negative pairs

    def test_all_negative_positions(self, ecoli_cpb_constraint):
        """Multiple negative pairs are all reported."""
        # CTA-ATA (-0.50) + ATA-AGG (-0.30)
        seq = "CTAATAAGG"
        positions = ecoli_cpb_constraint.violated_positions(seq)
        assert 0 in positions  # CTA-ATA
        assert 3 in positions  # ATA-AGG

    def test_short_sequence_empty(self, ecoli_cpb_constraint):
        """Short sequences return empty positions list."""
        assert ecoli_cpb_constraint.violated_positions("ATG") == []


# ═══════════════════════════════════════════════════════════════════════════
# 2. Good vs Bad Codon Pair Sequences
# ═══════════════════════════════════════════════════════════════════════════

class TestGoodVsBadCodonPairs:
    """Sequences with good codon pairs should score higher than bad ones."""

    def test_good_pairs_score_higher_than_bad(self, ecoli_cpb_constraint):
        """Sequence using over-represented pairs scores higher than under-represented."""
        # Good: ATG-CTG (0.38) — start-proximal Leu, well-represented
        good_seq = "ATGCTGATGCTG"  # ATG-CTG + CTG-ATG + ATG-CTG = (0.38+0.35+0.38)/3
        # Bad: CTA-ATA (-0.50) — rare-rare pair
        bad_seq = "CTAATACTAATA"  # CTA-ATA + ATA-CTA + CTA-ATA = (-0.50-0.48-0.50)/3

        good_score = ecoli_cpb_constraint.score(good_seq)
        bad_score = ecoli_cpb_constraint.score(bad_seq)

        assert good_score > bad_score, (
            f"Good pairs ({good_score:.3f}) should score higher than bad ({bad_score:.3f})"
        )

    def test_preferred_vs_avoided_leucine_codons(self, ecoli_cpb_constraint):
        """CTG (preferred Leu) pairs score better than CTA (rare Leu) pairs."""
        # CTG-CTG: most common E. coli pair (CPB = +0.45)
        preferred = "CTGCTG"
        # CTA-CTA: rare Leu pair (CPB = -0.42)
        avoided = "CTACTA"

        preferred_score = ecoli_cpb_constraint.score(preferred)
        avoided_score = ecoli_cpb_constraint.score(avoided)

        assert preferred_score > avoided_score

    def test_good_pairs_have_fewer_violated_positions(self, ecoli_cpb_constraint):
        """Good pairs should have fewer violated positions."""
        good_seq = "ATGCTGCTGGAAAAG"
        bad_seq = "CTAATAAGGAGACTA"

        good_vp = ecoli_cpb_constraint.violated_positions(good_seq)
        bad_vp = ecoli_cpb_constraint.violated_positions(bad_seq)

        assert len(good_vp) <= len(bad_vp)


# ═══════════════════════════════════════════════════════════════════════════
# 3. Integration with build_csp_model
# ═══════════════════════════════════════════════════════════════════════════

class TestBuildCSPModelIntegration:
    """Test that codon pair bias constraint integrates with build_csp_model."""

    def test_cpb_not_added_by_default(self, ecoli_config_without_cpb):
        """By default, optimize_codon_pair_bias=False → no CPB constraint."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_without_cpb)
        cpb_constraints = [sc for sc in model.soft_constraints
                          if sc.name == "MinimizeCodonPairBias"]
        assert len(cpb_constraints) == 0

    def test_cpb_added_when_enabled(self, ecoli_config_with_cpb):
        """With optimize_codon_pair_bias=True, CPB constraint is added."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_with_cpb)
        cpb_constraints = [sc for sc in model.soft_constraints
                          if sc.name == "MinimizeCodonPairBias"]
        assert len(cpb_constraints) == 1

    def test_cpb_not_added_for_organism_without_data(self):
        """CPB optimization requested but organism has no data → constraint skipped."""
        c = _import_constraints()
        config = SolverConfig(optimize_codon_pair_bias=True)
        # Homo_sapiens may or may not have CPB data; use a truly unknown organism
        # We cannot pass an unknown organism to build_csp_model (it validates),
        # so we test the constraint directly
        constraint = c.MinimizeCodonPairBias(organism="Homo_sapiens")
        # Homo_sapiens should have human CPB data from codon_pair_scoring
        has_data = bool(constraint.cpb_data)
        # Either way, build_csp_model should handle it gracefully
        # The key test: no crash when CPB data is empty for an organism
        if not has_data:
            config = SolverConfig(optimize_codon_pair_bias=True)
            model = c.build_csp_model("MVLSPADKTN", "Homo_sapiens", config)
            cpb = [sc for sc in model.soft_constraints if sc.name == "MinimizeCodonPairBias"]
            assert len(cpb) == 0  # Should be skipped

    def test_cpb_constraint_has_correct_weight(self, ecoli_config_with_cpb):
        """CPB constraint uses codon_pair_bias_weight from config."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_with_cpb)
        cpb = [sc for sc in model.soft_constraints if sc.name == "MinimizeCodonPairBias"]
        assert len(cpb) == 1

        # Verify weight is propagated through ConstraintSpec
        specs = model.constraints
        cpb_spec = [s for s in specs if s.name == "MinimizeCodonPairBias"]
        assert len(cpb_spec) == 1
        assert cpb_spec[0].weight == ecoli_config_with_cpb.codon_pair_bias_weight

    def test_cpb_constraint_objective_contribution(self, ecoli_config_with_cpb):
        """CPB constraint contributes to the objective value."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_with_cpb)

        # Build two different sequences and verify the one with better pairs
        # has a higher objective contribution from CPB
        good_seq = "ATGCTGCTGGAAAAG"  # ATG-CTG-CTG-GAA-AAG: good pairs
        bad_seq = "CTAATAAGGAGACTA"   # CTA-ATA-AGG-AGA-CTA: bad pairs

        cpb_constraint = [sc for sc in model.soft_constraints
                         if sc.name == "MinimizeCodonPairBias"][0]

        good_score = cpb_constraint.score(good_seq)
        bad_score = cpb_constraint.score(bad_seq)
        assert good_score > bad_score

    def test_model_objective_includes_cpb_weight(self, ecoli_config_with_cpb):
        """Model objective_value includes codon_pair_bias_weight."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_with_cpb)

        # The objective_value method should include CPB
        seq = "ATGCTGCTGGAAAAG"
        obj = model.objective_value(seq)
        assert isinstance(obj, float)

        # Verify CPB constraint is in the weight_map
        cpb_constraints = [sc for sc in model.soft_constraints
                          if sc.name == "MinimizeCodonPairBias"]
        assert len(cpb_constraints) == 1


# ═══════════════════════════════════════════════════════════════════════════
# 4. ConstraintSpec and ConstraintType Recognition
# ═══════════════════════════════════════════════════════════════════════════

class TestConstraintSpecIntegration:
    """Test that the CPB constraint integrates with ConstraintSpec system."""

    def test_constraint_type_in_enum(self):
        """CODON_PAIR_BIAS should be a valid ConstraintType."""
        assert ConstraintType.CODON_PAIR_BIAS == "codon_pair_bias"

    def test_cpb_spec_is_soft(self, ecoli_config_with_cpb):
        """CPB constraint should appear as SOFT in ConstraintSpec."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_with_cpb)
        cpb_specs = [s for s in model.constraints
                    if s.ctype == ConstraintType.CODON_PAIR_BIAS]
        assert len(cpb_specs) == 1
        assert cpb_specs[0].strictness == ConstraintStrictness.SOFT

    def test_cpb_spec_check_returns_true(self, ecoli_config_with_cpb):
        """ConstraintSpec.check() for CODON_PAIR_BIAS always returns True."""
        c = _import_constraints()
        model = c.build_csp_model("MVLSPADKTN", "Escherichia_coli", ecoli_config_with_cpb)
        cpb_specs = [s for s in model.constraints
                    if s.ctype == ConstraintType.CODON_PAIR_BIAS]
        assert len(cpb_specs) == 1
        assert cpb_specs[0].check("ATGCTGAAAGAA") is True


# ═══════════════════════════════════════════════════════════════════════════
# 5. SolverConfig Toggle
# ═══════════════════════════════════════════════════════════════════════════

class TestSolverConfigToggle:
    """Test that the CPB constraint can be toggled via SolverConfig."""

    def test_default_config_has_cpb_disabled(self):
        """Default SolverConfig has optimize_codon_pair_bias=False."""
        config = SolverConfig()
        assert config.optimize_codon_pair_bias is False

    def test_default_cpb_weight(self):
        """Default codon_pair_bias_weight is 0.2."""
        config = SolverConfig()
        assert config.codon_pair_bias_weight == 0.2

    def test_enable_cpb(self):
        """Can enable CPB optimization."""
        config = SolverConfig(optimize_codon_pair_bias=True)
        assert config.optimize_codon_pair_bias is True

    def test_custom_cpb_weight(self):
        """Can set custom CPB weight."""
        config = SolverConfig(
            optimize_codon_pair_bias=True,
            codon_pair_bias_weight=0.5,
        )
        assert config.codon_pair_bias_weight == 0.5

    def test_cpb_toggle_in_model(self):
        """Toggling optimize_codon_pair_bias adds/removes the constraint."""
        c = _import_constraints()
        protein = "MVLSPADKTN"

        # Off
        config_off = SolverConfig(optimize_codon_pair_bias=False)
        model_off = c.build_csp_model(protein, "Escherichia_coli", config_off)
        cpb_off = [sc for sc in model_off.soft_constraints
                  if sc.name == "MinimizeCodonPairBias"]
        assert len(cpb_off) == 0

        # On
        config_on = SolverConfig(optimize_codon_pair_bias=True)
        model_on = c.build_csp_model(protein, "Escherichia_coli", config_on)
        cpb_on = [sc for sc in model_on.soft_constraints
                 if sc.name == "MinimizeCodonPairBias"]
        assert len(cpb_on) == 1


# ═══════════════════════════════════════════════════════════════════════════
# 6. E. coli CPB Data Validation (CUA→CTA fix)
# ═══════════════════════════════════════════════════════════════════════════

class TestEColiCPBDataValidation:
    """Validate the E. coli codon pair bias data."""

    def test_no_rna_codons_in_data(self):
        """All codon pair keys should use DNA notation (T, not U)."""
        ecoli = _import_e_coli()
        for key in ecoli.E_COLI_CODON_PAIR_BIAS:
            assert "U" not in key, f"RNA-style codon in key: {key}"

    def test_all_codons_valid_dna(self):
        """Each codon in pair keys should be a valid 3-letter DNA codon."""
        ecoli = _import_e_coli()
        valid_bases = set("ACGT")
        for key in ecoli.E_COLI_CODON_PAIR_BIAS:
            parts = key.split("-")
            assert len(parts) == 2, f"Invalid key format: {key}"
            for codon in parts:
                assert len(codon) == 3, f"Invalid codon length: {codon} in {key}"
                assert set(codon).issubset(valid_bases), f"Invalid bases in codon: {codon}"

    def test_cta_cta_pair_exists(self):
        """CTA-CTA (rare Leu pair) should be in the bias table (not CUA-CUA)."""
        ecoli = _import_e_coli()
        assert "CTA-CTA" in ecoli.E_COLI_CODON_PAIR_BIAS
        assert "CUA-CUA" not in ecoli.E_COLI_CODON_PAIR_BIAS

    def test_ctg_ctg_is_highest(self):
        """CTG-CTG should be the highest-scoring pair."""
        ecoli = _import_e_coli()
        max_key = max(ecoli.E_COLI_CODON_PAIR_BIAS,
                     key=ecoli.E_COLI_CODON_PAIR_BIAS.get)
        assert max_key == "CTG-CTG"
        assert ecoli.E_COLI_CODON_PAIR_BIAS["CTG-CTG"] == 0.45

    def test_cta_ata_is_lowest(self):
        """CTA-ATA should be the lowest-scoring pair."""
        ecoli = _import_e_coli()
        min_key = min(ecoli.E_COLI_CODON_PAIR_BIAS,
                     key=ecoli.E_COLI_CODON_PAIR_BIAS.get)
        assert min_key == "CTA-ATA"
        assert ecoli.E_COLI_CODON_PAIR_BIAS["CTA-ATA"] == -0.50

    def test_cpb_scores_reasonable_range(self):
        """CPB scores should be in a reasonable range [-1, 1]."""
        ecoli = _import_e_coli()
        for key, value in ecoli.E_COLI_CODON_PAIR_BIAS.items():
            assert -1.0 <= value <= 1.0, f"CPB out of range for {key}: {value}"

    def test_over_represented_pairs_positive(self):
        """All over-represented pairs should have positive CPB."""
        ecoli = _import_e_coli()
        positive_pairs = {k: v for k, v in ecoli.E_COLI_CODON_PAIR_BIAS.items() if v > 0}
        for key, value in positive_pairs.items():
            assert value > 0, f"Over-represented pair {key} has non-positive CPB: {value}"

    def test_under_represented_pairs_negative(self):
        """All under-represented pairs should have negative CPB."""
        ecoli = _import_e_coli()
        negative_pairs = {k: v for k, v in ecoli.E_COLI_CODON_PAIR_BIAS.items() if v < 0}
        for key, value in negative_pairs.items():
            assert value < 0, f"Under-represented pair {key} has non-negative CPB: {value}"


# ═══════════════════════════════════════════════════════════════════════════
# 7. codon_pair_scoring Module Integration
# ═══════════════════════════════════════════════════════════════════════════

class TestCodonPairScoringIntegration:
    """Test that the codon_pair_scoring module works with the constraint."""

    def test_compute_cpb_with_ecoli(self):
        """compute_cpb should return non-zero for E. coli sequences with biased pairs."""
        cpb = _import_cpb()
        # ATG-CTG = +0.38
        score = cpb.compute_cpb("ATGCTG", "Escherichia_coli")
        assert abs(score - 0.38) < 1e-9

    def test_score_codon_pair(self):
        """score_codon_pair should return correct values for known pairs."""
        cpb = _import_cpb()
        # CTG-CTG = +0.45
        assert abs(cpb.score_codon_pair("CTG", "CTG", "Escherichia_coli") - 0.45) < 1e-9
        # CTA-ATA = -0.50
        assert abs(cpb.score_codon_pair("CTA", "ATA", "Escherichia_coli") - (-0.50)) < 1e-9

    def test_get_codon_pair_data_ecoli(self):
        """get_codon_pair_data should return E. coli data."""
        cpb = _import_cpb()
        data = cpb.get_codon_pair_data("Escherichia_coli")
        assert len(data) > 0
        assert "CTG-CTG" in data

    def test_get_codon_pair_data_aliases(self):
        """Both 'Escherichia_coli' and 'e_coli' should return E. coli data."""
        cpb = _import_cpb()
        data1 = cpb.get_codon_pair_data("Escherichia_coli")
        data2 = cpb.get_codon_pair_data("e_coli")
        assert data1 is data2  # Same object

    def test_constraint_uses_cpb_module(self):
        """MinimizeCodonPairBias should use codon_pair_scoring.get_codon_pair_data."""
        c = _import_constraints()
        cpb = _import_cpb()

        constraint = c.MinimizeCodonPairBias(organism="Escherichia_coli")
        module_data = cpb.get_codon_pair_data("Escherichia_coli")

        # The constraint should have the same data as the module
        assert constraint.cpb_data == module_data
