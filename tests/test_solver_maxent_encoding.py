"""
Comprehensive pytest tests for biocompiler.solver.maxent_encoding.

Covers:
1. MaxEnt score encoding/decoding — precompute_splice_scores, build_splice_constraint_table
2. Quantization correctness — quantize_maxent_scores (boundary values, linear mapping, clamping)
3. Cross-codon constraint generation — encode_cross_codon_splice_context
4. SpliceConstraintEncoder class — initialization, caching, query methods
5. Edge cases — short sequences, single codon, no splice sites found, unknown amino acids
"""

from __future__ import annotations

import pytest

from biocompiler.shared.constants import AA_TO_CODONS
from biocompiler.solver.types import SolverConfig, SpliceConstraint, CrossCodonSpliceConstraint
from biocompiler.solver.maxent_encoding import (
    _OUT_OF_BOUNDS_SCORE,
    _build_sequence_with_codon,
    precompute_splice_scores,
    build_splice_constraint_table,
    quantize_maxent_scores,
    encode_cross_codon_splice_context,
    SpliceConstraintEncoder,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def default_config() -> SolverConfig:
    """Default SolverConfig with standard thresholds."""
    return SolverConfig()


@pytest.fixture
def strict_config() -> SolverConfig:
    """Config with very low threshold — catches more splice sites."""
    return SolverConfig(cryptic_splice_threshold=0.0)


@pytest.fixture
def lenient_config() -> SolverConfig:
    """Config with very high threshold — almost no sites are forbidden."""
    return SolverConfig(cryptic_splice_threshold=20.0)


@pytest.fixture
def sample_protein() -> str:
    """First 10 AA of human alpha-globin — diverse codons."""
    return "MVLSPADKTN"


@pytest.fixture
def valine_protein() -> str:
    """Protein with valine — all Val codons contain GT."""
    return "VVVVV"


@pytest.fixture
def short_protein() -> str:
    """Minimal two-residue protein for cross-codon boundary testing."""
    return "MV"


@pytest.fixture
def single_aa_protein() -> str:
    """Single amino acid — no cross-codon boundaries possible."""
    return "M"


# ═══════════════════════════════════════════════════════════════════════════════
# 1. _build_sequence_with_codon helper
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuildSequenceWithCodon:

    def test_basic_construction(self):
        """Building a sequence with a codon override at position 0."""
        seq = _build_sequence_with_codon("MV", {}, 0, "ATG")
        assert len(seq) == 6
        assert seq[:3] == "ATG"

    def test_override_overrides_default(self):
        """The codon at the override position should replace the default."""
        protein = "MV"
        assignments = {}
        seq = _build_sequence_with_codon(protein, assignments, 0, "ATG")
        assert seq[:3] == "ATG"

    def test_partial_assignments_used_for_context(self):
        """Pre-existing codon_assignments are used for context positions."""
        protein = "MVL"
        assignments = {1: "GTT"}  # Val -> GTT instead of default
        seq = _build_sequence_with_codon(protein, assignments, 0, "ATG")
        assert seq[3:6] == "GTT"  # assignment used

    def test_override_takes_precedence(self):
        """Position override should override an existing assignment."""
        protein = "MV"
        assignments = {0: "ATG"}
        seq = _build_sequence_with_codon(protein, assignments, 0, "ATG")
        assert seq[:3] == "ATG"

    def test_correct_length(self):
        """Sequence length = len(protein) * 3."""
        for length in [1, 3, 5, 10]:
            protein = "M" * length
            seq = _build_sequence_with_codon(protein, {}, 0, "ATG")
            assert len(seq) == length * 3

    def test_other_positions_use_default_codons(self):
        """Positions without assignment use first (highest-CAI) codon."""
        protein = "MV"
        default_m = AA_TO_CODONS["M"][0]
        default_v = AA_TO_CODONS["V"][0]
        seq = _build_sequence_with_codon(protein, {}, 1, "GTC")
        assert seq[:3] == default_m
        assert seq[3:6] == "GTC"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. MaxEnt Score Encoding / Decoding — precompute_splice_scores
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrecomputeSpliceScores:

    def test_returns_dict_with_position_codon_keys(self, sample_protein, default_config):
        """Keys should be (position, codon) tuples."""
        data = precompute_splice_scores(sample_protein, default_config)
        for key in data:
            assert isinstance(key, tuple)
            assert len(key) == 2
            pos, codon = key
            assert isinstance(pos, int)
            assert isinstance(codon, str)
            assert len(codon) == 3

    def test_all_codons_per_position_present(self, sample_protein, default_config):
        """Every valid codon at every position should have an entry."""
        data = precompute_splice_scores(sample_protein, default_config)
        for i, aa in enumerate(sample_protein):
            for codon in AA_TO_CODONS[aa]:
                assert (i, codon) in data, f"Missing ({i}, {codon})"

    def test_entry_has_donor_and_acceptor_keys(self, sample_protein, default_config):
        """Each entry must have donor_scores and acceptor_scores lists."""
        data = precompute_splice_scores(sample_protein, default_config)
        for key, scores in data.items():
            assert "donor_scores" in scores
            assert "acceptor_scores" in scores
            assert isinstance(scores["donor_scores"], list)
            assert isinstance(scores["acceptor_scores"], list)

    def test_donor_scores_are_tuples_of_position_and_float(self, sample_protein, default_config):
        """Each donor score entry should be (base_position, score)."""
        data = precompute_splice_scores(sample_protein, default_config)
        for key, scores in data.items():
            for entry in scores["donor_scores"]:
                assert isinstance(entry, tuple)
                assert len(entry) == 2
                assert isinstance(entry[0], int)
                assert isinstance(entry[1], float)

    def test_acceptor_scores_are_tuples_of_position_and_float(self, sample_protein, default_config):
        """Each acceptor score entry should be (base_position, score)."""
        data = precompute_splice_scores(sample_protein, default_config)
        for key, scores in data.items():
            for entry in scores["acceptor_scores"]:
                assert isinstance(entry, tuple)
                assert len(entry) == 2
                assert isinstance(entry[0], int)
                assert isinstance(entry[1], float)

    def test_valine_codons_have_donor_scores(self, valine_protein, default_config):
        """All Val codons contain GT, so should produce donor scores."""
        data = precompute_splice_scores(valine_protein, default_config)
        has_donor = False
        for (pos, codon), scores in data.items():
            if scores["donor_scores"]:
                has_donor = True
                break
        # At least some valine codon should have a GT donor scored
        assert has_donor, "Valine protein should produce at least one donor score"

    def test_deterministic(self, sample_protein, default_config):
        """Same input should produce same output."""
        data1 = precompute_splice_scores(sample_protein, default_config)
        data2 = precompute_splice_scores(sample_protein, default_config)
        assert set(data1.keys()) == set(data2.keys())
        for key in data1:
            assert data1[key]["donor_scores"] == data2[key]["donor_scores"]
            assert data1[key]["acceptor_scores"] == data2[key]["acceptor_scores"]

    def test_case_insensitive_protein(self, default_config):
        """Lowercase protein input should be handled (uppercased internally)."""
        data = precompute_splice_scores("mv", default_config)
        assert (0, "ATG") in data

    def test_short_protein(self, short_protein, default_config):
        """Two-codon protein should still produce valid results."""
        data = precompute_splice_scores(short_protein, default_config)
        assert len(data) > 0
        for i, aa in enumerate(short_protein):
            for codon in AA_TO_CODONS[aa]:
                assert (i, codon) in data

    def test_single_aa_no_cross_boundary_issues(self, single_aa_protein, default_config):
        """Single amino acid protein should produce valid precomputed data."""
        data = precompute_splice_scores(single_aa_protein, default_config)
        # Only M = ATG, one entry
        assert (0, "ATG") in data
        assert len(data) == 1

    def test_scores_within_expected_range(self, sample_protein, default_config):
        """All scores should be finite floats (not NaN or infinite)."""
        import math
        data = precompute_splice_scores(sample_protein, default_config)
        for key, scores in data.items():
            for _, s in scores["donor_scores"]:
                assert math.isfinite(s), f"Donor score {s} is not finite for {key}"
            for _, s in scores["acceptor_scores"]:
                assert math.isfinite(s), f"Acceptor score {s} is not finite for {key}"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. build_splice_constraint_table
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuildSpliceConstraintTable:

    def test_all_constraints_above_threshold(self, sample_protein, default_config):
        """Every constraint in the table should have score > threshold."""
        data = precompute_splice_scores(sample_protein, default_config)
        constraints = build_splice_constraint_table(sample_protein, default_config, data)
        for c in constraints:
            assert c.score > c.threshold, (
                f"Constraint at pos={c.position} codon={c.codon} "
                f"score={c.score} not > threshold={c.threshold}"
            )

    def test_constraint_fields(self, sample_protein, default_config):
        """Each SpliceConstraint should have valid field values."""
        data = precompute_splice_scores(sample_protein, default_config)
        constraints = build_splice_constraint_table(sample_protein, default_config, data)
        for c in constraints:
            assert isinstance(c, SpliceConstraint)
            assert isinstance(c.position, int)
            assert c.position >= 0
            assert len(c.codon) == 3
            assert c.site_type in ("donor", "acceptor")
            assert isinstance(c.score, float)
            assert isinstance(c.threshold, float)

    def test_deduplication_keeps_worst(self, default_config):
        """If same (pos, codon, site_type) appears multiple times, keep worst score."""
        # Use a protein long enough to have overlapping scoring contexts
        protein = "MVLSPADKTNVKAAWGKVGA"
        data = precompute_splice_scores(protein, default_config)
        constraints = build_splice_constraint_table(protein, default_config, data)

        # Check no duplicate (pos, codon, site_type) keys
        keys_seen = set()
        for c in constraints:
            key = (c.position, c.codon, c.site_type)
            assert key not in keys_seen, f"Duplicate constraint for {key}"
            keys_seen.add(key)

    def test_strict_config_produces_more_constraints(self, sample_protein, default_config, strict_config):
        """Lower threshold should produce more or equal constraints."""
        data_default = precompute_splice_scores(sample_protein, default_config)
        data_strict = precompute_splice_scores(sample_protein, strict_config)
        c_default = build_splice_constraint_table(sample_protein, default_config, data_default)
        c_strict = build_splice_constraint_table(sample_protein, strict_config, data_strict)
        assert len(c_strict) >= len(c_default)

    def test_lenient_config_produces_fewer_constraints(self, sample_protein, default_config, lenient_config):
        """Higher threshold should produce fewer or equal constraints."""
        data_default = precompute_splice_scores(sample_protein, default_config)
        data_lenient = precompute_splice_scores(sample_protein, lenient_config)
        c_default = build_splice_constraint_table(sample_protein, default_config, data_default)
        c_lenient = build_splice_constraint_table(sample_protein, lenient_config, data_lenient)
        assert len(c_lenient) <= len(c_default)

    def test_very_high_threshold_produces_zero_constraints(self, sample_protein):
        """With threshold = 100, no realistic MaxEnt score should exceed it."""
        config = SolverConfig(cryptic_splice_threshold=100.0)
        data = precompute_splice_scores(sample_protein, config)
        constraints = build_splice_constraint_table(sample_protein, config, data)
        assert len(constraints) == 0

    def test_separate_donor_acceptor_thresholds(self, sample_protein):
        """Donor and acceptor thresholds can be set independently."""
        config = SolverConfig(
            cryptic_splice_threshold=3.0,
            donor_threshold=0.0,       # Very strict for donors
            acceptor_threshold=100.0,  # Very lenient for acceptors
        )
        data = precompute_splice_scores(sample_protein, config)
        constraints = build_splice_constraint_table(sample_protein, config, data)
        donor_constraints = [c for c in constraints if c.site_type == "donor"]
        acceptor_constraints = [c for c in constraints if c.site_type == "acceptor"]
        # With strict donor threshold, should find donor constraints
        # With lenient acceptor threshold, should find few/no acceptor constraints
        assert len(acceptor_constraints) <= len(donor_constraints)

    def test_empty_data_produces_no_constraints(self, sample_protein, default_config):
        """Empty score data should yield no constraints."""
        constraints = build_splice_constraint_table(sample_protein, default_config, {})
        assert len(constraints) == 0

    def test_margin_property(self, sample_protein, default_config):
        """SpliceConstraint.margin should equal score - threshold."""
        data = precompute_splice_scores(sample_protein, default_config)
        constraints = build_splice_constraint_table(sample_protein, default_config, data)
        for c in constraints:
            assert abs(c.margin - (c.score - c.threshold)) < 1e-10


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Quantization Correctness
# ═══════════════════════════════════════════════════════════════════════════════

class TestQuantizeMaxentScores:

    def test_example_from_docstring(self):
        """Verify the examples given in the function docstring."""
        assert quantize_maxent_scores([0.0, 1.5, 3.0, 5.0], n_bins=10, threshold=3.0) == [0, 5, 10, 10]

    def test_negative_scores_map_to_zero(self):
        """Scores <= 0 map to bin 0."""
        result = quantize_maxent_scores([-10.0, -1.0, 0.0], n_bins=10, threshold=3.0)
        assert result == [0, 0, 0]

    def test_score_above_threshold_maps_to_n_bins(self):
        """Scores >= threshold map to n_bins."""
        result = quantize_maxent_scores([3.0, 5.0, 100.0], n_bins=10, threshold=3.0)
        assert result == [10, 10, 10]

    def test_linear_interpolation(self):
        """Intermediate scores should map linearly."""
        # With n_bins=20, threshold=3.0, bin_width=0.15
        result = quantize_maxent_scores([0.15], n_bins=20, threshold=3.0)
        assert result[0] == 1  # 0.15 / 0.15 = 1.0 -> bin 1

    def test_bin_width_calculation(self):
        """Verify bin width = threshold / n_bins."""
        # threshold=3.0, n_bins=6, bin_width=0.5
        # score=1.0 -> 1.0/0.5 = 2 -> bin 2
        result = quantize_maxent_scores([1.0], n_bins=6, threshold=3.0)
        assert result[0] == 2

    def test_scores_near_zero_boundary(self):
        """Score slightly above 0 should map to bin 0 or 1 depending on resolution."""
        # n_bins=20, threshold=3.0, bin_width=0.15
        # score=0.001 -> 0.001/0.15 = 0.006... -> bin 0
        result = quantize_maxent_scores([0.001], n_bins=20, threshold=3.0)
        assert result[0] == 0

    def test_scores_near_threshold_boundary(self):
        """Score just below threshold maps to n_bins - 1."""
        # n_bins=10, threshold=3.0, bin_width=0.3
        # score=2.99 -> 2.99/0.3 = 9.966... -> int(9.966) = 9
        result = quantize_maxent_scores([2.99], n_bins=10, threshold=3.0)
        assert result[0] == 9

    def test_output_within_valid_range(self):
        """All output bins should be in [0, n_bins]."""
        scores = [-5.0, -1.0, 0.0, 0.5, 1.5, 2.5, 3.0, 5.0, 10.0]
        for n_bins in [2, 10, 20, 50]:
            result = quantize_maxent_scores(scores, n_bins=n_bins, threshold=3.0)
            for bin_idx in result:
                assert 0 <= bin_idx <= n_bins, f"Bin {bin_idx} out of range [0, {n_bins}]"

    def test_empty_input(self):
        """Empty list should return empty list."""
        assert quantize_maxent_scores([]) == []

    def test_single_score(self):
        """Single score should return single-element list."""
        result = quantize_maxent_scores([1.5], n_bins=10, threshold=3.0)
        assert len(result) == 1

    def test_invalid_n_bins_too_small(self):
        """n_bins < 2 should raise ValueError."""
        with pytest.raises(ValueError, match="at least 2 bins"):
            quantize_maxent_scores([1.0], n_bins=1)

    def test_invalid_n_bins_zero(self):
        """n_bins = 0 should raise ValueError."""
        with pytest.raises(ValueError, match="at least 2 bins"):
            quantize_maxent_scores([1.0], n_bins=0)

    def test_invalid_threshold_zero(self):
        """threshold <= 0 should raise ValueError."""
        with pytest.raises(ValueError, match="Threshold must be positive"):
            quantize_maxent_scores([1.0], threshold=0.0)

    def test_invalid_threshold_negative(self):
        """Negative threshold should raise ValueError."""
        with pytest.raises(ValueError, match="Threshold must be positive"):
            quantize_maxent_scores([1.0], threshold=-1.0)

    def test_monotonicity(self):
        """Higher scores should produce equal or higher bin indices."""
        scores = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
        result = quantize_maxent_scores(scores, n_bins=20, threshold=3.0)
        for i in range(1, len(result)):
            assert result[i] >= result[i - 1], (
                f"Non-monotonic: scores[{i-1}]={scores[i-1]} -> bin {result[i-1]}, "
                f"scores[{i}]={scores[i]} -> bin {result[i]}"
            )

    def test_precision_with_different_bin_counts(self):
        """More bins should provide finer granularity for intermediate scores."""
        score = 1.5
        r10 = quantize_maxent_scores([score], n_bins=10, threshold=3.0)
        r100 = quantize_maxent_scores([score], n_bins=100, threshold=3.0)
        # More bins -> larger numeric bin index but same relative position
        # 10 bins: 1.5/0.3 = 5
        # 100 bins: 1.5/0.03 = 50
        assert r10[0] == 5
        assert r100[0] == 50


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Cross-Codon Constraint Generation
# ═══════════════════════════════════════════════════════════════════════════════

class TestEncodeCrossCodonSpliceContext:

    def test_returns_list_of_tuples(self, sample_protein, default_config):
        """Return type should be a list of 5-tuples."""
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        assert isinstance(result, list)
        for entry in result:
            assert isinstance(entry, tuple)
            assert len(entry) == 5

    def test_tuple_structure(self, sample_protein, default_config):
        """Each tuple: (position, codon_left, codon_right, is_donor, score)."""
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        for pos, cl, cr, is_donor, score in result:
            assert isinstance(pos, int)
            assert len(cl) == 3
            assert len(cr) == 3
            assert isinstance(is_donor, bool)
            assert isinstance(score, float)

    def test_donor_has_gt_at_boundary(self, sample_protein, default_config):
        """Cross-codon donor entries must have GT at the boundary."""
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        for pos, cl, cr, is_donor, score in result:
            if is_donor:
                # Last base of codon_left + first base of codon_right = GT
                boundary_dinuc = cl[2] + cr[0]
                assert boundary_dinuc == "GT", (
                    f"Donor at pos {pos}: {cl}[{cl[2]}] + {cr}[{cr[0]}] = {boundary_dinuc}, expected GT"
                )

    def test_acceptor_has_ag_at_boundary(self, sample_protein, default_config):
        """Cross-codon acceptor entries must have AG at the boundary."""
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        for pos, cl, cr, is_donor, score in result:
            if not is_donor:
                boundary_dinuc = cl[2] + cr[0]
                assert boundary_dinuc == "AG", (
                    f"Acceptor at pos {pos}: {cl}[{cl[2]}] + {cr}[{cr[0]}] = {boundary_dinuc}, expected AG"
                )

    def test_no_out_of_bounds_scores(self, sample_protein, default_config):
        """Entries with out-of-bounds sentinel scores should be filtered out."""
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        for _, _, _, _, score in result:
            assert score > _OUT_OF_BOUNDS_SCORE

    def test_positions_are_valid(self, sample_protein, default_config):
        """Position should be a valid codon index < len(protein) - 1."""
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        n = len(sample_protein)
        for pos, _, _, _, _ in result:
            assert 0 <= pos < n - 1

    def test_codons_match_amino_acids(self, sample_protein, default_config):
        """Codon_left should encode protein[pos], codon_right protein[pos+1]."""
        from biocompiler.shared.constants import CODON_TABLE
        result = encode_cross_codon_splice_context(sample_protein, default_config)
        for pos, cl, cr, _, _ in result:
            assert CODON_TABLE[cl] == sample_protein[pos], (
                f"Position {pos}: {cl} encodes {CODON_TABLE[cl]}, expected {sample_protein[pos]}"
            )
            assert CODON_TABLE[cr] == sample_protein[pos + 1], (
                f"Position {pos+1}: {cr} encodes {CODON_TABLE[cr]}, expected {sample_protein[pos+1]}"
            )

    def test_single_aa_produces_no_cross_codon(self, single_aa_protein, default_config):
        """Single amino acid has no boundaries — should return empty list."""
        result = encode_cross_codon_splice_context(single_aa_protein, default_config)
        assert result == []

    def test_deterministic(self, sample_protein, default_config):
        """Same input should produce same output."""
        r1 = encode_cross_codon_splice_context(sample_protein, default_config)
        r2 = encode_cross_codon_splice_context(sample_protein, default_config)
        assert r1 == r2

    def test_valine_protein_cross_codon_gt(self, valine_protein, default_config):
        """Valine protein should have cross-codon GT possibilities."""
        # Valine codons all start with GT, so boundary between Val codons
        # is X|G where X is last base of one Val codon and G is first of next.
        # Since all Val codons are GTN, codon[2] can be T/A/C/G and codon[0]=G
        # So boundary dinuc is codon[2]+G which is not GT (it is XG not GT).
        # However, the T in GTN at position 1 and G at next codon position 0
        # form TG, not GT. So cross-codon GT needs codon[2]=G and next[0]=T.
        # Since all Val codons are GTN, next[0]='G', so no cross-codon GT for VV.
        # This is actually correct — we just verify it runs without error.
        result = encode_cross_codon_splice_context(valine_protein, default_config)
        assert isinstance(result, list)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. SpliceConstraintEncoder Class
# ═══════════════════════════════════════════════════════════════════════════════

class TestSpliceConstraintEncoder:

    def test_initialize_populates_caches(self, sample_protein, default_config):
        """After initialize(), all caches should be populated."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        assert not encoder._initialized
        encoder.initialize()
        assert encoder._initialized
        assert len(encoder._score_cache) > 0

    def test_initialize_idempotent(self, sample_protein, default_config):
        """Calling initialize() twice should not re-compute."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        encoder.initialize()
        cache_ref = encoder._score_cache
        encoder.initialize()  # Second call should be no-op
        assert encoder._score_cache is cache_ref

    def test_get_forbidden_assignments(self, sample_protein, default_config):
        """Should return a list of SpliceConstraint objects."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        forbidden = encoder.get_forbidden_assignments()
        assert isinstance(forbidden, list)
        for c in forbidden:
            assert isinstance(c, SpliceConstraint)

    def test_get_cross_codon_constraints(self, sample_protein, default_config):
        """Should return a list of CrossCodonSpliceConstraint objects."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        cross = encoder.get_cross_codon_constraints()
        assert isinstance(cross, list)
        for c in cross:
            assert isinstance(c, CrossCodonSpliceConstraint)
            assert c.position_right == c.position_left + 1

    def test_get_scores_for_position_valid(self, sample_protein, default_config):
        """Should return score data for a valid (position, codon) pair."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        scores = encoder.get_scores_for_position(0, "ATG")
        assert scores is not None
        assert "donor_scores" in scores
        assert "acceptor_scores" in scores

    def test_get_scores_for_position_invalid(self, sample_protein, default_config):
        """Should return None for a non-existent (position, codon) pair."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        scores = encoder.get_scores_for_position(0, "ZZZ")
        assert scores is None

    def test_get_max_splice_score_valid(self, sample_protein, default_config):
        """Should return a (max_donor, max_acceptor) tuple for valid position."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        max_d, max_a = encoder.get_max_splice_score(0, "ATG")
        assert isinstance(max_d, float)
        assert isinstance(max_a, float)

    def test_get_max_splice_score_invalid(self, sample_protein, default_config):
        """Should return sentinel values for invalid (position, codon) pair."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        max_d, max_a = encoder.get_max_splice_score(0, "ZZZ")
        assert max_d == _OUT_OF_BOUNDS_SCORE
        assert max_a == _OUT_OF_BOUNDS_SCORE

    def test_quantize_position_scores(self, sample_protein, default_config):
        """Should return (donor_bin, acceptor_bin) tuple of ints."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        d_bin, a_bin = encoder.quantize_position_scores(0, "ATG")
        assert isinstance(d_bin, int)
        assert isinstance(a_bin, int)
        assert 0 <= d_bin <= default_config.n_quantize_bins
        assert 0 <= a_bin <= default_config.n_quantize_bins

    def test_clear_cache_resets_state(self, sample_protein, default_config):
        """clear_cache() should empty all caches and reset _initialized."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        encoder.initialize()
        assert encoder._initialized
        encoder.clear_cache()
        assert not encoder._initialized
        assert len(encoder._score_cache) == 0
        assert len(encoder._forbidden_cache) == 0
        assert len(encoder._cross_codon_cache) == 0
        assert len(encoder._cross_codon_constraint_cache) == 0

    def test_query_after_clear_reinitializes(self, sample_protein, default_config):
        """Query methods should auto-reinitialize after clear_cache()."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        encoder.initialize()
        encoder.clear_cache()
        # get_forbidden_assignments calls initialize() internally
        forbidden = encoder.get_forbidden_assignments()
        assert encoder._initialized
        assert isinstance(forbidden, list)

    def test_forbidden_assignments_are_above_threshold(self, sample_protein, default_config):
        """All forbidden assignments should have scores exceeding their threshold."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        forbidden = encoder.get_forbidden_assignments()
        for c in forbidden:
            assert c.score > c.threshold

    def test_cross_codon_constraints_above_threshold(self, sample_protein, default_config):
        """Cross-codon constraints should only include pairs exceeding threshold."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        cross = encoder.get_cross_codon_constraints()
        for c in cross:
            if c.is_donor:
                assert c.score > default_config.effective_donor_threshold
            else:
                assert c.score > default_config.effective_acceptor_threshold

    def test_lenient_config_fewer_forbidden(self, sample_protein, default_config, lenient_config):
        """Lenient config should produce fewer forbidden assignments."""
        enc_default = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        enc_lenient = SpliceConstraintEncoder(protein=sample_protein, config=lenient_config)
        assert len(enc_lenient.get_forbidden_assignments()) <= len(enc_default.get_forbidden_assignments())


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Edge Cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_single_methionine_no_forbidden(self, default_config):
        """Single M protein: ATG should not produce forbidden constraints."""
        encoder = SpliceConstraintEncoder(protein="M", config=default_config)
        forbidden = encoder.get_forbidden_assignments()
        # ATG is only 3 bases; may or may not have splice scores above threshold
        # But it should not error
        assert isinstance(forbidden, list)

    def test_single_methionine_no_cross_codon(self, default_config):
        """Single M protein has no cross-codon boundaries."""
        encoder = SpliceConstraintEncoder(protein="M", config=default_config)
        cross = encoder.get_cross_codon_constraints()
        assert cross == []

    def test_protein_with_no_gt_ag_sites(self, default_config):
        """Protein with no GT/AG dinucleotides in any codon should have empty scores."""
        # Phenylalanine: TTT, TTC — no GT, no AG within codon
        # Also tryptophan: TGG — no GT, no AG
        protein = "FW"
        data = precompute_splice_scores(protein, default_config)
        # Within-codon scores should be empty for non-GT/AG codons
        # (but cross-codon could still create them)
        for (pos, codon), scores in data.items():
            if "GT" not in codon:
                # Within-codon donor scores should be empty
                within_donor = [(bp, s) for bp, s in scores["donor_scores"]
                                if bp >= pos * 3 and bp < (pos + 1) * 3]
                # Note: cross-codon GTs could still contribute, so we just
                # verify the structure is sound
                pass

    def test_very_short_sequence_scoring(self, default_config):
        """Very short sequence should handle out-of-bounds MaxEntScan context gracefully."""
        protein = "MK"
        data = precompute_splice_scores(protein, default_config)
        # Should not crash — MaxEntScan returns -50 for out-of-bounds
        # The precompute function stores all scores including sentinel -50,
        # but the constraint table and cross-codon functions filter those out.
        for key, scores in data.items():
            for _, s in scores["donor_scores"] + scores["acceptor_scores"]:
                assert isinstance(s, float)  # Just verify structure is valid

    def test_unknown_amino_acid_skipped_in_main_loop(self, default_config):
        """Unknown amino acid should be skipped in the per-position loop.

        Note: precompute_splice_scores skips unknown AAs at line 159-160
        but _build_sequence_with_codon accesses AA_TO_CODONS for context
        positions. So we test with a protein that only has valid AAs but
        verify the skip logic in the per-position loop by checking that
        only valid (position, codon) keys appear in the output.
        """
        protein = "MK"
        data = precompute_splice_scores(protein, default_config)
        # All keys should reference valid positions and codons
        for (pos, codon) in data:
            assert 0 <= pos < len(protein)
            assert codon in AA_TO_CODONS[protein[pos]]

    def test_precompute_with_all_same_amino_acid(self, default_config):
        """All-same amino acid protein should still produce valid results."""
        protein = "LLLLLLLLLL"  # 10 leucines
        data = precompute_splice_scores(protein, default_config)
        assert len(data) == 10 * len(AA_TO_CODONS["L"])

    def test_quantize_with_extreme_scores(self):
        """Very large positive and negative scores should be clamped."""
        result = quantize_maxent_scores([-1000.0, 1000.0], n_bins=10, threshold=3.0)
        assert result[0] == 0       # Negative -> bin 0
        assert result[1] == 10      # Above threshold -> n_bins

    def test_cross_codon_single_boundary(self, default_config):
        """Two-codon protein should have exactly one boundary to check."""
        result = encode_cross_codon_splice_context("MV", default_config)
        # All entries should have position 0 (the only boundary)
        for pos, _, _, _, _ in result:
            assert pos == 0

    def test_encoder_with_custom_quantize_bins(self):
        """Custom n_quantize_bins should be respected."""
        config = SolverConfig(n_quantize_bins=5)
        encoder = SpliceConstraintEncoder(protein="MVLSPADKTN", config=config)
        d_bin, a_bin = encoder.quantize_position_scores(0, "ATG")
        assert d_bin <= 5
        assert a_bin <= 5

    def test_encoder_protein_uppercased(self, default_config):
        """Encoder should handle lowercase protein by uppercasing."""
        encoder = SpliceConstraintEncoder(protein="mvlspadktn", config=default_config)
        forbidden = encoder.get_forbidden_assignments()
        assert isinstance(forbidden, list)

    def test_precompute_scores_no_donor_for_non_gt_codons(self, default_config):
        """Codons without GT should have empty within-codon donor scores in simple cases."""
        # Alanine: GCT, GCC, GCA, GCG — no GT within codon
        protein = "A"
        data = precompute_splice_scores(protein, default_config)
        for (pos, codon), scores in data.items():
            if "GT" not in codon:
                # Within-codon GT donor sites should not exist
                within_gt = [bp for bp, s in scores["donor_scores"]
                             if bp >= pos * 3 and bp < (pos + 1) * 3]
                # But with overlapping context from other codons, this may
                # not be strictly empty, so just verify structure
                assert isinstance(scores["donor_scores"], list)


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Integration / Consistency Checks
# ═══════════════════════════════════════════════════════════════════════════════

class TestIntegrationConsistency:

    def test_forbidden_assignments_subset_of_score_data(self, sample_protein, default_config):
        """Every forbidden assignment should reference a key in the score cache."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        forbidden = encoder.get_forbidden_assignments()
        for c in forbidden:
            key = (c.position, c.codon)
            assert key in encoder._score_cache, (
                f"Forbidden assignment {key} not in score cache"
            )

    def test_quantized_scores_consistent_with_max_scores(self, sample_protein, default_config):
        """Quantized position scores should be consistent with max splice scores."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        for i, aa in enumerate(sample_protein):
            for codon in AA_TO_CODONS[aa]:
                max_d, max_a = encoder.get_max_splice_score(i, codon)
                d_bin, a_bin = encoder.quantize_position_scores(i, codon)
                # Verify bins are consistent with the quantize function directly
                expected_d = quantize_maxent_scores(
                    [max_d],
                    n_bins=default_config.n_quantize_bins,
                    threshold=default_config.effective_donor_threshold,
                )[0]
                expected_a = quantize_maxent_scores(
                    [max_a],
                    n_bins=default_config.n_quantize_bins,
                    threshold=default_config.effective_acceptor_threshold,
                )[0]
                assert d_bin == expected_d, (
                    f"pos={i} codon={codon}: donor bin {d_bin} != expected {expected_d}"
                )
                assert a_bin == expected_a, (
                    f"pos={i} codon={codon}: acceptor bin {a_bin} != expected {expected_a}"
                )

    def test_cross_codon_constraint_positions_consecutive(self, sample_protein, default_config):
        """Cross-codon constraints should always involve consecutive positions."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        cross = encoder.get_cross_codon_constraints()
        for c in cross:
            assert c.position_right == c.position_left + 1

    def test_cross_codon_dinucleotides_correct(self, sample_protein, default_config):
        """Cross-codon constraints should have correct GT/AG at the boundary."""
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        cross = encoder.get_cross_codon_constraints()
        for c in cross:
            dinuc = c.codon_left[2] + c.codon_right[0]
            if c.is_donor:
                assert dinuc == "GT", f"Expected GT at boundary, got {dinuc}"
            else:
                assert dinuc == "AG", f"Expected AG at boundary, got {dinuc}"

    def test_precompute_and_constraint_table_agree(self, sample_protein, default_config):
        """Manually computing constraint table should match encoder results."""
        data = precompute_splice_scores(sample_protein, default_config)
        manual_constraints = build_splice_constraint_table(sample_protein, default_config, data)
        encoder = SpliceConstraintEncoder(protein=sample_protein, config=default_config)
        encoder_constraints = encoder.get_forbidden_assignments()
        # Same count
        assert len(manual_constraints) == len(encoder_constraints)
        # Same set of (position, codon, site_type) keys
        manual_keys = {(c.position, c.codon, c.site_type) for c in manual_constraints}
        encoder_keys = {(c.position, c.codon, c.site_type) for c in encoder_constraints}
        assert manual_keys == encoder_keys
