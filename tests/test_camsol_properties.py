"""
Property-based tests for camsol.py using Hypothesis
====================================================
Covers three core properties:

  1. Solubility scores are bounded in [-3, +3]
  2. Same sequence always produces the same score (determinism)
  3. Valid protein sequences do not crash compute_intrinsic_solubility

Import from biocompiler.camsol.
"""

import pytest
pytest.importorskip("hypothesis")
from hypothesis import given, settings, assume
from hypothesis import strategies as st
pytest.importorskip("hypothesis")

from biocompiler.engines.camsol import (
    CamSolResult,
    SolubilityResult,
    compute_intrinsic_solubility,
    compute_solubility,
    classify_solubility,
    predict_idp,
    select_hydropathy_scale,
    clear_cache,
    CAMSOL_HYDROPATHY,
    URRY_HYDROPATHY,
    CAMSOL_CHARGE,
    CAMSOL_ALPHA_HELIX,
    CAMSOL_BETA_STRAND,
    HydrophobicityScale,
)
from biocompiler.shared.exceptions import CamSolError


# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

STANDARD_AAS = sorted(CAMSOL_HYDROPATHY.keys())
"""The 20 standard amino acids (from CAMSOL_HYDROPATHY keys)."""

VALID_SOLUBILITY_CLASSES = frozenset(
    {"highly_soluble", "soluble", "marginally_soluble", "insoluble"}
)

# CamSol documented score bounds
SCORE_MIN = -3.0
SCORE_MAX = 3.0


# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

amino_acid = st.sampled_from(STANDARD_AAS)
"""Strategy for a single standard amino acid."""


@st.composite
def protein_sequence(draw, min_size: int = 1, max_size: int = 200) -> str:
    """Generate a valid protein sequence of standard amino acids.

    The minimum size defaults to 1 since compute_intrinsic_solubility
    requires at least 1 residue (empty string raises CamSolError).
    """
    length = draw(st.integers(min_value=min_size, max_value=max_size))
    residues = draw(st.lists(amino_acid, min_size=length, max_size=length))
    return "".join(residues)


@st.composite
def hydrophobic_protein(draw, min_size: int = 10, max_size: int = 100) -> str:
    """Generate a protein rich in hydrophobic residues (V, I, L, F, W, M)."""
    hydrophobic_aas = st.sampled_from(list("VLIFWM"))
    length = draw(st.integers(min_value=min_size, max_value=max_size))
    residues = draw(st.lists(hydrophobic_aas, min_size=length, max_size=length))
    return "".join(residues)


@st.composite
def hydrophilic_protein(draw, min_size: int = 10, max_size: int = 100) -> str:
    """Generate a protein rich in charged/hydrophilic residues (K, E, D, R, Q, N)."""
    hydrophilic_aas = st.sampled_from(list("KEDRQN"))
    length = draw(st.integers(min_value=min_size, max_value=max_size))
    residues = draw(st.lists(hydrophilic_aas, min_size=length, max_size=length))
    return "".join(residues)


@st.composite
def idp_like_protein(draw, min_size: int = 20, max_size: int = 100) -> str:
    """Generate a protein likely predicted as intrinsically disordered.

    Uses disorder-promoting residues (P, E, S, Q, K, R, D, G) with
    few order-promoting residues to satisfy predict_idp thresholds.
    """
    disorder_aas = st.sampled_from(list("PESQKRDG"))
    length = draw(st.integers(min_value=min_size, max_value=max_size))
    residues = draw(st.lists(disorder_aas, min_size=length, max_size=length))
    return "".join(residues)


hydrophobicity_scale_name = st.sampled_from(["auto", "wimley_white", "urry"])
"""Strategy for the hydrophobicity_scale parameter."""


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _clear_camsol_cache():
    """Clear the CamSol cache before and after each test to avoid
    cached results masking non-determinism."""
    clear_cache()
    yield
    clear_cache()


# ===========================================================================
# Property 1: Solubility scores are bounded in [-3, +3]
# ===========================================================================

class TestSolubilityScoresBounded:
    """All solubility scores returned by compute_intrinsic_solubility
    must lie within the documented range of [-3, +3]."""

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_intrinsic_score_bounded(self, seq):
        """compute_intrinsic_solubility returns intrinsic_score in [-3, +3]
        for any valid protein sequence."""
        result = compute_intrinsic_solubility(seq)
        assert SCORE_MIN <= result.intrinsic_score <= SCORE_MAX, (
            f"intrinsic_score {result.intrinsic_score} out of bounds "
            f"[{SCORE_MIN}, {SCORE_MAX}] for sequence starting with {seq[:10]!r}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_primary_score_bounded(self, seq):
        """The primary_score (overall_score alias) is also in [-3, +3]."""
        result = compute_intrinsic_solubility(seq)
        assert SCORE_MIN <= result.primary_score <= SCORE_MAX, (
            f"primary_score {result.primary_score} out of bounds "
            f"[{SCORE_MIN}, {SCORE_MAX}]"
        )

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_per_residue_scores_bounded(self, seq):
        """Each per-residue score should be a finite float.
        (Per-residue scores may exceed [-3, +3] before averaging,
        but they must be finite numbers.)"""
        result = compute_intrinsic_solubility(seq)
        for i, score in enumerate(result.per_residue_scores):
            assert isinstance(score, float), (
                f"Per-residue score at position {i} is {type(score)}, not float"
            )
            assert abs(score) != float("inf"), (
                f"Per-residue score at position {i} is infinite: {score}"
            )
            assert score == score, (  # NaN check
                f"Per-residue score at position {i} is NaN"
            )

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_compute_solubility_score_bounded(self, seq):
        """compute_solubility (unified API) also returns scores in [-3, +3]."""
        result = compute_solubility(seq)
        assert SCORE_MIN <= result.primary_score <= SCORE_MAX, (
            f"compute_solubility primary_score {result.primary_score} "
            f"out of bounds [{SCORE_MIN}, {SCORE_MAX}]"
        )

    @given(seq=hydrophobic_protein())
    @settings(max_examples=50, deadline=5000)
    def test_hydrophobic_score_bounded(self, seq):
        """Hydrophobic protein scores are still within bounds."""
        result = compute_intrinsic_solubility(seq)
        assert SCORE_MIN <= result.intrinsic_score <= SCORE_MAX, (
            f"Hydrophobic protein score {result.intrinsic_score} out of bounds"
        )

    @given(seq=hydrophilic_protein())
    @settings(max_examples=50, deadline=5000)
    def test_hydrophilic_score_bounded(self, seq):
        """Hydrophilic protein scores are still within bounds."""
        result = compute_intrinsic_solubility(seq)
        assert SCORE_MIN <= result.intrinsic_score <= SCORE_MAX, (
            f"Hydrophilic protein score {result.intrinsic_score} out of bounds"
        )

    @given(score=st.floats(min_value=-5.0, max_value=5.0))
    @settings(max_examples=200)
    def test_classify_solubility_returns_valid_class(self, score):
        """classify_solubility always returns a valid class string
        regardless of the input score (even outside [-3, +3])."""
        result = classify_solubility(score)
        assert result in VALID_SOLUBILITY_CLASSES, (
            f"classify_solubility({score}) = {result!r}, "
            f"expected one of {VALID_SOLUBILITY_CLASSES}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_classification_matches_score(self, seq):
        """The solubility_class returned must be consistent with the score
        via classify_solubility."""
        result = compute_intrinsic_solubility(seq)
        expected_class = classify_solubility(result.intrinsic_score)
        assert result.classification == expected_class, (
            f"classification={result.classification!r} != "
            f"classify_solubility({result.intrinsic_score})={expected_class!r}"
        )


# ===========================================================================
# Property 2: Same sequence always produces same score (determinism)
# ===========================================================================

class TestSolubilityDeterminism:
    """Calling compute_intrinsic_solubility with the same inputs must
    always produce the same result.  This holds even across cache
    clears (pure-function property)."""

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_same_sequence_same_score(self, seq):
        """Two calls with the same sequence produce the same intrinsic_score."""
        result1 = compute_intrinsic_solubility(seq)
        clear_cache()  # Force recomputation, not just cache hit
        result2 = compute_intrinsic_solubility(seq)
        assert result1.intrinsic_score == result2.intrinsic_score, (
            f"Non-deterministic score for {seq[:10]!r}: "
            f"{result1.intrinsic_score} != {result2.intrinsic_score}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_same_sequence_same_per_residue(self, seq):
        """Two calls produce identical per_residue_scores."""
        result1 = compute_intrinsic_solubility(seq)
        clear_cache()
        result2 = compute_intrinsic_solubility(seq)
        assert result1.per_residue_scores == result2.per_residue_scores, (
            f"Non-deterministic per-residue scores for {seq[:10]!r}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_same_sequence_same_aggregation_regions(self, seq):
        """Two calls produce identical aggregation_prone_regions."""
        result1 = compute_intrinsic_solubility(seq)
        clear_cache()
        result2 = compute_intrinsic_solubility(seq)
        assert result1.aggregation_prone_regions == result2.aggregation_prone_regions, (
            f"Non-deterministic aggregation regions for {seq[:10]!r}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_same_sequence_same_class(self, seq):
        """Two calls produce the same solubility classification."""
        result1 = compute_intrinsic_solubility(seq)
        clear_cache()
        result2 = compute_intrinsic_solubility(seq)
        assert result1.classification == result2.classification, (
            f"Non-deterministic classification for {seq[:10]!r}: "
            f"{result1.classification} != {result2.classification}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_same_sequence_same_recommendations(self, seq):
        """Two calls produce identical recommendations."""
        result1 = compute_intrinsic_solubility(seq)
        clear_cache()
        result2 = compute_intrinsic_solubility(seq)
        assert result1.recommendations == result2.recommendations, (
            f"Non-deterministic recommendations for {seq[:10]!r}"
        )

    @given(
        seq=protein_sequence(min_size=1, max_size=50),
        scale=hydrophobicity_scale_name,
    )
    @settings(max_examples=50, deadline=5000)
    def test_determinism_with_explicit_scale(self, seq, scale):
        """Determinism holds even when specifying a hydrophobicity scale."""
        result1 = compute_intrinsic_solubility(seq, hydrophobicity_scale=scale)
        clear_cache()
        result2 = compute_intrinsic_solubility(seq, hydrophobicity_scale=scale)
        assert result1.intrinsic_score == result2.intrinsic_score, (
            f"Non-deterministic with scale={scale!r}: "
            f"{result1.intrinsic_score} != {result2.intrinsic_score}"
        )

    @given(seq=protein_sequence(min_size=1, max_size=50))
    @settings(max_examples=30, deadline=5000)
    def test_compute_solubility_determinism(self, seq):
        """The unified compute_solubility API is also deterministic."""
        result1 = compute_solubility(seq)
        clear_cache()
        result2 = compute_solubility(seq)
        assert result1.primary_score == result2.primary_score, (
            f"compute_solubility non-deterministic: "
            f"{result1.primary_score} != {result2.primary_score}"
        )

    # Concrete (non-Hypothesis) cross-check

    def test_determinism_concrete_known_proteins(self):
        """Determinism on concrete known protein sequences."""
        proteins = {
            "soluble": "MSEKKDKKEKEKKDEKKDEEKKDEEKKDESKKDEEKKDEEKKDESKKDEEKKDEEKK",
            "insoluble": "MVVVIIVVVLLLFLLLLFFFFWWWAAAIIIMMM",
            "balanced": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAVDILSKKGDVQVIK",
        }
        for label, seq in proteins.items():
            r1 = compute_intrinsic_solubility(seq)
            clear_cache()
            r2 = compute_intrinsic_solubility(seq)
            assert r1.intrinsic_score == r2.intrinsic_score, (
                f"Non-deterministic for {label} protein"
            )


# ===========================================================================
# Property 3: Valid protein sequences do not crash
# ===========================================================================

class TestValidSequencesNoCrash:
    """compute_intrinsic_solubility must not raise exceptions for any
    valid protein sequence (standard 20 amino acids only)."""

    @given(seq=protein_sequence(min_size=1, max_size=300))
    @settings(max_examples=200, deadline=10000)
    def test_no_crash_random_sequence(self, seq):
        """Any valid protein sequence of 1–300 residues completes without
        raising an exception."""
        result = compute_intrinsic_solubility(seq)
        assert result.success is True, (
            f"compute_intrinsic_solubility failed for {seq[:10]!r}...: "
            f"error={result.error!r}"
        )
        assert isinstance(result, CamSolResult)

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_no_crash_with_each_scale(self, seq):
        """No crash regardless of which hydrophobicity scale is selected."""
        for scale in ("auto", "wimley_white", "urry"):
            result = compute_intrinsic_solubility(
                seq, hydrophobicity_scale=scale
            )
            assert result.success is True, (
                f"Failed with scale={scale!r}: {result.error!r}"
            )

    @given(seq=protein_sequence(min_size=1, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_no_crash_hydrophobicity_enum(self, seq):
        """No crash when passing HydrophobicityScale enum members."""
        for scale in HydrophobicityScale:
            result = compute_intrinsic_solubility(
                seq, hydrophobicity_scale=scale
            )
            assert result.success is True

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_result_fields_well_formed(self, seq):
        """Result has all expected fields with correct types."""
        result = compute_intrinsic_solubility(seq)
        # Sequence matches input
        assert result.sequence == seq
        # Per-residue scores length matches sequence length
        assert len(result.per_residue_scores) == len(seq), (
            f"per_residue_scores length {len(result.per_residue_scores)} "
            f"!= sequence length {len(seq)}"
        )
        # Aggregation regions have valid coordinates
        for start, end, avg in result.aggregation_prone_regions:
            assert 0 <= start < end <= len(seq), (
                f"Invalid aggregation region ({start}, {end}) "
                f"for sequence length {len(seq)}"
            )
        # Classification is a valid string
        assert result.classification in VALID_SOLUBILITY_CLASSES
        # Method is correct
        assert result.method == "camsol_intrinsic"
        # Engine name
        assert result.engine_name == "camsol"
        # Structural score is None for intrinsic computation
        assert result.structural_score is None
        # Hydrophobicity scale used is valid
        assert result.hydrophobicity_scale_used in ("wimley_white", "urry")

    @given(seq=protein_sequence(min_size=1, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_compute_solubility_no_crash(self, seq):
        """The unified compute_solubility API also never crashes on
        valid sequences (without PDB)."""
        result = compute_solubility(seq)
        assert result.success is True
        assert isinstance(result, CamSolResult)

    @given(seq=idp_like_protein(min_size=15, max_size=100))
    @settings(max_examples=50, deadline=5000)
    def test_idp_sequences_no_crash(self, seq):
        """Sequences likely predicted as IDPs complete successfully
        (exercise the Urry scale auto-selection path)."""
        result = compute_intrinsic_solubility(seq)
        assert result.success is True
        # If predicted as IDP, the Urry scale should have been used
        if predict_idp(seq):
            assert result.hydrophobicity_scale_used == "urry", (
                f"IDP sequence got scale={result.hydrophobicity_scale_used!r} "
                f"instead of 'urry'"
            )

    @given(seq=protein_sequence(min_size=1, max_size=5))
    @settings(max_examples=50, deadline=5000)
    def test_very_short_sequences_no_crash(self, seq):
        """Even very short sequences (1–5 residues) do not crash."""
        result = compute_intrinsic_solubility(seq)
        assert result.success is True

    @given(single_aa=amino_acid)
    @settings(max_examples=20, deadline=5000)
    def test_single_residue_no_crash(self, single_aa):
        """A single amino acid does not crash."""
        result = compute_intrinsic_solubility(single_aa)
        assert result.success is True
        assert len(result.per_residue_scores) == 1

    @given(aa=amino_acid, count=st.integers(min_value=1, max_value=50))
    @settings(max_examples=40, deadline=5000)
    def test_homopolymer_no_crash(self, aa, count):
        """Homopolymers (repeated single AA) never crash."""
        seq = aa * count
        result = compute_intrinsic_solubility(seq)
        assert result.success is True

    def test_empty_sequence_raises_camsol_error(self):
        """Empty sequence must raise CamSolError (not a bare exception)."""
        with pytest.raises(CamSolError):
            compute_intrinsic_solubility("")

    @given(invalid_char=st.sampled_from(list("BJOUXZ")))
    @settings(max_examples=6)
    def test_invalid_aa_raises_camsol_error(self, invalid_char):
        """Non-standard amino acid codes must raise CamSolError."""
        seq = "MAK" + invalid_char + "EL"
        with pytest.raises(CamSolError):
            compute_intrinsic_solubility(seq)

    def test_invalid_scale_raises_value_error(self):
        """Invalid hydrophobicity_scale value raises ValueError."""
        with pytest.raises(ValueError, match="Unknown hydrophobicity scale"):
            compute_intrinsic_solubility("MAKEL", hydrophobicity_scale="invalid")


# ===========================================================================
# Additional property: predict_idp is consistent with select_hydropathy_scale
# ===========================================================================

class TestIDPPredictionConsistency:
    """predict_idp and select_hydropathy_scale must be consistent:
    when predict_idp returns True for "auto" mode, the Urry scale
    must be selected."""

    @given(seq=protein_sequence(min_size=10, max_size=200))
    @settings(max_examples=100, deadline=5000)
    def test_idp_consistent_with_scale_selection(self, seq):
        """If predict_idp(seq) is True, select_hydropathy_scale with
        'auto' must choose the Urry scale."""
        is_idp = predict_idp(seq)
        _, scale_name = select_hydropathy_scale(seq, "auto")
        if is_idp:
            assert scale_name == "urry", (
                f"predict_idp=True but scale={scale_name!r} "
                f"(expected 'urry') for {seq[:15]!r}..."
            )
        else:
            assert scale_name == "wimley_white", (
                f"predict_idp=False but scale={scale_name!r} "
                f"(expected 'wimley_white') for {seq[:15]!r}..."
            )

    @given(seq=protein_sequence(min_size=1, max_size=9))
    @settings(max_examples=30, deadline=5000)
    def test_short_sequence_not_idp(self, seq):
        """Sequences shorter than the IDP minimum length (10) are
        never predicted as IDP."""
        is_idp = predict_idp(seq)
        assert is_idp is False, (
            f"Sequence of length {len(seq)} predicted as IDP "
            f"(min length is 10)"
        )

    @given(seq=protein_sequence(min_size=10, max_size=200))
    @settings(max_examples=50, deadline=5000)
    def test_explicit_scale_overrides_auto(self, seq):
        """When an explicit scale is given, it overrides the auto selection."""
        _, scale = select_hydropathy_scale(seq, "urry")
        assert scale == "urry"
        _, scale = select_hydropathy_scale(seq, "wimley_white")
        assert scale == "wimley_white"
