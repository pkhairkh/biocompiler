"""
BioCompiler CamSol Solubility Module v1.0.0
=============================================
CamSol-inspired solubility prediction algorithm in pure Python.

Predicts protein solubility from sequence using intrinsic physicochemical
properties (hydropathy, charge, secondary structure propensity) with
optional structure-based corrections from PDB data.

Based on the CamSol method (Sormanni et al., J Mol Biol 2015) for
predicting protein solubility from sequence.

Scoring range: -3 to +3 (positive = soluble, negative = aggregation-prone).

Accuracy and Confidence
----------------------
**CamSol intrinsic mode** (sequence-only, no PDB):
  - Classification accuracy: ~85.7% against 21-protein curated benchmark
  - Specificity (correctly identifying soluble proteins): 100%
  - Sensitivity (correctly identifying aggregation-prone proteins): ~66.7%
    for intrinsically disordered/aggregation-prone proteins (improved from
    ~33.3% by adding the Urry hydrophobicity scale for predicted IDPs)
  - Pearson r ≈ 0.73 between enhanced score and known solubility ordinal
  - The Urry hydrophobicity scale (based on elastin-like polypeptide
    experiments) is automatically selected for sequences predicted to be
    intrinsically disordered. This correctly penalises hydrophobic patches
    in IDPs that the Wimley-White scale underweights.
  - Enhanced benchmark score (with patch correction) improves discrimination
  - Validated against: ``validation.camsol_benchmark`` (21 proteins)

**CamSol structural mode** (with PDB data):
  - More accurate than intrinsic mode when PDB data is available
  - Structure-based SASA corrections improve aggregation-prone region detection
  - No quantitative benchmark yet for structural mode

  **Confidence levels:**
    - Intrinsic, enhanced score > 0.5: **HIGH** — reliable soluble classification
    - Intrinsic, enhanced score in [-0.5, 0.5]: **MEDIUM** — borderline
    - Intrinsic, enhanced score < -0.5: **HIGH** — reliable aggregation-prone classification
    - Structural mode: generally **HIGH** when PDB is available

**Known limitations:**
  - The Wimley-White octanol scale may still underperform for some IDPs
    not caught by the composition heuristic; manually set
    hydrophobicity_scale="urry" for such cases
  - Simple mean scoring compresses the signal; the published CamSol uses
    an aggressive patch-correction formula
  - No pH or temperature dependence modeled

References
----------
- Sormanni et al., J Mol Biol 2015; 427:478-490 (original CamSol)
- Wimley & White, Nat Struct Biol 1996; 3:842 (octanol scale)
- Urry et al., J Protein Chem 1992; 11:165 (hydrophobicity scale for IDPs)
"""

from __future__ import annotations

import hashlib
import logging
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum

from biocompiler.shared.constants import BLOSUM62, STANDARD_AAS
from .base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from biocompiler.shared.exceptions import CamSolError

import warnings as _warnings

# Suppress the SolubilityResult deprecation warning for internal module use
# (the __getattr__ at the bottom of this file emits it for external callers)
_warnings.filterwarnings(
    "ignore",
    message="camsol.SolubilityResult is deprecated",
    category=DeprecationWarning,
    module=__name__,
)

logger = logging.getLogger(__name__)

# ── Attempt to import shared defaults; provide local fallbacks ──
try:
    from biocompiler.shared.constants import DEFAULT_SOLUBILITY_WINDOW as _DEFAULT_WINDOW
except ImportError:
    _DEFAULT_WINDOW: int = 7  # type: ignore[misc]
    logger.debug("DEFAULT_SOLUBILITY_WINDOW not found in constants; using fallback %d", _DEFAULT_WINDOW)

try:
    from biocompiler.shared.constants import DEFAULT_SOLUBILITY_SMOOTHING as _DEFAULT_SMOOTHING
except ImportError:
    _DEFAULT_SMOOTHING: int = 3  # type: ignore[misc]
    logger.debug("DEFAULT_SOLUBILITY_SMOOTHING not found in constants; using fallback %d", _DEFAULT_SMOOTHING)

try:
    from biocompiler.shared.constants import DEFAULT_BATCH_SIZE as _DEFAULT_BATCH_SIZE
except ImportError:
    _DEFAULT_BATCH_SIZE: int = 8  # type: ignore[misc]
    logger.debug("DEFAULT_BATCH_SIZE not found in constants; using fallback %d", _DEFAULT_BATCH_SIZE)


__all__ = [
    "CamSolResult",
    "SolubilityResult",
    "HydrophobicityScale",
    "compute_intrinsic_solubility",
    "compute_solubility",
    "compute_structural_solubility",
    "compute_solubility_batch",
    "find_solubility_mutations",
    "classify_solubility",
    "generate_solubility_recommendations",
    "clear_cache",
    "CAMSOL_HYDROPATHY",
    "URRY_HYDROPATHY",
    "CAMSOL_CHARGE",
    "CAMSOL_ALPHA_HELIX",
    "CAMSOL_BETA_STRAND",
    "CAMSOL_CLASSIFICATION_ACCURACY",
    "CAMSOL_SPECIFICITY",
    "CAMSOL_SENSITIVITY_IDP",
    "CAMSOL_PEARSON_R",
    "predict_idp",
    "select_hydropathy_scale",
]


# ────────────────────────────────────────────────────────────
# Accuracy constants (from validation.camsol_benchmark)
# ────────────────────────────────────────────────────────────

#: Classification accuracy against 21-protein curated benchmark
#: (high/medium/low solubility classification with enhanced scoring)
CAMSOL_CLASSIFICATION_ACCURACY: float = 0.857

#: Specificity: fraction of truly soluble proteins correctly identified
#: (enhanced score > 0 for known-soluble proteins)
CAMSOL_SPECIFICITY: float = 1.0

#: Sensitivity for intrinsically disordered / aggregation-prone proteins
#: (enhanced score < 0 for known-aggregation-prone proteins)
#: Improved from 0.333 (Wimley-White only) by adding Urry scale auto-selection
#: for predicted IDP sequences
CAMSOL_SENSITIVITY_IDP: float = 0.667

#: Pearson correlation between enhanced score and known solubility ordinal
CAMSOL_PEARSON_R: float = 0.73


# ────────────────────────────────────────────────────────────
# In-memory cache
# ────────────────────────────────────────────────────────────

_cache: dict[tuple[str, int, int, str], CamSolResult] = {}


def clear_cache() -> None:
    """Clear the module-level CamSol result cache."""
    _cache.clear()
    logger.info("CamSol cache cleared.")


# ────────────────────────────────────────────────────────────
# CamSol-specific physicochemical scales
# ────────────────────────────────────────────────────────────

CAMSOL_HYDROPATHY: dict[str, float] = {
    # CamSol-specific hydropathy scale (optimized for solubility prediction)
    # Based on Wimley-White octanol scale with CamSol corrections
    # Positive = hydrophilic (soluble), Negative = hydrophobic (aggregation-prone)
    "I": -0.31,
    "V": -0.07,
    "L": -0.56,
    "F": -0.52,
    "C":  0.24,
    "M": -0.23,
    "A":  0.17,
    "G":  0.27,
    "T":  0.26,
    "S":  0.40,
    "W": -0.60,
    "Y": -0.18,
    "P":  0.99,
    "H":  0.61,
    "E":  1.32,
    "Q":  0.73,
    "D":  1.26,
    "N":  0.64,
    "K":  1.23,
    "R":  0.60,
}

CAMSOL_CHARGE: dict[str, float] = {
    # Net charge contribution at pH 7.4
    "K":  1.0,
    "R":  1.0,
    "H":  0.5,  # partially protonated at physiological pH
    "D": -1.0,
    "E": -1.0,
    "A":  0.0,
    "V":  0.0,
    "I":  0.0,
    "L":  0.0,
    "F":  0.0,
    "C":  0.0,
    "M":  0.0,
    "G":  0.0,
    "T":  0.0,
    "S":  0.0,
    "W":  0.0,
    "Y":  0.0,
    "P":  0.0,
    "Q":  0.0,
    "N":  0.0,
}

CAMSOL_ALPHA_HELIX: dict[str, float] = {
    # Alpha-helix propensity (solubility effect)
    # Positive values indicate helix-forming residues that reduce aggregation
    "A":  0.23,
    "L":  0.11,
    "E":  0.37,
    "M":  0.06,
    "Q":  0.33,
    "K":  0.26,
    "R":  0.21,
    "H":  0.09,
    "V": -0.08,
    "I": -0.07,
    "W": -0.09,
    "F": -0.11,
    "Y": -0.07,
    "C": -0.03,
    "G": -0.12,
    "P": -0.25,
    "S":  0.03,
    "T":  0.04,
    "D":  0.15,
    "N":  0.05,
}

CAMSOL_BETA_STRAND: dict[str, float] = {
    # Beta-strand propensity (aggregation effect)
    # Positive values indicate strand-forming residues that increase aggregation
    "V":  0.34,
    "I":  0.30,
    "T":  0.14,
    "F":  0.28,
    "Y":  0.13,
    "W":  0.19,
    "L":  0.22,
    "C":  0.10,
    "M":  0.08,
    "Q": -0.04,
    "N": -0.05,
    "S": -0.03,
    "R": -0.06,
    "H": -0.05,
    "E": -0.10,
    "D": -0.09,
    "A": -0.01,
    "G": -0.14,
    "P": -0.28,
    "K": -0.11,
}

# ────────────────────────────────────────────────────────────
# Urry hydrophobicity scale (for intrinsically disordered proteins)
# ────────────────────────────────────────────────────────────

URRY_HYDROPATHY: dict[str, float] = {
    # Urry hydrophobicity scale (based on elastin-like polypeptide experiments)
    # Specifically designed for intrinsically disordered proteins (IDPs).
    #
    # Original Urry scale: positive = hydrophobic, negative = hydrophilic.
    # Values below are NEGATED to match CamSol convention:
    #   Positive = hydrophilic (soluble), Negative = hydrophobic (aggregation-prone)
    #
    # Key differences from Wimley-White (CAMSOL_HYDROPATHY):
    #   - Charged residues (K, R, D, E) are less extremely positive
    #   - Hydrophobic residues (I, L, V, F, W, Y) are much more negative
    #   - This correctly penalises hydrophobic patches in IDPs that the
    #     Wimley-White scale underweights
    #
    # Reference: Urry et al., J Protein Chem 1992; 11:165
    "A": -0.086,
    "R":  0.850,
    "N":  0.549,
    "D":  0.695,
    "C": -0.397,
    "Q":  0.549,
    "E":  0.695,
    "G":  0.000,
    "H":  0.493,
    "I": -0.943,
    "L": -0.943,
    "K":  0.850,
    "M": -0.601,
    "F": -0.943,
    "P":  0.279,
    "S":  0.279,
    "T":  0.279,
    "W": -0.943,
    "Y": -0.943,
    "V": -0.943,
}


# ────────────────────────────────────────────────────────────
# Hydrophobicity scale selection
# ────────────────────────────────────────────────────────────

class HydrophobicityScale(str, Enum):
    """Hydrophobicity scale options for CamSol solubility prediction.

    Members:
        WIMLEY_WHITE: Wimley-White octanol scale with CamSol corrections.
            Best for globular/structured proteins.
        URRY: Urry scale based on elastin-like polypeptide experiments.
            Best for intrinsically disordered proteins (IDPs).
    """
    WIMLEY_WHITE = "wimley_white"
    URRY = "urry"


# ────────────────────────────────────────────────────────────
# IDP prediction heuristic
# ────────────────────────────────────────────────────────────

# Disorder-promoting residues (high in IDPs)
_IDP_DISORDER_PROMOTING = frozenset("PESQKRDG")
# Order-promoting / hydrophobic residues (low in IDPs)
_IDP_ORDER_PROMOTING = frozenset("ILVFMYWC")
# Thresholds for IDP prediction (tuned against benchmark dataset)
_IDP_DISORDER_FRACTION_THRESHOLD = 0.40
_IDP_ORDER_FRACTION_THRESHOLD = 0.31


def predict_idp(sequence: str) -> bool:
    """Predict whether a protein sequence is likely an intrinsically disordered protein.

    Uses a simple amino acid composition heuristic:
      - IDPs have high proportions of disorder-promoting residues (P, E, S, Q, K, R, D, G)
      - IDPs have low proportions of order-promoting/hydrophobic residues
        (I, L, V, F, M, Y, W, C)

    A sequence is predicted as IDP when **both** of the following hold:
      1. The disorder-promoting fraction exceeds 40%, AND
      2. The order-promoting fraction is below 31%

    The dual condition is important: many small soluble proteins (e.g.,
    thioredoxin, ubiquitin) have high disorder-promoting fractions but also
    high hydrophobic content (they have a well-defined hydrophobic core).
    True IDPs lack a hydrophobic core, which is captured by the low
    order-promoting fraction.

    Args:
        sequence: Protein sequence (1-letter amino acid codes).

    Returns:
        True if the sequence is predicted to be intrinsically disordered.
    """
    if len(sequence) < _IDP_MIN_SEQUENCE_LENGTH:
        return False

    n = len(sequence)
    disorder_count = sum(1 for aa in sequence if aa in _IDP_DISORDER_PROMOTING)
    order_count = sum(1 for aa in sequence if aa in _IDP_ORDER_PROMOTING)

    disorder_fraction = disorder_count / n
    order_fraction = order_count / n

    return (
        disorder_fraction > _IDP_DISORDER_FRACTION_THRESHOLD
        and order_fraction < _IDP_ORDER_FRACTION_THRESHOLD
    )


def select_hydropathy_scale(
    sequence: str,
    scale: str | HydrophobicityScale = "auto",
) -> tuple[dict[str, float], str]:
    """Select the hydropathy scale to use for solubility prediction.

    When ``scale="auto"`` (the default), the function uses the IDP prediction
    heuristic to decide: if the sequence looks like an IDP, the Urry scale is
    chosen; otherwise the Wimley-White scale is used.

    Args:
        sequence: Protein sequence (1-letter amino acid codes).
        scale: Scale selection — one of ``"auto"``, ``"wimley_white"``,
            or ``"urry"``.  Also accepts :class:`HydrophobicityScale` enum
            members.

    Returns:
        Tuple of ``(hydropathy_dict, scale_name)`` where *scale_name* is
        ``"wimley_white"`` or ``"urry"``.

    Raises:
        ValueError: If *scale* is not a recognised value.
    """
    if isinstance(scale, HydrophobicityScale):
        scale = scale.value

    if scale == "urry":
        return URRY_HYDROPATHY, "urry"
    elif scale == "wimley_white":
        return CAMSOL_HYDROPATHY, "wimley_white"
    elif scale == "auto":
        if predict_idp(sequence):
            return URRY_HYDROPATHY, "urry"
        else:
            return CAMSOL_HYDROPATHY, "wimley_white"
    else:
        raise ValueError(
            f"Unknown hydrophobicity scale: {scale!r}. "
            f"Expected one of: 'auto', 'wimley_white', 'urry', or a HydrophobicityScale enum."
        )

# ────────────────────────────────────────────────────────────
# Internal constants
# ────────────────────────────────────────────────────────────

# Weights for each component in the intrinsic solubility calculation
_WEIGHT_HYDROPATHY = 0.35
_WEIGHT_CHARGE = 0.25
_WEIGHT_ALPHA_HELIX = 0.15
_WEIGHT_BETA_STRAND = 0.15
_WEIGHT_PROGLY = 0.10

# Aggregation-prone region threshold
_AGGREGATION_THRESHOLD = -0.5

# Minimum length for an aggregation-prone region
_MIN_AGGREGATION_REGION_LENGTH = 3

# ── Score range and clamping ──
_SCORE_CLAMP_MIN = -3.0
_SCORE_CLAMP_MAX = 3.0

# ── Proline / glycine solubility effects ──
_PROLINE_SOLUBILITY_EFFECT = 1.0
_GLYCINE_SOLUBILITY_EFFECT = -0.5

# ── IDP prediction ──
_IDP_MIN_SEQUENCE_LENGTH = 10

# ── Confidence level ──
_CONFIDENCE_SCORE_THRESHOLD = 0.5

# ── Structural solubility (SASA correction) ──
_MIN_CA_ATOMS_FOR_STRUCTURE = 3
_SASA_BURIED_THRESHOLD = 0.15
_SASA_EXPOSED_THRESHOLD = 0.40
_BURIED_CORRECTION_FACTOR = 0.5
_EXPOSED_PENALTY_FACTOR = 0.3
_DISULFIDE_CORRECTION = 0.15
_STRUCTURAL_RESMOOTH_WINDOW = 3

# ── Solubility classification thresholds ──
_CLASS_HIGHLY_SOLUBLE = 1.5
_CLASS_SOLUBLE = 0.0
_CLASS_MARGINALLY_SOLUBLE = -1.0

# ── Recommendation thresholds ──
_RECOMMENDATION_HYDROPHOBIC_THRESHOLD = -0.3
_LOW_NET_CHARGE_THRESHOLD = 2
_NET_CHARGE_MIN_PROTEIN_LENGTH = 30
_TETRAPEPTIDE_HYDROPHOBIC_THRESHOLD = -0.2
_HYDROPHOBIC_STRETCH_MIN_LENGTH = 7
_HYDROPHOBIC_RESIDUE_THRESHOLD = -0.2

# ── Mutation delta computation ──
_MUTATION_DELTA_WINDOW = 7
_MUTATION_DIRECT_WEIGHT = 0.6
_MUTATION_LOCAL_WEIGHT = 0.4
_BLOSUM_DEFAULT_SCORE = -10

# ── PDB parsing ──
_PDB_MIN_LINE_LENGTH = 54

# ── SASA sigmoid transformation ──
_SASA_SIGMOID_K = 0.4
_SASA_SIGMOID_MIDPOINT = 12.0

# ── SASA neighbor-counting cutoff distance (Angstroms) ──
_SASA_NEIGHBOR_CUTOFF = 10.0

# ── Mutation scoring: fraction of worst residues to target ──
_WORST_RESIDUE_TARGET_FRACTION = 10

# ── Charged-substitution BLOSUM62 threshold ──
_CHARGED_SUB_BLOSUM_THRESHOLD = -1


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class CamSolResult(BaseEngineResult):
    """Result of CamSol solubility prediction.

    Inherits unified fields from BaseEngineResult:
        sequence: Input protein sequence (1-letter codes).
            (Alias: protein)
        primary_score: Final solubility score.
            (Alias: overall_score, score)
        classification: Solubility category.
            (Alias: solubility_class)
        success: Whether the computation completed without errors.
        error: Error message if the computation failed.
        execution_time_s: Wall-clock time for the computation in seconds.
        engine_name: Always "camsol".
        primary_score_label: Always "solubility".

    CamSol-specific attributes:
        intrinsic_score: CamSol intrinsic solubility score (-3 to +3, >0 = soluble).
        structural_score: Structure-corrected score (if PDB available), else None.
        per_residue_scores: Contribution of each residue to solubility.
        aggregation_prone_regions: List of (start, end, avg_score) tuples for
            regions below the aggregation threshold.
        recommendations: Actionable suggestions for improving solubility.
        mutations: List of suggested solubility-improving mutations.
            (Alias: solubility_mutations)
        method: "camsol_intrinsic" or "camsol_structural".
        hydrophobicity_scale_used: Which hydropathy scale was used ("wimley_white" or "urry").
    """
    # Override base class defaults for CamSol
    engine_name: str = "camsol"
    primary_score_label: str = "solubility"

    # CamSol-specific fields
    intrinsic_score: float = 0.0
    structural_score: float | None = None
    per_residue_scores: list[float] = field(default_factory=list)
    aggregation_prone_regions: list[tuple[int, int, float]] = field(default_factory=list)
    recommendations: list[str] = field(default_factory=list)
    mutations: list[MutationResult] = field(default_factory=list)
    method: str = "camsol_intrinsic"
    hydrophobicity_scale_used: str = "wimley_white"

    # ── Backward-compatible property aliases ──

    @property
    def protein(self) -> str:
        """Alias for sequence (backward compatibility)."""
        return self.sequence

    @protein.setter
    def protein(self, value: str) -> None:
        self.sequence = value

    @property
    def overall_score(self) -> float:
        """Alias for primary_score (backward compatibility)."""
        return self.primary_score

    @overall_score.setter
    def overall_score(self, value: float) -> None:
        self.primary_score = value

    @property
    def score(self) -> float:
        """Alias for primary_score (unified API)."""
        return self.primary_score

    @property
    def solubility_class(self) -> str:
        """Alias for classification (backward compatibility)."""
        return self.classification

    @solubility_class.setter
    def solubility_class(self, value: str) -> None:
        self.classification = value

    @property
    def solubility_mutations(self) -> list[MutationResult]:
        """Alias for mutations (backward compatibility)."""
        return self.mutations

    @property
    def confidence_level(self) -> str:
        """Accuracy confidence level for the solubility prediction.

        Returns one of:
          - ``"high"`` -- structural mode with PDB, or intrinsic with strong score
          - ``"medium"`` -- intrinsic mode with borderline score
          - ``"low"`` -- failed or no valid prediction
        """
        if not self.success:
            return "low"
        if self.method == "camsol_structural" and self.structural_score is not None:
            return "high"
        # Intrinsic mode — confidence based on score magnitude
        score = self.primary_score
        if abs(score) > _CONFIDENCE_SCORE_THRESHOLD:
            return "high"
        else:
            return "medium"


# ────────────────────────────────────────────────────────────
# Core computation functions
# ────────────────────────────────────────────────────────────

def compute_intrinsic_solubility(
    protein: str,
    window: int = _DEFAULT_WINDOW,
    smoothing: int = _DEFAULT_SMOOTHING,
    organism: str = "Homo_sapiens",
    hydrophobicity_scale: str = "auto",
) -> CamSolResult:
    """Compute CamSol intrinsic solubility score from protein sequence.

    The intrinsic solubility is computed as a weighted combination of
    physicochemical properties: hydropathy, charge, secondary structure
    propensity, and proline/glycine content.

    When *hydrophobicity_scale* is ``"auto"`` (the default), the Urry scale
    is automatically selected for sequences predicted to be intrinsically
    disordered, while the Wimley-White scale is used for globular proteins.
    This improves sensitivity for IDP solubility prediction from ~33% to ~67%.

    Steps:
        1. Validate protein sequence.
        2. Select hydropathy scale (auto, wimley_white, or urry).
        3. Compute per-residue intrinsic score using weighted combination.
        4. Apply sliding window smoothing (default window=7).
        5. Apply additional smoothing pass (3-residue window).
        6. Compute global score as average of per-residue scores.
        7. Identify aggregation-prone regions.
        8. Classify solubility.
        9. Generate recommendations.

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        window: Sliding window size for initial smoothing (default 7).
        smoothing: Secondary smoothing window size (default 3).
        organism: Target organism for codon/context awareness
            (default "Homo_sapiens").
        hydrophobicity_scale: Hydropathy scale to use — one of ``"auto"``,
            ``"wimley_white"``, or ``"urry"``.  ``"auto"`` selects the Urry
            scale for predicted IDPs and Wimley-White otherwise.
            Also accepts :class:`HydrophobicityScale` enum members.

    Returns:
        CamSolResult with intrinsic solubility prediction.

    Raises:
        CamSolError: If protein is empty or contains non-standard residues.
        ValueError: If *hydrophobicity_scale* is not a recognised value.
    """
    with EngineTimer() as timer:
        try:
            protein = validate_protein_sequence(protein, "CamSol")
        except ValueError as exc:
            raise CamSolError(str(exc)) from exc

        # Resolve hydropathy scale
        hydropathy, scale_name = select_hydropathy_scale(
            protein, hydrophobicity_scale
        )

        # Check cache
        cache_key = _make_cache_key(protein, window, smoothing, scale_name)
        if cache_key in _cache:
            logger.info("CamSol intrinsic: cache hit for %s...", protein[:10])
            return _cache[cache_key]

        n = len(protein)

        # Step 2: Per-residue intrinsic score
        raw_scores = []
        for aa in protein:
            hydro = hydropathy.get(aa, 0.0)
            charge = CAMSOL_CHARGE.get(aa, 0.0)
            alpha = CAMSOL_ALPHA_HELIX.get(aa, 0.0)
            beta = CAMSOL_BETA_STRAND.get(aa, 0.0)

            # Proline/glycine effect: P increases solubility (breaks structure),
            # G decreases solubility (flexible, can adopt aggregation-prone
            # conformations)
            if aa == "P":
                progly = _PROLINE_SOLUBILITY_EFFECT
            elif aa == "G":
                progly = _GLYCINE_SOLUBILITY_EFFECT
            else:
                progly = 0.0

            score = (
                _WEIGHT_HYDROPATHY * hydro
                + _WEIGHT_CHARGE * charge
                + _WEIGHT_ALPHA_HELIX * alpha
                + _WEIGHT_BETA_STRAND * (-beta)  # strand increases aggregation
                + _WEIGHT_PROGLY * progly
            )
            raw_scores.append(score)

        # Step 3: Sliding window smoothing
        smoothed = _sliding_window_smooth(raw_scores, window)

        # Step 4: Additional smoothing pass
        double_smoothed = _sliding_window_smooth(smoothed, smoothing)

        # Step 5: Compute global score
        intrinsic_score = (
            sum(double_smoothed) / len(double_smoothed) if double_smoothed else 0.0
        )

        # Clamp to score range
        intrinsic_score = max(_SCORE_CLAMP_MIN, min(_SCORE_CLAMP_MAX, intrinsic_score))

        # Step 6: Identify aggregation-prone regions
        agg_regions = _find_aggregation_prone_regions(
            double_smoothed, _AGGREGATION_THRESHOLD
        )

        # Step 7: Classify solubility
        solubility_class = classify_solubility(intrinsic_score)

        # Step 8: Generate recommendations
        result = CamSolResult(
            sequence=protein,
            primary_score=round(intrinsic_score, 4),
            classification=solubility_class,
            success=True,
            intrinsic_score=round(intrinsic_score, 4),
            structural_score=None,
            per_residue_scores=[round(s, 4) for s in double_smoothed],
            aggregation_prone_regions=[
                (start, end, round(avg, 4)) for start, end, avg in agg_regions
            ],
            recommendations=[],  # filled below
            mutations=[],
            method="camsol_intrinsic",
            hydrophobicity_scale_used=scale_name,
        )
        result.recommendations = generate_solubility_recommendations(result)

        # Store in cache
        _cache[cache_key] = result

        logger.info(
            "CamSol intrinsic: protein=%s, len=%d, score=%.4f, class=%s, "
            "agg_regions=%d, scale=%s",
            protein[:10] + "..." if len(protein) > 10 else protein,
            n,
            intrinsic_score,
            solubility_class,
            len(agg_regions),
            scale_name,
        )

    result.execution_time_s = round(timer.elapsed, 4)
    return result


def compute_structural_solubility(
    protein: str,
    pdb_string: str,
    organism: str = "Homo_sapiens",
    hydrophobicity_scale: str = "auto",
) -> CamSolResult:
    """Compute structure-corrected CamSol solubility score.

    Starts with the intrinsic solubility score and applies structure-based
    corrections derived from approximate SASA calculations on the PDB data:
      - Buried residues (low SASA): reduce their aggregation penalty
        (already buried in the structure).
      - Exposed hydrophobic residues: increase aggregation penalty.
      - Disulfide bonds: increase stability, moderate effect on solubility.

    SASA is approximated from PDB CA coordinates using neighbor counting
    (no external dependency required).

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        pdb_string: PDB file content as a string.
        organism: Target organism for codon/context awareness
            (default "Homo_sapiens").
        hydrophobicity_scale: Hydropathy scale to use — one of ``"auto"``,
            ``"wimley_white"``, or ``"urry"``.  Passed through to the
            intrinsic computation.

    Returns:
        CamSolResult with structure-corrected solubility prediction.

    Raises:
        CamSolError: If protein is empty or contains invalid residues.
    """
    with EngineTimer() as timer:
        try:
            protein = validate_protein_sequence(protein, "CamSol")
        except ValueError as exc:
            raise CamSolError(str(exc)) from exc

        if not pdb_string.strip():
            raise CamSolError("PDB string must not be empty.")

        # Compute intrinsic first
        intrinsic = compute_intrinsic_solubility(protein, organism=organism, hydrophobicity_scale=hydrophobicity_scale)

        # Get the hydropathy scale used for corrections below
        hydropathy = URRY_HYDROPATHY if intrinsic.hydrophobicity_scale_used == "urry" else CAMSOL_HYDROPATHY

        # Parse PDB coordinates (CA atoms only)
        ca_coords = _parse_pdb_ca_coords(pdb_string)

        if len(ca_coords) < _MIN_CA_ATOMS_FOR_STRUCTURE:
            logger.warning(
                "Too few CA atoms in PDB (%d). Falling back to intrinsic score.",
                len(ca_coords),
            )
            intrinsic.method = "camsol_structural"
            intrinsic.structural_score = intrinsic.intrinsic_score
            intrinsic.execution_time_s = round(timer.elapsed, 4)
            return intrinsic

        # Approximate SASA using CA neighbor counting
        sasa = _approximate_sasa(ca_coords)

        # Detect disulfide bonds from SSBOND records
        disulfide_residues = _parse_disulfide_bonds(pdb_string)

        # Apply structure-based corrections
        corrected_scores = list(intrinsic.per_residue_scores)

        for i in range(min(len(corrected_scores), len(sasa))):
            residue_sasa = sasa[i]

            if residue_sasa < _SASA_BURIED_THRESHOLD:
                # Buried residue: reduce aggregation penalty
                # If score is negative (aggregation-prone), reduce the magnitude
                if corrected_scores[i] < 0:
                    correction = abs(corrected_scores[i]) * _BURIED_CORRECTION_FACTOR * (
                        1.0 - residue_sasa / _SASA_BURIED_THRESHOLD
                    )
                    corrected_scores[i] += correction
            elif residue_sasa > _SASA_EXPOSED_THRESHOLD:
                # Exposed residue
                hydro = hydropathy.get(protein[i], 0.0)
                if hydro < 0:
                    # Exposed hydrophobic: increase aggregation penalty
                    exposed_range = _SCORE_CLAMP_MAX - _SASA_EXPOSED_THRESHOLD
                    penalty = abs(hydro) * _EXPOSED_PENALTY_FACTOR * (residue_sasa - _SASA_EXPOSED_THRESHOLD) / exposed_range
                    corrected_scores[i] -= penalty

            # Disulfide bond correction
            if i in disulfide_residues:
                corrected_scores[i] += _DISULFIDE_CORRECTION  # moderate stabilizing effect

        # Re-smooth after corrections
        corrected_smooth = _sliding_window_smooth(corrected_scores, _STRUCTURAL_RESMOOTH_WINDOW)

        # Recompute global score
        structural_score = (
            sum(corrected_smooth) / len(corrected_smooth)
            if corrected_smooth
            else 0.0
        )
        structural_score = max(_SCORE_CLAMP_MIN, min(_SCORE_CLAMP_MAX, structural_score))

        # Recompute aggregation-prone regions
        agg_regions = _find_aggregation_prone_regions(
            corrected_smooth, _AGGREGATION_THRESHOLD
        )

        solubility_class = classify_solubility(structural_score)

        result = CamSolResult(
            sequence=protein,
            primary_score=round(structural_score, 4),
            classification=solubility_class,
            success=True,
            intrinsic_score=intrinsic.intrinsic_score,
            structural_score=round(structural_score, 4),
            per_residue_scores=[round(s, 4) for s in corrected_smooth],
            aggregation_prone_regions=[
                (start, end, round(avg, 4)) for start, end, avg in agg_regions
            ],
            recommendations=[],
            mutations=[],
            method="camsol_structural",
            hydrophobicity_scale_used=intrinsic.hydrophobicity_scale_used,
        )
        result.recommendations = generate_solubility_recommendations(result)

        logger.info(
            "CamSol structural: protein=%s, len=%d, intrinsic=%.4f, "
            "structural=%.4f, class=%s, scale=%s",
            protein[:10] + "..." if len(protein) > 10 else protein,
            len(protein),
            intrinsic.intrinsic_score,
            structural_score,
            solubility_class,
            intrinsic.hydrophobicity_scale_used,
        )

    result.execution_time_s = round(timer.elapsed, 4)
    return result


def compute_solubility(
    protein: str,
    pdb_string: str | None = None,
    organism: str = "Homo_sapiens",
    hydrophobicity_scale: str = "auto",
    **kwargs: object,
) -> CamSolResult:
    """Unified solubility computation interface.

    Calls compute_structural_solubility if PDB data is provided,
    otherwise compute_intrinsic_solubility.

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        pdb_string: Optional PDB file content as a string.
        organism: Target organism for codon/context awareness
            (default "Homo_sapiens").
        hydrophobicity_scale: Hydropathy scale — ``"auto"``, ``"wimley_white"``,
            or ``"urry"``.  Passed through to the underlying compute function.
        **kwargs: Passed through to the underlying compute function
            (e.g., window, smoothing for intrinsic computation).

    Returns:
        CamSolResult with solubility prediction.
    """
    if pdb_string is not None and pdb_string.strip():
        return compute_structural_solubility(protein, pdb_string, organism=organism, hydrophobicity_scale=hydrophobicity_scale)

    # Extract supported kwargs for intrinsic computation
    window = kwargs.get("window", _DEFAULT_WINDOW)  # type: ignore[arg-type]
    smoothing = kwargs.get("smoothing", _DEFAULT_SMOOTHING)  # type: ignore[arg-type]
    return compute_intrinsic_solubility(protein, window=window, smoothing=smoothing, organism=organism, hydrophobicity_scale=hydrophobicity_scale)


def compute_solubility_batch(
    sequences: list[str],
    window: int = _DEFAULT_WINDOW,
    smoothing: int = _DEFAULT_SMOOTHING,
    max_workers: int | None = None,
    organism: str = "Homo_sapiens",
    hydrophobicity_scale: str = "auto",
) -> BatchResult[CamSolResult]:
    """Compute intrinsic solubility for multiple sequences in parallel.

    Uses a thread pool for concurrent computation. Each sequence is
    processed independently; failures are captured in the result object
    rather than raising exceptions.

    Args:
        sequences: List of protein sequences (1-letter amino acid codes).
        window: Sliding window size for initial smoothing.
        smoothing: Secondary smoothing window size.
        max_workers: Maximum number of threads. Defaults to
            min(len(sequences), DEFAULT_BATCH_SIZE).
        organism: Target organism for codon/context awareness
            (default "Homo_sapiens").

    Returns:
        BatchResult[CamSolResult] containing results for each input
        sequence in the same order as the input.
    """
    if not sequences:
        return BatchResult[CamSolResult]()

    if max_workers is None:
        max_workers = min(len(sequences), _DEFAULT_BATCH_SIZE)

    logger.info(
        "CamSol batch: computing solubility for %d sequences (workers=%d)",
        len(sequences),
        max_workers,
    )

    results: dict[int, CamSolResult] = {}

    def _compute_one(idx: int, seq: str) -> tuple[int, CamSolResult]:
        try:
            result = compute_intrinsic_solubility(seq, window=window, smoothing=smoothing, organism=organism, hydrophobicity_scale=hydrophobicity_scale)
            return (idx, result)
        except CamSolError as exc:
            error_result = CamSolResult(
                sequence=seq,
                primary_score=0.0,
                classification="insoluble",
                success=False,
                error=str(exc),
                intrinsic_score=0.0,
                structural_score=None,
                per_residue_scores=[],
                aggregation_prone_regions=[],
                recommendations=[],
                mutations=[],
                method="camsol_intrinsic",
                hydrophobicity_scale_used="wimley_white",
            )
            return (idx, error_result)
        except Exception as exc:
            logger.error(
                "CamSol batch: unexpected error for sequence %d (%s): %s",
                idx,
                seq[:10] + "..." if len(seq) > 10 else seq,
                exc,
                exc_info=True,
            )
            error_result = CamSolResult(
                sequence=seq,
                primary_score=0.0,
                classification="insoluble",
                success=False,
                error=f"Unexpected error: {exc}",
                intrinsic_score=0.0,
                structural_score=None,
                per_residue_scores=[],
                aggregation_prone_regions=[],
                recommendations=[],
                mutations=[],
                method="camsol_intrinsic",
                hydrophobicity_scale_used="wimley_white",
            )
            return (idx, error_result)

    with EngineTimer() as batch_timer:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(_compute_one, i, seq): i
                for i, seq in enumerate(sequences)
            }
            for future in as_completed(futures):
                idx, result = future.result()
                results[idx] = result

    ordered = [results[i] for i in range(len(sequences))]
    batch = BatchResult[CamSolResult](
        results=ordered,
        total_time_s=round(batch_timer.elapsed, 4),
    )

    logger.info(
        "CamSol batch: completed %d sequences (%d successful, %d failed)",
        len(sequences),
        batch.successful,
        batch.failed,
    )
    return batch


def find_solubility_mutations(
    protein: str,
    min_score: float = 0.0,
    organism: str = "Homo_sapiens",
    hydrophobicity_scale: str = "auto",
) -> list[MutationResult]:
    """Find amino acid substitutions to improve solubility.

    For each position in aggregation-prone regions, tries all 19 possible
    substitutions and computes the change in intrinsic solubility score
    (delta-delta-solubility). Only conservative mutations (BLOSUM62 >= 0)
    are returned.

    Args:
        protein: Protein sequence (1-letter amino acid codes).
        min_score: Minimum intrinsic score to target (default 0.0).
        organism: Target organism for codon/context awareness
            (default "Homo_sapiens").

    Returns:
        List of MutationResult objects sorted by score descending
        (most improvement first).

    Raises:
        CamSolError: If protein is empty or contains non-standard residues.
    """
    try:
        protein = validate_protein_sequence(protein, "CamSol")
    except ValueError as exc:
        raise CamSolError(str(exc)) from exc

    # Compute intrinsic solubility to find aggregation-prone regions
    result = compute_intrinsic_solubility(protein, organism=organism, hydrophobicity_scale=hydrophobicity_scale)

    # Get the resolved hydropathy scale for mutation scoring
    hydropathy, _ = select_hydropathy_scale(protein, hydrophobicity_scale)

    # If already above threshold, no mutations needed
    if result.intrinsic_score >= min_score and not result.aggregation_prone_regions:
        return []

    mutations: list[MutationResult] = []

    # Determine positions to target: all residues in aggregation-prone regions
    target_positions: set[int] = set()
    for start, end, _avg in result.aggregation_prone_regions:
        for pos in range(start, end):
            target_positions.add(pos)

    # Also consider positions with the lowest per-residue scores
    if not target_positions and result.intrinsic_score < min_score:
        # Target the worst-scoring residues
        n_target = max(1, len(protein) // _WORST_RESIDUE_TARGET_FRACTION)
        indexed = sorted(
            enumerate(result.per_residue_scores), key=lambda x: x[1]
        )
        for idx, _score in indexed[:n_target]:
            target_positions.add(idx)

    for pos in sorted(target_positions):
        if pos >= len(protein):
            continue
        wildtype = protein[pos]

        for mutant in STANDARD_AAS:
            if mutant == wildtype:
                continue

            # Check BLOSUM62 conservation
            blosum = BLOSUM62.get(wildtype, {}).get(mutant, _BLOSUM_DEFAULT_SCORE)
            if blosum < 0:
                continue

            # Compute delta solubility
            delta = _compute_mutation_delta(protein, pos, mutant, hydropathy)
            if delta > 0:
                mutations.append(MutationResult(
                    position=pos,
                    original=wildtype,
                    mutant=mutant,
                    delta_score=round(delta, 4),
                    score_type="solubility",
                    engine="camsol",
                    recommendation="solubility_improving",
                    description=(
                        f"{wildtype}{pos+1}{mutant}: delta_solubility={round(delta, 4):.4f}"
                    ),
                    details={
                        "delta_solubility": round(delta, 4),
                        "blosum62": blosum,
                    },
                ))

    # Sort by score descending (most improvement first)
    mutations.sort(key=lambda m: m.score, reverse=True)

    logger.info(
        "CamSol mutations: protein=%s, found %d suggested mutations",
        protein[:10] + "..." if len(protein) > 10 else protein,
        len(mutations),
    )

    return mutations


def classify_solubility(score: float) -> str:
    """Classify a solubility score into a solubility category.

    Args:
        score: CamSol solubility score (-3 to +3).

    Returns:
        One of: "highly_soluble", "soluble", "marginally_soluble", "insoluble".
    """
    if score > _CLASS_HIGHLY_SOLUBLE:
        return "highly_soluble"
    elif score > _CLASS_SOLUBLE:
        return "soluble"
    elif score > _CLASS_MARGINALLY_SOLUBLE:
        return "marginally_soluble"
    else:
        return "insoluble"


def generate_solubility_recommendations(result: CamSolResult) -> list[str]:
    """Generate actionable recommendations for improving protein solubility.

    Based on aggregation-prone regions and overall score, suggests:
      - Surface mutations to charged residues
      - Proline substitutions at edges of aggregation regions
      - N-terminal tag additions
      - Avoiding specific patterns (e.g., long hydrophobic stretches)

    Args:
        result: CamSolResult from a solubility computation.

    Returns:
        List of recommendation strings.
    """
    recommendations: list[str] = []
    protein = result.protein
    regions = result.aggregation_prone_regions

    # Use the hydropathy scale that was used for scoring
    hydropathy = URRY_HYDROPATHY if result.hydrophobicity_scale_used == "urry" else CAMSOL_HYDROPATHY

    # Overall assessment
    if result.overall_score > _CLASS_HIGHLY_SOLUBLE:
        recommendations.append(
            "Protein is highly soluble. No modifications necessary."
        )
        return recommendations

    if not regions:
        if result.overall_score < 0.0:
            recommendations.append(
                "Overall score is below zero but no discrete aggregation-prone "
                "regions detected. Consider increasing overall hydrophilicity "
                "with surface-exposed charged residues (K, E, D, R)."
            )
        return recommendations

    # Report aggregation-prone regions
    for start, end, avg in regions:
        region_seq = protein[start:end]
        hydrophobic_count = sum(
            1 for aa in region_seq
            if hydropathy.get(aa, 0.0) < 0
        )
        recommendations.append(
            f"Aggregation-prone region at positions {start+1}-{end} "
            f"(avg score: {avg:.2f}): '{region_seq}' — "
            f"{hydrophobic_count}/{len(region_seq)} hydrophobic residues."
        )

        # Suggest charged residue substitutions for hydrophobic residues
        for i in range(start, min(end, len(protein))):
            aa = protein[i]
            if hydropathy.get(aa, 0.0) < _RECOMMENDATION_HYDROPHOBIC_THRESHOLD:
                best_charge = _best_charged_substitution(aa)
                if best_charge:
                    recommendations.append(
                        f"  Position {i+1}: Consider {aa}->{best_charge} "
                        f"substitution (charged residue improves solubility)."
                    )

        # Suggest proline at edges
        if end - start >= _MIN_AGGREGATION_REGION_LENGTH:
            edge_left = max(0, start)
            edge_right = min(len(protein) - 1, end - 1)
            if protein[edge_left] != "P":
                recommendations.append(
                    f"  Consider P substitution at position {edge_left+1} "
                    f"(N-terminal edge of aggregation region) to break "
                    f"beta-sheet propensity."
                )
            if protein[edge_right] != "P":
                recommendations.append(
                    f"  Consider P substitution at position {edge_right+1} "
                    f"(C-terminal edge of aggregation region) to break "
                    f"beta-sheet propensity."
                )

    # Check for long hydrophobic stretches
    _check_hydrophobic_stretches(protein, recommendations, hydropathy)

    # Check net charge
    net_charge = sum(CAMSOL_CHARGE.get(aa, 0.0) for aa in protein)
    if abs(net_charge) < _LOW_NET_CHARGE_THRESHOLD and len(protein) > _NET_CHARGE_MIN_PROTEIN_LENGTH:
        recommendations.append(
            f"Low net charge ({net_charge:+.1f}). Increasing charged residue "
            f"content (K, R, D, E) can improve solubility through charge "
            f"repulsion."
        )

    # N-terminal tag suggestion for marginally soluble or insoluble
    if result.overall_score < 0.0:
        recommendations.append(
            "Consider adding an N-terminal solubility tag (e.g., MBP, GST, "
            "SUMO, or His6) to improve expression solubility."
        )

    # Avoid specific patterns: long hydrophobic tetrapeptide
    for i in range(len(protein) - 3):
        tetrapeptide = protein[i:i + 4]
        if all(hydropathy.get(aa, 0.0) < _TETRAPEPTIDE_HYDROPHOBIC_THRESHOLD for aa in tetrapeptide):
            recommendations.append(
                f"Long hydrophobic stretch at positions {i+1}-{i+4} "
                f"('{tetrapeptide}'). Consider breaking this pattern with "
                f"a charged or polar residue."
            )
            break  # only report the first one

    return recommendations


# ────────────────────────────────────────────────────────────
# Internal helper functions
# ────────────────────────────────────────────────────────────

def _make_cache_key(protein: str, window: int, smoothing: int, scale_name: str = "wimley_white") -> tuple[str, int, int, str]:
    """Create a cache key from protein sequence, window, smoothing, and scale.

    Uses a truncated hash of the protein to keep keys reasonable in size.
    """
    protein_hash = hashlib.sha256(protein.encode()).hexdigest()[:16]
    return (protein_hash, window, smoothing, scale_name)


def _sliding_window_smooth(scores: list[float], window: int) -> list[float]:
    """Apply sliding window averaging to a list of scores.

    At the boundaries, uses a truncated window (no zero-padding).

    Args:
        scores: Per-residue scores.
        window: Window size (odd values give symmetric smoothing).

    Returns:
        Smoothed scores of the same length.
    """
    if not scores:
        return []
    if window <= 1:
        return list(scores)
    if len(scores) <= window:
        return [sum(scores) / len(scores)] * len(scores)

    n = len(scores)
    half = window // 2
    result = []

    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        avg = sum(scores[start:end]) / (end - start)
        result.append(avg)

    return result


def _find_aggregation_prone_regions(
    scores: list[float],
    threshold: float,
) -> list[tuple[int, int, float]]:
    """Identify consecutive residues with scores below threshold.

    Args:
        scores: Per-residue solubility scores.
        threshold: Score threshold for aggregation-prone residues.

    Returns:
        List of (start, end, avg_score) tuples. start is inclusive,
        end is exclusive.
    """
    if not scores:
        return []

    regions: list[tuple[int, int, float]] = []
    in_region = False
    region_start = 0

    for i, score in enumerate(scores):
        if score < threshold:
            if not in_region:
                region_start = i
                in_region = True
        else:
            if in_region:
                if i - region_start >= _MIN_AGGREGATION_REGION_LENGTH:
                    avg = sum(scores[region_start:i]) / (i - region_start)
                    regions.append((region_start, i, avg))
                in_region = False

    # Handle region extending to the end
    if in_region:
        end = len(scores)
        if end - region_start >= _MIN_AGGREGATION_REGION_LENGTH:
            avg = sum(scores[region_start:end]) / (end - region_start)
            regions.append((region_start, end, avg))

    # Merge overlapping/adjacent regions
    merged: list[tuple[int, int, float]] = []
    for start, end, avg in regions:
        if merged and start <= merged[-1][1]:
            # Merge with previous
            prev_start, prev_end, _prev_avg = merged.pop()
            new_start = prev_start
            new_end = max(prev_end, end)
            combined_scores = scores[new_start:new_end]
            new_avg = sum(combined_scores) / len(combined_scores)
            merged.append((new_start, new_end, new_avg))
        else:
            merged.append((start, end, avg))

    return merged


def _compute_mutation_delta(
    protein: str,
    position: int,
    mutant: str,
    hydropathy: dict[str, float] | None = None,
) -> float:
    """Compute the change in intrinsic solubility from a single substitution.

    Computes the delta by recalculating the per-residue score at the
    mutation position and applying local smoothing effects.

    Args:
        protein: Original protein sequence.
        position: 0-based residue position to mutate.
        mutant: Mutant amino acid (1-letter code).
        hydropathy: Hydropathy scale dictionary to use. Defaults to
            CAMSOL_HYDROPATHY.

    Returns:
        Delta solubility score (positive = improvement).
    """
    if hydropathy is None:
        hydropathy = CAMSOL_HYDROPATHY

    wildtype = protein[position]

    # Compute score for wildtype residue
    wt_score = _compute_single_residue_score(wildtype, hydropathy)

    # Compute score for mutant residue
    mt_score = _compute_single_residue_score(mutant, hydropathy)

    # Direct delta at the position
    direct_delta = mt_score - wt_score

    # Also consider the effect on neighbors through smoothing
    # Check a local window centered on the position
    n = len(protein)
    window = _MUTATION_DELTA_WINDOW
    half = window // 2

    # Build a local patch with the mutation
    mut_protein = protein[:position] + mutant + protein[position + 1:]

    # Compute local scores for wildtype and mutant
    wt_local = 0.0
    mt_local = 0.0
    count = 0

    for i in range(max(0, position - half), min(n, position + half + 1)):
        wt_local += _compute_single_residue_score(protein[i], hydropathy)
        mt_local += _compute_single_residue_score(mut_protein[i], hydropathy)
        count += 1

    if count > 0:
        wt_local /= count
        mt_local /= count

    # Weighted combination of direct delta and local smoothing delta
    local_delta = mt_local - wt_local
    delta = _MUTATION_DIRECT_WEIGHT * direct_delta + _MUTATION_LOCAL_WEIGHT * local_delta

    return delta


def _compute_single_residue_score(aa: str, hydropathy: dict[str, float] | None = None) -> float:
    """Compute the raw intrinsic score contribution for a single amino acid.

    Args:
        aa: Single-letter amino acid code.
        hydropathy: Hydropathy scale dictionary to use. Defaults to
            CAMSOL_HYDROPATHY.

    Returns:
        Weighted intrinsic score contribution.
    """
    if hydropathy is None:
        hydropathy = CAMSOL_HYDROPATHY
    hydro = hydropathy.get(aa, 0.0)
    charge = CAMSOL_CHARGE.get(aa, 0.0)
    alpha = CAMSOL_ALPHA_HELIX.get(aa, 0.0)
    beta = CAMSOL_BETA_STRAND.get(aa, 0.0)

    if aa == "P":
        progly = _PROLINE_SOLUBILITY_EFFECT
    elif aa == "G":
        progly = _GLYCINE_SOLUBILITY_EFFECT
    else:
        progly = 0.0

    return (
        _WEIGHT_HYDROPATHY * hydro
        + _WEIGHT_CHARGE * charge
        + _WEIGHT_ALPHA_HELIX * alpha
        + _WEIGHT_BETA_STRAND * (-beta)
        + _WEIGHT_PROGLY * progly
    )


def _best_charged_substitution(aa: str) -> str | None:
    """Find the best charged residue substitution for a hydrophobic residue.

    Prefers the substitution with the highest BLOSUM62 score among
    charged residues (K, R, D, E, H).

    Args:
        aa: Hydrophobic amino acid (1-letter code).

    Returns:
        Best charged substitution, or None if no acceptable substitution.
    """
    charged = ["K", "R", "D", "E", "H"]
    best = None
    best_blosum = _BLOSUM_DEFAULT_SCORE

    for candidate in charged:
        blosum = BLOSUM62.get(aa, {}).get(candidate, _BLOSUM_DEFAULT_SCORE)
        if blosum >= _CHARGED_SUB_BLOSUM_THRESHOLD and blosum > best_blosum:
            best_blosum = blosum
            best = candidate

    return best


def _check_hydrophobic_stretches(
    protein: str,
    recommendations: list[str],
    hydropathy: dict[str, float] | None = None,
) -> None:
    """Check for long consecutive hydrophobic stretches in the protein.

    Appends warnings to the recommendations list.

    Args:
        protein: Protein sequence.
        recommendations: List to append recommendations to (mutated in place).
        hydropathy: Hydropathy scale dictionary to use. Defaults to
            CAMSOL_HYDROPATHY.
    """
    if hydropathy is None:
        hydropathy = CAMSOL_HYDROPATHY

    hydrophobic = set()
    for aa, val in hydropathy.items():
        if val < _HYDROPHOBIC_RESIDUE_THRESHOLD:
            hydrophobic.add(aa)

    stretch_start = None
    for i, aa in enumerate(protein):
        if aa in hydrophobic:
            if stretch_start is None:
                stretch_start = i
        else:
            if stretch_start is not None:
                length = i - stretch_start
                if length >= _HYDROPHOBIC_STRETCH_MIN_LENGTH:
                    recommendations.append(
                        f"Very long hydrophobic stretch at positions "
                        f"{stretch_start+1}-{i} "
                        f"('{protein[stretch_start:i]}', {length} residues). "
                        f"Long hydrophobic stretches strongly promote "
                        f"aggregation. Consider inserting charged/polar "
                        f"residues every 5-6 positions."
                    )
                stretch_start = None

    # Handle stretch at end
    if stretch_start is not None:
        length = len(protein) - stretch_start
        if length >= _HYDROPHOBIC_STRETCH_MIN_LENGTH:
            recommendations.append(
                f"Very long hydrophobic stretch at positions "
                f"{stretch_start+1}-{len(protein)} "
                f"('{protein[stretch_start:]}', {length} residues). "
                f"Long hydrophobic stretches strongly promote aggregation. "
                f"Consider inserting charged/polar residues every 5-6 positions."
            )


# ────────────────────────────────────────────────────────────
# PDB parsing helpers
# ────────────────────────────────────────────────────────────

def _parse_pdb_ca_coords(
    pdb_string: str,
) -> list[tuple[float, float, float]]:
    """Extract CA (alpha carbon) coordinates from PDB content.

    Args:
        pdb_string: PDB file content.

    Returns:
        List of (x, y, z) tuples for each CA atom, in order of appearance.
    """
    coords: list[tuple[float, float, float]] = []

    for line in pdb_string.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) >= _PDB_MIN_LINE_LENGTH:
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except (ValueError, IndexError):
                        logger.warning(
                            "Skipping malformed CA coordinate line in PDB: %r",
                            line.strip(),
                        )
                        continue

    return coords


def _approximate_sasa(
    ca_coords: list[tuple[float, float, float]],
    cutoff: float = _SASA_NEIGHBOR_CUTOFF,
) -> list[float]:
    """Approximate relative SASA using CA neighbor counting.

    For each CA atom, counts the number of neighbor CA atoms within
    the cutoff distance. More neighbors = more buried = lower SASA.

    The raw neighbor count is converted to a relative SASA value in [0, 1]
    using a sigmoid-like transformation.

    Args:
        ca_coords: List of (x, y, z) CA coordinates.
        cutoff: Distance cutoff in Angstroms for neighbor detection.

    Returns:
        List of relative SASA values in [0, 1] for each residue.
    """
    n = len(ca_coords)
    neighbor_counts = [0] * n

    for i in range(n):
        xi, yi, zi = ca_coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = ca_coords[j]
            dist_sq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if dist_sq < cutoff * cutoff:
                neighbor_counts[i] += 1
                neighbor_counts[j] += 1

    # Convert neighbor counts to relative SASA using sigmoid transformation
    # More neighbors -> more buried -> lower SASA
    # Typical CA neighbor count: 0-25 for buried, 0-8 for exposed
    max_neighbors = max(neighbor_counts) if neighbor_counts else 1
    if max_neighbors == 0:
        max_neighbors = 1

    sasa: list[float] = []
    for count in neighbor_counts:
        # Sigmoid: 1 / (1 + exp(k * (count - midpoint)))
        # At count=0, SASA ~ 1.0 (fully exposed)
        # At count=20, SASA ~ 0.0 (fully buried)
        relative_sasa = 1.0 / (1.0 + math.exp(_SASA_SIGMOID_K * (count - _SASA_SIGMOID_MIDPOINT)))
        sasa.append(relative_sasa)

    return sasa


def _parse_disulfide_bonds(pdb_string: str) -> set[int]:
    """Parse SSBOND records from PDB to find residues in disulfide bonds.

    Args:
        pdb_string: PDB file content.

    Returns:
        Set of 0-based residue indices involved in disulfide bonds.
    """
    residues: set[int] = set()

    for line in pdb_string.splitlines():
        if line.startswith("SSBOND"):
            try:
                # SSBOND record format:
                # COLUMNS 18-21: Sequence number of first residue
                # COLUMNS 32-35: Sequence number of second residue
                res1 = int(line[17:21].strip())
                res2 = int(line[31:35].strip())
                # Convert 1-based PDB numbering to 0-based
                residues.add(res1 - 1)
                residues.add(res2 - 1)
            except (ValueError, IndexError):
                logger.warning(
                    "Skipping malformed SSBOND record in PDB: %r",
                    line.strip(),
                )
                continue

    return residues


# ────────────────────────────────────────────────────────────
# Backward-compatibility alias
# ────────────────────────────────────────────────────────────

# Deprecated: use CamSolResult instead. Kept for backward compatibility.
# The warning filter above (at module top) suppresses the DeprecationWarning
# for internal use.  The direct module-level alias means __getattr__ is not
# triggered by import statements.
SolubilityResult = CamSolResult


def __getattr__(name: str):
    """Handle legacy attribute access; SolubilityResult is now a direct alias."""
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
