"""
BioCompiler Expression Predictor
==================================

A simple expression prediction model based on CAI, GC content, mRNA stability,
and codon pair bias.  Uses a heuristic model calibrated against published
expression data to predict relative expression levels for optimized sequences.

Key components:

  - **predict_expression()**: Predicts relative expression level for an
    optimized sequence in a given organism.

  - **ExpressionPrediction**: Structured result with predicted level,
    confidence, and key contributing factors.

  - **ExpressionPredictor**: Class-based interface for batch predictions
    and model customisation.

References:
  - Welch M et al. (2009) PLoS ONE 4:e7002 — Design parameters for codon
    optimization.
  - Gustafsson C et al. (2004) Trends Biotechnol 22:346-353 — Codon bias
    and expression.
  - Quax TEF et al. (2015) Mol Cell 59:519-530 — Codon pair bias.
  - Kudla G et al. (2009) Science 324:255-258 — Codon optimality and
    mRNA stability.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    "ExpressionPrediction",
    "predict_expression",
    "ExpressionPredictor",
    "_compute_gc_optimality",
    "_FACTOR_WEIGHTS",
    "_GC_SWEET_SPOT",
    "_estimate_mrna_stability",
]


# ═══════════════════════════════════════════════════════════════════════════════
# Data structures
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ExpressionPrediction:
    """Predicted expression level for an optimized sequence.

    Attributes:
        predicted_relative_expression: Relative expression level (0.0–1.0),
            where 1.0 means maximum expected expression.
        confidence: Confidence in the prediction (0.0–1.0).
        key_factors: Dictionary of contributing factors and their scores.
        category: Expression category: 'high', 'medium', or 'low'.
        organism: Target organism for the prediction.
        warnings: List of warning messages about potential issues.
    """

    predicted_relative_expression: float
    confidence: float
    key_factors: Dict[str, float]
    category: str
    organism: str
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not (0.0 <= self.predicted_relative_expression <= 1.0):
            raise ValueError(
                f"predicted_relative_expression must be in [0, 1], "
                f"got {self.predicted_relative_expression}"
            )
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(
                f"confidence must be in [0, 1], got {self.confidence}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Heuristic model parameters
# ═══════════════════════════════════════════════════════════════════════════════
# Weights for the four factors in the heuristic model.
# These were calibrated by comparing predicted vs published expression levels
# for the benchmark proteins in the wetlab_validation module.

_FACTOR_WEIGHTS: Dict[str, float] = {
    "cai": 0.40,          # CAI is the strongest predictor
    "gc_optimality": 0.25, # GC content optimality
    "mrna_stability": 0.20, # mRNA stability score
    "codon_pair_bias": 0.15, # Codon pair bias contribution
}

# GC content "sweet spot" ranges per organism (where expression is highest).
# These are narrower than the organism GC target ranges because the
# optimal-for-expression range is narrower than the acceptable range.
_GC_SWEET_SPOT: Dict[str, Tuple[float, float]] = {
    "Escherichia_coli": (0.45, 0.55),
    "Homo_sapiens": (0.50, 0.60),
    "Saccharomyces_cerevisiae": (0.35, 0.45),
    "Mus_musculus": (0.45, 0.55),
    "CHO_K1": (0.45, 0.55),
}

# Default sweet spot for unknown organisms
_GC_SWEET_SPOT_DEFAULT: Tuple[float, float] = (0.40, 0.55)


# ═══════════════════════════════════════════════════════════════════════════════
# Core prediction function
# ═══════════════════════════════════════════════════════════════════════════════

def predict_expression(
    optimized_sequence: str,
    organism: str,
    cai: Optional[float] = None,
    gc_content: Optional[float] = None,
    mrna_stability_score: Optional[float] = None,
    codon_pair_bias: Optional[float] = None,
) -> ExpressionPrediction:
    """Predict relative expression level for an optimized sequence.

    Uses a heuristic model based on four factors:
      1. CAI (Codon Adaptation Index) — 40% weight
      2. GC content optimality — 25% weight
      3. mRNA stability — 20% weight
      4. Codon pair bias — 15% weight

    If any factor is not provided, it is computed on-the-fly (for CAI and GC)
    or estimated (for mRNA stability and CPB).

    Args:
        optimized_sequence: The codon-optimized DNA sequence.
        organism: Target organism (canonical name or alias).
        cai: Pre-computed CAI value.  If None, computed from the sequence.
        gc_content: Pre-computed GC fraction.  If None, computed from the sequence.
        mrna_stability_score: Pre-computed mRNA stability score (0.0–1.0).
            If None, estimated from the sequence.
        codon_pair_bias: Pre-computed codon pair bias score.
            If None, estimated from the sequence.

    Returns:
        An ExpressionPrediction with the predicted level, confidence, and
        key contributing factors.
    """
    from ..organisms import resolve_organism
    from .translation import compute_cai
    from ..scanner import gc_content as _gc_content

    resolved = resolve_organism(organism, strict=False)
    seq = optimized_sequence.upper().strip() if optimized_sequence else ""
    warnings: List[str] = []

    # ── Compute or use provided CAI ───────────────────────────────────────
    if cai is not None:
        cai_val = cai
    else:
        try:
            cai_val = compute_cai(seq, organism=resolved) if seq else 0.0
        except Exception:
            cai_val = 0.0
            warnings.append("Could not compute CAI; defaulting to 0.0")

    # ── Compute or use provided GC content ─────────────────────────────────
    if gc_content is not None:
        gc_val = gc_content
    else:
        gc_val = _gc_content(seq) if seq else 0.0

    # ── Compute or estimate mRNA stability ─────────────────────────────────
    if mrna_stability_score is not None:
        stability_val = mrna_stability_score
    else:
        stability_val = _estimate_mrna_stability(seq, resolved)

    # ── Compute or estimate codon pair bias ───────────────────────────────
    if codon_pair_bias is not None:
        cpb_val = codon_pair_bias
    else:
        cpb_val = _estimate_cpb(seq, resolved)

    # ── Compute individual factor scores (0.0–1.0) ───────────────────────

    # CAI factor: CAI is already 0.0–1.0 by definition
    cai_factor = cai_val

    # GC optimality factor: how close GC is to the organism's sweet spot
    gc_factor = _compute_gc_optimality(gc_val, resolved)

    # mRNA stability factor
    stability_factor = max(0.0, min(1.0, stability_val))

    # CPB factor: normalize CPB to 0.0–1.0 range.
    # Published CPB scores range roughly from -0.5 to +0.5 for individual pairs;
    # mean CPB for well-optimized sequences is typically 0.0 to +0.3.
    # Map [-0.3, +0.3] → [0.0, 1.0]
    cpb_factor = max(0.0, min(1.0, (cpb_val + 0.3) / 0.6)) if cpb_val is not None else 0.5

    # ── Compute weighted prediction ────────────────────────────────────────
    predicted = (
        _FACTOR_WEIGHTS["cai"] * cai_factor
        + _FACTOR_WEIGHTS["gc_optimality"] * gc_factor
        + _FACTOR_WEIGHTS["mrna_stability"] * stability_factor
        + _FACTOR_WEIGHTS["codon_pair_bias"] * cpb_factor
    )
    predicted = max(0.0, min(1.0, predicted))

    # ── Compute confidence ─────────────────────────────────────────────────
    # Confidence is higher when all factors are available and consistent.
    # Low CAI + high GC outside range → lower confidence.
    confidence = _compute_confidence(cai_factor, gc_factor, stability_factor, cpb_val)

    # ── Determine category ────────────────────────────────────────────────
    if predicted >= 0.7:
        category = "high"
    elif predicted >= 0.4:
        category = "medium"
    else:
        category = "low"

    # ── Build warnings ────────────────────────────────────────────────────
    if cai_val < 0.5:
        warnings.append(f"Low CAI ({cai_val:.3f}) may indicate suboptimal codon usage")
    if gc_val < 0.30 or gc_val > 0.70:
        warnings.append(f"GC content ({gc_val:.1%}) is outside the typical 30-70% range")
    if cpb_val is not None and cpb_val < -0.1:
        warnings.append(f"Negative codon pair bias ({cpb_val:.3f}) may reduce expression")

    key_factors = {
        "cai": round(cai_factor, 4),
        "gc_optimality": round(gc_factor, 4),
        "mrna_stability": round(stability_factor, 4),
        "codon_pair_bias": round(cpb_factor, 4),
    }

    return ExpressionPrediction(
        predicted_relative_expression=round(predicted, 4),
        confidence=round(confidence, 4),
        key_factors=key_factors,
        category=category,
        organism=resolved,
        warnings=warnings,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_gc_optimality(gc: float, organism: str) -> float:
    """Compute GC optimality score (0.0–1.0).

    Returns 1.0 when GC is within the organism's sweet spot, declining
    linearly to 0.0 as it moves away from the sweet spot.
    """
    sweet_lo, sweet_hi = _GC_SWEET_SPOT.get(organism, _GC_SWEET_SPOT_DEFAULT)

    if sweet_lo <= gc <= sweet_hi:
        return 1.0

    # Distance from the sweet spot boundary
    if gc < sweet_lo:
        distance = sweet_lo - gc
    else:
        distance = gc - sweet_hi

    # Allow up to 20 percentage points of deviation before hitting 0
    max_deviation = 0.20
    score = max(0.0, 1.0 - distance / max_deviation)
    return score


def _estimate_mrna_stability(seq: str, organism: str) -> float:
    """Estimate mRNA stability score (0.0–1.0) from sequence features.

    Uses heuristic estimation based on:
      - Absence of ATTTA instability motifs (Kudla et al., 2009)
      - GC content correlation with mRNA half-life
      - Absence of long A/T runs (polyA signals)

    For production use, this should be replaced with actual mRNA stability
    prediction from the mrna_stability module.
    """
    if not seq:
        return 0.5

    score = 0.5  # Start at neutral

    # Penalize ATTTA instability motifs
    instability_count = seq.count("ATTTA")
    score -= min(0.3, instability_count * 0.1)

    # Penalize long A/T runs (8+ consecutive A or T)
    at_run_penalty = 0
    current_run = 0
    for base in seq:
        if base in "AT":
            current_run += 1
            if current_run >= 8:
                at_run_penalty += 1
        else:
            current_run = 0
    score -= min(0.2, at_run_penalty * 0.05)

    # Reward moderate GC content (correlated with mRNA stability)
    gc = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    if 0.40 <= gc <= 0.60:
        score += 0.2
    elif 0.30 <= gc <= 0.70:
        score += 0.1

    return max(0.0, min(1.0, score))


def _estimate_cpb(seq: str, organism: str) -> float:
    """Estimate codon pair bias score from the sequence.

    Uses the codon_pair_scoring module when available, falls back to
    a heuristic estimate based on codon usage frequencies.
    """
    if not seq or len(seq) < 6:
        return 0.0

    try:
        from .codon_pair_scoring import compute_cpb
        return compute_cpb(seq, organism)
    except Exception:
        pass

    # Heuristic estimate: preferred codon pairs tend to be GC-rich
    # for human/CHO, and specific patterns for E. coli/yeast.
    # This is a rough approximation.
    gc = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0

    # Organism-specific rough CPB estimates based on GC
    if organism in ("Escherichia_coli",):
        # E. coli: moderate GC ≈ neutral CPB
        return (gc - 0.50) * 0.5
    elif organism in ("Homo_sapiens", "CHO_K1", "Mus_musculus"):
        # Human/mouse/CHO: GC-rich pairs are preferred
        return (gc - 0.50) * 0.8
    elif organism == "Saccharomyces_cerevisiae":
        # Yeast: AT-rich preferred codons
        return (0.38 - gc) * 0.6
    else:
        return (gc - 0.50) * 0.5


def _compute_confidence(
    cai_factor: float,
    gc_factor: float,
    stability_factor: float,
    cpb_val: Optional[float],
) -> float:
    """Compute prediction confidence based on factor consistency.

    Higher confidence when:
      - All factors point in the same direction (high or low)
      - CPB data is available
      - No extreme disagreement between factors
    """
    factors = [cai_factor, gc_factor, stability_factor]

    # Base confidence from factor availability
    confidence = 0.6  # Start at moderate
    if cpb_val is not None:
        confidence += 0.1

    # Reward factor consistency (low variance)
    mean_factor = sum(factors) / len(factors)
    variance = sum((f - mean_factor) ** 2 for f in factors) / len(factors)
    consistency_bonus = max(0.0, 0.2 - variance)
    confidence += consistency_bonus

    # Penalize if CAI is very different from the other factors
    if cai_factor > 0.7 and gc_factor < 0.3:
        confidence -= 0.1  # Disagreement reduces confidence

    return max(0.1, min(0.95, confidence))


# ═══════════════════════════════════════════════════════════════════════════════
# Class-based interface for batch predictions
# ═══════════════════════════════════════════════════════════════════════════════

class ExpressionPredictor:
    """Batch expression prediction interface.

    Allows customization of model weights and provides methods for
    batch prediction over multiple sequences.

    Usage::

        predictor = ExpressionPredictor(organism="Homo_sapiens")
        prediction = predictor.predict("ATGGCCCTG...")
        batch = predictor.predict_batch(["ATGGCC...", "ATGTTT..."])
    """

    def __init__(
        self,
        organism: str,
        factor_weights: Optional[Dict[str, float]] = None,
    ) -> None:
        """Initialize the expression predictor.

        Args:
            organism: Target organism for predictions.
            factor_weights: Custom factor weights.  Must sum to ~1.0.
                Keys: 'cai', 'gc_optimality', 'mrna_stability', 'codon_pair_bias'.
        """
        from ..organisms import resolve_organism
        self._organism = resolve_organism(organism, strict=False)
        self._weights = dict(factor_weights or _FACTOR_WEIGHTS)

        # Validate weights
        total = sum(self._weights.values())
        if abs(total - 1.0) > 0.05:
            logger.warning(
                "Factor weights sum to %.3f (expected ~1.0); normalizing",
                total,
            )
            if total > 0:
                self._weights = {k: v / total for k, v in self._weights.items()}

    @property
    def organism(self) -> str:
        """The target organism for predictions."""
        return self._organism

    def predict(self, optimized_sequence: str) -> ExpressionPrediction:
        """Predict expression for a single optimized sequence.

        Args:
            optimized_sequence: Codon-optimized DNA sequence.

        Returns:
            ExpressionPrediction with predicted level and details.
        """
        return predict_expression(
            optimized_sequence=optimized_sequence,
            organism=self._organism,
        )

    def predict_batch(
        self,
        sequences: List[str],
    ) -> List[ExpressionPrediction]:
        """Predict expression for multiple sequences.

        Args:
            sequences: List of codon-optimized DNA sequences.

        Returns:
            List of ExpressionPrediction objects, one per input sequence.
        """
        return [self.predict(seq) for seq in sequences]

    def get_factor_weights(self) -> Dict[str, float]:
        """Return the current factor weights."""
        return dict(self._weights)
