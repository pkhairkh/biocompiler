"""
BioCompiler Protein Design & Engineering Helpers v1.0.0
========================================================
Ties together sequence optimization and mutagenesis capabilities
for goal-directed protein engineering: thermostability, solubility,
deimmunization, and multi-objective design.

Design philosophy:
  - Each design function iteratively proposes mutations, verifies them
    against user-specified constraints (stability, solubility, immunogenicity,
    BLOSUM62 conservation), and stops when the target is reached or the
    mutation budget is exhausted.
  - Heavy dependencies (optimization, mutagenesis, species CAI tables) are
    lazily imported to keep the import graph lightweight.
  - All predictions are heuristic estimates — not substitutes for
    experimental validation.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, TypedDict

from biocompiler.shared.constants import BLOSUM62, HYDROPATHY, HYDROPHOBIC_AAS
from .base import BaseEngineResult, EngineTimer, validate_protein_sequence

logger = logging.getLogger(__name__)

__all__ = [
    "DesignResult",
    "DesignConstraints",
    "MutationScore",
    "score_mutation",
    "find_disulfide_opportunities",
    "find_proline_substitution_sites",
    "design_thermostable",
    "design_soluble",
    "design_low_immunogenicity",
    "design_multi_objective",
    # Re-exported constants for convenience
    "STABILITY_BLOSUM_WEIGHT",
    "STABILITY_HYDROPATHY_WEIGHT",
    "STABILITY_HYDROPATHY_COEFF",
    "STABILITY_PROLINE_COEFF",
    "STABILITY_GLYCINE_COEFF",
    "STABILITY_DISULFIDE_COEFF",
    "SOLUBILITY_HYDROPATHY_WEIGHT",
    "SOLUBILITY_CHARGED_WEIGHT",
    "SOLUBILITY_AGGREGATION_WEIGHT",
    "SOLUBILITY_HYDROPATHY_COEFF",
    "SOLUBILITY_CHARGED_COEFF",
    "SOLUBILITY_STRETCH_PENALTY",
    "MHC_II_WINDOW",
    "IMMUNOGENICITY_SCALE",
    "IMMUNOGENICITY_NORM_OFFSET",
    "IMMUNOGENICITY_NORM_RANGE",
    "HYDROPHOBIC_STRETCH_THRESHOLD",
    "SOFT_CONSTRAINT_TOLERANCE",
    "DISULFIDE_OPPORTUNITY_LIMIT",
    "BLOSUM62_MISSING_SCORE",
    "RESIDUE_CA_DISTANCE_ANGSTROMS",
    "MAX_EFFECTIVE_DISTANCE_ANGSTROMS",
    "DISULFIDE_MIN_RESIDUE_DISTANCE",
    "DISULFIDE_SHORT_LOOP_MAX_DISTANCE",
    "DISULFIDE_BASE_STABILIZATION",
    "DISULFIDE_DISTANCE_COEFF",
    "DISULFIDE_REFERENCE_DISTANCE",
    "DISULFIDE_MIN_STABILIZATION",
    "DISULFIDE_SHORT_LOOP_BONUS",
    "DEFAULT_WEIGHT_STABILITY",
    "DEFAULT_WEIGHT_SOLUBILITY",
    "DEFAULT_WEIGHT_IMMUNOGENICITY",
    "MULTI_OBJECTIVE_SOLUBILITY_TOLERANCE",
    "_base_immunogenicity",
    "_estimate_ddg",
    "_base_solubility",
    "_base_stability",
    "_check_constraints",
    "_estimate_immunogenicity_delta",
    "_estimate_solubility_delta",
    "_is_preserved",
    "_predict_secondary_structure_simple",
]

# ────────────────────────────────────────────────────────────
# Named constants for heuristic weights & thresholds
# ────────────────────────────────────────────────────────────

# Stability (ΔΔG) estimation weights
STABILITY_BLOSUM_WEIGHT: float = -0.15
STABILITY_HYDROPATHY_WEIGHT: float = 0.05
PROLINE_BONUS: float = -0.3
GLYCINE_BONUS: float = 0.3

# Base stability estimation coefficients
STABILITY_HYDROPATHY_COEFF: float = 0.5
STABILITY_PROLINE_COEFF: float = 2.0
STABILITY_GLYCINE_COEFF: float = 1.5
STABILITY_DISULFIDE_COEFF: float = 2.0

# Solubility estimation weights
SOLUBILITY_HYDROPATHY_WEIGHT: float = 0.2
SOLUBILITY_CHARGED_WEIGHT: float = 0.3
SOLUBILITY_AGGREGATION_WEIGHT: float = 0.2

# Base solubility estimation coefficients
SOLUBILITY_HYDROPATHY_COEFF: float = 0.3
SOLUBILITY_CHARGED_COEFF: float = 1.5
SOLUBILITY_STRETCH_PENALTY: float = 0.1

# Immunogenicity estimation
MHC_II_WINDOW: int = 9
IMMUNOGENICITY_SCALE: float = 0.05

# Base immunogenicity normalization constants
IMMUNOGENICITY_NORM_OFFSET: float = 4.5
IMMUNOGENICITY_NORM_RANGE: float = 9.0

# Hydrophobic stretch detection
HYDROPHOBIC_STRETCH_THRESHOLD: int = 5

# Soft constraint tolerance (allows gradual improvement)
SOFT_CONSTRAINT_TOLERANCE: float = 0.3

# Disulfide bond search limit
DISULFIDE_OPPORTUNITY_LIMIT: int = 20

# BLOSUM62 missing score fallback (used when a substitution is not in the matrix)
BLOSUM62_MISSING_SCORE: int = -4

# Disulfide bond geometric & stabilization model constants
RESIDUE_CA_DISTANCE_ANGSTROMS: float = 3.5   # Å per residue in extended chain
MAX_EFFECTIVE_DISTANCE_ANGSTROMS: float = 15.0  # Cap for effective Cα-Cα distance (Å)
DISULFIDE_MIN_RESIDUE_DISTANCE: int = 5        # Minimum sequence separation for disulfide
DISULFIDE_SHORT_LOOP_MAX_DISTANCE: int = 15    # Max distance for "short loop" bonus
DISULFIDE_BASE_STABILIZATION: float = -2.0     # Base stabilization (kcal/mol) for disulfide bond
DISULFIDE_DISTANCE_COEFF: float = 0.02         # Δstabilization per residue beyond reference distance
DISULFIDE_REFERENCE_DISTANCE: int = 10         # Reference distance for stabilization model
DISULFIDE_MIN_STABILIZATION: float = -3.0      # Floor for stabilization estimate (kcal/mol)
DISULFIDE_SHORT_LOOP_BONUS: float = -0.5       # Extra stabilization for short loops (kcal/mol)

# Default weights for multi-objective scoring
DEFAULT_WEIGHT_STABILITY: float = 0.4
DEFAULT_WEIGHT_SOLUBILITY: float = 0.3
DEFAULT_WEIGHT_IMMUNOGENICITY: float = 0.3

# Multi-objective soft solubility tolerance (wider than SOFT_CONSTRAINT_TOLERANCE)
MULTI_OBJECTIVE_SOLUBILITY_TOLERANCE: float = 0.5

# ────────────────────────────────────────────────────────────
# Convenience sets for mutation strategies
# ────────────────────────────────────────────────────────────

# HYDROPHOBIC_AAS imported from constants (Kyte-Doolittle > 1.0): {A, C, F, I, L, M, V}
# _AGGREGATION_PRONE_AAS is a strict subset of HYDROPHOBIC_AAS, limited to residues
# most strongly associated with β-aggregation (I, V, L, F, Y, W).  It excludes
# alanine, cysteine, and methionine which are hydrophobic but less aggregation-prone.
_CHARGED_AAS = set("DEKR")
_POLAR_AAS = set("STNQ")
_SURFACE_FAVORED_AAS = set("DEKRQN")   # charged + polar — good on surface
_AGGREGATION_PRONE_AAS = set("IVLFYW") # subset of HYDROPHOBIC_AAS most prone to aggregation

# Standard amino acid list (BLOSUM62 order)
_BLOSUM_SCORES = list("ARNDCQEGHILKMFPSTWYV")


# ────────────────────────────────────────────────────────────
# Data Classes
# ────────────────────────────────────────────────────────────

class MutationScore(TypedDict):
    """Structured return type for score_mutation().

    Keys:
        stability_ddg: Predicted ΔΔG for the mutation (kcal/mol, negative = stabilizing).
        solubility_delta: Predicted change in solubility score (positive = improved).
        immunogenicity_delta: Predicted change in immunogenicity (negative = deimmunized).
        blosum62: BLOSUM62 substitution score for wildtype → mutant.
        weighted_score: Combined weighted score across all dimensions.
    """
    stability_ddg: float
    solubility_delta: float
    immunogenicity_delta: float
    blosum62: int
    weighted_score: float


@dataclass
class DesignResult(BaseEngineResult):
    """Result of a protein design run.

    Inherits from BaseEngineResult for unified API compatibility.
    Domain-specific fields are preserved for backward compatibility.

    Unified field mapping:
      sequence → designed_protein
      primary_score → stability_change
      classification → 'design_success' | 'design_partial' | 'design_failed'
      engine_name → 'protein_design'
      primary_score_label → 'ddg'
    """

    # Override base class required fields with defaults for keyword-arg compat
    sequence: str = ""
    primary_score: float = 0.0
    classification: str = ""
    success: bool = False
    error: str | None = None
    execution_time_s: float = 0.0
    engine_name: str = "protein_design"
    primary_score_label: str = "ddg"

    # Domain-specific fields (backward-compatible)
    original_protein: str = ""
    designed_protein: str = ""
    mutations: list[dict[str, Any]] = field(default_factory=list)  # all mutations applied
    stability_change: float = 0.0           # ΔΔG (kcal/mol, negative = stabilizing)
    solubility_change: float = 0.0          # ΔCamSol score
    immunogenicity_change: float = 0.0      # Δimmunogenicity (negative = deimmunized)
    cai: float | None = None                # CAI of designed sequence (if DNA available)
    iterations: int = 0
    constraints_satisfied: list[str] = field(default_factory=list)  # names of satisfied constraints
    constraints_violated: list[str] = field(default_factory=list)   # names of violated constraints

    def __post_init__(self):
        # Sync unified base fields from domain-specific fields
        if not self.sequence and self.designed_protein:
            object.__setattr__(self, 'sequence', self.designed_protein)
        if self.primary_score == 0.0 and self.stability_change != 0.0:
            object.__setattr__(self, 'primary_score', self.stability_change)
        elif self.stability_change == 0.0 and self.primary_score != 0.0:
            object.__setattr__(self, 'stability_change', self.primary_score)
        if not self.classification:
            if not self.constraints_violated:
                label = "design_success"
            elif self.success:
                label = "design_partial"
            else:
                label = "design_failed"
            object.__setattr__(self, 'classification', label)

    @property
    def designed_sequence(self) -> str:
        """Unified API alias for designed_protein."""
        return self.designed_protein


@dataclass
class DesignConstraints:
    """User-specified constraints for protein design."""

    min_stability_kcal: float = -5.0
    min_solubility_score: float = 0.0
    max_immunogenicity: float = 0.5
    max_mutations: int = 10
    blosum62_min: int = 0
    max_ddg_per_mutation: float = 2.0
    preserve_positions: list[int] | None = None
    preserve_residues: list[str] | None = None  # e.g. ["C"] to preserve all cysteines


# ────────────────────────────────────────────────────────────
# Internal helpers — heuristic predictors
# ────────────────────────────────────────────────────────────

def _estimate_ddg(wildtype: str, mutant: str) -> float:
    """Heuristic ΔΔG estimate for a single-point mutation.

    Uses a simplified model:
      - BLOSUM62 score contributes: STABILITY_BLOSUM_WEIGHT * blosum
      - Hydropathy change contributes: STABILITY_HYDROPATHY_WEIGHT * Δhydropathy
      - Proline in a non-Pro context gets a stability bonus of PROLINE_BONUS
      - Glycine introduction gets a penalty of GLYCINE_BONUS

    Negative ΔΔG = stabilizing.
    This is a rough heuristic and should NOT be treated as a physics-based prediction.
    """
    blosum = BLOSUM62.get(wildtype, {}).get(mutant, BLOSUM62_MISSING_SCORE)
    dh = HYDROPATHY.get(mutant, 0.0) - HYDROPATHY.get(wildtype, 0.0)
    ddg = STABILITY_BLOSUM_WEIGHT * blosum + STABILITY_HYDROPATHY_WEIGHT * dh
    # Proline stabilization bonus (loop rigidification)
    if mutant == "P" and wildtype != "P":
        ddg += PROLINE_BONUS
    # Glycine flexibility penalty
    if mutant == "G" and wildtype != "G":
        ddg += GLYCINE_BONUS
    return round(ddg, 3)


def _estimate_solubility_delta(wildtype: str, mutant: str) -> float:
    """Heuristic change in CamSol-like solubility score for a mutation.

    Replacing hydrophobic surface residues with charged/polar improves solubility.
    Positive delta = improved solubility.
    """
    w_h = HYDROPATHY.get(wildtype, 0.0)
    m_h = HYDROPATHY.get(mutant, 0.0)
    delta = (w_h - m_h) * SOLUBILITY_HYDROPATHY_WEIGHT  # hydrophobic → hydrophilic = positive

    # Bonus for charged residues (D, E, K, R)
    if mutant in _CHARGED_AAS and wildtype not in _CHARGED_AAS:
        delta += SOLUBILITY_CHARGED_WEIGHT
    # Penalty for introducing aggregation-prone residue
    if mutant in _AGGREGATION_PRONE_AAS and wildtype not in _AGGREGATION_PRONE_AAS:
        delta -= SOLUBILITY_AGGREGATION_WEIGHT
    return round(delta, 3)


def _estimate_immunogenicity_delta(protein: str, position: int, mutant: str) -> float:
    """Heuristic change in immunogenicity for a mutation at *position*.

    Simplified model:
      - MHC-II binding is roughly correlated with hydrophobicity of the
        MHC_II_WINDOW-mer peptide centered on the position.
      - Replacing a hydrophobic anchor residue with polar/charged reduces
        predicted binding affinity.
    Negative delta = reduced immunogenicity.
    """
    half = MHC_II_WINDOW // 2
    start = max(0, position - half)
    end = min(len(protein), position + half + 1)

    old_hydro_sum = sum(HYDROPATHY.get(protein[i], 0.0) for i in range(start, end))
    new_protein = protein[:position] + mutant + protein[position + 1:]
    new_hydro_sum = sum(HYDROPATHY.get(new_protein[i], 0.0) for i in range(start, end))

    delta = (new_hydro_sum - old_hydro_sum) * IMMUNOGENICITY_SCALE
    return round(delta, 3)  # positive = more hydrophobic = more immunogenic


def _predict_secondary_structure_simple(protein: str) -> list[str]:
    """Very simple secondary-structure prediction (Chou-Fasman-inspired).

    Returns a list of structure codes, one per residue:
      'H' = helix, 'E' = strand, 'L' = loop/coil

    This is a heuristic placeholder — in production, use PSIPRED or similar.
    """
    HELIX_FORMERS = set("AELM")
    STRAND_FORMERS = set("VIY")
    MIN_HELIX = 4
    MIN_STRAND = 3

    n = len(protein)
    structure = ["L"] * n

    # Nucleation: find runs of helix/strand formers
    i = 0
    while i < n:
        # Try helix
        if protein[i] in HELIX_FORMERS:
            run = 0
            j = i
            while j < n and protein[j] in HELIX_FORMERS:
                run += 1
                j += 1
            if run >= MIN_HELIX:
                for k in range(i, j):
                    structure[k] = "H"
                i = j
                continue
        # Try strand
        if protein[i] in STRAND_FORMERS:
            run = 0
            j = i
            while j < n and protein[j] in STRAND_FORMERS:
                run += 1
                j += 1
            if run >= MIN_STRAND:
                for k in range(i, j):
                    structure[k] = "E"
                i = j
                continue
        i += 1

    return structure


def _is_preserved(position: int, wildtype: str, constraints: DesignConstraints) -> bool:
    """Check whether a position or residue type is protected from mutation."""
    if constraints.preserve_positions is not None and position in constraints.preserve_positions:
        return True
    if constraints.preserve_residues is not None and wildtype in constraints.preserve_residues:
        return True
    return False


def _check_constraints(
    protein: str,
    constraints: DesignConstraints,
    total_ddg: float,
    solubility: float,
    immunogenicity: float,
) -> tuple[list[str], list[str]]:
    """Return (satisfied, violated) constraint names for current state."""
    satisfied: list[str] = []
    violated: list[str] = []

    if total_ddg <= constraints.min_stability_kcal:
        satisfied.append("min_stability")
    else:
        violated.append("min_stability")

    if solubility >= constraints.min_solubility_score:
        satisfied.append("min_solubility")
    else:
        violated.append("min_solubility")

    if immunogenicity <= constraints.max_immunogenicity:
        satisfied.append("max_immunogenicity")
    else:
        violated.append("max_immunogenicity")

    return satisfied, violated


def _compute_cai_for_protein(protein: str, organism: str) -> float | None:
    """Try to compute CAI for a protein by optimizing it first. Returns None on failure."""
    try:
        from biocompiler.optimizer import optimize_sequence
        result = optimize_sequence(protein, organism=organism)
        return result.cai
    except Exception:
        return None


def _base_solubility(protein: str) -> float:
    """Compute a baseline CamSol-like intrinsic solubility score.

    Simplified model: average hydropathy, with penalties for long
    hydrophobic stretches and bonuses for charged content.
    """
    if not protein:
        return 0.0
    n = len(protein)
    avg_hydro = sum(HYDROPATHY.get(aa, 0.0) for aa in protein) / n
    charged_frac = sum(1 for aa in protein if aa in _CHARGED_AAS) / n

    # Penalty for long hydrophobic stretches (≥HYDROPHOBIC_STRETCH_THRESHOLD residues)
    stretch_penalty = 0.0
    run = 0
    for aa in protein:
        if aa in _AGGREGATION_PRONE_AAS:
            run += 1
            if run >= HYDROPHOBIC_STRETCH_THRESHOLD:
                stretch_penalty += SOLUBILITY_STRETCH_PENALTY
        else:
            run = 0

    score = -avg_hydro * SOLUBILITY_HYDROPATHY_COEFF + charged_frac * SOLUBILITY_CHARGED_COEFF - stretch_penalty
    return round(score, 3)


def _base_immunogenicity(protein: str) -> float:
    """Compute a baseline immunogenicity score (0-1 scale).

    Uses average hydrophobicity of MHC_II_WINDOW-mer windows as a proxy for MHC-II
    binding propensity. Higher = more immunogenic.
    """
    if len(protein) < MHC_II_WINDOW:
        return 0.0
    scores = []
    for i in range(len(protein) - MHC_II_WINDOW + 1):
        window = protein[i:i + MHC_II_WINDOW]
        avg_h = sum(HYDROPATHY.get(aa, 0.0) for aa in window) / MHC_II_WINDOW
        scores.append(avg_h)
    if not scores:
        return 0.0
    # Normalize: typical range is -IMMUNOGENICITY_NORM_OFFSET to +IMMUNOGENICITY_NORM_OFFSET; map to 0-1
    max_score = max(scores)
    normalized = (max_score + IMMUNOGENICITY_NORM_OFFSET) / IMMUNOGENICITY_NORM_RANGE
    return round(max(0.0, min(1.0, normalized)), 3)


def _base_stability(protein: str) -> float:
    """Compute a baseline stability estimate (ΔG in kcal/mol).

    Simplified: based on average hydrophobicity (core packing proxy),
    proline content (rigidification), and glycine content (flexibility).
    More negative = more stable.
    """
    if not protein:
        return 0.0
    n = len(protein)
    avg_hydro = sum(HYDROPATHY.get(aa, 0.0) for aa in protein) / n
    pro_frac = protein.count("P") / n
    gly_frac = protein.count("G") / n
    cys_count = protein.count("C")
    # Disulfide pairs contribute ~STABILITY_DISULFIDE_COEFF kcal/mol each
    disulfide_pairs = cys_count // 2
    ddg = -avg_hydro * STABILITY_HYDROPATHY_COEFF - pro_frac * STABILITY_PROLINE_COEFF + gly_frac * STABILITY_GLYCINE_COEFF - disulfide_pairs * STABILITY_DISULFIDE_COEFF
    return round(ddg, 3)


# ────────────────────────────────────────────────────────────
# Public API: Disulfide bond & proline substitution scanners
# ────────────────────────────────────────────────────────────

def find_disulfide_opportunities(protein: str) -> list[dict[str, Any]]:
    """Find positions where introducing cysteine pairs could form disulfide bonds.

    Criteria:
      - Both positions should be in predicted loop/coil regions (between
        secondary structure elements).
      - Positions must be at least DISULFIDE_MIN_RESIDUE_DISTANCE residues apart in sequence.
      - Neither position is already a cysteine.

    Returns:
        List of dicts with keys: position1, position2, distance_estimate,
        stabilizing_estimate.
    """
    ss = _predict_secondary_structure_simple(protein)
    n = len(protein)
    opportunities: list[dict[str, Any]] = []

    # Find loop positions (not already C)
    loop_positions = [i for i in range(n) if ss[i] == "L" and protein[i] != "C"]

    for idx_a in range(len(loop_positions)):
        for idx_b in range(idx_a + 1, len(loop_positions)):
            i = loop_positions[idx_a]
            j = loop_positions[idx_b]
            distance = abs(j - i)
            if distance < DISULFIDE_MIN_RESIDUE_DISTANCE:
                continue
            # Estimate Cα-Cα distance from sequence separation
            distance_estimate = distance * RESIDUE_CA_DISTANCE_ANGSTROMS
            # In a loop, the effective distance is much shorter than sequence distance
            effective_distance = min(distance_estimate, MAX_EFFECTIVE_DISTANCE_ANGSTROMS)

            # Stabilizing estimate: depends on how constrained the loop currently is
            # and the sequence distance (shorter loops with disulfides are more stabilizing)
            stabilizing_estimate = max(
                DISULFIDE_MIN_STABILIZATION,
                DISULFIDE_BASE_STABILIZATION + DISULFIDE_DISTANCE_COEFF * (distance - DISULFIDE_REFERENCE_DISTANCE),
            )
            # Short loops get more stabilization
            if distance <= DISULFIDE_SHORT_LOOP_MAX_DISTANCE:
                stabilizing_estimate += DISULFIDE_SHORT_LOOP_BONUS

            opportunities.append({
                "position1": i,
                "position2": j,
                "distance_estimate": round(effective_distance, 1),
                "stabilizing_estimate": round(stabilizing_estimate, 3),
            })

    # Sort by most stabilizing first
    opportunities.sort(key=lambda x: x["stabilizing_estimate"])
    return opportunities


def find_proline_substitution_sites(protein: str) -> list[dict[str, Any]]:
    """Find positions where proline substitution would stabilize (loop → rigid).

    Criteria:
      - Position is in a predicted loop/coil region.
      - Original residue is not already proline.
      - Original residue is not glycine (G→P is highly destabilizing).
      - BLOSUM62 score for the substitution is >= -1 (moderately conservative).

    Returns:
        List of dicts with keys: position, wildtype, ddg_estimate, in_loop.
    """
    ss = _predict_secondary_structure_simple(protein)
    n = len(protein)
    sites: list[dict[str, Any]] = []

    for i in range(n):
        if ss[i] != "L":
            continue
        wt = protein[i]
        if wt == "P" or wt == "G":
            continue
        blosum = BLOSUM62.get(wt, {}).get("P", BLOSUM62_MISSING_SCORE)
        if blosum < -1:
            continue
        ddg_est = _estimate_ddg(wt, "P")
        sites.append({
            "position": i,
            "wildtype": wt,
            "ddg_estimate": ddg_est,
            "in_loop": True,
        })

    # Sort by most stabilizing first
    sites.sort(key=lambda x: x["ddg_estimate"])
    return sites


# ────────────────────────────────────────────────────────────
# Public API: Mutation scoring
# ────────────────────────────────────────────────────────────

def score_mutation(
    protein: str,
    position: int,
    mutant: str,
    organism: str = "Homo_sapiens",
    weights: dict | None = None,
) -> MutationScore:
    """Score a single mutation across all dimensions.

    Args:
        protein: Original protein sequence (1-letter codes).
        position: 0-based residue position.
        mutant: Mutant amino acid (1-letter code).
        organism: Target organism for CAI lookup.
        weights: Optional weight dict with keys 'stability', 'solubility',
                 'immunogenicity'. Defaults to DEFAULT_WEIGHT_* constants.

    Returns:
        MutationScore with keys: stability_ddg, solubility_delta, immunogenicity_delta,
        blosum62, weighted_score.
    """
    if weights is None:
        weights = {
            "stability": DEFAULT_WEIGHT_STABILITY,
            "solubility": DEFAULT_WEIGHT_SOLUBILITY,
            "immunogenicity": DEFAULT_WEIGHT_IMMUNOGENICITY,
        }

    wildtype = protein[position]
    ddg = _estimate_ddg(wildtype, mutant)
    sol_delta = _estimate_solubility_delta(wildtype, mutant)
    imm_delta = _estimate_immunogenicity_delta(protein, position, mutant)
    blosum = BLOSUM62.get(wildtype, {}).get(mutant, BLOSUM62_MISSING_SCORE)

    # Normalize components for weighted score:
    #   stability: lower ddg is better → score = -ddg (positive = good)
    #   solubility: higher delta is better → score = delta
    #   immunogenicity: lower delta is better → score = -delta
    w_stab = weights.get("stability", DEFAULT_WEIGHT_STABILITY)
    w_sol = weights.get("solubility", DEFAULT_WEIGHT_SOLUBILITY)
    w_imm = weights.get("immunogenicity", DEFAULT_WEIGHT_IMMUNOGENICITY)

    weighted_score = (
        w_stab * (-ddg) +
        w_sol * sol_delta +
        w_imm * (-imm_delta)
    )

    return {
        "stability_ddg": ddg,
        "solubility_delta": sol_delta,
        "immunogenicity_delta": imm_delta,
        "blosum62": blosum,
        "weighted_score": round(weighted_score, 4),
    }


# ────────────────────────────────────────────────────────────
# Internal: Mutation proposal callbacks
# ────────────────────────────────────────────────────────────

def _propose_thermostable(
    current: list[str],
    current_protein: str,
    current_sol: float,
    current_imm: float,
    current_stab: float,
    cons: DesignConstraints,
) -> dict[str, Any] | None:
    """Propose the best stabilizing mutation using three strategies.

    Strategies (in order of consideration):
      1. Single-point mutations that lower ΔΔG.
      2. Proline substitutions in loop regions.
      3. Disulfide bond introduction (Cys pairs in loops).
    """
    best_mutation = None
    best_ddg = 0.0  # only accept stabilizing (negative ddg)

    # --- Strategy 1: Scan all single-point mutations for best stabilizing one ---
    for pos in range(len(current)):
        wt = current[pos]
        if _is_preserved(pos, wt, cons):
            continue
        for mutant_aa in _BLOSUM_SCORES:
            if mutant_aa == wt:
                continue
            blosum = BLOSUM62.get(wt, {}).get(mutant_aa, BLOSUM62_MISSING_SCORE)
            if blosum < cons.blosum62_min:
                continue
            ddg = _estimate_ddg(wt, mutant_aa)
            if ddg >= 0:
                continue  # not stabilizing
            if abs(ddg) > cons.max_ddg_per_mutation:
                continue  # too large a change per step
            # Soft constraint checks: do not make things significantly worse
            sol_delta = _estimate_solubility_delta(wt, mutant_aa)
            imm_delta = _estimate_immunogenicity_delta(
                current_protein, pos, mutant_aa
            )
            new_sol = current_sol + sol_delta
            new_imm = current_imm + imm_delta
            if new_sol < cons.min_solubility_score - SOFT_CONSTRAINT_TOLERANCE:
                continue
            if new_imm > cons.max_immunogenicity + SOFT_CONSTRAINT_TOLERANCE:
                continue
            if ddg < best_ddg:
                best_ddg = ddg
                best_mutation = {
                    "position": pos,
                    "wildtype": wt,
                    "mutant": mutant_aa,
                    "ddg": ddg,
                    "solubility_delta": sol_delta,
                    "immunogenicity_delta": imm_delta,
                    "blosum62": blosum,
                    "strategy": "single_point",
                }

    # --- Strategy 2: Proline substitutions in loops ---
    proline_sites = find_proline_substitution_sites(current_protein)
    for site in proline_sites:
        pos = site["position"]
        wt = current[pos]
        if _is_preserved(pos, wt, cons):
            continue
        ddg = site["ddg_estimate"]
        if ddg >= 0 or ddg >= best_ddg:
            continue
        blosum = BLOSUM62.get(wt, {}).get("P", BLOSUM62_MISSING_SCORE)
        if blosum < cons.blosum62_min:
            continue
        sol_delta = _estimate_solubility_delta(wt, "P")
        imm_delta = _estimate_immunogenicity_delta(
            current_protein, pos, "P"
        )
        new_sol = current_sol + sol_delta
        new_imm = current_imm + imm_delta
        if new_sol < cons.min_solubility_score - SOFT_CONSTRAINT_TOLERANCE:
            continue
        if new_imm > cons.max_immunogenicity + SOFT_CONSTRAINT_TOLERANCE:
            continue
        best_ddg = ddg
        best_mutation = {
            "position": pos,
            "wildtype": wt,
            "mutant": "P",
            "ddg": ddg,
            "solubility_delta": sol_delta,
            "immunogenicity_delta": imm_delta,
            "blosum62": blosum,
            "strategy": "proline_in_loop",
        }

    # --- Strategy 3: Disulfide bond introduction ---
    disulfide_ops = find_disulfide_opportunities(current_protein)
    for opp in disulfide_ops[:DISULFIDE_OPPORTUNITY_LIMIT]:
        pos1 = opp["position1"]
        pos2 = opp["position2"]
        wt1 = current[pos1]
        wt2 = current[pos2]
        if _is_preserved(pos1, wt1, cons):
            continue
        if _is_preserved(pos2, wt2, cons):
            continue
        ddg_pair = opp["stabilizing_estimate"]
        if ddg_pair >= best_ddg:
            continue
        # Check constraints for both mutations combined
        sol1 = _estimate_solubility_delta(wt1, "C")
        sol2 = _estimate_solubility_delta(wt2, "C")
        imm1 = _estimate_immunogenicity_delta(current_protein, pos1, "C")
        imm2 = _estimate_immunogenicity_delta(current_protein, pos2, "C")
        new_sol = current_sol + sol1 + sol2
        new_imm = current_imm + imm1 + imm2
        if new_sol < cons.min_solubility_score - SOFT_CONSTRAINT_TOLERANCE:
            continue
        if new_imm > cons.max_immunogenicity + SOFT_CONSTRAINT_TOLERANCE:
            continue
        best_ddg = ddg_pair
        best_mutation = {
            "position": pos1,
            "wildtype": wt1,
            "mutant": "C",
            "ddg": ddg_pair / 2,
            "solubility_delta": sol1,
            "immunogenicity_delta": imm1,
            "blosum62": BLOSUM62.get(wt1, {}).get("C", BLOSUM62_MISSING_SCORE),
            "strategy": "disulfide_pair",
            "pair_position": pos2,
            "pair_wildtype": wt2,
        }

    return best_mutation


def _apply_disulfide_partner(
    current: list[str],
    mutation: dict[str, Any],
    total_ddg: float,
) -> tuple[float, list[dict[str, Any]]]:
    """Apply the partner cysteine for a disulfide pair mutation."""
    if mutation.get("strategy") != "disulfide_pair":
        return total_ddg, []
    pair_pos = mutation["pair_position"]
    pair_wt = mutation["pair_wildtype"]
    current[pair_pos] = "C"
    total_ddg += mutation["ddg"]  # symmetric contribution
    partner = {
        "position": pair_pos,
        "wildtype": pair_wt,
        "mutant": "C",
        "ddg": mutation["ddg"],
        "solubility_delta": _estimate_solubility_delta(pair_wt, "C"),
        "immunogenicity_delta": _estimate_immunogenicity_delta(
            "".join(current), pair_pos, "C"
        ),
        "blosum62": BLOSUM62.get(pair_wt, {}).get("C", BLOSUM62_MISSING_SCORE),
        "strategy": "disulfide_pair_partner",
    }
    return total_ddg, [partner]


def _propose_soluble(
    current: list[str],
    current_protein: str,
    current_sol: float,
    current_imm: float,
    current_stab: float,
    cons: DesignConstraints,
) -> dict[str, Any] | None:
    """Propose the best solubility-improving mutation.

    Focuses on replacing hydrophobic/aggregation-prone surface residues
    with charged or polar residues.
    """
    best_mutation = None
    best_sol_delta = 0.0

    for pos in range(len(current)):
        wt = current[pos]
        if _is_preserved(pos, wt, cons):
            continue
        # Focus on hydrophobic / aggregation-prone residues
        if wt not in HYDROPHOBIC_AAS and wt not in _AGGREGATION_PRONE_AAS:
            continue
        for mutant_aa in _SURFACE_FAVORED_AAS | _POLAR_AAS:
            if mutant_aa == wt:
                continue
            blosum = BLOSUM62.get(wt, {}).get(mutant_aa, BLOSUM62_MISSING_SCORE)
            if blosum < cons.blosum62_min:
                continue
            ddg = _estimate_ddg(wt, mutant_aa)
            if abs(ddg) > cons.max_ddg_per_mutation:
                continue
            # Soft stability check: only block if making stability
            # significantly worse than current (not vs. threshold)
            if ddg > cons.max_ddg_per_mutation:
                continue  # strongly destabilizing
            sol_delta = _estimate_solubility_delta(wt, mutant_aa)
            if sol_delta <= best_sol_delta:
                continue
            imm_delta = _estimate_immunogenicity_delta(current_protein, pos, mutant_aa)
            new_imm = current_imm + imm_delta
            # Soft immunogenicity check: allow some tolerance
            if new_imm > cons.max_immunogenicity + SOFT_CONSTRAINT_TOLERANCE:
                continue
            best_sol_delta = sol_delta
            best_mutation = {
                "position": pos,
                "wildtype": wt,
                "mutant": mutant_aa,
                "ddg": ddg,
                "solubility_delta": sol_delta,
                "immunogenicity_delta": imm_delta,
                "blosum62": blosum,
                "strategy": "surface_hydrophilic",
            }

    return best_mutation


def _propose_deimmunize(
    current: list[str],
    current_protein: str,
    current_sol: float,
    current_imm: float,
    current_stab: float,
    cons: DesignConstraints,
) -> dict[str, Any] | None:
    """Propose the best immunogenicity-reducing mutation.

    Finds the most immunogenic MHC_II_WINDOW-mer window and proposes
    a substitution within it that reduces hydrophobicity.
    """
    # Find the most immunogenic MHC_II_WINDOW-mer window
    best_window_start = 0
    best_window_score = -999.0
    for i in range(len(current_protein) - MHC_II_WINDOW + 1):
        window = current_protein[i:i + MHC_II_WINDOW]
        score = sum(HYDROPATHY.get(aa, 0.0) for aa in window) / MHC_II_WINDOW
        if score > best_window_score:
            best_window_score = score
            best_window_start = i

    # Within that window, find the best mutation to reduce hydrophobicity
    best_mutation = None
    best_imm_delta = 0.0  # negative = reducing immunogenicity

    for pos in range(best_window_start, best_window_start + MHC_II_WINDOW):
        if pos >= len(current):
            break
        wt = current[pos]
        if _is_preserved(pos, wt, cons):
            continue
        for mutant_aa in _POLAR_AAS | _CHARGED_AAS | set("AST"):
            if mutant_aa == wt:
                continue
            blosum = BLOSUM62.get(wt, {}).get(mutant_aa, BLOSUM62_MISSING_SCORE)
            if blosum < cons.blosum62_min:
                continue
            ddg = _estimate_ddg(wt, mutant_aa)
            if abs(ddg) > cons.max_ddg_per_mutation:
                continue
            # Soft stability check: do not make things much worse
            if ddg > cons.max_ddg_per_mutation:
                continue
            imm_delta = _estimate_immunogenicity_delta(current_protein, pos, mutant_aa)
            if imm_delta >= best_imm_delta:
                continue  # want negative (reducing immunogenicity)
            sol_delta = _estimate_solubility_delta(wt, mutant_aa)
            new_sol = current_sol + sol_delta
            # Soft solubility check: allow tolerance
            if new_sol < cons.min_solubility_score - SOFT_CONSTRAINT_TOLERANCE:
                continue
            best_imm_delta = imm_delta
            best_mutation = {
                "position": pos,
                "wildtype": wt,
                "mutant": mutant_aa,
                "ddg": ddg,
                "solubility_delta": sol_delta,
                "immunogenicity_delta": imm_delta,
                "blosum62": blosum,
                "strategy": "deimmunize_window",
            }

    return best_mutation


def _propose_multi_objective(
    current: list[str],
    current_protein: str,
    current_sol: float,
    current_imm: float,
    current_stab: float,
    cons: DesignConstraints,
    weights: dict[str, float],
) -> dict[str, Any] | None:
    """Propose the best multi-objective mutation (weighted score improvement).

    Returns None if no mutation improves the weighted score (i.e., best_weighted <= 0).
    """
    best_mutation = None
    best_weighted = -999.0

    w_stab = weights.get("stability", DEFAULT_WEIGHT_STABILITY)
    w_sol = weights.get("solubility", DEFAULT_WEIGHT_SOLUBILITY)
    w_imm = weights.get("immunogenicity", DEFAULT_WEIGHT_IMMUNOGENICITY)

    for pos in range(len(current)):
        wt = current[pos]
        if _is_preserved(pos, wt, cons):
            continue
        for mutant_aa in _BLOSUM_SCORES:
            if mutant_aa == wt:
                continue
            blosum = BLOSUM62.get(wt, {}).get(mutant_aa, BLOSUM62_MISSING_SCORE)
            if blosum < cons.blosum62_min:
                continue
            ddg = _estimate_ddg(wt, mutant_aa)
            if abs(ddg) > cons.max_ddg_per_mutation:
                continue
            # Soft stability check: only block if mutation makes stability
            # significantly worse than current state (not vs. threshold)
            if ddg > cons.max_ddg_per_mutation:
                continue  # strongly destabilizing

            sol_delta = _estimate_solubility_delta(wt, mutant_aa)
            imm_delta = _estimate_immunogenicity_delta(current_protein, pos, mutant_aa)

            new_sol = current_sol + sol_delta
            new_imm = current_imm + imm_delta

            # Soft constraints: prefer improvements but do not hard-block
            # unless the resulting value is clearly worse than threshold
            if new_sol < cons.min_solubility_score - MULTI_OBJECTIVE_SOLUBILITY_TOLERANCE:
                continue
            if new_imm > cons.max_immunogenicity + SOFT_CONSTRAINT_TOLERANCE:
                continue

            # Compute weighted score
            weighted = (
                w_stab * (-ddg) +
                w_sol * sol_delta +
                w_imm * (-imm_delta)
            )

            if weighted > best_weighted:
                best_weighted = weighted
                best_mutation = {
                    "position": pos,
                    "wildtype": wt,
                    "mutant": mutant_aa,
                    "ddg": ddg,
                    "solubility_delta": sol_delta,
                    "immunogenicity_delta": imm_delta,
                    "blosum62": blosum,
                    "weighted_score": round(weighted, 4),
                    "strategy": "multi_objective",
                }

    # Only return if the best mutation actually improves things
    if best_weighted <= 0:
        return None

    return best_mutation


# ────────────────────────────────────────────────────────────
# Internal: Common iterative design loop
# ────────────────────────────────────────────────────────────

def _design_iterative(
    protein: str,
    organism: str,
    constraints: DesignConstraints,
    propose_mutation: Callable[
        [list[str], str, float, float, float, DesignConstraints],
        dict[str, Any] | None,
    ],
    is_objective_met: Callable[[float, float, float], bool],
    is_success: Callable[[float, float, float], bool] | None = None,
    apply_extra: Callable[
        [list[str], dict[str, Any], float],
        tuple[float, list[dict[str, Any]]],
    ] | None = None,
) -> DesignResult:
    """Common iterative design loop shared by design_thermostable/soluble/low_immunogenicity.

    Handles timer setup/teardown, state initialization, the iteration loop,
    final evaluation, and DesignResult construction.  Each caller provides
    its own ``propose_mutation`` callback that encapsulates the domain-specific
    mutation search strategy.

    Args:
        protein: Validated protein sequence.
        organism: Target organism for CAI lookup.
        constraints: Design constraints.
        propose_mutation: Callback that searches for the best mutation given
            current state.  Signature:
            ``(current_list, current_protein_str, current_sol, current_imm,
              current_stab, constraints) -> mutation dict | None``.
        is_objective_met: Returns True when the design objective is reached.
            Signature: ``(current_stab, current_sol, current_imm) -> bool``.
        is_success: Returns True if the final design is successful.
            Defaults to *is_objective_met* if not provided.
            Signature: ``(current_stab, current_sol, current_imm) -> bool``.
        apply_extra: Optional callback after applying the primary mutation, for
            additional logic (e.g., disulfide pair partner).  Signature:
            ``(current_list, mutation, total_ddg) -> (new_total_ddg, extra_mutations)``.
    """
    if is_success is None:
        is_success = is_objective_met

    _timer = EngineTimer()
    _timer.__enter__()

    current = list(protein)
    base_stab = _base_stability(protein)
    current_stab = base_stab
    total_ddg = 0.0
    mutations: list[dict[str, Any]] = []
    iterations = 0

    for iteration in range(constraints.max_mutations):
        iterations += 1
        current_protein = "".join(current)
        current_sol = _base_solubility(current_protein)
        current_imm = _base_immunogenicity(current_protein)

        if is_objective_met(current_stab, current_sol, current_imm):
            break

        best_mutation = propose_mutation(
            current, current_protein, current_sol, current_imm,
            current_stab, constraints,
        )

        if best_mutation is None:
            logger.info("No improving mutation found at iteration %d", iteration)
            break

        # Apply the best mutation
        pos = best_mutation["position"]
        current[pos] = best_mutation["mutant"]
        total_ddg += best_mutation["ddg"]
        current_stab = base_stab + total_ddg
        mutations.append(best_mutation)

        # Apply extra logic (e.g., disulfide pair partner)
        if apply_extra is not None:
            total_ddg, extra_mutations = apply_extra(current, best_mutation, total_ddg)
            current_stab = base_stab + total_ddg
            mutations.extend(extra_mutations)

    designed = "".join(current)
    final_sol = _base_solubility(designed)
    final_imm = _base_immunogenicity(designed)
    cai = _compute_cai_for_protein(designed, organism)

    satisfied, violated = _check_constraints(
        designed, constraints, current_stab, final_sol, final_imm,
    )

    _timer.__exit__(None, None, None)

    return DesignResult(
        original_protein=protein,
        designed_protein=designed,
        mutations=mutations,
        stability_change=total_ddg,
        solubility_change=final_sol - _base_solubility(protein),
        immunogenicity_change=final_imm - _base_immunogenicity(protein),
        cai=cai,
        iterations=iterations,
        constraints_satisfied=satisfied,
        constraints_violated=violated,
        success=is_success(current_stab, final_sol, final_imm),
        execution_time_s=round(_timer.elapsed, 4),
    )


# ────────────────────────────────────────────────────────────
# Public API: Design functions
# ────────────────────────────────────────────────────────────

def design_thermostable(
    protein: str,
    organism: str = "Homo_sapiens",
    target_stability: float = -10.0,
    constraints: DesignConstraints | None = None,
) -> DesignResult:
    """Find mutations that increase thermostability.

    Strategy:
      1. Search for stabilizing mutations (negative ΔΔG).
      2. Consider disulfide bond introduction (Cys pairs in loops).
      3. Consider proline substitutions in loops.
      4. Apply mutations iteratively until target_stability is reached
         or max_mutations is exhausted.
      5. Verify each mutation does not violate solubility/immunogenicity
         constraints.

    Args:
        protein: Input protein sequence (1-letter codes).
        organism: Target organism.
        target_stability: Target ΔG in kcal/mol (more negative = more stable).
        constraints: Design constraints.

    Returns:
        DesignResult with all mutations and metrics.
    """
    protein = validate_protein_sequence(protein, "protein_design")
    if constraints is None:
        constraints = DesignConstraints()

    return _design_iterative(
        protein=protein,
        organism=organism,
        constraints=constraints,
        propose_mutation=_propose_thermostable,
        is_objective_met=lambda stab, _sol, _imm: stab <= target_stability,
        apply_extra=_apply_disulfide_partner,
    )


def design_soluble(
    protein: str,
    organism: str = "Homo_sapiens",
    min_solubility: float = 0.5,
    constraints: DesignConstraints | None = None,
) -> DesignResult:
    """Find mutations that improve solubility.

    Strategy:
      1. Replace hydrophobic surface residues with charged/polar residues.
      2. Break aggregation-prone regions (long hydrophobic stretches).
      3. Apply mutations iteratively until min_solubility is reached
         or max_mutations is exhausted.
      4. Verify stability is preserved after each mutation.

    Args:
        protein: Input protein sequence (1-letter codes).
        organism: Target organism.
        min_solubility: Target solubility score.
        constraints: Design constraints.

    Returns:
        DesignResult with all mutations and metrics.
    """
    protein = validate_protein_sequence(protein, "protein_design")
    if constraints is None:
        constraints = DesignConstraints()

    return _design_iterative(
        protein=protein,
        organism=organism,
        constraints=constraints,
        propose_mutation=_propose_soluble,
        is_objective_met=lambda _stab, sol, _imm: sol >= min_solubility,
    )


def design_low_immunogenicity(
    protein: str,
    organism: str = "Homo_sapiens",
    max_immunogenicity: float = 0.3,
    constraints: DesignConstraints | None = None,
) -> DesignResult:
    """Reduce immunogenicity via amino acid substitution.

    Strategy:
      - Identify high-immunogenicity MHC_II_WINDOW-mer windows (hydrophobic-rich).
      - Propose substitutions that reduce window hydrophobicity while
        preserving BLOSUM62 conservation.
      - Apply iteratively until max_immunogenicity threshold is met
        or max_mutations exhausted.

    Args:
        protein: Input protein sequence.
        organism: Target organism.
        max_immunogenicity: Target immunogenicity score.
        constraints: Design constraints.

    Returns:
        DesignResult with all mutations and metrics.
    """
    protein = validate_protein_sequence(protein, "protein_design")
    if constraints is None:
        constraints = DesignConstraints()

    return _design_iterative(
        protein=protein,
        organism=organism,
        constraints=constraints,
        propose_mutation=_propose_deimmunize,
        is_objective_met=lambda _stab, _sol, imm: imm <= max_immunogenicity,
    )


def design_multi_objective(
    protein: str,
    organism: str = "Homo_sapiens",
    constraints: DesignConstraints | None = None,
    weights: dict | None = None,
) -> DesignResult:
    """Multi-objective optimization: balance stability, solubility, and immunogenicity.

    For each candidate mutation, compute a weighted score improvement.
    Apply the mutation with the best combined improvement.
    Iterate until all constraints are satisfied or max_mutations is reached.

    Default weights: {"stability": 0.4, "solubility": 0.3, "immunogenicity": 0.3}

    Args:
        protein: Input protein sequence.
        organism: Target organism.
        constraints: Design constraints.
        weights: Objective weights dict.

    Returns:
        DesignResult with all mutations and metrics.
    """
    protein = validate_protein_sequence(protein, "protein_design")
    if constraints is None:
        constraints = DesignConstraints()
    if weights is None:
        weights = {
            "stability": DEFAULT_WEIGHT_STABILITY,
            "solubility": DEFAULT_WEIGHT_SOLUBILITY,
            "immunogenicity": DEFAULT_WEIGHT_IMMUNOGENICITY,
        }

    def _propose_mo(
        current: list[str],
        current_protein: str,
        current_sol: float,
        current_imm: float,
        current_stab: float,
        cons: DesignConstraints,
    ) -> dict[str, Any] | None:
        return _propose_multi_objective(
            current, current_protein, current_sol, current_imm,
            current_stab, cons, weights,
        )

    def _all_constraints_met(stab: float, sol: float, imm: float) -> bool:
        return (
            stab <= constraints.min_stability_kcal
            and sol >= constraints.min_solubility_score
            and imm <= constraints.max_immunogenicity
        )

    return _design_iterative(
        protein=protein,
        organism=organism,
        constraints=constraints,
        propose_mutation=_propose_mo,
        is_objective_met=_all_constraints_met,
    )
