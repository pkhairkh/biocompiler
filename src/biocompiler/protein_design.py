"""
BioCompiler Protein Design & Engineering Helpers v7.2.0
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

from .constants import BLOSUM62, HYDROPATHY

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Convenience sets for mutation strategies
# ────────────────────────────────────────────────────────────

_HYDROPHOBIC_AAS = set("AILMFWV")
_CHARGED_AAS = set("DEKR")
_POLAR_AAS = set("STNQ")
_SURFACE_FAVORED = set("DEKRQN")   # charged + polar — good on surface
_AGGREGATION_PRONE = set("IVLFYW") # hydrophobic stretches

# Standard amino acid index (BLOSUM62 order)
_BLOSUM_INDEX = list("ARNDCQEGHILKMFPSTWYV")


# ────────────────────────────────────────────────────────────
# Data Classes
# ────────────────────────────────────────────────────────────

@dataclass
class DesignResult:
    """Result of a protein design run."""

    original_protein: str
    designed_protein: str
    mutations: list[dict]             # all mutations applied
    stability_change: float           # ΔΔG (kcal/mol, negative = stabilizing)
    solubility_change: float          # ΔCamSol score
    immunogenicity_change: float      # Δimmunogenicity (negative = deimmunized)
    cai: float | None                 # CAI of designed sequence (if DNA available)
    iterations: int
    constraints_satisfied: list[str]  # names of satisfied constraints
    constraints_violated: list[str]   # names of violated constraints
    success: bool


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
      - BLOSUM62 score contributes: -0.15 * blosum (higher conservation → less destabilizing)
      - Hydropathy change contributes: -0.05 * Δhydropathy
      - Proline in a non-Pro context gets a stability bonus of -0.3 (rigidifies loops)
      - Glycine introduction gets a penalty of +0.3 (increases backbone flexibility)

    Negative ΔΔG = stabilizing.
    This is a rough heuristic and should NOT be treated as a physics-based prediction.
    """
    blosum = BLOSUM62.get(wildtype, {}).get(mutant, -4)
    dh = HYDROPATHY.get(mutant, 0.0) - HYDROPATHY.get(wildtype, 0.0)
    ddg = -0.15 * blosum + 0.05 * dh
    # Proline stabilization bonus (loop rigidification)
    if mutant == "P" and wildtype != "P":
        ddg -= 0.3
    # Glycine flexibility penalty
    if mutant == "G" and wildtype != "G":
        ddg += 0.3
    return round(ddg, 3)


def _estimate_solubility_delta(wildtype: str, mutant: str) -> float:
    """Heuristic change in CamSol-like solubility score for a mutation.

    Replacing hydrophobic surface residues with charged/polar improves solubility.
    Positive delta = improved solubility.
    """
    w_h = HYDROPATHY.get(wildtype, 0.0)
    m_h = HYDROPATHY.get(mutant, 0.0)
    delta = (w_h - m_h) * 0.2  # hydrophobic → hydrophilic = positive

    # Bonus for charged residues (D, E, K, R)
    if mutant in _CHARGED_AAS and wildtype not in _CHARGED_AAS:
        delta += 0.3
    # Penalty for introducing aggregation-prone residue
    if mutant in _AGGREGATION_PRONE and wildtype not in _AGGREGATION_PRONE:
        delta -= 0.2
    return round(delta, 3)


def _estimate_immunogenicity_delta(protein: str, position: int, mutant: str) -> float:
    """Heuristic change in immunogenicity for a mutation at *position*.

    Simplified model:
      - MHC-II binding is roughly correlated with hydrophobicity of the
        9-mer peptide centered on the position.
      - Replacing a hydrophobic anchor residue with polar/charged reduces
        predicted binding affinity.
    Negative delta = reduced immunogenicity.
    """
    window = 9
    half = window // 2
    start = max(0, position - half)
    end = min(len(protein), position + half + 1)

    old_hydro_sum = sum(HYDROPATHY.get(protein[i], 0.0) for i in range(start, end))
    new_protein = protein[:position] + mutant + protein[position + 1:]
    new_hydro_sum = sum(HYDROPATHY.get(new_protein[i], 0.0) for i in range(start, end))

    delta = (new_hydro_sum - old_hydro_sum) * 0.05
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
        from .optimization import optimize_sequence
        result = optimize_sequence(protein, organism=organism)
        return result.cai
    except Exception:
        return None


def _get_cai_weights(organism: str) -> dict[str, float]:
    """Return CAI adaptiveness table for organism, or empty dict."""
    try:
        from .organisms import CODON_ADAPTIVENESS_TABLES
        return CODON_ADAPTIVENESS_TABLES.get(organism, {})
    except ImportError:
        return {}


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

    # Penalty for long hydrophobic stretches (≥5 residues)
    stretch_penalty = 0.0
    run = 0
    for aa in protein:
        if aa in _AGGREGATION_PRONE:
            run += 1
            if run >= 5:
                stretch_penalty += 0.1
        else:
            run = 0

    score = -avg_hydro * 0.3 + charged_frac * 1.5 - stretch_penalty
    return round(score, 3)


def _base_immunogenicity(protein: str) -> float:
    """Compute a baseline immunogenicity score (0-1 scale).

    Uses average hydrophobicity of 9-mer windows as a proxy for MHC-II
    binding propensity. Higher = more immunogenic.
    """
    if len(protein) < 9:
        return 0.0
    scores = []
    for i in range(len(protein) - 8):
        window = protein[i:i + 9]
        avg_h = sum(HYDROPATHY.get(aa, 0.0) for aa in window) / 9.0
        scores.append(avg_h)
    if not scores:
        return 0.0
    # Normalize: typical range is -4.5 to +4.5; map to 0-1
    max_score = max(scores)
    normalized = (max_score + 4.5) / 9.0
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
    # Disulfide pairs contribute ~-2.0 kcal/mol each
    disulfide_pairs = cys_count // 2
    ddg = -avg_hydro * 0.5 - pro_frac * 2.0 + gly_frac * 1.5 - disulfide_pairs * 2.0
    return round(ddg, 3)


# ────────────────────────────────────────────────────────────
# Public API: Disulfide bond & proline substitution scanners
# ────────────────────────────────────────────────────────────

def find_disulfide_opportunities(protein: str) -> list[dict]:
    """Find positions where introducing cysteine pairs could form disulfide bonds.

    Criteria:
      - Both positions should be in predicted loop/coil regions (between
        secondary structure elements).
      - Positions must be at least 5 residues apart in sequence.
      - Neither position is already a cysteine.

    Returns:
        List of dicts with keys: position1, position2, distance_estimate,
        stabilizing_estimate.
    """
    ss = _predict_secondary_structure_simple(protein)
    n = len(protein)
    opportunities: list[dict] = []

    # Find loop positions (not already C)
    loop_positions = [i for i in range(n) if ss[i] == "L" and protein[i] != "C"]

    for idx_a in range(len(loop_positions)):
        for idx_b in range(idx_a + 1, len(loop_positions)):
            i = loop_positions[idx_a]
            j = loop_positions[idx_b]
            distance = abs(j - i)
            if distance < 5:
                continue
            # Estimate Cα-Cα distance from sequence separation (rough:
            # ~3.8 Å per residue for extended chain, loops are shorter)
            distance_estimate = distance * 3.5  # Å
            # Disulfide bonds are feasible for Cα-Cα distances of ~4-7 Å
            # In a loop, the effective distance is much shorter than sequence distance
            effective_distance = min(distance_estimate, 15.0)

            # Stabilizing estimate: depends on how constrained the loop currently is
            # and the sequence distance (shorter loops with disulfides are more stabilizing)
            stabilizing_estimate = max(-3.0, -2.0 + 0.02 * (distance - 10))
            # Short loops get more stabilization
            if distance <= 15:
                stabilizing_estimate -= 0.5

            opportunities.append({
                "position1": i,
                "position2": j,
                "distance_estimate": round(effective_distance, 1),
                "stabilizing_estimate": round(stabilizing_estimate, 3),
            })

    # Sort by most stabilizing first
    opportunities.sort(key=lambda x: x["stabilizing_estimate"])
    return opportunities


def find_proline_substitution_sites(protein: str) -> list[dict]:
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
    sites: list[dict] = []

    for i in range(n):
        if ss[i] != "L":
            continue
        wt = protein[i]
        if wt == "P" or wt == "G":
            continue
        blosum = BLOSUM62.get(wt, {}).get("P", -4)
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
) -> dict:
    """Score a single mutation across all dimensions.

    Args:
        protein: Original protein sequence (1-letter codes).
        position: 0-based residue position.
        mutant: Mutant amino acid (1-letter code).
        organism: Target organism for CAI lookup.
        weights: Optional weight dict with keys 'stability', 'solubility',
                 'immunogenicity'. Defaults to {0.4, 0.3, 0.3}.

    Returns:
        Dict with keys: stability_ddg, solubility_delta, immunogenicity_delta,
        blosum62, weighted_score.
    """
    if weights is None:
        weights = {"stability": 0.4, "solubility": 0.3, "immunogenicity": 0.3}

    wildtype = protein[position]
    ddg = _estimate_ddg(wildtype, mutant)
    sol_delta = _estimate_solubility_delta(wildtype, mutant)
    imm_delta = _estimate_immunogenicity_delta(protein, position, mutant)
    blosum = BLOSUM62.get(wildtype, {}).get(mutant, -4)

    # Normalize components for weighted score:
    #   stability: lower ddg is better → score = -ddg (positive = good)
    #   solubility: higher delta is better → score = delta
    #   immunogenicity: lower delta is better → score = -delta
    w_stab = weights.get("stability", 0.4)
    w_sol = weights.get("solubility", 0.3)
    w_imm = weights.get("immunogenicity", 0.3)

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
      5. Verify each mutation doesn't violate solubility/immunogenicity
         constraints.

    Args:
        protein: Input protein sequence (1-letter codes).
        organism: Target organism.
        target_stability: Target ΔG in kcal/mol (more negative = more stable).
        constraints: Design constraints.

    Returns:
        DesignResult with all mutations and metrics.
    """
    if constraints is None:
        constraints = DesignConstraints()

    current = list(protein)
    total_ddg = 0.0
    base_stab = _base_stability(protein)
    current_stab = base_stab
    mutations: list[dict] = []
    iterations = 0

    for iteration in range(constraints.max_mutations):
        iterations += 1
        if current_stab <= target_stability:
            break

        best_mutation = None
        best_ddg = 0.0  # only accept stabilizing (negative ddg)

        current_protein_str = "".join(current)
        current_sol = _base_solubility(current_protein_str)
        current_imm = _base_immunogenicity(current_protein_str)

        # --- Strategy 1: Scan all single-point mutations for best stabilizing one ---
        for pos in range(len(current)):
            wt = current[pos]
            if _is_preserved(pos, wt, constraints):
                continue
            for mutant_aa in _BLOSUM_INDEX:
                if mutant_aa == wt:
                    continue
                blosum = BLOSUM62.get(wt, {}).get(mutant_aa, -4)
                if blosum < constraints.blosum62_min:
                    continue
                ddg = _estimate_ddg(wt, mutant_aa)
                if ddg >= 0:
                    continue  # not stabilizing
                if abs(ddg) > constraints.max_ddg_per_mutation:
                    continue  # too large a change per step
                # Soft constraint checks: don't make things significantly worse
                sol_delta = _estimate_solubility_delta(wt, mutant_aa)
                imm_delta = _estimate_immunogenicity_delta(
                    current_protein_str, pos, mutant_aa
                )
                new_sol = current_sol + sol_delta
                new_imm = current_imm + imm_delta
                # Block only if the mutation worsens a constraint below a
                # tolerance margin (0.3 units) — allows gradual improvement
                if new_sol < constraints.min_solubility_score - 0.3:
                    continue
                if new_imm > constraints.max_immunogenicity + 0.3:
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
        proline_sites = find_proline_substitution_sites("".join(current))
        for site in proline_sites:
            pos = site["position"]
            wt = current[pos]
            if _is_preserved(pos, wt, constraints):
                continue
            ddg = site["ddg_estimate"]
            if ddg >= 0 or ddg >= best_ddg:
                continue
            blosum = BLOSUM62.get(wt, {}).get("P", -4)
            if blosum < constraints.blosum62_min:
                continue
            sol_delta = _estimate_solubility_delta(wt, "P")
            imm_delta = _estimate_immunogenicity_delta(
                current_protein_str, pos, "P"
            )
            new_sol = current_sol + sol_delta
            new_imm = current_imm + imm_delta
            if new_sol < constraints.min_solubility_score - 0.3:
                continue
            if new_imm > constraints.max_immunogenicity + 0.3:
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
        disulfide_ops = find_disulfide_opportunities("".join(current))
        for opp in disulfide_ops[:20]:  # limit search
            pos1 = opp["position1"]
            pos2 = opp["position2"]
            wt1 = current[pos1]
            wt2 = current[pos2]
            if _is_preserved(pos1, wt1, constraints):
                continue
            if _is_preserved(pos2, wt2, constraints):
                continue
            ddg_pair = opp["stabilizing_estimate"]
            if ddg_pair >= best_ddg:
                continue
            # Check constraints for both mutations combined
            sol1 = _estimate_solubility_delta(wt1, "C")
            sol2 = _estimate_solubility_delta(wt2, "C")
            imm1 = _estimate_immunogenicity_delta(current_protein_str, pos1, "C")
            imm2 = _estimate_immunogenicity_delta(current_protein_str, pos2, "C")
            new_sol = current_sol + sol1 + sol2
            new_imm = current_imm + imm1 + imm2
            if new_sol < constraints.min_solubility_score - 0.3:
                continue
            if new_imm > constraints.max_immunogenicity + 0.3:
                continue
            best_ddg = ddg_pair
            best_mutation = {
                "position": pos1,
                "wildtype": wt1,
                "mutant": "C",
                "ddg": ddg_pair / 2,
                "solubility_delta": sol1,
                "immunogenicity_delta": imm1,
                "blosum62": BLOSUM62.get(wt1, {}).get("C", -4),
                "strategy": "disulfide_pair",
                "pair_position": pos2,
                "pair_wildtype": wt2,
            }

        if best_mutation is None:
            logger.info("No more stabilizing mutations found at iteration %d", iteration)
            break

        # Apply the best mutation
        pos = best_mutation["position"]
        current[pos] = best_mutation["mutant"]
        total_ddg += best_mutation["ddg"]
        current_stab = base_stab + total_ddg
        mutations.append(best_mutation)

        # If disulfide pair, apply second Cysteine too
        if best_mutation.get("strategy") == "disulfide_pair":
            pair_pos = best_mutation["pair_position"]
            pair_wt = best_mutation["pair_wildtype"]
            current[pair_pos] = "C"
            total_ddg += best_mutation["ddg"]  # symmetric contribution
            current_stab = base_stab + total_ddg
            mutations.append({
                "position": pair_pos,
                "wildtype": pair_wt,
                "mutant": "C",
                "ddg": best_mutation["ddg"],
                "solubility_delta": _estimate_solubility_delta(pair_wt, "C"),
                "immunogenicity_delta": _estimate_immunogenicity_delta(
                    "".join(current), pair_pos, "C"
                ),
                "blosum62": BLOSUM62.get(pair_wt, {}).get("C", -4),
                "strategy": "disulfide_pair_partner",
            })

    designed = "".join(current)
    final_sol = _base_solubility(designed)
    final_imm = _base_immunogenicity(designed)
    cai = _compute_cai_for_protein(designed, organism)

    satisfied, violated = _check_constraints(
        designed, constraints, current_stab, final_sol, final_imm,
    )

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
        success=current_stab <= target_stability,
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
    if constraints is None:
        constraints = DesignConstraints()

    current = list(protein)
    base_stab = _base_stability(protein)
    current_stab = base_stab
    total_ddg = 0.0
    mutations: list[dict] = []
    iterations = 0

    for iteration in range(constraints.max_mutations):
        iterations += 1
        current_protein = "".join(current)
        current_sol = _base_solubility(current_protein)
        current_imm = _base_immunogenicity(current_protein)
        if current_sol >= min_solubility:
            break

        best_mutation = None
        best_sol_delta = 0.0

        for pos in range(len(current)):
            wt = current[pos]
            if _is_preserved(pos, wt, constraints):
                continue
            # Focus on hydrophobic / aggregation-prone residues
            if wt not in _HYDROPHOBIC_AAS and wt not in _AGGREGATION_PRONE:
                continue
            for mutant_aa in _SURFACE_FAVORED | _POLAR_AAS:
                if mutant_aa == wt:
                    continue
                blosum = BLOSUM62.get(wt, {}).get(mutant_aa, -4)
                if blosum < constraints.blosum62_min:
                    continue
                ddg = _estimate_ddg(wt, mutant_aa)
                if abs(ddg) > constraints.max_ddg_per_mutation:
                    continue
                # Soft stability check: only block if making stability
                # significantly worse than current (not vs. threshold)
                if ddg > constraints.max_ddg_per_mutation:
                    continue  # strongly destabilizing
                sol_delta = _estimate_solubility_delta(wt, mutant_aa)
                if sol_delta <= best_sol_delta:
                    continue
                imm_delta = _estimate_immunogenicity_delta(current_protein, pos, mutant_aa)
                new_imm = current_imm + imm_delta
                # Soft immunogenicity check: allow some tolerance
                if new_imm > constraints.max_immunogenicity + 0.3:
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

        if best_mutation is None:
            logger.info("No more solubility-improving mutations at iteration %d", iteration)
            break

        pos = best_mutation["position"]
        current[pos] = best_mutation["mutant"]
        total_ddg += best_mutation["ddg"]
        current_stab = base_stab + total_ddg
        mutations.append(best_mutation)

    designed = "".join(current)
    final_sol = _base_solubility(designed)
    final_imm = _base_immunogenicity(designed)
    cai = _compute_cai_for_protein(designed, organism)

    satisfied, violated = _check_constraints(
        designed, constraints, current_stab, final_sol, final_imm,
    )

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
        success=final_sol >= min_solubility,
    )


def design_low_immunogenicity(
    protein: str,
    organism: str = "Homo_sapiens",
    max_immunogenicity: float = 0.3,
    constraints: DesignConstraints | None = None,
) -> DesignResult:
    """Reduce immunogenicity via amino acid substitution.

    Wrapper around deimmunization logic with constraint checking.
    Additionally verifies stability and solubility aren't compromised.

    Strategy:
      - Identify high-immunogenicity 9-mer windows (hydrophobic-rich).
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
    if constraints is None:
        constraints = DesignConstraints()

    current = list(protein)
    base_stab = _base_stability(protein)
    current_stab = base_stab
    total_ddg = 0.0
    mutations: list[dict] = []
    iterations = 0

    for iteration in range(constraints.max_mutations):
        iterations += 1
        current_protein = "".join(current)
        current_imm = _base_immunogenicity(current_protein)
        if current_imm <= max_immunogenicity:
            break

        # Find the most immunogenic 9-mer window
        best_window_start = 0
        best_window_score = -999.0
        for i in range(len(current_protein) - 8):
            window = current_protein[i:i + 9]
            score = sum(HYDROPATHY.get(aa, 0.0) for aa in window) / 9.0
            if score > best_window_score:
                best_window_score = score
                best_window_start = i

        # Within that window, find the best mutation to reduce hydrophobicity
        best_mutation = None
        best_imm_delta = 0.0  # negative = reducing immunogenicity

        for pos in range(best_window_start, best_window_start + 9):
            if pos >= len(current):
                break
            wt = current[pos]
            if _is_preserved(pos, wt, constraints):
                continue
            for mutant_aa in _POLAR_AAS | _CHARGED_AAS | set("AST"):
                if mutant_aa == wt:
                    continue
                blosum = BLOSUM62.get(wt, {}).get(mutant_aa, -4)
                if blosum < constraints.blosum62_min:
                    continue
                ddg = _estimate_ddg(wt, mutant_aa)
                if abs(ddg) > constraints.max_ddg_per_mutation:
                    continue
                # Soft stability check: don't make things much worse
                if ddg > constraints.max_ddg_per_mutation:
                    continue
                imm_delta = _estimate_immunogenicity_delta(current_protein, pos, mutant_aa)
                if imm_delta >= best_imm_delta:
                    continue  # want negative (reducing immunogenicity)
                sol_delta = _estimate_solubility_delta(wt, mutant_aa)
                new_sol = _base_solubility(current_protein) + sol_delta
                # Soft solubility check: allow tolerance
                if new_sol < constraints.min_solubility_score - 0.3:
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

        if best_mutation is None:
            logger.info("No more deimmunizing mutations at iteration %d", iteration)
            break

        pos = best_mutation["position"]
        current[pos] = best_mutation["mutant"]
        total_ddg += best_mutation["ddg"]
        current_stab = base_stab + total_ddg
        mutations.append(best_mutation)

    designed = "".join(current)
    final_sol = _base_solubility(designed)
    final_imm = _base_immunogenicity(designed)
    cai = _compute_cai_for_protein(designed, organism)

    satisfied, violated = _check_constraints(
        designed, constraints, current_stab, final_sol, final_imm,
    )

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
        success=final_imm <= max_immunogenicity,
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
    if constraints is None:
        constraints = DesignConstraints()
    if weights is None:
        weights = {"stability": 0.4, "solubility": 0.3, "immunogenicity": 0.3}

    current = list(protein)
    base_stab = _base_stability(protein)
    current_stab = base_stab
    total_ddg = 0.0
    mutations: list[dict] = []
    iterations = 0

    for iteration in range(constraints.max_mutations):
        iterations += 1
        current_protein = "".join(current)
        current_sol = _base_solubility(current_protein)
        current_imm = _base_immunogenicity(current_protein)

        satisfied, violated = _check_constraints(
            current_protein, constraints, current_stab, current_sol, current_imm,
        )
        if not violated:
            break  # all constraints satisfied

        best_mutation = None
        best_weighted = -999.0

        for pos in range(len(current)):
            wt = current[pos]
            if _is_preserved(pos, wt, constraints):
                continue
            for mutant_aa in _BLOSUM_INDEX:
                if mutant_aa == wt:
                    continue
                blosum = BLOSUM62.get(wt, {}).get(mutant_aa, -4)
                if blosum < constraints.blosum62_min:
                    continue
                ddg = _estimate_ddg(wt, mutant_aa)
                if abs(ddg) > constraints.max_ddg_per_mutation:
                    continue
                # Soft stability check: only block if mutation makes stability
                # significantly worse than current state (not vs. threshold)
                if ddg > constraints.max_ddg_per_mutation:
                    continue  # strongly destabilizing

                sol_delta = _estimate_solubility_delta(wt, mutant_aa)
                imm_delta = _estimate_immunogenicity_delta(current_protein, pos, mutant_aa)

                new_sol = current_sol + sol_delta
                new_imm = current_imm + imm_delta

                # Soft constraints: prefer improvements but don't hard-block
                # unless the resulting value is clearly worse than threshold
                if new_sol < constraints.min_solubility_score - 0.5:
                    continue
                if new_imm > constraints.max_immunogenicity + 0.3:
                    continue

                # Compute weighted score
                w_stab = weights.get("stability", 0.4)
                w_sol = weights.get("solubility", 0.3)
                w_imm = weights.get("immunogenicity", 0.3)

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

        if best_mutation is None or best_weighted <= 0:
            logger.info(
                "No improving multi-objective mutation at iteration %d "
                "(best_weighted=%.4f)", iteration, best_weighted,
            )
            break

        pos = best_mutation["position"]
        current[pos] = best_mutation["mutant"]
        total_ddg += best_mutation["ddg"]
        current_stab = base_stab + total_ddg
        mutations.append(best_mutation)

    designed = "".join(current)
    final_sol = _base_solubility(designed)
    final_imm = _base_immunogenicity(designed)
    cai = _compute_cai_for_protein(designed, organism)

    satisfied, violated = _check_constraints(
        designed, constraints, current_stab, final_sol, final_imm,
    )

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
        success=len(violated) == 0,
    )
