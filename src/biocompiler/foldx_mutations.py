"""
BioCompiler FoldX Mutation Scanning Module v7.2.0
===================================================
Systematic mutation scanning and stability landscape analysis.

Provides empirical ΔΔG estimation for all point mutations in a protein
using BLOSUM62 substitution scores, hydrophobicity changes, and volume
changes.  Identifies stabilizing/destabilizing mutations, conserved
positions, compensatory mutations, and structural/functional hotspots.

Formula (empirical):
    ddg ≈ -0.1 * BLOSUM62(wt, mut) + 0.5 * |Δhydro| + 0.3 * |Δvolume|
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# BLOSUM62 substitution matrix  (20×20, symmetric)
# ────────────────────────────────────────────────────────────

BLOSUM62: dict[str, dict[str, int]] = {
    "A": {"A":  4, "R": -1, "N": -2, "D": -2, "C":  0, "Q": -1, "E": -1, "G":  0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "T":  0, "W": -3, "Y": -2, "V":  0},
    "R": {"A": -1, "R":  5, "N":  0, "D": -2, "C": -3, "Q":  1, "E":  0, "G": -2, "H":  0, "I": -3, "L": -2, "K":  2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"A": -2, "R":  0, "N":  6, "D":  1, "C": -3, "Q":  0, "E":  0, "G":  0, "H":  1, "I": -3, "L": -3, "K":  0, "M": -2, "F": -3, "P": -2, "S":  1, "T":  0, "W": -4, "Y": -2, "V": -3},
    "D": {"A": -2, "R": -2, "N":  1, "D":  6, "C": -3, "Q":  0, "E":  2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S":  0, "T": -1, "W": -4, "Y": -3, "V": -3},
    "C": {"A":  0, "R": -3, "N": -3, "D": -3, "C":  9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"A": -1, "R":  1, "N":  0, "D":  0, "C": -3, "Q":  5, "E":  2, "G": -2, "H":  0, "I": -3, "L": -2, "K":  1, "M":  0, "F": -3, "P": -1, "S":  0, "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"A": -1, "R":  0, "N":  0, "D":  2, "C": -4, "Q":  2, "E":  5, "G": -2, "H":  0, "I": -3, "L": -3, "K":  1, "M": -2, "F": -3, "P": -1, "S":  0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"A":  0, "R": -2, "N":  0, "D": -1, "C": -3, "Q": -2, "E": -2, "G":  6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S":  0, "T": -2, "W": -2, "Y": -3, "V": -3},
    "H": {"A": -2, "R":  0, "N":  1, "D": -1, "C": -3, "Q":  0, "E":  0, "G": -2, "H":  8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y":  2, "V": -3},
    "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I":  4, "L":  2, "K": -3, "M":  1, "F":  0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V":  3},
    "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I":  2, "L":  4, "K": -2, "M":  2, "F":  0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V":  1},
    "K": {"A": -1, "R":  2, "N":  0, "D": -1, "C": -3, "Q":  1, "E":  1, "G": -2, "H": -1, "I": -3, "L": -2, "K":  5, "M": -1, "F": -3, "P": -1, "S":  0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q":  0, "E": -2, "G": -3, "H": -2, "I":  1, "L":  2, "K": -1, "M":  5, "F":  0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V":  1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I":  0, "L":  0, "K": -3, "M":  0, "F":  6, "P": -4, "S": -2, "T": -2, "W":  1, "Y":  3, "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P":  7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A":  1, "R": -1, "N":  1, "D":  0, "C": -1, "Q":  0, "E":  0, "G":  0, "H": -1, "I": -2, "L": -2, "K":  0, "M": -1, "F": -2, "P": -1, "S":  4, "T":  1, "W": -3, "Y": -2, "V": -2},
    "T": {"A":  0, "R": -1, "N":  0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "T":  5, "W": -2, "Y": -2, "V":  0},
    "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F":  1, "P": -4, "S": -3, "T": -2, "W": 11, "Y":  2, "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H":  2, "I": -1, "L": -1, "K": -2, "M": -1, "F":  3, "P": -3, "S": -2, "T": -2, "W":  2, "Y":  7, "V": -1},
    "V": {"A":  0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I":  3, "L":  1, "K": -2, "M":  1, "F": -1, "P": -2, "S": -2, "T":  0, "W": -3, "Y": -1, "V":  4},
}

# ────────────────────────────────────────────────────────────
# Van der Waals volumes (Å³) — Creighton, 1993
# ────────────────────────────────────────────────────────────

AA_VOLUME: dict[str, float] = {
    "A":  88.6,  "R": 173.4,  "N": 114.1,  "D": 111.1,  "C": 108.5,
    "Q": 143.8,  "E": 138.4,  "G":  60.1,  "H": 153.2,  "I": 166.7,
    "L": 166.7,  "K": 168.6,  "M": 162.9,  "F": 189.9,  "P": 112.7,
    "S":  89.0,  "T": 116.1,  "W": 227.8,  "Y": 193.6,  "V": 140.0,
}

# ────────────────────────────────────────────────────────────
# Kyte-Doolittle hydropathy scale
# ────────────────────────────────────────────────────────────

HYDROPATHY: dict[str, float] = {
    "A":  1.8,  "R": -4.5,  "N": -3.5,  "D": -3.5,  "C":  2.5,
    "Q": -3.5,  "E": -3.5,  "G": -0.4,  "H": -3.2,  "I":  4.5,
    "L":  3.8,  "K": -3.9,  "M":  1.9,  "F":  2.8,  "P": -1.6,
    "S": -0.8,  "T": -0.7,  "W": -0.9,  "Y": -1.3,  "V":  4.2,
}

# ────────────────────────────────────────────────────────────
# Standard 20 amino acids
# ────────────────────────────────────────────────────────────

_AMINO_ACIDS: list[str] = list("ACDEFGHIKLMNPQRSTVWY")

# ΔΔG category thresholds (kcal/mol)
_STABILIZING_THRESHOLD = -0.5
_DESTABILIZING_THRESHOLD = 0.5

# Volume normalization factor: raw volumes are in Å³ (60–228), producing
# |Δvolume| up to ~168.  Dividing by 100 keeps the volume term in the
# same order of magnitude as the BLOSUM and hydropathy terms, yielding
# physically reasonable ΔΔG estimates (roughly -1 to +6 kcal/mol).
_VOLUME_SCALE = 100.0


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class StabilityLandscape:
    """Complete stability landscape for a protein across all point mutations."""

    protein: str
    wildtype_stability: float                # ΔG of wildtype (kcal/mol)
    mutations: list[dict]                     # [{position, wildtype, mutant, ddg, stabilizing}]
    stabilizing_count: int
    destabilizing_count: int
    neutral_count: int
    most_stabilizing: dict | None             # best mutation
    most_destabilizing: dict | None           # worst mutation
    positions_scanned: list[int]
    method: str                               # "empirical" or "foldx"


@dataclass
class ConservationScore:
    """Conservation analysis for a single protein position."""

    position: int
    wildtype: str
    conservation: float                       # 0 (variable) to 1 (fully conserved)
    substitution_tolerance: float             # average ΔΔG for all 19 substitutions
    critical: bool                            # conservation > 0.8 or avg_ddg > 3.0


# ────────────────────────────────────────────────────────────
# Core ΔΔG estimation
# ────────────────────────────────────────────────────────────

def _estimate_ddg(wt: str, mut: str) -> float:
    """Estimate ΔΔG for a single substitution using empirical formula.

    ddg ≈ -0.1 * BLOSUM62(wt, mut) + 0.5 * |Δhydro| + 0.3 * |Δvolume|/100

    The volume change is normalized by 100 (from Å³ to a unit that keeps
    the term in the same order of magnitude as BLOSUM and hydropathy),
    yielding physically reasonable ΔΔG estimates in kcal/mol.

    Positive ΔΔG → destabilizing; negative → stabilizing.
    """
    blosum = BLOSUM62.get(wt, {}).get(mut, -4)
    delta_hydro = abs(HYDROPATHY.get(wt, 0.0) - HYDROPATHY.get(mut, 0.0))
    delta_volume = abs(AA_VOLUME.get(wt, 0.0) - AA_VOLUME.get(mut, 0.0)) / _VOLUME_SCALE
    return -0.1 * blosum + 0.5 * delta_hydro + 0.3 * delta_volume


# ────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────

def scan_all_mutations(
    protein: str,
    method: str = "empirical",
) -> StabilityLandscape:
    """Scan all 19 substitutions at every position in *protein*.

    Args:
        protein: Amino-acid sequence (1-letter codes).
        method:  "empirical" for formula-based ΔΔG; "foldx" reserved
                 for future integration (falls back to empirical).

    Returns:
        A :class:`StabilityLandscape` with every mutation scored and
        categorized.
    """
    if not protein:
        return StabilityLandscape(
            protein="",
            wildtype_stability=0.0,
            mutations=[],
            stabilizing_count=0,
            destabilizing_count=0,
            neutral_count=0,
            most_stabilizing=None,
            most_destabilizing=None,
            positions_scanned=[],
            method=method,
        )

    protein = protein.upper()
    effective_method = method
    if method == "foldx":
        logger.warning(
            "FoldX backend not available; falling back to empirical estimation"
        )
        effective_method = "empirical"

    mutations: list[dict] = []
    positions_scanned: list[int] = []
    stabilizing = 0
    destabilizing = 0
    neutral = 0
    best: dict | None = None
    worst: dict | None = None

    for pos, wt in enumerate(protein):
        if wt not in _AMINO_ACIDS:
            continue
        positions_scanned.append(pos)
        for mut in _AMINO_ACIDS:
            if mut == wt:
                continue
            ddg = _estimate_ddg(wt, mut)
            is_stabilizing = ddg < _STABILIZING_THRESHOLD
            entry = {
                "position": pos,
                "wildtype": wt,
                "mutant": mut,
                "ddg": round(ddg, 4),
                "stabilizing": is_stabilizing,
            }
            mutations.append(entry)

            if is_stabilizing:
                stabilizing += 1
            elif ddg > _DESTABILIZING_THRESHOLD:
                destabilizing += 1
            else:
                neutral += 1

            if best is None or ddg < best["ddg"]:
                best = entry
            if worst is None or ddg > worst["ddg"]:
                worst = entry

    # Wildtype stability is set to 0.0 as the reference state
    return StabilityLandscape(
        protein=protein,
        wildtype_stability=0.0,
        mutations=mutations,
        stabilizing_count=stabilizing,
        destabilizing_count=destabilizing,
        neutral_count=neutral,
        most_stabilizing=best,
        most_destabilizing=worst,
        positions_scanned=positions_scanned,
        method=effective_method,
    )


def scan_position(
    protein: str,
    position: int,
    method: str = "empirical",
) -> list[dict]:
    """Scan all 19 substitutions at a single *position*.

    Args:
        protein:  Amino-acid sequence (1-letter codes).
        position: 0-based residue index.
        method:   "empirical" or "foldx" (falls back to empirical).

    Returns:
        List of mutation dicts sorted by ΔΔG ascending (most
        stabilizing first).
    """
    protein = protein.upper()
    if position < 0 or position >= len(protein):
        logger.warning("Position %d out of range for protein of length %d", position, len(protein))
        return []

    wt = protein[position]
    if wt not in _AMINO_ACIDS:
        logger.warning("Wildtype residue '%s' at position %d is not a standard AA", wt, position)
        return []

    results: list[dict] = []
    for mut in _AMINO_ACIDS:
        if mut == wt:
            continue
        ddg = _estimate_ddg(wt, mut)
        is_stabilizing = ddg < _STABILIZING_THRESHOLD
        results.append({
            "position": position,
            "wildtype": wt,
            "mutant": mut,
            "ddg": round(ddg, 4),
            "stabilizing": is_stabilizing,
        })

    results.sort(key=lambda m: m["ddg"])
    return results


def compute_conservation(
    protein: str,
    method: str = "empirical",
) -> list[ConservationScore]:
    """Compute conservation score for every position.

    Conservation = 1 - (num_tolerated_substitutions / 19)
    A tolerated substitution has ΔΔG < 1.0 kcal/mol.
    A position is *critical* if conservation > 0.8 **or** average ΔΔG > 3.0.

    Args:
        protein: Amino-acid sequence (1-letter codes).
        method:  "empirical" or "foldx" (falls back to empirical).

    Returns:
        List of :class:`ConservationScore` objects, one per residue.
    """
    protein = protein.upper()
    scores: list[ConservationScore] = []

    for pos, wt in enumerate(protein):
        if wt not in _AMINO_ACIDS:
            continue

        ddgs: list[float] = []
        tolerated = 0
        for mut in _AMINO_ACIDS:
            if mut == wt:
                continue
            ddg = _estimate_ddg(wt, mut)
            ddgs.append(ddg)
            if ddg < 1.0:
                tolerated += 1

        conservation = 1.0 - (tolerated / 19.0)
        avg_ddg = sum(ddgs) / len(ddgs) if ddgs else 0.0
        critical = conservation > 0.8 or avg_ddg > 3.0

        scores.append(ConservationScore(
            position=pos,
            wildtype=wt,
            conservation=round(conservation, 4),
            substitution_tolerance=round(avg_ddg, 4),
            critical=critical,
        ))

    return scores


def find_compensatory_mutations(
    protein: str,
    destabilizing_mutations: list[dict],
) -> list[dict]:
    """Find second-site compensatory mutations for destabilizing variants.

    A compensatory mutation is one that, when combined with the original
    destabilizing mutation, reduces the total ΔΔG.  The heuristic looks
    for stabilizing mutations at positions within ±5 residues of the
    original mutation.

    Args:
        protein:               Amino-acid sequence.
        destabilizing_mutations: List of dicts each with keys
                                ``position``, ``wildtype``, ``mutant``,
                                ``ddg``.

    Returns:
        List of dicts with keys ``position``, ``original_mutation``,
        ``compensatory_mutation``, ``combined_ddg``.
    """
    if not protein or not destabilizing_mutations:
        return []

    protein = protein.upper()

    # Pre-compute per-position mutation lists for efficiency
    position_mutations: dict[int, list[dict]] = {}
    for pos, wt in enumerate(protein):
        if wt not in _AMINO_ACIDS:
            continue
        muts: list[dict] = []
        for mut in _AMINO_ACIDS:
            if mut == wt:
                continue
            ddg = _estimate_ddg(wt, mut)
            muts.append({
                "position": pos,
                "wildtype": wt,
                "mutant": mut,
                "ddg": round(ddg, 4),
                "stabilizing": ddg < _STABILIZING_THRESHOLD,
            })
        position_mutations[pos] = muts

    results: list[dict] = []

    for dm in destabilizing_mutations:
        dm_pos = dm.get("position", -1)
        dm_ddg = dm.get("ddg", 0.0)

        if dm_pos < 0 or dm_pos >= len(protein):
            continue

        # Search nearby positions (within 5 residues)
        best_comp: dict | None = None
        best_combined = dm_ddg  # start with no compensation

        for offset in range(-5, 6):
            if offset == 0:
                continue
            nearby = dm_pos + offset
            if nearby < 0 or nearby >= len(protein):
                continue
            if nearby not in position_mutations:
                continue

            for cm in position_mutations[nearby]:
                combined = dm_ddg + cm["ddg"]
                # A compensatory mutation must reduce total ΔΔG and be at
                # least mildly stabilising (ddg < 0) on its own.
                if combined < best_combined and cm["ddg"] < 0:
                    best_combined = combined
                    best_comp = cm

        if best_comp is not None:
            results.append({
                "position": best_comp["position"],
                "original_mutation": {
                    "position": dm_pos,
                    "wildtype": dm.get("wildtype", ""),
                    "mutant": dm.get("mutant", ""),
                    "ddg": dm_ddg,
                },
                "compensatory_mutation": {
                    "position": best_comp["position"],
                    "wildtype": best_comp["wildtype"],
                    "mutant": best_comp["mutant"],
                    "ddg": best_comp["ddg"],
                },
                "combined_ddg": round(best_combined, 4),
            })

    return results


def rank_positions_by_mutability(
    protein: str,
) -> list[tuple[int, float]]:
    """Rank positions from most to least mutable.

    Mutability score = average ΔΔG across all 19 substitutions.
    Lower scores → more mutable (easier to change without destabilizing).

    Args:
        protein: Amino-acid sequence.

    Returns:
        List of ``(position, avg_ddg)`` sorted by *avg_ddg* ascending
        (most mutable first).
    """
    protein = protein.upper()
    rankings: list[tuple[int, float]] = []

    for pos, wt in enumerate(protein):
        if wt not in _AMINO_ACIDS:
            continue

        ddgs: list[float] = []
        for mut in _AMINO_ACIDS:
            if mut == wt:
                continue
            ddgs.append(_estimate_ddg(wt, mut))

        avg_ddg = sum(ddgs) / len(ddgs) if ddgs else 0.0
        rankings.append((pos, round(avg_ddg, 4)))

    rankings.sort(key=lambda x: x[1])
    return rankings


def identify_hotspot_regions(
    protein: str,
    window: int = 5,
    threshold: float = 2.0,
) -> list[tuple[int, int]]:
    """Find contiguous regions where average ΔΔG exceeds *threshold*.

    Hotspots are structural/functional regions that are hard to mutate
    without destabilizing the protein.  A sliding window of *window*
    residues is used; if the average ΔΔG of all 19 substitutions across
    all positions in the window exceeds *threshold*, the window is
    flagged.  Overlapping flagged windows are merged into contiguous
    (start, end) intervals.

    Args:
        protein:   Amino-acid sequence.
        window:    Sliding-window size (number of residues).
        threshold: Average ΔΔG above which a window is a hotspot.

    Returns:
        List of ``(start, end)`` position tuples (0-based, inclusive).
    """
    protein = protein.upper()
    if not protein or window < 1:
        return []

    # Compute per-position average ΔΔG
    pos_avg: dict[int, float] = {}
    for pos, wt in enumerate(protein):
        if wt not in _AMINO_ACIDS:
            continue
        ddgs = [_estimate_ddg(wt, mut) for mut in _AMINO_ACIDS if mut != wt]
        pos_avg[pos] = sum(ddgs) / len(ddgs) if ddgs else 0.0

    # Sliding window scan
    hot_positions: set[int] = set()
    for start in range(len(protein) - window + 1):
        window_positions = [p for p in range(start, start + window) if p in pos_avg]
        if len(window_positions) < window:
            continue
        window_avg = sum(pos_avg[p] for p in window_positions) / len(window_positions)
        if window_avg > threshold:
            hot_positions.update(window_positions)

    # Merge contiguous positions into (start, end) intervals
    if not hot_positions:
        return []

    sorted_pos = sorted(hot_positions)
    regions: list[tuple[int, int]] = []
    region_start = sorted_pos[0]
    region_end = sorted_pos[0]

    for p in sorted_pos[1:]:
        if p == region_end + 1:
            region_end = p
        else:
            regions.append((region_start, region_end))
            region_start = p
            region_end = p
    regions.append((region_start, region_end))

    return regions
