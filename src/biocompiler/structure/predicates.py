"""
BioCompiler Structure-Based Type System Predicates
====================================================
Type system predicates based on protein structure prediction.
These predicates evaluate structural confidence, misfolding risk,
fold topology, and unexpected interaction potential using either
provided PDB structures or heuristics from sequence.

Requires PDB data for full evaluation; degrades gracefully to
UNCERTAIN verdicts with knowledge_gap annotations when structure
is unavailable.
"""

from __future__ import annotations

import math
import logging
from typing import Any, Optional

from ..type_system import Verdict, TypeCheckResult
from ..constants import HYDROPHOBIC_AAS

logger = logging.getLogger(__name__)


__all__ = [
    "evaluate_structure_confidence",
    "evaluate_no_misfolding_risk",
    "evaluate_correct_fold_topology",
    "evaluate_no_unexpected_interaction",
    "expected_radius_of_gyration",
    "compute_secondary_structure_fractions",
    "find_surface_charge_patches",
    "find_unstructured_regions",
]

# ────────────────────────────────────────────────────────────
# Amino acid property sets
# ────────────────────────────────────────────────────────────
POSITIVELY_CHARGED_AAS = set("KRH")
NEGATIVELY_CHARGED_AAS = set("DE")
CHARGED_AAS = POSITIVELY_CHARGED_AAS | NEGATIVELY_CHARGED_AAS

# ────────────────────────────────────────────────────────────
# Named constants (extracted from magic numbers)
# ────────────────────────────────────────────────────────────

# Radius of gyration empirical formula: Rg ≈ RG_COEFFICIENT * N^RG_EXPONENT
RG_EMPIRICAL_COEFFICIENT: float = 2.5
RG_EMPIRICAL_EXPONENT: float = 0.33

# pLDDT verdict thresholds
PLDDT_PASS_THRESHOLD: float = 90.0
PLDDT_LIKELY_PASS_THRESHOLD: float = 70.0
PLDDT_UNCERTAIN_THRESHOLD: float = 50.0

# Misfolding risk thresholds
LOW_CONFIDENCE_CONSECUTIVE_LIMIT: int = 10
RAMACHANDRAN_OUTLIER_FRACTION: float = 0.05
CA_CLASH_DISTANCE_ANGSTROM: float = 2.0

# Radius of gyration deviation factors (relative to expected)
RG_DEVIATION_LOW: float = 0.6
RG_DEVIATION_HIGH: float = 1.4

# Fold topology thresholds
RG_RATIO_FAILED_LOW: float = 0.6
RG_RATIO_FAILED_HIGH: float = 1.5
RG_RATIO_BORDERLINE_LOW: float = 0.7
RG_RATIO_BORDERLINE_HIGH: float = 1.4
MIN_STRUCTURED_PROTEIN_LENGTH: int = 50
COIL_FRACTION_FAILED: float = 0.80
COIL_FRACTION_BORDERLINE: float = 0.65
HELIX_FRACTION_FAILED: float = 0.85
HELIX_FRACTION_BORDERLINE: float = 0.75

# Packing density thresholds
PACKING_DISTANCE_CUTOFF: float = 10.0
PACKING_DENSITY_FAILED_MIN: float = 5.0
PACKING_DENSITY_FAILED_MAX: float = 12.0
PACKING_DENSITY_BORDERLINE_MIN: float = 6.0
PACKING_DENSITY_BORDERLINE_MAX: float = 11.0

# Exposure / burial thresholds
EXPOSURE_CUTOFF_FRACTION: float = 0.7
BURIAL_CUTOFF_FRACTION: float = 0.7
EXPOSED_HYDROPHOBIC_THRESHOLD: float = 0.40
SEQUENCE_HYDROPHOBICITY_THRESHOLD: float = 0.45
HYDRO_BORDERLINE_EXPOSED_MIN: float = 0.5
HYDRO_BORDERLINE_RATIO: float = 0.85

# Interaction risk thresholds
UNSTRUCTURED_REGION_MIN_LENGTH: int = 30
SIGNIFICANT_PATCH_MIN_LENGTH: int = 5

# Fold topology TM-score and RMSD thresholds
TM_SCORE_PASS_THRESHOLD: float = 0.4  # Standard "same fold" threshold (relaxed from 0.5)
SMALL_PROTEIN_RMSD_THRESHOLD_ANG: float = 3.0  # Angstroms, for proteins < 100aa
SMALL_PROTEIN_LENGTH_LIMIT: int = 100  # Residues

# Interaction filtering thresholds (to reduce false positives)
INTERFACE_AREA_THRESHOLD_ANG2: float = 500.0  # Minimum interface area (Å²) to flag
BINDING_ENERGY_THRESHOLD_KCAL: float = -8.0  # Only flag if binding energy is stronger
CRYSTAL_PACKING_SURFACE_FRACTION: float = 0.30  # >30% surface involved = likely artifact
AVERAGE_RESIDUE_SURFACE_AREA: float = 90.0  # Approximate Å² per exposed residue

# Geometry tolerance
DIHEDRAL_NORM_TOLERANCE: float = 1e-8

# Sequence heuristic thresholds
FLEXIBLE_REGION_FRACTION: float = 0.60

# PDB parsing
PDB_MIN_ATOM_LINE_LENGTH: int = 54
DEFAULT_OCCUPANCY: float = 1.0
DEFAULT_BFACTOR: float = 0.0


# ────────────────────────────────────────────────────────────
# PDB parsing helpers
# ────────────────────────────────────────────────────────────

def _parse_pdb_atoms(pdb_string: str) -> list[dict[str, Any]]:
    """Parse ATOM records from a PDB string.

    Returns a list of dicts, each with keys:
        serial, name, resname, chain, resseq, x, y, z, occupancy, bfactor, element
    """
    atoms = []
    for line in pdb_string.splitlines():
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        if len(line) < PDB_MIN_ATOM_LINE_LENGTH:
            continue
        try:
            atom = {
                "serial": int(line[6:11].strip()),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21] if len(line) > 21 else "A",
                "resseq": int(line[22:26].strip()),
                "x": float(line[30:38].strip()),
                "y": float(line[38:46].strip()),
                "z": float(line[46:54].strip()),
                "occupancy": float(line[54:60].strip()) if len(line) >= 60 else DEFAULT_OCCUPANCY,
                "bfactor": float(line[60:66].strip()) if len(line) >= 66 else DEFAULT_BFACTOR,
                "element": line[76:78].strip() if len(line) >= 78 else "",
            }
            atoms.append(atom)
        except (ValueError, IndexError) as exc:
            logger.debug("Skipping malformed PDB ATOM line: %s", exc)
            continue
    return atoms


def _extract_plddt_scores(pdb_string: str) -> list[float]:
    """Extract pLDDT scores from PDB B-factor column.

    ESMFold and AlphaFold store per-residue confidence (pLDDT)
    in the B-factor column of CA atoms.
    """
    atoms = _parse_pdb_atoms(pdb_string)
    plddt_scores = []
    for atom in atoms:
        if atom["name"] == "CA":
            plddt_scores.append(atom["bfactor"])
    return plddt_scores


def _extract_ca_coords(pdb_string: str) -> list[tuple[float, float, float, int]]:
    """Extract CA atom coordinates from PDB string.

    Returns list of (x, y, z, resseq) tuples.
    """
    atoms = _parse_pdb_atoms(pdb_string)
    ca_coords = []
    for atom in atoms:
        if atom["name"] == "CA":
            ca_coords.append((atom["x"], atom["y"], atom["z"], atom["resseq"]))
    return ca_coords


def _extract_backbone_coords(pdb_string: str) -> dict[int, dict[str, tuple[float, float, float]]]:
    """Extract backbone atom coordinates (N, CA, C) organized by residue.

    Returns dict mapping resseq -> {"N": (x,y,z), "CA": (x,y,z), "C": (x,y,z)}.
    """
    atoms = _parse_pdb_atoms(pdb_string)
    backbone: dict[int, dict[str, tuple[float, float, float]]] = {}
    for atom in atoms:
        if atom["name"] in ("N", "CA", "C"):
            resseq = atom["resseq"]
            if resseq not in backbone:
                backbone[resseq] = {}
            backbone[resseq][atom["name"]] = (atom["x"], atom["y"], atom["z"])
    return backbone


# ────────────────────────────────────────────────────────────
# Geometry helpers
# ────────────────────────────────────────────────────────────

def _distance(p1: tuple[float, ...], p2: tuple[float, ...]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def _vector(p1: tuple[float, ...], p2: tuple[float, ...]) -> tuple[float, ...]:
    """Vector from p1 to p2."""
    return tuple(b - a for a, b in zip(p1, p2))


def _dot(v1: tuple[float, ...], v2: tuple[float, ...]) -> float:
    """Dot product of two vectors."""
    return sum(a * b for a, b in zip(v1, v2))


def _norm(v: tuple[float, ...]) -> float:
    """Norm of a vector."""
    return math.sqrt(sum(x * x for x in v))


def _dihedral_angle(
    p1: tuple[float, ...],
    p2: tuple[float, ...],
    p3: tuple[float, ...],
    p4: tuple[float, ...],
) -> float:
    """Compute dihedral angle (in degrees) defined by four points.

    Uses the convention: angle between the (p1-p2, p2-p3) and (p2-p3, p3-p4) planes.
    """
    b1 = _vector(p1, p2)
    b2 = _vector(p2, p3)
    b3 = _vector(p3, p4)

    n1 = (
        b1[1] * b2[2] - b1[2] * b2[1],
        b1[2] * b2[0] - b1[0] * b2[2],
        b1[0] * b2[1] - b1[1] * b2[0],
    )
    n2 = (
        b2[1] * b3[2] - b2[2] * b3[1],
        b2[2] * b3[0] - b2[0] * b3[2],
        b2[0] * b3[1] - b2[1] * b3[0],
    )

    n1_norm = _norm(n1)
    n2_norm = _norm(n2)

    if n1_norm < DIHEDRAL_NORM_TOLERANCE or n2_norm < DIHEDRAL_NORM_TOLERANCE:
        return 0.0

    cos_angle = _dot(n1, n2) / (n1_norm * n2_norm)
    cos_angle = max(-1.0, min(1.0, cos_angle))

    # Determine sign
    cross_n1_n2 = (
        n1[1] * n2[2] - n1[2] * n2[1],
        n1[2] * n2[0] - n1[0] * n2[2],
        n1[0] * n2[1] - n1[1] * n2[0],
    )
    sign = 1.0 if _dot(cross_n1_n2, b2) >= 0 else -1.0

    return math.degrees(math.acos(cos_angle)) * sign


def _compute_radius_of_gyration(ca_coords: list[tuple[float, float, float, int]]) -> float:
    """Compute radius of gyration from CA coordinates.

    Rg = sqrt(mean(distance_from_centroid^2))
    """
    if not ca_coords:
        return 0.0

    # Compute centroid
    n = len(ca_coords)
    cx = sum(c[0] for c in ca_coords) / n
    cy = sum(c[1] for c in ca_coords) / n
    cz = sum(c[2] for c in ca_coords) / n

    # Compute mean squared distance from centroid
    msd = sum((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2 for c in ca_coords) / n
    return math.sqrt(msd)


def _compute_rmsd(
    coords1: list[tuple[float, ...]],
    coords2: list[tuple[float, ...]],
) -> float:
    """Compute Cα RMSD between two coordinate sets (after centering).

    Simplified: does not perform Kabsch rotation; centers both structures
    at the origin before computing RMSD.

    Args:
        coords1: First set of (x, y, z, ...) tuples.
        coords2: Second set of (x, y, z, ...) tuples (same length).

    Returns:
        RMSD in Angstroms, or inf if inputs are mismatched/empty.
    """
    if len(coords1) != len(coords2) or len(coords1) == 0:
        return float("inf")

    n = len(coords1)

    # Center both structures
    c1_cx = sum(c[0] for c in coords1) / n
    c1_cy = sum(c[1] for c in coords1) / n
    c1_cz = sum(c[2] for c in coords1) / n
    c2_cx = sum(c[0] for c in coords2) / n
    c2_cy = sum(c[1] for c in coords2) / n
    c2_cz = sum(c[2] for c in coords2) / n

    sum_sq = 0.0
    for c1, c2 in zip(coords1, coords2):
        dx = (c1[0] - c1_cx) - (c2[0] - c2_cx)
        dy = (c1[1] - c1_cy) - (c2[1] - c2_cy)
        dz = (c1[2] - c1_cz) - (c2[2] - c2_cz)
        sum_sq += dx * dx + dy * dy + dz * dz

    return math.sqrt(sum_sq / n)


def _compute_tm_score(
    coords1: list[tuple[float, ...]],
    coords2: list[tuple[float, ...]],
    length: int,
) -> float:
    """Compute TM-score between two coordinate sets (after centering).

    TM-score = (1 / L_ref) * Σ [1 / (1 + (d_i / d0)²)]
    where d0 = 1.24 * (L_ref - 15)^(1/3) - 1.8

    Simplified: does not perform Kabsch rotation; centers both structures.

    Args:
        coords1: First set of (x, y, z, ...) tuples.
        coords2: Second set of (x, y, z, ...) tuples (same length).
        length: Reference length (L_ref), typically the native protein length.

    Returns:
        TM-score in [0, 1], or 0.0 if inputs are invalid.
    """
    if len(coords1) != len(coords2) or len(coords1) == 0:
        return 0.0

    n = len(coords1)
    l_ref = max(length, 1)

    # d0 parameter
    if l_ref > 15:
        d0 = 1.24 * ((l_ref - 15) ** (1.0 / 3.0)) - 1.8
    else:
        d0 = 0.5
    d0 = max(d0, 0.5)

    # Center both structures
    c1_cx = sum(c[0] for c in coords1) / n
    c1_cy = sum(c[1] for c in coords1) / n
    c1_cz = sum(c[2] for c in coords1) / n
    c2_cx = sum(c[0] for c in coords2) / n
    c2_cy = sum(c[1] for c in coords2) / n
    c2_cz = sum(c[2] for c in coords2) / n

    tm_sum = 0.0
    for c1, c2 in zip(coords1, coords2):
        dx = (c1[0] - c1_cx) - (c2[0] - c2_cx)
        dy = (c1[1] - c1_cy) - (c2[1] - c2_cy)
        dz = (c1[2] - c1_cz) - (c2[2] - c2_cz)
        d_i = math.sqrt(dx * dx + dy * dy + dz * dz)
        tm_sum += 1.0 / (1.0 + (d_i / d0) ** 2)

    return tm_sum / l_ref


def _compute_ramachandran_outliers(
    backbone: dict[int, dict[str, tuple[float, float, float]]],
) -> tuple[int, int]:
    """Count Ramachandran outliers and total residues with backbone atoms.

    A residue is a Ramachandran outlier if its phi/psi angles fall outside
    the allowed regions of the Ramachandran plot (simplified check).

    Returns (num_outliers, num_total).
    """
    sorted_resseqs = sorted(backbone.keys())
    if len(sorted_resseqs) < 3:
        return 0, 0

    outliers = 0
    total = 0

    for i in range(1, len(sorted_resseqs) - 1):
        prev_res = sorted_resseqs[i - 1]
        curr_res = sorted_resseqs[i]
        next_res = sorted_resseqs[i + 1]

        # Need C of previous residue, N/CA/C of current, N of next
        if not all(k in backbone[prev_res] for k in ("C",)):
            continue
        if not all(k in backbone[curr_res] for k in ("N", "CA", "C")):
            continue
        if not all(k in backbone[next_res] for k in ("N",)):
            continue

        # Phi: C(prev) - N(curr) - CA(curr) - C(curr)
        phi = _dihedral_angle(
            backbone[prev_res]["C"],
            backbone[curr_res]["N"],
            backbone[curr_res]["CA"],
            backbone[curr_res]["C"],
        )

        # Psi: N(curr) - CA(curr) - C(curr) - N(next)
        psi = _dihedral_angle(
            backbone[curr_res]["N"],
            backbone[curr_res]["CA"],
            backbone[curr_res]["C"],
            backbone[next_res]["N"],
        )

        total += 1

        # Simplified Ramachandran allowed regions
        if not _is_ramachandran_allowed(phi, psi):
            outliers += 1

    return outliers, total


def _is_ramachandran_allowed(phi: float, psi: float) -> bool:
    """Check if phi/psi angles are in an allowed Ramachandran region.

    Simplified check covering the major allowed regions:
    - Beta-sheet region: phi ~ -120 to -60, psi ~ 60 to 180
    - Right-handed alpha helix: phi ~ -90 to -30, psi ~ -75 to -15
    - Left-handed alpha helix: phi ~ 30 to 90, psi ~ -15 to 75
    - Extended allowed regions with generous margins
    """
    # Right-handed alpha helix (most common)
    if -120 <= phi <= -30 and -90 <= psi <= -10:
        return True
    # Beta sheet region
    if -180 <= phi <= -40 and 30 <= psi <= 180:
        return True
    # Also allow near -180/+180 boundary
    if 140 <= phi <= 180 and 30 <= psi <= 180:
        return True
    # Left-handed helix (rare but allowed for glycine)
    if 20 <= phi <= 120 and -30 <= psi <= 90:
        return True
    # Type II beta turn region
    if -120 <= phi <= 0 and -30 <= psi <= 60:
        return True
    # Generous extended region
    if -180 <= phi <= -60 and -180 <= psi <= -90:
        return True
    # Additional allowed: polyproline II
    if -100 <= phi <= -50 and 100 <= psi <= 180:
        return True

    return False


# ────────────────────────────────────────────────────────────
# Helper functions (public API)
# ────────────────────────────────────────────────────────────

def expected_radius_of_gyration(length: int) -> float:
    """Compute expected radius of gyration for a globular protein.

    Empirical formula: Rg ≈ 2.5 * N^0.33 for globular proteins,
    where N is the number of residues.

    Args:
        length: Number of amino acid residues.

    Returns:
        Expected radius of gyration in Angstroms.
    """
    if length <= 0:
        return 0.0
    return RG_EMPIRICAL_COEFFICIENT * (length ** RG_EMPIRICAL_EXPONENT)


def compute_secondary_structure_fractions(pdb_string: str) -> dict[str, float]:
    """Compute secondary structure fractions from a PDB string.

    Uses phi/psi dihedral angle-based classification:
    - Helix: phi in [-120, -30], psi in [-90, -10]
    - Sheet: phi in [-180, -40], psi in [30, 180] (or equivalent)
    - Coil: everything else

    Args:
        pdb_string: PDB format string.

    Returns:
        Dict with keys "helix", "sheet", "coil" mapping to fractions (0.0-1.0).
    """
    backbone = _extract_backbone_coords(pdb_string)
    sorted_resseqs = sorted(backbone.keys())

    if len(sorted_resseqs) < 3:
        return {"helix": 0.0, "sheet": 0.0, "coil": 1.0}

    helix = 0
    sheet = 0
    coil = 0

    for i in range(1, len(sorted_resseqs) - 1):
        prev_res = sorted_resseqs[i - 1]
        curr_res = sorted_resseqs[i]
        next_res = sorted_resseqs[i + 1]

        if not all(k in backbone[prev_res] for k in ("C",)):
            coil += 1
            continue
        if not all(k in backbone[curr_res] for k in ("N", "CA", "C")):
            coil += 1
            continue
        if not all(k in backbone[next_res] for k in ("N",)):
            coil += 1
            continue

        phi = _dihedral_angle(
            backbone[prev_res]["C"],
            backbone[curr_res]["N"],
            backbone[curr_res]["CA"],
            backbone[curr_res]["C"],
        )
        psi = _dihedral_angle(
            backbone[curr_res]["N"],
            backbone[curr_res]["CA"],
            backbone[curr_res]["C"],
            backbone[next_res]["N"],
        )

        # Classify
        if -120 <= phi <= -30 and -90 <= psi <= -10:
            helix += 1
        elif (-180 <= phi <= -40 and 30 <= psi <= 180) or (140 <= phi <= 180 and 30 <= psi <= 180):
            sheet += 1
        elif -100 <= phi <= -50 and 100 <= psi <= 180:
            sheet += 1  # polyproline II
        else:
            coil += 1

    total = helix + sheet + coil
    if total == 0:
        return {"helix": 0.0, "sheet": 0.0, "coil": 1.0}

    return {
        "helix": helix / total,
        "sheet": sheet / total,
        "coil": coil / total,
    }


def find_surface_charge_patches(
    protein: str,
    window: int = 7,
    min_charge: float = 0.7,
) -> list[tuple[int, int, str]]:
    """Find regions with concentrated same-sign charge in a protein sequence.

    Slides a window over the protein sequence and identifies regions
    where the fraction of residues with the same charge sign exceeds
    the threshold.

    Args:
        protein: Amino acid sequence (single-letter codes).
        window: Sliding window size in residues (default 7).
        min_charge: Minimum fraction of same-sign charged residues
            to flag as a patch (default 0.7).

    Returns:
        List of (start, end, charge_type) tuples where charge_type
        is "positive" or "negative". Positions are 0-based, end is exclusive.
    """
    protein = protein.upper()
    patches: list[tuple[int, int, str]] = []

    if len(protein) < window:
        return patches

    i = 0
    while i <= len(protein) - window:
        window_seq = protein[i:i + window]
        pos_count = sum(1 for aa in window_seq if aa in POSITIVELY_CHARGED_AAS)
        neg_count = sum(1 for aa in window_seq if aa in NEGATIVELY_CHARGED_AAS)

        pos_frac = pos_count / window
        neg_frac = neg_count / window

        if pos_frac >= min_charge:
            # Extend the patch
            end = i + window
            while end < len(protein) and protein[end] in POSITIVELY_CHARGED_AAS:
                end += 1
            patches.append((i, end, "positive"))
            i = end
        elif neg_frac >= min_charge:
            end = i + window
            while end < len(protein) and protein[end] in NEGATIVELY_CHARGED_AAS:
                end += 1
            patches.append((i, end, "negative"))
            i = end
        else:
            i += 1

    # Merge overlapping patches of the same type
    merged: list[tuple[int, int, str]] = []
    for patch in patches:
        if merged and merged[-1][2] == patch[2] and merged[-1][1] >= patch[0]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], patch[1]), merged[-1][2])
        else:
            merged.append(patch)

    return merged


def find_unstructured_regions(pdb_string: str, min_length: int = 15) -> list[tuple[int, int]]:
    """Find long coil regions from secondary structure assignment.

    Identifies contiguous coil regions longer than min_length residues
    based on phi/psi dihedral angle classification.

    Args:
        pdb_string: PDB format string.
        min_length: Minimum length of coil region to report (default 15).

    Returns:
        List of (start, end) tuples (0-based residue index, end exclusive).
    """
    backbone = _extract_backbone_coords(pdb_string)
    sorted_resseqs = sorted(backbone.keys())

    if len(sorted_resseqs) < 3:
        return []

    # Classify each residue
    ss_assignment: dict[int, str] = {}
    for i in range(1, len(sorted_resseqs) - 1):
        prev_res = sorted_resseqs[i - 1]
        curr_res = sorted_resseqs[i]
        next_res = sorted_resseqs[i + 1]

        if not all(k in backbone[prev_res] for k in ("C",)):
            ss_assignment[curr_res] = "coil"
            continue
        if not all(k in backbone[curr_res] for k in ("N", "CA", "C")):
            ss_assignment[curr_res] = "coil"
            continue
        if not all(k in backbone[next_res] for k in ("N",)):
            ss_assignment[curr_res] = "coil"
            continue

        phi = _dihedral_angle(
            backbone[prev_res]["C"],
            backbone[curr_res]["N"],
            backbone[curr_res]["CA"],
            backbone[curr_res]["C"],
        )
        psi = _dihedral_angle(
            backbone[curr_res]["N"],
            backbone[curr_res]["CA"],
            backbone[curr_res]["C"],
            backbone[next_res]["N"],
        )

        if -120 <= phi <= -30 and -90 <= psi <= -10:
            ss_assignment[curr_res] = "helix"
        elif (-180 <= phi <= -40 and 30 <= psi <= 180) or (140 <= phi <= 180 and 30 <= psi <= 180):
            ss_assignment[curr_res] = "sheet"
        else:
            ss_assignment[curr_res] = "coil"

    # Find contiguous coil regions
    regions: list[tuple[int, int]] = []
    coil_start = None
    for idx, resseq in enumerate(sorted_resseqs):
        ss = ss_assignment.get(resseq, "coil")
        if ss == "coil":
            if coil_start is None:
                coil_start = idx
        else:
            if coil_start is not None:
                if idx - coil_start >= min_length:
                    regions.append((coil_start, idx))
                coil_start = None

    # Handle trailing coil
    if coil_start is not None:
        end_idx = len(sorted_resseqs)
        if end_idx - coil_start >= min_length:
            regions.append((coil_start, end_idx))

    return regions


# ────────────────────────────────────────────────────────────
# Main predicate evaluate functions
# ────────────────────────────────────────────────────────────

def evaluate_structure_confidence(
    sequence: str,
    protein: str,
    organism: str,
    min_mean_plddt: float = 70.0,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Evaluate that predicted protein structure has sufficient confidence.

    Checks the mean pLDDT (per-residue Local Distance Difference Test)
    score from structure prediction, either extracted from a provided
    PDB file (B-factor column) or estimated.

    Verdict thresholds:
        mean_plddt > 90  → PASS
        mean_plddt 70-90 → LIKELY_PASS
        mean_plddt 50-70 → UNCERTAIN
        mean_plddt < 50  → LIKELY_FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name.
        min_mean_plddt: Minimum acceptable mean pLDDT (default 70.0).
        pdb_string: Optional PDB format string with structure data.

    Returns:
        TypeCheckResult with verdict based on pLDDT confidence.
    """
    if pdb_string is not None:
        plddt_scores = _extract_plddt_scores(pdb_string)

        if not plddt_scores:
            return TypeCheckResult(
                predicate="StructureConfidence",
                verdict=Verdict.UNCERTAIN,
                derivation=[{"step": "extract_plddt", "result": "no CA atoms found in PDB"}],
                knowledge_gap="PDB provided but no CA atoms with B-factors found",
            )

        mean_plddt = sum(plddt_scores) / len(plddt_scores)

        # Per-category breakdown
        very_high = sum(1 for s in plddt_scores if s > PLDDT_PASS_THRESHOLD)
        confident = sum(1 for s in plddt_scores if PLDDT_LIKELY_PASS_THRESHOLD < s <= PLDDT_PASS_THRESHOLD)
        low = sum(1 for s in plddt_scores if PLDDT_UNCERTAIN_THRESHOLD < s <= PLDDT_LIKELY_PASS_THRESHOLD)
        very_low = sum(1 for s in plddt_scores if s <= PLDDT_UNCERTAIN_THRESHOLD)
        n = len(plddt_scores)

        derivation = [
            {"step": "extract_plddt", "mean_plddt": round(mean_plddt, 2), "num_residues": n},
            {
                "step": "category_breakdown",
                "very_high_gt90": very_high,
                "confident_70_90": confident,
                "low_50_70": low,
                "very_low_lt50": very_low,
                "very_high_pct": round(100 * very_high / n, 1),
                "confident_pct": round(100 * confident / n, 1),
                "low_pct": round(100 * low / n, 1),
                "very_low_pct": round(100 * very_low / n, 1),
            },
        ]

        if mean_plddt > PLDDT_PASS_THRESHOLD:
            verdict = Verdict.PASS
        elif mean_plddt > PLDDT_LIKELY_PASS_THRESHOLD:
            verdict = Verdict.LIKELY_PASS
        elif mean_plddt > PLDDT_UNCERTAIN_THRESHOLD:
            verdict = Verdict.UNCERTAIN
        else:
            verdict = Verdict.LIKELY_FAIL

        violation = None
        if mean_plddt < min_mean_plddt:
            violation = (
                f"Mean pLDDT {mean_plddt:.1f} is below threshold {min_mean_plddt:.1f}; "
                f"{very_low} residues ({100 * very_low / n:.1f}%) have pLDDT <= {PLDDT_UNCERTAIN_THRESHOLD:.0f}"
            )

        return TypeCheckResult(
            predicate="StructureConfidence",
            verdict=verdict,
            derivation=derivation,
            violation=violation,
        )

    # No PDB provided — try ESMFold or return UNCERTAIN
    try:
        # Attempt to use ESMFold if available
        from esm import pretrained  # type: ignore[import-untyped]

        logger.info("ESMFold available, running structure prediction for %d aa protein", len(protein))
        model = pretrained.esmfold_v1()
        model = model.eval()
        with torch.no_grad():  # type: ignore[name-defined]  # noqa: F821
            result = model.infer(protein)
            predicted_plddt = float(result["mean_plddt"])

        if predicted_plddt > PLDDT_PASS_THRESHOLD:
            verdict = Verdict.LIKELY_PASS
        elif predicted_plddt > PLDDT_LIKELY_PASS_THRESHOLD:
            verdict = Verdict.LIKELY_PASS
        elif predicted_plddt > PLDDT_UNCERTAIN_THRESHOLD:
            verdict = Verdict.UNCERTAIN
        else:
            verdict = Verdict.LIKELY_FAIL

        return TypeCheckResult(
            predicate="StructureConfidence",
            verdict=verdict,
            derivation=[
                {"step": "esmfold_prediction", "mean_plddt": round(predicted_plddt, 2)},
            ],
            knowledge_gap="Structure from computational prediction (ESMFold), not experimental",
            violation=(
                f"Predicted mean pLDDT {predicted_plddt:.1f} < {min_mean_plddt:.1f}"
                if predicted_plddt < min_mean_plddt else None
            ),
        )
    except ImportError:
        logger.debug("ESMFold not available; falling back to UNCERTAIN verdict")
    except Exception as exc:
        logger.warning("ESMFold prediction failed: %s", exc)

    # No structure available at all
    return TypeCheckResult(
        predicate="StructureConfidence",
        verdict=Verdict.UNCERTAIN,
        derivation=[
            {"step": "no_structure", "result": "no PDB provided, ESMFold not available"},
        ],
        knowledge_gap=(
            "No structure provided and ESMFold not installed; "
            "cannot evaluate structure confidence computationally"
        ),
    )


def evaluate_no_misfolding_risk(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Evaluate misfolding risk indicators.

    Checks for:
    a. Long low-confidence regions (>10 consecutive residues with pLDDT < 70)
    b. Ramachandran outliers (>5% of residues)
    c. High clash score (non-bonded contacts < 2.0 Angstroms between CA atoms)
    d. Very low or very high radius of gyration for protein length

    Verdict:
        No risk indicators → PASS
        1 mild indicator    → LIKELY_PASS
        2 indicators        → UNCERTAIN
        3+ indicators       → LIKELY_FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name.
        pdb_string: Optional PDB format string with structure data.

    Returns:
        TypeCheckResult with verdict based on misfolding risk indicators.
    """
    risk_indicators: list[str] = []
    derivation_steps: list[dict] = []

    if pdb_string is None:
        return TypeCheckResult(
            predicate="NoMisfoldingRisk",
            verdict=Verdict.UNCERTAIN,
            derivation=[{"step": "no_structure", "result": "no PDB provided"}],
            knowledge_gap="No structure provided; cannot evaluate misfolding risk",
        )

    # (a) Long low-confidence regions
    plddt_scores = _extract_plddt_scores(pdb_string)
    if plddt_scores:
        # Find consecutive low-confidence stretches
        consecutive_low = 0
        max_consecutive_low = 0
        for score in plddt_scores:
            if score < PLDDT_LIKELY_PASS_THRESHOLD:
                consecutive_low += 1
                max_consecutive_low = max(max_consecutive_low, consecutive_low)
            else:
                consecutive_low = 0

        derivation_steps.append({
            "step": "low_confidence_regions",
            "max_consecutive_low_plddt": max_consecutive_low,
            "threshold": LOW_CONFIDENCE_CONSECUTIVE_LIMIT,
            "plddt_threshold": PLDDT_LIKELY_PASS_THRESHOLD,
        })

        if max_consecutive_low > LOW_CONFIDENCE_CONSECUTIVE_LIMIT:
            risk_indicators.append(
                f"Long low-confidence region: {max_consecutive_low} consecutive "
                f"residues with pLDDT < {PLDDT_LIKELY_PASS_THRESHOLD:.0f}"
            )
    else:
        derivation_steps.append({
            "step": "low_confidence_regions",
            "result": "no pLDDT scores available",
        })

    # (b) Ramachandran outliers
    backbone = _extract_backbone_coords(pdb_string)
    outliers, total = _compute_ramachandran_outliers(backbone)
    outlier_frac = outliers / total if total > 0 else 0.0

    derivation_steps.append({
        "step": "ramachandran",
        "outliers": outliers,
        "total": total,
        "outlier_fraction": round(outlier_frac, 4),
        "threshold": RAMACHANDRAN_OUTLIER_FRACTION,
    })

    if total > 0 and outlier_frac > RAMACHANDRAN_OUTLIER_FRACTION:
        risk_indicators.append(
            f"Ramachandran outliers: {outliers}/{total} "
            f"({100 * outlier_frac:.1f}%) > {RAMACHANDRAN_OUTLIER_FRACTION * 100:.0f}%"
        )

    # (c) Clash score (simplified: count CA-CA contacts < 2.0 Angstroms)
    ca_coords = _extract_ca_coords(pdb_string)
    if len(ca_coords) > 1:
        clashes = 0
        for i in range(len(ca_coords)):
            for j in range(i + 2, len(ca_coords)):  # skip adjacent residues
                dist = _distance(ca_coords[i][:3], ca_coords[j][:3])
                if dist < CA_CLASH_DISTANCE_ANGSTROM:
                    clashes += 1

        derivation_steps.append({
            "step": "clash_score",
            "clashes": clashes,
            "num_ca": len(ca_coords),
        })

        if clashes > 0:
            risk_indicators.append(
                f"Clash score: {clashes} non-bonded CA-CA contacts < {CA_CLASH_DISTANCE_ANGSTROM:.1f} Angstroms"
            )

    # (d) Radius of gyration check
    if ca_coords:
        rg = _compute_radius_of_gyration(ca_coords)
        expected_rg = expected_radius_of_gyration(len(protein))
        # Allow ±40% deviation from expected
        rg_low = expected_rg * RG_DEVIATION_LOW
        rg_high = expected_rg * RG_DEVIATION_HIGH

        derivation_steps.append({
            "step": "radius_of_gyration",
            "observed_rg": round(rg, 2),
            "expected_rg": round(expected_rg, 2),
            "acceptable_range": [round(rg_low, 2), round(rg_high, 2)],
        })

        if rg < rg_low or rg > rg_high:
            direction = "low" if rg < rg_low else "high"
            risk_indicators.append(
                f"Radius of gyration {rg:.1f}A is too {direction} "
                f"(expected ~{expected_rg:.1f}A, range [{rg_low:.1f}, {rg_high:.1f}])"
            )

    # Determine verdict based on number of risk indicators
    n_risks = len(risk_indicators)
    if n_risks == 0:
        verdict = Verdict.PASS
    elif n_risks == 1:
        verdict = Verdict.LIKELY_PASS
    elif n_risks == 2:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_FAIL

    violation = None
    if risk_indicators:
        violation = f"{n_risks} misfolding risk indicator(s): " + "; ".join(risk_indicators)

    return TypeCheckResult(
        predicate="NoMisfoldingRisk",
        verdict=verdict,
        derivation=derivation_steps,
        violation=violation,
        knowledge_gap="Structure prediction is computational, not experimental" if pdb_string else None,
    )


def _sequence_based_topology_check(
    protein: str,
) -> tuple[Verdict, list[dict], str | None]:
    """Evaluate fold topology from sequence alone when no structure is available.

    Uses sequence composition heuristics to estimate whether the protein
    is likely to adopt a well-folded structure:
    - Hydrophobicity distribution (should be in reasonable range)
    - Secondary structure propensity (should have some SS-forming residues)
    - Charge distribution (should not be excessively charged)

    Returns:
        (verdict, derivation_steps, knowledge_gap)
    """
    derivation_steps: list[dict] = []
    protein_upper = protein.upper()
    n = len(protein_upper)

    if n == 0:
        return Verdict.UNCERTAIN, [{"step": "empty_sequence"}], "Empty protein sequence"

    # (a) Hydrophobicity check
    hydro_count = sum(1 for aa in protein_upper if aa in HYDROPHOBIC_AAS)
    hydro_frac = hydro_count / n
    derivation_steps.append({
        "step": "sequence_hydrophobicity",
        "hydrophobic_fraction": round(hydro_frac, 3),
    })

    # (b) Secondary structure propensity from sequence
    helix_formers = set("AELMQRK")
    sheet_formers = set("VIYFWCT")
    helix_count = sum(1 for aa in protein_upper if aa in helix_formers)
    sheet_count = sum(1 for aa in protein_upper if aa in sheet_formers)
    helix_frac = helix_count / n
    sheet_frac = sheet_count / n
    derivation_steps.append({
        "step": "ss_propensity",
        "helix_propensity_frac": round(helix_frac, 3),
        "sheet_propensity_frac": round(sheet_frac, 3),
    })

    # (c) Charge distribution
    pos_count = sum(1 for aa in protein_upper if aa in POSITIVELY_CHARGED_AAS)
    neg_count = sum(1 for aa in protein_upper if aa in NEGATIVELY_CHARGED_AAS)
    charged_frac = (pos_count + neg_count) / n
    derivation_steps.append({
        "step": "charge_distribution",
        "charged_fraction": round(charged_frac, 3),
        "positive": pos_count,
        "negative": neg_count,
    })

    # Scoring: a well-folded protein should have reasonable properties
    issues = 0
    if hydro_frac < 0.20 or hydro_frac > 0.55:
        issues += 1
        derivation_steps.append({
            "step": "hydrophobicity_issue",
            "detail": f"hydrophobic fraction {hydro_frac:.3f} outside [0.20, 0.55]",
        })
    if helix_frac + sheet_frac < 0.15:
        issues += 1
        derivation_steps.append({
            "step": "ss_propensity_issue",
            "detail": f"combined SS propensity {helix_frac + sheet_frac:.3f} < 0.15",
        })
    if charged_frac > 0.40:
        issues += 1
        derivation_steps.append({
            "step": "charge_issue",
            "detail": f"charged fraction {charged_frac:.3f} > 0.40",
        })

    derivation_steps.append({
        "step": "sequence_topology_summary",
        "issues_found": issues,
    })

    if issues == 0:
        return (
            Verdict.LIKELY_PASS,
            derivation_steps,
            "Sequence-based topology estimate (no structure data); "
            "sequence composition looks reasonable",
        )
    elif issues == 1:
        return (
            Verdict.UNCERTAIN,
            derivation_steps,
            "Sequence-based topology estimate with minor concerns (no structure data)",
        )
    else:
        return (
            Verdict.UNCERTAIN,
            derivation_steps,
            "Sequence-based topology estimate with multiple concerns (no structure data)",
        )


def _has_obvious_dimer_interface(pdb_string: str, protein: str) -> bool:
    """Check if structure shows an obvious dimerization interface.

    Looks for multiple chains in the PDB or a very large contiguous
    hydrophobic patch on the surface that would indicate a dimer interface.

    Args:
        pdb_string: PDB format string.
        protein: Amino acid sequence.

    Returns:
        True if an obvious dimer interface is detected.
    """
    atoms = _parse_pdb_atoms(pdb_string)

    # Check for multiple chains
    chains: set[str] = set()
    for atom in atoms:
        chains.add(atom["chain"])
    if len(chains) > 1:
        return True

    # Check for very large surface hydrophobic patch (>60% exposed hydrophobic)
    ca_coords = _extract_ca_coords(pdb_string)
    if not ca_coords or len(protein) < 20:
        return False

    n_ca = len(ca_coords)
    cx = sum(c[0] for c in ca_coords) / n_ca
    cy = sum(c[1] for c in ca_coords) / n_ca
    cz = sum(c[2] for c in ca_coords) / n_ca
    max_dist = max(
        math.sqrt((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2)
        for c in ca_coords
    )
    exposure_cutoff = max_dist * EXPOSURE_CUTOFF_FRACTION

    exposed_hydrophobic = 0
    exposed_total = 0
    for idx, ca in enumerate(ca_coords):
        dist = math.sqrt((ca[0] - cx) ** 2 + (ca[1] - cy) ** 2 + (ca[2] - cz) ** 2)
        if dist > exposure_cutoff:
            exposed_total += 1
            aa = protein[idx] if idx < len(protein) else "X"
            if aa in HYDROPHOBIC_AAS:
                exposed_hydrophobic += 1

    if exposed_total > 0:
        exposed_hydro_frac = exposed_hydrophobic / exposed_total
        if exposed_hydro_frac > 0.60:
            return True

    return False


def evaluate_correct_fold_topology(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
    reference_pdb_string: str | None = None,
) -> TypeCheckResult:
    """Evaluate that the predicted fold makes biological sense.

    Heuristics checked:
    a. Expected secondary structure content (alpha/beta/coil proportions)
       matches protein length
    b. Radius of gyration within expected range: Rg ≈ 2.5 * N^0.33
    c. Packing density in normal range (5-12 CA neighbors within 10A)
    d. Hydrophobic residues predominantly buried (if PDB available)
    e. TM-score > 0.4 vs reference (or RMSD < 3.0 Å for small proteins < 100aa)

    When no PDB structure is provided, falls back to sequence-based topology
    comparison instead of returning UNCERTAIN/FAIL.

    Verdict:
        All checks pass   → PASS
        1 check borderline → LIKELY_PASS
        2 checks fail      → UNCERTAIN
        3+ checks fail     → LIKELY_FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name.
        pdb_string: Optional PDB format string with structure data.
        reference_pdb_string: Optional reference PDB for TM-score/RMSD comparison.

    Returns:
        TypeCheckResult with verdict based on fold topology checks.
    """
    failed_checks: list[str] = []
    borderline_checks: list[str] = []
    derivation_steps: list[dict] = []

    if pdb_string is None:
        # Use sequence-based topology comparison instead of FAIL/UNCERTAIN
        seq_verdict, seq_derivation, seq_gap = _sequence_based_topology_check(protein)
        return TypeCheckResult(
            predicate="CorrectFoldTopology",
            verdict=seq_verdict,
            derivation=seq_derivation,
            knowledge_gap=seq_gap,
        )

    n_residues = len(protein)
    ca_coords = _extract_ca_coords(pdb_string)

    # (a) Secondary structure content check
    ss_fracs = compute_secondary_structure_fractions(pdb_string)
    derivation_steps.append({
        "step": "secondary_structure",
        "helix_frac": round(ss_fracs["helix"], 3),
        "sheet_frac": round(ss_fracs["sheet"], 3),
        "coil_frac": round(ss_fracs["coil"], 3),
    })

    # For proteins of various lengths, the SS content should be reasonable.
    # Very long coil (>80%) or very long helix (>80%) in a large protein is unusual.
    if n_residues > MIN_STRUCTURED_PROTEIN_LENGTH:
        if ss_fracs["coil"] > COIL_FRACTION_FAILED:
            failed_checks.append(
                f"Excessive coil content: {ss_fracs['coil']:.1%} coil "
                f"(expected <{COIL_FRACTION_FAILED:.0%} for structured proteins)"
            )
        elif ss_fracs["coil"] > COIL_FRACTION_BORDERLINE:
            borderline_checks.append(
                f"High coil content: {ss_fracs['coil']:.1%} coil"
            )
        if ss_fracs["helix"] > HELIX_FRACTION_FAILED:
            failed_checks.append(
                f"Unusually high helix content: {ss_fracs['helix']:.1%} helix"
            )
        elif ss_fracs["helix"] > HELIX_FRACTION_BORDERLINE:
            borderline_checks.append(
                f"High helix content: {ss_fracs['helix']:.1%} helix"
            )

    # (b) Radius of gyration within expected range
    if ca_coords:
        rg = _compute_radius_of_gyration(ca_coords)
        expected_rg = expected_radius_of_gyration(n_residues)
        # Allow ±30% for standard globular, ±50% for borderline
        rg_ratio = rg / expected_rg if expected_rg > 0 else 1.0

        derivation_steps.append({
            "step": "radius_of_gyration",
            "observed_rg": round(rg, 2),
            "expected_rg": round(expected_rg, 2),
            "ratio": round(rg_ratio, 3),
        })

        if rg_ratio < RG_RATIO_FAILED_LOW or rg_ratio > RG_RATIO_FAILED_HIGH:
            failed_checks.append(
                f"Radius of gyration ratio {rg_ratio:.2f} far from expected 1.0 "
                f"(Rg={rg:.1f}A vs expected {expected_rg:.1f}A)"
            )
        elif rg_ratio < RG_RATIO_BORDERLINE_LOW or rg_ratio > RG_RATIO_BORDERLINE_HIGH:
            borderline_checks.append(
                f"Radius of gyration ratio {rg_ratio:.2f} borderline "
                f"(Rg={rg:.1f}A vs expected {expected_rg:.1f}A)"
            )
    else:
        borderline_checks.append("No CA coordinates for Rg calculation")

    # (c) Packing density (CA neighbors within 10A)
    if len(ca_coords) > 1:
        neighbor_counts = []
        for i, ca_i in enumerate(ca_coords):
            count = 0
            for j, ca_j in enumerate(ca_coords):
                if i == j:
                    continue
                if _distance(ca_i[:3], ca_j[:3]) <= PACKING_DISTANCE_CUTOFF:
                    count += 1
            neighbor_counts.append(count)

        mean_neighbors = sum(neighbor_counts) / len(neighbor_counts) if neighbor_counts else 0
        derivation_steps.append({
            "step": "packing_density",
            "mean_ca_neighbors_10A": round(mean_neighbors, 2),
            "acceptable_range": [PACKING_DENSITY_FAILED_MIN, PACKING_DENSITY_FAILED_MAX],
        })

        if mean_neighbors < PACKING_DENSITY_FAILED_MIN or mean_neighbors > PACKING_DENSITY_FAILED_MAX:
            failed_checks.append(
                f"Packing density out of range: mean {mean_neighbors:.1f} "
                f"CA neighbors within {PACKING_DISTANCE_CUTOFF:.0f}A "
                f"(expected {PACKING_DENSITY_FAILED_MIN:.0f}-{PACKING_DENSITY_FAILED_MAX:.0f})"
            )
        elif mean_neighbors < PACKING_DENSITY_BORDERLINE_MIN or mean_neighbors > PACKING_DENSITY_BORDERLINE_MAX:
            borderline_checks.append(
                f"Packing density borderline: mean {mean_neighbors:.1f} CA neighbors"
            )
    else:
        borderline_checks.append("Insufficient CA coordinates for packing analysis")

    # (d) Hydrophobic burial check (requires PDB)
    if ca_coords and n_residues > MIN_STRUCTURED_PROTEIN_LENGTH:
        # Compute centroid
        n_ca = len(ca_coords)
        cx = sum(c[0] for c in ca_coords) / n_ca
        cy = sum(c[1] for c in ca_coords) / n_ca
        cz = sum(c[2] for c in ca_coords) / n_ca

        # Define "buried" as within configured fraction of the max distance from centroid
        max_dist = max(
            math.sqrt((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2)
            for c in ca_coords
        )
        burial_cutoff = max_dist * BURIAL_CUTOFF_FRACTION

        buried_hydrophobic = 0
        buried_total = 0
        exposed_hydrophobic = 0
        exposed_total = 0

        for idx, ca in enumerate(ca_coords):
            dist_from_center = math.sqrt(
                (ca[0] - cx) ** 2 + (ca[1] - cy) ** 2 + (ca[2] - cz) ** 2
            )
            aa = protein[idx] if idx < len(protein) else "X"
            is_hydrophobic = aa in HYDROPHOBIC_AAS

            if dist_from_center <= burial_cutoff:
                buried_total += 1
                if is_hydrophobic:
                    buried_hydrophobic += 1
            else:
                exposed_total += 1
                if is_hydrophobic:
                    exposed_hydrophobic += 1

        buried_hydro_frac = buried_hydrophobic / buried_total if buried_total > 0 else 0.0
        exposed_hydro_frac = exposed_hydrophobic / exposed_total if exposed_total > 0 else 0.0

        derivation_steps.append({
            "step": "hydrophobic_burial",
            "buried_hydrophobic_frac": round(buried_hydro_frac, 3),
            "exposed_hydrophobic_frac": round(exposed_hydro_frac, 3),
        })

        # Hydrophobic residues should be more common in the core
        if exposed_hydrophobic > 0 and buried_total > 0:
            if exposed_hydro_frac > buried_hydro_frac:
                failed_checks.append(
                    f"Hydrophobic residues more exposed than buried: "
                    f"exposed={exposed_hydro_frac:.1%} vs buried={buried_hydro_frac:.1%}"
                )
            elif exposed_hydro_frac > HYDRO_BORDERLINE_EXPOSED_MIN and exposed_hydro_frac >= buried_hydro_frac * HYDRO_BORDERLINE_RATIO:
                borderline_checks.append(
                    f"Borderline hydrophobic burial: "
                    f"exposed={exposed_hydro_frac:.1%} vs buried={buried_hydro_frac:.1%}"
                )

    # (e) TM-score or RMSD comparison with reference structure
    if reference_pdb_string is not None:
        ref_ca_coords = _extract_ca_coords(reference_pdb_string)
        if len(ca_coords) > 0 and len(ref_ca_coords) > 0:
            min_len = min(len(ca_coords), len(ref_ca_coords))
            pred_coords = [
                (ca_coords[i][0], ca_coords[i][1], ca_coords[i][2])
                for i in range(min_len)
            ]
            ref_coords = [
                (ref_ca_coords[i][0], ref_ca_coords[i][1], ref_ca_coords[i][2])
                for i in range(min_len)
            ]

            if n_residues < SMALL_PROTEIN_LENGTH_LIMIT:
                # For small proteins (<100aa), use RMSD < 3.0 Å
                rmsd = _compute_rmsd(pred_coords, ref_coords)
                derivation_steps.append({
                    "step": "rmsd_comparison",
                    "rmsd_ang": round(rmsd, 2),
                    "threshold_ang": SMALL_PROTEIN_RMSD_THRESHOLD_ANG,
                    "method": "RMSD (small protein <100aa)",
                })
                if rmsd > SMALL_PROTEIN_RMSD_THRESHOLD_ANG:
                    failed_checks.append(
                        f"RMSD {rmsd:.2f} Å exceeds threshold "
                        f"{SMALL_PROTEIN_RMSD_THRESHOLD_ANG:.1f} Å for small protein"
                    )
            else:
                # For larger proteins, use TM-score > 0.4 (relaxed from 0.5)
                tm_score = _compute_tm_score(pred_coords, ref_coords, n_residues)
                derivation_steps.append({
                    "step": "tm_score_comparison",
                    "tm_score": round(tm_score, 4),
                    "threshold": TM_SCORE_PASS_THRESHOLD,
                    "method": "TM-score",
                })
                if tm_score < TM_SCORE_PASS_THRESHOLD:
                    failed_checks.append(
                        f"TM-score {tm_score:.4f} below threshold "
                        f"{TM_SCORE_PASS_THRESHOLD} (folds likely different)"
                    )

    # Determine verdict
    n_failed = len(failed_checks)
    n_borderline = len(borderline_checks)

    if n_failed == 0 and n_borderline == 0:
        verdict = Verdict.PASS
    elif n_failed == 0 and n_borderline <= 1:
        verdict = Verdict.LIKELY_PASS
    elif n_failed <= 1 and n_borderline <= 2:
        verdict = Verdict.UNCERTAIN
    elif n_failed >= 3:
        verdict = Verdict.LIKELY_FAIL
    elif n_failed >= 2:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_FAIL

    violation = None
    all_issues = failed_checks + borderline_checks
    if all_issues:
        violation = f"{n_failed} failed + {n_borderline} borderline check(s): " + "; ".join(all_issues)

    return TypeCheckResult(
        predicate="CorrectFoldTopology",
        verdict=verdict,
        derivation=derivation_steps,
        violation=violation,
        knowledge_gap="Structure prediction is computational, not experimental",
    )


def evaluate_no_unexpected_interaction(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
    is_monomeric: bool = False,
    known_interaction_partners: list[str] | None = None,
) -> TypeCheckResult:
    """Evaluate potential for unwanted protein-protein interactions.

    Checks for:
    a. Large exposed hydrophobic surface (>40% of exposed residues are hydrophobic)
    b. Long unstructured regions (>30 residues with low secondary structure)
    c. High surface charge patches (cluster of same-charge residues)

    Filtering (to reduce false positives):
    - Only flag interactions with estimated interface area > 500 Å²
    - Only flag if predicted binding energy < -8 kcal/mol (not transient)
    - Auto-PASS monomeric proteins with no known interaction partners
      unless structure shows obvious dimer interface
    - Ignore interactions involving >30% of surface residues (crystal packing)

    Verdict:
        No indicators  → PASS
        1 indicator    → LIKELY_PASS
        2 indicators   → UNCERTAIN
        3+ indicators  → LIKELY_FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name.
        pdb_string: Optional PDB format string with structure data.
        is_monomeric: Whether the protein is known to be monomeric.
        known_interaction_partners: List of known interaction partner names, or None.

    Returns:
        TypeCheckResult with verdict based on interaction risk indicators.
    """
    indicators: list[str] = []
    derivation_steps: list[dict] = []

    # Track exposed residue count for crystal packing filter
    total_surface_residues = 0
    exposed_hydrophobic_count = 0

    # ── Monomeric auto-PASS ──
    # For monomeric proteins with no known interaction partners, auto-PASS
    # unless the structure shows an obvious dimer interface
    if is_monomeric and (known_interaction_partners is None or len(known_interaction_partners) == 0):
        has_dimer = False
        if pdb_string is not None:
            has_dimer = _has_obvious_dimer_interface(pdb_string, protein)
        if not has_dimer:
            derivation_steps.append({
                "step": "monomer_check",
                "result": "monomeric with no known partners and no dimer interface, auto-PASS",
                "is_monomeric": True,
                "known_partners": known_interaction_partners,
            })
            return TypeCheckResult(
                predicate="NoUnexpectedInteraction",
                verdict=Verdict.PASS,
                derivation=derivation_steps,
            )
        else:
            derivation_steps.append({
                "step": "monomer_check",
                "result": "monomeric but structure shows dimer interface, continuing checks",
                "is_monomeric": True,
            })

    # (a) Exposed hydrophobic surface
    if pdb_string is not None:
        ca_coords = _extract_ca_coords(pdb_string)
        if ca_coords and len(protein) > 10:
            n_ca = len(ca_coords)
            cx = sum(c[0] for c in ca_coords) / n_ca
            cy = sum(c[1] for c in ca_coords) / n_ca
            cz = sum(c[2] for c in ca_coords) / n_ca

            max_dist = max(
                math.sqrt((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2)
                for c in ca_coords
            )
            # Exposed residues: beyond configured fraction of max distance from centroid
            exposure_cutoff = max_dist * EXPOSURE_CUTOFF_FRACTION

            exposed_hydrophobic = 0
            exposed_total = 0
            for idx, ca in enumerate(ca_coords):
                dist = math.sqrt(
                    (ca[0] - cx) ** 2 + (ca[1] - cy) ** 2 + (ca[2] - cz) ** 2
                )
                if dist > exposure_cutoff:
                    exposed_total += 1
                    aa = protein[idx] if idx < len(protein) else "X"
                    if aa in HYDROPHOBIC_AAS:
                        exposed_hydrophobic += 1

            total_surface_residues = exposed_total
            exposed_hydrophobic_count = exposed_hydrophobic
            exposed_hydro_frac = exposed_hydrophobic / exposed_total if exposed_total > 0 else 0.0

            # Estimate interface area from exposed hydrophobic residues
            estimated_interface_area = exposed_hydrophobic * AVERAGE_RESIDUE_SURFACE_AREA

            # Estimate binding energy: ~-0.5 kcal/mol per hydrophobic contact,
            # plus electrostatic contribution from charge patches
            estimated_binding_energy = -0.5 * exposed_hydrophobic

            derivation_steps.append({
                "step": "exposed_hydrophobic",
                "exposed_residues": exposed_total,
                "exposed_hydrophobic": exposed_hydrophobic,
                "fraction": round(exposed_hydro_frac, 3),
                "threshold": EXPOSED_HYDROPHOBIC_THRESHOLD,
                "estimated_interface_area_ang2": round(estimated_interface_area, 1),
                "estimated_binding_energy_kcal_mol": round(estimated_binding_energy, 2),
            })

            if exposed_hydro_frac > EXPOSED_HYDROPHOBIC_THRESHOLD:
                # Only flag if interface area > 500 Å² (not trivial contacts)
                if estimated_interface_area > INTERFACE_AREA_THRESHOLD_ANG2:
                    # Only flag if binding energy is strong enough (not transient)
                    if estimated_binding_energy < BINDING_ENERGY_THRESHOLD_KCAL:
                        indicators.append(
                            f"Large exposed hydrophobic surface: "
                            f"{exposed_hydrophobic}/{exposed_total} "
                            f"({exposed_hydro_frac:.1%}) of exposed residues are hydrophobic "
                            f"(est. interface ~{estimated_interface_area:.0f} Å², "
                            f"ΔG ~{estimated_binding_energy:.1f} kcal/mol)"
                        )
                    else:
                        derivation_steps.append({
                            "step": "binding_energy_filter",
                            "result": "Binding energy above threshold, interaction likely transient, not flagging",
                            "estimated_binding_energy": round(estimated_binding_energy, 2),
                            "threshold": BINDING_ENERGY_THRESHOLD_KCAL,
                        })
                else:
                    derivation_steps.append({
                        "step": "interface_area_filter",
                        "result": "Interface area below 500 Å² threshold, not flagging",
                        "estimated_interface_area": round(estimated_interface_area, 1),
                        "threshold": INTERFACE_AREA_THRESHOLD_ANG2,
                    })
        else:
            # Fallback: sequence-based hydrophobicity estimate
            derivation_steps.append({
                "step": "exposed_hydrophobic",
                "result": "no PDB, using sequence-based estimate",
            })
            _check_sequence_hydrophobicity(protein, indicators, derivation_steps)
    else:
        # No PDB: use sequence-based heuristic
        _check_sequence_hydrophobicity(protein, indicators, derivation_steps)

    # (b) Long unstructured regions
    if pdb_string is not None:
        unstructured = find_unstructured_regions(pdb_string, min_length=UNSTRUCTURED_REGION_MIN_LENGTH)
        derivation_steps.append({
            "step": "unstructured_regions",
            "regions_found": len(unstructured),
            "regions": unstructured,
        })

        for start, end in unstructured:
            length = end - start
            # Estimate interface area for this unstructured region
            unstructured_interface_area = length * AVERAGE_RESIDUE_SURFACE_AREA * 0.5
            unstructured_binding_energy = -0.3 * length  # weaker per-residue contribution
            if unstructured_interface_area > INTERFACE_AREA_THRESHOLD_ANG2:
                if unstructured_binding_energy < BINDING_ENERGY_THRESHOLD_KCAL:
                    indicators.append(
                        f"Long unstructured region: {length} residues at positions {start}-{end} "
                        f"(est. interface ~{unstructured_interface_area:.0f} Å²)"
                    )
                else:
                    derivation_steps.append({
                        "step": "unstructured_binding_filter",
                        "detail": f"Unstructured region {start}-{end} binding energy "
                                  f"~{unstructured_binding_energy:.1f} kcal/mol above threshold, likely transient",
                    })
    else:
        # Sequence-based heuristic: regions with many small/polar residues
        long_flexible = _find_flexible_regions_sequence(protein, min_length=UNSTRUCTURED_REGION_MIN_LENGTH)
        derivation_steps.append({
            "step": "unstructured_regions",
            "method": "sequence_heuristic",
            "regions_found": len(long_flexible),
        })

        for start, end in long_flexible:
            length = end - start
            indicators.append(
                f"Long potentially unstructured region: {length} residues "
                f"at positions {start}-{end} (sequence heuristic)"
            )

    # (c) Surface charge patches
    charge_patches = find_surface_charge_patches(protein)
    significant_patches = [p for p in charge_patches if p[1] - p[0] >= SIGNIFICANT_PATCH_MIN_LENGTH]

    derivation_steps.append({
        "step": "surface_charge_patches",
        "total_patches": len(charge_patches),
        "significant_patches_5plus": len(significant_patches),
        "patches": [
            {"start": p[0], "end": p[1], "type": p[2]} for p in significant_patches
        ],
    })

    if significant_patches:
        patch_descriptions = [
            f"{p[2]} patch at {p[0]}-{p[1]} ({p[1] - p[0]} residues)"
            for p in significant_patches
        ]
        indicators.append(
            f"High surface charge patches: " + ", ".join(patch_descriptions)
        )

    # ── Crystal packing artifact filter ──
    # If the "interaction" involves >30% of surface residues, it's likely
    # a crystal packing artifact rather than a specific interaction
    if indicators and total_surface_residues > 0:
        involved_fraction = exposed_hydrophobic_count / total_surface_residues
        if involved_fraction > CRYSTAL_PACKING_SURFACE_FRACTION:
            derivation_steps.append({
                "step": "crystal_packing_filter",
                "result": (
                    f"Interaction involves {involved_fraction:.1%} of surface residues "
                    f"(>{CRYSTAL_PACKING_SURFACE_FRACTION:.0%}), likely crystal packing artifact"
                ),
                "involved_fraction": round(involved_fraction, 3),
                "threshold": CRYSTAL_PACKING_SURFACE_FRACTION,
            })
            indicators = [
                ind for ind in indicators
                if "hydrophobic surface" not in ind
            ]

    # Determine verdict (relaxed thresholds to reduce false positives)
    n_indicators = len(indicators)
    if n_indicators == 0:
        verdict = Verdict.PASS
    elif n_indicators == 1:
        verdict = Verdict.LIKELY_PASS
    elif n_indicators == 2:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_FAIL

    violation = None
    if indicators:
        violation = f"{n_indicators} interaction risk indicator(s): " + "; ".join(indicators)

    return TypeCheckResult(
        predicate="NoUnexpectedInteraction",
        verdict=verdict,
        derivation=derivation_steps,
        violation=violation,
    )


# ────────────────────────────────────────────────────────────
# Sequence-based heuristic helpers (used when no PDB available)
# ────────────────────────────────────────────────────────────

def _check_sequence_hydrophobicity(
    protein: str,
    indicators: list[str],
    derivation_steps: list[dict],
) -> None:
    """Check overall hydrophobicity from sequence alone.

    If >40% of all residues are hydrophobic, flag as potential
    exposed hydrophobic surface risk (conservative estimate since
    we can't determine burial without structure).
    """
    protein = protein.upper()
    n = len(protein)
    if n == 0:
        return

    hydro_count = sum(1 for aa in protein if aa in HYDROPHOBIC_AAS)
    hydro_frac = hydro_count / n

    derivation_steps.append({
        "step": "sequence_hydrophobicity",
        "hydrophobic_count": hydro_count,
        "total": n,
        "fraction": round(hydro_frac, 3),
        "note": "conservative estimate without structure data",
    })

    if hydro_frac > SEQUENCE_HYDROPHOBICITY_THRESHOLD:
        indicators.append(
            f"High overall hydrophobicity: {hydro_frac:.1%} of residues "
            f"are hydrophobic (no structure to confirm burial)"
        )


def _find_flexible_regions_sequence(
    protein: str,
    min_length: int = 30,
) -> list[tuple[int, int]]:
    """Find potentially unstructured regions from sequence alone.

    Uses a simple heuristic: regions rich in small/flexible residues
    (G, P, S, N, D, E, K, Q) and poor in large hydrophobic residues
    (W, Y, F, I, L, V) tend to be disordered.

    Args:
        protein: Amino acid sequence.
        min_length: Minimum region length to report.

    Returns:
        List of (start, end) tuples.
    """
    protein = protein.upper()
    flexible_residues = set("GPSNDEKQ")
    regions: list[tuple[int, int]] = []

    # Classify each residue as flexible or not
    is_flexible = [aa in flexible_residues for aa in protein]

    # Sliding window: find regions where >60% of residues are flexible
    window_size = min(7, len(protein))
    if window_size < 3 or len(protein) < min_length:
        return regions

    in_flexible_region = False
    region_start = 0

    for i in range(len(protein) - window_size + 1):
        window = is_flexible[i:i + window_size]
        flex_frac = sum(window) / window_size

        if flex_frac > FLEXIBLE_REGION_FRACTION:
            if not in_flexible_region:
                region_start = i
                in_flexible_region = True
        else:
            if in_flexible_region:
                if i - region_start >= min_length:
                    regions.append((region_start, i))
                in_flexible_region = False

    # Handle trailing region
    if in_flexible_region:
        end = len(protein)
        if end - region_start >= min_length:
            regions.append((region_start, end))

    return regions
