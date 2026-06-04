"""
BioCompiler Structure Quality Assessment
==========================================
Comprehensive quality metrics for predicted protein structures,
including pLDDT assessment, Ramachandran analysis, clash detection,
packing density, solvent accessibility, and hydrophobic exposure.

All metrics are computed from PDB coordinate data using simplified
geometric approaches — no external dependencies required.
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass
from typing import Any

from ..constants import HYDROPATHY, HYDROPHOBIC_AAS

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

# KYTE_DOOLITTLE: backward-compatible alias for constants.HYDROPATHY
KYTE_DOOLITTLE = HYDROPATHY

VDW_RADII: dict[str, float] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "H": 1.20,
    "P": 1.80,
    "SE": 1.80,
}

# Maximum theoretical SASA for each amino acid (Tien et al. 2013)
MAX_SASA: dict[str, float] = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "E": 223.0,
    "Q": 225.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}

# HYDROPHOBIC_AAS now imported from ..constants (standardized in v7.5.0)

# ────────────────────────────────────────────────────────────
# pLDDT thresholds
# ────────────────────────────────────────────────────────────
PLDDT_VERY_HIGH_THRESHOLD: float = 90.0
PLDDT_CONFIDENT_THRESHOLD: float = 70.0
PLDDT_LOW_THRESHOLD: float = 50.0
PLDDT_RUNNING_AVG_WINDOW: int = 5

# ────────────────────────────────────────────────────────────
# Clash detection thresholds
# ────────────────────────────────────────────────────────────
CLASH_OVERLAP_FACTOR: float = 0.4
CLASH_BONDED_SKIP_RANGE: int = 4

# ────────────────────────────────────────────────────────────
# SASA / exposure thresholds
# ────────────────────────────────────────────────────────────
SASA_NEIGHBOR_CUTOFF: float = 12.0
SASA_EXPOSURE_THRESHOLD: int = 15
SASA_DECAY_CONSTANT: float = 0.05
SASA_DEFAULT_MAX: float = 180.0

# ────────────────────────────────────────────────────────────
# Overall quality thresholds
# ────────────────────────────────────────────────────────────
EXCELLENT_PLDDT: float = 85.0
EXCELLENT_RAMA: float = 90.0
EXCELLENT_CLASH: float = 10.0
GOOD_PLDDT: float = 70.0
GOOD_RAMA: float = 80.0
GOOD_CLASH: float = 20.0
ACCEPTABLE_PLDDT: float = 50.0
ACCEPTABLE_RAMA: float = 70.0
FAIL_PLDDT_THRESHOLD: float = 30.0

# ────────────────────────────────────────────────────────────
# CA-dihedral scaling
# ────────────────────────────────────────────────────────────
CA_DIHEDRAL_SCALE: float = 1.6


# ────────────────────────────────────────────────────────────
# Data Classes
# ────────────────────────────────────────────────────────────

@dataclass
class StructureQualityReport:
    """Comprehensive quality report for a predicted protein structure.

    Attributes:
        mean_plddt: Average pLDDT confidence score across all residues.
        plddt_categories: Count of residues per confidence category.
            Keys: "very_high", "confident", "low", "very_low".
        ramachandran_favored: Percentage of residues in Ramachandran
            favored regions.
        ramachandran_allowed: Percentage of residues in Ramachandran
            allowed regions.
        ramachandran_outliers: Percentage of residues in Ramachandran
            outlier regions.
        clash_score: Steric clashes per 1000 atoms.
        molprobity_score: Approximate MolProbity score combining
            clash, Ramachandran, and rotamer metrics.
        radius_of_gyration: Radius of gyration in Angstroms.
        packing_density: Mean number of CA neighbors within cutoff.
        exposed_hydrophobic_fraction: Fraction of exposed residues
            that are hydrophobic (aggregation risk indicator).
        overall_quality: One of "excellent", "good", "acceptable", "poor".
        verdict: One of "PASS", "LIKELY_PASS", "UNCERTAIN",
            "LIKELY_FAIL", "FAIL".
    """

    mean_plddt: float
    plddt_categories: dict[str, int]
    ramachandran_favored: float
    ramachandran_allowed: float
    ramachandran_outliers: float
    clash_score: float
    molprobity_score: float
    radius_of_gyration: float
    packing_density: float
    exposed_hydrophobic_fraction: float
    overall_quality: str
    verdict: str


# ────────────────────────────────────────────────────────────
# pLDDT Assessment
# ────────────────────────────────────────────────────────────

def assess_plddt(plddt_scores: list[float]) -> dict[str, Any]:
    """Assess per-residue pLDDT confidence scores.

    Categorizes each residue into confidence bands:
      - very_high: pLDDT > 90
      - confident: 70 <= pLDDT <= 90
      - low: 50 <= pLDDT < 70
      - very_low: pLDDT < 50

    Also computes a running average over windows of 5 residues for
    local quality assessment, and identifies low-confidence regions
    (consecutive residues with pLDDT < 70).

    Args:
        plddt_scores: Per-residue pLDDT scores (0-100 scale).

    Returns:
        Dictionary with keys:
          - "counts": dict with category -> count
          - "percentages": dict with category -> percentage
          - "mean": mean pLDDT
          - "running_average": list of smoothed pLDDT values
          - "low_confidence_regions": list of (start, end) tuples
    """
    if not plddt_scores:
        return {
            "counts": {"very_high": 0, "confident": 0, "low": 0, "very_low": 0},
            "percentages": {"very_high": 0.0, "confident": 0.0, "low": 0.0, "very_low": 0.0},
            "mean": 0.0,
            "running_average": [],
            "low_confidence_regions": [],
        }

    n = len(plddt_scores)
    counts: dict[str, int] = {"very_high": 0, "confident": 0, "low": 0, "very_low": 0}

    for score in plddt_scores:
        if score > PLDDT_VERY_HIGH_THRESHOLD:
            counts["very_high"] += 1
        elif score >= PLDDT_CONFIDENT_THRESHOLD:
            counts["confident"] += 1
        elif score >= PLDDT_LOW_THRESHOLD:
            counts["low"] += 1
        else:
            counts["very_low"] += 1

    percentages = {k: (v / n) * 100.0 for k, v in counts.items()}
    mean_val = sum(plddt_scores) / n

    # Running average over configured window
    half_w = PLDDT_RUNNING_AVG_WINDOW // 2
    running_avg: list[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        running_avg.append(sum(plddt_scores[start:end]) / (end - start))

    # Identify low-confidence regions (consecutive residues < 70)
    low_conf_regions: list[tuple[int, int]] = []
    in_region = False
    region_start = 0
    for i, score in enumerate(plddt_scores):
        if score < PLDDT_CONFIDENT_THRESHOLD:
            if not in_region:
                region_start = i
                in_region = True
        else:
            if in_region:
                low_conf_regions.append((region_start, i - 1))
                in_region = False
    if in_region:
        low_conf_regions.append((region_start, n - 1))

    return {
        "counts": counts,
        "percentages": percentages,
        "mean": mean_val,
        "running_average": running_avg,
        "low_confidence_regions": low_conf_regions,
    }


# ────────────────────────────────────────────────────────────
# Ramachandran Assessment
# ────────────────────────────────────────────────────────────

def assess_ramachandran(phi_psi: list[tuple[float, float]]) -> dict[str, Any]:
    """Classify phi/psi angle pairs into Ramachandran regions.

    Favored regions (approximate):
      - alpha-helix: phi in [-100, -30], psi in [-80, -10]
      - beta-sheet: phi in [-180, -60], psi in [60, 180]
      - left-handed helix: phi in [30, 100], psi in [-10, 70]

    Allowed regions are broader bands around favored regions.

    Args:
        phi_psi: List of (phi, psi) angle pairs in degrees.

    Returns:
        Dictionary with:
          - "favored": percentage in favored regions
          - "allowed": percentage in allowed (but not favored) regions
          - "outliers": percentage in outlier regions
          - "classifications": per-residue classification list
    """
    if not phi_psi:
        return {
            "favored": 0.0,
            "allowed": 0.0,
            "outliers": 0.0,
            "classifications": [],
        }

    n = len(phi_psi)
    classifications: list[str] = []

    for phi, psi in phi_psi:
        if _in_favored_region(phi, psi):
            classifications.append("favored")
        elif _in_allowed_region(phi, psi):
            classifications.append("allowed")
        else:
            classifications.append("outlier")

    fav_count = classifications.count("favored")
    allow_count = classifications.count("allowed")
    outlier_count = classifications.count("outlier")

    return {
        "favored": (fav_count / n) * 100.0,
        "allowed": (allow_count / n) * 100.0,
        "outliers": (outlier_count / n) * 100.0,
        "classifications": classifications,
    }


def _in_favored_region(phi: float, psi: float) -> bool:
    """Check if (phi, psi) falls in a favored Ramachandran region."""
    # Alpha-helix favored
    if -100 <= phi <= -30 and -80 <= psi <= -10:
        return True
    # Beta-sheet favored
    if -180 <= phi <= -60 and 60 <= psi <= 180:
        return True
    # Left-handed helix favored
    if 30 <= phi <= 100 and -10 <= psi <= 70:
        return True
    return False


def _in_allowed_region(phi: float, psi: float) -> bool:
    """Check if (phi, psi) falls in an allowed Ramachandran region.

    Allowed regions are broader bands around the favored regions.
    """
    # Alpha-helix allowed (broader)
    if -150 <= phi <= -10 and -150 <= psi <= 50:
        return True
    # Beta-sheet allowed (broader)
    if -180 <= phi <= -30 and 30 <= psi <= 180:
        return True
    # Also allow negative psi with moderate phi (beta-region extension)
    if -180 <= phi <= -30 and -180 <= psi <= -120:
        return True
    # Left-handed helix allowed (broader)
    if 0 <= phi <= 130 and -40 <= psi <= 100:
        return True
    # Alpha-L region extension
    if 40 <= phi <= 100 and 0 <= psi <= 50:
        return True
    return False


# ────────────────────────────────────────────────────────────
# Clash Score Computation
# ────────────────────────────────────────────────────────────

def compute_clash_score(atoms: list[dict[str, Any]]) -> float:
    """Compute steric clash score (clashes per 1000 atoms).

    A steric clash is defined as a non-bonded atom pair whose distance
    is less than 0.4 times the sum of their van der Waals radii.
    Only inter-residue pairs are checked; bonded neighbors within
    4 residue sequence positions are skipped.

    Uses simplified van der Waals radii: C=1.70, N=1.55, O=1.52, S=1.80.

    Args:
        atoms: List of atom dicts, each with keys:
            - "element": str (e.g. "C", "N", "O", "S")
            - "x", "y", "z": float coordinates in Angstroms
            - "residue_index": int (0-based residue position)

    Returns:
        Clashes per 1000 atoms.
    """
    n = len(atoms)
    if n < 2:
        return 0.0

    clash_count = 0

    for i in range(n):
        elem_i = atoms[i].get("element", "C")
        r_i = VDW_RADII.get(elem_i, 1.70)
        ri_idx = atoms[i].get("residue_index", 0)
        xi = atoms[i]["x"]
        yi = atoms[i]["y"]
        zi = atoms[i]["z"]

        for j in range(i + 1, n):
            ri_jdx = atoms[j].get("residue_index", 0)

            # Skip bonded neighbors within configured range
            if abs(ri_idx - ri_jdx) <= CLASH_BONDED_SKIP_RANGE:
                continue

            elem_j = atoms[j].get("element", "C")
            r_j = VDW_RADII.get(elem_j, 1.70)

            min_dist = CLASH_OVERLAP_FACTOR * (r_i + r_j)

            dx = xi - atoms[j]["x"]
            dy = yi - atoms[j]["y"]
            dz = zi - atoms[j]["z"]
            dist_sq = dx * dx + dy * dy + dz * dz

            if dist_sq < min_dist * min_dist:
                clash_count += 1

    return (clash_count / n) * 1000.0


# ────────────────────────────────────────────────────────────
# Packing Density
# ────────────────────────────────────────────────────────────

def compute_packing_density(
    ca_coords: list[tuple[float, float, float]],
    distance_cutoff: float = 10.0,
) -> float:
    """Compute mean packing density from CA coordinates.

    For each CA atom, count the number of neighboring CA atoms
    within the distance cutoff. The mean across all residues is
    the packing density. Higher values indicate more compact structures.

    Args:
        ca_coords: List of (x, y, z) tuples for CA atoms.
        distance_cutoff: Distance threshold in Angstroms.

    Returns:
        Mean number of neighbors per residue.
    """
    n = len(ca_coords)
    if n == 0:
        return 0.0

    cutoff_sq = distance_cutoff * distance_cutoff
    total_neighbors = 0

    for i in range(n):
        xi, yi, zi = ca_coords[i]
        neighbors = 0
        for j in range(n):
            if i == j:
                continue
            dx = xi - ca_coords[j][0]
            dy = yi - ca_coords[j][1]
            dz = zi - ca_coords[j][2]
            if dx * dx + dy * dy + dz * dz <= cutoff_sq:
                neighbors += 1
        total_neighbors += neighbors

    return total_neighbors / n


# ────────────────────────────────────────────────────────────
# Exposed Hydrophobic Fraction
# ────────────────────────────────────────────────────────────

def compute_exposed_hydrophobic(
    ca_coords: list[tuple[float, float, float]],
    residues: list[str],
    probe_radius: float = 1.4,
) -> float:
    """Compute the fraction of exposed residues that are hydrophobic.

    Approximate solvent accessibility: a CA atom is considered
    "exposed" if it has fewer than 15 neighbors within 12 Angstroms.

    Hydrophobic amino acids: defined by constants.HYDROPHOBIC_AAS.

    A high value indicates aggregation risk — hydrophobic residues
    on the surface may drive non-specific interactions.

    Args:
        ca_coords: List of (x, y, z) tuples for CA atoms.
        residues: List of single-letter amino acid codes.
        probe_radius: Probe radius in Angstroms (used to adjust
            the neighbor cutoff, default 1.4).

    Returns:
        Fraction of exposed residues that are hydrophobic (0.0-1.0).
    """
    n = len(ca_coords)
    if n == 0:
        return 0.0

    # Effective neighbor cutoff incorporates the probe radius
    neighbor_cutoff = SASA_NEIGHBOR_CUTOFF - 1.4 + probe_radius
    exposure_threshold = SASA_EXPOSURE_THRESHOLD

    exposed_count = 0
    exposed_hydrophobic = 0

    for i in range(n):
        xi, yi, zi = ca_coords[i]
        neighbors = 0
        for j in range(n):
            if i == j:
                continue
            dx = xi - ca_coords[j][0]
            dy = yi - ca_coords[j][1]
            dz = zi - ca_coords[j][2]
            if dx * dx + dy * dy + dz * dz <= neighbor_cutoff * neighbor_cutoff:
                neighbors += 1

        if neighbors < exposure_threshold:
            exposed_count += 1
            if i < len(residues) and residues[i] in HYDROPHOBIC_AAS:
                exposed_hydrophobic += 1

    if exposed_count == 0:
        return 0.0

    return exposed_hydrophobic / exposed_count


# ────────────────────────────────────────────────────────────
# Low Confidence Region Detection
# ────────────────────────────────────────────────────────────

def find_low_confidence_regions(
    plddt_scores: list[float],
    window: int = 10,
    threshold: float = 70.0,
) -> list[tuple[int, int]]:
    """Find regions with mean pLDDT below threshold using sliding window.

    Scans the sequence with the specified window size. If the mean
    pLDDT in a window falls below the threshold, that window is
    flagged. Overlapping and adjacent flagged windows are merged
    into contiguous regions.

    Args:
        plddt_scores: Per-residue pLDDT scores.
        window: Sliding window size in residues.
        threshold: Mean pLDDT threshold below which regions are
            considered low confidence.

    Returns:
        List of (start, end) index tuples (inclusive) for
        low-confidence regions.
    """
    n = len(plddt_scores)
    if n == 0 or window <= 0:
        return []

    # Determine which residues are in low-confidence windows
    flagged = [False] * n
    half_w = window // 2

    for center in range(n):
        start = max(0, center - half_w)
        end = min(n, center + half_w + 1)
        if end - start < 1:
            continue
        window_mean = sum(plddt_scores[start:end]) / (end - start)
        if window_mean < threshold:
            for idx in range(start, end):
                flagged[idx] = True

    # Merge contiguous flagged regions
    regions: list[tuple[int, int]] = []
    in_region = False
    region_start = 0

    for i in range(n):
        if flagged[i]:
            if not in_region:
                region_start = i
                in_region = True
        else:
            if in_region:
                regions.append((region_start, i - 1))
                in_region = False

    if in_region:
        regions.append((region_start, n - 1))

    return regions


# ────────────────────────────────────────────────────────────
# SASA Approximation
# ────────────────────────────────────────────────────────────

def compute_sasa_approximation(
    ca_coords: list[tuple[float, float, float]],
    residue_names: list[str],
) -> list[float]:
    """Approximate relative solvent-accessible surface area per residue.

    Uses a simple neighbor-counting approach: for each CA atom, count
    the number of other CA atoms within 12 Angstroms. This count is
    mapped to an approximate absolute SASA, then normalized by the
    maximum theoretical SASA for the amino acid type (Tien et al. 2013).

    Args:
        ca_coords: List of (x, y, z) tuples for CA atoms.
        residue_names: List of single-letter amino acid codes.

    Returns:
        Per-residue relative SASA values (0.0 to 1.0).
    """
    n = len(ca_coords)
    if n == 0:
        return []

    cutoff_sq = SASA_NEIGHBOR_CUTOFF * SASA_NEIGHBOR_CUTOFF

    # Count neighbors for each CA
    neighbor_counts: list[int] = []
    for i in range(n):
        xi, yi, zi = ca_coords[i]
        count = 0
        for j in range(n):
            if i == j:
                continue
            dx = xi - ca_coords[j][0]
            dy = yi - ca_coords[j][1]
            dz = zi - ca_coords[j][2]
            if dx * dx + dy * dy + dz * dz <= cutoff_sq:
                count += 1
        neighbor_counts.append(count)

    # Map neighbor count to approximate absolute SASA
    # Using a simple decay model: max_sasa * exp(-k * neighbors)
    # where k is calibrated so that ~15 neighbors ~ 50% burial

    relative_sasa: list[float] = []
    for i in range(n):
        aa = residue_names[i] if i < len(residue_names) else "A"
        max_sasa = MAX_SASA.get(aa, SASA_DEFAULT_MAX)
        approx_abs_sasa = max_sasa * math.exp(-SASA_DECAY_CONSTANT * neighbor_counts[i])
        rel = approx_abs_sasa / max_sasa
        # Clamp to [0.0, 1.0]
        relative_sasa.append(max(0.0, min(1.0, rel)))

    return relative_sasa


# ────────────────────────────────────────────────────────────
# Radius of Gyration
# ────────────────────────────────────────────────────────────

def _compute_radius_of_gyration(
    ca_coords: list[tuple[float, float, float]],
) -> float:
    """Compute the radius of gyration from CA coordinates.

    Rg = sqrt( sum_i |r_i - r_centroid|^2 / N )

    Args:
        ca_coords: List of (x, y, z) tuples for CA atoms.

    Returns:
        Radius of gyration in Angstroms.
    """
    n = len(ca_coords)
    if n == 0:
        return 0.0

    # Compute centroid
    cx = sum(c[0] for c in ca_coords) / n
    cy = sum(c[1] for c in ca_coords) / n
    cz = sum(c[2] for c in ca_coords) / n

    # Sum of squared distances from centroid
    sum_sq = 0.0
    for x, y, z in ca_coords:
        dx = x - cx
        dy = y - cy
        dz = z - cz
        sum_sq += dx * dx + dy * dy + dz * dz

    return math.sqrt(sum_sq / n)


# ────────────────────────────────────────────────────────────
# MolProbity Score Approximation
# ────────────────────────────────────────────────────────────

def _approximate_molprobity_score(
    clash_score: float,
    ramachandran_outliers: float,
    rota_outliers: float = 0.0,
) -> float:
    """Approximate MolProbity score from component metrics.

    The MolProbity score combines clash score, Ramachandran outliers,
    and rotamer outliers into a single quality metric. Lower is better.

    Approximation: 0.5 * log10(clash_score + 1) + 0.3 * (rama_outliers / 10)
                   + 0.2 * (rota_outliers / 10)

    A score < 2.0 is excellent, 2.0-3.0 is good, > 3.0 is problematic.

    Args:
        clash_score: Steric clashes per 1000 atoms.
        ramachandran_outliers: Percentage of Ramachandran outliers.
        rota_outliers: Percentage of rotamer outliers (estimated as 0
            if not available).

    Returns:
        Approximate MolProbity score.
    """
    clash_component = 0.5 * math.log10(clash_score + 1.0) if clash_score > 0 else 0.0
    rama_component = 0.3 * (ramachandran_outliers / 10.0)
    rota_component = 0.2 * (rota_outliers / 10.0)
    return clash_component + rama_component + rota_component


# ────────────────────────────────────────────────────────────
# PDB Parsing
# ────────────────────────────────────────────────────────────

def _parse_pdb_string(pdb_string: str) -> dict[str, Any]:
    """Parse a PDB string to extract atom information.

    Simple inline parser that reads ATOM and HETATM records.
    Does not import from structure.py — self-contained.

    Args:
        pdb_string: PDB file content as a string.

    Returns:
        Dictionary with:
          - "atoms": list of atom dicts (element, x, y, z, residue_index,
            residue_name, atom_name, b_factor)
          - "ca_coords": list of (x, y, z) for CA atoms
          - "b_factors": list of B-factor values for CA atoms
          - "residue_names": list of 1-letter AA codes for CA atoms
          - "residues": list of 3-letter AA codes for CA atoms
    """
    THREE_TO_ONE: dict[str, str] = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M",  # Selenomethionine -> M
    }

    atoms: list[dict] = []
    ca_coords: list[tuple[float, float, float]] = []
    b_factors: list[float] = []
    residue_names: list[str] = []
    residue_three: list[str] = []

    seen_residues: set[tuple[str, int]] = set()
    current_res_idx = -1

    for line in pdb_string.splitlines():
        record = line[:6].strip()
        if record not in ("ATOM", "HETATM"):
            continue

        # PDB column definitions (1-based):
        # 1-6: Record name
        # 7-11: Atom serial number
        # 13-16: Atom name
        # 17: Alternate location indicator
        # 18-20: Residue name
        # 22: Chain ID
        # 23-26: Residue sequence number
        # 27: Code for insertion of residues
        # 31-38: X coordinate
        # 39-46: Y coordinate
        # 47-54: Z coordinate
        # 55-60: Occupancy
        # 61-66: Temperature factor (B-factor)

        try:
            atom_name = line[12:16].strip()
            alt_loc = line[16] if len(line) > 16 else " "
            res_name_3 = line[17:20].strip()
            chain_id = line[21] if len(line) > 21 else " "
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            b_factor = float(line[60:66].strip()) if len(line) > 60 else 0.0
        except (ValueError, IndexError):
            continue

        # Skip alternate conformations (keep only first or 'A')
        if alt_loc not in (" ", "A", ""):
            continue

        # Determine element
        element = line[76:78].strip() if len(line) > 77 else atom_name[0]
        if not element or element[0].isdigit():
            element = atom_name[0] if atom_name else "C"
        element = element.upper()

        # Map residue name to 1-letter code
        res_name_1 = THREE_TO_ONE.get(res_name_3, "X")

        # Track residue index
        res_key = (chain_id, res_seq)
        if res_key not in seen_residues:
            seen_residues.add(res_key)
            current_res_idx += 1

        atom_dict = {
            "element": element,
            "x": x,
            "y": y,
            "z": z,
            "residue_index": current_res_idx,
            "residue_name": res_name_3,
            "atom_name": atom_name,
            "b_factor": b_factor,
        }
        atoms.append(atom_dict)

        # Extract CA-specific data
        if atom_name == "CA":
            ca_coords.append((x, y, z))
            b_factors.append(b_factor)
            residue_names.append(res_name_1)
            residue_three.append(res_name_3)

    return {
        "atoms": atoms,
        "ca_coords": ca_coords,
        "b_factors": b_factors,
        "residue_names": residue_names,
        "residues": residue_three,
    }


# ────────────────────────────────────────────────────────────
# Main Quality Computation
# ────────────────────────────────────────────────────────────

def compute_structure_quality(pdb_string: str) -> StructureQualityReport:
    """Compute comprehensive structure quality metrics from a PDB string.

    Parses the PDB string inline (simple ATOM line parsing), extracts
    CA coordinates and B-factors, then computes all quality metrics.

    Overall quality determination:
      - excellent: mean_plddt > 85, ramachandran_favored > 90%, clash_score < 10
      - good: mean_plddt > 70, ramachandran_favored > 80%, clash_score < 20
      - acceptable: mean_plddt > 50, ramachandran_favored > 70%
      - poor: anything else

    Verdict mapping:
      - excellent -> PASS
      - good -> LIKELY_PASS
      - acceptable -> UNCERTAIN
      - poor -> LIKELY_FAIL (if mean_plddt < 30) or FAIL

    Args:
        pdb_string: PDB file content as a string.

    Returns:
        StructureQualityReport with all computed metrics.
    """
    parsed = _parse_pdb_string(pdb_string)

    atoms = parsed["atoms"]
    ca_coords = parsed["ca_coords"]
    b_factors = parsed["b_factors"]
    residue_names = parsed["residue_names"]

    # Handle empty PDB
    if not ca_coords:
        return StructureQualityReport(
            mean_plddt=0.0,
            plddt_categories={"very_high": 0, "confident": 0, "low": 0, "very_low": 0},
            ramachandran_favored=0.0,
            ramachandran_allowed=0.0,
            ramachandran_outliers=0.0,
            clash_score=0.0,
            molprobity_score=0.0,
            radius_of_gyration=0.0,
            packing_density=0.0,
            exposed_hydrophobic_fraction=0.0,
            overall_quality="poor",
            verdict="FAIL",
        )

    # pLDDT assessment (B-factors used as proxy for pLDDT in predicted structures)
    plddt_result = assess_plddt(b_factors)
    mean_plddt = plddt_result["mean"]
    plddt_categories = plddt_result["counts"]

    # Ramachandran assessment (approximate from CA-only data)
    # Since we only have CA coordinates, we approximate phi/psi from
    # the backbone geometry. For a more accurate assessment, full
    # backbone coordinates would be needed. We use a simplified approach.
    phi_psi = _approximate_phi_psi(ca_coords)
    rama_result = assess_ramachandran(phi_psi)

    # Clash score
    clash = compute_clash_score(atoms)

    # Packing density
    packing = compute_packing_density(ca_coords)

    # Radius of gyration
    rg = _compute_radius_of_gyration(ca_coords)

    # Exposed hydrophobic fraction
    exposed_hydro = compute_exposed_hydrophobic(ca_coords, residue_names)

    # MolProbity score approximation
    molprobity = _approximate_molprobity_score(clash, rama_result["outliers"])

    # Determine overall quality
    overall_quality = _determine_overall_quality(
        mean_plddt,
        rama_result["favored"],
        clash,
    )

    # Determine verdict
    verdict = _quality_to_verdict(overall_quality, mean_plddt)

    return StructureQualityReport(
        mean_plddt=mean_plddt,
        plddt_categories=plddt_categories,
        ramachandran_favored=rama_result["favored"],
        ramachandran_allowed=rama_result["allowed"],
        ramachandran_outliers=rama_result["outliers"],
        clash_score=clash,
        molprobity_score=molprobity,
        radius_of_gyration=rg,
        packing_density=packing,
        exposed_hydrophobic_fraction=exposed_hydro,
        overall_quality=overall_quality,
        verdict=verdict,
    )


def _determine_overall_quality(
    mean_plddt: float,
    ramachandran_favored: float,
    clash_score: float,
) -> str:
    """Determine overall quality category from key metrics.

    Args:
        mean_plddt: Mean pLDDT score.
        ramachandran_favored: Percentage in favored Ramachandran regions.
        clash_score: Steric clashes per 1000 atoms.

    Returns:
        One of "excellent", "good", "acceptable", "poor".
    """
    if mean_plddt > EXCELLENT_PLDDT and ramachandran_favored > EXCELLENT_RAMA and clash_score < EXCELLENT_CLASH:
        return "excellent"
    elif mean_plddt > GOOD_PLDDT and ramachandran_favored > GOOD_RAMA and clash_score < GOOD_CLASH:
        return "good"
    elif mean_plddt > ACCEPTABLE_PLDDT and ramachandran_favored > ACCEPTABLE_RAMA:
        return "acceptable"
    else:
        return "poor"


def _quality_to_verdict(overall_quality: str, mean_plddt: float) -> str:
    """Map overall quality to a verdict string.

    Args:
        overall_quality: One of "excellent", "good", "acceptable", "poor".
        mean_plddt: Mean pLDDT score (used to distinguish LIKELY_FAIL/FAIL).

    Returns:
        One of "PASS", "LIKELY_PASS", "UNCERTAIN", "LIKELY_FAIL", "FAIL".
    """
    if overall_quality == "excellent":
        return "PASS"
    elif overall_quality == "good":
        return "LIKELY_PASS"
    elif overall_quality == "acceptable":
        return "UNCERTAIN"
    elif overall_quality == "poor":
        if mean_plddt < FAIL_PLDDT_THRESHOLD:
            return "FAIL"
        else:
            return "LIKELY_FAIL"
    else:
        return "FAIL"


def _approximate_phi_psi(
    ca_coords: list[tuple[float, float, float]],
) -> list[tuple[float, float]]:
    """Approximate phi/psi dihedral angles from CA-only coordinates.

    This is a rough approximation since proper phi/psi calculation
    requires N, CA, and C backbone atoms. From CA coordinates alone,
    we estimate backbone geometry using virtual bond angles and
    dihedrals formed by consecutive CA atoms.

    The approximation uses 4 consecutive CA atoms to estimate each
    dihedral angle, which captures the overall backbone conformation
    but is less accurate than full-atom calculations.

    Args:
        ca_coords: List of (x, y, z) tuples for CA atoms.

    Returns:
        List of (phi, psi) angle pairs in degrees. The first and last
        residues will have approximate values since they lack flanking
        CA atoms on one side.
    """
    n = len(ca_coords)
    if n < 4:
        # Not enough residues for dihedral estimation
        # Return rough estimates based on secondary structure heuristics
        return [(0.0, 0.0)] * n

    def _dihedral(
        p0: tuple[float, float, float],
        p1: tuple[float, float, float],
        p2: tuple[float, float, float],
        p3: tuple[float, float, float],
    ) -> float:
        """Compute dihedral angle between 4 points in degrees."""
        b1 = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
        b2 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
        b3 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])

        # Normal to plane (b1, b2)
        n1 = (
            b1[1] * b2[2] - b1[2] * b2[1],
            b1[2] * b2[0] - b1[0] * b2[2],
            b1[0] * b2[1] - b1[1] * b2[0],
        )
        # Normal to plane (b2, b3)
        n2 = (
            b2[1] * b3[2] - b2[2] * b3[1],
            b2[2] * b3[0] - b2[0] * b3[2],
            b2[0] * b3[1] - b2[1] * b3[0],
        )

        # Normalize
        n1_len = math.sqrt(n1[0] ** 2 + n1[1] ** 2 + n1[2] ** 2)
        n2_len = math.sqrt(n2[0] ** 2 + n2[1] ** 2 + n2[2] ** 2)

        if n1_len < 1e-10 or n2_len < 1e-10:
            return 0.0

        n1 = (n1[0] / n1_len, n1[1] / n1_len, n1[2] / n1_len)
        n2 = (n2[0] / n2_len, n2[1] / n2_len, n2[2] / n2_len)

        # Cosine of dihedral
        cos_angle = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]
        cos_angle = max(-1.0, min(1.0, cos_angle))

        # Cross product for sign
        m1 = (
            n1[1] * n2[2] - n1[2] * n2[1],
            n1[2] * n2[0] - n1[0] * n2[2],
            n1[0] * n2[1] - n1[1] * n2[0],
        )
        # Sign from dot(m1, b2)
        sign = m1[0] * b2[0] + m1[1] * b2[1] + m1[2] * b2[2]
        if sign < 0:
            angle = -math.degrees(math.acos(cos_angle))
        else:
            angle = math.degrees(math.acos(cos_angle))

        return angle

    phi_psi: list[tuple[float, float]] = []

    for i in range(n):
        # Approximate phi from CA[i-1], CA[i], CA[i+1], CA[i+2]
        # Approximate psi from CA[i-2], CA[i-1], CA[i], CA[i+1]
        phi = 0.0
        psi = 0.0

        if i >= 1 and i + 2 < n:
            phi = _dihedral(ca_coords[i - 1], ca_coords[i], ca_coords[i + 1], ca_coords[i + 2])

        if i >= 2 and i < n:
            psi = _dihedral(ca_coords[i - 2], ca_coords[i - 1], ca_coords[i], ca_coords[min(i + 1, n - 1)])

        # Scale approximation: CA-based dihedrals are roughly 60% of
        # actual backbone dihedrals in magnitude
        phi *= CA_DIHEDRAL_SCALE
        psi *= CA_DIHEDRAL_SCALE

        phi_psi.append((phi, psi))

    return phi_psi
