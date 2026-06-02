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
from typing import Optional

from ..type_system import Verdict, TypeCheckResult

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Amino acid property sets
# ────────────────────────────────────────────────────────────
HYDROPHOBIC_AAS = set("AVILMFWY")
POSITIVELY_CHARGED_AAS = set("KRH")
NEGATIVELY_CHARGED_AAS = set("DE")
CHARGED_AAS = POSITIVELY_CHARGED_AAS | NEGATIVELY_CHARGED_AAS


# ────────────────────────────────────────────────────────────
# PDB parsing helpers
# ────────────────────────────────────────────────────────────

def _parse_pdb_atoms(pdb_string: str) -> list[dict]:
    """Parse ATOM records from a PDB string.

    Returns a list of dicts, each with keys:
        serial, name, resname, chain, resseq, x, y, z, occupancy, bfactor, element
    """
    atoms = []
    for line in pdb_string.splitlines():
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        if len(line) < 54:
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
                "occupancy": float(line[54:60].strip()) if len(line) >= 60 else 1.0,
                "bfactor": float(line[60:66].strip()) if len(line) >= 66 else 0.0,
                "element": line[76:78].strip() if len(line) >= 78 else "",
            }
            atoms.append(atom)
        except (ValueError, IndexError):
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

    if n1_norm < 1e-8 or n2_norm < 1e-8:
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
    return 2.5 * (length ** 0.33)


def compute_secondary_structure_fractions(pdb_string: str) -> dict:
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
        very_high = sum(1 for s in plddt_scores if s > 90)
        confident = sum(1 for s in plddt_scores if 70 < s <= 90)
        low = sum(1 for s in plddt_scores if 50 < s <= 70)
        very_low = sum(1 for s in plddt_scores if s <= 50)
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

        if mean_plddt > 90:
            verdict = Verdict.PASS
        elif mean_plddt > 70:
            verdict = Verdict.LIKELY_PASS
        elif mean_plddt > 50:
            verdict = Verdict.UNCERTAIN
        else:
            verdict = Verdict.LIKELY_FAIL

        violation = None
        if mean_plddt < min_mean_plddt:
            violation = (
                f"Mean pLDDT {mean_plddt:.1f} is below threshold {min_mean_plddt:.1f}; "
                f"{very_low} residues ({100 * very_low / n:.1f}%) have pLDDT <= 50"
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

        if predicted_plddt > 90:
            verdict = Verdict.LIKELY_PASS
        elif predicted_plddt > 70:
            verdict = Verdict.LIKELY_PASS
        elif predicted_plddt > 50:
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
        pass

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
            if score < 70:
                consecutive_low += 1
                max_consecutive_low = max(max_consecutive_low, consecutive_low)
            else:
                consecutive_low = 0

        derivation_steps.append({
            "step": "low_confidence_regions",
            "max_consecutive_low_plddt": max_consecutive_low,
            "threshold": 10,
            "plddt_threshold": 70,
        })

        if max_consecutive_low > 10:
            risk_indicators.append(
                f"Long low-confidence region: {max_consecutive_low} consecutive "
                f"residues with pLDDT < 70"
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
        "threshold": 0.05,
    })

    if total > 0 and outlier_frac > 0.05:
        risk_indicators.append(
            f"Ramachandran outliers: {outliers}/{total} "
            f"({100 * outlier_frac:.1f}%) > 5%"
        )

    # (c) Clash score (simplified: count CA-CA contacts < 2.0 Angstroms)
    ca_coords = _extract_ca_coords(pdb_string)
    if len(ca_coords) > 1:
        clashes = 0
        for i in range(len(ca_coords)):
            for j in range(i + 2, len(ca_coords)):  # skip adjacent residues
                dist = _distance(ca_coords[i][:3], ca_coords[j][:3])
                if dist < 2.0:
                    clashes += 1

        derivation_steps.append({
            "step": "clash_score",
            "clashes": clashes,
            "num_ca": len(ca_coords),
        })

        if clashes > 0:
            risk_indicators.append(
                f"Clash score: {clashes} non-bonded CA-CA contacts < 2.0 Angstroms"
            )

    # (d) Radius of gyration check
    if ca_coords:
        rg = _compute_radius_of_gyration(ca_coords)
        expected_rg = expected_radius_of_gyration(len(protein))
        # Allow ±40% deviation from expected
        rg_low = expected_rg * 0.6
        rg_high = expected_rg * 1.4

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


def evaluate_correct_fold_topology(
    sequence: str,
    protein: str,
    organism: str,
    pdb_string: str | None = None,
) -> TypeCheckResult:
    """Evaluate that the predicted fold makes biological sense.

    Heuristics checked:
    a. Expected secondary structure content (alpha/beta/coil proportions)
       matches protein length
    b. Radius of gyration within expected range: Rg ≈ 2.5 * N^0.33
    c. Packing density in normal range (5-12 CA neighbors within 10A)
    d. Hydrophobic residues predominantly buried (if PDB available)

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

    Returns:
        TypeCheckResult with verdict based on fold topology checks.
    """
    failed_checks: list[str] = []
    borderline_checks: list[str] = []
    derivation_steps: list[dict] = []

    if pdb_string is None:
        return TypeCheckResult(
            predicate="CorrectFoldTopology",
            verdict=Verdict.UNCERTAIN,
            derivation=[{"step": "no_structure", "result": "no PDB provided"}],
            knowledge_gap="No structure provided; cannot evaluate fold topology",
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
    if n_residues > 50:
        if ss_fracs["coil"] > 0.80:
            failed_checks.append(
                f"Excessive coil content: {ss_fracs['coil']:.1%} coil "
                f"(expected <80% for structured proteins)"
            )
        elif ss_fracs["coil"] > 0.65:
            borderline_checks.append(
                f"High coil content: {ss_fracs['coil']:.1%} coil"
            )
        if ss_fracs["helix"] > 0.85:
            failed_checks.append(
                f"Unusually high helix content: {ss_fracs['helix']:.1%} helix"
            )
        elif ss_fracs["helix"] > 0.75:
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

        if rg_ratio < 0.6 or rg_ratio > 1.5:
            failed_checks.append(
                f"Radius of gyration ratio {rg_ratio:.2f} far from expected 1.0 "
                f"(Rg={rg:.1f}A vs expected {expected_rg:.1f}A)"
            )
        elif rg_ratio < 0.7 or rg_ratio > 1.4:
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
                if _distance(ca_i[:3], ca_j[:3]) <= 10.0:
                    count += 1
            neighbor_counts.append(count)

        mean_neighbors = sum(neighbor_counts) / len(neighbor_counts) if neighbor_counts else 0
        derivation_steps.append({
            "step": "packing_density",
            "mean_ca_neighbors_10A": round(mean_neighbors, 2),
            "acceptable_range": [5, 12],
        })

        if mean_neighbors < 5 or mean_neighbors > 12:
            failed_checks.append(
                f"Packing density out of range: mean {mean_neighbors:.1f} "
                f"CA neighbors within 10A (expected 5-12)"
            )
        elif mean_neighbors < 6 or mean_neighbors > 11:
            borderline_checks.append(
                f"Packing density borderline: mean {mean_neighbors:.1f} CA neighbors"
            )
    else:
        borderline_checks.append("Insufficient CA coordinates for packing analysis")

    # (d) Hydrophobic burial check (requires PDB)
    if ca_coords and n_residues > 10:
        # Compute centroid
        n_ca = len(ca_coords)
        cx = sum(c[0] for c in ca_coords) / n_ca
        cy = sum(c[1] for c in ca_coords) / n_ca
        cz = sum(c[2] for c in ca_coords) / n_ca

        # Define "buried" as within 70% of the max distance from centroid
        max_dist = max(
            math.sqrt((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2)
            for c in ca_coords
        )
        burial_cutoff = max_dist * 0.7

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
            elif exposed_hydro_frac > 0.5 and exposed_hydro_frac >= buried_hydro_frac * 0.85:
                borderline_checks.append(
                    f"Borderline hydrophobic burial: "
                    f"exposed={exposed_hydro_frac:.1%} vs buried={buried_hydro_frac:.1%}"
                )

    # Determine verdict
    n_failed = len(failed_checks)
    n_borderline = len(borderline_checks)

    if n_failed == 0 and n_borderline == 0:
        verdict = Verdict.PASS
    elif n_failed == 0 and n_borderline == 1:
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
) -> TypeCheckResult:
    """Evaluate potential for unwanted protein-protein interactions.

    Checks for:
    a. Large exposed hydrophobic surface (>40% of exposed residues are hydrophobic)
    b. Long unstructured regions (>30 residues with low secondary structure)
    c. High surface charge patches (cluster of same-charge residues)

    Verdict:
        No indicators  → PASS
        1 indicator    → UNCERTAIN
        2+ indicators  → LIKELY_FAIL

    Args:
        sequence: DNA coding sequence.
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism name.
        pdb_string: Optional PDB format string with structure data.

    Returns:
        TypeCheckResult with verdict based on interaction risk indicators.
    """
    indicators: list[str] = []
    derivation_steps: list[dict] = []

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
            # Exposed residues: beyond 70% of max distance from centroid
            exposure_cutoff = max_dist * 0.7

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

            exposed_hydro_frac = exposed_hydrophobic / exposed_total if exposed_total > 0 else 0.0

            derivation_steps.append({
                "step": "exposed_hydrophobic",
                "exposed_residues": exposed_total,
                "exposed_hydrophobic": exposed_hydrophobic,
                "fraction": round(exposed_hydro_frac, 3),
                "threshold": 0.40,
            })

            if exposed_hydro_frac > 0.40:
                indicators.append(
                    f"Large exposed hydrophobic surface: "
                    f"{exposed_hydrophobic}/{exposed_total} "
                    f"({exposed_hydro_frac:.1%}) of exposed residues are hydrophobic"
                )
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
        unstructured = find_unstructured_regions(pdb_string, min_length=30)
        derivation_steps.append({
            "step": "unstructured_regions",
            "regions_found": len(unstructured),
            "regions": unstructured,
        })

        for start, end in unstructured:
            length = end - start
            indicators.append(
                f"Long unstructured region: {length} residues at positions {start}-{end}"
            )
    else:
        # Sequence-based heuristic: regions with many small/polar residues
        long_flexible = _find_flexible_regions_sequence(protein, min_length=30)
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
    significant_patches = [p for p in charge_patches if p[1] - p[0] >= 5]

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

    # Determine verdict
    n_indicators = len(indicators)
    if n_indicators == 0:
        verdict = Verdict.PASS
    elif n_indicators == 1:
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

    if hydro_frac > 0.45:
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

        if flex_frac > 0.60:
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
