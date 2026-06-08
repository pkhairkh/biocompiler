"""
BioCompiler ProLIF Integration Module
======================================
Comprehensive protein-ligand interaction fingerprint analysis using ProLIF
(Protein-Ligand Interaction Fingerprints).

This module provides structured interaction fingerprint computation,
comparison, and binding affinity estimation from protein-ligand complexes.

ProLIF detects the following interaction types:
- Hydrophobic: Non-polar contacts between hydrophobic groups
- HBAcceptor: Hydrogen bond where the protein accepts from the ligand
- HBDonor: Hydrogen bond where the protein donates to the ligand
- Cationic: Cationic group interaction (protein cation)
- Anionic: Anionic group interaction (protein anion)
- CationPi: Cation-pi interaction (cation near aromatic ring)
- PiCation: Pi-cation interaction (aromatic ring near cation)
- PiStacking: Pi-pi stacking between aromatic rings
- VdWContact: Van der Waals contact

All functions gracefully handle missing ProLIF with clear error messages,
returning empty/zero results rather than raising ImportError.

References
----------
- Bouysset C & Fiorucci S (2021) ProLIF: interaction fingerprints for
  protein-ligand complexes. J Cheminform 13:72.
  https://doi.org/10.1186/s13321-021-00548-5
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Optional dependency detection
# ────────────────────────────────────────────────────────────

try:
    import prolif as plf
    _HAS_PROLIF = True
except ImportError:
    _HAS_PROLIF = False

try:
    from rdkit import Chem
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

#: ProLIF-supported interaction types
SUPPORTED_INTERACTION_TYPES: list[str] = [
    "Hydrophobic",
    "HBAcceptor",
    "HBDonor",
    "Cationic",
    "Anionic",
    "CationPi",
    "PiCation",
    "PiStacking",
    "VdWContact",
]

#: Interaction type weights for binding affinity estimation.
#: Stronger / more specific interactions receive higher weights.
_INTERACTION_WEIGHTS: dict[str, float] = {
    "Hydrophobic": 0.5,
    "HBAcceptor": 1.5,
    "HBDonor": 1.5,
    "Cationic": 1.8,
    "Anionic": 1.8,
    "CationPi": 2.0,
    "PiCation": 2.0,
    "PiStacking": 1.8,
    "VdWContact": 0.3,
}

#: Scaling constant for binding affinity estimation (kcal/mol per weighted unit).
#: Derived from rough correlation: each weighted interaction contributes
#: approximately -0.7 kcal/mol to binding free energy.
_AFFINITY_SCALE: float = -0.7


# ────────────────────────────────────────────────────────────
# Dataclasses
# ────────────────────────────────────────────────────────────

@dataclass
class InteractionFingerprint:
    """Interaction fingerprint for a single residue.

    Attributes:
        residue_id: Unique residue identifier (e.g. "ALA45.A" for
            alanine at position 45 on chain A).
        residue_name: 3-letter amino acid code.
        interaction_types: List of interaction type strings detected
            for this residue (e.g. ["Hydrophobic", "HBAcceptor"]).
    """
    residue_id: str
    residue_name: str
    interaction_types: list[str] = field(default_factory=list)


@dataclass
class InteractionReport:
    """Comprehensive interaction report for a protein-ligand complex.

    Attributes:
        interactions: List of per-residue interaction fingerprints.
        n_total: Total number of individual interactions detected.
        n_residues_involved: Number of distinct residues with at least
            one interaction.
        binding_affinity_estimate: Rough binding affinity estimate
            in kcal/mol (negative values indicate favorable binding).
    """
    interactions: list[InteractionFingerprint] = field(default_factory=list)
    n_total: int = 0
    n_residues_involved: int = 0
    binding_affinity_estimate: float = 0.0


# ────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────

def is_prolif_available() -> bool:
    """Check whether ProLIF is installed and importable.

    Returns:
        True if the ``prolif`` package can be imported, False otherwise.
    """
    return _HAS_PROLIF


def compute_interaction_fingerprint(
    protein_pdb: str,
    ligand_pdb: str,
) -> InteractionReport:
    """Compute interaction fingerprint from PDB file paths.

    Loads protein and ligand structures from PDB files and uses ProLIF
    to detect all intermolecular interactions.

    Args:
        protein_pdb: Path to the protein PDB file.
        ligand_pdb: Path to the ligand PDB file.

    Returns:
        An :class:`InteractionReport` with per-residue fingerprints,
        total counts, and a rough binding affinity estimate.

    Raises:
        RuntimeError: If ProLIF is not installed (check with
            :func:`is_prolif_available` beforehand).
    """
    if not _HAS_PROLIF:
        raise RuntimeError(
            "ProLIF is required for interaction fingerprint computation "
            "but is not installed. Install it with: pip install prolif"
        )
    if not _HAS_RDKIT:
        raise RuntimeError(
            "RDKit is required to load PDB files for ProLIF analysis "
            "but is not installed. Install it with: pip install rdkit"
        )

    protein_mol = Chem.MolFromPDBFile(protein_pdb, removeHs=False)
    ligand_mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)

    if protein_mol is None:
        raise ValueError(
            f"Failed to load protein PDB file: {protein_pdb!r}. "
            "Ensure the file exists and contains valid ATOM records."
        )
    if ligand_mol is None:
        raise ValueError(
            f"Failed to load ligand PDB file: {ligand_pdb!r}. "
            "Ensure the file exists and contains valid HETATM/ATOM records."
        )

    return compute_interaction_fingerprint_from_mols(protein_mol, ligand_mol)


def compute_interaction_fingerprint_from_mols(
    protein_mol: Any,
    ligand_mol: Any,
) -> InteractionReport:
    """Compute interaction fingerprint from RDKit Mol objects.

    Wraps ProLIF's :class:`~prolif.Fingerprint` to detect interactions
    between the protein and ligand molecules.

    Args:
        protein_mol: RDKit Mol object for the protein (with 3D coordinates).
        ligand_mol: RDKit Mol object for the ligand (with 3D coordinates).

    Returns:
        An :class:`InteractionReport` with per-residue fingerprints,
        total counts, and a rough binding affinity estimate.

    Raises:
        RuntimeError: If ProLIF is not installed.
    """
    if not _HAS_PROLIF:
        raise RuntimeError(
            "ProLIF is required for interaction fingerprint computation "
            "but is not installed. Install it with: pip install prolif"
        )

    try:
        prot = plf.Molecule(protein_mol)
        lig = plf.Molecule(ligand_mol)

        fp = plf.Fingerprint()
        fp.run_from_iterable([lig], prot)

        ifp_df = fp.to_dataframe()

        # Aggregate interactions by residue
        residue_interactions: dict[str, dict[str, Any]] = {}
        for col in ifp_df.columns:
            residue_label, interaction_type = col
            residue_key = str(residue_label)

            if residue_key not in residue_interactions:
                # Parse residue info from ProLIF label (e.g. "ALA45.A")
                residue_name, residue_id = _parse_prolif_residue(residue_key)
                residue_interactions[residue_key] = {
                    "residue_id": residue_id,
                    "residue_name": residue_name,
                    "interaction_types": set(),
                }

            residue_interactions[residue_key]["interaction_types"].add(
                str(interaction_type)
            )

        # Build InteractionFingerprint list
        fingerprints: list[InteractionFingerprint] = []
        n_total = 0
        for info in residue_interactions.values():
            itypes = sorted(info["interaction_types"])
            fingerprints.append(InteractionFingerprint(
                residue_id=info["residue_id"],
                residue_name=info["residue_name"],
                interaction_types=itypes,
            ))
            n_total += len(itypes)

        report = InteractionReport(
            interactions=fingerprints,
            n_total=n_total,
            n_residues_involved=len(fingerprints),
            binding_affinity_estimate=0.0,
        )
        report.binding_affinity_estimate = estimate_binding_affinity(report)

        return report

    except Exception as exc:
        logger.error("ProLIF interaction fingerprint computation failed: %s", exc)
        return InteractionReport(
            interactions=[],
            n_total=0,
            n_residues_involved=0,
            binding_affinity_estimate=0.0,
        )


def detect_interaction_types(
    protein_pdb: str,
    ligand_pdb: str,
) -> dict[str, int]:
    """Count each interaction type in a protein-ligand complex.

    Args:
        protein_pdb: Path to the protein PDB file.
        ligand_pdb: Path to the ligand PDB file.

    Returns:
        Dictionary mapping interaction type names to their counts.
        Returns an empty dict if ProLIF is not available or
        computation fails.
    """
    try:
        report = compute_interaction_fingerprint(protein_pdb, ligand_pdb)
    except (RuntimeError, ValueError) as exc:
        logger.warning("Could not compute interaction fingerprint: %s", exc)
        return {}

    counts: dict[str, int] = {itype: 0 for itype in SUPPORTED_INTERACTION_TYPES}
    for fp in report.interactions:
        for itype in fp.interaction_types:
            counts[itype] = counts.get(itype, 0) + 1

    return counts


def estimate_binding_affinity(report: InteractionReport) -> float:
    """Estimate binding affinity from an interaction report.

    Uses a simple weighted-count model: each interaction type is assigned
    a weight reflecting its typical contribution to binding free energy.
    The total weighted score is scaled by :data:`_AFFINITY_SCALE` to
    produce a rough ΔG estimate in kcal/mol.

    This is an *approximation* — actual binding affinity depends on
    entropic contributions, desolvation penalties, and other factors
    not captured by a simple count-based model.

    Args:
        report: An :class:`InteractionReport` with computed interactions.

    Returns:
        Estimated binding free energy in kcal/mol (negative = favorable).
    """
    weighted_sum = 0.0
    for fp in report.interactions:
        for itype in fp.interaction_types:
            weight = _INTERACTION_WEIGHTS.get(itype, 0.5)
            weighted_sum += weight

    return weighted_sum * _AFFINITY_SCALE


def compare_interaction_patterns(
    report1: InteractionReport,
    report2: InteractionReport,
) -> float:
    """Compute Tanimoto similarity between two interaction fingerprint reports.

    Each report is converted to a binary bit vector over the space of
    (residue_id, interaction_type) pairs. The Tanimoto coefficient
    (Jaccard index) is then computed:

        T = |A ∩ B| / |A ∪ B|

    where A and B are the sets of (residue_id, interaction_type) pairs
    from each report.

    Args:
        report1: First interaction report.
        report2: Second interaction report.

    Returns:
        Tanimoto similarity in [0, 1]. Returns 0.0 if both reports
        are empty (undefined similarity).
    """
    set1 = _report_to_bit_set(report1)
    set2 = _report_to_bit_set(report2)

    if not set1 and not set2:
        return 0.0

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    if union == 0:
        return 0.0

    return intersection / union


# ────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────

def _parse_prolif_residue(label: str) -> tuple[str, str]:
    """Parse a ProLIF residue label into (residue_name, residue_id).

    ProLIF labels typically look like "ALA45.A" (name + number + chain).
    This function extracts the 3-letter code and the full identifier.

    Args:
        label: ProLIF residue label string.

    Returns:
        Tuple of (residue_name, residue_id) where residue_name is the
        3-letter amino acid code and residue_id is the full label.
    """
    # Try to extract 3-letter code from the beginning
    residue_name = label[:3] if len(label) >= 3 else label
    return residue_name, label


def _report_to_bit_set(report: InteractionReport) -> set[tuple[str, str]]:
    """Convert an InteractionReport to a set of (residue_id, interaction_type) pairs.

    Args:
        report: Interaction report to convert.

    Returns:
        Set of (residue_id, interaction_type) tuples representing all
        detected interactions in the report.
    """
    bit_set: set[tuple[str, str]] = set()
    for fp in report.interactions:
        for itype in fp.interaction_types:
            bit_set.add((fp.residue_id, itype))
    return bit_set


# ────────────────────────────────────────────────────────────
# Module-level exports
# ────────────────────────────────────────────────────────────

__all__ = [
    # Dataclasses
    "InteractionFingerprint",
    "InteractionReport",
    # Functions
    "compute_interaction_fingerprint",
    "compute_interaction_fingerprint_from_mols",
    "detect_interaction_types",
    "estimate_binding_affinity",
    "compare_interaction_patterns",
    "is_prolif_available",
    # Constants
    "SUPPORTED_INTERACTION_TYPES",
]
