"""
BioCompiler Ligand Binding Analysis (v2)
==========================================
Advanced protein-ligand binding site prediction, scoring, and analysis
module with optional RDKit, AutoDock Vina, and ProLIF integration.

This module provides:
- Pharmacophore-based binding site detection
- Per-residue energy decomposition
- SMILES feature extraction (regex fallback + RDKit)
- 3D conformer generation (RDKit ETKDGv3)
- AutoDock Vina docking wrapper
- ProLIF interaction fingerprint analysis

References
----------
- Deng W, et al. (2023) J Chem Inf Model 63:1716-1734 (pharmacophore review)
- Sink R, et al. (2023) Nat Rev Drug Discov 22:1-20 (docking review)
- Trott O & Olson AJ (2010) J Comput Chem 31:455-461 (AutoDock Vina)
- Landrum G (2023) RDKit: Open-source cheminformatics. rdkit.org
- Bouysset C & Fiorucci S (2021) J Cheminform 13:72 (ProLIF)
"""

from __future__ import annotations

import math
import re
import logging
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Optional dependency detection
# ────────────────────────────────────────────────────────────

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Lipinski, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from vina import Vina as VinaDock
    HAS_VINA = True
except ImportError:
    HAS_VINA = False

try:
    import prolif as plf
    HAS_PROLIF = True
except ImportError:
    HAS_PROLIF = False


__all__ = [
    # Data classes
    "BindingSite",
    "PharmacophoreFeature",
    "LigandInfo",
    "DockingResult",
    # Core analysis
    "detect_binding_sites",
    "score_binding_site",
    "decompose_per_residue_energy",
    # SMILES parsing
    "parse_smiles_features_rdkit",
    "_parse_smiles_features",
    # 3D conformer generation
    "generate_3d_conformer",
    # Docking
    "dock_ligand_vina",
    # Interaction fingerprints
    "compute_interaction_fingerprint",
    # Feature flags
    "HAS_RDKIT",
    "HAS_VINA",
    "HAS_PROLIF",
]


# ────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────

# H-bond distance cutoff (Angstroms)
HBOND_DISTANCE_CUTOFF: float = 3.5

# Van der Waals contact distance cutoff
VDW_CONTACT_CUTOFF: float = 4.0

# Aromatic ring centroid distance cutoff
AROMATIC_DISTANCE_CUTOFF: float = 5.5

# Halogen bond distance cutoff
HALOGEN_BOND_CUTOFF: float = 3.5

# Ionic interaction distance cutoff
IONIC_DISTANCE_CUTOFF: float = 6.0

# Default binding site search radius from any protein atom
DEFAULT_SITE_RADIUS: float = 6.0

# Minimum number of contacts to define a binding site
MIN_SITE_CONTACTS: int = 3

# Donor atom names (protein)
PROTEIN_DONOR_ATOMS: set[str] = {"N", "NE", "NH1", "NH2", "NZ", "ND1", "NE2", "OG", "OG1", "OH", "NE1"}

# Acceptor atom names (protein)
PROTEIN_ACCEPTOR_ATOMS: set[str] = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD"}

# Aromatic residue names
AROMATIC_RESIDUES: set[str] = {"PHE", "TYR", "TRP", "HIS"}

# Charged residue names (positive)
POSITIVELY_CHARGED_RESIDUES: set[str] = {"ARG", "LYS", "HIS"}

# Charged residue names (negative)
NEGATIVELY_CHARGED_RESIDUES: set[str] = {"ASP", "GLU"}

# Halogen elements
HALOGEN_ELEMENTS: set[str] = {"F", "Cl", "Br", "I"}

# Scoring weights for binding site features
SCORING_WEIGHTS: dict[str, float] = {
    "hbond": 2.0,
    "hydrophobic": 1.0,
    "aromatic": 1.5,
    "ionic": 2.5,
    "halogen": 1.8,
    "vdw": 0.5,
}


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class PharmacophoreFeature:
    """A pharmacophore feature detected in a binding site.

    Attributes:
        feature_type: Type of pharmacophore feature (hbond_donor, hbond_acceptor,
            aromatic, hydrophobic, ionic_positive, ionic_negative, halogen).
        residue_name: 3-letter amino acid code.
        chain_id: Chain identifier.
        residue_seq: Residue sequence number.
        atom_name: Atom name.
        coordinates: (x, y, z) coordinates.
        score: Feature importance score.
    """
    feature_type: str
    residue_name: str
    chain_id: str
    residue_seq: int
    atom_name: str
    coordinates: tuple[float, float, float]
    score: float = 1.0


@dataclass
class BindingSite:
    """A predicted ligand binding site.

    Attributes:
        residues: List of (chain_id, residue_seq, residue_name) tuples.
        features: List of pharmacophore features.
        score: Overall binding site score.
        center: Geometric center of the site (x, y, z).
        radius: Effective radius of the site in Angstroms.
    """
    residues: list[tuple[str, int, str]] = field(default_factory=list)
    features: list[PharmacophoreFeature] = field(default_factory=list)
    score: float = 0.0
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    radius: float = 0.0


@dataclass
class LigandInfo:
    """Information about a ligand molecule.

    Attributes:
        name: Ligand name/identifier.
        smiles: SMILES string (if available).
        atoms: List of atom dicts with keys: name, element, x, y, z.
        features: Dict of extracted features.
    """
    name: str = ""
    smiles: str = ""
    atoms: list[dict[str, Any]] = field(default_factory=list)
    features: dict[str, Any] = field(default_factory=dict)


@dataclass
class DockingResult:
    """Result of a docking calculation.

    Attributes:
        best_energy: Best (lowest) binding energy in kcal/mol.
        all_energies: All binding energies from generated poses.
        n_poses: Number of poses generated.
        error: Error message if docking failed.
        interaction_fingerprint: ProLIF interaction fingerprint dict (if available).
    """
    best_energy: float | None = None
    all_energies: list[float] = field(default_factory=list)
    n_poses: int = 0
    error: str | None = None
    interaction_fingerprint: dict | None = None


# ────────────────────────────────────────────────────────────
# Internal helper functions
# ────────────────────────────────────────────────────────────

def _distance(p1: tuple[float, ...], p2: tuple[float, ...]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def _l_is_donor(latom_element: str, latom_name: str, latom: Any = None) -> bool:
    """Check if a ligand atom is an H-bond donor.

    An atom is a donor only if it has at least one attached hydrogen.

    Args:
        latom_element: Element symbol of the ligand atom.
        latom_name: Atom name of the ligand atom.
        latom: Optional RDKit atom object for precise H-count check.

    Returns:
        True if the atom can donate a hydrogen bond.
    """
    # If RDKit atom object is available, check hydrogen count
    if latom is not None:
        try:
            if latom.GetTotalNumHs() > 0:
                return latom.GetAtomicNum() in (7, 8, 16)  # N, O, S with H
            return False
        except Exception:
            pass

    # Fallback: check typical donor patterns by atom name
    e = latom_element.upper()
    n = latom_name.upper()

    # Nitrogen donors: must have H attached
    if e == 'N':
        # Typical donor nitrogen names in PDB
        if any(x in n for x in ['NH', 'NE', 'NZ', 'ND', 'NG', 'N1', 'N2', 'N3', 'N4', 'N6']):
            return True
        # Negative: aromatic/imine nitrogen (no H)
        if any(x in n for x in ['NA', 'NB', 'NC', 'ND2']):
            return False
        # Ambiguous: assume donor for N with H likely
        return True  # Conservative

    # Oxygen donors: OH, water
    if e == 'O':
        if 'OH' in n or 'OW' in n:
            return True
        return False  # Most oxygens are acceptors, not donors

    # Sulfur donors: SH
    if e == 'S':
        if 'SH' in n or 'SG' in n:
            return True
        return False

    return False


def _p_is_acceptor(patom_element: str, patom_name: str) -> bool:
    """Check if a protein atom is a hydrogen bond acceptor.

    Args:
        patom_element: Element symbol of the protein atom.
        patom_name: Atom name of the protein atom.

    Returns:
        True if the atom can accept a hydrogen bond.
    """
    if patom_name in PROTEIN_ACCEPTOR_ATOMS:
        return True
    if patom_element in ("O", "N", "S"):
        return True
    return False


def _p_is_donor(patom_element: str, patom_name: str) -> bool:
    """Check if a protein atom is a hydrogen bond donor.

    Args:
        patom_element: Element symbol of the protein atom.
        patom_name: Atom name of the protein atom.

    Returns:
        True if the atom can donate a hydrogen bond.
    """
    if patom_name in PROTEIN_DONOR_ATOMS:
        return True
    if patom_element in ("N", "O"):
        return True
    return False


def _l_is_acceptor(latom_element: str, latom_name: str) -> bool:
    """Check if a ligand atom is a hydrogen bond acceptor.

    Args:
        latom_element: Element symbol of the ligand atom.
        latom_name: Atom name of the ligand atom.

    Returns:
        True if the atom can accept a hydrogen bond.
    """
    if latom_element in ("O", "N", "S"):
        return True
    return False


def _is_aromatic_ring(residue_name: str) -> bool:
    """Check if a residue has an aromatic ring system."""
    return residue_name in AROMATIC_RESIDUES


def _is_charged_residue(residue_name: str) -> tuple[bool, str]:
    """Check if a residue is charged and return the charge type.

    Returns:
        Tuple of (is_charged, charge_type) where charge_type is
        'positive', 'negative', or ''.
    """
    if residue_name in POSITIVELY_CHARGED_RESIDUES:
        return True, "positive"
    if residue_name in NEGATIVELY_CHARGED_RESIDUES:
        return True, "negative"
    return False, ""


def _is_hydrophobic(residue_name: str) -> bool:
    """Check if a residue is hydrophobic."""
    hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"}
    return residue_name in hydrophobic


def _parse_pdb_atoms(pdb_string: str, include_het: bool = True) -> list[dict[str, Any]]:
    """Parse ATOM/HETATM records from a PDB string.

    Returns a list of dicts with keys:
        serial, name, resname, chain, resseq, x, y, z, occupancy, bfactor, element
    """
    atoms = []
    for line in pdb_string.splitlines():
        record_type = line[0:6].strip() if len(line) >= 6 else ""
        if record_type == "ATOM" or (include_het and record_type == "HETATM"):
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


# ────────────────────────────────────────────────────────────
# Regex-based SMILES feature parsing (fallback)
# ────────────────────────────────────────────────────────────

def _parse_smiles_features(smiles: str) -> dict:
    """Parse SMILES string using regex-based heuristics.

    This is a fragile fallback when RDKit is not available.
    Provides rough estimates of molecular features.

    Args:
        smiles: SMILES string.

    Returns:
        Dict with estimated features.
    """
    if not smiles:
        return {
            "hbond_donors": 0,
            "hbond_acceptors": 0,
            "aromatic_rings": 0,
            "halogens": 0,
            "charged_groups": 0,
            "rotatable_bonds": 0,
            "molecular_weight": 0.0,
            "logp": 0.0,
        }

    # Count N-H and O-H patterns (rough donor estimate)
    donor_pattern = re.compile(r"[NO][Hh]")
    hbond_donors = len(donor_pattern.findall(smiles))

    # Count N and O atoms (rough acceptor estimate)
    acceptor_pattern = re.compile(r"[nNoO]")
    hbond_acceptors = len(acceptor_pattern.findall(smiles))

    # Count aromatic rings (lowercase letters in SMILES)
    aromatic_atoms = len(re.findall(r"[cnops]", smiles))
    # Rough estimate: ~5-6 aromatic atoms per ring
    aromatic_rings = max(0, aromatic_atoms // 5)

    # Count halogens
    halogen_pattern = re.compile(r"[Ff]|[Cc]l|[Bb]r|[Ii]")
    halogens = len(halogen_pattern.findall(smiles))

    # Count charged groups
    charge_pattern = re.compile(r"[+-]")
    charged_groups = len(charge_pattern.findall(smiles))

    # Estimate rotatable bonds (single bonds not in rings)
    # Rough: count non-ring single bonds
    rotatable_bonds = max(0, smiles.count("-") + smiles.count("/") + smiles.count("\\"))

    # Very rough MW estimate from atom count
    atom_count = len(re.findall(r"[A-Za-z]", smiles))
    molecular_weight = atom_count * 12.0  # Rough average

    # Very rough logP estimate
    logp = 0.0

    return {
        "hbond_donors": hbond_donors,
        "hbond_acceptors": hbond_acceptors,
        "aromatic_rings": aromatic_rings,
        "halogens": halogens,
        "charged_groups": charged_groups,
        "rotatable_bonds": rotatable_bonds,
        "molecular_weight": molecular_weight,
        "logp": logp,
    }


# ────────────────────────────────────────────────────────────
# RDKit SMILES feature parsing (Upgrade 1)
# ────────────────────────────────────────────────────────────

def parse_smiles_features_rdkit(smiles: str) -> dict:
    """Parse SMILES using RDKit for accurate feature extraction.

    Falls back to regex-based _parse_smiles_features() when RDKit unavailable.

    Args:
        smiles: SMILES string

    Returns:
        Dict with features: hbond_donors, hbond_acceptors, aromatic_rings,
        halogens, charged_groups, rotatable_bonds, molecular_weight, logp
    """
    if not HAS_RDKIT:
        return _parse_smiles_features(smiles)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return _parse_smiles_features(smiles)

    return {
        "hbond_donors": Lipinski.NumHDonors(mol),
        "hbond_acceptors": Lipinski.NumHAcceptors(mol),
        "aromatic_rings": Lipinski.NumAromaticRings(mol),
        "halogens": sum(1 for atom in mol.GetAtoms()
                       if atom.GetAtomicNum() in (9, 17, 35, 53)),
        "charged_groups": sum(1 for atom in mol.GetAtoms()
                            if atom.GetFormalCharge() != 0),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
    }


# ────────────────────────────────────────────────────────────
# RDKit 3D conformer generation (Upgrade 2)
# ────────────────────────────────────────────────────────────

def generate_3d_conformer(smiles: str, n_conformers: int = 10) -> np.ndarray | None:
    """Generate 3D conformers using RDKit ETKDGv3 + MMFF94.

    Args:
        smiles: SMILES string
        n_conformers: Number of conformers to generate

    Returns:
        Nx3 numpy array of lowest-energy conformer coordinates, or None
    """
    if not HAS_RDKIT:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Generate conformers using ETKDGv3
    params = AllChem.ETKDGv3()
    params.numConformers = n_conformers
    params.randomSeed = 42

    cids = AllChem.EmbedMultipleConfs(mol, params)
    if not cids:
        return None

    # Minimize each conformer with MMFF94
    energies = []
    for cid in cids:
        try:
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cid
            )
            if ff is not None:
                ff.Minimize(maxIters=500)
                energies.append((cid, ff.CalcEnergy()))
        except Exception:
            continue

    if not energies:
        # Fallback to first conformer without minimization
        conf = mol.GetConformer(cids[0])
        return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

    # Select lowest energy conformer
    best_cid = min(energies, key=lambda x: x[1])[0]
    conf = mol.GetConformer(best_cid)
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

    return coords


# ────────────────────────────────────────────────────────────
# AutoDock Vina docking wrapper (Upgrade 3)
# ────────────────────────────────────────────────────────────

def dock_ligand_vina(
    receptor_pdbqt: str,
    ligand_smiles: str,
    box_center: tuple[float, float, float] = (0, 0, 0),
    box_size: tuple[int, int, int] = (20, 20, 20),
    exhaustiveness: int = 32,
    n_poses: int = 20,
) -> dict:
    """Dock ligand against receptor using AutoDock Vina.

    Args:
        receptor_pdbqt: Path to receptor PDBQT file
        ligand_smiles: SMILES string of ligand
        box_center: Center of search box (x,y,z)
        box_size: Size of search box (nx,ny,nz) in Angstroms
        exhaustiveness: Search thoroughness (default 32)
        n_poses: Number of poses to generate

    Returns:
        Dict with best_energy, all_energies, n_poses, poses_pdbqt
    """
    if not HAS_VINA:
        return {"best_energy": None, "error": "Vina not installed"}

    try:
        # Convert SMILES to PDBQT using meeko
        try:
            from meeko import MoleculePreparation
            from rdkit import Chem
            mol = Chem.MolFromSmiles(ligand_smiles)
            if mol is not None:
                mol = Chem.AddHs(mol)
                preparator = MoleculePreparation()
                preparator.prepare(mol)
                ligand_pdbqt = preparator.write_pdbqt_file()
            else:
                return {"best_energy": None, "error": "Invalid SMILES"}
        except ImportError:
            return {"best_energy": None, "error": "meeko not installed for SMILES→PDBQT conversion"}

        v = VinaDock(sf_name='vina')
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_string(ligand_pdbqt)
        v.compute_vina_maps(center=box_center, box_size=box_size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        energies = v.energies()
        result = {
            "best_energy": energies[0][0] if len(energies) > 0 else None,
            "all_energies": [e[0] for e in energies],
            "n_poses": len(energies),
        }

        # Chain ProLIF interaction fingerprint after successful docking
        if result["best_energy"] is not None:
            try:
                from .prolif_integration import compute_interaction_fingerprint_from_mols
                # Compute IFP from docked pose
                docked_pose_sdf = v.write_pose(format='sdf')
                ifp = compute_interaction_fingerprint_from_mols(receptor_pdbqt, docked_pose_sdf)
                result["interaction_fingerprint"] = ifp
            except (ImportError, Exception):
                pass  # ProLIF not available

        return result
    except Exception as e:
        return {"best_energy": None, "error": str(e)}


# ────────────────────────────────────────────────────────────
# ProLIF interaction fingerprint (Upgrade 4)
# ────────────────────────────────────────────────────────────

def compute_interaction_fingerprint(protein_pdb: str, ligand_pdb: str) -> dict:
    """Compute protein-ligand interaction fingerprint using ProLIF.

    Args:
        protein_pdb: Path to protein PDB file or RDKit Mol
        ligand_pdb: Path to ligand PDB file or RDKit Mol

    Returns:
        Dict mapping residue_id -> set of interaction types
    """
    if not HAS_PROLIF:
        return {"error": "ProLIF not installed. pip install prolif"}

    try:
        from rdkit import Chem
        protein = Chem.MolFromPDBFile(protein_pdb, removeHs=False)
        ligand = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)

        if protein is None or ligand is None:
            return {"error": "Failed to load PDB files"}

        fp = plf.Fingerprint()
        fp.run_from_iterable([plf.Molecule(ligand)], plf.Molecule(protein))

        interactions: dict[str, set[str]] = {}
        ifp = fp.to_dataframe()
        for col in ifp.columns:
            residue, interaction = col
            if str(residue) not in interactions:
                interactions[str(residue)] = set()
            interactions[str(residue)].add(interaction)

        # Convert sets to lists for JSON serialization
        return {k: list(v) for k, v in interactions.items()}
    except Exception as e:
        return {"error": str(e)}


# ────────────────────────────────────────────────────────────
# Core binding site detection
# ────────────────────────────────────────────────────────────

def detect_binding_sites(
    pdb_string: str,
    ligand_atoms: list[dict[str, Any]] | None = None,
    ligand_smiles: str | None = None,
    site_radius: float = DEFAULT_SITE_RADIUS,
    min_contacts: int = MIN_SITE_CONTACTS,
) -> list[BindingSite]:
    """Detect ligand binding sites in a protein structure.

    If ligand atoms are provided (from HETATM records or a SMILES-derived
    conformer), finds the protein residues within site_radius of any
    ligand atom and scores them.

    If no ligand is provided, attempts to detect pockets using
    cavity-detection heuristics.

    Args:
        pdb_string: PDB format string of the protein (with or without ligand).
        ligand_atoms: Optional list of ligand atom dicts with x, y, z keys.
        ligand_smiles: Optional SMILES string for feature extraction.
        site_radius: Search radius around ligand atoms (Angstroms).
        min_contacts: Minimum contacts to define a binding site.

    Returns:
        List of BindingSite objects, sorted by score (descending).
    """
    protein_atoms = _parse_pdb_atoms(pdb_string, include_het=False)

    if not protein_atoms:
        return []

    # If no ligand provided, try to detect pockets heuristically
    if ligand_atoms is None or len(ligand_atoms) == 0:
        return _detect_pockets_heuristic(protein_atoms, site_radius, min_contacts)

    # Find protein atoms within site_radius of any ligand atom
    contact_atoms: list[dict[str, Any]] = []
    for patom in protein_atoms:
        for latom in ligand_atoms:
            dist = _distance(
                (patom["x"], patom["y"], patom["z"]),
                (latom["x"], latom["y"], latom["z"]),
            )
            if dist <= site_radius:
                contact_atoms.append(patom)
                break

    if not contact_atoms:
        return []

    # Group contact atoms by residue
    residue_map: dict[tuple[str, int, str], list[dict[str, Any]]] = {}
    for atom in contact_atoms:
        key = (atom["chain"], atom["resseq"], atom["resname"])
        if key not in residue_map:
            residue_map[key] = []
        residue_map[key].append(atom)

    # Extract pharmacophore features
    features = _extract_pharmacophore_features(
        protein_atoms, ligand_atoms, ligand_smiles, site_radius
    )

    # Group features by residue
    feature_by_residue: dict[tuple[str, int, str], list[PharmacophoreFeature]] = {}
    for feat in features:
        key = (feat.chain_id, feat.residue_seq, feat.residue_name)
        if key not in feature_by_residue:
            feature_by_residue[key] = []
        feature_by_residue[key].append(feat)

    # Build binding site
    residues = list(residue_map.keys())
    site_features = []
    for res_key in residues:
        if res_key in feature_by_residue:
            site_features.extend(feature_by_residue[res_key])

    # Compute geometric center
    if contact_atoms:
        cx = sum(a["x"] for a in contact_atoms) / len(contact_atoms)
        cy = sum(a["y"] for a in contact_atoms) / len(contact_atoms)
        cz = sum(a["z"] for a in contact_atoms) / len(contact_atoms)
    else:
        cx, cy, cz = 0.0, 0.0, 0.0

    # Compute effective radius
    max_dist = 0.0
    for a in contact_atoms:
        d = _distance((a["x"], a["y"], a["z"]), (cx, cy, cz))
        max_dist = max(max_dist, d)
    radius = max_dist + 2.0  # Add buffer

    # Score the site
    score = score_binding_site(site_features, ligand_smiles)

    site = BindingSite(
        residues=residues,
        features=site_features,
        score=score,
        center=(cx, cy, cz),
        radius=radius,
    )

    return [site]


def _detect_pockets_heuristic(
    protein_atoms: list[dict[str, Any]],
    site_radius: float,
    min_contacts: int,
) -> list[BindingSite]:
    """Detect potential binding pockets using surface curvature heuristics.

    Finds clusters of exposed residues with high feature density.

    Args:
        protein_atoms: List of protein atom dicts.
        site_radius: Search radius.
        min_contacts: Minimum contacts.

    Returns:
        List of BindingSite objects.
    """
    # Simple approach: find residues that are potential pharmacophore points
    feature_residues: list[dict[str, Any]] = []

    for atom in protein_atoms:
        resname = atom["resname"]
        atom_name = atom["name"]
        is_feature = False

        if atom_name in PROTEIN_DONOR_ATOMS or atom_name in PROTEIN_ACCEPTOR_ATOMS:
            is_feature = True
        elif _is_aromatic_ring(resname) and atom_name in ("CG", "CZ", "CD1", "CE1", "CD2", "CE2"):
            is_feature = True
        elif resname in POSITIVELY_CHARGED_RESIDUES and atom_name in ("CZ", "NH1", "NH2", "NZ", "NE"):
            is_feature = True
        elif resname in NEGATIVELY_CHARGED_RESIDUES and atom_name in ("CG", "OD1", "OD2", "CD", "OE1", "OE2"):
            is_feature = True

        if is_feature:
            feature_residues.append(atom)

    if not feature_residues:
        return []

    # Cluster feature residues by spatial proximity
    clusters = _cluster_atoms_spatial(feature_residues, site_radius)

    sites = []
    for cluster in clusters:
        if len(cluster) < min_contacts:
            continue

        residue_map: dict[tuple[str, int, str], list[dict[str, Any]]] = {}
        for atom in cluster:
            key = (atom["chain"], atom["resseq"], atom["resname"])
            if key not in residue_map:
                residue_map[key] = []
            residue_map[key].append(atom)

        residues = list(residue_map.keys())

        # Build features for this cluster
        features = []
        for atom in cluster:
            feat_type = _classify_atom_feature(atom)
            if feat_type:
                features.append(PharmacophoreFeature(
                    feature_type=feat_type,
                    residue_name=atom["resname"],
                    chain_id=atom["chain"],
                    residue_seq=atom["resseq"],
                    atom_name=atom["name"],
                    coordinates=(atom["x"], atom["y"], atom["z"]),
                ))

        # Geometric center
        cx = sum(a["x"] for a in cluster) / len(cluster)
        cy = sum(a["y"] for a in cluster) / len(cluster)
        cz = sum(a["z"] for a in cluster) / len(cluster)

        max_dist = max(_distance((a["x"], a["y"], a["z"]), (cx, cy, cz)) for a in cluster)
        radius = max_dist + 2.0

        score = score_binding_site(features)

        sites.append(BindingSite(
            residues=residues,
            features=features,
            score=score,
            center=(cx, cy, cz),
            radius=radius,
        ))

    # Sort by score descending
    sites.sort(key=lambda s: s.score, reverse=True)
    return sites


def _classify_atom_feature(atom: dict[str, Any]) -> str | None:
    """Classify a protein atom into a pharmacophore feature type.

    Args:
        atom: Atom dict with keys: name, resname, element.

    Returns:
        Feature type string or None.
    """
    name = atom["name"]
    resname = atom["resname"]

    if name in PROTEIN_DONOR_ATOMS:
        return "hbond_donor"
    if name in PROTEIN_ACCEPTOR_ATOMS:
        return "hbond_acceptor"
    if _is_aromatic_ring(resname):
        return "aromatic"
    is_charged, charge_type = _is_charged_residue(resname)
    if is_charged:
        return f"ionic_{charge_type}"
    if _is_hydrophobic(resname):
        return "hydrophobic"
    return None


def _cluster_atoms_spatial(
    atoms: list[dict[str, Any]],
    radius: float,
) -> list[list[dict[str, Any]]]:
    """Simple spatial clustering of atoms using distance-based grouping.

    Uses single-linkage clustering: atoms within radius of any
    cluster member are added to that cluster.

    Args:
        atoms: List of atom dicts with x, y, z keys.
        radius: Clustering distance cutoff.

    Returns:
        List of clusters (each cluster is a list of atom dicts).
    """
    if not atoms:
        return []

    # Assign each atom to a cluster via union-find
    n = len(atoms)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            dist = _distance(
                (atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]),
                (atoms[j]["x"], atoms[j]["y"], atoms[j]["z"]),
            )
            if dist <= radius:
                union(i, j)

    # Group by root
    clusters: dict[int, list[dict[str, Any]]] = {}
    for i in range(n):
        root = find(i)
        if root not in clusters:
            clusters[root] = []
        clusters[root].append(atoms[i])

    return list(clusters.values())


def _extract_pharmacophore_features(
    protein_atoms: list[dict[str, Any]],
    ligand_atoms: list[dict[str, Any]],
    ligand_smiles: str | None,
    site_radius: float,
) -> list[PharmacophoreFeature]:
    """Extract pharmacophore features at the protein-ligand interface.

    Args:
        protein_atoms: List of protein atom dicts.
        ligand_atoms: List of ligand atom dicts.
        ligand_smiles: Optional SMILES string for ligand feature extraction.
        site_radius: Distance cutoff for interface detection.

    Returns:
        List of PharmacophoreFeature objects.
    """
    features: list[PharmacophoreFeature] = []

    # Get ligand features if SMILES available
    ligand_features: dict[str, Any] = {}
    if ligand_smiles:
        ligand_features = parse_smiles_features_rdkit(ligand_smiles)

    for patom in protein_atoms:
        patom_coord = (patom["x"], patom["y"], patom["z"])
        patom_name = patom["name"]
        patom_element = patom.get("element", "")
        resname = patom["resname"]

        for latom in ligand_atoms:
            latom_coord = (latom["x"], latom["y"], latom["z"])
            latom_element = latom.get("element", "")
            latom_name = latom.get("name", "")

            dist = _distance(patom_coord, latom_coord)

            if dist > site_radius:
                continue

            # Check H-bond donor/acceptor interactions
            # Bug fix: use _l_is_donor and _p_is_acceptor with correct signatures
            if dist <= HBOND_DISTANCE_CUTOFF:
                # Ligand donor → Protein acceptor
                if _l_is_donor(latom_element, latom_name) and _p_is_acceptor(patom_element, patom_name):
                    features.append(PharmacophoreFeature(
                        feature_type="hbond",
                        residue_name=resname,
                        chain_id=patom["chain"],
                        residue_seq=patom["resseq"],
                        atom_name=patom_name,
                        coordinates=patom_coord,
                        score=SCORING_WEIGHTS["hbond"] * (1.0 - dist / HBOND_DISTANCE_CUTOFF),
                    ))
                # Protein donor → Ligand acceptor
                if _p_is_donor(patom_element, patom_name) and _l_is_acceptor(latom_element, latom_name):
                    features.append(PharmacophoreFeature(
                        feature_type="hbond",
                        residue_name=resname,
                        chain_id=patom["chain"],
                        residue_seq=patom["resseq"],
                        atom_name=patom_name,
                        coordinates=patom_coord,
                        score=SCORING_WEIGHTS["hbond"] * (1.0 - dist / HBOND_DISTANCE_CUTOFF),
                    ))

            # Aromatic stacking
            if dist <= AROMATIC_DISTANCE_CUTOFF and _is_aromatic_ring(resname):
                if latom_element in ("C", "N") and latom.get("is_aromatic", False):
                    features.append(PharmacophoreFeature(
                        feature_type="aromatic",
                        residue_name=resname,
                        chain_id=patom["chain"],
                        residue_seq=patom["resseq"],
                        atom_name=patom_name,
                        coordinates=patom_coord,
                        score=SCORING_WEIGHTS["aromatic"] * (1.0 - dist / AROMATIC_DISTANCE_CUTOFF),
                    ))

            # Ionic interactions
            if dist <= IONIC_DISTANCE_CUTOFF:
                is_charged, charge_type = _is_charged_residue(resname)
                if is_charged:
                    # Check for complementary charge on ligand
                    ligand_charged = ligand_features.get("charged_groups", 0) > 0
                    if ligand_charged:
                        features.append(PharmacophoreFeature(
                            feature_type="ionic",
                            residue_name=resname,
                            chain_id=patom["chain"],
                            residue_seq=patom["resseq"],
                            atom_name=patom_name,
                            coordinates=patom_coord,
                            score=SCORING_WEIGHTS["ionic"] * (1.0 - dist / IONIC_DISTANCE_CUTOFF),
                        ))

            # Halogen bonds
            if dist <= HALOGEN_BOND_CUTOFF and latom_element in ("F", "Cl", "Br", "I"):
                features.append(PharmacophoreFeature(
                    feature_type="halogen",
                    residue_name=resname,
                    chain_id=patom["chain"],
                    residue_seq=patom["resseq"],
                    atom_name=patom_name,
                    coordinates=patom_coord,
                    score=SCORING_WEIGHTS["halogen"] * (1.0 - dist / HALOGEN_BOND_CUTOFF),
                ))

            # Hydrophobic contacts
            if dist <= VDW_CONTACT_CUTOFF and _is_hydrophobic(resname):
                features.append(PharmacophoreFeature(
                    feature_type="hydrophobic",
                    residue_name=resname,
                    chain_id=patom["chain"],
                    residue_seq=patom["resseq"],
                    atom_name=patom_name,
                    coordinates=patom_coord,
                    score=SCORING_WEIGHTS["hydrophobic"] * (1.0 - dist / VDW_CONTACT_CUTOFF),
                ))

    return features


# ────────────────────────────────────────────────────────────
# Binding site scoring
# ────────────────────────────────────────────────────────────

def score_binding_site(
    features: list[PharmacophoreFeature],
    ligand_smiles: str | None = None,
) -> float:
    """Score a binding site based on its pharmacophore features.

    Combines feature scores with ligand complementarity.

    Args:
        features: List of pharmacophore features at the site.
        ligand_smiles: Optional SMILES string for ligand complementarity.

    Returns:
        Binding site score (higher is better).
    """
    if not features:
        return 0.0

    # Sum weighted feature scores
    total_score = 0.0
    feature_type_counts: dict[str, int] = {}

    for feat in features:
        total_score += feat.score
        feature_type_counts[feat.feature_type] = feature_type_counts.get(feat.feature_type, 0) + 1

    # Diversity bonus: having multiple feature types is better
    n_feature_types = len(feature_type_counts)
    diversity_bonus = math.log1p(n_feature_types) * 0.5

    # Ligand complementarity bonus
    complementarity_bonus = 0.0
    if ligand_smiles:
        ligand_features = parse_smiles_features_rdkit(ligand_smiles)
        if ligand_features.get("hbond_donors", 0) > 0 and "hbond" in feature_type_counts:
            complementarity_bonus += 1.0
        if ligand_features.get("hbond_acceptors", 0) > 0 and "hbond" in feature_type_counts:
            complementarity_bonus += 1.0
        if ligand_features.get("aromatic_rings", 0) > 0 and "aromatic" in feature_type_counts:
            complementarity_bonus += 1.5
        if ligand_features.get("charged_groups", 0) > 0 and "ionic" in feature_type_counts:
            complementarity_bonus += 2.0
        if ligand_features.get("halogens", 0) > 0 and "halogen" in feature_type_counts:
            complementarity_bonus += 1.0

    return total_score + diversity_bonus + complementarity_bonus


# ────────────────────────────────────────────────────────────
# Per-residue energy decomposition
# ────────────────────────────────────────────────────────────

def decompose_per_residue_energy(
    pdb_string: str,
    ligand_atoms: list[dict[str, Any]],
    ligand_smiles: str | None = None,
) -> dict[tuple[str, int, str], dict[str, float]]:
    """Decompose binding energy into per-residue contributions.

    Uses a simplified force-field approach with empirical scoring
    functions for hydrogen bonds, van der Waals contacts, ionic
    interactions, aromatic stacking, and hydrophobic effects.

    Args:
        pdb_string: PDB format string of the protein.
        ligand_atoms: List of ligand atom dicts with x, y, z, element keys.
        ligand_smiles: Optional SMILES string for feature extraction.

    Returns:
        Dict mapping (chain_id, residue_seq, residue_name) →
        dict with energy components: total, hbond, vdw, ionic, aromatic, hydrophobic.
    """
    protein_atoms = _parse_pdb_atoms(pdb_string, include_het=False)

    if not protein_atoms or not ligand_atoms:
        return {}

    # Initialize per-residue energy accumulators
    energy_map: dict[tuple[str, int, str], dict[str, float]] = {}

    # Bug fix: Use list[str] = [] instead of list[str] = set()
    # (set() is not compatible with list[str] type annotation)
    features: list[str] = []

    for patom in protein_atoms:
        patom_coord = (patom["x"], patom["y"], patom["z"])
        patom_name = patom["name"]
        patom_element = patom.get("element", "")
        resname = patom["resname"]
        resseq = patom["resseq"]
        chain = patom["chain"]
        res_key = (chain, resseq, resname)

        if res_key not in energy_map:
            energy_map[res_key] = {
                "total": 0.0,
                "hbond": 0.0,
                "vdw": 0.0,
                "ionic": 0.0,
                "aromatic": 0.0,
                "hydrophobic": 0.0,
            }

        for latom in ligand_atoms:
            latom_coord = (latom["x"], latom["y"], latom["z"])
            latom_element = latom.get("element", "")
            latom_name = latom.get("name", "")

            dist = _distance(patom_coord, latom_coord)

            if dist > VDW_CONTACT_CUTOFF * 1.5:
                continue

            # Hydrogen bond energy (distance-dependent)
            if dist <= HBOND_DISTANCE_CUTOFF:
                # Bug fix: use (element, name) argument order consistently
                # _l_is_donor(latom_element, latom_name) and _p_is_acceptor(patom_element, patom_name)
                if _l_is_donor(latom_element, latom_name) and _p_is_acceptor(patom_element, patom_name):
                    energy = _hbond_energy(dist)
                    energy_map[res_key]["hbond"] += energy
                    if "hbond_donor" not in features:
                        features.append("hbond_donor")

                if _p_is_donor(patom_element, patom_name) and _l_is_acceptor(latom_element, latom_name):
                    energy = _hbond_energy(dist)
                    energy_map[res_key]["hbond"] += energy
                    if "hbond_acceptor" not in features:
                        features.append("hbond_acceptor")

            # Van der Waals energy
            if dist <= VDW_CONTACT_CUTOFF and dist > 1.0:
                energy = _vdw_energy(dist, patom_element, latom_element)
                energy_map[res_key]["vdw"] += energy

            # Ionic interaction energy
            if dist <= IONIC_DISTANCE_CUTOFF:
                is_charged, charge_type = _is_charged_residue(resname)
                if is_charged:
                    energy = _ionic_energy(dist, charge_type)
                    energy_map[res_key]["ionic"] += energy
                    if "ionic" not in features:
                        features.append("ionic")

            # Aromatic stacking energy
            if dist <= AROMATIC_DISTANCE_CUTOFF and _is_aromatic_ring(resname):
                energy = _aromatic_energy(dist)
                energy_map[res_key]["aromatic"] += energy
                if "aromatic" not in features:
                    features.append("aromatic")

            # Hydrophobic contact energy
            if dist <= VDW_CONTACT_CUTOFF and _is_hydrophobic(resname):
                energy = _hydrophobic_energy(dist)
                energy_map[res_key]["hydrophobic"] += energy
                if "hydrophobic" not in features:
                    features.append("hydrophobic")

    # Compute total per-residue energies
    for res_key, energies in energy_map.items():
        energies["total"] = (
            energies["hbond"]
            + energies["vdw"]
            + energies["ionic"]
            + energies["aromatic"]
            + energies["hydrophobic"]
        )

    return energy_map


# ────────────────────────────────────────────────────────────
# Empirical energy functions
# ────────────────────────────────────────────────────────────

def _hbond_energy(distance: float) -> float:
    """Compute hydrogen bond energy using a 12-10 potential.

    E = 5 * (r0/r)^12 - 6 * (r0/r)^10
    where r0 = 2.9 Angstroms (optimal H-bond distance).

    Args:
        distance: Distance in Angstroms.

    Returns:
        Energy in kcal/mol (negative = favorable).
    """
    r0 = 2.9
    if distance < 0.5:
        distance = 0.5  # Prevent numerical overflow
    r_ratio = r0 / distance
    energy = 5.0 * (r_ratio ** 12) - 6.0 * (r_ratio ** 10)
    return -abs(energy)  # H-bonds are favorable


def _vdw_energy(distance: float, element1: str, element2: str) -> float:
    """Compute van der Waals energy using a 6-12 Lennard-Jones potential.

    Args:
        distance: Distance in Angstroms.
        element1: Element symbol of first atom.
        element2: Element symbol of second atom.

    Returns:
        Energy in kcal/mol (negative = favorable at equilibrium).
    """
    # Simplified sigma values by element
    sigma_map = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "H": 1.2}
    sigma1 = sigma_map.get(element1, 1.7)
    sigma2 = sigma_map.get(element2, 1.7)
    sigma = (sigma1 + sigma2) / 2.0
    epsilon = 0.1  # kcal/mol (simplified)

    if distance < 0.5:
        distance = 0.5

    r_ratio = sigma / distance
    energy = 4.0 * epsilon * (r_ratio ** 12 - r_ratio ** 6)
    return energy


def _ionic_energy(distance: float, charge_type: str) -> float:
    """Compute ionic interaction energy.

    Uses Coulomb's law with a distance-dependent dielectric constant.

    Args:
        distance: Distance in Angstroms.
        charge_type: 'positive' or 'negative'.

    Returns:
        Energy in kcal/mol (negative = favorable for complementary charges).
    """
    if distance < 1.0:
        distance = 1.0

    # Distance-dependent dielectric: epsilon = 4r
    dielectric = 4.0 * distance
    # Coulomb constant in kcal*Angstrom/(mol*e^2)
    k = 332.0

    # Assuming opposite charges (favorable)
    energy = -k / (dielectric * distance)
    return energy


def _aromatic_energy(distance: float) -> float:
    """Compute aromatic stacking energy.

    Simplified model: favorable at ~3.5-5.0 Angstroms.

    Args:
        distance: Distance in Angstroms.

    Returns:
        Energy in kcal/mol (negative = favorable).
    """
    r0 = 4.5  # Optimal stacking distance
    sigma = 0.5
    energy = -1.5 * math.exp(-0.5 * ((distance - r0) / sigma) ** 2)
    return energy


def _hydrophobic_energy(distance: float) -> float:
    """Compute hydrophobic contact energy.

    Simplified: short-range favorable interaction.

    Args:
        distance: Distance in Angstroms.

    Returns:
        Energy in kcal/mol (negative = favorable).
    """
    r0 = 3.5
    if distance > VDW_CONTACT_CUTOFF:
        return 0.0
    energy = -0.3 * (1.0 - (distance - r0) / (VDW_CONTACT_CUTOFF - r0))
    return min(0.0, energy)
