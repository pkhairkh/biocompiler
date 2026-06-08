"""
BioCompiler RDKit Integration Module
======================================
RDKit-based cheminformatics functions for ligand binding analysis,
with graceful fallback when RDKit is not available.

This module provides:
- SMILES parsing to RDKit Mol objects
- Molecular descriptor computation (MW, LogP, TPSA, H-bond donors/acceptors,
  rotatable bonds, aromatic rings, fraction sp3)
- 3D conformer generation with energy ranking (ETKDGv3 / ETKDGv2 / KDG)
- Partial charge computation (Gasteiger)
- 3D pharmacophore feature extraction (HBD, HBA, aromatic, hydrophobic,
  positive, negative ionizable)
- RMSD computation between conformers
- Conformer clustering by RMSD using DBSCAN

All functions return ``None`` or empty results when RDKit is not installed,
so calling code can safely degrade without import errors.

References
----------
- Landrum G (2023) RDKit: Open-source cheminformatics.
  https://www.rdkit.org
- Riniker S, Landrum GA (2015) Better Informed Distance Geometry:
  Using What We Know To Improve Conformation Generation.
  J Chem Inf Model 55(12):2562-2574. doi:10.1021/acs.jcim.5b00654
- Wang R, Gao Y, Lai L (2000) Calculating partition coefficient by
  atom-additive method. Perspect Drug Discov Des 19:47-66.
- Gasteiger J, Marsili M (1980) Iterative partial equalization of
  orbital electronegativity. Tetrahedron 36(22):3219-3228.
- Ester M, Kriegel HP, Sander J, Xu X (1996) A density-based algorithm
  for discovering clusters in large spatial databases with noise.
  KDD-96:226-231.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional dependency detection
# ---------------------------------------------------------------------------

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig

    HAS_RDKIT: bool = True
except ImportError:
    HAS_RDKIT: bool = False

# Try to import numpy for coordinate arrays (optional but preferred)
try:
    import numpy as np

    HAS_NUMPY: bool = True
except ImportError:
    HAS_NUMPY: bool = False


# ---------------------------------------------------------------------------
# Public API flag
# ---------------------------------------------------------------------------

__all__ = [
    # Dataclasses
    "ConformerInfo",
    "PharmacophoreFeature3D",
    # Functions
    "is_rdkit_available",
    "parse_smiles",
    "compute_molecular_descriptors",
    "generate_conformers",
    "compute_partial_charges",
    "compute_rotatable_bonds",
    "extract_pharmacophore_features_3d",
    "compute_rmsd",
    "cluster_poses_by_rmsd",
]


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ConformerInfo:
    """Information about a single 3D conformer.

    Attributes:
        conformer_id: Conformer index within the RDKit molecule.
        energy: MMFF94 minimised energy (kcal/mol).  ``None`` if
            minimisation was not performed or failed.
        coordinates: Nx3 array of atomic coordinates (Angstroms).
            Stored as a list of (x, y, z) tuples when numpy is
            unavailable, or as an ndarray when numpy is present.
    """

    conformer_id: int
    energy: float | None = None
    coordinates: Any = field(default_factory=list)


@dataclass
class PharmacophoreFeature3D:
    """A 3D pharmacophore feature extracted from a ligand conformer.

    Attributes:
        feature_type: One of ``'HBD'``, ``'HBA'``, ``'aromatic'``,
            ``'hydrophobic'``, ``'positive'``, ``'negative'``.
        position: (x, y, z) centroid of the feature in Angstroms.
        direction: Optional direction vector (dx, dy, dz) for
            directional features such as H-bond donors/acceptors.
            ``(0.0, 0.0, 0.0)`` when direction is undefined.
        atom_indices: Indices of the atoms that define this feature.
    """

    feature_type: str
    position: tuple[float, float, float] = (0.0, 0.0, 0.0)
    direction: tuple[float, float, float] = (0.0, 0.0, 0.0)
    atom_indices: list[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

def is_rdkit_available() -> bool:
    """Check whether RDKit is installed and importable.

    Returns:
        ``True`` if ``from rdkit import Chem`` succeeds, ``False`` otherwise.
    """
    return HAS_RDKIT


# ---------------------------------------------------------------------------
# SMILES parsing
# ---------------------------------------------------------------------------

def parse_smiles(smiles: str) -> Any:
    """Parse a SMILES string into an RDKit Mol object.

    Args:
        smiles: A valid SMILES string (e.g. ``'c1ccccc1'``).

    Returns:
        An ``rdkit.Chem.rdchem.Mol`` instance, or ``None`` if RDKit is
        unavailable or the SMILES cannot be parsed.
    """
    if not HAS_RDKIT:
        logger.debug("RDKit not available; parse_smiles returning None")
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Failed to parse SMILES: %s", smiles)
    return mol


# ---------------------------------------------------------------------------
# Molecular descriptors
# ---------------------------------------------------------------------------

def compute_molecular_descriptors(smiles: str) -> dict[str, Any]:
    """Compute a standard set of drug-likeness molecular descriptors.

    Descriptors computed (when RDKit is available):
        - **molecular_weight**: Exact molecular weight (Da)
        - **logp**: Wildman-Crippen LogP (XLogP3-like)
        - **tpsa**: Topological polar surface area (Å²)
        - **hbond_donors**: Number of hydrogen-bond donors (Lipinski)
        - **hbond_acceptors**: Number of hydrogen-bond acceptors (Lipinski)
        - **rotatable_bonds**: Number of rotatable bonds
        - **aromatic_rings**: Number of aromatic rings
        - **fraction_sp3**: Fraction of sp³-hybridised carbons (Fsp³)

    When RDKit is unavailable all values are ``None``.

    Args:
        smiles: SMILES string.

    Returns:
        Dictionary mapping descriptor names to their computed values.
    """
    empty: dict[str, Any] = {
        "molecular_weight": None,
        "logp": None,
        "tpsa": None,
        "hbond_donors": None,
        "hbond_acceptors": None,
        "rotatable_bonds": None,
        "aromatic_rings": None,
        "fraction_sp3": None,
    }

    if not HAS_RDKIT:
        return empty

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return empty

    return {
        "molecular_weight": Descriptors.ExactMolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "hbond_donors": Lipinski.NumHDonors(mol),
        "hbond_acceptors": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "aromatic_rings": Lipinski.NumAromaticRings(mol),
        "fraction_sp3": Lipinski.FractionCSP3(mol),
    }


# ---------------------------------------------------------------------------
# 3D conformer generation
# ---------------------------------------------------------------------------

# Mapping from user-facing method names to RDKit embedding parameter classes
_METHOD_MAP: dict[str, str] = {
    "ETKDGv3": "ETKDGv3",
    "ETKDGv2": "ETKDGv2",
    "ETKDG": "ETKDGv3",   # alias – default to v3
    "KDG": "KDG",
}


def generate_conformers(
    smiles: str,
    n_conformers: int = 10,
    method: str = "ETKDGv3",
) -> list[dict]:
    """Generate 3D conformers with MMFF94 energy ranking.

    Uses RDKit's distance-geometry embedding followed by MMFF94
    force-field minimisation.  Conformers are returned sorted by
    ascending energy.

    Args:
        smiles: SMILES string.
        n_conformers: Number of conformers to generate (default 10).
        method: Embedding method — ``'ETKDGv3'`` (default, Riniker &
            Landrum 2015), ``'ETKDGv2'``, or ``'KDG'``.

    Returns:
        List of dicts, each with keys:
            - ``conformer_id`` (int): RDKit conformer index
            - ``energy`` (float | None): MMFF94 energy in kcal/mol
            - ``coordinates`` (list[tuple[float,float,float]]): Nx3 coords
        Returns an empty list when RDKit is unavailable or the SMILES
        is invalid.
    """
    if not HAS_RDKIT:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    mol = Chem.AddHs(mol)

    # Select embedding parameters
    method_key = _METHOD_MAP.get(method, "ETKDGv3")
    if method_key == "ETKDGv3":
        params = AllChem.ETKDGv3()
    elif method_key == "ETKDGv2":
        params = AllChem.ETKDGv2()
    elif method_key == "KDG":
        params = AllChem.KDG()
    else:
        params = AllChem.ETKDGv3()

    params.numConformers = n_conformers
    params.randomSeed = 42

    cids = AllChem.EmbedMultipleConfs(mol, params)
    if not cids:
        logger.warning("EmbedMultipleConfs returned no conformers for: %s", smiles)
        return []

    # Minimise each conformer with MMFF94 and collect energies
    results: list[dict] = []
    for cid in cids:
        energy: float | None = None
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                if ff is not None:
                    ff.Minimize(maxIters=500)
                    energy = ff.CalcEnergy()
        except Exception:
            logger.debug("MMFF minimisation failed for conformer %d", cid)

        conf = mol.GetConformer(cid)
        n_atoms = mol.GetNumAtoms()
        coords = [tuple(conf.GetAtomPosition(i)) for i in range(n_atoms)]

        results.append({
            "conformer_id": cid,
            "energy": energy,
            "coordinates": coords,
        })

    # Sort by energy (None values last)
    results.sort(key=lambda r: (r["energy"] is None, r["energy"] or 0.0))

    return results


# ---------------------------------------------------------------------------
# Partial charges
# ---------------------------------------------------------------------------

def compute_partial_charges(
    smiles: str,
    method: str = "Gasteiger",
) -> list[float]:
    """Compute partial atomic charges.

    Currently only the **Gasteiger** method is supported, which is
    fast and does not require a 3D geometry.

    Args:
        smiles: SMILES string.
        method: Charge computation method (default ``'Gasteiger'``).

    Returns:
        List of partial charges, one per heavy atom.  Empty list when
        RDKit is unavailable or the SMILES is invalid.
    """
    if not HAS_RDKIT:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    method_lower = method.lower()
    if method_lower == "gasteiger":
        AllChem.ComputeGasteigerCharges(mol)
        charges = []
        for atom in mol.GetAtoms():
            try:
                val = atom.GetDoubleProp("_GasteigerCharge")
                # RDKit sometimes assigns NaN for problematic atoms
                if math.isnan(val):
                    charges.append(0.0)
                else:
                    charges.append(val)
            except Exception:
                charges.append(0.0)
        return charges
    else:
        logger.warning("Unsupported charge method '%s'; returning empty list", method)
        return []


# ---------------------------------------------------------------------------
# Rotatable bonds
# ---------------------------------------------------------------------------

def compute_rotatable_bonds(smiles: str) -> int:
    """Count the number of rotatable bonds (Lipinski definition).

    Args:
        smiles: SMILES string.

    Returns:
        Number of rotatable bonds, or ``0`` if RDKit is unavailable
        or the SMILES is invalid.
    """
    if not HAS_RDKIT:
        return 0

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0

    return Lipinski.NumRotatableBonds(mol)


# ---------------------------------------------------------------------------
# 3D pharmacophore feature extraction
# ---------------------------------------------------------------------------

def extract_pharmacophore_features_3d(
    smiles: str,
    conformer_id: int = 0,
) -> list[dict]:
    """Extract 3D pharmacophore features from a ligand conformer.

    Feature types extracted:
        - **HBD** — hydrogen bond donor
        - **HBA** — hydrogen bond acceptor
        - **aromatic** — aromatic ring centroid
        - **hydrophobic** — hydrophobic region
        - **positive** — positively ionisable group
        - **negative** — negatively ionisable group

    When RDKit's ``BaseFeatures.fdef`` cannot be loaded (e.g. minimal
    installs), a heuristic atom-based fallback is used.

    Args:
        smiles: SMILES string.
        conformer_id: Index of the conformer to use (default 0).

    Returns:
        List of dicts with keys matching :class:`PharmacophoreFeature3D`
        fields (``feature_type``, ``position``, ``direction``,
        ``atom_indices``).  Empty list when RDKit is unavailable.
    """
    if not HAS_RDKIT:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    mol = Chem.AddHs(mol)

    # Generate a 3D conformer if one does not already exist
    if mol.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol, params)
        if status == -1:
            logger.warning("Could not embed molecule for pharmacophore extraction: %s", smiles)
            return []
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                if ff is not None:
                    ff.Minimize(maxIters=500)
        except Exception:
            pass

    # Ensure requested conformer_id exists
    if conformer_id >= mol.GetNumConformers():
        conformer_id = 0

    features: list[dict] = []

    # Attempt feature-based pharmacophore extraction via BaseFeatures.fdef
    try:
        import os
        fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
        if os.path.isfile(fdef_path):
            feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
            mol_feats = feat_factory.GetFeaturesForMol(mol, confId=conformer_id)

            for mf in mol_feats:
                feat_type = mf.GetFamily()
                # Map RDKit families to our canonical names
                type_map = {
                    "Donor": "HBD",
                    "Acceptor": "HBA",
                    "Aromatic": "aromatic",
                    "Hydrophobe": "hydrophobic",
                    "PosIonizable": "positive",
                    "NegIonizable": "negative",
                    "LumpedHydrophobe": "hydrophobic",
                }
                canonical = type_map.get(feat_type)
                if canonical is None:
                    continue

                atom_indices = list(mf.GetAtomIds())
                pos = mf.GetPos()
                position = (pos.x, pos.y, pos.z)

                # Direction (available for some features)
                direction = (0.0, 0.0, 0.0)
                try:
                    dir_vec = mf.GetDirection()
                    if dir_vec is not None:
                        direction = (dir_vec.x, dir_vec.y, dir_vec.z)
                except Exception:
                    pass

                features.append({
                    "feature_type": canonical,
                    "position": position,
                    "direction": direction,
                    "atom_indices": atom_indices,
                })

            if features:
                return features
    except Exception:
        logger.debug("Feature-based pharmacophore extraction failed; using heuristic fallback")

    # Heuristic atom-based fallback
    features = _extract_pharmacophore_heuristic(mol, conformer_id)
    return features


def _extract_pharmacophore_heuristic(mol, conformer_id: int = 0) -> list[dict]:
    """Heuristic pharmacophore extraction when BaseFeatures.fdef is unavailable.

    Uses simple SMARTS patterns and atom properties to identify
    pharmacophore features and computes centroids from conformer
    coordinates.

    Args:
        mol: RDKit Mol with at least one conformer and explicit Hs.
        conformer_id: Conformer index.

    Returns:
        List of pharmacophore feature dicts.
    """
    conf = mol.GetConformer(conformer_id)
    features: list[dict] = []

    def _centroid(atom_ids: list[int]) -> tuple[float, float, float]:
        """Compute centroid of atoms."""
        cx, cy, cz = 0.0, 0.0, 0.0
        for aid in atom_ids:
            pos = conf.GetAtomPosition(aid)
            cx += pos.x
            cy += pos.y
            cz += pos.z
        n = len(atom_ids) or 1
        return (cx / n, cy / n, cz / n)

    # HBD — N-H or O-H donors
    hbd_pattern = Chem.MolFromSmarts("[N,O;!H0]")
    if hbd_pattern is not None:
        matches = mol.GetSubstructMatches(hbd_pattern)
        seen: set[int] = set()
        for match in matches:
            aid = match[0]
            if aid not in seen:
                seen.add(aid)
                # Direction: from donor heavy atom toward its H
                atom = mol.GetAtomWithIdx(aid)
                direction = (0.0, 0.0, 0.0)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        h_pos = conf.GetAtomPosition(neighbor.GetIdx())
                        d_pos = conf.GetAtomPosition(aid)
                        dx = h_pos.x - d_pos.x
                        dy = h_pos.y - d_pos.y
                        dz = h_pos.z - d_pos.z
                        norm = math.sqrt(dx * dx + dy * dy + dz * dz) or 1.0
                        direction = (dx / norm, dy / norm, dz / norm)
                        break
                features.append({
                    "feature_type": "HBD",
                    "position": _centroid([aid]),
                    "direction": direction,
                    "atom_indices": [aid],
                })

    # HBA — N or O with lone pairs
    hba_pattern = Chem.MolFromSmarts("[N,O;!+1]")
    if hba_pattern is not None:
        matches = mol.GetSubstructMatches(hba_pattern)
        seen_hba: set[int] = set()
        for match in matches:
            aid = match[0]
            if aid not in seen_hba:
                seen_hba.add(aid)
                features.append({
                    "feature_type": "HBA",
                    "position": _centroid([aid]),
                    "direction": (0.0, 0.0, 0.0),
                    "atom_indices": [aid],
                })

    # Aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_ring_atoms: list[list[int]] = []
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(aid).GetIsAromatic() for aid in ring):
            aromatic_ring_atoms.append(list(ring))

    for ring_atoms in aromatic_ring_atoms:
        features.append({
            "feature_type": "aromatic",
            "position": _centroid(ring_atoms),
            "direction": (0.0, 0.0, 0.0),
            "atom_indices": ring_atoms,
        })

    # Hydrophobic — carbon atoms not in polar environments
    hydrophobic_pattern = Chem.MolFromSmarts("[C;!+1;!$(C=[N,O,S])]")
    if hydrophobic_pattern is not None:
        matches = mol.GetSubstructMatches(hydrophobic_pattern)
        if matches:
            hydro_atoms = sorted(set(m[0] for m in matches))
            features.append({
                "feature_type": "hydrophobic",
                "position": _centroid(hydro_atoms),
                "direction": (0.0, 0.0, 0.0),
                "atom_indices": hydro_atoms,
            })

    # Positive ionisable — amines, amidines, guanidines
    pos_patterns = [
        Chem.MolFromSmarts("[N;+;!H0]"),
        Chem.MolFromSmarts("[n;+]"),
        Chem.MolFromSmarts("[NH2,NH1,NH0]C(=[NH2+,NH1+,NH0+])NH2"),
    ]
    pos_atoms: set[int] = set()
    for pat in pos_patterns:
        if pat is not None:
            for match in mol.GetSubstructMatches(pat):
                pos_atoms.update(match)
    if pos_atoms:
        aids = sorted(pos_atoms)
        features.append({
            "feature_type": "positive",
            "position": _centroid(aids),
            "direction": (0.0, 0.0, 0.0),
            "atom_indices": aids,
        })

    # Negative ionisable — carboxylates, phosphates, sulfonates
    neg_patterns = [
        Chem.MolFromSmarts("[O-]C(=O)"),
        Chem.MolFromSmarts("[O-]P(=O)"),
        Chem.MolFromSmarts("[O-]S(=O)(=O)"),
    ]
    neg_atoms: set[int] = set()
    for pat in neg_patterns:
        if pat is not None:
            for match in mol.GetSubstructMatches(pat):
                neg_atoms.update(match)
    if neg_atoms:
        aids = sorted(neg_atoms)
        features.append({
            "feature_type": "negative",
            "position": _centroid(aids),
            "direction": (0.0, 0.0, 0.0),
            "atom_indices": aids,
        })

    return features


# ---------------------------------------------------------------------------
# RMSD between conformers
# ---------------------------------------------------------------------------

def compute_rmsd(
    smiles: str,
    conformer_i: int,
    conformer_j: int,
) -> float | None:
    """Compute the RMSD between two conformers of the same molecule.

    The conformers are generated on-the-fly using ETKDGv3 + MMFF94
    minimisation.  If the requested conformer indices do not exist,
    ``None`` is returned.

    Args:
        smiles: SMILES string.
        conformer_i: Index of the first conformer.
        conformer_j: Index of the second conformer.

    Returns:
        RMSD in Ångströms, or ``None`` if RDKit is unavailable,
        the SMILES is invalid, or the conformer indices are out of range.
    """
    if not HAS_RDKIT:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Generate enough conformers to cover both indices
    n_needed = max(conformer_i, conformer_j) + 1
    params = AllChem.ETKDGv3()
    params.numConformers = n_needed
    params.randomSeed = 42

    cids = AllChem.EmbedMultipleConfs(mol, params)
    if not cids or len(cids) < 2:
        return None

    # Minimise
    try:
        for cid in cids:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                if ff is not None:
                    ff.Minimize(maxIters=500)
    except Exception:
        pass

    # Map requested indices to actual conformer IDs
    if conformer_i >= len(cids) or conformer_j >= len(cids):
        return None

    cid_i = cids[conformer_i]
    cid_j = cids[conformer_j]

    try:
        rmsd = AllChem.GetConformerRMS(mol, cid_i, cid_j)
        return rmsd
    except Exception:
        logger.debug("RMSD computation failed for conformers %d, %d", cid_i, cid_j)
        return None


# ---------------------------------------------------------------------------
# Conformer clustering by RMSD (DBSCAN)
# ---------------------------------------------------------------------------

def cluster_poses_by_rmsd(
    smiles: str,
    n_conformers: int = 20,
    rmsd_threshold: float = 2.0,
) -> list[list[int]]:
    """Cluster conformers by pairwise RMSD using DBSCAN.

    Generates *n_conformers* 3D conformers, computes the full pairwise
    RMSD matrix, and applies DBSCAN clustering with the given distance
    threshold (``eps``).

    DBSCAN is chosen because it does not require specifying the number
    of clusters in advance and naturally handles noise (outlier conformers)
    by assigning them to cluster ``-1``.

    Args:
        smiles: SMILES string.
        n_conformers: Number of conformers to generate (default 20).
        rmsd_threshold: DBSCAN ``eps`` parameter — maximum RMSD (Å)
            for two conformers to be considered neighbours (default 2.0).

    Returns:
        List of clusters, each cluster being a list of conformer indices.
        Conformers classified as noise are placed in a cluster labelled
        ``[-1]``.  Returns an empty list when RDKit is unavailable or
        the SMILES is invalid.
    """
    if not HAS_RDKIT:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    mol = Chem.AddHs(mol)

    # Generate conformers
    params = AllChem.ETKDGv3()
    params.numConformers = n_conformers
    params.randomSeed = 42

    cids = AllChem.EmbedMultipleConfs(mol, params)
    if not cids or len(cids) < 2:
        return [[0]] if cids else []

    # Minimise each conformer
    for cid in cids:
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                if ff is not None:
                    ff.Minimize(maxIters=500)
        except Exception:
            pass

    # Compute pairwise RMSD matrix
    n = len(cids)
    dist_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            try:
                rmsd = AllChem.GetConformerRMS(mol, cids[i], cids[j])
            except Exception:
                rmsd = float("inf")
            dist_matrix[i][j] = rmsd
            dist_matrix[j][i] = rmsd

    # DBSCAN implementation (no sklearn dependency required)
    clusters = _dbscan(dist_matrix, eps=rmsd_threshold, min_samples=1)

    # Organise results: list of clusters, each cluster is a list of indices
    cluster_map: dict[int, list[int]] = {}
    for idx, label in enumerate(clusters):
        if label not in cluster_map:
            cluster_map[label] = []
        cluster_map[label].append(idx)

    # Sort clusters by size (largest first); noise cluster (-1) last
    sorted_clusters = sorted(
        cluster_map.items(),
        key=lambda kv: (kv[0] == -1, -len(kv[1])),
    )

    return [indices for _, indices in sorted_clusters]


def _dbscan(
    dist_matrix: list[list[float]],
    eps: float,
    min_samples: int = 1,
) -> list[int]:
    """Simple DBSCAN clustering on a precomputed distance matrix.

    Args:
        dist_matrix: NxN symmetric distance matrix.
        eps: Neighbourhood radius.
        min_samples: Minimum number of points to form a dense region.

    Returns:
        List of cluster labels (int).  Noise points are labelled ``-1``.
    """
    n = len(dist_matrix)
    labels = [-1] * n  # -1 = unvisited / noise
    visited = [False] * n
    cluster_id = 0

    def _region_query(point: int) -> list[int]:
        """Return indices of points within eps of *point*."""
        neighbours = []
        for other in range(n):
            if dist_matrix[point][other] <= eps:
                neighbours.append(other)
        return neighbours

    for i in range(n):
        if visited[i]:
            continue
        visited[i] = True
        neighbours = _region_query(i)

        if len(neighbours) < min_samples:
            # Noise (may be later claimed by another cluster)
            labels[i] = -1
        else:
            # Start a new cluster
            labels[i] = cluster_id
            seed_set = list(neighbours)
            seed_set = [s for s in seed_set if s != i]

            j = 0
            while j < len(seed_set):
                q = seed_set[j]
                if not visited[q]:
                    visited[q] = True
                    q_neighbours = _region_query(q)
                    if len(q_neighbours) >= min_samples:
                        for nb in q_neighbours:
                            if nb not in seed_set:
                                seed_set.append(nb)
                if labels[q] == -1:
                    labels[q] = cluster_id
                j += 1

            cluster_id += 1

    return labels
