"""
BioCompiler AutoDock Vina Wrapper
==================================
Dedicated wrapper module for AutoDock Vina molecular docking with
graceful dependency handling and BioCompiler integration.

This module provides:
- VinaResult / DockingConfig dataclasses for structured docking I/O
- SMILES→PDBQT conversion (meeko/RDKit)
- PDB→PDBQT conversion (OpenBabel / MGLTools prepare_receptor)
- Auto-detection of search box from receptor geometry
- Full docking workflow via dock_smiles() and dock_pdbqt()
- Integrated binding assessment via score_ligand_binding()

All functions degrade gracefully when optional dependencies
(vina, meeko, rdkit, openbabel) are unavailable, returning
error-bearing VinaResult objects or informative dicts rather
than raising ImportError.

References
----------
- Trott O & Olson AJ (2010) AutoDock Vina: improving the speed and
  accuracy of docking with a new scoring function, efficient
  optimization, and multithreading. J Comput Chem 31:455-461.
  doi:10.1002/jcc.21334
- Eberhardt J, Santos-Martins D, Tillack AF & Forli S (2021)
  AutoDock Vina 1.2.0: New Docking Methods, Enhanced Scoring,
  and GPU Acceleration. J Chem Inf Model 61(8):3891-3898.
  doi:10.1021/acs.jcim.1c01178
- Corso G, Stark H, Jing B, Barzilay R, Jaakkola T (2024)
  DiffDock: Diffusion Steps, Twists, and Turns for Molecular
  Docking. ICLR 2024. (Referenced for future integration.)
"""

from __future__ import annotations

import logging
import math
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Any, Optional

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# Optional dependency detection
# ────────────────────────────────────────────────────────────

try:
    from vina import Vina as VinaDock
    HAS_VINA = True
except ImportError:
    HAS_VINA = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from meeko import MoleculePreparation
    HAS_MEEKO = True
except ImportError:
    HAS_MEEKO = False

try:
    from openbabel import openbabel as ob
    HAS_OPENBABEL = True
except ImportError:
    HAS_OPENBABEL = False


__all__ = [
    # Data classes
    "VinaResult",
    "DockingConfig",
    # Core docking
    "dock_smiles",
    "dock_pdbqt",
    # Format conversion
    "smiles_to_pdbqt",
    "pdb_to_pdbqt",
    # Box auto-detection
    "compute_box_from_receptor",
    # Availability check
    "is_vina_available",
    # BioCompiler integration
    "score_ligand_binding",
    # Feature flags
    "HAS_VINA",
    "HAS_RDKIT",
    "HAS_MEEKO",
    "HAS_OPENBABEL",
]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class VinaResult:
    """Structured result from an AutoDock Vina docking calculation.

    Attributes:
        best_energy: Best (lowest) binding affinity in kcal/mol.
            None if docking failed or Vina unavailable.
        all_energies: List of binding affinities for all generated
            poses, ordered best-first (lowest first).
        n_poses: Number of poses actually generated.
        method: Docking method used (e.g. "vina", "vinardo",
            "ad4", or "unavailable").
        error: Error message if docking failed; None on success.
    """
    best_energy: float | None = None
    all_energies: list[float] = field(default_factory=list)
    n_poses: int = 0
    method: str = "unavailable"
    error: str | None = None


@dataclass
class DockingConfig:
    """Configuration for an AutoDock Vina docking run.

    Attributes:
        exhaustiveness: Search thoroughness. Higher values yield
            more reliable results at the cost of runtime.
            Default 32 (Trott & Olson 2010 recommend ≥8).
        n_poses: Maximum number of poses to generate.
            Default 20.
        energy_range: Maximum energy difference (kcal/mol) between
            the best and worst pose reported. Default 3.
        box_center: (x, y, z) coordinates of the search box
            center in Angstroms. Required for docking.
        box_size: (sx, sy, sz) dimensions of the search box
            in Angstroms. Must be ≥10.75 per dimension.
    """
    exhaustiveness: int = 32
    n_poses: int = 20
    energy_range: float = 3.0
    box_center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0)


# ────────────────────────────────────────────────────────────
# Internal PDBQT parsing helper
# ────────────────────────────────────────────────────────────

def _parse_pdbqt_atoms(pdbqt_path: str) -> list[dict[str, Any]]:
    """Parse atom records from a PDBQT file for box computation.

    Returns list of dicts with keys: x, y, z, atom_type.
    Only ATOM and HETATM records are parsed.
    """
    atoms: list[dict[str, Any]] = []
    try:
        with open(pdbqt_path, "r") as fh:
            for line in fh:
                rec = line[:6].strip() if len(line) >= 6 else ""
                if rec in ("ATOM", "HETATM") and len(line) >= 54:
                    try:
                        atoms.append({
                            "x": float(line[30:38]),
                            "y": float(line[38:46]),
                            "z": float(line[46:54]),
                            "atom_type": line[77:80].strip() if len(line) >= 80 else "",
                        })
                    except (ValueError, IndexError):
                        continue
    except OSError as exc:
        logger.warning("Cannot read PDBQT file %s: %s", pdbqt_path, exc)
    return atoms


# ────────────────────────────────────────────────────────────
# Availability check
# ────────────────────────────────────────────────────────────

def is_vina_available() -> bool:
    """Check whether the AutoDock Vina Python bindings are installed.

    Returns True if ``from vina import Vina`` succeeds, False otherwise.

    Examples
    --------
    >>> is_vina_available()
    False  # unless vina package is installed
    """
    return HAS_VINA


# ────────────────────────────────────────────────────────────
# Format conversion: SMILES → PDBQT
# ────────────────────────────────────────────────────────────

def smiles_to_pdbqt(
    smiles: str,
    output_path: str | None = None,
) -> str:
    """Convert a SMILES string to PDBQT format for Vina docking.

    Uses meeko for PDBQT preparation when available (recommended),
    with a fallback to RDKit + OpenBabel.

    Parameters
    ----------
    smiles : str
        Input SMILES string (e.g. ``"CC(=O)OC1=CC=CC=C1C(=O)O"``
        for aspirin).
    output_path : str or None
        Path to write the PDBQT file. If None, a temporary file
        is created and its path is returned.

    Returns
    -------
    str
        Path to the generated PDBQT file.

    Raises
    ------
    RuntimeError
        If neither meeko nor (RDKit + OpenBabel) are available,
        or if the SMILES string cannot be parsed.
    """
    if not smiles:
        raise RuntimeError("Empty SMILES string provided")

    pdbqt_string: str | None = None

    # Strategy 1: meeko (preferred — produces Vina-compatible PDBQT)
    if HAS_MEEKO and HAS_RDKIT:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise RuntimeError(f"RDKit cannot parse SMILES: {smiles}")
            mol = Chem.AddHs(mol)
            preparator = MoleculePreparation()
            preparator.prepare(mol)
            pdbqt_string = preparator.write_pdbqt_string()
            logger.debug("Converted SMILES→PDBQT via meeko for %s", smiles[:40])
        except Exception as exc:
            logger.warning("meeko conversion failed: %s; trying fallback", exc)

    # Strategy 2: RDKit 3D embedding + OpenBabel conversion
    if pdbqt_string is None and HAS_RDKIT and HAS_OPENBABEL:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise RuntimeError(f"RDKit cannot parse SMILES: {smiles}")
            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                result = AllChem.EmbedMolecule(
                    mol, randomSeed=42, useRandomCoords=True
                )
            if result == -1:
                raise RuntimeError(
                    f"Cannot generate 3D coordinates for SMILES: {smiles}"
                )
            # Minimise with MMFF94
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass  # non-fatal

            # Write to SDF, convert with OpenBabel
            with tempfile.NamedTemporaryFile(
                suffix=".sdf", delete=False, mode="w"
            ) as sdf_tmp:
                writer = Chem.SDWriter(sdf_tmp.name)
                writer.write(mol)
                writer.close()
                sdf_path = sdf_tmp.name

            try:
                ob_conversion = ob.OBConversion()
                ob_conversion.SetInAndOutFormats("sdf", "pdbqt")
                ob_mol = ob.OBMol()
                ob_conversion.ReadFile(ob_mol, sdf_path)
                pdbqt_string = ob_conversion.WriteString(ob_mol)
                logger.debug(
                    "Converted SMILES→PDBQT via RDKit+OpenBabel for %s",
                    smiles[:40],
                )
            finally:
                try:
                    os.unlink(sdf_path)
                except OSError:
                    pass
        except Exception as exc:
            logger.warning("RDKit+OpenBabel conversion failed: %s", exc)

    # Strategy 3: OpenBabel CLI fallback
    if pdbqt_string is None and HAS_OPENBABEL:
        try:
            with tempfile.NamedTemporaryFile(
                suffix=".smi", delete=False, mode="w"
            ) as smi_tmp:
                smi_tmp.write(smiles)
                smi_path = smi_tmp.name
            pdbqt_tmp_path = smi_path.replace(".smi", ".pdbqt")
            try:
                subprocess.run(
                    ["obabel", smi_path, "-O", pdbqt_tmp_path, "--gen3d", "-h"],
                    capture_output=True,
                    timeout=60,
                    check=True,
                )
                with open(pdbqt_tmp_path, "r") as fh:
                    pdbqt_string = fh.read()
                logger.debug(
                    "Converted SMILES→PDBQT via OpenBabel CLI for %s",
                    smiles[:40],
                )
            finally:
                for p in (smi_path, pdbqt_tmp_path):
                    try:
                        os.unlink(p)
                    except OSError:
                        pass
        except Exception as exc:
            logger.warning("OpenBabel CLI conversion failed: %s", exc)

    if pdbqt_string is None:
        raise RuntimeError(
            "Cannot convert SMILES to PDBQT: install meeko and rdkit "
            "(pip install meeko rdkit), or openbabel "
            "(conda install -c conda-forge openbabel)"
        )

    # Write output
    if output_path is None:
        fd, output_path = tempfile.mkstemp(suffix="_ligand.pdbqt")
        with os.fdopen(fd, "w") as fh:
            fh.write(pdbqt_string)
    else:
        with open(output_path, "w") as fh:
            fh.write(pdbqt_string)

    return output_path


# ────────────────────────────────────────────────────────────
# Format conversion: PDB → PDBQT
# ────────────────────────────────────────────────────────────

def pdb_to_pdbqt(
    pdb_path: str,
    output_path: str | None = None,
) -> str:
    """Convert a PDB file to PDBQT format for Vina receptor input.

    Tries OpenBabel first (preserves protonation), then falls back
    to the MGLTools ``prepare_receptor4.py`` script.

    Parameters
    ----------
    pdb_path : str
        Path to the input PDB file.
    output_path : str or None
        Path to write the PDBQT file. If None, a temporary file
        is created.

    Returns
    -------
    str
        Path to the generated PDBQT file.

    Raises
    ------
    FileNotFoundError
        If the input PDB file does not exist.
    RuntimeError
        If neither OpenBabel nor MGLTools are available.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if output_path is None:
        fd, output_path = tempfile.mkstemp(suffix="_receptor.pdbqt")
        os.close(fd)

    # Strategy 1: OpenBabel Python bindings
    if HAS_OPENBABEL:
        try:
            ob_conversion = ob.OBConversion()
            ob_conversion.SetInAndOutFormats("pdb", "pdbqt")
            ob_mol = ob.OBMol()
            if ob_conversion.ReadFile(ob_mol, pdb_path):
                # Add polar hydrogens if missing
                ob_mol.AddPolarHydrogens()
                ob_conversion.WriteFile(ob_mol, output_path)
                logger.debug(
                    "Converted PDB→PDBQT via OpenBabel for %s", pdb_path
                )
                return output_path
        except Exception as exc:
            logger.warning("OpenBabel PDB→PDBQT conversion failed: %s", exc)

    # Strategy 2: OpenBabel CLI
    obabel_bin = shutil.which("obabel")
    if obabel_bin is not None:
        try:
            subprocess.run(
                [obabel_bin, pdb_path, "-O", output_path, "-h", "-xr"],
                capture_output=True,
                timeout=60,
                check=True,
            )
            logger.debug("Converted PDB→PDBQT via obabel CLI for %s", pdb_path)
            return output_path
        except Exception as exc:
            logger.warning("obabel CLI conversion failed: %s", exc)

    # Strategy 3: MGLTools prepare_receptor4.py
    prepare_script = shutil.which("prepare_receptor4.py")
    if prepare_script is None:
        # Check common MGLTools install locations
        mgl_candidates = [
            "/usr/local/MGLTools-1.5.7/bin/prepare_receptor4.py",
            "/opt/mgltools/bin/prepare_receptor4.py",
        ]
        for candidate in mgl_candidates:
            if os.path.isfile(candidate):
                prepare_script = candidate
                break

    if prepare_script is not None:
        try:
            subprocess.run(
                [
                    "python",
                    prepare_script,
                    "-r", pdb_path,
                    "-o", output_path,
                ],
                capture_output=True,
                timeout=60,
                check=True,
            )
            logger.debug(
                "Converted PDB→PDBQT via MGLTools for %s", pdb_path
            )
            return output_path
        except Exception as exc:
            logger.warning("MGLTools conversion failed: %s", exc)

    raise RuntimeError(
        "Cannot convert PDB to PDBQT: install OpenBabel "
        "(conda install -c conda-forge openbabel) or MGLTools"
    )


# ────────────────────────────────────────────────────────────
# Auto-detect search box from receptor
# ────────────────────────────────────────────────────────────

def compute_box_from_receptor(
    receptor_pdbqt: str,
    padding: float = 10.0,
) -> DockingConfig:
    """Auto-detect a Vina search box from the receptor geometry.

    Parses the receptor PDBQT file, computes the axis-aligned
    bounding box of all atoms, and expands it by ``padding``
    Angstroms on each side.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    padding : float
        Padding in Angstroms to add around the bounding box
        on each side. Default 10.0.

    Returns
    -------
    DockingConfig
        Configuration with box_center and box_size set from
        the receptor geometry. If the PDBQT file cannot be
        read, returns a default DockingConfig.
    """
    atoms = _parse_pdbqt_atoms(receptor_pdbqt)

    if not atoms:
        logger.warning(
            "No atoms found in %s; returning default DockingConfig",
            receptor_pdbqt,
        )
        return DockingConfig()

    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    zs = [a["z"] for a in atoms]

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    center_z = (min_z + max_z) / 2.0

    size_x = (max_x - min_x) + 2.0 * padding
    size_y = (max_y - min_y) + 2.0 * padding
    size_z = (max_z - min_z) + 2.0 * padding

    return DockingConfig(
        box_center=(center_x, center_y, center_z),
        box_size=(size_x, size_y, size_z),
    )


# ────────────────────────────────────────────────────────────
# Core docking: PDBQT ligand
# ────────────────────────────────────────────────────────────

def dock_pdbqt(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    config: DockingConfig | None = None,
) -> VinaResult:
    """Dock a PDBQT-format ligand against a PDBQT receptor.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    ligand_pdbqt : str
        Path to the ligand PDBQT file.
    config : DockingConfig or None
        Docking configuration. If None, auto-detects the search
        box from the receptor and uses default parameters.

    Returns
    -------
    VinaResult
        Structured docking result. On failure, ``error`` is set
        and ``best_energy`` is None.
    """
    if not HAS_VINA:
        return VinaResult(
            error="AutoDock Vina Python bindings not installed. "
            "Install with: pip install vina",
            method="unavailable",
        )

    if config is None:
        config = compute_box_from_receptor(receptor_pdbqt)

    try:
        v = VinaDock(sf_name="vina")
        v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)

        v.compute_vina_maps(
            center=config.box_center,
            box_size=config.box_size,
        )
        v.dock(
            exhaustiveness=config.exhaustiveness,
            n_poses=config.n_poses,
            min_rmsd=1.0,
            max_evals=0,
        )

        energies = v.energies()
        # energies() returns list of [affinity, rmsd_lb, rmsd_ub]
        all_e = [row[0] for row in energies] if energies else []

        return VinaResult(
            best_energy=all_e[0] if all_e else None,
            all_energies=all_e,
            n_poses=len(all_e),
            method="vina",
        )
    except Exception as exc:
        return VinaResult(
            error=f"Vina docking failed: {exc}",
            method="vina",
        )


# ────────────────────────────────────────────────────────────
# Core docking: SMILES ligand
# ────────────────────────────────────────────────────────────

def dock_smiles(
    receptor_pdbqt: str,
    ligand_smiles: str,
    config: DockingConfig | None = None,
) -> VinaResult:
    """Dock a SMILES ligand against a receptor PDBQT file.

    Handles SMILES→PDBQT conversion automatically (using meeko/RDKit
    when available), then delegates to :func:`dock_pdbqt`.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    ligand_smiles : str
        SMILES string of the ligand.
    config : DockingConfig or None
        Docking configuration. If None, auto-detects the search
        box from the receptor and uses default parameters.

    Returns
    -------
    VinaResult
        Structured docking result. On failure (including SMILES
        conversion failure), ``error`` is set.
    """
    # Convert SMILES to PDBQT
    ligand_pdbqt_path: str | None = None
    try:
        ligand_pdbqt_path = smiles_to_pdbqt(ligand_smiles)
    except RuntimeError as exc:
        return VinaResult(
            error=f"SMILES→PDBQT conversion failed: {exc}",
            method="vina",
        )

    if ligand_pdbqt_path is None:
        return VinaResult(
            error="SMILES→PDBQT conversion returned no output path",
            method="vina",
        )

    # Run docking
    try:
        result = dock_pdbqt(receptor_pdbqt, ligand_pdbqt_path, config)
        return result
    finally:
        # Clean up temporary ligand PDBQT file
        if ligand_pdbqt_path and os.path.isfile(ligand_pdbqt_path):
            try:
                os.unlink(ligand_pdbqt_path)
            except OSError:
                pass


# ────────────────────────────────────────────────────────────
# BioCompiler integration: full binding assessment
# ────────────────────────────────────────────────────────────

def score_ligand_binding(
    protein_seq: str,
    ligand_smiles: str,
    receptor_pdb_path: str | None = None,
) -> dict[str, Any]:
    """Full binding assessment combining docking with pharmacophore scoring.

    This is the primary BioCompiler integration point for ligand
    binding evaluation. It combines:

    1. AutoDock Vina docking (if receptor structure available)
    2. Pharmacophore-based binding site scoring
    3. Ligand feature analysis (via RDKit or regex fallback)

    Parameters
    ----------
    protein_seq : str
        Amino acid sequence of the protein (single-letter codes).
        Used for pharmacophore scoring when PDB structure is
        unavailable.
    ligand_smiles : str
        SMILES string of the ligand.
    receptor_pdb_path : str or None
        Path to the receptor PDB file. If provided, the receptor
        will be converted to PDBQT and used for docking. If None,
        only sequence-based pharmacophore scoring is performed.

    Returns
    -------
    dict
        Binding assessment dictionary with keys:

        - ``"docking"``: VinaResult dict (or None if undockable)
        - ``"pharmacophore_score"``: float from binding site scoring
        - ``"ligand_features"``: dict of ligand molecular features
        - ``"binding_verdict"``: str — "strong", "moderate", "weak",
          or "none"
        - ``"error"``: str or None
    """
    result: dict[str, Any] = {
        "docking": None,
        "pharmacophore_score": 0.0,
        "ligand_features": {},
        "binding_verdict": "none",
        "error": None,
    }

    # ── Ligand feature extraction ──────────────────────────
    try:
        from .ligand_binding_v2 import parse_smiles_features_rdkit
        result["ligand_features"] = parse_smiles_features_rdkit(ligand_smiles)
    except (ImportError, Exception) as exc:
        logger.debug("Ligand feature extraction failed: %s", exc)
        # Minimal fallback
        result["ligand_features"] = {"smiles_length": len(ligand_smiles)}

    # ── Pharmacophore scoring ──────────────────────────────
    pharma_score = 0.0
    if receptor_pdb_path is not None and os.path.isfile(receptor_pdb_path):
        try:
            from .ligand_binding_v2 import (
                detect_binding_sites,
                score_binding_site,
            )
            with open(receptor_pdb_path, "r") as fh:
                pdb_string = fh.read()
            sites = detect_binding_sites(
                pdb_string,
                ligand_smiles=ligand_smiles,
            )
            if sites:
                best_site = max(sites, key=lambda s: s.score)
                pharma_score = best_site.score
        except (ImportError, Exception) as exc:
            logger.debug("Pharmacophore scoring failed: %s", exc)
            # Sequence-based heuristic fallback
            pharma_score = _pharmacophore_score_from_sequence(
                protein_seq, ligand_smiles
            )
    else:
        # No receptor PDB — use sequence-based heuristic
        pharma_score = _pharmacophore_score_from_sequence(
            protein_seq, ligand_smiles
        )

    result["pharmacophore_score"] = pharma_score

    # ── Docking ────────────────────────────────────────────
    vina_result: VinaResult | None = None
    if receptor_pdb_path is not None and os.path.isfile(receptor_pdb_path):
        try:
            # Convert receptor PDB → PDBQT
            receptor_pdbqt_path = pdb_to_pdbqt(receptor_pdb_path)

            # Compute search box from receptor
            config = compute_box_from_receptor(receptor_pdbqt_path)

            # Dock
            vina_result = dock_smiles(
                receptor_pdbqt_path, ligand_smiles, config
            )
            result["docking"] = {
                "best_energy": vina_result.best_energy,
                "all_energies": vina_result.all_energies,
                "n_poses": vina_result.n_poses,
                "method": vina_result.method,
                "error": vina_result.error,
            }

            # Clean up temporary receptor PDBQT
            try:
                os.unlink(receptor_pdbqt_path)
            except OSError:
                pass
        except (RuntimeError, FileNotFoundError, Exception) as exc:
            result["error"] = f"Docking pipeline failed: {exc}"
            logger.warning("Docking pipeline failed: %s", exc)

    # ── Verdict ────────────────────────────────────────────
    result["binding_verdict"] = _compute_binding_verdict(
        vina_result, pharma_score
    )

    return result


# ────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────

def _pharmacophore_score_from_sequence(
    protein_seq: str,
    ligand_smiles: str,
) -> float:
    """Estimate pharmacophore compatibility from sequence alone.

    This is a heuristic that scores complementarity between
    protein sequence composition and ligand features without
    3D structural information.

    Parameters
    ----------
    protein_seq : str
        Amino acid sequence (single-letter codes).
    ligand_smiles : str
        SMILES string of the ligand.

    Returns
    -------
    float
        Heuristic pharmacophore compatibility score (≥0).
    """
    score = 0.0
    seq = protein_seq.upper()
    n = len(seq) or 1

    # Ligand feature extraction
    lig_features: dict[str, Any] = {}
    try:
        from .ligand_binding_v2 import parse_smiles_features_rdkit
        lig_features = parse_smiles_features_rdkit(ligand_smiles)
    except (ImportError, Exception):
        pass

    # Charged residue complementarity
    pos_count = sum(1 for aa in seq if aa in "KRH")
    neg_count = sum(1 for aa in seq if aa in "DE")
    ligand_charged = lig_features.get("charged_groups", 0)
    if ligand_charged > 0 and (pos_count > 0 or neg_count > 0):
        charged_frac = (pos_count + neg_count) / n
        score += min(charged_frac * 10.0, 2.0)  # cap at 2.0

    # Hydrogen bond donors/acceptors
    ligand_donors = lig_features.get("hbond_donors", 0)
    ligand_acceptors = lig_features.get("hbond_acceptors", 0)
    donor_residues = sum(1 for aa in seq if aa in "STCNQ")
    acceptor_residues = sum(1 for aa in seq if aa in "DEKRH")
    if ligand_donors > 0 and acceptor_residues > 0:
        score += min(acceptor_residues / n * 5.0, 1.5)
    if ligand_acceptors > 0 and donor_residues > 0:
        score += min(donor_residues / n * 5.0, 1.5)

    # Aromatic complementarity
    aromatic_residues = sum(1 for aa in seq if aa in "FWYH")
    ligand_aromatic = lig_features.get("aromatic_rings", 0)
    if ligand_aromatic > 0 and aromatic_residues > 0:
        score += min(aromatic_residues / n * 8.0, 2.0)

    # Hydrophobic complementarity
    hydrophobic_residues = sum(1 for aa in seq if aa in "AVLIMFPW")
    if lig_features.get("logp", 0) > 2.0 and hydrophobic_residues / n > 0.3:
        score += 1.0

    return score


def _compute_binding_verdict(
    vina_result: VinaResult | None,
    pharma_score: float,
) -> str:
    """Compute an overall binding verdict from docking and pharmacophore.

    Parameters
    ----------
    vina_result : VinaResult or None
        Docking result (None if docking was not performed).
    pharma_score : float
        Pharmacophore compatibility score.

    Returns
    -------
    str
        One of "strong", "moderate", "weak", "none".
    """
    # Docking-based verdict (binding energy in kcal/mol)
    energy_score = 0.0
    if vina_result is not None and vina_result.best_energy is not None:
        be = vina_result.best_energy
        if be <= -9.0:
            energy_score = 4.0
        elif be <= -7.0:
            energy_score = 3.0
        elif be <= -5.0:
            energy_score = 2.0
        elif be <= -3.0:
            energy_score = 1.0
        else:
            energy_score = 0.0

    # Pharmacophore-based score contribution
    pharma_contribution = 0.0
    if pharma_score >= 5.0:
        pharma_contribution = 2.0
    elif pharma_score >= 3.0:
        pharma_contribution = 1.5
    elif pharma_score >= 1.0:
        pharma_contribution = 1.0
    elif pharma_score > 0.0:
        pharma_contribution = 0.5

    combined = energy_score + pharma_contribution

    if combined >= 5.0:
        return "strong"
    elif combined >= 3.0:
        return "moderate"
    elif combined >= 1.5:
        return "weak"
    else:
        return "none"
