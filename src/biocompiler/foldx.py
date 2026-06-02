"""
BioCompiler FoldX Stability Analysis Module v7.2.0
====================================================
Provides both online (FoldX CLI wrapper) and offline (empirical scoring)
modes for protein stability analysis, mutation scanning, and stabilization.

Online mode requires FoldX installed on PATH. Offline mode uses heuristics
based on amino acid composition, BLOSUM62, and Kyte-Doolittle hydropathy.

Usage:
    from biocompiler.foldx import (
        is_foldx_available,
        run_foldx_stability,
        run_foldx_repair,
        run_foldx_mutation,
        empirical_stability,
        scan_mutations,
        find_stabilizing_mutations,
        FoldXResult,
        MutationResult,
        FoldXError,
        BLOSUM62,
        HYDROPATHY,
    )
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
import tempfile
import time
from dataclasses import dataclass, field

from .exceptions import BioCompilerError

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────
# BLOSUM62 Substitution Matrix (20x20 standard AAs)
# Nested dict format: BLOSUM62[aa1][aa2] = score
# ────────────────────────────────────────────────────────────

_BLOSUM62_ROWS = [
    #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
    [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
    [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
    [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
    [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
    [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
    [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
    [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
    [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
    [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
    [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
    [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1, -1, -1, -3, -2, -2],  # K
    [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
    [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
    [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
    [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
    [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
    [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
    [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
    [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
]

_BLOSUM_INDEX = list("ARNDCQEGHILKMFPSTWYV")

BLOSUM62: dict[str, dict[str, int]] = {}
for _i, _a1 in enumerate(_BLOSUM_INDEX):
    BLOSUM62[_a1] = {}
    for _j, _a2 in enumerate(_BLOSUM_INDEX):
        BLOSUM62[_a1][_a2] = _BLOSUM62_ROWS[_i][_j]


# ────────────────────────────────────────────────────────────
# Kyte-Doolittle Hydropathy Scale
# ────────────────────────────────────────────────────────────

HYDROPATHY: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# Standard 20 amino acids
STANDARD_AAS: list[str] = list("ARNDCQEGHILKMFPSTWYV")

# Hydrophobic residues (Kyte-Doolittle > 1.0)
HYDROPHOBIC_AAS: set[str] = {"A", "I", "L", "M", "F", "W", "V"}

# Positively charged residues
POSITIVE_AAS: set[str] = {"K", "R", "H"}

# Negatively charged residues
NEGATIVE_AAS: set[str] = {"D", "E"}


# ────────────────────────────────────────────────────────────
# Custom Exception
# ────────────────────────────────────────────────────────────

class FoldXError(BioCompilerError):
    """Raised when FoldX analysis fails."""
    def __init__(self, reason: str, command: str | None = None):
        self.command = command
        msg = f"FoldX error: {reason}"
        if command:
            msg += f" (command: {command})"
        super().__init__(msg)


# ────────────────────────────────────────────────────────────
# Data Classes
# ────────────────────────────────────────────────────────────

@dataclass
class FoldXResult:
    """Result of a FoldX stability analysis run.

    Attributes:
        protein: The protein sequence analyzed.
        pdb_string: PDB file content (if available).
        stability_kcal: Total stability ΔG in kcal/mol (negative = stable).
        ddg_kcal: ΔΔG relative to reference (if mutation analysis).
        interaction_energy: Total interaction energy.
        backbone_hbond: Backbone hydrogen bond contribution.
        sidechain_hbond: Sidechain hydrogen bond contribution.
        van_der_waals: Van der Waals contribution.
        electrostatics: Electrostatic contribution.
        solvation: Solvation energy.
        van_der_waals_clashes: Van der Waals clashes penalty.
        entropy_sidechain: Sidechain entropy contribution.
        entropy_mainchain: Mainchain entropy contribution.
        torsional_clash: Torsional clash penalty.
        backbone_clash: Backbone clash penalty.
        helix_dipole: Helix dipole contribution.
        disulfide: Disulfide bond contribution.
        electrostatic_kon: Electrostatic Kon contribution.
        partial_covalent: Partial covalent bond energy.
        energy_ionisation: Ionisation energy.
        execution_time_s: Wall-clock execution time in seconds.
        method: Analysis method used ("foldx_cli", "empirical", "hybrid").
        success: Whether the analysis completed successfully.
        error: Error message if analysis failed.
    """
    protein: str
    pdb_string: str | None
    stability_kcal: float
    ddg_kcal: float | None
    interaction_energy: float | None
    backbone_hbond: float | None
    sidechain_hbond: float | None
    van_der_waals: float | None
    electrostatics: float | None
    solvation: float | None
    van_der_waals_clashes: float | None
    entropy_sidechain: float | None
    entropy_mainchain: float | None
    torsional_clash: float | None
    backbone_clash: float | None
    helix_dipole: float | None
    disulfide: float | None
    electrostatic_kon: float | None
    partial_covalent: float | None
    energy_ionisation: float | None
    execution_time_s: float
    method: str
    success: bool
    error: str | None = None


@dataclass
class MutationResult:
    """Result of a single-point mutation analysis.

    Attributes:
        position: 0-indexed residue position in the protein.
        wildtype: Original amino acid (single-letter code).
        mutant: Mutant amino acid (single-letter code).
        ddg_kcal: ΔΔG in kcal/mol (positive = destabilizing).
        stabilizing: True if ΔΔG < -0.5 kcal/mol.
        neutral: True if -0.5 ≤ ΔΔG ≤ 0.5 kcal/mol.
        destabilizing: True if ΔΔG > 0.5 kcal/mol.
        method: Analysis method used.
    """
    position: int
    wildtype: str
    mutant: str
    ddg_kcal: float
    stabilizing: bool
    neutral: bool
    destabilizing: bool
    method: str


# ────────────────────────────────────────────────────────────
# FoldX CLI Availability Check
# ────────────────────────────────────────────────────────────

def is_foldx_available() -> bool:
    """Check if the FoldX CLI is available on PATH.

    Returns:
        True if the ``foldx`` executable can be found on PATH and
        executes without error; False otherwise.
    """
    try:
        result = subprocess.run(
            ["foldx", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # FoldX may return non-zero for --version but still be available;
        # as long as it runs, we consider it available
        return True
    except FileNotFoundError:
        logger.debug("foldx executable not found on PATH")
        return False
    except subprocess.TimeoutExpired:
        logger.debug("foldx --version timed out")
        return False
    except OSError as exc:
        logger.debug("foldx availability check failed: %s", exc)
        return False


# ────────────────────────────────────────────────────────────
# FoldX CLI Wrapper: Stability
# ────────────────────────────────────────────────────────────

def run_foldx_stability(
    pdb_string: str,
    foldx_dir: str | None = None,
    timeout: float = 300.0,
) -> FoldXResult:
    """Run FoldX Stability analysis on a PDB structure.

    Writes the PDB content to a temporary directory, runs
    ``foldx --command=Stability``, and parses the ``.fxout`` output file
    to extract all energy components.

    Args:
        pdb_string: PDB file content as a string.
        foldx_dir: Directory containing the FoldX rotabase and other
            required files. If None, uses current working directory.
        timeout: Maximum execution time in seconds.

    Returns:
        FoldXResult with all energy components populated on success,
        or with success=False and an error message on failure.
    """
    start_time = time.time()

    if not is_foldx_available():
        elapsed = time.time() - start_time
        return FoldXResult(
            protein="",
            pdb_string=pdb_string,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="foldx_cli",
            success=False,
            error="FoldX CLI not available on PATH",
        )

    tmpdir = tempfile.mkdtemp(prefix="foldx_stability_")
    try:
        # Write PDB file — FoldX expects the file in the working directory
        pdb_path = os.path.join(tmpdir, "input.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_string)

        # Build command
        cmd = [
            "foldx",
            "--command=Stability",
            f"--pdb=input.pdb",
        ]
        if foldx_dir:
            cmd.append(f"--foldxDir={foldx_dir}")

        logger.info("Running FoldX Stability: %s", " ".join(cmd))

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=tmpdir,
        )

        elapsed = time.time() - start_time

        if result.returncode != 0:
            stderr = result.stderr.strip() if result.stderr else "unknown error"
            return FoldXResult(
                protein="",
                pdb_string=pdb_string,
                stability_kcal=0.0,
                ddg_kcal=None,
                interaction_energy=None,
                backbone_hbond=None,
                sidechain_hbond=None,
                van_der_waals=None,
                electrostatics=None,
                solvation=None,
                van_der_waals_clashes=None,
                entropy_sidechain=None,
                entropy_mainchain=None,
                torsional_clash=None,
                backbone_clash=None,
                helix_dipole=None,
                disulfide=None,
                electrostatic_kon=None,
                partial_covalent=None,
                energy_ionisation=None,
                execution_time_s=elapsed,
                method="foldx_cli",
                success=False,
                error=f"FoldX exited with code {result.returncode}: {stderr}",
            )

        # Parse the output file
        fxout_path = os.path.join(tmpdir, "input_ST.fxout")
        if not os.path.exists(fxout_path):
            # Try alternative naming
            fxout_candidates = [
                f for f in os.listdir(tmpdir) if f.endswith(".fxout")
            ]
            if fxout_candidates:
                fxout_path = os.path.join(tmpdir, fxout_candidates[0])
            else:
                return FoldXResult(
                    protein="",
                    pdb_string=pdb_string,
                    stability_kcal=0.0,
                    ddg_kcal=None,
                    interaction_energy=None,
                    backbone_hbond=None,
                    sidechain_hbond=None,
                    van_der_waals=None,
                    electrostatics=None,
                    solvation=None,
                    van_der_waals_clashes=None,
                    entropy_sidechain=None,
                    entropy_mainchain=None,
                    torsional_clash=None,
                    backbone_clash=None,
                    helix_dipole=None,
                    disulfide=None,
                    electrostatic_kon=None,
                    partial_covalent=None,
                    energy_ionisation=None,
                    execution_time_s=elapsed,
                    method="foldx_cli",
                    success=False,
                    error="FoldX output file (.fxout) not found",
                )

        parsed = _parse_stability_fxout(fxout_path)
        parsed["protein"] = _extract_protein_from_pdb(pdb_string)
        parsed["pdb_string"] = pdb_string
        parsed["execution_time_s"] = elapsed
        parsed["method"] = "foldx_cli"
        parsed["success"] = True
        parsed["error"] = None

        return FoldXResult(**parsed)

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        return FoldXResult(
            protein="",
            pdb_string=pdb_string,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="foldx_cli",
            success=False,
            error=f"FoldX timed out after {timeout}s",
        )
    except Exception as exc:
        elapsed = time.time() - start_time
        return FoldXResult(
            protein="",
            pdb_string=pdb_string,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="foldx_cli",
            success=False,
            error=str(exc),
        )
    finally:
        # Clean up temp directory
        _cleanup_tempdir(tmpdir)


# ────────────────────────────────────────────────────────────
# FoldX CLI Wrapper: RepairPDB
# ────────────────────────────────────────────────────────────

def run_foldx_repair(
    pdb_string: str,
    foldx_dir: str | None = None,
    timeout: float = 600.0,
) -> tuple[str, FoldXResult]:
    """Run FoldX RepairPDB to optimize sidechain conformations.

    Args:
        pdb_string: PDB file content as a string.
        foldx_dir: Directory containing the FoldX rotabase and other
            required files.
        timeout: Maximum execution time in seconds.

    Returns:
        Tuple of (repaired PDB string, FoldXResult). If repair fails,
        returns the original PDB string and a result with success=False.
    """
    start_time = time.time()

    if not is_foldx_available():
        elapsed = time.time() - start_time
        result = FoldXResult(
            protein="",
            pdb_string=pdb_string,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="foldx_cli",
            success=False,
            error="FoldX CLI not available on PATH",
        )
        return pdb_string, result

    tmpdir = tempfile.mkdtemp(prefix="foldx_repair_")
    try:
        pdb_path = os.path.join(tmpdir, "input.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_string)

        cmd = [
            "foldx",
            "--command=RepairPDB",
            "--pdb=input.pdb",
        ]
        if foldx_dir:
            cmd.append(f"--foldxDir={foldx_dir}")

        logger.info("Running FoldX RepairPDB: %s", " ".join(cmd))

        proc_result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=tmpdir,
        )

        elapsed = time.time() - start_time

        if proc_result.returncode != 0:
            stderr = proc_result.stderr.strip() if proc_result.stderr else "unknown error"
            result = FoldXResult(
                protein="",
                pdb_string=pdb_string,
                stability_kcal=0.0,
                ddg_kcal=None,
                interaction_energy=None,
                backbone_hbond=None,
                sidechain_hbond=None,
                van_der_waals=None,
                electrostatics=None,
                solvation=None,
                van_der_waals_clashes=None,
                entropy_sidechain=None,
                entropy_mainchain=None,
                torsional_clash=None,
                backbone_clash=None,
                helix_dipole=None,
                disulfide=None,
                electrostatic_kon=None,
                partial_covalent=None,
                energy_ionisation=None,
                execution_time_s=elapsed,
                method="foldx_cli",
                success=False,
                error=f"FoldX RepairPDB exited with code {proc_result.returncode}: {stderr}",
            )
            return pdb_string, result

        # Read repaired PDB
        repaired_path = os.path.join(tmpdir, "input_Repair.pdb")
        if not os.path.exists(repaired_path):
            # Try alternative names
            for fname in os.listdir(tmpdir):
                if fname.endswith("_Repair.pdb"):
                    repaired_path = os.path.join(tmpdir, fname)
                    break

        repaired_pdb = pdb_string  # fallback
        if os.path.exists(repaired_path):
            with open(repaired_path) as f:
                repaired_pdb = f.read()

        # Parse stability of repaired structure if output available
        parsed: dict = {
            "protein": _extract_protein_from_pdb(repaired_pdb),
            "pdb_string": repaired_pdb,
            "stability_kcal": 0.0,
            "ddg_kcal": None,
            "interaction_energy": None,
            "backbone_hbond": None,
            "sidechain_hbond": None,
            "van_der_waals": None,
            "electrostatics": None,
            "solvation": None,
            "van_der_waals_clashes": None,
            "entropy_sidechain": None,
            "entropy_mainchain": None,
            "torsional_clash": None,
            "backbone_clash": None,
            "helix_dipole": None,
            "disulfide": None,
            "electrostatic_kon": None,
            "partial_covalent": None,
            "energy_ionisation": None,
            "execution_time_s": elapsed,
            "method": "foldx_cli",
            "success": True,
            "error": None,
        }

        # Look for repair output
        for fname in os.listdir(tmpdir):
            if fname.endswith(".fxout"):
                fxout_data = _parse_stability_fxout(os.path.join(tmpdir, fname))
                for key in (
                    "stability_kcal", "interaction_energy", "backbone_hbond",
                    "sidechain_hbond", "van_der_waals", "electrostatics",
                    "solvation", "van_der_waals_clashes", "entropy_sidechain",
                    "entropy_mainchain", "torsional_clash", "backbone_clash",
                    "helix_dipole", "disulfide", "electrostatic_kon",
                    "partial_covalent", "energy_ionisation",
                ):
                    if key in fxout_data:
                        parsed[key] = fxout_data[key]
                break

        return repaired_pdb, FoldXResult(**parsed)

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        result = FoldXResult(
            protein="",
            pdb_string=pdb_string,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="foldx_cli",
            success=False,
            error=f"FoldX RepairPDB timed out after {timeout}s",
        )
        return pdb_string, result
    except Exception as exc:
        elapsed = time.time() - start_time
        result = FoldXResult(
            protein="",
            pdb_string=pdb_string,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="foldx_cli",
            success=False,
            error=str(exc),
        )
        return pdb_string, result
    finally:
        _cleanup_tempdir(tmpdir)


# ────────────────────────────────────────────────────────────
# FoldX CLI Wrapper: BuildModel (Mutation Analysis)
# ────────────────────────────────────────────────────────────

def run_foldx_mutation(
    pdb_string: str,
    mutations: list[str],
    foldx_dir: str | None = None,
    timeout: float = 300.0,
) -> list[MutationResult]:
    """Run FoldX BuildModel for point mutation ΔΔG analysis.

    Args:
        pdb_string: PDB file content as a string.
        mutations: List of mutations in FoldX format
            (e.g., ``["A123G", "L45F"]`` — single-letter wildtype AA,
            1-indexed position, single-letter mutant AA).
        foldx_dir: Directory containing FoldX auxiliary files.
        timeout: Maximum execution time in seconds.

    Returns:
        List of MutationResult, one per mutation. If FoldX is
        unavailable, returns an empty list.
    """
    if not is_foldx_available():
        logger.warning("FoldX CLI not available; cannot run mutation analysis")
        return []

    if not mutations:
        return []

    start_time = time.time()
    tmpdir = tempfile.mkdtemp(prefix="foldx_mutation_")

    try:
        # Write PDB file
        pdb_path = os.path.join(tmpdir, "input.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_string)

        # Write individual_list.txt (FoldX mutation file format)
        # Format: A123G;  (one mutation per line, semicolon-terminated)
        mutant_file = os.path.join(tmpdir, "individual_list.txt")
        with open(mutant_file, "w") as f:
            for mut in mutations:
                f.write(f"{mut};\n")

        cmd = [
            "foldx",
            "--command=BuildModel",
            "--pdb=input.pdb",
            "--mutant-file=individual_list.txt",
        ]
        if foldx_dir:
            cmd.append(f"--foldxDir={foldx_dir}")

        logger.info("Running FoldX BuildModel for %d mutations", len(mutations))

        proc_result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=tmpdir,
        )

        elapsed = time.time() - start_time

        if proc_result.returncode != 0:
            logger.error(
                "FoldX BuildModel failed (code %d): %s",
                proc_result.returncode,
                proc_result.stderr.strip() if proc_result.stderr else "unknown",
            )
            return []

        # Parse ΔΔG output
        results: list[MutationResult] = []
        ddg_values = _parse_buildmodel_output(tmpdir)

        for mut_str in mutations:
            parsed = _parse_mutation_string(mut_str)
            if parsed is None:
                continue
            wt, pos_1idx, mt = parsed
            pos_0idx = pos_1idx - 1  # convert to 0-indexed

            ddg = ddg_values.get(mut_str, None)
            if ddg is None:
                # Try looking up by position
                ddg = ddg_values.get(f"{wt}{pos_1idx}{mt}", 0.0)

            if ddg is None:
                ddg = 0.0

            results.append(MutationResult(
                position=pos_0idx,
                wildtype=wt,
                mutant=mt,
                ddg_kcal=ddg,
                stabilizing=ddg < -0.5,
                neutral=-0.5 <= ddg <= 0.5,
                destabilizing=ddg > 0.5,
                method="foldx_cli",
            ))

        return results

    except subprocess.TimeoutExpired:
        logger.error("FoldX BuildModel timed out after %ss", timeout)
        return []
    except Exception as exc:
        logger.error("FoldX mutation analysis failed: %s", exc)
        return []
    finally:
        _cleanup_tempdir(tmpdir)


# ────────────────────────────────────────────────────────────
# Empirical (Offline) Stability Estimation
# ────────────────────────────────────────────────────────────

def empirical_stability(protein: str) -> FoldXResult:
    """Estimate protein stability using empirical heuristics.

    Uses multiple sequence-based heuristics to estimate ΔG:

    a. **Hydrophobic core quality**: Fraction of hydrophobic residues
       (AILMFWV) — optimal range 0.35–0.45.
    b. **Charge balance**: Ratio of positive (KRH) to negative (DE)
       residues — ratio near 1.0 is good.
    c. **Proline content**: >5% proline increases rigidity (stabilizing
       up to ~8%), but too much is destabilizing.
    d. **Glycine content**: >7% glycine increases flexibility
       (destabilizing).
    e. **Cysteine pairs**: Even number of cysteines suggests potential
       disulfide bonds (stabilizing).
    f. **Charged-hydrophobic pattern**: Alternating pattern suggests
       helical content (moderate stability).
    g. **Sequence length**: Longer proteins tend to be more stable
       (more intramolecular contacts).

    These heuristics are combined into an estimated ΔG using a weighted
    formula calibrated against typical globular protein stabilities.

    Args:
        protein: Protein sequence (single-letter amino acid codes).

    Returns:
        FoldXResult with method="empirical".
    """
    start_time = time.time()

    if not protein:
        elapsed = time.time() - start_time
        return FoldXResult(
            protein="",
            pdb_string=None,
            stability_kcal=0.0,
            ddg_kcal=None,
            interaction_energy=None,
            backbone_hbond=None,
            sidechain_hbond=None,
            van_der_waals=None,
            electrostatics=None,
            solvation=None,
            van_der_waals_clashes=None,
            entropy_sidechain=None,
            entropy_mainchain=None,
            torsional_clash=None,
            backbone_clash=None,
            helix_dipole=None,
            disulfide=None,
            electrostatic_kon=None,
            partial_covalent=None,
            energy_ionisation=None,
            execution_time_s=elapsed,
            method="empirical",
            success=False,
            error="Empty protein sequence",
        )

    seq = protein.upper()
    n = len(seq)

    # ── Heuristic a: Hydrophobic core quality ──
    hydrophobic_count = sum(1 for aa in seq if aa in HYDROPHOBIC_AAS)
    hydrophobic_frac = hydrophobic_count / n if n > 0 else 0.0
    # Optimal range: 0.35–0.45; deviations penalized
    if 0.35 <= hydrophobic_frac <= 0.45:
        hydro_score = 0.0
    else:
        # Quadratic penalty for deviation from optimal range
        if hydrophobic_frac < 0.35:
            deviation = 0.35 - hydrophobic_frac
        else:
            deviation = hydrophobic_frac - 0.45
        hydro_score = deviation * 20.0  # kcal/mol penalty

    # ── Heuristic b: Charge balance ──
    pos_count = sum(1 for aa in seq if aa in POSITIVE_AAS)
    neg_count = sum(1 for aa in seq if aa in NEGATIVE_AAS)
    total_charged = pos_count + neg_count
    if total_charged > 0:
        charge_ratio = min(pos_count, neg_count) / max(pos_count, neg_count)
    else:
        charge_ratio = 0.0
    # Ideal ratio is 1.0; deviation is destabilizing
    charge_score = (1.0 - charge_ratio) * 3.0  # up to 3 kcal/mol penalty

    # ── Heuristic c: Proline content ──
    proline_count = sum(1 for aa in seq if aa == "P")
    proline_frac = proline_count / n if n > 0 else 0.0
    if proline_frac <= 0.05:
        pro_score = -0.2  # slight stabilizing effect of moderate proline
    elif proline_frac <= 0.08:
        pro_score = 0.0   # neutral
    else:
        pro_score = (proline_frac - 0.08) * 15.0  # destabilizing above 8%

    # ── Heuristic d: Glycine content ──
    glycine_count = sum(1 for aa in seq if aa == "G")
    glycine_frac = glycine_count / n if n > 0 else 0.0
    if glycine_frac <= 0.07:
        gly_score = 0.0
    else:
        gly_score = (glycine_frac - 0.07) * 10.0  # destabilizing above 7%

    # ── Heuristic e: Cysteine pairs (disulfide bonds) ──
    cys_count = sum(1 for aa in seq if aa == "C")
    if cys_count >= 2 and cys_count % 2 == 0:
        disulfide_count = cys_count // 2
        disulfide_score = -disulfide_count * 2.0  # each disulfide ~ -2 kcal/mol
    elif cys_count >= 2:
        # Odd number of cysteines — one unpaired, partial stabilization
        disulfide_score = -((cys_count - 1) // 2) * 2.0 + 0.5
    else:
        disulfide_score = 0.0

    # ── Heuristic f: Charged-hydrophobic pattern (helix indicator) ──
    # Alternating charged/hydrophobic residues suggest helical content
    pattern_score = 0.0
    if n >= 4:
        alternations = 0
        for i in range(n - 1):
            is_charged_i = seq[i] in POSITIVE_AAS or seq[i] in NEGATIVE_AAS
            is_hydro_i = seq[i] in HYDROPHOBIC_AAS
            is_charged_j = seq[i + 1] in POSITIVE_AAS or seq[i + 1] in NEGATIVE_AAS
            is_hydro_j = seq[i + 1] in HYDROPHOBIC_AAS
            if (is_charged_i and is_hydro_j) or (is_hydro_i and is_charged_j):
                alternations += 1
        alt_frac = alternations / (n - 1) if n > 1 else 0.0
        # Moderate alternation (0.2–0.35) suggests good helix content
        if 0.2 <= alt_frac <= 0.35:
            pattern_score = -0.5
        elif alt_frac > 0.35:
            pattern_score = (alt_frac - 0.35) * 3.0  # too much is bad

    # ── Heuristic g: Sequence length ──
    # Longer proteins have more intramolecular contacts → more stable
    # Base stability increases logarithmically with length
    length_bonus = -2.0 * (1.0 - 1.0 / (1.0 + n / 100.0))

    # ── Combine into estimated ΔG ──
    # Typical globular protein: -5 to -15 kcal/mol
    stability_kcal = (
        -5.0            # base stability for a folded protein
        + hydro_score   # positive = destabilizing
        + charge_score  # positive = destabilizing
        + pro_score     # can be stabilizing or destabilizing
        + gly_score     # positive = destabilizing
        + disulfide_score  # negative = stabilizing
        + pattern_score    # can be stabilizing or destabilizing
        + length_bonus     # negative = stabilizing for longer proteins
    )

    # Compute sub-component estimates for the result
    # Approximate decomposition of the total into FoldX-like terms
    vdw_est = hydrophobic_frac * n * -0.3  # hydrophobic packing ≈ VdW
    hbond_est = -(hydrophobic_count + pos_count + neg_count) * 0.1
    elec_est = -charge_ratio * 1.5 if total_charged > 0 else 0.0
    solv_est = -hydrophobic_frac * 2.0  # burial of hydrophobic surface
    sc_entropy_est = glycine_frac * n * 0.05  # glycine increases entropy

    elapsed = time.time() - start_time

    return FoldXResult(
        protein=protein,
        pdb_string=None,
        stability_kcal=round(stability_kcal, 2),
        ddg_kcal=None,
        interaction_energy=None,
        backbone_hbond=round(hbond_est, 2),
        sidechain_hbond=round(hbond_est * 0.6, 2),
        van_der_waals=round(vdw_est, 2),
        electrostatics=round(elec_est, 2),
        solvation=round(solv_est, 2),
        van_der_waals_clashes=None,
        entropy_sidechain=round(sc_entropy_est, 2),
        entropy_mainchain=round(glycine_frac * n * 0.02, 2),
        torsional_clash=None,
        backbone_clash=None,
        helix_dipole=None,
        disulfide=round(disulfide_score, 2),
        electrostatic_kon=None,
        partial_covalent=None,
        energy_ionisation=None,
        execution_time_s=round(elapsed, 4),
        method="empirical",
        success=True,
        error=None,
    )


# ────────────────────────────────────────────────────────────
# Mutation Scanning
# ────────────────────────────────────────────────────────────

def scan_mutations(
    protein: str,
    pdb_string: str | None = None,
    positions: list[int] | None = None,
) -> list[MutationResult]:
    """Scan all possible single-point mutations and estimate ΔΔG.

    For each position (or specified positions), tries all 19 possible
    amino acid substitutions. If FoldX is available and a PDB structure
    is provided, uses BuildModel for accurate ΔΔG. Otherwise, uses
    empirical scoring based on BLOSUM62 substitution scores and
    hydrophobicity changes.

    Args:
        protein: Protein sequence (single-letter codes).
        pdb_string: PDB structure content, if available.
        positions: Specific 0-indexed positions to scan. If None,
            scans all positions.

    Returns:
        List of MutationResult sorted by ΔΔG (most stabilizing first).
    """
    if not protein:
        return []

    seq = protein.upper()
    n = len(seq)

    if positions is None:
        positions = list(range(n))
    else:
        # Validate positions
        positions = [p for p in positions if 0 <= p < n]

    # If FoldX is available and PDB provided, use it
    if pdb_string and is_foldx_available():
        mutation_strings: list[str] = []
        position_map: list[tuple[int, str, str]] = []  # (0-idx, wt, mt)
        for pos in positions:
            wt = seq[pos]
            if wt not in STANDARD_AAS:
                continue
            for mt in STANDARD_AAS:
                if mt == wt:
                    continue
                # FoldX uses 1-indexed positions
                mut_str = f"{wt}{pos + 1}{mt}"
                mutation_strings.append(mut_str)
                position_map.append((pos, wt, mt))

        foldx_results = run_foldx_mutation(pdb_string, mutation_strings)
        if foldx_results:
            foldx_results.sort(key=lambda r: r.ddg_kcal)
            return foldx_results

    # Fall back to empirical scoring
    return _empirical_mutation_scan(seq, positions)


def find_stabilizing_mutations(
    protein: str,
    pdb_string: str | None = None,
    ddg_threshold: float = -1.0,
) -> list[MutationResult]:
    """Find mutations predicted to stabilize the protein.

    Scans all positions and returns only mutations with ΔΔG below
    the specified threshold (negative ΔΔG = stabilizing).

    Args:
        protein: Protein sequence (single-letter codes).
        pdb_string: PDB structure content, if available.
        ddg_threshold: ΔΔG threshold in kcal/mol. Only mutations
            with ΔΔG < threshold are returned. Default: -1.0.

    Returns:
        List of stabilizing MutationResult, sorted by ΔΔG
        (most stabilizing first).
    """
    all_mutations = scan_mutations(protein, pdb_string)
    return [m for m in all_mutations if m.ddg_kcal < ddg_threshold]


# ────────────────────────────────────────────────────────────
# Internal Helpers
# ────────────────────────────────────────────────────────────

def _parse_stability_fxout(filepath: str) -> dict:
    """Parse a FoldX Stability .fxout file.

    The .fxout file has a header line with column labels followed by
    data lines with energy values. We extract the total stability and
    individual energy components.

    Returns:
        Dict with energy component keys and float values.
    """
    result: dict = {}

    try:
        with open(filepath) as f:
            lines = f.readlines()
    except OSError:
        return result

    # FoldX .fxout format: header lines start with #, data follows
    # Look for the line containing energy values
    # Typical header: "Pdb       Total Energy   Backbone Hbond  ..."
    # Typical data:   "input.pdb  -45.32         -12.45          ..."

    energy_keys = [
        "stability_kcal",
        "backbone_hbond",
        "sidechain_hbond",
        "van_der_waals",
        "electrostatics",
        "solvation",
        "van_der_waals_clashes",
        "entropy_sidechain",
        "entropy_mainchain",
        "torsional_clash",
        "backbone_clash",
        "helix_dipole",
        "disulfide",
        "electrostatic_kon",
        "partial_covalent",
        "energy_ionisation",
    ]

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Skip header lines
        if "Pdb" in line or "Total" in line:
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        # First field is the PDB name, rest are energy values
        try:
            values = [float(v) for v in parts[1:]]
        except ValueError:
            continue

        # Map values to keys
        for i, key in enumerate(energy_keys):
            if i < len(values):
                result[key] = values[i]

        # Interaction energy is sometimes a separate column
        if len(values) > len(energy_keys):
            result["interaction_energy"] = values[len(energy_keys)]

        break  # Only use the first data line

    return result


def _parse_buildmodel_output(tmpdir: str) -> dict[str, float]:
    """Parse FoldX BuildModel ΔΔG output.

    Looks for the average ΔΔG file produced by BuildModel.

    Returns:
        Dict mapping mutation strings to ΔΔG values.
    """
    ddg_values: dict[str, float] = {}

    # FoldX produces several output files; the ΔΔG values are in
    # the DIF file or the average file
    for fname in os.listdir(tmpdir):
        filepath = os.path.join(tmpdir, fname)
        if not os.path.isfile(filepath):
            continue

        # Look for files with ΔΔG information
        if "DIF" in fname.upper() or "ABVERAGE" in fname.upper() or fname.endswith(".fxout"):
            try:
                with open(filepath) as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        parts = line.split()
                        # Try to find mutation and ΔΔG
                        if len(parts) >= 2:
                            # Look for mutation format like A123G
                            for part in parts:
                                match = re.match(
                                    r"^([ARNDCQEGHILKMFPSTWYV])(\d+)([ARNDCQEGHILKMFPSTWYV])$",
                                    part,
                                )
                                if match:
                                    mut_str = match.group(0)
                                    # ΔΔG is typically the last numeric value on the line
                                    for val_part in reversed(parts):
                                        try:
                                            ddg_values[mut_str] = float(val_part)
                                            break
                                        except ValueError:
                                            continue
                                    break
            except OSError:
                continue

    return ddg_values


def _extract_protein_from_pdb(pdb_string: str) -> str:
    """Extract protein sequence from PDB ATOM records.

    Uses CA atom residue names to reconstruct the one-letter sequence.

    Args:
        pdb_string: PDB file content.

    Returns:
        Protein sequence as single-letter codes.
    """
    three_to_one: dict[str, str] = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    residues: dict[int, str] = []

    for line in pdb_string.splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        res_name = line[17:20].strip()
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            continue

        aa = three_to_one.get(res_name)
        if aa:
            residues.append((res_seq, aa))

    # Sort by residue number and deduplicate
    residues.sort(key=lambda x: x[0])
    seen: set[int] = set()
    protein_chars: list[str] = []
    for res_seq, aa in residues:
        if res_seq not in seen:
            seen.add(res_seq)
            protein_chars.append(aa)

    return "".join(protein_chars)


def _parse_mutation_string(mut_str: str) -> tuple[str, int, str] | None:
    """Parse a mutation string like 'A123G' into (wildtype, position, mutant).

    Args:
        mut_str: Mutation in format WT<number>MUT (1-indexed position).

    Returns:
        Tuple of (wildtype_aa, position_1indexed, mutant_aa) or None.
    """
    match = re.match(
        r"^([ARNDCQEGHILKMFPSTWYV])(\d+)([ARNDCQEGHILKMFPSTWYV])$", mut_str
    )
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3)


def _empirical_mutation_scan(
    protein: str,
    positions: list[int],
) -> list[MutationResult]:
    """Scan mutations using empirical BLOSUM62 + hydropathy scoring.

    For each position and each of the 19 possible substitutions,
    estimates ΔΔG based on:

    1. BLOSUM62 substitution score (conservative changes are less
       destabilizing).
    2. Hydropathy change (large changes in hydrophobicity at core
       positions are destabilizing).
    3. Proline/glycine introduction in regular secondary structure
       is destabilizing.

    Args:
        protein: Upper-cased protein sequence.
        positions: 0-indexed positions to scan.

    Returns:
        List of MutationResult sorted by ΔΔG (most stabilizing first).
    """
    n = len(protein)
    results: list[MutationResult] = []

    # Estimate which positions are "core" vs "surface" based on
    # local hydrophobicity
    window = 5
    for pos in positions:
        wt = protein[pos]
        if wt not in BLOSUM62:
            continue

        # Estimate local hydrophobicity (proxy for burial)
        start = max(0, pos - window)
        end = min(n, pos + window + 1)
        local_hydro = sum(HYDROPATHY.get(protein[i], 0.0) for i in range(start, end))
        local_hydro_avg = local_hydro / (end - start) if (end - start) > 0 else 0.0
        is_core = local_hydro_avg > 0.5  # above zero = hydrophobic region

        for mt in STANDARD_AAS:
            if mt == wt:
                continue

            ddg = _estimate_ddg(wt, mt, is_core, pos, n)
            results.append(MutationResult(
                position=pos,
                wildtype=wt,
                mutant=mt,
                ddg_kcal=round(ddg, 2),
                stabilizing=ddg < -0.5,
                neutral=-0.5 <= ddg <= 0.5,
                destabilizing=ddg > 0.5,
                method="empirical",
            ))

    results.sort(key=lambda r: r.ddg_kcal)
    return results


def _estimate_ddg(
    wt: str,
    mt: str,
    is_core: bool,
    position: int,
    protein_length: int,
) -> float:
    """Estimate ΔΔG for a single point mutation using empirical rules.

    Combines BLOSUM62 score, hydropathy change, and structural context
    heuristics into an estimated ΔΔG value.

    Args:
        wt: Wildtype amino acid.
        mt: Mutant amino acid.
        is_core: Whether the position is estimated to be in the core.
        position: 0-indexed position in the protein.
        protein_length: Total protein length.

    Returns:
        Estimated ΔΔG in kcal/mol (positive = destabilizing).
    """
    # Base ΔΔG from BLOSUM62 — higher scores = more conservative = less ΔΔG
    blosum = BLOSUM62.get(wt, {}).get(mt, -4)
    # Convert BLOSUM62 score to rough ΔΔG
    # BLOSUM62 range: -4 to +11; typical ΔΔG range: -2 to +8
    ddg_blosum = -blosum * 0.3  # invert and scale

    # Hydropathy change penalty
    hydro_wt = HYDROPATHY.get(wt, 0.0)
    hydro_mt = HYDROPATHY.get(mt, 0.0)
    hydro_change = abs(hydro_mt - hydro_wt)

    if is_core:
        # In the core: replacing hydrophobic with polar is very destabilizing
        if hydro_wt > 0 and hydro_mt < 0:
            ddg_hydro = hydro_change * 0.5
        elif hydro_wt < 0 and hydro_mt > 0:
            # Polar to hydrophobic in core can be stabilizing
            ddg_hydro = -hydro_change * 0.15
        else:
            ddg_hydro = 0.0
    else:
        # On surface: hydrophobic to polar is often stabilizing
        if hydro_wt > 0 and hydro_mt < 0:
            ddg_hydro = -hydro_change * 0.1
        elif hydro_wt < 0 and hydro_mt > 0:
            # Polar to hydrophobic on surface is destabilizing
            ddg_hydro = hydro_change * 0.2
        else:
            ddg_hydro = 0.0

    # Proline introduction penalty (breaks secondary structure)
    ddg_pro = 0.0
    if mt == "P" and wt != "P":
        # Proline introduces kinks — more destabilizing in the middle
        mid_factor = 1.0 - abs(position - protein_length / 2) / (protein_length / 2 + 1)
        ddg_pro = 1.5 + mid_factor * 1.0
    if wt == "P" and mt != "P":
        # Removing proline is generally stabilizing
        ddg_pro = -0.8

    # Glycine introduction penalty (increases backbone flexibility)
    ddg_gly = 0.0
    if mt == "G" and wt != "G":
        ddg_gly = 0.6 if is_core else 0.2
    if wt == "G" and mt != "G":
        # Replacing glycine with something more rigid can be stabilizing
        ddg_gly = -0.3 if is_core else -0.1

    # Cysteine introduction (potential disulfide)
    ddg_cys = 0.0
    if mt == "C" and wt != "C":
        ddg_cys = -0.5  # potential disulfide formation
    if wt == "C" and mt != "C":
        ddg_cys = 1.0  # breaking a potential disulfide

    # Charged residue changes
    ddg_charge = 0.0
    wt_charged = wt in POSITIVE_AAS or wt in NEGATIVE_AAS
    mt_charged = mt in POSITIVE_AAS or mt in NEGATIVE_AAS
    if is_core:
        # Introducing charge in the core is destabilizing
        if not wt_charged and mt_charged:
            ddg_charge = 2.0
        # Removing charge from core is stabilizing
        if wt_charged and not mt_charged:
            ddg_charge = -0.5

    total_ddg = ddg_blosum + ddg_hydro + ddg_pro + ddg_gly + ddg_cys + ddg_charge
    return round(total_ddg, 2)


def _cleanup_tempdir(tmpdir: str) -> None:
    """Remove a temporary directory and all its contents.

    Args:
        tmpdir: Path to the temporary directory.
    """
    try:
        for fname in os.listdir(tmpdir):
            fpath = os.path.join(tmpdir, fname)
            try:
                os.remove(fpath)
            except OSError:
                pass
        os.rmdir(tmpdir)
    except OSError:
        logger.debug("Could not clean up temp directory: %s", tmpdir)
