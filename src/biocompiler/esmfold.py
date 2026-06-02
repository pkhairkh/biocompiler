"""
BioCompiler ESMFold Client — Protein Structure Prediction

Provides a client for the ESMFold protein structure prediction API
(ESM Atlas) with local esm package fallback and offline graceful
degradation.  Parses PDB output to extract pLDDT confidence scores,
backbone dihedrals, and contact maps.

Pipeline integration point:
  After translation (protein sequence) → predict 3-D structure →
  confidence assessment → certificate enrichment.

All network calls include retry logic with exponential backoff.
Individual batch failures are isolated — they never crash sibling
predictions.
"""

from __future__ import annotations

import logging
import math
import re
import time
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Semaphore
from typing import Any

from .constants import CODON_TABLE, AA_TO_CODONS
from .exceptions import BioCompilerError

logger = logging.getLogger(__name__)

# ==============================================================================
# Constants
# ==============================================================================

STANDARD_AMINO_ACIDS: set[str] = set("ACDEFGHIKLMNPQRSTVWY")

DEFAULT_API_URL = "https://api.esmatlas.com/fetchPredictedStructure"
DEFAULT_TIMEOUT = 120.0
MAX_RETRIES = 3
RETRY_BASE_DELAY = 2.0  # seconds, doubled each attempt

# PDB ATOM record regex — captures key fields from fixed-width format
# Columns:  1-6  Record name  7-11  Serial  13-16  Atom name
#           17    Alt loc     18-20  Res name  22     Chain
#           23-26 Res seq     27     iCode
#           31-38 X           39-46  Y         47-54  Z
#           55-60 Occupancy  61-66  B-factor (pLDDT)
_PDB_ATOM_RE = re.compile(
    r"^ATOM\s+"
    r"(\d+)\s+"           # serial
    r"(.{1,4})\s+"        # atom name (may have leading space)
    r"(\w)\s+"            # alt loc
    r"(\w{3})\s+"         # residue name
    r"(\w)\s+"            # chain ID
    r"(\d+)\s+"           # residue sequence number
    r"(\w)?\s*"           # insertion code
    r"([-]?\d+\.\d+)\s+"  # x
    r"([-]?\d+\.\d+)\s+"  # y
    r"([-]?\d+\.\d+)\s+"  # z
    r"(\d+\.\d+)\s+"      # occupancy
    r"(\d+\.\d+)"         # b-factor (pLDDT)
)


# ==============================================================================
# Exceptions
# ==============================================================================

class ESMFoldError(BioCompilerError):
    """Raised when ESMFold prediction fails."""

    def __init__(self, reason: str, protein: str | None = None):
        self.reason = reason
        self.protein = protein
        msg = f"ESMFold prediction failed: {reason}"
        if protein:
            msg += f" (protein length={len(protein)})"
        super().__init__(msg)


# ==============================================================================
# Data classes
# ==============================================================================

@dataclass
class ESMFoldResult:
    """Result of an ESMFold protein structure prediction.

    Attributes:
        pdb_string:       Predicted structure in PDB format.
        plddt_scores:     Per-residue pLDDT confidence scores (0-100).
        mean_plddt:       Average pLDDT across all residues.
        pae_matrix:       Predicted Aligned Error matrix (residue×residue).
                          None when the API does not return PAE.
        protein:          Input protein sequence (single-letter codes).
        model_name:       ESMFold model identifier (e.g. "esmfold_v1").
        execution_time_s: Wall-clock time for the prediction in seconds.
        success:          Whether prediction completed without error.
        error:            Error message if success is False, else None.
    """

    pdb_string: str
    plddt_scores: list[float]
    mean_plddt: float
    pae_matrix: list[list[float]] | None
    protein: str
    model_name: str
    execution_time_s: float
    success: bool
    error: str | None = None


# ==============================================================================
# Input validation
# ==============================================================================

def _validate_protein(protein: str) -> None:
    """Raise ESMFoldError if *protein* contains non-standard amino acids.

    Only the 20 canonical single-letter codes are accepted.
    """
    if not protein:
        raise ESMFoldError("Empty protein sequence", protein=protein)
    invalid = set(protein) - STANDARD_AMINO_ACIDS
    if invalid:
        raise ESMFoldError(
            f"Invalid amino acid(s) in protein: {invalid}. "
            f"Only standard single-letter codes {sorted(STANDARD_AMINO_ACIDS)} are permitted.",
            protein=protein,
        )


# ==============================================================================
# Availability check
# ==============================================================================

def is_esmfold_available() -> bool:
    """Check whether ESMFold is available (API reachable or local esm installed).

    Returns True if either:
      - The ESM Atlas API responds to a HEAD request, or
      - The ``esm`` Python package is importable.

    Returns False otherwise (network error, no local esm, etc.).
    """
    # 1. Try the remote API
    try:
        import urllib.request
        import urllib.error

        req = urllib.request.Request(
            DEFAULT_API_URL, method="HEAD"
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            if resp.status < 500:
                logger.debug("ESM Atlas API reachable (status %d)", resp.status)
                return True
    except Exception as exc:
        logger.debug("ESM Atlas API not reachable: %s", exc)

    # 2. Try local esm package
    try:
        import esm  # noqa: F401
        logger.debug("Local esm package is importable")
        return True
    except ImportError:
        pass

    return False


# ==============================================================================
# PDB parsing helpers
# ==============================================================================

def parse_pdb(pdb_string: str) -> dict[str, Any]:
    """Parse a PDB format string into a structured dictionary.

    Extracts atom records, residue information, chain identifiers, CA
    coordinates, and B-factor (pLDDT) scores.

    Returns:
        {
            "atoms":     list[dict]  — one per ATOM record,
            "residues":  list[dict]  — one per unique residue,
            "chains":    list[str]   — unique chain identifiers,
            "plddt_scores": list[float] — per-residue pLDDT (from CA B-factor),
        }
    """
    atoms: list[dict[str, Any]] = []
    residues: list[dict[str, Any]] = []
    chains: set[str] = set()
    plddt_scores: list[float] = []

    seen_residues: dict[tuple[str, int, str], dict[str, Any]] = {}

    for line in pdb_string.splitlines():
        if not line.startswith("ATOM"):
            continue

        # Fixed-width PDB parsing
        try:
            atom_name = line[12:16].strip()
            alt_loc = line[16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21].strip() if len(line) > 21 else ""
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            occupancy = float(line[54:60].strip()) if len(line) > 59 else 1.0
            b_factor = float(line[60:66].strip()) if len(line) > 65 else 0.0
        except (ValueError, IndexError) as exc:
            logger.debug("Skipping unparseable ATOM line: %r (%s)", line[:66], exc)
            continue

        atom_record = {
            "atom_name": atom_name,
            "alt_loc": alt_loc,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_seq": res_seq,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "b_factor": b_factor,
        }
        atoms.append(atom_record)
        chains.add(chain_id)

        # Track unique residues
        res_key = (chain_id, res_seq, res_name)
        if res_key not in seen_residues:
            res_dict: dict[str, Any] = {
                "chain_id": chain_id,
                "res_seq": res_seq,
                "res_name": res_name,
                "ca": None,  # will be filled from CA atom
                "b_factor": b_factor,
                "atom_names": [atom_name],
            }
            seen_residues[res_key] = res_dict
        else:
            seen_residues[res_key]["atom_names"].append(atom_name)

        # Capture CA atom B-factor as pLDDT
        if atom_name == "CA" and alt_loc in ("", "A"):
            if seen_residues[res_key]["ca"] is None:
                seen_residues[res_key]["ca"] = (x, y, z)
                seen_residues[res_key]["b_factor"] = b_factor

    # Build residue list and pLDDT scores in sequence order
    for res_key in sorted(seen_residues.keys()):
        res = seen_residues[res_key]
        residues.append(res)
        plddt_scores.append(res["b_factor"])

    return {
        "atoms": atoms,
        "residues": residues,
        "chains": sorted(chains),
        "plddt_scores": plddt_scores,
    }


# ==============================================================================
# Backbone dihedral computation
# ==============================================================================

def compute_backbone_dihedrals(pdb_string: str) -> dict[str, list[float | None]]:
    """Compute phi (φ) and psi (ψ) backbone dihedral angles from a PDB structure.

    Uses CA coordinates as an approximation (ideal phi/psi requires N, CA, C
    atoms; here we approximate from CA traces for speed and robustness).

    Returns:
        {"phi": [float|None, ...], "psi": [float|None, ...]}
        Angles are in degrees.  First phi and last psi are None
        (undefined at chain termini).
    """
    parsed = parse_pdb(pdb_string)
    residues = parsed["residues"]

    # Collect CA coordinates
    ca_coords: list[tuple[float, float, float]] = []
    for res in residues:
        if res["ca"] is not None:
            ca_coords.append(res["ca"])

    n = len(ca_coords)
    phi: list[float | None] = [None] * n
    psi: list[float | None] = [None] * n

    if n < 4:
        # Not enough residues for dihedrals
        return {"phi": phi, "psi": psi}

    for i in range(1, n - 2):
        # Approximate phi(i) from CA_{i-1}, CA_i, CA_{i+1}, CA_{i+2}
        angle = _dihedral_angle(
            ca_coords[i - 1], ca_coords[i], ca_coords[i + 1], ca_coords[i + 2]
        )
        if angle is not None:
            phi[i] = math.degrees(angle)

    for i in range(1, n - 2):
        # Approximate psi(i) from CA_i, CA_{i+1}, CA_{i+2}, CA_{i+3}
        # but we use the more conventional offset:
        # psi(i) ~ dihedral(CA_{i}, CA_{i+1}, CA_{i+2}) mapped
        angle = _dihedral_angle(
            ca_coords[i], ca_coords[i + 1],
            ca_coords[min(i + 2, n - 1)],
            ca_coords[min(i + 3, n - 1)],
        )
        if angle is not None and i < n - 1:
            psi[i] = math.degrees(angle)

    return {"phi": phi, "psi": psi}


def _dihedral_angle(
    p0: tuple[float, float, float],
    p1: tuple[float, float, float],
    p2: tuple[float, float, float],
    p3: tuple[float, float, float],
) -> float | None:
    """Compute the dihedral angle defined by four 3-D points (radians).

    Returns None if any two points are coincident (degenerate).
    """
    import numpy as np  # type: ignore[import-untyped]

    b0 = np.array(p1) - np.array(p0)
    b1 = np.array(p2) - np.array(p1)
    b2 = np.array(p3) - np.array(p2)

    # Normalise b1 so that it does not influence magnitude of vector
    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-8:
        return None
    b1_unit = b1 / b1_norm

    # v = projection of b0 onto plane perpendicular to b1
    v = b0 - np.dot(b0, b1_unit) * b1_unit
    w = b2 - np.dot(b2, b1_unit) * b1_unit

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_unit, v), w)

    return math.atan2(y, x)


# ==============================================================================
# pLDDT classification
# ==============================================================================

def classify_plddt(mean_plddt: float) -> str:
    """Classify a mean pLDDT score into a confidence category.

    Thresholds follow the AlphaFold / ESMFold convention:

        > 90 : "Very high (experimental)"
        70-90: "Confident"
        50-70: "Low confidence"
        < 50 : "Very low"

    Args:
        mean_plddt: Average per-residue pLDDT score (0-100).

    Returns:
        Human-readable confidence classification string.
    """
    if mean_plddt > 90:
        return "Very high (experimental)"
    elif mean_plddt > 70:
        return "Confident"
    elif mean_plddt > 50:
        return "Low confidence"
    else:
        return "Very low"


# ==============================================================================
# Contact map estimation
# ==============================================================================

def estimate_contact_map(
    pdb_string: str,
    distance_threshold: float = 8.0,
) -> list[list[int]]:
    """Estimate a binary residue contact map from CA coordinates.

    Two residues are in contact if the Euclidean distance between their
    Cα atoms is below *distance_threshold* (default 8.0 Å).

    Args:
        pdb_string:        PDB format structure string.
        distance_threshold: Distance cutoff in Ångströms.

    Returns:
        Square binary matrix (list of lists).  entry[i][j] == 1 if
        residues i and j are in contact, 0 otherwise.  Diagonal is 0.
    """
    parsed = parse_pdb(pdb_string)
    ca_coords: list[tuple[float, float, float]] = []
    for res in parsed["residues"]:
        if res["ca"] is not None:
            ca_coords.append(res["ca"])

    n = len(ca_coords)
    contact_map: list[list[int]] = [[0] * n for _ in range(n)]

    for i in range(n):
        xi, yi, zi = ca_coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = ca_coords[j]
            dist = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            if dist < distance_threshold:
                contact_map[i][j] = 1
                contact_map[j][i] = 1

    return contact_map


# ==============================================================================
# Core prediction — single protein
# ==============================================================================

def predict_structure(
    protein: str,
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> ESMFoldResult:
    """Predict the 3-D structure of a protein using ESMFold.

    Strategy priority:
      1. **API** — POST the sequence to ESM Atlas and parse the PDB response.
      2. **Local esm** — ``import esm`` and run ESMFold locally (if installed).
      3. **Offline** — Return ``ESMFoldResult(success=False)`` with an error.

    Retry logic: up to 3 attempts with exponential backoff on transient
    API errors (network, 5xx, 429 rate-limit).  Invalid protein input
    (non-standard AA) causes immediate failure with no retry.

    Args:
        protein: Protein sequence (single-letter amino acid codes).
        use_api: If True, try the remote API first.
        api_url: ESM Atlas API endpoint.
        timeout: HTTP request timeout in seconds.

    Returns:
        ESMFoldResult with PDB string, pLDDT scores, and metadata.

    Raises:
        ESMFoldError: If the protein contains invalid amino acid codes.
    """
    _validate_protein(protein)

    t0 = time.monotonic()

    # --- Strategy 1: Remote API ------------------------------------------------
    if use_api:
        result = _predict_via_api(protein, api_url, timeout)
        if result is not None:
            result.execution_time_s = time.monotonic() - t0
            return result

    # --- Strategy 2: Local esm package -----------------------------------------
    result = _predict_via_local_esm(protein)
    if result is not None:
        result.execution_time_s = time.monotonic() - t0
        return result

    # --- Strategy 3: Offline fallback ------------------------------------------
    elapsed = time.monotonic() - t0
    return ESMFoldResult(
        pdb_string="",
        plddt_scores=[],
        mean_plddt=0.0,
        pae_matrix=None,
        protein=protein,
        model_name="esmfold_v1",
        execution_time_s=elapsed,
        success=False,
        error="ESMFold unavailable: API unreachable and local esm package not installed",
    )


def _predict_via_api(
    protein: str,
    api_url: str,
    timeout: float,
) -> ESMFoldResult | None:
    """Attempt ESMFold prediction via the ESM Atlas API with retries.

    Returns None on exhaustion so the caller can try the next strategy.
    """
    import urllib.request
    import urllib.error
    import json

    last_error: str | None = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            data = protein.encode("utf-8")
            req = urllib.request.Request(
                api_url,
                data=data,
                headers={"Content-Type": "application/x-fasta"},
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                pdb_string = resp.read().decode("utf-8")

            if not pdb_string or len(pdb_string) < 50:
                last_error = "API returned empty or truncated PDB"
                logger.warning("Attempt %d/%d: %s", attempt, MAX_RETRIES, last_error)
                continue

            return _build_result_from_pdb(protein, pdb_string, "esmfold_v1")

        except urllib.error.HTTPError as exc:
            if exc.code == 429:
                # Rate limited — respect Retry-After or use backoff
                retry_after = exc.headers.get("Retry-After")
                wait = float(retry_after) if retry_after else RETRY_BASE_DELAY * (2 ** (attempt - 1))
                last_error = f"API rate limited (429), waiting {wait:.1f}s"
                logger.warning("Attempt %d/%d: %s", attempt, MAX_RETRIES, last_error)
                time.sleep(wait)
                continue
            elif exc.code >= 500:
                wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
                last_error = f"API server error ({exc.code}), retrying in {wait:.1f}s"
                logger.warning("Attempt %d/%d: %s", attempt, MAX_RETRIES, last_error)
                time.sleep(wait)
                continue
            elif exc.code >= 400:
                # Client error (bad request, etc.) — do not retry
                last_error = f"API client error ({exc.code}): {exc.read().decode('utf-8', errors='replace')[:200]}"
                logger.error("API client error, not retrying: %s", last_error)
                break
            else:
                last_error = f"HTTP error {exc.code}"
                logger.warning("Attempt %d/%d: %s", attempt, MAX_RETRIES, last_error)
                continue

        except urllib.error.URLError as exc:
            wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
            last_error = f"Network error: {exc.reason}"
            logger.warning("Attempt %d/%d: %s (retry in {wait:.1f}s)", attempt, MAX_RETRIES, last_error)
            time.sleep(wait)
            continue

        except TimeoutError:
            wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
            last_error = f"Request timed out after {timeout}s"
            logger.warning("Attempt %d/%d: %s", attempt, MAX_RETRIES, last_error)
            time.sleep(wait)
            continue

        except Exception as exc:
            last_error = f"Unexpected error: {exc}"
            logger.error("Unexpected API error: %s", exc)
            break

    # All retries exhausted
    if last_error is not None:
        logger.error("All %d API attempts failed. Last error: %s", MAX_RETRIES, last_error)
    return None


def _predict_via_local_esm(protein: str) -> ESMFoldResult | None:
    """Attempt ESMFold prediction using the locally-installed esm package.

    Returns None if esm is not installed, or if prediction fails.
    """
    try:
        import esm
        import torch  # type: ignore[import-untyped]
    except ImportError:
        logger.debug("Local esm/torch not available")
        return None

    try:
        logger.info("Running ESMFold locally via esm package")
        model, alphabet = esm.pretrained.esmfold_v1()
        model = model.eval()
        if torch.cuda.is_available():
            model = model.cuda()

        batch_converter = alphabet.get_batch_converter()
        data = [("protein", protein)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        if torch.cuda.is_available():
            batch_tokens = batch_tokens.cuda()

        with torch.no_grad():
            output = model(batch_tokens)

        pdb_string = model.output_to_pdb(output)[0]
        mean_plddt = float(output["mean_plddt"][0])

        # Per-residue pLDDT from atom-level outputs
        plddt_per_residue = output.get("plddt")
        if plddt_per_residue is not None:
            plddt_scores = plddt_per_residue[0, :len(protein)].cpu().tolist()
        else:
            plddt_scores = [mean_plddt] * len(protein)

        return ESMFoldResult(
            pdb_string=pdb_string,
            plddt_scores=plddt_scores,
            mean_plddt=mean_plddt,
            pae_matrix=None,
            protein=protein,
            model_name="esmfold_v1",
            execution_time_s=0.0,  # filled by caller
            success=True,
        )

    except Exception as exc:
        logger.error("Local ESMFold prediction failed: %s", exc)
        return None


def _build_result_from_pdb(
    protein: str,
    pdb_string: str,
    model_name: str,
) -> ESMFoldResult:
    """Construct an ESMFoldResult by parsing a PDB string."""
    parsed = parse_pdb(pdb_string)
    plddt_scores = parsed["plddt_scores"]
    mean_plddt = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0

    return ESMFoldResult(
        pdb_string=pdb_string,
        plddt_scores=plddt_scores,
        mean_plddt=mean_plddt,
        pae_matrix=None,
        protein=protein,
        model_name=model_name,
        execution_time_s=0.0,  # filled by caller
        success=True,
    )


# ==============================================================================
# Batch prediction
# ==============================================================================

def predict_structure_batch(
    proteins: list[str],
    max_concurrent: int = 3,
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> list[ESMFoldResult]:
    """Predict structures for multiple proteins concurrently.

    Uses a :class:`ThreadPoolExecutor` with a semaphore for rate
    limiting.  Individual prediction failures are captured in the
    returned ``ESMFoldResult`` objects — they never crash sibling
    predictions.

    Args:
        proteins:      List of protein sequences (single-letter codes).
        max_concurrent: Maximum number of concurrent API requests.
        use_api:       If True, try the remote API first.
        api_url:       ESM Atlas API endpoint.
        timeout:       Per-protein request timeout in seconds.

    Returns:
        List of ESMFoldResult objects, one per input protein, in the
        same order as *proteins*.
    """
    semaphore = Semaphore(max_concurrent)
    results: dict[int, ESMFoldResult] = {}

    def _predict_with_semaphore(idx: int, protein: str) -> tuple[int, ESMFoldResult]:
        with semaphore:
            try:
                result = predict_structure(
                    protein,
                    use_api=use_api,
                    api_url=api_url,
                    timeout=timeout,
                )
            except ESMFoldError as exc:
                result = ESMFoldResult(
                    pdb_string="",
                    plddt_scores=[],
                    mean_plddt=0.0,
                    pae_matrix=None,
                    protein=protein,
                    model_name="esmfold_v1",
                    execution_time_s=0.0,
                    success=False,
                    error=str(exc),
                )
            except Exception as exc:
                result = ESMFoldResult(
                    pdb_string="",
                    plddt_scores=[],
                    mean_plddt=0.0,
                    pae_matrix=None,
                    protein=protein,
                    model_name="esmfold_v1",
                    execution_time_s=0.0,
                    success=False,
                    error=f"Unexpected error: {exc}",
                )
            return idx, result

    with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
        futures = {
            executor.submit(_predict_with_semaphore, i, prot): i
            for i, prot in enumerate(proteins)
        }
        for future in as_completed(futures):
            idx, result = future.result()
            results[idx] = result

    # Preserve input order
    return [results[i] for i in range(len(proteins))]


# ==============================================================================
# Convenience: full structure analysis
# ==============================================================================

def analyze_structure(
    protein: str,
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> dict[str, Any]:
    """Run a full structure analysis on a single protein.

    Combines prediction, PDB parsing, dihedral computation, contact
    map estimation, and pLDDT classification into one call.

    Args:
        protein: Protein sequence.
        use_api: Try the remote API first.
        api_url: ESM Atlas endpoint.
        timeout: Request timeout.

    Returns:
        Dict with keys: result (ESMFoldResult), parsed (dict),
        dihedrals (dict), contact_map (list[list[int]]),
        confidence (str).
    """
    result = predict_structure(protein, use_api=use_api, api_url=api_url, timeout=timeout)

    if not result.success:
        return {
            "result": result,
            "parsed": {},
            "dihedrals": {"phi": [], "psi": []},
            "contact_map": [],
            "confidence": "N/A",
        }

    parsed = parse_pdb(result.pdb_string)
    dihedrals = compute_backbone_dihedrals(result.pdb_string)
    contact_map = estimate_contact_map(result.pdb_string)
    confidence = classify_plddt(result.mean_plddt)

    return {
        "result": result,
        "parsed": parsed,
        "dihedrals": dihedrals,
        "contact_map": contact_map,
        "confidence": confidence,
    }
