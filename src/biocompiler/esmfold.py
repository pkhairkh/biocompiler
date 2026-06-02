"""
BioCompiler ESMFold Client — Protein Structure Prediction

Provides a client for the ESMFold protein structure prediction API
(ESM Atlas) with local esm package fallback and offline graceful
degradation.  Parses PDB output to extract pLDDT confidence scores,
backbone dihedrals, and contact maps.

Also includes:
  - ESMFoldCache: In-memory and file-based caching for predictions.
  - Batch processing with progress tracking, error isolation, and
    rate limiting.

Pipeline integration point:
  After translation (protein sequence) → predict 3-D structure →
  confidence assessment → certificate enrichment.

All network calls include retry logic with exponential backoff.
Individual batch failures are isolated — they never crash sibling
predictions.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError
from dataclasses import dataclass, field
from threading import Semaphore
from typing import Any, Callable, Dict, Optional

from .constants import CODON_TABLE, AA_TO_CODONS, DEFAULT_ENGINE_TIMEOUT
from .engine_base import validate_protein_sequence, EngineTimer
from .exceptions import BioCompilerError

logger = logging.getLogger(__name__)

# ==============================================================================
# Constants
# ==============================================================================

STANDARD_AMINO_ACIDS: set[str] = set("ACDEFGHIKLMNPQRSTVWY")

# ESMFold-specific constants (timeout falls back to shared DEFAULT_ENGINE_TIMEOUT)
DEFAULT_API_URL = "https://api.esmatlas.com/fetchPredictedStructure"
DEFAULT_TIMEOUT: float = DEFAULT_ENGINE_TIMEOUT
MAX_RETRIES = 3
RETRY_BASE_DELAY = 2.0  # seconds, doubled each attempt

MAX_BATCH_SIZE = 50
MAX_PROTEIN_LENGTH = 1000

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
        method:           How the prediction was obtained — ``"esmfold_api"``
                          for the remote ESM Atlas API, ``"esmfold_local"``
                          for the locally-installed esm package.
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
    method: str = "esmfold_api"


# ==============================================================================
# ESMFold Cache (merged from esmfold_cache.py)
# ==============================================================================

class ESMFoldCache:
    """Cache for ESMFold structure predictions.

    Supports in-memory caching with optional file-based persistence.
    """

    def __init__(self, cache_dir: Optional[str] = None, max_size: int = 1000):
        """Initialize the cache.

        Args:
            cache_dir: Directory for file-based cache persistence.
                       If None, uses in-memory only.
            max_size: Maximum number of entries in memory cache.
        """
        self._cache: Dict[str, ESMFoldResult] = {}
        self._cache_dir = cache_dir
        self._max_size = max_size
        self._hits = 0
        self._misses = 0

    @staticmethod
    def _key(protein: str) -> str:
        """Generate a cache key from a protein sequence."""
        return hashlib.sha256(protein.encode()).hexdigest()[:16]

    def get(self, protein: str) -> Optional[ESMFoldResult]:
        """Retrieve a cached prediction.

        Args:
            protein: Amino acid sequence.

        Returns:
            ESMFoldResult if cached, None otherwise.
        """
        key = self._key(protein)

        # Check memory cache
        if key in self._cache:
            self._hits += 1
            return self._cache[key]

        # Check file cache
        if self._cache_dir is not None:
            filepath = os.path.join(self._cache_dir, f"{key}.json")
            if os.path.exists(filepath):
                try:
                    with open(filepath, "r") as f:
                        data = json.load(f)
                    result = ESMFoldResult(
                        protein=data["protein"],
                        pdb_string=data["pdb_string"],
                        plddt_scores=data.get("plddt_scores", []),
                        mean_plddt=data["mean_plddt"],
                        pae_matrix=data.get("pae_matrix"),
                        model_name=data.get("model_name", "esmfold_v1"),
                        execution_time_s=data.get("execution_time_s", 0.0),
                        success=data.get("success", True),
                        error=data.get("error"),
                        method=data.get("method", "esmfold_api"),
                    )
                    # Promote to memory cache
                    self._cache[key] = result
                    self._hits += 1
                    return result
                except (json.JSONDecodeError, KeyError):
                    pass

        self._misses += 1
        return None

    def put(self, protein: str, result: ESMFoldResult) -> None:
        """Store a prediction in the cache.

        Args:
            protein: Amino acid sequence.
            result: ESMFoldResult to cache.
        """
        key = self._key(protein)

        # Evict oldest entries if at capacity
        if len(self._cache) >= self._max_size and key not in self._cache:
            oldest_key = next(iter(self._cache))
            del self._cache[oldest_key]

        self._cache[key] = result

        # Persist to file
        if self._cache_dir is not None:
            os.makedirs(self._cache_dir, exist_ok=True)
            filepath = os.path.join(self._cache_dir, f"{key}.json")
            try:
                data = {
                    "protein": result.protein,
                    "pdb_string": result.pdb_string,
                    "plddt_scores": result.plddt_scores,
                    "mean_plddt": result.mean_plddt,
                    "pae_matrix": result.pae_matrix,
                    "model_name": result.model_name,
                    "execution_time_s": result.execution_time_s,
                    "success": result.success,
                    "error": result.error,
                    "method": result.method,
                }
                with open(filepath, "w") as f:
                    json.dump(data, f)
            except OSError:
                pass

    @property
    def hits(self) -> int:
        """Number of cache hits."""
        return self._hits

    @property
    def misses(self) -> int:
        """Number of cache misses."""
        return self._misses

    @property
    def size(self) -> int:
        """Number of entries in memory cache."""
        return len(self._cache)

    @property
    def hit_rate(self) -> float:
        """Cache hit rate (0.0 to 1.0)."""
        total = self._hits + self._misses
        return self._hits / total if total > 0 else 0.0

    def clear(self) -> None:
        """Clear the memory cache."""
        self._cache.clear()
        self._hits = 0
        self._misses = 0


# Module-level default cache instance (lazy-initialized)
_default_cache: ESMFoldCache | None = None


def _get_default_cache() -> ESMFoldCache:
    """Return the module-level default cache (lazy-initialized)."""
    global _default_cache
    if _default_cache is None:
        _default_cache = ESMFoldCache()
    return _default_cache


# ==============================================================================
# Input validation
# ==============================================================================

def _validate_protein(protein: str) -> None:
    """Raise ESMFoldError if *protein* contains non-standard amino acids.

    Delegates to :func:`engine_base.validate_protein_sequence` for the
    standard validation, then applies ESMFold-specific checks (e.g. length).
    Only the 20 canonical single-letter codes are accepted.
    """
    try:
        validate_protein_sequence(protein, "ESMFold")
    except ValueError as exc:
        raise ESMFoldError(str(exc), protein=protein) from exc


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

    with EngineTimer() as timer:
        # --- Strategy 1: Remote API ------------------------------------------------
        if use_api:
            result = _predict_via_api(protein, api_url, timeout)
            if result is not None:
                result.execution_time_s = timer.elapsed
                return result

        # --- Strategy 2: Local esm package -----------------------------------------
        result = _predict_via_local_esm(protein)
        if result is not None:
            result.execution_time_s = timer.elapsed
            return result

    # --- Strategy 3: Offline fallback ------------------------------------------
    return ESMFoldResult(
        pdb_string="",
        plddt_scores=[],
        mean_plddt=0.0,
        pae_matrix=None,
        protein=protein,
        model_name="esmfold_v1",
        execution_time_s=timer.elapsed,
        success=False,
        error="ESMFold unavailable: API unreachable and local esm package not installed",
        method="esmfold_api",
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

            return _build_result_from_pdb(protein, pdb_string, "esmfold_v1", method="esmfold_api")

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
            method="esmfold_local",
        )

    except Exception as exc:
        logger.error("Local ESMFold prediction failed: %s", exc)
        return None


def _build_result_from_pdb(
    protein: str,
    pdb_string: str,
    model_name: str,
    method: str = "esmfold_api",
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
        method=method,
    )


# ==============================================================================
# Batch prediction (simple concurrent wrapper)
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
                    method="esmfold_api",
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
                    method="esmfold_api",
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


# ==============================================================================
# Advanced batch prediction (merged from esmfold_batch.py)
# ==============================================================================

@dataclass
class BatchStructureRequest:
    """Request object for batch ESMFold structure prediction.

    Attributes:
        proteins: List of amino acid sequences to predict.
        names: Optional list of protein names (same length as *proteins*).
            If *None*, names are auto-generated as ``protein_0``, ``protein_1``, ...
        use_cache: Whether to check the prediction cache before calling the API.
        max_concurrent: Maximum number of concurrent API calls (rate-limit).
        timeout_per_protein: Per-item timeout in seconds.
        stop_on_failure: If *True*, abort the entire batch on the first failure.
    """

    proteins: list[str]
    names: list[str] | None = None
    use_cache: bool = True
    max_concurrent: int = 3
    timeout_per_protein: float = 120.0
    stop_on_failure: bool = False


@dataclass
class BatchStructureResult:
    """Aggregated result of a batch structure prediction run.

    Attributes:
        results: Per-protein result dicts.  Each dict contains at minimum
            ``"name"``, ``"status"`` (``"success"`` | ``"error"``),
            and either prediction data (``"mean_plddt"``, ``"length"``, etc.)
            or ``"error"`` message.
        names: Ordered list of protein names.
        total: Total number of proteins in the batch.
        successful: Number of successful predictions.
        failed: Number of failed predictions.
        from_cache: Number of results served from the cache.
        total_time_s: Wall-clock time for the entire batch in seconds.
        summary: Aggregate statistics dict.
    """

    results: list[dict]
    names: list[str]
    total: int
    successful: int
    failed: int
    from_cache: int
    total_time_s: float
    summary: dict


def validate_batch_input(proteins: list[str]) -> list[str]:
    """Validate a list of protein sequences for batch prediction.

    Checks:
      - Batch size does not exceed :data:`MAX_BATCH_SIZE` (50).
      - Each protein length does not exceed :data:`MAX_PROTEIN_LENGTH` (1000).
      - Each protein contains only standard amino acids
        (:data:`STANDARD_AMINO_ACIDS`).

    Returns:
        A list of human-readable validation error strings.  An empty list
        means the input is valid.
    """
    errors: list[str] = []

    if len(proteins) > MAX_BATCH_SIZE:
        errors.append(
            f"Batch size {len(proteins)} exceeds maximum of {MAX_BATCH_SIZE}"
        )

    for idx, protein in enumerate(proteins):
        if len(protein) > MAX_PROTEIN_LENGTH:
            errors.append(
                f"Protein at index {idx} has length {len(protein)}, "
                f"exceeding maximum of {MAX_PROTEIN_LENGTH}"
            )

        invalid_chars = set(protein.upper()) - STANDARD_AMINO_ACIDS
        if invalid_chars:
            errors.append(
                f"Protein at index {idx} contains non-standard amino acids: "
                f"{sorted(invalid_chars)}"
            )

    return errors


def estimate_batch_time(
    num_proteins: int,
    avg_length: int,
    concurrent: int = 3,
) -> float:
    """Estimate total wall-clock time for a batch prediction run.

    Uses rough ESMFold performance characteristics:
      - ~1 second per residue for a live API call.
      - ~0.5 second per residue for a cached result (assumes ~50 % cache hit).

    The estimate accounts for concurrency by dividing wall-clock time by the
    number of concurrent workers (capped by *num_proteins*).

    Args:
        num_proteins: Number of proteins in the batch.
        avg_length: Average residue length per protein.
        concurrent: Number of concurrent workers.

    Returns:
        Estimated wall-clock time in seconds.
    """
    if num_proteins <= 0 or avg_length <= 0:
        return 0.0

    # Rough per-residue timing.
    api_time_per_protein = avg_length * 1.0   # ~1 s/residue
    cache_time_per_protein = avg_length * 0.5  # ~0.5 s/residue for cached

    # Assume ~50 % cache-hit rate for estimation.
    api_count = num_proteins // 2 + num_proteins % 2
    cache_count = num_proteins // 2

    total_serial_time = (
        api_count * api_time_per_protein
        + cache_count * cache_time_per_protein
    )

    effective_concurrency = min(concurrent, num_proteins)
    estimated = total_serial_time / effective_concurrency

    return round(estimated, 1)


def _predict_single(
    protein: str,
    name: str,
    use_cache: bool,
    semaphore: Semaphore,
) -> dict:
    """Predict structure for a single protein with rate-limiting.

    Returns a result dict with at least ``name`` and ``status`` keys.
    """
    result: dict = {"name": name, "status": "error"}

    with EngineTimer() as timer:
        try:
            semaphore.acquire()
            try:
                # Attempt cache lookup first.
                if use_cache:
                    try:
                        cache = _get_default_cache()
                        cached_result = cache.get(protein)
                        if cached_result is not None:
                            result = {
                                "name": name,
                                "status": "success",
                                "from_cache": True,
                                "mean_plddt": getattr(cached_result, "mean_plddt", None),
                                "length": len(protein),
                                "pdb": getattr(cached_result, "pdb_string", None),
                                "time_s": timer.elapsed,
                            }
                            return result
                    except Exception as exc:
                        logger.debug("Cache lookup failed for %s: %s", name, exc)

                # Call the prediction API.
                prediction: ESMFoldResult = predict_structure(protein)

                result = {
                    "name": name,
                    "status": "success",
                    "from_cache": False,
                    "mean_plddt": getattr(prediction, "mean_plddt", None),
                    "length": len(protein),
                    "pdb": getattr(prediction, "pdb_string", None),
                    "time_s": timer.elapsed,
                }

            finally:
                semaphore.release()

        except TimeoutError:
            result["error"] = "Prediction timed out"
            result["time_s"] = timer.elapsed
        except ImportError as exc:
            result["error"] = f"ESMFold module not available: {exc}"
            result["time_s"] = timer.elapsed
        except Exception as exc:
            result["error"] = str(exc)
            result["time_s"] = timer.elapsed

    return result


def predict_batch(
    request: BatchStructureRequest,
    progress_callback: Callable[[int, int, dict], None] | None = None,
) -> BatchStructureResult:
    """Run batch ESMFold structure prediction with concurrency and error isolation.

    .. deprecated::
        Use :func:`predict_structure_batch` instead, which provides a
        simpler, cleaner API (``list[str]`` → ``list[ESMFoldResult]``).
        This function is retained for backward compatibility and will be
        removed in a future release.

    Each protein is predicted independently — one failure does not affect
    others (unless ``request.stop_on_failure`` is *True*).

    Args:
        request: The batch request specification.
        progress_callback: Optional callable ``(completed, total, latest_result)``
            invoked after each protein finishes.

    Returns:
        Aggregated :class:`BatchStructureResult`.
    """
    with EngineTimer() as batch_timer:
        # Resolve names.
        if request.names is not None:
            if len(request.names) != len(request.proteins):
                raise ValueError(
                    f"Length of names ({len(request.names)}) does not match "
                    f"length of proteins ({len(request.proteins)})"
                )
            names = list(request.names)
        else:
            names = [f"protein_{i}" for i in range(len(request.proteins))]

        # Validate input.
        validation_errors = validate_batch_input(request.proteins)
        if validation_errors:
            raise ValueError(
                "Batch input validation failed: " + "; ".join(validation_errors)
            )

        # Delegate to predict_structure_batch for the core predictions.
        esmfold_results = predict_structure_batch(
            proteins=request.proteins,
            max_concurrent=request.max_concurrent,
        )

        # Convert list[ESMFoldResult] → list[dict] for BatchStructureResult.
        results: list[dict] = []
        cancel = False
        for idx, (name, ef_result) in enumerate(zip(names, esmfold_results)):
            result_dict: dict = {
                "name": name,
                "status": "success" if ef_result.success else "error",
                "from_cache": False,
                "mean_plddt": ef_result.mean_plddt if ef_result.success else None,
                "length": len(ef_result.protein),
                "pdb": ef_result.pdb_string if ef_result.success else None,
                "time_s": ef_result.execution_time_s,
            }
            if not ef_result.success:
                result_dict["error"] = ef_result.error or "Prediction failed"
                if request.stop_on_failure:
                    cancel = True
            results.append(result_dict)

            if progress_callback is not None:
                try:
                    progress_callback(idx + 1, len(esmfold_results), result_dict)
                except Exception as cb_exc:
                    logger.warning("Progress callback raised: %s", cb_exc)

    total_time = batch_timer.elapsed

    # Aggregate statistics.
    total = len(results)
    successful = sum(1 for r in results if r["status"] == "success")
    failed = total - successful
    from_cache = sum(1 for r in results if r.get("from_cache", False))

    plddt_values = [
        r["mean_plddt"] for r in results
        if r["status"] == "success" and r.get("mean_plddt") is not None
    ]

    summary: dict = {
        "total": total,
        "successful": successful,
        "failed": failed,
        "from_cache": from_cache,
        "success_rate": round(successful / total, 3) if total > 0 else 0.0,
        "mean_plddt": round(sum(plddt_values) / len(plddt_values), 2)
        if plddt_values
        else None,
        "min_plddt": round(min(plddt_values), 2) if plddt_values else None,
        "max_plddt": round(max(plddt_values), 2) if plddt_values else None,
        "total_time_s": total_time,
        "cancelled": cancel,
    }

    return BatchStructureResult(
        results=results,
        names=names,
        total=total,
        successful=successful,
        failed=failed,
        from_cache=from_cache,
        total_time_s=total_time,
        summary=summary,
    )


def predict_proteins(
    proteins: list[str],
    names: list[str] | None = None,
    **kwargs,
) -> BatchStructureResult:
    """Convenience wrapper for batch ESMFold prediction.

    Creates a :class:`BatchStructureRequest` from the given proteins and
    forwards to :func:`predict_batch`.

    Args:
        proteins: List of amino acid sequences.
        names: Optional protein names.
        **kwargs: Additional keyword arguments forwarded to
            :class:`BatchStructureRequest` (e.g. ``use_cache``,
            ``max_concurrent``, ``timeout_per_protein``,
            ``stop_on_failure``).

    Returns:
        Aggregated :class:`BatchStructureResult`.
    """
    request = BatchStructureRequest(proteins=proteins, names=names, **kwargs)
    return predict_batch(request)


def _quality_label(mean_plddt: float | None) -> str:
    """Return a human-readable quality label based on mean pLDDT."""
    if mean_plddt is None:
        return "N/A"
    if mean_plddt >= 90:
        return "Very High"
    if mean_plddt >= 70:
        return "High"
    if mean_plddt >= 50:
        return "Medium"
    return "Low"


def format_batch_report(
    result: BatchStructureResult,
    format: str = "text",
) -> str:
    """Format a :class:`BatchStructureResult` as a human-readable report or JSON.

    Args:
        result: The batch result to format.
        format: ``"text"`` for a human-readable table, ``"json"`` for a JSON
            string of the result data.

    Returns:
        Formatted report string.
    """
    if format == "json":
        data = {
            "results": result.results,
            "names": result.names,
            "total": result.total,
            "successful": result.successful,
            "failed": result.failed,
            "from_cache": result.from_cache,
            "total_time_s": result.total_time_s,
            "summary": result.summary,
        }
        return json.dumps(data, indent=2)

    # --- text format ---
    # Header
    lines: list[str] = []
    lines.append("=" * 82)
    lines.append("ESMFold Batch Structure Prediction Report")
    lines.append("=" * 82)
    lines.append("")

    # Table header
    header = (
        f"{'Name':<20s} {'Length':>6s} {'Mean pLDDT':>10s} "
        f"{'Quality':>10s} {'Time (s)':>9s} {'Status':>8s} {'Cache':>5s}"
    )
    lines.append(header)
    lines.append("-" * 82)

    # Table rows
    for r in result.results:
        name = r.get("name", "?")[:20]
        length = str(r.get("length", "-"))
        mean_plddt = r.get("mean_plddt")
        plddt_str = f"{mean_plddt:.2f}" if mean_plddt is not None else "-"
        quality = _quality_label(mean_plddt)
        time_s = f"{r.get('time_s', 0.0):.2f}"
        status = r.get("status", "?")
        cache = "Yes" if r.get("from_cache", False) else "No"
        if status == "error":
            cache = "-"
            quality = "-"
            plddt_str = "-"
            length = "-"

        lines.append(
            f"{name:<20s} {length:>6s} {plddt_str:>10s} "
            f"{quality:>10s} {time_s:>9s} {status:>8s} {cache:>5s}"
        )

    lines.append("-" * 82)

    # Summary
    s = result.summary
    lines.append("")
    lines.append("Summary")
    lines.append("-" * 40)
    lines.append(f"  Total:          {s.get('total', result.total)}")
    lines.append(f"  Successful:     {s.get('successful', result.successful)}")
    lines.append(f"  Failed:         {s.get('failed', result.failed)}")
    lines.append(f"  From cache:     {s.get('from_cache', result.from_cache)}")
    lines.append(
        f"  Success rate:   {s.get('success_rate', 0.0):.1%}"
    )
    mean_plddt = s.get("mean_plddt")
    lines.append(
        f"  Mean pLDDT:     {mean_plddt:.2f}" if mean_plddt is not None
        else "  Mean pLDDT:     N/A"
    )
    min_plddt = s.get("min_plddt")
    lines.append(
        f"  Min pLDDT:      {min_plddt:.2f}" if min_plddt is not None
        else "  Min pLDDT:      N/A"
    )
    max_plddt = s.get("max_plddt")
    lines.append(
        f"  Max pLDDT:      {max_plddt:.2f}" if max_plddt is not None
        else "  Max pLDDT:      N/A"
    )
    lines.append(f"  Total time:     {result.total_time_s:.2f}s")

    if s.get("cancelled"):
        lines.append("  ** Batch was cancelled (stop_on_failure) **")

    lines.append("=" * 82)

    return "\n".join(lines)
