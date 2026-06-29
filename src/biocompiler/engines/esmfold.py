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

Accuracy and Confidence
----------------------
**ESMFold pLDDT scores** (from real ESMFold model):
  - Per-residue pLDDT has ~0.8 Pearson correlation with experimental
    structure accuracy (Lin et al., Science 2023).
  - pLDDT >= 90: Very high confidence (backbone accurate to ~1 Å)
  - pLDDT 70-90: Confident (generally correct backbone)
  - pLDDT 50-70: Low confidence (may have domain-level errors)
  - pLDDT < 50: Very low (likely disordered or mispredicted)

**Our implementation** depends on API availability:
  - **API mode** (esmfold_api): Uses the ESM Atlas remote API, which
    runs the real ESMFold v1 model.  pLDDT scores have full accuracy.
  - **Local esm mode** (esmfold_local): Uses the locally-installed
    ``esm`` Python package.  Same model accuracy as API mode.
  - **Heuristic fallback** (heuristic_fallback): When both API and local
    esm are unavailable, uses sequence-based heuristics (hydrophobicity,
    charge distribution, secondary structure propensity) to produce a
    low-confidence estimate.  Mean pLDDT is calibrated per-residue based on
    predicted secondary structure (helix: 45-55, sheet: 40-50, coil: 25-35)
    with overall mean typically 35-50, capped at 55.0.  Confidence is always
    < 0.5.  This is significantly less accurate than ESMFold and should only
    be used as a last resort to avoid UNCERTAIN verdicts.
  - **Offline mode** (complete failure): Returns success=False with no
    prediction.  Only reached if the heuristic fallback itself fails.

  **Confidence levels:**
    - API/local mode, mean pLDDT >= 70: **HIGH**
    - API/local mode, mean pLDDT 50-70: **MEDIUM**
    - API/local mode, mean pLDDT < 50: **LOW**
    - Heuristic fallback: **VERY LOW** (confidence < 0.5, per-residue pLDDT
      based on SS prediction: helix 45-55, sheet 40-50, coil 25-35)
    - Offline mode (no prediction): **NONE**

**Known limitations:**
  - API availability is not guaranteed (ESM Atlas may be down or rate-limited)
  - Local esm requires GPU for reasonable speed
  - Heuristic fallback is significantly less accurate than ESMFold;
    it uses simple sequence-based rules and cannot capture long-range
    interactions.  Results should be treated as tentative estimates only.
  - ESMFold does not predict PAE (Predicted Aligned Error) via the API;
    this field is always None in our results

References
----------
- Lin et al., Science 2023; 379:1043 (ESMFold / ESM-2)
- Jumper et al., Nature 2021; 596:583 (AlphaFold2 pLDDT methodology)
"""

from __future__ import annotations

import functools
import hashlib
import json
import logging
import math
import os
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError
from dataclasses import dataclass, field
from threading import Semaphore
from typing import Any, Callable, Optional

from biocompiler.shared.constants import DEFAULT_ENGINE_TIMEOUT
from .base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    classify_score,
    validate_protein_sequence,
)
try:
    from biocompiler.shared.exceptions import ESMFoldError
except ImportError:
    from biocompiler.shared.exceptions import BioCompilerError

    class ESMFoldError(BioCompilerError):
        """Raised when ESMFold prediction fails."""

        def __init__(self, reason: str, protein: str | None = None):
            self.reason = reason
            self.protein = protein
            msg = f"ESMFold prediction failed: {reason}"
            if protein:
                msg += f" (protein length={len(protein)})"
            super().__init__(msg)

logger = logging.getLogger(__name__)

__all__ = [
    "ESMFoldError",
    "ESMFoldResult",
    "ESMFoldCache",
    "BatchStructureRequest",
    "BatchStructureResult",
    "StructureComparison",
    "STANDARD_AMINO_ACIDS",
    "DEFAULT_API_URL",
    "DEFAULT_TIMEOUT",
    "MAX_RETRIES",
    "RETRY_BASE_DELAY",
    "MAX_BATCH_SIZE",
    "MAX_PROTEIN_LENGTH",
    "clear_cache",
    "predict_structure",
    "predict_structure_batch",
    "predict_batch",
    "predict_proteins",
    "predict_pair",
    "analyze_structure",
    "is_esmfold_available",
    "parse_pdb",
    "compute_backbone_dihedrals",
    "classify_plddt",
    "estimate_contact_map",
    "validate_batch_input",
    "estimate_batch_time",
    "format_batch_report",
    "ESMFOLD_PLDDT_CORRELATION",
    "_validate_protein",
    "_build_result_from_pdb",
    "_get_default_cache",
]

# ==============================================================================
# Constants
# ==============================================================================

STANDARD_AMINO_ACIDS: set[str] = set("ACDEFGHIKLMNPQRSTVWY")

# ────────────────────────────────────────────────────────────
# Accuracy constants
# ────────────────────────────────────────────────────────────

#: Pearson correlation between ESMFold pLDDT and experimental accuracy
#: (Lin et al., Science 2023)
ESMFOLD_PLDDT_CORRELATION: float = 0.8

#: ESMFold pLDDT classification thresholds (from AlphaFold convention)
#: pLDDT >= 90: Very high (experimental accuracy)
#: pLDDT 70-90: Confident
#: pLDDT 50-70: Low confidence
#: pLDDT < 50: Very low

# ESMFold-specific constants (timeout falls back to shared DEFAULT_ENGINE_TIMEOUT)
DEFAULT_API_URL: str = "https://api.esmatlas.com/fetchPredictedStructure"
DEFAULT_TIMEOUT: float = DEFAULT_ENGINE_TIMEOUT
MAX_RETRIES: int = 3
RETRY_BASE_DELAY: float = 2.0  # seconds, doubled each attempt

MAX_BATCH_SIZE: int = 50
MAX_PROTEIN_LENGTH: int = 1000

# Named constants for previously-magic numbers
_AVAILABILITY_TIMEOUT: float = 10.0       # seconds — HEAD request timeout for is_esmfold_available
_MIN_PDB_LENGTH: int = 50                  # minimum bytes for a valid PDB response
_ERROR_TRUNCATE: int = 200                 # max chars for truncated API error messages
_CACHE_KEY_LENGTH: int = 16                # hex chars kept from SHA-256 digest for cache keys
_DEGENERATE_BOND_THRESHOLD: float = 1e-8   # below this, bonds are considered coincident in dihedral calc
_DEFAULT_BATCH_TIMEOUT: float = 120.0      # seconds — default per-protein timeout for batch requests
_API_TIME_PER_RESIDUE: float = 1.0         # seconds — estimated API time per residue for batch estimates
_CACHE_TIME_PER_RESIDUE: float = 0.5       # seconds — estimated cached-result time per residue
_REPORT_LINE_WIDTH: int = 82               # characters — width of formatted batch report lines

# NOTE: PDB ATOM record parsing uses fixed-width column slicing
# (see parse_pdb).  A regex-based parser was previously defined here
# but was never invoked; it has been removed to avoid dead code.


# ==============================================================================
# Data classes
# ==============================================================================

@dataclass
class ESMFoldResult(BaseEngineResult):
    """Result of an ESMFold protein structure prediction.

    Inherits unified fields from :class:`BaseEngineResult`:
      - sequence:       input protein sequence (alias: ``protein``)
      - primary_score:  mean pLDDT score (aliases: ``plddt``, ``mean_plddt``)
      - classification: confidence category (alias: ``confidence_class``)
      - success, error, execution_time_s, engine_name, primary_score_label

    ESMFold-specific attributes:
        plddt_scores:     Per-residue pLDDT confidence scores (0-100).
        pae_matrix:       Predicted Aligned Error matrix (residue×residue).
                          None when the API does not return PAE.
                          (alias: ``pae``)
        pdb_string:       Predicted structure in PDB format.
        model_name:       ESMFold model identifier (e.g. "esmfold_v1").
        method:           How the prediction was obtained — ``"esmfold_api"``
                          for the remote ESM Atlas API, ``"esmfold_local"``
                          for the locally-installed esm package.

    Backward compatibility:
        The constructor accepts both unified names (``sequence``,
        ``primary_score``, ``classification``) and legacy aliases
        (``protein``, ``mean_plddt``, ``plddt``, ``confidence_class``).
        Legacy names are mapped to unified fields automatically.
    """

    # --- Override base class defaults for ESMFold ---
    sequence: str = ""
    primary_score: float = 0.0
    classification: str = ""
    success: bool = True
    error: Optional[str] = None
    execution_time_s: float = 0.0
    engine_name: str = "esmfold"
    primary_score_label: str = "pLDDT"

    # --- ESMFold-specific fields ---
    plddt_scores: list[float] = field(default_factory=list)
    pae_matrix: list[list[float]] | None = None
    pdb_string: str = ""
    model_name: str = "esmfold_v1"
    method: str = "esmfold_api"

    def __init__(
        self,
        # Unified API fields (from BaseEngineResult)
        sequence: str = "",
        primary_score: float = 0.0,
        classification: str = "",
        success: bool = True,
        error: str | None = None,
        execution_time_s: float = 0.0,
        engine_name: str = "esmfold",
        primary_score_label: str = "pLDDT",
        # ESMFold-specific fields
        plddt_scores: list[float] | None = None,
        pae_matrix: list[list[float]] | None = None,
        pdb_string: str = "",
        model_name: str = "esmfold_v1",
        method: str = "esmfold_api",
        # Backward-compat constructor aliases
        protein: str | None = None,
        mean_plddt: float | None = None,
        plddt: float | None = None,
        confidence_class: str | None = None,
    ) -> None:
        # Resolve backward-compat aliases to unified fields
        if protein is not None and not sequence:
            sequence = protein
        if mean_plddt is not None and primary_score == 0.0:
            primary_score = mean_plddt
        if plddt is not None and primary_score == 0.0:
            primary_score = plddt
        if confidence_class is not None and not classification:
            classification = confidence_class

        # Auto-compute classification if not provided and score is available
        if not classification and primary_score > 0:
            classification = classify_plddt(primary_score)

        # Set all fields
        self.sequence = sequence
        self.primary_score = primary_score
        self.classification = classification
        self.success = success
        self.error = error
        self.execution_time_s = execution_time_s
        self.engine_name = engine_name
        self.primary_score_label = primary_score_label
        self.plddt_scores = plddt_scores if plddt_scores is not None else []
        self.pae_matrix = pae_matrix
        self.pdb_string = pdb_string
        self.model_name = model_name
        self.method = method

    # --- Property aliases for backward compatibility ---

    @property
    def protein(self) -> str:
        """Alias for :attr:`sequence` (backward compatibility)."""
        return self.sequence

    @protein.setter
    def protein(self, value: str) -> None:
        self.sequence = value

    @property
    def mean_plddt(self) -> float:
        """Alias for :attr:`primary_score` (backward compatibility)."""
        return self.primary_score

    @mean_plddt.setter
    def mean_plddt(self, value: float) -> None:
        self.primary_score = value

    @property
    def plddt(self) -> float:
        """Alias for :attr:`primary_score` (unified API)."""
        return self.primary_score

    @plddt.setter
    def plddt(self, value: float) -> None:
        self.primary_score = value

    @property
    def confidence_class(self) -> str:
        """Alias for :attr:`classification` (unified API)."""
        return self.classification

    @confidence_class.setter
    def confidence_class(self, value: str) -> None:
        self.classification = value

    @property
    def pae(self) -> list[list[float]] | None:
        """Alias for :attr:`pae_matrix`."""
        return self.pae_matrix

    @pae.setter
    def pae(self, value: list[list[float]] | None) -> None:
        self.pae_matrix = value

    @property
    def confidence_level(self) -> str:
        """Accuracy confidence level for the structure prediction.

        Returns one of:
          - ``"high"`` -- successful prediction with mean pLDDT >= 70
          - ``"medium"`` -- successful prediction with mean pLDDT 50-70
          - ``"low"`` -- successful prediction with mean pLDDT < 50
          - ``"very_low"`` -- heuristic fallback prediction (pLDDT < 50,
            method is ``"heuristic_fallback"``)
          - ``"none"`` -- failed prediction (no structure obtained)
        """
        if not self.success:
            return "none"
        if self.method == "heuristic_fallback":
            return "very_low"
        score = self.primary_score  # mean pLDDT
        if score >= 70:
            return "high"
        elif score >= 50:
            return "medium"
        else:
            return "low"

    @property
    def normalized_confidence(self) -> float:
        """Normalized confidence score in the 0–1 range.

        Converts the mean pLDDT score (0–100 scale) to a 0–1 confidence
        value.  For heuristic fallback predictions, the confidence is
        capped at 0.5 since these estimates are inherently low-confidence.

        The normalization follows the pLDDT interpretation convention:
          - pLDDT >= 90 → confidence >= 0.9  (very high)
          - pLDDT 70–90 → confidence 0.7–0.9 (confident)
          - pLDDT 50–70 → confidence 0.5–0.7 (low)
          - pLDDT < 50  → confidence < 0.5   (very low / disordered)

        For failed predictions (``success=False``), returns 0.0.

        Robustness: if ``primary_score`` is NaN or infinite, the property
        returns 0.0 instead of propagating the invalid value.

        Returns:
            Float in the range [0.0, 1.0].
        """
        if not self.success:
            return 0.0
        # Guard against NaN / infinite primary_score
        try:
            raw = self.primary_score / 100.0
            if not math.isfinite(raw):
                return 0.0
        except (TypeError, ZeroDivisionError):
            return 0.0
        # Clamp to [0, 1] for safety
        normalized = max(0.0, min(1.0, raw))
        # Cap heuristic fallback at 0.5
        if self.method == "heuristic_fallback":
            normalized = min(normalized, 0.5)
        return round(normalized, 4)


# ==============================================================================
# ESMFold Cache (merged from esmfold_cache.py)
# ==============================================================================

class ESMFoldCache:
    """Cache for ESMFold structure predictions.

    Supports in-memory LRU (Least Recently Used) caching with optional
    file-based persistence.  When the cache reaches *max_size*, the
    least-recently-used entry is evicted first.
    """

    def __init__(self, cache_dir: Optional[str] = None, max_size: int = 256):
        """Initialize the cache.

        Args:
            cache_dir: Directory for file-based cache persistence.
                       If None, uses in-memory only.
            max_size: Maximum number of entries in memory cache.
        """
        self._cache: OrderedDict[str, ESMFoldResult] = OrderedDict()
        self._cache_dir = cache_dir
        self._max_size = max_size
        self._hits = 0
        self._misses = 0

    @staticmethod
    def _key(protein: str) -> str:
        """Generate a cache key from a protein sequence."""
        return hashlib.sha256(protein.encode()).hexdigest()[:_CACHE_KEY_LENGTH]

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
            # Move to end (most recently used) for LRU ordering
            self._cache.move_to_end(key)
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
                    # Promote to memory cache (at end = most recently used)
                    self._cache[key] = result
                    self._hits += 1
                    return result
                except (json.JSONDecodeError, KeyError) as exc:
                    logger.warning(
                        "Corrupt cache file %s, skipping: %s", filepath, exc
                    )

        self._misses += 1
        return None

    def put(self, protein: str, result: ESMFoldResult) -> None:
        """Store a prediction in the cache.

        Args:
            protein: Amino acid sequence.
            result: ESMFoldResult to cache.
        """
        key = self._key(protein)

        # If key already exists, remove it first (will be re-added at end)
        if key in self._cache:
            del self._cache[key]

        # Evict least-recently-used entries if at capacity
        while len(self._cache) >= self._max_size:
            # popitem(last=False) removes the first (oldest / least recently used) item
            self._cache.popitem(last=False)

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
            except OSError as exc:
                logger.warning(
                    "Failed to write cache file %s: %s", filepath, exc
                )

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


def clear_cache() -> None:
    """Clear the ESMFold prediction cache.

    Clears both the ESMFoldCache (in-memory LRU + file persistence) and
    the ``functools.lru_cache`` memoization layer.
    """
    _get_default_cache().clear()
    _predict_structure_cached.cache_clear()


@functools.lru_cache(maxsize=256)
def _predict_structure_cached(protein: str) -> ESMFoldResult:
    """LRU-cached wrapper for structure prediction.

    This provides a fast memoization layer on top of the full prediction
    pipeline.  Calls are keyed on the protein sequence string only;
    parameters like ``use_api`` and ``timeout`` are NOT part of the key
    — the assumption is that the same sequence should always yield the
    same structure prediction.

    Call :func:`clear_cache` to invalidate.
    """
    return _predict_structure_uncached(protein)


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
        with urllib.request.urlopen(req, timeout=_AVAILABILITY_TIMEOUT) as resp:
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
        logger.debug("Local esm package not importable")

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
    if b1_norm < _DEGENERATE_BOND_THRESHOLD:
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

# pLDDT classification thresholds for use with classify_score()
# NOTE: classify_score uses >= comparison, but the original classify_plddt
# uses strict >.  We add a tiny epsilon to thresholds to emulate strict >
# behaviour for all practical pLDDT values (0-100, typically 2 decimal places).
_EPS = 1e-9

_PLDDT_THRESHOLDS: list[tuple[float, str]] = [
    (90 + _EPS, "Very high (experimental)"),
    (70 + _EPS, "Confident"),
    (50 + _EPS, "Low confidence"),
]


def classify_plddt(mean_plddt: float) -> str:
    """Classify a mean pLDDT score into a confidence category.

    Thresholds follow the AlphaFold / ESMFold convention:

        >= 90 : "Very high (experimental)"
        70-90: "Confident"
        50-70: "Low confidence"
        < 50 : "Very low"

    Delegates to :func:`engine_base.classify_score` for unified
    classification logic.

    Args:
        mean_plddt: Average per-residue pLDDT score (0-100).

    Returns:
        Human-readable confidence classification string.
    """
    return classify_score(
        mean_plddt,
        thresholds=_PLDDT_THRESHOLDS,
        fallback="Very low",
    )


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
    organism: str = "Homo_sapiens",
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> ESMFoldResult:
    """Predict the 3-D structure of a protein using ESMFold.

    This function is the main public API for structure prediction.  It
    automatically checks the LRU cache before performing a prediction;
    repeated calls with the same sequence return the cached result.

    Strategy priority:
      1. **LRU cache** — Return cached result if this sequence was predicted
         before (via ``functools.lru_cache``).
      2. **ESMFoldCache** — Check the module-level ESMFoldCache instance
         (in-memory LRU + optional file persistence).
      3. **API** — POST the sequence to ESM Atlas and parse the PDB response.
      4. **Local esm** — ``import esm`` and run ESMFold locally (if installed).
      5. **Heuristic fallback** — Use sequence-based heuristics (hydrophobicity,
         charge distribution, secondary structure propensity) to produce a
         low-confidence estimate.  Per-residue pLDDT is calibrated based on
         predicted secondary structure (helix: 45-55, sheet: 40-50, coil: 25-35)
         with mean typically 35-50, capped at 55.0.  Confidence is always < 0.5.
         The ``method`` field is set to ``"heuristic_fallback"``.
         This is NOT a substitute for ESMFold — it exists so that structure
         predicates can return a tentative verdict instead of UNCERTAIN.
      6. **Complete failure** — Return ``ESMFoldResult(success=False)`` with
         an error (only if the heuristic fallback itself raises an exception).

    Retry logic: up to 3 attempts with exponential backoff on transient
    API errors (network, 5xx, 429 rate-limit).  Invalid protein input
    (non-standard AA) causes immediate failure with no retry.

    Args:
        protein: Protein sequence (single-letter amino acid codes).
        organism: Target organism name (for API consistency across engines;
            currently unused by ESMFold but reserved for future use).
        use_api: If True, try the remote API first.
        api_url: ESM Atlas API endpoint.
        timeout: HTTP request timeout in seconds.

    Returns:
        ESMFoldResult with PDB string, pLDDT scores, and metadata.

    Raises:
        ESMFoldError: If the protein contains invalid amino acid codes.
    """
    # Check ESMFoldCache first (in-memory LRU + file persistence)
    cache = _get_default_cache()
    cached = cache.get(protein)
    if cached is not None:
        return cached

    # Perform the actual prediction
    result = _predict_structure_uncached(
        protein, organism=organism, use_api=use_api,
        api_url=api_url, timeout=timeout,
    )

    # Store in ESMFoldCache for future lookups
    cache.put(protein, result)

    return result


def _predict_structure_uncached(
    protein: str,
    organism: str = "Homo_sapiens",
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> ESMFoldResult:
    """Internal uncached implementation of predict_structure.

    See :func:`predict_structure` for documentation.  This function
    performs the actual prediction without checking caches.
    """
    _validate_protein(protein)

    with EngineTimer() as timer:
        # --- Strategy 1: Remote API ------------------------------------------------
        if use_api:
            result = _predict_via_api(protein, api_url, timeout)
            if result is not None:
                result.execution_time_s = round(timer.elapsed, 4)
                return result

        # --- Strategy 2: Local esm package -----------------------------------------
        result = _predict_via_local_esm(protein)
        if result is not None:
            result.execution_time_s = round(timer.elapsed, 4)
            return result

    # --- Strategy 3: Heuristic fallback (offline) ------------------------------
    try:
        from .esmfold_fallback import predict_structure_heuristic
        heuristic = predict_structure_heuristic(protein)
        result = ESMFoldResult(
            pdb_string="",
            plddt_scores=heuristic["plddt_scores"],
            mean_plddt=heuristic["mean_plddt"],
            pae_matrix=None,
            protein=protein,
            model_name=heuristic["model_name"],
            execution_time_s=round(timer.elapsed, 4),
            success=True,
            error=None,
            method=heuristic["method"],
        )
        logger.info(
            "ESMFold API/local unavailable; used heuristic fallback "
            "(estimated pLDDT=%.1f, confidence=%.2f)",
            heuristic["mean_plddt"],
            heuristic["confidence"],
        )
        return result
    except Exception as exc:
        logger.warning("Heuristic fallback also failed: %s", exc)

    # --- Strategy 4: Complete failure -----------------------------------------
    return ESMFoldResult(
        pdb_string="",
        plddt_scores=[],
        mean_plddt=0.0,
        pae_matrix=None,
        protein=protein,
        model_name="esmfold_v1",
        execution_time_s=round(timer.elapsed, 4),
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

    # SSRF protection: only allow ESM Atlas domain
    from urllib.parse import urlparse
    parsed = urlparse(api_url)
    if parsed.scheme != "https" or parsed.netloc not in (
        "api.esmatlas.com",
        "www.api.esmatlas.com",
    ):
        raise ValueError(
            f"SSRF protection: api_url must point to api.esmatlas.com, got {parsed.netloc!r}"
        )

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

            if not pdb_string or len(pdb_string) < _MIN_PDB_LENGTH:
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
                last_error = f"API client error ({exc.code}): {exc.read().decode('utf-8', errors='replace')[:_ERROR_TRUNCATE]}"
                logger.error("API client error, not retrying: %s", last_error)
                break
            else:
                last_error = f"HTTP error {exc.code}"
                logger.warning("Attempt %d/%d: %s", attempt, MAX_RETRIES, last_error)
                continue

        except urllib.error.URLError as exc:
            wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
            last_error = f"Network error: {exc.reason}"
            logger.warning("Attempt %d/%d: %s (retry in %.1fs)", attempt, MAX_RETRIES, last_error, wait)
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
# Pair prediction (wild-type vs mutant / optimized)
# ==============================================================================

@dataclass
class StructureComparison:
    """Result of comparing two protein structure predictions.

    Typically used for wild-type vs. mutant (optimized) comparisons.

    Attributes:
        wildtype_result:    ESMFoldResult for the wild-type sequence.
        mutant_result:      ESMFoldResult for the mutant / optimized sequence.
        plddt_delta:        Difference in mean pLDDT (mutant − wildtype).
                            Positive values indicate the mutant has higher
                            predicted confidence.
        plddt_per_residue_delta: Per-residue pLDDT difference
                                 (mutant − wildtype).  Length matches the
                                 shorter of the two sequences; beyond that,
                                 values are None.
        confidence_delta:   Difference in normalized_confidence
                            (mutant − wildtype).
        improved:           True if the mutant has equal or higher mean pLDDT
                            than the wildtype.
    """

    wildtype_result: ESMFoldResult
    mutant_result: ESMFoldResult
    plddt_delta: float
    plddt_per_residue_delta: list[float | None]
    confidence_delta: float
    improved: bool


def predict_pair(
    wildtype: str,
    mutant: str,
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> StructureComparison:
    """Predict structures for wild-type and mutant sequences in the same batch.

    Convenience method for comparing wild-type vs. optimized structures.
    Both sequences are predicted in the same batch call for efficiency,
    then compared on pLDDT and normalized confidence.

    This is particularly useful for evaluating whether a sequence
    optimization (e.g., codon optimization, stabilizing mutations)
    is predicted to improve structural confidence.

    Args:
        wildtype: Wild-type protein sequence (single-letter codes).
        mutant:   Mutant / optimized protein sequence (single-letter codes).
        use_api:  If True, try the remote API first.
        api_url:  ESM Atlas API endpoint.
        timeout:  Per-protein request timeout in seconds.

    Returns:
        :class:`StructureComparison` with both results and delta metrics.

    Raises:
        ESMFoldError: If either sequence contains invalid amino acid codes.
    """
    # Predict both in one batch call for efficiency
    batch = predict_structure_batch(
        [wildtype, mutant],
        max_concurrent=2,
        use_api=use_api,
        api_url=api_url,
        timeout=timeout,
    )

    # Handle the case where batch results are incomplete or predictions failed.
    # Each prediction always returns an ESMFoldResult (even on failure), but
    # we guard against edge cases where the batch might return fewer entries.
    if len(batch.results) < 2:
        # If we got zero or one result, create a failure placeholder for
        # whichever prediction is missing so the comparison can still be
        # constructed (callers can check .success on each result).
        wt_result = batch.results[0] if len(batch.results) >= 1 else ESMFoldResult(
            protein=wildtype, success=False,
            error="Prediction unavailable: batch returned insufficient results",
        )
        mt_result = batch.results[1] if len(batch.results) >= 2 else ESMFoldResult(
            protein=mutant, success=False,
            error="Prediction unavailable: batch returned insufficient results",
        )
    else:
        wt_result = batch.results[0]
        mt_result = batch.results[1]

    # Compute delta metrics — use 0.0 for failed predictions so that
    # the delta reflects the asymmetric outcome (e.g. if only the mutant
    # failed, plddt_delta will be negative).
    wt_score = wt_result.primary_score if wt_result.success else 0.0
    mt_score = mt_result.primary_score if mt_result.success else 0.0
    wt_conf = wt_result.normalized_confidence if wt_result.success else 0.0
    mt_conf = mt_result.normalized_confidence if mt_result.success else 0.0

    plddt_delta = mt_score - wt_score
    confidence_delta = mt_conf - wt_conf

    # Per-residue delta — only meaningful when both predictions succeeded
    # and have per-residue scores.
    if wt_result.success and mt_result.success:
        wt_scores = wt_result.plddt_scores
        mt_scores = mt_result.plddt_scores
        min_len = min(len(wt_scores), len(mt_scores))
        max_len = max(len(wt_scores), len(mt_scores))
        per_residue_delta: list[float | None] = []
        for i in range(max_len):
            if i < min_len:
                per_residue_delta.append(round(mt_scores[i] - wt_scores[i], 2))
            else:
                per_residue_delta.append(None)
    else:
        # Cannot compute meaningful per-residue deltas when a prediction failed
        per_residue_delta = []

    # A comparison is "improved" only when both predictions succeeded and
    # the mutant has equal or higher pLDDT.  If either prediction failed,
    # we conservatively report not improved.
    improved = wt_result.success and mt_result.success and plddt_delta >= 0

    return StructureComparison(
        wildtype_result=wt_result,
        mutant_result=mt_result,
        plddt_delta=round(plddt_delta, 2),
        plddt_per_residue_delta=per_residue_delta,
        confidence_delta=round(confidence_delta, 4),
        improved=improved,
    )


# ==============================================================================
# Batch prediction (simple concurrent wrapper)
# ==============================================================================

def predict_structure_batch(
    sequences: list[str],
    max_concurrent: int = 3,
    max_workers: int | None = None,
    use_api: bool = True,
    api_url: str = DEFAULT_API_URL,
    timeout: float = DEFAULT_TIMEOUT,
) -> BatchResult[ESMFoldResult]:
    """Predict structures for multiple proteins concurrently.

    Uses a :class:`ThreadPoolExecutor` with a semaphore for rate
    limiting.  Individual prediction failures are captured in the
    returned ``ESMFoldResult`` objects — they never crash sibling
    predictions.

    Args:
        sequences:      List of protein sequences (single-letter codes).
        max_concurrent: Maximum number of concurrent API requests.
        max_workers:    Maximum number of worker threads.  If *None*,
            defaults to *max_concurrent*.
        use_api:       If True, try the remote API first.
        api_url:       ESM Atlas API endpoint.
        timeout:       Per-protein request timeout in seconds.

    Returns:
        :class:`BatchResult` [ESMFoldResult] containing per-protein
        results in the same order as *sequences*, along with aggregate
        timing and success/failure counts.
    """
    workers = max_workers if max_workers is not None else max_concurrent
    semaphore = Semaphore(max_concurrent)
    results: dict[int, ESMFoldResult] = {}

    with EngineTimer() as batch_timer:
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

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_predict_with_semaphore, i, prot): i
                for i, prot in enumerate(sequences)
            }
            for future in as_completed(futures):
                idx, result = future.result()
                results[idx] = result

    # Preserve input order and wrap in BatchResult
    result_list = [results[i] for i in range(len(sequences))]
    return BatchResult(
        results=result_list,
        total_time_s=round(batch_timer.elapsed, 4),
    )


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
    timeout_per_protein: float = _DEFAULT_BATCH_TIMEOUT
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

    results: list[dict[str, Any]]
    names: list[str]
    total: int
    successful: int
    failed: int
    from_cache: int
    total_time_s: float
    summary: dict[str, Any]


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
    api_time_per_protein = avg_length * _API_TIME_PER_RESIDUE
    cache_time_per_protein = avg_length * _CACHE_TIME_PER_RESIDUE

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
) -> dict[str, Any]:
    """Predict structure for a single protein with rate-limiting.

    Returns a result dict with at least ``name`` and ``status`` keys.
    """
    result: dict[str, Any] = {"name": name, "status": "error"}

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
                                "time_s": round(timer.elapsed, 4),
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
                    "time_s": round(timer.elapsed, 4),
                }

            finally:
                semaphore.release()

        except TimeoutError:
            result["error"] = "Prediction timed out"
            result["time_s"] = round(timer.elapsed, 4)
        except ImportError as exc:
            result["error"] = f"ESMFold module not available: {exc}"
            result["time_s"] = round(timer.elapsed, 4)
        except Exception as exc:
            result["error"] = str(exc)
            result["time_s"] = round(timer.elapsed, 4)

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
    import warnings as _warnings
    _warnings.warn(
        "esmfold.predict_batch() is deprecated — use "
        "predict_structure_batch() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
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
        batch_result = predict_structure_batch(
            sequences=request.proteins,
            max_concurrent=request.max_concurrent,
        )

        # Convert BatchResult[ESMFoldResult] → list[dict[str, Any]] for BatchStructureResult.
        results: list[dict[str, Any]] = []
        cancel = False
        for idx, (name, ef_result) in enumerate(zip(names, batch_result.results)):
            result_dict: dict[str, Any] = {
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
                    progress_callback(idx + 1, len(names), result_dict)
                except Exception as cb_exc:
                    logger.warning("Progress callback raised: %s", cb_exc)

    total_time = round(batch_timer.elapsed, 4)

    # Aggregate statistics.
    total = len(results)
    successful = sum(1 for r in results if r["status"] == "success")
    failed = total - successful
    from_cache = sum(1 for r in results if r.get("from_cache", False))

    plddt_values = [
        r["mean_plddt"] for r in results
        if r["status"] == "success" and r.get("mean_plddt") is not None
    ]

    summary: dict[str, Any] = {
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
    **kwargs: Any,
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
    lines.append("=" * _REPORT_LINE_WIDTH)
    lines.append("ESMFold Batch Structure Prediction Report")
    lines.append("=" * _REPORT_LINE_WIDTH)
    lines.append("")

    # Table header
    header = (
        f"{'Name':<20s} {'Length':>6s} {'Mean pLDDT':>10s} "
        f"{'Quality':>10s} {'Time (s)':>9s} {'Status':>8s} {'Cache':>5s}"
    )
    lines.append(header)
    lines.append("-" * _REPORT_LINE_WIDTH)

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

    lines.append("-" * _REPORT_LINE_WIDTH)

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

    lines.append("=" * _REPORT_LINE_WIDTH)

    return "\n".join(lines)
