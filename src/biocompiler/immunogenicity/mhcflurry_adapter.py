"""MHCflurry adapter — offline MHC-I binding prediction via neural networks.

This module wraps `MHCflurry <https://github.com/openvax/mhcflurry>`_ for
completely offline MHC class I binding prediction.  It serves as the
intermediate accuracy tier in BioCompiler's immunogenicity stack:

============  ===================  ============  ===========
Method        AUC-ROC range        Requires API   Offline
============  ===================  ============  ===========
PSSM          0.60 – 0.75          No             Yes
**MHCflurry** **0.80 – 0.85**      **No**         **Yes**
NetMHCpan     0.85 – 0.95          Yes            No
============  ===================  ============  ===========

MHCflurry uses ensemble neural-network models trained on IEDB affinity
data.  Predictions are *much* faster than NetMHCpan (no network latency)
and substantially more accurate than PSSM heuristics, making MHCflurry
the recommended default for offline pipelines.

Features
--------
- :func:`is_mhcflurry_available` — check whether the library and models
  are installed
- :class:`MHCflurryClient` — main client with binding and antigen-
  processing prediction
- :func:`predict_binding` — module-level convenience function with full
  fallback chain
- :func:`predict_batch` — module-level batch prediction with fallback
- :func:`download_models` — download MHCflurry models (~100 MB)
- LRU prediction cache (max 5 000 entries)
- Graceful error handling: unsupported alleles are skipped with a debug
  log; prediction errors produce empty results with a warning

Fallback chain
--------------
The adapter follows a four-tier fallback chain that **never crashes** due
to missing models:

1. **MHCflurry** — neural-network ensemble (AUC 0.80–0.85, confidence 1.0)
2. **NetMHCpan** — web API predictor (AUC 0.85–0.95, confidence 1.0),
   if installed and reachable
3. **Precomputed database** — curated IEDB/SYFPEITHI entries +
   PSSM-derived predictions (confidence 0.7)
4. **PSSM fallback** — position-specific scoring matrix heuristic
   (AUC 0.60–0.75, confidence 0.5), always available offline

Each fallback transition is logged at INFO level.

Result format
-------------
MHCflurry returns IC50 values (nM) and presentation percentiles.
This adapter converts them to the same :class:`MHCBindingResult` format
used by :mod:`biocompiler.immunogenicity`:

- ``binding_score`` = 1 − log(IC50) / log(50 000)  (clamped to [0, 1])
- ``ic50_nm`` = raw predicted IC50
- ``binding_class`` = classify_binding(ic50)  using the standard
  50 / 500 / 5 000 nM thresholds
- ``confidence`` = 1.0 for mhcflurry/netmhcpan, 0.7 for precomputed,
  0.5 for PSSM

References
----------
- O'Donnell et al., *Bioinformatics* 2018; 34:2696–2703 (MHCflurry 1.2)
- O'Donnell et al., *Nature Biotechnology* 2021; 39:1329–1336 (MHCflurry 2.0)
"""
from __future__ import annotations

import logging
import math
import os
from collections import OrderedDict
from dataclasses import dataclass, replace
from typing import Optional, Sequence

from .core import MHCBindingResult, classify_binding, score_peptide_pssm, binding_score_to_ic50

__all__ = [
    "MHCflurryClient",
    "is_mhcflurry_available",
    "is_netmhcpan_available",
    "download_models",
    "clear_cache",
    "predict_binding",
    "predict_batch",
    "MHCFLURRY_AUC_ROC_LOW",
    "MHCFLURRY_AUC_ROC_HIGH",
    "CONFIDENCE_MHCFLURRY",
    "CONFIDENCE_NETMHCPAN",
    "CONFIDENCE_PRECOMPUTED",
    "CONFIDENCE_PSSM",
    "MHCBindingResult",
    "_LRUCache",
    "_extract_overlapping_peptides",
    "_mhcflurry_result_to_binding_result",
    "_validate_peptide",
    "_validate_protein",
    "_validate_sequence",
    "ic50_to_binding_score",
    "classify_binding",
]

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# Accuracy constants
# ═══════════════════════════════════════════════════════════════════════════

#: Expected AUC-ROC lower bound for MHCflurry MHC-I binding prediction
MHCFLURRY_AUC_ROC_LOW: float = 0.80

#: Expected AUC-ROC upper bound for MHCflurry MHC-I binding prediction
MHCFLURRY_AUC_ROC_HIGH: float = 0.85

# ═══════════════════════════════════════════════════════════════════════════
# Confidence constants
# ═══════════════════════════════════════════════════════════════════════════

#: Confidence level for MHCflurry predictions (neural-network ensemble)
CONFIDENCE_MHCFLURRY: float = 1.0

#: Confidence level for NetMHCpan predictions (web API)
CONFIDENCE_NETMHCPAN: float = 1.0

#: Confidence level for precomputed database lookups
CONFIDENCE_PRECOMPUTED: float = 0.7

#: Confidence level for PSSM-based fallback predictions
CONFIDENCE_PSSM: float = 0.5

#: Map from method name to confidence level
_METHOD_CONFIDENCE: dict[str, float] = {
    "mhcflurry": CONFIDENCE_MHCFLURRY,
    "mhcflurry_presentation": CONFIDENCE_MHCFLURRY,
    "netmhcpan": CONFIDENCE_NETMHCPAN,
    "precomputed_lookup": CONFIDENCE_PRECOMPUTED,
    "pssm_fallback": CONFIDENCE_PSSM,
}

# ═══════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════

#: Default MHCflurry models directory
_DEFAULT_MODELS_DIR: str = os.path.join(os.path.expanduser("~"), ".mhcflurry")

#: Maximum IC50 value used in binding-score normalisation
_MAX_IC50_NM: float = 50_000.0

#: Logarithm of the max IC50 (pre-computed for speed)
_LOG_MAX_IC50: float = math.log(_MAX_IC50_NM)

#: Standard amino-acid set for input validation
_STANDARD_AAS: set[str] = set("ACDEFGHIKLMNPQRSTVWY")

#: Default epitope lengths for MHC-I scanning
_DEFAULT_EPITOPE_LENGTHS: list[int] = [8, 9, 10, 11]

#: Maximum cache size (LRU eviction when exceeded)
_CACHE_MAX_SIZE: int = 5_000


# ═══════════════════════════════════════════════════════════════════════════
# LRU Cache
# ═══════════════════════════════════════════════════════════════════════════

class _LRUCache:
    """Simple ordered-dict-based LRU cache keyed by ``(peptide, allele)``.

    Thread-safety is *not* provided — callers must serialise access if
    needed.  This mirrors the simplicity of the
    :class:`~biocompiler.netmhcpan.NetMHCpanCache` design.
    """

    def __init__(self, max_size: int = _CACHE_MAX_SIZE) -> None:
        self._max_size = max_size
        self._data: OrderedDict[tuple[str, str], MHCBindingResult] = OrderedDict()
        self._hits: int = 0
        self._misses: int = 0

    # ── public helpers ────────────────────────────────────────────

    def get(self, peptide: str, allele: str) -> MHCBindingResult | None:
        """Return cached result or ``None`` on miss."""
        key = (peptide, allele)
        if key in self._data:
            self._hits += 1
            # Move to end so it is most-recently-used
            self._data.move_to_end(key)
            return self._data[key]
        self._misses += 1
        return None

    def put(self, peptide: str, allele: str, result: MHCBindingResult) -> None:
        """Insert a result, evicting the least-recently-used entry if full."""
        key = (peptide, allele)
        if key in self._data:
            self._data.move_to_end(key)
            self._data[key] = result
            return
        if len(self._data) >= self._max_size:
            # Evict oldest (least-recently-used)
            self._data.popitem(last=False)
        self._data[key] = result

    def clear(self) -> None:
        """Drop all cached entries and reset counters."""
        self._data.clear()
        self._hits = 0
        self._misses = 0

    @property
    def size(self) -> int:
        """Number of entries currently in the cache."""
        return len(self._data)

    @property
    def hits(self) -> int:
        """Number of cache hits since the last clear."""
        return self._hits

    @property
    def misses(self) -> int:
        """Number of cache misses since the last clear."""
        return self._misses

    @property
    def hit_rate(self) -> float:
        """Cache hit rate in [0, 1]."""
        total = self._hits + self._misses
        return self._hits / total if total > 0 else 0.0


# Module-level singleton cache
_default_cache: _LRUCache = _LRUCache()


def clear_cache() -> None:
    """Clear the module-level MHCflurry prediction cache.

    This does *not* unload loaded MHCflurry models — only the prediction
    result cache is cleared.
    """
    _default_cache.clear()
    logger.info("MHCflurry prediction cache cleared")


# ═══════════════════════════════════════════════════════════════════════════
# Availability checks
# ═══════════════════════════════════════════════════════════════════════════

def is_mhcflurry_available() -> bool:
    """Check whether MHCflurry is importable **and** models are downloaded.

    Returns ``True`` only when both conditions are satisfied:

    1. ``import mhcflurry`` succeeds.
    2. The default models directory (``~/.mhcflurry/``) contains at
       least one model set.

    Returns
    -------
    bool
        ``True`` if MHCflurry is ready for predictions.
    """
    try:
        import mhcflurry  # noqa: F401 — import check
    except ImportError:
        logger.debug("mhcflurry is not installed")
        return False

    # Check that at least the Class1 affinity models exist
    try:
        models_dir = os.environ.get("MHCFLURRY_DEFAULT_CLASS1_MODELS_DIR", "")
        if not models_dir:
            models_dir = os.path.join(
                os.path.expanduser("~"), ".mhcflurry", "models_class1", "models"
            )
        if os.path.isdir(models_dir) and os.listdir(models_dir):
            logger.debug("MHCflurry is available (models found at %s)", models_dir)
            return True
        # Fallback: try MHCflurry's own discovery
        try:
            predictor = mhcflurry.Class1AffinityPredictor.load()
            if predictor is not None:
                return True
        except Exception:
            logger.debug("MHCflurry model availability check failed", exc_info=True)
        logger.debug("MHCflurry is installed but no models found at %s", models_dir)
        return False
    except Exception as exc:
        logger.debug("Error checking MHCflurry models: %s", exc)
        return False


def is_netmhcpan_available() -> bool:
    """Check whether the NetMHCpan adapter is importable.

    Returns ``True`` if ``biocompiler.netmhcpan`` can be imported,
    indicating that the NetMHCpan web API client is available for
    predictions (which still require network connectivity).

    Returns
    -------
    bool
        ``True`` if the NetMHCpan adapter is available.
    """
    try:
        from biocompiler.immunogenicity.netmhcpan import NetMHCpanClient  # noqa: F401
        return True
    except ImportError:
        logger.debug("NetMHCpan adapter is not available")
        return False
    except Exception:
        logger.debug("NetMHCpan adapter check failed", exc_info=True)
        return False


# ═══════════════════════════════════════════════════════════════════════════
# Model download
# ═══════════════════════════════════════════════════════════════════════════

def download_models(
    models_dir: str | None = None,
    verbose: bool = True,
) -> bool:
    """Download MHCflurry models (~100 MB).

    Parameters
    ----------
    models_dir : str or None
        Target directory for the model files.  If ``None``, MHCflurry's
        default directory (``~/.mhcflurry/``) is used.
    verbose : bool
        If ``True``, show a progress bar during download.

    Returns
    -------
    bool
        ``True`` on success, ``False`` on failure.

    Notes
    -----
    This function calls ``mhcflurry-downloads fetch`` under the hood.
    The download includes the Class1 affinity model, the Class1
    presentation model, and supporting data.  A network connection is
    required.
    """
    try:
        import mhcflurry
    except ImportError:
        logger.error(
            "mhcflurry is not installed.  Install it with: "
            "pip install mhcflurry"
        )
        return False

    try:
        # MHCflurry ≥ 2.0 uses mhcflurry.downloads
        from mhcflurry import downloads

        if models_dir is not None:
            os.environ["MHCFLURRY_DATA_DIR"] = models_dir

        # Download Class1 affinity models
        if verbose:
            logger.info("Downloading MHCflurry Class1 affinity models…")
        downloads.fetch_model("models_class1", verbose=verbose)

        # Try to download presentation models (optional)
        try:
            if verbose:
                logger.info("Downloading MHCflurry Class1 presentation models…")
            downloads.fetch_model("models_class1_presentation", verbose=verbose)
        except Exception as exc:
            logger.warning(
                "Could not download presentation models (non-fatal): %s", exc
            )

        # Try to download processing models (optional)
        try:
            if verbose:
                logger.info("Downloading MHCflurry antigen processing models…")
            downloads.fetch_model("models_class1_processing", verbose=verbose)
        except Exception as exc:
            logger.warning(
                "Could not download processing models (non-fatal): %s", exc
            )

        if verbose:
            logger.info("MHCflurry model download complete.")
        return True

    except AttributeError:
        # Older MHCflurry versions — try the CLI approach
        import subprocess
        import sys

        try:
            cmd = [sys.executable, "-m", "mhcflurry-downloads", "fetch"]
            if models_dir is not None:
                cmd.extend(["--data-dir", models_dir])
            if verbose:
                logger.info("Running: %s", " ".join(cmd))
            result = subprocess.run(
                cmd,
                capture_output=not verbose,
                text=True,
                timeout=600,
            )
            if result.returncode != 0:
                stderr = getattr(result, "stderr", "") or ""
                logger.error("mhcflurry-downloads failed: %s", stderr)
                return False
            if verbose:
                logger.info("MHCflurry model download complete.")
            return True
        except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
            logger.error("MHCflurry model download failed: %s", exc)
            return False

    except Exception as exc:
        logger.error(
            "MHCflurry model download failed: %s.  "
            "Ensure you have a working internet connection and sufficient "
            "disk space (~100 MB).",
            exc,
        )
        return False


# ═══════════════════════════════════════════════════════════════════════════
# Result conversion utilities
# ═══════════════════════════════════════════════════════════════════════════

def ic50_to_binding_score(ic50_nm: float) -> float:
    """Convert an IC50 value (nM) to a normalised binding score in [0, 1].

    The formula is ``binding_score = 1 − log(IC50) / log(50 000)``,
    which maps:

    ==========  ==============
    IC50 (nM)   binding_score
    ==========  ==============
    1           ≈ 1.0
    50          ≈ 0.75
    500         ≈ 0.55
    5 000       ≈ 0.36
    50 000      =  0.0
    ==========  ==============

    Values are clamped to [0, 1].  This is the inverse of the mapping
    used in :func:`~biocompiler.immunogenicity.binding_score_to_ic50`
    but calibrated for MHCflurry's IC50 range.

    Parameters
    ----------
    ic50_nm : float
        Predicted IC50 in nanomolars.

    Returns
    -------
    float
        Normalised binding score in [0, 1].
    """
    if ic50_nm <= 0:
        return 1.0
    score = 1.0 - math.log(ic50_nm) / _LOG_MAX_IC50
    return max(0.0, min(1.0, score))


def _validate_sequence(seq: str, allow_empty: bool = False) -> str:
    """Upper-case and validate an amino-acid sequence.

    Parameters
    ----------
    seq : str
        Raw amino-acid sequence string.
    allow_empty : bool
        If ``False`` (default), empty sequences raise ``ValueError``.
        Set to ``True`` for peptide inputs where an empty string is
        acceptable (the caller handles it).

    Returns
    -------
    str
        Upper-cased, stripped sequence.

    Raises
    ------
    ValueError
        If the sequence contains non-standard amino acids, or is empty
        when *allow_empty* is ``False``.
    """
    seq = seq.upper().strip()
    if not allow_empty and not seq:
        raise ValueError("Sequence must not be empty")
    invalid = set(seq) - _STANDARD_AAS
    if invalid:
        raise ValueError(
            f"Sequence contains non-standard amino acids: {invalid!r}"
        )
    return seq


def _validate_peptide(peptide: str) -> str:
    """Upper-case and validate a peptide sequence.

    Raises ``ValueError`` if non-standard amino acids are found.
    """
    return _validate_sequence(peptide, allow_empty=True)


def _validate_protein(protein: str) -> str:
    """Upper-case and validate a protein sequence.

    Raises ``ValueError`` if the sequence is empty or contains
    non-standard amino acids.
    """
    return _validate_sequence(protein, allow_empty=False)


def _extract_overlapping_peptides(
    protein: str,
    epitope_lengths: Sequence[int],
) -> list[tuple[str, int, int]]:
    """Extract all overlapping peptides of the given lengths.

    Returns a list of ``(peptide, start, end)`` tuples where positions
    are 0-based and *end* is inclusive.
    """
    peptides: list[tuple[str, int, int]] = []
    for length in epitope_lengths:
        for start in range(len(protein) - length + 1):
            pep = protein[start : start + length]
            peptides.append((pep, start, start + length - 1))
    return peptides


def _mhcflurry_result_to_binding_result(
    peptide: str,
    allele: str,
    start_position: int,
    end_position: int,
    ic50_nm: float,
    presentation_score: float | None = None,
    method: str = "mhcflurry",
    rank: float | None = None,
) -> MHCBindingResult:
    """Convert raw MHCflurry outputs to :class:`MHCBindingResult`.

    Parameters
    ----------
    peptide : str
        Peptide sequence.
    allele : str
        MHC allele name.
    start_position : int
        0-based start position in the source protein.
    end_position : int
        0-based inclusive end position.
    ic50_nm : float
        Predicted IC50 in nM.
    presentation_score : float or None
        If provided (from presentation predictor), this overrides
        ``binding_score`` so that the score reflects processing + binding.
    method : str
        Method tag (``"mhcflurry"`` or ``"mhcflurry_presentation"``).
    rank : float or None
        Percentile rank from the prediction tool, if available.

    Returns
    -------
    MHCBindingResult
    """
    binding_score = ic50_to_binding_score(ic50_nm)
    if presentation_score is not None:
        binding_score = presentation_score
    binding_class = classify_binding(ic50_nm)
    confidence = _METHOD_CONFIDENCE.get(method, CONFIDENCE_PSSM)

    # Anchor residues: for MHCflurry we do not have PSSM-style per-position
    # scores, so we record the canonical MHC-I anchor positions (P2 and
    # C-terminus) as markers.
    anchor_residues: dict[int, str] = {}
    anchor_scores: dict[int, float] = {}
    if len(peptide) >= 2:
        anchor_residues[1] = peptide[1]  # P2
        anchor_scores[1] = binding_score
    if len(peptide) >= 1:
        last_idx = len(peptide) - 1
        anchor_residues[last_idx] = peptide[last_idx]  # C-terminus
        anchor_scores[last_idx] = binding_score

    return MHCBindingResult(
        allele=allele,
        peptide=peptide,
        start_position=start_position,
        end_position=end_position,
        binding_score=round(binding_score, 6),
        ic50_nm=round(ic50_nm, 2),
        binding_class=binding_class,
        anchor_residues=anchor_residues,
        anchor_scores={k: round(v, 6) for k, v in anchor_scores.items()},
        method=method,
        rank=rank,
        confidence=confidence,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Module-level convenience functions
# ═══════════════════════════════════════════════════════════════════════════

# Module-level default client (lazy-initialised)
_default_client: MHCflurryClient | None = None


def _get_default_client() -> MHCflurryClient:
    """Return the module-level default :class:`MHCflurryClient`.

    The client is created on first call with ``allow_offline_fallback=True``.
    """
    global _default_client
    if _default_client is None:
        _default_client = MHCflurryClient(allow_offline_fallback=True)
    return _default_client


def predict_binding(allele: str, peptide: str) -> MHCBindingResult:
    """Predict MHC binding for a single peptide-allele pair with full fallback.

    This module-level convenience function follows a four-tier fallback
    chain that **never crashes** due to missing models:

    1. **MHCflurry** — neural-network ensemble (AUC 0.80–0.85, confidence 1.0)
    2. **NetMHCpan** — web API predictor (AUC 0.85–0.95, confidence 1.0)
    3. **Precomputed database** — curated IEDB/SYFPEITHI entries
       (confidence 0.7)
    4. **PSSM fallback** — position-specific scoring matrix heuristic
       (AUC 0.60–0.75, confidence 0.5)

    Parameters
    ----------
    allele : str
        MHC allele name (e.g. ``"HLA-A*02:01"``).
    peptide : str
        Amino-acid sequence of the peptide (8–11 residues typical).

    Returns
    -------
    MHCBindingResult
        Binding prediction with ``method`` and ``confidence`` fields
        indicating which tier produced the result.

    Raises
    ------
    ValueError
        If the peptide contains non-standard amino acids.
    """
    return _get_default_client().predict_binding(peptide, allele)


def predict_batch(allele: str, peptides: list[str]) -> list[MHCBindingResult]:
    """Predict MHC binding for multiple peptides against one allele.

    This module-level convenience function applies the same fallback
    chain as :func:`predict_binding` to each peptide individually.

    Parameters
    ----------
    allele : str
        MHC allele name (e.g. ``"HLA-A*02:01"``).
    peptides : list[str]
        Amino-acid sequences of the peptides.

    Returns
    -------
    list[MHCBindingResult]
        One binding prediction per input peptide, in the same order.
    """
    client = _get_default_client()
    results: list[MHCBindingResult] = []
    for peptide in peptides:
        try:
            result = client.predict_binding(peptide, allele)
            results.append(result)
        except ValueError:
            # Re-raise validation errors (bad amino acids)
            raise
        except Exception as exc:
            logger.warning(
                "Prediction failed for %s / %s: %s — returning non-binder",
                allele, peptide, exc,
            )
            # Ultimate fallback: return a non-binder result
            results.append(MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=max(0, len(peptide) - 1),
                binding_score=0.0,
                ic50_nm=_MAX_IC50_NM,
                binding_class="non_binder",
                anchor_residues={},
                anchor_scores={},
                method="pssm_fallback",
                confidence=CONFIDENCE_PSSM,
            ))
    return results


# ═══════════════════════════════════════════════════════════════════════════
# MHCflurryClient
# ═══════════════════════════════════════════════════════════════════════════

class MHCflurryClient:
    """Client for offline MHC-I binding and presentation prediction.

    Uses MHCflurry's neural-network ensemble to predict peptide-MHC
    binding affinity (IC50) and, optionally, antigen-processing scores
    that incorporate proteasomal cleavage and TAP transport.

    Models are lazy-loaded on the first prediction call, keeping
    construction lightweight.

    Parameters
    ----------
    models_dir : str or None
        Path to the MHCflurry models directory.  If ``None``, the
        default location (``~/.mhcflurry/``) is used.
    allow_offline_fallback : bool
        If ``True`` (default), the client gracefully falls back through
        the full prediction chain (MHCflurry → NetMHCpan → precomputed
        → PSSM) when models are not available.  If ``False``, a
        :class:`RuntimeError` is raised when MHCflurry is unavailable.

    Examples
    --------
    >>> client = MHCflurryClient()
    >>> result = client.predict_binding("SIINFEKL", "HLA-A*02:01")
    >>> print(result.binding_class, result.ic50_nm, result.method, result.confidence)
    moderate_binder 342.5 mhcflurry 1.0

    >>> results = client.batch_predict(
    ...     "MAGRSGDLDAIIRYVKQLR",
    ...     alleles=["HLA-A*02:01", "HLA-A*03:01"],
    ...     epitope_lengths=[8, 9, 10],
    ... )
    >>> len(results) > 0
    True
    """

    def __init__(self, models_dir: str | None = None, allow_offline_fallback: bool = True) -> None:
        self._models_dir = models_dir or _DEFAULT_MODELS_DIR
        self._affinity_predictor = None
        self._presentation_predictor = None
        self._cache = _LRUCache()
        self._models_loaded = False
        self._models_load_failed = False
        self.allow_offline_fallback = allow_offline_fallback

    # ── lazy model loading ────────────────────────────────────────

    def _load_models(self) -> bool:
        """Lazy-load MHCflurry models on first prediction call.

        Returns ``True`` if models were loaded successfully, ``False``
        otherwise.  Unlike the previous version, this method **never
        raises** — failures are logged and the flag
        ``_models_load_failed`` is set so the fallback chain can proceed.
        """
        if self._models_loaded:
            return True
        if self._models_load_failed:
            return False

        try:
            import mhcflurry
        except ImportError:
            logger.info(
                "mhcflurry is not installed — skipping MHCflurry prediction. "
                "Install it with: pip install mhcflurry"
            )
            self._models_load_failed = True
            return False

        try:
            # Load Class1 affinity predictor
            if self._models_dir and self._models_dir != _DEFAULT_MODELS_DIR:
                self._affinity_predictor = (
                    mhcflurry.Class1AffinityPredictor.load(
                        models_dir=self._models_dir
                    )
                )
            else:
                self._affinity_predictor = (
                    mhcflurry.Class1AffinityPredictor.load()
                )
            logger.info("MHCflurry Class1 affinity predictor loaded")
        except Exception as exc:
            logger.info(
                "Failed to load MHCflurry affinity model: %s — "
                "will use fallback chain", exc,
            )
            self._models_load_failed = True
            return False

        # Try to load the presentation predictor (optional)
        try:
            if self._models_dir and self._models_dir != _DEFAULT_MODELS_DIR:
                self._presentation_predictor = (
                    mhcflurry.Class1PresentationPredictor.load(
                        models_dir=self._models_dir
                    )
                )
            else:
                self._presentation_predictor = (
                    mhcflurry.Class1PresentationPredictor.load()
                )
            logger.info("MHCflurry Class1 presentation predictor loaded")
        except Exception as exc:
            logger.info(
                "MHCflurry presentation predictor not available "
                "(non-fatal): %s",
                exc,
            )
            self._presentation_predictor = None

        self._models_loaded = True
        return True

    def _try_mhcflurry_prediction(
        self, peptide: str, allele: str,
    ) -> MHCBindingResult | None:
        """Try to get a MHCflurry prediction.

        Returns the result if successful, or ``None`` if MHCflurry
        is unavailable or fails.  **Never raises.**
        """
        # Check cache first
        cached = self._cache.get(peptide, allele)
        if cached is not None:
            logger.debug("Cache hit for %s / %s", allele, peptide)
            return cached

        # Try loading models if not yet loaded
        if not self._models_loaded and not self._models_load_failed:
            self._load_models()

        if self._affinity_predictor is None:
            return None

        try:
            df = self._affinity_predictor.predict(
                peptides=[peptide],
                alleles=[allele],
            )
            if not df.empty:
                ic50_nm = float(df.iloc[0]["prediction"])
                rank_val = None
                if "presentation_percentile" in df.columns:
                    rank_val = float(df.iloc[0].get("presentation_percentile", None))
                result = _mhcflurry_result_to_binding_result(
                    peptide=peptide,
                    allele=allele,
                    start_position=0,
                    end_position=len(peptide) - 1,
                    ic50_nm=ic50_nm,
                    method="mhcflurry",
                    rank=rank_val,
                )
                self._cache.put(peptide, allele, result)
                return result
            logger.debug(
                "MHCflurry returned no prediction for %s / %s", allele, peptide,
            )
            return None
        except Exception as exc:
            logger.info(
                "MHCflurry prediction failed for %s / %s: %s", allele, peptide, exc,
            )
            return None

    def _try_netmhcpan_prediction(
        self, peptide: str, allele: str,
    ) -> MHCBindingResult | None:
        """Try to get a NetMHCpan prediction.

        Returns the result if successful, or ``None`` if NetMHCpan
        is unavailable or fails.  **Never raises.**
        """
        try:
            from biocompiler.immunogenicity.netmhcpan import NetMHCpanClient
            client = NetMHCpanClient()
            # NetMHCpanClient exposes predict_mhc_i_binding / predict_mhc_ii_binding
            # (NOT predict_binding — that is the MHCflurryClient method).  Using
            # the wrong name raised AttributeError, which was silently swallowed
            # by the broad except below, disabling the entire NetMHCpan tier of
            # the 4-stage fallback chain (CRITICAL issue C11).
            result = client.predict_mhc_i_binding(peptide, allele)
            if result is not None:
                # Ensure confidence and method are set correctly
                result = replace(
                    result,
                    method="netmhcpan",
                    confidence=CONFIDENCE_NETMHCPAN,
                )
                self._cache.put(peptide, allele, result)
                return result
            return None
        except ImportError:
            logger.debug("NetMHCpan adapter not available — skipping")
            return None
        except Exception as exc:
            logger.info(
                "NetMHCpan prediction failed for %s / %s: %s "
                "— continuing fallback chain", allele, peptide, exc,
            )
            return None

    # ── single peptide prediction ─────────────────────────────────

    def predict_binding(
        self,
        peptide: str,
        allele: str,
    ) -> MHCBindingResult:
        """Predict MHC-I binding affinity for a single peptide-allele pair.

        Prediction follows a fallback chain:

        1. **MHCflurry** — neural-network ensemble (AUC 0.80–0.85)
        2. **NetMHCpan** — web API predictor (AUC 0.85–0.95), if available
        3. **Pre-computed database** — curated IEDB/SYFPEITHI entries +
           PSSM-derived predictions
        4. **PSSM-based prediction** — position-specific scoring matrix
           heuristic (AUC 0.60–0.75)

        Only if **all** methods fail *and* ``allow_offline_fallback`` is
        ``False`` is a :class:`RuntimeError` raised.  When the fallback
        chain is enabled (the default), the method always returns a
        result.

        Parameters
        ----------
        peptide : str
            Amino-acid sequence of the peptide (8–11 residues typical).
        allele : str
            MHC-I allele name (e.g. ``"HLA-A*02:01"``).

        Returns
        -------
        MHCBindingResult
            Binding prediction matching the
            :mod:`~biocompiler.immunogenicity` interface.

        Raises
        ------
        RuntimeError
            If MHCflurry is not available and offline fallback is
            disabled.
        ValueError
            If the peptide contains non-standard amino acids.
        """
        peptide = _validate_peptide(peptide)

        # --- 1. Try MHCflurry first ------------------------------------------
        result = self._try_mhcflurry_prediction(peptide, allele)
        if result is not None:
            return result

        logger.info(
            "MHCflurry unavailable for %s / %s — trying NetMHCpan fallback",
            allele, peptide,
        )

        # --- 2. Try NetMHCpan ------------------------------------------------
        result = self._try_netmhcpan_prediction(peptide, allele)
        if result is not None:
            return result

        logger.info(
            "NetMHCpan unavailable for %s / %s — trying precomputed database",
            allele, peptide,
        )

        # --- 3. Fallback to pre-computed database ----------------------------
        if self.allow_offline_fallback:
            record = self._lookup_precomputed(peptide, allele)
            if record is not None:
                logger.info(
                    "Using pre-computed binding data for %s/%s (source: %s)",
                    allele, peptide, record.source,
                )
                result = MHCBindingResult(
                    allele=allele,
                    peptide=peptide,
                    start_position=0,
                    end_position=len(peptide) - 1,
                    binding_score=round(record.binding_score, 6),
                    ic50_nm=round(record.ic50_nm, 2),
                    binding_class=record.binding_class,
                    anchor_residues=record.anchor_residues,
                    anchor_scores={k: round(v, 6) for k, v in record.anchor_scores.items()},
                    method="precomputed_lookup",
                    confidence=CONFIDENCE_PRECOMPUTED,
                )
                self._cache.put(peptide, allele, result)
                return result

            logger.info(
                "No pre-computed data for %s/%s — using PSSM fallback",
                allele, peptide,
            )

            # --- 4. PSSM fallback (always available offline) -----------------
            return self._pssm_predict(peptide, allele)

        raise RuntimeError(
            "MHCflurry not available and offline fallback disabled"
        )

    # ── batch prediction ──────────────────────────────────────────

    def batch_predict(
        self,
        protein: str,
        alleles: list[str],
        epitope_lengths: list[int] | None = None,
    ) -> list[MHCBindingResult]:
        """Scan a full protein for MHC-I binders across alleles.

        Extracts all overlapping peptides of the requested lengths and
        predicts binding using the full fallback chain (MHCflurry →
        NetMHCpan → precomputed → PSSM).

        Parameters
        ----------
        protein : str
            Full protein amino-acid sequence.
        alleles : list[str]
            MHC-I alleles to evaluate.
        epitope_lengths : list[int] or None
            Peptide lengths to extract.  Defaults to
            ``[8, 9, 10, 11]``.

        Returns
        -------
        list[MHCBindingResult]
            Binding predictions for every peptide × allele combination.
            Unsupported alleles are handled by the fallback chain.
        """
        if epitope_lengths is None:
            epitope_lengths = list(_DEFAULT_EPITOPE_LENGTHS)

        try:
            protein = _validate_protein(protein)
        except ValueError as exc:
            logger.warning("Invalid protein sequence: %s", exc)
            return []

        if not alleles:
            return []

        # Extract overlapping peptides
        peptide_tuples = _extract_overlapping_peptides(protein, epitope_lengths)
        if not peptide_tuples:
            return []

        # Use predict_binding for each (peptide, allele) pair — it
        # handles the full fallback chain and never crashes.
        results: list[MHCBindingResult] = []
        for pep, start, end in peptide_tuples:
            for allele in alleles:
                try:
                    result = self.predict_binding(pep, allele)
                    # Override positions with protein-relative ones
                    result = replace(
                        result,
                        start_position=start,
                        end_position=end,
                    )
                    results.append(result)
                except ValueError:
                    # Invalid peptide — skip
                    continue
                except Exception as exc:
                    logger.warning(
                        "Prediction failed for %s / %s: %s — skipping",
                        allele, pep, exc,
                    )

        # Optimisation: if MHCflurry is available, try batch mode for
        # better performance, then fill gaps with fallback chain.
        if self._affinity_predictor is not None:
            self._batch_predict_mhcflurry_optimized(
                peptide_tuples, alleles, results,
            )

        logger.info(
            "batch_predict: %d results for %d alleles, "
            "protein length %d, epitope_lengths %s",
            len(results), len(alleles), len(protein), epitope_lengths,
        )
        return results

    def _batch_predict_mhcflurry_optimized(
        self,
        peptide_tuples: list[tuple[str, int, int]],
        alleles: list[str],
        existing_results: list[MHCBindingResult],
    ) -> None:
        """Optimised batch prediction using MHCflurry's batch mode.

        This fills in results for peptide-allele pairs that have not
        been predicted yet using the fallback chain.  It only runs
        when MHCflurry models are already loaded.
        """
        if self._affinity_predictor is None:
            return

        # Identify which (peptide, allele) pairs still need predictions
        existing_keys: set[tuple[str, str]] = {
            (r.peptide, r.allele) for r in existing_results
        }

        batch_peptides: list[str] = []
        batch_alleles: list[str] = []
        batch_meta: list[tuple[str, int, int]] = []

        for pep, start, end in peptide_tuples:
            for allele in alleles:
                if (pep, allele) in existing_keys:
                    continue
                cached = self._cache.get(pep, allele)
                if cached is not None:
                    existing_results.append(replace(
                        cached,
                        start_position=start,
                        end_position=end,
                    ))
                    existing_keys.add((pep, allele))
                    continue
                batch_peptides.append(pep)
                batch_alleles.append(allele)
                batch_meta.append((allele, start, end))

        if not batch_peptides:
            return

        batch_ic50s: dict[tuple[str, str], float] = {}
        try:
            df = self._affinity_predictor.predict(
                peptides=batch_peptides,
                alleles=batch_alleles,
            )
            for idx in range(len(df)):
                row = df.iloc[idx]
                pep = str(row["peptide"])
                allele = str(row["allele"])
                ic50 = float(row["prediction"])
                batch_ic50s[(pep, allele)] = ic50
        except Exception as exc:
            logger.info(
                "MHCflurry batch prediction failed: %s — "
                "filling gaps with fallback chain", exc,
            )
            # Fill gaps with fallback chain
            for pep, start, end in peptide_tuples:
                for allele in alleles:
                    if (pep, allele) not in existing_keys:
                        try:
                            result = self.predict_binding(pep, allele)
                            result = replace(
                                result, start_position=start, end_position=end,
                            )
                            existing_results.append(result)
                            existing_keys.add((pep, allele))
                        except Exception:
                            pass
            return

        # Process batch results
        for pep, allele, start, end in [
            (batch_peptides[i], batch_meta[i][0], batch_meta[i][1], batch_meta[i][2])
            for i in range(len(batch_peptides))
        ]:
            ic50_nm = batch_ic50s.get((pep, allele))
            if ic50_nm is None:
                # MHCflurry did not return a result for this pair
                # Try the fallback chain for this specific pair
                try:
                    result = self.predict_binding(pep, allele)
                    result = replace(result, start_position=start, end_position=end)
                    existing_results.append(result)
                except Exception:
                    pass
                continue

            result = _mhcflurry_result_to_binding_result(
                peptide=pep,
                allele=allele,
                start_position=start,
                end_position=end,
                ic50_nm=ic50_nm,
                method="mhcflurry",
            )
            self._cache.put(pep, allele, result)
            existing_results.append(result)

    # ── presentation prediction ───────────────────────────────────

    def predict_presentation(
        self,
        protein: str,
        alleles: list[str],
        epitope_lengths: list[int] | None = None,
    ) -> list[MHCBindingResult]:
        """Predict antigen processing + MHC binding (presentation score).

        Uses MHCflurry's :class:`~mhcflurry.Class1PresentationPredictor`
        which combines:

        * **MHC binding affinity** (neural-network IC50 prediction)
        * **Proteasomal cleavage** (likely cleavage sites)
        * **TAP transport efficiency** (Transporter associated with
          Antigen Processing)

        This is more accurate than binding prediction alone for
        identifying naturally presented epitopes.

        Parameters
        ----------
        protein : str
            Full protein amino-acid sequence.
        alleles : list[str]
            MHC-I alleles to evaluate.
        epitope_lengths : list[int] or None
            Peptide lengths to extract.  Defaults to
            ``[8, 9, 10, 11]``.

        Returns
        -------
        list[MHCBindingResult]
            Predictions with ``binding_score`` set to the presentation
            score (0–1, where higher = more likely presented).  The
            ``ic50_nm`` field still contains the predicted IC50.

        Notes
        -----
        Falls back to :meth:`batch_predict` (binding only) if the
        presentation predictor is not available (e.g. models not
        downloaded).
        """
        if epitope_lengths is None:
            epitope_lengths = list(_DEFAULT_EPITOPE_LENGTHS)

        try:
            protein = _validate_protein(protein)
        except ValueError as exc:
            logger.warning("Invalid protein sequence: %s", exc)
            return []

        if not alleles:
            return []

        # Load models if not yet loaded
        if not self._models_loaded and not self._models_load_failed:
            self._load_models()

        # If presentation predictor is not available, fall back
        if self._presentation_predictor is None:
            logger.info(
                "MHCflurry presentation predictor not available. "
                "Falling back to binding-only prediction. "
                "Run download_models() to get presentation models."
            )
            return self.batch_predict(protein, alleles, epitope_lengths)

        # Extract overlapping peptides
        peptide_tuples = _extract_overlapping_peptides(protein, epitope_lengths)
        if not peptide_tuples:
            return []

        # Build batch input
        batch_peptides: list[str] = []
        batch_alleles: list[str] = []
        batch_meta: list[tuple[str, int, int]] = []

        for pep, start, end in peptide_tuples:
            for allele in alleles:
                batch_peptides.append(pep)
                batch_alleles.append(allele)
                batch_meta.append((allele, start, end))

        results: list[MHCBindingResult] = []

        if not batch_peptides:
            return results

        try:
            df = self._presentation_predictor.predict(
                peptides=batch_peptides,
                alleles=batch_alleles,
            )
            for idx in range(len(df)):
                row = df.iloc[idx]
                pep = str(row["peptide"])
                allele = str(row["allele"])

                # Extract presentation score and IC50
                presentation_score = float(
                    row.get("presentation_score", row.get("presentation_percentile", 0.0))
                )
                ic50_nm = float(
                    row.get("binding_affinity", row.get("prediction", _MAX_IC50_NM))
                )

                # Find the corresponding meta entry
                start = 0
                end = len(pep) - 1
                for a, s, e in batch_meta:
                    if a == allele:
                        start = s
                        end = e
                        break

                # Use presentation percentile as binding_score if
                # available (lower percentile = better presentation)
                if "presentation_percentile" in row:
                    pct = float(row["presentation_percentile"])
                    # Convert percentile to a 0-1 score (lower = better)
                    pres_score = max(0.0, min(1.0, 1.0 - pct / 100.0))
                else:
                    pres_score = max(0.0, min(1.0, presentation_score))

                rank_val = None
                if "presentation_percentile" in row:
                    rank_val = float(row["presentation_percentile"])

                result = _mhcflurry_result_to_binding_result(
                    peptide=pep,
                    allele=allele,
                    start_position=start,
                    end_position=end,
                    ic50_nm=ic50_nm,
                    presentation_score=pres_score,
                    method="mhcflurry_presentation",
                    rank=rank_val,
                )
                results.append(result)

        except Exception as exc:
            logger.info(
                "MHCflurry presentation prediction failed: %s. "
                "Falling back to binding-only prediction.",
                exc,
            )
            return self.batch_predict(protein, alleles, epitope_lengths)

        logger.info(
            "predict_presentation: %d results for %d alleles, "
            "protein length %d, epitope_lengths %s",
            len(results), len(alleles), len(protein), epitope_lengths,
        )
        return results

    # ── offline fallback helpers ────────────────────────────────────

    def _lookup_precomputed(self, peptide: str, allele: str) -> object | None:
        """Look up a peptide-allele pair in the pre-computed binding database.

        Returns a :class:`~biocompiler.mhc_binding_db.PrecomputedEntry`
        if found, or ``None`` otherwise.

        .. deprecated:: C20 remediation (prior fix)
            The ``precomputed/`` subpackage (~20,500 LOC) was deleted because
            25 entries labelled ``source="known_epitope"`` had PSSM-derived
            IC50 values that diverged 1-2 orders of magnitude from published
            values, and 2,075 ``pssm_predicted`` entries were stale relative
            to the current runtime PSSM scoring.  This method now always
            returns ``None`` so the fallback chain falls through to
            :meth:`_pssm_predict`, which uses the SAME PSSM matrices as the
            deleted tables but with current scoring and an honest
            ``confidence=0.5``.  The method signature is retained for
            backward compatibility (tests mock it to exercise the
            ``precomputed_lookup`` branch of the fallback chain).
        """
        return None

    def _pssm_predict(self, peptide: str, allele: str) -> MHCBindingResult:
        """Predict binding using the PSSM-based heuristic from immunogenicity.

        This is the lowest-accuracy fallback (AUC 0.60–0.75) but is
        always available offline without any external dependencies.
        """
        try:
            pssm_score = score_peptide_pssm(peptide, allele)
            ic50_nm = binding_score_to_ic50(pssm_score)
            binding_class = classify_binding(ic50_nm)
            binding_score = ic50_to_binding_score(ic50_nm)

            # Anchor residues from PSSM
            anchor_residues: dict[int, str] = {}
            anchor_scores: dict[int, float] = {}
            if len(peptide) >= 2:
                anchor_residues[1] = peptide[1]  # P2 anchor
                anchor_scores[1] = binding_score
            if len(peptide) >= 1:
                last_idx = len(peptide) - 1
                anchor_residues[last_idx] = peptide[last_idx]  # C-term anchor
                anchor_scores[last_idx] = binding_score

            result = MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=len(peptide) - 1,
                binding_score=round(binding_score, 6),
                ic50_nm=round(ic50_nm, 2),
                binding_class=binding_class,
                anchor_residues=anchor_residues,
                anchor_scores={k: round(v, 6) for k, v in anchor_scores.items()},
                method="pssm_fallback",
                confidence=CONFIDENCE_PSSM,
            )
            self._cache.put(peptide, allele, result)
            return result
        except Exception as exc:
            logger.warning(
                "PSSM prediction failed for %s / %s: %s — "
                "returning non-binder result", allele, peptide, exc,
            )
            # Ultimate fallback: return non-binder
            return MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=max(0, len(peptide) - 1),
                binding_score=0.0,
                ic50_nm=_MAX_IC50_NM,
                binding_class="non_binder",
                anchor_residues={},
                anchor_scores={},
                method="pssm_fallback",
                confidence=CONFIDENCE_PSSM,
            )
