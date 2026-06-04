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
- :func:`download_models` — download MHCflurry models (~100 MB)
- LRU prediction cache (max 5 000 entries)
- Graceful error handling: unsupported alleles are skipped with a debug
  log; prediction errors produce empty results with a warning

Result format
-------------
MHCflurry returns IC50 values (nM) and presentation percentiles.
This adapter converts them to the same :class:`MHCBindingResult` format
used by :mod:`biocompiler.immunogenicity`:

- ``binding_score`` = 1 − log(IC50) / log(50 000)  (clamped to [0, 1])
- ``ic50_nm`` = raw predicted IC50
- ``binding_class`` = classify_binding(ic50)  using the standard
  50 / 500 / 5 000 nM thresholds

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
from typing import Sequence

from .immunogenicity import MHCBindingResult, classify_binding, score_peptide_pssm, binding_score_to_ic50

__all__ = [
    "MHCflurryClient",
    "is_mhcflurry_available",
    "download_models",
    "clear_cache",
    "MHCFLURRY_AUC_ROC_LOW",
    "MHCFLURRY_AUC_ROC_HIGH",
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
# Availability check
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

    Returns
    -------
    MHCBindingResult
    """
    binding_score = ic50_to_binding_score(ic50_nm)
    if presentation_score is not None:
        binding_score = presentation_score
    binding_class = classify_binding(ic50_nm)

    # Anchor residues: for MHCflurry we don't have PSSM-style per-position
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
    )


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

    Examples
    --------
    >>> client = MHCflurryClient()
    >>> result = client.predict_binding("SIINFEKL", "HLA-A*02:01")
    >>> print(result.binding_class, result.ic50_nm)
    moderate_binder 342.5

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
        self.allow_offline_fallback = allow_offline_fallback

    # ── lazy model loading ────────────────────────────────────────

    def _load_models(self) -> None:
        """Lazy-load MHCflurry models on first prediction call.

        Sets ``self._models_loaded = True`` on success.  If loading
        fails, logs an error and leaves ``_models_loaded`` as ``False``
        so subsequent calls will retry.
        """
        if self._models_loaded:
            return

        try:
            import mhcflurry
        except ImportError as exc:
            logger.error(
                "mhcflurry is not installed. Install it with: "
                "pip install mhcflurry"
            )
            raise

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
            logger.error("Failed to load MHCflurry affinity model: %s", exc)
            raise

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
            logger.warning(
                "MHCflurry presentation predictor not available "
                "(non-fatal): %s",
                exc,
            )
            self._presentation_predictor = None

        self._models_loaded = True

    def _ensure_models(self) -> None:
        """Ensure models are loaded, raising on failure."""
        if not self._models_loaded:
            self._load_models()
        if self._affinity_predictor is None:
            raise RuntimeError(
                "MHCflurry affinity predictor is not available. "
                "Run `download_models()` first or check your installation."
            )

    # ── single peptide prediction ─────────────────────────────────

    def predict_binding(
        self,
        peptide: str,
        allele: str,
    ) -> MHCBindingResult:
        """Predict MHC-I binding affinity for a single peptide-allele pair.

        Prediction follows a fallback chain:

        1. **MHCflurry** — neural-network ensemble (AUC 0.80–0.85)
        2. **Pre-computed database** — curated IEDB/SYFPEITHI entries +
           PSSM-derived predictions
        3. **PSSM-based prediction** — position-specific scoring matrix
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
            disabled, or if *all* prediction methods fail.
        ValueError
            If the peptide contains non-standard amino acids.
        """
        peptide = _validate_peptide(peptide)

        # Check cache
        cached = self._cache.get(peptide, allele)
        if cached is not None:
            logger.debug("Cache hit for %s / %s", allele, peptide)
            return cached

        # --- 1. Try MHCflurry first ------------------------------------------
        if is_mhcflurry_available() and self._affinity_predictor is not None:
            try:
                self._ensure_models()
                df = self._affinity_predictor.predict(
                    peptides=[peptide],
                    alleles=[allele],
                )
                if not df.empty:
                    ic50_nm = float(df.iloc[0]["prediction"])
                    result = _mhcflurry_result_to_binding_result(
                        peptide=peptide,
                        allele=allele,
                        start_position=0,
                        end_position=len(peptide) - 1,
                        ic50_nm=ic50_nm,
                    )
                    self._cache.put(peptide, allele, result)
                    return result
                logger.debug(
                    "MHCflurry returned no prediction for %s / %s "
                    "— trying offline fallback", allele, peptide,
                )
            except RuntimeError:
                raise
            except Exception as exc:
                logger.debug(
                    "MHCflurry prediction failed for %s / %s: %s "
                    "— trying offline fallback", allele, peptide, exc,
                )
        elif is_mhcflurry_available():
            # MHCflurry is installed but models are not yet loaded;
            # try loading them.
            try:
                self._ensure_models()
                # Retry after loading
                return self.predict_binding(peptide, allele)
            except Exception as exc:
                logger.debug(
                    "MHCflurry model loading failed: %s "
                    "— trying offline fallback", exc,
                )
        else:
            logger.debug(
                "MHCflurry not available — using offline fallback "
                "for %s / %s", allele, peptide,
            )

        # --- 2. Fallback to pre-computed database ----------------------------
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
                )
                self._cache.put(peptide, allele, result)
                return result
            else:
                logger.warning(
                    "No pre-computed data for %s/%s, using PSSM prediction",
                    allele, peptide,
                )
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
        runs MHCflurry in batch mode (significantly faster than calling
        :meth:`predict_binding` individually).

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
            Unsupported alleles are silently skipped.
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

        self._ensure_models()

        # Extract overlapping peptides
        peptide_tuples = _extract_overlapping_peptides(protein, epitope_lengths)
        if not peptide_tuples:
            return []

        # Build batch input: one entry per (peptide, allele) pair
        batch_peptides: list[str] = []
        batch_alleles: list[str] = []
        batch_meta: list[tuple[str, int, int]] = []  # (allele, start, end)

        for pep, start, end in peptide_tuples:
            for allele in alleles:
                # Check cache first
                cached = self._cache.get(pep, allele)
                if cached is not None:
                    continue  # will be added from cache later
                batch_peptides.append(pep)
                batch_alleles.append(allele)
                batch_meta.append((allele, start, end))

        # Run batch prediction
        batch_ic50s: dict[tuple[str, str], float] = {}
        if batch_peptides:
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
                logger.warning(
                    "MHCflurry batch prediction failed: %s. "
                    "Returning cached results only.",
                    exc,
                )

        # Assemble results
        results: list[MHCBindingResult] = []

        # Collect cached results
        for pep, start, end in peptide_tuples:
            for allele in alleles:
                cached = self._cache.get(pep, allele)
                if cached is not None:
                    # Override positions with current protein-relative ones
                    results.append(replace(
                        cached,
                        start_position=start,
                        end_position=end,
                    ))

        # Collect batch results
        for pep, start, end in peptide_tuples:
            for allele in alleles:
                # Skip if already added from cache
                if self._cache.get(pep, allele) is not None:
                    continue

                ic50_nm = batch_ic50s.get((pep, allele))

                if ic50_nm is None:
                    # Allele may not be supported by MHCflurry
                    logger.debug(
                        "No MHCflurry prediction for %s / %s — skipping",
                        allele, pep,
                    )
                    continue

                result = _mhcflurry_result_to_binding_result(
                    peptide=pep,
                    allele=allele,
                    start_position=start,
                    end_position=end,
                    ic50_nm=ic50_nm,
                )
                self._cache.put(pep, allele, result)
                results.append(result)

        logger.info(
            "batch_predict: %d results for %d alleles, "
            "protein length %d, epitope_lengths %s",
            len(results), len(alleles), len(protein), epitope_lengths,
        )
        return results

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

        # If presentation predictor is not available, fall back
        if self._presentation_predictor is None:
            logger.warning(
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

                result = _mhcflurry_result_to_binding_result(
                    peptide=pep,
                    allele=allele,
                    start_position=start,
                    end_position=end,
                    ic50_nm=ic50_nm,
                    presentation_score=pres_score,
                    method="mhcflurry_presentation",
                )
                results.append(result)

        except Exception as exc:
            logger.warning(
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
        """
        try:
            from .mhc_binding_db import get_database
            db = get_database(allele)
            if db is None:
                return None
            return db.lookup(peptide)
        except Exception as exc:
            logger.debug(
                "Pre-computed database lookup failed for %s/%s: %s",
                allele, peptide, exc,
            )
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

            # Compute anchor residues using the PSSM
            anchor_residues: dict[int, str] = {}
            anchor_scores: dict[int, float] = {}
            from .immunogenicity import _get_mhc_i_pssms, _identify_anchor_positions
            pssm_data = _get_mhc_i_pssms().get(allele)
            if pssm_data is not None and len(peptide) == len(pssm_data):
                anchor_residues, anchor_scores = _identify_anchor_positions(
                    peptide, pssm_data
                )
            else:
                # Fallback: use canonical MHC-I anchor positions
                if len(peptide) >= 2:
                    anchor_residues[1] = peptide[1]  # P2
                    anchor_scores[1] = pssm_score
                if len(peptide) >= 1:
                    last_idx = len(peptide) - 1
                    anchor_residues[last_idx] = peptide[last_idx]
                    anchor_scores[last_idx] = pssm_score

            result = MHCBindingResult(
                allele=allele,
                peptide=peptide,
                start_position=0,
                end_position=len(peptide) - 1,
                binding_score=round(pssm_score, 6),
                ic50_nm=round(ic50_nm, 2),
                binding_class=binding_class,
                anchor_residues=anchor_residues,
                anchor_scores={k: round(v, 6) for k, v in anchor_scores.items()},
                method="pssm_fallback",
            )
            self._cache.put(peptide, allele, result)
            return result
        except Exception as exc:
            raise RuntimeError(
                f"All prediction methods failed for allele "
                f"{allele!r} / peptide {peptide!r}: MHCflurry unavailable, "
                f"pre-computed data not found, and PSSM fallback error: {exc}"
            ) from exc

    # ── cache management ──────────────────────────────────────────

    def clear_cache(self) -> None:
        """Clear the prediction cache for this client instance."""
        self._cache.clear()
        logger.debug("MHCflurryClient cache cleared")

    @property
    def cache_size(self) -> int:
        """Number of entries in the client's prediction cache."""
        return self._cache.size

    @property
    def cache_hit_rate(self) -> float:
        """Cache hit rate (0.0 – 1.0) for this client instance."""
        return self._cache.hit_rate

    # ── supported alleles ─────────────────────────────────────────

    def supported_alleles(self) -> list[str]:
        """Return the list of alleles supported by the loaded model.

        Returns
        -------
        list[str]
            Allele names supported by the MHCflurry affinity predictor.

        Raises
        ------
        RuntimeError
            If models have not been loaded yet.
        """
        self._ensure_models()
        try:
            return sorted(self._affinity_predictor.supported_alleles)
        except AttributeError:
            # Older MHCflurry versions may not expose this
            logger.debug(
                "MHCflurry predictor does not expose supported_alleles"
            )
            return []

    def is_allele_supported(self, allele: str) -> bool:
        """Check whether a specific allele is supported.

        Parameters
        ----------
        allele : str
            MHC-I allele name.

        Returns
        -------
        bool
        """
        self._ensure_models()
        try:
            return allele in self._affinity_predictor.supported_alleles
        except AttributeError:
            # Older MHCflurry versions may not expose supported_alleles.
            # Attempt a real prediction with a test peptide to check allele support.
            try:
                df = self._affinity_predictor.predict(
                    peptides=["AAAAAAAAA"],
                    alleles=[allele],
                )
                return not df.empty
            except Exception:
                logger.debug(
                    "Cannot determine allele support for %r via prediction probe",
                    allele,
                )
                return False
