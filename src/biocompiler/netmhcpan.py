"""BioCompiler NetMHCpan Client — MHC Binding Prediction via NetMHCpan API.

Provides a client for the NetMHCpan 4.1 web API at
https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
for peptide-MHC binding affinity prediction.

This module provides an upgrade path from the PSSM-based heuristics in
:mod:`biocompiler.immunogenicity`.  When the NetMHCpan API is available,
predictions are more accurate; when it is unavailable, the system
gracefully falls back to PSSM-based scoring with a warning.

Features:
  - NetMHCpanClient: Class that calls the NetMHCpan web API
  - predict_mhc_i_binding: MHC-I binding prediction for a single peptide/allele
  - predict_mhc_ii_binding: MHC-II binding prediction (via NetMHCIIpan)
  - batch_predict: Batch prediction for a full protein across multiple alleles
  - NetMHCpanCache: In-memory caching for API results
  - Retry logic with exponential backoff
  - Timeout handling
  - Output parsing for NetMHCpan text format
  - Binding classification: Strong binder (rank < 0.5%),
    Weak binder (rank < 2%), No binding

All network calls use ``urllib.request`` (consistent with :mod:`esmfold`).
Individual prediction failures are isolated — they never crash sibling
predictions in batch mode.

References
----------
- NetMHCpan 4.1: Jurtz et al., J Immunol 2017; 199:3360-3368
- NetMHCIIpan 4.0: Reynisson et al., Nucleic Acids Res 2020; 48:W449-W456

MHC Prediction Fallback Chain
-----------------------------
When predicting MHC binding affinity, the system follows this priority:

1. **MHCflurry** (if installed) — fastest local predictor, good accuracy
2. **NetMHCpan** (if installed locally or web API reachable) — gold-standard
   predictor from DTU Health Tech
3. **Precomputed database** — lookup from pre-computed binding affinity
   databases (e.g. IEDB)
4. **PSSM fallback** — position-specific scoring matrix heuristic,
   always available as a last resort

The adapter functions :func:`is_netmhcpan_available`,
:func:`predict_binding_netmhcpan`, and
:func:`batch_predict_binding_netmhcpan` provide a clean interface for
the fallback chain.  They return ``None`` or empty results when
NetMHCpan is not available, allowing the caller to try the next
fallback without raising exceptions.
"""
from __future__ import annotations

import hashlib
import logging
import re
import shutil
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Optional
from .constants import DEFAULT_ENGINE_TIMEOUT
from .engine_base import validate_protein_sequence
from .exceptions import ImmunogenicityError

__all__ = [
    "NetMHCpanError",
    "MHCBindingResult",
    "NetMHCpanCache",
    "NetMHCpanClient",
    "DEFAULT_API_URL",
    "DEFAULT_MHCII_API_URL",
    "DEFAULT_TIMEOUT",
    "MAX_RETRIES",
    "RETRY_BASE_DELAY",
    "POLL_INTERVAL",
    "MAX_POLL_ATTEMPTS",
    "STRONG_BINDER_RANK_THRESHOLD",
    "WEAK_BINDER_RANK_THRESHOLD",
    "MHC_II_EPITOPE_LENGTH",
    "DEFAULT_MHC_I_EPITOPE_LENGTHS",
    "clear_cache",
    "is_netmhcpan_installed",
    "is_netmhcpan_available",
    "predict_binding_netmhcpan",
    "batch_predict_binding_netmhcpan",
    "predict_mhc_i_binding",
    "predict_mhc_ii_binding",
    "batch_predict",
    "parse_netmhcpan_output",
    "classify_binding_rank",
]

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════

DEFAULT_API_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
DEFAULT_MHCII_API_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
DEFAULT_TIMEOUT: float = DEFAULT_ENGINE_TIMEOUT
MAX_RETRIES: int = 3
RETRY_BASE_DELAY: float = 2.0  # seconds, doubled each attempt
POLL_INTERVAL: float = 5.0  # seconds between polling for job completion
MAX_POLL_ATTEMPTS: int = 120  # 120 * 5s = 10 minutes max wait

# Binding classification thresholds (rank %)
STRONG_BINDER_RANK_THRESHOLD: float = 0.5  # rank < 0.5% → strong binder
WEAK_BINDER_RANK_THRESHOLD: float = 2.0  # rank < 2% → weak binder

# Epitope length defaults
MHC_II_EPITOPE_LENGTH: int = 15  # 15-mer peptides for MHC-II prediction
DEFAULT_MHC_I_EPITOPE_LENGTHS: list[int] = [8, 9, 10, 11]  # Standard MHC-I peptide lengths

# NetMHCpan IC50 conversion constant (50000 nM = 50 µM, the log50k reference)
_NETMHCPAN_LOG50K_CONSTANT: float = 50000.0

# Rank-to-score mapping divisor (log10 scale: rank 0.01→1.0, rank 100→0.0)
_RANK_SCORE_LOG10_DIVISOR: float = 4.0

# Fallback rank used when parsing yields no rank value
_FALLBACK_RANK: float = 100.0

# Minimum number of whitespace-separated tokens required to attempt parsing
_MIN_DATA_LINE_TOKENS: int = 6

# Maximum characters to read from an HTTP error response body
_ERROR_BODY_MAX_CHARS: int = 500

# NetMHCpan config file identifiers for the webface2 CGI
# NOTE: These look like local filesystem paths but are actually external service
# identifiers required by the DTU NetMHCpan webface2 CGI API. They reference
# config paths on the remote DTU server, NOT local files.
_NETMHCPAN_CONFIG: str = "/usr/opt/www/pub/CBS/services/NetMHCpan-4.1/webface.NetMHCpan-4.1.cfg"
_NETMHCIIPAN_CONFIG: str = "/usr/opt/www/pub/CBS/services/NetMHCIIpan-4.0/webface.NetMHCIIpan-4.0.cfg"


# ═══════════════════════════════════════════════════════════════════════════
# Exceptions
# ═══════════════════════════════════════════════════════════════════════════

class NetMHCpanError(ImmunogenicityError):
    """Raised when NetMHCpan API call fails.

    This covers network errors, API response parsing failures,
    timeout conditions, and invalid input errors.
    """

    def __init__(self, reason: str, allele: str | None = None):
        self.reason = reason
        self.allele = allele
        self._message = f"NetMHCpan error: {reason}"
        if allele:
            self._message += f" (allele={allele})"
        super().__init__(self._message)

    def __str__(self) -> str:
        return self._message


# ═══════════════════════════════════════════════════════════════════════════
# Data classes
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class MHCBindingResult:
    """Result of a single peptide-MHC binding prediction from NetMHCpan.

    Compatible with the ``MHCBindingResult`` from :mod:`immunogenicity`
    but extended with NetMHCpan-specific fields.

    Attributes
    ----------
    allele : str
        MHC allele name (e.g. ``"HLA-A*02:01"``).
    peptide : str
        Amino acid sequence of the predicted peptide.
    start_position : int
        0-based start position in the source protein.
    end_position : int
        0-based end position (inclusive) in the source protein.
    binding_score : float
        Normalised binding score in [0, 1] (1 = strongest).
    ic50_nm : float or None
        Predicted IC50 in nM, if available.
    binding_class : str
        One of ``"strong_binder"``, ``"weak_binder"``, ``"non_binder"``.
    rank : float or None
        Percentile rank from NetMHCpan, if available.
    anchor_residues : dict
        Position → amino acid at anchor positions.
    anchor_scores : dict
        Position → binding contribution score.
    method : str
        ``"netmhcpan"`` or ``"pssm_fallback"``.
    """
    allele: str
    peptide: str
    start_position: int
    end_position: int
    binding_score: float = 0.0
    ic50_nm: float | None = None
    binding_class: str = "non_binder"
    rank: float | None = None
    anchor_residues: dict = field(default_factory=dict)
    anchor_scores: dict = field(default_factory=dict)
    method: str = "netmhcpan"


# ═══════════════════════════════════════════════════════════════════════════
# Cache
# ═══════════════════════════════════════════════════════════════════════════

class NetMHCpanCache:
    """In-memory cache for NetMHCpan prediction results.

    Follows the same pattern as :class:`ESMFoldCache` for consistency.
    Cache key is derived from allele + peptide + length.
    """

    def __init__(self, max_size: int = 5000):
        self._cache: dict[str, MHCBindingResult] = {}
        self._max_size = max_size
        self._hits = 0
        self._misses = 0

    @staticmethod
    def _key(allele: str, peptide: str, epitope_length: int) -> str:
        """Generate a cache key from allele, peptide, and length."""
        raw = f"{allele}:{peptide}:{epitope_length}"
        return hashlib.sha256(raw.encode()).hexdigest()[:16]

    def get(self, allele: str, peptide: str, epitope_length: int = 9) -> MHCBindingResult | None:
        """Retrieve a cached prediction result.

        Returns None on cache miss.
        """
        key = self._key(allele, peptide, epitope_length)
        if key in self._cache:
            self._hits += 1
            return self._cache[key]
        self._misses += 1
        return None

    def put(self, allele: str, peptide: str, result: MHCBindingResult, epitope_length: int = 9) -> None:
        """Store a prediction result in the cache."""
        key = self._key(allele, peptide, epitope_length)
        if len(self._cache) >= self._max_size and key not in self._cache:
            # Evict oldest entry (FIFO)
            oldest_key = next(iter(self._cache))
            del self._cache[oldest_key]
        self._cache[key] = result

    def put_batch(self, results: list[MHCBindingResult], epitope_length: int = 9) -> None:
        """Store a batch of results in the cache."""
        for r in results:
            self.put(r.allele, r.peptide, r, epitope_length)

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
        """Number of entries in the cache."""
        return len(self._cache)

    @property
    def hit_rate(self) -> float:
        """Cache hit rate (0.0 to 1.0)."""
        total = self._hits + self._misses
        return self._hits / total if total > 0 else 0.0

    def clear(self) -> None:
        """Clear the cache."""
        self._cache.clear()
        self._hits = 0
        self._misses = 0


# Module-level default cache instance (lazy-initialized)
_default_cache: NetMHCpanCache | None = None


def _get_default_cache() -> NetMHCpanCache:
    """Return the module-level default cache (lazy-initialized)."""
    global _default_cache
    if _default_cache is None:
        _default_cache = NetMHCpanCache()
    return _default_cache


def clear_cache() -> None:
    """Clear the NetMHCpan prediction cache."""
    _get_default_cache().clear()
    logger.info("NetMHCpan prediction cache cleared")


# ═══════════════════════════════════════════════════════════════════════════
# Binding classification
# ═══════════════════════════════════════════════════════════════════════════

def classify_binding_rank(rank: float) -> str:
    """Classify a peptide by its NetMHCpan rank score.

    Parameters
    ----------
    rank : float
        Percentile rank from NetMHCpan.

    Returns
    -------
    str
        One of ``"strong_binder"`` (rank < 0.5%),
        ``"weak_binder"`` (rank < 2%), or ``"non_binder"``.
    """
    if rank < STRONG_BINDER_RANK_THRESHOLD:
        return "strong_binder"
    elif rank < WEAK_BINDER_RANK_THRESHOLD:
        return "weak_binder"
    else:
        return "non_binder"


def _rank_to_binding_score(rank: float) -> float:
    """Convert a percentile rank to a normalised binding score in [0, 1].

    Lower rank → higher binding score.
    Uses a log-scale mapping for better dynamic range:
      rank ~0.01% → score ~1.0
      rank ~0.5%  → score ~0.8
      rank ~2%    → score ~0.5
      rank ~10%   → score ~0.2
      rank >50%   → score ~0.0
    """
    if rank <= 0:
        return 1.0
    import math
    # log10 mapping: rank 0.01 → 1.0, rank 100 → 0.0
    score = max(0.0, min(1.0, 1.0 - math.log10(rank) / _RANK_SCORE_LOG10_DIVISOR))
    return round(score, 6)


# ═══════════════════════════════════════════════════════════════════════════
# Output parsing
# ═══════════════════════════════════════════════════════════════════════════

# NetMHCpan output line pattern
# Example lines (whitespace-separated columns):
#   1  HLA-A*02:01  AAAAAAAAV  ...  0.025  <0.5%  97.5
#   1  HLA-A*02:01  AAAAAAAAV  ...  0.025  WB     1.5
#   1  HLA-A*02:01  AAAAAAAAV  ...  0.025         5.0
#
# The output format varies slightly between versions but the key columns are:
#   - Position (1-based)
#   - HLA allele
#   - Peptide sequence
#   - 1-log50k(IC50) or raw score
#   - Binding level (SB/<0.5%, WB/<2%, or blank)
#   - %Rank

# Regex to match a data line in NetMHCpan output
# Columns: Pos HLA Peptide [Core Of Gp Gl Ip Il Ic] Identity Score [BindLevel] Rank
_NETMHCPAN_DATA_RE = re.compile(
    r"^\s*(\d+)\s+"              # Position (1-based)
    r"(\S+)\s+"                   # HLA allele
    r"([A-Z]+)\s+"                # Peptide sequence
    r".*?"                         # Optional core/offset columns
    r"(\d+\.\d+)\s+"              # Score (1-log50k(IC50))
    r"(\S+)?\s*"                   # BindLevel (SB, WB, or blank) — optional
    r"(\d+\.\d+)"                  # %Rank
)

# Simpler pattern for parsing: capture position, allele, peptide, score, rank
# The BindLevel column may be missing or combined with rank
_NETMHCPAN_LINE_RE = re.compile(
    r"^\s*(\d+)\s+"              # Position
    r"([\w\*:\-]+)\s+"           # HLA allele (e.g. HLA-A*02:01 or H-2-Db)
    r"([A-Z]{6,11})\s+"          # Peptide (6-11 aa for MHC-I)
)


def parse_netmhcpan_output(output_text: str, allele: str) -> list[dict]:
    """Parse NetMHCpan text output into structured results.

    NetMHCpan output is a space-delimited table with a header section
    followed by data lines.  Each data line contains:
      Position, Allele, Peptide, [core columns], Identity,
      Score (1-log50k IC50), BindLevel, %Rank

    Parameters
    ----------
    output_text : str
        Raw text output from NetMHCpan.
    allele : str
        Expected allele name (for validation).

    Returns
    -------
    list[dict]
        Each dict has keys: position, allele, peptide, score, ic50_nm,
        rank, binding_class.
    """
    results: list[dict] = []
    lines = output_text.strip().splitlines()

    for line in lines:
        # Skip header, comment, and separator lines
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("-"):
            continue
        if stripped.startswith("Pos") or stripped.startswith("  Pos"):
            continue

        # Try to parse as a data line
        parsed = _parse_data_line(line)
        if parsed is not None:
            results.append(parsed)

    return results


def _parse_data_line(line: str) -> dict | None:
    """Parse a single NetMHCpan output data line.

    Returns None if the line cannot be parsed.
    """
    # Split on whitespace and filter empty tokens
    tokens = line.split()
    if len(tokens) < _MIN_DATA_LINE_TOKENS:
        return None

    try:
        # Find the peptide (a sequence of uppercase letters, 6-15 chars)
        peptide_idx = None
        for i, tok in enumerate(tokens):
            if re.match(r"^[A-Z]{6,15}$", tok):
                peptide_idx = i
                break

        if peptide_idx is None:
            return None

        # Position is the first token (1-based)
        position = int(tokens[0])

        # Allele is after position, before peptide
        allele = tokens[1]

        # Peptide
        peptide = tokens[peptide_idx]

        # Look for numeric values after the peptide/identity columns
        # The score and rank are the last two numeric columns
        # BindLevel (SB/WB/<0.5%/etc.) may appear between them
        numeric_tokens = []
        bind_level = None

        # Scan from the right: the last token should be the rank
        # The second-to-last numeric should be the score
        for tok in reversed(tokens):
            try:
                val = float(tok)
                numeric_tokens.append(val)
            except ValueError:
                # Check for bind level markers
                tok_upper = tok.upper()
                if tok_upper in ("SB", "WB", "<0.5%", "<2%", "<0.5", "<2", "<=0.5%", "<=2%"):
                    bind_level = tok_upper

        if len(numeric_tokens) < 2:
            return None

        # Last numeric = rank, second-to-last = score
        rank = numeric_tokens[0]  # reversed, so first = last in original
        score = numeric_tokens[1]

        # Convert NetMHCpan score (1-log50k IC50) to IC50 in nM
        # 1-log50k(IC50) where IC50 is in nM
        # So IC50 = _NETMHCPAN_LOG50K_CONSTANT^(1-score)
        import math
        if score > 0:
            ic50 = _NETMHCPAN_LOG50K_CONSTANT ** (1.0 - score)
        else:
            ic50 = _NETMHCPAN_LOG50K_CONSTANT  # maximum

        # Determine binding class from rank (preferred) or bind_level
        if bind_level in ("SB", "<0.5%", "<0.5", "<=0.5%"):
            binding_class = "strong_binder"
        elif bind_level in ("WB", "<2%", "<2", "<=2%"):
            binding_class = "weak_binder"
        else:
            binding_class = classify_binding_rank(rank)

        return {
            "position": position - 1,  # Convert to 0-based
            "allele": allele,
            "peptide": peptide,
            "score": round(score, 6),
            "ic50_nm": round(ic50, 2),
            "rank": round(rank, 4),
            "binding_class": binding_class,
        }

    except (ValueError, IndexError):
        logger.debug("Failed to parse NetMHCpan data line", exc_info=True)
        return None


# ═══════════════════════════════════════════════════════════════════════════
# Availability check
# ═══════════════════════════════════════════════════════════════════════════

def is_netmhcpan_available(timeout: float = 10.0) -> bool:
    """Check whether NetMHCpan is available (locally installed or web API reachable).

    Returns True if a local NetMHCpan binary is found on PATH, or if the
    DTU Health Tech server responds to a HEAD request.
    Returns False otherwise (no local binary and network error/timeout).
    """
    # Check local installation first (fast, no network)
    if is_netmhcpan_installed():
        return True

    # Fall back to web API check
    try:
        req = urllib.request.Request(
            "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/",
            method="HEAD",
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if resp.status < 500:
                logger.debug("NetMHCpan API reachable (status %d)", resp.status)
                return True
    except Exception as exc:
        logger.debug("NetMHCpan API not reachable: %s", exc)

    return False


# ═══════════════════════════════════════════════════════════════════════════
# NetMHCpan Client
# ═══════════════════════════════════════════════════════════════════════════

class NetMHCpanClient:
    """Client for the NetMHCpan web API.

    Provides methods for MHC-I and MHC-II binding prediction via
    the DTU Health Tech NetMHCpan service.  Includes caching,
    retry logic with exponential backoff, and graceful fallback.

    Usage
    -----
    >>> client = NetMHCpanClient()
    >>> result = client.predict_mhc_i_binding("SIINFEKL", "HLA-A*02:01")
    >>> print(result.binding_class, result.rank)
    weak_binder 1.23

    Parameters
    ----------
    api_url : str
        NetMHCpan API endpoint URL.
    timeout : float
        HTTP request timeout in seconds.
    max_retries : int
        Maximum number of retry attempts per request.
    cache : NetMHCpanCache or None
        Cache instance. If None, uses the module-level default cache.
    use_cache : bool
        Whether to use caching.
    """

    def __init__(
        self,
        api_url: str = DEFAULT_API_URL,
        timeout: float = DEFAULT_TIMEOUT,
        max_retries: int = MAX_RETRIES,
        cache: NetMHCpanCache | None = None,
        use_cache: bool = True,
    ):
        self.api_url = api_url
        self.timeout = timeout
        self.max_retries = max_retries
        self._cache = cache if cache is not None else _get_default_cache()
        self.use_cache = use_cache

    def predict_mhc_i_binding(
        self,
        peptide_sequence: str,
        allele: str,
        epitope_length: int = 9,
    ) -> MHCBindingResult:
        """Predict MHC class I binding for a single peptide.

        Calls the NetMHCpan 4.1 API with the given peptide and allele.

        Parameters
        ----------
        peptide_sequence : str
            Amino acid sequence of the peptide.
        allele : str
            MHC-I allele name (e.g. ``"HLA-A*02:01"``).
        epitope_length : int
            Expected peptide length (default 9).

        Returns
        -------
        MHCBindingResult
            Binding prediction result.

        Raises
        ------
        NetMHCpanError
            If the API call fails after all retries.
        """
        peptide_sequence = peptide_sequence.upper().strip()

        # Check cache
        if self.use_cache:
            cached = self._cache.get(allele, peptide_sequence, epitope_length)
            if cached is not None:
                logger.debug("Cache hit for %s/%s (len=%d)", allele, peptide_sequence, epitope_length)
                return cached

        # Validate input
        if not peptide_sequence:
            raise NetMHCpanError("Peptide sequence must not be empty", allele=allele)
        if len(peptide_sequence) != epitope_length:
            logger.warning(
                "Peptide length %d != epitope_length %d for allele %s",
                len(peptide_sequence), epitope_length, allele,
            )

        # Submit to API
        output_text = self._submit_mhc_i(peptide_sequence, allele, epitope_length)

        # Parse results
        parsed = parse_netmhcpan_output(output_text, allele)

        if not parsed:
            # No binding data returned — treat as non-binder
            result = MHCBindingResult(
                allele=allele,
                peptide=peptide_sequence,
                start_position=0,
                end_position=len(peptide_sequence) - 1,
                binding_score=0.0,
                ic50_nm=_NETMHCPAN_LOG50K_CONSTANT,
                binding_class="non_binder",
                rank=None,
                method="netmhcpan",
            )
        else:
            # Take the first (and usually only) result
            p = parsed[0]
            rank = p.get("rank", _FALLBACK_RANK)
            binding_class = p.get("binding_class", classify_binding_rank(rank if rank is not None else _FALLBACK_RANK))
            binding_score = _rank_to_binding_score(rank) if rank is not None else 0.0

            result = MHCBindingResult(
                allele=allele,
                peptide=peptide_sequence,
                start_position=p.get("position", 0),
                end_position=p.get("position", 0) + len(peptide_sequence) - 1,
                binding_score=binding_score,
                ic50_nm=p.get("ic50_nm"),
                binding_class=binding_class,
                rank=rank,
                method="netmhcpan",
            )

        # Cache result
        if self.use_cache:
            self._cache.put(allele, peptide_sequence, result, epitope_length)

        return result

    def predict_mhc_ii_binding(
        self,
        peptide_sequence: str,
        allele: str,
    ) -> MHCBindingResult:
        """Predict MHC class II binding for a single peptide.

        Calls the NetMHCIIpan 4.0 API with the given peptide and allele.

        Parameters
        ----------
        peptide_sequence : str
            Amino acid sequence of the peptide (15-mer recommended).
        allele : str
            MHC-II allele name (e.g. ``"HLA-DRB1*01:01"``).

        Returns
        -------
        MHCBindingResult
            Binding prediction result.

        Raises
        ------
        NetMHCpanError
            If the API call fails after all retries.
        """
        peptide_sequence = peptide_sequence.upper().strip()

        # Check cache
        if self.use_cache:
            cached = self._cache.get(allele, peptide_sequence, MHC_II_EPITOPE_LENGTH)
            if cached is not None:
                logger.debug("Cache hit for %s/%s (MHC-II)", allele, peptide_sequence)
                return cached

        if not peptide_sequence:
            raise NetMHCpanError("Peptide sequence must not be empty", allele=allele)

        # Submit to MHC-II API
        output_text = self._submit_mhc_ii(peptide_sequence, allele)

        # Parse results
        parsed = parse_netmhcpan_output(output_text, allele)

        if not parsed:
            result = MHCBindingResult(
                allele=allele,
                peptide=peptide_sequence,
                start_position=0,
                end_position=len(peptide_sequence) - 1,
                binding_score=0.0,
                ic50_nm=_NETMHCPAN_LOG50K_CONSTANT,
                binding_class="non_binder",
                rank=None,
                method="netmhcpan",
            )
        else:
            # For MHC-II, take the best-scoring core
            best = min(parsed, key=lambda p: p.get("rank", _FALLBACK_RANK))
            rank = best.get("rank", _FALLBACK_RANK)
            binding_class = best.get("binding_class", classify_binding_rank(rank))
            binding_score = _rank_to_binding_score(rank)

            result = MHCBindingResult(
                allele=allele,
                peptide=peptide_sequence,
                start_position=best.get("position", 0),
                end_position=best.get("position", 0) + len(peptide_sequence) - 1,
                binding_score=binding_score,
                ic50_nm=best.get("ic50_nm"),
                binding_class=binding_class,
                rank=rank,
                method="netmhcpan",
            )

        # Cache result
        if self.use_cache:
            self._cache.put(allele, peptide_sequence, result, MHC_II_EPITOPE_LENGTH)

        return result

    def batch_predict(
        self,
        protein_sequence: str,
        alleles: list[str],
        epitope_lengths: list[int] | None = None,
    ) -> list[MHCBindingResult]:
        """Batch predict MHC binding for a full protein across multiple alleles.

        Extracts all overlapping peptides of the specified lengths from the
        protein sequence and predicts binding for each peptide x allele
        combination.

        Parameters
        ----------
        protein_sequence : str
            Full protein amino acid sequence.
        alleles : list[str]
            MHC alleles to evaluate.
        epitope_lengths : list[int] or None
            Peptide lengths to extract (default [8, 9, 10, 11]).

        Returns
        -------
        list[MHCBindingResult]
            Binding predictions for every peptide x allele combination.
        """
        if epitope_lengths is None:
            epitope_lengths = list(DEFAULT_MHC_I_EPITOPE_LENGTHS)

        protein_sequence = protein_sequence.upper().strip()

        # Validate
        try:
            protein_sequence = validate_protein_sequence(protein_sequence, "NetMHCpan")
        except ValueError as exc:
            raise NetMHCpanError(str(exc)) from exc

        results: list[MHCBindingResult] = []

        for allele in alleles:
            # Determine if this is an MHC-I or MHC-II allele
            is_mhc_ii = _is_mhc_ii_allele(allele)

            if is_mhc_ii:
                # MHC-II: use 15-mer peptides
                peptide_length = MHC_II_EPITOPE_LENGTH
                for start in range(len(protein_sequence) - peptide_length + 1):
                    peptide = protein_sequence[start : start + peptide_length]
                    try:
                        result = self.predict_mhc_ii_binding(peptide, allele)
                        # Adjust positions to protein-relative
                        result.start_position = start
                        result.end_position = start + peptide_length - 1
                        results.append(result)
                    except NetMHCpanError as exc:
                        logger.warning(
                            "NetMHCpan MHC-II prediction failed for %s/%s: %s",
                            allele, peptide, exc,
                        )
                        # Append a non-binder result as fallback
                        results.append(MHCBindingResult(
                            allele=allele,
                            peptide=peptide,
                            start_position=start,
                            end_position=start + peptide_length - 1,
                            binding_score=0.0,
                            ic50_nm=_NETMHCPAN_LOG50K_CONSTANT,
                            binding_class="non_binder",
                            method="netmhcpan_failed",
                        ))
            else:
                # MHC-I: try each epitope length
                for epi_len in epitope_lengths:
                    for start in range(len(protein_sequence) - epi_len + 1):
                        peptide = protein_sequence[start : start + epi_len]
                        try:
                            result = self.predict_mhc_i_binding(peptide, allele, epi_len)
                            # Adjust positions to protein-relative
                            result.start_position = start
                            result.end_position = start + epi_len - 1
                            results.append(result)
                        except NetMHCpanError as exc:
                            logger.warning(
                                "NetMHCpan MHC-I prediction failed for %s/%s: %s",
                                allele, peptide, exc,
                            )
                            results.append(MHCBindingResult(
                                allele=allele,
                                peptide=peptide,
                                start_position=start,
                                end_position=start + epi_len - 1,
                                binding_score=0.0,
                                ic50_nm=_NETMHCPAN_LOG50K_CONSTANT,
                                binding_class="non_binder",
                                method="netmhcpan_failed",
                            ))

        logger.info(
            "batch_predict: %d results for %d alleles, protein length %d",
            len(results), len(alleles), len(protein_sequence),
        )
        return results

    # ────────────────────────────────────────────────────────────
    # Private methods: API submission
    # ────────────────────────────────────────────────────────────

    def _submit_mhc_i(
        self,
        peptide: str,
        allele: str,
        epitope_length: int,
    ) -> str:
        """Submit an MHC-I prediction job to the NetMHCpan API.

        Returns the raw text output from NetMHCpan.
        """
        form_data = urllib.parse.urlencode({
            "configfile": _NETMHCPAN_CONFIG,
            "seq": peptide,
            "allele": allele,
            "len": str(epitope_length),
            "sort": "0",  # Sort by position
        }).encode("utf-8")

        return self._submit_job(form_data)

    def _submit_mhc_ii(
        self,
        peptide: str,
        allele: str,
    ) -> str:
        """Submit an MHC-II prediction job to the NetMHCIIpan API.

        Returns the raw text output from NetMHCIIpan.
        """
        form_data = urllib.parse.urlencode({
            "configfile": _NETMHCIIPAN_CONFIG,
            "seq": peptide,
            "allele": allele,
        }).encode("utf-8")

        return self._submit_job(form_data)

    def _submit_job(self, form_data: bytes) -> str:
        """Submit a job to the NetMHCpan webface2 CGI endpoint.

        Handles the submit → poll → retrieve workflow with retry
        logic and exponential backoff.

        Returns
        -------
        str
            Raw text output from the API.

        Raises
        ------
        NetMHCpanError
            If the API call fails after all retries.
        """
        last_error: str | None = None

        for attempt in range(1, self.max_retries + 1):
            try:
                # Submit the job
                req = urllib.request.Request(
                    self.api_url,
                    data=form_data,
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    method="POST",
                )
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    response_text = resp.read().decode("utf-8")

                # Check if the response contains a redirect/job page
                # The webface2 CGI typically returns either:
                # 1. Direct results (for small inputs)
                # 2. A page with a job ID that we need to poll
                result_url = _extract_result_url(response_text)

                if result_url is not None:
                    # Poll for results
                    output = self._poll_for_results(result_url)
                    return output
                else:
                    # The response might contain the results directly
                    # Check if it looks like NetMHCpan output
                    if _looks_like_netmhcpan_output(response_text):
                        return response_text
                    else:
                        # Try to extract error message
                        error_msg = _extract_error_message(response_text)
                        if error_msg:
                            last_error = f"API returned error: {error_msg}"
                            logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)
                        else:
                            last_error = "API returned unexpected response format"
                            logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)

            except urllib.error.HTTPError as exc:
                if exc.code == 429:
                    # Rate limited
                    retry_after = exc.headers.get("Retry-After")
                    wait = float(retry_after) if retry_after else RETRY_BASE_DELAY * (2 ** (attempt - 1))
                    last_error = f"API rate limited (429), waiting {wait:.1f}s"
                    logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)
                    time.sleep(wait)
                    continue
                elif exc.code >= 500:
                    wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
                    last_error = f"API server error ({exc.code}), retrying in {wait:.1f}s"
                    logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)
                    time.sleep(wait)
                    continue
                elif exc.code >= 400:
                    # Client error — do not retry
                    error_body = exc.read().decode("utf-8", errors="replace")[:_ERROR_BODY_MAX_CHARS]
                    last_error = f"API client error ({exc.code}): {error_body}"
                    logger.error("API client error, not retrying: %s", last_error)
                    break
                else:
                    last_error = f"HTTP error {exc.code}"
                    logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)
                    continue

            except urllib.error.URLError as exc:
                wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
                last_error = f"Network error: {exc.reason}"
                logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)
                time.sleep(wait)
                continue

            except TimeoutError:
                wait = RETRY_BASE_DELAY * (2 ** (attempt - 1))
                last_error = f"Request timed out after {self.timeout}s"
                logger.warning("Attempt %d/%d: %s", attempt, self.max_retries, last_error)
                time.sleep(wait)
                continue

            except Exception as exc:
                last_error = f"Unexpected error: {exc}"
                logger.error("Unexpected NetMHCpan API error: %s", exc)
                break

        # All retries exhausted
        raise NetMHCpanError(
            f"API call failed after {self.max_retries} attempts. Last error: {last_error}"
        )

    def _poll_for_results(self, result_url: str) -> str:
        """Poll a NetMHCpan job URL until results are available.

        Parameters
        ----------
        result_url : str
            URL to poll for completed results.

        Returns
        -------
        str
            Raw text output from NetMHCpan.

        Raises
        ------
        NetMHCpanError
            If polling exceeds MAX_POLL_ATTEMPTS.
        """
        for poll_attempt in range(MAX_POLL_ATTEMPTS):
            try:
                req = urllib.request.Request(result_url, method="GET")
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    response_text = resp.read().decode("utf-8")

                # Check if results are ready
                if _looks_like_netmhcpan_output(response_text):
                    return response_text

                # Check if still processing
                if "still processing" in response_text.lower() or "job is running" in response_text.lower():
                    time.sleep(POLL_INTERVAL)
                    continue

                # Check for errors
                error_msg = _extract_error_message(response_text)
                if error_msg:
                    raise NetMHCpanError(f"Job failed: {error_msg}")

                # If we can't determine the state, wait and retry
                time.sleep(POLL_INTERVAL)

            except urllib.error.URLError as exc:
                logger.warning("Poll attempt %d failed: %s", poll_attempt, exc)
                time.sleep(POLL_INTERVAL)
                continue

        raise NetMHCpanError(
            f"Polling for results timed out after {MAX_POLL_ATTEMPTS * POLL_INTERVAL:.0f}s"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════

def _is_mhc_ii_allele(allele: str) -> bool:
    """Determine if an allele is MHC class II based on naming convention."""
    allele_upper = allele.upper()
    # MHC-II alleles: HLA-DR, HLA-DQ, HLA-DP, H2-I, H2-Ab, H2-Aa
    return any(
        allele_upper.startswith(prefix)
        for prefix in ("HLA-DR", "HLA-DQ", "HLA-DP", "H2-I", "H2-A")
    )


def _extract_result_url(html: str) -> str | None:
    """Extract the result URL from a webface2 job submission response.

    The webface2 CGI returns HTML with a meta-refresh or a link to the
    results page.
    """
    # Look for meta refresh: <meta http-equiv="refresh" content="0;url=...">
    meta_match = re.search(
        r'<meta\s+http-equiv=["\']refresh["\']\s+content=["\']\d+;\s*url=([^"\']+)["\']',
        html, re.IGNORECASE,
    )
    if meta_match:
        return meta_match.group(1)

    # Look for a link to the result page
    link_match = re.search(
        r'href=["\']([^"\']*(?:result|output)[^"\']*)["\']',
        html, re.IGNORECASE,
    )
    if link_match:
        return link_match.group(1)

    # Look for a job ID that we can construct a URL from
    job_match = re.search(r'jobid[=:]\s*["\']?(\w+)', html, re.IGNORECASE)
    if job_match:
        job_id = job_match.group(1)
        return f"https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid={job_id}"

    return None


def _looks_like_netmhcpan_output(text: str) -> bool:
    """Check if text looks like NetMHCpan output."""
    # NetMHCpan output contains characteristic header lines
    markers = [
        "NetMHCpan",
        "NetMHCIIpan",
        "#  Pos",
        "1-log50k",
        "%Rank",
        "BindLevel",
    ]
    return sum(1 for m in markers if m in text) >= 2


def _extract_error_message(html: str) -> str | None:
    """Extract an error message from an HTML response."""
    # Look for error messages in HTML
    error_match = re.search(
        r'(?:error|Error|ERROR)[:\s]+([^<\n]{10,200})',
        html,
    )
    if error_match:
        return error_match.group(1).strip()

    # Look for common error indicators
    if "unrecognized allele" in html.lower():
        return "Unrecognized allele"
    if "invalid sequence" in html.lower():
        return "Invalid sequence"

    return None


# ═══════════════════════════════════════════════════════════════════════════
# Convenience module-level functions
# ═══════════════════════════════════════════════════════════════════════════

_default_client: NetMHCpanClient | None = None


def _get_default_client() -> NetMHCpanClient:
    """Return the module-level default client (lazy-initialized)."""
    global _default_client
    if _default_client is None:
        _default_client = NetMHCpanClient()
    return _default_client


def predict_mhc_i_binding(
    peptide_sequence: str,
    allele: str,
    epitope_length: int = 9,
) -> MHCBindingResult:
    """Predict MHC-I binding using the module-level default client.

    Convenience wrapper around :meth:`NetMHCpanClient.predict_mhc_i_binding`.
    """
    return _get_default_client().predict_mhc_i_binding(peptide_sequence, allele, epitope_length)


def predict_mhc_ii_binding(
    peptide_sequence: str,
    allele: str,
) -> MHCBindingResult:
    """Predict MHC-II binding using the module-level default client.

    Convenience wrapper around :meth:`NetMHCpanClient.predict_mhc_ii_binding`.
    """
    return _get_default_client().predict_mhc_ii_binding(peptide_sequence, allele)


def batch_predict(
    protein_sequence: str,
    alleles: list[str],
    epitope_lengths: list[int] | None = None,
) -> list[MHCBindingResult]:
    """Batch predict MHC binding using the module-level default client.

    Convenience wrapper around :meth:`NetMHCpanClient.batch_predict`.
    """
    return _get_default_client().batch_predict(protein_sequence, alleles, epitope_lengths)


# ═══════════════════════════════════════════════════════════════════════════
# Fallback chain adapter functions
#
# These provide a clean interface for the MHC prediction fallback chain.
# They return None when NetMHCpan is not available, allowing the caller
# to try the next fallback without raising exceptions.
# ═══════════════════════════════════════════════════════════════════════════

# Names of NetMHCpan binaries to search for on PATH
_NETMHCPAN_BINARIES: list[str] = ["netMHCpan", "NetMHCpan", "netmhcpan"]

# Cached result for is_netmhcpan_installed
_netmhcpan_installed_cache: bool | None = None


def is_netmhcpan_installed() -> bool:
    """Check whether a NetMHCpan binary is installed locally.

    Searches the system PATH for any of the known NetMHCpan binary names.
    The result is cached after the first call.

    Returns
    -------
    bool
        True if a NetMHCpan binary is found on PATH.
    """
    global _netmhcpan_installed_cache
    if _netmhcpan_installed_cache is not None:
        return _netmhcpan_installed_cache
    for binary in _NETMHCPAN_BINARIES:
        if shutil.which(binary) is not None:
            _netmhcpan_installed_cache = True
            return True
    _netmhcpan_installed_cache = False
    return False


def _predict_binding_local(
    allele: str,
    peptide: str,
    epitope_length: int = 9,
    timeout: float = 120.0,
) -> MHCBindingResult:
    """Run a local NetMHCpan binary to predict MHC binding.

    Parameters
    ----------
    allele : str
        MHC allele name.
    peptide : str
        Peptide sequence (must not be empty).
    epitope_length : int
        Expected peptide length (default 9).
    timeout : float
        Maximum runtime in seconds (default 120).

    Returns
    -------
    MHCBindingResult
        Binding prediction result.

    Raises
    ------
    NetMHCpanError
        If the binary is not found, the peptide is empty, or
        the binary returns a non-zero exit code.
    """
    if not peptide:
        raise NetMHCpanError("Peptide must not be empty", allele=allele)

    # Find binary
    binary_path = None
    for binary in _NETMHCPAN_BINARIES:
        path = shutil.which(binary)
        if path is not None:
            binary_path = path
            break

    if binary_path is None:
        raise NetMHCpanError(
            "NetMHCpan binary not found on PATH", allele=allele
        )

    # Run the binary
    try:
        result = subprocess.run(
            [binary_path, "-a", allele, "-p", peptide, "-l", str(epitope_length)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise NetMHCpanError(
            f"NetMHCpan timed out after {timeout}s", allele=allele
        ) from exc

    if result.returncode != 0:
        raise NetMHCpanError(
            f"NetMHCpan exited with code {result.returncode}", allele=allele
        )

    # Parse the output
    parsed = parse_netmhcpan_output(result.stdout, allele)
    if parsed:
        p = parsed[0]
        rank = p.get("rank", _FALLBACK_RANK)
        binding_class = p.get("binding_class", classify_binding_rank(rank))
        binding_score = _rank_to_binding_score(rank)
        return MHCBindingResult(
            allele=allele,
            peptide=peptide,
            start_position=p.get("position", 0),
            end_position=p.get("position", 0) + len(peptide) - 1,
            binding_score=binding_score,
            ic50_nm=p.get("ic50_nm"),
            binding_class=binding_class,
            rank=rank,
            method="netmhcpan_local",
        )

    # No parseable data — return non-binder
    return MHCBindingResult(
        allele=allele,
        peptide=peptide,
        start_position=0,
        end_position=len(peptide) - 1,
        binding_score=0.0,
        binding_class="non_binder",
        method="netmhcpan_local",
    )


def predict_binding_netmhcpan(
    allele: str,
    peptide: str,
    epitope_length: int = 9,
    timeout: float = 10.0,
) -> MHCBindingResult | None:
    """Adapter for the MHC prediction fallback chain — single peptide.

    Returns None when NetMHCpan is not available (neither locally
    installed nor web API reachable), allowing the caller to try
    the next fallback predictor without catching exceptions.

    Parameters
    ----------
    allele : str
        MHC allele name.
    peptide : str
        Peptide amino acid sequence.
    epitope_length : int
        Expected peptide length (default 9).
    timeout : float
        Timeout for availability check (default 10s).

    Returns
    -------
    MHCBindingResult or None
        Binding prediction result, or None if NetMHCpan is unavailable.
    """
    if not is_netmhcpan_available(timeout=timeout):
        return None

    # Try local binary first, then fall back to web API
    if is_netmhcpan_installed():
        try:
            return _predict_binding_local(allele, peptide, epitope_length)
        except Exception:
            return None

    # Use web API
    try:
        return predict_mhc_i_binding(peptide, allele, epitope_length)
    except Exception:
        return None


def batch_predict_binding_netmhcpan(
    protein_sequence: str,
    alleles: list[str],
    epitope_lengths: list[int] | None = None,
    timeout: float = 10.0,
) -> list[MHCBindingResult] | None:
    """Adapter for the MHC prediction fallback chain — batch.

    Returns None when NetMHCpan is not available, allowing the
    caller to try the next fallback predictor.

    Parameters
    ----------
    protein_sequence : str
        Full protein amino acid sequence.
    alleles : list[str]
        MHC alleles to evaluate.
    epitope_lengths : list[int] or None
        Peptide lengths to extract (default [8, 9, 10, 11]).
    timeout : float
        Timeout for availability check (default 10s).

    Returns
    -------
    list[MHCBindingResult] or None
        Batch binding predictions, or None if NetMHCpan is unavailable.
    """
    if not is_netmhcpan_available(timeout=timeout):
        return None

    try:
        return batch_predict(protein_sequence, alleles, epitope_lengths)
    except Exception:
        return None
