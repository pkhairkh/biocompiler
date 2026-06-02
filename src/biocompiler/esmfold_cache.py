"""
BioCompiler ESMFold Cache
=========================

In-memory and file-based caching for ESMFold structure predictions.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass, field
from typing import Dict, Optional

from .esmfold import ESMFoldResult


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
                        mean_plddt=data["mean_plddt"],
                        ptdmt=data.get("ptdmt", 0.0),
                        num_residues=data.get("num_residues", len(data["protein"])),
                        model_name=data.get("model_name", "esmfold-offline"),
                        cached=True,
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
                    "mean_plddt": result.mean_plddt,
                    "model_name": result.model_name,
                }
                # Add optional fields if they exist
                if hasattr(result, 'plddt_scores'):
                    data["plddt_scores"] = result.plddt_scores
                if hasattr(result, 'execution_time_s'):
                    data["execution_time_s"] = result.execution_time_s
                if hasattr(result, 'success'):
                    data["success"] = result.success
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
