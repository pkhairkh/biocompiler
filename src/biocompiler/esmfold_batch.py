"""Batch ESMFold structure prediction with progress and error isolation.

Provides concurrent batch prediction of protein structures using ESMFold,
with per-item error isolation, caching, rate limiting, and progress tracking.

Example usage::

    from biocompiler.esmfold_batch import predict_proteins, format_batch_report

    proteins = ["MVHLTPEEKSAVTALWGKVNV", "MVLSPADKTNVKAAWGKVGA"]
    result = predict_proteins(proteins, names=["HBB_frag", "Alpha_frag"])
    print(format_batch_report(result))
"""

from __future__ import annotations

import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError
from dataclasses import dataclass, field
from threading import Semaphore
from typing import Callable

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STANDARD_AMINO_ACIDS: set[str] = {
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
}

MAX_BATCH_SIZE = 50
MAX_PROTEIN_LENGTH = 1000


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Time estimation
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Batch prediction
# ---------------------------------------------------------------------------

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
    start = time.monotonic()

    try:
        semaphore.acquire()
        try:
            # Attempt cache lookup first.
            if use_cache:
                try:
                    from biocompiler.esmfold_cache import cached_predict

                    cached_result = cached_predict(protein)
                    if cached_result is not None:
                        elapsed = round(time.monotonic() - start, 3)
                        result = {
                            "name": name,
                            "status": "success",
                            "from_cache": True,
                            "mean_plddt": getattr(cached_result, "mean_plddt", None),
                            "length": len(protein),
                            "pdb": getattr(cached_result, "pdb", None),
                            "time_s": elapsed,
                        }
                        return result
                except ImportError:
                    logger.debug("esmfold_cache not available, skipping cache")
                except Exception as exc:
                    logger.debug("Cache lookup failed for %s: %s", name, exc)

            # Call the prediction API.
            from biocompiler.esmfold import predict_structure, ESMFoldResult

            prediction: ESMFoldResult = predict_structure(protein)
            elapsed = round(time.monotonic() - start, 3)

            result = {
                "name": name,
                "status": "success",
                "from_cache": False,
                "mean_plddt": getattr(prediction, "mean_plddt", None),
                "length": len(protein),
                "pdb": getattr(prediction, "pdb", None),
                "time_s": elapsed,
            }

        finally:
            semaphore.release()

    except TimeoutError:
        result["error"] = f"Prediction timed out"
        result["time_s"] = round(time.monotonic() - start, 3)
    except ImportError as exc:
        result["error"] = f"ESMFold module not available: {exc}"
        result["time_s"] = round(time.monotonic() - start, 3)
    except Exception as exc:
        result["error"] = str(exc)
        result["time_s"] = round(time.monotonic() - start, 3)

    return result


def predict_batch(
    request: BatchStructureRequest,
    progress_callback: Callable[[int, int, dict], None] | None = None,
) -> BatchStructureResult:
    """Run batch ESMFold structure prediction with concurrency and error isolation.

    Each protein is predicted independently — one failure does not affect
    others (unless ``request.stop_on_failure`` is *True*).

    Args:
        request: The batch request specification.
        progress_callback: Optional callable ``(completed, total, latest_result)``
            invoked after each protein finishes.

    Returns:
        Aggregated :class:`BatchStructureResult`.
    """
    t0 = time.monotonic()

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

    # Prepare concurrency primitives.
    semaphore = Semaphore(request.max_concurrent)
    total = len(request.proteins)
    completed = 0
    results: list[dict] = [None] * total  # type: ignore[list-item]
    cancel = False

    # Submit all tasks.
    with ThreadPoolExecutor(max_workers=request.max_concurrent) as executor:
        future_to_index: dict = {}
        for idx, (protein, name) in enumerate(zip(request.proteins, names)):
            future = executor.submit(
                _predict_single,
                protein,
                name,
                request.use_cache,
                semaphore,
            )
            future_to_index[future] = idx

        # Collect results as they complete.
        for future in as_completed(future_to_index, timeout=None):
            idx = future_to_index[future]
            try:
                result = future.result(timeout=request.timeout_per_protein)
            except TimeoutError:
                result = {
                    "name": names[idx],
                    "status": "error",
                    "error": "Prediction timed out",
                    "time_s": request.timeout_per_protein,
                }
            except Exception as exc:
                result = {
                    "name": names[idx],
                    "status": "error",
                    "error": str(exc),
                    "time_s": 0.0,
                }

            results[idx] = result
            completed += 1

            # Progress callback.
            if progress_callback is not None:
                try:
                    progress_callback(completed, total, result)
                except Exception as cb_exc:
                    logger.warning("Progress callback raised: %s", cb_exc)

            # Stop-on-failure: cancel remaining work.
            if request.stop_on_failure and result["status"] == "error":
                cancel = True
                for pending_future in future_to_index:
                    pending_future.cancel()
                # Fill remaining slots with cancelled results.
                for j in range(total):
                    if results[j] is None:
                        results[j] = {
                            "name": names[j],
                            "status": "error",
                            "error": "Cancelled due to stop_on_failure",
                            "time_s": 0.0,
                        }
                break

    # If not cancelled, fill any remaining None slots (shouldn't happen normally).
    for j in range(total):
        if results[j] is None:
            results[j] = {
                "name": names[j],
                "status": "error",
                "error": "Result not collected (unknown reason)",
                "time_s": 0.0,
            }

    total_time = round(time.monotonic() - t0, 3)

    # Aggregate statistics.
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


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

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
