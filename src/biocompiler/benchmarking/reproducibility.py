"""
Benchmark reproducibility utilities.
Ensures benchmarks are deterministic, versioned, and reproducible.

Provides:
  - ``capture_environment``: Record the full execution environment.
  - ``save_benchmark_results``: Persist results with environment metadata.
  - ``load_benchmark_results``: Load previously saved results.
  - ``validate_benchmark_reproducibility``: Compare new results against baseline.

Usage::

    from biocompiler.benchmarking.reproducibility import (
        capture_environment,
        save_benchmark_results,
        load_benchmark_results,
        validate_benchmark_reproducibility,
    )

    env = capture_environment()

    # After running benchmarks:
    save_benchmark_results(results, "benchmark_results/v12_baseline.json")

    # Later, validate reproducibility:
    ok = validate_benchmark_reproducibility(
        "benchmark_results/v12_baseline.json",
        new_results=current_results,
        tolerance=0.01,
    )
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "capture_environment",
    "deterministic_seed",
    "validate_benchmark_reproducibility",
    "save_benchmark_results",
    "load_benchmark_results",
]


# ---------------------------------------------------------------------------
# Environment capture
# ---------------------------------------------------------------------------


def capture_environment() -> dict[str, Any]:
    """Capture the full execution environment for reproducibility.

    Records Python version, installed package versions, OS details,
    CPU information, and git commit hash so that benchmark runs
    can be fully reproduced.

    Returns:
        Dictionary with the following keys:

        - ``python_version``: Full Python version string.
        - ``python_executable``: Path to the Python interpreter.
        - ``platform_system``: OS name (e.g., "Linux").
        - ``platform_release``: OS release version.
        - ``platform_machine``: CPU architecture (e.g., "x86_64").
        - ``platform_processor``: Processor identifier.
        - ``cpu_count``: Number of logical CPUs.
        - ``git_commit``: Short commit hash, or ``"unknown"``.
        - ``git_branch``: Current git branch, or ``"unknown"``.
        - ``git_dirty``: Whether the working tree has uncommitted changes.
        - ``package_versions``: Dict of installed package → version.
        - ``biocompiler_version``: BioCompiler version string.
        - ``timestamp``: ISO 8601 UTC timestamp.
        - ``hostname``: Machine hostname.
        - ``env_vars``: Selected environment variables affecting benchmarks.
    """
    env: dict[str, Any] = {}

    # Python details
    env["python_version"] = sys.version
    env["python_executable"] = sys.executable

    # Platform details
    env["platform_system"] = platform.system()
    env["platform_release"] = platform.release()
    env["platform_machine"] = platform.machine()
    env["platform_processor"] = platform.processor() or "unknown"
    env["cpu_count"] = os.cpu_count() or 0

    # Git details
    git_info = _capture_git_info()
    env["git_commit"] = git_info["commit"]
    env["git_branch"] = git_info["branch"]
    env["git_dirty"] = git_info["dirty"]

    # Package versions
    env["package_versions"] = _capture_package_versions()

    # BioCompiler version
    env["biocompiler_version"] = _capture_biocompiler_version()

    # Timestamp
    env["timestamp"] = datetime.now(timezone.utc).isoformat()

    # Hostname
    env["hostname"] = platform.node()

    # Selected environment variables
    env["env_vars"] = {
        k: os.environ.get(k, "")
        for k in [
            "BIOCOMPILER_SEED",
            "BIOCOMPILER_NUMBA_DISABLE_JIT",
            "NUMBA_DISABLE_JIT",
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
        ]
        if os.environ.get(k) is not None
    }

    return env


def _capture_git_info() -> dict[str, Any]:
    """Capture git commit, branch, and dirty state.

    Returns:
        Dict with keys ``commit``, ``branch``, ``dirty``.
    """
    result: dict[str, Any] = {
        "commit": "unknown",
        "branch": "unknown",
        "dirty": True,
    }

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        result["commit"] = commit
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    try:
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        result["branch"] = branch
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    try:
        subprocess.check_call(
            ["git", "diff", "--quiet"],
            stderr=subprocess.DEVNULL,
        )
        result["dirty"] = False
    except (subprocess.CalledProcessError, FileNotFoundError):
        result["dirty"] = True

    return result


def _capture_package_versions() -> dict[str, str]:
    """Capture versions of key packages.

    Returns:
        Dict mapping package name to version string.
    """
    key_packages = [
        "biocompiler",
        "dnachisel",
        "numpy",
        "scipy",
        "matplotlib",
        "biopython",
        "numba",
    ]
    versions: dict[str, str] = {}

    for pkg in key_packages:
        try:
            import importlib
            mod = importlib.import_module(pkg)
            ver = getattr(mod, "__version__", "installed")
            versions[pkg] = ver
        except ImportError:
            versions[pkg] = "not_installed"

    return versions


def _capture_biocompiler_version() -> str:
    """Capture the BioCompiler version string.

    Returns:
        Version string, or ``"unknown"`` if unavailable.
    """
    try:
        from .. import __version__
        return __version__
    except (ImportError, AttributeError):
        return "unknown"


# ---------------------------------------------------------------------------
# Deterministic seed
# ---------------------------------------------------------------------------


def deterministic_seed() -> int:
    """Generate a deterministic seed based on the current commit hash.

    The seed is derived from the SHA-1 of the current git commit,
    ensuring that the same code always produces the same seed.
    If git is unavailable, falls back to a fixed seed.

    Returns:
        Integer seed suitable for ``random.seed()`` or ``numpy.random``.
    """
    git_info = _capture_git_info()
    commit = git_info["commit"]

    if commit == "unknown":
        # Fixed fallback seed
        return 42

    # Hash the commit to get a stable integer
    digest = hashlib.sha256(commit.encode()).hexdigest()
    return int(digest[:8], 16)


# ---------------------------------------------------------------------------
# Results persistence
# ---------------------------------------------------------------------------


def save_benchmark_results(results: dict[str, Any], output_path: str) -> None:
    """Save benchmark results to a JSON file with environment metadata.

    Args:
        results: Benchmark results dictionary.
        output_path: File path to save results.
    """
    output: dict[str, Any] = {"results": results}
    output["environment"] = capture_environment()

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, default=str)

    logger.info("Benchmark results saved to %s", path)


def load_benchmark_results(path: str) -> dict[str, Any]:
    """Load benchmark results from a JSON file.

    Args:
        path: File path to load results from.

    Returns:
        Dictionary with ``results`` and ``environment`` keys.
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    logger.info("Benchmark results loaded from %s", path)
    return data


# ---------------------------------------------------------------------------
# Reproducibility validation
# ---------------------------------------------------------------------------


def validate_benchmark_reproducibility(
    results_path: str,
    new_results: dict[str, Any],
    tolerance: float = 0.01,
) -> bool:
    """Compare new benchmark results against a saved baseline for reproducibility.

    Checks that key numerical metrics are within ``tolerance`` of the
    baseline values. This is essential for verifying that benchmark
    results are stable across runs, machines, or Python versions.

    Args:
        results_path: Path to the saved baseline results JSON file.
        new_results: Current results to compare against the baseline.
        tolerance: Maximum allowed relative difference (default 1%).

    Returns:
        ``True`` if all metrics are within tolerance, ``False`` otherwise.
    """
    baseline = load_benchmark_results(results_path)
    baseline_results = baseline.get("results", {})

    metric_keys = [
        "mean_cai_biocompiler",
        "mean_cai_dnachisel",
        "mean_speed_biocompiler_ms",
        "mean_speed_dnachisel_ms",
        "constraint_satisfaction_rate_biocompiler",
        "constraint_satisfaction_rate_dnachisel",
    ]

    # Also check summary sub-dict if present
    if "summary" in baseline_results:
        for key in list(metric_keys):
            metric_keys.append(f"summary.{key}")

    all_within_tolerance = True

    for key in metric_keys:
        baseline_val = _extract_nested(baseline_results, key)
        current_val = _extract_nested(new_results, key)

        if baseline_val is None or current_val is None:
            continue

        if baseline_val == 0 and current_val == 0:
            continue

        if baseline_val == 0:
            logger.warning(
                "Metric '%s': baseline is 0 but current is %s", key, current_val
            )
            all_within_tolerance = False
            continue

        rel_diff = abs(current_val - baseline_val) / abs(baseline_val)

        if rel_diff > tolerance:
            logger.warning(
                "Metric '%s' exceeds tolerance: "
                "baseline=%s, current=%s, rel_diff=%.4f > %.4f",
                key, baseline_val, current_val, rel_diff, tolerance,
            )
            all_within_tolerance = False
        else:
            logger.info(
                "Metric '%s' within tolerance: "
                "baseline=%s, current=%s, rel_diff=%.4f <= %.4f",
                key, baseline_val, current_val, rel_diff, tolerance,
            )

    return all_within_tolerance


def _extract_nested(data: dict, key: str) -> float | None:
    """Extract a value from a possibly nested dict using dot notation.

    Args:
        data: Dictionary to search.
        key: Key, possibly with dots for nested access (e.g., "cai.mean").

    Returns:
        Float value, or ``None`` if not found.
    """
    parts = key.split(".")
    current: Any = data

    for part in parts:
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            return None

    try:
        return float(current)
    except (TypeError, ValueError):
        return None
