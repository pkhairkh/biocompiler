#!/usr/bin/env python3
"""Helper script to run mutation testing with mutmut.

BioCompiler mutation testing runner — wraps mutmut with sensible defaults
and provides per-module and full-suite modes.

Usage:
    # Run mutation tests on all configured modules
    python scripts/run_mutation_tests.py

    # Run on a single module only
    python scripts/run_mutation_tests.py --module greedy

    # Run with a custom timeout per mutant (seconds)
    python scripts/run_mutation_tests.py --timeout 30

    # Show results from the last run (no re-run)
    python scripts/run_mutation_tests.py --results

    # Clean up mutation testing artifacts
    python scripts/run_mutation_tests.py --clean

Requirements:
    pip install mutmut
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Project root (one level up from scripts/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Modules available for targeted mutation testing
MODULE_ALIASES: dict[str, str] = {
    "greedy": "src/biocompiler/optimizer/greedy.py",
    "pipeline": "src/biocompiler/optimizer/pipeline.py",
    "biosecurity": "src/biocompiler/biosecurity.py",
}

DEFAULT_TEST_COMMAND = (
    "python -m pytest tests/ -x -q -m 'not slow' --timeout=60"
)

ARTIFACTS = [".mutmut-cache", ".mutmut-coverage"]


def _run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Run a command in the project root."""
    print(f"[mutmut] $ {' '.join(cmd)}")
    return subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=check)


def run_full(timeout: int | None = None) -> None:
    """Run mutation testing on all paths configured in pyproject.toml."""
    cmd = [sys.executable, "-m", "mutmut", "run"]
    if timeout:
        cmd += ["--timeout", str(timeout)]
    _run(cmd)


def run_module(module: str, timeout: int | None = None) -> None:
    """Run mutation testing on a single module by alias."""
    path = MODULE_ALIASES.get(module)
    if path is None:
        available = ", ".join(sorted(MODULE_ALIASES))
        print(f"Unknown module '{module}'. Available: {available}", file=sys.stderr)
        sys.exit(1)

    cmd = [
        sys.executable,
        "-m",
        "mutmut",
        "run",
        "--paths-to-mutate",
        path,
    ]
    if timeout:
        cmd += ["--timeout", str(timeout)]
    _run(cmd)


def show_results() -> None:
    """Display results from the last mutation testing run."""
    _run([sys.executable, "-m", "mutmut", "results"])


def show_html() -> None:
    """Generate and open HTML report."""
    _run([sys.executable, "-m", "mutmut", "html"])
    html_path = PROJECT_ROOT / "mutmut-report.html"
    if html_path.exists():
        print(f"[mutmut] HTML report generated: {html_path}")


def clean() -> None:
    """Remove mutation testing cache and coverage artifacts."""
    for artifact in ARTIFACTS:
        p = PROJECT_ROOT / artifact
        if p.is_dir():
            shutil.rmtree(p)
            print(f"[mutmut] Removed directory: {p}")
        elif p.exists():
            p.unlink()
            print(f"[mutmut] Removed file: {p}")
    # Also remove generated HTML report
    html_path = PROJECT_ROOT / "mutmut-report.html"
    if html_path.exists():
        html_path.unlink()
        print(f"[mutmut] Removed: {html_path}")
    print("[mutmut] Clean complete.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BioCompiler mutation testing helper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--module",
        choices=sorted(MODULE_ALIASES),
        help="Run mutation tests on a single module instead of all configured paths",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Timeout per mutant in seconds (default: mutmut built-in)",
    )
    parser.add_argument(
        "--results",
        action="store_true",
        help="Show results from the last run without re-running",
    )
    parser.add_argument(
        "--html",
        action="store_true",
        help="Generate an HTML report of the last run",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove mutation testing cache and artifacts",
    )

    args = parser.parse_args()

    if args.clean:
        clean()
        return

    # Ensure mutmut is installed
    try:
        subprocess.run(
            [sys.executable, "-m", "mutmut", "--version"],
            capture_output=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(
            "[mutmut] mutmut is not installed. Install with: pip install mutmut",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.results:
        show_results()
        return

    if args.html:
        show_html()
        return

    if args.module:
        run_module(args.module, timeout=args.timeout)
    else:
        run_full(timeout=args.timeout)

    # Show a summary after the run
    show_results()


if __name__ == "__main__":
    main()
