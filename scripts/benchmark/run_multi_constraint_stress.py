#!/usr/bin/env python3
"""
Run the Multi-Constraint Stress Test
=====================================

Executes all 4 stress-test scenarios (restriction sites, no CpG,
everything-conflicts, no cryptic splice) on HBB (human β-globin)
and GFP, then prints the provenance-driven CAI tradeoff report.

Usage::

    python scripts/run_multi_constraint_stress.py
"""

from __future__ import annotations

import logging
import sys

# Configure logging to show progress
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)


def main() -> None:
    """Entry point for the multi-constraint stress test runner."""
    from biocompiler.benchmarking.multi_constraint_stress import (
        run_all_stress_tests,
        print_stress_test_report,
    )

    print("Starting multi-constraint stress tests...")
    print("This may take a minute as each scenario runs full optimisation.")
    print()

    results = run_all_stress_tests()
    print_stress_test_report(results)

    # Print a brief summary to stderr for CI pipelines
    n_ok = sum(1 for r in results if r.cai_constrained > 0)
    n_total = len(results)
    print(
        f"\nCompleted {n_ok}/{n_total} scenarios successfully.",
        file=sys.stderr,
    )

    if n_ok < n_total:
        print("WARNING: Some scenarios did not produce valid results.", file=sys.stderr)
        sys.exit(1)
    else:
        print("All scenarios completed successfully.", file=sys.stderr)
        sys.exit(0)


if __name__ == "__main__":
    main()
