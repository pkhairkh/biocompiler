#!/usr/bin/env python3
"""Run the organism-aware CAI recovery benchmark and print results.

Demonstrates that disabling eukaryotic constraints (splice-site avoidance,
CpG-island avoidance) for prokaryotic targets such as E. coli recovers
the ~0.27 CAI gap caused by applying irrelevant constraints.

Usage::

    python -m biocompiler.scripts.run_organism_aware_benchmark
    python scripts/run_organism_aware_benchmark.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the organism-aware CAI recovery benchmark",
    )
    parser.add_argument(
        "--organism",
        default="Escherichia_coli",
        help="Target organism (default: Escherichia_coli)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Output results as JSON instead of formatted text",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    from biocompiler.benchmarking.organism_aware_benchmark import (
        benchmark_organism_aware_cai,
        print_organism_aware_report,
    )

    results = benchmark_organism_aware_cai(organism=args.organism)

    if args.json_output:
        print(json.dumps(results, indent=2))
    else:
        print_organism_aware_report(results)

    # Exit with a hint about the recovery
    recovery = results.get("mean_cai_recovery", 0.0)
    if recovery > 0.01:
        print(f"\n[PASS] CAI recovery of +{recovery:.4f} confirmed for {args.organism}")
    else:
        print(f"\n[FAIL] Negligible CAI recovery ({recovery:+.4f}) for {args.organism}")

    sys.exit(0)


if __name__ == "__main__":
    main()
