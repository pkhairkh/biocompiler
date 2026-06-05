#!/usr/bin/env python3
"""Run the BioCompiler vs DNAchisel head-to-head benchmark.

Usage::

    python scripts/run_head_to_head.py [--organism Escherichia_coli] [--genes lacZ,recA,groEL]

This script runs the head-to-head benchmark comparing BioCompiler (with
organism-aware constraints) against DNAchisel, measuring CAI, speed,
and constraint satisfaction across the E. coli gene panel.
"""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    """Entry point for the head-to-head benchmark script."""
    parser = argparse.ArgumentParser(
        description="BioCompiler vs DNAchisel head-to-head CAI optimizer benchmark",
    )
    parser.add_argument(
        "--organism",
        default="Escherichia_coli",
        help="Target organism for codon optimization (default: Escherichia_coli)",
    )
    parser.add_argument(
        "--genes",
        default=None,
        help="Comma-separated list of gene names to benchmark (default: all E. coli genes)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Configure logging
    import logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # Parse gene list
    genes: list[str] | None = None
    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]

    # Run benchmark
    from biocompiler.benchmarking.head_to_head_benchmark import (
        run_head_to_head,
        print_head_to_head_report,
    )

    print(f"Running head-to-head benchmark for {args.organism}...")
    result = run_head_to_head(genes=genes, organism=args.organism)

    # Print report
    print_head_to_head_report(result)

    # Exit code based on whether benchmark actually ran
    if result.num_genes == 0:
        print("ERROR: No genes were benchmarked.", file=sys.stderr)
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
