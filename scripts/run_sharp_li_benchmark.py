#!/usr/bin/env python3
"""
Run the Sharp-Li vs Kazusa CAI benchmark.

Compares CAI values computed with two different reference sets against
published values from Sharp & Li (1987) and Puigbo et al. (2008).

Usage::

    python -m biocompiler.scripts.run_sharp_li_benchmark

References
----------
1. Sharp, P.M. & Li, W.-H. (1987). Nucleic Acids Research 15:1281-1295.
2. Puigbo, P., Bravo, I.G. & Garcia-Vallve, S. (2008). BMC Bioinformatics 9:65.
"""

from __future__ import annotations

import sys
import os

# Add project src to path so the script can be run directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from biocompiler.benchmarking.sharp_li_benchmark import (
    benchmark_sharp_li_cai,
    print_benchmark_report,
)


def main() -> int:
    """Run the Sharp-Li vs Kazusa CAI benchmark and print the report.

    Returns:
        Exit code: 0 if Sharp-Li is closer to published values, 1 otherwise.
    """
    results = benchmark_sharp_li_cai()
    print_benchmark_report(results)

    # Exit code indicates whether Sharp-Li is closer
    if results.get("sharp_li_is_closer", False):
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
