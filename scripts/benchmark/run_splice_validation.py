#!/usr/bin/env python3
"""
Run Splice Scoring Validation
================================

Script that runs the splice scoring validation comparing the deprecated
splicing.maxent_score() against the proper maxentscan.score_donor()/
score_acceptor() implementations, and prints a detailed report.

Usage::

    python scripts/run_splice_validation.py

Task ID: F5.4
"""

from __future__ import annotations

import sys


def main() -> int:
    """Run the splice scoring validation and print the report."""
    try:
        from biocompiler.benchmarking.splice_scoring_validation import (
            validate_splice_scoring,
            print_splice_validation_report,
        )
    except ImportError:
        # Fallback for development / editable installs
        from src.biocompiler.benchmarking.splice_scoring_validation import (
            validate_splice_scoring,
            print_splice_validation_report,
        )

    result = validate_splice_scoring()
    print_splice_validation_report(result)

    # Exit with non-zero code if anti-correlation is not detected
    # (this would indicate the validation itself has a problem)
    if not result.is_anti_correlated:
        print(
            "WARNING: Anti-correlation was NOT detected. "
            "This is unexpected — the deprecated maxent_score() should be "
            "anti-correlated with the proper MaxEntScan scoring.",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
