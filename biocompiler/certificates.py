"""
BioCompiler Certificate Generation v7.0.0
============================================
Generates optimization certificates with GOLD/SILVER/BRONZE levels.
"""

from typing import List, Optional
from .type_system import CertLevel, PredicateResult


def compute_certificate(results: List[PredicateResult]) -> CertLevel:
    """Compute certificate level from predicate results.

    GOLD:   All predicates satisfied by optimization alone
    SILVER: All predicates satisfied, some required mutagenesis
    BRONZE: Some predicates could not be fully satisfied
    """
    has_mutagenesis = False
    has_unsatisfied = False

    for r in results:
        if not r.passed:
            has_unsatisfied = True
        if "mutagenesis" in r.details.lower():
            has_mutagenesis = True

    if has_unsatisfied:
        return CertLevel.BRONZE
    elif has_mutagenesis:
        return CertLevel.SILVER
    else:
        return CertLevel.GOLD


def format_certificate(results: List[PredicateResult], seq: str, species: str) -> str:
    """Format a human-readable certificate report."""
    cert = compute_certificate(results)
    lines = [
        "=" * 60,
        "  BioCompiler v7.0.0 — Optimization Certificate",
        "=" * 60,
        f"  Sequence length: {len(seq)} bp",
        f"  Species:         {species}",
        f"  Certificate:     {cert.value}",
        "-" * 60,
        "  Predicate Results:",
    ]
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        verdict_str = f" [{r.verdict.name}]" if r.verdict else ""
        lines.append(f"    [{status}{verdict_str}] {r.predicate}: {r.details}")
    lines.append("=" * 60)
    return "\n".join(lines)
