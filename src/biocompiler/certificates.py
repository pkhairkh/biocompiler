"""
BioCompiler Certificate Generation v7.0.0
============================================
Generates optimization certificates with GOLD/SILVER/BRONZE levels.

GOLD:   All predicates satisfied by optimization alone
SILVER: All predicates satisfied, some required mutagenesis or have
        unavoidable constraints (e.g., Valine GT dinucleotides)
BRONZE: Some predicates could not be fully satisfied
"""

from typing import List, Optional
from .type_system import CertLevel, PredicateResult


def compute_certificate(results: List[PredicateResult]) -> CertLevel:
    """Compute certificate level from predicate results.

    GOLD:   All predicates satisfied by optimization alone
    SILVER: All predicates satisfied, some required mutagenesis or have
            unavoidable constraints (e.g., Valine GT dinucleotides that
            can only be removed by AA substitution)
    BRONZE: Some predicates could not be fully satisfied
    """
    has_mutagenesis = False
    has_unavoidable = False
    has_unsatisfied = False

    for r in results:
        if not r.passed:
            has_unsatisfied = True
        if "mutagenesis" in r.details.lower():
            has_mutagenesis = True
        if "unavoidable" in r.details.lower():
            has_unavoidable = True

    if has_unsatisfied:
        return CertLevel.BRONZE
    elif has_mutagenesis or has_unavoidable:
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
        # Mark mutagenesis-resolved constraints
        mutagenesis_marker = ""
        if "mutagenesis" in r.details.lower():
            mutagenesis_marker = " [MUTAGENESIS]"
        elif "unavoidable" in r.details.lower() and r.passed:
            mutagenesis_marker = " [UNAVOIDABLE]"
        lines.append(f"    [{status}{verdict_str}{mutagenesis_marker}] {r.predicate}: {r.details}")
    lines.append("=" * 60)
    lines.append("")
    lines.append("  Certificate Levels:")
    lines.append("    GOLD   — All constraints satisfied by synonymous optimization")
    lines.append("    SILVER — All constraints satisfied (some required AA substitution)")
    lines.append("    BRONZE — Some constraints could not be satisfied")
    lines.append("=" * 60)
    return "\n".join(lines)
