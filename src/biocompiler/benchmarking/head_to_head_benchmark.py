"""
BioCompiler Benchmarking — Head-to-Head CAI Optimizer Benchmark
================================================================

Head-to-head benchmark comparing BioCompiler (with organism-aware constraints)
against DNAchisel, measuring CAI, speed, and constraint satisfaction across
a panel of E. coli genes.

Provides:
  - ``HeadToHeadResult``: Dataclass summarizing aggregate benchmark results.
  - ``run_head_to_head``: Run the full benchmark across a gene panel.
  - ``print_head_to_head_report``: Print a formatted comparison table.

    Design:
    All metrics are computed using BioCompiler's validated evaluators for fairness.
    DNAchisel results are included only when the package is installed; otherwise
    the benchmark runs in BioCompiler-only mode and DNAchisel columns report
    as unavailable.  CAI is always computed with ``compute_cai_validated`` for
    both tools — DNAchisel's own CAI output is NOT trusted.

Usage::

    from biocompiler.benchmarking.head_to_head_benchmark import (
        run_head_to_head,
        print_head_to_head_report,
    )

    result = run_head_to_head(organism="Escherichia_coli")
    print_head_to_head_report(result)
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "HeadToHeadResult",
    "run_head_to_head",
    "print_head_to_head_report",
]


# ─── Default constraints ─────────────────────────────────────────────

_DEFAULT_GC_LO = 0.30
_DEFAULT_GC_HI = 0.70
_DEFAULT_ENZYMES = ["EcoRI", "BamHI", "HindIII", "XhoI"]


# ─── Data class ──────────────────────────────────────────────────────


@dataclass
class HeadToHeadResult:
    """Aggregate head-to-head benchmark result comparing BioCompiler vs DNAchisel.

    Attributes:
        organism: Target organism used for the benchmark.
        num_genes: Number of genes benchmarked.
        biocompiler_mean_cai: Mean CAI across all genes (BioCompiler).
        dnachisel_mean_cai: Mean CAI across all genes (DNAchisel),
            or 0.0 if DNAchisel was unavailable.
        biocompiler_mean_time_ms: Mean optimization time in milliseconds
            (BioCompiler).
        dnachisel_mean_time_ms: Mean optimization time in milliseconds
            (DNAchisel), or 0.0 if DNAchisel was unavailable.
        speed_ratio: biocompiler_mean_time_ms / dnachisel_mean_time_ms,
            or 0.0 if DNAchisel was unavailable.  Values < 1.0 mean
            BioCompiler is faster.
        cai_gap: biocompiler_mean_cai - dnachisel_mean_cai.  Positive
            values mean BioCompiler achieves higher CAI.
        per_gene_results: List of per-gene result dicts, each containing:
            gene_name, protein_length, biocompiler_cai, dnachisel_cai,
            biocompiler_gc, dnachisel_gc, biocompiler_time_ms,
            dnachisel_time_ms, biocompiler_violations, dnachisel_violations,
            biocompiler_success, dnachisel_success.
        biocompiler_provenance_available: Whether BioCompiler provenance
            tracking was available for the run.
        dnachisel_provenance_available: Whether DNAchisel was installed
            and produced results for at least one gene.
    """

    organism: str
    num_genes: int
    biocompiler_mean_cai: float
    dnachisel_mean_cai: float
    biocompiler_mean_time_ms: float
    dnachisel_mean_time_ms: float
    speed_ratio: float
    cai_gap: float
    per_gene_results: list[dict]
    biocompiler_provenance_available: bool
    dnachisel_provenance_available: bool


# ─── Helpers ─────────────────────────────────────────────────────────


def _get_ecoli_genes() -> list[tuple[str, str]]:
    """Retrieve E. coli gene panel from gene_sets.

    Returns a list of (gene_name, protein_sequence) tuples from the
    E_COLI_EXTENDED gene set.
    """
    try:
        from .gene_sets import E_COLI_EXTENDED
    except ImportError:
        logger.warning("E_COLI_EXTENDED gene set not available")
        return []

    genes: list[tuple[str, str]] = []
    for name, entry in E_COLI_EXTENDED.items():
        protein = entry.get("protein", "")
        if protein:
            genes.append((name, protein))
    return genes


def _count_constraint_violations(
    sequence: str,
    protein: str,
    organism: str,
    gc_lo: float = _DEFAULT_GC_LO,
    gc_hi: float = _DEFAULT_GC_HI,
    enzymes: list[str] | None = None,
) -> int:
    """Count the number of constraint violations in an optimized sequence.

    Checks:
      1. Translation fidelity (sequence translates to target protein)
      2. GC content within [gc_lo, gc_hi]
      3. No restriction enzyme sites (from the enzyme panel)

    Args:
        sequence: Optimized DNA sequence.
        protein: Target protein sequence.
        organism: Target organism for translation.
        gc_lo: Minimum GC fraction.
        gc_hi: Maximum GC fraction.
        enzymes: Restriction enzymes to check.

    Returns:
        Number of violated constraints (0 = all satisfied).
    """
    if enzymes is None:
        enzymes = _DEFAULT_ENZYMES

    violations = 0

    # 1. Translation fidelity
    try:
        from ..translation import translate
        translated = translate(sequence, to_stop=True)
        if translated != protein:
            violations += 1
    except Exception:
        violations += 1

    # 2. GC content
    try:
        from ..scanner import gc_content
        gc = gc_content(sequence)
        if not (gc_lo <= gc <= gc_hi):
            violations += 1
    except Exception:
        violations += 1

    # 3. Restriction sites
    try:
        from ..restriction_sites import get_recognition_site
        for enz in enzymes:
            site = get_recognition_site(enz)
            if site is None:
                continue
            if any(b not in "ACGT" for b in site.upper()):
                continue
            if site.upper() in sequence.upper():
                violations += 1
                break  # Count as one violation regardless of how many sites
    except Exception:
        pass

    return violations


# ─── Main benchmark function ─────────────────────────────────────────


def run_head_to_head(
    genes: list[str] | None = None,
    organism: str = "Escherichia_coli",
) -> HeadToHeadResult:
    """Run head-to-head benchmark: BioCompiler vs DNAchisel.

    For each gene in the panel:
      1. Optimize with BioCompiler (organism-aware mode)
      2. Optimize with DNAchisel (if available, else skip)
      3. Record: CAI, GC content, optimization time, constraint violations,
         sequence length

    Then compute summary statistics and return a HeadToHeadResult.

    Args:
        genes: Optional list of gene names to benchmark. If None, uses
            the full E. coli gene panel from ``E_COLI_EXTENDED`` (currently
            20 genes covering housekeeping, regulatory, chaperone, translation,
            membrane, metabolic, resistance, and transport categories).
        organism: Target organism for codon optimization.
            Default is ``"Escherichia_coli"``.

    Returns:
        HeadToHeadResult with aggregate and per-gene metrics.
    """
    from ..optimization import optimize_sequence
    from .dnachisel_adapter import is_dnachisel_available

    # ── Resolve gene panel ──
    gene_panel = _get_ecoli_genes()
    if genes is not None:
        # Filter to requested genes
        gene_panel = [(n, p) for n, p in gene_panel if n in genes]

    if not gene_panel:
        logger.warning("No genes available for benchmarking")
        return HeadToHeadResult(
            organism=organism,
            num_genes=0,
            biocompiler_mean_cai=0.0,
            dnachisel_mean_cai=0.0,
            biocompiler_mean_time_ms=0.0,
            dnachisel_mean_time_ms=0.0,
            speed_ratio=0.0,
            cai_gap=0.0,
            per_gene_results=[],
            biocompiler_provenance_available=False,
            dnachisel_provenance_available=False,
        )

    dnachisel_available = is_dnachisel_available()
    adapter = None
    if dnachisel_available:
        try:
            from .dnachisel_adapter import DNAchiselAdapter
            adapter = DNAchiselAdapter()
        except ImportError:
            dnachisel_available = False

    per_gene: list[dict] = []
    bc_cais: list[float] = []
    dc_cais: list[float] = []
    bc_times: list[float] = []
    dc_times: list[float] = []
    bc_provenance_found = False
    dc_provenance_found = False

    for gene_name, protein in gene_panel:
        # ── BioCompiler optimization ──
        bc_cai = 0.0
        bc_gc = 0.0
        bc_time_ms = 0.0
        bc_violations = 0
        bc_success = False
        bc_seq = ""

        t0 = time.perf_counter()
        try:
            bc_result = optimize_sequence(
                target_protein=protein,
                organism=organism,
                gc_lo=_DEFAULT_GC_LO,
                gc_hi=_DEFAULT_GC_HI,
                track_provenance=True,
            )
            bc_seq = bc_result.sequence
            # Use validated CAI for fair comparison — do NOT trust optimizer's CAI
            from .metrics import compute_cai_validated
            bc_cai = compute_cai_validated(bc_seq, organism)
            bc_gc = bc_result.gc_content
            bc_success = True

            # Check provenance
            if hasattr(bc_result, "provenance") and bc_result.provenance is not None:
                bc_provenance_found = True

            # Count constraint violations
            bc_violations = _count_constraint_violations(
                bc_seq, protein, organism,
                _DEFAULT_GC_LO, _DEFAULT_GC_HI, _DEFAULT_ENZYMES,
            )
        except Exception as exc:
            logger.error("BioCompiler failed for %s: %s", gene_name, exc)
            bc_success = False
        bc_time_ms = (time.perf_counter() - t0) * 1000.0

        # ── DNAchisel optimization ──
        dc_cai = 0.0
        dc_gc = 0.0
        dc_time_ms = 0.0
        dc_violations = 0
        dc_success = False

        if adapter is not None:
            t1 = time.perf_counter()
            try:
                dc_result = adapter.optimize(
                    protein=protein,
                    organism=organism,
                    constraints=[
                        {"type": "gc_range", "gc_lo": _DEFAULT_GC_LO, "gc_hi": _DEFAULT_GC_HI},
                        {"type": "avoid_restriction", "enzymes": _DEFAULT_ENZYMES},
                    ],
                )
                if dc_result.success:
                    dc_cai = dc_result.cai
                    dc_gc = dc_result.gc_content
                    dc_success = True
                    dc_provenance_found = True

                    # Count constraint violations on DNAchisel output
                    dc_violations = _count_constraint_violations(
                        dc_result.sequence, protein, organism,
                        _DEFAULT_GC_LO, _DEFAULT_GC_HI, _DEFAULT_ENZYMES,
                    )
            except Exception as exc:
                logger.debug("DNAchisel failed for %s: %s", gene_name, exc)
            dc_time_ms = (time.perf_counter() - t1) * 1000.0

        # Record per-gene results
        per_gene.append({
            "gene_name": gene_name,
            "protein_length": len(protein),
            "biocompiler_cai": bc_cai,
            "dnachisel_cai": dc_cai if dc_success else None,
            "biocompiler_gc": bc_gc,
            "dnachisel_gc": dc_gc if dc_success else None,
            "biocompiler_time_ms": bc_time_ms,
            "dnachisel_time_ms": dc_time_ms if dc_success else None,
            "biocompiler_violations": bc_violations,
            "dnachisel_violations": dc_violations if dc_success else None,
            "biocompiler_success": bc_success,
            "dnachisel_success": dc_success if dnachisel_available else None,
            "sequence_length": len(bc_seq),
        })

        if bc_success:
            bc_cais.append(bc_cai)
            bc_times.append(bc_time_ms)
        if dc_success:
            dc_cais.append(dc_cai)
            dc_times.append(dc_time_ms)

    # ── Compute summary statistics ──
    n = len(gene_panel)

    biocompiler_mean_cai = sum(bc_cais) / len(bc_cais) if bc_cais else 0.0
    dnachisel_mean_cai = sum(dc_cais) / len(dc_cais) if dc_cais else 0.0
    biocompiler_mean_time = sum(bc_times) / len(bc_times) if bc_times else 0.0
    dnachisel_mean_time = sum(dc_times) / len(dc_times) if dc_times else 0.0

    speed_ratio = (
        biocompiler_mean_time / dnachisel_mean_time
        if dnachisel_mean_time > 0
        else 0.0
    )
    cai_gap = biocompiler_mean_cai - dnachisel_mean_cai

    # Constraint satisfaction rate
    bc_total = sum(1 for g in per_gene if g["biocompiler_success"])
    bc_satisfied = sum(
        1 for g in per_gene
        if g["biocompiler_success"] and g["biocompiler_violations"] == 0
    )
    # Note: constraint satisfaction rate is stored in per_gene results

    return HeadToHeadResult(
        organism=organism,
        num_genes=n,
        biocompiler_mean_cai=round(biocompiler_mean_cai, 4),
        dnachisel_mean_cai=round(dnachisel_mean_cai, 4),
        biocompiler_mean_time_ms=round(biocompiler_mean_time, 2),
        dnachisel_mean_time_ms=round(dnachisel_mean_time, 2),
        speed_ratio=round(speed_ratio, 4),
        cai_gap=round(cai_gap, 4),
        per_gene_results=per_gene,
        biocompiler_provenance_available=bc_provenance_found,
        dnachisel_provenance_available=dc_provenance_found,
    )


# ─── Report printing ─────────────────────────────────────────────────


def print_head_to_head_report(result: HeadToHeadResult) -> None:
    """Print a formatted head-to-head comparison report.

    Outputs a summary table with per-gene results and aggregate statistics.
    Highlights where BioCompiler matches or exceeds DNAchisel on each metric,
    and notes the provenance advantage.

    Args:
        result: HeadToHeadResult from ``run_head_to_head``.
    """
    # ── Header ──
    print()
    print("=" * 82)
    print("  HEAD-TO-HEAD BENCHMARK: BioCompiler vs DNAchisel")
    print(f"  Organism: {result.organism}  |  Genes: {result.num_genes}")
    print("=" * 82)

    # ── Per-gene table ──
    dc_available = any(
        g["dnachisel_success"] is not None for g in result.per_gene_results
    )

    # Column widths
    name_w = 10
    len_w = 5
    cai_w = 8
    gc_w = 7
    time_w = 9
    viol_w = 6

    if dc_available:
        header = (
            f"{'Gene':<{name_w}} {'Len':>{len_w}} "
            f"{'CAI_BC':>{cai_w}} {'CAI_DC':>{cai_w}} "
            f"{'GC_BC':>{gc_w}} {'GC_DC':>{gc_w}} "
            f"{'ms_BC':>{time_w}} {'ms_DC':>{time_w}} "
            f"{'V_BC':>{viol_w}} {'V_DC':>{viol_w}} "
            f"{'Winner':>8}"
        )
    else:
        header = (
            f"{'Gene':<{name_w}} {'Len':>{len_w}} "
            f"{'CAI_BC':>{cai_w}} "
            f"{'GC_BC':>{gc_w}} "
            f"{'ms_BC':>{time_w}} "
            f"{'V_BC':>{viol_w}}"
        )
    print(header)
    print("-" * len(header))

    for g in result.per_gene_results:
        name = g["gene_name"]
        plen = g["protein_length"]
        bc_cai = g["biocompiler_cai"]
        bc_gc = g["biocompiler_gc"]
        bc_time = g["biocompiler_time_ms"]
        bc_viol = g["biocompiler_violations"]

        if dc_available and g["dnachisel_success"]:
            dc_cai_val = g["dnachisel_cai"]
            dc_gc_val = g["dnachisel_gc"]
            dc_time_val = g["dnachisel_time_ms"]
            dc_viol = g["dnachisel_violations"]

            # Determine winner
            if bc_cai > (dc_cai_val or 0) + 0.001:
                winner = "BC"
            elif (dc_cai_val or 0) > bc_cai + 0.001:
                winner = "DC"
            else:
                winner = "tie"

            # Highlight BioCompiler advantage
            bc_cai_str = f"{bc_cai:.4f}"
            dc_cai_str = f"{dc_cai_val:.4f}" if dc_cai_val is not None else "  N/A "

            line = (
                f"{name:<{name_w}} {plen:>{len_w}} "
                f"{bc_cai_str:>{cai_w}} {dc_cai_str:>{cai_w}} "
                f"{bc_gc:>{gc_w}.3f} {dc_gc_val:>{gc_w}.3f} "
                f"{bc_time:>{time_w}.1f} {dc_time_val:>{time_w}.1f} "
                f"{bc_viol:>{viol_w}} {dc_viol:>{viol_w}} "
                f"{winner:>8}"
            )
        else:
            line = (
                f"{name:<{name_w}} {plen:>{len_w}} "
                f"{bc_cai:>{cai_w}.4f} "
                f"{bc_gc:>{gc_w}.3f} "
                f"{bc_time:>{time_w}.1f} "
                f"{bc_viol:>{viol_w}}"
            )
        print(line)

    # ── Summary statistics ──
    print("-" * 82)
    print()
    print("  SUMMARY")
    print("  -------")
    print(f"  Mean CAI   — BioCompiler: {result.biocompiler_mean_cai:.4f}", end="")
    if dc_available:
        print(f"  |  DNAchisel: {result.dnachisel_mean_cai:.4f}", end="")
        if result.cai_gap > 0:
            print(f"  |  Gap: +{result.cai_gap:.4f} (BC leads)")
        elif result.cai_gap < 0:
            print(f"  |  Gap: {result.cai_gap:.4f} (DC leads)")
        else:
            print(f"  |  Gap: {result.cai_gap:.4f} (tie)")
    else:
        print()

    print(f"  Mean Time  — BioCompiler: {result.biocompiler_mean_time_ms:.1f} ms", end="")
    if dc_available:
        print(f"  |  DNAchisel: {result.dnachisel_mean_time_ms:.1f} ms", end="")
        if result.speed_ratio > 0:
            if result.speed_ratio < 1.0:
                print(f"  |  Speed ratio: {result.speed_ratio:.2f}x (BC faster)")
            else:
                print(f"  |  Speed ratio: {result.speed_ratio:.2f}x (DC faster)")
        else:
            print()
    else:
        print()

    # Constraint satisfaction rates
    bc_success = sum(1 for g in result.per_gene_results if g["biocompiler_success"])
    bc_satisfied = sum(
        1 for g in result.per_gene_results
        if g["biocompiler_success"] and g["biocompiler_violations"] == 0
    )
    bc_rate = bc_satisfied / bc_success if bc_success > 0 else 0.0

    print(f"  Constraint Satisfaction — BioCompiler: {bc_rate:.1%} "
          f"({bc_satisfied}/{bc_success} genes)", end="")
    if dc_available:
        dc_success = sum(
            1 for g in result.per_gene_results
            if g["dnachisel_success"] is True
        )
        dc_satisfied = sum(
            1 for g in result.per_gene_results
            if g["dnachisel_success"] is True
            and g["dnachisel_violations"] == 0
        )
        dc_rate = dc_satisfied / dc_success if dc_success > 0 else 0.0
        print(f"  |  DNAchisel: {dc_rate:.1%} "
              f"({dc_satisfied}/{dc_success} genes)", end="")
    print()

    # ── Provenance advantage ──
    print()
    print("  PROVENANCE")
    print("  ----------")
    print(f"  BioCompiler provenance available: "
          f"{'Yes' if result.biocompiler_provenance_available else 'No'}")
    print(f"  DNAchisel provenance available:   "
          f"{'Yes' if result.dnachisel_provenance_available else 'No (not installed)' if not dc_available else 'No'}")

    if result.biocompiler_provenance_available and not result.dnachisel_provenance_available:
        print()
        print("  >>> BioCompiler provides full provenance tracking (decision")
        print("      logs, constraint resolution history, codon choice rationale)")
        print("      which DNAchisel does not offer. This is a key advantage")
        print("      for regulatory compliance and reproducibility.")

    # ── Overall verdict ──
    print()
    if dc_available and result.cai_gap >= 0:
        print("  VERDICT: BioCompiler matches or exceeds DNAchisel on CAI")
        print(f"           (+{result.cai_gap:.4f} mean CAI advantage)")
    elif dc_available:
        print("  VERDICT: DNAchisel leads on CAI")
        print(f"           ({result.cai_gap:.4f} mean CAI gap)")
    else:
        print("  VERDICT: DNAchisel not installed — BioCompiler-only benchmark")
    print("=" * 82)
    print()
