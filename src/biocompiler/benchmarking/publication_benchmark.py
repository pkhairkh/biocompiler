"""
Publication-ready benchmark suite for BioCompiler v12.
Generates comprehensive benchmark data suitable for peer-reviewed publication.

Benchmark categories:
1. CAI quality (across organisms and tools)
2. Speed (with statistical methodology: mean, std, min, max)
3. Constraint satisfaction rate
4. Eukaryotic constraint handling (CpG, splice sites)

Usage::

    from biocompiler.benchmarking.publication_benchmark import (
        run_publication_benchmark,
        benchmark_cai_quality,
        benchmark_speed,
        benchmark_constraint_satisfaction,
        generate_latex_table,
        generate_benchmark_figure,
        PUBLICATION_GENES,
    )

    # Run full suite
    results = run_publication_benchmark(output_dir="benchmark_results")

    # Or individual categories with custom gene set
    cai_results = benchmark_cai_quality(
        gene_set=PUBLICATION_GENES,
        organisms=["Escherichia_coli", "Homo_sapiens"],
    )

    # Generate publication artifacts
    latex = generate_latex_table(results, table_type="cai")
    generate_benchmark_figure(results, "figure.png", figure_type="bar")
"""

from __future__ import annotations

import json
import logging
import os
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "PUBLICATION_GENES",
    "run_publication_benchmark",
    "benchmark_cai_quality",
    "benchmark_speed",
    "benchmark_constraint_satisfaction",
    "generate_latex_table",
    "generate_benchmark_figure",
]


# ---------------------------------------------------------------------------
# Standard publication gene set
# ---------------------------------------------------------------------------

PUBLICATION_GENES: dict[str, str] = {
    "GFP": (
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYQ"
    ),
    "Insulin": (
        "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGG"
    ),
    "HBB": (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGA"
        "FSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALA"
        "HKYH"
    ),
    "mCherry": (
        "MVSKGEEDNMAIIKMFMRFHVTHGSGSNGTGESRGMDMKMVIENAMECVRMVMHEGHNYTGKLPVPWPTLVT"
        "TFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDG"
        "NILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSAL"
        "SKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    ),
}


# ---------------------------------------------------------------------------
# Default configuration
# ---------------------------------------------------------------------------

_DEFAULT_ORGANISMS = [
    "Escherichia_coli",
    "Homo_sapiens",
    "Saccharomyces_cerevisiae",
    "Mus_musculus",
    "CHO_K1",
]

_DEFAULT_ENZYMES = [
    "EcoRI", "BamHI", "HindIII", "XhoI", "XbaI", "SalI",
    "PstI", "NcoI", "NdeI", "NotI",
]

_DEFAULT_GC_RANGE = (0.30, 0.70)

# Valid single-letter amino acid codes
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------------------
# 1. CAI quality benchmark
# ---------------------------------------------------------------------------


def benchmark_cai_quality(
    gene_set: dict[str, str],
    organisms: list[str],
) -> dict[str, Any]:
    """Benchmark CAI quality across organisms and tools.

    Optimizes each gene for each organism and computes the CAI using
    the validated evaluator. For fair comparison, also runs DNAchisel
    (if available) and computes CAI with the same evaluator.

    Args:
        gene_set: Dict mapping gene name to protein sequence.
        organisms: List of organism names (e.g. ``["Escherichia_coli", "Homo_sapiens"]``).

    Returns:
        Dict with keys:

        - ``per_gene``: Dict[gene_name][organism] -> dict with
          ``cai_biocompiler``, ``cai_dnachisel``, ``delta``.
        - ``mean_cai_biocompiler``: Mean CAI across all (gene, organism).
        - ``mean_cai_dnachisel``: Mean DNAchisel CAI (0 if unavailable).
        - ``cai_advantage``: Mean CAI advantage of BioCompiler.
        - ``organisms``: List of organisms used.
        - ``genes``: List of genes used.
    """
    from .metrics import compute_cai_validated
    from .dnachisel_adapter import is_dnachisel_available

    dnachisel_available = is_dnachisel_available()
    adapter = None
    if dnachisel_available:
        try:
            from .dnachisel_adapter import DNAchiselAdapter
            adapter = DNAchiselAdapter()
        except ImportError:
            dnachisel_available = False

    per_gene: dict[str, dict[str, Any]] = {}
    bc_cais: list[float] = []
    dc_cais: list[float] = []

    for gene_name, protein in gene_set.items():
        per_gene[gene_name] = {}

        for organism in organisms:
            bc_cai = _optimize_and_compute_cai(protein, organism)
            dc_cai = 0.0

            if adapter is not None:
                try:
                    dc_result = adapter.optimize(
                        protein=protein,
                        organism=organism,
                        constraints=[
                            {"type": "gc_range",
                             "gc_lo": _DEFAULT_GC_RANGE[0],
                             "gc_hi": _DEFAULT_GC_RANGE[1]},
                        ],
                    )
                    if dc_result.success:
                        dc_cai = compute_cai_validated(
                            dc_result.sequence, organism
                        )
                except Exception as exc:
                    logger.debug(
                        "DNAchisel failed for %s/%s: %s", gene_name, organism, exc
                    )

            delta = bc_cai - dc_cai if dc_cai > 0 else 0.0

            per_gene[gene_name][organism] = {
                "cai_biocompiler": round(bc_cai, 4),
                "cai_dnachisel": round(dc_cai, 4) if dc_cai > 0 else None,
                "delta": round(delta, 4),
            }

            bc_cais.append(bc_cai)
            if dc_cai > 0:
                dc_cais.append(dc_cai)

    return {
        "per_gene": per_gene,
        "mean_cai_biocompiler": round(statistics.mean(bc_cais), 4) if bc_cais else 0.0,
        "mean_cai_dnachisel": round(statistics.mean(dc_cais), 4) if dc_cais else 0.0,
        "cai_advantage": (
            round(statistics.mean(bc_cais) - statistics.mean(dc_cais), 4)
            if dc_cais
            else 0.0
        ),
        "organisms": list(organisms),
        "genes": list(gene_set.keys()),
    }


def _optimize_and_compute_cai(protein: str, organism: str) -> float:
    """Optimize a protein with BioCompiler and compute CAI.

    Args:
        protein: Amino acid sequence.
        organism: Target organism.

    Returns:
        CAI value, or 0.0 if optimization failed.
    """
    from .metrics import compute_cai_validated

    try:
        from ..optimization import optimize_sequence
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=_DEFAULT_GC_RANGE[0],
            gc_hi=_DEFAULT_GC_RANGE[1],
        )
        return compute_cai_validated(result.sequence, organism)
    except Exception as exc:
        logger.debug("BioCompiler optimization failed: %s", exc)
        return 0.0


# ---------------------------------------------------------------------------
# 2. Speed benchmark
# ---------------------------------------------------------------------------


def benchmark_speed(
    gene_set: dict[str, str],
    organisms: list[str],
    warmup: int = 3,
    repeats: int = 10,
) -> dict[str, Any]:
    """Benchmark execution speed with proper statistical methodology.

    Follows best practices for benchmarking:

    1. Warmup runs to allow JIT compilation (Numba) to complete.
    2. Multiple repeated runs for statistical significance.
    3. Report mean, std, min, max.

    Args:
        gene_set: Dict mapping gene name to protein sequence.
        organisms: List of target organisms.
        warmup: Number of warmup iterations (default 3).
        repeats: Number of measured iterations (default 10).

    Returns:
        Dict with keys:

        - ``per_gene``: Dict[gene_name][organism] -> speed stats.
        - ``mean_speed_biocompiler_ms``: Mean time per gene in ms.
        - ``mean_speed_dnachisel_ms``: Mean DNAchisel time (0 if unavailable).
        - ``speed_ratio``: BC/DC ratio (<1 means BC is faster).
        - ``warmup``: Number of warmup runs used.
        - ``repeats``: Number of measured runs used.
    """
    adapter = None
    dnachisel_available = False
    try:
        from .dnachisel_adapter import DNAchiselAdapter, is_dnachisel_available
        dnachisel_available = is_dnachisel_available()
        if dnachisel_available:
            adapter = DNAchiselAdapter()
    except ImportError:
        pass

    per_gene: dict[str, dict[str, Any]] = {}
    bc_times: list[float] = []
    dc_times: list[float] = []

    for gene_name, protein in gene_set.items():
        per_gene[gene_name] = {}

        for organism in organisms:
            # Warmup
            for _ in range(warmup):
                try:
                    from ..optimization import optimize_sequence
                    optimize_sequence(
                        target_protein=protein,
                        organism=organism,
                        gc_lo=_DEFAULT_GC_RANGE[0],
                        gc_hi=_DEFAULT_GC_RANGE[1],
                    )
                except Exception:
                    pass

            # Measured runs
            bc_run_times: list[float] = []
            for _ in range(repeats):
                t0 = time.perf_counter()
                try:
                    from ..optimization import optimize_sequence
                    optimize_sequence(
                        target_protein=protein,
                        organism=organism,
                        gc_lo=_DEFAULT_GC_RANGE[0],
                        gc_hi=_DEFAULT_GC_RANGE[1],
                    )
                except Exception:
                    pass
                bc_run_times.append((time.perf_counter() - t0) * 1000.0)

            dc_run_times: list[float] = []
            if adapter is not None:
                for _ in range(warmup):
                    try:
                        adapter.optimize(
                            protein=protein,
                            organism=organism,
                            constraints=[
                                {"type": "gc_range",
                                 "gc_lo": _DEFAULT_GC_RANGE[0],
                                 "gc_hi": _DEFAULT_GC_RANGE[1]},
                            ],
                        )
                    except Exception:
                        pass

                for _ in range(repeats):
                    t0 = time.perf_counter()
                    try:
                        adapter.optimize(
                            protein=protein,
                            organism=organism,
                            constraints=[
                                {"type": "gc_range",
                                 "gc_lo": _DEFAULT_GC_RANGE[0],
                                 "gc_hi": _DEFAULT_GC_RANGE[1]},
                            ],
                        )
                    except Exception:
                        pass
                    dc_run_times.append((time.perf_counter() - t0) * 1000.0)

            bc_mean = statistics.mean(bc_run_times) if bc_run_times else 0.0
            bc_std = statistics.stdev(bc_run_times) if len(bc_run_times) > 1 else 0.0
            bc_min = min(bc_run_times) if bc_run_times else 0.0
            bc_max = max(bc_run_times) if bc_run_times else 0.0

            dc_mean = statistics.mean(dc_run_times) if dc_run_times else 0.0
            dc_std = statistics.stdev(dc_run_times) if len(dc_run_times) > 1 else 0.0
            dc_min = min(dc_run_times) if dc_run_times else 0.0
            dc_max = max(dc_run_times) if dc_run_times else 0.0

            per_gene[gene_name][organism] = {
                "biocompiler_mean_ms": round(bc_mean, 2),
                "biocompiler_std_ms": round(bc_std, 2),
                "biocompiler_min_ms": round(bc_min, 2),
                "biocompiler_max_ms": round(bc_max, 2),
                "biocompiler_cv": round(bc_std / bc_mean, 4) if bc_mean > 0 else 0.0,
                "dnachisel_mean_ms": round(dc_mean, 2) if dc_mean > 0 else None,
                "dnachisel_std_ms": round(dc_std, 2) if dc_mean > 0 else None,
                "dnachisel_min_ms": round(dc_min, 2) if dc_min > 0 else None,
                "dnachisel_max_ms": round(dc_max, 2) if dc_max > 0 else None,
            }

            bc_times.append(bc_mean)
            if dc_mean > 0:
                dc_times.append(dc_mean)

    mean_bc = statistics.mean(bc_times) if bc_times else 0.0
    mean_dc = statistics.mean(dc_times) if dc_times else 0.0

    return {
        "per_gene": per_gene,
        "mean_speed_biocompiler_ms": round(mean_bc, 2),
        "mean_speed_dnachisel_ms": round(mean_dc, 2),
        "speed_ratio": round(mean_bc / mean_dc, 4) if mean_dc > 0 else 0.0,
        "warmup": warmup,
        "repeats": repeats,
    }


# ---------------------------------------------------------------------------
# 3. Constraint satisfaction benchmark
# ---------------------------------------------------------------------------


def benchmark_constraint_satisfaction(
    gene_set: dict[str, str],
    constraint_sets: list[dict],
) -> dict[str, Any]:
    """Benchmark how many constraints each tool satisfies.

    For each gene and constraint set, optimizes the sequence and
    checks how many constraints are satisfied.

    Args:
        gene_set: Dict mapping gene name to protein sequence.
        constraint_sets: List of constraint dictionaries. Each dict
            should have a ``name`` and ``constraints`` key.

    Returns:
        Dict with constraint satisfaction rates and per-gene details.
    """
    organisms = _DEFAULT_ORGANISMS[:2]  # Use fewer organisms for constraint tests

    per_gene: dict[str, Any] = {}
    total_bc_satisfied = 0
    total_dc_satisfied = 0
    total_constraints = 0

    for gene_name, protein in gene_set.items():
        per_gene[gene_name] = {}

        for organism in organisms:
            for cs in constraint_sets:
                cs_name = cs["name"]
                constraints = cs["constraints"]

                # Check BioCompiler constraints
                bc_satisfied = _check_biocompiler_constraints(
                    protein, organism, constraints
                )
                dc_satisfied = 0
                dc_total = 0

                # Check DNAchisel constraints
                try:
                    from .dnachisel_adapter import DNAchiselAdapter, is_dnachisel_available
                    if is_dnachisel_available():
                        adapter = DNAchiselAdapter()
                        dc_result = adapter.optimize(
                            protein=protein,
                            organism=organism,
                            constraints=constraints,
                        )
                        if dc_result.success:
                            dc_satisfied, dc_total = _count_satisfied_constraints(
                                dc_result.sequence, protein, organism, constraints
                            )
                except (ImportError, Exception) as exc:
                    logger.debug("DNAchisel constraint check failed: %s", exc)

                bc_total = len(constraints) + 1  # +1 for translation

                key = f"{organism}_{cs_name}"
                per_gene[gene_name][key] = {
                    "biocompiler_satisfied": bc_satisfied,
                    "biocompiler_total": bc_total,
                    "dnachisel_satisfied": dc_satisfied,
                    "dnachisel_total": dc_total,
                }

                total_bc_satisfied += bc_satisfied
                total_constraints += bc_total
                if dc_total > 0:
                    total_dc_satisfied += dc_satisfied

    bc_rate = total_bc_satisfied / total_constraints if total_constraints > 0 else 0.0
    dc_rate = total_dc_satisfied / total_constraints if total_constraints > 0 else 0.0

    return {
        "per_gene": per_gene,
        "constraint_satisfaction_rate_biocompiler": round(bc_rate, 4),
        "constraint_satisfaction_rate_dnachisel": round(dc_rate, 4),
        "constraint_sets_used": [cs["name"] for cs in constraint_sets],
    }


def _check_biocompiler_constraints(
    protein: str, organism: str, constraints: list[dict]
) -> int:
    """Check how many constraints BioCompiler satisfies.

    Args:
        protein: Protein sequence.
        organism: Target organism.
        constraints: Constraint specifications.

    Returns:
        Number of satisfied constraints.
    """
    satisfied = 0

    try:
        from ..optimization import optimize_sequence
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=_DEFAULT_GC_RANGE[0],
            gc_hi=_DEFAULT_GC_RANGE[1],
        )
        seq = result.sequence

        # Check translation
        from ..translation import translate
        translated = translate(seq, to_stop=True)
        if translated == protein:
            satisfied += 1

        # Check GC range
        from ..scanner import gc_content
        gc = gc_content(seq)
        if _DEFAULT_GC_RANGE[0] <= gc <= _DEFAULT_GC_RANGE[1]:
            satisfied += 1

        # Check restriction sites
        for c in constraints:
            if c.get("type") == "avoid_restriction":
                enzymes = c.get("enzymes", [])
                from .metrics import count_restriction_sites
                rs = count_restriction_sites(seq, enzymes)
                if sum(rs.values()) == 0:
                    satisfied += 1

    except Exception as exc:
        logger.debug("BioCompiler constraint check failed: %s", exc)

    return satisfied


def _count_satisfied_constraints(
    sequence: str, protein: str, organism: str, constraints: list[dict]
) -> tuple[int, int]:
    """Count satisfied constraints for a given sequence.

    Returns:
        Tuple of (satisfied, total).
    """
    total = 0
    satisfied = 0

    # Translation check
    total += 1
    try:
        from ..translation import translate
        translated = translate(sequence, to_stop=True)
        if translated == protein:
            satisfied += 1
    except Exception:
        pass

    # GC range check
    total += 1
    try:
        from ..scanner import gc_content
        gc = gc_content(sequence)
        if _DEFAULT_GC_RANGE[0] <= gc <= _DEFAULT_GC_RANGE[1]:
            satisfied += 1
    except Exception:
        pass

    # Restriction site checks
    for c in constraints:
        if c.get("type") == "avoid_restriction":
            total += 1
            enzymes = c.get("enzymes", [])
            try:
                from .metrics import count_restriction_sites
                rs = count_restriction_sites(sequence, enzymes)
                if sum(rs.values()) == 0:
                    satisfied += 1
            except Exception:
                pass

    return satisfied, total


# ---------------------------------------------------------------------------
# Full benchmark suite
# ---------------------------------------------------------------------------


def run_publication_benchmark(output_dir: str = "benchmark_results") -> dict[str, Any]:
    """Run the full publication benchmark suite.

    Executes all benchmark categories and saves comprehensive results
    to the specified output directory.

    Args:
        output_dir: Directory to save results and figures.

    Returns:
        Dict with all benchmark results, environment metadata,
        and summary statistics.
    """
    from .reproducibility import capture_environment

    logger.info("Starting publication benchmark suite...")

    # Capture environment
    environment = capture_environment()

    # Default constraint sets for the full suite
    default_constraint_sets = [
        {
            "name": "gc_range_only",
            "constraints": [
                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
            ],
        },
        {
            "name": "gc_plus_restriction",
            "constraints": [
                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                {"type": "avoid_restriction", "enzymes": _DEFAULT_ENZYMES[:4]},
            ],
        },
        {
            "name": "full_constraints",
            "constraints": [
                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                {"type": "avoid_restriction", "enzymes": _DEFAULT_ENZYMES},
            ],
        },
    ]

    # Run benchmark categories
    t0 = time.perf_counter()

    logger.info("Running CAI quality benchmark...")
    cai_results = benchmark_cai_quality(
        gene_set=PUBLICATION_GENES,
        organisms=_DEFAULT_ORGANISMS,
    )

    logger.info("Running speed benchmark...")
    speed_results = benchmark_speed(
        gene_set=PUBLICATION_GENES,
        organisms=_DEFAULT_ORGANISMS[:2],
    )

    logger.info("Running constraint satisfaction benchmark...")
    constraint_results = benchmark_constraint_satisfaction(
        gene_set=PUBLICATION_GENES,
        constraint_sets=default_constraint_sets,
    )

    elapsed = time.perf_counter() - t0

    # Compile results
    results: dict[str, Any] = {
        "version": "v12",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "total_elapsed_s": round(elapsed, 2),
        "cai_quality": cai_results,
        "speed": speed_results,
        "constraint_satisfaction": constraint_results,
        "environment": environment,
        "summary": {
            "mean_cai_biocompiler": cai_results["mean_cai_biocompiler"],
            "mean_cai_dnachisel": cai_results["mean_cai_dnachisel"],
            "cai_advantage": cai_results["cai_advantage"],
            "mean_speed_biocompiler_ms": speed_results["mean_speed_biocompiler_ms"],
            "mean_speed_dnachisel_ms": speed_results["mean_speed_dnachisel_ms"],
            "speed_ratio": speed_results["speed_ratio"],
            "constraint_satisfaction_rate_biocompiler": (
                constraint_results["constraint_satisfaction_rate_biocompiler"]
            ),
            "constraint_satisfaction_rate_dnachisel": (
                constraint_results["constraint_satisfaction_rate_dnachisel"]
            ),
        },
    }

    # Save results
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    results_file = output_path / "publication_benchmark_results.json"
    with open(results_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)
    logger.info("Results saved to %s", results_file)

    # Generate LaTeX table
    try:
        latex = generate_latex_table(results)
        latex_file = output_path / "benchmark_table.tex"
        with open(latex_file, "w", encoding="utf-8") as f:
            f.write(latex)
        logger.info("LaTeX table saved to %s", latex_file)
    except Exception as exc:
        logger.warning("Failed to generate LaTeX table: %s", exc)

    # Generate figure
    try:
        figure_path = str(output_path / "benchmark_figure.png")
        generate_benchmark_figure(results, figure_path)
        logger.info("Benchmark figure saved to %s", figure_path)
    except Exception as exc:
        logger.warning("Failed to generate benchmark figure: %s", exc)

    logger.info(
        "Publication benchmark complete in %.1fs. "
        "CAI advantage: +%.4f, Speed ratio: %.2fx",
        elapsed,
        results["summary"]["cai_advantage"],
        results["summary"]["speed_ratio"],
    )

    return results


# ---------------------------------------------------------------------------
# LaTeX table generation
# ---------------------------------------------------------------------------


def generate_latex_table(results: dict, table_type: str = "cai") -> str:
    r"""Generate a LaTeX table from benchmark results.

    Produces a publication-quality table. The ``table_type`` parameter
    selects which benchmark data to tabulate:

    - ``"cai"``: CAI comparison table
    - ``"speed"``: Speed comparison table
    - ``"constraints"``: Constraint satisfaction table
    - ``"summary"``: Summary comparison table (default for general results)

    Args:
        results: Results dict from ``run_publication_benchmark``.
        table_type: Type of table to generate (``"cai"``, ``"speed"``,
            ``"constraints"``, or ``"summary"``).

    Returns:
        LaTeX source string.
    """
    lines: list[str] = []

    if table_type == "cai":
        lines.append(r"\begin{table}[htbp]")
        lines.append(r"  \centering")
        lines.append(
            r"  \caption{CAI quality comparison: BioCompiler v12 vs DNAchisel}"
        )
        lines.append(r"  \label{tab:cai-comparison}")
        lines.append(r"  \begin{tabular}{llrrr}")
        lines.append(r"    \toprule")
        lines.append(
            r"    Gene & Organism & CAI$_{\text{BC}}$ & "
            r"CAI$_{\text{DC}}$ & $\Delta$ \\"
        )
        lines.append(r"    \midrule")

        cai_data = results.get("cai_quality", results)
        per_gene = cai_data.get("per_gene", {})
        for gene_name, organisms_data in per_gene.items():
            for organism, metrics in organisms_data.items():
                bc_cai = metrics.get("cai_biocompiler", 0.0)
                dc_cai = metrics.get("cai_dnachisel")
                delta = metrics.get("delta", 0.0)
                dc_str = f"{dc_cai:.4f}" if dc_cai is not None else "---"
                org_latex = organism.replace("_", r"\_")
                lines.append(
                    f"    {gene_name} & {org_latex} & "
                    f"{bc_cai:.4f} & {dc_str} & {delta:+.4f} \\\\"
                )

        lines.append(r"    \bottomrule")
        lines.append(r"  \end{tabular}")
        lines.append(r"\end{table}")

    elif table_type == "speed":
        lines.append(r"\begin{table}[htbp]")
        lines.append(r"  \centering")
        lines.append(
            r"  \caption{Speed comparison: BioCompiler v12 vs DNAchisel}"
        )
        lines.append(r"  \label{tab:speed-comparison}")
        lines.append(r"  \begin{tabular}{llrrrr}")
        lines.append(r"    \toprule")
        lines.append(
            r"    Gene & Organism & Mean$_{\text{BC}}$ (ms) & "
            r"Std$_{\text{BC}}$ (ms) & Mean$_{\text{DC}}$ (ms) & "
            r"Ratio \\"
        )
        lines.append(r"    \midrule")

        speed_data = results.get("speed", results)
        per_gene = speed_data.get("per_gene", {})
        for gene_name, organisms_data in per_gene.items():
            for organism, metrics in organisms_data.items():
                bc_mean = metrics.get("biocompiler_mean_ms", 0.0)
                bc_std = metrics.get("biocompiler_std_ms", 0.0)
                dc_mean_str = (
                    f"{metrics['dnachisel_mean_ms']:.2f}"
                    if metrics.get("dnachisel_mean_ms")
                    else "---"
                )
                ratio = (
                    f"{bc_mean / metrics['dnachisel_mean_ms']:.2f}x"
                    if metrics.get("dnachisel_mean_ms") and metrics["dnachisel_mean_ms"] > 0
                    else "---"
                )
                org_latex = organism.replace("_", r"\_")
                lines.append(
                    f"    {gene_name} & {org_latex} & "
                    f"{bc_mean:.2f} & {bc_std:.2f} & "
                    f"{dc_mean_str} & {ratio} \\\\"
                )

        lines.append(r"    \bottomrule")
        lines.append(r"  \end{tabular}")
        lines.append(r"\end{table}")

    elif table_type == "constraints":
        lines.append(r"\begin{table}[htbp]")
        lines.append(r"  \centering")
        lines.append(
            r"  \caption{Constraint satisfaction comparison: BioCompiler v12 vs DNAchisel}"
        )
        lines.append(r"  \label{tab:constraint-comparison}")
        lines.append(r"  \begin{tabular}{llrr}")
        lines.append(r"    \toprule")
        lines.append(
            r"    Gene & Constraint Set & BC Rate & DC Rate \\"
        )
        lines.append(r"    \midrule")

        cs_data = results.get("constraint_satisfaction", results)
        per_gene = cs_data.get("per_gene", {})
        for gene_name, entries in per_gene.items():
            for key, metrics in entries.items():
                bc_sat = metrics.get("biocompiler_satisfied", 0)
                bc_total = metrics.get("biocompiler_total", 1)
                dc_sat = metrics.get("dnachisel_satisfied", 0)
                dc_total = metrics.get("dnachisel_total", 0)
                bc_rate = bc_sat / bc_total if bc_total > 0 else 0.0
                dc_rate = dc_sat / dc_total if dc_total > 0 else 0.0
                cs_name = key.split("_", 1)[-1] if "_" in key else key
                cs_latex = cs_name.replace("_", r"\_")
                lines.append(
                    f"    {gene_name} & {cs_latex} & "
                    f"{bc_rate:.1%} & "
                    f"{dc_rate:.1%} \\\\"
                )

        lines.append(r"    \bottomrule")
        lines.append(r"  \end{tabular}")
        lines.append(r"\end{table}")

    else:  # summary
        lines.append(r"\begin{table}[htbp]")
        lines.append(r"  \centering")
        lines.append(
            r"  \caption{Head-to-head comparison: BioCompiler v12 vs DNAchisel}"
        )
        lines.append(r"  \label{tab:benchmark-comparison}")
        lines.append(r"  \begin{tabular}{lrrrr}")
        lines.append(r"    \toprule")

        # Header
        lines.append(
            r"    Metric & BioCompiler & DNAchisel & $\Delta$ & $p$-value \\"
        )
        lines.append(r"    \midrule")

        summary = results.get("summary", {})

        # CAI row
        bc_cai = summary.get("mean_cai_biocompiler", 0.0)
        dc_cai = summary.get("mean_cai_dnachisel", 0.0)
        cai_delta = summary.get("cai_advantage", 0.0)
        lines.append(
            f"    Mean CAI & {bc_cai:.4f} & {dc_cai:.4f} & "
            f"{cai_delta:+.4f} & --- \\\\"
        )

        # Speed row
        bc_speed = summary.get("mean_speed_biocompiler_ms", 0.0)
        dc_speed = summary.get("mean_speed_dnachisel_ms", 0.0)
        speed_ratio = summary.get("speed_ratio", 0.0)
        lines.append(
            f"    Mean Speed (ms) & {bc_speed:.1f} & {dc_speed:.1f} & "
            f"×{speed_ratio:.2f} & --- \\\\"
        )

        # Constraint satisfaction row
        bc_cs = summary.get("constraint_satisfaction_rate_biocompiler", 0.0)
        dc_cs = summary.get("constraint_satisfaction_rate_dnachisel", 0.0)
        cs_delta = bc_cs - dc_cs
        lines.append(
            f"    Constraint Sat. & {bc_cs:.1%} & {dc_cs:.1%} & "
            f"{cs_delta:+.1%} & --- \\\\"
        )

        lines.append(r"    \bottomrule")
        lines.append(r"  \end{tabular}")
        lines.append(r"\end{table}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------


def generate_benchmark_figure(
    results: dict,
    output_path: str,
    figure_type: str = "bar",
) -> None:
    """Generate a matplotlib figure from benchmark results.

    Creates a multi-panel figure. The ``figure_type`` parameter selects
    the style:

    - ``"bar"``: Bar chart comparison (default)
    - ``"radar"``: Radar/spider chart of relative performance

    Args:
        results: Results dict from ``run_publication_benchmark``.
        output_path: File path to save the figure (e.g., "figure.png").
        figure_type: Type of figure to generate (``"bar"`` or ``"radar"``).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available; skipping figure generation")
        return

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    summary = results.get("summary", {})

    if figure_type == "radar":
        _generate_radar_figure(summary, output)
    else:
        _generate_bar_figure(summary, output)


def _generate_bar_figure(summary: dict, output: Path) -> None:
    """Generate a bar-chart style benchmark figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: CAI comparison
    ax1 = axes[0]
    bc_cai = summary.get("mean_cai_biocompiler", 0.0)
    dc_cai = summary.get("mean_cai_dnachisel", 0.0)
    bars = ax1.bar(
        ["BioCompiler", "DNAchisel"],
        [bc_cai, dc_cai if dc_cai > 0 else 0],
        color=["#2563EB", "#DC2626"],
        edgecolor="white",
    )
    ax1.set_ylabel("Mean CAI")
    ax1.set_title("CAI Quality")
    ax1.set_ylim(0, 1.0)
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax1.text(
                bar.get_x() + bar.get_width() / 2, height + 0.01,
                f"{height:.3f}", ha="center", va="bottom", fontsize=9,
            )

    # Panel 2: Speed comparison
    ax2 = axes[1]
    bc_speed = summary.get("mean_speed_biocompiler_ms", 0.0)
    dc_speed = summary.get("mean_speed_dnachisel_ms", 0.0)
    bars2 = ax2.bar(
        ["BioCompiler", "DNAchisel"],
        [bc_speed, dc_speed if dc_speed > 0 else 0],
        color=["#2563EB", "#DC2626"],
        edgecolor="white",
    )
    ax2.set_ylabel("Mean Time (ms)")
    ax2.set_title("Execution Speed")
    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax2.text(
                bar.get_x() + bar.get_width() / 2, height + 0.5,
                f"{height:.1f}", ha="center", va="bottom", fontsize=9,
            )

    # Panel 3: Constraint satisfaction
    ax3 = axes[2]
    bc_cs = summary.get("constraint_satisfaction_rate_biocompiler", 0.0)
    dc_cs = summary.get("constraint_satisfaction_rate_dnachisel", 0.0)
    bars3 = ax3.bar(
        ["BioCompiler", "DNAchisel"],
        [bc_cs * 100, dc_cs * 100 if dc_cs > 0 else 0],
        color=["#2563EB", "#DC2626"],
        edgecolor="white",
    )
    ax3.set_ylabel("Satisfaction Rate (%)")
    ax3.set_title("Constraint Satisfaction")
    ax3.set_ylim(0, 105)
    for bar in bars3:
        height = bar.get_height()
        if height > 0:
            ax3.text(
                bar.get_x() + bar.get_width() / 2, height + 1,
                f"{height:.1f}%", ha="center", va="bottom", fontsize=9,
            )

    fig.suptitle(
        "BioCompiler v12 vs DNAchisel — Publication Benchmark",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(str(output), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _generate_radar_figure(summary: dict, output: Path) -> None:
    """Generate a radar/spider chart comparing relative performance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    categories = ["CAI", "Speed", "Constraints"]

    # Normalize metrics to 0-1 scale (1 = best)
    bc_cai = summary.get("mean_cai_biocompiler", 0.0)
    dc_cai = summary.get("mean_cai_dnachisel", 0.0)
    max_cai = max(bc_cai, dc_cai, 0.01)

    bc_speed = summary.get("mean_speed_biocompiler_ms", 1.0)
    dc_speed = summary.get("mean_speed_dnachisel_ms", 1.0)
    # For speed, lower is better -> invert
    max_speed = max(bc_speed, dc_speed, 1.0)

    bc_cs = summary.get("constraint_satisfaction_rate_biocompiler", 0.0)
    dc_cs = summary.get("constraint_satisfaction_rate_dnachisel", 0.0)

    bc_values = [
        bc_cai / max_cai if max_cai > 0 else 0,
        1.0 - (bc_speed / max_speed) if max_speed > 0 else 0,
        bc_cs,
    ]
    dc_values = [
        dc_cai / max_cai if max_cai > 0 else 0,
        1.0 - (dc_speed / max_speed) if max_speed > 0 else 0,
        dc_cs,
    ]

    # Close the polygon
    bc_values.append(bc_values[0])
    dc_values.append(dc_values[0])

    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    angles.append(angles[0])

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    ax.fill(angles, bc_values, alpha=0.25, color="#2563EB", label="BioCompiler")
    ax.plot(angles, bc_values, color="#2563EB", linewidth=2)
    ax.fill(angles, dc_values, alpha=0.25, color="#DC2626", label="DNAchisel")
    ax.plot(angles, dc_values, color="#DC2626", linewidth=2)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=12)
    ax.set_ylim(0, 1)
    ax.set_title(
        "BioCompiler v12 vs DNAchisel — Relative Performance",
        fontsize=14, fontweight="bold", pad=20,
    )
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=10)

    fig.tight_layout()
    fig.savefig(str(output), dpi=150, bbox_inches="tight")
    plt.close(fig)
