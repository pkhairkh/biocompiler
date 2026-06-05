"""
BioCompiler Benchmarking Sub-Package
======================================

Gene sets, cross-validation tools, and comparison adapters for BioCompiler
performance evaluation.

This sub-package provides the infrastructure needed to benchmark BioCompiler
against alternative tools (e.g., DNAchisel) and to validate its metrics
against published implementations.

Provides:
  - BenchmarkRunner: Convenience class for running benchmark comparisons
  - BenchmarkResult: Per-gene comparison dataclass
  - BenchmarkReport: Aggregate benchmark report
  - ComparisonSummary: Aggregate statistical analysis
  - compute_cai_sharp_li: Independent Sharp & Li (1987) CAI implementation
  - DNAchiselAdapter: Wraps DNAchisel's API for comparable results
  - run_comparison: Single-gene head-to-head comparison
  - compare_results: Multi-gene statistical analysis
  - benchmark_sharp_li_cai: Sharp-Li vs Kazusa CAI reference set benchmark
  - benchmark_organism_aware_cai: Organism-aware constraint selection benchmark
  - BENCHMARK_SUITES: Dict mapping benchmark names to runner functions
  - HUMAN_THERAPEUTIC_GENES: 10 human therapeutic proteins for real-world benchmarking
  - VACCINE_ANTIGEN_GENES: 10 vaccine antigen proteins for real-world benchmarking
  - STRESS_TEST_GENES: 5 synthetic stress-test sequences for edge-case evaluation
  - get_all_gene_sets: Returns all gene sets merged

Usage::

    from biocompiler.benchmarking import BenchmarkRunner, BenchmarkResult

    # Run benchmarks with the convenience runner
    runner = BenchmarkRunner(organism="Homo_sapiens")
    report = runner.run_genes(["HBB", "INS"])

    # Or use the functional API
    from biocompiler.benchmarking import run_comparison, compare_results
    result = run_comparison("HBB", "MVLSPADKTNVKAAWGKVGA", "Homo_sapiens")
    summary = compare_results([result])

    # Run a named benchmark suite
    from biocompiler.benchmarking import run_benchmark_by_name
    result = run_benchmark_by_name("sharp_li_cai")
"""

from __future__ import annotations

import logging
import time as _time
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

from .visualization import (  # noqa: F401
    plot_cai_comparison,
    plot_gc_comparison,
    plot_runtime_comparison,
    plot_constraint_satisfaction,
    plot_summary_dashboard,
    generate_latex_table,
)

from .metrics import compute_cai_validated  # noqa: F401

__all__ = [
    # Runner & report
    "BenchmarkRunner",
    "BenchmarkReport",
    # Comparison types
    "BenchmarkResult",
    "ComparisonSummary",
    "run_comparison",
    "compare_results",
    # Validated CAI
    "compute_cai_sharp_li",
    "compute_cai_validated",
    # DNAchisel adapter
    "DNAchiselAdapter",
    # Reference-set benchmarks
    "benchmark_sharp_li_cai",
    "benchmark_organism_aware_cai",
    "BENCHMARK_SUITES",
    "run_benchmark_by_name",
    # Visualization
    "plot_cai_comparison",
    "plot_gc_comparison",
    "plot_runtime_comparison",
    "plot_constraint_satisfaction",
    "plot_summary_dashboard",
    "generate_latex_table",
]

# ---------------------------------------------------------------------------
# Re-export from sub-modules
# ---------------------------------------------------------------------------

from .comparison import (
    BenchmarkResult,
    ComparisonSummary,
    run_comparison,
    compare_results,
)

from .cai_validated import compute_cai_sharp_li

from .dnachisel_adapter import DNAchiselAdapter

from .sharp_li_benchmark import benchmark_sharp_li_cai

from .organism_aware_benchmark import benchmark_organism_aware_cai


# ---------------------------------------------------------------------------
# BENCHMARK_SUITES — registry of named benchmark suites
# ---------------------------------------------------------------------------

BENCHMARK_SUITES: dict[str, callable] = {
    "sharp_li_cai": benchmark_sharp_li_cai,
    "organism_aware_cai": benchmark_organism_aware_cai,
}
"""Registry mapping benchmark names to their runner functions.

Each value is a callable that runs the benchmark and returns a dict of
results.  Use :func:`run_benchmark_by_name` for convenient dispatch.

Available benchmarks:

  - ``"sharp_li_cai"`` — Compares CAI computed with the Kazusa vs Sharp-Li
    reference sets against published values from Sharp & Li (1987) and
    Puigbo et al. (2008).  Validates both the algorithm implementation and
    the reference-set provenance.
  - ``"organism_aware_cai"`` — Demonstrates CAI recovery when
    eukaryotic-specific constraints (cryptic splice-site avoidance, CpG-island
    avoidance) are disabled for prokaryotic targets such as *E. coli*.
"""


def run_benchmark_by_name(name: str) -> dict:
    """Run a benchmark suite by name.

    Convenience function that looks up the benchmark in
    :data:`BENCHMARK_SUITES` and executes it, returning the results dict.

    Args:
        name: Benchmark name. Must be a key in :data:`BENCHMARK_SUITES`.

    Returns:
        Dict of benchmark results (structure depends on the benchmark).

    Raises:
        ValueError: If ``name`` is not a recognised benchmark.

    Examples::

        from biocompiler.benchmarking import run_benchmark_by_name

        # Run the Sharp-Li vs Kazusa CAI benchmark
        result = run_benchmark_by_name("sharp_li_cai")
        print(f"Sharp-Li closer: {result['sharp_li_is_closer']}")

        # Run the organism-aware constraint benchmark
        result = run_benchmark_by_name("organism_aware_cai")
        print(f"Mean CAI recovery: {result['mean_cai_recovery']:+.4f}")
    """
    if name not in BENCHMARK_SUITES:
        available = ", ".join(sorted(BENCHMARK_SUITES.keys()))
        raise ValueError(
            f"Unknown benchmark '{name}'. Available: {available}"
        )
    return BENCHMARK_SUITES[name]()

# Lazy imports for gene sets — module may not exist in all installations
try:
    from .gene_sets import (
        HUMAN_THERAPEUTIC_GENES,
        VACCINE_ANTIGEN_GENES,
        STRESS_TEST_GENES,
        get_all_gene_sets,
    )
    __all__ += [
        "HUMAN_THERAPEUTIC_GENES",
        "VACCINE_ANTIGEN_GENES",
        "STRESS_TEST_GENES",
        "get_all_gene_sets",
    ]
except ImportError:
    pass

# Also try importing extended gene sets
try:
    from .gene_sets import (
        E_COLI_EXTENDED,
        HUMAN_THERAPEUTIC,
        HUMAN_SIGNALING,
        YEAST_INDUSTRIAL,
        MOUSE_MODEL,
    )
    __all__ += [
        "E_COLI_EXTENDED",
        "HUMAN_THERAPEUTIC",
        "HUMAN_SIGNALING",
        "YEAST_INDUSTRIAL",
        "MOUSE_MODEL",
    ]
except ImportError:
    pass


# ---------------------------------------------------------------------------
# BenchmarkReport — aggregate report for multiple benchmark runs
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkReport:
    """Aggregate benchmark report for a set of gene comparisons.

    Collects the results from running multiple genes through the comparison
    pipeline and provides summary statistics.

    Attributes:
        organism: Target organism used for all comparisons.
        gene_names: List of gene names that were benchmarked.
        results: List of BenchmarkResult objects from each gene comparison.
        total_time_s: Total wall-clock time for all comparisons.
        timestamp: ISO 8601 timestamp of when the report was generated.
        summary: ComparisonSummary (populated after calling compute_summary).
    """

    organism: str
    gene_names: list[str] = field(default_factory=list)
    results: list[BenchmarkResult] = field(default_factory=list)
    total_time_s: float = 0.0
    timestamp: str = ""
    summary: ComparisonSummary | None = None

    def compute_summary(self) -> ComparisonSummary:
        """Compute and store the aggregate statistical summary.

        Returns:
            ComparisonSummary with statistical analysis across all results.
        """
        self.summary = compare_results(self.results)
        return self.summary

    def to_dict(self) -> dict[str, Any]:
        """Serialize the report to a JSON-compatible dict."""
        return {
            "organism": self.organism,
            "gene_names": list(self.gene_names),
            "num_results": len(self.results),
            "total_time_s": self.total_time_s,
            "timestamp": self.timestamp,
            "summary": self.summary.__dict__ if self.summary else None,
        }


# ---------------------------------------------------------------------------
# BenchmarkRunner — convenience class for running benchmark comparisons
# ---------------------------------------------------------------------------

class BenchmarkRunner:
    """Convenience runner for benchmarking BioCompiler on gene sets.

    Wraps the functional API (``run_comparison``, ``compare_results``)
    in a reusable class that can be configured once and run multiple times.

    Example::

        runner = BenchmarkRunner(
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        report = runner.run_genes(["HBB", "INS", "EGFP"])
        print(report.summary)
    """

    def __init__(
        self,
        organism: str = "Homo_sapiens",
        enzymes: list[str] | None = None,
        gc_lo: float = 0.30,
        gc_hi: float = 0.70,
    ) -> None:
        """Initialize the benchmark runner.

        Args:
            organism: Target organism for codon usage.
            enzymes: Restriction enzymes to avoid. If None, uses defaults.
            gc_lo: Minimum GC content fraction.
            gc_hi: Maximum GC content fraction.
        """
        self.organism = organism
        self.enzymes = enzymes
        self.gc_lo = gc_lo
        self.gc_hi = gc_hi

    def run_gene(
        self,
        gene_name: str,
        protein: str,
    ) -> BenchmarkResult:
        """Run a head-to-head comparison for a single gene.

        Args:
            gene_name: Name of the gene (e.g., "HBB").
            protein: Amino acid sequence.

        Returns:
            BenchmarkResult with metrics from both optimizers.
        """
        return run_comparison(
            gene_name=gene_name,
            protein=protein,
            organism=self.organism,
            enzymes=self.enzymes,
            gc_lo=self.gc_lo,
            gc_hi=self.gc_hi,
        )

    def run_genes(
        self,
        gene_names: list[str],
        protein_sequences: dict[str, str] | None = None,
    ) -> BenchmarkReport:
        """Run head-to-head comparisons for multiple genes.

        If ``protein_sequences`` is provided, uses those sequences.
        Otherwise, attempts to look up genes from the built-in gene sets.

        Args:
            gene_names: List of gene names to benchmark.
            protein_sequences: Optional dict mapping gene name → protein
                sequence. If None, built-in gene sets are used.

        Returns:
            BenchmarkReport with all results and aggregate summary.
        """
        from datetime import datetime, timezone

        results: list[BenchmarkResult] = []
        t0 = _time.perf_counter()

        for gene_name in gene_names:
            protein = ""
            if protein_sequences and gene_name in protein_sequences:
                protein = protein_sequences[gene_name]
            else:
                # Try built-in gene sets
                protein = self._lookup_gene(gene_name)

            if not protein:
                logger.warning(
                    "No protein sequence found for gene '%s'; skipping",
                    gene_name,
                )
                continue

            result = self.run_gene(gene_name, protein)
            results.append(result)

        elapsed = _time.perf_counter() - t0

        report = BenchmarkReport(
            organism=self.organism,
            gene_names=list(gene_names),
            results=results,
            total_time_s=elapsed,
            timestamp=datetime.now(timezone.utc).isoformat(),
        )
        if results:
            report.compute_summary()

        return report

    @staticmethod
    def _lookup_gene(gene_name: str) -> str:
        """Look up a protein sequence from built-in gene sets.

        Args:
            gene_name: Gene name to look up.

        Returns:
            Protein sequence string, or empty string if not found.
        """
        try:
            from .gene_sets import HUMAN_THERAPEUTIC_GENES
            entry = HUMAN_THERAPEUTIC_GENES.get(gene_name)
            if entry:
                return entry.get("protein_sequence", "")
        except ImportError:
            pass

        # Try vaccine antigen gene set
        try:
            from .gene_sets import VACCINE_ANTIGEN_GENES
            entry = VACCINE_ANTIGEN_GENES.get(gene_name)
            if entry:
                return entry.get("protein_sequence", "")
        except ImportError:
            pass

        # Try stress-test gene set
        try:
            from .gene_sets import STRESS_TEST_GENES
            entry = STRESS_TEST_GENES.get(gene_name)
            if entry:
                return entry.get("protein_sequence", "")
        except ImportError:
            pass

        # Common vaccine antigen genes (fallback)
        _COMMON_GENES: dict[str, str] = {
            "INS": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
            "HBB": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
            "EGFP": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        }
        return _COMMON_GENES.get(gene_name, "")
