"""
BioCompiler Benchmarking — Head-to-Head Runner
================================================

Orchestrates head-to-head benchmark comparisons between BioCompiler and
DNAchisel across curated gene sets.  Provides a high-level ``BenchmarkRunner``
API that iterates over genes, invokes both optimizers, computes shared
metrics, determines winners, and generates reports.

When DNAchisel is not installed, the runner still executes BioCompiler and
records results with DNAchisel fields set to ``None`` — emitting a warning
instead of crashing.

Usage::

    from biocompiler.benchmarking.runner import BenchmarkRunner, BenchmarkReport
    from biocompiler.benchmarking.gene_sets import VACCINE_ANTIGEN_GENES

    runner = BenchmarkRunner(gene_sets=VACCINE_ANTIGEN_GENES)
    results = runner.run_all()
    print(BenchmarkReport.generate(results))
    BenchmarkReport.to_csv(results, "benchmark_results.csv")
    stats = BenchmarkReport.summary_stats(results)

Design:
    All metrics are computed using BioCompiler's own evaluators (via
    :func:`biocompiler.benchmarking.metrics.compute_all_metrics`) for
    fairness — both tools' outputs are evaluated with the same CAI
    computation, GC measurement, and restriction site scanner.
"""

from __future__ import annotations

import csv
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .metrics import BenchmarkMetrics, compute_all_metrics
from .dnachisel_adapter import is_dnachisel_available

logger = logging.getLogger(__name__)

__all__ = [
    "BenchmarkResult",
    "BenchmarkRunner",
    "BenchmarkReport",
    "DEFAULT_ENZYME_PANEL",
    "AVAILABLE_BENCHMARKS",
    "run_benchmark_by_name",
]

# ---------------------------------------------------------------------------
# Default enzyme panel
# ---------------------------------------------------------------------------

DEFAULT_ENZYME_PANEL: list[str] = [
    "EcoRI", "BamHI", "HindIII", "XhoI", "NotI", "SalI", "PstI", "KpnI",
]
"""Standard 8-enzyme panel for benchmarking restriction-site avoidance."""

# GC range defaults for constraint checking
_DEFAULT_GC_LO: float = 0.30
_DEFAULT_GC_HI: float = 0.70

# ---------------------------------------------------------------------------
# Available benchmark registry
# ---------------------------------------------------------------------------

AVAILABLE_BENCHMARKS: dict[str, str] = {
    "sharp_li_cai": (
        "Compare CAI computed with Kazusa vs Sharp-Li reference sets "
        "against published values from Sharp & Li (1987) and Puigbo et al. (2008)"
    ),
    "organism_aware_cai": (
        "Demonstrate CAI recovery when eukaryotic-specific constraints "
        "(splice/CpG) are disabled for prokaryotic targets such as E. coli"
    ),
}
"""Registry of available named benchmarks with descriptions.

Keys are benchmark names that can be passed to :func:`run_benchmark_by_name`.
Values are human-readable descriptions of what each benchmark does.
"""

# Epsilon for floating-point comparisons
_EPSILON: float = 0.001

# Organism name normalisation (gene_sets use "Homo sapiens" with space;
# the optimizer expects "Homo_sapiens" with underscore).
_ORGANISM_ALIASES: dict[str, str] = {
    "homo sapiens": "Homo_sapiens",
    "escherichia coli": "Escherichia_coli",
    "mus musculus": "Mus_musculus",
    "saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
    "cricetulus griseus": "CHO_K1",
    "cho k1": "CHO_K1",
}


def _normalise_organism(organism: str) -> str:
    """Normalise organism name to the format expected by the optimizer.

    Handles both "Homo sapiens" (space) and "Homo_sapiens" (underscore)
    formats, as well as common lowercase aliases.
    """
    # Already in correct format?
    from ..organisms import SUPPORTED_ORGANISMS

    if organism in SUPPORTED_ORGANISMS:
        return organism
    # Try alias lookup (case-insensitive)
    alias = _ORGANISM_ALIASES.get(organism.lower())
    if alias is not None:
        return alias
    # Try replacing spaces with underscores
    candidate = organism.replace(" ", "_")
    if candidate in SUPPORTED_ORGANISMS:
        return candidate
    # Fallback — return as-is; the optimizer will raise if unsupported
    return candidate


# ---------------------------------------------------------------------------
# BenchmarkResult dataclass
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkResult:
    """Per-gene benchmark result comparing BioCompiler vs DNAchisel.

    Captures optimized sequences, computed metrics, runtimes, constraint
    satisfaction, and the determined winner for a single gene.

    Attributes
    ----------
    gene_name : str
        Gene symbol (e.g. ``"INS"``, ``"SARS2_S"``).
    organism : str
        Target organism name (normalised, e.g. ``"Homo_sapiens"``).
    protein_length : int
        Length of the protein in amino acids.
    biocompiler_sequence : str
        Optimized DNA sequence from BioCompiler.
    dnachisel_sequence : str | None
        Optimized DNA sequence from DNAchisel, or ``None`` if unavailable.
    biocompiler_metrics : BenchmarkMetrics
        Computed metrics for the BioCompiler output.
    dnachisel_metrics : BenchmarkMetrics | None
        Computed metrics for the DNAchisel output, or ``None`` if unavailable.
    biocompiler_runtime_s : float
        Wall-clock optimisation time for BioCompiler (seconds).
    dnachisel_runtime_s : float | None
        Wall-clock optimisation time for DNAchisel (seconds), or ``None``
        if unavailable.
    constraint_satisfaction : dict[str, bool]
        Constraint satisfaction flags for the BioCompiler output.
        Keys include ``"gc_in_range"``, ``"no_restriction_sites"``, etc.
    winner : str
        One of ``"biocompiler"``, ``"dnachisel"``, ``"tie"``, or
        ``"dnachisel_unavailable"``.
    """

    gene_name: str
    organism: str
    protein_length: int
    biocompiler_sequence: str
    dnachisel_sequence: str | None
    biocompiler_metrics: BenchmarkMetrics
    dnachisel_metrics: BenchmarkMetrics | None
    biocompiler_runtime_s: float
    dnachisel_runtime_s: float | None
    constraint_satisfaction: dict[str, bool] = field(default_factory=dict)
    winner: str = "dnachisel_unavailable"

    # ── Convenience properties ──

    @property
    def dnachisel_available(self) -> bool:
        """Whether DNAchisel results are available for this gene."""
        return self.dnachisel_metrics is not None

    @property
    def all_constraints_satisfied(self) -> bool:
        """Whether all evaluated constraints are satisfied."""
        return all(self.constraint_satisfaction.values())


# ---------------------------------------------------------------------------
# Constraint evaluation helper
# ---------------------------------------------------------------------------

def _evaluate_constraint_satisfaction(
    metrics: BenchmarkMetrics,
    gc_lo: float = _DEFAULT_GC_LO,
    gc_hi: float = _DEFAULT_GC_HI,
) -> dict[str, bool]:
    """Evaluate constraint satisfaction from a BenchmarkMetrics object.

    Returns a dict mapping constraint name → bool (True = satisfied).
    """
    return {
        "gc_in_range": gc_lo <= metrics.gc_profile.mean <= gc_hi,
        "no_restriction_sites": metrics.restriction_site_total == 0,
        "low_cryptic_splice_sites": metrics.cryptic_splice_sites == 0,
        "no_cpg_islands": metrics.cpg_islands == 0,
        "cai_above_threshold": metrics.cai >= 0.5,
        "good_mrna_stability": metrics.mrna_stability >= 0.5,
    }


# ---------------------------------------------------------------------------
# Winner determination
# ---------------------------------------------------------------------------

def _determine_winner(
    bc_metrics: BenchmarkMetrics,
    dc_metrics: BenchmarkMetrics | None,
) -> str:
    """Determine the winner based on composite metric comparison.

    Scoring: each metric category awards one point to the tool that is
    better (or zero for a tie).  The tool with the most points wins.

    Categories:
      1. CAI (higher is better)
      2. GC in range + closest to 0.50 midpoint
      3. Fewer restriction sites
      4. Fewer cryptic splice sites
      5. Higher mRNA stability

    Returns ``"biocompiler"``, ``"dnachisel"``, ``"tie"``, or
    ``"dnachisel_unavailable"``.
    """
    if dc_metrics is None:
        return "dnachisel_unavailable"

    bc_score = 0
    dc_score = 0

    # 1. CAI (higher is better)
    if bc_metrics.cai > dc_metrics.cai + _EPSILON:
        bc_score += 1
    elif dc_metrics.cai > bc_metrics.cai + _EPSILON:
        dc_score += 1

    # 2. GC closeness to 0.50 midpoint (closer is better)
    _GC_MID = 0.50
    bc_gc_diff = abs(bc_metrics.gc_profile.mean - _GC_MID)
    dc_gc_diff = abs(dc_metrics.gc_profile.mean - _GC_MID)
    if bc_gc_diff < dc_gc_diff - _EPSILON:
        bc_score += 1
    elif dc_gc_diff < bc_gc_diff - _EPSILON:
        dc_score += 1

    # 3. Restriction sites (fewer is better)
    if bc_metrics.restriction_site_total < dc_metrics.restriction_site_total:
        bc_score += 1
    elif dc_metrics.restriction_site_total < bc_metrics.restriction_site_total:
        dc_score += 1

    # 4. Cryptic splice sites (fewer is better)
    if bc_metrics.cryptic_splice_sites < dc_metrics.cryptic_splice_sites:
        bc_score += 1
    elif dc_metrics.cryptic_splice_sites < bc_metrics.cryptic_splice_sites:
        dc_score += 1

    # 5. mRNA stability (higher is better)
    if bc_metrics.mrna_stability > dc_metrics.mrna_stability + _EPSILON:
        bc_score += 1
    elif dc_metrics.mrna_stability > bc_metrics.mrna_stability + _EPSILON:
        dc_score += 1

    if bc_score > dc_score:
        return "biocompiler"
    if dc_score > bc_score:
        return "dnachisel"
    return "tie"


# ---------------------------------------------------------------------------
# BenchmarkRunner
# ---------------------------------------------------------------------------

class BenchmarkRunner:
    """Orchestrates head-to-head benchmark comparisons.

    Runs both BioCompiler and DNAchisel on each gene in the provided
    gene set, computes shared metrics, and returns ``BenchmarkResult``
    objects.

    Parameters
    ----------
    gene_sets : dict[str, dict]
        Mapping from gene name to gene data dict.  Each dict must contain
        at least ``"protein_sequence"`` and ``"organism"`` keys.
    enzymes : list[str] | None
        Restriction enzymes to avoid and to scan for.  Defaults to
        :data:`DEFAULT_ENZYME_PANEL`.

    Examples
    --------
    >>> from biocompiler.benchmarking.runner import BenchmarkRunner
    >>> from biocompiler.benchmarking.gene_sets import VACCINE_ANTIGEN_GENES
    >>> runner = BenchmarkRunner(gene_sets=VACCINE_ANTIGEN_GENES)
    >>> results = runner.run_all()
    """

    def __init__(
        self,
        gene_sets: dict[str, dict],
        enzymes: list[str] | None = None,
    ) -> None:
        self.gene_sets = gene_sets
        self.enzymes = enzymes if enzymes is not None else list(DEFAULT_ENZYME_PANEL)

        # Check DNAchisel availability once at init
        self._dnachisel_available = is_dnachisel_available()
        if not self._dnachisel_available:
            logger.warning(
                "DNAchisel is not installed — benchmark will run "
                "BioCompiler only.  Install with: pip install dnachisel"
            )

    # ── Single gene ──────────────────────────────────────────────────

    def run_gene(self, gene_name: str, gene_data: dict) -> BenchmarkResult:
        """Run both BioCompiler and DNAchisel on a single gene.

        Parameters
        ----------
        gene_name : str
            Gene symbol (e.g. ``"INS"``).
        gene_data : dict
            Dict with keys ``"protein_sequence"`` and ``"organism"``.

        Returns
        -------
        BenchmarkResult
            Comparison result for this gene.  DNAchisel fields will be
            ``None`` if DNAchisel is unavailable.
        """
        protein = gene_data["protein_sequence"]
        raw_organism = gene_data.get("organism", "Homo_sapiens")
        organism = _normalise_organism(raw_organism)

        logger.info("Benchmarking gene %s (%d aa, %s)", gene_name, len(protein), organism)

        # ── BioCompiler ──
        bc_seq: str = ""
        bc_metrics: BenchmarkMetrics | None = None
        bc_runtime: float = 0.0

        t0 = time.perf_counter()
        try:
            from ..optimization import optimize_sequence

            bc_result = optimize_sequence(
                target_protein=protein,
                organism=organism,
                gc_lo=_DEFAULT_GC_LO,
                gc_hi=_DEFAULT_GC_HI,
            )
            bc_seq = bc_result.sequence
            bc_metrics = compute_all_metrics(
                dna=bc_seq,
                protein=protein,
                organism=organism,
                enzymes=self.enzymes,
            )
        except Exception as exc:
            logger.error("BioCompiler failed for %s: %s", gene_name, exc)
            # Create a zeroed-metrics placeholder so we can still return
            bc_metrics = self._placeholder_metrics()
        bc_runtime = time.perf_counter() - t0

        # ── DNAchisel ──
        dc_seq: str | None = None
        dc_metrics: BenchmarkMetrics | None = None
        dc_runtime: float | None = None

        if self._dnachisel_available:
            t1 = time.perf_counter()
            try:
                from .dnachisel_adapter import DNAchiselAdapter

                adapter = DNAchiselAdapter()
                dc_result = adapter.optimize(
                    protein=protein,
                    organism=organism,
                    constraints=[
                        {"type": "gc_range", "gc_lo": _DEFAULT_GC_LO, "gc_hi": _DEFAULT_GC_HI},
                        {"type": "avoid_restriction", "enzymes": self.enzymes},
                    ],
                )
                if dc_result.success:
                    dc_seq = dc_result.sequence
                    dc_metrics = compute_all_metrics(
                        dna=dc_seq,
                        protein=protein,
                        organism=organism,
                        enzymes=self.enzymes,
                    )
                else:
                    logger.warning(
                        "DNAchisel optimization failed for %s: %s",
                        gene_name,
                        dc_result.error,
                    )
            except ImportError:
                logger.warning(
                    "DNAchisel import failed for %s — marking as unavailable",
                    gene_name,
                )
            except Exception as exc:
                logger.error("DNAchisel failed for %s: %s", gene_name, exc)
            dc_runtime = time.perf_counter() - t1

        # ── Assemble result ──
        assert bc_metrics is not None  # guaranteed by logic above

        constraint_satisfaction = _evaluate_constraint_satisfaction(
            bc_metrics, _DEFAULT_GC_LO, _DEFAULT_GC_HI,
        )
        winner = _determine_winner(bc_metrics, dc_metrics)

        return BenchmarkResult(
            gene_name=gene_name,
            organism=organism,
            protein_length=len(protein),
            biocompiler_sequence=bc_seq,
            dnachisel_sequence=dc_seq,
            biocompiler_metrics=bc_metrics,
            dnachisel_metrics=dc_metrics,
            biocompiler_runtime_s=bc_runtime,
            dnachisel_runtime_s=dc_runtime,
            constraint_satisfaction=constraint_satisfaction,
            winner=winner,
        )

    # ── All genes ────────────────────────────────────────────────────

    def run_all(self) -> list[BenchmarkResult]:
        """Run benchmarks for all genes in the gene set.

        Returns
        -------
        list[BenchmarkResult]
            Results for every gene, in iteration order.
        """
        results: list[BenchmarkResult] = []
        for gene_name, gene_data in self.gene_sets.items():
            result = self.run_gene(gene_name, gene_data)
            results.append(result)
        return results

    # ── Filter by organism ───────────────────────────────────────────

    def run_organism(self, organism: str) -> list[BenchmarkResult]:
        """Run benchmarks for genes matching a specific organism.

        Organism matching is case-insensitive and normalises spaces to
        underscores (e.g. ``"Homo sapiens"`` matches ``"Homo_sapiens"``).

        Parameters
        ----------
        organism : str
            Organism name to filter by.

        Returns
        -------
        list[BenchmarkResult]
            Results for genes from the specified organism.
        """
        normalised = _normalise_organism(organism)
        results: list[BenchmarkResult] = []
        for gene_name, gene_data in self.gene_sets.items():
            gene_org = _normalise_organism(gene_data.get("organism", ""))
            if gene_org == normalised:
                result = self.run_gene(gene_name, gene_data)
                results.append(result)
        return results

    # ── Placeholder metrics ──────────────────────────────────────────

    @staticmethod
    def _placeholder_metrics() -> BenchmarkMetrics:
        """Create a zeroed-out BenchmarkMetrics for error cases."""
        from .metrics import GCProfile

        return BenchmarkMetrics(
            cai=0.0,
            gc_profile=GCProfile(
                mean=0.0, std=0.0, min_=0.0, max_=0.0, window=50,
            ),
            restriction_sites={},
            restriction_site_total=0,
            cryptic_splice_sites=0,
            cpg_islands=0,
            codon_pair_bias=0.0,
            mrna_stability=0.0,
            sequence_identity=0.0,
            runtime_s=0.0,
        )


# ---------------------------------------------------------------------------
# run_benchmark_by_name — dispatch to the correct benchmark
# ---------------------------------------------------------------------------

def run_benchmark_by_name(name: str) -> dict:
    """Run a benchmark suite by name and return its results.

    Dispatches to the appropriate benchmark function based on *name*.
    This is the runner-level entry point that complements the higher-level
    :func:`biocompiler.benchmarking.run_benchmark_by_name` (which imports
    from this module and from ``sharp_li_benchmark`` / ``organism_aware_benchmark``).

    Parameters
    ----------
    name : str
        Benchmark name.  Must be a key in :data:`AVAILABLE_BENCHMARKS`.

    Returns
    -------
    dict
        Benchmark results (structure depends on the benchmark).

    Raises
    ------
    ValueError
        If *name* is not a recognised benchmark.

    Examples
    --------
    >>> from biocompiler.benchmarking.runner import run_benchmark_by_name
    >>> result = run_benchmark_by_name("sharp_li_cai")
    >>> print(result["sharp_li_is_closer"])
    True
    """
    if name not in AVAILABLE_BENCHMARKS:
        available = ", ".join(sorted(AVAILABLE_BENCHMARKS.keys()))
        raise ValueError(
            f"Unknown benchmark '{name}'. Available: {available}"
        )

    if name == "sharp_li_cai":
        from .sharp_li_benchmark import benchmark_sharp_li_cai
        return benchmark_sharp_li_cai()
    elif name == "organism_aware_cai":
        from .organism_aware_benchmark import benchmark_organism_aware_cai
        return benchmark_organism_aware_cai()
    else:
        # Should not reach here if AVAILABLE_BENCHMARKS is in sync
        raise ValueError(f"Benchmark '{name}' is registered but has no runner.")


# ---------------------------------------------------------------------------
# BenchmarkReport
# ---------------------------------------------------------------------------

class BenchmarkReport:
    """Static methods for generating human-readable and CSV benchmark reports.

    All methods are class methods (no instantiation needed).
    """

    @staticmethod
    def generate(results: list[BenchmarkResult]) -> str:
        """Generate a human-readable summary report.

        Parameters
        ----------
        results : list[BenchmarkResult]
            Benchmark results to summarise.

        Returns
        -------
        str
            Multi-line text report suitable for terminal output or log files.
        """
        if not results:
            return "No benchmark results to report.\n"

        lines: list[str] = []
        lines.append("=" * 72)
        lines.append("  BioCompiler Benchmark Report")
        lines.append("=" * 72)
        lines.append("")

        # Per-gene details
        for r in results:
            lines.append(f"Gene: {r.gene_name}  |  Organism: {r.organism}  |  "
                         f"Protein: {r.protein_length} aa")
            lines.append("-" * 72)

            # BioCompiler metrics
            bc = r.biocompiler_metrics
            lines.append(f"  BioCompiler  CAI={bc.cai:.4f}  "
                         f"GC={bc.gc_profile.mean:.3f}  "
                         f"RS={bc.restriction_site_total}  "
                         f"Splice={bc.cryptic_splice_sites}  "
                         f"CpG={bc.cpg_islands}  "
                         f"mRNA={bc.mrna_stability:.3f}  "
                         f"Time={r.biocompiler_runtime_s:.3f}s")

            # DNAchisel metrics
            if r.dnachisel_metrics is not None:
                dc = r.dnachisel_metrics
                dc_time = r.dnachisel_runtime_s if r.dnachisel_runtime_s is not None else 0.0
                lines.append(f"  DNAchisel    CAI={dc.cai:.4f}  "
                             f"GC={dc.gc_profile.mean:.3f}  "
                             f"RS={dc.restriction_site_total}  "
                             f"Splice={dc.cryptic_splice_sites}  "
                             f"CpG={dc.cpg_islands}  "
                             f"mRNA={dc.mrna_stability:.3f}  "
                             f"Time={dc_time:.3f}s")
            else:
                lines.append("  DNAchisel    [unavailable]")

            # Constraint satisfaction
            sat_count = sum(1 for v in r.constraint_satisfaction.values() if v)
            total_count = len(r.constraint_satisfaction)
            lines.append(f"  Constraints: {sat_count}/{total_count} satisfied  "
                         f"Winner: {r.winner}")
            lines.append("")

        # Aggregate summary
        stats = BenchmarkReport.summary_stats(results)
        lines.append("=" * 72)
        lines.append("  Summary")
        lines.append("=" * 72)
        lines.append(f"  Total genes:            {stats['total_genes']}")
        lines.append(f"  DNAchisel available:     {stats['genes_with_dnachisel']}")
        lines.append(f"  BioCompiler wins:       {stats['biocompiler_wins']}")
        lines.append(f"  DNAchisel wins:         {stats['dnachisel_wins']}")
        lines.append(f"  Ties:                   {stats['ties']}")
        lines.append(f"  Mean CAI (BioCompiler): {stats['mean_cai_biocompiler']:.4f}")
        if stats['mean_cai_dnachisel'] is not None:
            lines.append(f"  Mean CAI (DNAchisel):   {stats['mean_cai_dnachisel']:.4f}")
        lines.append(f"  Mean GC  (BioCompiler): {stats['mean_gc_biocompiler']:.4f}")
        if stats['mean_gc_dnachisel'] is not None:
            lines.append(f"  Mean GC  (DNAchisel):   {stats['mean_gc_dnachisel']:.4f}")
        lines.append(f"  Mean Time (BioCompiler):{stats['mean_runtime_biocompiler']:.4f}s")
        if stats['mean_runtime_dnachisel'] is not None:
            lines.append(f"  Mean Time (DNAchisel):  {stats['mean_runtime_dnachisel']:.4f}s")
        lines.append(f"  Constraint satisfaction: {stats['constraint_satisfaction_rate']:.1%}")
        lines.append("")

        return "\n".join(lines)

    @staticmethod
    def to_csv(
        results: list[BenchmarkResult],
        filepath: str | Path,
    ) -> None:
        """Write benchmark results to a CSV file.

        The CSV contains one row per gene with columns for gene name,
        organism, protein length, metrics from both tools, and the winner.

        Parameters
        ----------
        results : list[BenchmarkResult]
            Benchmark results to export.
        filepath : str | Path
            Output file path.  Parent directories are created if needed.
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        fieldnames = [
            "gene_name",
            "organism",
            "protein_length",
            "biocompiler_cai",
            "biocompiler_gc",
            "biocompiler_restriction_sites",
            "biocompiler_cryptic_splice_sites",
            "biocompiler_cpg_islands",
            "biocompiler_mrna_stability",
            "biocompiler_runtime_s",
            "dnachisel_cai",
            "dnachisel_gc",
            "dnachisel_restriction_sites",
            "dnachisel_cryptic_splice_sites",
            "dnachisel_cpg_islands",
            "dnachisel_mrna_stability",
            "dnachisel_runtime_s",
            "gc_in_range",
            "no_restriction_sites",
            "low_cryptic_splice_sites",
            "no_cpg_islands",
            "cai_above_threshold",
            "good_mrna_stability",
            "winner",
        ]

        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for r in results:
                bc = r.biocompiler_metrics
                dc = r.dnachisel_metrics

                row: dict[str, Any] = {
                    "gene_name": r.gene_name,
                    "organism": r.organism,
                    "protein_length": r.protein_length,
                    # BioCompiler metrics
                    "biocompiler_cai": f"{bc.cai:.4f}",
                    "biocompiler_gc": f"{bc.gc_profile.mean:.4f}",
                    "biocompiler_restriction_sites": bc.restriction_site_total,
                    "biocompiler_cryptic_splice_sites": bc.cryptic_splice_sites,
                    "biocompiler_cpg_islands": bc.cpg_islands,
                    "biocompiler_mrna_stability": f"{bc.mrna_stability:.4f}",
                    "biocompiler_runtime_s": f"{r.biocompiler_runtime_s:.4f}",
                    # DNAchisel metrics
                    "dnachisel_cai": f"{dc.cai:.4f}" if dc is not None else "",
                    "dnachisel_gc": f"{dc.gc_profile.mean:.4f}" if dc is not None else "",
                    "dnachisel_restriction_sites": dc.restriction_site_total if dc is not None else "",
                    "dnachisel_cryptic_splice_sites": dc.cryptic_splice_sites if dc is not None else "",
                    "dnachisel_cpg_islands": dc.cpg_islands if dc is not None else "",
                    "dnachisel_mrna_stability": f"{dc.mrna_stability:.4f}" if dc is not None else "",
                    "dnachisel_runtime_s": (f"{r.dnachisel_runtime_s:.4f}"
                                            if r.dnachisel_runtime_s is not None else ""),
                    # Constraints
                    "gc_in_range": r.constraint_satisfaction.get("gc_in_range", ""),
                    "no_restriction_sites": r.constraint_satisfaction.get("no_restriction_sites", ""),
                    "low_cryptic_splice_sites": r.constraint_satisfaction.get("low_cryptic_splice_sites", ""),
                    "no_cpg_islands": r.constraint_satisfaction.get("no_cpg_islands", ""),
                    "cai_above_threshold": r.constraint_satisfaction.get("cai_above_threshold", ""),
                    "good_mrna_stability": r.constraint_satisfaction.get("good_mrna_stability", ""),
                    # Winner
                    "winner": r.winner,
                }
                writer.writerow(row)

        logger.info("Benchmark results written to %s", filepath)

    @staticmethod
    def summary_stats(results: list[BenchmarkResult]) -> dict[str, Any]:
        """Compute aggregate statistics across benchmark results.

        Parameters
        ----------
        results : list[BenchmarkResult]
            Benchmark results to aggregate.

        Returns
        -------
        dict
            Aggregate statistics including means, win counts, and
            constraint satisfaction rates.
        """
        if not results:
            return {
                "total_genes": 0,
                "genes_with_dnachisel": 0,
                "biocompiler_wins": 0,
                "dnachisel_wins": 0,
                "ties": 0,
                "mean_cai_biocompiler": 0.0,
                "mean_cai_dnachisel": None,
                "mean_gc_biocompiler": 0.0,
                "mean_gc_dnachisel": None,
                "mean_runtime_biocompiler": 0.0,
                "mean_runtime_dnachisel": None,
                "constraint_satisfaction_rate": 0.0,
            }

        n = len(results)

        # Win counts
        biocompiler_wins = sum(1 for r in results if r.winner == "biocompiler")
        dnachisel_wins = sum(1 for r in results if r.winner == "dnachisel")
        ties = sum(1 for r in results if r.winner == "tie")

        # DNAchisel availability
        genes_with_dc = sum(1 for r in results if r.dnachisel_available)

        # Mean CAI
        mean_cai_bc = sum(r.biocompiler_metrics.cai for r in results) / n
        dc_cais = [r.dnachisel_metrics.cai for r in results if r.dnachisel_metrics is not None]
        mean_cai_dc = sum(dc_cais) / len(dc_cais) if dc_cais else None

        # Mean GC
        mean_gc_bc = sum(r.biocompiler_metrics.gc_profile.mean for r in results) / n
        dc_gcs = [r.dnachisel_metrics.gc_profile.mean for r in results if r.dnachisel_metrics is not None]
        mean_gc_dc = sum(dc_gcs) / len(dc_gcs) if dc_gcs else None

        # Mean runtime
        mean_rt_bc = sum(r.biocompiler_runtime_s for r in results) / n
        dc_rts = [r.dnachisel_runtime_s for r in results if r.dnachisel_runtime_s is not None]
        mean_rt_dc = sum(dc_rts) / len(dc_rts) if dc_rts else None

        # Constraint satisfaction rate
        total_constraints = sum(len(r.constraint_satisfaction) for r in results)
        satisfied_constraints = sum(
            sum(1 for v in r.constraint_satisfaction.values() if v)
            for r in results
        )
        constraint_rate = (satisfied_constraints / total_constraints
                          if total_constraints > 0 else 0.0)

        return {
            "total_genes": n,
            "genes_with_dnachisel": genes_with_dc,
            "biocompiler_wins": biocompiler_wins,
            "dnachisel_wins": dnachisel_wins,
            "ties": ties,
            "mean_cai_biocompiler": round(mean_cai_bc, 4),
            "mean_cai_dnachisel": round(mean_cai_dc, 4) if mean_cai_dc is not None else None,
            "mean_gc_biocompiler": round(mean_gc_bc, 4),
            "mean_gc_dnachisel": round(mean_gc_dc, 4) if mean_gc_dc is not None else None,
            "mean_runtime_biocompiler": round(mean_rt_bc, 4),
            "mean_runtime_dnachisel": round(mean_rt_dc, 4) if mean_rt_dc is not None else None,
            "constraint_satisfaction_rate": constraint_rate,
        }
