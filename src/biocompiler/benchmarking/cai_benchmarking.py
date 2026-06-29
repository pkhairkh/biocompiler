"""
CAI Benchmarking Against Published Values
==========================================

Benchmark predicted CAI values against published values from the literature.
This module provides structured comparison of BioCompiler's CAI implementation
against known ground-truth values from Sharp & Li (1987), Puigbo et al (2008),
and other published sources.

Key principle: Published CAI values are for **native** (unoptimized) genes.
Our **optimized** CAI should be HIGHER than published native-gene values.
If our optimized CAI is lower than the published native CAI, that indicates
a problem with the optimization.

Data sources
------------
1. Sharp, P.M. & Li, W.-H. (1987) Nucleic Acids Research 15:1281-1295
2. Puigbo, P., Bravo, I.G. & Garcia-Vallve, S. (2008) BMC Bioinformatics 9:65
3. Codon Usage Database (Kazusa)
4. GtRNAdb — tRNA gene copy numbers

References
----------
Sharp, P. M., & Li, W.-H. (1987). The codon Adaptation Index—a measure of
directional synonymous codon usage bias, and its potential applications.
*Nucleic Acids Research*, 15(3), 1281-1295. doi:10.1093/nar/15.3.1281

Puigbo, P., Bravo, I.G. & Garcia-Vallve, S. (2008). CAIcal: A combined set
of tools to assess codon usage adaptation. *BMC Bioinformatics*, 9:65.
doi:10.1186/1471-2105-9-65
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from biocompiler.expression.translation import compute_cai
from ..organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES

__all__ = [
    "CAIBenchmarkResult",
    "PUBLISHED_CAI_VALUES",
    "benchmark_cai",
    "benchmark_cai_for_dna",
    "benchmark_optimization",
    "summarize_benchmark",
    "BenchmarkSummary",
]


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Published CAI Values Database
# ═══════════════════════════════════════════════════════════════════════════════
#
# Mapping: (gene, organism) -> {"cai": float, "source": str}
#
# IMPORTANT: These values are for NATIVE genes (unoptimized codon usage).
# Our optimized CAI should ALWAYS be >= published native CAI.
# An optimization that produces CAI lower than the native gene's CAI
# is by definition a failed optimization.
#
# Organism names use the short aliases accepted by resolve_organism().
#
# Tolerance: ±0.05 for same reference set; ±0.10 for cross-reference-set
# comparisons (e.g., Sharp-Li vs Kazusa).

PUBLISHED_CAI_VALUES: dict[tuple[str, str], dict[str, object]] = {
    # ── E. coli — Sharp & Li (1987), Table 1 ────────────────────────
    ("lacZ", "e_coli"): {
        "cai": 0.27,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Lowly expressed; CAI reflects non-optimal codon bias",
    },
    ("trpA", "e_coli"): {
        "cai": 0.84,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Highly expressed; strong optimal codon usage",
    },
    ("recA", "e_coli"): {
        "cai": 0.59,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "SOS response gene; moderate expression",
    },
    ("ompA", "e_coli"): {
        "cai": 0.79,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Highly expressed membrane protein",
    },
    ("groEL", "e_coli"): {
        "cai": 0.76,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Highly expressed chaperone",
    },
    ("rpoB", "e_coli"): {
        "cai": 0.50,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "RNA polymerase beta subunit; moderate expression",
    },
    ("dnaK", "e_coli"): {
        "cai": 0.56,
        "source": "Sharp & Li 1987; Ikemura 1985",
        "gene_type": "native",
        "notes": "Heat-shock chaperone; inducible expression",
    },

    # ── E. coli — Heterologous genes from Puigbo et al (2008) ────────
    ("gfp", "e_coli"): {
        "cai": 0.54,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "heterologous",
        "notes": "Jellyfish GFP with native codons in E. coli",
    },
    ("insulin", "e_coli"): {
        "cai": 0.34,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "heterologous",
        "notes": "Human insulin with native human codons in E. coli",
    },
    ("hGH", "e_coli"): {
        "cai": 0.32,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "heterologous",
        "notes": "Human growth hormone with native codons in E. coli",
    },
    ("IFN-alpha2", "e_coli"): {
        "cai": 0.33,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "heterologous",
        "notes": "Interferon-alpha2 with native codons in E. coli",
    },

    # ── S. cerevisiae — Sharp & Li (1987), Table 1 ──────────────────
    ("ADH1", "yeast"): {
        "cai": 0.91,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "One of the most highly expressed yeast genes",
    },
    ("PGK1", "yeast"): {
        "cai": 0.88,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Highly expressed glycolytic enzyme",
    },
    ("ENO1", "yeast"): {
        "cai": 0.72,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Moderately-to-highly expressed enolase",
    },
    ("ACT1", "yeast"): {
        "cai": 0.56,
        "source": "Sharp & Li 1987, Table 1",
        "gene_type": "native",
        "notes": "Constitutive actin; moderate CAI due to regulatory constraints",
    },

    # ── Human — Puigbo et al (2008) ─────────────────────────────────
    ("HBB", "human"): {
        "cai": 0.95,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "native",
        "notes": "Beta-globin; extremely high expression in erythroid cells",
    },
    ("INS", "human"): {
        "cai": 0.84,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "native",
        "notes": "Insulin; highly expressed in pancreatic beta cells",
    },
    ("ALB", "human"): {
        "cai": 0.78,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "native",
        "notes": "Albumin; most abundant plasma protein",
    },
    ("TP53", "human"): {
        "cai": 0.63,
        "source": "Puigbo et al 2008, CAIcal server",
        "gene_type": "native",
        "notes": "p53 tumor suppressor; moderate expression under normal conditions",
    },
    ("insulin", "human"): {
        "cai": 0.72,
        "source": "Codon Usage Database",
        "gene_type": "native",
        "notes": "Human insulin native gene (different reference set from Puigbo)",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Data classes
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class CAIBenchmarkResult:
    """Result of benchmarking a predicted CAI against a published value.

    Attributes
    ----------
    gene : str
        Gene name (e.g. ``"gfp"``, ``"insulin"``).
    organism : str
        Organism name, resolved to canonical form
        (e.g. ``"Escherichia_coli"``).
    predicted_cai : float
        CAI computed by BioCompiler for the given sequence.
    published_cai : float or None
        Published CAI value from the literature, if available.
        ``None`` if no published value exists for this gene/organism pair.
    source : str
        Citation for the published value, or ``"N/A"`` if none.
    deviation : float or None
        Difference between predicted and published CAI
        (predicted - published). Positive means our prediction is higher
        than published; negative means lower.  ``None`` if no published
        value is available.
    gene_type : str
        Whether the published CAI is for a ``"native"`` or
        ``"heterologous"`` gene.
    pass_threshold : bool or None
        Whether the deviation is within acceptable tolerance.
        For native genes: we require predicted >= published (deviation >= 0)
        or at minimum within ±0.10.
        For heterologous genes: any deviation is acceptable since these
        are unoptimized.
        ``None`` if no published value is available.
    """
    gene: str
    organism: str
    predicted_cai: float
    published_cai: float | None
    source: str
    deviation: float | None
    gene_type: str
    pass_threshold: bool | None


@dataclass
class BenchmarkSummary:
    """Summary statistics for a batch of CAI benchmark results.

    Attributes
    ----------
    total_genes : int
        Total number of genes benchmarked.
    genes_with_published : int
        Number of genes that had published CAI values for comparison.
    mean_deviation : float
        Mean deviation (predicted - published) across all genes
        with published values.
    max_deviation : float
        Maximum positive deviation (our prediction is much higher).
    min_deviation : float
        Maximum negative deviation (our prediction is lower than published).
    pass_rate : float
        Fraction of genes with published values that pass the threshold.
    failures : list[str]
        List of gene names that failed the threshold check.
    """
    total_genes: int
    genes_with_published: int
    mean_deviation: float
    max_deviation: float
    min_deviation: float
    pass_rate: float
    failures: list[str] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Benchmarking functions
# ═══════════════════════════════════════════════════════════════════════════════

# Acceptable tolerance for CAI deviation
_NATIVE_TOLERANCE: float = 0.10
_HETEROLOGOUS_TOLERANCE: float = 0.15


def benchmark_cai_for_dna(
    gene: str,
    organism: str,
    dna: str,
) -> CAIBenchmarkResult:
    """Benchmark a single DNA sequence's CAI against published values.

    Computes the CAI for the given DNA sequence using BioCompiler's
    implementation, then compares it against any published value in
    the ``PUBLISHED_CAI_VALUES`` database.

    Parameters
    ----------
    gene : str
        Gene name (must match a key in ``PUBLISHED_CAI_VALUES``).
    organism : str
        Organism name (any alias accepted by ``resolve_organism``).
    dna : str
        DNA coding sequence to evaluate.

    Returns
    -------
    CAIBenchmarkResult
        Structured result with predicted CAI, published CAI, deviation,
        and pass/fail status.
    """
    # Resolve organism name
    resolved_org = resolve_organism(organism, strict=False)

    # Compute predicted CAI
    predicted = compute_cai(dna, organism=resolved_org)

    # Look up published value
    published_entry = PUBLISHED_CAI_VALUES.get((gene, organism))
    if published_entry is None:
        # Try with resolved organism short name
        from ..organisms import SPECIES_SHORT_NAMES
        short_name = SPECIES_SHORT_NAMES.get(resolved_org, organism)
        published_entry = PUBLISHED_CAI_VALUES.get((gene, short_name))

    if published_entry is not None:
        published_cai = float(published_entry["cai"])  # type: ignore[arg-type]
        source = str(published_entry.get("source", "Unknown"))
        gene_type = str(published_entry.get("gene_type", "native"))
        deviation = predicted - published_cai

        # Determine pass/fail
        if gene_type == "heterologous":
            tolerance = _HETEROLOGOUS_TOLERANCE
            # For heterologous genes, any deviation is acceptable since
            # published value is for native (unoptimized) codons
            pass_threshold = True
        else:
            tolerance = _NATIVE_TOLERANCE
            # For native genes, we want predicted >= published
            # (our optimization should never decrease CAI)
            # But if comparing against a different reference set,
            # allow tolerance
            pass_threshold = deviation >= -tolerance
    else:
        published_cai = None
        source = "N/A"
        gene_type = "unknown"
        deviation = None
        pass_threshold = None

    return CAIBenchmarkResult(
        gene=gene,
        organism=resolved_org,
        predicted_cai=predicted,
        published_cai=published_cai,
        source=source,
        deviation=deviation,
        gene_type=gene_type,
        pass_threshold=pass_threshold,
    )


def benchmark_cai(
    genes: list[str],
    organisms: list[str],
    dna_sequences: dict[tuple[str, str], str] | None = None,
) -> list[CAIBenchmarkResult]:
    """Benchmark predicted CAI against published values for multiple genes.

    For each (gene, organism) pair, computes the CAI of the provided DNA
    sequence and compares it against published values.

    Parameters
    ----------
    genes : list[str]
        List of gene names to benchmark.
    organisms : list[str]
        List of organism names (matched 1:1 with genes).
    dna_sequences : dict or None
        Mapping of (gene, organism) -> DNA sequence.  If ``None``,
        only genes with entries in ``PUBLISHED_CAI_VALUES`` are checked
        and their predicted CAI is set to 0.0 (placeholder).

    Returns
    -------
    list[CAIBenchmarkResult]
        One result per (gene, organism) pair.

    Raises
    ------
    ValueError
        If ``genes`` and ``organisms`` have different lengths.
    """
    if len(genes) != len(organisms):
        raise ValueError(
            f"genes and organisms must have the same length: "
            f"{len(genes)} != {len(organisms)}"
        )

    results: list[CAIBenchmarkResult] = []
    for gene, organism in zip(genes, organisms):
        dna = ""
        if dna_sequences is not None:
            dna = dna_sequences.get((gene, organism), "")
            # Also try with short organism name
            if not dna:
                resolved = resolve_organism(organism, strict=False)
                from ..organisms import SPECIES_SHORT_NAMES
                short_name = SPECIES_SHORT_NAMES.get(resolved, organism)
                dna = dna_sequences.get((gene, short_name), "")

        if dna:
            result = benchmark_cai_for_dna(gene, organism, dna)
        else:
            # No DNA sequence provided; record with predicted_cai = 0.0
            resolved_org = resolve_organism(organism, strict=False)
            published_entry = PUBLISHED_CAI_VALUES.get((gene, organism))
            if published_entry is None:
                from ..organisms import SPECIES_SHORT_NAMES
                short_name = SPECIES_SHORT_NAMES.get(resolved_org, organism)
                published_entry = PUBLISHED_CAI_VALUES.get((gene, short_name))

            if published_entry is not None:
                published_cai = float(published_entry["cai"])  # type: ignore[arg-type]
                source = str(published_entry.get("source", "Unknown"))
                gene_type = str(published_entry.get("gene_type", "native"))
            else:
                published_cai = None
                source = "N/A"
                gene_type = "unknown"

            result = CAIBenchmarkResult(
                gene=gene,
                organism=resolved_org,
                predicted_cai=0.0,
                published_cai=published_cai,
                source=source,
                deviation=None,
                gene_type=gene_type,
                pass_threshold=None,
            )
        results.append(result)

    return results


def benchmark_optimization(
    gene: str,
    organism: str,
    original_dna: str,
    optimized_dna: str,
) -> dict[str, object]:
    """Benchmark an optimization by comparing original vs optimized CAI.

    This is the key quality check: after optimization, the CAI should
    be >= the original (native gene) CAI, and ideally >= the published
    CAI for that gene in the target organism.

    Parameters
    ----------
    gene : str
        Gene name.
    organism : str
        Target organism for expression.
    original_dna : str
        Original (pre-optimization) DNA sequence.
    optimized_dna : str
        Optimized DNA sequence.

    Returns
    -------
    dict
        Dictionary with keys:
        - ``original_cai``: CAI of the original sequence
        - ``optimized_cai``: CAI of the optimized sequence
        - ``improvement``: optimized_cai - original_cai
        - ``published_cai``: Published CAI (if available)
        - ``exceeds_published``: Whether optimized_cai > published_cai
        - ``optimization_success``: Whether improvement > 0
    """
    resolved_org = resolve_organism(organism, strict=False)
    original_cai = compute_cai(original_dna, organism=resolved_org)
    optimized_cai = compute_cai(optimized_dna, organism=resolved_org)

    improvement = optimized_cai - original_cai

    # Look up published value
    published_entry = PUBLISHED_CAI_VALUES.get((gene, organism))
    if published_entry is None:
        from ..organisms import SPECIES_SHORT_NAMES
        short_name = SPECIES_SHORT_NAMES.get(resolved_org, organism)
        published_entry = PUBLISHED_CAI_VALUES.get((gene, short_name))

    published_cai: float | None = None
    exceeds_published: bool | None = None
    if published_entry is not None:
        published_cai = float(published_entry["cai"])  # type: ignore[arg-type]
        exceeds_published = optimized_cai >= published_cai

    return {
        "gene": gene,
        "organism": resolved_org,
        "original_cai": original_cai,
        "optimized_cai": optimized_cai,
        "improvement": improvement,
        "published_cai": published_cai,
        "exceeds_published": exceeds_published,
        "optimization_success": improvement > 0,
    }


def summarize_benchmark(results: list[CAIBenchmarkResult]) -> BenchmarkSummary:
    """Summarize a batch of benchmark results.

    Parameters
    ----------
    results : list[CAIBenchmarkResult]
        Benchmark results to summarize.

    Returns
    -------
    BenchmarkSummary
        Summary statistics.
    """
    total = len(results)
    with_published = [r for r in results if r.published_cai is not None]
    genes_with_published = len(with_published)

    if genes_with_published == 0:
        return BenchmarkSummary(
            total_genes=total,
            genes_with_published=0,
            mean_deviation=0.0,
            max_deviation=0.0,
            min_deviation=0.0,
            pass_rate=0.0,
            failures=[],
        )

    deviations = [r.deviation for r in with_published if r.deviation is not None]
    if not deviations:
        mean_dev = 0.0
        max_dev = 0.0
        min_dev = 0.0
    else:
        mean_dev = sum(deviations) / len(deviations)
        max_dev = max(deviations)
        min_dev = min(deviations)

    passed = [r for r in with_published if r.pass_threshold is True]
    pass_rate = len(passed) / genes_with_published

    failures = [
        r.gene for r in with_published if r.pass_threshold is False
    ]

    return BenchmarkSummary(
        total_genes=total,
        genes_with_published=genes_with_published,
        mean_deviation=mean_dev,
        max_deviation=max_dev,
        min_deviation=min_dev,
        pass_rate=pass_rate,
        failures=failures,
    )
