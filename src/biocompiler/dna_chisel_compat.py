"""
BioCompiler DNA Chisel Compatibility Layer — Comparative Benchmarking

Optional integration with DNA Chisel (https://github.com/Edinburgh-Genome-Foundry/DnaChisel)
for comparative benchmarking and alternative optimization. This module provides:

- ``optimize_with_dna_chisel`` — Use DNA Chisel's constraint solver as an alternative
- ``compare_optimizers`` — Run both BioCompiler and DNA Chisel side-by-side
- ``run_comparative_benchmark`` — Full benchmark suite comparing both tools

DNA Chisel is an OPTIONAL dependency. All public functions degrade gracefully
when it is not installed, returning informative error dicts instead of raising.

Design Philosophy:
    This module demonstrates that BioCompiler can work *alongside* DNA Chisel,
    not just compete with it. Users may prefer DNA Chisel's declarative
    constraint language for some tasks, and BioCompiler's type-system/certificate
    approach for others. Comparative benchmarking helps users choose.

DNA Chisel API Mapping:
    BioCompiler constraints -> DNA Chisel specifications:
        - GC range (gc_lo, gc_hi) -> EnforceGCContent
        - Restriction site avoidance -> AvoidPattern
        - Translation fidelity -> EnforceTranslation
        - CAI threshold -> (no direct equivalent; post-hoc check)

    DNA Chisel uses a different amino acid representation. We convert
    BioCompiler's single-letter codes to DNA Chisel's expected format.
"""

import logging
import time
from dataclasses import dataclass, field
from typing import Optional

from .benchmark import REFERENCE_GENES
from .optimization import optimize_sequence, OptimizationResult
from .translation import compute_cai
from .scanner import gc_content
from .constants import RESTRICTION_ENZYMES, CODON_TABLE, AA_TO_CODONS

logger = logging.getLogger(__name__)

# ─── DNA Chisel availability check ───────────────────────────────────

_DNA_CHISEL_AVAILABLE: bool = False
_DNA_CHISEL_ERROR: str = ""

try:
    from dnachisel import (
        DnaOptimizationProblem,
        AvoidPattern,
        EnforceGCContent,
        EnforceTranslation,
    )
    _DNA_CHISEL_AVAILABLE = True
    logger.debug("DNA Chisel is available for comparative benchmarking")
except ImportError as exc:
    _DNA_CHISEL_ERROR = str(exc)
    logger.debug(
        "DNA Chisel not installed — comparative benchmarking will be limited. "
        "Install with: pip install 'biocompiler[compare]'"
    )


def is_dna_chisel_available() -> bool:
    """Return True if DNA Chisel is installed and importable."""
    return _DNA_CHISEL_AVAILABLE


# ─── Amino acid conversion utilities ─────────────────────────────────

# DNA Chisel uses standard single-letter amino acid codes, same as BioCompiler.
# However, DNA Chisel's EnforceTranslation expects the amino acid sequence
# to be provided in a specific format. We need to ensure compatibility.

# Standard 1-letter to 3-letter amino acid mapping (for DNA Chisel interop)
_AA_ONE_TO_THREE: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}

# BioCompiler organism name -> DNA Chisel organism name mapping
_ORGANISM_MAP: dict[str, str] = {
    "Homo_sapiens": "h_sapiens",
    "Escherichia_coli": "e_coli",
    "Mus_musculus": "m_musculus",
    "CHO_K1": "c_griseus",
    "Saccharomyces_cerevisiae": "s_cerevisiae",
}


def _build_initial_sequence(protein: str, organism: str = "Homo_sapiens") -> str:
    """
    Build an initial DNA sequence from a protein using highest-CAI codons.

    DNA Chisel requires an initial sequence to optimize from. We seed it
    with BioCompiler's best-codon-per-position sequence, giving DNA Chisel
    a strong starting point.

    Args:
        protein: Amino acid sequence (single-letter codes)
        organism: Target organism for codon usage

    Returns:
        DNA sequence using preferred codons for each amino acid
    """
    from .organisms import CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS

    if organism not in SUPPORTED_ORGANISMS:
        # Fall back to human if organism not supported
        organism = "Homo_sapiens"

    usage = CODON_ADAPTIVENESS_TABLES.get(organism, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])

    # Sort codons by adaptiveness for each amino acid
    sorted_codons: dict[str, list[str]] = {}
    for aa in set(protein):
        codons = AA_TO_CODONS.get(aa, [])
        codons_sorted = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        sorted_codons[aa] = codons_sorted

    # Build sequence using best codon per position
    sequence_chars: list[str] = []
    for aa in protein:
        codons = sorted_codons.get(aa, [])
        if codons:
            sequence_chars.append(codons[0])
        else:
            # Unknown amino acid — use a neutral codon
            logger.warning("Unknown amino acid '%s', using NNN placeholder", aa)
            sequence_chars.append("NNN")

    return "".join(sequence_chars)


def _build_dna_chisel_spec(
    protein: str,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_enzymes: list[str] | None = None,
) -> list:
    """
    Convert BioCompiler constraints to DNA Chisel specification format.

    Creates a list of DNA Chisel constraint objects that correspond to
    BioCompiler's type predicates:

    - EnforceTranslation: ensure the sequence encodes the target protein
    - EnforceGCContent: keep GC content within bounds
    - AvoidPattern: avoid restriction enzyme recognition sites

    Note: DNA Chisel does not have a direct CAI optimization constraint.
    CAI is handled indirectly through codon usage preferences in the
    initial sequence seeding, and checked post-hoc.

    Args:
        protein: Target protein sequence (single-letter codes)
        gc_lo: Minimum GC content fraction
        gc_hi: Maximum GC content fraction
        restriction_enzymes: List of enzyme names to avoid

    Returns:
        List of DNA Chisel Specification objects

    Raises:
        ImportError: If DNA Chisel is not installed
    """
    if not _DNA_CHISEL_AVAILABLE:
        raise ImportError(
            f"DNA Chisel is required for specification building: {_DNA_CHISEL_ERROR}. "
            "Install with: pip install 'biocompiler[compare]'"
        )

    constraints = []

    # Enforce that the sequence translates to our target protein
    # DNA Chisel's EnforceTranslation takes 'translation' keyword
    constraints.append(EnforceTranslation(translation=protein))

    # Enforce GC content range
    # DNA Chisel uses 'mini'/'maxi' (not 'gc_min'/'gc_max')
    constraints.append(
        EnforceGCContent(
            mini=gc_lo,
            maxi=gc_hi,
            window=50,  # Use a 50bp sliding window for local GC enforcement
        )
    )

    # Avoid restriction enzyme sites
    if restriction_enzymes:
        for enz_name in restriction_enzymes:
            site = RESTRICTION_ENZYMES.get(enz_name)
            if site:
                # Avoid the pattern on both strands
                constraints.append(AvoidPattern(site))
                # DNA Chisel's AvoidPattern handles reverse complement automatically

    return constraints


def _count_restriction_sites(sequence: str, restriction_enzymes: list[str] | None = None) -> int:
    """
    Count restriction enzyme recognition sites in a sequence.

    Checks both forward and reverse complement strands.

    Args:
        sequence: DNA sequence to check
        restriction_enzymes: List of enzyme names to check for

    Returns:
        Total count of restriction sites found
    """
    from .constants import reverse_complement

    if not restriction_enzymes:
        restriction_enzymes = list(RESTRICTION_ENZYMES.keys())

    count = 0
    seq_upper = sequence.upper()
    for enz_name in restriction_enzymes:
        site = RESTRICTION_ENZYMES.get(enz_name, "")
        if not site:
            continue
        site_upper = site.upper()
        # Only check concrete sites (no IUPAC ambiguity)
        if any(b not in "ACGT" for b in site_upper):
            continue
        # Forward strand
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
        # Reverse complement strand
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:  # Avoid double-counting palindromes
            start = 0
            while True:
                pos = seq_upper.find(site_rc, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1

    return count


# ─── Public API ───────────────────────────────────────────────────────

@dataclass
class ChiselResult:
    """Result of DNA Chisel optimization."""
    sequence: str
    protein: str
    cai: float
    gc_content: float
    restriction_site_count: int
    execution_time_s: float
    success: bool
    error: str | None = None


@dataclass
class ComparisonResult:
    """Side-by-side comparison of BioCompiler vs DNA Chisel."""
    protein: str
    organism: str
    biocompiler: dict
    dna_chisel: dict | None
    dna_chisel_available: bool
    winner: dict  # Per-metric comparison


def optimize_with_dna_chisel(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_enzymes: list[str] | None = None,
) -> ChiselResult:
    """
    Optimize a protein sequence using DNA Chisel's constraint solver.

    This is an alternative to BioCompiler's own optimizer. DNA Chisel uses
    a different algorithm (random mutation + constraint propagation) which
    may produce different results, especially for long sequences.

    Args:
        protein: Target protein sequence (single-letter amino acid codes)
        organism: Target organism for codon usage
        gc_lo: Minimum GC content fraction
        gc_hi: Maximum GC content fraction
        restriction_enzymes: List of enzyme names whose sites should be avoided

    Returns:
        ChiselResult with optimized sequence and metrics

    Note:
        Requires DNA Chisel to be installed. Returns a ChiselResult with
        success=False if DNA Chisel is not available.
    """
    if not _DNA_CHISEL_AVAILABLE:
        return ChiselResult(
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=0.0,
            success=False,
            error=f"DNA Chisel not installed: {_DNA_CHISEL_ERROR}. "
                  "Install with: pip install 'biocompiler[compare]'",
        )

    enzyme_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]  # Limit for performance

    t0 = time.perf_counter()
    try:
        # Build initial sequence from preferred codons
        initial_seq = _build_initial_sequence(protein, organism)

        # Build constraint specification
        constraints = _build_dna_chisel_spec(
            protein, gc_lo, gc_hi, enzyme_names
        )

        # Create and solve the optimization problem
        problem = DnaOptimizationProblem(
            sequence=initial_seq,
            constraints=constraints,
        )

        # Resolve constraints
        problem.resolve_constraints()

        # Extract the optimized sequence
        optimized = str(problem.sequence)

        elapsed = time.perf_counter() - t0

        # Compute metrics using BioCompiler's own evaluators for fair comparison
        cai = compute_cai(optimized, organism)
        gc = gc_content(optimized)
        rs_count = _count_restriction_sites(optimized, enzyme_names)

        return ChiselResult(
            sequence=optimized,
            protein=protein,
            cai=cai,
            gc_content=gc,
            restriction_site_count=rs_count,
            execution_time_s=elapsed,
            success=True,
        )

    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("DNA Chisel optimization failed: %s", exc, exc_info=True)
        return ChiselResult(
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=elapsed,
            success=False,
            error=str(exc),
        )


def compare_optimizers(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_enzymes: list[str] | None = None,
    cai_threshold: float = 0.2,
) -> ComparisonResult:
    """
    Run both BioCompiler and DNA Chisel optimization on the same protein,
    returning a structured comparison.

    This is the primary entry point for users who want to evaluate which
    optimizer produces better results for their specific protein. All
    metrics are computed using BioCompiler's evaluators for fairness.

    Args:
        protein: Target protein sequence (single-letter amino acid codes)
        organism: Target organism for codon usage
        gc_lo: Minimum GC content fraction
        gc_hi: Maximum GC content fraction
        restriction_enzymes: List of enzyme names whose sites should be avoided
        cai_threshold: Minimum CAI score for BioCompiler's optimizer

    Returns:
        ComparisonResult with metrics for both optimizers and per-metric winner
    """
    enzyme_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]

    # ── Run BioCompiler optimization ──
    t0_bc = time.perf_counter()
    try:
        bc_result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
        )
        bc_elapsed = time.perf_counter() - t0_bc
        bc_metrics = {
            "sequence": bc_result.sequence,
            "sequence_length": len(bc_result.sequence),
            "cai": bc_result.cai,
            "gc_content": bc_result.gc_content,
            "restriction_site_count": _count_restriction_sites(
                bc_result.sequence, enzyme_names
            ),
            "execution_time_s": round(bc_elapsed, 4),
            "success": True,
            "satisfied_predicates": bc_result.satisfied_predicates,
            "failed_predicates": bc_result.failed_predicates,
            "fallback_used": bc_result.fallback_used,
        }
    except Exception as exc:
        bc_elapsed = time.perf_counter() - t0_bc
        logger.error("BioCompiler optimization failed: %s", exc, exc_info=True)
        bc_metrics = {
            "sequence": "",
            "sequence_length": 0,
            "cai": 0.0,
            "gc_content": 0.0,
            "restriction_site_count": 0,
            "execution_time_s": round(bc_elapsed, 4),
            "success": False,
            "error": str(exc),
        }

    # ── Run DNA Chisel optimization ──
    chisel_result = optimize_with_dna_chisel(
        protein=protein,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        restriction_enzymes=enzyme_names,
    )

    dc_metrics: dict | None = None
    if chisel_result.success:
        dc_metrics = {
            "sequence": chisel_result.sequence,
            "sequence_length": len(chisel_result.sequence),
            "cai": chisel_result.cai,
            "gc_content": chisel_result.gc_content,
            "restriction_site_count": chisel_result.restriction_site_count,
            "execution_time_s": round(chisel_result.execution_time_s, 4),
            "success": True,
        }
    elif not _DNA_CHISEL_AVAILABLE:
        dc_metrics = None  # Not available at all
    else:
        dc_metrics = {
            "sequence": "",
            "sequence_length": 0,
            "cai": 0.0,
            "gc_content": 0.0,
            "restriction_site_count": 0,
            "execution_time_s": round(chisel_result.execution_time_s, 4),
            "success": False,
            "error": chisel_result.error,
        }

    # ── Determine winners ──
    winner = _compute_winners(bc_metrics, dc_metrics, gc_lo, gc_hi)

    return ComparisonResult(
        protein=protein,
        organism=organism,
        biocompiler=bc_metrics,
        dna_chisel=dc_metrics,
        dna_chisel_available=_DNA_CHISEL_AVAILABLE,
        winner=winner,
    )


def _compute_winners(
    bc: dict, dc: dict | None, gc_lo: float, gc_hi: float
) -> dict:
    """
    Compare BioCompiler and DNA Chisel metrics and determine per-metric winners.

    Higher-is-better: CAI
    Lower-is-better: restriction_site_count, execution_time_s
    Target-range: gc_content (closer to midpoint is better)

    Args:
        bc: BioCompiler metrics dict
        dc: DNA Chisel metrics dict (None if unavailable)
        gc_lo: Minimum GC content
        gc_hi: Maximum GC content

    Returns:
        Dict with per-metric winner and overall assessment
    """
    result: dict = {
        "metrics": {},
        "overall": "biocompiler",  # Default if no comparison possible
    }

    if dc is None or not dc.get("success"):
        # No valid DNA Chisel result — BioCompiler wins by default
        result["overall"] = "biocompiler (dna_chisel_unavailable)"
        for metric in ("cai", "gc_content", "restriction_site_count", "execution_time_s"):
            result["metrics"][metric] = {
                "biocompiler": bc.get(metric),
                "dna_chisel": None,
                "winner": "biocompiler",
            }
        return result

    if not bc.get("success"):
        result["overall"] = "dna_chisel (biocompiler_failed)"
        for metric in ("cai", "gc_content", "restriction_site_count", "execution_time_s"):
            result["metrics"][metric] = {
                "biocompiler": None,
                "dna_chisel": dc.get(metric),
                "winner": "dna_chisel",
            }
        return result

    # Both succeeded — fair comparison
    bc_wins = 0
    dc_wins = 0

    # CAI: higher is better
    bc_cai = bc.get("cai", 0.0)
    dc_cai = dc.get("cai", 0.0)
    cai_winner = "tie"
    if bc_cai > dc_cai + 0.001:
        cai_winner = "biocompiler"
        bc_wins += 1
    elif dc_cai > bc_cai + 0.001:
        cai_winner = "dna_chisel"
        dc_wins += 1
    result["metrics"]["cai"] = {
        "biocompiler": bc_cai,
        "dna_chisel": dc_cai,
        "winner": cai_winner,
    }

    # GC content: closer to midpoint is better
    gc_mid = (gc_lo + gc_hi) / 2.0
    bc_gc_diff = abs(bc.get("gc_content", 0.0) - gc_mid)
    dc_gc_diff = abs(dc.get("gc_content", 0.0) - gc_mid)
    gc_winner = "tie"
    if bc_gc_diff < dc_gc_diff - 0.001:
        gc_winner = "biocompiler"
        bc_wins += 1
    elif dc_gc_diff < bc_gc_diff - 0.001:
        gc_winner = "dna_chisel"
        dc_wins += 1
    result["metrics"]["gc_content"] = {
        "biocompiler": bc.get("gc_content"),
        "dna_chisel": dc.get("gc_content"),
        "winner": gc_winner,
        "note": f"target_midpoint={gc_mid:.3f}",
    }

    # Restriction sites: fewer is better
    bc_rs = bc.get("restriction_site_count", 999)
    dc_rs = dc.get("restriction_site_count", 999)
    rs_winner = "tie"
    if bc_rs < dc_rs:
        rs_winner = "biocompiler"
        bc_wins += 1
    elif dc_rs < bc_rs:
        rs_winner = "dna_chisel"
        dc_wins += 1
    result["metrics"]["restriction_site_count"] = {
        "biocompiler": bc_rs,
        "dna_chisel": dc_rs,
        "winner": rs_winner,
    }

    # Execution time: faster is better
    bc_time = bc.get("execution_time_s", float("inf"))
    dc_time = dc.get("execution_time_s", float("inf"))
    time_winner = "tie"
    if bc_time < dc_time * 0.9:  # 10% margin to avoid counting noise
        time_winner = "biocompiler"
        bc_wins += 1
    elif dc_time < bc_time * 0.9:
        time_winner = "dna_chisel"
        dc_wins += 1
    result["metrics"]["execution_time_s"] = {
        "biocompiler": bc_time,
        "dna_chisel": dc_time,
        "winner": time_winner,
    }

    # Overall winner
    if bc_wins > dc_wins:
        result["overall"] = "biocompiler"
    elif dc_wins > bc_wins:
        result["overall"] = "dna_chisel"
    else:
        result["overall"] = "tie"

    return result


# ─── Comparative Benchmark Suite ──────────────────────────────────────

@dataclass
class ComparativeBenchmarkReport:
    """Complete comparative benchmark report."""
    timestamp: str
    dna_chisel_available: bool
    gene_results: list[dict] = field(default_factory=list)
    summary: dict = field(default_factory=dict)

    @property
    def total_genes(self) -> int:
        return len(self.gene_results)


def run_comparative_benchmark(
    genes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
) -> ComparativeBenchmarkReport:
    """
    Benchmark BioCompiler vs DNA Chisel on the reference gene set.

    For each reference gene, extracts the protein, runs both optimizers,
    and records metrics for comparison. The report shows which tool
    performs better on each metric across all genes.

    This function works correctly even without DNA Chisel installed —
    it will simply report BioCompiler's results with a note that DNA
    Chisel was unavailable.

    Args:
        genes: Subset of gene names to benchmark (None = all REFERENCE_GENES)
        gc_lo: Minimum GC content for optimization
        gc_hi: Maximum GC content for optimization
        cai_threshold: CAI threshold for BioCompiler's optimizer

    Returns:
        ComparativeBenchmarkReport with per-gene results and summary
    """
    from datetime import datetime, timezone

    gene_names = genes or list(REFERENCE_GENES.keys())
    gene_results: list[dict] = []

    for gene_name in gene_names:
        gene_data = REFERENCE_GENES.get(gene_name)
        if not gene_data:
            logger.warning("Unknown gene: %s, skipping", gene_name)
            continue

        logger.info("Comparative benchmark: %s", gene_name)

        # Extract protein from reference gene
        seq = gene_data["pre_mrna"].replace(" ", "")
        exons = gene_data["exon_boundaries"]
        organism = gene_data["organism"]

        # Get the coding sequence and translate
        coding_seq = "".join(seq[start:end] for start, end in exons)
        from .translation import translate
        protein = translate(coding_seq).rstrip("*")

        if not protein:
            gene_results.append({
                "gene": gene_name,
                "error": "Empty protein translation",
                "biocompiler": None,
                "dna_chisel": None,
            })
            continue

        # Run comparison
        comparison = compare_optimizers(
            protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
        )

        gene_results.append({
            "gene": gene_name,
            "description": gene_data["description"],
            "organism": organism,
            "protein_length": len(protein),
            "biocompiler": comparison.biocompiler,
            "dna_chisel": comparison.dna_chisel,
            "winner": comparison.winner,
        })

    # Compute summary
    summary = _compute_comparative_summary(gene_results)

    return ComparativeBenchmarkReport(
        timestamp=datetime.now(timezone.utc).isoformat(),
        dna_chisel_available=_DNA_CHISEL_AVAILABLE,
        gene_results=gene_results,
        summary=summary,
    )


def _compute_comparative_summary(gene_results: list[dict]) -> dict:
    """
    Compute aggregate summary statistics from per-gene benchmark results.

    Determines how many genes each optimizer wins on each metric,
    and produces an overall assessment.

    Args:
        gene_results: List of per-gene comparison dicts

    Returns:
        Summary dict with per-metric win counts and overall assessment
    """
    summary: dict = {
        "total_genes": len(gene_results),
        "genes_with_errors": 0,
        "metric_wins": {
            "cai": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
            "gc_content": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
            "restriction_site_count": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
            "execution_time_s": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
        },
        "overall_wins": {"biocompiler": 0, "dna_chisel": 0, "tie": 0, "unavailable": 0},
        "avg_cai": {"biocompiler": 0.0, "dna_chisel": 0.0},
        "avg_gc": {"biocompiler": 0.0, "dna_chisel": 0.0},
        "avg_restriction_sites": {"biocompiler": 0.0, "dna_chisel": 0.0},
        "avg_execution_time_s": {"biocompiler": 0.0, "dna_chisel": 0.0},
    }

    bc_cais: list[float] = []
    dc_cais: list[float] = []
    bc_gcs: list[float] = []
    dc_gcs: list[float] = []
    bc_rss: list[int] = []
    dc_rss: list[int] = []
    bc_times: list[float] = []
    dc_times: list[float] = []

    for result in gene_results:
        if result.get("error"):
            summary["genes_with_errors"] += 1
            continue

        winner_info = result.get("winner", {})
        bc_data = result.get("biocompiler", {})
        dc_data = result.get("dna_chisel")

        # Collect averages
        if bc_data and bc_data.get("success"):
            bc_cais.append(bc_data.get("cai", 0.0))
            bc_gcs.append(bc_data.get("gc_content", 0.0))
            bc_rss.append(bc_data.get("restriction_site_count", 0))
            bc_times.append(bc_data.get("execution_time_s", 0.0))

        if dc_data and dc_data.get("success"):
            dc_cais.append(dc_data.get("cai", 0.0))
            dc_gcs.append(dc_data.get("gc_content", 0.0))
            dc_rss.append(dc_data.get("restriction_site_count", 0))
            dc_times.append(dc_data.get("execution_time_s", 0.0))

        # Count metric wins
        metrics = winner_info.get("metrics", {})
        for metric_name in ("cai", "gc_content", "restriction_site_count", "execution_time_s"):
            metric_winner = metrics.get(metric_name, {}).get("winner", "unavailable")
            if metric_winner in summary["metric_wins"][metric_name]:
                summary["metric_wins"][metric_name][metric_winner] += 1
            else:
                summary["metric_wins"][metric_name]["unavailable"] += 1

        # Count overall wins
        overall = winner_info.get("overall", "unavailable")
        if overall in summary["overall_wins"]:
            summary["overall_wins"][overall] += 1
        else:
            summary["overall_wins"]["unavailable"] += 1

    # Compute averages
    if bc_cais:
        summary["avg_cai"]["biocompiler"] = round(sum(bc_cais) / len(bc_cais), 4)
    if dc_cais:
        summary["avg_cai"]["dna_chisel"] = round(sum(dc_cais) / len(dc_cais), 4)
    if bc_gcs:
        summary["avg_gc"]["biocompiler"] = round(sum(bc_gcs) / len(bc_gcs), 4)
    if dc_gcs:
        summary["avg_gc"]["dna_chisel"] = round(sum(dc_gcs) / len(dc_gcs), 4)
    if bc_rss:
        summary["avg_restriction_sites"]["biocompiler"] = round(sum(bc_rss) / len(bc_rss), 2)
    if dc_rss:
        summary["avg_restriction_sites"]["dna_chisel"] = round(sum(dc_rss) / len(dc_rss), 2)
    if bc_times:
        summary["avg_execution_time_s"]["biocompiler"] = round(sum(bc_times) / len(bc_times), 4)
    if dc_times:
        summary["avg_execution_time_s"]["dna_chisel"] = round(sum(dc_times) / len(dc_times), 4)

    return summary


def format_comparative_report_text(report: ComparativeBenchmarkReport) -> str:
    """
    Format a comparative benchmark report as human-readable text.

    Args:
        report: ComparativeBenchmarkReport to format

    Returns:
        Multi-line string with formatted report
    """
    lines = [
        "BioCompiler vs DNA Chisel — Comparative Benchmark Report",
        f"Timestamp: {report.timestamp}",
        f"DNA Chisel Available: {'Yes' if report.dna_chisel_available else 'No'}",
        "",
    ]

    if not report.dna_chisel_available:
        lines.append(
            "NOTE: DNA Chisel is not installed. Only BioCompiler results are shown."
        )
        lines.append(
            "Install with: pip install 'biocompiler[compare]'"
        )
        lines.append("")

    for gr in report.gene_results:
        gene = gr.get("gene", "unknown")
        desc = gr.get("description", "")
        lines.append(f"  Gene: {gene} — {desc}")

        bc = gr.get("biocompiler")
        if bc and bc.get("success"):
            lines.append(
                f"    BioCompiler: CAI={bc['cai']:.4f}, GC={bc['gc_content']:.3f}, "
                f"RS={bc['restriction_site_count']}, "
                f"Time={bc['execution_time_s']:.3f}s"
            )
        else:
            lines.append(f"    BioCompiler: FAILED — {bc.get('error', 'unknown')}")

        dc = gr.get("dna_chisel")
        if dc is None:
            lines.append("    DNA Chisel:   NOT AVAILABLE")
        elif dc.get("success"):
            lines.append(
                f"    DNA Chisel:   CAI={dc['cai']:.4f}, GC={dc['gc_content']:.3f}, "
                f"RS={dc['restriction_site_count']}, "
                f"Time={dc['execution_time_s']:.3f}s"
            )
        else:
            lines.append(f"    DNA Chisel:   FAILED — {dc.get('error', 'unknown')}")

        winner = gr.get("winner", {})
        overall = winner.get("overall", "N/A")
        lines.append(f"    Winner: {overall}")
        lines.append("")

    # Summary
    summary = report.summary
    lines.append("Summary:")
    lines.append(f"  Total genes: {summary.get('total_genes', 0)}")
    lines.append(f"  Genes with errors: {summary.get('genes_with_errors', 0)}")

    avg = summary.get("avg_cai", {})
    lines.append(
        f"  Average CAI: BioCompiler={avg.get('biocompiler', 0.0):.4f}, "
        f"DNA Chisel={avg.get('dna_chisel', 0.0):.4f}"
    )
    avg_gc = summary.get("avg_gc", {})
    lines.append(
        f"  Average GC:  BioCompiler={avg_gc.get('biocompiler', 0.0):.3f}, "
        f"DNA Chisel={avg_gc.get('dna_chisel', 0.0):.3f}"
    )
    avg_rs = summary.get("avg_restriction_sites", {})
    lines.append(
        f"  Average RS:  BioCompiler={avg_rs.get('biocompiler', 0.0):.1f}, "
        f"DNA Chisel={avg_rs.get('dna_chisel', 0.0):.1f}"
    )
    avg_time = summary.get("avg_execution_time_s", {})
    lines.append(
        f"  Average Time: BioCompiler={avg_time.get('biocompiler', 0.0):.3f}s, "
        f"DNA Chisel={avg_time.get('dna_chisel', 0.0):.3f}s"
    )

    # Per-metric wins
    lines.append("")
    lines.append("Per-Metric Wins:")
    metric_wins = summary.get("metric_wins", {})
    for metric, wins in metric_wins.items():
        bc_w = wins.get("biocompiler", 0)
        dc_w = wins.get("dna_chisel", 0)
        ties = wins.get("tie", 0)
        lines.append(f"  {metric}: BC={bc_w}, DC={dc_w}, Tie={ties}")

    overall_wins = summary.get("overall_wins", {})
    lines.append("")
    lines.append("Overall Wins:")
    for opt, count in overall_wins.items():
        lines.append(f"  {opt}: {count}")

    return "\n".join(lines)
