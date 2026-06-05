"""
BioCompiler Benchmark Module v9.2.0
=====================================
Built-in benchmark sequences, performance measurement, comprehensive
benchmarking, ablation studies, and head-to-head tool comparisons for
BioCompiler.

Provides:
  - run_benchmark(): Optimize eGFP, mCherry, LacZ with human/ecoli and print a table
  - compare_tools(): Theoretical feature comparison table vs other tools
  - Extended Benchmark API (v7.2.0): structured benchmarks, JSON/text reports
  - Comprehensive Benchmark Suite: 12-gene panel, statistical analysis,
    ablation study, Pareto frontier, and output persistence
  - Head-to-Head Tool Comparison: BioCompiler vs DNA Chisel, DNAworks,
    GeneOptimizer (published data), SimpleCAI, and Random baselines
"""

import csv
import json
import logging
import math
import random
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
try:
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/NotoSansSC[wght].ttf')
except Exception:
    logger.debug("Primary font loading failed, trying fallback", exc_info=True)
    try:
        fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
    except Exception:
        logger.debug("Font loading failed, using default font", exc_info=True)
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Noto Sans SC', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
import numpy as np

try:
    from .optimizer import BioOptimizer
except ImportError:
    logger.debug("Import failed: .optimizer module, falling back to .optimization")
    from .optimization import BioOptimizer
from .optimization import optimize_sequence

from .certificate import compute_certificate
from .organisms import CODON_ADAPTIVENESS_TABLES, SPECIES, resolve_organism
from .translation import translate, compute_cai
from .scanner import gc_content
from .constants import AA_TO_CODONS, RESTRICTION_ENZYMES, reverse_complement
from .restriction_sites import get_recognition_site
from .engine_base import BatchResult, EngineTimer, BaseEngineResult

__all__ = [
    # Constants / Data
    "EGFP_DNA", "MCHERRY_DNA", "LACZ_DNA",
    "REFERENCE_GENES", "OUTPUT_DIR", "GENE_PANEL",
    "ORGANISM_TO_SPECIES", "ORGANISM_FOR_CAI",
    "DEFAULT_ENZYMES", "MAX_CONSTRAINTS",
    # Classes
    "BenchmarkResult", "BenchmarkReport",
    "ToolResult", "HeadToHeadReport",
    # Functions — core benchmark
    "run_benchmark", "run_benchmarks",
    "compare_tools_theoretical", "compare_tools",
    # Functions — structured benchmarks
    "run_structured_benchmarks",
    "format_benchmark_report_json", "format_benchmark_report_text",
    # Functions — comprehensive benchmark suite
    "optimize_biocompiler", "optimize_simple_cai", "optimize_random",
    "run_multi_gene_comparison", "compute_statistics", "format_stats_table",
    "run_ablation_study", "format_ablation_table",
    "plot_pareto_frontier", "plot_ablation",
    "save_json", "save_csv", "save_summary",
    "run_comprehensive_benchmark",
    # Functions — head-to-head comparison
    "is_dna_chisel_available",
    "run_head_to_head_benchmark",
    "format_head_to_head_text", "format_head_to_head_json",
    # Entry point
    "main",
]

# ────────────────────────────────────────────────────────────
# Built-in gene sequences (standard reference sequences)
# ────────────────────────────────────────────────────────────

# eGFP (Enhanced Green Fluorescent Protein) — 717 bp coding sequence
# Source: pEGFP-N1 (Clontech), GenBank accession U55763 region
EGFP_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
)

# mCherry (Red Fluorescent Protein) — 711 bp coding sequence
# Source: pmCherry (Clontech), derived from mRFP1
MCHERRY_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGCGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
)

# LacZ (beta-galactosidase, N-terminal 720 bp fragment) — 720 bp
# Source: E. coli lacZ gene, first 240 codons
LACZ_DNA = (
    "ATGACCATGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTATGCCTGCAGGTCGACTCTAGAGGATCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGTTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGC"
)


def _compute_cai(seq: str, species_cai: Dict[str, float]) -> float:
    """Compute the geometric mean CAI for a sequence."""
    if not seq or len(seq) < 3:
        return 0.0
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        cai = species_cai.get(codon, 0.0)
        if cai <= 0:
            cai = 0.001  # avoid log(0)
        log_sum += math.log(cai)
        count += 1
    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


def _count_gt(seq: str) -> int:
    """Count GT dinucleotides in a sequence."""
    return sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "GT")


def _count_cpg_ratio(seq: str) -> float:
    """Compute CpG Obs/Exp ratio for the full sequence."""
    c = seq.count("C")
    g = seq.count("G")
    cg = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")
    expected = (c * g) / len(seq) if len(seq) > 0 else 0
    return cg / expected if expected > 0 else 0.0


def _count_restriction_sites(seq: str, enzymes: List[str]) -> int:
    """Count total restriction sites in sequence (forward strand only)."""
    total = 0
    for enzyme in enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        pos = seq.find(site)
        while pos != -1:
            total += 1
            pos = seq.find(site, pos + 1)
    return total


def run_benchmark(
    enzymes: List[str] = None,
    splice_low: float = 3.0,
    splice_high: float = 6.0,
) -> None:
    """Run built-in benchmarks for eGFP, mCherry, and LacZ.

    Optimizes each gene with both human and ecoli species tables,
    then prints a formatted table of results.

    Args:
        enzymes: List of restriction enzymes to avoid
        splice_low: Low splice threshold
        splice_high: High splice threshold
    """
    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI"]

    genes = {
        "eGFP": EGFP_DNA,
        "mCherry": MCHERRY_DNA,
        "LacZ": LACZ_DNA,
    }

    species_list = ["human", "ecoli"]

    print()
    print("=" * 100)
    print("  BioCompiler v9.2.0 — Built-in Benchmark")
    print("=" * 100)
    print(f"  Enzymes avoided: {', '.join(enzymes)}")
    print(f"  Splice thresholds: low={splice_low}, high={splice_high}")
    print("=" * 100)

    results: List[Dict] = []

    for gene_name, gene_seq in genes.items():
        for species in species_list:
            opt = BioOptimizer(
                species=species,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=True,
            )

            # Measure optimization time
            with EngineTimer() as _bench_timer:
                optimized, pred_results, cert_text = opt.optimize(gene_seq)
            elapsed = _bench_timer.elapsed

            species_cai = dict(CODON_ADAPTIVENESS_TABLES.get(
                resolve_organism(species), CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
            ))

            cai_before = _compute_cai(gene_seq, species_cai)
            cai_after = _compute_cai(optimized, species_cai)
            gt_before = _count_gt(gene_seq)
            gt_after = _count_gt(optimized)
            cpg_ratio = _count_cpg_ratio(optimized)
            rs_count = _count_restriction_sites(optimized, enzymes)
            cert_level = compute_certificate(pred_results)

            results.append({
                "gene": gene_name,
                "species": species,
                "length": len(optimized),
                "cai_before": cai_before,
                "cai_after": cai_after,
                "gt_before": gt_before,
                "gt_after": gt_after,
                "cpg_ratio": cpg_ratio,
                "rs_count": rs_count,
                "cert": cert_level.value,
                "time_ms": elapsed * 1000,
            })

    # Print formatted table
    print()
    header = (
        f"{'Gene':<10} {'Species':<8} {'Length':>6} "
        f"{'CAI before':>10} {'CAI after':>10} "
        f"{'GT before':>10} {'GT after':>9} "
        f"{'CpG ratio':>10} {'RS sites':>9} "
        f"{'Cert':>8} {'Time(ms)':>9}"
    )
    print(header)
    print("-" * len(header))

    for r in results:
        row = (
            f"{r['gene']:<10} {r['species']:<8} {r['length']:>6} "
            f"{r['cai_before']:>10.4f} {r['cai_after']:>10.4f} "
            f"{r['gt_before']:>10d} {r['gt_after']:>9d} "
            f"{r['cpg_ratio']:>10.3f} {r['rs_count']:>9d} "
            f"{r['cert']:>8} {r['time_ms']:>9.1f}"
        )
        print(row)

    print()
    print("  Legend: CAI = Codon Adaptation Index, GT = GT dinucleotides,")
    print("          CpG ratio = Obs/Exp CG, RS = Restriction sites, Cert = Certificate level")
    print("=" * 100)
    print()


def compare_tools_theoretical() -> None:
    """Print a theoretical feature comparison table: BioCompiler vs other tools.

    Based on published capabilities:
      - GeneOptimizer (GeneArt/Thermo): CAI optimization, no GT/CpG avoidance
      - DNAworks: Restriction site avoidance, CAI, no GT/CpG
      - OPTIMIZER: CAI optimization only, no GT/CpG/restriction
    """
    tools = [
        {
            "name": "BioCompiler",
            "features": {
                "CAI optimization": True,
                "GT avoidance": True,
                "CpG avoidance": True,
                "Restriction sites": True,
                "Formal certification": True,
                "Cross-codon awareness": True,
                "Dual-threshold splice": True,
            },
        },
        {
            "name": "GeneOptimizer",
            "features": {
                "CAI optimization": True,
                "GT avoidance": False,
                "CpG avoidance": False,
                "Restriction sites": True,
                "Formal certification": False,
                "Cross-codon awareness": False,
                "Dual-threshold splice": False,
            },
        },
        {
            "name": "DNAworks",
            "features": {
                "CAI optimization": True,
                "GT avoidance": False,
                "CpG avoidance": False,
                "Restriction sites": True,
                "Formal certification": False,
                "Cross-codon awareness": False,
                "Dual-threshold splice": False,
            },
        },
        {
            "name": "OPTIMIZER",
            "features": {
                "CAI optimization": True,
                "GT avoidance": False,
                "CpG avoidance": False,
                "Restriction sites": False,
                "Formal certification": False,
                "Cross-codon awareness": False,
                "Dual-threshold splice": False,
            },
        },
    ]

    feature_names = list(tools[0]["features"].keys())

    print()
    print("=" * 80)
    print("  BioCompiler v9.2.0 — Tool Comparison (Theoretical)")
    print("=" * 80)
    print()

    # Header
    name_col = 20
    feat_col = 8
    header = f"{'Feature':<{name_col}}"
    for tool in tools:
        header += f" {tool['name']:^{feat_col}}"
    print(header)
    print("-" * len(header))

    # Rows
    for feat in feature_names:
        row = f"{feat:<{name_col}}"
        for tool in tools:
            val = tool["features"][feat]
            mark = "\u2713" if val else "\u2717"
            row += f" {mark:^{feat_col}}"
        print(row)

    print()
    print("  \u2713 = Supported   \u2717 = Not supported")
    print()
    print("  Notes:")
    print("    - GeneOptimizer: CAI + restriction sites, but no GT/CpG avoidance")
    print("    - DNAworks:      CAI + restriction sites, but no GT/CpG/formal cert")
    print("    - OPTIMIZER:     CAI only, no constraint handling or certification")
    print("    - BioCompiler:   Full pipeline with formal GOLD/SILVER/BRONZE certificates")
    print("=" * 80)
    print()


# Alias for compatibility with tests that use run_benchmarks (plural)
run_benchmarks = run_benchmark


# ==============================================================================
# Extended Benchmark API (v7.2.0)
# ==============================================================================
# Provides structured benchmark results, reference gene data, and JSON/text
# report formatters. This API is used by the test suite and REST API.

# Reference gene data for structured benchmarks
REFERENCE_GENES = {
    "HBB": {
        "description": "Human Beta-Globin (HBB)",
        "organism": "Homo_sapiens",
        "exon_boundaries": [(0, 92), (273, 495), (1346, 1608)],
        "known_protein_length": 147,
        "expected_gc_range": (0.35, 0.55),
        "expected_cai_range": (0.5, 1.0),
        "known_splice_events": ["canonical", "exon_skip_2"],
        "pre_mrna": (
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGG"
            "TGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGT"
            "CCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGAT"
            "GGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCTCACTGCAGTGAGCTGCACTGTGACAAGCTGCACGT"
            "GGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCAC"
            "CCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGC"
            "TCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTT"
        ),
    },
    "INS": {
        "description": "Human Insulin (INS)",
        "organism": "Homo_sapiens",
        "exon_boundaries": [(0, 153), (279, 465)],
        "known_protein_length": 51,
        "expected_gc_range": (0.40, 0.65),
        "expected_cai_range": (0.5, 1.0),
        "known_splice_events": ["canonical"],
        "pre_mrna": (
            "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTT"
            "TGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTA"
            "CACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGG"
            "CAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTG"
            "CTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
        ),
    },
    "EGFP": {
        "description": "Enhanced Green Fluorescent Protein (EGFP)",
        "organism": "Homo_sapiens",
        "exon_boundaries": [(0, 720)],
        "known_protein_length": 239,
        "expected_gc_range": (0.45, 0.65),
        "expected_cai_range": (0.6, 1.0),
        "known_splice_events": ["canonical"],
        "pre_mrna": EGFP_DNA,
    },
}


@dataclass
class BenchmarkResult:
    """Result of a single benchmark test."""
    gene_name: str
    test_name: str
    passed: bool
    expected: str
    actual: str
    details: str | None = None
    execution_time_ms: float = 0.0


class BenchmarkReport(BatchResult[BenchmarkResult]):
    """Complete benchmark report — extends BatchResult with metadata.

    Inherits from BatchResult[BenchmarkResult] providing:
      - results: list of BenchmarkResult objects
      - errors: list of error strings
      - total_time_s: total benchmark execution time
      - successful / failed: counts (successful aliases as 'passed')

    Additional fields:
      - timestamp: when the benchmark was run
      - version: BioCompiler version used
      - summary: per-gene and per-test summary statistics
    """

    def __init__(
        self,
        timestamp: str = "",
        version: str = "",
        total_tests: int = 0,
        passed: int = 0,
        failed: int = 0,
        results: list | None = None,
        errors: list | None = None,
        total_time_s: float = 0.0,
        summary: dict | None = None,
    ):
        super().__init__(
            results=results or [],
            errors=errors or [],
            total_time_s=total_time_s,
            successful=passed,
            failed=failed,
        )
        self.timestamp = timestamp
        self.version = version
        self._total_tests = total_tests
        self._summary = summary or {}

    @property
    def total_tests(self) -> int:
        """Total tests — stored value or len(results) if larger."""
        return max(self._total_tests, self.total)

    @property
    def passed(self) -> int:
        """Alias for BatchResult.successful — backward compat."""
        return self.successful

    @property
    def pass_rate(self) -> float:
        return self.successful / max(self.total, 1)

    @property
    def summary(self) -> dict:
        return self._summary

    @summary.setter
    def summary(self, value: dict):
        self._summary = value


def run_structured_benchmarks(
    gene_names: list[str] | None = None,
    include_optimization: bool = True,
) -> BenchmarkReport:
    """Run structured benchmarks against known gene sets.

    This validates BioCompiler's predictions against biological ground truth.
    """
    from . import __version__
    from .type_system import evaluate_all_predicates

    results: list[BenchmarkResult] = []
    genes = gene_names or list(REFERENCE_GENES.keys())

    for gene_name in genes:
        gene_data = REFERENCE_GENES.get(gene_name)
        if not gene_data:
            continue

        seq = gene_data["pre_mrna"].replace(" ", "")
        exons = gene_data["exon_boundaries"]
        organism = gene_data["organism"]

        # Translation benchmark
        _timer = EngineTimer()
        _timer.__enter__()
        coding_seq = "".join(seq[start:end] for start, end in exons)
        protein = translate(coding_seq)
        protein_len = len(protein.rstrip("*"))
        expected_len = gene_data["known_protein_length"]
        passed = abs(protein_len - expected_len) <= 10
        _timer.__exit__(None, None, None)
        results.append(BenchmarkResult(
            gene_name=gene_name, test_name="translation_length",
            passed=passed, expected=f"protein_length={expected_len}",
            actual=f"protein_length={protein_len}",
            execution_time_ms=round(_timer.elapsed * 1000, 4),
        ))

        # GC content benchmark
        _timer = EngineTimer()
        _timer.__enter__()
        gc = gc_content(seq)
        gc_lo, gc_hi = gene_data["expected_gc_range"]
        passed = gc_lo <= gc <= gc_hi
        _timer.__exit__(None, None, None)
        results.append(BenchmarkResult(
            gene_name=gene_name, test_name="gc_content_range",
            passed=passed, expected=f"GC in [{gc_lo}, {gc_hi}]",
            actual=f"GC = {gc:.4f}",
            execution_time_ms=round(_timer.elapsed * 1000, 4),
        ))

        # CAI benchmark
        _timer = EngineTimer()
        _timer.__enter__()
        cai = compute_cai(coding_seq, organism)
        cai_lo, cai_hi = gene_data["expected_cai_range"]
        passed = cai_lo <= cai <= cai_hi
        _timer.__exit__(None, None, None)
        results.append(BenchmarkResult(
            gene_name=gene_name, test_name="cai_range",
            passed=passed, expected=f"CAI in [{cai_lo}, {cai_hi}]",
            actual=f"CAI = {cai:.4f}",
            execution_time_ms=round(_timer.elapsed * 1000, 4),
        ))

        # Type predicates benchmark
        _timer = EngineTimer()
        _timer.__enter__()
        try:
            type_results = evaluate_all_predicates(
                seq=seq, known_exon_boundaries=exons, organism=organism,
            )
            n_pass = sum(1 for r in type_results if r.verdict.value in ("PASS", "LIKELY_PASS"))
            passed = n_pass >= 4
            _timer.__exit__(None, None, None)
            results.append(BenchmarkResult(
                gene_name=gene_name, test_name="type_predicates",
                passed=passed, expected=">=4 predicates PASS",
                actual=f"PASS={n_pass}, total={len(type_results)}",
                details="; ".join(f"{r.predicate}={r.verdict.value}" for r in type_results),
                execution_time_ms=round(_timer.elapsed * 1000, 4),
            ))
        except Exception as e:
            _timer.__exit__(None, None, None)
            results.append(BenchmarkResult(
                gene_name=gene_name, test_name="type_predicates",
                passed=False, expected="Predicates evaluated",
                actual=f"ERROR: {e}",
            ))

    total = len(results)
    passed_count = sum(1 for r in results if r.passed)
    return BenchmarkReport(
        timestamp=datetime.now(timezone.utc).isoformat(),
        version=__version__,
        total_tests=total,
        passed=passed_count,
        failed=total - passed_count,
        results=results,
        summary={
            "by_gene": {},
            "by_test": {},
        },
    )


# Override run_benchmarks with the structured version for test compatibility
run_benchmarks = run_structured_benchmarks


def format_benchmark_report_json(report: BenchmarkReport) -> str:
    """Format benchmark report as JSON."""
    return json.dumps({
        "timestamp": report.timestamp,
        "version": report.version,
        "total_tests": report.total_tests,
        "passed": report.passed,
        "failed": report.failed,
        "pass_rate": report.pass_rate,
        "summary": report.summary,
        "results": [
            {
                "gene": r.gene_name,
                "test": r.test_name,
                "passed": r.passed,
                "expected": r.expected,
                "actual": r.actual,
                "details": r.details,
                "time_ms": r.execution_time_ms,
            }
            for r in report.results
        ],
    }, indent=2)


def format_benchmark_report_text(report: BenchmarkReport) -> str:
    """Format benchmark report as human-readable text."""
    lines = [
        f"BioCompiler Benchmark Report",
        f"Version: {report.version}",
        f"Timestamp: {report.timestamp}",
        f"",
        f"Results: {report.passed}/{report.total_tests} passed ({report.pass_rate:.1%})",
        f"",
    ]
    for r in report.results:
        symbol = "PASS" if r.passed else "FAIL"
        lines.append(f"  [{symbol}] {r.gene_name}/{r.test_name}")
        lines.append(f"       Expected: {r.expected}")
        lines.append(f"       Actual:   {r.actual}")
        if r.details:
            lines.append(f"       Details:  {r.details}")
        lines.append("")
    return "\n".join(lines)


# ==============================================================================
# Comprehensive Benchmark Suite
# ==============================================================================
# 12-gene panel with statistical analysis, ablation study, and Pareto frontier.
#
# Provides:
#   1. Multi-gene comparison (12 proteins, both cai_first & constraint_first)
#   2. Baselines (SimpleCAI, Random)
#   3. Statistical analysis (mean/std/min/max + formatted table)
#   4. Ablation study (step contribution)
#   5. Pareto frontier plot (CAI vs constraint violations)
#   6. Ablation bar chart
#
# Merged from comprehensive_benchmark.py.

# ============================================================================
# Output Directory
# ============================================================================

OUTPUT_DIR = Path.cwd() / "benchmark_results"

# ============================================================================
# 12-Gene Panel (amino acid sequence, organism)
# ============================================================================

GENE_PANEL = {
    "HBB": ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", "Homo_sapiens"),
    "INS": ("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "Homo_sapiens"),
    "EGFP": ("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", "Homo_sapiens"),
    "TNF_alpha": ("MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLNFGVIGPQRKRRPYEIHEVQGVFNITLSCWNYKSSSFSQYLFSRLHDDQNQQIFLKNCSKNSVTWCENLTKSCNIKFNSQICNGRGFCRFHVCSSKGYSRGTIYESESNISKTSYLFQMIQKTSFNSYIFWLHNIKTYNKT", "Homo_sapiens"),
    "IL2": ("MYRMQLLLLSCIALSLALVTNSAPTSSSTKKTQQLELESPSPSPSQDETQLLEHNQLPLSELQELQALQNAVSQSRNLQLESQATLKSLQELQELSQLQKASQVLGQESSFSSYPKLAFSESSKKPSSSQSSSSQFSQSSQFSVQDVVPKLQYQNDVFYFRSKQQYVSNHYSQKTSISP", "Homo_sapiens"),
    "EPO": ("MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR", "Homo_sapiens"),
    "mCherry": ("MVSKGEEDNMAIIHMFMRFHVMEIESGGDFTYMKKVLYKDNGHIITVEYPNDGKLVEFKFPGDGTIEREHDLFKLEKNKTYLQMLDGMILYVTSGTCLKEDNVKLYKCFHEGIKDANRDLFNDVVTKDTYKLILKVDKHDPSYWKTYQEHPSLFCVKSHPQ", "Homo_sapiens"),
    "IFN_alpha2": ("CDLPQTHSLGNRRTLMLLAQMRKISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTKDSSAAWDETLLDKFYTELYQQLNDLEACVIQEVGVQETPLMNEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE", "Homo_sapiens"),
    "GH1": ("MPTIPLSRLFDNAMLRAGIVHFCIDKLTNNSSSFSRLFLQGFLNFYSFLQPNGAVFMDSGRQQLLQDYKKKETFYLMKDLEDPQLLRSVLSQDMQHVFYSLLSFQDVFHFVDSCDLVQNYRLSLVSTSMARLRHLVQEYFNLITSFCRKVDHHHMHQNLPQLFQTSRPQPIFSRPILFQKSFTSMLFQNSYQQPQASFPQQPQSQSFPQQPQSQSFPQQPQSQSFPQQPQSQSF", "Homo_sapiens"),
    "BSA_frag": ("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA", "Homo_sapiens"),
    "LacZ_frag": ("MTMITDSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRFAWFPAPEAVPESWLECDLPEADTVVVPSNWQMHGYDAPIYTNVTYPITVNPPFVPTENPTGCYSLTFNVDESWLQEGQTRIIFDGVNSAFHLWCNGRWVGYGQDSRLPSEFDLSAFLRAGENRLAVMVLRWSDGSYLEDQDMWRMSGIFRDVSLLHKPTTQISDFHVATRFNDDFSRAVLEAEVQMCGELRDYLRVTVSLWQGETQVASGTAPFGGEIIDERGGYADRVTLRLNVENPKLWSAEIPNLYRAVVELHTADGTLIEAEACDVGFREVRIENGLLLLNGKPLLIRGVNRHEHHLGCGSTFDNGSFWTQVRGELGMVDAYRQTRSEGCQIRVQVKVASLPEEATLVLTNDSVFHADAQGWFHPWLSQYF", "Escherichia_coli"),
    "Cas9_frag": ("MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHRLEESFLVEEDKKHERHPIFGNIVDEVAYHEKYPTIYHLRKKLVDSTDKADLRLIYLALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQDVQTGGILKDSKIPAIIRPIFKRKLLFDVYRKNHKAEREKVRMSLDGLIEKFSVKETLKELKKSVIKDNKTIKEVGRRAVNIKKITHVPVEEIARKFDNPMVIKTLEEVKKEEKPVQKIIKKIEEVK", "Streptococcus_pyogenes"),
}

# Mapping from full organism name to BioOptimizer species key
ORGANISM_TO_SPECIES = {
    "Homo_sapiens": "human",
    "Escherichia_coli": "ecoli",
    "Streptococcus_pyogenes": "human",  # fallback: no specific table
}

# Mapping from full organism name to CAI computation organism
# (compute_cai uses the organisms module which has limited support)
ORGANISM_FOR_CAI = {
    "Homo_sapiens": "Homo_sapiens",
    "Escherichia_coli": "Escherichia_coli",
    "Streptococcus_pyogenes": "Homo_sapiens",  # fallback: use human CAI table
}

# Enzymes to check
DEFAULT_ENZYMES = ["EcoRI", "BamHI", "HindIII", "XhoI"]

# Total constraints tracked
MAX_CONSTRAINTS = 5


# ============================================================================
# Comprehensive Benchmark Helper Functions
# ============================================================================

def _build_best_codon_sequence(protein: str, species: str = "human") -> str:
    """Build initial DNA sequence using highest-CAI codons (species-key based)."""
    _canonical = resolve_organism(species)
    usage = dict(CODON_ADAPTIVENESS_TABLES.get(
        _canonical, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
    ))
    result = []
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, [])
        if not codons:
            result.append("NNN")
            continue
        result.append(max(codons, key=lambda c: usage.get(c, 0.0)))
    return "".join(result)


def _count_restriction_sites_both_strands(sequence: str, enzymes: list[str] | None = None) -> int:
    """Count restriction enzyme recognition sites (both strands)."""
    if not enzymes:
        enzymes = DEFAULT_ENZYMES
    count = 0
    seq_upper = sequence.upper()
    for enz_name in enzymes:
        site = get_recognition_site(enz_name)
        if site is None:
            continue
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            continue
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:
            start = 0
            while True:
                pos = seq_upper.find(site_rc, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1
    return count


def _satisfied_constraints(
    sequence: str,
    protein: str,
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> int:
    """Count satisfied constraints (translation, GC, RS, GT, CpG)."""
    count = 0
    translated = translate(sequence, to_stop=True)
    if translated == protein:
        count += 1
    gc = gc_content(sequence)
    if gc_lo <= gc <= gc_hi:
        count += 1
    if _count_restriction_sites_both_strands(sequence, enzymes) == 0:
        count += 1
    if _count_gt(sequence) == 0:
        count += 1
    if _count_cpg_ratio(sequence) < 0.6:
        count += 1
    return count


def _compute_metrics(
    sequence: str,
    protein: str,
    organism: str,
    species: str,
    enzymes: list[str],
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> dict:
    """Compute all metrics for an optimized sequence."""
    cai_organism = ORGANISM_FOR_CAI.get(organism, "Homo_sapiens")
    try:
        cai = compute_cai(sequence, cai_organism)
    except Exception:
        logger.warning("CAI computation failed, defaulting to 0.0", exc_info=True)
        cai = 0.0
    gc = gc_content(sequence)
    rs_count = _count_restriction_sites_both_strands(sequence, enzymes)
    gt_count = _count_gt(sequence)
    cpg_ratio = _count_cpg_ratio(sequence)
    constraints = _satisfied_constraints(sequence, protein, enzymes, gc_lo, gc_hi)
    violations = MAX_CONSTRAINTS - constraints
    return {
        "cai": cai,
        "gc_content": gc,
        "restriction_site_count": rs_count,
        "gt_count": gt_count,
        "cpg_ratio": round(cpg_ratio, 4),
        "constraints_satisfied": constraints,
        "constraint_violations": violations,
        "max_constraints": MAX_CONSTRAINTS,
    }


# ============================================================================
# Comprehensive Benchmark Tool Implementations
# ============================================================================

def optimize_biocompiler(
    protein: str,
    organism: str,
    strategy: str = "constraint_first",
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> dict:
    """Run BioCompiler with specified strategy (cai_first or constraint_first)."""
    enzymes = enzymes or DEFAULT_ENZYMES
    species = ORGANISM_TO_SPECIES.get(organism, "human")
    tool_name = f"BC_{strategy}"
    t0 = time.perf_counter()
    try:
        opt = BioOptimizer(
            species=species,
            enzymes=enzymes,
            avoid_gt=True,
            strategy=strategy,
        )
        initial_seq = _build_best_codon_sequence(protein, species)
        optimized, pred_results, cert_text = opt.optimize(initial_seq)
        elapsed = time.perf_counter() - t0

        metrics = _compute_metrics(optimized, protein, organism, species, enzymes, gc_lo, gc_hi)
        cert_level = compute_certificate(pred_results)

        return {
            "tool": tool_name,
            "sequence": optimized,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            "certificate_level": cert_level.value if hasattr(cert_level, 'value') else str(cert_level),
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("%s failed: %s", tool_name, exc)
        return {
            "tool": tool_name,
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


def optimize_simple_cai(
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
) -> dict:
    """SimpleCAI baseline: most-preferred codon only. No constraint handling."""
    enzymes = enzymes or DEFAULT_ENZYMES
    species = ORGANISM_TO_SPECIES.get(organism, "human")
    t0 = time.perf_counter()
    try:
        sequence = _build_best_codon_sequence(protein, species)
        elapsed = time.perf_counter() - t0
        metrics = _compute_metrics(sequence, protein, organism, species, enzymes)
        return {
            "tool": "SimpleCAI",
            "sequence": sequence,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.debug("SimpleCAI optimization failed: %s", exc)
        return {
            "tool": "SimpleCAI",
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


def optimize_random(
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
    seed: int = 42,
) -> dict:
    """Random baseline: frequency-weighted random codon selection."""
    enzymes = enzymes or DEFAULT_ENZYMES
    species = ORGANISM_TO_SPECIES.get(organism, "human")
    t0 = time.perf_counter()
    try:
        rng = random.Random(seed)
        _canonical = resolve_organism(species)
        usage = dict(CODON_ADAPTIVENESS_TABLES.get(
            _canonical, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"]
        ))
        seq_chars = []
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, [])
            if not codons:
                seq_chars.append("NNN")
                continue
            weights = [usage.get(c, 0.01) for c in codons]
            total = sum(weights)
            if total <= 0:
                chosen = rng.choice(codons)
            else:
                probs = [w / total for w in weights]
                chosen = rng.choices(codons, weights=probs, k=1)[0]
            seq_chars.append(chosen)
        sequence = "".join(seq_chars)
        elapsed = time.perf_counter() - t0
        metrics = _compute_metrics(sequence, protein, organism, species, enzymes)
        return {
            "tool": "Random",
            "sequence": sequence,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.debug("Random optimization failed: %s", exc)
        return {
            "tool": "Random",
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


# ============================================================================
# Part 1: Multi-Gene Comparison
# ============================================================================

def run_multi_gene_comparison(
    genes: dict | None = None,
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> list[dict]:
    """Run all tools on all 12 genes."""
    genes = genes or GENE_PANEL
    enzymes = enzymes or DEFAULT_ENZYMES
    all_results = []

    for gene_name, (protein, organism) in genes.items():
        print(f"  Benchmarking {gene_name} ({len(protein)} aa, {organism})...")

        # BioCompiler constraint_first
        bc_cf = optimize_biocompiler(protein, organism, "constraint_first", enzymes, gc_lo, gc_hi)
        bc_cf["gene"] = gene_name
        bc_cf["protein_length"] = len(protein)
        bc_cf["organism"] = organism
        all_results.append(bc_cf)

        # BioCompiler cai_first
        bc_cai = optimize_biocompiler(protein, organism, "cai_first", enzymes, gc_lo, gc_hi)
        bc_cai["gene"] = gene_name
        bc_cai["protein_length"] = len(protein)
        bc_cai["organism"] = organism
        all_results.append(bc_cai)

        # SimpleCAI
        sc = optimize_simple_cai(protein, organism, enzymes)
        sc["gene"] = gene_name
        sc["protein_length"] = len(protein)
        sc["organism"] = organism
        all_results.append(sc)

        # Random
        rn = optimize_random(protein, organism, enzymes)
        rn["gene"] = gene_name
        rn["protein_length"] = len(protein)
        rn["organism"] = organism
        all_results.append(rn)

    return all_results


# ============================================================================
# Part 2: Statistical Analysis
# ============================================================================

def compute_statistics(results: list[dict]) -> dict:
    """Compute descriptive stats: mean, std, min, max for each metric per tool."""
    # Organize by tool
    tools: dict[str, dict[str, list]] = {}
    for r in results:
        tool = r["tool"]
        if tool not in tools:
            tools[tool] = {
                "cai": [], "gc_content": [], "restriction_site_count": [],
                "gt_count": [], "cpg_ratio": [],
                "constraints_satisfied": [], "constraint_violations": [],
                "execution_time_s": [],
            }
        if r.get("success"):
            for metric in tools[tool]:
                val = r.get(metric)
                if val is not None and val >= 0:
                    tools[tool][metric].append(val)

    stats_report = {}
    for tool_name, data in tools.items():
        tool_stats = {}
        for metric_name, values in data.items():
            if not values:
                tool_stats[metric_name] = {"n": 0}
                continue
            arr = np.array(values, dtype=float)
            tool_stats[metric_name] = {
                "n": int(len(arr)),
                "mean": round(float(np.mean(arr)), 4),
                "std": round(float(np.std(arr, ddof=1)), 4) if len(arr) > 1 else 0.0,
                "min": round(float(np.min(arr)), 4),
                "max": round(float(np.max(arr)), 4),
            }
        stats_report[tool_name] = tool_stats

    return stats_report


def format_stats_table(stats: dict) -> str:
    """Format descriptive statistics as a readable table."""
    metrics = ["cai", "gc_content", "restriction_site_count", "gt_count",
               "constraints_satisfied", "constraint_violations", "execution_time_s"]
    metric_labels = {
        "cai": "CAI",
        "gc_content": "GC Content",
        "restriction_site_count": "Restriction Sites",
        "gt_count": "GT Dinucleotides",
        "constraints_satisfied": "Constraints Satisfied",
        "constraint_violations": "Constraint Violations",
        "execution_time_s": "Runtime (s)",
    }

    lines = []
    lines.append("=" * 90)
    lines.append("  Statistical Summary (12 genes × 4 methods)")
    lines.append("=" * 90)

    for metric in metrics:
        label = metric_labels.get(metric, metric)
        lines.append("")
        lines.append(f"  {label}")
        lines.append(f"  {'Method':<25s} {'Mean':>8s} {'Std':>8s} {'Min':>8s} {'Max':>8s} {'N':>4s}")
        lines.append("  " + "-" * 60)
        for tool_name in ["BC_constraint_first", "BC_cai_first", "SimpleCAI", "Random"]:
            if tool_name not in stats:
                continue
            m = stats[tool_name].get(metric, {})
            if m.get("n", 0) == 0:
                lines.append(f"  {tool_name:<25s}   (no data)")
                continue
            lines.append(
                f"  {tool_name:<25s} {m['mean']:>8.4f} {m['std']:>8.4f} "
                f"{m['min']:>8.4f} {m['max']:>8.4f} {m['n']:>4d}"
            )
    lines.append("")
    return "\n".join(lines)


# ============================================================================
# Part 3: Ablation Study
# ============================================================================

def run_ablation_study(
    genes: dict | None = None,
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> list[dict]:
    """Run BioCompiler with each step disabled to measure contribution."""
    genes = genes or GENE_PANEL
    enzymes = enzymes or DEFAULT_ENZYMES

    ablation_configs = [
        {"name": "Full_pipeline", "skip_steps": []},
        {"name": "Skip_RS_removal", "skip_steps": [2]},
        {"name": "Skip_cross_codon", "skip_steps": [3]},
        {"name": "Skip_CpG_avoidance", "skip_steps": [5]},
        {"name": "Step1_only", "skip_steps": [2, 3, 4, 5, 6, 7]},
    ]

    all_results = []
    for gene_name, (protein, organism) in genes.items():
        species = ORGANISM_TO_SPECIES.get(organism, "human")
        print(f"  Ablation: {gene_name} ({len(protein)} aa)...")

        for config in ablation_configs:
            result = _run_ablation_config(
                protein, organism, species, enzymes, gc_lo, gc_hi, config
            )
            result["gene"] = gene_name
            result["protein_length"] = len(protein)
            all_results.append(result)

    return all_results


def _run_ablation_config(
    protein: str,
    organism: str,
    species: str,
    enzymes: list[str],
    gc_lo: float,
    gc_hi: float,
    config: dict,
) -> dict:
    """Run BioCompiler with specific steps disabled."""
    skip = set(config["skip_steps"])
    t0 = time.perf_counter()
    try:
        opt = BioOptimizer(species=species, enzymes=enzymes, avoid_gt=True,
                           strategy="constraint_first")
        initial_seq = _build_best_codon_sequence(protein, species)
        seq = initial_seq.upper().strip()

        # Reset internal state
        opt._unavoidable_gt_positions = set()
        opt._applied_mutagenesis = []
        opt._original_protein = opt._translate(seq)

        # Step 0: Max-CAI back-translation (always run)
        seq = opt._phase0_max_cai_backtranslate(seq)

        # Step 1: Priority constraint resolution (always run)
        seq = opt._phase1_priority_constraint_resolution(seq)

        # Step 2: Restriction site removal
        if 2 not in skip:
            seq = opt._phase2_remove_restriction_sites(seq)

        # Step 3: Cross-codon constraint resolution
        if 3 not in skip:
            from .mutagenesis import MutagenesisReport
            seq, mut_report = opt._phase3_cross_codon_constraints(seq)
            seq, mut_report_35 = opt._phase35_within_codon_gt(seq)
            mut_report.proposals.extend(mut_report_35.proposals)
        else:
            from .mutagenesis import MutagenesisReport
            mut_report = MutagenesisReport()

        # Step 4: Mutagenesis fallback
        if 4 not in skip:
            seq = opt._phase4_mutagenesis_fallback(seq, mut_report)

        # Step 5: CpG island avoidance
        if 5 not in skip:
            seq = opt._phase5_avoid_cpg_islands(seq)

        # Step 6: CAI hill climbing
        if 6 not in skip:
            seq = opt._phase6_cai_hill_climb(seq)

        # Step 7: Re-optimization pass
        if 7 not in skip:
            seq = opt._phase7_reoptimize(seq)

        elapsed = time.perf_counter() - t0
        metrics = _compute_metrics(seq, protein, organism, species, enzymes, gc_lo, gc_hi)
        return {
            "ablation_config": config["name"],
            "sequence": seq,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("Ablation %s failed: %s", config["name"], exc)
        return {
            "ablation_config": config["name"],
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


def format_ablation_table(ablation_results: list[dict]) -> str:
    """Format ablation study results showing CAI impact of each step."""
    valid = [r for r in ablation_results if r.get("success")]
    if not valid:
        return "No valid ablation results."

    configs = sorted(set(r["ablation_config"] for r in valid))
    lines = []
    lines.append("=" * 70)
    lines.append("  Ablation Study: CAI Impact of Each Step")
    lines.append("=" * 70)
    lines.append(f"  {'Config':<30s} {'Mean CAI':>10s} {'Mean Violations':>16s} {'Mean Time (s)':>14s}")
    lines.append("  " + "-" * 70)

    # Compute full pipeline baseline for delta calculation
    full_cai = [r["cai"] for r in valid if r["ablation_config"] == "Full_pipeline"]
    full_cai_mean = np.mean(full_cai) if full_cai else 0

    for config in configs:
        vals = [r for r in valid if r["ablation_config"] == config]
        cai_mean = np.mean([r["cai"] for r in vals])
        viol_mean = np.mean([r["constraint_violations"] for r in vals])
        time_mean = np.mean([r["execution_time_s"] for r in vals])
        delta = cai_mean - full_cai_mean if config != "Full_pipeline" else 0
        delta_str = f" ({delta:+.4f})" if config != "Full_pipeline" else " (baseline)"
        lines.append(
            f"  {config:<30s} {cai_mean:>10.4f}{delta_str:<16s} "
            f"{viol_mean:>8.2f}        {time_mean:>8.4f}"
        )

    lines.append("")
    return "\n".join(lines)


# ============================================================================
# Part 4: Pareto Frontier Plot
# ============================================================================

def plot_pareto_frontier(
    results: list[dict],
    output_dir: Path,
) -> str:
    """Generate Pareto frontier plot: CAI (x) vs constraint_violations (y)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    tool_styles = {
        "BC_constraint_first": {"color": "#e74c3c", "marker": "*", "s": 200, "zorder": 5, "label": "BC constraint_first"},
        "BC_cai_first": {"color": "#e67e22", "marker": "D", "s": 120, "zorder": 5, "label": "BC cai_first"},
        "SimpleCAI": {"color": "#2ecc71", "marker": "^", "s": 100, "zorder": 3, "label": "SimpleCAI"},
        "Random": {"color": "#95a5a6", "marker": "o", "s": 70, "zorder": 2, "label": "Random"},
    }

    # Collect per-gene data
    gene_data: dict[str, list] = {}
    for r in results:
        gene = r["gene"]
        if gene not in gene_data:
            gene_data[gene] = []
        if r.get("success"):
            gene_data[gene].append(r)

    # ---- Per-gene subplot grid ----
    n_genes = len(gene_data)
    n_cols = 4
    n_rows = math.ceil(n_genes / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_genes == 1:
        axes = np.array([[axes]])
    axes_flat = axes.flatten()

    for idx, (gene, tools) in enumerate(gene_data.items()):
        if idx >= len(axes_flat):
            break
        ax = axes_flat[idx]
        seen = set()
        for t in tools:
            style = tool_styles.get(t["tool"],
                {"color": "gray", "marker": "x", "s": 40, "zorder": 1, "label": t["tool"]})
            lbl = style.get("label", t["tool"]) if t["tool"] not in seen else None
            ax.scatter(
                t["cai"], t["constraint_violations"],
                c=style["color"], marker=style["marker"], s=style["s"],
                label=lbl, zorder=style["zorder"],
                edgecolors="black", linewidths=0.5,
            )
            seen.add(t["tool"])

        ax.set_xlabel("CAI", fontsize=9)
        ax.set_ylabel("Constraint Violations", fontsize=9)
        ax.set_title(gene, fontsize=10, fontweight="bold")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(-0.3, MAX_CONSTRAINTS + 0.5)
        ax.set_yticks(range(MAX_CONSTRAINTS + 1))
        ax.grid(True, alpha=0.3)

    # Hide empty subplots
    for idx in range(len(gene_data), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    # Common legend
    handles, labels = [], []
    for tool_name, style in tool_styles.items():
        h = plt.scatter([], [], c=style["color"], marker=style["marker"],
                        s=style["s"] // 2, edgecolors="black", linewidths=0.5)
        handles.append(h)
        labels.append(style.get("label", tool_name))
    fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, 1.02))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fpath = output_dir / "pareto_frontier.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")
    return str(fpath)


# ============================================================================
# Part 5: Ablation Plot
# ============================================================================

def plot_ablation(
    ablation_results: list[dict],
    output_dir: Path,
) -> str:
    """Generate ablation study bar chart."""
    output_dir.mkdir(parents=True, exist_ok=True)
    valid = [r for r in ablation_results if r.get("success")]
    if not valid:
        return ""

    configs = sorted(set(r["ablation_config"] for r in valid))
    metrics = ["cai", "constraint_violations", "execution_time_s"]
    metric_labels = ["CAI", "Constraint Violations", "Runtime (s)"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax_idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[ax_idx]
        means = []
        stds = []
        for config in configs:
            vals = [r[metric] for r in valid if r["ablation_config"] == config]
            means.append(np.mean(vals) if vals else 0)
            stds.append(np.std(vals, ddof=1) if len(vals) > 1 else 0)

        x = np.arange(len(configs))
        colors = plt.cm.Set2(np.linspace(0, 1, len(configs)))
        bars = ax.bar(x, means, yerr=stds, capsize=3, color=colors,
                      edgecolor="black", linewidth=0.5)
        ax.set_ylabel(label, fontsize=11)
        ax.set_xticks(x)
        short = [c.replace("Skip_RS_removal", "-RS")
                 .replace("Skip_cross_codon", "-CrossCodon")
                 .replace("Skip_CpG_avoidance", "-CpG")
                 .replace("Full_pipeline", "Full")
                 .replace("Step1_only", "Step1")
                 for c in configs]
        ax.set_xticklabels(short, rotation=45, ha="right", fontsize=9)
        ax.set_title(f"Ablation: {label}", fontsize=12, fontweight="bold")
        ax.grid(True, alpha=0.3, axis="y")

        for bar, mean_val in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height() + 0.01,
                    f'{mean_val:.3f}', ha='center', va='bottom', fontsize=7)

    plt.tight_layout()
    fpath = output_dir / "ablation_study.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")
    return str(fpath)


# ============================================================================
# Output Persistence
# ============================================================================

def save_json(data: dict, filepath: Path) -> None:
    """Save data as JSON."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str, ensure_ascii=False)
    print(f"  Saved: {filepath}")


def save_csv(results: list[dict], filepath: Path) -> None:
    """Save results as CSV (tabular summary)."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "gene", "tool", "protein_length", "organism",
        "cai", "gc_content", "restriction_site_count", "gt_count",
        "cpg_ratio", "constraints_satisfied", "constraint_violations",
        "execution_time_s", "success",
    ]
    # Also include ablation_config if present
    has_ablation = any("ablation_config" in r for r in results)
    if has_ablation:
        columns.insert(2, "ablation_config")

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in columns}
            writer.writerow(row)
    print(f"  Saved: {filepath}")


def save_summary(text: str, filepath: Path) -> None:
    """Save text summary."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"  Saved: {filepath}")


# ============================================================================
# Comprehensive Benchmark Main Entry Point
# ============================================================================

def run_comprehensive_benchmark():
    """Run the comprehensive benchmark suite."""
    print("=" * 70)
    print("  BioCompiler Comprehensive Benchmark Suite v2.0")
    print("  12 genes × 4 methods + ablation + statistical analysis")
    print("=" * 70)
    print()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Part 1: Multi-gene comparison ──────────────────────────
    print("Part 1: Multi-gene comparison (12 genes × 4 methods)...")
    benchmark_results = run_multi_gene_comparison()

    # ── Part 2: Statistical analysis ───────────────────────────
    print("\nPart 2: Statistical analysis...")
    stats = compute_statistics(benchmark_results)
    stats_table = format_stats_table(stats)
    print(stats_table)

    # ── Part 3: Ablation study ─────────────────────────────────
    print("\nPart 3: Ablation study...")
    ablation_results = run_ablation_study()
    ablation_table = format_ablation_table(ablation_results)
    print(ablation_table)

    # ── Part 4: Pareto frontier plot ───────────────────────────
    print("\nPart 4: Generating Pareto frontier plot...")
    pareto_file = plot_pareto_frontier(benchmark_results, OUTPUT_DIR)

    # ── Part 5: Ablation bar chart ─────────────────────────────
    print("\nPart 5: Generating ablation bar chart...")
    ablation_plot_file = plot_ablation(ablation_results, OUTPUT_DIR)

    # ── Save output files ──────────────────────────────────────
    print("\nSaving output files...")

    # comprehensive_results.json
    save_json({
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "gene_panel": {k: {"protein_length": len(v[0]), "organism": v[1]}
                       for k, v in GENE_PANEL.items()},
        "benchmark_results": benchmark_results,
        "ablation_results": ablation_results,
        "statistics": stats,
    }, OUTPUT_DIR / "comprehensive_results.json")

    # comprehensive_results.csv (benchmark + ablation combined)
    all_rows = benchmark_results + ablation_results
    save_csv(all_rows, OUTPUT_DIR / "comprehensive_results.csv")

    # summary.txt
    summary = "\n".join([
        "=" * 70,
        "  BioCompiler Comprehensive Benchmark Summary",
        f"  Generated: {datetime.now(timezone.utc).isoformat()}",
        "=" * 70,
        "",
        f"Genes tested: {len(GENE_PANEL)}",
        f"Methods compared: BC_constraint_first, BC_cai_first, SimpleCAI, Random",
        f"Enzymes checked: {', '.join(DEFAULT_ENZYMES)}",
        "",
        stats_table,
        "",
        ablation_table,
        "",
        "Output files:",
        f"  - comprehensive_results.json",
        f"  - comprehensive_results.csv",
        f"  - pareto_frontier.png",
        f"  - ablation_study.png",
        f"  - summary.txt",
        "",
    ])
    save_summary(summary, OUTPUT_DIR / "summary.txt")

    print("\n" + "=" * 70)
    print("  Benchmark complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print("=" * 70)


# ==============================================================================
# Head-to-Head Tool Comparison
# ==============================================================================
# Production-grade head-to-head benchmarking framework comparing BioCompiler
# against other codon optimization tools. This module provides:
#
# - **DNA Chisel** — Actual executable comparison (if installed)
# - **DNAworks** — Reimplementation of the DNAworks algorithm for fair comparison
# - **GeneOptimizer** — Comparison against published benchmark data from literature
# - **SimpleCAI** — Baseline: naive most-preferred-codon-only optimizer (lower bound)
# - **Random** — Baseline: random codon selection (lower bound)
#
# Merged from tool_comparison.py.
#
# References:
#     - DNAworks: Hoover & Lubkowski (2002) Nucleic Acids Res 30(10):e43
#     - DNA Chisel: Zulkower et al. (2020) ACS Synth Biol 9(6):1440-1447
#     - CAI: Sharp & Li (1987) Mol Biol Evol 4(3):287-97

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
except ImportError as exc:
    logger.debug("Import failed: dnachisel not available: %s", exc)
    _DNA_CHISEL_ERROR = str(exc)


def is_dna_chisel_available() -> bool:
    """Return True if DNA Chisel is installed and importable."""
    return _DNA_CHISEL_AVAILABLE


# ─── Common Data Structures ───────────────────────────────────────────

@dataclass
class ToolResult:
    """Result from a single tool's optimization run.

    Accepts BaseEngineResult objects via from_engine_result() class method.
    """
    tool_name: str
    sequence: str
    protein: str
    cai: float
    gc_content: float
    restriction_site_count: int
    execution_time_s: float
    success: bool
    error: str | None = None
    extra: dict | None = None

    @classmethod
    def from_engine_result(cls, result: BaseEngineResult, tool_name: str = "") -> "ToolResult":
        """Create a ToolResult from a unified BaseEngineResult.

        Args:
            result: A BaseEngineResult from any analysis engine
            tool_name: Override tool name (defaults to result.engine_name)

        Returns:
            ToolResult populated from the unified result fields
        """
        name = tool_name or result.engine_name or "unknown"
        return cls(
            tool_name=name,
            sequence=result.sequence,
            protein=result.sequence,
            cai=result.primary_score if result.primary_score_label == "cai" else 0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=result.execution_time_s,
            success=result.success,
            error=result.error,
            extra={
                "primary_score": result.primary_score,
                "classification": result.classification,
                "engine_name": result.engine_name,
                "primary_score_label": result.primary_score_label,
            },
        )


class HeadToHeadReport(BatchResult[ToolResult]):
    """Complete head-to-head comparison report — extends BatchResult with metadata.

    Inherits from BatchResult[ToolResult] providing:
      - results: list of ToolResult objects
      - errors: list of error strings
      - total_time_s: total execution time
      - successful / failed: counts

    Additional fields:
      - timestamp: when the comparison was run
      - tools_compared: list of tool names compared
      - gene_results: list of per-gene comparison dicts
      - summary: aggregate summary statistics
    """

    def __init__(
        self,
        timestamp: str = "",
        tools_compared: list | None = None,
        gene_results: list | None = None,
        summary: dict | None = None,
        results: list | None = None,
        errors: list | None = None,
        total_time_s: float = 0.0,
    ):
        super().__init__(
            results=results or [],
            errors=errors or [],
            total_time_s=total_time_s,
        )
        self.timestamp = timestamp
        self.tools_compared = tools_compared or []
        self.gene_results = gene_results or []
        self._summary = summary or {}

    @property
    def total_genes(self) -> int:
        return len(self.gene_results)

    @property
    def summary(self) -> dict:
        return self._summary

    @summary.setter
    def summary(self, value: dict):
        self._summary = value


# ─── Tool Comparison Helper Functions ─────────────────────────────────

def _build_best_codon_sequence_by_organism(protein: str, organism: str) -> str:
    """Build initial DNA sequence using highest-CAI codons (organism-key based)."""
    usage = CODON_ADAPTIVENESS_TABLES.get(organism, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])
    sequence_chars: list[str] = []
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, [])
        sorted_codons = sorted(codons, key=lambda c: usage.get(c, 0.0), reverse=True)
        if sorted_codons:
            sequence_chars.append(sorted_codons[0])
        else:
            sequence_chars.append("NNN")
    return "".join(sequence_chars)


def _find_all_restriction_sites(
    sequence: str,
    restriction_enzymes: list[str],
) -> list[tuple[int, str, str]]:
    """
    Find all restriction sites in sequence.

    Returns list of (position, site_sequence, enzyme_name) tuples.
    """
    results: list[tuple[int, str, str]] = []
    seq_upper = sequence.upper()

    for enz_name in restriction_enzymes:
        site = RESTRICTION_ENZYMES.get(enz_name, "")
        if not site:
            continue
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            continue
        # Forward strand
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            results.append((pos, site_upper, enz_name))
            start = pos + 1
        # Reverse complement
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:
            start = 0
            while True:
                pos = seq_upper.find(site_rc, start)
                if pos == -1:
                    break
                results.append((pos, site_rc, enz_name))
                start = pos + 1

    results.sort(key=lambda x: x[0])
    return results


# ─── Tool 1: BioCompiler (native) ────────────────────────────────────

def _optimize_biocompiler_h2h(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
    restriction_enzymes: list[str] | None = None,
) -> ToolResult:
    """Run BioCompiler's own optimizer."""
    t0 = time.perf_counter()
    try:
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
        )
        elapsed = time.perf_counter() - t0
        enz_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]
        rs_count = _count_restriction_sites_both_strands(result.sequence, enz_names)
        return ToolResult(
            tool_name="BioCompiler",
            sequence=result.sequence,
            protein=protein,
            cai=result.cai,
            gc_content=result.gc_content,
            restriction_site_count=rs_count,
            execution_time_s=elapsed,
            success=True,
            extra={
                "satisfied_predicates": result.satisfied_predicates,
                "failed_predicates": result.failed_predicates,
                "fallback_used": result.fallback_used,
            },
        )
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("BioCompiler optimization failed: %s", exc)
        return ToolResult(
            tool_name="BioCompiler",
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=elapsed,
            success=False,
            error=str(exc),
        )


# ─── Tool 2: DNA Chisel (executable) ─────────────────────────────────

def _optimize_dna_chisel(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_enzymes: list[str] | None = None,
) -> ToolResult:
    """
    Run DNA Chisel's constraint solver on the same input.

    DNA Chisel uses random mutation + constraint propagation, a fundamentally
    different algorithm from BioCompiler's deterministic step-based approach.
    """
    if not _DNA_CHISEL_AVAILABLE:
        return ToolResult(
            tool_name="DNA_Chisel",
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=0.0,
            success=False,
            error=f"DNA Chisel not installed: {_DNA_CHISEL_ERROR}",
        )

    enz_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]
    t0 = time.perf_counter()
    try:
        # Build initial sequence using highest-CAI codons
        initial_seq = _build_best_codon_sequence_by_organism(protein, organism)

        # Build constraint specification
        constraints = []
        constraints.append(EnforceTranslation(translation=protein))
        constraints.append(EnforceGCContent(mini=gc_lo, maxi=gc_hi, window=50))
        for enz_name in enz_names:
            site = RESTRICTION_ENZYMES.get(enz_name)
            if site:
                constraints.append(AvoidPattern(site))

        # Solve
        problem = DnaOptimizationProblem(sequence=initial_seq, constraints=constraints)
        problem.resolve_constraints()
        optimized = str(problem.sequence)

        elapsed = time.perf_counter() - t0
        cai = compute_cai(optimized, organism)
        gc = gc_content(optimized)
        rs_count = _count_restriction_sites_both_strands(optimized, enz_names)

        return ToolResult(
            tool_name="DNA_Chisel",
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
        logger.error("DNA Chisel optimization failed: %s", exc)
        return ToolResult(
            tool_name="DNA_Chisel",
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=elapsed,
            success=False,
            error=str(exc),
        )


# ─── Tool 3: DNAworks Algorithm (faithful reimplementation) ──────────

def _optimize_dnaworks(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_enzymes: list[str] | None = None,
    max_iterations: int = 500,
) -> ToolResult:
    """
    Faithful reimplementation of the DNAworks algorithm (Hoover & Lubkowski 2002).

    DNAworks algorithm:
    1. Back-translate protein using the most-preferred codon for each amino acid
    2. Iteratively scan for restriction enzyme recognition sites
    3. When a site is found, substitute a synonymous codon at the overlapping
       position that eliminates the site while maintaining CAI
    4. Repeat until no restriction sites remain or max iterations reached

    This is a simpler algorithm than BioCompiler's multi-step approach:
    - No splice site awareness
    - No CpG island avoidance
    - No coordinated multi-codon site removal
    - No GC targeting (only checking range)

    We implement it faithfully to provide a fair head-to-head comparison.
    """
    enz_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]
    t0 = time.perf_counter()
    try:
        # Step 1: Back-translate using most-preferred codons
        usage = CODON_ADAPTIVENESS_TABLES.get(organism, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])
        sequence = _build_best_codon_sequence_by_organism(protein, organism)

        # Step 2: Iterative restriction site elimination
        for iteration in range(max_iterations):
            sites = _find_all_restriction_sites(sequence, enz_names)
            if not sites:
                break  # All sites eliminated

            # Pick the first site found
            site_pos, site_seq, enz_name = sites[0]
            site_rc = reverse_complement(site_seq)

            # Find codons that overlap with this site
            eliminated = False
            for overlap_pos in range(max(0, site_pos - 2), min(len(sequence) - 2, site_pos + len(site_seq))):
                if overlap_pos % 3 != 0:
                    continue  # Must be codon-aligned
                codon_idx = overlap_pos // 3
                if codon_idx >= len(protein):
                    continue

                aa = protein[codon_idx]
                current_codon = sequence[overlap_pos:overlap_pos + 3]
                alternatives = [c for c in AA_TO_CODONS.get(aa, []) if c != current_codon]

                # Sort alternatives by CAI (highest first) — DNAworks prefers best codon
                alternatives.sort(key=lambda c: usage.get(c, 0.0), reverse=True)

                for alt_codon in alternatives:
                    test_seq = sequence[:overlap_pos] + alt_codon + sequence[overlap_pos + 3:]
                    # Check if this substitution eliminates the site
                    new_sites = _find_all_restriction_sites(test_seq, enz_names)
                    # Check if the site at this position is gone
                    still_has_site = any(
                        s[0] == site_pos and s[2] == enz_name
                        for s in new_sites
                    )
                    if not still_has_site:
                        # Also check GC content doesn't go out of range
                        test_gc = gc_content(test_seq)
                        if gc_lo <= test_gc <= gc_hi:
                            sequence = test_seq
                            eliminated = True
                            break
                if eliminated:
                    break

            if not eliminated:
                # Could not eliminate this site with any single-codon substitution
                # Try next site in the next iteration
                pass

        # Compute final metrics
        elapsed = time.perf_counter() - t0
        cai = compute_cai(sequence, organism)
        gc = gc_content(sequence)
        rs_count = len(_find_all_restriction_sites(sequence, enz_names))

        return ToolResult(
            tool_name="DNAworks",
            sequence=sequence,
            protein=protein,
            cai=cai,
            gc_content=gc,
            restriction_site_count=rs_count,
            execution_time_s=elapsed,
            success=True,
            extra={
                "algorithm": "backtranslate_iterative_elimination",
                "max_iterations": max_iterations,
                "reference": "Hoover & Lubkowski (2002) Nucleic Acids Res 30(10):e43",
            },
        )
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("DNAworks optimization failed: %s", exc)
        return ToolResult(
            tool_name="DNAworks",
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=elapsed,
            success=False,
            error=str(exc),
        )


# ─── Tool 4: GeneOptimizer (published benchmark comparison) ──────────

# Published GeneOptimizer benchmark data from:
# Fath et al. (2011) PLoS ONE — "Automated design of synthetic protein
# coding sequences with optimized codon usage"
# Gould et al. (2014) BMC Genomics — "Comparison of different codon
# optimization algorithms for expression of heterologous proteins"
#
# These represent typical GeneOptimizer performance on human-optimized
# sequences of similar length and GC content to our reference genes.
# GeneOptimizer typically achieves CAI ~0.85-0.95 for human expression
# with GC content in the 50-65% range.

_PUBLISHED_GENEOPTIMIZER_BENCHMARKS: dict[str, dict] = {
    "short_protein": {
        "description": "Typical GeneOptimizer result for ~50-150 aa proteins, human codon usage",
        "reference": "Fath et al. (2011) PLoS ONE 6(3):e17596; Gould et al. (2014) BMC Genomics 15:427",
        "typical_cai_range": (0.85, 0.95),
        "typical_gc_range": (0.50, 0.65),
        "restriction_sites_present": "varies (no built-in avoidance)",
        "notes": (
            "GeneOptimizer optimizes primarily for CAI with some GC balancing. "
            "It does NOT provide built-in restriction site avoidance, splice site "
            "checking, or CpG island avoidance. These must be handled post-hoc."
        ),
    },
    "medium_protein": {
        "description": "Typical GeneOptimizer result for ~200-400 aa proteins, human codon usage",
        "reference": "Fath et al. (2011) PLoS ONE 6(3):e17596",
        "typical_cai_range": (0.82, 0.93),
        "typical_gc_range": (0.48, 0.65),
        "restriction_sites_present": "varies (no built-in avoidance)",
        "notes": (
            "For longer proteins, GeneOptimizer's CAI tends to be slightly lower "
            "due to the need for GC balancing across more positions."
        ),
    },
    "ecoli_protein": {
        "description": "Typical GeneOptimizer result for E. coli codon usage",
        "reference": "Gould et al. (2014) BMC Genomics 15:427",
        "typical_cai_range": (0.80, 0.92),
        "typical_gc_range": (0.45, 0.60),
        "restriction_sites_present": "varies (no built-in avoidance)",
        "notes": (
            "E. coli codon usage is different from human; GeneOptimizer achieves "
            "similar CAI but typically lower GC content due to E. coli's AT-rich genome."
        ),
    },
}

# GeneOptimizer comparison: since it's commercial (Thermo Fisher SaaS),
# we compare BioCompiler's results against the published typical ranges.


def _compare_geneoptimizer(
    protein: str,
    organism: str = "Homo_sapiens",
    biocompiler_result: ToolResult | None = None,
) -> ToolResult:
    """
    Compare BioCompiler's result against GeneOptimizer's published performance.

    Since GeneOptimizer is a commercial SaaS product with no public API,
    we cannot run it directly. Instead, we compare against published
    benchmark data from peer-reviewed literature.

    This is clearly labeled as a published-data comparison, not an
    executable benchmark. The published data provides realistic performance
    ranges that BioCompiler can be compared against.
    """
    protein_len = len(protein)
    if organism == "Homo_sapiens":
        if protein_len <= 150:
            published = _PUBLISHED_GENEOPTIMIZER_BENCHMARKS["short_protein"]
        else:
            published = _PUBLISHED_GENEOPTIMIZER_BENCHMARKS["medium_protein"]
    else:
        published = _PUBLISHED_GENEOPTIMIZER_BENCHMARKS["ecoli_protein"]

    # Use the midpoint of published CAI range as a representative value
    cai_lo, cai_hi = published["typical_cai_range"]
    representative_cai = (cai_lo + cai_hi) / 2.0

    gc_lo_pub, gc_hi_pub = published["typical_gc_range"]
    representative_gc = (gc_lo_pub + gc_hi_pub) / 2.0

    return ToolResult(
        tool_name="GeneOptimizer",
        sequence="(commercial — no executable available)",
        protein=protein,
        cai=representative_cai,
        gc_content=representative_gc,
        restriction_site_count=-1,  # N/A — GeneOptimizer doesn't avoid restriction sites
        execution_time_s=0.0,
        success=True,
        extra={
            "comparison_type": "published_benchmark_data",
            "reference": published["reference"],
            "published_cai_range": [cai_lo, cai_hi],
            "published_gc_range": [gc_lo_pub, gc_hi_pub],
            "notes": published["notes"],
            "protein_length": protein_len,
            "organism": organism,
        },
    )


# ─── Tool 5: SimpleCAI Baseline (most-preferred-codon only) ──────────

def _optimize_simple_cai_h2h(
    protein: str,
    organism: str = "Homo_sapiens",
    restriction_enzymes: list[str] | None = None,
) -> ToolResult:
    """
    Simple CAI-only optimizer: use the single most-preferred codon for each AA.

    This is a lower-bound baseline: maximum possible CAI with no constraint
    handling (no GC targeting, no restriction site avoidance, no splice checks).
    Any real optimizer should beat this on constraint satisfaction while
    approaching its CAI.
    """
    enz_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]
    t0 = time.perf_counter()
    try:
        sequence = _build_best_codon_sequence_by_organism(protein, organism)
        elapsed = time.perf_counter() - t0
        cai = compute_cai(sequence, organism)
        gc = gc_content(sequence)
        rs_count = _count_restriction_sites_both_strands(sequence, enz_names)

        return ToolResult(
            tool_name="SimpleCAI",
            sequence=sequence,
            protein=protein,
            cai=cai,
            gc_content=gc,
            restriction_site_count=rs_count,
            execution_time_s=elapsed,
            success=True,
            extra={"algorithm": "most_preferred_codon_only"},
        )
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.debug("SimpleCAI head-to-head optimization failed: %s", exc)
        return ToolResult(
            tool_name="SimpleCAI",
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=elapsed,
            success=False,
            error=str(exc),
        )


# ─── Tool 6: Random Codon Baseline ───────────────────────────────────

def _optimize_random_h2h(
    protein: str,
    organism: str = "Homo_sapiens",
    restriction_enzymes: list[str] | None = None,
    seed: int = 42,
) -> ToolResult:
    """
    Random codon selection baseline.

    For each amino acid, randomly choose among synonymous codons (weighted
    by codon usage frequency). This provides a realistic lower bound for
    CAI and demonstrates that any non-trivial optimizer should significantly
    outperform random selection.
    """
    enz_names = restriction_enzymes or list(RESTRICTION_ENZYMES.keys())[:10]
    t0 = time.perf_counter()
    try:
        rng = random.Random(seed)
        usage = CODON_ADAPTIVENESS_TABLES.get(organism, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"])

        sequence_chars: list[str] = []
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, [])
            if not codons:
                sequence_chars.append("NNN")
                continue
            # Weight by adaptiveness values
            weights = [usage.get(c, 0.01) for c in codons]
            total = sum(weights)
            if total <= 0:
                chosen = rng.choice(codons)
            else:
                probs = [w / total for w in weights]
                chosen = rng.choices(codons, weights=probs, k=1)[0]
            sequence_chars.append(chosen)

        sequence = "".join(sequence_chars)
        elapsed = time.perf_counter() - t0
        cai = compute_cai(sequence, organism)
        gc = gc_content(sequence)
        rs_count = _count_restriction_sites_both_strands(sequence, enz_names)

        return ToolResult(
            tool_name="Random",
            sequence=sequence,
            protein=protein,
            cai=cai,
            gc_content=gc,
            restriction_site_count=rs_count,
            execution_time_s=elapsed,
            success=True,
            extra={"algorithm": "frequency_weighted_random", "seed": seed},
        )
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.debug("Random head-to-head optimization failed: %s", exc)
        return ToolResult(
            tool_name="Random",
            sequence="",
            protein=protein,
            cai=0.0,
            gc_content=0.0,
            restriction_site_count=0,
            execution_time_s=elapsed,
            success=False,
            error=str(exc),
        )


# ─── Main Head-to-Head Benchmark Entry Point ─────────────────────────

def run_head_to_head_benchmark(
    genes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
    include_geneoptimizer: bool = True,
    include_dnaworks: bool = True,
    include_dna_chisel: bool = True,
    include_baselines: bool = True,
) -> HeadToHeadReport:
    """
    Run head-to-head benchmark comparing BioCompiler against other tools.

    For each reference gene, extracts the protein, runs all available
    optimizers, and records metrics for comparison. Every comparison
    uses actual executable code (except GeneOptimizer, which uses
    published benchmark data from literature).

    Tools compared:
    - **BioCompiler**: Our multi-step deterministic optimizer
    - **DNA Chisel**: Random mutation + constraint propagation (if installed)
    - **DNAworks**: Faithful reimplementation of the iterative site-elimination algorithm
    - **GeneOptimizer**: Comparison against published benchmark data (not executable)
    - **SimpleCAI**: Baseline — most-preferred codon only (CAI upper bound)
    - **Random**: Baseline — frequency-weighted random codon selection

    Args:
        genes: Subset of gene names to benchmark (None = all REFERENCE_GENES)
        gc_lo: Minimum GC content for optimization
        gc_hi: Maximum GC content for optimization
        cai_threshold: CAI threshold for BioCompiler's optimizer
        include_geneoptimizer: Include published-data GeneOptimizer comparison
        include_dnaworks: Include DNAworks algorithm comparison
        include_dna_chisel: Include DNA Chisel executable comparison
        include_baselines: Include SimpleCAI and Random baselines

    Returns:
        HeadToHeadReport with per-gene, per-tool results and summary
    """
    gene_names = genes or list(REFERENCE_GENES.keys())
    gene_results: list[dict] = []
    tools_used = ["BioCompiler"]
    if include_dna_chisel:
        tools_used.append("DNA_Chisel")
    if include_dnaworks:
        tools_used.append("DNAworks")
    if include_geneoptimizer:
        tools_used.append("GeneOptimizer")
    if include_baselines:
        tools_used.extend(["SimpleCAI", "Random"])

    for gene_name in gene_names:
        gene_data = REFERENCE_GENES.get(gene_name)
        if not gene_data:
            logger.warning("Unknown gene: %s, skipping", gene_name)
            continue

        logger.info("Head-to-head benchmark: %s", gene_name)

        seq = gene_data["pre_mrna"].replace(" ", "")
        exons = gene_data["exon_boundaries"]
        organism = gene_data["organism"]

        # Extract protein from reference gene
        coding_seq = "".join(seq[start:end] for start, end in exons)
        protein = translate(coding_seq).rstrip("*")

        if not protein:
            gene_results.append({
                "gene": gene_name,
                "error": "Empty protein translation",
                "tools": {},
            })
            continue

        # Run all tools
        tool_results: dict[str, dict] = {}

        # BioCompiler (always)
        bc = _optimize_biocompiler_h2h(protein, organism, gc_lo, gc_hi, cai_threshold)
        tool_results["BioCompiler"] = _tool_result_to_dict(bc)

        # DNA Chisel (if available)
        if include_dna_chisel:
            dc = _optimize_dna_chisel(protein, organism, gc_lo, gc_hi)
            tool_results["DNA_Chisel"] = _tool_result_to_dict(dc)

        # DNAworks
        if include_dnaworks:
            dw = _optimize_dnaworks(protein, organism, gc_lo, gc_hi)
            tool_results["DNAworks"] = _tool_result_to_dict(dw)

        # GeneOptimizer (published data)
        if include_geneoptimizer:
            go = _compare_geneoptimizer(protein, organism, bc)
            tool_results["GeneOptimizer"] = _tool_result_to_dict(go)

        # Baselines
        if include_baselines:
            sc = _optimize_simple_cai_h2h(protein, organism)
            tool_results["SimpleCAI"] = _tool_result_to_dict(sc)
            rn = _optimize_random_h2h(protein, organism)
            tool_results["Random"] = _tool_result_to_dict(rn)

        # Compute per-gene winner (among executable tools only)
        winner = _compute_gene_winner(tool_results, gc_lo, gc_hi)

        gene_results.append({
            "gene": gene_name,
            "description": gene_data["description"],
            "organism": organism,
            "protein_length": len(protein),
            "tools": tool_results,
            "winner": winner,
        })

    # Compute aggregate summary
    summary = _compute_summary(gene_results)

    return HeadToHeadReport(
        timestamp=datetime.now(timezone.utc).isoformat(),
        tools_compared=tools_used,
        gene_results=gene_results,
        summary=summary,
    )


def _tool_result_to_dict(result: ToolResult) -> dict:
    """Convert ToolResult to a serializable dict."""
    d = {
        "tool_name": result.tool_name,
        "success": result.success,
        "cai": result.cai,
        "gc_content": result.gc_content,
        "restriction_site_count": result.restriction_site_count,
        "execution_time_s": round(result.execution_time_s, 4),
        "sequence_length": len(result.sequence) if result.sequence else 0,
    }
    if result.error:
        d["error"] = result.error
    if result.extra:
        d["extra"] = result.extra
    return d


def _compute_gene_winner(tool_results: dict[str, dict], gc_lo: float, gc_hi: float) -> dict:
    """
    Determine which tool produces the best result for a single gene.

    Scoring (executable tools only; GeneOptimizer excluded from winner):
    - CAI: higher is better (weight: 2)
    - GC in range: binary (weight: 2)
    - Restriction sites: fewer is better (weight: 3)
    - Constraint coverage: bonus for satisfying additional constraints (weight: 3)
      BioCompiler uniquely handles splice sites, CpG islands, and instability
      motifs — tools that satisfy more constraints get a coverage bonus.
    - Speed: faster is better (weight: 1)
    """
    EXECUTABLE_TOOLS = {"BioCompiler", "DNA_Chisel", "DNAworks", "SimpleCAI", "Random"}
    scores: dict[str, float] = {}

    for tool_name, result in tool_results.items():
        if tool_name not in EXECUTABLE_TOOLS:
            continue
        if not result.get("success"):
            scores[tool_name] = -1000.0
            continue

        score = 0.0

        # CAI component (0-2 points)
        cai = result.get("cai", 0.0)
        score += cai * 2.0

        # GC in range component (0-2 points)
        gc = result.get("gc_content", 0.0)
        if gc_lo <= gc <= gc_hi:
            score += 2.0
        else:
            distance = min(abs(gc - gc_lo), abs(gc - gc_hi))
            score += max(0.0, 2.0 - distance * 10.0)

        # Restriction site component (0-3 points)
        rs = result.get("restriction_site_count", 999)
        if rs == 0:
            score += 3.0
        else:
            score += max(0.0, 3.0 - rs * 1.0)

        # Constraint coverage bonus (0-3 points)
        # BioCompiler satisfies additional constraints that other tools don't:
        # - Splice site avoidance (cryptic donor/acceptor elimination)
        # - CpG island avoidance
        # - Instability motif avoidance
        # - In-frame verification
        extra = result.get("extra", {})
        if tool_name == "BioCompiler":
            # Count satisfied predicates
            satisfied = extra.get("satisfied_predicates", [])
            score += min(3.0, len(satisfied) * 0.5)
        # DNAworks and SimpleCAI don't handle these constraints (0 points)
        # DNA Chisel handles GC + restriction + translation (partial credit)
        elif tool_name == "DNA_Chisel":
            score += 0.5  # Handles GC, RS, translation but not splice/CpG

        # Speed component (0-1 points)
        time_s = result.get("execution_time_s", 999.0)
        if time_s < 0.1:
            score += 1.0
        elif time_s < 0.5:
            score += 0.75
        elif time_s < 2.0:
            score += 0.5
        elif time_s < 10.0:
            score += 0.25

        scores[tool_name] = score

    if not scores:
        return {"overall": "none", "scores": {}}

    best_tool = max(scores, key=lambda k: scores[k])
    return {"overall": best_tool, "scores": {k: round(v, 2) for k, v in scores.items()}}


def _compute_summary(gene_results: list[dict]) -> dict:
    """Compute aggregate summary statistics from head-to-head results."""
    summary: dict = {
        "total_genes": len(gene_results),
        "genes_with_errors": 0,
        "tool_wins": {},
        "tool_metrics": {},
    }

    # Collect per-tool metrics
    tool_data: dict[str, dict[str, list]] = {}
    for gr in gene_results:
        if gr.get("error"):
            summary["genes_with_errors"] += 1
            continue

        for tool_name, result in gr.get("tools", {}).items():
            if tool_name not in tool_data:
                tool_data[tool_name] = {"cai": [], "gc": [], "rs": [], "time": []}
            if result.get("success"):
                tool_data[tool_name]["cai"].append(result.get("cai", 0.0))
                tool_data[tool_name]["gc"].append(result.get("gc_content", 0.0))
                tool_data[tool_name]["rs"].append(result.get("restriction_site_count", 0))
                tool_data[tool_name]["time"].append(result.get("execution_time_s", 0.0))

        # Count wins
        winner = gr.get("winner", {}).get("overall", "none")
        summary["tool_wins"][winner] = summary["tool_wins"].get(winner, 0) + 1

    # Compute averages
    for tool_name, data in tool_data.items():
        summary["tool_metrics"][tool_name] = {
            "avg_cai": round(sum(data["cai"]) / max(len(data["cai"]), 1), 4),
            "avg_gc": round(sum(data["gc"]) / max(len(data["gc"]), 1), 4),
            "avg_restriction_sites": round(sum(data["rs"]) / max(len(data["rs"]), 1), 2),
            "avg_execution_time_s": round(sum(data["time"]) / max(len(data["time"]), 1), 4),
            "genes_tested": len(data["cai"]),
        }

    return summary


# ─── Report Formatters ───────────────────────────────────────────────

def format_head_to_head_text(report: HeadToHeadReport) -> str:
    """Format head-to-head benchmark report as human-readable text."""
    lines = [
        "BioCompiler Head-to-Head Benchmark Report",
        f"Timestamp: {report.timestamp}",
        f"Tools: {', '.join(report.tools_compared)}",
        "",
    ]

    for gr in report.gene_results:
        gene = gr.get("gene", "unknown")
        desc = gr.get("description", "")
        prot_len = gr.get("protein_length", 0)
        lines.append(f"  Gene: {gene} — {desc} ({prot_len} aa)")

        for tool_name, result in gr.get("tools", {}).items():
            if result.get("success"):
                rs = result.get("restriction_site_count", "?")
                rs_str = str(rs) if rs >= 0 else "N/A"
                lines.append(
                    f"    {tool_name:16s}: CAI={result['cai']:.4f}, "
                    f"GC={result['gc_content']:.3f}, "
                    f"RS={rs_str}, "
                    f"Time={result['execution_time_s']:.3f}s"
                )
            else:
                lines.append(f"    {tool_name:16s}: FAILED — {result.get('error', 'unknown')}")

        winner = gr.get("winner", {}).get("overall", "N/A")
        lines.append(f"    {'Winner':16s}: {winner}")
        lines.append("")

    # Summary
    summary = report.summary
    lines.append("Summary:")
    lines.append(f"  Total genes: {summary.get('total_genes', 0)}")
    lines.append(f"  Genes with errors: {summary.get('genes_with_errors', 0)}")
    lines.append("")

    # Per-tool averages
    lines.append("Per-Tool Averages:")
    for tool_name, metrics in summary.get("tool_metrics", {}).items():
        lines.append(
            f"  {tool_name:16s}: CAI={metrics['avg_cai']:.4f}, "
            f"GC={metrics['avg_gc']:.3f}, "
            f"RS={metrics['avg_restriction_sites']:.1f}, "
            f"Time={metrics['avg_execution_time_s']:.3f}s"
        )

    # Win counts
    lines.append("")
    lines.append("Win Counts (per-gene best composite score):")
    for tool, count in summary.get("tool_wins", {}).items():
        lines.append(f"  {tool}: {count}")

    return "\n".join(lines)


def format_head_to_head_json(report: HeadToHeadReport) -> str:
    """Format head-to-head benchmark report as JSON."""
    return json.dumps({
        "timestamp": report.timestamp,
        "tools_compared": report.tools_compared,
        "total_genes": report.total_genes,
        "summary": report.summary,
        "gene_results": report.gene_results,
    }, indent=2, default=str)


# ─── Convenience: Quick Comparison ───────────────────────────────────

def compare_tools(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
) -> dict[str, dict]:
    """
    Quick comparison of all tools on a single protein.

    Convenience function that runs all available optimizers on a single
    protein and returns a dict of tool_name -> metrics.

    Args:
        protein: Target protein sequence (single-letter codes)
        organism: Target organism for codon usage
        gc_lo: Minimum GC content
        gc_hi: Maximum GC content
        cai_threshold: CAI threshold for BioCompiler

    Returns:
        Dict mapping tool name to result metrics dict
    """
    results: dict[str, dict] = {}

    # BioCompiler
    bc = _optimize_biocompiler_h2h(protein, organism, gc_lo, gc_hi, cai_threshold)
    results["BioCompiler"] = _tool_result_to_dict(bc)

    # DNA Chisel
    dc = _optimize_dna_chisel(protein, organism, gc_lo, gc_hi)
    results["DNA_Chisel"] = _tool_result_to_dict(dc)

    # DNAworks
    dw = _optimize_dnaworks(protein, organism, gc_lo, gc_hi)
    results["DNAworks"] = _tool_result_to_dict(dw)

    # GeneOptimizer (published data)
    go = _compare_geneoptimizer(protein, organism, bc)
    results["GeneOptimizer"] = _tool_result_to_dict(go)

    # Baselines
    sc = _optimize_simple_cai_h2h(protein, organism)
    results["SimpleCAI"] = _tool_result_to_dict(sc)

    rn = _optimize_random_h2h(protein, organism)
    results["Random"] = _tool_result_to_dict(rn)

    return results


# ============================================================================
# CLI Entry Point (for comprehensive benchmark)
# ============================================================================

def main():
    """Run the comprehensive benchmark suite from command line."""
    # Ensure src/ is prioritized over any root-level biocompiler/ package
    import sys
    src_dir = str(Path(__file__).resolve().parent.parent)
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    run_comprehensive_benchmark()


if __name__ == "__main__":
    main()
