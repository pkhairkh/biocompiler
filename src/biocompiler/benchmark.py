"""
BioCompiler Benchmark Module v7.0.0
=====================================
Built-in benchmark sequences and performance measurement for BioCompiler.

Provides:
  - run_benchmark(): Optimize eGFP, mCherry, LacZ with human/ecoli and print a table
  - compare_tools(): Theoretical feature comparison table vs other tools
"""

import time
from typing import Dict, List, Tuple

try:
    from .optimizer import BioOptimizer
except ImportError:
    from .optimization import BioOptimizer
from .type_system import (
    CODON_TABLE,
    CertLevel,
    PredicateResult,
    check_no_avoidable_gt,
    check_no_cpg_island,
    check_no_restriction_site,
)
from .certificates import compute_certificate
from .species import SPECIES

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
    import math
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
    """Count total restriction sites in sequence."""
    from .restriction_sites import get_recognition_site
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
    print("  BioCompiler v7.0.0 — Built-in Benchmark")
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
            t0 = time.perf_counter()
            optimized, pred_results, cert_text = opt.optimize(gene_seq)
            t1 = time.perf_counter()
            elapsed = t1 - t0

            species_cai = SPECIES.get(species, SPECIES["ecoli"])

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


def compare_tools() -> None:
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
    print("  BioCompiler v7.0.0 — Tool Comparison (Theoretical)")
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

import json as _json
from dataclasses import dataclass as _dataclass, field as _field
from datetime import datetime as _dt, timezone as _tz


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


@_dataclass
class BenchmarkResult:
    """Result of a single benchmark test."""
    gene_name: str
    test_name: str
    passed: bool
    expected: str
    actual: str
    details: str | None = None
    execution_time_ms: float = 0.0


@_dataclass
class BenchmarkReport:
    """Complete benchmark report."""
    timestamp: str
    version: str
    total_tests: int
    passed: int
    failed: int
    results: list = _field(default_factory=list)
    summary: dict = _field(default_factory=dict)

    @property
    def pass_rate(self) -> float:
        return self.passed / max(self.total_tests, 1)


def run_structured_benchmarks(
    gene_names: list[str] | None = None,
    include_optimization: bool = True,
) -> BenchmarkReport:
    """Run structured benchmarks against known gene sets.

    This validates BioCompiler's predictions against biological ground truth.
    """
    from . import __version__
    from .scanner import gc_content
    from .translation import translate, compute_cai
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
        t0 = time.perf_counter()
        coding_seq = "".join(seq[start:end] for start, end in exons)
        protein = translate(coding_seq)
        protein_len = len(protein.rstrip("*"))
        expected_len = gene_data["known_protein_length"]
        passed = abs(protein_len - expected_len) <= 10
        results.append(BenchmarkResult(
            gene_name=gene_name, test_name="translation_length",
            passed=passed, expected=f"protein_length={expected_len}",
            actual=f"protein_length={protein_len}",
            execution_time_ms=(time.perf_counter() - t0) * 1000,
        ))

        # GC content benchmark
        t0 = time.perf_counter()
        gc = gc_content(seq)
        gc_lo, gc_hi = gene_data["expected_gc_range"]
        passed = gc_lo <= gc <= gc_hi
        results.append(BenchmarkResult(
            gene_name=gene_name, test_name="gc_content_range",
            passed=passed, expected=f"GC in [{gc_lo}, {gc_hi}]",
            actual=f"GC = {gc:.4f}",
            execution_time_ms=(time.perf_counter() - t0) * 1000,
        ))

        # CAI benchmark
        t0 = time.perf_counter()
        cai = compute_cai(coding_seq, organism)
        cai_lo, cai_hi = gene_data["expected_cai_range"]
        passed = cai_lo <= cai <= cai_hi
        results.append(BenchmarkResult(
            gene_name=gene_name, test_name="cai_range",
            passed=passed, expected=f"CAI in [{cai_lo}, {cai_hi}]",
            actual=f"CAI = {cai:.4f}",
            execution_time_ms=(time.perf_counter() - t0) * 1000,
        ))

        # Type predicates benchmark
        t0 = time.perf_counter()
        try:
            type_results = evaluate_all_predicates(
                seq=seq, known_exon_boundaries=exons, organism=organism,
            )
            n_pass = sum(1 for r in type_results if r.verdict.value in ("PASS", "LIKELY_PASS"))
            passed = n_pass >= 4
            results.append(BenchmarkResult(
                gene_name=gene_name, test_name="type_predicates",
                passed=passed, expected=">=4 predicates PASS",
                actual=f"PASS={n_pass}, total={len(type_results)}",
                details="; ".join(f"{r.predicate}={r.verdict.value}" for r in type_results),
                execution_time_ms=(time.perf_counter() - t0) * 1000,
            ))
        except Exception as e:
            results.append(BenchmarkResult(
                gene_name=gene_name, test_name="type_predicates",
                passed=False, expected="Predicates evaluated",
                actual=f"ERROR: {e}",
            ))

    total = len(results)
    passed_count = sum(1 for r in results if r.passed)
    return BenchmarkReport(
        timestamp=_dt.now(_tz.utc).isoformat(),
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
    return _json.dumps({
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
