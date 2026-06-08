"""
NUMBA vs Pure-Python Benchmarking Utility.

Compares the performance of NUMBA-accelerated kernels against their
pure-Python fallbacks for CAI computation, GC counting, dinucleotide
counting, sliding-window GC, and batch codon swap scoring.

Usage::

    from biocompiler.optimizer.benchmark_numba import run_numba_benchmark

    report = run_numba_benchmark(protein="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    print(report.summary())
"""

from __future__ import annotations

import logging
import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class KernelBenchmarkResult:
    """Result of benchmarking a single NUMBA kernel vs pure Python."""
    kernel_name: str
    numba_time_ms: float
    python_time_ms: float
    speedup: float
    numba_result_valid: bool
    results_match: bool
    n_iterations: int = 100

    def summary(self) -> str:
        match_str = "MATCH" if self.results_match else "MISMATCH"
        valid_str = "valid" if self.numba_result_valid else "INVALID"
        return (
            f"{self.kernel_name:30s}  "
            f"NUMBA={self.numba_time_ms:8.3f}ms  "
            f"Python={self.python_time_ms:8.3f}ms  "
            f"speedup={self.speedup:6.2f}x  "
            f"results={match_str}  NUMBA={valid_str}"
        )


@dataclass
class NumbaBenchmarkReport:
    """Full benchmark report comparing NUMBA vs pure-Python paths."""
    numba_available: bool
    use_numba: bool
    kernel_results: List[KernelBenchmarkResult] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"NUMBA Benchmark Report",
            f"  NUMBA available: {self.numba_available}",
            f"  USE_NUMBA flag:   {self.use_numba}",
            f"",
            f"  {'Kernel':30s}  {'NUMBA':>10s}  {'Python':>10s}  {'Speedup':>8s}  {'Results':>8s}",
            f"  {'-'*30}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*8}",
        ]
        for kr in self.kernel_results:
            lines.append(f"  {kr.summary()}")
        return "\n".join(lines)


def _bench_cai_computation(protein: str, organism: str, n_iter: int = 100) -> KernelBenchmarkResult:
    """Benchmark CAI computation: NUMBA kernel vs pure Python."""
    from ..translation import translate
    from ..organisms import CODON_ADAPTIVENESS_TABLES
    from ..type_system import AA_TO_CODONS, CODON_TABLE
    from .cai import _compute_cai_fast

    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
    if adaptiveness is None:
        return KernelBenchmarkResult("CAI computation", 0, 0, 0, False, False, n_iter)

    # Build a test DNA sequence using optimal codons
    dna_parts = []
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, ["ATG"])
        best = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
        dna_parts.append(best)
    seq = "".join(dna_parts)

    # Warm up
    _compute_cai_fast(seq, adaptiveness)

    # Time NUMBA path
    from ..numba_kernels import USE_NUMBA as orig_use_numba
    try:
        from ..numba_kernels import USE_NUMBA
        # Force NUMBA on
        import biocompiler.numba_kernels as _nk
        _nk.USE_NUMBA = True
        t0 = time.perf_counter()
        for _ in range(n_iter):
            numba_cai = _compute_cai_fast(seq, adaptiveness)
        numba_time = (time.perf_counter() - t0) / n_iter * 1000
        _nk.USE_NUMBA = orig_use_numba
    except Exception:
        numba_time = 0.0
        numba_cai = 0.0

    # Time pure-Python path
    import biocompiler.numba_kernels as _nk
    _nk.USE_NUMBA = False
    t0 = time.perf_counter()
    for _ in range(n_iter):
        py_cai = _compute_cai_fast(seq, adaptiveness)
    py_time = (time.perf_counter() - t0) / n_iter * 1000
    _nk.USE_NUMBA = orig_use_numba

    results_match = abs(numba_cai - py_cai) < 0.01
    numba_valid = numba_cai > 0.0
    speedup = py_time / numba_time if numba_time > 0 else float("inf")

    return KernelBenchmarkResult(
        "CAI computation",
        numba_time, py_time, speedup,
        numba_valid, results_match, n_iter,
    )


def _bench_gc_counting(seq: str, n_iter: int = 200) -> KernelBenchmarkResult:
    """Benchmark GC counting: NUMBA count_gc vs Python str.count."""
    from ..numba_kernels import count_gc, seq_to_bytes, USE_NUMBA as orig_use_numba, HAS_NUMBA

    if not HAS_NUMBA:
        return KernelBenchmarkResult("GC counting", 0, 0, 0, False, False, n_iter)

    seq_bytes = seq_to_bytes(seq)

    # Warm up NUMBA
    count_gc(seq_bytes)

    # Time NUMBA
    t0 = time.perf_counter()
    for _ in range(n_iter):
        numba_gc = count_gc(seq_bytes)
    numba_time = (time.perf_counter() - t0) / n_iter * 1000

    # Time Python
    t0 = time.perf_counter()
    for _ in range(n_iter):
        py_gc = seq.count("G") + seq.count("C")
    py_time = (time.perf_counter() - t0) / n_iter * 1000

    results_match = numba_gc == py_gc
    speedup = py_time / numba_time if numba_time > 0 else float("inf")

    return KernelBenchmarkResult(
        "GC counting",
        numba_time, py_time, speedup,
        numba_gc >= 0, results_match, n_iter,
    )


def _bench_dinucleotide_counting(seq: str, n_iter: int = 200) -> KernelBenchmarkResult:
    """Benchmark multi-dinucleotide counting: NUMBA vs Python."""
    from .cai import _count_dinucs_fast
    from ..numba_kernels import USE_NUMBA as orig_use_numba, HAS_NUMBA

    dinucs = ("GT", "CG", "AG")

    # Warm up
    _count_dinucs_fast(seq, *dinucs)

    # Time NUMBA path
    import biocompiler.numba_kernels as _nk
    _nk.USE_NUMBA = True
    t0 = time.perf_counter()
    for _ in range(n_iter):
        numba_counts = _count_dinucs_fast(seq, *dinucs)
    numba_time = (time.perf_counter() - t0) / n_iter * 1000
    _nk.USE_NUMBA = orig_use_numba

    # Time Python path
    _nk.USE_NUMBA = False
    t0 = time.perf_counter()
    for _ in range(n_iter):
        py_counts = _count_dinucs_fast(seq, *dinucs)
    py_time = (time.perf_counter() - t0) / n_iter * 1000
    _nk.USE_NUMBA = orig_use_numba

    results_match = numba_counts == py_counts
    speedup = py_time / numba_time if numba_time > 0 else float("inf")

    return KernelBenchmarkResult(
        "Dinucleotide counting (GT/CG/AG)",
        numba_time, py_time, speedup,
        True, results_match, n_iter,
    )


def _bench_gc_window(seq: str, n_iter: int = 100) -> KernelBenchmarkResult:
    """Benchmark sliding-window GC: NUMBA fast_gc_window vs Python."""
    from ..numba_kernels import fast_gc_window, seq_to_bytes, HAS_NUMBA, USE_NUMBA as orig_use_numba

    if not HAS_NUMBA:
        return KernelBenchmarkResult("GC window (sliding)", 0, 0, 0, False, False, n_iter)

    seq_bytes = seq_to_bytes(seq)
    window_size = 50

    # Warm up
    fast_gc_window(seq_bytes, window_size)

    # Time NUMBA
    t0 = time.perf_counter()
    for _ in range(n_iter):
        numba_result = fast_gc_window(seq_bytes, window_size)
    numba_time = (time.perf_counter() - t0) / n_iter * 1000

    # Time Python
    t0 = time.perf_counter()
    for _ in range(n_iter):
        n = len(seq)
        gc_count = sum(1 for b in seq[:window_size] if b in "GC")
        py_result = [gc_count / window_size]
        for i in range(1, n - window_size + 1):
            if seq[i - 1] in "GC":
                gc_count -= 1
            if seq[i + window_size - 1] in "GC":
                gc_count += 1
            py_result.append(gc_count / window_size)
    py_time = (time.perf_counter() - t0) / n_iter * 1000

    # Verify results match
    results_match = len(numba_result) == len(py_result) and all(
        abs(a - b) < 1e-10 for a, b in zip(numba_result, py_result)
    )
    speedup = py_time / numba_time if numba_time > 0 else float("inf")

    return KernelBenchmarkResult(
        "GC window (sliding)",
        numba_time, py_time, speedup,
        len(numba_result) > 0, results_match, n_iter,
    )


def _bench_batch_codon_swap(protein: str, organism: str, n_iter: int = 50) -> KernelBenchmarkResult:
    """Benchmark batch codon swap scoring: NUMBA vs Python."""
    from ..organisms import CODON_ADAPTIVENESS_TABLES
    from ..type_system import AA_TO_CODONS
    from .cai import _BatchSwapScorer
    from ..numba_kernels import USE_NUMBA as orig_use_numba

    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism)
    if adaptiveness is None:
        return KernelBenchmarkResult("Batch codon swap scoring", 0, 0, 0, False, False, n_iter)

    # Build test sequence
    dna_parts = []
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, ["ATG"])
        best = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
        dna_parts.append(best)
    seq = "".join(dna_parts)
    seq_codons = [seq[i:i+3] for i in range(0, len(seq), 3)]

    # Pick a position with alternatives
    test_pos = 1
    aa = protein[test_pos] if test_pos < len(protein) else protein[0]
    candidates = AA_TO_CODONS.get(aa, ["GCT"])

    # NUMBA path
    import biocompiler.numba_kernels as _nk
    scorer = _BatchSwapScorer(adaptiveness)
    scorer.reset_incremental_state(seq_codons)

    _nk.USE_NUMBA = True
    scorer._np_built = False  # Reset lazy build
    t0 = time.perf_counter()
    for _ in range(n_iter):
        numba_scores = scorer.score_candidates(seq_codons, test_pos, candidates)
    numba_time = (time.perf_counter() - t0) / n_iter * 1000
    _nk.USE_NUMBA = orig_use_numba

    # Python path
    scorer2 = _BatchSwapScorer(adaptiveness)
    scorer2.reset_incremental_state(seq_codons)
    _nk.USE_NUMBA = False
    scorer2._np_built = False
    t0 = time.perf_counter()
    for _ in range(n_iter):
        py_scores = scorer2.score_candidates(seq_codons, test_pos, candidates)
    py_time = (time.perf_counter() - t0) / n_iter * 1000
    _nk.USE_NUMBA = orig_use_numba

    results_match = len(numba_scores) == len(py_scores) and all(
        abs(a - b) < 0.01 for a, b in zip(numba_scores, py_scores)
    )
    speedup = py_time / numba_time if numba_time > 0 else float("inf")

    return KernelBenchmarkResult(
        "Batch codon swap scoring",
        numba_time, py_time, speedup,
        len(numba_scores) > 0, results_match, n_iter,
    )


def run_numba_benchmark(
    protein: str = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    organism: str = "Escherichia_coli",
    n_iter: int = 100,
) -> NumbaBenchmarkReport:
    """Run a comprehensive NUMBA vs Python benchmark.

    Benchmarks all NUMBA v2 kernels against their pure-Python fallbacks,
    verifying that results match and measuring speedup.

    Args:
        protein: Protein sequence to use for benchmarking.
        organism: Organism for codon adaptiveness tables.
        n_iter: Number of iterations per benchmark (higher = more stable).

    Returns:
        NumbaBenchmarkReport with detailed results per kernel.
    """
    from ..numba_kernels import HAS_NUMBA, USE_NUMBA
    from ..type_system import AA_TO_CODONS
    from ..organisms import CODON_ADAPTIVENESS_TABLES

    report = NumbaBenchmarkReport(
        numba_available=HAS_NUMBA,
        use_numba=USE_NUMBA,
    )

    if not HAS_NUMBA:
        logger.info("NUMBA is not installed; benchmark skipped")
        return report

    # Build test DNA sequence
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(organism, {})
    dna_parts = []
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, ["ATG"])
        best = max(codons, key=lambda c: adaptiveness.get(c, 0.0))
        dna_parts.append(best)
    seq = "".join(dna_parts)

    # Run benchmarks
    report.kernel_results.append(_bench_cai_computation(protein, organism, n_iter))
    report.kernel_results.append(_bench_gc_counting(seq, n_iter * 2))
    report.kernel_results.append(_bench_dinucleotide_counting(seq, n_iter * 2))
    report.kernel_results.append(_bench_gc_window(seq, n_iter))
    report.kernel_results.append(_bench_batch_codon_swap(protein, organism, max(10, n_iter // 2)))

    return report


__all__ = [
    "KernelBenchmarkResult",
    "NumbaBenchmarkReport",
    "run_numba_benchmark",
]
