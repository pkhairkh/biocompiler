"""BioCompiler Benchmarks Package.

Provides accuracy benchmarking against published experimental data and
validated ground-truth references for each BioCompiler module.

Public API::

    from biocompiler.benchmarks.benchmark_suite import (
        BenchmarkResult,
        BenchmarkReport,
        benchmark_thermal_stability,
        benchmark_rna_degradation,
        benchmark_dna_damage,
        benchmark_nucleosome,
        benchmark_ribosome,
        benchmark_mfe,
        benchmark_ligand_binding,
        run_all_benchmarks,
        print_benchmark_report,
        save_benchmark_report,
    )
"""

# Use lazy imports to avoid triggering biocompiler.__init__ circular imports.
# The benchmark functions themselves handle ImportError/AttributeError
# gracefully at call time, so deferring the actual module imports is safe.

def __getattr__(name: str):
    """Lazy re-export from benchmark_suite to avoid circular import."""
    from . import benchmark_suite as _bs
    if name in _bs.__all__:
        return getattr(_bs, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "BenchmarkResult",
    "BenchmarkReport",
    "benchmark_thermal_stability",
    "benchmark_rna_degradation",
    "benchmark_dna_damage",
    "benchmark_nucleosome",
    "benchmark_ribosome",
    "benchmark_mfe",
    "benchmark_ligand_binding",
    "run_all_benchmarks",
    "print_benchmark_report",
    "save_benchmark_report",
]
