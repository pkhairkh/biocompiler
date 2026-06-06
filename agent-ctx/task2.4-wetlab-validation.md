# Task 2.4: Wet-Lab Validation Framework with In-Silico Benchmarks

## Summary

Implemented a comprehensive wet-lab validation framework for validating that optimized sequences would actually work in the lab. The framework includes published benchmark datasets, in-silico validation runners, quantitative metrics, regression detection, and an expression predictor.

## Files Created/Modified

### New Files

1. **`src/biocompiler/wetlab_validation.py`** — Complete rewrite with:
   - `BenchmarkEntry` dataclass for curated published gene expression results
   - `BENCHMARK_DATASET` — 14 curated entries covering 10+ proteins (GFP, Insulin, HBB, EPO, TNF-alpha, IL-2, IFN-alpha, hGH, Albumin, Lysozyme) across 4 organisms (E. coli, human, CHO, yeast)
   - `run_validation_suite()` — In-silico validation runner that optimizes each benchmark protein, checks CAI, GC, restriction sites, and protein fidelity
   - `ValidationSuiteResult` with aggregate metrics: pass_rate, cai_mean/std, gc_mean/std, protein_fidelity_rate, constraint_violation_rate, comparison_vs_dnachisel
   - `check_regression()` — Regression detection that flags >5% degradation with detailed per-protein reports
   - Retained original `WetLabProtocol`, `WetLabResult`, `compare_insilico_vs_wetlab` from previous module

2. **`src/biocompiler/expression_predictor.py`** — New module with:
   - `predict_expression()` — Heuristic expression prediction based on 4 factors (CAI 40%, GC optimality 25%, mRNA stability 20%, CPB 15%)
   - `ExpressionPrediction` dataclass with predicted level, confidence, key factors, and warnings
   - `ExpressionPredictor` class for batch predictions with custom weights
   - Organism-specific GC "sweet spot" ranges for expression optimality
   - mRNA stability estimation from ATTTA motifs and GC content

3. **`tests/test_wetlab_validation.py`** — 38 tests covering:
   - Benchmark dataset integrity (minimum entries, required fields, valid proteins, CAI/GC ranges, required proteins, multiple organisms)
   - Validation suite execution (default config, custom config, metrics, DNAchisel comparison, organism filter, empty list, per-protein fields)
   - Regression detection (identical results, pass_rate/CAI regression, small changes, constraint violation, per-protein regression)

4. **`tests/test_expression_predictor.py`** — 36 tests covering:
   - ExpressionPrediction dataclass validation
   - predict_expression() with various inputs and organisms
   - GC optimality scoring
   - mRNA stability estimation
   - Confidence computation
   - ExpressionPredictor class interface
   - Integration with benchmark proteins

### Modified Files

5. **`src/biocompiler/cli.py`** — Added `validate` CLI command:
   - `biocompiler validate --suite benchmark --organism human`
   - Options: `--suite {benchmark,quick}`, `--organism`, `--strategy`, `--gc-lo`, `--gc-hi`, `--json`, `--baseline`, `--verbose`
   - JSON output for machine-readable results
   - Baseline regression detection via `--baseline FILE`

## Test Results

- 74/74 new tests PASSED
- CLI `biocompiler validate` command verified working
- Quick suite validated with `--organism ecoli --json` producing correct output

## Key Design Decisions

1. **14 benchmark entries** instead of minimum 10, to cover multiple expression systems (E. coli, HEK293, CHO, Yeast) for the same protein, enabling cross-organism validation.

2. **Four-factor heuristic model** for expression prediction: CAI (strongest predictor), GC optimality (sweet-spot based), mRNA stability (motif-based), and codon pair bias (from codon_pair_scoring module or estimated).

3. **5% degradation threshold** for regression detection — this is standard in bioinformatics software regression testing and balances sensitivity vs false positive rate.

4. **Organism-specific GC sweet spots** are narrower than the organism GC target ranges, because the optimal-for-expression range is tighter than the acceptable-for-coding range.
