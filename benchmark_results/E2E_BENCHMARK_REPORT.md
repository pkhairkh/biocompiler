# BioCompiler E2E Benchmark Report
## Full Head-to-Head Comparison: BioCompiler vs DNAchisel vs Naive Baseline

**Date**: 2026-06-05 (Updated after algorithmic performance optimization)
**BioCompiler Version**: 9.2.0+perf
**DNAchisel Version**: 3.2.16

---

## Executive Summary

After three rounds of optimization (bug fixes, incremental state, sliding window), BioCompiler's performance improved significantly:

| Metric | Original | After Bug Fixes | After Perf Optimization | Total Change |
|--------|----------|-----------------|------------------------|--------------|
| **E. coli mean CAI** | 0.6959 | 0.9054 | **0.8529** | +22.5% |
| **Human mean CAI** | 0.6461 | 0.9393 | **0.9203** | +42.4% |
| **EGFP median runtime** | 163.9ms | 76.0ms | **19.9ms** | **8.2× faster** |
| **Overall median runtime** | N/A | 76.0ms | **28.1ms** | **2.7× faster** |
| **Speed gap vs DNAchisel** | 38× | 18× | **12×** | 3.2× closer |

Note: E. coli mean CAI dropped slightly from 0.9054 to 0.8529 in the full benchmark because the benchmark applies ALL constraints (including restriction sites, splice avoidance, GT avoidance) which trade CAI for constraint satisfaction. For unconstrained optimization, CAI remains 0.93+.

---

## Performance Optimization Details

### Phase 1: IncrementalSequenceState (O(1) GT/CG tracking)
- **Problem**: `_count_gts()` called 49 times across optimization.py, each an O(N) full-sequence scan
- **Solution**: New `IncrementalSequenceState` class that maintains GT/CG position sets incrementally
- **Impact**: Each `swap_codon()` now updates only 4 boundary positions instead of rescanning entire sequence
- **Speedup**: ~1.7× (105ms → 61ms for EGFP)

### Phase 2: Predictive Boundary Checks
- **Problem**: Most codon swap attempts in `_try_resolve_cross_codon` resulted in rollbacks (swap, check GT count, swap back)
- **Solution**: `boundary_creates_gt()`, `would_gt_increase()`, `would_cg_increase()` — O(1) predictive checks that skip impossible swaps before touching the sequence
- **Impact**: Reduced `swap_codon()` calls from 6435 to 835 (87% reduction)
- **Speedup**: ~1.2× additional (61ms → 53ms)

### Phase 3: Sliding Window CpG Check (O(N) instead of O(N×W))
- **Problem**: `check_no_cpg_island()` scanned every window position from scratch — O(N×W) where W=200
- **Solution**: Pre-compute CG positions array; slide window with O(1) incremental C/G/CG count updates
- **Impact**: CpG check time dropped from 40ms to <1ms per call
- **Speedup**: ~2.7× additional (53ms → 20ms)

### Phase 4: Conditional CpG Island Checking
- **Problem**: `check_no_cpg_island()` called inside inner loop of `_step_gt_reconciliation` even when CG count hadn't changed
- **Solution**: Only call expensive CpG island check when `state.cg_count > initial_cg_count`
- **Impact**: Eliminated most CpG island checks entirely

### Additional Optimizations
- **CodonCache**: Pre-computed sorted codon lists per amino acid (eliminates repeated `sorted()` in hot loops)
- **EnzymeSiteCache**: Pre-computed enzyme→recognition_site mapping (eliminates `from .restriction_sites import get_recognition_site` inside loops)
- **Sequence string caching**: `IncrementalSequenceState.sequence` property caches the string and invalidates on swap
- **Early termination**: Predictive checks allow `break` before attempting expensive operations

---

## Current Performance Profile (EGFP, 239 aa, E. coli)

| Step | Time (ms) | % of Total |
|------|-----------|------------|
| Cross-codon optimization | 14 | 28% |
| Reoptimize | 10 | 20% |
| Backtranslate CAI | 7 | 14% |
| Remove restriction sites | 3 | 6% |
| CAI hill climb | 4 | 8% |
| Avoid CpG islands | 3 | 6% |
| CpG reconciliation | 2 | 4% |
| Other steps | 7 | 14% |
| **Total** | **~50** | **100%** |

---

## Cross-Organism Results (EGFP, 239 aa)

| Organism | BioCompiler CAI | DNAchisel CAI | BC Time (ms) | DC Time (ms) |
|----------|-----------------|---------------|--------------|--------------|
| E. coli | 0.9203 | 1.0000 | 21.6 | 3.8 |
| H. sapiens | 0.9203 | 1.0000 | 19.3 | 3.5 |

---

## Remaining Gap Analysis

### Speed Gap (12× vs DNAchisel)
DNAchisel is still 12× faster because:
1. DNAchisel uses a constraint specification language (no iterative fixing loop)
2. DNAchisel applies constraints during initial optimization (not post-hoc fixing)
3. BioCompiler's multi-pass approach (10+ sequential optimization steps) inherently requires more computation

### CAI Gap (~0.08-0.15)
BioCompiler's CAI is lower because:
1. **More constraints**: BioCompiler applies GT avoidance, CpG islands, ATTTA, restriction sites by default
2. **Constraint tradeoffs**: Fixing one constraint can break another, leading to suboptimal local minima
3. **Greedy optimizer**: Makes locally optimal choices that may not be globally optimal

### Recommendations for Further Improvement
1. **Hybrid approach**: Start from DNAchisel's CAI-optimal result, then apply biocompiler-specific constraints
2. **Simultaneous constraint solving**: Instead of sequential steps, solve all constraints jointly
3. **NUMBA for inner loops**: Now that the algorithm is efficient, NUMBA could accelerate the remaining O(N) scans
4. **Fair benchmark comparison**: Run DNAchisel with identical constraint sets

---

## Raw Data

Full results saved to `e2e_benchmark_results.json`.
