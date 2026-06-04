# BioCompiler E2E Benchmark Report
## Full Head-to-Head Comparison: BioCompiler vs DNAchisel vs Naive Baseline

**Date**: 2026-06-05 (Updated after bug fixes)
**BioCompiler Version**: 9.2.0+fixes
**DNAchisel Version**: 3.2.16

---

## Executive Summary

After fixing three critical bugs (SPECIES data access, organism-aware predicates, CSP solver), BioCompiler's performance improved dramatically:

| Metric | Before Fixes | After Fixes | Change |
|--------|-------------|-------------|--------|
| **E. coli mean CAI** | 0.6959 | **0.9054** | +30.1% |
| **Human mean CAI** | 0.6461 | **0.9393** | +45.4% |
| **Overall mean CAI** | 0.6857 | **0.9142** | +33.3% |
| **Head-to-head wins** | 0/34 | **3/34** | +3 |
| **Mean runtime/gene** | 163.9ms | **76.0ms** | 2.16× faster |
| **CodonOptimality violations** | 100% | **0%** | Fixed |
| **NoCpGIsland on prokaryotes** | 71% | **0%** | Fixed |

---

## Bugs Fixed

### Bug #1: `_back_translate_protein()` SPECIES data-access bug
- **Root cause**: SPECIES dict was restructured to `{"cai_weights": {...}, "codon_usage_validation": True}`, but all consumer code accessed it as a flat `codon→float` dict, getting 0.0 for every codon weight.
- **Fix**: Added `get_species_cai_weights()` helper function; updated 9 access sites across 6 source files.
- **Impact**: Starting sequence CAI went from 0.73 → 1.0 for E. coli.

### Bug #2: Organism-aware constraint evaluation not working for prokaryotes
- **Root cause**: `_evaluate_all_predicates()` unconditionally ran eukaryote-specific checks (CpG, splice, GT) even for prokaryotes.
- **Fix**: Added `organism_name`/`organism_domain` attributes to BioOptimizer; skip eukaryotic predicates for prokaryotes; added `organism` parameter to `check_no_cryptic_splice()` and `check_no_avoidable_gt()`.
- **Impact**: Eliminated false constraint violations for prokaryotes; recovered ~0.05 CAI.

### Bug #3: CSP Solver returning empty sequences
- **Root cause**: OR-Tools engine had no diagnostic logging; Z3 engine constructor raised ImportError (caught silently); greedy fallback never attempted; SPECIES namespace mismatch.
- **Fix**: Added diagnostic logging; Z3 constructor no longer raises; greedy fallback now attempted; fixed organism key mapping; added `organism` field to types.CSPModel.
- **Impact**: All three backends (OR-Tools, Z3, Greedy) now produce valid results.

---

## Post-Fix Cross-Organism Results (EGFP, 239 aa)

| Organism | BioCompiler CAI | DNAchisel CAI | GC% | Violations |
|----------|-----------------|---------------|-----|------------|
| E. coli | 0.9144 | 1.0000 | 0.5384 | None |
| S. cerevisiae | 0.8923 | 0.9398 | 0.3445 | NoGTDinucleotide |
| H. sapiens | 0.9350 | 1.0000 | 0.5551 | None |
| M. musculus | 0.9299 | 0.9934 | 0.5356 | None |
| CHO-K1 | 0.9215 | N/A | 0.5439 | None |

---

## Remaining Gap Analysis

BioCompiler still trails DNAchisel by ~0.06-0.09 CAI on most genes. This gap is caused by:

1. **Constraint set mismatch**: BioCompiler applies more constraints (ATTTA motifs, T-runs, restriction sites by default) than DNAchisel in the benchmark
2. **Optimization algorithm**: DNAchisel uses a more efficient local search; BioCompiler's greedy optimizer sometimes makes suboptimal codon substitutions
3. **Speed**: BioCompiler is still ~18× slower than DNAchisel

### Recommendations for Next Steps

1. Add CAI-aware constraint resolver that minimizes CAI loss per constraint fix
2. Benchmark with identical constraint sets for fair comparison
3. Profile and optimize the greedy optimizer's hot path
4. Consider hybrid approach: start from DNAchisel's CAI-optimal, then apply biocompiler-specific constraints

---

## Raw Data

Full results saved to `e2e_benchmark_results.json`.
