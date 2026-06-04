# BioCompiler E2E Benchmark Report
## Full Head-to-Head Comparison: BioCompiler vs DNAchisel vs Naive Baseline

**Date**: 2026-06-05
**BioCompiler Version**: 9.2.0
**DNAchisel Version**: 3.2.16
**Test Genes**: 20 E. coli + 5 human therapeutic + 3 vaccine antigens + 5 stress-test + 5-organism EGFP

---

## Executive Summary

**BioCompiler is fundamentally broken as a codon optimizer.** It achieves:
- **Mean CAI: 0.6857** vs DNAchisel's **0.9788** vs Naive baseline's **1.0000**
- **0 wins out of 34** head-to-head comparisons against DNAchisel
- **38× slower** than DNAchisel (163.9ms vs 4.3ms per gene)
- **100% of genes** still have CodonOptimality violations after "optimization"
- **71% of genes** have NoCpGIsland violations (including prokaryotes where CpG is irrelevant)

The root cause is a **critical data-access bug** in `_back_translate_protein()` that causes the optimizer to start from the worst possible codons instead of the best.

---

## 1. Head-to-Head CAI Comparison

### E. coli Prokaryotic Optimization (20 genes)

| Metric | Naive | BioCompiler | DNAchisel | DNAchisel+RE |
|--------|-------|-------------|-----------|--------------|
| Mean CAI | 1.0000 | 0.6959 | 0.9728 | 0.9997 |
| CAI Range | [1.0, 1.0] | [0.50, 0.80] | [0.68, 1.0] | [1.0, 1.0] |
| Mean GC% | 0.5350 | 0.4198 | 0.5350 | 0.5287 |
| Mean Time | 0.0ms | 125.6ms | 2.4ms | 5.5ms |

### Human Eukaryotic Optimization (9 genes)

| Metric | Naive | BioCompiler | DNAchisel |
|--------|-------|-------------|-----------|
| Mean CAI | 1.0000 | 0.6461 | 0.9989 |
| CAI Range | [1.0, 1.0] | [0.56, 0.72] | [0.99, 1.0] |
| Mean GC% | 0.6444 | 0.3604 | 0.6325 |
| Mean Time | 0.0ms | 134.0ms | 2.0ms |

### Cross-Organism EGFP (5 organisms)

| Organism | Naive CAI | BioCompiler CAI | DNAchisel CAI |
|----------|-----------|-----------------|---------------|
| E. coli | 1.0000 | 0.7455 | 1.0000 |
| S. cerevisiae | 1.0000 | 0.8370 | 0.9398 |
| H. sapiens | 1.0000 | 0.6736 | 1.0000 |
| M. musculus | 1.0000 | 0.6674 | 0.9934 |
| CHO-K1 | 1.0000 | 0.6460 | ERROR |

---

## 2. Speed Comparison

| Tool | Total Time | Mean per Gene | Median | Slowdown Factor |
|------|-----------|---------------|--------|----------------|
| Naive | 0.00s | 0.0ms | 0.0ms | 1× (baseline) |
| DNAchisel | 0.16s | 4.3ms | 2.4ms | 38× faster than BC |
| DNAchisel+RE | 0.11s | 5.6ms | 2.4ms | 29× faster than BC |
| BioCompiler | 6.23s | 163.9ms | 108.7ms | 1× (slowest) |

---

## 3. Constraint Violations (BioCompiler)

| Violation | Count | % of Genes | Relevance |
|-----------|-------|-----------|-----------|
| CodonOptimality | 38/38 | 100% | **Critical**: Primary objective always fails |
| NoGTDinucleotide | 28/38 | 74% | Misguided: GT is normal in coding sequences |
| NoCpGIsland | 27/38 | 71% | **Wrong for prokaryotes**: CpG irrelevant in bacteria |
| NoRestrictionSite | 19/38 | 50% | Expected when default enzymes are applied |
| NoInstabilityMotif | 8/38 | 21% | Minor |
| GCInRange | 7/38 | 18% | **Critical**: GC constraint violated after optimization |

---

## 4. Root Cause Analysis

### Bug #1: CRITICAL — `_back_translate_protein()` uses wrong data structure

**File**: `src/biocompiler/optimization.py:2049-2063`

The function `_back_translate_protein()` accesses `SPECIES["ecoli"]` as a flat dictionary of codon→weight, but the actual structure is:
```python
SPECIES["ecoli"] = {
    "cai_weights": { "TTT": 0.867, "TTC": 1.0, ... },
    "codon_usage_validation": True
}
```

When the function does `species_cai.get(codon, 0.0)`, it gets `0.0` for **every codon** because codon keys are nested under `"cai_weights"`. The `max()` function then picks the first codon in `AA_TO_CODONS[aa]`, which is often the **worst** codon:

| Amino Acid | Buggy Pick | CAI | Correct Pick | CAI |
|-----------|------------|-----|-------------|-----|
| Leucine (L) | TTA | 0.264 | CTG | 1.000 |
| Proline (P) | CCT | 0.340 | CCG | 1.000 |
| Valine (V) | GTT | 0.686 | GTG | 1.000 |
| Alanine (A) | GCT | 0.498 | GCG | 1.000 |
| Threonine (T) | ACT | 0.540 | ACC | 1.000 |
| Glycine (G) | GGT | 0.830 | GGC | 1.000 |

**11 out of 20 amino acids** have wrong codon selection. This single bug causes the optimizer to start from CAI≈0.73 instead of CAI=1.0, and the subsequent optimization cannot recover.

### Bug #2: Organism-aware constraints not fully effective

Even though `organism_domain='prokaryote'` sets `avoid_gt=False` and `splice_low=999.0`, the BioOptimizer still:
- Reports NoCpGIsland violations (CpG avoidance active even for prokaryotes)
- Reports NoGTDinucleotide violations (GT avoidance active even for prokaryotes)
- Modifies 46.4% of codons away from the optimal starting sequence

### Bug #3: CSP Solver completely non-functional

The CSP solver (both OR-Tools and Z3 backends) returns empty sequences:
```
backend_used: NONE
cai: 0.0
sequence: ""
metadata: {'reason': 'All CSP backends unavailable or infeasible'}
```

Despite both OR-Tools and Z3 being installed and importable, the solver cannot find feasible solutions.

---

## 5. CAI Delta Analysis (BioCompiler - DNAchisel)

- **Mean delta**: -0.2910 (BioCompiler is 0.29 CAI points below DNAchisel)
- **Median delta**: -0.3092
- **Range**: [-0.4287, -0.0738]
- **Wins**: 0/34
- **Losses**: 34/34

The delta is remarkably consistent (~0.27-0.33) across almost all genes, confirming it's a systematic bug rather than edge-case behavior.

---

## 6. Impact of the Bug Fix

If `_back_translate_protein()` used `CODON_ADAPTIVENESS_TABLES` directly:

| Metric | Current | After Fix | DNAchisel |
|--------|---------|-----------|-----------|
| E. coli EGFP CAI | 0.7455 | **1.0000** | 1.0000 |
| E. coli EGFP GC% | 0.4198 | **0.4895** | 0.4895 |
| Initial sequence CAI | 0.7293 | **1.0000** | N/A |

The fix would immediately bring the starting sequence to CAI=1.0, and the optimizer would only need to make minimal changes for other constraints (restriction sites, GC range).

---

## 7. Recommendations (Priority Order)

1. **FIX `_back_translate_protein()`** — Use `CODON_ADAPTIVENESS_TABLES` or `PREFERRED_CODON_TABLES` instead of broken `SPECIES` dict. This is a one-line fix that recovers ~0.27 CAI.

2. **FIX organism-aware constraint evaluation** — `_evaluate_all_predicates()` unconditionally checks CpG/GT regardless of organism. Skip eukaryotic-only predicates for prokaryotes.

3. **FIX CSP Solver** — The OR-Tools/Z3 backends return empty sequences. Debug the model construction and feasibility.

4. **ADD CAI-aware constraint resolver** — When applying constraint fixes, minimize CAI loss per codon change. Current optimizer makes changes without considering CAI impact.

5. **Benchmark against DNAchisel WITH the same constraint set** — After fixing the above, run a fair comparison where both tools optimize for the same constraints.

---

## 8. Charts Generated

1. `chart1_ecoli_cai_comparison.png` — Bar chart: E. coli CAI by gene and tool
2. `chart2_cai_delta.png` — Bar chart: CAI delta (BC - DC) per gene
3. `chart3_speed_comparison.png` — Box plot: Runtime comparison
4. `chart4_violations.png` — Horizontal bar: Constraint violation frequencies
5. `chart5_cross_organism.png` — Bar chart: Cross-organism EGFP CAI
6. `chart6_gc_distribution.png` — Histogram: GC content distributions
7. `chart7_dashboard.png` — 4-panel summary dashboard

---

## Raw Data

Full results saved to `e2e_benchmark_results.json` (38 genes × 4 tools = 152 data points).
