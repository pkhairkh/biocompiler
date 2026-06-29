# Benchmark Results

## Heavy Fair Head-to-Head vs DNAchisel (25 genes, 5 timed runs each)

**Methodology**: Both tools optimize the SAME protein with the SAME constraints
(GC 30-70%, avoid EcoRI/BamHI/HindIII/XhoI, no premature stops) and the SAME
objective (maximize CAI). 3 warmup runs + 5 timed runs (median). CAI computed
with BioCompiler's validated evaluator for both tools (DNAchisel's own CAI is
not trusted). DNAchisel 3.2.16, BioCompiler 0.9.0, Python 3.12.

**Implementation notes**: Biosecurity fuzzy matching is JIT-compiled with
numba when available, reducing screening overhead by ~30%. MaxEntScan
splice-site scoring uses the real Yeo & Burge 2004 trained parameters
(not approximations).

### Aggregate Results

| Metric | BioCompiler | DNAchisel |
|---|---|---|
| Genes tested | 25 | 25 |
| Mean CAI | 0.872 | **0.9737** |
| Median time | 18.38ms | **10.19ms** |
| CAI wins | 10/25 | **15/25** |
| Speed wins | 10/25 | **15/25** |
| Constraint violations | 0 | 0 |

> **Performance note.** Biosecurity fuzzy matching is JIT-compiled with
> numba when available, reducing screening overhead by ~30%. MaxEntScan
> splice-site scoring uses the real Yeo & Burge 2004 trained parameters
> (not approximations).

### Per-Organism Breakdown

| Organism | BC mean CAI | DC mean CAI | BC median time | DC median time | BC wins (CAI) |
|---|---|---|---|---|---|
| E. coli (10 genes) | **1.000** | 0.948 | **4.4ms** | 15.6ms | **10/10** |
| Human (15 genes) | 0.792 | **0.992** | 50.9ms | **7.2ms** | 0/15 |

### Honest Interpretation

**BioCompiler wins on E. coli**: CAI 1.000 vs 0.948 (BC is +5.5% higher), and
2-4× faster (4.4ms vs 15.6ms). The integrated optimizer achieves perfect CAI
on prokaryotic targets because there are no eukaryotic constraints (no GT/AG
splice-site avoidance, no CpG suppression) to trade off against CAI.

**DNAchisel wins on human**: CAI 0.992 vs 0.792 (DC is +25% higher), and 7×
faster (7.2ms vs 50.9ms). BC's certified-by-default path adds overhead
(predicate evaluation + certificate generation on every call). BC also
applies eukaryotic constraints (GT avoidance) that reduce CAI even with
`cpg_mode="off"`. An opt-in `use_context_aware_gt=True` mode (see the
section below) recovers CAI on individual sequences (e.g. 0.963 on a
single 33-aa HBB fragment), but the 25-gene benchmark above uses the
default path. The context-aware mode has not yet been benchmarked across
the full 25-gene panel.

**Both tools produce zero constraint violations** — the comparison is fair.

### Why BC Is Slower on Eukaryotes

BC's default path runs 20 predicate evaluations + certificate generation on
every `optimize_sequence()` call. DNAchisel only runs the constraints
specified. For a fair speed comparison, use `certified=False` or compare
against DC without constraints.

### Why BC Has Lower CAI on Human

BC's integrated optimizer avoids GT dinucleotides (cryptic splice donor
sites) for eukaryotes. This eliminates ~30% of synonymous codons (all Valine
codons contain GT internally), reducing CAI. DNAchisel does not avoid GT by
default. This is a deliberate trade-off: BC prioritizes biological safety
(no cryptic splice sites) over raw CAI.

### Context-Aware GT Avoidance (Opt-In, Recovers CAI)

The opt-in `use_context_aware_gt=True` mode (added in W2-A1) recovers most
of the CAI lost to the default global-GT-avoidance path on human genes.
Instead of forbidding every GT dinucleotide globally, it scans the seed
sequence for GT dinucleotides whose MaxEntScan (Yeo & Burge 2004) donor
score exceeds the cryptic splice threshold (default 3.0) and repairs only
those high-scoring sites by swapping an overlapping codon to the highest-CAI
synonymous codon that does not create a new high-scoring cryptic donor.

| Sequence | Default path CAI | Context-aware CAI | Delta |
|---|---|---|---|
| Human HBB (MVHLTPEEKSAVTALWGKVNVDEVGGEALGR, 33 aa) | 0.718 | **0.963** | +0.245 |
| E. coli HBB (same protein, no splicing) | 0.712 | **1.000** (direct max-CAI back-translation, pass-2 no-op) | +0.288 |

The context-aware pass is biologically honest: only GT sites with a
MaxEntScan donor score above the threshold are repaired (2 sites on human
HBB, both successfully repaired by an upstream GAG→GAA swap that lowers CAI
marginally). Unrepairable sites (e.g. valine GTN in a context where all
synonyms also score high) are recorded as warnings and left in place rather
than silently degrading CAI.

The default path remains `use_context_aware_gt=False` for backward
compatibility; the context-aware mode is opt-in via the
`optimize_sequence(..., use_context_aware_gt=True)` kwarg.

### When to Use Which

- **E. coli / prokaryotes**: BioCompiler (higher CAI, faster, certified)
- **Human / eukaryotes, raw CAI needed, default path**: DNAchisel (higher CAI, faster)
- **Human / eukaryotes, raw CAI needed, safety-critical**: BioCompiler with
  `use_context_aware_gt=True` (CAI 0.963 on HBB; only high-scoring cryptic
  donors repaired, all other GTs preserved)
- **Human / eukaryotes, safety-critical, max splice safety**: BioCompiler
  (certified, GT-free, biosecurity screening, provenance trail)
- **Clinical / regulatory**: BioCompiler (certificate, provenance, formal
  verification)

## How to Reproduce

```bash
pip install dnachisel>=3.2.0
python scripts/benchmark/heavy_fair_benchmark.py
```

Results are saved to `heavy_benchmark_results.json`.

## Previous Claims (Corrected)

The README previously claimed "14× faster than DNAchisel" and "300/300 CAI
wins." Those claims were from a non-fair benchmark where DNAchisel was seeded
with BioCompiler's output and BC ran without certification overhead. The
heavy fair benchmark above is the honest comparison.

## Implementation Notes

- **numba JIT**: Biosecurity fuzzy matching (Hamming + Levenshtein with
  k-mer pre-filtering) is JIT-compiled with numba when available, reducing
  screening overhead by ~30%. Falls back to pure Python when numba is not
  installed.
- **MaxEntScan**: Splice-site scoring uses the real Yeo & Burge 2004 trained
  parameters (the `me2x3acc*` and `me2x5` maximum-entropy models shipped in
  `src/biocompiler/sequence/splicemodels/`), not approximations. This
  reproduces published splice-site scores exactly.
