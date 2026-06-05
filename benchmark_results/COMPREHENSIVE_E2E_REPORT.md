# Comprehensive E2E Benchmark: BioCompiler vs DNAchisel

## Fairness Guarantee

Both tools are evaluated with `compute_cai_validated` — the same validated
CAI function following Sharp & Li (1987). DNAchisel's own CAI output is
**NOT trusted**. DNAchisel uses `CodonOptimize(species=...)` as an objective,
but CAI is always recomputed with our validated evaluator.

## Per-Gene Results

| Gene | Organism | BC CAI | DC CAI | BC Time (ms) | DC Time (ms) | Speed Ratio | BC Wins |
|------|----------|--------|--------|-------------|-------------|-------------|---------|
| ADH1_Yeast | Yeast | 1.0000 | 0.9871 | 25.8 | 10.5 | 2.45x | ✓ |
| Albumin | Human | 1.0000 | 1.0000 | 51.2 | 5.2 | 9.80x | tie |
| Erythropoietin | Human | 1.0000 | 1.0000 | 12.5 | 1.8 | 6.78x | tie |
| GFP | E. coli | 0.9989 | 0.9559 | 1.5 | 22.5 | 0.07x | ✓ |
| GFP_Human | Human | 1.0000 | 1.0000 | 20.9 | 2.3 | 9.13x | tie |
| GFP_Yeast | Yeast | 0.9108 | 0.9393 | 114.3 | 6.6 | 17.32x |  |
| Growth_Hormone | Human | 1.0000 | 1.0000 | 10.6 | 1.9 | 5.55x | tie |
| HBB_Ecoli | E. coli | 0.9888 | 0.9650 | 1.1 | 5.2 | 0.21x | ✓ |
| HBB_beta_globin | Human | 1.0000 | 1.0000 | 14.7 | 1.7 | 8.88x | tie |
| Insulin_A_chain | Human | 1.0000 | 1.0000 | 1.1 | 1.0 | 1.13x | tie |
| Insulin_B_chain | Human | 1.0000 | 1.0000 | 1.3 | 0.8 | 1.58x | tie |
| PGK1_Yeast | Yeast | 0.9968 | 0.9247 | 16.8 | 6.8 | 2.47x | ✓ |
| T4_lysozyme | E. coli | 1.0000 | 0.9699 | 1.0 | 8.9 | 0.11x | ✓ |
| TDH3_Yeast | Yeast | 0.9242 | 0.9765 | 11.7 | 6.5 | 1.79x |  |
| Taq_polymerase | E. coli | 0.9990 | 0.9612 | 5.9 | 89.5 | 0.07x | ✓ |
| mCherry | E. coli | 0.9916 | 0.9547 | 1.5 | 11.1 | 0.13x | ✓ |
| p53_DBD | Human | 1.0000 | 1.0000 | 14.8 | 2.0 | 7.30x | tie |

## Aggregate Statistics

- **Mean CAI across all genes**: BioCompiler = 0.9888, DNAchisel = 0.9785
- **Mean speed ratio**: 1.66x (DC faster)
- **Head-to-head wins**: BC = 7, DC = 2, Ties = 8
- **Paired t-test p-value**: 0.157013 (not significant)

## Per-Organism Breakdown

### E. coli

- BC mean CAI: 0.9957
- DC mean CAI: 0.9613
- BC mean time: 2.2 ms
- DC mean time: 27.5 ms
- BC wins: 5, DC wins: 0

### Human

- BC mean CAI: 1.0000
- DC mean CAI: 1.0000
- BC mean time: 15.9 ms
- DC mean time: 2.1 ms
- BC wins: 0, DC wins: 0

### Yeast

- BC mean CAI: 0.9579
- DC mean CAI: 0.9569
- BC mean time: 42.1 ms
- DC mean time: 7.6 ms
- BC wins: 2, DC wins: 2
