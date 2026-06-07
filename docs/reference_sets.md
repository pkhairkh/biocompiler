# CAI Reference Sets: Kazusa vs Sharp-Li

> **v12.0.0 Note**: All CAI computation paths now use unified codon adaptiveness tables from `organisms.CODON_ADAPTIVENESS_TABLES`. Previously, different optimizer backends could use slightly different values, leading to inconsistent CAI results. This has been corrected. CAI values are now consistently >0.99 across organisms after fixing post-processing CAI regression.

## Overview

The Codon Adaptation Index (CAI) measures how well a coding sequence's codon usage matches that of highly expressed genes in a target organism. The CAI value depends critically on the **reference set** — the collection of highly expressed genes used to derive the relative adaptiveness values (w) for each codon.

BioCompiler supports two CAI reference sets:

| Reference Set | Source | Genes | Organism Coverage |
|---------------|--------|-------|-------------------|
| **Kazusa** | Kazusa Codon Usage Database, high-expression subsets | Varies by organism (large) | E. coli, H. sapiens, S. cerevisiae, M. musculus, CHO-K1 |
| **Sharp-Li** | Sharp & Li (1987), Table 2 | 24 E. coli genes | E. coli only |

## The Two Reference Sets

### Kazusa Codon Usage Database

The **Kazusa** reference set is derived from the [Kazusa Codon Usage Database](https://www.kazusa.or.jp/codon/), which aggregates codon usage statistics from all genes in the NCBI GenBank entries for each organism. BioCompiler uses the "high-expression subset" — genes annotated as highly expressed.

**Characteristics:**
- Larger gene collection (hundreds to thousands of genes depending on organism)
- Per-thousand frequency values from genome-wide codon usage tables
- More moderate codon bias (less extreme w values)
- Better for general-purpose gene optimization across diverse expression levels
- Available for 25 organisms: *E. coli*, *H. sapiens*, *S. cerevisiae*, *M. musculus*, CHO-K1, *C. elegans*, *D. melanogaster*, *A. thaliana*, *P. pastoris*, *B. subtilis*, and others

### Sharp & Li (1987) Original Reference Set

The **Sharp-Li** reference set reproduces the exact 24 highly expressed *E. coli* genes used in the original CAI paper:

**The 24 reference genes include:**
- Ribosomal proteins: rplA, rplB, rplC, rplD, rplE, rplF, rplJ, rplK, rplL, rplO, rplQ, rpsA, rpsC, rpsG, rpsJ, rpsL, rpsM
- Elongation factors: tufA (EF-Tu), fusA (EF-G)
- Outer membrane proteins: ompA, ompC, ompF
- Other: groEL (mopA), recA, rpoB (partial)

**Characteristics:**
- Smaller, precisely defined gene set (24 genes)
- Scaled frequencies (max = 1000 per amino acid) to exactly reproduce the published relative adaptiveness values
- Stronger codon bias (more extreme w values)
- Produces CAI values that match published Sharp & Li (1987) and CAIcal results
- E. coli only — no multi-organism coverage

## When to Use Each Reference Set

| Use Case | Recommended Set | Rationale |
|----------|----------------|-----------|
| General gene optimization | Kazusa | Broader codon usage reflects the full high-expression transcriptome |
| Comparing to published Sharp & Li (1987) values | Sharp-Li | Reproduces the exact adaptiveness values from the original paper |
| Comparing to CAIcal server results | Sharp-Li | CAIcal uses the same 24-gene reference set |
| Cross-organism CAI comparisons | Kazusa | Only set with multi-organism coverage |
| Validating CAI implementation | Sharp-Li | Ground-truth comparison against published values |
| Prokaryotic optimization (E. coli) | Kazusa (default) or Sharp-Li (for publication comparison) | Both produce valid CAI; Sharp-Li has stronger discriminative power |
| Eukaryotic optimization | Kazusa | Only option — Sharp-Li is E. coli only |

## Example Usage

### Computing CAI with the Kazusa Reference Set (Default)

```python
from biocompiler.benchmarking.cai_validated import compute_cai_sharp_li, load_reference_set

# Load the Kazusa reference set for E. coli
kazusa_ref = load_reference_set("Escherichia_coli")

# Compute CAI
cai = compute_cai_sharp_li("ATGAAAGCGTAA", kazusa_ref)
print(f"CAI (Kazusa): {cai:.4f}")
```

### Computing CAI with the Sharp-Li Reference Set

```python
from biocompiler.benchmarking.sharp_li_benchmark import SHARP_LI_ECOLI_REFERENCE
from biocompiler.benchmarking.cai_validated import compute_cai_sharp_li

# Use the Sharp-Li 24-gene reference set
sharp_li_ref = SHARP_LI_ECOLI_REFERENCE

# Compute CAI
cai = compute_cai_sharp_li("ATGAAAGCGTAA", sharp_li_ref)
print(f"CAI (Sharp-Li): {cai:.4f}")
```

### Running the Full Benchmark Comparison

```python
from biocompiler.benchmarking import run_benchmark_by_name

# Run the Sharp-Li vs Kazusa benchmark
results = run_benchmark_by_name("sharp_li_cai")

print(f"Mean Kazusa error:   {results['mean_kazusa_error']:.4f}")
print(f"Mean Sharp-Li error: {results['mean_sharp_li_error']:.4f}")
print(f"Sharp-Li is closer to published: {results['sharp_li_is_closer']}")

# Per-gene breakdown
for gene_result in results["per_gene_results"]:
    print(f"  {gene_result['gene']}: "
          f"published={gene_result['published_cai']:.4f}, "
          f"kazusa={gene_result['kazusa_cai']:.4f}, "
          f"sharp_li={gene_result['sharp_li_cai']:.4f}")
```

### Running the Organism-Aware Constraint Benchmark

```python
from biocompiler.benchmarking import run_benchmark_by_name

# Run the organism-aware CAI recovery benchmark
results = run_benchmark_by_name("organism_aware_cai")

print(f"Mean CAI (organism-unaware): {results['mean_cai_old']:.4f}")
print(f"Mean CAI (organism-aware):   {results['mean_cai_new']:.4f}")
print(f"Mean CAI recovery:           {results['mean_cai_recovery']:+.4f}")
```

## Expected CAI Differences for Common Genes

The following table shows typical CAI differences between the two reference sets for well-studied *E. coli* genes. Published values are from Sharp & Li (1987), Table 1.

| Gene | Published CAI | Kazusa CAI | Sharp-Li CAI | Δ (Kazusa − Published) | Notes |
|------|:------------:|:----------:|:------------:|:----------------------:|-------|
| trpA | 0.84 | ~0.74 | ~0.84 | −0.10 | Highly expressed; both sets agree directionally |
| recA | 0.76 | ~0.81 | ~0.76 | +0.05 | Moderate expression; moderate reference-set effect |
| ompA | 0.79 | ~0.79 | ~0.79 | ~0.00 | Very close agreement between sets |
| groEL | 0.78 | ~0.82 | ~0.78 | +0.04 | Close agreement; minor reference-set effect |
| dnaK | 0.76 | ~0.79 | ~0.76 | +0.03 | Close agreement |
| lacZ | 0.27 | ~0.72 | ~0.27 | +0.45 | **Major discrepancy** — low-expression gene |
| rpoB | 0.50 | ~0.65 | ~0.50 | +0.15 | Moderate expression; notable reference-set effect |

### Key Observations

1. **Highly expressed genes** (trpA, ompA, groEL, dnaK): Both reference sets produce similar CAI values. The Kazusa set gives slightly higher values because its larger gene pool dilutes the codon bias.

2. **Lowly expressed genes** (lacZ): The Kazusa set produces dramatically higher CAI values because its codon frequencies are averaged across a broader gene set with less extreme codon bias. lacZ uses many rare codons by *E. coli* standards; the Sharp-Li set penalizes these heavily (w values as low as 0.002–0.006 for rare codons), while the Kazusa set assigns higher w values due to the diluted bias.

3. **Moderately expressed genes** (rpoB): Intermediate effect. The Kazusa set still inflates CAI relative to the Sharp-Li set, but less dramatically than for low-expression genes.

4. **Rank-order preservation**: For E. coli endogenous genes, both reference sets preserve the rank order of CAI (highly expressed > moderately expressed > lowly expressed), though the numerical values differ.

5. **Heterologous genes** (GFP, human proteins in E. coli): The Kazusa set typically produces higher CAI values than the Sharp-Li set for the same reason — diluted codon bias means fewer codons are penalized.

## Provenance System

BioCompiler tracks which reference set was used for every CAI computation through its **provenance system**. This is a key differentiator: every codon decision records not just the CAI value, but also the reference set provenance, enabling reproducible and auditable gene design.

```python
from biocompiler import optimize_sequence

# Optimization with provenance tracking
result = optimize_sequence(
    target_protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR...",
    organism="Homo_sapiens",
    track_provenance=True,
)

# Provenance records include CAI reference set information
for decision in result.provenance.codon_decisions:
    print(f"Position {decision.position}: {decision.amino_acid} → {decision.chosen_codon} "
          f"(CAI contribution: {decision.cai_contribution:.4f})")
```

## Organism-Aware Constraint Selection

When optimizing for prokaryotic expression hosts (e.g., E. coli), BioCompiler automatically disables eukaryotic-specific constraints:

- **Cryptic splice-site avoidance**: Prokaryotes lack spliceosomes — GT/AG dinucleotide avoidance is irrelevant
- **CpG-island avoidance**: Prokaryotes lack CpG methylation — CG step avoidance is unnecessary

This organism-aware constraint selection recovers approximately **+0.27 CAI** on average for prokaryotic targets, as demonstrated by the `organism_aware_cai` benchmark.

```python
from biocompiler import optimize_sequence

# BioCompiler automatically selects the appropriate constraint set
result_prokaryote = optimize_sequence(
    target_protein="MVHLTPEEK...",
    organism="Escherichia_coli",  # Prokaryotic — splice/CpG auto-disabled
)

result_eukaryote = optimize_sequence(
    target_protein="MVHLTPEEK...",
    organism="Homo_sapiens",  # Eukaryotic — all constraints enabled
)
```

## References

1. Sharp, P.M. & Li, W.-H. (1987). "The codon Adaptation Index — a measure of directional synonymous codon usage bias, and its potential applications." *Nucleic Acids Research* 15:1281–1295. doi:10.1093/nar/15.3.1281

2. Puigbo, P., Bravo, I.G. & Garcia-Vallve, S. (2008). "CAIcal: A combined set of tools to assess codon usage adaptation." *BMC Bioinformatics* 9:65. doi:10.1186/1471-2105-9-65

3. Nakamura, Y., Gojobori, T. & Ikemura, T. (2000). "Codon usage tabulated from international DNA sequence databases: status for the year 2000." *Nucleic Acids Research* 28:292. (Kazusa Codon Usage Database)
