# ADR-0018: tRNA Adaptation Index (tAI)

## Status: Accepted

## Date

2026-03-05

## Context

BioCompiler uses the Codon Adaptation Index (CAI) as its sole metric for translation efficiency. CAI, introduced by Sharp & Li (1987), measures the relative adaptiveness of codons based on their frequency in a set of highly expressed genes. While CAI is well-established and computationally simple, it has well-documented limitations as a predictor of actual translation efficiency:

1. **Reference set dependency**: CAI requires a carefully curated reference set of highly expressed genes. The quality and representativeness of this set directly affects CAI's predictive power. For organisms with limited gene expression data (e.g., CHO cells, many non-model organisms), the reference set may be small or biased, producing unreliable CAI values.

2. **No tRNA consideration**: CAI is purely based on codon usage frequencies. It does not account for tRNA gene copy numbers, tRNA expression levels, or wobble base-pairing rules. A codon may be "frequent" in highly expressed genes but correspond to a low-abundance tRNA, leading to ribosomal stalling.

3. **Wobble base-pairing ignored**: The genetic code allows non-standard (wobble) base pairing at the third position of the codon-anticodon interaction. A single tRNA species can decode multiple synonymous codons with different efficiencies. CAI treats each codon independently and cannot model this.

4. **Organism-specific tRNA pools**: Different organisms have different tRNA gene copy numbers and modification patterns. CAI captures this indirectly through codon usage, but tRNA abundance data provides a more direct and biologically grounded measurement of translation efficiency.

The tRNA Adaptation Index (tAI), introduced by dos Reis et al. (2004), addresses these limitations by directly modeling the efficiency of codon-anticodon pairing based on tRNA gene copy numbers and wobble pairing rules. tAI has been shown to correlate better with protein expression levels than CAI in multiple organisms (Tuller et al., 2010; Chu et al., 2011).

**How tAI works:**

For each codon, tAI computes a weight based on:
- The gene copy number of each tRNA that can decode the codon (via Watson-Crick or wobble pairing)
- Wobble pairing efficiency weights (e.g., G-U wobble pairs have lower efficiency than Watson-Crick pairs)
- The geometric mean of per-codon relative adaptiveness values (similar structure to CAI, but using tRNA-based weights instead of reference gene frequencies)

The wobble rules are:
- Watson-Crick (A-U, G-C): efficiency = 1.0
- G-U wobble: efficiency = 0.6 (empirically determined)
- Modified wobble (e.g., inosine-I pairing with A, C, U): efficiency varies by modification type
- Non-pairing: efficiency = 0.0

**Alternatives Considered:**

1. **Keep CAI only** — Continue using CAI as the sole translation efficiency metric. But: CAI's limitations are well-documented, and users working with non-model organisms or optimizing for maximal protein expression have requested a more biologically grounded metric.

2. **Replace CAI with tAI entirely** — Use tAI as the default optimization objective instead of CAI. But: this would be a breaking change with significant migration impact. CAI is the de facto standard in the field, and most published codon optimization work uses CAI. Dropping CAI would make it harder to compare BioCompiler results with published benchmarks.

3. **Add tAI as an alternative objective alongside CAI** — Implement tAI as a new module that users can select as an optimization objective. Keep CAI as the default for backward compatibility. Provide both metrics in optimization results so users can compare them. This follows the same pattern as the custom objectives framework (ADR-0017).

## Decision

Adopt Alternative 3: add tAI as an alternative objective alongside CAI.

The `biocompiler.tai` module provides:

1. **tAI computation**: `compute_tai(dna, organism)` — computes the tRNA Adaptation Index for a DNA sequence and organism. Uses the same geometric-mean structure as CAI, but with tRNA-based weights.

2. **tRNA gene copy number database**: Per-organism tRNA gene copy numbers derived from the Genomic tRNA Database (GtRNAdb). Supported organisms: E. coli, human, mouse, CHO, yeast. Users can provide custom tRNA copy number data for organisms not in the built-in database.

3. **Wobble pairing rules**: Configurable wobble efficiency weights based on dos Reis et al. (2004). Default weights:
   - Watson-Crick pairing: 1.0
   - G:U wobble: 0.635
   - I:A, I:C, I:U (inosine wobble): 0.0, 0.0, 0.0 (prokaryotic); 0.29, 0.29, 0.41 (eukaryotic)
   - Other modified wobble: organism-specific

4. **Per-codon tRNA efficiency weights**: `codon_trna_weights(organism)` — returns the tRNA-based weight for each codon, analogous to codon adaptiveness values but derived from tRNA data.

5. **Integration with optimization**: Users can specify `objective="tai"` or `objective=tai_objective` to optimize for tAI instead of CAI. The result includes both `cai` and `tai` fields when tAI is computed.

6. **Result field**: `OptimizationResult.tai` — the tAI value of the optimized sequence (populated when tAI is computed, None otherwise).

**Default behavior**: CAI remains the default optimization objective. tAI is opt-in. This ensures full backward compatibility while providing the more biologically grounded metric for users who need it.

**When to use tAI instead of CAI:**
- Optimizing for maximal protein expression in a well-characterized organism with reliable tRNA gene copy number data
- Working with non-model organisms where the CAI reference set is small or unreliable, but tRNA gene copy numbers are available from genome sequencing
- Detecting ribosomal stalling risk: codons with high CAI but low tAI may indicate tRNA bottlenecks
- Validating CAI results: if tAI and CAI disagree significantly, the sequence may have codon-tRNA mismatches that warrant investigation

## Consequences

- **Positive**: Users now have access to a more biologically meaningful translation efficiency metric that accounts for tRNA availability and wobble pairing. This is particularly valuable for organisms where CAI reference sets are unreliable.
- **Positive**: The dual-metric approach (CAI + tAI) allows users to detect discrepancies between codon usage-based and tRNA-based optimization, flagging potential ribosomal stalling sites.
- **Positive**: tAI can be used as a custom objective via the existing `objectives` framework, requiring no special optimizer changes.
- **Positive**: Backward compatible — CAI remains the default, and existing code continues to work unchanged.
- **Negative**: tAI requires per-organism tRNA gene copy number data. For organisms not in the built-in database, users must provide this data manually. The quality of tAI depends on the accuracy and completeness of the tRNA gene copy number data.
- **Negative**: Wobble pairing efficiency weights are empirically determined and may not be accurate for all organisms or cell types. The dos Reis et al. (2004) weights were calibrated on a limited set of organisms.
- **Negative**: Adding tAI increases the module surface area and the number of metrics users need to understand. Documentation must clearly explain when to use tAI vs. CAI.
- **Trade-off**: tAI is more biologically grounded but less widely used than CAI. Published benchmarks and comparisons predominantly use CAI, making it harder to validate tAI results against the literature. We accept this trade-off because tAI provides complementary information that CAI cannot capture.

## References

- Sharp, P.M., & Li, W.H. (1987). The codon Adaptation Index - a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281-1295.
- dos Reis, M., Savva, R., & Wernisch, L. (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036-5044.
- Tuller, T., Carmi, A., Vestsigian, K., et al. (2010). An evolutionarily conserved mechanism for controlling the efficiency of protein translation. *Cell*, 141(2), 344-354.
- Chu, D., Kazana, E., Bellanger, N., et al. (2011). Translation elongation can control translation initiation on eukaryotic mRNAs. *EMBO Journal*, 30(1), 115-126.
- Genomic tRNA Database (GtRNAdb): http://gtrnadb.ucsc.edu/
- ADR-0017: Feature Parity with DNAchisel (custom objectives framework)
- ADR-0008: Greedy Multi-Objective Optimizer as Default
