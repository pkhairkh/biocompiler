# ADR-0012: CpG Avoidance in Greedy Optimizer

## Status: Accepted

## Date

2026-03-05

## Context

The NoCpGIsland predicate was failing on many optimized sequences because the
optimizer had no step to avoid CG dinucleotides. High CAI in human often requires
GC-rich codons, which naturally create CpG dinucleotides that can trigger epigenetic
silencing.

CpG dinucleotides (5'-CG-3') are significant in mammalian gene expression because
they are targets for DNA methylation. Methylation of CpG dinucleotides in promoter
regions and gene bodies is associated with transcriptional silencing through
chromatin remodeling. The NoCpGIsland predicate checks whether the designed sequence
contains a CpG island — a region of at least 200bp with GC content > 50% and an
observed/expected CpG ratio > 0.6 — which could trigger unwanted epigenetic
silencing of the transgene.

The greedy optimizer's CAI maximization step selects the highest-CAI codon at
each position, which in human often favors GC-rich codons (e.g., GCC for Alanine,
GGC for Glycine, CGC for Arginine). These high-CAI codons create CG dinucleotides
at codon boundaries (e.g., ...GCC-CGC... creates a CG across the boundary) and
within codons (CGC, CGG, CGT, CGA all contain CG). Without an explicit step to
disrupt these CpG dinucleotides, the optimizer's output frequently fails the
NoCpGIsland check.

**Alternatives Considered:**

1. **No CpG step** — Rely on the NoCpGIsland predicate to catch violations and
   report them. But: without a fix step, these violations cannot be repaired
   post-hoc without manual intervention.

2. **CpG avoidance in CAI maximization step** — Incorporate CpG avoidance into the initial
   codon selection. But: this would require a multi-objective optimization in
   the CAI maximization step, significantly complicating the greedy algorithm and potentially
   sacrificing CAI unnecessarily (many CpG positions can be fixed without CAI loss).

3. **Dedicated CpG disruption step** — Add a post-restriction-site step that
   specifically targets CG dinucleotides for disruption. This keeps the steps
   separate and allows each to be reasoned about independently.

## Decision

Add CpG dinucleotide disruption step. This step runs after cryptic splice
elimination and before final reconciliation. For each CG
dinucleotide found in the sequence, this step tries codon swaps that eliminate the
CG while preserving:

1. No reintroduction of restriction sites (already removed in restriction site removal step)
2. No creation of new cryptic splice donors above threshold (already fixed in
   cryptic splice elimination step)
3. No creation of ATTTA motifs (already removed in ATTTA motif removal step)

Post-CpG reconciliation step adds reconciliation to ensure CpG fixes don't undo restriction site
removal. This is necessary because a codon swap that eliminates a CG might
introduce a new restriction site that was previously removed.

This step processes CG dinucleotides in order of impact on the CpG island metric
(Those in the densest CG regions are prioritized). For each CG, it tries the
highest-CAI alternative codon that eliminates the CG and doesn't violate constraints
1–3 above. If no such codon exists, the CG position is left unmodified and reported
as unrepairable.

## Consequences

- **Positive**: NoCpGIsland pass rate improves significantly. The dedicated step
  ensures that CG dinucleotides are actively disrupted, rather than passively
  reported as violations.
- **Negative**: Some CAI reduction is expected — CpG-free codons may be less
  optimal than their CG-containing alternatives. For example, replacing CGC
  (Arginine, highest CAI) with AGA (Arginine, lower CAI) eliminates the CG but
  reduces CAI.
- **Negative**: This step is best-effort — not all CG positions can be eliminated
  without changing the amino acid sequence. Some amino acids (e.g., Arginine) have
  codons where CG is unavoidable (CGN codons), and the alternative AGA/AGG codons
  may be lower-CAI or may create other constraint violations.
- **Positive**: Post-CpG reconciliation step ensures that CpG fixes don't introduce
  regressions in other constraint categories, maintaining the invariant that all
  previous step fixes are preserved.

## References

- ADR-0014: Predicate Checking Delegation (NoCpGIsland in type system)
- ADR-0008: Greedy Multi-Objective Optimizer as Default (step architecture)
