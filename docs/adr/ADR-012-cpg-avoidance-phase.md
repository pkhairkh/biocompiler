# ADR-012: CpG Avoidance Phase in Greedy Optimizer

## Status: Accepted

## Date

2026-03-05

## Context

The NoCpGIsland predicate was failing on many optimized sequences because the
optimizer had no phase to avoid CG dinucleotides. High CAI in human often requires
GC-rich codons, which naturally create CpG dinucleotides that can trigger epigenetic
silencing.

CpG dinucleotides (5'-CG-3') are significant in mammalian gene expression because
they are targets for DNA methylation. Methylation of CpG dinucleotides in promoter
regions and gene bodies is associated with transcriptional silencing through
chromatin remodeling. The NoCpGIsland predicate checks whether the designed sequence
contains a CpG island — a region of at least 200bp with GC content > 50% and an
observed/expected CpG ratio > 0.6 — which could trigger unwanted epigenetic
silencing of the transgene.

The greedy optimizer's Phase 1 (CAI maximization) selects the highest-CAI codon at
each position, which in human often favors GC-rich codons (e.g., GCC for Alanine,
GGC for Glycine, CGC for Arginine). These high-CAI codons create CG dinucleotides
at codon boundaries (e.g., ...GCC-CGC... creates a CG across the boundary) and
within codons (CGC, CGG, CGT, CGA all contain CG). Without an explicit phase to
disrupt these CpG dinucleotides, the optimizer's output frequently fails the
NoCpGIsland check.

**Alternatives Considered:**

1. **No CpG phase** — Rely on the NoCpGIsland predicate to catch violations and
   report them. But: without a fix phase, these violations cannot be repaired
   post-hoc without manual intervention.

2. **CpG avoidance in CAI phase** — Incorporate CpG avoidance into the initial
   codon selection. But: this would require a multi-objective optimization in
   Phase 1, significantly complicating the greedy algorithm and potentially
   sacrificing CAI unnecessarily (many CpG positions can be fixed without CAI loss).

3. **Dedicated CpG disruption phase** — Add a post-restriction-site phase that
   specifically targets CG dinucleotides for disruption. This keeps the phases
   separate and allows each to be reasoned about independently.

## Decision

Add Phase 7.5: CpG dinucleotide disruption. This phase runs after Phase 7 (cryptic
splice elimination) and before Phase 8 (final reconciliation). For each CG
dinucleotide found in the sequence, the phase tries codon swaps that eliminate the
CG while preserving:

1. No reintroduction of restriction sites (already removed in Phase 5)
2. No creation of new cryptic splice donors above threshold (already fixed in
   Phase 7)
3. No creation of ATTTA motifs (already removed in Phase 6)

Phase 8.5 adds reconciliation to ensure CpG fixes don't undo restriction site
removal. This is necessary because a codon swap that eliminates a CG might
introduce a new restriction site that was previously removed.

The phase processes CG dinucleotides in order of impact on the CpG island metric
(Those in the densest CG regions are prioritized). For each CG, it tries the
highest-CAI alternative codon that eliminates the CG and doesn't violate constraints
1–3 above. If no such codon exists, the CG position is left unmodified and reported
as unrepairable.

## Consequences

- **Positive**: NoCpGIsland pass rate improves significantly. The dedicated phase
  ensures that CG dinucleotides are actively disrupted, rather than passively
  reported as violations.
- **Negative**: Some CAI reduction is expected — CpG-free codons may be less
  optimal than their CG-containing alternatives. For example, replacing CGC
  (Arginine, highest CAI) with AGA (Arginine, lower CAI) eliminates the CG but
  reduces CAI.
- **Negative**: The phase is best-effort — not all CG positions can be eliminated
  without changing the amino acid sequence. Some amino acids (e.g., Arginine) have
  codons where CG is unavoidable (CGN codons), and the alternative AGA/AGG codons
  may be lower-CAI or may create other constraint violations.
- **Positive**: Phase 8.5 reconciliation ensures that CpG fixes don't introduce
  regressions in other constraint categories, maintaining the invariant that all
  previous phase fixes are preserved.

## References

- ADR-0011: Predicate Checking Delegation (NoCpGIsland in type system)
- ADR-0008: Greedy Multi-Objective Optimizer as Default (phase architecture)
