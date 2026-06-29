# ADR-0003: Non-Deterministic Finite-State Transducers for Splicing

> **Note**: BDDs were evaluated but NOT implemented. Python sets are used for isoform set representation instead.

## Status

Accepted

## Date

2026-05-30

## Context

The splicing stage must model alternative splicing — the process by which a single pre-mRNA sequence can produce multiple distinct mRNA isoforms through different combinations of exon inclusion and exclusion. The model must capture the essential non-determinism of the process (multiple valid isoforms for the same input) while maintaining deterministic computation (the same input always produces the same isoform set). The model must also support cellular context parameterization (different cell types produce different isoform sets due to different splicing factor concentrations).

**Alternatives Considered:**

1. **Probabilistic Context-Free Grammar (PCFG)** — Assigns probabilities to alternative splicing outcomes based on training data. Well-understood formalism; can rank isoforms by likelihood. But: requires training data and parameter estimation; probabilities need calibration; composition requires independence assumptions; no formal guarantees about correctness.

2. **Hidden Markov Model (HMM)** — Well-established for gene finding and sequence analysis. Can model position-specific emission probabilities. But: produces a single most-likely output (Viterbi decoding) rather than a set of possible outputs; does not naturally model the transduction aspect (input→output transformation); probabilistic outputs are not compositional without independence assumptions.

3. **Neural network (e.g., SpliceAI-style)** — State-of-the-art splice site prediction accuracy. But: black-box model with no interpretable rules; no formal guarantees; requires large training datasets; non-deterministic across runs without fixed seeds; cannot be composed with the deterministic type system.

4. **Non-Deterministic Finite-State Transducer (NDFST)** — Models splicing as a set-valued function: given an input, produces the set of all possible outputs. Computation is deterministic (same set for same input). No probabilities. Composes via set union. Cellular context parameterizes transitions.

## Decision

Model splicing as a Non-Deterministic Finite-State Transducer (Alternative 4). The NDFST captures the essential non-determinism: alternative splicing produces a set of possible isoforms, not a probability distribution over isoforms. It avoids the need for probability estimates: the system computes the set of possible isoforms and checks whether all of them satisfy the desired properties. It composes: the combined isoform set for a circuit is the Cartesian product of individual isoform sets, with no independence assumptions. The cellular context parameterizes the NDFST by enabling or disabling transitions based on splicing factor concentration thresholds.

## Consequences

- Positive: (1) Formal guarantees: the isoform set is complete (all valid parse paths explored) and sound (every isoform satisfies the grammar). (2) No training data required: the grammar is constructed from curated biological knowledge. (3) Deterministic computation: same input always produces same isoform set, enabling reproducibility and certificate verification. (4) Compositional: NDFSTs compose via set operations without probability. (5) Cellular context parameterization is natural (enable/disable transitions).
- Negative: (1) Cannot rank isoforms by likelihood: the set is unranked, which is less informative than a probability distribution for exploratory analysis. (2) Conservative: the isoform set may include rare isoforms that are biologically unlikely but grammatically valid, leading to more UNCERTAIN verdicts. (3) Regular language limitation: finite-state transducers cannot express arbitrary long-range dependencies, requiring approximations for exon definition interactions across large introns. (4) Set size can grow exponentially for genes with many alternative exons (mitigated by hard cap on isoform count (default 100) to prevent combinatorial explosion).
