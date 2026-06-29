# ADR-0004: Constraint Satisfaction for Gene Design

## Status

Superseded by ADR-0008

## Date

2026-05-30

## Supersession Note

This ADR specified z3 CSP as the primary optimization strategy. In practice, z3 does not scale beyond ~100 amino acids and handles multi-codon constraints poorly. ADR-0008 replaces z3 with a greedy multi-step optimizer as default, keeping z3 as an optional backend for short proteins. The type-directed mutagenesis engine (ADR-0009) extends the approach to handle mathematical impossibility at the codon level.

## Context

The gene design problem is to find an mRNA sequence that translates to a target protein and satisfies multiple biological correctness constraints: splicing correctness (no cryptic splice sites, correct splicing pattern), codon adaptation (CAI above a threshold), GC content (within a specified range), restriction site absence, reading frame consistency, and instability motif absence. The problem is an inverse problem: given a target protein and a set of constraints, find an mRNA sequence that satisfies all constraints.

**Alternatives Considered:**

1. **Multi-objective weighted optimization** — Combine constraints into a single weighted score (e.g., 0.3×CAI + 0.2×GC_score + 0.5×splice_safety_score) and maximize. Flexible; can trade off between objectives. But: weights are arbitrary and require tuning; no hard guarantees — a high-scoring solution may violate individual constraints; cannot provide formal guarantees needed for safety-critical applications.

2. **Genetic algorithm** — Evolve a population of mRNA sequences, selecting for fitness. Handles non-convex search spaces. But: no convergence guarantee; no formal guarantees on solution quality; non-deterministic across runs; cannot diagnose infeasibility.

3. **Simulated annealing** — Stochastic search that avoids local optima. But: same problems as genetic algorithm: no guarantees, non-deterministic, no diagnosis.

4. **Constraint Satisfaction Problem (CSP) with MUS diagnosis** — Formulate gene design as a CSP: variables are codon positions, domains are synonymous codons, constraints are biological correctness properties. Solver finds a feasible assignment or reports INFEASIBLE with a Minimal Unsatisfiable Subset (MUS) explaining the conflict. Deterministic, sound, complete (for tractable instances).

## Decision

Formulate gene design as a Constraint Satisfaction Problem (Alternative 4). Safety-critical gene design requires hard guarantees, not soft scores. A gene designed for therapeutic use must be guaranteed to splice correctly, not merely likely to splice correctly. The CSP formulation provides three properties that no optimization approach can match: (1) Soundness: every solution satisfies all constraints by construction. (2) Completeness: if a feasible assignment exists, the solver finds one (for tractable instances with small domain sizes). (3) Diagnosis: the MUS computation tells the user exactly which constraints conflict, enabling targeted constraint relaxation rather than blind parameter tuning.

## Consequences

- Positive: (1) Hard guarantees: every returned solution satisfies all constraints. (2) Completeness for feasible instances. (3) Diagnostic MUS when infeasible — tells user exactly what is wrong. (4) Deterministic: same input always produces same solution or same INFEASIBLE report. (5) No training data, no parameter tuning, no arbitrary weights.
- Negative: (1) May report INFEASIBLE when a "good enough" solution exists — the hard-constraint formulation is less flexible than soft scoring. (2) CSP is NP-hard in general; proteins longer than ~1000 amino acids may exceed the 60-second time budget. (3) The feasible set may be large, and selecting among feasible solutions requires a secondary objective (currently CAI maximization), which adds complexity. (4) Not suitable for exploratory design where users want to see a range of trade-offs — the CSP produces one (or a few) feasible solutions, not a Pareto front.
