# ADR-0008: Greedy Multi-Objective Optimizer as Default

## Status

Accepted

## Date

2026-05-31

## Context

The original architecture (ADR-0004) specified Constraint Satisfaction Problem (CSP) with z3 as the primary solving strategy. However, in practice, the z3 CSP solver faces fundamental limitations for gene design:

1. **Scalability**: z3's SAT-based solver has exponential worst-case complexity. Proteins longer than ~100 amino acids regularly exceed the 60-second time budget, making z3 impractical for typical therapeutic genes (200-1500 amino acids).

2. **Multi-codon interactions**: Restriction sites that straddle codon boundaries (e.g., PstI CTG|CAG) require coordinated multi-codon changes. The z3 formulation encodes each codon independently, making cross-codon constraints expensive to express.

3. **Cryptic splice site elimination**: MaxEntScan scoring is a continuous function that doesn't decompose cleanly into Boolean constraints. Forcing it into z3's Boolean framework requires approximations that lose precision.

4. **Step ordering advantage**: The greedy optimizer's sequential step ordering (CAI → restriction sites → ATTTA → T-runs → GC → cryptic splice → reconciliation) exploits the structure of the problem — hard constraints are solved first, soft constraints second, with reconciliation passes to prevent regression.

**Alternatives Considered:**

1. **z3 as default** — Theoretical soundness and completeness guarantees. But: doesn't scale beyond ~100 amino acids; poor handling of continuous objectives (CAI, GC); cryptic splice scoring doesn't map to Boolean constraints; 10-100x slower than greedy on real proteins.

2. **Genetic algorithm** — Handles non-convex search spaces; can optimize CAI directly. But: non-deterministic; no convergence guarantee; no infeasibility diagnosis; can't provide formal guarantees.

3. **Greedy multi-step with coordinated solving** — Deterministic; fast (O(n) per step); handles multi-codon sites; scales to any protein length; reconciliation passes prevent step interference. But: not complete — may miss feasible solutions; no formal optimality guarantee.

## Decision

Use the greedy multi-step optimizer as default (Alternative 3), with z3 available as an optional backend for short proteins where formal completeness is desired. The greedy optimizer's architecture is:

- **CAI maximization**: Best codon per position (maximize CAI)
- **Restriction site removal**: Remove restriction sites (multi-codon coordinated solving)
- **ATTTA motif removal**: Remove ATTTA instability motifs
- **T-run fixing**: Fix 6+ consecutive T runs
- **GC content adjustment**: Adjust GC content (organism-specific targeting)
- **Reconciliation**: Ensure GC adjustment didn't reintroduce restriction sites
- **Cryptic splice elimination**: Eliminate cryptic splice donor/acceptor sites (MaxEntScan context disruption)
- **Final reconciliation**: Ensure splice fixes didn't reintroduce restriction sites

When the greedy optimizer fails to satisfy all predicates, the type-directed mutagenesis engine (ADR-0009) proposes conservative amino acid substitutions to make constraint satisfaction possible.

## Consequences

- Positive: (1) Scales to any protein length — 1500 AA proteins optimize in seconds. (2) Deterministic — same input always produces same output. (3) Multi-codon coordinated solving handles cross-boundary restriction sites correctly. (4) Fast — typically 100-1000x faster than z3 for real proteins. (5) Good practical CAI — typically >0.90 for human codon usage. (6) Reconciliation passes prevent step interference. (7) When it fails, failure is diagnosable via warnings.
- Negative: (1) Not complete — may miss feasible solutions that require non-greedy codon choices. (2) No formal optimality guarantee on CAI — the greedy strategy maximizes CAI in the CAI maximization step but subsequent steps may reduce it. (3) Step ordering is fixed — different orderings might produce better results for some proteins. (4) The reconciliation passes add complexity but may not catch all interactions.
