# ADR-0014: Predicate Checking Delegation (Separation of Concerns)

## Status

Accepted

## Date

2026-05-31

## Context

The optimizer module originally contained a `_check_predicates` function that re-implemented predicate logic inline — checking GC content, restriction sites, instability motifs, CAI, and cryptic splice sites using direct computation rather than delegating to the type system. This violated the Single Responsibility Principle: the optimizer is responsible for *finding* sequences, not *verifying* them.

This duplication caused three problems:

1. **Divergent semantics**: The optimizer's `_check_predicates` reported `"NoCrypticSplice"` as a simple pass/fail based on max donor score, while the type system's `evaluate_no_cryptic_splice` evaluated donor-acceptor pairs. This led to the optimizer claiming NoCrypticSplice was satisfied when the type system would disagree, or vice versa.

2. **Maintenance burden**: Every new predicate or change to existing predicate logic required updating both the type system AND the optimizer. This led to inconsistencies (e.g., NoCpGIsland was in the type system but not in the optimizer's check).

3. **Inconsistent reporting**: The optimizer reported predicate names differently from the type system (e.g., `"CodonAdapted"` vs. `"CodonAdapted(Homo_sapiens, 0.2)"`), making it impossible to reconcile results between the two modules.

Additionally, the NoCrypticSplice predicate was only checking donor-acceptor *pairs*, missing standalone strong cryptic donors. A strong GT site above the MaxEntScan threshold is problematic even without a paired acceptor, because it can pair with downstream genomic acceptors in vivo.

**Alternatives Considered:**

1. **Keep duplicate logic** — Simpler optimizer; no cross-module dependency. But: semantic divergence; maintenance burden; inconsistent results.

2. **Move all checking to optimizer** — Optimizer is the single source of truth. But: type system becomes vestigial; certificate verification would depend on optimizer; circular dependency risk.

3. **Delegate to type system** — Optimizer calls type system's `evaluate_all_predicates` for checking. Type system is the single source of truth for all predicate logic. Optimizer focuses solely on sequence generation.

## Decision

Delegate all predicate checking to the type system (Alternative 3). The optimizer's `_check_predicates` is replaced with `_check_predicates_via_type_system`, which calls `evaluate_all_predicates` and converts results to (satisfied, failed) lists.

Additionally, the NoCrypticSplice predicate is updated to flag standalone strong donors and acceptors (not just pairs), aligning it with the optimizer's cryptic splice elimination logic.

## Consequences

- Positive: (1) Single source of truth — predicate semantics are defined once in the type system. (2) NoCpGIsland is now checked in the optimization pipeline. (3) Consistent predicate naming between optimizer and type system. (4) NoCrypticSplice now catches standalone strong donors (the V-codon problem). (5) Certificate verification naturally uses the same logic as the optimizer.
- Negative: (1) Cross-module dependency — optimizer now imports from type_system. (2) Type system evaluation is slower than inline checking (overhead of creating TypeCheckResult objects). (3) The type system's parameterized predicate names (e.g., `"CodonAdapted(Homo_sapiens, 0.2)"`) are more verbose in results.
