# ADR-0016: Default Safety Measures

## Status: Accepted

## Date

2026-03-05

## Context

Previously, BioCompiler's safety features were opt-in. Provenance tracking, strict mode (refusing results with failed predicates), and translation verification each required explicit user configuration to enable. The default behavior produced results even when predicates failed, with no warning, no verification, and no audit trail.

This created several problems:

1. **Silent predicate failures**: A user could receive an optimized sequence with failed predicates (e.g., cryptic splice sites, CpG islands, restriction sites) without any indication that the result was unsafe. The `result.predicates_failed` list was populated but never surfaced unless the user explicitly checked it.

2. **No audit trail**: Without provenance tracking, there was no record of which optimization steps were applied, what parameters were used, or what constraints were active. This made it impossible to reproduce results or trace the origin of a problematic sequence.

3. **Unverified translations**: The optimizer could produce a DNA sequence that does not translate back to the input protein — a catastrophic correctness failure — without any automated check. This could happen due to bugs in the optimization pipeline or unexpected interactions between constraint steps.

4. **Missing biosafety context**: Exports (GenBank, FASTA, JSON) contained no biosafety metadata, meaning that downstream consumers (synthesis providers, lab information systems) had no way to know whether the sequence had been screened or verified.

The principle of **secure by default** — well-established in software engineering (e.g., HTTPS by default, secure defaults in cloud IAM) — should apply equally to biological design tools. A tool that produces DNA sequences should be safe, verified, and auditable without requiring the user to remember to enable each safety feature individually.

**Alternatives Considered:**

1. **Keep opt-in safety** — Continue requiring explicit configuration for safety features. But: the current opt-in approach has demonstrably failed — users routinely omit safety settings, and the tool has produced unsafe results in production use.

2. **Deprecation period with warnings** — Emit deprecation warnings when safety features are disabled, then change defaults in a later version. But: this delays the fix and still requires users to act. The silent-failure risk persists during the deprecation period.

3. **Safe by default with explicit opt-out** — Enable all safety features by default. Users who need the previous behavior must explicitly disable each one. This follows the secure-by-default principle and makes unsafe usage a conscious decision.

## Decision

Adopt safe-by-default with explicit opt-out (Alternative 3). The following safety measures are now enabled by default:

1. **Provenance tracking ON by default**: Every optimization result includes a `provenance` record containing:
   - Input protein sequence and parameters
   - Optimization steps applied (in order, with per-step parameters)
   - Constraint configuration
   - Timestamp and BioCompiler version
   - Hash of the input for tamper detection

   The `provenance` field is always populated. Previously it was populated only when `track_provenance=True` was set. Users who wish to disable provenance tracking (e.g., for performance-sensitive batch jobs) can set `track_provenance=False`.

2. **Strict mode ON by default**: The optimizer refuses to return results where any predicate has failed. Previously, `strict_mode=False` (the default) would return results with `predicates_failed` populated but no error raised. Now, `strict_mode=True` by default means that `optimize_sequence()` raises `OptimizationError` when predicates fail. Users who need partial results can set `strict_mode=False`.

3. **Translation verification automatic post-optimization**: After optimization completes, the result DNA sequence is automatically translated and compared to the input protein. If the translation does not match, an `OptimizationError` is raised. Previously, this check was only performed when `verify_translation=True` was set. The check is now always performed and cannot be disabled — producing a sequence that doesn't translate to the input protein is never acceptable.

4. **Biosafety annotations in all exports**: Every export format (GenBank, FASTA, JSON) includes biosafety metadata:
   - Biosecurity screen result (from ADR-0015)
   - Predicate pass/fail status
   - Provenance hash
   - Optimization parameters summary

   These annotations are always included. There is no option to omit them.

## Consequences

- **Positive**: Every design produced by BioCompiler is auditable (provenance), verified (translation check), and safe (strict mode + biosafety annotations). This aligns with responsible biodesign principles and the WHO Laboratory Biorisk Management guidance.
- **Positive**: Downstream consumers of BioCompiler output (synthesis providers, lab LIMS) have the metadata they need to make informed decisions about the sequence, without requiring the designer to remember to include it.
- **Positive**: Translation verification as an always-on check provides a critical correctness guarantee. A sequence that doesn't translate to the input protein is a bug, not a valid result — this check catches such bugs at the output boundary.
- **Positive**: Strict mode by default means that silent predicate failures are no longer possible. Users must consciously opt into receiving partial results.
- **Negative**: This is a **breaking change** — users who relied on receiving partial results (with failed predicates) will now get `OptimizationError`. They must explicitly set `strict_mode=False` to get the previous behavior. This justifies a major version bump.
- **Migration**: This change targets **v12.0.0**. The major version bump is justified because:
  - `strict_mode` default change breaks any code that calls `optimize_sequence()` and doesn't handle `OptimizationError`
  - `track_provenance` default change adds overhead to every call
  - Translation verification adds a small but non-zero latency
  - Export format changes (biosafety annotations) may break parsers that expect the old format
- **Negative**: Provenance tracking adds memory and compute overhead. For large batch jobs, the provenance records can be significant. The `track_provenance=False` opt-out provides an escape valve for performance-critical use cases, but the default is now slower.
- **Negative**: Some legitimate workflows (e.g., exploratory optimization where partial results are informative) will be disrupted by strict mode. Users must explicitly opt out, which adds friction.

## References

- ADR-0015: Biosecurity Sequence Screening (screening that feeds into biosafety annotations)
- ADR-0010: Graduated Certificates (predicate verification framework)
- ADR-0014: Predicate Checking Delegation (type system as single source of truth for predicate evaluation)
- OWASP Secure by Design principles
- WHO Laboratory Biorisk Management Framework (2024)
