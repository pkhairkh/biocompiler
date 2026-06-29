# ADR-0010: Graduated Certificates

## Status

Accepted

## Date

2026-05-31

## Context

The original certificate system generated certificates only when ALL predicates passed (strict mode). This created a binary outcome: either the design is fully certified, or it gets nothing. In practice, many real-world designs have partial compliance — a gene with CAI=0.99 that still has one PstI site is still useful, and the biologist can decide whether PstI is relevant for their cloning workflow.

Additionally, with the introduction of type-directed mutagenesis (ADR-0009), certificates need to document the substitutions that were applied. A certificate for a mutagenized protein must record: (a) that mutagenesis was applied, (b) which substitutions were made, and (c) the BLOSUM62 scores and reasons for each substitution. This provenance is essential for reproducibility and for biologists evaluating whether the substitutions are acceptable.

**Alternatives Considered:**

1. **Strict-only certificates** — Only generate when all predicates pass. Simple; no ambiguity. But: rejects useful partial designs; no documentation for near-feasible designs; biologist cannot make informed trade-offs; too conservative for practical use.

2. **Soft scoring certificates** — Assign a numerical score (0-100) based on how many predicates pass. Quantitative. But: scores are arbitrary; 80/100 does not tell the user WHICH predicates failed; cannot be independently verified in a meaningful way.

3. **Graduated certificates with full documentation** — Generate certificates for any design, documenting all predicate results (PASS, FAIL, or UNCERTAIN) and any mutagenesis applied. The certificate is a complete record of what was verified and what was not. Downstream users can make informed decisions about partial compliance.

## Decision

Adopt graduated certificates with full documentation (Alternative 3). The design is:

1. **Default mode (graduated)**: Certificate is generated even when predicates fail. The `overall_status` field indicates the fraction of passing predicates (e.g., `PARTIAL_7/8`).

2. **Strict mode (optional)**: `require_all_pass=True` raises `CertificateGenerationError` if any predicate fails. This preserves backward compatibility.

3. **Mutagenesis metadata**: When type-directed mutagenesis is applied, the certificate's `provenance` section includes a `mutagenesis` object documenting:
   - `applied`: boolean
   - `n_substitutions`: count
   - `substitutions`: list of {position, from, to, blosum62, reason, predicate}
   - `description`: human-readable summary

4. **Independent verification**: The certificate verification procedure re-evaluates every predicate from scratch using only the sequence and parameters embedded in the certificate. It checks that the claimed verdicts are consistent with re-evaluation, but does NOT require all predicates to pass for graduated certificates.

## Consequences

- Positive: (1) Practical — biologists get documentation even for partial designs. (2) Informative — the certificate tells users exactly which predicates pass and which fail. (3) Actionable — failed predicates guide the user toward targeted fixes. (4) Mutagenesis provenance — certificates document all amino acid substitutions with rationale. (5) Backward compatible — strict mode preserves the old behavior. (6) Verifiable — graduated certificates can still be independently verified; verification checks consistency, not universal passing.
- Negative: (1) More complex certificate format — the `mutagenesis` section adds schema complexity. (2) Risk of misuse — users might treat a PARTIAL certificate as a FULL_PASS, not realizing that some predicates failed. (3) Certificate size increases with mutagenesis metadata. (4) Verification must handle both graduated and strict certificates, adding implementation complexity.
