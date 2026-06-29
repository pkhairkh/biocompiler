# ADR-0005: Three-Valued Logic for Verdicts

> **Note**: The system uses five-valued logic (PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL). This ADR describes the three-valued subset used in the Lean4 formal model (PASS/FAIL/UNCERTAIN).

## Status

Superseded

## Date

2026-05-30

## Context

The type system must produce verdicts for each biological correctness property checked against an mRNA sequence. The verdict must honestly represent the system's state of knowledge: guaranteed correct, guaranteed incorrect, or cannot determine. The verdict must also compose: when multiple properties are checked, the combined verdict must accurately represent the combined state of knowledge without requiring independence assumptions or probability estimates.

**Alternatives Considered:**

1. **Probability scores (0.0–1.0)** — Fine-grained, quantitative, familiar. But: requires calibration (predicted probabilities must match observed frequencies); does not compose without independence assumptions (P(A and B) requires knowing P(A|B)); training data needed for probability estimation; a "0.97 safe" verdict is less useful than a "guaranteed safe" verdict for gene therapy.

2. **Binary pass/fail** — Simple, unambiguous. But: hides uncertainty — a sequence that "probably" satisfies a property but cannot be proven to do so gets the same FAIL verdict as a sequence that definitely violates it. This is too conservative and leads to low acceptance rates for valid designs.

3. **Confidence intervals** — Quantitative representation of uncertainty. But: still requires calibration and training data; still requires independence assumptions for composition; the interval width is itself an estimate.

4. **Three-valued logic (PASS/FAIL/UNCERTAIN)** — Honest representation of three epistemic states: guaranteed correct, guaranteed incorrect, or cannot determine. Composes without independence assumptions. No calibration or training data needed.

## Decision

Use three-valued logic PASS/FAIL/UNCERTAIN (Alternative 4). Three values capture exactly the information we have: guaranteed correct (PASS), guaranteed incorrect (FAIL), or cannot determine (UNCERTAIN). The logic composes cleanly: PASS ∧ PASS = PASS; PASS ∧ UNCERTAIN = UNCERTAIN; FAIL ∧ anything = FAIL. No independence assumptions, no calibration, no training data. An UNCERTAIN verdict is more interpretable than P = 0.63: it tells the user "the framework cannot determine the answer from the available information; you need more specific constraints or more knowledge about the cellular context."

## Consequences

- Positive: (1) Honest about uncertainty — UNCERTAIN is a legitimate, meaningful verdict, not a failure. (2) Composable without independence assumptions. (3) No calibration or training data needed. (4) Interpretable — users understand PASS/FAIL/UNCERTAIN more intuitively than probability values. (5) Maps directly to abstract interpretation's YES/NO/MAYBE.
- Negative: (1) Less actionable than probabilities — UNCERTAIN provides less guidance than P = 0.63 for prioritizing designs. (2) Conservative — may produce UNCERTAIN when the true answer is PASS or FAIL, leading to unnecessary additional analysis. (3) Cannot compare designs quantitatively — "Design A: PASS" and "Design B: PASS" look equivalent even if Design A is more robust. (4) Users may find UNCERTAIN frustrating if it occurs frequently for their use case.

> **Note:** The system now uses five-valued logic (PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL). See FiveValued.lean for formal proofs. The three-valued model is retained as the Lean4 formal foundation.
