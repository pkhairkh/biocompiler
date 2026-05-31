# ADR-0007: Declarative Grammar Configuration

## Status

Accepted

## Date

2026-05-30

## Context

The splicing grammar rules (splice site consensus sequences, regulatory element thresholds, branch point motifs, polypyrimidine tract requirements) are scientific knowledge that evolves as new experimental data becomes available. The system must allow these rules to be updated without modifying the source code of the splicing engine. Additionally, different organisms and cell types require different grammar parameters, which must be configurable.

**Alternatives Considered:**

1. **Hardcoded rules in source code** — Simplest implementation; no parsing overhead; rules are part of the codebase. But: updating rules requires modifying source code and redeploying; non-programmers cannot contribute rules; version control of rules is tied to source code commits rather than scientific publication dates; conflicts with the principle that scientific knowledge should be configurable by domain experts.

2. **Learned rules (grammar induction)** — Automatically discover rules from data. But: grammar induction is undecidable for the required grammar class (Flaw #4); ADIOS learns statistics, not physics (Flaw #3); no formal guarantees can be provided for learned rules; conflicts with REQ-CON-002 (no grammar induction).

3. **Declarative configuration files (YAML)** — Rules specified in human-readable YAML files loaded at pipeline initialization. Domain experts can edit rules without touching source code. Grammar changes trigger full regression testing. Different configuration files for different organisms and cell types.

## Decision

Use declarative configuration files in YAML format (Alternative 3). Splicing rules are scientific knowledge that should be updatable by domain experts without modifying source code. YAML is human-readable, supports comments, and is widely used for configuration. The splicing engine loads the grammar rules at initialization time and validates them against a schema. Different configuration files support different organisms (human, mouse, Drosophila) and cell types (HEK293T, HepG2, neuron).

## Consequences

- Positive: (1) Domain experts can update grammar rules without programming. (2) Grammar changes are tracked independently from source code changes. (3) Different configurations for different organisms and cell types. (4) Community contributions: researchers can submit grammar rule updates as pull requests on the configuration files. (5) Schema validation catches errors before pipeline execution.
- Negative: (1) Runtime parsing overhead for YAML files (negligible for the expected file sizes). (2) Configuration validation must be thorough to catch errors (e.g., invalid consensus sequences, inconsistent thresholds). (3) Configuration drift: different users may use different grammar versions, leading to different isoform sets for the same input sequence. (4) YAML is less precise than code — there is no type checking at configuration time (mitigated by JSON Schema validation).
