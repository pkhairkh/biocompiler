# ADR-0006: Foreign Function Interface for Folding and PTMs

## Status

Accepted

## Date

2026-05-30

## Context

Protein folding and post-translational modification (PTM) prediction are essential for a complete gene-to-protein pipeline, but they are not formalizable as string transformations (Flaw #1, Flaw #2 from the critical analysis). Folding is a thermodynamic process; PTMs are co-determined by the cellular environment. The system must provide access to these capabilities without compromising the deterministic guarantees of the core pipeline.

**Alternatives Considered:**

1. **Internal ML model** — Train or integrate a neural network (e.g., AlphaFold-style) as a core pipeline component. But: introduces probabilistic computation into the deterministic core; requires GPU resources for all pipeline runs; model updates require re-validation of the entire pipeline; conflicts with REQ-CON-001 (no probabilistic models for internal stages).

2. **Internal physics simulation** — Implement molecular dynamics (MD) for folding. But: computationally infeasible for routine use (hours to days per protein); requires extensive force field parameterization; still approximate (no guarantee of correct structure); conflicts with the pipeline's performance requirements.

3. **Omit folding/PTM entirely** — Restrict the system to splicing + translation + verification. Simpler; no external dependencies. But: severely limits the system's utility — many users need structural information; misses the opportunity to provide typed IR for structure data.

4. **Foreign Function Interface (FFI) with SLOT-filling** — Invoke external tools (AlphaFold, NetPhos, etc.) through a defined adapter interface. The core pipeline treats FFI outputs as non-deterministic SLOT fields in the IR. The FFI boundary isolates non-determinism from the deterministic core.

## Decision

Implement a Foreign Function Interface with SLOT-filling (Alternative 4). The FFI boundary isolates the non-determinism of external tools from the deterministic core pipeline. The IR schema defines empty SLOT fields (e.g., secondary_structure_pred, ptm_sites) that are filled by FFI adapters. The core pipeline does not model the internal computation of any external tool; it only guarantees correct input formatting, correct output parsing, and provenance metadata preservation. FFI invocations are explicitly treated as non-deterministic (REQ-FUNC-053).

## Consequences

- Positive: (1) Uses best-in-class external tools (AlphaFold, NetPhos) without reimplementing them. (2) Clean separation of paradigms — deterministic core and non-deterministic periphery are isolated. (3) No GPU requirement for core pipeline — FFI stages are optional (REQ-NFR-042). (4) Extensible — new external tools can be added by implementing the adapter interface. (5) Provenance tracking for all FFI outputs.
- Negative: (1) Dependency on external tools — if AlphaFold changes its API or output format, the adapter must be updated. (2) No formal guarantees on FFI output — the system can only guarantee correct input formatting and output parsing, not correctness of the prediction. (3) FFI stages are non-deterministic, which means pipeline output is not fully reproducible when FFI stages are included. (4) External tool installation and configuration add deployment complexity.
