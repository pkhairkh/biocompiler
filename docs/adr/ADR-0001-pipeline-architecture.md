# ADR-0001: Pipeline Architecture

## Status

Accepted

## Date

2026-05-30

## Context

The BioCompiler system needs to process gene sequences through multiple transformation stages — scanning, splicing, translation, type checking, optimization, and certification — with well-defined interfaces between stages. The stages have different computational characteristics: some are purely symbolic (scanning, translation), some are non-deterministic (splicing), and some require external tools (folding, PTM prediction). The architecture must support independent development, testing, and replacement of stages while maintaining end-to-end correctness guarantees.

**Alternatives Considered:**

1. **Monolithic tool** — All stages implemented as functions in a single program with no formal boundaries. Simpler initial development; no IR serialization overhead. But: stages cannot be developed or tested independently; changes to one stage risk breaking others; no natural boundary for isolating formalizable from non-formalizable components; does not support compositional analysis.

2. **Microservices** — Each stage as an independent service communicating over HTTP/gRPC. Maximum isolation and independent deployment. But: massive over-engineering for a batch processing pipeline; network latency and serialization overhead for every stage transition; operational complexity (service discovery, failure handling, monitoring); no benefit over in-process pipeline for the expected workload (single-user, batch processing).

3. **Pipeline with typed IR** — Staged transformation with typed intermediate representations mediating between stages. Each stage is a pure function that consumes an IR record and produces an IR record. The IR schemas enforce data contracts. Mirrors LLVM's proven architecture.

## Decision

Adopt the pipeline architecture with typed IR (Alternative 3). The pipeline architecture mirrors LLVM's proven multi-pass design (Lattner & Adve, 2004), which demonstrated that a staged pipeline with typed IR enables: (1) independent development of passes by different teams, (2) independent testing of each pass against its input/output contract, (3) replacement of passes without affecting others, (4) compositional analysis where each pass adds information to the IR, and (5) formal verification at each stage via IR invariants. These benefits are directly applicable to BioCompiler. The pipeline architecture also provides the right abstraction for isolating formalizable stages (splicing, translation) from non-formalizable ones (folding, PTMs) behind the FFI boundary, which is essential given Flaw #1 and Flaw #2 from the critical analysis.

## Consequences

- Positive: (1) Modularity — each stage can be developed, tested, and replaced independently. (2) Testability — each stage has a defined input/output contract that can be tested in isolation. (3) Extensibility — new passes can be inserted without modifying existing ones. (4) Formal verification — IR invariants can be checked at each stage boundary. (5) Separation of paradigms — the FFI boundary cleanly isolates symbolic from continuous computation.
- Negative: (1) More upfront design effort is required for IR schemas — each IR level must be specified before stages can be developed in parallel. (2) IR serialization and deserialization add overhead compared to in-process function calls, though this is negligible for the data sizes involved. (3) Schema evolution requires careful version management (see ADR-0002). (4) The pipeline architecture implies a fixed ordering of stages, which may not be optimal for all use cases — the optimizer loop is a special case handled by COMP-06.
