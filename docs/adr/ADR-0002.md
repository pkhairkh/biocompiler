# ADR-0002: Protocol Buffers for IR Schemas

## Status

Accepted

## Date

2026-05-30

## Context

The IR Bus requires schemas that define the structure of each IR level (IR-Seq, IR-Peptide, IR-Structure, IR-Circuit). These schemas must be enforced at runtime, support backward-compatible evolution, support efficient serialization for persistence and inter-process communication, and ideally support code generation for multiple programming languages.

**Alternatives Considered:**

1. **JSON Schema only** — Human-readable, widely supported, no compilation step. But: no type enforcement at compile time; no efficient binary serialization; no code generation; schema validation is slow for large datasets; JSON parsing is significantly slower than protobuf deserialization.

2. **HDF5** — Optimized for large numerical arrays (e.g., 3D structure coordinates). But: poor support for structured, record-oriented data; no code generation; complex C library dependency; overkill for the data sizes involved (sequences of thousands to millions of nucleotides, not terabytes of array data).

3. **Custom binary format** — Maximum flexibility and performance. But: maximum maintenance burden; no standard tooling; no backward compatibility guarantees; every consumer must implement its own parser; version management is entirely manual.

4. **Protocol Buffers (proto3)** — Schema enforcement, backward compatibility, efficient binary serialization, multi-language code generation (Python, C++, Java, Go), well-supported by Google and the open-source community.

## Decision

Use Protocol Buffers (proto3) for IR schema definitions (Alternative 4). Protocol Buffers provide the best combination of schema enforcement, backward compatibility, efficient serialization, and multi-language code generation. Schema enforcement is automatic: the generated Python classes have typed fields and raise errors on schema violations. Backward compatibility is built into the format: new fields can be added without breaking existing consumers (unknown fields are preserved). Efficient binary serialization reduces IR size by 3–10x compared to JSON. Code generation for Python eliminates manual parsing and validation.

## Consequences

- Positive: (1) Compile-time type safety from generated code. (2) Backward-compatible schema evolution via field numbering. (3) Efficient binary serialization (3–10x smaller than JSON). (4) Multi-language code generation for future ports. (5) Well-understood tooling (protoc, grpc).
- Negative: (1) Build step required for protoc compilation before Python code can run. (2) Proto3 limitations: no required fields (all are optional), no custom default values, limited support for one-of types. (3) Binary format is not human-readable (requires protoc --decode for inspection). (4) Learning curve for developers unfamiliar with protobuf.
