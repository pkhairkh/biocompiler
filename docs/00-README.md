# BioCompiler — Documentation Index

| Field | Value |
|---|---|
| **Project** | BioCompiler: A Compiler Framework for Deterministic Gene Design |
| **Version** | 1.0.0-draft |
| **Status** | ROUGH DRAFT — Not reviewed, not baselined |
| **Date** | 2026-05-30 |

## Specification Documents

| Doc | File | Purpose |
|---|---|---|
| SRS | [01-SRS.md](01-SRS.md) | What the system must do, shall-not-do, and under what constraints (IEEE 830 / ISO/IEC/IEEE 29148) |
| SAD | [02-SAD.md](02-SAD.md) | Component decomposition, communication, and architecture rationale (ISO/IEC/IEEE 42010) |
| SDD | [03-SDD.md](03-SDD.md) | Detailed design: algorithms, data structures, invariants per component (IEEE 1016) |
| ICD | [04-ICD.md](04-ICD.md) | Exact contracts for every IR schema, API, and FFI boundary (MIL-STD-2549 / IEEE 1320) |

## Verification & Validation

| Doc | File | Purpose |
|---|---|---|
| SVVP | [05-SVVP.md](05-SVVP.md) | How we prove the system meets every requirement (IEEE 1012) |
| TM | [08-Traceability-Matrix.md](08-Traceability-Matrix.md) | Bidirectional traceability: requirements ↔ architecture ↔ design ↔ interface ↔ test ↔ risk (ISO/IEC/IEEE 15289) |

## Design & Rationale

| Doc | File | Purpose |
|---|---|---|
| DR | [06-Design-Rationale.md](06-Design-Rationale.md) | Why each design decision was made, alternatives considered, and what was rejected (ISO/IEC/IEEE 42010 §7) |
| PP | [07-Project-Plan.md](07-Project-Plan.md) | Work breakdown, schedule, risk register, staffing (IEEE 16326 / PMBOK) |

## Technical Analysis

| Doc | File | Purpose |
|---|---|---|
| Critical Analysis | [09-Critical-Analysis.md](09-Critical-Analysis.md) | Nine fatal flaws in the original framework concept and salvageable value |
| Deterministic Methods | [10-Deterministic-Methods.md](10-Deterministic-Methods.md) | Six formal methods yielding deterministic answers from non-deterministic biology without probability |
| Refinement Mapping | [11-Refinement-Mapping.md](11-Refinement-Mapping.md) | Stepwise refinement from abstract specification to implementation |
| Engine Accuracy | [12-Engine-Accuracy.md](12-Engine-Accuracy.md) | Accuracy analysis and verification of the compilation engine |

## Validation Results

| Doc | File | Purpose |
|---|---|---|
| Retrospective Validation | [13-Retrospective-Validation.md](13-Retrospective-Validation.md) | Retrospective validation against known biological data |
| SLOT Proof-Implementation Gap | [14-SLOT-Proof-Implementation-Gap.md](14-SLOT-Proof-Implementation-Gap.md) | Analysis of gaps between SLOT proofs and running implementation |
| Technical Reference | [15-Reference.md](15-Reference.md) | 28-predicate tables, unified engine API, TCB, honest limitations |

## Architecture Decision Records

Nygard-format ADRs documenting each major architectural decision: [adr/](adr/)

| ADR | File |
|---|---|
| ADR-0001 | [ADR-0001.md](adr/ADR-0001.md) |
| ADR-0002 | [ADR-0002.md](adr/ADR-0002.md) |
| ADR-0003 | [ADR-0003.md](adr/ADR-0003.md) |
| ADR-0004 | [ADR-0004.md](adr/ADR-0004.md) |
| ADR-0005 | [ADR-0005.md](adr/ADR-0005.md) |
| ADR-0006 | [ADR-0006.md](adr/ADR-0006.md) |
| ADR-0007 | [ADR-0007.md](adr/ADR-0007.md) |
| ADR-0008 | [ADR-0008.md](adr/ADR-0008.md) |
| ADR-0009 | [ADR-0009.md](adr/ADR-0009.md) |
| ADR-0010 | [ADR-0010.md](adr/ADR-0010.md) |
| ADR-0011 | [ADR-0011-gt-free-codon-prioritization.md](adr/ADR-0011-gt-free-codon-prioritization.md) |
| ADR-0012 | [ADR-0012-cpg-avoidance.md](adr/ADR-0012-cpg-avoidance.md) |
| ADR-0013 | [ADR-0013-mutagenesis-gt-mandatory.md](adr/ADR-0013-mutagenesis-gt-mandatory.md) |
| ADR-0014 | [ADR-0014.md](adr/ADR-0014.md) |

## Conventions

| Convention | Meaning |
|---|---|
| **SHALL** | Mandatory requirement |
| **SHALL NOT** | Mandatory prohibition |
| **SHOULD** | Recommended but waivable with documented rationale |
| **MAY** | Optional capability |
| `REQ-FUNC-XXX` | Functional requirement ID |
| `REQ-NFR-XXX` | Non-functional requirement ID |
| `REQ-CON-XXX` | Constraint (design or environmental) |
| `COMP-XX` | Component ID in architecture |
| `IF-XX` | Interface ID in ICD |
| `TC-XX` | Test case ID in SVVP |
| `RISK-XX` | Risk ID in Project Plan |
| `AD-XX` | Architectural decision ID |
| `INV-XXX` | Invariant ID |
| `[DRAFT: ...]` | Known gap requiring further elicitation |
| `[TBD]` | To be determined at baseline review |
