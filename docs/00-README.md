# BioCompiler — Software Engineering Documentation Set

## Document Identification

| Field | Value |
|---|---|
| **Project** | BioCompiler: A Compiler Framework for Deterministic Gene Design |
| **Document Set Version** | 1.0.0-draft |
| **Status** | ROUGH DRAFT — Not reviewed, not baselined |
| **Classification** | Internal / Research |
| **Date** | 2026-05-30 |

## Purpose

This document set specifies a software system that applies compiler-engineering design patterns to the formalizable stages of gene-to-protein processing (splicing and translation). The system provides deterministic guarantees — formal, machine-checkable proofs of correctness — for designed mRNA sequences, without relying on probabilistic models.

This is NOT an academic paper. It is an engineering specification: every requirement is traceable, every interface is defined, every verification criterion is testable, and every design decision is justified against alternatives.

## Document Map

| ID | Document | Standard | Purpose |
|---|---|---|---|
| `DOC-00` | README (this document) | — | Document map, conventions, reading order |
| `DOC-01` | Software Requirements Specification (SRS) | IEEE 830 / ISO/IEC/IEEE 29148 | What the system must do, shall-not-do, and under what constraints |
| `DOC-02` | Software Architecture Document (SAD) | ISO/IEC/IEEE 42010 | How the system is decomposed, what communicates with what, and why |
| `DOC-03` | Software Design Document (SDD) | IEEE 1016 | Detailed design of each component: algorithms, data structures, invariants |
| `DOC-04` | Interface Control Document (ICD) | MIL-STD-2549 / IEEE 1320 | Exact contracts for every IR schema, API, and FFI boundary |
| `DOC-05` | Software Verification & Validation Plan (SVVP) | IEEE 1012 | How we prove the system meets every requirement |
| `DOC-06` | Design Rationale (DR) | ISO/IEC/IEEE 42010 §7 | Why each design decision was made, what alternatives were considered, and what was rejected |
| `DOC-07` | Project Plan (PP) | IEEE 16326 / PMBOK | Work breakdown, schedule, risk register, staffing |
| `DOC-08` | Traceability Matrix (TM) | ISO/IEC/IEEE 15289 | Complete bidirectional traceability: requirements ↔ architecture ↔ design ↔ interface ↔ test ↔ risk |
| `DOC-09` | Critical Analysis of Original Framework | — | Identifies the nine fatal flaws in the original "Compiler for Protein Synthesis" concept and extracts the salvageable value |
| `DOC-10` | Deterministic Methods for Non-Deterministic Biology | — | Six formal methods that produce deterministic answers from non-deterministic biological systems without probability |
| `ADR/` | Architecture Decision Records | Nygard ADR | 7 standalone ADRs in Nygard format (ADR-0001 through ADR-0007) |

## Reading Order

### For Developers Building the System

1. `adr/` (Architecture Decision Records) — Read the formal Nygard-format ADRs for each major decision.
2. `DOC-06` (Design Rationale) — Understand why the system is designed this way (extended narrative with full alternatives analysis).
3. `DOC-01` (SRS) — Understand every SHALL and SHALL NOT.
4. `DOC-02` (SAD) — Understand the component decomposition and data flow.
5. `DOC-03` (SDD) — Understand the detailed algorithms and invariants.
6. `DOC-04` (ICD) — Implement against the exact interface contracts.
7. `DOC-05` (SVVP) — Write tests that map to requirement IDs.

### For Reviewers / Stakeholders

1. `DOC-09` (Critical Analysis) — Understand what was wrong with the original idea and what remains.
2. `DOC-10` (Deterministic Methods) — Understand the theoretical foundations.
3. `DOC-01` (SRS) Sections 1–3 — Scope, stakeholders, and functional requirements.
4. `DOC-02` (SAD) Section 2 — Architecture rationale (why this design).
5. `DOC-07` (Project Plan) — Schedule, risk, and resource estimates.

### For Quality Assurance

1. `DOC-05` (SVVP) — The test plan, with traceability to every requirement.
2. `DOC-08` (Traceability Matrix) — Full bidirectional traceability.
3. `DOC-01` (SRS) Section 5 — Acceptance criteria per requirement.

### For Anyone Evaluating the Idea

1. `DOC-09` (Critical Analysis) — What's wrong with the original and what's salvageable.
2. `DOC-10` (Deterministic Methods) — How to get deterministic guarantees without probabilistic models.
3. `DOC-06` (Design Rationale) — Why each design decision was made.

## Requirement Traceability

Every requirement in `DOC-01` has a unique ID (e.g., `REQ-FUNC-001`). These IDs are traced forward into:

- **Architecture**: Which component satisfies the requirement (`DOC-02`)
- **Design**: Which algorithm and data structure implement the requirement (`DOC-03`)
- **Interface**: Which API/IR field implements the requirement (`DOC-04`)
- **Verification**: Which test case validates the requirement (`DOC-05`)
- **Rationale**: Why the requirement exists and what alternatives were considered (`DOC-06`)
- **Risk**: Which risks threaten the requirement (`DOC-07`)

The complete traceability matrix is in `DOC-08`.

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

## Status of This Draft

This is a **rough draft**. It has not been through:

- [ ] Technical review
- [ ] Stakeholder review
- [ ] Formal baseline
- [ ] Configuration control

Known gaps are marked with `[DRAFT: ...]` inline. Requirements marked `[TBD]` need further elicitation.
