# Architecture Decision Records

This directory contains the Architectural Decision Records (ADRs) for the BioCompiler project, following the Michael Nygard ADR format.

## Index

| ADR | Title | Status | Date |
|---|---|---|---|
| ADR-0001 | Pipeline Architecture | Accepted | 2026-05-30 |
| ADR-0002 | Protocol Buffers for IR Schemas | Accepted | 2026-05-30 |
| ADR-0003 | Non-Deterministic Finite-State Transducers for Splicing | Accepted | 2026-05-30 |
| ADR-0004 | Constraint Satisfaction for Gene Design | Superseded by ADR-0008 | 2026-05-30 |
| ADR-0005 | Three-Valued Logic for Verdicts | Accepted | 2026-05-30 |
| ADR-0006 | Foreign Function Interface for Folding and PTMs | Accepted | 2026-05-30 |
| ADR-0007 | Declarative Grammar Configuration | Accepted | 2026-05-30 |
| ADR-0008 | Greedy Multi-Objective Optimizer as Default | Accepted | 2026-05-31 |
| ADR-0009 | Type-Directed Protein Mutagenesis | Accepted | 2026-05-31 |
| ADR-0010 | Graduated Certificates | Accepted | 2026-05-31 |
| ADR-0011 | GT-Free Codon Prioritization in Cryptic Splice Elimination | Accepted | 2026-03-05 |
| ADR-0012 | CpG Avoidance in Greedy Optimizer | Accepted | 2026-03-05 |
| ADR-0013 | Mutagenesis GT-Mandatory vs Optimizer Weakness Distinction | Accepted | 2026-03-05 |
| ADR-0014 | Predicate Checking Delegation (Separation of Concerns) | Accepted | 2026-05-31 |
| ADR-0015 | Biosecurity Sequence Screening | Accepted | 2026-03-05 |
| ADR-0016 | Default Safety Measures | Accepted | 2026-03-05 |
| ADR-0017 | Feature Parity with DNAchisel | Accepted | 2026-03-05 |
| ADR-0018 | tRNA Adaptation Index (tAI) | Accepted | 2026-03-05 |

## Relationship to DOC-06 (Design Rationale)

The full design rationale — including the critical analysis of the original proposal, the six deterministic methods, and the detailed reasoning behind each decision — is in `../06-Design-Rationale.md`. The ADRs in this directory provide the formal decision records in the Nygard format; DOC-06 provides the extended narrative rationale.

## Format

Each ADR follows the Nygard format:
- **Title**: A short noun phrase describing the decision
- **Status**: Proposed, Accepted, Deprecated, or Superseded
- **Date**: The date the decision was made
- **Context**: The forces at play, including technological, political, social, and project-local
- **Decision**: The decision that was made
- **Consequences**: The resulting context after applying the decision, including positive and negative
