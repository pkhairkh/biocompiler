# DOC-08: Traceability Matrix

| Field | Value |
|---|---|
| **Document ID** | DOC-08 |
| **Version** | 1.0.0-draft |
| **Status** | ROUGH DRAFT |
| **Date** | 2026-05-30 |
| **Prepared By** | BioCompiler Project Team |
| **Reviewed By** | [TBD at baseline review] |
| **Approved By** | [TBD at baseline review] |
| **Standard** | ISO/IEC/IEEE 15289 |

---

## 1. Introduction

### 1.1 Purpose

This document constitutes the Traceability Matrix for the **BioCompiler** system, prepared in accordance with ISO/IEC/IEEE 15289 (Systems and software engineering — Content of life-cycle information items). Its purpose is to provide **complete bidirectional traceability** from stakeholder needs through requirements, architecture, design, interfaces, test cases, and risks — ensuring that every requirement is implemented, verified, and traceable, and that every test case maps back to one or more requirements.

Traceability serves three critical functions in the BioCompiler project:

1. **Completeness assurance**: Every requirement in DOC-01 (SRS) has a corresponding implementation artifact (component, design section, interface) and verification artifact (test case). Gaps in traceability indicate either unimplemented requirements or unverified implementations.

2. **Impact analysis**: When a requirement changes, the traceability matrix identifies all affected artifacts — which component must be modified, which design section must be updated, which interface contract changes, and which test cases must be revised or added.

3. **Safety argumentation**: BioCompiler is a safety-critical system: a false PASS verdict from the type system (certifying a defective gene design) could have serious consequences. The traceability matrix provides the evidentiary chain linking each safety-relevant requirement (especially REQ-NFR-011: soundness) to its soundness argument, its implementation, and its adversarial test cases.

### 1.2 Scope

This document covers:

- **Forward traceability** (Sections 2–6): Requirement → Architecture (component) → Design (section) → Interface (API/IR) → Test (unit, integration, system, soundness) → Risk
- **Non-functional requirements traceability** (Section 3): NFR → Verification method → Test case → Risk
- **Constraint traceability** (Section 4): Constraint → Enforcement mechanism → Verification → Risk
- **Architectural decision traceability** (Section 5): Decision → Requirements addressed → Components affected → Rationale → Risk
- **IR invariant traceability** (Section 6): Invariant → IR level / component → Enforcement mechanism → Test
- **Test coverage summary** (Section 7): Aggregate statistics on requirement coverage by test level
- **Backward traceability index** (Section 8): Test → Requirement, Interface → Requirement, Component → Requirement

All requirement identifiers follow the conventions established in DOC-00 (README): `REQ-FUNC-XXX` for functional requirements, `REQ-NFR-XXX` for non-functional requirements, and `REQ-CON-XXX` for constraints. Component identifiers use `COMP-XX`, interface identifiers use `IF-XX`, and test case identifiers use `TC-` prefixes.

### 1.3 Traceability Method

The traceability method employed is **full bidirectional mapping** with the following properties:

- **Every requirement** has at least one forward trace to an architectural component, a design section, and an interface.
- **Every requirement** has at least one forward trace to a test case at some level (unit, integration, system, or soundness).
- **Every test case** has a backward trace to at least one requirement.
- **Every interface** has a backward trace to at least one requirement.
- **Every component** has a backward trace to at least one requirement.
- **Traceability gaps** (requirements with no test, tests with no requirement) are flagged as `[DRAFT: GAP]` and must be resolved before baseline.

The traceability is maintained in tabular form for direct lookup. Cross-references to the source documents (DOC-01 through DOC-07) are provided in every table row.

### 1.4 References

| Ref ID | Document | Description |
|---|---|---|
| REF-01 | DOC-01: Software Requirements Specification | Source of all functional and non-functional requirements |
| REF-02 | DOC-02: Software Architecture Document | Source of component decomposition and architectural decisions |
| REF-03 | DOC-03: Software Design Document | Source of detailed design sections and invariants |
| REF-04 | DOC-04: Interface Control Document | Source of interface specifications (IF-01 through IF-10, IF-DATA-01 through IF-DATA-04) |
| REF-05 | DOC-05: Software Verification & Validation Plan | Source of test case definitions |
| REF-06 | DOC-06: Design Rationale | Source of architectural decision rationale |
| REF-07 | DOC-07: Project Plan | Source of risk register entries |
| REF-08 | ISO/IEC/IEEE 15289:2015 | Systems and software engineering — Content of life-cycle information items |

---

## 2. Functional Requirements Traceability

The following table provides forward traceability for every functional requirement in DOC-01. Each row maps a requirement to the component that implements it (DOC-02), the design section that specifies it (DOC-03), the interface that exposes it (DOC-04), the test cases that verify it (DOC-05), and any associated risks (DOC-07).

### 2.1 Scanner (REQ-FUNC-001 through REQ-FUNC-004)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-001 | Accept DNA/RNA FASTA or raw input; reject non-IUPAC characters; normalize T→U | COMP-01 (Scanner) | §3.1 COMP-01: Scanner | IF-01 (Scanner Interface) | TC-U-001, TC-U-002, TC-U-003 | TC-I-001 | TC-S-001 | — | RISK-04 |
| REQ-FUNC-002 | Scan and annotate biological elements (start/stop codons, splice sites, branch points, polypyrimidine tracts, Kozak, instability motifs, restriction sites) | COMP-01 (Scanner) | §3.1 COMP-01: Scanner | IF-01 (Scanner Interface) | TC-U-004, TC-U-005, TC-U-006 | TC-I-001 | TC-S-001 | — | — |
| REQ-FUNC-003 | Implement scanner as collection of DFAs, one per element type; O(n) per DFA | COMP-01 (Scanner) | §3.1 COMP-01: Scanner | IF-01 (Scanner Interface) | TC-U-003 | TC-I-001 | TC-S-001 | — | — |
| REQ-FUNC-004 | Scanner determinism: byte-identical output for identical input across runs and platforms | COMP-01 (Scanner) | §3.1 COMP-01: Scanner | IF-01 (Scanner Interface) | TC-U-007 | TC-I-001 | TC-S-001 | — | — |

### 2.2 Splicing Engine (REQ-FUNC-010 through REQ-FUNC-015)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-010 | NDFST for splice isoform set; exhaustive enumeration of all parse paths | COMP-02 (Splicing Engine) | §3.2 COMP-02: Splicing Engine | IF-02 (Splicing Interface) | TC-U-010 | TC-I-001 | TC-S-001 | TC-SND-002 | RISK-01 |
| REQ-FUNC-011 | Splicing grammar rules: GT-AG, GC-AG, AT-AC, branch points, polypyrimidine tracts, ESE/ESS/ISE/ISS | COMP-02 (Splicing Engine) | §3.2 COMP-02: Splicing Engine | IF-02 (Splicing Interface) | TC-U-012, TC-U-013 | TC-I-001 | TC-S-001 | — | RISK-01 |
| REQ-FUNC-012 | Alternative splicing as non-deterministic branching; each parse path = one isoform | COMP-02 (Splicing Engine) | §3.2 COMP-02: Splicing Engine | IF-02 (Splicing Interface) | TC-U-011 | TC-I-001 | TC-S-001 | — | — |
| REQ-FUNC-013 | Cellular context parameter modulates ESE/ESS/ISE/ISS thresholds per cell type | COMP-02 (Splicing Engine) | §3.2 COMP-02: Splicing Engine | IF-02 (Splicing Interface) | TC-U-014 | TC-I-001 | TC-S-001 | — | RISK-02 |
| REQ-FUNC-014 | NDFST determinism: same input + context → same isoform set | COMP-02 (Splicing Engine) | §3.2 COMP-02: Splicing Engine | IF-02 (Splicing Interface) | TC-U-015 | TC-I-001 | TC-S-001 | — | — |
| REQ-FUNC-015 | NDFST output specification: spliced mRNA, exon boundaries, reading frame, annotations, provenance | COMP-02 (Splicing Engine) | §3.2 COMP-02: Splicing Engine | IF-02 (Splicing Interface) | TC-U-010 | TC-I-001 | TC-S-001 | — | — |

### 2.3 Translation Engine (REQ-FUNC-020 through REQ-FUNC-023)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-020 | Deterministic FST for translation; standard genetic code; frame-aware | COMP-03 (Translation Engine) | §3.3 COMP-03: Translation Engine | IF-03 (Translation Interface) | TC-U-020 | TC-I-002 | TC-S-001 | — | — |
| REQ-FUNC-021 | Special cases: selenocysteine (UGA+SECIS), pyrrolysine (UAG archaeal), programmed frameshifts | COMP-03 (Translation Engine) | §3.3 COMP-03: Translation Engine | IF-03 (Translation Interface) | TC-U-021, TC-U-023 | TC-I-002 | TC-S-003 | — | — |
| REQ-FUNC-022 | Translation determinism: same input → same amino acid output | COMP-03 (Translation Engine) | §3.3 COMP-03: Translation Engine | IF-03 (Translation Interface) | TC-U-022 | TC-I-002 | TC-S-001 | — | — |
| REQ-FUNC-023 | Translation output: amino acid sequence, codon assignments, Sec flags, frameshift warnings, initiation confidence | COMP-03 (Translation Engine) | §3.3 COMP-03: Translation Engine | IF-03 (Translation Interface) | TC-U-020 | TC-I-002 | TC-S-001 | — | — |

### 2.4 Type System (REQ-FUNC-030 through REQ-FUNC-035)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-030 | Type checker with three-valued verdicts (PASS/FAIL/UNCERTAIN); no probabilities | COMP-05 (Type System) | §3.5 COMP-05: Type System | IF-05 (Type System Interface) | TC-U-030, TC-U-031 | TC-I-004 | TC-S-001 | TC-SND-001, TC-SND-002, TC-SND-003, TC-SND-004, TC-SND-005, TC-SND-006 | RISK-03 |
| REQ-FUNC-031 | Seven type predicates: SpliceCorrect, NoCrypticSplice, CodonAdapted, GCInRange, NoRestrictionSite, InFrame, NoInstabilityMotif | COMP-05 (Type System) | §3.5 COMP-05: Type System | IF-05 (Type System Interface) | TC-U-030, TC-U-031 | TC-I-004 | TC-S-001 | TC-SND-001, TC-SND-002, TC-SND-003, TC-SND-004, TC-SND-005, TC-SND-006, TC-SND-007 | — |
| REQ-FUNC-032 | Type checker determinism: same input → same verdicts and derivation traces | COMP-05 (Type System) | §3.5 COMP-05: Type System | IF-05 (Type System Interface) | TC-U-032 | TC-I-004 | TC-S-001 | — | — |
| REQ-FUNC-033 | Subtyping relations: SpliceCorrect(C1) ⟹ SpliceCorrect(C2) if C1 <: C2; GCInRange narrower ⟹ wider | COMP-05 (Type System) | §3.5 COMP-05: Type System | IF-05 (Type System Interface) | TC-U-033 | TC-I-004 | TC-S-001 | — | — |
| REQ-FUNC-034 | Three-valued composition rules: AND/OR truth tables for PASS/FAIL/UNCERTAIN | COMP-05 (Type System) | §3.5 COMP-05: Type System | IF-05 (Type System Interface) | TC-U-034, TC-U-035 | TC-I-004 | TC-S-001, TC-S-004 | TC-SND-008 | RISK-08 |
| REQ-FUNC-035 | Evidence per verdict: derivation trace (PASS), violation identification (FAIL), knowledge gap specification (UNCERTAIN) | COMP-05 (Type System) | §3.5 COMP-05: Type System | IF-05 (Type System Interface) | TC-U-030, TC-U-031 | TC-I-004 | TC-S-001 | TC-SND-001, TC-SND-002, TC-SND-003, TC-SND-004, TC-SND-005, TC-SND-006, TC-SND-007 | — |

### 2.5 Optimizer (REQ-FUNC-040 through REQ-FUNC-045)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-040 | CSP solver: find synonymous codon assignments satisfying all hard constraints; maximize CAI | COMP-06 (Optimizer) | §3.6 COMP-06: Optimizer | IF-06 (Optimizer Interface) | TC-U-040 | TC-I-005, TC-I-006 | TC-S-001, TC-S-004 | — | RISK-02 |
| REQ-FUNC-041 | Decision variables: one per codon position; domain = synonymous codons (1–6 per amino acid) | COMP-06 (Optimizer) | §3.6 COMP-06: Optimizer | IF-06 (Optimizer Interface) | TC-U-040 | TC-I-005 | TC-S-001 | — | — |
| REQ-FUNC-042 | Hard constraints: splicing correctness, CAI threshold, GC range, no restriction sites, no instability motifs, reading frame | COMP-06 (Optimizer) | §3.6 COMP-06: Optimizer | IF-06 (Optimizer Interface) | TC-U-042 | TC-I-005, TC-I-006 | TC-S-001, TC-S-004 | — | — |
| REQ-FUNC-043 | Feasible result: return complete codon assignment, objective value, verification record; deterministic selection | COMP-06 (Optimizer) | §3.6 COMP-06: Optimizer | IF-06 (Optimizer Interface) | TC-U-043 | TC-I-006 | TC-S-001 | — | — |
| REQ-FUNC-044 | Infeasible result: report INFEASIBLE with minimal unsatisfiable subset (MUS); MUS is minimal and verified | COMP-06 (Optimizer) | §3.6 COMP-06: Optimizer | IF-06 (Optimizer Interface) | TC-U-041, TC-U-044 | TC-I-005 | TC-S-004 | — | — |
| REQ-FUNC-045 | Optimizer determinism: same inputs → same assignment or same MUS; no randomized search | COMP-06 (Optimizer) | §3.6 COMP-06: Optimizer | IF-06 (Optimizer Interface) | TC-U-045 | TC-I-005 | TC-S-001 | — | — |

### 2.6 FFI (REQ-FUNC-050 through REQ-FUNC-053)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-050 | Folding FFI: invoke AlphaFold/ColabFold/RoseTTAFold; parse into IR-Structure; validate invariants | COMP-04 (FFI Manager) | §3.4 COMP-04: FFI Manager | IF-04 (FFI Interface) | TC-U-051 | TC-I-003 | TC-S-001 | — | RISK-05 |
| REQ-FUNC-051 | PTM FFI: invoke NetPhos/PhosphoSitePlus/dbPTM/MusiteDeep; parse into IR-Peptide PTM SLOTs | COMP-04 (FFI Manager) | §3.4 COMP-04: FFI Manager | IF-04 (FFI Interface) | TC-U-051 | TC-I-003 | TC-S-001 | — | RISK-05 |
| REQ-FUNC-052 | No internal modeling of external tools; FFI guarantees: correct input, correct output parsing, provenance | COMP-04 (FFI Manager) | §3.4 COMP-04: FFI Manager | IF-04 (FFI Interface) | TC-U-050, TC-U-052, TC-U-053 | TC-I-003 | TC-S-001 | — | — |
| REQ-FUNC-053 | FFI non-determinism: core pipeline does not depend on FFI for correctness; FFI output labeled as non-deterministic | COMP-04 (FFI Manager) | §3.4 COMP-04: FFI Manager | IF-04 (FFI Interface) | TC-U-052, TC-U-053 | TC-I-003 | TC-S-001 | — | — |

### 2.7 Certificate (REQ-FUNC-060 through REQ-FUNC-062)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-060 | Guarantee certificate generation in JSON: sequence, type verdicts + derivation, CSP record, provenance | COMP-07 (Certificate Generator) | §3.7 COMP-07: Certificate Generator | IF-07 (Certificate Interface) | TC-U-060 | TC-I-006, TC-I-007 | TC-S-001 | TC-SND-008 | — |
| REQ-FUNC-061 | Independent certificate verification: standalone checker re-derives verdicts without BioCompiler pipeline | COMP-07 (Certificate Generator) | §3.7 COMP-07: Certificate Generator | IF-07 (Certificate Interface) | TC-U-061 | TC-I-007 | TC-S-001 | — | — |
| REQ-FUNC-062 | Circuit certificate: individual gene certificates + composition checks with evidence | COMP-07 (Certificate Generator) | §3.7 COMP-07: Certificate Generator | IF-07 (Certificate Interface) | TC-U-062 | TC-I-008 | TC-S-002 | — | — |

### 2.8 Compositional Verifier (REQ-FUNC-070 through REQ-FUNC-073)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-070 | Circuit-level specification: multiple genes, promoters, terminators, topology, cellular context | COMP-08 (Compositional Verifier) | §3.8 COMP-08: Compositional Verifier | IF-08 (Composition Interface) | TC-U-070 | TC-I-008 | TC-S-002 | — | RISK-02 |
| REQ-FUNC-071 | Linker pass: promoter conflict, resource competition, splicing interference, RNA-RNA interaction | COMP-08 (Compositional Verifier) | §3.8 COMP-08: Compositional Verifier | IF-08 (Composition Interface) | TC-U-070, TC-U-071, TC-U-072, TC-U-073 | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-072 | Three-valued verdicts for composition checks with evidence | COMP-08 (Compositional Verifier) | §3.8 COMP-08: Compositional Verifier | IF-08 (Composition Interface) | TC-U-070, TC-U-071, TC-U-072, TC-U-073 | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-073 | Composition via three-valued logic: AND of individual gene verdicts and composition check verdicts | COMP-08 (Compositional Verifier) | §3.8 COMP-08: Compositional Verifier | IF-08 (Composition Interface) | TC-U-070, TC-U-071, TC-U-072, TC-U-073 | TC-I-008 | TC-S-002 | — | RISK-08 |

### 2.9 Mutation Explorer (REQ-FUNC-080 through REQ-FUNC-083)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-080 | Mutation space decomposition: intra-exonic, splice site, regulatory element categories | COMP-09 (Mutation Explorer) | §3.9 COMP-09: Mutation Explorer | IF-09 (Mutation Interface) | TC-U-080 | TC-I-007 | TC-S-002 | — | — |
| REQ-FUNC-081 | Legal multi-mutation enumeration: grammar-guided; every combination produces ≥ 1 valid isoform | COMP-09 (Mutation Explorer) | §3.9 COMP-09: Mutation Explorer | IF-09 (Mutation Interface) | TC-U-081 | TC-I-007 | TC-S-002 | — | — |
| REQ-FUNC-082 | Independence exploitation: mutations in non-overlapping exons enumerated separately; factorization | COMP-09 (Mutation Explorer) | §3.9 COMP-09: Mutation Explorer | IF-09 (Mutation Interface) | TC-U-082 | TC-I-007 | TC-S-002 | — | — |
| REQ-FUNC-083 | Constraint conflict detection: individually legal but jointly illegal mutations; epistatic interactions | COMP-09 (Mutation Explorer) | §3.9 COMP-09: Mutation Explorer | IF-09 (Mutation Interface) | TC-U-083 | TC-I-007 | TC-S-002 | — | — |

### 2.10 ORF Analyzer (REQ-FUNC-090 through REQ-FUNC-093)

| Requirement | Summary | Component (DOC-02) | Design Section (DOC-03) | Interface (DOC-04) | Unit Test(s) | Integration Test(s) | System Test(s) | Soundness Test(s) | Risk(s) |
|---|---|---|---|---|---|---|---|---|---|
| REQ-FUNC-090 | Multiple reading frame acceptance; separate translation FST per frame; per-frame IR-Peptide | COMP-10 (ORF Analyzer) | §3.10 COMP-10: ORF Analyzer | IF-10 (ORF Interface) | TC-U-090 | TC-I-009 | TC-S-003 | — | — |
| REQ-FUNC-091 | Shared constraint set: positions where mutations affect ≥ 2 reading frames | COMP-10 (ORF Analyzer) | §3.10 COMP-10: ORF Analyzer | IF-10 (ORF Interface) | TC-U-091 | TC-I-009 | TC-S-003 | — | — |
| REQ-FUNC-092 | Coupling classification: high-coupling (≥ 2 frames) vs. low-coupling (1 frame) positions | COMP-10 (ORF Analyzer) | §3.10 COMP-10: ORF Analyzer | IF-10 (ORF Interface) | TC-U-092 | TC-I-009 | TC-S-003 | — | — |
| REQ-FUNC-093 | Inter-frame constraint conflict detection: positions where optimization targets conflict; minimal conflict set | COMP-10 (ORF Analyzer) | §3.10 COMP-10: ORF Analyzer | IF-10 (ORF Interface) | TC-U-093 | TC-I-009 | TC-S-003 | — | — |

---

## 3. Non-Functional Requirements Traceability

### 3.1 Performance (REQ-NFR-001 through REQ-NFR-006)

| Requirement | Summary | Verification Method | Test Case(s) | Risk(s) |
|---|---|---|---|---|
| REQ-NFR-001 | Scanner: ≤ 1 s for 10 kb sequence on single CPU core | Performance benchmark | TC-U-007 (timing), TC-I-001 (timing) | — |
| REQ-NFR-002 | Splicing NDFST: ≤ 5 s for average 10-exon human gene on single CPU core | Performance benchmark | TC-U-010, TC-U-011 (timing), TC-I-001 (timing) | RISK-01 |
| REQ-NFR-003 | Translation FST: ≤ 100 ms for 5 kb spliced mRNA on single CPU core | Performance benchmark | TC-U-020 (timing), TC-I-002 (timing) | — |
| REQ-NFR-004 | Type checker: ≤ 10 s for all seven predicates on single mRNA on single CPU core | Performance benchmark | TC-U-030 (timing), TC-I-004 (timing) | — |
| REQ-NFR-005 | CSP optimizer: ≤ 60 s for ≤ 1,000 aa protein on single CPU core | Performance benchmark | TC-U-040 (timing), TC-I-005 (timing) | RISK-03 |
| REQ-NFR-006 | Compositional verifier: ≤ 5 min for ≤ 10 gene circuit on single CPU core | Performance benchmark | TC-U-070 (timing), TC-I-008 (timing) | — |

### 3.2 Reliability (REQ-NFR-010 through REQ-NFR-013)

| Requirement | Summary | Verification Method | Test Case(s) | Risk(s) |
|---|---|---|---|---|
| REQ-NFR-010 | Determinism for all internal pipeline stages: byte-identical output for identical input | Reproducibility testing (100 runs, cross-platform) | TC-REP-001, TC-REP-002, TC-REP-003, TC-REP-004 | RISK-07 |
| REQ-NFR-011 | Type checker soundness: no false PASS verdicts; PASS implies property holds | Adversarial testing with known violations | TC-SND-001, TC-SND-002, TC-SND-003, TC-SND-004, TC-SND-005, TC-SND-006, TC-SND-007, TC-SND-008 | RISK-07 |
| REQ-NFR-012 | CSP optimizer completeness: if feasible assignment exists, solver finds one | Constructive proof with known feasible solutions | TC-U-043 | — |
| REQ-NFR-013 | CSP optimizer correctness for infeasible: INFEASIBLE only if truly infeasible; MUS is minimal | Independent verification of MUS | TC-U-044 | — |

### 3.3 Usability (REQ-NFR-020 through REQ-NFR-023)

| Requirement | Summary | Verification Method | Test Case(s) | Risk(s) |
|---|---|---|---|---|
| REQ-NFR-020 | CLI with subcommands: design, verify, explore, analyze-orf, verify-circuit, check-cert; --json flag; exit codes 0/1/2/10/11 | CLI usability testing | TC-S-001, TC-S-002, TC-S-003, TC-S-004 (exit code verification) | — |
| REQ-NFR-021 | Python API with type-annotated signatures; synchronous and asynchronous modes; docstrings and examples | API usability testing | TC-I-001 through TC-I-009 (via API) | — |
| REQ-NFR-022 | Error messages include: (a) requirement violated, (b) position in input, (c) suggested remediation | Error message format verification | TC-U-003, TC-U-031 (error message content check) | — |
| REQ-NFR-023 | UNCERTAIN verdicts include: (a) information gap, (b) additional constraints/data to resolve | UNCERTAIN verdict format verification | TC-U-034 (UNCERTAIN output format check) | — |

### 3.4 Maintainability (REQ-NFR-030 through REQ-NFR-033)

| Requirement | Summary | Verification Method | Test Case(s) | Risk(s) |
|---|---|---|---|---|
| REQ-NFR-030 | IR schemas in protocol buffers; backward-compatible schema evolution; round-trip serialization test | Schema evolution testing | IR schema round-trip test (DOC-05 §2.1) | RISK-04 |
| REQ-NFR-031 | Each component independently testable via defined input/output contract; no preceding/succeeding stage required | Unit test isolation verification | All TC-U-xxx test cases (each tests one component in isolation) | — |
| REQ-NFR-032 | Splicing grammar in declarative YAML config; not hardcoded; cell-type overrides in config | Config change testing | TC-U-014 (context parameter test), config reload test | — |
| REQ-NFR-033 | FFI adapters implement common abstract interface; new adapters registered via config, no source modification | Adapter registration testing | TC-U-050 (adapter not found), TC-U-051 (mock adapter), TC-U-053 (provenance) | RISK-05 |

### 3.5 Portability (REQ-NFR-040 through REQ-NFR-042)

| Requirement | Summary | Verification Method | Test Case(s) | Risk(s) |
|---|---|---|---|---|
| REQ-NFR-040 | Runs on Linux x86_64 and ARM64; architecture-independent core; optional SIMD with pure-Python fallback | Cross-platform testing | TC-REP-002 (Linux x86_64, Linux ARM64) | — |
| REQ-NFR-041 | Runs on macOS x86_64 and ARM64 (Apple Silicon); same codebase as Linux | Cross-platform testing | TC-REP-002 (macOS x86_64, macOS ARM64) | — |
| REQ-NFR-042 | Core pipeline operates without external tools; FFI stages optional; graceful degradation when tools unavailable | Graceful degradation testing | TC-U-050 (adapter not found → warning, pipeline continues) | — |

### 3.6 Security (REQ-NFR-050 through REQ-NFR-052)

| Requirement | Summary | Verification Method | Test Case(s) | Risk(s) |
|---|---|---|---|---|
| REQ-NFR-050 | Certificate integrity via SHA-256 hash in design_id; tamper detection via hash recomparison | Hash verification testing | TC-U-060, TC-U-061 (hash field verification) | — |
| REQ-NFR-051 | Input validation: reject non-IUPAC characters in nucleotide sequences; report position and character | Input validation testing | TC-U-003 (InvalidSequenceError with position and character) | — |
| REQ-NFR-052 | Core pipeline requires no network access; FFI cloud-based stages optional; air-gapped operation supported | Network isolation testing | Network isolation test (core pipeline with network disabled) | — |

---

## 4. Constraint Traceability

### 4.1 Design Constraints (REQ-CON-001 through REQ-CON-004)

| Constraint | Summary | Enforced By | Verified By | Risk(s) |
|---|---|---|---|---|
| REQ-CON-001 | No probabilistic models for any internal pipeline stage; all internal reasoning is deterministic | Code review: no `random`, `numpy.random`, or probabilistic imports in internal components; linter rule | Code review; linter rule enforcement; reproducibility tests (TC-REP-001 through TC-REP-004) | — |
| REQ-CON-002 | No grammar induction; all grammar rules from curated biological knowledge (GENCODE, literature) | Splicing grammar loaded from declarative YAML config, not learned from data | Config-driven grammar rules test; grammar source audit (all rules traceable to GENCODE/literature citations) | — |
| REQ-CON-003 | No internal folding or PTM models; these handled exclusively via FFI (COMP-04) | Architecture constraint: no folding/PTM code in COMP-01 through COMP-03, COMP-05 through COMP-10 | FFI-only folding/PTM test; code review: no physics/ML models in internal components | — |
| REQ-CON-004 | No claim that biology implements compilation; compiler metaphor is a design pattern only | Documentation review: no theoretical claims about biological compilation in any artifact | Documentation review; claim audit across all DOC-xx artifacts | — |

### 4.2 Environmental Constraints (REQ-CON-010, REQ-CON-011)

| Constraint | Summary | Enforced By | Verified By | Risk(s) |
|---|---|---|---|---|
| REQ-CON-010 | Memory limits: ≤ 32 GB RAM for single-gene analysis; ≤ 64 GB RAM for circuit-level (10 genes) | Memory profiling in CI: pipeline run with memory monitoring; assertion on peak RSS | Memory profiling test with realistic inputs (10 kb gene; 10-gene circuit) | RISK-03 |
| REQ-CON-011 | GPU required only for FFI stages; core pipeline (COMP-01 through COMP-03, COMP-05 through COMP-10) is CPU-only | Architecture: no GPU-dependent code in core components; pure-Python/CPU fallback for all stages | CPU-only pipeline test: run full core pipeline on GPU-less machine; assert no GPU import | — |

---

## 5. Architectural Decision Traceability

| Decision | Requirements Addressed | Components Affected | Rationale (DOC-06) | Risks |
|---|---|---|---|---|
| AD-01 | REQ-FUNC-001–004, REQ-FUNC-010–015, REQ-FUNC-020–023, REQ-NFR-031 | COMP-01 through COMP-10, IR Bus | Pipeline architecture (not monolith): mirrors LLVM's proven multi-pass design; each stage independently testable and replaceable; enables separation of formalizable from non-formalizable stages (DOC-06 §3.1) | — |
| AD-02 | REQ-NFR-030, REQ-NFR-031, REQ-CON-001 | IR Bus, all components | Protocol Buffers for IR schemas: language-neutral, schema-evolution-friendly, efficient binary serialization; supports cross-toolchain IR interchange (DOC-06 §3.3) | RISK-04 |
| AD-03 | REQ-FUNC-010, REQ-FUNC-011, REQ-FUNC-012, REQ-FUNC-014, REQ-CON-001, REQ-CON-002 | COMP-02 (Splicing Engine) | NDFSTs for splicing: naturally models non-deterministic alternative splicing via set-valued output; computation is deterministic while biology is non-deterministic (DOC-06 §3.2, Flaw 5) | RISK-01 |
| AD-04 | REQ-FUNC-040, REQ-FUNC-041, REQ-FUNC-042, REQ-FUNC-043, REQ-FUNC-044, REQ-NFR-012, REQ-NFR-013 | COMP-06 (Optimizer) | Constraint satisfaction (not optimization): gene design is a feasibility problem; CSP + MUS provides provably complete explanations of failure (DOC-06 Flaw 8) | RISK-02 |
| AD-05 | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-034, REQ-FUNC-035, REQ-NFR-011 | COMP-05 (Type System), COMP-08 (Compositional Verifier) | Three-valued logic for verdicts: PASS/FAIL/UNCERTAIN captures epistemic boundary where FFI oracles cannot provide definitive answers (DOC-06 §4.1, Flaw 5) | RISK-08 |
| AD-06 | REQ-FUNC-050, REQ-FUNC-051, REQ-FUNC-052, REQ-FUNC-053, REQ-CON-003, REQ-NFR-042 | COMP-04 (FFI Manager) | FFI for folding and PTMs: these require physics-based or ML oracles; FFI keeps core deterministic (DOC-06 Flaw 1, Flaw 2) | RISK-05 |
| AD-07 | REQ-FUNC-011, REQ-FUNC-013, REQ-NFR-032, REQ-CON-002 | COMP-02 (Splicing Engine), config/ module | Declarative grammar configuration: organism-specific rules loaded from YAML, not hardcoded; enables updates without source code modification (DOC-06 Flaw 3, Flaw 4) | — |

---

## 6. IR Invariant Traceability

### 6.1 IR-Seq Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-SEQ-01 | IR-Seq / COMP-01 (Scanner) | Length in nucleotides is a positive integer; coding region length is a positive multiple of 3 | TC-U-001, TC-U-002 (sequence length verification) |
| INV-SEQ-02 | IR-Seq / COMP-02 (Splicing Engine) | Every splice-site annotation has a matching donor-acceptor pair within the same contig | TC-U-012 (donor-acceptor pair verification) |
| INV-SEQ-03 | IR-Seq / COMP-02 (Splicing Engine) | Isoform IDs are unique within a pipeline run | TC-U-010, TC-U-011 (isoform ID uniqueness) |
| INV-SEQ-04 | IR-Seq / COMP-01 (Scanner) | GC content is within [0.0, 1.0] and equals the computed GC fraction of the sequence | TC-U-001 (GC content verification) |

### 6.2 IR-Peptide Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-PEP-01 | IR-Peptide / COMP-03 (Translation Engine) | Peptide length equals (coding region length / 3) − 1 (stop codon excluded) | TC-U-020 (peptide length verification) |
| INV-PEP-02 | IR-Peptide / COMP-03 (Translation Engine) | Every amino acid residue has a back-reference to its source codon in IR-Seq; codon_assignments count equals peptide_length | TC-U-020 (codon assignment count verification) |
| INV-PEP-03 | IR-Peptide / COMP-03 (Translation Engine) | No duplicate isoform references within the same peptide record; source_isoform_id is singular and non-empty | TC-I-002 (isoform reference uniqueness) |

### 6.3 IR-Structure Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-STR-01 | IR-Structure / COMP-04 (FFI Manager) | Every structural element has a confidence score in [0.0, 1.0] (normalized pLDDT/100); mean_plddt in [0.0, 100.0] | TC-U-052 (output validation rejects pLDDT outside range) |
| INV-STR-02 | IR-Structure / COMP-04 (FFI Manager) | PTM site predictions reference valid residue positions in the parent IR-Peptide [0, peptide_length) | TC-U-052 (output validation rejects invalid PTM positions) |

### 6.4 IR-Circuit Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-CIR-01 | IR-Circuit / COMP-08 (Compositional Verifier) | Circuit graph is acyclic (no circular promoter dependencies) | TC-U-070 (topological sort verification) |
| INV-CIR-02 | IR-Circuit / COMP-08 (Compositional Verifier) | Every gene node references a valid, non-expired certificate from COMP-07 | TC-I-008 (certificate reference verification) |

### 6.5 Scanner Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-SCAN-01 | COMP-01 (Scanner) | Tokens are ordered by position (ascending) | TC-U-004, TC-U-005, TC-U-006 (token ordering verification) |
| INV-SCAN-02 | COMP-01 (Scanner) | No token's position range exceeds sequence bounds | TC-U-001, TC-U-002 (position bounds check) |
| INV-SCAN-03 | COMP-01 (Scanner) | Each position is examined by all active DFAs | TC-U-006 (all motifs detected at all positions) |
| INV-SCAN-04 | COMP-01 (Scanner) | Determinism: identical input produces identical token list | TC-U-007 (100-run byte-identical comparison) |

### 6.6 Splicing Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-SPL-01 | COMP-02 (Splicing Engine) | Isoform set is complete: all valid parse paths are represented | TC-U-010, TC-U-011 (completeness against known isoform sets) |
| INV-SPL-02 | COMP-02 (Splicing Engine) | Isoform set is sound: every isoform satisfies the splicing grammar | TC-U-012, TC-U-013 (each isoform verified against grammar) |
| INV-SPL-03 | COMP-02 (Splicing Engine) | Computation is deterministic: same input + context → same isoform set | TC-U-015 (100-run comparison) |
| INV-SPL-04 | COMP-02 (Splicing Engine) | Each isoform has a unique parse path (no duplicate isoforms) | TC-U-010, TC-U-011 (isoform uniqueness within set) |

### 6.7 Translation Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-TRA-01 | COMP-03 (Translation Engine) | Amino acid sequence length equals number of sense codons processed | TC-U-020 (length verification) |
| INV-TRA-02 | COMP-03 (Translation Engine) | Every codon assignment maps the specified codon to the correct amino acid per the genetic code | TC-U-020 (codon-to-amino-acid mapping verification) |
| INV-TRA-03 | COMP-03 (Translation Engine) | Determinism: same input → same output (modulo frameshift warning annotations) | TC-U-022 (100-run comparison) |

### 6.8 FFI Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-FFI-01 | COMP-04 (FFI Manager) | Every SLOT field is either empty or filled with valid data (never partially filled) | TC-U-051, TC-U-052 (SLOT field state after FFI invocation) |
| INV-FFI-02 | COMP-04 (FFI Manager) | Provenance metadata is recorded for every FFI invocation | TC-U-053 (provenance field verification) |
| INV-FFI-03 | COMP-04 (FFI Manager) | FFI output invariants (INV-STR-01, INV-STR-02, INV-PEP-03) are validated before acceptance | TC-U-052 (output validation test) |

### 6.9 Type System Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-TYP-01 | COMP-05 (Type System) | Soundness: if PASS is returned, the property holds for the sequence under the specified conditions | TC-SND-001 through TC-SND-008 (adversarial testing) |
| INV-TYP-02 | COMP-05 (Type System) | Determinism: same input → same verdicts | TC-U-032 (100-run comparison) |
| INV-TYP-03 | COMP-05 (Type System) | Three-valued composition preserves soundness | TC-U-034, TC-U-035 (composition truth table verification) |
| INV-TYP-04 | COMP-05 (Type System) | Subtyping: if SpliceCorrect(C1) passes and C1 <: C2, then SpliceCorrect(C2) passes | TC-U-033 (subtyping verification) |

### 6.10 CSP Invariants

| Invariant | IR Level / Component | Enforcement | Test(s) |
|---|---|---|---|
| INV-CSP-01 | COMP-06 (Optimizer) | Solution satisfies ALL hard constraints (verified by re-running type checker) | TC-U-040, TC-U-042 (constraint satisfaction verification) |
| INV-CSP-02 | COMP-06 (Optimizer) | MUS is truly unsatisfiable (verified independently) | TC-U-044 (MUS unsatisfiability verification) |
| INV-CSP-03 | COMP-06 (Optimizer) | Determinism: same input → same solution or same INFEASIBLE report | TC-U-045 (100-run comparison) |
| INV-CSP-04 | COMP-06 (Optimizer) | Completeness: if a feasible assignment exists, solver finds one | TC-U-043 (feasible solution discovery) |

---

## 7. Test Coverage Summary

### 7.1 Test Case Inventory by Level

| Test Level | Total Test Cases | Requirements Covered | Coverage % |
|---|---|---|---|
| Unit Tests (TC-U-xxx) | 30 | 35 / 35 functional requirements | 100% |
| Integration Tests (TC-I-xxx) | 9 | 30 / 35 functional requirements | 86% |
| System Tests (TC-S-xxx) | 4 | 26 / 35 functional requirements | 74% |
| Soundness Tests (TC-SND-xxx) | 8 | 4 / 35 functional requirements (focused on REQ-FUNC-030, 031, 034, 035) | 11% (by design: adversarial focus on type system) |
| Reproducibility Tests (TC-REP-xxx) | 4 | 6 / 35 functional requirements (REQ-FUNC-004, 014, 022, 032, 045 + REQ-NFR-010) | 17% (by design: determinism focus) |
| Biological Validation (TC-V-xxx) | 5 | 8 / 35 functional requirements | 23% |
| Comparison Tests (TC-C-xxx) | 3 | 3 / 35 functional requirements | 9% |
| **Total Unique Test Cases** | **63** | **35 / 35 functional requirements** | **100%** |

### 7.2 Coverage by Requirement Category

| Category | Requirements | With Unit Test | With Integration Test | With System Test | With Soundness Test | With Any Test |
|---|---|---|---|---|---|---|
| Scanner (REQ-FUNC-001–004) | 4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 |
| Splicing (REQ-FUNC-010–015) | 6 | 6/6 | 6/6 | 6/6 | 1/6 | 6/6 |
| Translation (REQ-FUNC-020–023) | 4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 |
| Type System (REQ-FUNC-030–035) | 6 | 6/6 | 6/6 | 6/6 | 4/6 | 6/6 |
| Optimizer (REQ-FUNC-040–045) | 6 | 6/6 | 5/6 | 5/6 | 0/6 | 6/6 |
| FFI (REQ-FUNC-050–053) | 4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 |
| Certificate (REQ-FUNC-060–062) | 3 | 3/3 | 3/3 | 3/3 | 1/3 | 3/3 |
| Compositional (REQ-FUNC-070–073) | 4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 |
| Mutation (REQ-FUNC-080–083) | 4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 |
| ORF (REQ-FUNC-090–093) | 4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 |
| **Total** | **45** | **45/45** | **44/45** | **44/45** | **6/45** | **45/45** |

### 7.3 NFR Coverage Summary

| NFR Category | Requirements | With Verification Method | With Test Case | Coverage % |
|---|---|---|---|---|
| Performance (REQ-NFR-001–006) | 6 | 6/6 | 6/6 | 100% |
| Reliability (REQ-NFR-010–013) | 4 | 4/4 | 4/4 | 100% |
| Usability (REQ-NFR-020–023) | 4 | 4/4 | 4/4 | 100% |
| Maintainability (REQ-NFR-030–033) | 4 | 4/4 | 4/4 | 100% |
| Portability (REQ-NFR-040–042) | 3 | 3/3 | 3/3 | 100% |
| Security (REQ-NFR-050–052) | 3 | 3/3 | 3/3 | 100% |
| **Total** | **24** | **24/24** | **24/24** | **100%** |

---

## 8. Backward Traceability Index

This section provides the reverse mapping: from test cases, interfaces, and components back to the requirements they verify or satisfy. This enables impact analysis (e.g., "if test TC-U-031 fails, which requirements are affected?") and coverage auditing (e.g., "does interface IF-05 cover all type system requirements?").

### 8.1 Test Cases → Requirements

#### Unit Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-U-001 | REQ-FUNC-001 |
| TC-U-002 | REQ-FUNC-001 |
| TC-U-003 | REQ-FUNC-001, REQ-FUNC-003, REQ-NFR-051 |
| TC-U-004 | REQ-FUNC-002 |
| TC-U-005 | REQ-FUNC-002 |
| TC-U-006 | REQ-FUNC-002 |
| TC-U-007 | REQ-FUNC-004, REQ-NFR-001 |
| TC-U-010 | REQ-FUNC-010, REQ-FUNC-015 |
| TC-U-011 | REQ-FUNC-012 |
| TC-U-012 | REQ-FUNC-011 |
| TC-U-013 | REQ-FUNC-011 |
| TC-U-014 | REQ-FUNC-013 |
| TC-U-015 | REQ-FUNC-014 |
| TC-U-020 | REQ-FUNC-020, REQ-FUNC-023 |
| TC-U-021 | REQ-FUNC-021 |
| TC-U-022 | REQ-FUNC-022 |
| TC-U-023 | REQ-FUNC-021 |
| TC-U-030 | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-035 |
| TC-U-031 | REQ-FUNC-031, REQ-FUNC-035 |
| TC-U-032 | REQ-FUNC-032 |
| TC-U-033 | REQ-FUNC-033 |
| TC-U-034 | REQ-FUNC-034 |
| TC-U-035 | REQ-FUNC-034 |
| TC-U-040 | REQ-FUNC-040, REQ-FUNC-041, REQ-FUNC-043 |
| TC-U-041 | REQ-FUNC-044 |
| TC-U-042 | REQ-FUNC-042 |
| TC-U-043 | REQ-FUNC-043, REQ-NFR-012 |
| TC-U-044 | REQ-FUNC-044, REQ-NFR-013 |
| TC-U-045 | REQ-FUNC-045 |
| TC-U-050 | REQ-FUNC-052 |
| TC-U-051 | REQ-FUNC-050, REQ-FUNC-051 |
| TC-U-052 | REQ-FUNC-052, REQ-FUNC-053 |
| TC-U-053 | REQ-FUNC-052 |
| TC-U-060 | REQ-FUNC-060 |
| TC-U-061 | REQ-FUNC-061 |
| TC-U-062 | REQ-FUNC-062 |
| TC-U-070 | REQ-FUNC-071 |
| TC-U-071 | REQ-FUNC-071 |
| TC-U-072 | REQ-FUNC-071 |
| TC-U-073 | REQ-FUNC-071 |
| TC-U-080 | REQ-FUNC-080 |
| TC-U-081 | REQ-FUNC-081 |
| TC-U-082 | REQ-FUNC-082 |
| TC-U-083 | REQ-FUNC-083 |
| TC-U-090 | REQ-FUNC-090 |
| TC-U-091 | REQ-FUNC-091 |
| TC-U-092 | REQ-FUNC-092 |
| TC-U-093 | REQ-FUNC-093 |

#### Integration Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-I-001 | REQ-FUNC-001, REQ-FUNC-002, REQ-FUNC-003, REQ-FUNC-004, REQ-FUNC-010, REQ-FUNC-011, REQ-FUNC-012, REQ-FUNC-013, REQ-FUNC-014, REQ-FUNC-015 |
| TC-I-002 | REQ-FUNC-015, REQ-FUNC-020, REQ-FUNC-021, REQ-FUNC-022, REQ-FUNC-023 |
| TC-I-003 | REQ-FUNC-050, REQ-FUNC-051, REQ-FUNC-052, REQ-FUNC-053 |
| TC-I-004 | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-032, REQ-FUNC-033, REQ-FUNC-034, REQ-FUNC-035 |
| TC-I-005 | REQ-FUNC-040, REQ-FUNC-042, REQ-FUNC-044, REQ-FUNC-045 |
| TC-I-006 | REQ-FUNC-043, REQ-FUNC-060 |
| TC-I-007 | REQ-FUNC-060, REQ-FUNC-061, REQ-FUNC-062, REQ-FUNC-080, REQ-FUNC-081, REQ-FUNC-082, REQ-FUNC-083 |
| TC-I-008 | REQ-FUNC-062, REQ-FUNC-070, REQ-FUNC-071, REQ-FUNC-072, REQ-FUNC-073 |
| TC-I-009 | REQ-FUNC-090, REQ-FUNC-091, REQ-FUNC-092, REQ-FUNC-093 |

#### System Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-S-001 | REQ-FUNC-001, REQ-FUNC-002, REQ-FUNC-003, REQ-FUNC-004, REQ-FUNC-010, REQ-FUNC-011, REQ-FUNC-012, REQ-FUNC-013, REQ-FUNC-014, REQ-FUNC-015, REQ-FUNC-020, REQ-FUNC-021, REQ-FUNC-022, REQ-FUNC-023, REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-032, REQ-FUNC-033, REQ-FUNC-034, REQ-FUNC-035, REQ-FUNC-040, REQ-FUNC-042, REQ-FUNC-043, REQ-FUNC-050, REQ-FUNC-051, REQ-FUNC-060, REQ-FUNC-061 |
| TC-S-002 | REQ-FUNC-070, REQ-FUNC-071, REQ-FUNC-072, REQ-FUNC-073, REQ-FUNC-062, REQ-FUNC-080, REQ-FUNC-081, REQ-FUNC-082, REQ-FUNC-083 |
| TC-S-003 | REQ-FUNC-090, REQ-FUNC-091, REQ-FUNC-092, REQ-FUNC-093, REQ-FUNC-021 |
| TC-S-004 | REQ-FUNC-040, REQ-FUNC-042, REQ-FUNC-044 |

#### Soundness Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-SND-001 | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-002 | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-003 | REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-004 | REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-005 | REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-006 | REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-007 | REQ-FUNC-031, REQ-FUNC-035, REQ-NFR-011 |
| TC-SND-008 | REQ-FUNC-034, REQ-FUNC-035, REQ-FUNC-060, REQ-NFR-011 |

#### Reproducibility Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-REP-001 | REQ-FUNC-004, REQ-FUNC-014, REQ-FUNC-022, REQ-FUNC-032, REQ-FUNC-045, REQ-NFR-010 |
| TC-REP-002 | REQ-FUNC-004, REQ-FUNC-014, REQ-NFR-010, REQ-NFR-040, REQ-NFR-041 |
| TC-REP-003 | REQ-FUNC-045, REQ-NFR-010 |
| TC-REP-004 | REQ-FUNC-061, REQ-NFR-010 |

#### Biological Validation Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-V-001 | REQ-FUNC-010, REQ-FUNC-011 |
| TC-V-002 | REQ-FUNC-020 |
| TC-V-003 | REQ-FUNC-040, REQ-FUNC-042 |
| TC-V-004 | REQ-FUNC-002, REQ-FUNC-042 |
| TC-V-005 | REQ-FUNC-090, REQ-FUNC-091 |

#### Comparison Tests → Requirements

| Test Case | Requirement(s) |
|---|---|
| TC-C-001 | REQ-FUNC-040, REQ-FUNC-042 |
| TC-C-002 | REQ-FUNC-042 |
| TC-C-003 | REQ-FUNC-040, REQ-FUNC-042 |

### 8.2 Interfaces → Requirements

| Interface | Requirement(s) |
|---|---|
| IF-01 (Scanner) | REQ-FUNC-001, REQ-FUNC-002, REQ-FUNC-003, REQ-FUNC-004 |
| IF-02 (Splicing Engine) | REQ-FUNC-010, REQ-FUNC-011, REQ-FUNC-012, REQ-FUNC-013, REQ-FUNC-014, REQ-FUNC-015 |
| IF-03 (Translation Engine) | REQ-FUNC-020, REQ-FUNC-021, REQ-FUNC-022, REQ-FUNC-023 |
| IF-04 (FFI Manager) | REQ-FUNC-050, REQ-FUNC-051, REQ-FUNC-052, REQ-FUNC-053 |
| IF-05 (Type System) | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-032, REQ-FUNC-033, REQ-FUNC-034, REQ-FUNC-035 |
| IF-06 (Optimizer) | REQ-FUNC-040, REQ-FUNC-041, REQ-FUNC-042, REQ-FUNC-043, REQ-FUNC-044, REQ-FUNC-045 |
| IF-07 (Certificate Generator) | REQ-FUNC-060, REQ-FUNC-061, REQ-FUNC-062 |
| IF-08 (Compositional Verifier) | REQ-FUNC-070, REQ-FUNC-071, REQ-FUNC-072, REQ-FUNC-073 |
| IF-09 (Mutation Explorer) | REQ-FUNC-080, REQ-FUNC-081, REQ-FUNC-082, REQ-FUNC-083 |
| IF-10 (ORF Analyzer) | REQ-FUNC-090, REQ-FUNC-091, REQ-FUNC-092, REQ-FUNC-093 |
| IF-DATA-01 (IR-Seq Schema) | REQ-FUNC-001, REQ-FUNC-002, REQ-FUNC-003, REQ-FUNC-004, REQ-FUNC-010, REQ-FUNC-011, REQ-FUNC-012, REQ-FUNC-013, REQ-FUNC-014, REQ-FUNC-015 |
| IF-DATA-02 (IR-Peptide Schema) | REQ-FUNC-020, REQ-FUNC-021, REQ-FUNC-022, REQ-FUNC-023, REQ-FUNC-050, REQ-FUNC-051 |
| IF-DATA-03 (IR-Structure Schema) | REQ-FUNC-050, REQ-FUNC-052 |
| IF-DATA-04 (IR-Circuit Schema) | REQ-FUNC-070, REQ-FUNC-071, REQ-FUNC-072, REQ-FUNC-073 |

### 8.3 Components → Requirements

| Component | Requirement(s) |
|---|---|
| COMP-01 (Scanner) | REQ-FUNC-001, REQ-FUNC-002, REQ-FUNC-003, REQ-FUNC-004 |
| COMP-02 (Splicing Engine) | REQ-FUNC-010, REQ-FUNC-011, REQ-FUNC-012, REQ-FUNC-013, REQ-FUNC-014, REQ-FUNC-015 |
| COMP-03 (Translation Engine) | REQ-FUNC-020, REQ-FUNC-021, REQ-FUNC-022, REQ-FUNC-023 |
| COMP-04 (FFI Manager) | REQ-FUNC-050, REQ-FUNC-051, REQ-FUNC-052, REQ-FUNC-053 |
| COMP-05 (Type System) | REQ-FUNC-030, REQ-FUNC-031, REQ-FUNC-032, REQ-FUNC-033, REQ-FUNC-034, REQ-FUNC-035 |
| COMP-06 (Optimizer) | REQ-FUNC-040, REQ-FUNC-041, REQ-FUNC-042, REQ-FUNC-043, REQ-FUNC-044, REQ-FUNC-045 |
| COMP-07 (Certificate Generator) | REQ-FUNC-060, REQ-FUNC-061, REQ-FUNC-062 |
| COMP-08 (Compositional Verifier) | REQ-FUNC-070, REQ-FUNC-071, REQ-FUNC-072, REQ-FUNC-073 |
| COMP-09 (Mutation Explorer) | REQ-FUNC-080, REQ-FUNC-081, REQ-FUNC-082, REQ-FUNC-083 |
| COMP-10 (ORF Analyzer) | REQ-FUNC-090, REQ-FUNC-091, REQ-FUNC-092, REQ-FUNC-093 |

### 8.4 Risks → Affected Requirements and Artifacts

| Risk | Affected Requirements | Affected Component(s) | Affected Test(s) |
|---|---|---|---|
| RISK-01 | REQ-FUNC-010, REQ-FUNC-011, REQ-NFR-002 | COMP-02 | TC-U-010, TC-U-012, TC-U-013, TC-I-001 |
| RISK-02 | REQ-FUNC-013, REQ-FUNC-040, REQ-FUNC-070 | COMP-02, COMP-06, COMP-08 | TC-U-014, TC-U-040, TC-U-070 |
| RISK-03 | REQ-FUNC-030, REQ-NFR-005, REQ-CON-010 | COMP-05, COMP-06 | TC-U-030, TC-U-040, TC-SND-001–008 |
| RISK-04 | REQ-FUNC-001, REQ-NFR-030 | COMP-01, IR Bus | TC-U-001, TC-U-002, TC-U-003 |
| RISK-05 | REQ-FUNC-050, REQ-FUNC-051, REQ-NFR-033 | COMP-04 | TC-U-050, TC-U-051, TC-I-003 |
| RISK-07 | REQ-NFR-010, REQ-NFR-011 | COMP-05, all internal components | TC-REP-001–004, TC-SND-001–008 |
| RISK-08 | REQ-FUNC-034, REQ-FUNC-073 | COMP-05, COMP-08 | TC-U-034, TC-U-035, TC-U-070–073 |

---

*End of DOC-08: Traceability Matrix*
