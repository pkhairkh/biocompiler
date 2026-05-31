# DOC-07: Project Plan (PP)

| Field | Value |
|---|---|
| **Document ID** | DOC-07 |
| **Version** | 1.0.0-draft |
| **Status** | ROUGH DRAFT |
| **Date** | 2026-05-30 |
| **Prepared By** | BioCompiler Project Team |
| **Reviewed By** | [TBD at baseline review] |
| **Approved By** | [TBD at baseline review] |
| **Standard** | IEEE 16326-2019 / PMBOK 7th Edition |

---

## 1. Project Overview

### 1.1 Objective — Build BioCompiler v1.0

The objective of the BioCompiler project is to design, implement, verify, and deliver **BioCompiler v1.0**, a compiler-framework software system that applies compiler-engineering design patterns — staged transformation, typed intermediate representations, composable passes, and type-based verification — to the formalizable stages of gene-to-protein processing, specifically splicing and translation. The system SHALL provide deterministic, machine-checkable guarantees for designed mRNA sequences without relying on probabilistic models.

BioCompiler v1.0 SHALL deliver the following capabilities at release:

- **Splicing grammar engine** that parses pre-mRNA sequences against a parameterized non-deterministic finite-state transducer (NDFST) to compute the set of all possible splice isoforms under a specified cellular context, without assigning probabilities to isoforms.
- **Translation engine** that maps codons to amino acids via a deterministic finite-state transducer (FST) implementing the standard genetic code, with support for selenocysteine insertion, pyrrolysine incorporation, and programmed ribosomal frameshift detection.
- **Type system** that checks mRNA sequences against biological correctness properties, producing three-valued verdicts (PASS / FAIL / UNCERTAIN) with derivation traces, violation identification, and knowledge-gap specification.
- **Constraint-satisfaction-based codon optimizer** that eliminates cryptic splice sites while maximizing codon adaptation index (CAI), reporting INFEASIBLE with a minimal unsatisfiable subset (MUS) when no feasible assignment exists.
- **Compositional verifier** for multi-gene circuits, checking promoter conflicts, resource competition, splicing interference, and RNA-RNA interactions, with three-valued verdicts composed via a defined algebra.
- **Certificate generator** that produces machine-checkable guarantee certificates independently verifiable by a separate checker program without access to the BioCompiler pipeline.
- **Mutation explorer** that enumerates legal multi-mutation combinations, exploits independence across exons, and detects constraint conflicts between individually legal but jointly illegal mutations.
- **ORF analyzer** that computes shared constraint sets for overlapping reading frames, classifies positions by coupling strength, and detects constraint conflicts between frames.

The project SHALL produce a system that satisfies all requirements specified in DOC-01 (SRS), verified and validated per DOC-05 (SVVP), with complete traceability documented in DOC-08 (TM).

### 1.2 Scope — Refer to DOC-01 SRS Section 1.2

The scope of the BioCompiler project is defined by the Software Requirements Specification (DOC-01), Section 1.2. The following summary references that specification; the authoritative scope definition resides in DOC-01.

**In Scope** (as specified in DOC-01 Section 1.2):

- Parsing pre-mRNA sequences against a splicing grammar to compute the set of all possible splice isoforms using an NDFST.
- Translating spliced mRNA to amino acid sequences via the standard genetic code implemented as a deterministic FST.
- Type-checking mRNA sequences against biological correctness properties, producing three-valued verdicts (PASS / FAIL / UNCERTAIN).
- Constraint-satisfaction-based codon optimization that eliminates cryptic splice sites while maximizing CAI.
- Compositional verification of multi-gene circuits with three-valued verdicts composed via a defined algebra.
- Generation of machine-checkable guarantee certificates (proof-carrying gene designs).
- Grammar-guided enumeration of legal multi-mutation combinations with conflict detection.
- Analysis of overlapping reading frames with shared constraint set computation and coupling classification.

**Out of Scope** (as specified in DOC-01 Section 1.2):

- Protein structure prediction (handled exclusively via FFI to external tools).
- Post-translational modification prediction (handled exclusively via FFI).
- Cellular dynamics simulation.
- Probabilistic predictions of any kind.
- Grammar induction from biological sequence data.
- Experimental laboratory validation.

### 1.3 Success Criteria — Refer to DOC-05 SVVP Section 4

The success criteria for the BioCompiler project are defined by the acceptance criteria in the Software Verification & Validation Plan (DOC-05), Section 4. The following summary references that specification; the authoritative success criteria reside in DOC-05.

**System-Level Acceptance** (DOC-05 Section 4.1):

| Acceptance ID | Criterion | Threshold |
|---|---|---|
| AC-01 | All unit tests pass | 100% pass rate; zero failures |
| AC-02 | Unit test coverage ≥ 90% overall; 100% for COMP-05 and COMP-06 | ≥ 90% line coverage overall; 100% for COMP-05, COMP-06 |
| AC-03 | All integration tests pass | 100% pass rate; zero failures |
| AC-04 | All system tests pass (TC-S-001 through TC-S-004) | 100% pass rate; zero failures |
| AC-05 | All soundness tests pass — zero false PASS verdicts | 100% pass rate; any false PASS blocks release |
| AC-06 | All reproducibility tests pass | 100% byte-identical on same platform; ±1 ULP across platforms |
| AC-07 | All biological validation tests pass | All thresholds met per DOC-05 Section 3.1 |

**Per-Phase Acceptance** (DOC-05 Section 4.2):

Each phase has its own acceptance gate with specific test case pass conditions. Phase transitions require passing the quality gate defined in Section 7 of this document before subsequent phase work may begin. The per-phase acceptance criteria in DOC-05 Section 4.2 are incorporated by reference and serve as the mandatory checkpoint definitions for milestones M1 through M4.

**Pilot Validation Studies** (DOC-05 Section 3.3):

Three pilot studies SHALL be completed and reviewed by domain experts as part of Phase 4 acceptance:

1. **GFP Gene Design for Mammalian Expression** — Design an optimized EGFP gene for HEK293T cells with full splicing correctness guarantees.
2. **Toggle Switch Circuit Verification** — Verify a 2-gene toggle switch circuit (LacI + TetR) for *E. coli* deployment.
3. **SARS-CoV-2 ORF1a/ORF1b Overlapping Frame Analysis** — Analyze overlapping reading frames in the SARS-CoV-2 ORF1a/ORF1b region.

---

## 2. Work Breakdown Structure (WBS)

### WBS Level 1: Phases

```
1.0 BioCompiler Project
  1.1 Phase 1: Foundation (Months 1–4)
  1.2 Phase 2: Splicing-Aware Gene Design (Months 5–9)
  1.3 Phase 3: Compositional Verification (Months 10–15)
  1.4 Phase 4: Advanced Capabilities (Months 16–21)
  1.5 Project Management (Ongoing)
```

### WBS Level 2: Work Packages per Phase

#### Phase 1: Foundation (18 person-weeks total)

Phase 1 establishes the core infrastructure: the Intermediate Representation (IR) schemas that all subsequent components depend on, the splicing grammar for the human genome, the type system prototype, and the FFI adapter for AlphaFold. The phase culminates in an end-to-end integration test for a single gene.

| WBS | Work Package | Deliverable | Effort (PW) | Dependencies | Requirements |
|---|---|---|---|---|---|
| 1.1.1 | IR Schema Definition | `ir_schemas/` with `.proto` files for IR-Seq, IR-Peptide, IR-Structure, IR-Circuit | 3 PW | none | REQ-FUNC-001–023 |
| 1.1.2 | Splicing Grammar (Human) | `splicing/` with NDFST implementation, GENCODE-derived PWMs, BDD-based isoform set computation | 6 PW | 1.1.1 | REQ-FUNC-010–015 |
| 1.1.3 | Type System Prototype | `type_system/` with three-valued verdict logic, seven type predicates, subtyping, conjunction algebra | 4 PW | 1.1.2 | REQ-FUNC-030–035 |
| 1.1.4 | AlphaFold FFI Adapter | `ffi/adapters/alphafold.py` with adapter contract, SLOT field population, provenance tracking, output validation | 3 PW | 1.1.1 | REQ-FUNC-050–053 |
| 1.1.5 | Phase 1 Integration Test | End-to-end test for one gene: FASTA input → scan → splice → translate → FFI call (mock) → output | 2 PW | 1.1.1–1.1.4 | — |

**Phase 1 Risk Focus**: RISK-01 (splicing grammar precision) and RISK-04 (IR schema design mistakes) are the primary risks during this phase. The integration test (1.1.5) serves as an early detector for schema design errors. The splicing grammar work package is the largest single effort at 6 PW, reflecting the complexity of constructing a biologically accurate NDFST from GENCODE annotations.

#### Phase 2: Splicing-Aware Gene Design (17 person-weeks total)

Phase 2 builds the core gene design capability: the CSP-based codon optimizer, splicing constraint integration into the optimizer, the certificate generator, and the GFP pilot study. The phase culminates in a complete single-gene design pipeline with formal verification and guarantee certificates.

| WBS | Work Package | Deliverable | Effort (PW) | Dependencies | Requirements |
|---|---|---|---|---|---|
| 1.2.1 | CSP Codon Optimizer | `optimizer/` with Z3-backed CSP solver, CAI objective, hard constraints (splicing, GC, restriction sites, instability motifs), MUS computation for infeasible cases | 5 PW | 1.1.1 | REQ-FUNC-040–045 |
| 1.2.2 | Splicing Constraint Integration | Splicing correctness constraints embedded in CSP solver; cryptic splice site elimination as hard constraint; cellular context parameterization | 4 PW | 1.2.1, 1.1.2 | REQ-FUNC-042 |
| 1.2.3 | Certificate Generator | `certificate/` with JSON certificate format, standalone checker, derivation trace inclusion, provenance metadata | 3 PW | 1.1.3 | REQ-FUNC-060–062 |
| 1.2.4 | GFP Pilot Study | End-to-end GFP gene design for HEK293T cells; validation against DOC-05 Pilot 1 criteria | 2 PW | 1.2.2, 1.2.3 | — |
| 1.2.5 | Phase 2 Integration & V&V | Full single-gene pipeline integration; unit, integration, and soundness testing per DOC-05 Section 4.2 Phase 2 criteria | 3 PW | 1.2.1–1.2.4 | REQ-NFR-011 |

**Phase 2 Risk Focus**: RISK-02 (type system too conservative) and RISK-03 (CSP solver too slow for long proteins) are the primary risks during this phase. The GFP pilot study (1.2.4) provides early feedback on practical performance and conservativeness. The V&V work package (1.2.5) is the gatekeeper for the soundness property (REQ-NFR-011).

#### Phase 3: Compositional Verification (14 person-weeks total)

Phase 3 extends the system from single-gene design to multi-gene circuit verification, adds the mutation explorer for single-gene mutation analysis, and integrates the ORF analyzer for overlapping reading frames. The phase culminates in the toggle switch pilot study.

| WBS | Work Package | Deliverable | Effort (PW) | Dependencies | Requirements |
|---|---|---|---|---|---|
| 1.3.1 | Circuit IR Definition | `ir_circuit.proto` extending IR schemas for multi-gene circuits; promoter, terminator, regulatory element representation; circuit topology (linear/circular) | 3 PW | 1.1.1 | REQ-FUNC-070 |
| 1.3.2 | Linker Passes | Four composition check passes: promoter conflict, resource competition, splicing interference, RNA-RNA interaction; three-valued verdict composition across components | 6 PW | 1.3.1, 1.1.2 | REQ-FUNC-071–073 |
| 1.3.3 | Toggle Switch Pilot | End-to-end toggle switch circuit verification; validation against DOC-05 Pilot 2 criteria | 2 PW | 1.3.2, 1.2.3 | — |
| 1.3.4 | Phase 3 Integration & V&V | Multi-gene pipeline integration; circuit certificate generation; testing per DOC-05 Section 4.2 Phase 3 criteria | 3 PW | 1.3.1–1.3.3 | REQ-NFR-011 |

**Phase 3 Risk Focus**: RISK-08 (three-valued logic composes incorrectly) is the primary risk during this phase. The composition of three-valued verdicts across circuit components must preserve the soundness guarantee: a circuit-level PASS must imply that every individual component PASS and every composition check PASS is sound. The linker passes (1.3.2) are the largest work package at 6 PW, reflecting the complexity of cross-component constraint checking.

#### Phase 4: Advanced Capabilities (15 person-weeks total)

Phase 4 delivers the mutation explorer and ORF analyzer, completing the full BioCompiler capability set. The phase culminates in the SARS-CoV-2 pilot study and full system V&V against all acceptance criteria.

| WBS | Work Package | Deliverable | Effort (PW) | Dependencies | Requirements |
|---|---|---|---|---|---|
| 1.4.1 | Mutation Explorer | `mutation/` with grammar-guided mutation decomposition, independence exploitation, legal combination enumeration, constraint conflict detection | 5 PW | 1.1.2 | REQ-FUNC-080–083 |
| 1.4.2 | ORF Analyzer | `orf/` with multi-frame construction, shared constraint set computation, coupling classification, inter-frame conflict detection | 4 PW | 1.1.3 | REQ-FUNC-090–093 |
| 1.4.3 | SARS-CoV-2 Pilot | End-to-end ORF1a/ORF1b overlapping frame analysis; validation against DOC-05 Pilot 3 criteria | 2 PW | 1.4.2 | — |
| 1.4.4 | Full System V&V | Complete system verification and validation per DOC-05 Section 4; all acceptance criteria AC-01 through AC-07; cross-platform reproducibility testing; biological validation suite execution | 4 PW | 1.4.1–1.4.3 | DOC-05 Section 4 |

**Phase 4 Risk Focus**: RISK-07 (soundness violation post-release) is mitigated by the comprehensive V&V work package (1.4.4), which includes adversarial soundness testing, biological validation, cross-platform reproducibility testing, and all pilot studies. The full system V&V is the final gate before the v1.0.0 release.

#### Total Effort Summary

| Phase | Person-Weeks | Calendar Months | Key Milestone |
|---|---|---|---|
| Phase 1: Foundation | 18 PW | 4 | M1: Phase 1 complete |
| Phase 2: Splicing-Aware Gene Design | 17 PW | 5 | M2: Phase 2 complete |
| Phase 3: Compositional Verification | 14 PW | 6 | M3: Phase 3 complete |
| Phase 4: Advanced Capabilities | 15 PW | 6 | M4: Phase 4 complete |
| **Total** | **64 PW** | **21** | **M4: v1.0.0 release** |

---

## 3. Schedule

### 3.1 Gantt Summary — ASCII Gantt Chart (21 Months)

```
Month        1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
             |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |

Phase 1:     ████████████████▏
Foundation   1.1.1 ██▏
             1.1.2   ██████▏
             1.1.3       ████▏
             1.1.4     ███▏
             1.1.5         ██▏
                                  M1

Phase 2:                    █████████████████████
Splicing                    1.2.1 █████▏
Design                      1.2.2     ████▏
                            1.2.3   ███▏
                            1.2.4       ██▏
                            1.2.5        ███▏
                                                    M2

Phase 3:                                     ████████████████████████████
Composition                                  1.3.1 ███▏
                                             1.3.2   ██████▏
                                             1.3.3         ██▏
                                             1.3.4          ███▏
                                                                        M3

Phase 4:                                                          ████████████████████████████
Advanced                                                          1.4.1 █████▏
                                                                  1.4.2   ████▏
                                                                  1.4.3      ██▏
                                                                  1.4.4        ████▏
                                                                                            M4

Proj Mgmt:   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
(Ongoing)                                                                                  
```

**Legend**: `█` = active work, `▏` = boundary, `░` = project management overhead

**Phase Overlap Policy**: Phases do not overlap. Each phase transition requires passing the quality gate defined in Section 7 before subsequent phase work may begin. The only exception is Project Management (WBS 1.5), which runs continuously across all phases.

### 3.2 Milestone Definitions

| Milestone | Month | Phase | Criteria | Deliverables |
|---|---|---|---|---|
| **M1** | Month 4 | Phase 1 complete | Quality Gate QG-1 passed (Section 7); all Phase 1 work packages (1.1.1–1.1.5) completed; Phase 1 acceptance criteria in DOC-05 Section 4.2 met; Phase 1 integration test passing | IR schemas (`.proto` files), NDFST splicing grammar for human, type system prototype, AlphaFold FFI adapter, Phase 1 integration test report |
| **M2** | Month 9 | Phase 2 complete | Quality Gate QG-2 passed (Section 7); all Phase 2 work packages (1.2.1–1.2.5) completed; Phase 2 acceptance criteria in DOC-05 Section 4.2 met; all Phase 2 soundness tests passing (zero false PASS); GFP pilot study reviewed by domain expert | CSP codon optimizer, splicing constraint integration, certificate generator with standalone checker, GFP pilot study report, Phase 2 V&V report |
| **M3** | Month 15 | Phase 3 complete | Quality Gate QG-3 passed (Section 7); all Phase 3 work packages (1.3.1–1.3.4) completed; Phase 3 acceptance criteria in DOC-05 Section 4.2 met; toggle switch pilot study reviewed by domain expert | Circuit IR schema, four linker passes, circuit certificate generation, toggle switch pilot study report, Phase 3 V&V report |
| **M4** | Month 21 | Phase 4 complete — all acceptance criteria met | Quality Gate QG-4 passed (Section 7); all Phase 4 work packages (1.4.1–1.4.4) completed; all system-level acceptance criteria AC-01 through AC-07 met (DOC-05 Section 4.1); SARS-CoV-2 pilot study reviewed by domain expert; v1.0.0 release tagged | Mutation explorer, ORF analyzer, SARS-CoV-2 pilot study report, full system V&V report, v1.0.0 release package with all documentation |

**Milestone Review Process**: Each milestone review SHALL be conducted as a formal review with the project team and at least one external domain expert. The review SHALL examine: (1) all deliverables against acceptance criteria, (2) risk register status and any new risks, (3) schedule variance and effort variance, and (4) go/no-go decision for the next phase. The review produces a signed milestone acceptance record.

---

## 4. Staffing

### 4.1 Required Roles

| Role | Count | Skills Required | Phase(s) |
|---|---|---|---|
| **Compiler Engineer** | 1 FTE | Finite-state transducer (FST/NDFST) design and implementation; protocol buffer schema design; compiler pipeline architecture with staged transformations; DFA construction and minimization; Python 3.10+ with type hints and `match` statements; BDD-based symbolic set representation; git workflow with feature branches and code review | Phase 1–4 |
| **Formal Methods Engineer** | 1 FTE | Type system design with subtyping and three-valued logic; constraint satisfaction problems (CSP) with Z3 or equivalent SMT solver; abstract interpretation for sound over-approximation; minimal unsatisfiable subset (MUS) computation; proof-carrying code and certificate generation; three-valued logic algebra and composition laws; property-based testing with Hypothesis | Phase 1–4 |
| **Computational Biologist** | 0.5 FTE | Splicing biology: splice site consensus sequences, branch point motifs, polypyrimidine tracts, ESE/ESS/ISE/ISS regulatory elements; GENCODE GFF3 annotation parsing and splice site PWM construction; codon optimization and CAI computation using organism-specific codon usage tables; AlphaFold2/3 input preparation and output interpretation; cellular context modeling for alternative splicing; SARS-CoV-2 genome organization and overlapping reading frames | Phase 1–4 |

**Role Allocation by Work Package**:

| Work Package | Compiler Engineer | Formal Methods Engineer | Computational Biologist |
|---|---|---|---|
| 1.1.1 IR Schema Definition | Lead | Review | Consult |
| 1.1.2 Splicing Grammar | Lead | Review | Co-lead |
| 1.1.3 Type System Prototype | Consult | Lead | Review |
| 1.1.4 AlphaFold FFI Adapter | Implement | Review | Co-lead |
| 1.1.5 Phase 1 Integration Test | Lead | Co-lead | Review |
| 1.2.1 CSP Codon Optimizer | Consult | Lead | Review |
| 1.2.2 Splicing Constraint Integration | Consult | Lead | Co-lead |
| 1.2.3 Certificate Generator | Review | Lead | Review |
| 1.2.4 GFP Pilot Study | Review | Review | Lead |
| 1.2.5 Phase 2 Integration & V&V | Co-lead | Lead | Review |
| 1.3.1 Circuit IR Definition | Lead | Review | Consult |
| 1.3.2 Linker Passes | Review | Lead | Co-lead |
| 1.3.3 Toggle Switch Pilot | Review | Review | Lead |
| 1.3.4 Phase 3 Integration & V&V | Co-lead | Lead | Review |
| 1.4.1 Mutation Explorer | Lead | Co-lead | Consult |
| 1.4.2 ORF Analyzer | Review | Lead | Co-lead |
| 1.4.3 SARS-CoV-2 Pilot | Review | Review | Lead |
| 1.4.4 Full System V&V | Co-lead | Lead | Co-lead |

### 4.2 Team Configurations

Three staffing configurations are provided, trading team size against calendar duration. The recommended configuration balances cost efficiency with schedule risk.

| Configuration | FTE | Calendar Duration | Monthly Cost Factor | Notes |
|---|---|---|---|---|
| **Minimum** | 2.0 FTE (Compiler Engineer + Formal Methods Engineer; Computational Biologist on contract) | 21 months | 1.0× (baseline) | Computational Biologist engaged at 0.5 FTE via contract for biology-specific work packages (1.1.2, 1.1.4, 1.2.2, 1.2.4, 1.3.2, 1.3.3, 1.4.2, 1.4.3). Phases execute sequentially with no overlap. Lowest total cost but longest schedule. Risk: if either FTE is unavailable, the critical path is immediately impacted since each FTE is the sole contributor to major work packages. |
| **Recommended** | 2.5 FTE (Compiler Engineer + Formal Methods Engineer + Computational Biologist at 0.5 FTE) | 18 months | 1.1× | Computational Biologist is a persistent team member at 0.5 FTE, enabling continuous biology-domain input and faster review cycles. Some work packages within phases can overlap when they have no mutual dependency (e.g., 1.2.1 and 1.2.3 can proceed in parallel). Three-month schedule reduction versus Minimum. |
| **Accelerated** | 3.5 FTE (Compiler Engineer + Formal Methods Engineer + Computational Biologist at 0.5 FTE + additional Compiler Engineer + additional Formal Methods Engineer at 0.5 FTE each) | 14 months | 1.5× | Additional staff enable significant within-phase parallelism: independent work packages (e.g., 1.1.2 and 1.1.3, or 1.4.1 and 1.4.2) execute simultaneously. Requires careful coordination to avoid merge conflicts and ensure IR schema consistency. Highest cost but shortest schedule. Risk: coordination overhead may reduce effective parallelism below theoretical maximum. |

**Decision**: The Recommended configuration (2.5 FTE, 18 months) is assumed for the schedule in Section 3.1. The Gantt chart is drawn for the 21-month minimum configuration; the recommended configuration compresses Phase 2–4 by approximately 3 months through within-phase parallelism.

---

## 5. Risk Register

| ID | Risk | Probability | Impact | Severity | Mitigation | Contingency |
|---|---|---|---|---|---|---|
| RISK-01 | **Splicing grammar precision insufficient** — The NDFST constructed from GENCODE annotations may produce isoform sets that either miss biologically real splice isoforms (false negatives, leading to incomplete safety guarantees) or produce spurious isoforms that do not occur in vivo (false positives, leading to over-conservative designs that reject feasible gene optimizations). Both cases undermine the system's value proposition: false negatives violate soundness, while false positives reduce practical utility. | Medium — GENCODE annotations are high-quality for canonical splice sites but coverage of non-canonical sites (GC-AG, AT-AC) and tissue-specific regulatory elements (ESE/ESS/ISE/ISS) is less complete. Empirical data from ENCODE suggests that approximately 10–20% of tissue-specific splice sites may be absent from current annotations. | High — If the grammar misses real splice sites, the NoCrypticSplice guarantee is unsound (a designed gene could contain a cryptic splice site that the grammar fails to detect). If the grammar is over-conservative, the optimizer will frequently report INFEASIBLE for designs that are biologically feasible, reducing user trust and adoption. | High | (1) Use the most comprehensive annotation source available (GENCODE v44 for human) and supplement with ENCODE tissue-specific splice junction data for cellular context parameterization. (2) Implement conservative (over-approximating) splice site detection: the NDFST SHALL include weak splice sites that are below the functional threshold but above the detection threshold, reporting them as UNCERTAIN rather than FAIL, so that soundness is preserved even with incomplete data. (3) Validate against MaxEntScan reference scores (DOC-05 TC-V-005) with Pearson correlation ≥ 0.95 to ensure scoring accuracy. (4) Perform biological validation against 50 GENCODE-annotated genes (DOC-05 TC-V-001) to measure recall and precision empirically. | If the grammar precision is insufficient for reliable guarantees, (1) downgrade the SpliceCorrect predicate to produce UNCERTAIN verdicts for tissue-specific contexts where annotation coverage is below a validated threshold, preserving soundness at the cost of reduced informativeness; (2) introduce a grammar extension API that allows domain experts to manually add splice site rules for specific cell types, reducing dependency on GENCODE coverage; (3) engage a second computational biologist to conduct a manual review of the grammar for the most commonly used cell types (HEK293T, HepG2, K562). |
| RISK-02 | **Type system too conservative** — The three-valued type system may produce FAIL or UNCERTAIN verdicts for sequences that are biologically correct but fall outside the grammar's formal characterization. This over-conservativeness manifests in two ways: (a) the type system rejects sequences that would function correctly in vivo because the grammar's sound over-approximation includes spurious isoforms, and (b) the type system reports UNCERTAIN for properties that could be resolved with additional biological knowledge not captured in the grammar. Over-conservativeness reduces the system's practical value by generating too many INFEASIBLE reports and too many UNCERTAIN verdicts, discouraging users who expect actionable results. | Medium — The tension between soundness (no false PASS) and completeness (no false FAIL/UNCERTAIN) is inherent in any formal verification system. The BioCompiler design explicitly prioritizes soundness (REQ-NFR-011), which inevitably introduces some over-approximation. The practical impact depends on the gap between the grammar's over-approximation and biological reality, which cannot be fully predicted before implementation. | Medium — Over-conservative verdicts reduce user adoption and trust but do not compromise the system's core safety property. Users who encounter frequent INFEASIBLE or UNCERTAIN results may abandon the tool in favor of less rigorous but more permissive alternatives. | Medium | (1) Design the type system with a clear three-valued semantics where UNCERTAIN is a first-class verdict, not a degenerate case, so that users receive meaningful information about what is known and what is not. (2) Provide detailed knowledge-gap specifications for UNCERTAIN verdicts, enabling users to understand what additional information would resolve the uncertainty. (3) Implement the subtyping mechanism (REQ-FUNC-033) so that more restrictive cellular contexts produce more specific (less conservative) verdicts, giving users a lever to reduce over-approximation by specifying more precise contexts. (4) Include the GFP pilot study (WBS 1.2.4) as an early test of practical conservativeness. | If the type system proves too conservative for practical use, (1) introduce a "relaxed mode" flag that permits UNCERTAIN verdicts in certificates with explicit disclaimers, allowing users who accept residual uncertainty to proceed with designs that are formally sound but incomplete; (2) add a "grammar tuning" interface that allows users to adjust splice site detection thresholds for their specific use case, at the cost of weakened formal guarantees; (3) invest additional effort in WBS 1.1.2 to improve grammar precision, reducing the over-approximation gap. |
| RISK-03 | **CSP solver too slow for long proteins** — The Z3-based CSP solver may require unacceptable computation time for proteins with many amino acid residues (e.g., > 500 aa), where the number of codon variables and constraints grows linearly with protein length. Each amino acid position introduces up to 6 synonymous codon choices (for leucine and arginine) and multiple constraints (CAI contribution, GC contribution, restriction site avoidance, splice site avoidance), creating a large search space. If the solver requires hours or days for a single optimization, the system becomes impractical for real-world gene design, which often targets proteins in the 300–2000 aa range. | Low — Z3 is a highly optimized SMT solver capable of handling constraint problems with thousands of variables. The BioCompiler CSP formulation uses linear arithmetic constraints (CAI, GC percentage) and Boolean constraints (restriction site avoidance, splice site avoidance), which are within Z3's core competency. However, the interaction between many constraints may create difficult search spaces, and performance cannot be guaranteed without empirical measurement. | High — If the solver cannot deliver results within the required time bounds (REQ-NFR-010: single-gene analysis in ≤ 10 minutes for ≤ 5 kb genes), the system fails its non-functional performance requirements and cannot be used for interactive gene design. Long solver runtimes also impede V&V activities, as the test suite must execute within practical time limits. | Medium | (1) Profile Z3 performance on representative problem sizes (200 aa, 500 aa, 1000 aa) during WBS 1.2.1 to establish empirical performance baselines early. (2) Decompose the CSP by codon position: exploit the fact that many constraints are local (restriction sites span ≤ 8 nt, splice sites span ≤ 30 nt), enabling the solver to reason about independent subproblems. (3) Use incremental solving: start with a relaxed constraint set and progressively tighten, enabling early termination if the relaxed problem is already infeasible. (4) Implement a solver timeout (default 300 seconds) with partial result reporting: if the solver times out, report the best feasible solution found so far (if any) with an UNCERTAIN verdict on optimality. | If the Z3 solver is too slow for long proteins, (1) implement a hybrid solver strategy: use Z3 for the core constraint satisfaction and a greedy heuristic for the optimization objective, reducing the solver's workload to feasibility only; (2) switch to a dedicated integer programming solver (e.g., Gurobi, CPLEX) that may handle large problems more efficiently than a general-purpose SMT solver; (3) partition long proteins into overlapping windows and solve each window independently, merging results with conflict resolution at boundaries; (4) provide a "quick mode" flag that uses heuristic codon assignment with subsequent verification (type-check the heuristic result rather than optimizing from scratch). |
| RISK-04 | **IR schema design mistakes propagate** — The Intermediate Representation schemas (IR-Seq, IR-Peptide, IR-Structure, IR-Circuit) are the shared data contracts between all pipeline components. If a schema is designed incorrectly — missing a required field, choosing an inappropriate data type, or defining an invariant that is too restrictive or too permissive — the mistake propagates to every component that depends on that schema. Because IR schemas are defined in protocol buffers with backward-compatible evolution rules, fixing a fundamental design mistake may require a breaking schema change that affects all downstream components simultaneously, causing cascading rework across multiple work packages. | Low — The schema design is guided by the detailed requirements in DOC-01 and the interface contracts in DOC-04. The schema author (Compiler Engineer) has experience with protocol buffer schema design. The schemas are small (four IR levels, each with < 20 message types), limiting the surface area for design errors. | High — A schema redesign during Phase 2 or later would require updating all components that produce or consume the affected IR level, potentially invalidating completed work packages and requiring re-implementation of serialization/deserialization logic, invariant checks, and integration tests. The cost of late schema changes increases quadratically with the number of dependent components. | Medium | (1) Invest adequate time in WBS 1.1.1 (3 PW, the largest non-splicing work package in Phase 1) to produce thorough, reviewed schema designs before downstream components begin implementation. (2) Conduct a formal schema review with both the Compiler Engineer and Formal Methods Engineer before closing WBS 1.1.1, using the interface contracts in DOC-04 as review criteria. (3) Implement schema validation tools early: a schema linter that checks for common design errors (missing required fields, inconsistent naming, invariant violations) and a fuzz tester that generates random IR instances to verify that invariants hold under all valid inputs. (4) Use protocol buffer `reserved` tags and field numbers to prevent accidental reuse of deprecated fields. | If a schema design mistake is discovered after Phase 1, (1) assess the blast radius: identify all components that consume or produce the affected IR level and estimate rework effort; (2) implement a schema migration path using protocol buffer's built-in backward compatibility, adding new fields to replace deprecated ones while maintaining read support for the old format; (3) if the mistake is fundamental (e.g., a wrong data type that cannot be migrated), schedule a focused "schema fix" sprint that updates all affected components in parallel, accepting the schedule impact; (4) add a post-mortem item to the Project Management process to identify why the review process failed to catch the error. |
| RISK-05 | **External tool API changes** — The BioCompiler FFI adapter framework depends on external tools (AlphaFold2/3, ColabFold, RoseTTAFold, NetPhos, PhosphoSitePlus, dbPTM, MusiteDeep) that are independently developed and may change their input/output formats, command-line interfaces, or Python APIs between BioCompiler releases. An API change that breaks an FFI adapter renders the corresponding pipeline feature non-functional until the adapter is updated. Because FFI output is treated as non-deterministic (REQ-FUNC-053) and SLOT fields are optional, a broken FFI adapter does not compromise the core pipeline's determinism or soundness, but it reduces the system's utility for users who rely on structure prediction or PTM analysis. | Medium — AlphaFold2 and ColabFold are relatively stable (AlphaFold2 has not had a breaking API change since its open-source release), but AlphaFold3 is newer and may undergo API evolution. NetPhos and other PTM prediction tools are web services that may change their REST API or output format without notice. The probability of at least one breaking API change during the 21-month development period is moderate. | Medium — A broken FFI adapter is isolated to the affected tool and does not propagate to other components. The core pipeline (scan → splice → translate → type-check → optimize → certificate) operates independently of FFI stages. However, users who need structure validation or PTM analysis will be unable to use those features until the adapter is fixed, reducing the system's value for those use cases. | Medium | (1) Design the FFI adapter framework with a versioned adapter interface: each adapter declares the external tool version it supports, and the adapter registry rejects adapters for unsupported versions. (2) Implement comprehensive output validation in every adapter (REQ-FUNC-052): if the external tool's output format changes, the validation layer catches the discrepancy and raises `OutputValidationError` rather than silently producing incorrect results. (3) Pin external tool versions in the build and test environment, and do not upgrade without explicit adapter testing. (4) Provide mock adapters for all FFI tools, enabling the core pipeline to be tested and used without any external tool dependencies (REQ-NFR-042). | If an external tool API changes, (1) the output validation layer will immediately detect the change and raise `OutputValidationError`, preventing silent data corruption; (2) update the affected adapter in a dedicated bug-fix sprint, targeting a turnaround of ≤ 5 business days for non-breaking format changes and ≤ 15 business days for breaking API changes; (3) if the tool has changed fundamentally (e.g., a different output schema), implement a new adapter version alongside the deprecated one, supporting both old and new formats during a transition period; (4) communicate the issue and fix timeline to users via the project communication channels. |
| RISK-06 | **No community adoption** — BioCompiler introduces a fundamentally different approach to gene design (deterministic formal verification with three-valued logic) compared to existing tools (heuristic optimization with probabilistic scoring). The target user community (synthetic biologists, circuit designers, virologists) may be unfamiliar with or resistant to formal methods concepts, preferring tools that produce binary accept/reject answers rather than three-valued PASS/FAIL/UNCERTAIN verdicts. If the community does not adopt BioCompiler, the project fails to achieve its impact goals regardless of technical success. | Medium — The synthetic biology community has limited exposure to formal verification methods. Existing gene design tools (IDT, Genscript, GeneDesign, DNAWorks) use heuristic approaches that are familiar to users. The three-valued logic paradigm, while more rigorous, adds conceptual complexity. However, regulatory pressure for formally verified gene designs (especially for gene therapy and mRNA vaccines) is increasing, creating a potential market pull. | Medium — Lack of adoption does not affect the technical correctness of the system but represents a strategic failure. The effort invested in developing BioCompiler yields limited return if the tool is not used. Adoption failure also reduces the likelihood of community contributions (bug reports, feature requests, new organism grammars) that would improve the system over time. | Medium | (1) Invest in user-facing documentation that explains three-valued logic in biological terms (PASS = "guaranteed correct," FAIL = "guaranteed incorrect," UNCERTAIN = "cannot determine with current knowledge") rather than formal methods jargon. (2) Provide a "quick start" tutorial that walks through the GFP pilot study as a concrete example, demonstrating that the tool produces practical, actionable results. (3) Publish the guarantee certificate format as an open standard, enabling third-party tools to consume and verify BioCompiler certificates without depending on the BioCompiler codebase. (4) Engage with regulatory agencies (FDA, EMA) early to understand their requirements for formally verified gene designs, positioning BioCompiler as a compliance tool. (5) Present results at synthetic biology conferences (SB, IWBDA, SEED) to build awareness. | If community adoption is low, (1) conduct user interviews to identify specific barriers (conceptual complexity, missing features, performance issues, integration difficulties) and address them in a post-v1.0 maintenance release; (2) develop a web-based GUI that hides the three-valued logic behind a simpler interface (green/yellow/red indicators) to reduce conceptual barrier; (3) partner with a gene synthesis company to offer BioCompiler verification as a value-added service alongside gene synthesis, integrating the tool into an existing workflow rather than requiring standalone adoption; (4) pivot to a "verification-as-a-service" model where users submit gene designs for verification without installing the tool locally. |
| RISK-07 | **Soundness violation post-release** — A soundness violation occurs when the type system produces a PASS verdict for a sequence that actually violates a declared constraint (a false PASS). This is the most severe possible failure for the BioCompiler system: it violates REQ-NFR-011, undermines the guarantee certificate's value, and could lead to the deployment of a gene design that contains a biological error (e.g., a cryptic splice site) that the user was assured was absent. Soundness violations are particularly dangerous because they are invisible to users who trust the certificate. | Low — The type system is designed with soundness as its primary invariant, and the adversarial soundness test suite (DOC-05 Section 2.4) specifically targets potential soundness violations with known-violation inputs. The three-valued logic design explicitly sacrifices completeness for soundness: the system will produce FAIL or UNCERTAIN verdicts for sequences where it cannot guarantee correctness, but it will never produce a PASS verdict for an incorrect sequence (assuming the grammar is sound). | Critical — A soundness violation destroys the system's core value proposition: formal guarantees. If a user deploys a gene design based on a false PASS certificate, the resulting biological failure (e.g., unintended splicing, incorrect translation) could have serious consequences in therapeutic applications. Even a single documented soundness violation would severely damage the system's credibility and could not be remediated by patching alone — all previously issued certificates would be suspect. | High | (1) Maintain the comprehensive soundness test suite (TC-SND-001 through TC-SND-008) as a mandatory gate for every release. Any test failure blocks the release unconditionally. (2) Implement property-based testing with Hypothesis that generates random sequences with injected violations, providing coverage beyond the hand-crafted adversarial tests. (3) Require 100% line coverage for COMP-05 (Type System) and COMP-06 (Optimizer), ensuring that every code path in the verification-critical components is exercised by at least one test. (4) Conduct an independent formal review of the type system's soundness proof by a qualified formal methods researcher, separate from the development team. (5) Implement the certificate standalone checker as a completely independent codebase with no shared implementation with the type system, providing a second implementation for cross-validation. | If a soundness violation is discovered post-release, (1) immediately issue a critical security advisory to all known users, identifying the specific predicate and input conditions that trigger the false PASS; (2) revoke all certificates that could be affected by the violation (using the certificate's predicate and parameter metadata to determine the blast radius); (3) implement and release a patch within 72 hours, adding the violating input to the soundness test suite as a regression test; (4) conduct a root cause analysis and post-mortem, examining why the existing test suite failed to catch the violation; (5) if the violation is in the splicing grammar (RISK-01 interaction), also update the grammar rules to prevent the specific false PASS condition; (6) engage an external auditor to review the type system for additional potential violations before the next release. |
| RISK-08 | **Three-valued logic composes incorrectly** — The BioCompiler type system uses three-valued logic (PASS / FAIL / UNCERTAIN) for individual predicate verdicts, and these verdicts are composed across multiple predicates (conjunction within a single gene) and across multiple components (circuit-level composition via linker passes). If the composition rules are defined incorrectly — for example, if the conjunction of two UNCERTAIN verdicts is defined as PASS instead of UNCERTAIN, or if the composition of a PASS in one circuit component with an UNCERTAIN in another produces PASS instead of UNCERTAIN — the resulting circuit-level verdict will be unsound. This risk is particularly insidious because the error may not be detectable by single-component testing; it only manifests when multiple components are verified together in a circuit. | Low — The three-valued logic conjunction algebra is well-defined in formal logic (Kleene's strong three-valued logic or Łukasiewicz's three-valued logic), and the BioCompiler design uses a defined algebra that prioritizes soundness (FAIL dominates, UNCERTAIN is absorbed by FAIL but propagates past PASS). The composition rules are documented in DOC-03 and DOC-05, and the test suite includes specific tests for verdict composition (TC-U-034, TC-U-035). However, the interaction between single-gene verdict composition and circuit-level verdict composition across linker passes introduces additional complexity that may harbor composition errors. | High — An incorrect composition rule would violate the soundness guarantee (REQ-NFR-011) at the circuit level: a circuit certificate could claim PASS when one or more composition checks should have returned FAIL or UNCERTAIN. This is a variant of RISK-07 but specific to the composition logic, making it harder to detect because it only manifests in multi-gene circuits. | Medium | (1) Formally specify the three-valued conjunction algebra and prove its soundness property before implementation. The proof SHALL demonstrate that for any combination of verdicts, the composed verdict is at least as strong as the weakest individual verdict (FAIL > UNCERTAIN > PASS in terms of restrictiveness). (2) Implement the composition logic as a single, well-tested function (not duplicated across components), ensuring that all composition operations use the same algebra. (3) Include specific tests for all 9 possible two-predicate conjunction combinations (PASS∧PASS, PASS∧FAIL, PASS∧UNCERTAIN, FAIL∧PASS, FAIL∧FAIL, FAIL∧UNCERTAIN, UNCERTAIN∧PASS, UNCERTAIN∧FAIL, UNCERTAIN∧UNCERTAIN) and all 27 possible three-predicate conjunction combinations. (4) Test circuit-level composition with the toggle switch pilot study (WBS 1.3.3), which exercises composition across two genes with multiple composition check types. | If an incorrect composition rule is discovered, (1) the fix is typically localized to the composition function and does not require changes to individual predicate logic; (2) re-run all soundness tests and integration tests to verify that the corrected composition rule does not introduce new violations; (3) re-generate all certificates from the toggle switch pilot study and GFP pilot study with the corrected composition rule; (4) if the incorrect rule was present in a released version, follow the soundness violation contingency from RISK-07, as this constitutes a false PASS at the circuit level. |

---

## 6. Configuration Management

### 6.1 Version Control

**Repository**: A single Git repository (`biocompiler`) SHALL contain all source code, IR schemas, test suites, documentation, and build configuration.

**Branching Strategy**: The project SHALL use a feature-branch workflow:

| Branch Type | Naming Convention | Purpose | Merge Policy |
|---|---|---|---|
| `main` | `main` | Stable, releasable code. Every commit on `main` SHALL pass the full CI suite. | Merge only via pull request with ≥ 1 approving review. Squash-merge to maintain clean history. |
| Feature branch | `feature/<WBS-ID>-<short-description>` (e.g., `feature/1.1.2-splicing-grammar`) | Development of a specific work package. Created from `main` at the start of the work package. | Merge to `main` via pull request after all work package acceptance criteria are met. Delete after merge. |
| Bug-fix branch | `fix/<issue-number>-<short-description>` | Correction of a reported defect. Created from `main`. | Merge to `main` via pull request with regression test. |
| Release branch | `release/vX.Y.Z` | Preparation of a specific release. Created from `main` when all phase acceptance criteria are met. Allows last-minute fixes without blocking ongoing development. | Merge to `main` and tag with version. |

**Main Branch Protection Rules**:

- No direct pushes to `main`; all changes via pull request.
- Pull requests require ≥ 1 approving review from a team member other than the author.
- Pull requests must pass the full CI suite (lint, unit tests, integration tests) before merge.
- The Compiler Engineer and Formal Methods Engineer SHALL cross-review each other's pull requests. The Computational Biologist SHALL review biology-domain changes (splicing grammar, codon tables, cellular context models).

**Commit Conventions**:

- Commit messages SHALL follow the Conventional Commits format: `type(scope): description` (e.g., `feat(splicing): implement NDFST isoform set computation`, `fix(optimizer): correct MUS minimality check`).
- Each commit SHALL reference the relevant WBS ID in the commit message footer (e.g., `WBS: 1.1.2`).

### 6.2 Release Strategy

| Version | Phase | Description | Criteria for Tag |
|---|---|---|---|
| `v0.1.0-alpha` | End of Phase 1 | First internal release with core pipeline (scanner, splicing, translation, FFI). Not suitable for external use. | M1 passed; Phase 1 integration test passing; basic CLI operational |
| `v0.2.0-alpha` | Mid Phase 2 | Core pipeline + type system prototype. Internal use only. | Type system producing three-valued verdicts; unit tests for COMP-05 passing |
| `v0.3.0-alpha` | End of Phase 2 | Core pipeline + optimizer + certificate generator. First release with formal verification capabilities. Internal use and select external reviewers. | M2 passed; GFP pilot study complete; all Phase 2 soundness tests passing |
| `v0.4.0-beta` | End of Phase 3 | Full single-gene pipeline + circuit verification. First release suitable for external beta testing. | M3 passed; toggle switch pilot study complete; circuit certificates passing standalone verification |
| `v0.5.0-beta` | Mid Phase 4 | All components implemented. Pre-release for final V&V. | Mutation explorer and ORF analyzer operational; SARS-CoV-2 pilot study in progress |
| `v1.0.0` | End of Phase 4 | Full BioCompiler release. All acceptance criteria met. | M4 passed; all acceptance criteria AC-01 through AC-07 met; all three pilot studies reviewed by domain experts; full system V&V complete |

**Version Numbering**: The project follows Semantic Versioning 2.0.0. Pre-release versions use `-alpha` (internal use, no stability guarantees) and `-beta` (feature-complete, external review, no stability guarantees for APIs). The `v1.0.0` release commits to backward compatibility for the IR schema, CLI interface, and Python API.

**Release Artifacts**: Each release SHALL include:

- Source code archive (tagged Git commit)
- Built Python package (wheel and sdist)
- IR schema files (compiled `.proto` descriptors)
- Test suite execution report with coverage
- Certificate checker binary (standalone)
- Release notes with change log, known issues, and migration guide (for breaking changes)

### 6.3 Baseline Definition

A **baseline** is a formally approved snapshot of the project's configuration items at a milestone. Each baseline is established at the corresponding milestone review and placed under configuration control: changes to baselined items require a formal change request.

| Baseline | Milestone | Configuration Items |
|---|---|---|
| **BL-1** (Requirements Baseline) | Project start | DOC-01 (SRS), DOC-05 (SVVP), DOC-07 (this document) — all at version 1.0.0-draft, baselined after initial review |
| **BL-2** (Architecture Baseline) | M1 (Month 4) | DOC-02 (SAD), DOC-04 (ICD), DOC-06 (Design Rationale), IR schemas (`.proto` files), Phase 1 source code and tests |
| **BL-3** (Design Baseline) | M2 (Month 9) | DOC-03 (SDD), all Phase 1 and Phase 2 source code and tests, `v0.3.0-alpha` release package |
| **BL-4** (Product Baseline) | M4 (Month 21) | All documentation (DOC-00 through DOC-10), all source code and tests, `v1.0.0` release package, all pilot study reports, full V&V report |

**Change Control**: Changes to items in the current baseline require a Change Request (CR) that documents: (1) the proposed change, (2) the affected configuration items, (3) the impact on schedule, effort, and risk, (4) the traceability to affected requirements. CRs for BL-1 and BL-2 items require approval from both the Compiler Engineer and Formal Methods Engineer. CRs for BL-3 and BL-4 items additionally require Computational Biologist approval if biology-domain items are affected. Emergency CRs (e.g., soundness violation fixes) may be approved by a single team member with retroactive review within 48 hours.

---

## 7. Quality Gates

Quality gates define the mandatory criteria that MUST be satisfied before a phase transition is authorized. Each quality gate is evaluated at the corresponding milestone review. A quality gate that is not fully passed blocks the phase transition; partial passes are not permitted.

| Gate ID | Phase Transition | Criteria | Evaluation Method |
|---|---|---|---|
| **QG-1** | Phase 1 → Phase 2 | (1) All Phase 1 work packages (1.1.1–1.1.5) completed and reviewed. (2) IR schemas validated against DOC-04 interface contracts with zero discrepancies. (3) NDFST splicing grammar produces correct isoform sets for at least 5 test genes (constitutive, cassette exon, mutually exclusive exons, non-canonical GC-AG, tissue-specific). (4) Type system prototype produces three-valued verdicts for all seven type predicates. (5) AlphaFold FFI adapter produces valid IR-Structure with SLOT fields populated and provenance metadata recorded (using mock adapter). (6) Phase 1 integration test passes: end-to-end pipeline for one gene from FASTA input to translated output. (7) Phase 1 unit test coverage ≥ 90%. (8) No open Critical or High severity defects. (9) Phase 1 risk register reviewed; no new unmitigated risks with severity ≥ High. | Formal milestone review with demonstration of all criteria. IR schema review against DOC-04. Test suite execution with coverage report. |
| **QG-2** | Phase 2 → Phase 3 | (1) All Phase 2 work packages (1.2.1–1.2.5) completed and reviewed. (2) CSP optimizer produces correct results for feasible and infeasible test cases (TC-U-040 through TC-U-045). (3) Splicing constraint integration verified: optimized sequences pass NoCrypticSplice type check. (4) Certificate generator produces valid certificates that pass standalone verification (TC-U-060, TC-U-061). (5) All Phase 2 soundness tests pass (TC-SND-001 through TC-SND-008) with zero false PASS verdicts. (6) GFP pilot study completed and reviewed by domain expert; all success criteria met. (7) Phase 2 unit test coverage: ≥ 90% overall, 100% for COMP-05 and COMP-06. (8) No open Critical or High severity defects. (9) `v0.3.0-alpha` release package produced. | Formal milestone review with demonstration of all criteria. Soundness test suite execution with zero-false-PASS verification. GFP pilot study report reviewed by Computational Biologist and external domain expert. |
| **QG-3** | Phase 3 → Phase 4 | (1) All Phase 3 work packages (1.3.1–1.3.4) completed and reviewed. (2) Circuit IR schema validated against DOC-04 interface contracts. (3) All four linker passes produce correct verdicts: promoter conflict (TC-U-070), resource competition (TC-U-071), splicing interference (TC-U-072), RNA interaction (TC-U-073). (4) Circuit certificate generation verified (TC-U-062). (5) Multi-gene pipeline integration test passes (TC-I-008). (6) Toggle switch pilot study completed and reviewed by domain expert; all success criteria met. (7) Three-valued logic composition verified for all 9 two-predicate conjunction combinations and all 27 three-predicate combinations. (8) No open Critical or High severity defects. (9) `v0.4.0-beta` release package produced. | Formal milestone review with demonstration of all criteria. Composition logic audit by Formal Methods Engineer. Toggle switch pilot study report reviewed by Computational Biologist and external domain expert. |
| **QG-4** | Phase 4 → Release | (1) All Phase 4 work packages (1.4.1–1.4.4) completed and reviewed. (2) Mutation explorer produces correct decomposition, enumeration, independence exploitation, and conflict detection (TC-U-080 through TC-U-083). (3) ORF analyzer produces correct multi-frame translations, shared constraint sets, coupling classifications, and frame conflict detection (TC-U-090 through TC-U-093). (4) All system-level acceptance criteria AC-01 through AC-07 met (DOC-05 Section 4.1). (5) SARS-CoV-2 pilot study completed and reviewed by domain expert; all success criteria met. (6) All biological validation tests pass (TC-V-001 through TC-V-005). (7) Cross-platform reproducibility verified (TC-REP-002). (8) No open defects of any severity. (9) All documentation (DOC-00 through DOC-10) at version 1.0.0 and baselined. (10) `v1.0.0` release package produced and verified. | Formal release review with full V&V report. All three pilot study reports reviewed by external domain experts. Documentation completeness review. Release package validation on all target platforms. |

---

## 8. Assumptions and Constraints

| ID | Assumption / Constraint | Type | Impact if Violated | Mitigation / Contingency |
|---|---|---|---|---|
| AC-01 | **The standard genetic code is correct and complete** — The 64 codon-to-amino-acid mappings (including three stop codons) are treated as ground truth (SRS ASSUME-01). The translation engine's correctness depends entirely on the accuracy of the genetic code table. If a previously unknown codon recoding mechanism is discovered (analogous to the historical discovery of selenocysteine and pyrrolysine), the translation engine would produce incorrect amino acid sequences for genes utilizing the new mechanism. | Assumption | Translation engine produces incorrect amino acid sequences for genes utilizing the unknown recoding mechanism. The InFrame type predicate may produce false PASS verdicts if the unknown mechanism creates a reading frame disruption that the standard code does not detect. | The standard genetic code is universally accepted and extremely well-characterized; the probability of a new recoding mechanism that affects the standard 61 sense codons is negligible. The system already handles selenocysteine and pyrrolysine as special cases; a new recoding mechanism could be added via the same extension mechanism (new flag in the translation FST). |
| AC-02 | **GENCODE splice site annotations are biologically valid** — The splicing grammar's correctness depends on the accuracy and completeness of splice site annotations derived from GENCODE (SRS ASSUME-02). If GENCODE annotations contain errors (incorrect splice site positions, missing tissue-specific annotations) or are incomplete (missing rare splice sites, missing non-canonical splice sites), the NDFST will produce isoform sets that do not accurately reflect biological reality. | Assumption | Splicing grammar produces incorrect isoform sets: either missing real isoforms (false negatives, leading to unsound guarantees) or including spurious isoforms (false positives, leading to over-conservative designs). The SpliceCorrect and NoCrypticSplice predicates are directly affected. | Use the most current GENCODE release (v44) and validate against ENCODE tissue-specific junction data. Implement conservative over-approximation (include weak sites as UNCERTAIN) to maintain soundness. Biological validation against 50 GENCODE-annotated genes (TC-V-001) provides empirical recall measurement. |
| AC-03 | **External tools produce correct output given valid input** — The FFI adapter framework does not verify the biological correctness of external tool output; it validates only schema conformance (SRS ASSUME-04). If an external tool (e.g., AlphaFold) produces an incorrect structure prediction, the BioCompiler IR-Structure will contain incorrect SLOT fields. However, because FFI output is non-deterministic (REQ-FUNC-053) and SLOT fields are optional, this does not affect the core pipeline's soundness. | Assumption | IR-Structure contains incorrect protein structure data (e.g., incorrect pLDDT scores, incorrect C-alpha coordinates). Users who rely on structure validation through the FFI may receive misleading information. The guarantee certificate is not affected because certificate validity does not depend on FFI output. | FFI output is clearly labeled as non-deterministic in the IR and certificate. The type system treats FFI-derived predicates as UNCERTAIN by default. Users are advised to validate FFI output independently for safety-critical applications. |
| AC-04 | **Project staffing remains stable throughout development** — The schedule and effort estimates assume that the assigned team members (Compiler Engineer, Formal Methods Engineer, Computational Biologist) are available for the duration of their assigned work packages. If a team member becomes unavailable (resignation, illness, reassignment), the critical path is immediately impacted because each role has specialized skills that cannot be easily replaced. | Constraint | Schedule delay proportional to the replacement ramp-up time (estimated 4–8 weeks for a new Compiler Engineer or Formal Methods Engineer to become productive on the BioCompiler codebase). If the Computational Biologist becomes unavailable, biology-domain work packages (1.1.2, 1.2.2, 1.2.4, 1.3.2, 1.3.3, 1.4.2, 1.4.3) are blocked. | Cross-training: the Compiler Engineer and Formal Methods Engineer SHALL review each other's code, ensuring that either can perform basic maintenance on the other's components in an emergency. Maintain a list of potential replacement candidates. For the Computational Biologist, maintain a consulting relationship with at least one additional computational biologist who can provide emergency coverage. |
| AC-05 | **Python 3.10+ and Z3 remain available and supported** — The BioCompiler system depends on Python 3.10+ (for type hint features and `match` statements) and Z3 (for the CSP optimizer). If Python introduces a breaking change that affects the BioCompiler codebase, or if Z3 becomes unavailable or unsupported, significant rework may be required. | Constraint | If Python introduces a breaking change, the codebase must be updated to conform, potentially affecting all components. If Z3 becomes unavailable, the optimizer (COMP-06) must be re-implemented using an alternative solver (e.g., Google OR-Tools, Gurobi), requiring an estimated 3–5 PW of rework. | Pin the Python version in the build environment (currently 3.10) and do not upgrade without comprehensive regression testing. Pin the Z3 version in requirements.txt. Evaluate alternative solvers during WBS 1.2.1 (CSP Codon Optimizer) to ensure that the optimizer's solver interface is abstract enough to support a backend swap with minimal code changes. |

---

## 9. Communication Plan

| Stakeholder | Frequency | Method | Content |
|---|---|---|---|
| **Project Team** (Compiler Engineer, Formal Methods Engineer, Computational Biologist) | Daily | Asynchronous standup via project chat channel (e.g., Slack, Zulip) | What was accomplished yesterday, what will be worked on today, any blockers. Each team member posts a brief update by 10:00 local time. |
| **Project Team** | Weekly | Synchronous video meeting (60 minutes) | Sprint review: work package progress against plan, risk register updates, upcoming milestones, technical decisions requiring team consensus. Meeting notes published to shared document space within 24 hours. |
| **Project Team** | Per-milestone | Formal milestone review meeting (3–4 hours) | Comprehensive review of all milestone deliverables against quality gate criteria. Includes live demonstration of functionality, test suite execution, and pilot study results. Produces signed milestone acceptance record. External domain expert participates for M2, M3, and M4 reviews. |
| **Technical Lead / PI** | Biweekly | Synchronous meeting or written report (30 minutes / 2 pages) | Executive summary: schedule status (on track / at risk / delayed), effort consumed vs. planned, risk register summary (new risks, changed risks, closed risks), key decisions made, decisions needed. |
| **External Domain Experts** (computational biology, formal methods) | Per-milestone | Email with milestone report attachment; optional video debrief (60 minutes) | Milestone report: phase deliverables, pilot study results, V&V summary, risk register, and specific questions for expert review. Experts are expected to provide written feedback within 10 business days. |
| **Regulatory Stakeholders** (FDA, EMA, or equivalent) | Quarterly | Written report (5–10 pages) | Regulatory briefing: system capabilities relevant to gene design verification, guarantee certificate format and standalone verification process, pilot study results, biological validation results. Intended to build familiarity with the BioCompiler approach for future regulatory submissions. |
| **Open-Source Community** (post-v0.3.0-alpha) | Monthly | GitHub releases, discussion forum posts, blog articles | Release notes, new features, known issues, contribution guidelines, roadmap update. Community feedback collected via GitHub Issues and Discussions. |
| **User Community** (post-v0.4.0-beta) | As needed | User group meetings, tutorials, documentation updates | Feature demonstrations, "how-to" guides for common use cases (GFP design, toggle switch verification, ORF analysis), Q&A sessions. Feedback collected via structured user survey after each pilot study. |

---

## 10. Metrics and Reporting

| Metric | Target | Measurement Method | Frequency |
|---|---|---|---|
| **Schedule Performance Index (SPI)** | ≥ 0.95 (within 5% of planned schedule) | SPI = Earned Value (EV) / Planned Value (PV), where EV is the sum of planned effort for completed work packages and PV is the sum of planned effort for work packages that should be complete by the measurement date. Calculated per phase. | Biweekly (with Technical Lead report) |
| **Effort Variance** | ≤ ±10% of planned effort per work package | Actual effort (tracked via time tracking in project management tool) minus planned effort per work package, expressed as a percentage of planned effort. | Biweekly |
| **Defect Density** | ≤ 5 defects per KLOC (thousand lines of code) at Phase 4 | Total defects reported divided by total lines of production code (excluding tests, comments, blank lines), measured at each phase completion. | Per-milestone |
| **Unit Test Coverage** | ≥ 90% overall; 100% for COMP-05 and COMP-06 | `pytest --cov` line coverage report executed on the full test suite. | Every commit (CI-enforced); reported at milestone reviews |
| **Soundness Test Pass Rate** | 100% (zero false PASS verdicts) | Execution of adversarial soundness test suite (TC-SND-001 through TC-SND-008). Any failure is a Critical defect that blocks release. | Every commit (CI-enforced); reported at milestone reviews |
| **Open Defect Count (by severity)** | Zero Critical or High at phase transitions; zero of any severity at release | Defect tracking system query: count of open defects grouped by severity (Critical / High / Medium / Low). | Weekly |
| **Risk Register Currency** | All risks reviewed within last 30 days | Date of last review for each risk entry in the risk register. Risks not reviewed in > 30 days are flagged as stale. | Weekly |
| **Code Review Turnaround** | ≤ 3 business days from pull request creation to merge | Median time from PR creation to merge, measured across all PRs in the reporting period. | Monthly |
| **CSP Solver Performance** | ≤ 300 seconds for proteins ≤ 500 aa; ≤ 600 seconds for proteins ≤ 1000 aa | Wall-clock time for CSP optimizer execution on benchmark problems, measured on the reference hardware (8-core x86_64, 32 GB RAM). Benchmarks run nightly on the main branch. | Nightly (automated); reported at milestone reviews |
| **Biological Validation Metrics** | GENCODE recall ≥ 100% for annotated isoforms; UniProt translation match ≥ 99%; MaxEntScan correlation ≥ 0.95 | Execution of biological validation test suite (TC-V-001 through TC-V-005). Metrics computed from test results. | Per-milestone (Phase 2 onward) |
| **Certificate Verification Rate** | 100% of generated certificates pass standalone verification | For each generated certificate, run the independent checker and verify CERTIFICATE_VALID result. Tracked as a test case (TC-U-061). | Every commit (CI-enforced); reported at milestone reviews |

---

## Appendix A: References

| Ref ID | Document / Source | Description |
|---|---|---|
| REF-01 | DOC-01: Software Requirements Specification (SRS) | Defines all requirements (functional, non-functional, constraints) for the BioCompiler system. |
| REF-02 | DOC-02: Software Architecture Document (SAD) | Defines the component decomposition, data flow, and architectural decisions. |
| REF-03 | DOC-03: Software Design Document (SDD) | Provides detailed algorithmic specifications, data structures, and invariants. |
| REF-04 | DOC-04: Interface Control Document (ICD) | Defines the exact contracts for every IR schema, API, and FFI boundary. |
| REF-05 | DOC-05: Software Verification & Validation Plan (SVVP) | Specifies how each requirement is verified and validated; defines acceptance criteria. |
| REF-06 | DOC-06: Design Rationale (DR) | Justifies each design decision and documents alternatives considered and rejected. |
| REF-07 | DOC-08: Traceability Matrix (TM) | Complete bidirectional traceability: requirements ↔ architecture ↔ design ↔ interface ↔ test ↔ risk. |
| REF-08 | IEEE 16326-2019 | IEEE Standard for Software Project Management Plans. |
| REF-09 | PMBOK 7th Edition | Project Management Institute, Project Management Body of Knowledge. |
| REF-10 | IEEE 1012-2016 | Standard for System, Software, and Hardware Verification and Validation. |
| REF-11 | ISO/IEC/IEEE 29148:2018 | Systems and software engineering — Life cycle processes — Requirements engineering. |

## Appendix B: Acronyms

| Acronym | Definition |
|---|---|
| BDD | Binary Decision Diagram |
| CAI | Codon Adaptation Index |
| CR | Change Request |
| CSP | Constraint Satisfaction Problem |
| DFA | Deterministic Finite Automaton |
| ESE | Exonic Splicing Enhancer |
| ESS | Exonic Splicing Silencer |
| EV | Earned Value |
| FFI | Foreign Function Interface |
| FST | Finite-State Transducer |
| FTE | Full-Time Equivalent |
| ISE | Intronic Splicing Enhancer |
| ISS | Intronic Splicing Silencer |
| IR | Intermediate Representation |
| KLOC | Thousand Lines of Code |
| MUS | Minimal Unsatisfiable Subset |
| NDFST | Non-Deterministic Finite-State Transducer |
| ORF | Open Reading Frame |
| PW | Person-Week |
| PV | Planned Value |
| QG | Quality Gate |
| SECIS | Selenocysteine Insertion Sequence |
| SMT | Satisfiability Modulo Theories |
| SPI | Schedule Performance Index |
| ULP | Unit in the Last Place |
| V&V | Verification and Validation |
| WBS | Work Breakdown Structure |
