# DOC-01: Software Requirements Specification (SRS)

| Field | Value |
|---|---|
| **Document ID** | DOC-01 |
| **Version** | 12.0.0 |
| **Status** | Current |
| **Date** | 2026-06-07 |
| **Prepared By** | BioCompiler Project Team |
| **Reviewed By** | [TBD at baseline review] |
| **Approved By** | [TBD at baseline review] |

---

## 1. Introduction

### 1.1 Purpose

This document constitutes the Software Requirements Specification (SRS) for the **BioCompiler** system, prepared in accordance with IEEE 830-1993 and ISO/IEC/IEEE 29148:2018. Its purpose is to define, in complete, unambiguous, and testable terms, what the BioCompiler system SHALL do, what it SHALL NOT do, and under what constraints it SHALL operate.

BioCompiler is a software system that applies compiler-engineering design patterns — staged transformation, typed intermediate representations, composable passes, and type-based verification — to the formalizable stages of gene-to-protein processing, specifically splicing and translation. The system provides deterministic, machine-checkable guarantees for designed mRNA sequences without relying on probabilistic models. It is intended for use by synthetic biologists, circuit designers, virologists, bioinformaticians, and regulatory reviewers who require formally verified gene designs rather than heuristic approximations.

This SRS serves as the authoritative reference for all subsequent engineering artifacts: the Software Architecture Document (DOC-02), the Software Design Document (DOC-03), the Interface Control Document (DOC-04), the Software Verification & Validation Plan (DOC-05), and the Traceability Matrix (DOC-08). Every requirement herein carries a unique identifier and is traceable forward into architecture, design, interface, and test artifacts, and backward into stakeholder needs and design rationale.

### 1.2 Scope

**In Scope:**

The BioCompiler system SHALL provide the following capabilities:

- **Parsing pre-mRNA sequences against a splicing grammar** to compute the set of all possible splice isoforms, using a non-deterministic finite-state transducer (NDFST) that captures alternative splicing as non-deterministic branching without assigning probabilities to isoforms.
- **Translating spliced mRNA to amino acid sequences** via the standard genetic code, implemented as a deterministic finite-state transducer (FST), with support for selenocysteine insertion, pyrrolysine incorporation, and detection of programmed ribosomal frameshifting.
- **Type-checking mRNA sequences against biological correctness properties** — including splicing correctness, codon adaptation, GC content compliance, restriction site absence, reading frame consistency, and instability motif absence — producing three-valued verdicts (PASS / FAIL / UNCERTAIN) for each property.
- **Constraint-satisfaction-based codon optimization** that eliminates cryptic splice sites while maximizing codon adaptation index (CAI), using a CSP solver that reports INFEASIBLE with a minimal unsatisfiable subset (MUS) when no feasible assignment exists.
- **Compositional verification of multi-gene circuits** encompassing promoter conflict detection, resource competition analysis, splicing interference detection, and RNA-RNA interaction screening, with three-valued verdicts composed via a defined algebra.
- **Generation of machine-checkable guarantee certificates** (proof-carrying gene designs) that are independently verifiable by a separate checker program without access to the BioCompiler pipeline.
- **Grammar-guided enumeration of legal multi-mutation combinations**, decomposing mutation space by splicing grammar nonterminals and exploiting independence across exons, with detection of constraint conflicts between individually legal but jointly illegal mutations.
- **Analysis of overlapping reading frames** with shared constraint set computation, high-coupling versus low-coupling position classification, and detection of constraint conflicts between frames.

**Out of Scope:**

The following are explicitly excluded from the BioCompiler system scope:

- **Protein structure prediction**: Handled exclusively as a foreign function interface (FFI) call to external tools (AlphaFold2, AlphaFold3, ColabFold, RoseTTAFold). The system SHALL NOT contain internal models for protein folding.
- **Post-translational modification (PTM) prediction**: Handled exclusively as an FFI call to external tools (NetPhos, PhosphoSitePlus, dbPTM, MusiteDeep). The system SHALL NOT contain internal models for PTM prediction.
- **Cellular dynamics simulation**: The system does not model transcription, translation kinetics, cellular growth, or population dynamics.
- **Probabilistic predictions of any kind**: No internal pipeline stage SHALL produce probability estimates, confidence intervals, or calibrated scores. The system uses three-valued logic (PASS/FAIL/UNCERTAIN), not probabilistic inference.
- **Grammar induction from biological sequence data**: All splicing grammar rules SHALL be specified from known biological knowledge. The system SHALL NOT learn grammar rules from data.
- **Experimental laboratory validation**: The system produces computational designs and formal guarantees; wet-lab validation of those designs is outside the system's scope.

### 1.3 Definitions, Acronyms, Abbreviations

| Term / Acronym | Definition |
|---|---|
| **CAI** | Codon Adaptation Index — a deterministic metric measuring the similarity of codon usage in a gene to the codon usage of highly expressed genes in a target organism. Values range from 0 to 1, with higher values indicating better adaptation. |
| **CSP** | Constraint Satisfaction Problem — a mathematical problem defined by a set of variables, each with a domain of possible values, and a set of constraints that restrict the values the variables can simultaneously take. The solver finds assignments satisfying all constraints or proves that no such assignment exists. |
| **DFA** | Deterministic Finite Automaton — a finite-state machine that accepts or rejects a given string of symbols by running through a uniquely determined sequence of states. Used in the BioCompiler scanner for motif detection. |
| **FFI** | Foreign Function Interface — a mechanism by which the BioCompiler system invokes external tools (folding predictors, PTM predictors) through a defined adapter contract, treating the external tool as a black box with non-deterministic output. |
| **FST** | Finite-State Transducer — a finite automaton that produces an output string for each input string, implementing a mapping from input sequences to output sequences. The translation engine uses a deterministic FST. |
| **IR** | Intermediate Representation — a typed, schema-enforced data structure (defined in protocol buffers) that mediates data flow between pipeline stages. Each IR level (IR-Seq, IR-Peptide, IR-Structure, IR-Circuit) has defined invariants. |
| **MUS** | Minimal Unsatisfiable Subset — the smallest subset of constraints from an infeasible CSP that remains unsatisfiable when considered in isolation. Used to diagnose infeasibility by pinpointing the exact set of conflicting constraints. |
| **NDFST** | Non-Deterministic Finite-State Transducer — a finite-state transducer that maps each input string to a set of possible output strings (not a probability distribution). The set-valued computation is deterministic: the same input always produces the same set of outputs. Used for splicing grammar parsing. |
| **ORF** | Open Reading Frame — a contiguous sequence of codons from a start codon (AUG) to a stop codon (UAA/UAG/UGA) that potentially encodes a protein. |
| **pLDDT** | Predicted Local Distance Difference Test — a per-residue confidence score (0–100 scale) produced by AlphaFold and related structure prediction tools. Higher scores indicate higher confidence in the predicted structure. |
| **PTM** | Post-Translational Modification — a chemical modification applied to a protein after translation, such as phosphorylation, glycosylation, acetylation, or ubiquitination. Predicted by external tools via the FFI. |
| **PWM** | Position Weight Matrix — a matrix representing the conservation of each nucleotide at each position of a sequence motif, derived from observed frequencies in known instances. Used for splice site scoring. |
| **SLOT** | A field in the Intermediate Representation (IR) that is initially empty and is filled by a specific FFI adapter. SLOT fields enable the separation of deterministic core computation from non-deterministic external tool output. |
| **ESE** | Exonic Splicing Enhancer — a short nucleotide motif within an exon that promotes the inclusion of that exon in the mature mRNA by recruiting splicing activator proteins (SR proteins). |
| **ESS** | Exonic Splicing Silencer — a short nucleotide motif within an exon that suppresses the inclusion of that exon in the mature mRNA by recruiting splicing repressor proteins (hnRNPs). |
| **ISE** | Intronic Splicing Enhancer — a short nucleotide motif within an intron that promotes the inclusion of the adjacent exon in the mature mRNA. |
| **ISS** | Intronic Splicing Silencer — a short nucleotide motif within an intron that suppresses the inclusion of the adjacent exon in the mature mRNA. |
| **SECIS** | Selenocysteine Insertion Sequence — a stem-loop structure in the 3' UTR of selenoprotein mRNAs that directs the ribosome to recode a UGA stop codon as selenocysteine (Sec, U) instead of terminating translation. |
| **UTR** | Untranslated Region — the portions of mRNA before the start codon (5' UTR) and after the stop codon (3' UTR) that are not translated into protein but contain regulatory elements. |
| **GFF3** | Generic Feature Format version 3 — a standard file format (defined by the Sequence Ontology project) for representing genomic features and annotations, used as the source format for GENCODE splice site data. |
| **MSA** | Multiple Sequence Alignment — a sequence alignment of three or more biological sequences (protein or nucleic acid), used as optional input for structure prediction FFI. |
| **SMT** | Satisfiability Modulo Theories — a decision problem for logical formulas with respect to combinations of background theories (e.g., linear arithmetic, bit-vectors). Z3 is an SMT solver used for the CSP optimizer. |
| **BDD** | Binary Decision Diagram — a data structure used to represent and manipulate Boolean functions compactly. Used for symbolic set representation of large isoform sets in the splicing grammar engine. |

### 1.4 References

| Ref ID | Document / Source | Description |
|---|---|---|
| **REF-01** | DOC-02: Software Architecture Document (SAD) | Defines the component decomposition, data flow, and architectural decisions implementing the requirements in this SRS. |
| **REF-02** | DOC-03: Software Design Document (SDD) | Provides detailed algorithmic specifications, data structures, and invariants for each component. |
| **REF-03** | DOC-04: Interface Control Document (ICD) | Defines the exact contracts for every IR schema, API, and FFI boundary. |
| **REF-04** | DOC-05: Software Verification & Validation Plan (SVVP) | Specifies how each requirement in this SRS is verified and validated. |
| **REF-05** | DOC-06: Design Rationale (DR) | Justifies each design decision and documents alternatives considered and rejected. |
| **REF-06** | DOC-07: Project Plan (PP) | Work breakdown, schedule, risk register, and staffing for system development. |
| **REF-07** | DOC-08: Traceability Matrix (TM) | Complete bidirectional traceability: requirements ↔ architecture ↔ design ↔ interface ↔ test ↔ risk. |
| **REF-08** | DOC-09: Critical Analysis of Original Framework | Identifies fatal flaws in the original concept and extracts salvageable value. |
| **REF-09** | DOC-10: Deterministic Methods for Non-Deterministic Biology | Formal methods that produce deterministic answers from non-deterministic biological systems. |
| **REF-10** | GENCODE v44, human genome annotation (GRCh38) | Source of splice site annotations, exon-intron boundaries, and gene models for splicing grammar construction. |
| **REF-11** | Codon Usage Database (Kazusa) | Source of organism-specific codon usage tables for CAI computation and codon optimization. |
| **REF-12** | REBASE, restriction enzyme database | Source of restriction enzyme recognition site sequences for the scanner and constraint solver. |
| **REF-13** | MaxEntScan splice site scoring model | Reference model for splice site strength scoring; used as a validation benchmark. |
| **REF-14** | Cousot & Cousot, "Abstract Interpretation," POPL 1977 | Foundational paper on abstract interpretation: sound over-approximation of program behavior for deterministic analysis. |
| **REF-15** | Milner, "A Theory of Type Polymorphism in Programming," JCSS 1978 | Foundational paper on type systems: "Well-typed programs don't go wrong" — the theoretical basis for the BioCompiler type system. |
| **REF-16** | Angluin, "Learning Regular Sets from Queries and Counterexamples," Information and Computation, 1987 | The L* algorithm for grammar induction; referenced as a rejected alternative (grammar induction is out of scope per REQ-CON-002). |
| **REF-17** | IEEE 830-1993 | IEEE Standard for Software Requirements Specifications — the governing standard for this document. |
| **REF-18** | ISO/IEC/IEEE 29148:2018 | Systems and software engineering — Life cycle processes — Requirements engineering — the updated international standard supplanting IEEE 830. |

### 1.5 Overview of SRS Structure

This SRS is organized as follows:

- **Section 1 (Introduction)** establishes the purpose, scope, terminology, references, and structure of this document.
- **Section 2 (Overall Description)** provides a high-level view of the product, including its operational context (product perspective), the major functional areas (product functions), the intended user classes, constraints, assumptions and dependencies, and detailed use cases.
- **Section 3 (Functional Requirements)** specifies, in complete and testable terms, every functional capability the system SHALL provide. Requirements are organized by pipeline stage: lexical analysis (Section 3.1), splicing grammar (Section 3.2), translation (Section 3.3), type system (Section 3.4), constraint-based optimization (Section 3.5), foreign function interface (Section 3.6), guarantee certificate generation (Section 3.7), compositional verification (Section 3.8), mutation exploration (Section 3.9), and overlapping reading frame analysis (Section 3.10). Each requirement has a unique identifier (REQ-FUNC-XXX).
- **Section 4 (Non-Functional Requirements)** specifies quality attributes: performance (Section 4.1), reliability (Section 4.2), usability (Section 4.3), maintainability (Section 4.4), portability (Section 4.5), and security (Section 4.6). Each requirement has a unique identifier (REQ-NFR-XXX).
- **Section 5 (Constraints)** specifies design constraints (Section 5.1) and environmental constraints (Section 5.2) that bound the solution space. Each constraint has a unique identifier (REQ-CON-XXX).
- **Section 6 (Requirement Prioritization)** maps each functional requirement to a MoSCoW priority (Must / Should / Could / Won't) to guide development phasing.
- **Section 7 (Traceability Matrix)** provides a summary mapping of every requirement to its architectural component, interface, test case, and associated risk, with a cross-reference to the full traceability in DOC-08.

---

## 2. Overall Description

### 2.1 Product Perspective

BioCompiler is a standalone software system with a pipeline architecture. It uses the compiler design pattern — staged transformation with typed intermediate representations, composable passes, and type-based verification — as an engineering methodology for building bioinformatics tools with formal verification capabilities. It is critically important to understand that the system does NOT claim that biology implements compilation. The compiler metaphor is a design pattern that provides a proven architectural template (separation of concerns, composability, formal verification at each stage), not a theoretical claim about biological information processing.

The system operates as a batch pipeline: it accepts a specification (target protein sequence plus design constraints and cellular context), executes a sequence of deterministic passes (scanning, parsing, translation, type checking, optimization), and produces one of two outcomes: (a) a verified design accompanied by a machine-checkable guarantee certificate, or (b) a failure report with a precise diagnosis identifying which constraints were violated or which knowledge gaps prevented resolution.

External tools for protein structure prediction (AlphaFold2/3, ColabFold, RoseTTAFold) and post-translational modification prediction (NetPhos, PhosphoSitePlus, dbPTM, MusiteDeep) are accessed through a Foreign Function Interface (FFI). The FFI boundary is strictly enforced: the deterministic core pipeline treats every FFI invocation as a non-deterministic black-box operation, and no guarantee certificate depends on FFI output for its core validity. FFI stages enrich the IR with optional SLOT fields that the type system may reference but does not require for its fundamental correctness guarantees.

The following context diagram illustrates the system's operational relationships:

```
                              +-------------------+
                              |   User / Client   |
                              |  (CLI or Python   |
                              |       API)        |
                              +--------+----------+
                                       |
                        Input: Target protein, constraints,
                               cellular context, grammar rules
                                       |
                                       v
                    +--------------------------------------+
                    |                                      |
                    |        BioCompiler Pipeline          |
                    |                                      |
                    |  +----------+   +-----------+        |
                    |  | COMP-01  |   | COMP-02   |        |
                    |  | Scanner  |-->| Splicing  |        |
                    |  +----------+   | Engine    |        |
                    |                 +-----------+        |
                    |                      |               |
                    |                      v               |
                    |                 +------------+       |
                    |                 | COMP-03    |       |
                    |                 | Translation|       |
                    |                 | Engine     |       |
                    |                 +-----+------+       |
                    |                       |              |
                    |          +------------+------------+ |
                    |          |                         | |
                    |          v                         v |
                    |   +-----------+            +---------+|
                    |   | COMP-04   |            | COMP-05  ||
                    |   | FFI Mgr   |            | Type     ||
                    |   | (folding, |            | System   ||
                    |   |  PTM)     |            |          ||
                    |   +-----------+            +---------+|
                    |          |                      |     |
                    |          |              +-------+----+|
                    |          |              |            ||
                    |          |              v            v ||
                    |          |       +----------+  +-------+|
                    |          |       | COMP-06  |  |COMP-07||
                    |          |       |Optimizer |  |Cert   ||
                    |          |       | (CSP)    |  |Gen    ||
                    |          |       +----------+  +-------+|
                    |          |                           |  |
                    |          |                           v  |
                    +----------|------------------------------+
                               |
                  Output: Verified design + guarantee certificate
                  OR:     Failure report + diagnosis
                               |
                    +----------+----------+
                    |  External Tools     |
                    |  (via FFI)          |
                    |                     |
                    |  AlphaFold2/3       |
                    |  ColabFold          |
                    |  RoseTTAFold        |
                    |  NetPhos            |
                    |  PhosphoSitePlus    |
                    |  dbPTM              |
                    |  MusiteDeep         |
                    +---------------------+
```

### 2.2 Product Functions (Summary)

The BioCompiler system provides the following eight major functional areas:

1. **Splicing Grammar Engine (COMP-02)**: Parse pre-mRNA sequences against a parameterized splicing grammar — implemented as a non-deterministic finite-state transducer (NDFST) — to compute the set of all possible splice isoforms consistent with known splice site consensus sequences, branch point motifs, polypyrimidine tracts, exon length constraints, and exonic/intronic splicing regulatory elements (ESE/ESS/ISE/ISS). The cellular context parameter modulates regulatory element thresholds, producing different isoform sets for different cell types.

2. **Translation Engine (COMP-03)**: Map each codon in a spliced mRNA sequence to its corresponding amino acid via a deterministic finite-state transducer (FST) implementing the standard genetic code. The engine handles all 61 sense codons, selenocysteine insertion (UGA recoding with SECIS element), pyrrolysine incorporation (UAG recoding in archaeal contexts), and detection of programmed ribosomal frameshifting motifs.

3. **Type System (COMP-05)**: Check mRNA sequences against declared biological correctness properties — SpliceCorrect(CellType), NoCrypticSplice, CodonAdapted(Organism, threshold), GCInRange(lo, hi), NoRestrictionSite(EnzymeSet), InFrame, NoInstabilityMotif — producing three-valued verdicts (PASS / FAIL / UNCERTAIN) for each property, with derivation traces for PASS, violation identification for FAIL, and knowledge gap specification for UNCERTAIN.

4. **Constraint-Based Optimizer (COMP-06)**: Find synonymous codon assignments satisfying all hard constraints (splicing correctness, CAI threshold, GC range, no restriction sites, no instability motifs, reading frame preservation) while maximizing a scalar objective (CAI), using a CSP solver that reports INFEASIBLE with a minimal unsatisfiable subset (MUS) when no feasible assignment exists.

5. **Compositional Verifier (COMP-08)**: Check cross-component constraints in multi-gene circuits, including promoter conflict (unintentional transcription factor regulation), resource competition (ribosome demand exceeding capacity), splicing interference (cryptic splice sites in one gene's transcript affecting another's), and RNA-RNA interaction (complementary transcript regions forming dsRNA).

6. **Certificate Generator (COMP-07)**: Produce machine-checkable guarantee certificates in JSON format for designs passing all type checks. Certificates include the verified sequence, each type predicate with its verdict and derivation trace, the CSP constraint set and assignment, and provenance metadata. Circuit certificates additionally include individual gene certificates and composition check results.

7. **Mutation Explorer (COMP-09)**: Decompose the mutation space of a gene into categories based on splicing grammar nonterminals (intra-exonic, splice site, regulatory element), enumerate legal multi-mutation combinations, exploit independence across exons, and report constraint conflicts (mutations that are individually legal but jointly violate the splicing grammar).

8. **ORF Analyzer (COMP-10)**: Compute shared constraint sets for overlapping reading frames in viral genomes and compact genomes, classify nucleotide positions as high-coupling (affecting multiple proteins) or low-coupling (affecting one protein), and detect constraint conflicts between frames where the optimization target for one frame contradicts the target for another.

### 2.3 User Characteristics

| User Class | Description | Needs |
|---|---|---|
| **Synthetic Biologist** | Designs genetic constructs for heterologous expression in target organisms (mammalian, yeast, bacterial). Typically has deep domain knowledge of molecular biology but limited formal methods expertise. Uses gene design tools (IDT, Genscript) but lacks formal verification of designed sequences. | Splicing-aware gene design with formal guarantees that cryptic splice sites are absent and the intended isoform is the only one produced. Deterministic verification rather than heuristic scoring. Machine-checkable certificates for regulatory submissions. |
| **Circuit Designer** | Assembles multi-gene circuits with regulatory interactions (toggle switches, oscillators, logic gates) for synthetic biology applications. Needs to verify that circuit components do not interfere with each other at the transcriptional, translational, or post-transcriptional level. | Compositional verification of circuit correctness: promoter conflict detection, resource competition analysis, splicing interference screening, and RNA-RNA interaction detection. Composable three-valued verdicts across all circuit components. |
| **Virologist** | Works with viral genomes containing overlapping reading frames (HIV, SARS-CoV-2, hepatitis B virus) where a single nucleotide position can encode amino acids in two or three different proteins simultaneously. Needs to understand coupling between frames for rational vaccine and drug design. | Shared constraint set analysis for overlapping reading frames. Classification of positions as high-coupling versus low-coupling. Detection of constraint conflicts between frames that single-frame analysis would miss. |
| **Bioinformatician** | Integrates specialized tools into custom analysis pipelines, often combining multiple software packages with glue scripts. Values well-defined interfaces, typed data formats, and programmatic access over graphical interfaces. | Typed Intermediate Representation (IR) as a universal interchange format between pipeline stages. Python API for programmatic access. Protocol buffer schemas for schema enforcement and backward-compatible evolution. |
| **Regulatory Reviewer** | Evaluates the safety of gene therapy constructs, mRNA vaccines, and other synthetic nucleic acid products for regulatory agencies (FDA, EMA). Requires auditable, verifiable evidence of safety properties rather than heuristic scores. | Machine-checkable guarantee certificates that can be independently verified by a separate checker program. Certificates that include derivation traces, constraint sets, and provenance metadata, enabling full audit trails. |

### 2.4 Constraints

Detailed constraints are specified in Section 5. In summary:

- **Design constraints**: No probabilistic models for internal stages (REQ-CON-001); no grammar induction (REQ-CON-002); no internal folding/PTM models (REQ-CON-003); no claim that biology implements compilation (REQ-CON-004).
- **Environmental constraints**: Memory limits of ≤ 32 GB RAM for single-gene analysis and ≤ 64 GB for circuit-level analysis (REQ-CON-010); GPU required only for FFI stages, not for core pipeline (REQ-CON-011).

### 2.5 Assumptions and Dependencies

**Assumptions:**

- **ASSUME-01**: The standard genetic code is correct and complete for all supported organisms. The 64 codon-to-amino-acid mappings (including the three stop codons) are treated as ground truth; the system does not validate the genetic code against experimental data.
- **ASSUME-02**: Splice site consensus sequences (the GT-AG rule for canonical splicing, the GC-AG and AT-AC rules for non-canonical splicing, branch point motifs, and polypyrimidine tract requirements) derived from GENCODE annotations are biologically valid for the target organism. The system's isoform set computation is only as complete as the grammar rules it uses; unknown or unannotated splice sites will not appear in the output.
- **ASSUME-03**: Codon usage tables from the Codon Usage Database (Kazusa) are representative of the target expression system's codon usage in highly expressed genes. CAI computations are only as meaningful as the underlying codon usage data.
- **ASSUME-04**: External tools (AlphaFold, ColabFold, NetPhos, PhosphoSitePlus, dbPTM, MusiteDeep) produce correct output when given valid input. The system does not verify external tool output; it validates only that the output conforms to the expected schema and invariants. FFI output is treated as non-deterministic (REQ-FUNC-053).
- **ASSUME-05**: The user-specified cellular context (cell type, organism) accurately describes the target expression environment. The system's splicing grammar produces different isoform sets for different cellular contexts (REQ-FUNC-013), but if the user specifies the wrong cell type, the computed isoform set will not match the biological reality.

**Dependencies:**

- **DEPEND-01**: Python 3.10+ runtime. The system is implemented in Python and requires Python 3.10 or later for type hint features (union types with `|`, `match` statements) and standard library modules.
- **DEPEND-02**: Protocol buffer compiler (`protoc`) v3.x+ for IR schema compilation. All IR schemas are defined in `.proto` files; the Python bindings are generated at build time.
- **DEPEND-03**: SAT/SMT solver library (Z3 or equivalent) for constraint satisfaction. The CSP optimizer (COMP-06) requires a solver capable of handling integer and real-valued constraints with optimization objectives.
- **DEPEND-04**: External tool installations for FFI stages. AlphaFold2/3 or ColabFold must be installed for folding FFI; NetPhos or equivalent must be installed for PTM FFI. These are optional — the core pipeline operates without them (REQ-NFR-042).
- **DEPEND-05**: GENCODE annotation files (GFF3 format) for splice site PWM construction. The splicing grammar PWMs are derived from annotated splice sites in GENCODE; without these files, the grammar cannot be constructed for the target organism.

### 2.6 Use Cases

#### UC-01: Design a Splicing-Safe Gene for Mammalian Expression

| Field | Value |
|---|---|
| **Use Case ID** | UC-01 |
| **Name** | Design a Splicing-Safe Gene for Mammalian Expression |
| **Primary Actor** | Synthetic Biologist |
| **Goal** | Produce an optimized mRNA sequence for a target protein that is guaranteed to splice correctly, free of cryptic splice sites, codon-adapted, and accompanied by a machine-checkable guarantee certificate. |
| **Preconditions** | (1) User has a target protein sequence. (2) Target organism and cell type are specified. (3) Constraint specification (CAI threshold, GC range, restriction enzyme avoidance set) is available. (4) Splicing grammar rules and codon usage tables for the target organism are installed. |
| **Postconditions** | (1) An optimized mRNA sequence is produced that translates to the target protein. (2) All type predicates return PASS. (3) A guarantee certificate is generated and passes standalone verification. |
| **Failure Postconditions** | (1) If INFEASIBLE: a failure report with MUS identifying the minimal conflicting constraint set is produced. (2) If UNCERTAIN: a report with knowledge gaps and required additional information is produced. |

**Main Success Scenario:**

| Step | Action | Component |
|---|---|---|
| 1 | User provides: target protein sequence (FASTA), target organism (e.g., *Homo sapiens*), cell type (e.g., HEK293T), and constraint specification (CAI ≥ 0.8, GC ∈ [40%, 60%], restriction enzymes to avoid: {EcoRI, BamHI, XhoI, HindIII, NotI}). | CLI / API |
| 2 | Scanner (COMP-01) tokenizes the initial codon assignment (generated from the protein sequence using a default codon table), annotating start codons, stop codons, splice donor/acceptor motifs, branch points, polypyrimidine tracts, Kozak consensus, RNA instability motifs, and restriction enzyme sites. | COMP-01 |
| 3 | Splicing Engine (COMP-02) runs the NDFST against the token stream under the specified cellular context, computing the set of all possible splice isoforms. | COMP-02 |
| 4 | Translation Engine (COMP-03) translates each splice isoform via the deterministic FST, producing IR-Peptide records with amino acid sequences, codon assignments, and any selenocysteine flags or frameshift warnings. | COMP-03 |
| 5 | Type System (COMP-05) checks all declared type predicates against the mRNA. In the typical first-pass case, cryptic splice sites are detected (NoCrypticSplice returns FAIL), or an unintended alternative isoform is found (SpliceCorrect returns FAIL). | COMP-05 |
| 6 | Optimizer (COMP-06) searches for synonymous codon substitutions that break cryptic splice sites, maintain or improve CAI, preserve GC range, and avoid introducing new violations. The CSP solver explores the space of synonymous codon assignments. | COMP-06 |
| 7 | Type System (COMP-05) re-checks all predicates against the optimized mRNA. All predicates now return PASS. | COMP-05 |
| 8 | Certificate Generator (COMP-07) produces a guarantee certificate in JSON format containing the verified sequence, all type predicates with PASS verdicts and derivation traces, the CSP constraint set and assignment, and provenance metadata. | COMP-07 |
| 9 | User receives: optimized mRNA sequence, guarantee certificate (independently verifiable), and verification report with all PASS verdicts and derivation traces. | Output |

**Alternative Flows:**

- **3a. No valid isoforms**: If the NDFST produces an empty isoform set (no parse path satisfies the grammar), the system reports a failure with the specific grammar constraint(s) that prevent any valid parse.
- **5a. Multiple FAIL verdicts**: The type system identifies all violations simultaneously. The optimizer attempts to resolve all violations in a single CSP solve. If some violations cannot be resolved by synonymous substitution alone (e.g., a premature stop codon in the protein sequence itself), the system reports INFEASIBLE.
- **6a. INFEASIBLE**: If the CSP solver finds no feasible assignment, it reports INFEASIBLE with a minimal unsatisfiable subset (MUS) identifying the smallest set of constraints that cannot be simultaneously satisfied.
- **6b. UNCERTAIN**: If some type predicates return UNCERTAIN (e.g., the splicing grammar cannot determine whether a weak splice site is functional in the specified cell type), the optimizer proceeds with the PASS/FAIL predicates but the final certificate includes UNCERTAIN verdicts.

#### UC-02: Verify a Multi-Gene Circuit

| Field | Value |
|---|---|
| **Use Case ID** | UC-02 |
| **Name** | Verify a Multi-Gene Circuit |
| **Primary Actor** | Circuit Designer |
| **Goal** | Verify that a multi-gene circuit is free of cross-component conflicts (promoter interference, resource competition, splicing interference, RNA-RNA interaction) and produce a circuit-level guarantee certificate. |
| **Preconditions** | (1) User has a circuit specification with multiple genes, their promoters, terminators, and regulatory elements. (2) Each gene has been individually designed and type-checked (UC-01). (3) Circuit topology (linear or circular) and organism are specified. |
| **Postconditions** | (1) All cross-component checks produce verdicts. (2) If all verdicts are PASS, a circuit-level guarantee certificate is generated. (3) If any verdict is FAIL or UNCERTAIN, a diagnostic report identifies the specific conflict. |

**Main Success Scenario:**

| Step | Action | Component |
|---|---|---|
| 1 | User provides: circuit specification (genes, promoters, terminators, regulatory elements, topology), organism, and cellular context. | CLI / API |
| 2 | For each gene in the circuit, the single-gene pipeline (UC-01, steps 2–5) is executed, producing per-gene IR-Seq, IR-Peptide, and type-check results. If any individual gene has FAIL verdicts that cannot be resolved, the circuit verification halts with a diagnostic report. | COMP-01..05 |
| 3 | Compositional Verifier (COMP-08) executes the linker pass, checking four categories of cross-component constraints: (a) promoter conflict — no transcription factor produced by one gene unintentionally regulates another gene's promoter; (b) resource competition — total ribosome demand does not exceed estimated cellular ribosome capacity; (c) splicing interference — no cryptic splice site in one gene's transcript interferes with splicing of another; (d) RNA-RNA interaction — no complementary regions between transcripts that could form dsRNA triggering silencing. | COMP-08 |
| 4 | If all composition verdicts are PASS, the Certificate Generator (COMP-07) produces a circuit certificate containing individual gene certificates plus composition check results with evidence. | COMP-07 |
| 5 | User receives: circuit guarantee certificate (independently verifiable), with individual gene certificates and composition check results. | Output |

**Alternative Flows:**

- **2a. Individual gene failure**: If any gene in the circuit fails individual type checking, circuit verification is not attempted. The user receives a per-gene failure report and must fix individual gene issues before re-attempting circuit verification.
- **3a. Composition FAIL**: If any composition check returns FAIL, the diagnostic report identifies the specific conflict (e.g., "Gene A produces transcription factor LacI, which binds to the promoter of Gene B, creating an unintended repression loop"). The certificate is not generated.
- **3b. Composition UNCERTAIN**: If any composition check returns UNCERTAIN (e.g., RNA-RNA interaction check cannot determine whether complementary regions form stable dsRNA under cellular conditions), the circuit certificate includes UNCERTAIN verdicts with specific knowledge gaps.

#### UC-03: Analyze Overlapping Reading Frames

| Field | Value |
|---|---|
| **Use Case ID** | UC-03 |
| **Name** | Analyze Overlapping Reading Frames |
| **Primary Actor** | Virologist |
| **Goal** | Identify nucleotide positions where mutations affect multiple proteins simultaneously, classify positions by coupling strength, detect constraint conflicts between frames, and compute the shared constraint set for rational vaccine or drug design. |
| **Preconditions** | (1) User has a nucleotide sequence with multiple annotated reading frames (e.g., a viral genome segment). (2) Reading frame boundaries (start, end, frame offset) are specified for each frame. |
| **Postconditions** | (1) Shared constraint set is computed, identifying positions where mutations have multi-protein effects. (2) All positions are classified as high-coupling or low-coupling. (3) Constraint conflicts between frames are detected and reported. |

**Main Success Scenario:**

| Step | Action | Component |
|---|---|---|
| 1 | User provides: viral genome (or compact genome segment) in FASTA format, with multiple annotated reading frames specified (frame offset, start position, end position, frame name). | CLI / API |
| 2 | Translation Engine (COMP-03) constructs a separate deterministic FST per reading frame, translating each frame's codons to amino acids independently. | COMP-03 |
| 3 | ORF Analyzer (COMP-10) computes the shared constraint set: for each nucleotide position, it determines which reading frames are affected by a mutation at that position (i.e., which frames include that nucleotide in a codon). Positions affecting multiple frames are added to the shared constraint set. | COMP-10 |
| 4 | ORF Analyzer classifies each nucleotide position as high-coupling (affects the amino acid sequence in two or more reading frames) or low-coupling (affects only one reading frame, or resides in a wobble position that does not change any amino acid). | COMP-10 |
| 5 | ORF Analyzer detects constraint conflicts: positions where the optimization target for one frame (e.g., a specific synonymous codon to avoid a cryptic splice site in Frame 1) conflicts with the optimization target for another frame (e.g., the same nucleotide position is part of a codon in Frame 2 that would introduce a premature stop codon under the Frame 1 optimization). | COMP-10 |
| 6 | User receives: shared constraint set (positions and affected frames), coupling classification (high/low per position), and conflict report (specific position-level conflicts between frames with affected amino acids and constraint descriptions). | Output |

**Alternative Flows:**

- **2a. Frameshift in one frame**: If one reading frame contains a programmed frameshift, the FST flags the position and the ORF Analyzer splits the frame into pre-frameshift and post-frameshift segments for independent analysis.
- **5a. No conflicts**: If no constraint conflicts are detected (all overlapping positions allow simultaneous optimization), the report confirms that single-frame optimization is sufficient.
- **5b. Unresolvable conflicts**: If conflicts are detected that cannot be resolved by any codon assignment, the report identifies the minimal set of proteins whose optimization targets are mutually exclusive.

---

## 3. Functional Requirements

### 3.1 Lexical Analysis (Scanner)

**REQ-FUNC-001**: The system SHALL accept a DNA or RNA nucleotide sequence as input, in either FASTA format (with one or more sequence records, each preceded by a header line beginning with '>') or raw format (a bare string of IUPAC nucleotide characters). The system SHALL normalize DNA input (T → U) for internal processing when operating in RNA mode. The system SHALL reject sequences containing non-IUPAC characters (i.e., any character other than A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N) and SHALL report the position and offending character in the error message. The maximum supported sequence length SHALL be 1,000,000 nucleotides (1 Mb).

**REQ-FUNC-002**: The system SHALL scan the input sequence and annotate the position and type of each of the following biological elements, producing an ordered list of (position, element_type, match_sequence, score) tuples:

- **Start codons** (ATG / AUG): every occurrence with its reading frame (0, 1, or 2). The scanner SHALL report all ATG occurrences, not just those in annotated ORFs.
- **Stop codons** (TAA, TAG, TGA / UAA, UAG, UGA): every occurrence with its reading frame. The scanner SHALL distinguish between in-frame and out-of-frame stop codons relative to annotated ORFs.
- **Splice donor motifs**: consensus GT at the 5' intron boundary, with extended consensus (C/A)AG|GURAGU scored against a position weight matrix (PWM). The score SHALL represent the match quality against the PWM, enabling downstream thresholding.
- **Splice acceptor motifs**: consensus AG at the 3' intron boundary, with extended consensus (C/U)AG|G scored against a PWM, including upstream polypyrimidine tract quality.
- **Branch point motifs**: consensus YNYRAY, located 18–40 nucleotides upstream of the acceptor site. The scanner SHALL search the entire 18–40 nt window and report the best-scoring match.
- **Polypyrimidine tracts**: runs of C and T/U upstream of the acceptor site. The scanner SHALL report the tract length, position, and pyrimidine purity (fraction of C + T/U).
- **Kozak consensus sequences**: translation initiation context (GCCRCCATGG, where R = A or G) around start codons. The scanner SHALL score each start codon's Kozak context for initiation strength.
- **RNA instability motifs**: AUUUA motifs (AU-rich elements in 3' UTR) and U-rich elements. The scanner SHALL report the position and motif type for each occurrence.
- **Restriction enzyme recognition sites**: all occurrences of recognition sequences for enzymes in a user-specified set, looked up from the REBASE database. The scanner SHALL report the position, enzyme name, and exact recognition sequence.

**REQ-FUNC-003**: The system SHALL implement the scanner as a collection of deterministic finite automata (DFAs), one per element type listed in REQ-FUNC-002. Each DFA SHALL be compiled from a motif specification (consensus sequence, PWM, or regular expression) at initialization time, and SHALL process the input sequence in a single left-to-right pass. The output of each DFA SHALL be a list of (position, match_sequence, score) tuples; the scanner SHALL merge the outputs of all DFAs into a single ordered token stream sorted by position. The DFA-based design ensures O(n) scanning time per element type, where n is the sequence length.

**REQ-FUNC-004**: The scanner SHALL be fully deterministic: given the same input sequence and the same DFA definitions (motif specifications and scoring parameters), the output SHALL be byte-identical across runs, including identical tuple ordering, identical scores (to floating-point precision), and identical position assignments. This determinism requirement extends across platforms (Linux x86_64, Linux ARM64, macOS x86_64, macOS ARM64) within floating-point tolerance (±1 ULP for score values).

### 3.2 Splicing Grammar (Parsing)

**REQ-FUNC-010**: The system SHALL implement a non-deterministic finite-state transducer (NDFST) that takes the token stream produced by the scanner (REQ-FUNC-002) as input and produces the set of all possible splice isoforms consistent with the splicing grammar. The NDFST captures alternative splicing as non-deterministic branching: each valid parse path through the transducer corresponds to a distinct splice isoform. The output is a set-valued result — a set of mRNA sequences, not a probability distribution over sequences. The set computation SHALL be exhaustive: the NDFST SHALL enumerate all parse paths that satisfy the grammar constraints, not merely a representative sample.

**REQ-FUNC-011**: The splicing grammar SHALL encode the following rules as NDFST transitions, with each rule contributing to the set of possible isoforms:

- **Canonical splice donor consensus**: (C/A)AG|GURAGU, where '|' marks the exon-intron boundary and R = A or G. The donor site SHALL be scored against a PWM derived from GENCODE-annotated donor sites. Sites exceeding a configurable score threshold SHALL be included as potential donor sites in the NDFST.
- **Canonical splice acceptor consensus**: (C/U)AG|G, where '|' marks the intron-exon boundary. The acceptor site SHALL be scored against a PWM derived from GENCODE-annotated acceptor sites, including upstream polypyrimidine tract quality. Sites exceeding a configurable score threshold SHALL be included as potential acceptor sites.
- **Non-canonical splice sites**: GC-AG (approximately 0.5% of human introns) and AT-AC (rare, U12-type introns) SHALL be included as rare alternatives with lower confidence scores. The NDFST SHALL include transitions for non-canonical sites, but the type system (REQ-FUNC-031, NoCrypticSplice) SHALL distinguish them from canonical sites.
- **Branch point motif**: YNYRAY (where Y = C or T/U, R = A or G, N = any), located 18–40 nucleotides upstream of the acceptor site. The grammar SHALL require a branch point with a score above a minimum threshold for each intron; the absence of a viable branch point SHALL eliminate the corresponding parse path.
- **Polypyrimidine tract**: minimum length and purity thresholds for the pyrimidine-rich region upstream of the acceptor site. The grammar SHALL include a soft constraint on polypyrimidine tract quality; tracts below the minimum threshold SHALL reduce the score of the corresponding parse path but SHALL NOT eliminate it entirely (reflecting the biological reality that weak tracts reduce but do not abolish splicing).
- **Exon length constraints**: minimum and maximum exon lengths SHALL be encoded as soft constraints (penalties for outlier lengths) rather than hard constraints, reflecting the biological reality that exon length distributions have long tails.
- **Exonic and intronic splicing regulatory elements**: ESEs, ESSs, ISEs, and ISSs SHALL modulate the inclusion or exclusion of exons. The grammar SHALL include transitions that are activated or deactivated based on the presence of these regulatory elements, with their effect modulated by the cellular context parameter (REQ-FUNC-013).

**REQ-FUNC-012**: The NDFST SHALL capture alternative splicing as non-deterministic branching. At each position where the grammar permits multiple valid transitions (e.g., a potential exon that can be either included or skipped, or a competing donor/acceptor pair that produces alternative exon boundaries), the NDFST SHALL branch into multiple parse paths. Each complete parse path from the start state to an accept state produces one splice isoform. The set of all complete parse paths constitutes the set of all possible splice isoforms. This non-determinism is biological (multiple real outcomes are possible) but the computation is deterministic (the same input always produces the same set of isoforms, per REQ-FUNC-014).

**REQ-FUNC-013**: The NDFST SHALL accept a cellular context specification parameter that modulates ESE, ESS, ISE, and ISS threshold values. Different cell types express different repertoires of splicing factors (SR proteins, hnRNPs), and the same mRNA sequence can produce different isoform sets in different cellular contexts. The cellular context parameter SHALL include: cell type identifier (e.g., "HEK293T", "HepG2", "neuron"), organism identifier (e.g., "Homo sapiens", "Mus musculus"), and optional signaling state (e.g., "hypoxia", "growth_factor_stimulation"). The grammar SHALL adjust the strength of regulatory element effects based on the known expression profiles of splicing factors in the specified cell type, producing a different isoform set for different cellular contexts. The cellular context SHALL be loaded from declarative configuration files (not hardcoded), enabling updates without source code modification.

**REQ-FUNC-014**: The computation SHALL be deterministic: given the same input sequence, the same cellular context parameters, and the same grammar rules (motif specifications, PWMs, regulatory element thresholds), the NDFST SHALL produce the same set of isoforms. This means: (a) the cardinality of the isoform set SHALL be identical across runs; (b) each isoform's spliced mRNA sequence, exon boundaries, and reading frame SHALL be identical across runs; (c) the ordering of isoforms in the output SHALL be deterministic (sorted by a defined criterion, e.g., lexicographic order of the spliced sequence). Determinism is critical for reproducibility, regression testing, and certificate verification.

**REQ-FUNC-015**: The NDFST output SHALL include, for each isoform in the set: (a) the spliced mRNA sequence (concatenation of exon sequences with introns removed); (b) exon boundaries (start and end positions in the pre-mRNA, in the coordinate system of the input sequence); (c) the reading frame (0, 1, or 2) established by the start codon and maintained across exon-exon junctions; (d) all annotations inherited from the scanner token stream, re-mapped to the coordinates of the spliced sequence (e.g., a restriction site that spans an exon-exon junction in the pre-mRNA may appear as a newly created site in the spliced mRNA); and (e) a provenance record identifying the grammar rules and cellular context that produced this isoform.

### 3.3 Translation

**REQ-FUNC-020**: The system SHALL implement translation as a deterministic finite-state transducer (FST) that maps each codon in a spliced mRNA sequence to its corresponding amino acid according to the standard genetic code. The FST SHALL consume the spliced mRNA sequence (output of REQ-FUNC-015) three nucleotides at a time, in the reading frame established by the start codon, and produce an amino acid for each codon. The FST SHALL terminate upon encountering the first in-frame stop codon. The FST states SHALL encode the current reading frame offset (0, 1, or 2), enabling frame-aware translation.

**REQ-FUNC-021**: The translation FST SHALL handle the following special cases:

- **All 61 sense codons**: Each of the 61 sense codons SHALL be mapped to the corresponding amino acid per the standard genetic code. The codon-to-amino-acid mapping SHALL be loaded from a codon table file (not hardcoded), supporting alternative genetic codes (e.g., vertebrate mitochondrial, yeast mitochondrial) in future extensions.
- **Selenocysteine (Sec, U) insertion**: The UGA codon, which is normally a stop codon, SHALL be recoded as selenocysteine when a SECIS (Selenocysteine Insertion Sequence) element is present in the 3' UTR of the mRNA. The FST SHALL detect the presence of a SECIS element (identified by the scanner in REQ-FUNC-002 or provided as an annotation in the input) and, when present, SHALL flag the UGA codon as selenocysteine rather than stop, producing a `SelenocysteineFlag` annotation at the corresponding position.
- **Pyrrolysine (Pyl, O) insertion**: The UAG codon, which is normally a stop codon, SHALL be recoded as pyrrolysine in archaeal contexts (when the organism parameter indicates an archaeon with known pyrrolysine incorporation). The FST SHALL flag the UAG codon as pyrrolysine in these contexts, producing a corresponding annotation.
- **Ribosomal frameshifting**: Known programmed frameshift motifs (e.g., the HIV-1 -1 frameshift motif U UUU UUUA, the coronavirus -1 frameshift motif U UUA AAC) SHALL be detected by the scanner (REQ-FUNC-002) and flagged by the FST as annotations rather than resolved. The FST SHALL NOT attempt to simulate frameshifting; instead, it SHALL produce a `FrameshiftWarning` annotation at the motif position, indicating the direction (-1 or +1) and the motif sequence. Downstream consumers (type system, certificate generator) SHALL treat frameshift warnings as requiring manual review.

**REQ-FUNC-022**: The translation FST SHALL be deterministic: given the same input mRNA sequence, the same codon table, and the same annotation set (SECIS presence, organism context), the same amino acid output SHALL be produced. This determinism holds modulo the explicitly flagged frameshift ambiguities (REQ-FUNC-021), which are not resolved by the FST but are reported as warnings. For sequences without programmed frameshift motifs, the translation SHALL be fully deterministic with no ambiguous output.

**REQ-FUNC-023**: The translation output SHALL include: (a) the amino acid sequence (using the IUPAC protein alphabet, including U for selenocysteine and O for pyrrolysine when applicable); (b) per-position codon assignment (mapping each amino acid position to the specific codon used, enabling downstream CAI computation and codon optimization); (c) selenocysteine flags (positions where UGA was recoded as Sec, with SECIS element presence confirmation); (d) frameshift warnings (positions where a programmed frameshift motif was detected, with motif sequence and direction); and (e) initiation site confidence scores (based on Kozak consensus strength from the scanner).

### 3.4 Type System (Semantic Analysis)

**REQ-FUNC-030**: The system SHALL implement a type checker that evaluates an mRNA sequence against a set of declared biological correctness types, producing for each type a verdict of PASS, FAIL, or UNCERTAIN. The three-valued verdict system reflects the epistemic status of the system's knowledge: PASS means the property is guaranteed to hold (the system has a proof); FAIL means the property is guaranteed not to hold (the system has a counterexample); UNCERTAIN means the system cannot determine whether the property holds or not (the available knowledge is insufficient to resolve the question). The type checker SHALL NOT produce probabilities, confidence scores, or rankings — only the three-valued verdict.

**REQ-FUNC-031**: The system SHALL support 33 type predicates (13 core + 20 SLOT-dependent), each of which is a function from an mRNA sequence and its annotations to a three-valued verdict. The predicates are organized into five categories: DNA-level (12), structure (4), stability (4), solubility (4), and immunogenicity (4). The 13 core predicates are those that can be evaluated without FFI-derived SLOT fields; the 20 SLOT-dependent predicates require SLOT fields filled by FFI adapters or external oracles.

**DNA-level predicates (12):**

- `SpliceCorrect(CellType)`: The mRNA produces exactly the intended splice isoform under the specified cellular context. PASS if the NDFST output set is a singleton containing only the intended isoform; FAIL if the NDFST output set contains an isoform other than the intended one; UNCERTAIN if the NDFST output set contains the intended isoform plus additional isoforms whose production depends on cellular context factors that are not fully characterized.
- `NoCrypticSplice`: The mRNA contains no sequences matching the splice site grammar above a configurable strength threshold. PASS if no potential splice donor or acceptor site exceeds the threshold in the coding sequence; FAIL if at least one site exceeds the threshold with sufficient evidence; UNCERTAIN if weak sites are present that might or might not function as splice sites depending on the cellular context.
- `CodonAdapted(Organism, threshold)`: The Codon Adaptation Index (CAI) of the mRNA meets or exceeds the specified threshold for the target organism. PASS if CAI ≥ threshold; FAIL if CAI < threshold and the codon table is complete; UNCERTAIN if the codon usage data for the organism is incomplete or based on too few gene samples.
- `GCInRange(lo, hi)`: The GC percentage of the mRNA falls within the inclusive range [lo, hi]. PASS if lo ≤ GC% ≤ hi; FAIL if GC% < lo or GC% > hi (this predicate is never UNCERTAIN, as GC% is a deterministic computation).
- `NoRestrictionSite(EnzymeSet)`: No recognition sites for any enzyme in the specified set are present in the mRNA sequence. PASS if no sites are found; FAIL if at least one site is found (this predicate is never UNCERTAIN, as site presence is a deterministic DFA match).
- `InFrame`: The reading frame is consistent throughout the coding region, with no premature stop codons. PASS if the start codon is in-frame with the stop codon and no in-frame stop codons occur between them; FAIL if a premature in-frame stop codon is detected; UNCERTAIN if frameshift motifs are present that make frame determination ambiguous.
- `NoInstabilityMotif`: No known mRNA destabilizing motifs (AUUUA, U-rich elements) are present in the 3' UTR. PASS if no motifs are found; FAIL if at least one motif is found in the 3' UTR; UNCERTAIN if motifs are found in regions whose UTR status is uncertain (e.g., alternative polyadenylation sites).
- `NoCpGIsland`: The mRNA contains no CpG islands (regions of elevated CpG density) that could trigger epigenetic silencing. PASS if no CpG islands are detected; FAIL if a CpG island exceeding the density threshold is found; UNCERTAIN if borderline CpG density is present.
- `NoCrypticPromoter(Organism, threshold)`: The mRNA contains no cryptic promoter sequences that could drive unintended transcription initiation. PASS if no promoter-like sequences exceed the threshold; FAIL if a cryptic promoter is detected; UNCERTAIN if weak promoter signals are present.
- `NoUnexpectedTMDomain(is_cytosolic, threshold)`: A cytosolic protein has not gained unexpected transmembrane domains after codon optimization. PASS if no hydrophobic stretch exceeds the TM threshold for cytosolic proteins; FAIL if an unexpected TM domain is detected; UNCERTAIN if borderline hydrophobicity is present.
- `mRNASecondaryStructure(window_start, window_end, dg_threshold)`: The mRNA secondary structure around the RBS/start codon is not so stable as to inhibit ribosome binding. PASS if the estimated free energy is above the threshold; FAIL if excessively stable secondary structure is detected; UNCERTAIN if the simplified model cannot resolve the prediction.
- `CoTranslationalFolding(Organism, min_pause_cai)`: Codon optimization preserves the co-translational folding profile, including the codon ramp and pause sites at domain boundaries. PASS if the ramp and pause sites are preserved; FAIL if the ramp is destroyed or domain boundary pause sites are disrupted; UNCERTAIN if moderate disruption is detected.

**Structure predicates (4, SLOT-dependent):**

- `StructureConfidence`: The predicted protein structure has sufficient confidence (pLDDT) for meaningful analysis. PASS if mean pLDDT exceeds threshold; FAIL if confidence is too low for reliable analysis; UNCERTAIN if FFI folding results are unavailable (SLOT unfilled).
- `NoMisfoldingRisk`: The predicted structure shows no signs of misfolding (clashes, Ramachandran outliers, rotamer outliers). PASS if no misfolding indicators are present; FAIL if significant structural anomalies are detected; UNCERTAIN if FFI results are unavailable.
- `CorrectFoldTopology`: The predicted protein fold matches the expected topology class. PASS if the fold topology is consistent with the expected architecture; FAIL if the predicted topology is incompatible; UNCERTAIN if FFI results are unavailable or the expected topology is not specified.
- `NoUnexpectedInteraction`: The predicted structure does not contain unexpected inter-domain or inter-chain interactions. PASS if no unexpected contacts are detected; FAIL if spurious interactions are found; UNCERTAIN if FFI results are unavailable.

**Stability predicates (4, SLOT-dependent):**

- `StableFolding(stability_threshold)`: The predicted protein structure is thermodynamically stable above the specified threshold. PASS if the folding energy estimate is favorable; FAIL if the structure is predicted to be unstable; UNCERTAIN if FFI-derived energy estimates are unavailable.
- `NoDestabilizingMutation(max_ddg)`: No mutation in the designed sequence destabilizes the protein beyond the specified ΔΔG threshold. PASS if all predicted ΔΔG values are within tolerance; FAIL if a destabilizing mutation is detected; UNCERTAIN if FFI results are unavailable.
- `DisulfideBondIntegrity`: Disulfide bonds in the native structure are preserved after codon optimization (no cysteine-disrupting mutations at bond positions). PASS if all disulfide-bonded cysteines are preserved; FAIL if a disulfide bond is disrupted; UNCERTAIN if FFI results are unavailable.
- `HydrophobicCoreQuality`: The hydrophobic core of the predicted structure is well-packed with no buried polar residues or voids. PASS if core quality metrics are satisfactory; FAIL if core packing defects are detected; UNCERTAIN if FFI results are unavailable.

**Solubility predicates (4, SLOT-dependent):**

- `SolubleExpression(min_solubility_score)`: The predicted solubility of the protein meets or exceeds the specified threshold. PASS if the solubility score is above threshold; FAIL if the protein is predicted to be insoluble; UNCERTAIN if FFI-derived solubility predictions are unavailable.
- `NoAggregationProneRegion`: The protein sequence contains no aggregation-prone regions (e.g., long stretches of hydrophobic or amyloidogenic residues). PASS if no aggregation-prone regions are found; FAIL if aggregation-prone regions are detected; UNCERTAIN if FFI results are unavailable.
- `ChargeComposition`: The overall charge composition of the protein is within acceptable bounds for soluble expression. PASS if the net charge and charge distribution are favorable; FAIL if charge composition predicts insolubility; UNCERTAIN if FFI results are unavailable.
- `NoLongHydrophobicStretch`: The protein contains no excessively long hydrophobic stretches that could cause membrane insertion or aggregation. PASS if all hydrophobic stretches are within tolerance; FAIL if a problematic stretch is detected; UNCERTAIN if FFI results are unavailable.

**Immunogenicity predicates (4, SLOT-dependent):**

- `LowImmunogenicity(max_immunogenicity_score)`: The overall immunogenicity score of the protein is below the specified threshold. PASS if the score is below threshold; FAIL if the protein is predicted to be immunogenic; UNCERTAIN if FFI-derived immunogenicity predictions are unavailable.
- `NoStrongTCellEpitope(mhc_alleles)`: The protein contains no strong T-cell epitopes for the specified MHC allele set. PASS if no strong epitopes are found; FAIL if a strong T-cell epitope is detected; UNCERTAIN if FFI results or MHC binding predictions are unavailable.
- `NoDominantBCellEpitope`: The protein contains no dominant B-cell epitopes on its surface. PASS if no dominant epitopes are found; FAIL if a dominant B-cell epitope is detected; UNCERTAIN if FFI results are unavailable.
- `PopulationCoverageSafe(mhc_alleles)`: The protein's epitope profile is safe across the specified population's MHC allele distribution. PASS if population coverage is within safe bounds; FAIL if a significant subpopulation is at risk; UNCERTAIN if FFI results or population allele frequency data are unavailable.

**REQ-FUNC-032**: The type checker SHALL be deterministic: given the same input mRNA sequence, the same set of type predicates, the same grammar rules, and the same cellular context parameters, the type checker SHALL produce the same set of verdicts. This means: (a) each predicate SHALL produce the same verdict (PASS, FAIL, or UNCERTAIN) across runs; (b) the derivation traces, violation identifications, and knowledge gap descriptions SHALL be identical across runs; (c) the ordering of verdicts in the output SHALL be deterministic (matching the ordering of predicates in the input).

**REQ-FUNC-033**: The system SHALL implement the following subtyping relations, which enable hierarchical verification:

- `SpliceCorrect(SpecificCellType) <: SpliceCorrect(GeneralCellType)`: When the specific cell type is a member of the general category (e.g., HEK293T <: HumanCell <: MammalianCell), a PASS verdict for the specific type implies a PASS verdict for the general type. This enables verification at the finest available granularity with automatic propagation to coarser categories.
- `CodonAdapted(O, t1) <: CodonAdapted(O, t2)`: When t1 ≥ t2, a PASS verdict for the higher threshold implies a PASS verdict for the lower threshold. A design verified at CAI ≥ 0.9 is automatically verified at CAI ≥ 0.8.
- `GCInRange(lo1, hi1) <: GCInRange(lo2, hi2)`: When lo2 ≤ lo1 and hi1 ≤ hi2 (the stricter range is a subset of the more permissive range), a PASS verdict for the stricter range implies a PASS verdict for the more permissive range.

The subtyping relation SHALL be used in two ways: (a) to automatically promote verified designs from specific to general contexts without re-verification, and (b) to detect when a user's constraint specification contains redundant predicates (e.g., both `CodonAdapted(Human, 0.8)` and `CodonAdapted(Human, 0.9)` — the latter subsumes the former).

**REQ-FUNC-034**: The system SHALL implement three-valued composition rules for combining verdicts from multiple type predicates. These rules ensure that the overall verdict for a conjunction (all predicates must pass) or disjunction (at least one predicate must pass) of predicates is well-defined and composable without probabilistic assumptions about predicate independence:

- **Conjunction (AND)**: PASS ∧ PASS = PASS; PASS ∧ UNCERTAIN = UNCERTAIN; PASS ∧ FAIL = FAIL; UNCERTAIN ∧ UNCERTAIN = UNCERTAIN; UNCERTAIN ∧ FAIL = FAIL; FAIL ∧ FAIL = FAIL. The key principle: a single FAIL dominates (the design cannot be accepted); UNCERTAIN propagates (the overall design cannot be fully verified but is not known to be incorrect).
- **Disjunction (OR)**: PASS ∨ PASS = PASS; PASS ∨ UNCERTAIN = PASS; PASS ∨ FAIL = PASS; UNCERTAIN ∨ UNCERTAIN = UNCERTAIN; UNCERTAIN ∨ FAIL = UNCERTAIN; FAIL ∨ FAIL = FAIL. The key principle: a single PASS dominates (the design satisfies at least one required property); UNCERTAIN is better than FAIL but not as good as PASS.

These composition rules SHALL be applied when: (a) combining the verdicts of multiple type predicates for a single gene; (b) combining the verdicts of multiple composition checks in a multi-gene circuit; and (c) combining per-gene verdicts with circuit-level composition verdicts.

**REQ-FUNC-035**: For each verdict, the type checker SHALL produce evidence appropriate to the verdict type:

- **For each PASS verdict**: The system SHALL produce a derivation trace — a chain of reasoning steps, each citing a specific grammar rule, constraint, or computation that justifies the PASS verdict. The derivation trace SHALL be detailed enough to enable independent verification: given the trace, the sequence, and the rule definitions, a human or machine reviewer SHALL be able to confirm that the PASS verdict is justified. Example: for `NoCrypticSplice`, the derivation trace SHALL list every position where a splice donor or acceptor motif was detected by the scanner, the score of each such motif, and the fact that each score falls below the cryptic splice site threshold.
- **For each FAIL verdict**: The system SHALL produce a specific violation identification including: (a) the position in the sequence where the violation occurs; (b) the rule or constraint that is violated; (c) the evidence for the violation (e.g., the exact cryptic splice site sequence, the CAI value and the threshold it fails to meet, the restriction enzyme site found). Example: for a `NoCrypticSplice` FAIL, the violation identification SHALL include the position of the cryptic donor site, the donor consensus match, the PWM score, and the threshold that was exceeded.
- **For each UNCERTAIN verdict**: The system SHALL identify the specific knowledge gap that prevents resolution: (a) what additional information is needed to convert the verdict to PASS or FAIL; (b) what additional constraints or data the user could provide to resolve the uncertainty; (c) which cellular context factors are uncharacterized. Example: for an UNCERTAIN `SpliceCorrect(HEK293T)`, the knowledge gap description SHALL identify the specific alternative isoform(s) that might or might not be produced, the weak splice site(s) whose functionality depends on unknown splicing factor concentrations, and what experimental data (e.g., RT-PCR under the specified cellular conditions) would resolve the verdict.

### 3.5 Constraint-Based Optimization

**REQ-FUNC-040**: The system SHALL implement a constraint satisfaction problem (CSP) solver that finds synonymous codon assignments satisfying all hard constraints while maximizing a scalar objective function (by default, the Codon Adaptation Index, CAI). The CSP formulation is as follows: given a target protein of length L amino acids, there are L decision variables (one per codon position); the domain of each variable is the set of synonymous codons for the amino acid at that position (1–6 codons depending on the amino acid). The solver SHALL search for an assignment of codons to positions such that all hard constraints are simultaneously satisfied, and among all such feasible assignments, return one (or more) that maximizes the objective.

**REQ-FUNC-041**: Decision variables SHALL be defined as one per codon position in the coding region of the mRNA. The domain of each decision variable SHALL be the set of synonymous codons for the amino acid at that position, as defined by the standard genetic code. For amino acids with a single codon (Met: AUG; Trp: UGG), the domain is a singleton and the variable is fixed. For amino acids with multiple codons (e.g., Leu: 6 codons; Ser: 6 codons; Arg: 6 codons), the domain size ranges from 2 to 6. The total search space is the Cartesian product of all domains, which can be as large as 6^L for a protein of length L (though constraint propagation dramatically reduces this in practice).

**REQ-FUNC-042**: The following hard constraints SHALL be supported by the CSP solver. A hard constraint is a requirement that MUST be satisfied by any feasible solution; violations are not tolerated:

- **Splicing correctness**: No synonymous substitution SHALL create a cryptic splice site (a sequence matching the splice site grammar above the NoCrypticSplice threshold). This constraint is implemented by checking, for each candidate codon at each position, whether the resulting local sequence (the codon plus its flanking context) creates a new donor or acceptor motif above the threshold.
- **CAI threshold**: The overall Codon Adaptation Index of the resulting mRNA SHALL meet or exceed the specified threshold value. This constraint is global (it depends on all codon choices simultaneously), not decomposable by position.
- **GC content range**: The overall GC percentage of the resulting mRNA SHALL fall within the inclusive range [lo, hi]. This constraint is global (it depends on all codon choices simultaneously).
- **No restriction sites**: No synonymous substitution SHALL introduce a restriction enzyme recognition site for any enzyme in the specified avoidance set. This constraint is local (it depends on the codon choice and its flanking context), but a single codon change can create a restriction site that spans the boundary between two codon positions.
- **No instability motifs**: No synonymous substitution SHALL introduce a known mRNA destabilizing motif (AUUUA, U-rich element) in the 3' UTR. This constraint is local to the 3' UTR region.
- **Reading frame preservation**: This constraint is satisfied by construction — all substitutions are synonymous (same amino acid, different codon), so the reading frame and protein sequence are preserved by definition. The constraint is included in the specification for completeness and to make explicit that the optimizer SHALL NOT consider non-synonymous mutations.

**REQ-FUNC-043**: When the CSP is feasible (at least one assignment satisfies all constraints), the solver SHALL return at least one such assignment. The returned assignment SHALL include: (a) the complete codon assignment for every position; (b) the objective value (CAI) of the assignment; (c) a verification record confirming that every hard constraint is satisfied (enabling independent verification). When multiple optimal or near-optimal solutions exist, the solver SHALL return either all solutions (when the number is small, configurable threshold ≤ 10) or a representative sample (when the number is large), sorted by objective value in descending order. The solver SHALL be deterministic in its selection: given the same problem, it SHALL return the same set of solutions in the same order.

**REQ-FUNC-044**: When the CSP is infeasible (no assignment satisfies all constraints), the solver SHALL report INFEASIBLE and provide a minimal unsatisfiable subset (MUS) — the smallest set of constraints that cannot be simultaneously satisfied. The MUS SHALL identify the precise source of infeasibility, enabling the user to relax or remove specific constraints. The MUS computation SHALL satisfy the following properties: (a) the MUS is unsatisfiable (no assignment satisfies all constraints in the MUS simultaneously); (b) the MUS is minimal (removing any single constraint from the MUS makes the remaining set satisfiable); (c) the MUS is reported with human-readable descriptions of each constraint in the conflict set (e.g., "CAI ≥ 0.9 for Homo sapiens" conflicts with "NoCrypticSplice at position 347" because all synonymous codons at position 347 that would break the cryptic splice site have low CAI scores). The solver MAY also provide relaxation suggestions: specific constraint modifications that would make the problem feasible.

**REQ-FUNC-045**: The optimizer SHALL be deterministic: given the same inputs (protein sequence, constraint set, objective function, codon tables, grammar rules, cellular context) and the same solver parameters (timeout, search strategy), the optimizer SHALL produce the same output. For feasible problems, the same assignment(s) SHALL be returned in the same order. For infeasible problems, the same MUS SHALL be reported. This determinism requirement is essential for reproducibility, regression testing, and certificate verification. The optimizer SHALL NOT use randomized search strategies (e.g., random restarts, stochastic local search) unless the random seed is deterministically derived from the input.

### 3.6 Foreign Function Interface (FFI)

**REQ-FUNC-050**: The system SHALL invoke external protein structure prediction tools (AlphaFold2, AlphaFold3, ColabFold, RoseTTAFold) through a defined FFI adapter interface. The FFI invocation SHALL provide: (a) the amino acid sequence (from IR-Peptide); (b) an optional multiple sequence alignment (MSA) in A3M format; and (c) adapter-specific configuration parameters. The FFI invocation SHALL consume: (a) 3D atomic coordinates (at minimum, C-alpha coordinates for every residue); (b) per-residue confidence scores (pLDDT, scale 0–100); (c) predicted aligned error (PAE) matrix; and (d) any validation flags from the external tool. The FFI adapter SHALL parse the external tool's output into the IR-Structure schema, validating that the output satisfies the schema invariants (INV-STR-01, INV-STR-02). If the output fails validation, the adapter SHALL raise an `OutputParseError` with a description of the invariant violation.

**REQ-FUNC-051**: The system SHALL invoke external PTM prediction tools (NetPhos, PhosphoSitePlus, dbPTM, MusiteDeep) through a defined FFI adapter interface. The FFI invocation SHALL provide: (a) the amino acid sequence (from IR-Peptide); (b) the cellular context specification (cell type, organism, signaling state); and (c) adapter-specific configuration parameters. The FFI invocation SHALL consume: (a) predicted PTM sites with modification type (phosphorylation, glycosylation, acetylation, ubiquitination, etc.); (b) confidence scores for each predicted site; and (c) contextual information (kinase specificity, sequence motif). The FFI adapter SHALL parse the external tool's output into the IR-Peptide PTM SLOT fields, validating that each PTM site references a valid position in the amino acid sequence.

**REQ-FUNC-052**: The system SHALL NOT model the internal computation of any external tool. The guarantee provided by the FFI is limited to three properties: (a) correct input formatting — the system SHALL construct the input to the external tool in the format expected by the tool's API or command-line interface; (b) correct output parsing — the system SHALL parse the external tool's output into the appropriate IR schema, validating schema invariants; and (c) provenance metadata preservation — the system SHALL record the external tool's name, version, invocation parameters, and output checksum in the IR provenance field, enabling full audit trails. The system SHALL NOT guarantee the correctness of the external tool's predictions, the consistency of the external tool's output across invocations, or the scientific validity of the external tool's models.

**REQ-FUNC-053**: FFI invocations SHALL be treated as non-deterministic: the system SHALL NOT assume that the same output is produced for the same input across runs of an external tool. This non-determinism has the following implications: (a) the core pipeline SHALL NOT depend on FFI output for its fundamental correctness guarantees (type-check verdicts and certificate generation SHALL be valid without FFI output); (b) when FFI output is present, the type system MAY reference it for additional checks (e.g., structure-aware validation), but these checks SHALL be clearly separated from the core type predicates and SHALL be labeled as FFI-dependent in the verification report; (c) the certificate SHALL clearly distinguish between guarantees derived from deterministic internal computation and annotations derived from non-deterministic FFI output; (d) reproducibility testing (REQ-NFR-010) SHALL NOT apply to FFI stages.

### 3.7 Guarantee Certificate Generation

**REQ-FUNC-060**: For every design that passes all type checks (all verdicts = PASS, or all verdicts ∈ {PASS, UNCERTAIN} with explicit user acceptance of UNCERTAIN verdicts), the system SHALL generate a guarantee certificate in JSON format containing the following fields:

- **version**: Certificate schema version identifier (e.g., "1.0.0").
- **design_id**: SHA-256 hash of the verified nucleotide sequence, providing a tamper-evident binding between the certificate and the sequence it certifies.
- **sequence**: The verified nucleotide sequence (the exact sequence to which the guarantee applies).
- **types**: An array of type predicate results, each containing: (a) the predicate name and parameters; (b) the verdict (PASS or UNCERTAIN; FAIL verdicts SHALL prevent certificate generation); (c) the derivation trace (for PASS verdicts) or the knowledge gap description (for UNCERTAIN verdicts).
- **constraints**: The CSP constraint set and the assignment found, enabling independent verification that the optimized codon assignment satisfies all stated constraints.
- **provenance**: Metadata including: (a) tool name and version; (b) timestamp (ISO 8601); (c) all input parameters; (d) SHA-256 hash of the input protein sequence and constraint specification; (e) version identifiers for all grammar rules, codon tables, and PWM files used; (f) FFI adapter outputs referenced (if any), with their own provenance metadata.

**REQ-FUNC-061**: The guarantee certificate SHALL be independently verifiable: a separate checker program, given only the certificate JSON file and the system's rule definitions (grammar rules, codon tables, PWMs, restriction enzyme sequences), SHALL be able to confirm that every PASS verdict is justified by re-executing the relevant type checks and verifying that the derivation traces are consistent with the re-executed checks. The independent checker SHALL NOT require the BioCompiler pipeline or any external tools; it SHALL be a standalone program that reads the certificate, re-derives the expected verdicts, and confirms or denies the certificate's claims. The independent checker SHALL produce a clear verification result: CERTIFICATE_VALID or CERTIFICATE_INVALID with specific reasons for any invalidity.

**REQ-FUNC-062**: For compositional verification (multi-gene circuits), the certificate SHALL additionally contain:

- **individual_gene_certificates**: One certificate per gene in the circuit (as specified in REQ-FUNC-060), each independently verifiable.
- **composition_checks**: An array of composition check results, each containing: (a) the check type (promoter_conflict, resource_competition, splicing_interference, rna_interaction); (b) the genes involved; (c) the verdict (PASS, FAIL, or UNCERTAIN); (d) the evidence supporting the verdict.
- **circuit_context**: The cellular context specification under which the composition was verified, including organism, cell type, and signaling state.
- **circuit_topology**: Whether the construct is linear or circular, and the gene ordering.

The circuit certificate as a whole is valid if and only if: (a) all individual gene certificates are valid; (b) all composition check verdicts are PASS (or UNCERTAIN with explicit user acceptance); and (c) the circuit context is specified and non-empty.

### 3.8 Compositional Verification (Multi-Gene Circuits)

**REQ-FUNC-070**: The system SHALL accept a circuit-level specification containing multiple genes with their promoters, terminators, and regulatory elements. The circuit specification SHALL include: (a) the nucleotide sequence for each gene (in FASTA or raw format); (b) the promoter associated with each gene (type: constitutive, inducible, repressible; strength; transcription factor activators and repressors); (c) the terminator associated with each gene (type, efficiency); (d) the construct topology (linear or circular); (e) the organism and cellular context; and (f) the relative positioning of genes in the construct (enabling adjacency-based interference analysis).

**REQ-FUNC-071**: The system SHALL implement a linker pass that checks the following compositional constraints across genes in the circuit:

- **Promoter conflict**: No transcription factor produced by one gene SHALL unintentionally regulate another gene's promoter. The checker SHALL compare the list of transcription factors produced by each gene (identified by the gene's protein product annotation) against the list of transcription factor activators and repressors for each promoter. If a match is found, the checker SHALL flag a conflict.
- **Resource competition**: The total ribosome demand for all genes in the circuit SHALL NOT exceed the estimated cellular ribosome capacity for the specified organism and cell type. Ribosome demand for each gene SHALL be estimated from the transcript's coding sequence length, codon adaptation index, and ribosome binding site strength. The capacity estimate SHALL be derived from literature values for the specified cell type. If the total demand exceeds capacity, the checker SHALL flag a conflict with the estimated demand distribution.
- **Splicing interference**: No cryptic splice site in one gene's transcript SHALL interfere with the splicing of another gene's transcript. Specifically, the checker SHALL verify that: (a) no splice donor from gene A pairs with a splice acceptor from gene B to create a cross-gene splicing event; (b) no regulatory element (ESE/ESS/ISE/ISS) from one gene's transcript falls within the splicing regulatory window of another gene's exon. This check is relevant only when transcripts are co-localized (e.g., polycistronic constructs or read-through transcription).
- **RNA-RNA interaction**: No complementary regions between transcripts SHALL have the potential to form stable double-stranded RNA (dsRNA) that could trigger RNA interference (RNAi) or the innate immune response. The checker SHALL scan all pairs of transcripts for complementary stretches exceeding a configurable length threshold (default: 19 bp for potential siRNA generation, 30 bp for potential dsRNA formation) and compute a simple base-pairing score. Regions exceeding the threshold SHALL be flagged as potential interactions.

**REQ-FUNC-072**: Each compositional check SHALL produce a three-valued verdict with evidence:

- **PASS**: No conflict detected. The evidence SHALL describe the search performed (e.g., "Scanned all 45 gene pairs for complementary regions ≥ 19 bp; no complementary stretches found").
- **FAIL**: A conflict is detected. The evidence SHALL identify the specific conflict: the genes involved, the positions, the nature of the conflict, and the biological consequence (e.g., "Gene LacI produces transcription factor LacI, which binds to the promoter of Gene TetR (repressible promoter with LacI as repressor), creating an unintended repression loop").
- **UNCERTAIN**: The check is inconclusive. The evidence SHALL identify the specific knowledge gap (e.g., "RNA-RNA interaction check found a 22 bp complementary region between Gene A and Gene C, but whether this forms stable dsRNA under cellular conditions depends on RNA secondary structure predictions that are not available").

**REQ-FUNC-073**: The linker SHALL compose individual gene type verdicts with composition verdicts using the three-valued logic defined in REQ-FUNC-034. The overall circuit verdict is the conjunction (AND) of: (a) each individual gene's overall verdict (itself a conjunction of all type predicates for that gene); and (b) each composition check's verdict. A circuit certificate is generated only when the overall verdict is PASS (or PASS ∧ UNCERTAIN with explicit user acceptance). If any component of the conjunction is FAIL, the circuit certificate is not generated and a failure report identifies the FAIL source(s). If any component is UNCERTAIN, the circuit certificate includes the UNCERTAIN verdicts with knowledge gaps.

### 3.9 Mutation Explorer

**REQ-FUNC-080**: The system SHALL decompose the mutation space of a gene into categories based on splicing grammar nonterminals: (a) intra-exonic mutations (point mutations within an exon that do not affect splice sites); (b) splice site mutations (mutations at or near donor/acceptor sites that may affect splicing); and (c) regulatory element mutations (mutations in ESE/ESS/ISE/ISS regions that may affect exon inclusion or exclusion). The decomposition SHALL reflect the structure of the splicing grammar: mutations in different nonterminal categories have different effects on the set of possible isoforms, and mutations in the same category have analogous effects.

**REQ-FUNC-081**: The system SHALL enumerate legal multi-mutation combinations — combinations where the splicing grammar still produces at least one valid isoform after all mutations are applied simultaneously. A "legal" mutation combination is one where: (a) every individual mutation is compatible with the grammar (i.e., the mutated sequence at each position still satisfies the grammar's terminal symbols); and (b) the combined effect of all mutations does not create a grammar violation (e.g., eliminating a required splice donor site). The enumeration SHALL be grammar-guided: instead of enumerating all possible nucleotide changes and testing each combination, the system SHALL enumerate only those combinations that are consistent with the grammar's nonterminal structure, dramatically reducing the search space.

**REQ-FUNC-082**: The system SHALL exploit independence: mutations in different exons that do not affect shared splice sites SHALL be enumerated separately, with the total count being the product of per-exon counts. Two mutations are independent if and only if they occur in different exons AND neither mutation affects a splice site or regulatory element shared between those exons. For independent mutations, the set of legal combinations factorizes: the legal combinations of mutations in exon A × the legal combinations of mutations in exon B = the legal combinations of mutations in both exons. This factorization avoids the combinatorial explosion of enumerating all cross-exon combinations.

**REQ-FUNC-083**: The system SHALL report constraint conflicts: mutation combinations that are individually legal (each mutation, considered alone, preserves the grammar) but jointly illegal (the combination of mutations creates a grammar violation). For example: mutation M1 in exon 2 eliminates a cryptic splice site but changes a codon to a rare synonymous codon; mutation M2 in exon 3 also eliminates a cryptic splice site but changes a different codon to a rare synonymous codon; individually, each mutation preserves the grammar and maintains CAI above the threshold, but jointly, the two rare codons reduce the overall CAI below the threshold, violating the `CodonAdapted` constraint. The system SHALL detect such epistatic interactions and report: (a) the set of jointly illegal mutations; (b) the constraint that is jointly violated; (c) which individual constraints remain satisfied.

### 3.10 Overlapping Reading Frame Analysis

**REQ-FUNC-090**: The system SHALL accept a nucleotide sequence with multiple annotated reading frames and construct a separate translation FST per frame. Each reading frame specification SHALL include: (a) the frame offset (0, 1, or 2, indicating which nucleotide is the first in the codon); (b) the start position and end position of the coding region within that frame; (c) a name or identifier for the protein product of that frame. The system SHALL validate that the reading frame specifications are internally consistent (start codon in the specified frame, no frame offsets that place the start codon out of frame). For each frame, the system SHALL produce a full IR-Peptide record including the amino acid sequence, codon assignments, and any special-case annotations (selenocysteine, frameshifts).

**REQ-FUNC-091**: The system SHALL compute the shared constraint set: the set of nucleotide positions where a mutation affects the amino acid sequence in more than one reading frame. For each nucleotide position i, the system SHALL determine: (a) which reading frames include position i in a codon (a nucleotide at position i is part of a codon in frame f if i mod 3 = f); (b) for each such frame, which codon position within its triplet the nucleotide occupies (first, second, or third position of the codon); (c) whether a mutation at position i changes the amino acid in each frame (a mutation at the third wobble position of a codon may or may not change the amino acid, depending on the codon and the substitution). Positions where a single mutation changes amino acids in two or more frames are added to the shared constraint set with annotations identifying the affected frames and amino acids.

**REQ-FUNC-092**: The system SHALL classify each nucleotide position as high-coupling or low-coupling:

- **High-coupling**: A position where a mutation changes the amino acid sequence in two or more reading frames simultaneously. These positions represent the strongest constraints on sequence design: any change at a high-coupling position must satisfy the optimization targets of all affected frames simultaneously. The system SHALL report the set of affected frames and the specific amino acid changes for each frame.
- **Low-coupling**: A position where a mutation changes the amino acid sequence in at most one reading frame. These positions are relatively unconstrained and can be optimized for a single frame without affecting other frames. Low-coupling positions that do not change any amino acid (e.g., synonymous wobble positions) are the least constrained.

The classification SHALL enable downstream users to prioritize: high-coupling positions require the most careful optimization and are the most likely sources of constraint conflicts; low-coupling positions offer the most design flexibility.

**REQ-FUNC-093**: The system SHALL detect constraint conflicts between frames: positions where the optimization target for one frame conflicts with the optimization target for another frame. A constraint conflict exists at position i when: (a) position i is in the shared constraint set (high-coupling); (b) the optimal codon assignment for Frame 1 at position i requires a specific nucleotide at position i; (c) the optimal codon assignment for Frame 2 at position i requires a different nucleotide at position i; and (d) no codon assignment satisfies both frames' constraints simultaneously. The system SHALL report: (a) the conflicting position(s); (b) the frames involved; (c) the amino acid requirements in each frame; (d) the set of possible codon assignments and which constraints each assignment satisfies or violates; and (e) a minimal conflict set (the smallest set of frames whose optimization targets are mutually exclusive at this position).

---

## 4. Non-Functional Requirements

### 4.1 Performance

**REQ-NFR-001**: The lexical scanner (COMP-01, implementing REQ-FUNC-001 through REQ-FUNC-004) SHALL process a 10,000 nucleotide (10 kb) input sequence, including all motif annotations (start codons, stop codons, splice donor/acceptor motifs, branch point motifs, polypyrimidine tracts, Kozak consensus, RNA instability motifs, and up to 20 restriction enzyme recognition sites), in under 1 second on a single CPU core (baseline: Intel Core i7-12700 or Apple M1). The measurement SHALL include input parsing, all DFA executions, and token stream construction. This performance target ensures that scanning is never the bottleneck in the pipeline.

**REQ-NFR-002**: The splicing grammar NDFST (COMP-02, implementing REQ-FUNC-010 through REQ-FUNC-015) SHALL compute the complete set of splice isoforms for a single human gene with an average of 10 exons (including alternative splicing variants) in under 5 seconds on a single CPU core. This performance target accounts for the combinatorial nature of isoform enumeration: genes with many alternative exons may produce large isoform sets, and the 5-second target applies to the average case. For genes with unusually complex alternative splicing patterns (> 20 exons with multiple alternative regions), the system MAY exceed this target but SHALL provide progress indication and a configurable timeout.

**REQ-NFR-003**: The translation FST (COMP-03, implementing REQ-FUNC-020 through REQ-FUNC-023) SHALL translate a 5,000 nucleotide (5 kb) spliced mRNA sequence in under 100 milliseconds on a single CPU core. This includes codon-by-codon transduction, SECIS element detection, selenocysteine and pyrrolysine flagging, and frameshift motif annotation. Translation is the simplest pipeline stage (a single-pass FST with a fixed codon table) and SHALL be the fastest.

**REQ-NFR-004**: The type checker (COMP-05, implementing REQ-FUNC-030 through REQ-FUNC-035) SHALL type-check a single mRNA sequence against all 33 type predicates (13 core + 20 SLOT-dependent) in under 10 seconds on a single CPU core. The prokaryote fast path (13 core predicates only, no splicing or FFI-dependent checks) SHALL complete in under 2 milliseconds. This includes the re-execution of the NDFST for SpliceCorrect, the re-scanning for NoCrypticSplice and NoRestrictionSite, the CAI computation for CodonAdapted, SLOT field evaluation for structure/stability/solubility/immunogenicity predicates, and the derivation trace generation for all predicates.

**REQ-NFR-005**: The CSP optimizer (COMP-06, implementing REQ-FUNC-040 through REQ-FUNC-045) SHALL find a feasible codon assignment for a protein of up to 1,000 amino acids (corresponding to an mRNA of approximately 3,000 nucleotides) in under 60 seconds on a single CPU core. This performance target assumes a typical constraint set (splicing correctness, CAI ≥ 0.8, GC ∈ [40%, 60%], ≤ 5 restriction enzymes, no instability motifs). For unusually tight constraint sets that require extensive backtracking, the system MAY exceed this target but SHALL provide progress indication and a configurable timeout (default: 60 seconds).

**REQ-NFR-006**: The compositional verifier (COMP-08, implementing REQ-FUNC-070 through REQ-FUNC-073) SHALL verify a circuit of up to 10 genes (each gene having been individually type-checked) in under 5 minutes on a single CPU core. This includes all four composition checks (promoter conflict, resource competition, splicing interference, RNA-RNA interaction) across all gene pairs. The 10-gene circuit produces 45 gene pairs for pairwise checks and 10 individual gene verdicts to compose.

### 4.2 Reliability

**REQ-NFR-010**: The system SHALL be deterministic for all internal pipeline stages: given the same inputs and parameters, the same outputs SHALL be produced. This requirement applies to all components except FFI stages (COMP-04, REQ-FUNC-050 and REQ-FUNC-051), which are explicitly treated as non-deterministic. Determinism SHALL be verified by reproducibility testing: running the same input through the pipeline 100 times and confirming byte-identical output (excluding FFI stages). This requirement is essential for: (a) regression testing (expected outputs can be stored and compared); (b) certificate verification (the independent checker must reproduce the same verdicts); and (c) scientific reproducibility (other researchers must be able to reproduce the same results).

**REQ-NFR-011**: The type checker SHALL be sound: if a PASS verdict is produced for a type predicate, the property SHALL be satisfied for the verified sequence under the specified conditions. The system SHALL NOT produce a PASS verdict for a sequence that violates the corresponding constraint. Soundness is the most critical reliability property: a false PASS (certifying a defective design) is a safety violation. Soundness SHALL be verified by adversarial testing: constructing sequences with known violations and confirming that the type checker always returns FAIL. Any soundness violation discovered post-release SHALL be treated as a Critical defect (see DOC-05, Section 6).

**REQ-NFR-012**: The CSP optimizer SHALL be complete for feasible problems: if a feasible codon assignment exists (i.e., at least one assignment satisfies all hard constraints), the solver SHALL find one. The solver SHALL NOT report INFEASIBLE when a feasible assignment exists. Completeness is verified by constructing problems with known feasible solutions and confirming that the solver finds a solution.

**REQ-NFR-013**: The CSP optimizer SHALL be correct for infeasible problems: if INFEASIBLE is reported, no feasible assignment exists under the stated constraints. The reported MUS SHALL be verified as genuinely unsatisfiable. Correctness for infeasible problems is verified by: (a) independently verifying that the MUS is unsatisfiable (using a different solver or manual analysis); and (b) confirming that adding any single constraint from the MUS to the remaining constraints creates a satisfiable problem.

### 4.3 Usability

**REQ-NFR-020**: The system SHALL provide a command-line interface (CLI) for pipeline invocation. The CLI SHALL support the following subcommands: `design` (full gene design pipeline), `verify` (type-check an existing mRNA sequence), `explore` (mutation space exploration), `analyze-orf` (overlapping reading frame analysis), `verify-circuit` (multi-gene circuit verification), and `check-cert` (independent certificate verification). The CLI SHALL provide clear, human-readable output by default and machine-readable output (JSON) via a `--json` flag. Exit codes SHALL distinguish between success (0), infeasible (1), uncertain (2), input error (10), and internal error (11).

**REQ-NFR-021**: The system SHALL provide a Python API for programmatic access to all pipeline stages. The API SHALL expose each component as a callable function with type-annotated signatures, enabling integration into custom pipelines and Jupyter notebooks. The API SHALL be documented with docstrings, usage examples, and type stubs. The API SHALL support both synchronous (blocking) and asynchronous (non-blocking, for FFI stages) invocation modes.

**REQ-NFR-022**: Error messages SHALL include three components: (a) the requirement violated — the specific REQ-FUNC or REQ-NFR identifier that is not satisfied; (b) the position in the input where the violation occurs — the nucleotide position, codon position, or gene identifier; and (c) a suggested remediation — a specific action the user can take to resolve the error. For example: "REQ-FUNC-042 (splicing correctness constraint): Cryptic splice donor site detected at position 347 (sequence: GUAAGU, score: 8.2, threshold: 7.5). Remediation: Use the optimizer to find a synonymous codon at position 115 that disrupts the donor motif."

**REQ-NFR-023**: UNCERTAIN verdicts SHALL include two components: (a) the specific information gap that prevents resolution — what the system needs to know but does not; and (b) a description of what additional constraints or data would convert the verdict to PASS or FAIL. For example: "UNCERTAIN: SpliceCorrect(HEK293T) — Weak splice donor site at position 234 (score: 6.8, threshold: 7.5) may or may not be functional in HEK293T cells depending on SRSF1 expression levels. To resolve: (1) Provide RT-PCR data for this transcript in HEK293T cells; or (2) Add a constraint eliminating all donor sites with score > 5.0 via synonymous substitution."

### 4.4 Maintainability

**REQ-NFR-030**: The IR schemas SHALL be defined in protocol buffers (.proto files), enabling backward-compatible schema evolution through protocol buffer's built-in versioning mechanisms (field numbers, optional fields, default values). New fields SHALL be added with new field numbers; existing fields SHALL NOT be removed or renumbered. Schema changes SHALL be validated by round-trip testing: serializing an IR instance to binary, deserializing it, and confirming byte-identical content. The schema SHALL be compiled to Python bindings at build time using `protoc`.

**REQ-NFR-031**: Each pipeline stage (COMP-01 through COMP-10) SHALL be independently testable through a defined input/output contract. A stage SHALL accept its defined input type (e.g., IR-Seq for COMP-02) and produce its defined output type (e.g., list[IR-Seq] for COMP-02) without requiring the execution of any preceding or succeeding stage. This enables: (a) unit testing of each stage in isolation; (b) substitution of mock implementations for testing; (c) independent development of stages by different team members.

**REQ-NFR-032**: The splicing grammar rules SHALL be specified in a declarative configuration file (YAML format), not hardcoded in source code. The configuration file SHALL contain: (a) donor and acceptor consensus sequences and PWM thresholds; (b) branch point motif specification and search window; (c) polypyrimidine tract minimum length and purity; (d) exon length soft constraints; (e) ESE/ESS/ISE/ISS motif lists and their regulatory effect strengths; (f) cell-type-specific parameter overrides. This enables updates to grammar rules (reflecting new biological knowledge) without modifying the NDFST implementation, recompiling, or redeploying the system.

**REQ-NFR-033**: External tool adapters SHALL implement a common abstract interface (`FFIAdapter`), enabling addition of new tools without modifying the pipeline core. The abstract interface SHALL define: (a) `name()` — the adapter's identifier; (b) `slot_fields()` — which IR SLOT fields the adapter fills; (c) `invoke(ir_peptide, config)` — the method that calls the external tool and populates SLOT fields; (d) `validate_output(output)` — the method that verifies schema invariants. New adapters SHALL be registered via a configuration file or entry point mechanism, without modifying any existing source code.

### 4.5 Portability

**REQ-NFR-040**: The system SHALL run on Linux (x86_64 and ARM64 architectures). All core pipeline functionality SHALL be architecture-independent; any architecture-specific code (e.g., SIMD optimizations for scanning) SHALL be optional with pure-Python fallbacks.

**REQ-NFR-041**: The system SHALL run on macOS (x86_64 and ARM64 / Apple Silicon architectures). The same codebase SHALL support both Linux and macOS without platform-specific branches (except for optional performance optimizations with fallbacks).

**REQ-NFR-042**: External tool dependencies SHALL be isolated behind the FFI interface (COMP-04), so that the core pipeline operates without any external tool installed. When an FFI stage is requested but the external tool is not available, the system SHALL: (a) log a warning indicating the tool is unavailable; (b) skip the FFI stage gracefully; (c) leave the corresponding IR SLOT fields empty; (d) continue with the remaining pipeline stages; and (e) indicate in the output (verification report, certificate) which FFI stages were skipped. The guarantee certificate SHALL be valid without FFI output; FFI-dependent annotations SHALL be clearly labeled as optional and non-deterministic.

### 4.6 Security

**REQ-NFR-050**: Certificate integrity SHALL be ensured by including a SHA-256 hash of the verified nucleotide sequence in the certificate's `design_id` field. Any modification to the sequence after certificate generation SHALL be detectable by recomputing the hash and comparing it to the certificate's `design_id`. The independent checker (REQ-FUNC-061) SHALL verify the hash as part of its validation procedure. The system SHALL NOT sign certificates with cryptographic keys (this is not a trusted computing system); the integrity guarantee is limited to tamper detection via hash comparison.

**REQ-NFR-051**: Input validation SHALL reject non-IUPAC characters in nucleotide sequences before processing. The scanner (COMP-01) SHALL validate every character in the input sequence against the IUPAC nucleotide alphabet (A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N) and raise an `InvalidSequenceError` identifying the position and offending character. This prevents injection of malicious or malformed data into the pipeline. For protein sequences, the system SHALL validate against the IUPAC protein alphabet (all 20 standard amino acids plus U, O, and ambiguity codes).

**REQ-NFR-052**: The core pipeline SHALL NOT require network access for its operation. All internal pipeline stages (COMP-01, COMP-02, COMP-03, COMP-05, COMP-06, COMP-07, COMP-08, COMP-09, COMP-10) SHALL operate entirely on local data (input sequences, local configuration files, local grammar rules). Network access is required ONLY for FFI stages that invoke cloud-based tools (e.g., ColabFold via API), and these stages are optional (REQ-NFR-042). This network isolation ensures that the core pipeline can be run in air-gapped environments (e.g., secure research networks, regulatory review environments).

---

## 5. Constraints

### 5.1 Design Constraints

**REQ-CON-001**: The system SHALL NOT use probabilistic models for any internal pipeline stage. All internal reasoning — scanning, splicing grammar parsing, translation, type checking, constraint satisfaction, certificate generation, compositional verification, mutation exploration, and ORF analysis — SHALL be deterministic. This constraint distinguishes BioCompiler from all existing gene design tools that use probabilistic scoring functions, hidden Markov models, or machine learning models for internal computation. The only exception is the FFI (COMP-04), where external tools MAY use probabilistic models internally, but the BioCompiler system treats their output as non-deterministic annotations rather than probabilistic predictions.

**REQ-CON-002**: The system SHALL NOT attempt grammar induction — the automatic learning of splicing grammar rules from biological sequence data. All grammar rules SHALL be specified from known biological knowledge: consensus sequences from the literature, PWMs derived from annotated splice sites in GENCODE, and regulatory element motifs from curated databases (ESEfinder, RESCUE-ESE). This constraint is motivated by two factors: (a) grammar induction algorithms (e.g., the Angluin L* algorithm, inside-outside re-estimation for PCFGs) require large training datasets and produce probabilistic grammars, which would violate REQ-CON-001; (b) the system's guarantee certificates (REQ-FUNC-060) require that every grammar rule has a known biological justification that can be cited in the derivation trace, which is incompatible with learned rules whose interpretation may be opaque.

**REQ-CON-003**: The system SHALL NOT model protein folding or post-translational modifications internally. These processes are not formalizable as string-to-string transformations (they depend on 3D physical chemistry and cellular signaling pathways, respectively) and are handled as foreign function calls (REQ-FUNC-050, REQ-FUNC-051). This constraint ensures that the system's internal computation remains within the domain of formal string transformations, where deterministic guarantees are achievable. The FFI boundary (COMP-04) isolates the non-deterministic, non-formalizable computation from the deterministic core.

**REQ-CON-004**: The system SHALL NOT claim that biology implements compilation. The compiler metaphor (staged transformation, typed IRs, type checking, composable passes) is used as a design pattern — a proven architectural template from software engineering — NOT as a theoretical claim about biological information processing. The system does not assert that splicing is "parsing" in the linguistic sense, or that translation is "code generation" in the compiler sense. It asserts only that the compiler design pattern provides a useful architectural framework for building bioinformatics tools with formal verification capabilities. This constraint is included to prevent over-interpretation of the metaphor and to maintain intellectual honesty about the system's scope and limitations.

### 5.2 Environmental Constraints

**REQ-CON-010**: The system SHALL operate within the memory limits of a standard scientific workstation: ≤ 32 GB RAM for single-gene analysis (all pipeline stages for a single gene, including the NDFST's isoform set representation and the CSP solver's search state) and ≤ 64 GB RAM for circuit-level analysis (10 genes processed in parallel with compositional verification). The memory constraint ensures that the system is usable on commonly available hardware without requiring specialized high-memory servers.

**REQ-CON-011**: GPU resources are required ONLY for FFI stages (specifically, AlphaFold/ColabFold structure prediction, which requires a GPU for reasonable performance). The core pipeline (COMP-01, COMP-02, COMP-03, COMP-05, COMP-06, COMP-07, COMP-08, COMP-09, COMP-10) SHALL NOT require a GPU and SHALL operate entirely on the CPU. This constraint ensures that the core pipeline is accessible on any workstation, including those without GPUs, and that the system degrades gracefully (skip FFI stages) when a GPU is not available (REQ-NFR-042).

---

## 6. Requirement Prioritization

> **Note**: Development phase numbers (Phase 1–4) referenced in the rationale column refer to project milestones and are separate from the optimizer's internal processing steps.

The following table maps each functional requirement to a MoSCoW priority category:

| Requirement ID | Description | Priority | Rationale |
|---|---|---|---|
| REQ-FUNC-001 | Accept DNA/RNA FASTA or raw input | **Must** | Foundational: all pipeline stages depend on input parsing |
| REQ-FUNC-002 | Scan and annotate biological elements | **Must** | Core scanner functionality: required by splicing engine and type system |
| REQ-FUNC-003 | Implement scanner as DFAs | **Must** | Ensures O(n) scanning and determinism (REQ-FUNC-004) |
| REQ-FUNC-004 | Scanner determinism | **Must** | Critical reliability requirement: reproducibility and certificate verification depend on it |
| REQ-FUNC-010 | NDFST for splice isoform set | **Must** | Core splicing engine: the system's primary differentiator |
| REQ-FUNC-011 | Splicing grammar rules | **Must** | Defines the biological knowledge encoded in the grammar |
| REQ-FUNC-012 | Alternative splicing as non-deterministic branching | **Must** | Captures the key biological phenomenon the system models |
| REQ-FUNC-013 | Cellular context parameter | **Should** | Enables cell-type-specific analysis; system works with default context |
| REQ-FUNC-014 | NDFST determinism | **Must** | Critical reliability requirement |
| REQ-FUNC-015 | NDFST output specification | **Must** | Defines the interface contract for downstream consumers |
| REQ-FUNC-020 | Deterministic FST for translation | **Must** | Core translation engine |
| REQ-FUNC-021 | Special cases (Sec, Pyl, frameshifts) | **Should** | Required for selenoprotein and archaeal gene design; not needed for standard human genes |
| REQ-FUNC-022 | Translation determinism | **Must** | Critical reliability requirement |
| REQ-FUNC-023 | Translation output specification | **Must** | Defines the interface contract for downstream consumers |
| REQ-FUNC-030 | Type checker with three-valued verdicts | **Must** | Core verification capability: the system's primary value proposition |
| REQ-FUNC-031 | Supported type predicates | **Must** | Defines the verification language |
| REQ-FUNC-032 | Type checker determinism | **Must** | Critical reliability requirement |
| REQ-FUNC-033 | Subtyping relations | **Should** | Enables hierarchical verification; system works without it but with redundant checks |
| REQ-FUNC-034 | Three-valued composition rules | **Must** | Required for composable verification of multi-predicate and multi-gene designs |
| REQ-FUNC-035 | Evidence per verdict (derivation trace / violation / knowledge gap) | **Must** | Required for certificate generation and user understanding |
| REQ-FUNC-040 | CSP solver for codon optimization | **Must** | Core optimization capability |
| REQ-FUNC-041 | Decision variables and domains | **Must** | Defines the CSP formulation |
| REQ-FUNC-042 | Hard constraints | **Must** | Defines the guarantee conditions |
| REQ-FUNC-043 | Feasible: return assignment(s) | **Must** | Required output for feasible problems |
| REQ-FUNC-044 | Infeasible: report MUS | **Must** | Critical for usability: tells the user why no solution exists |
| REQ-FUNC-045 | Optimizer determinism | **Must** | Critical reliability requirement |
| REQ-FUNC-050 | Folding FFI | **Should** | Enhances verification with structural data; core pipeline works without it |
| REQ-FUNC-051 | PTM prediction FFI | **Could** | Useful for therapeutic protein design; not needed for basic gene design |
| REQ-FUNC-052 | No internal modeling of external tools | **Must** | Design constraint ensuring FFI boundary integrity |
| REQ-FUNC-053 | FFI non-determinism | **Must** | Design constraint ensuring correct treatment of external tool output |
| REQ-FUNC-060 | Guarantee certificate generation | **Must** | Core value proposition: proof-carrying gene designs |
| REQ-FUNC-061 | Independent certificate verifiability | **Must** | Critical for trust: certificates must be checkable without the pipeline |
| REQ-FUNC-062 | Circuit certificate format | **Should** | Required for multi-gene circuit verification (Phase 3) |
| REQ-FUNC-070 | Circuit-level specification acceptance | **Should** | Required for multi-gene circuit verification (Phase 3) |
| REQ-FUNC-071 | Linker pass composition checks | **Should** | Core compositional verification capability |
| REQ-FUNC-072 | Three-valued verdicts for composition checks | **Should** | Consistent with the system's verification framework |
| REQ-FUNC-073 | Composition via three-valued logic | **Should** | Required for composable circuit verification |
| REQ-FUNC-080 | Mutation space decomposition | **Could** | Advanced capability; enables rational mutagenesis design |
| REQ-FUNC-081 | Legal multi-mutation enumeration | **Could** | Advanced capability; useful for directed evolution planning |
| REQ-FUNC-082 | Independence exploitation | **Could** | Performance optimization for mutation enumeration |
| REQ-FUNC-083 | Constraint conflict detection | **Could** | Advanced capability; detects epistatic interactions |
| REQ-FUNC-090 | Multiple reading frame acceptance | **Could** | Advanced capability; primarily for virology applications |
| REQ-FUNC-091 | Shared constraint set computation | **Could** | Core ORF analysis capability |
| REQ-FUNC-092 | Coupling classification | **Could** | Enables prioritized optimization of overlapping frames |
| REQ-FUNC-093 | Inter-frame constraint conflict detection | **Could** | Detects design infeasibility that single-frame analysis misses |

---

## 7. Traceability Matrix (Summary)

The following table maps every requirement to its architectural component, interface, test case, and associated risk. Full bidirectional traceability is maintained in DOC-08 (Traceability Matrix).

### 7.1 Functional Requirements Traceability

| Requirement | Architecture Component | Interface | Test Case(s) | Risk |
|---|---|---|---|---|
| REQ-FUNC-001 | COMP-01 (Scanner) | IF-01 | TC-U-001, TC-U-002, TC-U-003 | RISK-04 |
| REQ-FUNC-002 | COMP-01 (Scanner) | IF-01 | TC-U-004, TC-U-005, TC-U-006 | — |
| REQ-FUNC-003 | COMP-01 (Scanner) | IF-01 | TC-U-003 | — |
| REQ-FUNC-004 | COMP-01 (Scanner) | IF-01 | TC-U-007, TC-REP-001 | — |
| REQ-FUNC-010 | COMP-02 (Splicing Engine) | IF-02 | TC-U-010, TC-U-011, TC-V-001, TC-V-002 | RISK-01 |
| REQ-FUNC-011 | COMP-02 (Splicing Engine) | IF-02 | TC-U-012, TC-U-013 | RISK-01 |
| REQ-FUNC-012 | COMP-02 (Splicing Engine) | IF-02 | TC-U-011 | — |
| REQ-FUNC-013 | COMP-02 (Splicing Engine) | IF-02 | TC-U-014 | RISK-02 |
| REQ-FUNC-014 | COMP-02 (Splicing Engine) | IF-02 | TC-U-015, TC-REP-001, TC-REP-002 | — |
| REQ-FUNC-015 | COMP-02 (Splicing Engine) | IF-02 | TC-I-001 | — |
| REQ-FUNC-020 | COMP-03 (Translation Engine) | IF-03 | TC-U-020, TC-V-003 | — |
| REQ-FUNC-021 | COMP-03 (Translation Engine) | IF-03 | TC-U-021, TC-U-023 | — |
| REQ-FUNC-022 | COMP-03 (Translation Engine) | IF-03 | TC-U-022, TC-REP-001 | — |
| REQ-FUNC-023 | COMP-03 (Translation Engine) | IF-03 | TC-I-002 | — |
| REQ-FUNC-030 | COMP-05 (Type System) | IF-05 | TC-U-030, TC-U-031, TC-SND-001–006 | RISK-03 |
| REQ-FUNC-031 | COMP-05 (Type System) | IF-05 | TC-U-030, TC-U-031 | — |
| REQ-FUNC-032 | COMP-05 (Type System) | IF-05 | TC-U-032, TC-REP-001 | — |
| REQ-FUNC-033 | COMP-05 (Type System) | IF-05 | TC-U-033 | — |
| REQ-FUNC-034 | COMP-05 (Type System) | IF-05 | TC-U-034, TC-U-035 | RISK-08 |
| REQ-FUNC-035 | COMP-05 (Type System) | IF-05 | TC-U-030, TC-U-031, TC-SND-001–006 | — |
| REQ-FUNC-040 | COMP-06 (Optimizer) | IF-06 | TC-U-040, TC-U-041, TC-SND-007 | RISK-02 |
| REQ-FUNC-041 | COMP-06 (Optimizer) | IF-06 | TC-U-040 | — |
| REQ-FUNC-042 | COMP-06 (Optimizer) | IF-06 | TC-U-042, TC-V-004 | — |
| REQ-FUNC-043 | COMP-06 (Optimizer) | IF-06 | TC-U-043 | — |
| REQ-FUNC-044 | COMP-06 (Optimizer) | IF-06 | TC-U-044 | — |
| REQ-FUNC-045 | COMP-06 (Optimizer) | IF-06 | TC-U-045, TC-REP-001, TC-REP-003 | — |
| REQ-FUNC-050 | COMP-04 (FFI Manager) | IF-04 | TC-I-003 | RISK-05 |
| REQ-FUNC-051 | COMP-04 (FFI Manager) | IF-04 | TC-I-003 | RISK-05 |
| REQ-FUNC-052 | COMP-04 (FFI Manager) | IF-04 | TC-I-003 | — |
| REQ-FUNC-053 | COMP-04 (FFI Manager) | IF-04 | TC-I-003 | — |
| REQ-FUNC-060 | COMP-07 (Certificate Gen) | IF-07 | TC-I-006, TC-I-007, TC-SND-008 | — |
| REQ-FUNC-061 | COMP-07 (Certificate Gen) | IF-07 | TC-I-007 | — |
| REQ-FUNC-062 | COMP-07 (Certificate Gen) | IF-07 | TC-I-007 | — |
| REQ-FUNC-070 | COMP-08 (Comp. Verifier) | IF-08 | TC-I-008, TC-S-003 | RISK-02 |
| REQ-FUNC-071 | COMP-08 (Comp. Verifier) | IF-08 | TC-I-008 | — |
| REQ-FUNC-072 | COMP-08 (Comp. Verifier) | IF-08 | TC-I-008 | — |
| REQ-FUNC-073 | COMP-08 (Comp. Verifier) | IF-08 | TC-I-008 | RISK-08 |
| REQ-FUNC-080 | COMP-09 (Mutation Explorer) | IF-09 | TC-S-002 | — |
| REQ-FUNC-081 | COMP-09 (Mutation Explorer) | IF-09 | TC-S-002 | — |
| REQ-FUNC-082 | COMP-09 (Mutation Explorer) | IF-09 | TC-S-002 | — |
| REQ-FUNC-083 | COMP-09 (Mutation Explorer) | IF-09 | TC-S-002 | — |
| REQ-FUNC-090 | COMP-10 (ORF Analyzer) | IF-10 | TC-I-009, TC-S-004, TC-V-005 | — |
| REQ-FUNC-091 | COMP-10 (ORF Analyzer) | IF-10 | TC-I-009, TC-V-005 | — |
| REQ-FUNC-092 | COMP-10 (ORF Analyzer) | IF-10 | TC-S-004 | — |
| REQ-FUNC-093 | COMP-10 (ORF Analyzer) | IF-10 | TC-S-004 | — |

### 7.2 Non-Functional Requirements Traceability

| Requirement | Architecture Component | Interface | Test Case(s) | Risk |
|---|---|---|---|---|
| REQ-NFR-001 | COMP-01 | IF-01 | Performance benchmark | — |
| REQ-NFR-002 | COMP-02 | IF-02 | Performance benchmark | RISK-01 |
| REQ-NFR-003 | COMP-03 | IF-03 | Performance benchmark | — |
| REQ-NFR-004 | COMP-05 | IF-05 | Performance benchmark | — |
| REQ-NFR-005 | COMP-06 | IF-06 | Performance benchmark | RISK-03 |
| REQ-NFR-006 | COMP-08 | IF-08 | Performance benchmark | — |
| REQ-NFR-010 | All internal components | All internal IFs | TC-REP-001–004 | RISK-07 |
| REQ-NFR-011 | COMP-05 | IF-05 | TC-SND-001–008 | RISK-07 |
| REQ-NFR-012 | COMP-06 | IF-06 | TC-U-043 | — |
| REQ-NFR-013 | COMP-06 | IF-06 | TC-U-044 | — |
| REQ-NFR-020 | CLI | — | Usability test | — |
| REQ-NFR-021 | Python API | — | Usability test | — |
| REQ-NFR-022 | All components | — | Usability test | — |
| REQ-NFR-023 | COMP-05 | IF-05 | TC-U-031 | — |
| REQ-NFR-030 | IR Bus | — | Round-trip test | RISK-04 |
| REQ-NFR-031 | All components | All IFs | Unit test isolation | — |
| REQ-NFR-032 | COMP-02 | IF-02 | Config change test | — |
| REQ-NFR-033 | COMP-04 | IF-04 | Adapter registration test | RISK-05 |
| REQ-NFR-040 | All components | — | Platform test (Linux x86_64, ARM64) | — |
| REQ-NFR-041 | All components | — | Platform test (macOS x86_64, ARM64) | — |
| REQ-NFR-042 | COMP-04 | IF-04 | Graceful degradation test | — |
| REQ-NFR-050 | COMP-07 | IF-07 | Hash verification test | — |
| REQ-NFR-051 | COMP-01 | IF-01 | TC-U-003 | — |
| REQ-NFR-052 | All internal components | — | Network isolation test | — |

### 7.3 Constraint Traceability

| Requirement | Architecture Component | Interface | Test Case(s) | Risk |
|---|---|---|---|---|
| REQ-CON-001 | All internal components | — | Code review; no random modules in internal code | — |
| REQ-CON-002 | COMP-02 | IF-02 | Config-driven grammar rules test | — |
| REQ-CON-003 | COMP-04 | IF-04 | FFI-only folding/PTM test | — |
| REQ-CON-004 | — | — | Documentation review; no theoretical claims | — |
| REQ-CON-010 | All components | — | Memory profiling test | RISK-03 |
| REQ-CON-011 | COMP-04 | IF-04 | CPU-only pipeline test | — |

---

*End of DOC-01: Software Requirements Specification (SRS), Version 12.0.0*

*Full bidirectional traceability is maintained in DOC-08 (Traceability Matrix). Requirement IDs are stable across all documents in this set.*
