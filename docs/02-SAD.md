# DOC-02: Software Architecture Document (SAD)

> **WARNING:** **NOTE: The module structure described in §4.1 represents the original design. The actual implementation uses a DDD-style Python package structure. See `src/biocompiler/` for the current layout. Key differences: Protocol Buffer IR schemas were NOT implemented (Python dataclasses in `shared/types.py` are used instead); the FFI Manager was replaced by the `engines/` package; `certificate/` → `provenance/`; `scanner/` → `sequence/scanner.py`; `splicing/` → `sequence/splicing.py`; `translation/` → `expression/translation.py`; `config/` with YAML → `grammars/` and `organisms/config.py`. The `three_valued.py` module was replaced by `shared/five_valued_logic.py` (five-valued logic).**


| Field    | Value              |
|----------|--------------------|
| ID       | DOC-02             |
| Version  | 0.9.0        |
| Status   | Current        |
| Date     | 2026-06-07         |
| Project  | BioCompiler        |
| Standard | ISO/IEC/IEEE 42010 |

---

## 1. Architecture Overview

### 1.1 Purpose

This Software Architecture Document (SAD) specifies the architectural design of the
BioCompiler system — a compiler-engineering framework that applies deterministic,
formally grounded design patterns to the formalizable stages of gene-to-protein
processing. The document serves as the authoritative reference for:

- **Stakeholders** who need to understand how the system is structured and why.
- **Developers** who implement, extend, or maintain individual components.
- **Verifiers** who must confirm that the architecture satisfies its quality
  requirements, especially determinism and soundness.
- **Integrators** who compose BioCompiler outputs into larger synthetic-biology
  toolchains.

The architecture is expressed using the **4+1 view model** (Kruchten, 1995) as
recommended by ISO/IEC/IEEE 42010:2022. Each view addresses a different concern
set; together they provide a complete picture of the system's structure and
behavior.

### 1.2 Scope — 4+1 View Model

| View              | Concern                                              | Section |
|-------------------|------------------------------------------------------|---------|
| Logical View      | Component decomposition, interfaces, IR bus          | §2      |
| Process View      | Runtime behavior, concurrency, data flow             | §3      |
| Development View  | Module structure, build, coding standards            | §4      |
| Physical View     | Deployment, hardware, storage, network               | §5      |
| Scenario View     | Use-case realizations that validate the architecture | §6      |

The architecture is **pipeline-oriented**, mirroring LLVM's multi-pass design.
Each stage consumes and produces well-typed intermediate representations (IRs)
that flow through an IR Bus. Stages that require non-deterministic models
(splicing) or external oracles (folding, PTMs) are explicitly isolated behind
well-defined boundaries, ensuring that the core pipeline remains deterministic
and testable.

### 1.3 Architectural Drivers

| # | Driver                   | Priority  | Description                                                                                        |
|---|--------------------------|-----------|----------------------------------------------------------------------------------------------------|
| D1| Determinism              | Critical  | Given identical inputs, the pipeline must produce bit-identical outputs across runs and platforms. |
| D2| Soundness                | Critical  | No PASS verdict may be issued for a gene design that violates a stated constraint.                 |
| D3| Separation of Paradigms  | High      | Deterministic compiler stages must not depend on probabilistic models; FFI isolates the boundary.  |
| D4| Composability            | High      | Multi-gene circuits must be verifiable from per-gene certificates without re-running the pipeline. |
| D5| Verifiability            | High      | Every PASS/FAIL verdict must be accompanied by a certificate that is independently checkable.      |
| D6| Extensibility            | Medium    | New organisms, folding algorithms, and constraint classes must be addable without modifying core.   |

### 1.4 Architectural Decisions

| ID    | Decision                                       | Rationale                                                                                                  | Alternatives Considered                             | Status   |
|-------|------------------------------------------------|------------------------------------------------------------------------------------------------------------|-----------------------------------------------------|----------|
| AD-01 | Pipeline architecture (not monolith)           | Mirrors LLVM's proven multi-pass design; each stage is independently testable and replaceable.             | Monolithic class hierarchy, actor model             | Accepted |
| AD-02 | Protocol Buffers for IR schemas                | Language-neutral, schema-evolution–friendly, efficient binary serialization; supports cross-toolchain IR.  | JSON, ASN.1, Cap'n Proto, flatbuffers              | Accepted |
| AD-03 | Non-deterministic FSTs (NDFST) for splicing    | Splicing is inherently non-deterministic (multiple isoforms); NDFSTs model this naturally via set-valued output. | Deterministic approximations, Monte Carlo sampling | Accepted |
| AD-04 | Constraint satisfaction (not optimization)     | Gene design is a feasibility problem; CSP + MUS diagnosis provides provably complete explanations of failure. | Genetic algorithms, Bayesian optimization           | Accepted |
| AD-05 | Five-valued logic for verdicts                | PASS/FAIL/UNCERTAIN captures the epistemic boundary where FFI oracles cannot provide definitive answers.   | Boolean (pass/fail), probabilistic scores           | Accepted |
| AD-06 | FFI for folding and PTMs                       | Folding and post-translational modifications require physics-based or ML oracles; FFI keeps the core deterministic. | Reimplement in-core, REST microservices             | Accepted |
| AD-07 | Declarative grammar configuration              | Organism-specific genetic grammars (splice-site rules, codon tables) are loaded from declarative configs, not hard-coded. | Hard-coded tables, database-backed configs          | Accepted |

---

## 2. Logical View

### 2.1 Component Decomposition

```
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │                              BioCompiler Core                                │
 │                                                                              │
 │  ┌──────────┐   ┌──────────────┐   ┌───────────────┐   ┌────────────────┐   │
 │  │ COMP-01  │   │   COMP-02    │   │   COMP-03     │   │   COMP-04      │   │
 │  │ Scanner  │──▶│  Splicing    │──▶│  Translation  │──▶│  FFI Manager   │   │
 │  │ (DFA)    │   │  Engine      │   │  Engine       │   │  (Folding/PTM) │   │
 │  └──────────┘   │  (NDFST)     │   │  (Det. FST)   │   │                │   │
 │       │         └──────────────┘   └───────────────┘   └───────┬────────┘   │
 │       │                │                   │                     │            │
 │       ▼                ▼                   ▼                     ▼            │
 │  ┌─────────────────────────────────────────────────────────────────────┐     │
 │  │                         IR  Bus                                     │     │
 │  │   IR-Seq        IR-Peptide       IR-Structure       IR-Circuit     │     │
 │  └──────┬──────────────┬─────────────────┬────────────────┬───────────┘     │
 │         │              │                 │                │                  │
 │         ▼              ▼                 ▼                ▼                  │
 │  ┌───────────┐  ┌────────────┐  ┌──────────────┐  ┌───────────────────┐    │
 │  │ COMP-10   │  │  COMP-05   │  │   COMP-06    │  │    COMP-08        │    │
 │  │ ORF       │  │  Type      │  │   Optimizer  │  │    Compositional  │    │
 │  │ Analyzer  │  │  System    │  │   (CSP+MUS)  │  │    Verifier       │    │
 │  └───────────┘  └────────────┘  └──────────────┘  └───────────────────┘    │
 │                       │                 │                     │              │
 │                       ▼                 ▼                     ▼              │
 │                ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐     │
 │                │  COMP-07     │  │  COMP-09     │  │                  │     │
 │                │  Certificate │  │  Mutation    │  │                  │     │
 │                │  Generator   │  │  Explorer    │  │                  │     │
 │                └──────────────┘  └──────────────┘  └──────────────────┘     │
 └──────────────────────────────────────────────────────────────────────────────┘
```

**Data flow summary:** Raw DNA sequence enters COMP-01 (Scanner) and is
tokenized into IR-Seq. COMP-02 (Splicing Engine) transforms IR-Seq into one
or more IR-Seq variants (one per isoform). COMP-03 (Translation Engine)
converts each IR-Seq variant into IR-Peptide. COMP-04 (FFI Manager) invokes
external oracles to annotate IR-Peptide with structural data, producing
IR-Structure. COMP-05 (Type System) assigns PASS/FAIL/UNCERTAIN verdicts
at each level. COMP-06 (Optimizer) solves constraint-satisfaction problems
over IR-Structure. COMP-07 (Certificate Generator) emits independently
verifiable JSON certificates. COMP-08 (Compositional Verifier) links
per-gene certificates into circuit-level proofs. COMP-09 (Mutation Explorer)
enumerates grammar-guided variants. COMP-10 (ORF Analyzer) identifies and
characterizes overlapping reading frames.

### 2.2 Component Catalog

| ID      | Name                     | Responsibility                                                                 | Key Requirements                          | Dependencies               |
|---------|--------------------------|--------------------------------------------------------------------------------|-------------------------------------------|----------------------------|
| COMP-01 | Scanner                  | DFA-based lexical analysis of raw DNA sequences; tokenization into codons, regulatory motifs, junction signals. | REQ-01 (Determinism), REQ-02 (Soundness) | None (entry point)         |
| COMP-02 | Splicing Engine          | NDFST-driven simulation of pre-mRNA splicing; produces set-valued isoform output (all possible mRNA transcripts). | REQ-03 (Completeness), REQ-05 (Separation) | COMP-01, IR-Seq            |
| COMP-03 | Translation Engine       | Deterministic FST mapping from mRNA codons to amino acids; handles alternative codon tables per organism. | REQ-01 (Determinism), REQ-02 (Soundness) | COMP-02, IR-Seq            |
| COMP-04 | FFI Manager              | Adapter layer for external folding solvers and PTM prediction tools; serializes calls, handles timeouts and fallbacks. | REQ-05 (Separation), REQ-09 (Extensibility) | COMP-03, IR-Peptide        |
| COMP-05 | Type System              | Three-valued (PASS/FAIL/UNCERTAIN) verdict engine; checks structural and functional constraints at each IR level. | REQ-02 (Soundness), REQ-04 (Verifiability) | All IR levels              |
| COMP-06 | Optimizer                | Constraint satisfaction solver for gene design; MUS (Minimal Unsatisfiable Subset) diagnosis on failure. | REQ-06 (Completeness), REQ-02 (Soundness) | COMP-05, IR-Structure      |
| COMP-07 | Certificate Generator    | Emits JSON certificates for every verdict; certificates are independently verifiable without BioCompiler. | REQ-04 (Verifiability), REQ-07 (Composability) | COMP-05, COMP-06           |
| COMP-08 | Compositional Verifier   | Linker for multi-gene circuits; composes per-gene certificates into circuit-level proof without re-running pipeline. | REQ-07 (Composability), REQ-02 (Soundness) | COMP-07, IR-Circuit        |
| COMP-09 | Mutation Explorer        | Grammar-guided enumeration of sequence variants; generates candidates within edit-distance bounds. | REQ-08 (Extensibility), REQ-10 (Exploration) | COMP-01, COMP-05           |
| COMP-10 | ORF Analyzer             | Identifies overlapping reading frames in IR-Seq; computes frame-shifted translations and constraint interactions. | REQ-11 (ORF Analysis), REQ-02 (Soundness) | COMP-01, IR-Seq            |

### 2.3 IR Bus

The IR Bus is the backbone of the pipeline. Every component communicates
exclusively through well-typed IR records; no component directly accesses
another component's internal state. All IR schemas are defined in Protocol
Buffers (`.proto`) and versioned independently.

| IR Level    | Schema File             | Produced By          | Consumed By                                          | Description                                                                 |
|-------------|-------------------------|----------------------|------------------------------------------------------|-----------------------------------------------------------------------------|
| IR-Seq      | `ir_seq_v1.proto`       | COMP-01, COMP-02     | COMP-02, COMP-03, COMP-05, COMP-09, COMP-10          | Nucleotide sequence with annotated regions: exons, introns, splice sites, promoters, start/stop codons. Each isoform from COMP-02 is a distinct IR-Seq record with a unique isoform ID. |
| IR-Peptide  | `ir_peptide_v1.proto`   | COMP-03              | COMP-04, COMP-05, COMP-09                            | Amino acid chain with codon provenance, signal peptide annotations, and domain boundaries. Includes back-reference to source IR-Seq isoform. |
| IR-Structure| `ir_structure_v1.proto` | COMP-04              | COMP-05, COMP-06, COMP-07                            | Peptide with secondary/tertiary structure annotations from FFI oracles. Includes confidence scores (0.0–1.0) for each structural element, PTM site predictions, and folding energy estimates. |
| IR-Circuit  | `ir_circuit_v1.proto`   | COMP-08              | COMP-05, COMP-07, COMP-08                            | Directed graph of gene nodes with interaction edges. Each node references a per-gene certificate. Edges carry interaction constraints (e.g., promoter interference, metabolic burden). |

**IR versioning rules:**

- Schema versions follow semver (MAJOR.MINOR.PATCH).
- Minor and patch changes are backward-compatible; consumers ignore unknown fields.
- Major changes require a new schema file (e.g., `ir_seq_v2.proto`) and an explicit migration path.
- Every IR record carries a `schema_version` field; the pipeline rejects records with unknown major versions.

**IR Bus interface contract:**

```python
class IRBus(Protocol):
    def publish(self, ir_level: str, record: bytes, metadata: dict) -> None: ...
    def subscribe(self, ir_level: str, consumer_id: str) -> Iterator[IRRecord]: ...
    def query(self, ir_level: str, filter: IRFilter) -> Sequence[IRRecord]: ...
```

The bus guarantees:

1. **Ordering:** Records for a given pipeline run are delivered in publication order.
2. **Immutability:** Once published, an IR record cannot be modified; corrections are new records with a `supersedes` field.
3. **Completeness:** A consumer can request a replay of all records for a given `ir_level` and `run_id`.

### 2.4 IR Invariants

Invariants are predicates that must hold for every IR record published to the
bus. Violation of any invariant is a fatal error that halts the pipeline and
generates a diagnostic report.

| ID          | IR Level    | Invariant                                                                                   | Enforced By |
|-------------|-------------|---------------------------------------------------------------------------------------------|-------------|
| INV-SEQ-01  | IR-Seq      | Length in nucleotides is a positive multiple of 3 for coding regions.                       | COMP-01     |
| INV-SEQ-02  | IR-Seq      | Every splice-site annotation has a matching donor-acceptor pair within the same contig.     | COMP-02     |
| INV-SEQ-03  | IR-Seq      | Isoform IDs are unique within a pipeline run.                                               | COMP-02     |
| INV-PEP-01  | IR-Peptide  | Peptide length equals (coding-region length / 3) − 1 (stop codon excluded).                 | COMP-03     |
| INV-PEP-02  | IR-Peptide  | Every amino acid residue has a back-reference to its source codon in IR-Seq.                 | COMP-03     |
| INV-PEP-03  | IR-Peptide  | No duplicate isoform references within the same peptide record.                              | COMP-03     |
| INV-STR-01  | IR-Structure| Every structural element has a confidence score in [0.0, 1.0].                               | COMP-04     |
| INV-STR-02  | IR-Structure| PTM site predictions reference valid residue positions in the parent IR-Peptide.             | COMP-04     |
| INV-STR-03  | IR-Structure| Folding energy estimates are finite IEEE 754 doubles.                                        | COMP-04     |
| INV-CIR-01  | IR-Circuit  | The circuit graph is acyclic (no circular promoter dependencies).                           | COMP-08     |
| INV-CIR-02  | IR-Circuit  | Every gene node references a valid, non-expired certificate from COMP-07.                   | COMP-08     |

---

## 3. Process View

### 3.1 Single-Gene Pipeline

```
 ┌───────────┐     ┌──────────────┐     ┌───────────────┐     ┌───────────────┐
 │  DNA Seq  │────▶│   COMP-01    │────▶│   COMP-02     │────▶│   COMP-03     │
 │  (input)  │     │   Scanner    │     │  Splicing     │     │  Translation  │
 └───────────┘     └──────────────┘     │  Engine       │     │  Engine       │
                         │              └───────┬───────┘     └───────┬───────┘
                         │                      │                     │
                         ▼                      ▼                     ▼
                    ┌──────────┐          ┌──────────┐          ┌──────────┐
                    │  IR-Seq  │          │  IR-Seq  │          │IR-Peptide│
                    │ (tokens) │          │(isoforms)│          │          │
                    └──────────┘          └──────────┘          └──────────┘
                                                                       │
                                              ┌────────────────────────┘
                                              ▼
                                       ┌───────────────┐     ┌───────────────┐
                                       │   COMP-04     │────▶│  COMP-05      │
                                       │  FFI Manager  │     │  Type System  │
                                       └───────┬───────┘     └───────┬───────┘
                                               │                     │
                                               ▼                     ▼
                                        ┌────────────┐        ┌───────────┐
                                        │IR-Structure│        │  Verdict  │
                                        │            │        │  (3-val)  │
                                        └────────────┘        └───────────┘
                                              │                     │
                                              └────────┬────────────┘
                                                       ▼
                                              ┌───────────────┐
                                              │   COMP-06     │──────┐
                                              │  Optimizer    │      │
                                              │  (CSP + MUS)  │      │
                                              └───────┬───────┘      │
                                                      │              │
                                                      ▼              ▼
                                              ┌───────────────┐ ┌──────────────┐
                                              │  COMP-07      │ │  Optimized   │
                                              │  Certificate  │ │  IR-Structure│
                                              │  Generator    │ │  (if feasible│
                                              └───────┬───────┘ │  solution)   │
                                                      │         └──────────────┘
                                                      ▼
                                               ┌─────────────┐
                                               │ Certificate │
                                               │  (JSON)     │
                                               └─────────────┘
```

**Stage execution rules:**

1. COMP-01 → COMP-02 is a strict sequential dependency; splicing requires
   complete tokenization.
2. COMP-02 may produce **multiple** IR-Seq records (one per isoform). Each
   isoform is processed independently through COMP-03 and COMP-04.
3. COMP-04 calls to FFI oracles are **serialized per isoform** but may run
   **concurrently across isoforms** (see §3.3).
4. COMP-05 runs **after** each stage, checking invariants and assigning
   verdicts. A FAIL verdict at any stage halts downstream processing for that
   isoform.
5. COMP-06 receives the IR-Structure and the type-system verdict. If the
   verdict is PASS, COMP-06 attempts constraint optimization. If FAIL, COMP-06
   performs MUS diagnosis to identify the minimal set of unsatisfiable
   constraints.
6. COMP-07 generates the certificate regardless of verdict — the certificate
   encodes the verdict, evidence, and (on failure) the MUS.

### 3.2 Multi-Gene Circuit Pipeline

```
 ┌───────────┐  ┌───────────┐       ┌───────────┐
 │  Gene A   │  │  Gene B   │  ...  │  Gene N   │
 └─────┬─────┘  └─────┬─────┘       └─────┬─────┘
       │              │                    │
       ▼              ▼                    ▼
 ┌─────────────────────────────────────────────────┐
 │          Single-Gene Pipeline (§3.1)             │
 │          (parallel across genes)                 │
 └────┬───────────────┬───────────────────┬────────┘
      │               │                   │
      ▼               ▼                   ▼
 ┌──────────┐  ┌──────────┐        ┌──────────┐
 │ Cert A   │  │ Cert B   │  ...   │ Cert N   │
 └────┬─────┘  └────┬─────┘        └────┬─────┘
      │              │                    │
      └──────────────┴────────────────────┘
                     │
                     ▼
            ┌──────────────────┐
            │    COMP-08       │
            │   Compositional  │
            │   Verifier       │
            └────────┬─────────┘
                     │
              ┌──────┴──────┐
              ▼             ▼
      ┌──────────────┐ ┌──────────────┐
      │  IR-Circuit  │ │  Circuit     │
      │              │ │  Certificate │
      └──────────────┘ └──────────────┘
```

**Circuit verification rules:**

1. Each gene's single-gene pipeline runs **independently and in parallel**.
2. COMP-08 collects all per-gene certificates and the circuit topology
   (interaction edges).
3. COMP-08 checks **cross-gene constraints**: promoter interference, metabolic
   burden, resource competition, and regulatory crosstalk.
4. If all per-gene certificates are PASS and all cross-gene constraints are
   satisfied, COMP-08 issues a circuit-level PASS certificate.
5. If any per-gene certificate is FAIL, COMP-08 propagates the failure and
   identifies which cross-gene interactions are affected.
6. The circuit certificate includes the full dependency graph, enabling
   downstream consumers to re-verify without re-running the pipeline.

### 3.3 Concurrency Model

The BioCompiler pipeline employs a **structured concurrency** model designed to
preserve determinism while allowing limited parallelism where it is safe.

| Concern                    | Strategy                                                                                 |
|----------------------------|------------------------------------------------------------------------------------------|
| Pure-function stages       | COMP-01, COMP-02, COMP-03, and COMP-05 are pure functions. They are trivially parallelizable across independent inputs (e.g., different genes or isoforms). No synchronization is needed beyond input ordering. |
| FFI serialization          | COMP-04 (FFI Manager) serializes calls to each external oracle. Calls to **different** oracles may proceed concurrently. Calls to the **same** oracle are queued to prevent race conditions in the oracle's internal state. |
| Synchronization barriers   | Between IR levels, a synchronization barrier ensures all producers for a given level have completed before any consumer begins. For example, all isoforms of a gene must finish translation before COMP-04 starts folding any of them. |
| Shared state               | No component holds mutable shared state. All communication is through the IR Bus, which is append-only. Reads are snapshot-consistent. |
| Thread pool                | A configurable `ThreadPoolExecutor` (default: `min(cpu_count, 8)`) manages concurrent stage invocations. FFI calls use a separate `ThreadPoolExecutor` (default: `min(gpu_count + 1, 4)`) to avoid blocking the main pool. |

**Concurrency invariant:** For a single gene with a single isoform, the
pipeline is strictly sequential and bit-deterministic. Non-determinism can
only arise from (a) FFI oracle non-determinism, which is captured as
UNCERTAIN verdicts, or (b) scheduling order of independent genes, which does
not affect final results because gene pipelines are independent.

### 3.4 Error Propagation

Errors in BioCompiler are classified into three categories, each with a
distinct propagation strategy:

| Category        | Example                                   | Propagation                                                                                   |
|-----------------|-------------------------------------------|-----------------------------------------------------------------------------------------------|
| **Fatal**       | Schema violation, invariant breach, OOM   | Immediate halt of the entire pipeline. A diagnostic report is written to stderr and the IR Bus. No certificate is generated. |
| **Constraint**  | Type-system FAIL, unsatisfiable CSP       | The failing isoform/gene is marked with a FAIL verdict. Downstream stages receive the FAIL verdict and skip processing. Other isoforms/genes continue. MUS diagnosis is triggered automatically. |
| **Uncertain**   | FFI timeout, low-confidence folding score | The affected annotation is marked with confidence < threshold. COMP-05 issues an UNCERTAIN verdict. Downstream stages proceed but propagate the uncertainty flag. The certificate includes a disclaimer. |

**Error context chain:** Every error carries a chain of context:

```
ErrorContext {
    stage:        str     # e.g., "COMP-04/FFI-Manager"
    ir_level:     str     # e.g., "IR-Structure"
    record_id:    UUID    # The IR record that triggered the error
    isoform_id:   UUID    # The isoform, if applicable
    gene_id:      UUID    # The gene, if applicable
    run_id:       UUID    # The pipeline run
    message:      str     # Human-readable description
    cause:        ErrorContext | None  # Chained cause
}
```

This context chain is included in the certificate (for Constraint/Uncertain
errors) or the diagnostic report (for Fatal errors), enabling precise root-
cause analysis.

---

## 4. Development View

### 4.1 Module Structure

```
biocompiler/
├── __init__.py                  # Package root; version string
├── scanner/                     # COMP-01: Scanner
│   ├── __init__.py
│   ├── dfa.py                   # DFA engine (table-driven)
│   ├── tokenizer.py             # Sequence tokenizer
│   ├── motifs.py                # Regulatory motif definitions
│   └── tests/
│       ├── test_dfa.py
│       ├── test_tokenizer.py
│       └── test_motifs.py
├── splicing/                    # COMP-02: Splicing Engine
│   ├── __init__.py
│   ├── ndfst.py                 # Non-deterministic FST core
│   ├── isoform_builder.py       # Set-valued isoform construction
│   ├── grammar.py               # Splice-site grammar rules
│   └── tests/
│       ├── test_ndfst.py
│       ├── test_isoform_builder.py
│       └── test_grammar.py
├── translation/                 # COMP-03: Translation Engine
│   ├── __init__.py
│   ├── fst.py                   # Deterministic FST core
│   ├── codon_table.py           # Organism-specific codon tables
│   ├── signal_peptide.py        # Signal peptide detection
│   └── tests/
│       ├── test_fst.py
│       ├── test_codon_table.py
│       └── test_signal_peptide.py
├── ffi/                         # COMP-04: FFI Manager
│   ├── __init__.py
│   ├── manager.py               # FFI orchestration and serialization
│   ├── folding_adapter.py       # Adapter for folding solvers (AlphaFold, Rosetta)
│   ├── ptm_adapter.py           # Adapter for PTM prediction tools
│   ├── fallback.py              # Fallback strategies on timeout/error
│   └── tests/
│       ├── test_manager.py
│       ├── test_folding_adapter.py
│       ├── test_ptm_adapter.py
│       └── test_fallback.py
├── type_system/                 # COMP-05: Type System
│   ├── __init__.py
│   ├── shared/five_valued_logic.py          # PASS / FAIL / UNCERTAIN logic
│   ├── checker.py               # Constraint checker engine
│   ├── constraints/             # Constraint definitions by IR level
│   │   ├── seq_constraints.py
│   │   ├── peptide_constraints.py
│   │   ├── structure_constraints.py
│   │   └── circuit_constraints.py
│   └── tests/
│       ├── test_shared/five_valued_logic.py
│       ├── test_checker.py
│       └── test_constraints.py
├── optimizer/                   # COMP-06: Optimizer
│   ├── __init__.py
│   ├── csp_solver.py            # Constraint satisfaction solver
│   ├── mus_diagnosis.py         # Minimal Unsatisfiable Subset finder
│   ├── variable_domain.py       # Domain definitions for CSP variables
│   └── tests/
│       ├── test_csp_solver.py
│       ├── test_mus_diagnosis.py
│       └── test_variable_domain.py
├── certificate/                 # COMP-07: Certificate Generator
│   ├── __init__.py
│   ├── generator.py             # Certificate construction
│   ├── schema.py                # JSON schema definition
│   ├── verifier.py              # Standalone certificate verifier
│   └── tests/
│       ├── test_generator.py
│       ├── test_schema.py
│       └── test_verifier.py
├── compositional/               # COMP-08: Compositional Verifier
│   ├── __init__.py
│   ├── linker.py                # Certificate linking engine
│   ├── circuit_graph.py         # Circuit DAG construction
│   ├── cross_gene.py            # Cross-gene constraint checker
│   └── tests/
│       ├── test_linker.py
│       ├── test_circuit_graph.py
│       └── test_cross_gene.py
├── mutation/                    # COMP-09: Mutation Explorer
│   ├── __init__.py
│   ├── enumerator.py            # Grammar-guided variant enumeration
│   ├── edit_distance.py         # Edit-distance bounding
│   ├── candidate_ranker.py      # Candidate ranking by constraint score
│   └── tests/
│       ├── test_enumerator.py
│       ├── test_edit_distance.py
│       └── test_candidate_ranker.py
├── orf/                         # COMP-10: ORF Analyzer
│   ├── __init__.py
│   ├── frame_scanner.py         # Six-frame scanning
│   ├── overlap_detector.py      # Overlapping ORF detection
│   ├── interaction.py           # Frame-shift constraint interactions
│   └── tests/
│       ├── test_frame_scanner.py
│       ├── test_overlap_detector.py
│       └── test_interaction.py
├── ir/                          # IR Bus and schemas
│   ├── __init__.py
│   ├── bus.py                   # IR Bus implementation
│   ├── protocol.py              # Bus interface protocol
│   ├── proto/                   # Protocol Buffer schema files
│   │   ├── ir_seq_v1.proto
│   │   ├── ir_peptide_v1.proto
│   │   ├── ir_structure_v1.proto
│   │   └── ir_circuit_v1.proto
│   └── tests/
│       ├── test_bus.py
│       └── test_protocol.py
├── config/                      # Declarative configuration
│   ├── __init__.py
│   ├── loader.py                # Config file loader (YAML/JSON)
│   ├── organisms/               # Organism-specific grammar configs
│   │   ├── homo_sapiens.yaml
│   │   ├── e_coli.yaml
│   │   └── saccharomyces.yaml
│   └── defaults.yaml            # Default configuration
├── cli/                         # Command-line interface
│   ├── __init__.py
│   ├── main.py                  # Click/Typer CLI entry point
│   ├── commands/                # Subcommand implementations
│   │   ├── scan.py
│   │   ├── splice.py
│   │   ├── translate.py
│   │   ├── verify.py
│   │   ├── optimize.py
│   │   ├── certify.py
│   │   ├── compose.py
│   │   ├── explore.py
│   │   └── analyze_orf.py
│   └── tests/
│       └── test_main.py
├── api/                         # Programmatic API (FastAPI)
│   ├── __init__.py
│   ├── app.py                   # FastAPI app + route registration
│   ├── routes.py                # REST endpoints
│   ├── models.py                # Pydantic models for API I/O
│   └── auth.py                  # API key auth
├── proto_out/                   # Generated protobuf Python files (git-ignored)
├── pyproject.toml               # Project metadata, dependencies, tool config
├── Makefile                     # Build and CI targets
└── README.md                    # Developer onboarding
```

### 4.2 Dependency Graph

```
                    ┌─────────┐
                    │  cli/   │
                    └────┬────┘
                         │
                    ┌────▼────┐
                    │  api/   │
                    └────┬────┘
                         │
           ┌─────────────┼─────────────────┐
           │             │                 │
      ┌────▼────┐  ┌─────▼──────┐   ┌──────▼──────┐
      │pipeline │  │  config/   │   │    ir/      │
      └────┬────┘  └────────────┘   └──────┬──────┘
           │                                 │
    ┌──────┼──────┬───────┬──────────┬───────┼──────┬──────────┐
    │      │      │       │          │       │      │          │
┌───▼──┐┌──▼───┐┌─▼────┐┌▼──────┐┌──▼───┐┌──▼───┐┌▼───────┐┌▼──────┐
│scanner││splic.││transl.││  ffi  ││ type ││optim.││certif. ││compos.│
└───┬──┘└──┬───┘└──┬───┘└──┬───┘└──┬───┘└──┬───┘└──┬─────┘└──┬───┘
    │      │       │       │       │       │       │         │
    │      │       │       │       │       │       │         │
┌───▼──────▼───────▼───────▼───────▼───────▼───────▼─────────▼───┐
│                          ir/ (IR Bus)                           │
└────────────────────────────────────────────────────────────────┘
    │                                                             │
┌───▼─────────────────────────────────────────────────────────────▼───┐
│  mutation/  │  orf/                                                │
└────────────────────────────────────────────────────────────────────┘
```

**Dependency rules (enforced by linter):**

1. No component may import another component's internal modules; all
   inter-component communication is through the IR Bus.
2. `ir/` has no dependencies on any component; it is the lowest layer.
3. `config/` is imported only by components that need grammar rules; it has no
   dependency on `ir/`.
4. `cli/` and `api/` are the top-level entry points; they may import any
   component but not vice versa.
5. Circular imports are forbidden; the build system enforces this with a
   static dependency analysis step.

### 4.3 Build and Test

| Artifact                | Command                            | Description                                                             |
|-------------------------|------------------------------------|-------------------------------------------------------------------------|
| Proto stubs             | `make proto`                       | Runs `protoc` on `ir/proto/*.proto`; outputs to `ir/proto_out/`        |
| Type stubs              | `make stubs`                       | Runs `mypy --stubgen` on public API surfaces                           |
| Unit tests              | `make test`                        | Runs `pytest biocompiler/ -v --cov=biocompiler`                        |
| Integration tests       | `make test-integration`            | Runs end-to-end pipeline tests with fixture organisms                  |
| Type checking           | `make typecheck`                   | Runs `mypy biocompiler/ --strict`                                      |
| Linting                 | `make lint`                        | Runs `ruff check biocompiler/`                                         |
| Formatting              | `make format`                      | Runs `ruff format biocompiler/`                                        |
| Build package           | `make build`                       | Runs `python -m build`; produces wheel and sdist                       |
| Certificate test vectors| `make test-vectors`                | Generates and verifies golden certificate test vectors                  |
| Dependency audit        | `make audit`                       | Runs `pip-audit` on pinned dependencies                                |

**CI pipeline:** On every pull request, the following must pass:

1. `make proto` (no stale proto stubs)
2. `make lint` (zero violations)
3. `make typecheck` (zero errors)
4. `make test` (100% pass, ≥90% line coverage)
5. `make test-integration` (100% pass)
6. `make audit` (no known vulnerabilities)

### 4.4 Coding Standards

| Aspect             | Standard                                                                                     |
|--------------------|----------------------------------------------------------------------------------------------|
| Language           | Python 3.12+ with type hints on all public APIs                                             |
| Type system        | Strict `mypy`; no `Any` in public interfaces; `Protocol` for structural typing               |
| Docstrings         | Google-style docstrings on all public classes and functions                                  |
| Naming             | `snake_case` for functions/variables; `PascalCase` for classes; `UPPER_SNAKE` for constants  |
| Immutability       | Prefer `@dataclass(frozen=True)` for IR records; no in-place mutation of IR data             |
| Error handling     | Use custom exception hierarchy; never catch bare `Exception`; always chain causes            |
| Logging            | `structlog` with JSON output; no `print()` in library code                                  |
| Imports            | `__all__` in every `__init__.py`; no wildcard imports                                        |
| Testing            | pytest; property-based testing via Hypothesis for DFA/FST engines; ≥90% line coverage        |
| Security           | No `eval()` / `exec()`; no shell injection in FFI calls; all file I/O via `pathlib`         |

---

## 5. Physical View

### 5.1 Deployment Model

```
 ┌─────────────────────────────────────────────────────────────────┐
 │                        Workstation                               │
 │                                                                  │
 │  ┌───────────────────────────────────────────────────────────┐  │
 │  │                  BioCompiler Core Pipeline                 │  │
 │  │                                                           │  │
 │  │   COMP-01 → COMP-02 → COMP-03 → COMP-05 → COMP-06 →     │  │
 │  │   COMP-07 → COMP-08 → COMP-09 → COMP-10                  │  │
 │  │                                                           │  │
 │  │   (CPU-bound: pure Python, deterministic, no GPU needed)  │  │
 │  └───────────────────────────┬───────────────────────────────┘  │
 │                              │                                   │
 │                              │ FFI calls (subprocess / gRPC)     │
 │                              │                                   │
 │  ┌───────────────────────────▼───────────────────────────────┐  │
 │  │              COMP-04: FFI Manager                          │  │
 │  │                                                           │  │
 │  │  ┌─────────────────┐  ┌─────────────────┐                 │  │
 │  │  │ Folding Oracle  │  │   PTM Oracle    │                 │  │
 │  │  │ (GPU-accelerated│  │  (CPU-based ML) │                 │  │
 │  │  │  e.g. AlphaFold)│  │                 │                 │  │
 │  │  └────────┬────────┘  └─────────────────┘                 │  │
 │  │           │                                                │  │
 │  └───────────┼────────────────────────────────────────────────┘  │
 │              │                                                   │
 │  ┌───────────▼───────────────────────────────────────────────┐  │
 │  │                      GPU (optional)                        │  │
 │  │                                                           │  │
 │  │   CUDA / ROCm device for folding oracle inference         │  │
 │  │   Min: 8 GB VRAM; Recommended: 24 GB VRAM                │  │
 │  └───────────────────────────────────────────────────────────┘  │
 │                                                                  │
 │  ┌───────────────────────────────────────────────────────────┐  │
 │  │                    Local Storage                           │  │
 │  │                                                           │  │
 │  │   ~/.biocompiler/                                         │  │
 │  │   ├── configs/        # Organism grammar configs          │  │
 │  │   ├── cache/          # FFI result cache (SQLite)         │  │
 │  │   ├── certificates/   # Generated certificates (JSON)     │  │
 │  │   └── logs/           # Structured JSON logs              │  │
 │  └───────────────────────────────────────────────────────────┘  │
 └─────────────────────────────────────────────────────────────────┘
```

**Deployment configurations:**

| Configuration   | CPU           | RAM    | GPU              | Use Case                              |
|-----------------|---------------|--------|------------------|---------------------------------------|
| Minimal         | 4 cores       | 8 GB   | None (CPU fallback for folding) | Development, small genes |
| Standard        | 8 cores       | 16 GB  | 1× 8 GB VRAM     | Single-gene production runs           |
| High-throughput | 32 cores      | 64 GB  | 2× 24 GB VRAM    | Multi-gene circuits, batch processing |

### 5.2 Data Storage

| Data                  | Format          | Location                            | Size Estimate (per gene)   | Retention   |
|-----------------------|-----------------|-------------------------------------|----------------------------|-------------|
| Input DNA sequence    | FASTA / plain   | CLI argument or file path           | 1 KB – 10 MB               | User-managed|
| IR-Seq records        | Protobuf binary | In-memory (IR Bus) + temp file      | 10 KB – 100 MB (multi-iso) | Run lifetime|
| IR-Peptide records    | Protobuf binary | In-memory (IR Bus) + temp file      | 5 KB – 50 MB               | Run lifetime|
| IR-Structure records  | Protobuf binary | In-memory (IR Bus) + temp file      | 100 KB – 500 MB            | Run lifetime|
| IR-Circuit records    | Protobuf binary | In-memory (IR Bus) + temp file      | 1 KB – 10 MB (per edge)    | Run lifetime|
| Certificates          | JSON            | `~/.biocompiler/certificates/`      | 10 KB – 1 MB               | Indefinite  |
| FFI cache             | SQLite          | `~/.biocompiler/cache/ffi_cache.db` | 100 MB – 10 GB (cumulative)| 30-day TTL  |
| Organism configs      | YAML            | `biocompiler/config/organisms/`     | 1 KB – 50 KB per organism  | Version-controlled|
| Logs                  | JSON lines      | `~/.biocompiler/logs/`              | 1 MB – 100 MB (cumulative) | 7-day rotation|

### 5.3 Network Requirements

| Scenario                        | Protocol  | Data Volume          | Latency Requirement        |
|---------------------------------|-----------|----------------------|----------------------------|
| FFI call to local folding oracle| gRPC      | 1–50 MB per call     | 10 s – 10 min (batching OK)|
| FFI call to local PTM oracle    | gRPC      | 100 KB – 5 MB/call   | 1 s – 30 s                 |
| FFI call to remote folding API  | HTTPS/REST| 1–50 MB per call     | 30 s – 30 min              |
| Certificate download            | File I/O  | 10 KB – 1 MB         | N/A (local)                |
| Config update check             | HTTPS/REST| < 100 KB             | < 5 s                      |

**Network isolation:** The core pipeline (COMP-01 through COMP-03, COMP-05
through COMP-10) requires **no network access**. Only COMP-04 (FFI Manager)
initiates network connections. This design ensures:

- Air-gapped deployments are possible for sensitive sequences.
- Firewall rules can restrict egress to only the FFI oracle endpoints.
- FFI calls can be replaced with cached results for fully offline operation.

---

## 6. Scenario View (Use Cases)

### UC-01: Design a Splicing-Safe Gene

**Actor:** Synthetic biologist
**Goal:** Design a gene whose splicing produces no unwanted isoforms and that
satisfies all structural and functional constraints.

| Step | Stage           | Action                                                                                         | Output                           |
|------|-----------------|------------------------------------------------------------------------------------------------|----------------------------------|
| 1    | COMP-01 Scanner | Submit DNA sequence to scanner; receive tokenized IR-Seq with annotated regulatory motifs.    | IR-Seq (tokenized)               |
| 2    | COMP-02 Splicing| Run NDFST simulation; enumerate all possible isoforms.                                        | Set of IR-Seq records (isoforms) |
| 3    | COMP-05 Type    | Check splicing constraints on each isoform: are cryptic splice sites present? Unwanted exons?  | Verdict per isoform              |
| 4    | COMP-03 Transl. | Translate each PASS-ing isoform to IR-Peptide.                                                | Set of IR-Peptide records        |
| 5    | COMP-04 FFI     | Submit peptides to folding oracle; retrieve structure predictions.                             | Set of IR-Structure records      |
| 6    | COMP-05 Type    | Check structural constraints: stability, aggregation propensity, active-site geometry.         | Verdict per isoform              |
| 7    | COMP-06 Optim.  | If any isoform FAILs, run CSP solver to find sequence modifications that eliminate the unwanted isoform or fix the structural issue. MUS diagnosis explains failure. | Optimized IR-Structure or MUS report |
| 8    | COMP-07 Certif. | Generate JSON certificate encoding the final verdict, all evidence, and (if applicable) MUS.  | Certificate (JSON)               |
| 9    | Review          | Synthetic biologist reviews certificate; if PASS, proceeds to lab synthesis. If FAIL, iterates with COMP-09 (Mutation Explorer). | Decision                         |

### UC-02: Verify a Multi-Gene Circuit

**Actor:** Circuit designer
**Goal:** Compose multiple pre-verified genes into a circuit and verify that
cross-gene interactions do not violate system-level constraints.

| Step | Stage               | Action                                                                                                    | Output                        |
|------|---------------------|-----------------------------------------------------------------------------------------------------------|-------------------------------|
| 1    | Collect             | Gather per-gene certificates for all genes in the circuit. Each certificate is independently verified.   | Set of valid certificates     |
| 2    | COMP-08 Circuit     | Build circuit DAG from topology specification; annotate edges with interaction types.                     | IR-Circuit (graph)            |
| 3    | COMP-08 Cross-gene  | Check cross-gene constraints: promoter interference, metabolic burden, resource competition, regulatory crosstalk. | Cross-gene verdicts           |
| 4    | COMP-07 Certificate | Compose per-gene certificates and cross-gene verdicts into a circuit-level certificate.                   | Circuit certificate (JSON)    |
| 5    | Review              | Circuit designer reviews certificate; if PASS, submits to foundry. If FAIL, uses MUS to identify minimal set of problematic interactions. | Decision                      |

### UC-03: Analyze Overlapping Reading Frames

**Actor:** Molecular biologist
**Goal:** Identify and characterize overlapping ORFs in a viral or compact
genome to understand frame-shifted coding potential and constraint interactions.

| Step | Stage               | Action                                                                                                          | Output                          |
|------|---------------------|-----------------------------------------------------------------------------------------------------------------|----------------------------------|
| 1    | COMP-01 Scanner     | Submit genome sequence; receive tokenized IR-Seq with all six reading frames annotated.                         | IR-Seq (six-frame)               |
| 2    | COMP-10 ORF Scanner | Scan all six frames for start/stop codon pairs; identify candidate ORFs.                                        | ORF candidate list               |
| 3    | COMP-10 Overlap     | Detect overlapping ORFs across frames; compute overlap regions and frame-shift relationships.                   | Overlap map                      |
| 4    | COMP-10 Interaction | Analyze constraint interactions between overlapping ORFs: shared nucleotides impose coupled constraints.         | Interaction constraints          |
| 5    | COMP-05 Type        | Check whether overlapping ORFs are jointly satisfiable; assign five-valued verdicts per ORF and per overlap.   | Verdicts per ORF and overlap     |
| 6    | COMP-07 Certificate | Generate certificate encoding ORF map, overlap relationships, and verdicts.                                     | ORF analysis certificate (JSON)  |

---

## 7. Architecture Quality Attributes

| Attribute      | How Achieved                                                                                                         | How Verified                                                                                    |
|----------------|----------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| Determinism    | All core stages (COMP-01, COMP-02, COMP-03, COMP-05, COMP-06, COMP-07, COMP-08, COMP-09, COMP-10) are pure functions with no hidden state. FFI non-determinism is captured as UNCERTAIN verdicts. | Property-based tests (Hypothsis) with fixed seeds assert bit-identical outputs across repeated runs. Golden test vectors for every pipeline stage. |
| Soundness      | Type system (COMP-05) enforces that no PASS verdict is issued when a constraint is violated. MUS diagnosis (COMP-06) guarantees that every failure has a provably minimal explanation. Certificate verifier (COMP-07) independently re-checks all verdicts. | Formal proof of type-system soundness for five-valued logic. Mutation testing: inject constraint violations and verify FAIL verdicts are always produced. |
| Completeness   | NDFST (COMP-02) exhaustively enumerates all isoforms, not just the most probable. CSP solver (COMP-06) is complete: if a solution exists, it will be found. | Isoform enumeration: compare NDFST output against brute-force enumeration on small inputs. CSP: verify solver finds all solutions on benchmark constraint sets. |
| Composability  | Compositional verifier (COMP-08) links per-gene certificates without re-running the pipeline. Circuit certificates are independently verifiable. | Integration tests with circuits of 2–20 genes; verify circuit certificate can be validated by standalone verifier with no access to pipeline. |
| Extensibility  | Declarative grammar configs (AD-07) for new organisms. FFI adapter pattern for new oracles. Constraint definitions are data-driven, not hard-coded. | Add a new organism config without modifying any source code; run full test suite. Add a mock FFI adapter; verify it integrates via the adapter protocol. |
| Performance    | Structured concurrency allows parallel processing of independent genes and isoforms. FFI result caching avoids redundant oracle calls. Protocol Buffers provide efficient serialization. | Benchmark suite: single gene (1 kb, 10 kb, 100 kb), multi-gene circuit (5, 20, 100 genes). Profile FFI cache hit rates. Target: <60 s for a 10 kb gene on Standard deployment. |

---

*End of DOC-02: Software Architecture Document (SAD)*


---

## Addendum: Corrected Module Tree (Actual Implementation)

The following reflects the actual Python package structure in `src/biocompiler/`:

```
src/biocompiler/
├── __init__.py
├── shared/
│   ├── types.py              # Core data structures (Token, SpliceIsoform, TypeCheckResult, Certificate, Verdict)
│   ├── five_valued_logic.py  # Five-valued verdict logic (PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL)
│   └── constants.py
├── sequence/
│   ├── scanner.py            # Multi-DFA motif detection (was scanner/ package)
│   └── splicing.py           # NDFST isoform computation (was splicing/ package)
├── expression/
│   └── translation.py        # Codon-to-amino-acid FST + CAI/tAI computation (was translation/ package)
├── type_system/              # Predicate registry + evaluator functions
│   ├── __init__.py
│   ├── codon_tables.py
│   ├── checks.py
│   ├── predicates.py
│   └── registry.py
├── engines/                  # Analysis engine adapters (was ffi/ package)
│   ├── base.py               # BaseEngineResult, MutationResult, BatchResult
│   ├── foldx.py              # FoldX protein stability
│   ├── camsol.py             # CamSol solubility
│   ├── esmfold.py            # ESMFold structure prediction
│   └── immunogenicity.py     # MHC binding + B-cell epitope
├── provenance/               # Certificate + SLOT verification (was certificate/ package)
│   ├── __init__.py
│   ├── certificate.py        # Graduated certificate generation + verification
│   └── slot_verification.py  # Three-mode SLOT verification
├── optimizer/                # Integrated constraint-solving optimizer + mutagenesis loop
│   ├── __init__.py
│   ├── pipeline_core.py        # Main optimization pipeline (default fast path)
│   ├── integrated_optimizer.py # Single-pass integrated constraint-solving optimizer (the only optimizer since 0.9.1)
│   ├── pipeline_certification.py # Certified-by-default predicate evaluation
│   ├── pipeline_paths.py       # Extended predicate evaluation
│   ├── mutagenesis.py        # Type-directed mutagenesis
│   └── assembly.py
├── organisms/                # Organism-specific data (25 organisms + tAI)
│   ├── config.py
│   └── ...
├── biosecurity/              # Biosecurity screening
├── grammars/                 # Declarative grammar configuration (was config/ with YAML)
│   └── ...
├── tai.py                    # tRNA Adaptation Index
├── api.py                    # REST API (FastAPI)
└── cli.py                    # Command-line interface
```

**Not implemented:** Protocol Buffer IR schemas (`ir/proto/`), `orf/` package, `compositional/` package. Python dataclasses in `shared/types.py` replace the planned IR Bus and Protocol Buffer schemas.
