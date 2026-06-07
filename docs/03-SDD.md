# DOC-03: Software Design Document (SDD)

| Field | Value |
|---|---|
| **Document ID** | DOC-03 |
| **Version** | 12.0.0 |
| **Status** | Current |
| **Date** | 2026-06-07 |

---

## 1. Introduction

### 1.1 Purpose

This document specifies the detailed software design for each component of the BioCompiler system. It provides sufficient detail for a developer to implement each component correctly: algorithms, data structures, invariants, error handling, and inter-component contracts. This document bridges the gap between the architectural decomposition (DOC-02) and the interface contracts (DOC-04).

### 1.2 Scope

Covers the detailed design of all ten components (COMP-01 through COMP-10) and the IR Bus. Does not cover the design of external tools invoked via FFI (AlphaFold, NetPhos, etc.) — these are treated as black boxes with defined input/output contracts.

### 1.3 References

| Ref ID | Document |
|---|---|
| REF-01 | DOC-01: Software Requirements Specification |
| REF-02 | DOC-02: Software Architecture Document |
| REF-03 | DOC-04: Interface Control Document |
| REF-04 | DOC-06: Design Rationale |

---

## 2. Design Methodology

### 2.1 Design Principles

| Principle | Description | Enforcement |
|---|---|---|
| **Determinism by Construction** | No randomness, no sampling, no probabilistic inference in any internal component | Code review + reproducibility tests |
| **Soundness by Proof** | Every PASS verdict is backed by a derivable proof | Soundness argument per type predicate + adversarial testing |
| **Separation of Paradigms** | Symbolic (grammar-based) and continuous (ML/physics) computations never mix in the same component | Architecture constraint: FFI boundary isolates them |
| **Fail Explicitly** | Every error condition has a named exception type with diagnostic information | Exception hierarchy + error propagation tests |
| **Design for Testability** | Every component is a pure function with defined input/output contract | Unit test coverage ≥ 90%; 100% for COMP-05, COMP-06 |

### 2.2 Design Patterns

| Pattern | Application | Component(s) |
|---|---|---|
| **Pipeline** | Staged transformation with typed IR | Core architecture |
| **Strategy** | FFI adapters are interchangeable strategies for external tool invocation | COMP-04 |
| **Visitor** | Type predicates visit the IR to produce verdicts | COMP-05 |
| **Builder** | IR objects are constructed incrementally through the pipeline | IR Bus |
| **Factory Method** | Certificate generator creates different certificate types | COMP-07 |
| **Composite** | Circuit verification composes individual gene verifications | COMP-08 |

---

## 3. Detailed Component Designs

### 3.1 COMP-01: Scanner

#### 3.1.1 Overview

The Scanner performs lexical analysis of nucleotide sequences. It implements a collection of deterministic finite automata (DFAs), one per motif type, that scan the input sequence simultaneously and produce an ordered stream of annotated tokens. The scanner is the entry point of the pipeline; its output feeds the Splicing Engine (COMP-02) and populates the IR-Seq.

#### 3.1.2 Class Structure

```
Scanner
  ├── DFAEngine
  │     ├── active_dfas: list[DFA]
  │     ├── run(sequence: str) -> list[Token]
  │     └── resolve_overlaps(tokens: list[Token]) -> list[Token]
  ├── MotifDefinition
  │     ├── pattern: str           # IUPAC regex
  │     ├── motif_type: str        # "start_codon" | "stop_codon" | "splice_donor" | ...
  │     ├── threshold: float       # Minimum match score (for PWMs)
  │     └── to_dfa() -> DFA
  └── Token
        ├── position: int
        ├── element_type: str
        ├── match_sequence: str
        └── score: float
```

#### 3.1.3 Algorithm: Multi-DFA Scanning

```
ALGORITHM MultiDFAScan:
  INPUT:  nucleotide sequence S, list of MotifDefinitions M
  OUTPUT: ordered list of Tokens T

  1. For each MotifDefinition m in M:
       Construct DFA_d = m.to_dfa()
       Add DFA_d to active_dfas

  2. Initialize: position p = 0, token_list T = []

  3. For p = 0 to len(S) - 1:
       For each DFA_d in active_dfas:
         Feed S[p] to DFA_d
         If DFA_d is in accepting state:
           Emit Token(position=p-start, type=d.motif_type,
                      sequence=S[p-start:p+1], score=d.current_score)
         If DFA_d is in dead state:
           Reset DFA_d to start state

  4. Resolve overlapping tokens (longest-match rule):
     If two tokens overlap in position range:
       Keep the longer match; discard the shorter
     If same length: keep higher score

  5. Sort T by position (ascending), then by element_type priority
  6. Return T
```

**Time Complexity:** O(n × d × k) where n = sequence length, d = number of DFAs, k = average DFA size. For typical inputs (n ≈ 10,000, d ≈ 9, k ≈ 20), this is approximately 1.8M operations — well within the 1-second requirement.

**Space Complexity:** O(n + t) where t = number of tokens emitted.

#### 3.1.4 Motif Definitions

| Motif Type | Pattern (IUPAC) | Notes |
|---|---|---|
| Start codon | `ATG` | Reading frame must be tracked |
| Stop codon | `TAA\|TAG\|TGA` | Reading frame must be tracked |
| Splice donor | `(C\|A)AGGT(A\|G)AGT` | GT-AG rule; consensus extends beyond GT |
| Splice acceptor | `(C\|T)AG\|G` | AG at 3' intron boundary |
| Branch point | `YNYRAY` | 18–40 nt upstream of acceptor; fuzzy match |
| Polypyrimidine tract | `[CT]{8,}` | Minimum 8 nt of C/T upstream of acceptor |
| Kozak consensus | `GCCACCATGG` | Translation initiation context |
| AUUUA instability | `ATTTA` | mRNA destabilizing motif |
| Restriction site | User-specified | From REBASE; exact match |

#### 3.1.5 Invariants

| ID | Invariant |
|---|---|
| INV-SCAN-01 | Tokens are ordered by position (ascending) |
| INV-SCAN-02 | No token's position range exceeds sequence bounds |
| INV-SCAN-03 | Each position is examined by all active DFAs |
| INV-SCAN-04 | Determinism: identical input produces identical token list |

#### 3.1.6 Error Handling

| Error | Condition | Diagnostic |
|---|---|---|
| `InvalidSequenceError` | Non-IUPAC character in input | Position + character |
| `FormatError` | Malformed FASTA header or multi-sequence FASTA | Line number |
| `EmptySequenceError` | Zero-length input | — |

---

### 3.2 COMP-02: Splicing Engine

#### 3.2.1 Overview

The Splicing Engine implements a Non-Deterministic Finite-State Transducer (NDFST) that parses the token stream from the Scanner against the splicing grammar and computes the set of all possible splice isoforms. The NDFST captures alternative splicing as non-deterministic branching: each valid parse path through the transducer corresponds to a distinct splice isoform. The computation is deterministic (same isoform set for same input), but the output is set-valued (multiple isoforms possible).

#### 3.2.2 Class Structure

```
SplicingEngine
  ├── NDFST
  │     ├── states: set[State]
  │     ├── transitions: dict[State, list[Transition]]
  │     ├── initial_state: State
  │     ├── accepting_states: set[State]
  │     ├── execute(token_stream, context) -> list[Isoform]
  │     └── subset_construction(token_stream) -> list[ParsePath]
  ├── GrammarRule
  │     ├── rule_type: str        # "donor" | "acceptor" | "branch_point" | ...
  │     ├── pattern: str          # Consensus sequence
  │     ├── context_requirements: dict  # Cell-type-dependent thresholds
  │     ├── confidence: float
  │     └── is_applicable(context) -> bool
  ├── CellularContext
  │     ├── cell_type: str
  │     ├── ese_threshold: float
  │     ├── ess_threshold: float
  │     ├── ise_threshold: float
  │     ├── iss_threshold: float
  │     └── splicing_factor_concentrations: dict
  └── Isoform
        ├── sequence: str
        ├── exon_boundaries: list[PositionRange]
        ├── reading_frame: ReadingFrame
        ├── confidence: float
        └── parse_path: list[GrammarRule]
```

#### 3.2.3 Algorithm: NDFST Splicing Parse

```
ALGORITHM NDFSTSplice:
  INPUT:  token_stream T, cellular_context C, grammar_rules G
  OUTPUT: set of Isoforms I

  1. Construct NDFST from grammar rules G:
     - Each donor site token creates a "splice out" transition
     - Each acceptor site token creates a "splice in" transition
     - Exon regions are "emit" transitions (output the exon sequence)
     - ESE/ESS/ISE/ISS tokens modulate transition probabilities (as thresholds)

  2. Execute NDFST via subset construction:
     current_states = {initial_state}
     parse_paths = [{path: [], output: ""}]

  3. For each token t in T:
     next_states = {}
     next_paths = []

     For each state s in current_states:
       For each transition (s -> s', output) applicable to t under context C:
         Add s' to next_states
         For each parse_path ending in s:
           Create new path: append (s', output, t)
           Add to next_paths

     current_states = next_states
     parse_paths = next_paths

     If current_states is empty:
       BACKTRACK or ABORT (no valid parse)

  4. Filter parse_paths: keep only those ending in accepting states
  5. For each surviving parse_path:
       Construct Isoform: concatenate emitted exon sequences,
       record exon boundaries, reading frame, confidence
  6. Return set of Isoforms I
```

**Key Property:** The NDFST explores ALL valid parse paths, so the isoform set is complete (no valid isoform is missed) and sound (every isoform in the set satisfies the grammar).

#### 3.2.4 Cellular Context Parameterization

The cellular context modulates which transitions are enabled in the NDFST:

| Context Parameter | Effect on NDFST |
|---|---|
| ESE threshold | Below threshold: ESE-enhanced donor sites are downgraded (less likely to be used as splice donors) |
| ESS threshold | Below threshold: ESS-silenced sites are upgraded (more likely to be used) |
| Splicing factor concentrations | High SF2/ASF concentration → favors exon inclusion; High hnRNP A1 → favors exon skipping |

These modulations do not eliminate transitions (which would risk missing valid isoforms) but adjust the confidence scores of the resulting isoforms. The isoform set remains complete; the confidence scores help users prioritize.

#### 3.2.5 Invariants

| ID | Invariant |
|---|---|
| INV-SPL-01 | Isoform set is complete: all valid parse paths are represented |
| INV-SPL-02 | Isoform set is sound: every isoform satisfies the splicing grammar |
| INV-SPL-03 | Computation is deterministic: same input + context → same isoform set |
| INV-SPL-04 | Each isoform has a unique parse path (no duplicate isoforms) |

---

### 3.3 COMP-03: Translation Engine

#### 3.3.1 Overview

The Translation Engine implements a deterministic finite-state transducer (FST) that maps each codon in a spliced mRNA sequence to its corresponding amino acid according to the standard genetic code. The FST handles three special cases: selenocysteine insertion, pyrrolysine insertion, and programmed ribosomal frameshifting.

#### 3.3.2 Algorithm: Codon-by-Codon Translation

```
ALGORITHM Translate:
  INPUT:  IR-Seq (single splice isoform)
  OUTPUT: IR-Peptide

  1. Extract coding sequence from IR-Seq:
     CDS = IR-Seq.sequence[IR-Seq.exon_boundaries concatenated]
     frame = IR-Seq.reading_frame

  2. Initialize:
     amino_acids = []
     codon_assignments = []
     sec_flags = []
     frameshift_warnings = []
     position = frame  # Start at reading frame offset

  3. While position + 2 < len(CDS):
     codon = CDS[position:position+3]

     # Special case: Selenocysteine
     If codon == "TGA" and SECIS_element_present(IR-Seq):
       aa = "U"  # Selenocysteine
       sec_flags.append(SelenocysteineFlag(position, True))
     # Special case: Pyrrolysine (archaeal context only)
     Else if codon == "TAG" and context.organism in ARCHAEAL_ORGANISMS:
       aa = "O"  # Pyrrolysine
     # Standard genetic code
     Else:
       aa = CODON_TABLE[codon]

     amino_acids.append(aa)
     codon_assignments.append(CodonAssignment(position, codon, aa))

     # Check for frameshift motifs at this position
     If matches_frameshift_motif(CDS, position):
       frameshift_warnings.append(FrameshiftWarning(position, motif, direction))

     position += 3

  4. Construct IR-Peptide:
     - amino_acid_sequence = "".join(amino_acids)
     - codon_assignment = codon_assignments
     - selenocysteine_flags = sec_flags
     - frameshift_warnings = frameshift_warnings
     - reading_frame = frame
     - source_isoform_id = IR-Seq.id

  5. Return IR-Peptide
```

**Time Complexity:** O(n) where n = length of coding sequence. Trivial for any biological sequence.

#### 3.3.3 Invariants

| ID | Invariant |
|---|---|
| INV-TRA-01 | Amino acid sequence length = number of sense codons processed |
| INV-TRA-02 | Every codon assignment maps the specified codon to the correct amino acid per the genetic code |
| INV-TRA-03 | Determinism: same input → same output (modulo frameshift warnings, which are deterministic annotations) |

---

### 3.4 COMP-04: FFI Manager

#### 3.4.1 Overview

The FFI Manager orchestrates invocations of external tools (protein structure predictors, PTM predictors) through a common adapter interface. It is responsible for input formatting, output parsing, invariant validation, provenance tracking, and error handling. The FFI Manager treats every external tool as a non-deterministic black box.

#### 3.4.2 Algorithm: SLOT-Fill Protocol

```
ALGORITHM SlotFill:
  INPUT:  IR-Peptide (with empty SLOT fields), adapter_name, adapter_config
  OUTPUT: IR-Peptide or IR-Structure (with SLOT fields filled)

  1. Select adapter A from registry by adapter_name
  2. If adapter not found: raise AdapterNotFoundError

  3. Format input for external tool:
     input_data = A.format_input(ir_peptide, adapter_config)

  4. Invoke external tool:
     try:
       raw_output = A.invoke(input_data, timeout=adapter_config.timeout)
     except TimeoutError:
       raise ExternalToolError("Tool exceeded timeout")
     except SubprocessError as e:
       raise ExternalToolError(f"Tool returned non-zero exit: {e}")

  5. Parse output:
     try:
       parsed = A.parse_output(raw_output)
     except ParseError:
       raise OutputParseError("Output format does not match expected schema")

  6. Validate invariants:
     If not A.validate_output(parsed):
       raise OutputValidationError("Parsed output violates IR invariants")

  7. Fill SLOT fields in IR:
     A.fill_slots(ir_record, parsed)

  8. Record provenance:
     ir_record.provenance = Provenance(
       tool=A.name(), version=A.version(),
       timestamp=now(), parameters=adapter_config,
       input_hash=sha256(ir_peptide.amino_acid_sequence))

  9. Return enriched IR record
```

#### 3.4.3 Adapter Contract

Every FFI adapter MUST implement:

```python
class FFIAdapter(ABC):
    @abstractmethod
    def name(self) -> str: ...

    @abstractmethod
    def slot_fields(self) -> list[str]: ...

    @abstractmethod
    def format_input(self, ir_peptide: IRPeptide, config: dict) -> bytes: ...

    @abstractmethod
    def invoke(self, input_data: bytes, timeout: float) -> bytes: ...

    @abstractmethod
    def parse_output(self, raw_output: bytes) -> dict: ...

    @abstractmethod
    def fill_slots(self, ir_record, parsed: dict) -> None: ...

    @abstractmethod
    def validate_output(self, output) -> bool: ...
```

#### 3.4.4 Invariants

| ID | Invariant |
|---|---|
| INV-FFI-01 | Every SLOT field is either empty or filled with valid data (never partially filled) |
| INV-FFI-02 | Provenance metadata is recorded for every FFI invocation |
| INV-FFI-03 | FFI output invariants (INV-STR-01, INV-STR-02, INV-PEP-03) are validated before acceptance |

---

### 3.5 COMP-05: Type System

#### 3.5.1 Overview

The Type System performs static verification of biological correctness properties on mRNA sequences. It evaluates each declared type predicate against the mRNA and produces a three-valued verdict: PASS (guaranteed correct), FAIL (guaranteed incorrect), or UNCERTAIN (cannot determine). The Type System is the safety-critical component of BioCompiler; its soundness (no false PASS) is the most important invariant of the entire system.

#### 3.5.2 Algorithm: Type Checking

```
ALGORITHM TypeCheck:
  INPUT:  IR-Seq, list of TypePredicates, target_isoform (optional)
  OUTPUT: list of TypeCheckResults

  1. Initialize results = []

  2. For each TypePredicate P in predicates:
     verdict = P.evaluate(ir_seq, target_isoform)

     If verdict == PASS:
       derivation = P.get_derivation()  # Chain of reasoning
       results.append(TypeCheckResult(P, PASS, derivation, None, None))

     Else if verdict == FAIL:
       violation = P.get_violation()    # Position + rule + evidence
       results.append(TypeCheckResult(P, FAIL, None, violation, None))

     Else:  # UNCERTAIN
       gap = P.get_knowledge_gap()      # What's needed to resolve
       results.append(TypeCheckResult(P, UNCERTAIN, None, None, gap))

  3. Compose results using three-valued logic:
     overall = PASS
     For each result in results:
       overall = three_valued_and(overall, result.verdict)

  4. Return results
```

#### 3.5.3 Three-Valued Logic Truth Tables

**Conjunction (AND):**

| AND | PASS | UNCERTAIN | FAIL |
|---|---|---|---|
| **PASS** | PASS | UNCERTAIN | FAIL |
| **UNCERTAIN** | UNCERTAIN | UNCERTAIN | FAIL |
| **FAIL** | FAIL | FAIL | FAIL |

**Disjunction (OR):**

| OR | PASS | UNCERTAIN | FAIL |
|---|---|---|---|
| **PASS** | PASS | PASS | PASS |
| **UNCERTAIN** | PASS | UNCERTAIN | UNCERTAIN |
| **FAIL** | PASS | UNCERTAIN | FAIL |

#### 3.5.4 Soundness Arguments per Predicate

| Predicate | Soundness Argument |
|---|---|
| `SpliceCorrect(C)` | PASS only if NDFST output set is a singleton containing exactly the target isoform. Sound because NDFST is complete (all valid parse paths explored) and the isoform is verified by re-parsing. |
| `NoCrypticSplice` | PASS only if no DFA scan match exceeds the splice site strength threshold. Sound because the scanner is exhaustive (all positions scanned) and the threshold is a conservative upper bound for non-functional sites. |
| `CodonAdapted(O, t)` | PASS only if computed CAI ≥ t. Sound because CAI is a deterministic function of the sequence and the organism's codon usage table. |
| `GCInRange(lo, hi)` | PASS only if lo ≤ GC% ≤ hi. Sound because GC% is a countable property of the sequence. |
| `NoRestrictionSite(S)` | PASS only if no exact match for any enzyme in S. Sound because exact string matching is deterministic and complete. |
| `InFrame` | PASS only if reading frame is consistent and no premature stop codons exist. Sound because frame consistency and stop codon detection are deterministic. |
| `NoInstabilityMotif` | PASS only if no AUUUA or U-rich motifs detected. Sound because motif scanning is exhaustive. |

#### 3.5.5 Invariants

| ID | Invariant |
|---|---|
| INV-TYP-01 | **Soundness**: If PASS is returned, the property holds for the sequence under the specified conditions |
| INV-TYP-02 | Determinism: same input → same verdicts |
| INV-TYP-03 | Three-valued composition preserves soundness |
| INV-TYP-04 | Subtyping: if SpliceCorrect(C1) passes and C1 <: C2, then SpliceCorrect(C2) passes |

---

### 3.6 COMP-06: Optimizer (CSP)

#### 3.6.1 Overview

The Optimizer solves a Constraint Satisfaction Problem (CSP) to find synonymous codon assignments that satisfy all hard constraints while maximizing a scalar objective (typically CAI). When the CSP is infeasible, it computes a Minimal Unsatisfiable Subset (MUS) to diagnose which constraints conflict.

#### 3.6.2 Algorithm: CSP Solving with Constraint Propagation

```
ALGORITHM CSPOptimize:
  INPUT:  IR-Seq, list of Constraints, Objective
  OUTPUT: FeasibleResult or InfeasibleResult

  1. Initialize CSP:
     For each codon position i (1..n):
       variables[i] = CSPVariable(
         position=i,
         domain=SYNONYMOUS_CODONS[amino_acid_at_i])

  2. Apply constraint propagation (AC-3):
     changed = True
     While changed:
       changed = False
       For each constraint C:
         For each variable v affected by C:
           original_domain = v.domain
           v.domain = C.prune(v, variables)
           If v.domain is empty:
             Return InfeasibleResult(mus=compute_mus(C, variables))
           If v.domain != original_domain:
             changed = True

  3. If all domains are singletons:
     assignment = [v.domain[0] for v in variables]
     If satisfies_all_constraints(assignment, constraints):
       Return FeasibleResult(assignment, objective_value)
     Else:
       Return InfeasibleResult(mus=compute_mus_all(constraints))

  4. Backtracking search:
     result = backtrack(variables, constraints, 0)
     If result is None:
       Return InfeasibleResult(mus=compute_mus_all(constraints))
     Else:
       Return FeasibleResult(result, compute_objective(result, objective))

ALGORITHM Backtrack(variables, constraints, index):
  If index == len(variables):
    If satisfies_all_constraints(assignment, constraints):
      Return assignment
    Else:
      Return None

  For each codon c in variables[index].domain (ordered by CAI, descending):
    variables[index].assignment = c
    If local_constraints_satisfied(index, variables, constraints):
      result = Backtrack(variables, constraints, index + 1)
      If result is not None:
        Return result
    variables[index].assignment = None

  Return None
```

#### 3.6.3 MUS Computation

```
ALGORITHM ComputeMUS:
  INPUT:  set of constraints C (all unsatisfiable together)
  OUTPUT: minimal unsatisfiable subset of C

  1. mus = copy(C)
  2. For each constraint c in C:
       temp = mus - {c}
       If not satisfiable(temp):  # Still unsatisfiable without c
         mus = temp              # c is not in MUS; remove it
       # Else: c is necessary for unsatisfiability; keep it
  3. Return mus
```

#### 3.6.4 Invariants

| ID | Invariant |
|---|---|
| INV-CSP-01 | Solution satisfies ALL hard constraints (verified by re-running type checker) |
| INV-CSP-02 | MUS is truly unsatisfiable (verified independently) |
| INV-CSP-03 | Determinism: same input → same solution or same INFEASIBLE report |
| INV-CSP-04 | Completeness: if a feasible assignment exists, solver finds one |

---

### 3.7 COMP-07: Certificate Generator

#### 3.7.1 Algorithm: Certificate Generation

```
ALGORITHM GenerateCertificate:
  INPUT:  IR-Seq, TypeCheckResults, OptimizationResult, CompositionResults (optional)
  OUTPUT: GuaranteeCertificate (JSON)

  1. Verify preconditions:
     Assert all TypeCheckResults have verdict PASS
     Assert OptimizationResult is FeasibleResult
     If CompositionResults provided:
       Assert all have verdict PASS

  2. Construct certificate:
     cert = {
       "version": "1.0.0",
       "design_id": sha256(ir_seq.sequence),
       "sequence": ir_seq.sequence,
       "types": [
         {
           "predicate": result.predicate.name(),
           "verdict": "PASS",
           "derivation": [{"step": s, "evidence": e} for s, e in result.derivation]
         }
         for result in type_results
       ],
       "optimization": {
         "objective": optimization.objective.name(),
         "value": optimization.objective_value,
         "constraints_satisfied": [c.name() for c in optimization.constraints]
       },
       "composition": [
         {
           "check_type": check.type,
           "genes_involved": check.genes,
           "verdict": "PASS",
           "evidence": check.evidence
         }
         for check in (composition_results or [])
       ],
       "provenance": {
         "tool": "BioCompiler",
         "version": "1.0.0",
         "timestamp": iso8601_now(),
         "parameters": input_parameters,
         "input_hash": sha256(original_input)
       }
     }

  3. Return cert
```

#### 3.7.2 Algorithm: Standalone Certificate Verification

```
ALGORITHM VerifyCertificate:
  INPUT:  GuaranteeCertificate JSON, system rule definitions
  OUTPUT: VERIFIED or REJECTED with specific failure

  1. Verify design_id: sha256(cert.sequence) == cert.design_id
  2. For each type entry in cert.types:
     Re-evaluate the predicate against cert.sequence
     Assert verdict == PASS
     Assert derivation trace is valid (each step follows from the previous)
  3. If composition checks present:
     For each composition check:
       Re-verify the cross-gene constraint
       Assert verdict == PASS
  4. Verify provenance: all required fields present
  5. Return VERIFIED if all checks pass; REJECTED otherwise
```

---

### 3.8 COMP-08: Compositional Verifier

#### 3.8.1 Algorithm: Four Composition Checks

```
ALGORITHM VerifyComposition:
  INPUT:  IRCircuit, gene_results dict
  OUTPUT: list of CompositionChecks

  checks = []

  # Check 1: Promoter Conflict
  For each pair (gene_i, gene_j) where i != j:
    tfs_produced_i = gene_i.promoter.tf_activators + gene_i.promoter.tf_repressors
    tfs_regulating_j = gene_j.promoter.tf_activators + gene_j.promoter.tf_repressors
    overlap = tfs_produced_i ∩ tfs_regulating_j
    If overlap is non-empty:
      checks.append(CompositionCheck("promoter_conflict", [i, j],
                     FAIL if unintended else UNCERTAIN, evidence))

  # Check 2: Resource Competition
  total_ribosome_demand = sum(estimate_ribosome_demand(g) for g in circuit.genes)
  capacity = RIBOSOME_CAPACITY[circuit.organism][circuit.cellular_context]
  If total_ribosome_demand > capacity:
    checks.append(CompositionCheck("resource_competition", ALL_GENES,
                   FAIL, f"Demand {total} exceeds capacity {capacity}"))

  # Check 3: Splicing Interference
  For each pair (gene_i, gene_j) where i != j:
    cryptic_sites_i = gene_i.type_results["NoCrypticSplice"].violations
    functional_sites_j = gene_j.splicing_sites
    For each cryptic in cryptic_sites_i:
      If cryptic.overlaps(functional_sites_j):
        checks.append(CompositionCheck("splicing_interference", [i, j],
                       FAIL, evidence))

  # Check 4: RNA-RNA Interaction
  For each pair of transcripts (t_i, t_j):
    complementary_regions = find_complementary(t_i, t_j, min_length=15)
    If complementary_regions:
      checks.append(CompositionCheck("rna_interaction", [i, j],
                     UNCERTAIN, f"Complementary region of length {len}"))

  Return checks
```

---

### 3.9 COMP-09: Mutation Explorer

#### 3.9.1 Algorithm: Grammar-Guided Mutation Decomposition

```
ALGORITHM ExploreMutations:
  INPUT:  IR-Seq, max_mutations, mutation_categories
  OUTPUT: MutationReport

  1. Decompose gene into grammar nonterminals:
     nonterminals = []
     For each exon in IR-Seq.exon_boundaries:
       nonterminals.append(("intra_exonic", exon.start, exon.end))
     For each splice_site in IR-Seq.splice_sites:
       nonterminals.append(("splice_site", site.position, site.position+2))
     For each regulatory_element in IR-Seq.regulatory_elements:
       nonterminals.append(("regulatory", elem.position, elem.position+len(elem)))

  2. Enumerate mutations per nonterminal:
     For each nt in nonterminals:
       If nt.type == "intra_exonic":
         nt.mutations = synonymous_substitutions(nt.region)
       If nt.type == "splice_site":
         nt.mutations = point_mutations(nt.region)
       If nt.type == "regulatory":
         nt.mutations = point_mutations(nt.region)

  3. Exploit independence:
     independent_groups = partition_by_nonoverlap(nonterminals)
     For each group:
       group.combinations = product(member.mutations for member in group)

  4. Filter for grammar validity:
     legal_combinations = []
     For each combination in cartesian_product(group.combinations):
       mutated_seq = apply_mutations(IR-Seq.sequence, combination)
       isoform_set = NDFST_splice(mutated_seq)
       If isoform_set is non-empty:
         legal_combinations.append(combination)

  5. Detect conflicts:
     conflicts = find_conflicts(legal_combinations, constraints)

  6. Return MutationReport(nonterminals, legal_combinations, conflicts)
```

---

### 3.10 COMP-10: ORF Analyzer

#### 3.10.1 Algorithm: Overlapping Reading Frame Analysis

```
ALGORITHM AnalyzeOverlappingORFs:
  INPUT:  sequence, list of ReadingFrameSpecs
  OUTPUT: ORFAnalysisReport

  1. Construct translation FST per frame:
     For each frame_spec in reading_frames:
       frame_fsts[frame_spec.name] = build_fst(sequence, frame_spec)

  2. Compute per-position frame membership:
     For position p = 0 to len(sequence) - 1:
       affecting_frames[p] = [f for f in reading_frames
                               if f.start <= p < f.end
                               and (p - f.start) % 3 == f.frame]

  3. Compute shared constraint set:
     shared = {p for p, frames in affecting_frames.items() if len(frames) >= 2}

  4. Coupling classification:
     For each position p:
       If len(affecting_frames[p]) >= 2:
         coupling[p] = "high"
       Else:
         coupling[p] = "low"

  5. Detect constraint conflicts:
     conflicts = []
     For each position p in shared:
       optimal_codons_per_frame = {}
       For each frame f in affecting_frames[p]:
         codon = get_codon_at(sequence, p, f)
         optimal = optimal_codon_for_position(p, f, objective)
         optimal_codons_per_frame[f] = optimal
       If not all_same_codon(optimal_codons_per_frame.values()):
         conflicts.append(Conflict(p, optimal_codons_per_frame))

  6. Return ORFAnalysisReport(shared, coupling, conflicts)
```

---

## 4. Cross-Cutting Design Concerns

### 4.1 Determinism Guarantee

The system achieves byte-identical output for identical input through the following mechanisms:

1. **No randomness**: No use of `random`, `numpy.random`, or any non-deterministic library in the core pipeline.
2. **No floating-point ambiguity**: All comparisons use exact integer arithmetic where possible. Floating-point comparisons (e.g., GC content) use a tolerance of ±0.01.
3. **No hash randomization**: Python's `PYTHONHASHSEED=0` is set for deterministic dictionary ordering.
4. **No parallel non-determinism**: All internal pipeline stages are pure functions. FFI stages are explicitly marked as non-deterministic.
5. **No external state**: No file I/O, no network access, no environment variable reads during core pipeline execution.

### 4.2 Error Handling Strategy

```
BioCompilerError (base)
  ├── InputError
  │     ├── InvalidSequenceError
  │     ├── FormatError
  │     └── EmptySequenceError
  ├── PipelineError
  │     ├── NoValidIsoformError
  │     ├── FrameError
  │     ├── NoStartCodonError
  │     └── UnknownCellTypeError
  ├── FFIError
  │     ├── AdapterNotFoundError
  │     ├── ExternalToolError
  │     ├── OutputParseError
  │     └── OutputValidationError
  ├── OptimizationError
  │     ├── InfeasibleError (contains MUS)
  │     └── SolverTimeoutError
  └── CertificateError
        └── VerificationError
```

### 4.3 Logging and Observability

- **Structured logging**: JSON-formatted log entries with component ID, pipeline run ID, timestamp, and message.
- **Provenance tracking**: Every IR record carries a Provenance object recording the tool, version, timestamp, and parameters that produced it.
- **Pipeline run ID**: Each pipeline invocation gets a unique UUID. All log entries and IR records for a run share this ID.

### 4.4 Configuration Management

| Config File | Content | Format |
|---|---|---|
| `config/splicing_rules.yaml` | Splicing grammar rules, consensus sequences, thresholds | YAML |
| `config/codon_tables/*.csv` | Per-organism codon usage tables | CSV |
| `config/cell_contexts/*.yaml` | Per-cell-type splicing factor concentrations, thresholds | YAML |
| `config/pwm/*.npy` | Position weight matrices for splice sites | NumPy arrays |

Configuration files are loaded at pipeline initialization time and are immutable during a pipeline run. Changes to configuration files require restarting the pipeline.

---

## 5. Data Design

### 5.1 IR Data Flow

```
Input FASTA
    │
    v
[COMP-01: Scanner] → IR-Seq (with tokens)
    │
    v
[COMP-02: Splicing Engine] → IR-Seq (per isoform, with exon boundaries)
    │
    v
[COMP-03: Translation Engine] → IR-Peptide (with codon assignments)
    │
    ├──→ [COMP-04: Folding FFI] → IR-Structure
    ├──→ [COMP-04: PTM FFI] → IR-Peptide (with PTM SLOTs filled)
    │
    v
[COMP-05: Type System] → Verification Report
    │
    v (if not all PASS → loop back through COMP-06)
[COMP-06: Optimizer] → Optimized IR-Seq
    │
    v
[COMP-07: Certificate Generator] → Guarantee Certificate (JSON)
```

### 5.2 Persistence Strategy

- **In-memory processing**: All IR objects live in memory during pipeline execution. No intermediate disk writes for core pipeline stages.
- **Output serialization**: Final outputs (optimized IR-Seq, guarantee certificates) are serialized to disk using protocol buffers (IR) and JSON (certificates).
- **Checkpointing (future)**: For long-running FFI stages, the FFI Manager writes intermediate results to disk so that pipeline runs can be resumed after FFI failures.

---

## 6. Design Decisions Log

| ID | Decision | Rationale | Alternatives Considered | Date |
|---|---|---|---|---|
| DD-01 | AC-3 for constraint propagation | Well-understood, efficient for CSPs with small domains (1–6 per variable) | Forward checking only; full look-ahead | 2026-06-07 |
| DD-02 | Backtracking with CAI-ordered domain values | Prioritizes high-CAI codons, finding near-optimal solutions faster | Random ordering; least-constraining value | 2026-06-07 |
| DD-03 | Deletion-based MUS computation | Simple, correct, sufficient for constraint sets of size ≤ 10 | QuickXPlain; SAT-based MUS extraction | 2026-06-07 |
| DD-04 | Complementary region detection via simple string matching | O(n²) but n ≤ 10,000 bp per transcript; fast enough | Smith-Waterman; RNAhybrid | 2026-06-07 |
| DD-05 | JSON for certificates (not protobuf) | Human-readable; easily inspectable by regulators | Protocol buffers; XML | 2026-06-07 |
| DD-06 | Protocol buffers for IR (not JSON) | Schema enforcement, efficient serialization, code generation | JSON Schema; HDF5; pickle | 2026-06-07 |
| DD-07 | SHA-256 for certificate design_id | Cryptographic integrity; universally available | MD5; BLAKE2; no hash | 2026-06-07 |
