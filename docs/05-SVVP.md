# DOC-05: Software Verification & Validation Plan (SVVP)

> **WARNING:** **CORRECTIONS (v1.0.0):** (1) "IR-Seq [NOT IMPLEMENTED — use Python dataclasses]" and other IR class references — NOT IMPLEMENTED; Python dataclasses (Token, SpliceIsoform, TypeCheckResult, Certificate) are used instead. (2) "five-valued AND truth table" → "five-valued AND truth table"; the system uses five-valued verdicts (PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL). The Lean4 formal model uses a three-valued subset (PASS/FAIL/UNCERTAIN) that the five-valued logic refines to. (3) 43 type predicates (Lean4 model: 17 core + 19 SLOT = 36; Python implementation adds 7 extended diagnostic predicates for 43 total). (4) Test cases COMP-08, COMP-09, COMP-10 are [NOT IMPLEMENTED]. (5) "# protobuf ≥ 4.24 [REMOVED — not used]" is NOT a dependency — Protocol Buffers were not implemented. (6) Non-existent exception references (MissingPrerequisiteError, AdapterNotFoundError, ExternalToolError, OutputParseError, OutputValidationError) should be disregarded.**


| Field | Value |
|---|---|
| **Document ID** | DOC-05 |
| **Version** | 0.9.0 |
| **Status** | Current |
| **Date** | 2026-06-07 |
| **Prepared By** | BioCompiler Project Team |
| **Reviewed By** | [TBD at baseline review] |
| **Approved By** | [TBD at baseline review] |

---

## 1. Introduction

### 1.1 Purpose

This document constitutes the Software Verification & Validation Plan (SVVP) for the **BioCompiler** system, prepared in accordance with IEEE 1012-2016 (Standard for System, Software, and Hardware Verification and Validation). Its purpose is to define the complete set of verification and validation activities that SHALL be performed to confirm that the BioCompiler system satisfies every requirement specified in the Software Requirements Specification (DOC-01) and that the final system is fit for its intended use by synthetic biologists, circuit designers, virologists, bioinformaticians, and regulatory reviewers.

Verification confirms that each development artifact (design, code, test) correctly implements the requirements of the preceding artifact ("are we building the system right?"). Validation confirms that the system, when placed in its intended operational environment, fulfills its intended purpose ("are we building the right system?").

### 1.2 Scope

This SVVP covers all verification and validation activities for the BioCompiler system from requirements analysis through final acceptance testing. It applies to all ten components (COMP-01 through COMP-10), the IR Bus, the CLI, the Python API, the FFI adapter framework, and the certificate generation and verification subsystem.

**In scope:**

- Unit testing of every component (Section 2.1)
- Integration testing of component pairs and the full pipeline (Section 2.2)
- System testing against use cases UC-01, UC-02, UC-03 (Section 2.3)
- Soundness testing — adversarial inputs with known violations (Section 2.4)
- Reproducibility testing leveraging BioCompiler's determinism (Section 2.5)
- Biological validation against reference databases (Section 3.1)
- Comparison against existing gene design tools (Section 3.2)
- Pilot validation studies (Section 3.3)
- Acceptance criteria (Section 4)
- Test environment specification (Section 5)
- Defect classification and SLAs (Section 6)
- Requirements-to-test-case traceability (Section 7)

**Out of scope:**

- Verification of external tools invoked via FFI (AlphaFold, NetPhos, etc.) — these are treated as black boxes per REQ-FUNC-052
- Wet-lab experimental validation of designed sequences
- Certification of the development process (ISO 9001, CMMI) — this is a verification plan for the product, not the process

### 1.3 V&V Strategy — Exploiting Determinism

BioCompiler's defining architectural property — **full determinism for all internal pipeline stages** (REQ-NFR-010) — enables three categories of testing that are impractical for non-deterministic systems:

1. **Exact Reproducibility Testing**: Because the same input always produces bit-identical output (excluding FFI stages), we can run the entire pipeline 100 times on the same input and assert byte-identical results. Any non-determinism in a core pipeline stage is an immediate failure. This provides a powerful regression guard: any change that alters output for the same input is flagged.

2. **Snapshot / Regression Testing**: Because output is deterministic, we can store expected outputs (golden files) for a suite of reference inputs and compare every run against those snapshots. Snapshot tests are exact — not approximate — so even a single-bit drift is detected. This enables continuous regression testing on every commit with zero tolerance for output variation.

3. **Adversarial / Soundness Testing**: Because the type system is sound (REQ-NFR-011: no false PASS), we can construct sequences with *known* biological violations (cryptic splice sites, restriction sites, out-of-range GC content, reading frame disruptions) and verify that the system *always* returns FAIL. A single false PASS on an adversarial input is a Critical defect. This is the most important test category: it directly validates the system's core safety property.

These three categories form the backbone of the V&V strategy. Reproducibility testing guards against accidental non-determinism; snapshot testing guards against regressions; and adversarial testing guards against unsoundness. Together, they provide coverage that is both broader and deeper than what is achievable for probabilistic systems.

### 1.4 References

| Ref ID | Document |
|---|---|
| REF-01 | DOC-01: Software Requirements Specification (SRS) |
| REF-02 | DOC-02: Software Architecture Document (SAD) |
| REF-03 | DOC-03: Software Design Document (SDD) |
| REF-04 | DOC-04: Interface Control Document (ICD) |
| REF-05 | IEEE 1012-2016: Standard for System, Software, and Hardware Verification and Validation |
| REF-06 | GENCODE v44, human genome annotation (GRCh38) |
| REF-07 | UniProt Knowledgebase |
| REF-08 | IDT Codon Optimization Tool |
| REF-09 | Kazusa Codon Usage Database |

---

## 2. Verification Activities

### 2.1 Unit Testing

**Framework**: pytest with pytest-cov for coverage measurement, Hypothesis for property-based testing.

**Coverage targets**:
- Overall: ≥ 90% line coverage
- COMP-05 (Type System): 100% line coverage
- COMP-06 (Optimizer): 100% line coverage

**CI execution**: Every commit triggers the full unit test suite via `make test`. A failing unit test blocks the pull request.

**Property-based testing**: Hypothesis strategies generate random nucleotide sequences, random codon assignments, and random constraint sets for DFA/FST engines and the CSP solver. These complement the explicit test cases below.

#### Complete Unit Test Case Catalog

| Test Case ID | Requirement | Component | Test Description | Pass Criterion |
|---|---|---|---|---|
| TC-U-001 | REQ-FUNC-001 | COMP-01 (Scanner) | Accept valid FASTA-formatted DNA input with header line and nucleotide sequence | Scanner produces valid IR-Seq [NOT IMPLEMENTED — use Python dataclasses] with correct sequence length and no errors |
| TC-U-002 | REQ-FUNC-001 | COMP-01 (Scanner) | Accept raw nucleotide string (no FASTA header) as input | Scanner produces valid IR-Seq [NOT IMPLEMENTED — use Python dataclasses]; T→U normalization applied in RNA mode |
| TC-U-003 | REQ-FUNC-001 | COMP-01 (Scanner) | Reject input containing non-IUPAC characters (e.g., 'Z', '0', whitespace within sequence) | `InvalidSequenceError` raised with correct position and offending character in message |
| TC-U-004 | REQ-FUNC-002 | COMP-01 (Scanner) | Detect and annotate start codons (AUG) in the input sequence | Token list contains start_codon entries at correct positions with match_sequence "AUG" |
| TC-U-005 | REQ-FUNC-002 | COMP-01 (Scanner) | Detect splice donor sites (GT dinucleotide with surrounding consensus) in the input sequence | Token list contains splice_donor entries at correct positions with score ≥ threshold |
| TC-U-006 | REQ-FUNC-002 | COMP-01 (Scanner) | Detect restriction enzyme recognition sites from a user-specified enzyme set (EcoRI: GAATTC, BamHI: GGATCC, XhoI: CTCGAG) | Token list contains restriction_site entries at all positions where each enzyme pattern occurs |
| TC-U-007 | REQ-FUNC-004 | COMP-01 (Scanner) | Determinism: scan the same 10 kb input sequence 100 times and compare all outputs | All 100 token lists are byte-identical (same positions, same scores, same ordering) |
| TC-U-010 | REQ-FUNC-010 | COMP-02 (Splicing) | Constitutive exon gene (no alternative splicing sites): NDFST produces exactly 1 isoform | Isoform set cardinality = 1; isoform matches expected spliced sequence |
| TC-U-011 | REQ-FUNC-012 | COMP-02 (Splicing) | Gene with one alternatively spliced exon (cassette exon): NDFST produces ≥ 2 isoforms | Isoform set cardinality ≥ 2; one isoform includes the cassette exon, one excludes it |
| TC-U-012 | REQ-FUNC-011 | COMP-02 (Splicing) | Gene with canonical GT-AG splice sites: verify donor/acceptor pair annotations in each isoform | Every isoform's exon boundaries align with GT donor / AG acceptor positions from the scanner token stream |
| TC-U-013 | REQ-FUNC-011 | COMP-02 (Splicing) | Gene with rare GC-AG non-canonical splice site: verify NDFST includes a parse path using GC-AG | Isoform set includes at least one isoform utilizing the GC-AG splice site; GC-AG isoform has lower confidence score than GT-AG isoform |
| TC-U-014 | REQ-FUNC-013 | COMP-02 (Splicing) | Same gene scanned under two different cellular contexts (HEK293T vs. HepG2): verify different isoform sets | Isoform set cardinalities differ OR at least one isoform present in one context is absent in the other; both sets are non-empty |
| TC-U-015 | REQ-FUNC-014 | COMP-02 (Splicing) | Determinism: run NDFST on the same gene + context 100 times and compare all isoform sets | All 100 isoform sets are identical (same cardinality, same sequences, same exon boundaries, same ordering) |
| TC-U-020 | REQ-FUNC-020 | COMP-03 (Translation) | Translate a standard coding sequence using the standard genetic code | Amino acid sequence matches expected protein; each codon maps to the correct amino acid per the standard table |
| TC-U-021 | REQ-FUNC-021 | COMP-03 (Translation) | Translate an mRNA with a UGA codon and a SECIS element in the 3' UTR: selenocysteine insertion | UGA codon is flagged as selenocysteine (amino acid "U") with SelenocysteineFlag; not treated as stop |
| TC-U-022 | REQ-FUNC-022 | COMP-03 (Translation) | Determinism: translate the same spliced mRNA 100 times | All 100 amino acid sequences are identical; all codon assignments, SEC flags, and frameshift warnings are identical |
| TC-U-023 | REQ-FUNC-021 | COMP-03 (Translation) | Translate an mRNA containing a known programmed ribosomal frameshift motif (e.g., HIV-1 slippery sequence U UUU UUUA) | FrameshiftWarning emitted at the motif position with correct direction (-1) and motif sequence; translation continues in the original frame (FST does not resolve frameshift) |
| TC-U-030 | REQ-FUNC-030, REQ-FUNC-031 | COMP-05 (Type System) | Type-check an mRNA that satisfies all declared predicates (correct splicing, no cryptic sites, CAI ≥ 0.8, GC ∈ [40%,60%], no restriction sites, in frame, no instability motifs) | All verdicts = PASS; each result includes a derivation trace |
| TC-U-031 | REQ-FUNC-031 | COMP-05 (Type System) | Type-check an mRNA with a known cryptic splice site (synthetic GT-AG-like motif inside an exon) | NoCrypticSplice verdict = FAIL; violation identifies the position and score of the cryptic site |
| TC-U-032 | REQ-FUNC-032 | COMP-05 (Type System) | Determinism: type-check the same mRNA with the same predicates and context 100 times | All 100 verdict lists are identical (same verdicts, same derivation traces, same violation positions) |
| TC-U-033 | REQ-FUNC-033 | COMP-05 (Type System) | Verify subtyping: if SpliceCorrect(HEK293T) passes, then SpliceCorrect(liver) passes where HEK293T is a stricter context | SpliceCorrect(liver) returns PASS when SpliceCorrect(HEK293T) returns PASS and HEK293T is a subtype of liver context |
| TC-U-034 | REQ-FUNC-034 | COMP-05 (Type System) | Conjunction of PASS and UNCERTAIN: type-check with one predicate PASS and another UNCERTAIN | Overall verdict = UNCERTAIN (per five-valued AND truth table) |
| TC-U-035 | REQ-FUNC-034 | COMP-05 (Type System) | Conjunction of FAIL with any other verdict: type-check with one predicate FAIL and another PASS or UNCERTAIN | Overall verdict = FAIL (per five-valued AND truth table: FAIL dominates) |
| TC-U-040 | REQ-FUNC-040, REQ-FUNC-043 | COMP-06 (Optimizer) | Feasible optimization problem: protein of 200 aa, constraints CAI ≥ 0.8, GC ∈ [40%,60%], no restriction sites (3 enzymes) | Solver returns FeasibleResult with a complete codon assignment; CAI ≥ 0.8 verified; all hard constraints satisfied |
| TC-U-041 | REQ-FUNC-044 | COMP-06 (Optimizer) | Infeasible optimization problem: impossible constraints (CAI ≥ 0.99 AND GC ∈ [70%,80%] for the same protein) | Solver returns InfeasibleResult with a non-empty MUS; MUS is verified as genuinely unsatisfiable |
| TC-U-042 | REQ-FUNC-042 | COMP-06 (Optimizer) | Optimization with a splicing correctness constraint: ensure optimized sequence has no cryptic splice sites | Optimized codon assignment produces a sequence where NoCrypticSplice = PASS; verify by re-running type checker |
| TC-U-043 | REQ-FUNC-040 | COMP-06 (Optimizer) | Completeness: construct a problem with a known feasible solution, verify solver finds one | Solver returns FeasibleResult; returned assignment is verified to satisfy all constraints (independent re-check) |
| TC-U-044 | REQ-FUNC-044 | COMP-06 (Optimizer) | MUS correctness: for an infeasible problem, verify the MUS is minimal (removing any single constraint makes it satisfiable) | For each constraint c in MUS, MUS \ {c} is satisfiable (verified by independent solve) |
| TC-U-045 | REQ-FUNC-045 | COMP-06 (Optimizer) | Determinism: solve the same CSP instance 100 times with identical parameters | All 100 results are identical: same assignment (feasible) or same MUS (infeasible), same objective value, same ordering |
| TC-U-050 | REQ-FUNC-052 | COMP-04 (FFI) | Request an FFI adapter that is not registered: adapter name does not exist in the registry | `AdapterNotFoundError` raised with the adapter name in the error message |
| TC-U-051 | REQ-FUNC-050, REQ-FUNC-051 | COMP-04 (FFI) | Invoke an external tool that returns a non-zero exit code (simulated via mock) | `ExternalToolError` raised with exit code and stderr content; SLOT fields remain empty; provenance records the failure |
| TC-U-052 | REQ-FUNC-052 | COMP-04 (FFI) | FFI output validation: external tool returns output that violates IR invariants (e.g., pLDDT score of 150, outside [0,100] range) | `OutputValidationError` raised; invalid output is rejected; SLOT fields remain empty |
| TC-U-053 | REQ-FUNC-052 | COMP-04 (FFI) | Provenance tracking: invoke a mock folding adapter and verify provenance metadata is recorded | Provenance includes tool name, version, timestamp, input_hash (SHA-256 of amino acid sequence), and adapter parameters |
| TC-U-060 | REQ-FUNC-060 | COMP-07 (Certificate) | Generate a guarantee certificate for a design that passes all type checks with PASS verdicts | Certificate JSON contains: version, design_id (SHA-256), sequence, all type entries with PASS verdicts and derivation traces, optimization record, provenance metadata |
| TC-U-061 | REQ-FUNC-061 | COMP-07 (Certificate) | Standalone verification of a valid certificate using the independent checker | Independent checker returns CERTIFICATE_VALID; re-derived verdicts match certificate claims; SHA-256 hash matches design_id |
| TC-U-062 | REQ-FUNC-062 | COMP-07 (Certificate) | Generate a circuit-level certificate for a multi-gene circuit that passes all composition checks | Circuit certificate contains individual_gene_certificates array, composition_checks array with PASS verdicts and evidence |
| TC-U-070 | REQ-FUNC-071 | COMP-08 (Compositional) | Promoter conflict detection: two-gene circuit where Gene A produces a TF that binds Gene B's promoter | Composition check of type "promoter_conflict" returns FAIL with evidence identifying the overlapping TF |
| TC-U-071 | REQ-FUNC-071 | COMP-08 (Compositional) | Resource competition: circuit with 10 genes where total ribosome demand exceeds estimated capacity | Composition check of type "resource_competition" returns FAIL with total demand and capacity values |
| TC-U-072 | REQ-FUNC-071 | COMP-08 (Compositional) | Splicing interference: Gene A's transcript contains a cryptic splice site that overlaps a functional site in Gene B | Composition check of type "splicing_interference" returns FAIL with the overlapping site positions |
| TC-U-073 | REQ-FUNC-071 | COMP-08 (Compositional) | RNA interaction: two transcripts with complementary regions ≥ 15 nt detected | Composition check of type "rna_interaction" returns UNCERTAIN with complementary region lengths and positions |
| TC-U-080 | REQ-FUNC-080 | COMP-09 (Mutation) | Decompose a 3-exon gene into mutation categories: intra-exonic, splice site, regulatory element | Mutation report contains three non-empty categories with correct position ranges matching splicing grammar nonterminals |
| TC-U-081 | REQ-FUNC-081 | COMP-09 (Mutation) | Enumerate legal single-mutation combinations for a gene with known mutation space | All enumerated combinations produce at least one valid isoform when re-run through the NDFST; no combination produces zero isoforms |
| TC-U-082 | REQ-FUNC-082 | COMP-09 (Mutation) | Independence: two mutations in different non-overlapping exons — verify the combination count equals the product of per-exon counts | Total combinations for both exons = (count_exon1 × count_exon2); factorization holds |
| TC-U-083 | REQ-FUNC-083 | COMP-09 (Mutation) | Conflict detection: two individually legal mutations that jointly violate the CodonAdapted constraint (each reduces CAI slightly; together they drop CAI below threshold) | Conflict report identifies the two mutations, the jointly violated constraint (CodonAdapted), and which individual constraints remain satisfied |
| TC-U-090 | REQ-FUNC-090 | COMP-10 (ORF) | Multi-frame construction: sequence with two overlapping reading frames (frame 0 and frame +1) — construct translation FST per frame | Two IR-Peptide [NOT IMPLEMENTED] records produced, one per frame; each has correct amino acid sequence for its frame offset |
| TC-U-091 | REQ-FUNC-091 | COMP-10 (ORF) | Shared constraint set computation: two overlapping frames on a 300 nt sequence — identify positions where mutations affect both frames | Shared constraint set is non-empty; each position in the set affects ≥ 2 frames; positions affecting only one frame are excluded |
| TC-U-092 | REQ-FUNC-092 | COMP-10 (ORF) | Coupling classification: classify all positions in a dual-frame region as high-coupling or low-coupling | Positions in shared constraint set classified as high-coupling; positions in only one frame classified as low-coupling; all positions classified |
| TC-U-093 | REQ-FUNC-093 | COMP-10 (ORF) | Conflict detection: two frames where optimal codon assignment for Frame 1 requires nucleotide X at position p, but Frame 2 requires nucleotide Y ≠ X at the same position | Conflict report identifies position p, both frames, the conflicting amino acid requirements, and the minimal conflict set |

---

### 2.2 Integration Testing

Integration tests verify that component pairs and multi-component sequences interact correctly through the IR Bus. Each test exercises at least two components in sequence, using real (not mocked) inter-component communication.

| Test Case ID | Components Integrated | Test Description | Pass Criterion |
|---|---|---|---|
| TC-I-001 | COMP-01 → COMP-02 (Scanner → Splicing) | Scan a 5-exon pre-mRNA gene, pass token stream to Splicing Engine, verify NDFST receives correct token stream and produces valid isoforms | Splicing Engine receives IR-Seq [NOT IMPLEMENTED — use Python dataclasses] with all tokens (start codons, splice sites, etc.); produces ≥ 1 isoform with correct exon boundaries |
| TC-I-002 | COMP-02 → COMP-03 (Splicing → Translation) | Splice a gene with 2 isoforms, translate each isoform, verify each produces a distinct IR-Peptide [NOT IMPLEMENTED] with correct amino acid sequences | Two IR-Peptide [NOT IMPLEMENTED] records produced; each has a different amino acid sequence; each references its source isoform ID |
| TC-I-003 | COMP-03 → COMP-04 (Translation → FFI) | Translate an mRNA, invoke a mock folding adapter, verify IR-Structure [NOT IMPLEMENTED] is produced with SLOT fields filled | IR-Structure [NOT IMPLEMENTED] contains C-alpha coordinates and pLDDT scores; provenance metadata recorded; INV-STR-01 and INV-STR-02 invariants hold |
| TC-I-004 | COMP-03 → COMP-05 (Translation → Type System) | Translate an mRNA, then type-check the resulting IR-Peptide [NOT IMPLEMENTED] against all 43 type predicates (17 core + 19 SLOT-dependent + 7 extended diagnostic) | Type checker receives IR-Peptide [NOT IMPLEMENTED] with valid codon assignments; produces verdicts for all 43 predicates; determinism holds |
| TC-I-005 | COMP-05 → COMP-06 (Type System → Optimizer) | Type-check a sequence that produces at least one FAIL verdict (e.g., NoCrypticSplice = FAIL), pass to Optimizer | Optimizer receives FAIL verdicts and constraint definitions; searches for feasible codon assignment that resolves the FAIL |
| TC-I-006 | COMP-06 → COMP-07 (Optimizer → Certificate) | Optimize a feasible gene design, then generate a guarantee certificate | Certificate includes optimized sequence, all PASS verdicts with derivation traces, CSP constraint set and assignment, and provenance |
| TC-I-007 | Full Pipeline (COMP-01 through COMP-07) | Execute the complete single-gene pipeline (UC-01 main success scenario): input protein + constraints → optimized mRNA + certificate | Pipeline produces optimized mRNA sequence; all type predicates return PASS; certificate is generated and passes standalone verification |
| TC-I-008 | Multi-Gene Pipeline (COMP-01×N → COMP-08) | Execute the multi-gene pipeline for a 3-gene circuit with no conflicts; each gene individually passes type checking; COMP-08 verifies composition | Circuit certificate generated; all individual gene certificates valid; all composition checks return PASS; independent checker confirms CERTIFICATE_VALID |
| TC-I-009 | Multi-Frame Pipeline (COMP-01 → COMP-03 → COMP-10) | Input a viral genome segment with 2 overlapping reading frames; execute scanning, translation per frame, and ORF analysis | ORF analysis report contains shared constraint set, coupling classification, and no constraint conflicts (for a compatible dual-frame input) |

---

### 2.3 System Testing

System tests exercise the complete BioCompiler system end-to-end against the use cases defined in the SRS (DOC-01, Section 2.6). Each test uses realistic input data and verifies the full output, including exit codes, output format, and certificate validity.

| Test Case ID | Use Case | Test Description | Pass Criterion |
|---|---|---|---|
| TC-S-001 | UC-01 | Design a splicing-safe gene for mammalian expression: input GFP protein (238 aa), target organism *Homo sapiens*, cell type HEK293T, constraints CAI ≥ 0.8, GC ∈ [40%, 60%], avoid {EcoRI, BamHI, XhoI, HindIII, NotI} | Optimized mRNA produced; all type predicates PASS; certificate generated and passes standalone verification; CLI exit code 0; output includes human-readable verification report and JSON certificate |
| TC-S-002 | UC-02 | Verify a multi-gene toggle switch circuit (2 genes: LacI repressor + TetR repressor with mutual inhibition promoters): verify promoter conflict detection, resource competition, splicing interference, and RNA interaction | Composition checks detect the intended mutual inhibition (PASS, since it is intentional) and verify no unintended conflicts; circuit certificate generated; CLI exit code 0 |
| TC-S-003 | UC-03 | Analyze overlapping reading frames in SARS-CoV-2 ORF1a/ORF1b region (programmed -1 frameshift): input nucleotide sequence with two frame specs (frame 0: ORF1a, frame -1: ORF1b) | ORF analysis report produced; shared constraint set computed; high-coupling positions identified near the frameshift site; frameshift warning emitted at the slippery sequence; CLI exit code 0 |
| TC-S-004 | N/A (infeasible case) | Design a gene with impossible constraints: CAI ≥ 0.95 AND GC ∈ [70%, 80%] AND no restriction sites (5 enzymes) for a 500 aa protein | Pipeline reports INFEASIBLE with MUS; MUS identifies the minimal conflicting constraint set; CLI exit code 1; no certificate generated; failure report includes human-readable MUS descriptions |

---

### 2.4 Soundness Testing (Critical)

Soundness is the most critical property of the BioCompiler system: the type system MUST NOT produce a PASS verdict for a sequence that violates a declared constraint (REQ-NFR-011). Soundness tests are adversarial: they construct sequences with *known* biological violations and verify that the type system *always* returns FAIL. A single false PASS on any soundness test is a **Critical** defect (see Section 6) that blocks release.

| Test Case ID | Predicate Under Test | Adversarial Input Description | Known Violation | Pass Criterion |
|---|---|---|---|---|
| TC-SND-001 | NoCrypticSplice | mRNA sequence with a synthetic GT-AG dinucleotide pair embedded inside an exon, surrounded by consensus-like flanking bases, creating a strong cryptic splice donor | Cryptic splice donor site with MaxEntScan-like score above the threshold | NoCrypticSplice verdict = FAIL; violation identifies the cryptic donor position and score |
| TC-SND-002 | SpliceCorrect(HEK293T) | Pre-mRNA sequence with an alternative exon that has functional ESE motifs for HEK293T context, producing a second isoform | Multiple isoforms exist (not a singleton set) | SpliceCorrect(HEK293T) verdict = FAIL; violation identifies the unintended isoform and the alternative splice path |
| TC-SND-003 | CodonAdapted(*H. sapiens*, 0.8) | mRNA where 40% of codons are rare codons (low-frequency synonymous codons for *H. sapiens*), producing CAI ≈ 0.45 | CAI < 0.8 | CodonAdapted verdict = FAIL; violation identifies the CAI value and positions of rare codons |
| TC-SND-004 | GCInRange(40, 60) | mRNA with 72% GC content (synthetic sequence rich in G and C nucleotides) | GC% = 72, outside [40, 60] | GCInRange verdict = FAIL; violation identifies the GC percentage |
| TC-SND-005 | NoRestrictionSite({EcoRI}) | mRNA containing the EcoRI recognition site GAATTC at position 450 | EcoRI site present | NoRestrictionSite verdict = FAIL; violation identifies the site position and match sequence |
| TC-SND-006 | InFrame | mRNA where a premature stop codon (UAA) appears at codon position 50 out of 200 | Premature stop codon disrupts reading frame | InFrame verdict = FAIL; violation identifies the premature stop codon position |
| TC-SND-007 | NoInstabilityMotif | mRNA containing three AUUUA instability motifs in the 3' UTR region | Instability motifs present | NoInstabilityMotif verdict = FAIL; violation identifies all three motif positions |
| TC-SND-008 | Combined | mRNA with two simultaneous violations: a cryptic splice site AND an EcoRI restriction site | Two distinct constraint violations | Both NoCrypticSplice and NoRestrictionSite verdicts = FAIL; overall verdict = FAIL (conjunction of FAIL with anything = FAIL); certificate is NOT generated |

---

### 2.5 Reproducibility Testing

BioCompiler's determinism guarantee (REQ-NFR-010) requires byte-identical output for identical input across runs and platforms. Reproducibility tests verify this property directly.

| Test Case ID | Scope | Test Description | Pass Criterion |
|---|---|---|---|
| TC-REP-001 | Single-platform, 100 runs | Run the full single-gene pipeline (COMP-01 through COMP-07) on the same input 100 times on the same machine (Linux x86_64); compare all outputs | All 100 outputs are byte-identical: same optimized sequence, same verdicts, same certificate JSON (excluding timestamp field if present) |
| TC-REP-002 | Cross-platform | Run the full single-gene pipeline on the same input on 4 platform configurations: Linux x86_64, Linux ARM64, macOS x86_64, macOS ARM64 | All 4 outputs are identical within floating-point tolerance (±1 ULP for score values); identical verdicts, identical sequences, identical certificate structure |
| TC-REP-003 | CSP solver reproducibility | Run the CSP optimizer on a feasible problem with 200 codon positions 100 times; compare assignments and objective values | All 100 results are identical: same codon assignment, same objective value, same ordering of solutions |
| TC-REP-004 | Cross-run certificate verification | Generate a certificate in Run 1; independently verify it in Run 2 (separate process invocation) | Independent checker returns CERTIFICATE_VALID in Run 2; re-derived verdicts match the certificate from Run 1 exactly |

---

## 3. Validation Activities

Validation activities confirm that the BioCompiler system produces biologically meaningful and scientifically valid results, not merely that it satisfies its formal requirements. These activities compare BioCompiler's outputs against real biological data and existing tools.

### 3.1 Biological Validation

Biological validation compares BioCompiler's computational outputs against established reference databases. Each test verifies that BioCompiler's deterministic analyses agree with known biological facts.

| Test Case ID | Reference Source | Test Description | Pass Criterion |
|---|---|---|---|
| TC-V-001 | GENCODE v44 (GRCh38) | Process 50 human genes with annotated alternative splicing from GENCODE; compare BioCompiler's computed isoform sets against GENCODE's annotated isoforms | For each gene, BioCompiler's computed isoform set is a superset of the GENCODE-annotated isoforms (BioCompiler may find additional isoforms not in GENCODE, but must not miss any annotated isoform) |
| TC-V-002 | UniProt Knowledgebase | Translate 100 human mRNA sequences from GENCODE; compare BioCompiler's amino acid outputs against UniProt canonical protein sequences | ≥ 99% exact match: at least 99 of 100 translations exactly match the UniProt canonical sequence; mismatches must be explained by known biological factors (selenocysteine, known sequence variants) |
| TC-V-003 | IDT Codon Optimization Tool | For 20 target proteins, compare BioCompiler's optimized codon assignments (under identical constraints) against IDT's codon optimization output | BioCompiler's CAI values are within ±0.05 of IDT's reported CAI; all BioCompiler constraints verified as satisfied; BioCompiler may differ in specific codon choices but achieves comparable CAI |
| TC-V-004 | REBASE (restriction enzyme database) | Scan 30 sequences with 10 restriction enzymes each; compare BioCompiler's detected restriction sites against REBASE's known site positions | 100% recall: all sites listed in REBASE are detected by BioCompiler; 100% precision: all sites detected by BioCompiler correspond to valid REBASE entries |
| TC-V-005 | MaxEntScan reference scores | Score 200 known splice sites (100 genuine, 100 decoy) using BioCompiler's splice site scoring; compare against MaxEntScan reference scores | Pearson correlation ≥ 0.95 between BioCompiler scores and MaxEntScan scores on the same sites; BioCompiler's accept/reject classification matches MaxEntScan on ≥ 95% of sites |

### 3.2 Comparison Against Existing Tools

| Test Case ID | Comparison Target | Test Description | Pass Criterion |
|---|---|---|---|
| TC-C-001 | GeneDesign (web tool) | Design 10 genes for *E. coli* expression using both BioCompiler and GeneDesign with equivalent constraints; compare codon usage statistics, GC content, and CAI | BioCompiler achieves CAI within ±0.03 of GeneDesign; GC content matches within ±2%; all BioCompiler hard constraints satisfied (GeneDesign does not guarantee splicing correctness) |
| TC-C-002 | DNAWorks | Optimize 10 genes for mammalian expression using both BioCompiler and DNAWorks; compare restriction site avoidance success and codon adaptation | BioCompiler successfully avoids all specified restriction sites in ≥ 10/10 cases; DNAWorks success rate compared; BioCompiler provides formal certificate of avoidance |
| TC-C-003 | GenSmart™ Codon Optimization | Optimize 5 genes with complex constraint sets (CAI + GC + restriction sites + no cryptic splice sites); compare feasibility determination and solution quality | BioCompiler produces a feasible solution for all 5 cases where GenSmart reports success; BioCompiler additionally provides a guarantee certificate verifying all constraints including NoCrypticSplice (which GenSmart does not check) |

### 3.3 Pilot Validation Studies

Three pilot studies provide end-to-end validation of the BioCompiler system on realistic use cases, with results reviewed by domain experts.

#### Pilot 1: GFP Gene Design for Mammalian Expression

| Field | Value |
|---|---|
| **Objective** | Design an optimized GFP (Green Fluorescent Protein) gene for expression in HEK293T cells with full splicing correctness guarantees |
| **Target protein** | EGFP, 239 amino acids (UniProt: C5MKY7) |
| **Organism** | *Homo sapiens* |
| **Cell type** | HEK293T |
| **Constraints** | CAI ≥ 0.8, GC ∈ [40%, 60%], avoid {EcoRI, BamHI, XhoI, HindIII, NotI}, NoCrypticSplice, NoInstabilityMotif |
| **Success criteria** | (1) Optimized mRNA produced with all type predicates PASS; (2) Certificate passes standalone verification; (3) CAI ≥ 0.8 confirmed; (4) No cryptic splice sites detected by independent review using MaxEntScan; (5) Comparison against commercial gene synthesis company's default optimization shows comparable or better CAI |

#### Pilot 2: Toggle Switch Circuit Verification

| Field | Value |
|---|---|
| **Objective** | Verify a 2-gene toggle switch circuit (LacI + TetR mutual inhibition) for synthetic biology deployment in *E. coli* |
| **Circuit** | Gene A: LacI under TetR-repressible promoter (pTet); Gene B: TetR under LacI-repressible promoter (pLac) |
| **Organism** | *Escherichia coli* K-12 |
| **Success criteria** | (1) Both genes individually pass all type predicates; (2) Composition check detects intentional mutual inhibition (PASS with annotation); (3) No unintended promoter conflicts detected; (4) No splicing interference (prokaryotic context — no splicing, verified as PASS); (5) Resource competition analysis reports total ribosome demand within *E. coli* capacity; (6) Circuit certificate generated and passes standalone verification |

#### Pilot 3: SARS-CoV-2 ORF1a/ORF1b Overlapping Frame Analysis

| Field | Value |
|---|---|
| **Objective** | Analyze the overlapping reading frames in the SARS-CoV-2 ORF1a/ORF1b region, where a programmed -1 ribosomal frameshift produces ORF1b from the same nucleotide sequence |
| **Input** | SARS-CoV-2 reference genome (NC_045512.2) nucleotides 266–13,483 (ORF1a + ORF1b region) |
| **Frames** | Frame 0: ORF1a (nucleotides 266–13,483); Frame -1: ORF1b (nucleotides 13,468–21,555, accessed via frameshift) |
| **Success criteria** | (1) ORF1a translation matches UniProt P0DTC1 (replicase polyprotein 1a) with ≥ 99% identity; (2) ORF1b translation matches UniProt P0DTC1 (replicase polyprotein 1ab) downstream region with ≥ 99% identity; (3) Frameshift warning emitted at the correct slippery sequence position (nucleotide ~13,467); (4) Shared constraint set correctly identifies positions where the frameshift causes overlap; (5) High-coupling positions near the frameshift site are correctly classified; (6) Constraint conflicts (if any) between frames are identified with correct minimal conflict sets |

---

## 4. Acceptance Criteria

### 4.1 System-Level Acceptance

The BioCompiler system is accepted only when ALL of the following system-level acceptance criteria are met:

| Acceptance ID | Criterion | Measurement Method | Threshold |
|---|---|---|---|
| AC-01 | All unit tests pass | `make test` execution | 100% pass rate; zero failures |
| AC-02 | Unit test coverage ≥ 90% overall; 100% for COMP-05 and COMP-06 | `pytest --cov` report | ≥ 90% line coverage overall; 100% for COMP-05, COMP-06 |
| AC-03 | All integration tests pass | `make test-integration` execution | 100% pass rate; zero failures |
| AC-04 | All system tests pass (TC-S-001 through TC-S-004) | Manual or CI execution of system test suite | 100% pass rate; zero failures |
| AC-05 | All soundness tests pass (TC-SND-001 through TC-SND-008) — zero false PASS verdicts | Adversarial test suite execution | 100% pass rate; zero false PASS verdicts (any false PASS blocks release) |
| AC-06 | All reproducibility tests pass (TC-REP-001 through TC-REP-004) | Reproducibility test suite execution on all target platforms | 100% byte-identical output on same platform; identical within ±1 ULP across platforms |
| AC-07 | All biological validation tests pass (TC-V-001 through TC-V-005) | Biological validation suite execution | All thresholds met (see Section 3.1 pass criteria) |

### 4.2 Per-Phase Acceptance

> **Note**: Development phase numbers (Phase 1–4) refer to project milestones and are separate from the optimizer's internal processing steps.

The BioCompiler system is developed in four phases, each with its own acceptance gate:

#### Phase 1: Core Pipeline (COMP-01, COMP-02, COMP-03, IR Bus)

| Criterion | Test Cases | Pass Condition |
|---|---|---|
| Scanner accepts FASTA and raw input; rejects non-IUPAC | TC-U-001, TC-U-002, TC-U-003 | All pass |
| Scanner detects all motif types (start codons, splice sites, restriction sites) | TC-U-004, TC-U-005, TC-U-006 | All pass |
| Scanner is deterministic (100 runs) | TC-U-007 | Byte-identical output |
| Splicing Engine produces correct isoform sets | TC-U-010, TC-U-011, TC-U-012, TC-U-013 | All pass |
| Splicing Engine respects cellular context | TC-U-014 | Different isoform sets for different contexts |
| Splicing Engine is deterministic | TC-U-015 | Byte-identical output across 100 runs |
| Translation Engine handles standard and special cases | TC-U-020, TC-U-021, TC-U-023 | All pass |
| Translation Engine is deterministic | TC-U-022 | Byte-identical output across 100 runs |
| Scanner→Splicing integration | TC-I-001 | Pass |
| Splicing→Translation integration | TC-I-002 | Pass |
| Phase 1 unit test coverage ≥ 90% | All TC-U-0xx | ≥ 90% line coverage |

#### Phase 2: Verification & Optimization (COMP-05, COMP-06, COMP-07)

| Criterion | Test Cases | Pass Condition |
|---|---|---|
| Type System produces correct five-valued verdicts | TC-U-030, TC-U-031, TC-U-034, TC-U-035 | All pass |
| Type System implements subtyping correctly | TC-U-033 | Pass |
| Type System is deterministic | TC-U-032 | Byte-identical verdicts across 100 runs |
| Optimizer solves feasible problems | TC-U-040, TC-U-043 | FeasibleResult returned; constraints verified |
| Optimizer detects infeasible problems with correct MUS | TC-U-041, TC-U-044 | InfeasibleResult with verified minimal MUS |
| Optimizer respects splicing constraints | TC-U-042 | Optimized sequence passes NoCrypticSplice |
| Optimizer is deterministic | TC-U-045 | Identical results across 100 runs |
| Certificate generation produces valid certificates | TC-U-060 | Valid JSON with all required fields |
| Certificate standalone verification works | TC-U-061 | CERTIFICATE_VALID |
| Type System→Optimizer integration | TC-I-004, TC-I-005 | Pass |
| Optimizer→Certificate integration | TC-I-006 | Pass |
| All Phase 2 soundness tests pass | TC-SND-001 through TC-SND-008 | Zero false PASS |
| Phase 2 unit test coverage: 100% for COMP-05 and COMP-06 | All TC-U-03x, TC-U-04x | 100% line coverage |

#### Phase 3: Composition & Exploration (COMP-08, COMP-09, COMP-10)

| Criterion | Test Cases | Pass Condition |
|---|---|---|
| Compositional Verifier detects promoter conflicts | TC-U-070 | FAIL verdict with correct evidence |
| Compositional Verifier detects resource competition | TC-U-071 | FAIL verdict with demand/capacity values |
| Compositional Verifier detects splicing interference | TC-U-072 | FAIL verdict with overlapping site positions |
| Compositional Verifier detects RNA interactions | TC-U-073 | UNCERTAIN verdict with complementary regions |
| Mutation Explorer decomposes mutation space | TC-U-080 | Three non-empty categories |
| Mutation Explorer enumerates legal combinations | TC-U-081 | All combinations produce ≥ 1 valid isoform |
| Mutation Explorer exploits independence | TC-U-082 | Factorization holds |
| Mutation Explorer detects constraint conflicts | TC-U-083 | Conflict reported with jointly violated constraint |
| ORF Analyzer constructs multi-frame translations | TC-U-090 | Two correct IR-Peptide [NOT IMPLEMENTED] records |
| ORF Analyzer computes shared constraint set | TC-U-091 | Non-empty shared set with correct positions |
| ORF Analyzer classifies coupling | TC-U-092 | All positions classified |
| ORF Analyzer detects frame conflicts | TC-U-093 | Conflict with minimal conflict set |
| Circuit certificate generation | TC-U-062 | Valid circuit certificate |
| Multi-gene pipeline integration | TC-I-008 | Pass |
| Multi-frame pipeline integration | TC-I-009 | Pass |
| Phase 3 unit test coverage ≥ 90% | All TC-U-07x, TC-U-08x, TC-U-09x | ≥ 90% line coverage |

#### Phase 4: Full System Validation & Release

| Criterion | Test Cases | Pass Condition |
|---|---|---|
| All system tests pass | TC-S-001 through TC-S-004 | 100% pass |
| All soundness tests pass | TC-SND-001 through TC-SND-008 | Zero false PASS |
| All reproducibility tests pass | TC-REP-001 through TC-REP-004 | Byte-identical outputs |
| All biological validation tests pass | TC-V-001 through TC-V-005 | All thresholds met |
| All comparison tests pass | TC-C-001 through TC-C-003 | All criteria met |
| Pilot 1 (GFP) passes | GFP pilot | All 5 success criteria met |
| Pilot 2 (Toggle Switch) passes | Toggle switch pilot | All 6 success criteria met |
| Pilot 3 (SARS-CoV-2) passes | SARS-CoV-2 pilot | All 6 success criteria met |
| Full pipeline integration test passes | TC-I-007 | Pass |
| CLI exit codes correct for all test scenarios | TC-S-001 through TC-S-004 | Exit codes: 0 (success), 1 (infeasible), 2 (uncertain), 10 (input error) |
| No Critical or High severity open defects | Defect tracker | Zero open Critical/High defects |
| Documentation complete and consistent | All DOC-00 through DOC-10 | All documents at baseline status; no TBD items |

---

## 5. Test Environment

| Parameter | Specification |
|---|---|
| **Operating Systems** | Ubuntu 22.04 LTS (x86_64); Ubuntu 22.04 LTS (ARM64); macOS 14 Sonoma (x86_64); macOS 14 Sonoma (Apple Silicon / M2) |
| **Python** | 3.12.x (exact version pinned in pyproject.toml) |
| **CPU** | Intel Core i7-12700 (baseline for performance tests) or Apple M1/M2 (equivalent) |
| **RAM** | 32 GB (single-gene tests); 64 GB (circuit-level tests) |
| **GPU** | Optional: NVIDIA A100 or Apple M2 GPU (FFI folding tests only; not required for core pipeline tests) |
| **Test Data** | GENCODE v44 GFF3 annotation files (human GRCh38); Kazusa codon usage tables (human, *E. coli*, yeast); REBASE restriction enzyme database (curated subset of 50 common enzymes); 50 human gene fixtures (FASTA) with known splicing patterns; SARS-CoV-2 reference genome NC_045512.2; 200 synthetic test sequences with known ground-truth violations for soundness testing |
| **CI/CD** | GitHub Actions; on every pull request: `make proto && make lint && make typecheck && make test && make test-integration && make audit`; on merge to main: nightly full test suite including reproducibility tests and biological validation |
| **Dependency Versions** | Z3 solver ≥ 4.12; # protobuf ≥ 4.24 [REMOVED — not used]; pytest ≥ 8.0; hypothesis ≥ 6.90; ruff ≥ 0.4; mypy ≥ 1.10 (all pinned in pyproject.toml) |

---

## 6. Defect Classification

Every defect discovered during V&V activities is classified by severity and assigned a Service Level Agreement (SLA) for resolution. The classification follows the impact on BioCompiler's safety properties (soundness > determinism > functionality > usability).

| Severity | Definition | Examples | Resolution SLA | Release Blocking? |
|---|---|---|---|---|
| **Critical** | Violation of a core safety property: soundness (false PASS), determinism failure in a core pipeline stage, or data corruption | Type system returns PASS for a sequence with a known cryptic splice site; non-deterministic output from COMP-05 on identical input; certificate claims PASS but independent checker finds FAIL | 24 hours for triage; fix within 72 hours; mandatory re-execution of all soundness and reproducibility tests before merge | **Yes** — no release with an open Critical defect |
| **High** | Incorrect output that does not violate soundness; wrong FAIL verdict (false alarm); incorrect MUS; incorrect isoform set (missing a valid isoform); performance regression exceeding 2× the target | Type system returns FAIL for a sequence that should pass; optimizer reports INFEASIBLE when a feasible solution exists; NDFST misses an annotated GENCODE isoform; scanning takes 3 seconds instead of 1 second | 72 hours for triage; fix within 1 week; targeted regression tests added | **Yes** — no release with an open High defect |
| **Medium** | Functional defect with workaround; UNCERTAIN verdict where PASS or FAIL is achievable; cosmetic certificate formatting issue; CLI error message missing a suggested remediation | Type system returns UNCERTAIN for a property where enough information exists to return PASS; certificate JSON has incorrect indentation; error message lacks remediation suggestion | 1 week for triage; fix within 2 weeks; test added to prevent recurrence | **No** — but must be documented in release notes |
| **Low** | Minor usability issue; documentation error; non-user-facing logging issue; CLI help text typo | Help text for `analyze-orf` subcommand has a typo; debug log message has incorrect component ID; SRS cross-reference in an error message points to wrong section | 2 weeks for triage; fix within 1 month | **No** |

**Defect lifecycle:**

1. **Discovery**: Tester or developer identifies a deviation from expected behavior.
2. **Classification**: Defect is assigned a severity by the test lead, with input from the developer who authored the affected component.
3. **Assignment**: Defect is assigned to the responsible developer.
4. **Resolution**: Developer implements a fix and adds a regression test.
5. **Verification**: Test lead verifies the fix resolves the defect and the regression test passes.
6. **Closure**: Defect is closed if verification passes; reopened if it fails.

**Escalation**: Any defect that is not resolved within its SLA is escalated to the project lead. Critical defects unresolved after 72 hours trigger a project-wide review.

---

## 7. Traceability: Requirements → Test Cases

The following table provides complete traceability from every functional requirement (REQ-FUNC-XXX), non-functional requirement (REQ-NFR-XXX), and constraint (REQ-CON-XXX) in the SRS (DOC-01) to the test cases in this SVVP that verify it. Every requirement must be covered by at least one test case.

| Requirement ID | Requirement Summary | Unit Tests | Integration Tests | System Tests | Soundness Tests | Repro. Tests |
|---|---|---|---|---|---|---|
| REQ-FUNC-001 | Accept FASTA/raw input; reject non-IUPAC | TC-U-001, TC-U-002, TC-U-003 | — | — | — | — |
| REQ-FUNC-002 | Scan and annotate biological elements | TC-U-004, TC-U-005, TC-U-006 | TC-I-001 | — | — | — |
| REQ-FUNC-003 | Scanner implemented as DFAs | TC-U-004, TC-U-005, TC-U-006 | TC-I-001 | — | — | — |
| REQ-FUNC-004 | Scanner determinism | TC-U-007 | — | — | — | TC-REP-001, TC-REP-002 |
| REQ-FUNC-010 | NDFST for splice isoform set | TC-U-010 | TC-I-001 | — | — | — |
| REQ-FUNC-011 | Splicing grammar rules (GT-AG, GC-AG, AT-AC) | TC-U-012, TC-U-013 | TC-I-001 | — | — | — |
| REQ-FUNC-012 | Alternative splicing as non-deterministic branching | TC-U-011 | TC-I-002 | — | — | — |
| REQ-FUNC-013 | Cellular context parameterization | TC-U-014 | TC-I-001 | — | — | — |
| REQ-FUNC-014 | Splicing Engine determinism | TC-U-015 | — | — | — | TC-REP-001, TC-REP-002 |
| REQ-FUNC-015 | NDFST output fields (sequence, boundaries, frame, annotations, provenance) | TC-U-010, TC-U-011 | TC-I-001, TC-I-002 | — | — | — |
| REQ-FUNC-020 | Translation as deterministic FST | TC-U-020 | TC-I-002 | — | — | — |
| REQ-FUNC-021 | Special cases: selenocysteine, frameshifting | TC-U-021, TC-U-023 | TC-I-003 | TC-S-003 | — | — |
| REQ-FUNC-022 | Translation FST determinism | TC-U-022 | — | — | — | TC-REP-001 |
| REQ-FUNC-023 | Translation output fields | TC-U-020, TC-U-021, TC-U-023 | TC-I-002, TC-I-004 | — | — | — |
| REQ-FUNC-030 | Three-valued type checker (PASS/FAIL/UNCERTAIN) | TC-U-030, TC-U-031 | TC-I-004 | — | TC-SND-001 through TC-SND-008 | — |
| REQ-FUNC-031 | Type predicates (SpliceCorrect, NoCrypticSplice, etc.) | TC-U-030, TC-U-031 | TC-I-004, TC-I-005 | TC-S-001 | TC-SND-001 through TC-SND-007 | — |
| REQ-FUNC-032 | Type checker determinism | TC-U-032 | — | — | — | TC-REP-001 |
| REQ-FUNC-033 | Subtyping relations | TC-U-033 | — | — | — | — |
| REQ-FUNC-034 | Three-valued composition rules | TC-U-034, TC-U-035 | — | TC-S-008 (TC-SND-008) | TC-SND-008 | — |
| REQ-FUNC-035 | Verdict evidence (derivation, violation, knowledge gap) | TC-U-030, TC-U-031 | TC-I-004 | — | TC-SND-001 through TC-SND-007 | — |
| REQ-FUNC-040 | CSP solver for codon optimization | TC-U-040 | TC-I-005, TC-I-006 | TC-S-001 | — | — |
| REQ-FUNC-041 | Decision variable domains | TC-U-040, TC-U-043 | TC-I-005 | — | — | — |
| REQ-FUNC-042 | Hard constraints (splicing, CAI, GC, restriction, frame, instability) | TC-U-042 | TC-I-005 | TC-S-001 | — | — |
| REQ-FUNC-043 | Feasible problem: return assignment + verification | TC-U-040, TC-U-043 | TC-I-005, TC-I-006 | TC-S-001 | — | — |
| REQ-FUNC-044 | Infeasible problem: MUS with minimality | TC-U-041, TC-U-044 | TC-I-005 | TC-S-004 | — | — |
| REQ-FUNC-045 | Optimizer determinism | TC-U-045 | — | — | — | TC-REP-003 |
| REQ-FUNC-050 | FFI: folding adapter interface | TC-U-051, TC-U-052, TC-U-053 | TC-I-003 | — | — | — |
| REQ-FUNC-051 | FFI: PTM adapter interface | TC-U-051, TC-U-052, TC-U-053 | TC-I-003 | — | — | — |
| REQ-FUNC-052 | FFI: correct input formatting, output parsing, provenance | TC-U-050, TC-U-051, TC-U-052, TC-U-053 | TC-I-003 | — | — | — |
| REQ-FUNC-053 | FFI non-determinism treatment | TC-U-051, TC-U-052 | TC-I-003 | — | — | — |
| REQ-FUNC-060 | Guarantee certificate generation | TC-U-060 | TC-I-006, TC-I-007 | TC-S-001 | — | — |
| REQ-FUNC-061 | Independent certificate verification | TC-U-061 | TC-I-006, TC-I-007 | TC-S-001 | — | TC-REP-004 |
| REQ-FUNC-062 | Circuit-level certificate | TC-U-062 | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-070 | Circuit specification input | — | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-071 | Composition checks (promoter, resource, splicing, RNA) | TC-U-070, TC-U-071, TC-U-072, TC-U-073 | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-072 | Composition verdict evidence | TC-U-070, TC-U-071, TC-U-072, TC-U-073 | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-073 | Three-valued composition of circuit verdicts | TC-U-073 | TC-I-008 | TC-S-002 | — | — |
| REQ-FUNC-080 | Mutation space decomposition | TC-U-080 | — | — | — | — |
| REQ-FUNC-081 | Legal mutation enumeration | TC-U-081 | — | — | — | — |
| REQ-FUNC-082 | Independence exploitation | TC-U-082 | — | — | — | — |
| REQ-FUNC-083 | Constraint conflict detection | TC-U-083 | — | — | — | — |
| REQ-FUNC-090 | Multi-frame translation FST construction | TC-U-090 | TC-I-009 | TC-S-003 | — | — |
| REQ-FUNC-091 | Shared constraint set computation | TC-U-091 | TC-I-009 | TC-S-003 | — | — |
| REQ-FUNC-092 | Coupling classification (high/low) | TC-U-092 | TC-I-009 | TC-S-003 | — | — |
| REQ-FUNC-093 | Inter-frame constraint conflict detection | TC-U-093 | TC-I-009 | TC-S-003 | — | — |
| REQ-NFR-001 | Scanner performance: ≤ 1 s for 10 kb | TC-U-007 (timed) | — | — | — | — |
| REQ-NFR-002 | Splicing performance: ≤ 5 s for average gene | TC-U-015 (timed) | — | — | — | — |
| REQ-NFR-003 | Translation performance: ≤ 100 ms for 5 kb | TC-U-022 (timed) | — | — | — | — |
| REQ-NFR-004 | Type checker performance: ≤ 10 s for all predicates | TC-U-032 (timed) | — | — | — | — |
| REQ-NFR-005 | Optimizer performance: ≤ 60 s for 1000 aa | TC-U-045 (timed) | — | — | — | — |
| REQ-NFR-006 | Composition verifier performance: ≤ 5 min for 10 genes | TC-U-071 (timed) | — | — | — | — |
| REQ-NFR-010 | Determinism for all internal stages | TC-U-007, TC-U-015, TC-U-022, TC-U-032, TC-U-045 | — | — | — | TC-REP-001, TC-REP-002, TC-REP-003, TC-REP-004 |
| REQ-NFR-011 | Soundness: no false PASS | TC-U-031 | — | — | TC-SND-001 through TC-SND-008 | — |
| REQ-NFR-012 | Completeness: solver finds feasible solution | TC-U-040, TC-U-043 | TC-I-005 | TC-S-001 | — | — |
| REQ-NFR-013 | Correctness: INFEASIBLE is truly infeasible; MUS valid | TC-U-041, TC-U-044 | TC-I-005 | TC-S-004 | — | — |
| REQ-NFR-020 | CLI with subcommands and exit codes | — | — | TC-S-001, TC-S-002, TC-S-003, TC-S-004 | — | — |
| REQ-NFR-021 | Python API with type-annotated signatures | TC-U-001 through TC-U-093 (all invoked via Python API) | TC-I-001 through TC-I-009 | — | — | — |
| REQ-NFR-022 | Error messages with REQ ID, position, remediation | TC-U-003, TC-U-050, TC-U-051, TC-U-052 | — | TC-S-004 | — | — |
| REQ-NFR-023 | UNCERTAIN verdicts with knowledge gap description | TC-U-034, TC-U-073 | — | — | — | — |
| REQ-NFR-030 | IR schemas in protocol buffers; schema evolution | — | TC-I-001 through TC-I-009 (all use protobuf IR) | — | — | — |
| REQ-NFR-031 | Each stage independently testable | TC-U-001 through TC-U-093 (each tests one component) | TC-I-001 through TC-I-009 | — | — | — |
| REQ-NFR-032 | Declarative grammar configuration | TC-U-014, TC-U-015 | TC-I-001 | — | — | — |
| REQ-NFR-033 | Common FFI adapter interface | TC-U-050, TC-U-051, TC-U-052, TC-U-053 | TC-I-003 | — | — | — |
| REQ-NFR-040 | Linux x86_64 and ARM64 support | — | — | — | — | TC-REP-002 |
| REQ-NFR-041 | macOS x86_64 and Apple Silicon support | — | — | — | — | TC-REP-002 |
| REQ-NFR-042 | Core pipeline works without external tools | TC-U-001 through TC-U-093 (none require FFI tools installed) | TC-I-007 (excluding TC-I-003) | TC-S-001, TC-S-004 | — | — |
| REQ-NFR-050 | Certificate integrity via SHA-256 hash | TC-U-060, TC-U-061 | TC-I-006, TC-I-007 | TC-S-001 | — | — |
| REQ-NFR-051 | Input validation: reject non-IUPAC characters | TC-U-003 | — | — | — | — |
| REQ-NFR-052 | No network access for core pipeline | TC-U-001 through TC-U-093 (all run locally) | TC-I-001, TC-I-002, TC-I-004 through TC-I-009 | TC-S-001 through TC-S-004 | — | — |
| REQ-CON-001 | No probabilistic models in internal stages | TC-U-007, TC-U-015, TC-U-022, TC-U-032, TC-U-045 (determinism implies no randomness) | — | — | — | TC-REP-001 |
| REQ-CON-002 | No grammar induction | — | — | — | — | — |
| REQ-CON-003 | No internal folding/PTM models | TC-U-050, TC-U-051, TC-U-052, TC-U-053 (FFI tests verify external-only handling) | TC-I-003 | — | — | — |
| REQ-CON-004 | No claim biology implements compilation | — | — | — | — | — |
| REQ-CON-010 | Memory ≤ 32 GB (single gene), ≤ 64 GB (circuit) | TC-U-040 (memory profiled) | TC-I-008 (memory profiled) | TC-S-002 (memory profiled) | — | — |
| REQ-CON-011 | GPU only for FFI; core pipeline CPU-only | TC-U-001 through TC-U-093 (no GPU required) | TC-I-007 (no GPU required) | TC-S-001, TC-S-004 (no GPU required) | — | — |

**Coverage summary**: All 30 functional requirements (REQ-FUNC-001 through REQ-FUNC-093), 17 non-functional requirements (REQ-NFR-001 through REQ-NFR-052), and 6 constraints (REQ-CON-001 through REQ-CON-011) are covered by at least one test case. Requirements that are design constraints rather than testable behaviors (REQ-CON-002: no grammar induction; REQ-CON-004: no claim biology implements compilation) are verified by code review and architectural inspection rather than automated testing.

---

*End of DOC-05: Software Verification & Validation Plan*
