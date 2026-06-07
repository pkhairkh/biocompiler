# BioCompiler Refinement Mapping: Lean4 → Python

## Purpose

This document maps each Lean4 theorem to its Python implementation, establishing a traceable connection between the formal model and the running code. The central caveat of the BioCompiler project is that "the Lean4 proofs are about a simplified model, not the actual Python implementation." This refinement mapping makes that gap explicit, theorem by theorem, and describes how each gap could be closed.

**Scope**: This document covers all 11+ proof modules in `proof/BioCompiler/`, and the corresponding Python modules in `src/biocompiler/`. (The former `lean/BioCompiler/` library modules were removed as redundant; their content is superseded by the proof/ modules.)

---

## 1. Type Correspondence Table

The foundational challenge of the refinement mapping is that the Lean4 formal model and the Python implementation use different type systems. The table below establishes the correspondence.

| Lean4 Type | Python Counterpart | Refinement Notes |
|---|---|---|
| `BioCompiler.Verdict` (3 values: PASS, FAIL, UNCERTAIN) | `biocompiler.types.Verdict` (5 values: PASS, LIKELY_PASS, UNCERTAIN, LIKELY_FAIL, FAIL) | Python has 2 extra values (LIKELY_PASS, LIKELY_FAIL). The 3-valued Lean4 model is a conservative abstraction: every Python PASS/LIKELY_PASS maps to Lean4 PASS; every Python FAIL/LIKELY_FAIL maps to Lean4 FAIL. See Gap §3.1. |
| `BioCompiler.Sequence` (= `List Nucleotide`) | `str` | Lean4 uses a typed list of `Nucleotide` constructors; Python uses raw strings. The refinement requires that every `str` character is in `{A, C, G, T}` — enforced at runtime by `scanner.validate_dna_sequence()`. |
| `BioCompiler.Nucleotide` | Single character in `str` | `Nucleotide.A` ↔ `'A'`, etc. Python has no compile-time guarantee that strings contain only valid bases. |
| `BioCompiler.CellularContext` | Organism parameters in `type_system/` package + `organisms/*.py` | Lean4 bundles cellType, ESE/ESS/ISE/ISS thresholds into one record. Python distributes these across species config dicts (`SPECIES`, `CODON_ADAPTIVENESS_TABLES`). |
| `BioCompiler.TypePredicate` (33 constructors) | `PREDICATE_NAMES` list + `check_*` functions | Lean4 uses an inductive type with 33 constructors (13 core + 20 SLOT). Python uses string-based dispatch and a registry pattern. |
| `BioCompiler.CertLevel` (GOLD/SILVER/BRONZE) | `type_system.CertLevel` (GOLD/SILVER/BRONZE) | Direct 1:1 correspondence. Both classify optimization quality. |
| `BioCompiler.SpliceVerdict` (PASS/UNCERTAIN/FAIL) | `type_system.SpliceVerdict` (PASS/UNCERTAIN/FAIL) | Direct 1:1 correspondence for the dual-threshold splice model. |
| `BioCompiler.SatisfactionMethod` | `MutagenesisProposal.chose_poorly` / `.impossible` flags | Lean4 uses a 3-constructor inductive type; Python uses boolean flags on a dataclass. |
| `BioCompiler.PredicateRecord` | `type_system.PredicateResult` | Lean4: name + passed + method + optional verdict. Python: predicate + passed + verdict + details + positions. Python carries more diagnostic data. |
| `BioCompiler.NDFST State` | NDFST logic in `splicing.py` | Lean4 models the NDFST as a parametric state machine with `ndfstRun`/`ndfstOutputSet`. Python uses a procedural implementation with the same semantics but no type-level state parameter. |
| `BioCompiler.SpliceIsoform` | `types.SpliceIsoform` | Lean4: `sequence` + `exonBoundaries`. Python: `sequence` + `exon_boundaries` + `parse_path` + `score`. Python has extra fields. |
| `BioCompiler.CodonAdaptationIndex` (typeclass) | `organisms.SPECIES` + `CODON_ADAPTIVENESS_TABLES` | Lean4 axiomatizes `computeCAI` as deterministic. Python uses empirical CAI tables from organism databases. |
| `BioCompiler.SLOTValues` | Individual FFI adapter outputs (ESMFold, NetPhos, etc.) | Lean4 models SLOT as a record of optional values. Python fills these via external tool calls at runtime. |
| `Rat` (Lean4 arbitrary precision) | `float` (Python IEEE 754) | Critical gap: Lean4 uses exact rational arithmetic; Python uses floating point. See Gap §3.4. |

---

## 2. Theorem-by-Theorem Mapping

### 2.1 Three-Valued Logic Module (`ThreeValued.lean`)

This module defines the algebraic properties of the three-valued verdict logic. It has 12 theorems, 0 sorry, 0 axioms.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 1 | `Verdict.and_PASS` — `Verdict.and PASS v = v` | `five_valued_and(Verdict.PASS, v)` returns `v` | Python has 5 values; the identity holds for all 5 when the first argument is PASS (highest in ordering) | Property-based test: ∀ v ∈ Verdict, `five_valued_and(PASS, v) == v` |
| 2 | `Verdict.and_FAIL_left` — `Verdict.and FAIL v = FAIL` | `five_valued_and(Verdict.FAIL, v)` returns `FAIL` | FAIL is lowest in ordering, so min(FAIL, anything) = FAIL holds for 5-valued logic too | Property-based test: ∀ v ∈ Verdict, `five_valued_and(FAIL, v) == FAIL` |
| 3 | `Verdict.and_FAIL_right` — `Verdict.and v FAIL = FAIL` | `five_valued_and(v, Verdict.FAIL)` returns `FAIL` | Same as above but commuted | Property-based test |
| 4 | `and_eq_PASS_iff` — `and init hd = PASS ↔ init = PASS ∧ hd = PASS` | Implicit in `combined_verdict()` logic | The 5-valued AND still has this property: only PASS ∧ PASS = PASS | Property-based test: ∀ a, b, `five_valued_and(a, b) == PASS ↔ a == PASS and b == PASS` |
| 5 | `and_comm` | Commutativity of `five_valued_and` | Holds for min-based ordering | Property-based test |
| 6 | `or_comm` | Commutativity of `five_valued_or` | Holds for max-based ordering | Property-based test |
| 7 | `and_assoc` | Associativity of `five_valued_and` | Holds for min | Property-based test |
| 8 | `or_assoc` | Associativity of `five_valued_or` | Holds for max | Property-based test |
| 9 | `and_idem` | `five_valued_and(v, v) == v` | Holds for min | Property-based test |
| 10 | `and_pass_pass` — `v₁ = PASS → v₂ = PASS → and v₁ v₂ = PASS` | Trivially holds | No gap | N/A |
| 11 | `foldl_and_pass_implies_all_pass` | `combined_verdict(vs) == PASS → ∀ v ∈ vs, v == PASS` | The 5-valued version must also have this: if the fold of min is PASS (4), every element must be 4 | Property-based test |
| 12 | `foldl_ne_pass_of_ne_pass` | If initial ≠ PASS, fold never reaches PASS | Holds: min of non-max with anything stays non-max | Property-based test |

**Refinement summary for ThreeValued**: The 3-valued logic is a conservative abstraction of the 5-valued logic. Every theorem that holds for 3 values also holds for the corresponding subset of 5 values. The 5-valued logic adds `LIKELY_PASS` (3) and `LIKELY_FAIL` (1) between PASS (4) and UNCERTAIN (2), and between UNCERTAIN (2) and FAIL (0). All lattice properties are preserved because the ordering is total and AND = min / OR = max.

---

### 2.2 Sequence Module (`Sequence.lean`)

This module defines nucleotide sequences and pattern matching. 3 theorems, 0 sorry, 0 axioms.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 13 | `matchesAt_spec` — `matchesAt seq pattern pos = true ↔ matchesAtProp seq pattern pos` | `seq[pos:pos+len(pattern)] == pattern` in Python | Lean4 uses `List Nucleotide` equality; Python uses string slicing. Semantically equivalent when string is valid DNA. | Property-based test with random DNA sequences |
| 14 | `containsPattern_complete` — if pattern appears, scanner returns true | `pattern in seq` or `seq.find(pattern)` in Python | Python's `str.find` has identical completeness: if the pattern exists at some position, `find` returns that position. | Unit test: ∀ pattern, seq, `pattern in seq ↔ _find_pattern(seq, pattern) != -1` |
| 15 | `containsPattern_sound` — if scanner returns true, pattern appears | Reverse of completeness | Same | Same |

**Additional functions**:

- `gcContent` → `scanner.gc_content()` / `type_system` GC checks. Lean4 uses `Rat`; Python uses `float`. The algorithm is identical: `(G+C)/total`. Gap: floating-point rounding.
- `hasPrematureStop` → `type_system.check_no_stop_codons()`. Lean4 checks all codons in the reading frame; Python does the same with `seq[i:i+3] in ("TAA", "TAG", "TGA")`.
- `stopCodons` → `{"TAA", "TAG", "TGA"}` in Python. Direct correspondence.

---

### 2.3 NDFST Module (`NDFST.lean`)

This module defines non-deterministic finite-state transducer semantics. 4 main theorems, 0 sorry, 0 axioms.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 16 | `ndfstRun_sound` — every output in `ndfstRun` corresponds to a valid path | `splicing.py` NDFST execution | Python's NDFST follows the same foldl semantics. Gap: Python uses mutable lists; Lean4 uses purely functional lists. | Integration test: run NDFST on known gene sequences, verify output matches Lean4 model |
| 17 | `ndfstRun_complete` — every valid path produces an output in `ndfstRun` | Same | Same | Same |
| 18 | `ndfst_deterministic` — `ndfstOutputSet` is deterministic (reflexivity) | Python NDFST is deterministic given the same input | No semantic gap; this is a computational determinism guarantee | Test: same input → same output (idempotence) |
| 19 | `ndfst_unique_deterministic` — `ndfstUniqueOutputSet` is deterministic | Same | Same | Same |

**Key structure mappings**:
- `CellularContext` → organism parameters in Python (cell type, splice thresholds)
- `SpliceIsoform` → `types.SpliceIsoform`
- `SplicingNDFST` typeclass → `splicing.py` NDFST implementation

---

### 2.4 Scanners Module (`Scanners.lean`)

This module implements concrete scanner functions with completeness proofs. Contains both concrete scanners (fully proved) and abstract scanner interfaces (axiomatized).

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 20 | `hasPattern_complete` — pattern found → `hasPattern` returns true | `seq.find(pattern)` / `pattern in seq` | Direct correspondence | Unit test |
| 21 | `hasPattern_sound` — `hasPattern` true → pattern found | Same | Same | Same |
| 22 | `hasAnyRestrictionSite_complete` | `check_no_restriction_site()` | Lean4: list of Sequence patterns. Python: enzyme name → site lookup + `str.find`. Python also checks reverse complement strand. | Integration test: verify every restriction site found by Python scanner is real |
| 23 | `hasAnyRestrictionSite_sound` | Same | Same | Same |
| 24 | `hasInstabilityMotif_attta_complete` | `scanner.scan_sequence()` for "ATTTA" motifs | Direct correspondence; Python scans for `INSTABILITY_MOTIF = "ATTTA"` | Unit test |
| 25 | `hasInstabilityMotif_urich_complete` | Python scans for T-runs ≥ 6 | Lean4: `uRichMotif = tRun 6`. Python: consecutive T detection in optimizer. | Unit test |
| 26 | `hasInstabilityMotif_complete` (contrapositive) | N/A (Python doesn't have a negative completeness check) | Python has no explicit "no motif exists" verification | Could add property-based test |
| 27 | `hasInstabilityMotif_sound` | Same | Same | Same |

**Abstract scanner axioms** (these are PARAMETERS of the proof, not gaps):

| # | Axiom | Python Implementation | Gap Assessment |
|---|---|---|---|
| A1 | `SpliceSiteScanner.scanner_completeness` | `maxentscan.score_donor()` / `score_acceptor()` | Python's MaxEntScan scoring is a well-validated computational model. Gap: the MaxEntScan model itself is heuristic, not formally verified. |
| A2 | `SpliceSiteScanner.scanner_soundness` | Same | Same |
| A3 | `SpliceSiteScanner.borderline_completeness` | Dual-threshold logic in `check_no_cryptic_splice()` | Python implements `low_thresh` and `high_thresh` exactly as the formal model specifies |
| A4 | `CpGIslandScanner.scanner_completeness` | `check_no_cpg_island()` sliding window | Python uses identical window/threshold logic |
| A5 | `CpGIslandScanner.scanner_soundness` | Same | Same |
| A6 | `PromoterScanner.scanner_completeness` | `check_no_cryptic_promoter()` | Python uses IUPAC-aware consensus scoring with position weight matrices — more sophisticated than the abstract model |
| A7 | `PromoterScanner.scanner_soundness` | Same | Same |
| A8 | `PromoterScanner.borderline_completeness` | Dual-threshold `threshold * 0.8` | Python implements exactly this: UNCERTAIN if `score >= threshold * 0.8` but `< threshold` |
| A9 | `TMDomainScanner.scanner_completeness` | `check_no_unexpected_tm_domain()` | Python uses sliding window hydrophobic fraction, identical to the formal model's `hydroFraction ≥ threshold` |
| A10 | `TMDomainScanner.scanner_soundness` | Same | Same |
| A11 | `TMDomainScanner.borderline_completeness` | `threshold * 0.85` borderline | Python implements: UNCERTAIN if `frac > threshold * 0.85` but `frac ≤ threshold` |
| A12 | `mRNAStructureOracle.oracle_completeness` | `check_mrna_secondary_structure()` | Python uses simplified nearest-neighbor ΔG estimation. Gap: not a full ViennaRNA computation. |
| A13 | `mRNAStructureOracle.borderline_completeness` | `dg_threshold * 0.7` borderline | Python implements this threshold |
| A14 | `CoTranslationalFoldingOracle.oracle_completeness` | `check_co_translational_folding()` | Python uses codon ramp analysis + pause site detection. Gap: heuristic, not a full co-translational folding simulation. |
| A15 | `CoTranslationalFoldingOracle.borderline_completeness` | Same | Same |
| A16 | `SplicingNDFST.output_is_valid` | NDFST in `splicing.py` | Python's NDFST produces the same valid isoforms |
| A17 | `SplicingNDFST.all_isoforms_produced` | Same | Same |
| A18 | `CodonAdaptationIndex.computeCAI` | `organisms.SPECIES` CAI tables | Python uses empirical CAI values; the deterministic computation property holds |

---

### 2.5 Type System Module (`TypeSystem.lean`)

The central module: defines all 33 type predicates and proves `type_soundness`. 1 main theorem + 1 per-predicate proof, 0 sorry.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 28 | **`type_soundness`** — `evaluate P seq ctx = PASS → propertyHolds P seq ctx` | Each `check_*` function returns `PredicateResult` with verdict | **This is the central theorem.** Python's `check_*` functions implement the same evaluation logic but as imperative code. The soundness guarantee is: if a Python check returns PASS, the property truly holds. | Property-based testing: for each predicate, verify that PASS result implies the property holds (by independent verification) |

**Per-predicate soundness proofs** (each is a case in `type_soundness`):

| Predicate | Lean4 evaluate logic | Python `check_*` function | Key Gap |
|---|---|---|---|
| `SpliceCorrect cellType` | `ctx.cellType != cellType → UNCERTAIN; else unique output length = 1 → PASS` | Not directly implemented as a standalone check in Python | Python's splice correctness is verified implicitly via `check_no_cryptic_splice` |
| `NoCrypticSplice` | `hasCrypticSpliceSite → FAIL; hasBorderlineSpliceSite → UNCERTAIN; else PASS` | `check_no_cryptic_splice(seq, low, high)` | Python uses MaxEntScan scoring; Lean4 abstracts via SpliceSiteScanner. Same dual-threshold logic. |
| `CodonAdapted org threshold` | `computeCAI ≥ threshold → PASS` | `check_codon_optimality()` | Python uses per-codon CAI; Lean4 uses whole-sequence CAI. Threshold semantics match. |
| `GCInRange lo hi` | `lo ≤ gcContent ≤ hi → PASS` | `gc_content()` used in optimizer constraints | Python uses `float`; Lean4 uses `Rat`. Potential rounding difference at boundaries. |
| `NoRestrictionSite enzymeSites` | `hasAnyRestrictionSite → FAIL` | `check_no_restriction_site(seq, enzymes)` | Python also checks reverse complement; Lean4 model does not explicitly handle revcomp |
| `InFrame rf boundaries` | `readingFrameConsistent ∧ ¬hasPrematureStop → PASS` | `check_valid_coding_seq()` partially covers this | Python doesn't have a separate InFrame check; frame consistency is implicit |
| `NoInstabilityMotif` | `hasInstabilityMotif → FAIL` | `scanner.scan_sequence()` detects ATTTA | Direct correspondence |
| `NoCpGIsland` | `hasCpGIsland → FAIL` | `check_no_cpg_island(seq, window, threshold)` | Same window/threshold logic; Python uses float division |
| `NoGTDinucleotide` | `hasPattern seq spliceDonorConsensus → FAIL` | `check_no_gt_dinucleotide(seq)` and `check_no_avoidable_gt(seq)` | Python has both strict (any GT → FAIL) and relaxed (only avoidable GT → FAIL) versions. Lean4 models only the strict version. |
| `NoStopCodons` | `hasPrematureStop → FAIL` | `check_no_stop_codons(seq)` | Direct correspondence; both skip the last codon |
| `ValidCodingSeq` | `length % 3 = 0 ∧ ¬hasPrematureStop` | `check_valid_coding_seq(seq)` | Direct correspondence |
| `CodonOptimality org threshold` | `computeCAI ≥ threshold → PASS` | `check_codon_optimality()` | Same as CodonAdapted |
| `NoCrypticPromoter org threshold` | `hasCrypticPromoter → FAIL; hasBorderlinePromoter → UNCERTAIN; else PASS` | `check_no_cryptic_promoter(seq, org, threshold)` | Python uses IUPAC-aware consensus scoring; Lean4 abstracts via PromoterScanner. Same dual-threshold. |
| 20 SLOT predicates | Always `UNCERTAIN` | Each has a Python `check_*` function that may return PASS/FAIL | **Major gap**: Lean4 models SLOT predicates as always UNCERTAIN (conservative). Python actually evaluates them with heuristic engines. See Gap §3.3. |

---

### 2.6 Compositional Soundness Module (`Compositional.lean`)

Proves that the overall PASS verdict implies all individual properties hold.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 29 | **`compositional_soundness`** — `evaluateAll = PASS → ∀ P, propertyHolds P` | `combined_verdict(results)` returns PASS only if all pass | Python's `combined_verdict` uses `five_valued_and` (min over 5-valued ordering). If the result is PASS (4), every element must be PASS (4). | Property-based test: `combined_verdict(vs) == PASS ↔ all(v == PASS for v in vs)` |
| 30 | `slot_predicates_uncertain` — SLOT predicates never produce PASS | In Python, SLOT predicates CAN produce PASS (e.g., `check_conservation_score` returns PASS/FAIL) | **Critical gap**: Lean4 says SLOT predicates always return UNCERTAIN; Python evaluates them and returns PASS/FAIL. This means Python certificates may claim PASS for SLOT predicates that the formal model says are uncertain. | This is a deliberate design choice: the formal model is conservative. Python provides "best effort" evaluation. The guarantee certificate should treat SLOT PASS as UNCERTAIN for formal purposes. |
| 31 | `slot_predicates_dont_affect_pass` — SLOT in list → no overall PASS | If Python evaluates a SLOT predicate to PASS, it contributes to overall PASS | Same as above | Certificate generation in Python should have a mode that treats SLOT predicates as UNCERTAIN |
| 32 | `all_core_if_pass` — evaluateAll = PASS → all predicates are core | Python does not enforce this | Same as above | Add validation in `generate_certificate()` |
| 33 | `dinucleotide_no_compose` — constraints don't compose (counterexample) | `_remove_site_multicodon()` handles junction-crossing sites | Python correctly handles this via multi-codon coordinated solving | Integration test: verify junction-crossing sites are detected |
| 34 | `restriction_site_no_compose` — same for restriction sites | Same | Same | Same |

---

### 2.7 Splicing Resolution Module (`SplicingResolution.lean`)

Proves key results about splicing resolution correctness.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 35 | `pass_implies_no_cryptic_sites` | `check_no_cryptic_splice()` returning PASS implies no cryptic sites | Python uses MaxEntScan which is a heuristic; the formal model assumes scanner completeness as an axiom | Cross-validate MaxEntScan results against known splice site databases |
| 36 | `fail_iff_cryptic_exists` | FAIL iff the scanner finds a cryptic site | Same | Same |
| 37 | `uncertain_iff_borderline_no_cryptic` | UNCERTAIN iff borderline but no cryptic | Same | Same |
| 38 | **`splice_resolution_deterministic`** — both PASS → unique isoform | Python's NDFST produces a deterministic set of isoforms | Python's NDFST implementation should produce exactly one isoform when splice is correct | Integration test: verify NDFST produces exactly 1 isoform for known single-isoform genes |
| 39 | `canonical_donor_has_gt` (axiom) | Splice donor sites always have GT at the boundary | Biological fact, not a computational guarantee | Validate against GENCODE annotation |
| 40 | `canonical_acceptor_has_ag` (axiom) | Same for AG | Same | Same |
| 41 | `no_cryptic_splice_verdicts_exclusive` — PASS/UNCERTAIN/FAIL are exclusive and exhaustive | Python's dual-threshold logic produces the same exclusive verdicts | Direct correspondence | Property-based test |
| 42 | `pass_both_implies_canonical_only` | Python checks this implicitly | Same as #38 | Same |
| 43 | `extension_cannot_remove_gt` | Python's optimization doesn't rely on this property explicitly | The theorem states a monotonicity property: adding nucleotides can only add GT sites, never remove them | Could be used to optimize the Python optimizer (early termination) |

---

### 2.8 Mutagenesis Module (`Mutagenesis.lean`)

Proves key results about type-directed mutagenesis. 15+ theorems, 0 sorry.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 44 | `synonymous_preserves_translation` | Optimization always preserves the amino acid sequence (invariant in optimizer) | Python's optimizer has assertions for this (`len(sequence) == len(protein) * 3`) | Property-based test: after optimization, translate back and verify same protein |
| 45 | `synonymous_gc_counterexample` — synonymous sub can change GC | `_find_gt_free_codons()` and GC adjustment in optimizer | Python correctly handles this: GC is re-checked after every substitution | Unit test: verify AAA→AAG changes GC count |
| 46 | `synonymous_restriction_counterexample` — synonymous sub can create restriction sites | `_remove_site_multicodon()` handles this | Python handles this via reconciliation passes | Integration test |
| 47 | `all_valine_codons_have_gt` | `GT_MANDATORY_AAS = {"V"}` in `mutagenesis.py` | Direct correspondence | Unit test: verify all 4 Valine codons contain "GT" |
| 48 | `mandatory_gt_has_gt` | `find_unrepairable_cryptic_donors()` identifies Valine positions | Direct correspondence | Unit test |
| 49 | `unrepairable_cryptic_donor_exists` | `find_unrepairable_cryptic_donors()` reports `fixable=False` for Valine positions | Direct correspondence | Integration test with known Valine-containing sequences |
| 50 | `synonymous_safe_for_SLOT_predicates` | SLOT predicates are trivially preserved (both sides UNCERTAIN in Lean4; in Python, SLOT checks are independent of sequence) | Gap: Python SLOT predicates DO depend on sequence, but Lean4 says they're always UNCERTAIN. In practice, the Python SLOT predicates that depend on external tools (ESMFold, etc.) ARE non-deterministic, so the Lean4 model is conservative. | N/A (conservative abstraction) |
| 51 | `synonymous_unsafe_for_dna_predicates` | `_greedy_optimize()` shows that GT elimination can introduce new GTs at other positions | Python handles this via iterative re-scanning after each substitution | Integration test: verify optimizer doesn't create new violations |
| 52 | `valine_only_mandatory_gt_aa` | `GT_MANDATORY_AAS = {"V"}` | Direct correspondence | Unit test for all 20 amino acids |
| 53 | `every_aa_has_ag_free_codon` | `_find_ag_free_codons()` works for all AAs | Direct correspondence | Unit test for all 20 amino acids |

---

### 2.9 Certificates Module (`Certificates.lean`)

Proves certificate soundness. 1 main theorem, 0 sorry.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 54 | **`certificate_soundness`** — valid certificate → all properties hold | `verify_certificate()` re-evaluates each predicate | Python's verification re-runs the same `check_*` functions. Gap: if a `check_*` function has a bug, both generation and verification will be wrong. Lean4's proof is about the abstract model, not the Python code. | Independent re-implementation of critical checks for cross-validation |

---

### 2.10 SLOT Independence Module (`SLOTIndependence.lean`)

Proves that certificate validity is independent of FFI output. 6 main theorems.

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 55 | `evaluate_slot_independent` — evaluation doesn't take SLOT values as arguments | Python's `check_*` functions don't take SLOT values directly | Direct correspondence for core predicates | N/A |
| 56 | `property_slot_independent` | Same | Same | Same |
| 57 | `certificate_slot_independent` | `verify_certificate()` doesn't use SLOT values | Direct correspondence | N/A |
| 58 | `ffi_never_pass` — FFI-dependent predicates never produce PASS | Python's SLOT-dependent predicates CAN produce PASS | **Critical gap**: see §3.3 | Conservative certificate mode |
| 59 | `full_slot_independence` — combined guarantee | Python's certificate system doesn't have this guarantee explicitly | Need to add SLOT-aware certificate mode | Implement `formal_mode=True` in certificate generation that treats SLOT predicates as UNCERTAIN |
| 60 | `predicate_is_core_or_slot` | `PREDICATE_NAMES` list categorizes predicates | Python has 33 predicates; Lean4 has 33. The predicate lists are aligned. |

---

### 2.11 Lean4 Library Modules (now superseded by proof/ modules)

> **Note**: The former `lean/BioCompiler/` directory (CodonTable.lean, Predicates.lean, Optimization.lean, Certificate.lean) was a simplified v7.0.0 model that has been superseded by the `proof/BioCompiler/` modules. The theorem mappings below remain valid as they describe concepts now covered more rigorously in the proof/ modules.

#### CodonTable.lean (superseded by Sequence.lean + TypeSystem.lean)

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 61 | `translateCodon` — standard genetic code | `CODON_TABLE` dict in `type_system.py` | Direct correspondence: both implement the standard genetic code | Unit test: compare all 64 codon translations |
| 62 | `exists_gt_free_codon` — every non-Stop, non-V AA has a GT-free codon | `_find_gt_free_codons()` returns empty for V, non-empty for all others | Direct correspondence | Unit test for all 20 AAs |
| 63 | `exists_cg_free_codon` — every non-Stop AA has a CG-free codon | Python doesn't have an explicit CG-free codon finder, but `AvoidCpG` in the optimizer does this | Correspondence through the optimizer's CpG avoidance step | Integration test |
| 64 | `cross_codon_gt_resolvable` — cross-codon GT can be resolved when not M/W | Cross-codon GT resolution in optimizer | Python handles M/W as special cases (methionine ATG, tryptophan TGG both end in G) | Integration test |
| 65 | `cross_codon_cg_resolvable` — cross-codon CG always resolvable | CpG avoidance in optimizer | Direct correspondence | Integration test |

#### Predicates.lean (superseded by TypeSystem.lean)

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 66 | `noGT_subsumes_crypticSplice_FAIL` — NoGT ⟹ NoCrypticSplice(FAIL) | If no GT dinucleotides exist, `check_no_cryptic_splice` returns PASS | Direct correspondence: no GT means no splice donors, so NoCrypticSplice vacuously passes | Unit test |
| 67 | `validCoding_implies_noInternalStops` | `check_valid_coding_seq` passing implies `check_no_stop_codons` passing | Direct correspondence | Unit test |
| 68 | `dual_threshold_monotonicity` — PASS at strict threshold → PASS at permissive | `classifySplice` threshold comparison | Direct correspondence for float-based thresholds | Property-based test |

#### Optimization.lean (superseded by TypeSystem.lean + Mutagenesis.lean)

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 69 | `phase1_preserves_translation` — GT-aware optimization preserves protein | `_greedy_optimize()` preserves amino acid sequence | Python's optimizer has assertions for this; Lean4 proves it | Property-based test: optimize random proteins, verify translation preserved |
| 70 | `phase2_preserves_translation` — restriction site removal preserves protein | `_remove_site_multicodon()` preserves amino acid sequence | Same | Same |
| 71 | `phase5_preserves_translation` — CpG avoidance preserves protein | CpG avoidance step in optimizer | Same | Same |
| 72 | `pipeline_preserves_protein` — end-to-end pipeline preserves protein | `optimize_sequence()` overall pipeline | Same | End-to-end property-based test |

#### Certificate.lean (library, superseded by Certificates.lean)

| # | Lean4 Theorem | Python Counterpart | Refinement Gap | Gap Closure Strategy |
|---|---|---|---|---|
| 73 | `gold_implies_all_optimization` | `compute_certificate()` returns GOLD only if no mutagenesis and no unsatisfied | Direct correspondence | Unit test |
| 74 | `silver_implies_all_passed` | SILVER means all passed but some required mutagenesis | Direct correspondence | Unit test |
| 75 | `bronze_implies_unsatisfied` | BRONZE means at least one predicate failed | Direct correspondence | Unit test |
| 76 | `gold_certificate_soundness` | Python certificate generation doesn't have a formal soundness proof | **Gap**: Python relies on testing, not proof | Cross-validation with independent re-implementation |

---

## 3. Known Refinement Gaps

### 3.1 Three-Valued vs. Five-Valued Logic

**Description**: The Lean4 model uses 3 verdict values (PASS, UNCERTAIN, FAIL). Python uses 5 values (PASS, LIKELY_PASS, UNCERTAIN, LIKELY_FAIL, FAIL).

**Impact**: The 3-valued model is a conservative abstraction. Every Python PASS maps to Lean4 PASS. Every Python FAIL maps to Lean4 FAIL. LIKELY_PASS and LIKELY_FAIL are "promoted" to UNCERTAIN in the 3-valued model. This means the formal model is more conservative: it may claim UNCERTAIN where Python claims LIKELY_PASS.

**Closure Strategy**: Define a formal `refine5to3` function:
- `PASS → PASS`, `LIKELY_PASS → UNCERTAIN`, `UNCERTAIN → UNCERTAIN`, `LIKELY_FAIL → UNCERTAIN`, `FAIL → FAIL`
Prove that all 3-valued theorems still hold after this refinement. Add property-based tests that verify `combined_verdict` is monotone under this refinement.

### 3.2 Lean4 Sequence vs. Python str

**Description**: Lean4 uses `List Nucleotide` where each element is one of 4 constructors. Python uses `str` which can contain any Unicode character.

**Impact**: Python code must validate that strings contain only A/C/G/T. This is done at runtime by `scanner.validate_dna_sequence()`, but there's no compile-time guarantee.

**Closure Strategy**:
1. Add a `ValidDNA` newtype or type hint that constrains strings to DNA alphabet
2. Add property-based tests that verify all internal functions maintain the DNA-alphabet invariant
3. Consider using a `DNASequence` class with validation in `__init__`

### 3.3 Axiom Assumptions vs. Python Heuristics

**Description**: The Lean4 proof takes 18 scanner/oracle axioms as parameters. Python implements these as heuristic functions (MaxEntScan, Kozak scoring, sliding windows, etc.).

**Impact**: The soundness theorem says "IF the scanners are correct, THEN the type system is sound." The Python scanners are well-validated but not formally proved correct. This is the standard approach in formal methods: prove the system correct conditional on oracles, then validate the oracles independently.

**Closure Strategy**:
1. **Independent validation**: For each scanner, run against gold-standard databases (e.g., GENCODE for splice sites, REBASE for restriction enzymes)
2. **Property-based testing**: Generate random sequences, verify scanner completeness (if a real site exists, the scanner finds it) and soundness (if the scanner reports a site, it's real)
3. **Cross-validation**: Compare Python scanner results against alternative implementations (e.g., ViennaRNA for mRNA structure, TMHMM for TM domains)
4. **Empirical coverage**: Measure false positive/negative rates on benchmark datasets

### 3.4 Float vs. Rat Differences

**Description**: Lean4 uses `Rat` (exact rational arithmetic) for GC content, CAI thresholds, splice scores, etc. Python uses IEEE 754 `float`.

**Impact**: Boundary comparisons may differ due to floating-point rounding. For example, `gcContent seq = 0.6` might be `0.600000000001` in Python, causing `lo ≤ gcContent` to pass when the exact value would be `0.6`.

**Closure Strategy**:
1. Add epsilon tolerance in Python comparisons (e.g., `gc_content >= gc_lo - 1e-10`)
2. Property-based test: for random sequences, verify that `float(gc_content(seq))` agrees with `Rat(gc_content(seq))` within tolerance
3. Document which comparisons are boundary-sensitive

### 3.5 SLOT Predicates: Formal UNCERTAIN vs. Python PASS/FAIL

**Description**: In the Lean4 model, SLOT-dependent predicates always evaluate to UNCERTAIN. In Python, they are evaluated with heuristic engines and can return PASS or FAIL.

**Impact**: This is the largest semantic gap. The formal model says "we cannot guarantee SLOT predicates," but Python makes a best-effort evaluation. If a Python SLOT predicate returns PASS, the formal model considers this result unverified.

**Closure Strategy**:
1. Add a `formal_mode` flag to certificate generation that treats SLOT predicates as UNCERTAIN
2. Document which predicates are SLOT-dependent in the certificate itself
3. For guarantee certificates (GOLD/SILVER), require that only core predicates contribute to the PASS verdict
4. Add a `slot_confidence` field to `TypeCheckResult` that indicates the reliability of SLOT verdicts

### 3.6 Reverse Complement Handling

**Description**: The Lean4 model doesn't explicitly model reverse complement strands for restriction site checking. Python checks both strands.

**Impact**: Python is more thorough than the formal model. The formal model's soundness guarantee holds a fortiori for the Python implementation (if it passes the stricter check, it would also pass the weaker one).

**Closure Strategy**: This is a case where the Python implementation is STRONGER than the formal model. No gap closure needed; consider extending the formal model to include reverse complement.

### 3.7 Number of Predicates: 33 (Lean4) vs. 33 (Python)

**Description**: Lean4 defines 33 `TypePredicate` constructors (13 core + 20 SLOT). Python's `PREDICATE_NAMES` list has 33 entries.

**Impact**: The Lean4 and Python predicate counts are now aligned (33 each). Previously, 4 predicates (`StructureConfidence`, `NoMisfoldingRisk`, `CorrectFoldTopology`, `NoUnexpectedInteraction`) were in the Lean4 model but grouped or handled differently in Python; these have since been consolidated.

**Closure Strategy**: Verify that all 33 Lean4 predicates have corresponding Python checks, even if grouped differently. Add any missing checks.

---

## 4. Testing Strategy

### 4.1 Property-Based Tests for Refinement

The following property-based tests would verify the refinement mapping:

| Property | Test Description | Covers Theorems |
|---|---|---|
| `and_pass_implies_both_pass` | `five_valued_and(a, b) == PASS ↔ a == PASS ∧ b == PASS` | #4, #11 |
| `and_fail_absorbs` | `five_valued_and(FAIL, v) == FAIL` for all v | #2, #3 |
| `combined_verdict_pass_iff_all_pass` | `combined_verdict(vs) == PASS ↔ all(v == PASS for v in vs)` | #11, #29 |
| `gc_content_rational_agreement` | `abs(float_gc - rational_gc) < 1e-10` for all sequences | Float/Rat gap |
| `translation_preserved_after_optimization` | After `optimize_sequence(protein)`, translate back to same protein | #44, #69-72 |
| `no_false_negatives_scanner` | For known sites in databases, scanner always finds them | Axioms A1-A15 |
| `no_false_positives_scanner` | Scanner reports only real sites | Axioms A1-A15 |
| `valine_gt_mandatory` | All 4 Valine codons contain "GT" | #47, #52 |
| `gt_free_codons_exist_for_non_v` | Every non-V, non-Stop AA has a GT-free codon | #62 |
| `certificate_verification_idempotent` | `verify_certificate(cert.to_dict())` returns VERIFIED for valid certs | #54 |
| `dual_threshold_exclusive` | For any score, exactly one of PASS/UNCERTAIN/FAIL holds | #41 |
| `refinement_5to3_monotone` | `refine5to3(five_valued_and(a, b)) == three_valued_and(refine5to3(a), refine5to3(b))` | #1-12 |
| `restriction_site_no_compose` | Junction-crossing sites are detected | #33, #34 |

### 4.2 Which Theorems Have Corresponding Property Tests

| Theorem Category | Property Test Coverage | Priority |
|---|---|---|
| Three-valued logic (12 theorems) | High — 5-valued extension is straightforward | High |
| Pattern matching (3 theorems) | High — direct string correspondence | High |
| NDFST semantics (4 theorems) | Medium — requires NDFST test infrastructure | Medium |
| Scanner completeness/soundness (7+ axioms) | High — cross-validation against databases | Critical |
| Type system soundness (32 predicates) | High — each predicate needs independent verification | Critical |
| Compositional soundness (4 theorems) | High — `combined_verdict` is the key function | High |
| Splicing resolution (8 theorems) | Medium — requires splice site test data | Medium |
| Mutagenesis (15+ theorems) | High — Valine/GT analysis is directly testable | High |
| Certificate soundness (1 theorem) | High — verification must be independently correct | Critical |
| SLOT independence (6 theorems) | Medium — needs formal_mode implementation | Medium |
| Pipeline preservation (4 theorems) | High — end-to-end testing | High |

---

## 5. Summary Statistics

| Metric | Count |
|---|---|
| Total Lean4 theorems mapped | 76 |
| Total Lean4 axioms (proof parameters) | 20 |
| Total theorem-to-Python mappings | 76 |
| Critical refinement gaps | 7 |
| Property-based tests proposed | 13 |
| Proof modules analyzed | 11 |
| Library modules analyzed | 4 |
| Python modules analyzed | 8+ |

---

## 6. Conclusion

The Lean4 formal model and the Python implementation share the same architectural intent: a type system for biological sequence optimization that provides soundness guarantees. The refinement mapping reveals 7 categories of gaps, the most critical being:

1. **3-valued vs. 5-valued logic** — solvable via conservative abstraction mapping
2. **SLOT predicate evaluation** — the formal model says UNCERTAIN; Python evaluates them — solvable via `formal_mode` flag
3. **Float vs. Rat** — solvable via epsilon tolerance
4. **Scanner/oracle correctness** — not formally verified but well-validated empirically

The 76 theorem-to-Python mappings provide a complete traceability chain from formal proof to running code. Each gap has a concrete closure strategy, and the proposed property-based tests would provide empirical evidence that the Python implementation refines the Lean4 model.
