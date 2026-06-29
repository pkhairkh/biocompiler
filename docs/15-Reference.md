# DOC-15: Technical Reference

**Document ID:** DOC-15
**Version:** 0.9.3
**Date:** 2026-06-07
**Classification:** Reference Documentation

---

This document consolidates detailed technical reference material previously in the project README. For architecture, design rationale, and formal methods, see the [docs/](00-README.md) document map.

---

## 1. Type Predicate System

BioCompiler features a type predicate system organized into five domains. The Lean4 soundness proof in `TypeSystem.lean` defines **36 type predicates** (17 core + 19 SLOT-dependent). The Python implementation includes 7 additional extended diagnostic predicates for a total of **43 type predicates**, evaluating a context-dependent subset by default (typically 12 for short sequences, up to 43 for full analysis with all engines available).

The core/SLOT classification follows the Lean4 formal model: core predicates evaluate deterministically from the sequence and context alone, while SLOT-dependent predicates require external tool output (FFI) and return UNCERTAIN in conservative mode.

### DNA-Level Predicates (13)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 1 | `SpliceCorrect` | Splice isoform correctness for cell type | Core |
| 2 | `NoCrypticSplice` | No splice-site-like motifs exceeding MaxEntScan threshold | Core |
| 3 | `CodonAdapted` | CAI ≥ threshold for target organism | Core |
| 4 | `GCInRange` | Global GC% within [lo, hi] bounds | Core |
| 5 | `NoRestrictionSite` | No enzyme recognition sites present | Core |
| 6 | `InFrame` | Reading frame and exon boundary integrity | Core |
| 7 | `NoInstabilityMotif` | No known instability motifs (e.g., AT-rich repeats) | Core |
| 8 | `NoCpGIsland` | No CpG islands (sliding window GC + Obs/Exp CG ratio) | Core |
| 9 | `NoGTDinucleotide` | No avoidable GT dinucleotides (cross-codon aware) | Core |
| 10 | `NoStopCodons` | No internal stop codons | Core |
| 11 | `ValidCodingSeq` | In-frame, valid codons only | Core |
| 12 | `CodonOptimality` | Geometric mean CAI ≥ threshold | Core |
| 13 | `NoCrypticPromoter` | Cryptic promoter avoidance | Core |
| 14 | `ConservationScore` | BLOSUM62-based AA conservation | SLOT-dependent |
| 15 | `NoUnexpectedTMDomain` | Unexpected transmembrane domain detection | SLOT-dependent |
| 16 | `mRNASecondaryStructure` | mRNA secondary structure around RBS | SLOT-dependent |
| 17 | `CoTranslationalFolding` | Co-translational folding pause-site preservation | SLOT-dependent |

### Structure Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 18 | `StructureConfidence` | ESMFold structure quality confidence | SLOT-dependent |
| 19 | `NoMisfoldingRisk` | Misfolding risk indicators | SLOT-dependent |
| 20 | `CorrectFoldTopology` | Fold topology validation | SLOT-dependent |
| 21 | `NoUnexpectedInteraction` | Unwanted protein-protein interactions | SLOT-dependent |

### Stability Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 22 | `StableFolding` | Thermodynamic stability (ΔG) | SLOT-dependent |
| 23 | `NoDestabilizingMutation` | No high-ΔΔG mutations | SLOT-dependent |
| 24 | `DisulfideBondIntegrity` | Cysteine pairing check | SLOT-dependent |
| 25 | `HydrophobicCoreQuality` | Hydrophobic core composition | SLOT-dependent |

### Solubility Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 26 | `SolubleExpression` | CamSol solubility score | SLOT-dependent |
| 27 | `NoAggregationProneRegion` | Aggregation-prone region detection | SLOT-dependent |
| 28 | `ChargeComposition` | Charge balance and pI | SLOT-dependent |
| 29 | `NoLongHydrophobicStretch` | Long hydrophobic stretch detection | SLOT-dependent |

### Immunogenicity Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 30 | `LowImmunogenicity` | Overall immunogenicity score | SLOT-dependent |
| 31 | `NoStrongTCellEpitope` | MHC binding epitope detection | SLOT-dependent |
| 32 | `NoDominantBCellEpitope` | B-cell epitope coverage | SLOT-dependent |
| 33 | `PopulationCoverageSafe` | MHC allele population coverage | SLOT-dependent |

> **Note on SLOT-dependent predicates:** SLOT-dependent predicates rely on external tool output (FFI) and cannot produce PASS in the formal model. They always return UNCERTAIN, which is the correct and safe behavior — the system declines to guarantee what it cannot verify deterministically. The Lean4 proof formally establishes that SLOT-dependent predicates never produce PASS (`ffi_never_pass`), so they cannot compromise compositional soundness. Certificate validity is independent of SLOT values. See [DOC-14: SLOT Predicate Proof-Implementation Gap](14-SLOT-Proof-Implementation-Gap.md).

---

## 2. SLOT Architecture

The 36-predicate type system is partitioned into **core** and **SLOT-dependent** predicates based on whether evaluation requires external tool output (FFI).

**Core predicates** (14) evaluate deterministically from the nucleotide sequence and grammar rules alone. They produce PASS or FAIL verdicts and are individually proved sound in Lean4.

**SLOT-dependent predicates** (22) require external tool output — structure prediction (AlphaFold/ESMFold), stability calculation (FoldX), solubility scoring (CamSol), immunogenicity prediction (MHC binding), etc. In the formal model, these always return UNCERTAIN because the external tool is treated as a non-deterministic black box.

The Lean4 proof establishes three key properties:

1. **`ffi_never_pass`**: SLOT-dependent predicates never produce PASS, regardless of SLOT values
2. **`slot_predicates_dont_affect_pass`**: If any SLOT-dependent predicate is in the evaluation list, `evaluateAll` cannot return PASS
3. **`certificate_slot_independent`**: Certificate validity does not change when SLOT values are modified

This means the soundness guarantee is preserved even when the type system includes SLOT-dependent predicates — they are safely isolated from the PASS/FAIL logic. For details on the three SLOT modes (CONSERVATIVE, VERIFIED, PERMISSIVE) and the proof-implementation gap, see [DOC-14](14-SLOT-Proof-Implementation-Gap.md).

---

## 3. Integrated Optimizer API

All six analysis engines (ESMFold, FoldX, CamSol, Immunogenicity, Deimmunization, Protein Design) share a unified result type hierarchy:

- **`BaseEngineResult`** — base class with `sequence`, `primary_score`, `classification`, `engine_name`, `primary_score_label`, `passed`, `success`, `execution_time_s`
- **`MutationResult`** — unified mutation representation with `position`, `original`, `mutant`, `delta_score`, `score_type`, `engine`, `recommendation`
- **`BatchResult[T]`** — batch result container with `results`, `errors`, `successful`, `failed`, `total_time_s`
- **`EngineTimer`** — context manager for execution timing
- **`EngineConfig`** — shared configuration (`use_cache`, `timeout_s`, `verbose`, `max_workers`)
- **`classify_score()`** — threshold-based score classification shared across engines

Each engine result type inherits from `BaseEngineResult` and provides backward-compatible property aliases (e.g., `FoldXResult.ddg` → `primary_score`, `ESMFoldResult.plddt` → `primary_score`, `CamSolResult.score` → `primary_score`).

```python
from biocompiler.engines.base import BaseEngineResult, MutationResult, BatchResult, EngineTimer

# All engine results share the unified interface
result = empirical_stability("MVHLTPEEK...")  # Returns FoldXResult(BaseEngineResult)
print(result.primary_score)    # ΔΔG value
print(result.classification)   # "stable", "unstable", etc.
print(result.passed)           # True if success
print(result.engine_name)      # "foldx"

# Batch operations return BatchResult
batch = predict_structure_batch(sequences)  # Returns BatchResult[ESMFoldResult]
print(f"{batch.successful}/{batch.total} succeeded")
```

---

## 4. Type-Directed Protein Mutagenesis

The key innovation: the type predicate does not just **verify** — it **directs design** across the central dogma boundary (DNA→RNA→Protein).

When the type system proves that NO codon assignment can satisfy all predicates (e.g., Valine's codons ALL contain GT, making cryptic splice donor elimination impossible), the mutagenesis engine proposes conservative amino acid substitutions ranked by BLOSUM62 score.

**Key improvement**: The mutagenesis engine distinguishes between GT-mandatory positions (Valine only — all codons contain GT) and optimizer weaknesses (GT-free codons exist but were not used). Mutagenesis is only proposed for GT-mandatory positions, preventing unnecessary protein modifications and exposing optimizer bugs for repair.

**HBB proof of concept**: 15 V→I substitutions (BLOSUM62=+3 each) turn an impossible constraint (5/6 predicates failing) into a solvable one, at only 0.2% CAI cost and 99.3% protein identity.

```
Optimizer → Type System → [FAIL: NoCrypticSplice at V positions]
                           |
                    Mutagenesis Engine
                    (V→I, BLOSUM62=+3, GT-free codons)
                           |
                    Modified Protein → Optimizer → Type System → [PASS]
```

See also: [ADR-0013: Mutagenesis GT-Mandatory Distinction](adr/ADR-0013-mutagenesis-gt-mandatory.md), [DOC-06: Design Rationale](06-Design-Rationale.md).

---

## 5. Trusted Computing Base (TCB)

The soundness proof rests on 3 class-field axioms in the SpliceSiteScanner type class (15 explicit `axiom` declarations (tool-soundness contracts in SLOTVerification.lean)). These are field axioms within the class — not standalone `axiom` declarations using Lean4's `axiom` keyword, but unproved assumptions that any SpliceSiteScanner instance must satisfy. These are the boundaries beyond which the proof does not extend:

| TCB | Assumption | Rationale |
|-----|-----------|-----------|
| TCB-1 | `SpliceSiteScanner.scanner_completeness`: the scanner finds all cryptic splice sites with score ≥ crypticThreshold | Adversarial testing; comparison with MaxEntScan reference scores |
| TCB-2 | `SpliceSiteScanner.scanner_soundness`: the scanner only reports real splice sites | Validation against GENCODE-annotated splice sites |
| TCB-3 | `SpliceSiteScanner.borderline_completeness`: the borderline scanner finds all sites with score in [uncertainLoThreshold, crypticThreshold) | Validation against borderline splice site predictions |

Every guarantee is conditional on these assumptions. The proof is honest about where the modeling ends and the biology begins.

---

## 6. Honest Limitations

1. **UNCERTAIN is the common case for complex genes.** For genes with many cryptic splice sites or GC content near boundaries, the type system may return UNCERTAIN. This is correct behavior — the system declines to guarantee what it cannot verify.

2. **The proof guarantees correctness of the model, not of biology.** If the NDFST grammar is wrong (TCB-1), the guarantee is about a wrong model. Wet-lab validation is needed to close this gap.

3. **No quantitative predictions.** The framework cannot answer "what fraction of transcripts will include exon 5?" For those, probability is unavoidable — but the system explicitly declines to answer rather than giving misleading deterministic answers.

4. **Grammar curation requires domain expertise.** The YAML grammar configuration (ADR-0007) must be curated by someone who understands the splicing biology of the target gene and cell type.

5. **Mutagenesis changes the protein.** Type-directed mutagenesis proposes conservative substitutions (BLOSUM62 ≥ 0 by default), but these are not guaranteed to preserve protein function. The biologist must evaluate whether each substitution is acceptable for their application.

---

## 7. Repository Structure

```
biocompiler/
├── proof/                        # Lean4 machine-verified soundness proof
│   ├── lakefile.lean             # Build configuration
│   └── BioCompiler/              # Proof modules (Soundness, TypeSystem, ThreeValued, etc.)
├── src/biocompiler/              # Production Python package
│   ├── types.py                  # Core data structures (Verdict, Token, Certificate)
│   ├── scanner.py                # Multi-DFA motif detection
│   ├── splicing.py               # NDFST isoform computation
│   ├── translation.py            # Codon-to-amino-acid FST + CAI computation
│   ├── type_system/              # Predicate registry + evaluator functions (decomposed)
│   ├── optimizer/                 # Greedy multi-phase optimizer + mutagenesis loop (decomposed)
│   ├── mutagenesis.py            # Type-directed mutagenesis
│   ├── certificate.py            # Graduated certificate generation + verification
│   ├── engine_base.py            # Unified engine API (BaseEngineResult, MutationResult)
│   ├── organisms/                # Organism-specific data (25 organisms + tAI)
│   ├── biosecurity/              # Biosecurity screening (v12)
│   ├── tai.py                    # tRNA Adaptation Index (v12)
│   ├── api.py                    # REST API (FastAPI)
│   └── cli.py                    # Command-line interface
├── tests/                        # Test suite (14,600+ tests)
├── docs/                         # Complete SE specification (19 docs + 18 ADRs)
└── paper/                        # LaTeX manuscript
```
