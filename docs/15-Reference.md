# DOC-15: Technical Reference

**Document ID:** DOC-15
**Version:** 1.0.0
**Date:** 2026-03-04
**Classification:** Reference Documentation

---

This document consolidates detailed technical reference material previously in the project README. For architecture, design rationale, and formal methods, see the [docs/](00-README.md) document map.

---

## 1. 28-Predicate Type System (v9.0.0)

BioCompiler v9.0.0 extends the type system from 8 to 28 predicates, organized into five domains. The Lean4 soundness proof covers all 28 predicates: 13 core predicates produce PASS/FAIL verdicts (individually proved sound), and 19 SLOT-dependent predicates always return UNCERTAIN (vacuously sound — they never produce PASS, so the implication `evaluate(P) = PASS → propertyHolds(P)` holds trivially).

### DNA-Level Predicates (12)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 1 | `NoStopCodons` | No internal stop codons | Core |
| 2 | `NoCrypticSplice` | No splice-site-like motifs exceeding MaxEntScan threshold | Core |
| 3 | `NoCpGIsland` | No CpG islands (sliding window GC + Obs/Exp CG ratio) | Core |
| 4 | `NoRestrictionSite` | No enzyme recognition sites present | Core |
| 5 | `NoGTDinucleotide` | No avoidable GT dinucleotides (cross-codon aware) | Core |
| 6 | `ValidCodingSeq` | In-frame, valid codons only | Core |
| 7 | `ConservationScore` | BLOSUM62-based AA conservation | Core |
| 8 | `CodonOptimality` | Geometric mean CAI ≥ threshold | Core |
| 9 | `NoCrypticPromoter` | Cryptic promoter avoidance | SLOT-dependent |
| 10 | `NoUnexpectedTMDomain` | Unexpected transmembrane domain detection | SLOT-dependent |
| 11 | `mRNASecondaryStructure` | mRNA secondary structure around RBS | SLOT-dependent |
| 12 | `CoTranslationalFolding` | Co-translational folding pause-site preservation | SLOT-dependent |

### Structure Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 13 | `StructureConfidence` | ESMFold structure quality confidence | SLOT-dependent |
| 14 | `NoMisfoldingRisk` | Misfolding risk indicators | SLOT-dependent |
| 15 | `CorrectFoldTopology` | Fold topology validation | SLOT-dependent |
| 16 | `NoUnexpectedInteraction` | Unwanted protein-protein interactions | SLOT-dependent |

### Stability Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 17 | `StableFolding` | Thermodynamic stability (ΔG) | SLOT-dependent |
| 18 | `NoDestabilizingMutation` | No high-ΔΔG mutations | SLOT-dependent |
| 19 | `DisulfideBondIntegrity` | Cysteine pairing check | SLOT-dependent |
| 20 | `HydrophobicCoreQuality` | Hydrophobic core composition | SLOT-dependent |

### Solubility Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 21 | `SolubleExpression` | CamSol solubility score | SLOT-dependent |
| 22 | `NoAggregationProneRegion` | Aggregation-prone region detection | SLOT-dependent |
| 23 | `ChargeComposition` | Charge balance and pI | Core |
| 24 | `NoLongHydrophobicStretch` | Long hydrophobic stretch detection | Core |

### Immunogenicity Predicates (4)

| # | Predicate | Checks | Category |
|---|-----------|--------|----------|
| 25 | `LowImmunogenicity` | Overall immunogenicity score | SLOT-dependent |
| 26 | `NoStrongTCellEpitope` | MHC binding epitope detection | SLOT-dependent |
| 27 | `NoDominantBCellEpitope` | B-cell epitope coverage | SLOT-dependent |
| 28 | `PopulationCoverageSafe` | MHC allele population coverage | SLOT-dependent |

> **Note on SLOT-dependent predicates:** SLOT-dependent predicates rely on external tool output (FFI) and cannot produce PASS in the formal model. They always return UNCERTAIN, which is the correct and safe behavior — the system declines to guarantee what it cannot verify deterministically. The Lean4 proof formally establishes that SLOT-dependent predicates never produce PASS (`ffi_never_pass`), so they cannot compromise compositional soundness. Certificate validity is independent of SLOT values. See [DOC-14: SLOT Predicate Proof-Implementation Gap](14-SLOT-Proof-Implementation-Gap.md).

---

## 2. SLOT Architecture

The 28-predicate type system is partitioned into **core** and **SLOT-dependent** predicates based on whether evaluation requires external tool output (FFI).

**Core predicates** (13) evaluate deterministically from the nucleotide sequence and grammar rules alone. They produce PASS or FAIL verdicts and are individually proved sound in Lean4.

**SLOT-dependent predicates** (19) require external tool output — structure prediction (AlphaFold/ESMFold), stability calculation (FoldX), solubility scoring (CamSol), immunogenicity prediction (MHC binding), etc. In the formal model, these always return UNCERTAIN because the external tool is treated as a non-deterministic black box.

The Lean4 proof establishes three key properties:

1. **`ffi_never_pass`**: SLOT-dependent predicates never produce PASS, regardless of SLOT values
2. **`slot_predicates_dont_affect_pass`**: If any SLOT-dependent predicate is in the evaluation list, `evaluateAll` cannot return PASS
3. **`certificate_slot_independent`**: Certificate validity does not change when SLOT values are modified

This means the soundness guarantee is preserved even when the type system includes SLOT-dependent predicates — they are safely isolated from the PASS/FAIL logic. For details on the three SLOT modes (CONSERVATIVE, VERIFIED, PERMISSIVE) and the proof-implementation gap, see [DOC-14](14-SLOT-Proof-Implementation-Gap.md).

---

## 3. Unified Engine API (v9.0.0)

All six analysis engines (ESMFold, FoldX, CamSol, Immunogenicity, Deimmunization, Protein Design) share a unified result type hierarchy:

- **`BaseEngineResult`** — base class with `sequence`, `primary_score`, `classification`, `engine_name`, `primary_score_label`, `passed`, `success`, `execution_time_s`
- **`MutationResult`** — unified mutation representation with `position`, `original`, `mutant`, `delta_score`, `score_type`, `engine`, `recommendation`
- **`BatchResult[T]`** — batch result container with `results`, `errors`, `successful`, `failed`, `total_time_s`
- **`EngineTimer`** — context manager for execution timing
- **`EngineConfig`** — shared configuration (`use_cache`, `timeout_s`, `verbose`, `max_workers`)
- **`classify_score()`** — threshold-based score classification shared across engines

Each engine result type inherits from `BaseEngineResult` and provides backward-compatible property aliases (e.g., `FoldXResult.ddg` → `primary_score`, `ESMFoldResult.plddt` → `primary_score`, `CamSolResult.score` → `primary_score`).

```python
from biocompiler.engine_base import BaseEngineResult, MutationResult, BatchResult, EngineTimer

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

The key innovation: the type predicate doesn't just **verify** — it **directs design** across the central dogma boundary (DNA→RNA→Protein).

When the type system proves that NO codon assignment can satisfy all predicates (e.g., Valine's codons ALL contain GT, making cryptic splice donor elimination impossible), the mutagenesis engine proposes conservative amino acid substitutions ranked by BLOSUM62 score.

**v7.1 improvement**: The mutagenesis engine distinguishes between GT-mandatory positions (Valine only — all codons contain GT) and optimizer weaknesses (GT-free codons exist but weren't used). Mutagenesis is only proposed for GT-mandatory positions, preventing unnecessary protein modifications and exposing optimizer bugs for repair.

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

The soundness proof rests on 5 explicitly stated assumptions. These are the boundaries beyond which the proof does not extend:

| TCB | Assumption | Rationale |
|-----|-----------|-----------|
| TCB-1 | The NDFST grammar faithfully models the splicing process for the target cell type | Grammar curation from experimental data (ENCODE/GTEx) |
| TCB-2 | The codon usage table accurately reflects translation efficiency in the target organism | Standard biochemistry data (Kazusa) |
| TCB-3 | The scanner thresholds are conservative upper bounds for non-functional sites | Domain-expert curation |
| TCB-4 | The restriction enzyme recognition sequences are correct | REBASE database |
| TCB-5 | Pattern matching is implemented correctly | Standard algorithm; could itself be verified |

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
│   ├── type_system.py            # Predicate registry + evaluator functions
│   ├── optimization.py           # Greedy multi-phase optimizer + mutagenesis loop
│   ├── mutagenesis.py            # Type-directed mutagenesis
│   ├── certificate.py            # Graduated certificate generation + verification
│   ├── engine_base.py            # Unified engine API (BaseEngineResult, MutationResult)
│   ├── organisms/                # Organism-specific data (5 organisms)
│   ├── api.py                    # REST API (FastAPI)
│   └── cli.py                    # Command-line interface
├── tests/                        # Test suite (420+ tests)
├── docs/                         # Complete SE specification (14 docs + ADRs)
└── paper/                        # LaTeX manuscript
```
