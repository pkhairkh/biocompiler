# BioCompiler Formal Proof Structure

This directory contains the Lean4 formal proof of the BioCompiler type system
soundness. The proof establishes that "well-typed genes don't go wrong" — if
the type system returns PASS for all predicates, then all semantic properties
hold.

## Prerequisites

- [Lean4](https://leanprover.org/) (v4.30.0, as specified in `lean-toolchain`)
- [Lake](https://github.com/leanprover/lake) build system (bundled with Lean4)
- [Elan](https://github.com/leanprover/elan) (Lean version manager, recommended)

## Building

```bash
# Install elan (Lean version manager)
curl -sSfL https://github.com/leanprover/elan/releases/latest/download/elan-x86_64-unknown-linux-gnu.tar.gz | tar xz
./elan-init -y --default-toolchain none
source "$HOME/.elan/env"

# Build all proofs (from the proof/ directory)
cd proof/
lake build

# Fetch Mathlib cache (if applicable)
lake exe cache get || true
```

## Proof Architecture

```
Theorem 4: Compositional Soundness (Compositional.lean)
  evaluateAll = PASS → all hold
  + slot_predicates_uncertain: SLOT predicates never PASS
  + slot_predicates_dont_affect_pass: SLOT in list → no PASS
                    │
Theorem 3: Per-Predicate Soundness (TypeSystem.lean)
  evaluate P = PASS → holds(P)  (33 predicates: 13 core + 20 SLOT)
                    │
      ┌─────────────┼─────────────┐
      │             │             │
Theorem 1:     Theorem 2:     Theorem 5:
3-Valued       NDFST          SLOT-Indep.
PASS⊓PASS     Deterministic  Certs don't
FAIL sticky   computation    depend on FFI
```

## Module Overview

### Core Logic

| Module | Description |
|--------|-------------|
| `ThreeValued.lean` | Three-valued logic (PASS/UNCERTAIN/FAIL) with proved algebraic properties |
| `FiveValued.lean` | Five-valued logic extension (PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL) with refinement proof |
| `Sequence.lean` | Nucleotide sequences, pattern matching with completeness and soundness |
| `NDFST.lean` | Non-deterministic finite-state transducers for splicing |
| `CellularContext` | Context with cell type, ESE/ESS/ISE/ISS thresholds (in NDFST.lean) |

### Scanner Layer

| Module | Description |
|--------|-------------|
| `Scanners.lean` | Abstract scanner interfaces + concrete implementations |
| `ScannerProofs.lean` | Completeness/soundness proofs for CpG, Promoter, TM scanners |
| `OracleProofs.lean` | Proofs for mRNA, co-translational folding, splicing, CAI oracles |

### Type System

| Module | Description |
|--------|-------------|
| `TypeSystem.lean` | 33 type predicates (13 core + 20 SLOT), evaluation, and per-predicate soundness |
| `Compositional.lean` | Compositional soundness: overall PASS implies all properties hold |
| `Certificates.lean` | Guarantee certificate structure and certificate soundness |
| `SLOTIndependence.lean` | SLOT independence: certificates don't depend on FFI output |

### SLOT Verification

| Module | Description |
|--------|-------------|
| `SLOTVerification.lean` | SLOT predicate verification conditions, modes (conservative/verified/permissive) |
| `Refinement.lean` | VERIFIED mode refines CONSERVATIVE mode; simulation theorem |

### Specialized Proofs

| Module | Description |
|--------|-------------|
| `SplicingResolution.lean` | Splice site resolution and canonical donor/acceptor proofs |
| `Mutagenesis.lean` | Synonymous mutation theorems, GT/AG analysis, codon degeneracy |

### Top-Level Entry Point

| Module | Description |
|--------|-------------|
| `Soundness.lean` | Main theorem module — imports all components, re-exports theorems |

## Proof Status

All modules are **FULLY PROVED** with 0 sorry. The Trusted Computing Base (TCB)
consists of 3 scanner axioms + 4 SLOT verification axioms:

**Scanner axioms (3):**
1. `SpliceSiteScanner.scanner_completeness` — scanner finds all cryptic sites
2. `SpliceSiteScanner.scanner_soundness` — scanner only reports real sites
3. `SpliceSiteScanner.borderline_completeness` — borderline scanner completeness

**SLOT verification axioms (4, added in Task 1.10):**
4. `vc_imply_no_unexpected_tm` — TMHMM VCs imply no unexpected TM domain
5. `vc_imply_stable_folding` — FoldX VCs imply stable folding
6. `vc_imply_soluble_expression` — ProteinSol VCs imply adequate solubility
7. `vc_imply_low_immunogenicity` — NetMHC VCs imply low immunogenicity

These are **parameters of the proof**, not gaps: the soundness theorem says
"ASSUMING the remaining scanners and tools are correct, the type system is
sound." This follows the standard approach in formal methods.

## Five-Valued Logic Extension

The `FiveValued.lean` module extends the three-valued logic with intermediate
verdicts:

```
Ordering: PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL
```

Key theorems:
- `and_project_refines`: `project (five_valued_and a b) = three_valued_and (project a) (project b)`
- `and_pass_pass`: `five_valued_and PASS PASS = PASS`
- `and_fail_absorb`: `five_valued_and FAIL x = FAIL` for all x
- `and_comm`, `and_assoc`: Commutativity and associativity

The projection maps LIKELY_PASS and LIKELY_FAIL to UNCERTAIN, making the
5-valued logic a conservative refinement of the 3-valued logic.

## SLOT Property Semantics

Four SLOT predicates have non-vacuous semantic definitions:

| Predicate | Property Semantics |
|-----------|-------------------|
| `NoUnexpectedTMDomain` | ∀ pos, tmHydrophobicFraction < threshold |
| `StableFolding` | predictedStabilityScore ≤ ddgThreshold |
| `SolubleExpression` | camSolScore ≥ minScore |
| `LowImmunogenicity` | ∀ pos, mhcBindingAffinity > maxScore |

The remaining 16 SLOT predicates use vacuous `True` semantics, awaiting
progressive strengthening.

## CI Integration

The `.github/workflows/ci.yml` includes a `proof-check` job that:
1. Installs elan and Lean4 (v4.30.0)
2. Runs `lake build` in the `proof/` directory
3. Checks for `sorry` in proof files (must be sorry-free)
4. Is a separate job from Python tests and can be skipped if Lean is not needed

## File Structure

```
proof/
├── lakefile.lean          # Lake build configuration
├── lean-toolchain         # Lean4 version (leanprover/lean4:v4.30.0)
├── lake-manifest.json     # Dependency manifest
├── BioCompiler/
│   ├── Soundness.lean     # Main theorem module (entry point)
│   ├── ThreeValued.lean   # 3-valued logic
│   ├── FiveValued.lean    # 5-valued logic extension
│   ├── Sequence.lean      # Nucleotide sequences
│   ├── NDFST.lean         # Non-deterministic FSTs
│   ├── Scanners.lean      # Scanner interfaces
│   ├── ScannerProofs.lean # Scanner completeness/soundness proofs
│   ├── OracleProofs.lean  # Oracle proofs
│   ├── TypeSystem.lean    # Type predicates and soundness
│   ├── Compositional.lean # Compositional soundness
│   ├── Certificates.lean  # Guarantee certificates
│   ├── SLOTIndependence.lean  # SLOT independence theorem
│   ├── SLOTVerification.lean  # SLOT verification conditions
│   ├── Refinement.lean    # Mode refinement theorem
│   ├── SplicingResolution.lean # Splice resolution proofs
│   └── Mutagenesis.lean   # Mutation analysis proofs
└── doc/
    └── DOC-11-Formal-Soundness-Proof.md
```
