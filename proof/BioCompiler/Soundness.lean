/-
  BioCompiler.Soundness — Main Theorem Module

  This is the top-level module that imports all components and states
  the complete soundness theorem for the BioCompiler type system.

  ╔══════════════════════════════════════════════════════════════════════════════╗
  ║                                                                            ║
  ║   THEOREM (BioCompiler Type System Soundness):                             ║
  ║                                                                            ║
  ║   ∀ (predicates : List TypePredicate)                                      ║
  ║     (seq : Sequence)                                                       ║
  ║     (ctx : CellularContext),                                               ║
  ║                                                                            ║
  ║     evaluateAll predicates seq ctx = Verdict.PASS                           ║
  ║       → ∀ P ∈ predicates, propertyHolds P seq ctx                         ║
  ║                                                                            ║
  ║   "Well-typed genes don't go wrong."                                       ║
  ║                                                                            ║
  ║   Corollary (SLOT-Independence):                                           ║
  ║     The soundness guarantee is independent of FFI output.                  ║
  ║     Guarantee certificates remain valid regardless of external tool         ║
  ║     behavior.                                                              ║
  ║                                                                            ║
  ╚══════════════════════════════════════════════════════════════════════════════╝

  Proof Architecture:
  ┌─────────────────────────────────────────────────────────────────────────┐
  │  Theorem 4: Compositional Soundness (Compositional.lean)              │
  │  evaluateAll = PASS → all hold                                        │
  └───────────────────────────┬─────────────────────────────────────────────┘
                              │
  ┌───────────────────────────┴─────────────────────────────────────────────┐
  │  Theorem 3: Per-Predicate Soundness (TypeSystem.lean)                 │
  │  evaluate P = PASS → holds(P)                                         │
  └───────────────────────────┬─────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
  ┌─────┴──────┐  ┌──────────┴──────┐  ┌──────────┴──────┐
  │ Theorem 1: │  │ Theorem 2:      │  │ Theorem 5:      │
  │ 3-Valued   │  │ NDFST Semantics │  │ SLOT-Indep.     │
  │ PASS⊓PASS │  │ Deterministic   │  │ Certs don't     │
  │ FAIL sticky│  │ computation     │  │ depend on FFI   │
  └────────────┘  └─────────────────┘  └─────────────────┘

  Proof Status:
  ┌─────────────────────────┬────────────┬────────────────────────────────────────┐
  │ Module                  │ Status     │ Notes                                  │
  ├─────────────────────────┼────────────┼────────────────────────────────────────┤
  │ ThreeValued.lean        │ FULLY      │ 12 theorems, 0 sorry, 0 axioms        │
  │ Sequence.lean           │ FULLY      │ Pattern matching with proved           │
  │                         │            │ completeness and soundness             │
  │ NDFST.lean              │ FULLY      │ ndfstRun_sound + ndfstRun_complete     │
  │                         │            │ both proved; ConsumesInput replaces    │
  │                         │            │ ValidPath for completeness             │
  │ Scanners.lean           │ FULLY      │ Concrete scanner implementations       │
  │                         │            │ with proved completeness; includes     │
  │                         │            │ dual-threshold borderline scanner      │
  │ TypeSystem.lean         │ FULLY      │ All 8 predicates proved, 0 sorry       │
  │                         │            │ (incl. NoCpGIsland); NoCrypticSplice   │
  │                         │            │ uses dual-threshold PASS/UNCERTAIN/FAIL│
  │ Compositional.lean      │ FULLY      │ UNCERTAIN propagation proved via       │
  │                         │            │ foldl_uncertain_ne_pass lemma          │
  │ SLOTIndependence.lean   │ FULLY      │ All theorems proved, 0 sorry           │
  └─────────────────────────┴────────────┴────────────────────────────────────────┘

  Trusted Computing Base (axioms that are NOT proved within Lean4):
  1. SpliceSiteScanner.scanner_completeness — scanner finds all cryptic sites
  2. SpliceSiteScanner.scanner_soundness — scanner only reports real sites
  3. SpliceSiteScanner.borderline_completeness — borderline scanner finds all
     sites with score in [uncertainLoThreshold, crypticThreshold)
  4. CpGIslandScanner.scanner_completeness — CpG island scanner finds all islands
  5. CpGIslandScanner.scanner_soundness — CpG island scanner only reports real islands
  6. SplicingNDFST.output_is_valid — NDFST outputs are valid isoforms
  7. SplicingNDFST.all_isoforms_produced — NDFST is complete
  8. CodonAdaptationIndex.computeCAI — deterministic CAI computation

  These are PARAMETERS of the proof, not gaps. The soundness theorem says:
  "ASSUMING the scanners are complete and the NDFST is correct,
   the type system is sound." This is the standard approach in formal methods:
  prove the TYPE SYSTEM sound conditional on correct oracles, then validate
  the oracles independently.

  DUAL-THRESHOLD BACKWARD COMPATIBILITY:
  The dual-threshold NoCrypticSplice (PASS/UNCERTAIN/FAIL) is backward
  compatible with the original binary (PASS/FAIL) behavior. If no borderline
  scanner is provided, the old behavior is recovered by treating all sites
  with score in [uncertainLoThreshold, crypticThreshold) as PASS (since they
  are below the cryptic threshold). The PASS guarantee is now stronger:
  all sites must have score < uncertainLoThreshold, not just < crypticThreshold.
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners
import BioCompiler.TypeSystem
import BioCompiler.Compositional
import BioCompiler.SLOTIndependence

namespace BioCompiler

-- Re-export the main theorems for easy access

/-- The central soundness theorem: PASS implies property holds. -/
def soundness := @type_soundness

/-- Compositional soundness: overall PASS implies all properties hold. -/
def comp_soundness := @compositional_soundness

/-- Certificate soundness: valid certificate implies all properties hold. -/
def cert_soundness := @certificate_soundness

/-- SLOT independence: certificate validity is FFI-output independent. -/
def slot_indep := @full_slot_independence

/-- Three-valued conjunction preserves PASS. -/
def and_sound := @and_pass_pass

/-- FFI-dependent predicates never produce PASS. -/
def ffi_no_pass := @ffi_never_pass

end BioCompiler
