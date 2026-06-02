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
  ║   Corollary (SLOT Predicates Don't Affect PASS):                           ║
  ║     If any SLOT-dependent predicate is in the list, evaluateAll            ║
  ║     cannot return PASS. Only core predicates contribute PASS.              ║
  ║                                                                            ║
  ╚══════════════════════════════════════════════════════════════════════════════╝

  Proof Architecture:
  ┌─────────────────────────────────────────────────────────────────────────┐
  │  Theorem 4: Compositional Soundness (Compositional.lean)              │
  │  evaluateAll = PASS → all hold                                        │
  │  + slot_predicates_uncertain: SLOT predicates never PASS              │
  │  + slot_predicates_dont_affect_pass: SLOT in list → no PASS           │
  └───────────────────────────┬─────────────────────────────────────────────┘
                              │
  ┌───────────────────────────┴─────────────────────────────────────────────┐
  │  Theorem 3: Per-Predicate Soundness (TypeSystem.lean)                 │
  │  evaluate P = PASS → holds(P)  (32 predicates: 16 core + 16 SLOT)    │
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
  │                         │            │ dual-threshold borderline scanner;     │
  │                         │            │ PromoterScanner, TMDomainScanner,      │
  │                         │            │ mRNAStructureOracle, CoTranslational-  │
  │                         │            │ FoldingOracle added                    │
  │ TypeSystem.lean         │ FULLY      │ All 32 predicates proved, 0 sorry      │
  │                         │            │ (16 core + 16 SLOT-dependent);         │
  │                         │            │ NoStopCodons, NoGTDinucleotide,        │
  │                         │            │ ValidCodingSeq, ConservationScore,     │
  │                         │            │ CodonOptimality, NoCrypticPromoter,    │
  │                         │            │ NoUnexpectedTMDomain,                  │
  │                         │            │ mRNASecondaryStructure,                │
  │                         │            │ CoTranslationalFolding added;          │
  │                         │            │ isSLOT classification function;        │
  │                         │            │ SLOT predicates have vacuously true     │
  │                         │            │ soundness                              │
  │ Compositional.lean      │ FULLY      │ UNCERTAIN propagation proved;          │
  │                         │            │ slot_predicates_uncertain lemma added; │
  │                         │            │ slot_predicates_dont_affect_pass;      │
  │                         │            │ all_core_if_pass corollary             │
  │ SLOTIndependence.lean   │ FULLY      │ All theorems updated for 32 predicates │
  │                         │            │ isCorePredicate extended to 32;        │
  │                         │            │ isCore_iff_not_slot theorem added      │
  └─────────────────────────┴────────────┴────────────────────────────────────────┘

  Trusted Computing Base (axioms that are NOT proved within Lean4):
  1. SpliceSiteScanner.scanner_completeness — scanner finds all cryptic sites
  2. SpliceSiteScanner.scanner_soundness — scanner only reports real sites
  3. SpliceSiteScanner.borderline_completeness — borderline scanner finds all
     sites with score in [uncertainLoThreshold, crypticThreshold)
  4. CpGIslandScanner.scanner_completeness — CpG island scanner finds all islands
  5. CpGIslandScanner.scanner_soundness — CpG island scanner only reports real islands
  6. PromoterScanner.scanner_completeness — promoter scanner finds all cryptic promoters
  7. PromoterScanner.scanner_soundness — promoter scanner only reports real promoters
  8. PromoterScanner.borderline_completeness — borderline promoter scanner finds all
     sites with score in [threshold*0.8, threshold)
  9. TMDomainScanner.scanner_completeness — TM domain scanner finds all domains
  10. TMDomainScanner.scanner_soundness — TM domain scanner only reports real domains
  11. TMDomainScanner.borderline_completeness — borderline TM scanner finds all
      borderline domains
  12. mRNAStructureOracle.oracle_completeness — mRNA structure oracle finds all
      strong structures
  13. mRNAStructureOracle.borderline_completeness — borderline mRNA structure oracle
  14. CoTranslationalFoldingOracle.oracle_completeness — co-translational folding
      oracle finds all disruptions
  15. CoTranslationalFoldingOracle.borderline_completeness — borderline folding oracle
  16. SplicingNDFST.output_is_valid — NDFST outputs are valid isoforms
  17. SplicingNDFST.all_isoforms_produced — NDFST is complete
  18. CodonAdaptationIndex.computeCAI — deterministic CAI computation

  These are PARAMETERS of the proof, not gaps. The soundness theorem says:
  "ASSUMING the scanners are complete and the NDFST is correct,
   the type system is sound." This is the standard approach in formal methods:
  prove the TYPE SYSTEM sound conditional on correct oracles, then validate
  the oracles independently.
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners
import BioCompiler.TypeSystem
import BioCompiler.Compositional
import BioCompiler.Certificates
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
def and_sound := @Verdict.and_pass_pass

/-- FFI-dependent predicates never produce PASS. -/
def ffi_no_pass := @ffi_never_pass

/-- SLOT-dependent type predicates never produce PASS. -/
def slot_no_pass := @slot_predicates_uncertain

/-- SLOT predicates in the list prevent overall PASS. -/
def slot_dont_affect := @slot_predicates_dont_affect_pass

end BioCompiler
