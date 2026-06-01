/-
  BioCompiler.Certificates — Certificate Validity and Soundness

  This module defines certificate validity and proves that a valid
  guarantee certificate implies all claimed properties hold.

  Reference: DOC-03 (SDD) §3.5.5, DOC-01 (SRS) INV-TYP-01
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners
import BioCompiler.TypeSystem
import BioCompiler.Compositional

namespace BioCompiler

open Verdict Sequence

-- ==============================================================================
-- Certificate Validity
-- ==============================================================================

/-- A guarantee certificate is valid if and only if the composed evaluation
    of all claimed predicates yields PASS. This is the formal definition that
    ties the type system to the certificate mechanism. -/
def certificateValid [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) : Prop :=
  @evaluateAll inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst predicates seq ctx = PASS

-- ==============================================================================
-- Certificate Soundness Theorem
-- ==============================================================================

/-- THEOREM (Certificate Soundness): A valid guarantee certificate implies
    all claimed properties hold.

    Proof: certificateValid ↔ evaluateAll = PASS (by definition),
    then compositional_soundness gives ∀ P ∈ predicates, propertyHolds P. -/
theorem certificate_soundness [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    @certificateValid inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst predicates seq ctx →
    ∀ P ∈ predicates, @propertyHolds inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx := by
  intro h_cert P hP
  exact @compositional_soundness inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst
    predicates seq ctx h_cert P hP

end BioCompiler
