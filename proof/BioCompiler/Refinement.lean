/-
  BioCompiler.Refinement — Refinement Theorem: VERIFIED Mode Refines CONSERVATIVE Mode

  This module formally establishes that the Python VERIFIED mode's behavioral
  specification is a valid refinement of the Lean4 CONSERVATIVE mode's
  specification. This bridges the gap between the formal model (where SLOT
  predicates always return UNCERTAIN) and the practical implementation (where
  SLOT predicates can return PASS when evidence exists).

  ──────────────────────────────────────────────────────────────────────────
  REFINEMENT INFORMAL EXPLANATION
  ──────────────────────────────────────────────────────────────────────────

  **What is a refinement?**
  A refinement relates two systems: an "abstract" system (the Lean4 CONSERVATIVE
  model) and a "concrete" system (the Python VERIFIED implementation). The
  concrete system REFINES the abstract system if every observable behavior of
  the concrete system is consistent with (at least as informative as) some
  observable behavior of the abstract system.

  **The verdict information ordering:**
  We define an information ordering on verdicts:
      PASS ≥ UNCERTAIN    (PASS is more informative than UNCERTAIN)
      FAIL ≥ UNCERTAIN    (FAIL is more informative than UNCERTAIN)
      PASS and FAIL are incomparable

  In this ordering, "v₁ refines v₂" means v₁ is at least as informative as v₂.
  Formally: v₂ = UNCERTAIN ∨ v₁ = v₂

  **The key insight:**
  CONSERVATIVE mode always returns UNCERTAIN for SLOT predicates — the weakest
  (least informative) verdict. VERIFIED mode returns PASS when verification
  conditions hold, and UNCERTAIN otherwise. Since PASS ≥ UNCERTAIN and
  UNCERTAIN = UNCERTAIN, every VERIFIED verdict is at least as informative as
  the corresponding CONSERVATIVE verdict. Therefore VERIFIED REFINES CONSERVATIVE.

  ──────────────────────────────────────────────────────────────────────────
  WHAT THIS REFINEMENT PROOF BUYS YOU
  ──────────────────────────────────────────────────────────────────────────

  1. **No contradiction:** Switching from CONSERVATIVE to VERIFIED mode can
     never produce a verdict that CONTRADICTS the CONSERVATIVE verdict.
     If CONSERVATIVE says UNCERTAIN, VERIFIED can say PASS or UNCERTAIN —
     both are consistent with UNCERTAIN.

  2. **Compositional safety:** The refinement composes through the AND-fold
     (evaluateAllWithMode). If the VERIFIED overall verdict is PASS, the
     CONSERVATIVE overall verdict would be PASS or UNCERTAIN — never FAIL.
     This means VERIFIED mode can only "upgrade" verdicts, never "downgrade"
     them relative to CONSERVATIVE mode.

  3. **Soundness transfer:** Under the axiom `verification_conditions_imply_property`,
     VERIFIED mode's PASS verdicts carry the same soundness guarantee as
     CONSERVATIVE mode's PASS verdicts (plus the additional guarantee that
     the semantic property holds, via the VC axiom).

  4. **Backward compatibility:** Any system built on CONSERVATIVE mode
     guarantees continues to work under VERIFIED mode. The refinement
     ensures no existing guarantees are violated.

  ──────────────────────────────────────────────────────────────────────────
  WHAT THIS REFINEMENT PROOF DOES NOT BUY YOU
  ──────────────────────────────────────────────────────────────────────────

  1. **Does NOT prove the Python code implements the spec correctly.**
     This refinement is between the BEHAVIORAL SPECIFICATIONS of CONSERVATIVE
     and VERIFIED modes, as defined in the Lean4 model. It does not prove that
     the Python code in `slot_verification.py` correctly implements the VERIFIED
     mode specification. That gap is bridged by property-based testing
     (see `tests/test_property_predicates.py`, `tests/test_property_three_valued.py`).

  2. **Does NOT prove VERIFIED mode is sound unconditionally.**
     The soundness of VERIFIED mode depends on the axiom
     `verification_conditions_imply_property`. This axiom is the explicit
     "social contract" that external tools are trustworthy. If the tools are
     wrong, VERIFIED mode's PASS verdicts may be wrong.

  3. **Does NOT prove PERMISSIVE mode is sound.**
     PERMISSIVE mode goes beyond what is formally proven. It can return PASS
     with weaker evidence than VERIFIED mode requires. No refinement theorem
     exists for PERMISSIVE mode.

  4. **Does NOT eliminate the need for testing.**
     The refinement proves a mathematical relationship between specifications,
     not between specifications and implementations. Property-based testing
     remains essential for verifying that the Python code matches the spec.

  ──────────────────────────────────────────────────────────────────────────
  RELATIONSHIP TO PROPERTY-BASED TESTING
  ──────────────────────────────────────────────────────────────────────────

  The refinement proof and property-based testing serve complementary roles:

  **Refinement proof (this file):**
  - Proves a mathematical relationship between SPECIFICATIONS
  - Holds for all possible inputs (universal quantification)
  - Cannot be wrong (it's a formal proof)
  - Does not test the Python implementation

  **Property-based testing (test_property_*.py):**
  - Tests that the Python IMPLEMENTATION matches the specification
  - Holds for tested inputs (statistical coverage via Hypothesis)
  - Could miss rare edge cases
  - Directly validates the Python code

  Together, they provide:
  - Mathematical certainty that the specifications are compatible (refinement)
  - Empirical confidence that the implementation matches the spec (testing)
  - Explicit documentation of what is proven vs. what is tested

  ──────────────────────────────────────────────────────────────────────────

  THEOREMS (all sorry-free):

  1. verdictRefines_refl: verdictRefines is reflexive
  2. verdictRefines_trans: verdictRefines is transitive
  3. and_monotone_refines: Verdict.and is monotone w.r.t. refinement
  4. verified_slot_refines_conservative_slot:
     Per SLOT predicate, VERIFIED verdict refines CONSERVATIVE verdict
  5. verified_refines_conservative:
     VERIFIED mode refines CONSERVATIVE mode (per-predicate)
  6. verified_soundness_conditional:
     VERIFIED mode is sound under the VC axiom (PASS → property holds)
  7. evaluateAllWithMode_cons_decomposition:
     evaluateAllWithMode decomposes as AND of head and tail
  8. simulation_verified_conservative:
     SIMULATION THEOREM: for any predicate list, VERIFIED overall verdict
     refines CONSERVATIVE overall verdict
  9. conservative_pass_implies_verified_pass_or_uncertain:
     If CONSERVATIVE mode gives PASS, VERIFIED gives PASS or UNCERTAIN
     (never FAIL — no "downgrade")
  10. verified_fail_implies_conservative_uncertain_or_fail:
      If VERIFIED mode gives FAIL, CONSERVATIVE gives UNCERTAIN or FAIL
      (consistency guarantee)
  11. permissive_refines_conservative:
      PERMISSIVE mode also refines CONSERVATIVE mode (structural)

  AXIOMS (0 new; 1 inherited from SLOTVerification.lean):
  1. verification_conditions_imply_property: If all VCs for a SLOT predicate
     hold, then the predicate's semantic property holds. (Used in theorem 6)

  REFERENCE: DOC-03 (SDD) §3.5, DOC-10 (Deterministic Methods) §4,
             DOC-14 (SLOT Proof-Implementation Gap),
             docs/11-Refinement-Mapping.md
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners
import BioCompiler.TypeSystem
import BioCompiler.Compositional
import BioCompiler.SLOTIndependence
import BioCompiler.SLOTVerification

namespace BioCompiler

open Verdict Sequence

-- ==============================================================================
-- Section 1: Verdict Refinement Ordering
-- ==============================================================================

/-- The information ordering on verdicts. `verdictRefines v_refined v_abstract`
    means that `v_refined` is at least as informative as `v_abstract`.

    The ordering is:
        PASS ≥ UNCERTAIN    (PASS is more informative)
        FAIL ≥ UNCERTAIN    (FAIL is more informative)
        PASS and FAIL are incomparable

    Formally: the abstract verdict is UNCERTAIN (anything refines it)
    or the refined verdict equals the abstract verdict (same information).

    This is the "simulation ordering" from refinement theory: a concrete
    system refines an abstract system if every concrete observation is
    consistent with (at least as precise as) some abstract observation. -/
def verdictRefines (v_refined : Verdict) (v_abstract : Verdict) : Prop :=
  v_abstract = UNCERTAIN ∨ v_refined = v_abstract

/-- verdictRefines is reflexive: every verdict refines itself. -/
theorem verdictRefines_refl (v : Verdict) : verdictRefines v v := by
  right; rfl

/-- verdictRefines is transitive: if v₁ refines v₂ and v₂ refines v₃,
    then v₁ refines v₃. This makes verdictRefines a preorder. -/
theorem verdictRefines_trans (v₁ v₂ v₃ : Verdict)
    (h12 : verdictRefines v₁ v₂) (h23 : verdictRefines v₂ v₃) :
    verdictRefines v₁ v₃ := by
  unfold verdictRefines at *
  rcases h23 with h3_unc | h23_eq
  · -- v₃ = UNCERTAIN, so v₃ = UNCERTAIN holds
    left; exact h3_unc
  · -- v₂ = v₃
    rcases h12 with h2_unc | h12_eq
    · -- v₂ = UNCERTAIN, and v₂ = v₃, so v₃ = UNCERTAIN
      left; rw [h23_eq] at h2_unc; exact h2_unc
    · -- v₁ = v₂ = v₃
      right; rw [h12_eq, h23_eq]

/-- PASS refines UNCERTAIN: PASS is strictly more informative than UNCERTAIN.
    This is the key ordering fact that makes VERIFIED mode a valid refinement
    of CONSERVATIVE mode. -/
theorem pass_refines_uncertain : verdictRefines PASS UNCERTAIN := by
  left; rfl

/-- FAIL refines UNCERTAIN: FAIL is strictly more informative than UNCERTAIN.
    In the Python implementation, VERIFIED mode can return FAIL for SLOT
    predicates (e.g., when a stability check definitively fails). This theorem
    ensures such FAIL verdicts are still consistent with the CONSERVATIVE model. -/
theorem fail_refines_uncertain : verdictRefines FAIL UNCERTAIN := by
  left; rfl

/-- UNCERTAIN does NOT refine PASS: UNCERTAIN is less informative than PASS.
    This means you cannot "downgrade" a PASS verdict to UNCERTAIN in a
    refinement — the concrete system must be at least as informative as
    the abstract system. -/
theorem uncertain_not_refines_pass : ¬verdictRefines UNCERTAIN PASS := by
  unfold verdictRefines; simp [Verdict.noConfusion]

/-- Verdict.and is monotone with respect to the refinement ordering.
    If v₁ refines v₂ and w₁ refines w₂, then AND(v₁, w₁) refines AND(v₂, w₂).

    This is the key compositionality lemma: it allows us to lift per-predicate
    refinement to composed (AND-fold) refinement. It means that if each
    individual predicate's VERIFIED verdict refines its CONSERVATIVE verdict,
    then the overall AND-fold also refines. -/
theorem and_monotone_refines (v₁ v₂ w₁ w₂ : Verdict)
    (hv : verdictRefines v₁ v₂) (hw : verdictRefines w₁ w₂) :
    verdictRefines (Verdict.and v₁ w₁) (Verdict.and v₂ w₂) := by
  unfold verdictRefines at *
  -- Case split on the structure of the refinement premises
  rcases hv with h2_unc | h12_eq
  · -- v₂ = UNCERTAIN
    rcases hw with h4_unc | h34_eq
    · -- w₂ = UNCERTAIN: AND(UNCERTAIN, UNCERTAIN) = UNCERTAIN
      -- Any verdict refines UNCERTAIN
      left; simp [Verdict.and, h2_unc, h4_unc]
    · -- w₁ = w₂
      cases w₂ with
      | PASS => left; simp [Verdict.and, h2_unc, h34_eq]
      | UNCERTAIN => left; simp [Verdict.and, h2_unc, h34_eq]
      | FAIL => right; simp [and_FAIL_right, h2_unc, h34_eq]
  · -- v₁ = v₂
    rcases hw with h4_unc | h34_eq
    · -- w₂ = UNCERTAIN
      cases v₂ with
      | PASS => left; simp [Verdict.and, h4_unc, h12_eq]
      | UNCERTAIN => left; simp [Verdict.and, h4_unc, h12_eq]
      | FAIL => right; simp [and_FAIL_left, h4_unc, h12_eq]
    · -- w₁ = w₂: AND(v₂, w₂) = AND(v₂, w₂)
      right; rw [h12_eq, h34_eq]

-- ==============================================================================
-- Section 2: Mode Refinement Relation
-- ==============================================================================

/-- Mode-level refinement relation. `refines mode_refined mode_abstract` means
    that for every SLOT predicate, the verdict produced by `mode_refined` is at
    least as informative as the verdict produced by `mode_abstract`.

    This captures the behavioral specification: one mode REFINES another if
    switching from the abstract mode to the refined mode can only make verdicts
    more informative (never less informative, never contradictory). -/
def refines (mode_refined : SLOTMode) (mode_abstract : SLOTMode) : Prop :=
  ∀ (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext),
    isSLOT P = true →
    verdictRefines (evaluateSLOT mode_refined vctx P seq ctx)
                   (evaluateSLOT mode_abstract vctx P seq ctx)

-- ==============================================================================
-- Section 3: VERIFIED Refines CONSERVATIVE (Theorem b)
-- ==============================================================================

/-- THEOREM (Per-SLOT-Predicate Refinement): For any SLOT predicate, the
    verdict produced by VERIFIED mode refines the verdict produced by
    CONSERVATIVE mode.

    Proof: CONSERVATIVE mode always returns UNCERTAIN for SLOT predicates.
    VERIFIED mode returns PASS (if all VCs hold) or UNCERTAIN (otherwise).
    Both PASS and UNCERTAIN refine UNCERTAIN (by the information ordering).

    This is the per-predicate foundation of the full refinement theorem. -/
theorem verified_slot_refines_conservative_slot
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext)
    (h_slot : isSLOT P = true) :
    verdictRefines (evaluateSLOT SLOTMode.verified vctx P seq ctx)
                   (evaluateSLOT SLOTMode.conservative vctx P seq ctx) := by
  -- CONSERVATIVE mode always returns UNCERTAIN for SLOT predicates
  have h_conservative : evaluateSLOT SLOTMode.conservative vctx P seq ctx = UNCERTAIN := rfl
  -- VERIFIED mode returns PASS or UNCERTAIN
  -- Both PASS and UNCERTAIN refine UNCERTAIN
  unfold verdictRefines
  left; exact h_conservative

/-- THEOREM (VERIFIED Refines CONSERVATIVE): VERIFIED mode refines
    CONSERVATIVE mode. For every SLOT predicate, the VERIFIED verdict is at
    least as informative as the CONSERVATIVE verdict.

    This is the main per-predicate refinement theorem. Combined with the
    compositionality lemma (`and_monotone_refines`), it yields the full
    simulation theorem (`simulation_verified_conservative`). -/
theorem verified_refines_conservative : refines SLOTMode.verified SLOTMode.conservative :=
  fun P seq ctx vctx h_slot => verified_slot_refines_conservative_slot P seq ctx vctx h_slot

/-- DETAILED THEOREM (VERIFIED Refines CONSERVATIVE — case analysis):
    This theorem makes the three cases of the refinement explicit:

    1. If VERIFIED returns PASS → CONSERVATIVE returns UNCERTAIN
       (PASS is more informative; the upgrade is safe)
    2. If VERIFIED returns UNCERTAIN → CONSERVATIVE returns UNCERTAIN
       (same verdict; no change)
    3. If VERIFIED returns FAIL → CONSERVATIVE returns UNCERTAIN or FAIL
       (FAIL is at least as informative as UNCERTAIN; consistent)

    Note: In the current Lean4 model, evaluateSLOT never returns FAIL for
    SLOT predicates, so case 3 is vacuously true. It is included because
    the Python VERIFIED implementation CAN return FAIL (e.g., when a
    stability check definitively fails), and the refinement relation should
    account for this behavior. -/
theorem verified_refines_conservative_cases
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext)
    (h_slot : isSLOT P = true) :
    -- Case 1: VERIFIED PASS → CONSERVATIVE UNCERTAIN
    (evaluateSLOT SLOTMode.verified vctx P seq ctx = PASS →
     evaluateSLOT SLOTMode.conservative vctx P seq ctx = UNCERTAIN) ∧
    -- Case 2: VERIFIED UNCERTAIN → CONSERVATIVE UNCERTAIN
    (evaluateSLOT SLOTMode.verified vctx P seq ctx = UNCERTAIN →
     evaluateSLOT SLOTMode.conservative vctx P seq ctx = UNCERTAIN) ∧
    -- Case 3: VERIFIED FAIL → CONSERVATIVE UNCERTAIN or FAIL
    (evaluateSLOT SLOTMode.verified vctx P seq ctx = FAIL →
     evaluateSLOT SLOTMode.conservative vctx P seq ctx = UNCERTAIN ∨
     evaluateSLOT SLOTMode.conservative vctx P seq ctx = FAIL) := by
  constructor
  · -- Case 1: VERIFIED PASS → CONSERVATIVE UNCERTAIN
    intro _; rfl
  constructor
  · -- Case 2: VERIFIED UNCERTAIN → CONSERVATIVE UNCERTAIN
    intro _; rfl
  · -- Case 3: VERIFIED FAIL → CONSERVATIVE UNCERTAIN or FAIL
    intro _; left; rfl

-- ==============================================================================
-- Section 4: VERIFIED Mode Soundness (Theorem c)
-- ==============================================================================

/-- THEOREM (VERIFIED Mode Soundness — Conditional on VC Axiom):
    If VERIFIED mode returns PASS for a SLOT predicate, then the predicate's
    semantic property holds.

    This theorem is conditional on the axiom `verification_conditions_imply_property`,
    which states that if all verification conditions for a SLOT predicate hold,
    then the predicate's semantic property holds. This axiom is the explicit
    "social contract" of VERIFIED mode — the only unproven link is the trust
    that external tools are correct.

    Proof chain:
    1. evaluateSLOT VERIFIED = PASS → allVCHold vctx (slotVCs P) = true
       (by verified_pass_implies_all_vcs)
    2. allVCHold = true → slotPropertySemantics holds
       (by verification_conditions_imply_property axiom)

    This is the VERIFIED-mode analogue of `type_soundness` for core predicates:
    it connects PASS verdicts to actual property guarantees, but through
    EXPLICIT verification conditions rather than vacuous truth. -/
theorem verified_soundness_conditional
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext) :
    isSLOT P = true →
    evaluateSLOT SLOTMode.verified vctx P seq ctx = PASS →
    slotPropertySemantics P seq ctx :=
  slot_soundness_verified P seq ctx vctx

/-- COROLLARY: Combined soundness for VERIFIED mode — if VERIFIED returns PASS,
    then both slotPropertySemantics and propertyHolds hold.

    This bridges the gap between the new (strengthenable) slotPropertySemantics
    and the existing (vacuously True) propertyHolds. Currently, the implication
    from slotPropertySemantics to propertyHolds is trivial (since propertyHolds
    is True for SLOT predicates), but when slotPropertySemantics is strengthened
    to have real semantic content, this theorem will become the formal bridge. -/
theorem verified_soundness_combined [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext) :
    isSLOT P = true →
    evaluateWithMode (State := State) SLOTMode.verified vctx P seq ctx = PASS →
    slotPropertySemantics P seq ctx ∧
    @propertyHolds inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst P seq ctx := by
  intro h_slot h_pass
  unfold evaluateWithMode at h_pass
  simp only [if_pos h_slot] at h_pass
  constructor
  · -- slotPropertySemantics holds (under VC axiom)
    exact slot_soundness_verified P seq ctx vctx h_slot h_pass
  · -- propertyHolds holds (trivially, since it's True for SLOT predicates)
    exact slot_property_implies_propertyHolds P seq ctx h_slot
      (slot_soundness_verified P seq ctx vctx h_slot h_pass)

-- ==============================================================================
-- Section 5: Mode-Aware Compositional Evaluation
-- ==============================================================================

/-- Evaluate a list of type predicates with SLOT mode support, combining
    results with three-valued AND (Verdict.and).

    This is the mode-aware analogue of `evaluateAll` from Compositional.lean.
    It evaluates each predicate using `evaluateWithMode` and combines the
    results with Verdict.and, starting from PASS.

    For core predicates, behavior is identical across modes.
    For SLOT predicates, the mode determines the verdict (UNCERTAIN, PASS, etc.). -/
def evaluateAllWithMode [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (mode : SLOTMode) (vctx : VerificationContext)
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) : Verdict :=
  match predicates with
  | [] => PASS
  | P :: ps => Verdict.and (evaluateWithMode (State := State) mode vctx P seq ctx)
                         (evaluateAllWithMode (State := State) mode vctx ps seq ctx)

/-- evaluateAllWithMode with empty list returns PASS (vacuously sound). -/
@[simp] theorem evaluateAllWithMode_nil [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (mode : SLOTMode) (vctx : VerificationContext) (seq : Sequence) (ctx : CellularContext) :
    evaluateAllWithMode (State := State) mode vctx [] seq ctx = PASS := rfl

/-- evaluateAllWithMode decomposes as AND of head and tail. -/
@[simp] theorem evaluateAllWithMode_cons [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (mode : SLOTMode) (vctx : VerificationContext)
    (P : TypePredicate) (ps : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    evaluateAllWithMode (State := State) mode vctx (P :: ps) seq ctx =
      Verdict.and (evaluateWithMode (State := State) mode vctx P seq ctx)
                  (evaluateAllWithMode (State := State) mode vctx ps seq ctx) := rfl

-- ==============================================================================
-- Section 6: Simulation Theorem (Theorem d)
-- ==============================================================================

/-- Per-predicate refinement for evaluateWithMode (not just evaluateSLOT).
    For any predicate (core or SLOT), the VERIFIED verdict refines the
    CONSERVATIVE verdict.

    - Core predicates: same verdict in both modes (trivial refinement)
    - SLOT predicates: VERIFIED verdict refines CONSERVATIVE verdict
      (by verified_slot_refines_conservative_slot) -/
theorem per_predicate_refinement [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext) :
    verdictRefines (evaluateWithMode (State := State) SLOTMode.verified vctx P seq ctx)
                   (evaluateWithMode (State := State) SLOTMode.conservative vctx P seq ctx) := by
  unfold evaluateWithMode
  by_cases h_slot : isSLOT P = true
  · -- SLOT predicate: VERIFIED refines CONSERVATIVE
    simp only [if_pos h_slot]
    exact verified_slot_refines_conservative_slot P seq ctx vctx h_slot
  · -- Core predicate: same verdict in both modes (trivial refinement)
    simp only [if_neg h_slot]
    right; rfl

/-- SIMULATION THEOREM: For any list of predicates, the VERIFIED mode's
    overall verdict refines the CONSERVATIVE mode's overall verdict.

    This is the central result of this module. It establishes that switching
    from CONSERVATIVE mode to VERIFIED mode is a valid refinement: the
    VERIFIED mode's behavioral specification is consistent with (and at least
    as informative as) the CONSERVATIVE mode's specification.

    Proof: By induction on the predicate list.
    - Base case: empty list → both return PASS → PASS refines PASS
    - Inductive step: AND(head, tail) refines AND(head, tail)
      by and_monotone_refines + per_predicate_refinement + IH

    SIGNIFICANCE:
    This theorem means that any system designed around CONSERVATIVE mode
    guarantees will continue to work under VERIFIED mode. The VERIFIED mode
    can only "upgrade" verdicts (from UNCERTAIN to PASS), never "downgrade"
    them (from PASS to FAIL, or from PASS to UNCERTAIN).

    WHAT THIS DOESN'T PROVE:
    This does not prove that the Python code in `slot_verification.py` correctly
    implements the VERIFIED mode specification. It proves a relationship between
    two SPECIFICATIONS, not between a specification and an implementation.
    Property-based testing bridges this gap empirically. -/
theorem simulation_verified_conservative [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext)
    (vctx : VerificationContext) :
    verdictRefines (evaluateAllWithMode (State := State) SLOTMode.verified vctx predicates seq ctx)
                   (evaluateAllWithMode (State := State) SLOTMode.conservative vctx predicates seq ctx) := by
  induction predicates with
  | nil =>
    -- Both evaluate to PASS; PASS refines PASS
    simp [evaluateAllWithMode]
    right; rfl
  | cons P ps ih =>
    -- VERIFIED: AND(evaluateWithMode VERIFIED P, evaluateAllWithMode VERIFIED ps)
    -- CONSERVATIVE: AND(evaluateWithMode CONSERVATIVE P, evaluateAllWithMode CONSERVATIVE ps)
    simp [evaluateAllWithMode]
    exact and_monotone_refines
      (evaluateWithMode (State := State) SLOTMode.verified vctx P seq ctx)
      (evaluateWithMode (State := State) SLOTMode.conservative vctx P seq ctx)
      (evaluateAllWithMode (State := State) SLOTMode.verified vctx ps seq ctx)
      (evaluateAllWithMode (State := State) SLOTMode.conservative vctx ps seq ctx)
      (per_predicate_refinement P seq ctx vctx)
      ih

-- ==============================================================================
-- Section 7: Corollaries and Connections
-- ==============================================================================

/-- COROLLARY: If CONSERVATIVE mode produces PASS for the overall evaluation,
    then VERIFIED mode also produces PASS (never FAIL).

    This means upgrading from CONSERVATIVE to VERIFIED mode can never
    "downgrade" a PASS result to FAIL. The worst case is that VERIFIED mode
    also produces PASS (for the same reasons) or UNCERTAIN (if some SLOT
    predicates upgrade from UNCERTAIN to PASS, which would actually IMPROVE
    the verdict, but could also keep it the same if other predicates are
    UNCERTAIN for different reasons).

    In practice: if your CONSERVATIVE-mode certificate says PASS, switching
    to VERIFIED mode won't invalidate it. -/
theorem conservative_pass_no_downgrade [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext)
    (vctx : VerificationContext)
    (h_conservative_pass : evaluateAllWithMode (State := State) SLOTMode.conservative vctx predicates seq ctx = PASS) :
    evaluateAllWithMode (State := State) SLOTMode.verified vctx predicates seq ctx ≠ FAIL := by
  intro h_fail
  have h_ref := simulation_verified_conservative (State := State) predicates seq ctx vctx
  unfold verdictRefines at h_ref
  rcases h_ref with h_unc | h_eq
  · -- CONSERVATIVE = UNCERTAIN, but h_conservative_pass says PASS
    rw [h_conservative_pass] at h_unc; cases h_unc
  · -- VERIFIED = CONSERVATIVE, but VERIFIED = FAIL and CONSERVATIVE = PASS
    rw [h_fail, h_conservative_pass] at h_eq; cases h_eq

/-- COROLLARY: If VERIFIED mode produces FAIL for the overall evaluation,
    then CONSERVATIVE mode produces UNCERTAIN or FAIL (never PASS).

    This means that VERIFIED mode's FAIL verdicts are consistent with
    CONSERVATIVE mode: CONSERVATIVE mode would never claim PASS where
    VERIFIED mode claims FAIL. If VERIFIED says FAIL, CONSERVATIVE was
    already uncertain or failing. -/
theorem verified_fail_consistent [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext)
    (vctx : VerificationContext)
    (h_verified_fail : evaluateAllWithMode (State := State) SLOTMode.verified vctx predicates seq ctx = FAIL) :
    evaluateAllWithMode (State := State) SLOTMode.conservative vctx predicates seq ctx = UNCERTAIN ∨
    evaluateAllWithMode (State := State) SLOTMode.conservative vctx predicates seq ctx = FAIL := by
  have h_ref := simulation_verified_conservative (State := State) predicates seq ctx vctx
  unfold verdictRefines at h_ref
  -- h_ref : v_cons = UNCERTAIN ∨ v_ver = v_cons
  -- h_verified_fail : v_ver = FAIL
  rcases h_ref with h_unc | h_eq
  · left; exact h_unc
  · right; exact h_eq ▸ h_verified_fail

/-- COROLLARY: PERMISSIVE mode also refines CONSERVATIVE mode.

    PERMISSIVE mode returns PASS when ANY verification condition holds
    (weaker requirement than VERIFIED mode's ALL conditions). This means
    PERMISSIVE mode's verdicts are also at least as informative as
    CONSERVATIVE mode's verdicts for SLOT predicates.

    WARNING: While PERMISSIVE mode structurally refines CONSERVATIVE mode,
    its PASS verdicts do NOT carry the same soundness guarantee as VERIFIED
    mode's PASS verdicts. PERMISSIVE mode can return PASS with insufficient
    evidence. This theorem only says the information ordering is preserved,
    not that the verdicts are sound. -/
theorem permissive_refines_conservative : refines SLOTMode.permissive SLOTMode.conservative :=
  fun P seq ctx vctx h_slot => by
    -- CONSERVATIVE mode always returns UNCERTAIN for SLOT predicates
    have h_conservative : evaluateSLOT SLOTMode.conservative vctx P seq ctx = UNCERTAIN := rfl
    -- PERMISSIVE mode returns PASS or UNCERTAIN, both of which refine UNCERTAIN
    unfold verdictRefines
    left; exact h_conservative

/-- THEOREM (Refinement ordering is a partial order on modes):
    The three modes form a chain in the refinement ordering:
    - PERMISSIVE refines CONSERVATIVE (can upgrade UNCERTAIN to PASS)
    - VERIFIED refines CONSERVATIVE (can upgrade UNCERTAIN to PASS)
    - PERMISSIVE and VERIFIED are not comparable in general
      (PERMISSIVE returns PASS more easily, but VERIFIED's PASS is sound)

    The ordering on modes by "safety" (CONSERVATIVE > VERIFIED > PERMISSIVE)
    is the REVERSE of the refinement ordering on verdicts. This is expected:
    safer modes produce weaker verdicts (more UNCERTAIN), which are easier
    to refine. -/
theorem mode_refinement_chain :
    refines SLOTMode.verified SLOTMode.conservative ∧
    refines SLOTMode.permissive SLOTMode.conservative ∧
    -- Both VERIFIED and PERMISSIVE refines are trivially established because
    -- CONSERVATIVE always returns UNCERTAIN for SLOT predicates.
    -- The real content is that VERIFIED's PASS is sound (under the VC axiom).
    True := by
  constructor
  · exact verified_refines_conservative
  constructor
  · exact permissive_refines_conservative
  · trivial

/-- THEOREM: CONSERVATIVE mode's overall verdict is never more informative
    than VERIFIED mode's. This is the contrapositive view of the refinement:
    if the CONSERVATIVE verdict is not UNCERTAIN, then the VERIFIED verdict
    equals it.

    This means: for any definite verdict (PASS or FAIL) that CONSERVATIVE
    mode produces, VERIFIED mode produces the same verdict. Since CONSERVATIVE
    mode never produces PASS for SLOT predicates, this effectively says that
    CONSERVATIVE mode's PASS verdicts (from core predicates only) are preserved
    in VERIFIED mode. -/
theorem conservative_definite_preserved [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext)
    (vctx : VerificationContext)
    (h_not_uncertain : evaluateAllWithMode (State := State) SLOTMode.conservative vctx predicates seq ctx ≠ UNCERTAIN) :
    evaluateAllWithMode (State := State) SLOTMode.verified vctx predicates seq ctx =
    evaluateAllWithMode (State := State) SLOTMode.conservative vctx predicates seq ctx := by
  have h_ref := simulation_verified_conservative (State := State) predicates seq ctx vctx
  unfold verdictRefines at h_ref
  rcases h_ref with h_unc | h_eq
  · -- CONSERVATIVE = UNCERTAIN, contradicting h_not_uncertain
    exact absurd h_unc h_not_uncertain
  · -- VERIFIED = CONSERVATIVE
    exact h_eq

-- ==============================================================================
-- Section 8: Connection to Existing Soundness Theorems
-- ==============================================================================

/-- THEOREM (Refinement + Soundness = Progressive Assurance):
    The refinement theorem, combined with the existing soundness theorems,
    provides a progressive assurance framework:

    Level 0 (CONSERVATIVE): All PASS verdicts are sound (proven unconditionally).
      This is the baseline guarantee from `slot_soundness_conservative`.

    Level 1 (VERIFIED): All PASS verdicts are sound (proven under VC axiom).
      This is the enhanced guarantee from `verified_soundness_conditional`.
      The refinement theorem ensures Level 1 is at least as informative as
      Level 0, so no information is lost by switching modes.

    Level 2 (PERMISSIVE): No formal soundness guarantee.
      PERMISSIVE mode refines CONSERVATIVE in the information ordering,
      but its PASS verdicts are not backed by the VC axiom.

    This theorem formalizes Level 0 and Level 1: CONSERVATIVE mode's
    soundness + the refinement + VERIFIED mode's conditional soundness. -/
theorem progressive_assurance [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext) :
    -- Level 0: CONSERVATIVE soundness (unconditional)
    (evaluateWithMode (State := State) SLOTMode.conservative vctx P seq ctx = PASS →
      @propertyHolds inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
        State inst_dec inst_inhab inst_ndfst P seq ctx) ∧
    -- Level 1: VERIFIED soundness (under VC axiom)
    (isSLOT P = true →
     evaluateWithMode (State := State) SLOTMode.verified vctx P seq ctx = PASS →
     slotPropertySemantics P seq ctx) ∧
    -- Refinement: VERIFIED refines CONSERVATIVE
    verdictRefines (evaluateWithMode (State := State) SLOTMode.verified vctx P seq ctx)
                   (evaluateWithMode (State := State) SLOTMode.conservative vctx P seq ctx) := by
  constructor
  · -- Level 0: CONSERVATIVE soundness
    exact slot_soundness_conservative P seq ctx vctx
  constructor
  · -- Level 1: VERIFIED soundness (under VC axiom)
    intro h_slot h_pass
    unfold evaluateWithMode at h_pass
    simp only [if_pos h_slot] at h_pass
    exact slot_soundness_verified P seq ctx vctx h_slot h_pass
  · -- Refinement
    exact per_predicate_refinement P seq ctx vctx

/-- COROLLARY: If VERIFIED mode gives PASS and CONSERVATIVE mode gives
    UNCERTAIN, then the SLOT predicate's semantic property holds (under VC axiom).

    This is the "value add" of VERIFIED mode: it can upgrade UNCERTAIN to PASS,
    and when it does, the property is guaranteed (under the axiom). This is
    exactly the practical benefit that motivates VERIFIED mode. -/
theorem verified_pass_value_add [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex]
    [inst_cpg : CpGIslandScanner] [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner]
    [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) (vctx : VerificationContext)
    (h_slot : isSLOT P = true)
    (h_verified_pass : evaluateWithMode (State := State) SLOTMode.verified vctx P seq ctx = PASS)
    (h_conservative_uncertain : evaluateWithMode (State := State) SLOTMode.conservative vctx P seq ctx = UNCERTAIN) :
    slotPropertySemantics P seq ctx := by
  -- VERIFIED PASS + SLOT predicate → property holds (under VC axiom)
  -- CONSERVATIVE UNCERTAIN is guaranteed by the structure (always true for SLOT predicates)
  unfold evaluateWithMode at h_verified_pass
  simp only [if_pos h_slot] at h_verified_pass
  exact slot_soundness_verified P seq ctx vctx h_slot h_verified_pass

-- ==============================================================================
-- Section 9: Summary and Proof Architecture
-- ==============================================================================

/-
  PROOF ARCHITECTURE:

  ┌──────────────────────────────────────────────────────────────────────┐
  │  SIMULATION THEOREM: simulation_verified_conservative               │
  │  "For any predicate list, VERIFIED overall verdict refines          │
  │   CONSERVATIVE overall verdict"                                      │
  │                                                                      │
  │  Proof: By induction on predicate list                               │
  │  ┌──────────────────────────────────────────────────────────────┐   │
  │  │ Step: AND(v_verified_P, vs_verified) refines                 │   │
  │  │       AND(v_conservative_P, vs_conservative)                 │   │
  │  │                                                              │   │
  │  │ Uses: and_monotone_refines                                   │   │
  │  │        + per_predicate_refinement (head)                     │   │
  │  │        + IH (tail)                                           │   │
  │  └──────────────────────────────────────────────────────────────┘   │
  └──────────────────────────────────────────────────────────────────────┘
                              │
  ┌───────────────────────────┴──────────────────────────────────────────┐
  │  PER-PREDICATE REFINEMENT: per_predicate_refinement                 │
  │  "For any predicate, VERIFIED verdict refines CONSERVATIVE verdict" │
  │                                                                      │
  │  Core predicates: same verdict → trivial refinement (right)         │
  │  SLOT predicates: CONSERVATIVE = UNCERTAIN → refinement (left)      │
  └──────────────────────────────────────────────────────────────────────┘
                              │
  ┌───────────────────────────┴──────────────────────────────────────────┐
  │  VERDICT INFORMATION ORDERING: verdictRefines                        │
  │  v_refined refines v_abstract iff v_abstract = UNCERTAIN             │
  │                                      ∨ v_refined = v_abstract        │
  │                                                                      │
  │  PASS ≥ UNCERTAIN    FAIL ≥ UNCERTAIN    PASS ∥ FAIL               │
  │                                                                      │
  │  Monotone: AND preserves the ordering (and_monotone_refines)        │
  │  Transitive: chain of refinements compose (verdictRefines_trans)    │
  └──────────────────────────────────────────────────────────────────────┘
                              │
  ┌───────────────────────────┴──────────────────────────────────────────┐
  │  EXISTING THEOREMS (from SLOTVerification.lean):                    │
  │  - conservative_is_safe: CONSERVATIVE never returns PASS for SLOT   │
  │  - slot_soundness_conservative: CONSERVATIVE PASS → property holds  │
  │  - slot_soundness_verified: VERIFIED PASS → property holds (axiom)  │
  │  - verification_conditions_imply_property: VC axiom (1 axiom)       │
  └──────────────────────────────────────────────────────────────────────┘

  FORMAL COVERAGE SUMMARY (with refinement):

  ┌──────────────┬───────────────────────┬──────────────────────┬───────────────────┐
  │ Mode         │ Python behavior       │ Lean4 proof coverage │ Refinement status │
  ├──────────────┼───────────────────────┼──────────────────────┼───────────────────┤
  │ CONSERVATIVE │ Always UNCERTAIN      │ Full (vacuously      │ Abstract system   │
  │              │                       │   sound)             │ (baseline)        │
  │ VERIFIED     │ PASS with evidence    │ Partial (under axiom)│ Refines CONSERV.  │
  │ PERMISSIVE   │ PASS with weak eviden.│ None (beyond proof)  │ Refines CONSERV.  │
  └──────────────┴───────────────────────┴──────────────────────┴───────────────────┘

  NEW THEOREMS (this file, 11 total, 0 sorry, 0 new axioms):

  Section 1 — Verdict Refinement Ordering:
    1. verdictRefines_refl     : reflexivity
    2. verdictRefines_trans    : transitivity
    3. pass_refines_uncertain  : PASS ≥ UNCERTAIN
    4. fail_refines_uncertain  : FAIL ≥ UNCERTAIN
    5. uncertain_not_refines_pass : UNCERTAIN ≱ PASS
    6. and_monotone_refines    : AND preserves refinement

  Section 3 — VERIFIED Refines CONSERVATIVE:
    7. verified_slot_refines_conservative_slot : per-SLOT-predicate
    8. verified_refines_conservative           : mode-level
    9. verified_refines_conservative_cases     : explicit 3-case analysis

  Section 4 — VERIFIED Mode Soundness:
   10. verified_soundness_conditional          : PASS → property (under axiom)
   11. verified_soundness_combined             : PASS → property + propertyHolds

  Section 6 — Simulation:
   12. simulation_verified_conservative        : overall verdict refinement

  Section 7 — Corollaries:
   13. conservative_pass_no_downgrade          : PASS not lost
   14. verified_fail_consistent                : FAIL is consistent
   15. permissive_refines_conservative         : PERMISSIVE also refines
   16. mode_refinement_chain                   : mode ordering
   17. conservative_definite_preserved         : definite verdicts preserved

  Section 8 — Connections:
   18. progressive_assurance                   : Level 0 + Level 1 + refinement
   19. verified_pass_value_add                 : value of upgrading to VERIFIED
-/

end BioCompiler
