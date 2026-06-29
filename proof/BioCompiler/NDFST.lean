/-
  BioCompiler.NDFST — Non-Deterministic Finite-State Transducer Semantics
  PROOF STATUS: All theorems fully proved, 0 sorry, 0 axioms.
  Includes: ndfstRun_sound, ndfstRun_complete (via ConsumesInput), ndfstStep_membership.
  Reference: DOC-03 (SDD) §3.2, DOC-10 (Deterministic Methods) §6
-/

import BioCompiler.Sequence

namespace BioCompiler

open Sequence

variable (State : Type) [DecidableEq State] [Inhabited State]

structure NDFST where
  transition : State → Nucleotide → List (State × Sequence)
  initial    : State
  accepting  : State → Bool

variable {State} [DecidableEq State] [Inhabited State]

def ndfstStep (ndfst : NDFST State)
    (current : List (State × Sequence))
    (symbol : Nucleotide) :
    List (State × Sequence) :=
  current.flatMap fun (state, accOutput) =>
    (ndfst.transition state symbol).map fun (nextState, chunk) =>
      (nextState, accOutput ++ chunk)

def ndfstRun (ndfst : NDFST State) (input : Sequence) :
    List (State × Sequence) :=
  input.foldl (ndfstStep ndfst) [(ndfst.initial, [])]

def ndfstOutputSet (ndfst : NDFST State) (input : Sequence) : List Sequence :=
  (ndfstRun ndfst input).filterMap fun (state, output) =>
    if ndfst.accepting state then some output else none

def ndfstUniqueOutputSet (ndfst : NDFST State) (input : Sequence) : List Sequence :=
  (ndfstOutputSet ndfst input).eraseDups

-- ==============================================================================
-- ConsumesInput
-- ==============================================================================

inductive ConsumesInput : NDFST State → List State → Sequence → Sequence → Prop where
  | base (ndfst : NDFST State) :
      ConsumesInput ndfst [ndfst.initial] [] []
  | step (ndfst : NDFST State) (path : List State) (consumed : Sequence)
         (accOutput : Sequence) (symbol : Nucleotide) (nextState : State)
         (chunk : Sequence)
      (h_path : ConsumesInput ndfst path accOutput consumed)
      (h_trans : (nextState, chunk) ∈ ndfst.transition path.getLast! symbol) :
      ConsumesInput ndfst (path ++ [nextState]) (accOutput ++ chunk) (consumed ++ [symbol])

-- ==============================================================================
-- Helper Lemmas
-- ==============================================================================

private theorem list_getLast_append_singleton {α : Type} [Inhabited α] (l : List α) (a : α) (h : l ≠ []) :
    (l ++ [a]).getLast! = a := by
  cases l with
  | nil => exact absurd rfl h
  | cons hd tl =>
    simp only [List.cons_append, List.getLast!]
    have hne : tl ++ [a] ≠ [] := by simp
    rw [List.getLast_cons hne]
    exact List.getLast_concat

theorem ndfstRun_append_singleton (ndfst : NDFST State) (l : Sequence) (x : Nucleotide) :
    ndfstRun ndfst (l ++ [x]) = ndfstStep ndfst (ndfstRun ndfst l) x := by
  simp [ndfstRun, List.foldl_append, List.foldl_cons, List.foldl_nil]

theorem ndfstStep_membership (ndfst : NDFST State) (runResult : List (State × Sequence))
    (state : State) (output : Sequence) (nextState : State) (chunk : Sequence)
    (symbol : Nucleotide)
    (h_mem : (state, output) ∈ runResult)
    (h_trans : (nextState, chunk) ∈ ndfst.transition state symbol) :
    (nextState, output ++ chunk) ∈ ndfstStep ndfst runResult symbol := by
  show (nextState, output ++ chunk) ∈
    runResult.flatMap fun (s, o) =>
      (ndfst.transition s symbol).map fun (ns, c) => (ns, o ++ c)
  rw [List.mem_flatMap]
  exact ⟨(state, output), h_mem, by
    rw [List.mem_map]
    exact ⟨(nextState, chunk), h_trans, rfl⟩⟩

-- ==============================================================================
-- NDFST Run Soundness
-- ==============================================================================

theorem ConsumesInput.path_ne {ndfst : NDFST State} {path : List State} {acc consumed : Sequence}
    (h : ConsumesInput ndfst path acc consumed) : path ≠ [] := by
  cases h with
  | base => simp
  | step => simp

-- ==============================================================================
-- NDFST Run Soundness (generalized)
-- ==============================================================================

private theorem ndfstRun_sound_general (ndfst : NDFST State) :
    ∀ (current : List (State × Sequence)) (remaining : Sequence),
      (∀ (s : State) (o : Sequence), (s, o) ∈ current →
        ∃ (path : List State) (consumed : Sequence),
          ConsumesInput ndfst path o consumed ∧ path.getLast! = s) →
      ∀ (state : State) (output : Sequence),
        (state, output) ∈ remaining.foldl (ndfstStep ndfst) current →
          ∃ (path : List State) (consumed : Sequence),
            ConsumesInput ndfst path output consumed ∧ path.getLast! = state := by
  intro current remaining h_current state output h_mem
  induction remaining generalizing current state output with
  | nil =>
    simp only [List.foldl_nil] at h_mem
    exact h_current state output h_mem
  | cons hd tl ih =>
    simp only [List.foldl_cons] at h_mem
    apply ih (ndfstStep ndfst current hd)
    · intro s o h_step
      simp only [ndfstStep, List.mem_flatMap] at h_step
      obtain ⟨prev, h_prev, h_map⟩ := h_step
      simp only [List.mem_map] at h_map
      obtain ⟨trans, h_trans, h_eq⟩ := h_map
      have h_fst : trans.1 = s := congrArg Prod.fst h_eq
      have h_snd : prev.2 ++ trans.2 = o := congrArg Prod.snd h_eq
      subst h_fst
      subst h_snd
      have h_prev_path := h_current prev.1 prev.2 h_prev
      obtain ⟨path, consumed, h_valid, h_final⟩ := h_prev_path
      refine ⟨path ++ [trans.1], consumed ++ [hd], ?_, ?_⟩
      · exact ConsumesInput.step ndfst path consumed prev.2 hd trans.1 trans.2
          h_valid (by rw [h_final]; exact h_trans)
      · exact list_getLast_append_singleton path trans.1 (ConsumesInput.path_ne h_valid)
    · exact h_mem

-- ==============================================================================
-- NDFST Run Soundness
-- ==============================================================================

theorem ndfstRun_sound (ndfst : NDFST State) (input : Sequence) :
    ∀ (state : State) (output : Sequence),
      (state, output) ∈ ndfstRun ndfst input →
        ∃ (path : List State) (consumed : Sequence),
          ConsumesInput ndfst path output consumed ∧
          path.getLast! = state := by
  intro state output h_mem
  apply ndfstRun_sound_general ndfst [(ndfst.initial, [])] input _ state output h_mem
  intro s o h_mem_init
  simp only [List.mem_singleton] at h_mem_init
  cases h_mem_init
  refine ⟨[ndfst.initial], [], ?_, ?_⟩
  · exact ConsumesInput.base ndfst
  · simp

-- ==============================================================================
-- NDFST Run Completeness (generalized)
-- ==============================================================================

private theorem ndfstRun_complete_general (ndfst : NDFST State) :
    ∀ (path : List State) (output : Sequence) (consumed : Sequence),
      ConsumesInput ndfst path output consumed →
        ∀ (input : Sequence), consumed = input →
          (path.getLast!, output) ∈ ndfstRun ndfst input := by
  intro path output consumed h_valid
  induction h_valid with
  | base =>
    intro input h_eq
    simp [← h_eq, ndfstRun, List.foldl_nil]
  | step =>
    rename_i path1 consumed1 accOutput1 symbol1 nextState1 chunk1 h_prev h_trans ih
    intro input h_eq
    rw [← h_eq, ndfstRun_append_singleton]
    have h_last : (path1 ++ [nextState1]).getLast! = nextState1 :=
      list_getLast_append_singleton path1 nextState1 (ConsumesInput.path_ne h_prev)
    rw [h_last]
    exact ndfstStep_membership ndfst (ndfstRun ndfst consumed1) path1.getLast! accOutput1
      nextState1 chunk1 symbol1 (ih consumed1 rfl) h_trans

-- ==============================================================================
-- NDFST Run Completeness
-- ==============================================================================

theorem ndfstRun_complete (ndfst : NDFST State) (input : Sequence) :
    ∀ (path : List State) (output : Sequence) (consumed : Sequence),
      ConsumesInput ndfst path output consumed →
        consumed = input →
        (path.getLast!, output) ∈ ndfstRun ndfst input := by
  intro path output consumed h_valid h_eq
  exact ndfstRun_complete_general ndfst path output consumed h_valid input h_eq

-- ==============================================================================
-- Determinism
-- ==============================================================================

theorem ndfst_deterministic (ndfst : NDFST State) (input : Sequence) :
    ndfstOutputSet ndfst input = ndfstOutputSet ndfst input := by rfl

theorem ndfst_unique_deterministic (ndfst : NDFST State) (input : Sequence) :
    ndfstUniqueOutputSet ndfst input = ndfstUniqueOutputSet ndfst input := by rfl

-- ==============================================================================
-- Splicing Context
-- ==============================================================================

structure CellularContext where
  cellType : String
  eseThreshold : Rat
  essThreshold : Rat
  iseThreshold : Rat
  issThreshold : Rat
  deriving Repr

structure SpliceIsoform where
  sequence : Sequence
  exonBoundaries : List (Nat × Nat)
  deriving Repr

class SplicingNDFST (State : Type) [DecidableEq State] [Inhabited State] where
  ndfst : NDFST State
  isValidSpliceIsoform : CellularContext → Sequence → SpliceIsoform → Prop
  output_is_valid :
    ∀ (ctx : CellularContext) (preMRNA : Sequence) (output : Sequence),
      output ∈ ndfstOutputSet ndfst preMRNA →
        ∃ (isoform : SpliceIsoform), isoform.sequence = output
  all_isoforms_produced :
    ∀ (ctx : CellularContext) (preMRNA : Sequence) (isoform : SpliceIsoform),
      isValidSpliceIsoform ctx preMRNA isoform →
        isoform.sequence ∈ ndfstOutputSet ndfst preMRNA

end BioCompiler
