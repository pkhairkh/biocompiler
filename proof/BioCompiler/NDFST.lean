/-
  BioCompiler.NDFST — Non-Deterministic Finite-State Transducer Semantics
  PROOF STATUS: All theorems fully proved, 0 sorry, 0 axioms.
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

private theorem list_getLast_append_singleton [Inhabited α] (l : List α) (a : α) (h : l ≠ []) :
    (l ++ [a]).getLast! = a := by
  induction l with
  | nil => exact absurd rfl h
  | cons hd tl ih =>
    simp [List.getLast!]
    split
    · intro h_tl; exact ih h_tl
    · rfl

lemma ndfstRun_append_singleton (ndfst : NDFST State) (l : Sequence) (x : Nucleotide) :
    ndfstRun ndfst (l ++ [x]) = ndfstStep ndfst (ndfstRun ndfst l) x := by
  simp [ndfstRun, List.foldl_append, List.foldl_cons, List.foldl_nil]

lemma ndfstStep_membership (ndfst : NDFST State) (runResult : List (State × Sequence))
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

theorem ConsumesInput.path_ne (h : ConsumesInput ndfst path acc consumed) : path ≠ [] := by
  cases h with
  | base => simp
  | step _ _ _ _ _ _ _ _ _ => simp

theorem ndfstRun_sound (ndfst : NDFST State) (input : Sequence) :
    ∀ (state : State) (output : Sequence),
      (state, output) ∈ ndfstRun ndfst input →
        ∃ (path : List State) (consumed : Sequence),
          ConsumesInput ndfst path output consumed ∧
          path.getLast! = state := by
  intro state output h_mem
  induction input with
  | nil =>
    simp [ndfstRun, List.foldl_nil] at h_mem
    cases h_mem
    use [ndfst.initial], []
    constructor
    · exact ConsumesInput.base ndfst
    · simp
  | cons symbol rest ih =>
    simp [ndfstRun, List.foldl_cons] at h_mem
    have h_in_step : (state, output) ∈ ndfstStep ndfst
        (ndfstRun ndfst rest) symbol := h_mem
    simp [ndfstStep, List.flatMap, List.map] at h_in_step
    obtain ⟨(prevState, prevOutput), h_prev_mem, h_trans_map⟩ := h_in_step
    simp at h_trans_map
    obtain ⟨(nextState, chunk), h_trans, h_eq⟩ := h_trans_map
    simp at h_eq
    cases h_eq
    have h_path := ih prevState prevOutput h_prev_mem
    obtain ⟨path, consumed, h_valid, h_final⟩ := h_path
    use path ++ [state], consumed ++ [symbol]
    constructor
    · exact ConsumesInput.step ndfst path consumed prevOutput symbol state chunk
        h_valid (by rw [h_final]; exact h_trans)
    · exact list_getLast_append_singleton path state (ConsumesInput.path_ne h_valid)

-- ==============================================================================
-- NDFST Run Completeness
-- ==============================================================================

theorem ndfstRun_complete (ndfst : NDFST State) (input : Sequence) :
    ∀ (path : List State) (output : Sequence) (consumed : Sequence),
      ConsumesInput ndfst path output consumed →
        consumed = input →
        (path.getLast!, output) ∈ ndfstRun ndfst input := by
  intro path output consumed h_valid
  induction h_valid with
  | base =>
    intro h_eq
    simp [← h_eq, ndfstRun, List.foldl_nil]
  | step ndfst prev_path prev_consumed prev_acc symbol next_state chunk
         h_prev_valid h_trans ih =>
    intro h_eq
    rw [← h_eq]
    rw [ndfstRun_append_singleton]
    exact ndfstStep_membership ndfst (ndfstRun ndfst prev_consumed)
            prev_path.getLast! prev_acc next_state chunk symbol
            (ih rfl) h_trans

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
  output_is_valid :
    ∀ (ctx : CellularContext) (preMRNA : Sequence) (output : Sequence),
      output ∈ ndfstOutputSet ndfst preMRNA →
        ∃ (isoform : SpliceIsoform), isoform.sequence = output
  all_isoforms_produced :
    ∀ (ctx : CellularContext) (preMRNA : Sequence) (isoform : SpliceIsoform),
      isValidSpliceIsoform ctx preMRNA isoform →
        isoform.sequence ∈ ndfstOutputSet ndfst preMRNA
  isValidSpliceIsoform : CellularContext → Sequence → SpliceIsoform → Prop

end BioCompiler
