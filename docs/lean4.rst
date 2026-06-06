Formal Proofs (Lean4)
=====================

BioCompiler includes a machine-checked formal proof of its type system
soundness, written in Lean4. The proof establishes that **"well-typed genes
don't go wrong"** — if the type system returns PASS for all predicates, then
all semantic properties hold.


Proof Architecture
------------------

The proof is structured as a hierarchy of theorems:

::

  Theorem 4: Compositional Soundness (Compositional.lean)
    evaluateAll = PASS → all hold
    + slot_predicates_uncertain: SLOT predicates never PASS
    + slot_predicates_dont_affect_pass: SLOT in list → no PASS
                      │
  Theorem 3: Per-Predicate Soundness (TypeSystem.lean)
    evaluate P = PASS → holds(P)  (32 predicates: 15 core + 17 SLOT)
                      │
        ┌─────────────┼─────────────┐
        │             │             │
  Theorem 1:     Theorem 2:     Theorem 5:
  3-Valued       NDFST          SLOT-Indep.
  PASS⊓PASS     Deterministic  Certs don't
  FAIL sticky   computation    depend on FFI


Module Overview
---------------

Core Logic
^^^^^^^^^^

=========================  =====================================================
Module                     Description
=========================  =====================================================
ThreeValued.lean           Three-valued logic (PASS/UNCERTAIN/FAIL) with algebraic properties
FiveValued.lean            Five-valued extension with refinement proof
Sequence.lean              Nucleotide sequences and pattern matching
NDFST.lean                 Non-deterministic finite-state transducers for splicing
=========================  =====================================================

Scanner Layer
^^^^^^^^^^^^^

=========================  =====================================================
Module                     Description
=========================  =====================================================
Scanners.lean              Abstract scanner interfaces + concrete implementations
ScannerProofs.lean         Completeness/soundness proofs for CpG, Promoter, TM scanners
OracleProofs.lean          Proofs for mRNA, co-translational folding, splicing, CAI oracles
=========================  =====================================================

Type System
^^^^^^^^^^^

=========================  =====================================================
Module                     Description
=========================  =====================================================
TypeSystem.lean            32 type predicates and per-predicate soundness
Compositional.lean         Compositional soundness: overall PASS implies all properties hold
Certificates.lean          Guarantee certificate structure and certificate soundness
SLOTIndependence.lean      SLOT independence: certificates don't depend on FFI output
=========================  =====================================================

SLOT Verification
^^^^^^^^^^^^^^^^^

=========================  =====================================================
Module                     Description
=========================  =====================================================
SLOTVerification.lean      SLOT predicate verification conditions and modes
Refinement.lean            VERIFIED mode refines CONSERVATIVE mode; simulation theorem
=========================  =====================================================

Specialized Proofs
^^^^^^^^^^^^^^^^^^

=========================  =====================================================
Module                     Description
=========================  =====================================================
SplicingResolution.lean    Splice site resolution and canonical donor/acceptor proofs
Mutagenesis.lean           Synonymous mutation theorems, GT/AG analysis, codon degeneracy
Soundness.lean             Main theorem module (entry point, re-exports all theorems)
=========================  =====================================================


Proof Status
------------

All modules are **FULLY PROVED** with 0 ``sorry``.

The Trusted Computing Base (TCB) consists of **3 scanner axioms + 4 SLOT
verification axioms**:

Scanner axioms (3):

1. ``SpliceSiteScanner.scanner_completeness`` — scanner finds all cryptic sites
2. ``SpliceSiteScanner.scanner_soundness`` — scanner only reports real sites
3. ``SpliceSiteScanner.borderline_completeness`` — borderline scanner completeness

SLOT verification axioms (4):

4. ``vc_imply_no_unexpected_tm`` — TMHMM VCs imply no unexpected TM domain
5. ``vc_imply_stable_folding`` — FoldX VCs imply stable folding
6. ``vc_imply_soluble_expression`` — ProteinSol VCs imply adequate solubility
7. ``vc_imply_low_immunogenicity`` — NetMHC VCs imply low immunogenicity

These are **parameters of the proof**, not gaps: the soundness theorem says
"ASSUMING the remaining scanners and tools are correct, the type system is
sound." This follows the standard approach in formal methods.


Building the Proofs
-------------------

Prerequisites: Lean4 (v4.30.0), Lake build system, Elan version manager.

.. code-block:: bash

   # Install elan
   curl -sSfL https://github.com/leanprover/elan/releases/latest/download/elan-x86_64-unknown-linux-gnu.tar.gz | tar xz
   ./elan-init -y --default-toolchain none
   source "$HOME/.elan/env"

   # Build all proofs
   cd proof/
   lake build

   # Fetch Mathlib cache
   lake exe cache get || true


Five-Valued Logic Extension
----------------------------

The ``FiveValued.lean`` module extends the three-valued logic with intermediate
verdicts (LIKELY_PASS, LIKELY_FAIL) while preserving soundness through a
conservative refinement mapping.

Key theorems:

- ``and_project_refines``: Projection of 5-valued AND equals 3-valued AND
- ``and_pass_pass``: ``five_valued_and PASS PASS = PASS``
- ``and_fail_absorb``: ``five_valued_and FAIL x = FAIL`` for all x
- ``and_comm``, ``and_assoc``: Commutativity and associativity


SLOT Property Semantics
-----------------------

Four SLOT predicates have non-vacuous semantic definitions:

=========================  =====================================================
Predicate                  Property Semantics
=========================  =====================================================
NoUnexpectedTMDomain       ∀ pos, tmHydrophobicFraction < threshold
StableFolding              predictedStabilityScore ≤ ddgThreshold
SolubleExpression          camSolScore ≥ minScore
LowImmunogenicity          ∀ pos, mhcBindingAffinity > maxScore
=========================  =====================================================

The remaining SLOT predicates use vacuous ``True`` semantics, awaiting
progressive strengthening.


CI Integration
--------------

The CI pipeline includes a ``proof-check`` job that:

1. Installs elan and Lean4 (v4.30.0)
2. Runs ``lake build`` in the ``proof/`` directory
3. Checks for ``sorry`` in proof files (must be sorry-free)
4. Runs as a separate job from Python tests
