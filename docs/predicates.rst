Predicate System
================

BioCompiler implements a comprehensive **28-predicate** type system for
certified gene optimization. Each predicate evaluates a specific biological
constraint and returns a three-valued verdict: **PASS**, **FAIL**, or
**UNCERTAIN**.

The predicate system is the core of BioCompiler's type-theoretic approach to
gene design: "well-typed genes don't go wrong." If all predicates return PASS,
the optimized sequence is guaranteed to satisfy all specified constraints.


Three-Valued Logic
------------------

Each predicate returns a verdict from a three-valued logic:

========  =================================================================
Verdict   Meaning
========  =================================================================
PASS      The constraint is definitively satisfied.
FAIL      The constraint is definitively violated.
UNCERTAIN Insufficient evidence to determine pass or fail.
========  =================================================================

Combined verdicts follow Kleene-style conjunction (weakest-link principle):
the overall result is determined by the worst individual verdict.

PASS ⊓ PASS = PASS
PASS ⊓ UNCERTAIN = UNCERTAIN
PASS ⊓ FAIL = FAIL
UNCERTAIN ⊓ FAIL = FAIL


Five-Valued Logic Extension
----------------------------

The runtime extends the three-valued logic with two intermediate verdicts
for cases where a tool is available but the result is not formally verified:

=============  =====  ======================================================
Verdict        Score  Meaning
=============  =====  ======================================================
PASS           1.0    Formally verified pass
LIKELY_PASS    0.75   High confidence but not formally verified
UNCERTAIN      0.5    Insufficient evidence
LIKELY_FAIL    0.25   High confidence of failure
FAIL           0.0    Formally verified failure
=============  =====  ======================================================

The five-valued logic refines conservatively to the three-valued model:
LIKELY_PASS → UNCERTAIN and LIKELY_FAIL → FAIL. This ensures that the
Lean4 soundness proofs remain valid when applied to five-valued verdicts.


28 Predicates
-------------

DNA-Level Predicates (12)
^^^^^^^^^^^^^^^^^^^^^^^^^

======  ============================  ==============================================
No.     Predicate                     Description
======  ============================  ==============================================
1       NoStopCodons                  No internal stop codons in reading frame
2       NoCrypticSplice               No cryptic splice sites (MaxEntScan scoring)
3       NoCpGIsland                   CpG island avoidance (Obs/Exp ratio check)
4       NoRestrictionSite             No restriction enzyme recognition sites
5       NoGTDinucleotide              GT dinucleotide avoidance (splice donor mimic)
6       ValidCodingSeq                In-frame, valid codons only
7       ConservationScore             BLOSUM62-based amino acid conservation
8       CodonOptimality               CAI-based codon quality
9       NoCrypticPromoter             Cryptic promoter avoidance
10      NoUnexpectedTMDomain          Unexpected transmembrane domain detection
11      mRNASecondaryStructure        mRNA secondary structure around RBS
12      CoTranslationalFolding        Co-translational folding pause-site preservation
======  ============================  ==============================================

Structure Predicates (4)
^^^^^^^^^^^^^^^^^^^^^^^^

======  ============================  ==============================================
No.     Predicate                     Description
======  ============================  ==============================================
13      StructureConfidence           ESMFold structure quality confidence (pLDDT)
14      NoMisfoldingRisk              Misfolding risk indicators
15      CorrectFoldTopology           Fold topology validation
16      NoUnexpectedInteraction       Unwanted protein-protein interactions
======  ============================  ==============================================

Stability Predicates (4)
^^^^^^^^^^^^^^^^^^^^^^^^

======  ============================  ==============================================
No.     Predicate                     Description
======  ============================  ==============================================
17      StableFolding                 Thermodynamic stability (ΔG from FoldX)
18      NoDestabilizingMutation       No high-ΔΔG mutations
19      DisulfideBondIntegrity        Cysteine pairing check
20      HydrophobicCoreQuality        Hydrophobic core composition
======  ============================  ==============================================

Solubility Predicates (4)
^^^^^^^^^^^^^^^^^^^^^^^^^

======  ============================  ==============================================
No.     Predicate                     Description
======  ============================  ==============================================
21      SolubleExpression             CamSol solubility score
22      NoAggregationProneRegion      Aggregation-prone region detection
23      ChargeComposition             Charge balance and isoelectric point
24      NoLongHydrophobicStretch      Long hydrophobic stretch detection
======  ============================  ==============================================

Immunogenicity Predicates (4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======  ============================  ==============================================
No.     Predicate                     Description
======  ============================  ==============================================
25      LowImmunogenicity             Overall immunogenicity score
26      NoStrongTCellEpitope          MHC binding epitope detection
27      NoDominantBCellEpitope        B-cell epitope coverage
28      PopulationCoverageSafe        MHC allele population coverage
======  ============================  ==============================================


Predicate Registry
------------------

All predicates are registered in a central ``PredicateRegistry`` that supports
runtime lookup, evaluation, and extensibility:

.. code-block:: python

   from biocompiler.type_system import registry

   # List all registered predicates
   names = registry.names()

   # Evaluate a specific predicate
   result = registry.evaluate("NoCrypticSplice", sequence=dna_seq, organism="Homo_sapiens")

   # Evaluate all predicates at once
   from biocompiler import evaluate_all_predicates
   results = evaluate_all_predicates(seq=dna_seq, organism="Homo_sapiens")


Organism Awareness
------------------

Many predicates are organism-aware and adjust their behavior accordingly:

- **Prokaryotes** (e.g., *E. coli*): Splice-site checks, CpG island checks,
  and GT dinucleotide checks are skipped since these are eukaryote-specific
  concerns.
- **Eukaryotes** (e.g., *H. sapiens*): All predicates are active with
  organism-specific thresholds.

This is handled automatically by the ``organism_config`` module:

.. code-block:: python

   from biocompiler.organism_config import is_eukaryotic_organism, auto_detect_organism_domain

   is_euk = is_eukaryotic_organism("E_coli")  # False
   domain = auto_detect_organism_domain("Saccharomyces_cerevisiae")  # "eukaryote"


SLOT Predicates
---------------

Some predicates (StructureConfidence, StableFolding, SolubleExpression,
LowImmunogenicity) rely on external tools (ESMFold, FoldX, ProteinSol,
NetMHC). These SLOT (Safe Lazy Oracles with Trust) predicates operate in
three modes:

=============  ================================================================
Mode           Behavior
=============  ================================================================
conservative   Returns UNCERTAIN when the tool is unavailable
verified       Requires tool output; fails if tool unavailable
permissive     Uses tool output when available, otherwise UNCERTAIN
=============  ================================================================
