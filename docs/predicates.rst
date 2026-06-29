Predicate System
================

BioCompiler implements a comprehensive **43-predicate** type system for
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


43 Predicates
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

Extended Diagnostic Predicates (15)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======  =================================  ==========================================
No.     Predicate                          Description
======  =================================  ==========================================
29      SpliceCorrect                      Splice correctness with MaxEntScan
30      GCInRange                          Global GC content in valid range
31      CodonAdapted                       CAI threshold compliance
32      InFrame                            Reading frame preservation
33      NoInstabilityMotif                 AU-rich element avoidance (ARE Class I/II/III)
34      NoBlastMatches                     Biosecurity BLAST match avoidance
35      PrimerCompatibility                Primer design compatibility
36      NoCrypticORF                       Out-of-frame open reading frames
37      NoRQCTrigger                       Ribosome quality control triggers
38      NoAluRepeat                        Alu/SINE repetitive elements
39      NoMiRNABindingSite                 miRNA seed match avoidance (multi-organism)
40      NoM6ASite                          m6A RRACH motif avoidance
41      NoPolyASignal                      Premature polyadenylation signals
42      NucleosideModificationGuidance     Nucleoside modification recommendations
43      SlidingGC                          Sliding-window GC constraint
======  =================================  ==========================================


NoMiRNABindingSite — Detailed Reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``NoMiRNABindingSite`` predicate (predicate #39) scans DNA coding
sequences for seed matches against known miRNA families.  miRNA binding
to mRNA coding regions can trigger translational repression, mRNA
destabilization, or off-target silencing.

**Parameters:**

===========  ========  ==========================================================
Parameter    Default   Description
===========  ========  ==========================================================
organism     required  Target organism (``"Homo_sapiens"``, ``"Mus_musculus"``,
                       ``"CHO_K1"``, etc.).  Determines the seed database.
tissue       None      Optional tissue filter.  When specified, only seeds
                       expressed in the given tissue *or* ubiquitously
                       expressed seeds are checked.  Valid values include:
                       ``"liver"``, ``"muscle"``, ``"immune"``, ``"neural"``,
                       ``"blood"``, ``"endothelial"``, ``"kidney"``,
                       ``"colon"``, ``"fibrotic"``, ``"epithelial"``,
                       ``"lymphoid"``, ``"smooth_muscle"``,
                       ``"hematopoietic"``, ``"tumor_suppressor"``.
min_seed_match  7      Minimum seed match length (6 or 7).
===========  ========  ==========================================================

**Multi-Organism Support:**

The predicate now supports organism-specific miRNA seed databases:

- **Human** (*Homo sapiens*): 30 miRNA seeds from miRBase v22 + Ludwig et al. 2016
- **Mouse** (*Mus musculus*): 22 miRNA seeds (conserved + mouse-specific)
- **CHO-K1** (Chinese Hamster Ovary): 10 miRNA seeds (limited dataset, primarily
  human orthologs)
- **Prokaryotes** (*E. coli*, *S. cerevisiae*): No miRNA checks (miRNA machinery
  is absent)

For organisms without a dedicated seed database, the human seed database is used
as a fallback (closest available).  This can be overridden by providing a custom
seed set via the ``get_mirna_seeds()`` API.

**Tissue Filtering:**

When the ``tissue`` parameter is provided, only seeds expressed in that tissue
or marked as ``"ubiquitous"`` are considered.  This dramatically reduces false
positives for tissue-specific therapeutic applications (e.g., liver-targeted
mRNA vaccines should primarily avoid liver-expressed miRNA sites).

.. code-block:: python

   from biocompiler.type_system import check_no_mirna_binding_site, get_mirna_seeds

   # Check all human miRNA seeds (default)
   result = check_no_mirna_binding_site(seq, organism="Homo_sapiens")

   # Check only liver-relevant seeds
   result = check_no_mirna_binding_site(seq, organism="Homo_sapiens", tissue="liver")

   # Browse available seeds for mouse
   seeds = get_mirna_seeds("Mus_musculus")

   # Browse liver-specific seeds for human
   seeds = get_mirna_seeds("Homo_sapiens", tissue="liver")

**Seed Match Types:**

Match type scoring follows the canonical miRNA target classification:

=========  ==============================  =====
Type       Description                     Score
=========  ==============================  =====
8mer       Positions 2–8 match + A at pos1  1.0
7mer-m8    Positions 2–8 match              0.9
7mer-A1    Positions 2–7 match + A at pos1  0.85
6mer       Positions 2–7 match              0.7
=========  ==============================  =====

**Verdict Logic:**

- **FAIL**: 8mer or 7mer-m8 match (score ≥ 0.9) for any top-10 abundantly
  expressed miRNA (tier 1)
- **UNCERTAIN**: 7mer-A1 or 6mer match, or match to lower-abundance miRNA
  (tier 2–3)
- **PASS**: No significant miRNA seed matches

**Elimination:**

The optimizer provides ``eliminate_mirna_binding_sites()`` to disrupt miRNA
seed matches via synonymous codon substitution:

.. code-block:: python

   from biocompiler.optimizer import eliminate_mirna_binding_sites

   new_seq, warnings = eliminate_mirna_binding_sites(
       sequence=seq,
       aas=amino_acids,
       sorted_codons=codon_table,
       usage=cai_usage,
       organism="Homo_sapiens",
       tissue="liver",
   )


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

   from biocompiler.organisms.config import is_eukaryotic_organism, auto_detect_organism_domain

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
