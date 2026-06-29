Optimizer
=========

BioCompiler provides multiple optimization strategies for gene design, all
producing deterministic, reproducible results.


Optimization Strategies
-----------------------

Integrated Optimizer (Default)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The integrated optimizer (``integrated_optimize``) runs the full greedy
multi-step sequence — CAI maximization, restriction-site removal, ATTTA/T-run
cleanup, GC adjustment, cryptic-splice elimination, and reconciliation —
with CAI hill climbing. This is the default strategy invoked by
``optimize_sequence`` and produces ~14× faster optimization than DNAchisel
with higher CAI.

.. code-block:: python

   from biocompiler import optimize_sequence

   result = optimize_sequence(
       "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
       organism="Homo_sapiens",
       gc_lo=0.30,
       gc_hi=0.70,
   )
   print(f"CAI: {result.cai:.4f}")
   print(f"GC: {result.gc_content:.4f}")
   print(f"Sequence: {result.sequence}")


Constraint-First Strategy
^^^^^^^^^^^^^^^^^^^^^^^^^

Resolves all constraints first, then maximizes CAI. This is useful when
constraint satisfaction is the priority:

.. code-block:: python

   result = optimize_sequence(
       protein_seq,
       organism="Homo_sapiens",
       strategy="constraint_first",
   )


CSP Solver Strategy
^^^^^^^^^^^^^^^^^^^

Uses constraint satisfaction problem (CSP) solvers — OR-Tools CP-SAT —
for optimal codon assignment:

.. code-block:: python

   result = optimize_sequence(
       protein_seq,
       organism="E_coli",
       strategy="csp",
   )

Requires the ``solver`` extra (OR-Tools CP-SAT): ``pip install biocompiler[solver]``


BioOptimizer Class
------------------

The ``BioOptimizer`` class provides a lower-level interface for optimization
with fine-grained control:

.. code-block:: python

   from biocompiler.optimizer import BioOptimizer

   opt = BioOptimizer(
       organism="Homo_sapiens",
       enzymes=["EcoRI", "BamHI", "XhoI"],
       splice_low=3.0,
       splice_high=6.0,
       avoid_gt=True,
       organism_domain="eukaryote",
   )
   optimized_seq, pred_results, cert_text = opt.optimize(dna_sequence)


Batch Optimization
------------------

Optimize multiple proteins in a single call:

.. code-block:: python

   from biocompiler import batch_optimize

   proteins = ["MVLSPADKTNVKAAWGKVGA", "MSKGEELFTGVVPILVELDG"]
   results = batch_optimize(proteins, organism="Homo_sapiens")


Large Sequence Support
----------------------

For proteins of any length (including multi-kb sequences), the unified
``optimize_sequence()`` entry point automatically selects an appropriate
strategy. There is no separate large-sequence module — chunk-based
optimization, when needed, is handled internally by the integrated
optimizer:

.. code-block:: python

   from biocompiler import optimize_sequence

   # Works for any protein length; the integrated optimizer handles
   # internal chunking and reconciliation automatically.
   result = optimize_sequence(long_protein, organism="CHO_K1")


Incremental Constraint Checking
-------------------------------

The integrated optimizer maintains incremental constraint state internally
for O(1) constraint re-checking after codon changes. This is not a separate
public module — it is part of ``integrated_optimize`` and is used
automatically by ``optimize_sequence``:

.. code-block:: python

   from biocompiler.optimizer.integrated_optimizer import integrated_optimize

   # Incremental re-checking is handled inside integrated_optimize; callers
   # do not need to (and cannot) construct an IncrementalSequenceState.
   seq, notes, secis_positions = integrated_optimize(
       protein_seq, organism="Homo_sapiens",
   )


OptimizationResult
------------------

The ``optimize_sequence()`` function returns an ``OptimizationResult`` object:

=============  ================================================================
Field          Description
=============  ================================================================
sequence       Optimized DNA sequence
cai            Codon Adaptation Index
gc_content     GC fraction of the optimized sequence
certificate    Certificate level (GOLD, SILVER, BRONZE)
predicate_results  List of PredicateResult objects
certificate_text   Formatted certificate string
=============  ================================================================


Objective Functions
-------------------

BioCompiler supports multiple objective functions for optimization:

.. code-block:: python

   from biocompiler.optimizer.objectives import resolve_objective, OBJECTIVE_REGISTRY

   # Available objectives
   print(list(OBJECTIVE_REGISTRY.keys()))
   # ['cai', 'cai_gc_balanced', 'codon_pair', 'min_max_gc']

   # Use a specific objective
   result = optimize_sequence(
       protein_seq,
       organism="Homo_sapiens",
       objective="cai_gc_balanced",
   )


Performance
-----------

Benchmark results for E. coli GFP optimization:

==============  ===========  ===========
Metric          Before v10   After v10
==============  ===========  ===========
CAI             0.67         0.999
Time            20 ms        2 ms
==============  ===========  ===========
