Optimizer
=========

BioCompiler provides multiple optimization strategies for gene design, all
producing deterministic, reproducible results.


Optimization Strategies
-----------------------

Hybrid Optimizer (Default)
^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``HybridOptimizer`` combines greedy initialization with priority-based
constraint satisfaction and CAI hill climbing. This is the default strategy
and provides 3–5× speed improvement over the legacy pipeline with higher CAI.

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

Uses constraint satisfaction problem (CSP) solvers — either OR-Tools CP-SAT
or Z3 SMT — for optimal codon assignment:

.. code-block:: python

   result = optimize_sequence(
       protein_seq,
       organism="E_coli",
       strategy="csp",
   )

Requires the ``solver`` extra: ``pip install biocompiler[solver]``


BioOptimizer Class
------------------

The ``BioOptimizer`` class provides a lower-level interface for optimization
with fine-grained control:

.. code-block:: python

   from biocompiler.optimizer import BioOptimizer

   opt = BioOptimizer(
       species="Homo_sapiens",
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

For proteins longer than 10 kb, BioCompiler uses chunk-based optimization:

.. code-block:: python

   from biocompiler.large_sequence import optimize_large_sequence

   result = optimize_large_sequence(long_protein, organism="CHO_K1")


Incremental Constraint Checking
-------------------------------

BioCompiler uses ``IncrementalSequenceState`` for O(1) constraint re-checking
after codon changes, providing 2–2000× faster constraint evaluation:

.. code-block:: python

   from biocompiler.incremental import IncrementalSequenceState

   state = IncrementalSequenceState(dna_seq, organism="Homo_sapiens")
   # Swap a codon and re-check only affected constraints
   state.swap_codon(position, new_codon)


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

   from biocompiler.objectives import resolve_objective, OBJECTIVE_REGISTRY

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
