Command-Line Interface
======================

BioCompiler provides a comprehensive CLI for certified gene optimization
and protein analysis. The CLI is accessible via the ``biocompiler`` command.


Commands
--------

optimize
^^^^^^^^

Optimize a protein sequence for a target organism:

.. code-block:: bash

   # Optimize a protein for human expression
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism Homo_sapiens

   # Optimize from a FASTA file
   biocompiler optimize --input gene.fasta --organism E_coli

   # Choose optimization strategy
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism human --strategy hybrid
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism human --strategy csp
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism human --strategy constraint_first

   # JSON output for pipeline integration
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism human --json

   # Verbose output with timing
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism human --verbose

**Key options:**

- ``--organism`` / ``--species``: Target organism (default: ``Homo_sapiens``)
- ``--strategy``: Optimization backend: ``hybrid`` (default), ``constraint_first``, or ``csp``
- ``--no-splice-check``: Skip eukaryotic splice-site constraints (for prokaryotes)
- ``--codon-pair-bias``: Optimize codon-pair bias during the run
- ``--gc-lo`` / ``--gc-hi``: GC content bounds (default: 0.30–0.70)
- ``--output``: Output FASTA file path
- ``--json``: Machine-readable JSON output
- ``--verbose``: Detailed optimization trace with timing
- ``--seed``: Deterministic random seed for reproducible optimization


batch
^^^^^

Batch-optimize proteins from a file:

.. code-block:: bash

   biocompiler batch proteins.txt --organism human --json

The input file should contain one protein per line, optionally with a name:

.. code-block:: text

   # Optional header/comments
   my_protein  MVLSPADKTNVKAAWGKVGA
   gfp         MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE


check
^^^^^

Read a FASTA file, evaluate all registered predicates, and print a certificate:

.. code-block:: bash

   # Check all predicates
   biocompiler check --input gene.fasta

   # List available predicates
   biocompiler check --list-predicates

   # Filter by category
   biocompiler check --input gene.fasta --category dna
   biocompiler check --input gene.fasta --category immunogenicity

   # Filter by specific predicates
   biocompiler check --input gene.fasta --predicates NoCrypticSplice,GCInRange


benchmark
^^^^^^^^^

Run built-in benchmarks:

.. code-block:: bash

   biocompiler benchmark --genes eGFP,mCherry,LacZ --organism E_coli


scan
^^^^

Scan a DNA sequence for features:

.. code-block:: bash

   biocompiler scan --input gene.fasta


serve
^^^^^

Start the REST API server:

.. code-block:: bash

   biocompiler serve --host 0.0.0.0 --port 8000


structure
^^^^^^^^^

Predict and assess protein structure using ESMFold:

.. code-block:: bash

   biocompiler structure --protein MVLSPADKTNVKAAWGKVGA


stability
^^^^^^^^^

Analyze protein stability using FoldX:

.. code-block:: bash

   biocompiler stability --protein MVLSPADKTNVKAAWGKVGA


solubility
^^^^^^^^^^

Analyze protein solubility using CamSol:

.. code-block:: bash

   biocompiler solubility --protein MVLSPADKTNVKAAWGKVGA


immunogenicity
^^^^^^^^^^^^^^

Analyze and reduce immunogenicity:

.. code-block:: bash

   biocompiler immunogenicity --protein MVLSPADKTNVKAAWGKVGA --organism human


assess
^^^^^^

Full protein assessment (structure + stability + solubility + immunogenicity):

.. code-block:: bash

   biocompiler assess --protein MVLSPADKTNVKAAWGKVGA --organism human


whatif
^^^^^^

Run what-if analysis on a protein sequence:

.. code-block:: bash

   biocompiler whatif --protein MVLSPADKTNVKAAWGKVGA --organism human


Predicate Categories
--------------------

Predicates can be filtered by category using the ``--category`` flag:

============  ================================================================
Category      Predicates
============  ================================================================
biosecurity   NoCrypticSplice, SpliceCorrect, NoCpGIsland, NoRestrictionSite,
              NoGTDinucleotide, NoCrypticPromoter, NoInstabilityMotif
dna           NoCrypticSplice, SpliceCorrect, GCInRange, SlidingGC,
              CodonAdapted, NoRestrictionSite, InFrame, NoInstabilityMotif,
              NoCpGIsland, NoStopCodons, NoGTDinucleotide, ValidCodingSeq,
              ConservationScore, CodonOptimality, NoCrypticPromoter,
              CoTranslationalFolding, NoUnexpectedTMDomain,
              mRNASecondaryStructure
codon_usage   CodonAdapted, CodonOptimality
structural    StructureConfidence, NoMisfoldingRisk, CorrectFoldTopology,
              NoUnexpectedInteraction
stability     StableFolding, NoDestabilizingMutation, DisulfideBondIntegrity,
              HydrophobicCoreQuality
solubility    SolubleExpression, NoAggregationProneRegion, ChargeComposition,
              NoLongHydrophobicStretch
immunogenicity LowImmunogenicity, NoStrongTCellEpitope, NoDominantBCellEpitope,
              PopulationCoverageSafe
============  ================================================================


Exit Codes
----------

======  ================================================================
Code    Meaning
======  ================================================================
0       Success
1       General error / invalid input
2       Biosecurity screening failure (enforce mode)
======  ================================================================
