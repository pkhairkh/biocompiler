BioCompiler API Reference
=========================

BioCompiler is a compiler framework for machine-verified gene design with
a comprehensive constraint system, unified CAI tables, hybrid optimization,
and deterministic type checking.

This section documents the public API for each top-level module.

.. toctree::
   :maxdepth: 2
   :caption: Core Modules

   biocompiler.optimization
   biocompiler.biosecurity
   biocompiler.protein_verification
   biocompiler.objectives

.. toctree::
   :maxdepth: 2
   :caption: Constraint Modules

   biocompiler.multigene
   biocompiler.pattern_enforcement
   biocompiler.iupac
   biocompiler.sliding_gc
   biocompiler.local_gc

.. toctree::
   :maxdepth: 2
   :caption: Design & Assembly

   biocompiler.parts
   biocompiler.assembly
   biocompiler.sbol_export
   biocompiler.sbol_import

.. toctree::
   :maxdepth: 2
   :caption: Integration & Export

   biocompiler.lims
   biocompiler.annotation
   biocompiler.genbank_roundtrip

.. toctree::
   :maxdepth: 2
   :caption: Scoring

   biocompiler.tai
   biocompiler.wetlab_validation
