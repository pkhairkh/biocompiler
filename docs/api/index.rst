BioCompiler API Reference
=========================

BioCompiler is a compiler framework for machine-verified gene design with
a comprehensive constraint system, unified CAI tables, hybrid optimization,
and deterministic type checking.

This section documents the public API for each top-level module.

.. toctree::
   :maxdepth: 2
   :caption: Core Modules

   biocompiler.optimizer
   biocompiler.biosecurity
   biocompiler.validation.protein_verification
   biocompiler.optimizer.objectives

.. toctree::
   :maxdepth: 2
   :caption: Constraint Modules

   biocompiler.optimizer.multigene
   biocompiler.sequence.pattern_enforcement
   biocompiler.sequence.iupac
   biocompiler.sequence.sliding_gc
   biocompiler.sequence.local_gc

.. toctree::
   :maxdepth: 2
   :caption: Design & Assembly

   biocompiler.optimizer.parts
   biocompiler.optimizer.assembly
   biocompiler.export.sbol_export
   biocompiler.export.sbol_import

.. toctree::
   :maxdepth: 2
   :caption: Integration & Export

   biocompiler.lims
   biocompiler.export.annotation
   biocompiler.export.genbank_roundtrip

.. toctree::
   :maxdepth: 2
   :caption: Scoring

   biocompiler.expression.tai
   biocompiler.validation.wetlab_validation
