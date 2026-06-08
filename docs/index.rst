BioCompiler Documentation
=========================

BioCompiler is a compiler framework for machine-verified gene design with a
comprehensive 33-predicate system, unified CAI tables, hybrid optimization,
and deterministic type checking.

The compilation pipeline follows a deterministic flow:

  **Scanner** → **NDFST Splicing** → **Translation** → **Type Check** → **Certificate** → **Verify**

All computation is deterministic: the same input always produces identical output.


.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   cli


.. toctree::
   :maxdepth: 2
   :caption: Core Concepts

   predicates
   optimizer


.. toctree::
   :maxdepth: 2
   :caption: Safety & Verification

   biosecurity
   lean4


.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
