BioCompiler Documentation
=========================

BioCompiler is a compiler framework for machine-verified gene design with a
comprehensive 43-predicate system, unified CAI tables, hybrid optimization,
and deterministic type checking.

The compilation pipeline follows a deterministic flow:

  **Scanner** → **NDFST Splicing** → **Translation** → **Type Check** → **Certificate** → **Verify**

All computation is deterministic: the same input always produces identical output.

.. tip::

   New to BioCompiler? Start with the **User Guide** and **Examples** — they
   cover the four most common tasks (optimize, verify, screen, YAML compile)
   with copy-paste-runnable code. The technical sections below go deeper.


Practical Guides (Markdown)
---------------------------

These standalone guides are written in Markdown for easy reading on GitHub
or any Markdown viewer. They are not part of the Sphinx toctree (which is
RestructuredText-only); they are linked here so the index points at them.

- `User Guide <USER_GUIDE.md>`_ — installation, the four core use cases
  (optimize, verify, screen, YAML compile), the full organism table, the
  CLI reference, and troubleshooting.
- `Examples <EXAMPLES.md>`_ — five runnable, end-to-end scripts: HBB for
  human, insulin for E. coli, gene therapy construct verification, batch
  optimization, and a head-to-head DNAchisel comparison.


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
