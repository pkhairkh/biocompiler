Installation
============

BioCompiler requires Python 3.10 or later.


Quick Install
-------------

Install from PyPI (when available):

.. code-block:: bash

   pip install biocompiler

Or install from source:

.. code-block:: bash

   git clone https://github.com/pkhairkh/biocompiler.git
   cd biocompiler
   pip install .


Optional Dependencies
---------------------

BioCompiler has several optional dependency groups that extend functionality:

.. code-block:: bash

   # CSP solver backends (OR-Tools + Z3)
   pip install biocompiler[solver]

   # Z3 solver only
   pip install biocompiler[optimizer]

   # ViennaRNA mRNA structure prediction
   pip install biocompiler[viennarna]

   # MHCflurry immunogenicity prediction
   pip install biocompiler[mhcflurry]

   # BioPython integration
   pip install biocompiler[biopython]

   # Jupyter notebook integration
   pip install biocompiler[jupyter]

   # Documentation building
   pip install biocompiler[docs]

   # All optional dependencies
   pip install biocompiler[all]

   # Development tools (linting, testing, docs)
   pip install biocompiler[dev]


Verifying the Installation
---------------------------

After installation, verify the CLI is accessible:

.. code-block:: bash

   biocompiler --help
   biocompiler --version

You can also verify the Python API:

.. code-block:: python

   import biocompiler
   print(biocompiler.__version__)

   # Quick optimization
   from biocompiler import optimize_sequence
   result = optimize_sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", organism="Homo_sapiens")
   print(f"CAI: {result.cai:.4f}")


Docker Installation
-------------------

A Dockerfile is provided for containerized usage:

.. code-block:: bash

   docker build -t biocompiler .
   docker run biocompiler biocompiler --version


Building Documentation
----------------------

To build the documentation locally:

.. code-block:: bash

   pip install biocompiler[docs]
   cd docs/
   make html

The built HTML documentation will be in ``docs/_build/html/``.
