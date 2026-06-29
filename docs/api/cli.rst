CLI Reference
=============

BioCompiler provides a comprehensive command-line interface for certified gene
optimization and protein analysis.

.. module:: biocompiler.cli
   :synopsis: Command-line interface for gene optimization and protein analysis.

Basic Usage
-----------

::

    biocompiler <command> [options]

Use ``biocompiler --help`` to see all commands, and ``biocompiler <command> --help``
for command-specific options.

Global Options
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Option
     - Description
   * - ``--version``
     - Show BioCompiler version and exit.
   * - ``--help``
     - Show help message and exit.
   * - ``--verbose``
     - Enable verbose output with timing information.
   * - ``--json``
     - Output results in machine-readable JSON format.

Commands
--------

optimize
^^^^^^^^

Optimize a protein sequence for a target organism.

Synopsis::

    biocompiler optimize PROTEIN [options]
    biocompiler optimize --input FASTA_FILE [options]

The ``optimize`` command accepts either a protein sequence string as a
positional argument, or a FASTA file via ``--input``. It back-translates
the protein to DNA using organism-specific codon tables and optimizes
the resulting sequence to satisfy all constraint predicates.

Input Options
"""""""""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``PROTEIN``
     - *(positional)*
     - Protein sequence (single-letter amino acid codes).
   * - ``--input PATH``
     - *(none)*
     - Input FASTA file (legacy mode). Mutually exclusive with positional PROTEIN.
   * - ``--organism NAME``
     - ``Homo_sapiens``
     - Target organism. Accepts canonical names, short keys, abbreviated binomials, or display names.
   * - ``--species NAME``
     - *(none)*
     - Deprecated alias for ``--organism``.

Optimization Strategy
"""""""""""""""""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--strategy``
     - ``hybrid``
     - Optimization backend: ``hybrid`` (greedy + hill climbing), ``constraint_first`` (constraints before CAI), or ``csp`` (constraint satisfaction problem solver via OR-Tools/Z3).
   * - ``--no-splice-check``
     - *(off)*
     - Skip eukaryotic splice-site constraints. For prokaryotic targets.
   * - ``--codon-pair-bias``
     - *(off)*
     - Optimize codon-pair bias during the run.
   * - ``--organism-domain``
     - ``auto``
     - Force organism domain: ``auto``, ``eukaryote``, or ``prokaryote``.

Constraint Parameters
"""""""""""""""""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--gc-lo``
     - 0.30
     - Minimum GC content fraction.
   * - ``--gc-hi``
     - 0.70
     - Maximum GC content fraction.
   * - ``--enzymes ENZYMES``
     - *(none)*
     - Comma-separated list of restriction enzymes to avoid (e.g. ``EcoRI,BamHI,XhoI``).
   * - ``--splice-low``
     - 3.0
     - Minimum MaxEntScan score for splice site detection.
   * - ``--splice-high``
     - 6.0
     - Maximum MaxEntScan score for cryptic splice sites.
   * - ``--avoid-gt``
     - ``true``
     - Avoid GT dinucleotides that can create donor splice sites.

Immunogenicity Options
""""""""""""""""""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--source-organism NAME``
     - *(none)*
     - Organism the protein originates from. Used for immunogenicity self-protein detection.
   * - ``--therapeutic``
     - *(off)*
     - Apply stricter immunogenicity thresholds for therapeutic proteins.

Output Options
""""""""""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--output PATH``
     - *(none)*
     - Output file path. Default extension depends on ``--format``.
   * - ``--format FORMAT``
     - ``fasta``
     - Output format: ``fasta``, ``genbank``, or ``sbol3``.
   * - ``--certificate PATH``
     - *(auto)*
     - Certificate output file path (legacy FASTA mode only).
   * - ``--json``
     - *(off)*
     - Output results as JSON.
   * - ``--provenance``
     - *(off)*
     - Track and save provenance data (legacy FASTA mode only).
   * - ``--biosecurity-report``
     - *(off)*
     - Include biosecurity screening report.
   * - ``--seed N``
     - *(none)*
     - Set random seed for reproducible optimization.

Examples
""""""""

Optimize eGFP for E. coli::

    biocompiler optimize MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE \
        --organism ecoli --json

Optimize from FASTA with GenBank output::

    biocompiler optimize --input gfp.fasta \
        --organism Homo_sapiens --format genbank \
        --output gfp_optimized.gb

Optimize with CSP solver and codon-pair bias::

    biocompiler optimize MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE \
        --organism ecoli --strategy csp --codon-pair-bias

Optimize a therapeutic protein with source organism::

    biocompiler optimize MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR \
        --organism human --source-organism ecoli --therapeutic \
        --output hbb_optimized.fasta

batch
^^^^^

Batch-optimize multiple proteins from a text file.

Synopsis::

    biocompiler batch PROTEINS_FILE [options]

The input file should contain one protein per line. Lines starting with
``#`` are treated as comments. Each line can be either:

- A plain protein sequence: ``MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE``
- A named protein: ``GFP MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE``

Options
"""""""

Supports the same optimization options as ``optimize`` (``--organism``,
``--strategy``, ``--gc-lo``, ``--gc-hi``, ``--codon-pair-bias``,
``--source-organism``, ``--therapeutic``, ``--seed``, ``--verbose``,
``--json``, ``--output``).

Example
"""""""

::

    # proteins.txt
    GFP  MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE
    Insulin MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVE
    # This is a comment

    biocompiler batch proteins.txt --organism ecoli --json --output results.fasta

check
^^^^^

Type-check a DNA sequence against all registered predicates without optimizing.

Synopsis::

    biocompiler check --input FASTA_FILE [options]

Evaluates every registered predicate (currently 12+), producing a certificate
that records which constraints pass and which fail.

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--input PATH``
     - *(required)*
     - Input FASTA file.
   * - ``--organism NAME``
     - ``Homo_sapiens``
     - Target organism.
   * - ``--species NAME``
     - *(none)*
     - Deprecated alias for ``--organism``.
   * - ``--enzymes ENZYMES``
     - *(none)*
     - Comma-separated restriction enzymes.
   * - ``--organism-domain``
     - ``auto``
     - Force organism domain.
   * - ``--source-organism``
     - *(none)*
     - Source organism for immunogenicity.
   * - ``--therapeutic``
     - *(off)*
     - Therapeutic mode.
   * - ``--strict-mode``
     - *(off)*
     - Strict mode.
   * - ``--json``
     - *(off)*
     - JSON output.

Example
"""""""

::

    biocompiler check --input hbb.fasta --organism Homo_sapiens --json

scan
^^^^

Scan a DNA sequence for motifs, restriction sites, and open reading frames.

Synopsis::

    biocompiler scan SEQUENCE [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``SEQUENCE``
     - *(positional)*
     - DNA sequence to scan.
   * - ``--enzymes ENZYMES``
     - *(none)*
     - Comma-separated restriction enzymes to scan for.
   * - ``--find-orfs``
     - *(off)*
     - Detect open reading frames.

Example
"""""""

::

    biocompiler scan ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGG \
        --enzymes EcoRI,BamHI --find-orfs

benchmark
^^^^^^^^^

Run built-in benchmarks or named gene sets.

Synopsis::

    biocompiler benchmark [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--gene-set NAME``
     - *(none)*
     - Named gene set: ``REFERENCE_GENES``, ``HUMAN_THERAPEUTIC``, or ``GENE_PANEL``.
   * - ``--list-gene-sets``
     - *(off)*
     - List available gene sets and exit.
   * - ``--enzymes ENZYMES``
     - *(none)*
     - Restriction enzymes.
   * - ``--splice-low``
     - 3.0
     - Splice site detection threshold.
   * - ``--splice-high``
     - 6.0
     - Cryptic splice threshold.
   * - ``--seed N``
     - *(none)*
     - Random seed for reproducibility.
   * - ``--output PATH``
     - *(none)*
     - Output file (``.json`` for JSON, else CSV).

Examples
""""""""

::

    biocompiler benchmark
    biocompiler benchmark --gene-set REFERENCE_GENES --output results.json
    biocompiler benchmark --list-gene-sets

serve
^^^^^

Start the REST API server.

Synopsis::

    biocompiler serve [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--host HOST``
     - ``0.0.0.0``
     - Bind host.
   * - ``--port PORT``
     - ``8000``
     - Bind port.
   * - ``--no-auth``
     - *(off)*
     - Disable API authentication (dangerous — use only for local development).
   * - ``--reload``
     - *(off)*
     - Enable auto-reload for development.

Example
"""""""

::

    biocompiler serve --port 8080 --no-auth
    biocompiler serve --host 127.0.0.1 --port 8000

structure
^^^^^^^^^

Predict and assess protein structure using ESMFold.

Synopsis::

    biocompiler structure --protein SEQUENCE [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--protein SEQUENCE``
     - *(required)*
     - Protein sequence.
   * - ``--organism NAME``
     - ``Homo_sapiens``
     - Target organism.
   * - ``--output PATH``
     - *(none)*
     - Output PDB file path.

stability
^^^^^^^^^

Analyze protein stability.

Synopsis::

    biocompiler stability --protein SEQUENCE [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--protein SEQUENCE``
     - *(required)*
     - Protein sequence.
   * - ``--organism NAME``
     - ``Homo_sapiens``
     - Target organism.
   * - ``--pdb PATH``
     - *(none)*
     - Optional PDB structure file.
   * - ``--json``
     - *(off)*
     - JSON output.

solubility
^^^^^^^^^^

Analyze protein solubility.

Synopsis::

    biocompiler solubility --protein SEQUENCE [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--protein SEQUENCE``
     - *(required)*
     - Protein sequence.
   * - ``--pdb PATH``
     - *(none)*
     - Optional PDB structure file.

immunogenicity
^^^^^^^^^^^^^^

Analyze and reduce protein immunogenicity.

Synopsis::

    biocompiler immunogenicity --protein SEQUENCE [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--protein SEQUENCE``
     - *(required)*
     - Protein sequence.
   * - ``--organism NAME``
     - ``Homo_sapiens``
     - Host organism.
   * - ``--source-organism``
     - *(none)*
     - Source organism for self-protein detection.
   * - ``--therapeutic``
     - *(off)*
     - Therapeutic mode.
   * - ``--deimmunize``
     - *(off)*
     - Run deimmunization.
   * - ``--target-score``
     - 0.3
     - Target immunogenicity score for deimmunization.
   * - ``--max-mutations``
     - 10
     - Maximum mutations for deimmunization.

assess
^^^^^^

Run a full protein assessment combining structure, stability, solubility,
and immunogenicity analyses.

Synopsis::

    biocompiler assess --protein SEQUENCE [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--protein SEQUENCE``
     - *(required)*
     - Protein sequence.
   * - ``--organism NAME``
     - ``Homo_sapiens``
     - Target organism.
   * - ``--pdb PATH``
     - *(none)*
     - Optional PDB structure file.
   * - ``--no-structure``
     - *(off)*
     - Skip structure prediction.
   * - ``--no-stability``
     - *(off)*
     - Skip stability analysis.
   * - ``--no-solubility``
     - *(off)*
     - Skip solubility analysis.
   * - ``--no-immunogenicity``
     - *(off)*
     - Skip immunogenicity analysis.

validate-cai
^^^^^^^^^^^^

Validate CAI computation against published ground-truth values.

Synopsis::

    biocompiler validate-cai [options]

This command runs the CAI validation benchmark, comparing BioCompiler's
CAI values against published reference values from Sharp & Li (1987) and
other sources.

validate-maxentscan
^^^^^^^^^^^^^^^^^^^

Validate MaxEntScan scores against published values.

Synopsis::

    biocompiler validate-maxentscan [options]

whatif
^^^^^^

Run what-if analysis on a protein sequence.

Synopsis::

    biocompiler whatif --protein SEQUENCE [options]

Explores how changes to optimization parameters would affect the result
without actually running a full optimization.

explain
^^^^^^^

Explain the optimization decisions for a specific position.

Synopsis::

    biocompiler explain --provenance-id ID --position N [options]

report
^^^^^^

Generate a provenance report for an optimization run.

Synopsis::

    biocompiler report --provenance-id ID [options]

Options
"""""""

.. list-table::
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--provenance-id ID``
     - *(required)*
     - Provenance trail ID.
   * - ``--format FORMAT``
     - ``text``
     - Report format: ``text``, ``markdown``, or ``json``.
   * - ``--output PATH``
     - *(none)*
     - Output file path.

Organism Name Resolution
-------------------------

All CLI commands that accept an ``--organism`` parameter support multiple
input forms that are automatically resolved to the canonical name:

.. list-table::
   :header-rows: 1

   * - Input Form
     - Example
     - Resolved To
   * - Short key
     - ``ecoli``
     - ``Escherichia_coli``
   * - Abbreviated binomial
     - ``E_coli``
     - ``Escherichia_coli``
   * - Display name
     - ``E. coli``
     - ``Escherichia_coli``
   * - Canonical name
     - ``Escherichia_coli``
     - ``Escherichia_coli``

Supported organisms:

.. list-table::
   :header-rows: 1

   * - Organism
     - Canonical Name
     - Short Key
     - Domain
   * - *Escherichia coli*
     - ``Escherichia_coli``
     - ``ecoli``
     - Prokaryote
   * - *Homo sapiens*
     - ``Homo_sapiens``
     - ``human``
     - Eukaryote
   * - *Mus musculus*
     - ``Mus_musculus``
     - ``mouse``
     - Eukaryote
   * - *Saccharomyces cerevisiae*
     - ``Saccharomyces_cerevisiae``
     - ``yeast``
     - Eukaryote
   * - CHO-K1
     - ``CHO_K1``
     - ``cho``
     - Eukaryote

Output Formats
--------------

The CLI supports multiple output formats:

- **Text** (default): Human-readable output with colour support on TTY.
- **JSON** (``--json``): Machine-readable JSON output for scripting and pipelines.
- **FASTA**: DNA sequence in FASTA format (``--format fasta``).
- **GenBank**: GenBank flat-file format (``--format genbank``).
- **SBOL3**: SBOL3 XML or JSON-LD (``--format sbol3``).

Colour Support
--------------

The CLI automatically enables ANSI colour output when running on a TTY.
When piped or redirected, colour codes are suppressed automatically.

Public API
----------

.. autofunction:: main

.. autofunction:: build_parser

.. autofunction:: colorize

.. autofunction:: cmd_optimize

.. autofunction:: cmd_batch

.. autofunction:: cmd_check

.. autofunction:: cmd_benchmark

.. autofunction:: cmd_scan

.. autofunction:: cmd_structure

.. autofunction:: cmd_stability

.. autofunction:: cmd_solubility

.. autofunction:: cmd_immunogenicity

.. autofunction:: cmd_assess

.. autofunction:: cmd_validate_cai

.. autofunction:: cmd_validate_maxentscan

.. autofunction:: cmd_whatif

.. autofunction:: cmd_explain

.. autofunction:: cmd_report
