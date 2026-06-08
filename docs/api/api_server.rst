REST API Server Reference
=========================

BioCompiler provides a production-grade FastAPI REST API for integration with
other bioinformatics tools, pipelines, and web applications.

.. module:: biocompiler.api
   :synopsis: FastAPI REST API for gene optimization and protein analysis.

Starting the Server
-------------------

Using the CLI::

    biocompiler serve --port 8000

Using uvicorn directly::

    uvicorn biocompiler.api:app --host 0.0.0.0 --port 8000

Environment Variables
---------------------

The API server behaviour can be configured via environment variables:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Variable
     - Default
     - Description
   * - ``BIOCOMPILER_API_KEY``
     - *(auto-generated)*
     - API key for authentication. Set to a specific string for production, ``disabled`` to turn off auth, or leave unset to auto-generate and persist to ``~/.biocompiler/api_key``.
   * - ``BIOCOMPILER_API_KEYS``
     - *(empty)*
     - Comma-separated list of API keys for key rotation. Takes precedence over ``BIOCOMPILER_API_KEY``.
   * - ``BIOCOMPILER_AUTH_MODE``
     - ``required``
     - Authentication enforcement mode: ``required``, ``optional``, or ``disabled``.
   * - ``BIOCOMPILER_RATE_LIMIT``
     - ``60``
     - Maximum requests per minute per client.
   * - ``BIOCOMPILER_RATE_LIMIT_DB``
     - ``~/.biocompiler/rate_limits.db``
     - Path to the SQLite database used for persistent rate-limit tracking.
   * - ``BIOCOMPILER_MAX_PROTEIN_LENGTH``
     - ``10000``
     - Maximum protein sequence length in amino acids.
   * - ``BIOCOMPILER_MAX_BATCH_SIZE``
     - ``50``
     - Maximum number of items per batch request.
   * - ``BIOCOMPILER_MAX_REQUEST_SIZE``
     - ``10000000``
     - Maximum request body size in bytes.
   * - ``BIOCOMPILER_OPTIMIZE_TIMEOUT``
     - ``300``
     - Optimization timeout in seconds.
   * - ``BIOCOMPILER_BATCH_ITEM_TIMEOUT``
     - ``30``
     - Per-item timeout in batch requests (seconds).
   * - ``BIOCOMPILER_PROVENANCE_DIR``
     - *(none)*
     - Directory for persistent provenance storage.

Authentication
--------------

The API uses key-based authentication via the ``X-API-Key`` header.

Auth Modes
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Mode
     - Behaviour
   * - ``required``
     - *(default)* Unauthenticated requests receive HTTP 401. This is the safe default for production.
   * - ``optional``
     - Unauthenticated requests are allowed but receive an ``X-Auth-Warning`` response header. Useful during migration.
   * - ``disabled``
     - Authentication is completely disabled. **DANGEROUS** for production. Use only with ``BIOCOMPILER_API_KEY=disabled`` or ``--no-auth`` CLI flag.

Key Management
^^^^^^^^^^^^^^

- **Auto-generated key**: If ``BIOCOMPILER_API_KEY`` is not set, a random 64-character hex key is generated, printed to the console on first startup, and persisted to ``~/.biocompiler/api_key`` for reuse.
- **Explicit key**: Set ``BIOCOMPILER_API_KEY=your-secret-key`` for production deployments.
- **Key rotation**: Set ``BIOCOMPILER_API_KEYS=key1,key2,key3`` for zero-downtime key rotation. The server accepts any key in the list.
- **Constant-time comparison**: All key comparisons use ``secrets.compare_digest`` to prevent timing attacks.

Rate Limiting
-------------

The API enforces per-client rate limiting using a sliding-window algorithm backed
by SQLite.

- **Default**: 60 requests per minute per client.
- **Client identification**: By API key (authenticated) or IP address (unauthenticated/optional mode).
- **Batch rate limiting**: Each item in a batch request consumes one rate-limit unit. A batch of 20 items uses 20 units.
- **Headers**: Rate-limit status is enforced via HTTP 429 responses with ``Retry-After`` guidance.
- **Configuration**: Set ``BIOCOMPILER_RATE_LIMIT`` to change the default.

Core Endpoints
--------------

Type-Check a Sequence
^^^^^^^^^^^^^^^^^^^^^^

.. http:post:: /check

   Evaluate all registered predicates against a DNA sequence and return a
   type-check verdict with an optional certificate.

   **Request Body** — :class:`SequenceInput`:

   .. list-table::
      :header-rows: 1
      :widths: 20 15 65

      * - Field
        - Type
        - Description
      * - ``sequence``
        - ``str``
        - DNA sequence (ACGTN). Required.
      * - ``organism``
        - ``str``
        - Target organism. Default: ``Homo_sapiens``. Accepts canonical names, short keys, abbreviated binomials, or display names.
      * - ``species``
        - ``str | None``
        - Deprecated alias for ``organism``. If both provided, ``species`` takes precedence.
      * - ``exon_boundaries``
        - ``list[tuple[int, int]] | None``
        - Exon boundaries as ``[(start, end), ...]``.
      * - ``gc_lo``
        - ``float``
        - Minimum GC content. Default: 0.30.
      * - ``gc_hi``
        - ``float``
        - Maximum GC content. Default: 0.70.
      * - ``cai_threshold``
        - ``float``
        - Minimum CAI threshold. Default: 0.5.
      * - ``enzymes``
        - ``list[str] | None``
        - Restriction enzymes to check for.
      * - ``cellular_context``
        - ``str``
        - Cellular context for splicing. Default: ``HEK293T``.

   **Response** — :class:`TypeCheckResponse`:

   .. list-table::
      :header-rows: 1
      :widths: 25 75

      * - Field
        - Description
      * - ``sequence_length``
        - Length of the input sequence in base pairs.
      * - ``gc_content``
        - Global GC content as a fraction (0.0–1.0).
      * - ``protein``
        - Translated protein sequence (truncated to 100 aa in summary).
      * - ``results``
        - Per-predicate evaluation results (predicate name, verdict, violation, knowledge_gap).
      * - ``overall_verdict``
        - ``PASS``, ``FAIL``, or ``UNCERTAIN``.
      * - ``certificate``
        - Certificate dict if overall verdict is PASS.

   :statuscode 200: Sequence checked successfully.
   :statuscode 401: Missing or invalid API key.
   :statuscode 422: Validation error in input.

Optimize a Sequence
^^^^^^^^^^^^^^^^^^^

.. http:post:: /optimize

   Optimize a protein sequence for a target organism, producing a codon-optimized
   DNA sequence that satisfies all constraint predicates.

   **Request Body** — :class:`ProteinInput`:

   .. list-table::
      :header-rows: 1
      :widths: 25 15 60

      * - Field
        - Type
        - Description
      * - ``protein``
        - ``str``
        - Target protein sequence (single-letter amino acid codes). Required.
      * - ``organism``
        - ``str``
        - Target organism. Default: ``Homo_sapiens``.
      * - ``species``
        - ``str | None``
        - Deprecated alias for ``organism``.
      * - ``gc_lo`` / ``gc_hi``
        - ``float``
        - GC content bounds. Default: 0.30 / 0.70.
      * - ``cai_threshold``
        - ``float``
        - Minimum CAI threshold. Default: 0.2.
      * - ``enzymes``
        - ``list[str] | None``
        - Restriction enzymes to avoid.
      * - ``cryptic_splice_threshold``
        - ``float``
        - MaxEntScan threshold for cryptic splice site detection. Default: 3.0.
      * - ``organism_domain``
        - ``str``
        - ``auto``, ``eukaryote``, or ``prokaryote``. Default: ``auto``.
      * - ``source_organism``
        - ``str | None``
        - Organism the protein originates from (for immunogenicity predicates). Default: same as ``organism`` (self-protein).
      * - ``therapeutic``
        - ``bool``
        - Apply stricter immunogenicity thresholds. Default: ``false``.
      * - ``self_protein``
        - ``bool | None``
        - Override self-protein detection. Default: auto-detect from ``source_organism``.
      * - ``strict_mode``
        - ``bool``
        - If true, refuse to return sequences with failed predicates (HTTP 422). Default: ``true``.
      * - ``track_provenance``
        - ``bool``
        - Track provenance for this optimization. Default: ``true``.

   **Response** — :class:`OptimizeResponse`:

   .. list-table::
      :header-rows: 1
      :widths: 25 75

      * - Field
        - Description
      * - ``sequence``
        - Optimized DNA sequence.
      * - ``protein``
        - Target protein sequence.
      * - ``cai``
        - Codon Adaptation Index of the optimized sequence.
      * - ``gc_content``
        - GC content fraction.
      * - ``satisfied_predicates``
        - List of predicate names that passed.
      * - ``failed_predicates``
        - List of predicate names that failed.
      * - ``fallback_used``
        - Whether a fallback optimizer was used.
      * - ``provenance_id``
        - Provenance trail ID if tracking was enabled.
      * - ``organism_domain``
        - Resolved organism domain (``eukaryote`` or ``prokaryote``).
      * - ``source_organism``
        - Resolved source organism.
      * - ``therapeutic``
        - Whether therapeutic-mode thresholds were applied.
      * - ``self_protein``
        - Self-protein status used for immunogenicity.

   :statuscode 200: Optimization succeeded.
   :statuscode 422: Optimization failed all predicates (strict mode).

Verify a Certificate
^^^^^^^^^^^^^^^^^^^^

.. http:post:: /verify

   Verify a previously generated certificate.

   **Request Body** — :class:`CertificateInput`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Type
        - Description
      * - ``certificate``
        - ``dict``
        - Certificate as a JSON dict.

   **Response** — :class:`VerifyResponse`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Description
      * - ``status``
        - ``valid`` or ``invalid``.
      * - ``failure_reasons``
        - List of reasons the certificate is invalid (empty if valid).

Scan a Sequence
^^^^^^^^^^^^^^^^

.. http:post:: /scan

   Scan a DNA sequence for motifs, restriction sites, and open reading frames.

   **Request Body** — :class:`ScanInput`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Type
        - Description
      * - ``sequence``
        - ``str``
        - DNA sequence. Required.
      * - ``enzymes``
        - ``list[str] | None``
        - Restriction enzymes to scan for.
      * - ``find_orfs``
        - ``bool``
        - Whether to find open reading frames. Default: ``false``.

   **Response** — :class:`ScanResponse`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Description
      * - ``sequence_length``
        - Length in base pairs.
      * - ``tokens``
        - List of detected motifs/sites.
      * - ``orfs``
        - List of detected ORFs (if requested).

Export Endpoints
----------------

Export to FASTA
^^^^^^^^^^^^^^^

.. http:post:: /export/fasta

   Export a DNA sequence in FASTA format.

   **Request Body** — :class:`ExportFastaInput`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Type
        - Default
        - Description
      * - ``sequence``
        - ``str``
        - *(required)*
        - DNA sequence.
      * - ``identifier``
        - ``str``
        - ``BioCompiler_design``
        - Sequence identifier in the FASTA header.
      * - ``description``
        - ``str``
        - ``""``
        - Description line.
      * - ``organism``
        - ``str``
        - ``Homo_sapiens``
        - Source organism.

   **Response** — :class:`ExportFastaResponse`:

   - ``format``: ``"fasta"``
   - ``content``: FASTA-formatted string.

Export to GenBank
^^^^^^^^^^^^^^^^^

.. http:post:: /export/genbank

   Export a DNA sequence in GenBank flat-file format.

   **Request Body** — :class:`ExportGenbankInput`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Type
        - Default
        - Description
      * - ``sequence``
        - ``str``
        - *(required)*
        - DNA sequence.
      * - ``locus_name``
        - ``str``
        - ``BIOCOMPILER``
        - LOCUS name (max 16 characters).
      * - ``definition``
        - ``str``
        - ``"BioCompiler designed sequence"``
        - DEFINITION line.
      * - ``organism``
        - ``str``
        - ``Homo_sapiens``
        - Source organism.
      * - ``gene_name``
        - ``str | None``
        - ``None``
        - Gene name annotation.
      * - ``exon_boundaries``
        - ``list | None``
        - ``None``
        - Exon boundary annotations.
      * - ``certificate``
        - ``dict | None``
        - ``None``
        - Certificate dict to embed in the GenBank file.

   **Response** — :class:`ExportGenbankResponse`:

   - ``format``: ``"genbank"``
   - ``content``: GenBank-formatted string.

Export to SBOL3
^^^^^^^^^^^^^^^

.. http:post:: /export/sbol3

   Export a DNA sequence as an SBOL3 (Synthetic Biology Open Language) document.

   **Request Body** — :class:`ExportSbol3Input`:

   .. list-table::
      :header-rows: 1

      * - Field
        - Type
        - Default
        - Description
      * - ``sequence``
        - ``str``
        - *(required)*
        - DNA sequence.
      * - ``organism``
        - ``str``
        - ``Homo_sapiens``
        - Target organism.
      * - ``gene_name``
        - ``str``
        - ``optimized_gene``
        - Gene name for ``display_id``.
      * - ``base_uri``
        - ``str``
        - ``https://biocompiler.org/sbol3``
        - Base URI for SBOL identities.
      * - ``format``
        - ``str``
        - ``sbol3``
        - ``sbol3`` (XML) or ``sbol3json`` (JSON-LD).

   **Response** — :class:`ExportSbol3Response`:

   - ``format``: Format used (``sbol3`` or ``sbol3json``).
   - ``content``: SBOL3 document content (XML or JSON-LD).

   See :doc:`sbol3_export` for detailed SBOL3 documentation.

Information Endpoints
---------------------

Health Check
^^^^^^^^^^^^

.. http:get:: /health

   Returns server health status, version, and configuration.

   **Response** — :class:`HealthResponse`:

   - ``status``: ``"healthy"``
   - ``version``: BioCompiler version string.
   - ``timestamp``: ISO 8601 timestamp.
   - ``auth_enabled``: Whether authentication is currently enabled.
   - ``rate_limit_rpm``: Current rate limit setting.

Server Info
^^^^^^^^^^^

.. http:get:: /info

   Returns server configuration limits and supported organisms.

   **Response** — :class:`InfoResponse`:

   - ``max_protein_length``, ``max_batch_size``, ``max_request_size``, ``optimize_timeout_s``
   - ``supported_organisms``: List of canonical organism names.
   - ``api_version``, ``safety_version``

List Organisms
^^^^^^^^^^^^^^

.. http:get:: /organisms

   Returns all supported organisms with their codon usage and GC target data.

   **Response** — :class:`OrganismResponse`:

   - ``organisms``: List of organism detail dicts.

List Predicates
^^^^^^^^^^^^^^^

.. http:get:: /predicates

   Returns all registered type-system predicates.

   **Response** — :class:`PredicateResponse`:

   - ``predicates``: List of predicate name strings.

List Enzymes
^^^^^^^^^^^^

.. http:get:: /enzymes

   Returns all supported restriction enzymes and their recognition sites.

   **Response** — :class:`EnzymeListResponse`:

   - ``enzymes``: Dict mapping enzyme name to recognition site sequence.

Batch Endpoints
---------------

BioCompiler supports batch operations for high-throughput workflows. Each item
in a batch is processed independently — one failure does not affect others.

Batch Type-Check
^^^^^^^^^^^^^^^^

.. http:post:: /batch/check

   Type-check up to **50** DNA sequences in a single request.

   **Request Body** — :class:`BatchCheckInput`:

   - ``sequences``: List of :class:`BatchCheckItem` objects (max 50).

   **Response** — :class:`BatchCheckResponse`:

   - ``results``: Per-item :class:`TypeCheckResponse` dicts.
   - ``summary``: :class:`BatchCheckSummary` with counts (``total``, ``pass``, ``fail``, ``uncertain``, ``errors``).

Batch Optimize
^^^^^^^^^^^^^^

.. http:post:: /batch/optimize

   Optimize up to **20** proteins in a single request.

   **Request Body** — :class:`BatchOptimizeInput`:

   - ``proteins``: List of :class:`BatchOptimizeItem` objects (max 20). Each item can have its own organism and parameters.

   **Fast Batch** — :class:`FastBatchOptimizeInput`:

   If all proteins share the same organism and parameters, use the simplified
   fast batch input with a shared ``organism``, ``gc_lo``, ``gc_hi``, etc.
   This reuses a single optimizer instance and is significantly faster.

   **Response** — :class:`BatchOptimizeResponse`:

   - ``results``: Per-item :class:`OptimizeResponse` dicts.
   - ``summary``: :class:`BatchOptimizeSummary` with counts (``total``, ``all_satisfied``, ``partial``, ``errors``).

Batch Export
^^^^^^^^^^^^

.. http:post:: /batch/export

   Export up to **50** sequences in a single request.

   **Request Body** — :class:`BatchExportInput`:

   - ``sequences``: List of :class:`BatchExportItem` objects (max 50). Each item can specify ``fasta`` or ``genbank`` format independently.

   **Response** — :class:`BatchExportResponse`:

   - ``results``: List of :class:`BatchExportResultItem` objects.

Protein Analysis Endpoints
--------------------------

All protein analysis endpoints are mounted under ``/protein/``.

Structure Prediction
^^^^^^^^^^^^^^^^^^^^

.. http:post:: /protein/structure/predict

   Predict protein structure using ESMFold.

   **Input**: :class:`StructurePredictInput` (``protein``, ``organism``, ``use_cache``).

   **Response**: :class:`StructurePredictResponse` (``pdb_string``, ``mean_plddt``, ``plddt_scores``, ``quality_class``, ``execution_time_s``).

.. http:post:: /protein/structure/batch

   Batch structure prediction (max 20).

.. http:post:: /protein/structure/quality

   Assess structure quality from a PDB string.

   **Input**: :class:`QualityAssessInput` (``pdb_string``).

   **Response**: :class:`QualityAssessResponse` (``mean_plddt``, ``ramachandran_favored``, ``clash_score``, ``overall_quality``, ``verdict``).

Stability Analysis
^^^^^^^^^^^^^^^^^^

.. http:post:: /protein/stability/analyze

   Analyze protein stability.

   **Input**: :class:`StabilityInput` (``protein``, ``organism``, ``pdb_string``).

.. http:post:: /protein/stability/mutations

   Scan mutations for stability effects.

   **Input**: :class:`MutationScanInput` (``protein``, ``positions``, ``method``).

.. http:post:: /protein/stability/batch

   Batch stability analysis (max 20).

Solubility Analysis
^^^^^^^^^^^^^^^^^^^

.. http:post:: /protein/solubility/analyze

   Analyze protein solubility.

   **Input**: :class:`SolubilityInput` (``protein``, ``pdb_string``).

.. http:post:: /protein/solubility/mutations

   Find solubility-improving mutations.

.. http:post:: /protein/solubility/batch

   Batch solubility analysis (max 20).

Immunogenicity Analysis
^^^^^^^^^^^^^^^^^^^^^^^

.. http:post:: /protein/immunogenicity/analyze

   Analyze protein immunogenicity.

   **Input**: :class:`ImmunogenicityInput` (``protein``, ``organism``, ``mhc_alleles``, ``source_organism``, ``therapeutic``, ``self_protein``).

.. http:post:: /protein/immunogenicity/deimmunize

   Deimmunize a protein by introducing mutations that reduce MHC binding.

   **Input**: :class:`DeimmunizeInput` (``protein``, ``organism``, ``target_score``, ``max_mutations``, ``blosum62_min``, ``source_organism``, ``therapeutic``).

.. http:post:: /protein/immunogenicity/batch

   Batch immunogenicity analysis (max 20).

Full Assessment
^^^^^^^^^^^^^^^

.. http:post:: /protein/assessment/full

   Run a comprehensive protein assessment combining structure, stability,
   solubility, and immunogenicity analyses.

   **Input**: :class:`FullAssessmentInput` (``protein``, ``organism``, ``pdb_string``, ``run_structure``, ``run_stability``, ``run_solubility``, ``run_immunogenicity``).

Provenance Endpoints
--------------------

List Provenance Records
^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /provenance

   Query and list provenance records.

   **Response** — :class:`ProvenanceListResponse`:

   - ``count``: Total number of matching records.
   - ``records``: List of :class:`ProvenanceRecordSummary` objects.

Get Provenance Detail
^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /provenance/{record_id}

   Retrieve a specific provenance record by ID.

   **Response** — :class:`ProvenanceDetailResponse`:

   - ``id``: Record ID.
   - ``trail``: Full decision-level provenance trail.

Organism Domain Resolution
--------------------------

.. autofunction:: resolve_organism_domain

Input Validation Helpers
------------------------

.. autofunction:: validate_protein_input

.. autofunction:: validate_organism_input

Pydantic Request Models
-----------------------

Input Models
^^^^^^^^^^^^

.. autoclass:: SequenceInput
   :members:
   :undoc-members:

.. autoclass:: ProteinInput
   :members:
   :undoc-members:

.. autoclass:: CertificateInput
   :members:

.. autoclass:: ExportFastaInput
   :members:

.. autoclass:: ExportGenbankInput
   :members:

.. autoclass:: ExportSbol3Input
   :members:

.. autoclass:: ScanInput
   :members:

Batch Input Models
^^^^^^^^^^^^^^^^^^

.. autoclass:: BatchCheckItem
   :members:

.. autoclass:: BatchCheckInput
   :members:

.. autoclass:: BatchOptimizeItem
   :members:

.. autoclass:: BatchOptimizeInput
   :members:

.. autoclass:: FastBatchOptimizeInput
   :members:

.. autoclass:: BatchExportItem
   :members:

.. autoclass:: BatchExportInput
   :members:

Protein Analysis Input Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: StructurePredictInput
   :members:

.. autoclass:: QualityAssessInput
   :members:

.. autoclass:: StabilityInput
   :members:

.. autoclass:: MutationScanInput
   :members:

.. autoclass:: SolubilityInput
   :members:

.. autoclass:: ImmunogenicityInput
   :members:

.. autoclass:: DeimmunizeInput
   :members:

.. autoclass:: FullAssessmentInput
   :members:

Pydantic Response Models
------------------------

.. autoclass:: TypeCheckResponse
   :members:

.. autoclass:: OptimizeResponse
   :members:

.. autoclass:: VerifyResponse
   :members:

.. autoclass:: ScanResponse
   :members:

.. autoclass:: OrganismResponse
   :members:

.. autoclass:: PredicateResponse
   :members:

.. autoclass:: HealthResponse
   :members:

.. autoclass:: InfoResponse
   :members:

.. autoclass:: EnzymeListResponse
   :members:

.. autoclass:: ExportFastaResponse
   :members:

.. autoclass:: ExportGenbankResponse
   :members:

.. autoclass:: ExportSbol3Response
   :members:

Batch Response Models
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: BatchCheckSummary
   :members:

.. autoclass:: BatchCheckResponse
   :members:

.. autoclass:: BatchOptimizeSummary
   :members:

.. autoclass:: BatchOptimizeResponse
   :members:

.. autoclass:: BatchExportResultItem
   :members:

.. autoclass:: BatchExportResponse
   :members:

Protein Analysis Response Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: StructurePredictResponse
   :members:

.. autoclass:: QualityAssessResponse
   :members:

.. autoclass:: StabilityResponse
   :members:

.. autoclass:: MutationScanResponse
   :members:

.. autoclass:: SolubilityResponse
   :members:

.. autoclass:: ImmunogenicityResponse
   :members:

.. autoclass:: DeimmunizeResponse
   :members:

.. autoclass:: FullAssessmentResponse
   :members:

Provenance Response Models
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ProvenanceResponse
   :members:

.. autoclass:: ProvenanceExplainResponse
   :members:

.. autoclass:: ProvenanceReportResponse
   :members:

.. autoclass:: ProvenanceRecordSummary
   :members:

.. autoclass:: ProvenanceListResponse
   :members:

.. autoclass:: ProvenanceDetailResponse
   :members:

Authentication Functions
------------------------

.. autofunction:: verify_api_key

.. autofunction:: set_no_auth_flag

.. autofunction:: get_auth_mode

.. autofunction:: get_configured_api_keys

.. autofunction:: is_auth_enabled

Constants
---------

.. data:: API_KEY_NAME

   The HTTP header name for API key authentication (``X-API-Key``).

.. data:: RATE_LIMIT_RPM

   Current rate limit in requests per minute (default: 60).

.. data:: BATCH_CHECK_MAX

   Maximum items per batch type-check request (50).

.. data:: BATCH_OPTIMIZE_MAX

   Maximum items per batch optimize request (20).

.. data:: BATCH_EXPORT_MAX

   Maximum items per batch export request (50).

.. data:: BATCH_ITEM_TIMEOUT_S

   Per-item timeout for batch requests in seconds (default: 30).

.. data:: ESMFOLD_TIMEOUT_S

   Timeout for ESMFold structure prediction in seconds (default: 120).

.. data:: MAX_PROTEIN_LENGTH

   Maximum protein sequence length in amino acids (default: 10000).

.. data:: MAX_BATCH_SIZE

   Maximum items per batch request (default: 50).

.. data:: MAX_REQUEST_SIZE

   Maximum request body size in bytes (default: 10,000,000).

.. data:: OPTIMIZE_TIMEOUT_S

   Optimization timeout in seconds (default: 300).
