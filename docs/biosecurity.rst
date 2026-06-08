Biosecurity Screening
=====================

BioCompiler includes a built-in biosecurity screening module that detects
hazardous biological sequences before optimization, preventing accidental or
intentional design of harmful constructs.


Screening Categories
--------------------

Select Agent Toxins
^^^^^^^^^^^^^^^^^^^

Detects catalytic and binding motifs from regulated toxins including:

- **Ricin** A-chain and B-chain (ribosome-inactivating protein)
- **Abrin** A-chain (Type II RIP)
- **Botulinum** neurotoxin (zinc endopeptidase)
- **Shiga** toxin A and B subunits
- **Diphtheria** toxin (ADP-ribosyltransferase)
- **Tetanus** toxin (zinc endopeptidase)
- **Cholera** toxin (ADP-ribosyltransferase)
- **Anthrax** EF, LF, and PA (edema factor, lethal factor, protective antigen)
- **SEB** (staphylococcal enterotoxin B, superantigen)

Risk level: **critical** for active-site motifs, **high** for binding motifs.

Viral Surface Proteins
^^^^^^^^^^^^^^^^^^^^^^

Detects fusion peptides, receptor-binding domains, and cleavage sites from:

- **Influenza** HA (fusion peptide, receptor binding, cleavage site) and NA (active site)
- **SARS-CoV-2** spike (RBD, fusion peptide, furin cleavage, heptad repeat)
- **HIV-1** Env (V3 loop, gp41 fusion core, CD4 binding)
- **Ebola** GP (receptor binding, fusion peptide, mucin-like domain)
- **Variola** (smallpox) envelope protein

Risk level: **high** for fusion/cleavage motifs, **medium** for receptor-binding motifs.

Antibiotic Resistance Markers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Detects both protein motifs and nucleotide patterns for:

- **blaTEM** (TEM beta-lactamase)
- **nptII** (kanamycin resistance)
- **aac(6')** (aminoglycoside acetyltransferase)
- **cat** (chloramphenicol acetyltransferase)
- **tetA** and **tetM/tetO** (tetracycline resistance)
- **vanA** (vancomycin resistance)
- **mecA** (methicillin resistance, PBP2a)
- **ctx-m** (extended-spectrum beta-lactamase)
- **ndm-1** (New Delhi metallo-beta-lactamase)

Risk level: **high** for ESBL/carbapenemase markers, **medium** for others.

Oncogenes and Growth Factors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Detects oncogenic protein motifs:

- **MYC** (transcriptional activation, bHLH domain)
- **RAS** (GTP binding, switch I/II regions)
- **EGFR** and **HER2** (kinase domain)
- **BRAF** (activation segment, V600E region)
- **VEGF**, **PDGF**, **TGF-beta** (receptor binding)
- **p53** (DNA binding domain mutations)

Risk level: **medium** for kinase/activation motifs, **low** for receptor-binding.


Screening API
-------------

.. code-block:: python

   from biocompiler.biosecurity import screen_hazardous_sequence, BiosecurityReport

   # Screen a protein sequence
   report = screen_hazardous_sequence(protein="NIRVGLPIIS...", dna="ATGGCG...")

   # Check results
   print(f"Hazardous: {report.is_hazardous}")
   print(f"Risk level: {report.risk_level}")
   print(f"Categories: {report.flagged_categories}")

   for match in report.matches:
       print(f"  {match.category}/{match.name} at pos {match.position} "
             f"(confidence={match.confidence:.2f}, type={match.match_type})")


Fuzzy Matching
--------------

The biosecurity module supports fuzzy matching to catch near-matches:

- **Hamming distance**: For short peptide motifs (<15 aa), allows 1–2 substitutions
- **Levenshtein distance**: For motifs with insertions/deletions, allows 1 edit
- **Reverse complement**: DNA screening also checks the reverse strand

Fuzzy matches receive reduced confidence scores based on edit distance.


Biosecurity Modes
-----------------

The screening behavior is controlled by the ``BIOCOMPILER_BIOSECURITY_MODE``
environment variable:

=========  ================================================================
Mode       Behavior
=========  ================================================================
enforce    (Default) Block optimization if hazardous sequences detected
warn       Log warnings but allow optimization to proceed
off        Disable biosecurity screening entirely
=========  ================================================================

.. code-block:: bash

   # Enable warning mode
   export BIOCOMPILER_BIOSECURITY_MODE=warn
   biocompiler optimize MVLSPADKTNVKAAWGKVGA --organism human

   # Disable screening (not recommended for production)
   export BIOCOMPILER_BIOSECURITY_MODE=off


Risk Levels
-----------

=========  ================================================================
Level      Meaning
=========  ================================================================
none       No hazardous sequences detected
low        Minor concern, unlikely to be functional
medium     Moderate concern, may be functional
high       Significant concern, likely functional
critical   Severe concern, select agent or equivalent
=========  ================================================================


Regulatory References
---------------------

- CDC Select Agent Program (42 CFR Part 73)
- Australia Group Common Control List
- WHO Laboratory Biosafety Manual, 4th ed. (2020)
- CARD: Comprehensive Antibiotic Resistance Database (https://card.mcmaster.ca)
