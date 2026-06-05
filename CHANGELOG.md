# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 12.0.0 (2026-03-05)

### BREAKING CHANGES
- **Version bumped to 12.0.0**: Major version reflects the cumulative scope of new features and safety infrastructure added since v11.1.0.

### New Features
- **Biosecurity screening pipeline**: Full `biosecurity` module with hazard signature database covering select agent toxins (ricin, abrin, botulinum, shiga, diphtheria, tetanus, cholera, anthrax EF/LF), viral surface proteins (influenza HA/NA, SARS-CoV-2 spike, HIV env, Ebola GP), antibiotic resistance markers (blaTEM, nptII, aac(6'), cat, tetA/M/O, vanA, mecA, ctx-m, ndm-1), and oncogenes/growth factors (MYC, RAS, EGFR, VEGF, BRAF, p53, HER2, PDGF, TGF-beta). Risk-level classification (none/low/medium/high/critical). Biosecurity gate blocks optimization for high/critical risk sequences.
- **Translation verification module** (`protein_verification`): Post-optimization safety check verifying that optimized DNA encodes the expected protein. Position-level mismatch reporting, premature stop detection, terminal stop codon handling. `verify_and_raise()` integration hook for the optimizer.
- **Optimization strict mode**: `OptimizationConstraintError` raised when strict mode (default) encounters unsatisfied predicates. Provides `failed_predicates` list and `partial_result` for inspection.
- **5-organism support**: Full optimization and CAI tables for E. coli, Human, Yeast, Mouse, and CHO. Organism-aware constraint profiles: prokaryotes skip splice/CpG constraints; eukaryotes get full constraint set.
- **Organism configuration system** (`organism_config`): `OrganismConfig` dataclass with GC targets, preferred codons, avoided motifs, max homopolymer runs, mRNA degradation models, and domain classification. Constraint profiles per organism.
- **Provenance and decision audit trail**: `DecisionRecord`, `ProvenanceTracker`, `OptimizationProvenance`, `OptimizationRecord` for full decision-level tracking during optimization. `decision_provenance` module with `CodonDecision`, `ConstraintDecision`, `DecisionProvenanceCollector`, `OptimizationDecisionTrail`.
- **UTR models**: `utr_models` module with organism-specific 5' and 3' UTR suggestions and scoring (UTRConfig, suggest_5utr, suggest_3utr, score_5utr, score_3utr).
- **Multi-gene / Operon support**: `multigene` module with `GeneSpec`, `MultiGeneResult`, `OperonConfig`, `optimize_multigene`, `optimize_operon`.
- **CSP/SMT Solver integration**: `solver` sub-package with `CSPSolver`, `SolverConfig`, `SolverResult`, `SolverBackend`, `CAIAwareConstraintResolver`, `solve_with_csp`, `csp_optimize`. OR-Tools CP-SAT and Z3 SMT backends.
- **MHC binding database**: `mhc_binding_db` with precomputed PSSMs for HLA-A*01:01, A*02:01, A*03:01, B*07:02, B*08:01, DRB1*01:01, DRB1*04:01, mouse H2-Kb, H2-Db.
- **ESMFold batch structure prediction**: Batch prediction API, caching, pLDDT classification, contact map estimation.
- **Immunogenicity engine**: T-cell epitope prediction (MHC-I, MHC-II PSSM), B-cell epitope prediction (Kolaskar-Tongaonkar, Parker, Chou-Fasman, Emini, BepiPred-like, conformational), population coverage, deimmunization pipeline.
- **Structure quality assessment**: `structure` sub-package with PDB parser, quality report, Ramachandran assessment, clash score, packing density, SASA approximation, structure predicates.
- **FoldX stability integration**: `foldx` module with stability landscape scanning, conservation scoring, compensatory mutation finding, hotspot identification.
- **CamSol solubility prediction**: Intrinsic and structural solubility, solubility-guided optimization, batch computation, `CamSolResult` (with backward-compatible `SolubilityResult` alias).
- **Protein design engine**: `protein_design` module with thermostable, soluble, low-immunogenicity, and multi-objective design modes.
- **mRNA stability module**: `mrna_stability` with destabilizing motif detection and removal.
- **Numba-accelerated kernels**: `numba_kernels` with fast CAI computation, dinucleotide counting, batch codon swap scoring, incremental CAI updates.
- **Aho-Corasick multi-pattern scanner**: `aho_corasick` for efficient restriction site detection.
- **Sliding window GC content**: `sliding_gc` module for local GC content checking and fixing.
- **MaxEntScan splice scoring**: `maxentscan` and `maxentscan_fast` with donor/acceptor scoring, splice site scanning.
- **Codon pair bias scoring**: `codon_pair_scoring` module.
- **Large sequence chunk-based optimization**: `large_sequence` module for proteins >10kb.
- **Incremental constraint tracking**: `incremental` module with `IncrementalSequenceState`, O(1) GC updates.
- **Sequence annotation enrichment**: `annotation` module with `SequenceAnnotation`, `annotate_sequence`, `annotate_to_genbank`.
- **GenBank round-trip**: `genbank_roundtrip` for import/export round-trip testing.
- **IUPAC ambiguity code support**: `iupac` module.
- **Rate limiter**: `rate_limiter` for API rate limiting.
- **LIMS integration**: `lims` module for laboratory information management.
- **What-if analysis**: `whatif_analysis` for exploring optimization alternatives.
- **Pattern enforcement**: `pattern_enforcement` for sequence pattern constraints.
- **Benchmarking sub-package**: `benchmarking` with structured head-to-head comparison, CAI-validated benchmarks, multi-constraint stress testing, organism-aware benchmarks, Sharp-Li reference benchmark.
- **BioPython deep integration**: Codon usage tables, alignment, phylogenetic distance, ORF detection, BLAST, back-translation.
- **SBOL3 export/import**: `sbol_export` and `sbol_import` for synthetic biology open language interchange.
- **Grammar-driven optimization**: `grammar_loader` with YAML grammar definitions for HBB/HEK293T and EGFP/HEK293T.
- **Formal verification proofs**: Lean4 proofs for soundness, SLOT verification, NDFST, compositional reasoning, mutagenesis, type system, certificates, oracle proofs, three-valued logic, refinement mapping, and scanner correctness.

### Bug Fixes
- All fixes from v10.0.0 and v11.1.0 carried forward
- Biosecurity gate correctly blocks optimization for high/critical risk sequences
- Translation verification catches codon table mismatches after optimization
- Strict mode prevents returning sequences with failed predicates

---

## 10.0.0 (2026-06-05)

### BREAKING CHANGES
- **CAI table unification**: The optimizer now uses the same CODON_ADAPTIVENESS_TABLES for both codon selection and CAI evaluation. Previously, the optimizer used SPECIES tables which disagreed with the evaluation tables, causing CAI to be incorrectly reported. CAI values will now be different (correct) from previous versions.
- **E. coli optimal codons corrected**: Five amino acids had incorrect optimal codons in the per_thousand data (Phe, Ile, Tyr, His, Arg). The adaptiveness table now correctly identifies the most frequent codon per amino acid.
- **Organism name resolution**: New `resolve_organism()` function provides centralized name resolution. Both `species` and `organism` parameters now accepted everywhere.

### New Features
- **HybridOptimizer**: New hybrid optimization strategy combining greedy initialization + priority-based constraint satisfaction + CAI hill climbing. 3-5× faster than the legacy pipeline with higher CAI.
- **CSP Solver fixed**: Both OR-Tools CP-SAT and Z3 SMT solver backends now produce valid sequences with correct constraints.
- **Prokaryote fast path**: E. coli and other prokaryotes skip eukaryotic constraints (splice sites, CpG islands), recovering ~0.3 CAI and 3× speed improvement.
- **DP back-translation**: Optional dynamic programming back-translation that avoids cross-codon restriction sites and GT dinucleotides.
- **Incremental constraint checking**: O(1) GC updates, 2-2000× faster constraint re-checking after codon changes.
- **CAI-aware constraint resolution**: All constraint resolution steps now prefer higher-CAI alternatives, minimizing CAI loss during constraint fixing.
- **CAI recovery pass**: Post-optimization pass that upgrades suboptimal codons to optimal synonyms without violating constraints.

### Performance
- E. coli GFP: CAI 0.67 → 0.999, Time 20ms → 2ms (10× faster)
- Speed gap vs DNAchisel: 38× → 2-3× 
- CAI gap vs DNAchisel: 0.31 → 0.001

### Bug Fixes
- Fixed _back_translate_protein using wrong CAI table (SPECIES vs CODON_ADAPTIVENESS_TABLES)
- Fixed 5 E. coli amino acids with incorrect optimal codons
- Fixed OR-Tools and Z3 CSP solvers returning empty sequences
- Fixed prokaryote optimization applying eukaryotic splice/CpG constraints
- Fixed MaxEntScan deprecated model anti-correlation
- Fixed mRNA stability pass degrading CAI by accepting suboptimal codon swaps

---

## 9.2.0 (2025-12-15)

### New Features
- Unified engine API: BaseEngineResult, MutationResult, BatchResult shared across all 6 analysis engines (ESMFold, FoldX, CamSol, Immunogenicity, Deimmunization, Protein Design)
- 28-predicate type system: 12 DNA + 4 structure + 4 stability + 4 solubility + 4 immunogenicity predicates
- SLOT architecture: 13 core predicates (PASS/FAIL) + 19 SLOT-dependent predicates (always UNCERTAIN); Lean4 proof covers all 28 predicates
- HBB full pass: all 8 optimizer predicates pass simultaneously
- CpG reconciliation, CAI reconciliation, cross-codon coordination

### Bug Fixes
- Fixed CpG reconciliation undoing restriction site removals
- Fixed CAI reconciliation reintroducing cryptic splice sites
- Fixed cross-codon coordination oscillation between GT and CG constraints
