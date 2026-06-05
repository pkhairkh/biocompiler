# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
