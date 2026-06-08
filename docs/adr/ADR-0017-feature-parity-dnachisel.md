# ADR-0017: Feature Parity with DNAchisel

## Status: Accepted

## Date

2026-03-05

## Context

DNAchisel is an established Python library for codon optimization and DNA sequence design that supports a rich set of constraint types and specification formats. BioCompiler was originally designed with only three core constraints (no restriction sites, no cryptic splice sites, GC content range). As users migrated from DNAchisel, they expected equivalent functionality and were blocked by missing features.

The gap between BioCompiler and DNAchisel fell into six categories:

1. **Constraint breadth**: DNAchisel supports 11+ constraint types (AvoidPattern, EnforcePattern, AvoidChanges, EnforceTranslation, AvoidBlastMatches, UniquifyAllKmers, etc.), while BioCompiler had only 3 (restriction sites, splice sites, GC). Users who needed CpG avoidance, IUPAC ambiguous base handling, sliding-window GC, or pattern enforcement could not use BioCompiler.

2. **IUPAC ambiguous base support**: DNAchisel handles sequences with IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) natively. BioCompiler required all inputs to be pure ACGT, forcing users to pre-process sequences manually with no guidance on resolution strategy.

3. **Custom objectives**: DNAchisel allows users to define custom optimization objectives (e.g., minimize GC deviation from target, maximize codon pair score). BioCompiler only maximized CAI — users who wanted to balance CAI against GC content or codon pair bias had no mechanism to do so.

4. **Sliding-window and local GC**: DNAchisel provides local/sliding-window GC constraints that catch extreme GC regions missed by global GC checks. BioCompiler only had global GC constraints, meaning locally extreme GC regions (which cause polymerase processivity issues, secondary structure, or synthesis difficulty) went undetected.

5. **Part libraries and assembly planning**: DNAchisel's ecosystem (via DNA Chisel Plus and related tools) supports standard biological parts, assembly method planning (Gibson, Golden Gate), and construct design. BioCompiler had no parts registry or assembly planning capability, requiring users to manually design constructs and check assembly compatibility.

6. **Interoperability**: DNAchisel users expected SBOL export/import, GenBank round-trip verification, LIMS integration, and automatic annotation. BioCompiler had limited export (GenBank and FASTA only) and no import or round-trip verification.

**Alternatives Considered:**

1. **Wrapper around DNAchisel** — Use DNAchisel as a backend for constraint types that BioCompiler doesn't implement natively. But: this would make BioCompiler dependent on DNAchisel's release cycle, would not integrate with BioCompiler's type system and certificate framework, and would duplicate the optimization pipeline (DNAchisel is a constraint solver, not a compiler with type-directed verification).

2. **Incremental feature addition with no parity goal** — Add features one at a time as users request them. But: this results in an unpredictable roadmap, inconsistent API design (each feature added ad hoc), and continued user frustration during the gap period.

3. **Systematic feature parity with BioCompiler-native implementations** — Implement each DNAchisel feature natively in BioCompiler, using BioCompiler's type system, certificate framework, and optimization pipeline. This ensures consistent API design, integration with the predicate/certificate system, and independence from DNAchisel's codebase.

## Decision

Adopt systematic feature parity (Alternative 3), implementing each DNAchisel feature natively in BioCompiler with the following module assignments:

1. **IUPAC ambiguous base support** (`biocompiler.sequence.iupac`):
   - `resolve_ambiguous()` with strategies: most_common, cai_optimal, gc_balanced, first
   - `expand_ambiguous()` for enumeration of all concrete sequences
   - `validate_iupac_sequence()` for input validation
   - Integration with the optimizer: IUPAC inputs are resolved before optimization using the `cai_optimal` strategy by default

2. **Pattern enforcement** (`biocompiler.sequence.pattern_enforcement`):
   - `PatternConstraint` dataclass with action (enforce/avoid), scope (dna/protein), strand (both/forward/reverse)
   - `check_pattern()` and `enforce_pattern()` for single constraints
   - `enforce_patterns()` for iterative multi-constraint enforcement
   - IUPAC pattern expansion and Aho-Corasick multi-pattern scanning
   - Supersedes the restriction-site-only avoidance; restriction sites are now a special case of avoid-pattern constraints

3. **Custom objectives** (`biocompiler.optimizer.objectives`):
   - `ObjectiveFunction` Protocol: `(dna, protein, organism) -> float`
   - Built-in objectives: cai, cai_gc_balanced, codon_pair, min_max_gc
   - `resolve_objective()` for name/callable resolution
   - Objective refinement pass after standard optimization
   - `objective_score` field in `OptimizationResult`

4. **Sliding-window GC** (`biocompiler.sequence.sliding_gc`):
   - `check_sliding_gc()` with configurable window size, min/max GC, step
   - `fix_sliding_gc_violations()` using CAI-aware synonymous codon substitution
   - `evaluate_sliding_gc()` as a type-system predicate
   - NUMBA-accelerated kernel with pure-Python fallback

5. **Local GC constraints** (`biocompiler.sequence.local_gc`):
   - `LocalGCConstraint` for region-specific GC bounds
   - `check_local_gc()` and `optimize_local_gc()` for fine-grained regional control

6. **Part libraries** (`biocompiler.optimizer.parts`):
   - `Part` dataclass and `PartLibrary` registry
   - Built-in default parts (promoters, RBS, terminators, linkers)
   - YAML/JSON loading from user-provided files
   - Search by type and organism

7. **Assembly planning** (`biocompiler.optimizer.assembly`):
   - `plan_golden_gate()` for Type IIS restriction enzyme assembly
   - `plan_gibson()` for overlap-based seamless assembly
   - Internal restriction site checking for Golden Gate compatibility

8. **SBOL3 export/import** (`biocompiler.export.sbol_export`, `biocompiler.export.sbol_import`):
   - Pure-Python SBOL3 RDF/XML generator (no external SBOL library dependency)
   - Component, Sequence, Measure, Activity, Plan elements
   - JSON-LD alternative serialization
   - Import and conversion to GeneSpec for re-optimization

9. **GenBank round-trip verification** (`biocompiler.export.genbank_roundtrip`):
   - `verify_genbank_roundtrip()` for export→import→compare verification
   - Sequence integrity and annotation preservation checks

10. **LIMS integration** (`biocompiler.lims`):
    - `LIMSIntegration` abstract base class
    - `BenchlingExporter` and `LabGuruExporter` concrete implementations
    - Convenience functions `export_to_benchling()`, `export_to_labguru()`

11. **Sequence annotation** (`biocompiler.export.annotation`):
    - Automatic feature detection: ORFs, restriction sites, CpG islands, splice sites, repeats, GC/AT-rich regions, RBS
    - `annotate_to_genbank()` for fully-annotated GenBank output

## Consequences

- **Positive**: BioCompiler now covers the full feature set that DNAchisel users expect, removing the primary migration barrier. Users can transition from DNAchisel without losing functionality.
- **Positive**: All new features integrate with BioCompiler's type system, certificate framework, and provenance tracking — unlike a DNAchisel wrapper, which would bypass these guarantees.
- **Positive**: IUPAC support enables direct processing of sequences from databases (GenBank, UniProt) that may contain ambiguous bases, without manual pre-processing.
- **Positive**: Custom objectives allow domain-specific optimization (e.g., balancing CAI with codon pair bias for viral vectors, targeting specific GC for synthesis compatibility).
- **Positive**: Sliding-window GC catches local extremes that global GC misses — this is a common source of synthesis failures in practice.
- **Positive**: Part libraries and assembly planning support the full construct design workflow, not just codon optimization.
- **Negative**: The feature expansion significantly increases the module surface area (18+ modules vs. the original ~6). This increases maintenance burden and the learning curve for new users.
- **Negative**: Some features (pattern enforcement, assembly planning) are less mature than their DNAchisel equivalents. Edge cases and performance may differ.
- **Trade-off**: We chose native implementations over a DNAchisel wrapper to maintain type-system integration and independence, but this means we must independently verify correctness for each feature rather than inheriting DNAchisel's test coverage.

## References

- DNAchisel: https://github.com/Edinburgh-Genome-Foundry/DnaChisel
- ADR-0018: tRNA Adaptation Index (tAI) — additional scoring metric beyond CAI
- ADR-0008: Greedy Multi-Objective Optimizer as Default
- ADR-0012: CpG Avoidance in Greedy Optimizer
- ADR-0015: Biosecurity Sequence Screening
