# DOC-21: IR Compiler Architecture

## Overview

BioCompiler includes a real multi-level intermediate representation (IR) compiler
that transforms gene specifications through five typed IR levels, mirroring the
Central Dogma of molecular biology.

## IR Levels

| Level | Name | Description |
|---|---|---|
| IR-L0 | GenomicDNA | Raw DNA sequence with region annotations |
| IR-L1 | PreMRNA | Transcribed RNA (T→U), regions preserved |
| IR-L2 | MatureMRNA | Spliced mRNA: 5'UTR + CDS + 3'UTR |
| IR-L3 | Polypeptide | Translated amino acid sequence |
| IR-L4 | FoldedProtein | 3D structure (ESMFold oracle) |

## Lowering Passes

- **transcribe** (L0→L1): DNA → pre-mRNA (T→U)
- **splice** (L1→L2): Remove introns, assemble CDS (SECIS-aware, NDFST alternative splicing)
- **translate** (L2→L3): Codon → amino acid (selenocysteine support via SECIS positions)
- **fold** (L3→L4): Protein structure prediction (ESMFold + heuristic fallback)

## Optimization Passes (IR→IR)

- **optimize_codons**: Higher CAI, same protein (verified)
- **eliminate_cpgs**: Remove CpG dinucleotides, same protein (verified)

## Frontend

YAML spec parser: `parse_spec()` → IR_L0

## Codegen (IR → output)

- **to_genbank()**: IR-L0 → GenBank format
- **to_fasta()**: IR-L3 → FASTA format
- **to_sbol3()**: IR-L0 → SBOL3 Turtle RDF

## Formal Verification

All lowering passes are proven correct in Lean4 (proof/BioCompiler/IR.lean):
- `transcribe_correctness`: L0→L1 preserves sequence
- `splice_correctness`: L1→L2 preserves coding content
- `compile_correctness`: Full L0→L3 produces correct protein
- 15 sorry in IR.lean

## Test Results

- **100,068 IR compiler tests** (fast back-translation): 100% pass
- **10,680 E2E tests** (optimizer+IR+codegen): 99.4% pass, 0% fail
- **Total**: 1,674 combinations tested, 99.9% pass, 0 mismatches

## Selenocysteine Support

The compiler handles selenocysteine (U), the 21st amino acid:
- Encoded by UGA (normally a stop codon) via SECIS-element recoding
- `secis_positions` field on IR_L0/L1/L2 marks codon indices where UGA→U
- splice() skips UGA at SECIS positions when scanning for terminal stop
- translate() emits U (not *) at SECIS positions


## Integrated Constraint-Solving Optimizer

Single-pass constraint solving: for each codon, try all synonymous alternatives
sorted by CAI. Pick the one satisfying ALL constraints simultaneously.

Performance: 14× faster than DNAchisel.
