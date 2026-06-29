# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.9.3 (2026-06-29) — Initial public release

### Core system

- **Five-level IR compiler**: Genomic DNA → pre-mRNA → mature mRNA → polypeptide → folded protein, with formally verified lowering passes between each level (`transcribe`, `splice`, `translate`, `fold`).
- **Lean 4 formalization** (`proof/BioCompiler/`, 17 modules, 8,746 LOC, **267 theorems**):
  - 17 core deterministic predicates fully proved (**0 `sorry`** in the whole development).
  - 19 SLOT predicates conservatively bounded (return `UNCERTAIN` in conservative mode; cannot produce a false `PASS`).
  - **0 `sorry` remain** in `SLOTVerification.lean` (the 2 former BLOSUM62-related `sorry` were discharged in W1-A5 by formalizing the BLOSUM62 substitution matrix in `proof/BioCompiler/BLOSUM62.lean`, 11 theorems).
  - The former 15 broad `axiom` declarations have been **narrowed to 34 specific, independently-testable `axiom` declarations** in `SLOTVerification.lean`, each a single-property tool-soundness contract for an external ML model (TMHMM, ViennaRNA, ESMFold/AlphaFold, FoldX, CamSol, Aggrescan, ExPASy, NetMHCpan/MHCflurry, BePiPred, IEDB). Each narrowed axiom is backed by a runtime evidence check in `src/biocompiler/provenance/runtime_evidence.py`.
  - TCB: 3 class-field axioms in `SpliceSiteScanner` (scanner completeness, soundness, borderline completeness).
  - Central theorem: `compile_correctness` — the full L0→L3 pipeline produces the correct protein.
- **Integrated constraint-solving optimizer**: greedy forward pass with bounded multi-pass cleanup. Achieves CAI 1.000 on E. coli targets (vs DNAchisel's 0.948), 0 constraint violations. With the opt-in **context-aware GT-avoidance mode** (`use_context_aware_gt=True`), achieves CAI 0.963 on human HBB (vs. 0.718 for the default global-GT-avoidance path on the same sequence) by repairing only high-scoring cryptic splice donors rather than forbidding every GT dinucleotide.
- **Biosecurity screening**: 89 hazard signatures across 5 categories (select agent toxins, oncogenes, viral surface proteins, antibiotic resistance markers, Australia Group pathogens). Pre-optimization exact matching (<1ms) and post-optimization fuzzy matching (Hamming + Levenshtein with k-mer pre-filtering).
- **SECIS-aware selenocysteine translation**: UGA codons at SECIS positions emit U (selenocysteine) rather than stop.
- **Standalone certificate verifier** (`scripts/standalone_verifier.py`, 473 LOC, stdlib-only) for independent re-checking of optimization certificates.
- **Codegen targets**: GenBank, FASTA (protein/mRNA), SBOL3.

### Test suite

- **14,129 tests collected** across 270 test files; **652 deselected by default** (slow/benchmark/external-tool markers). **389 tests marked `requires_external`** (require external tools like NetMHCpan, ViennaRNA, ESMFold that are not installed in CI).
- **IR + certificate + provenance stack**: 308/308 pass.
- **25-selenoprotein end-to-end test suite** (`tests/test_selenoproteins_e2e.py`): 130 tests covering all 25 known human selenoproteins. The catalog in `src/biocompiler/selenoproteins.py` records verified Sec positions from Kryukov et al. 2003, Gladyshev et al. 2016, and UniProtKB feature tables. Includes multi-Sec proteins SELENOP (10 Sec residues) and SELENOU (24 Sec residues).
- **10 SECIS unit tests** (`tests/test_secis.py`): synthetic test cases covering UGA recoding, multiple SECIS positions, and end-to-end back-translation.
- **500 `@given` decorators** across 18 test files using Hypothesis for property-based testing.
- **Heavy fair benchmark** (`heavy_benchmark_results.json`): 25-gene head-to-head comparison vs DNAchisel. BC wins E. coli (CAI 1.000 vs 0.948, 2-4× faster); DNAchisel wins human on the default path (CAI 0.992 vs 0.792, 7× faster). Both achieve 0 constraint violations. With the opt-in context-aware GT-avoidance mode, BC achieves CAI 0.963 on human HBB (vs. 0.718 default-path on the same sequence). See `BENCHMARKS.md` for full methodology.
- **Runtime evidence checks** (`tests/test_runtime_evidence.py`): 76 tests (34 positive + 34 negative + 8 module/aggregator) covering the 34 narrowed tool-soundness axioms in `SLOTVerification.lean`. Each check verifies that the corresponding tool's output is self-consistent and within expected ranges at runtime; this catches tool malfunctions even though it does not prove the axioms.

### Verification

- **Lean 4 proofs build cleanly** with `lake build` using Lean 4.30.0 (elan-managed). All 19 build jobs succeed; only linter warnings (unused simp arguments, unused variables), no errors.
- **Proof statistics verified against docs**:
  - 267 theorems (`grep -c "^theorem " proof/BioCompiler/*.lean`)
  - 34 narrowed `axiom` declarations in SLOTVerification.lean (narrowed from 15 former broad axioms); plus 0 in BLOSUM62.lean (BLOSUM62 case fully discharged, no axiom needed) (`grep -c "^axiom " proof/BioCompiler/*.lean`)
  - 0 `sorry` in SLOTVerification.lean (the 2 former BLOSUM62-related sorry were discharged in W1-A5 via BLOSUM62.lean)

### Limitations (stated honestly)

- **Not validated in wet-lab or clinical settings.** BioCompiler outputs are design candidates to be reviewed by a qualified molecular biologist before synthesis.
- **Formal verification is partial.** 17 core predicates are fully proved; 19 SLOT predicates are conservatively bounded but not fully proved (require formalizing external ML models, which is out of scope).
- **Benchmark performance is mixed.** BioCompiler wins on E. coli (prokaryotic) targets but loses to DNAchisel on human (eukaryotic) targets. The value proposition is certified-by-default predicate evaluation and biosecurity screening, not raw speed/CAI.
