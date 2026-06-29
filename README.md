# BioCompiler 0.9.3

**Status: Beta — not yet validated in wet-lab or clinical settings.**

A formally-verified gene design compiler with multi-level intermediate
representations. BioCompiler applies compiler-engineering design patterns —
staged transformation, typed IR, composable passes, type-based verification —
to the formalizable stages of gene-to-protein processing (transcription,
splicing, translation, folding) and to a constraint-driven codon optimizer
that feeds the IR.

The compiler is real: a YAML spec goes in one end, and a GenBank / FASTA /
SBOL3 file comes out the other, after passing through five typed IR levels
that mirror the Central Dogma. The verification is real: 17 core
predicates are fully proved in Lean4, 19 SLOT predicates are vacuously sound
in conservative mode, and **0 `sorry`** remain in `SLOTVerification.lean`
(the 2 former BLOSUM62-related `sorry` were discharged by formalizing the
BLOSUM62 substitution matrix in `proof/BioCompiler/BLOSUM62.lean`, 11
theorems). The former 15 broad tool-soundness `axiom` declarations have been
**narrowed to 34 specific, independently-testable contracts** — each asserts
a single property of an external tool's output (window-size matching,
threshold-range validity, proxy correctness, or a soundness contract) and is
backed by a runtime evidence check in
`src/biocompiler/provenance/runtime_evidence.py` (34 checks, one per narrowed
axiom). The Lean4 development totals 267 theorems across 17 modules (8,746
LOC).

> **The 17/19 separation, stated honestly.** The 17 core predicates are the
> **deterministic layer** — every predicate whose verdict can be decided by
> sequence computation alone (no external tools). The 19 SLOT predicates are
> **oracle-dependent** — they require external tools (ESMFold, NetMHCpan,
> ViennaRNA) and are handled by SLOT, which proves they **cannot produce a
> false PASS** without external tool evidence. This separation is the
> contribution: deterministic safety is formally guaranteed; oracle-dependent
> safety is conservatively bounded.

> **The 0 open `sorry`, stated honestly.** 0 `sorry` remain in
> `SLOTVerification.lean` — the 2 former BLOSUM62-related `sorry` were
> discharged by formalizing the BLOSUM62 substitution matrix in
> `proof/BioCompiler/BLOSUM62.lean` (11 theorems). The former 15 broad
> tool-soundness `axiom` declarations have been narrowed to **34 specific,
> independently-testable contracts**, each asserting a single property of an
> external tool's output and each backed by a runtime evidence check in
> `src/biocompiler/provenance/runtime_evidence.py`. These narrowed axioms
> are **tool-interface correctness properties** — they assert that external
> tool outputs faithfully reflect biological properties. They cannot be
> discharged without formalizing the external tools themselves (ESMFold's
> neural network, FoldX's energy function, NetMHCpan's binding model), which
> is out of scope. Each is documented as `NEEDS_TOOL_AXIOM`. The proof
> **degrades gracefully**: in conservative mode, the system never produces a
> false PASS even with these obligations open.

---

## The IR Pipeline

This is the heart of the compiler. A gene specification flows through five
typed IR levels, each produced by a pure, invariant-checked lowering pass:

```
YAML spec -> IR-L0 (GenomicDNA) -> transcribe -> IR-L1 (PreMRNA) -> splice -> IR-L2 (MatureMRNA) -> translate -> IR-L3 (Polypeptide) -> fold -> IR-L4 (FoldedProtein)
                                    |                                            |
                                    v                                            v
                          NDFST alternative splicing                       GenBank / FASTA / SBOL3 codegen
                   (exon-skip / intron-retention / no-splice
                    isoforms, MaxEntScan-scored; SECIS-aware)
```

Each arrow is a pure, invariant-checked lowering pass. The NDFST branch
means a single IR-L1 can yield multiple IR-L2 isoforms. Codegen can emit
from IR-L0 (GenBank/SBOL3), IR-L2 (FASTA mRNA), or IR-L3 (FASTA protein).

### IR levels (`biocompiler.ir.types`)

| Level | Type | Carries |
|-------|------|---------|
| L0 | `IR_L0_GenomicDNA` | DNA sequence, region annotations, SECIS positions |
| L1 | `IR_L1_PreMRNA` | RNA sequence (T->U), regions preserved |
| L2 | `IR_L2_MatureMRNA` | `5'UTR` + `CDS` + `3'UTR` (spliced) |
| L3 | `IR_L3_Polypeptide` | Amino-acid sequence (with `*` stop, `U` Sec) |
| L4 | `IR_L4_FoldedProtein` | 3-D coordinates (oracle), pLDDT, SS string |

### Lowering passes (`biocompiler.ir.passes`)

- **`transcribe`** (L0->L1): DNA -> pre-mRNA. Pure, validated.
- **`splice`** (L1->L2): intron removal, CDS assembly. SECIS-aware (UGA at a
  SECIS position is selenocysteine, not stop). Falls back to start/stop scan
  when no region annotations are present.
- **`translate`** (L2->L3): standard genetic code. UGA recoded to `U`
  (selenocysteine) at SECIS positions; unknown codons -> `X`.
- **`fold`** (L3->L4): ESMFold oracle with Chou-Fasman heuristic fallback.
  Always returns a result; `metadata["oracle"]` records which backend ran.

### Optimization passes (`biocompiler.ir.optimization`)

IR->IR passes that **preserve the translated protein** (verified by
re-translation after each pass - `IRError` is raised on any drift):

- **`optimize_codons`** - codon optimization for the target organism (CAI).
- **`eliminate_cpgs`** - remove CpG dinucleotides from the CDS (relevant for
  gene therapy, where CpGs trigger TLR9 responses).
- **`run_optimization_pipeline`** - chain the above in order.

### Alternative splicing (`biocompiler.ir.splicing`)

Splicing is modeled as a Non-deterministic Finite State Transducer (NDFST):
one pre-mRNA can produce multiple mature mRNAs. `splice_ndfst(ir_l1)` returns
the primary isoform plus alternative isoforms (exon skipping, intron
retention, no-splicing), each scored by MaxEntScan splice-site strength.

### Frontend (`biocompiler.ir.frontend`)

```yaml
# gene.yaml
gene_name: HBB
organism: human
sequence: ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA
regions:
  - type: cds
    start: 0
    end: 96
```

### Codegen (`biocompiler.ir.codegen`)

- `to_genbank(ir_l0)` - GenBank flat file with region annotations.
- `to_fasta_protein(ir_l3)` / `to_fasta_mrna(ir_l2)` - FASTA (60-col wrap).
- `to_sbol3(ir_l0)` - SBOL3 (synthetic biology standard).

Coordinates convert from IR's 0-indexed half-open `[start, end)` to GenBank's
1-indexed inclusive `start..end` (`gb_start = ir_start + 1`, `gb_end = ir_end`).

---

## Formal Verification Status (honest)

The Lean4 development lives in [`proof/`](proof/). It builds cleanly with
`lake build` (Lean toolchain managed via `elan`).

### What is proved

- **17 core predicates** - fully proved sound in Lean4 (0 `sorry` in core;
  0 `sorry` in the whole development after the W1-A5 BLOSUM62 discharge).
  Soundness theorem (for these 17): if `evaluate(P, seq, C) = PASS`, then
  `propertyHolds(P, seq, C)`. This is a real, machine-checked guarantee -
  scoped to those 17 predicates, not to the whole type system. These 17 are
  the **deterministic layer** — every predicate whose verdict can be decided
  by sequence computation alone (no external tools).
  The full Lean4 development totals **267 theorems** across 17 modules
  (8,746 LOC), including the new `BLOSUM62.lean` module (11 theorems
  formalizing the BLOSUM62 substitution matrix).
- **19 SLOT predicates** - vacuously sound in `CONSERVATIVE` mode: by design
  they never return `PASS` (they return `UNCERTAIN` and require external tool
  evidence to escalate). They cannot produce a false PASS. These 19 are
  **oracle-dependent** — they require external tools (ESMFold, NetMHCpan,
  ViennaRNA) and are handled by SLOT, which proves they cannot produce a
  false PASS without external tool evidence.
- **All lowering passes** (`transcribe`, `splice`, `translate`) - proved
  correct in Lean4 (`proof/BioCompiler/IR.lean`,
  `proof/BioCompiler/SplicingResolution.lean`).
- **Compositional soundness** - composing proved predicates preserves
  soundness under five-valued logic.
- **SLOT independence** - FFI / external-tool calls never produce a `PASS`
  verdict by construction.

> In short: **deterministic safety is formally guaranteed; oracle-dependent
> safety is conservatively bounded.** This separation is the contribution.

### What is not yet proved

- **0 `sorry` remain in `SLOTVerification.lean`.** The 2 former
  BLOSUM62-related `sorry` were discharged in W1-A5 by formalizing the
  BLOSUM62 substitution matrix in `proof/BioCompiler/BLOSUM62.lean`
  (11 theorems proving symmetry, diagonal non-negativity, and range bounds).
- **34 narrowed `axiom` declarations in `SLOTVerification.lean`** (narrowed
  from 15 former broad axioms, one per tool). These are tool-interface
  soundness contracts for external ML models (TMHMM, ViennaRNA, ESMFold,
  FoldX, CamSol, Aggrescan, ExPASy, NetMHCpan, BePiPred, IEDB). Each asserts
  a single specific testable property of an external tool's output
  (window-size matching, threshold-range validity, proxy correctness, or a
  soundness contract). Each is backed by a runtime evidence check in
  `src/biocompiler/provenance/runtime_evidence.py` (34 checks total, one per
  narrowed axiom). Discharging them as theorems would require formalizing
  the external tools themselves, which is out of scope.
- **The trusted computing base (TCB)** is 3 class-field axioms in
  `SpliceSiteScanner` (down from an initial 18) — documented in
  [`proof/README.md`](proof/README.md).

### Predicate counts (corrected)

| Tier | Lean4 model | Python implementation |
|------|-------------|------------------------|
| Core (fully proved) | 17 | 17 |
| SLOT (vacuously sound in conservative mode) | 19 | 19 |
| Extended diagnostic (Python only, not in Lean4) | - | 7 |
| **Total** | **36** | **43** |

> Canonical counts: 17 core + 19 SLOT = 36 in Lean4; 43 in Python (with 7
> extended diagnostic predicates that are runtime-only and not formalized).

> **The 17/19 separation, stated honestly.** The 17 core predicates are the
> **deterministic layer** — every predicate whose verdict can be decided by
> sequence computation alone (no external tools). The 19 SLOT predicates are
> **oracle-dependent** — they require external tools (ESMFold, NetMHCpan,
> ViennaRNA) and are handled by SLOT, which proves they **cannot produce a
> false PASS** without external tool evidence. This separation is the
> contribution: deterministic safety is formally guaranteed; oracle-dependent
> safety is conservatively bounded.

### Build

```bash
cd proof/
lake build    # Requires Lean4 (via elan)
```

A successful build confirms all proofs machine-check with **0 `sorry`**
(after the W1-A5 BLOSUM62 discharge), 3 TCB class-field axioms, and **34
narrowed tool-soundness `axiom`s** in `SLOTVerification.lean` for
oracle-dependent predicates (each backed by a runtime evidence check in
`src/biocompiler/provenance/runtime_evidence.py`).

---

## Benchmark Results

### Heavy Fair Benchmark vs DNAchisel (25 genes, 5 timed runs each)

Both tools optimize the SAME protein with the SAME constraints (GC 30-70%,
avoid 4 restriction enzymes, no premature stops) and SAME objective (maximize
CAI). 3 warmup + 5 timed runs (median). See [`BENCHMARKS.md`](BENCHMARKS.md)
for full methodology and per-gene results.

| Metric | BioCompiler | DNAchisel |
|---|---|---|
| Mean CAI | 0.872 | **0.9737** |
| Median time | 18.38ms | **10.19ms** |
| CAI wins | 10/25 | **15/25** |
| Speed wins | 10/25 | **15/25** |
| Violations | 0 | 0 |

**Honest summary**: DNAchisel wins on raw CAI (0.9737 vs 0.872) and speed
(10.19ms vs 18.38ms) overall. BioCompiler wins on E. coli (CAI 1.000 vs 0.948,
2-4× faster) but loses on human (CAI 0.792 vs 0.992, 7× slower) on the
default path. **With the opt-in context-aware GT-avoidance mode**
(`use_context_aware_gt=True`), BioCompiler achieves CAI 0.963 on human HBB
(vs. 0.718 for the default path on the same sequence) by repairing only
high-scoring cryptic splice donors rather than forbidding every GT
dinucleotide. BC's value is not raw speed/CAI — it's the certified-by-default
path (predicate evaluation, certificate generation, biosecurity screening,
provenance trail) and eukaryotic safety constraints (context-aware GT
avoidance). For raw throughput without certification, use DNAchisel.

> **Performance note.** Biosecurity fuzzy matching is JIT-compiled with
> numba when available, reducing screening overhead by ~30%. MaxEntScan
> splice-site scoring uses the real Yeo & Burge 2004 trained parameters
> (not approximations).

### Large-Scale E2E Test (16,616 UniProt proteins × 6 organisms = 99,696 combos)

This is the E2E slice of the 200,064-combination test sweep (see Test Results
below for the full breakdown).

| Metric | Value |
|---|---|
| Proteins tested | 16,616 (all reviewed human proteins, 50-800aa) |
| Organisms | 6 (human, mouse, yeast, E. coli, CHO, Pichia) |
| Total combinations (this slice) | **99,696** |
| PASS | **99,696 (100.0%)** |
| FAIL | 0 |
| TIMEOUT | 0 |
| ERROR | 0 |

All 25 known human selenoproteins are verified end-to-end in
`tests/test_selenoproteins_e2e.py` (130 tests). The catalog in
`src/biocompiler/selenoproteins.py` records the verified Sec positions from
the canonical selenoprotein literature (Kryukov 2003; Gladyshev 2016;
UniProtKB feature tables). For each selenoprotein, the test constructs a
fragment containing U at a verified position, back-translates to DNA with
TGA at the Sec codon, runs the full L0→L3 pipeline with SECIS annotation,
and verifies that U is preserved (not converted to stop). A negative
control verifies that without SECIS annotation, the same TGA correctly
translates to stop. An additional 10 SECIS unit tests in `tests/test_secis.py`
cover UGA recoding, multiple SECIS positions, and end-to-end back-translation
on synthetic cases.

### Integrated Constraint-Solving Optimizer

BioCompiler uses a single-pass integrated constraint-solving optimizer that selects the best codon at each position satisfying all constraints. On E. coli it achieves CAI 1.000 (vs DNAchisel's 0.948) and is 2-4× faster. On human proteins it is slower (certified-by-default overhead) and has lower CAI (GT-avoidance constraint). See BENCHMARKS.md for the full fair comparison.

## Test Results

BioCompiler is tested at scale. Three test regimes cover the IR compiler,
the end-to-end optimizer->IR->codegen pipeline, and the head-to-head
comparison.

| Regime | Combinations tested | Result |
|--------|--------------------:|--------|
| IR compiler (L0->L1->L2->L3, including selenocysteine) | 100,068 | **100% pass** |
| E2E (optimizer + IR + codegen, 16,616 proteins × 6 organisms) | 99,696 | **100% pass** |
| Head-to-head vs DNAchisel (50 proteins × 6 organisms) | 300 | **100% pass** |
| **Total** | **200,064** | **100% pass (200,064/200,064 combinations, including selenocysteine)** |

All 150 selenocysteine-containing combinations (25 proteins × 6 organisms)
now pass after the integrated-optimizer SeC fix (commit b9f079d). There are
zero confirmed correctness failures across the 200,064 combinations.

Additional coverage:

- **Biosecurity screening** - 8 of 8 reference hazard fragments correctly flagged (ricin A, shiga, diphtheria, tetanus, cholera, abrin, botulinum, anthrax LF). The botulinum and anthrax LF signatures were added by GAP-1; all 8 select-agent toxins now block optimization. Verified by `tests/test_biosecurity_coverage.py`.
- **Immunogenicity honesty** - the immunogenicity layer returns `UNCERTAIN`
  when its underlying scores cannot be backed by a real tool call (e.g.
  NetMHCpan unavailable). It does not fabricate confidence. Immunogenicity
  analysis is available as a **separate API** (`compute_immunogenicity()`,
  `deimmunize()`) that integrates NetMHCpan/MHCflurry when installed. It is
  **NOT part of the default `optimize_sequence()` path** because it requires
  external tools and adds ~seconds of latency. The `LowImmunogenicity`
  predicate is a SLOT predicate that returns `UNCERTAIN` when tools are
  unavailable.
- **Standalone verifier** - a 473-line, stdlib-only Python verifier
  (`scripts/standalone_verifier.py`) re-checks certificate claims without
  importing the rest of BioCompiler. This is distinct from the
  `biocompiler-verify` console entry point, which is the *in-process*
  verifier (it imports BioCompiler and runs the full `check` pipeline).
  Use `python scripts/standalone_verifier.py <cert.json>` when you need a
  TCB-independent re-check; use `biocompiler-verify <cert.json>` for the
  convenience in-process re-check.

---

## Features

### Compiler

- **Multi-level IR** - five typed levels (L0-L4) mirroring the Central
  Dogma, with invariant checks at each boundary.
- **YAML frontend** - `parse_spec` / `compile_from_spec` accept a file path
  or inline YAML string.
- **Codegen backends** - GenBank, FASTA, SBOL3.
- **Selenocysteine (U) support** - `secis_positions` on IR-L0/L1/L2 cause
  `UGA` to be translated as `U` rather than stop, matching biology of
  selenoproteins.
- **NDFST alternative splicing** - enumerates exon-skip / intron-retention /
  no-splice isoforms from one pre-mRNA, scored by MaxEntScan.
- **ESMFold integration** - `fold` (L3->L4) calls ESMFold when available and
  falls back to a pure-Python Chou-Fasman heuristic; the oracle used is
  recorded on the IR-L4 object.
- **IR optimization passes** - `optimize_codons`, `eliminate_cpgs`, with
  protein-preservation verification after every pass.

### Optimizer (pre-IR)

- **Integrated constraint-solving optimizer** - single-pass, selects the
  best codon at each position satisfying ALL constraints simultaneously.
  Operates on raw protein->DNA; the IR optimization passes are thin adapters
  over this. Reach it via `optimize_sequence(...)` (the only optimizer since
  0.9.1; the legacy slow-path stack was removed).
- **30 curated organisms** - codon usage tables, organism-aware constraint
  profiles (prokaryotes skip splice/CpG constraints that don't apply).
- **Two CAI reference sets** - Kazusa (25 organisms) and Sharp-Li (E. coli,
  reproduces the original 1987 values for literature comparison).
- **tAI metric** - tRNA-adaptation index for 10 organisms.
- **Custom objectives, multi-gene constructs, part libraries, assembly
  planning (Golden Gate / Gibson), pattern enforcement, IUPAC ambiguity
  resolution.**

### Safety

- **Biosecurity screening on by default** - input sequences screened against
  pathogen / toxin signature databases before any optimization work.
- **Protein verification on by default** - optimized DNA is back-translated
  and checked against the target protein.
- **Strict mode (opt-in)** - every warning becomes a hard failure; intended
  for clinical / production pipelines.
- **Sliding-window GC** - detects local GC extremes that global averages
  hide (default 50 bp window).
- **Provenance tracking** - every optimization decision is recorded for
  audit and reproducibility.

### Verification & certificates

- **Graduated certificates** - Gold / Silver / Bronze tiers, with the tier
  capped by the worst verdict in the predicate set (any `UNCERTAIN` -> not
  Gold; multiple `UNCERTAIN` -> not Silver).
- **Standalone verifier** - 473 LOC, stdlib-only, re-checks a certificate
  without trusting the BioCompiler package. Invoke as
  `python scripts/standalone_verifier.py <cert.json>`; this is distinct
  from `biocompiler-verify`, the in-process verifier that imports the
  BioCompiler package.
- **Five-valued verdicts** - `PASS / LIKELY_PASS / UNCERTAIN / LIKELY_FAIL
  / FAIL`, with `PASS/FAIL/UNCERTAIN` modelled in Lean4.

---

## What It Is NOT

Honesty matters in safety-adjacent software. The following are explicitly
**not** claimed:

- **Not clinically validated.** No wet-lab, no animal model, no human trial.
  The 200,064-test sweep is a software correctness argument, not a
  biological validation.
- **Not the fastest optimizer in every scenario.** On the heavy fair benchmark (25 genes, same constraints), DNAchisel wins overall: CAI 0.973 vs 0.872, 10.4ms vs 17.6ms. BioCompiler wins on E. coli (CAI 1.000, 2-4× faster) but loses on human on the default path (BC applies GT-avoidance that reduces CAI). The opt-in context-aware GT-avoidance mode (`use_context_aware_gt=True`) recovers most of the lost CAI — on human HBB it achieves CAI 0.963 vs. 0.718 for the default path — by repairing only high-scoring cryptic splice donors rather than forbidding every GT dinucleotide. The certified-by-default path adds ~10-50ms overhead. For raw throughput, use DNAchisel.
- **Not a general-purpose bioinformatics tool.** BioCompiler does not do
  alignment, variant calling, RNA-seq, phylogenetics, or assembly. It is a
  gene-design compiler: it takes a target protein and emits a gene design
  with a provenance trail and a certificate.
- **Not a complete formal proof.** 17 core predicates are proved; 19 SLOT
  predicates are vacuously sound; **0 `sorry` remain in `SLOTVerification.lean`**
  (the 2 former BLOSUM62-related `sorry` were discharged in W1-A5 by
  formalizing the BLOSUM62 substitution matrix in
  `proof/BioCompiler/BLOSUM62.lean`, 11 theorems). The former 15 broad
  tool-soundness `axiom` declarations have been **narrowed to 34 specific,
  independently-testable contracts**, each backed by a runtime evidence check
  in `src/biocompiler/provenance/runtime_evidence.py` (34 checks total).
  These narrowed axioms serve as the closed-form TCB extension for
  oracle-dependent predicates (ESMFold, FoldX, NetMHCpan correctness
  contracts). Each is tracked as `NEEDS_TOOL_AXIOM`, not silently skipped.
  The proof **degrades gracefully**: in conservative mode, the system never
  produces a false PASS even with these obligations open. Do not claim
  "machine-verified correctness" for the whole system — claim it for the 17
  core predicates (the **deterministic layer**), and be explicit about the 34
  narrowed tool-soundness axioms (the **oracle-dependent layer**, which is
  conservatively bounded but not fully formally guaranteed).

---

## Related Work

BioCompiler builds on a substantial body of prior art. We acknowledge it
here rather than claim novelty we don't have:

- **Beal & BBN - "Proto BioCompiler" (2011).** A compiler-style framework
  for genetic engineering with typed IRs and verification. BioCompiler's
  overall shape - typed IR, staged lowering, verification - is in the same
  family. We do **not** claim to be the first formal verification of a type
  system for biological sequence analysis; that prior work exists.
- **Nielsen et al. - Cello (2016).** A digital-circuit-to-DNA compiler for
  prokaryotic gene circuits. Cello targets regulatory logic, not
  protein-coding gene optimization, but the compiler analogy is the same.
- **Zulkower & Rosser - DNAchisel (2020).** A flexible, fast Python library
  for codon optimization with a constraint-solving API. DNAchisel is the
  performance baseline for BioCompiler's optimizer (see "What it is NOT").
- **Lin et al. - ESMFold (2023).** The structure-prediction oracle used at
  IR-L4. Treated as an external tool whose output is evidence, not a proven
  property (SLOT-style).
- **Yeo & Burge - MaxEntScan (2004).** First-order Markov model for splice
  site scoring, used by the NDFST splicing module.
- **Sharp & Li - CAI (1987).** The Codon Adaptation Index; BioCompiler
  reproduces the original E. coli reference values for literature
  comparison.

---

## Quick Start

### Install

```bash
pip install -e .

# Optional extras:
pip install -e ".[all]"          # All optional dependencies
pip install -e ".[compare]"      # DNAchisel comparison harness
```

### Compile a gene through the IR pipeline

```python
from biocompiler.ir.frontend import compile_from_spec
from biocompiler.ir.types import IRLevel
from biocompiler.ir.codegen import to_genbank, to_fasta_protein
from biocompiler.ir.optimization import run_optimization_pipeline

# 1. Parse the YAML frontend -> IR-L0, then lower all the way to IR-L3 (protein).
ir_l3 = compile_from_spec("src/biocompiler/ir/example_specs/hbb.yaml",
                          target_level=IRLevel.L3)
print(f"Protein: {ir_l3.sequence}")
# -> MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*

# 2. Stop at IR-L2 (mature mRNA) and run IR optimization passes on it.
ir_l2 = compile_from_spec("src/biocompiler/ir/example_specs/hbb.yaml",
                          target_level=IRLevel.L2)
optimized = run_optimization_pipeline(ir_l2)
print(f"Passes applied: {optimized.metadata['passes_applied']}")
print(f"CAI after: {optimized.metadata['cai_after']:.4f}")

# 3. Emit standard file formats from any IR level.
ir_l0 = compile_from_spec("src/biocompiler/ir/example_specs/hbb.yaml",
                          target_level=IRLevel.L0)
print(to_genbank(ir_l0)[:200])          # GenBank flat file
print(to_fasta_protein(ir_l3))          # FASTA protein
```

The YAML spec can also be passed inline as a string:

```python
ir_l3 = compile_from_spec(
    "gene_name: test\norganism: e_coli\nsequence: ATGGCTAAATGGCGTTAA\n",
    target_level=IRLevel.L3,
)
# ir_l3.sequence == "MAKWR*"
```

### Standalone optimizer (pre-IR, raw protein -> DNA)

```python
from biocompiler.optimizer.pipeline_core import optimize_sequence

result = optimize_sequence(
    protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGR",
    organism="Homo_sapiens",
    strict_mode=True,           # warnings become hard failures
)
print(f"DNA:      {result.sequence}")
print(f"CAI:      {result.cai:.4f}")
print(f"GC:       {result.gc_content:.2%}")
print(f"Verified: {result.verification_passed}")
print(f"Biosec:   {result.biosecurity_screening_result}")
```

### CLI / REST API

```bash
biocompiler optimize --protein "MVHLTPEEK..." --organism Homo_sapiens --enable-mutagenesis
uvicorn biocompiler.api:app --host 0.0.0.0 --port 8000    # API docs at /docs
```

### Building the proof

```bash
cd proof/
lake build    # Requires Lean4 (via elan)
```

---

## Supported Organisms

BioCompiler ships with 30 curated organisms with full codon usage tables and
organism-aware constraint profiles. The most common:

| Organism | Canonical name | Short key | Domain |
|----------|---------------|-----------|--------|
| *E. coli* | `Escherichia_coli` | `ecoli` | Prokaryote |
| *B. subtilis* | `Bacillus_subtilis` | `bacillus` | Prokaryote |
| *H. sapiens* | `Homo_sapiens` | `human` | Eukaryote |
| *M. musculus* | `Mus_musculus` | `mouse` | Eukaryote |
| *S. cerevisiae* | `Saccharomyces_cerevisiae` | `yeast` | Eukaryote |
| *K. phaffii* (Pichia) | `Komagataella_phaffii` | `pichia` | Eukaryote |
| CHO-K1 | `CHO_K1` | `cho` | Eukaryote |
| HEK-293T | `HEK293T` | `hek293` | Eukaryote |
| *D. rerio* | `Danio_rerio` | `danio` | Eukaryote |
| *D. melanogaster* | `D_melanogaster` | `drosophila` | Eukaryote |
| *A. thaliana* | `Arabidopsis_thaliana` | `arabidopsis` | Eukaryote |

Full list (30 organisms, including NS0, PER.C6, Cricetulus, Spodoptera,
Trichoplusia, and 8 more) in
[`docs/15-Reference.md`](docs/15-Reference.md). The domain is auto-detected
from the organism name; override with `organism_domain='auto'|'eukaryote'|'prokaryote'`.
Additional organisms can be loaded dynamically from the Kazusa Codon Usage
Database:

```python
from biocompiler.organisms import resolve_or_download_organism
resolve_or_download_organism('Thermus_thermophilus')
```

tAI is available for 10 organisms: E. coli, human, mouse, yeast, CHO,
C. elegans, Drosophila, Arabidopsis, Pichia, Bacillus.

---

## CAI Reference Sets

BioCompiler supports two CAI reference sets:

| Set | Source | Coverage | Best for |
|------|--------|----------|----------|
| **Kazusa** | Kazusa Codon Usage Database (high-expression subsets) | 25 organisms | General optimization, cross-organism comparisons |
| **Sharp-Li** | Sharp & Li (1987), 24 highly expressed E. coli genes | E. coli only | Comparing against published CAI values |

The Sharp-Li set reproduces the exact relative-adaptiveness values from the
original CAI paper. For low-expression genes like lacZ, the two sets can
differ by >0.4 CAI - see [`docs/reference_sets.md`](docs/reference_sets.md)
for guidance.

---

## Performance

The integrated optimizer achieves CAI >= 0.99 on standard benchmark genes
(GFP, HBB, insulin) for both prokaryotic and eukaryotic targets. It is a
single-pass O(n×k) algorithm: it selects the best codon at each position
satisfying ALL constraints simultaneously, with no distinct phases.

The integrated optimizer achieves CAI 1.000 on E. coli (vs DNAchisel's 0.948) and is 2-4× faster on prokaryotes. On eukaryotes, DNAchisel is faster and achieves higher CAI (no GT-avoidance constraint). See BENCHMARKS.md for the full fair
benchmark (see Benchmark Results above). The full IR pipeline, invariant
checking, provenance, and certificate generation add real overhead to the
end-to-end pipeline. The trade is deliberate.

---

## Architecture

```
YAML spec --> frontend --> IR-L0 --> transcribe --> IR-L1 --> splice --> IR-L2 --> translate --> IR-L3 --> fold --> IR-L4
                              |                  |                 |                  |               |
                          GenBank            (NDFST           (SECIS-aware      (ESMFold      FASTA protein
                          SBOL3          alternative splicing)  UGA->U at Sec)   + heuristic)   (mRNA + protein)
```

The optimizer is the integrated constraint-solving optimizer
(single-pass, O(n×k)); since 0.9.1 it is the only optimizer. The optimizer
sits upstream of the IR: it takes a target protein and produces a DNA
sequence, which then enters the IR as IR-L0 for verification, certificate
generation, and codegen.

### Repository Layout

```
biocompiler/
+-- proof/                     # Lean4 formalization (17 modules, 0 sorry, 3 TCB + 34 narrowed tool-soundness axioms, 267 theorems)
+-- src/biocompiler/
|   +-- ir/                    # Multi-level IR (types, passes, codegen, frontend, splicing, folding, optimization)
|   +-- optimizer/             # Integrated constraint-solving optimizer (protein -> DNA) — the only optimizer since 0.9.1
|   +-- type_system/           # 43 predicates (17 core proved, 19 SLOT, 7 extended diagnostic)
|   +-- biosecurity/           # Pathogen / toxin screening
|   +-- immunogenicity/        # NetMHCpan / MHCflurry / deimmunization
|   +-- engines/               # External tool integration (ESMFold, FoldX, CamSol, ViennaRNA)
|   +-- organisms/             # 30 curated organism configs + codon tables
|   +-- export/                # GenBank / SBOL2 / SBOL3 exporters
|   +-- provenance/            # Decision tracking + graduated certificates
|   +-- solver/                # Constraint satisfaction infrastructure (conflict resolution, MUS diagnosis)
|   +-- sequence/              # Scanners, splicing, sliding GC, IUPAC, MaxEntScan
|   +-- expression/            # CAI / tAI / Sharp-Li / mRNA stability
|   +-- api/                   # FastAPI REST API
|   +-- cli/                   # Command-line interface
|   +-- shared/                # Five-valued logic, types, constants, exceptions
+-- scripts/
|   +-- standalone_verifier.py # 473-LOC stdlib-only certificate verifier
+-- tests/                     # 200,064-combination test sweep
+-- docs/                      # SE spec (20 docs + 18 ADRs)
+-- papers/                    # Manuscripts
```

---

## Documentation

| Document | Content |
|----------|---------|
| [DOC-00](docs/00-README.md) | Document map and reading order |
| [DOC-02](docs/02-SAD.md) | Software Architecture (ISO 42010) |
| [DOC-06](docs/06-Design-Rationale.md) | Design rationale and critical analysis |
| [DOC-10](docs/10-Deterministic-Methods.md) | Six deterministic methods for non-deterministic biology |
| [DOC-14](docs/14-SLOT-Proof-Implementation-Gap.md) | SLOT predicate proof-implementation gap (0 sorry + 34 narrowed tool-soundness axioms) |
| [DOC-15](docs/15-Reference.md) | 43-predicate tables, engine API, TCB, limitations |
| [Reference Sets](docs/reference_sets.md) | Kazusa vs Sharp-Li CAI reference sets |
| [docs/adr/](docs/adr/) | 18 Architecture Decision Records |
| [proof/README.md](proof/README.md) | Lean4 TCB inventory (3 class-field axioms + 34 narrowed tool-soundness axioms) |

---

## Citation

```bibtex
@misc{biocompiler2026,
  title     = {BioCompiler: A Formally-Verified Gene Design Compiler
               with Multi-Level Intermediate Representations},
  author    = {Khairkhah, Parham},
  year      = {2026},
  note      = {Beta. 17/36 Lean4 predicates fully proved; 0 sorry (all discharged,
              BLOSUM62 formalized in BLOSUM62.lean); 34 narrowed tool-soundness
              axioms in SLOTVerification.lean; 267 theorems total;
               279-gene reference test sweep; not clinically validated},
  url       = {https://github.com/pkhairkh/biocompiler}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).
