# BioCompiler v12.0.0

A compiler framework for human protein synthesis using intermediate representations — with a machine-verified soundness proof.

---

## What This Is

BioCompiler applies compiler-engineering design patterns — staged transformation, typed intermediate representations, composable passes, and type-based verification — to the formalizable stages of gene-to-protein processing: RNA splicing and translation.

The central result is a **machine-verified soundness proof** in Lean4: the type system never produces a false PASS verdict. If the type checker says a gene design is correct, it *is* correct — under explicitly stated assumptions.

```
type_soundness: ∀ P, ∀ seq, ∀ C, evaluate(P, seq, C) = PASS → propertyHolds(P, seq, C)
```

This is, to our knowledge, the first formal verification of a type system governing biological sequence analysis.

Where existing gene design tools provide scores, heuristics, or ML predictions, BioCompiler provides a categorical answer — not a better heuristic, but a different *kind* of answer — using three-valued logic (PASS / FAIL / UNCERTAIN), type predicates enforcing biological constraints, and compositional soundness that preserves guarantees without independence assumptions.

For safety-critical gene design — gene therapy, clinical-grade synthetic biology — a false positive can kill a patient. No existing gene design tool provides a machine-verified guarantee that its output is correct. BioCompiler does.

---

## What's New in v12

### Biosecurity Screening

Automated screening of input protein sequences against known pathogen-associated and toxin-related databases. Every optimization now runs biosecurity checks by default, flagging sequences with homology to dangerous agents before any gene synthesis work begins.

```python
result = optimize_sequence(protein, organism='Escherichia_coli')
print(f"Biosecurity flags: {result.biosecurity_flags}")
```

### Protein Verification

Post-optimization protein verification ensures the optimized DNA sequence translates back to the exact target protein. This catches silent bugs caused by codon-table edge cases, IUPAC ambiguity resolution, or mutagenesis side effects.

```python
result = optimize_sequence(protein, organism='Homo_sapiens')
assert result.protein_verified  # Always True unless mutagenesis was applied
```

### Strict Mode

A new `strict_mode` parameter makes every warning a hard failure. In strict mode, any unresolved constraint violation, uncertain predicate, or biosecurity flag causes the optimizer to raise an exception rather than return a potentially unsafe result. Recommended for clinical and production pipelines.

```python
result = optimize_sequence(protein, organism='Homo_sapiens', strict_mode=True)
```

### Sliding-Window GC Content

In addition to global GC%, BioCompiler now computes GC content over a configurable sliding window (default 50 bp). This catches local GC extremes — AT-rich or GC-rich islands — that global averages hide and that can cause expression failures in vivo.

```python
result = optimize_sequence(protein, organism='Escherichia_coli', gc_window_size=50)
print(f"Max local GC: {result.max_local_gc:.2%}")
print(f"Min local GC: {result.min_local_gc:.2%}")
```

### Custom Objectives

Users can now define custom optimization objectives beyond the built-in CAI, GC, and constraint predicates. Custom objectives are scored alongside native objectives and influence codon selection during the hill-climbing phase.

```python
from biocompiler.objectives import CustomObjective

my_obj = CustomObjective(name="avoid_homopolymers", weight=0.5, evaluate_fn=...)
result = optimize_sequence(protein, organism='Escherichia_coli', custom_objectives=[my_obj])
```

### Multi-Gene Constructs

First-class support for optimizing multiple genes as a single construct — with shared promoter/terminator regions, per-gene constraint tracking, and construct-level GC and restriction-site management.

```python
from biocompiler.multigene import optimize_construct

result = optimize_construct(
    genes={"gene_a": protein_a, "gene_b": protein_b},
    organism='Homo_sapiens',
)
```

### IUPAC Ambiguity Support

Input sequences may now contain IUPAC ambiguity codes (e.g., `N`, `R`, `Y`). BioCompiler resolves ambiguity codes deterministically during back-translation, selecting the highest-CAI codon compatible with each ambiguous position.

### Pattern Enforcement

Arbitrary sequence patterns can be enforced or avoided in the optimized output. This generalizes restriction-site avoidance to any user-defined motif — Kozak sequences, TATA boxes, cryptic promoter elements, etc.

```python
result = optimize_sequence(
    protein, organism='Homo_sapiens',
    avoid_patterns=['TATAAA', 'TATATA'],
    enforce_patterns=[('Kozak', 'GCCACC')],
)
```

### Part Libraries

A composable part library system for standard biological parts (promoters, RBS, terminators, tags). Parts are assembled into constructs with automatic constraint propagation across junction boundaries.

```python
from biocompiler.parts import PartLibrary, AssemblySpec

lib = PartLibrary.from_default()
parts = [lib.get("T7_promoter"), lib.get("His_tag"), lib.get("T7_terminator")]
```

### Assembly Planning

Automated assembly plan generation given a set of parts and a target construct. BioCompiler selects optimal restriction enzymes, generates fragment boundaries, and produces step-by-step assembly protocols compatible with Golden Gate, Gibson, and traditional restriction-ligation methods.

```python
from biocompiler.assembly import plan_assembly

plan = plan_assembly(parts, method='golden_gate')
print(plan.protocol_steps)
```

### SBOL Support

Import and export designs in SBOL (Synthetic Biology Open Language) format. This enables interoperability with SBOL-compliant tools and registries.

```python
from biocompiler.sbol_export import export_sbol
from biocompiler.sbol_import import import_sbol

# Export
sbol_doc = export_sbol(result, gene_name="GFP")

# Import
design = import_sbol("design.xml")
```

### LIMS Integration

A built-in LIMS (Laboratory Information Management System) integration module for tracking optimization runs, storing results, and managing sample metadata. Supports SQLite (default) and PostgreSQL backends.

```python
from biocompiler.lims import LIMSClient

lims = LIMSClient()  # SQLite by default
run_id = lims.record_optimization(result, metadata={"project": "gene_therapy_v2"})
```

### tAI Metric

Transfer RNA Adaptation Index (tAI) is now available as an alternative to CAI. tAI accounts for tRNA gene copy numbers and wobble pairing rules, providing a more biologically grounded measure of codon optimality for organisms with well-characterized tRNA pools.

```python
result = optimize_sequence(protein, organism='Homo_sapiens', metric='tai')
print(f"tAI: {result.tai:.4f}")
```

---

## What Was New in v10

### HybridOptimizer — 3-Phase Gene Optimization

The **HybridOptimizer** is a high-performance optimizer combining greedy initialization, priority-based local search, and CAI hill climbing:

| Phase | Strategy | Purpose |
|-------|----------|---------|
| **Phase 1** | Greedy Initialization | Back-translate with highest-CAI codons |
| **Phase 2** | Priority-Queue Local Search | Fix constraint violations by severity with incremental re-evaluation |
| **Phase 3** | CAI Hill Climbing | Aggressively recover CAI while maintaining all constraints |

Key innovation: instead of sequential constraint resolution that can undo previous fixes, the HybridOptimizer uses a **priority queue** and **incremental constraint evaluation** to avoid the "fix A → break B → fix B → break A" oscillation loop. After each fix, only affected positions are re-evaluated.

Performance target: **<1.5ms** for GFP (714bp), **CAI > 0.98**.

### CAI Table Unification (Breaking Change in v10)

v10 unified the CAI computation tables across all optimizer backends. Previously, different optimizers could use slightly different codon adaptiveness values, leading to inconsistent CAI results. All paths now use the same unified table from `organisms.CODON_ADAPTIVENESS_TABLES`.

**Migration impact**: CAI values may differ slightly from v9. In particular, v9 CAI values for some genes were **inflated** due to a table mismatch — v10 values are the correct ones.

### `species` → `organism` Parameter Migration

The `species` parameter (e.g., `'ecoli'`, `'human'`) is deprecated in favor of the more explicit `organism` parameter (e.g., `'Escherichia_coli'`, `'Homo_sapiens'`). Both forms still work, and `species` will continue to be supported for backward compatibility.

```python
# Both forms work:
result = optimize_sequence(protein, species='ecoli')
result = optimize_sequence(protein, organism='Escherichia_coli')
```

The `organism` parameter also accepts short keys, abbreviated binomials, and display names, which are resolved automatically:

| Input Form | Example | Resolved To |
|-----------|---------|-------------|
| Short key | `'ecoli'` | `Escherichia_coli` |
| Abbreviated binomial | `'E_coli'` | `Escherichia_coli` |
| Display name | `'E. coli'` | `Escherichia_coli` |
| Canonical name | `'Escherichia_coli'` | `Escherichia_coli` |

### E. coli Codon Usage Data Correction

Codon usage data for 5 amino acids in E. coli has been corrected to match the latest Kazusa Codon Database values. This affects CAI computation for Ala, Arg, Gly, Leu, and Val.

---

## Performance

### Optimization Benchmarks

| Gene | Organism | CAI | Time | vs DNAchisel |
|------|----------|:---:|:----:|:------------:|
| GFP | E. coli | 0.999 | 2ms | 2.4× |
| HBB | E. coli | 1.000 | 3ms | — |
| GFP | Human | 0.983 | 15ms | — |
| Insulin | E. coli | 0.998 | 1ms | 3.0× |
| HBB | Human | 0.97 | 18ms | — |

The HybridOptimizer achieves **2–3× speedup** over DNAchisel on standard benchmark genes while producing equal or higher CAI values.

### HybridOptimizer Phase Breakdown

```
GFP (714bp, E. coli):
  Phase 1 (Greedy):     CAI 1.000
  Phase 2 (Constraints): CAI 0.999  (1 violation fixed)
  Phase 3 (Hill Climb):  CAI 0.999  (0 improvements needed)
  Total: 2ms

GFP (714bp, Human):
  Phase 1 (Greedy):     CAI 0.997
  Phase 2 (Constraints): CAI 0.975  (8 violations fixed)
  Phase 3 (Hill Climb):  CAI 0.983  (3 improvements)
  Total: 15ms
```

---

## Quick Start

### Install

```bash
pip install -e .
```

For optional dependencies:

```bash
pip install -e ".[optimizer]"    # Z3 constraint solver
pip install -e ".[all]"          # All optional dependencies
pip install -e ".[compare]"      # DNAchisel comparison
```

### Python API

```python
from biocompiler.api import optimize_sequence

# Recommended: strict_mode + biosecurity (default in v12)
result = optimize_sequence(
    protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR...",
    organism='Escherichia_coli',
    strict_mode=True,          # Fail fast on any constraint violation
)

print(f"Optimized DNA: {result.sequence}")
print(f"CAI: {result.cai:.4f}")
print(f"GC: {result.gc_content:.2%}")
print(f"Protein verified: {result.protein_verified}")
print(f"Biosecurity flags: {result.biosecurity_flags}")

# Legacy parameter forms still work:
result = optimize_sequence(protein, species='ecoli')
result = optimize_sequence(protein, organism='Escherichia_coli')
```

### Full Optimization with Certificates

```python
from biocompiler import optimize_sequence, generate_certificate, verify_certificate
from biocompiler import evaluate_all_predicates

# Optimize HBB for human expression with mutagenesis
result = optimize_sequence(
    target_protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR...",
    organism="Homo_sapiens",
    enable_mutagenesis=True,
)

print(f"CAI: {result.cai:.4f}")
print(f"Satisfied: {result.satisfied_predicates}")
print(f"Mutagenesis: {result.mutagenesis_applied}")

# Generate and verify a certificate
results = evaluate_all_predicates(seq=result.sequence, known_exon_boundaries=[(0, len(result.sequence))])
cert = generate_certificate(result.sequence, results, {"gene": "HBB", "organism": "Homo_sapiens"})
status, failures = verify_certificate(cert.to_dict())
print(f"Certificate: {status}")
```

### Using the HybridOptimizer Directly

```python
from biocompiler import HybridOptimizer

optimizer = HybridOptimizer(
    species='ecoli',                   # or organism='Escherichia_coli'
    enzymes=['EcoRI', 'BamHI', 'XhoI'],
    gc_lo=0.30,
    gc_hi=0.70,
    avoid_gt=True,                     # Skip for prokaryotic targets
)

result = optimizer.optimize(protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR")

print(f"CAI: {result.cai:.4f}")
print(f"GC:  {result.gc_content:.2%}")
print(f"Phase 1 CAI: {result.phase1_cai:.4f}")
print(f"Phase 2 CAI: {result.phase2_cai:.4f}")
print(f"Phase 3 CAI: {result.phase3_cai:.4f}")
print(f"Violations fixed: {result.violations_fixed}")
print(f"Hill climb improvements: {result.hill_climb_improvements}")
```

### CLI / REST API

```bash
biocompiler optimize --protein "MVHLTPEEK..." --organism Homo_sapiens --enable-mutagenesis
uvicorn biocompiler.api:app --host 0.0.0.0 --port 8000
```

### Building the Proof

```bash
cd proof/
lake build    # Requires Lean4 (via elan)
```

A successful build confirms all proofs machine-check with 0 `sorry` and 0 axioms.

---

## Safety Features

BioCompiler v12 ships with **safety-by-default** — new features that protect against common gene-design hazards are enabled out of the box:

| Feature | Default | What It Does |
|---------|---------|-------------|
| **Biosecurity screening** | ✅ On | Screens input sequences against pathogen/toxin databases before optimization |
| **Protein verification** | ✅ On | Verifies optimized DNA translates back to the target protein |
| **Strict mode** | Off (opt-in) | Makes every warning a hard failure — ideal for clinical pipelines |
| **Sliding-window GC** | ✅ On (50 bp) | Detects local GC extremes that global averages miss |
| **Restriction-site avoidance** | ✅ On | Removes recognition sites for common restriction enzymes |
| **Cryptic splice-site avoidance** | ✅ On (eukaryotes) | Avoids GT/AG donor/acceptor motifs in eukaryotic targets |
| **CpG island avoidance** | ✅ On (eukaryotes) | Reduces CpG dinucleotides to prevent epigenetic silencing |
| **Provenance tracking** | ✅ On | Records every optimization decision for audit and reproducibility |

### Why Safety-By-Default?

In gene therapy and clinical synthetic biology, a silently unsafe optimization — wrong protein, pathogen homology, hidden splice site — can have catastrophic consequences. BioCompiler's type system already provides a machine-verified soundness guarantee; the v12 safety defaults extend that philosophy from the formal model into the engineering defaults.

```python
# Safety defaults are active even in the simplest call:
result = optimize_sequence(protein, organism='Homo_sapiens')

# For maximum safety in clinical use:
result = optimize_sequence(
    protein, organism='Homo_sapiens',
    strict_mode=True,           # Fail on any unresolved issue
    gc_window_size=30,          # Tighter local GC window
)
```

---

## Breaking Changes in v10

### 1. CAI Values Now Correctly Computed

CAI values were inflated in v9 due to a table mismatch between the optimizer and the CAI computation module. v10 unifies the tables, so reported CAI values are now accurate. **If you relied on v9 CAI values, expect small decreases** — the new values are the correct ones.

### 2. `species` Parameter Deprecated in Favor of `organism`

The `species` parameter (e.g., `'ecoli'`) still works but is deprecated. Use `organism` (e.g., `'Escherichia_coli'`) for new code. If both are provided, `species` takes precedence for backward compatibility, but a `DeprecationWarning` is emitted.

### 3. E. coli Codon Usage Data Corrected for 5 Amino Acids

Codon usage frequencies for Ala, Arg, Gly, Leu, and Val in E. coli have been updated to match the latest Kazusa Codon Database. This may change optimization results for E. coli targets compared to v9.

---

## Key Results

| Result | Status |
|--------|--------|
| **Soundness proof** — 0 `sorry`, 0 axioms, 5 TCB assumptions | Lean4 formalization |
| **Compositional soundness** — composing predicates preserves soundness under three-valued logic | Proved in Lean4 |
| **SLOT independence** — FFI calls never produce PASS | Proved in Lean4 |
| **Per-predicate soundness** — 13 core proved sound + 19 SLOT-dependent proved UNCERTAIN | Proved in Lean4 |
| **Production pipeline** — Scanner → NDFST → Translation → TypeCheck → Certificate → Verify | Python (420+ tests) |
| **HBB full pass** — 8 optimizer predicates pass simultaneously | HBB CAI=0.97 |
| **28-predicate type system** — 12 DNA + 4 structure + 4 stability + 4 solubility + 4 immunogenicity | 13 core + 19 SLOT |
| **HybridOptimizer** — 3-phase optimizer with priority-based local search | CAI 0.999, 2ms (GFP) |
| **Biosecurity screening** — pathogen/toxin homology checks on by default | v12 safety-by-default |
| **Protein verification** — post-opt back-translation check on by default | v12 safety-by-default |
| **Strict mode** — warnings become hard failures | v12 opt-in |
| **Sliding-window GC** — local GC extremes detected automatically | v12 on by default |
| **tAI metric** — tRNA-based codon optimality alternative to CAI | v12 |
| **Multi-gene constructs** — per-gene constraint tracking in shared construct | v12 |
| **IUPAC support** — ambiguity codes resolved deterministically | v12 |
| **Part libraries + assembly planning** — Golden Gate, Gibson, restriction-ligation | v12 |
| **SBOL import/export** — interoperability with SBOL tools and registries | v12 |
| **LIMS integration** — SQLite and PostgreSQL backends for run tracking | v12 |
| **Type-directed mutagenesis** — V→I substitutions make HBB feasible (BLOSUM62=+3) | Proof of concept |
| **SE specification** — 11 IEEE/ISO-standard documents + 16 ADRs | Complete |

---

## Architecture

```
DNA Sequence → Scanner → NDFST (Splicing) → Translation → TypeCheck → Certificate → Verify
                  ↓                          ↓              ↓
              IR-Seq tokens            IR-Peptide     PASS/FAIL/UNCERTAIN
```

### HybridOptimizer Pipeline

```
Protein Sequence
       ↓
  Phase 1: Greedy Init ──── Back-translate with highest-CAI codons
       ↓
  Phase 2: Constraint Fix ── Priority queue: fix most severe violations first
       ↓                       Incremental re-evaluation after each fix
  Phase 3: CAI Hill Climb ── Upgrade codons while maintaining constraints
       ↓
  Optimized DNA Sequence
```

The pipeline processes gene sequences through typed intermediate representations. Core predicates (13) evaluate deterministically; SLOT-dependent predicates (19) delegate to external tools and return UNCERTAIN in the formal model. Type-directed mutagenesis proposes conservative amino acid substitutions when no codon assignment satisfies all predicates.

### Repository Structure

```
biocompiler/
├── proof/              # Lean4 soundness proof
├── src/biocompiler/    # Production Python package
│   ├── hybrid_optimizer.py   # HybridOptimizer (3-phase)
│   ├── optimization.py       # BioOptimizer + optimize_sequence API
│   ├── biosecurity.py        # Pathogen/toxin screening (v12)
│   ├── protein_verification.py # Post-opt verification (v12)
│   ├── sliding_gc.py         # Sliding-window GC (v12)
│   ├── objectives.py         # Custom objectives (v12)
│   ├── multigene.py          # Multi-gene constructs (v12)
│   ├── iupac.py              # IUPAC ambiguity support (v12)
│   ├── pattern_enforcement.py # Pattern enforce/avoid (v12)
│   ├── parts.py              # Part libraries (v12)
│   ├── assembly.py           # Assembly planning (v12)
│   ├── sbol_export.py        # SBOL export (v12)
│   ├── sbol_import.py        # SBOL import (v12)
│   ├── lims.py               # LIMS integration (v12)
│   ├── api.py                # FastAPI REST API
│   └── ...
├── scripts/            # Benchmark and utility scripts
├── tests/              # Test suite (420+ tests)
├── docs/               # Full SE specification (14 docs + 16 ADRs)
└── paper/              # LaTeX manuscript
```

---

## CAI Reference Sets

BioCompiler supports two CAI (Codon Adaptation Index) reference sets for computing codon optimality:

| Reference Set | Source | Coverage | Best For |
|---------------|--------|----------|----------|
| **Kazusa** | Kazusa Codon Usage Database (high-expression subsets) | 5 organisms (E. coli, H. sapiens, S. cerevisiae, M. musculus, CHO-K1) | General-purpose optimization, cross-organism comparisons |
| **Sharp-Li** | Sharp & Li (1987), 24 highly expressed E. coli genes | E. coli only | Comparing against published literature values, CAIcal validation |

The Sharp-Li reference set reproduces the exact relative adaptiveness values from the original CAI paper, enabling direct comparison with published CAI values. For low-expression genes like lacZ, the two sets can differ by >0.4 CAI — see [`docs/reference_sets.md`](docs/reference_sets.md) for detailed guidance on when to use each.

```python
from biocompiler.benchmarking import run_benchmark_by_name

# Validate CAI against published Sharp & Li (1987) values
result = run_benchmark_by_name("sharp_li_cai")
print(f"Sharp-Li closer to published: {result['sharp_li_is_closer']}")
```

### Organism-Aware Constraint Selection

BioCompiler automatically disables biologically inappropriate constraints based on the target organism's domain:

- **Prokaryotic targets** (E. coli): Cryptic splice-site avoidance and CpG-island avoidance are disabled — prokaryotes lack spliceosomes and CpG methylation
- **Eukaryotic targets** (H. sapiens, CHO-K1, etc.): All constraints are enabled

This recovers ~0.27 CAI on average for prokaryotic targets:

```python
from biocompiler.benchmarking import run_benchmark_by_name

result = run_benchmark_by_name("organism_aware_cai")
print(f"Mean CAI recovery: {result['mean_cai_recovery']:+.4f}")
```

### Provenance System

Every optimization decision is tracked with full provenance — including which reference set was used for CAI computation, which constraints were active, and what alternatives were considered. This makes BioCompiler's gene designs reproducible and auditable, a key differentiator for safety-critical applications.

---

## Supported Organisms

| Organism | Canonical Name | Short Key | Domain |
|----------|---------------|-----------|--------|
| *Escherichia coli* | `Escherichia_coli` | `ecoli` | Prokaryote |
| *Homo sapiens* | `Homo_sapiens` | `human` | Eukaryote |
| *Mus musculus* | `Mus_musculus` | `mouse` | Eukaryote |
| *Saccharomyces cerevisiae* | `Saccharomyces_cerevisiae` | `yeast` | Eukaryote |
| CHO-K1 | `CHO_K1` | `cho` | Eukaryote |

The domain is auto-detected from the organism name. You can override it with the `organism_domain` parameter (`'auto'`, `'eukaryote'`, `'prokaryote'`).

---

## Documentation

Full technical documentation is in [`docs/`](docs/):

| Document | Content |
|----------|---------|
| [DOC-00](docs/00-README.md) | Document map and reading order |
| [DOC-02](docs/02-SAD.md) | Software Architecture (ISO 42010) |
| [DOC-06](docs/06-Design-Rationale.md) | Design Rationale and critical analysis of the original proposal |
| [DOC-10](docs/10-Deterministic-Methods.md) | Six deterministic methods for non-deterministic biology |
| [DOC-14](docs/14-SLOT-Proof-Implementation-Gap.md) | SLOT predicate proof-implementation gap |
| [DOC-15](docs/15-Reference.md) | Technical reference: 28-predicate tables, engine API, TCB, limitations |
| [Reference Sets](docs/reference_sets.md) | CAI reference sets (Kazusa vs Sharp-Li): when to use each, example usage, expected differences |
| [docs/adr/](docs/adr/) | 14 Architecture Decision Records |

---

## Citation

```bibtex
@misc{biocompiler2026,
  title     = {BioCompiler: A Compiler Framework for Human Protein Synthesis
               Using Intermediate Representations with a Machine-Verified Soundness Proof},
  author    = {Khairkhah, Pouya},
  year      = {2026},
  note      = {v12.0.0 — Safety-by-default, biosecurity screening, protein verification, strict mode, multi-gene, tAI},
  url       = {https://github.com/pkhairkh/biocompiler}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).
