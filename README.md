# BioCompiler v10.0.0

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

## What's New in v10

### HybridOptimizer — 3-Phase Gene Optimization

The **HybridOptimizer** is a new high-performance optimizer combining greedy initialization, priority-based local search, and CAI hill climbing:

| Phase | Strategy | Purpose |
|-------|----------|---------|
| **Phase 1** | Greedy Initialization | Back-translate with highest-CAI codons |
| **Phase 2** | Priority-Queue Local Search | Fix constraint violations by severity with incremental re-evaluation |
| **Phase 3** | CAI Hill Climbing | Aggressively recover CAI while maintaining all constraints |

Key innovation: instead of sequential constraint resolution that can undo previous fixes, the HybridOptimizer uses a **priority queue** and **incremental constraint evaluation** to avoid the "fix A → break B → fix B → break A" oscillation loop. After each fix, only affected positions are re-evaluated.

Performance target: **<1.5ms** for GFP (714bp), **CAI > 0.98**.

### CAI Table Unification (Breaking Change)

v10 unifies the CAI computation tables across all optimizer backends. Previously, different optimizers could use slightly different codon adaptiveness values, leading to inconsistent CAI results. Now, all paths use the same unified table from `organisms.CODON_ADAPTIVENESS_TABLES`.

**Migration impact**: CAI values may differ slightly from v9. In particular, v9 CAI values for some genes were **inflated** due to a table mismatch — v10 values are the correct ones.

### `species` → `organism` Parameter Migration

The `species` parameter (e.g., `'ecoli'`, `'human'`) is now deprecated in favor of the more explicit `organism` parameter (e.g., `'Escherichia_coli'`, `'Homo_sapiens'`). Both forms still work, and `species` will continue to be supported for backward compatibility.

```python
# Both forms work in v10:
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

# Both parameter forms work:
result = optimize_sequence(protein, species='ecoli')
result = optimize_sequence(protein, organism='Escherichia_coli')

print(f"Optimized DNA: {result.sequence}")
print(f"CAI: {result.cai:.4f}")
print(f"GC: {result.gc_content:.2%}")
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
| **Type-directed mutagenesis** — V→I substitutions make HBB feasible (BLOSUM62=+3) | Proof of concept |
| **SE specification** — 11 IEEE/ISO-standard documents + 14 ADRs | Complete |

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
│   ├── hybrid_optimizer.py   # HybridOptimizer (new in v10)
│   ├── optimization.py       # BioOptimizer + optimize_sequence API
│   ├── api.py                # FastAPI REST API
│   └── ...
├── tests/              # Test suite (420+ tests)
├── docs/               # Full SE specification (14 docs + 14 ADRs)
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
  note      = {v10.0.0 — HybridOptimizer, unified CAI tables, organism parameter, 28-predicate type system},
  url       = {https://github.com/pkhairkh/biocompiler}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).
