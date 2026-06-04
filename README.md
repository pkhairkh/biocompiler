# BioCompiler v9.0.0

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
| **Type-directed mutagenesis** — V→I substitutions make HBB feasible (BLOSUM62=+3) | Proof of concept |
| **SE specification** — 11 IEEE/ISO-standard documents + 14 ADRs | Complete |

---

## Quick Start

### Install

```bash
pip install -e .
```

### Python API

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

## Architecture

```
DNA Sequence → Scanner → NDFST (Splicing) → Translation → TypeCheck → Certificate → Verify
                  ↓                          ↓              ↓
              IR-Seq tokens            IR-Peptide     PASS/FAIL/UNCERTAIN
```

The pipeline processes gene sequences through typed intermediate representations. Core predicates (13) evaluate deterministically; SLOT-dependent predicates (19) delegate to external tools and return UNCERTAIN in the formal model. Type-directed mutagenesis proposes conservative amino acid substitutions when no codon assignment satisfies all predicates.

### Repository Structure

```
biocompiler/
├── proof/              # Lean4 soundness proof
├── src/biocompiler/    # Production Python package
├── tests/              # Test suite (420+ tests)
├── docs/               # Full SE specification (14 docs + 14 ADRs)
└── paper/              # LaTeX manuscript
```

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
| [docs/adr/](docs/adr/) | 14 Architecture Decision Records |

---

## Citation

```bibtex
@misc{biocompiler2026,
  title     = {BioCompiler: A Compiler Framework for Human Protein Synthesis
               Using Intermediate Representations with a Machine-Verified Soundness Proof},
  author    = {Khairkhah, Pouya},
  year      = {2026},
  note      = {v9.0.0 — unified engine API, 28-predicate type system, SLOT architecture, HBB full pass},
  url       = {https://github.com/pkhairkh/biocompiler}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).
