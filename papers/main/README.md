# Paper Directory — BioCompiler (VSTTE submission)

## Current Draft

**`main.tex`** — Canonical VSTTE submission draft.

- Document class: `acmart` (ACM format, anonymous review)
- Title: "Well-Typed Genes Don't Go Wrong: A Type System for Machine-Verified Gene Design"
- 43 type predicates across five biological domains (DNA-level, structure, stability, solubility, immunogenicity)
- SLOT (Sealed Local Oracle Theory) framework with refinement theorem
- Lean 4 soundness proof: 17 modules, 8,746 LOC, 2 sorry in SLOTVerification.lean (NEEDS_TOOL_AXIOM), 0 sorry in other 16 modules
- Self-contained: all sections inline (no `\input`/`\include`)

## Venue Rationale

This paper targets **VSTTE** (Verified Software: Theories, Tools, Experiments),
not POPL. The contribution is a machine-checked soundness proof in Lean 4 plus a
five-valued type system for codon-level gene design constraints — a formal-methods
contribution rather than a programming-language theory contribution. There is no
frontend programming language (the user-facing input is a YAML config), no type
inference algorithm, and no code generation to a machine target (GenBank/FASTA/SBOL3
are data interchange formats, not machine code). VSTTE is the natural home for
verified-tools papers of this shape.

## Section Fragments

The following files contain extracted sections for parallel editing. They are **not** included by `main.tex` (which is self-contained). They duplicate content already present in `main.tex` and may be used for collaborative review or split editing.

| File | Sections | Notes |
|------|----------|-------|
| `sections_1_3.tex` | 1 (Introduction), 2 (Motivating Example), 3 (The Type System) | Fragment — requires preamble from main document |
| `sections_4_6.tex` | 4 (SLOT Framework & Refinement), 5 (Formal Verification), 6 (Certificates & Mutagenesis) | Standalone compilable — has own preamble |
| `sections_7_9.tex` | 7 (Implementation), 8 (Evaluation), 9 (Related Work) | Fragment — requires preamble from main document |

## Supporting Files

| File | Purpose |
|------|---------|
| `references.bib` | BibTeX bibliography (paper-specific entries; self-contained bibliography) |
| `texmf/acmart.cls` | ACM document class |
| `texmf/mathpartir.sty` | Inference rule typesetting |

## Build Instructions

```bash
cd papers/main/
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Requires `acmart.cls` and `mathpartir.sty` (both in `texmf/`).

## Consolidation History

This directory was consolidated from the former `paper/` directory at the
project root, then renamed from `popl2027/` to `main/` when the venue target
was changed from POPL 2027 to VSTTE. File renames during consolidation:

- `paper/popl2027.tex` → `papers/popl2027/popl2027.tex` → `papers/main/main.tex`
- `paper/popl_sections_1_3.tex` → `papers/popl2027/sections_1_3.tex` → `papers/main/sections_1_3.tex`
- `paper/popl_sections_4_6.tex` → `papers/popl2027/sections_4_6.tex` → `papers/main/sections_4_6.tex`
- `paper/popl_sections_7_9.tex` → `papers/popl2027/sections_7_9.tex` → `papers/main/sections_7_9.tex`
