# Paper Directory — BioCompiler POPL 2027

## Current Draft

**`popl2027.tex`** — Canonical POPL 2027 submission draft.

- Document class: `acmart` (ACM format, anonymous review)
- Title: "Well-Typed Genes Don't Go Wrong: A Type System for Machine-Verified Gene Design"
- 33 type predicates across five biological domains (DNA-level, structure, stability, solubility, immunogenicity)
- SLOT (Sealed Local Oracle Theory) framework with refinement theorem
- Lean 4 soundness proof: 18 modules, 7,930 LOC, zero sorry, zero axioms
- Self-contained: all sections inline (no `\input`/`\include`)

## Section Fragments

The following files contain extracted sections for parallel editing. They are **not** included by `popl2027.tex` (which is self-contained). They duplicate content already present in `popl2027.tex` and may be used for collaborative review or split editing.

| File | Sections | Notes |
|------|----------|-------|
| `popl_sections_1_3.tex` | 1 (Introduction), 2 (Motivating Example), 3 (The Type System) | Fragment — requires preamble from main document |
| `popl_sections_4_6.tex` | 4 (SLOT Framework & Refinement), 5 (Formal Verification), 6 (Certificates & Mutagenesis) | Standalone compilable — has own preamble |
| `popl_sections_7_9.tex` | 7 (Implementation), 8 (Evaluation), 9 (Related Work) | Fragment — requires preamble from main document |

## Supporting Files

| File | Purpose |
|------|---------|
| `references.bib` | BibTeX bibliography (all citations) |
| `texmf/acmart.cls` | ACM document class (POPL format) |
| `texmf/mathpartir.sty` | Inference rule typesetting |

## Removed Files (2025-03-04 cleanup)

| File | Reason |
|------|--------|
| `icfp_draft.tex` | Old ICFP 2026 draft (deadline passed), superseded by POPL 2027 draft |
| `main.tex` | Older draft with 7 predicates and llncs class, superseded by popl2027.tex |
| `icfp_draft.pdf` | Stale compiled output |
| `main.pdf` | Stale compiled output |
| `icfp_draft.aux/.log/.out` | LaTeX build artifacts |

## Build Instructions

```bash
cd paper/
pdflatex popl2027.tex
bibtex popl2027
pdflatex popl2027.tex
pdflatex popl2027.tex
```

Requires `acmart.cls` and `mathpartir.sty` (both in `texmf/`).
