# BioCompiler — Papers

This directory contains the two BioCompiler papers.

## Paper Inventory

| Directory | Title | Target Venue | Status |
|-----------|-------|-------------|--------|
| `technical_report/` | BioCompiler: A Formally-Verified Gene Design Compiler with Multi-Level IR | arXiv | Complete (13 pages, VLM-verified) |
| `main/` | Well-Typed Genes Don't Go Wrong | VSTTE | Draft (most complete, 50 sections) |

## `technical_report/` — arXiv Technical Report

The honest, scoped claim document. 438 lines of LaTeX, compiles to a 13-page
PDF (`main.pdf`). Covers: five-level IR pipeline, integrated constraint-solving
optimizer, Lean4 formal verification (17 modules, 15 open sorry), 1,674-combination
test sweep, biosecurity screening. Referenced in the project README as the
arXiv technical report.

```bash
cd papers/technical_report/
tectonic main.tex   # or: pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## `main/` — VSTTE Submission Draft

The flagship paper: "Well-Typed Genes Don't Go Wrong." 3,592 lines of LaTeX
across 4 files (`main.tex` + 3 section fragments). Targets VSTTE (Verified
Software: Theories, Tools, Experiments) — a formal-methods venue appropriate
for the Lean4 soundness proof contribution. Compiles to a 398KB PDF.

```bash
cd papers/main/
tectonic main.tex   # or: pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## `demonstrations/`

Five demonstration scripts (`demo1`–`demo5`) that exercise the BioCompiler API
on real proteins (therapeutic autopsy, adversarial protein, multi-organism
atlas, COVID vaccine, PCR primer compatibility). All 5 pass. See
`demonstrations/run_all_demonstrations.py`.
