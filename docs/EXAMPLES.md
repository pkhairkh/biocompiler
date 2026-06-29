# BioCompiler Examples

Five runnable, end-to-end scripts. Each one was tested against
BioCompiler 0.9.0 before it was written down — copy, paste, and run.

> **Beta disclaimer.** BioCompiler is not clinically validated. Treat
> every output as a design candidate to be reviewed by a qualified
> molecular biologist before synthesis.

---

## Example 1: Optimize HBB for Human Expression

Human hemoglobin beta (HBB), N-terminal 31 residues, optimized for
*Homo sapiens* expression. Uses non-strict mode because the input is a
short fragment (see the [User Guide troubleshooting section](USER_GUIDE.md#norqctrigger-failure-strict-mode)).

```python
"""
Example 1 — Optimize HBB (human hemoglobin beta) for human expression.
"""
from biocompiler import optimize_sequence

PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
ORGANISM = "Homo_sapiens"

result = optimize_sequence(
    target_protein=PROTEIN,
    organism=ORGANISM,
    strict_mode=False,          # short fragment — accept partial result
    track_provenance=True,
)

# --- Quality metrics -------------------------------------------------------
print("=== HBB / Homo_sapiens ===")
print(f"DNA length:     {len(result.sequence)} nt")
print(f"CAI:            {result.cai:.3f}    (1.0 = perfect adaptation)")
print(f"GC content:     {result.gc_content:.1%}")
print(f"IR verified:    {result.ir_verified}  (True = protein preserved)")
print(f"Timed out:      {result.timed_out}")

# --- Predicate results -----------------------------------------------------
print("\n=== Predicate Results ===")
for p in result.predicate_results:
    status = "PASS" if p.passed else "FAIL"
    print(f"  [{status:4}] {p.predicate}: {p.details[:80]}")

if result.failed_predicates:
    print(f"\nFailed predicates: {result.failed_predicates}")
if result.warnings:
    print(f"\nWarnings:")
    for w in result.warnings:
        print(f"  - {w}")

# --- File-format exports ---------------------------------------------------
print("\n=== Exports ===")
print(f"GenBank:  {len(result.genbank)} chars")
print(f"FASTA:    {len(result.fasta)} chars")
print(f"SBOL3:    {len(result.sbol3)} chars")

# Write them to disk if you want to hand them to a synthesis provider:
# with open("hbb_optimized.gb", "w") as f: f.write(result.genbank)
# with open("hbb_optimized.fasta", "w") as f: f.write(result.fasta)
# with open("hbb_optimized.sbol3.ttl", "w") as f: f.write(result.sbol3)

# --- Sanity check ----------------------------------------------------------
# `optimize_sequence()` already back-translates the optimized DNA through
# the IR pipeline (L0->L1->L2->L3) and confirms the translation matches
# the input protein. The result is exposed on the OptimizationResult:
assert result.ir_verified is True, "optimizer failed IR verification"
# `result.ir_l3.sequence` is the protein the optimized DNA translates to.
# It always starts with the input protein and ends with '*' (stop codon).
print(f"\nBack-translation matches input protein: {result.ir_l3.sequence!r}")
assert result.ir_l3.sequence.startswith(PROTEIN)
```

---

## Example 2: Optimize Insulin for E. coli

The human preproinsulin signal + B chain, optimized for high-yield *E. coli*
expression. Uses non-strict mode with a tighter-than-default GC band
(40%-60%): the prokaryote fast path skips eukaryote-only constraints, but
the tight GC band makes the sliding-window check harder to satisfy on a
90 bp sequence, so we accept the partial result and inspect
`failed_predicates`.

```python
"""
Example 2 — Optimize human insulin B chain for E. coli expression.
"""
from biocompiler import optimize_sequence

# Insulin B chain (UniProt P01313, residues 25-54 of preproinsulin)
INSULIN_B = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"

# Default constraints (GC 30%-70%, default enzyme panel): strict_mode passes
# cleanly. Tightening the GC band to 40%-60% leaves less room for the
# optimizer to satisfy the sliding-window GC check on short sequences, so
# we use strict_mode=False here and inspect failed_predicates afterwards.
result = optimize_sequence(
    target_protein=INSULIN_B,
    organism="e_coli",
    strict_mode=False,
    gc_lo=0.40, gc_hi=0.60,     # E. coli optimal GC band (tighter than default)
    enzymes=["NdeI", "XhoI", "BamHI", "HindIII"],  # common cloning enzymes
)

print("=== Insulin B / E. coli ===")
print(f"DNA:        {result.sequence}")
print(f"Protein:    {result.protein}")
print(f"CAI:        {result.cai:.3f}")
print(f"GC:         {result.gc_content:.1%}")
print(f"Verified:   {result.ir_verified}")
print(f"Failed:     {result.failed_predicates}")
# Typical output: CAI ~0.79, GC ~62%, ir_verified=True, failed=['SlidingGC']
# (the tight 40-60% band is hard to satisfy in every 50 bp window for a
# 90 bp sequence — accept the result and inspect failed_predicates)

# Restriction-site check: confirm the optimized DNA has zero sites for the
# enzymes we asked it to avoid.
from biocompiler.type_system import check_no_restriction_site
sites_check = check_no_restriction_site(
    result.sequence, enzymes=["NdeI", "XhoI", "BamHI", "HindIII"],
)
print(f"\nRestriction sites: {sites_check.passed}  ({sites_check.details})")

# Write the FASTA for the synthesis provider.
fasta_block = result.fasta
print(f"\n--- FASTA ---\n{fasta_block}")
```

---

## Example 3: Verify a Gene Therapy Construct

A full verification pipeline for a designed gene-therapy construct:
(1) compile the DNA through the IR to confirm the translation,
(2) run the 43-predicate type system,
(3) screen for biosecurity hazards.

```python
"""
Example 3 — Verify a gene therapy construct end-to-end.

The example DNA is the natural HBB N-terminus (96 nt). In a real workflow
you would substitute your own designed DNA and your expected protein.
"""
from biocompiler.ir.types import IR_L0_GenomicDNA, IRLevel
from biocompiler.ir.passes import compile_gene
from biocompiler.ir.codegen import to_genbank, to_sbol3
from biocompiler.type_system import evaluate_all_predicates
from biocompiler.biosecurity import screen_hazardous_sequence

# --- Inputs ----------------------------------------------------------------
DESIGNED_DNA = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG"
    "AACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA"
)
EXPECTED_PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"
ORGANISM = "human"

# --- Step 1: IR pipeline (DNA -> protein) ----------------------------------
ir_l0 = IR_L0_GenomicDNA(
    sequence=DESIGNED_DNA,
    regions=[],                 # no intron annotations — start/stop scan
    organism=ORGANISM,
    gene_name="HBB-therapy",
)
ir_l3 = compile_gene(ir_l0, IRLevel.L3)

print("=== Step 1: IR compilation ===")
print(f"Translated: {ir_l3.sequence}")
if ir_l3.sequence != EXPECTED_PROTEIN:
    print(f"MISMATCH! Expected: {EXPECTED_PROTEIN}")
    raise SystemExit(1)
print("OK — translation matches expected protein.")

# --- Step 2: 43-predicate type system -------------------------------------
print("\n=== Step 2: Type system (43 predicates) ===")
results = evaluate_all_predicates(DESIGNED_DNA, organism=ORGANISM)

verdict_counts = {}
for r in results:
    v = r.verdict.value
    verdict_counts[v] = verdict_counts.get(v, 0) + 1

for v in ("PASS", "LIKELY_PASS", "UNCERTAIN", "LIKELY_FAIL", "FAIL"):
    if v in verdict_counts:
        print(f"  {v:<12} {verdict_counts[v]:>3}")

failures = [r for r in results if r.verdict.value in ("FAIL", "LIKELY_FAIL")]
if failures:
    print("\nFailures (investigate before synthesis):")
    for r in failures:
        print(f"  - {r.predicate}: {r.violation}")
else:
    print("\nNo FAIL / LIKELY_FAIL verdicts.")

# --- Step 3: Biosecurity screening ----------------------------------------
print("\n=== Step 3: Biosecurity screening ===")
report = screen_hazardous_sequence(
    protein=ir_l3.sequence.rstrip("*"),   # screener expects AA-only input
    dna=DESIGNED_DNA,
)
if report.is_hazardous:
    print(f"  !!! HAZARD DETECTED — risk level: {report.risk_level}")
    for m in report.matches:
        print(f"      [{m.category}] {m.name} (pos {m.position}, "
              f"conf={m.confidence:.2f}, type={m.match_type})")
    for rec in report.recommendations:
        print(f"      -> {rec}")
else:
    print(f"  OK — no signatures matched (risk_level={report.risk_level})")
    print("  Note: 'no signatures matched' is NOT the same as 'safe'.")

# --- Step 4: Emit synthesis-ready files ------------------------------------
print("\n=== Step 4: File exports ===")
print(f"GenBank: {len(to_genbank(ir_l0))} chars")
print(f"SBOL3:   {len(to_sbol3(ir_l0))} chars")

# Save the GenBank file for your cloning / synthesis provider:
# with open("hbb_therapy.gb", "w") as f: f.write(to_genbank(ir_l0))
```

---

## Example 4: Batch Optimize Multiple Genes

Optimize a panel of proteins against the same host organism and produce
a single multi-FASTA output. We iterate over `optimize_sequence()` rather
than `batch_optimize()` because the per-protein API also populates the
IR objects, the GenBank / FASTA / SBOL3 exports, and the IR-verification
flag — the batch API does not.

```python
"""
Example 4 — Batch-optimize a panel of proteins for E. coli expression.
"""
import time
from biocompiler import optimize_sequence

GENE_PANEL = [
    # (gene_name, uniprot_id, protein_sequence)
    ("HBB",       "P68871", "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"),
    ("INS-B",     "P01313", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"),
    ("GFP",       "P42212", "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYG"),
    ("mCherry",   "X5DSL3", "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGR"),
    ("SUMO1",     "P63165", "MSDQEAKPSTEDLGDKKEGEYIKLKVIGQDSSEIHFKVKMT"),
]
ORGANISM = "e_coli"

results = []
t0 = time.time()
for name, uniprot, protein in GENE_PANEL:
    r = optimize_sequence(
        target_protein=protein,
        organism=ORGANISM,
        strict_mode=True,
    )
    results.append((name, uniprot, protein, r))
    print(f"  {name:<10}  CAI={r.cai:.3f}  GC={r.gc_content:.1%}  "
          f"ir_ok={r.ir_verified}  failed={r.failed_predicates}")
print(f"\nTotal wall-clock: {time.time() - t0:.1f}s "
      f"({len(GENE_PANEL)} proteins)")

# --- Write a multi-FASTA for the synthesis provider -----------------------
with open("gene_panel_optimized.fasta", "w") as fh:
    for name, uniprot, protein, r in results:
        fh.write(f">{name}|{uniprot}|{ORGANISM}|len={len(r.sequence) // 3}\n")
        # Wrap DNA at 60 cols (FASTA convention).
        for i in range(0, len(r.sequence), 60):
            fh.write(r.sequence[i:i + 60] + "\n")
print("Wrote gene_panel_optimized.fasta")

# --- Summary table ---------------------------------------------------------
print("\n=== Summary ===")
print(f"{'Gene':<10} {'CAI':>6} {'GC%':>6} {'IR ok':>6} {'Failed':>20}")
for name, _, _, r in results:
    print(f"{name:<10} {r.cai:>6.3f} {r.gc_content*100:>5.1f}% "
          f"{str(r.ir_verified):>6} {','.join(r.failed_predicates) or '-':>20}")
```

---

## Example 5: Compare with DNAchisel

Head-to-head comparison of BioCompiler and [DNAchisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)
on the same protein and constraint set. Requires the optional
`compare` extra:

```bash
pip install biocompiler[compare]
```

If DNAchisel is not installed, the script runs BioCompiler only and prints
a reminder instead of crashing.

```python
"""
Example 5 — Head-to-head comparison: BioCompiler vs DNAchisel.

The adapter translates BioCompiler's constraint model into DNAchisel's
specification, runs both tools on the same input, and reports CAI / GC /
restriction-site count side-by-side. CAI is always recomputed with
BioCompiler's evaluator so the metric is comparable across tools.
"""
import time
from biocompiler import optimize_sequence
from biocompiler.expression.translation import compute_cai
from biocompiler.benchmarking.dnachisel_adapter import (
    DNAchiselAdapter,
    is_dnachisel_available,
)

PROTEIN = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
ORGANISM = "Homo_sapiens"
ENZYMES = ["EcoRI", "BamHI", "XhoI", "HindIII"]

# --- BioCompiler -----------------------------------------------------------
t0 = time.time()
bc_result = optimize_sequence(
    target_protein=PROTEIN,
    organism=ORGANISM,
    enzymes=ENZYMES,
    strict_mode=False,
)
bc_time = time.time() - t0

print("=== BioCompiler ===")
print(f"  CAI:        {bc_result.cai:.4f}")
print(f"  GC:         {bc_result.gc_content:.1%}")
print(f"  IR verified:{bc_result.ir_verified}")
print(f"  Time:       {bc_time:.3f}s")

# --- DNAchisel (optional) --------------------------------------------------
if not is_dnachisel_available():
    print("\n=== DNAchisel ===")
    print("  DNAchisel is not installed. Install with: pip install dnachisel")
    print("  (or: pip install biocompiler[compare])")
else:
    adapter = DNAchiselAdapter()
    t0 = time.time()
    dc_result = adapter.optimize(
        protein=PROTEIN,
        organism=ORGANISM,
        constraints=[
            {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
            {"type": "avoid_restriction", "enzymes": ENZYMES},
        ],
    )
    dc_time = time.time() - t0

    print("\n=== DNAchisel ===")
    print(f"  CAI:        {dc_result.cai:.4f}   (recomputed with BioCompiler)")
    print(f"  GC:         {dc_result.gc_content:.1%}")
    print(f"  Rest sites: {dc_result.restriction_site_count}")
    print(f"  Time:       {dc_result.execution_time_s:.3f}s")

    # --- Side-by-side summary ---------------------------------------------
    print("\n=== Comparison ===")
    print(f"{'Metric':<20} {'BioCompiler':>14} {'DNAchisel':>14}")
    print(f"{'CAI':<20} {bc_result.cai:>14.4f} {dc_result.cai:>14.4f}")
    print(f"{'GC content':<20} {bc_result.gc_content:>13.1%} "
          f"{dc_result.gc_content:>13.1%}")
    print(f"{'Wall-clock (s)':<20} {bc_time:>14.3f} "
          f"{dc_result.execution_time_s:>14.3f}")
    print()
    print("Note: BioCompiler is roughly 14x faster than DNAchisel on")
    print("standard benchmark genes, while additionally providing IR invariant")
    print("checking, provenance tracking, and certificate generation. See the")
    print("README section 'What It Is NOT' for the full performance disclaimer.")
```

---

## Where to Go Next

- **[USER_GUIDE.md](USER_GUIDE.md)** — the practical user guide with
  installation, the four core use cases, the full organism table, and
  troubleshooting.
- **[../README.md](../README.md)** — the technical README: IR pipeline,
  formal verification status, performance numbers, "What It Is NOT".
- **[21-IR-Compiler.md](21-IR-Compiler.md)** — the IR pipeline in depth
  (L0->L4, lowering passes, codegen, alternative splicing).
- **[15-Reference.md](15-Reference.md)** — the 43-predicate reference,
  engine API, TCB inventory, known limitations.
