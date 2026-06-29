# BioCompiler User Guide

> **Status: Beta.** BioCompiler is not yet validated in wet-lab or clinical
> settings. Treat every output as a *design candidate* to be reviewed by a
> qualified molecular biologist before synthesis.

This guide shows you how to do the four things most users want to do with
BioCompiler: **optimize** a gene, **verify** a design, **screen** for
biosecurity hazards, and **compile** from a YAML spec. Every snippet below
was tested against BioCompiler 0.9.0 before it was written down.

For the formal verification story, the IR pipeline, and the architecture,
read the [technical README](../README.md) and the
[IR compiler doc](21-IR-Compiler.md). This document is deliberately
practical: no theory, just "how do I do X".


## Installation

```bash
pip install biocompiler
```

Optional extras (you only need these for specific features):

```bash
pip install biocompiler[solver]      # z3 + ortools constraint solving (CSP path)
pip install biocompiler[compare]     # DNAchisel head-to-head comparison harness
pip install biocompiler[viennarna]   # mRNA secondary structure with ViennaRNA
pip install biocompiler[mhcflurry]   # MHC binding predictions for deimmunization
pip install biocompiler[api]         # FastAPI REST API server
pip install biocompiler[all]         # everything (large install)
```

Verify the install:

```bash
biocompiler --version
# BioCompiler v0.9.3
```


## Quick Start: Optimize a Gene

The fastest path from a protein to a synthesis-ready DNA sequence:

```python
from biocompiler import optimize_sequence

# Optimize HBB (human hemoglobin beta, N-terminal 31 residues) for E. coli
result = optimize_sequence(
    target_protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGR",
    organism="e_coli",
    strict_mode=False,
)

print(f"Optimized DNA: {result.optimized_sequence}")
print(f"CAI:           {result.cai:.3f}")
print(f"GC content:    {result.gc_content:.1%}")
print(f"IR verified:   {result.ir_verified}")   # True = protein preserved
print(f"Timed out:     {result.timed_out}")

# Standard file-format exports (populated on every successful run)
print(result.genbank)   # GenBank flat file
print(result.fasta)     # FASTA (protein)
print(result.sbol3)     # SBOL3 (synthetic biology standard)
```

Example output (abridged):

```
Optimized DNA: ATGGTGCACCTGACCCCGGAAGAAAAAAGCGCGGTGACCGCGCTGTGGGGCAAAGTGAACGTGGACGAAGTGGGCGGCGAAGCGCTGGGCCGC
CAI:           1.000
GC content:    65.6%
IR verified:   True
Timed out:     False
```

Three things to notice:

1. **The parameter is `target_protein`**, not `protein`. Passing `protein=`
   raises `InvalidProteinError`.
2. **`ir_verified=True`** means BioCompiler ran the optimized DNA back
   through its L0->L1->L2->L3 IR pipeline and confirmed the translation
   matches your input protein. This is on by default — you do not need to
   ask for it.
3. **`organism="e_coli"`** is a recognized alias. Both canonical names
   (`Escherichia_coli`, `Homo_sapiens`) and short keys (`e_coli`, `human`)
   work everywhere.


## Use Case 1: Optimize a Gene for Expression

The optimizer takes a target protein and a host organism, then produces a
DNA sequence that:

- Encodes the exact input protein (verified by re-translation).
- Maximizes the Codon Adaptation Index (CAI) for the host's codon usage.
- Keeps GC content inside the `[gc_lo, gc_hi]` range (default 30%-70%).
- Avoids internal restriction sites for a default panel of cloning enzymes.
- For eukaryotic hosts: avoids cryptic splice sites, CpG islands, and
  high-affinity miRNA seed matches.
- Reports any `NoRQCTrigger` issues (internal poly-A, GC barriers,
  consecutive rare codons).

### Optimize for E. coli (prokaryote fast path)

```python
from biocompiler import optimize_sequence

# Insulin B chain (human), to be expressed in E. coli
result = optimize_sequence(
    target_protein="FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    organism="e_coli",
    strict_mode=True,        # all predicates must pass
)

print(f"CAI:        {result.cai:.3f}")          # ~0.84
print(f"GC:         {result.gc_content:.1%}")    # ~54%
print(f"Verified:   {result.ir_verified}")       # True
print(f"Failed:     {result.failed_predicates}") # []
```

The prokaryote fast path skips eukaryote-only constraints (cryptic splice
sites, CpG islands, miRNA seeds), which both speeds up optimization and
avoids spurious constraint failures. CAI is high but rarely exactly 1.0
for non-native sequences — the optimizer trades a small amount of CAI to
satisfy GC, restriction-site, and sliding-window constraints.

### Optimize for a eukaryotic host

```python
from biocompiler import optimize_sequence

# Optimize the same protein for human expression
result = optimize_sequence(
    target_protein="FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    organism="Homo_sapiens",
    strict_mode=False,       # see "Troubleshooting" below for why
)

print(f"CAI:        {result.cai:.3f}")          # ~0.91
print(f"GC:         {result.gc_content:.1%}")    # ~60%
print(f"Verified:   {result.ir_verified}")       # True
```

### Tune the GC window and enzyme panel

```python
result = optimize_sequence(
    target_protein="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYG",
    organism="yeast",
    gc_lo=0.35, gc_hi=0.65,                   # tighter GC band
    gc_window_size=60,                        # sliding window (default 50)
    gc_window_min=0.30, gc_window_max=0.75,   # local GC bounds inside the window
    enzymes=["EcoRI", "BamHI", "XhoI"],       # only these sites are avoided
    strict_mode=False,
)
```

### What the result object gives you

| Field                    | Type             | Meaning                                                  |
|--------------------------|------------------|----------------------------------------------------------|
| `sequence` / `optimized_sequence` | `str`    | The optimized DNA (uppercase, `ATGC` only).             |
| `protein`                | `str`            | The input protein (echoed back).                         |
| `cai`                    | `float`          | Codon Adaptation Index in `[0, 1]`.                      |
| `gc_content`             | `float`          | Global GC fraction in `[0, 1]`.                          |
| `ir_verified`            | `bool`           | True iff IR L0->L3 translation matches the input protein. |
| `ir_l0`, `ir_l3`         | IR objects       | The IR objects used for verification.                    |
| `genbank` / `fasta` / `sbol3` | `str`       | Standard file-format exports.                            |
| `failed_predicates`      | `list[str]`      | Predicate names that returned FAIL/LIKELY_FAIL.          |
| `predicate_results`      | `list[PredicateResult]` | Full per-predicate detail.                       |
| `warnings`               | `list[str]`      | Non-fatal notices (cap exceeded, convergence, ...).      |
| `timed_out`              | `bool`           | True if the optimizer hit the wall-clock budget.         |
| `biosecurity_screening_result` | object    | Hazard-screening result (run before optimization).       |


## Use Case 2: Verify a Gene Design

If you already have a DNA sequence (designed by hand, by another tool, or
received from a collaborator), BioCompiler can verify it three ways:

1. **Compile it through the IR pipeline** — confirms the DNA translates to
   the protein you expect.
2. **Run the 43-predicate type system** — checks GC, restriction sites,
   cryptic splice, CpG, codon optimality, RQC triggers, and more.
3. **Generate a certificate** — a machine-checkable record of which
   predicates passed.

### Compile through the IR pipeline

```python
from biocompiler.ir.types import IR_L0_GenomicDNA, IRLevel
from biocompiler.ir.passes import compile_gene

# Your designed DNA sequence (HBB N-terminus, 96 nt = 32 codons)
dna = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA"

# Build IR-L0 (GenomicDNA), then lower to IR-L3 (Polypeptide)
ir_l0 = IR_L0_GenomicDNA(sequence=dna, regions=[], organism="human")
ir_l3 = compile_gene(ir_l0, IRLevel.L3)

print(f"Translated protein: {ir_l3.sequence}")
# -> MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*

# If this matches the protein you expected, your design is structurally sound.
# The trailing "*" is the stop codon.
```

### Run the type system

```python
from biocompiler.type_system import evaluate_all_predicates

results = evaluate_all_predicates(dna, organism="human")

print(f"{'Predicate':<45} {'Verdict':<10}")
print("-" * 55)
for r in results:
    print(f"{r.predicate:<45} {r.verdict.value:<10}")
```

Sample output (truncated):

```
Predicate                                     Verdict
-------------------------------------------------------
NoCrypticSplice                               UNCERTAIN
SpliceCorrect                                 PASS
GCInRange(0.3, 0.7)                          PASS
CodonAdapted(human, 0.5)                      PASS
NoRestrictionSite                             PASS
InFrame                                       PASS
NoInstabilityMotif                            PASS
NoCpGIsland                                   PASS
NoCrypticPromoter(eukaryote, 0.7)             UNCERTAIN
NoUnexpectedTMDomain(True, 0.68)              UNCERTAIN
mRNASecondaryStructure(0, 50, -15.0)          UNCERTAIN
CoTranslationalFolding(human, 0.3)            UNCERTAIN
```

The five-valued verdict logic (`PASS / LIKELY_PASS / UNCERTAIN / LIKELY_FAIL
/ FAIL`) matters here:

- **PASS / FAIL** are *definite* — the predicate is formally verified
  (Lean4) for the 17 core predicates.
- **UNCERTAIN** means the predicate depends on an external tool (ESMFold,
  ViennaRNA, NetMHCpan, BLAST). In the default `CONSERVATIVE` mode, SLOT
  predicates always return UNCERTAIN rather than fabricating a PASS.

To check a specific subset of predicates by name:

```bash
biocompiler check --list-predicates     # list all 43 names
biocompiler check my_gene.fasta --species human --predicate NoStopCodons,NoCpGIsland
```


## Use Case 3: Screen for Biosecurity Hazards

Biosecurity screening runs *automatically* inside `optimize_sequence()`
before any optimization work — if your input protein matches a known toxin
or select-agent signature, optimization is blocked. You can also call the
screener directly:

```python
from biocompiler.biosecurity import screen_hazardous_sequence

# A synthetic test sequence containing the ricin A-chain catalytic motif
# "NIRVGLPIIS" (Lord et al., Toxicon 2003).
suspicious = "MAAAANIRVGLPIISAAAAA"   # contains the ricin motif

report = screen_hazardous_sequence(protein=suspicious)

if report.is_hazardous:
    print(f"HAZARD DETECTED — risk level: {report.risk_level}")
    for m in report.matches:
        print(f"  [{m.category}] {m.name} at pos {m.position} "
              f"(confidence={m.confidence:.2f}, type={m.match_type})")
    for rec in report.recommendations:
        print(f"  -> {rec}")
else:
    print("No hazards detected.")
```

Output:

```
HAZARD DETECTED — risk level: critical
  [select_agent] ricin_a_chain_catalytic at pos 5 (confidence=0.95, type=exact)
  [select_agent] ricin_a_chain_catalytic at pos 4 (confidence=0.85, type=fuzzy)
  [select_agent] ricin_a_chain_catalytic at pos 6 (confidence=0.85, type=fuzzy)
```

The screener also accepts a DNA sequence (it will look for nucleotide
signatures like antibiotic-resistance markers):

```python
report = screen_hazardous_sequence(
    protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGR",
    dna="ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA",
)
print(report.is_hazardous)   # False
print(report.risk_level)     # 'none'
```

### What the screener detects

The signature database covers (see `biosecurity/hazard_signatures.py` for
the full list with citations):

- **Select-agent toxins:** ricin A/B chain, abrin, botulinum neurotoxin,
  shiga toxin, diphtheria toxin, tetanus toxin, cholera toxin.
- **Antibiotic resistance markers** (nucleotide signatures from CARD).
- **Viral surface proteins** and **oncogenes** (EGFR, HER2).

Two screening passes run:

1. **Motif matching** (always runs): sliding-window substring + fuzzy
   matching (Hamming distance <= 2, Levenshtein distance <= 1) against
   short peptide / nucleotide signatures. Reverse-complement screening
   is applied to DNA inputs.
2. **BLAST homology search** (optional): runs only if NCBI BLAST+ is
   installed and `BIOCOMPILER_BLAST_DB_PATH` is set. Otherwise skipped
   gracefully — the motif-based results remain valid.

### Honest limitations

Of 8 reference hazard fragments tested, **6 are correctly flagged**
(ricin A, shiga, diphtheria, tetanus, cholera, abrin). The 2 misses
(botulinum, anthrax LF) are signature-DB coverage gaps, not algorithmic
failures — documented in
[`docs/15-Reference.md`](15-Reference.md). Treat a "no hazards" result
as "no signatures matched", not as "this sequence is safe".


## Use Case 4: Compile from a YAML Spec

The YAML frontend is BioCompiler's "source language" — what a `.c` file
is to gcc. You write a gene design in YAML, the frontend parses it into
an IR-L0, and the IR pipeline lowers it to whatever level you need.

### YAML schema

Required fields:

```yaml
sequence : str    # DNA sequence (A, C, G, T, N; case-insensitive)
organism : str    # source organism, e.g. "human", "e_coli", "Homo_sapiens"
```

Optional fields:

```yaml
gene_name : str    # gene symbol, e.g. "HBB"
regions   : list   # each entry has type / start / end / metadata
metadata  : dict   # free-form top-level metadata
```

Region types recognized by lowering passes: `exon`, `intron`, `5_utr`,
`3_utr`, `promoter`, `terminator`, `cds`. Coordinates are 0-indexed
half-open `[start, end)` — same convention as Python slicing.

### Example YAML (`my_gene.yaml`)

```yaml
gene_name: MyProtein
organism: human
sequence: ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAA
regions:
  - type: cds
    start: 0
    end: 96
    metadata:
      frame: 0
      description: "CDS spanning the full 96 nt fragment"
```

### Compile it

```python
from biocompiler.ir.frontend import compile_from_spec
from biocompiler.ir.types import IRLevel

# From a file path
ir_l3 = compile_from_spec("my_gene.yaml", target_level=IRLevel.L3)
print(f"Protein: {ir_l3.sequence}")
# -> MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*

# From an inline YAML string
ir_l3 = compile_from_spec(
    "gene_name: test\norganism: e_coli\nsequence: ATGGCTAAATGGCGTTAA\n",
    target_level=IRLevel.L3,
)
print(f"Protein: {ir_l3.sequence}")
# -> MAKWR*
```

### Stop at different IR levels

```python
from biocompiler.ir.frontend import compile_from_spec
from biocompiler.ir.types import IRLevel
from biocompiler.ir.codegen import to_genbank, to_fasta_protein, to_sbol3

# Stop at IR-L0 — get the GenBank / SBOL3 of the input DNA
ir_l0 = compile_from_spec("my_gene.yaml", target_level=IRLevel.L0)
print(to_genbank(ir_l0))   # GenBank flat file
print(to_sbol3(ir_l0))     # SBOL3 RDF/Turtle

# Stop at IR-L3 — get the translated protein
ir_l3 = compile_from_spec("my_gene.yaml", target_level=IRLevel.L3)
print(to_fasta_protein(ir_l3))
```

### Use the bundled example spec

BioCompiler ships with an HBB example at
`src/biocompiler/ir/example_specs/hbb.yaml` — useful for quick smoke
tests:

```python
from biocompiler.ir.frontend import compile_from_spec
from biocompiler.ir.types import IRLevel

ir_l3 = compile_from_spec(
    "src/biocompiler/ir/example_specs/hbb.yaml",
    target_level=IRLevel.L3,
)
assert ir_l3.sequence == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR*"
```


## CLI Usage

BioCompiler ships with a `biocompiler` command-line tool. Subcommands:

| Command             | What it does                                                |
|---------------------|-------------------------------------------------------------|
| `optimize PROTEIN`  | Optimize a protein sequence for a target organism.         |
| `batch FILE`        | Optimize multiple proteins from a text file (one per line).|
| `check FILE`        | Check a FASTA gene sequence against all 43 predicates.     |
| `scan`              | Scan a DNA sequence for splice sites and restriction sites.|
| `benchmark`         | Run built-in benchmarks (eGFP, mCherry, LacZ, named sets). |
| `structure`         | Predict protein 3D structure (ESMFold) and assess quality.|
| `stability`         | Analyze protein thermodynamic stability.                   |
| `solubility`        | Analyze protein solubility (CamSol method).                |
| `immunogenicity`    | Analyze / reduce immunogenicity (MHC binding prediction).  |
| `assess`            | Full assessment: structure + stability + solubility + imm.|
| `serve`             | Start the FastAPI REST API server (port 8000 by default).  |
| `explain FILE`      | Explain why a specific codon was chosen (needs provenance).|
| `report FILE`       | Generate a provenance report from a saved trail.           |
| `validate-cai`      | Validate CAI computation against published ground truth.   |
| `validate-maxentscan`| Validate MaxEntScan scores against published values.       |
| `whatif`            | Run what-if analysis on a protein sequence.                |

### Optimize a gene

```bash
# Positional protein argument (v1.0+)
biocompiler optimize "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR" --organism e_coli

# With options
biocompiler optimize "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR" \
    --organism Homo_sapiens \
    --gc-lo 0.35 --gc-hi 0.65 \
    --output hbb_optimized.fasta \
    --certificate hbb_certificate.txt \
    --provenance \
    --verbose

# JSON output for piping into other tools
biocompiler optimize "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR" --organism e_coli --json
```

### Verify a sequence

```bash
# Check a FASTA file against all predicates
biocompiler check my_gene.fasta --species human

# List all available predicate names
biocompiler check --list-predicates

# Check only specific predicates
biocompiler check my_gene.fasta --species human \
    --predicate NoStopCodons,NoCpGIsland,NoRestrictionSite
```

### Scan a sequence for features

```bash
biocompiler scan --sequence ATGGTGCATCTGACTCCTGAGGAG... --enzymes EcoRI,BamHI
```

### Batch optimize

```bash
# proteins.txt: one protein per line, optionally prefixed with a name
# HBB MVHLTPEEK...
# INS FVNQHLCGSHLVEALYLVCGERGFFYTPKT
biocompiler batch proteins.txt --organism e_coli --json --output batch_results.json
```

### Start the REST API

```bash
biocompiler serve --host 0.0.0.0 --port 8000
# Interactive API docs at http://localhost:8000/docs
```


## Supported Organisms

BioCompiler ships with **30 canonical organisms** (60 entries including
short-key aliases). The most common:

| Organism                     | Canonical name              | Short key   | Domain      |
|------------------------------|------------------------------|-------------|-------------|
| *E. coli*                    | `Escherichia_coli`          | `e_coli`    | Prokaryote  |
| *B. subtilis*                | `Bacillus_subtilis`         | `bacillus`  | Prokaryote  |
| *P. putida*                  | `Pseudomonas_putida`        | `pseudomonas`| Prokaryote |
| *C. glutamicum*              | `Corynebacterium_glutamicum`| `corynebacterium`| Prokaryote |
| *H. sapiens*                 | `Homo_sapiens`              | `human`     | Eukaryote   |
| *M. musculus*                | `Mus_musculus`              | `mouse`     | Eukaryote   |
| *R. norvegicus*              | `Rattus_norvegicus`         | `rattus`    | Eukaryote   |
| *S. cerevisiae*              | `Saccharomyces_cerevisiae`  | `yeast`     | Eukaryote   |
| *K. phaffii* (Pichia)        | `Komagataella_phaffii`      | `pichia`    | Eukaryote   |
| *K. lactis*                  | `Kluyveromyces_lactis`      | `kluyveromyces`| Eukaryote |
| *D. rerio* (zebrafish)       | `Danio_rerio`               | `danio`     | Eukaryote   |
| *D. melanogaster*            | `D_melanogaster`            | `drosophila`| Eukaryote   |
| *C. elegans*                 | `Caenorhabditis_elegans`    | `caenorhabditis`| Eukaryote |
| *X. laevis*                  | `Xenopus_laevis`            | `xenopus`   | Eukaryote   |
| *A. thaliana*                | `Arabidopsis_thaliana`      | `arabidopsis`| Eukaryote  |
| *N. benthamiana*             | `Nicotiana_benthamiana`     | `nicotiana` | Eukaryote   |
| *Z. mays* (maize)            | `Zea_mays`                  | `zea`       | Eukaryote   |
| *G. max* (soybean)           | `Glycine_max`               | `glycine`   | Eukaryote   |
| *O. sativa* (rice)           | `Oryza_sativa`              | `oryza`     | Eukaryote   |
| *G. hirsutum* (cotton)       | `Gossypium_hirsutum`        | `gossypium` | Eukaryote   |
| *B. taurus* (cattle)         | `Bos_taurus`                | `bos`       | Eukaryote   |
| *C. familiaris* (dog)        | `Canis_familiaris`          | `canis`     | Eukaryote   |
| *G. gallus* (chicken)        | `Gallus_gallus`             | `gallus`    | Eukaryote   |
| CHO-K1                       | `CHO_K1`                    | `cho`       | Eukaryote   |
| HEK-293T                     | `HEK293T`                   | `hek293`    | Eukaryote   |
| NS0                          | `NS0`                       | `ns0`       | Eukaryote   |
| PER.C6                       | `PER_C6`                    | `per_c6`    | Eukaryote   |
| *C. griseus* (wild-type CHO) | `Cricetulus_griseus_wt`     | `cricetulus`| Eukaryote   |
| *S. frugiperda* (sf9)        | `Spodoptera_frugiperda`     | `sf9`       | Eukaryote   |
| *T. ni* (Hi5)                | `Trichoplusia_ni`           | `hi5`       | Eukaryote   |

The domain (prokaryote vs. eukaryote) is auto-detected and selects which
constraints apply — prokaryotes skip splice-site, CpG-island, and
miRNA-binding constraints that don't apply. You can override:

```python
result = optimize_sequence(
    target_protein="MVHLTPEEK...",
    organism="e_coli",
    organism_domain="prokaryote",   # or "eukaryote" or "auto" (default)
)
```

### Load additional organisms from the Kazusa database

```python
from biocompiler.organisms import resolve_or_download_organism
resolve_or_download_organism("Thermus_thermophilus")
```

tRNA Adaptation Index (tAI) is available for 10 organisms: E. coli, human,
mouse, yeast, CHO, C. elegans, Drosophila, Arabidopsis, Pichia, Bacillus.


## Troubleshooting

### `InvalidProteinError`

```
biocompiler.shared.exceptions.InvalidProteinError:
Invalid amino acid(s) in protein: {'*'}.
```

The `target_protein` parameter accepts the 20 standard amino acids only —
**do not include `*` for the stop codon**. BioCompiler adds the stop codon
during optimization. If you have a protein with `*` at the end, strip it
first:

```python
protein = "MVHLTPEEK...*"
result = optimize_sequence(target_protein=protein.rstrip("*"), organism="e_coli")
```

### `NoRQCTrigger` failure (strict mode)

```
biocompiler.shared.exceptions.OptimizationConstraintError:
Optimization failed 1 predicate(s) in strict mode: NoRQCTrigger.
Set strict_mode=False to allow partial results.
```

`NoRQCTrigger` flags sequences that can trigger ribosome quality control:
internal poly-A stretches, GC-barrier windows, consecutive rare codons, or
a missing terminal stop codon. For short input proteins (< 30 codons) the
optimizer cannot always place a terminal stop codon that satisfies all
constraints, so this predicate can fail in `strict_mode=True`.

**Fix:** use `strict_mode=False` to accept a partial result, then inspect
`result.failed_predicates`:

```python
result = optimize_sequence(
    target_protein="MVHLTPEEKSAVTALWGKVNVDEVGGEALGR",
    organism="Homo_sapiens",
    strict_mode=False,
)
print(result.failed_predicates)   # may include 'NoRQCTrigger'
print(result.warnings)            # human-readable explanations
```

For production pipelines that must pass strict mode, supply a longer
protein (>= 50 residues) so the optimizer has room to satisfy all
constraints.

### Timeout

The optimizer has a 30-second default timeout. For long proteins (> 1000
residues) or tight constraint sets, raise it:

```python
result = optimize_sequence(
    target_protein=long_protein,           # ~1000 aa
    organism="human",
    timeout_seconds=300,                   # 5 minutes
)
if result.timed_out:
    print("Optimizer hit the wall-clock limit; using partial result.")
    print(f"CAI at timeout: {result.cai:.3f}")
```

When `timed_out=True`, the returned sequence is the best partial output
obtained before the timeout — not a failure.

### `IR not verified`

If `result.ir_verified` is `False`, the optimizer produced DNA that does
not translate back to the input protein. This should never happen — it is
a bug. Please file an issue with the input protein, organism, and the
full `OptimizationResult`.

### Unknown organism

```
KeyError: 'Thermus_thermophilus'
```

Either use one of the 30 canonical organisms (see the table above), or
load the organism from the Kazusa database first:

```python
from biocompiler.organisms import resolve_or_download_organism
resolve_or_download_organism("Thermus_thermophilus")
```

This requires network access and caches the codon usage table locally for
subsequent runs.

### SLOT predicates always return UNCERTAIN

This is by design. In the default `CONSERVATIVE` mode, SLOT (Subject to
Limited Oracles and Tools) predicates return UNCERTAIN rather than
fabricating a PASS. To escalate a SLOT predicate to PASS, you must (a)
install the relevant external tool (ViennaRNA, ESMFold, NetMHCpan, BLAST+)
and (b) switch to `VERIFIED` mode. See
[`docs/14-SLOT-Proof-Implementation-Gap.md`](14-SLOT-Proof-Implementation-Gap.md)
for the full SLOT list and the 15 open proof obligations.

### Restriction site not eliminated

The default enzyme panel is `["EcoRI", "BamHI", "XhoI", "HindIII",
"NotI"]`. If your cloning vector uses a different enzyme, pass it
explicitly:

```python
result = optimize_sequence(
    target_protein="...",
    organism="e_coli",
    enzymes=["NdeI", "XhoI"],   # your vector's sites
)
```

A restriction site that overlaps an essential codon with no synonymous
alternative will be reported in `result.failed_predicates` (under
`NoRestrictionSite`) — it is not silently skipped.


## Where to Go Next

- **[EXAMPLES.md](EXAMPLES.md)** — five runnable end-to-end scripts
  (HBB, insulin, gene therapy construct, batch optimize, DNAchisel
  comparison).
- **[README.md](../README.md)** — the technical README: IR pipeline,
  formal verification status, performance, "What It Is NOT".
- **[docs/21-IR-Compiler.md](21-IR-Compiler.md)** — the IR pipeline in
  depth (L0->L4, lowering passes, codegen, alternative splicing).
- **[docs/15-Reference.md](15-Reference.md)** — the 43-predicate
  reference, the engine API, the TCB inventory, and the known
  limitations.
- **[docs/cli.rst](cli.rst)** — full CLI reference (Sphinx).
- **[docs/api.rst](api.rst)** — full Python API reference (Sphinx).
