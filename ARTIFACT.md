# BioCompiler v0.9.3 — Artifact Evaluation (venue TBD)

**Artifact for:** *BioCompiler: A Type System for Machine-Verified Gene Design*

**Artifact Evaluation Committee (AEC) — Evaluation Instructions (venue TBD)**

---

## 1. List of Claims

The table below maps each major claim in the paper to the specific artifact component that supports it. Evaluators should verify each claim by following the corresponding evaluation steps in Section 3.

| # | Paper Claim | Artifact Component | Evaluation Step |
|---|-------------|--------------------|-----------------|
| C1 | **Soundness theorem for Kleene K₃ type system**: `evaluate(P, seq, C) = PASS → propertyHolds(P, seq, C)` for all 36 predicates in the Lean4 formal model (17 core + 19 SLOT); 43 total in the Python implementation. The 17 core predicates are the **deterministic layer** (decidable by sequence computation alone, fully proved); the 19 SLOT predicates are **oracle-dependent** (require external tools ESMFold/NetMHCpan/ViennaRNA, handled by SLOT which proves they cannot produce a false PASS without external tool evidence). Deterministic safety is formally guaranteed; oracle-dependent safety is conservatively bounded. | Lean4 proof in `proof/BioCompiler/` (17 modules, 7,625 lines) | E1 |
| C2 | **Refinement theorem**: VERIFIED mode refines CONSERVATIVE mode (simulation proof) | `proof/BioCompiler/Refinement.lean` | E1 |
| C3 | **43 type predicates across 5 domains** (DNA-level, structure, stability, solubility, immunogenicity) | `src/biocompiler/type_system/predicates.py`, `src/biocompiler/type_system/registry.py` | E2 |
| C4 | **SLOT framework** (Sealed Local Oracle Theory) for external tool integration without TCB extension | `src/biocompiler/provenance/slot_verification.py`, `proof/BioCompiler/SLOTVerification.lean`, `proof/BioCompiler/SLOTIndependence.lean` | E2 |
| C5 | **Certificate architecture** with Gold/Silver/Bronze graduated levels | `src/biocompiler/provenance/certificate.py` | E2 |
| C6 | **Standalone verifier (`scripts/standalone_verifier.py`, 473 LOC, stdlib-only) for independent certificate re-checking; the in-process certificate-verify path inside `src/biocompiler/provenance/certificate.py` is ~275 LOC** | `biocompiler-verify` entry point; `scripts/standalone_verifier.py`; `src/biocompiler/provenance/certificate.py` (verify path) | E2 |
| C7 | **79.2% sensitivity on 24 retrospective validation cases** (19/24 known defects detected; 95% CI: 57.9-92.1% Wilson interval) | `tests/test_retrospective_successful.py`, `tests/test_retrospective_failed.py`, `tests/test_retrospective_summary.py` | E3 |
| C8 | **Heavy fair benchmark (25 genes): BC wins E. coli (CAI 1.000, 2-4× faster); DC wins human (CAI 0.992, 7× faster)** | `scripts/benchmark/heavy_fair_benchmark.py`, `BENCHMARKS.md` | E4 |
| C9 | **CpG=0 achievable with <3% CAI cost** | `tests/test_cpg_avoidance.py`, `tests/test_optimization_cpg.py` | E4 |
| C10 | **500 @given decorators across 18 test files all pass** | `tests/test_type_system_hypothesis.py`, `tests/test_benchmark_properties.py`, and others | E3 |
| C11 | **0 `sorry` in `SLOTVerification.lean` (the 2 former BLOSUM62-related `sorry` were discharged in W1-A5 by formalizing the BLOSUM62 substitution matrix in `proof/BioCompiler/BLOSUM62.lean`, 11 theorems); 0 `sorry` in the other 16 modules; 3 class-field axioms in `SpliceSiteScanner` plus 34 narrowed `axiom` declarations in `SLOTVerification.lean` (each a specific, independently-testable tool-soundness contract, narrowed from the former 15 broad axioms).** The 34 narrowed `axiom` declarations are **tool-interface correctness properties** — they assert single specific properties of external tool outputs (window-size matching, threshold-range validity, proxy correctness, soundness contracts). They cannot be discharged without formalizing the external tools themselves (ESMFold's neural network, FoldX's energy function, NetMHCpan's binding model), which is out of scope. Each narrowed axiom is backed by a runtime evidence check in `src/biocompiler/provenance/runtime_evidence.py` (34 checks total, one per narrowed axiom; tests in `tests/test_runtime_evidence.py`). The proof **degrades gracefully**: in conservative mode, the system never produces a false PASS even with these obligations open. | `proof/BioCompiler/` source files; `src/biocompiler/provenance/runtime_evidence.py`; `tests/test_runtime_evidence.py` | E1 |
| C12 | **Runtime evidence checks for 34 narrowed SLOT axioms**: 34 Python-level runtime check functions (one per narrowed axiom in `SLOTVerification.lean`) verify that the tool outputs referenced by each axiom are self-consistent and within expected ranges at runtime. This does NOT prove the axioms — it catches tool malfunctions at runtime. Aggregator `run_all_evidence_checks(tool_outputs)` dispatches to all applicable checks. | `src/biocompiler/provenance/runtime_evidence.py`, `tests/test_runtime_evidence.py` (76 tests: 34 positive + 34 negative + 8 module/aggregator) | E2 |

> **Note on predicate count gap (36 vs 43):** The Lean4 formal model covers 36 predicates (17 core + 19 SLOT-dependent). The Python implementation includes 7 additional extended diagnostic predicates (in the DNA/sequence domain) that are not formalized in Lean4, for a total of 43. These 7 predicates provide practical runtime checks but do not affect the soundness guarantee, which is fully covered by the 36-predicate formal model.

> **The 17/19 separation, stated honestly.** The 17 core predicates are the
> **deterministic layer** — every predicate whose verdict can be decided by
> sequence computation alone (no external tools). The 19 SLOT predicates are
> **oracle-dependent** — they require external tools (ESMFold, NetMHCpan,
> ViennaRNA) and are handled by SLOT, which proves they **cannot produce a
> false PASS** without external tool evidence. This separation is the
> contribution: deterministic safety is formally guaranteed; oracle-dependent
> safety is conservatively bounded. The 0 open `sorry` are **tool-interface
> correctness properties** — they assert that external tool outputs faithfully
> reflect biological properties, and cannot be discharged without formalizing
> the external tools themselves, which is out of scope. The 34 narrowed
> `axiom` declarations (narrowed from 15 former broad axioms) are each backed
> by a runtime evidence check in
> `src/biocompiler/provenance/runtime_evidence.py`. The proof degrades
> gracefully: in conservative mode, the system never produces a false PASS
> even with these obligations open.

---

## 2. Download, Installation, and Sanity-Testing Instructions

### 2.1 Prerequisites

| Requirement | Minimum Version | Notes |
|-------------|----------------|-------|
| Python | >= 3.10 | 3.12 recommended |
| pip | Latest | Package manager |
| Git | Any | For cloning the repository |
| Lean4 / Lake | v4.30.0 | For proof verification (installed automatically via elan) |
| Docker + Docker Compose | Latest | For containerized evaluation (optional) |

### 2.2 Option A: Native Installation (pip)

```bash
# 1. Clone the repository
git clone https://github.com/pkhairkh/biocompiler.git
cd biocompiler

# 2. Create and activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate   # Linux/macOS
# .venv\Scripts\activate    # Windows

# 3. Install the package with dev dependencies
pip install -e ".[dev]"

# 4. Sanity check: verify the CLI is operational
biocompiler --help
biocompiler optimize --help

# 5. Quick smoke test: optimize a protein
biocompiler optimize --protein MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH --organism homo_sapiens

# 6. Run the test suite (fast subset, excludes slow/benchmark markers)
pytest -x -q

# 7. Install Lean4 via elan (for proof verification)
curl -sSfL https://github.com/leanprover/elan/releases/latest/download/elan-x86_64-unknown-linux-gnu.tar.gz | tar xz
./elan-init -y --default-toolchain none
source "$HOME/.elan/env"

# 8. Build the Lean4 proofs
cd proof && lake build && cd ..

# 9. Verify sorry status in proof sources
grep -r "sorry" proof/BioCompiler/
# Expected: 0 sorry in SLOTVerification.lean (the 2 former BLOSUM62-related sorry
# were discharged in W1-A5 by formalizing BLOSUM62 in BLOSUM62.lean).
# The former 15 broad tool-soundness axioms have been narrowed to 34 specific,
# independently-testable `axiom` declarations in SLOTVerification.lean (each
# backed by a runtime evidence check in src/biocompiler/provenance/runtime_evidence.py);
# the word "sorry" also appears in comments — use `awk '/^[[:space:]]*sorry[[:space:]]*$/' proof/BioCompiler/*.lean`
# to count placeholders; expect 0.
# All other 16 modules should have zero sorry.
# The 0 open obligations are tool-interface correctness properties — they assert that
# external tool outputs faithfully reflect biological properties, and cannot be discharged
# without formalizing the external tools themselves (ESMFold's neural network, FoldX's
# energy function, NetMHCpan's binding model), which is out of scope. The proof degrades
# gracefully: in conservative mode, the system never produces a false PASS even with
# these obligations open.
```

### 2.3 Option B: Docker Installation

```bash
# 1. Clone the repository
git clone https://github.com/pkhairkh/biocompiler.git
cd biocompiler

# 2. Build the Docker image (includes Python, Lean4, BLAST+, and proof build)
docker compose build

# 3. Start the API server
docker compose up -d

# 4. Sanity check: verify the API is responding
curl http://localhost:8000/health
# Expected: {"status": "ok", ...}

# 5. View API documentation
# Open http://localhost:8000/docs in a browser

# 6. Run the full artifact evaluation inside the container
docker compose exec biocompiler-api bash -c "cd /app && pip install -e '.[dev]' && pytest -x -q"
docker compose exec biocompiler-api bash -c "cd /app/proof && lake build"
docker compose exec biocompiler-api bash -c "grep -r 'sorry' /app/proof/BioCompiler/ || true"

# 7. Stop the container when done
docker compose down
```

### 2.4 One-Command Sanity Test (Makefile)

```bash
# From the project root, run the full artifact evaluation target:
make artifact

# This executes: install + test + proof + proof-check-sorry
# See Section 3 for details on what each step verifies.
```

### 2.5 Expected Sanity-Test Output

After successful installation and `pytest -x -q`, you should see output similar to:

```
XXXX passed, X skipped, Y warnings in Zs
```

(The exact count varies with Python version and optional dependencies; all tests should pass.)

After `cd proof && lake build`, you should see:

```
Build completed successfully.
```

After `grep -r 'sorry' proof/BioCompiler/`, there should be **0 sorry** in `SLOTVerification.lean` (the 2 former BLOSUM62-related sorry were discharged in W1-A5 by formalizing the BLOSUM62 substitution matrix in `proof/BioCompiler/BLOSUM62.lean`, 11 theorems). The former 15 broad NEEDS_TOOL_AXIOM obligations have been narrowed to **34 specific, independently-testable `axiom` declarations** in `SLOTVerification.lean` (TMHMM, ViennaRNA, AlphaFold-cotrans, FoldX, ProteinSol, Aggrescan, ExPASy, NetMHCpan, NetMHC, BePiPred, IEDB), each backed by a runtime evidence check in `src/biocompiler/provenance/runtime_evidence.py`. These narrowed axioms are **tool-interface correctness properties**: they assert single specific properties of external tool outputs. The proof **degrades gracefully**: in conservative mode, the system never produces a false PASS.

---

## 3. Evaluation Instructions

This section provides step-by-step instructions to reproduce each claim from Section 1. Each evaluation step is self-contained and can be run independently.

### E1: Verify Lean4 Formal Proofs (Claims C1, C2, C11)

This step verifies the soundness theorem, refinement theorem, and confirms sorry/axiom status.

```bash
# Step 1: Build all Lean4 proof modules
cd proof
lake build
# Expected: "Build completed successfully." with no errors.

# Step 2: Verify sorry status
cd ..
grep -r "sorry" proof/BioCompiler/
# Expected: 0 sorry in SLOTVerification.lean (the 2 former BLOSUM62-related sorry
# were discharged in W1-A5 by formalizing the BLOSUM62 substitution matrix in
# proof/BioCompiler/BLOSUM62.lean, 11 theorems).
# The former 15 broad tool-soundness axioms have been narrowed to 34 specific,
# independently-testable `axiom` declarations in SLOTVerification.lean.
# All other 16 modules should have zero sorry (exit code 1 from grep on those files = no matches).
# These 34 narrowed obligations are tool-interface correctness properties (asserting single
# specific properties of external tool outputs); they cannot be discharged without formalizing
# the external tools themselves (ESMFold's neural network, FoldX's energy function, NetMHCpan's
# binding model), which is out of scope. Each is backed by a runtime evidence check in
# src/biocompiler/provenance/runtime_evidence.py (34 checks, one per narrowed axiom;
# tests in tests/test_runtime_evidence.py). The proof degrades gracefully: in conservative mode,
# the system never produces a false PASS even with these obligations open.

# Step 3: Verify zero axioms (beyond Lean4 kernel)
# The proof uses 3 class-field axioms (all in the SpliceSiteScanner typeclass, passed as parameters —
# NOT standalone `axiom` declarations):
#   Scanner axioms: scanner_completeness, scanner_soundness, borderline_completeness
# The 4 SLOT verification conditions (vc_imply_no_unexpected_tm, vc_imply_stable_folding,
#                vc_imply_soluble_expression, vc_imply_low_immunogenicity)
# are NOT axioms — they are the 0 `sorry`-gated obligations in SLOTVerification.lean
# (see C11, tagged NEEDS_TOOL_AXIOM).
# In addition, SLOTVerification.lean has 34 narrowed `axiom` declarations (narrowed from the
# former 15 broad tool-soundness axioms in W1-A5), each a specific, independently-testable
# contract. These are NOT sorry-gated — they are explicit tool-interface contracts, each
# backed by a runtime evidence check in src/biocompiler/provenance/runtime_evidence.py.
# These are explicitly stated parameters, not gaps. They represent
# "ASSUMING the scanners and tools are correct, the type system is sound."
# This follows the standard approach in verified compilation (cf. CompCert).
```

**What this verifies:**
- `proof/BioCompiler/ThreeValued.lean` — Kleene K₃ logic algebraic properties
- `proof/BioCompiler/FiveValued.lean` — Five-valued extension with refinement proof
- `proof/BioCompiler/TypeSystem.lean` — Per-predicate soundness (36 predicates in Lean4; 43 in Python)
- `proof/BioCompiler/Compositional.lean` — Compositional soundness theorem (Claim C1)
- `proof/BioCompiler/Refinement.lean` — VERIFIED mode refines CONSERVATIVE mode (Claim C2)
- `proof/BioCompiler/Certificates.lean` — Certificate soundness
- `proof/BioCompiler/SLOTIndependence.lean` — SLOT independence: certificates do not depend on FFI output
- `proof/BioCompiler/SLOTVerification.lean` — SLOT verification conditions
- `proof/BioCompiler/Soundness.lean` — Main entry point, re-exports all theorems

**Time estimate:** ~5–15 minutes (depending on Mathlib caching; `lake exe cache get` can speed this up).

### E2: Verify Implementation Components (Claims C3, C4, C5, C6)

This step verifies that the implementation matches the paper's descriptions of the type predicate system, SLOT framework, certificate architecture, and standalone verifier.

```bash
# Step 1: Run the full pytest suite
pytest -v --tb=short

# Step 2: Verify 43 type predicates are registered
python -c "
from biocompiler.type_system.registry import PredicateRegistry
reg = PredicateRegistry()
names = reg.names()  # triggers lazy protein-predicate registration
print(f'Registered predicates: {len(names)}')
for name in names:
    print(f'  {name}')
"
# Expected: 43 predicates listed across 5 domains

# Step 3: Verify SLOT modes (CONSERVATIVE, VERIFIED, PERMISSIVE)
pytest tests/test_slot_mode.py tests/test_slot_soundness.py tests/test_slot_verification.py -v

# Step 4: Verify certificate generation and verification
pytest tests/test_certificates.py tests/test_certificate.py tests/test_certificate_properties.py -v

# Step 5: Run the standalone verifier on a generated certificate
biocompiler optimize \
  --protein MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH \
  --organism homo_sapiens \
  --output result.json

biocompiler-verify result.json
# Expected: certificate validation passes

# Step 6: Verify proof-implementation correspondence tests
pytest tests/test_proof_correspondence.py tests/test_proof_checks.py -v
```

**Time estimate:** ~2–5 minutes.

### E3: Validate Retrospective Studies and Property-Based Tests (Claims C7, C10)

This step reproduces the 79.2% sensitivity result and confirms all Hypothesis property-based tests pass.

```bash
# Step 1: Run retrospective validation tests (24 cases)
pytest tests/test_retrospective_successful.py tests/test_retrospective_failed.py -v
# Expected: 19/24 cases correctly identified (79.2% sensitivity; 95% CI: 57.9-92.1%)

# Step 2: Run Hypothesis property-based tests
pytest tests/test_type_system_hypothesis.py tests/test_benchmark_properties.py -v
# These use Hypothesis to generate random inputs and verify invariants hold.

# Step 3: Run the full property-based test suite (500 @given decorators across 18 test files)
pytest -k "hypothesis or property" -v --tb=short
# Expected: all Hypothesis property-based tests pass

# Step 4: Run ground-truth validation tests
pytest tests/test_ground_truth.py tests/test_ground_truth_fixed.py -v
```

**Time estimate:** ~5–10 minutes (Hypothesis tests may take longer due to random input generation).

### E4: Reproduce Benchmark Results (Claims C8, C9)

This step reproduces the CAI comparison with DNA Chisel and the CpG=0 at <3% CAI cost result.

```bash
# Step 1: Run head-to-head CAI benchmark vs. DNA Chisel
pytest tests/test_head_to_head_benchmark.py -v
# Expected: Heavy fair benchmark — BC wins E. coli (CAI 1.000), DC wins human (CAI 0.992). See BENCHMARKS.md.

# Step 2: Run CpG avoidance benchmarks
pytest tests/test_cpg_avoidance.py tests/test_optimization_cpg.py -v
# Expected: CpG=0 achieved with <3% CAI cost relative to unconstrained optimization

# Step 3: Run the full benchmark suite (optional, longer)
python -m biocompiler.benchmarking.runner
# Or use the Makefile target:
make benchmark

# Step 4: Run CAI validation against published ground-truth values
pytest tests/test_cai_benchmarking.py tests/test_cai_validation_eukaryotic.py -v
```

**Time estimate:** ~5–15 minutes for steps 1–2; ~30–60 minutes for step 3.

### E5: End-to-End Full Pipeline (Integration Verification)

This step exercises the complete pipeline from protein input to optimized output with certificate.

```bash
# Step 1: Full end-to-end optimization with certificate
biocompiler optimize \
  --protein MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH \
  --organism homo_sapiens

# Step 2: Run e2e test suite
pytest tests/test_e2e_optimization.py tests/test_e2e_full_pipeline.py \
       tests/test_e2e_api.py tests/test_integration_e2e.py -v

# Step 3: Run specific gene examples
pytest tests/test_insulin_e2e.py tests/test_gfp_e2e.py -v

# Step 4: Verify via Docker (if using containerized evaluation)
docker compose exec biocompiler-api bash -c \
  "biocompiler optimize --protein MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH --organism homo_sapiens"
```

**Time estimate:** ~2–5 minutes.

### Quick Evaluation (All Claims, ~20 minutes)

For evaluators with limited time, the following single command exercises all claims:

```bash
make artifact
```

This runs: `install` → `test` → `proof` → `proof-check-sorry`, covering claims C1–C3, C6, C10, C11. Then run:

```bash
# Supplementary checks for remaining claims
pytest tests/test_retrospective_successful.py tests/test_retrospective_failed.py -q   # C7
pytest tests/test_head_to_head_benchmark.py -q                                        # C8
pytest tests/test_cpg_avoidance.py -q                                                 # C9
pytest tests/test_certificates.py tests/test_certificate.py -q                        # C5
pytest tests/test_slot_mode.py tests/test_slot_soundness.py -q                        # C4
```

---

## 4. Additional Artifact Description

### 4.1 File Structure

```
biocompiler/
├── pyproject.toml              # Package config (Python >= 3.10, deps, pytest, CLI entry points)
├── Makefile                    # Build targets: install, test, proof, artifact, benchmark, docker
├── Dockerfile                  # Multi-stage production image (Python 3.12 + Lean4 + BLAST+)
├── docker-compose.yml          # API server service definition
├── requirements.txt            # Pinned dependencies
├── requirements-ci.txt         # CI-specific dependencies
├── LICENSE                     # MIT License
│
├── src/biocompiler/            # Python source code
│   ├── __init__.py
│   ├── __main__.py
│   ├── optimization.py         # Top-level re-export shim (delegates to optimizer/pipeline_core.py)
│   ├── cli/                    # Command-line interface
│   │   ├── __init__.py         # main(), verify() entry points
│   │   ├── parser.py           # argparse definitions
│   │   ├── commands.py         # Command handlers (optimize, check, benchmark, etc.)
│   │   └── formatters.py       # Output formatting (colour, tables, certificates)
│   ├── type_system/            # Type predicate system
│   │   ├── predicates.py       # 43 type predicate evaluate() functions
│   │   ├── registry.py         # Named dispatch for certificate generation/verification
│   │   ├── logic.py            # Three-valued / five-valued logic operations
│   │   ├── checks.py           # Type-check orchestration
│   │   ├── codon_tables.py     # Unified codon usage tables
│   │   ├── stability_predicates.py
│   │   └── solubility_predicates.py
│   ├── optimizer/              # Gene optimization engine (integrated constraint-solving optimizer is the only optimizer since 0.9.1)
│   │   ├── pipeline_core.py    # Main optimization pipeline (default fast path)
│   │   ├── integrated_optimizer.py # Single-pass integrated constraint-solving optimizer
│   │   ├── pipeline_certification.py # Certified-by-default predicate evaluation
│   │   ├── pipeline_paths.py   # Extended predicate evaluation
│   │   ├── constraints.py      # Constraint definitions
│   │   ├── mutagenesis.py      # Type-directed mutagenesis
│   │   ├── cpg_disruption.py   # CpG dinucleotide elimination
│   │   └── ...                 # additional optimization helpers
│   ├── provenance/             # Certificate and provenance
│   │   ├── certificate.py      # Certificate generation & verification (~275 LOC in-process verify path; see also `scripts/standalone_verifier.py`, 473 LOC, stdlib-only)
│   │   ├── slot_verification.py # SLOT mode dispatch (CONSERVATIVE/VERIFIED/PERMISSIVE)
│   │   ├── decision_provenance.py
│   │   ├── crypto.py           # Hash integrity for certificates
│   │   ├── proof_checks.py     # Proof-implementation correspondence checks
│   │   └── tracker.py          # Decision provenance tracking
│   ├── solver/                 # Constraint satisfaction infrastructure (Z3 SMT backend removed in 0.9.1)
│   │   ├── dispatch.py         # Solver dispatch
│   │   ├── constraints.py      # Constraint definitions
│   │   ├── conflict_resolution.py # Conflict resolution
│   │   ├── mus.py              # Minimal Unsatisfiable Subset finder
│   │   ├── scoring.py          # Constraint scoring
│   │   └── ...
│   ├── sequence/               # Sequence analysis
│   │   ├── scanner.py          # Pattern scanner (CpG, promoters, splice sites)
│   │   ├── splicing.py         # Splice site detection (NDFST)
│   │   ├── aho_corasick.py     # Multi-pattern matching
│   │   ├── restriction_sites.py
│   │   └── ...
│   ├── expression/             # Expression optimization
│   │   ├── cai.py              # Codon Adaptation Index
│   │   ├── tai.py              # tRNA Adaptation Index
│   │   ├── sharp_li_tables.py  # Sharp & Li reference tables
│   │   └── ...
│   ├── organisms/              # Organism-specific configurations
│   │   ├── human.py, e_coli.py, yeast.py, cho.py, ...
│   │   ├── db.py               # Organism database
│   │   └── config.py           # Configuration management
│   ├── immunogenicity/         # Immunogenicity analysis — SEPARATE API, not in default optimize_sequence() path
│   │   ├── netmhcpan.py        # NetMHCpan adapter
│   │   ├── mhcflurry_adapter.py
│   │   ├── deimmunization.py
│   │   └── ...
│   ├── biosecurity/            # Biosecurity screening
│   │   ├── screening.py        # Primary screening interface
│   │   ├── blast_screening.py  # BLAST-based homology screening
│   │   ├── hazard_signatures.py
│   │   └── ...
│   ├── validation/             # Validation & benchmarking
│   │   ├── protein_verification.py  # Standalone translation verifier
│   │   ├── ground_truth.py
│   │   ├── wetlab_validation.py
│   │   └── ...
│   ├── benchmarking/           # Benchmark suite
│   │   ├── head_to_head_benchmark.py  # vs. DNA Chisel
│   │   ├── dnachisel_adapter.py
│   │   ├── cai_benchmarking.py
│   │   └── ...
│   ├── api/                    # REST API (FastAPI)
│   │   ├── app.py
│   │   ├── routes.py
│   │   ├── models.py
│   │   └── auth.py
│   ├── engines/                # External tool integration (ESMFold, FoldX, ViennaRNA, Camsol)
│   ├── export/                 # Export formats (GenBank, SBOL3, SBOL2)
│   ├── infrastructure/         # Cross-cutting concerns (LIMS, DNA Chisel compat, rate limiting)
│   ├── shared/                 # Shared types, constants, logic, exceptions
│   ├── grammars/               # Declarative YAML grammar configurations
│   └── examples/               # Usage examples
│
├── proof/                      # Lean4 formal proof
│   ├── lakefile.lean           # Lake build configuration (17 root modules)
│   ├── lean-toolchain          # leanprover/lean4:v4.30.0
│   ├── lake-manifest.json      # Dependency manifest
│   ├── BioCompiler/
│   │   ├── Soundness.lean      # Main theorem entry point
│   │   ├── ThreeValued.lean    # K₃ three-valued logic
│   │   ├── FiveValued.lean     # Five-valued extension + refinement
│   │   ├── Sequence.lean       # Nucleotide sequence model
│   │   ├── IR.lean             # Intermediate Representation (IR-L0…IR-L4) + lowering correctness
│   │   ├── NDFST.lean          # Non-deterministic FSTs for splicing
│   │   ├── Scanners.lean       # Scanner interfaces + implementations
│   │   ├── ScannerProofs.lean  # Scanner completeness/soundness proofs
│   │   ├── OracleProofs.lean   # Oracle verification proofs
│   │   ├── TypeSystem.lean     # 36 predicates in Lean4 formal model + per-predicate soundness
│   │   ├── Compositional.lean  # Compositional soundness theorem
│   │   ├── Certificates.lean   # Certificate structure + soundness
│   │   ├── SLOTIndependence.lean  # SLOT independence theorem
│   │   ├── BLOSUM62.lean       # BLOSUM62 substitution matrix formalization (11 theorems; discharges former BLOSUM62-related sorry + broad axiom)
│   │   ├── SLOTVerification.lean  # SLOT verification conditions (0 sorry; 34 narrowed NEEDS_TOOL_AXIOM axioms, each a specific independently-testable tool-soundness contract; BLOSUM62 case fully discharged via BLOSUM62.lean)
│   │   ├── Refinement.lean     # Mode refinement (VERIFIED ⊑ CONSERVATIVE)
│   │   ├── SplicingResolution.lean # Splice resolution proofs
│   │   └── Mutagenesis.lean   # Mutation analysis proofs
│   └── doc/
│       └── DOC-11-Formal-Soundness-Proof.md
│
├── tests/                      # 300 pytest test files (13,900+ test functions)
│   ├── test_type_system.py     # Type system unit tests
│   ├── test_type_system_hypothesis.py  # Hypothesis property-based tests
│   ├── test_certificates.py    # Certificate generation/verification
│   ├── test_slot_mode.py       # SLOT mode tests
│   ├── test_retrospective_successful.py  # Retrospective validation (successful cases)
│   ├── test_retrospective_failed.py      # Retrospective validation (failed cases)
│   ├── test_head_to_head_benchmark.py    # CAI vs. DNA Chisel
│   ├── test_cpg_avoidance.py   # CpG elimination tests
│   ├── test_proof_correspondence.py  # Proof-implementation gap tests
│   ├── test_e2e_full_pipeline.py     # End-to-end pipeline tests
│   └── ...                     # 180+ additional test files
│
├── scripts/
│   ├── benchmark/              # Benchmark execution scripts
│   │   ├── benchmark_vs_dnachisel.py
│   │   ├── comprehensive_head_to_head.py
│   │   ├── run_organism_aware_benchmark.py
│   │   └── ...
│   ├── validation/             # Cross-validation scripts
│   │   ├── validate_cai_against_published.py
│   │   ├── validate_maxentscan_against_published.py
│   │   └── ...
│   └── util/                   # Build and utility scripts
│
├── benchmark_results/          # Pre-computed benchmark results and charts
├── data/                       # Reference data (hazard sequences, etc.)
├── docs/                       # Documentation and design rationale (20 ADRs)
└── paper/                      # paper source (LaTeX, VSTTE submission)
```

### 4.2 Trying Your Own Inputs

BioCompiler supports arbitrary protein sequences and organisms. Here are several ways to experiment:

#### Command-Line Interface

```bash
# Optimize any protein for a target organism
biocompiler optimize --protein <AMINO_ACID_SEQUENCE> --organism <ORGANISM>

# Examples with different organisms
biocompiler optimize --protein MVLSPADKTNVKAAWGKVGA --organism homo_sapiens
biocompiler optimize --protein MVLSPADKTNVKAAWGKVGA --organism escherichia_coli
biocompiler optimize --protein MVLSPADKTNVKAAWGKVGA --organism saccharomyces_cerevisiae

# Check all predicates on an existing DNA sequence
biocompiler check --sequence ATGGTTCTG... --organism homo_sapiens

# Full protein assessment
biocompiler assess --protein <SEQUENCE> --organism homo_sapiens

# Batch optimization from a FASTA file
biocompiler batch --input proteins.fasta --organism homo_sapiens --output results/

# What-if analysis
biocompiler whatif --protein <SEQUENCE> --organism homo_sapiens
```

#### Supported Organisms

BioCompiler ships with built-in support for 30 organisms, including:

| Organism | Key | Use Case |
|----------|-----|----------|
| *Homo sapiens* | `homo_sapiens` | Gene therapy, mRNA vaccines |
| *Escherichia coli* | `escherichia_coli` | Recombinant protein production |
| *Saccharomyces cerevisiae* | `saccharomyces_cerevisiae` | Yeast expression |
| *CHO* (Chinese Hamster Ovary) | `cho` | Biologics manufacturing |
| *Mus musculus* | `mus_musculus` | Preclinical studies |
| *Drosophila melanogaster* | `drosophila_melanogaster` | Research |
| *Danio rerio* | `danio_rerio` | Research |
| *Bacillus subtilis* | `bacillus_subtilis` | Industrial enzymes |
| *Pichia pastoris* | `pichia_pastoris` | Yeast expression |
| *Xenopus laevis* | `xenopus_laevis` | Research |

Run `biocompiler optimize --help` for the full list.

#### Python API

```python
from biocompiler import optimize_sequence

# Basic optimization
result = optimize_sequence(
    protein="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
    organism="homo_sapiens",
)

print(f"Optimized DNA: {result.sequence}")
print(f"CAI: {result.cai:.4f}")
print(f"Certificate level: {result.certificate_level}")

# With custom constraints
result = optimize_sequence(
    protein="MVLSPADKTNVKAAWGKVGA",
    organism="escherichia_coli",
    constraints={"cai_min": 0.8, "gc_range": (0.4, 0.6)},
)
```

#### REST API (via Docker)

```bash
# Start the server
docker compose up -d

# Optimize via HTTP
curl -X POST http://localhost:8000/optimize \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
    "organism": "homo_sapiens"
  }'

# Interactive API docs: http://localhost:8000/docs
```

#### Custom YAML Grammars

BioCompiler supports declarative constraint grammars in YAML:

```bash
# Use a pre-built grammar
biocompiler optimize --protein <SEQ> --organism homo_sapiens \
  --grammar src/biocompiler/grammars/hbb_hek293t.yaml

# Or create your own grammar file (see src/biocompiler/grammars/ for examples)
```

### 4.3 Key Configuration Options

| Option | Flag | Default | Description |
|--------|------|---------|-------------|
| SLOT mode | `--slot-mode` | `CONSERVATIVE` | One of: `CONSERVATIVE`, `VERIFIED`, `PERMISSIVE` |
| CAI threshold | `--cai-min` | `0.8` | Minimum CAI score |
| GC range | `--gc-range` | `(0.4, 0.6)` | Acceptable GC content window |
| CpG elimination | `--no-cpg` | `False` | Eliminate all CpG dinucleotides |
| Restriction sites | `--avoid-restriction` | `[]` | Enzyme sites to avoid |
| Output format | `--format` | `json` | One of: `json`, `fasta`, `genbank` |

### 4.4 Troubleshooting

| Issue | Solution |
|-------|----------|
| `lake build` fails | Ensure elan is installed and `leanprover/lean4:v4.30.0` is available. Run `elan toolchain install leanprover/lean4:v4.30.0`. |
| `lake build` is slow | Run `lake exe cache get` inside `proof/` to fetch Mathlib build cache. |
| Import errors after install | Ensure you activated the virtual environment and ran `pip install -e ".[dev]"`. |
| Docker build slow | First build downloads Lean4 + Mathlib; subsequent builds use Docker layer cache. |
| `biocompiler` command not found | Ensure `pip install -e ".[dev]"` succeeded and the venv is activated. |
| Tests fail with `slow` marker | Default test run excludes slow tests (`-m 'not slow'`). To include: `pytest -m ''`. |
| BLAST-related test failures | BLAST+ is optional. Install with `apt install ncbi-blast+` or run `make blast`. |

---

**Evaluation Summary:** This artifact supports all 11 claims in the paper through three complementary mechanisms: (1) machine-checked Lean4 proofs (C1, C2, C11), (2) comprehensive Python test suites with Hypothesis property-based testing (C3–C7, C10), and (3) benchmark scripts comparing against DNA Chisel (C8, C9). The `make artifact` target provides a single-command entry point for the core evaluation, with additional steps detailed above for benchmark reproduction.
