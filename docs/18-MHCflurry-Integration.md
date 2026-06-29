# 18 — MHCflurry Integration

> **Document ID:** 18
> **Date:** 2026-06-04
> **Addresses:** Caveat 3 (Heuristic Engines) — Immunogenicity Upgrade
> **Related:** `biocompiler/immunogenicity/mhcflurry_adapter.py`, `biocompiler/immunogenicity/mhcflurry_population.py`, `12-Engine-Accuracy.md`

---

## 1. Overview

BioCompiler's immunogenicity module previously offered two prediction tiers:

| Tier | Method | AUC-ROC | Offline? |
|------|--------|---------|----------|
| Default | PSSM (position-specific scoring matrices) | 0.60–0.75 | Yes |
| Optional | NetMHCpan 4.1 (remote API) | 0.85–0.95 | No |

The gap between these tiers is enormous — both in accuracy and operational
characteristics.  PSSMs are unreliable (AUC as low as 0.60); NetMHCpan
is accurate but requires a network connection and is subject to rate limits
and downtime.

**MHCflurry fills this gap.** It is an open-source, offline neural-network
predictor that:

- Runs **completely offline** — no API calls, no rate limits, no network
  dependency.
- Achieves AUC-ROC of **0.80–0.85**, a substantial improvement over PSSM.
- Supports **15,000+ MHC-I alleles** (compared to 9 in the PSSM tier).
- Includes **antigen processing prediction** (proteasomal cleavage + TAP
  transport) via the presentation model.
- Downloads models once (~100 MB), cached in `~/.mhcflurry/`.

MHCflurry replaces PSSM as the recommended default for offline
immunogenicity screening, establishing a four-tier prediction hierarchy:

```
NetMHCpan (0.85–0.95)  ←  Best accuracy, requires API
        │
        ▼ fallback
MHCflurry (0.80–0.85)  ←  Best offline, recommended default
        │
        ▼ fallback
Precomputed database   ←  Pre-computed binding scores for common alleles
        │
        ▼ fallback
PSSM (0.60–0.75)       ←  Always available, lowest accuracy
```

---

## 2. Installation

### 2.1 pip install

```bash
pip install mhcflurry
```

MHCflurry requires Python ≥ 3.8 and TensorFlow (or TensorFlow-compatible
backend).  The pip install pulls in all dependencies automatically.

### 2.2 Model download

On first use, MHCflurry needs its trained neural-network models (~100 MB).
Download them explicitly:

```python
from biocompiler.immunogenicity.mhcflurry_adapter import download_models

download_models()  # downloads to ~/.mhcflurry/ by default
```

Or from the command line:

```bash
python -m mhcflurry-downloads fetch
```

The download includes three model sets:

| Model Set | Purpose | Size | Required |
|-----------|---------|------|----------|
| `models_class1` | MHC-I binding affinity (IC50) | ~60 MB | **Yes** |
| `models_class1_presentation` | Antigen processing + binding | ~25 MB | Optional |
| `models_class1_processing` | Proteasomal cleavage + TAP | ~15 MB | Optional |

Models are cached in `~/.mhcflurry/` and only need to be downloaded once.
Subsequent runs load from disk with no network access.

### 2.3 Verification

```python
from biocompiler.immunogenicity.mhcflurry_adapter import is_mhcflurry_available

if is_mhcflurry_available():
    print("MHCflurry is ready for predictions")
else:
    print("MHCflurry not available — install + download models")
```

`is_mhcflurry_available()` returns `True` only when **both** conditions are
met: (1) `import mhcflurry` succeeds, and (2) at least one model set is
present in the models directory.

---

## 3. Prediction Hierarchy

BioCompiler's immunogenicity stack now follows a four-tier fallback chain
for MHC-I binding prediction:

```
predict_t_cell_epitopes()
    │
    ├─ 1. NetMHCpan  (if use_netmhcpan=True AND API reachable)
    │     AUC-ROC: 0.85–0.95
    │     Requires: network, NetMHCpan service
    │
    ├─ 2. MHCflurry  (if mhcflurry available AND enabled)
    │     AUC-ROC: 0.80–0.85
    │     Requires: mhcflurry package + downloaded models
    │
    ├─ 3. Precomputed database  (if available)
    │     Pre-computed binding scores for common allele-peptide pairs
    │     Requires: database file
    │
    └─ 4. PSSM       (always available)
          AUC-ROC: 0.60–0.75
          Requires: nothing
```

**Priority rules:**

1. If NetMHCpan is explicitly requested and available → use NetMHCpan.
2. If MHCflurry is enabled and models are present → use MHCflurry.
3. If precomputed database is available → use precomputed scores.
4. Fall back to PSSM if none of the above are available.

This hierarchy ensures the best available accuracy is always used, with
graceful degradation when external resources are unavailable.

---

## 4. API Reference

### 4.1 Module-level functions

```python
from biocompiler.immunogenicity.mhcflurry_adapter import (
    is_mhcflurry_available,   # Check if MHCflurry is ready
    download_models,          # Download neural-network models
    clear_cache,              # Clear the module-level prediction cache
    MHCFLURRY_AUC_ROC_LOW,   # 0.80
    MHCFLURRY_AUC_ROC_HIGH,  # 0.85
)
```

#### `is_mhcflurry_available() -> bool`

Check whether MHCflurry is importable **and** models are downloaded.
Returns `True` only when both conditions are satisfied.

#### `download_models(models_dir=None, verbose=True) -> bool`

Download MHCflurry models (~100 MB).  If `models_dir` is `None`, uses
`~/.mhcflurry/`.  Returns `True` on success, `False` on failure.
A network connection is required for the download itself.

#### `clear_cache() -> None`

Clear the module-level LRU prediction cache.  This does **not** unload
loaded MHCflurry models — only the prediction result cache is cleared.

### 4.2 `MHCflurryClient`

The main client for offline MHC-I binding and presentation prediction.

```python
from biocompiler.immunogenicity.mhcflurry_adapter import MHCflurryClient

client = MHCflurryClient(models_dir=None)  # uses ~/.mhcflurry/ by default
```

Models are **lazy-loaded** on the first prediction call, keeping
construction lightweight.

#### `predict_binding(peptide, allele) -> MHCBindingResult`

Predict MHC-I binding affinity for a single peptide–allele pair.

```python
result = client.predict_binding("SIINFEKL", "HLA-A*02:01")
print(result.binding_class)  # "strong_binder", "moderate_binder", etc.
print(result.ic50_nm)        # predicted IC50 in nM
print(result.binding_score)  # normalised 0–1 score
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `peptide` | `str` | Amino-acid sequence (8–11 residues typical) |
| `allele` | `str` | MHC-I allele name (e.g. `"HLA-A*02:01"`) |

Returns an `MHCBindingResult` dataclass matching the
`biocompiler.immunogenicity` interface.

#### `batch_predict(protein, alleles, epitope_lengths=None) -> list[MHCBindingResult]`

Scan a full protein for MHC-I binders across multiple alleles.
Extracts all overlapping peptides of the requested lengths and runs
MHCflurry in batch mode (significantly faster than calling
`predict_binding` individually).

```python
results = client.batch_predict(
    protein="MAGRSGDLDAIIRYVKQLRYLENGKETLQRT",
    alleles=["HLA-A*02:01", "HLA-A*03:01", "HLA-B*07:02"],
    epitope_lengths=[8, 9, 10],
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protein` | `str` | — | Full protein amino-acid sequence |
| `alleles` | `list[str]` | — | MHC-I alleles to evaluate |
| `epitope_lengths` | `list[int] \| None` | `[8,9,10,11]` | Peptide lengths to extract |

Unsupported alleles are silently skipped with a debug log.

#### `predict_presentation(protein, alleles, epitope_lengths=None) -> list[MHCBindingResult]`

Predict antigen processing **plus** MHC binding (presentation score).
Uses MHCflurry's `Class1PresentationPredictor` which combines:

1. **MHC binding affinity** (neural-network IC50 prediction)
2. **Proteasomal cleavage** (likely cleavage sites)
3. **TAP transport efficiency** (Transporter associated with Antigen Processing)

This is more accurate than binding prediction alone for identifying
naturally presented epitopes.

```python
results = client.predict_presentation(
    protein="MAGRSGDLDAIIRYVKQLRYLENGKETLQRT",
    alleles=["HLA-A*02:01"],
)
# result.binding_score now reflects processing + binding
```

Falls back to `batch_predict` (binding only) if the presentation model
is not available.



#### `clear_cache() -> None`

Clear the prediction cache for this client instance.

### 4.3 Result format — `MHCBindingResult`

MHCflurry returns IC50 values (nM) and presentation percentiles.  The
adapter converts them to the standard `MHCBindingResult` format used
throughout BioCompiler's immunogenicity module:

```python
@dataclass
class MHCBindingResult:
    allele: str              # e.g. "HLA-A*02:01"
    peptide: str             # e.g. "SIINFEKL"
    start_position: int      # 0-based start in source protein
    end_position: int        # 0-based inclusive end
    binding_score: float     # normalised 0–1 (higher = stronger binding)
    ic50_nm: float           # predicted IC50 in nanomolars
    binding_class: str       # "strong_binder", "moderate_binder", etc.
    anchor_residues: dict    # {position: amino_acid}
    anchor_scores: dict      # {position: score}
```

**Binding score conversion:**

```
binding_score = 1 − log(IC50) / log(50,000)
```

| IC50 (nM) | binding_score | Classification |
|-----------|---------------|----------------|
| 1 | ≈ 1.00 | Strong binder |
| 50 | ≈ 0.75 | Strong binder |
| 500 | ≈ 0.55 | Moderate binder |
| 5,000 | ≈ 0.36 | Weak binder |
| 50,000 | = 0.00 | Non-binder |

**Binding class thresholds** (standard across all BioCompiler
immunogenicity methods):

| IC50 (nM) | Class |
|-----------|-------|
| < 50 | `strong_binder` |
| 50 – 500 | `moderate_binder` |
| 500 – 5,000 | `weak_binder` |
| > 5,000 | `non_binder` |

---

## 5. Allele Coverage

### 5.1 Before MHCflurry: 9 alleles (PSSM)

The original PSSM implementation in `immunogenicity.py` provides matrices
for only:

| Class | Alleles | Count |
|-------|---------|-------|
| MHC-I | `HLA-A*02:01`, `HLA-A*01:01`, `HLA-A*03:01`, `HLA-A*24:02`, `HLA-B*07:02`, `HLA-B*08:01` | 6 |
| MHC-II | `HLA-DRB1*01:01`, `HLA-DRB1*04:01`, `HLA-DRB1*07:01` | 3 |
| **Total** | | **9** |

This covers only the most common HLA alleles in Caucasian populations
and provides very poor representation for African, Asian, and Indigenous
populations.

### 5.2 After MHCflurry: 15,000+ alleles

MHCflurry 2.0 models are trained on IEDB data covering over 15,000
MHC-I alleles.  The `mhcflurry_population.py` module ships with
frequency data for the 50 most common alleles across 6 population groups:

| Locus | Alleles with frequency data | MHCflurry-supported |
|-------|-----------------------------|---------------------|
| HLA-A | 20 | ~60+ |
| HLA-B | 18 | ~90+ |
| HLA-C | 12 | ~25+ |
| Non-human | 6 (mouse H-2) | 6 |

The `SUPPORTED_MHCFLURRY_ALLELES` constant in `mhcflurry_population.py`
lists ~200 of the most commonly used alleles with explicit classification.
The actual MHCflurry model supports far more — check the MHCflurry package
documentation for the full list of supported alleles.

### 5.3 Allele classification

Each allele is classified by both MHC class and data source:

```python
from biocompiler.immunogenicity.mhcflurry_population import ALLELE_CLASSIFICATION

ALLELE_CLASSIFICATION["HLA-A*02:01"]  # "I:both"    (PSSM + MHCflurry)
ALLELE_CLASSIFICATION["HLA-A*68:01"]  # "I:mhcflurry" (MHCflurry only)
ALLELE_CLASSIFICATION["HLA-DRB1*01:01"]  # "II:pssm"  (PSSM only)
```

Format: `"<class>:<source>"` where:
- `<class>`: `I` or `II`
- `<source>`: `pssm`, `mhcflurry`, or `both`

---

## 6. Accuracy Comparison

### 6.1 AUC-ROC benchmarks (MHC-I binding prediction)

| Method | AUC-ROC Range | Offline? | Speed (100 peptides × 5 alleles) |
|--------|--------------|----------|----------------------------------|
| **NetMHCpan 4.1** | **0.85–0.95** | No | ~10 s (network latency dominates) |
| **MHCflurry 2.0** | **0.80–0.85** | **Yes** | ~2 s (local computation) |
| **PSSM** | **0.60–0.75** | Yes | ~0.1 s (simple matrix scan) |

**Key observations:**

1. MHCflurry closes ~60% of the accuracy gap between PSSM and NetMHCpan,
   while remaining fully offline.
2. PSSM captures anchor positions but misses peptide conformation, TCR
   contact residues, and MHC-peptide interaction dynamics.
3. NetMHCpan remains the gold standard for high-stakes predictions, but
   requires network access and is subject to rate limits.

### 6.2 MHC-II accuracy

| Method | AUC-ROC Range | Notes |
|--------|--------------|-------|
| NetMHCIIpan 4.0 | 0.80–0.90 | Gold standard for MHC-II |
| PSSM | 0.55–0.70 | Even lower for MHC-II (open groove) |
| **MHCflurry** | **N/A** | **MHC-II not supported** |

### 6.3 Presentation prediction accuracy

MHCflurry's presentation predictor (binding + processing) improves
positive predictive value for naturally presented epitopes by filtering
out peptides that bind MHC but are unlikely to be generated by the
antigen processing pathway.  Published benchmarks (O'Donnell et al.,
*Nat Biotechnol* 2021) show:

- Presentation percentile ≤ 2.0 identifies ~80% of experimentally
  validated epitopes at ~5% false positive rate.
- Binding-only prediction at IC50 ≤ 50 nM identifies ~90% of epitopes
  but at ~20% false positive rate.

### 6.4 Module accuracy constants

```python
from biocompiler.immunogenicity.mhcflurry_adapter import (
    MHCFLURRY_AUC_ROC_LOW,   # 0.80
    MHCFLURRY_AUC_ROC_HIGH,  # 0.85
)

from biocompiler.immunogenicity import (
    IMMUNOGENICITY_PSSM_AUC_ROC_LOW,       # 0.60
    IMMUNOGENICITY_PSSM_AUC_ROC_HIGH,      # 0.75
    IMMUNOGENICITY_NETMHCPAN_AUC_ROC_LOW,  # 0.85
    IMMUNOGENICITY_NETMHCPAN_AUC_ROC_HIGH, # 0.95
)
```

---

## 7. Population Coverage

### 7.1 Why population coverage matters

Immunogenicity screening is only meaningful if the MHC alleles tested
represent the target patient population.  A therapeutic protein that
appears non-immunogenic when tested against 6 common alleles may still
trigger immune responses in populations carrying rarer alleles.

### 7.2 Coverage improvement

| Allele set | Alleles | Global coverage | Caucasian | African | Asian |
|------------|---------|-----------------|-----------|---------|-------|
| PSSM MHC-I (6) | 6 | ~40% | ~55% | ~20% | ~25% |
| PSSM MHC-I + MHC-II (9) | 9 | ~45% | ~60% | ~25% | ~30% |
| Top 12 MHCflurry (optimized) | 12 | ~75% | ~85% | ~60% | ~70% |
| Top 25 MHCflurry | 25 | ~90% | ~95% | ~85% | ~88% |
| Top 50 MHCflurry | 50 | **>95%** | **>98%** | **>90%** | **>95%** |

### 7.3 Population coverage API

The `mhcflurry_population.py` module provides utilities for computing
and optimising population coverage:

```python
from biocompiler.immunogenicity.mhcflurry_population import (
    compute_population_coverage,
    find_coverage_optimizing_alleles,
    get_allele_frequency,
    EXPANDED_POPULATION_COVERAGE,
    POPULATION_GROUPS,
    POPULATION_WEIGHTS,
    SUPPORTED_MHCFLURRY_ALLELES,
)
```

#### `compute_population_coverage(alleles, population="global") -> float`

Compute the fraction of a population covered by a set of alleles.

```python
# Coverage with the 6 PSSM MHC-I alleles
coverage = compute_population_coverage([
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01",
    "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01",
], population="global")
# ≈ 0.40

# Coverage with 50 optimised alleles
best_50 = find_coverage_optimizing_alleles(n_alleles=50)
coverage = compute_population_coverage(best_50, population="global")
# ≈ 0.95+
```

Coverage formula (Hardy-Weinberg, independence assumption):

```
P(covered) = 1 − ∏_i (1 − freq_i / 100)
```

The independence assumption is an approximation — true HLA alleles are
in linkage disequilibrium — but provides a useful upper-bound estimate
for screening.

#### `find_coverage_optimizing_alleles(n_alleles=6, population="global") -> list[str]`

Greedy selection of alleles that maximise population coverage.

```python
# Find the best 12 alleles for global coverage
best_12 = find_coverage_optimizing_alleles(n_alleles=12, population="global")

# Find the best 6 alleles for Asian populations
best_asian_6 = find_coverage_optimizing_alleles(n_alleles=6, population="Asian")
```

#### `get_allele_frequency(allele, population) -> float`

Look up the phenotype frequency of an allele in a population.

```python
freq = get_allele_frequency("HLA-A*02:01", "Caucasian")  # 28.0 (%)
freq = get_allele_frequency("HLA-A*02:01", "Asian")       # 10.0 (%)
```

### 7.4 Population groups and data sources

Six population groups are supported:

| Group | Weight (millions) | Data source |
|-------|-------------------|-------------|
| Caucasian | 1,200 | AFND, NMDP |
| African | 1,400 | AFND, Middleton & Williams 2000 |
| Asian | 4,700 | AFND, Mori et al. 1997 |
| Hispanic | 650 | AFND, NMDP |
| Native American | 50 | AFND (limited data) |
| Oceanian | 45 | AFND (limited data) |

**All frequency values are approximate estimates** suitable for
computational screening but **not for clinical decision-making**.
Validate against the current AFND release before clinical use.

---

## 8. Configuration

### 8.1 Enabling MHCflurry

MHCflurry is automatically used when available and enabled.  Configure
via the `compute_immunogenicity()` function parameters or the API:

```python
from biocompiler.immunogenicity import compute_immunogenicity

# Default: uses MHCflurry if available, falls back to PSSM
result = compute_immunogenicity(protein, mhc_alleles=my_alleles)

# Force PSSM only (disable MHCflurry)
result = compute_immunogenicity(
    protein,
    mhc_alleles=my_alleles,
    use_mhcflurry=False,
)

# Force NetMHCpan (overrides MHCflurry)
result = compute_immunogenicity(
    protein,
    mhc_alleles=my_alleles,
    use_netmhcpan=True,
)
```

### 8.2 Custom models directory

By default, MHCflurry models are stored in `~/.mhcflurry/`.  To use a
custom location:

```python
from biocompiler.immunogenicity.mhcflurry_adapter import MHCflurryClient, download_models

# Download to custom directory
download_models(models_dir="/data/mhcflurry_models")

# Use custom directory
client = MHCflurryClient(models_dir="/data/mhcflurry_models")
```

Or set the environment variable:

```bash
export MHCFLURRY_DEFAULT_CLASS1_MODELS_DIR=/data/mhcflurry_models/models_class1/models
```

### 8.3 Cache configuration

The LRU cache holds up to 5,000 entries by default.  To clear the cache:

```python
from biocompiler.immunogenicity.mhcflurry_adapter import clear_cache

clear_cache()  # module-level cache

# Or per-client
client = MHCflurryClient()
client.clear_cache()
```

### 8.4 Epitope lengths

Default epitope lengths for MHC-I scanning are `[8, 9, 10, 11]`.
Override per call:

```python
results = client.batch_predict(
    protein=protein,
    alleles=alleles,
    epitope_lengths=[9, 10],  # common MHC-I lengths only
)
```

---

## 9. Integration Points

MHCflurry plugs into BioCompiler at several levels:

### 9.1 `immunogenicity.py` — T-cell epitope prediction

The `predict_t_cell_epitopes()` function is the primary integration
point.  When MHCflurry is available, MHC-I binding predictions use the
neural-network predictor instead of PSSM matrices.

```
predict_t_cell_epitopes()
    └── predict_mhc_i_binding()
            ├── MHCflurryClient.batch_predict()  ← if mhcflurry available
            └── _pssm_scan()                     ← fallback
```

MHC-II predictions continue to use PSSM or NetMHCIIpan (MHCflurry does
not support MHC-II).

### 9.2 `compute_immunogenicity()` — Combined score

The overall immunogenicity computation delegates to
`predict_t_cell_epitopes()`, so MHCflurry results flow through
automatically:

```
compute_immunogenicity()
    ├── predict_t_cell_epitopes()  ← uses MHCflurry for MHC-I
    │       └── MHCflurryClient.batch_predict()
    ├── predict_epitopes()         ← B-cell (unchanged)
    └── find_deimmunization_mutations()  ← uses epitope results
```

### 9.3 `mhcflurry_population.py` — Population coverage

The expanded allele frequency data is used by:

- `compute_immunogenicity()` — population coverage metric
- `PopulationCoverageSafe` predicate (type system predicate #28)
- Coverage-optimising allele selection for screening panels

### 9.4 `api.py` — REST API

The FastAPI endpoints for immunogenicity analysis accept the same
parameters and follow the same prediction hierarchy:

```
POST /api/v1/immunogenicity/analyze
    └── use_mhcflurry: bool = True  (default: auto-detect)

POST /api/v1/immunogenicity/deimmunize
    └── uses MHCflurry predictions for mutation prioritisation
```

### 9.5 `slot_verification.py` — SLOT predicates

MHCflurry predictions feed into the immunogenicity SLOT predicates:

- `LowImmunogenicity` — overall immunogenicity score
- `NoStrongTCellEpitope` — MHC binding epitope detection
- `PopulationCoverageSafe` — MHC allele population coverage

These predicates always return `UNCERTAIN` in the formal model
(SLOT-dependent), but the MHCflurry-backed predictions provide
substantially more accurate SLOT values for downstream review.

---

## 10. Examples

### 10.1 Basic binding prediction

```python
from biocompiler.immunogenicity.mhcflurry_adapter import MHCflurryClient

client = MHCflurryClient()

# Predict binding for a single peptide
result = client.predict_binding("SIINFEKL", "HLA-A*02:01")
print(f"IC50: {result.ic50_nm:.1f} nM")
print(f"Binding class: {result.binding_class}")
print(f"Binding score: {result.binding_score:.3f}")
```

### 10.2 Full protein scan

```python
from biocompiler.immunogenicity.mhcflurry_adapter import MHCflurryClient

client = MHCflurryClient()

protein = "MAGRSGDLDAIIRYVKQLRYLENGKETLQRTDAPKSVPLPKGLSNYDYTR"
alleles = ["HLA-A*02:01", "HLA-A*03:01", "HLA-B*07:02"]

results = client.batch_predict(protein, alleles, epitope_lengths=[8, 9, 10])

# Filter for strong binders
strong = [r for r in results if r.binding_class == "strong_binder"]
print(f"Found {len(strong)} strong binders out of {len(results)} predictions")

for r in sorted(strong, key=lambda x: x.ic50_nm):
    print(f"  {r.allele}: {r.peptide} (IC50={r.ic50_nm:.1f} nM, "
          f"pos={r.start_position}-{r.end_position})")
```

### 10.3 Antigen processing prediction

```python
from biocompiler.immunogenicity.mhcflurry_adapter import MHCflurryClient

client = MHCflurryClient()

protein = "MAGRSGDLDAIIRYVKQLRYLENGKETLQRTDAPKSVPLPKGLSNYDYTR"
alleles = ["HLA-A*02:01"]

# Predict presentation (binding + proteasomal cleavage + TAP transport)
results = client.predict_presentation(protein, alleles)

for r in results[:5]:  # top 5
    print(f"  {r.peptide} — presentation score: {r.binding_score:.4f}, "
          f"IC50: {r.ic50_nm:.1f} nM")
```

### 10.4 Population coverage optimisation

```python
from biocompiler.immunogenicity.mhcflurry_population import (
    find_coverage_optimizing_alleles,
    compute_population_coverage,
)

# Find the 12 alleles that maximise global population coverage
best_12 = find_coverage_optimizing_alleles(n_alleles=12, population="global")
coverage = compute_population_coverage(best_12, population="global")
print(f"Global coverage with 12 alleles: {coverage:.1%}")

# Compare per-population
for pop in ["Caucasian", "African", "Asian", "Hispanic"]:
    cov = compute_population_coverage(best_12, population=pop)
    print(f"  {pop}: {cov:.1%}")
```

### 10.5 Integration with compute_immunogenicity()

```python
from biocompiler.immunogenicity import compute_immunogenicity
from biocompiler.immunogenicity.mhcflurry_population import find_coverage_optimizing_alleles

# Use optimised allele panel for immunogenicity screening
best_alleles = find_coverage_optimizing_alleles(n_alleles=25)

protein = "MAGRSGDLDAIIRYVKQLRYLENGKETLQRTDAPKSVPLPKGLSNYDYTR"
result = compute_immunogenicity(protein, mhc_alleles=best_alleles)

print(f"Overall immunogenicity: {result.immunogenicity_score:.3f}")
print(f"Classification: {result.immunogenicity_class}")
print(f"T-cell epitopes: {len(result.t_cell_epitopes)}")
```

### 10.6 Availability check and graceful fallback

```python
from biocompiler.immunogenicity.mhcflurry_adapter import is_mhcflurry_available, MHCflurryClient
from biocompiler.immunogenicity import compute_immunogenicity

if is_mhcflurry_available():
    client = MHCflurryClient()
    print(f"MHCflurry available — using MHCflurry for predictions")
else:
    # Falls back to PSSM automatically
    alleles = None  # uses default PSSM alleles
    print("MHCflurry not available — using PSSM fallback")

result = compute_immunogenicity(protein, mhc_alleles=alleles)
```

---

## 11. Known Limitations

### 11.1 MHC-II not supported

MHCflurry 2.0 does **not** natively support MHC class II binding
prediction.  MHC-II alleles continue to use PSSM matrices or
NetMHCIIpan (when available).

**Impact:** For therapeutic proteins where MHC-II-mediated CD4+ T-cell
responses are the primary concern, NetMHCIIpan remains the only
high-accuracy option and requires network access.

**Future:** The MHCflurry-MHCII extension is under development.  When
released, it will be integrated as a fourth tier in the prediction
hierarchy.

### 11.2 Model size and download

MHCflurry models require ~100 MB of disk space for the initial download.
This is a one-time cost — models are cached in `~/.mhcflurry/` and
loaded from disk on subsequent runs.

**Impact:** Environments with limited disk space or no network access
during setup cannot use MHCflurry.  The PSSM fallback remains available.

### 11.3 Memory usage

Loading MHCflurry models into memory requires ~500 MB RAM (TensorFlow
runtime + model weights).  This is negligible on modern workstations but
may be significant in constrained environments (e.g., CI runners,
serverless functions).

### 11.4 Prediction speed

MHCflurry batch prediction is fast (~2 s for 100 peptides × 5 alleles
on CPU), but significantly slower than PSSM (~0.1 s for the same input).
For very large proteins (>1000 aa) with many alleles (>50), the
computation can take 10–30 seconds on CPU.

**Mitigation:** Use `epitope_lengths=[9, 10]` instead of the default
`[8, 9, 10, 11]` to reduce the peptide × allele search space by ~25%.

### 11.5 Population coverage data accuracy

The allele frequency data in `mhcflurry_population.py` is compiled from
published sources (AFND, NMDP) but consists of **approximate estimates**
rounded to the nearest 0.5%.  The independence assumption in the coverage
formula ignores linkage disequilibrium between HLA loci.

**Impact:** Coverage estimates are upper bounds.  True coverage may be
5–10% lower due to LD effects.  Do not use for clinical decision-making
without validation.

### 11.6 GPU not required but beneficial

MHCflurry runs on CPU by default.  GPU acceleration is available when
TensorFlow detects a CUDA-compatible GPU, providing ~5–10× speedup for
large batch predictions.  However, GPU is not required and adds
significant installation complexity.

### 11.7 No B-cell epitope prediction

MHCflurry is an MHC binding / antigen processing predictor only.  It
does not provide B-cell epitope prediction.  B-cell epitopes continue
to use the existing Kolaskar-Tongaonkar, Parker, and Chou-Fasman scales.

**Future:** Integration of BepiPred-2.0 (LSTM-based) is planned as a
separate enhancement.

---

## References

1. O'Donnell TJ et al. (2018) MHCflurry: open-source class I MHC binding
   affinity prediction. *Bioinformatics* 34:2696–2703.
2. O'Donnell TJ et al. (2021) The MHCflurry 2.0 open-source major
   histocompatibility complex class I peptide binding prediction model
   and antigen processing model. *Nature Biotechnology* 39:1329–1336.
3. Reynisson B et al. (2020) NetMHCpan-4.1 and NetMHCIIpan-4.0:
   improved predictions of MHC antigen presentation by concurrent
   motif deconvolution and integration of MS MHC eluted ligand data.
   *Nucleic Acids Res* 48:W449–W454.
4. Gonzalez-Galarza FF et al. (2020) Allele Frequency Net Database
   (AFND) 2020 update: gold-standard data classification, open access
   genotype data and new query tools. *Nucleic Acids Res* 48:D946–D953.
5. Mori A et al. (1997) HLA allele and haplotype frequencies in the
   Japanese population. *Tissue Antigens* 50:354–362.
6. Middleton D & Williams F (2000) HLA allele and haplotype frequencies
   in African populations. *Eur J Immunogenet* 27:295–312.
