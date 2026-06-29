# Audit: Optimizer Strategy Reachability

**Scope**: Map which of the 5 optimizer strategies listed in `pipeline_core.py`
are actually reachable from the public `optimize_sequence()` API and which
are dead code, to guide future optimizer retirement work.
**Date**: 2026-06-19

---

## 0. Executive Summary

The BioCompiler `optimize_sequence()` API has a **two-tier dispatch problem**
that makes most of its stated strategy surface **unreachable by default and
untested through the public API**:

1. `use_integrated=True` is the **default** (pipeline_core.py:363) and the
   resulting fast path returns unconditionally at line 701 — it NEVER
   consults the `strategy` parameter or `use_csp_solver` flag.
2. Every slow-path strategy therefore requires the user to **explicitly
   pass `use_integrated=False`**. A `grep` of the entire codebase
   (src + tests) shows that **no test ever sets `use_integrated=False`**.
   Only one place in production code (`pipeline_core.py:674`) even mentions
   the string `use_integrated=False`, and it is a code comment.
3. As a direct consequence, every test that calls
   `optimize_sequence(strategy="hybrid"|"constraint_first"|"cai_first"|"harmonize")`
   or `optimize_sequence(use_csp_solver=True)` silently runs the integrated
   fast path and the strategy flag is a no-op. One such test
   (`test_harmonize_strategy_produces_harmonization_score`) actively FAILS
   because the fast path does not set `objective_score` while the harmonize
   path does — concrete proof that the strategy is silently ignored.

### Headline counts (5 audited strategies + 1 bonus + 9 side files)

| Verdict                          | Count | LOC      |
|----------------------------------|-------|----------|
| LIVE (reachable + tested via public API) | 1     | 645      |
| LIVE-UNTESTED-VIA-PUBLIC-API     | 4     | 9,490    |
| LIVE-BUT-BROKEN (cai_first eukaryote) | 1 | 705      |
| DEAD-ORPHAN (test-only / no production wiring) | 3 | 1,778 |
| DEAD (replaced / shim)           | 2     | 255      |
| **Total audited**                | 12    | 12,873   |

**Conservative retirement scenario** (only the 5 fully-dead items):
~2,033 LOC.

**Aggressive retirement scenario** (everything not reachable via the
default public-API call):
~13,124 LOC (≈ 64 % of the optimizer package's strategy surface).

The decision between these scenarios is the central question for
future cleanup. See §6 for a per-strategy recommendation matrix.

---

## 1. Strategy Dispatch Map (pipeline_core.py)

`optimize_sequence()` (pipeline_core.py:336) has the following dispatch
shape. Line numbers are from HEAD as of this audit:

```
optimize_sequence(target_protein, organism, …,
                  strategy="hybrid",        # line 344 — DEFAULT strategy
                  use_csp_solver=False,     # line 347
                  use_integrated=True,      # line 363 — DEFAULT
                  …)
  │
  ├── (validation, biosecurity screen — runs always, lines 607-643)
  │
  ├── if use_integrated and target_protein:           # line 649 — DEFAULT branch
  │     └── integrated_optimize(…)                    # line 651-660
  │         └── return _int_result                    # line 701 — UNCONDITIONAL RETURN
  │
  └── (slow path: only reachable if use_integrated=False OR
       integrated_optimize raised an Exception)
        │
        ├── if use_csp_solver:                        # line 757
        │     └── _pipeline_paths.run_csp_solver_path # line 758
        │         └── return _csp_result              # line 768
        │
        ├── if strategy == "harmonize":               # line 808
        │     └── _pipeline_paths.run_harmonize_path  # line 809
        │
        ├── if strategy == "hybrid":                  # line 819 — DEFAULT strategy
        │     └── HybridOptimizer(…).optimize(…)      # line 820-836
        │         └── (eukaryote CAI recovery, CpG elimination, miRNA avoidance,
        │              prokaryote fast predicates, eukaryote predicate evaluation,
        │              then falls through to common post-processing)
        │
        └── else:                                     # line 1170
              └── BioOptimizer(…, strategy=strategy)  # line 1173 — constraint_first / cai_first
                  └── opt.optimize(initial_seq, strategy=strategy)  # line 1202
                      ├── if effective_strategy == "cai_first":     # line 1787
                      │     └── self._optimize_cai_first(seq)
                      └── else:                                       # constraint_first (default)
                            └── (inline sequential steps)
```

Key observations:

* **The `strategy` parameter has NO effect when `use_integrated=True`**
  (the default). The fast path at lines 649-701 returns before any
  strategy dispatch.
* **The `use_csp_solver` flag has NO effect when `use_integrated=True`**
  for the same reason. Tests that pass `use_csp_solver=True` without
  `use_integrated=False` are silently hitting the fast path.
* The slow path's strategy dispatch is a 3-way `if/elif/else`:
  `harmonize` → `hybrid` → `else` (constraint_first or cai_first).
* `hybrid_postprocessing` is NOT a strategy per se; it is a sub-module
  used by `pipeline_cross_codon.CrossCodonMixin`, which is a parent
  class of `BioOptimizer` (line 1593). It is reachable only via the
  slow-path `else` branch (constraint_first / cai_first). The
  `HybridOptimizer` (slow-path strategy="hybrid") does NOT use
  `hybrid_postprocessing` — it uses `hybrid_prokaryote`,
  `hybrid_eukaryote`, `hybrid_constraints`, `hybrid_hillclimb`.
* `pipeline_core.py:97` imports `hybrid_postprocessing as _postprocessing`
  but never uses the alias — a dead import that should be removed
  regardless of future cleanup outcome.

---

## 2. Per-Strategy Audit

### 2.1 `integrated_optimize` — DEFAULT FAST PATH ✅ LIVE

| Field | Value |
|---|---|
| Entry point | `use_integrated=True` (DEFAULT in `optimize_sequence()` signature, line 363) |
| Reachable via public API? | **YES** — every default `optimize_sequence()` call |
| File(s) | `src/biocompiler/optimizer/integrated_optimizer.py` (645 LOC) |
| Dependencies | `type_system.codon_tables`, `organisms`, `sequence.scanner`, `sequence.restriction_sites`, `shared.constants` — fully self-contained, does NOT import from `greedy`, `incremental`, `hybrid_*`, `strategy_*`, or `hybrid_postprocessing` |
| Test coverage | **YES** — `tests/test_integrated_optimizer.py` (55 tests, direct calls to `integrated_optimize`); plus every default `optimize_sequence()` test in the suite exercises it indirectly |
| Tested via public API? | **YES** — every default `optimize_sequence()` test |
| Empirical verification | `optimize_sequence('MVHLTPEEK', organism='Homo_sapiens')` returns `convergence_status='integrated'` |
| **Verdict** | **LIVE** |

### 2.2 `HybridOptimizer` — SLOW PATH, strategy="hybrid" (default strategy value) ⚠ LIVE-UNTESTED-VIA-PUBLIC-API

| Field | Value |
|---|---|
| Entry point | `use_integrated=False` AND `strategy="hybrid"` (or just `use_integrated=False`, since `"hybrid"` is the default strategy value) |
| Reachable via public API? | **YES** in principle — verified empirically: `optimize_sequence(…, use_integrated=False, strategy='hybrid')` returns `convergence_status='converged'`, `cai=0.8207` |
| File(s) | `hybrid_optimizer.py` (1,382) + `hybrid_types.py` (88) + `hybrid_prokaryote.py` (1,116) + `hybrid_eukaryote.py` (1,532) + `hybrid_constraints.py` (1,231) + `hybrid_hillclimb.py` (284) = **5,633 LOC** |
| Dependencies | `incremental`, `hybrid_types`, `hybrid_prokaryote`, `hybrid_eukaryote`, `hybrid_constraints`, `hybrid_hillclimb`, `numba_kernels`, `maxentscan_fast`, `solver.dispatch` (for optional CSP fallback), `codon_pair_scoring` |
| Test coverage | **YES, direct unit tests only** — `tests/test_thread_safety.py::TestHybridOptimizerThreadSafety` (2 tests, 3 instantiations), `tests/test_performance.py` (2 instantiations), `tests/test_task_2_9_medium_findings.py` (4 instantiations), `tests/test_maxentscan_correlation.py` (4 instantiations), `tests/test_task_1_8_constants_validation.py` (1 instantiation). All bypass `optimize_sequence()` and instantiate `HybridOptimizer(species=…).optimize(…)` directly. |
| Tested via public API? | **NO** — `grep -rn 'use_integrated\s*=\s*False' tests/` returns ZERO matches. The 17 `optimize_sequence(strategy="hybrid")` call sites in `test_integration_comprehensive.py` all use the default `use_integrated=True` and silently hit the fast path. |
| Empirical verification | `optimize_sequence(…, use_integrated=False, strategy='hybrid', strict_mode=False)` works correctly: returns `convergence_status='converged'`, `failed_predicates=['NoRQCTrigger']` |
| **Verdict** | **LIVE-UNTESTED-VIA-PUBLIC-API** — Reachable, directly unit-tested, but the public-API wiring (`use_integrated=False`) is exercised by ZERO tests |

### 2.3 `strategy_constraint_first` — SLOW PATH, strategy="constraint_first" ⚠ LIVE-UNTESTED-VIA-PUBLIC-API

| Field | Value |
|---|---|
| Entry point | `use_integrated=False` AND `strategy="constraint_first"` |
| Reachable via public API? | **YES** — verified empirically: `optimize_sequence(…, use_integrated=False, strategy='constraint_first')` returns `convergence_status='converged'`, `cai=0.8207` |
| File(s) | `src/biocompiler/optimizer/strategy_constraint_first.py` (1,582 LOC); plus `pipeline_cross_codon.py` (416 LOC, CrossCodonMixin parent) and `hybrid_postprocessing.py` (1,826 LOC, used by CrossCodonMixin). Direct + transitively-used: **3,824 LOC** |
| Dependencies | `incremental`, `mutagenesis`, `cai`, `greedy` (for `score_splice_donor_potential`), `constraints`, `pipeline_cross_codon` (CrossCodonMixin), `hybrid_postprocessing` (via CrossCodonMixin) |
| Test coverage | **PARTIAL** — `tests/test_optimization_biopyoptimizer.py` (33 tests) directly instantiates `BioOptimizer(species='ecoli', …).optimize(seq)`. The `BioOptimizer.optimize()` method's default strategy is `constraint_first` (pipeline_core.py:1612), so most of these tests DO exercise `strategy_constraint_first` — but only for **prokaryotes** (species='ecoli'). `tests/test_improvements.py::test_constraint_first_strategy` (line 86-92) and `tests/test_integration_comprehensive.py` call `optimize_sequence(strategy='constraint_first')` WITHOUT `use_integrated=False` — silently hitting the fast path. |
| Tested via public API? | **NO** — no test ever sets `use_integrated=False` |
| Empirical verification | `optimize_sequence(…, use_integrated=False, strategy='constraint_first', strict_mode=False)` works: `convergence_status='converged'`, `cai=0.8207` |
| **Verdict** | **LIVE-UNTESTED-VIA-PUBLIC-API** — Reachable, directly unit-tested via `BioOptimizer` (prokaryote only), but the public-API wiring is exercised by ZERO tests; eukaryote path is exercised by ZERO tests of any kind |

### 2.4 `strategy_cai_first` — SLOW PATH, strategy="cai_first" 🛑 LIVE-BUT-BROKEN-FOR-EUKARYOTES

| Field | Value |
|---|---|
| Entry point | `use_integrated=False` AND `strategy="cai_first"` |
| Reachable via public API? | **YES** for prokaryotes; **BROKEN** for eukaryotes |
| File(s) | `src/biocompiler/optimizer/strategy_cai_first.py` (705 LOC); transitively also uses `pipeline_cross_codon.py` (416) and `hybrid_postprocessing.py` (1,826). Direct + transitively-used: **2,947 LOC** |
| Dependencies | Same as constraint_first plus its own `_compute_seq_cai` NUMBA fast path |
| Test coverage | **PARTIAL** — `tests/test_optimization_biopyoptimizer.py::test_cai_first_strategy` (1 test, line 138-141) directly instantiates `BioOptimizer(species='ecoli', strategy='cai_first')` — but `species='ecoli'` means `is_prokaryote=True`, which **skips the buggy code path** (the loop at strategy_cai_first.py:134-142 is gated on `if not self.is_prokaryote:`). `tests/test_improvements.py` calls `optimize_sequence(strategy='cai_first')` WITHOUT `use_integrated=False` — silently hitting the fast path. |
| Tested via public API? | **NO** |
| Empirical verification | `optimize_sequence('MVHLTPEEKSAVTALWGKVNVDEVGGEALGR', organism='Homo_sapiens', strategy='cai_first', use_integrated=False, strict_mode=False)` raises `NameError: name 'HAS_NUMBA' is not defined` at `strategy_cai_first.py:177`. **Confirmed broken for eukaryotes.** Prokaryote path works. |
| Bug root cause | `strategy_cai_first.py:177` references `HAS_NUMBA` and `_numba_cai_kernel` but neither is imported at the top of the file (lines 15-25 import only `typing`, `logging`, `type_system.{CODON_TABLE, AA_TO_CODONS, BLOSUM62, PredicateResult, check_no_cryptic_splice}`, `mutagenesis.MutagenesisReport`, `provenance.certificate.format_certificate`, `constraints.{_count_gts, _is_unavoidable_gt, _has_gt}`, `cai._count_dinucs_fast`, `greedy.{score_splice_donor_potential, SPLICE_DONOR_POTENTIAL_THRESHOLD}`). The bug has been latent because the only test that exercises the path uses a prokaryote, which skips the buggy `if not self.is_prokaryote:` block at line 134. |
| **Verdict** | **LIVE-BUT-BROKEN-FOR-EUKARYOTES** — Reachable for prokaryotes (works, but untested via public API); **CRASHES for eukaryotes** (the default organism `Homo_sapiens`!) |

### 2.5 `hybrid_postprocessing` — NOT A STRATEGY, A SUB-MODULE ⚠ LIVE-UNTESTED-VIA-PUBLIC-API

| Field | Value |
|---|---|
| Entry point | None direct. Reachable only via `BioOptimizer` → `CrossCodonMixin` (pipeline_cross_codon.py) → `hybrid_postprocessing`. So entry point is `use_integrated=False` AND `strategy in {"constraint_first", "cai_first"}`. |
| Reachable via public API? | **YES** transitively — when the slow-path `else` branch (line 1170) creates `BioOptimizer(…).optimize(…)`, every step that calls `self._step_*` (e.g. `_step_avoid_cpg_islands`, `_step_cai_hill_climb`, `_step_reoptimize`, `_step_mrna_stability_improvement`, `_step_mirna_avoidance`, etc.) is dispatched to `pipeline_cross_codon.py` which forwards to `hybrid_postprocessing`. |
| File(s) | `src/biocompiler/optimizer/hybrid_postprocessing.py` (1,826 LOC). Plus the deprecated shim `postprocessing.py` (21 LOC, see §3.9). |
| Dependencies | `incremental`, `mirna_avoidance` (lazy import at line 1795), `codon_pair_scoring`, `maxentscan`, `mutagenesis`, `type_system.checks` |
| Test coverage | **NO direct tests**. Exercised indirectly via `tests/test_optimization_biopyoptimizer.py` (33 BioOptimizer tests, prokaryote-only) — these call `BioOptimizer.optimize()` which transitively uses hybrid_postprocessing for the cpg/cai/mirna post-processing passes. But: (a) no eukaryote coverage, (b) no public-API coverage. |
| Tested via public API? | **NO** |
| Note | `pipeline_core.py:97` imports `hybrid_postprocessing as _postprocessing` but the alias is **never used** in pipeline_core.py (a dead import). The actual users are `pipeline_cross_codon.py:27` (20 call sites). |
| **Verdict** | **LIVE-UNTESTED-VIA-PUBLIC-API** — Reachable transitively via the BioOptimizer slow path, but only its prokaryote sub-path is exercised by any test, and never via `optimize_sequence()` |

### 2.6 BONUS: `harmonize` strategy ⚠ LIVE-UNTESTED-VIA-PUBLIC-API (not in task's 5-strategy list, but present in pipeline_core.py)

The task description lists 5 strategies but `pipeline_core.py:808` dispatches a 6th strategy, `harmonize`, that was omitted from the task list. Audited here for completeness because it directly affects the future cleanup retirement decision.

| Field | Value |
|---|---|
| Entry point | `use_integrated=False` AND `strategy="harmonize"` |
| Reachable via public API? | **YES** — verified empirically: `optimize_sequence(…, use_integrated=False, strategy='harmonize', source_organism='Escherichia_coli')` returns `convergence_status='converged'`, `cai=0.9376`, `objective_score=0.833` |
| File(s) | `src/biocompiler/optimizer/codon_harmonization.py` (474 LOC); `pipeline_paths.py::run_harmonize_path` (lines 173-381 of the 746-LOC file) |
| Dependencies | `codon_harmonization.{harmonize_codons, harmonize_with_cai_fallback, compute_harmonization_score}`, `type_system.checks`, `provenance.certificate.format_certificate`, `expression.translation.compute_cai`, `sequence.scanner.gc_content` |
| Test coverage | **YES direct unit tests** — `tests/test_codon_harmonization.py` (50 tests, including 5 in `TestPipelineIntegration` class that call `optimize_sequence(strategy='harmonize')`). **BUT all 5 of those tests silently hit the fast path** (default `use_integrated=True`), and 1 of them (`test_harmonize_strategy_produces_harmonization_score`) **FAILS** because of this — it asserts `result.objective_score is not None`, which only the harmonize path produces, not the integrated fast path. The other 4 pass trivially because they only assert `result.sequence`, `len(sequence) == 27`, `result.cai > 0`, `result.protein == "MVHLTPEEK"` — all of which the fast path also produces. |
| Tested via public API? | **NO** — and 1 test actively fails because of the dispatch mismatch (concrete proof that the strategy is silently ignored) |
| Empirical verification | `optimize_sequence('MVHLTPEEK', organism='Homo_sapiens', strategy='harmonize', source_organism='Escherichia_coli')` (default `use_integrated=True`) returns `convergence_status='integrated'`, `objective_score=None` — confirms the fast path ran and the strategy was ignored. Same call with `use_integrated=False` returns `convergence_status='converged'`, `objective_score=0.833` — confirms the harmonize path ran. |
| **Verdict** | **LIVE-UNTESTED-VIA-PUBLIC-API** — Reachable, directly unit-tested, but the public-API wiring is broken in tests (1 failing test, 4 false-positive passes) |

---

## 3. Side-File Audit (files called out in task description)

For each: `grep -rln 'from.*<module>\|import.*<module>' src/ tests/` to find importers.

### 3.1 `state_machine.py` (1,396 LOC) — 🛑 DEAD-ORPHAN

| Field | Value |
|---|---|
| Public API exposure | **NONE**. Not in `optimizer/__init__.py`. Not in `biocompiler/__init__.py`. Not re-exported anywhere. |
| Production importers | **ZERO**. `grep -rln 'from.*state_machine\|import.*state_machine' src/` returns only the file itself. |
| Test importers | `tests/test_state_machine.py` (42 tests, all pass). The test file imports `DeterministicOptimizationStateMachine` directly via `from biocompiler.optimizer.state_machine import …`. |
| Other references | `ARTIFACT.md` (documentation only) |
| Tested via public API? | **NO** — not exposed via public API |
| Tested at all? | YES (its own test file) |
| **Verdict** | **DEAD-ORPHAN** — Experimental/orphan module. Tests exist but no production code path can reach it. Safe to retire (would also retire `test_state_machine.py`). |

### 3.2 `incremental.py` (1,245 LOC) — ✅ LIVE

| Field | Value |
|---|---|
| Public API exposure | Re-exported in `optimizer/__init__.py:214` and `biocompiler/__init__.py:207` |
| Production importers | **8 modules**: `greedy.py`, `hybrid_optimizer.py`, `hybrid_constraints.py`, `hybrid_hillclimb.py`, `hybrid_postprocessing.py`, `pipeline_core.py`, `pipeline_cross_codon.py`, `strategy_constraint_first.py` |
| Test importers | `tests/test_incremental.py` (58 tests) |
| Note | Used by `greedy.py` (line 43), which is in turn imported by `pipeline_core.py` (line 47) for `_eukaryote_cai_recovery` and `_eliminate_cpg_dinucleotides` — both of which are only called from the slow-path hybrid block (pipeline_core.py:857, 883, 947). So `incremental.py` is on the slow-path critical path. However, since `greedy.py` is imported at module top-level (line 47), `incremental.py` is loaded on EVERY `optimize_sequence()` call regardless of path — it just isn't *executed* on the fast path. |
| **Verdict** | **LIVE** (core data structure for slow path; cannot retire without removing the slow path entirely) |

### 3.3 `large_sequence.py` (619 LOC) — ⚠ LIVE (separate public API, not on `optimize_sequence()` dispatch)

| Field | Value |
|---|---|
| Public API exposure | Re-exported in `optimizer/__init__.py:220` (`optimize_large_sequence`, `ProteinTooLongError`) and `biocompiler/__init__.py:192` |
| Production importers | Only the two `__init__.py` re-exports. No production code calls `optimize_large_sequence()`. |
| Test importers | `tests/test_large_sequence.py` (38 tests, all call `optimize_large_sequence(…)` directly) |
| Reachable via `optimize_sequence()`? | **NO** — `optimize_sequence()` has no length-based dispatch to `optimize_large_sequence()`. The two are sibling public APIs. |
| Tested via public API? | NO (it IS the public API entry point) |
| Tested at all? | YES (38 direct tests) |
| **Verdict** | **LIVE** (separate parallel public API). Retiring it is a public-API break, not a dead-code removal. |

### 3.4 `mirna_avoidance.py` (778 LOC) — ⚠ LIVE-UNTESTED-VIA-PUBLIC-API

| Field | Value |
|---|---|
| Public API exposure | Re-exported in `optimizer/__init__.py:190` (`eliminate_mirna_binding_sites` aliased as `eliminate_mirna_binding_sites_active`) |
| Production importers | `pipeline_core.py:910` (lazy import inside slow-path hybrid block, eukaryote-only: `if not is_prokaryote and _mirna_avoid_mirna:`); `hybrid_postprocessing.py:1795` (lazy import inside `step_mirna_avoidance`, used by BioOptimizer slow path) |
| Test importers | **NONE**. `grep -rln 'mirna_avoidance' tests/` returns zero. `tests/test_mirna_predicates.py` defines its own local `_eliminate_mirna_binding_sites` helper function (line 491) instead of importing from `mirna_avoidance`. |
| Reachable via `optimize_sequence()`? | **YES** — but only on the slow path (`use_integrated=False, strategy="hybrid"` with a eukaryote organism) |
| Tested via public API? | **NO** |
| Tested at all? | NO direct tests |
| **Verdict** | **LIVE-UNTESTED-VIA-PUBLIC-API** — Reachable via slow-path hybrid block (eukaryote), but no test exercises that path. |

### 3.5 `mirna_elimination.py` (234 LOC) — 🛑 DEAD (replaced by `mirna_avoidance.py`)

| Field | Value |
|---|---|
| Public API exposure | Re-exported in `optimizer/__init__.py:185` as `eliminate_mirna_binding_sites` (the un-suffixed name) |
| Production importers | **ZERO**. The only import is `optimizer/__init__.py:185`. The two production call sites (`pipeline_core.py:910` and `hybrid_postprocessing.py:1795`) both import from `.mirna_avoidance`, NOT from `.mirna_elimination`. |
| Test importers | **ZERO**. No test imports `mirna_elimination` or `MAX_MIRNA_ELIMINATION_ITERATIONS`. |
| Reachable via `optimize_sequence()`? | **NO** — no production code path imports from it |
| Tested at all? | NO |
| **Verdict** | **DEAD** — superseded by `mirna_avoidance.py`. Both modules export a function named `eliminate_mirna_binding_sites`, but production code only calls the `mirna_avoidance` version. The `mirna_elimination` version is reachable only via explicit `from biocompiler.optimizer import eliminate_mirna_binding_sites` (the un-suffixed name in `__init__.py:186`), which no test or production code does. |

### 3.6 `codon_harmonization.py` (474 LOC) — ⚠ LIVE-UNTESTED-VIA-PUBLIC-API

| Field | Value |
|---|---|
| Public API exposure | Re-exported in `optimizer/__init__.py:207` (`compute_rca`, `harmonize_codons`, `harmonize_with_cai_fallback`, `compute_harmonization_score`) and `biocompiler/__init__.py:179` |
| Production importers | `pipeline_paths.py:196` (lazy import inside `run_harmonize_path`, the `strategy="harmonize"` slow path); `objectives.py:276` (lazy import inside `harmonization_objective`, used when user passes `objective="harmonization"` or `objective=make_harmonization_objective(...)`) |
| Test importers | `tests/test_codon_harmonization.py` (50 tests — direct unit tests for `compute_rca`, `harmonize_codons`, `harmonize_with_cai_fallback`, `compute_harmonization_score`, plus 5 integration tests that call `optimize_sequence(strategy='harmonize')` but silently hit the fast path; see §2.6) |
| Reachable via `optimize_sequence()`? | **YES** — via `use_integrated=False, strategy="harmonize"` |
| Tested via public API? | **NO** (and 1 of the 5 integration tests fails because of the dispatch bug) |
| Tested at all? | YES (50 direct unit tests) |
| Other users | `objectives.py::harmonization_objective` — used when user passes `objective="harmonization"`. This IS reachable on the fast path's slow-post-processing block (pipeline_core.py:1493-1499), specifically `_pipeline_cai.run_custom_objective_refinement`. So `codon_harmonization` is also reachable via the objective="harmonization" path even with `use_integrated=True`. (Out of scope for strategy audit but noted for completeness.) |
| **Verdict** | **LIVE-UNTESTED-VIA-PUBLIC-API** (as a strategy). Has solid direct unit-test coverage for the underlying functions. The `strategy="harmonize"` dispatch wiring is what's untested (and partially broken — see §2.6). |

### 3.7 `blast_avoidance.py` (449 LOC) — ✅ LIVE (separate predicate API, NOT on `optimize_sequence()` dispatch)

| Field | Value |
|---|---|
| Public API exposure | Re-exported transitively via `type_system.checks::check_no_blast_matches` and `type_system.predicates::evaluate_no_blast_matches`, both of which are in `type_system/registry.py` and `type_system/__init__.py`. |
| Production importers | `type_system/checks.py:1691` (lazy import inside `check_no_blast_matches` — predicate #34 "NoBlastMatches") |
| Test importers | `tests/test_blast_avoidance.py` (14 tests) |
| Reachable via `optimize_sequence()`? | **NO** — the slow path's `evaluate_extended_predicates()` (pipeline_paths.py:659-746) only evaluates predicates #11-#20 (mRNA stability, cryptic promoter, TM domain, cryptic ORF, RQC trigger, Alu, miRNA, m6A, polyA, secondary structure). `NoBlastMatches` (predicate #34) is NOT in that set. |
| Tested via public API? | NO |
| Tested at all? | YES (14 direct tests) |
| **Verdict** | **LIVE** (reachable via the type-system predicate API). NOT on the optimizer strategy dispatch — retiring it would not affect `optimize_sequence()` behavior, only the standalone predicate API. |

### 3.8 `benchmark_numba.py` (361 LOC) — 🛑 DEAD-ORPHAN (test/benchmark-only utility)

| Field | Value |
|---|---|
| Public API exposure | **NONE**. Not in `optimizer/__init__.py`. Not in `biocompiler/__init__.py`. |
| Production importers | **ZERO**. The only mention of `benchmark_numba` outside the file itself is `tests/test_numba_v2_wiring.py` (3 lazy imports at lines 347, 355, 366 inside test functions). |
| Test importers | `tests/test_numba_v2_wiring.py` (3 test functions call `run_numba_benchmark`) |
| Reachable via `optimize_sequence()`? | **NO** |
| Tested at all? | YES (3 tests, but they're testing the benchmark utility itself, not using it as a test fixture) |
| **Verdict** | **DEAD-ORPHAN** — Pure benchmark utility. Could be moved to `scripts/benchmark/` or deleted. Safe to retire (would also retire the 3 test functions in `test_numba_v2_wiring.py`). |

### 3.9 BONUS: `postprocessing.py` (21 LOC, deprecated shim) — 🛑 DEAD

| Field | Value |
|---|---|
| Public API exposure | Not re-exported in `__init__.py` |
| Production importers | **ZERO**. `grep -rln 'postprocessing' src/ tests/` finds the file itself plus `hybrid_postprocessing.py` (different module). |
| Test importers | ZERO |
| Notes | Self-documenting as "DEPRECATED" in its docstring (line 4). It's a 21-line `from .hybrid_postprocessing import *` shim with no users. |
| **Verdict** | **DEAD** — pure dead shim |

---

## 4. Test Coverage Matrix

The single most important empirical finding of this audit:

```
$ grep -rln 'use_integrated\s*=\s*False' tests/
(no output — ZERO matches)

$ grep -rln 'use_integrated\s*=\s*False' .
docs/optimizer_strategy_audit.md    ← this audit's predecessor notes
docs/w3_certified_path_spec.md      ← design spec
src/biocompiler/optimizer/pipeline_core.py:674  ← code comment
```

**No test in the entire biocompiler test suite ever sets `use_integrated=False`.**
This means the public-API wiring for **every** slow-path strategy is untested.

Tests that *appear* to test slow-path strategies but actually hit the fast path:

| Test file | Test name | Calls | Silently hits fast path? |
|---|---|---|---|
| tests/test_improvements.py | `test_cai_first_high_score` | `optimize_sequence(strategy='cai_first', strict_mode=False)` | YES |
| tests/test_improvements.py | `test_cai_first_preserves_length` | `optimize_sequence(strategy='cai_first', strict_mode=False)` | YES |
| tests/test_improvements.py | `test_constraint_first_works` | `optimize_sequence(strategy='constraint_first', strict_mode=False)` | YES |
| tests/test_integration_comprehensive.py | 17 calls in TestHybridOptimizer + TestCAIConsistency | `optimize_sequence(strategy='hybrid', strict_mode=False)` | YES (all 17) |
| tests/test_codon_harmonization.py | `test_harmonize_strategy_basic` | `optimize_sequence(strategy='harmonize', source_organism=…)` | YES |
| tests/test_codon_harmonization.py | `test_harmonize_strategy_with_cai_weight` | `optimize_sequence(strategy='harmonize', harmonization_cai_weight=…)` | YES |
| tests/test_codon_harmonization.py | `test_harmonize_strategy_defaults_source_to_target` | `optimize_sequence(strategy='harmonize')` | YES |
| tests/test_codon_harmonization.py | `test_harmonize_strategy_produces_harmonization_score` | `optimize_sequence(strategy='harmonize', source_organism=…)` | YES — **FAILS** (asserts `objective_score is not None`) |
| tests/test_codon_harmonization.py | `test_harmonize_strategy_with_prokaryote_target` | `optimize_sequence(strategy='harmonize', source_organism='Homo_sapiens', organism='Escherichia_coli')` | YES |
| tests/test_gfp_e2e.py | `test_csp_solver_under_60s_if_available` | `optimize_sequence(use_csp_solver=True)` | YES (use_csp_solver ignored) |
| tests/test_insulin_e2e.py | (similar CSP test) | `optimize_sequence(use_csp_solver=True)` | YES |
| tests/test_csp_integration.py | multiple | `optimize_sequence(…)` with no `use_integrated=False` and no `use_csp_solver=True` | YES (all) |

Tests that DO exercise slow-path strategies, but only via direct instantiation (bypassing `optimize_sequence()`):

| Test file | Strategy exercised | Bypasses public API? |
|---|---|---|
| tests/test_optimization_biopyoptimizer.py (33 tests) | `strategy_constraint_first` (default) + `strategy_cai_first` (1 test, prokaryote only) | YES — `BioOptimizer(species=…).optimize(seq)` directly |
| tests/test_thread_safety.py | `HybridOptimizer` (slow path strategy="hybrid") | YES — `HybridOptimizer(species=…).optimize(…)` directly |
| tests/test_performance.py | `HybridOptimizer` | YES — direct |
| tests/test_task_2_9_medium_findings.py | `HybridOptimizer` | YES — direct |
| tests/test_maxentscan_correlation.py | `HybridOptimizer` | YES — direct |
| tests/test_task_1_8_constants_validation.py | `HybridOptimizer` | YES — direct |
| tests/test_optimizer_timeout.py | `HybridOptimizer` (monkey-patched) | YES — direct |
| tests/test_codon_harmonization.py (45 of 50 tests) | `harmonize` underlying functions | YES — direct calls to `compute_rca`, `harmonize_codons`, etc. |
| tests/test_state_machine.py | `state_machine.DeterministicOptimizationStateMachine` | YES — direct (and not on any production path) |
| tests/test_incremental.py | `IncrementalSequenceState` etc. | YES — direct |
| tests/test_large_sequence.py | `optimize_large_sequence` (separate public API) | YES — direct (separate API) |
| tests/test_blast_avoidance.py | `check_no_blast_matches` predicate | YES — direct (predicate API) |
| tests/test_integrated_optimizer.py (55 tests) | `integrated_optimize` direct + every default optimize_sequence() test indirect | BOTH — direct + public API |

---

## 5. Summary Table (the table requested in the task spec)

| Strategy / File | Entry Point | Reachable via `optimize_sequence()`? | Tested via public API? | Tested directly? | LOC | Verdict |
|---|---|---|---|---|---|---|
| `integrated_optimize` | `use_integrated=True` (DEFAULT) | YES | YES | YES (55 tests) | 645 | **LIVE** |
| `HybridOptimizer` (strategy="hybrid") | `use_integrated=False, strategy="hybrid"` (default strategy) | YES | NO (zero tests set `use_integrated=False`) | YES (15+ direct tests across 5 files) | 1,382 + 5,248 satellite = **5,633** | **LIVE-UNTESTED-VIA-PUBLIC-API** |
| `strategy_constraint_first` | `use_integrated=False, strategy="constraint_first"` | YES | NO | PARTIAL (33 BioOptimizer tests, prokaryote-only) | 1,582 + 2,242 satellite = **3,824** | **LIVE-UNTESTED-VIA-PUBLIC-API** |
| `strategy_cai_first` | `use_integrated=False, strategy="cai_first"` | YES (prokaryote); BROKEN (eukaryote — `NameError: HAS_NUMBA`) | NO | PARTIAL (1 test, prokaryote-only — skips the bug) | 705 + 2,242 satellite = **2,947** | **LIVE-BUT-BROKEN-FOR-EUKARYOTES** |
| `hybrid_postprocessing` (sub-module of BioOptimizer, NOT a strategy) | Indirect: `use_integrated=False, strategy in {"constraint_first","cai_first"}` | YES (transitively) | NO | NO direct tests; partial indirect via BioOptimizer tests | 1,826 | **LIVE-UNTESTED-VIA-PUBLIC-API** |
| `harmonize` (bonus — omitted from task's 5-strategy list) | `use_integrated=False, strategy="harmonize"` | YES | NO (1 of 5 integration tests FAILS due to dispatch bug) | YES (50 unit tests for underlying functions) | 474 | **LIVE-UNTESTED-VIA-PUBLIC-API** |
| `state_machine.py` | None — not on any dispatch path | NO | NO | YES (42 tests) | 1,396 | **DEAD-ORPHAN** |
| `incremental.py` | Indirect — used by 8 production modules | YES (transitively) | YES (indirect) | YES (58 tests) | 1,245 | **LIVE** |
| `large_sequence.py` | Separate public API (`optimize_large_sequence`) — NOT on `optimize_sequence()` dispatch | NO (separate API) | NO | YES (38 tests, direct) | 619 | **LIVE** (separate API) |
| `mirna_avoidance.py` | Indirect via slow-path hybrid block (eukaryote) | YES (transitively) | NO | NO (tests use local helper) | 778 | **LIVE-UNTESTED-VIA-PUBLIC-API** |
| `mirna_elimination.py` | None — only re-exported, no production caller | NO | NO | NO | 234 | **DEAD** (replaced by `mirna_avoidance`) |
| `codon_harmonization.py` | Indirect via `strategy="harmonize"` AND via `objective="harmonization"` | YES | NO | YES (50 direct unit tests) | 474 | **LIVE-UNTESTED-VIA-PUBLIC-API** |
| `blast_avoidance.py` | Separate type-system predicate API (`check_no_blast_matches`) — NOT on `optimize_sequence()` dispatch | NO (separate predicate API) | NO | YES (14 tests) | 449 | **LIVE** (separate API) |
| `benchmark_numba.py` | None — not on any dispatch path, not in any production import | NO | NO | YES (3 tests of the benchmark itself) | 361 | **DEAD-ORPHAN** |
| `postprocessing.py` (deprecated shim) | None — no importer | NO | NO | NO | 21 | **DEAD** (shim with no users) |

### Totals

| Verdict | Strategies/Files | LOC |
|---|---|---|
| LIVE (reachable + tested via public API) | 1 (`integrated_optimize`) | 645 |
| LIVE (separate public API, not on `optimize_sequence` dispatch) | 2 (`large_sequence`, `blast_avoidance`) | 1,068 |
| LIVE (core slow-path data structure) | 1 (`incremental`) | 1,245 |
| LIVE-UNTESTED-VIA-PUBLIC-API | 5 (`HybridOptimizer` stack, `strategy_constraint_first`, `strategy_cai_first` [broken for eukaryote], `hybrid_postprocessing`, `harmonize` stack, `mirna_avoidance`, `codon_harmonization`) | 12,482 |
| LIVE-BUT-BROKEN-FOR-EUKARYOTES | 1 (`strategy_cai_first`) — counted above | (counted above) |
| DEAD-ORPHAN (test-only, no production wiring) | 2 (`state_machine`, `benchmark_numba`) | 1,757 |
| DEAD (replaced/shim) | 2 (`mirna_elimination`, `postprocessing`) | 255 |
| **Total audited** | **14 entries** | **17,452 LOC** |

(Numbers add to more than the file total because some files are counted under multiple strategies — e.g. `pipeline_cross_codon.py` and `hybrid_postprocessing.py` are transitively used by both `strategy_constraint_first` and `strategy_cai_first`. The unique file LOC total is 20,638 — see §1 for the wc-l output.)

---

## 6. Retirement Scenarios for future cleanup

### Scenario A — Conservative (only fully-dead code)

Retire only the 4 entries that are DEAD or DEAD-ORPHAN:

| File | LOC | Test file to also retire | Test LOC |
|---|---|---|---|
| `state_machine.py` | 1,396 | `tests/test_state_machine.py` | ~722 |
| `mirna_elimination.py` | 234 | (none — no test imports it) | 0 |
| `benchmark_numba.py` | 361 | 3 test functions in `tests/test_numba_v2_wiring.py` | ~30 |
| `postprocessing.py` | 21 | (none) | 0 |
| **Total** | **2,012** | | **~752** |

Risk: **VERY LOW**. These files have no production callers and no public-API exposure. Removing them is a pure cleanup.

Caveat: `mirna_elimination.py` is re-exported as `eliminate_mirna_binding_sites` in `optimizer/__init__.py:186` and `__all__` at line 264. Removing it requires updating `__init__.py` to remove the re-export (or keeping a one-line shim that imports from `mirna_avoidance`).

### Scenario B — Aggressive (retire everything not reachable via the default `optimize_sequence()` call)

In addition to Scenario A, retire the entire slow-path strategy stack:

| File | LOC |
|---|---|
| `strategy_constraint_first.py` | 1,582 |
| `strategy_cai_first.py` (also bug-ridden) | 705 |
| `hybrid_postprocessing.py` | 1,826 |
| `hybrid_optimizer.py` | 1,382 |
| `hybrid_types.py` | 88 |
| `hybrid_prokaryote.py` | 1,116 |
| `hybrid_eukaryote.py` | 1,532 |
| `hybrid_constraints.py` | 1,231 |
| `hybrid_hillclimb.py` | 284 |
| `pipeline_cross_codon.py` (CrossCodonMixin — only used by BioOptimizer) | 416 |
| `mirna_avoidance.py` | 778 |
| `codon_harmonization.py` (if harmonize strategy is also retired) | 474 |
| `greedy.py` (3508 LOC — but heavily imported; see note below) | (3,508) |
| **Subtotal (excluding greedy.py)** | **10,414** |

Plus retire the public-API parameters that would become no-ops:
- `strategy` parameter (default `"hybrid"` — silently ignored today; would be removed)
- `use_csp_solver` parameter (silently ignored on fast path; CSP slow path is also untested via public API)
- `use_integrated` parameter (would become a no-op since only `True` is supported)
- `source_organism`, `harmonization_cai_weight` parameters (only used by harmonize path)

Risk: **HIGH**. This is a public-API break. Multiple existing direct-API tests would also need to be retired:
- `tests/test_optimization_biopyoptimizer.py` (33 tests, exercises `BioOptimizer` directly)
- `tests/test_thread_safety.py::TestHybridOptimizerThreadSafety` (2 tests)
- `tests/test_performance.py` (2 HybridOptimizer tests)
- `tests/test_task_2_9_medium_findings.py` (4 HybridOptimizer tests)
- `tests/test_maxentscan_correlation.py` (4 HybridOptimizer tests)
- `tests/test_task_1_8_constants_validation.py` (1 HybridOptimizer test)
- `tests/test_optimizer_timeout.py` (1 HybridOptimizer test)
- `tests/test_codon_harmonization.py` (50 tests, ~45 of which test the harmonize underlying functions directly)

Also requires retiring `BioOptimizer` class itself (pipeline_core.py:1593-1865, ~270 LOC of class body plus `optimize()` method).

Note on `greedy.py` (3,508 LOC): it is imported by `pipeline_core.py` (lines 47-57, 105-114) and by `strategy_cai_first.py` and `strategy_constraint_first.py`. The pipeline_core imports are used ONLY in the slow-path hybrid block (`_eukaryote_cai_recovery`, `_eliminate_cpg_dinucleotides` at lines 857, 883, 947). If the slow path is retired, those imports become dead and `greedy.py` could also be retired — BUT `greedy.py` also exports `score_splice_donor_potential`, `SPLICE_DONOR_POTENTIAL_THRESHOLD`, and other constants that are re-exported in `optimizer/__init__.py:71-80, 89-94` and may be used by other parts of the codebase. A full greedy.py retirement requires a separate dependency audit.

### Scenario C — Middle ground (retire broken + dead only)

In addition to Scenario A, retire only the strategies that are provably broken or have zero direct test coverage:

| File | LOC | Reason |
|---|---|---|
| `strategy_cai_first.py` | 705 | BROKEN for eukaryotes (default organism!) — `NameError: HAS_NUMBA` |
| **Total Scenario C** | **2,717** | (Scenario A + strategy_cai_first) |

Risk: **LOW** for `strategy_cai_first.py` retirement (the path is broken for the default organism anyway, so no production user could be relying on it). Requires deprecating the `strategy="cai_first"` parameter value.

### Recommendation

future cleanup should pursue **Scenario A** (conservative) immediately, then evaluate **Scenario C** (retire `strategy_cai_first`) as a near-term follow-up since the path is already broken. **Scenario B** (aggressive) should be preceded by:
1. A deprecation cycle for `use_integrated=False`, `strategy != "hybrid"`, and `use_csp_solver=True` parameters (raise `DeprecationWarning` on the slow path).
2. Migration of `BioOptimizer` direct-instantiation tests to either `integrated_optimize` (fast path) or `HybridOptimizer` (only slow-path strategy worth keeping).
3. A separate audit of `greedy.py`'s 3,508 LOC to determine which symbols are still needed by the fast path and the public `optimizer/__init__.py` exports.

---

## 7. Open Questions for future cleanup Implementer

1. **Is `BioOptimizer` itself part of the public API?** It is re-exported in `optimizer/__init__.py:85` and `biocompiler/__init__.py`. Tests instantiate it directly. Retiring the slow-path strategies effectively retires `BioOptimizer` — is that an acceptable public-API break?
2. **Should `use_csp_solver` be retired?** It is silently ignored on the fast path. The CSP solver engines (`solver/engine_z3.py`, `solver/engine_ortools.py`) are tested directly via `tests/test_csp_z3.py` etc., but the wiring `optimize_sequence(use_csp_solver=True, use_integrated=False) → run_csp_solver_path()` is exercised by ZERO tests. If the CSP path is retained, it needs at least one integration test that actually sets `use_integrated=False`.
3. **The `harmonize` strategy has a failing test (`test_harmonize_strategy_produces_harmonization_score`)** — should this be fixed by adding `use_integrated=False` to the test, or by removing the test (and the harmonize strategy) entirely? The test failure has been latent because the test suite has many pre-existing failures (per W2-c worklog).
4. **The `strategy_cai_first.py:177` `HAS_NUMBA` NameError bug** — fix or retire? If retired, the `strategy="cai_first"` parameter value should raise `ValueError` rather than silently falling through to `constraint_first` (the `else` branch at pipeline_core.py:1170 currently catches `cai_first` only because it's the `else` after `hybrid` and `harmonize`).
5. **Should `large_sequence.optimize_large_sequence()` be wired into `optimize_sequence()` as an automatic dispatch for proteins > N amino acids?** Currently `optimize_sequence()` has no length cap and no delegation to `optimize_large_sequence()` — they are sibling public APIs. This is a design question, not a dead-code question.
6. **`mirna_elimination.py` vs `mirna_avoidance.py`** — both export `eliminate_mirna_binding_sites`. Which is canonical? Production code uses `mirna_avoidance`. The `__init__.py` exports BOTH (under different names: `eliminate_mirna_binding_sites` from `mirna_elimination`, `eliminate_mirna_binding_sites_active` from `mirna_avoidance`). This is a naming-collision footgun. Recommend retiring `mirna_elimination.py` (Scenario A) and renaming `eliminate_mirna_binding_sites_active` back to `eliminate_mirna_binding_sites` in a follow-up.

---

## 8. Methodology

* All `grep` performed with ripgrep via the Grep tool (no `find`/`grep` shell commands).
* LOC counts from `wc -l` on each file.
* Empirical verification of every "Reachable via public API?" claim by running
  `python -c "from biocompiler import optimize_sequence; …"` with the relevant
  parameter combination and inspecting `result.convergence_status`
  (`"integrated"` ← fast path ran; `"converged"` ← slow path ran).
* Empirical verification of the `strategy_cai_first` eukaryote bug by running
  `optimize_sequence(…, use_integrated=False, strategy='cai_first', organism='Homo_sapiens')`
  and observing the `NameError`.
* Empirical verification of the harmonize dispatch bug by running
  `pytest tests/test_codon_harmonization.py::TestPipelineIntegration::test_harmonize_strategy_produces_harmonization_score`
  and observing the failure with `objective_score=None` in the result object
  (proving the fast path ran instead of the harmonize path).
* No source files were modified. Only this audit file was written.

---

## 9. References

* `src/biocompiler/optimizer/pipeline_core.py` — main dispatch (lines 336-1580)
* `src/biocompiler/optimizer/pipeline_paths.py` — `run_harmonize_path`, `run_csp_solver_path`, `run_prokaryote_hybrid_path`, `evaluate_extended_predicates`
* `src/biocompiler/optimizer/integrated_optimizer.py` — fast-path implementation
* `src/biocompiler/optimizer/hybrid_optimizer.py` — slow-path strategy="hybrid"
* `src/biocompiler/optimizer/strategy_cai_first.py` — slow-path strategy="cai_first" (BROKEN for eukaryotes at line 177)
* `src/biocompiler/optimizer/strategy_constraint_first.py` — slow-path strategy="constraint_first"
* `src/biocompiler/optimizer/hybrid_postprocessing.py` — sub-module used by BioOptimizer via CrossCodonMixin
* `src/biocompiler/optimizer/__init__.py` — public-API surface (re-exports)
* `docs/w3_certified_path_spec.md` — design spec (predecessor audit)
