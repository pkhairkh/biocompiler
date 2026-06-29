# DOC-16: CSP/SMT Solver for Codon Optimization

| Field    | Value            |
|----------|------------------|
| ID       | DOC-16           |
| Version  | 0.9.0            |
| Status   | Final            |
| Date     | 2026-03-05       |
| Parent   | DOC-10 §5        |
| Replaces | ADR-0008 (partially) |

---

## 1. Overview

### 1.1 Why a CSP Solver?

The BioCompiler's codon optimizer must select, for each amino acid position, one
of several synonymous codons such that the resulting mRNA satisfies *all*
biological correctness constraints simultaneously — no cryptic splice sites, GC
content in range, no restriction enzyme sites, no instability motifs, and CAI
above a user-specified threshold. This is a **constraint satisfaction problem**
(CSP): find an assignment of values (codons) to variables (positions) subject to
constraints, or prove that no such assignment exists.

As argued in DOC-10 §5, the CSP formulation provides three properties that no
optimization approach can match:

1. **Soundness.** Every returned solution satisfies all hard constraints by
   construction — not with high probability, but necessarily.
2. **Completeness.** If a feasible assignment exists, the solver finds one.
3. **Diagnosis.** When infeasible, the Minimal Unsatisfiable Subset (MUS)
   identifies exactly which constraints conflict, enabling targeted relaxation.

### 1.2 What It Replaces

The previous production solver was the **greedy multi-step optimizer**
(ADR-0008), which processes constraints sequentially:

```
CAI maximization → restriction sites → ATTTA → T-runs → GC → reconciliation
→ cryptic splice → final reconciliation
```

The greedy optimizer has known problems:

| Problem | Description |
|---|---|
| **Step ordering matters** | Running GC adjustment before cryptic-splice elimination can lock in codons that create splice sites; running them in reverse can violate GC. The fixed ordering is a heuristic that works for most proteins but fails for some. |
| **No global optimality** | The CAI-maximization step picks the highest-CAI codon per position, but subsequent steps may force substitutions that reduce CAI below what a globally optimal assignment would achieve. |
| **Reconciliation brittleness** | Reconciliation passes fix regressions introduced by earlier steps, but two passes may not suffice for deeply interacting constraints. In rare cases, the optimizer oscillates without converging. |
| **No infeasibility proof** | When the greedy optimizer fails, it emits warnings but cannot prove that *no* feasible assignment exists. The user does not know whether to relax constraints or try a different algorithm. |

The CSP solver addresses all four problems. Constraints are considered
simultaneously, eliminating step-ordering artifacts. The solver explores the
full search space (within time bounds), guaranteeing optimal CAI within the
feasible region. There is no reconciliation — constraints never regress because
they are never applied sequentially. And infeasibility is reported with a MUS,
not a list of warnings.

### 1.3 Solver Hierarchy

No single solver backend handles all constraint types efficiently. The
BioCompiler therefore uses a **two-tier solver hierarchy**:

```
┌─────────────────────────────────┐
│  Tier 1: OR-Tools CP-SAT        │  ← Primary: fast, automaton constraints
│  (default, all constraint types) │
├─────────────────────────────────┤
│  Tier 2: Greedy optimizer       │  ← Last resort: timeout / unavailable
│  (with warning, no guarantees)  │
└─────────────────────────────────┘
```

The primary backend (OR-Tools CP-SAT) handles 99% of real-world problems in
under 10 seconds. The greedy optimizer is the emergency fallback when CP-SAT
times out or is unavailable at runtime; it produces a best-effort result with
an explicit `SOLVER_FALLBACK` warning.

> **Note (Z3 backend removed):** Earlier versions of BioCompiler shipped a Z3
> SMT backend as an intermediate "Tier 2" fallback for continuous (MaxEntScan)
> constraints. The `engine_z3.py` module was **removed** in second-pass
> cleanup (see `CHANGELOG.md` and `solver/dispatch.py`). Z3's SAT engine had
> exponential worst-case complexity on codon-optimization problems and timed
> out for proteins longer than ~300 amino acids, so the tier added latency
> without practical benefit. An explicit `SolverBackend.Z3` request is now
> treated identically to the default OR-Tools path. Historical Z3 benchmark
> data is retained in §10.2 for reference only.

---

## 2. Theoretical Foundation

### 2.1 CSP Formulation (DOC-10 §5)

The gene design CSP is defined as follows:

- **Variables**: One integer variable `x[i]` per codon position `i ∈ {0, …, n-1}`.
- **Domain**: `Dom(x[i]) = {c₁, c₂, …, c_k}` where each `c_j` is a synonymous
  codon for the amino acid at position `i`. The domain size ranges from 1
  (Met, Trp) to 6 (Leu, Arg, Ser).
- **Hard constraints**: All constraints that must be satisfied for the solution
  to be valid (see §3.1).
- **Objective**: Maximize CAI (see §3.2).

The decision problem ("does a feasible assignment exist?") is NP-complete in
general, but the BioCompiler's constraint structure — small per-variable domains
(max 6 values) and mostly local constraints (spanning ≤ 10 codon positions) —
makes it highly tractable in practice.

### 2.2 CAI as a Linearizable Objective

The Codon Adaptation Index is defined as the geometric mean of relative
synonymous codon usage (RSCU) values:

```
CAI(m) = (∏_{i=0}^{n-1} w(x[i])) ^ (1/n)
```

where `w(c)` is the relative adaptiveness of codon `c`. Since geometric means
do not linearize directly, we maximize the log-CAI instead:

```
log CAI(m) = (1/n) × Σ_{i=0}^{n-1} log w(x[i])
```

This is a weighted sum where each codon choice contributes a fixed bonus
`log w(c)`. CP-SAT handles this natively via integer-weighted objectives; Z3
handles it via real-arithmetic optimization.

### 2.3 Relationship to DOC-10 §5

DOC-10 §5 establishes the principle of constraint satisfaction over optimization
and introduces MUS-based diagnosis. This document specifies the *implementation*
of that principle: the concrete constraint model, the solver backends, the
automaton encoding for forbidden substrings, and the MUS extraction algorithm.

---

## 3. Constraint Model

### 3.1 Hard Constraints

Hard constraints are Boolean predicates that must hold for a solution to be
accepted. The solver treats them as inviolable; any violation renders the
assignment infeasible.

| ID | Constraint | Scope | Description |
|---|---|---|---|
| H1 | `CorrectTranslation` | Per-position | `x[i]` must encode the amino acid at position `i`. Enforced by domain construction — only synonymous codons appear in `Dom(x[i])`. |
| H2 | `NoRestrictionSite(S)` | Sliding window (6–8 nt) | No subsequence of the assembled mRNA matches any restriction enzyme recognition sequence in set `S`. Cross-codon: spans up to 3 codon positions. |
| H3 | `NoInstabilityMotif` | Sliding window (5 nt) | No subsequence matches AUUUA (ARE) or contains ≥6 consecutive T/U residues. Cross-codon: spans 2 codon positions. |
| H4 | `GCInRange(lo, hi)` | Global | GC content of the assembled mRNA must lie in `[lo, hi]`. Encoded as a linear inequality on the sum of per-codon GC counts. |
| H5 | `NoCrypticSpliceDonor` | Sliding window (9 nt) | No 9-mer in the assembled mRNA has a MaxEntScan donor score above the cryptic threshold. Cross-codon: spans 4 codon positions. (See §5–6.) |
| H6 | `NoCrypticSpliceAcceptor` | Sliding window (23 nt) | No 23-mer has a MaxEntScan acceptor score above the cryptic threshold. Cross-codon: spans 8 codon positions. (See §5–6.) |
| H7 | `NoPrematureStop` | Per-position | No in-frame codon (except the last) is a stop codon. Trivially enforced by domain construction for most amino acids; relevant for selenocysteine contexts. |

### 3.2 Soft Objective: Maximize CAI

CAI maximization is the soft objective. Among all feasible assignments (those
satisfying H1–H7), the solver selects the one with the highest log-CAI. This
provides a principled tie-breaking rule: if multiple codon assignments satisfy
all constraints, the one that maximizes expression potential is preferred.

### 3.3 Constraint Interaction Matrix

Constraints interact when they share codon variables. The following matrix
shows which constraint pairs can conflict (i.e., tightening one makes the
other harder to satisfy):

| | H2 (Restriction) | H3 (Instability) | H4 (GC) | H5 (Splice Donor) | H6 (Splice Acceptor) |
|---|---|---|---|---|---|
| **H2** | — | Low | Medium | High | High |
| **H3** | Low | — | Low | Medium | Low |
| **H4** | Medium | Low | — | **Very High** | **Very High** |
| **H5** | High | Medium | **Very High** | — | High |
| **H6** | High | Low | **Very High** | High | — |

The **H4–H5/H6 interaction** is the most common source of infeasibility.
Cryptic splice site elimination often requires substituting a GC-rich codon
with an AT-rich synonym, which conflicts with a minimum GC constraint. This
is the interaction most frequently flagged by MUS analysis (see §7).

### 3.4 Cross-Codon Constraints

Constraints H2, H3, H5, and H6 are **cross-codon**: they depend on the
sequence of nucleotides spanning two or more adjacent codon positions. For
example, a restriction enzyme recognition site `CTGCAG` (PstI) straddles a
codon boundary:

```
position i:   ...CTG | CAG...   position i+1
              ←codon i→ ←codon i+1→
```

The codon at position `i` must end with `CTG` and the codon at position `i+1`
must start with `CAG`. Neither codon alone determines whether the restriction
site exists; the constraint involves a **joint assignment** to `x[i]` and
`x[i+1]`.

Cross-codon constraints are the primary reason why a simple per-position greedy
strategy fails: fixing `x[i]` without considering `x[i+1]` can create or
destroy a restriction site at the boundary. The CSP solver handles these
natively by posting constraints over multiple variables simultaneously.

---

## 4. Solver Backends

### 4.1 OR-Tools CP-SAT (Primary)

**Google OR-Tools CP-SAT** is a constraint programming solver using a
conflict-driven clause-learning (CDCL) search engine with linear arithmetic
and Boolean reasoning. It is the default backend for all CSP invocations.

**Why CP-SAT is primary:**

| Property | CP-SAT | Z3 | Greedy |
|---|---|---|---|
| Automaton/table constraints | Native (`AddAutomaton`) | Not native | Not applicable |
| Solve time (500 AA) | 2–8 seconds | 30–300 seconds | <1 second |
| Solve time (1000 AA) | 5–20 seconds | Often timeout | <1 second |
| Optimality guarantee | Yes (within tolerance) | Yes | No |
| Infeasibility detection | Yes | Yes | No |
| MUS extraction | Via deletion-based | Native | Not available |
| Determinism | Fully deterministic | Fully deterministic | Fully deterministic |
| Cross-codon constraints | Via sliding-window tables | Via quantifier-free formulas | Sequential heuristic |

**CP-SAT-specific encoding:**

- Each codon variable `x[i]` is an `IntVar` with domain `{0, 1, …, k-1}` where
  `k = |Dom(x[i])|`. A lookup table maps integer values to actual codons.
- GC content is encoded as a linear constraint: `Σ gc_weight[x[i]] ∈ [lo*n/3, hi*n/3]`.
- Restriction sites and instability motifs are encoded via **forbidden assignment
  tables** (Boolean literals asserting that certain value combinations for
  adjacent positions are prohibited).
- Cryptic splice site avoidance uses the **automaton constraint** (§5).
- The objective is `Maximize(Σ cai_weight[x[i]])` where `cai_weight[v]` is the
  integer-rounded `log w(c)` scaled by 10⁴ for integer arithmetic.

### 4.2 Z3 SMT (Removed)

> **The Z3 SMT backend (`engine_z3.py`) was removed in second-pass cleanup.**
> This section is retained for historical reference only. The current solver
> hierarchy is OR-Tools CP-SAT → greedy fallback (see §1.3). An explicit
> `SolverBackend.Z3` request is treated identically to the default OR-Tools
> path.

**Z3** was an SMT solver supporting quantifier-free nonlinear real arithmetic
(QF_NRA) and quantifier-free bit-vector logic (QF_BV). It previously served
as the fallback backend for problems where continuous constraints dominate.

**When Z3 was preferred (historical):**

- Problems with dense MaxEntScan threshold constraints where the CP-SAT
  discretization (§6) loses too much precision.
- Debugging: Z3's native MUS extraction (`z3.unsat_core()`) is often faster
  than CP-SAT's deletion-based approach for small unsatisfiable cores.
- Research: users who want to experiment with custom SMT-LIB constraints.

**Z3-specific encoding:**

- Each codon variable is an `IntSort` constant with domain constraints.
- Nucleotide composition is modeled via auxiliary `Bool` variables: for each
  codon position `i` and each nucleotide slot `j ∈ {0,1,2}`, `nt[i][j]` is a
  function of `x[i]` mapping codon index → nucleotide.
- MaxEntScan scores are expressed as `RealSort` arithmetic expressions over
  nucleotide indicator variables.
- Forbidden substrings are encoded via explicit transition constraints that
  simulate a DFA (less efficient than CP-SAT's `AddAutomaton` but more flexible).

**Z3 limitations for gene design:**

1. **Scalability.** Z3's SAT engine has exponential worst-case complexity on
   problems with many small-domain variables and local constraints — exactly
   the structure of codon optimization. Proteins longer than ~300 amino acids
   regularly exceed the 60-second budget (ADR-0008).

2. **No native automaton constraint.** Forbidden-substring constraints must be
   encoded as large conjunctions of Boolean implications, which blow up the
   clause database for long patterns (e.g., 23-mer splice acceptor contexts).

3. **Real arithmetic overhead.** MaxEntScan constraints involving continuous
  scores trigger Z3's NRA solver, which is significantly slower than its
   Boolean/integer modes.

### 4.3 Greedy Optimizer (Last Resort)

When CP-SAT fails (timeout, library unavailable, or internal error),
the system falls back to the greedy multi-step optimizer (formerly ADR-0008,
now folded into `integrated_optimize`).

**Fallback behavior:**

- The result object carries `solver_used = "greedy_fallback"` and a
  `SOLVER_FALLBACK` warning.
- The greedy result is *not* guaranteed to satisfy all constraints. The caller
  (typically COMP-05 Type System) re-checks all predicates and may still issue
  FAIL verdicts.
- If the greedy result satisfies all constraints, the certificate includes a
  note that optimality is not guaranteed.

**When fallback is triggered:**

| Condition | Action |
|---|---|
| CP-SAT times out (default 60s) | Fall back to greedy |
| CP-SAT not installed (`ImportError`) | Greedy with `SOLVER_UNAVAILABLE` warning |
| CP-SAT returns internal error | Fall back to greedy |
| Explicit `SolverBackend.Z3` requested | Treated as default OR-Tools path (Z3 removed) |

---

## 5. Automaton Constraints

### 5.1 Forbidden Substrings as DFA Acceptance

A key advantage of CP-SAT over Z3 for gene design is the **automaton
constraint**: `model.AddAutomaton(vars, start_state, accept_states, transitions)`.
This constraint asserts that the sequence of values taken by `vars` is accepted
by a deterministic finite automaton (DFA).

For forbidden substring constraints, we construct the complement: a DFA that
accepts all sequences *except* those containing a forbidden substring. The
construction is a standard Aho-Corasick automaton:

1. Build a trie of all forbidden substrings (e.g., all 6-8 nt restriction
   enzyme recognition sequences).
2. Add failure transitions (Aho-Corasick extension) to make the trie a DFA.
3. Mark all states that represent the end of a forbidden substring as
   **rejecting** (not in the accept set).
4. All other states are **accepting**.

The solver then enforces that the sequence of nucleotide values (derived from
codon assignments) traces a path through the DFA that ends in an accepting
state — guaranteeing that no forbidden substring appears.

### 5.2 Encoding Codon Sequences as Symbol Sequences

The automaton constraint operates on a sequence of integer variables. Since
our primary variables are codon-level (`x[i]` ∈ codon indices), we introduce
auxiliary **nucleotide-level variables**:

```
For each codon position i and nucleotide slot j ∈ {0, 1, 2}:
    nt[3*i + j] = nucleotide_lookup[x[i]][j]
```

where `nucleotide_lookup[c][j]` maps codon index `c` to its `j`-th nucleotide
(encoded as an integer 0–3 for {A, C, G, T}).

The automaton constraint is posted over the `nt` variables, not the `x`
variables directly. CP-SAT handles the channeling between `x[i]` and `nt[3*i+j]`
automatically via element constraints.

### 5.3 Combining Multiple Forbidden Patterns

Multiple forbidden substrings (restriction sites, instability motifs, cryptic
splice site 9-mers, cryptic splice acceptor 23-mers) are combined into a
**single Aho-Corasick automaton**. This is more efficient than posting separate
automaton constraints for each pattern because:

- The solver maintains one transition table instead of many.
- Shared prefixes between patterns (e.g., `GUAAGU` and `GUAAGC` for splice
  donors) are merged in the trie, reducing the state count.
- The Aho-Corasick construction runs in O(total pattern length) time.

For a typical problem with 10 restriction sites + ~50 forbidden 9-mers +
~20 forbidden 23-mers, the combined automaton has ~500–2000 states — well within
CP-SAT's capacity.

### 5.4 Handling Valine and GT-Containing Codons

All four Valine codons (GTT, GTC, GTA, GTG) contain the GT dinucleotide, which
is the splice donor consensus. When a Valine position falls in a context that
creates a strong cryptic splice donor, *no* synonymous substitution can remove
the GT. This is the well-known Valine impossibility (ADR-0009, ADR-0011).

The automaton constraint handles this correctly: if the DFA has no accepting
path through a Valine position's nucleotide triplets, the constraint is
locally infeasible, and MUS analysis will identify the Valine position and the
conflicting splice-site constraint. The user is then directed to the
type-directed mutagenesis engine (ADR-0009) for conservative amino acid
substitution proposals.

---

## 6. MaxEntScan Encoding

### 6.1 The Challenge: Continuous Scores as Discrete Constraints

MaxEntScan is a maximum-entropy model that assigns a continuous score to each
9-mer (donor) or 23-mer (acceptor) window. A high score indicates a strong
splice site. The constraint `NoCrypticSpliceDonor` requires that no 9-mer in
the mRNA has a MaxEntScan score above a threshold (default: 3.0 for donors,
3.0 for acceptors, configurable).

The challenge is that MaxEntScan scores are computed by a lookup table over
4⁹ = 262,144 possible 9-mers (donor) or 4²³ ≈ 70 billion possible 23-mers
(acceptor). The acceptor space is far too large for explicit enumeration.

### 6.2 Discretization Strategy

The BioCompiler uses a **threshold-based discretization** approach:

1. **Pre-compute** the set of 9-mers whose MaxEntScan donor score exceeds the
   threshold. For a threshold of 3.0, this set typically contains ~200–500
   9-mers out of 262,144. These are the **forbidden 9-mers**.

2. **Pre-compute** the set of 23-mers whose MaxEntScan acceptor score exceeds
   the threshold. To avoid enumerating 4²³ possibilities, we use a
   **window-decomposition** heuristic: compute the score for each 23-mer that
   actually appears in the candidate mRNA's neighborhood (i.e., 23-mers
   reachable by synonymous substitution within a sliding window). This reduces
   the enumeration to O(n × d⁸) where `d` is the max domain size — typically
   ~10,000–100,000 candidate 23-mers per problem, of which ~50–200 exceed
   the threshold.

3. **Feed forbidden subsequences** to the Aho-Corasick automaton (§5.3). The
   solver then avoids all nucleotide sequences that would produce a forbidden
   9-mer or 23-mer.

### 6.3 Dynamic Threshold Adjustment

When the solver reports infeasibility due to a MaxEntScan constraint in the
MUS, the user can relax the threshold. The system supports **dynamic
re-thresholding**: the caller re-invokes the solver with a higher threshold,
which reduces the set of forbidden 9-mers/23-mers, potentially restoring
feasibility.

```python
# Example: Relax MaxEntScan donor threshold from 3.0 to 5.0
config = SolverConfig(
    maxentscan_donor_threshold=5.0,   # was 3.0
    maxentscan_acceptor_threshold=5.0, # was 3.0
)
result = solve_with_csp(ir_seq, constraints, config=config)
```

The relaxation is a conscious engineering trade-off: a higher threshold means
weaker splice sites are tolerated, increasing the risk of cryptic splicing in
vivo. The type system (COMP-05) will issue an UNCERTAIN verdict for positions
whose MaxEntScan score is between the old and new thresholds, preserving
soundness while acknowledging the uncertainty.

### 6.4 Z3 Encoding for MaxEntScan

When Z3 is the backend, MaxEntScan constraints are encoded more directly:
each 9-mer window is constrained via a real-arithmetic expression that computes
the score as a sum of per-position log-odds terms. This avoids the
discretization step and captures the continuous nature of the score.

However, this encoding is significantly slower (Z3 must reason about real
arithmetic) and is only used when the user explicitly requests the Z3 backend
or when the discretized CP-SAT formulation produces a false-negative
infeasibility (i.e., CP-SAT says INFEASIBLE but Z3 finds a solution — possible
when the threshold is near a 9-mer boundary).

---

## 7. MUS Analysis

### 7.1 What Is a MUS?

A **Minimal Unsatisfiable Subset** (MUS) is a subset of constraints that is
(1) unsatisfiable — no assignment satisfies all constraints in the subset — and
(2) minimal — removing any single constraint from the subset restores
satisfiability.

When the CSP is infeasible, the MUS tells the user *exactly which constraints
conflict*. This is far more useful than a generic "no solution exists" message,
because it enables targeted relaxation.

### 7.2 MUS Example

Suppose the solver reports the following MUS for a 500-amino-acid protein:

```
MUS (4 constraints):
  1. H4: GCInRange(0.40, 0.60)
  2. H5: NoCrypticSpliceDonor (MaxEntScan threshold 3.0)
  3. H1: CorrectTranslation at position 247 (Valine)
  4. H5: NoCrypticSpliceDonor at window spanning positions 246-248
```

This means: the Valine at position 247 creates a GT dinucleotide that, in
combination with its flanking codons (positions 246 and 248), forms a 9-mer
with MaxEntScan score > 3.0. All four Valine codons contain GT, so
`CorrectTranslation` and `NoCrypticSpliceDonor` are inherently in conflict at
this position. The GC constraint is also in the MUS because the only codons
available for the flanking positions that would disrupt the splice site are
AT-rich, pushing GC below 40%.

The user can then make an informed decision:

- **Option A**: Relax the MaxEntScan threshold to 5.0 (accepting a weaker
  guarantee on splice avoidance).
- **Option B**: Allow a conservative amino acid substitution at position 247
  (e.g., Val → Ile, BLOSUM62 score +3 — see ADR-0009).
- **Option C**: Widen the GC range to [35%, 60%].

### 7.3 MUS Extraction Algorithms

| Backend | Algorithm | Typical Time |
|---|---|---|
| CP-SAT | Deletion-based: iteratively remove constraints and re-solve. If infeasibility persists, the constraint is not in the MUS. | 2–10× single-solve time |
| Z3 | Native: `z3.unsat_core()` extracts the MUS from the proof tree. | Near-instant (already computed) |

For CP-SAT, the deletion-based algorithm is:

```
MUS = all_constraints
for c in all_constraints:
    if solve(MUS \ {c}) == INFEASIBLE:
        MUS = MUS \ {c}      # c is not needed for infeasibility
    else:
        MUS = MUS ∪ {c}      # c is essential; keep it
return MUS
```

This requires at most `|constraints|` solver invocations. For typical problems
with 5–20 constraints, this adds 5–20 seconds to the solve time.

### 7.4 MUS-Based Relaxation Suggestions

The solver can optionally compute **relaxation suggestions** for each constraint
in the MUS:

| Constraint Type | Relaxation |
|---|---|
| `GCInRange(lo, hi)` | Widen the range by ±5% increments until feasible |
| `NoCrypticSpliceDonor` | Raise MaxEntScan threshold by 1.0 increments |
| `NoRestrictionSite(S)` | Remove the most restrictive enzyme from S |
| `CorrectTranslation` | Propose amino acid substitution via mutagenesis engine |

The user receives a structured report:

```json
{
  "mus": ["GCInRange(0.40,0.60)", "NoCrypticSpliceDonor(threshold=3.0)", ...],
  "relaxation_suggestions": [
    {"constraint": "GCInRange(0.40,0.60)", "action": "Widen to [0.35, 0.60]", "feasible_after": true},
    {"constraint": "NoCrypticSpliceDonor(threshold=3.0)", "action": "Raise to 5.0", "feasible_after": true}
  ]
}
```

---

## 8. API Reference

### 8.1 `solve_with_csp`

```python
def solve_with_csp(
    ir_seq: IRSeq,
    constraints: list[ConstraintSpec],
    config: SolverConfig | None = None,
) -> SolverResult:
    """Solve the codon optimization CSP using the best available backend.

    Tries backends in order: CP-SAT → Z3 → Greedy. Returns as soon as a
    backend produces a result (feasible or infeasible). Falls back to the
    next backend only on timeout or internal error.

    Args:
        ir_seq: IR-Seq record to optimize.
        constraints: Hard constraints that must be satisfied.
        config: Solver configuration. If None, uses SolverConfig() defaults.

    Returns:
        SolverResult containing either a feasible assignment or an
        infeasible report with MUS.

    Raises:
        SolverTimeoutError: All backends exceeded their time budgets.
        ConstraintDefinitionError: A constraint is malformed.
    """
    ...
```

### 8.2 `SolverConfig`

```python
@dataclass(frozen=True)
class SolverConfig:
    """Configuration for the CSP solver."""

    # ── Backend selection ──────────────────────────────────────────────
    backend: Literal["auto", "cp_sat", "z3", "greedy"] = "auto"
    # "auto" tries CP-SAT → Z3 → Greedy in order.

    # ── Time budgets ───────────────────────────────────────────────────
    cp_sat_timeout: float = 60.0     # seconds
    z3_timeout: float = 60.0         # seconds

    # ── MaxEntScan thresholds ──────────────────────────────────────────
    maxentscan_donor_threshold: float = 3.0
    maxentscan_acceptor_threshold: float = 3.0

    # ── MUS analysis ───────────────────────────────────────────────────
    compute_mus: bool = True          # Compute MUS on infeasibility
    compute_relaxation_suggestions: bool = True

    # ── Search parameters ──────────────────────────────────────────────
    num_workers: int = 8              # CP-SAT parallel workers
    log_search_progress: bool = False # Emit CP-SAT search logs

    # ── Fallback behavior ──────────────────────────────────────────────
    allow_greedy_fallback: bool = True
```

### 8.3 `SolverResult`

```python
@dataclass(frozen=True)
class SolverResult:
    """Result from the CSP solver."""

    # ── Status ─────────────────────────────────────────────────────────
    status: Literal["FEASIBLE", "INFEASIBLE", "TIMEOUT", "ERROR"]
    solver_used: str                  # "cp_sat", "z3", "greedy_fallback"

    # ── Feasible result (status == "FEASIBLE") ─────────────────────────
    assignments: list[CodonAssignment] | None = None
    optimized_sequence: str | None = None
    cai: float | None = None
    gc_content: float | None = None
    objective_value: float | None = None

    # ── Infeasible result (status == "INFEASIBLE") ─────────────────────
    mus: list[ConstraintSpec] | None = None
    mus_description: str | None = None
    relaxation_suggestions: list[RelaxationSuggestion] | None = None

    # ── Metadata ───────────────────────────────────────────────────────
    solve_time_seconds: float = 0.0
    warnings: list[str] = field(default_factory=list)

    @property
    def is_feasible(self) -> bool:
        return self.status == "FEASIBLE"
```

### 8.4 Supporting Types

```python
@dataclass(frozen=True)
class CodonAssignment:
    """An assignment of a specific codon to a position."""
    position: int          # 0-based codon position
    original_codon: str    # Original codon in the input sequence
    assigned_codon: str    # Optimized codon assignment
    amino_acid: str        # Amino acid encoded

@dataclass(frozen=True)
class ConstraintSpec:
    """A hard constraint for the CSP optimizer."""
    name: str                      # e.g., "NoCrypticSpliceDonor"
    constraint_type: str           # e.g., "splicing", "gc_content"
    parameters: dict[str, Any]     # e.g., {"threshold": 3.0}

@dataclass(frozen=True)
class RelaxationSuggestion:
    """A suggested constraint relaxation to restore feasibility."""
    constraint: ConstraintSpec
    action: str                    # Human-readable action description
    feasible_after: bool           # Whether the relaxed problem was verified feasible
```

---

## 9. Migration Guide: From Greedy to CSP

### 9.1 Why Migrate?

The CSP solver is the default as of v1.0.0. Existing code using the greedy
optimizer directly should migrate to `solve_with_csp` for the following
benefits:

| Benefit | Greedy | CSP |
|---|---|---|
| Global CAI optimality | No (greedy per-position) | Yes (within feasible region) |
| Infeasibility proof | No | Yes, with MUS |
| Cross-codon constraints | Sequential (may regress) | Simultaneous (no regression) |
| Determinism | Yes | Yes |
| Fallback | N/A | Falls back to greedy automatically |

### 9.2 API Migration

**Before (legacy default optimizer — module now deleted):**

The standalone `biocompiler.optimizer.greedy` module was removed when the
greedy strategy was folded into the unified `integrated_optimize` pipeline.
The equivalent current call is `optimize_sequence`, which runs the same
greedy multi-step sequence (CAI → restriction sites → ATTTA → T-runs → GC
→ cryptic splice → reconciliation) internally:

```python
from biocompiler import optimize_sequence

result = optimize_sequence(
    protein_seq,
    organism="Homo_sapiens",
    enzymes=["EcoRI", "XhoI"],
    gc_lo=0.40, gc_hi=0.60,
    cryptic_splice_threshold=3.0,
)
# result.sequence, result.warnings
```

**After (CSP solver):**

```python
from biocompiler.solver.dispatch import solve_with_csp, SolverConfig
from biocompiler.solver.types import ConstraintSpec

constraints = [
    ConstraintSpec(name="NoRestrictionSite", constraint_type="restriction",
                   parameters={"enzymes": ["EcoRI", "XhoI"]}),
    ConstraintSpec(name="GCInRange", constraint_type="gc_content",
                   parameters={"lo": 0.40, "hi": 0.60}),
    ConstraintSpec(name="NoCrypticSpliceDonor", constraint_type="splicing",
                   parameters={"threshold": 3.0}),
    ConstraintSpec(name="NoCrypticSpliceAcceptor", constraint_type="splicing",
                   parameters={"threshold": 3.0}),
]

result = solve_with_csp(
    ir_seq=ir_seq,
    constraints=constraints,
    config=SolverConfig(),
)

if result.is_feasible:
    print(f"Optimized: CAI={result.cai:.3f}, GC={result.gc_content:.2f}")
else:
    print(f"Infeasible. MUS: {[c.name for c in result.mus]}")
```

### 9.3 Behavioral Differences

| Aspect | Greedy | CSP |
|---|---|---|
| Result on feasible input | Always produces a sequence (may not satisfy all constraints) | Produces a sequence satisfying *all* hard constraints |
| Result on infeasible input | Best-effort sequence + warnings | `INFEASIBLE` status + MUS |
| CAI on feasible input | Typically 0.85–0.92 | Typically 0.88–0.95 (higher due to global optimization) |
| Solve time | <1 second | 2–20 seconds |
| When CP-SAT/Z3 unavailable | N/A (greedy is always available) | Falls back to greedy with `SOLVER_FALLBACK` warning |

### 9.4 Preserving Greedy as Explicit Backend

Users who prefer the greedy optimizer (e.g., for interactive use where speed
matters more than optimality) can select it explicitly:

```python
result = solve_with_csp(
    ir_seq=ir_seq,
    constraints=constraints,
    config=SolverConfig(backend="greedy"),
)
```

This bypasses CP-SAT and Z3 entirely, producing the same result as the old
`greedy_optimize` function but wrapped in the `SolverResult` interface.

---

## 10. Performance Benchmarks

### 10.1 Benchmark Setup

All benchmarks were run on a single machine with the following specifications:

- **CPU**: AMD EPYC 7763 (64 cores, 2.45 GHz)
- **RAM**: 256 GB DDR4
- **Python**: 3.12.2
- **OR-Tools**: 9.10.4027
- **Z3**: 4.13.0
- **OS**: Ubuntu 22.04 LTS

Proteins were selected from UniProt to cover a range of lengths and amino acid
compositions. All runs used default `SolverConfig` (60s timeout per backend).

### 10.2 Solve Times

> **Historical note:** The Z3 columns below are retained from before the Z3
> backend was removed (see §1.3). They are no longer reproducible with the
> current codebase and are kept only to document why Z3 was dropped (it timed
> out on proteins longer than ~500 AA).

| Protein | Length (AA) | CP-SAT Time | Z3 Time *(historical)* | Greedy Time | CAI (CSP) | CAI (Greedy) |
|---|---|---|---|---|---|---|
| Insulin | 51 | 0.3 s | 1.2 s | 0.01 s | 0.94 | 0.91 |
| GFP | 238 | 1.8 s | 45 s | 0.04 s | 0.92 | 0.89 |
| Cas9 | 1368 | 12 s | TIMEOUT | 0.15 s | 0.89 | 0.84 |
| Titin fragment | 3000 | 38 s | TIMEOUT | 0.30 s | 0.87 | 0.82 |
| p53 | 393 | 2.5 s | 120 s | 0.05 s | 0.93 | 0.90 |
| mCherry | 236 | 1.5 s | 38 s | 0.03 s | 0.91 | 0.88 |
| Antibody HC | 450 | 4.2 s | TIMEOUT | 0.07 s | 0.90 | 0.86 |

**Key observations:**

1. CP-SAT solves all test proteins within 60 seconds, including the 3000-AA
   titin fragment.
2. Z3 times out for proteins longer than ~500 amino acids, consistent with
   ADR-0008's findings.
3. CSP-optimized CAI is consistently 2–5% higher than greedy-optimized CAI,
   because the CSP solver can make globally optimal codon choices that the
   greedy per-position strategy misses.
4. Greedy is 50–300× faster than CP-SAT but produces lower-quality results.

### 10.3 Infeasibility Detection Time

| Scenario | Constraints | CP-SAT MUS Time | Z3 MUS Time |
|---|---|---|---|
| Overconstrained GC + splice | 6 | 1.2 s | 0.8 s |
| Impossible Valine position | 4 | 0.5 s | 0.3 s |
| Multiple conflicting restriction sites | 8 | 2.1 s | 1.5 s |
| Dense splice constraints (50+ forbidden 9-mers) | 12 | 4.8 s | 3.2 s |

Z3's native `unsat_core()` is faster for MUS extraction because it avoids
repeated solving. However, the total time (solve + MUS) favors CP-SAT because
CP-SAT's initial solve is much faster.

### 10.4 Scaling Behavior

```
Solve Time (seconds)
  40 ┤                                          ■ Titin (3000 AA)
     │
  30 ┤
     │
  20 ┤
     │                           ■ Cas9 (1368 AA)
  10 ┤
     │
   5 ┤           ■ Ab HC (450 AA)
     │
   2 ┤   ■ p53   ■ GFP  ■ mCherry
     │
   1 ┤   ■ Insulin
     │
   0 ┼───┬────┬────┬────┬────┬────┬────┬────┬───
       0  500  1000 1500 2000 2500 3000 3500 4000
                   Protein Length (AA)
```

CP-SAT solve time is approximately linear in protein length for the constraint
densities typical of gene design (O(n) with a small constant). This is because
most constraints are local (spanning ≤ 8 codon positions), and CP-SAT's
constraint propagation effectively decomposes the problem into loosely coupled
subproblems.

---

## 11. Known Limitations

### 11.1 Very Long Sequences (> 5000 AA)

For proteins exceeding 5000 amino acids, CP-SAT solve times may exceed the
60-second default timeout. Two mitigations are available:

1. **Windowed solving**: Partition the protein into overlapping windows of
   1000 AA (with 50 AA overlap), solve each window independently, and merge
   results at the boundaries. This sacrifices global optimality but completes
   within O(window_count × 10s).

2. **Increased timeout**: Set `cp_sat_timeout = 300` (5 minutes). For most
   proteins up to 5000 AA, CP-SAT solves within 2–3 minutes.

### 11.2 Timeout Behavior

When the solver times out:

- `SolverResult.status = "TIMEOUT"`
- `SolverResult.solver_used` indicates which backend timed out.
- If `allow_greedy_fallback = True`, the system automatically tries the next
  backend in the hierarchy. If all backends time out, the final status is
  `TIMEOUT` and no solution is returned.
- If `allow_greedy_fallback = False`, a `SolverTimeoutError` is raised.

### 11.3 MaxEntScan Precision

The CP-SAT discretization of MaxEntScan constraints (§6.2) is conservative: it
forbids any 9-mer whose score exceeds the threshold. This means the solver may
report INFEASIBLE when a solution exists that places the score *at* the
threshold boundary (within floating-point precision). In practice, this is rare
(affects < 0.1% of problems) and can be resolved by raising the threshold by
0.1.

### 11.4 Non-Determinism in Parallel CP-SAT

CP-SAT uses multiple workers by default (`num_workers = 8`). While the search
is deterministic in single-worker mode, parallel search may explore different
parts of the search space first on different hardware. The **optimal objective
value** is deterministic, but the **specific codon assignment** may differ
across runs. To guarantee bit-identical results, set `num_workers = 1`.

### 11.5 Z3 Memory Usage

For large proteins, Z3's internal clause database can grow to several GB before
timing out. The system monitors memory usage and aborts Z3 if it exceeds a
configurable limit (default: 4 GB). This prevents out-of-memory conditions on
machines with limited RAM.

### 11.6 Missing Solver Dependencies

If neither `ortools` nor `z3-solver` is installed, the system falls back to the
greedy optimizer. The `SolverResult.warnings` list will contain
`"SOLVER_UNAVAILABLE: Neither OR-Tools nor Z3 is installed. Using greedy
fallback with no optimality guarantee."` Users who want formal guarantees must
install at least one solver:

```bash
pip install ortools    # Recommended (CP-SAT)
# or
pip install z3-solver  # Alternative (Z3)
```

### 11.7 Constraint Expression Limits

The current constraint model (§3) does not support:

- **Tissue-specific splicing constraints** that depend on cell type — these
  require a parameterized MaxEntScan model per cell type, which is future work.
- **RNA secondary structure constraints** — these require ViennaRNA integration
  (planned) and a fundamentally different constraint encoding.
- **Multi-gene circuit constraints** — these are handled by COMP-08
  (Compositional Verifier), not by the per-gene CSP solver.

---

## Appendix A: Decision Flow

```
                ┌──────────────────┐
                │ solve_with_csp() │
                └────────┬─────────┘
                         │
                ┌────────▼─────────┐
                │ backend == auto? │
                └────┬────────┬────┘
                     │ Yes    │ No
          ┌──────────▼──┐  ┌──▼──────────────┐
          │ Try CP-SAT   │  │ Use specified    │
          │              │  │ backend directly │
          └──────┬───────┘  └──┬──────────────┘
                 │              │
          ┌──────▼───────┐     │
          │ Solved?       │     │
          └──┬───────┬───┘     │
        Yes  │       │ No      │
     ┌───────▼──┐ ┌──▼──────┐ │
     │ Return    │ │ Try Z3  │ │
     │ FEASIBLE  │ │         │ │
     └───────────┘ └──┬──────┘ │
                     │        │
              ┌──────▼──────┐  │
              │ Solved?      │  │
              └──┬──────┬───┘  │
            Yes  │      │ No   │
         ┌───────▼──┐ ┌─▼─────┴──────┐
         │ Return    │ │ Greedy       │
         │ FEASIBLE  │ │ fallback?    │
         └───────────┘ └──┬──────┬────┘
                       Yes │      │ No
                  ┌────────▼──┐ ┌─▼──────────┐
                  │ Return    │ │ Return      │
                  │ FALLBACK  │ │ TIMEOUT     │
                  └───────────┘ └─────────────┘
```


## Appendix B: Constraint-to-Backend Compatibility

Not every constraint type is supported equally by every solver backend. The following compatibility matrix documents which backends can handle each constraint natively, which require encoding transformations, and which are unsupported entirely. This matrix should be consulted when selecting a backend for a specific problem instance or when debugging solver failures.

| Constraint | CP-SAT | Z3 | Greedy |
|---|---|---|---|
| `CorrectTranslation` | Native (domain construction) | Native (domain construction) | Native (codon table lookup) |
| `NoRestrictionSite` | Native (`AddAutomaton` / forbidden assignment tables) | Encoded (DFA simulation via Boolean implications) | Heuristic (sequential scan + repair) |
| `NoInstabilityMotif` | Native (`AddAutomaton`) | Encoded (DFA simulation) | Heuristic (ATTTA/T-run scan) |
| `GCInRange` | Native (linear inequality on GC counts) | Native (real-arithmetic constraint) | Heuristic (post-hoc GC adjustment) |
| `NoCrypticSpliceDonor` | Native (forbidden 9-mer table + `AddAutomaton`) | Native (real-arithmetic MaxEntScan expression) | Heuristic (MaxEntScan scan + codon substitution) |
| `NoCrypticSpliceAcceptor` | Encoded (window-decomposition heuristic, see Section 6.2) | Native (real-arithmetic MaxEntScan expression) | Heuristic (MaxEntScan scan) |
| `NoPrematureStop` | Native (domain construction excludes stop codons) | Native (domain construction) | Native (codon table lookup) |
| `Maximize CAI` | Native (integer-weighted objective) | Native (real-arithmetic optimization) | Heuristic (per-position greedy CAI) |

**Legend:**

- **Native**: The constraint is expressed directly in the backend's native constraint language with no loss of precision. The solver can reason about it efficiently using purpose-built algorithms.
- **Encoded**: The constraint requires a transformation (e.g., DFA simulation via Boolean clauses, discretization of continuous scores) that may introduce minor precision loss. The encoding is sound but may be less efficient than a native implementation.
- **Heuristic**: The constraint is handled by a procedural scan-and-repair strategy that provides no formal guarantee of finding a solution even when one exists. The greedy optimizer may miss feasible solutions or produce assignments that violate the constraint.

**Implications for backend selection:**

- For problems dominated by `NoCrypticSpliceDonor` and `NoCrypticSpliceAcceptor` constraints with tight MaxEntScan thresholds, Z3 may produce more precise results because it avoids the discretization step. However, Z3's real-arithmetic reasoning is substantially slower than CP-SAT's integer-based approach for typical gene design problems, so this precision advantage only materializes for small proteins or when the user is willing to accept longer solve times.
- For problems with many `NoRestrictionSite` constraints spanning multiple enzymes, CP-SAT's native `AddAutomaton` support is significantly more efficient than Z3's Boolean-encoded DFA simulation. The Aho-Corasick automaton built by CP-SAT merges shared prefixes across all forbidden patterns, reducing the state count and accelerating constraint propagation.
- The greedy backend provides no formal guarantee for any cross-codon constraint. It should only be used as a last resort when formal solver backends are unavailable or have timed out. Any solution produced by the greedy optimizer must be independently verified by the type system before it can be trusted.
