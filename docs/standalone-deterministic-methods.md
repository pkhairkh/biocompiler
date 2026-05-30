# Avoiding Probability: Deterministic Methods for a Non-Deterministic Biology

## How to get deterministic guarantees without pretending biology is deterministic

---

## The Problem

The revised framework embraced probability: PCFGs for splicing, HMMs for translation, calibration-based validation. This works, but it creates a tension. The compiler metaphor is powerful precisely because compilers provide **deterministic guarantees** (type safety, memory safety, well-typed programs don't go wrong). A probabilistic compiler is a category error — compilers prove things, they don't estimate them.

The original paper's mistake was trying to make *biology* deterministic (by averaging away noise). That fails because stochasticity is not noise — it's the biology. But there's a third option: **keep the biology non-deterministic, but make the *framework's answers* deterministic**.

This is not a word game. It's a well-established move in computer science, with concrete mathematical machinery behind it. The key insight:

> **You don't need a deterministic system to get deterministic answers. You need to ask the right questions — questions whose answers are deterministic even when the system is not.**

---

## The Core Move: Change the Question

| Probabilistic Question | Deterministic Equivalent |
|------------------------|--------------------------|
| "What is the probability that mRNA *m* produces protein *p*?" | "Is protein *p* in the **set of all possible outputs** of mRNA *m*?" |
| "What is the probability that this splice site is functional?" | "Is this splice site **guaranteed to be functional** / **guaranteed to be non-functional** / **uncertain**?" |
| "Find the mRNA that maximizes P(correct splicing) × P(high expression)" | "Find an mRNA that **satisfies all hard constraints**: correct splicing, high expression, no cryptic sites" |
| "What is the distribution of folding outcomes?" | "What is the **range** of possible folding energies?" |

Same biology, different questions. The answers to the right-column questions are deterministic (yes/no, set membership, bounds) even though the underlying process is stochastic.

This move has a name in computer science: **abstract interpretation**.

---

## Method 1: Abstract Interpretation

### What It Is

Abstract interpretation is the mathematical foundation of all static program analysis. It was developed by Patrick and Radhia Cousot in 1977 to answer questions about programs *without running them*. The key idea: instead of computing the exact (possibly non-deterministic) behavior of a program, compute a **sound over-approximation** — an "abstract" behavior that captures everything the program *might* do, without computing exact probabilities.

The classic example: a compiler doesn't need to know the exact value of variable `x` at runtime to prove that array access `a[x]` is safe. It only needs to know that `x` is in the range `[0, len(a)-1]`. The range is an **abstract value** — a sound over-approximation of the concrete value. The analysis is deterministic: for any program, the abstract interpreter produces a unique abstract value for each variable.

### Applied to Protein Synthesis

Instead of computing P(splice isoform *i* | mRNA *m*, cell type *c*), compute the **abstract domain** of all possible splicing outcomes:

| Abstract Domain | What It Captures | Example |
|----------------|-----------------|---------|
| **Set of possible isoforms** | `{isoform₁, isoform₂, isoform₃}` | "This mRNA can produce one of these three proteins" |
| **Interval of inclusion rates** | `[0.85, 1.0]` for exon 5 | "Exon 5 is included in 85-100% of transcripts" |
| **Three-valued classification** | `{INCLUDED, EXCLUDED, UNCERTAIN}` per exon | "Exon 5 is definitely included; exon 7 is uncertain" |
| **Bounds on folding energy** | `[-45.2, -38.7]` kcal/mol | "The protein will fold with energy in this range" |

Each of these is a **deterministic answer**. Given an mRNA and a cell type, the abstract interpreter produces a unique abstract value. No probabilities involved.

### Why This Is Not Just Renaming Probability

The critical distinction: abstract interpretation is **sound by construction**. If the abstract interpreter says "exon 5 is definitely included," then exon 5 is included in *every possible* execution of the splicing process — not just with high probability, but *necessarily*. If it says "exon 5 is uncertain," then there exists at least one execution where exon 5 is excluded. The three-valued answer (yes/no/maybe) is **complete**: it covers all cases, and each case has a precise meaning.

This is fundamentally different from a probabilistic model that says P(exon 5 included) = 0.97. That probability is an estimate, and it could be wrong. The abstract interpretation's answer is a guarantee, and it cannot be wrong (if the abstract domain is correctly designed).

### The Price

Sound over-approximation means the framework is **conservative**. It will sometimes say "uncertain" when the true answer is "yes" or "no." This is the price of determinism without probabilities. The framework never gives a wrong answer, but it sometimes gives an incomplete answer.

For many applications, this is acceptable — even preferable:
- **Gene design**: You want *guarantees* that your synthetic gene won't be mis-spliced, not a 97% probability.
- **Safety-critical applications** (gene therapy, vaccine design): A conservative answer ("uncertain — don't use this mRNA") is better than a confident-but-wrong answer.
- **Compositional verification**: When combining multiple components, you need *guarantees*, not probabilities. You cannot compose "97% safe" components and get a meaningful result for the combination.

---

## Method 2: Type Systems for Biological Correctness

### What It Is

Programming language type systems are the most successful example of deterministic reasoning about non-deterministic systems. A type checker proves properties of programs *without executing them* — and without knowing the inputs. "Well-typed programs don't go wrong" (Milner, 1978) is a deterministic guarantee about a system whose runtime behavior is non-deterministic (it depends on inputs).

### Applied to mRNA Design

Define a **type system** for mRNA sequences. The types encode biological correctness properties:

| Type | Meaning | Checkable How? |
|------|---------|----------------|
| `SpliceCorrect(C)` | Guaranteed to splice correctly in cell type C | Splice site grammar + known splicing factor concentrations |
| `NoCrypticSplice` | Contains no sequences matching splice site consensus | Scan against splice site PWMs |
| `CodonAdapted(O, θ)` | Codon adaptation index ≥ θ for organism O | Lookup codon usage table, compute CAI |
| `GCInRange(lo, hi)` | GC content in [lo%, hi%] | Count G/C nucleotides |
| `NoRestrictionSite(S)` | No restriction enzyme cut sites from set S | Scan against restriction enzyme database |
| `InFrame` | All ORFs are in the correct reading frame | Grammar check: start codon alignment |

The type checker is **deterministic**: given an mRNA sequence and a type, it either passes or fails. No probabilities. If an mRNA has type `SpliceCorrect(HEK293T) ∧ NoCrypticSplice ∧ CodonAdapted(H_sapiens, 0.8)`, then it is *guaranteed* to satisfy all three properties.

### Subtyping as Biological Refinement

Subtyping captures the relationship between cell types:

```
SpliceCorrect(HEK293T) <: SpliceCorrect(HumanCell)
SpliceCorrect(HumanCell) <: SpliceCorrect(MammalianCell)
```

If an mRNA is guaranteed to splice correctly in HEK293T cells, it is also guaranteed to splice correctly in any supertype (more general cell type). This enables **hierarchical verification**: verify the mRNA for a specific cell type, and you get guarantees for all more general cell types for free.

### Type Inference for Gene Design

Instead of checking whether a given mRNA has a desired type, the compiler can **infer** an mRNA that has the desired type. This is the inverse problem, reframed as type inference:

**Input**: Target protein + type constraints (cell type, organism, CAI threshold, etc.)
**Output**: An mRNA sequence that is well-typed — i.e., guaranteed to satisfy all constraints.

This is exactly what a compiler does when it infers types for a program. The inference is deterministic: either a well-typed solution exists (and the compiler finds one) or it doesn't (and the compiler reports a type error). No probabilities needed.

### What Makes This Novel

No current gene design tool provides type-system guarantees. Codon optimization tools are heuristic: they produce "good" designs but cannot *prove* they are correct. A type-checked mRNA design is **proven correct by construction** — not with high probability, but necessarily.

---

## Method 3: Constraint Satisfaction, Not Optimization

### The Move

Instead of maximizing P(correct splicing) × P(high expression) × P(no immunogenicity), solve a **constraint satisfaction problem (CSP)**:

```
FIND: mRNA sequence m
SUCH THAT:
    m translates to target protein P
    m has no cryptic splice sites (hard constraint)
    m has CAI ≥ 0.8 for Homo sapiens (hard constraint)
    m has GC content in [40%, 60%] (hard constraint)
    m has no restriction sites from set S (hard constraint)
    m has no RNA instability motifs (hard constraint)
```

If a solution exists, the solver finds one. If no solution exists, the solver reports infeasibility. Either way, the answer is deterministic.

### Comparison

| Probabilistic Optimization | Constraint Satisfaction |
|---------------------------|------------------------|
| Maximize expected quality | Guarantee all constraints are met |
| Solution quality is a score | Solution quality is binary: satisfies or doesn't |
| No hard guarantees | Hard guarantees by construction |
| May produce designs that fail with low probability | Cannot produce designs that violate any constraint |
| Requires probability estimates for all terms | Requires only yes/no predicates for all terms |

### When Constraints Are Too Tight

If no solution satisfies all constraints, the framework does **constraint relaxation** — not by converting to probabilities, but by identifying the *minimal* set of constraints that must be relaxed:

```
INFEASIBLE. Minimal conflicting constraints:
  - CAI ≥ 0.9 (can be relaxed to 0.78 without violating other constraints)
  - No cryptic splice sites (cannot be relaxed: 2 sites found, both in exon 3)
```

This is deterministic debugging — the compiler tells you exactly which constraints are in conflict and what the options are. No probabilities needed.

### What Makes This Novel

Current codon optimization tools use weighted scoring functions (which implicitly treat constraints as soft probabilities). A CSP-based approach treats constraints as **hard**, which is what you want for safety-critical gene design. No publicly available gene design tool uses CSP.

---

## Method 4: Non-Deterministic Automata (Without Probabilities)

### The Move

Non-deterministic finite-state transducers (NFSTs) are *non-deterministic* but not *probabilistic*. An NFST maps each input to a **set** of possible outputs, not a probability distribution over outputs. The set computation is deterministic: for any input, the NFST produces a unique set of outputs.

### Applied to Splicing

A splicing NFST:
- Input: pre-mRNA sequence
- Output: **set of all possible splice isoforms** (given known splice sites and alternative splicing rules)
- Computation: deterministic set enumeration

No probabilities. The NFST tells you "these are the possible isoforms" without ranking them. If you want to rank them, you need additional information — but the framework doesn't require it.

### Composition Without Probability

NFSTs compose deterministically: if NFST₁ maps input *x* to set S₁, and NFST₂ maps each element of S₁ to a set S₂(y), then the composition NFST₁ ∘ NFST₂ maps *x* to ⋃{S₂(y) | y ∈ S₁}. This is a deterministic computation.

This means the entire pipeline (splicing NFST ∘ translation NFST ∘ [folding as external call]) composes without probabilities. Each stage produces a set of possible outputs; the composition produces the set of all possible final outputs.

### When Sets Get Large

The set of possible isoforms for a gene with many alternative exons can be large (hundreds to thousands). But:
1. For gene design, you typically want to *avoid* alternative splicing, so the set should be a singleton (one correct isoform). The framework verifies this.
2. For analysis, you can use **symbolic set representation** (e.g., BDDs — binary decision diagrams) to represent large sets compactly. BDDs have been used for decades in hardware verification to represent sets of billions of states.

---

## Method 5: Game Semantics — Worst-Case Guarantees

### The Move

Model the cellular context as an **adversary** that chooses the worst-case outcome. The compiler guarantees correctness *even in the worst case*.

| Question | Probabilistic Answer | Game-Semantic Answer |
|----------|---------------------|---------------------|
| Will this mRNA splice correctly? | "With probability 0.97" | "Yes, even in the worst-case cellular context in class C" |
| Will this protein fold correctly? | "With probability 0.85" | "Yes, for all cellular conditions in class C" |

The universal quantifier ("for all") replaces the probability. The answer is deterministic.

### When This Is Too Conservative

Worst-case analysis is conservative: it may reject mRNAs that would work in practice. But this is the right behavior for:
- **Gene therapy**: You're injecting DNA into a patient. You want guarantees, not probabilities.
- **Vaccine design**: You're giving a synthetic gene to millions of people. Conservative is correct.
- **Industrial bioproduction**: You're running a bioreactor for weeks. A mis-splicing event that happens 1% of the time will occur thousands of times.

For less conservative guarantees, you can **refine the adversary class**:
- Adversary class "all mammalian cells" → very conservative
- Adversary class "HEK293T cells" → less conservative
- Adversary class "HEK293T cells under standard culture conditions" → even less conservative

The more precisely you specify the cellular context, the less conservative the guarantee. But the guarantee is always deterministic.

---

## Method 6: Three-Valued Logic

### The Move

Replace probabilities with a three-valued logic: **YES**, **NO**, **MAYBE**.

| Three-Valued Answer | Meaning | Action |
|---------------------|---------|--------|
| YES | Property is guaranteed to hold | Proceed |
| NO | Property is guaranteed not to hold | Reject |
| MAYBE | Property might or might not hold | Investigate further / add more constraints |

This is exactly what abstract interpretation produces. It's also what a type checker produces (well-typed, ill-typed, or type-error-that-can-be-resolved-with-more-information).

### Why Three Values Beat Probabilities

1. **Composability**: Three-valued logic composes cleanly. P(A and B) requires knowing P(A) and P(B) and their dependence. Three-valued logic: YES ∧ YES = YES, YES ∧ MAYBE = MAYBE, MAYBE ∧ MAYBE = MAYBE. No dependence modeling needed.

2. **No calibration needed**: Probabilistic models must be calibrated (the predicted probabilities must match observed frequencies). Three-valued logic requires no calibration — the answers are correct by construction (sound abstract interpretation).

3. **No training data needed**: Probabilistic models require training data to estimate parameters. Three-valued logic requires only the *rules* (splice site consensus, genetic code, codon usage tables) — all of which are already known.

4. **Interpretable**: "MAYBE" is more interpretable than "P = 0.63." It tells the user: "The framework cannot determine the answer from the available information. You need more specific constraints or more knowledge about the cellular context."

---

## The Unified Framework: Deterministic by Construction

All six methods share a common structure:

```
Input:  mRNA + specification (type constraints, correctness properties, cellular context class)
        │
        ▼
┌─────────────────────────┐
│  Abstract Interpreter    │  (deterministic computation)
│  / Type Checker          │
│  / Constraint Solver     │
│  / NFST Composition      │
│  / Game Analyzer          │
│  / Three-Valued Evaluator│
└────────┬────────────────┘
         │
         ▼
Output: YES / NO / MAYBE  +  (if YES: a verified design)
                           +  (if NO: a counterexample / conflict set)
                           +  (if MAYBE: a set of uncertain regions)
```

The computation is deterministic. The input produces a unique output. No randomness, no probabilities, no sampling.

The biology underneath is still non-deterministic — but the framework doesn't need to model the non-determinism probabilistically. It models it **set-theoretically** (what are all possible outcomes?) or **logically** (what properties are guaranteed?) or **type-theoretically** (what types are satisfied?).

---

## Comparison: Three Approaches

| Aspect | Original Paper (Fake Determinism) | Revised Framework (Probability) | This Document (True Determinism) |
|--------|-----------------------------------|--------------------------------|----------------------------------|
| Model of biology | Deterministic (averaged) | Probabilistic | Non-deterministic (sets of outcomes) |
| Framework's answers | Single output | Probability distribution | Deterministic (YES/NO/MAYBE) |
| Guarantees | Wrong (averaging destroys info) | Soft (calibrated probabilities) | Hard (sound by construction) |
| Composition | Trivial (single output) | Complex (dependence modeling) | Clean (set union / logical AND) |
| Validation | Circular (parsing success) | Calibration | Proof (abstract interpretation is sound) |
| When wrong | Silently (wrong average) | Measurably (miscalibration) | Never (sound over-approximation) |
| When uncertain | Pretends not to be | Quantifies uncertainty | Reports MAYBE explicitly |
| Training data | Not needed (but needed for induction) | Needed (parameter estimation) | Not needed (only known rules) |
| Novelty | Low (codon opt. is solved) | Medium (PCFGs for splicing) | High (no one has done abstract interpretation for gene design) |

---

## Concrete Example: Splicing-Aware Gene Design, Deterministically

**Problem**: Design an mRNA for protein P that is guaranteed to splice correctly in HEK293T cells.

**Step 1: Define the type**
```
TargetType = SpliceCorrect(HEK293T) ∧ NoCrypticSplice ∧ CodonAdapted(H_sapiens, 0.8) ∧ GCInRange(40, 60)
```

**Step 2: Enumerate candidate mRNAs** (synonymous codon substitutions)
- There are ~3ⁿ candidates for a protein of length n
- Use constraint propagation to prune: if a codon choice creates a cryptic splice site, eliminate it immediately
- This is standard CSP solving (arc consistency, forward checking)

**Step 3: Type-check each surviving candidate**
- SpliceCorrect(HEK293T): Run the splicing NFST. If the output set is a singleton {P}, the candidate is splice-correct. Otherwise, MAYBE or NO.
- NoCrypticSplice: Scan against splice site PWMs. If no match exceeds the threshold, YES.
- CodonAdapted: Compute CAI. If ≥ 0.8, YES.
- GCInRange: Count. If in [40%, 60%], YES.

**Step 4: Output**
- If a candidate passes all checks → **YES**. Deliver the verified design with a formal guarantee.
- If no candidate passes → **NO**. Report the conflicting constraints.
- If some candidates pass some checks but not others → **MAYBE**. Report the uncertain regions and suggest constraint relaxations.

No probabilities. No training data beyond known rules. No calibration. Deterministic by construction.

---

## What This Enables That Probability Cannot

1. **Formal verification of gene therapy constructs**: Before injecting DNA into a patient, you want a *proof* that the mRNA will splice correctly — not a probability estimate.

2. **Compositional guarantees for genetic circuits**: When combining N circuit components, you need the *intersection* of their correctness properties, not the product of N probabilities (which decays exponentially).

3. **Regulatory compliance**: The FDA does not approve gene therapies based on probability estimates. It requires evidence of safety. A formal guarantee (type-checked, abstractly interpreted) is a stronger form of evidence than a calibrated probability.

4. **Zero-shot gene design for novel organisms**: Probabilistic models need training data from the target organism. Abstract interpretation needs only the known rules (splice site consensus, codon usage table), which can be obtained from a single genome sequence — no experimental data required.

5. **Proof-carrying gene designs**: The compiler can produce a **certificate** — a machine-checkable proof that the mRNA satisfies all stated properties. This is analogous to proof-carrying code in computer security. No one has done this for gene design.

---

## The Honest Limitation

The deterministic framework cannot answer questions that require quantitative predictions:
- "What fraction of transcripts will include exon 5?" → Needs probability
- "What is the expected protein yield?" → Needs probability
- "Which of two designs is better?" → Needs a ranking, which implies quantification

For these questions, probability is unavoidable. But for the questions that matter most in gene design — "Is this guaranteed to work?", "Is this guaranteed to be safe?", "Does this satisfy all constraints?" — determinism is not only possible but preferable.

The framework should be honest about this boundary. It provides deterministic guarantees for the questions it can answer, and it explicitly declines to answer questions that require probability — rather than giving a probabilistic answer that might be wrong.

---

## Summary

| Method | Key Idea | Deterministic Answer | Best For |
|--------|----------|---------------------|----------|
| Abstract interpretation | Sound over-approximation of all possible outcomes | Set of outcomes / three-valued (YES/NO/MAYBE) | Splicing analysis, mutation landscape exploration |
| Type systems | Static verification of correctness properties | Well-typed or type error | Gene design, safety-critical applications |
| Constraint satisfaction | Hard constraints, no soft scores | Satisfied or infeasible (with conflict set) | Codon optimization, gene synthesis |
| Non-deterministic automata | Set-valued functions, not distributions | Set of possible outputs | Pipeline composition |
| Game semantics | Worst-case analysis | Guaranteed for all contexts in class C | Safety-critical design |
| Three-valued logic | YES/NO/MAYBE instead of probabilities | Composable, no calibration needed | All of the above |

**The unifying principle**: Replace the question "What is the probability?" with "What is guaranteed?" The biology remains non-deterministic, but the framework's answers are deterministic, sound, and composable — the three properties that make compiler technology so powerful.
