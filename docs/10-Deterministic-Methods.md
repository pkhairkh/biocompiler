# DOC-10: Deterministic Methods for Non-Deterministic Biology

| Field    | Value            |
|----------|------------------|
| ID       | DOC-10           |
| Version  | 12.0.0      |
| Status   | Current      |
| Date     | 2026-06-07       |

---

## 1. The Problem

The compiler metaphor is powerful precisely because compilers provide deterministic guarantees. A type checker does not say "this program is probably type-safe." A memory-safety analysis does not say "there is a 97% chance this access is in bounds." Compilers prove things. They do not estimate them.

A probabilistic compiler is, therefore, a category error. It collapses the distinction between engineering and statistics — between knowing something will hold and believing it probably will. No one would trust a compiler that produces executables which "probably won't segfault." Yet this is exactly what a probabilistic gene design framework offers: sequences that "probably" splice correctly, "probably" fold well, "probably" express.

The original BioCompiler paper made an understandable but fundamental mistake: it tried to make biology deterministic by averaging over outcomes. If splicing is stochastic, the reasoning went, then we should model the distribution of splicing outcomes and optimize for the expected case. This turns the compiler into an optimizer — a tool that finds the best probabilistic compromise rather than one that guarantees correctness.

But there is a third option, one that neither forces biology to be deterministic nor abandons determinism in the framework. The key insight is this:

**"You don't need a deterministic system to get deterministic answers. You need to ask the right questions — questions whose answers are deterministic even when the system is not."**

Consider a simple analogy. A non-deterministic finite automaton (NFA) does not produce a single output for a given input — it may follow many possible paths. Yet the question "does this NFA accept string w?" has a perfectly deterministic answer: yes or no. The automaton is non-deterministic; the answer is not. The same principle applies to biology. A gene may produce many possible splicing outcomes, but the question "can this gene produce an out-of-frame isoform?" has a deterministic answer. The biology is non-deterministic; the answer is not.

This document describes six formal methods that exploit this principle. Each method takes a biological system that is inherently non-deterministic — subject to cellular context, stochastic molecular events, and environmental variation — and produces deterministic answers about that system. None of these methods rely on probability. All of them provide guarantees that hold regardless of which non-deterministic outcome actually occurs.

These six methods are not alternatives to be chosen among. They are complementary perspectives on the same fundamental move: changing the question from "what will happen?" to "what can we guarantee?" Together, they form the mathematical foundation of a BioCompiler that is a true compiler — a tool that proves properties of designs, not one that estimates their likely performance.

---

## 2. The Core Move: Change the Question

The shift from probabilistic to deterministic answers is not a change in the biology. The biology is what it is: messy, context-dependent, stochastic. The change is in what we ask about the biology.

Consider the following table, which rephrases common probabilistic questions as their deterministic equivalents:

| Probabilistic Question | Deterministic Equivalent |
|---|---|
| "What is P(mRNA → protein)?" | "Is protein p in the set of all possible outputs?" |
| "What is P(splice site functional)?" | "Guaranteed functional / non-functional / uncertain?" |
| "Find mRNA maximizing P(correct) × P(expression)" | "Find mRNA satisfying all hard constraints" |
| "Distribution of folding outcomes?" | "Range of possible folding energies?" |

Each row asks about the same biological phenomenon. The left column asks for a probability — a number that requires distributional assumptions, training data, and calibration. The right column asks for a guarantee — a statement that can be derived from first principles, known rules, and logical inference.

The right-column answers are deterministic. They do not vary with repeated sampling. They do not depend on model calibration. They are either true or false, provably so.

This is not merely a philosophical preference. It has concrete engineering consequences:

- **Composability**: Deterministic answers compose without modeling dependencies. If component A is guaranteed correct and component B is guaranteed correct, then the composition is guaranteed correct. Probabilistic answers require joint distributions to compose, which are almost never available.
- **Verification**: A deterministic guarantee can be independently verified. A probability estimate cannot — it requires replicating the same experimental conditions, which is impossible in a living cell.
- **Regulation**: The FDA does not approve drugs that "probably work." It requires evidence of safety and efficacy. Deterministic guarantees are evidence. Probabilities are predictions.
- **Debugging**: When a deterministic framework says "NO, this design violates constraint X," you know exactly what went wrong. When a probabilistic framework says "this design has a 3% failure rate," you do not know which 3% or why.

The rest of this document elaborates six methods that operationalize this core move. Each method provides a different way to extract deterministic answers from non-deterministic biological systems.

---

## 3. Method 1: Abstract Interpretation

### 3.1 What It Is

Abstract interpretation is the mathematical foundation of static program analysis, introduced by Cousot and Cousot in 1977. The core idea is simple but profound: instead of computing the exact behavior of a program (which may be undecidable or intractable), compute a **sound over-approximation** of its behavior.

"Sound" means that every behavior the real system can exhibit is captured by the abstraction. "Over-approximation" means the abstraction may include behaviors the real system cannot exhibit. The abstraction is conservative: it never misses a real behavior, but it may include spurious ones.

The classic example from compiler construction illustrates this beautifully. A compiler does not need to know the exact value of variable `x` at runtime to prove that `a[x]` is a safe array access. It only needs to know that `x ∈ [0, len(a)-1]` — that is, that every possible value of `x` falls within the valid index range. This interval is an abstract domain that soundly approximates the concrete values `x` can take. The analysis is fast, deterministic, and gives a hard guarantee: if the abstract interpreter says the access is safe, it is safe in every possible execution.

Abstract interpretation provides a **Galois connection** between the concrete domain (actual behaviors) and the abstract domain (approximated behaviors). This Galois connection ensures soundness by construction. The abstraction function α maps concrete values to their abstract representations; the concretization function γ maps abstract values back to the sets of concrete values they represent. Soundness is the invariant: for every concrete value c, c ∈ γ(α(c)).

### 3.2 Applied to Protein Synthesis

In the context of the BioCompiler, the "program" is an mRNA sequence, and the "execution" is the process of protein synthesis — transcription, splicing, translation, folding. Each step is non-deterministic: splicing may skip exons, ribosomes may stall, proteins may misfold. Abstract interpretation allows us to reason about all possible executions simultaneously by computing in an abstract domain.

The following table defines abstract domains relevant to protein synthesis:

| Abstract Domain | What It Captures | Example |
|---|---|---|
| Set of possible isoforms | All proteins that can result from any splicing of this mRNA | {isoform₁, isoform₂, isoform₃} — "This mRNA can produce one of these three proteins" |
| Interval of inclusion rates | Range of possible exon inclusion rates across cellular contexts | [0.85, 1.0] for exon 5 — "Exon 5 included in 85–100% of transcripts" |
| Three-valued classification | Definite inclusion/exclusion/uncertainty per exon | {INCLUDED, EXCLUDED, UNCERTAIN} — "Exon 5 definitely included; exon 7 uncertain" |
| Bounds on folding energy | Interval of possible free energies for any fold the protein can adopt | [-45.2, -38.7] kcal/mol — "Folding energy in this range" |

Each of these is a deterministic answer about a non-deterministic process. The set of possible isoforms is computed once; it does not change with repeated sampling. The interval of inclusion rates is derived from known binding affinities and concentration ranges; it is not estimated from data. The three-valued classification is a logical inference from splicing rules; it is not a statistical prediction.

Consider the isoform set more concretely. An mRNA with three alternatively spliced exons can, in principle, produce 2³ = 8 possible isoforms. But not all combinations are biologically possible — splicing enhancers and silencers create dependencies between exons. The abstract interpreter computes the set of **feasible** isoforms by propagating constraints through the splicing grammar. The result is a set: for example, {001, 011, 111} (where each bit indicates inclusion). This set is a deterministic characterization of the mRNA's splicing behavior. It does not tell you how likely each isoform is, but it tells you exactly which isoforms are possible and — crucially — which are impossible.

### 3.3 Why This Is Not Just Renaming Probability

A skeptic might object: "You're just computing intervals instead of point estimates. An interval [0.85, 1.0] for exon inclusion is really just a confidence interval by another name."

This objection misses the fundamental distinction. A confidence interval is a statistical construct: it says that, with some confidence level (say 95%), the true parameter lies in the interval. It depends on the sample size, the distributional assumptions, and the calibration of the statistical model. Change the confidence level, and the interval changes.

An abstract interpretation interval is a **logical** construct. It says that, in **every possible execution** of the system — not 95% of executions, not 99%, but **all of them** — the exon inclusion rate falls within [0.85, 1.0]. This is not a statement about likelihood; it is a statement about necessity. If the abstract interpreter says "exon 5 definitely included," it is included in **every** possible execution — not just with high probability, but **necessarily**.

This is fundamentally different from P(exon 5 included) = 0.97. The probabilistic statement says "probably yes." The abstract interpretation says "necessarily yes." The difference matters when you are designing a gene therapy construct that will be injected into a patient. "Probably safe" and "provably safe" are not the same thing.

The soundness guarantee also means that abstract interpretation is **correct by construction**. If the abstract interpreter says a property holds, it holds in every concrete execution. This is not an empirical claim that needs to be validated by experiment — it is a mathematical theorem that follows from the Galois connection between the concrete and abstract domains.

### 3.4 The Price

Abstract interpretation is conservative. Because it soundly over-approximates, it sometimes says "uncertain" when the true answer is "yes" or "no." An exon might be included in 99.99% of executions, but if there exists a pathological cellular context (extreme stress, unusual splicing factor concentrations) in which it is skipped, the abstract interpreter will classify it as UNCERTAIN.

This conservatism is not a flaw — it is the price of soundness. And in the context of gene design, it is an acceptable price for three reasons:

1. **Gene design wants guarantees.** The purpose of a gene design tool is not to find the best possible design, but to find a design that is guaranteed to work. A conservative tool that sometimes rejects good designs is preferable to an optimistic tool that sometimes approves bad ones.

2. **Safety-critical applications demand conservatism.** In gene therapy, a false negative (rejecting a safe design) is far less costly than a false positive (approving an unsafe design). Abstract interpretation is biased toward false negatives, which is exactly the right direction for safety.

3. **Compositional verification requires conservatism.** When composing multiple analyses (splicing, folding, expression), a sound analysis of each component yields a sound analysis of the composition. A probabilistic analysis of each component does **not** yield a valid probabilistic analysis of the composition without modeling the dependencies between components — dependencies that are typically unknown.

The conservatism of abstract interpretation is not arbitrary. It can be systematically reduced by refining the abstract domain. A three-valued classification (INCLUDED/EXCLUDED/UNCERTAIN) is coarse. An interval classification ([0.85, 1.0]) is finer. A set of possible isoforms is finer still. The choice of abstract domain is an engineering decision that trades precision for computational cost.

---

## 4. Method 2: Type Systems for Biological Correctness

### 4.1 What It Is

Robin Milner's famous slogan — "Well-typed programs don't go wrong" — captures the essence of type theory as a verification method. A type checker proves properties of programs without executing them, without knowing their inputs, and without making any probabilistic assumptions. The guarantee is absolute: if the type checker says the program is well-typed, then no execution of the program will ever produce a type error.

This is a deterministic guarantee about a non-deterministic system. The program's execution depends on its inputs, which are unknown at compile time. The program's behavior is therefore non-deterministic from the compiler's perspective. Yet the type checker provides a deterministic answer: well-typed or not.

The power of type systems comes from their compositionality. The type of a compound expression is determined by the types of its sub-expressions, without knowing their concrete values. This means type checking scales: the cost of checking a program is linear (or near-linear) in its size, regardless of the size of its input space.

### 4.2 Applied to mRNA Design

In the BioCompiler, an mRNA sequence is analogous to a program, and the cellular machinery that processes it (polymerase, spliceosome, ribosome, chaperones) is analogous to the runtime system. A "type error" in this analogy is a failure mode: a cryptic splice site that causes aberrant splicing, a frameshift that produces a truncated protein, a destabilizing motif that causes mRNA degradation.

The following table defines types for mRNA design:

| Type | Meaning | Checkable How? |
|---|---|---|
| `SpliceCorrect(C)` | Guaranteed to splice correctly in cell type C | Splice site grammar + splicing factor concentrations for C |
| `NoCrypticSplice` | No sequences matching splice site consensus | Scan against splice site position weight matrices (PWMs) |
| `CodonAdapted(O, θ)` | Codon Adaptation Index ≥ θ for organism O | Lookup codon usage table for O, compute CAI |
| `GCInRange(lo, hi)` | GC content in [lo%, hi%] | Count G/C nucleotides, divide by total length |
| `NoRestrictionSite(S)` | No restriction enzyme cut sites from set S | Scan against REBASE recognition sequences |
| `InFrame` | All open reading frames in correct reading frame | Grammar check: verify exon boundaries preserve reading frame |
| `NoInstabilityMotif` | No mRNA destabilizing motifs (e.g., AU-rich elements, premature stop codons) | Motif scan against known destabilizing sequences |

Each type is checked by a deterministic procedure. The type checker either passes or fails — there is no "probably passes." The checking procedure is specified precisely enough that two independent implementations would produce identical results for any input.

Consider `NoCrypticSplice` in more detail. A cryptic splice site is a sequence within an exon or intron that sufficiently resembles the splice site consensus (GT...AG for the major spliceosome) to be mistakenly recognized by the spliceosome. The type checker scans the entire mRNA sequence against position weight matrices derived from known splice sites. If the match score exceeds a threshold, the sequence fails the type check. The threshold is chosen conservatively — it may flag sequences that would never be recognized as splice sites in vivo — but the check is deterministic and sound.

Consider `SpliceCorrect(C)`. This is a more sophisticated type that depends on the cellular context C. To check this type, the compiler needs a model of splicing in cell type C — specifically, which splicing factors are expressed and at what concentrations. This model defines the "runtime environment" for the mRNA. The type checker then verifies that, given this environment, the mRNA's splice sites are correctly recognized and its splicing pattern produces the intended isoform. This is equivalent to a type-and-effect system in programming language theory, where the effect depends on the runtime environment.

### 4.3 Subtyping as Biological Refinement

Type systems support subtyping, which enables hierarchical verification. The subtyping relation `A <: B` means "every value of type A is also of type B" — or, in biological terms, "every mRNA that satisfies specification A also satisfies specification B."

For the `SpliceCorrect` type, a natural subtyping hierarchy emerges:

```
SpliceCorrect(HEK293T) <: SpliceCorrect(HumanCell) <: SpliceCorrect(MammalianCell)
```

An mRNA that is guaranteed to splice correctly in HEK293T cells (a specific human cell line) is also guaranteed to splice correctly in any human cell (because HEK293T expresses a subset of the splicing factors found in all human cells). And an mRNA guaranteed to splice correctly in any human cell is also guaranteed to splice correctly in any mammalian cell (because human splicing machinery is a subset of mammalian splicing machinery).

Wait — this direction seems wrong. If HEK293T expresses a **subset** of human splicing factors, then an mRNA that splices correctly in HEK293T might encounter additional splicing factors in other human cells, potentially creating new splice variants. The subtyping direction depends on whether additional splicing factors can only **help** correct splicing (in which case the subtype direction is as stated) or can also **disrupt** it (in which case the direction reverses).

This tension reveals an important design principle: the subtyping hierarchy must be derived from the biological model, not assumed. If the model says that additional splicing factors can only promote intended splicing (a conservative assumption in some contexts), then the hierarchy is as stated. If the model says additional factors can create cryptic splicing, then a more refined type is needed — for example, `SpliceCorrect(HEK293T, factorSubset)` where the factor subset is explicit.

In any case, the subtyping relation is a deterministic, checkable relation. Given two types, the compiler can determine whether one is a subtype of the other. This enables hierarchical verification: verify against the most specific type first, and inherit the guarantees of all supertypes.

### 4.4 Type Inference for Gene Design

Type checking verifies that a given mRNA has a desired type. But in gene design, we often face the inverse problem: given a desired type (a specification), find an mRNA that has that type.

This is type inference — the process of deducing a type for an expression, or, in the inverse direction, constructing an expression that has a desired type. In programming language theory, type inference is what allows programmers to write code without explicit type annotations; the compiler infers the types. In gene design, type inference is what allows designers to specify what they want (e.g., "an mRNA encoding protein P that splices correctly in HEK293T cells") and have the compiler find a sequence that satisfies the specification.

The key insight is that type inference is deterministic: either a well-typed solution exists or it doesn't. There is no "maybe a solution exists with probability 0.7." The type inference algorithm either finds a solution (and provides a witness — the mRNA sequence) or proves that no solution exists (and provides a reason — the conflicting constraints).

This reframes gene design as a type inference problem. Instead of searching for a sequence that maximizes a probabilistic objective, we search for a sequence that satisfies a deterministic specification. The search space is still enormous (4^n for an mRNA of length n), but constraint propagation and type-directed synthesis can prune it dramatically.

---

## 5. Method 3: Constraint Satisfaction, Not Optimization

### 5.1 The Move

The standard approach to gene design is optimization: find the mRNA that maximizes some objective function, typically a product of probabilities:

```
maximize  P(correct splicing) × P(high expression) × P(no immunogenicity)
```

This formulation has several problems. The probabilities are estimated from limited data and may be poorly calibrated. The product assumes independence between factors, which is almost never true. The optimum may lie on a boundary where small perturbations cause failure. And the result is a single point estimate with no guarantee that it will work in practice.

The deterministic alternative is constraint satisfaction:

```
FIND:    mRNA m
SUCH THAT:
         m translates to protein P
         m has no cryptic splice sites
         CAI(m) ≥ 0.8
         GC(m) ∈ [40%, 60%]
         m has no restriction sites from set S
         m has no instability motifs
```

This is a constraint satisfaction problem (CSP). The answer is deterministic: either a feasible mRNA exists (and the solver finds one) or no feasible mRNA exists (and the solver reports that the constraints are inconsistent, along with a minimal unsatisfiable subset).

### 5.2 Comparison Table

The following table contrasts the two approaches across several dimensions:

| Dimension | Probabilistic Optimization | Constraint Satisfaction |
|---|---|---|
| Goal | Maximize expected quality | Guarantee all constraints met |
| Solution quality | A score (continuous) | Binary: satisfies or doesn't |
| Guarantees | None (only expected performance) | Hard guarantees by construction |
| Failure mode | May produce designs that fail with low probability | Cannot produce constraint-violating designs |
| When no solution exists | Returns "best" infeasible design | Reports infeasibility with explanation |
| Composition | Requires joint distributions | Constraints compose by conjunction |
| Calibration | Requires training data | Requires only known rules |
| Interpretability | "Why this design?" — "It has the highest score" | "Why this design?" — "It satisfies all constraints" |
| Robustness | Optimal design may be fragile | Any feasible design meets all constraints |

The constraint satisfaction approach is not inherently better in all situations. It is better when guarantees matter — when the cost of failure is high, when the system must work in diverse conditions, or when regulatory approval requires evidence of correctness rather than predictions of performance.

### 5.3 When Constraints Are Too Tight

What happens when the CSP is infeasible — when no mRNA satisfies all the constraints? This is not a failure of the method; it is a feature. Infeasibility is a diagnostic signal. It means the specification is overconstrained, and the designer needs to relax one or more constraints.

Modern CSP solvers provide a powerful debugging tool: the **minimal unsatisfiable subset (MUS)**. An MUS is a minimal set of constraints that, taken together, are inconsistent. Removing any single constraint from the MUS restores feasibility. The MUS tells the designer exactly which constraints are in conflict, enabling targeted relaxation.

For example, suppose the solver reports that the following constraints form an MUS:

1. CAI ≥ 0.9
2. GC ∈ [40%, 60%]
3. No cryptic splice sites
4. Encodes protein P

This means no mRNA encoding protein P can simultaneously have CAI ≥ 0.9, GC content in the specified range, and no cryptic splice sites. The designer can then make an informed decision: relax the CAI threshold to 0.8, widen the GC range, or accept a higher risk of cryptic splicing. Each relaxation is a conscious engineering trade-off, not an opaque compromise buried in an optimization objective.

This process — constraint relaxation guided by MUS analysis — is deterministic debugging. It replaces the probabilistic question "which constraint should I relax to maximize expected quality?" with the deterministic question "which constraint must I relax to make the problem feasible?" The latter has a definitive answer; the former does not.

---

## 6. Method 4: Non-Deterministic Automata (Without Probabilities)

### 6.1 The Move

A non-deterministic finite-state transducer (NDFST) maps each input to a **set** of possible outputs, not a distribution over outputs. This is a crucial distinction. Non-determinism means multiple outcomes are possible; it does not mean those outcomes have probabilities.

The set of possible outputs is computed deterministically. Given an input string and an NDFST, the set of outputs is uniquely determined — there is no randomness in the computation. The non-determinism is in the system being modeled (the biological process), not in the framework that models it.

This is precisely the situation in biology. Splicing is non-deterministic — the same pre-mRNA can produce different isoforms in different cells or even in the same cell. But this non-determinism is not probabilistic in any useful sense: the "probabilities" depend on cellular context, which varies across cells and over time. What we can determine, however, is the **set** of possible isoforms. This set is a well-defined, deterministically computable property of the mRNA sequence and the splicing machinery.

### 6.2 Applied to Splicing

We can model splicing as an NDFST:

- **Input**: pre-mRNA sequence
- **States**: positions in the pre-mRNA, splicing decisions made so far
- **Transitions**: at each potential splice site, either splice or don't splice (subject to biological constraints — splice sites must be recognized, exons must have minimum length, etc.)
- **Output**: the set of all possible splice isoforms

The NDFST computes the set of all isoforms that are biologically possible given the mRNA sequence and the splicing machinery. No probabilities are assigned to individual isoforms. The output is a set: {isoform₁, isoform₂, isoform₃, ...}.

This set is a deterministic answer to a non-deterministic process. It tells you exactly which outcomes are possible and which are impossible. It does not tell you how likely each outcome is, but for many design questions, this is exactly the right information:

- "Can this mRNA produce a truncated protein?" — Check whether any isoform in the set is truncated.
- "Are all possible isoforms in-frame?" — Check whether every isoform in the set preserves the reading frame.
- "Is the intended isoform the only possible one?" — Check whether the set is a singleton.

Each of these questions has a deterministic answer, derived from the set of possible isoforms.

### 6.3 Composition Without Probability

A key advantage of NDFSTs is that they compose without probability. If splicing is modeled as NDFST₁ and translation is modeled as NDFST₂, then the composition splicing-then-translation is NDFST₂ ∘ NDFST₁. The output set of the composition is the union of NDFST₂'s output sets over all inputs in NDFST₁'s output set:

```
Output(NDFST₂ ∘ NDFST₁, input) = ⋃{s ∈ Output(NDFST₁, input)} Output(NDFST₂, s)
```

This is a deterministic set computation. No joint distributions are needed. No independence assumptions are required. The composition is sound by construction: every output of the composed transducer is reachable by some execution of the pipeline, and every reachable output is in the set.

This is in stark contrast to probabilistic composition, where computing the output distribution of a pipeline requires the joint distribution of all intermediate variables — information that is typically unavailable and, even if available, would require modeling complex dependencies.

### 6.4 When Sets Get Large

A practical concern with set-based representations is that the sets can become very large. An mRNA with n alternatively spliced exons can produce up to 2^n possible isoforms. For large n, explicitly enumerating all isoforms is infeasible.

The solution is **symbolic representation**. Binary Decision Diagrams (BDDs) and their variants can compactly encode large sets of isoforms. A BDD represents a set as a directed acyclic graph where each node corresponds to a decision (splice or don't splice at a given site) and each path from root to leaf corresponds to an isoform. BDDs exploit regularity in the set — if many isoforms share common prefixes, the BDD merges them, resulting in a compact representation.

Symbolic set operations (union, intersection, membership testing) are efficient on BDDs. This means that NDFST composition, which requires computing unions of output sets, can be performed efficiently even when the underlying sets are exponentially large.

Other symbolic representations are also possible: ZDDs (zero-suppressed decision diagrams) for sparse sets, SAT/SMT solvers for constraint-based representations, and regular expressions for sets with repetitive structure. The choice of representation is an implementation detail; the key principle is that set-based computation is feasible even for large sets, provided the right symbolic representation is used.

---

## 7. Method 5: Game Semantics — Worst-Case Guarantees

### 7.1 The Move

Game semantics models computation as a game between two players: the system (which tries to satisfy a specification) and the environment (which tries to violate it). The system wins if it can satisfy the specification regardless of the environment's moves. The environment wins if it can force a violation.

In the biological context, the "system" is the gene design, and the "environment" is the cellular context — the concentrations of splicing factors, the temperature, the stress level, the availability of ribosomes, and all the other variables that affect gene expression. The cellular context is non-deterministic: it varies across cells and over time, and we cannot predict it precisely.

Game semantics replaces probabilistic reasoning about the environment with adversarial reasoning. Instead of asking "what is the probability that the cellular context will cause this mRNA to splice incorrectly?", we ask "is there **any** cellular context in class C that can cause this mRNA to splice incorrectly?" If the answer is no, then the mRNA is guaranteed to splice correctly in every context in C — not with high probability, but necessarily.

| Question | Probabilistic Answer | Game-Semantic Answer |
|---|---|---|
| Will this mRNA splice correctly? | "P = 0.97" | "Yes, even in worst-case context in class C" |
| Will this protein fold? | "P = 0.85" | "Yes, for all temperatures in [37°C, 39°C]" |
| Is this mRNA stable? | "Expected half-life = 6 hours" | "Half-life ≥ 2 hours in all contexts in C" |

The universal quantifier ("for all contexts") replaces the probability. The answer is deterministic: either the mRNA is correct in every context in the class, or there exists a context in which it fails.

### 7.2 Refining the Adversary Class

The conservativeness of game semantics depends on the adversary class — the set of contexts the environment can choose from. A broad adversary class (e.g., "all possible mammalian cells") yields very conservative guarantees: the design must work in every possible context, including pathological ones. A narrow adversary class (e.g., "HEK293T cells under standard culture conditions") yields tighter guarantees but is applicable to fewer situations.

The adversary class can be refined to make guarantees more precise:

1. **Broadest class**: All possible cellular contexts. This yields the most conservative guarantees and is appropriate for gene therapies that must work in diverse patient populations.

2. **Organism class**: All cells of a given organism. This is less conservative because it restricts the range of splicing factor concentrations to those found in the organism.

3. **Tissue class**: All cells of a given tissue type. Even less conservative, because cells within a tissue share a more similar expression profile.

4. **Cell line class**: A specific cell line under defined culture conditions. The least conservative, because the cellular context is most constrained.

More precise context specification leads to less conservative guarantees. But the guarantees are always deterministic. The refinement does not introduce probability — it simply narrows the set of contexts over which the universal quantifier ranges.

This approach naturally supports a hierarchy of guarantees: a design that is guaranteed correct for all mammalian cells is also guaranteed correct for all human cells, which is also guaranteed correct for HEK293T cells. The hierarchy goes in the direction of increasing strength: broader adversary classes imply stronger guarantees.

---

## 8. Method 6: Three-Valued Logic

### 8.1 The Move

Classical logic is two-valued: every proposition is either TRUE or FALSE. But in the presence of non-determinism and incomplete information, two values are not enough. Some propositions are genuinely uncertain — we cannot determine their truth value from the available information.

Three-valued logic adds a third value: MAYBE (also called UNKNOWN or UNDEFINED). The three values correspond to three epistemic states:

| Answer | Meaning | Action |
|---|---|---|
| YES | Guaranteed to hold in every possible execution | Proceed with confidence |
| NO | Guaranteed not to hold in any possible execution | Reject the design |
| MAYBE | Might hold in some executions and not in others | Investigate further; refine the analysis |

Three-valued logic is not a weakening of classical logic — it is an honest representation of the limits of our knowledge. When we do not have enough information to determine whether a property holds, we say MAYBE rather than guessing YES or NO.

The logical connectives of three-valued logic are defined to preserve soundness:

- **Conjunction**: YES ∧ YES = YES; YES ∧ MAYBE = MAYBE; YES ∧ NO = NO; MAYBE ∧ MAYBE = MAYBE; MAYBE ∧ NO = NO; NO ∧ NO = NO.
- **Disjunction**: YES ∨ YES = YES; YES ∨ MAYBE = YES; YES ∨ NO = YES; MAYBE ∨ MAYBE = MAYBE; MAYBE ∨ NO = MAYBE; NO ∨ NO = NO.
- **Negation**: ¬YES = NO; ¬NO = YES; ¬MAYBE = MAYBE.

These definitions ensure that YES and NO propagate soundly. If any conjunct is NO, the conjunction is NO. If any disjunct is YES, the disjunction is YES. MAYBE arises only when the information is genuinely insufficient.

### 8.2 Why Three Values Beat Probabilities

Three-valued logic has several practical advantages over probabilistic reasoning:

1. **Composability without dependence modeling.** Probabilistic conjunction requires P(A ∧ B) = P(A) × P(B | A), which depends on the conditional probability P(B | A). If A and B are independent, this simplifies to P(A) × P(B), but independence is almost never guaranteed in biology. Three-valued conjunction requires no dependence information: YES ∧ YES = YES, regardless of any dependence between the conjuncts.

2. **No calibration needed.** Probabilities must be calibrated — P = 0.97 must mean that the event occurs 97% of the time, or the probability is meaningless. Calibration requires extensive experimental data and careful statistical modeling. Three-valued logic requires no calibration: YES means "guaranteed," NO means "impossible," and MAYBE means "unknown." No data is needed to interpret these values.

3. **No training data needed.** Probabilistic models are typically learned from training data. If the training data is biased, incomplete, or unrepresentative, the model's predictions will be wrong. Three-valued logic requires only **known rules** — biological constraints that have been established by the literature. It makes no claims about phenomena that are not covered by the rules, instead returning MAYBE.

4. **Interpretable.** A probability of 0.97 is difficult to interpret in context. Does it mean the design will fail 3% of the time? In 3% of cells? In 3% of patients? Under 3% of conditions? A YES/NO/MAYBE answer is immediately interpretable: the design is guaranteed to work, guaranteed to fail, or its behavior is uncertain and needs further investigation.

The primary limitation of three-valued logic is its coarseness. It cannot distinguish between "almost certainly yes" and "barely maybe." For quantitative questions — "what fraction of cells will express this gene?" — probabilities are necessary. But for the yes/no questions that arise in design verification — "will this gene ever produce a truncated protein?" — three-valued logic is both sufficient and more honest than probability.

---

## 9. The Unified Framework: Deterministic by Construction

The six methods described in this document are not independent alternatives. They are different facets of a single framework for extracting deterministic answers from non-deterministic biology. The following diagram shows how they combine into a unified pipeline:

```
                              ┌──────────────────────────────────────┐
                              │                                      │
    Input:                    │   Deterministic Analysis Engine       │
                              │                                      │
    mRNA sequence             │   ┌─────────────────────┐            │
    +                         │   │ Abstract Interpreter │            │
    specification             │   └──────────┬──────────┘            │
         │                    │              │                       │
         │                    │   ┌──────────┴──────────┐            │
         └────────────────────┤   │ Type Checker        │            │
                              │   └──────────┬──────────┘            │
                              │              │                       │
                              │   ┌──────────┴──────────┐            │
                              │   │ CSP Solver           │            │
                              │   └──────────┬──────────┘            │
                              │              │                       │
                              │   ┌──────────┴──────────┐            │
                              │   │ NDFST Engine         │            │
                              │   └──────────┬──────────┘            │
                              │              │                       │
                              │   ┌──────────┴──────────┐            │
                              │   │ Game Analyzer        │            │
                              │   └──────────┬──────────┘            │
                              │              │                       │
                              │   ┌──────────┴──────────┐            │
                              │   │ Three-Valued Evaluator│           │
                              │   └──────────┬──────────┘            │
                              │              │                       │
                              └──────────────┼───────────────────────┘
                                             │
                                             ▼
                              ┌──────────────────────────────────────┐
                              │                                      │
                              │   Output:                            │
                              │                                      │
                              │   If YES:  verified design +         │
                              │            proof certificate         │
                              │                                      │
                              │   If NO:   counterexample +          │
                              │            conflict set / MUS        │
                              │                                      │
                              │   If MAYBE: uncertain regions +      │
                              │             suggested refinements     │
                              │                                      │
                              └──────────────────────────────────────┘
```

The pipeline works as follows:

1. The user provides an mRNA sequence and a specification (the desired properties the mRNA should satisfy).

2. The deterministic analysis engine applies all six methods, each contributing a different perspective:
   - **Abstract interpretation** computes over-approximations of the mRNA's behavior (possible isoforms, folding energy ranges, etc.).
   - **Type checking** verifies that the mRNA satisfies the specified types (no cryptic splice sites, correct reading frame, etc.).
   - **CSP solving** checks whether the design constraints are simultaneously satisfiable and, if not, identifies the minimal unsatisfiable subset.
   - **NDFST analysis** computes the set of all possible outputs of the mRNA processing pipeline (splicing, translation, folding).
   - **Game analysis** checks whether the mRNA satisfies the specification in the worst-case cellular context.
   - **Three-valued evaluation** combines the results of all analyses into a single YES/NO/MAYBE verdict.

3. The output is deterministic:
   - **YES**: The mRNA is guaranteed to satisfy the specification. A proof certificate is provided, documenting which analyses passed and why.
   - **NO**: The mRNA is guaranteed to violate the specification. A counterexample is provided, showing a specific violation, along with the conflict set or MUS identifying the root cause.
   - **MAYBE**: The available information is insufficient to determine whether the mRNA satisfies the specification. The uncertain regions are identified, and suggestions for refining the analysis (e.g., narrowing the adversary class, using a finer abstract domain) are provided.

The computation is deterministic. The same input always produces the same output. The biology is non-deterministic — the mRNA may behave differently in different cells — but the framework models this non-determinism set-theoretically (as sets of possible outcomes) or logically (as universal quantification over contexts). No probabilities are needed.

---

## 10. Comparison: Three Approaches

The following table compares three approaches to gene design across thirteen dimensions:

| Dimension | Original Paper (Fake Determinism) | Revised Framework (Probability) | This Document (True Determinism) |
|---|---|---|---|
| **Model of biology** | Deterministic (averaged) | Stochastic (distributional) | Non-deterministic (set-theoretic) |
| **Framework's answers** | Deterministic but wrong | Probabilistic | Deterministic and correct |
| **Guarantees** | False precision | Expected values, no hard guarantees | Hard guarantees (sound by construction) |
| **Composition** | Trivial (assumes independence) | Requires joint distributions | Set union / logical conjunction |
| **Validation** | Matches averaged data | Matches distribution | Subsumes all possible outcomes |
| **When wrong** | Systematically misleading | Poorly calibrated | Over-conservative (never under-conservative) |
| **When uncertain** | Reports false precision | Reports wide confidence intervals | Reports MAYBE with diagnostic info |
| **Training data** | Required for averaging | Required for distribution estimation | Not required (only known rules) |
| **Novelty** | Cannot generalize beyond training distribution | Generalizes via distributional assumptions | Generalizes to any context within the adversary class |
| **Interpretability** | "This design is optimal" (misleading) | "This design has the highest expected score" | "This design satisfies all constraints" or "No such design exists" |
| **Debugging** | Difficult (opaque optimization) | Difficult (many contributing factors) | Structured (MUS, counterexamples, uncertain regions) |
| **Regulatory applicability** | Weak (averages hide failures) | Moderate (probabilities require justification) | Strong (guarantees are evidence) |
| **Philosophical status** | Category error (deterministic model of non-deterministic system) | Honest but incomplete (probabilities are not guarantees) | Correct (deterministic answers about non-deterministic system) |

The three approaches are not mutually exclusive. A practical gene design tool may use all three: true determinism for safety-critical properties, probability for quantitative predictions, and averaged models for quick initial screening. The key is to use each approach for what it is good at and to be honest about its limitations.

The original paper's approach — fake determinism — is the least defensible. Averaging over non-deterministic outcomes and treating the average as the outcome is not just an approximation; it is a misrepresentation. It hides the variability that is the most important feature of biological systems.

The probabilistic framework is honest about variability but cannot provide guarantees. It is appropriate for quantitative questions ("what yield can I expect?") but inappropriate for safety questions ("will this gene therapy ever produce a toxic protein?").

The true deterministic approach described in this document is honest about variability and provides guarantees. It is appropriate for safety questions and design verification, but it cannot answer quantitative questions. For those, probability is unavoidable.

---

## 11. Concrete Example: Splicing-Aware Gene Design, Deterministically

This section walks through a concrete example of splicing-aware gene design using the deterministic framework. The goal is to design an mRNA encoding a therapeutic protein (e.g., Factor VIII for hemophilia A) that is guaranteed to splice correctly in human liver cells.

### Step 1: Define the Type

The specification is expressed as a type:

```
mRNA m : SpliceCorrect(HumanHepatocyte) ∧ NoCrypticSplice ∧ CodonAdapted(Human, 0.8) ∧ GCInRange(40, 60) ∧ NoRestrictionSite({EcoRI, BamHI, XhoI}) ∧ InFrame
```

This type specifies that the mRNA must splice correctly in human hepatocytes, have no cryptic splice sites, have a codon adaptation index of at least 0.8 for humans, have GC content between 40% and 60%, contain no restriction sites for the specified enzymes, and preserve the correct reading frame.

### Step 2: Enumerate Candidates with Constraint Propagation

The codon space for a protein of typical therapeutic size (~1500 amino acids for Factor VIII) is astronomically large (roughly 3^1500, given ~3 synonymous codons per amino acid on average). Exhaustive enumeration is impossible, but constraint propagation dramatically reduces the search space:

1. **CodonAdapted(Human, 0.8)**: Eliminate low-CAI codons. For each amino acid, keep only codons with relative adaptiveness above a threshold. This typically eliminates 1–2 codons per amino acid, reducing the space by a factor of ~2^1500.

2. **NoCrypticSplice**: Scan the remaining codon choices for sequences matching the splice site consensus (GT...AG with sufficient PWM score). Where a cryptic splice site is found, eliminate the codon choices that create it. This is a local constraint that eliminates specific codons at specific positions.

3. **GCInRange(40, 60)**: Track the running GC count. If adding a GC-rich codon would make the GC content exceed 60%, eliminate that codon. If avoiding GC-rich codons would make the GC content fall below 40%, eliminate the low-GC alternatives. This is a global constraint that interacts with the other choices.

4. **NoRestrictionSite({EcoRI, BamHI, XhoI})**: Scan for restriction enzyme recognition sequences (GAATTC for EcoRI, GGATCC for BamHI, CTCGAG for XhoI). Eliminate codon choices that create these sequences.

5. **InFrame**: Ensure exon boundaries preserve the reading frame. This is a structural constraint that depends on the exon structure of the gene.

After constraint propagation, the remaining search space is typically orders of magnitude smaller than the original. In many cases, constraint propagation alone is sufficient to find a solution — no search is needed.

### Step 3: Type-Check Survivors

For each candidate mRNA that survives constraint propagation, run the full type checker:

1. **SpliceCorrect(HumanHepatocyte)**: Model splicing in human hepatocytes using the NDFST splicing model. Compute the set of all possible splice isoforms. Verify that the set contains only the intended isoform (or, more conservatively, that every isoform in the set produces a functional protein).

2. **NoCrypticSplice**: Rescan with a more sensitive splice site model (e.g., including auxiliary sequences like exonic splicing enhancers and intronic splicing silencers) to catch cryptic sites that the simple PWM scan missed.

3. **CodonAdapted(Human, 0.8)**: Compute the CAI exactly and verify it meets the threshold.

4. **GCInRange(40, 60)**: Compute the exact GC content and verify it is in range.

5. **NoRestrictionSite({EcoRI, BamHI, XhoI})**: Verify no restriction sites exist.

6. **InFrame**: Verify all reading frames are correct.

Each check is deterministic. Pass or fail. No probabilities.

### Step 4: Output YES/NO/MAYBE with Evidence

- **If YES**: The mRNA satisfies all constraints. Output the verified design along with a proof certificate documenting each check that passed, the abstract domains used, the adversary class considered, and the reasoning steps.

- **If NO**: The mRNA violates one or more constraints. Output a counterexample (the specific violation) and a conflict set / MUS (the minimal set of constraints that are mutually inconsistent). The designer can then relax the appropriate constraints and try again.

- **If MAYBE**: Some constraints could not be definitively verified. Output the uncertain regions (e.g., "Exon 7 splicing is uncertain in hepatocytes under stress conditions") and suggested refinements (e.g., "Narrow the adversary class to hepatocytes under standard conditions" or "Use a finer splicing model that includes stress-responsive splicing factors").

The entire process is deterministic. The same mRNA and specification always produce the same verdict. The biology is non-deterministic (the mRNA may splice differently in different cells), but the framework accounts for this by considering all possible outcomes and providing guarantees that hold for every outcome.

---

## 12. What This Enables That Probability Cannot

The deterministic framework enables several capabilities that are fundamentally impossible with a probabilistic approach:

### 1. Formal Verification of Gene Therapy Constructs

Gene therapy constructs are safety-critical medical devices. They must be verified to be safe before being administered to patients. A probabilistic guarantee ("97% chance of correct splicing") is not sufficient for regulatory approval — the FDA requires evidence that the construct will not produce harmful proteins, not merely that it probably won't. The deterministic framework provides exactly this kind of evidence: a proof that the construct satisfies safety properties in all possible cellular contexts within a specified class.

### 2. Compositional Guarantees for Genetic Circuits

Genetic circuits — engineered systems of interacting genes — are composed of multiple genetic parts, each of which may be non-deterministic. Verifying the correctness of the composition requires reasoning about all possible interactions between parts. Probabilistic composition requires modeling the joint distribution of part behaviors, which is intractable for realistic circuits. Deterministic composition — set union for possible outcomes, logical conjunction for properties — is tractable and sound.

### 3. Regulatory Compliance

Regulatory agencies (FDA, EMA, PMDA) require evidence of safety and efficacy for gene therapy products. Probabilistic predictions are not evidence — they are estimates that may be wrong. Deterministic guarantees are evidence: they are logical consequences of established biological rules and can be independently verified by the regulator's own computational tools. The deterministic framework can produce machine-checkable proof certificates that regulators can verify without trusting the BioCompiler.

### 4. Zero-Shot Gene Design for Novel Organisms

Probabilistic models require training data from the target organism. If no data is available (e.g., for a newly sequenced organism), the model's predictions are unreliable. The deterministic framework requires only **known rules** — conserved biological mechanisms that are shared across organisms. For example, the splice site consensus (GT...AG) is conserved across all eukaryotes. A design that is guaranteed to have no cryptic splice sites in any organism that uses the major spliceosome is a zero-shot guarantee — it requires no organism-specific data.

### 5. Proof-Carrying Gene Designs

A proof-carrying gene design is an mRNA sequence accompanied by a machine-checkable certificate that proves the sequence satisfies its specification. The certificate can be independently verified by any party — the designer, the regulator, the manufacturer, the clinician — without trusting the tool that produced it. This is analogous to proof-carrying code in computer science, where a program is accompanied by a proof of memory safety that can be checked by the runtime system.

Proof-carrying gene designs are only possible with a deterministic framework. A probabilistic "proof" would be a statistical argument that depends on the validity of the probabilistic model, the quality of the training data, and the correctness of the calibration — all of which must be trusted. A deterministic proof depends only on the logical rules used by the verifier, which can be independently inspected and validated.

---

## 13. The Honest Limitation

The deterministic framework described in this document has a fundamental limitation: **it cannot answer quantitative questions.**

Questions like the following are inherently probabilistic and cannot be answered deterministically:

- "What fraction of cells will express this gene?"
- "What is the expected protein yield?"
- "Which of these two designs will produce more protein?"
- "What is the probability that this gene therapy will be effective?"

For these questions, probability is unavoidable. The deterministic framework should decline to answer them rather than giving misleading deterministic answers.

This is not a flaw in the framework — it is an honest acknowledgment of its scope. The framework provides guarantees, not predictions. It tells you what will always happen and what can never happen. It does not tell you what will usually happen or what is most likely to happen.

For many design questions, guarantees are sufficient. If a gene therapy construct is guaranteed to produce only the intended protein (no toxic byproducts), guaranteed to express in the target tissue (even in the worst case), and guaranteed to be stable (no degradation motifs), then the quantitative details of expression level can be determined experimentally. The guarantee eliminates the most important risks; the remaining quantitative optimization can be done empirically.

But for other questions — particularly those involving cost optimization, yield maximization, or comparative evaluation of alternative designs — quantitative predictions are essential. The deterministic framework should not pretend to answer these questions. Instead, it should be complemented by a probabilistic framework that handles quantitative prediction, with the understanding that probabilistic predictions are estimates, not guarantees.

The honest division of labor is:

- **Deterministic framework**: Safety, correctness, feasibility. "Will this design work?" → YES/NO/MAYBE.
- **Probabilistic framework**: Performance, yield, optimization. "How well will this design work?" → Point estimate + confidence interval.

Neither framework subsumes the other. Both are needed. But the deterministic framework should come first, because a design that is not guaranteed to be correct should not be optimized for performance.

---

## 14. Summary Table

| Method | Key Idea | Deterministic Answer | Best For |
|---|---|---|---|
| **Abstract Interpretation** | Compute sound over-approximations of behavior in abstract domains | Over-approximated properties (e.g., set of possible isoforms, interval of folding energies) | Reasoning about all possible outcomes simultaneously; proving safety properties; scalable analysis of large sequences |
| **Type Systems** | Classify sequences by the properties they guarantee; check types without execution | Type judgment: well-typed or not | Fast rejection of incorrect designs; compositional verification; hierarchical specification via subtyping |
| **Constraint Satisfaction** | Replace optimization with feasibility; find designs satisfying all hard constraints | Feasible/infeasible + witness/MUS | Design problems with multiple hard constraints; debugging overconstrained specifications; systematic exploration of the design space |
| **Non-Deterministic Automata** | Model biological processes as NDFSTs; compute sets of possible outputs without probabilities | Set of all possible outputs | Modeling splicing and other alternative processing; composing multi-step biological pipelines; checking properties of all possible outcomes |
| **Game Semantics** | Model cellular context as adversary; guarantee correctness in worst case | "Correct for all contexts in class C" | Safety-critical applications; gene therapy verification; designs that must work across diverse cellular environments |
| **Three-Valued Logic** | Replace probabilities with YES/NO/MAYBE; honest representation of uncertainty | YES / NO / MAYBE with diagnostic information | Combining results from multiple analyses; reporting verification results to users; identifying regions that need further investigation |

Each method addresses a different aspect of the fundamental challenge: extracting deterministic answers from non-deterministic biology. Together, they form a coherent framework that provides the kinds of guarantees that compilers provide for software — not by ignoring non-determinism, not by averaging over it, but by asking questions whose answers are deterministic even when the system is not.

The framework is conservative: it sometimes says MAYBE when the answer is really YES, and it sometimes rejects designs that would work in practice. But it is never wrong: when it says YES, the design is guaranteed to work; when it says NO, the design is guaranteed to fail. In the context of gene design — where failure can mean a toxic protein in a patient's cells — this conservatism is not just acceptable. It is essential.

---

*End of DOC-10*
