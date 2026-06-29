# BioCompiler Case Studies — "Types Compose, Constraints Don't"

Three compelling case studies demonstrating why BioCompiler's type system
approach succeeds where checklist-based approaches (like DNA Chisel) fail.

---

## Case Study 1: "The Impossible Design" — Type System Detects Unresolvable Conflicts

### The Problem

Consider a gene design where position 14 requires a Valine (V) amino acid.
All four Valine codons — GTT, GTC, GTA, GTG — start with the dinucleotide
"GT", which is the universal 5' splice donor signal. This creates a
fundamental conflict:

| Constraint | Requirement | Status |
|---|---|---|
| NoGTDinucleotide | No GT dinucleotides anywhere | Requires avoiding Valine codons |
| EnforceTranslation | Encode the target protein | Requires Valine at position 14 |
| ConservationScore | Preserve amino acid identity | Blocks V→I substitution (BLOSUM62=3) |

### Why Checklist Approaches Fail

DNA Chisel and similar tools treat constraints as an independent checklist.
They attempt to resolve each constraint in isolation:

1. **EnforceTranslation**: [OK] The sequence encodes the target protein
2. **AvoidPattern('GT')**: [FAIL] Valine codons contain GT
3. **Resolution attempt**: Try synonymous substitution...
   → All Val codons have GT. No synonymous fix exists.
4. **Try again?** Same result. And again. And again...

DNA Chisel would either:
- (a) Loop indefinitely trying to resolve the unresolvable
- (b) Hit a timeout and return PARTIAL with no explanation WHY
- (c) Silently relax the GT constraint without informing the user

### Why BioCompiler's Type System Wins

BioCompiler COMPOSES constraints via logical AND and detects when they
form an unsatisfiable conjunction:

```
NoGTDinucleotide ∧ ConservationScore(min=4) ∧ ValineAt(14) = UNSATISFIABLE
```

The type system returns **FAIL** with a certificate explaining:
- Which predicates failed (NoGTDinucleotide, NoCrypticSplice)
- The exact positions of violations
- WHY the conflict is unresolvable (all V codons contain GT, conservation blocks substitution)

This transforms an opaque failure into **actionable information**. The
biologist can now decide:
1. Accept the GT (if splice risk is tolerable for this context)
2. Allow mutagenesis V→I (if BLOSUM62=3 is acceptable)
3. Redesign the protein (if neither compromise is acceptable)

### Key Insight

**Types compose, constraints do not.** The checklist approach checks each
constraint independently and assumes PASS ∧ PASS → PASS. The type system
checks the CONJUNCTION of all constraints and detects when they are
mutually exclusive.

---

## Case Study 2: "Compositional Failure" — Where Constraints Don't Compose

### The Problem

Two gene fragments individually satisfy all constraints, but their
concatenation violates a constraint at the junction:

| Fragment | Sequence | Last/First Codon | Sau3AI-free? |
|---|---|---|---|
| A | ...GAT | GAT (Asp, D) | [OK] PASS |
| B | CTT... | CTT (Leu, L) | [OK] PASS |
| A+B | ...GAT**C**TT... | Junction: GATC | [FAIL] **FAIL** |

The Sau3AI restriction site **GATC** is formed at the junction between
the last base of Fragment A (T from GAT) and the first base of Fragment
B (C from CTT). Neither fragment contains GATC internally — it only
appears at the boundary.

### Why Checklist Approaches Miss This

```
Checklist approach:
1. Check Fragment A against [no_restriction_sites] → PASS [OK]
2. Check Fragment B against [no_restriction_sites] → PASS [OK]
3. Concatenate A + B
4. ??? No re-check step ???
5. Result: GATC site exists but was never detected
```

The fundamental issue: checklist items are checked independently. There is
no compositional logic — "PASS ∧ PASS → PASS" is ASSUMED but NOT
GUARANTEED. Junction regions are invisible to per-fragment checking.

### How BioCompiler Catches It

BioCompiler's `find_cross_codon_restriction()` function scans for
restriction sites that span codon boundaries. When fragments are
composed, it detects the cross-boundary GATC site at position 15.

The fix is trivial: change Fragment A's last codon from GAT to GAC
(both encode Aspartic acid, D). The junction becomes GACCTT — no
Sau3AI site.

### Key Insight

**Composition requires re-verification.** You cannot just AND two
constraint-satisfaction results. The type system enforces this
automatically by checking the composed artifact, not just the individual
components.

This is directly analogous to type checking in programming languages:
`f: A → B` and `g: B → C` compose to `g ∘ f: A → C`, but you must
verify the intermediate type B matches. BioCompiler verifies that
the "junction type" between fragments is consistent.

---

## Case Study 3: "The Certificate Saves the Day" — Verification Catches Optimizer Bug

### The Problem

An optimizer produces a sequence that appears "good" based on per-codon
checks, but a subtle cross-codon bug was introduced during optimization:

| Sequence | Codon 1 | Codon 2 | Boundary | GT? |
|---|---|---|---|---|
| Good | AAA (K) | TTT (F) | A\|T | No |
| Buggy | AAG (K) | TTT (F) | **G\|T** | **Yes** |

The optimizer changed AAA→AAG thinking it was merely a synonymous
substitution for Lysine (K). But AAG ends with G, and the next codon
TTT starts with T, creating a cross-codon GT dinucleotide.

### Why Naive Optimizers Miss This

A naive per-codon check examines each codon in isolation:

```
Codon 0 (ATG, M): GT within codon? False
Codon 1 (AAG, K): GT within codon? False   ← No GT within AAG
Codon 2 (TTT, F): GT within codon? False
...
→ All clear. No GT found.
```

The cross-codon GT at position 5-6 (G from AAG + T from TTT) is
invisible to per-codon checking. The optimizer says "all clear" but
the sequence has a cryptic splice donor signal.

### How the Certificate Verifier Catches It

BioCompiler's certificate system provides **three layers of protection**:

1. **GENERATION**: Every optimization produces a certificate documenting
   what was verified and at what level (GOLD/SILVER/BRONZE).

2. **VERIFICATION**: An independent verifier RE-EVALUATES every predicate
   from scratch. It does NOT trust the optimizer's claims.

3. **GRADUATED LEVELS**:
   - **GOLD**: All constraints satisfied by synonymous optimization
   - **SILVER**: All constraints met, but some have unavoidable trade-offs
   - **BRONZE**: Some constraints could not be satisfied — manual review needed

In this case, the good sequence gets a **GOLD** certificate (all predicates
pass). The buggy sequence gets a **BRONZE** certificate (NoGTDinucleotide
fails due to an avoidable GT at position 5).

The verifier independently re-runs `check_no_avoidable_gt()` and
detects the cross-codon GT, even though the optimizer missed it.

### Key Insight

**Independent verification is essential.** The optimizer is both the
producer and the checker — a conflict of interest. The certificate
system separates these roles: the optimizer produces, the verifier
checks. You do not have to trust the optimizer — you can VERIFY its
output.

This is directly analogous to formal verification in software engineering:
the compiler produces code, but the type checker (or proof checker)
independently verifies correctness. BioCompiler applies the same
principle to gene design.

---

## Running the Case Studies

> **Note:** These case studies are thought experiments illustrating BioCompiler's design philosophy, not runnable scripts. No `case_studies.py` script exists in the codebase.

```bash
# These case studies are conceptual demonstrations.
# They cannot be run as scripts.
```

Analysis results, if produced interactively, would be saved to a location of your choosing (not `download/case_studies/`).
- `case_study_1.json` — The Impossible Design
- `case_study_2.json` — Compositional Failure
- `case_study_3.json` — The Certificate Saves the Day
- `all_case_studies.json` — Combined results

## Code Location

The case studies are conceptual demonstrations embedded in this document.
There is no separate `case_studies.py` script in the codebase.
