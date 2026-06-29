# ADR-0009: Type-Directed Protein Mutagenesis

## Status

Accepted

## Date

2026-05-31

## Context

The type system proves properties about DNA sequences, but some properties are *mathematically impossible* to satisfy at the codon level for certain amino acids. The canonical example is Valine (V): all four Valine codons (GTT, GTC, GTA, GTG) contain the GT dinucleotide, which is the splice donor consensus. When a Valine position falls in a context that creates a strong cryptic splice donor (MaxEntScan score above threshold), NO codon swap can fix it — the GT is mandatory.

This is not an optimizer failure; it is a fundamental mathematical impossibility at the codon level. The only solution is to change the amino acid itself — crossing the central dogma boundary from DNA design to protein design.

**Key insight**: The type system's derivation (the detailed chain of reasoning for each FAIL verdict) is not just documentation — it is a *repair directive*. When the type system proves a predicate unsatisfiable at the codon level, the derivation tells us exactly which amino acid positions are responsible and why, enabling targeted substitutions.

**Alternatives Considered:**

1. **Accept infeasibility** — Report UNSAT and let the user decide. Simple; no risk of changing the protein. But: provides no guidance; the biologist does not know what to change; most biologists would prefer a concrete proposal they can evaluate.

2. **Brute-force substitution search** — Try every possible amino acid substitution at every position and re-optimize. Thorough. But: O(20^n) search space; computationally infeasible; most substitutions are biochemically unreasonable.

3. **Type-directed mutagenesis with BLOSUM62 ranking** — When the type system proves a predicate unsatisfiable at the codon level, analyze *why* it is unsatisfiable (which amino acid positions have mandatory problematic dinucleotides), propose conservative amino acid substitutions ranked by BLOSUM62 score, apply the best substitution, and re-optimize in a loop until all predicates pass or no more substitutions are possible. This is a targeted, biologically informed search.

## Decision

Implement type-directed mutagenesis with BLOSUM62 ranking (Alternative 3). The architecture is:

```
Optimizer → Type System → [FAIL]
                           |
                    Mutagenesis Engine
                    (analyzes FAIL derivation)
                           |
                    Propose V→I (BLOSUM62=+3)
                           |
                    Modified Protein → Optimizer → Type System → [PASS]
```

The key design choices:

1. **BLOSUM62 for substitution scoring**: The BLOSUM62 matrix (Henikoff & Henikoff, 1992) scores substitutions based on observed frequencies in aligned protein families. Positive scores indicate substitutions observed more often than expected by chance — i.e., evolutionarily acceptable changes. V→I scores +3, meaning this substitution is very common in nature.

2. **Conservative threshold**: Only substitutions with BLOSUM62 >= -1 (default) are proposed. This ensures all proposed substitutions are at least neutrally acceptable from a structural standpoint, while allowing slightly conservative substitutions that can resolve constraints.

3. **GT-mandatory amino acid identification**: The engine pre-computes which amino acids have the GT (or AG) dinucleotide in ALL their codons. Only Valine (V) has GT in all codons, making it the primary target for donor-elimination mutagenesis.

4. **Confidence levels**: Each proposed substitution is labeled as high/medium/low confidence based on whether the substitution definitively resolves the issue (e.g., V→I eliminates GT entirely = high confidence) or only partially helps.

> **Note:** Confidence levels are described here as designed but are not yet implemented in the current codebase.

5. **Loop with re-optimization**: After applying substitutions, the full optimization pipeline is re-run, because substitutions may change which codons are optimal and which constraints are binding.

**Proof of concept**: The HBB (hemoglobin beta) gene contains 15 Valine positions. Without mutagenesis, NoCrypticSplice fails (max_donor=5.71). With V→I substitutions (BLOSUM62=+3 each), the protein achieves 99.3% identity, CAI drops only ~0.2% (from 0.94 to 0.93), and NoCrypticSplice passes.

## Consequences

- Positive: (1) Makes previously impossible designs feasible — the HBB proof of concept turns a 5/6 predicate failure into 7/8 predicates passing. (2) Conservative — only proposes substitutions with strong BLOSUM62 support. (3) Transparent — every substitution is documented with position, original AA, substitute AA, BLOSUM62 score, reason, and predicate addressed. (4) Integrated with certificates — substitutions are recorded in certificate provenance metadata. (5) The derivation-as-repair-directive paradigm extends the type system's power beyond verification to design.
- Negative: (1) Changes the protein — not a pure codon optimization tool anymore. (2) BLOSUM62 is a heuristic — a +3 substitution is likely conservative but not guaranteed to preserve function. (3) The iterative loop may propose many substitutions before finding a feasible design, potentially reducing protein identity below acceptable thresholds. (4) Only addresses splice-related impossibility currently; other types of codon-level impossibility would need new analysis functions. (5) The user must explicitly opt in to enable mutagenesis — this is NOT applied by default, preserving backward compatibility.
