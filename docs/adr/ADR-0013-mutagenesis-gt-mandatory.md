# ADR-0013: Mutagenesis GT-Mandatory vs Optimizer Weakness Distinction

## Status: Accepted

## Date

2026-03-05

## Context

The mutagenesis engine was proposing amino acid substitutions for positions where
the OPTIMIZER should have fixed the problem (by choosing a GT-free codon). This
conflates two fundamentally different issues:

1. **GT-mandatory** (Valine): No codon swap can eliminate GT from any Valine codon
   (GTT, GTC, GTA, GTG all contain GT). Mutagenesis IS needed — the only way to
   remove the GT is to substitute Valine with a different amino acid whose codons
   don't contain GT (e.g., Isoleucine: ATT, ATC, ATA).

2. **Optimizer weakness**: The amino acid has GT-free codons but the optimizer
   didn't use them. For example, Alanine has GT-free codons (GCA, GCC, GCT) but
   the optimizer may have chosen GCG or GGT (if Glycine) which contain GT. This
   is a bug or limitation in the optimizer, not a fundamental impossibility.
   Mutagenesis is NOT the right tool — the optimizer should be fixed.

Before this distinction was made, the mutagenesis engine would propose substitutions
for both categories, leading to:
- Unnecessary protein modifications (the optimizer could have fixed it with a codon
  swap)
- Masked optimizer bugs (the mutagenesis "fix" hid the fact that the optimizer
  failed to use available GT-free codons)
- Reduced protein identity (unnecessary AA substitutions lower the similarity to
  the original protein)

The type system's `evaluate_no_cryptic_splice` function already produces rich
derivations including `responsible_codon_index`, `responsible_aa`,
`problematic_dinucleotide`, and `all_codons_for_aa`. This information is sufficient
to determine whether a failing position is GT-mandatory or an optimizer weakness.

## Decision

1. Add `is_gt_mandatory(aa)` function to the mutagenesis module. This function
   returns `True` if and only if ALL synonymous codons for the given amino acid
   contain the GT dinucleotide. Currently, this is True only for Valine (V).

2. Extend `find_unrepairable_cryptic_donors` to include a `gt_mandatory` field in
   each returned position. This field is computed by checking whether the amino
   acid at the cryptic donor position has any GT-free codons.

3. Only propose AA substitutions for GT-mandatory positions in
   `type_directed_mutagenesis`. Positions where `gt_mandatory=False` are skipped
   — these represent optimizer weaknesses that should be addressed by improving
   the optimizer (cryptic splice elimination step GT-free codon prioritization, ADR-0011), not by changing
   the protein.

4. Add `diagnose_optimizer_weakness()` function that identifies positions where
   the optimizer failed to use available GT-free codons. This function returns a
   list of `(position, current_codon, gt_free_alternatives)` tuples, providing
   actionable diagnostic information for optimizer debugging.

5. Extend the type system's derivation to include `gt_mandatory` and
   `gt_free_alternatives` fields. When `evaluate_no_cryptic_splice` returns a FAIL
   verdict, the derivation now includes:
   - `gt_mandatory: bool` — True if all codons for this AA contain GT
   - `gt_free_alternatives: list[str]` — The GT-free codons available (empty if
     GT-mandatory)

## Consequences

- **Positive**: Mutagenesis only proposes substitutions that are truly necessary.
  This reduces the number of AA substitutions and preserves more of the original
  protein sequence.
- **Positive**: Optimizer bugs are clearly identified rather than masked by
  mutagenesis. The `diagnose_optimizer_weakness()` function provides specific
  information about which positions the optimizer failed to fix and what GT-free
  alternatives were available.
- **Positive**: The V→I substitution (BLOSUM62=+3) remains the primary mutagenesis
  proposal, since Valine is the only GT-mandatory amino acid. This is a highly
  conservative substitution with well-characterized structural effects.
- **Negative**: Some genes that previously "passed" via unnecessary mutagenesis may
  now fail because the optimizer weakness is exposed rather than masked. This is
  the correct behavior — the optimizer should be fixed, not the protein.
- **Positive**: The type system's enhanced derivation provides richer diagnostic
  information, enabling both the mutagenesis engine and human reviewers to
  understand exactly why a position fails and what the options are.

## References

- ADR-0009: Type-Directed Protein Mutagenesis (BLOSUM62 substitution engine)
- ADR-0011: GT-Free Codon Prioritization in Cryptic Splice Elimination (fixes
  optimizer weaknesses for non-Valine amino acids)
- ADR-0014: Predicate Checking Delegation (rich derivations from type system)
