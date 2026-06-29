# ADR-0011: GT-Free Codon Prioritization in Cryptic Splice Elimination

## Status: Accepted

## Date

2026-03-05

## Context

The cryptic splice elimination step was failing on ~80% of genes because it did not
prioritize GT-free codons for amino acids that have them. Valine is the ONLY amino
acid whose ALL codons contain GT. For C, G, R, S, GT-free alternatives exist and
should be tried first (guaranteed to eliminate the GT dinucleotide).

The previous implementation of the cryptic splice elimination step attempted codon swaps without considering
whether GT-free alternatives existed for a given amino acid. This meant that the
optimizer would often cycle through all synonymous codons for Cysteine (TGC, TGT),
both of which contain GT, when it could have simply chosen a GT-free codon — except
Cysteine has no GT-free codons. However, for amino acids like Alanine (GCA, GCC, GCG,
GCT), Glycine (GGA, GGC, GGG, GGT), Arginine (CGT-containing codons vs. AGA/AGG),
and Serine (AGT-containing codons vs. AGC/TCN), GT-free alternatives exist and
provide a guaranteed path to eliminating the GT dinucleotide.

The key biological insight is that the GT dinucleotide is the core recognition
sequence for splice donors. Any codon containing GT is a potential cryptic splice
donor. For amino acids where all synonymous codons contain GT (only Valine), no
codon-level fix is possible — this requires mutagenesis (V→I substitution). For all
other amino acids with GT-containing codons, a codon swap to a GT-free alternative
is guaranteed to eliminate the cryptic donor.

## Decision

Rewrite the cryptic splice elimination step with a 3-strategy approach:

1. **GT-free codon swap** (highest priority): For amino acids with GT-free alternatives,
   swap to the highest-CAI GT-free codon. This GUARANTEES elimination of the GT
   dinucleotide from that codon position. The swap is selected to maximize CAI among
   all GT-free synonymous codons, so the impact on codon adaptation is minimal.

2. **Context disruption** (for Valine): Since all V codons contain GT (GTT, GTC, GTA,
   GTG), no codon swap can eliminate the GT. Two sub-strategies are employed:
   - Select the V codon that produces the lowest MaxEntScan score when the GT is
     embedded in its 9-mer context (the surrounding nucleotides modulate splice
     site strength).
   - Disrupt the 9-mer context by swapping neighboring codons to break the extended
     consensus sequence that makes a GT look like a functional donor.

3. **Accept unrepairable**: Some Valine positions will remain above threshold even
   after context disruption. These are flagged for type-directed mutagenesis
   (V→I substitution, BLOSUM62=+3). The type system's derivation explicitly records
   these positions as `gt_mandatory=True`, distinguishing them from positions where
   the optimizer failed to use available GT-free codons.

The implementation adds a `_find_gt_free_codons(aa)` helper that returns the set
of synonymous codons for an amino acid that do not contain the GT dinucleotide, and
an `_is_unavoidable_gt_aa(aa)` function that returns True only for Valine (the sole amino
acid where all codons contain GT).

## Consequences

- **Positive**: NoCrypticSplice pass rate dramatically improves (from ~20% to expected
  ~60%+). The GT-free codon swap strategy provides a deterministic, guaranteed fix
  for the majority of cryptic splice positions.
- **Positive**: CAI impact is minimal — GT-free codons are often high-CAI themselves.
  For example, GCC (Alanine, GT-free) is the highest-CAI Alanine codon in human.
- **Positive**: The optimizer now properly distinguishes between "fixable" positions
  (where GT-free codons exist) and "unrepairable" positions (Valine only). This
  classification is propagated through the type system's derivation information,
  enabling the mutagenesis engine to act only on truly necessary positions.
- **Negative**: Context disruption for Valine is best-effort and may not always
  reduce the MaxEntScan score below threshold. These positions require mutagenesis.
- **Negative**: The 3-strategy approach adds complexity to the cryptic splice elimination step, but the
  improvement in pass rate justifies it.

## References

- ADR-0009: Type-Directed Protein Mutagenesis (V→I substitution)
- ADR-0014: Predicate Checking Delegation (standalone strong donor detection)
- ADR-0013: Mutagenesis GT-Mandatory vs Optimizer Weakness Distinction
