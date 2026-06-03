/-
BioCompiler v7.0.0 — Formal Verification
==========================================
Entry point for the Lean formal proof of certified gene optimization.

This library proves:
  1. The genetic code is a total, deterministic function from codons to amino acids
  2. 8 certification predicates are well-defined and compose correctly
  3. The dual-threshold NoCrypticSplice system is monotone
  4. NoGTDinucleotide subsumes NoCrypticSplice(FAIL)
  5. Mutagenesis preserves conservation when BLOSUM62 >= threshold
  6. GOLD certificates imply all 8 predicates are satisfied
  7. SILVER certificates imply all predicates passed (some via mutagenesis)
  8. BRONZE certificates imply at least one predicate is unsatisfied
-/

import BioCompiler.CodonTable
import BioCompiler.Predicates
import BioCompiler.Optimization
import BioCompiler.Certificate
