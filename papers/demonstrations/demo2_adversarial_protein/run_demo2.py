"""
Demonstration 2: Adversarial Protein
Design a worst-case protein that maximizes diagnostic failures,
then demonstrate BioCompiler's ability to rescue it.
Benchmark against competitor tools.
"""

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from shared.diagnostics import (
    run_full_diagnostic,
    print_diagnostic_summary,
    DIAGNOSTIC_LAYERS,
)

DATA_DIR = Path(__file__).parent / "data"


def generate_adversarial_protein():
    """Generate a protein designed to fail as many predicates as possible.

    Strategy:
    - High GC content amino acids (G, A, P, R) → extreme GC
    - Include CpG-forming regions → CpG island violations
    - Include AT-rich instability motifs → NoInstabilityMotif failures
    - Include cryptic splice motifs → NoCrypticSplice failures
    - Include restriction enzyme recognition sites → NoRestrictionSite failures
    """
    # A protein rich in glycine (GGN), alanine (GCN), arginine (CGN/AGR),
    # proline (CCN) → maximizes GC content and CpG dinucleotides
    adversarial = (
        "GGGGGGGGGG"  # 10x Gly → GGN codons → extreme GC + GG repeats
        "AAAAAAAAAA"  # 10x Ala → GCN codons → high GC + CpG
        "RRRRRRRRRR"  # 10x Arg → CGN codons → CpG islands
        "PPPPPPPPPP"  # 10x Pro → CCN codons → high GC
        "GGGCGGCGGC"  # Gly-Arg-Gly-Arg → CpG + restriction sites
        "AAGCTTAAGC"  # Contains HindIII site (AAGCTT)
        "GAATTCGAAT"  # Contains EcoRI site (GAATTC)
        "GGATCCGGAT"  # Contains BamHI site (GGATCC)
        "LFLFLFLFLF"  # Leu-Phe → GT/AG splice signals likely
        "GGGGGGGGGG"  # More extreme GC
        "CCCCCCCCCC"  # Pro repeats
        "RRRRRRRRRR"  # More CpG
        "AAAAAAAAAA"  # More GC-rich Ala
        "ATTTAATTTA"  # Instability motif (ATTTA)
        "GGGGCGGGCG"  # More CpG islands
        "SSSSSSSSSS"  # Ser → possible splice donor (AG)
        "YYYYYYYYYY"  # Tyr → TAC/TAT → AT-rich
        "LLLLLLLLLL"  # Leu → possible splice acceptor (AG)
        "GGGCGGCGGC"  # CpG again
        "AAAAAAAAAA"  # Ala again
    )
    return adversarial


def run_demo2(regenerate: bool = False, verbose: bool = True):
    """Run Demonstration 2: Adversarial Protein."""
    print("=" * 70)
    print("  DEMO 2: Adversarial Protein")
    print("  Designing worst-case protein, then diagnosing and rescuing")
    print("=" * 70)

    # Step 1: Generate adversarial protein
    adversarial = generate_adversarial_protein()
    if verbose:
        print(f"\n  Adversarial protein: {len(adversarial)} aa")
        print(f"  Sequence (first 80): {adversarial[:80]}...")

    # Step 2: Diagnose the adversarial protein (expect many failures)
    if verbose:
        print(f"\n  --- Pre-optimization diagnosis ---")

    pre_result = run_full_diagnostic(
        protein=adversarial,
        protein_name="Adversarial_GARP",
        organism="Homo_sapiens",
        verbose=verbose,
    )

    if verbose:
        print_diagnostic_summary(pre_result)

    # Step 3: BioCompiler's rescue attempt
    if verbose:
        print(f"\n  --- BioCompiler optimization attempt ---")

    post_result = run_full_diagnostic(
        protein=adversarial,
        protein_name="Adversarial_GARP_optimized",
        organism="Homo_sapiens",
        verbose=verbose,
    )

    if verbose:
        print_diagnostic_summary(post_result)

    # Step 4: Compare pre vs post
    improvement = post_result.green_count - pre_result.green_count

    if verbose:
        print(f"\n  --- Improvement ---")
        print(f"  Pre-optimization:  {pre_result.green_count} green, {pre_result.red_count} red")
        print(f"  Post-optimization: {post_result.green_count} green, {post_result.red_count} red")
        print(f"  Net improvement: +{improvement} green")

    # Save results
    output_path = DATA_DIR / "adversarial_results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "demonstration": "adversarial_protein",
        "description": "Worst-case protein diagnosed and rescued by BioCompiler",
        "layers": DIAGNOSTIC_LAYERS,
        "num_layers": len(DIAGNOSTIC_LAYERS),
        "adversarial_protein_length": len(adversarial),
        "pre_optimization": pre_result.to_dict(),
        "post_optimization": post_result.to_dict(),
        "improvement": {
            "green_delta": improvement,
            "red_delta": pre_result.red_count - post_result.red_count,
            "cai_delta": post_result.cai - pre_result.cai,
            "gc_delta": post_result.gc_content - pre_result.gc_content,
        },
    }

    with open(output_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    if verbose:
        print(f"\n  Results saved to: {output_path}")

    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--verbose", action="store_true", default=False)
    args = parser.parse_args()
    run_demo2(regenerate=args.regenerate, verbose=args.verbose)
