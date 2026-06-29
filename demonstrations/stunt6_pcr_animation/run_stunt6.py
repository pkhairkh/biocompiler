"""
<<<<<<<< HEAD:papers/demonstrations/demo5_pcr_animation/run_demo5.py
Demonstration 5: Primer Compatibility Demo
===========================================
========
Demo 6: Primer Compatibility Demo
==================================
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/stunt6_pcr_animation/run_stunt6.py
Demonstrate the primer compatibility diagnostic layer —
a critical but often overlooked aspect of mRNA vaccine design.
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

# Test sequences with different GC profiles
TEST_SEQUENCES = {
    "Normal_GC_Insulin": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    "High_GC_Custom": "GGGGGGGGGGAAAAAAAAAARRRRRRRRRRPPPPPPPPPPGGGGCGGGCGAAGCTTAAGCGAATTCGAATGGATCCGGATCCCCCCCCCCRRRRRRRRRRAAAAAAAAAAGGGGGGGGGG",
    "Low_GC_Collagen": "GPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPP",
}


<<<<<<<< HEAD:papers/demonstrations/demo5_pcr_animation/run_demo5.py
def run_demo5(regenerate: bool = False, verbose: bool = True):
    """Run Demonstration 5: Primer Compatibility Demo."""
    print("=" * 70)
    print("  DEMO 5: Primer Compatibility Demo")
========
def run_stunt6(regenerate: bool = False, verbose: bool = True):
    """Run Demo 6: Primer Compatibility Demo."""
    print("=" * 70)
    print("  DEMO 6: Primer Compatibility Demo")
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/stunt6_pcr_animation/run_stunt6.py
    print("  Testing PCR compatibility across different GC profiles")
    print("=" * 70)

    results = []

    for name, protein in TEST_SEQUENCES.items():
        if verbose:
            print(f"\n>>> {name} ({len(protein)} aa)")

        result = run_full_diagnostic(
            protein=protein,
            protein_name=name,
            organism="Homo_sapiens",
            verbose=verbose,
        )
        results.append(result)

        if verbose:
            print(f"  GC: {result.gc_content:.1%}, Primer compat: {'YES' if result.primer_compatibility else 'NO'}")

    # Save results
    output_path = DATA_DIR / "primer_compat_results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "demonstration": "primer_compatibility",
        "description": "PCR/primer compatibility across different GC profiles",
        "results": [r.to_dict() for r in results],
        "primer_compat_rules": {
            "gc_range": "35-65% for standard PCR",
            "amplicon_length": "<5 kb for standard Taq polymerase",
            "primer_design": "Tm 55-65°C, GC clamp, no self-complementarity",
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
    run_demo5(regenerate=args.regenerate, verbose=args.verbose)
