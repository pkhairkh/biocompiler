"""
<<<<<<<< HEAD:papers/demonstrations/demo3_atlas/run_demo3.py
Demonstration 3: Atlas — Large-Scale Compilation Survey
========================================================
========
Demo 3: Atlas — Large-Scale Compilation Survey
================================================
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/stunt3_atlas/run_stunt3.py
Optimize a representative set of proteins across multiple organisms
to measure pass rates per predicate and identify hardest constraints.
"""

import json
import sys
import time
import random
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent))

from shared.diagnostics import (
    run_full_diagnostic,
    DiagnosticResult,
    CANONICAL_PREDICATES,
    DIAGNOSTIC_LAYERS,
)

DATA_DIR = Path(__file__).parent / "data"

# Representative proteins (shorter for speed)
PROTEIN_PANEL = {
    "Insulin": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    "GH": "FPTIPLSRLFDNAMLRAGIVHFCIDRAGRKFDLQKTEGHVVVLAEAMKMDGVVVRSLRLNEHKVDAYLRNLKEFHEKFQKLPSPEYKEFTKRINLLQDSMLSVSRLGHAIEKLTSEAPKLLDENQK",
    "HBB": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
    "EGFP": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "BSA": "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
}

# Organisms to test
ORGANISMS = {
    "mammalian": ["Homo_sapiens", "Mus_musculus", "CHO"],
    "bacterial": ["E_coli", "Bacillus_subtilis"],
}


def run_atlas_sample(
    num_proteins_per_organism: int = 3,
    verbose: bool = True,
) -> dict:
    """Run a scaled-down Atlas compilation survey."""
    random.seed(42)

    protein_names = list(PROTEIN_PANEL.keys())
    all_organisms = []
    for group, orgs in ORGANISMS.items():
        all_organisms.extend(orgs)

    results = []
    pass_rates = defaultdict(lambda: {"pass": 0, "fail": 0, "uncertain": 0, "total": 0})

    total_runs = len(all_organisms) * min(num_proteins_per_organism, len(protein_names))
    done = 0

    for organism in all_organisms:
        # Select proteins for this organism
        selected = protein_names[:num_proteins_per_organism]

        for prot_name in selected:
            protein = PROTEIN_PANEL[prot_name]
            done += 1

            if verbose:
                print(f"  [{done}/{total_runs}] {prot_name} × {organism}...", end="", flush=True)

            try:
                result = run_full_diagnostic(
                    protein=protein,
                    protein_name=f"{prot_name}_{organism}",
                    organism=organism,
                    verbose=False,
                )

                # Record per-predicate pass rates
                for p in result.canonical_predicates:
                    pred_name = p["predicate"]
                    pass_rates[pred_name]["total"] += 1
                    color = {"green": "pass", "red": "fail", "amber": "uncertain"}[
                        {"PASS": "green", "LIKELY_PASS": "green",
                         "FAIL": "red", "LIKELY_FAIL": "red"}.get(
                            p["verdict"].upper(), "amber"
                        )
                    ]
                    pass_rates[pred_name][color] += 1

                # Extended predicates
                pass_rates["BiosecurityScreening"]["total"] += 1
                pass_rates["BiosecurityScreening"]["pass" if result.biosecurity_passed else "fail"] += 1

                pass_rates["Immunogenicity"]["total"] += 1
                imm_color = "pass" if result.immunogenicity_risk in ("low", "moderate") else "fail"
                pass_rates["Immunogenicity"][imm_color] += 1

                pass_rates["PrimerCompatibility"]["total"] += 1
                pass_rates["PrimerCompatibility"]["pass" if result.primer_compatibility else "fail"] += 1

                results.append(result.to_dict())

                if verbose:
                    print(f" {result.green_count}G/{result.red_count}R")

            except Exception as e:
                if verbose:
                    print(f" ERROR: {e}")

    # Compute summary
    layer_pass_rates = {}
    for layer_name, counts in pass_rates.items():
        if counts["total"] > 0:
            rate = counts["pass"] / counts["total"]
            layer_pass_rates[layer_name] = round(rate, 4)

    # Sort by pass rate (hardest first)
    sorted_layers = sorted(layer_pass_rates.items(), key=lambda x: x[1])

    # All-pass rate
    all_pass_count = sum(1 for r in results if r.get("all_predicates_pass", False))

    summary = {
        "demonstration": "atlas",
        "description": f"Large-scale compilation survey: {total_runs} runs",
        "layers": DIAGNOSTIC_LAYERS,
        "num_layers": len(DIAGNOSTIC_LAYERS),
        "total_compilations": total_runs,
        "all_pass_rate": round(all_pass_count / max(total_runs, 1), 4),
        "all_pass_count": all_pass_count,
        "layer_pass_rates": layer_pass_rates,
        "hardest_layers": [name for name, rate in sorted_layers[:5]],
        "easiest_layers": [name for name, rate in sorted_layers[-5:]],
        "results": results,
    }

    return summary


<<<<<<<< HEAD:papers/demonstrations/demo3_atlas/run_demo3.py
def run_demo3(regenerate: bool = False, verbose: bool = True):
    """Run Demonstration 3: Atlas."""
========
def run_stunt3(regenerate: bool = False, verbose: bool = True):
    """Run Demo 3: Atlas."""
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/stunt3_atlas/run_stunt3.py
    print("=" * 70)
    print("  DEMO 3: Atlas — Large-Scale Compilation Survey")
    print("  Measuring predicate pass rates across organisms and proteins")
    print("=" * 70)

    summary = run_atlas_sample(
        num_proteins_per_organism=3,
        verbose=verbose,
    )

    # Save results
    output_path = DATA_DIR / "atlas_results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    if verbose:
        print(f"\n{'='*70}")
        print(f"  DEMO 3 SUMMARY")
        print(f"{'='*70}")
        print(f"  Total compilations: {summary['total_compilations']}")
        print(f"  All-pass rate: {summary['all_pass_rate']:.1%}")
        print(f"  Hardest layers: {', '.join(summary['hardest_layers'])}")
        print(f"  Easiest layers: {', '.join(summary['easiest_layers'])}")
        print(f"\n  Per-layer pass rates:")
        for layer, rate in sorted(summary["layer_pass_rates"].items(), key=lambda x: x[1]):
            bar = "█" * int(rate * 30)
            print(f"    {layer:30s} {rate:6.1%} {bar}")
        print(f"\n  Results saved to: {output_path}")

    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--verbose", action="store_true", default=False)
    args = parser.parse_args()
    run_demo3(regenerate=args.regenerate, verbose=args.verbose)
