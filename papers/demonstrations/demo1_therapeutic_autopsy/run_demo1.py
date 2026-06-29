"""
Demonstration 1: Therapeutic Autopsy
Diagnose 5 FDA-approved therapeutic proteins through BioCompiler's
full 15-diagnostic-layer pipeline (12 canonical predicates +
biosecurity, immunogenicity, primer compatibility).

Paper narrative: "If these approved drugs were designed today with
BioCompiler, which hidden risks would we catch?"
"""

import json
import sys
import time
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from shared.diagnostics import (
    run_full_diagnostic,
    print_diagnostic_summary,
    DiagnosticResult,
    CANONICAL_PREDICATES,
    DIAGNOSTIC_LAYERS,
)

DATA_DIR = Path(__file__).parent / "data"

# FDA-approved therapeutic proteins
FDA_DRUGS = {
    "Insulin (Humalin)": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    "EPO (Epogen)": "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRTLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR",
    "Factor VIII (Advate)": "MQKVCLNQFSLYKSSFYLMQANLRLTDSWQMMKNDLSPRHCQSLESIRTFCSFWNSLFKEYKFYVDVSQDGKRQVFQNSNNQGFNKTIVNNLIQYMGYYITGFHLHELRQNGKLISFTNFVKLYHFTKQNIITHLIKSSYGDVAVGLGASHVSIFKAVNPKFVNFLAFLNNEPSLAKVIENLHQLPSSQFQKTSYDYTFNHNLTMMHSFAVHKRIYVNFIDNKTITYNENLNNLTFTFKHSYVLQRDLITYVDKDKKSDQIHEQYVSVNVNQSTVDLKLSFQSYINSCLVPVDNNTTTIWKENLSLPVNKETVQYLSFTNVISYSVSGVDSFQISSKILKSLENQISNQVKNSKILTFNVTSQNVDNVSISTLTSKLCIFKNSYSKNISYTSVLVNWVNKQVTITKLNKTIQGLNKTISNSYLFQNFNITTKDVKVNLTGLNVVTNNIFQNKKSISTKLNVTVNISYSKLFQNISSYITDSNLNFNITSIENNTLSKFKDLVTNLKIQSHPNLSAIFESLKQNISNLSISLFNISNVSKSQFNITISYITNDLTLKNKNLSFIAENNKFIKNKLNSHLINILNFSNKINSFSINLKISNLSKPKIDINNINISLISNKFSIENDIKKNINLSYINNKNIKQFSINCSLNIQINNSKINISLSYITNNKIPLNFNINLILNNSNKINLSKFKNDLFINKININSILNINKNISISTIKSLLISNISDLINELNLSNFLISNLSIPTIKQINHSISNITKISFSINKLFNLSNNLINIFNSLKISNKNINYFKINSFSISNLNYLSIDNINLTNIKNYLSNFNKSIYISNTIKNYLSIKNYLSLFNKKISSTINLFNITNNLINNFNFNINSINSKINKNLSIQNINSIN",
    "G-CSF (Neupogen)": "MTPLGPASSLPQSFLLKCLEQVRKIQGDGAALQEKLCATYKLCHPEELVLLGHSLGIPWAPLSSCPSQALQLAGCLSQLHSGLFLYQGLLQALEGISPELGPTLDTLQLDVADFATTIWQQMEELGMAPALQPTQGAMPAFASAFQRRAGGVLVASHLQSFLEVSYRVLRHLAQP",
    "Interferon-alpha (Intron A)": "MALTFGLLVASLVTLSSSVSGCDLPQTHSLGNRRAFSLLTQYRRLQLLMSNLKLDKCSESPSPEPIFSLGERLLQLSLQSKRRKDLCRAFHEFLQEVSSLAQRSLDFQNELTAVAEKMHQFSLVNKESDLKKIIQQTVDFLQENKDTGALSRSADVLFSFTQLLDSLSASGLSFDQVLQKADSLIPESLESIRSLTFQQQLLRKLPQNQQFFSTLNNNFKDKDKDKKDKKDKKDKDKK",
}


def run_demo1(regenerate: bool = False, verbose: bool = True):
    """Run Demonstration 1: Therapeutic Autopsy."""
    print("=" * 70)
    print("  DEMO 1: Therapeutic Autopsy")
    print("  Diagnosing 5 FDA-approved drugs through 15 diagnostic layers")
    print("=" * 70)

    results = []

    for name, protein in FDA_DRUGS.items():
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
            print_diagnostic_summary(result)

    # Save results
    output_path = DATA_DIR / "diagnostic_results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "demonstration": "therapeutic_autopsy",
        "description": "5 FDA-approved drugs diagnosed through 15 diagnostic layers",
        "layers": DIAGNOSTIC_LAYERS,
        "num_layers": len(DIAGNOSTIC_LAYERS),
        "results": [r.to_dict() for r in results],
        "summary": {
            "total_drugs": len(results),
            "all_pass": sum(1 for r in results if r.all_predicates_pass),
            "avg_green": sum(r.green_count for r in results) / len(results),
            "avg_red": sum(r.red_count for r in results) / len(results),
            "avg_time_s": sum(r.optimization_time_s for r in results) / len(results),
        },
    }

    with open(output_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    if verbose:
        print(f"\n{'='*70}")
        print(f"  DEMO 1 SUMMARY")
        print(f"{'='*70}")
        print(f"  Drugs diagnosed: {len(results)}")
        print(f"  Diagnostic layers: {len(DIAGNOSTIC_LAYERS)}")
        print(f"  All-pass: {summary['summary']['all_pass']}/{len(results)}")
        print(f"  Avg green/total: {summary['summary']['avg_green']:.1f}/{len(DIAGNOSTIC_LAYERS)}")
        print(f"  Avg time: {summary['summary']['avg_time_s']:.1f}s")
        print(f"  Results saved to: {output_path}")

    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--verbose", action="store_true", default=False)
    args = parser.parse_args()
    run_demo1(regenerate=args.regenerate, verbose=args.verbose)
