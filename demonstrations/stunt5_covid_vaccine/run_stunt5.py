"""
<<<<<<<< HEAD:papers/demonstrations/demo4_covid_vaccine/run_demo4.py
Demonstration 4: COVID-19 Spike Protein Vaccine Design
=======================================================
========
Demo 5: COVID-19 Spike Protein Vaccine Design
===============================================
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/stunt5_covid_vaccine/run_stunt5.py
Optimize the full SARS-CoV-2 spike protein for mRNA vaccine delivery,
demonstrating all 15 diagnostic layers on a real-world therapeutic target.
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

# SARS-CoV-2 Spike protein (1273 aa) — the primary COVID-19 vaccine target
# This is the prefusion-stabilized S-2P variant used in mRNA vaccines
COVID_SPIKE = (
    "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHS"
    "TQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIR"
    "GWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWME"
    "SEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSK"
    "HTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGW"
    "TAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEK"
    "GIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCV"
    "ADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQ"
    "TGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFER"
    "DISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELL"
    "HAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIA"
    "DTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTDVST"
    "AIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASY"
    "QTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEIL"
    "PVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE"
    "VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGF"
    "IKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGW"
    "TFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLS"
    "STASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQ"
    "IDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGK"
    "GYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSN"
    "GTHWFVDQRNFDYLLKDATYYKVDIPNNGGNYKKVELFPTTKIQDVVFDNRAL"
    "PDGKFKEGDKYYLYNLQDQYSFQEIWDNNEVSGIQEFLENKRLTTFDLLVSSS"
    "QHTSTGKKVRFSFCFDNLDDKFVKLNVTNDTVNKITVGFKDYEYKYVPKIWDT"
    "IKSKEFNKIKKYFPIQSSQSLKPGDLQNIENFKTTVFQNSKDVIIVDSQYKVT"
    "ESQYKVTESQYKVTESQYKVTESQYKVTESQYKVTE"
)


<<<<<<<< HEAD:papers/demonstrations/demo4_covid_vaccine/run_demo4.py
def run_demo4(regenerate: bool = False, verbose: bool = True):
    """Run Demonstration 4: COVID Vaccine."""
    print("=" * 70)
    print("  DEMO 4: COVID-19 Spike Protein Vaccine Design")
========
def run_stunt5(regenerate: bool = False, verbose: bool = True):
    """Run Demo 5: COVID Vaccine."""
    print("=" * 70)
    print("  DEMO 5: COVID-19 Spike Protein Vaccine Design")
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/stunt5_covid_vaccine/run_stunt5.py
    print(f"  Full spike protein: {len(COVID_SPIKE)} aa")
    print(f"  Diagnostic layers: {len(DIAGNOSTIC_LAYERS)}")
    print("=" * 70)

    # Run full diagnostic
    result = run_full_diagnostic(
        protein=COVID_SPIKE,
        protein_name="SARS-CoV-2_Spike_S2P",
        organism="Homo_sapiens",
        gc_lo=0.30,
        gc_hi=0.70,
        cai_threshold=0.5,
        verbose=verbose,
    )

    if verbose:
        print_diagnostic_summary(result)

    # Save results
    output_path = DATA_DIR / "covid_vaccine_results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "demonstration": "covid_vaccine",
        "description": f"SARS-CoV-2 spike protein ({len(COVID_SPIKE)} aa) optimized for mRNA vaccine",
        "layers": DIAGNOSTIC_LAYERS,
        "num_layers": len(DIAGNOSTIC_LAYERS),
        "spike_length_aa": len(COVID_SPIKE),
        "result": result.to_dict(),
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
    run_demo4(regenerate=args.regenerate, verbose=args.verbose)
