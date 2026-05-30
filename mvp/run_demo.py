"""
BioCompiler MVP — Entry Point

Runs the complete pipeline demo for both target proteins:
  1. Human Insulin (51 aa) — small, fast to optimize
  2. eGFP (239 aa) — larger, more challenging

All output is deterministic and reproducible.
"""

from biocompiler_mvp import run_pipeline, HUMAN_INSULIN, EGFP


def main() -> None:
    """Run the BioCompiler MVP demo."""
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║                    BioCompiler MVP Demo                             ║")
    print("║  MaxEntScan + CSP Optimization with Three-Valued Type System       ║")
    print("║                                                                     ║")
    print("║  Pipeline:                                                          ║")
    print("║    1. CSP optimization (z3) → DNA sequence                          ║")
    print("║    2. MaxEntScan → splice site scoring                              ║")
    print("║    3. Type checking → 7 predicates (PASS/FAIL/UNCERTAIN)            ║")
    print("║    4. Certificate generation → JSON                                 ║")
    print("║    5. Certificate verification → independent re-check               ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    # Demo 1: Human Insulin (51 aa)
    print("╭──────────────────────────────────────────────────────────────────────╮")
    print("│ DEMO 1: Human Insulin (51 amino acids)                              │")
    print("╰──────────────────────────────────────────────────────────────────────╯")
    run_pipeline(
        target_protein=HUMAN_INSULIN,
        organism="human",
        protein_name="Human_Insulin",
    )

    print("\n\n")

    # Demo 2: eGFP (239 aa)
    print("╭──────────────────────────────────────────────────────────────────────╮")
    print("│ DEMO 2: eGFP (239 amino acids)                                      │")
    print("╰──────────────────────────────────────────────────────────────────────╯")
    run_pipeline(
        target_protein=EGFP,
        organism="human",
        protein_name="eGFP",
    )

    print("\n")
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║                      Demo Complete                                  ║")
    print("║  Certificates saved as JSON files in:                               ║")
    print("║    /home/z/my-project/download/biocompiler-mvp/                     ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")


if __name__ == "__main__":
    main()
