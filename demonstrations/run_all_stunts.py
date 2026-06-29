"""
BioCompiler Demonstrations — Main Runner
=========================================
Run all paper demonstrations for the BioCompiler POPL 2027 paper.

Demonstrations:
  1: Therapeutic Autopsy — diagnose 5 FDA drugs
  2: Adversarial Protein — worst-case protein rescue
  3: Atlas — large-scale compilation survey
  5: COVID Vaccine — SARS-CoV-2 spike optimization
  6: Primer Compatibility — PCR demo

Usage:
  python run_all_stunts.py              # Run all
  python run_all_stunts.py --stunt 1    # Run specific demonstration
  python run_all_stunts.py --regenerate # Force regeneration of data
  python run_all_stunts.py --verbose    # Verbose output
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def run_stunt(stunt_num: int, regenerate: bool, verbose: bool):
    """Run a specific demonstration by number."""
    if stunt_num == 1:
        from stunt1_therapeutic_autopsy.run_stunt1 import run_stunt1
        return run_stunt1(regenerate=regenerate, verbose=verbose)
    elif stunt_num == 2:
        from stunt2_adversarial_protein.run_stunt2 import run_stunt2
        return run_stunt2(regenerate=regenerate, verbose=verbose)
    elif stunt_num == 3:
        from stunt3_atlas.run_stunt3 import run_stunt3
        return run_stunt3(regenerate=regenerate, verbose=verbose)
    elif stunt_num == 5:
        from stunt5_covid_vaccine.run_stunt5 import run_stunt5
        return run_stunt5(regenerate=regenerate, verbose=verbose)
    elif stunt_num == 6:
        from stunt6_pcr_animation.run_stunt6 import run_stunt6
        return run_stunt6(regenerate=regenerate, verbose=verbose)
    else:
        print(f"Unknown demonstration: {stunt_num}. Available: 1, 2, 3, 5, 6")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="BioCompiler Demonstrations Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--stunt", type=int, nargs="+",
        help="Demonstration number(s) to run (1, 2, 3, 5, 6). Default: all"
    )
    parser.add_argument(
        "--regenerate", action="store_true",
        help="Force regeneration of output data"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=True,
        help="Verbose output (default: True)"
    )

    args = parser.parse_args()

    stunts = args.stunt or [1, 2, 3, 5, 6]

    print("=" * 70)
    print("  BioCompiler Demonstrations")
    print(f"  Running: {', '.join(f'Demo {s}' for s in stunts)}")
    print(f"  Regenerate: {args.regenerate}")
    print("=" * 70)

    t_start = time.time()
    results = {}

    for s in stunts:
        print(f"\n{'#'*70}")
        print(f"  Running Demo {s}...")
        print(f"{'#'*70}")

        try:
            t0 = time.time()
            result = run_stunt(s, args.regenerate, args.verbose)
            elapsed = time.time() - t0
            results[s] = {"status": "PASS", "time_s": elapsed, "result": result}
            print(f"\n  Demo {s}: PASS ({elapsed:.1f}s)")
        except Exception as e:
            results[s] = {"status": "FAIL", "error": str(e)}
            print(f"\n  Demo {s}: FAIL ({e})")

    total_time = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"  ALL DEMONSTRATIONS COMPLETE")
    print(f"{'='*70}")

    for s, r in results.items():
        status = r["status"]
        time_s = f"{r.get('time_s', 0):.1f}s" if status == "PASS" else "N/A"
        print(f"  Demo {s}: {status} ({time_s})")

    passed = sum(1 for r in results.values() if r["status"] == "PASS")
    print(f"\n  Total: {passed}/{len(results)} passed in {total_time:.1f}s")

    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
