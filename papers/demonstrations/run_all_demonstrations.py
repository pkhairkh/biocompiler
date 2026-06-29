"""
BioCompiler Paper Demonstrations — Main Runner
Run all paper demonstrations for the BioCompiler POPL 2027 paper.

Demonstrations:
  1: Therapeutic Autopsy — diagnose 5 FDA drugs
  2: Adversarial Protein — worst-case protein rescue
  3: Atlas — large-scale compilation survey
  4: COVID Vaccine — SARS-CoV-2 spike optimization
  5: Primer Compatibility — PCR demo

Usage:
  python run_all_demonstrations.py              # Run all
  python run_all_demonstrations.py --demo 1     # Run specific demonstration
  python run_all_demonstrations.py --regenerate # Force regeneration of data
  python run_all_demonstrations.py --verbose    # Verbose output
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def run_demo(demo_num: int, regenerate: bool, verbose: bool):
    """Run a specific demonstration by number."""
    if demo_num == 1:
        from demo1_therapeutic_autopsy.run_demo1 import run_demo1
        return run_demo1(regenerate=regenerate, verbose=verbose)
    elif demo_num == 2:
        from demo2_adversarial_protein.run_demo2 import run_demo2
        return run_demo2(regenerate=regenerate, verbose=verbose)
    elif demo_num == 3:
        from demo3_atlas.run_demo3 import run_demo3
        return run_demo3(regenerate=regenerate, verbose=verbose)
    elif demo_num == 4:
        from demo4_covid_vaccine.run_demo4 import run_demo4
        return run_demo4(regenerate=regenerate, verbose=verbose)
    elif demo_num == 5:
        from demo5_pcr_animation.run_demo5 import run_demo5
        return run_demo5(regenerate=regenerate, verbose=verbose)
    else:
        print(f"Unknown demonstration: {demo_num}. Available: 1, 2, 3, 4, 5")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="BioCompiler Paper Demonstrations Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--demo", type=int, nargs="+",
        help="Demonstration number(s) to run (1, 2, 3, 4, 5). Default: all"
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

    demos = args.demo or [1, 2, 3, 4, 5]

    print("=" * 70)
    print("  BioCompiler Paper Demonstrations")
    print(f"  Running: {', '.join(f'Demo {s}' for s in demos)}")
    print(f"  Regenerate: {args.regenerate}")
    print("=" * 70)

    t_start = time.time()
    results = {}

    for s in demos:
        print(f"\n{'#'*70}")
        print(f"  Running Demo {s}...")
        print(f"{'#'*70}")

        try:
            t0 = time.time()
            result = run_demo(s, args.regenerate, args.verbose)
            elapsed = time.time() - t0
            results[s] = {"status": "PASS", "time_s": elapsed, "result": result}
            print(f"\n  Demo {s}: PASS ({elapsed:.1f}s)")
        except Exception as e:
            results[s] = {"status": "FAIL", "error": str(e)}
            print(f"\n  Demo {s}: FAIL ({e})")

    total_time = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"  ALL DEMOS COMPLETE")
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
