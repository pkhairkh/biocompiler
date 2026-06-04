"""
BioCompiler CSP Optimizer Demo
================================
Compares the CSP (Constraint Satisfaction Problem) solver against the
greedy optimizer on eGFP and HBB protein sequences.

The CSP solver uses Z3 or OR-Tools (via biocompiler.solver.dispatch) to
find globally optimal codon assignments that satisfy all constraints
simultaneously.  If no CSP backend is installed, the demo uses the
greedy optimizer with the "cai_first" strategy as a simulated CSP
baseline and compares it against the default greedy approach.

Usage:
    python -m biocompiler.examples.csp_optimize_demo
    python src/biocompiler/examples/csp_optimize_demo.py
"""

from __future__ import annotations

import importlib.util
import logging
import sys
import time
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# Allow running as ``python src/biocompiler/examples/csp_optimize_demo.py``
# from the repo root by adding src/ to sys.path when the package is not
# installed in development mode.
# ---------------------------------------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import os
    _repo_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    _src_dir = os.path.join(_repo_root, "src")
    if _src_dir not in sys.path:
        sys.path.insert(0, _src_dir)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Hardcoded protein sequences
# ---------------------------------------------------------------------------

# Enhanced Green Fluorescent Protein (eGFP) -- 239 aa
EGFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTT"
    "LTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKG"
    "IDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDG"
    "PVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

# Human Hemoglobin Beta (HBB) -- 147 aa
HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAH"
    "GKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQ"
    "AAYQKVVAGVANALAHKYH"
)

# ---------------------------------------------------------------------------
# CSP solver availability check
# ---------------------------------------------------------------------------

def _check_csp_available() -> dict[str, bool]:
    """Probe which CSP backends are available without importing them.

    Uses ``importlib.util.find_spec`` so heavy packages are never loaded
    just to test availability.
    """
    ortools_ok = importlib.util.find_spec("ortools") is not None
    z3_ok = importlib.util.find_spec("z3") is not None
    return {
        "ortools": ortools_ok,
        "z3": z3_ok,
        "any": ortools_ok or z3_ok,
    }


CSP_BACKENDS = _check_csp_available()

# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class OptimizationRun:
    """Captures the outcome of one optimizer run for comparison."""
    name: str
    protein_name: str
    protein: str
    sequence: str
    gc_content: float
    cai: float
    elapsed_s: float
    passed_predicates: list[str] = field(default_factory=list)
    failed_predicates: list[str] = field(default_factory=list)
    solver_used: str = "greedy"

    @property
    def total_predicates(self) -> int:
        return len(self.passed_predicates) + len(self.failed_predicates)

    @property
    def pass_rate(self) -> float:
        total = self.total_predicates
        return len(self.passed_predicates) / total if total else 0.0


# ---------------------------------------------------------------------------
# Greedy optimizer wrapper
# ---------------------------------------------------------------------------

def run_greedy(protein: str, protein_name: str) -> OptimizationRun:
    """Run the built-in greedy optimizer and collect results."""
    from biocompiler.optimization import optimize_sequence

    t0 = time.perf_counter()
    result = optimize_sequence(
        target_protein=protein,
        organism="Homo_sapiens",
    )
    elapsed = time.perf_counter() - t0

    return OptimizationRun(
        name="Greedy",
        protein_name=protein_name,
        protein=protein,
        sequence=result.sequence,
        gc_content=result.gc_content,
        cai=result.cai,
        elapsed_s=elapsed,
        passed_predicates=list(result.satisfied_predicates),
        failed_predicates=list(result.failed_predicates),
        solver_used="greedy",
    )


# ---------------------------------------------------------------------------
# CSP optimizer -- uses biocompiler.solver.dispatch when available
# ---------------------------------------------------------------------------

def _try_solver_dispatch(protein: str) -> Optional[str]:
    """Attempt to solve via biocompiler.solver.dispatch.csp_optimize.

    Returns the optimized DNA sequence, or None if the solver package
    is unavailable or fails.
    """
    try:
        from biocompiler.solver.dispatch import csp_optimize, is_csp_available
    except (ImportError, Exception):
        return None

    avail = is_csp_available()
    if not avail.get("any", False):
        return None

    try:
        result = csp_optimize(
            protein=protein,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
        )
        if result and hasattr(result, "sequence") and result.sequence:
            return result.sequence
    except Exception:
        logger.warning("CSP solver dispatch failed, skipping", exc_info=True)

    return None


def _fallback_csp(protein: str) -> str:
    """Greedy with cai_first strategy as a simulated CSP baseline.

    Used when no CSP backend is installed.  This uses a different
    optimisation strategy (maximise CAI first, then fix constraints)
    which approximates the global-optimisation behaviour of a true
    CSP solver.
    """
    from biocompiler.optimization import optimize_sequence

    result = optimize_sequence(
        target_protein=protein,
        organism="Homo_sapiens",
        strategy="cai_first",
    )
    return result.sequence


def run_csp(protein: str, protein_name: str) -> OptimizationRun:
    """Run the CSP optimizer and collect results.

    Tries the real CSP solver first, then falls back to the greedy
    cai_first strategy as a simulated baseline.
    """
    from biocompiler.translation import compute_cai
    from biocompiler.type_system import evaluate_all_predicates

    # Determine which approach to use
    solver_label = "none"
    t0 = time.perf_counter()

    sequence = _try_solver_dispatch(protein)

    if sequence is not None:
        solver_label = "csp-dispatch"
    else:
        # No real CSP backend -- use greedy cai_first as baseline
        sequence = _fallback_csp(protein)
        solver_label = "simulated (cai_first)"

    elapsed = time.perf_counter() - t0

    # Compute metrics
    gc = (sequence.count("G") + sequence.count("C")) / max(len(sequence), 1)
    gc = round(gc, 4)

    try:
        cai_val = compute_cai(sequence, "Homo_sapiens")
    except Exception:
        cai_val = 0.0

    # Evaluate predicates
    boundaries = [(0, len(sequence))]
    pred_results = evaluate_all_predicates(
        seq=sequence,
        boundaries=boundaries,
        organism="Homo_sapiens",
        gc_lo=0.30,
        gc_hi=0.70,
    )

    passed = [r.predicate for r in pred_results if getattr(r, "passed", True)]
    failed = [r.predicate for r in pred_results if not getattr(r, "passed", True)]

    return OptimizationRun(
        name=f"CSP ({solver_label})",
        protein_name=protein_name,
        protein=protein,
        sequence=sequence,
        gc_content=gc,
        cai=cai_val,
        elapsed_s=elapsed,
        passed_predicates=passed,
        failed_predicates=failed,
        solver_used=solver_label,
    )


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _divider(char: str = "-", width: int = 72) -> str:
    return char * width


def _format_predicate_table(run: OptimizationRun) -> str:
    """Format pass/fail summary for a single run."""
    lines: list[str] = []
    all_preds = sorted(
        [(p, "PASS") for p in run.passed_predicates]
        + [(p, "FAIL") for p in run.failed_predicates]
    )
    for name, status in all_preds:
        tag = "[PASS]" if status == "PASS" else "[FAIL]"
        lines.append(f"  {tag:8s} {name}")
    if not all_preds:
        lines.append("  (no predicate results)")
    return "\n".join(lines)


def _format_comparison(greedy: OptimizationRun, csp: OptimizationRun) -> str:
    """Format a side-by-side comparison of two runs."""
    lines: list[str] = []
    lines.append(_divider("=", 72))
    lines.append(f"  Comparison: {greedy.name} vs {csp.name}")
    lines.append(f"  Protein: {greedy.protein_name} ({len(greedy.protein)} aa)")
    lines.append(_divider("-", 72))

    # Metrics comparison
    lines.append("")
    lines.append("  Metric           | Greedy          | CSP")
    lines.append("  -----------------+-----------------+-----------------")

    def _row(label: str, g_val: str, c_val: str) -> str:
        return f"  {label:17s}| {g_val:15s} | {c_val:15s}"

    lines.append(_row("GC content", f"{greedy.gc_content:.4f}",
                       f"{csp.gc_content:.4f}"))
    lines.append(_row("CAI", f"{greedy.cai:.4f}", f"{csp.cai:.4f}"))
    lines.append(_row("Time (s)", f"{greedy.elapsed_s:.3f}",
                       f"{csp.elapsed_s:.3f}"))
    lines.append(_row("Passed preds", str(len(greedy.passed_predicates)),
                       str(len(csp.passed_predicates))))
    lines.append(_row("Failed preds", str(len(greedy.failed_predicates)),
                       str(len(csp.failed_predicates))))
    lines.append(_row("Pass rate", f"{greedy.pass_rate:.1%}",
                       f"{csp.pass_rate:.1%}"))
    lines.append(_row("Solver", greedy.solver_used, csp.solver_used))
    lines.append(_row("Seq length (bp)", str(len(greedy.sequence)),
                       str(len(csp.sequence))))

    # Delta calculations
    lines.append("")
    lines.append("  Deltas (CSP - Greedy):")

    gc_delta = csp.gc_content - greedy.gc_content
    cai_delta = csp.cai - greedy.cai
    time_delta = csp.elapsed_s - greedy.elapsed_s
    pass_delta = len(csp.passed_predicates) - len(greedy.passed_predicates)

    lines.append(f"    GC content:  {gc_delta:+.4f}")
    lines.append(f"    CAI:         {cai_delta:+.4f}")
    lines.append(f"    Time:        {time_delta:+.3f} s")
    lines.append(f"    Passed pred: {pass_delta:+d}")

    return "\n".join(lines)


def _format_run_detail(run: OptimizationRun) -> str:
    """Format detailed output for a single optimization run."""
    lines: list[str] = []
    lines.append(_divider("-", 72))
    lines.append(f"  {run.name} Optimization - {run.protein_name}")
    lines.append(_divider("-", 72))
    lines.append(f"  Sequence length: {len(run.sequence)} bp")
    lines.append(f"  GC content:      {run.gc_content:.4f}")
    lines.append(f"  CAI:             {run.cai:.4f}")
    lines.append(f"  Elapsed time:    {run.elapsed_s:.3f} s")
    lines.append(f"  Solver:          {run.solver_used}")
    lines.append(f"  Passed:          {len(run.passed_predicates)}"
                 f"/{run.total_predicates}")
    lines.append(f"  Failed:          {len(run.failed_predicates)}"
                 f"/{run.total_predicates}")
    lines.append("")
    lines.append("  Predicate Results:")
    lines.append(_format_predicate_table(run))
    lines.append("")

    # Show first 80 bp of sequence
    preview = run.sequence[:80]
    suffix = "..." if len(run.sequence) > 80 else ""
    lines.append(f"  Sequence preview: {preview}{suffix}")
    lines.append(_divider("-", 72))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main demo
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the CSP vs Greedy optimization comparison demo."""
    print()
    print("=" * 72)
    print("  BioCompiler CSP Optimizer Demo")
    print("  Comparing CSP solver vs Greedy optimizer on eGFP and HBB")
    print("=" * 72)
    print()

    # Report backend availability
    if CSP_BACKENDS["ortools"]:
        print("  [OK] OR-Tools CP-SAT solver is available")
    else:
        print("  [--] OR-Tools CP-SAT not installed")

    if CSP_BACKENDS["z3"]:
        print("  [OK] Z3 SMT solver is available")
    else:
        print("  [--] Z3 SMT solver not installed")

    if CSP_BACKENDS["any"]:
        print("  CSP backend available -- will use real CSP solver")
    else:
        print("  No CSP backend installed -- using simulated CSP fallback")
        print("  Install with: pip install z3-solver  OR  pip install ortools")
    print()

    proteins = [
        ("eGFP", EGFP_PROTEIN),
        ("HBB", HBB_PROTEIN),
    ]

    all_comparisons: list[tuple[OptimizationRun, OptimizationRun]] = []

    for protein_name, protein in proteins:
        print(_divider("=", 72))
        print(f"  Optimizing: {protein_name} ({len(protein)} aa)")
        print(_divider("=", 72))
        print()

        # --- Greedy run ---
        print(f"  Running greedy optimizer on {protein_name}...")
        try:
            greedy_run = run_greedy(protein, protein_name)
            print(f"  Greedy done in {greedy_run.elapsed_s:.3f} s")
        except Exception as exc:
            print(f"  Greedy optimizer FAILED: {exc}")
            continue

        # --- CSP run ---
        print(f"  Running CSP optimizer on {protein_name}...")
        try:
            csp_run = run_csp(protein, protein_name)
            print(f"  CSP done in {csp_run.elapsed_s:.3f} s")
        except Exception as exc:
            print(f"  CSP optimizer FAILED: {exc}")
            print(f"  Showing greedy results only.")
            print()
            print(_format_run_detail(greedy_run))
            continue

        all_comparisons.append((greedy_run, csp_run))

        # --- Detailed results ---
        print()
        print(_format_run_detail(greedy_run))
        print()
        print(_format_run_detail(csp_run))
        print()
        print(_format_comparison(greedy_run, csp_run))
        print()

    # --- Grand summary ---
    if all_comparisons:
        print()
        print("=" * 72)
        print("  Grand Summary")
        print("=" * 72)
        print()

        header = (
            f"  {'Protein':<10s} {'Solver':<20s} "
            f"{'CAI':>7s} {'GC':>7s} {'Pass':>5s} "
            f"{'Fail':>5s} {'Time':>8s}"
        )
        print(header)
        print("  " + "-" * (len(header) - 2))

        for greedy, csp in all_comparisons:
            for run in (greedy, csp):
                print(
                    f"  {run.protein_name:<10s} "
                    f"{run.solver_used:<20s} "
                    f"{run.cai:>7.4f} "
                    f"{run.gc_content:>7.4f} "
                    f"{len(run.passed_predicates):>5d} "
                    f"{len(run.failed_predicates):>5d} "
                    f"{run.elapsed_s:>7.3f}s"
                )
            print()

        # Verdict
        print(_divider("-", 72))
        csp_wins = 0
        greedy_wins = 0
        for greedy, csp in all_comparisons:
            # Score: weighted sum of predicate pass count + CAI bonus
            g_score = len(greedy.passed_predicates) + greedy.cai
            c_score = len(csp.passed_predicates) + csp.cai
            if c_score > g_score:
                csp_wins += 1
            elif g_score > c_score:
                greedy_wins += 1

        if csp_wins > greedy_wins:
            print("  CSP optimizer produces better results overall.")
        elif greedy_wins > csp_wins:
            print("  Greedy optimizer produces better results overall.")
        else:
            print("  Both optimizers produce comparable results.")

        if not CSP_BACKENDS["any"]:
            print()
            print("  NOTE: Install a CSP backend for true constraint-based")
            print("  optimisation (globally optimal codon assignment):")
            print("    pip install z3-solver")
            print("    pip install ortools")
            print()
            print("  The simulated CSP fallback uses the greedy optimizer")
            print("  with a different strategy (cai_first), so differences")
            print("  may be minimal without a real CSP backend.")

        print()
        print("=" * 72)
        print("  Demo complete.")
        print("=" * 72)


if __name__ == "__main__":
    main()
