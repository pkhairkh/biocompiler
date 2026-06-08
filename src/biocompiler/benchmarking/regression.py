"""
BioCompiler Benchmark Regression Detection
============================================
Detects performance and quality regressions by comparing current benchmark
results against a baseline from a previous release.

Regression thresholds:
  - Performance regression: > 10% slowdown in optimization time
  - Quality regression:     > 5% CAI drop
  - Constraint regression:  > 5% drop in constraint satisfaction rate
  - GC regression:          > 10% deviation from target GC

The module also generates human-readable regression reports suitable for
CI comments and release notes.

Usage::

    from biocompiler.benchmarking.regression import (
        compare_results,
        generate_regression_report,
        RegressionReport,
    )

    comparison = compare_results(current_results, baseline_results)
    report = generate_regression_report(comparison)
    if comparison["has_regressions"]:
        print(report)
        # Fail CI, post PR comment, etc.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from .registry import BENCHMARK_METRICS, BenchmarkMetric, MetricDirection

logger = logging.getLogger(__name__)

__all__ = [
    "Regression",
    "RegressionSeverity",
    "compare_results",
    "generate_regression_report",
    "compare_single_metric",
    "RegressionReport",
]


# ---------------------------------------------------------------------------
# Regression severity levels
# ---------------------------------------------------------------------------

class RegressionSeverity(str):
    """Severity level for a detected regression."""
    NONE = "none"
    WARNING = "warning"
    CRITICAL = "critical"


# ---------------------------------------------------------------------------
# Regression — a single detected regression
# ---------------------------------------------------------------------------

@dataclass
class Regression:
    """A single detected regression.

    Attributes
    ----------
    protein_name : str
        Gene symbol.
    organism_name : str
        Target organism.
    metric_name : str
        Name of the metric that regressed.
    baseline_value : float
        Value in the baseline (previous release).
    current_value : float
        Value in the current run.
    change_fraction : float
        Fractional change from baseline (negative = improvement for
        higher-is-better, positive = regression).
    threshold : float
        Regression threshold that was exceeded.
    severity : str
        ``"warning"`` or ``"critical"``.
    message : str
        Human-readable description.
    """

    protein_name: str
    organism_name: str
    metric_name: str
    baseline_value: float
    current_value: float
    change_fraction: float
    threshold: float
    severity: str = RegressionSeverity.WARNING
    message: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "protein_name": self.protein_name,
            "organism_name": self.organism_name,
            "metric_name": self.metric_name,
            "baseline_value": round(self.baseline_value, 6),
            "current_value": round(self.current_value, 6),
            "change_fraction": round(self.change_fraction, 4),
            "threshold": self.threshold,
            "severity": self.severity,
            "message": self.message,
        }


# ---------------------------------------------------------------------------
# compare_single_metric — compare one metric between baseline and current
# ---------------------------------------------------------------------------

def compare_single_metric(
    metric_name: str,
    baseline_value: float,
    current_value: float,
    metric_definition: BenchmarkMetric | None = None,
) -> Regression | None:
    """Compare a single metric value against the baseline.

    Parameters
    ----------
    metric_name : str
        Name of the metric.
    baseline_value : float
        Value in the baseline.
    current_value : float
        Value in the current run.
    metric_definition : BenchmarkMetric | None
        Metric definition from the registry, or ``None`` to use defaults.

    Returns
    -------
    Regression | None
        A Regression if the change exceeds the threshold, else ``None``.
        Note: protein_name and organism_name must be set by the caller.
    """
    if metric_definition is None:
        metric_definition = BENCHMARK_METRICS.get(metric_name)

    # Default thresholds if metric is not in registry
    threshold = 0.10  # 10% default
    direction = MetricDirection.HIGHER_IS_BETTER
    target_value = None

    if metric_definition is not None:
        threshold = metric_definition.regression_threshold
        direction = metric_definition.direction
        target_value = metric_definition.target_value

    # Skip comparison if baseline is zero or very small
    if abs(baseline_value) < 1e-10:
        return None

    if direction == MetricDirection.HIGHER_IS_BETTER:
        # Regression if current is significantly lower
        change = (current_value - baseline_value) / abs(baseline_value)
        if change < -threshold:
            severity = RegressionSeverity.CRITICAL if abs(change) > 2 * threshold else RegressionSeverity.WARNING
            return Regression(
                protein_name="",  # Set by caller
                organism_name="",  # Set by caller
                metric_name=metric_name,
                baseline_value=baseline_value,
                current_value=current_value,
                change_fraction=change,
                threshold=threshold,
                severity=severity,
                message=(
                    f"{metric_name}: {current_value:.4f} vs baseline {baseline_value:.4f} "
                    f"({change:+.1%}, threshold: -{threshold:.0%})"
                ),
            )

    elif direction == MetricDirection.LOWER_IS_BETTER:
        # Regression if current is significantly higher
        change = (current_value - baseline_value) / abs(baseline_value)
        if change > threshold:
            severity = RegressionSeverity.CRITICAL if abs(change) > 2 * threshold else RegressionSeverity.WARNING
            return Regression(
                protein_name="",
                organism_name="",
                metric_name=metric_name,
                baseline_value=baseline_value,
                current_value=current_value,
                change_fraction=change,
                threshold=threshold,
                severity=severity,
                message=(
                    f"{metric_name}: {current_value:.4f} vs baseline {baseline_value:.4f} "
                    f"({change:+.1%}, threshold: +{threshold:.0%})"
                ),
            )

    elif direction == MetricDirection.CLOSER_TO_TARGET:
        # Regression if current is further from target than baseline
        if target_value is None:
            return None
        baseline_deviation = abs(baseline_value - target_value)
        current_deviation = abs(current_value - target_value)
        if baseline_deviation < 1e-10:
            return None
        deviation_change = (current_deviation - baseline_deviation) / baseline_deviation
        if deviation_change > threshold:
            severity = RegressionSeverity.CRITICAL if deviation_change > 2 * threshold else RegressionSeverity.WARNING
            return Regression(
                protein_name="",
                organism_name="",
                metric_name=metric_name,
                baseline_value=baseline_value,
                current_value=current_value,
                change_fraction=deviation_change,
                threshold=threshold,
                severity=severity,
                message=(
                    f"{metric_name}: {current_value:.4f} vs baseline {baseline_value:.4f} "
                    f"(deviation from target {target_value}: {current_deviation:.4f} vs "
                    f"{baseline_deviation:.4f}, +{deviation_change:.1%})"
                ),
            )

    return None


# ---------------------------------------------------------------------------
# compare_results — compare full benchmark run against baseline
# ---------------------------------------------------------------------------

def compare_results(
    current: Any,
    baseline: Any,
) -> dict[str, Any]:
    """Compare current benchmark results against a baseline.

    Parameters
    ----------
    current : BenchmarkRunResult
        Current benchmark results (from ``perf_runner``).
    baseline : BenchmarkRunResult
        Baseline results from a previous release.

    Returns
    -------
    dict
        Comparison report with:
        - ``"has_regressions"``: bool
        - ``"regressions"``: list[dict]
        - ``"improvements"``: list[dict]
        - ``"unchanged"``: list[dict]
        - ``"summary"``: dict
    """
    from .perf_runner import BenchmarkRunResult

    # Build lookup maps: (protein, organism) -> result
    current_map = {
        (r.protein_name, r.organism_name): r
        for r in current.results
    }
    baseline_map = {
        (r.protein_name, r.organism_name): r
        for r in baseline.results
    }

    # Find common keys
    common_keys = set(current_map.keys()) & set(baseline_map.keys())
    new_keys = set(current_map.keys()) - set(baseline_map.keys())
    removed_keys = set(baseline_map.keys()) - set(current_map.keys())

    regressions: list[Regression] = []
    improvements: list[dict[str, Any]] = []
    unchanged: list[dict[str, Any]] = []

    # Metrics to compare
    compare_metrics = [
        ("cai", "cai"),
        ("gc_mean", "gc_mean"),
        ("optimization_time_s", "optimization_time"),
        ("constraint_satisfaction_rate", "constraint_satisfaction_rate"),
        ("restriction_site_total", "restriction_site_total"),
        ("cryptic_splice_sites", "cryptic_splice_sites"),
        ("cpg_islands", "cpg_islands"),
        ("mrna_stability", "mrna_stability"),
    ]

    for key in sorted(common_keys):
        protein_name, organism_name = key
        curr = current_map[key]
        base = baseline_map[key]

        # Skip failed runs
        if curr.error is not None or base.error is not None:
            continue

        for attr_name, metric_name in compare_metrics:
            curr_val = getattr(curr, attr_name, None)
            base_val = getattr(base, attr_name, None)

            if curr_val is None or base_val is None:
                continue

            reg = compare_single_metric(metric_name, base_val, curr_val)
            if reg is not None:
                reg.protein_name = protein_name
                reg.organism_name = organism_name
                regressions.append(reg)
            else:
                # Check for improvement
                change = (curr_val - base_val) / abs(base_val) if abs(base_val) > 1e-10 else 0.0
                metric_def = BENCHMARK_METRICS.get(metric_name)
                if metric_def is not None:
                    direction = metric_def.direction
                    if direction == MetricDirection.HIGHER_IS_BETTER and change > metric_def.regression_threshold:
                        improvements.append({
                            "protein_name": protein_name,
                            "organism_name": organism_name,
                            "metric_name": metric_name,
                            "baseline_value": base_val,
                            "current_value": curr_val,
                            "change_fraction": change,
                            "message": f"{metric_name}: {curr_val:.4f} vs {base_val:.4f} ({change:+.1%})",
                        })
                    elif direction == MetricDirection.LOWER_IS_BETTER and change < -metric_def.regression_threshold:
                        improvements.append({
                            "protein_name": protein_name,
                            "organism_name": organism_name,
                            "metric_name": metric_name,
                            "baseline_value": base_val,
                            "current_value": curr_val,
                            "change_fraction": change,
                            "message": f"{metric_name}: {curr_val:.4f} vs {base_val:.4f} ({change:+.1%})",
                        })
                    else:
                        unchanged.append({
                            "protein_name": protein_name,
                            "organism_name": organism_name,
                            "metric_name": metric_name,
                            "change_fraction": change,
                        })

    # Summary
    critical_count = sum(1 for r in regressions if r.severity == RegressionSeverity.CRITICAL)
    warning_count = sum(1 for r in regressions if r.severity == RegressionSeverity.WARNING)

    summary = {
        "current_version": current.version,
        "baseline_version": baseline.version,
        "common_benchmarks": len(common_keys),
        "new_benchmarks": len(new_keys),
        "removed_benchmarks": len(removed_keys),
        "total_regressions": len(regressions),
        "critical_regressions": critical_count,
        "warning_regressions": warning_count,
        "improvements": len(improvements),
        "unchanged_metrics": len(unchanged),
        "has_regressions": len(regressions) > 0,
        "has_critical_regressions": critical_count > 0,
    }

    return {
        "has_regressions": len(regressions) > 0,
        "regressions": [r.to_dict() for r in regressions],
        "improvements": improvements,
        "unchanged": unchanged,
        "summary": summary,
        "new_benchmarks": [f"{p}/{o}" for p, o in new_keys],
        "removed_benchmarks": [f"{p}/{o}" for p, o in removed_keys],
    }


# ---------------------------------------------------------------------------
# RegressionReport — structured report for CI / release notes
# ---------------------------------------------------------------------------

@dataclass
class RegressionReport:
    """Structured regression report suitable for CI comments and release notes.

    Attributes
    ----------
    has_regressions : bool
        Whether any regressions were detected.
    has_critical : bool
        Whether any critical regressions were detected.
    text : str
        Human-readable markdown report.
    comparison : dict
        Full comparison data.
    """

    has_regressions: bool
    has_critical: bool
    text: str
    comparison: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "has_regressions": self.has_regressions,
            "has_critical": self.has_critical,
            "text": self.text,
            "comparison": self.comparison,
        }


def generate_regression_report(comparison: dict[str, Any]) -> RegressionReport:
    """Generate a human-readable regression report.

    Parameters
    ----------
    comparison : dict
        Output from :func:`compare_results`.

    Returns
    -------
    RegressionReport
        Structured report with markdown text.
    """
    summary = comparison.get("summary", {})
    regressions = comparison.get("regressions", [])
    improvements = comparison.get("improvements", [])

    lines: list[str] = []
    lines.append("## BioCompiler Benchmark Regression Report")
    lines.append("")

    # Header
    lines.append(f"| Metric | Value |")
    lines.append(f"|--------|-------|")
    lines.append(f"| Current version | {summary.get('current_version', 'unknown')} |")
    lines.append(f"| Baseline version | {summary.get('baseline_version', 'unknown')} |")
    lines.append(f"| Common benchmarks | {summary.get('common_benchmarks', 0)} |")
    lines.append(f"| New benchmarks | {summary.get('new_benchmarks', 0)} |")
    lines.append(f"| **Total regressions** | **{summary.get('total_regressions', 0)}** |")
    lines.append(f"| Critical regressions | {summary.get('critical_regressions', 0)} |")
    lines.append(f"| Warning regressions | {summary.get('warning_regressions', 0)} |")
    lines.append(f"| Improvements | {summary.get('improvements', 0)} |")
    lines.append("")

    # Regressions
    if regressions:
        lines.append("### :rotating_light: Regressions")
        lines.append("")
        lines.append("| Protein | Organism | Metric | Baseline | Current | Change | Severity |")
        lines.append("|---------|----------|--------|----------|---------|--------|----------|")
        for r in sorted(regressions, key=lambda x: (x.get("severity", ""), x.get("metric_name", ""))):
            severity_icon = ":red_circle:" if r.get("severity") == "critical" else ":yellow_circle:"
            lines.append(
                f"| {r['protein_name']} | {r['organism_name']} | {r['metric_name']} | "
                f"{r['baseline_value']:.4f} | {r['current_value']:.4f} | "
                f"{r['change_fraction']:+.1%} | {severity_icon} {r['severity']} |"
            )
        lines.append("")

    # Improvements
    if improvements:
        lines.append("### :white_check_mark: Improvements")
        lines.append("")
        lines.append("| Protein | Organism | Metric | Baseline | Current | Change |")
        lines.append("|---------|----------|--------|----------|---------|--------|")
        for imp in sorted(improvements, key=lambda x: abs(x.get("change_fraction", 0)), reverse=True)[:10]:
            lines.append(
                f"| {imp['protein_name']} | {imp['organism_name']} | {imp['metric_name']} | "
                f"{imp['baseline_value']:.4f} | {imp['current_value']:.4f} | "
                f"{imp['change_fraction']:+.1%} |"
            )
        if len(improvements) > 10:
            lines.append(f"| ... | ... | ... | ... | ... | _{len(improvements) - 10} more_ |")
        lines.append("")

    # New / removed benchmarks
    new_benchmarks = comparison.get("new_benchmarks", [])
    removed_benchmarks = comparison.get("removed_benchmarks", [])
    if new_benchmarks:
        lines.append(f"### New benchmarks: {', '.join(new_benchmarks)}")
        lines.append("")
    if removed_benchmarks:
        lines.append(f"### Removed benchmarks: {', '.join(removed_benchmarks)}")
        lines.append("")

    # Verdict
    if comparison.get("has_regressions"):
        lines.append("### :x: VERDICT: Regressions detected")
        lines.append("")
        lines.append("This PR introduces performance or quality regressions.")
        lines.append("Please review the table above and address before merging.")
    else:
        lines.append("### :white_check_mark: VERDICT: No regressions detected")
        lines.append("")
        lines.append("All benchmarks are within acceptable thresholds compared to baseline.")

    text = "\n".join(lines)

    return RegressionReport(
        has_regressions=comparison.get("has_regressions", False),
        has_critical=summary.get("has_critical_regressions", False),
        text=text,
        comparison=comparison,
    )
