"""
BioCompiler Wet-Lab Validation Framework

This module provides tools for validating computational predictions against
experimental (wet-lab) data. It enables closing the loop between in-silico
gene design optimization and actual experimental measurements.

Key classes:
  - ExperimentalResult: A single experimental data point
  - ValidationComparison: Comparison between prediction and experiment
  - WetLabValidator: Orchestrates validation across multiple data points

Statistical methods:
  - Pearson correlation (linear relationship between predicted CAI and measured expression)
  - Spearman rank correlation (monotonic relationship, robust to outliers)
  - Rank-order matching (whether optimization predictions preserve experimental ordering)

Usage:
    validator = WetLabValidator()
    validator.add_experimental_result(ExperimentalResult(
        gene_name="INS",
        organism="Escherichia_coli",
        measured_expression_level=1.2e6,
        measured_cai=0.95,
        sequence_used="ATGGCT...",
        notes="Shake flask, 37°C, 16h",
    ))
    # ... add more results ...

    # Compare with optimization predictions
    comparison = validator.compare_with_prediction(opt_result, "INS")

    # Get overall correlation
    r, p = validator.compute_correlation()
    rho, p_rho = validator.compute_rank_correlation()

    # Generate validation report
    report = validator.validation_report()
    validator.save_report("validation_report.json")
"""

from __future__ import annotations

import csv
import json
import math
import logging
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .optimization import OptimizationResult

logger = logging.getLogger(__name__)


@dataclass
class ExperimentalResult:
    """A single wet-lab experimental data point.

    Attributes:
        gene_name: Name/identifier of the gene (e.g., "INS", "GFP").
        organism: Target organism (canonical name, e.g., "Escherichia_coli").
        measured_expression_level: Measured expression level in arbitrary units
            (e.g., fluorescence intensity, protein concentration in mg/L).
        measured_cai: Measured codon adaptation index (if available from sequencing).
            This may be None if CAI was not measured experimentally.
        sequence_used: DNA sequence used in the experiment.
        notes: Free-text notes about experimental conditions.
    """
    gene_name: str
    organism: str
    measured_expression_level: float
    measured_cai: Optional[float] = None
    sequence_used: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a plain dict."""
        d = asdict(self)
        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ExperimentalResult":
        """Deserialize from a plain dict."""
        return cls(
            gene_name=data["gene_name"],
            organism=data["organism"],
            measured_expression_level=float(data["measured_expression_level"]),
            measured_cai=float(data["measured_cai"]) if data.get("measured_cai") is not None else None,
            sequence_used=data.get("sequence_used", ""),
            notes=data.get("notes", ""),
        )


@dataclass
class ValidationComparison:
    """Comparison between a computational prediction and experimental result.

    Attributes:
        gene_name: Name/identifier of the gene.
        predicted_cai: CAI predicted by the optimizer.
        measured_expression: Measured expression level from the experiment.
        correlation: Pearson correlation between predicted CAI and measured expression
            for this data point (1.0 if perfectly correlated, -1.0 if anti-correlated).
            This is the per-point contribution to the overall correlation.
        rank_order_match: Whether the rank ordering of this gene matches between
            prediction and experiment. True if the gene's rank by CAI equals its
            rank by measured expression.
    """
    gene_name: str
    predicted_cai: float
    measured_expression: float
    correlation: float
    rank_order_match: bool

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a plain dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ValidationComparison":
        """Deserialize from a plain dict."""
        return cls(
            gene_name=data["gene_name"],
            predicted_cai=float(data["predicted_cai"]),
            measured_expression=float(data["measured_expression"]),
            correlation=float(data["correlation"]),
            rank_order_match=bool(data["rank_order_match"]),
        )


class WetLabValidator:
    """Validates computational predictions against experimental data.

    This class manages a collection of experimental results and provides
    methods to compare them with optimization predictions, compute
    statistical correlations, and generate validation reports.

    Usage:
        validator = WetLabValidator()
        validator.add_experimental_result(result)
        comparison = validator.compare_with_prediction(opt_result, "INS")
        r, p = validator.compute_correlation()
        report = validator.validation_report()
    """

    def __init__(self) -> None:
        self._results: list[ExperimentalResult] = []
        self._comparisons: list[ValidationComparison] = []

    @property
    def results(self) -> list[ExperimentalResult]:
        """Return a copy of the stored experimental results."""
        return list(self._results)

    @property
    def comparisons(self) -> list[ValidationComparison]:
        """Return a copy of the stored comparisons."""
        return list(self._comparisons)

    def add_experimental_result(self, result: ExperimentalResult) -> None:
        """Add an experimental data point to the validator.

        Args:
            result: The experimental result to add.

        Raises:
            ValueError: If a result with the same gene_name and organism
                already exists.
        """
        # Check for duplicates
        for existing in self._results:
            if existing.gene_name == result.gene_name and existing.organism == result.organism:
                raise ValueError(
                    f"Experimental result for gene '{result.gene_name}' "
                    f"in organism '{result.organism}' already exists. "
                    f"Remove the existing result before adding a new one."
                )
        self._results.append(result)

    def remove_experimental_result(self, gene_name: str, organism: str) -> bool:
        """Remove an experimental result by gene_name and organism.

        Returns:
            True if a result was removed, False if not found.
        """
        for i, r in enumerate(self._results):
            if r.gene_name == gene_name and r.organism == organism:
                self._results.pop(i)
                # Also remove associated comparisons
                self._comparisons = [
                    c for c in self._comparisons if c.gene_name != gene_name
                ]
                return True
        return False

    def compare_with_prediction(
        self,
        optimization_result: "OptimizationResult",
        gene_name: str,
    ) -> ValidationComparison:
        """Compare an optimization prediction with experimental data.

        Args:
            optimization_result: The OptimizationResult from BioCompiler.
            gene_name: Gene name to match with experimental data.

        Returns:
            A ValidationComparison with correlation and rank-order information.

        Raises:
            ValueError: If no experimental result is found for the gene.
        """
        # Find the matching experimental result
        exp_result = None
        for r in self._results:
            if r.gene_name == gene_name:
                exp_result = r
                break

        if exp_result is None:
            raise ValueError(
                f"No experimental result found for gene '{gene_name}'. "
                f"Available genes: {[r.gene_name for r in self._results]}"
            )

        predicted_cai = optimization_result.cai
        measured_expression = exp_result.measured_expression_level

        # Compute per-point correlation contribution
        # For a single point, we use a simple signed deviation measure
        # The actual correlation is computed across all data points
        correlation = self._compute_point_correlation(predicted_cai, measured_expression)

        # Compute rank-order match
        rank_order_match = self._compute_rank_order_match(gene_name, predicted_cai, measured_expression)

        comparison = ValidationComparison(
            gene_name=gene_name,
            predicted_cai=predicted_cai,
            measured_expression=measured_expression,
            correlation=correlation,
            rank_order_match=rank_order_match,
        )

        # Update existing comparison or add new one
        existing_idx = None
        for i, c in enumerate(self._comparisons):
            if c.gene_name == gene_name:
                existing_idx = i
                break

        if existing_idx is not None:
            self._comparisons[existing_idx] = comparison
        else:
            self._comparisons.append(comparison)

        return comparison

    def _compute_point_correlation(self, predicted_cai: float, measured_expression: float) -> float:
        """Compute a per-point correlation contribution.

        For a single point, we can't compute a meaningful Pearson r.
        Instead, we compute the z-score of this point relative to all
        comparisons, or return 1.0 if this is the first point.
        """
        if len(self._comparisons) < 1:
            return 1.0  # First point, trivially correlated

        # Check if this point is consistent with the trend
        existing_cais = [c.predicted_cai for c in self._comparisons]
        existing_exprs = [c.measured_expression for c in self._comparisons]

        # Simple heuristic: if higher CAI tends to give higher expression
        # among existing points, check if this new point follows the trend
        if len(existing_cais) < 2:
            return 1.0

        # Compute slope of best-fit line for existing points
        mean_cai = sum(existing_cais) / len(existing_cais)
        mean_expr = sum(existing_exprs) / len(existing_exprs)

        numerator = sum(
            (cai - mean_cai) * (expr - mean_expr)
            for cai, expr in zip(existing_cais, existing_exprs)
        )
        denominator_cai = sum((cai - mean_cai) ** 2 for cai in existing_cais)
        denominator_expr = sum((expr - mean_expr) ** 2 for expr in existing_exprs)

        if denominator_cai == 0 or denominator_expr == 0:
            return 0.0

        slope = numerator / denominator_cai

        # Check if the new point follows the trend
        if slope > 0:
            # Positive trend: higher CAI → higher expression
            if predicted_cai > mean_cai and measured_expression > mean_expr:
                return 1.0
            elif predicted_cai < mean_cai and measured_expression < mean_expr:
                return 1.0
            else:
                return -1.0
        elif slope < 0:
            # Negative trend: higher CAI → lower expression
            if predicted_cai > mean_cai and measured_expression < mean_expr:
                return 1.0
            elif predicted_cai < mean_cai and measured_expression > mean_expr:
                return 1.0
            else:
                return -1.0
        else:
            return 0.0

    def _compute_rank_order_match(
        self, gene_name: str, predicted_cai: float, measured_expression: float,
    ) -> bool:
        """Check if the rank ordering of this gene matches between prediction and experiment.

        A gene's rank order matches if, among all compared genes, the ordering
        by predicted CAI is the same as the ordering by measured expression.
        """
        if len(self._comparisons) < 1:
            return True  # First point, trivially matches

        # Collect all data points including the new one
        all_cais = [(c.gene_name, c.predicted_cai) for c in self._comparisons]
        all_cais.append((gene_name, predicted_cai))
        all_exprs = [(c.gene_name, c.measured_expression) for c in self._comparisons]
        all_exprs.append((gene_name, measured_expression))

        # Sort by CAI and expression separately
        cai_ranks = {name: rank for rank, (name, _) in enumerate(
            sorted(all_cais, key=lambda x: x[1], reverse=True)
        )}
        expr_ranks = {name: rank for rank, (name, _) in enumerate(
            sorted(all_exprs, key=lambda x: x[1], reverse=True)
        )}

        # Check if this gene's rank matches
        return cai_ranks.get(gene_name, -1) == expr_ranks.get(gene_name, -1)

    def compute_correlation(self) -> tuple[float, float]:
        """Compute Pearson correlation between predicted CAI and measured expression.

        Returns:
            Tuple of (Pearson r, p-value). If fewer than 3 data points are
            available, returns (0.0, 1.0).

        Raises:
            ValueError: If no comparisons are available.
        """
        if not self._comparisons:
            raise ValueError("No comparisons available. Call compare_with_prediction first.")

        if len(self._comparisons) < 3:
            # Not enough data points for meaningful correlation
            if len(self._comparisons) == 1:
                return (1.0, 1.0)  # Trivially correlated
            # Two points: compute simple correlation
            c0, c1 = self._comparisons[0], self._comparisons[1]
            if c0.predicted_cai == c1.predicted_cai or c0.measured_expression == c1.measured_expression:
                return (0.0, 1.0)
            # Two points: r is either +1 or -1
            direction = (
                (c1.predicted_cai - c0.predicted_cai) *
                (c1.measured_expression - c0.measured_expression)
            )
            r = 1.0 if direction > 0 else -1.0
            return (r, 1.0)

        # Compute Pearson r from scratch (no scipy dependency)
        cais = [c.predicted_cai for c in self._comparisons]
        exprs = [c.measured_expression for c in self._comparisons]
        n = len(cais)

        mean_cai = sum(cais) / n
        mean_expr = sum(exprs) / n

        numerator = sum(
            (cai - mean_cai) * (expr - mean_expr)
            for cai, expr in zip(cais, exprs)
        )
        ss_cai = sum((cai - mean_cai) ** 2 for cai in cais)
        ss_expr = sum((expr - mean_expr) ** 2 for expr in exprs)

        if ss_cai == 0 or ss_expr == 0:
            return (0.0, 1.0)

        r = numerator / math.sqrt(ss_cai * ss_expr)

        # Clamp to [-1, 1] for floating point safety
        r = max(-1.0, min(1.0, r))

        # Approximate p-value using t-distribution
        # t = r * sqrt((n-2) / (1 - r^2))
        if abs(r) >= 1.0:
            p_value = 0.0
        else:
            t_stat = r * math.sqrt((n - 2) / (1 - r ** 2))
            # Approximate two-tailed p-value
            # Using the incomplete beta function approximation
            p_value = self._approximate_p_value(t_stat, n - 2)

        return (r, p_value)

    def _approximate_p_value(self, t_stat: float, df: int) -> float:
        """Approximate two-tailed p-value from t-statistic using normal approximation.

        For df >= 30, the t-distribution is well-approximated by the normal.
        For small df, we use a simple correction.
        """
        if df <= 0:
            return 1.0

        # Normal approximation to t-distribution
        x = abs(t_stat) / math.sqrt(df)

        # Simple approximation using the error function
        # P(T > |t|) ≈ 2 * (1 - Φ(|t| * correction))
        # For the normal CDF, we use the Abramowitz & Stegun approximation
        z = abs(t_stat)

        # For small df, apply a correction
        if df < 30:
            correction = df - 2
            if correction <= 0:
                # For very small df, just use the normal approximation without correction
                pass
            else:
                z = z * math.sqrt(df / correction)

        # Normal CDF approximation (Abramowitz & Stegun)
        if z > 10:
            return 0.0

        # Standard normal survival function approximation
        t_val = 1.0 / (1.0 + 0.2316419 * z)
        d = 0.3989422804014327  # 1/sqrt(2*pi)
        poly = t_val * (
            0.319381530 +
            t_val * (-0.356563782 +
            t_val * (1.781477937 +
            t_val * (-1.821255978 +
            t_val * 1.330274429)))
        )
        p_one_tail = d * math.exp(-0.5 * z * z) * poly
        p_value = 2.0 * max(0.0, min(1.0, p_one_tail))

        return p_value

    def compute_rank_correlation(self) -> tuple[float, float]:
        """Compute Spearman rank correlation between predicted CAI and measured expression.

        Returns:
            Tuple of (Spearman rho, p-value). If fewer than 3 data points,
            returns (0.0, 1.0).

        Raises:
            ValueError: If no comparisons are available.
        """
        if not self._comparisons:
            raise ValueError("No comparisons available. Call compare_with_prediction first.")

        if len(self._comparisons) < 3:
            if len(self._comparisons) == 1:
                return (1.0, 1.0)
            c0, c1 = self._comparisons[0], self._comparisons[1]
            if c0.predicted_cai == c1.predicted_cai or c0.measured_expression == c1.measured_expression:
                return (0.0, 1.0)
            direction = (
                (c1.predicted_cai - c0.predicted_cai) *
                (c1.measured_expression - c0.measured_expression)
            )
            rho = 1.0 if direction > 0 else -1.0
            return (rho, 1.0)

        cais = [c.predicted_cai for c in self._comparisons]
        exprs = [c.measured_expression for c in self._comparisons]
        n = len(cais)

        # Compute ranks
        cai_ranks = self._compute_ranks(cais)
        expr_ranks = self._compute_ranks(exprs)

        # Spearman rho = 1 - 6 * sum(d_i^2) / (n * (n^2 - 1))
        d_squared_sum = sum(
            (cr - er) ** 2
            for cr, er in zip(cai_ranks, expr_ranks)
        )

        denominator = n * (n ** 2 - 1)
        if denominator == 0:
            return (0.0, 1.0)

        rho = 1.0 - (6.0 * d_squared_sum) / denominator

        # Clamp
        rho = max(-1.0, min(1.0, rho))

        # Approximate p-value using t-distribution
        if abs(rho) >= 1.0:
            p_value = 0.0
        else:
            t_stat = rho * math.sqrt((n - 2) / (1 - rho ** 2))
            p_value = self._approximate_p_value(t_stat, n - 2)

        return (rho, p_value)

    def _compute_ranks(self, values: list[float]) -> list[float]:
        """Compute ranks for a list of values, handling ties with average rank."""
        n = len(values)
        indexed = sorted(enumerate(values), key=lambda x: x[1])

        ranks = [0.0] * n
        i = 0
        while i < n:
            # Find all tied values
            j = i
            while j < n and indexed[j][1] == indexed[i][1]:
                j += 1
            # Average rank for tied values
            avg_rank = (i + j - 1) / 2.0 + 1.0  # 1-based
            for k in range(i, j):
                ranks[indexed[k][0]] = avg_rank
            i = j

        return ranks

    def validation_report(self) -> dict[str, Any]:
        """Generate a comprehensive validation report.

        Returns:
            Dict with validation statistics and per-gene comparisons.

        Raises:
            ValueError: If no comparisons are available.
        """
        if not self._comparisons:
            return {
                "status": "no_data",
                "message": "No comparisons available. Add experimental results and compare with predictions.",
                "num_results": len(self._results),
                "num_comparisons": 0,
            }

        cais = [c.predicted_cai for c in self._comparisons]
        exprs = [c.measured_expression for c in self._comparisons]

        # Compute correlations
        try:
            pearson_r, pearson_p = self.compute_correlation()
        except ValueError:
            pearson_r, pearson_p = 0.0, 1.0

        try:
            spearman_rho, spearman_p = self.compute_rank_correlation()
        except ValueError:
            spearman_rho, spearman_p = 0.0, 1.0

        # Rank order match percentage
        rank_matches = sum(1 for c in self._comparisons if c.rank_order_match)
        rank_match_pct = rank_matches / len(self._comparisons) * 100 if self._comparisons else 0

        # Summary statistics
        report: dict[str, Any] = {
            "status": "ok",
            "num_results": len(self._results),
            "num_comparisons": len(self._comparisons),
            "pearson_correlation": {
                "r": pearson_r,
                "p_value": pearson_p,
                "significant_at_005": pearson_p < 0.05,
            },
            "spearman_rank_correlation": {
                "rho": spearman_rho,
                "p_value": spearman_p,
                "significant_at_005": spearman_p < 0.05,
            },
            "rank_order_match_pct": rank_match_pct,
            "predicted_cai_stats": {
                "mean": sum(cais) / len(cais) if cais else 0,
                "min": min(cais) if cais else 0,
                "max": max(cais) if cais else 0,
            },
            "measured_expression_stats": {
                "mean": sum(exprs) / len(exprs) if exprs else 0,
                "min": min(exprs) if exprs else 0,
                "max": max(exprs) if exprs else 0,
            },
            "comparisons": [c.to_dict() for c in self._comparisons],
            "experimental_results": [r.to_dict() for r in self._results],
        }

        return report

    def load_from_csv(self, path: str) -> int:
        """Load experimental results from a CSV file.

        Expected columns:
            gene_name, organism, measured_expression_level,
            measured_cai (optional), sequence_used (optional), notes (optional)

        Args:
            path: Path to the CSV file.

        Returns:
            Number of results loaded.

        Raises:
            FileNotFoundError: If the file doesn't exist.
            ValueError: If required columns are missing.
        """
        csv_path = Path(path)
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {path}")

        required_columns = {"gene_name", "organism", "measured_expression_level"}
        count = 0

        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)

            # Validate columns
            if reader.fieldnames is None:
                raise ValueError("CSV file has no header row")
            missing = required_columns - set(reader.fieldnames)
            if missing:
                raise ValueError(f"Missing required columns: {missing}")

            for row in reader:
                try:
                    result = ExperimentalResult(
                        gene_name=row["gene_name"],
                        organism=row["organism"],
                        measured_expression_level=float(row["measured_expression_level"]),
                        measured_cai=(
                            float(row["measured_cai"])
                            if row.get("measured_cai") and row["measured_cai"].strip()
                            else None
                        ),
                        sequence_used=row.get("sequence_used", ""),
                        notes=row.get("notes", ""),
                    )
                    self.add_experimental_result(result)
                    count += 1
                except (ValueError, KeyError) as e:
                    logger.warning(f"Skipping row {count + 1}: {e}")
                    continue

        return count

    def save_report(self, path: str) -> None:
        """Save the validation report as a JSON file.

        Args:
            path: Path to save the JSON report.
        """
        report = self.validation_report()
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, default=str)


__all__ = [
    "ExperimentalResult",
    "ValidationComparison",
    "WetLabValidator",
]

# Global validator instance for API usage.
# In production, this would be replaced with a database-backed store.
_global_validator = WetLabValidator()
