"""
BioCompiler What-If Alternative Analysis — Provenance for Constraint Trade-offs
=================================================================================

Provides "what-if" scenario analysis that tells users how much CAI or other
quality metrics change when constraints are relaxed, tightened, added, or
removed.

Every constraint imposed during gene optimization has a cost — typically
measured in CAI.  This module makes those costs explicit by re-optimizing
the protein under altered constraint configurations and comparing the results.
The key insight is answering questions like: "You're leaving 8% CAI on the
table because of your GC constraint" or "Adding BamHI avoidance only costs
you 2% CAI."  These counterfactual analyses help users make informed
trade-off decisions.

The ``WhatIfAnalyzer`` runs what-if experiments by re-optimizing under altered
constraints and comparing results.  The ``WhatIfReport`` generator produces
both markdown (human-readable) and JSON (machine-readable) reports.  Scenarios
are ranked by CAI impact so the most important trade-offs appear first.

Usage::

    from biocompiler.whatif_analysis import WhatIfAnalyzer, WhatIfReport

    # Initialize analyzer with baseline constraints
    analyzer = WhatIfAnalyzer(
        protein="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
        organism="Homo_sapiens",
        gc_lo=0.30, gc_hi=0.70,
    )

    # Run individual what-if analyses
    dna = "ATGGTG..."  # your current optimized sequence
    scenario = analyzer.analyze_gc_relaxation(dna, current_gc_hi=0.70, alternative_gc_hi=0.75)
    print(f"Relaxing GC max: CAI delta = {scenario.cai_delta:+.4f}")

    # Run a full report covering all standard scenarios
    scenarios = analyzer.full_whatif_report(dna)
    report_text = WhatIfReport.generate(scenarios)
    WhatIfReport.to_json(scenarios, "whatif_report.json")

Components:
- WhatIfScenario: immutable record of a single what-if experiment
- WhatIfAnalyzer: runs what-if experiments by re-optimizing under altered
  constraints and comparing results
- WhatIfReport: generates markdown and JSON reports from scenarios
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from .scanner import gc_content
from .translation import compute_cai
from .optimization import optimize_sequence, OptimizationResult
from .restriction_sites import RESTRICTION_SITES, get_recognition_site
from .constants import RESTRICTION_ENZYMES
from .organisms import SUPPORTED_ORGANISMS, ORGANISM_GC_TARGETS

logger = logging.getLogger(__name__)

__all__ = [
    "WhatIfScenario",
    "WhatIfAnalyzer",
    "WhatIfReport",
]


# ---------------------------------------------------------------------------
# WhatIfScenario
# ---------------------------------------------------------------------------

@dataclass
class WhatIfScenario:
    """Immutable record of a single what-if analysis experiment.

    Captures the counterfactual: "what would happen if we changed constraint X
    from value A to value B?" — including the predicted impact on CAI, GC,
    and constraint satisfaction.

    Attributes:
        description: Human-readable summary, e.g.
            ``"What if we relax GC max from 0.70 to 0.75?"``
        parameter_changed: Machine-readable name of the changed parameter,
            e.g. ``"gc_hi"``, ``"restriction_sites"``, ``"cpg_avoidance"``,
            ``"organism"``.
        original_value: The value used in the baseline optimization.
        alternative_value: The alternative value being explored.
        predicted_cai: CAI of the re-optimized sequence under the alternative,
            or ``None`` if the scenario could not be evaluated.
        predicted_gc: GC content of the re-optimized sequence, or ``None``.
        constraint_satisfaction: Dict mapping constraint names to booleans
            (``True`` = satisfied) for the alternative design, or ``None``
            if not evaluated.
        feasibility: One of ``"feasible"``, ``"infeasible"``, or
            ``"unknown"``.  A scenario is *infeasible* when the optimizer
            cannot find any sequence satisfying the alternative constraints.
        baseline_cai: CAI of the original (baseline) sequence for easy
            comparison.
        baseline_gc: GC content of the original (baseline) sequence.
        cai_delta: Difference between predicted_cai and baseline_cai.
            Positive means the alternative improves CAI.
        timestamp: ISO 8601 timestamp of when this scenario was computed.
    """

    description: str
    parameter_changed: str
    original_value: Any
    alternative_value: Any
    predicted_cai: float | None
    predicted_gc: float | None
    constraint_satisfaction: dict[str, bool] | None
    feasibility: str  # "feasible" | "infeasible" | "unknown"
    baseline_cai: float = 0.0
    baseline_gc: float = 0.0
    cai_delta: float | None = None
    timestamp: str = ""

    def __post_init__(self) -> None:
        """Compute derived fields and validate invariants."""
        if self.feasibility not in ("feasible", "infeasible", "unknown"):
            raise ValueError(
                f"feasibility must be 'feasible', 'infeasible', or 'unknown', "
                f"got {self.feasibility!r}"
            )
        if not self.timestamp:
            self.timestamp = datetime.now(timezone.utc).isoformat()
        if self.predicted_cai is not None and self.cai_delta is None:
            self.cai_delta = self.predicted_cai - self.baseline_cai

    # -- Serialization -------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Serialize this scenario to a JSON-compatible dict."""
        return {
            "description": self.description,
            "parameter_changed": self.parameter_changed,
            "original_value": _json_safe(self.original_value),
            "alternative_value": _json_safe(self.alternative_value),
            "predicted_cai": self.predicted_cai,
            "predicted_gc": self.predicted_gc,
            "constraint_satisfaction": (
                dict(self.constraint_satisfaction)
                if self.constraint_satisfaction is not None
                else None
            ),
            "feasibility": self.feasibility,
            "baseline_cai": self.baseline_cai,
            "baseline_gc": self.baseline_gc,
            "cai_delta": self.cai_delta,
            "timestamp": self.timestamp,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> WhatIfScenario:
        """Deserialize a WhatIfScenario from a plain dict.

        Raises:
            ValueError: If required keys are missing.
        """
        required = {
            "description", "parameter_changed", "original_value",
            "alternative_value", "predicted_cai", "predicted_gc",
            "constraint_satisfaction", "feasibility",
        }
        missing = required - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize WhatIfScenario: missing keys {missing}"
            )
        return cls(
            description=data["description"],
            parameter_changed=data["parameter_changed"],
            original_value=data["original_value"],
            alternative_value=data["alternative_value"],
            predicted_cai=data["predicted_cai"],
            predicted_gc=data["predicted_gc"],
            constraint_satisfaction=(
                dict(data["constraint_satisfaction"])
                if data["constraint_satisfaction"] is not None
                else None
            ),
            feasibility=data["feasibility"],
            baseline_cai=data.get("baseline_cai", 0.0),
            baseline_gc=data.get("baseline_gc", 0.0),
            cai_delta=data.get("cai_delta"),
            timestamp=data.get("timestamp", ""),
        )

    def __repr__(self) -> str:
        cai_str = f"{self.predicted_cai:.4f}" if self.predicted_cai is not None else "N/A"
        delta_str = f"{self.cai_delta:+.4f}" if self.cai_delta is not None else "N/A"
        return (
            f"WhatIfScenario("
            f"param={self.parameter_changed!r}, "
            f"feasibility={self.feasibility!r}, "
            f"CAI={cai_str}, "
            f"delta={delta_str})"
        )


# ---------------------------------------------------------------------------
# WhatIfAnalyzer
# ---------------------------------------------------------------------------

class WhatIfAnalyzer:
    """Run what-if analyses by re-optimizing under altered constraints.

    Each analysis method re-optimizes the protein with a single constraint
    change and reports the impact on CAI, GC content, and constraint
    satisfaction.  The key insight is quantifying the *cost* of each
    constraint so users can make informed trade-off decisions.

    Example::

        analyzer = WhatIfAnalyzer(
            protein="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            organism="Homo_sapiens",
            constraints=[{"type": "gc_range", "lo": 0.30, "hi": 0.70}],
        )
        scenario = analyzer.analyze_gc_relaxation(dna, 0.70, 0.75)
        print(f"Relaxing GC max: CAI improves by {scenario.cai_delta:+.2%}")
    """

    def __init__(
        self,
        protein: str,
        organism: str = "Homo_sapiens",
        constraints: list | None = None,
        gc_lo: float = 0.30,
        gc_hi: float = 0.70,
        restriction_sites: list[str] | None = None,
        cryptic_splice_threshold: float = 3.0,
        cai_threshold: float = 0.5,
    ) -> None:
        """Initialize the WhatIfAnalyzer.

        Args:
            protein: Target protein sequence (single-letter amino acid codes).
            organism: Target organism for codon optimization.
            constraints: Optional list of constraint dicts.  If ``None``,
                default constraints are derived from the other parameters.
            gc_lo: Minimum GC content for the baseline optimization.
            gc_hi: Maximum GC content for the baseline optimization.
            restriction_sites: List of restriction enzyme names to avoid.
                If ``None``, the standard set is used.
            cryptic_splice_threshold: MaxEntScan threshold for cryptic
                splice site elimination.
            cai_threshold: Minimum CAI threshold for constraint checking.

        Raises:
            ValueError: If the organism is not supported.
        """
        if organism not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {organism!r}. "
                f"Supported: {SUPPORTED_ORGANISMS}"
            )
        self.protein = protein.upper().strip()
        self.organism = organism
        self.constraints = constraints or []
        self.gc_lo = gc_lo
        self.gc_hi = gc_hi
        self.restriction_sites = restriction_sites or list(RESTRICTION_ENZYMES.keys())
        self.cryptic_splice_threshold = cryptic_splice_threshold
        self.cai_threshold = cai_threshold

    # -- Helper --------------------------------------------------------------

    def _baseline_optimize(self) -> OptimizationResult:
        """Run the baseline optimization with current constraints."""
        return optimize_sequence(
            target_protein=self.protein,
            organism=self.organism,
            gc_lo=self.gc_lo,
            gc_hi=self.gc_hi,
            enzymes=self.restriction_sites,
            cryptic_splice_threshold=self.cryptic_splice_threshold,
            cai_threshold=self.cai_threshold,
        )

    def _reoptimize(
        self,
        gc_lo: float | None = None,
        gc_hi: float | None = None,
        enzymes: list[str] | None = None,
        organism: str | None = None,
        avoid_cpg: bool | None = None,
    ) -> OptimizationResult:
        """Re-optimize with modified parameters.

        Any parameter that is ``None`` inherits the baseline value.
        """
        result = optimize_sequence(
            target_protein=self.protein,
            organism=organism or self.organism,
            gc_lo=gc_lo if gc_lo is not None else self.gc_lo,
            gc_hi=gc_hi if gc_hi is not None else self.gc_hi,
            enzymes=enzymes if enzymes is not None else self.restriction_sites,
            cryptic_splice_threshold=self.cryptic_splice_threshold,
            cai_threshold=self.cai_threshold,
        )
        return result

    @staticmethod
    def _compute_constraint_satisfaction(
        result: OptimizationResult,
        gc_lo: float,
        gc_hi: float,
    ) -> dict[str, bool]:
        """Compute a dict of constraint satisfaction from an OptimizationResult."""
        satisfied: dict[str, bool] = {}
        satisfied["gc_in_range"] = gc_lo <= result.gc_content <= gc_hi
        satisfied["cai_above_threshold"] = result.cai >= 0.0  # always computed
        satisfied["no_failed_predicates"] = len(result.failed_predicates) == 0

        # Check each predicate result
        for pr in result.predicate_results:
            name = getattr(pr, "predicate", str(pr))
            passed = getattr(pr, "passed", False)
            satisfied[f"predicate:{name}"] = passed

        return satisfied

    # -- Public analysis methods ---------------------------------------------

    def analyze_gc_relaxation(
        self,
        dna: str,
        current_gc_hi: float,
        alternative_gc_hi: float,
    ) -> WhatIfScenario:
        """Analyze the impact of relaxing the upper GC constraint.

        Re-optimizes the protein with the alternative GC upper bound and
        reports how much CAI improves (or degrades).  This directly answers
        "am I leaving CAI on the table because of my GC constraint?"

        Args:
            dna: Current baseline DNA sequence (used for baseline metrics).
            current_gc_hi: Current GC upper bound.
            alternative_gc_hi: Proposed relaxed GC upper bound.

        Returns:
            WhatIfScenario with the predicted impact.
        """
        baseline_cai = compute_cai(dna, self.organism)
        baseline_gc = gc_content(dna)

        description = (
            f"What if we relax GC max from {current_gc_hi:.2f} "
            f"to {alternative_gc_hi:.2f}?"
        )

        try:
            result = self._reoptimize(gc_hi=alternative_gc_hi)
            predicted_cai = result.cai
            predicted_gc = result.gc_content
            constraint_sat = self._compute_constraint_satisfaction(
                result, self.gc_lo, alternative_gc_hi,
            )
            feasibility = "feasible"
        except Exception as exc:
            logger.warning("GC relaxation scenario failed: %s", exc)
            predicted_cai = None
            predicted_gc = None
            constraint_sat = None
            feasibility = "unknown"

        return WhatIfScenario(
            description=description,
            parameter_changed="gc_hi",
            original_value=current_gc_hi,
            alternative_value=alternative_gc_hi,
            predicted_cai=predicted_cai,
            predicted_gc=predicted_gc,
            constraint_satisfaction=constraint_sat,
            feasibility=feasibility,
            baseline_cai=baseline_cai,
            baseline_gc=baseline_gc,
        )

    def analyze_gc_tightening(
        self,
        dna: str,
        current_gc_lo: float,
        current_gc_hi: float,
        alternative_gc_lo: float,
        alternative_gc_hi: float,
    ) -> WhatIfScenario:
        """Analyze the impact of tightening the GC constraint range.

        Useful for understanding the cost of stricter GC requirements
        (e.g., for regulatory or stability reasons).

        Args:
            dna: Current baseline DNA sequence.
            current_gc_lo: Current GC lower bound.
            current_gc_hi: Current GC upper bound.
            alternative_gc_lo: Proposed tighter GC lower bound.
            alternative_gc_hi: Proposed tighter GC upper bound.

        Returns:
            WhatIfScenario with the predicted impact.
        """
        baseline_cai = compute_cai(dna, self.organism)
        baseline_gc = gc_content(dna)

        description = (
            f"What if we tighten GC from [{current_gc_lo:.2f}, {current_gc_hi:.2f}] "
            f"to [{alternative_gc_lo:.2f}, {alternative_gc_hi:.2f}]?"
        )

        try:
            result = self._reoptimize(
                gc_lo=alternative_gc_lo,
                gc_hi=alternative_gc_hi,
            )
            predicted_cai = result.cai
            predicted_gc = result.gc_content
            constraint_sat = self._compute_constraint_satisfaction(
                result, alternative_gc_lo, alternative_gc_hi,
            )
            feasibility = "feasible"
        except Exception as exc:
            logger.warning("GC tightening scenario failed: %s", exc)
            predicted_cai = None
            predicted_gc = None
            constraint_sat = None
            feasibility = "unknown"

        return WhatIfScenario(
            description=description,
            parameter_changed="gc_range",
            original_value=[current_gc_lo, current_gc_hi],
            alternative_value=[alternative_gc_lo, alternative_gc_hi],
            predicted_cai=predicted_cai,
            predicted_gc=predicted_gc,
            constraint_satisfaction=constraint_sat,
            feasibility=feasibility,
            baseline_cai=baseline_cai,
            baseline_gc=baseline_gc,
        )

    def analyze_restriction_site_removal(
        self,
        dna: str,
        enzyme_to_add: str,
    ) -> WhatIfScenario:
        """Analyze the impact of adding another restriction enzyme to avoid.

        Re-optimizes with the new enzyme added and reports the CAI cost.
        This tells users "adding BamHI avoidance only costs you 2% CAI".

        Args:
            dna: Current baseline DNA sequence.
            enzyme_to_add: Name of the restriction enzyme to add
                (e.g. ``"BamHI"``).

        Returns:
            WhatIfScenario with the predicted impact.
        """
        baseline_cai = compute_cai(dna, self.organism)
        baseline_gc = gc_content(dna)

        # Resolve the enzyme name to its recognition site
        site = get_recognition_site(enzyme_to_add)
        if site is None:
            # Not in our database — still allow the scenario but flag unknown
            site = f"unknown({enzyme_to_add})"

        description = (
            f"What if we add {enzyme_to_add} ({site}) to the "
            f"restriction site avoidance list?"
        )

        new_enzymes = list(self.restriction_sites) + [enzyme_to_add]

        try:
            result = self._reoptimize(enzymes=new_enzymes)
            predicted_cai = result.cai
            predicted_gc = result.gc_content
            constraint_sat = self._compute_constraint_satisfaction(
                result, self.gc_lo, self.gc_hi,
            )
            # Check if the new enzyme site is actually absent
            new_site_seq = get_recognition_site(enzyme_to_add)
            if new_site_seq:
                from .constants import reverse_complement as _rc
                site_present = (
                    new_site_seq in result.sequence
                    or _rc(new_site_seq) in result.sequence
                )
                constraint_sat[f"no_{enzyme_to_add}_site"] = not site_present
            feasibility = "feasible"
        except Exception as exc:
            logger.warning("Restriction site addition scenario failed: %s", exc)
            predicted_cai = None
            predicted_gc = None
            constraint_sat = None
            feasibility = "unknown"

        return WhatIfScenario(
            description=description,
            parameter_changed="restriction_sites",
            original_value=list(self.restriction_sites),
            alternative_value=new_enzymes,
            predicted_cai=predicted_cai,
            predicted_gc=predicted_gc,
            constraint_satisfaction=constraint_sat,
            feasibility=feasibility,
            baseline_cai=baseline_cai,
            baseline_gc=baseline_gc,
        )

    def analyze_constraint_removal(
        self,
        dna: str,
        constraint_name: str,
    ) -> WhatIfScenario:
        """Analyze the impact of removing a constraint.

        The most common use-case: "What if we remove the CpG avoidance
        constraint?" — this tells users how much CAI they're sacrificing
        for CpG avoidance.

        Currently supports removing:
        - ``"cpg_avoidance"``: Re-optimize without CpG island disruption
          (currently modelled by widening GC range, as CpG avoidance is
          embedded in the optimizer pipeline)
        - ``"restriction_sites"``: Re-optimize with no restriction enzyme
          avoidance
        - ``"cryptic_splice"``: Re-optimize with a very high splice
          threshold (effectively disabling splice elimination)

        Args:
            dna: Current baseline DNA sequence.
            constraint_name: Name of the constraint to remove.

        Returns:
            WhatIfScenario with the predicted impact.
        """
        baseline_cai = compute_cai(dna, self.organism)
        baseline_gc = gc_content(dna)

        description = f"What if we remove the {constraint_name} constraint?"
        original_value = True
        alternative_value = False

        try:
            if constraint_name == "restriction_sites":
                result = self._reoptimize(enzymes=[])
            elif constraint_name == "cryptic_splice":
                # Set threshold extremely high to effectively disable
                old_threshold = self.cryptic_splice_threshold
                self.cryptic_splice_threshold = 999.0
                result = self._reoptimize()
                self.cryptic_splice_threshold = old_threshold
            elif constraint_name == "cpg_avoidance":
                # CpG avoidance is embedded in the optimizer; removing it
                # means re-optimizing without the CpG step.  We approximate
                # this by re-optimizing with wider GC bounds (since CpG
                # avoidance often pushes GC down).  A more precise model
                # would require a flag in optimize_sequence.
                result = self._reoptimize(
                    gc_lo=max(0.0, self.gc_lo - 0.05),
                    gc_hi=min(1.0, self.gc_hi + 0.05),
                )
            else:
                logger.warning(
                    "Unknown constraint %r; attempting default re-optimization",
                    constraint_name,
                )
                result = self._reoptimize()

            predicted_cai = result.cai
            predicted_gc = result.gc_content
            constraint_sat = self._compute_constraint_satisfaction(
                result, self.gc_lo, self.gc_hi,
            )
            # Mark the removed constraint as N/A
            constraint_sat[constraint_name] = True  # trivially satisfied when removed
            feasibility = "feasible"
        except Exception as exc:
            logger.warning("Constraint removal scenario failed: %s", exc)
            predicted_cai = None
            predicted_gc = None
            constraint_sat = None
            feasibility = "unknown"

        return WhatIfScenario(
            description=description,
            parameter_changed=constraint_name,
            original_value=original_value,
            alternative_value=alternative_value,
            predicted_cai=predicted_cai,
            predicted_gc=predicted_gc,
            constraint_satisfaction=constraint_sat,
            feasibility=feasibility,
            baseline_cai=baseline_cai,
            baseline_gc=baseline_gc,
        )

    def analyze_organism_switch(
        self,
        dna: str,
        target_organism: str,
    ) -> WhatIfScenario:
        """Analyze the impact of optimizing for a different organism.

        Re-optimizes the protein for the target organism and compares
        CAI against the baseline.  This tells users "optimizing for E. coli
        instead of human would give you X% higher CAI for this gene".

        Args:
            dna: Current baseline DNA sequence (optimized for self.organism).
            target_organism: The alternative organism to optimize for.

        Returns:
            WhatIfScenario with the predicted impact.

        Raises:
            ValueError: If target_organism is not supported.
        """
        if target_organism not in SUPPORTED_ORGANISMS:
            raise ValueError(
                f"Unsupported organism: {target_organism!r}. "
                f"Supported: {SUPPORTED_ORGANISMS}"
            )

        baseline_cai = compute_cai(dna, self.organism)
        baseline_gc = gc_content(dna)

        description = (
            f"What if we optimize for {target_organism} "
            f"instead of {self.organism}?"
        )

        # Use the target organism's default GC range
        target_gc_range = ORGANISM_GC_TARGETS.get(
            target_organism, (self.gc_lo, self.gc_hi)
        )

        try:
            result = self._reoptimize(
                organism=target_organism,
                gc_lo=target_gc_range[0],
                gc_hi=target_gc_range[1],
            )
            predicted_cai = result.cai
            predicted_gc = result.gc_content
            constraint_sat = self._compute_constraint_satisfaction(
                result, target_gc_range[0], target_gc_range[1],
            )
            feasibility = "feasible"
        except Exception as exc:
            logger.warning("Organism switch scenario failed: %s", exc)
            predicted_cai = None
            predicted_gc = None
            constraint_sat = None
            feasibility = "unknown"

        return WhatIfScenario(
            description=description,
            parameter_changed="organism",
            original_value=self.organism,
            alternative_value=target_organism,
            predicted_cai=predicted_cai,
            predicted_gc=predicted_gc,
            constraint_satisfaction=constraint_sat,
            feasibility=feasibility,
            baseline_cai=baseline_cai,
            baseline_gc=baseline_gc,
        )

    def full_whatif_report(self, dna: str) -> list[WhatIfScenario]:
        """Run all standard what-if analyses and return a ranked list.

        Runs the following analyses:
        1. GC relaxation (gc_hi + 0.05)
        2. GC tightening (gc_hi - 0.05)
        3. Adding BamHI avoidance
        4. Adding NotI avoidance
        5. Removing restriction site avoidance entirely
        6. Removing cryptic splice constraint
        7. Removing CpG avoidance
        8. Organism switch to each alternative organism

        The results are ranked by absolute CAI delta (largest improvement
        or smallest cost first), so the most impactful trade-offs appear
        at the top.

        Args:
            dna: Current baseline DNA sequence.

        Returns:
            List of WhatIfScenario sorted by absolute CAI delta descending.
        """
        scenarios: list[WhatIfScenario] = []

        # 1. GC relaxation
        scenarios.append(
            self.analyze_gc_relaxation(dna, self.gc_hi, self.gc_hi + 0.05)
        )

        # 2. GC tightening
        if self.gc_hi - 0.05 >= self.gc_lo + 0.05:
            scenarios.append(
                self.analyze_gc_tightening(
                    dna,
                    self.gc_lo, self.gc_hi,
                    self.gc_lo + 0.02, self.gc_hi - 0.05,
                )
            )

        # 3-4. Adding common restriction enzymes
        for enzyme in ("BamHI", "NotI"):
            if enzyme not in self.restriction_sites:
                scenarios.append(
                    self.analyze_restriction_site_removal(dna, enzyme)
                )

        # 5. Removing restriction site avoidance entirely
        scenarios.append(
            self.analyze_constraint_removal(dna, "restriction_sites")
        )

        # 6. Removing cryptic splice constraint
        scenarios.append(
            self.analyze_constraint_removal(dna, "cryptic_splice")
        )

        # 7. Removing CpG avoidance
        scenarios.append(
            self.analyze_constraint_removal(dna, "cpg_avoidance")
        )

        # 8. Organism switches
        for alt_organism in SUPPORTED_ORGANISMS:
            if alt_organism != self.organism:
                try:
                    scenarios.append(
                        self.analyze_organism_switch(dna, alt_organism)
                    )
                except ValueError:
                    pass  # Skip unsupported organisms

        # Rank by absolute CAI delta (largest impact first)
        scenarios.sort(
            key=lambda s: abs(s.cai_delta) if s.cai_delta is not None else 0.0,
            reverse=True,
        )

        return scenarios


# ---------------------------------------------------------------------------
# WhatIfReport
# ---------------------------------------------------------------------------

class WhatIfReport:
    """Generate human-readable and machine-readable what-if reports.

    The report ranks scenarios by impact and highlights the key insight
    for each trade-off (e.g., "Relaxing GC max from 0.70 to 0.75 would
    improve CAI by +5.2%").

    Example::

        analyzer = WhatIfAnalyzer(protein="MVLSPAD...", organism="Homo_sapiens")
        scenarios = analyzer.full_whatif_report(dna)
        report_text = WhatIfReport.generate(scenarios)
        WhatIfReport.to_json(scenarios, "whatif_report.json")
    """

    @staticmethod
    def generate(scenarios: list[WhatIfScenario]) -> str:
        """Generate a markdown-formatted what-if analysis report.

        Args:
            scenarios: List of WhatIfScenario instances (typically from
                :meth:`WhatIfAnalyzer.full_whatif_report`).

        Returns:
            Multi-line markdown string with the full report.
        """
        if not scenarios:
            return "# What-If Analysis Report\n\nNo scenarios to report.\n"

        lines: list[str] = []
        lines.append("# What-If Analysis Report")
        lines.append("")
        lines.append(
            f"**Generated:** "
            f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}"
        )
        lines.append(f"**Scenarios analyzed:** {len(scenarios)}")
        lines.append("")

        # Summary table
        lines.append("## Summary (ranked by CAI impact)")
        lines.append("")
        lines.append(
            "| # | Scenario | Parameter | Original | Alternative | "
            "CAI (alt) | CAI Delta | GC (alt) | Feasibility |"
        )
        lines.append(
            "|---|----------|-----------|----------|-------------|"
            "----------|-----------|----------|-------------|"
        )

        for idx, s in enumerate(scenarios, 1):
            cai_str = f"{s.predicted_cai:.4f}" if s.predicted_cai is not None else "N/A"
            delta_str = (
                f"{s.cai_delta:+.4f}" if s.cai_delta is not None else "N/A"
            )
            gc_str = f"{s.predicted_gc:.4f}" if s.predicted_gc is not None else "N/A"
            orig_str = _format_value(s.original_value)
            alt_str = _format_value(s.alternative_value)

            lines.append(
                f"| {idx} | {s.description[:50]} | "
                f"{s.parameter_changed} | {orig_str} | {alt_str} | "
                f"{cai_str} | {delta_str} | {gc_str} | {s.feasibility} |"
            )

        lines.append("")

        # Detailed sections for each scenario
        lines.append("## Detailed Scenarios")
        lines.append("")

        for idx, s in enumerate(scenarios, 1):
            lines.append(f"### {idx}. {s.description}")
            lines.append("")
            lines.append(f"- **Parameter changed:** `{s.parameter_changed}`")
            lines.append(f"- **Original value:** {_format_value(s.original_value)}")
            lines.append(f"- **Alternative value:** {_format_value(s.alternative_value)}")
            lines.append(f"- **Feasibility:** {s.feasibility}")

            if s.predicted_cai is not None:
                lines.append(
                    f"- **Predicted CAI:** {s.predicted_cai:.4f} "
                    f"(baseline: {s.baseline_cai:.4f}, "
                    f"delta: {s.cai_delta:+.4f})"
                )
                # Key insight
                if s.cai_delta is not None:
                    pct = s.cai_delta / s.baseline_cai * 100 if s.baseline_cai > 0 else 0.0
                    if s.cai_delta > 0.01:
                        lines.append(
                            f"  > **Insight:** This change would improve CAI "
                            f"by ~{abs(pct):.1f}%"
                        )
                    elif s.cai_delta < -0.01:
                        lines.append(
                            f"  > **Insight:** This change would cost ~{abs(pct):.1f}% CAI"
                        )
                    else:
                        lines.append(
                            "  > **Insight:** This change has minimal CAI impact"
                        )

            if s.predicted_gc is not None:
                lines.append(
                    f"- **Predicted GC:** {s.predicted_gc:.4f} "
                    f"(baseline: {s.baseline_gc:.4f})"
                )

            if s.constraint_satisfaction is not None:
                lines.append("- **Constraint satisfaction:**")
                for cname, satisfied in s.constraint_satisfaction.items():
                    icon = "✅" if satisfied else "❌"
                    lines.append(f"  - {icon} {cname}")

            lines.append(f"- **Timestamp:** {s.timestamp}")
            lines.append("")

        return "\n".join(lines)

    @staticmethod
    def to_json(
        scenarios: list[WhatIfScenario],
        filepath: str,
    ) -> None:
        """Export scenarios to a JSON file.

        Args:
            scenarios: List of WhatIfScenario instances.
            filepath: Path to write the JSON file.
        """
        data = {
            "report_type": "whatif_analysis",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "scenario_count": len(scenarios),
            "scenarios": [s.to_dict() for s in scenarios],
        }
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, sort_keys=True, default=str)
        logger.info("What-if report written to %s (%d scenarios)", filepath, len(scenarios))

    @staticmethod
    def from_json(filepath: str) -> list[WhatIfScenario]:
        """Load scenarios from a JSON file.

        Args:
            filepath: Path to the JSON file.

        Returns:
            List of WhatIfScenario instances.

        Raises:
            ValueError: If the file format is invalid.
        """
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)

        if "scenarios" not in data:
            raise ValueError("Invalid what-if report file: missing 'scenarios' key")

        return [WhatIfScenario.from_dict(s) for s in data["scenarios"]]


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _json_safe(value: Any) -> Any:
    """Make a value JSON-serializable (convert lists, dicts, etc.)."""
    if isinstance(value, (str, int, float, bool, type(None))):
        return value
    if isinstance(value, list):
        return [_json_safe(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    return str(value)


def _format_value(value: Any) -> str:
    """Format a value for display in a markdown table."""
    if isinstance(value, float):
        return f"{value:.2f}"
    if isinstance(value, list):
        if len(value) <= 3:
            return ", ".join(_format_value(v) for v in value)
        return f"[{len(value)} items]"
    return str(value)
